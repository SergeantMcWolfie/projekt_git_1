import sys
from math import sin, cos, tan,  sqrt, atan, degrees, radians
from numpy import array


class Transformacje:
    def __init__(self, model: str = "wgs84"):
        """
        Parametry elipsoid:
            a - duża półoś elipsoidy - promień równikowy
            b - mała półoś elipsoidy - promień południkowy
            flat - spłaszczenie
            ecc2 - mimośród^2
        + WGS84: https://en.wikipedia.org/wiki/World_Geodetic_System#WGS84
        + Inne powierzchnie odniesienia: https://en.wikibooks.org/wiki/PROJ.4#Spheroid
        + Parametry planet: https://nssdc.gsfc.nasa.gov/planetary/factsheet/index.html
        """
        if model == "wgs84":
            self.a = 6378137.0 # semimajor_axis
            self.b = 6356752.31424518 # semiminor_axis
        elif model == "grs80":
            self.a = 6378137.0
            self.b = 6356752.31414036
        elif model == "mars":
            self.a = 3396900.0
            self.b = 3376097.80585952
        else:
            raise NotImplementedError(f"{model} model not implemented")
        self.flat = (self.a - self.b) / self.a
        self.ecc = sqrt(2 * self.flat - self.flat ** 2) # eccentricity  WGS84:0.0818191910428 
        self.ecc2 = (2 * self.flat - self.flat ** 2) # eccentricity**2
        
    def xyz2plh(self, X, Y, Z, output = 'dec_degree'):
        """
        Algorytm Hirvonena - algorytm transformacji współrzędnych ortokartezjańskich (x, y, z)
        na współrzędne geodezyjne długość szerokość i wysokośc elipsoidalna (phi, lam, h). Jest to proces iteracyjny. 
        W wyniku 3-4-krotneej iteracji wyznaczenia wsp. phi można przeliczyć współrzędne z dokładnoscią ok 1 cm.     
        Parameters
        ----------
        X, Y, Z : FLOAT
             współrzędne w układzie orto-kartezjańskim, 
    
        Returns
        -------
        lat
            [stopnie dziesiętne] - szerokość geodezyjna
        lon
            [stopnie dziesiętne] - długośc geodezyjna.
        h : TYPE
            [metry] - wysokość elipsoidalna
        output [STR] - optional, defoulf 
            dec_degree - decimal degree
            dms - degree, minutes, sec
        """
        r   = sqrt(X**2 + Y**2)           # promień
        lat_prev = atan(Z / (r * (1 - self.ecc2)))    # pierwsze przybliilizenie
        lat = 0
        while abs(lat_prev - lat) > 0.000001/206265:    
            lat_prev = lat
            N = self.a / sqrt(1 - self.ecc2 * sin(lat_prev)**2)
            h = r / cos(lat_prev) - N
            lat = atan((Z/r) * (((1 - self.ecc2 * N/(N + h))**(-1))))
        lon = atan(Y/X)
        N = self.a / sqrt(1 - self.ecc2 * (sin(lat))**2);
        h = r / cos(lat) - N       
        if output == "dec_degree":
            return degrees(lat), degrees(lon), h 
        elif output == "dms":
            lat = self.deg2dms(degrees(lat))
            lon = self.deg2dms(degrees(lon))
            return f"{lat[0]:02d}:{lat[1]:02d}:{lat[2]:.2f}", f"{lon[0]:02d}:{lon[1]:02d}:{lon[2]:.2f}", f"{h:.3f}"
        else:
            raise NotImplementedError(f"{output} - output format not defined")
    
    def plh2xyz(self, f, l, h):
        '''
        Algorytm odwrotny do tego zaproponowanego przez Hirvonena - transformuje współrzędne geocentryczne
        do współrzędnych ortokartezjańskich w postaci X, Y, Z.
    
        Parameters
        ----------
        f : float
            Wspolrzedna fi punktu (radiany)
        l : float
            Wspolrzedna lambda punktu (radiany)
        h : float
            Wysokosc elipsoidalna punktu (metry)
    
        Returns
        -------
        X, Y, Z : FLOAT
             [metry] współrzędne w układzie orto-kartezjańskim.
    
        '''
        N = self.a / sqrt((1 - self.ecc2 * sin(f)**2))
        X = (N + h) * cos(f) * cos(l)
        Y = (N + h) * cos(f) * sin(l)
        Z = ((N * (1 - self.ecc2)) + h) * sin(f) 
        return X, Y, Z
    
    def xyz2neu(self, X_st, Y_st, Z_st, X_end, Y_end, Z_end):
        '''
        Funkcja przekształca współrzędne orto-karrtezjańskie XYZ punktu początkowego
        do współrzędnych topocentrycznych NEU (northing, easting, up).
        
        Parameters
        ----------
        X_st, Y_st, Z_st : FLOAT
             [metry] współrzędne punktu początkowego (np. anteny) w układzie orto-kartezjańskim.
             
        X_end, Y_end, Z_end : FLOAT
             [metry] współrzędne punktu końcowego (np. staelity) w układzie orto-kartezjańskim.
    
        Returns
        -------
        n, e, u : FLOAT
             współrzędne w układzie topocentrycznym.
    
        '''
        r   = sqrt(X_st**2 + Y_st**2)           # promień
        lat_prev = atan(Z_st / (r * (1 - self.ecc2)))    # pierwsze przybliilizenie
        lat = 0
        while abs(lat_prev - lat) > 0.000001/206265:    
            lat_prev = lat
            N = self.a / sqrt(1 - self.ecc2 * sin(lat_prev)**2)
            h = r / cos(lat_prev) - N
            lat = atan((Z_st/r) * (((1 - self.ecc2 * N/(N + h))**(-1))))
        lon = atan(Y_st/X_st)
        N = self.a / sqrt(1 - self.ecc2 * (sin(lat))**2);
        h = r / cos(lat) - N
        dX = [X_st - X_end, Y_st - Y_end, Z_st - Z_end]
        R = array([[-sin(lat)*cos(lon), -sin(lat)*sin(lon), cos(lat)],
                      [-sin(lon), cos(lon), 0],
                      [cos(lat)*cos(lon), cos(lat)*sin(lon), sin(lat)]])
        dx = R @ dX
        n, e, u = dx[0], dx[1], dx[2]
        return n, e, u
    
    def bl2two(self, f, l, lb0):
        '''
        Funkcja zwraca wartosci wspolrzednych X, Y  w odwzorowaniu PL2000.
    
        Parameters
        ----------
        f : TYPE
            DESCRIPTION.
        l : TYPE
            DESCRIPTION.
        lb0 : TYPE
            DESCRIPTION.
    
        Returns
        -------
        x2000, y2000 : TYPE
            DESCRIPTION.
    
        '''
        if lb0 == 15:
            strefa = 5
        elif lb0 == 18:
            strefa = 6
        elif lb0 == 21:
            strefa = 7
        elif lb0 == 24:
            strefa = 8
        b2 = self.a**2 * (1 - self.ecc2)
        eeprim = (self.a**2 - b2) / b2
        deltalb = l - lb0
        t = tan(f)
        eta2 = eeprim * cos(f)**2
        N = self.a / sqrt((1 - self.ecc2 * sin(f)**2))
        A0 = 1 - self.ecc2/4 - (3*(self.ecc2**2))/64 - (5*(self.ecc2**3))/256
        A2 = (3/8)*(self.ecc2 + self.ecc2**2/4 + (15*(self.ecc2**3))/128)
        A4 = (15/256)*(self.ecc2**2 + (3*(self.ecc2**3))/4)
        A6 = (36*(self.ecc2**3))/3072
        sig = self.a*(A0*f - A2*sin(2*f) + A4*sin(4*f) - A6*sin(6*f))
        xgk = sig + (deltalb**2 / 2) * N * sin(f) * cos(f) * (1 + (deltalb**2 / 12) * cos(f)**2 * (5 - t**2 + 9*eta2 + 4*eta2**2)
                                                                    + (deltalb**4 / 360) * cos(f)**4 * (61 - 58*t**2 + t**4 + 270*eta2 - 330 * eta2 * t**2))
        ygk = deltalb * N * cos(f) * (1 + (deltalb**2 / 6) * cos(f)**2 * (1 - t**2 + eta2)
                                         + (deltalb**4 / 120) * cos(f)**4 * (5 - 18*t**2 + t**4 + 14*eta2 - 58 * eta2 * t**2))
        m2000 = 0.999923 
        x2000 = xgk * m2000
        y2000 = ygk * m2000 + strefa * 1000000 + 500000
        return x2000, y2000
    
    def bl2nine(self, f, l, lb0 = radians(19)):
        '''
        
    
        Parameters
        ----------
        f : TYPE
            DESCRIPTION.
        l : TYPE
            DESCRIPTION.
        lb0 : TYPE, optional
            DESCRIPTION. The default is radians(19).
    
        Returns
        -------
        x92 : TYPE
            DESCRIPTION.
        y92 : TYPE
            DESCRIPTION.
    
        '''
        b2 = self.a**2 * (1 - self.ecc2)
        eeprim = (self.a**2 - b2) / b2
        deltalb = l - lb0
        t = tan(f)
        eta2 = eeprim * cos(f)**2
        N = self.a / sqrt((1 - self.ecc2 * sin(f)**2))
        A0 = 1 - self.ecc2/4 - (3*(self.ecc2**2))/64 - (5*(self.ecc2**3))/256
        A2 = (3/8)*(self.ecc2 + self.ecc2**2/4 + (15*(self.ecc2**3))/128)
        A4 = (15/256)*(self.ecc2**2 + (3*(self.ecc2**3))/4)
        A6 = (36*(self.ecc2**3))/3072
        sig = self.a*(A0*f - A2*sin(2*f) + A4*sin(4*f) - A6*sin(6*f))
        xgk = sig + (deltalb**2 / 2) * N * sin(f) * cos(f) * (1 + (deltalb**2 / 12) * cos(f)**2 * (5 - t**2 + 9*eta2 + 4*eta2**2)
                                                                    + (deltalb**4 / 360) * cos(f)**4 * (61 - 58*t**2 + t**4 + 270*eta2 - 330 * eta2 * t**2))
        ygk = deltalb * N * cos(f) * (1 + (deltalb**2 / 6) * cos(f)**2 * (1 - t**2 + eta2)
                                         + (deltalb**4 / 120) * cos(f)**4 * (5 - 18*t**2 + t**4 + 14*eta2 - 58 * eta2 * t**2))
        m92 = 0.9993
        x92 = xgk * m92 - 5300000
        y92 = ygk * m92 + 500000
        
        return x92, y92


if __name__ == "__main__":
    if 'wgs84' in sys.argv[-4]:
        geo = Transformacje(model = "wgs84")
    elif 'grs80' in sys.argv[-4]:
        geo = Transformacje(model = "grs80")
    elif 'mars' in sys.argv[-4]:
        geo = Transformacje(model = "mars")
    else:
        raise 'Podaj model elipsoidy!'
    imp_file_path = sys.argv[-1]
    header_lines = int(sys.argv[-2])
    
    if "xyz2blh" in sys.argv[-3]:
        
    # XYZ TO BLH
        coords_plh = []
        with open(imp_file_path) as f:
            lines = f.readlines()
            lines = lines[header_lines:]
            for line in lines:
                line = line.strip()
                x_str, y_str, z_str = line.split(',')
                x, y, z = (float(x_str), float(y_str), float(z_str))
                p, l, h = geo.xyz2plh(x, y, z)
                coords_plh.append([p, l, h])
        with open('result_xyz2blh.txt', 'w') as f:
            # HEADER
            f.write('Podane współrzędne kartezjańskie XYZ zamienione na geocentryczne phi, lambda, height.\n')
            f.write(' phi[deg]   lam[deg]   h[m] \n')
            f.write('#---------------------------\n')
            # CONTENT
            for coords in coords_plh:
                line = ','.join([str(coord) for coord in coords])
                f.write(line + '\n')
                
    elif "xyz2neu" in sys.argv[-3]:
        
    # XYZ TO NEU
        coords_neu = []
        with open(imp_file_path) as f:
                lines = f.readlines()
                lines = lines[header_lines:]
                for line in lines:
                    line = line.strip()
                    x_str, y_str, z_str, x2_str, y2_str, z2_str = line.split(',')
                    x_st, y_st, z_st, x_end, y_end, z_end = (float(x_str), float(y_str), float(z_str), float(x2_str), float(y2_str), float(z2_str))
                    n, e, u = geo.xyz2neu(x_st, y_st, z_st, x_end, y_end, z_end)
                    coords_neu.append([n, e, u])
        with open('result_xyz2neu.txt', 'w') as f:
            # HEADER
            f.write('Podane współrzędne kartezjańskie XYZ punktu początkowego\n zamienione na topocentryczne NEU - northing, easting, up.\n')
            f.write(' Northing   Easting   Up \n')
            f.write('#---------------------------\n')
            # CONTENT
            for coords in coords_neu:
                line = ','.join([str(coord) for coord in coords])
                f.write(line + '\n')
            
    elif "blh2xyz" in sys.argv[-3]:
    # BLH TO XYZ
    
        coords_xyz = []
        with open(imp_file_path) as f:
                lines = f.readlines()
                lines = lines[header_lines:]
                for line in lines:
                    line = line.strip()
                    fi_str, lb_str, h_str = line.split(',')
                    fi, lb, h = (float(fi_str), float(lb_str), float(h_str))
                    x, y, z = geo.plh2xyz(radians(fi), radians(lb), h)
                    x, y, z = round(x,3), round(y, 3), round(z, 3)
                    coords_xyz.append([x, y, z])
        with open('result_blh2xyz.txt', 'w') as f:
            # HEADER
            f.write('Podane współrzędne geocentryczne lat, lon, h zamienione na kartezjańskie XYZ.\n')
            f.write(' X[m]   Y[m]   Z[m] \n')
            f.write('#---------------------------\n')
            # CONTENT
            for coords in coords_xyz:
                line = ','.join([str(coord) for coord in coords])
                f.write(line + '\n')
                
    elif "bl2two" in sys.argv[-3]:
    # BL TO XY2000
        
        coords_two = []
        with open(imp_file_path) as f:
                lines = f.readlines()
                lines = lines[header_lines:]
                for line in lines:
                    line = line.strip()
                    b_str, l_str = line.split(',')
                    b, l = (float(b_str), float(l_str))
                    x, y = geo.bl2two(b, l, lb0 = 15)
                    x, y = round(x, 3), round(y, 3)
                    coords_two.append([x, y])
        with open('result_bl2two.txt', 'w') as f:
            # HEADER
            f.write('Podane współrzędne kartezjańskie XYZ punktu początkowego\n zamienione na topocentryczne NEU - northing, easting, up.\n')
            f.write(' Northing   Easting   Up \n')
            f.write('#---------------------------\n')
            # CONTENT
            for coords in coords_two:
                line = ','.join([str(coord) for coord in coords])
                f.write(line + '\n')

    elif "bl2nine" in sys.argv[-3]:
    # BL TO XY1992    
        
        coords_two = []
        with open(imp_file_path) as f:
                lines = f.readlines()
                lines = lines[header_lines:]
                for line in lines:
                    line = line.strip()
                    b_str, l_str = line.split(',')
                    b, l = (float(b_str), float(l_str))
                    x, y = geo.bl2nine(b, l)
                    x, y = round(x, 3), round(y, 3)
                    coords_two.append([x, y])
        with open('result_bl2nine.txt', 'w') as f:
            # HEADER
            f.write('Podane współrzędne kartezjańskie XYZ punktu początkowego\n zamienione na topocentryczne NEU - northing, easting, up.\n')
            f.write(' Northing   Easting   Up \n')
            f.write('#---------------------------\n')
            # CONTENT
            for coords in coords_two:
                line = ','.join([str(coord) for coord in coords])
                f.write(line + '\n')
    
    
    else:
        raise "Musisz podać rodzaj transformacji współrzędnych!"
            
print(sys.argv)