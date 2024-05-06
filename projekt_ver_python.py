import sys
from math import sin, cos, tan,  sqrt, atan, atan2, degrees, radians
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
        
        def plh2xyz(self, p, l, h):
            '''
            Algorytm odwrotny do tego zaproponowanego przez Hirvonena - transformuje współrzędne geocentryczne
            do współrzędnych ortokartezjańskich w postaci X, Y, Z.

            Parameters
            ----------
            p : float
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
        
        def xyz2neu(self, X, Y, Z):
            '''
            Funkcja przekształca współrzędne orto-karrtezjańskie XYZ do współrzędnych topocentrycznych NEU (northing, easting, up).
            
            Parameters
            ----------
            X, Y, Z : FLOAT
                 [metry] współrzędne w układzie orto-kartezjańskim.

            Returns
            -------
            n, e, u : FLOAT
                 współrzędne w układzie topocentrycznym, 

            '''
            f, l, h = xyz2plh(self, X, Y, Z)
            R = array([[-sin(f)*cos(l), -sin(f)*sin(l), cos(f)],
                          [-sin(l), cos(l), 0],
                          [cos(f)*cos(l), cos(f)*sin(l), sin(f)]])
            dX = [X, Y, Z]
            dx = R @ dX
            n, e, u = dx[0], dx[1], dx[2]
            return n, e, u
        
        def bl2two(self, f, l, lb0):
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
        
        def bl2nine(self, b, l, lb0 = radians(19)):
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
    # utworzenie obiektu
    geo = Transformacje(model = "wgs84")
    imp_file_path = sys.argv[-1]
    coords_plh = []
    with open(imp_file_path) as f:
        lines = f.readlines()
        lines = lines[4:]
        for line in lines:
            line = line.strip()
            x_str, y_str, z_str = line.split(',')
            x, y, z = (float(x_str), float(y_str), float(z_str))
            p, l, h = geo.xyz2plh(x, y, z)
            coords_plh.append([p, l, h])
    with open('result_xyz2plh.txt', 'w') as f:
        f.write('phi[deg], lam[deg], h[m]\n')
        for coords in coords_plh:
            line = ','.join([str(coord) for coord in coords])
            f.write(line + '\n')