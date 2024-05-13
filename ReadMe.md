Transformacje_wsp_v5


Dostępne transformacje:

   XYZ ---> BLH
   BLH ---> XYZ
   XYZ ---> NEU
   BL ---> PL1992
   BL ---> PL2000

 
Obsługiwane modele elipsoidy:

   GRS80
   WGS84
   Elipsoida Krasowskiego
   Mars
 

Wymagania:

   system operacyjny Windows 11 
   python 3.12 
   biblioteka numpy
  


Opis programu:
 
 Plik przyjmuje argumenty podane za pomocą flag oddzielonych spacją w następującej kolejności:


'python' ---> 'Transformacje_wsp_v5' ---> transformacja ---> model elipsoidy --->  plik_wejsciowy.txt ---> liczba_wierszy_nagłówka_pliku_wejściowego   


  Wybór transformacji możliwy jest poprzez wpisanie jednej z poniższych nazw:

   'xyzblh' dla transformacji XYZ_BLH
   'blh2xyz' dla transformacji BLH_XYZ
   'xyz2neu' dla transformacji XYZ_NEU
   'bl2two' dla transformacji BL_PL1992
   'bl2nine' dla transformacji BL_PL2000
  

  Wybór modelu elipsoidy możliwy jest poprzez wpisanie jednej z poniższych nazw:

   'wgs84' dla elipsoidy WGS80
   'grs80' dla elipsoidy GRS84
   'krass' dla elipsoidy Krasowskiego
   'mars'  dla Mars
  
  Po wyborze parametrów i załadowaniu pliku z danymi utworzy się plik tekstowy zawierający wyniki wykonanych obliczeń.

  Plik z wynikami zapisuje się pod nazwą:

  "result_{funkcja}.txt"
  (gdzie {funkcja} oznacza nazwę transformacji, którą chcemy wykonać)


  Dodatkowo dla niektórych funkcji program zażąda podanie dodatkowych danych
 
 - po wywołaniu transformacji 'bl2two' należy podać jeden z południków osiowych: 15/18/21/24
 - po wywołaniu transformacji 'xyz2neu' należy podać współrzędne odbiornika w metrach, w następującej kolejności: X,Y,Z 
  (należy pamiętać aby oddzielić dane przecinkiem!)
  
  
Przykładowe transformacje (dla elipsoidy GRS80):
  

 XYZ ---> BLH
  Algorytm Hirvonena - algorytm transformacji współrzędnych ortokartezjańskich (X[m], Y[m], Z[m])
   na współrzędne geodezyjne długość szerokość i wysokośc elipsoidalna (phi, lam, h). Jest to proces iteracyjny. 
   W wyniku 3-4-krotneej iteracji wyznaczenia wsp. phi można przeliczyć współrzędne z dokładnoscią ok 1 cm. 

  dla danych z pliku 'wsp_inp_XYZ.txt' (kolejno X[m], Y[m], Z[m])

    3664940.500,1409153.590,5009571.170
    3664940.510,1409153.580,5009571.167

  otrzymujemy wyniki (kolejno fi, lambda, h)

    52.097272,21.031533,141.399
    52.097272,21.031533,141.4

  
 BLH ---> XYZ
  Algorytm odwrotny do tego zaproponowanego przez Hirvonena - transformuje współrzędne geocentryczne
   do współrzędnych ortokartezjańskich w postaci X, Y, Z.

  dla danych z pliku 'wsp_inp_BL.txt' (kolejno fi, lambda, h)

    52.097272,21.031533,141.399
    52.097272,21.031533,141.4

  otrzymujemy wyniki (kolejno X[m], Y[m], Z[m])

    3664940.526,1409153.576,5009571.155
    3664940.527,1409153.576,5009571.156


 XYZ ---> neu
  Funkcja przekształca współrzędne orto-kartezjańskie XYZ punktu początkowego
   do współrzędnych topocentrycznych NEU (northing, easting, up).

  dla danych z pliku 'wsp_XYZ.txt' (kolejno X[m], Y[m], Z[m])

    3664940.500,1409153.590,5009571.170
    3664940.510,1409153.580,5009571.167
  
  otrzymujemy wyniki (kolejno n, e, u)
  
    537311.974,613451.991,-881637.196
    537312.049,613451.933,-881637.196
  
 
 BL ---> XY PL1992
  Funkcja zwraca wartości współrzędnych X, Y  w odwzorowaniu PL2000.

  dla danych z pliku 'wsp_inp_BL.txt' (kolejno fi, lambda)
  
    52.097272,21.031533
    52.097272,21.031533

  otrzymujemy wyniki (kolejno x92[m], y92[m])

    472071.316,639114.469
    472071.316,639114.469

  
 BL ---> XY PL2000
  Funkcja zwraca wartości współrzędnych X, Y  w odwzorowaniu PL1992.
  
  dla danych z pliku 'wsp_inp_BL.txt' (kolejno fi, lambda)
  
    52.097272,21.031533
    52.097272,21.031533
  
  otrzymujemy wyniki (kolejno x2000[m], y2000[m])

    5773722.697,7502160.76
    5773722.697,7502160.76 

   
   



