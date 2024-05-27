Transformacje_wsp_v5


Dostępne transformacje:<br>

   XYZ ---> BLH<br>
   BLH ---> XYZ<br>
   XYZ ---> NEU<br>
   BL ---> PL1992<br>
   BL ---> PL2000<br>

 
Obsługiwane modele elipsoidy:<br>

   GRS80<br>
   WGS84<br>
   Elipsoida Krasowskiego<br>
   Mars<br>
 

Wymagania:<br>

   system operacyjny Windows 11 <br>
   python 3.12 lub 3.11 <br>
   biblioteka numpy<br>
  


Opis programu:<br>
 
 Plik przyjmuje argumenty podane za pomocą flag oddzielonych spacją w następującej kolejności:<br>


'python' ---> 'Transformacje_wsp_v5' ---> transformacja ---> model elipsoidy --->  plik_wejsciowy.txt ---> liczba_wierszy_nagłówka_pliku_wejściowego   <br>


  Wybór transformacji możliwy jest poprzez wpisanie jednej z poniższych nazw:<br>

   'xyzblh' dla transformacji XYZ_BLH<br>
   'blh2xyz' dla transformacji BLH_XYZ<br>
   'xyz2neu' dla transformacji XYZ_NEU<br>
   'bl2nine' dla transformacji BL_PL1992<br>
   'bl2two' dla transformacji BL_PL2000<br>
  

  Wybór modelu elipsoidy możliwy jest poprzez wpisanie jednej z poniższych nazw:<br>

   'wgs84' dla elipsoidy WGS84<br>
   'grs80' dla elipsoidy GRS80<br>
   'krass' dla elipsoidy Krasowskiego<br>
   'mars'  dla Mars<br>

   Przy czym użycie elipsoidy Krasowskiego możliwe jest jedynie dla transformacji: 'xyzblh', 'blh2xyz', 'xyz2neu'<br>
   (program nie obsługuje w tym wypadku transformacji: 'bl2nine', 'bl2two' - zwróci wtedy błędne wyniki!)<br>
    
  
  Po wyborze parametrów i załadowaniu pliku z danymi utworzy się plik tekstowy zawierający wyniki wykonanych obliczeń.<br>

  Plik z wynikami zapisuje się pod nazwą:<br>

  "result_{funkcja}.txt"<br>
  (gdzie {funkcja} oznacza nazwę transformacji, którą chcemy wykonać)<br>


  Dodatkowo dla niektórych funkcji program zażąda podanie dodatkowych danych<br>
 
 - po wywołaniu transformacji 'bl2two' należy podać jeden z południków osiowych: 15/18/21/24<br>
 - po wywołaniu transformacji 'xyz2neu' należy podać współrzędne odbiornika w metrach, w następującej kolejności: X,Y,Z <br>
  (należy pamiętać aby oddzielić dane przecinkiem!)<br>
  
  
Przykładowe transformacje (dla elipsoidy GRS80):<br>
  

 XYZ ---> BLH<br>
  Algorytm Hirvonena - algorytm transformacji współrzędnych ortokartezjańskich (X[m], Y[m], Z[m])
   na współrzędne geodezyjne długość szerokość i wysokośc elipsoidalna (phi, lam, h). Jest to proces iteracyjny. 
   W wyniku 3-4-krotneej iteracji wyznaczenia wsp. phi można przeliczyć współrzędne z dokładnoscią ok 1 cm. <br>

  dla danych z pliku 'wsp_inp_XYZ.txt' (kolejno X[m], Y[m], Z[m])

    3664940.500,1409153.590,5009571.170
    3664940.510,1409153.580,5009571.167

  otrzymujemy wyniki (kolejno fi, lambda, h)

    52.097272,21.031533,141.399
    52.097272,21.031533,141.4

  
 BLH ---> XYZ<br>
  Algorytm odwrotny do tego zaproponowanego przez Hirvonena - transformuje współrzędne geocentryczne
   do współrzędnych ortokartezjańskich w postaci X, Y, Z.<br>

  dla danych z pliku 'wsp_inp_BL.txt' (kolejno fi, lambda, h)

    52.097272,21.031533,141.399
    52.097272,21.031533,141.4

  otrzymujemy wyniki (kolejno X[m], Y[m], Z[m])

    3664940.526,1409153.576,5009571.155
    3664940.527,1409153.576,5009571.156


 XYZ ---> neu<br>
  Funkcja przekształca współrzędne orto-kartezjańskie XYZ punktu początkowego
   do współrzędnych topocentrycznych NEU (northing, easting, up).<br>

  dla danych z pliku 'wsp_XYZ.txt' (kolejno X[m], Y[m], Z[m])

    3664940.500,1409153.590,5009571.170
    3664940.510,1409153.580,5009571.167
  
  otrzymujemy wyniki (kolejno n, e, u)
  
    537311.974,613451.991,-881637.196
    537312.049,613451.933,-881637.196
  
 
 BL ---> XY PL1992<br>
  Funkcja zwraca wartości współrzędnych X, Y  w odwzorowaniu PL1992.<br>

  dla danych z pliku 'wsp_inp_BL.txt' (kolejno fi, lambda)
  
    52.097272,21.031533
    52.097272,21.031533

  otrzymujemy wyniki (kolejno x92[m], y92[m])

    472071.316,639114.469
    472071.316,639114.469

  
 BL ---> XY PL2000<br>
  Funkcja zwraca wartości współrzędnych X, Y  w odwzorowaniu PL2000.<br>
  
  dla danych z pliku 'wsp_inp_BL.txt' (kolejno fi, lambda)
  
    52.097272,21.031533
    52.097272,21.031533
  
  otrzymujemy wyniki (kolejno x2000[m], y2000[m])

    5773722.697,7502160.76
    5773722.697,7502160.76 

   
   



