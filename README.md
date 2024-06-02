PROJEKT NR 1: TRANSFORMACJE

Obsługiwane transformacje
XYZ ---> BLH
BLH ---> XYZ
XYZ ---> NEU
BL ---> PL1992
BL ---> PL2000

Obsługiwane elipsoidy
GRS80
WGS84
Elipsoida Krasowskiego

Wymagania programu
system operacyjny Window 11 lub macOS 11.7.6
Python 3.11.5
biblioteka numpy
bibliteka argparse
biblioteka math (sin, cos, sqrt, atan, atan2, degrees, radians, pi, tan)
biblioteka sys

Na początku funkcja \texttt{wczytanie\_zapisanie\_pliku} odczytuje dane z pliku wejściowego, który zawiera współrzędne X, Y, Z. Następnie używamy funkcji \texttt{xyz2plh} do przeliczenia tych współrzędnych na długość geograficzną, szerokość geograficzną i wysokość elipsoidalną. W przypadku, gdy długość i szerokość mieszczą się w odpowiednich zakresach dla układu PL1992 lub PL2000, obliczamy również współrzędne w tych układach przy użyciu funkcji \texttt{flh2PL92} i \texttt{flh2PL00}.

Po obliczeniach funkcja przekształca współrzędne XYZ na współrzędne NEU za pomocą funkcji \texttt{xyz2neu}. Następnie wyniki są zapisywane do dwóch oddzielnych plików tekstowych: jeden zawierający współrzędne XYZ, długość i szerokość geograficzną, wysokość elipsoidalną oraz współrzędne w układach PL1992 i PL2000, a drugi zawierający współrzędne NEU.

Obsługa programu w terminalu:
Uruchomienie programu wymaga podania parametrów przez wiersz poleceń.


Program przyjmuje argumenty podane za pomocą następujących flag
-m <nazwa_elipsoidy>: Nazwa modelu elipsoidy, na której chcemy dokonać transformacji
-trans <nazwa_transformacji>: Nazwa transformacji, którą chcemy wykonać

Wybór elipsoidy możliwy jest poprzez wpisanie jednej z poniższych nazw:
'wgs84'
'grs80'
'Krasowski'


Natępnie przkazujemy parametry w kolejności, oznaczone nastepującymi flagami:
-fa, -la, -ha, -fb, -lb, -hb
*fa - szerokość geodezyjna punktu A
*la - długość geodezyjna punktu A
*ha - wysokość punktu A
*fb - szerokość geodezyjna punktu B
*lb - długość geodezyjna punktu B
*hb - wysokość punktu B


PRZYKŁADOWY POPRAWNY ZAPIS WYWOŁANIA:
python geo_v1.py -m grs80 -trans wynik.txt -fa 53 -la 67 -ha 54 -fb 56 -lb 68 -hb 115

PRZYKŁADOWY PLIK WYNIKOWY:
Wyniki_obliczen_Geodezyjnych; n, e, u.

|                        n                         |                        e                         |                        u                         |

|                                0.0116965954423486|                               -0.0072036099442109|                                0.0041590449204018|
|                                0.0142304947739723|                               -0.0024513568876594|                                0.0006956011613583|
|                                0.0134563734931296|                               -0.0024513598954034|                                0.0035940022406434|
|                                0.0131983327569383|                               -0.0024513621404213|                                0.0045601348965502|
|                               -0.0297971735089245|                               -0.0090253286550346|                               -0.0058882843057203|
|                                0.0024531660178090|                                0.0083568357704472|                               -0.0034850066630017|
|                                0.0132195394281337|                               -0.0038441650789937|                                0.0035307466986917|
|                                0.0123646853515949|                               -0.0002395847820951|                                0.0064075853461758|
|                               -0.0308093993992371|                               -0.0001782219649152|                               -0.0071936874040200|
|                                0.0044510530891616|                                0.0042619814678192|                                0.0001537584044234|
|                               -0.0081491917165731|                               -0.0071311561061653|                               -0.0021765307194158|
|                                0.0070787082016427|                                0.0079501283073926|                                0.0039607261207125|


Po wpisaniu poprawnie wyśtwiela się komunikat informujący nas o:
-wybranej elipsoidzie 
-wynikach: w formie Wyniki_z_flh2neu; n = ..., e = ..., u = ...
- "Dziękujemy za skorzystanie z naszego programu"

ZNANE BŁĘDY
- Zawarte w kodzie transformacje nie działają, gdy użytkownik próbuje wykonać je na elipsoidzie Krasowskiego.
