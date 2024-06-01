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

Program przyjmuje argumenty podane za pomocą następujących flag
-m <nazwa_elipsoidy>: Nazwa modelu elipsoidy, na której chcemy dokonać transformacji
-trans <nazwa_transformacji>: Nazwa transformacji, którą chcemy wykonać

Wybór elipsoidy możliwy jest poprzez wpisanie jednej z poniższych nazw:
'wgs84'
'grs80'
'Krasowski'


Wybór transformacji możliwy jest poprzez wpisanie:
'neu'

Natępnie przkazujemy parametry w kolejności, oznaczone nastepującymi flagami:
-fa, -la, -ha, -fb, -lb, -hb
*fa - szerokość geodezyjna punktu A
*la - długość geodezyjna punktu A
*ha - wysokość punktu A
*fb - szerokość geodezyjna punktu B
*lb - długość geodezyjna punktu B
*hb - wysokość punktu B


PRZYKŁADOWY POPRAWNY ZAPIS WYWOŁANIA:
python geo_v1.py -m grs80 -trans neu wynik.txt -fa 53 -la 67 -ha 54 -fb 56 -lb 68 -hb 115

Po wpisaniu poprawnie wyśtwiela się komunikat informujący nas o:
-wybranej elipsoidzie 
-wynikach: w formie Wyniki_z_flh2neu; n = ..., e = ..., u = ...
- "Dziękujemy za skorzystanie z naszego programu"

ZNANE BŁĘDY
- Zawarte w kodzie transformacje nie działają, gdy użytkownik próbuje wykonać je na elipsoidzie Krasowskiego.
