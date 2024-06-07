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

Obsługa programu w terminalu:
Uruchomienie programu wymaga podania parametrów przez wiersz poleceń.


Program przyjmuje argumenty podane za pomocą następujących flag
-m <nazwa_elipsoidy>: Nazwa modelu elipsoidy, na której chcemy dokonać transformacji
-trans <nazwa_transformacji>: Nazwa transformacji, którą chcemy wykonać
-i <input, podejmy nazwę pliku tekstowego, z którego program ma wczytać dane>

Po fladze -m wpisujemy jedną z poniższych nazw:
'wgs84'
'grs80'
'Krasowski'

Po fladze -trans wpisujemy jedną z poniższych nazw (transformacje, która chcemy wykonać)
'xyz2plh;
'neu'
'flh22000_92'
'plh2xyz'


PRZYKŁADOWY POPRAWNY ZAPIS WYWOŁANIA:
python geo_v1.py -m grs80 -trans plh2xyz -i wsp_inp.txt


Po wpisaniu poprawnie wyśtwiela się komunikat informujący nas o:
- "Dziękujemy za skorzystanie z naszego programu"
W przypadku wyboru neu: 
"Plik ze współrzędnymi neu został zapisany Wyniki_neu.txt
Dziękujemy za skorzystanie z naszego programu"

Dla każdej transformacji zapisywany jest plik tekstowy z wynikami.

ZNANE BŁĘDY
- Zawarte w kodzie transformacje nie działają, gdy użytkownik próbuje wykonać je na elipsoidzie Krasowskiego.
- Wykonanie transformacji xyz2plh zapisuje plik z odpowiednimi wynikami oraz dodatkowo wyniki transformacji dla neu.
- Program działa tylko na wczytywanie pliku wsp_inp.txt
