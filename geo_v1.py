from math import sin, cos, sqrt, atan, atan2, degrees, radians, pi, tan
import sys
import argparse 
from argparse import ArgumentParser
import os
import tkinter as tk
from numpy import array as np


o = object()

class Transformacje:
    def __init__(self, model: str = "wgs84"):
        """
        Parametry elipsoid:
            a - duża półoś elipsoidy - promień równikowy
            b - mała półoś elipsoidy - promień południkowy
            flat - spłaszczenie
            ecc2 - mimośród^2
        
        """
        if model == "wgs84":
            self.a = 6378137.0  # semimajor_axis
            self.b = 6356752.31424518  # semiminor_axis
        elif model == "grs80":
            self.a = 6378137.0
            self.b = 6356752.31414036
        elif model == "Krasowski":
            self.a = 6378245.000
            self.e2 = 0.00669342162296
        else:
            raise NotImplementedError(f"{model} Ta powierzchnia odniesienia jest nieobługiwana")
        self.flat = (self.a - self.b) / self.a
        self.ecc = sqrt(2 * self.flat - self.flat ** 2)  # eccentricity  WGS84:0.0818191910428 
        self.ecc2 = (2 * self.flat - self.flat ** 2)  # eccentricity**2
        
    def Np(self,f):
        """
        Obliczenie promienia przekroju w pierwszym wertykale

        Parameters
        ----------
        f : FLOAT
            Szerokosc geodezyjna,[radiany]

        Returns
        -------
        N : FLOAT
            Promień przekroju w pierwszym wertykale, [metry]

        """
        N=self.a/sqrt(1-self.ep2*sin(f)**2)
        return N
  
    def sigma(self, f):
        """
        

        Parameters
        ----------
        f : FLOAT 
            Szerokosc geodezyjna.

        Returns
        sig: FLOAT

        """
        A0 = 1 - self.ep2/4 - 3 * self.ep2**2/64 - 5 * self.ep2**3/256
        A2 = (3/8) * (self.ep2 + self.ep2**2/4 + 15 * self.ep2**3/128)
        A4 = (15/256) * (self.ep2**2 + 3 * self.ep2**3/4)
        A6 = 35 * self.ep2**3/3072
        sig = self.a * (A0 * f - A2 * sin(2 * f) + A4 * sin(4 * f) - A6 * sin(6 * f))
        return(sig)
    
    def xyz2plh(self, X, Y, Z, output = 'dec_degree'):
        """
        Algorytm Hirvonena - algorytm transformacji współrzędnych ortokartezjańskich (x, y, z)
        na współrzędne geodezyjne długość szerokość i wysokośc elipsoidalna (phi, lam, h). Jest to proces iteracyjny. 
        W wyniku 3-4-krotnej iteracji wyznaczenia wsp. phi można przeliczyć współrzędne z dokładnoscią ok 1 cm.     
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
        r   = sqrt(X**2 + Y**2)  # promień
        lat_prev = atan(Z / (r * (1 - self.ecc2)))  # pierwsze przybliilizenie
        lat = 0
        while abs(lat_prev - lat) > 0.000001/206265:    
            lat_prev = lat
            N = self.a / sqrt(1 - self.ecc2 * sin(lat_prev)**2)
            h = r / cos(lat_prev) - N
            lat = atan((Z/r) * (((1 - self.ecc2 * N/(N + h))**(-1))))
        lon = atan(Y/X)
        N = self.a / sqrt(1 - self.ecc2 * (sin(lat))**2)
        h = r / cos(lat) - N       
        if output == "dec_degree":
            return degrees(lat), degrees(lon), h 
        elif output == "dms":
            lat = self.deg2dms(degrees(lat))
            lon = self.deg2dms(degrees(lon))
            return f"{lat[0]:02d}:{lat[1]:02d}:{lat[2]:.2f}", f"{lon[0]:02d}:{lon[1]:02d}:{lon[2]:.2f}", f"{h:.3f}"
        else:
            raise NotImplementedError(f"{output} - output format not defined")
            

    def plh2XYZ(self,phi,lam,h):
        """
        Przeliczenie współrzędnych prostokątnych na współrzędne ortokartezjańskie.

        Parameters
        ----------
        phi : FLOAT
            Szerokosc, [stopnie dziesietne]
        lam : FLOAT
            Dlugosc, [stopnie dziesietne]
        h : FLOAT
            Wysokosć, [metry]

        Returns
        result: LIST
        Lista zawierająca współrzędne ortokartezjańskie
        """
        result=[]
        phi=radians(phi)
        lam=radians(lam)
        Rn=self.a/sqrt(1-self.ecc2*sin(phi)**2)
        q=Rn*self.ecc2 *sin(phi)
        x=(Rn+h)*cos(phi)*cos(lam)
        y=(Rn+h)*cos(phi)*sin(lam)
        z=(Rn+h)*sin(phi)-q
        result.append([x, y, z])
        return result 
    
    def xyz2neu(self, x, y, z, x_0, y_0, z_0):
        '''
        Przeliczenie X,Y,Z na neu (układ współrzędnych horyzontaknych).

        Parameters
        ----------
        x, y, z, x_0, y_0, z_0: FLOAT [metry], współrzędne w układzie ortokartezjańskim          

        Returns
        -------
        result : LIST
            Funkcja zwraca listę ze współrzędnymi horyzontalnymi.

        '''
        result=[]
        phi, lam, h = [radians(coord) for coord in self.xyz2plh(x_0, y_0, z_0)]
        R = np.array([[-sin(lam), -sin(phi)*cos(lam), cos(phi)*cos(lam)],
                      [cos(lam), -sin(phi)*sin(lam), cos(phi)*sin(lam)],
                      [0, cos(phi), sin(phi)]])
        xyz_t = np.array([[x-x_0],
                          [y-y_0],
                          [z-z_0]])
        [[E],[N],[U]] = R.T @ xyz_t
        result.append(R.T @ xyz_t)
        return result
    
    
    def GK2000(self, f, l, m=0.999923):
        """
        Przeliczenie współrzędnych  do układu PL2000.

        Parameters
        ----------
        f : float
            DESCRIPTION. Szerokosc geograficzna.[stopnie dziesiętne]
        l : float
            DESCRIPTION. Długosc geograficzna. [stopnie dziesiętne]
        m : float, optional
            DESCRIPTION. The default is 0.999923. Skala odwzorowawcza. [bez jednostki]

        Returns
        -------
        result : list
            DESCRIPTION. Program zwraca listę zawierająca współrzędne X,Y w układzie PL2000.

        """
        result = []
        for f, l in zip(f,l):
            l0 = 0 
            strefa = 0
            if l > radians(13.5) and l < radians(16.5):
                strefa = 5
                l0 = radians(15)
            elif l > radians(16.5) and l < radians(19.5):
                strefa = 6
                l0 = radians(18)
            elif l > radians(19.5) and l < radians(22.5):
                strefa =7
                l0 = radians(21)
            elif l > radians(22.5) and l < radians(25.5):
                strefa = 8
                l0 = radians(24)
            b2 = (self.a**2) * (1-self.ep2)   #krotsza polos
            e2p = ( self.a**2 - b2 ) / b2   #drugi mimosrod elipsy
            dl = l - l0
            t = tan(f)
            ni = sqrt(e2p * (cos(f))**2)
            N = self.Np(f)
            sigma = self.sigma(f)
            XGK20 = sigma + ((dl**2)/2)*N*sin(f)*cos(f) * ( 1 + ((dl**2)/12)*(cos(f))**2 * ( 5 - (t**2)+9*(ni**2) + 4*(ni**4)     )  + ((dl**4)/360)*(cos(f)**4) * (61-58*(t**2)+(t**4) + 270*(ni**2) - 330*(ni**2)*(t**2))  )
            YGK20 = (dl*N* cos(f)) * (1+(((dl)**2/6)*(cos(f))**2) *(1-(t**2)+(ni**2))+((dl**4)/120)*(cos(f)**4)*(5-18*(t**2)+(t**4)+14*(ni**2)-58*(ni**2)*(t**2)) )
            X2000 = XGK20 * m 
            Y2000 = YGK20 * m + strefa*1000000 + 500000
            result.append([X2000, Y2000])
        
        return result
    
    
    
    def GK1992(self, f, l, m = 0.9993):
        """
        Przeliczenie współrzędnych do układu PL1992.

        Parameters
        ----------
        f : FLOAT
            Szerokosc geograficzna. [stopnie dziesiętne]
        l : FLOAT
            Długosc geograficzna. [stopnie dziesiętne]
        m : FLOAT, optional
            The default is 0.9993. Skala odwzorowawcza, [brak jednostki].

        Returns
        -------
        result : LIST
            Zwraca listę zawierjącą współrzędne X,Y w układzie PL1992.

        """
        result = []
        lam0 = radians(19)
        for f, l in zip(f,l):
            b2 = (self.a**2) * (1-self.ep2)   #krotsza polos
            e2p = ( self.a**2 - b2 ) / b2   #drugi mimosrod elipsy
            dlam = l - lam0
            t = tan(f)
            ni = sqrt(e2p * (cos(f))**2)
            N = self.Np(f)

            sigma = self.sigma(f)

            xgk = sigma + ((dlam**2)/2)*N*sin(f)*cos(f) * ( 1+ ((dlam**2)/12)*(cos(f))**2 * ( 5 - (t**2)+9*(ni**2) + 4*(ni**4)) + ((dlam**4)/360)*(cos(f)**4) * (61-58*(t**2)+(t**4) + 270*(ni**2) - 330*(ni**2)*(t**2)))
            ygk = (dlam*N* cos(f)) * (1+(((dlam)**2/6)*(cos(f))**2) *(1-(t**2)+(ni**2))+((dlam**4)/120)*(cos(f)**4)*(5-18*(t**2)+(t**4)+14*(ni**2)-58*(ni**2)*(t**2)) )
            
            x92 = xgk*m - 5300000
            y92 = ygk*m + 500000
            
            result.append([x92, y92])
        
        return result 

    def plik(self, file, transf):
        dane = np.genfromtxt(file,delimiter = ' ')
        if transf == 'XYZ2FLH':
            result  = self.XYZ2FLH(dane[:,0], dane[:,1], dane[:,2])
            np.savetxt(f"result_{transf}_{args.elip}.txt", result, delimiter=' ', fmt='%0.10f %0.10f %0.3f')
        elif transf == 'FLH2XYZ':
            result  = self.FLH2XYZ(np.deg2rad((dane[:,0])), np.deg2rad(dane[:,1]), dane[:,2])
            np.savetxt(f"result_{transf}_{args.elip}.txt", result, delimiter =' ', fmt ='%0.3f %0.3f %0.3f' )
        elif transf == 'XYZ2NEU':
            result  = self.XYZ2NEU(dane[1:,0], dane[1:,1], dane[1:,2], dane[0,0], dane[0,1], dane[0,2])
            np.savetxt(f"result_{transf}_{args.elip}.txt", result, delimiter =' ', fmt ='%0.3f %0.3f %0.3f' )
        elif transf == 'GK2000':
            result  = self.GK2000(np.deg2rad(dane[:,0]), np.deg2rad(dane[:,1]))
            np.savetxt(f"result_{transf}_{args.elip}.txt", result, delimiter=' ', fmt='%0.3f %0.3f')
        elif transf == 'GK1992':
            result  = self.GK1992(np.deg2rad(dane[:,0]), np.deg2rad(dane[:,1]))
            np.savetxt(f"result_{transf}_{args.elip}.txt", result, delimiter=' ', fmt='%0.3f %0.3f')

if __name__ == '__main__':
    parser = ArgumentParse()
    parser.add_argument('-dane', type=str, help='Wpisz sciezke do pliku z danymi wejsciowymi')
    parser.add_argument('-elip',gt type=str, help='Wybierz elipsoide sposrod dostepnych: WGS84, GRS80, KRASOWSKI')
    parser.add_argument('-transf', type=str, help='Wybierz transformacje, z ktorej chcesz skorzystac, sposrod dostepnych: XYZ2flh, flh2XYZ, saz2neu, GK2000, GK1992, XYZ2NEU')
    args = parser.parse_args()
    elip = {'WGS84': [6378137.000, 0.00669438002290], 'GRS80': [6378137.000, 0.00669438002290], 'KRASOWSKI': [6378245.000, 0.00669342162296]}
    transf = {'XYZ2FLH': 'XYZ2FLH', 'FLH2XYZ': 'FLH2XYZ','XYZ2NEU': 'XYZ2NEU', 'GK2000': 'GK2000', 'GK1992': 'GK1992'}

    
    try:
        wsp = Transformacje(elip[args.elip])
        wczyt = wsp.plik(args.dane, transf[args.transf.upper()])
        print("Utworzono plik ze wspolrzednymi.")
    except AttributeError as e:
        print("Error:", e)    
    except FileNotFoundError:
        print("Nie znaleziono podanego pliku")
    except KeyError:
        print("Niepoprawna nazwa Elipsoidy lub Transformacji")
    except IndexError:
        print("Dane w podanym pliku sa w nieodpowiednim formacie")
    except ValueError:
        print("Dane w podanym pliku sa w nieodpowiednim formacie")
    finally:
        print("Mamy nadzieję, że nasz program był dla Ciebie użyteczny")



        

        
    
