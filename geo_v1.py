from math import sin, cos, sqrt, atan, atan2, degrees, radians, pi, tan
import sys
import argparse 
from argparse import ArgumentParser
import os
import tkinter as tk
import numpy as np



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
    
    def deg2dms(self, x):
        '''
        Funkcja dms służy do zamienienia wizualnego FLOAT (radiany) na stopnie minuty i sekundy.   

        Parametry
        ----------
        x : FLOAT
        [radiany].

        Returns
        x : STR
        [dms] - stopnie, minuty, sekundy
        '''
        sig = ' '
        if x < 0:
            sig = '-'
            x = abs(x)
        x = x * 180/np.pi
        d = int(x)
        m = int(60 * (x - d))
        s = (x - d - m/60)*3600
        if s > 59.999995:
            s = 0
            m = m + 1
        if m == 60:
            m = 0
            d = d + 1
            
        d = str(d)
        if len(d) == 1:
            d = "  " + d
        elif len(d) == 2:
            d = " " + d
        elif len(d) == 3:
            d = d
                
        if m < 10:
            m = str(m)
            m = "0" + m
                
        if s < 10:
            s = "%.5f"%s
            s = str(s)
            s= "0" + s
        else:
            s = ("%.5f"%s)
                
        x1=(f'{d}°{m}′{s}″')  
        return(x1)
        
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
        N=self.a/sqrt(1-self.ecc2*sin(f)**2)
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
            return (lat,lon,h)
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
        phi=radians(phi)
        lam=radians(lam)
        Rn=self.a/sqrt(1-self.ecc2*sin(phi)**2)
        q=Rn*self.ecc2 *sin(phi)
        x=(Rn+h)*cos(phi)*cos(lam)
        y=(Rn+h)*cos(phi)*sin(lam)
        z=(Rn+h)*sin(phi)-q
        return x,y,z
    
    def xyz2neu(self, f, l, xa, ya, za, xb, yb, zb):
        '''
        Układ współrzędnych horyzontalnych – układ współrzędnych astronomicznych, w którym oś główną stanowi 
        lokalny kierunek pionu, a płaszczyzną podstawową jest płaszczyzna horyzontu astronomicznego. 
        Biegunami układu są zenit i nadir. Ich położenie na sferze niebieskiej zależy od współrzędnych geograficznych 
        obserwatora oraz momentu obserwacji, tak więc współrzędne horyzontalne opisują jedynie chwilowe położenie ciała niebieskiego.

        Parametry
        ----------
        f : FLOAT
            [stopnie dziesiętne] - szerokość geodezyjna..
        l : FLOAT
            [stopnie dziesiętne] - długośc geodezyjna.
        XA, YA, ZA, XB, YB, ZB: FLOAT
             współrzędne w układzie orto-kartezjańskim, 

        Returns
        -------
        n , e, u : STR
            współrzędne horyzontalne
            

        '''
        dX = Transformacje.dXYZ(self, xa, ya, za, xb, yb, zb)
        R = np.array([[-np.sin(f)*np.cos(l), -np.sin(l), np.cos(f)*np.cos(l)],
                      [-np.sin(f)*np.sin(l),  np.cos(l), np.cos(f)*np.sin(l)],
                      [np.cos(f),             0,         np.sin(f)          ]])
        neu = R.T @ dX
        n = neu[0];   e = neu[1];   u = neu[2]
        n = "%.16f"%n; e = "%.16f"%e; u="%.16f"%u
        dlugosc = []
        xx = len(n); dlugosc.append(xx)
        yy = len(e); dlugosc.append(yy)
        zz = len(u); dlugosc.append(zz)
        P = 50
        
        while xx < P:
            n = str(" ") + n
            xx += 1
        
        while yy < P:
            e = str(" ") + e
            yy += 1
            
        while zz < P:
            u = str(" ") + u
            zz +=1
            
        return(n, e, u)
    
    def dXYZ(self, xa, ya, za, xb, yb, zb):
        '''
        funkcja liczy macierz różnicy współrzednych punktów A i B, która jest potrzebna do obliczenia macierzy neu

        Parametry
        ----------
        XA, YA, ZA, XB, YB, ZB: FLOAT
             współrzędne w układzie orto-kartezjańskim, 

        Returns
        -------
        dXYZ : ARRAY
            macierz różnicy współrzędnych

        '''
        dXYZ = np.array([xb-xa, yb-ya, zb-za])
        return(dXYZ)
    
    def flh2PL00(self, f, l):
        '''
        Układ współrzędnych 2000 – układ współrzędnych płaskich prostokątnych zwany układem „2000”, 
        powstały w wyniku zastosowania odwzorowania Gaussa-Krügera dla elipsoidy GRS 80 w czterech 
        trzystopniowych strefach o południkach osiowych 15°E, 18°E, 21°E i 24°E, oznaczone odpowiednio numerami – 5, 6, 7 i 8.

        Parametry
        ----------
        f : FLOAT
            [stopnie dziesiętne] - szerokość geodezyjna..
        l : FLOAT
            [stopnie dziesiętne] - długośc geodezyjna.

        Returns
        -------
         X2000, Y2000 : FLOAT
              [metry] - współrzędne w układzie 2000

        '''
          
        if l >= 13.5 and l < 16.5:
            l0 = np.radians(15)
        elif l >= 16.5 and l < 19.5:
            l0 = np.radians(18)
        elif l >= 19.5 and l < 22.5:
            l0 = np.radians(21)
        elif l >= 22.5 and l <= 25.5:
            l0 = np.radians(24)
        else:
            raise NotImplementedError(f"{Transformacje.dms(self, np.radians(l))} ten południk nie jest obsługiwany przez układ współrzędnych płaskich PL2000")
        
        if f > 55 or f < 48.9:
            raise NotImplementedError(f"{Transformacje.dms(self, np.radians(f))} ten równoleżnik nie jest obsługiwany przez układ współrzędnych płaskich PL2000")
            
        f = np.radians(f)
        l = np.radians(l)
        a2 = self.a**2
        b2 = a2 * (1 - self.ecc2)
        e_2 = (a2 - b2)/b2
        dl = l - l0
        dl2 = dl**2
        dl4 = dl**4
        t = np.tan(f)
        t2 = t**2
        t4 = t**4
        n2 = e_2 * (np.cos(f)**2)
        n4 = n2 ** 2
        N = Transformacje.Np(self, f)
        e4 = self.ecc2**2
        e6 = self.ecc2**3
        A0 = 1 - (self.ecc2/4) - ((3*e4)/64) - ((5*e6)/256)
        A2 = (3/8) * (self.ecc2 + e4/4 + (15*e6)/128)
        A4 = (15/256) * (e4 + (3*e6)/4)
        A6 = (35*e6)/3072
        sigma = self.a * ((A0 * f) - A2 * np.sin(2*f) + A4 * np.sin(4*f) - A6 * np.sin(6*f))
        xgk = sigma + ((dl**2)/2) * N * np.sin(f) * np.cos(f) * (1 + ((dl**2)/12)*(np.cos(f)**2)*(5 - t2 + 9 * n2 + 4 * n4) + (dl4/360) * (np.cos(f)**4)*(61 - (58 * t2) + t4 + (270 * n2) - (330 * n2 * t2)))
        ygk = dl * N * np.cos(f) * (1 + (dl2/6) * (np.cos(f)**2) * (1 - t2 + n2) + (dl4/120) * (np.cos(f)**4) * (5 - (18 * t2) + t4 + (14 * n2) - 58 * n2 * t2))
        strefa = round(l0 * 180/np.pi)/3
        x00 = xgk * 0.999923
        y00 = ygk * 0.999923 + strefa * 1000000 + 500000
        return(x00,y00)
    
    
    
    def flh2PL92(self, f, l):
        '''
        Układ współrzędnych 1992 (Państwowy Układ Współrzędnych Geodezyjnych 1992) – układ współrzędnych 
        płaskich prostokątnych oparty na odwzorowaniu Gaussa-Krügera dla elipsoidy GRS80 w jednej dziesięciostopniowej strefie.

        Parametry
        ----------
        f : FLOAT
            [stopnie dziesiętne] - szerokość geodezyjna..
        l : FLOAT
            [stopnie dziesiętne] - długośc geodezyjna.

        Returns
        -------
         X1992, Y1992 : FLOAT
              [metry] - współrzędne w układzie 1992

        '''
        
        if l > 25.5 or l < 13.5:
            raise NotImplementedError(f"{Transformacje.dms(self, np.radians(l))} ten południk nie jest obsługiwany przez układ współrzędnych płaskich PL1992")
            
        if f > 55 or f < 48.9:
            raise NotImplementedError(f"{Transformacje.dms(self, np.radians(f))} ten równoleżnik nie jest obsługiwany przez układ współrzędnych płaskich PL1992")
            
        f = np.radians(f)
        l = np.radians(l)
        a2 = self.a**2
        b2 = a2 * (1 - self.ecc2)
        e_2 = (a2 - b2)/b2
        l0 = np.radians(19)
        dl = l - l0
        dl2 = dl**2
        dl4 = dl**4
        t = np.tan(f)
        t2 = t**2
        t4 = t**4
        n2 = e_2 * (np.cos(f)**2)
        n4 = n2 ** 2
        N = Transformacje.Np(self, f)
        e4 = self.ecc2**2
        e6 = self.ecc2**3
        A0 = 1 - (self.ecc2/4) - ((3*e4)/64) - ((5*e6)/256)
        A2 = (3/8) * (self.ecc2 + e4/4 + (15*e6)/128)
        A4 = (15/256) * (e4 + (3*e6)/4)
        A6 = (35*e6)/3072
        sigma = self.a * ((A0 * f) - A2 * np.sin(2*f) + A4 * np.sin(4*f) - A6 * np.sin(6*f))
        xgk = sigma + ((dl**2)/2) * N * np.sin(f) * np.cos(f) * (1 + ((dl**2)/12)*(np.cos(f)**2)*(5 - t2 + 9 * n2 + 4 * n4) + (dl4/360) * (np.cos(f)**4)*(61 - (58 * t2) + t4 + (270 * n2) - (330 * n2 * t2)))
        ygk = dl * N * np.cos(f) * (1 + (dl2/6) * (np.cos(f)**2) * (1 - t2 + n2) + (dl4/120) * (np.cos(f)**4) * (5 - (18 * t2) + t4 + (14 * n2) - 58 * n2 * t2))
        x92 = xgk * 0.9993 - 5300000
        y92 = ygk * 0.9993 + 500000
        return(x92,y92)
    
    def zamiana_float2string_rad(self, liczba):
        '''
        zamienia z float na str z dokładnocia taką jaką musi mieć jednostka radiany. 
        Zrobione jest to po to żeby wynik się ładnie zapisywał w pliku wynikowym 

        Parameters
        ----------
        liczba : FLOAT
            Liczba w radianach.

        Returns
        -------
        liczba : STR
            string z okreloną stała iloscia znaków

        '''
        zm_liczba = "%.12f"%liczba
        P = 16
        xx = len(zm_liczba)
        while xx < P:
            zm_liczba = str(" ") + zm_liczba
            xx += 1
        return(zm_liczba)  

    
    def zamiana_float2string_fl(self, liczba):
        '''
        zamienia z float na str z dokładnocia taką jaką musi mieć jednostka stopnie dziesiętne. 
        Zrobione jest to po to żeby wynik się ładnie zapisywał w pliku wynikowym 

        Parameters
        ----------
        liczba : FLOAT
            Liczba w stopniach dziesiętnych.

        Returns
        -------
        liczba : STR
            string z okreloną stała iloscia znaków

        '''
        zm_liczba = "%.10f"%liczba
        P = 16
        xx = len(zm_liczba)
        while xx < P:
            zm_liczba = str(" ") + zm_liczba
            xx += 1
        return(zm_liczba)        
    
    
    def zamiana_float2string(self, liczba):
        '''
        zamienia z float na str z dokładnocia taką jaką musi mieć jednostka dms. 
        Zrobione jest to po to żeby wynik się ładnie zapisywał w pliku wynikowym 

        Parameters
        ----------
        liczba : FLOAT
            Liczba w dms.

        Returns
        -------
        liczbe : STR
            string z okreloną stała iloscia znaków


        '''
        zm_liczba = "%.3f"%liczba
        P = 21
        xx = len(zm_liczba)
        while xx < P:
            zm_liczba = str(" ") + zm_liczba
            xx += 1
        return(zm_liczba)
    '''
    def plik(self, file, transf):
        dane = np.genfromtxt(file,delimiter = ',')
        if transf == 'XYZ2FLH':
            result  = self.XYZ2FLH(dane[:,0], dane[:,1], dane[:,2])
            np.savetxt(f"result_{transf}_{args.elip}.txt", result, delimiter='', fmt='%0.10f %0.10f %0.3f')
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
            '''            
    def zapisanie_pliku(self, X, Y, Z, f, l, h, x92, y92, x00, y00, N, E, U, xyz_txt, neu_txt): 
        with open(xyz_txt, "w",  encoding="utf-8") as plik:
            plik.write("Wyniki obliczeń Geodezyjnych; X, Y, Z, fi, lambda, h, x1992, y1992, x2000, y2000.\n")
            plik.write("Znak '-' w koordynatach; x1992, y1992, x2000, y2000 oznacza, że dla podanych współrzędnych ortokartezjańskich (X, Y, Z) po obliczeniu współrzędnych geodezyjnych fi i lambda. fi i lambda nie należą do dozwolonych współrzędnych \ngeodezyjnych układów PL1992, PL2000.\n")
            plik.write("|          X          |          Y          |          Z          |          fi         |        lambda       |          h          |        x1992        |        y1992        |        x2000        |        y2000        |\n")
            for x, y, z, f, l, h, x92, y92, x00, y00 in zip(X, Y, Z, f, l, h, x92, y92, x00, y00):
                plik.write(f"|{Transformacje.zamiana_float2string(self, x)}|{Transformacje.zamiana_float2string(self, y)}|{Transformacje.zamiana_float2string(self, z)}|     {f}|     {l}|{h}|{x92}|{y92}|{x00}|{y00}|\n")

        with open(neu_txt, "w", encoding="utf-8") as plik1:
            plik1.write("Wyniki obliczeń Geodezyjnych; n, e, u.\n")
            plik1.write("-" * 154 + "\n")
            plik1.write("|                        n                         |                        e                         |                        u                         |\n")
            for n, e, u in zip(N, E, U):
                plik1.write(f"|{n}|{e}|{u}|\n")
            
    def wczytaj_dane(self, nazwa_pliku):
        '''
        Funkcja wczytuje plik z danymi X, Y, Z i tworzy z nich listy posegregowanych X, Y i Z.

        Parametry
        ----------
        nazwa_pliku : str
            Nazwa pliku do wczytania wraz z rozszerzeniem txt.

        Returns
        -------
        X, Y, Z : list
            Listy danych X, Y i Z.
            ilosc_wierszy : int
            Liczba wierszy w pliku.
        '''

        with open(nazwa_pliku, "r") as plik:
            tablica = np.genfromtxt(plik, delimiter=",", dtype='<U20', skip_header=4)
            X = []
            Y = []
            Z = []
            for wiersz in tablica:
                X.append(float(wiersz[0]))
                Y.append(float(wiersz[1]))
                Z.append(float(wiersz[2]))
                ilosc_wierszy = len(X)
                return X, Y, Z, ilosc_wierszy
        
    def wczytanie_zapisanie_pliku(self, dane, output='dms', xyz_txt='Wyniki_transformacji_X_Y_Z_fi_lambda_h_x1992_y1992_x2000_y2000.txt', neu_txt='Wyniki_transformacji_n_e_u.txt'):
        '''
        Wczytuje i zapisuje plik za pomocą jednej funkcji.

        Parameters
        ----------
        dane : str
            Plik z danymi xyz.
            output : str, optional
            Sposób zapisywania współrzędnych f, l [dms, radiany, dec_degree] (domyślnie 'dms').
            xyz_txt : str, optional
            Nazwa pliku wynikowego na xyz, flh, PL1992, PL2000 (domyślnie 'Wyniki_transformacji_X_Y_Z_fi_lambda_h_x1992_y1992_x2000_y2000.txt').
            neu_txt : str, optional
            Nazwa pliku wynikowego dla współrzędnych horyzontalnych (domyślnie 'Wyniki_transformacji_n_e_u.txt').

        Returns
        -------
        Plik txt
        '''

        X, Y, Z, C = Transformacje.wczytaj_dane(self, dane)
        F = []
        L = []
        H = []
        X92 = []
        Y92 = []
        X00 = []
        Y00 = []
        N = []
        E = []
        U = []

        for x, y, z in zip(X, Y, Z):
            f, l, h = Transformacje.xyz2plh(self, x, y, z, output=output)
            if output == "dms":
                F.append(f)
                L.append(l)
            elif output == "radiany":
                f = Transformacje.zamiana_float2string_rad(self, f)
                l = Transformacje.zamiana_float2string_rad(self, l)
            else:
                f = Transformacje.zamiana_float2string_fl(self, f)
                l = Transformacje.zamiana_float2string_fl(self, l)
                F.append(f)
                L.append(l)
            H.append(Transformacje.zamiana_float2string(self, h))
            f,l,h = Transformacje.xyz2plh(self, x, y, z)

            if 13.5 <= l <= 25.5 and 48.9 <= f <= 55.0:
                x92, y92 = Transformacje.flh2PL92(self, f, l)
                X92.append(Transformacje.zamiana_float2string(self, x92))
                Y92.append(Transformacje.zamiana_float2string(self, y92))
                x00, y00 = Transformacje.flh2PL00(self, f, l)
                X00.append(Transformacje.zamiana_float2string(self, x00))
                Y00.append(Transformacje.zamiana_float2string(self, y00))
            else:
                X92.append("         '-'         ")
                Y92.append("         '-'         ")
                X00.append("         '-'         ")
                Y00.append("         '-'         ")

        f1, l1, h1 = Transformacje.xyz2plh(self, X[0], Y[0], Z[0])
        n1,e1,u1 = Transformacje.xyz2neu(self, f1, l1, X[0], Y[0], Z[0], X[-1], Y[-1], Z[-1])
        N.append(n1)
        E.append(e1)
        U.append(u1)

        for i in range(C-1):
            f, l, h = Transformacje.xyz2plh(self, X[i], Y[i], Z[i])
            n, e, u = Transformacje.xyz2neu(self, f, l, X[i], Y[i], Z[i], X[i+1], Y[i+1], Z[i+1])
            N.append(n)
            E.append(e)
            U.append(u)
            print(N)
            
        Transformacje.zapisanie_pliku(self, X, Y, Z, F, L, H, X92, Y92, X00, Y00, N, E, U, xyz_txt, neu_txt )

if __name__ == '__main__':
    
    prze = Transformacje("grs80")
    prze.wczytanie_zapisanie_pliku("wsp_inp.txt")
    
    parser = ArgumentParser()
    parser.add_argument('-m', '--m', type=str, help="Podaj jedną z wskazanych elipsoid: GRS80, WGS84, Krasowski")
    parser.add_argument('-neu', '--neu', type=str, help="Podaj nazwe pliku wynikiowego dla neu z rozszerzeniem txt")
    parser.add_argument('-fa', '--fa', type=float)
    parser.add_argument('-la', '--la', type=float)
    parser.add_argument('-ha', '--ha', type=float)
    parser.add_argument('-fb', '--fb', type=float)
    parser.add_argument('-lb', '--lb', type=float)
    parser.add_argument('-hb', '--hb', type=float)
    args = parser.parse_args()
    
    prze = Transformacje(model = args.m)
    
    
    xa, ya, za = prze.plh2XYZ(args.fa, args.la, args.ha)
    xb, yb, zb = prze.plh2XYZ(args.fb, args.lb, args.hb)
    
    print("")
    print("")
    n, e, u=prze.xyz2neu(args.fa, args.la, xa, ya, za, xb, yb, zb)
    
    n = float(n)
    e = float(e)
    u = float(u)


    print("Elipsida:", args.m)
    print(f"Wyniki_z_flh2neu; n = {n}, e = {e}, u = {u}")
    print("")
    print("")   

    
    
    
'''
    parser = ArgumentParser()
    parser.add_argument('-dane', type=str, help='Wpisz sciezke do pliku z danymi wejsciowymi')
    parser.add_argument('-elip',type=str, help='Wybierz elipsoide sposrod dostepnych: WGS84, GRS80, KRASOWSKI')
    parser.add_argument('-transf', type=str, help='Wybierz transformacje, z ktorej chcesz skorzystac, sposrod dostepnych: XYZ2flh, flh2XYZ, saz2neu, GK2000, GK1992, XYZ2NEU')
    args = parser.parse_args()
    elip = {'WGS84': [6378137.000, 0.00669438002290], 'GRS80': [6378137.000, 0.00669438002290], 'KRASOWSKI': [6378245.000, 0.00669342162296]}
    transf = {'XYZ2FLH': 'XYZ2FLH', 'FLH2XYZ': 'FLH2XYZ','XYZ2NEU': 'XYZ2NEU', 'GK2000': 'GK2000', 'GK1992': 'GK1992'}
'''
'''   
    try:
        wsp = Transformacje(args.elip)
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
   '''     

    



        

        
    
