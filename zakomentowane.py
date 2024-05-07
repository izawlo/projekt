# -*- coding: utf-8 -*-
"""
Created on Tue May  7 11:23:25 2024

@author: alawo
"""

'''
if __name__ == "__main__":
    # utworzenie obiektu
    geo = Transformacje(model = "wgs84")
    print(sys.argv)
    # dane XYZ geocentryczne
    X = 3664940.500; Y = 1409153.590; Z = 5009571.170
    phi, lam, h = geo.xyz2plh(X, Y, Z)
    print(phi, lam, h)
    inp_file_path=sys.argv[-1]
    # phi, lam, h = geo.xyz2plh2(X, Y, Z)
    # print(phi, lam, h)
    inp_file_path = sys.argv[-1]
    
    # if '--xyz2plh' in sys.argv and '--plh2xyz' in sys.argv:
    #     print('możesz podać tylko jednej flagę')
    # elif '--xyz2plh' in sys.argv:
    #     coords_plh=[]
    #     with open(inp_file_path) as f:
            
    






    if '--xyz2plh' in sys.argv[-1]
    coords_plh = []
    with open ('wsp_inp.txt') as f:
    lines = f.readlines()
    lines = lines[4:]
    for line in lines:
    line = line.strip()
    x_str, y_str, z_str = line.split(',')
    x, y, z =(float(x_str), float(y_str), float(z_str))
    p, l, h = geo.xyz2plh(x,y,z)
    coords_plh.append([p,l,h])

    with open('result_xyz2plh.txt', 'w') as f:
    f.write('phi[deg],lam[phi],h[m]\n')
    for coords in coords_plh:
    line = ','.join([str(coord) for coord in coords])
    f.write(line + '\n')
    
    elif '--plh2xyz' in sys.argv:
        coords_xyz=[]
        with open ('inp_file_path.txt') as f:
        lines = f.readlines()
        lines = lines[1:]
            for line in lines:
            line = line.strip()
            phi_str, lam_str, h_str = line.split(',')
            phi, lam, h =(float(phi_str), float(lam_str), float(h_str))
            x, y, z = geo.xyz2plh(phi,lam,h)
            coords_plh.append([x,y,z])

    with open('result_plh2xyz.txt', 'w') as f:
    f.write('x[m],y[m],z[m]\n')
        for coords in coords_plh:
        line = ','.join([str(coord) for coord in coords])
        f.write(line + '\n')
'''
    
"""Tranformacja współrzędnych geocentryczny do współrzędnych topocentrycznych
def Rneu(self,f,l):
    R = np.array([[-np.sin(f)*np.cos(l), -np.sin(l), np.cos(f)*np.cos(l)],
                  [-np.sin(f)*np.sin(l), np.cos(l), np.cos(f)*np.sin(l)],
                  [np.cos(f), 0., np.sin(f)]])
    return(R)
"""
