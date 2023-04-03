from math import *
import random as rd
from itertools import combinations

#Coordonnées GPS
PR = (49.415, 1.08434, 0)
S1 = (50.2, 28.8, 20148)
S2 = (45.1, 27, 20596)
S3 = (20.9, 38.8, 20479)
S4 = (28.2, 5.3, 20332)
S5 = (40, 63.2, 20136)

#Données de mesure
P1map = (40.7600000, -73.9840000, 6371)
P2map = (40.69001, -75.16777, 6371)
d12 = 100085.7 #m

def conversion(P):
    Lat, Lon, R = P
    R = R + 6371
    x = R * cos(Lat) * cos(Lon)
    y = R * cos(Lat) * sin(Lon)
    z = R * sin(Lat)
    return (x, y, z)

def rayon_sphère(P, P1):
    x = P[0]
    y = P[1]
    z = P[2]
    x1 = P1[0]
    y1 = P1[1]
    z1 = P1[2]
    d = sqrt(((x1 - x)**2)+((y1-y)**2)+((z1-z)**2))
    return d

def rayon_cercle(P1, P2, R1, R2):
    O1O2 = rayon_sphère(P1, P2)
    beta = acos(((O1O2 ** 2) + (R2 ** 2) - (R1 ** 2)) / (2 * O1O2 * R2))
    d2 = cos(beta) * R2
    h = tan(beta) * d2
    d1 = O1O2 - d2
    x1 = P1[0]
    y1 = P1[1]
    z1 = P1[2]
    x2 = P2[0]
    y2 = P2[1]
    z2 = P2[2]
    x = ((d1*x2)+(d2*x1))/(d1+d2)
    y = ((d1*y2)+(d2*y1))/(d1+d2)
    z = ((d1*z2)+(d2*z1))/(d1+d2)
    P = (x, y, z)
    return (h, P)

def interspheredroite1(P, R, vect, point):
    Sx = P[0]
    Sy = P[1]
    Sz = P[2]
    Sr = R**2
    DVa = vect[0]
    DVb = vect[1]
    DVc = vect[2]
    Da = point[0]
    Db = point[1]
    Dc = point[2]
    A1 = Da - Sx
    A2 = Db - Sy
    A3 = Dc - Sz
    B1 = (A1**2) + (A2**2) + (A3**2) - Sr
    B2 = (2*A1*DVa) + (2*A2*DVb) + (2*A3*DVc)
    B3 = (DVa**2) + (DVb**2) + (DVc**2)
    a = B3
    b = B2
    c = B1
    L = (a, b, c)
    return L

def delta(l:list):
    a = l[0]
    b = l[1]
    c = l[2]
    d = (b**2)-4*a*c
    if d > 0:
        r1 = (-b-sqrt(d))/(2*a)
        r2 = (-b+sqrt(d))/(2*a)
        return (r1, r2)
    elif d == 0:
        r0 = -b/(2*a)
        return r0
    else:
        return None

def interspheredroite2(vect, point, P):
    r1 = P[0]
    r2 = P[1]
    DVa = vect[0]
    DVb = vect[1]
    DVc = vect[2]
    Da = point[0]
    Db = point[1]
    Dc = point[2]
    V111 = DVa*r1
    V112 = DVb*r1
    V113 = DVc*r1
    V121 = V111 + Da
    V122 = V112 + Db
    V123 = V113 + Dc
    V211 = DVa*r2
    V212 = DVb*r2
    V213 = DVc*r2
    V221 = V211 + Da
    V222 = V212 + Db
    V223 = V213 + Dc
    P1 = (V121, V122, V123)
    P2 = (V221, V222, V223)
    return (P1, P2)

def Point_final(d, d1):
    if d > d1:
        return True
    elif d < d1:
        return False
    else:
        exit()

def one(PR, S1, S2, S3, S4, S5, P1, P2, d12):
    Lerreur = []
    Lpoint = [conversion(S1), conversion(S2), conversion(S3), conversion(S4), conversion(S5), conversion(PR)]
    SAT1 = conversion(S1)
    SAT2 = conversion(S2)
    SAT3 = conversion(S3)
    SAT4 = conversion(S4)
    SAT5 = conversion(S5)
    Pointfinalx = []
    Pointfinaly = []
    Pointfinalz = []
    Proba = []
    Pmesure = [conversion(P1), conversion(P2)]
    distP1MP2M = rayon_sphère(Pmesure[0], Pmesure[1])
    dinter = distP1MP2M / d12
    for i in range(5):
        eo = rd.uniform(2, 50)  # erreur orbite
        ei = rd.uniform(0.5, 100)  # erreur ionosphère
        et = rd.uniform(0.01, 0.5)  # erreur troposphère
        e = eo + ei + et
        e = e * dinter
        Lerreur.append(e)
    Rayon = [rayon_sphère(Lpoint[5], Lpoint[0]) + Lerreur[0], rayon_sphère(Lpoint[5], Lpoint[1]) + Lerreur[1],
             rayon_sphère(Lpoint[5], Lpoint[2]) + Lerreur[2], rayon_sphère(Lpoint[5], Lpoint[3]) + Lerreur[3],
             rayon_sphère(Lpoint[5], Lpoint[4])]
    for i in range(len(Lpoint)-1):
        r = rayon_sphère(Lpoint[-1], Lpoint[i])
        Rayon.append(r)
    A = [SAT1, SAT2, SAT3, SAT4, SAT5]
    temp = combinations(A, 3)
    for i in list(temp):
        Proba.append(i)
    #print(Proba)
    #print(len(Proba))
    n = 0
    while n < len(Proba):
        x = Proba[n]
        Rayonbis = []
        if SAT1 in x:
            Rayonbis.append(Rayon[0])
        if SAT2 in x:
            Rayonbis.append(Rayon[1])
        if SAT3 in x:
            Rayonbis.append(Rayon[2])
        if SAT4 in x:
            Rayonbis.append(Rayon[3])
        if SAT5 in x:
            Rayonbis.append(Rayon[4])
        #print(x)
        #print(Rayonbis)
        Rayon1 = rayon_cercle(x[0], x[1], Rayonbis[0], Rayonbis[1])
        PC1 = Rayon1[1]
        Rayon2 = rayon_cercle(x[1], x[2], Rayonbis[1], Rayonbis[2])
        PC2 = Rayon2[1]
        normal = ((x[0][0] - PC1[0]) / sqrt((x[0][0] - PC1[0]) ** 2 + (x[0][1] - PC1[1]) ** 2 + (x[0][2] - PC1[2]) ** 2),
                  (x[0][1] - PC1[1]) / sqrt((x[0][0] - PC1[0]) ** 2 + (x[0][1] - PC1[1]) ** 2 + (x[0][2] - PC1[2]) ** 2),
                  (x[0][2] - PC1[2]) / sqrt((x[0][0] - PC1[0]) ** 2 + (x[0][1] - PC1[1]) ** 2 + (x[0][2] - PC1[2]) ** 2))
        d = (-normal[0] * PC1[0] - normal[1] * PC1[1] - normal[2] * PC1[2])
        normal2 = ((x[2][0] - PC2[0]) / sqrt((x[2][0] - PC2[0]) ** 2 + (x[2][1] - PC2[1]) ** 2 + (x[2][2] - PC2[2]) ** 2),
                   (x[2][1] - PC2[1]) / sqrt((x[2][0] - PC2[0]) ** 2 + (x[2][1] - PC2[1]) ** 2 + (x[2][2] - PC2[2]) ** 2),
                   (x[2][2] - PC2[2]) / sqrt((x[2][0] - PC2[0]) ** 2 + (x[2][1] - PC2[1]) ** 2 + (x[2][2] - PC2[2]) ** 2))
        d2 = (-normal2[0] * PC2[0] - normal2[1] * PC2[1] - normal2[2] * PC2[2])

        # Droite
        zt = (-normal[0]) / normal[2]
        zy = (-normal[1]) / normal[2]
        zo = (-d) / normal[2]

        it = normal2[2] * zt
        iy = normal2[2] * zy
        io = normal2[2] * zo

        nt = normal2[0] + it
        ny = normal2[1] + iy
        no = d2 + io

        # y
        yt = (-nt) / ny
        yo = (-no) / ny

        nzt = (zy * yt) + zt
        nzo = (zy * yo) + zo

        vect = (1, yt, nzt)
        point = (0, yo, nzo)

        Pi1 = delta(interspheredroite1(PC1, Rayon1[0], vect, point))
        Pi2 = delta(interspheredroite1(PC2, Rayon2[0], vect, point))

        P1f = (interspheredroite2(vect, point, Pi1))[0]
        P2f = (interspheredroite2(vect, point, Pi1))[1]

        P1f2 = (interspheredroite2(vect, point, Pi2))[0]
        P2f2 = (interspheredroite2(vect, point, Pi2))[1]

        if ((round(P1f[0]) == round(P1f2[0]) and round(P1f[1]) == round(P1f2[1]) and round(P1f[1]) == round(P1f2[1]))
            and (round(P2f[0]) == round(P2f2[0]) and round(P2f[1]) == round(P2f2[1]) and round(P2f[1]) == round(P2f2[1]))) \
                or ((round(P1f[0]) == round(P2f2[0]) and round(P1f[1]) == round(P2f2[1]) and round(P1f[1]) == round(P2f2[1]))
                    and (round(P2f[0]) == round(P1f2[0]) and round(P2f[1]) == round(P1f2[1]) and round(P2f[1]) == round(P1f2[1]))):
            None
        else:
            exit()

        doups1 = rayon_sphère(P1f, Lpoint[5])
        doups2 = rayon_sphère(P2f, Lpoint[5])

        if Point_final(doups1, doups2) == False:
            # print('Les coordonnées cartésiennes du point final sont :', P1f)
            Pf = P1f
        elif Point_final(doups1, doups2) == True:
            # print('Les coordonnées cartésiennes du point final sont :', P2f)
            Pf = P2f

        Pointfinalx.append(Pf[0])
        Pointfinaly.append(Pf[1])
        Pointfinalz.append(Pf[2])

        n = n + 1

    Pfx = sum(Pointfinalx)/len(Pointfinalx)
    Pfy = sum(Pointfinaly)/len(Pointfinaly)
    Pfz = sum(Pointfinalz)/len(Pointfinalz)
    #print(Pfx, Pfy, Pfz)
    Pointfinal = (Pfx, Pfy, Pfz)
    #print(Pointfinal)
    dprecision = rayon_sphère(Pointfinal, Lpoint[5])
    dfinal = 1 / dinter
    dprecisionmetre = dfinal * dprecision
    return dprecisionmetre

#print(one(PR, S1, S2, S3, S4, S5, P1map, P2map, d12))

L = []
for i in range(10000):
    x = (one(PR, S1, S2, S3, S4, S5, P1map, P2map, d12))
    L.append(x)
#print(L)

m = max(L)

print("Le maximum de l'écart de précision GPS qu'on peut avoir lorsqu'on exécute le programme 10000 fois est d'environ", round(m), "mètres.")
