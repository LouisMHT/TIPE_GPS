# Programme Python du T.I.P.E.
from math import *
import matplotlib.pyplot as plt
import numpy as np
import random as rd

#Définition des 3 Satellites
# in-the-sky.org
#Information du 28 Janvier 2022 à 8h

#Point de référence Lycée Marcel Sembat
latitude = 49.415
longitude = 1.08434
altitude = 0

CoordRef = (latitude, longitude, altitude)
print("Coordonnées du Point de Référence :", CoordRef)

#1er Satellite NAVSTAR 47 ou GPS 20
latitude1 = 50.2
longitude1 = 28.8
altitude1 = 20148

#2ème Satellite NAVSTAR 62 ou GPS 7
latitude2 = 45.1
longitude2 = 27
altitude2 = 20596

#3ème Satellite NAVSTAR 56 ou GPS 2
latitude3 = 20.9
longitude3 = 38.8
altitude3 = 20479

#4ème Satellite Navstar 69 ou GPS 30
latitude4 = 28.2
longitude4 = 5.3
altitude4 = 20332

# Nouveau Programme de conversion en coordonnées cartésiennes
def conversion(R, Lat, Lon):
    x = R * cos(Lat) * cos(Lon)
    y = R * cos(Lat) * sin(Lon)
    z = R * sin(Lat)
    return x, y, z

#Conversion Point de Référence
PointRef = conversion(6371, latitude, longitude)
print("Coordonnées Cartésiennes du point de référence :", PointRef)

#Conversion Sat1
PointSat1 = conversion((6371+altitude1), latitude1, longitude1)
print("Coordonnées Cartésiennes de Sat1 :", PointSat1)

#Conversion Sat2
PointSat2 = conversion((6371+altitude2), latitude2, longitude2)
print("Coordonnées Cartésiennes de Sat2 :", PointSat2)

#Conversion Sat3
PointSat3 = conversion((6371+altitude3), latitude3, longitude3)
print("Coordonnées Cartésiennes de Sat3 :", PointSat3)

#Conversion Sat4
PointSat4 = conversion((6371+altitude4), latitude4, longitude4)
print("Coordonnées Cartésiennes de Sat4 :", PointSat4)

#Fonction de calcul de Rayon des sphères
def rayon_sphère(P, P1):
    x = P[0]
    y = P[1]
    z = P[2]
    x1 = P1[0]
    y1 = P1[1]
    z1 = P1[2]
    d = sqrt(((x1 - x)**2)+((y1-y)**2)+((z1-z)**2))
    return d

#Rayon 1er sphère
Rayonc1 = rayon_sphère(PointRef, PointSat1)
print("Le rayon de la 1ère sphère est :", Rayonc1)

#Rayon 2ème sphère
Rayonc2 = rayon_sphère(PointRef, PointSat2)
print("Le rayon de la 2ème sphère est :", Rayonc2)

#Rayon 3ème sphère
Rayonc3 = rayon_sphère(PointRef, PointSat3)
print("Le rayon de la 3ème sphère est :", Rayonc3)

#Rayon 4ème sphère
Rayonc4 = rayon_sphère(PointRef, PointSat4)
print("Le rayon de la 3ème sphère est :", Rayonc3)

#Représentation des points et sphères en 3D
fig = plt.figure()
ax = fig.add_subplot(projection='3d')

#ax.scatter([PointRef[0]], [PointRef[1]], [PointRef[2]], color='g')
ax.scatter([PointSat1[0]], [PointSat1[1]], [PointSat1[2]], color='r')
ax.scatter([PointSat2[0]], [PointSat2[1]], [PointSat2[2]], color='r')
ax.scatter([PointSat3[0]], [PointSat3[1]], [PointSat3[2]], color='r')
ax.scatter([PointSat4[0]], [PointSat4[1]], [PointSat4[2]], color='r')

#sphère1
u, v = np.mgrid[0:2*np.pi:100j, 0:np.pi:100j]
x = Rayonc1*np.cos(u)*np.sin(v) + PointSat1[0]
y = Rayonc1*np.sin(u)*np.sin(v) + PointSat1[1]
z = Rayonc1*np.cos(v) + PointSat1[2]
ax.plot_wireframe(x, y, z, color='b', linewidth=0.1)

#sphère2
u, v = np.mgrid[0:2*np.pi:100j, 0:np.pi:100j]
x = Rayonc2*np.cos(u)*np.sin(v) + PointSat2[0]
y = Rayonc2*np.sin(u)*np.sin(v) + PointSat2[1]
z = Rayonc2*np.cos(v) + PointSat2[2]
ax.plot_wireframe(x, y, z, color='r', linewidth=0.1)

#sphère3
u, v = np.mgrid[0:2*np.pi:100j, 0:np.pi:100j]
x = Rayonc3*np.cos(u)*np.sin(v) + PointSat3[0]
y = Rayonc3*np.sin(u)*np.sin(v) + PointSat3[1]
z = Rayonc3*np.cos(v) + PointSat3[2]
ax.plot_wireframe(x, y, z, color='k', linewidth=0.1)


#Segment entre 2 points
x1_valeurs = [PointSat1[0], PointSat2[0]]
y1_valeurs = [PointSat1[1], PointSat2[1]]
z1_valeurs = [PointSat1[2], PointSat2[2]]

plt.plot(x1_valeurs, y1_valeurs, z1_valeurs, color="r")

x2_valeurs = [PointSat2[0], PointSat3[0]]
y2_valeurs = [PointSat2[1], PointSat3[1]]
z2_valeurs = [PointSat2[2], PointSat3[2]]

plt.plot(x2_valeurs, y2_valeurs, z2_valeurs, color="r")

#Segment de vérification
x3_valeurs = [PointRef[0], PointSat4[0]]
y3_valeurs = [PointRef[1], PointSat4[1]]
z3_valeurs = [PointRef[2], PointSat4[2]]

plt.plot(x3_valeurs, y3_valeurs, z3_valeurs, color="g")

#Intersection entre la sphère 1 et la sphère 2
#Utilisation de la formule de Euclide
#Distance entre les deux centres

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

Rayon1 = rayon_cercle(PointSat1, PointSat2, Rayonc1, Rayonc2)
print("Le rayon du cercle entre la 1er et la 2ème sphère est :", Rayon1[0])
PC1 = Rayon1[1]
print("Coordonnées Cartésiennes du centre du 1er cercle :", PC1)
ax.scatter([PC1[0]], [PC1[1]], [PC1[2]], color='b')

Rayon2 = rayon_cercle(PointSat2, PointSat3, Rayonc2, Rayonc3)
print("Le rayon du cercle entre la 2ème et la 3ème Sphère est :", Rayon2[0])
PC2 = Rayon2[1]
print("Coordonnées Cartésiennes du centre du 2ème cercle :", PC2)
ax.scatter([PC2[0]], [PC2[1]], [PC2[2]], color='b')

#plan
#point PC1
normal = ((PointSat1[0] - PC1[0])/sqrt((PointSat1[0] - PC1[0])**2+(PointSat1[1] - PC1[1])**2+(PointSat1[2] - PC1[2])**2), (PointSat1[1] - PC1[1])/sqrt((PointSat1[0] - PC1[0])**2+(PointSat1[1] - PC1[1])**2+(PointSat1[2] - PC1[2])**2), (PointSat1[2] - PC1[2])/sqrt((PointSat1[0] - PC1[0])**2+(PointSat1[1] - PC1[1])**2+(PointSat1[2] - PC1[2])**2))
print(normal)
d = (-normal[0]*PC1[0]-normal[1]*PC1[1]-normal[2]*PC1[2])
print(d)

#point PC2
normal2 = ((PointSat3[0] - PC2[0])/sqrt((PointSat3[0] - PC2[0])**2+(PointSat3[1] - PC2[1])**2+(PointSat3[2] - PC2[2])**2), (PointSat3[1] - PC2[1])/sqrt((PointSat3[0] - PC2[0])**2+(PointSat3[1] - PC2[1])**2+(PointSat3[2] - PC2[2])**2), (PointSat3[2] - PC2[2])/sqrt((PointSat3[0] - PC2[0])**2+(PointSat3[1] - PC2[1])**2+(PointSat3[2] - PC2[2])**2))
print(normal2)
d2 = (-normal2[0]*PC2[0]-normal2[1]*PC2[1]-normal2[2]*PC2[2])
print(d2)

#Représentation graphique pour les cercles
u1 = (normal[1]/sqrt((normal[0]**2)+(normal[1]**2)), -normal[0]/sqrt((normal[0]**2)+(normal[1]**2)), 0)
print(u1)
u2 = (normal2[1]/sqrt((normal2[0]**2)+(normal2[1]**2)), -normal2[0]/sqrt((normal2[0]**2)+(normal2[1]**2)), 0)
print(u2)

def produit_vect(n, u):
    x = n[1]*u[2] - n[2]*u[1]
    y = n[2]*u[0] - n[0]*u[2]
    z = n[0]*u[1] - n[1]*u[0]
    return (x, y, z)

v1 = produit_vect(normal, u1)
v2 = produit_vect(normal2, u2)

#cercle1
phi = [i for i in np.linspace(0, 2*np.pi, 100)]
a1 = PC1[0] + Rayon1[0]*np.cos(phi)*u1[0] + Rayon1[0]*np.sin(phi)*v1[0]
b1 = PC1[1] + Rayon1[0]*np.cos(phi)*u1[1] + Rayon1[0]*np.sin(phi)*v1[1]
c1 = PC1[2] + Rayon1[0]*np.cos(phi)*u1[2] + Rayon1[0]*np.sin(phi)*v1[2]
ax.plot(a1, b1, c1, color='b')

#cercle2
phi = [i for i in np.linspace(0, 2*np.pi, 100)]
a2 = PC2[0] + Rayon2[0]*np.cos(phi)*u2[0] + Rayon2[0]*np.sin(phi)*v2[0]
b2 = PC2[1] + Rayon2[0]*np.cos(phi)*u2[1] + Rayon2[0]*np.sin(phi)*v2[1]
c2 = PC2[2] + Rayon2[0]*np.cos(phi)*u2[2] + Rayon2[0]*np.sin(phi)*v2[2]
ax.plot(a2, b2, c2, color='b')

#droite sous forme paramétrique
#x = t

zt = (-normal[0])/normal[2]
zy = (-normal[1])/normal[2]
zo = (-d)/normal[2]

it = normal2[2]*zt
iy = normal2[2]*zy
io = normal2[2]*zo

nt = normal2[0]+it
ny = normal2[1]+iy
no = d2+io

#y
yt = (-nt)/ny
yo = (-no)/ny
print("Droite sous forme paramétrique :")
print("x = t")
print("y = ", yt, "t +", "(", yo, ")")

nzt = (zy*yt)+zt
nzo = (zy*yo)+zo

print("z = ", nzt, "t +", "(", nzo, ")")

#Verif final
print("t =", PointRef[0])
#PointRef[0] = yt*t+yo
ty = (PointRef[1] - yo)/yt
print("t =", ty)
#PointRef[2] = nzt*t+nzo
tz = (PointRef[2] - nzo)/nzt
print("t =", tz)

if round(PointRef[0]) == round(ty) == round(tz):
    print("Vérification correcte !")
else:
    print("nope")

vect = (1, yt, nzt)
point = (0, yo, nzo)

print("équation de la sphère issu de l'intersection des deux sphères :")
print("(x-",PC1[0],")^2+(y-",PC1[1],")^2+(z-",PC1[2],")^2=",Rayon1[0],"^2")

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

Pi1 = delta(interspheredroite1(PC1, Rayon1[0], vect, point))
Pi2 = delta(interspheredroite1(PC2, Rayon2[0], vect, point))

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

P1f = (interspheredroite2(vect, point, Pi1))[0]
P2f = (interspheredroite2(vect, point, Pi1))[1]

P1f2 = (interspheredroite2(vect, point, Pi2))[0]
P2f2 = (interspheredroite2(vect, point, Pi2))[1]

print(P1f)
print(P2f)
print(P1f2)
print(P2f2)

if ((round(P1f[0]) == round(P1f2[0]) and round(P1f[1]) == round(P1f2[1]) and round(P1f[1]) == round(P1f2[1]))\
        and (round(P2f[0]) == round(P2f2[0]) and round(P2f[1]) == round(P2f2[1]) and round(P2f[1]) == round(P2f2[1])))\
        or ((round(P1f[0]) == round(P2f2[0]) and round(P1f[1]) == round(P2f2[1]) and round(P1f[1]) == round(P2f2[1]))\
        and (round(P2f[0]) == round(P1f2[0]) and round(P2f[1]) == round(P1f2[1]) and round(P2f[1]) == round(P1f2[1]))):
    print('La vérification des cercles est ok')
else:
    print('ERROR 404')
    exit()

ax.scatter([P1f[0]], [P1f[1]], [P1f[2]], color='magenta')
ax.scatter([P2f[0]], [P2f[1]], [P2f[2]], color='magenta')

def Point_final(P):
    if round(rayon_sphère(P, PointSat4)) == round(rayon_sphère(PointRef, PointSat4)):
        return True
    else:
        return False

if Point_final(P1f) == True:
    print('Les coordonnées cartésiennes du point final sont :', P1f)
    Pf = P1f
elif Point_final(P2f) == True:
    print('Les coordonnées cartésiennes du point final sont :', P2f)
    Pf = P2f
else:
    print('ERROR 404')

if round(PointRef[0]) == round(Pf[0]) and round(PointRef[1]) == round(Pf[1]) and round(PointRef[2]) == round(Pf[2]):
    print('Le point final est bien le même point que le point de référence !')
else:
    print('ERROR 404')

ax.set_box_aspect([60000, 60000, 60000])
plt.axis("off")
plt.show()
