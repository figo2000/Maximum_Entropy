import numpy as np
import matplotlib.pyplot as plt
import random
from numpy import linalg as LA

d=1000
esp=1/d
a=range(d)
b=range(d)
puntos=[]
for i in range(d):
	 for j in range(d):
	   puntos.append((i*esp, j*esp))


centros=[[0.9048693291718051, 0.6890027262351821],
 [0.9287697231008333, 0.8417041058004717],
 [0.0718836319351367, 0.7765797612573372],
 [0.31713419632009066, 0.17100664982252867],
 [0.3961623437597509, 0.5519265434295347],
 [0.6595792958755342, 0.12658371121037293]]


l=[0.25,0.05,-.08,-.150,0.085,0.035]
lamb=[0,0,0,0,0,0]

def distancia(x,y,z,w):
    return (np.sqrt((x-z)*(x-z)+(y-w)*(y-w)))
    
  


#Retorna un arreglo de 6 arreglos con los puntos de cada una de las regiones
def regiones(lamb):
    regiones=[]
    for i in range(6):
        regiones.append([])

    region_punto=[]
    solo_regiones=[]
    for punto in puntos:
        distancia_cen=[]
        for centro in centros:
            distancia_cen.append(distancia(punto[0],punto[1], centro[0], centro[1])-lamb[centros.index(centro)])
        centroide=min(distancia_cen)
        indice=distancia_cen.index(centroide)
        region_punto.append((punto[0],punto[1],indice))
        solo_regiones.append(indice)
        regiones[indice].append([punto[0],punto[1]])
    return(regiones)   

def solo_reg(lamb):
    regiones=[]
    for i in range(6):
        regiones.append([])

    region_punto=[]
    solo_regiones=[]
    for punto in puntos:
        distancia_cen=[]
        for centro in centros:
            distancia_cen.append(distancia(punto[0],punto[1], centro[0], centro[1])-lamb[centros.index(centro)])
        centroide=min(distancia_cen)
        indice=distancia_cen.index(centroide)
        region_punto.append((punto[0],punto[1],indice))
        solo_regiones.append(indice)
        regiones[indice].append([punto[0],punto[1]])
    return(solo_regiones) 

solo_r=solo_reg(l)  
base_voro=[0,0,0,0,0,0]
solo_r_voro=solo_reg(base_voro)



data1 = np.array(puntos)
x, y = data1.T
plt.scatter(x,y, c=solo_r_voro)
plt.show()
