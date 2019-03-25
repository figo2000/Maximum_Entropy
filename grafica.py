import numpy as np
import matplotlib.pyplot as plt
import random
from numpy import linalg as LA
from mpl_toolkits.mplot3d import Axes3D

#Los centros de la distribucion empirica
centros=[[0.44232470625375619, 0.82280622459131059],
 [0.35513191769163754, 0.66307976190303786],
 [0.096687945406576659, 0.96880545347250746],
 [0.16432801880909642, 0.060775949774398041],
 [0.78029051537435767, 0.40552617150737336],
 [0.0097722975376862697, 0.60431616552166056]]


centros=[[0.9048693291718051, 0.6890027262351821],
 [0.9287697231008333, 0.8417041058004717],
 [0.0718836319351367, 0.7765797612573372],
 [0.31713419632009066, 0.17100664982252867],
 [0.3961623437597509, 0.5519265434295347],
 [0.6595792958755342, 0.12658371121037293]]


l=[0.25,0.05,-.08,-.150,0.085,0.035]
lamb=[0,0,0,0,0,0]

d=100
esp=0.01
a=range(d)
b=range(d)
puntos=[]
for i in range(d):
	 for j in range(d):
	   puntos.append((i*esp, j*esp))
#Distancia entre dos puntos
def distancia(x,y,z,w):
    return (np.sqrt((x-z)*(x-z)+(y-w)*(y-w)))
    

def gradiente(u, v, lam, t):
    integral=0
    for punto in puntos:
        region=region_punto(punto,lam)
        la=lam[region]
        mini=distancia(punto[0],punto[1], centros[region][0], centros[region][1])-la
        integral+=np.exp(-1-u*mini-v)*mini*esp*esp
    return(t-integral, 1-integral_u_v_entro(u,v,lam)) 

#Retorna la region a la que pertenece el punto dado
def region_punto(punto,la):
    distancia_cen=[]
    for centro in centros:
        resta=centros.index(centro)
        distancia_cen.append(distancia(punto[0],punto[1], centro[0], centro[1])-la[resta])
        centroide=min(distancia_cen)
        indice=distancia_cen.index(centroide)
    return indice

#dada una distribucion dada por (lamb, u,v) la evalua en (x,y) (Esto es para graficar)
def fun_puntual_entrop(lamb, u, v, x, y):
    region=region_punto([x,y],lamb)
    l=lamb[region]
    mini=distancia(x,y, centros[region][0], centros[region][1])
    expon=-1-(u*(mini-l)+v)
    return (np.exp(expon))
# dados u,v, lamb calcula \iint_{R}e^{-1-u*min{...}-v}dA    
def integral_u_v_entro(u,v,lam):
    integral=0
    for punto in puntos:
        region=region_punto(punto,lam)
        la=lam[region]
        mini=distancia(punto[0],punto[1], centros[region][0], centros[region][1])-la
        integral+=np.exp(-1-u*mini-v)*esp*esp
    return(integral)     


minimizador=  4.14750228, -1.5251104
minimizador=[6.97841773, -1.41446287]
l=[0.25,0.05,-.08,-.150,0.085,0.035]
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
x = y = np.arange(0, 1.0, 0.01)
X, Y = np.meshgrid(x, y)
zs = np.array([fun_puntual_entrop(l,    minimizador[0],  minimizador[1],x,y) for x,y in zip(np.ravel(X), np.ravel(Y))])
Z = zs.reshape(X.shape)

ax.plot_surface(X, Y, Z)

ax.set_xlabel('X Label')
ax.set_ylabel('Y Label')
ax.set_zlabel('Z Label')

plt.show()

print(gradiente(minimizador[0],minimizador[1],lamb, 0.1))


#Ejemplo 5 pasos  [10.16837388,  -2.63724589]