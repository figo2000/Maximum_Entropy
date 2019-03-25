import numpy as np
import matplotlib.pyplot as plt
import random
from numpy import linalg as LA

d=100
esp=0.01
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
    
  




#Retorna la region a la que pertenece el punto dado con el lambda dado
def region_punto(punto,la):
    distancia_cen=[]
    for centro in centros:
        resta=centros.index(centro)
        distancia_cen.append(distancia(punto[0],punto[1], centro[0], centro[1])-la[resta])
        centroide=min(distancia_cen)
        indice=distancia_cen.index(centroide)
    return indice
        
#Evalua la funcion objetivo de la entropia maxima en el punto (x,y) dado lambda, u, v
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

# halla int_{R}e^{-1-u*min_{i}{||x-x_i||-lambda_i}+v} dados u,v y lambda [LA FUNCION A MINIMIZAR]
def pa_min_u_v(u,v,lam,t):
    integral=0
    for punto in puntos:
        region=region_punto(punto,lam)
        la=lam[region]
        mini=distancia(punto[0],punto[1], centros[region][0], centros[region][1])-la
        integral+=np.exp(-1-u*mini-v)*esp*esp
    return(integral+u*t+v)   




#MINIMIZACION U,V: se hallan los u,v optimos y se tiene la funcion F(u,v,\lambda) y calculamos el valor de la entropia de esta funcion.(ESTO SE HACE EN 15 MINUTOS)


#Calcula el gradiente de la funcion a minimizar
def gradiente(u, v, lam, t):
    integral=0
    for punto in puntos:
        region=region_punto(punto,lam)
        la=lam[region]
        mini=distancia(punto[0],punto[1], centros[region][0], centros[region][1])-la
        integral+=np.exp(-1-u*mini-v)*mini*esp*esp
    return(t-integral, 1-integral_u_v_entro(u,v,lam)) 

#Calcula la longitud del paso en el descenso del gradiente con BACKTRACKING
def function_backtracking(p,alph,bet,lam,t):
    t = 1
    while (integral_u_v_entro((np.array(p) -t*np.array(gradiente(p[0],p[1], lam, t)))[0],(np.array(p) -t*np.array(gradiente(p[0],p[1], lam, t)))[1], lam) > integral_u_v_entro(p[0],p[1],lam) + alph*t*LA.norm(gradiente(p[0],p[1], lam, t))**2): 
        t = bet*t
    return t

#Calcula el siguiente paso del descenso con la restriccion u>=0
def desc_back_step(p, alph, bet, lam, t):
    rta = np.array(p)-function_backtracking(p,alph,bet,lam,t)*np.array(gradiente(p[0],p[1],lam,t))
    if rta[0]>=0:
        return (rta)
    else: 
        rta[0]=0
        return (rta)
#halla un minimo local empezando en un punto p. con 10 pasos  backtracking
def descenso_gr_back(p,alph, bet, lam, t):
    p_actual=p
    #gradiente_act=1
    conteo=0
    while conteo <5  :
        conteo+=1
        p_actual = desc_back_step(p_actual, alph, bet, lam, t)
        #gradiente_act =LA.norm(gradiente(p_actual[0],p_actual[1],l , 0.001))
    return(p_actual)             


#Halla el minimo empezando con una grilla
def minimos_u_v(lamb,t):
    puntos_s=[[1,-8],[1,-3],[1,2],[1,7],[4,-8],[4,-3],[4,2],[4,7], [7,-8],[7,-3],[7,2],[7,7], [10,-8],[10,-3],[10,2],[10,7]]
    minimos_locales=[]
    minimizadores_local=[]
    gradientes=[]
    for i in puntos_s:
        local=descenso_gr_back([i[0],  i[1]],0.5,0.5,lamb,t)
        minimizadores_local.append(local)        
    gradientes=[]
    minimos_locales=[]
    for i in minimizadores_local:
        gradientes.append(LA.norm(gradiente(i[0],i[1],lamb, t)))
        minimos_locales.append(pa_min_u_v(i[0],i[1],lamb,t))
    minimo=min(minimos_locales)
    indice=minimos_locales.index(minimo)
    minimizador=minimizadores_local[indice]
    return(minimo, minimizador)       


print(minimos_u_v(l,0.))



def g(u,v,l):
	rta=[0,0,0,0,0,0]
	for punto in puntos:
		region=region_punto(punto,l)
		rta[region]-=fun_puntual_entrop(l,u,v,punto[0],punto[1])*esp*esp
        
	return(rta)


	
v=np.dot(l, [-0.094991399392880496, -0.0378595075341773, -0.0044523034955456864, -0.011950884202663973, -0.072596007029297652, -0.030531496828136218])
mat=[[1,1,1,1,1,1],[-1,-1,-1,-1,-1,-1], [-0.094991399392880496, -0.0378595075341773, -0.0044523034955456864, -0.011950884202663973, -0.072596007029297652, -0.030531496828136218]]
b=[0,0, v]

def escoge_lam(matriz, vector):
    c = [1,0,1,0,0,3]
    A = matriz
    b = vector
    x0_bounds = (None, np.sqrt(2))
    x1_bounds = (None, np.sqrt(2))
    x2_bounds = (None, np.sqrt(2))
    x3_bounds = (None, np.sqrt(2))
    x4_bounds = (None, np.sqrt(2))
    x5_bounds = (None, np.sqrt(2))
    from scipy.optimize import linprog
    res = linprog(c, A_ub=A, b_ub=b, bounds=(x0_bounds, x1_bounds, x2_bounds, x3_bounds, x4_bounds, x5_bounds))
    return(res.x)

















ll=[[0.44232470625375619, 0.82280622459131059],
 [0.35513191769163754, 0.66307976190303786],
 [0.096687945406576659, 0.96880545347250746],
 [0.16432801880909642, 0.060775949774398041],
 [0.78029051537435767, 0.40552617150737336],
 [0.0097722975376862697, 0.60431616552166056]]










                     
