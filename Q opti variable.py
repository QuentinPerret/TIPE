import matplotlib.pyplot as plt
import  numpy as np
from scipy.integrate import odeint

#Initialisation des constantes du problèmes
duree_s = 20. #durée de l'expérience (calcul de x et y)
duree_t = 5000. #durée de l'expérience (calcul de z)
dt = 1.  #durée d'un interval de temps entre 2 mesures
m = 1. #constante de la loi de création de la biomasse
Vr = 1. #volume du bioréacteur
Vl = 100. # volume del'étendu d'eau 

def model(var,t,Q):	
	x = var[0]
	y = var[1]
	z = var[2]	
	dxdt =  m*x*y - Q/Vr*x
	dydt = Q/Vr*(z - y) - m*x*y
	dzdt =  Q * (y - z) / Vl	
	return [dxdt , dydt , dzdt]

#Conditions Initiales
x0 = 0.1 #valeur initiale de x
y0 = 1. #valeur initiale de y
z0 = 1. #valeur initiale de z
var0= [x0,y0,z0]

#Définition du temps 
n = int(duree_t//dt)
t = np.linspace(0,duree_t,n)

x = np.zeros(n)
y = np.zeros(n)
z = np.zeros(n)


def Qconstant(var0,obj,esp):	
	""" Entrée : var0 constantes initiales du problème , obj objectif de dépollution , eps précision de la valeur du débit 
		Sortie : débit contant optimal pour atteindre obj le plus rapidement possible Q est précis à eps près et l'indice dans le tableau t du temps mis pour atteindre l'objectif obj	"""	
	q = np.linspace(0,1,int(1/esp)) #on définit toutes les valeurs que peut prendre le débit q
	Q_min , m = 0 , n #on crée les variables pour contenir les résultats à renvoyer	
	for Q in q: #on parcours tous les débits possibles	
		Z = odeint(model,var0,t,args=(Q,))[:,2] #pour chaque débit on calcul les valeurs de z à chaque instant 
		j = 0 		
		while j < n and Z[j] > obj: #tant que la valeur de Z est plus grande que l'objectif de dépollution j += 1 			
			j += 1 		
		#Ici on a j l'indice du temps dans t du temùps nécéssaire pour atteindre l'objectif obj		
		if j < m: #si j est plus petit que mù soit le temps mis par Q inférieur à celui par Q_min pour atteindre le plus rapidement l'obkectif de dépollution 			
			Q_min = Q
			m = j #alors on change les valeurs de Q_mùin et mù
	return  Q_min , m #on renvoie le débit contant optimal pour atteindre obj le plus rapidement possible Q est précis à eps près et l'indice dans le tableau t du temps mis pour atteindre l'objectif obj
	
	
def Q_variable_opti(obj,eps):	
	'''renvoie un tableau res avec toutes les valeurs de Q à la précision eps 
	   pour chaque intervalle de temps dt afin d'atteindre l'objectif de dépollution obj le plus rapidement possible  '''	
	N = int(1/eps)
	t = np.linspace(0,duree_t,n)
	z = z0
	var_i = var0	
	res = np.zeros([N,2])	
	palier = np.linspace(z0,obj,N)	
	i = 0	
	while  i < N:		
		if palier[i] > var_i[2]:			
			i+= 1			
		else:			
			Q_opti , T_I_opti = Qconstant(var_i,palier[i],eps)
			Q = Q_opti			
			var = odeint( model,var_i,t,args=(Q,))
			var_i = var[T_I_opti]			
			t = t[T_I_opti:]
			res[i,0] , res[i,1] = Q_opti , t[0]			
			i += 1			
	return res

def Q(Q_variable,temps):
	''' Renvois le débit Q optimal à l'instant t pour un débit variable '''
	i = 0	
	while i < len(Q_variable)-1 and temps > Q_variable[i,1]:		
		i += 1
	return Q_variable[i,0] 

def dessin_Z_opti(obj,eps):	
	'''Fonction dessinant z en fonction de t pour Q(t) donné'''	
	Q_variable = Q_variable_opti(obj,eps)	
	var0  = [0.1,1.,1.]	
	x[0] = var0[0]
	y[0] = var0[1]
	z[0] = var0[2]		
	n = int(duree_t//dt)
	t = np.linspace(0,duree_t,n)	
	for i in range(1,n):		
		tspan =[ t[i-1] , t[i] ]		
		var = odeint(model,var0,t,args=(Q(Q_variable,t[i]),))		
		var0 = var[1]
		x[i] = var0[0]
		y[i] = var0[1]
		z[i] = var0[2]		
	plt.plot(t , z , "b" , linewidth=1, label = 'Concentration en polluant dans le lac, Q = Q(t)*' )
	plt.xlabel('Temps (s)',fontsize=14)
	plt.ylabel('Concentration',fontsize=14)	
	var0= [x0,y0,z0]
	q = Qconstant(var0,obj,eps)[0]
	var = odeint( model,var0,t,args=(q,) )
	Z = var[:,2]	
	plt.plot(t , Z  , "r" , linewidth=1, label = 'Concentration en polluant dans le lac, Q = Q constant*')
	plt.title('Concentration en polluant dans le lac en fonction du temps')
	plt.legend()
	plt.show()
	
dessin_Z_opti(.1,0.01)
