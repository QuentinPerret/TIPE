import matplotlib.pyplot as plt
import  numpy as np
from scipy.integrate import odeint

#Initialisation des constantes du problèmes
duree_s = 2000. #durée de l'expérience (calcul de x et y)
duree_t = 5000. #durée de l'expérience (calcul de z)
dt = 1.  #durée d'un interval de temps entre 2 mesures
m = 1. #constante de la loi de création de la biomasse
Vr = 1. #volume du bioréacteur
Vl = 100. # volume del'étendu d'eau 

#Définition des foonction x' , y' et z'
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

#Définition des tableaux des variables
x = np.zeros(n)
y = np.zeros(n)
z = np.zeros(n)

def Qconstant(var0,obj,esp):	
	"""	Entrée : var0 constantes initiales du problème , obj objectif de dépollution , eps précision de la valeur du débit 
		Sortie : débit contant optimal pour atteindre obj le plus rapidement possible Q est précis à eps près et l'indice dans le tableau t du temps mis pour atteindre l'objectif obj """		
	q = np.linspace(0,1,int(1/esp)) #on définit toutes les valeurs que peut prendre le débit q
	Q_min , m = 0 , n #on crée les variables pour contenir les résultats à renvoyer	
	for Q in q: #on parcours tous les débits possibles		
		Z = odeint(model,var0,t,args=(Q,))[:,2] #pour chaque débit on calcul les valeurs de z à chaque instant 
		j = 0 
		while j < n and Z[j] > obj: #tant que la valeur de Z est plus grande que l'objectif de dépollution j += 1 	
			j += 1 		
		#Ici on a j l'indice du temps dans t du temps nécéssaire pour atteindre l'objectif obj		
		if j < m: #si j est plus petit que mù soit le temps mis par Q inférieur à celui par Q_min pour atteindre le plus rapidement l'obkectif de dépollution 			
			Q_min = Q
			m = j #alors on change les valeurs de Q_min et m
	return  Q_min , m #on renvoie le débit contant optimal pour atteindre obj le plus rapidement possible Q est précis à eps près et l'indice dans le tableau t du temps mis pour atteindre l'objectif obj

def dessin_Z_opti(obj,eps):
	'''Fonction dessinant z en fonction de t pour Q contant optimal'''
	Q = Qconstant(var0,obj,eps)[0]
	var = odeint( model,var0,t,args=(Q,) )
	Z = var[:,2]
	label3 = 'Concentration en polluant dans le lac, Q = ' + str(Q)[:len(str(eps))]
	plt.plot(t , Z  , linewidth=1, label = label3)
	plt.xlabel('Temps (s)',fontsize=14)
	plt.ylabel('Concentration',fontsize=14)
	plt.title('Concentration en polluant dans le lac en fonction du temps \navec un débit constant,\n pour un objectif de dépollution à '+ str(obj))	
	plt.legend()
	plt.show()

dessin_Z_opti(0.1,0.01)
