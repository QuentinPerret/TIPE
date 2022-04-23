import matplotlib.pyplot as plt
import  numpy as np
from scipy.integrate import odeint

#Initialisation des constantes du problèmes
duree_s = 20. #durée de l'expérience (calcul de x et y)
duree_t = 5000. #durée de l'expérience (calcul de z)
dt = 0.05  #durée d'un interval de temps entre 2 mesures
m = 1. #constante de la loi de création de la biomasse
Vr = 1. #volume du bioréacteur
Vl = 100. # volume del'étendu d'eau 
z_cst = 1
z0 = 1

#Définition des foonction x' , y' et z'
def model(var,t,Q):	
	x = var[0]
	y = var[1]	
	dxdt =  m*x*y - Q/Vr*x
	dydt = Q/Vr*(z_cst - y) - m*x*y	
	return [dxdt , dydt]
	
#On doit définir deux model puisque les équations de x' et y' 
#ne sont pas couplées avec celle de z' dans l'approximation lent-rapide

def model_z(var,t,Q,y_cst):	
	z = var[0]	
	dzdt = Q / Vl  * (y_cst - z)	
	return [dzdt]

#Conditions Initiales
x0 = 0.1 #valeur initiale de x
y0 = 1. #valeur initiale de y
var0= [x0,y0]

def dessin_X_Y(Q):	
	'''Fonction dessinant les courbes x et y en fonction de t pour Q donné'''	
	#Création du tableau de temps 
	n = int(duree_s//dt)
	t = np.linspace(0,duree_s,n)	
	#Création des tableaux contenant x et y
	x = np.zeros(n)
	y = np.zeros(n)	
	#Calcul de toutes les valeurs de x et y pour chaque élément de t 
	var = odeint(model,var0,t,args= (Q,))	
	X,Y  = var[:,0] , var[:,1]	
	#Affichage des courbes de x et y en fonction du temps 		
	plt.plot(t , X ,'r', linewidth=1, label = 'Concentration en biomasse, Q = ' + str(Q) )
	plt.plot(t , Y , 'b', linewidth=1, label = 'Concentration en polluant, Q = ' + str(Q) )
	plt.xlabel('Temps (s)',fontsize=14)
	plt.ylabel('Concentration',fontsize=14)
	plt.title('Modèle approximation lent_rapide')
	plt.legend()
	plt.show()

def dessin_Z(Q):	
	'''Fonction dessinant z en fonction de t pour Q donné'''
	#Création du tableau de temps 
	n = int(duree_t//dt)
	t = np.linspace(0,duree_t,n)
	#création des tableaux contenant x, y et z 
	x = np.zeros(n)
	y = np.zeros(n)
	z = np.zeros(n)
	#Calcul de x et y pour chaque élément de t 
	xy = odeint(model,var0,t, args = (Q,))
	#Calcul de z pour chaque élément de t 
	y_cst = xy[int(duree_s/dt),1]	#on définit y_inf en prenant la dernière valeur de y pour t = duree_s
	var = odeint(model_z,z0,t,args=(Q,y_cst))	
	#Affichage de courbe de z en fonction du temps  
	plt.plot(t , var  , linewidth=1, label = 'Concentration en polluant dans le lac, Q = ' + str(Q))
	plt.xlabel('Temps (s)',fontsize=14)
	plt.ylabel('Concentration',fontsize=14)
	plt.title("Modèle approximation lent_rapide")
	plt.legend()
	plt.show()

dessin_X_Y(0.05)
dessin_Z(0.05)
