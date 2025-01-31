import math, random
import numpy as np
import matplotlib.pyplot as plt

cylindres =  np.loadtxt("maps/donnees-map-1.txt")
fileOut = open("scripts/script-1.txt", "w")

# Variables globales
b = 3
b0 = 100
V0 = 1
a = 0.0698
etat = [0.0, 0.0, 90]
qteFuel = 10000.0

# Paramètres de cout 
a0 = 1 
a1 = 1

# Paramètres algo fourmie
gamma = 0.05
alpha = 1
beta = 2.0
rho = 0.5
nb_fourmis = 20
Q = 0.3

# Ecriture dans le fichier de sortie
def ecrireAvancer(distance):
    fileOut.write("GO %f\r" % (distance))
def ecrireTourner(angle):
    fileOut.write("TURN %f\r" % (angle))
def ecrireFinish():
    fileOut.write("FINISH\r\n")

# Calcul de la distance euclidienne
def distanceEuclidienne(p1, p2):
    return math.sqrt(math.pow((p1[0]-p2[0]), 2)+math.pow((p1[1]-p2[1]), 2))

# Permet de savoir si un point se situe dans un segment
def is_between(a,c,b):
    return round(distanceEuclidienne(a,c) + distanceEuclidienne(c,b), 3) == round(distanceEuclidienne(a,b), 3)

# Création du graphe des possibilités entre les cylindres
def graphePossibilites(cylindres, taille):
    res=np.zeros((taille,taille))
    for i in range(taille):
        for j in range(i, taille):
            if (i==j):
                res[i][j]=0
            else: 
                bool = True
                for k in range(taille):
                    if(k!=i and k!=j):
                        vectDirect = np.array([cylindres[j][0]-cylindres[i][0], cylindres[j][1]-cylindres[i][1]])
                        vectProj = np.array([cylindres[k][0]-cylindres[i][0], cylindres[k][1]-cylindres[i][1]])
                        proj = (np.dot(vectProj, vectDirect)/np.dot(vectDirect, vectDirect)) * vectDirect
                        point_proj = (cylindres[i][0], cylindres[i][1])+proj
                        if(is_between((cylindres[i][0], cylindres[i][1]), point_proj, (cylindres[j][0], cylindres[j][1]))):
                            if(distanceEuclidienne(point_proj, (cylindres[k][0], cylindres[k][1]))<0.5):
                                bool=False
                res[i][j]=bool 
    return res

pos0=[0,0]
orientation0 = 0
def CylindreToText(pos0,orientation0,cylindreOrdre):
    global cylindres
    pos=pos0
    orientation=orientation0
    for i in cylindreOrdre:
        Vecteur = [cylindres[i][0]-pos[0],cylindres[i][1]-pos[1]]
        Distance = np.sqrt(Vecteur[0]**2+Vecteur[1]**2)
        Vecteur = Vecteur / Distance
        a=np.arcsin(Vecteur[1])
        b=np.arccos(Vecteur[0])
        if Vecteur[1]<0:
            NouvelleOrientation = -b
        else:
            NouvelleOrientation = b
        NouvelleOrientation = NouvelleOrientation*180/np.pi
        OP = NouvelleOrientation-orientation
        if OP>180:
            OP=OP-360
        elif OP<-180:
            OP=OP+360
        ecrireTourner(OP)
        ecrireAvancer(Distance)
        orientation=NouvelleOrientation
        pos=[cylindres[i][0],cylindres[i][1]]
    ecrireFinish()
    return None

def coutCarb(mu, cylindres):
    global b, b0
    q = b*mu + b0
    res = np.zeros((taille, taille))
    for i in range(taille):
        for j in range(taille):
            res[i, j] = q*distanceEuclidienne(cylindres[i,0:2], cylindres[j,0:2])
    return res

def coutCarbPoint(pt1, pt2, mu, cylindres):
    global b, b0
    q = b*mu + b0
    return q*distanceEuclidienne(cylindres[pt1,0:2], cylindres[pt2,0:2])

def coutTemps(graphe, mu, num, cylindres):
    global a, V0
    Vmax = V0*math.exp(-a*mu)
    res = np.zeros((taille))
    for i in range(taille):
        if(graphe[i]==1):
            res[i] = distanceEuclidienne(cylindres[i,0:2], cylindres[num,0:2])/Vmax
    return res

def coutTempsPoint(pt1, pt2, mu, cylindres):
    global a, V0
    Vmax = V0*math.exp(-a*mu)
    return distanceEuclidienne(cylindres[pt1,0:2], cylindres[pt2,0:2])/Vmax

def coutTotPoint(pt1, pt2, mu, cylindres):
    return a0*coutTempsPoint(pt1, pt2, mu, cylindres) + a1* coutCarbPoint(pt1, pt2, mu, cylindres)

def tracer_chemin(coords, chemin): 
    plt.figure(figsize = (30,30))
    for i in range(len(chemin)-1): 
        p1=coords[chemin[i]]
        p2=coords[chemin[i+1]]
        plt.plot([p1[0], p2[0]], [p1[1], p2[1]], 'b-o')

    plt.scatter(coords[:, 0], coords[:, 1], c='red', label='Points')
    for i,(x,y) in enumerate(coords): 
        plt.text(x,y,f'{x}', fontsize = 12, ha='right')
    plt.title('Chemin opti')
    plt.xlabel('X')
    plt.ylabel('Y')
    plt.legend()
    plt.grid(True)
    plt.show()


def conversionTypeCylindre(typeC): 
    if(typeC==1.0):
        return [1, 1]
    if(typeC==2.0):
        return [2, 2]
    if(typeC==3.0):
        return [2, 3]

taille=len(cylindres)
cylindres = np.concatenate((cylindres, [[0.0, 0.0, 0.0]]), axis=0)
graphePossible = graphePossibilites(cylindres, taille+1)
graphePossibleInit = graphePossible.T + graphePossible

sum_value = 0
pb = 0
nb_cylindres = taille +1 

pheromones = np.ones((taille+1, taille+1), dtype=float)

def nextChoice(point, vuuuuu, pheromones, couts, mu):
    global alpha, beta
    nb_options = len(couts)
    proba = []

    for i in range(nb_options):
        if i in vuuuuu :
            proba.append(0)
        else:
            pheromone = math.pow(pheromones[point, i], alpha)
            desirabilite = math.pow(1/coutTotPoint(point, i, mu, cylindres), beta)
            proba.append(pheromone*desirabilite)
    proba = np.array(proba)
    proba /= np.sum(proba)
    return np.random.choice(range(nb_options), p=proba)

def fourmis(iterations):
    global pheromones, alpha, beta, rho, nb_fourmis, cylindres, taille
    meilleur_chemin = None
    meilleur_cout = float('inf')
    meilleur_score = None
    mu = 0 
    couts = coutCarb(mu, cylindres)
    nb_cylindres = len(couts)

    for _ in range(iterations):
        chemins = []
        coutCarb_totals = []
        coutTemps_totals = []
        scores = []

        for _ in range(nb_fourmis):
            mu = 0 
            value = 0
            act = taille 
            vuuus = [act]
            coutCarbTotal = 0 
            coutTempsTotal = 0 

            while(len(vuuus)<nb_cylindres+1 ): #and coutCarbTotal<qteFuel
                choix = nextChoice(act, vuuus, pheromones, couts, mu)
                coutCarbTotal += coutCarbPoint(act, choix, mu, cylindres)
                coutTempsTotal += coutTempsPoint(act, choix, mu, cylindres)
                vuuus.append(choix)
                act = choix
                mu += conversionTypeCylindre(cylindres[act, 2])[0]
                value += conversionTypeCylindre(cylindres[act, 2])[1]
            
            chemins.append(vuuus)
            scores.append(value)
            coutCarb_totals.append(coutCarbTotal)
            coutTemps_totals.append(coutTempsTotal)

        for i in range(nb_fourmis):
            if(coutCarb_totals[i]<meilleur_cout):
                meilleur_chemin = chemins[i]
                meilleur_score = scores[i]
                meilleur_cout = coutCarb_totals[i]
        
        pheromones *= (1 - rho)
        for i in range(nb_fourmis):
            chemin = chemins[i]
            cout = coutCarb_totals[i]
            for j in range(len(chemin) - 1):
                u, v = chemin[j], chemin[j+1]
                pheromones[u, v] += 1/cout
                pheromones[v, u] += 1/cout

    return meilleur_chemin, meilleur_cout, meilleur_score
    
coords = cylindres[:, :2]
chemin, cout, score = fourmis(200)

#tracer_chemin(coords, chemin)

CylindreToText(pos0, orientation0, np.delete(chemin, 0))

print(chemin)
print(cout)
print(score)
fileOut.close()