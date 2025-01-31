# challenge_robotique

Création d'un algorithme d'optimisation.

## Contexte 

Dans le cadre de l'UV challenge à IMT Nord Europe, nous avons du créer un algorithme permettant d'optimiser le parcours d'un robot.
L'objectif du robot est de récupérer des cylindres de valeur sur une carte en minimisant le coût en carburant dans un temps imparti.
Ce parcours est simulé dans un environnement unity et nous avions le droit à python avec la librairie numpy pour créer l'algorithme.

## Algorithme

Nous nous sommes basés sur un algorithme de colonnies de fourmis avec le dépôt de phéromones pour les trajets les plus valorisés. 
Ceci nous permet d'obtenir des résultats très satisfaisants en une centaine d'itérations seulement. 

L'algorithme calcule à chaque itération et pour chaque foumis les probabilités qu'elle se déplace vers chaque nouveau cylindre selon le coût du trajet et la quantité de phéromones déposée sur ce trajet. 
Il sélectionne ensuite un des cylindres et recommence jusqu'à avoir tout visiter.
L'algorithme sélectionne la plus performante des fourmis et ajoute des phéromones en fonction du score fourni par celle-ci.