# Exemple de projet CGAL

## Compilation

* se rendre dans le répertoire ```src```, et lancer la commande ```cgal_create_CMakeLists``` pour générer le CMakeLists.txt correspondant aux sources ```.cpp``` disponibles dans ce dossier
* se rendre dans le répertoire ```build```, puis préparer la compilation avec ```cmake ..```
* compiler avec la commande ```make```

## Utilisation

Le seul programme pour l'instant implémenté est un programme qui calcule le genre d'une surface sans bord à une seul composante connexe. On peut le tester avec la commande suivante, depuis le répertoire ```build```:

```
src/genre ../data/cube.off
```

Ce qui produit la sortie suivante:

```
Nombre de sommets: 8
Nombre d'arêtes: 12
Nombre de faces: 6
En supposant que le maillage contienne une unique surface sans bord, alors son genre est de 0
````

## Ajout d'un nouveau programme

Tous les fichiers au format c++ présents dans le répertoire ```src``` seront considérés comme des nouveaux programmes CGAL par la commande ```cgal_create_CMakeLists``` lancée depuis ```src```. Pensez à utiliser un sous-répertoire si vous voulez faire des ```include```.
