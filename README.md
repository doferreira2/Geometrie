# Projet CGAL : Segementation Mesh

## Compilation

Pour compiler le programme, il faut se rendre dans le répertoire `build` et exécuter la commande `cmake ..` pour préparer la compilation. Ensuite, il faut compiler le programme avec la commande `make`.

## Utilisation

Le programme implémente la segmentation d'un fichier .OFF. Pour tester le programme avec un fichier d'entrée `.off`, on peut utiliser la commande suivante à partir du répertoire `build` :

```
./color ../data/cube.off
```

On peut également ajouter deux paramètres optionnels :

Le premier paramètre est le type de mesure à utiliser pour la segmentation. On peut utiliser soit `P` pour la mesure de périmètre, soit `A` pour la mesure d'aire.
Le deuxième paramètre est le nombre de classes initiales, qui définit également le seuil pour la segmentation complexe.

Voici un exemple d'utilisation avec les deux paramètres optionnels :

```
./color ../data/pig.off P 15
```

Ce qui produit la sortie suivante:

```
Nb Class avec segmentation simple : 6
Le résultat a été exporté dans result.off !
Nb Class avec segmentation CC : 301
Le résultat a été exporté dans resultCC.off !
Nb Class avec segmentation complex : 3
Le résultat a été exporté dans resultcompl.off !
```

Cette commande affiche le nombre de classes obtenues avec la segmentation simple, la segmentation en composantes connexes (CC) et la segmentation complexe. Le résultat est également exporté dans les fichiers `result.off`, `resultCC.off` et `resultcompl.off`, qui sont situés dans le dossier `build`.

