# Alignement de séquences

Le but de ce projet est de pouvoir rechercher efficacement des sous-séquences dans un génome
en suivant la stratégie de l’outil bowtie. Elle consiste à calculer un index de la séquence de référence
en utilisant la transformée de Burrow-Wheeler (BW) puis d’utiliser une procédure de recherche
efficace sur la séquence transformée.

## Exécution du code

Ce projet est fait en *Python*, pour lancer le code il suffit de le lancer de cette façon :

```sh
./script.py args
```
ou
```sh
python script.py args
```

Tous les scripts *Python* se trouve dans le dossier `src`.
Les bibliothèques utilisées font parties des bibliothèques standard de Python, il n'y a rien d'autres à installer mis
à part Python3 pour pouvoir utiliser ce projet.

## Utilisation des scripts

### bw-build

Ce script prend en entrée un fichier contenant une séquence sur laquelle on veut appliquer la transformée de
Burrow-Wheeler, un fichier de sortie, une fréquence de création d'index et des paramètres optionnels `--compress` pour
compresser le fichier de sortie et `--progressive k` pour faire une transformée progressive et ainsi pouvoir appliquer
cette transformée à de grandes séquences.

Il est possible d'avoir plus d'informations avec l'option `--help`.

### bw-search

Ce script prend en entrée un fichier contenant une transformée de Burrow-Wheeler (sortie du script précédent) et une
séquence à chercher dans la séquence initiale.
Les indexes où la séquence recherchée commence seront ensuite affichés sur la sortie standard. Il est également possible
d'ajouter l'option `--count-only` pour ne pas afficher toutes les séquences mais juste le nombre de séquences trouvées.

Il est possible d'avoir plus d'informations avec l'option `--help`.