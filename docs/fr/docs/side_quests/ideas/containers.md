# Partie 1 : Approfondissement des Conteneurs

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Traduction assistée par IA - [en savoir plus et suggérer des améliorations](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

[TODO]

---

## 1. Comment trouver ou créer des images de conteneurs

Certain·es développeur·ses de logiciels fournissent des images de conteneurs pour leurs logiciels disponibles sur des registres de conteneurs comme Docker Hub, mais beaucoup ne le font pas.
Dans cette section optionnelle, nous vous montrerons deux façons d'obtenir une image de conteneur pour les outils que vous souhaitez utiliser dans vos pipelines Nextflow : en utilisant Seqera Containers et en construisant vous-même l'image de conteneur.

Vous allez obtenir/construire une image de conteneur pour le package pip `quote`, qui sera utilisé dans l'exercice à la fin de cette section.

### 1.1. Obtenir une image de conteneur depuis Seqera Containers

Seqera Containers est un service gratuit qui construit des images de conteneurs pour les outils installables via pip et conda (y compris bioconda).
Accédez à [Seqera Containers](https://www.seqera.io/containers/) et recherchez le package pip `quote`.

![Seqera Containers](img/seqera-containers-1.png)

Cliquez sur "+Add" puis "Get Container" pour demander une image de conteneur pour le package pip `quote`.

![Seqera Containers](img/seqera-containers-2.png)

Si c'est la première fois qu'un conteneur communautaire est construit pour cette version du package, cela peut prendre quelques minutes.
Cliquez pour copier l'URI (par exemple `community.wave.seqera.io/library/pip_quote:ae07804021465ee9`) de l'image de conteneur qui a été créée pour vous.

Vous pouvez maintenant utiliser l'image de conteneur pour exécuter la commande `quote` et obtenir une citation aléatoire de Grace Hopper.

```bash
docker run --rm community.wave.seqera.io/library/pip_quote:ae07804021465ee9 quote "Grace Hopper"
```

Sortie :

```console title="Sortie"
Humans are allergic to change. They love to say, 'We've always done it
this way.' I try to fight that. That's why I have a clock on my wall
that runs counter-clockwise.
```

### 1.2. Construire l'image de conteneur vous-même

Utilisons quelques détails de construction du site web Seqera Containers pour construire nous-mêmes l'image de conteneur pour le package pip `quote`.
Retournez sur le site web Seqera Containers et cliquez sur le bouton "Build Details".

Le premier élément que nous allons examiner est le `Dockerfile`, un type de fichier script qui contient toutes les commandes nécessaires pour construire l'image de conteneur.
Nous avons ajouté quelques commentaires explicatifs au Dockerfile ci-dessous pour vous aider à comprendre ce que fait chaque partie.

```Dockerfile title="Dockerfile"
# Partir de l'image docker de base micromamba
FROM mambaorg/micromamba:1.5.10-noble
# Copier le fichier conda.yml dans le conteneur
COPY --chown=$MAMBA_USER:$MAMBA_USER conda.yml /tmp/conda.yml
# Installer divers utilitaires pour que Nextflow les utilise et les packages dans le fichier conda.yml
RUN micromamba install -y -n base -f /tmp/conda.yml \
    && micromamba install -y -n base conda-forge::procps-ng \
    && micromamba env export --name base --explicit > environment.lock \
    && echo ">> CONDA_LOCK_START" \
    && cat environment.lock \
    && echo "<< CONDA_LOCK_END" \
    && micromamba clean -a -y
# Exécuter le conteneur en tant qu'utilisateur root
USER root
# Définir la variable d'environnement PATH pour inclure le répertoire d'installation de micromamba
ENV PATH="$MAMBA_ROOT_PREFIX/bin:$PATH"
```

Le deuxième élément que nous allons examiner est le fichier `conda.yml`, qui contient la liste des packages qui doivent être installés dans l'image de conteneur.

```conda.yml title="conda.yml"
channels:
- conda-forge
- bioconda
dependencies:
- pip
- pip:
  - quote==3.0.0 #
```

Copiez le contenu de ces fichiers dans les emplacements situés dans le répertoire `containers/build`, puis exécutez la commande suivante pour construire l'image de conteneur vous-même.

!!! Note

    Nous utilisons l'option `-t quote:latest` pour étiqueter l'image de conteneur avec le nom `quote` et l'étiquette `latest`.
    Nous pourrons utiliser cette étiquette pour faire référence à l'image de conteneur lors de son exécution sur ce système.

```bash
docker build -t quote:latest containers/build
```

Une fois la construction terminée, vous pouvez exécuter l'image de conteneur que vous venez de construire.

```bash
docker run --rm quote:latest quote "Margaret Oakley Dayhoff"
```

### À retenir

Vous avez appris deux façons différentes d'obtenir une image de conteneur pour un outil que vous souhaitez utiliser dans vos pipelines Nextflow : en utilisant Seqera Containers et en construisant l'image de conteneur vous-même.

### Et ensuite ?

Vous avez tout ce dont vous avez besoin pour continuer au [chapitre suivant](./04_hello_genomics.md) de cette série de formation.
Vous pouvez également poursuivre avec un exercice optionnel pour récupérer des citations sur les pionniers de l'informatique/biologie en utilisant le conteneur `quote` et les afficher en utilisant le conteneur `cowsay`.

---

## 2. Faire dire des citations de scientifiques célèbres à la vache

Cette section contient quelques exercices avancés, pour pratiquer ce que vous avez appris jusqu'à présent.
Faire ces exercices n'est _pas requis_ pour comprendre les parties suivantes de la formation, mais ils offrent une façon amusante de renforcer vos acquis en déterminant comment faire dire des citations de scientifiques célèbres à la vache.

```console title="cowsay-output-Grace-Hopper.txt"
  _________________________________________________
 /                                                 \
| Humans are allergic to change. They love to       |
| say, 'We've always done it this way.' I try to fi |
| ght that. That's why I have a clock on my wall th |
| at runs counter-clockwise.                        |
| -Grace Hopper                                     |
 \                                                 /
  =================================================
                                                 \
                                                  \
                                                    ^__^
                                                    (oo)\_______
                                                    (__)\       )\/\
                                                        ||----w |
                                                        ||     ||
```

### 2.1. Modifier le script `hello-containers.nf` pour utiliser un processus getQuote

Nous avons une liste de pionniers de l'informatique et de la biologie dans le fichier `containers/data/pioneers.csv`.
À un niveau élevé, pour compléter cet exercice vous devrez :

- Modifier le `params.input_file` par défaut pour pointer vers le fichier `pioneers.csv`.
- Créer un processus `getQuote` qui utilise le conteneur `quote` pour récupérer une citation pour chaque entrée.
- Connecter la sortie du processus `getQuote` au processus `cowsay` pour afficher la citation.

Pour l'image de conteneur `quote`, vous pouvez soit utiliser celle que vous avez construite vous-même dans l'exercice avancé précédent, soit utiliser celle que vous avez obtenue depuis Seqera Containers.

!!! Hint "Indice"

    Un bon choix pour le bloc `script` de votre processus getQuote pourrait être :
        ```groovy
        script:
            def safe_author = author.tokenize(' ').join('-')
            """
            quote "$author" > quote-${safe_author}.txt
            echo "-${author}" >> quote-${safe_author}.txt
            """
        ```

Vous pouvez trouver une solution à cet exercice dans `containers/solutions/hello-containers-4.1.nf`.

### 2.2. Modifier votre pipeline Nextflow pour lui permettre de s'exécuter en modes `quote` et `sayHello`.

Ajoutez une logique de branchement à votre pipeline pour lui permettre d'accepter des entrées destinées à la fois à `quote` et `sayHello`.
Voici un exemple de la façon d'utiliser une instruction `if` dans un workflow Nextflow :

```groovy title="hello-containers.nf"
workflow {
    if (params.quote) {
        ...
    }
    else {
        ...
    }
    cowSay(text_ch)
}
```

!!! Hint "Indice"

    Vous pouvez utiliser `new_ch = processName.out` pour assigner un nom au canal de sortie d'un processus.

Vous pouvez trouver une solution à cet exercice dans `containers/solutions/hello-containers-4.2.nf`.

### À retenir

Vous savez comment utiliser les conteneurs dans Nextflow pour exécuter des processus, et comment construire une logique de branchement dans vos pipelines !

### Et ensuite ?

Félicitations, prenez une pause et buvez de l'eau !

Lorsque vous êtes prêt·e, passez à la Partie 3 de cette série de formation pour apprendre comment appliquer ce que vous avez appris jusqu'à présent à un cas d'utilisation d'analyse de données plus réaliste.
