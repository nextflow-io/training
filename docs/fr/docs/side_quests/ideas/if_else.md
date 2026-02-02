# Partie 2 : If - Else

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Traduction assistée par IA - [en savoir plus et suggérer des améliorations](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

[TODO]

---

## 1. Faire citer des scientifiques célèbres à la vache

Cette section contient quelques exercices supplémentaires pour pratiquer ce que vous avez appris jusqu'à présent.
Réaliser ces exercices n'est _pas obligatoire_ pour comprendre les parties suivantes de la formation, mais ils constituent une façon amusante de renforcer vos apprentissages en découvrant comment faire citer des scientifiques célèbres à la vache.

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

### 1.1. Modifier le script `hello-containers.nf` pour utiliser un processus getQuote

Nous avons une liste de pionniers de l'informatique et de la biologie dans le fichier `containers/data/pioneers.csv`.
De manière générale, pour compléter cet exercice, vous devrez :

- Modifier le paramètre par défaut `params.input_file` pour qu'il pointe vers le fichier `pioneers.csv`.
- Créer un processus `getQuote` qui utilise le conteneur `quote` pour récupérer une citation pour chaque entrée.
- Connecter la sortie du processus `getQuote` au processus `cowsay` pour afficher la citation.

Pour l'image de conteneur `quote`, vous pouvez soit utiliser celle que vous avez construite vous-même dans l'exercice supplémentaire précédent, soit utiliser celle que vous avez obtenue depuis Seqera Containers.

!!! Hint

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

### 1.2. Modifier votre pipeline Nextflow pour lui permettre de s'exécuter en modes `quote` et `sayHello`.

Ajoutez une logique de branchement à votre pipeline pour lui permettre d'accepter des entrées destinées à la fois à `quote` et `sayHello`.
Voici un exemple d'utilisation d'une instruction `if` dans un workflow Nextflow :

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

!!! Hint

    Vous pouvez utiliser `new_ch = processName.out` pour attribuer un nom au canal de sortie d'un processus.

Vous pouvez trouver une solution à cet exercice dans `containers/solutions/hello-containers-4.2.nf`.

### À retenir

Vous savez comment utiliser des conteneurs dans Nextflow pour exécuter des processus, et comment construire une logique de branchement dans vos pipelines !

### Et ensuite ?

Félicitations, prenez une pause pour vous étirer et boire de l'eau !

Lorsque vous êtes prêt·e, passez à la Partie 3 de cette série de formation pour apprendre à appliquer ce que vous avez appris jusqu'à présent à un cas d'usage d'analyse de données plus réaliste.
