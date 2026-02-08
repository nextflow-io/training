# Partie 2 : Hello Channels

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Traduction assistée par IA - [en savoir plus et suggérer des améliorations](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

<div class="video-wrapper">
  <iframe width="560" height="315" src="https://www.youtube.com/embed/yDR66fzAMOg?si=xCItHLiOQWqoqBB9&amp;list=PLPZ8WHdZGxmWKozQuzr27jyMGqp9kElVK&amp;cc_load_policy=1&amp;cc_lang_pref=fr" title="YouTube video player" frameborder="0" allow="accelerometer; autoplay; clipboard-write; encrypted-media; gyroscope; picture-in-picture; web-share" referrerpolicy="strict-origin-when-cross-origin" allowfullscreen></iframe>
</div>

/// caption
:fontawesome-brands-youtube:{ .youtube } Voir [la playlist complète](https://youtube.com/playlist?list=PLPZ8WHdZGxmWKozQuzr27jyMGqp9kElVK&si=eF7cLR62goy-lc6n) sur la chaîne YouTube Nextflow.

:green_book: La transcription de la vidéo est disponible [ici](./transcripts/02_hello_channels.md).
///

Dans la Partie 1 de ce cours (Hello World), nous vous avons montré comment fournir une entrée variable à un processus en fournissant l'entrée directement dans l'appel du processus : `sayHello(params.input)`.
C'était une approche délibérément simplifiée.
En pratique, cette approche a des limitations majeures ; à savoir qu'elle ne fonctionne que pour des cas très simples où nous voulons seulement exécuter le processus une fois, sur une seule valeur.
Dans la plupart des cas d'utilisation de workflow réalistes, nous voulons traiter plusieurs valeurs (données expérimentales pour plusieurs échantillons, par exemple), nous avons donc besoin d'une façon plus sophistiquée de gérer les entrées.

C'est à cela que servent les [**canaux**](https://nextflow.io/docs/latest/channel.html) Nextflow.
Les canaux sont des files d'attente conçues pour gérer les entrées efficacement et les faire passer d'une étape à l'autre dans des workflows à plusieurs étapes, tout en fournissant un parallélisme intégré et de nombreux autres avantages.

Dans cette partie du cours, vous apprendrez à utiliser un canal pour gérer plusieurs entrées provenant de diverses sources différentes.
Vous apprendrez également à utiliser des [**opérateurs**](https://nextflow.io/docs/latest/reference/operator.html) pour transformer le contenu des canaux selon les besoins.

??? info "Comment commencer à partir de cette section"

    Cette section du cours suppose que vous avez terminé la Partie 1 du cours [Hello Nextflow](./index.md), mais si vous êtes à l'aise avec les bases couvertes dans cette section, vous pouvez commencer ici sans rien faire de spécial.

---

## 0. Échauffement : Exécuter `hello-channels.nf`

Nous allons utiliser le script de workflow `hello-channels.nf` comme point de départ.
Il est équivalent au script produit en suivant la Partie 1 de ce cours de formation, sauf que nous avons changé la destination de sortie :

```groovy title="hello-channels.nf" linenums="37" hl_lines="3"
output {
    first_output {
        path 'hello_channels'
        mode 'copy'
    }
}
```

Juste pour s'assurer que tout fonctionne, exécutez le script une fois avant d'apporter des modifications :

```bash
nextflow run hello-channels.nf --input 'Hello Channels!'
```

??? success "Sortie de la commande"

    ```console
     N E X T F L O W   ~  version 25.10.2

    Launching `hello-channels.nf` [wise_jennings] DSL2 - revision: b24f4902d6

    executor >  local (1)
    [6f/824bc1] process > sayHello [100%] 1 of 1 ✔
    ```

Comme précédemment, vous trouverez le fichier de sortie nommé `output.txt` dans le répertoire `results/hello_channels` (comme spécifié dans le bloc `output` du script de workflow, montré ci-dessus).

??? abstract "Contenu du répertoire"

    ```console title="results/hello_channels" hl_lines="2-3"
    results
    ├── hello_channels
    │   └── output.txt
    ├── hello_world
    │   └── output.txt
    └── output.txt -> /workspaces/training/hello-nextflow/work/8c/79499c11beea6e9d43605141f2817f/output.txt
    ```

??? abstract "Contenu du fichier"

    ```console title="results/hello_channels/output.txt"
    Hello Channels!
    ```

Si cela a fonctionné pour vous, vous êtes prêt·e à apprendre les canaux.

---

## 1. Fournir des entrées variables via un canal explicitement

Nous allons créer un **canal** pour passer l'entrée variable au processus `sayHello()` au lieu de compter sur la gestion implicite, qui a certaines limitations.

### 1.1. Créer un canal d'entrée

Il existe une variété de [**fabriques de canaux**](https://nextflow.io/docs/latest/reference/channel.html) que nous pouvons utiliser pour configurer un canal.
Pour garder les choses simples pour l'instant, nous allons utiliser la fabrique de canaux la plus basique, appelée [`channel.of`](https://nextflow.io/docs/latest/reference/channel.html#of), qui créera un canal contenant une seule valeur.
Fonctionnellement, ce sera similaire à la configuration que nous avions avant, mais au lieu de laisser Nextflow créer un canal implicitement, nous le faisons explicitement maintenant.

Voici la ligne de code que nous allons utiliser :

```console title="Syntaxe"
greeting_ch = channel.of('Hello Channels!')
```

Cela crée un canal appelé `greeting_ch` en utilisant la fabrique de canaux `channel.of()`, qui configure un simple canal de file d'attente, et charge la chaîne `'Hello Channels!'` à utiliser comme valeur de salutation.

<figure class="excalidraw">
--8<-- "docs/en/docs/hello_nextflow/img/hello-pipeline-channel.svg"
</figure>

!!! note "Note"

    Nous revenons temporairement aux chaînes codées en dur au lieu d'utiliser un paramètre CLI pour des raisons de lisibilité. Nous reviendrons à l'utilisation des paramètres CLI une fois que nous aurons couvert ce qui se passe au niveau du canal.

Dans le bloc workflow, ajoutez le code de la fabrique de canaux :

=== "Après"

    ```groovy title="hello-channels.nf" linenums="27" hl_lines="4 5"
    workflow {

        main:
        // créer un canal pour les entrées
        greeting_ch = channel.of('Hello Channels!')
        // émettre une salutation
        sayHello(params.input)

        publish:
        first_output = sayHello.out
    }
    ```

=== "Avant"

    ```groovy title="hello-channels.nf" linenums="27"
    workflow {

        main:
        // émettre une salutation
        sayHello(params.input)

        publish:
        first_output = sayHello.out
    }
    ```

Ce n'est pas encore fonctionnel puisque nous n'avons pas encore changé l'entrée vers l'appel du processus.

### 1.2. Ajouter le canal comme entrée à l'appel du processus

Maintenant nous devons effectivement brancher notre canal nouvellement créé dans l'appel du processus `sayHello()`, remplaçant le paramètre CLI que nous fournissions directement avant.

Dans le bloc workflow, faites la modification de code suivante :

=== "Après"

    ```groovy title="hello-channels.nf" linenums="27" hl_lines="7"
    workflow {

        main:
        // créer un canal pour les entrées
        greeting_ch = channel.of('Hello Channels!')
        // émettre une salutation
        sayHello(greeting_ch)

        publish:
        first_output = sayHello.out
    }
    ```

=== "Avant"

    ```groovy title="hello-channels.nf" linenums="27" hl_lines="7"
    workflow {

        main:
        // créer un canal pour les entrées
        greeting_ch = channel.of('Hello Channels!')
        // émettre une salutation
        sayHello(params.input)

        publish:
        first_output = sayHello.out
    }
    ```

Cela indique à Nextflow d'exécuter le processus `sayHello` sur le contenu du canal `greeting_ch`.

Maintenant notre workflow est correctement fonctionnel ; c'est l'équivalent explicite d'écrire `sayHello('Hello Channels!')`.

### 1.3. Exécuter le workflow

Exécutons-le !

```bash
nextflow run hello-channels.nf
```

??? success "Sortie de la commande"

    ```console
     N E X T F L O W   ~  version 25.10.2

    Launching `hello-channels.nf` [fabulous_crick] DSL2 - revision: 23e20f76e8

    executor >  local (1)
    [c0/4f1872] process > sayHello (1) [100%] 1 of 1 ✔
    ```

Si vous avez fait les deux modifications correctement, vous devriez obtenir une exécution réussie.
Vous pouvez vérifier le répertoire des résultats pour vous assurer que le résultat est toujours le même qu'avant.

??? abstract "Contenu du fichier"

    ```console title="results/hello_channels/output.txt"
    Hello Channels!
    ```

Nous avons donc augmenté la flexibilité de notre workflow tout en obtenant le même résultat final.
Cela peut sembler que nous écrivons plus de code sans bénéfice tangible, mais la valeur deviendra claire dès que nous commencerons à gérer plus d'entrées.

En avant-première, regardons une chose de plus avant de passer à la suite : un petit avantage pratique de l'utilisation d'un canal explicite pour gérer l'entrée de données.

### 1.4. Utiliser `view()` pour inspecter le contenu du canal

Les canaux Nextflow sont construits d'une manière qui nous permet d'opérer sur leur contenu en utilisant des opérateurs, que nous couvrirons en détail plus tard dans ce chapitre.

Pour l'instant, nous allons juste vous montrer comment utiliser un opérateur super simple appelé [`view()`](https://www.nextflow.io/docs/latest/reference/operator.html#view) pour inspecter le contenu d'un canal.
Vous pouvez considérer `view()` comme un outil de débogage, comme une instruction `print()` en Python, ou son équivalent dans d'autres langages.

Ajoutez cette petite ligne au bloc workflow :

=== "Après"

    ```groovy title="hello-channels.nf" linenums="27" hl_lines="7"
    workflow {

        main:
        // créer un canal pour les entrées
        greeting_ch = channel.of('Hello Channels!')
                             .view()
        // émettre une salutation
        sayHello(greeting_ch)

        publish:
        first_output = sayHello.out
    }
    ```

=== "Avant"

    ```groovy title="hello-channels.nf" linenums="27"
    workflow {

        main:
        // créer un canal pour les entrées
        greeting_ch = channel.of('Hello Channels!')
        // émettre une salutation
        sayHello(greeting_ch)

        publish:
        first_output = sayHello.out
    }
    ```

La quantité exacte d'espaces n'a pas d'importance tant que c'est un multiple de 4 ; nous visons juste à aligner le début de l'instruction `.view()` avec la partie `.of()` de la construction du canal.

Maintenant exécutez à nouveau le workflow :

```bash
nextflow run hello-channels.nf
```

??? success "Sortie de la commande"

    ```console hl_lines="7"
     N E X T F L O W   ~  version 25.10.2

    Launching `hello-channels.nf` [scruffy_shaw] DSL2 - revision: 2ede41e14a

    executor >  local (1)
    [ef/f7e40a] sayHello (1) [100%] 1 of 1 ✔
    Hello Channels!
    ```

Comme vous pouvez le voir, cela affiche le contenu du canal dans la console.
Ici, nous n'avons qu'un élément, mais quand nous commencerons à charger plusieurs valeurs dans le canal dans la section suivante, vous verrez que cela est configuré pour afficher un élément par ligne.

### À retenir

Vous savez comment utiliser une fabrique de canaux de base pour fournir une entrée à un processus.

### Et ensuite ?

Apprendre à utiliser les canaux pour faire itérer le workflow sur plusieurs valeurs d'entrée.

---

## 2. Modifier le workflow pour s'exécuter sur plusieurs valeurs d'entrée

Les workflows s'exécutent généralement sur des lots d'entrées destinés à être traités en masse, nous voulons donc améliorer le workflow pour accepter plusieurs valeurs d'entrée.

### 2.1. Charger plusieurs salutations dans le canal d'entrée

De façon pratique, la fabrique de canaux `channel.of()` que nous utilisons est tout à fait capable d'accepter plus d'une valeur, donc nous n'avons pas besoin de la modifier du tout.
Nous pouvons simplement charger plusieurs valeurs dans le canal.

Mettons `'Hello'`, `'Bonjour'` et `'Holà'`.

#### 2.1.1. Ajouter plus de salutations

Avant le bloc workflow, faites la modification de code suivante :

=== "Après"

    ```groovy title="hello-channels.nf" linenums="30" hl_lines="2"
    // créer un canal pour les entrées
    greeting_ch = channel.of('Hello','Bonjour','Holà')
                         .view()
    ```

=== "Avant"

    ```groovy title="hello-channels.nf" linenums="30" hl_lines="2"
    // créer un canal pour les entrées
    greeting_ch = channel.of('Hello Channels')
                         .view()
    ```

La documentation nous dit que cela devrait fonctionner. Cela peut-il vraiment être aussi simple ?

#### 2.1.2. Exécuter la commande et regarder la sortie du log

Essayons.

```bash
nextflow run hello-channels.nf
```

??? success "Sortie de la commande"

    ```console hl_lines="6"
     N E X T F L O W   ~  version 25.10.2

    Launching `hello-channels.nf` [amazing_crick] DSL2 - revision: 59a9a5888a

    executor >  local (3)
    [f4/c9962c] process > sayHello (1) [100%] 3 of 3 ✔
    Hello
    Bonjour
    Holà
    ```

Cela semble avoir bien fonctionné.
Le moniteur d'exécution montre que `3 of 3` appels ont été effectués pour le processus `sayHello`, et nous voyons les trois salutations énumérées par l'instruction `view()`, une par ligne comme promis.

Cependant, il n'y a toujours qu'une seule sortie dans le répertoire des résultats :

??? abstract "Contenu du répertoire"

    ```console title="results/hello_channels" hl_lines="3"
    results
    ├── hello_channels
    │   └── output.txt
    ├── hello_world
    │   └── output.txt
    └── output.txt -> /workspaces/training/hello-nextflow/work/8c/79499c11beea6e9d43605141f2817f/output.txt
    ```

??? abstract "Contenu du fichier"

    ```console title="results/hello_channels/output.txt"
    Holà
    ```

Vous devriez voir l'une des trois salutations dedans, bien que celle que vous avez obtenue pourrait être différente de ce qui est montré ici.
Pouvez-vous deviner pourquoi ?

<figure class="excalidraw">
--8<-- "docs/en/docs/hello_nextflow/img/hello-pipeline-channel-multi.svg"
</figure>

_Dans le diagramme, le canal est représenté en vert, et l'ordre des éléments est représenté comme des billes dans un tuyau : le premier chargé est à droite, puis le deuxième au milieu, puis le troisième est à gauche._

En regardant le moniteur d'exécution, il nous a donné seulement un chemin de sous-répertoire (`f4/c9962c`).
Jetons un coup d'œil dedans.

??? abstract "Contenu du répertoire"

    ```console hl_lines="9"
    work/f4/c9962ce91ef87480babcb86b2b9042/
    ├── .command.begin
    ├── .command.err
    ├── .command.log
    ├── .command.out
    ├── .command.run
    ├── .command.sh
    ├── .exitcode
    └── output.txt
    ```

??? abstract "Contenu du fichier"

    ```console title="work/f4/c9962ce91ef87480babcb86b2b9042/output.txt"
    Hello
    ```

Ce n'est même pas la même salutation que nous avons obtenue dans le répertoire des résultats ! Que se passe-t-il ?

À ce stade, nous devons vous dire que par défaut, le système de journalisation ANSI écrit les logs de plusieurs appels au même processus sur la même ligne.
Donc le statut des trois appels au processus sayHello() atterrissent au même endroit.

Heureusement, nous pouvons désactiver ce comportement pour voir la liste complète des appels de processus.

#### 2.1.3. Exécuter à nouveau la commande avec l'option `-ansi-log false`

Pour développer la journalisation afin d'afficher une ligne par appel de processus, ajoutez `-ansi-log false` à la commande.

```bash
nextflow run hello-channels.nf -ansi-log false
```

??? success "Sortie de la commande"

    ```console
     N E X T F L O W  ~  version 25.10.2
    Launching `hello-channels.nf` [desperate_monod] DSL2 - revision: 59a9a5888a
    Hello
    Bonjour
    Holà
    [23/871c7e] Submitted process > sayHello (2)
    [7f/21e2c2] Submitted process > sayHello (1)
    [f4/ea10a6] Submitted process > sayHello (3)
    ```

Cette fois nous voyons les trois exécutions de processus et leurs sous-répertoires de travail associés listés dans la sortie.

C'est beaucoup mieux, au moins pour un workflow simple.
Pour un workflow complexe, ou un grand nombre d'entrées, avoir la liste complète affichée dans le terminal deviendrait un peu accablant.
C'est pourquoi `-ansi-log false` n'est pas le comportement par défaut.

!!! tip "Astuce"

    La façon dont le statut est rapporté est un peu différente entre les deux modes de journalisation.
    Dans le mode condensé, Nextflow rapporte si les appels se sont terminés avec succès ou non.
    Dans ce mode étendu, il rapporte seulement qu'ils ont été soumis.

De toute façon, maintenant que nous avons les sous-répertoires de chaque appel de processus, nous pouvons chercher leurs logs et sorties.

??? abstract "Contenu du répertoire"

    ```console
    work/23/871c7ec3642a898ecd5e6090d21300/
    ├── .command.begin
    ├── .command.err
    ├── .command.log
    ├── .command.out
    ├── .command.run
    ├── .command.sh
    ├── .exitcode
    └── output.txt
    ```

    ```console
    work/7f/21e2c2f3cc8833ef3858b236e5575c/
    ├── .command.begin
    ├── .command.err
    ├── .command.log
    ├── .command.out
    ├── .command.run
    ├── .command.sh
    ├── .exitcode
    └── output.txt
    ```

    ```console
    work/f4/ea10a680d5687596d3eaa3fcf69272/
    ├── .command.begin
    ├── .command.err
    ├── .command.log
    ├── .command.out
    ├── .command.run
    ├── .command.sh
    ├── .exitcode
    └── output.txt
    ```

??? abstract "Contenu des fichiers"

    ```txt title="work/23/871c7ec3642a898ecd5e6090d21300/output.txt"
    Bonjour
    ```

    ```txt title="work/7f/21e2c2f3cc8833ef3858b236e5575c/output.txt"
    Hello
    ```

    ```txt title="work/f4/ea10a680d5687596d3eaa3fcf69272/output.txt"
    Holà
    ```

Cela montre que les trois processus se sont exécutés avec succès (youpi).

Cela dit, nous avons toujours le problème qu'il n'y a qu'un seul fichier de sortie dans le répertoire des résultats.

Vous vous souvenez peut-être que nous avions codé en dur le nom du fichier de sortie pour le processus `sayHello`, donc les trois appels ont produit un fichier appelé `output.txt`.

<figure class="excalidraw">
--8<-- "docs/en/docs/hello_nextflow/img/hello-channels-task-dirs.svg"
</figure>

Tant que les fichiers de sortie restent dans les sous-répertoires de travail, isolés des autres processus, c'est correct.
Mais quand ils sont publiés dans le même répertoire de résultats, celui qui a été copié là en premier est écrasé par le suivant, et ainsi de suite.

### 2.2. S'assurer que les noms de fichiers de sortie seront uniques

Nous pouvons continuer à publier toutes les sorties dans le même répertoire de résultats, mais nous devons nous assurer qu'elles auront des noms uniques.
Spécifiquement, nous devons modifier le premier processus pour générer un nom de fichier dynamiquement afin que les noms de fichiers finaux soient uniques.

Alors comment rendre les noms de fichiers uniques ?
Une façon courante de le faire est d'utiliser une pièce unique de métadonnées des entrées (reçues du canal d'entrée) comme partie du nom du fichier de sortie.
Ici, pour plus de commodité, nous utiliserons simplement la salutation elle-même puisque c'est juste une courte chaîne, et la préfixerons au nom de fichier de sortie de base.

<figure class="excalidraw">
--8<-- "docs/en/docs/hello_nextflow/img/hello-pipeline-channel-multi-unique.svg"
</figure>

#### 2.2.1. Construire un nom de fichier de sortie dynamique

Dans le bloc process, faites les modifications de code suivantes :

=== "Après"

    ```groovy title="hello-channels.nf" linenums="6" hl_lines="7 11"
    process sayHello {

        input:
        val greeting

        output:
        path "${greeting}-output.txt"

        script:
        """
        echo '${greeting}' > '${greeting}-output.txt'
        """
    }
    ```

=== "Avant"

    ```groovy title="hello-channels.nf" linenums="6" hl_lines="7 11"
    process sayHello {

        input:
        val greeting

        output:
        path 'output.txt'

        script:
        """
        echo '${greeting}' > output.txt
        """
    }
    ```

Assurez-vous de remplacer `output.txt` à la fois dans la définition de sortie et dans le bloc de commande `script:`.

!!! tip "Astuce"

    Dans la définition de sortie, vous DEVEZ utiliser des guillemets doubles autour de l'expression du nom de fichier de sortie (PAS des guillemets simples), sinon cela échouera.

Cela devrait produire un nom de fichier de sortie unique chaque fois que le processus est appelé, afin qu'il puisse être distingué des sorties d'autres appels au même processus dans le répertoire de sortie.

#### 2.2.2. Exécuter le workflow

Exécutons-le. Notez que nous revenons à l'exécution avec les paramètres de log ANSI par défaut.

```bash
nextflow run hello-channels.nf
```

??? success "Sortie de la commande"

    ```console
     N E X T F L O W   ~  version 25.10.2

    Launching `hello-channels.nf` [sharp_minsky] DSL2 - revision: 16a291febe

    executor >  local (3)
    [e8/33ee64] sayHello (2) [100%] 3 of 3 ✔
    Hello
    Bonjour
    Holà
    ```

En revenant à la vue résumée, la sortie est à nouveau résumée sur une ligne.
Jetez un œil au répertoire `results` pour voir si toutes les salutations de sortie y sont.

??? abstract "Contenu du répertoire"

    ```console
    results/hello_channels/
    ├── Bonjour-output.txt
    ├── Hello-output.txt
    ├── Holà-output.txt
    └── output.txt
    ```

Oui ! Et chacune a le contenu attendu.

??? abstract "Contenu des fichiers"

    ```console title="Bonjour-output.txt"
    Bonjour
    ```

    ```console title="Hello-output.txt"
    Hello
    ```

    ```console title="Holà-output.txt"
    Holà
    ```

Succès ! Maintenant nous pouvons ajouter autant de salutations que nous voulons sans nous soucier de l'écrasement des fichiers de sortie.

!!! tip "Astuce"

    En pratique, nommer les fichiers en fonction des données d'entrée elles-mêmes est presque toujours peu pratique.
    La meilleure façon de générer des noms de fichiers dynamiques est de passer des métadonnées à un processus avec les fichiers d'entrée.
    Les métadonnées sont généralement fournies via une « feuille d'échantillons » ou équivalents.
    Vous apprendrez comment faire cela plus tard dans votre formation Nextflow (voir la [quête annexe sur les métadonnées](../side_quests/metadata.md)).

### À retenir

Vous savez comment alimenter plusieurs éléments d'entrée à travers un canal.

### Et ensuite ?

Apprendre à utiliser un opérateur pour transformer le contenu d'un canal.

---

## 3. Fournir plusieurs entrées via un tableau

Nous venons de vous montrer comment gérer plusieurs éléments d'entrée qui étaient codés en dur directement dans la fabrique de canaux.
Et si nous voulions fournir ces entrées multiples d'une manière différente ?

Par exemple, imaginez que nous configurions une variable d'entrée contenant un tableau d'éléments comme ceci :

`greetings_array = ['Hello','Bonjour','Holà']`

Pouvons-nous charger cela dans notre canal de sortie et nous attendre à ce que ça fonctionne ?

<figure class="excalidraw">
--8<-- "docs/en/docs/hello_nextflow/img/hello-pipeline-multi-inputs-array.svg"
</figure>

Découvrons-le.

### 3.1. Fournir un tableau de valeurs comme entrée au canal

Le bon sens suggère que nous devrions pouvoir simplement passer un tableau de valeurs au lieu d'une seule valeur.
Essayons ; nous devrons configurer la variable d'entrée et la charger dans la fabrique de canaux.

#### 3.1.1. Configurer la variable d'entrée

Prenons la variable `greetings_array` que nous venons d'imaginer et rendons-la réelle en l'ajoutant au bloc workflow :

=== "Après"

    ```groovy title="hello-channels.nf" linenums="27" hl_lines="4 5"
    workflow {

        main:
        // déclarer un tableau de salutations d'entrée
        greetings_array = ['Hello','Bonjour','Holà']
        // créer un canal pour les entrées
        greeting_ch = channel.of('Hello','Bonjour','Holà')
                             .view()
        // émettre une salutation
        sayHello(greeting_ch)

        publish:
        first_output = sayHello.out
    }
    ```

=== "Avant"

    ```groovy title="hello-channels.nf" linenums="27"
    workflow {

        main:
        // créer un canal pour les entrées
        greeting_ch = channel.of('Hello','Bonjour','Holà')
                             .view()
        // émettre une salutation
        sayHello(greeting_ch)

        publish:
        first_output = sayHello.out
    }
    ```

Ce n'est pas encore fonctionnel, nous avons juste ajouté une déclaration pour le tableau.

#### 3.1.2. Définir le tableau de salutations comme entrée de la fabrique de canaux

Maintenant nous allons remplacer les valeurs `'Hello','Bonjour','Holà'` actuellement codées en dur dans la fabrique de canaux par le `greetings_array` que nous venons de créer.

Dans le bloc workflow, faites la modification suivante :

=== "Après"

    ```groovy title="hello-channels.nf" linenums="27" hl_lines="7"
    workflow {

        main:
        // déclarer un tableau de salutations d'entrée
        greetings_array = ['Hello','Bonjour','Holà']
        // créer un canal pour les entrées
        greeting_ch = channel.of(greetings_array)
                             .view()
        // émettre une salutation
        sayHello(greeting_ch)

        publish:
        first_output = sayHello.out
    }
    ```

=== "Avant"

    ```groovy title="hello-channels.nf" linenums="27" hl_lines="7"
    workflow {

        main:
        // déclarer un tableau de salutations d'entrée
        greetings_array = ['Hello','Bonjour','Holà']
        // créer un canal pour les entrées
        greeting_ch = channel.of('Hello','Bonjour','Holà')
                             .view()
        // émettre une salutation
        sayHello(greeting_ch)

        publish:
        first_output = sayHello.out
    }
    ```

Cela devrait être fonctionnel maintenant.

#### 3.1.3. Exécuter le workflow

Essayons de l'exécuter :

```bash
nextflow run hello-channels.nf
```

??? failure "Sortie de la commande"

    ```console hl_lines="7 11 16"
    N E X T F L O W   ~  version 25.10.2

    Launching `hello-channels.nf` [friendly_koch] DSL2 - revision: 97256837a7

    executor >  local (1)
    [a8/1f6ead] sayHello (1) | 0 of 1
    [Hello, Bonjour, Holà]
    ERROR ~ Error executing process > 'sayHello (1)'

    Caused by:
      Missing output file(s) `[Hello, Bonjour, Holà]-output.txt` expected by process `sayHello (1)`


    Command executed:

      echo '[Hello, Bonjour, Holà]' > '[Hello, Bonjour, Holà]-output.txt'

    Command exit status:
      0

    Command output:
      (empty)

    Work dir:
      /workspaces/training/hello-nextflow/work/a8/1f6ead5f3fa30a3c508e2e7cf83ffb

    Tip: you can replicate the issue by changing to the process work dir and entering the command `bash .command.run`

    -- Check '.nextflow.log' file for details
    ```

Oh non ! Il y a une erreur !

Regardez la sortie de `view()` et les messages d'erreur.

On dirait que Nextflow a essayé d'exécuter un seul appel de processus, en utilisant `[Hello, Bonjour, Holà]` comme valeur de chaîne, au lieu d'utiliser les trois chaînes dans le tableau comme valeurs séparées.

<figure class="excalidraw">
--8<-- "docs/en/docs/hello_nextflow/img/hello-channels-array-fail.svg"
</figure>

Donc c'est l'« emballage » qui cause le problème.
Comment faire pour que Nextflow déballe le tableau et charge les chaînes individuelles dans le canal ?

### 3.2. Utiliser un opérateur pour transformer le contenu du canal

C'est là que les [**opérateurs**](https://nextflow.io/docs/latest/reference/operator.html) entrent en jeu.
Vous avez déjà utilisé l'opérateur `.view()`, qui regarde juste ce qu'il y a à l'intérieur.
Maintenant nous allons regarder les opérateurs qui nous permettent d'agir sur le contenu d'un canal.

Si vous parcourez la [liste des opérateurs](https://nextflow.io/docs/latest/reference/operator.html) dans la documentation Nextflow, vous trouverez [`flatten()`](https://nextflow.io/docs/latest/reference/operator.html#flatten), qui fait exactement ce dont nous avons besoin : déballer le contenu d'un tableau et les émettre comme éléments individuels.

#### 3.2.1. Ajouter l'opérateur `flatten()`

Pour appliquer l'opérateur `flatten()` à notre canal d'entrée, nous l'ajoutons à la déclaration de la fabrique de canaux.

Dans le bloc workflow, faites la modification de code suivante :

=== "Après"

    ```groovy title="hello-channels.nf" linenums="27" hl_lines="9"
    workflow {

        main:
        // déclarer un tableau de salutations d'entrée
        greetings_array = ['Hello','Bonjour','Holà']
        // créer un canal pour les entrées
        greeting_ch = channel.of(greetings_array)
                             .view()
                             .flatten()
        // émettre une salutation
        sayHello(greeting_ch)

        publish:
        first_output = sayHello.out
    }
    ```

=== "Avant"

    ```groovy title="hello-channels.nf" linenums="27"
    workflow {

        main:
        // déclarer un tableau de salutations d'entrée
        greetings_array = ['Hello','Bonjour','Holà']
        // créer un canal pour les entrées
        greeting_ch = channel.of(greetings_array)
                             .view()
        // émettre une salutation
        sayHello(greeting_ch)

        publish:
        first_output = sayHello.out
    }
    ```

Ici nous avons ajouté l'opérateur sur la ligne suivante pour la lisibilité, mais vous pouvez ajouter des opérateurs sur la même ligne que la fabrique de canaux si vous préférez, comme ceci :
`greeting_ch = channel.of(greetings_array).view().flatten()`

<figure class="excalidraw">
--8<-- "docs/en/docs/hello_nextflow/img/hello-channels-array-success.svg"
</figure>

#### 3.2.2. Affiner les instructions `view()`

Nous pourrions exécuter cela tout de suite pour tester si ça fonctionne, mais pendant qu'on y est, nous allons affiner la façon dont nous inspectons le contenu du canal.

Nous voulons pouvoir contraster ce à quoi ressemble le contenu avant et après l'application de l'opérateur `flatten()`, donc nous allons en ajouter un second, ET nous allons ajouter un peu de code pour les étiqueter plus clairement dans la sortie.

Dans le bloc workflow, faites la modification de code suivante :

=== "Après"

    ```groovy title="hello-channels.nf" linenums="27" hl_lines="8-10"
    workflow {

        main:
        // déclarer un tableau de salutations d'entrée
        greetings_array = ['Hello','Bonjour','Holà']
        // créer un canal pour les entrées
        greeting_ch = channel.of(greetings_array)
                             .view { greeting -> "Before flatten: $greeting" }
                             .flatten()
                             .view { greeting -> "After flatten: $greeting" }
        // émettre une salutation
        sayHello(greeting_ch)

        publish:
        first_output = sayHello.out
    }
    ```

=== "Avant"

    ```groovy title="hello-channels.nf" linenums="27" hl_lines="8-9"
    workflow {

        main:
        // déclarer un tableau de salutations d'entrée
        greetings_array = ['Hello','Bonjour','Holà']
        // créer un canal pour les entrées
        greeting_ch = channel.of(greetings_array)
                             .view()
                             .flatten()
        // émettre une salutation
        sayHello(greeting_ch)

        publish:
        first_output = sayHello.out
    }
    ```

Vous voyez que nous avons ajouté une deuxième instruction `.view`, et pour chacune d'elles, nous avons remplacé les parenthèses vides (`()`) par des accolades contenant du code, tel que `{ greeting -> "Before flatten: $greeting" }`.

Ce sont des _closures_. Le code qu'elles contiennent sera exécuté pour chaque élément dans le canal.
Nous définissons une variable temporaire pour la valeur interne, ici appelée `greeting` (mais ce pourrait être n'importe quel nom arbitraire), qui n'est utilisée que dans la portée de cette closure.

Dans cet exemple, `$greeting` représente chaque élément individuel chargé dans le canal.
Cela donnera une sortie console bien étiquetée.

!!! info "Information"

    Dans certains pipelines, vous pourriez voir une variable spéciale appelée `$it` utilisée à l'intérieur des closures d'opérateurs.
    C'est une variable _implicite_ qui permet un accès raccourci à la variable interne,
    sans avoir besoin de la définir avec un `->`.

    Nous préférons être explicites pour aider à la clarté du code, ainsi la syntaxe `$it` est découragée et sera progressivement supprimée du langage Nextflow.

#### 3.2.3. Exécuter le workflow

Finalement, vous pouvez essayer d'exécuter à nouveau le workflow !

```bash
nextflow run hello-channels.nf
```

??? success "Sortie de la commande"

    ```console hl_lines="7-10"
     N E X T F L O W   ~  version 25.10.2

    Launching `hello-channels.nf` [sleepy_gutenberg] DSL2 - revision: 1db4f760ee

    executor >  local (3)
    [b1/6a1e15] sayHello (2) [100%] 3 of 3 ✔
    Before flatten: [Hello, Bonjour, Holà]
    After flatten: Hello
    After flatten: Bonjour
    After flatten: Holà
    ```

Cette fois ça fonctionne ET nous donne l'aperçu supplémentaire de ce à quoi ressemble le contenu du canal avant et après l'exécution de l'opérateur `flatten()`.

- Une seule instruction `Before flatten:` parce qu'à ce moment le canal contient un élément, le tableau original.
- Trois instructions `After flatten:` séparées, une pour chaque salutation, qui sont maintenant des éléments individuels dans le canal.

De manière importante, cela signifie que chaque élément peut maintenant être traité séparément par le workflow.

!!! tip "Astuce"

    Il est techniquement possible d'obtenir les mêmes résultats en utilisant une fabrique de canaux différente, [`channel.fromList`](https://nextflow.io/docs/latest/reference/channel.html#fromlist), qui inclut une étape de mapping implicite dans son fonctionnement.
    Ici nous avons choisi de ne pas l'utiliser afin de démontrer l'utilisation d'un opérateur sur un cas d'utilisation simple.

### À retenir

Vous savez comment utiliser un opérateur comme `flatten()` pour transformer le contenu d'un canal, et comment utiliser l'opérateur `view()` pour inspecter le contenu du canal avant et après l'application d'un opérateur.

### Et ensuite ?

Apprendre à faire prendre au workflow un fichier comme source de valeurs d'entrée.

---

## 4. Lire les valeurs d'entrée depuis un fichier CSV

De façon réaliste, nous allons rarement, voire jamais, commencer à partir d'un tableau de valeurs.
Le plus probable est que nous aurons un ou plusieurs fichiers contenant les données qui doivent être traitées, dans un format structuré quelconque.

Nous avons préparé un fichier CSV appelé `greetings.csv` qui contient plusieurs salutations d'entrée, imitant le type de données tabulaires que vous pourriez vouloir traiter dans une analyse de données réelle, stocké sous `data/`.
(Les nombres ne sont pas significatifs, ils sont juste là à des fins d'illustration.)

```csv title="data/greetings.csv" linenums="1"
Hello,English,123
Bonjour,French,456
Holà,Spanish,789
```

Notre prochaine tâche est d'adapter notre workflow pour lire les valeurs de ce fichier.

<figure class="excalidraw">
--8<-- "docs/en/docs/hello_nextflow/img/hello-pipeline-multi-inputs-csv.svg"
</figure>

Voyons comment nous pouvons faire cela.

### 4.1. Modifier le script pour attendre un fichier CSV comme source de salutations

Pour commencer, nous allons devoir faire deux changements clés au script :

- Changer le paramètre d'entrée pour pointer vers le fichier CSV
- Changer la fabrique de canaux pour une conçue pour gérer un fichier

#### 4.1.1. Changer le paramètre d'entrée pour pointer vers le fichier CSV

Vous vous souvenez du paramètre `params.input` que nous avons configuré dans la Partie 1 ?
Nous allons le mettre à jour pour pointer vers le fichier CSV contenant nos salutations.

Faites la modification suivante à la déclaration du paramètre :

=== "Après"

    ```groovy title="hello-channels.nf" linenums="20" hl_lines="5"
    /*
    * Paramètres du pipeline
    */
    params {
        input: Path = 'data/greetings.csv'
    }
    ```

=== "Avant"

    ```groovy title="hello-channels.nf" linenums="20" hl_lines="5"
    /*
     * Paramètres du pipeline
     */
    input: String = 'Holà mundo!'
    ```

Cela suppose que le fichier est co-localisé avec le code du workflow.
Vous apprendrez comment gérer d'autres emplacements de données plus tard dans votre parcours Nextflow.

#### 4.1.2. Passer à une fabrique de canaux conçue pour gérer un fichier

Puisque nous voulons maintenant utiliser un fichier au lieu de simples chaînes comme entrée, nous ne pouvons pas utiliser la fabrique de canaux `channel.of()` d'avant.
Nous devons passer à l'utilisation d'une nouvelle fabrique de canaux, [`channel.fromPath()`](https://www.nextflow.io/docs/latest/reference/channel.html#frompath), qui a une fonctionnalité intégrée pour gérer les chemins de fichiers.

Dans le bloc workflow, faites la modification de code suivante :

=== "Après"

    ```groovy title="hello-channels.nf" linenums="27" hl_lines="4-8"
    workflow {

        main:
        // créer un canal pour les entrées depuis un fichier CSV
        greeting_ch = channel.fromPath(params.input)
                             .view { greeting -> "Before flatten: $greeting" }
                             // .flatten()
                             // .view { greeting -> "After flatten: $greeting" }
        // émettre une salutation
        sayHello(greeting_ch)

        publish:
        first_output = sayHello.out
    }
    ```

=== "Avant"

    ```groovy title="hello-channels.nf" linenums="27" hl_lines="4-8"
    workflow {

        main:
        // déclarer un tableau de salutations d'entrée
        greetings_array = ['Hello','Bonjour','Holà']
        // créer un canal pour les entrées
        greeting_ch = channel.of(greetings_array)
                             .view { greeting -> "Before flatten: $greeting" }
                             .flatten()
                             .view { greeting -> "After flatten: $greeting" }
        // émettre une salutation
        sayHello(greeting_ch)

        publish:
        first_output = sayHello.out
    }
    ```

Vous remarquerez que nous avons rechangé l'entrée du canal en `param.input`, et supprimé la déclaration `greetings_array` puisque nous n'en aurons plus besoin.
Nous avons également commenté le `flatten()` et la deuxième instruction `view()`.

#### 4.1.3. Exécuter le workflow

Essayons d'exécuter le workflow avec la nouvelle fabrique de canaux et le fichier d'entrée.

```bash
nextflow run hello-channels.nf
```

??? failure "Sortie de la commande"

    ```console hl_lines="5 6 9 14"
     N E X T F L O W   ~  version 25.10.2

    Launching `hello-channels.nf` [peaceful_poisson] DSL2 - revision: a286c08ad5

    [-        ] sayHello [  0%] 0 of 1
    Before flatten: /workspaces/training/hello-nextflow/data/greetings.csv
    ERROR ~ Error executing process > 'sayHello (1)'

    Caused by:
      File `/workspaces/training/hello-nextflow/data/greetings.csv-output.txt` is outside the scope of the process work directory: /workspaces/training/hello-nextflow/work/30/e610cb4ea5ae8693f456ac3329c92f


    Command executed:

      echo '/workspaces/training/hello-nextflow/data/greetings.csv' > '/workspaces/training/hello-nextflow/data/greetings.csv-output.txt'

    Command exit status:
      -

    Command output:
      (empty)

    Work dir:
      /workspaces/training/hello-nextflow/work/30/e610cb4ea5ae8693f456ac3329c92f

    Tip: when you have fixed the problem you can continue the execution adding the option `-resume` to the run command line

    -- Check '.nextflow.log' file for details
    ```

Oh non, ça ne fonctionne pas. Regardez le début de la sortie console et le message d'erreur.
La partie `Command executed:` est particulièrement utile ici.

Cela peut sembler un peu familier.
On dirait que Nextflow a essayé d'exécuter un seul appel de processus en utilisant le chemin du fichier lui-même comme valeur de chaîne.
Donc il a résolu le chemin du fichier correctement, mais il n'a pas réellement analysé son contenu, ce qui est ce que nous voulions.

Comment faire pour que Nextflow ouvre le fichier et charge son contenu dans le canal ?

On dirait qu'on a besoin d'un autre [opérateur](https://nextflow.io/docs/latest/reference/operator.html) !

### 4.2. Utiliser l'opérateur `splitCsv()` pour analyser le fichier

En parcourant à nouveau la liste des opérateurs, nous trouvons [`splitCsv()`](https://nextflow.io/docs/latest/reference/operator.html#splitcsv), qui est conçu pour analyser et diviser le texte au format CSV.

#### 4.2.1. Appliquer `splitCsv()` au canal

Pour appliquer l'opérateur, nous l'ajoutons à la ligne de la fabrique de canaux comme précédemment.

Dans le bloc workflow, faites la modification de code suivante pour remplacer `flatten()` par `splitcsv()` (décommenté) :

=== "Après"

    ```groovy title="hello-channels.nf" linenums="27" hl_lines="6-8"
    workflow {

        main:
        // créer un canal pour les entrées depuis un fichier CSV
        greeting_ch = channel.fromPath(params.input)
                             .view { csv -> "Before splitCsv: $csv" }
                             .splitCsv()
                             .view { csv -> "After splitCsv: $csv" }
        // émettre une salutation
        sayHello(greeting_ch)

        publish:
        first_output = sayHello.out
    }
    ```

=== "Avant"

    ```groovy title="hello-channels.nf" linenums="27" hl_lines="6-8"
    workflow {

        main:
        // créer un canal pour les entrées depuis un fichier CSV
        greeting_ch = channel.fromPath(params.input)
                             .view { greeting -> "Before flatten: $greeting" }
                             // .flatten()
                             // .view { greeting -> "After flatten: $greeting" }
        // émettre une salutation
        sayHello(greeting_ch)

        publish:
        first_output = sayHello.out
    }
    ```

Comme vous pouvez le voir, nous avons également mis à jour les instructions `view()` avant/après.
Techniquement nous aurions pu utiliser le même nom de variable (`greeting`) mais nous l'avons mis à jour en quelque chose de plus approprié (`csv`) pour rendre le code plus lisible par les autres.

#### 4.2.2. Exécuter à nouveau le workflow

Essayons d'exécuter le workflow avec la logique d'analyse CSV ajoutée.

```bash
nextflow run hello-channels.nf
```

??? failure "Sortie de la commande"

    ```console hl_lines="7-11 14 19"
     N E X T F L O W   ~  version 25.10.2

    Launching `hello-channels.nf` [insane_fermat] DSL2 - revision: 8e62fcbeb1

    executor >  local (3)
    [24/76da2f] sayHello (2) [  0%] 0 of 3 ✘
    Before splitCsv: /workspaces/training/hello-nextflow/data/greetings.csv
    After splitCsv: [Hello, English, 123]
    After splitCsv: [Bonjour, French, 456]
    After splitCsv: [Holà, Spanish, 789]
    ERROR ~ Error executing process > 'sayHello (2)'

    Caused by:
      Missing output file(s) `[Bonjour, French, 456]-output.txt` expected by process `sayHello (2)`


    Command executed:

      echo '[Bonjour, French, 456]' > '[Bonjour, French, 456]-output.txt'

    Command exit status:
      0

    Command output:
      (empty)

    Work dir:
      /workspaces/training/hello-nextflow/work/24/76da2fcc4876b61632749f99e26a50

    Tip: you can try to figure out what's wrong by changing to the process work dir and showing the script file named `.command.sh`

    -- Check '.nextflow.log' file for details
    ```

Intéressant, cela échoue aussi, mais avec une erreur différente.
Cette fois Nextflow a analysé le contenu du fichier (youpi !) mais il a chargé chaque ligne comme un tableau, et chaque tableau est un élément dans le canal.

<figure class="excalidraw">
--8<-- "docs/en/docs/hello_nextflow/img/hello-channels-split-fail.svg"
</figure>

Nous devons lui dire de ne prendre que la première colonne de chaque ligne.
Alors comment déballer cela ?

Nous avons précédemment utilisé `flatten()` pour déballer le contenu d'un canal, mais cela ne fonctionnerait pas ici parce que flatten déballe _tout_ (n'hésitez pas à l'essayer si vous voulez voir par vous-même).

À la place, nous utiliserons un autre opérateur appelé `map()` qui est vraiment utile et apparaît beaucoup dans les pipelines Nextflow.

### 4.3. Utiliser l'opérateur `map()` pour extraire les salutations

L'opérateur [`map()`](https://nextflow.io/docs/latest/reference/operator.html#map) est un petit outil très pratique qui nous permet de faire toutes sortes de mappings sur le contenu d'un canal.

Dans ce cas, nous allons l'utiliser pour extraire cet unique élément que nous voulons de chaque ligne dans notre fichier de données.
Voici à quoi ressemble la syntaxe :

```groovy title="Syntaxe"
.map { row -> row[0] }
```

Cela signifie « pour chaque ligne dans le canal, prendre le 0ème (premier) élément qu'elle contient ».

Appliquons donc cela à notre analyse CSV.

#### 4.3.1. Appliquer `map()` au canal

Dans le bloc workflow, faites la modification de code suivante :

=== "Après"

    ```groovy title="hello-channels.nf" linenums="27" hl_lines="9 10"
    workflow {

        main:
        // créer un canal pour les entrées depuis un fichier CSV
        greeting_ch = channel.fromPath(params.input)
                             .view { csv -> "Before splitCsv: $csv" }
                             .splitCsv()
                             .view { csv -> "After splitCsv: $csv" }
                             .map { item -> item[0] }
                             .view { csv -> "After map: $csv" }
        // émettre une salutation
        sayHello(greeting_ch)

        publish:
        first_output = sayHello.out
    }
    ```

=== "Avant"

    ```groovy title="hello-channels.nf" linenums="27"
    workflow {

        main:
        // créer un canal pour les entrées depuis un fichier CSV
        greeting_ch = channel.fromPath(params.input)
                             .view { csv -> "Before splitCsv: $csv" }
                             .splitCsv()
                             .view { csv -> "After splitCsv: $csv" }
        // émettre une salutation
        sayHello(greeting_ch)

        publish:
        first_output = sayHello.out
    }
    ```

Vous voyez que nous avons ajouté un autre appel `view()` pour confirmer que l'opérateur fait ce que nous attendons.

#### 4.3.2. Exécuter le workflow

Exécutons cela une fois de plus :

```bash
nextflow run hello-channels.nf
```

??? success "Sortie de la commande"

    ```console
     N E X T F L O W   ~  version 25.10.2

    Launching `hello-channels.nf` [focused_volhard] DSL2 - revision: de435e45be

    executor >  local (3)
    [54/6eebe3] sayHello (3) [100%] 3 of 3 ✔
    Before splitCsv: /workspaces/training/hello-nextflow/data/greetings.csv
    After splitCsv: [Hello, English, 123]
    After splitCsv: [Bonjour, French, 456]
    After splitCsv: [Holà, Spanish, 789]
    After map: Hello
    After map: Bonjour
    After map: Holà
    ```

Cette fois cela devrait s'exécuter sans erreur.

En regardant la sortie des instructions `view()`, vous voyez ce qui suit :

- Une seule instruction `Before splitCsv:` : à ce moment le canal contient un élément, le chemin du fichier original.
- Trois instructions `After splitCsv:` séparées : une pour chaque salutation, mais chacune est contenue dans un tableau qui correspond à cette ligne dans le fichier.
- Trois instructions `After map:` séparées : une pour chaque salutation, qui sont maintenant des éléments individuels dans le canal.

_Notez que les lignes peuvent apparaître dans un ordre différent dans votre sortie._

<figure class="excalidraw">
--8<-- "docs/en/docs/hello_nextflow/img/hello-channels-split-and-map.svg"
</figure>

Vous pouvez également regarder les fichiers de sortie pour vérifier que chaque salutation a été correctement extraite et traitée à travers le workflow.

Nous avons obtenu le même résultat qu'avant, mais maintenant nous avons beaucoup plus de flexibilité pour ajouter plus d'éléments au canal de salutations que nous voulons traiter en modifiant un fichier d'entrée, sans modifier aucun code.
Vous apprendrez des approches plus sophistiquées pour gérer des entrées complexes dans une formation ultérieure.

### À retenir

Vous savez comment utiliser le constructeur de canal `.fromPath()` et les opérateurs `splitCsv()` et `map()` pour lire un fichier de valeurs d'entrée et les gérer de manière appropriée.

Plus généralement, vous avez une compréhension de base de la façon dont Nextflow utilise les **canaux** pour gérer les entrées vers les processus et les **opérateurs** pour transformer leur contenu.
Vous avez également vu comment les canaux gèrent l'exécution parallèle de manière implicite.

<figure class="excalidraw">
--8<-- "docs/en/docs/hello_nextflow/img/hello-channels-parallel.svg"
</figure>

### Et ensuite ?

Prenez une grande pause, vous avez travaillé dur sur celui-ci !

Quand vous êtes prêt·e, passez à la [**Partie 3 : Hello Workflow**](./03_hello_workflow.md) pour apprendre comment ajouter plus d'étapes et les connecter ensemble dans un vrai workflow.

---

## Quiz

<quiz>
Qu'est-ce qu'un canal dans Nextflow ?
- [ ] Une spécification de chemin de fichier
- [ ] Une définition de processus
- [x] Une structure de type file d'attente pour passer des données entre processus
- [ ] Un paramètre de configuration

En savoir plus : [1.1. Créer un canal d'entrée](#11-creer-un-canal-dentree)
</quiz>

<quiz>
Que va produire ce code ?

```groovy
channel.of('Hello', 'Bonjour', 'Hola')
    .view()
```

- [ ] `['Hello', 'Bonjour', 'Hola']` (une seule liste)
- [x] Chaque élément sur une ligne séparée : `Hello`, `Bonjour`, `Hola`
- [ ] Rien (les canaux n'affichent pas par défaut)
- [ ] Une erreur (syntaxe invalide)

En savoir plus : [1.1. Créer un canal d'entrée](#11-creer-un-canal-dentree)
</quiz>

<quiz>
Quand un canal contient plusieurs valeurs, comment Nextflow gère-t-il l'exécution du processus ?
- [ ] Le processus s'exécute une fois avec toutes les valeurs
- [x] Le processus s'exécute une fois pour chaque valeur dans le canal
- [ ] Le processus s'exécute seulement avec la première valeur
- [ ] Le processus s'exécute seulement avec la dernière valeur

En savoir plus : [2. Modifier le workflow pour s'exécuter sur plusieurs valeurs d'entrée](#2-modifier-le-flux-de-travail-pour-sexecuter-sur-plusieurs-valeurs-dentree)
</quiz>

<quiz>
Que fait l'opérateur `flatten()` ?
- [ ] Combine plusieurs canaux en un seul
- [ ] Trie les éléments du canal
- [x] Déballe les tableaux en éléments individuels
- [ ] Supprime les éléments dupliqués

En savoir plus : [3.2.1. Ajouter l'opérateur `flatten()`](#321-ajouter-loperateur-flatten)
</quiz>

<quiz>
Quel est le but de l'opérateur `view()` ?
- [ ] Filtrer le contenu du canal
- [ ] Transformer les éléments du canal
- [x] Inspecter et déboguer le contenu du canal
- [ ] Sauvegarder le contenu du canal dans un fichier

En savoir plus : [1.4. Utiliser `view()` pour inspecter le contenu du canal](#14-utiliser-view-pour-inspecter-le-contenu-du-canal)
</quiz>

<quiz>
Que fait `splitCsv()` ?
- [ ] Crée un fichier CSV à partir du contenu du canal
- [ ] Divise une chaîne par des virgules
- [x] Analyse un fichier CSV en tableaux représentant chaque ligne
- [ ] Fusionne plusieurs fichiers CSV

En savoir plus : [4.2. Utiliser l'opérateur `splitCsv()` pour analyser le fichier](#42-utiliser-loperateur-splitcsv-pour-analyser-le-fichier)
</quiz>

<quiz>
Quel est le but de l'opérateur `map()` ?
- [ ] Filtrer les éléments d'un canal
- [ ] Combiner plusieurs canaux
- [x] Transformer chaque élément dans un canal
- [ ] Compter les éléments dans un canal

En savoir plus : [4.3. Utiliser l'opérateur `map()` pour extraire les salutations](#43-utiliser-loperateur-map-pour-extraire-les-salutations)
</quiz>

<quiz>
Pourquoi est-il important d'utiliser des noms de fichiers de sortie dynamiques lors du traitement de plusieurs entrées ?
- [ ] Pour améliorer les performances
- [ ] Pour réduire l'espace disque
- [x] Pour empêcher les fichiers de sortie de s'écraser mutuellement
- [ ] Pour activer la fonctionnalité de reprise

En savoir plus : [2.2. S'assurer que les noms de fichiers de sortie seront uniques](#22-sassurer-que-les-noms-de-fichiers-de-sortie-seront-uniques)
</quiz>
