# Partie 5 : Hello Containers

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Traduction assistée par IA - [en savoir plus et suggérer des améliorations](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

<div class="video-wrapper">
  <iframe width="560" height="315" src="https://www.youtube.com/embed/Xqr--bKEN9U?si=y8lAedhEHWaTV4zd&amp;list=PLPZ8WHdZGxmWKozQuzr27jyMGqp9kElVK&amp;cc_load_policy=1" title="YouTube video player" frameborder="0" allow="accelerometer; autoplay; clipboard-write; encrypted-media; gyroscope; picture-in-picture; web-share" referrerpolicy="strict-origin-when-cross-origin" allowfullscreen></iframe>
</div>

/// caption
:fontawesome-brands-youtube:{ .youtube } Consultez [la playlist complète](https://youtube.com/playlist?list=PLPZ8WHdZGxmWKozQuzr27jyMGqp9kElVK&si=eF7cLR62goy-lc6n) sur la chaîne YouTube Nextflow.

:green_book: La transcription vidéo est disponible [ici](./transcripts/05_hello_containers.md).
///

Dans les parties 1 à 4 de cette formation, vous avez appris à utiliser les éléments de base de Nextflow pour assembler un workflow simple capable de traiter du texte, de paralléliser l'exécution en présence de plusieurs entrées, et de collecter les résultats pour un traitement ultérieur.

Cependant, vous étiez limité·e aux outils UNIX de base disponibles dans votre environnement.
Les tâches réelles nécessitent souvent divers outils et paquets non inclus par défaut.
Typiquement, vous devriez installer ces outils, gérer leurs dépendances et résoudre les éventuels conflits.

Tout cela est très fastidieux et ennuyeux, nous allons donc vous montrer comment utiliser des **conteneurs** pour résoudre ce problème de manière beaucoup plus pratique.

Un **conteneur** est une unité logicielle légère, autonome et exécutable créée à partir d'une **image** de conteneur qui inclut tout le nécessaire pour exécuter une application, y compris le code, les bibliothèques système et les paramètres.
Comme vous pouvez l'imaginer, cela va être très utile pour rendre vos pipelines plus reproductibles.

Notez que nous enseignons cela en utilisant [Docker](https://www.docker.com/get-started/), mais gardez à l'esprit que Nextflow prend en charge [plusieurs autres technologies de conteneurs](https://nextflow.io/docs/latest/container.html) également.

??? info "Comment commencer à partir de cette section"

    Cette section de la formation suppose que vous avez terminé les parties 1 à 4 de la formation [Hello Nextflow](./index.md) et que vous disposez d'un pipeline complet et fonctionnel.

    Si vous commencez la formation à partir de ce point, vous devrez copier le répertoire `modules` depuis les solutions :

    ```bash
    cp -r solutions/4-hello-modules/modules .
    ```

---

## 0. Échauffement : Exécuter `hello-containers.nf`

Nous allons utiliser le script de workflow `hello-containers.nf` comme point de départ.
Il est équivalent au script produit en suivant la partie 4 de cette formation, sauf que nous avons modifié les destinations de sortie :

```groovy title="hello-containers.nf" linenums="37" hl_lines="3 7 11 15"
output {
    first_output {
        path 'hello_containers'
        mode 'copy'
    }
    uppercased {
        path 'hello_containers'
        mode 'copy'
    }
    collected {
        path 'hello_containers'
        mode 'copy'
    }
    batch_report {
        path 'hello_containers'
        mode 'copy'
    }
}
```

Juste pour vérifier que tout fonctionne, exécutez le script une fois avant d'apporter des modifications :

```bash
nextflow run hello-containers.nf
```

??? success "Sortie de la commande"

    ```console
     N E X T F L O W   ~  version 25.10.2

    Launching `hello-containers.nf` [nice_escher] DSL2 - revision: d5dfdc9872

    executor > local (7)
    [5a/ec1fa1] sayHello (2) [100%] 3 of 3 ✔
    [30/32b5b8] convertToUpper (3) [100%] 3 of 3 ✔
    [d3/be01bc] collectGreetings [100%] 1 of 1 ✔

    ```

Comme précédemment, vous trouverez les fichiers de sortie dans le répertoire spécifié dans le bloc `output` (`results/hello_containers/`).

??? abstract "Contenu du répertoire"

    ```console
    results/hello_containers/
    ├── Bonjour-output.txt
    ├── COLLECTED-batch-output.txt
    ├── Hello-output.txt
    ├── Holà-output.txt
    ├── batch-report.txt
    ├── UPPER-Bonjour-output.txt
    ├── UPPER-Hello-output.txt
    └── UPPER-Holà-output.txt
    ```

Si cela a fonctionné pour vous, vous êtes prêt·e à apprendre comment utiliser les conteneurs.

---

## 1. Utiliser un conteneur 'manuellement'

Ce que nous voulons faire, c'est ajouter une étape à notre workflow qui utilisera un conteneur pour l'exécution.

Cependant, nous allons d'abord passer en revue quelques concepts et opérations de base pour consolider votre compréhension de ce que sont les conteneurs avant de commencer à les utiliser dans Nextflow.

### 1.1. Télécharger l'image de conteneur

Pour utiliser un conteneur, vous téléchargez généralement ou _tirez_ une image de conteneur depuis un registre de conteneurs, puis vous exécutez l'image de conteneur pour créer une instance de conteneur.

La syntaxe générale est la suivante :

```bash title="Syntaxe"
docker pull '<container>'
```

La partie `docker pull` est l'instruction au système de conteneurs pour tirer une image de conteneur depuis un dépôt.

La partie `'<container>'` est l'adresse URI de l'image de conteneur.

À titre d'exemple, tirons une image de conteneur qui contient [cowpy](https://github.com/jeffbuttars/cowpy), une implémentation Python d'un outil appelé `cowsay` qui génère de l'art ASCII pour afficher des entrées de texte arbitraires de manière amusante.

```txt title="Exemple"
 ________________________
< Are we having fun yet? >
 ------------------------
    \                                  ___-------___
     \                             _-~~             ~~-_
      \                         _-~                    /~-_
             /^\__/^\         /~  \                   /    \
           /|  O|| O|        /      \_______________/        \
          | |___||__|      /       /                \          \
          |          \    /      /                    \          \
          |   (_______) /______/                        \_________ \
          |         / /         \                      /            \
           \         \^\\         \                  /               \     /
             \         ||           \______________/      _-_       //\__//
               \       ||------_-~~-_ ------------- \ --/~   ~\    || __/
                 ~-----||====/~     |==================|       |/~~~~~
                  (_(__/  ./     /                    \_\      \.
                         (_(___/                         \_____)_)
```

Il existe divers dépôts où vous pouvez trouver des conteneurs publiés.
Nous avons utilisé le service [Seqera Containers](https://seqera.io/containers/) pour générer cette image de conteneur Docker à partir du paquet Conda `cowpy` : `'community.wave.seqera.io/library/cowpy:1.1.5--3db457ae1977a273'`.

Exécutez la commande pull complète :

```bash
docker pull 'community.wave.seqera.io/library/cowpy:1.1.5--3db457ae1977a273'
```

??? success "Sortie de la commande"

    ```console
    1.1.5--3db457ae1977a273: Pulling from library/cowpy
    dafa2b0c44d2: Pull complete
    dec6b097362e: Pull complete
    f88da01cff0b: Pull complete
    4f4fb700ef54: Pull complete
    92dc97a3ef36: Pull complete
    403f74b0f85e: Pull complete
    10b8c00c10a5: Pull complete
    17dc7ea432cc: Pull complete
    bb36d6c3110d: Pull complete
    0ea1a16bbe82: Pull complete
    030a47592a0a: Pull complete
    c23bdb422167: Pull complete
    e1686ff32a11: Pull complete
    Digest: sha256:1ebc0043e8cafa61203bf42d29fd05bd14e7b4298e5e8cf986504c15f5aa4160
    Status: Downloaded newer image for community.wave.seqera.io/library/cowpy:1.1.5--3db457ae1977a273
    community.wave.seqera.io/library/cowpy:1.1.5--3db457ae1977a273
    ```

Si vous n'avez jamais téléchargé l'image auparavant, cela peut prendre une minute. Une fois terminé, vous disposez d'une copie locale de l'image de conteneur.

### 1.2. Utiliser le conteneur pour exécuter `cowpy` en tant que commande unique

Une façon très courante d'utiliser les conteneurs est de les exécuter directement, _c'est-à-dire_ de manière non interactive.
C'est idéal pour exécuter des commandes ponctuelles.

La syntaxe générale est la suivante :

```bash title="Syntaxe"
docker run --rm '<container>' [tool command]
```

La partie `docker run --rm '<container>'` est l'instruction au système de conteneurs de créer une instance de conteneur à partir d'une image de conteneur et d'y exécuter une commande.
Le drapeau `--rm` indique au système d'arrêter l'instance de conteneur une fois la commande terminée.

La syntaxe `[tool command]` dépend de l'outil que vous utilisez et de la façon dont le conteneur est configuré.
Commençons simplement avec `cowpy`.

Entièrement assemblée, la commande d'exécution du conteneur ressemble à ceci ; allez-y et exécutez-la.

```bash
docker run --rm 'community.wave.seqera.io/library/cowpy:1.1.5--3db457ae1977a273' cowpy
```

??? success "Sortie de la commande"

    ```console
    ______________________________________________________
    < Cowacter, eyes:default, tongue:False, thoughts:False >
    ------------------------------------------------------
        \   ^__^
          \  (oo)\_______
            (__)\       )\/\
              ||----w |
              ||     ||
    ```

Le système a créé le conteneur, exécuté la commande `cowpy` avec ses paramètres, envoyé la sortie vers la console et enfin, arrêté l'instance de conteneur.

### 1.3. Utiliser le conteneur pour exécuter `cowpy` de manière interactive

Vous pouvez également exécuter un conteneur de manière interactive, ce qui vous donne une invite shell à l'intérieur du conteneur et vous permet d'expérimenter avec la commande.

#### 1.3.1. Démarrer le conteneur

Pour exécuter de manière interactive, nous ajoutons simplement `-it` à la commande `docker run`.
Optionnellement, nous pouvons spécifier le shell que nous voulons utiliser à l'intérieur du conteneur en ajoutant _par exemple_ `/bin/bash` à la commande.

```bash
docker run --rm -it 'community.wave.seqera.io/library/cowpy:1.1.5--3db457ae1977a273' /bin/bash
```

Notez que votre invite change pour quelque chose comme `(base) root@b645838b3314:/tmp#`, ce qui indique que vous êtes maintenant à l'intérieur du conteneur.

Vous pouvez le vérifier en exécutant `ls /` pour lister le contenu du répertoire depuis la racine du système de fichiers :

```bash
ls /
```

??? abstract "Sortie de la commande"

    ```console
    bin  boot  dev  etc  home  lib  lib64  media  mnt  opt  proc  root  run  sbin  srv  sys  tmp  usr  var
    ```

Nous utilisons `ls` ici au lieu de `tree` car l'utilitaire `tree` n'est pas disponible dans ce conteneur.
Vous pouvez voir que le système de fichiers à l'intérieur du conteneur est différent du système de fichiers de votre système hôte.

Une limitation de ce que nous venons de faire est que le conteneur est complètement isolé du système hôte par défaut.
Cela signifie que le conteneur ne peut pas accéder aux fichiers du système hôte à moins que vous ne l'autorisiez explicitement à le faire.

Nous allons vous montrer comment faire cela dans un instant.

#### 1.3.2. Exécuter la ou les commandes de l'outil souhaité

Maintenant que vous êtes à l'intérieur du conteneur, vous pouvez exécuter la commande `cowpy` directement et lui donner des paramètres.
Par exemple, la documentation de l'outil indique que nous pouvons changer le personnage ('cowacter') avec `-c`.

```bash
cowpy "Hello Containers" -c tux
```

??? success "Sortie de la commande"

    ```console
    __________________
    < Hello Containers >
    ------------------
      \
        \
            .--.
          |o_o |
          |:_/ |
          //   \ \
        (|     | )
        /'\_   _/`\
        \___)=(___/
    ```

Maintenant la sortie montre le pingouin Linux, Tux, au lieu de la vache par défaut, parce que nous avons spécifié le paramètre `-c tux`.

Parce que vous êtes à l'intérieur du conteneur, vous pouvez exécuter la commande `cowpy` autant de fois que vous le souhaitez, en variant les paramètres d'entrée, sans avoir à vous soucier des commandes Docker.

!!! Tip "Astuce"

    Utilisez le drapeau '-c' pour choisir un personnage différent, notamment :
    `beavis`, `cheese`, `daemon`, `dragonandcow`, `ghostbusters`, `kitty`, `moose`, `milk`, `stegosaurus`, `turkey`, `turtle`, `tux`

C'est sympa. Ce qui serait encore plus sympa, c'est si nous pouvions fournir notre `greetings.csv` en entrée.
Mais comme nous n'avons pas accès au système de fichiers, nous ne pouvons pas.

Corrigeons cela.

#### 1.3.3. Quitter le conteneur

Pour quitter le conteneur, vous pouvez taper `exit` à l'invite ou utiliser le raccourci clavier ++ctrl+d++.

```bash
exit
```

Votre invite devrait maintenant être revenue à ce qu'elle était avant de démarrer le conteneur.

#### 1.3.4. Monter des données dans le conteneur

Comme indiqué précédemment, le conteneur est isolé du système hôte par défaut.

Pour permettre au conteneur d'accéder au système de fichiers hôte, vous pouvez **monter** un **volume** depuis le système hôte dans le conteneur en utilisant la syntaxe suivante :

```bash title="Syntaxe"
-v <outside_path>:<inside_path>
```

Dans notre cas, `<outside_path>` sera le répertoire de travail actuel, nous pouvons donc simplement utiliser un point (`.`), et `<inside_path>` est juste un alias que nous inventons ; appelons-le `/my_project` (le chemin interne doit être absolu).

Pour monter un volume, nous remplaçons les chemins et ajoutons l'argument de montage de volume à la commande docker run comme suit :

```bash
docker run --rm -it -v .:/my_project 'community.wave.seqera.io/library/cowpy:1.1.5--3db457ae1977a273' /bin/bash
```

Cela monte le répertoire de travail actuel en tant que volume qui sera accessible sous `/my_project` à l'intérieur du conteneur.

Vous pouvez vérifier que cela fonctionne en listant le contenu de `/my_project` :

```bash
ls /my_project
```

??? success "Sortie de la commande"

    ```console
    data               hello-config.nf      hello-modules.nf   hello-world.nf  nextflow.config  solutions         work
    hello-channels.nf  hello-containers.nf  hello-workflow.nf  modules         results          test-params.json
    ```

Vous pouvez maintenant voir le contenu du répertoire de travail depuis l'intérieur du conteneur, y compris le fichier `greetings.csv` sous `data/`.

Cela a effectivement établi un tunnel à travers la paroi du conteneur que vous pouvez utiliser pour accéder à cette partie de votre système de fichiers.

#### 1.3.5. Utiliser les données montées

Maintenant que nous avons monté le répertoire de travail dans le conteneur, nous pouvons utiliser la commande `cowpy` pour afficher le contenu du fichier `greetings.csv`.

Pour ce faire, nous utiliserons `cat /my_project/data/greetings.csv | ` pour rediriger le contenu du fichier CSV vers la commande `cowpy`.

```bash
cat /my_project/data/greetings.csv | cowpy -c turkey
```

??? success "Sortie de la commande"

    ```console title="data/greetings.csv"
     ____________________
    / Hello,English,123  \
    | Bonjour,French,456 |
    \ Holà,Spanish,789   /
    --------------------
      \                                  ,+*^^*+___+++_
      \                           ,*^^^^              )
        \                       _+*                     ^**+_
        \                    +^       _ _++*+_+++_,         )
                  _+^^*+_    (     ,+*^ ^          \+_        )
                {       )  (    ,(    ,_+--+--,      ^)      ^\
                { (\@)    } f   ,(  ,+-^ __*_*_  ^^\_   ^\       )
              {:;-/    (_+*-+^^^^^+*+*<_ _++_)_    )    )      /
              ( /  (    (        ,___    ^*+_+* )   <    <      \
              U _/     )    *--<  ) ^\-----++__)   )    )       )
                (      )  _(^)^^))  )  )\^^^^^))^*+/    /       /
              (      /  (_))_^)) )  )  ))^^^^^))^^^)__/     +^^
            (     ,/    (^))^))  )  ) ))^^^^^^^))^^)       _)
              *+__+*       (_))^)  ) ) ))^^^^^^))^^^^^)____*^
              \             \_)^)_)) ))^^^^^^^^^^))^^^^)
              (_             ^\__^^^^^^^^^^^^))^^^^^^^)
                ^\___            ^\__^^^^^^))^^^^^^^^)\\
                      ^^^^^\uuu/^^\uuu/^^^^\^\^\^\^\^\^\^\
                        ___) >____) >___   ^\_\_\_\_\_\_\)
                        ^^^//\\_^^//\\_^       ^(\_\_\_\)
                          ^^^ ^^ ^^^ ^
    ```

Cela produit l'art ASCII souhaité d'une dinde récitant nos exemples de salutations !
Sauf qu'ici la dinde répète les lignes complètes au lieu de juste les salutations.
Nous savons déjà que notre workflow Nextflow fera un meilleur travail !

N'hésitez pas à expérimenter avec cette commande.
Lorsque vous avez terminé, quittez le conteneur comme précédemment :

```bash
exit
```

Vous vous retrouverez dans votre shell normal.

### À retenir

Vous savez comment tirer un conteneur et l'exécuter soit en tant que commande unique soit de manière interactive. Vous savez également comment rendre vos données accessibles depuis l'intérieur de votre conteneur, ce qui vous permet d'essayer n'importe quel outil qui vous intéresse sur des données réelles sans avoir à installer de logiciel sur votre système.

### Et ensuite ?

Apprenez comment utiliser les conteneurs pour l'exécution des processus Nextflow.

---

## 2. Utiliser des conteneurs dans Nextflow

Nextflow dispose d'un support intégré pour exécuter des processus à l'intérieur de conteneurs afin de vous permettre d'exécuter des outils que vous n'avez pas installés dans votre environnement de calcul.
Cela signifie que vous pouvez utiliser n'importe quelle image de conteneur de votre choix pour exécuter vos processus, et Nextflow se chargera de tirer l'image, de monter les données et d'exécuter le processus à l'intérieur.

Pour démontrer cela, nous allons ajouter une étape `cowpy` au pipeline que nous avons développé, après l'étape `collectGreetings`.

<figure class="excalidraw">
--8<-- "docs/en/docs/hello_nextflow/img/hello-pipeline-cowpy.svg"
</figure>

### 2.1. Écrire un module `cowpy`

Tout d'abord, créons le module de processus `cowpy`.

#### 2.1.1. Créer un fichier squelette pour le nouveau module

Créez un fichier vide pour le module appelé `cowpy.nf`.

```bash
touch modules/cowpy.nf
```

Cela nous donne un endroit où mettre le code du processus.

#### 2.1.2. Copier le code du processus `cowpy` dans le fichier du module

Nous pouvons modeler notre processus `cowpy` sur les autres processus que nous avons écrits précédemment.

```groovy title="modules/cowpy.nf" linenums="1"
#!/usr/bin/env nextflow

// Générer de l'art ASCII avec cowpy
process cowpy {

    input:
    path input_file
    val character

    output:
    path "cowpy-${input_file}"

    script:
    """
    cat ${input_file} | cowpy -c "${character}" > cowpy-${input_file}
    """

}
```

Le processus attend un `input_file` contenant les salutations ainsi qu'une valeur `character`.

La sortie sera un nouveau fichier texte contenant l'art ASCII généré par l'outil `cowpy`.

### 2.2. Ajouter cowpy au workflow

Maintenant nous devons importer le module et appeler le processus.

#### 2.2.1. Importer le processus `cowpy` dans `hello-containers.nf`

Insérez la déclaration d'importation au-dessus du bloc workflow et remplissez-la de manière appropriée.

=== "Après"

    ```groovy title="hello-containers.nf" linenums="3" hl_lines="5"
    // Inclure les modules
    include { sayHello } from './modules/sayHello.nf'
    include { convertToUpper } from './modules/convertToUpper.nf'
    include { collectGreetings } from './modules/collectGreetings.nf'
    include { cowpy } from './modules/cowpy.nf'
    ```

=== "Avant"

    ```groovy title="hello-containers.nf" linenums="3"
    // Inclure les modules
    include { sayHello } from './modules/sayHello.nf'
    include { convertToUpper } from './modules/convertToUpper.nf'
    include { collectGreetings } from './modules/collectGreetings.nf'
    ```

Maintenant le module `cowpy` est disponible pour être utilisé dans le workflow.

#### 2.2.2. Ajouter un appel au processus `cowpy` dans le workflow

Connectons le processus `cowpy()` à la sortie du processus `collectGreetings()`, qui comme vous vous en souvenez peut-être produit deux sorties :

- `collectGreetings.out.outfile` contient le fichier de sortie <--_ce que nous voulons_
- `collectGreetings.out.report` contient le fichier de rapport avec le nombre de salutations par lot

Dans le bloc workflow, effectuez la modification de code suivante :

=== "Après"

    ```groovy title="hello-containers.nf" linenums="19" hl_lines="12-13"
        main:
        // créer un canal pour les entrées depuis un fichier CSV
        greeting_ch = channel.fromPath(params.input)
                            .splitCsv()
                            .map { line -> line[0] }
        // émettre une salutation
        sayHello(greeting_ch)
        // convertir la salutation en majuscules
        convertToUpper(sayHello.out)
        // collecter toutes les salutations dans un fichier
        collectGreetings(convertToUpper.out.collect(), params.batch)
        // générer de l'art ASCII des salutations avec cowpy
        cowpy(collectGreetings.out.outfile, params.character)
    ```

=== "Avant"

    ```groovy title="hello-containers.nf" linenums="19"
        main:
        // créer un canal pour les entrées depuis un fichier CSV
        greeting_ch = channel.fromPath(params.input)
                            .splitCsv()
                            .map { line -> line[0] }
        // émettre une salutation
        sayHello(greeting_ch)
        // convertir la salutation en majuscules
        convertToUpper(sayHello.out)
        // collecter toutes les salutations dans un fichier
        collectGreetings(convertToUpper.out.collect(), params.batch)
    ```

Notez que nous avons déclaré un nouveau paramètre CLI, `params.character`, afin de spécifier quel personnage nous voulons pour dire les salutations.

#### 2.2.3. Ajouter le paramètre `character` au bloc `params`

C'est techniquement optionnel mais c'est la pratique recommandée et c'est l'occasion de définir une valeur par défaut pour le personnage pendant que nous y sommes.

=== "Après"

    ```groovy title="hello-containers.nf" linenums="9" hl_lines="7"
    /*
    * Paramètres du pipeline
    */
    params {
        input: Path = 'data/greetings.csv'
        batch: String = 'batch'
        character: String = 'turkey'
    }
    ```

=== "Avant"

    ```groovy title="hello-containers.nf" linenums="9"
    /*
    * Paramètres du pipeline
    */
    params {
        input: Path = 'data/greetings.csv'
        batch: String = 'batch'
    }
    ```

Maintenant nous pouvons être paresseux·ses et éviter de taper le paramètre character dans nos lignes de commande.

#### 2.2.4. Mettre à jour les sorties du workflow

Nous devons mettre à jour les sorties du workflow pour publier la sortie du processus `cowpy`.

##### 2.2.4.1. Mettre à jour la section `publish:`

Dans le `bloc workflow`, effectuez la modification de code suivante :

=== "Après"

    ```groovy title="hello-containers.nf" linenums="34" hl_lines="6"
        publish:
        first_output = sayHello.out
        uppercased = convertToUpper.out
        collected = collectGreetings.out.outfile
        batch_report = collectGreetings.out.report
        cowpy_art = cowpy.out
    ```

=== "Avant"

    ```groovy title="hello-containers.nf" linenums="34"
        publish:
        first_output = sayHello.out
        uppercased = convertToUpper.out
        collected = collectGreetings.out.outfile
        batch_report = collectGreetings.out.report
    ```

Le processus `cowpy` ne produit qu'une seule sortie, nous pouvons donc y faire référence de la manière habituelle en ajoutant `.out`.

Mais pour l'instant, finissons de mettre à jour les sorties au niveau du workflow.

##### 2.2.4.2. Mettre à jour le bloc `output`

Nous devons ajouter la sortie finale `cowpy_art` au bloc `output`. Pendant que nous y sommes, modifions également les destinations de publication puisque maintenant notre pipeline est complet et nous savons quelles sorties nous intéressent vraiment.

Dans le bloc `output`, effectuez les modifications de code suivantes :

=== "Après"

    ```groovy title="hello-containers.nf" linenums="42" hl_lines="3 7 11 15 18-21"
    output {
        first_output {
            path 'hello_containers/intermediates'
            mode 'copy'
        }
        uppercased {
            path 'hello_containers/intermediates'
            mode 'copy'
        }
        collected {
            path 'hello_containers/intermediates'
            mode 'copy'
        }
        batch_report {
            path 'hello_containers'
            mode 'copy'
        }
        cowpy_art {
            path 'hello_containers'
            mode 'copy'
        }
    }
    ```

=== "Avant"

    ```groovy title="hello-containers.nf" linenums="42" hl_lines="3 7 11 15"
    output {
        first_output {
            path 'hello_containers'
            mode 'copy'
        }
        uppercased {
            path 'hello_containers'
            mode 'copy'
        }
        collected {
            path 'hello_containers'
            mode 'copy'
        }
        batch_report {
            path 'hello_containers'
            mode 'copy'
        }
    }
    ```

Maintenant les sorties publiées seront un peu plus organisées.

#### 2.2.5. Exécuter le workflow

Pour récapituler, voici ce que nous visons :

<figure class="excalidraw">
--8<-- "docs/en/docs/hello_nextflow/img/hello_pipeline_complete.svg"
</figure>

Pensez-vous que cela va fonctionner ?

Supprimons les sorties publiées précédentes pour avoir une table rase, et exécutons le workflow avec le drapeau `-resume`.

```bash
rm -r hello_containers/
nextflow run hello-containers.nf -resume
```

??? failure "Sortie de la commande (éditée pour plus de clarté)"

    ```console hl_lines="10 13 20-21 26-27"
     N E X T F L O W   ~  version 25.10.2

    Launching `hello-containers.nf` [lonely_woese] DSL2 - revision: abf1dccf7f

    executor >  local (1)
    [c9/f5c686] sayHello (3)       [100%] 3 of 3, cached: 3 ✔
    [ef/3135a8] convertToUpper (3) [100%] 3 of 3, cached: 3 ✔
    [7f/f435e3] collectGreetings   [100%] 1 of 1, cached: 1 ✔
    [9b/02e776] cowpy              [  0%] 0 of 1 ✘
    ERROR ~ Error executing process > 'cowpy'

    Caused by:
      Process `cowpy` terminated with an error exit status (127)


    Command executed:

      cat COLLECTED-batch-output.txt | cowpy -c "turkey" > cowpy-COLLECTED-batch-output.txt

    Command exit status:
      127

    Command output:
      (empty)

    Command error:
      .command.sh: line 2: cowpy: command not found

    Work dir:
      /workspaces/training/hello-nextflow/work/9b/02e7761db848f82db3c3e59ff3a9b6

    Tip: when you have fixed the problem you can continue the execution adding the option `-resume` to the run command line

    -- Check '.nextflow.log' file for details
    ERROR ~ Cannot access first() element from an empty List

    -- Check '.nextflow.log' file for details
    ```

Oh non, il y a une erreur !
Le code d'erreur donné par `error exit status (127)` signifie que l'exécutable que nous avons demandé n'a pas été trouvé.

C'est logique, puisque nous appelons l'outil `cowpy` mais nous n'avons pas encore spécifié de conteneur (oups).

### 2.3. Utiliser un conteneur pour exécuter le processus `cowpy`

Nous devons spécifier un conteneur et dire à Nextflow de l'utiliser pour le processus `cowpy()`.

#### 2.3.1. Spécifier un conteneur pour `cowpy`

Nous pouvons utiliser la même image que nous utilisions directement dans la première section de ce tutoriel.

Modifiez le module `cowpy.nf` pour ajouter la directive `container` à la définition du processus comme suit :

=== "Après"

    ```groovy title="modules/cowpy.nf" linenums="4" hl_lines="3"
    process cowpy {

        container 'community.wave.seqera.io/library/cowpy:1.1.5--3db457ae1977a273'

        input:
        path input_file
        val character

        output:
        path "cowpy-${input_file}"

        script:
        """
        cat ${input_file} | cowpy -c "${character}" > cowpy-${input_file}
        """
    }
    ```

=== "Avant"

    ```groovy title="modules/cowpy.nf" linenums="4"
    process cowpy {

        input:
        path input_file
        val character

        output:
        path "cowpy-${input_file}"

        script:
        """
        cat ${input_file} | cowpy -c "${character}" > cowpy-${input_file}
        """
    }
    ```

Cela indique à Nextflow que _si l'utilisation de Docker est activée_, il doit utiliser l'image de conteneur spécifiée ici pour exécuter le processus.

#### 2.3.2. Activer l'utilisation de Docker via le fichier `nextflow.config`

Notez que nous avons dit _'si l'utilisation de Docker est activée'_. Par défaut, elle ne l'est pas, nous devons donc dire à Nextflow qu'il est autorisé à utiliser Docker.
À cette fin, nous allons légèrement anticiper le sujet de la prochaine et dernière partie de cette formation (Partie 6), qui couvre la configuration.

L'un des principaux moyens que Nextflow offre pour configurer l'exécution du workflow est d'utiliser un fichier `nextflow.config`.
Lorsqu'un tel fichier est présent dans le répertoire actuel, Nextflow le chargera automatiquement et appliquera toute configuration qu'il contient.

Nous avons fourni un fichier `nextflow.config` avec une seule ligne de code qui désactive explicitement Docker : `docker.enabled = false`.

Maintenant, changeons cela en `true` pour activer Docker :

=== "Après"

    ```console title="nextflow.config" linenums="1" hl_lines="1"
    docker.enabled = true
    ```

=== "Avant"

    ```console title="nextflow.config" linenums="1" hl_lines="1"
    docker.enabled = false
    ```

!!! tip "Astuce"

    Il est possible d'activer l'exécution Docker depuis la ligne de commande, pour une exécution donnée, en utilisant le paramètre `-with-docker <container>`.
    Cependant, cela ne nous permet de spécifier qu'un seul conteneur pour l'ensemble du workflow, alors que l'approche que nous venons de vous montrer nous permet de spécifier un conteneur différent par processus.
    C'est meilleur pour la modularité, la maintenance du code et la reproductibilité.

#### 2.3.3. Exécuter le workflow avec Docker activé

Exécutez le workflow avec le drapeau `-resume` :

```bash
nextflow run hello-containers.nf -resume
```

??? success "Sortie de la commande"

    ```console
     N E X T F L O W   ~  version 25.10.2

    Launching `hello-containers.nf` [drunk_perlman] DSL2 - revision: abf1dccf7f

    executor >  local (1)
    [c9/f5c686] sayHello (3)       [100%] 3 of 3, cached: 3 ✔
    [ef/3135a8] convertToUpper (3) [100%] 3 of 3, cached: 3 ✔
    [7f/f435e3] collectGreetings   [100%] 1 of 1, cached: 1 ✔
    [98/656c6c] cowpy              [100%] 1 of 1 ✔
    ```

Cette fois, cela fonctionne effectivement !
Comme d'habitude, vous pouvez trouver les sorties du workflow dans le répertoire de résultats correspondant, bien que cette fois elles soient un peu plus soigneusement organisées, avec seulement le rapport et la sortie finale au niveau supérieur, et tous les fichiers intermédiaires rangés dans un sous-répertoire.

??? abstract "Contenu du répertoire"

    ```console
    results/hello_containers/
    ├── cowpy-COLLECTED-batch-output.txt
    ├── intermediates
    │   ├── Bonjour-output.txt
    │   ├── COLLECTED-batch-output.txt
    │   ├── Hello-output.txt
    │   ├── Holà-output.txt
    │   ├── UPPER-Bonjour-output.txt
    │   ├── UPPER-Hello-output.txt
    │   └── UPPER-Holà-output.txt
    └── batch-report.txt
    ```

La sortie finale d'art ASCII se trouve dans le répertoire `results/hello_containers/`, sous le nom `cowpy-COLLECTED-batch-output.txt`.

??? abstract "Contenu du fichier"

    ```console title="results/hello_containers/cowpy-COLLECTED-batch-output.txt"
    _________
    / HOLà    \
    | HELLO   |
    \ BONJOUR /
    ---------
      \                                  ,+*^^*+___+++_
      \                           ,*^^^^              )
        \                       _+*                     ^**+_
        \                    +^       _ _++*+_+++_,         )
                  _+^^*+_    (     ,+*^ ^          \+_        )
                {       )  (    ,(    ,_+--+--,      ^)      ^\
                { (\@)    } f   ,(  ,+-^ __*_*_  ^^\_   ^\       )
              {:;-/    (_+*-+^^^^^+*+*<_ _++_)_    )    )      /
              ( /  (    (        ,___    ^*+_+* )   <    <      \
              U _/     )    *--<  ) ^\-----++__)   )    )       )
                (      )  _(^)^^))  )  )\^^^^^))^*+/    /       /
              (      /  (_))_^)) )  )  ))^^^^^))^^^)__/     +^^
            (     ,/    (^))^))  )  ) ))^^^^^^^))^^)       _)
              *+__+*       (_))^)  ) ) ))^^^^^^))^^^^^)____*^
              \             \_)^)_)) ))^^^^^^^^^^))^^^^)
              (_             ^\__^^^^^^^^^^^^))^^^^^^^)
                ^\___            ^\__^^^^^^))^^^^^^^^)\\
                      ^^^^^\uuu/^^\uuu/^^^^\^\^\^\^\^\^\^\
                        ___) >____) >___   ^\_\_\_\_\_\_\)
                        ^^^//\\_^^//\\_^       ^(\_\_\_\)
                          ^^^ ^^ ^^^ ^
    ```

Et voilà, notre belle dinde disant les salutations comme souhaité.

#### 2.3.4. Inspecter comment Nextflow a lancé la tâche conteneurisée

En guise de coda final à cette section, jetons un coup d'œil au sous-répertoire de travail pour l'un des appels de processus `cowpy` afin d'obtenir un peu plus d'informations sur la façon dont Nextflow fonctionne avec les conteneurs en coulisses.

Vérifiez la sortie de votre commande `nextflow run` pour trouver le chemin vers le sous-répertoire de travail pour le processus `cowpy`.
En regardant ce que nous avons obtenu pour l'exécution montrée ci-dessus, la ligne de journal de console pour le processus `cowpy` commence par `[98/656c6c]`.
Cela correspond au chemin de répertoire tronqué suivant : `work/98/656c6c`.

Dans ce répertoire, vous trouverez le fichier `.command.run` qui contient toutes les commandes que Nextflow a exécutées en votre nom au cours de l'exécution du pipeline.

??? abstract "Contenu du fichier"

    ```console title="work/98/656c6c90cce1667c094d880f4b6dcc/.command.run"
    #!/bin/bash
    ### ---
    ### name: 'cowpy'
    ### container: 'community.wave.seqera.io/library/cowpy:1.1.5--3db457ae1977a273'
    ### outputs:
    ### - 'cowpy-COLLECTED-batch-output.txt'
    ### ...
    set -e
    set -u
    NXF_DEBUG=${NXF_DEBUG:=0}; [[ $NXF_DEBUG > 1 ]] && set -x
    NXF_ENTRY=${1:-nxf_main}


    nxf_sleep() {
      sleep $1 2>/dev/null || sleep 1;
    }

    nxf_date() {
        local ts=$(date +%s%3N);
        if [[ ${#ts} == 10 ]]; then echo ${ts}000
        elif [[ $ts == *%3N ]]; then echo ${ts/\%3N/000}
        elif [[ $ts == *3N ]]; then echo ${ts/3N/000}
        elif [[ ${#ts} == 13 ]]; then echo $ts
        else echo "Unexpected timestamp value: $ts"; exit 1
        fi
    }

    nxf_env() {
        echo '============= task environment ============='
        env | sort | sed "s/\(.*\)AWS\(.*\)=\(.\{6\}\).*/\1AWS\2=\3xxxxxxxxxxxxx/"
        echo '============= task output =================='
    }

    nxf_kill() {
        declare -a children
        while read P PP;do
            children[$PP]+=" $P"
        done < <(ps -e -o pid= -o ppid=)

        kill_all() {
            [[ $1 != $$ ]] && kill $1 2>/dev/null || true
            for i in ${children[$1]:=}; do kill_all $i; done
        }

        kill_all $1
    }

    nxf_mktemp() {
        local base=${1:-/tmp}
        mkdir -p "$base"
        if [[ $(uname) = Darwin ]]; then mktemp -d $base/nxf.XXXXXXXXXX
        else TMPDIR="$base" mktemp -d -t nxf.XXXXXXXXXX
        fi
    }

    nxf_fs_copy() {
      local source=$1
      local target=$2
      local basedir=$(dirname $1)
      mkdir -p $target/$basedir
      cp -fRL $source $target/$basedir
    }

    nxf_fs_move() {
      local source=$1
      local target=$2
      local basedir=$(dirname $1)
      mkdir -p $target/$basedir
      mv -f $source $target/$basedir
    }

    nxf_fs_rsync() {
      rsync -rRl $1 $2
    }

    nxf_fs_rclone() {
      rclone copyto $1 $2/$1
    }

    nxf_fs_fcp() {
      fcp $1 $2/$1
    }

    on_exit() {
        local last_err=$?
        local exit_status=${nxf_main_ret:=0}
        [[ ${exit_status} -eq 0 && ${nxf_unstage_ret:=0} -ne 0 ]] && exit_status=${nxf_unstage_ret:=0}
        [[ ${exit_status} -eq 0 && ${last_err} -ne 0 ]] && exit_status=${last_err}
        printf -- $exit_status > /workspaces/training/hello-nextflow/work/98/656c6c90cce1667c094d880f4b6dcc/.exitcode
        set +u
        docker rm $NXF_BOXID &>/dev/null || true
        exit $exit_status
    }

    on_term() {
        set +e
        docker stop $NXF_BOXID
    }

    nxf_launch() {
        docker run -i --cpu-shares 1024 -e "NXF_TASK_WORKDIR" -v /workspaces/training/hello-nextflow/work:/workspaces/training/hello-nextflow/work -w "$NXF_TASK_WORKDIR" --name $NXF_BOXID community.wave.seqera.io/library/cowpy:1.1.5--3db457ae1977a273 /bin/bash -ue /workspaces/training/hello-nextflow/work/98/656c6c90cce1667c094d880f4b6dcc/.command.sh
    }

    nxf_stage() {
        true
        # stage input files
        rm -f COLLECTED-batch-output.txt
        ln -s /workspaces/training/hello-nextflow/work/7f/f435e3f2cf95979b5f3d7647ae6696/COLLECTED-batch-output.txt COLLECTED-batch-output.txt
    }

    nxf_unstage_outputs() {
        true
    }

    nxf_unstage_controls() {
        true
    }

    nxf_unstage() {
        if [[ ${nxf_main_ret:=0} == 0 ]]; then
            (set -e -o pipefail; (nxf_unstage_outputs | tee -a .command.out) 3>&1 1>&2 2>&3 | tee -a .command.err)
            nxf_unstage_ret=$?
        fi
        nxf_unstage_controls
    }

    nxf_main() {
        trap on_exit EXIT
        trap on_term TERM INT USR2
        trap '' USR1

        [[ "${NXF_CHDIR:-}" ]] && cd "$NXF_CHDIR"
        export NXF_BOXID="nxf-$(dd bs=18 count=1 if=/dev/urandom 2>/dev/null | base64 | tr +/ 0A | tr -d '\r\n')"
        NXF_SCRATCH=''
        [[ $NXF_DEBUG > 0 ]] && nxf_env
        touch /workspaces/training/hello-nextflow/work/98/656c6c90cce1667c094d880f4b6dcc/.command.begin
        set +u
        set -u
        [[ $NXF_SCRATCH ]] && cd $NXF_SCRATCH
        export NXF_TASK_WORKDIR="$PWD"
        nxf_stage

        set +e
        (set -o pipefail; (nxf_launch | tee .command.out) 3>&1 1>&2 2>&3 | tee .command.err) &
        pid=$!
        wait $pid || nxf_main_ret=$?
        nxf_unstage
    }

    $NXF_ENTRY

    ```

Si vous recherchez `nxf_launch` dans ce fichier, vous devriez voir quelque chose comme ceci :

```console
nxf_launch() {
    docker run -i --cpu-shares 1024 -e "NXF_TASK_WORKDIR" -v /workspaces/training/hello-nextflow/work:/workspaces/training/hello-nextflow/work -w "$NXF_TASK_WORKDIR" --name $NXF_BOXID community.wave.seqera.io/library/cowpy:1.1.5--3db457ae1977a273 /bin/bash -ue /workspaces/training/hello-nextflow/work/98/656c6c90cce1667c094d880f4b6dcc/.command.sh
}
```

Comme vous pouvez le voir, Nextflow utilise la commande `docker run` pour lancer l'appel de processus.
Il monte également le sous-répertoire de travail correspondant dans le conteneur, définit le répertoire de travail à l'intérieur du conteneur en conséquence, et exécute notre script bash modélisé dans le fichier `.command.sh`.

Tout le travail difficile que nous devions faire manuellement dans la première section ? Nextflow le fait pour nous en coulisses !

```txt
 _______________________
< Hourra pour les robots...! >
 -----------------------
                                   ,-----.
                                   |     |
                                ,--|     |-.
                         __,----|  |     | |
                       ,;::     |  `_____' |
                       `._______|    i^i   |
                                `----| |---'| .
                           ,-------._| |== ||//
                           |       |_|P`.  /'/
                           `-------' 'Y Y/'/'
                                     .==\ /_\
   ^__^                             /   /'|  `i
   (oo)\_______                   /'   /  |   |
   (__)\       )\/\             /'    /   |   `i
       ||----w |           ___,;`----'.___L_,-'`\__
       ||     ||          i_____;----\.____i""\____\
```

### À retenir

Vous savez comment utiliser les conteneurs dans Nextflow pour exécuter des processus.

### Et ensuite ?

Faites une pause !

Lorsque vous êtes prêt·e, passez à [**Partie 6 : Hello Config**](./06_hello_config.md) pour apprendre comment configurer l'exécution de votre pipeline pour l'adapter à votre infrastructure ainsi que gérer la configuration des entrées et des paramètres.

C'est la toute dernière partie, et ensuite vous aurez terminé cette formation !

---

## Quiz

<quiz>
Qu'est-ce qu'un conteneur ?
- [ ] Un type de machine virtuelle
- [ ] Un format de compression de fichiers
- [x] Une unité exécutable légère et autonome qui inclut tout le nécessaire pour exécuter une application
- [ ] Un protocole réseau
</quiz>

<quiz>
Quelle est la différence entre une image de conteneur et une instance de conteneur ?
- [ ] Ce sont la même chose
- [x] Une image est un modèle ; une instance est un conteneur en cours d'exécution créé à partir de cette image
- [ ] Une instance est un modèle ; une image est un conteneur en cours d'exécution
- [ ] Les images sont pour Docker ; les instances sont pour Singularity
</quiz>

<quiz>
Que fait le drapeau `-v` dans une commande `docker run` ?
- [ ] Active la sortie détaillée
- [ ] Valide le conteneur
- [x] Monte un volume depuis le système hôte dans le conteneur
- [ ] Spécifie la version du conteneur

En savoir plus : [1.3.4. Monter des données dans le conteneur](#134-monter-des-données-dans-le-conteneur)
</quiz>

<quiz>
Pourquoi devez-vous monter des volumes lors de l'utilisation de conteneurs ?
- [ ] Pour améliorer les performances du conteneur
- [ ] Pour économiser de l'espace disque
- [x] Parce que les conteneurs sont isolés du système de fichiers hôte par défaut
- [ ] Pour activer le réseau

En savoir plus : [1.3.4. Monter des données dans le conteneur](#134-monter-des-données-dans-le-conteneur)
</quiz>

<quiz>
Comment spécifier un conteneur pour un processus Nextflow ?
- [ ] `docker 'container-uri'`
- [ ] `image 'container-uri'`
- [x] `container 'container-uri'`
- [ ] `use 'container-uri'`

En savoir plus : [2.3.1. Spécifier un conteneur pour cowpy](#231-spécifier-un-conteneur-pour-cowpy)
</quiz>

<quiz>
Quel paramètre `nextflow.config` active Docker pour votre workflow ?
- [ ] `#!groovy process.docker = true`
- [x] `#!groovy docker.enabled = true`
- [ ] `#!groovy container.engine = 'docker'`
- [ ] `#!groovy docker.activate = true`

En savoir plus : [2.3.2. Activer l'utilisation de Docker via le fichier `nextflow.config`](#232-activer-lutilisation-de-docker-via-le-fichier-nextflowconfig)
</quiz>

<quiz>
Que gère automatiquement Nextflow lors de l'exécution d'un processus dans un conteneur ? (Sélectionnez toutes les réponses applicables)
- [x] Tirer l'image de conteneur si nécessaire
- [x] Monter le répertoire de travail
- [x] Exécuter le script du processus à l'intérieur du conteneur
- [x] Nettoyer l'instance de conteneur après l'exécution

En savoir plus : [2.3.4. Inspecter comment Nextflow a lancé la tâche conteneurisée](#234-inspecter-comment-nextflow-a-lancé-la-tâche-conteneurisée)
</quiz>
