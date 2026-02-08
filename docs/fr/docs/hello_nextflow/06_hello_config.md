# Partie 6 : Hello Config

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Traduction assistée par IA - [en savoir plus et suggérer des améliorations](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

<div class="video-wrapper">
  <iframe width="560" height="315" src="https://www.youtube.com/embed/FcZTiE25TeA?si=tnXTi6mRkITY0zW_&amp;list=PLPZ8WHdZGxmWKozQuzr27jyMGqp9kElVK&amp;cc_load_policy=1&amp;cc_lang_pref=fr" title="YouTube video player" frameborder="0" allow="accelerometer; autoplay; clipboard-write; encrypted-media; gyroscope; picture-in-picture; web-share" referrerpolicy="strict-origin-when-cross-origin" allowfullscreen></iframe>
</div>

/// caption
:fontawesome-brands-youtube:{ .youtube } Voir [la playlist complète](https://youtube.com/playlist?list=PLPZ8WHdZGxmWKozQuzr27jyMGqp9kElVK&si=eF7cLR62goy-lc6n) sur la chaîne YouTube de Nextflow.

:green_book: La transcription de la vidéo est disponible [ici](./transcripts/06_hello_config.md).
///

Cette section explorera comment configurer et gérer la configuration de votre pipeline Nextflow afin que vous puissiez personnaliser son comportement, l'adapter à différents environnements et optimiser l'utilisation des ressources _sans modifier une seule ligne du code du workflow lui-même_.

Il existe plusieurs façons de le faire, qui peuvent être utilisées en combinaison et sont interprétées selon l'[ordre de priorité](https://nextflow.io/docs/latest/config.html) décrit dans la documentation de configuration.

Dans cette partie du cours, nous allons vous montrer le mécanisme de fichier de configuration le plus simple et le plus courant, le fichier [`nextflow.config`](https://nextflow.io/docs/latest/config.html), que vous avez déjà rencontré dans la Partie 5 : Hello Containers.

Nous passerons en revue les composants essentiels de la configuration Nextflow tels que les directives de processus, les exécuteurs, les profils et les fichiers de paramètres.
En apprenant à utiliser efficacement ces options de configuration, vous pouvez améliorer la flexibilité, l'évolutivité et les performances de vos pipelines.

??? info "Comment commencer à partir de cette section"

    Cette section du cours suppose que vous avez complété les Parties 1-5 du cours [Hello Nextflow](./index.md) et que vous avez un pipeline fonctionnel complet.

    Si vous commencez le cours à partir de ce point, vous devrez copier le répertoire `modules` et le fichier `nextflow.config` depuis les solutions :

    ```bash
    cp -r solutions/5-hello-containers/modules .
    cp solutions/5-hello-containers/nextflow.config .
    ```

    Le fichier `nextflow.config` contient la ligne `docker.enabled = true` qui active l'utilisation des conteneurs Docker.

    Si vous n'êtes pas familier·ère avec le pipeline Hello ou si vous avez besoin d'un rappel, consultez [cette page d'information](../info/hello_pipeline.md).

---

## 0. Échauffement : Exécuter `hello-config.nf`

Nous allons utiliser le script de workflow `hello-config.nf` comme point de départ.
Il est équivalent au script produit en travaillant à travers la Partie 5 de ce cours de formation, sauf que nous avons changé les destinations de sortie :

```groovy title="hello-config.nf" linenums="37" hl_lines="3 7 11 15"
output {
    first_output {
        path 'hello_config/intermediates'
        mode 'copy'
    }
    uppercased {
        path 'hello_config/intermediates'
        mode 'copy'
    }
    collected {
        path 'hello_config/intermediates'
        mode 'copy'
    }
    batch_report {
        path 'hello_config'
        mode 'copy'
    }
    cowpy_art {
        path 'hello_config'
        mode 'copy'
    }
}
```

Juste pour s'assurer que tout fonctionne, exécutez le script une fois avant d'effectuer des modifications :

```bash
nextflow run hello-config.nf
```

??? success "Sortie de la commande"

    ```console
     N E X T F L O W   ~  version 25.10.2

    Launching `hello-config.nf` [nice_escher] DSL2 - revision: d5dfdc9872

    executor > local (7)
    [6a/bc46a6] sayHello (2) [100%] 3 of 3 ✔
    [33/67bc48] convertToUpper (3) [100%] 3 of 3 ✔
    [b5/de03ba] collectGreetings [100%] 1 of 1 ✔
    [98/c6b57b] cowpy              | 1 of 1 ✔
    ```

Comme précédemment, vous trouverez les fichiers de sortie dans le répertoire spécifié dans le bloc `output` (`results/hello_config/`).

??? abstract "Contenu du répertoire"

    ```console
    results/hello_config/
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

La sortie finale d'art ASCII est dans le répertoire `results/hello_config/`, sous le nom `cowpy-COLLECTED-batch-output.txt`.

??? abstract "Contenu du fichier"

    ```console title="results/hello_config/cowpy-COLLECTED-batch-output.txt"
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

Si cela a fonctionné pour vous, vous êtes prêt·e à apprendre comment configurer vos pipelines.

---

## 1. Gérer les paramètres d'entrée du workflow

Nous allons commencer par un aspect de la configuration qui est simplement une extension de ce avec quoi nous avons travaillé jusqu'à présent : la gestion des paramètres d'entrée.

Actuellement, notre workflow est configuré pour accepter plusieurs valeurs de paramètres via la ligne de commande, avec des valeurs par défaut définies dans un bloc `params` dans le script de workflow lui-même.
Cependant, vous pourriez vouloir remplacer ces valeurs par défaut sans avoir à spécifier les paramètres sur la ligne de commande ou à modifier le fichier de script original.

Il existe plusieurs façons de le faire ; nous allons vous montrer trois façons de base qui sont très couramment utilisées.

### 1.1. Déplacer les valeurs par défaut vers `nextflow.config`

C'est l'approche la plus simple, bien qu'elle soit peut-être la moins flexible puisque le fichier `nextflow.config` principal n'est pas quelque chose que vous voulez éditer pour chaque exécution.
Mais elle a l'avantage de séparer les préoccupations de _déclarer_ les paramètres dans le workflow (ce qui appartient définitivement là) versus fournir des _valeurs par défaut_, qui sont plus à leur place dans un fichier de configuration.

Faisons cela en deux étapes.

#### 1.1.1. Créer un bloc `params` dans le fichier de configuration

Effectuez les modifications de code suivantes dans le fichier `nextflow.config` :

=== "Après"

    ```groovy title="nextflow.config" linenums="1" hl_lines="3-10"
    docker.enabled = true

    /*
    * Paramètres du pipeline
    */
    params {
        input = 'data/greetings.csv'
        batch = 'batch'
        character = 'turkey'
    }
    ```

=== "Avant"

    ```groovy title="nextflow.config" linenums="1"
    docker.enabled = true
    ```

Notez que nous n'avons pas simplement copié le bloc `params` du workflow vers le fichier de configuration.
La syntaxe est un peu différente.
Dans le fichier de workflow, ce sont des déclarations typées.
Dans la configuration, ce sont des affectations de valeurs.

Techniquement, cela suffit pour remplacer les valeurs par défaut encore spécifiées dans le fichier de workflow.
Vous pourriez modifier le personnage, par exemple, et exécuter le workflow pour vous assurer que la valeur définie dans le fichier de configuration remplace celle définie dans le fichier de workflow.

Mais dans l'esprit de déplacer complètement la configuration vers le fichier de configuration, supprimons entièrement ces valeurs du fichier de workflow.

#### 1.1.2. Supprimer les valeurs du bloc `params` dans le fichier de workflow

Effectuez les modifications de code suivantes dans le fichier de workflow `hello-config.nf` :

=== "Après"

    ```groovy title="hello-config.nf" linenums="9" hl_lines="5-7"
    /*
    * Paramètres du pipeline
    */
    params {
        input: Path
        batch: String
        character: String
    }
    ```

=== "Avant"

    ```groovy title="hello-config.nf" linenums="9" hl_lines="5-7"
    /*
    * Paramètres du pipeline
    */
    params {
        input: Path = 'data/greetings.csv'
        batch: String = 'batch'
        character: String = 'turkey'
    }
    ```

Maintenant, le fichier de workflow lui-même ne définit aucune valeur par défaut pour ces paramètres.

#### 1.1.3. Exécuter le pipeline

Testons que cela fonctionne correctement.

```bash
nextflow run hello-config.nf
```

??? success "Sortie de la commande"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `hello-config.nf` [disturbed_einstein] DSL2 - revision: ede9037d02

    executor >  local (8)
    [f0/35723c] sayHello (2)       | 3 of 3 ✔
    [40/3efd1a] convertToUpper (3) | 3 of 3 ✔
    [17/e97d32] collectGreetings   | 1 of 1 ✔
    [98/c6b57b] cowpy              | 1 of 1 ✔
    ```

Cela produit toujours la même sortie qu'avant.

La sortie finale d'art ASCII est dans le répertoire `results/hello_config/`, sous le nom `cowpy-COLLECTED-batch-output.txt`, comme avant.

??? abstract "Contenu du fichier"

    ```console title="results/hello_config/cowpy-COLLECTED-batch-output.txt"
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

Fonctionnellement, ce déplacement n'a rien changé, mais conceptuellement c'est un peu plus propre d'avoir les valeurs par défaut définies dans le fichier de configuration.

### 1.2. Utiliser un fichier de configuration spécifique à l'exécution

C'est très bien, mais parfois vous pourriez vouloir exécuter des expériences temporaires avec différentes valeurs par défaut sans toucher au fichier de configuration principal.
Vous pouvez le faire en créant un nouveau fichier `nextflow.config` dans un sous-répertoire que vous utiliserez comme répertoire de travail pour vos expériences.

#### 1.2.1. Créer le répertoire de travail avec une configuration vide

Commençons par créer un nouveau répertoire et nous y déplacer :

```bash
mkdir -p tux-run
cd tux-run
```

Ensuite, créez un fichier de configuration vide dans ce répertoire :

```bash
touch nextflow.config
```

Cela produit un fichier vide.

#### 1.2.2. Configurer la configuration expérimentale

Maintenant, ouvrez le nouveau fichier et ajoutez les paramètres que vous voulez personnaliser :

```groovy title="tux-run/nextflow.config" linenums="1"
params {
    input = '../data/greetings.csv'
    batch = 'experiment'
    character = 'tux'
}
```

Notez que le chemin vers le fichier d'entrée doit refléter la structure du répertoire.

#### 1.2.3. Exécuter le pipeline

Nous pouvons maintenant exécuter notre pipeline depuis notre nouveau répertoire de travail.
Assurez-vous d'adapter le chemin en conséquence !

```bash
nextflow run ../hello-config.nf
```

??? success "Sortie de la commande"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `../hello-config.nf` [trusting_escher] DSL2 - revision: 356df0818d

    executor >  local (8)
    [59/b66913] sayHello (2)       [100%] 3 of 3 ✔
    [ad/f06364] convertToUpper (3) [100%] 3 of 3 ✔
    [10/714895] collectGreetings   [100%] 1 of 1 ✔
    [88/3ece98] cowpy              [100%] 1 of 1 ✔
    ```

Cela créera un nouvel ensemble de répertoires sous `tux-run/` incluant `tux-run/work/` et `tux-run/results/`.

Dans cette exécution, Nextflow combine le `nextflow.config` dans notre répertoire actuel avec le `nextflow.config` dans le répertoire racine du pipeline, et remplace ainsi le personnage par défaut (turkey) par le personnage tux.

Le fichier de sortie final devrait contenir le personnage tux disant les messages de bienvenue.

??? abstract "Contenu du fichier"

    ```console title="tux-run/results/hello_config/cowpy-COLLECTED-experiment-output.txt"
    _________
    / HELLO   \
    | BONJOUR |
    \ HOLà    /
    ---------
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

Voilà ; maintenant vous avez un espace pour expérimenter sans modifier votre configuration « normale ».

!!! warning "Avertissement"

    Assurez-vous de revenir au répertoire précédent avant de passer à la section suivante !

    ```bash
    cd ..
    ```

Maintenant, regardons une autre façon utile de définir les valeurs des paramètres.

### 1.3. Utiliser un fichier de paramètres

L'approche du sous-répertoire fonctionne très bien pour l'expérimentation, mais elle implique un peu de configuration et nécessite que vous adaptiez les chemins en conséquence.
Il existe une approche plus simple lorsque vous voulez exécuter votre pipeline avec un ensemble spécifique de valeurs, ou permettre à quelqu'un d'autre de le faire avec un minimum d'effort.

Nextflow nous permet de spécifier des paramètres via un [fichier de paramètres](https://nextflow.io/docs/latest/config.html#params-file) au format YAML ou JSON, ce qui rend très pratique la gestion et la distribution d'ensembles alternatifs de valeurs par défaut, par exemple, ainsi que de valeurs de paramètres spécifiques à l'exécution.

#### 1.3.1. Examiner le fichier de paramètres exemple

Pour démontrer ceci, nous fournissons un fichier de paramètres exemple dans le répertoire actuel, appelé `test-params.yaml` :

```yaml title="test-params.yaml" linenums="1"
input: "data/greetings.csv"
batch: "yaml"
character: "stegosaurus"
```

Ce fichier de paramètres contient une paire clé-valeur pour chacune des entrées que nous voulons spécifier.
Notez l'utilisation de deux-points (`:`) au lieu de signes égal (`=`) si vous comparez la syntaxe avec le fichier de configuration.
Le fichier de configuration est écrit en Groovy, tandis que le fichier de paramètres est écrit en YAML.

!!! info "Info"

    Nous fournissons également une version JSON du fichier de paramètres comme exemple mais nous n'allons pas l'exécuter ici.
    N'hésitez pas à essayer celui-là par vous-même.

#### 1.3.2. Exécuter le pipeline

Pour exécuter le workflow avec ce fichier de paramètres, ajoutez simplement `-params-file <nom_de_fichier>` à la commande de base.

```bash
nextflow run hello-config.nf -params-file test-params.yaml
```

??? success "Sortie de la commande"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `hello-config.nf` [disturbed_sammet] DSL2 - revision: ede9037d02

    executor >  local (8)
    [f0/35723c] sayHello (2)       | 3 of 3 ✔
    [40/3efd1a] convertToUpper (3) | 3 of 3 ✔
    [17/e97d32] collectGreetings   | 1 of 1 ✔
    [98/c6b57b] cowpy              | 1 of 1 ✔
    ```

Le fichier de sortie final devrait contenir le personnage stegosaurus disant les messages de bienvenue.

??? abstract "Contenu du fichier"

    ```console title="results/hello_config/cowpy-COLLECTED-yaml-output.txt"
    _________
    / HELLO   \
    | HOLà    |
    \ BONJOUR /
    ---------
    \                             .       .
    \                           / `.   .' "
      \                  .---.  <    > <    >  .---.
      \                 |    \  \ - ~ ~ - /  /    |
            _____          ..-~             ~-..-~
            |     |   \~~~\.'                    `./~~~/
          ---------   \__/                        \__/
          .'  O    \     /               /       \  "
        (_____,    `._.'               |         }  \/~~~/
          `----.          /       }     |        /    \__/
                `-.      |       /      |       /      `. ,~~|
                    ~-.__|      /_ - ~ ^|      /- _      `..-'
                        |     /        |     /     ~-.     `-. _  _  _
                        |_____|        |_____|         ~ - . _ _ _ _ _>
    ```

Utiliser un fichier de paramètres peut sembler excessif lorsque vous n'avez que quelques paramètres à spécifier, mais certains pipelines attendent des dizaines de paramètres.
Dans ces cas, utiliser un fichier de paramètres nous permettra de fournir des valeurs de paramètres à l'exécution sans avoir à taper des lignes de commande massives et sans modifier le script de workflow.

Cela facilite également la distribution d'ensembles de paramètres aux collaborateurs·trices, ou comme information de support pour une publication, par exemple.
Cela rend votre travail plus reproductible par d'autres.

### À retenir

Vous savez comment tirer parti des options de configuration clés pour gérer les entrées du workflow.

### Et ensuite ?

Apprenez à gérer où et comment les sorties de votre workflow sont publiées.

---

## 2. Gérer les sorties du workflow

Jusqu'à présent, nous avons codé en dur tous les chemins pour les déclarations de sortie au niveau du workflow, et comme nous l'avons noté lorsque nous avons commencé à ajouter plusieurs sorties, il peut y avoir un peu de répétition impliquée.

Regardons quelques façons courantes de configurer cela pour être plus flexible.

### 2.1. Personnaliser le répertoire de sortie avec `-output-dir`

Lorsque nous contrôlons comment nos sorties « publiées » sont organisées, nous avons deux priorités distinctes :

- Le répertoire de sortie de niveau supérieur
- Comment les fichiers sont organisés dans ce répertoire

Nous avons utilisé le répertoire de niveau supérieur par défaut jusqu'à présent : `results`.
Commençons par personnaliser cela, en utilisant l'option CLI `-output-dir`.

#### 2.1.1. Exécuter le pipeline avec `-output-dir`

L'option `-output-dir` (raccourci : `-o`) remplace le répertoire de sortie par défaut (`results/`) pour toutes les sorties du workflow.
C'est la méthode recommandée pour contrôler le chemin racine où les sorties sont publiées.

```bash
nextflow run hello-config.nf -output-dir custom-outdir-cli/
```

??? success "Sortie de la commande"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `hello-config.nf` [prickly_kay] DSL2 - revision: 32ecc4fba2

    executor >  local (8)
    [9f/332636] sayHello (1)       [100%] 3 of 3 ✔
    [03/a55991] convertToUpper (3) [100%] 3 of 3 ✔
    [e5/ab7893] collectGreetings   [100%] 1 of 1 ✔
    [a8/97338e] cowpy              [100%] 1 of 1 ✔
    ```

Cela publie les sorties dans `custom-outdir-cli/` au lieu de `results/` :

??? abstract "Contenu du répertoire"

    ```console
    custom-outdir-cli/
    └── hello_config
        ├── batch-report.txt
        ├── cowpy-COLLECTED-batch-output.txt
        └── intermediates
            ├── Bonjour-output.txt
            ├── COLLECTED-batch-output.txt
            ├── Hello-output.txt
            ├── Holà-output.txt
            ├── UPPER-Bonjour-output.txt
            ├── UPPER-Hello-output.txt
            └── UPPER-Holà-output.txt
    ```

Notez que nous avons toujours le sous-répertoire `hello_config` provenant des déclarations `path` dans le bloc output.
Nettoyons cela.

#### 2.1.2. Supprimer les chemins codés en dur du bloc output

Le préfixe `hello_config/` était codé en dur dans les chapitres précédents, mais puisque nous apprenons maintenant à configurer les chemins de sortie de manière flexible, nous pouvons supprimer ce codage en dur.
Pour les sorties qui n'ont pas besoin d'un sous-répertoire, nous pouvons définir la directive `path` à une chaîne vide, ou la supprimer entièrement.

Effectuez les modifications de code suivantes dans le fichier de workflow :

=== "Après"

    ```groovy title="hello-config.nf" linenums="42" hl_lines="3 7 11 15 19"
    output {
        first_output {
            path 'intermediates'
            mode 'copy'
        }
        uppercased {
            path 'intermediates'
            mode 'copy'
        }
        collected {
            path 'intermediates'
            mode 'copy'
        }
        batch_report {
            path ''
            mode 'copy'
        }
        cowpy_art {
            path ''
            mode 'copy'
        }
    }
    ```

=== "Avant"

    ```groovy title="hello-config.nf" linenums="42" hl_lines="3 7 11 15 19"
    output {
        first_output {
            path 'hello_config/intermediates'
            mode 'copy'
        }
        uppercased {
            path 'hello_config/intermediates'
            mode 'copy'
        }
        collected {
            path 'hello_config/intermediates'
            mode 'copy'
        }
        batch_report {
            path 'hello_config'
            mode 'copy'
        }
        cowpy_art {
            path 'hello_config'
            mode 'copy'
        }
    }
    ```

Exécutez à nouveau le pipeline :

```bash
nextflow run hello-config.nf -output-dir custom-outdir-cli-2/
```

Maintenant les sorties sont publiées directement sous `custom-outdir-cli-2/`, sans le sous-répertoire `hello_config` :

??? abstract "Contenu du répertoire"

    ```console
    custom-outdir-cli-2/
    ├── batch-report.txt
    ├── cowpy-COLLECTED-batch-output.txt
    └── intermediates
        ├── Bonjour-output.txt
        ├── COLLECTED-batch-output.txt
        ├── Hello-output.txt
        ├── Holà-output.txt
        ├── UPPER-Bonjour-output.txt
        ├── UPPER-Hello-output.txt
        └── UPPER-Holà-output.txt
    ```

!!! tip "Astuce"

    L'option `-output-dir` est utilisée pour contrôler _où_ vont les sorties, tandis que la directive `path` dans le bloc output contrôle la _structure des sous-répertoires_.

### 2.2. Chemins de sortie dynamiques

En plus de changer le répertoire de sortie via le CLI, nous pouvons également définir une valeur par défaut personnalisée dans le fichier de configuration en utilisant `outputDir`.
Cela nous permet de définir le chemin du répertoire de manière dynamique - pas seulement en utilisant des chaînes statiques.

#### 2.2.1. Définir `outputDir` dans le fichier de configuration

Ajoutez le code suivant au fichier `nextflow.config` :

=== "Après"

    ```groovy title="nextflow.config" linenums="9" hl_lines="10-13"
    /*
    * Paramètres du pipeline
    */
    params {
        input = 'data/greetings.csv'
        batch = 'batch'
        character = 'turkey'
    }

    /*
    * Paramètres de sortie
    */
    outputDir = "custom-outdir-config/${params.batch}"
    ```

=== "Avant"

    ```groovy title="nextflow.config" linenums="9"
    /*
    * Paramètres du pipeline
    */
    params {
        input = 'data/greetings.csv'
        batch = 'batch'
        character = 'turkey'
    }
    ```

Cela définit le répertoire de sortie à `custom-outdir-config/` plus la valeur du paramètre `batch` comme sous-répertoire.
Maintenant vous pouvez changer l'emplacement de sortie en définissant le paramètre `--batch` :

```bash
nextflow run hello-config.nf --batch my_run
```

Cela publie les sorties dans `custom-outdir-config/my_run/`.

!!! note "Note"

    L'option CLI `-output-dir` a priorité sur le paramètre de configuration `outputDir`.
    Si elle est définie, l'option de configuration sera entièrement ignorée.

#### 2.2.2. Sous-répertoires avec noms de lot et de processus

Nous pouvons également définir des déclarations de `path` de sortie de sous-répertoire de manière dynamique, par sortie.

Par exemple, nous pouvons organiser nos sorties par processus en référençant `<processus>.name` dans la déclaration de chemin de sortie :

=== "Après"

    ```groovy title="hello-config.nf" linenums="42" hl_lines="3 7 11 15 19"
    output {
        first_output {
            path { sayHello.name }
            mode 'copy'
        }
        uppercased {
            path { convertToUpper.name }
            mode 'copy'
        }
        collected {
            path { collectGreetings.name }
            mode 'copy'
        }
        batch_report {
            path { collectGreetings.name }
            mode 'copy'
        }
        cowpy_art {
            path { cowpy.name }
            mode 'copy'
        }
    }
    ```

=== "Avant"

    ```groovy title="hello-config.nf" linenums="42" hl_lines="3 7 11 15 19"
    output {
        first_output {
            path 'intermediates'
            mode 'copy'
        }
        uppercased {
            path 'intermediates'
            mode 'copy'
        }
        collected {
            path 'intermediates'
            mode 'copy'
        }
        batch_report {
            path ''
            mode 'copy'
        }
        cowpy_art {
            path ''
            mode 'copy'
        }
    }
    ```

Nous pouvons aller plus loin et composer des chemins de sous-répertoires plus complexes.

Dans l'édition ci-dessus, nous avons effacé la distinction entre les `intermediates` versus les sorties finales étant au niveau supérieur.
Remettons cela, et plaçons également les fichiers dans un sous-répertoire `params.batch`.

!!! tip "Astuce"

    Inclure `params.batch` dans le `path` du bloc output, au lieu du `outputDir` de configuration, signifie qu'il ne sera pas écrasé avec `-output-dir` sur le CLI.

D'abord, mettez à jour le fichier de configuration pour supprimer `${params.batch}` du `outputDir` (puisque nous le déplaçons vers les déclarations de chemin) :

=== "Après"

    ```groovy title="nextflow.config" linenums="12" hl_lines="4"
    /*
    * Paramètres de sortie
    */
    outputDir = "custom-outdir-config/"
    ```

=== "Avant"

    ```groovy title="nextflow.config" linenums="12" hl_lines="4"
    /*
    * Paramètres de sortie
    */
    outputDir = "custom-outdir-config/${params.batch}"
    ```

Ensuite, effectuez les modifications suivantes dans le fichier de workflow :

=== "Après"

    ```groovy title="hello-config.nf" linenums="42" hl_lines="3 7 11 15 19"
    output {
        first_output {
            path { "${params.batch}/intermediates/${sayHello.name}" }
            mode 'copy'
        }
        uppercased {
            path { "${params.batch}/intermediates/${convertToUpper.name}" }
            mode 'copy'
        }
        collected {
            path { "${params.batch}/intermediates/${collectGreetings.name}" }
            mode 'copy'
        }
        batch_report {
            path { "${params.batch}/${collectGreetings.name}" }
            mode 'copy'
        }
        cowpy_art {
            path { "${params.batch}/${cowpy.name}" }
            mode 'copy'
        }
    }
    ```

=== "Avant"

    ```groovy title="hello-config.nf" linenums="42" hl_lines="3 7 11 15 19"
    output {
        first_output {
            path { sayHello.name }
            mode 'copy'
        }
        uppercased {
            path { convertToUpper.name }
            mode 'copy'
        }
        collected {
            path { collectGreetings.name }
            mode 'copy'
        }
        batch_report {
            path { collectGreetings.name }
            mode 'copy'
        }
        cowpy_art {
            path { cowpy.name }
            mode 'copy'
        }
    }
    ```

#### 2.2.3. Exécuter le pipeline

Voyons comment cela fonctionne en pratique, en définissant à la fois `-output-dir` (ou `-o` en abrégé) à `custom-outdir-config-2` et le nom du lot à `rep2` depuis la ligne de commande :

```bash
nextflow run hello-config.nf -output-dir custom-outdir-config-2 --batch rep2
```

??? success "Sortie de la commande"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `hello-config.nf` [mad_curry] DSL2 - revision: 668a98ccb9

    executor >  local (8)
    [9e/6095e0] sayHello (1)       [100%] 3 of 3 ✔
    [05/454d52] convertToUpper (3) [100%] 3 of 3 ✔
    [ed/e3ddfb] collectGreetings   [100%] 1 of 1 ✔
    [39/5e063a] cowpy              [100%] 1 of 1 ✔
    ```

Cela publie les sorties dans `custom-outdir-config-2/rep2/`, avec le chemin de base spécifié _et_ le sous-répertoire du nom de lot _et_ les résultats groupés par processus :

??? abstract "Contenu du répertoire"

    ```console
    custom-outdir-config-2
    └── rep2
        ├── collectGreetings
        │   └── rep2-report.txt
        ├── cowpy
        │   └── cowpy-COLLECTED-rep2-output.txt
        └── intermediates
            ├── collectGreetings
            │   └── COLLECTED-rep2-output.txt
            ├── convertToUpper
            │   ├── UPPER-Bonjour-output.txt
            │   ├── UPPER-Hello-output.txt
            │   └── UPPER-Holà-output.txt
            └── sayHello
                ├── Bonjour-output.txt
                ├── Hello-output.txt
                └── Holà-output.txt
    ```

### 2.3. Définir le mode de publication au niveau du workflow

Enfin, dans l'esprit de réduire la quantité de code répétitif, nous pouvons remplacer les déclarations `mode` par sortie par une seule ligne dans la configuration.

#### 2.3.1. Ajouter `workflow.output.mode` au fichier de configuration

Ajoutez le code suivant au fichier `nextflow.config` :

=== "Après"

    ```groovy title="nextflow.config" linenums="12" hl_lines="5"
    /*
    * Paramètres de sortie
    */
    outputDir = "custom-outdir-config/"
    workflow.output.mode = 'copy'
    ```

=== "Avant"

    ```groovy title="nextflow.config" linenums="12"
    /*
    * Paramètres de sortie
    */
    outputDir = "custom-outdir-config/"
    ```

Définir `workflow.output.mode` dans le fichier de configuration est suffisant pour remplacer ce qui est défini dans le fichier de workflow, mais supprimons quand même le code inutile.

#### 2.3.2. Supprimer le mode de sortie du fichier de workflow

Effectuez les modifications suivantes dans le fichier de workflow :

=== "Après"

    ```groovy title="hello-config.nf" linenums="42"
    output {
        first_output {
            path { "${params.batch}/intermediates/${sayHello.name}" }
        }
        uppercased {
            path { "${params.batch}/intermediates/${convertToUpper.name}" }
        }
        collected {
            path { "${params.batch}/intermediates/${collectGreetings.name}" }
        }
        batch_report {
            path { "${params.batch}/${collectGreetings.name}" }
        }
        cowpy_art {
            path { "${params.batch}/${cowpy.name}" }
        }
    }
    ```

=== "Avant"

    ```groovy title="hello-config.nf" linenums="42" hl_lines="4 8 12 16 20"
    output {
        first_output {
            path { "${params.batch}/intermediates/${sayHello.name}" }
            mode 'copy'
        }
        uppercased {
            path { "${params.batch}/intermediates/${convertToUpper.name}" }
            mode 'copy'
        }
        collected {
            path { "${params.batch}/intermediates/${collectGreetings.name}" }
            mode 'copy'
        }
        batch_report {
            path { "${params.batch}/${collectGreetings.name}" }
            mode 'copy'
        }
        cowpy_art {
            path { "${params.batch}/${cowpy.name}" }
            mode 'copy'
        }
    }
    ```

C'est plus concis, n'est-ce pas ?

#### 2.3.3. Exécuter le pipeline

Testons que cela fonctionne correctement :

```bash
nextflow run hello-config.nf -output-dir config-output-mode
```

??? success "Sortie de la commande"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `hello-config.nf` [small_stone] DSL2 - revision: 024d6361b5

    executor >  local (8)
    [e8/a0e93e] sayHello (1)       [100%] 3 of 3 ✔
    [14/176c9d] convertToUpper (3) [100%] 3 of 3 ✔
    [23/d667ca] collectGreetings   [100%] 1 of 1 ✔
    [e6/1dc80e] cowpy              [100%] 1 of 1 ✔
    ```

Cela publie les sorties dans `config-output-mode/`, et ce sont toujours toutes de vraies copies, pas des liens symboliques.

??? abstract "Contenu du répertoire"

    ```console
    config-output-mode
    └── batch
        ├── collectGreetings
        │   └── batch-report.txt
        ├── cowpy
        │   └── cowpy-COLLECTED-batch-output.txt
        └── intermediates
            ├── collectGreetings
            │   └── COLLECTED-batch-output.txt
            ├── convertToUpper
            │   ├── UPPER-Bonjour-output.txt
            │   ├── UPPER-Hello-output.txt
            │   └── UPPER-Holà-output.txt
            └── sayHello
                ├── Bonjour-output.txt
                ├── Hello-output.txt
                └── Holà-output.txt
    ```

La principale raison pour laquelle vous pourriez encore vouloir utiliser la façon par sortie de définir le mode est si vous voulez mélanger et assortir au sein du même workflow, _c.-à-d._ avoir certaines sorties copiées et d'autres liées symboliquement.

Il y a plein d'autres options que vous pouvez personnaliser de cette façon, mais espérons que cela vous donne une idée de l'éventail d'options et comment les utiliser efficacement selon vos préférences.

### À retenir

Vous savez comment contrôler le nommage et la structure des répertoires où vos sorties sont publiées, ainsi que le mode de publication des sorties du workflow.

### Et ensuite ?

Apprenez à adapter la configuration de votre workflow à votre environnement de calcul, en commençant par la technologie de packaging logiciel.

---

## 3. Sélectionner une technologie de packaging logiciel

Jusqu'à présent, nous avons examiné des éléments de configuration qui contrôlent comment les entrées entrent et où les sorties sortent. Maintenant, il est temps de nous concentrer plus spécifiquement sur l'adaptation de la configuration de votre workflow à votre environnement de calcul.

La première étape sur ce chemin est de spécifier d'où proviendront les packages logiciels qui seront exécutés à chaque étape.
Sont-ils déjà installés dans l'environnement de calcul local ?
Devons-nous récupérer des images et les exécuter via un système de conteneurs ?
Ou devons-nous récupérer des packages Conda et construire un environnement Conda local ?

Dans la toute première partie de ce cours de formation (Parties 1-4), nous avons simplement utilisé des logiciels installés localement dans notre workflow.
Puis dans la Partie 5, nous avons introduit les conteneurs Docker et le fichier `nextflow.config`, que nous avons utilisé pour activer l'utilisation des conteneurs Docker.

Maintenant, voyons comment nous pouvons configurer une option de packaging logiciel alternative via le fichier `nextflow.config`.

### 3.1. Désactiver Docker et activer Conda dans le fichier de configuration

Prétendons que nous travaillons sur un cluster HPC et que l'administrateur·trice n'autorise pas l'utilisation de Docker pour des raisons de sécurité.
Heureusement pour nous, Nextflow prend en charge plusieurs autres technologies de conteneurs comme Singularity (qui est plus largement utilisé sur HPC), et des gestionnaires de packages logiciels tels que Conda.

Nous pouvons changer notre fichier de configuration pour utiliser [Conda](https://nextflow.io/docs/latest/conda.html) au lieu de Docker.
Pour ce faire, changeons la valeur de `docker.enabled` à `false`, et ajoutons une directive activant l'utilisation de Conda :

=== "Après"

    ```groovy title="nextflow.config" linenums="1" hl_lines="1-2"
    docker.enabled = false
    conda.enabled = true
    ```

=== "Avant"

    ```groovy title="nextflow.config" linenums="1" hl_lines="1"
    docker.enabled = true
    ```

Cela permettra à Nextflow de créer et d'utiliser des environnements Conda pour les processus qui ont des packages Conda spécifiés.
Ce qui signifie que nous devons maintenant en ajouter un à notre processus `cowpy` !

### 3.2. Spécifier un package Conda dans la définition du processus

Nous avons déjà récupéré l'URI pour un package Conda contenant l'outil `cowpy` : `conda-forge::cowpy==1.1.5`

Maintenant nous ajoutons l'URI à la définition du processus `cowpy` en utilisant la directive `conda` :

=== "Après"

    ```groovy title="modules/cowpy.nf" linenums="4" hl_lines="4"
    process cowpy {

        container 'community.wave.seqera.io/library/cowpy:1.1.5--3db457ae1977a273'
        conda 'conda-forge::cowpy==1.1.5'

        input:
    ```

=== "Avant"

    ```groovy title="modules/cowpy.nf" linenums="4"
    process cowpy {

        container 'community.wave.seqera.io/library/cowpy:1.1.5--3db457ae1977a273'

        input:
    ```

Pour être clair, nous ne _remplaçons_ pas la directive `docker`, nous _ajoutons_ une option alternative.

!!! tip "Astuce"

    Il existe plusieurs façons d'obtenir l'URI pour un package conda donné.
    Nous recommandons d'utiliser la requête de recherche [Seqera Containers](https://seqera.io/containers/), qui vous donnera un URI que vous pouvez copier et coller, même si vous ne prévoyez pas de créer un conteneur à partir de celui-ci.

### 3.3. Exécuter le workflow pour vérifier qu'il peut utiliser Conda

Essayons.

```bash
nextflow run hello-config.nf --batch conda
```

??? success "Sortie de la commande"

    ```console title="Sortie"
    N E X T F L O W   ~  version 25.10.2

    Launching `hello-config.nf` [friendly_lamport] DSL2 - revision: 024d6361b5

    executor >  local (8)
    [e8/91c116] sayHello (2)       [100%] 3 of 3 ✔
    [fe/6a70ce] convertToUpper (3) [100%] 3 of 3 ✔
    [99/7cc493] collectGreetings   [100%] 1 of 1 ✔
    [3c/09fb59] cowpy              [100%] 1 of 1 ✔
    ```

Cela devrait fonctionner sans problème et produire les mêmes sorties qu'avant sous `custom-outdir-config/conda`.

En coulisses, Nextflow a récupéré les packages Conda et créé l'environnement, ce qui prend normalement un peu de travail ; c'est donc agréable que nous n'ayons pas à faire tout cela nous-mêmes !

!!! note "Note"

    Cela s'exécute rapidement parce que le package `cowpy` est assez petit, mais si vous travaillez avec de gros packages, cela peut prendre un peu plus de temps que d'habitude la première fois, et vous pourriez voir la sortie de la console rester « coincée » pendant une minute environ avant de se terminer.
    C'est normal et c'est dû au travail supplémentaire que fait Nextflow la première fois que vous utilisez un nouveau package.

De notre point de vue, cela semble fonctionner exactement de la même manière qu'avec Docker, même si en coulisses la mécanique est un peu différente.

Cela signifie que nous sommes prêts à exécuter avec des environnements Conda si nécessaire.

??? info "Mélanger et assortir Docker et Conda"

    Puisque ces directives sont assignées par processus, il est possible de « mélanger et assortir », _c.-à-d._ configurer certains des processus de votre workflow pour s'exécuter avec Docker et d'autres avec Conda, par exemple, si l'infrastructure de calcul que vous utilisez prend en charge les deux.
    Dans ce cas, vous activeriez à la fois Docker et Conda dans votre fichier de configuration.
    Si les deux sont disponibles pour un processus donné, Nextflow priorisera les conteneurs.

    Et comme noté précédemment, Nextflow prend en charge plusieurs autres technologies de packaging logiciel et de conteneurs, donc vous n'êtes pas limité·e à ces deux-là.

### À retenir

Vous savez comment configurer quel package logiciel chaque processus doit utiliser, et comment basculer entre les technologies.

### Et ensuite ?

Apprenez à changer la plateforme d'exécution utilisée par Nextflow pour effectivement faire le travail.

---

## 4. Sélectionner une plateforme d'exécution

Jusqu'à présent, nous avons exécuté notre pipeline avec l'exécuteur local.
Celui-ci exécute chaque tâche sur la machine sur laquelle Nextflow s'exécute.
Quand Nextflow démarre, il regarde les CPU et la mémoire disponibles.
Si les ressources des tâches prêtes à s'exécuter dépassent les ressources disponibles, Nextflow retiendra les dernières tâches de l'exécution jusqu'à ce qu'une ou plusieurs des tâches précédentes soient terminées, libérant les ressources nécessaires.

L'exécuteur local est pratique et efficace, mais il est limité à cette seule machine. Pour de très grandes charges de travail, vous pourriez découvrir que votre machine locale est un goulot d'étranglement, soit parce que vous avez une seule tâche qui nécessite plus de ressources que ce que vous avez disponible, soit parce que vous avez tellement de tâches qu'attendre qu'une seule machine les exécute prendrait trop de temps.

Nextflow prend en charge [de nombreux backends d'exécution différents](https://nextflow.io/docs/latest/executor.html), y compris les ordonnanceurs HPC (Slurm, LSF, SGE, PBS, Moab, OAR, Bridge, HTCondor et autres) ainsi que les backends d'exécution cloud tels que (AWS Batch, Google Cloud Batch, Azure Batch, Kubernetes et plus).

### 4.1. Cibler un backend différent

Le choix de l'exécuteur est défini par une directive de processus appelée `executor`.
Par défaut, elle est définie sur `local`, donc la configuration suivante est implicite :

```groovy title="Configuration intégrée"
process {
    executor = 'local'
}
```

Pour définir l'exécuteur pour cibler un backend différent, vous spécifieriez simplement l'exécuteur que vous voulez en utilisant une syntaxe similaire à celle décrite ci-dessus pour les allocations de ressources (voir la [documentation des exécuteurs](https://nextflow.io/docs/latest/executor.html) pour toutes les options).

```groovy title="nextflow.config"
process {
    executor = 'slurm'
}
```

!!! warning "Avertissement"

    Nous ne pouvons pas réellement tester ceci dans l'environnement de formation car il n'est pas configuré pour se connecter à un HPC.

### 4.2. Gérer la syntaxe spécifique au backend pour les paramètres d'exécution

La plupart des plateformes de calcul haute performance permettent (et parfois exigent) que vous spécifiiez certains paramètres tels que les demandes et limitations d'allocation de ressources (par exemple, nombre de CPU et mémoire) et le nom de la file d'attente de jobs à utiliser.

Malheureusement, chacun de ces systèmes utilise des technologies, syntaxes et configurations différentes pour définir comment un job doit être défini et soumis à l'ordonnanceur correspondant.

??? abstract "Exemples"

    Par exemple, le même job nécessitant 8 CPU et 4 Go de RAM à exécuter sur la file « my-science-work » doit être exprimé de différentes manières selon le backend.

    ```bash title="Config pour SLURM / soumettre avec sbatch"
    #SBATCH -o /path/to/my/task/directory/my-task-1.log
    #SBATCH --no-requeue
    #SBATCH -c 8
    #SBATCH --mem 4096M
    #SBATCH -p my-science-work
    ```

    ```bash title="Config pour PBS / soumettre avec qsub"
    #PBS -o /path/to/my/task/directory/my-task-1.log
    #PBS -j oe
    #PBS -q my-science-work
    #PBS -l nodes=1:ppn=5
    #PBS -l mem=4gb
    ```

    ```bash title="Config pour SGE / soumettre avec qsub"
    #$ -o /path/to/my/task/directory/my-task-1.log
    #$ -j y
    #$ -terse
    #$ -notify
    #$ -q my-science-work
    #$ -l slots=5
    #$ -l h_rss=4096M,mem_free=4096M
    ```

Heureusement, Nextflow simplifie tout cela.
Il fournit une syntaxe standardisée pour que vous puissiez spécifier les propriétés pertinentes telles que [`cpus`](https://nextflow.io/docs/latest/reference/process.html#cpus), [`memory`](https://nextflow.io/docs/latest/reference/process.html#memory) et [`queue`](https://nextflow.io/docs/latest/reference/process.html#queue) (voir les [directives de processus](https://nextflow.io/docs/latest/reference/process.html#process-directives) pour d'autres propriétés) une seule fois.
Ensuite, à l'exécution, Nextflow utilisera ces paramètres pour générer les scripts appropriés spécifiques au backend en fonction du paramètre d'exécuteur.

Nous couvrirons cette syntaxe standardisée dans la section suivante.

### À retenir

Vous savez maintenant comment changer l'exécuteur pour utiliser différents types d'infrastructure de calcul.

### Et ensuite ?

Apprenez à évaluer et exprimer les allocations et limitations de ressources dans Nextflow.

---

## 5. Contrôler les allocations de ressources de calcul

La plupart des plateformes de calcul haute performance permettent (et parfois exigent) que vous spécifiiez certains paramètres d'allocation de ressources tels que le nombre de CPU et la mémoire.

Par défaut, Nextflow utilisera un seul CPU et 2 Go de mémoire pour chaque processus.
Les directives de processus correspondantes sont appelées `cpus` et `memory`, donc la configuration suivante est implicite :

```groovy title="Configuration intégrée" linenums="1"
process {
    cpus = 1
    memory = 2.GB
}
```

Vous pouvez modifier ces valeurs, soit pour tous les processus, soit pour des processus nommés spécifiques, en utilisant des directives de processus supplémentaires dans votre fichier de configuration.
Nextflow les traduira en instructions appropriées pour l'exécuteur choisi.

Mais comment savez-vous quelles valeurs utiliser ?

### 5.1. Exécuter le workflow pour générer un rapport d'utilisation des ressources

Si vous ne savez pas à l'avance combien de CPU et de mémoire vos processus sont susceptibles d'avoir besoin, vous pouvez faire du profilage de ressources, ce qui signifie que vous exécutez le workflow avec des allocations par défaut, enregistrez combien chaque processus a utilisé, et à partir de là, estimez comment ajuster les allocations de base.

De manière pratique, Nextflow inclut des outils intégrés pour faire ceci, et générera volontiers un rapport pour vous sur demande.

Pour ce faire, ajoutez `-with-report <nom_de_fichier>.html` à votre ligne de commande.

```bash
nextflow run hello-config.nf -with-report report-config-1.html
```

Le rapport est un fichier html, que vous pouvez télécharger et ouvrir dans votre navigateur. Vous pouvez également faire un clic droit dessus dans l'explorateur de fichiers à gauche et cliquer sur `Show preview` pour le visualiser dans l'environnement de formation.

Prenez quelques minutes pour parcourir le rapport et voir si vous pouvez identifier des opportunités d'ajustement des ressources.
Assurez-vous de cliquer sur les onglets qui montrent les résultats d'utilisation en pourcentage de ce qui a été alloué.

Voir [Rapports](https://nextflow.io/docs/latest/reports.html) pour la documentation sur toutes les fonctionnalités disponibles.

### 5.2. Définir les allocations de ressources pour tous les processus

Le profilage montre que les processus de notre workflow de formation sont très légers, donc réduisons l'allocation de mémoire par défaut à 1 Go par processus.

Ajoutez ce qui suit à votre fichier `nextflow.config`, avant la section des paramètres du pipeline :

=== "Après"

    ```groovy title="nextflow.config" linenums="1" hl_lines="4-9"
    docker.enabled = false
    conda.enabled = true

    /*
    * Paramètres de processus
    */
    process {
        memory = 1.GB
    }

    /*
    * Paramètres du pipeline
    */
    params {
        input = 'data/greetings.csv'
        batch = 'batch'
        character = 'turkey'
    }
    ```

=== "Avant"

    ```groovy title="nextflow.config" linenums="1"
    docker.enabled = false
    conda.enabled = true

    /*
    * Paramètres du pipeline
    */
    params {
        input = 'data/greetings.csv'
        batch = 'batch'
        character = 'turkey'
    }
    ```

Cela aidera à réduire la quantité de calcul que nous consommons.

### 5.3. Définir les allocations de ressources pour un processus spécifique

En même temps, nous allons prétendre que le processus `cowpy` nécessite plus de ressources que les autres, juste pour pouvoir démontrer comment ajuster les allocations pour un processus individuel.

=== "Après"

    ```groovy title="nextflow.config" linenums="4" hl_lines="6-9"
    /*
    * Paramètres de processus
    */
    process {
        memory = 1.GB
        withName: 'cowpy' {
            memory = 2.GB
            cpus = 2
        }
    }
    ```

=== "Avant"

    ```groovy title="nextflow.config" linenums="4"
    /*
    * Paramètres de processus
    */
    process {
        memory = 1.GB
    }
    ```

Avec cette configuration, tous les processus demanderont 1 Go de mémoire et un seul CPU (la valeur par défaut implicite), sauf le processus `cowpy`, qui demandera 2 Go et 2 CPU.

!!! tip "Astuce"

    Si vous avez une machine avec peu de CPU et que vous en allouez un nombre élevé par processus, vous pourriez voir des appels de processus mis en file d'attente derrière d'autres.
    C'est parce que Nextflow s'assure que nous ne demandons pas plus de CPU qu'il n'y en a de disponibles.

### 5.4. Exécuter le workflow avec la configuration mise à jour

Essayons cela, en fournissant un nom de fichier différent pour le rapport de profilage afin que nous puissions comparer les performances avant et après les changements de configuration.

```bash
nextflow run hello-config.nf -with-report report-config-2.html
```

Vous ne remarquerez probablement aucune différence réelle puisque c'est une si petite charge de travail, mais c'est l'approche que vous utiliseriez pour analyser les performances et les besoins en ressources d'un workflow du monde réel.

C'est très utile lorsque vos processus ont des besoins en ressources différents. Cela vous permet de dimensionner correctement les allocations de ressources que vous configurez pour chaque processus en fonction de données réelles, pas de suppositions.

!!! tip "Astuce"

    Ce n'est qu'un petit aperçu de ce que vous pouvez faire pour optimiser votre utilisation des ressources.
    Nextflow lui-même a une [logique de réessai dynamique](https://nextflow.io/docs/latest/process.html#dynamic-task-resources) vraiment chouette intégrée pour réessayer les jobs qui échouent en raison de limitations de ressources.
    De plus, la Seqera Platform offre des outils pilotés par l'IA pour optimiser automatiquement vos allocations de ressources également.

### 5.5. Ajouter des limites de ressources

Selon l'exécuteur de calcul et l'infrastructure de calcul que vous utilisez, il peut y avoir des contraintes sur ce que vous pouvez (ou devez) allouer.
Par exemple, votre cluster peut exiger que vous restiez dans certaines limites.

Vous pouvez utiliser la directive `resourceLimits` pour définir les limitations pertinentes. La syntaxe ressemble à ceci quand elle est seule dans un bloc process :

```groovy title="Exemple de syntaxe"
process {
    resourceLimits = [
        memory: 750.GB,
        cpus: 200,
        time: 30.d
    ]
}
```

Nextflow traduira ces valeurs en instructions appropriées selon l'exécuteur que vous avez spécifié.

Nous n'allons pas exécuter ceci, puisque nous n'avons pas accès à l'infrastructure pertinente dans l'environnement de formation.
Cependant, si vous essayiez d'exécuter le workflow avec des allocations de ressources qui dépassent ces limites, puis que vous regardiez la commande `sbatch` dans le fichier de script `.command.run`, vous verriez que les demandes qui sont réellement envoyées à l'exécuteur sont plafonnées aux valeurs spécifiées par `resourceLimits`.

??? info "Configurations de référence institutionnelles"

    Le projet nf-core a compilé une [collection de fichiers de configuration](https://nf-co.re/configs/) partagés par diverses institutions à travers le monde, couvrant un large éventail d'exécuteurs HPC et cloud.

    Ces configurations partagées sont précieuses à la fois pour les personnes qui travaillent là-bas et peuvent donc simplement utiliser la configuration de leur institution prête à l'emploi, et comme modèle pour les personnes qui cherchent à développer une configuration pour leur propre infrastructure.

### À retenir

Vous savez comment générer un rapport de profilage pour évaluer l'utilisation des ressources et comment modifier les allocations de ressources pour tous les processus et/ou pour des processus individuels, ainsi que définir des limitations de ressources pour l'exécution sur HPC.

### Et ensuite ?

Apprenez à configurer des profils de configuration prédéfinis et à basculer entre eux à l'exécution.

---

## 6. Utiliser des profils pour basculer entre des configurations prédéfinies

Nous vous avons montré plusieurs façons de personnaliser la configuration de votre pipeline selon le projet sur lequel vous travaillez ou l'environnement de calcul que vous utilisez.

Vous pourriez vouloir basculer entre des paramètres alternatifs selon l'infrastructure de calcul que vous utilisez. Par exemple, vous pourriez vouloir développer et exécuter des tests à petite échelle localement sur votre ordinateur portable, puis exécuter des charges de travail à grande échelle sur HPC ou cloud.

Nextflow vous permet de configurer n'importe quel nombre de [profils](https://nextflow.io/docs/latest/config.html#config-profiles) qui décrivent différentes configurations, que vous pouvez ensuite sélectionner à l'exécution en utilisant un argument de ligne de commande, plutôt que d'avoir à modifier le fichier de configuration lui-même.

### 6.1. Créer des profils pour basculer entre le développement local et l'exécution sur HPC

Configurons deux profils alternatifs ; un pour exécuter de petites charges sur un ordinateur ordinaire, où nous utiliserons des conteneurs Docker, et un pour l'exécution sur un HPC universitaire avec un ordonnanceur Slurm, où nous utiliserons des packages Conda.

#### 6.1.1. Configurer les profils

Ajoutez ce qui suit à votre fichier `nextflow.config`, après la section des paramètres du pipeline mais avant les paramètres de sortie :

=== "Après"

    ```groovy title="nextflow.config" linenums="15" hl_lines="10-27"
    /*
    * Paramètres du pipeline
    */
    params {
        input = 'data/greetings.csv'
        batch = 'batch'
        character = 'turkey'
    }

    /*
    * Profils
    */
    profiles {
        my_laptop {
            process.executor = 'local'
            docker.enabled = true
        }
        univ_hpc {
            process.executor = 'slurm'
            conda.enabled = true
            process.resourceLimits = [
                memory: 750.GB,
                cpus: 200,
                time: 30.d
            ]
        }
    }

    /*
    * Paramètres de sortie
    */
    outputDir = "custom-outdir-config/"
    workflow.output.mode = 'copy'
    ```

=== "Avant"

    ```groovy title="nextflow.config" linenums="15"
    /*
    * Paramètres du pipeline
    */
    params {
        input = 'data/greetings.csv'
        batch = 'batch'
        character = 'turkey'
    }

    /*
    * Paramètres de sortie
    */
    outputDir = "custom-outdir-config/"
    workflow.output.mode = 'copy'
    ```

Vous voyez que pour le HPC universitaire, nous spécifions également des limitations de ressources.

#### 6.1.2. Exécuter le workflow avec un profil

Pour spécifier un profil dans notre ligne de commande Nextflow, nous utilisons l'argument `-profile`.

Essayons d'exécuter le workflow avec la configuration `my_laptop`.

```bash
nextflow run hello-config.nf -profile my_laptop
```

??? success "Sortie de la commande"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `hello-config.nf` [hungry_sanger] DSL2 - revision: 024d6361b5

    executor >  local (8)
    [b0/fb2ec9] sayHello (3)       [100%] 3 of 3 ✔
    [4a/e039f0] convertToUpper (3) [100%] 3 of 3 ✔
    [6f/408fa9] collectGreetings   [100%] 1 of 1 ✔
    [f1/fd6520] cowpy              [100%] 1 of 1 ✔
    ```

Comme vous pouvez le voir, cela nous permet de basculer entre les configurations très facilement à l'exécution.

!!! warning "Avertissement"

    Le profil `univ_hpc` ne fonctionnera pas correctement dans l'environnement de formation puisque nous n'avons pas accès à un ordonnanceur Slurm.

Si à l'avenir nous trouvons d'autres éléments de configuration qui co-surviennent toujours avec ceux-ci, nous pouvons simplement les ajouter au(x) profil(s) correspondant(s).
Nous pouvons également créer des profils supplémentaires s'il y a d'autres éléments de configuration que nous voulons regrouper.

### 6.2. Créer un profil de paramètres de test

Les profils ne sont pas seulement pour la configuration de l'infrastructure.
Nous pouvons également les utiliser pour définir des valeurs par défaut pour les paramètres du workflow, pour faciliter à d'autres d'essayer le workflow sans avoir à rassembler eux·elles-mêmes des valeurs d'entrée appropriées.
Vous pouvez considérer cela comme une alternative à l'utilisation d'un fichier de paramètres.

#### 6.2.1. Configurer le profil

La syntaxe pour exprimer les valeurs par défaut dans ce contexte ressemble à ceci, pour un profil que nous nommons `test` :

```groovy title="Exemple de syntaxe"
    test {
        params.<parameter1>
        params.<parameter2>
        ...
    }
```

Si nous ajoutons un profil de test pour notre workflow, le bloc `profiles` devient :

```groovy title="nextflow.config" linenums="24" hl_lines="18-22"
/*
* Profils
*/
profiles {
    my_laptop {
        process.executor = 'local'
        docker.enabled = true
    }
    univ_hpc {
        process.executor = 'slurm'
        conda.enabled = true
        process.resourceLimits = [
            memory: 750.GB,
            cpus: 200,
            time: 30.d
        ]
    }
    test {
        params.input = 'data/greetings.csv'
        params.batch = 'test'
        params.character = 'dragonandcow'
    }
}
```

Tout comme pour les profils de configuration technique, vous pouvez configurer plusieurs profils différents spécifiant des paramètres sous n'importe quel nom arbitraire que vous aimez.

#### 6.2.2. Exécuter le workflow localement avec le profil de test

De manière pratique, les profils ne sont pas mutuellement exclusifs, donc nous pouvons spécifier plusieurs profils dans notre ligne de commande en utilisant la syntaxe suivante `-profile <profile1>,<profile2>` (pour n'importe quel nombre de profils).

Si vous combinez des profils qui définissent des valeurs pour les mêmes éléments de configuration et sont décrits dans le même fichier de configuration, Nextflow résoudra le conflit en utilisant la valeur qu'il a lue en dernier (_c.-à-d._ ce qui vient plus tard dans le fichier).
Si les paramètres conflictuels sont définis dans différentes sources de configuration, l'[ordre de priorité](https://nextflow.io/docs/latest/config.html) par défaut s'applique.

Essayons d'ajouter le profil de test à notre commande précédente :

```bash
nextflow run hello-config.nf -profile my_laptop,test
```

??? success "Sortie de la commande"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `hello-config.nf` [modest_becquerel] DSL2 - revision: 024d6361b5

    executor >  local (8)
    [4c/fe2580] sayHello (1)       [100%] 3 of 3 ✔
    [fd/7d9017] convertToUpper (3) [100%] 3 of 3 ✔
    [13/1523bd] collectGreetings   [100%] 1 of 1 ✔
    [06/a1ee14] cowpy              [100%] 1 of 1 ✔
    ```

Cela utilisera Docker quand c'est possible et produira des sorties sous `custom-outdir-config/test`, et cette fois le personnage est le duo comique `dragonandcow`.

??? abstract "Contenu du fichier"

    ```console title="custom-outdir-config/test/cowpy/cowpy-COLLECTED-test-output.txt"
     _________
    / HOLà    \
    | HELLO   |
    \ BONJOUR /
    ---------
                \                    ^    /^
                  \                  / \  // \
                  \   |\___/|      /   \//  .\
                    \  /O  O  \__  /    //  | \ \           *----*
                      /     /  \/_/    //   |  \  \          \   |
                      \@___\@`    \/_   //    |   \   \         \/\ \
                    0/0/|       \/_ //     |    \    \         \ \
                0/0/0/0/|        \///      |     \     \       | |
              0/0/0/0/0/_|_ /   (  //       |      \     _\     |  /
          0/0/0/0/0/0/`/,_ _ _/  ) ; -.    |    _ _\.-~       /   /
                      ,-}        _      *-.|.-~-.           .~    ~
      \     \__/        `/\      /                 ~-. _ .-~      /
      \____(oo)           *.   }            {                   /
      (    (--)          .----~-.\        \-`                 .~
      //__\\  \__ Ack!   ///.----..<        \             _ -~
      //    \\               ///-._ _ _ _ _ _ _{^ - - - - ~
    ```

Cela signifie que tant que nous distribuons tous les fichiers de données de test avec le code du workflow, n'importe qui peut rapidement essayer le workflow sans avoir à fournir ses propres entrées via la ligne de commande ou un fichier de paramètres.

!!! tip "Astuce"

    Nous pouvons pointer vers des URLs pour des fichiers plus volumineux qui sont stockés à l'extérieur.
    Nextflow les téléchargera automatiquement tant qu'il y a une connexion ouverte.

    Pour plus de détails, voir la Quête secondaire [Working with Files](../side_quests/working_with_files.md)

### 6.3. Utiliser `nextflow config` pour voir la configuration résolue

Comme noté ci-dessus, parfois le même paramètre peut être défini à des valeurs différentes dans des profils que vous voulez combiner.
Et plus généralement, il existe de nombreux endroits où des éléments de configuration peuvent être stockés, et parfois les mêmes propriétés peuvent être définies à des valeurs différentes à différents endroits.

Nextflow applique un [ordre de priorité](https://nextflow.io/docs/latest/config.html) défini pour résoudre les conflits, mais cela peut être difficile à déterminer par vous-même.
Et même si rien n'est en conflit, il peut être fastidieux de chercher tous les endroits possibles où les choses pourraient être configurées.

Heureusement, Nextflow inclut un outil utilitaire pratique appelé `config` qui peut automatiser tout ce processus pour vous.

L'outil `config` explorera tout le contenu de votre répertoire de travail actuel, aspirera tous les fichiers de configuration, et produira la configuration entièrement résolue que Nextflow utiliserait pour exécuter le workflow.
Cela vous permet de savoir quels paramètres seront utilisés sans avoir à lancer quoi que ce soit.

#### 6.3.1. Résoudre la configuration par défaut

Exécutez cette commande pour résoudre la configuration qui serait appliquée par défaut.

```bash
nextflow config
```

??? success "Sortie de la commande"

    ```groovy
    params {
      input = 'data/greetings.csv'
      batch = 'batch'
      character = 'turkey'
    }

    docker {
      enabled = false
    }

    conda {
      enabled = true
    }

    process {
      memory = '1 GB'
      withName:cowpy {
          memory = '2 GB'
          cpus = 2
      }
    }

    outputDir = 'custom-outdir-config/'

    workflow {
      output {
          mode = 'copy'
      }
    }
    ```

Cela vous montre la configuration de base que vous obtenez si vous ne spécifiez rien de plus dans la ligne de commande.

#### 6.3.2. Résoudre la configuration avec des paramètres spécifiques activés

Si vous fournissez des paramètres de ligne de commande, par exemple en activant un ou plusieurs profils ou en chargeant un fichier de paramètres, la commande les prendra également en compte.

```bash
nextflow config -profile my_laptop,test
```

??? success "Sortie de la commande"

    ```groovy
    params {
      input = 'data/greetings.csv'
      batch = 'test'
      character = 'dragonandcow'
    }

    docker {
      enabled = true
    }

    conda {
      enabled = true
    }

    process {
      memory = '1 GB'
      withName:cowpy {
          memory = '2 GB'
          cpus = 2
      }
      executor = 'local'
    }

    outputDir = 'custom-outdir-config/'

    workflow {
      output {
          mode = 'copy'
      }
    }
    ```

Cela devient particulièrement utile pour des projets complexes qui impliquent plusieurs couches de configuration.

### À retenir

Vous savez comment utiliser des profils pour sélectionner une configuration prédéfinie à l'exécution avec un minimum de tracas.
Plus généralement, vous savez comment configurer les exécutions de votre workflow pour s'adapter à différentes plateformes de calcul et améliorer la reproductibilité de vos analyses.

### Et ensuite ?

Célébrez et donnez-vous une grande tape dans le dos ! Vous avez terminé votre tout premier cours de développeur·se Nextflow.

Rendez-vous au [résumé final du cours](./next_steps.md) pour revoir ce que vous avez appris et découvrir ce qui vient ensuite.

---

## Quiz

<quiz>
Quel est le nom du fichier de configuration que Nextflow charge automatiquement ?
- [ ] `config.nf`
- [ ] `pipeline.config`
- [x] `nextflow.config`
- [ ] `workflow.config`
</quiz>

<quiz>
Qu'est-ce qui a priorité lorsque le même paramètre est défini à la fois dans le fichier de configuration et sur la ligne de commande ?
- [ ] La valeur du fichier de configuration
- [x] La valeur de la ligne de commande
- [ ] La première valeur rencontrée
- [ ] Aucune ; cela provoque une erreur

En savoir plus : [1.1. Déplacer les valeurs par défaut vers `nextflow.config`](#11-deplacer-les-valeurs-par-defaut-vers-nextflowconfig)
</quiz>

<quiz>
Pouvez-vous avoir à la fois Docker et Conda activés dans la même configuration ?
- [x] Oui, Nextflow peut utiliser les deux selon les directives de processus
- [ ] Non, un seul peut être activé à la fois
- [ ] Oui, mais seulement dans les profils
- [ ] Non, ils sont mutuellement exclusifs
</quiz>

<quiz>
Si Docker et Conda sont tous deux activés et qu'un processus a les deux directives, lequel est priorisé ?
- [x] Docker (conteneurs)
- [ ] Conda
- [ ] Le premier défini
- [ ] Cela provoque une erreur

En savoir plus : [3. Sélectionner une technologie de packaging logiciel](#3-selectionner-une-technologie-de-packaging-logiciel)
</quiz>

<quiz>
Quelle est l'allocation de mémoire par défaut pour les processus Nextflow ?
- [ ] 1 Go
- [x] 2 Go
- [ ] 4 Go
- [ ] Pas de limite
</quiz>

<quiz>
Comment définissez-vous les besoins en ressources pour un processus spécifique dans le fichier de configuration ?
- [ ] `#!groovy processName.memory = '4 GB'`
- [ ] `#!groovy process.memory.processName = '4 GB'`
- [x] `#!groovy process { withName: 'processName' { memory = '4 GB' } }`
- [ ] `#!groovy resources.processName.memory = '4 GB'`

En savoir plus : [5.3. Définir les allocations de ressources pour un processus spécifique](#53-definir-les-allocations-de-ressources-pour-un-processus-specifique)
</quiz>

<quiz>
Quelle option de ligne de commande génère un rapport d'utilisation des ressources ?
- [ ] `-with-metrics`
- [ ] `-with-stats`
- [x] `-with-report`
- [ ] `-with-profile`

En savoir plus : [5.1. Exécuter le workflow pour générer un rapport d'utilisation des ressources](#51-executer-le-workflow-pour-generer-un-rapport-dutilisation-des-ressources)
</quiz>

<quiz>
Que fait la directive `resourceLimits` ?
- [ ] Définit les besoins minimaux en ressources
- [ ] Alloue des ressources aux processus
- [x] Plafonne le maximum de ressources qui peuvent être demandées
- [ ] Surveille l'utilisation des ressources

En savoir plus : [5.5. Ajouter des limites de ressources](#55-ajouter-des-limites-de-ressources)
</quiz>

<quiz>
Quel est l'exécuteur par défaut dans Nextflow ?
- [x] `local`
- [ ] `slurm`
- [ ] `kubernetes`
- [ ] `aws`

En savoir plus : [4. Sélectionner une plateforme d'exécution](#4-selectionner-une-plateforme-dexecution)
</quiz>

<quiz>
Comment spécifiez-vous un fichier de paramètres lors de l'exécution de Nextflow ?
- [ ] `--params params.json`
- [ ] `-config params.json`
- [x] `-params-file params.json`
- [ ] `--input params.json`

En savoir plus : [1.3. Utiliser un fichier de paramètres](#13-utiliser-un-fichier-de-parametres)
</quiz>

<quiz>
À quoi peuvent servir les profils ? (Sélectionnez toutes les réponses applicables)
- [x] Définir des paramètres spécifiques à l'infrastructure
- [x] Définir des limites de ressources pour différents environnements
- [x] Fournir des paramètres de test
- [ ] Définir de nouveaux processus

En savoir plus : [6. Utiliser des profils pour basculer entre des configurations prédéfinies](#6-utiliser-des-profils-pour-basculer-entre-des-configurations-predefinies)
</quiz>

<quiz>
Comment spécifiez-vous plusieurs profils dans une seule commande ?
- [ ] `-profile profile1 -profile profile2`
- [ ] `-profiles profile1,profile2`
- [x] `-profile profile1,profile2`
- [ ] `--profile profile1 --profile profile2`

En savoir plus : [6. Utiliser des profils pour basculer entre des configurations prédéfinies](#6-utiliser-des-profils-pour-basculer-entre-des-configurations-predefinies)
</quiz>
