# Partie 3 : Configuration d'exécution

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Traduction assistée par IA - [en savoir plus et suggérer des améliorations](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

Cette section explorera comment gérer la configuration d'un pipeline Nextflow afin de personnaliser son comportement, l'adapter à différents environnements et optimiser l'utilisation des ressources _sans modifier une seule ligne du code du workflow lui-même_.

Il existe plusieurs façons de le faire, qui peuvent être utilisées en combinaison et sont interprétées selon l'ordre de priorité décrit dans la documentation [Configuration](https://nextflow.io/docs/latest/config.html).

Dans cette partie de la formation, nous allons vous montrer le mécanisme de fichier de configuration le plus simple et le plus courant, le fichier `nextflow.config`, que vous avez déjà rencontré dans la section sur les conteneurs de la Partie 2.

Nous passerons en revue les composants essentiels de la configuration Nextflow tels que les directives de process, les executors, les profils et les fichiers de paramètres.
En apprenant à utiliser efficacement ces options de configuration, vous pourrez tirer pleinement parti de la flexibilité, de l'évolutivité et des performances des pipelines Nextflow.

Pour exercer ces éléments de configuration, nous allons exécuter une nouvelle copie du workflow que nous avons exécuté à la fin de la Partie 2 de cette formation, renommé `3-main.nf`.

Si vous n'êtes pas familier·ère avec le pipeline Hello ou si vous avez besoin d'un rappel, consultez [cette page d'information](../info/hello_pipeline.md).

---

## 1. Gérer les paramètres d'entrée du workflow

??? example "Scénario"

    Vous avez téléchargé un pipeline et souhaitez l'exécuter de manière répétée avec les mêmes fichiers d'entrée et paramètres, mais vous ne voulez pas taper tous les paramètres à chaque fois.
    Ou peut-être configurez-vous le pipeline pour un·e collègue qui n'est pas à l'aise avec les arguments en ligne de commande.

Nous allons commencer par un aspect de la configuration qui est simplement une extension de ce avec quoi nous avons travaillé jusqu'à présent : la gestion des paramètres d'entrée.

Actuellement, notre workflow est configuré pour accepter plusieurs valeurs de paramètres via la ligne de commande, déclarées dans un bloc `params` dans le script du workflow lui-même.
L'un a une valeur par défaut définie dans sa déclaration.

Cependant, vous pourriez vouloir définir des valeurs par défaut pour tous, ou remplacer la valeur par défaut existante sans avoir à spécifier les paramètres sur la ligne de commande, ou modifier le fichier de script original.

Il existe plusieurs façons de le faire ; nous allons vous montrer trois façons basiques qui sont très couramment utilisées.

### 1.1. Configurer des valeurs dans `nextflow.config`

C'est l'approche la plus simple, bien que ce soit peut-être la moins flexible puisque le fichier `nextflow.config` principal n'est pas quelque chose que vous voulez modifier pour chaque exécution.
Mais cela a l'avantage de séparer les préoccupations de _déclarer_ les paramètres dans le workflow (ce qui y a définitivement sa place) versus fournir des _valeurs par défaut_, qui sont plus à leur place dans un fichier de configuration.

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
Pour le paramètre `batch` qui avait déjà une valeur par défaut déclarée, la syntaxe est un peu différente.
Dans le fichier du workflow, c'est une déclaration typée.
Dans la configuration, ce sont des affectations de valeurs.

Techniquement, cela suffit pour remplacer les valeurs par défaut toujours spécifiées dans le fichier du workflow.
Vous pourriez modifier la valeur par défaut pour `batch` et exécuter le workflow pour vous assurer que la valeur définie dans le fichier de configuration remplace celle définie dans le fichier du workflow.

Mais dans l'esprit de déplacer complètement la configuration vers le fichier de configuration, supprimons entièrement cette valeur par défaut du fichier du workflow.

#### 1.1.2. Supprimer la valeur par défaut pour `batch` dans le fichier du workflow

Effectuez la modification de code suivante dans le fichier du workflow `3-main.nf` :

=== "Après"

    ```groovy title="3-main.nf" linenums="9" hl_lines="6"
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

    ```groovy title="3-main.nf" linenums="9" hl_lines="6"
    /*
    * Paramètres du pipeline
    */
    params {
        input: Path
        batch: String = 'batch'
        character: String
    }
    ```

Maintenant, le fichier du workflow lui-même ne définit aucune valeur par défaut pour ces paramètres.

#### 1.1.3. Exécuter le pipeline

Testons que cela fonctionne correctement sans spécifier de paramètres dans la ligne de commande.

```bash
nextflow run 3-main.nf
```

??? success "Sortie de la commande"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `3-main.nf` [disturbed_einstein] DSL2 - revision: ede9037d02

    executor >  local (8)
    [f0/35723c] sayHello (2)       | 3 of 3 ✔
    [40/3efd1a] convertToUpper (3) | 3 of 3 ✔
    [17/e97d32] collectGreetings   | 1 of 1 ✔
    [98/c6b57b] cowpy              | 1 of 1 ✔
    ```

Cela produit toujours la même sortie qu'auparavant.

La sortie finale d'art ASCII est dans le répertoire `results/3-main/`, sous le nom `cowpy-COLLECTED-batch-output.txt`, comme avant.

??? abstract "Contenu du fichier"

    ```console title="results/3-main/cowpy-COLLECTED-batch-output.txt"
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

??? example "Scénario"

    Vous voulez expérimenter avec différents paramètres sans modifier votre fichier de configuration principal.

Vous pouvez le faire en créant un nouveau fichier `nextflow.config` dans un sous-répertoire que vous utiliserez comme répertoire de travail pour vos expériences.

#### 1.2.1. Créer le répertoire de travail avec une configuration vierge

Commençons par créer un nouveau répertoire et nous y déplacer :

```bash
mkdir -p tux-run
cd tux-run
```

Ensuite, créez un fichier de configuration vierge dans ce répertoire :

```bash
touch nextflow.config
```

Cela produit un fichier vide.

#### 1.2.2. Configurer la configuration expérimentale

Ouvrez maintenant le nouveau fichier et ajoutez les paramètres que vous souhaitez personnaliser :

```groovy title="tux-run/nextflow.config" linenums="1"
params {
    input = '../data/greetings.csv'
    batch = 'experiment'
    character = 'tux'
}
```

Notez que le chemin vers le fichier d'entrée doit refléter la structure des répertoires.

#### 1.2.3. Exécuter le pipeline

Nous pouvons maintenant exécuter notre pipeline depuis notre nouveau répertoire de travail.
Assurez-vous d'adapter le chemin en conséquence !

```bash
nextflow run ../3-main.nf
```

??? success "Sortie de la commande"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `../3-main.nf` [trusting_escher] DSL2 - revision: 356df0818d

    executor >  local (8)
    [59/b66913] sayHello (2)       [100%] 3 of 3 ✔
    [ad/f06364] convertToUpper (3) [100%] 3 of 3 ✔
    [10/714895] collectGreetings   [100%] 1 of 1 ✔
    [88/3ece98] cowpy              [100%] 1 of 1 ✔
    ```

Cela créera un nouvel ensemble de répertoires sous `tux-run/` incluant `tux-run/work/` et `tux-run/results/`.

Dans cette exécution, Nextflow combine le `nextflow.config` dans notre répertoire actuel avec le `nextflow.config` dans le répertoire racine du pipeline, et remplace ainsi le personnage par défaut (turkey) par le personnage tux.

Le fichier de sortie final devrait contenir le personnage tux disant les salutations.

??? abstract "Contenu du fichier"

    ```console title="tux-run/results/3-main/cowpy-COLLECTED-experiment-output.txt"
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

Maintenant regardons une autre façon utile de définir des valeurs de paramètres.

### 1.3. Utiliser un fichier de paramètres

??? example "Scénario"

    Vous devez partager les paramètres exacts d'exécution avec un·e collaborateur·trice, ou les enregistrer pour une publication.

L'approche du sous-répertoire fonctionne très bien pour expérimenter, mais elle implique un peu de configuration et nécessite que vous adaptiez les chemins en conséquence.
Il existe une approche plus simple pour quand vous voulez exécuter votre pipeline avec un ensemble spécifique de valeurs, ou permettre à quelqu'un d'autre de le faire avec un effort minimal.

Nextflow nous permet de spécifier des paramètres via un [fichier de paramètres](https://nextflow.io/docs/latest/config.html#parameter-file) au format YAML ou JSON, ce qui rend très pratique la gestion et la distribution d'ensembles alternatifs de valeurs par défaut, par exemple, ainsi que de valeurs de paramètres spécifiques à l'exécution.

#### 1.3.1. Examiner le fichier de paramètres d'exemple

Pour démontrer cela, nous fournissons un fichier de paramètres d'exemple dans le répertoire actuel, appelé `test-params.yaml` :

```yaml title="test-params.yaml" linenums="1"
input: "data/greetings.csv"
batch: "yaml"
character: "stegosaurus"
```

Ce fichier de paramètres contient une paire clé-valeur pour chacune des entrées que nous voulons spécifier.
Notez l'utilisation de deux-points (`:`) au lieu de signes égal (`=`) si vous comparez la syntaxe au fichier de configuration.
Le fichier de config est écrit en Groovy, tandis que le fichier de paramètres est écrit en YAML.

!!! info "Information"

    Nous fournissons également une version JSON du fichier de paramètres comme exemple mais nous n'allons pas l'exécuter ici.
    N'hésitez pas à essayer celui-là par vous-même.

#### 1.3.2. Exécuter le pipeline

Pour exécuter le workflow avec ce fichier de paramètres, ajoutez simplement `-params-file <filename>` à la commande de base.

```bash
nextflow run 3-main.nf -params-file test-params.yaml
```

??? success "Sortie de la commande"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `3-main.nf` [disturbed_sammet] DSL2 - revision: ede9037d02

    executor >  local (8)
    [2b/9a7d1e] sayHello (2)       | 3 of 3 ✔
    [5c/8f3b2a] convertToUpper (3) | 3 of 3 ✔
    [a3/29d8fb] collectGreetings   | 1 of 1 ✔
    [b7/83ef12] cowpy              | 1 of 1 ✔
    ```

Le fichier de sortie final devrait contenir le personnage stegosaurus disant les salutations.

??? abstract "Contenu du fichier"

    ```console title="results/3-main/cowpy-COLLECTED-yaml-output.txt"
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

Utiliser un fichier de paramètres peut sembler excessif quand vous n'avez que quelques paramètres à spécifier, mais certains pipelines attendent des dizaines de paramètres.
Dans ces cas, utiliser un fichier de paramètres nous permettra de fournir des valeurs de paramètres au moment de l'exécution sans avoir à taper des lignes de commande massives et sans modifier le script du workflow.

Cela facilite également la distribution d'ensembles de paramètres à des collaborateur·trices, ou comme information supplémentaire pour une publication, par exemple.
Cela rend votre travail plus reproductible par d'autres.

### Récapitulatif

Vous savez comment tirer parti des options de configuration clés pour gérer les entrées du workflow.

### Et ensuite ?

Apprenez à gérer où et comment les sorties de votre workflow sont publiées.

---

## 2. Gérer les sorties du workflow

??? example "Scénario"

    Votre pipeline publie les sorties dans un répertoire codé en dur, mais vous voulez organiser les résultats par projet ou nom d'expérience sans modifier le code du workflow à chaque fois.

Le workflow que nous avons hérité utilise des chemins pour les déclarations de sortie au niveau du workflow, ce qui n'est pas terriblement flexible et implique beaucoup de répétition.

Regardons quelques façons courantes de configurer cela pour être plus flexible.

### 2.1. Personnaliser le nom du répertoire `outputDir`

Chaque version du workflow que nous avons exécutée jusqu'à présent a publié ses sorties dans un sous-répertoire différent codé en dur dans les définitions de sortie.

Nous avons changé l'emplacement de ce sous-répertoire dans la Partie 1 en utilisant le drapeau CLI `-output-dir`, mais ce n'est toujours qu'une chaîne statique.
Configurons plutôt cela dans un fichier de config, où nous pouvons définir des chemins dynamiques plus complexes.
Nous pourrions créer un tout nouveau paramètre pour cela, mais utilisons le paramètre `batch` puisqu'il est juste là.

#### 2.1.1. Définir une valeur pour `outputDir` dans le fichier de configuration

Le chemin que Nextflow utilise pour publier les sorties est contrôlé par l'option `outputDir`.
Pour changer le chemin pour toutes les sorties, vous pouvez définir une valeur pour cette option dans le fichier de configuration `nextflow.config`.

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
    outputDir = "results_config/${params.batch}"
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

Cela remplacera le chemin par défaut intégré, `results/`, par `results_config/` plus la valeur du paramètre `batch` comme sous-répertoire.

N'oubliez pas que vous pouvez également définir cette option depuis la ligne de commande en utilisant le paramètre `-output-dir` dans votre commande (`-o` en abrégé), mais alors vous ne pourriez pas utiliser la valeur du paramètre `batch`.
L'utilisation du drapeau CLI écrasera `outputDir` dans la config s'il est défini.

#### 2.1.2. Supprimer la partie répétée du chemin codé en dur

Nous avons encore un sous-répertoire codé en dur dans les options de sortie, alors débarrassons-nous de cela maintenant.

Effectuez les modifications de code suivantes dans le fichier du workflow :

=== "Après"

    ```groovy title="3-main.nf" linenums="42" hl_lines="3 7 11 15 19"
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

    ```groovy title="3-main.nf" linenums="42" hl_lines="3 7 11 15 19"
    output {
        first_output {
            path '3-main/intermediates'
            mode 'copy'
        }
        uppercased {
            path '3-main/intermediates'
            mode 'copy'
        }
        collected {
            path '3-main/intermediates'
            mode 'copy'
        }
        batch_report {
            path '3-main'
            mode 'copy'
        }
        cowpy_art {
            path '3-main'
            mode 'copy'
        }
    }
    ```

Nous aurions également pu simplement ajouter `${params.batch}` à chaque chemin au lieu de modifier la valeur par défaut de `outputDir`, mais c'est plus concis.

#### 2.1.3. Exécuter le pipeline

Testons que cela fonctionne correctement, en définissant le nom du lot sur `outdir` depuis la ligne de commande.

```bash
nextflow run 3-main.nf --batch outdir
```

??? success "Sortie de la commande"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `3-main.nf` [amazing_church] DSL2 - revision: 6e18cd130e

    executor >  local (8)
    [9c/6a03ea] sayHello (2)       [100%] 3 of 3 ✔
    [11/9e58a6] convertToUpper (3) [100%] 3 of 3 ✔
    [c8/1977e5] collectGreetings   [100%] 1 of 1 ✔
    [38/f01eda] cowpy              [100%] 1 of 1 ✔
    ```

Cela produit toujours la même sortie qu'auparavant, sauf que cette fois nous trouvons nos sorties sous `results_config/outdir/`.

??? abstract "Contenu du répertoire"

    ```console
    results_config/outdir
    ├── cowpy-COLLECTED-outdir-output.txt
    ├── intermediates
    │   ├── Bonjour-output.txt
    │   ├── COLLECTED-outdir-output.txt
    │   ├── Hello-output.txt
    │   ├── Holà-output.txt
    │   ├── UPPER-Bonjour-output.txt
    │   ├── UPPER-Hello-output.txt
    │   └── UPPER-Holà-output.txt
    └── outdir-report.txt
    ```

Vous pouvez combiner cette approche avec des définitions de chemins personnalisés pour construire n'importe quelle hiérarchie de répertoires que vous souhaitez.

### 2.2. Organiser les sorties par process

Une façon populaire d'organiser davantage les sorties est de le faire par process, _c.-à-d._ créer des sous-répertoires pour chaque process exécuté dans le pipeline.

#### 2.2.1. Remplacer les chemins de sortie par une référence aux noms de process

Tout ce que vous avez à faire est de référencer le nom du process comme `<process>.name` dans la déclaration du chemin de sortie.

Effectuez les modifications suivantes dans le fichier du workflow :

=== "Après"

    ```groovy title="3-main.nf" linenums="42" hl_lines="3 7 11 15 19"
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

    ```groovy title="3-main.nf" linenums="42" hl_lines="3 7 11 15 19"
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

Cela supprime les éléments codés en dur restants de la configuration du chemin de sortie.

#### 2.2.2. Exécuter le pipeline

Testons que cela fonctionne correctement, en définissant le nom du lot sur `pnames` depuis la ligne de commande.

```bash
nextflow run 3-main.nf --batch pnames
```

??? success "Sortie de la commande"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `3-main.nf` [jovial_mcclintock] DSL2 - revision: ede9037d02

    executor >  local (8)
    [4a/c2e6b8] sayHello (2)       | 3 of 3 ✔
    [6f/d4a172] convertToUpper (3) | 3 of 3 ✔
    [e8/4f19d7] collectGreetings   | 1 of 1 ✔
    [f2/a85c36] cowpy              | 1 of 1 ✔
    ```

Cela produit toujours la même sortie qu'auparavant, sauf que cette fois nous trouvons nos sorties sous `results_config/pnames/`, et elles sont groupées par process.

??? abstract "Contenu du répertoire"

    ```console
    results_config/pnames/
    ├── collectGreetings
    │   ├── COLLECTED-pnames-output.txt
    │   └── pnames-report.txt
    ├── convertToUpper
    │   ├── UPPER-Bonjour-output.txt
    │   ├── UPPER-Hello-output.txt
    │   └── UPPER-Holà-output.txt
    ├── cowpy
    │   └── cowpy-COLLECTED-pnames-output.txt
    └── sayHello
        ├── Bonjour-output.txt
        ├── Hello-output.txt
        └── Holà-output.txt
    ```

!!! note "Note"

    Notez qu'ici nous avons effacé la distinction entre `intermediates` versus les sorties finales au niveau supérieur.
    Vous pouvez mélanger ces approches et même inclure plusieurs variables, par exemple en définissant le chemin de la première sortie comme `#!groovy "${params.batch}/intermediates/${sayHello.name}"`

### 2.3. Définir le mode de publication au niveau du workflow

Enfin, dans l'esprit de réduire la quantité de code répétitif, nous pouvons remplacer les déclarations `mode` par sortie par une seule ligne dans la configuration.

#### 2.3.1. Ajouter `workflow.output.mode` au fichier de configuration

Ajoutez le code suivant au fichier `nextflow.config` :

=== "Après"

    ```groovy title="nextflow.config" linenums="2" hl_lines="5"
    /*
    * Paramètres de sortie
    */
    outputDir = "results_config/${params.batch}"
    workflow.output.mode = 'copy'
    ```

=== "Avant"

    ```groovy title="nextflow.config" linenums="12"
    /*
    * Paramètres de sortie
    */
    outputDir = "results_config/${params.batch}"
    ```

Tout comme l'option `outputDir`, donner une valeur à `workflow.output.mode` dans le fichier de configuration serait suffisant pour remplacer ce qui est défini dans le fichier du workflow, mais supprimons quand même le code inutile.

#### 2.3.2. Supprimer le mode de sortie du fichier du workflow

Effectuez les modifications suivantes dans le fichier du workflow :

=== "Après"

    ```groovy title="3-main.nf" linenums="42"
    output {
        first_output {
            path { sayHello.name }
        }
        uppercased {
            path { convertToUpper.name }
        }
        collected {
            path { collectGreetings.name }
        }
        batch_report {
            path { collectGreetings.name }
        }
        cowpy_art {
            path { cowpy.name }
        }
    }
    ```

=== "Avant"

    ```groovy title="3-main.nf" linenums="42" hl_lines="3 7 11 15 19"
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

C'est plus concis, n'est-ce pas ?

#### 2.3.3. Exécuter le pipeline

Testons que cela fonctionne correctement, en définissant le nom du lot sur `outmode` depuis la ligne de commande.

```bash
nextflow run 3-main.nf --batch outmode
```

??? success "Sortie de la commande"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `3-main.nf` [rowdy_sagan] DSL2 - revision: ede9037d02

    executor >  local (8)
    [5b/d91e3c] sayHello (2)       | 3 of 3 ✔
    [8a/f6c241] convertToUpper (3) | 3 of 3 ✔
    [89/cd3a48] collectGreetings   | 1 of 1 ✔
    [9e/71fb52] cowpy              | 1 of 1 ✔
    ```

Cela produit toujours la même sortie qu'auparavant, sauf que cette fois nous trouvons nos sorties sous `results_config/outmode/`.
Ce sont toujours toutes des copies propres, pas des liens symboliques.

??? abstract "Contenu du répertoire"

    ```console
    results_config/outmode/
    ├── collectGreetings
    │   ├── COLLECTED-outmode-output.txt
    │   └── outmode-report.txt
    ├── convertToUpper
    │   ├── UPPER-Bonjour-output.txt
    │   ├── UPPER-Hello-output.txt
    │   └── UPPER-Holà-output.txt
    ├── cowpy
    │   └── cowpy-COLLECTED-outmode-output.txt
    └── sayHello
        ├── Bonjour-output.txt
        ├── Hello-output.txt
        └── Holà-output.txt
    ```

La principale raison pour laquelle vous pourriez encore vouloir utiliser la façon par sortie de définir le mode est si vous voulez mélanger dans le même workflow, _c.-à-d._ avoir certaines sorties copiées et d'autres liées symboliquement.

Il y a beaucoup d'autres options que vous pouvez personnaliser de cette façon, mais nous espérons que cela vous donne une idée de la gamme d'options et comment les utiliser efficacement selon vos préférences.

### Récapitulatif

Vous savez comment contrôler le nommage et la structure des répertoires où vos sorties sont publiées, ainsi que le mode de publication des sorties du workflow.

### Et ensuite ?

Apprenez à adapter la configuration de votre workflow à votre environnement de calcul, en commençant par la technologie d'empaquetage logiciel.

---

## 3. Sélectionner une technologie d'empaquetage logiciel

Jusqu'à présent, nous avons examiné les éléments de configuration qui contrôlent comment les entrées entrent et d'où les sorties sortent. Maintenant il est temps de se concentrer plus spécifiquement sur l'adaptation de la configuration de votre workflow à votre environnement de calcul.

La première étape sur ce chemin est de spécifier d'où vont venir les packages logiciels qui seront exécutés à chaque étape.
Sont-ils déjà installés dans l'environnement de calcul local ?
Devons-nous récupérer des images et les exécuter via un système de conteneurs ?
Ou devons-nous récupérer des packages Conda et construire un environnement Conda local ?

Dans la toute première partie de cette formation (Parties 1-4) nous avons juste utilisé des logiciels installés localement dans notre workflow.
Ensuite dans la Partie 5, nous avons introduit les conteneurs Docker et le fichier `nextflow.config`, que nous avons utilisé pour activer l'utilisation des conteneurs Docker.

Maintenant voyons comment nous pouvons configurer une option d'empaquetage logiciel alternative via le fichier `nextflow.config`.

### 3.1. Désactiver Docker et activer Conda dans le fichier de config

??? example "Scénario"

    Vous déplacez votre pipeline vers un cluster HPC où Docker n'est pas autorisé pour des raisons de sécurité.
    Le cluster prend en charge Singularity et Conda, donc vous devez changer votre configuration en conséquence.

Comme noté précédemment, Nextflow prend en charge plusieurs technologies de conteneurs, y compris Singularity (qui est plus largement utilisé sur HPC), ainsi que des gestionnaires de packages logiciels tels que Conda.

Nous pouvons changer notre fichier de configuration pour utiliser Conda au lieu de Docker.
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

Cela permettra à Nextflow de créer et d'utiliser des environnements Conda pour les processes qui ont des packages Conda spécifiés.
Ce qui signifie que nous devons maintenant en ajouter un à notre process `cowpy` !

### 3.2. Spécifier un package Conda dans la définition du process

Nous avons déjà récupéré l'URI pour un package Conda contenant l'outil `cowpy` : `conda-forge::cowpy==1.1.5`

Maintenant nous ajoutons l'URI à la définition du process `cowpy` en utilisant la directive `conda` :

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

Essayons cela.

```bash
nextflow run 3-main.nf --batch conda
```

??? success "Sortie de la commande"

    ```console title="Sortie"
    N E X T F L O W   ~  version 25.10.2

    Launching `3-main.nf` [trusting_lovelace] DSL2 - revision: 028a841db1

    executor >  local (8)
    [ee/4ca1f2] sayHello (3)       | 3 of 3 ✔
    [20/2596a7] convertToUpper (1) | 3 of 3 ✔
    [b3/e15de5] collectGreetings   | 1 of 1 ✔
    [c5/af5f88] cowpy              | 1 of 1 ✔
    ```

Cela devrait fonctionner sans problème et produire les mêmes sorties qu'auparavant sous `results_config/conda`.

En coulisses, Nextflow a récupéré les packages Conda et créé l'environnement, ce qui prend normalement un peu de travail ; donc c'est agréable de ne pas avoir à faire tout cela nous-mêmes !

!!! info "Information"

    Cela s'exécute rapidement car le package `cowpy` est assez petit, mais si vous travaillez avec de gros packages, cela peut prendre un peu plus de temps que d'habitude la première fois, et vous pourriez voir la sortie console rester « bloquée » pendant une minute ou plus avant de se terminer.
    C'est normal et est dû au travail supplémentaire que Nextflow fait la première fois que vous utilisez un nouveau package.

De notre point de vue, cela semble fonctionner exactement de la même manière que l'exécution avec Docker, même si en arrière-plan les mécanismes sont un peu différents.

Cela signifie que nous sommes prêt·es à exécuter avec des environnements Conda si nécessaire.

??? info "Mélanger Docker et Conda"

    Puisque ces directives sont assignées par process, il est possible de « mélanger », _c.-à-d._ de configurer certains des processes de votre workflow pour s'exécuter avec Docker et d'autres avec Conda, par exemple, si l'infrastructure de calcul que vous utilisez prend en charge les deux.
    Dans ce cas, vous activeriez à la fois Docker et Conda dans votre fichier de configuration.
    Si les deux sont disponibles pour un process donné, Nextflow priorisera les conteneurs.

    Et comme noté précédemment, Nextflow prend en charge plusieurs autres technologies d'empaquetage logiciel et de conteneurs, donc vous n'êtes pas limité·e à ces deux-là.

### Récapitulatif

Vous savez comment configurer quel package logiciel chaque process doit utiliser, et comment basculer entre les technologies.

### Et ensuite ?

Apprenez à changer la plateforme d'exécution utilisée par Nextflow pour faire réellement le travail.

---

## 4. Sélectionner une plateforme d'exécution

??? example "Scénario"

    Vous avez développé et testé votre pipeline sur votre ordinateur portable, mais maintenant vous devez l'exécuter sur des milliers d'échantillons.
    Votre institution a un cluster HPC avec un ordonnanceur Slurm que vous aimeriez utiliser à la place.

Jusqu'à présent, nous avons exécuté notre pipeline avec l'executor local.
Cela exécute chaque tâche sur la machine sur laquelle Nextflow s'exécute.
Quand Nextflow démarre, il examine les CPU et la mémoire disponibles.
Si les ressources des tâches prêtes à s'exécuter dépassent les ressources disponibles, Nextflow retiendra les dernières tâches de l'exécution jusqu'à ce qu'une ou plusieurs des tâches précédentes soient terminées, libérant les ressources nécessaires.

L'executor local est pratique et efficace, mais il est limité à cette seule machine. Pour de très grandes charges de travail, vous pourriez découvrir que votre machine locale est un goulot d'étranglement, soit parce que vous avez une seule tâche qui nécessite plus de ressources que vous n'en avez disponibles, soit parce que vous avez tellement de tâches qu'attendre qu'une seule machine les exécute prendrait trop de temps.

Nextflow prend en charge [de nombreux backends d'exécution différents](https://nextflow.io/docs/latest/executor.html), y compris les ordonnanceurs HPC (Slurm, LSF, SGE, PBS, Moab, OAR, Bridge, HTCondor et d'autres) ainsi que les backends d'exécution cloud (AWS Batch, Google Cloud Batch, Azure Batch, Kubernetes et plus).

### 4.1. Cibler un backend différent

Le choix de l'executor est défini par une directive de process appelée `executor`.
Par défaut, elle est définie sur `local`, donc la configuration suivante est implicite :

```groovy title="Configuration intégrée"
process {
    executor = 'local'
}
```

Pour définir l'executor pour cibler un backend différent, vous spécifieriez simplement l'executor que vous voulez en utilisant une syntaxe similaire à celle décrite ci-dessus pour les allocations de ressources (voir [Executors](https://nextflow.io/docs/latest/executor.html) pour toutes les options).

```groovy title="nextflow.config"
process {
    executor = 'slurm'
}
```

!!! warning "Avertissement"

    Nous ne pouvons pas réellement tester cela dans l'environnement de formation car il n'est pas configuré pour se connecter à un HPC.

### 4.2. Gérer la syntaxe spécifique au backend pour les paramètres d'exécution

La plupart des plateformes de calcul haute performance permettent (et parfois exigent) que vous spécifiiez certains paramètres tels que les demandes et limitations d'allocation de ressources (par exemple le nombre de CPU et la mémoire) et le nom de la file d'attente des jobs à utiliser.

Malheureusement, chacun de ces systèmes utilise différentes technologies, syntaxes et configurations pour définir comment un job doit être défini et soumis à l'ordonnanceur concerné.

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
Il fournit une syntaxe standardisée pour que vous puissiez spécifier les propriétés pertinentes telles que `cpus`, `memory` et `queue` une seule fois (voir [Directives de process](https://nextflow.io/docs/latest/reference/process.html#process-directives) pour toutes les options disponibles).
Ensuite, au moment de l'exécution, Nextflow utilisera ces paramètres pour générer les scripts spécifiques au backend appropriés en fonction du paramètre d'executor.

Nous couvrirons cette syntaxe standardisée dans la section suivante.

### Récapitulatif

Vous savez maintenant comment changer l'executor pour utiliser différents types d'infrastructure de calcul.

### Et ensuite ?

Apprenez à évaluer et exprimer les allocations et limitations de ressources dans Nextflow.

---

## 5. Contrôler les allocations de ressources de calcul

??? example "Scénario"

    Votre pipeline continue d'échouer sur le cluster parce que les tâches sont tuées pour avoir dépassé les limites de mémoire.
    Ou peut-être êtes-vous facturé·e pour des ressources que vous n'utilisez pas et souhaitez-vous optimiser les coûts.

La plupart des plateformes de calcul haute performance permettent (et parfois exigent) que vous spécifiiez certains paramètres d'allocation de ressources tels que le nombre de CPU et la mémoire.

Par défaut, Nextflow utilisera un seul CPU et 2 Go de mémoire pour chaque process.
Les directives de process correspondantes sont appelées `cpus` et `memory`, donc la configuration suivante est implicite :

```groovy title="Configuration intégrée" linenums="1"
process {
    cpus = 1
    memory = 2.GB
}
```

Vous pouvez modifier ces valeurs, soit pour tous les processes, soit pour des processes nommés spécifiques, en utilisant des directives de process supplémentaires dans votre fichier de configuration.
Nextflow les traduira en instructions appropriées pour l'executor choisi.

Mais comment savez-vous quelles valeurs utiliser ?

### 5.1. Exécuter le workflow pour générer un rapport d'utilisation des ressources

??? example "Scénario"

    Vous ne savez pas combien de mémoire ou de CPU vos processes ont besoin et voulez éviter de gaspiller des ressources ou d'avoir des jobs tués.

Si vous ne savez pas à l'avance combien de CPU et de mémoire vos processes sont susceptibles d'avoir besoin, vous pouvez faire du profilage de ressources, ce qui signifie que vous exécutez le workflow avec quelques allocations par défaut, enregistrez combien chaque process a utilisé, et à partir de là, estimez comment ajuster les allocations de base.

Commodément, Nextflow inclut des outils intégrés pour faire cela, et générera volontiers un rapport pour vous sur demande.

Pour ce faire, ajoutez `-with-report <filename>.html` à votre ligne de commande.

```bash
nextflow run 3-main.nf -with-report report-config-1.html
```

Le rapport est un fichier html, que vous pouvez télécharger et ouvrir dans votre navigateur. Vous pouvez également faire un clic droit dessus dans l'explorateur de fichiers à gauche et cliquer sur `Show preview` pour le voir dans l'environnement de formation.

Prenez quelques minutes pour parcourir le rapport et voir si vous pouvez identifier des opportunités d'ajustement des ressources.
Assurez-vous de cliquer sur les onglets qui montrent les résultats d'utilisation en pourcentage de ce qui a été alloué.

Consultez [Rapports](https://nextflow.io/docs/latest/reports.html) pour la documentation sur toutes les fonctionnalités disponibles.

### 5.2. Définir les allocations de ressources pour tous les processes

Le profilage montre que les processes de notre workflow de formation sont très légers, alors réduisons l'allocation de mémoire par défaut à 1 Go par process.

Ajoutez ce qui suit à votre fichier `nextflow.config`, avant la section des paramètres du pipeline :

```groovy title="nextflow.config" linenums="4"
/*
* Paramètres de process
*/
process {
    memory = 1.GB
}
```

Cela aidera à réduire la quantité de calcul que nous consommons.

### 5.3. Définir les allocations de ressources pour un process spécifique

En même temps, nous allons prétendre que le process `cowpy` nécessite plus de ressources que les autres, juste pour que nous puissions démontrer comment ajuster les allocations pour un process individuel.

=== "Après"

    ```groovy title="nextflow.config" linenums="4" hl_lines="6-9"
    /*
    * Paramètres de process
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
    * Paramètres de process
    */
    process {
        memory = 1.GB
    }
    ```

Avec cette configuration, tous les processes demanderont 1 Go de mémoire et un seul CPU (la valeur par défaut implicite), sauf le process `cowpy`, qui demandera 2 Go et 2 CPU.

!!! info "Information"

    Si vous avez une machine avec peu de CPU et que vous en allouez un grand nombre par process, vous pourriez voir des appels de process mis en file d'attente les uns derrière les autres.
    C'est parce que Nextflow s'assure que nous ne demandons pas plus de CPU qu'il n'y en a de disponibles.

### 5.4. Exécuter le workflow avec la configuration mise à jour

Essayons cela, en fournissant un nom de fichier différent pour le rapport de profilage afin de pouvoir comparer les performances avant et après les modifications de configuration.

```bash
nextflow run 3-main.nf -with-report report-config-2.html
```

Vous ne remarquerez probablement pas de différence réelle puisque c'est une charge de travail si petite, mais c'est l'approche que vous utiliseriez pour analyser les performances et les besoins en ressources d'un workflow du monde réel.

C'est très utile lorsque vos processes ont des besoins en ressources différents. Cela vous permet de dimensionner correctement les allocations de ressources que vous définissez pour chaque process en fonction de données réelles, pas de suppositions.

!!! tip "Astuce"

    Ceci n'est qu'un petit aperçu de ce que vous pouvez faire pour optimiser votre utilisation des ressources.
    Nextflow lui-même a une [logique de réessai dynamique](https://nextflow.io/docs/latest/process.html#dynamic-task-resources) vraiment intéressante intégrée pour réessayer les jobs qui échouent en raison de limitations de ressources.
    De plus, la plateforme Seqera offre des outils pilotés par l'IA pour optimiser vos allocations de ressources automatiquement également.

### 5.5. Ajouter des limites de ressources

Selon l'executor de calcul et l'infrastructure de calcul que vous utilisez, il peut y avoir des contraintes sur ce que vous pouvez (ou devez) allouer.
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

Nextflow traduira ces valeurs en instructions appropriées selon l'executor que vous avez spécifié.

Nous n'allons pas exécuter cela, puisque nous n'avons pas accès à l'infrastructure pertinente dans l'environnement de formation.
Cependant, si vous deviez essayer d'exécuter le workflow avec des allocations de ressources qui dépassent ces limites, puis rechercher la commande `sbatch` dans le fichier de script `.command.run`, vous verriez que les demandes qui sont réellement envoyées à l'executor sont plafonnées aux valeurs spécifiées par `resourceLimits`.

??? info "Configurations de référence institutionnelles"

    Le projet nf-core a compilé une [collection de fichiers de configuration](https://nf-co.re/configs/) partagés par diverses institutions à travers le monde, couvrant une large gamme d'executors HPC et cloud.

    Ces configurations partagées sont précieuses à la fois pour les personnes qui y travaillent et peuvent donc simplement utiliser la configuration de leur institution immédiatement, et comme modèle pour les personnes qui cherchent à développer une configuration pour leur propre infrastructure.

### Récapitulatif

Vous savez comment générer un rapport de profilage pour évaluer l'utilisation des ressources et comment modifier les allocations de ressources pour tous les processes et/ou pour des processes individuels, ainsi que définir des limitations de ressources pour l'exécution sur HPC.

### Et ensuite ?

Apprenez à configurer des profils de configuration prédéfinis et à basculer entre eux au moment de l'exécution.

---

## 6. Utiliser des profils pour basculer entre des configurations prédéfinies

??? example "Scénario"

    Vous basculez régulièrement entre l'exécution de pipelines sur votre ordinateur portable pour le développement et sur le HPC de votre institution pour les exécutions de production.
    Vous en avez assez de changer manuellement les paramètres de configuration chaque fois que vous changez d'environnement.

Nous vous avons montré un certain nombre de façons de personnaliser la configuration de votre pipeline selon le projet sur lequel vous travaillez ou l'environnement de calcul que vous utilisez.

Vous pourriez vouloir basculer entre des paramètres alternatifs selon l'infrastructure de calcul que vous utilisez. Par exemple, vous pourriez vouloir développer et exécuter des tests à petite échelle localement sur votre ordinateur portable, puis exécuter des charges de travail à grande échelle sur HPC ou cloud.

Nextflow vous permet de configurer n'importe quel nombre de [**profils**](https://nextflow.io/docs/latest/config.html#profiles) qui décrivent différentes configurations, que vous pouvez ensuite sélectionner au moment de l'exécution en utilisant un argument de ligne de commande, plutôt que d'avoir à modifier le fichier de configuration lui-même.

### 6.1. Créer des profils pour basculer entre le développement local et l'exécution sur HPC

Configurons deux profils alternatifs ; un pour exécuter des charges à petite échelle sur un ordinateur ordinaire, où nous utiliserons des conteneurs Docker, et un pour exécuter sur un HPC universitaire avec un ordonnanceur Slurm, où nous utiliserons des packages Conda.

#### 6.1.1. Configurer les profils

Ajoutez ce qui suit à votre fichier `nextflow.config`, après la section des paramètres du pipeline mais avant les paramètres de sortie :

```groovy title="nextflow.config" linenums="24"
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
```

Vous voyez que pour le HPC universitaire, nous spécifions également des limitations de ressources.

#### 6.1.2. Exécuter le workflow avec un profil

Pour spécifier un profil dans notre ligne de commande Nextflow, nous utilisons l'argument `-profile`.

Essayons d'exécuter le workflow avec la configuration `my_laptop`.

```bash
nextflow run 3-main.nf -profile my_laptop
```

??? success "Sortie de la commande"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `3-main.nf` [gigantic_brazil] DSL2 - revision: ede9037d02

    executor >  local (8)
    [58/da9437] sayHello (3)       | 3 of 3 ✔
    [35/9cbe77] convertToUpper (2) | 3 of 3 ✔
    [67/857d05] collectGreetings   | 1 of 1 ✔
    [37/7b51b5] cowpy              | 1 of 1 ✔
    ```

Comme vous pouvez le voir, cela nous permet de basculer entre les configurations très commodément au moment de l'exécution.

!!! warning "Avertissement"

    Le profil `univ_hpc` ne fonctionnera pas correctement dans l'environnement de formation puisque nous n'avons pas accès à un ordonnanceur Slurm.

Si à l'avenir nous trouvons d'autres éléments de configuration qui sont toujours co-occurrents avec ceux-ci, nous pouvons simplement les ajouter au(x) profil(s) correspondant(s).
Nous pouvons également créer des profils supplémentaires s'il y a d'autres éléments de configuration que nous voulons regrouper.

### 6.2. Créer un profil de paramètres de test

??? example "Scénario"

    Vous voulez que d'autres puissent essayer rapidement votre pipeline sans avoir à rassembler leurs propres données d'entrée.

Les profils ne sont pas seulement pour la configuration d'infrastructure.
Nous pouvons également les utiliser pour définir des valeurs par défaut pour les paramètres du workflow, pour faciliter aux autres l'essai du workflow sans avoir à rassembler eux-mêmes des valeurs d'entrée appropriées.
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

```groovy title="nextflow.config" linenums="24"
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

Commodément, les profils ne sont pas mutuellement exclusifs, donc nous pouvons spécifier plusieurs profils dans notre ligne de commande en utilisant la syntaxe suivante `-profile <profile1>,<profile2>` (pour n'importe quel nombre de profils).

Si vous combinez des profils qui définissent des valeurs pour les mêmes éléments de configuration et sont décrits dans le même fichier de configuration, Nextflow résoudra le conflit en utilisant la valeur qu'il a lue en dernier (_c.-à-d._ ce qui vient plus tard dans le fichier).
Si les paramètres en conflit sont définis dans différentes sources de configuration, l'[ordre de priorité](https://www.nextflow.io/docs/latest/config.html#configuration-file) par défaut s'applique.

Essayons d'ajouter le profil de test à notre commande précédente :

```bash
nextflow run 3-main.nf -profile my_laptop,test
```

??? success "Sortie de la commande"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `3-main.nf` [jovial_coulomb] DSL2 - revision: 46a6763141

    executor >  local (8)
    [9b/687cdc] sayHello (2)       | 3 of 3 ✔
    [ca/552187] convertToUpper (3) | 3 of 3 ✔
    [e8/83e306] collectGreetings   | 1 of 1 ✔
    [fd/e84fa9] cowpy              | 1 of 1 ✔
    ```

Cela utilisera Docker où possible et produira des sorties sous `results_config/test`, et cette fois le personnage est le duo comique `dragonandcow`.

??? abstract "Contenu du fichier"

    ```console title="results_config/test/"
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

    Nous pouvons pointer vers des URL pour des fichiers plus grands qui sont stockés en externe.
    Nextflow les téléchargera automatiquement tant qu'il y a une connexion ouverte.

    Pour plus de détails, consultez la Quête secondaire [Travailler avec les fichiers](../side_quests/working_with_files.md)

### 6.3. Utiliser `nextflow config` pour voir la configuration résolue

Comme noté ci-dessus, parfois le même paramètre peut être défini sur différentes valeurs dans des profils que vous voulez combiner.
Et plus généralement, il y a de nombreux endroits où des éléments de configuration peuvent être stockés, et parfois les mêmes propriétés peuvent être définies sur différentes valeurs à différents endroits.

Nextflow applique un [ordre de priorité](https://nextflow.io/docs/latest/config.html#configuration-file) défini pour résoudre les conflits, mais cela peut être difficile à déterminer vous-même.
Et même si rien n'est en conflit, il peut être fastidieux de rechercher tous les endroits possibles où les choses pourraient être configurées.

Heureusement, Nextflow inclut un outil utilitaire pratique appelé `config` qui peut automatiser tout ce processus pour vous.

L'outil `config` explorera tout le contenu de votre répertoire de travail actuel, aspirera tous les fichiers de configuration, et produira la configuration entièrement résolue que Nextflow utiliserait pour exécuter le workflow.
Cela vous permet de découvrir quels paramètres seront utilisés sans avoir à lancer quoi que ce soit.

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

    outputDir = 'results_config/batch'

    workflow {
      output {
          mode = 'copy'
      }
    }
    ```

Cela vous montre la configuration de base que vous obtenez si vous ne spécifiez rien de supplémentaire dans la ligne de commande.

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

    outputDir = 'results_config/test'

    workflow {
      output {
          mode = 'copy'
      }
    }
    ```

Cela devient particulièrement utile pour les projets complexes qui impliquent plusieurs couches de configuration.

### Récapitulatif

Vous savez comment utiliser les profils pour sélectionner une configuration prédéfinie au moment de l'exécution avec un minimum de tracas.
Plus généralement, vous savez comment configurer vos exécutions de workflow pour convenir à différentes plateformes de calcul et améliorer la reproductibilité de vos analyses.

### Et ensuite ?

Apprenez à exécuter des pipelines directement depuis des dépôts distants comme GitHub.

---

## 7. Exécuter des pipelines depuis des dépôts distants

??? example "Scénario"

    Vous voulez exécuter un pipeline bien établi comme ceux de nf-core sans avoir à télécharger et gérer le code vous-même.

Jusqu'à présent, nous avons exécuté des scripts de workflow situés dans le répertoire actuel.
En pratique, vous voudrez souvent exécuter des pipelines stockés dans des dépôts distants, tels que GitHub.

Nextflow rend cela simple : vous pouvez exécuter n'importe quel pipeline directement depuis une URL de dépôt Git sans avoir à le télécharger manuellement d'abord.

### 7.1. Exécuter un pipeline depuis GitHub

La syntaxe de base pour exécuter un pipeline distant est `nextflow run <repository>`, où `<repository>` peut être un chemin de dépôt GitHub comme `nextflow-io/hello`, une URL complète, ou un chemin vers GitLab, Bitbucket, ou d'autres services d'hébergement Git.

Essayez d'exécuter le pipeline de démonstration officiel Nextflow « hello » :

```bash
nextflow run nextflow-io/hello
```

??? success "Sortie de la commande"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Pulling nextflow-io/hello ...
     downloaded from https://github.com/nextflow-io/hello.git
    Launching `https://github.com/nextflow-io/hello` [sleepy_swanson] DSL2 - revision: 2ce0b0e294 [master]

    executor >  local (4)
    [ba/08236d] sayHello (4) [100%] 4 of 4 ✔
    Ciao world!

    Hello world!

    Bonjour world!

    Hola world!
    ```

La première fois que vous exécutez un pipeline distant, Nextflow le télécharge et le met en cache localement.
Les exécutions suivantes utilisent la version mise en cache à moins que vous ne demandiez explicitement une mise à jour.

### 7.2. Spécifier une version pour la reproductibilité

Par défaut, Nextflow exécute la dernière version de la branche par défaut.
Vous pouvez spécifier une version particulière (tag), une branche ou un commit en utilisant le drapeau `-r` :

```bash
nextflow run nextflow-io/hello -r v1.3
```

??? success "Sortie de la commande"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `https://github.com/nextflow-io/hello` [sick_carson] DSL2 - revision: 2ce0b0e294 [v1.3]

    executor >  local (4)
    [61/e11f77] sayHello (4) [100%] 4 of 4 ✔
    Ciao world!

    Bonjour world!

    Hello world!

    Hola world!
    ```

Spécifier des versions exactes est essentiel pour la reproductibilité.

### Récapitulatif

Vous savez comment exécuter des pipelines directement depuis GitHub et d'autres dépôts distants, et comment spécifier des versions pour la reproductibilité.

### Et ensuite ?

Félicitez-vous bien !
Vous savez tout ce que vous devez savoir pour commencer à exécuter et gérer des pipelines Nextflow.

Cela conclut cette formation, mais si vous êtes impatient·e de continuer à apprendre, nous avons deux recommandations principales :

- Si vous voulez approfondir le développement de vos propres pipelines, consultez [Hello Nextflow](../hello_nextflow/index.md), une formation pour débutants qui couvre la même progression générale que celle-ci mais va beaucoup plus en détail sur les channels et les opérateurs.
- Si vous souhaitez continuer à apprendre comment exécuter des pipelines Nextflow sans aller plus profondément dans le code, consultez la première partie de [Hello nf-core](../hello_nf-core/index.md), qui introduit les outils pour trouver et exécuter des pipelines du projet très populaire [nf-core](https://nf-co.re/).

Amusez-vous bien !

---

## Quiz

<quiz>
Lorsque des valeurs de paramètres sont définies à la fois dans le fichier du workflow et dans `nextflow.config`, laquelle prévaut ?
- [ ] La valeur du fichier du workflow
- [x] La valeur du fichier de configuration
- [ ] La première valeur rencontrée
- [ ] Cela provoque une erreur

En savoir plus : [1.1. Configurer des valeurs dans `nextflow.config`](#11-configurer-des-valeurs-dans-nextflowconfig)
</quiz>

<quiz>
Quelle est la différence de syntaxe entre définir une valeur par défaut de paramètre dans un fichier de workflow vs. un fichier de config ?
- [ ] Ils utilisent la même syntaxe
- [x] Le workflow utilise une déclaration typée (`#!groovy param: Type = value`), la config utilise une affectation (`#!groovy param = value`)
- [ ] La config utilise une déclaration typée, le workflow utilise une affectation
- [ ] Seuls les fichiers de config peuvent définir des valeurs par défaut

En savoir plus : [1.1. Configurer des valeurs dans `nextflow.config`](#11-configurer-des-valeurs-dans-nextflowconfig)
</quiz>

<quiz>
Comment spécifiez-vous un fichier de paramètres lors de l'exécution d'un workflow ?
- [ ] `--params params.yaml`
- [ ] `-config params.yaml`
- [x] `-params-file params.yaml`
- [ ] `--input-params params.yaml`

En savoir plus : [1.3. Utiliser un fichier de paramètres](#13-utiliser-un-fichier-de-parametres)
</quiz>

<quiz>
Que contrôle l'option de configuration `outputDir` ?
- [ ] L'emplacement du répertoire de travail
- [x] Le chemin de base où les sorties du workflow sont publiées
- [ ] Le répertoire pour les fichiers de log
- [ ] L'emplacement des fichiers de module

En savoir plus : [2.1. Personnaliser le nom du répertoire outputDir](#21-personnaliser-le-nom-du-repertoire-outputdir)
</quiz>

<quiz>
Comment référencez-vous dynamiquement un nom de process dans la configuration du chemin de sortie ?
- [ ] `#!groovy ${processName}`
- [ ] `#!groovy path "<process>.name"`
- [x] `#!groovy path { <process>.name }`
- [ ] `@processName`

En savoir plus : [2.2. Organiser les sorties par process](#22-organiser-les-sorties-par-process)
</quiz>

<quiz>
Si Docker et Conda sont tous deux activés et qu'un process a les deux directives, lequel est priorisé ?
- [x] Docker (conteneurs)
- [ ] Conda
- [ ] Le premier défini dans le process
- [ ] Cela provoque une erreur

En savoir plus : [3. Sélectionner une technologie d'empaquetage logiciel](#3-selectionner-une-technologie-dempaquetage-logiciel)
</quiz>

<quiz>
Quel est l'executor par défaut dans Nextflow ?
- [x] `local`
- [ ] `slurm`
- [ ] `kubernetes`
- [ ] `aws`

En savoir plus : [4. Sélectionner une plateforme d'exécution](#4-selectionner-une-plateforme-dexecution)
</quiz>

<quiz>
Quelle commande génère un rapport d'utilisation des ressources ?
- [ ] `nextflow run workflow.nf -with-metrics`
- [ ] `nextflow run workflow.nf -with-stats`
- [x] `nextflow run workflow.nf -with-report report.html`
- [ ] `nextflow run workflow.nf -profile report`

En savoir plus : [5.1. Exécuter le workflow pour générer un rapport d'utilisation des ressources](#51-executer-le-workflow-pour-generer-un-rapport-dutilisation-des-ressources)
</quiz>

<quiz>
Comment définissez-vous les exigences de ressources pour un process spécifique nommé `cowpy` dans le fichier de config ?
- [ ] `#!groovy cowpy.memory = '2.GB'`
- [ ] `#!groovy process.cowpy.memory = '2.GB'`
- [x] `#!groovy process { withName: 'cowpy' { memory = '2.GB' } }`
- [ ] `#!groovy resources.cowpy.memory = '2.GB'`

En savoir plus : [5.3. Définir les allocations de ressources pour un process spécifique](#53-definir-les-allocations-de-ressources-pour-un-process-specifique)
</quiz>

<quiz>
Que fait la directive `resourceLimits` ?
- [ ] Définit les exigences minimales de ressources
- [ ] Alloue des ressources aux processes
- [x] Plafonne les ressources maximales qui peuvent être demandées
- [ ] Surveille l'utilisation des ressources en temps réel

En savoir plus : [5.5. Ajouter des limites de ressources](#55-ajouter-des-limites-de-ressources)
</quiz>

<quiz>
Comment spécifiez-vous plusieurs profils dans une seule commande ?
- [ ] `-profile profile1 -profile profile2`
- [ ] `-profiles profile1,profile2`
- [x] `-profile profile1,profile2`
- [ ] `--profile profile1 --profile profile2`

En savoir plus : [6. Utiliser des profils pour basculer entre des configurations prédéfinies](#6-utiliser-des-profils-pour-basculer-entre-des-configurations-predefinies)
</quiz>

<quiz>
Quelle commande affiche la configuration entièrement résolue que Nextflow utiliserait ?
- [ ] `nextflow show-config`
- [ ] `nextflow settings`
- [x] `nextflow config`
- [ ] `nextflow resolve`

En savoir plus : [6.3. Utiliser `nextflow config` pour voir la configuration résolue](#63-utiliser-nextflow-config-pour-voir-la-configuration-resolue)
</quiz>

<quiz>
À quoi peuvent servir les profils ? (Sélectionnez toutes les réponses applicables)
- [x] Définir des paramètres spécifiques à l'infrastructure (executors, conteneurs)
- [x] Définir des limites de ressources pour différents environnements
- [x] Fournir des paramètres de test pour faciliter le test du workflow
- [ ] Définir de nouveaux processes

En savoir plus : [6. Utiliser des profils pour basculer entre des configurations prédéfinies](#6-utiliser-des-profils-pour-basculer-entre-des-configurations-predefinies)
</quiz>
