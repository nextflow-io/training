# Partie 3 : Utiliser un module nf-core

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Traduction assistée par IA - [en savoir plus et suggérer des améliorations](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

Dans cette troisième partie du cours de formation Hello nf-core, nous vous montrons comment trouver, installer et utiliser un module nf-core existant dans votre pipeline.

L'un des grands avantages de travailler avec nf-core est la possibilité de tirer parti de modules pré-construits et testés du dépôt [nf-core/modules](https://github.com/nf-core/modules).
Plutôt que d'écrire chaque processus à partir de zéro, vous pouvez installer et utiliser des modules maintenus par la communauté qui suivent les meilleures pratiques.

Pour démontrer comment cela fonctionne, nous remplacerons le module personnalisé `collectGreetings` par le module `cat/cat` de nf-core/modules dans le pipeline `core-hello`.

??? info "Comment commencer à partir de cette section"

    Cette section du cours suppose que vous avez terminé la [Partie 2 : Réécrire Hello pour nf-core](./02_rewrite_hello.md) et que vous disposez d'un pipeline `core-hello` fonctionnel.

    Si vous n'avez pas terminé la Partie 2 ou si vous souhaitez repartir de zéro pour cette partie, vous pouvez utiliser la solution `core-hello-part2` comme point de départ.
    Exécutez cette commande depuis le répertoire `hello-nf-core/` :

    ```bash
    cp -r solutions/core-hello-part2 core-hello
    cd core-hello
    ```

    Cela vous donne un pipeline nf-core entièrement fonctionnel, prêt pour l'ajout de modules.
    Vous pouvez tester qu'il s'exécute avec succès en exécutant la commande suivante :

    ```bash
    nextflow run . --outdir core-hello-results -profile test,docker --validate_params false
    ```

---

## 1. Trouver et installer un module nf-core approprié

Tout d'abord, apprenons comment trouver un module nf-core existant et l'installer dans notre pipeline.

Nous allons chercher à remplacer le processus `collectGreetings`, qui utilise la commande Unix `cat` pour concaténer plusieurs fichiers de salutations en un seul.
La concaténation de fichiers est une opération très courante, il est donc raisonnable de penser qu'il pourrait déjà exister un module dans nf-core conçu à cet effet.

Plongeons-nous dans le sujet.

### 1.1. Parcourir les modules disponibles sur le site web nf-core

Le projet nf-core maintient un catalogue centralisé de modules sur [https://nf-co.re/modules](https://nf-co.re/modules).

Accédez à la page des modules dans votre navigateur web et utilisez la barre de recherche pour rechercher 'concatenate'.

![résultats de recherche de module](./img/module-search-results.png)

Comme vous pouvez le voir, il y a pas mal de résultats, dont beaucoup sont des modules conçus pour concaténer des types de fichiers très spécifiques.
Parmi eux, vous devriez voir un module appelé `cat_cat` qui est générique.

!!! note "Convention de nommage des modules"

    Le caractère de soulignement (`_`) est utilisé comme substitut du caractère barre oblique (`/`) dans les noms de modules.

    Les modules nf-core suivent la convention de nommage `logiciel/commande` lorsqu'un outil fournit plusieurs commandes, comme `samtools/view` (package samtools, commande view) ou `gatk/haplotypecaller` (package GATK, commande HaplotypeCaller).
    Pour les outils qui ne fournissent qu'une seule commande principale, les modules utilisent un seul niveau comme `fastqc` ou `multiqc`.

Cliquez sur la boîte du module `cat_cat` pour voir la documentation du module.

La page du module affiche :

- Une brève description : "A module for concatenation of gzipped or uncompressed files"
- La commande d'installation : `nf-core modules install cat/cat`
- La structure des canaux d'entrée et de sortie
- Les paramètres disponibles

### 1.2. Lister les modules disponibles depuis la ligne de commande

Alternativement, vous pouvez également rechercher des modules directement depuis la ligne de commande en utilisant les outils nf-core.

```bash
nf-core modules list remote
```

Cela affichera une liste de tous les modules disponibles dans le dépôt nf-core/modules, bien que ce soit un peu moins pratique si vous ne connaissez pas déjà le nom du module que vous recherchez.
Cependant, si vous le connaissez, vous pouvez rediriger la liste vers `grep` pour trouver des modules spécifiques :

```bash
nf-core modules list remote | grep 'cat/cat'
```

??? success "Sortie de la commande"

    ```console
    │ cat/cat
    ```

Gardez simplement à l'esprit que l'approche `grep` ne récupérera que les résultats avec le terme de recherche dans leur nom, ce qui ne fonctionnerait pas pour `cat_cat`.

### 1.3. Obtenir des informations détaillées sur le module

Pour voir des informations détaillées sur un module spécifique depuis la ligne de commande, utilisez la commande `info` :

```bash
nf-core modules info cat/cat
```

Cela affiche la documentation sur le module, y compris ses entrées, ses sorties et des informations d'utilisation de base.

??? success "Sortie de la commande"

    ```console

                                              ,--./,-.
              ___     __   __   __   ___     /,-._.--~\
        |\ | |__  __ /  ` /  \ |__) |__         }  {
        | \| |       \__, \__/ |  \ |___     \`-._,-`-,
                                              `._,._,'

        nf-core/tools version 3.5.2 - https://nf-co.re


    ╭─ Module: cat/cat  ─────────────────────────────────────────────────╮
    │ 🌐 Repository: https://github.com/nf-core/modules.git              │
    │ 🔧 Tools: cat                                                      │
    │ 📖 Description: A module for concatenation of gzipped or           │
    │ uncompressed files                                                 │
    ╰────────────────────────────────────────────────────────────────────╯
                      ╷                                          ╷
    📥 Inputs        │Description                               │Pattern
    ╺━━━━━━━━━━━━━━━━━┿━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┿━━━━━━━╸
    input[0]         │                                          │
    ╶─────────────────┼──────────────────────────────────────────┼───────╴
      meta  (map)     │Groovy Map containing sample information  │
                      │e.g. [ id:'test', single_end:false ]      │
    ╶─────────────────┼──────────────────────────────────────────┼───────╴
      files_in  (file)│List of compressed / uncompressed files   │      *
                      ╵                                          ╵
                          ╷                                 ╷
    📥 Outputs           │Description                      │     Pattern
    ╺━━━━━━━━━━━━━━━━━━━━━┿━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┿━━━━━━━━━━━━╸
    file_out             │                                 │
    ╶─────────────────────┼─────────────────────────────────┼────────────╴
      meta  (map)         │Groovy Map containing sample     │
                          │information                      │
    ╶─────────────────────┼─────────────────────────────────┼────────────╴
      ${prefix}  (file)   │Concatenated file. Will be       │ ${file_out}
                          │gzipped if file_out ends with    │
                          │".gz"                            │
    ╶─────────────────────┼─────────────────────────────────┼────────────╴
    versions_cat         │                                 │
    ╶─────────────────────┼─────────────────────────────────┼────────────╴
      versions_cat (tuple)│Software version information     │
                          ╵                                 ╵

    💻  Installation command: nf-core modules install cat/cat

    ```

Ce sont exactement les mêmes informations que vous pouvez trouver sur le site web.

### 1.4. Installer le module cat/cat

Maintenant que nous avons trouvé le module que nous voulons, nous devons l'ajouter au code source de notre pipeline.

La bonne nouvelle est que le projet nf-core inclut des outils pour faciliter cette partie.
Plus précisément, la commande `nf-core modules install` permet d'automatiser la récupération du code et de le rendre disponible à votre projet en une seule étape.

Accédez au répertoire de votre pipeline et exécutez la commande d'installation :

```bash
cd core-hello
nf-core modules install cat/cat
```

L'outil procédera à l'installation du module.

??? success "Sortie de la commande"

    ```console

                                          ,--./,-.
          ___     __   __   __   ___     /,-._.--~\
    |\ | |__  __ /  ` /  \ |__) |__         }  {
    | \| |       \__, \__/ |  \ |___     \`-._,-`-,
                                          `._,._,'

    nf-core/tools version 3.5.2 - https://nf-co.re


    INFO Installing 'cat/cat'
    INFO Use the following statement to include this module:

        include { CAT_CAT } from '../modules/nf-core/cat/cat/main'
    ```

La commande effectue automatiquement :

- Le téléchargement des fichiers du module dans `modules/nf-core/cat/cat/`
- La mise à jour de `modules.json` pour suivre le module installé
- La fourniture de l'instruction `include` correcte à utiliser dans votre workflow

!!! tip "Astuce"

    Assurez-vous toujours que votre répertoire de travail actuel est la racine de votre projet de pipeline avant d'exécuter la commande d'installation de module.

Vérifions que le module a été installé correctement :

```bash
tree -L 4 modules
```

??? abstract "Contenu du répertoire"

    ```console
    modules
    ├── local
    │   ├── collectGreetings.nf
    │   ├── convertToUpper.nf
    │   ├── cowpy.nf
    │   └── sayHello.nf
    └── nf-core
        └── cat
            └── cat
                ├── environment.yml
                ├── main.nf
                ├── meta.yml
                └── tests

    5 directories, 7 files
    ```

Vous pouvez également vérifier l'installation en demandant à l'utilitaire nf-core de lister les modules installés localement :

```bash
nf-core modules list local
```

??? success "Sortie de la commande"

    ```console
    INFO     Repository type: pipeline
    INFO     Modules installed in '.':

    ┏━━━━━━━━━━━━━┳━━━━━━━━━━━━━━━━━┳━━━━━━━━━━━━━┳━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┳━━━━━━━━━━━━┓
    ┃ Module Name ┃ Repository      ┃ Version SHA ┃ Message                                ┃ Date       ┃
    ┡━━━━━━━━━━━━━╇━━━━━━━━━━━━━━━━━╇━━━━━━━━━━━━━╇━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━╇━━━━━━━━━━━━┩
    │ cat/cat     │ nf-core/modules │ 41dfa3f     │ update meta.yml of all modules (#8747) │ 2025-07-07 │
    └─────────────┴─────────────────┴─────────────┴────────────────────────────────────────┴────────────┘
    ```

Cela confirme que le module `cat/cat` fait maintenant partie du code source de votre projet.

Cependant, pour utiliser réellement le nouveau module, nous devons l'importer dans notre pipeline.

### 1.5. Mettre à jour les importations de module

Remplaçons l'instruction `include` pour le module `collectGreetings` par celle pour `CAT_CAT` dans la section des importations du workflow `workflows/hello.nf`.

Pour rappel, l'outil d'installation de module nous a donné l'instruction exacte à utiliser :

```groovy title="Instruction d'importation produite par la commande d'installation"
include { CAT_CAT } from '../modules/nf-core/cat/cat/main'
```

Notez que la convention nf-core est d'utiliser des majuscules pour les noms de modules lors de leur importation.

Ouvrez [core-hello/workflows/hello.nf](core-hello/workflows/hello.nf) et effectuez la substitution suivante :

=== "Après"

    ```groovy title="core-hello/workflows/hello.nf" linenums="1" hl_lines="11"
    /*
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    */
    include { paramsSummaryMap       } from 'plugin/nf-schema'
    include { softwareVersionsToYAML } from '../subworkflows/nf-core/utils_nfcore_pipeline'
    include { sayHello               } from '../modules/local/sayHello.nf'
    include { convertToUpper         } from '../modules/local/convertToUpper.nf'
    include { cowpy                  } from '../modules/local/cowpy.nf'
    include { CAT_CAT                } from '../modules/nf-core/cat/cat/main'
    ```

=== "Avant"

    ```groovy title="core-hello/workflows/hello.nf" linenums="1" hl_lines="10"
    /*
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    */
    include { paramsSummaryMap       } from 'plugin/nf-schema'
    include { softwareVersionsToYAML } from '../subworkflows/nf-core/utils_nfcore_pipeline'
    include { sayHello               } from '../modules/local/sayHello.nf'
    include { convertToUpper         } from '../modules/local/convertToUpper.nf'
    include { collectGreetings       } from '../modules/local/collectGreetings.nf'
    include { cowpy                  } from '../modules/local/cowpy.nf'
    ```

Remarquez comment le chemin pour le module nf-core diffère des modules locaux :

- **Module nf-core** : `'../modules/nf-core/cat/cat/main'` (référence à `main.nf`)
- **Module local** : `'../modules/local/collectGreetings.nf'` (référence à un fichier unique)

Le module est maintenant disponible pour le workflow, donc tout ce que nous devons faire est de remplacer l'appel à `collectGreetings` pour utiliser `CAT_CAT`. N'est-ce pas ?

Pas si vite.

À ce stade, vous pourriez être tenté de vous lancer et de commencer à éditer le code, mais il vaut la peine de prendre un moment pour examiner attentivement ce que le nouveau module attend et ce qu'il produit.

Nous allons traiter cela comme une section séparée car cela implique un nouveau mécanisme que nous n'avons pas encore couvert : les métadonnées sous forme de map.

!!! note "Note"

    Vous pouvez éventuellement supprimer le fichier `collectGreetings.nf` :

    ```bash
    rm modules/local/collectGreetings.nf
    ```

    Cependant, vous pourriez vouloir le conserver comme référence pour comprendre les différences entre les modules locaux et nf-core.

### À retenir

Vous savez comment trouver un module nf-core et le rendre disponible pour votre projet.

### Et ensuite ?

Évaluer ce qu'un nouveau module requiert et identifier les changements importants nécessaires pour l'intégrer dans un pipeline.

---

## 2. Évaluer les exigences du nouveau module

Plus précisément, nous devons examiner l'**interface** du module, c'est-à-dire ses définitions d'entrée et de sortie, et la comparer à l'interface du module que nous cherchons à remplacer.
Cela nous permettra de déterminer si nous pouvons simplement traiter le nouveau module comme un remplacement direct ou si nous devrons adapter une partie du câblage.

Idéalement, c'est quelque chose que vous devriez faire _avant_ même d'installer le module, mais bon, mieux vaut tard que jamais.
(Pour information, il existe une commande `uninstall` pour se débarrasser des modules que vous décidez de ne plus vouloir.)

!!! note "Note"

    Le processus CAT_CAT inclut une gestion assez intelligente de différents types de compression, d'extensions de fichiers, etc. qui ne sont pas strictement pertinents pour ce que nous essayons de vous montrer ici, donc nous ignorerons la plupart de ces éléments et nous concentrerons uniquement sur les parties importantes.

### 2.1. Comparer les interfaces des deux modules

Pour rappel, voici à quoi ressemble l'interface de notre module `collectGreetings` :

```groovy title="modules/local/collectGreetings.nf (extrait)" linenums="1" hl_lines="6-7 10"
process collectGreetings {

    publishDir 'results', mode: 'copy'

    input:
        path input_files
        val batch_name

    output:
        path "COLLECTED-${batch_name}-output.txt" , emit: outfile
```

Le module `collectGreetings` prend deux entrées :

- `input_files` contient un ou plusieurs fichiers d'entrée à traiter ;
- `batch_name` est une valeur que nous utilisons pour attribuer un nom spécifique à l'exécution au fichier de sortie, qui est une forme de métadonnées.

Une fois terminé, `collectGreetings` produit un seul chemin de fichier, émis avec l'étiquette `outfile`.

En comparaison, l'interface du module `cat/cat` est plus complexe :

```groovy title="modules/nf-core/cat/cat/main.nf (extrait)" linenums="1" hl_lines="11 14"
process CAT_CAT {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/pigz:2.8' :
        'biocontainers/pigz:2.8' }"

    input:
    tuple val(meta), path(files_in)

    output:
    tuple val(meta), path("${prefix}"), emit: file_out
    tuple val("${task.process}"), val("pigz"), eval("pigz --version 2>&1 | sed 's/pigz //g'"), topic: versions, emit: versions_cat
```

Le module CAT_CAT prend une seule entrée, mais cette entrée est un tuple contenant deux choses :

- `meta` est une structure contenant des métadonnées, appelée metamap ;
- `files_in` contient un ou plusieurs fichiers d'entrée à traiter, équivalent à `input_files` de `collectGreetings`.

Une fois terminé, CAT_CAT livre ses sorties en deux parties :

- Un autre tuple contenant le metamap et le fichier de sortie concaténé, émis avec l'étiquette `file_out` ;
- Un tuple de version publié dans le canal topic `versions` pour le suivi des versions de logiciels.

Notez également que par défaut, le fichier de sortie sera nommé en fonction d'un identifiant qui fait partie des métadonnées (code non montré ici).

Cela peut sembler beaucoup à retenir en regardant simplement le code, voici donc un diagramme pour vous aider à visualiser comment tout s'articule.

<figure class="excalidraw">
--8<-- "docs/en/docs/hello_nf-core/img/module_comparison.svg"
</figure>

Vous pouvez voir que les deux modules ont des exigences d'entrée similaires en termes de contenu (un ensemble de fichiers d'entrée plus quelques métadonnées) mais des attentes très différentes quant à la façon dont ce contenu est empaqueté.
En ignorant la sortie versions pour l'instant, leur sortie principale est également équivalente (un fichier concaténé), sauf que CAT_CAT émet également le metamap conjointement avec le fichier de sortie.

Les différences d'empaquetage seront assez faciles à gérer, comme vous le verrez dans un instant.
Cependant, pour comprendre la partie metamap, nous devons vous présenter un contexte supplémentaire.

### 2.2. Comprendre les metamaps

Nous venons de vous dire que le module CAT_CAT attend une map de métadonnées comme partie de son tuple d'entrée.
Prenons quelques minutes pour examiner de plus près ce que c'est.

La **map de métadonnées**, souvent appelée **metamap** en abrégé, est une map de style Groovy contenant des informations sur des unités de données.
Dans le contexte des pipelines Nextflow, les unités de données peuvent être tout ce que vous voulez : des échantillons individuels, des lots d'échantillons ou des ensembles de données entiers.

Par convention, un metamap nf-core est nommé `meta` et contient le champ requis `id`, qui est utilisé pour nommer les sorties et suivre les unités de données.

Par exemple, une map de métadonnées typique pourrait ressembler à ceci :

```groovy title="Exemple de metamap au niveau de l'échantillon"
[id: 'sample1', single_end: false, strandedness: 'forward']
```

Ou dans un cas où les métadonnées sont attachées au niveau du lot :

```groovy title="Exemple de metamap au niveau du lot"
[id: 'batch1', date: '25.10.01']
```

Maintenant, plaçons cela dans le contexte du processus `CAT_CAT`, qui attend que les fichiers d'entrée soient empaquetés dans un tuple avec un metamap, et produit également le metamap comme partie du tuple de sortie.

```groovy title="modules/nf-core/cat/cat/main.nf (extrait)" linenums="1" hl_lines="2 5"
input:
tuple val(meta), path(files_in)

output:
tuple val(meta), path("${prefix}"), emit: file_out
```

En conséquence, chaque unité de données circule dans le pipeline avec les métadonnées pertinentes attachées.
Les processus suivants peuvent alors facilement accéder à ces métadonnées également.

Vous vous souvenez que nous vous avons dit que le fichier produit par `CAT_CAT` sera nommé en fonction d'un identifiant qui fait partie des métadonnées ?
Voici le code pertinent :

```groovy title="modules/nf-core/cat/cat/main.nf (extrait)" linenums="35"
prefix   = task.ext.prefix ?: "${meta.id}${getFileSuffix(file_list[0])}"
```

Cela se traduit approximativement comme suit : si un `prefix` est fourni via le système de paramètres de tâche externes (`task.ext`), utilisez-le pour nommer le fichier de sortie ; sinon créez-en un en utilisant `${meta.id}`, qui correspond au champ `id` dans le metamap.

Vous pouvez imaginer le canal d'entrée entrant dans ce module avec un contenu comme ceci :

```groovy title="Exemple de contenu de canal d'entrée"
ch_input = [[[id: 'batch1', date: '25.10.01'], ['file1A.txt', 'file1B.txt']],
            [[id: 'batch2', date: '25.10.26'], ['file2A.txt', 'file2B.txt']],
            [[id: 'batch3', date: '25.11.14'], ['file3A.txt', 'file3B.txt']]]
```

Puis le contenu du canal de sortie sortant comme ceci :

```groovy title="Exemple de contenu de canal de sortie"
ch_input = [[[id: 'batch1', date: '25.10.01'], 'batch1.txt'],
            [[id: 'batch2', date: '25.10.26'], 'batch2.txt'],
            [[id: 'batch3', date: '25.11.14'], 'batch3.txt']]
```

Comme mentionné précédemment, la configuration d'entrée `tuple val(meta), path(files_in)` est un modèle standard utilisé dans tous les modules nf-core.

Nous espérons que vous commencez à voir à quel point cela peut être utile.
Non seulement cela vous permet de nommer les sorties en fonction des métadonnées, mais vous pouvez également faire des choses comme l'utiliser pour appliquer différentes valeurs de paramètres, et en combinaison avec des opérateurs spécifiques, vous pouvez même regrouper, trier ou filtrer les données au fur et à mesure qu'elles circulent dans le pipeline.

!!! note "En savoir plus sur les métadonnées"

    Pour une introduction complète au travail avec les métadonnées dans les workflows Nextflow, y compris comment lire les métadonnées à partir de samplesheets et les utiliser pour personnaliser le traitement, consultez la quête secondaire [Métadonnées dans les workflows](../side_quests/metadata).

### 2.3. Résumer les changements à effectuer

Sur la base de ce que nous avons examiné, voici les changements majeurs que nous devons apporter à notre pipeline pour utiliser le module `cat/cat` :

- Créer un metamap contenant le nom du lot ;
- Empaqueter le metamap dans un tuple avec l'ensemble des fichiers d'entrée à concaténer (provenant de `convertToUpper`) ;
- Changer l'appel de `collectGreetings()` à `CAT_CAT` ;
- Extraire le fichier de sortie du tuple produit par le processus `CAT_CAT` avant de le passer à `cowpy`.

Cela devrait faire l'affaire ! Maintenant que nous avons un plan, nous sommes prêts à nous lancer.

### À retenir

Vous savez comment évaluer l'interface d'entrée et de sortie d'un nouveau module pour identifier ses exigences, et vous avez appris comment les metamaps sont utilisés par les pipelines nf-core pour garder les métadonnées étroitement associées aux données au fur et à mesure qu'elles circulent dans un pipeline.

### Et ensuite ?

Intégrer le nouveau module dans un workflow.

---

## 3. Intégrer CAT_CAT dans le workflow `hello.nf`

Maintenant que vous savez tout sur les metamaps (ou suffisamment pour les besoins de ce cours, en tout cas), il est temps de réellement implémenter les changements que nous avons décrits ci-dessus.

Par souci de clarté, nous allons décomposer cela et couvrir chaque étape séparément.

!!! note "Note"

    Tous les changements montrés ci-dessous sont apportés à la logique du workflow dans le bloc `main` dans le fichier de workflow `core-hello/workflows/hello.nf`.

### 3.1. Créer une map de métadonnées

Tout d'abord, nous devons créer une map de métadonnées pour `CAT_CAT`, en gardant à l'esprit que les modules nf-core exigent que le metamap contienne au moins un champ `id`.

Puisque nous n'avons besoin d'aucune autre métadonnée, nous pouvons rester simple et utiliser quelque chose comme ceci :

```groovy title="Exemple de syntaxe"
def cat_meta = [id: 'test']
```

Sauf que nous ne voulons pas coder en dur la valeur `id` ; nous voulons utiliser la valeur du paramètre `params.batch`.
Donc le code devient :

```groovy title="Exemple de syntaxe"
def cat_meta = [id: params.batch]
```

Oui, c'est littéralement aussi simple que cela de créer un metamap de base.

Ajoutons ces lignes après l'appel à `convertToUpper`, en supprimant l'appel à `collectGreetings` :

=== "Après"

    ```groovy title="core-hello/workflows/hello.nf" linenums="26" hl_lines="7-8"
        // émettre une salutation
        sayHello(ch_samplesheet)

        // convertir la salutation en majuscules
        convertToUpper(sayHello.out)

        // créer une map de métadonnées avec le nom du lot comme ID
        def cat_meta = [ id: params.batch ]

        // générer de l'art ASCII des salutations avec cowpy
        cowpy(collectGreetings.out.outfile, params.character)
    ```

=== "Avant"

    ```groovy title="core-hello/workflows/hello.nf" linenums="26" hl_lines="7-8"
        // émettre une salutation
        sayHello(ch_samplesheet)

        // convertir la salutation en majuscules
        convertToUpper(sayHello.out)

        // rassembler toutes les salutations dans un seul fichier
        collectGreetings(convertToUpper.out.collect(), params.batch)

        // générer de l'art ASCII des salutations avec cowpy
        cowpy(collectGreetings.out.outfile, params.character)
    ```

Cela crée une map de métadonnées simple où l'`id` est défini sur le nom de notre lot (qui sera `test` lors de l'utilisation du profil test).

### 3.2. Créer un canal avec des tuples de métadonnées

Ensuite, transformez le canal de fichiers en un canal de tuples contenant des métadonnées et des fichiers :

=== "Après"

    ```groovy title="core-hello/workflows/hello.nf" linenums="26" hl_lines="10-11"
        // émettre une salutation
        sayHello(ch_samplesheet)

        // convertir la salutation en majuscules
        convertToUpper(sayHello.out)

        // créer une map de métadonnées avec le nom du lot comme ID
        def cat_meta = [ id: params.batch ]

        // créer un canal avec des métadonnées et des fichiers au format tuple
        ch_for_cat = convertToUpper.out.collect().map { files -> tuple(cat_meta, files) }

        // générer de l'art ASCII des salutations avec cowpy
        cowpy(collectGreetings.out.outfile, params.character)
    ```

=== "Avant"

    ```groovy title="core-hello/workflows/hello.nf" linenums="26"
        // émettre une salutation
        sayHello(ch_samplesheet)

        // convertir la salutation en majuscules
        convertToUpper(sayHello.out)

        // créer une map de métadonnées avec le nom du lot comme ID
        def cat_meta = [ id: params.batch ]

        // générer de l'art ASCII des salutations avec cowpy
        cowpy(collectGreetings.out.outfile, params.character)
    ```

La ligne que nous avons ajoutée accomplit deux choses :

- `.collect()` rassemble tous les fichiers de la sortie `convertToUpper` dans une seule liste
- `#!groovy .map { files -> tuple(cat_meta, files) }` crée un tuple de `[métadonnées, fichiers]` dans le format attendu par `CAT_CAT`

C'est tout ce que nous devons faire pour configurer le tuple d'entrée pour `CAT_CAT`.

### 3.3. Appeler le module CAT_CAT

Maintenant, appelez `CAT_CAT` sur le canal nouvellement créé :

=== "Après"

    ```groovy title="core-hello/workflows/hello.nf" linenums="26" hl_lines="13-14"
        // émettre une salutation
        sayHello(ch_samplesheet)

        // convertir la salutation en majuscules
        convertToUpper(sayHello.out)

        // créer une map de métadonnées avec le nom du lot comme ID
        def cat_meta = [ id: params.batch ]

        // créer un canal avec des métadonnées et des fichiers au format tuple
        ch_for_cat = convertToUpper.out.collect().map { files -> tuple(cat_meta, files) }

        // concaténer les salutations
        CAT_CAT(ch_for_cat)

        // générer de l'art ASCII des salutations avec cowpy
        cowpy(collectGreetings.out.outfile, params.character)
    ```

=== "Avant"

    ```groovy title="core-hello/workflows/hello.nf" linenums="26"
        // émettre une salutation
        sayHello(ch_samplesheet)

        // convertir la salutation en majuscules
        convertToUpper(sayHello.out)

        // créer une map de métadonnées avec le nom du lot comme ID
        def cat_meta = [ id: params.batch ]

        // créer un canal avec des métadonnées et des fichiers au format tuple
        ch_for_cat = convertToUpper.out.collect().map { files -> tuple(cat_meta, files) }

        // générer de l'art ASCII des salutations avec cowpy
        cowpy(collectGreetings.out.outfile, params.character)
    ```

Cela complète la partie la plus délicate de cette substitution, mais nous n'avons pas tout à fait terminé : nous devons encore mettre à jour la façon dont nous passons la sortie concaténée au processus `cowpy`.

### 3.4. Extraire le fichier de sortie du tuple pour `cowpy`

Auparavant, le processus `collectGreetings` produisait simplement un fichier que nous pouvions passer directement à `cowpy`.
Cependant, le processus `CAT_CAT` produit un tuple qui inclut le metamap en plus du fichier de sortie.

Puisque `cowpy` n'accepte pas encore les tuples de métadonnées (nous corrigerons cela dans la prochaine partie du cours), nous devons extraire le fichier de sortie du tuple produit par `CAT_CAT` avant de le transmettre à `cowpy` :

=== "Après"

    ```groovy title="core-hello/workflows/hello.nf" linenums="26" hl_lines="16-17 20"
        // émettre une salutation
        sayHello(ch_samplesheet)

        // convertir la salutation en majuscules
        convertToUpper(sayHello.out)

        // créer une map de métadonnées avec le nom du lot comme ID
        def cat_meta = [ id: params.batch ]

        // créer un canal avec des métadonnées et des fichiers au format tuple
        ch_for_cat = convertToUpper.out.collect().map { files -> tuple(cat_meta, files) }

        // concaténer les salutations
        CAT_CAT(ch_for_cat)

        // extraire le fichier du tuple puisque cowpy n'utilise pas encore les métadonnées
        ch_for_cowpy = CAT_CAT.out.file_out.map{ meta, file -> file }

        // générer de l'art ASCII des salutations avec cowpy
        cowpy(ch_for_cowpy, params.character)
    ```

=== "Avant"

    ```groovy title="core-hello/workflows/hello.nf" linenums="26" hl_lines="17"
        // émettre une salutation
        sayHello(ch_samplesheet)

        // convertir la salutation en majuscules
        convertToUpper(sayHello.out)

        // créer une map de métadonnées avec le nom du lot comme ID
        def cat_meta = [ id: params.batch ]

        // créer un canal avec des métadonnées et des fichiers au format tuple
        ch_for_cat = convertToUpper.out.collect().map { files -> tuple(cat_meta, files) }

        // concaténer les salutations
        CAT_CAT(ch_for_cat)

        // générer de l'art ASCII des salutations avec cowpy
        cowpy(collectGreetings.out.outfile, params.character)
    ```

L'opération `#!groovy .map { meta, file -> file }` extrait le fichier du tuple `[métadonnées, fichier]` produit par `CAT_CAT` dans un nouveau canal, `ch_for_cowpy`.

Ensuite, il suffit de passer `ch_for_cowpy` à `cowpy` au lieu de `collectGreetings.out.outfile` dans cette dernière ligne.

!!! note "Note"

    Dans la prochaine partie du cours, nous mettrons à jour `cowpy` pour qu'il fonctionne directement avec les tuples de métadonnées, donc cette étape d'extraction ne sera plus nécessaire.

### 3.5. Tester le workflow

Testons que le workflow fonctionne avec le module `cat/cat` nouvellement intégré :

```bash
nextflow run . --outdir core-hello-results -profile test,docker --validate_params false
```

Cela devrait s'exécuter assez rapidement.

??? success "Sortie de la commande"

    ```console
    N E X T F L O W ~ version 25.04.3

        Launching `./main.nf` [evil_pike] DSL2 - revision: b9e9b3b8de

        Input/output options
          input                     : /workspaces/training/hello-nf-core/core-hello/assets/greetings.csv
          outdir                    : core-hello-results

        Institutional config options
          config_profile_name       : Test profile
          config_profile_description: Minimal test dataset to check pipeline function

        Generic options
          validate_params           : false
          trace_report_suffix       : 2025-10-30_18-50-58

        Core Nextflow options
          runName                   : evil_pike
          containerEngine           : docker
          launchDir                 : /workspaces/training/hello-nf-core/core-hello
          workDir                   : /workspaces/training/hello-nf-core/core-hello/work
          projectDir                : /workspaces/training/hello-nf-core/core-hello
          userName                  : root
          profile                   : test,docker
          configFiles               : /workspaces/training/hello-nf-core/core-hello/nextflow.config

        !! Only displaying parameters that differ from the pipeline defaults !!
        ------------------------------------------------------
        executor >  local (8)
        [b3/f005fd] CORE_HELLO:HELLO:sayHello (3)       [100%] 3 of 3 ✔
        [08/f923d0] CORE_HELLO:HELLO:convertToUpper (3) [100%] 3 of 3 ✔
        [34/3729a9] CORE_HELLO:HELLO:CAT_CAT (test)     [100%] 1 of 1 ✔
        [24/df918a] CORE_HELLO:HELLO:cowpy              [100%] 1 of 1 ✔
        -[core/hello] Pipeline completed successfully-
    ```

Remarquez que `CAT_CAT` apparaît maintenant dans la liste d'exécution des processus au lieu de `collectGreetings`.

Et voilà ! Nous utilisons maintenant un module robuste et maintenu par la communauté au lieu d'un code personnalisé de qualité prototype pour cette étape du pipeline.

### À retenir

Vous savez maintenant comment :

- Trouver et installer des modules nf-core
- Évaluer les exigences d'un module nf-core
- Créer une map de métadonnées simple pour l'utiliser avec un module nf-core
- Intégrer un module nf-core dans votre workflow

### Et ensuite ?

Apprendre à adapter vos modules locaux pour suivre les conventions nf-core.
Nous vous montrerons également comment créer de nouveaux modules nf-core à partir d'un modèle en utilisant les outils nf-core.
