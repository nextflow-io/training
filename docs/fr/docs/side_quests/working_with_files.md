# Traitement des fichiers en entrée

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Traduction assistée par IA - [en savoir plus et suggérer des améliorations](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

Les flux de travail d'analyse scientifique impliquent souvent le traitement d'un grand nombre de fichiers.
Nextflow fournit des outils puissants pour gérer les fichiers efficacement, vous aidant à organiser et traiter vos données avec un minimum de code.

### Objectifs d'apprentissage

Dans cette quête secondaire, nous allons explorer comment Nextflow gère les fichiers, des opérations de base aux techniques plus avancées pour travailler avec des collections de fichiers.
Vous apprendrez à extraire des métadonnées à partir de noms de fichiers, ce qui est une exigence courante dans les pipelines d'analyse scientifique.

À la fin de cette quête secondaire, vous serez capable de :

- Créer des objets Path à partir de chaînes de caractères représentant des chemins de fichiers en utilisant la méthode `file()` de Nextflow
- Accéder aux attributs de fichiers tels que le nom, l'extension et le répertoire parent
- Gérer de manière transparente les fichiers locaux et distants en utilisant des URI
- Utiliser les canaux pour automatiser la gestion des fichiers avec `channel.fromPath()` et `channel.fromFilePairs()`
- Extraire et structurer les métadonnées à partir de noms de fichiers en utilisant la manipulation de chaînes
- Grouper les fichiers associés en utilisant la correspondance de motifs et les expressions glob
- Intégrer les opérations sur fichiers dans les processus Nextflow avec une gestion appropriée des entrées
- Organiser les sorties de processus en utilisant des structures de répertoires pilotées par les métadonnées

Ces compétences vous aideront à construire des flux de travail capables de gérer différents types d'entrées de fichiers avec une grande flexibilité.

### Prérequis

Avant d'entreprendre cette quête secondaire, vous devriez :

- Avoir terminé le tutoriel [Hello Nextflow](../../hello_nextflow/) ou un cours équivalent pour débutants.
- Être à l'aise avec l'utilisation des concepts et mécanismes de base de Nextflow (processus, canaux, opérateurs)

<!-- I removed the suggestion to do the metamaps SQ first because that works more naturally after -->

---

## 0. Démarrage

#### Ouvrir l'espace de code de formation

Si vous ne l'avez pas encore fait, assurez-vous d'ouvrir l'environnement de formation comme décrit dans la [Configuration de l'Environnement](../envsetup/index.md).

[![Open in GitHub Codespaces](https://github.com/codespaces/badge.svg)](https://codespaces.new/nextflow-io/training?quickstart=1&ref=master)

#### Se déplacer dans le répertoire du projet

Déplaçons-nous dans le répertoire où se trouvent les fichiers pour ce tutoriel.

```bash
cd side-quests/working_with_files
```

Vous pouvez configurer VSCode pour se concentrer sur ce répertoire :

```bash
code .
```

#### Examiner les matériaux

Vous trouverez un fichier de flux de travail simple appelé `main.nf`, un répertoire `modules` contenant deux fichiers de modules, et un répertoire `data` contenant quelques fichiers de données d'exemple.

??? abstract "Contenu du répertoire"

    ```console
    .
    ├── data
    │   ├── patientA_rep1_normal_R1_001.fastq.gz
    │   ├── patientA_rep1_normal_R2_001.fastq.gz
    │   ├── patientA_rep1_tumor_R1_001.fastq.gz
    │   ├── patientA_rep1_tumor_R2_001.fastq.gz
    │   ├── patientA_rep2_normal_R1_001.fastq.gz
    │   ├── patientA_rep2_normal_R2_001.fastq.gz
    │   ├── patientA_rep2_tumor_R1_001.fastq.gz
    │   ├── patientA_rep2_tumor_R2_001.fastq.gz
    │   ├── patientB_rep1_normal_R1_001.fastq.gz
    │   ├── patientB_rep1_normal_R2_001.fastq.gz
    │   ├── patientB_rep1_tumor_R1_001.fastq.gz
    │   ├── patientB_rep1_tumor_R2_001.fastq.gz
    │   ├── patientC_rep1_normal_R1_001.fastq.gz
    │   ├── patientC_rep1_normal_R2_001.fastq.gz
    │   ├── patientC_rep1_tumor_R1_001.fastq.gz
    │   └── patientC_rep1_tumor_R2_001.fastq.gz
    ├── main.nf
    └── modules
        ├── analyze_reads.nf
        └── count_lines.nf
    ```

Ce répertoire contient des données de séquençage en paires provenant de trois patients (A, B, C).

Pour chaque patient, nous avons des échantillons de type `tumor` (provenant généralement de biopsies tumorales) ou `normal` (prélevés sur des tissus sains ou du sang).
Si vous n'êtes pas familier avec l'analyse du cancer, sachez simplement que cela correspond à un modèle expérimental qui utilise des échantillons appariés tumeur/normal pour effectuer des analyses contrastives.

Pour le patient A spécifiquement, nous avons deux ensembles de réplicats techniques (répétitions).

Les fichiers de données de séquençage sont nommés avec une convention typique `_R1_` et `_R2_` pour ce qu'on appelle les 'lectures avant' et 'lectures arrière'.

_Ne vous inquiétez pas si vous n'êtes pas familier avec ce plan expérimental, ce n'est pas critique pour comprendre ce tutoriel._

#### Examiner l'exercice

Votre défi est d'écrire un flux de travail Nextflow qui va :

1. **Charger** les fichiers d'entrée en utilisant les méthodes de gestion de fichiers de Nextflow
2. **Extraire** les métadonnées (ID patient, réplicat, type d'échantillon) de la structure du nom de fichier
3. **Grouper** les fichiers appariés (R1/R2) ensemble en utilisant `channel.fromFilePairs()`
4. **Traiter** les fichiers avec un module d'analyse fourni
5. **Organiser** les sorties dans une structure de répertoires basée sur les métadonnées extraites

#### Liste de vérification de préparation

Vous pensez être prêt·e à vous lancer ?

- [ ] Je comprends l'objectif de ce cours et ses prérequis
- [ ] Mon espace de code est opérationnel
- [ ] J'ai défini mon répertoire de travail de manière appropriée
- [ ] Je comprends l'exercice

Si vous pouvez cocher toutes les cases, vous êtes prêt·e.

---

## 1. Opérations de base sur les fichiers

### 1.1. Identifier le type d'un objet avec `.class`

Regardez le fichier de flux de travail `main.nf` :

```groovy title="main.nf" linenums="1"
#!/usr/bin/env nextflow

workflow {

    // Crée un objet Path à partir d'une chaîne de chemin
    myFile = 'data/patientA_rep1_normal_R1_001.fastq.gz'

    println "${myFile} is of class ${myFile.class}"
}
```

Il s'agit d'un mini-flux de travail (sans processus) qui fait référence à un seul chemin de fichier dans son workflow, puis l'imprime dans la console, avec sa classe.

??? info "Qu'est-ce que `.class` ?"

    Dans Nextflow, `.class` nous indique quel type d'objet nous manipulons. C'est comme demander "quel genre de chose est-ce ?" pour savoir s'il s'agit d'une chaîne de caractères, d'un nombre, d'un fichier ou d'autre chose.
    Cela nous aidera à illustrer la différence entre une simple chaîne de caractères et un objet Path dans les sections suivantes.

Exécutons le flux de travail :

```bash
nextflow run main.nf
```

??? success "Sortie de la commande"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `main.nf` [romantic_chandrasekhar] DSL2 - revision: 5a4a89bc3a

    data/patientA_rep1_normal_R1_001.fastq.gz is of class java.lang.String
    ```

Comme vous pouvez le voir, Nextflow a imprimé le chemin sous forme de chaîne exactement comme nous l'avons écrit.

Il s'agit simplement d'une sortie texte ; Nextflow n'a encore rien fait de spécial avec.
Nous avons également confirmé qu'en ce qui concerne Nextflow, il s'agit seulement d'une chaîne (de classe `java.lang.String`).
Cela a du sens, puisque nous n'avons pas encore indiqué à Nextflow qu'il correspond à un fichier.

### 1.2. Créer un objet Path avec file()

Nous pouvons indiquer à Nextflow comment gérer les fichiers en créant des [objets Path](https://www.nextflow.io/docs/latest/reference/stdlib-types.html#path) à partir de chaînes de chemins.

Dans notre flux de travail, nous pouvons convertir le chemin sous forme de chaîne `data/patientA_rep1_normal_R1_001.fastq.gz` en objet Path en utilisant la méthode `file()`, qui donne accès aux propriétés et opérations sur les fichiers.

Modifiez le fichier `main.nf` pour envelopper la chaîne avec `file()` comme suit :

=== "Après"

    ```groovy title="main.nf" linenums="5" hl_lines="2"
        // Crée un objet Path à partir d'une chaîne de chemin
        myFile = file('data/patientA_rep1_normal_R1_001.fastq.gz')

        println "${myFile} is of class ${myFile.class}"
    ```

=== "Avant"

    ```groovy title="main.nf" linenums="5" hl_lines="2"
        // Crée un objet Path à partir d'une chaîne de chemin
        myFile = 'data/patientA_rep1_normal_R1_001.fastq.gz'

        println "${myFile} is of class ${myFile.class}"
    ```

Maintenant, exécutez à nouveau le flux de travail :

```bash
nextflow run main.nf
```

??? success "Sortie de la commande"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `main.nf` [kickass_coulomb] DSL2 - revision: 5af44b1b59

    /workspaces/training/side-quests/working_with_files/data/patientA_rep1_normal_R1_001.fastq.gz is of class class sun.nio.fs.UnixPath
    ```

Cette fois, vous voyez le chemin absolu complet au lieu du chemin relatif que nous avons fourni en entrée.

Nextflow a converti notre chaîne en objet Path et l'a résolu vers l'emplacement réel du fichier sur le système.
Le chemin du fichier sera maintenant absolu, comme dans `/workspaces/training/side-quests/working_with_files/data/patientA_rep1_normal_R1_001.fastq.gz`.

Notez également que la classe de l'objet Path est `sun.nio.fs.UnixPath` : c'est la façon dont Nextflow représente les fichiers locaux.
Comme nous le verrons plus tard, les fichiers distants auront des noms de classe différents (comme `nextflow.file.http.XPath` pour les fichiers HTTP), mais ils fonctionnent tous exactement de la même manière et peuvent être utilisés de manière identique dans vos flux de travail.

!!! tip

    **La différence clé :**

    - **Chaîne de chemin** : Juste du texte que Nextflow traite comme des caractères
    - **Objet Path** : Une référence de fichier intelligente avec laquelle Nextflow peut travailler

    Pensez-y ainsi : une chaîne de chemin est comme écrire une adresse sur papier, tandis qu'un objet Path est comme avoir l'adresse chargée dans un appareil GPS qui sait comment naviguer jusqu'à là et peut vous donner des détails sur le trajet.

### 1.3. Accéder aux attributs de fichier

Pourquoi est-ce utile ? Eh bien, maintenant que Nextflow comprend que `myFile` est un objet Path et pas seulement une chaîne, nous pouvons accéder aux différents attributs de l'objet Path.

Mettons à jour notre flux de travail pour imprimer les attributs de fichier intégrés :

=== "Après"

    ```groovy title="main.nf" linenums="5" hl_lines="4-9"
        // Crée un objet Path à partir d'une chaîne de chemin
        myFile = file('data/patientA_rep1_normal_R1_001.fastq.gz')

        // Affiche les attributs du fichier
        println "File object class: ${myFile.class}"
        println "File name: ${myFile.name}"
        println "Simple name: ${myFile.simpleName}"
        println "Extension: ${myFile.extension}"
        println "Parent directory: ${myFile.parent}"
    ```

=== "Avant"

    ```groovy title="main.nf" linenums="5" hl_lines="4"
        // Crée un objet Path à partir d'une chaîne de chemin
        myFile = file('data/patientA_rep1_normal_R1_001.fastq.gz')

        println "${myFile} is of class ${myFile.class}"
    ```

Exécutez le flux de travail :

```bash
nextflow run main.nf
```

??? success "Sortie de la commande"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `main.nf` [ecstatic_ampere] DSL2 - revision: f3fa3dcb48

    File object class: sun.nio.fs.UnixPath
    File name: patientA_rep1_normal_R1_001.fastq.gz
    Simple name: patientA_rep1_normal_R1_001
    Extension: gz
    Parent directory: /workspaces/training/side-quests/working_with_files/data
    ```

Vous voyez les différents attributs de fichier imprimés dans la console ci-dessus.

### 1.4. Fournir le fichier à un processus

La différence entre les chaînes et les objets Path devient critique lorsque vous commencez à construire des flux de travail réels avec des processus.
Jusqu'à présent, nous avons vérifié que Nextflow traite maintenant notre fichier d'entrée comme un fichier, mais voyons si nous pouvons réellement exécuter quelque chose sur ce fichier dans un processus.

#### 1.4.1. Importer le processus et examiner le code

Nous vous fournissons un module de processus pré-écrit appelé `COUNT_LINES` qui prend un fichier en entrée et compte combien de lignes il contient.

Pour utiliser le processus dans le flux de travail, vous devez simplement ajouter une instruction include avant le bloc workflow :

=== "Après"

    ```groovy title="main.nf" linenums="1" hl_lines="3"
    #!/usr/bin/env nextflow

    include { COUNT_LINES } from './modules/count_lines.nf'

    workflow {
    ```

=== "Avant"

    ```groovy title="main.nf" linenums="1"
    #!/usr/bin/env nextflow

    workflow {
    ```

Vous pouvez ouvrir le fichier de module pour examiner son code :

```groovy title="modules/count_lines.nf" linenums="1"
#!/usr/bin/env nextflow

process COUNT_LINES {
    debug true

    input:
    path input_file

    script:
    """
    set -o pipefail
    echo "Processing file: $input_file"
    gzip -dc $input_file | wc -l
    """
}
```

Comme vous pouvez le voir, c'est un petit script assez simple qui décompresse le fichier et compte combien de lignes il contient.

??? info "Que fait `debug true` ?"

    La directive `debug true` dans la définition du processus fait en sorte que Nextflow imprime la sortie de votre script (comme le nombre de lignes "40") directement dans le journal d'exécution.
    Sans cela, vous ne verriez que le statut d'exécution du processus mais pas la sortie réelle de votre script.

    Pour plus d'informations sur le débogage des flux de travail Nextflow, consultez la quête secondaire [Débogage des Flux de Travail Nextflow](debugging.md).

#### 1.4.2. Ajouter un appel à `COUNT_LINES`

Maintenant que le processus est disponible pour le flux de travail, nous pouvons ajouter un appel au processus `COUNT_LINES` pour l'exécuter sur le fichier d'entrée.

Effectuez les modifications suivantes dans le flux de travail :

=== "Après"

    ```groovy title="main.nf" linenums="7" hl_lines="11-12"
        // Crée un objet Path à partir d'une chaîne de chemin
        myFile = file('data/patientA_rep1_normal_R1_001.fastq.gz')

        // Affiche les attributs du fichier
        println "File object class: ${myFile.class}"
        println "File name: ${myFile.name}"
        println "Simple name: ${myFile.simpleName}"
        println "Extension: ${myFile.extension}"
        println "Parent directory: ${myFile.parent}"

        // Compte les lignes dans le fichier
        COUNT_LINES(myFile)
    ```

=== "Avant"

    ```groovy title="main.nf" linenums="7" hl_lines="4-9"
        // Crée un objet Path à partir d'une chaîne de chemin
        myFile = file('data/patientA_rep1_normal_R1_001.fastq.gz')

        // Affiche les attributs du fichier
        println "File object class: ${myFile.class}"
        println "File name: ${myFile.name}"
        println "Simple name: ${myFile.simpleName}"
        println "Extension: ${myFile.extension}"
        println "Parent directory: ${myFile.parent}"
    ```

Et maintenant exécutez le flux de travail :

```bash
nextflow run main.nf
```

??? success "Sortie de la commande"

    ```console
     N E X T F L O W   ~  version 25.10.2

    Launching `main.nf` [cheeky_hypatia] DSL2 - revision: 281d13c414

    File object class: class sun.nio.fs.UnixPath
    File name: patientA_rep1_normal_R1_001.fastq.gz
    Simple name: patientA_rep1_normal_R1_001
    Extension: gz
    Parent directory: /workspaces/training/side-quests/working_with_files/data
    executor >  local (1)
    [e9/341c05] COUNT_LINES [100%] 1 of 1 ✔
    Processing file: /workspaces/training/side-quests/working_with_files/data/patientA_rep1_normal_R1_001.fastq.gz
    40
    ```

Cela montre que nous sommes capables d'opérer sur le fichier de manière appropriée à l'intérieur d'un processus.

Plus précisément, Nextflow a effectué avec succès les opérations suivantes :

- Mis en place le fichier dans le répertoire de travail
- Décompressé le fichier .gz
- Compté les lignes (40 lignes dans ce cas)
- Terminé sans erreur

La clé de cette opération fluide est que nous indiquons explicitement à Nextflow que notre entrée est un fichier et doit être traitée comme telle.

### 1.5. Dépanner les erreurs de base sur les entrées de fichiers

Cela déroute souvent les nouveaux venus à Nextflow, alors prenons quelques minutes pour examiner ce qui se passe lorsque vous le faites mal.

Il y a deux endroits principaux où vous pouvez vous tromper dans la gestion des fichiers : au niveau du flux de travail et au niveau du processus.

#### 1.5.1. Erreur au niveau du flux de travail

Voyons ce qui se passe si nous revenons à traiter le fichier comme une chaîne lorsque nous spécifions l'entrée dans le bloc workflow.

Effectuez les modifications suivantes dans le flux de travail, en veillant à commenter les instructions d'impression spécifiques aux chemins :

=== "Après"

    ```groovy title="main.nf" linenums="7" hl_lines="2 6-11"
        // Crée un objet Path à partir d'une chaîne de chemin
        myFile = 'data/patientA_rep1_normal_R1_001.fastq.gz'

        // Affiche les attributs du fichier
        println "File object class: ${myFile.class}"
        /*
        println "File name: ${myFile.name}"
        println "Simple name: ${myFile.simpleName}"
        println "Extension: ${myFile.extension}"
        println "Parent directory: ${myFile.parent}"
        */

        // Compte les lignes dans le fichier
        COUNT_LINES(myFile)
    ```

=== "Avant"

    ```groovy title="main.nf" linenums="7" hl_lines="4-9"
        // Crée un objet Path à partir d'une chaîne de chemin
        myFile = file('data/patientA_rep1_normal_R1_001.fastq.gz')

        // Affiche les attributs du fichier
        println "File object class: ${myFile.class}"
        println "File name: ${myFile.name}"
        println "Simple name: ${myFile.simpleName}"
        println "Extension: ${myFile.extension}"
        println "Parent directory: ${myFile.parent}"

        // Compte les lignes dans le fichier
        COUNT_LINES(myFile)
    ```

Et maintenant exécutez le flux de travail :

```bash
nextflow run main.nf
```

??? failure "Sortie de la commande"

    ```console
     N E X T F L O W   ~  version 25.10.2

    Launching `main.nf` [friendly_goodall] DSL2 - revision: ae50609b20

    [-        ] COUNT_LINES -
    ERROR ~ Error executing process > 'COUNT_LINES'

    Caused by:
      Not a valid path value: 'data/patientA_rep1_normal_R1_001.fastq.gz'



    Tip: view the complete command output by changing to the process work dir and entering the command `cat .command.out`

    -- Check '.nextflow.log' file for details
    ```

Voici la partie importante :

```console
Not a valid path value: 'data/patientA_rep1_normal_R1_001.fastq.gz'
```

Lorsque vous spécifiez une entrée `path`, Nextflow valide que vous passez des références de fichiers réelles, pas seulement des chaînes.
Cette erreur vous indique que `'data/patientA_rep1_normal_R1_001.fastq.gz'` n'est pas une valeur de chemin valide car c'est une chaîne, pas un objet Path.

Nextflow a immédiatement détecté le problème et s'est arrêté avant même de démarrer le processus.

#### 1.5.2. Erreur au niveau du processus

L'autre endroit où nous pourrions oublier de spécifier que nous voulons que Nextflow traite l'entrée comme un fichier est dans la définition du processus.

!!! warning "Conservez l'erreur du flux de travail de 1.5.1"

    Pour que ce test fonctionne correctement, conservez le flux de travail dans son état défaillant (utilisant une simple chaîne au lieu de `file()`).
    Lorsqu'il est combiné avec `val` dans le processus, cela produit l'erreur montrée ci-dessous.

Effectuez la modification suivante dans le module :

=== "Après"

    ```groovy title="modules/count_lines.nf" linenums="3" hl_lines="5"
    process COUNT_LINES {
        debug true

        input:
        val input_file
    ```

=== "Avant"

    ```groovy title="modules/count_lines.nf" linenums="3" hl_lines="5"
    process COUNT_LINES {
        debug true

        input:
        path input_file
    ```

Et maintenant exécutez à nouveau le flux de travail :

```bash
nextflow run main.nf
```

??? failure "Sortie de la commande"

    ```console
     N E X T F L O W   ~  version 25.10.2

    Launching `main.nf` [soggy_golick] DSL2 - revision: ae50609b20

    executor >  local (1)
    [b3/b3023c] COUNT_LINES [  0%] 0 of 1 ✘
    ERROR ~ Error executing process > 'COUNT_LINES'

    Caused by:
      Process `COUNT_LINES` terminated with an error exit status (1)


    Command executed:

      set -o pipefail
      echo "Processing file: data/patientA_rep1_normal_R1_001.fastq.gz"
      gzip -dc data/patientA_rep1_normal_R1_001.fastq.gz | wc -l

    Command exit status:
      1

    Command output:
      Processing file: data/patientA_rep1_normal_R1_001.fastq.gz
      0

    Command error:
      Processing file: data/patientA_rep1_normal_R1_001.fastq.gz
      gzip: data/patientA_rep1_normal_R1_001.fastq.gz: No such file or directory
      0

    Work dir:
      /workspaces/training/side-quests/working_with_files/work/b3/b3023cb2ccb986851301d8e369e79f

    Tip: when you have fixed the problem you can continue the execution adding the option `-resume` to the run command line

    -- Check '.nextflow.log' file for details
    ```

Cela montre beaucoup de détails sur l'erreur car le processus est configuré pour afficher des informations de débogage, comme indiqué précédemment.

Voici les sections les plus pertinentes :

```console
Command executed:

  set -o pipefail
  echo "Processing file: data/patientA_rep1_normal_R1_001.fastq.gz"
  gzip -dc data/patientA_rep1_normal_R1_001.fastq.gz | wc -l
```

```console
Command error:
  Processing file: data/patientA_rep1_normal_R1_001.fastq.gz
  gzip: data/patientA_rep1_normal_R1_001.fastq.gz: No such file or directory
  0
```

Cela indique que le système n'a pas pu trouver le fichier ; cependant, si vous regardez le chemin, il y a bien un fichier avec ce nom à cet emplacement.

Lorsque nous avons exécuté cela, Nextflow a transmis la valeur de chaîne au script, mais il n'a pas _mis en place_ le fichier réel dans le répertoire de travail.
Ainsi, le processus a essayé d'utiliser la chaîne relative, `data/patientA_rep1_normal_R1_001.fastq.gz`, mais ce fichier n'existe pas dans le répertoire de travail du processus.

Pris ensemble, ces deux exemples vous montrent à quel point il est important d'indiquer à Nextflow si une entrée doit être gérée comme un fichier.

!!! note

    Assurez-vous de revenir en arrière et de corriger les deux erreurs intentionnelles avant de continuer à la section suivante.

### Enseignements clés

- Chaînes de chemin vs objets Path : Les chaînes sont juste du texte, les objets Path sont des références de fichiers intelligentes
- La méthode `file()` convertit un chemin sous forme de chaîne en objet Path avec lequel Nextflow peut travailler
- Vous pouvez accéder aux propriétés de fichier comme `name`, `simpleName`, `extension` et `parent` [en utilisant les attributs de fichier](https://www.nextflow.io/docs/latest/working-with-files.html#getting-file-attributes)
- L'utilisation d'objets Path au lieu de chaînes permet à Nextflow de gérer correctement les fichiers dans votre flux de travail
- Résultats des entrées de processus : Une gestion appropriée des fichiers nécessite des objets Path, et non des chaînes, pour garantir que les fichiers sont correctement mis en place et accessibles pour être utilisés par les processus.

---

## 2. Utilisation de fichiers distants

L'une des fonctionnalités clés de Nextflow est la capacité de basculer de manière transparente entre les fichiers locaux (sur la même machine) et les fichiers distants accessibles via Internet.

Si vous le faites correctement, vous ne devriez jamais avoir besoin de modifier la logique de votre flux de travail pour accommoder des fichiers provenant de différents emplacements.
Tout ce que vous devez faire pour utiliser un fichier distant est de spécifier le préfixe approprié dans le chemin du fichier lorsque vous le fournissez au flux de travail.

Par exemple, `/path/to/data` n'a pas de préfixe, indiquant qu'il s'agit d'un chemin de fichier local 'normal', alors que `s3://path/to/data` inclut le préfixe `s3://`, indiquant qu'il est situé dans le stockage d'objets S3 d'Amazon.

De nombreux protocoles différents sont pris en charge :

- HTTP(S)/FTP (http://, https://, ftp://)
- Amazon S3 (s3://)
- Azure Blob Storage (az://)
- Google Cloud Storage (gs://)

Pour utiliser l'un de ces protocoles, spécifiez simplement le préfixe pertinent dans la chaîne, qui est alors techniquement appelée Uniform Resource Identifier (URI) au lieu de chemin de fichier.
Nextflow gérera l'authentification et la mise en place des fichiers au bon endroit, le téléchargement ou le chargement et toutes les autres opérations sur fichiers auxquelles vous vous attendez.

La force principale de ce système est qu'il nous permet de basculer entre les environnements sans modifier aucune logique de pipeline.
Par exemple, vous pouvez développer avec un petit ensemble de tests local avant de passer à un ensemble de tests à grande échelle situé dans un stockage distant simplement en changeant l'URI.

### 2.1. Utiliser un fichier depuis Internet

Testons cela en remplaçant le chemin local que nous fournissons à notre flux de travail par un chemin HTTPS pointant vers une copie des mêmes données stockées dans Github.

!!! warning

    Cela ne fonctionnera que si vous avez une connexion Internet active.

Ouvrez à nouveau `main.nf` et modifiez le chemin d'entrée comme suit :

=== "Après"

    ```groovy title="main.nf" linenums="2" hl_lines="2"
        // Utilise un fichier distant depuis Internet
        myFile = file('https://raw.github.com/nextflow-io/training/master/side-quests/working_with_files/data/patientA_rep1_normal_R1_001.fastq.gz')

        // Affiche les attributs du fichier
        println "File object class: ${myFile.class}"
        println "File name: ${myFile.name}"
        println "Simple name: ${myFile.simpleName}"
        println "Extension: ${myFile.extension}"
        println "Parent directory: ${myFile.parent}"
    ```

=== "Avant"

    ```groovy title="main.nf" linenums="2" hl_lines="2"
        // Crée un objet Path à partir d'une chaîne de chemin
        myFile = file('data/patientA_rep1_normal_R1_001.fastq.gz')

        // Affiche les attributs du fichier
        println "File object class: ${myFile.class}"
        println "File name: ${myFile.name}"
        println "Simple name: ${myFile.simpleName}"
        println "Extension: ${myFile.extension}"
        println "Parent directory: ${myFile.parent}"
    ```

Exécutons le flux de travail :

```bash
nextflow run main.nf
```

??? success "Sortie de la commande"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `main.nf` [insane_swartz] DSL2 - revision: fff18abe6d

    File object class: class nextflow.file.http.XPath
    File name: patientA_rep1_normal_R1_001.fastq.gz
    Simple name: patientA_rep1_normal_R1_001
    Extension: gz
    Parent directory: /nextflow-io/training/master/side-quests/working_with_files/data
    executor >  local (1)
    [8a/2ab7ca] COUNT_LINES [100%] 1 of 1 ✔
    Processing file: patientA_rep1_normal_R1_001.fastq.gz
    40
    ```

Ça marche ! Vous pouvez voir que très peu de choses ont changé.

La seule différence dans la sortie de la console est que la classe de l'objet path est maintenant `nextflow.file.http.XPath`, alors que pour le chemin local la classe était `sun.nio.fs.UnixPath`.
Vous n'avez pas besoin de mémoriser ces classes ; nous mentionnons simplement cela pour démontrer que Nextflow identifie et gère les différents emplacements de manière appropriée.

En coulisses, Nextflow a téléchargé le fichier dans un répertoire de mise en place situé dans le répertoire de travail.
Ce fichier mis en place peut ensuite être traité comme un fichier local et créé en lien symbolique dans le répertoire de processus pertinent.

Vous pouvez vérifier que cela s'est produit ici en regardant le contenu du répertoire de travail situé à la valeur de hachage du processus.

??? abstract "Contenu du répertoire de travail"

    Si le hachage du processus était `8a/2ab7ca`, vous pourriez explorer le répertoire de travail :

    ```console
    $ ls -la work/8a/2ab7ca*/
    total 16
    drwxr-xr-x  6 user  staff   192 Jan 28 10:00 .
    drwxr-xr-x  3 user  staff    96 Jan 28 10:00 ..
    -rw-r--r--  1 user  staff     0 Jan 28 10:00 .command.begin
    -rw-r--r--  1 user  staff   127 Jan 28 10:00 .command.sh
    lrwxr-xr-x  1 user  staff    89 Jan 28 10:00 patientA_rep1_normal_R1_001.fastq.gz -> /path/to/work/stage/.../patientA_rep1_normal_R1_001.fastq.gz
    ```

    Le lien symbolique pointe vers une copie mise en place du fichier distant que Nextflow a téléchargé automatiquement.

Notez que pour les fichiers plus volumineux, l'étape de téléchargement prendra du temps supplémentaire par rapport à l'exécution sur des fichiers locaux.
Cependant, Nextflow vérifie s'il a déjà une copie mise en place pour éviter les téléchargements inutiles.
Donc, si vous réexécutez sur le même fichier et que vous n'avez pas supprimé le fichier mis en place, Nextflow utilisera la copie mise en place.

Cela montre à quel point il est facile de basculer entre les données locales et distantes en utilisant Nextflow, ce qui est une fonctionnalité clé de Nextflow.

!!! note

    La seule exception importante à ce principe est que vous ne pouvez pas utiliser de motifs glob ou de chemins de répertoires avec HTTPS car HTTPS ne peut pas lister plusieurs fichiers, vous devez donc spécifier des URL de fichiers exactes.
    Cependant, d'autres protocoles de stockage tels que le stockage blob (`s3://`, `az://`, `gs://`) peuvent utiliser à la fois des globs et des chemins de répertoires.

    Voici comment vous pourriez utiliser des motifs glob avec le stockage cloud :

    ```groovy title="Exemples de stockage cloud (non exécutables dans cet environnement)"
    // S3 avec motifs glob - correspondrait à plusieurs fichiers
    ch_s3_files = channel.fromPath('s3://my-bucket/data/*.fastq.gz')

    // Azure Blob Storage avec motifs glob
    ch_azure_files = channel.fromPath('az://container/data/patient*_R{1,2}.fastq.gz')

    // Google Cloud Storage avec motifs glob
    ch_gcs_files = channel.fromPath('gs://bucket/data/sample_*.fastq.gz')
    ```

    Nous vous montrerons comment travailler avec des globs dans la pratique dans la section suivante.

### 2.2. Revenir au fichier local

Nous allons revenir à l'utilisation de nos fichiers d'exemple locaux pour le reste de cette quête secondaire, alors remettons l'entrée du flux de travail sur le fichier d'origine :

=== "Après"

    ```groovy title="main.nf" linenums="2" hl_lines="2"
        // Crée un objet Path à partir d'une chaîne de chemin
        myFile = file('data/patientA_rep1_normal_R1_001.fastq.gz')

        // Affiche les attributs du fichier
        println "File object class: ${myFile.class}"
        println "File name: ${myFile.name}"
        println "Simple name: ${myFile.simpleName}"
        println "Extension: ${myFile.extension}"
        println "Parent directory: ${myFile.parent}"
    ```

=== "Avant"

    ```groovy title="main.nf" linenums="2" hl_lines="2"
        // Crée un objet Path à partir d'une chaîne de chemin
        myFile = file('https://raw.github.com/nextflow-io/training/master/side-quests/working_with_files/data/patientA_rep1_normal_R1_001.fastq.gz')

        // Affiche les attributs du fichier
        println "File object class: ${myFile.class}"
        println "File name: ${myFile.name}"
        println "Simple name: ${myFile.simpleName}"
        println "Extension: ${myFile.extension}"
        println "Parent directory: ${myFile.parent}"
    ```

### Enseignements clés

- Les données distantes sont accessibles en utilisant un URI (HTTP, FTP, S3, Azure, Google Cloud)
- Nextflow téléchargera et mettra en place automatiquement les données au bon endroit, tant que ces chemins sont fournis aux processus
- N'écrivez pas de logique pour télécharger ou charger des fichiers distants !
- Les fichiers locaux et distants produisent différents types d'objets mais fonctionnent de manière identique
- **Important** : HTTP/HTTPS fonctionne uniquement avec des fichiers uniques (pas de motifs glob)
- Le stockage cloud (S3, Azure, GCS) prend en charge à la fois les fichiers uniques et les motifs glob
- Vous pouvez basculer de manière transparente entre les sources de données locales et distantes sans modifier la logique du code (tant que le protocole prend en charge vos opérations requises)

---

## 3. Utilisation de la fabrique de canal `fromPath()`

Jusqu'à présent, nous avons travaillé avec un seul fichier à la fois, mais dans Nextflow, nous allons généralement vouloir créer un canal d'entrée avec plusieurs fichiers d'entrée à traiter.

Une façon naïve de le faire serait de combiner la méthode `file()` avec [`channel.of()`](https://www.nextflow.io/docs/latest/reference/channel.html#of) comme ceci :

```groovy title="Exemple de syntaxe"
ch_files = channel.of([file('data/patientA_rep1_normal_R1_001.fastq.gz')],
                      [file('data/patientA_rep1_normal_R1_001.fastq.gz')])
```

Cela fonctionne, mais c'est maladroit.

!!! tip "Quand utiliser `file()` vs `channel.fromPath()`"

    - Utilisez `file()` lorsque vous avez besoin d'un seul objet Path pour une manipulation directe (vérifier si un fichier existe, lire ses attributs, ou le passer à un seul appel de processus)
    - Utilisez `channel.fromPath()` lorsque vous avez besoin d'un canal pouvant contenir plusieurs fichiers, surtout avec des motifs glob, ou lorsque les fichiers vont circuler à travers plusieurs processus

C'est là qu'intervient [`channel.fromPath()`](https://www.nextflow.io/docs/latest/reference/channel.html#frompath) : une fabrique de canal pratique qui regroupe toutes les fonctionnalités dont nous avons besoin pour générer un canal à partir d'une ou plusieurs chaînes de fichiers statiques ainsi que de motifs glob.

### 3.1. Ajouter la fabrique de canal

Mettons à jour notre flux de travail pour utiliser `channel.fromPath`.

=== "Après"

    ```groovy title="main.nf" linenums="7" hl_lines="1-3"
        // Charge les fichiers avec channel.fromPath
        ch_files = channel.fromPath('data/patientA_rep1_normal_R1_001.fastq.gz')
        ch_files.view { myFile -> "Found file: $myFile" }

        // Affiche les attributs du fichier
        /* Comment these out for now, we'll come back to them!
        println "File object class: ${myFile.class}"
        println "File name: ${myFile.name}"
        println "Simple name: ${myFile.simpleName}"
        println "Extension: ${myFile.extension}"
        println "Parent directory: ${myFile.parent}"
        */

        // Compte les lignes dans le fichier
        // COUNT_LINES(myFile)
    ```

=== "Avant"

    ```groovy title="main.nf" linenums="7" hl_lines="1-2"
        // Crée un objet Path à partir d'une chaîne de chemin
        myFile = file('data/patientA_rep1_normal_R1_001.fastq.gz')

        // Affiche les attributs du fichier
        println "File object class: ${myFile.class}"
        println "File name: ${myFile.name}"
        println "Simple name: ${myFile.simpleName}"
        println "Extension: ${myFile.extension}"
        println "Parent directory: ${myFile.parent}"

        // Compte les lignes dans le fichier
        COUNT_LINES(myFile)
    ```

Nous avons également commenté le code qui imprime les attributs pour le moment, et ajouté une instruction `.view` pour imprimer juste le nom de fichier à la place.

Exécutez le flux de travail :

```bash
nextflow run main.nf
```

??? success "Sortie de la commande"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `main.nf` [grave_meucci] DSL2 - revision: b09964a583

    Found file: /workspaces/training/side-quests/working_with_files/data/patientA_rep1_normal_R1_001.fastq.gz
    ```

Comme vous pouvez le voir, le chemin du fichier est chargé en tant qu'objet de type `Path` dans le canal.
C'est similaire à ce qu'aurait fait `file()`, sauf que maintenant nous avons un canal dans lequel nous pouvons charger plus de fichiers si nous le souhaitons.

L'utilisation de `channel.fromPath()` est un moyen pratique de créer un nouveau canal peuplé par une liste de fichiers.

### 3.2. Voir les attributs des fichiers dans le canal

Dans notre première utilisation de la fabrique de canal, nous avons simplifié le code et imprimé simplement le nom de fichier.

Revenons à l'impression des attributs de fichier complets :

=== "Après"

    ```groovy title="main.nf" linenums="7" hl_lines="3-9 12"
        // Charge les fichiers avec channel.fromPath
        ch_files = channel.fromPath('data/patientA_rep1_normal_R1_001.fastq.gz')
        ch_files.view { myFile ->
            println "File object class: ${myFile.class}"
            println "File name: ${myFile.name}"
            println "Simple name: ${myFile.simpleName}"
            println "Extension: ${myFile.extension}"
            println "Parent directory: ${myFile.parent}"
        }

        // Compte les lignes dans le fichier
        COUNT_LINES(ch_files)
    ```

=== "Avant"

    ```groovy title="main.nf" linenums="7" hl_lines="3"
        // Charge les fichiers avec channel.fromPath
        ch_files = channel.fromPath('data/patientA_rep1_normal_R1_001.fastq.gz')
        ch_files.view { myFile -> "Found file: $myFile" }

        // Compte les lignes dans le fichier
        // COUNT_LINES(ch_files)
    ```

Nous réactivons également l'appel au processus `COUNT_LINES` pour vérifier que le traitement des fichiers fonctionne toujours correctement avec notre approche basée sur les canaux.

Exécutez le flux de travail :

```bash
nextflow run main.nf
```

??? success "Sortie de la commande"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `main.nf` [furious_swanson] DSL2 - revision: c35c34950d

    executor >  local (1)
    [9d/6701a6] COUNT_LINES (1) [100%] 1 of 1 ✔
    File object class: sun.nio.fs.UnixPath
    File name: patientA_rep1_normal_R1_001.fastq.gz
    Simple name: patientA_rep1_normal_R1_001
    Extension: gz
    Parent directory: /workspaces/training/side-quests/working_with_files/data
    Processing file: patientA_rep1_normal_R1_001.fastq.gz
    40
    ```

Et voilà, mêmes résultats qu'avant mais maintenant nous avons le fichier dans un canal, donc nous pouvons en ajouter plus.

### 3.3. Utiliser un glob pour correspondre à plusieurs fichiers

Il y a plusieurs façons de charger plus de fichiers dans le canal.
Ici, nous allons vous montrer comment utiliser des motifs glob, qui sont un moyen pratique de faire correspondre et récupérer des noms de fichiers et de répertoires basés sur des caractères génériques.
Le processus de correspondance de ces motifs est appelé "globbing" ou "expansion de nom de fichier".

!!! note

    Comme indiqué précédemment, Nextflow prend en charge le globbing pour gérer les fichiers d'entrée et de sortie dans la majorité des cas, sauf avec les chemins de fichiers HTTPS car HTTPS ne peut pas lister plusieurs fichiers.

Disons que nous voulons récupérer les deux fichiers d'une paire de fichiers associés à un patient donné, `patientA` :

```console
patientA_rep1_normal_R1_001.fastq.gz
patientA_rep1_normal_R2_001.fastq.gz
```

Puisque la seule différence entre les noms de fichiers est le numéro de réplicat, _c'est-à-dire_ le numéro après `R`, nous pouvons utiliser le caractère générique `*` pour représenter le numéro comme suit :

```console
patientA_rep1_normal_R*_001.fastq.gz
```

C'est le motif glob dont nous avons besoin.

Maintenant, tout ce que nous devons faire est de mettre à jour le chemin de fichier dans la fabrique de canal pour utiliser ce motif glob comme suit :

=== "Après"

    ```groovy title="main.nf" linenums="7"
      // Charge les fichiers avec channel.fromPath
        ch_files = channel.fromPath('data/patientA_rep1_normal_R*_001.fastq.gz')
    ```

=== "Avant"

    ```groovy title="main.nf" linenums="7"
      // Charge les fichiers avec channel.fromPath
        ch_files = channel.fromPath('data/patientA_rep1_normal_R1_001.fastq.gz')
    ```

Nextflow reconnaîtra automatiquement qu'il s'agit d'un motif glob et le gérera de manière appropriée.

Exécutez le flux de travail pour tester cela :

```bash
nextflow run main.nf
```

??? success "Sortie de la commande"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `main.nf` [boring_sammet] DSL2 - revision: d2aa789c9a

    executor >  local (2)
    [3c/a65de5] COUNT_LINES (2) [100%] 2 of 2 ✔
    File object class: class sun.nio.fs.UnixPath
    File name: patientA_rep1_normal_R1_001.fastq.gz
    Simple name: patientA_rep1_normal_R1_001
    Extension: gz
    Parent directory: /workspaces/training/side-quests/working_with_files/data
    File object class: class sun.nio.fs.UnixPath
    File name: patientA_rep1_normal_R2_001.fastq.gz
    Simple name: patientA_rep1_normal_R2_001
    Extension: gz
    Parent directory: /workspaces/training/side-quests/working_with_files/data
    Processing file: patientA_rep1_normal_R1_001.fastq.gz
    40

    Processing file: patientA_rep1_normal_R2_001.fastq.gz
    40
    ```

Comme vous pouvez le voir, nous avons maintenant deux objets Path dans notre canal, ce qui montre que Nextflow a effectué l'expansion de nom de fichier correctement et a chargé et traité les deux fichiers comme prévu.

En utilisant cette méthode, nous pouvons récupérer autant ou aussi peu de fichiers que nous le souhaitons simplement en changeant le motif glob. Si nous le rendions plus généreux, par exemple en remplaçant toutes les parties variables des noms de fichiers par `*` (_par exemple_ `data/patient*_rep*_*_R*_001.fastq.gz`), nous pourrions récupérer tous les fichiers d'exemple dans le répertoire `data`.

### Enseignements clés

- `channel.fromPath()` crée un canal avec des fichiers correspondant à un motif
- Chaque fichier est émis en tant qu'élément séparé dans le canal
- Nous pouvons utiliser un motif glob pour correspondre à plusieurs fichiers
- Les fichiers sont automatiquement convertis en objets Path avec des attributs complets
- La méthode `.view()` permet l'inspection du contenu du canal

---

## 4. Extraction de métadonnées de base à partir de noms de fichiers

Dans la plupart des domaines scientifiques, il est très courant d'avoir des métadonnées encodées dans les noms des fichiers qui contiennent les données.
Par exemple, en bio-informatique, les fichiers contenant des données de séquençage sont souvent nommés de manière à encoder des informations sur l'échantillon, la condition, le réplicat et le numéro de lecture.

Si les noms de fichiers sont construits selon une convention cohérente, vous pouvez extraire ces métadonnées de manière standardisée et les utiliser dans le cours de votre analyse.
C'est un grand 'si', bien sûr, et vous devriez être très prudent chaque fois que vous vous fiez à la structure des noms de fichiers ; mais la réalité est que cette approche est très largement utilisée, alors jetons un œil à comment cela se fait dans Nextflow.

Dans le cas de nos données d'exemple, nous savons que les noms de fichiers incluent des métadonnées structurées de manière cohérente.
Par exemple, le nom de fichier `patientA_rep1_normal_R2_001` encode ce qui suit :

- ID patient : `patientA`
- ID réplicat : `rep1`
- type d'échantillon : `normal` (par opposition à `tumor`)
- ensemble de lectures : `R1` (par opposition à `R2`)

Nous allons modifier notre flux de travail pour récupérer ces informations en trois étapes :

1. Récupérer le `simpleName` du fichier, qui inclut les métadonnées
2. Séparer les métadonnées en utilisant une méthode appelée `tokenize()`
3. Utiliser une map pour organiser les métadonnées

!!! warning

    Vous ne devriez jamais encoder d'informations sensibles dans les noms de fichiers, telles que les noms de patients ou d'autres caractéristiques d'identification, car cela peut compromettre la confidentialité des patients ou d'autres restrictions de sécurité pertinentes.

### 4.1. Récupérer le `simpleName`

Le `simpleName` est un attribut de fichier qui correspond au nom de fichier dépouillé de son chemin et de son extension.

Effectuez les modifications suivantes dans le flux de travail :

=== "Après"

    ```groovy title="main.nf" linenums="7" hl_lines="3-6"
        // Charge les fichiers avec channel.fromPath
        ch_files = channel.fromPath('data/patientA_rep1_normal_R*_001.fastq.gz')
        ch_files.map { myFile ->
            [ myFile.simpleName, myFile ]
        }
        .view()
    ```

=== "Avant"

    ```groovy title="main.nf" linenums="7" hl_lines="3-9"
        // Charge les fichiers avec channel.fromPath
        ch_files = channel.fromPath('data/patientA_rep1_normal_R*_001.fastq.gz')
        ch_files.view { myFile ->
            println "File object class: ${myFile.class}"
            println "File name: ${myFile.name}"
            println "Simple name: ${myFile.simpleName}"
            println "Extension: ${myFile.extension}"
            println "Parent directory: ${myFile.parent}"
        }
    ```

Cela récupère le `simpleName` et l'associe à l'objet fichier complet en utilisant une opération `map()`.

Exécutez le flux de travail pour tester qu'il fonctionne :

```bash
nextflow run main.nf
```

??? success "Sortie de la commande"

    ```console hl_lines="7-8"
    N E X T F L O W   ~  version 25.10.2

    Launching `main.nf` [suspicious_mahavira] DSL2 - revision: ae8edc4e48

    executor >  local (2)
    [e9/55774b] COUNT_LINES (2) [100%] 2 of 2 ✔
    [patientA_rep1_normal_R2_001, /workspaces/training/side-quests/working_with_files/data/patientA_rep1_normal_R2_001.fastq.gz]
    [patientA_rep1_normal_R1_001, /workspaces/training/side-quests/working_with_files/data/patientA_rep1_normal_R1_001.fastq.gz]
    Processing file: patientA_rep1_normal_R1_001.fastq.gz
    40

    Processing file: patientA_rep1_normal_R2_001.fastq.gz
    40
    ```

Chaque élément dans le canal est maintenant un tuple contenant le `simpleName` et l'objet fichier d'origine.

### 4.2. Extraire les métadonnées du `simplename`

À ce stade, les métadonnées que nous voulons sont intégrées dans le `simplename`, mais nous ne pouvons pas accéder directement aux éléments individuels.
Nous devons donc diviser le `simplename` en ses composants.
Heureusement, ces composants sont simplement séparés par des traits de soulignement dans le nom de fichier d'origine, nous pouvons donc appliquer une méthode Nextflow courante appelée `tokenize()` qui est parfaite pour cette tâche.

Effectuez les modifications suivantes dans le flux de travail :

=== "Après"

    ```groovy title="main.nf" linenums="7" hl_lines="4"
        // Charge les fichiers avec channel.fromPath
        ch_files = channel.fromPath('data/patientA_rep1_normal_R*_001.fastq.gz')
        ch_files.map { myFile ->
            [ myFile.simpleName.tokenize('_'), myFile ]
        }
    ```

=== "Avant"

    ```groovy title="main.nf" linenums="7" hl_lines="4"
        // Charge les fichiers avec channel.fromPath
        ch_files = channel.fromPath('data/patientA_rep1_normal_R*_001.fastq.gz')
        ch_files.map { myFile ->
            [ myFile.simpleName, myFile ]
        }
    ```

La méthode `tokenize()` divisera la chaîne `simpleName` partout où elle trouve des traits de soulignement, et retournera une liste contenant les sous-chaînes.

Exécutez le flux de travail :

```bash
nextflow run main.nf
```

??? success "Sortie de la commande"

    ```console hl_lines="7-8"
    N E X T F L O W   ~  version 25.10.2

    Launching `main.nf` [gigantic_gauss] DSL2 - revision: a39baabb57

    executor >  local (2)
    [e7/da2f4b] COUNT_LINES (2) [100%] 2 of 2 ✔
    [[patientA, rep1, normal, R2, 001], /workspaces/training/side-quests/working_with_files/data/patientA_rep1_normal_R2_001.fastq.gz]
    [[patientA, rep1, normal, R1, 001], /workspaces/training/side-quests/working_with_files/data/patientA_rep1_normal_R1_001.fastq.gz]
    Processing file: patientA_rep1_normal_R2_001.fastq.gz
    40

    Processing file: patientA_rep1_normal_R1_001.fastq.gz
    40
    ```

Maintenant le tuple pour chaque élément dans notre canal contient la liste de métadonnées (_par exemple_ `[patientA, rep1, normal, R1, 001]`) et l'objet fichier d'origine.

C'est génial !
Nous avons décomposé les informations de notre patient d'une seule chaîne en une liste de chaînes.
Nous pouvons maintenant gérer chaque partie de l'information du patient séparément.

### 4.3. Utiliser une map pour organiser les métadonnées

Nos métadonnées ne sont qu'une liste plate pour le moment.
Elle est assez facile à utiliser mais difficile à lire.

```console
[patientA, rep1, normal, R1, 001]
```

Quel est l'élément à l'index 3 ? Pouvez-vous le dire sans vous référer à l'explication originale de la structure des métadonnées ?

C'est une excellente opportunité d'utiliser un stockage clé-valeur, où chaque élément a un ensemble de clés et leurs valeurs associées, vous pouvez donc facilement vous référer à chaque clé pour obtenir la valeur correspondante.

Dans notre exemple, cela signifie passer de cette organisation :

```groovy
data = [patientA, 1, normal, R1]

println data[3]
```

À celle-ci :

```groovy
data = [id: patientA, replicate: 1, type: normal, readNum: 1]

println data.readNum
```

Dans Nextflow, cela s'appelle une [map](https://nextflow.io/docs/latest/script.html#maps).

Convertissons maintenant notre liste plate en map.
Effectuez les modifications suivantes dans le flux de travail :

=== "Après"

    ```groovy title="main.nf" linenums="7" hl_lines="4-13"
        // Charge les fichiers avec channel.fromPath
        ch_files = channel.fromPath('data/patientA_rep1_normal_R*_001.fastq.gz')
        ch_files.map { myFile ->
            def (patient, replicate, type, readNum) = myFile.simpleName.tokenize('_')
            [
              [
                id: patient,
                replicate: replicate.replace('rep', ''),
                type: type,
                readNum: readNum.replace('R', ''),
              ],
              myFile
            ]
        }
    ```

===
