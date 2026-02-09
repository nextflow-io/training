# Traitement des fichiers en entrée

Les workflows d'analyse scientifique impliquent souvent le traitement d'un grand nombre de fichiers.
Nextflow fournit des outils puissants pour gérer les fichiers de manière efficace, vous aidant à organiser et traiter vos données avec un minimum de code.

### Objectifs d'apprentissage

Dans cette quête secondaire, nous explorerons comment Nextflow gère les fichiers, des opérations de base aux techniques plus avancées pour travailler avec des collections de fichiers.
Vous apprendrez à extraire des métadonnées à partir des noms de fichiers, ce qui est une exigence courante dans les pipelines d'analyse scientifique.

À la fin de cette quête secondaire, vous serez capable de :

- Créer des objets Path à partir de chaînes de chemins de fichiers en utilisant la méthode `file()` de Nextflow
- Accéder aux attributs de fichiers tels que le nom, l'extension et le répertoire parent
- Gérer de manière transparente les fichiers locaux et distants en utilisant des URI
- Utiliser des canaux pour automatiser la gestion des fichiers avec `channel.fromPath()` et `channel.fromFilePairs()`
- Extraire et structurer les métadonnées à partir des noms de fichiers en utilisant la manipulation de chaînes
- Regrouper les fichiers associés en utilisant la correspondance de motifs et les expressions glob
- Intégrer les opérations sur fichiers dans les processus Nextflow avec une gestion appropriée des entrées
- Organiser les sorties de processus en utilisant des structures de répertoires basées sur les métadonnées

Ces compétences vous aideront à construire des workflows capables de gérer différents types d'entrées de fichiers avec une grande flexibilité.

### Prérequis

Avant d'entreprendre cette quête secondaire, vous devriez :

- Avoir complété le tutoriel [Hello Nextflow](../../hello_nextflow/) ou un cours équivalent pour débutant·es.
- Être à l'aise avec l'utilisation des concepts et mécanismes de base de Nextflow (processus, canaux, opérateurs)

<!-- I removed the suggestion to do the metamaps SQ first because that works more naturally after -->

---

## 0. Commencer

#### Ouvrir l'espace de code de formation

Si vous ne l'avez pas encore fait, assurez-vous d'ouvrir l'environnement de formation comme décrit dans la [Configuration de l'environnement](../envsetup/index.md).

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

#### Examiner le matériel

Vous trouverez un fichier de workflow simple appelé `main.nf`, un répertoire `modules` contenant deux fichiers de modules, et un répertoire `data` contenant quelques fichiers de données d'exemple.

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

Ce répertoire contient des données de séquençage paired-end de trois patient·es (A, B, C).

Pour chaque patient·e, nous avons des échantillons de type `tumor` (provenant généralement de biopsies tumorales) ou `normal` (prélevés sur des tissus sains ou du sang).
Si vous n'êtes pas familier·ère avec l'analyse du cancer, sachez simplement que cela correspond à un modèle expérimental qui utilise des échantillons appariés tumeur/normal pour effectuer des analyses contrastives.

Pour le patient A spécifiquement, nous avons deux ensembles de réplicats techniques (répétitions).

Les fichiers de données de séquençage sont nommés avec une convention typique `_R1_` et `_R2_` pour ce qu'on appelle les 'lectures forward' et 'lectures reverse'.

_Ne vous inquiétez pas si vous n'êtes pas familier·ère avec ce plan expérimental, ce n'est pas critique pour comprendre ce tutoriel._

#### Examiner l'exercice

Votre défi est d'écrire un workflow Nextflow qui va :

1. **Charger** les fichiers d'entrée en utilisant les méthodes de gestion de fichiers de Nextflow
2. **Extraire** les métadonnées (ID patient, réplicat, type d'échantillon) de la structure du nom de fichier
3. **Regrouper** les fichiers appariés (R1/R2) ensemble en utilisant `channel.fromFilePairs()`
4. **Traiter** les fichiers avec un module d'analyse fourni
5. **Organiser** les sorties dans une structure de répertoires basée sur les métadonnées extraites

#### Liste de vérification de préparation

Pensez-vous être prêt·e à vous lancer ?

- [ ] Je comprends l'objectif de ce cours et ses prérequis
- [ ] Mon espace de code est opérationnel
- [ ] J'ai défini mon répertoire de travail de manière appropriée
- [ ] Je comprends l'exercice

Si vous pouvez cocher toutes les cases, vous êtes prêt·e à commencer.

---

## 1. Opérations de base sur les fichiers

### 1.1. Identifier le type d'un objet avec `.class`

Jetez un œil au fichier de workflow `main.nf` :

```groovy title="main.nf" linenums="1"
#!/usr/bin/env nextflow

workflow {

    // Create a Path object from a string path
    myFile = 'data/patientA_rep1_normal_R1_001.fastq.gz'

    println "${myFile} is of class ${myFile.class}"
}
```

Il s'agit d'un mini-workflow (sans aucun processus) qui fait référence à un seul chemin de fichier dans son workflow, puis l'affiche dans la console, avec sa classe.

??? info "Qu'est-ce que `.class` ?"

    Dans Nextflow, `.class` nous indique quel type d'objet nous manipulons. C'est comme demander "quel genre de chose est-ce ?" pour savoir s'il s'agit d'une chaîne, d'un nombre, d'un fichier ou d'autre chose.
    Cela nous aidera à illustrer la différence entre une simple chaîne et un objet Path dans les sections suivantes.

Exécutons le workflow :

```bash
nextflow run main.nf
```

??? success "Sortie de la commande"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `main.nf` [romantic_chandrasekhar] DSL2 - revision: 5a4a89bc3a

    data/patientA_rep1_normal_R1_001.fastq.gz is of class java.lang.String
    ```

Comme vous pouvez le voir, Nextflow a affiché le chemin de la chaîne exactement comme nous l'avons écrit.

Il s'agit simplement d'une sortie texte ; Nextflow n'a encore rien fait de spécial avec.
Nous avons également confirmé qu'en ce qui concerne Nextflow, il s'agit seulement d'une chaîne (de classe `java.lang.String`).
Cela a du sens, puisque nous n'avons pas encore dit à Nextflow qu'elle correspond à un fichier.

### 1.2. Créer un objet Path avec file()

Nous pouvons indiquer à Nextflow comment gérer les fichiers en créant des [objets Path](https://www.nextflow.io/docs/latest/reference/stdlib-types.html#path) à partir de chaînes de chemins.

Dans notre workflow, nous pouvons convertir le chemin de chaîne `data/patientA_rep1_normal_R1_001.fastq.gz` en un objet Path en utilisant la méthode `file()`, qui fournit un accès aux propriétés et opérations sur les fichiers.

Modifiez le `main.nf` pour envelopper la chaîne avec `file()` comme suit :

=== "Après"

    ```groovy title="main.nf" linenums="5" hl_lines="2"
        // Create a Path object from a string path
        myFile = file('data/patientA_rep1_normal_R1_001.fastq.gz')

        println "${myFile} is of class ${myFile.class}"
    ```

=== "Avant"

    ```groovy title="main.nf" linenums="5" hl_lines="2"
        // Create a Path object from a string path
        myFile = 'data/patientA_rep1_normal_R1_001.fastq.gz'

        println "${myFile} is of class ${myFile.class}"
    ```

Maintenant, exécutez à nouveau le workflow :

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

Nextflow a converti notre chaîne en un objet Path et l'a résolu vers l'emplacement réel du fichier sur le système.
Le chemin du fichier sera maintenant absolu, comme dans `/workspaces/training/side-quests/working_with_files/data/patientA_rep1_normal_R1_001.fastq.gz`.

Notez également que la classe de l'objet Path est `sun.nio.fs.UnixPath` : c'est la façon dont Nextflow représente les fichiers locaux.
Comme nous le verrons plus tard, les fichiers distants auront des noms de classe différents (comme `nextflow.file.http.XPath` pour les fichiers HTTP), mais ils fonctionnent tous exactement de la même manière et peuvent être utilisés de manière identique dans vos workflows.

!!! tip

    **La différence clé :**

    - **Chaîne de chemin** : Juste du texte que Nextflow traite comme des caractères
    - **Objet Path** : Une référence de fichier intelligente avec laquelle Nextflow peut travailler

    Pensez-y comme ceci : une chaîne de chemin est comme écrire une adresse sur papier, tandis qu'un objet Path est comme avoir l'adresse chargée dans un appareil GPS qui sait comment y naviguer et peut vous donner des détails sur le trajet.

### 1.3. Accéder aux attributs de fichier

Pourquoi est-ce utile ? Eh bien, maintenant que Nextflow comprend que `myFile` est un objet Path et pas seulement une chaîne, nous pouvons accéder aux différents attributs de l'objet Path.

Mettons à jour notre workflow pour afficher les attributs de fichier intégrés :

=== "Après"

    ```groovy title="main.nf" linenums="5" hl_lines="4-9"
        // Create a Path object from a string path
        myFile = file('data/patientA_rep1_normal_R1_001.fastq.gz')

        // Print file attributes
        println "File object class: ${myFile.class}"
        println "File name: ${myFile.name}"
        println "Simple name: ${myFile.simpleName}"
        println "Extension: ${myFile.extension}"
        println "Parent directory: ${myFile.parent}"
    ```

=== "Avant"

    ```groovy title="main.nf" linenums="5" hl_lines="4"
        // Create a Path object from a string path
        myFile = file('data/patientA_rep1_normal_R1_001.fastq.gz')

        println "${myFile} is of class ${myFile.class}"
    ```

Exécutez le workflow :

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

Vous voyez les différents attributs de fichier affichés dans la console ci-dessus.

### 1.4. Fournir le fichier à un processus

La différence entre les chaînes et les objets Path devient critique lorsque vous commencez à construire des workflows réels avec des processus.
Jusqu'à présent, nous avons vérifié que Nextflow traite maintenant notre fichier d'entrée comme un fichier, mais voyons si nous pouvons réellement exécuter quelque chose sur ce fichier dans un processus.

#### 1.4.1. Importer le processus et examiner le code

Nous vous fournissons un module de processus pré-écrit appelé `COUNT_LINES` qui prend un fichier en entrée et compte combien de lignes il contient.

Pour utiliser le processus dans le workflow, vous devez simplement ajouter une instruction include avant le bloc workflow :

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

    La directive `debug true` dans la définition du processus fait en sorte que Nextflow affiche la sortie de votre script (comme le nombre de lignes "40") directement dans le journal d'exécution.
    Sans cela, vous ne verriez que le statut d'exécution du processus mais pas la sortie réelle de votre script.

    Pour plus d'informations sur le débogage des workflows Nextflow, consultez la quête secondaire [Debugging Nextflow Workflows](debugging.md).

#### 1.4.2. Ajouter un appel à `COUNT_LINES`

Maintenant que le processus est disponible pour le workflow, nous pouvons ajouter un appel au processus `COUNT_LINES` pour l'exécuter sur le fichier d'entrée.

Effectuez les modifications suivantes dans le workflow :

=== "Après"

    ```groovy title="main.nf" linenums="7" hl_lines="11-12"
        // Create a Path object from a string path
        myFile = file('data/patientA_rep1_normal_R1_001.fastq.gz')

        // Print file attributes
        println "File object class: ${myFile.class}"
        println "File name: ${myFile.name}"
        println "Simple name: ${myFile.simpleName}"
        println "Extension: ${myFile.extension}"
        println "Parent directory: ${myFile.parent}"

        // Count the lines in the file
        COUNT_LINES(myFile)
    ```

=== "Avant"

    ```groovy title="main.nf" linenums="7" hl_lines="4-9"
        // Create a Path object from a string path
        myFile = file('data/patientA_rep1_normal_R1_001.fastq.gz')

        // Print file attributes
        println "File object class: ${myFile.class}"
        println "File name: ${myFile.name}"
        println "Simple name: ${myFile.simpleName}"
        println "Extension: ${myFile.extension}"
        println "Parent directory: ${myFile.parent}"
    ```

Et maintenant exécutez le workflow :

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

Spécifiquement, Nextflow a effectué les opérations suivantes avec succès :

- Mise en place du fichier dans le répertoire de travail
- Décompression du fichier .gz
- Comptage des lignes (40 lignes dans ce cas)
- Achèvement sans erreur

La clé de cette opération fluide est que nous disons explicitement à Nextflow que notre entrée est un fichier et doit être traitée comme tel.

### 1.5. Dépanner les erreurs de base d'entrée de fichier

Cela déroute souvent les nouveaux·elles venu·es à Nextflow, alors prenons quelques minutes pour voir ce qui se passe quand vous le faites mal.

Il y a deux endroits principaux où vous pouvez vous tromper dans la gestion des fichiers : au niveau du workflow, et au niveau du processus.

#### 1.5.1. Erreur au niveau du workflow

Voyons ce qui se passe si nous revenons à traiter le fichier comme une chaîne lorsque nous spécifions l'entrée dans le bloc workflow.

Effectuez les modifications suivantes dans le workflow, en vous assurant de commenter les instructions d'affichage spécifiques au chemin :

=== "Après"

    ```groovy title="main.nf" linenums="7" hl_lines="2 6-11"
        // Create a Path object from a string path
        myFile = 'data/patientA_rep1_normal_R1_001.fastq.gz'

        // Print file attributes
        println "File object class: ${myFile.class}"
        /*
        println "File name: ${myFile.name}"
        println "Simple name: ${myFile.simpleName}"
        println "Extension: ${myFile.extension}"
        println "Parent directory: ${myFile.parent}"
        */

        // Count the lines in the file
        COUNT_LINES(myFile)
    ```

=== "Avant"

    ```groovy title="main.nf" linenums="7" hl_lines="4-9"
        // Create a Path object from a string path
        myFile = file('data/patientA_rep1_normal_R1_001.fastq.gz')

        // Print file attributes
        println "File object class: ${myFile.class}"
        println "File name: ${myFile.name}"
        println "Simple name: ${myFile.simpleName}"
        println "Extension: ${myFile.extension}"
        println "Parent directory: ${myFile.parent}"

        // Count the lines in the file
        COUNT_LINES(myFile)
    ```

Et maintenant exécutez le workflow :

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

!!! warning "Conserver l'erreur du workflow de 1.5.1"

    Pour que ce test fonctionne correctement, gardez le workflow dans son état cassé (utilisant une simple chaîne au lieu de `file()`).
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

Et maintenant exécutez à nouveau le workflow :

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

Cela montre beaucoup de détails sur l'erreur car le processus est configuré pour afficher des informations de débogage, comme noté ci-dessus.

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

Cela indique que le système n'a pas pu trouver le fichier ; cependant, si vous vérifiez le chemin, il y a bien un fichier avec ce nom à cet emplacement.

Lorsque nous avons exécuté cela, Nextflow a transmis la valeur de chaîne au script, mais il n'a pas _mis en place_ le fichier réel dans le répertoire de travail.
Donc le processus a essayé d'utiliser la chaîne relative, `data/patientA_rep1_normal_R1_001.fastq.gz`, mais ce fichier n'existe pas dans le répertoire de travail du processus.

Pris ensemble, ces deux exemples vous montrent à quel point il est important de dire à Nextflow si une entrée doit être gérée comme un fichier.

!!! note

    Assurez-vous de revenir en arrière et de corriger les deux erreurs intentionnelles avant de continuer à la section suivante.

### À retenir

- Chaînes de chemin vs objets Path : Les chaînes sont juste du texte, les objets Path sont des références de fichiers intelligentes
- La méthode `file()` convertit un chemin de chaîne en un objet Path avec lequel Nextflow peut travailler
- Vous pouvez accéder aux propriétés de fichier comme `name`, `simpleName`, `extension`, et `parent` [en utilisant les attributs de fichier](https://www.nextflow.io/docs/latest/working-with-files.html#getting-file-attributes)
- L'utilisation d'objets Path au lieu de chaînes permet à Nextflow de gérer correctement les fichiers dans votre workflow
- Résultats d'entrée de processus : Une gestion appropriée des fichiers nécessite des objets Path, pas des chaînes, pour garantir que les fichiers sont correctement mis en place et accessibles pour être utilisés par les processus.

---

## 2. Utilisation de fichiers distants

L'une des fonctionnalités clés de Nextflow est la capacité de basculer de manière transparente entre les fichiers locaux (sur la même machine) et les fichiers distants accessibles via Internet.

Si vous le faites correctement, vous ne devriez jamais avoir besoin de changer la logique de votre workflow pour accommoder des fichiers provenant de différents emplacements.
Tout ce que vous devez faire pour utiliser un fichier distant est de spécifier le préfixe approprié dans le chemin du fichier lorsque vous le fournissez au workflow.

Par exemple, `/path/to/data` n'a pas de préfixe, indiquant qu'il s'agit d'un chemin de fichier local 'normal', alors que `s3://path/to/data` inclut le préfixe `s3://`, indiquant qu'il est situé dans le stockage d'objets S3 d'Amazon.

De nombreux protocoles différents sont pris en charge :

- HTTP(S)/FTP (http://, https://, ftp://)
- Amazon S3 (s3://)
- Azure Blob Storage (az://)
- Google Cloud Storage (gs://)

Pour utiliser l'un de ceux-ci, spécifiez simplement le préfixe pertinent dans la chaîne, qui est alors techniquement appelée un Uniform Resource Identifier (URI) au lieu d'un chemin de fichier.
Nextflow gérera l'authentification et la mise en place des fichiers au bon endroit, le téléchargement ou l'envoi et toutes les autres opérations sur fichiers que vous attendez.

La force clé de ce système est qu'il nous permet de basculer entre les environnements sans changer aucune logique de pipeline.
Par exemple, vous pouvez développer avec un petit ensemble de test local avant de passer à un ensemble de test complet situé dans un stockage distant simplement en changeant l'URI.

### 2.1. Utiliser un fichier depuis Internet

Testons cela en remplaçant le chemin local que nous fournissons à notre workflow par un chemin HTTPS pointant vers une copie des mêmes données stockées dans Github.

!!! warning

    Cela ne fonctionnera que si vous avez une connexion Internet active.

Ouvrez à nouveau `main.nf` et changez le chemin d'entrée comme suit :

=== "Après"

    ```groovy title="main.nf" linenums="2" hl_lines="2"
        // Using a remote file from the internet
        myFile = file('https://raw.github.com/nextflow-io/training/master/side-quests/working_with_files/data/patientA_rep1_normal_R1_001.fastq.gz')

        // Print file attributes
        println "File object class: ${myFile.class}"
        println "File name: ${myFile.name}"
        println "Simple name: ${myFile.simpleName}"
        println "Extension: ${myFile.extension}"
        println "Parent directory: ${myFile.parent}"
    ```

=== "Avant"

    ```groovy title="main.nf" linenums="2" hl_lines="2"
        // Create a Path object from a string path
        myFile = file('data/patientA_rep1_normal_R1_001.fastq.gz')

        // Print file attributes
        println "File object class: ${myFile.class}"
        println "File name: ${myFile.name}"
        println "Simple name: ${myFile.simpleName}"
        println "Extension: ${myFile.extension}"
        println "Parent directory: ${myFile.parent}"
    ```

Exécutons le workflow :

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

La seule différence dans la sortie de la console est que la classe de l'objet chemin est maintenant `nextflow.file.http.XPath`, alors que pour le chemin local la classe était `sun.nio.fs.UnixPath`.
Vous n'avez pas besoin de mémoriser ces classes ; nous mentionnons simplement cela pour démontrer que Nextflow identifie et gère les différents emplacements de manière appropriée.

En coulisses, Nextflow a téléchargé le fichier dans un répertoire de mise en place situé dans le répertoire de travail.
Ce fichier mis en place peut ensuite être traité comme un fichier local et lié symboliquement dans le répertoire de processus pertinent.

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
Donc si vous exécutez à nouveau sur le même fichier et que vous n'avez pas supprimé le fichier mis en place, Nextflow utilisera la copie mise en place.

Cela montre à quel point il est facile de basculer entre les données locales et distantes en utilisant Nextflow, ce qui est une fonctionnalité clé de Nextflow.

!!! note

    La seule exception importante à ce principe est que vous ne pouvez pas utiliser de motifs glob ou de chemins de répertoires avec HTTPS car HTTPS ne peut pas lister plusieurs fichiers, vous devez donc spécifier des URL de fichiers exactes.
    Cependant, d'autres protocoles de stockage tels que le stockage blob (`s3://`, `az://`, `gs://`) peuvent utiliser à la fois des globs et des chemins de répertoires.

    Voici comment vous pourriez utiliser des motifs glob avec le stockage cloud :

    ```groovy title="Cloud storage examples (not runnable in this environment)"
    // S3 with glob patterns - would match multiple files
    ch_s3_files = channel.fromPath('s3://my-bucket/data/*.fastq.gz')

    // Azure Blob Storage with glob patterns
    ch_azure_files = channel.fromPath('az://container/data/patient*_R{1,2}.fastq.gz')

    // Google Cloud Storage with glob patterns
    ch_gcs_files = channel.fromPath('gs://bucket/data/sample_*.fastq.gz')
    ```

    Nous vous montrerons comment travailler avec les globs en pratique dans la section suivante.

### 2.2. Revenir au fichier local

Nous allons revenir à l'utilisation de nos fichiers d'exemple locaux pour le reste de cette quête secondaire, alors revenons à l'entrée du workflow vers le fichier original :

=== "Après"

    ```groovy title="main.nf" linenums="2" hl_lines="2"
        // Create a Path object from a string path
        myFile = file('data/patientA_rep1_normal_R1_001.fastq.gz')

        // Print file attributes
        println "File object class: ${myFile.class}"
        println "File name: ${myFile.name}"
        println "Simple name: ${myFile.simpleName}"
        println "Extension: ${myFile.extension}"
        println "Parent directory: ${myFile.parent}"
    ```

=== "Avant"

    ```groovy title="main.nf" linenums="2" hl_lines="2"
        // Create a Path object from a string path
        myFile = file('https://raw.github.com/nextflow-io/training/master/side-quests/working_with_files/data/patientA_rep1_normal_R1_001.fastq.gz')

        // Print file attributes
        println "File object class: ${myFile.class}"
        println "File name: ${myFile.name}"
        println "Simple name: ${myFile.simpleName}"
        println "Extension: ${myFile.extension}"
        println "Parent directory: ${myFile.parent}"
    ```

### À retenir

- Les données distantes sont accessibles en utilisant un URI (HTTP, FTP, S3, Azure, Google Cloud)
- Nextflow téléchargera et mettra automatiquement en place les données au bon endroit, tant que ces chemins sont fournis aux processus
- N'écrivez pas de logique pour télécharger ou envoyer des fichiers distants !
- Les fichiers locaux et distants produisent différents types d'objets mais fonctionnent de manière identique
- **Important** : HTTP/HTTPS ne fonctionne qu'avec des fichiers uniques (pas de motifs glob)
- Le stockage cloud (S3, Azure, GCS) prend en charge à la fois les fichiers uniques et les motifs glob
- Vous pouvez basculer de manière transparente entre les sources de données locales et distantes sans changer la logique du code (tant que le protocole prend en charge vos opérations requises)

---

## 3. Utilisation du channel factory `fromPath()`

Jusqu'à présent, nous avons travaillé avec un seul fichier à la fois, mais dans Nextflow, nous voulons généralement créer un canal d'entrée avec plusieurs fichiers d'entrée à traiter.

Une façon naïve de le faire serait de combiner la méthode `file()` avec [`channel.of()`](https://www.nextflow.io/docs/latest/reference/channel.html#of) comme ceci :

```groovy title="Syntax example"
ch_files = channel.of([file('data/patientA_rep1_normal_R1_001.fastq.gz')],
                      [file('data/patientA_rep1_normal_R1_001.fastq.gz')])
```

Cela fonctionne, mais c'est maladroit.

!!! tip "Quand utiliser `file()` vs `channel.fromPath()`"

    - Utilisez `file()` lorsque vous avez besoin d'un seul objet Path pour une manipulation directe (vérifier si un fichier existe, lire ses attributs, ou passer à une seule invocation de processus)
    - Utilisez `channel.fromPath()` lorsque vous avez besoin d'un canal pouvant contenir plusieurs fichiers, en particulier avec des motifs glob, ou lorsque les fichiers circuleront à travers plusieurs processus

C'est là qu'intervient [`channel.fromPath()`](https://www.nextflow.io/docs/latest/reference/channel.html#frompath) : un channel factory pratique qui regroupe toutes les fonctionnalités dont nous avons besoin pour générer un canal à partir d'une ou plusieurs chaînes de fichiers statiques ainsi que de motifs glob.

### 3.1. Ajouter le channel factory

Mettons à jour notre workflow pour utiliser `channel.fromPath`.

=== "Après"

    ```groovy title="main.nf" linenums="7" hl_lines="1-3"
        // Load files with channel.fromPath
        ch_files = channel.fromPath('data/patientA_rep1_normal_R1_001.fastq.gz')
        ch_files.view { myFile -> "Found file: $myFile" }

        // Print file attributes
        /* Commentez-les pour l'instant, nous y reviendrons !
        println "File object class: ${myFile.class}"
        println "File name: ${myFile.name}"
        println "Simple name: ${myFile.simpleName}"
        println "Extension: ${myFile.extension}"
        println "Parent directory: ${myFile.parent}"
        */

        // Count the lines in the file
        // COUNT_LINES(myFile)
    ```

=== "Avant"

    ```groovy title="main.nf" linenums="7" hl_lines="1-2"
        // Create a Path object from a string path
        myFile = file('data/patientA_rep1_normal_R1_001.fastq.gz')

        // Print file attributes
        println "File object class: ${myFile.class}"
        println "File name: ${myFile.name}"
        println "Simple name: ${myFile.simpleName}"
        println "Extension: ${myFile.extension}"
        println "Parent directory: ${myFile.parent}"

        // Count the lines in the file
        COUNT_LINES(myFile)
    ```

Nous avons également commenté le code qui affiche les attributs pour l'instant, et ajouté une instruction `.view` pour afficher uniquement le nom de fichier à la place.

Exécutez le workflow :

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
C'est similaire à ce que `file()` aurait fait, sauf que maintenant nous avons un canal dans lequel nous pouvons charger plus de fichiers si nous le voulons.

L'utilisation de `channel.fromPath()` est un moyen pratique de créer un nouveau canal peuplé par une liste de fichiers.

### 3.2. Afficher les attributs des fichiers dans le canal

Dans notre première tentative d'utilisation du channel factory, nous avons simplifié le code et affiché simplement le nom du fichier.

Revenons à l'affichage des attributs complets du fichier :

=== "Après"

    ```groovy title="main.nf" linenums="7" hl_lines="3-9 12"
        // Load files with channel.fromPath
        ch_files = channel.fromPath('data/patientA_rep1_normal_R1_001.fastq.gz')
        ch_files.view { myFile ->
            println "File object class: ${myFile.class}"
            println "File name: ${myFile.name}"
            println "Simple name: ${myFile.simpleName}"
            println "Extension: ${myFile.extension}"
            println "Parent directory: ${myFile.parent}"
        }

        // Count the lines in the file
        COUNT_LINES(ch_files)
    ```

=== "Avant"

    ```groovy title="main.nf" linenums="7" hl_lines="3"
        // Load files with channel.fromPath
        ch_files = channel.fromPath('data/patientA_rep1_normal_R1_001.fastq.gz')
        ch_files.view { myFile -> "Found file: $myFile" }

        // Count the lines in the file
        // COUNT_LINES(ch_files)
    ```

Nous réactivons également l'appel au processus `COUNT_LINES` pour vérifier que le traitement des fichiers fonctionne toujours correctement avec notre approche basée sur les canaux.

Exécutez le workflow :

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

Et voilà, mêmes résultats qu'avant mais maintenant nous avons le fichier dans un canal, donc nous pouvons en ajouter d'autres.

### 3.3. Utiliser un glob pour correspondre à plusieurs fichiers

Il existe plusieurs façons de charger plus de fichiers dans le canal.
Ici, nous allons vous montrer comment utiliser des motifs glob, qui sont un moyen pratique de faire correspondre et récupérer des noms de fichiers et de répertoires basés sur des caractères génériques.
Le processus de correspondance de ces motifs est appelé "globbing" ou "expansion de nom de fichier".

!!! note

    Comme noté précédemment, Nextflow prend en charge le globbing pour gérer les fichiers d'entrée et de sortie dans la majorité des cas, sauf avec les chemins de fichiers HTTPS car HTTPS ne peut pas lister plusieurs fichiers.

Disons que nous voulons récupérer les deux fichiers d'une paire de fichiers associés à un patient donné, `patientA` :

```console
patientA_rep1_normal_R1_001.fastq.gz
patientA_rep1_normal_R2_001.fastq.gz
```

Puisque la seule différence entre les noms de fichiers est le numéro de réplicat, _c'est-à-dire_ le numéro après `R`, nous pouvons utiliser le caractère générique `*` pour remplacer le numéro comme suit :

```console
patientA_rep1_normal_R*_001.fastq.gz
```

C'est le motif glob dont nous avons besoin.

Maintenant, tout ce que nous devons faire est de mettre à jour le chemin du fichier dans le channel factory pour utiliser ce motif glob comme suit :

=== "Après"

    ```groovy title="main.nf" linenums="7"
      // Load files with channel.fromPath
        ch_files = channel.fromPath('data/patientA_rep1_normal_R*_001.fastq.gz')
    ```

=== "Avant"

    ```groovy title="main.nf" linenums="7"
      // Load files with channel.fromPath
        ch_files = channel.fromPath('data/patientA_rep1_normal_R1_001.fastq.gz')
    ```

Nextflow reconnaîtra automatiquement qu'il s'agit d'un motif glob et le gérera de manière appropriée.

Exécutez le workflow pour tester cela :

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

Comme vous pouvez le voir, nous avons maintenant deux objets Path dans notre canal, ce qui montre que Nextflow a effectué l'expansion de nom de fichier correctement, et a chargé et traité les deux fichiers comme prévu.

En utilisant cette méthode, nous pouvons récupérer autant ou aussi peu de fichiers que nous le voulons simplement en changeant le motif glob. Si nous le rendions plus généreux, par exemple en remplaçant toutes les parties variables des noms de fichiers par `*` (_par exemple_ `data/patient*_rep*_*_R*_001.fastq.gz`) nous pourrions récupérer tous les fichiers d'exemple dans le répertoire `data`.

### À retenir

- `channel.fromPath()` crée un canal avec des fichiers correspondant à un motif
- Chaque fichier est émis comme un élément séparé dans le canal
- Nous pouvons utiliser un motif glob pour correspondre à plusieurs fichiers
- Les fichiers sont automatiquement convertis en objets Path avec tous les attributs
- La méthode `.view()` permet l'inspection du contenu du canal

---

## 4. Extraction de métadonnées de base à partir des noms de fichiers

Dans la plupart des domaines scientifiques, il est très courant d'avoir des métadonnées encodées dans les noms des fichiers qui contiennent les données.
Par exemple, en bioinformatique, les fichiers contenant des données de séquençage sont souvent nommés d'une manière qui encode des informations sur l'échantillon, la condition, le réplicat et le numéro de lecture.

Si les noms de fichiers sont construits selon une convention cohérente, vous pouvez extraire ces métadonnées de manière standardisée et les utiliser au cours de votre analyse.
C'est un grand 'si', bien sûr, et vous devriez être très prudent·e chaque fois que vous vous fiez à la structure des noms de fichiers ; mais la réalité est que cette approche est très largement utilisée, alors voyons comment cela se fait dans Nextflow.

Dans le cas de nos données d'exemple, nous savons que les noms de fichiers incluent des métadonnées structurées de manière cohérente.
Par exemple, le nom de fichier `patientA_rep1_normal_R2_001` encode ce qui suit :

- ID patient : `patientA`
- ID réplicat : `rep1`
- type d'échantillon : `normal` (par opposition à `tumor`)
- ensemble de lectures : `R1` (par opposition à `R2`)

Nous allons modifier notre workflow pour récupérer ces informations en trois étapes :

1. Récupérer le `simpleName` du fichier, qui inclut les métadonnées
2. Séparer les métadonnées en utilisant une méthode appelée `tokenize()`
3. Utiliser une map pour organiser les métadonnées

!!! warning

    Vous ne devriez jamais encoder d'informations sensibles dans les noms de fichiers, telles que les noms de patient·es ou d'autres caractéristiques d'identification, car cela peut compromettre la confidentialité des patient·es ou d'autres restrictions de sécurité pertinentes.

### 4.1. Récupérer le `simpleName`

Le `simpleName` est un attribut de fichier qui correspond au nom de fichier dépouillé de son chemin et de son extension.

Effectuez les modifications suivantes dans le workflow :

=== "Après"

    ```groovy title="main.nf" linenums="7" hl_lines="3-6"
        // Load files with channel.fromPath
        ch_files = channel.fromPath('data/patientA_rep1_normal_R*_001.fastq.gz')
        ch_files.map { myFile ->
            [ myFile.simpleName, myFile ]
        }
        .view()
    ```

=== "Avant"

    ```groovy title="main.nf" linenums="7" hl_lines="3-9"
        // Load files with channel.fromPath
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

Exécutez le workflow pour tester que cela fonctionne :

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

Chaque élément dans le canal est maintenant un tuple contenant le `simpleName` et l'objet fichier original.

### 4.2. Extraire les métadonnées du `simplename`

À ce stade, les métadonnées que nous voulons sont intégrées dans le `simplename`, mais nous ne pouvons pas accéder directement aux éléments individuels.
Nous devons donc diviser le `simplename` en ses composants.
Heureusement, ces composants sont simplement séparés par des underscores dans le nom de fichier original, nous pouvons donc appliquer une méthode Nextflow courante appelée `tokenize()` qui est parfaite pour cette tâche.

Effectuez les modifications suivantes dans le workflow :

=== "Après"

    ```groovy title="main.nf" linenums="7" hl_lines="4"
        // Load files with channel.fromPath
        ch_files = channel.fromPath('data/patientA_rep1_normal_R*_001.fastq.gz')
        ch_files.map { myFile ->
            [ myFile.simpleName.tokenize('_'), myFile ]
        }
    ```

=== "Avant"

    ```groovy title="main.nf" linenums="7" hl_lines="4"
        // Load files with channel.fromPath
        ch_files = channel.fromPath('data/patientA_rep1_normal_R*_001.fastq.gz')
        ch_files.map { myFile ->
            [ myFile.simpleName, myFile ]
        }
    ```

La méthode `tokenize()` divisera la chaîne `simpleName` partout où elle trouve des underscores, et retournera une liste contenant les sous-chaînes.

Exécutez le workflow :

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

Maintenant le tuple pour chaque élément dans notre canal contient la liste de métadonnées (_par exemple_ `[patientA, rep1, normal, R1, 001]`) et l'objet fichier original.

C'est génial !
Nous avons décomposé nos informations patient d'une seule chaîne en une liste de chaînes.
Nous pouvons maintenant gérer chaque partie des informations patient séparément.

### 4.3. Utiliser une map pour organiser les métadonnées

Nos métadonnées ne sont qu'une liste plate pour le moment.
C'est assez facile à utiliser mais difficile à lire.

```console
[patientA, rep1, normal, R1, 001]
```

Quel est l'élément à l'index 3 ? Pouvez-vous le dire sans vous référer à l'explication originale de la structure des métadonnées ?

C'est une excellente opportunité d'utiliser un magasin clé-valeur, où chaque élément a un ensemble de clés et leurs valeurs associées, vous pouvez donc facilement vous référer à chaque clé pour obtenir la valeur correspondante.

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

Convertissons maintenant notre liste plate en une map.
Effectuez les modifications suivantes dans le workflow :

=== "Après"

    ```groovy title="main.nf" linenums="7" hl_lines="4-13"
        // Load files with channel.fromPath
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

=== "Avant"

    ```groovy title="main.nf" linenums="7" hl_lines="4"
        // Load files with channel.fromPath
        ch_files = channel.fromPath('data/patientA_rep1_normal_R*_001.fastq.gz')
        ch_files.map { myFile ->
            [ myFile.simpleName.tokenize('_'), myFile ]
        }
    ```

Les changements clés ici sont :

- **Affectation par déstructuration** : `def (patient, replicate, type, readNum) = ...` extrait les valeurs tokenisées dans des variables nommées en une seule ligne
- **Syntaxe littérale de map** : `[id: patient, replicate: ...]` crée une map où chaque clé (comme `id`) est associée à une valeur (comme `patient`)
- **Structure imbriquée** : La liste externe `[..., myFile]` apparie la map de métadonnées avec l'objet fichier original

Nous avons également simplifié quelques-unes des chaînes de métadonnées en utilisant une méthode de remplacement de chaîne appelée `replace()` pour supprimer certains caractères qui sont inutiles (_par exemple_ `replicate.replace('rep', '')` pour ne garder que le numéro des ID de réplicat).

Exécutons à nouveau le workflow :

```bash
nextflow run main.nf
```

??? success "Sortie de la commande"

    ```console hl_lines="7-8"
    N E X T F L O W   ~  version 25.10.2

    Launching `main.nf` [infallible_swartz] DSL2 - revision: 7f4e68c0cb

    executor >  local (2)
    [1b/e7fb27] COUNT_LINES (1) [100%] 2 of 2 ✔
    [[id:patientA, replicate:1, type:normal, readNum:2], /workspaces/training/side-quests/working_with_files/data/patientA_rep1_normal_R2_001.fastq.gz]
    [[id:patientA, replicate:1, type:normal, readNum:1], /workspaces/training/side-quests/working_with_files/data/patientA_rep1_normal_R1_001.fastq.gz]
    Processing file: patientA_rep1_normal_R2_001.fastq.gz
    40

    Processing file: patientA_rep1_normal_R1_001.fastq.gz
    40
    ```

Maintenant les métadonnées sont bien étiquetées (_par exemple_ `[id:patientA, replicate:1, type:normal, readNum:2]`) donc il est beaucoup plus facile de dire ce qui est quoi.

Il sera également beaucoup plus facile d'utiliser réellement des éléments de métadonnées dans le workflow, et rendra notre code plus facile à lire et plus maintenable.

### À retenir

- Nous pouvons gérer les noms de fichiers dans Nextflow avec la puissance d'un langage de programmation complet
- Nous pouvons traiter les noms de fichiers comme des chaînes pour extraire les informations pertinentes
- L'utilisation de méthodes comme `tokenize()` et `replace()` nous permet de manipuler les chaînes dans le nom de fichier
- L'opération `.map()` transforme les éléments du canal tout en préservant la structure
- Les métadonnées structurées (maps) rendent le code plus lisible et maintenable que les listes positionnelles

Ensuite, nous verrons comment gérer les fichiers de données appariés.

---

## 5. Gestion des fichiers de données appariés

De nombreux plans expérimentaux produisent des fichiers de données appariés qui bénéficient d'être gérés de manière explicitement appariée.
Par exemple, en bioinformatique, les données de séquençage sont souvent générées sous forme de lectures appariées, c'est-à-dire des chaînes de séquence qui proviennent du même fragment d'ADN (souvent appelées 'forward' et 'reverse' car elles sont lues à partir d'extrémités opposées).

C'est le cas de nos données d'exemple, où R1 et R2 font référence aux deux ensembles de lectures.

```console
data/patientA_rep1_normal_R1_001.fastq.gz
data/patientA_rep1_normal_R2_001.fastq.gz
```

Nextflow fournit un channel factory spécialisé pour travailler avec des fichiers appariés comme celui-ci appelé `channel.fromFilePairs()`, qui regroupe automatiquement les fichiers en fonction d'un motif de nommage partagé. Cela vous permet d'associer les fichiers appariés plus étroitement avec moins d'effort.

Nous allons modifier notre workflow pour en profiter.
Cela va prendre deux étapes :

1. Changer le channel factory pour `channel.fromFilePairs()`
2. Extraire et mapper les métadonnées

### 5.1. Changer le channel factory pour `channel.fromFilePairs()`

Pour utiliser `channel.fromFilePairs`, nous devons spécifier le motif que Nextflow devrait utiliser pour identifier les deux membres d'une paire.

En revenant à nos données d'exemple, nous pouvons formaliser le motif de nommage comme suit :

```console
data/patientA_rep1_normal_R{1,2}_001.fastq.gz
```

C'est similaire au motif glob que nous avons utilisé précédemment, sauf que cela énumère spécifiquement les sous-chaînes (soit `1` soit `2` venant juste après le R) qui identifient les deux membres de la paire.

Mettons à jour le workflow `main.nf` en conséquence :

=== "Après"

    ```groovy title="main.nf" linenums="7" hl_lines="1-2"
        // Load files with channel.fromFilePairs
        ch_files = channel.fromFilePairs('data/patientA_rep1_normal_R{1,2}_001.fastq.gz')
        /* Commentez le mapping pour l'instant, nous y reviendrons !
        ch_files.map { myFile ->
            def (sample, replicate, type, readNum) = myFile.simpleName.tokenize('_')
            [
                [
                    id: sample,
                    replicate: replicate.replace('rep', ''),
                    type: type,
                    readNum: readNum,
                ],
                myFile
            ]
        }
        */
        .view()
    ```

=== "Avant"

    ```groovy title="main.nf" linenums="7" hl_lines="1-2"
        // Load files with channel.fromFilePairs
        ch_files = channel.fromPath('data/patientA_rep1_normal_R*_001.fastq.gz')
        ch_files.map { myFile ->
            def (sample, replicate, type, readNum) = myFile.simpleName.tokenize('_')
            [
                [
                    id: sample,
                    replicate: replicate.replace('rep', ''),
                    type: type,
                    readNum: readNum,
                ],
                myFile
            ]
        }
        .view()
    ```

Nous avons changé le channel factory et adapté le motif de correspondance de fichiers, et pendant que nous y étions, nous avons commenté l'opération map.
Nous la rajouterons plus tard, avec quelques modifications.

Exécutez le workflow pour le tester :

```bash
nextflow run main.nf
```

??? failure "Sortie de la commande"

    ```console hl_lines="7-8"
     N E X T F L O W   ~  version 25.10.2

    Launching `main.nf` [angry_koch] DSL2 - revision: 44fdf66105

    [-        ] COUNT_LINES -
    [-        ] COUNT_LINES -
    [patientA_rep1_normal_R, [/workspaces/training/side-quests/working_with_files/data/patientA_rep1_normal_R1_001.fastq.gz, /workspaces/training/side-quests/working_with_files/data/patientA_rep1_normal_R2_001.fastq.gz]]
    ERROR ~ Error executing process > 'COUNT_LINES (1)'

    Caused by:
      Not a valid path value: 'patientA_rep1_normal_R'



    Tip: when you have fixed the problem you can continue the execution adding the option `-resume` to the run command line

    -- Check '.nextflow.log' file for details
    ```

Oh-oh, cette fois l'exécution a échoué !

La partie pertinente du message d'erreur est ici :

```console
Not a valid path value: 'patientA_rep1_normal_R'
```

C'est parce que nous avons changé le channel factory.
Jusqu'à présent, le canal d'entrée original ne contenait que les chemins de fichiers.
Toute la manipulation de métadonnées que nous avons faite n'affectait pas réellement le contenu du canal.

Maintenant que nous utilisons le channel factory `.fromFilePairs`, le contenu du canal résultant est différent.
Nous ne voyons qu'un seul élément de canal, composé d'un tuple contenant deux éléments : la partie du `simpleName` partagée par les deux fichiers, qui sert d'identifiant, et un tuple contenant les deux objets fichiers, au format `id, [ file1, file2 ]`.

C'est génial, car Nextflow a fait le travail difficile d'extraire le nom du patient en examinant le préfixe partagé et en l'utilisant comme identifiant de patient.

Cependant, cela casse notre workflow actuel.
Si nous voulions toujours exécuter `COUNT_LINES` de la même manière sans changer le processus, nous devrions appliquer une opération de mapping pour extraire les chemins de fichiers.
Mais nous n'allons pas faire cela, car notre objectif ultime est d'utiliser un processus différent, `ANALYZE_READS`, qui gère les paires de fichiers de manière appropriée.

Alors commentons simplement (ou supprimons) l'appel à `COUNT_LINES` et continuons.

=== "Après"

    ```groovy title="main.nf" linenums="26" hl_lines="2"
        // Count the lines in the file
        // COUNT_LINES(ch_files)
    ```

=== "Avant"

    ```groovy title="main.nf" linenums="26" hl_lines="2"
        // Count the lines in the file
        COUNT_LINES(ch_files)
    ```

Vous pouvez également commenter ou supprimer l'instruction include de `COUNT_LINES`, mais cela n'aura aucun effet fonctionnel.

Maintenant exécutons à nouveau le workflow :

```bash
nextflow run main.nf
```

??? success "Sortie de la commande"

    ```console hl_lines="5"
     N E X T F L O W   ~  version 25.10.2

    Launching `main.nf` [fabulous_davinci] DSL2 - revision: 22b53268dc

    [patientA_rep1_normal_R, [/workspaces/training/side-quests/working_with_files/data/patientA_rep1_normal_R1_001.fastq.gz, /workspaces/training/side-quests/working_with_files/data/patientA_rep1_normal_R2_001.fastq.gz]]
    ```

Yay, cette fois le workflow réussit !

Cependant, nous devons toujours extraire le reste des métadonnées du champ `id`.

### 5.2. Extraire et organiser les métadonnées des paires de fichiers

Notre opération `map` d'avant ne fonctionnera pas car elle ne correspond pas à la structure de données, mais nous pouvons la modifier pour qu'elle fonctionne.

Nous avons déjà accès à l'identifiant réel du patient dans la chaîne que `fromFilePairs()` a utilisée comme identifiant, nous pouvons donc l'utiliser pour extraire les métadonnées sans obtenir le `simpleName` de l'objet Path comme nous l'avons fait auparavant.

Décommentez l'opération map dans le workflow et effectuez les modifications suivantes :

=== "Après"

    ```groovy title="main.nf" linenums="7" hl_lines="3-4 9 11 13"
        // Load files with channel.fromFilePairs
        ch_files = channel.fromFilePairs('data/patientA_rep1_normal_R{1,2}_001.fastq.gz')
        ch_files.map { id, files ->
            def (sample, replicate, type) = id.tokenize('_')
            [
                [
                    id: sample,
                    replicate: replicate.replace('rep', ''),
                    type: type
                ],
                files
            ]
        }
        .view()
    ```

=== "Avant"

    ```groovy title="main.nf" linenums="7" hl_lines="3-5 11 13"
        // Load files with channel.fromFilePairs
        ch_files = channel.fromFilePairs('data/patientA_rep1_normal_R{1,2}_001.fastq.gz')
        /* Commentez le mapping pour l'instant, nous y reviendrons !
        ch_files.map { myFile ->
            def (sample, replicate, type, readNum) = myFile.simpleName.tokenize('_')
            [
                [
                    id: sample,
                    replicate: replicate.replace('rep', ''),
                    type: type,
                    readNum: readNum,
                ],
                myFile
            ]
        }
        */
        .view()
    ```

Cette fois la map commence par `id, files` au lieu de juste `myFile`, et `tokenize()` est appliqué à `id` au lieu de `myFile.simpleName`.

Notez également que nous avons supprimé `readNum` de la ligne `tokenize()` ; toutes les sous-chaînes que nous ne nommons pas spécifiquement (en partant de la gauche) seront silencieusement supprimées.
Nous pouvons faire cela car les fichiers appariés sont maintenant étroitement associés, nous n'avons donc plus besoin de `readNum` dans la map de métadonnées.

Exécutons le workflow :

```bash
nextflow run main.nf
```

??? success "Sortie de la commande"

    ```console

    N E X T F L O W   ~  version 25.10.2

    Launching `main.nf` [prickly_stonebraker] DSL2 - revision: f62ab10a3f

    [[id:patientA, replicate:1, type:normal], [/workspaces/training/side-quests/working_with_files/data/patientA_rep1_normal_R1_001.fastq.gz, /workspaces/training/side-quests/working_with_files/data/patientA_rep1_normal_R2_001.fastq.gz]]
    ```

Et voilà : nous avons la map de métadonnées (`[id:patientA, replicate:1, type:normal]`) en première position du tuple de sortie, suivie du tuple de fichiers appariés, comme prévu.

Bien sûr, cela ne récupérera et traitera que cette paire spécifique de fichiers.
Si vous voulez expérimenter le traitement de plusieurs paires, vous pouvez essayer d'ajouter des wildcards dans le motif d'entrée et voir ce qui se passe.
Par exemple, essayez d'utiliser `data/patientA_rep1_*_R{1,2}_001.fastq.gz`

### À retenir

- [`channel.fromFilePairs()` trouve et apparie automatiquement les fichiers associés](https://www.nextflow.io/docs/latest/reference/channel.html#fromfilepairs)
- Cela simplifie la gestion des lectures paired-end dans votre pipeline
- Les fichiers appariés peuvent être regroupés en tuples `[id, [file1, file2]]`
- L'extraction de métadonnées peut être faite à partir de l'ID de fichier apparié plutôt que des fichiers individuels

---

## 6. Utilisation des opérations sur fichiers dans les processus

Maintenant, rassemblons tout cela dans un processus simple pour renforcer comment utiliser les opérations sur fichiers à l'intérieur d'un processus Nextflow.

Nous vous fournissons un module de processus pré-écrit appelé `ANALYZE_READS` qui prend un tuple de métadonnées et une paire de fichiers d'entrée et les analyse.
Nous pourrions imaginer que cela fait de l'alignement de séquences, ou de l'appel de variants ou toute autre étape qui a du sens pour ce type de données.

Commençons.

### 6.1. Importer le processus et examiner le code

Pour utiliser ce processus dans le workflow, nous devons simplement ajouter une instruction d'inclusion de module avant le bloc workflow.

Effectuez la modification suivante dans le workflow :

=== "Après"

    ```groovy title="main.nf" linenums="1" hl_lines="3"
    #!/usr/bin/env nextflow

    include { ANALYZE_READS } from './modules/analyze_reads.nf'

    workflow {
    ```

=== "Avant"

    ```groovy title="main.nf" linenums="1"
    #!/usr/bin/env nextflow

    workflow {
    ```

Vous pouvez ouvrir le fichier de module pour examiner son code :

```groovy title="modules/analyze_reads.nf - process example" linenums="1"
#!/usr/bin/env nextflow

process ANALYZE_READS {
    tag { meta.id }

    publishDir { "results/${meta.id}" }, mode: 'copy'

    input:
    tuple val(meta), path(files)

    output:
    tuple val(meta.id), path("${meta.id}_stats.txt")

    script:
    """
    echo "Sample metadata: ${meta.id}" > ${meta.id}_stats.txt
    echo "Replicate: ${meta.replicate}" >> ${meta.id}_stats.txt
    echo "Type: ${meta.type}" >> ${meta.id}_stats.txt
    echo "Read 1: ${files[0]}" >> ${meta.id}_stats.txt
    echo "Read 2: ${files[1]}" >> ${meta.id}_stats.txt
    echo "Read 1 size: \$(gunzip -dc ${files[0]} | wc -l | awk '{print \$1/4}') reads" >> ${meta.id}_stats.txt
    echo "Read 2 size: \$(gunzip -dc ${files[1]} | wc -l | awk '{print \$1/4}') reads" >> ${meta.id}_stats.txt
    """
}
```

!!! note

    Les directives `tag` et `publishDir` utilisent la syntaxe de closure (`{ ... }`) au lieu de l'interpolation de chaîne (`"${...}"`).
    C'est parce que ces directives référencent des variables d'entrée (`meta`) qui ne sont pas disponibles avant l'exécution.
    La syntaxe de closure diffère l'évaluation jusqu'à ce que le processus s'exécute réellement.

!!! note

    Nous appelons notre map de métadonnées `meta` par convention.
    Pour une plongée plus profonde dans les meta maps, consultez la quête secondaire [Metadata and meta maps](./metadata.md).

### 6.2. Appeler le processus dans le workflow

Maintenant que le processus est disponible pour le workflow, nous pouvons ajouter un appel au processus `ANALYZE_READS` pour l'exécuter.

Pour l'exécuter sur nos données d'exemple, nous devrons faire deux choses :

1. Donner un nom au canal remappé
2. Ajouter un appel au processus

#### 6.2.1. Nommer le canal d'entrée remappé

Nous avons précédemment appliqué les manipulations de mapping directement au canal d'entrée.
Afin de fournir le contenu remappé au processus `ANALYZE_READS` (et de le faire d'une manière claire et facile à lire) nous voulons créer un nouveau canal nommé `ch_samples`.

Nous pouvons le faire en utilisant l'opérateur [`set`](https://www.nextflow.io/docs/latest/reference/operator.html#set).

Dans le workflow principal, remplacez l'opérateur `.view()` par `.set { ch_samples }`, et ajoutez une ligne testant que nous pouvons nous référer au canal par son nom.

=== "Après"

    ```groovy title="main.nf" linenums="7" hl_lines="14 16-17"
        // Load files with channel.fromFilePairs
        ch_files = channel.fromFilePairs('data/patientA_rep1_normal_R{1,2}_001.fastq.gz')
        ch_files.map { id,  files ->
           def (sample, replicate, type, readNum) = id.tokenize('_')
           [
               [
                   id: sample,
                   replicate: replicate.replace('rep', ''),
                   type: type
               ],
               files
           ]
        }
            .set { ch_samples }

        // Temporary: peek into ch_samples
        ch_samples.view()
    ```

=== "Avant"

    ```groovy title="main.nf" linenums="7" hl_lines="14"
        // Load files with channel.fromFilePairs
        ch_files = channel.fromFilePairs('data/patientA_rep1_normal_R{1,2}_001.fastq.gz')
        ch_files.map { id,  files ->
           def (sample, replicate, type, readNum) = id.tokenize('_')
           [
               [
                   id: sample,
                   replicate: replicate.replace('rep', ''),
                   type: type
               ],
               files
           ]
        }
        .view()
    }
    ```

Exécutons ceci :

```bash
nextflow run main.nf
```

??? success "Sortie de la commande"

    ```console
     N E X T F L O W   ~  version 25.10.2

    Launching `main.nf` [goofy_kirch] DSL2 - revision: 3313283e42

    [[id:patientA, replicate:1, type:normal], [/workspaces/training/side-quests/working_with_files/data/patientA_rep1_normal_R1_001.fastq.gz, /workspaces/training/side-quests/working_with_files/data/patientA_rep1_normal_R2_001.fastq.gz]]
    ```

Cela confirme que nous pouvons maintenant nous référer au canal par son nom.

#### 6.2.2. Appeler le processus sur les données

Maintenant appelons réellement le processus `ANALYZE_READS` sur le canal `ch_samples`.

Dans le workflow principal, effectuez les modifications de code suivantes :

=== "Après"

    ```groovy title="main.nf" linenums="23"
        // Run the analysis
        ANALYZE_READS(ch_samples)
    ```

=== "Avant"

    ```groovy title="main.nf" linenums="23"
        // Temporary: peek into ch_samples
        ch_samples.view()
    ```

Exécutons ceci :

```bash
nextflow run main.nf
```

??? success "Sortie de la commande"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `./main.nf` [shrivelled_cori] DSL2 - revision: b546a31769

    executor >  local (1)
    [b5/110360] process > ANALYZE_READS (patientA) [100%] 1 of 1 ✔
    ```

Ce processus est configuré pour publier ses sorties dans un répertoire `results`, alors jetez-y un œil.

??? abstract "Contenu du répertoire et du fichier"

    ```console
    results
    └── patientA
        └── patientA_stats.txt
    ```

    ```txt title="patientA_stats.txt"
    Sample metadata: patientA
    Replicate: 1
    Type: normal
    Read 1: patientA_rep1_normal_R1_001.fastq.gz
    Read 2: patientA_rep1_normal_R2_001.fastq.gz
    Read 1 size: 10 reads
    Read 2 size: 10 reads
    ```

Le processus a pris nos entrées et créé un nouveau fichier contenant les métadonnées du patient, comme prévu.
Splendide !

### 6.3. Inclure beaucoup plus de patients

Bien sûr, cela ne traite qu'une seule paire de fichiers pour un·e seul·e patient·e, ce qui n'est pas exactement le genre de haut débit que vous espérez obtenir avec Nextflow.
Vous voudrez probablement traiter beaucoup plus de données à la fois.

Rappelez-vous que `channel.fromPath()` accepte un _glob_ en entrée, ce qui signifie qu'il peut accepter n'importe quel nombre de fichiers qui correspondent au motif.
Par conséquent, si nous voulons inclure tous les patients, nous pouvons simplement modifier la chaîne d'entrée pour inclure plus de patients, comme noté en passant précédemment.

Faisons semblant de vouloir être aussi gourmand·es que possible.
Effectuez les modifications suivantes dans le workflow :

=== "Après"

    ```groovy title="main.nf" linenums="7" hl_lines="2"
        // Load files with channel.fromFilePairs
        ch_files = channel.fromFilePairs('data/*_R{1,2}_001.fastq.gz')
    ```

=== "Avant"

    ```groovy title="main.nf" linenums="7" hl_lines="2"
        // Load files with channel.fromFilePairs
        ch_files = channel.fromFilePairs('data/patientA_rep1_normal_R{1,2}_001.fastq.gz')
    ```

Exécutez à nouveau le pipeline :

```bash
nextflow run main.nf
```

??? success "Sortie de la commande"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `./main.nf` [big_stonebraker] DSL2 - revision: f7f9b8a76c

    executor >  local (8)
    [d5/441891] process > ANALYZE_READS (patientC) [100%] 8 of 8 ✔
    ```

Le répertoire results devrait maintenant contenir des résultats pour toutes les données disponibles.

??? abstract "Contenu du répertoire"

    ```console
    results
    ├── patientA
    │   └── patientA_stats.txt
    ├── patientB
    │   └── patientB_stats.txt
    └── patientC
        └── patientC_stats.txt
    ```

Succès ! Nous avons analysé tous les patients en une seule fois ! N'est-ce pas ?

Peut-être pas.
Si vous regardez de plus près, nous avons un problème : nous avons deux réplicats pour patientA, mais un seul fichier de sortie !
Nous écrasons le fichier de sortie à chaque fois.

### 6.4. Rendre les fichiers publiés uniques

Puisque nous avons accès aux métadonnées du patient, nous pouvons les utiliser pour rendre les fichiers publiés uniques en incluant des métadonnées différenciatrices, soit dans la structure de répertoires soit dans les noms de fichiers eux-mêmes.

Effectuez la modification suivante dans le workflow :

=== "Après"

    ```groovy title="modules/analyze_reads.nf" linenums="6"
        publishDir { "results/${meta.type}/${meta.id}/${meta.replicate}" }, mode: 'copy'
    ```

=== "Avant"

    ```groovy title="modules/analyze_reads.nf" linenums="6"
        publishDir { "results/${meta.id}" }, mode: 'copy'
    ```

Ici nous montrons l'option d'utiliser des niveaux de répertoires supplémentaires pour tenir compte des types d'échantillons et des réplicats, mais vous pourriez expérimenter en le faisant au niveau du nom de fichier également.

Maintenant exécutez le pipeline une dernière fois, mais assurez-vous de supprimer d'abord le répertoire results pour vous donner un espace de travail propre :

```bash
rm -r results
nextflow run main.nf
```

??? success "Sortie de la commande"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `./main.nf` [insane_swartz] DSL2 - revision: fff18abe6d

    executor >  local (8)
    [e3/449081] process > ANALYZE_READS (patientC) [100%] 8 of 8 ✔
    ```

Vérifiez maintenant le répertoire results :

??? abstract "Contenu du répertoire"

    ```console
    results/
    ├── normal
    │   ├── patientA
    │   │   ├── 1
    │   │   │   └── patientA_stats.txt
    │   │   └── 2
    │   │       └── patientA_stats.txt
    │   ├── patientB
    │   │   └── 1
    │   │       └── patientB_stats.txt
    │   └── patientC
    │       └── 1
    │           └── patientC_stats.txt
    └── tumor
        ├── patientA
        │   ├── 1
        │   │   └── patientA_stats.txt
        │   └── 2
        │       └── patientA_stats.txt
        ├── patientB
        │   └── 1
        │       └── patientB_stats.txt
        └── patientC
            └── 1
                └── patientC_stats.txt
    ```

Et voilà, toutes nos métadonnées, bien organisées. C'est un succès !

Il y a beaucoup plus que vous pouvez faire une fois que vous avez vos métadonnées chargées dans une map comme celle-ci :

1. Créer des répertoires de sortie organisés basés sur les attributs des patient·es
2. Prendre des décisions dans les processus basées sur les propriétés des patient·es
3. Diviser, joindre et recombiner les données basées sur les valeurs de métadonnées

Ce modèle de garder les métadonnées explicites et attachées aux données (plutôt qu'encodées dans les noms de fichiers) est une bonne pratique fondamentale dans Nextflow qui permet de construire des workflows d'analyse robustes et maintenables.
Vous pouvez en apprendre plus à ce sujet dans la quête secondaire [Metadata and meta maps](./metadata.md).

### À retenir

- La directive `publishDir` peut organiser les sorties basées sur les valeurs de métadonnées
- Les métadonnées dans les tuples permettent une organisation structurée des résultats
- Cette approche crée des workflows maintenables avec une provenance de données claire
- Les processus peuvent prendre des tuples de métadonnées et de fichiers en entrée
- La directive `tag` fournit l'identification du processus dans les journaux d'exécution
- La structure du workflow sépare la création de canaux de l'exécution des processus

---

## Résumé

Dans cette quête secondaire, vous avez appris à travailler avec des fichiers dans Nextflow, des opérations de base aux techniques plus avancées pour gérer des collections de fichiers.

L'application de ces techniques dans votre propre travail vous permettra de construire des workflows plus efficaces et maintenables, en particulier lorsque vous travaillez avec un grand nombre de fichiers avec des conventions de nommage complexes.

### Modèles clés

1.  **Opérations de base sur les fichiers :** Nous avons créé des objets Path avec `file()` et accédé aux attributs de fichiers comme le nom, l'extension et le répertoire parent, en apprenant la différence entre les chaînes et les objets Path.

    - Créer un objet Path avec `file()`

    ```groovy
    myFile = file('path/to/file.txt')
    ```

    - Obtenir les attributs du fichier

    ```groovy
    println myFile.name       // file.txt
    println myFile.baseName   // file
    println myFile.extension  // txt
    println myFile.parent     // path/to
    ```

2.  **Utilisation de fichiers distants** : Nous avons appris à basculer de manière transparente entre les fichiers locaux et distants en utilisant des URI, démontrant la capacité de Nextflow à gérer des fichiers de diverses sources sans changer la logique du workflow.

    - Fichier local

    ```groovy
    myFile = file('path/to/file.txt')
    ```

    - FTP

    ```groovy
    myFile = file('ftp://path/to/file.txt')
    ```

    - HTTPS

    ```groovy
    myFile = file('https://path/to/file.txt')
    ```

    - Amazon S3

    ```groovy
    myFile = file('s3://path/to/file.txt')
    ```

    - Azure Blob Storage

    ```groovy
    myFile = file('az://path/to/file.txt')
    ```

    - Google Cloud Storage

    ```groovy
    myFile = file('gs://path/to/file.txt')
    ```

3.  **Chargement de fichiers en utilisant le channel factory `fromPath()` :** Nous avons créé des canaux à partir de motifs de fichiers avec `channel.fromPath()` et affiché leurs attributs de fichiers, y compris les types d'objets.

    - Créer un canal à partir d'un motif de fichier

    ```groovy
     ch_files = channel.fromPath('data/*.fastq.gz')
    ```

    - Obtenir les attributs du fichier

    ```groovy
     ch_files.view { myFile ->
        println "File object class: ${myFile.class}"
        println "File name: ${myFile.name}"
        println "Simple name: ${myFile.simpleName}"
        println "Extension: ${myFile.extension}"
        println "Parent directory: ${myFile.parent}"
    }
    ```

4.  **Extraction de métadonnées patient à partir des noms de fichiers :** Nous avons utilisé `tokenize()` et `replace()` pour extraire et structurer les métadonnées à partir des noms de fichiers, en les convertissant en maps organisées.

    ```groovy
    def name = file.name.tokenize('_')
    def patientId = name[0]
    def replicate = name[1].replace('rep', '')
    def type = name[2]
    def readNum = name[3].replace('R', '')
    ```

5.  **Simplification avec channel.fromFilePairs :** Nous avons utilisé `channel.fromFilePairs()` pour apparier automatiquement les fichiers associés et extraire les métadonnées des ID de fichiers appariés.

    ```groovy
    ch_pairs = channel.fromFilePairs('data/*_R{1,2}_001.fastq.gz')
    ```

6.  **Utilisation des opérations sur fichiers dans les processus :** Nous avons intégré les opérations sur fichiers dans les processus Nextflow avec une gestion appropriée des entrées, en utilisant `publishDir` pour organiser les sorties basées sur les métadonnées.

    - Associer une meta map avec les entrées du processus

    ```groovy
    ch_files = channel.fromFilePairs('data/patientA_rep1_normal_R{1,2}_001.fastq.gz')
    ch_files.map { id,  files ->
        def (sample, replicate, type, readNum) = id.tokenize('_')
        [
            [
                id: sample,
                replicate: replicate.replace('rep', ''),
                type: type
            ],
             files
        ]
    }
        .set { ch_samples }

    ANALYZE_READS(ch_samples)
    ```

    - Organiser les sorties basées sur les métadonnées

    ```groovy
    publishDir { "results/${meta.type}/${meta.id}/${meta.replicate}" }, mode: 'copy'
    ```

### Ressources supplémentaires

- [Documentation Nextflow : Working with Files](https://www.nextflow.io/docs/latest/working-with-files.html)
- [channel.fromPath](https://www.nextflow.io/docs/latest/reference/channel.html#frompath)
- [channel.fromFilePairs](https://www.nextflow.io/docs/latest/reference/channel.html#fromfilepairs)

---

## Et ensuite ?

Retournez au [menu des Quêtes secondaires](./index.md) ou cliquez sur le bouton en bas à droite de la page pour passer au sujet suivant dans la liste.
