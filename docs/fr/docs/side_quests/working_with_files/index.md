# Traitement des fichiers en entrée

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Traduction assistée par IA - [en savoir plus et suggérer des améliorations](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

Les workflows d'analyse scientifique impliquent souvent le traitement d'un grand nombre de fichiers.
Nextflow fournit des outils puissants pour gérer les fichiers efficacement, vous aidant à organiser et traiter vos données avec un minimum de code.

### Objectifs d'apprentissage

Dans cette quête secondaire, nous allons explorer comment Nextflow gère les fichiers, des opérations de base aux techniques plus avancées pour travailler avec des collections de fichiers.
Vous apprendrez à extraire des métadonnées à partir des noms de fichiers, ce qui est une exigence courante dans les pipelines d'analyse scientifique.

À la fin de cette quête secondaire, vous serez en mesure de :

- Créer des objets Path à partir de chaînes de chemins de fichiers en utilisant la méthode `file()` de Nextflow
- Accéder aux attributs des fichiers tels que le nom, l'extension et le répertoire parent
- Gérer de manière transparente les fichiers locaux et distants en utilisant des URIs
- Utiliser des canaux pour automatiser la gestion des fichiers avec `channel.fromPath()` et `channel.fromFilePairs()`
- Extraire et structurer des métadonnées à partir des noms de fichiers en utilisant la manipulation de chaînes
- Regrouper des fichiers liés en utilisant la correspondance de motifs et les expressions glob
- Intégrer des opérations sur les fichiers dans des processus Nextflow avec une gestion appropriée des entrées
- Organiser les sorties des processus en utilisant des structures de répertoires basées sur les métadonnées

Ces compétences vous aideront à construire des workflows capables de gérer différents types d'entrées de fichiers avec une grande flexibilité.

### Prérequis

Avant de vous lancer dans cette quête secondaire, vous devriez :

- Avoir complété le tutoriel [Hello Nextflow](../../hello_nextflow/) ou un cours équivalent pour débutant·es.
- Être à l'aise avec les concepts et mécanismes de base de Nextflow (processus, canaux, opérateurs)

<!-- I removed the suggestion to do the metamaps SQ first because that works more naturally after -->

---

## 0. Premiers pas

#### Ouvrir l'environnement de formation

Si vous ne l'avez pas encore fait, assurez-vous d'ouvrir l'environnement de formation comme décrit dans la [Configuration de l'environnement](../envsetup/index.md).

[![Open in GitHub Codespaces](https://github.com/codespaces/badge.svg)](https://codespaces.new/nextflow-io/training?quickstart=1&ref=master)

#### Se déplacer dans le répertoire du projet

Déplaçons-nous dans le répertoire où se trouvent les fichiers de ce tutoriel.

```bash
cd side-quests/working_with_files
```

Vous pouvez configurer VSCode pour qu'il se concentre sur ce répertoire :

```bash
code .
```

#### Examiner les fichiers

Vous trouverez un fichier workflow simple appelé `main.nf`, un répertoire `modules` contenant deux fichiers de modules, et un répertoire `data` contenant quelques exemples de fichiers de données.

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

Ce répertoire contient des données de séquençage paired-end provenant de trois patient·es (A, B, C).

Pour chaque patient·e, nous avons des échantillons de type `tumor` (provenant généralement de biopsies tumorales) ou `normal` (prélevés sur des tissus sains ou du sang).
Si vous n'êtes pas familier·ère avec l'analyse du cancer, sachez simplement que cela correspond à un modèle expérimental qui utilise des paires d'échantillons tumeur/normal pour effectuer des analyses contrastives.

Pour le patient A spécifiquement, nous avons deux ensembles de réplicats techniques (répétitions).

Les fichiers de données de séquençage sont nommés selon la convention typique `_R1_` et `_R2_` pour ce que l'on appelle les « lectures directes » et les « lectures inverses ».

_Ne vous inquiétez pas si vous n'êtes pas familier·ère avec ce design expérimental, ce n'est pas essentiel pour comprendre ce tutoriel._

#### Examiner l'exercice

Votre défi est d'écrire un workflow Nextflow qui va :

1. **Charger** les fichiers d'entrée en utilisant les méthodes de gestion de fichiers de Nextflow
2. **Extraire** les métadonnées (identifiant du patient, réplicat, type d'échantillon) à partir de la structure du nom de fichier
3. **Regrouper** les fichiers appariés (R1/R2) en utilisant `channel.fromFilePairs()`
4. **Traiter** les fichiers avec un module d'analyse fourni
5. **Organiser** les sorties dans une structure de répertoires basée sur les métadonnées extraites

#### Liste de vérification

Vous pensez être prêt·e à vous lancer ?

- [ ] Je comprends l'objectif de ce cours et ses prérequis
- [ ] Mon codespace est opérationnel
- [ ] J'ai défini mon répertoire de travail de manière appropriée
- [ ] Je comprends l'exercice

Si vous pouvez cocher toutes les cases, vous êtes prêt·e à commencer.

---

## 1. Opérations de base sur les fichiers

### 1.1. Identifier le type d'un objet avec `.class`

Jetez un œil au fichier workflow `main.nf` :

```groovy title="main.nf" linenums="1"
#!/usr/bin/env nextflow

workflow {

    // Crée un objet Path à partir d'un chemin sous forme de chaîne
    myFile = 'data/patientA_rep1_normal_R1_001.fastq.gz'

    println "${myFile} is of class ${myFile.class}"
}
```

Il s'agit d'un mini-workflow (sans aucun processus) qui fait référence à un seul chemin de fichier dans son workflow, puis l'affiche dans la console, accompagné de sa classe.

??? info "Qu'est-ce que `.class` ?"

    Dans Nextflow, `.class` nous indique le type d'objet avec lequel nous travaillons. C'est comme demander « de quel type est cet objet ? » pour savoir s'il s'agit d'une chaîne, d'un nombre, d'un fichier ou d'autre chose.
    Cela nous aidera à illustrer la différence entre une chaîne simple et un objet Path dans les sections suivantes.

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

Comme vous pouvez le voir, Nextflow a affiché le chemin sous forme de chaîne exactement tel que nous l'avons écrit.

Il s'agit simplement d'une sortie texte ; Nextflow n'a encore rien fait de spécial avec elle.
Nous avons également confirmé que, du point de vue de Nextflow, il s'agit uniquement d'une chaîne (de classe `java.lang.String`).
C'est logique, puisque nous n'avons pas encore indiqué à Nextflow qu'elle correspond à un fichier.

### 1.2. Créer un objet Path avec file()

Nous pouvons indiquer à Nextflow comment gérer les fichiers en créant des [objets Path](https://www.nextflow.io/docs/latest/reference/stdlib-types.html#path) à partir de chaînes de chemins.

Dans notre workflow, nous pouvons convertir la chaîne de chemin `data/patientA_rep1_normal_R1_001.fastq.gz` en objet Path en utilisant la méthode `file()`, qui donne accès aux propriétés et opérations du fichier.

Modifiez `main.nf` pour envelopper la chaîne avec `file()` comme suit :

=== "Après"

    ```groovy title="main.nf" linenums="5" hl_lines="2"
        // Crée un objet Path à partir d'un chemin sous forme de chaîne
        myFile = file('data/patientA_rep1_normal_R1_001.fastq.gz')

        println "${myFile} is of class ${myFile.class}"
    ```

=== "Avant"

    ```groovy title="main.nf" linenums="5" hl_lines="2"
        // Crée un objet Path à partir d'un chemin sous forme de chaîne
        myFile = 'data/patientA_rep1_normal_R1_001.fastq.gz'

        println "${myFile} is of class ${myFile.class}"
    ```

Exécutez maintenant le workflow à nouveau :

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
Le chemin du fichier sera désormais absolu, comme dans `/workspaces/training/side-quests/working_with_files/data/patientA_rep1_normal_R1_001.fastq.gz`.

Notez également que la classe de l'objet Path est `sun.nio.fs.UnixPath` : c'est la façon dont Nextflow représente les fichiers locaux.
Comme nous le verrons plus tard, les fichiers distants auront des noms de classe différents (comme `nextflow.file.http.XPath` pour les fichiers HTTP), mais ils fonctionnent tous exactement de la même manière et peuvent être utilisés de façon identique dans vos workflows.

!!! tip "Astuce"

    **La différence essentielle :**

    - **Chaîne de chemin** : Simple texte que Nextflow traite comme des caractères
    - **Objet Path** : Une référence de fichier intelligente avec laquelle Nextflow peut travailler

    Pensez-y ainsi : une chaîne de chemin, c'est comme écrire une adresse sur du papier, tandis qu'un objet Path, c'est comme avoir l'adresse chargée dans un GPS qui sait comment y naviguer et peut vous donner des détails sur le trajet.

### 1.3. Accéder aux attributs des fichiers

En quoi est-ce utile ? Maintenant que Nextflow comprend que `myFile` est un objet Path et pas seulement une chaîne, nous pouvons accéder aux différents attributs de l'objet Path.

Mettons à jour notre workflow pour afficher les attributs de fichier intégrés :

=== "Après"

    ```groovy title="main.nf" linenums="5" hl_lines="4-9"
        // Crée un objet Path à partir d'un chemin sous forme de chaîne
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
        // Crée un objet Path à partir d'un chemin sous forme de chaîne
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

Vous voyez les différents attributs du fichier affichés dans la console ci-dessus.

### 1.4. Passer le fichier à un processus

La différence entre les chaînes et les objets Path devient critique lorsque vous commencez à construire des workflows réels avec des processus.
Jusqu'à présent, nous avons vérifié que Nextflow traite maintenant notre fichier d'entrée comme un fichier, mais voyons si nous pouvons réellement exécuter quelque chose sur ce fichier dans un processus.

#### 1.4.1. Importer le processus et examiner le code

Nous vous fournissons un module de processus pré-écrit appelé `COUNT_LINES` qui prend un fichier en entrée et compte le nombre de lignes qu'il contient.

Pour utiliser le processus dans le workflow, il vous suffit d'ajouter une instruction include avant le bloc workflow :

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

Comme vous pouvez le voir, il s'agit d'un petit script assez simple qui décompresse le fichier et compte le nombre de lignes qu'il contient.

??? info "Que fait `debug true` ?"

    La directive `debug true` dans la définition du processus amène Nextflow à afficher la sortie de votre script (comme le nombre de lignes « 40 ») directement dans le journal d'exécution.
    Sans cela, vous ne verriez que le statut d'exécution du processus, mais pas la sortie réelle de votre script.

    Pour plus d'informations sur le débogage des processus Nextflow, consultez la quête secondaire [Débogage des workflows Nextflow](debugging.md).

#### 1.4.2. Ajouter un appel à `COUNT_LINES`

Maintenant que le processus est disponible pour le workflow, nous pouvons ajouter un appel au processus `COUNT_LINES` pour l'exécuter sur le fichier d'entrée.

Effectuez les modifications suivantes dans le workflow :

=== "Après"

    ```groovy title="main.nf" linenums="7" hl_lines="11-12"
        // Crée un objet Path à partir d'un chemin sous forme de chaîne
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
        // Crée un objet Path à partir d'un chemin sous forme de chaîne
        myFile = file('data/patientA_rep1_normal_R1_001.fastq.gz')

        // Affiche les attributs du fichier
        println "File object class: ${myFile.class}"
        println "File name: ${myFile.name}"
        println "Simple name: ${myFile.simpleName}"
        println "Extension: ${myFile.extension}"
        println "Parent directory: ${myFile.parent}"
    ```

Exécutez maintenant le workflow :

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

Cela montre que nous sommes en mesure d'opérer correctement sur le fichier à l'intérieur d'un processus.

Plus précisément, Nextflow a effectué avec succès les opérations suivantes :

- Mis en place (staging) le fichier dans le répertoire de travail
- Décompressé le fichier .gz
- Compté les lignes (40 lignes dans ce cas)
- Terminé sans erreur

La clé de cette opération fluide est que nous indiquons explicitement à Nextflow que notre entrée est un fichier et doit être traitée comme tel.

### 1.5. Résoudre les erreurs de base liées aux fichiers en entrée

Cela déroute souvent les nouveaux·elles utilisateur·trices de Nextflow, alors prenons quelques minutes pour examiner ce qui se passe quand on fait une erreur.

Il y a deux endroits principaux où la gestion des fichiers peut être incorrecte : au niveau du workflow, et au niveau du processus.

#### 1.5.1. Erreur au niveau du workflow

Voyons ce qui se passe si nous revenons à traiter le fichier comme une chaîne lorsque nous spécifions l'entrée dans le bloc workflow.

Effectuez les modifications suivantes dans le workflow, en vous assurant de commenter les instructions d'affichage spécifiques aux chemins :

=== "Après"

    ```groovy title="main.nf" linenums="7" hl_lines="2 6-11"
        // Crée un objet Path à partir d'un chemin sous forme de chaîne
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
        // Crée un objet Path à partir d'un chemin sous forme de chaîne
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

Exécutez maintenant le workflow :

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

Lorsque vous spécifiez une entrée `path`, Nextflow vérifie que vous passez de véritables références de fichiers, et non de simples chaînes.
Cette erreur vous indique que `'data/patientA_rep1_normal_R1_001.fastq.gz'` n'est pas une valeur de chemin valide car c'est une chaîne, et non un objet Path.

Nextflow a immédiatement détecté le problème et s'est arrêté avant même de démarrer le processus.

#### 1.5.2. Erreur au niveau du processus

L'autre endroit où nous pourrions oublier de spécifier que nous voulons que Nextflow traite l'entrée comme un fichier se trouve dans la définition du processus.

!!! warning "Avertissement"

    **Conservez l'erreur du workflow de la section 1.5.1**

    Pour que ce test fonctionne correctement, laissez le workflow dans son état défectueux (en utilisant une chaîne simple au lieu de `file()`).
    Combiné avec `val` dans le processus, cela produit l'erreur illustrée ci-dessous.

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

Exécutez maintenant le workflow à nouveau :

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

Cela affiche de nombreux détails sur l'erreur car le processus est configuré pour afficher des informations de débogage, comme indiqué ci-dessus.

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

Cela indique que le système n'a pas pu trouver le fichier ; pourtant, si vous vérifiez le chemin, il existe bien un fichier portant ce nom à cet emplacement.

Lors de l'exécution, Nextflow a transmis la valeur de la chaîne au script, mais n'a pas _mis en place_ (staging) le fichier réel dans le répertoire de travail.
Le processus a donc essayé d'utiliser la chaîne relative, `data/patientA_rep1_normal_R1_001.fastq.gz`, mais ce fichier n'existe pas dans le répertoire de travail du processus.

Ensemble, ces deux exemples montrent à quel point il est important d'indiquer à Nextflow si une entrée doit être traitée comme un fichier.

!!! note "Note"

    Assurez-vous de corriger les deux erreurs intentionnelles avant de passer à la section suivante.

### À retenir

- Chaînes de chemin vs objets Path : les chaînes sont simplement du texte, les objets Path sont des références de fichiers intelligentes
- La méthode `file()` convertit une chaîne de chemin en objet Path avec lequel Nextflow peut travailler
- Vous pouvez accéder aux propriétés des fichiers comme `name`, `simpleName`, `extension` et `parent` [en utilisant les attributs de fichier](https://www.nextflow.io/docs/latest/working-with-files.html#getting-file-attributes)
- L'utilisation d'objets Path plutôt que de chaînes permet à Nextflow de gérer correctement les fichiers dans votre workflow
- Résultats de la gestion des entrées dans les processus : une gestion correcte des fichiers nécessite des objets Path, et non des chaînes, pour s'assurer que les fichiers sont correctement mis en place (staging) et accessibles pour être utilisés par les processus.

---

## 2. Utilisation de fichiers distants

L'une des fonctionnalités clés de Nextflow est la capacité de basculer de manière transparente entre des fichiers locaux (sur la même machine) et des fichiers distants accessibles via internet.

Si vous procédez correctement, vous ne devriez jamais avoir besoin de modifier la logique de votre workflow pour accommoder des fichiers provenant de différents emplacements.
Tout ce que vous avez à faire pour utiliser un fichier distant est de spécifier le préfixe approprié dans le chemin du fichier lorsque vous le fournissez au workflow.

Par exemple, `/path/to/data` n'a pas de préfixe, indiquant qu'il s'agit d'un chemin de fichier local « normal », tandis que `s3://path/to/data` inclut le préfixe `s3://`, indiquant qu'il est situé dans le stockage d'objets S3 d'Amazon.

De nombreux protocoles différents sont pris en charge :

- HTTP(S)/FTP (http://, https://, ftp://)
- Amazon S3 (s3://)
- Azure Blob Storage (az://)
- Google Cloud Storage (gs://)

Pour utiliser l'un d'eux, il suffit de spécifier le préfixe pertinent dans la chaîne, qui est alors techniquement appelée un Uniform Resource Identifier (URI) plutôt qu'un chemin de fichier.
Nextflow gérera l'authentification et la mise en place des fichiers au bon endroit, en téléchargeant ou en téléversant et en effectuant toutes les autres opérations sur les fichiers que vous attendez.

La force principale de ce système est qu'il nous permet de basculer entre les environnements sans modifier la logique du pipeline.
Par exemple, vous pouvez développer avec un petit ensemble de test local avant de passer à un ensemble de test à grande échelle situé dans un stockage distant, simplement en changeant l'URI.

### 2.1. Utiliser un fichier depuis internet

Testons cela en remplaçant le chemin local que nous fournissons à notre workflow par un chemin HTTPS pointant vers une copie des mêmes données stockées sur Github.

!!! warning "Avertissement"

    Cela ne fonctionnera que si vous avez une connexion internet active.

Ouvrez à nouveau `main.nf` et modifiez le chemin d'entrée comme suit :

=== "Après"

    ```groovy title="main.nf" linenums="2" hl_lines="2"
        // Utilisation d'un fichier distant depuis internet
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
        // Crée un objet Path à partir d'un chemin sous forme de chaîne
        myFile = file('data/patientA_rep1_normal_R1_001.fastq.gz')

        // Affiche les attributs du fichier
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

Ça fonctionne ! Vous pouvez voir que très peu de choses ont changé.

La seule différence dans la sortie de la console est que la classe de l'objet Path est maintenant `nextflow.file.http.XPath`, alors que pour le chemin local, la classe était `sun.nio.fs.UnixPath`.
Vous n'avez pas besoin de mémoriser ces classes ; nous les mentionnons simplement pour démontrer que Nextflow identifie et gère les différents emplacements de manière appropriée.

En coulisses, Nextflow a téléchargé le fichier dans un répertoire de mise en place (staging) situé dans le répertoire de travail.
Ce fichier mis en place peut ensuite être traité comme un fichier local et lié symboliquement dans le répertoire du processus concerné.

Vous pouvez vérifier que cela s'est produit en examinant le contenu du répertoire de travail situé à la valeur de hachage du processus.

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

Notez que pour les fichiers plus volumineux, l'étape de téléchargement prendra plus de temps par rapport à l'exécution sur des fichiers locaux.
Cependant, Nextflow vérifie s'il dispose déjà d'une copie mise en place pour éviter les téléchargements inutiles.
Donc, si vous exécutez à nouveau sur le même fichier et que vous n'avez pas supprimé le fichier mis en place, Nextflow utilisera la copie mise en place.

Cela montre à quel point il est facile de basculer entre des données locales et distantes avec Nextflow, ce qui est une fonctionnalité clé de Nextflow.

!!! note "Note"

    L'exception importante à ce principe est que vous ne pouvez pas utiliser des motifs glob ou des chemins de répertoires avec HTTPS car HTTPS ne peut pas lister plusieurs fichiers, vous devez donc spécifier des URLs de fichiers exacts.
    Cependant, d'autres protocoles de stockage tels que le stockage blob (`s3://`, `az://`, `gs://`) peuvent utiliser à la fois des globs et des chemins de répertoires.

    Voici comment vous pourriez utiliser des motifs glob avec le stockage cloud :

    ```groovy title="Cloud storage examples (not runnable in this environment)"
    // S3 avec des motifs glob - correspondrait à plusieurs fichiers
    ch_s3_files = channel.fromPath('s3://my-bucket/data/*.fastq.gz')

    // Azure Blob Storage avec des motifs glob
    ch_azure_files = channel.fromPath('az://container/data/patient*_R{1,2}.fastq.gz')

    // Google Cloud Storage avec des motifs glob
    ch_gcs_files = channel.fromPath('gs://bucket/data/sample_*.fastq.gz')
    ```

    Nous vous montrerons comment travailler avec des globs en pratique dans la section suivante.

### 2.2. Revenir au fichier local

Nous allons revenir à l'utilisation de nos exemples de fichiers locaux pour le reste de cette quête secondaire, alors revenons à l'entrée de fichier d'origine dans le workflow :

=== "Après"

    ```groovy title="main.nf" linenums="2" hl_lines="2"
        // Crée un objet Path à partir d'un chemin sous forme de chaîne
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
        // Crée un objet Path à partir d'un chemin sous forme de chaîne
        myFile = file('https://raw.github.com/nextflow-io/training/master/side-quests/working_with_files/data/patientA_rep1_normal_R1_001.fastq.gz')

        // Affiche les attributs du fichier
        println "File object class: ${myFile.class}"
        println "File name: ${myFile.name}"
        println "Simple name: ${myFile.simpleName}"
        println "Extension: ${myFile.extension}"
        println "Parent directory: ${myFile.parent}"
    ```

### À retenir

- Les données distantes sont accessibles via un URI (HTTP, FTP, S3, Azure, Google Cloud)
- Nextflow téléchargera et mettra en place automatiquement les données au bon endroit, à condition que ces chemins soient transmis aux processus
- N'écrivez pas de logique pour télécharger ou téléverser des fichiers distants !
- Les fichiers locaux et distants produisent des types d'objets différents mais fonctionnent de manière identique
- **Important** : HTTP/HTTPS ne fonctionne qu'avec des fichiers uniques (pas de motifs glob)
- Le stockage cloud (S3, Azure, GCS) prend en charge à la fois les fichiers uniques et les motifs glob
- Vous pouvez basculer de manière transparente entre des sources de données locales et distantes sans modifier la logique du code (à condition que le protocole prenne en charge les opérations requises)

---

## 3. Utilisation de la fabrique de canaux `fromPath()`

Jusqu'à présent, nous avons travaillé avec un seul fichier à la fois, mais dans Nextflow, nous allons généralement vouloir créer un canal d'entrée avec plusieurs fichiers d'entrée à traiter.

Une façon naïve de faire cela serait de combiner la méthode `file()` avec [`channel.of()`](https://www.nextflow.io/docs/latest/reference/channel.html#of) comme ceci :

```groovy title="Syntax example"
ch_files = channel.of([file('data/patientA_rep1_normal_R1_001.fastq.gz')],
                      [file('data/patientA_rep1_normal_R1_001.fastq.gz')])
```

Cela fonctionne, mais c'est maladroit.

!!! tip "Astuce : Quand utiliser `file()` vs `channel.fromPath()`"

    - Utilisez `file()` lorsque vous avez besoin d'un seul objet Path pour une manipulation directe (vérifier si un fichier existe, lire ses attributs, ou le passer à une seule invocation de processus)
    - Utilisez `channel.fromPath()` lorsque vous avez besoin d'un canal pouvant contenir plusieurs fichiers, notamment avec des motifs glob, ou lorsque les fichiers vont traverser plusieurs processus

C'est là qu'intervient [`channel.fromPath()`](https://www.nextflow.io/docs/latest/reference/channel.html#frompath) : une fabrique de canaux pratique qui regroupe toutes les fonctionnalités dont nous avons besoin pour générer un canal à partir d'une ou plusieurs chaînes de fichiers statiques ainsi que des motifs glob.

### 3.1. Ajouter la fabrique de canaux

Mettons à jour notre workflow pour utiliser `channel.fromPath`.

=== "Après"

    ```groovy title="main.nf" linenums="7" hl_lines="1-3"
        // Charge les fichiers avec channel.fromPath
        ch_files = channel.fromPath('data/patientA_rep1_normal_R1_001.fastq.gz')
        ch_files.view { myFile -> "Found file: $myFile" }

        // Affiche les attributs du fichier
        /* Commentons cela pour l'instant, nous y reviendrons !
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
        // Crée un objet Path à partir d'un chemin sous forme de chaîne
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

Nous avons également commenté le code qui affiche les attributs pour l'instant, et ajouté une instruction `.view` pour afficher uniquement le nom du fichier à la place.

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
C'est similaire à ce qu'aurait fait `file()`, sauf que nous avons maintenant un canal dans lequel nous pouvons charger davantage de fichiers si nous le souhaitons.

L'utilisation de `channel.fromPath()` est un moyen pratique de créer un nouveau canal peuplé d'une liste de fichiers.

### 3.2. Afficher les attributs des fichiers dans le canal

Dans notre première utilisation de la fabrique de canaux, nous avons simplifié le code et affiché uniquement le nom du fichier.

Revenons à l'affichage des attributs complets du fichier :

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

Et voilà, les mêmes résultats qu'avant, mais maintenant nous avons le fichier dans un canal, ce qui nous permet d'en ajouter d'autres.

### 3.3. Utiliser un glob pour correspondre à plusieurs fichiers

Il existe plusieurs façons de charger davantage de fichiers dans le canal.
Ici, nous allons vous montrer comment utiliser des motifs glob, qui sont un moyen pratique de faire correspondre et de récupérer des noms de fichiers et de répertoires basés sur des caractères génériques.
Le processus de correspondance de ces motifs est appelé « globbing » ou « expansion de noms de fichiers ».

!!! note "Note"

    Comme indiqué précédemment, Nextflow prend en charge le globbing pour gérer les fichiers d'entrée et de sortie dans la majorité des cas, sauf avec les chemins de fichiers HTTPS car HTTPS ne peut pas lister plusieurs fichiers.

Disons que nous voulons récupérer les deux fichiers d'une paire de fichiers associés à un patient donné, `patientA` :

```console
patientA_rep1_normal_R1_001.fastq.gz
patientA_rep1_normal_R2_001.fastq.gz
```

Puisque la seule différence entre les noms de fichiers est le numéro de réplicat, _c'est-à-dire_ le nombre après `R`, nous pouvons utiliser le caractère générique `*` pour remplacer le nombre comme suit :

```console
patientA_rep1_normal_R*_001.fastq.gz
```

C'est le motif glob dont nous avons besoin.

Maintenant, tout ce que nous avons à faire est de mettre à jour le chemin du fichier dans la fabrique de canaux pour utiliser ce motif glob comme suit :

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

Comme vous pouvez le voir, nous avons maintenant deux objets Path dans notre canal, ce qui montre que Nextflow a effectué l'expansion des noms de fichiers correctement, et a chargé et traité les deux fichiers comme prévu.

En utilisant cette méthode, nous pouvons récupérer autant ou aussi peu de fichiers que nous le souhaitons simplement en modifiant le motif glob. Si nous le rendions plus général, par exemple en remplaçant toutes les parties variables des noms de fichiers par `*` (_par ex._ `data/patient*_rep*_*_R*_001.fastq.gz`), nous pourrions récupérer tous les exemples de fichiers dans le répertoire `data`.

### À retenir

- `channel.fromPath()` crée un canal avec les fichiers correspondant à un motif
- Chaque fichier est émis comme un élément séparé dans le canal
- Nous pouvons utiliser un motif glob pour correspondre à plusieurs fichiers
- Les fichiers sont automatiquement convertis en objets Path avec tous leurs attributs
- La méthode `.view()` permet d'inspecter le contenu du canal

---

## 4. Extraction de métadonnées de base à partir des noms de fichiers

Dans la plupart des domaines scientifiques, il est très courant d'avoir des métadonnées encodées dans les noms des fichiers contenant les données.
Par exemple, en bioinformatique, les fichiers contenant des données de séquençage sont souvent nommés de manière à encoder des informations sur l'échantillon, la condition, le réplicat et le numéro de lecture.

Si les noms de fichiers sont construits selon une convention cohérente, vous pouvez extraire ces métadonnées de manière standardisée et les utiliser dans le cadre de votre analyse.
C'est un grand « si », bien sûr, et vous devriez être très prudent·e chaque fois que vous vous fiez à la structure des noms de fichiers ; mais la réalité est que cette approche est très largement utilisée, alors voyons comment cela se fait dans Nextflow.

Dans le cas de nos données d'exemple, nous savons que les noms de fichiers incluent des métadonnées structurées de manière cohérente.
Par exemple, le nom de fichier `patientA_rep1_normal_R2_001` encode les informations suivantes :

- identifiant du patient : `patientA`
- identifiant du réplicat : `rep1`
- type d'échantillon : `normal` (par opposition à `tumor`)
- ensemble de lectures : `R1` (par opposition à `R2`)

Nous allons modifier notre workflow pour récupérer ces informations en trois étapes :

1. Récupérer le `simpleName` du fichier, qui inclut les métadonnées
2. Séparer les métadonnées en utilisant une méthode appelée `tokenize()`
3. Utiliser une map pour organiser les métadonnées

!!! warning "Avertissement"

    Vous ne devez jamais encoder d'informations sensibles dans les noms de fichiers, comme les noms de patients ou d'autres caractéristiques d'identification, car cela peut compromettre la confidentialité des patients ou d'autres restrictions de sécurité pertinentes.

### 4.1. Récupérer le `simpleName`

Le `simpleName` est un attribut de fichier qui correspond au nom de fichier dépouillé de son chemin et de son extension.

Effectuez les modifications suivantes dans le workflow :

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

Chaque élément du canal est maintenant un tuple contenant le `simpleName` et l'objet fichier d'origine.

### 4.2. Extraire les métadonnées du `simpleName`

À ce stade, les métadonnées que nous voulons sont intégrées dans le `simpleName`, mais nous ne pouvons pas accéder directement aux éléments individuels.
Nous devons donc diviser le `simpleName` en ses composants.
Heureusement, ces composants sont simplement séparés par des underscores dans le nom de fichier d'origine, nous pouvons donc appliquer une méthode Nextflow courante appelée `tokenize()` qui est parfaite pour cette tâche.

Effectuez les modifications suivantes dans le workflow :

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

La méthode `tokenize()` va diviser la chaîne `simpleName` partout où elle trouve des underscores, et retournera une liste contenant les sous-chaînes.

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

Maintenant, le tuple pour chaque élément de notre canal contient la liste des métadonnées (_par ex._ `[patientA, rep1, normal, R1, 001]`) et l'objet fichier d'origine.

Excellent !
Nous avons décomposé les informations de notre patient d'une seule chaîne en une liste de chaînes.
Nous pouvons maintenant gérer chaque partie des informations du patient séparément.

### 4.3. Utiliser une map pour organiser les métadonnées

Nos métadonnées ne sont pour l'instant qu'une liste plate.
Elle est assez facile à utiliser mais difficile à lire.

```console
[patientA, rep1, normal, R1, 001]
```

Quel est l'élément à l'index 3 ? Pouvez-vous le dire sans vous référer à l'explication originale de la structure des métadonnées ?

C'est une excellente occasion d'utiliser un stockage clé-valeur, où chaque élément a un ensemble de clés et leurs valeurs associées, afin que vous puissiez facilement vous référer à chaque clé pour obtenir la valeur correspondante.

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
Effectuez les modifications suivantes dans le workflow :

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

=== "Avant"

    ```groovy title="main.nf" linenums="7" hl_lines="4"
        // Charge les fichiers avec channel.fromPath
        ch_files = channel.fromPath('data/patientA_rep1_normal_R*_001.fastq.gz')
        ch_files.map { myFile ->
            [ myFile.simpleName.tokenize('_'), myFile ]
        }
    ```

Les changements clés ici sont :

- **Affectation par déstructuration** : `def (patient, replicate, type, readNum) = ...` extrait les valeurs tokenisées dans des variables nommées en une seule ligne
- **Syntaxe littérale de map** : `[id: patient, replicate: ...]` crée une map où chaque clé (comme `id`) est associée à une valeur (comme `patient`)
- **Structure imbriquée** : La liste externe `[..., myFile]` associe la map de métadonnées à l'objet fichier d'origine

Nous avons également simplifié quelques-unes des chaînes de métadonnées en utilisant une méthode de remplacement de chaîne appelée `replace()` pour supprimer certains caractères inutiles (_par ex._ `replicate.replace('rep', '')` pour ne conserver que le numéro des identifiants de réplicat).

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

Maintenant, les métadonnées sont clairement étiquetées (_par ex._ `[id:patientA, replicate:1, type:normal, readNum:2]`), ce qui rend beaucoup plus facile de savoir ce qui est quoi.

Il sera également beaucoup plus facile d'utiliser réellement les éléments des métadonnées dans le workflow, et cela rendra notre code plus lisible et plus facile à maintenir.

### À retenir

- Nous pouvons gérer les noms de fichiers dans Nextflow avec la puissance d'un langage de programmation complet
- Nous pouvons traiter les noms de fichiers comme des chaînes pour en extraire des informations pertinentes
- L'utilisation de méthodes comme `tokenize()` et `replace()` nous permet de manipuler les chaînes dans le nom de fichier
- L'opération `.map()` transforme les éléments du canal tout en préservant la structure
- Les métadonnées structurées (maps) rendent le code plus lisible et plus facile à maintenir que les listes positionnelles

Ensuite, nous allons voir comment gérer les fichiers de données appariés.

---

## 5. Gestion des fichiers de données appariés

De nombreux designs expérimentaux produisent des fichiers de données appariés qui bénéficient d'être gérés de manière explicitement appariée.
Par exemple, en bioinformatique, les données de séquençage sont souvent générées sous forme de lectures appariées, c'est-à-dire des séquences de chaînes provenant du même fragment d'ADN (souvent appelées « directes » et « inverses » car elles sont lues depuis des extrémités opposées).

C'est le cas de nos données d'exemple, où R1 et R2 font référence aux deux ensembles de lectures.

```console
data/patientA_rep1_normal_R1_001.fastq.gz
data/patientA_rep1_normal_R2_001.fastq.gz
```

Nextflow fournit une fabrique de canaux spécialisée pour travailler avec des fichiers appariés comme ceux-ci, appelée `channel.fromFilePairs()`, qui regroupe automatiquement les fichiers en fonction d'un motif de nommage partagé. Cela vous permet d'associer les fichiers appariés plus étroitement avec moins d'effort.

Nous allons modifier notre workflow pour tirer parti de cela.
Cela va prendre deux étapes :

1. Changer la fabrique de canaux pour `channel.fromFilePairs()`
2. Extraire et mapper les métadonnées

### 5.1. Changer la fabrique de canaux pour `channel.fromFilePairs()`

Pour utiliser `channel.fromFilePairs`, nous devons spécifier le motif que Nextflow doit utiliser pour identifier les deux membres d'une paire.

En revenant à nos données d'exemple, nous pouvons formaliser le motif de nommage comme suit :

```console
data/patientA_rep1_normal_R{1,2}_001.fastq.gz
```

C'est similaire au motif glob que nous avons utilisé précédemment, sauf que celui-ci énumère spécifiquement les sous-chaînes (soit `1` soit `2` venant juste après le R) qui identifient les deux membres de la paire.

Mettons à jour le workflow `main.nf` en conséquence :

=== "Après"

    ```groovy title="main.nf" linenums="7" hl_lines="1-2"
        // Charge les fichiers avec channel.fromFilePairs
        ch_files = channel.fromFilePairs('data/patientA_rep1_normal_R{1,2}_001.fastq.gz')
        /* Commentons le mapping pour l'instant, nous y reviendrons !
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
        // Charge les fichiers avec channel.fromFilePairs
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

Nous avons changé la fabrique de canaux et adapté le motif de correspondance des fichiers, et pendant que nous y étions, nous avons commenté l'opération map.
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

Oups, cette fois l'exécution a échoué !

La partie pertinente du message d'erreur est ici :

```console
Not a valid path value: 'patientA_rep1_normal_R'
```

C'est parce que nous avons changé la fabrique de canaux.
Jusqu'à présent, le canal d'entrée d'origine ne contenait que les chemins de fichiers.
Toutes les manipulations de métadonnées que nous avons effectuées n'ont pas réellement affecté le contenu du canal.

Maintenant que nous utilisons la fabrique de canaux `.fromFilePairs`, le contenu du canal résultant est différent.
Nous ne voyons qu'un seul élément de canal, composé d'un tuple contenant deux éléments : la partie du `simpleName` partagée par les deux fichiers, qui sert d'identifiant, et un tuple contenant les deux objets fichiers, dans le format `id, [ file1, file2 ]`.

C'est excellent, car Nextflow a fait le travail difficile d'extraire le nom du patient en examinant le préfixe partagé et en l'utilisant comme identifiant du patient.

Cependant, cela casse notre workflow actuel.
Si nous voulions toujours exécuter `COUNT_LINES` de la même manière sans modifier le processus, nous devrions appliquer une opération de mapping pour extraire les chemins de fichiers.
Mais nous n'allons pas faire cela, car notre objectif final est d'utiliser un processus différent, `ANALYZE_READS`, qui gère les paires de fichiers de manière appropriée.

Alors commentons simplement (ou supprimons) l'appel à `COUNT_LINES` et passons à la suite.

=== "Après"

    ```groovy title="main.nf" linenums="26" hl_lines="2"
        // Compte les lignes dans le fichier
        // COUNT_LINES(ch_files)
    ```

=== "Avant"

    ```groovy title="main.nf" linenums="26" hl_lines="2"
        // Compte les lignes dans le fichier
        COUNT_LINES(ch_files)
    ```

Vous pouvez également commenter ou supprimer l'instruction include de `COUNT_LINES`, mais cela n'aura aucun effet fonctionnel.

Exécutons maintenant le workflow à nouveau :

```bash
nextflow run main.nf
```

??? success "Sortie de la commande"

    ```console hl_lines="5"
     N E X T F L O W   ~  version 25.10.2

    Launching `main.nf` [fabulous_davinci] DSL2 - revision: 22b53268dc

    [patientA_rep1_normal_R, [/workspaces/training/side-quests/working_with_files/data/patientA_rep1_normal_R1_001.fastq.gz, /workspaces/training/side-quests/working_with_files/data/patientA_rep1_normal_R2_001.fastq.gz]]
    ```

Super, cette fois le workflow réussit !

Cependant, nous devons encore extraire le reste des métadonnées du champ `id`.

### 5.2. Extraire et organiser les métadonnées à partir des paires de fichiers

Notre opération `map` précédente ne fonctionnera pas car elle ne correspond pas à la structure des données, mais nous pouvons la modifier pour qu'elle fonctionne.

Nous avons déjà accès à l'identifiant réel du patient dans la chaîne que `fromFilePairs()` a utilisée comme identifiant, nous pouvons donc l'utiliser pour extraire les métadonnées sans obtenir le `simpleName` de l'objet Path comme nous l'avons fait précédemment.

Décommentez l'opération map dans le workflow et effectuez les modifications suivantes :

=== "Après"

    ```groovy title="main.nf" linenums="7" hl_lines="3-4 9 11 13"
        // Charge les fichiers avec channel.fromFilePairs
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
        // Charge les fichiers avec channel.fromFilePairs
        ch_files = channel.fromFilePairs('data/patientA_rep1_normal_R{1,2}_001.fastq.gz')
        /* Commentons le mapping pour l'instant, nous y reviendrons !
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

Cette fois, la map commence par `id, files` au lieu de simplement `myFile`, et `tokenize()` est appliqué à `id` au lieu de `myFile.simpleName`.

Notez également que nous avons supprimé `readNum` de la ligne `tokenize()` ; toutes les sous-chaînes que nous ne nommons pas spécifiquement (en partant de la gauche) seront silencieusement ignorées.
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

Bien sûr, cela ne récupérera et ne traitera que cette paire spécifique de fichiers.
Si vous souhaitez expérimenter avec le traitement de plusieurs paires, vous pouvez essayer d'ajouter des caractères génériques dans le motif d'entrée et voir ce qui se passe.
Par exemple, essayez d'utiliser `data/patientA_rep1_*_R{1,2}_001.fastq.gz`

### À retenir

- [`channel.fromFilePairs()` trouve et apparie automatiquement les fichiers liés](https://www.nextflow.io/docs/latest/reference/channel.html#fromfilepairs)
- Cela simplifie la gestion des lectures paired-end dans votre pipeline
- Les fichiers appariés peuvent être regroupés sous forme de tuples `[id, [file1, file2]]`
- L'extraction des métadonnées peut être effectuée à partir de l'identifiant des fichiers appariés plutôt que des fichiers individuels

---

## 6. Utilisation des opérations sur les fichiers dans les processus

Maintenant, mettons tout cela ensemble dans un processus simple pour renforcer la façon d'utiliser les opérations sur les fichiers à l'intérieur d'un processus Nextflow.

Nous vous fournissons un module de processus pré-écrit appelé `ANALYZE_READS` qui prend un tuple de métadonnées et une paire de fichiers d'entrée et les analyse.
Nous pourrions imaginer que cela effectue un alignement de séquences, un appel de variants ou toute autre étape qui a du sens pour ce type de données.

Commençons.

### 6.1. Importer le processus et examiner le code

Pour utiliser ce processus dans le workflow, nous avons juste besoin d'ajouter une instruction include de module avant le bloc workflow.

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

!!! note "Note"

    Les directives `tag` et `publishDir` utilisent la syntaxe de closure (`{ ... }`) au lieu de l'interpolation de chaîne (`"${...}"`).
    C'est parce que ces directives font référence à des variables d'entrée (`meta`) qui ne sont pas disponibles avant l'exécution.
    La syntaxe de closure diffère l'évaluation jusqu'à ce que le processus s'exécute réellement.

!!! note "Note"

    Nous appelons notre map de métadonnées `meta` par convention.
    Pour une exploration plus approfondie des meta maps, consultez la quête secondaire [Métadonnées et meta maps](../metadata/).

### 6.2. Appeler le processus dans le workflow

Maintenant que le processus est disponible pour le workflow, nous pouvons ajouter un appel au processus `ANALYZE_READS` pour l'exécuter.

Pour l'exécuter sur nos données d'exemple, nous devrons faire deux choses :

1. Donner un nom au canal remappé
2. Ajouter un appel au processus

#### 6.2.1. Nommer le canal d'entrée remappé

Nous avons précédemment appliqué les manipulations de mapping directement au canal d'entrée.
Afin d'alimenter le contenu remappé au processus `ANALYZE_READS` (et de le faire d'une manière claire et facile à lire), nous voulons créer un nouveau canal nommé `ch_samples`.

Nous pouvons faire cela en utilisant l'opérateur [`set`](https://www.nextflow.io/docs/latest/reference/operator.html#set).

Dans le workflow principal, remplacez l'opérateur `.view()` par `.set { ch_samples }`, et ajoutez une ligne testant que nous pouvons faire référence au canal par son nom.

=== "Après"

    ```groovy title="main.nf" linenums="7" hl_lines="14 16-17"
        // Charge les fichiers avec channel.fromFilePairs
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

        // Temporaire : aperçu de ch_samples
        ch_samples.view()
    ```

=== "Avant"

    ```groovy title="main.nf" linenums="7" hl_lines="14"
        // Charge les fichiers avec channel.fromFilePairs
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

Cela confirme que nous pouvons maintenant faire référence au canal par son nom.

#### 6.2.2. Appeler le processus sur les données

Appelons maintenant réellement le processus `ANALYZE_READS` sur le canal `ch_samples`.

Dans le workflow principal, effectuez les modifications de code suivantes :

=== "Après"

    ```groovy title="main.nf" linenums="23"
        // Exécute l'analyse
        ANALYZE_READS(ch_samples)
    ```

=== "Avant"

    ```groovy title="main.nf" linenums="23"
        // Temporaire : aperçu de ch_samples
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

Ce processus est configuré pour publier ses sorties dans un répertoire `results`, alors jetez un œil là-dedans.

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

### 6.3. Inclure beaucoup plus de patient·es

Bien sûr, cela ne traite qu'une seule paire de fichiers pour un·e seul·e patient·e, ce qui n'est pas exactement le type de débit élevé que vous espérez obtenir avec Nextflow.
Vous voudrez probablement traiter beaucoup plus de données à la fois.

Rappelons que `channel.fromPath()` accepte un _glob_ en entrée, ce qui signifie qu'il peut accepter n'importe quel nombre de fichiers correspondant au motif.
Par conséquent, si nous voulons inclure tous les patient·es, nous pouvons simplement modifier la chaîne d'entrée pour inclure davantage de patient·es, comme mentionné en passant précédemment.

Faisons comme si nous voulions être aussi exhaustif·ves que possible.
Effectuez les modifications suivantes dans le workflow :

=== "Après"

    ```groovy title="main.nf" linenums="7" hl_lines="2"
        // Charge les fichiers avec channel.fromFilePairs
        ch_files = channel.fromFilePairs('data/*_R{1,2}_001.fastq.gz')
    ```

=== "Avant"

    ```groovy title="main.nf" linenums="7" hl_lines="2"
        // Charge les fichiers avec channel.fromFilePairs
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

Le répertoire de résultats devrait maintenant contenir des résultats pour toutes les données disponibles.

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

Succès ! Nous avons analysé tous les patient·es en une seule fois ! N'est-ce pas ?

Peut-être pas.
Si vous regardez de plus près, nous avons un problème : nous avons deux réplicats pour patientA, mais un seul fichier de sortie !
Nous écrasons le fichier de sortie à chaque fois.

### 6.4. Rendre les fichiers publiés uniques

Puisque nous avons accès aux métadonnées du patient, nous pouvons les utiliser pour rendre les fichiers publiés uniques en incluant des métadonnées différenciatrices, soit dans la structure de répertoires, soit dans les noms de fichiers eux-mêmes.

Effectuez la modification suivante dans le workflow :

=== "Après"

    ```groovy title="modules/analyze_reads.nf" linenums="6"
        publishDir { "results/${meta.type}/${meta.id}/${meta.replicate}" }, mode: 'copy'
    ```

=== "Avant"

    ```groovy title="modules/analyze_reads.nf" linenums="6"
        publishDir { "results/${meta.id}" }, mode: 'copy'
    ```

Ici, nous montrons l'option d'utiliser des niveaux de répertoires supplémentaires pour tenir compte des types d'échantillons et des réplicats, mais vous pourriez également expérimenter en le faisant au niveau des noms de fichiers.

Exécutez maintenant le pipeline une fois de plus, mais assurez-vous de supprimer d'abord le répertoire de résultats pour vous donner un espace de travail propre :

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

Vérifiez maintenant le répertoire de résultats :

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

Et voilà, toutes nos métadonnées, soigneusement organisées. C'est un succès !

Il y a beaucoup plus que vous pouvez faire une fois que vous avez vos métadonnées chargées dans une map comme celle-ci :

1. Créer des répertoires de sortie organisés basés sur les attributs du patient
2. Prendre des décisions dans les processus basées sur les propriétés du patient
3. Diviser, joindre et recombiner des données basées sur les valeurs des métadonnées

Ce modèle consistant à garder les métadonnées explicites et attachées aux données (plutôt qu'encodées dans les noms de fichiers) est une bonne pratique fondamentale dans Nextflow qui permet de construire des workflows d'analyse robustes et faciles à maintenir.
Vous pouvez en apprendre davantage à ce sujet dans la quête secondaire [Métadonnées et meta maps](../metadata/).

### À retenir

- La directive `publishDir` peut organiser les sorties en fonction des valeurs des métadonnées
- Les métadonnées dans les tuples permettent une organisation structurée des résultats
- Cette approche crée des workflows faciles à maintenir avec une provenance des données claire
- Les processus peuvent prendre des tuples de métadonnées et de fichiers en entrée
- La directive `tag` fournit l'identification du processus dans les journaux d'exécution
- La structure du workflow sépare la création des canaux de l'exécution des processus

---

## Résumé

Dans cette quête secondaire, vous avez appris à travailler avec des fichiers dans Nextflow, des opérations de base aux techniques plus avancées pour gérer des collections de fichiers.

L'application de ces techniques dans votre propre travail vous permettra de construire des workflows plus efficaces et plus faciles à maintenir, en particulier lorsque vous travaillez avec un grand nombre de fichiers ayant des conventions de nommage complexes.

### Modèles clés

1.  **Opérations de base sur les fichiers :** Nous avons créé des objets Path avec `file()` et accédé aux attributs des fichiers comme le nom, l'extension et le répertoire parent, en apprenant la différence entre les chaînes et les objets Path.

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

2.  **Utilisation de fichiers distants** : Nous avons appris à basculer de manière transparente entre des fichiers locaux et distants en utilisant des URIs, démontrant la capacité de Nextflow à gérer des fichiers provenant de diverses sources sans modifier la logique du workflow.

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

3.  **Chargement de fichiers avec la fabrique de canaux `fromPath()` :** Nous avons créé des canaux à partir de motifs de fichiers avec `channel.fromPath()` et affiché leurs attributs de fichier, y compris les types d'objets.

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

4.  **Extraction des métadonnées du patient à partir des noms de fichiers :** Nous avons utilisé `tokenize()` et `replace()` pour extraire et structurer les métadonnées à partir des noms de fichiers, en les convertissant en maps organisées.

    ```groovy
    def name = file.name.tokenize('_')
    def patientId = name[0]
    def replicate = name[1].replace('rep', '')
    def type = name[2]
    def readNum = name[3].replace('R', '')
    ```

5.  **Simplification avec channel.fromFilePairs :** Nous avons utilisé `channel.fromFilePairs()` pour apparier automatiquement les fichiers liés et extraire les métadonnées à partir des identifiants des fichiers appariés.

    ```groovy
    ch_pairs = channel.fromFilePairs('data/*_R{1,2}_001.fastq.gz')
    ```

6.  **Utilisation des opérations sur les fichiers dans les processus :** Nous avons intégré des opérations sur les fichiers dans des processus Nextflow avec une gestion appropriée des entrées, en utilisant `publishDir` pour organiser les sorties en fonction des métadonnées.

    - Associer une meta map aux entrées du processus

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

    - Organiser les sorties en fonction des métadonnées

    ```groovy
    publishDir { "results/${meta.type}/${meta.id}/${meta.replicate}" }, mode: 'copy'
    ```

### Ressources supplémentaires

- [Documentation Nextflow : Travailler avec des fichiers](https://www.nextflow.io/docs/latest/working-with-files.html)
- [channel.fromPath](https://www.nextflow.io/docs/latest/reference/channel.html#frompath)
- [channel.fromFilePairs](https://www.nextflow.io/docs/latest/reference/channel.html#fromfilepairs)

---

## Et ensuite ?

Retournez au [menu des Quêtes secondaires](../) ou cliquez sur le bouton en bas à droite de la page pour passer au sujet suivant dans la liste.
