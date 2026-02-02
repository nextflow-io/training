# Modèles de Script Nextflow Essentiels

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Traduction assistée par IA - [en savoir plus et suggérer des améliorations](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

Nextflow est un langage de programmation qui s'exécute sur la Java Virtual Machine. Bien que Nextflow soit construit sur [Groovy](http://groovy-lang.org/) et partage une grande partie de sa syntaxe, Nextflow est plus qu'un simple "Groovy avec des extensions" -- c'est un langage autonome avec une [syntaxe](https://nextflow.io/docs/latest/reference/syntax.html) et une [bibliothèque standard](https://nextflow.io/docs/latest/reference/stdlib.html) entièrement spécifiées.

Vous pouvez écrire beaucoup de code Nextflow sans aller au-delà de la syntaxe de base pour les variables, les maps et les listes. La plupart des tutoriels Nextflow se concentrent sur l'orchestration des workflows (channels, processes et flux de données), et vous pouvez aller étonnamment loin avec seulement cela.

Cependant, lorsque vous devez manipuler des données, analyser des noms de fichiers complexes, implémenter une logique conditionnelle ou construire des workflows de production robustes, il est utile de penser à deux aspects distincts de votre code : le **dataflow** (channels, opérateurs, processes et workflows) et le **scripting** (le code à l'intérieur des closures, des fonctions et des scripts de process). Bien que cette distinction soit quelque peu arbitraire—c'est tout du code Nextflow—elle fournit un modèle mental utile pour comprendre quand vous orchestrez votre pipeline par rapport à quand vous manipulez des données. Maîtriser les deux améliore considérablement votre capacité à écrire des workflows clairs et maintenables.

### Objectifs d'apprentissage

Cette quête secondaire vous emmène dans un voyage pratique des concepts de base aux modèles prêts pour la production.
Nous allons transformer un workflow simple de lecture CSV en un pipeline bioinformatique sophistiqué, en le faisant évoluer étape par étape à travers des défis réalistes :

- **Comprendre les frontières :** Distinguer entre les opérations de dataflow et le scripting, et comprendre comment ils fonctionnent ensemble
- **Manipulation de données :** Extraire, transformer et sélectionner des maps et des collections en utilisant des opérateurs puissants
- **Traitement de chaînes :** Analyser des schémas de nommage de fichiers complexes avec des motifs regex et maîtriser l'interpolation de variables
- **Fonctions réutilisables :** Extraire une logique complexe dans des fonctions nommées pour des workflows plus propres et plus maintenables
- **Logique dynamique :** Construire des processes qui s'adaptent à différents types d'entrées et utiliser des closures pour l'allocation dynamique de ressources
- **Routage conditionnel :** Router intelligemment les échantillons à travers différents processes en fonction de leurs caractéristiques de métadonnées
- **Opérations sûres :** Gérer les données manquantes avec élégance grâce aux opérateurs null-safe et valider les entrées avec des messages d'erreur clairs
- **Gestionnaires basés sur la configuration :** Utiliser les gestionnaires d'événements de workflow pour la journalisation, les notifications et la gestion du cycle de vie

### Prérequis

Avant d'entreprendre cette quête secondaire, vous devriez :

- Avoir complété le tutoriel [Hello Nextflow](../hello_nextflow/README.md) ou un cours équivalent pour débutants.
- Être à l'aise avec l'utilisation des concepts et mécanismes de base de Nextflow (processes, channels, opérateurs, travail avec des fichiers, métadonnées)
- Avoir une familiarité de base avec les constructions de programmation courantes (variables, maps, listes)

Ce tutoriel expliquera les concepts de programmation au fur et à mesure que nous les rencontrerons, vous n'avez donc pas besoin d'une vaste expérience en programmation.
Nous commencerons par des concepts fondamentaux et construirons jusqu'aux modèles avancés.

---

## 0. Commencer

#### Ouvrir le codespace de formation

Si vous ne l'avez pas encore fait, assurez-vous d'ouvrir l'environnement de formation comme décrit dans [Configuration de l'environnement](../envsetup/index.md).

[![Open in GitHub Codespaces](https://github.com/codespaces/badge.svg)](https://codespaces.new/nextflow-io/training?quickstart=1&ref=master)

#### Se déplacer dans le répertoire du projet

Déplaçons-nous dans le répertoire où se trouvent les fichiers pour ce tutoriel.

```bash
cd side-quests/essential_scripting_patterns
```

#### Examiner les matériaux

Vous trouverez un fichier de workflow principal et un répertoire `data` contenant des exemples de fichiers de données.

```console title="Contenu du répertoire"
.
├── collect.nf
├── data
│   ├── samples.csv
│   └── sequences
│       ├── SAMPLE_001_S1_L001_R1_001.fastq
│       ├── SAMPLE_002_S2_L001_R1_001.fastq
│       └── SAMPLE_003_S3_L001_R1_001.fastq
├── main.nf
├── modules
│   ├── fastp.nf
│   ├── generate_report.nf
│   └── trimgalore.nf
└── nextflow.config
```

Notre exemple de CSV contient des informations sur des échantillons biologiques qui nécessitent un traitement différent en fonction de leurs caractéristiques :

```console title="samples.csv"
sample_id,organism,tissue_type,sequencing_depth,file_path,quality_score
SAMPLE_001,human,liver,30000000,data/sequences/SAMPLE_001_S1_L001_R1_001.fastq,38.5
SAMPLE_002,mouse,brain,25000000,data/sequences/SAMPLE_002_S2_L001_R1_001.fastq,35.2
SAMPLE_003,human,kidney,45000000,data/sequences/SAMPLE_003_S3_L001_R1_001.fastq,42.1
```

Nous utiliserons cet ensemble de données réaliste pour explorer des techniques de programmation pratiques que vous rencontrerez dans des workflows bioinformatiques réels.

<!-- TODO: Can we make this more domain-agnostic? -->

<!-- TODO: add an assignment statement? #### Examiner l'assignation -->

#### Liste de vérification de la préparation

Pensez-vous être prêt·e à vous lancer ?

- [ ] Je comprends l'objectif de ce cours et ses prérequis
- [ ] Mon codespace est opérationnel
- [ ] J'ai défini mon répertoire de travail de manière appropriée
<!-- - [ ] Je comprends l'assignation -->

Si vous pouvez cocher toutes les cases, vous êtes prêt·e.

---

## 1. Dataflow vs Scripting : Comprendre les Frontières

### 1.1. Identifier Ce Qui Est Quoi

Lors de l'écriture de workflows Nextflow, il est important de distinguer entre le **dataflow** (comment les données se déplacent à travers les channels et les processes) et le **scripting** (le code qui manipule les données et prend des décisions). Construisons un workflow démontrant comment ils fonctionnent ensemble.

#### 1.1.1. Workflow Nextflow de Base

Commencez par un workflow simple qui lit simplement le fichier CSV (nous l'avons déjà fait pour vous dans `main.nf`) :

```groovy title="main.nf" linenums="1"
workflow {
    ch_samples = channel.fromPath("./data/samples.csv")
        .splitCsv(header: true)
        .view()
}
```

Le bloc `workflow` définit la structure de notre pipeline, tandis que `channel.fromPath()` crée un channel à partir d'un chemin de fichier. L'opérateur `.splitCsv()` traite le fichier CSV et convertit chaque ligne en une structure de données map.

Exécutez ce workflow pour voir les données CSV brutes :

```bash
nextflow run main.nf
```

??? success "Sortie de la commande"

    ```console
    Launching `main.nf` [marvelous_tuckerman] DSL2 - revision: 6113e05c17

    [sample_id:SAMPLE_001, organism:human, tissue_type:liver, sequencing_depth:30000000, file_path:data/sequences/SAMPLE_001_S1_L001_R1_001.fastq, quality_score:38.5]
    [sample_id:SAMPLE_002, organism:mouse, tissue_type:brain, sequencing_depth:25000000, file_path:data/sequences/SAMPLE_002_S2_L001_R1_001.fastq, quality_score:35.2]
    [sample_id:SAMPLE_003, organism:human, tissue_type:kidney, sequencing_depth:45000000, file_path:data/sequences/SAMPLE_003_S3_L001_R1_001.fastq, quality_score:42.1]
    ```

#### 1.1.2. Ajout de l'Opérateur Map

Maintenant, nous allons ajouter du scripting pour transformer les données, en utilisant l'opérateur `.map()` que vous connaissez probablement déjà. Cet opérateur prend une 'closure' où nous pouvons écrire du code pour transformer chaque élément.

!!! note

    Une **closure** est un bloc de code qui peut être passé et exécuté plus tard. Considérez-la comme une fonction que vous définissez en ligne. Les closures sont écrites avec des accolades `{ }` et peuvent prendre des paramètres. Elles sont fondamentales pour le fonctionnement des opérateurs Nextflow et si vous écrivez du Nextflow depuis un moment, vous les utilisez peut-être déjà sans vous en rendre compte !

Voici à quoi ressemble cette opération map :

=== "Après"

    ```groovy title="main.nf" linenums="2" hl_lines="3-6"
        ch_samples = channel.fromPath("./data/samples.csv")
            .splitCsv(header: true)
            .map { row ->
                return row
            }
            .view()
    ```

=== "Avant"

    ```groovy title="main.nf" linenums="2" hl_lines="3"
        ch_samples = channel.fromPath("./data/samples.csv")
            .splitCsv(header: true)
            .view()
    ```

C'est notre première **closure** - une fonction anonyme que vous pouvez passer comme argument (similaire aux lambdas en Python ou aux fonctions fléchées en JavaScript). Les closures sont essentielles pour travailler avec les opérateurs Nextflow.

La closure `{ row -> return row }` prend un paramètre `row` (qui pourrait avoir n'importe quel nom : `item`, `sample`, etc.).

Lorsque l'opérateur `.map()` traite chaque élément du channel, il passe cet élément à votre closure. Ici, `row` contient une ligne CSV à la fois.

Appliquez ce changement et exécutez le workflow :

```bash
nextflow run main.nf
```

Vous verrez la même sortie qu'avant, car nous retournons simplement l'entrée sans modification. Cela confirme que l'opérateur map fonctionne correctement. Maintenant, commençons à transformer les données.

#### 1.1.3. Création d'une Structure de Données Map

Maintenant, nous allons écrire une logique de **scripting** à l'intérieur de notre closure pour transformer chaque ligne de données. C'est ici que nous traitons des éléments de données individuels plutôt que d'orchestrer le flux de données.

=== "Après"

    ```groovy title="main.nf" linenums="2" hl_lines="4-12"
        ch_samples = channel.fromPath("./data/samples.csv")
            .splitCsv(header: true)
            .map { row ->
                // Scripting pour la transformation de données
                def sample_meta = [
                    id: row.sample_id.toLowerCase(),
                    organism: row.organism,
                    tissue: row.tissue_type.replaceAll('_', ' ').toLowerCase(),
                    depth: row.sequencing_depth.toInteger(),
                    quality: row.quality_score.toDouble()
                ]
                return sample_meta
            }
            .view()
    ```

=== "Avant"

    ```groovy title="main.nf" linenums="2" hl_lines="4"
        ch_samples = channel.fromPath("./data/samples.csv")
            .splitCsv(header: true)
            .map { row ->
                return row
            }
            .view()
    ```

La map `sample_meta` est une structure de données clé-valeur (comme les dictionnaires en Python, les objets en JavaScript ou les hashes en Ruby) stockant des informations liées : ID d'échantillon, organisme, type de tissu, profondeur de séquençage et score de qualité.

Nous utilisons des méthodes de manipulation de chaînes comme `.toLowerCase()` et `.replaceAll()` pour nettoyer nos données, et des méthodes de conversion de type comme `.toInteger()` et `.toDouble()` pour convertir les données de chaîne du CSV en types numériques appropriés.

Appliquez ce changement et exécutez le workflow :

```bash
nextflow run main.nf
```

??? success "Sortie de la commande"

    ```console
    [id:sample_001, organism:human, tissue:liver, depth:30000000, quality:38.5]
    [id:sample_002, organism:mouse, tissue:brain, depth:25000000, quality:35.2]
    [id:sample_003, organism:human, tissue:kidney, depth:45000000, quality:42.1]
    ```

#### 1.1.4. Ajout de Logique Conditionnelle

Maintenant, ajoutons plus de scripting - cette fois en utilisant un opérateur ternaire pour prendre des décisions basées sur les valeurs des données.

Effectuez le changement suivant :

=== "Après"

    ```groovy title="main.nf" linenums="2" hl_lines="11-12"
        ch_samples = channel.fromPath("./data/samples.csv")
            .splitCsv(header: true)
            .map { row ->
                def sample_meta = [
                    id: row.sample_id.toLowerCase(),
                    organism: row.organism,
                    tissue: row.tissue_type.replaceAll('_', ' ').toLowerCase(),
                    depth: row.sequencing_depth.toInteger(),
                    quality: row.quality_score.toDouble()
                ]
                def priority = sample_meta.quality > 40 ? 'high' : 'normal'
                return sample_meta + [priority: priority]
            }
            .view()
    ```

=== "Avant"

    ```groovy title="main.nf" linenums="2" hl_lines="11"
        ch_samples = channel.fromPath("./data/samples.csv")
            .splitCsv(header: true)
            .map { row ->
                def sample_meta = [
                    id: row.sample_id.toLowerCase(),
                    organism: row.organism,
                    tissue: row.tissue_type.replaceAll('_', ' ').toLowerCase(),
                    depth: row.sequencing_depth.toInteger(),
                    quality: row.quality_score.toDouble()
                ]
                return sample_meta
            }
            .view()
    ```

L'opérateur ternaire est un raccourci pour une instruction if/else qui suit le modèle `condition ? valeur_si_vrai : valeur_si_faux`. Cette ligne signifie : "Si la qualité est supérieure à 40, utiliser 'high', sinon utiliser 'normal'". Son cousin, l'**opérateur Elvis** (`?:`), fournit des valeurs par défaut lorsque quelque chose est null ou vide - nous explorerons ce modèle plus tard dans ce tutoriel.

L'opérateur d'addition de map `+` crée une **nouvelle map** plutôt que de modifier celle existante. Cette ligne crée une nouvelle map qui contient toutes les paires clé-valeur de `sample_meta` plus la nouvelle clé `priority`.

!!! Note

    Ne modifiez jamais les maps passées dans les closures - créez-en toujours de nouvelles en utilisant `+` (par exemple). Dans Nextflow, les mêmes données transitent souvent simultanément par plusieurs opérations. Modifier une map sur place peut provoquer des effets secondaires imprévisibles lorsque d'autres opérations référencent ce même objet. La création de nouvelles maps garantit que chaque opération dispose de sa propre copie propre.

Exécutez le workflow modifié :

```bash
nextflow run main.nf
```

??? success "Sortie de la commande"

    ```console
    [id:sample_001, organism:human, tissue:liver, depth:30000000, quality:38.5, priority:normal]
    [id:sample_002, organism:mouse, tissue:brain, depth:25000000, quality:35.2, priority:normal]
    [id:sample_003, organism:human, tissue:kidney, depth:45000000, quality:42.1, priority:high]
    ```

Nous avons ajouté avec succès une logique conditionnelle pour enrichir nos métadonnées avec un niveau de priorité basé sur les scores de qualité.

#### 1.1.5. Sélection de Maps avec `.subMap()`

Alors que l'opérateur `+` ajoute des clés à une map, vous devez parfois faire l'inverse - extraire uniquement des clés spécifiques. La méthode `.subMap()` est parfaite pour cela.

Ajoutons une ligne pour créer une version simplifiée de nos métadonnées qui ne contient que les champs d'identification :

=== "Après"

    ```groovy title="main.nf" linenums="2" hl_lines="12-15"
        ch_samples = channel.fromPath("./data/samples.csv")
            .splitCsv(header: true)
            .map { row ->
                // Scripting pour la transformation de données
                def sample_meta = [
                    id: row.sample_id.toLowerCase(),
                    organism: row.organism,
                    tissue: row.tissue_type.replaceAll('_', ' ').toLowerCase(),
                    depth: row.sequencing_depth.toInteger(),
                    quality: row.quality_score.toDouble()
                ]
                def id_only = sample_meta.subMap(['id', 'organism', 'tissue'])
                println "ID fields only: ${id_only}"

                def priority = sample_meta.quality > 40 ? 'high' : 'normal'
                return sample_meta + [priority: priority]
            }
            .view()
    ```

=== "Avant"

    ```groovy title="main.nf" linenums="2" hl_lines="12"
        ch_samples = channel.fromPath("./data/samples.csv")
            .splitCsv(header: true)
            .map { row ->
                // Scripting pour la transformation de données
                def sample_meta = [
                    id: row.sample_id.toLowerCase(),
                    organism: row.organism,
                    tissue: row.tissue_type.replaceAll('_', ' ').toLowerCase(),
                    depth: row.sequencing_depth.toInteger(),
                    quality: row.quality_score.toDouble()
                ]
                def priority = sample_meta.quality > 40 ? 'high' : 'normal'
                return sample_meta + [priority: priority]
            }
            .view()
    ```

Exécutez le workflow modifié :

```bash
nextflow run main.nf
```

??? success "Sortie de la commande"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `main.nf` [peaceful_cori] DSL2 - revision: 4cc4a8340f

    ID fields only: [id:sample_001, organism:human, tissue:liver]
    ID fields only: [id:sample_002, organism:mouse, tissue:brain]
    ID fields only: [id:sample_003, organism:human, tissue:kidney]
    [id:sample_001, organism:human, tissue:liver, depth:30000000, quality:38.5, priority:normal]
    [id:sample_002, organism:mouse, tissue:brain, depth:25000000, quality:35.2, priority:normal]
    [id:sample_003, organism:human, tissue:kidney, depth:45000000, quality:42.1, priority:high]
    ```

Cela montre à la fois les métadonnées complètes affichées par l'opération `view()` et le sous-ensemble extrait que nous avons imprimé avec `println`.

La méthode `.subMap()` prend une liste de clés et retourne une nouvelle map contenant uniquement ces clés. Si une clé n'existe pas dans la map d'origine, elle n'est tout simplement pas incluse dans le résultat.

Ceci est particulièrement utile lorsque vous devez créer différentes versions de métadonnées pour différents processes - certains peuvent avoir besoin de métadonnées complètes tandis que d'autres n'ont besoin que de champs d'identification minimaux.

Maintenant, supprimez ces instructions println pour restaurer votre workflow à son état précédent, car nous n'en avons pas besoin pour la suite.

!!! tip "Résumé des Opérations sur les Maps"

    - **Ajouter des clés** : `map1 + [new_key: value]` - Crée une nouvelle map avec des clés supplémentaires
    - **Extraire des clés** : `map1.subMap(['key1', 'key2'])` - Crée une nouvelle map avec uniquement les clés spécifiées
    - **Les deux opérations créent de nouvelles maps** - Les maps originales restent inchangées

#### 1.1.6. Combiner des Maps et Retourner des Résultats

Jusqu'à présent, nous n'avons retourné que ce que la communauté Nextflow appelle la 'meta map', et nous avons ignoré les fichiers auxquels ces métadonnées se rapportent. Mais si vous écrivez des workflows Nextflow, vous voulez probablement faire quelque chose avec ces fichiers.

Produisons une structure de channel comprenant un tuple de 2 éléments : la map de métadonnées enrichie et le chemin de fichier correspondant. C'est un modèle courant dans Nextflow pour passer des données aux processes.

=== "Après"

    ```groovy title="main.nf" linenums="2" hl_lines="12"
        ch_samples = channel.fromPath("./data/samples.csv")
            .splitCsv(header: true)
            .map { row ->
                def sample_meta = [
                    id: row.sample_id.toLowerCase(),
                    organism: row.organism,
                    tissue: row.tissue_type.replaceAll('_', ' ').toLowerCase(),
                    depth: row.sequencing_depth.toInteger(),
                    quality: row.quality_score.toDouble()
                ]
                def priority = sample_meta.quality > 40 ? 'high' : 'normal'
                return tuple( sample_meta + [priority: priority], file(row.file_path) )
            }
            .view()
    ```

=== "Avant"

    ```groovy title="main.nf" linenums="2" hl_lines="12"
        ch_samples = channel.fromPath("./data/samples.csv")
            .splitCsv(header: true)
            .map { row ->
                def sample_meta = [
                    id: row.sample_id.toLowerCase(),
                    organism: row.organism,
                    tissue: row.tissue_type.replaceAll('_', ' ').toLowerCase(),
                    depth: row.sequencing_depth.toInteger(),
                    quality: row.quality_score.toDouble()
                ]
                def priority = sample_meta.quality > 40 ? 'high' : 'normal'
                return sample_meta + [priority: priority]
            }
            .view()
    ```

Appliquez ce changement et exécutez le workflow :

```bash
nextflow run main.nf
```

??? success "Sortie de la commande"

    ```console
    [[id:sample_001, organism:human, tissue:liver, depth:30000000, quality:38.5, priority:normal], /workspaces/training/side-quests/essential_scripting_patterns/data/sequences/SAMPLE_001_S1_L001_R1_001.fastq]
    [[id:sample_002, organism:mouse, tissue:brain, depth:25000000, quality:35.2, priority:normal], /workspaces/training/side-quests/essential_scripting_patterns/data/sequences/SAMPLE_002_S2_L001_R1_001.fastq]
    [[id:sample_003, organism:human, tissue:kidney, depth:45000000, quality:42.1, priority:high], /workspaces/training/side-quests/essential_scripting_patterns/data/sequences/SAMPLE_003_S3_L001_R1_001.fastq]
    ```

Cette structure de tuple `[meta, file]` est un modèle courant dans Nextflow pour passer à la fois des métadonnées et des fichiers associés aux processes.

!!! note

    **Maps et Métadonnées** : Les maps sont fondamentales pour travailler avec les métadonnées dans Nextflow. Pour une explication plus détaillée du travail avec les maps de métadonnées, consultez la quête secondaire [Travailler avec les métadonnées](./metadata.md).

Notre workflow démontre le modèle de base : les **opérations de dataflow** (`workflow`, `channel.fromPath()`, `.splitCsv()`, `.map()`, `.view()`) orchestrent comment les données se déplacent à travers le pipeline, tandis que le **scripting** (maps `[key: value]`, méthodes de chaînes, conversions de type, opérateurs ternaires) à l'intérieur de la closure `.map()` gère la transformation des éléments de données individuels.

### 1.2. Comprendre les Différents Types : Channel vs List

Jusqu'à présent, nous avons pu distinguer entre les opérations de dataflow et le scripting. Mais qu'en est-il lorsque le même nom de méthode existe dans les deux contextes ?

Un exemple parfait est la méthode `collect`, qui existe à la fois pour les types de channel et les types List dans la bibliothèque standard Nextflow. La méthode `collect()` sur une List transforme chaque élément, tandis que l'opérateur `collect()` sur un channel rassemble toutes les émissions du channel en un channel à élément unique.

Démontrons cela avec quelques données d'exemple, en commençant par nous rappeler ce que fait l'opérateur `collect()` de channel. Consultez `collect.nf` :

```groovy title="collect.nf" linenums="1"
def sample_ids = ['sample_001', 'sample_002', 'sample_003']

// channel.collect() - regroupe plusieurs émissions de channel en une seule
ch_input = channel.fromList(sample_ids)
ch_input.view { sample -> "Individual channel item: ${sample}" }
ch_collected = ch_input.collect()
ch_collected.view { list -> "channel.collect() result: ${list} (${list.size()} items grouped into 1)" }
```

Étapes :

- Définir une List d'ID d'échantillons
- Créer un channel avec `fromList()` qui émet chaque ID d'échantillon séparément
- Imprimer chaque élément avec `view()` au fur et à mesure qu'il transite
- Rassembler tous les éléments en une seule liste avec l'opérateur `collect()` du channel
- Imprimer le résultat collecté (élément unique contenant tous les ID d'échantillons) avec un deuxième `view()`

Nous avons modifié la structure du channel, mais nous n'avons pas modifié les données elles-mêmes.

Exécutez le workflow pour confirmer cela :

```bash
nextflow run collect.nf
```

??? success "Sortie de la commande"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `collect.nf` [loving_mendel] DSL2 - revision: e8d054a46e

    Individual channel item: sample_001
    Individual channel item: sample_002
    Individual channel item: sample_003
    channel.collect() result: [sample_001, sample_002, sample_003] (3 items grouped into 1)
    ```

`view()` retourne une sortie pour chaque émission de channel, nous savons donc que cette sortie unique contient tous les 3 éléments originaux regroupés en une seule liste.

Maintenant, voyons la méthode `collect` sur une List en action. Modifiez `collect.nf` pour appliquer la méthode `collect` de la List à la liste originale d'ID d'échantillons :

=== "Après"

    ```groovy title="main.nf" linenums="1" hl_lines="9-13"
    def sample_ids = ['sample_001', 'sample_002', 'sample_003']

    // channel.collect() - regroupe plusieurs émissions de channel en une seule
    ch_input = channel.fromList(sample_ids)
    ch_input.view { sample -> "Individual channel item: ${sample}" }
    ch_collected = ch_input.collect()
    ch_collected.view { list -> "channel.collect() result: ${list} (${list.size()} items grouped into 1)" }

    // List.collect() - transforme chaque élément, préserve la structure
    def formatted_ids = sample_ids.collect { id ->
        id.toUpperCase().replace('SAMPLE_', 'SPECIMEN_')
    }
    println "List.collect() result: ${formatted_ids} (${sample_ids.size()} items transformed into ${formatted_ids.size()})"
    ```

=== "Avant"

    ```groovy title="main.nf" linenums="1"
    def sample_ids = ['sample_001', 'sample_002', 'sample_003']

    // channel.collect() - regroupe plusieurs émissions de channel en une seule
    ch_input = channel.fromList(sample_ids)
    ch_input.view { sample -> "Individual channel item: ${sample}" }
    ch_collected = ch_input.collect()
    ch_collected.view { list -> "channel.collect() result: ${list} (${list.size()} items grouped into 1)" }
    ```

Dans ce nouvel extrait, nous :

- Définissons une nouvelle variable `formatted_ids` qui utilise la méthode `collect` de la List pour transformer chaque ID d'échantillon dans la liste originale
- Imprimons le résultat en utilisant `println`

Exécutez le workflow modifié :

```bash
nextflow run collect.nf
```

??? success "Sortie de la commande"

    ```console hl_lines="5"
    N E X T F L O W   ~  version 25.10.2

    Launching `collect.nf` [cheeky_stonebraker] DSL2 - revision: 2d5039fb47

    List.collect() result: [SPECIMEN_001, SPECIMEN_002, SPECIMEN_003] (3 items transformed into 3)
    Individual channel item: sample_001
    Individual channel item: sample_002
    Individual channel item: sample_003
    channel.collect() result: [sample_001, sample_002, sample_003] (3 items grouped into 1)
    ```

Cette fois, nous n'avons PAS modifié la structure des données, nous avons toujours 3 éléments dans la liste, mais nous AVONS transformé chaque élément en utilisant la méthode `collect` de la List pour produire une nouvelle liste avec des valeurs modifiées. Ceci est similaire à l'utilisation de l'opérateur `map` sur un channel, mais il opère sur une structure de données List plutôt que sur un channel.

`collect` est un cas extrême que nous utilisons ici pour faire valoir un point. La leçon clé est que lorsque vous écrivez des workflows, distinguez toujours entre les **structures de données** (Lists, Maps, etc.) et les **channels** (constructions de dataflow). Les opérations peuvent partager des noms mais se comportent complètement différemment selon le type sur lequel elles sont appelées.

### 1.3. L'Opérateur Spread (`*.`) - Raccourci pour l'Extraction de Propriétés

Liée à la méthode `collect` de List est l'opérateur spread (`*.`), qui fournit un moyen concis d'extraire des propriétés de collections. C'est essentiellement du sucre syntaxique pour un modèle `collect` courant.

Ajoutons une démonstration à notre fichier `collect.nf` :

=== "Après"

    ```groovy title="collect.nf" linenums="1" hl_lines="15-18"
    def sample_ids = ['sample_001', 'sample_002', 'sample_003']

    // channel.collect() - regroupe plusieurs émissions de channel en une seule
    ch_input = channel.fromList(sample_ids)
    ch_input.view { sample -> "Individual channel item: ${sample}" }
    ch_collected = ch_input.collect()
    ch_collected.view { list -> "channel.collect() result: ${list} (${list.size()} items grouped into 1)" }

    // List.collect() - transforme chaque élément, préserve la structure
    def formatted_ids = sample_ids.collect { id ->
        id.toUpperCase().replace('SAMPLE_', 'SPECIMEN_')
    }
    println "List.collect() result: ${formatted_ids} (${sample_ids.size()} items transformed into ${formatted_ids.size()})"

    // Opérateur Spread - accès concis aux propriétés
    def sample_data = [[id: 's1', quality: 38.5], [id: 's2', quality: 42.1], [id: 's3', quality: 35.2]]
    def all_ids = sample_data*.id
    println "Spread operator result: ${all_ids}"
    ```

=== "Avant"

    ```groovy title="collect.nf" linenums="1"
    def sample_ids = ['sample_001', 'sample_002', 'sample_003']

    // channel.collect() - regroupe plusieurs émissions de channel en une seule
    ch_input = channel.fromList(sample_ids)
    ch_input.view { sample -> "Individual channel item: ${sample}" }
    ch_collected = ch_input.collect()
    ch_collected.view { list -> "channel.collect() result: ${list} (${list.size()} items grouped into 1)" }

    // List.collect() - transforme chaque élément, préserve la structure
    def formatted_ids = sample_ids.collect { id ->
        id.toUpperCase().replace('SAMPLE_', 'SPECIMEN_')
    }
    println "List.collect() result: ${formatted_ids} (${sample_ids.size()} items transformed into ${formatted_ids.size()})"
    ```

Exécutez le workflow mis à jour :

```bash title="Tester l'opérateur spread"
nextflow run collect.nf
```

??? success "Sortie de la commande"

    ```console hl_lines="6"
    N E X T F L O W   ~  version 25.10.2

    Launching `collect.nf` [cranky_galileo] DSL2 - revision: 5f3c8b2a91

    List.collect() result: [SPECIMEN_001, SPECIMEN_002, SPECIMEN_003] (3 items transformed into 3)
    Spread operator result: [s1, s2, s3]
    Individual channel item: sample_001
    Individual channel item: sample_002
    Individual channel item: sample_003
    channel.collect() result: [sample_001, sample_002, sample_003] (3 items grouped into 1)
    ```

L'opérateur spread `*.` est un raccourci pour un modèle collect courant :

```groovy
// Ceux-ci sont équivalents :
def ids = samples*.id
def ids = samples.collect { it.id }

// Fonctionne également avec les appels de méthodes :
def names = files*.getName()
def names = files.collect { it.getName() }
```

L'opérateur spread est particulièrement utile lorsque vous devez extraire une seule propriété d'une liste d'objets - c'est plus lisible que d'écrire la closure `collect` complète.

!!! tip "Quand Utiliser Spread vs Collect"

    - **Utilisez spread (`*.`)** pour un accès simple aux propriétés : `samples*.id`, `files*.name`
    - **Utilisez collect** pour des transformations ou une logique complexe : `samples.collect { it.id.toUpperCase() }`, `samples.collect { [it.id, it.quality > 40] }`

### Récapitulatif

Dans cette section, vous avez appris :

- **Dataflow vs scripting** : Les opérateurs de channel orchestrent la façon dont les données transitent par votre pipeline, tandis que le scripting transforme les éléments de données individuels
- **Comprendre les types** : La même méthode (comme `collect`) peut se comporter différemment selon le type sur lequel elle est appelée (Channel vs List)
- **Le contexte compte** : Soyez toujours conscient de savoir si vous travaillez avec des channels (dataflow) ou des structures de données (scripting)

Comprendre ces frontières est essentiel pour le débogage, la documentation et l'écriture de workflows maintenables.

Ensuite, nous plongerons plus profondément dans les capacités de traitement de chaînes, qui sont essentielles pour gérer des données du monde réel.

---

## 2. Traitement de Chaînes et Génération Dynamique de Scripts

Maîtriser le traitement de chaînes distingue les workflows fragiles des pipelines robustes. Cette section couvre l'analyse de noms de fichiers complexes, la génération dynamique de scripts et l'interpolation de variables.

### 2.1. Correspondance de Motifs et Expressions Régulières

Les fichiers bioinformatiques ont souvent des conventions de nommage complexes encodant des métadonnées. Extrayons cela automatiquement en utilisant la correspondance de motifs avec des expressions régulières.

Nous allons retourner à notre workflow `main.nf` et ajouter une logique de correspondance de motifs pour extraire des informations d'échantillon supplémentaires des noms de fichiers. Les fichiers FASTQ de notre ensemble de données suivent les conventions de nommage de style Illumina avec des noms comme `SAMPLE_001_S1_L001_R1_001.fastq.gz`. Ceux-ci peuvent sembler cryptiques, mais ils encodent en fait des métadonnées utiles comme l'ID d'échantillon, le numéro de voie et la direction de lecture. Nous allons utiliser les capacités regex pour analyser ces noms.

Effectuez le changement suivant à votre workflow `main.nf` existant :

=== "Après"

    ```groovy title="main.nf" linenums="4" hl_lines="10-21"
            .map { row ->
                // Scripting pour la transformation de données
                def sample_meta = [
                    id: row.sample_id.toLowerCase(),
                    organism: row.organism,
                    tissue: row.tissue_type.replaceAll('_', ' ').toLowerCase(),
                    depth: row.sequencing_depth.toInteger(),
                    quality: row.quality_score.toDouble()
                ]
                def fastq_path = file(row.file_path)

                def m = (fastq_path.name =~ /^(.+)_S(\d+)_L(\d{3})_(R[12])_(\d{3})\.fastq(?:\.gz)?$/)
                def file_meta = m ? [
                    sample_num: m[0][2].toInteger(),
                    lane: m[0][3],
                    read: m[0][4],
                    chunk: m[0][5]
                ] : [:]

                def priority = sample_meta.quality > 40 ? 'high' : 'normal'
                return tuple(sample_meta + file_meta + [priority: priority], fastq_path)
            }
    ```

=== "Avant"

    ```groovy title="main.nf" linenums="4" hl_lines="10-11"
            .map { row ->
                // Scripting pour la transformation de données
                def sample_meta = [
                    id: row.sample_id.toLowerCase(),
                    organism: row.organism,
                    tissue: row.tissue_type.replaceAll('_', ' ').toLowerCase(),
                    depth: row.sequencing_depth.toInteger(),
                    quality: row.quality_score.toDouble()
                ]
                def priority = sample_meta.quality > 40 ? 'high' : 'normal'
                return tuple(sample_meta + [priority: priority], file(row.file_path))
            }
    ```

Cela démontre les **concepts clés de traitement de chaînes** :

1. **Littéraux d'expressions régulières** utilisant la syntaxe `~/pattern/` - cela crée un motif regex sans avoir besoin d'échapper les barres obliques inverses
2. **Correspondance de motifs** avec l'opérateur `=~` - cela tente de faire correspondre une chaîne avec un motif regex
3. **Objets Matcher** qui capturent des groupes avec `[0][1]`, `[0][2]`, etc. - `[0]` fait référence à la correspondance entière, `[1]`, `[2]`, etc. font référence aux groupes capturés entre parenthèses

Décomposons le motif regex `^(.+)_S(\d+)_L(\d{3})_(R[12])_(\d{3})\.fastq(?:\.gz)?$` :

| Motif               | Correspond à                                 | Capture                            |
| ------------------- | -------------------------------------------- | ---------------------------------- |
| `^(.+)`             | Nom d'échantillon depuis le début            | Groupe 1 : nom d'échantillon       |
| `_S(\d+)`           | Numéro d'échantillon `_S1`, `_S2`, etc.      | Groupe 2 : numéro d'échantillon    |
| `_L(\d{3})`         | Numéro de voie `_L001`                       | Groupe 3 : voie (3 chiffres)       |
| `_(R[12])`          | Direction de lecture `_R1` ou `_R2`          | Groupe 4 : direction de lecture    |
| `_(\d{3})`          | Numéro de fragment `_001`                    | Groupe 5 : fragment (3 chiffres)   |
| `\.fastq(?:\.gz)?$` | Extension de fichier `.fastq` ou `.fastq.gz` | Non capturé (?: est non-capturant) |

Cela analyse les conventions de nommage de style Illumina pour extraire automatiquement les métadonnées.

Exécutez le workflow modifié :

```bash title="Tester la correspondance de motifs"
nextflow run main.nf
```

??? success "Sortie de la commande"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `main.nf` [clever_pauling] DSL2 - revision: 605d2058b4

    [[id:sample_001, organism:human, tissue:liver, depth:30000000, quality:38.5, sample_num:1, lane:001, read:R1, chunk:001, priority:normal], /workspaces/training/side-quests/essential_scripting_patterns/data/sequences/SAMPLE_001_S1_L001_R1_001.fastq]
    [[id:sample_002, organism:mouse, tissue:brain, depth:25000000, quality:35.2, sample_num:2, lane:001, read:R1, chunk:001, priority:normal], /workspaces/training/side-quests/essential_scripting_patterns/data/sequences/SAMPLE_002_S2_L001_R1_001.fastq]
    [[id:sample_003, organism:human, tissue:kidney, depth:45000000, quality:42.1, sample_num:3, lane:001, read:R1, chunk:001, priority:high], /workspaces/training/side-quests/essential_scripting_patterns/data/sequences/SAMPLE_003_S3_L001_R1_001.fastq]
    ```

Cela montre les métadonnées enrichies à partir des noms de fichiers.

### 2.2. Génération Dynamique de Scripts dans les Processes

Les blocs script de process sont essentiellement des chaînes multi-lignes qui sont passées au shell. Vous pouvez utiliser une **logique conditionnelle** (if/else, opérateurs ternaires) pour générer dynamiquement différentes chaînes de script en fonction des caractéristiques d'entrée. Ceci est essentiel pour gérer des types d'entrée divers—comme les lectures single-end vs paired-end—sans dupliquer les définitions de process.

Ajoutons un process à notre workflow qui démontre ce modèle. Ouvrez `modules/fastp.nf` et jetez un œil :

```groovy title="modules/fastp.nf" linenums="1"
process FASTP {
    container 'community.wave.seqera.io/library/fastp:0.24.0--62c97b06e8447690'

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("*_trimmed*.fastq.gz"), emit: reads

    script:
    """
    fastp \\
        --in1 ${reads[0]} \\
        --in2 ${reads[1]} \\
        --out1 ${meta.id}_trimmed_R1.fastq.gz \\
        --out2 ${meta.id}_trimmed_R2.fastq.gz \\
        --json ${meta.id}.fastp.json \\
        --html ${meta.id}.fastp.html \\
        --thread $task.cpus
    """
}
```

Le process prend des fichiers FASTQ en entrée et exécute l'outil `fastp` pour rogner les adaptateurs et filtrer les lectures de faible qualité. Malheureusement, la personne qui a écrit ce process n'a pas prévu les lectures single-end que nous avons dans notre ensemble de données d'exemple. Ajoutons-le à notre workflow et voyons ce qui se passe :

Tout d'abord, incluez le module à la toute première ligne de votre workflow `main.nf` :

```groovy title="main.nf" linenums="1"
include { FASTP } from './modules/fastp.nf'
```

Ensuite, modifiez le bloc `workflow` pour connecter le channel `ch_samples` au process `FASTP` :

=== "Après"

    ```groovy title="main.nf" linenums="25" hl_lines="27"
    workflow {

        ch_samples = channel.fromPath("./data/samples.csv")
            .splitCsv(header: true)
            .map { row ->
                def sample_meta = [
                    id: row.sample_id.toLowerCase(),
                    organism: row.organism,
                    tissue: row.tissue_type.replaceAll('_', ' ').toLowerCase(),
                    depth: row.sequencing_depth.toInteger(),
                    quality: row.quality_score.toDouble()
                ]
                def fastq_path = file(row.file_path)

                def m = (fastq_path.name =~ /^(.+)_S(\d+)_L(\d{3})_(R[12])_(\d{3})\.fastq(?:\.gz)?$/)
                def file_meta = m ? [
                    sample_num: m[0][2].toInteger(),
                    lane: m[0][3],
                    read: m[0][4],
                    chunk: m[0][5]
                ] : [:]

                def priority = sample_meta.quality > 40 ? 'high' : 'normal'
                return tuple(sample_meta + file_meta + [priority: priority], fastq_path)
            }

        ch_fastp = FASTP(ch_samples)
    }
    ```

=== "Avant"

    ```groovy title="main.nf" linenums="25" hl_lines="26"
    workflow {

        ch_samples = channel.fromPath("./data/samples.csv")
            .splitCsv(header: true)
            .map { row ->
                def sample_meta = [
                    id: row.sample_id.toLowerCase(),
                    organism: row.organism,
                    tissue: row.tissue_type.replaceAll('_', ' ').toLowerCase(),
                    depth: row.sequencing_depth.toInteger(),
                    quality: row.quality_score.toDouble()
                ]
                def fastq_path = file(row.file_path)

                def m = (fastq_path.name =~ /^(.+)_S(\d+)_L(\d{3})_(R[12])_(\d{3})\.fastq(?:\.gz)?$/)
                def file_meta = m ? [
                    sample_num: m[0][2].toInteger(),
                    lane: m[0][3],
                    read: m[0][4],
                    chunk: m[0][5]
                ] : [:]

                def priority = sample_meta.quality > 40 ? 'high' : 'normal'
                return [sample_meta + file_meta + [priority: priority], file(row.file_path)]
            }
            .view()
    }
    ```

Exécutez ce workflow modifié :

```bash
nextflow run main.nf
```

??? failure "Sortie de la commande"

    ```console
    ERROR ~ Error executing process > 'FASTP (3)'

    Caused by:
      Process `FASTP (3)` terminated with an error exit status (255)


    Command executed:

      fastp \
          --in1 SAMPLE_003_S3_L001_R1_001.fastq \
          --in2 null \
          --out1 sample_003_trimmed_R1.fastq.gz \
          --out2 sample_003_trimmed_R2.fastq.gz \
          --json sample_003.fastp.json \
          --html sample_003.fastp.html \
          --thread 2

    Command exit status:
      255

    Command output:
      (empty)
    ```

Vous pouvez voir que le process essaie d'exécuter `fastp` avec une valeur `null` pour le deuxième fichier d'entrée, ce qui provoque son échec. C'est parce que notre ensemble de données contient des lectures single-end, mais le process est codé en dur pour attendre des lectures paired-end (deux fichiers d'entrée à la fois).

Corrigez cela en ajoutant une logique conditionnelle au bloc `script:` du process `FASTP`. Une instruction if/else vérifie le nombre de fichiers de lecture et ajuste la commande en conséquence.

=== "Après"

    ```groovy title="main.nf" linenums="10" hl_lines="3-27"
        script:
        // Détection simple single-end vs paired-end
        def is_single = reads instanceof List ? reads.size() == 1 : true

        if (is_single) {
            def input_file = reads instanceof List ? reads[0] : reads
            """
            fastp \\
                --in1 ${input_file} \\
                --out1 ${meta.id}_trimmed.fastq.gz \\
                --json ${meta.id}.fastp.json \\
                --html ${meta.id}.fastp.html \\
                --thread $task.cpus
            """
        } else {
            """
            fastp \\
                --in1 ${reads[0]} \\
                --in2 ${reads[1]} \\
                --out1 ${meta.id}_trimmed_R1.fastq.gz \\
                --out2 ${meta.id}_trimmed_R2.fastq.gz \\
                --json ${meta.id}.fastp.json \\
                --html ${meta.id}.fastp.html \\
                --thread $task.cpus
            """
        }
    ```

=== "Avant"

    ```groovy title="main.nf" linenums="10" hl_lines="2-11"
            script:
            """
            fastp \\
                --in1 ${reads[0]} \\
                --in2 ${reads[1]} \\
                --out1 ${meta.id}_trimmed_R1.fastq.gz \\
                --out2 ${meta.id}_trimmed_R2.fastq.gz \\
                --json ${meta.id}.fastp.json \\
                --html ${meta.id}.fastp.html \\
                --thread $task.cpus
            """
        }
    ```

Maintenant, le workflow peut gérer gracieusement les lectures single-end et paired-end. La logique conditionnelle vérifie le nombre de fichiers d'entrée et construit la commande appropriée pour `fastp`. Voyons si cela fonctionne :

```bash
nextflow run main.nf
```

??? success "Sortie de la commande"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `main.nf` [adoring_rosalind] DSL2 - revision: 04b1cd93e9

    executor >  local (3)
    [31/a8ad4d] process > FASTP (3) [100%] 3 of 3 ✔
    ```

Ça a l'air bon ! Si nous vérifions les commandes réelles qui ont été exécutées (personnalisez pour votre hash de tâche) :

```console title="Vérifier les commandes exécutées"
cat work/31/a8ad4d95749e685a6d842d3007957f/.command.sh
```

Nous pouvons voir que Nextflow a correctement choisi la bonne commande pour les lectures single-end :

```bash title=".command.sh"
#!/bin/bash -ue
fastp \
    --in1 SAMPLE_003_S3_L001_R1_001.fastq \
    --out1 sample_003_trimmed.fastq.gz \
    --json sample_003.fastp.json \
    --html sample_003.fastp.html \
    --thread 2
```

Un autre usage courant de la logique de script dynamique peut être vu dans [le module Génomique de Nextflow pour la Science](../../nf4science/genomics/02_joint_calling). Dans ce module, le process GATK appelé peut prendre plusieurs fichiers d'entrée, mais chacun doit être préfixé par `-V` pour former une ligne de commande correcte. Le process utilise du scripting pour transformer une collection de fichiers d'entrée (`all_gvcfs`) en arguments de commande corrects :

```groovy title="manipulation de ligne de commande pour GATK" linenums="1" hl_lines="2 5"
    script:
    def gvcfs_line = all_gvcfs.collect { gvcf -> "-V ${gvcf}" }.join(' ')
    """
    gatk GenomicsDBImport \
        ${gvcfs_line} \
        -L ${interval_list} \
        --genomicsdb-workspace-path ${cohort_name}_gdb
    """
```

Ces modèles d'utilisation du scripting dans les blocs script de process sont extrêmement puissants et peuvent être appliqués dans de nombreux scénarios - de la gestion de types d'entrée variables à la construction d'arguments de ligne de commande complexes à partir de collections de fichiers, rendant vos processes vraiment adaptables aux exigences diverses des données du monde réel.

### 2.3. Interpolation de Variables : Variables Nextflow et Shell

Les scripts de process mélangent des variables Nextflow, des variables shell et des substitutions de commandes, chacune avec une syntaxe d'interpolation différente. L'utilisation de la mauvaise syntaxe provoque des erreurs. Explorons-les avec un process qui crée un rapport de traitement.

Jetez un œil au fichier de module `modules/generate_report.nf` :

```groovy title="modules/generate_report.nf" linenums="1"
process GENERATE_REPORT {

    publishDir 'results/reports', mode: 'copy'

    input:
    tuple val(meta), path(reads)

    output:
    path "${meta.id}_report.txt"

    script:
    """
    echo "Processing ${reads}" > ${meta.id}_report.txt
    echo "Sample: ${meta.id}" >> ${meta.id}_report.txt
    """
}
```

Ce process écrit un rapport simple avec l'ID d'échantillon et le nom de fichier. Maintenant, exécutons-le pour voir ce qui se passe lorsque nous devons mélanger différents types de variables.

Incluez le process dans votre `main.nf` et ajoutez-le au workflow :

=== "Après"

    ```groovy title="main.nf" linenums="1" hl_lines="2 30"
    include { FASTP } from './modules/fastp.nf'
    include { GENERATE_REPORT } from './modules/generate_report.nf'

    workflow {
        ch_samples = channel.fromPath("./data/samples.csv")
            .splitCsv(header: true)
            .map { row ->
                def sample_meta = [
                    id: row.sample_id.toLowerCase(),
                    organism: row.organism,
                    tissue: row.tissue_type.replaceAll('_', ' ').toLowerCase(),
                    depth: row.sequencing_depth.toInteger(),
                    quality: row.quality_score.toDouble()
                ]
                def fastq_path = file(row.file_path)

                def m = (fastq_path.name =~ /^(.+)_S(\d+)_L(\d{3})_(R[12])_(\d{3})\.fastq(?:\.gz)?$/)
                def file_meta = m ? [
                    sample_num: m[0][2].toInteger(),
                    lane: m[0][3],
                    read: m[0][4],
                    chunk: m[0][5]
                ] : [:]

                def priority = sample_meta.quality > 40 ? 'high' : 'normal'
                return tuple(sample_meta + file_meta + [priority: priority], fastq_path)
            }

        ch_fastp = FASTP(ch_samples)
        GENERATE_REPORT(ch_samples)
    }
    ```

=== "Avant"

    ```groovy title="main.nf" linenums="1" hl_lines="1 10-29"
    include { FASTP } from './modules/fastp.nf'

    workflow {
        ch_samples = channel.fromPath("./data/samples.csv")
            .splitCsv(header: true)
            .map { row ->
                def sample_meta = [
                    id: row.sample_id.toLowerCase(),
                    organism: row.organism,
                    tissue: row.tissue_type.replaceAll('_', ' ').toLowerCase(),
                    depth: row.sequencing_depth.toInteger(),
                    quality: row.quality_score.toDouble()
                ]
                def fastq_path = file(row.file_path)

                def m = (fastq_path.name =~ /^(.+)_S(\d+)_L(\d{3})_(R[12])_(\d{3})\.fastq(?:\.gz)?$/)
                def file_meta = m ? [
                    sample_num: m[0][2].toInteger(),
                    lane: m[0][3],
                    read: m[0][4],
                    chunk: m[0][5]
                ] : [:]

                def priority = sample_meta.quality > 40 ? 'high' : 'normal'
                return tuple(sample_meta + file_meta + [priority: priority], fastq_path)
            }

        ch_fastp = FASTP(ch_samples)
    }
    ```

Maintenant, exécutez le workflow et vérifiez les rapports générés dans `results/reports/`. Ils devraient contenir des informations de base sur chaque échantillon.

<!-- TODO: add the run command -->

??? success "Sortie de la commande"

    ```console
    <!-- TODO: output -->
    ```

Mais que se passe-t-il si nous voulons ajouter des informations sur quand et où le traitement a eu lieu ? Modifions le process pour utiliser des variables **shell** et un peu de substitution de commande pour inclure l'utilisateur actuel, le nom d'hôte et la date dans le rapport :

=== "Après"

    ```groovy title="modules/generate_report.nf" linenums="10" hl_lines="5-7"
        script:
        """
        echo "Processing ${reads}" > ${meta.id}_report.txt
        echo "Sample: ${meta.id}" >> ${meta.id}_report.txt
        echo "Processed by: ${USER}" >> ${meta.id}_report.txt
        echo "Hostname: $(hostname)" >> ${meta.id}_report.txt
        echo "Date: $(date)" >> ${meta.id}_report.txt
        """
    ```

=== "Avant"

    ```groovy title="modules/generate_report.nf" linenums="10"
        script:
        """
        echo "Processing ${reads}" > ${meta.id}_report.txt
        echo "Sample: ${meta.id}" >> ${meta.id}_report.txt
        """
    ```

Si vous exécutez cela, vous remarquerez une erreur - Nextflow essaie d'interpréter `${USER}` comme une variable Nextflow qui n'existe pas.

??? failure "Sortie de la commande"

    ```console
    Error modules/generate_report.nf:15:27: `USER` is not defined
    │  15 |     echo "Processed by: ${USER}" >> ${meta.id}_report.txt
    ╰     |                           ^^^^

    ERROR ~ Script compilation failed
    ```
