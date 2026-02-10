# Modèles de Script Nextflow Essentiels

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Traduction assistée par IA - [en savoir plus et suggérer des améliorations](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

Nextflow est un langage de programmation qui s'exécute sur la machine virtuelle Java. Bien que Nextflow soit construit sur [Groovy](http://groovy-lang.org/) et partage une grande partie de sa syntaxe, Nextflow est plus qu'un simple « Groovy avec des extensions » -- c'est un langage autonome avec une [syntaxe](https://nextflow.io/docs/latest/reference/syntax.html) et une [bibliothèque standard](https://nextflow.io/docs/latest/reference/stdlib.html) entièrement spécifiées.

Vous pouvez écrire beaucoup de code Nextflow sans vous aventurer au-delà de la syntaxe de base pour les variables, les maps et les listes. La plupart des tutoriels Nextflow se concentrent sur l'orchestration des workflows (canaux, processus et flux de données), et vous pouvez aller étonnamment loin avec seulement cela.

Cependant, lorsque vous devez manipuler des données, analyser des noms de fichiers complexes, implémenter une logique conditionnelle ou construire des workflows de production robustes, il est utile de penser à deux aspects distincts de votre code : le **dataflow** (canaux, opérateurs, processus et workflows) et le **scripting** (le code à l'intérieur des closures, fonctions et scripts de processus). Bien que cette distinction soit quelque peu arbitraire -- c'est tout du code Nextflow -- elle fournit un modèle mental utile pour comprendre quand vous orchestrez votre pipeline et quand vous manipulez des données. Maîtriser les deux améliore considérablement votre capacité à écrire des workflows clairs et maintenables.

### Objectifs d'apprentissage

Cette quête secondaire vous emmène dans un voyage pratique des concepts de base aux modèles prêts pour la production.
Nous allons transformer un simple workflow de lecture de CSV en un pipeline bioinformatique sophistiqué, en le faisant évoluer étape par étape à travers des défis réalistes :

- **Comprendre les frontières :** Distinguer les opérations de dataflow et le scripting, et comprendre comment ils fonctionnent ensemble
- **Manipulation de données :** Extraire, transformer et sous-ensembler des maps et des collections en utilisant des opérateurs puissants
- **Traitement de chaînes :** Analyser des schémas de nommage de fichiers complexes avec des motifs regex et maîtriser l'interpolation de variables
- **Fonctions réutilisables :** Extraire la logique complexe dans des fonctions nommées pour des workflows plus propres et plus maintenables
- **Logique dynamique :** Construire des processus qui s'adaptent à différents types d'entrées et utiliser des closures pour l'allocation dynamique de ressources
- **Routage conditionnel :** Router intelligemment les échantillons à travers différents processus en fonction de leurs caractéristiques de métadonnées
- **Opérations sûres :** Gérer les données manquantes avec élégance grâce aux opérateurs null-safe et valider les entrées avec des messages d'erreur clairs
- **Gestionnaires basés sur la configuration :** Utiliser les gestionnaires d'événements de workflow pour la journalisation, les notifications et la gestion du cycle de vie

### Prérequis

Avant d'entreprendre cette quête secondaire, vous devriez :

- Avoir terminé le tutoriel [Hello Nextflow](../hello_nextflow/README.md) ou un cours équivalent pour débutant·es.
- Être à l'aise avec l'utilisation des concepts et mécanismes de base de Nextflow (processus, canaux, opérateurs, travail avec des fichiers, métadonnées)
- Avoir une familiarité de base avec les constructions de programmation courantes (variables, maps, listes)

Ce tutoriel expliquera les concepts de programmation au fur et à mesure que nous les rencontrerons, vous n'avez donc pas besoin d'une expérience de programmation approfondie.
Nous commencerons par les concepts fondamentaux et construirons jusqu'aux modèles avancés.

---

## 0. Premiers pas

#### Ouvrir le codespace de formation

Si vous ne l'avez pas encore fait, assurez-vous d'ouvrir l'environnement de formation comme décrit dans la [Configuration de l'environnement](../envsetup/index.md).

[![Open in GitHub Codespaces](https://github.com/codespaces/badge.svg)](https://codespaces.new/nextflow-io/training?quickstart=1&ref=master)

#### Se déplacer dans le répertoire du projet

Déplaçons-nous dans le répertoire où se trouvent les fichiers de ce tutoriel.

```bash
cd side-quests/essential_scripting_patterns
```

#### Examiner le matériel

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

Nous utiliserons cet ensemble de données réaliste pour explorer des techniques de programmation pratiques que vous rencontrerez dans de vrais workflows bioinformatiques.

<!-- TODO: Can we make this more domain-agnostic? -->

<!-- TODO: add an assignment statement? #### Review the assignment -->

#### Liste de vérification de préparation

Pensez-vous être prêt·e à vous lancer ?

- [ ] Je comprends l'objectif de ce cours et ses prérequis
- [ ] Mon codespace est opérationnel
- [ ] J'ai défini mon répertoire de travail de manière appropriée
<!-- - [ ] I understand the assignment -->

Si vous pouvez cocher toutes les cases, vous êtes prêt·e à commencer.

---

## 1. Dataflow vs Scripting : Comprendre les Frontières

### 1.1. Identifier ce qui est quoi

Lors de l'écriture de workflows Nextflow, il est important de distinguer le **dataflow** (comment les données se déplacent à travers les canaux et les processus) et le **scripting** (le code qui manipule les données et prend des décisions). Construisons un workflow démontrant comment ils fonctionnent ensemble.

#### 1.1.1. Workflow Nextflow de base

Commençons par un simple workflow qui lit simplement le fichier CSV (nous l'avons déjà fait pour vous dans `main.nf`) :

```groovy title="main.nf" linenums="1"
workflow {
    ch_samples = channel.fromPath("./data/samples.csv")
        .splitCsv(header: true)
        .view()
}
```

Le bloc `workflow` définit la structure de notre pipeline, tandis que `channel.fromPath()` crée un canal à partir d'un chemin de fichier. L'opérateur `.splitCsv()` traite le fichier CSV et convertit chaque ligne en une structure de données map.

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

#### 1.1.2. Ajout de l'opérateur Map

Maintenant, nous allons ajouter du scripting pour transformer les données, en utilisant l'opérateur `.map()` que vous connaissez probablement déjà. Cet opérateur prend une 'closure' où nous pouvons écrire du code pour transformer chaque élément.

!!! note "Note"

    Une **closure** est un bloc de code qui peut être passé et exécuté plus tard. Pensez-y comme une fonction que vous définissez en ligne. Les closures sont écrites avec des accolades `{ }` et peuvent prendre des paramètres. Elles sont fondamentales pour le fonctionnement des opérateurs Nextflow et si vous écrivez du Nextflow depuis un moment, vous les utilisez peut-être déjà sans vous en rendre compte !

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

La closure `{ row -> return row }` prend un paramètre `row` (pourrait être n'importe quel nom : `item`, `sample`, etc.).

Lorsque l'opérateur `.map()` traite chaque élément du canal, il passe cet élément à votre closure. Ici, `row` contient une ligne CSV à la fois.

Appliquez cette modification et exécutez le workflow :

```bash
nextflow run main.nf
```

Vous verrez la même sortie qu'auparavant, car nous retournons simplement l'entrée inchangée. Cela confirme que l'opérateur map fonctionne correctement. Maintenant, commençons à transformer les données.

#### 1.1.3. Création d'une structure de données Map

Maintenant, nous allons écrire une logique de **scripting** à l'intérieur de notre closure pour transformer chaque ligne de données. C'est là que nous traitons des éléments de données individuels plutôt que d'orchestrer le flux de données.

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

Appliquez cette modification et exécutez le workflow :

```bash
nextflow run main.nf
```

??? success "Sortie de la commande"

    ```console
    [id:sample_001, organism:human, tissue:liver, depth:30000000, quality:38.5]
    [id:sample_002, organism:mouse, tissue:brain, depth:25000000, quality:35.2]
    [id:sample_003, organism:human, tissue:kidney, depth:45000000, quality:42.1]
    ```

#### 1.1.4. Ajout de logique conditionnelle

Maintenant, ajoutons plus de scripting - cette fois en utilisant un opérateur ternaire pour prendre des décisions basées sur les valeurs des données.

Effectuez la modification suivante :

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

L'opérateur ternaire est un raccourci pour une instruction if/else qui suit le modèle `condition ? valeur_si_vrai : valeur_si_faux`. Cette ligne signifie : « Si la qualité est supérieure à 40, utiliser 'high', sinon utiliser 'normal' ». Son cousin, l'**opérateur Elvis** (`?:`), fournit des valeurs par défaut lorsque quelque chose est null ou vide - nous explorerons ce modèle plus tard dans ce tutoriel.

L'opérateur d'addition de map `+` crée une **nouvelle map** plutôt que de modifier celle existante. Cette ligne crée une nouvelle map qui contient toutes les paires clé-valeur de `sample_meta` plus la nouvelle clé `priority`.

!!! Note "Note"

    Ne modifiez jamais les maps passées dans les closures - créez-en toujours de nouvelles en utilisant `+` (par exemple). Dans Nextflow, les mêmes données circulent souvent simultanément à travers plusieurs opérations. Modifier une map sur place peut causer des effets secondaires imprévisibles lorsque d'autres opérations référencent ce même objet. Créer de nouvelles maps garantit que chaque opération a sa propre copie propre.

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

Nous avons réussi à ajouter une logique conditionnelle pour enrichir nos métadonnées avec un niveau de priorité basé sur les scores de qualité.

#### 1.1.5. Sous-ensemble de Maps avec `.subMap()`

Alors que l'opérateur `+` ajoute des clés à une map, parfois vous devez faire l'inverse - extraire uniquement des clés spécifiques. La méthode `.subMap()` est parfaite pour cela.

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
                println "Champs ID uniquement : ${id_only}"

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

    Champs ID uniquement : [id:sample_001, organism:human, tissue:liver]
    Champs ID uniquement : [id:sample_002, organism:mouse, tissue:brain]
    Champs ID uniquement : [id:sample_003, organism:human, tissue:kidney]
    [id:sample_001, organism:human, tissue:liver, depth:30000000, quality:38.5, priority:normal]
    [id:sample_002, organism:mouse, tissue:brain, depth:25000000, quality:35.2, priority:normal]
    [id:sample_003, organism:human, tissue:kidney, depth:45000000, quality:42.1, priority:high]
    ```

Cela montre à la fois les métadonnées complètes affichées par l'opération `view()` et le sous-ensemble extrait que nous avons imprimé avec `println`.

La méthode `.subMap()` prend une liste de clés et retourne une nouvelle map contenant uniquement ces clés. Si une clé n'existe pas dans la map d'origine, elle n'est simplement pas incluse dans le résultat.

Ceci est particulièrement utile lorsque vous devez créer différentes versions de métadonnées pour différents processus - certains peuvent avoir besoin de métadonnées complètes tandis que d'autres n'ont besoin que de champs d'identification minimaux.

Maintenant, supprimez ces instructions println pour restaurer votre workflow à son état précédent, car nous n'en avons pas besoin pour la suite.

!!! tip "Résumé des opérations sur les Maps"

    - **Ajouter des clés** : `map1 + [new_key: value]` - Crée une nouvelle map avec des clés supplémentaires
    - **Extraire des clés** : `map1.subMap(['key1', 'key2'])` - Crée une nouvelle map avec uniquement les clés spécifiées
    - **Les deux opérations créent de nouvelles maps** - Les maps d'origine restent inchangées

#### 1.1.6. Combiner des Maps et retourner des résultats

Jusqu'à présent, nous n'avons retourné que ce que la communauté Nextflow appelle la 'meta map', et nous avons ignoré les fichiers auxquels ces métadonnées se rapportent. Mais si vous écrivez des workflows Nextflow, vous voulez probablement faire quelque chose avec ces fichiers.

Produisons une structure de canal comprenant un tuple de 2 éléments : la map de métadonnées enrichie et le chemin de fichier correspondant. C'est un modèle courant dans Nextflow pour passer des données aux processus.

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

Appliquez cette modification et exécutez le workflow :

```bash
nextflow run main.nf
```

??? success "Sortie de la commande"

    ```console
    [[id:sample_001, organism:human, tissue:liver, depth:30000000, quality:38.5, priority:normal], /workspaces/training/side-quests/essential_scripting_patterns/data/sequences/SAMPLE_001_S1_L001_R1_001.fastq]
    [[id:sample_002, organism:mouse, tissue:brain, depth:25000000, quality:35.2, priority:normal], /workspaces/training/side-quests/essential_scripting_patterns/data/sequences/SAMPLE_002_S2_L001_R1_001.fastq]
    [[id:sample_003, organism:human, tissue:kidney, depth:45000000, quality:42.1, priority:high], /workspaces/training/side-quests/essential_scripting_patterns/data/sequences/SAMPLE_003_S3_L001_R1_001.fastq]
    ```

Cette structure de tuple `[meta, file]` est un modèle courant dans Nextflow pour passer à la fois des métadonnées et des fichiers associés aux processus.

!!! note "Note"

    **Maps et métadonnées** : Les maps sont fondamentales pour travailler avec les métadonnées dans Nextflow. Pour une explication plus détaillée du travail avec les maps de métadonnées, consultez la quête secondaire [Travailler avec les métadonnées](./metadata.md).

Notre workflow démontre le modèle de base : les **opérations de dataflow** (`workflow`, `channel.fromPath()`, `.splitCsv()`, `.map()`, `.view()`) orchestrent comment les données se déplacent à travers le pipeline, tandis que le **scripting** (maps `[key: value]`, méthodes de chaînes, conversions de types, opérateurs ternaires) à l'intérieur de la closure `.map()` gère la transformation des éléments de données individuels.

### 1.2. Comprendre les différents types : Channel vs List

Jusqu'à présent, tout va bien, nous pouvons distinguer les opérations de dataflow et le scripting. Mais qu'en est-il lorsque le même nom de méthode existe dans les deux contextes ?

Un exemple parfait est la méthode `collect`, qui existe à la fois pour les types de canaux et les types List dans la bibliothèque standard Nextflow. La méthode `collect()` sur une List transforme chaque élément, tandis que l'opérateur `collect()` sur un canal rassemble toutes les émissions de canal en un canal à élément unique.

Démontrons cela avec quelques exemples de données, en commençant par nous rafraîchir la mémoire sur ce que fait l'opérateur `collect()` de canal. Consultez `collect.nf` :

```groovy title="collect.nf" linenums="1"
def sample_ids = ['sample_001', 'sample_002', 'sample_003']

// channel.collect() - regroupe plusieurs émissions de canal en une seule
ch_input = channel.fromList(sample_ids)
ch_input.view { sample -> "Élément de canal individuel : ${sample}" }
ch_collected = ch_input.collect()
ch_collected.view { list -> "Résultat de channel.collect() : ${list} (${list.size()} éléments regroupés en 1)" }
```

Étapes :

- Définir une List d'IDs d'échantillons
- Créer un canal avec `fromList()` qui émet chaque ID d'échantillon séparément
- Imprimer chaque élément avec `view()` au fur et à mesure qu'il circule
- Rassembler tous les éléments dans une seule liste avec l'opérateur `collect()` du canal
- Imprimer le résultat collecté (élément unique contenant tous les IDs d'échantillons) avec un second `view()`

Nous avons changé la structure du canal, mais nous n'avons pas changé les données elles-mêmes.

Exécutez le workflow pour confirmer cela :

```bash
nextflow run collect.nf
```

??? success "Sortie de la commande"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `collect.nf` [loving_mendel] DSL2 - revision: e8d054a46e

    Élément de canal individuel : sample_001
    Élément de canal individuel : sample_002
    Élément de canal individuel : sample_003
    Résultat de channel.collect() : [sample_001, sample_002, sample_003] (3 éléments regroupés en 1)
    ```

`view()` retourne une sortie pour chaque émission de canal, nous savons donc que cette sortie unique contient les 3 éléments d'origine regroupés en une seule liste.

Maintenant, voyons la méthode `collect` sur une List en action. Modifiez `collect.nf` pour appliquer la méthode `collect` de la List à la liste d'origine des IDs d'échantillons :

=== "Après"

    ```groovy title="main.nf" linenums="1" hl_lines="9-13"
    def sample_ids = ['sample_001', 'sample_002', 'sample_003']

    // channel.collect() - regroupe plusieurs émissions de canal en une seule
    ch_input = channel.fromList(sample_ids)
    ch_input.view { sample -> "Élément de canal individuel : ${sample}" }
    ch_collected = ch_input.collect()
    ch_collected.view { list -> "Résultat de channel.collect() : ${list} (${list.size()} éléments regroupés en 1)" }

    // List.collect() - transforme chaque élément, préserve la structure
    def formatted_ids = sample_ids.collect { id ->
        id.toUpperCase().replace('SAMPLE_', 'SPECIMEN_')
    }
    println "Résultat de List.collect() : ${formatted_ids} (${sample_ids.size()} éléments transformés en ${formatted_ids.size()})"
    ```

=== "Avant"

    ```groovy title="main.nf" linenums="1"
    def sample_ids = ['sample_001', 'sample_002', 'sample_003']

    // channel.collect() - regroupe plusieurs émissions de canal en une seule
    ch_input = channel.fromList(sample_ids)
    ch_input.view { sample -> "Élément de canal individuel : ${sample}" }
    ch_collected = ch_input.collect()
    ch_collected.view { list -> "Résultat de channel.collect() : ${list} (${list.size()} éléments regroupés en 1)" }
    ```

Dans ce nouveau snippet, nous :

- Définissons une nouvelle variable `formatted_ids` qui utilise la méthode `collect` de la List pour transformer chaque ID d'échantillon dans la liste d'origine
- Imprimons le résultat en utilisant `println`

Exécutez le workflow modifié :

```bash
nextflow run collect.nf
```

??? success "Sortie de la commande"

    ```console hl_lines="5"
    N E X T F L O W   ~  version 25.10.2

    Launching `collect.nf` [cheeky_stonebraker] DSL2 - revision: 2d5039fb47

    Résultat de List.collect() : [SPECIMEN_001, SPECIMEN_002, SPECIMEN_003] (3 éléments transformés en 3)
    Élément de canal individuel : sample_001
    Élément de canal individuel : sample_002
    Élément de canal individuel : sample_003
    Résultat de channel.collect() : [sample_001, sample_002, sample_003] (3 éléments regroupés en 1)
    ```

Cette fois, nous n'avons PAS changé la structure des données, nous avons toujours 3 éléments dans la liste, mais nous AVONS transformé chaque élément en utilisant la méthode `collect` de la List pour produire une nouvelle liste avec des valeurs modifiées. C'est similaire à l'utilisation de l'opérateur `map` sur un canal, mais il opère sur une structure de données List plutôt que sur un canal.

`collect` est un cas extrême que nous utilisons ici pour faire valoir un point. La leçon clé est que lorsque vous écrivez des workflows, distinguez toujours les **structures de données** (Lists, Maps, etc.) et les **canaux** (constructions de dataflow). Les opérations peuvent partager des noms mais se comporter complètement différemment selon le type sur lequel elles sont appelées.

### 1.3. L'opérateur Spread (`*.`) - Raccourci pour l'extraction de propriétés

Lié à la méthode `collect` de la List est l'opérateur spread (`*.`), qui fournit un moyen concis d'extraire des propriétés de collections. C'est essentiellement du sucre syntaxique pour un modèle `collect` courant.

Ajoutons une démonstration à notre fichier `collect.nf` :

=== "Après"

    ```groovy title="collect.nf" linenums="1" hl_lines="15-18"
    def sample_ids = ['sample_001', 'sample_002', 'sample_003']

    // channel.collect() - regroupe plusieurs émissions de canal en une seule
    ch_input = channel.fromList(sample_ids)
    ch_input.view { sample -> "Élément de canal individuel : ${sample}" }
    ch_collected = ch_input.collect()
    ch_collected.view { list -> "Résultat de channel.collect() : ${list} (${list.size()} éléments regroupés en 1)" }

    // List.collect() - transforme chaque élément, préserve la structure
    def formatted_ids = sample_ids.collect { id ->
        id.toUpperCase().replace('SAMPLE_', 'SPECIMEN_')
    }
    println "Résultat de List.collect() : ${formatted_ids} (${sample_ids.size()} éléments transformés en ${formatted_ids.size()})"

    // Opérateur spread - accès concis aux propriétés
    def sample_data = [[id: 's1', quality: 38.5], [id: 's2', quality: 42.1], [id: 's3', quality: 35.2]]
    def all_ids = sample_data*.id
    println "Résultat de l'opérateur spread : ${all_ids}"
    ```

=== "Avant"

    ```groovy title="collect.nf" linenums="1"
    def sample_ids = ['sample_001', 'sample_002', 'sample_003']

    // channel.collect() - regroupe plusieurs émissions de canal en une seule
    ch_input = channel.fromList(sample_ids)
    ch_input.view { sample -> "Élément de canal individuel : ${sample}" }
    ch_collected = ch_input.collect()
    ch_collected.view { list -> "Résultat de channel.collect() : ${list} (${list.size()} éléments regroupés en 1)" }

    // List.collect() - transforme chaque élément, préserve la structure
    def formatted_ids = sample_ids.collect { id ->
        id.toUpperCase().replace('SAMPLE_', 'SPECIMEN_')
    }
    println "Résultat de List.collect() : ${formatted_ids} (${sample_ids.size()} éléments transformés en ${formatted_ids.size()})"
    ```

Exécutez le workflow mis à jour :

```bash title="Tester l'opérateur spread"
nextflow run collect.nf
```

??? success "Sortie de la commande"

    ```console hl_lines="6"
    N E X T F L O W   ~  version 25.10.2

    Launching `collect.nf` [cranky_galileo] DSL2 - revision: 5f3c8b2a91

    Résultat de List.collect() : [SPECIMEN_001, SPECIMEN_002, SPECIMEN_003] (3 éléments transformés en 3)
    Résultat de l'opérateur spread : [s1, s2, s3]
    Élément de canal individuel : sample_001
    Élément de canal individuel : sample_002
    Élément de canal individuel : sample_003
    Résultat de channel.collect() : [sample_001, sample_002, sample_003] (3 éléments regroupés en 1)
    ```

L'opérateur spread `*.` est un raccourci pour un modèle collect courant :

```groovy
// Ces deux expressions sont équivalentes :
def ids = samples*.id
def ids = samples.collect { it.id }

// Fonctionne aussi avec les appels de méthode :
def names = files*.getName()
def names = files.collect { it.getName() }
```

L'opérateur spread est particulièrement utile lorsque vous devez extraire une seule propriété d'une liste d'objets - il est plus lisible que d'écrire la closure `collect` complète.

!!! tip "Quand utiliser Spread vs Collect"

    - **Utilisez spread (`*.`)** pour un accès simple aux propriétés : `samples*.id`, `files*.name`
    - **Utilisez collect** pour les transformations ou la logique complexe : `samples.collect { it.id.toUpperCase() }`, `samples.collect { [it.id, it.quality > 40] }`

### À retenir

Dans cette section, vous avez appris :

- **Dataflow vs scripting** : Les opérateurs de canal orchestrent comment les données circulent à travers votre pipeline, tandis que le scripting transforme les éléments de données individuels
- **Comprendre les types** : La même méthode (comme `collect`) peut se comporter différemment selon le type sur lequel elle est appelée (Channel vs List)
- **Le contexte compte** : Soyez toujours conscient·e de si vous travaillez avec des canaux (dataflow) ou des structures de données (scripting)

Comprendre ces frontières est essentiel pour le débogage, la documentation et l'écriture de workflows maintenables.

Ensuite, nous approfondirons les capacités de traitement de chaînes, qui sont essentielles pour gérer des données du monde réel.

---

## 2. Traitement de chaînes et génération de scripts dynamiques

Maîtriser le traitement de chaînes sépare les workflows fragiles des pipelines robustes. Cette section couvre l'analyse de noms de fichiers complexes, la génération de scripts dynamiques et l'interpolation de variables.

### 2.1. Correspondance de motifs et expressions régulières

Les fichiers bioinformatiques ont souvent des conventions de nommage complexes encodant des métadonnées. Extrayons cela automatiquement en utilisant la correspondance de motifs avec des expressions régulières.

Nous allons retourner à notre workflow `main.nf` et ajouter une logique de correspondance de motifs pour extraire des informations supplémentaires sur les échantillons à partir des noms de fichiers. Les fichiers FASTQ de notre ensemble de données suivent les conventions de nommage de style Illumina avec des noms comme `SAMPLE_001_S1_L001_R1_001.fastq.gz`. Ceux-ci peuvent sembler cryptiques, mais ils encodent en fait des métadonnées utiles comme l'ID d'échantillon, le numéro de lane et la direction de lecture. Nous allons utiliser les capacités regex pour analyser ces noms.

Effectuez la modification suivante dans votre workflow `main.nf` existant :

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

Cela démontre des **concepts clés de traitement de chaînes** :

1. **Littéraux d'expressions régulières** utilisant la syntaxe `~/pattern/` - cela crée un motif regex sans avoir besoin d'échapper les barres obliques inverses
2. **Correspondance de motifs** avec l'opérateur `=~` - cela tente de faire correspondre une chaîne à un motif regex
3. **Objets Matcher** qui capturent des groupes avec `[0][1]`, `[0][2]`, etc. - `[0]` fait référence à la correspondance entière, `[1]`, `[2]`, etc. font référence aux groupes capturés entre parenthèses

Décomposons le motif regex `^(.+)_S(\d+)_L(\d{3})_(R[12])_(\d{3})\.fastq(?:\.gz)?$` :

| Motif               | Correspond à                                 | Capture                            |
| ------------------- | -------------------------------------------- | ---------------------------------- |
| `^(.+)`             | Nom d'échantillon depuis le début            | Groupe 1 : nom d'échantillon       |
| `_S(\d+)`           | Numéro d'échantillon `_S1`, `_S2`, etc.      | Groupe 2 : numéro d'échantillon    |
| `_L(\d{3})`         | Numéro de lane `_L001`                       | Groupe 3 : lane (3 chiffres)       |
| `_(R[12])`          | Direction de lecture `_R1` ou `_R2`          | Groupe 4 : direction de lecture    |
| `_(\d{3})`          | Numéro de chunk `_001`                       | Groupe 5 : chunk (3 chiffres)      |
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

### 2.2. Génération de scripts dynamiques dans les processus

Les blocs script de processus sont essentiellement des chaînes multi-lignes qui sont passées au shell. Vous pouvez utiliser une **logique conditionnelle** (if/else, opérateurs ternaires) pour générer dynamiquement différentes chaînes de script en fonction des caractéristiques d'entrée. Ceci est essentiel pour gérer divers types d'entrées -- comme les lectures single-end vs paired-end -- sans dupliquer les définitions de processus.

Ajoutons un processus à notre workflow qui démontre ce modèle. Ouvrez `modules/fastp.nf` et jetez un œil :

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

Le processus prend des fichiers FASTQ en entrée et exécute l'outil `fastp` pour couper les adaptateurs et filtrer les lectures de faible qualité. Malheureusement, la personne qui a écrit ce processus n'a pas prévu les lectures single-end que nous avons dans notre ensemble de données d'exemple. Ajoutons-le à notre workflow et voyons ce qui se passe :

Tout d'abord, incluez le module à la toute première ligne de votre workflow `main.nf` :

```groovy title="main.nf" linenums="1"
include { FASTP } from './modules/fastp.nf'
```

Ensuite, modifiez le bloc `workflow` pour connecter le canal `ch_samples` au processus `FASTP` :

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

Vous pouvez voir que le processus essaie d'exécuter `fastp` avec une valeur `null` pour le deuxième fichier d'entrée, ce qui le fait échouer. C'est parce que notre ensemble de données contient des lectures single-end, mais le processus est codé en dur pour attendre des lectures paired-end (deux fichiers d'entrée à la fois).

Corrigez cela en ajoutant une logique conditionnelle au bloc `script:` du processus `FASTP`. Une instruction if/else vérifie le nombre de fichiers de lecture et ajuste la commande en conséquence.

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

Maintenant, le workflow peut gérer à la fois les lectures single-end et paired-end avec élégance. La logique conditionnelle vérifie le nombre de fichiers d'entrée et construit la commande appropriée pour `fastp`. Voyons si cela fonctionne :

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

Une autre utilisation courante de la logique de script dynamique peut être vue dans [le module Nextflow for Science Genomics](../../nf4science/genomics/02_joint_calling). Dans ce module, le processus GATK appelé peut prendre plusieurs fichiers d'entrée, mais chacun doit être préfixé par `-V` pour former une ligne de commande correcte. Le processus utilise le scripting pour transformer une collection de fichiers d'entrée (`all_gvcfs`) en arguments de commande corrects :

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

Ces modèles d'utilisation du scripting dans les blocs script de processus sont extrêmement puissants et peuvent être appliqués dans de nombreux scénarios - de la gestion de types d'entrée variables à la construction d'arguments de ligne de commande complexes à partir de collections de fichiers, rendant vos processus vraiment adaptables aux exigences diverses des données du monde réel.

### 2.3. Interpolation de variables : Variables Nextflow et Shell

Les scripts de processus mélangent des variables Nextflow, des variables shell et des substitutions de commandes, chacune avec une syntaxe d'interpolation différente. Utiliser la mauvaise syntaxe provoque des erreurs. Explorons cela avec un processus qui crée un rapport de traitement.

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
    echo "Traitement de ${reads}" > ${meta.id}_report.txt
    echo "Échantillon : ${meta.id}" >> ${meta.id}_report.txt
    """
}
```

Ce processus écrit un rapport simple avec l'ID d'échantillon et le nom de fichier. Maintenant, exécutons-le pour voir ce qui se passe lorsque nous devons mélanger différents types de variables.

Incluez le processus dans votre `main.nf` et ajoutez-le au workflow :

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

Mais que se passe-t-il si nous voulons ajouter des informations sur quand et où le traitement a eu lieu ? Modifions le processus pour utiliser des variables **shell** et un peu de substitution de commande pour inclure l'utilisateur·trice actuel·le, le nom d'hôte et la date dans le rapport :

=== "Après"

    ```groovy title="modules/generate_report.nf" linenums="10" hl_lines="5-7"
        script:
        """
        echo "Traitement de ${reads}" > ${meta.id}_report.txt
        echo "Échantillon : ${meta.id}" >> ${meta.id}_report.txt
        echo "Traité par : ${USER}" >> ${meta.id}_report.txt
        echo "Nom d'hôte : $(hostname)" >> ${meta.id}_report.txt
        echo "Date : $(date)" >> ${meta.id}_report.txt
        """
    ```

=== "Avant"

    ```groovy title="modules/generate_report.nf" linenums="10"
        script:
        """
        echo "Traitement de ${reads}" > ${meta.id}_report.txt
        echo "Échantillon : ${meta.id}" >> ${meta.id}_report.txt
        """
    ```

Si vous exécutez cela, vous remarquerez une erreur - Nextflow essaie d'interpréter `${USER}` comme une variable Nextflow qui n'existe pas.

??? failure "Sortie de la commande"

    ```console
    Error modules/generate_report.nf:15:27: `USER` is not defined
    │  15 |     echo "Traité par : ${USER}" >> ${meta.id}_report.txt
    ╰     |                           ^^^^

    ERROR ~ Script compilation failed
    ```

Nous devons l'échapper pour que Bash puisse le gérer à la place.

Corrigez cela en échappant les variables shell et les substitutions de commandes avec une barre oblique inverse (`\`) :

=== "Après"

    ```groovy title="modules/generate_report.nf" linenums="10" hl_lines="5-7"
        script:
        """
        echo "Traitement de ${reads}" > ${meta.id}_report.txt
        echo "Échantillon : ${meta.id}" >> ${meta.id}_report.txt
        echo "Traité par : \${USER}" >> ${meta.id}_report.txt
        echo "Nom d'hôte : \$(hostname)" >> ${meta.id}_report.txt
        echo "Date : \$(date)" >> ${meta.id}_report.txt
        """
    ```

=== "Avant"

    ```groovy title="modules/generate_report.nf" linenums="10"
        script:
        """
        echo "Traitement de ${reads}" > ${meta.id}_report.txt
        echo "Échantillon : ${meta.id}" >> ${meta.id}_report.txt
        echo "Traité par : ${USER}" >> ${meta.id}_report.txt
        echo "Nom d'hôte : $(hostname)" >> ${meta.id}_report.txt
        echo "Date : $(date)" >> ${meta.id}_report.txt
        """
    ```

Maintenant ça fonctionne ! La barre oblique inverse (`\`) dit à Nextflow « n'interprète pas ceci, passe-le à Bash ».

### À retenir

Dans cette section, vous avez appris des techniques de **traitement de chaînes** :

- **Expressions régulières pour l'analyse de fichiers** : Utilisation de l'opérateur `=~` et des motifs regex (`~/pattern/`) pour extraire des métadonnées de conventions de nommage de fichiers complexes
- **Génération de scripts dynamiques** : Utilisation de logique conditionnelle (if/else, opérateurs ternaires) pour générer différentes chaînes de script en fonction des caractéristiques d'entrée
- **Interpolation de variables** : Comprendre quand Nextflow interprète les chaînes vs quand le shell le fait
  - `${var}` - Variables Nextflow (interpolées par Nextflow au moment de la compilation du workflow)
  - `\${var}` - Variables d'environnement shell (échappées, passées à bash au moment de l'exécution)
  - `\$(cmd)` - Substitution de commande shell (échappée, exécutée par bash au moment de l'exécution)

Ces modèles de traitement et de génération de chaînes sont essentiels pour gérer les divers formats de fichiers et conventions de nommage que vous rencontrerez dans les workflows bioinformatiques du monde réel.

---

## 3. Création de fonctions réutilisables

La logique de workflow complexe en ligne dans les opérateurs de canal ou les définitions de processus réduit la lisibilité et la maintenabilité. Les **fonctions** vous permettent d'extraire cette logique dans des composants nommés et réutilisables.

Notre opération map est devenue longue et complexe. Extrayons-la dans une fonction réutilisable en utilisant le mot-clé `def`.

Pour illustrer à quoi cela ressemble avec notre workflow existant, effectuez la modification ci-dessous, en utilisant `def` pour définir une fonction réutilisable appelée `separateMetadata` :

=== "Après"

    ```groovy title="main.nf" linenums="1" hl_lines="4-24 29"
    include { FASTP } from './modules/fastp.nf'
    include { GENERATE_REPORT } from './modules/generate_report.nf'

    def separateMetadata(row) {
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

    workflow {
        ch_samples = channel.fromPath("./data/samples.csv")
            .splitCsv(header: true)
            .map{ row -> separateMetadata(row) }

        ch_fastp = FASTP(ch_samples)
        GENERATE_REPORT(ch_samples)
    }
    ```

=== "Avant"

    ```groovy title="main.nf" linenums="1" hl_lines="7-27"
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

En extrayant cette logique dans une fonction, nous avons réduit la logique réelle du workflow à quelque chose de beaucoup plus propre :

```groovy title="workflow minimal"
    ch_samples = channel.fromPath("./data/samples.csv")
        .splitCsv(header: true)
        .map{ row -> separateMetadata(row) }

    ch_fastp = FASTP(ch_samples)
    GENERATE_REPORT(ch_samples)
```

Cela rend la logique du workflow beaucoup plus facile à lire et à comprendre d'un coup d'œil. La fonction `separateMetadata` encapsule toute la logique complexe d'analyse et d'enrichissement des métadonnées, la rendant réutilisable et testable.

Exécutez le workflow pour vous assurer qu'il fonctionne toujours :

```bash
nextflow run main.nf
```

??? success "Sortie de la commande"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `main.nf` [admiring_panini] DSL2 - revision: 8cc832e32f

    executor >  local (6)
    [8c/2e3f91] process > FASTP (3)           [100%] 3 of 3 ✔
    [7a/1b4c92] process > GENERATE_REPORT (3) [100%] 3 of 3 ✔
    ```

La sortie devrait montrer les deux processus se terminant avec succès. Le workflow est maintenant beaucoup plus propre et plus facile à maintenir, avec toute la logique complexe de traitement des métadonnées encapsulée dans la fonction `separateMetadata`.

### À retenir

Dans cette section, vous avez appris la **création de fonctions** :

- **Définir des fonctions avec `def`** : Le mot-clé pour créer des fonctions nommées (comme `def` en Python ou `function` en JavaScript)
- **Portée des fonctions** : Les fonctions définies au niveau du script sont accessibles dans tout votre workflow Nextflow
- **Valeurs de retour** : Les fonctions retournent automatiquement la dernière expression, ou utilisent un `return` explicite
- **Code plus propre** : Extraire la logique complexe dans des fonctions est une pratique fondamentale d'ingénierie logicielle dans n'importe quel langage

Ensuite, nous explorerons comment utiliser les closures dans les directives de processus pour l'allocation dynamique de ressources.

---

## 4. Directives de ressources dynamiques avec des closures

Jusqu'à présent, nous avons utilisé le scripting dans le bloc `script` des processus. Mais les **closures** (introduites dans la Section 1.1) sont également incroyablement utiles dans les directives de processus, en particulier pour l'allocation dynamique de ressources. Ajoutons des directives de ressources à notre processus FASTP qui s'adaptent en fonction des caractéristiques de l'échantillon.

### 4.1. Allocation de ressources spécifique aux échantillons

Actuellement, notre processus FASTP utilise des ressources par défaut. Rendons-le plus intelligent en allouant plus de CPUs pour les échantillons à haute profondeur. Modifiez `modules/fastp.nf` pour inclure une directive `cpus` dynamique et une directive `memory` statique :

=== "Après"

    ```groovy title="modules/fastp.nf" linenums="1" hl_lines="4-5"
    process FASTP {
        container 'community.wave.seqera.io/library/fastp:0.24.0--62c97b06e8447690'

        cpus { meta.depth > 40000000 ? 2 : 1 }
        memory 2.GB

        input:
        tuple val(meta), path(reads)
    ```

=== "Avant"

    ```groovy title="modules/fastp.nf" linenums="1"
    process FASTP {
        container 'community.wave.seqera.io/library/fastp:0.24.0--62c97b06e8447690'

        input:
        tuple val(meta), path(reads)
    ```

La closure `{ meta.depth > 40000000 ? 2 : 1 }` utilise l'**opérateur ternaire** (couvert dans la Section 1.1) et est évaluée pour chaque tâche, permettant une allocation de ressources par échantillon. Les échantillons à haute profondeur (>40M lectures) obtiennent 2 CPUs, tandis que les autres obtiennent 1 CPU.

!!! note "Accès aux variables d'entrée dans les directives"

    La closure peut accéder à toutes les variables d'entrée (comme `meta` ici) car Nextflow évalue ces closures dans le contexte de chaque exécution de tâche.

Exécutez à nouveau le workflow avec l'option `-ansi-log false` pour faciliter la visualisation des hashs de tâches.

```bash
nextflow run main.nf -ansi-log false
```

??? success "Sortie de la commande"

    ```console
    N E X T F L O W  ~  version 25.10.2
    Launching `main.nf` [fervent_albattani] DSL2 - revision: fa8f249759
    [bd/ff3d41] Submitted process > FASTP (2)
    [a4/a3aab2] Submitted process > FASTP (1)
    [48/6db0c9] Submitted process > FASTP (3)
    [ec/83439d] Submitted process > GENERATE_REPORT (3)
    [bd/15d7cc] Submitted process > GENERATE_REPORT (2)
    [42/699357] Submitted process > GENERATE_REPORT (1)
    ```

Vous pouvez vérifier la commande `docker` exacte qui a été exécutée pour voir l'allocation de CPU pour une tâche donnée :

```console title="Vérifier la commande docker"
cat work/48/6db0c9e9d8aa65e4bb4936cd3bd59e/.command.run | grep "docker run"
```

Vous devriez voir quelque chose comme :

```bash title="commande docker"
    docker run -i --cpu-shares 4096 --memory 2048m -e "NXF_TASK_WORKDIR" -v /workspaces/training/side-quests/essential_scripting_patterns:/workspaces/training/side-quests/essential_scripting_patterns -w "$NXF_TASK_WORKDIR" --name $NXF_BOXID community.wave.seqera.io/library/fastp:0.24.0--62c97b06e8447690 /bin/bash -ue /workspaces/training/side-quests/essential_scripting_patterns/work/48/6db0c9e9d8aa65e4bb4936cd3bd59e/.command.sh
```

Dans cet exemple, nous avons choisi un exemple qui a demandé 2 CPUs (`--cpu-shares 2048`), car c'était un échantillon à haute profondeur, mais vous devriez voir différentes allocations de CPU selon la profondeur de l'échantillon. Essayez cela pour les autres tâches également.

### 4.2. Stratégies de nouvelle tentative

Un autre modèle puissant consiste à utiliser `task.attempt` pour les stratégies de nouvelle tentative. Pour montrer pourquoi c'est utile, nous allons commencer par réduire l'allocation de mémoire à FASTP à moins que ce dont il a besoin. Changez la directive `memory` dans `modules/fastp.nf` à `1.GB` :

=== "Après"

    ```groovy title="modules/fastp.nf" linenums="1" hl_lines="5"
    process FASTP {
        container 'community.wave.seqera.io/library/fastp:0.24.0--62c97b06e8447690'

        cpus { meta.depth > 40000000 ? 4 : 2 }
        memory 1.GB

        input:
        tuple val(meta), path(reads)
    ```

=== "Avant"

    ```groovy title="modules/fastp.nf" linenums="1" hl_lines="5"
    process FASTP {
        container 'community.wave.seqera.io/library/fastp:0.24.0--62c97b06e8447690'

        cpus { meta.depth > 40000000 ? 4 : 2 }
        memory 2.GB

        input:
        tuple val(meta), path(reads)
    ```

... et exécutez à nouveau le workflow :

```bash
nextflow run main.nf
```

??? failure "Sortie de la commande"

    ```console hl_lines="2 11"
    Command exit status:
      137

    Command output:
      (empty)

    Command error:
      Detecting adapter sequence for read1...
      No adapter detected for read1

      .command.sh: line 7:   101 Killed                  fastp --in1 SAMPLE_002_S2_L001_R1_001.fastq --out1 sample_002_trimmed.fastq.gz --json sample_002.fastp.json --html sample_002.fastp.html --thread 2
    ```

Cela indique que le processus a été tué pour avoir dépassé les limites de mémoire.

C'est un scénario très courant dans les workflows du monde réel - parfois vous ne savez tout simplement pas combien de mémoire une tâche aura besoin jusqu'à ce que vous l'exécutiez.

Pour rendre notre workflow plus robuste, nous pouvons implémenter une stratégie de nouvelle tentative qui augmente l'allocation de mémoire à chaque tentative, encore une fois en utilisant une closure Groovy. Modifiez la directive `memory` pour multiplier la mémoire de base par `task.attempt`, et ajoutez les directives `errorStrategy 'retry'` et `maxRetries 2` :

=== "Après"

    ```groovy title="modules/fastp.nf" linenums="1" hl_lines="5-7"
    process FASTP {
        container 'community.wave.seqera.io/library/fastp:0.24.0--62c97b06e8447690'

        cpus { meta.depth > 40000000 ? 4 : 2 }
        memory { 1.GB * task.attempt }
        errorStrategy 'retry'
        maxRetries 2

        input:
        tuple val(meta), path(reads)
    ```

=== "Avant"

    ```groovy title="modules/fastp.nf" linenums="1" hl_lines="5"
    process FASTP {
        container 'community.wave.seqera.io/library/fastp:0.24.0--62c97b06e8447690'

        cpus { meta.depth > 40000000 ? 4 : 2 }
        memory 2.GB

        input:
        tuple val(meta), path(reads)
    ```

Maintenant, si le processus échoue en raison d'une mémoire insuffisante, Nextflow réessaiera avec plus de mémoire :

- Première tentative : 1 GB (task.attempt = 1)
- Deuxième tentative : 2.GB (task.attempt = 2)

... et ainsi de suite, jusqu'à la limite `maxRetries`.

### À retenir

Les directives dynamiques avec des closures vous permettent de :

- Allouer des ressources en fonction des caractéristiques d'entrée
- Implémenter des stratégies de nouvelle tentative automatiques avec des ressources croissantes
- Combiner plusieurs facteurs (métadonnées, numéro de tentative, priorités)
- Utiliser une logique conditionnelle pour des calculs de ressources complexes

Cela rend vos workflows à la fois plus efficaces (pas de sur-allocation) et plus robustes (nouvelle tentative automatique avec plus de ressources).

---

## 5. Logique conditionnelle et contrôle de processus

Précédemment, nous avons utilisé `.map()` avec du scripting pour transformer les données de canal. Maintenant, nous allons utiliser la logique conditionnelle pour contrôler quels processus s'exécutent en fonction des données -- essentiel pour des workflows flexibles s'adaptant à différents types d'échantillons.

Les [opérateurs de dataflow](https://www.nextflow.io/docs/latest/reference/operator.html) de Nextflow prennent des closures évaluées au moment de l'exécution, permettant à la logique conditionnelle de piloter les décisions de workflow en fonction du contenu du canal.

### 5.1. Routage avec `.branch()`

Par exemple, supposons que nos échantillons de séquençage doivent être coupés avec FASTP uniquement s'il s'agit d'échantillons humains avec une couverture supérieure à un certain seuil. Les échantillons de souris ou les échantillons à faible couverture devraient être exécutés avec Trimgalore à la place (c'est un exemple artificiel, mais il illustre le point).

Nous avons fourni un processus Trimgalore simple dans `modules/trimgalore.nf`, jetez-y un œil si vous le souhaitez, mais les détails ne sont pas importants pour cet exercice. Le point clé est que nous voulons router les échantillons en fonction de leurs métadonnées.

Incluez le nouveau module depuis `modules/trimgalore.nf` :

=== "Après"

    ```groovy title="main.nf" linenums="1" hl_lines="2"
    include { FASTP } from './modules/fastp.nf'
    include { TRIMGALORE } from './modules/trimgalore.nf'
    ```

=== "Avant"

    ```groovy title="main.nf" linenums="1"
    include { FASTP } from './modules/fastp.nf'
    ```

... puis modifiez votre workflow `main.nf` pour brancher les échantillons en fonction de leurs métadonnées et les router à travers le processus de coupe approprié, comme ceci :

=== "Après"

    ```groovy title="main.nf" linenums="28" hl_lines="5-12"
        ch_samples = channel.fromPath("./data/samples.csv")
            .splitCsv(header: true)
            .map { row -> separateMetadata(row) }

        trim_branches = ch_samples
            .branch { meta, reads ->
                fastp: meta.organism == 'human' && meta.depth >= 30000000
                trimgalore: true
            }

        ch_fastp = FASTP(trim_branches.fastp)
        ch_trimgalore = TRIMGALORE(trim_branches.trimgalore)
        GENERATE_REPORT(ch_samples)
    ```

=== "Avant"

    ```groovy title="main.nf" linenums="28" hl_lines="5"
        ch_samples = channel.fromPath("./data/samples.csv")
            .splitCsv(header: true)
            .map { row -> separateMetadata(row) }

        ch_fastp = FASTP(ch_samples)
        GENERATE_REPORT(ch_samples)
    ```

Exécutez ce workflow modifié :

```bash
nextflow run main.nf
```

??? success "Sortie de la commande"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `main.nf` [adoring_galileo] DSL2 - revision: c9e83aaef1

    executor >  local (6)
    [1d/0747ac] process > FASTP (2)           [100%] 2 of 2 ✔
    [cc/c44caf] process > TRIMGALORE (1)      [100%] 1 of 1 ✔
    [34/bd5a9f] process > GENERATE_REPORT (1) [100%] 3 of 3 ✔
    ```

Ici, nous avons utilisé de petites mais puissantes expressions conditionnelles à l'intérieur de l'opérateur `.branch{}` pour router les échantillons en fonction de leurs métadonnées. Les échantillons humains avec une couverture élevée passent par `FASTP`, tandis que tous les autres échantillons passent par `TRIMGALORE`.

### 5.2. Utilisation de `.filter()` avec la véracité

Un autre modèle puissant pour contrôler l'exécution du workflow est l'opérateur `.filter()`, qui utilise une closure pour déterminer quels éléments doivent continuer dans le pipeline. À l'intérieur de la closure de filtre, vous écrirez des **expressions booléennes** qui décident quels éléments passent.

Nextflow (comme de nombreux langages dynamiques) a un concept de **"véracité"** qui détermine quelles valeurs s'évaluent à `true` ou `false` dans des contextes booléens :

- **Vrai (Truthy)** : Valeurs non nulles, chaînes non vides, nombres non nuls, collections non vides
- **Faux (Falsy)** : `null`, chaînes vides `""`, zéro `0`, collections vides `[]` ou `[:]`, `false`

Cela signifie que `meta.id` seul (sans `!= null` explicite) vérifie si l'ID existe et n'est pas vide. Utilisons cela pour filtrer les échantillons qui ne répondent pas à nos exigences de qualité.

Ajoutez ce qui suit avant l'opération de branchement :

=== "Après"

    ```groovy title="main.nf" linenums="28" hl_lines="5-11"
        ch_samples = channel.fromPath("./data/samples.csv")
            .splitCsv(header: true)
            .map { row -> separateMetadata(row) }

        // Filtrer les échantillons invalides ou de faible qualité
        ch_valid_samples = ch_samples
            .filter { meta, reads ->
                meta.id && meta.organism && meta.depth >= 25000000
            }

        trim_branches = ch_valid_samples
            .branch { meta, reads ->
                fastp: meta.organism == 'human' && meta.depth >= 30000000
                trimgalore: true
            }
    ```

=== "Avant"

    ```groovy title="main.nf" linenums="28" hl_lines="5"
        ch_samples = channel.fromPath("./data/samples.csv")
            .splitCsv(header: true)
            .map { row -> separateMetadata(row) }

        trim_branches = ch_samples
            .branch { meta, reads ->
                fastp: meta.organism == 'human' && meta.depth >= 30000000
                trimgalore: true
            }
    ```

Exécutez à nouveau le workflow :

```bash
nextflow run main.nf
```

??? success "Sortie de la commande"

    ```console hl_lines="5"
    N E X T F L O W  ~  version 25.10.2
    Launching `main.nf` [lonely_williams] DSL2 - revision: d0b3f121ec
    [94/b48eac] Submitted process > FASTP (2)
    [2c/d2b28f] Submitted process > GENERATE_REPORT (2)
    [65/2e3be4] Submitted process > GENERATE_REPORT (1)
    [94/b48eac] NOTE: Process `FASTP (2)` terminated with an error exit status (137) -- Execution is retried (1)
    [3e/0d8664] Submitted process > TRIMGALORE (1)
    [6a/9137b0] Submitted process > FASTP (1)
    [6a/9137b0] NOTE: Process `FASTP (1)` terminated with an error exit status (137) -- Execution is retried (1)
    [83/577ac0] Submitted process > GENERATE_REPORT (3)
    [a2/5117de] Re-submitted process > FASTP (1)
    [1f/a1a4ca] Re-submitted process > FASTP (2)
    ```

Parce que nous avons choisi un filtre qui exclut certains échantillons, moins de tâches ont été exécutées.

L'expression de filtre `meta.id && meta.organism && meta.depth >= 25000000` combine la véracité avec des comparaisons explicites :

- `meta.id && meta.organism` vérifie que les deux champs existent et ne sont pas vides (en utilisant la véracité)
- `meta.depth >= 25000000` assure une profondeur de séquençage suffisante avec une comparaison explicite

!!! note "Véracité en pratique"

    L'expression `meta.id && meta.organism` est plus concise que d'écrire :
    ```groovy
    meta.id != null && meta.id != '' && meta.organism != null && meta.organism != ''
    ```

    Cela rend la logique de filtrage beaucoup plus propre et plus facile à lire.

### À retenir

Dans cette section, vous avez appris à utiliser la logique conditionnelle pour contrôler l'exécution du workflow en utilisant les interfaces de closure des opérateurs Nextflow comme `.branch{}` et `.filter{}`, en tirant parti de la véracité pour écrire des expressions conditionnelles concises.

Notre pipeline route maintenant intelligemment les échantillons à travers les processus appropriés, mais les workflows de production doivent gérer les données invalides avec élégance. Rendons notre workflow robuste contre les valeurs manquantes ou nulles.

---

## 6. Navigation sûre et opérateurs Elvis

Notre fonction `separateMetadata` suppose actuellement que tous les champs CSV sont présents et valides. Mais que se passe-t-il avec des données incomplètes ? Découvrons-le.

### 6.1. Le problème : Accéder à des propriétés qui n'existent pas

Disons que nous voulons ajouter la prise en charge d'informations de séquençage optionnelles. Dans certains laboratoires, les échantillons peuvent avoir un champ supplémentaire pour l'ID de séquençage ou le numéro de lot, mais notre CSV actuel n'a pas cette colonne. Essayons d'y accéder quand même.

Modifiez la fonction `separateMetadata` pour inclure un champ run_id :

=== "Après"

    ```groovy title="main.nf" linenums="5" hl_lines="9"
    def separateMetadata(row) {
        def sample_meta = [
            id: row.sample_id.toLowerCase(),
            organism: row.organism,
            tissue: row.tissue_type.replaceAll('_', ' ').toLowerCase(),
            depth: row.sequencing_depth.toInteger(),
            quality: row.quality_score.toDouble()
        ]
        def run_id = row.run_id.toUpperCase()
    ```

=== "Avant"

    ```groovy title="main.nf" linenums="5"
    def separateMetadata(row) {
        def sample_meta = [
            id: row.sample_id.toLowerCase(),
            organism: row.organism,
            tissue: row.tissue_type.replaceAll('_', ' ').toLowerCase(),
            depth: row.sequencing_depth.toInteger(),
            quality: row.quality_score.toDouble()
        ]
    ```

Maintenant, exécutez le workflow :

```bash
nextflow run main.nf
```

??? failure "Sortie de la commande"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `main.nf` [trusting_torvalds] DSL2 - revision: b56fbfbce2

    ERROR ~ Cannot invoke method toUpperCase() on null object

    -- Check script 'main.nf' at line: 13 or see '.nextflow.log' file for more details
    ```

Cela plante avec une NullPointerException.

Le problème est que `row.run_id` retourne `null` car la colonne `run_id` n'existe pas dans notre CSV. Lorsque nous essayons d'appeler `.toUpperCase()` sur `null`, cela plante. C'est là que l'opérateur de navigation sûre sauve la mise.

### 6.2. Opérateur de navigation sûre (`?.`)

L'opérateur de navigation sûre (`?.`) retourne `null` au lieu de lever une exception lorsqu'il est appelé sur une valeur `null`. Si l'objet avant `?.` est `null`, l'expression entière s'évalue à `null` sans exécuter la méthode.

Mettez à jour la fonction pour utiliser la navigation sûre :

=== "Après"

    ```groovy title="main.nf" linenums="4" hl_lines="9"
    def separateMetadata(row) {
        def sample_meta = [
            id: row.sample_id.toLowerCase(),
            organism: row.organism,
            tissue: row.tissue_type.replaceAll('_', ' ').toLowerCase(),
            depth: row.sequencing_depth.toInteger(),
            quality: row.quality_score.toDouble()
        ]
        def run_id = row.run_id?.toUpperCase()
    ```

=== "Avant"

    ```groovy title="main.nf" linenums="4" hl_lines="9"
    def separateMetadata(row) {
        def sample_meta = [
            id: row.sample_id.toLowerCase(),
            organism: row.organism,
            tissue: row.tissue_type.replaceAll('_', ' ').toLowerCase(),
            depth: row.sequencing_depth.toInteger(),
            quality: row.quality_score.toDouble()
        ]
        def run_id = row.run_id.toUpperCase()
    ```

Exécutez à nouveau :

```bash
nextflow run main.nf
```

??? success "Sortie de la commande"

    ```console
    <!-- TODO: output -->
    ```

Pas de plantage ! Le workflow gère maintenant le champ manquant avec élégance. Lorsque `row.run_id` est `null`, l'opérateur `?.` empêche l'appel `.toUpperCase()`, et `run_id` devient `null` au lieu de causer une exception.

### 6.3. Opérateur Elvis (`?:`) pour les valeurs par défaut

L'opérateur Elvis (`?:`) fournit des valeurs par défaut lorsque le côté gauche est "faux" (comme expliqué précédemment). Il est nommé d'après Elvis Presley car `?:` ressemble à ses cheveux et ses yeux célèbres vus de côté !

Maintenant que nous utilisons la navigation sûre, `run_id` sera `null` pour les échantillons sans ce champ. Utilisons l'opérateur Elvis pour fournir une valeur par défaut et l'ajouter à notre map `sample_meta` :

=== "Après"

    ```groovy title="main.nf" linenums="5" hl_lines="9-10"
    def separateMetadata(row) {
        def sample_meta = [
            id: row.sample_id.toLowerCase(),
            organism: row.organism,
            tissue: row.tissue_type.replaceAll('_', ' ').toLowerCase(),
            depth: row.sequencing_depth.toInteger(),
            quality: row.quality_score.toDouble()
        ]
        def run_id = row.run_id?.toUpperCase() ?: 'UNSPECIFIED'
        sample_meta.run = run_id
    ```

=== "Avant"

    ```groovy title="main.nf" linenums="5" hl_lines="9"
    def separateMetadata(row) {
        def sample_meta = [
            id: row.sample_id.toLowerCase(),
            organism: row.organism,
            tissue: row.tissue_type.replaceAll('_', ' ').toLowerCase(),
            depth: row.sequencing_depth.toInteger(),
            quality: row.quality_score.toDouble()
        ]
        def run_id = row.run_id?.toUpperCase()
    ```

Ajoutez également un opérateur `view()` dans le workflow pour voir les résultats :

=== "Après"

    ```groovy title="main.nf" linenums="30" hl_lines="4"
        ch_samples = channel.fromPath("./data/samples.csv")
            .splitCsv(header: true)
            .map{ row -> separateMetadata(row) }
            .view()
    ```

=== "Avant"

    ```groovy title="main.nf" linenums="30"
        ch_samples = channel.fromPath("./data/samples.csv")
            .splitCsv(header: true)
            .map{ row -> separateMetadata(row) }
    ```

et exécutez le workflow :

```bash
nextflow run main.nf
```

??? success "Sortie de la commande"

    ```console
    [[id:sample_001, organism:human, tissue:liver, depth:30000000, quality:38.5, run:UNSPECIFIED, sample_num:1, lane:001, read:R1, chunk:001, priority:normal], /workspaces/training/side-quests/essential_scripting_patterns/data/sequences/SAMPLE_001_S1_L001_R1_001.fastq]
    [[id:sample_002, organism:mouse, tissue:brain, depth:25000000, quality:35.2, run:UNSPECIFIED, sample_num:2, lane:001, read:R1, chunk:001, priority:normal], /workspaces/training/side-quests/essential_scripting_patterns/data/sequences/SAMPLE_002_S2_L001_R1_001.fastq]
    [[id:sample_003, organism:human, tissue:kidney, depth:45000000, quality:42.1, run:UNSPECIFIED, sample_num:3, lane:001, read:R1, chunk:001, priority:high], /workspaces/training/side-quests/essential_scripting_patterns/data/sequences/SAMPLE_003_S3_L001_R1_001.fastq]
    ```

Parfait ! Maintenant, tous les échantillons ont un champ `run` avec soit leur ID de séquençage réel (en majuscules) soit la valeur par défaut 'UNSPECIFIED'. La combinaison de `?.` et `?:` fournit à la fois la sécurité (pas de plantages) et des valeurs par défaut sensées.

Retirez maintenant l'opérateur `.view()` maintenant que nous avons confirmé que cela fonctionne.

!!! tip "Combiner navigation sûre et Elvis"

    Le modèle `value?.method() ?: 'default'` est courant dans les workflows de production :

    - `value?.method()` - Appelle la méthode en toute sécurité, retourne `null` si `value` est `null`
    - `?: 'default'` - Fournit une valeur de secours si le résultat est `null`

    Ce modèle gère les données manquantes/incomplètes avec élégance.

Utilisez ces opérateurs de manière cohérente dans les fonctions, les closures d'opérateurs (`.map{}`, `.filter{}`), les scripts de processus et les fichiers de configuration. Ils empêchent les plantages lors de la gestion de données du monde réel.

### À retenir

- **Navigation sûre (`?.`)** : Empêche les plantages sur les valeurs null - retourne null au lieu de lever une exception
- **Opérateur Elvis (`?:`)** : Fournit des valeurs par défaut - `value ?: 'default'`
- **Combinaison** : `value?.method() ?: 'default'` est le modèle courant

Ces opérateurs rendent les workflows résilients aux données incomplètes - essentiel pour le travail du monde réel.

---

## 7. Validation avec `error()` et `log.warn`

Parfois, vous devez arrêter immédiatement le workflow si les paramètres d'entrée sont invalides. Dans Nextflow, vous pouvez utiliser des fonctions intégrées comme `error()` et `log.warn`, ainsi que des constructions de programmation standard comme les instructions `if` et la logique booléenne, pour implémenter une logique de validation. Ajoutons une validation à notre workflow.

Créez une fonction de validation avant votre bloc workflow, appelez-la depuis le workflow, et changez la création du canal pour utiliser un paramètre pour le chemin du fichier CSV. Si le paramètre est manquant ou si le fichier n'existe pas, appelez `error()` pour arrêter l'exécution avec un message clair.

=== "Après"

    ```groovy title="main.nf" linenums="1" hl_lines="5-20 23-24"
    include { FASTP } from './modules/fastp.nf'
    include { TRIMGALORE } from './modules/trimgalore.nf'
    include { GENERATE_REPORT } from './modules/generate_report.nf'

    def validateInputs() {
        // Vérifier que le paramètre d'entrée est fourni
        if (!params.input) {
            error("Chemin du fichier CSV d'entrée non fourni. Veuillez spécifier --input <file.csv>")
        }

        // Vérifier que le fichier CSV existe
        if (!file(params.input).exists()) {
            error("Fichier CSV d'entrée introuvable : ${params.input}")
        }
    }
    ...
    workflow {
        validateInputs()
        ch_samples = channel.fromPath(params.input)
    ```

=== "Avant"

    ```groovy title="main.nf" linenums="1"
    include { FASTP } from './modules/fastp.nf'
    include { TRIMGALORE } from './modules/trimgalore.nf'
    include { GENERATE_REPORT } from './modules/generate_report.nf'

    ...
    workflow {
        ch_samples = channel.fromPath("./data/samples.csv")
    ```

Maintenant, essayez d'exécuter sans le fichier CSV :

```bash
nextflow run main.nf
```

??? failure "Sortie de la commande"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `main.nf` [confident_coulomb] DSL2 - revision: 07059399ed

    WARN: Access to undefined parameter `input` -- Initialise it to a default value eg. `params.input = some_value`
    Chemin du fichier CSV d'entrée non fourni. Veuillez spécifier --input <file.csv>
    ```

Le workflow s'arrête immédiatement avec un message d'erreur clair au lieu d'échouer mystérieusement plus tard

Maintenant, exécutez-le avec un fichier inexistant :

```bash
nextflow run main.nf --input ./data/nonexistent.csv
```

??? failure "Sortie de la commande"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `main.nf` [cranky_gates] DSL2 - revision: 26839ae3eb

    Fichier CSV d'entrée introuvable : ./data/nonexistent.csv
    ```

Enfin, exécutez-le avec le fichier correct :

```bash
nextflow run main.nf --input ./data/samples.csv
```

??? success "Sortie de la commande"

    ```console
    <!-- TODO: output -->
    ```

Cette fois, il s'exécute avec succès.

Vous pouvez également ajouter une validation dans la fonction `separateMetadata`. Utilisons le non-fatal `log.warn` pour émettre des avertissements pour les échantillons avec une faible profondeur de séquençage, mais permettons quand même au workflow de continuer :

=== "Après"

    ```groovy title="main.nf" linenums="1" hl_lines="3-6"
        def priority = sample_meta.quality > 40 ? 'high' : 'normal'

        // Valider que les données ont du sens
        if (sample_meta.depth < 30000000) {
            log.warn "Faible profondeur de séquençage pour ${sample_meta.id} : ${sample_meta.depth}"
        }

        return tuple(sample_meta + file_meta + [priority: priority], fastq_path)
    }
    ```

=== "Avant"

    ```groovy title="main.nf" linenums="1"
        def priority = sample_meta.quality > 40 ? 'high' : 'normal'

        return tuple(sample_meta + file_meta + [priority: priority], fastq_path)
    }
    ```

Exécutez à nouveau le workflow avec le CSV d'origine :

```bash
nextflow run main.nf --input ./data/samples.csv
```

??? warning "Sortie de la commande"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `main.nf` [awesome_goldwasser] DSL2 - revision: a31662a7c1

    executor >  local (5)
    [ce/df5eeb] process > FASTP (2)           [100%] 2 of 2 ✔
    [-        ] process > TRIMGALORE          -
    [d1/7d2b4b] process > GENERATE_REPORT (3) [100%] 3 of 3 ✔
    WARN: Faible profondeur de séquençage pour sample_002 : 25000000
    ```

Nous voyons un avertissement concernant la faible profondeur de séquençage pour l'un des échantillons.

### À retenir

- **`error()`** : Arrête immédiatement le workflow avec un message clair
- **`log.warn`** : Émet des avertissements sans arrêter le workflow
- **Validation précoce** : Vérifier les entrées avant le traitement pour échouer rapidement avec des erreurs utiles
- **Fonctions de validation** : Créer une logique de validation réutilisable qui peut être appelée au démarrage du workflow

Une validation appropriée rend les workflows plus robustes et conviviaux en détectant les problèmes tôt avec des messages d'erreur clairs.

---

## 8. Gestionnaires d'événements de workflow

Jusqu'à présent, nous avons écrit du code dans nos scripts de workflow et nos définitions de processus. Mais il y a une autre fonctionnalité importante que vous devriez connaître : les gestionnaires d'événements de workflow.

Les gestionnaires d'événements sont des closures qui s'exécutent à des points spécifiques du cycle de vie de votre workflow. Ils sont parfaits pour ajouter de la journalisation, des notifications ou des opérations de nettoyage. Ces gestionnaires doivent être définis dans votre script de workflow aux côtés de votre définition de workflow.

### 8.1. Le gestionnaire `onComplete`

Le gestionnaire d'événements le plus couramment utilisé est `onComplete`, qui s'exécute lorsque votre workflow se termine (qu'il ait réussi ou échoué). Ajoutons-en un pour résumer les résultats de notre pipeline.

Ajoutez le gestionnaire d'événements à votre fichier `main.nf`, à l'intérieur de votre définition de workflow :

=== "Après"

    ```groovy title="main.nf" linenums="66" hl_lines="5-16"
        ch_fastp = FASTP(trim_branches.fastp)
        ch_trimgalore = TRIMGALORE(trim_branches.trimgalore)
        GENERATE_REPORT(ch_samples)

        workflow.onComplete = {
            println ""
            println "Résumé de l'exécution du pipeline :"
            println "=========================="
            println "Terminé à : ${workflow.complete}"
            println "Durée     : ${workflow.duration}"
            println "Succès    : ${workflow.success}"
            println "workDir   : ${workflow.workDir}"
            println "statut de sortie : ${workflow.exitStatus}"
            println ""
        }
    }
    ```

=== "Avant"

    ```groovy title="main.nf" linenums="66" hl_lines="4"
        ch_fastp = FASTP(trim_branches.fastp)
        ch_trimgalore = TRIMGALORE(trim_branches.trimgalore)
        GENERATE_REPORT(ch_samples)
    }
    ```

Cette closure s'exécute lorsque le workflow se termine. À l'intérieur, vous avez accès à l'objet `workflow` qui fournit des propriétés utiles sur l'exécution.

Exécutez votre workflow et vous verrez ce résumé apparaître à la fin !

```bash
nextflow run main.nf --input ./data/samples.csv -ansi-log false
```

??? success "Sortie de la commande"

    ```console
    N E X T F L O W  ~  version 25.10.2
    Launching `main.nf` [marvelous_boltzmann] DSL2 - revision: a31662a7c1
    WARN: Faible profondeur de séquençage pour sample_002 : 25000000
    [9b/d48e40] Submitted process > FASTP (2)
    [6a/73867a] Submitted process > GENERATE_REPORT (2)
    [79/ad0ac5] Submitted process > GENERATE_REPORT (1)
    [f3/bda6cb] Submitted process > FASTP (1)
    [34/d5b52f] Submitted process > GENERATE_REPORT (3)

    Résumé de l'exécution du pipeline :
    ==========================
    Terminé à : 2025-10-10T12:14:24.885384+01:00
    Durée     : 2.9s
    Succès    : true
    workDir   : /workspaces/training/side-quests/essential_scripting_patterns/work
    statut de sortie : 0
    ```

Rendons-le plus utile en ajoutant une logique conditionnelle :

=== "Après"

    ```groovy title="main.nf" linenums="66" hl_lines="5-22"
        ch_fastp = FASTP(trim_branches.fastp)
        ch_trimgalore = TRIMGALORE(trim_branches.trimgalore)
        GENERATE_REPORT(ch_samples)

        workflow.onComplete = {
            println ""
            println "Résumé de l'exécution du pipeline :"
            println "=========================="
            println "Terminé à : ${workflow.complete}"
            println "Durée     : ${workflow.duration}"
            println "Succès    : ${workflow.success}"
            println "workDir   : ${workflow.workDir}"
            println "statut de sortie : ${workflow.exitStatus}"
            println ""

            if (workflow.success) {
                println "✅ Pipeline terminé avec succès !"
            } else {
                println "❌ Le pipeline a échoué !"
                println "Erreur : ${workflow.errorMessage}"
            }
        }
    }
    ```

=== "Avant"

    ```groovy title="main.nf" linenums="66" hl_lines="5-16"
        ch_fastp = FASTP(trim_branches.fastp)
        ch_trimgalore = TRIMGALORE(trim_branches.trimgalore)
        GENERATE_REPORT(ch_samples)

        workflow.onComplete = {
            println ""
            println "Résumé de l'exécution du pipeline :"
            println "=========================="
            println "Terminé à : ${workflow.complete}"
            println "Durée     : ${workflow.duration}"
            println "Succès    : ${workflow.success}"
            println "workDir   : ${workflow.workDir}"
            println "statut de sortie : ${workflow.exitStatus}"
            println ""
        }
    }
    ```

Maintenant, nous obtenons un résumé encore plus informatif, incluant un message de succès/échec et le répertoire de sortie s'il est spécifié :

<!-- TODO: add run command -->

??? success "Sortie de la commande"

    ```console
    N E X T F L O W  ~  version 25.10.2
    Launching `main.nf` [boring_linnaeus] DSL2 - revision: a31662a7c1
    WARN: Faible profondeur de séquençage pour sample_002 : 25000000
    [e5/242efc] Submitted process > FASTP (2)
    [3b/74047c] Submitted process > GENERATE_REPORT (3)
    [8a/7a57e6] Submitted process > GENERATE_REPORT (1)
    [a8/b1a31f] Submitted process > GENERATE_REPORT (2)
    [40/648429] Submitted process > FASTP (1)

    Résumé de l'exécution du pipeline :
    ==========================
    Terminé à : 2025-10-10T12:16:00.522569+01:00
    Durée     : 3.6s
    Succès    : true
    workDir   : /workspaces/training/side-quests/essential_scripting_patterns/work
    statut de sortie : 0

    ✅ Pipeline terminé avec succès !
    ```

Vous pouvez également écrire le résumé dans un fichier en utilisant des opérations sur les fichiers :

```groovy title="main.nf - Écrire le résumé dans un fichier"
workflow {
    // ... votre code de workflow ...

    workflow.onComplete = {
        def summary = """
        Résumé de l'exécution du pipeline
        ===========================
        Terminé : ${workflow.complete}
        Durée   : ${workflow.duration}
        Succès  : ${workflow.success}
        Commande : ${workflow.commandLine}
        """

        println summary

        // Écrire dans un fichier journal
        def log_file = file("${workflow.launchDir}/pipeline_summary.txt")
        log_file.text = summary
    }
}
```

### 8.2. Le gestionnaire `onError`

Outre `onComplete`, il existe un autre gestionnaire d'événements que vous pouvez utiliser : `onError`, qui s'exécute uniquement si le workflow échoue :

```groovy title="main.nf - gestionnaire onError"
workflow {
    // ... votre code de workflow ...

    workflow.onError = {
        println "="* 50
        println "L'exécution du pipeline a échoué !"
        println "Message d'erreur : ${workflow.errorMessage}"
        println "="* 50

        // Écrire un journal d'erreur détaillé
        def error_file = file("${workflow.launchDir}/error.log")
        error_file.text = """
        Rapport d'erreur du workflow
        =====================
        Heure : ${new Date()}
        Erreur : ${workflow.errorMessage}
        Rapport d'erreur : ${workflow.errorReport ?: 'Aucun rapport détaillé disponible'}
        """

        println "Détails de l'erreur écrits dans : ${error_file}"
    }
}
```

Vous pouvez utiliser plusieurs gestionnaires ensemble dans votre script de workflow :

```groovy title="main.nf - Gestionnaires combinés"
workflow {
    // ... votre code de workflow ...

    workflow.onError = {
        println "Le workflow a échoué : ${workflow.errorMessage}"
    }

    workflow.onComplete = {
        def duration_mins = workflow.duration.toMinutes().round(2)
        def status = workflow.success ? "SUCCÈS ✅" : "ÉCHEC ❌"

        println """
        Pipeline terminé : ${status}
        Durée : ${duration_mins} minutes
        """
    }
}
```

### À retenir

Dans cette section, vous avez appris :

- **Closures de gestionnaires d'événements** : Closures dans votre script de workflow qui s'exécutent à différents points du cycle de vie
- **Gestionnaire `onComplete`** : Pour les résumés d'exécution et les rapports de résultats
- **Gestionnaire `onError`** : Pour la gestion des erreurs et la journalisation des échecs
- **Propriétés de l'objet workflow** : Accès à `workflow.success`, `workflow.duration`, `workflow.errorMessage`, etc.

Les gestionnaires d'événements montrent comment vous pouvez utiliser toute la puissance du langage Nextflow dans vos scripts de workflow pour ajouter des capacités sophistiquées de journalisation et de notification.

---

## Résumé

Félicitations, vous avez réussi !

Tout au long de cette quête secondaire, vous avez construit un pipeline complet de traitement d'échantillons qui a évolué d'une gestion de métadonnées de base à un workflow sophistiqué et prêt pour la production.
Chaque section s'est appuyée sur la précédente, démontrant comment les constructions de programmation transforment des workflows simples en systèmes puissants de traitement de données, avec les avantages suivants :

- **Code plus clair** : Comprendre le dataflow vs le scripting vous aide à écrire des workflows plus organisés
- **Gestion robuste** : La navigation sûre et les opérateurs Elvis rendent les workflows résilients aux données manquantes
- **Traitement flexible** : La logique conditionnelle permet à vos workflows de traiter différents types d'échantillons de manière appropriée
- **Ressources adaptatives** : Les directives dynamiques optimisent l'utilisation des ressources en fonction des caractéristiques d'entrée

Cette progression reflète l'évolution réelle des pipelines bioinformatiques, des prototypes de recherche traitant quelques échantillons aux systèmes de production traitant des milliers d'échantillons à travers les laboratoires et les institutions.
Chaque défi que vous avez résolu et chaque modèle que vous avez appris reflète des problèmes réels auxquels les développeur·ses sont confronté·es lors de la mise à l'échelle des workflows Nextflow.

L'application de ces modèles dans votre propre travail vous permettra de construire des workflows robustes et prêts pour la production.

### Modèles clés

1.  **Dataflow vs Scripting :** Vous avez appris à distinguer les opérations de dataflow (orchestration de canaux) et le scripting (code qui manipule les données), y compris les différences cruciales entre les opérations sur différents types comme `collect` sur Channel vs List.

    - Dataflow : orchestration de canaux

    ```groovy
    channel.fromPath('*.fastq').splitCsv(header: true)
    ```

    - Scripting : traitement de données sur les collections

    ```groovy
    sample_data.collect { it.toUpperCase() }
    ```

2.  **Traitement avancé de chaînes** : Vous avez maîtrisé les expressions régulières pour analyser les noms de fichiers, la génération de scripts dynamiques dans les processus, et l'interpolation de variables (Nextflow vs Bash vs Shell).

    - Correspondance de motifs

    ```groovy
    filename =~ ~/^(\w+)_(\w+)_(\d+)\.fastq$/
    ```

    - Fonction avec retour conditionnel

    ```groovy
    def parseSample(filename) {
        def matcher = filename =~ pattern
        return matcher ? [valid: true, data: matcher[0]] : [valid: false]
    }
    ```

    - Collection de fichiers vers arguments de commande (dans le bloc script de processus)

    ```groovy
    script:
    def file_args = input_files.collect { file -> "--input ${file}" }.join(' ')
    """
    analysis_tool ${file_args} --output results.txt
    """
    ```

3.  **Création de fonctions réutilisables** : Vous avez appris à extraire la logique complexe dans des fonctions nommées qui peuvent être appelées depuis les opérateurs de canal, rendant les workflows plus lisibles et maintenables.

    - Définir une fonction nommée

    ```groovy
    def separateMetadata(row) {
        def sample_meta = [ /* code caché pour la brièveté */ ]
        def fastq_path = file(row.file_path)
        def m = (fastq_path.name =~ /^(.+)_S(\d+)_L(\d{3})_(R[12])_(\d{3})\.fastq(?:\.gz)?$/)
        def file_meta = m ? [ /* code caché pour la brièveté */ ] : [:]
        def priority = sample_meta.quality > 40 ? 'high' : 'normal'

        return tuple(sample_meta + file_meta + [priority: priority], fastq_path)
    }
    ```

    - Appeler la fonction nommée dans un workflow

    ```groovy
    workflow {
        ch_samples = channel.fromPath("./data/samples.csv")
            .splitCsv(header: true)
            .map{ row -> separateMetadata(row) }

        ch_fastp = FASTP(ch_samples)
    }
    ```

4.  **Directives de ressources dynamiques avec des closures** : Vous avez exploré l'utilisation de closures dans les directives de processus pour une allocation de ressources adaptative basée sur les caractéristiques d'entrée.

    - Closures nommées et composition

    ```groovy
    def enrichData = normalizeId >> addQualityCategory >> addFlags
    def processor = generalFunction.curry(fixedParam)
    ```

    - Closures avec accès à la portée

    ```groovy
    def collectStats = { data -> stats.count++; return data }
    ```

5.  **Logique conditionnelle et contrôle de processus** : Vous avez ajouté un routage intelligent en utilisant les opérateurs `.branch()` et `.filter()`, en tirant parti de la véracité pour des expressions conditionnelles concises.

    - Utiliser `.branch()` pour router les données à travers différentes branches de workflow

    ```groovy
    trim_branches = ch_samples
    .branch { meta, reads ->
        fastp: meta.organism == 'human' && meta.depth >= 30000000
        trimgalore: true
    }

    ch_fastp = FASTP(trim_branches.fastp)
    ch_trimgalore = TRIMGALORE(trim_branches.trimgalore)
    ```

    - Évaluation booléenne avec la véracité Groovy

    ```groovy
    if (sample.files) println "A des fichiers"
    ```

    - Utiliser `filter()` pour sous-ensembler les données avec la 'véracité'

    ````groovy
    ch_valid_samples = ch```groovy
    ch_valid_samples = ch_samples
        .filter { meta, reads ->
            meta.id && meta.organism && meta.depth >= 25000000
        }
    ````

6.  **Navigation sûre et opérateurs Elvis** : Vous avez rendu le pipeline robuste contre les données manquantes en utilisant `?.` pour l'accès sûr aux propriétés et `?:` pour fournir des valeurs par défaut.

    ```groovy
    def id = data?.sample?.id ?: 'unknown'
    ```

7.  **Validation avec error() et log.warn** : Vous avez appris à valider les entrées tôt et à échouer rapidement avec des messages d'erreur clairs.

    ```groovy
    try {
        def errors = validateSample(sample)
        if (errors) throw new RuntimeException("Invalide : ${errors.join(', ')}")
    } catch (Exception e) {
        println "Erreur : ${e.message}"
    }
    ```

8.  **Gestionnaires d'événements de configuration** : Vous avez appris à utiliser les gestionnaires d'événements de workflow (`onComplete` et `onError`) pour la journalisation, les notifications et la gestion du cycle de vie.

    - Utiliser `onComplete` pour journaliser et notifier

    ```groovy
    workflow.onComplete = {
        println "Succès     : ${workflow.success}"
        println "statut de sortie : ${workflow.exitStatus}"

        if (workflow.success) {
            println "✅ Pipeline terminé avec succès !"
        } else {
            println "❌ Le pipeline a échoué !"
            println "Erreur : ${workflow.errorMessage}"
        }
    }
    ```

    - Utiliser `onError` pour prendre des mesures spécifiquement en cas d'échec

    ```groovy
    workflow.onError = {
        // Écrire un journal d'erreur détaillé
        def error_file = file("${workflow.launchDir}/error.log")
        error_file.text = """
        Heure : ${new Date()}
        Erreur : ${workflow.errorMessage}
        Rapport d'erreur : ${workflow.errorReport ?: 'Aucun rapport détaillé disponible'}
        """

        println "Détails de l'erreur écrits dans : ${error_file}"
    }
    ```

### Ressources supplémentaires

- [Référence du langage Nextflow](https://nextflow.io/docs/latest/reference/syntax.html)
- [Opérateurs Nextflow](https://www.nextflow.io/docs/latest/operator.html)
- [Syntaxe de script Nextflow](https://www.nextflow.io/docs/latest/script.html)
- [Bibliothèque standard Nextflow](https://nextflow.io/docs/latest/reference/stdlib.html)

Assurez-vous de consulter ces ressources lorsque vous devez explorer des fonctionnalités plus avancées.

Vous bénéficierez de la pratique et de l'expansion de vos compétences afin de :

- Écrire des workflows plus propres avec une séparation appropriée entre le dataflow et le scripting
- Maîtriser l'interpolation de variables pour éviter les pièges courants avec les variables Nextflow, Bash et shell
- Utiliser des directives de ressources dynamiques pour des workflows efficaces et adaptatifs
- Transformer des collections de fichiers en arguments de ligne de commande correctement formatés
- Gérer différentes conventions de nommage de fichiers et formats d'entrée avec élégance en utilisant regex et le traitement de chaînes
- Construire du code réutilisable et maintenable en utilisant des modèles de closure avancés et la programmation fonctionnelle
- Traiter et organiser des ensembles de données complexes en utilisant des opérations sur les collections
- Ajouter de la validation, de la gestion des erreurs et de la journalisation pour rendre vos workflows prêts pour la production
- Implémenter la gestion du cycle de vie du workflow avec des gestionnaires d'événements

---

## Et ensuite ?

Retournez au [menu des Quêtes secondaires](./index.md) ou cliquez sur le bouton en bas à droite de la page pour passer au sujet suivant de la liste.
