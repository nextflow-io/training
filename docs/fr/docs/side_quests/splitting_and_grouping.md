# Division et regroupement

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Traduction assistée par IA - [en savoir plus et suggérer des améliorations](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

Nextflow fournit des outils puissants pour travailler avec les données de manière flexible. Une capacité clé est la division des données en différents flux, puis le regroupement d'éléments connexes. Ceci est particulièrement précieux dans les workflows de bioinformatique où vous devez traiter différents types d'échantillons séparément avant de combiner les résultats pour l'analyse.

Imaginez cela comme le tri du courrier : vous séparez les lettres par destination, traitez chaque pile différemment, puis recombinez les éléments allant à la même personne. Nextflow utilise des opérateurs spéciaux pour accomplir cela avec des données scientifiques. Cette approche est également communément appelée le modèle **scatter/gather** dans l'informatique distribuée et les workflows de bioinformatique.

Le système de canaux de Nextflow est au cœur de cette flexibilité. Les canaux connectent différentes parties de votre workflow, permettant aux données de circuler à travers votre analyse. Vous pouvez créer plusieurs canaux à partir d'une seule source de données, traiter chaque canal différemment, puis fusionner les canaux lorsque nécessaire. Cette approche vous permet de concevoir des workflows qui reflètent naturellement les chemins de branchement et de convergence d'analyses bioinformatiques complexes.

### Objectifs d'apprentissage

Dans cette quête annexe, vous apprendrez à diviser et regrouper des données en utilisant les opérateurs de canaux de Nextflow.
Nous commencerons avec un fichier CSV contenant des informations sur les échantillons et les fichiers de données associés, puis manipulerons et réorganiserons ces données.

À la fin de cette quête annexe, vous serez capable de séparer et combiner efficacement des flux de données, en utilisant les techniques suivantes :

- Lire des données à partir de fichiers en utilisant `splitCsv`
- Filtrer et transformer des données avec `filter` et `map`
- Combiner des données connexes en utilisant `join` et `groupTuple`
- Créer des combinaisons de données avec `combine` pour un traitement parallèle
- Optimiser la structure des données en utilisant `subMap` et des stratégies de déduplication
- Construire des fonctions réutilisables avec des closures nommées pour vous aider à manipuler les structures de canaux

Ces compétences vous aideront à construire des workflows capables de gérer efficacement plusieurs fichiers d'entrée et différents types de données, tout en maintenant une structure de code propre et maintenable.

### Prérequis

Avant d'entreprendre cette quête annexe, vous devez :

- Avoir terminé le tutoriel [Hello Nextflow](../hello_nextflow/README.md) ou un cours équivalent pour débutants.
- Être à l'aise avec les concepts et mécanismes de base de Nextflow (processes, channels, operators, travailler avec des fichiers, méta données)

**Optionnel :** Nous recommandons de compléter d'abord la quête annexe [Metadata in workflows](./metadata.md).
Elle couvre les fondamentaux de la lecture de fichiers CSV avec `splitCsv` et de la création de meta maps, que nous utiliserons beaucoup ici.

---

## 0. Commencer

#### Ouvrir le codespace de formation

Si vous ne l'avez pas encore fait, assurez-vous d'ouvrir l'environnement de formation comme décrit dans la [Configuration de l'environnement](../envsetup/index.md).

[![Open in GitHub Codespaces](https://github.com/codespaces/badge.svg)](https://codespaces.new/nextflow-io/training?quickstart=1&ref=master)

#### Se déplacer dans le répertoire du projet

Déplaçons-nous dans le répertoire où se trouvent les fichiers pour ce tutoriel.

```bash
cd side-quests/splitting_and_grouping
```

Vous pouvez configurer VSCode pour se concentrer sur ce répertoire :

```bash
code .
```

#### Examiner les matériaux

Vous trouverez un fichier de workflow principal et un répertoire `data` contenant une feuille d'échantillons nommée `samplesheet.csv`.

```console title="Contenu du répertoire"
.
├── data
│   └── samplesheet.csv
└── main.nf
```

La feuille d'échantillons contient des informations sur des échantillons de différents patients, incluant l'identifiant du patient, le numéro de répétition de l'échantillon, le type (normal ou tumeur), et les chemins vers des fichiers de données hypothétiques (qui n'existent pas réellement, mais nous ferons comme s'ils existaient).

```console title="samplesheet.csv"
id,repeat,type,bam
patientA,1,normal,patientA_rep1_normal.bam
patientA,1,tumor,patientA_rep1_tumor.bam
patientA,2,normal,patientA_rep2_normal.bam
patientA,2,tumor,patientA_rep2_tumor.bam
patientB,1,normal,patientB_rep1_normal.bam
patientB,1,tumor,patientB_rep1_tumor.bam
patientC,1,normal,patientC_rep1_normal.bam
patientC,1,tumor,patientC_rep1_tumor.bam
```

Cette feuille d'échantillons répertorie huit échantillons provenant de trois patients (A, B, C).

Pour chaque patient, nous avons des échantillons de type `tumor` (provenant généralement de biopsies tumorales) ou `normal` (prélevés sur du tissu sain ou du sang).
Si vous n'êtes pas familier avec l'analyse du cancer, sachez simplement que cela correspond à un modèle expérimental qui utilise des échantillons appariés tumeur/normal pour effectuer des analyses contrastives.

Pour le patient A spécifiquement, nous avons deux ensembles de réplicats techniques (répétitions).

!!! note

    Ne vous inquiétez pas si vous n'êtes pas familier avec cette conception expérimentale, ce n'est pas critique pour comprendre ce tutoriel.

#### Examiner l'assignation

Votre défi est d'écrire un workflow Nextflow qui va :

1. **Lire** les données d'échantillons à partir d'un fichier CSV et les structurer avec des meta maps
2. **Séparer** les échantillons dans différents canaux en fonction du type (normal vs tumeur)
3. **Joindre** les paires appariées tumeur/normal par identifiant de patient et numéro de réplicat
4. **Distribuer** les échantillons sur des intervalles génomiques pour un traitement parallèle
5. **Regrouper** les échantillons connexes pour une analyse en aval

Ceci représente un modèle bioinformatique courant où vous devez diviser les données pour un traitement indépendant, puis recombiner les éléments connexes pour une analyse comparative.

#### Liste de vérification de préparation

Pensez-vous être prêt·e à vous lancer ?

- [ ] Je comprends l'objectif de ce cours et ses prérequis
- [ ] Mon codespace est opérationnel
- [ ] J'ai configuré mon répertoire de travail de manière appropriée
- [ ] Je comprends l'assignation

Si vous pouvez cocher toutes les cases, vous êtes prêt·e à commencer.

---

## 1. Lire les données d'échantillons

### 1.1. Lire les données d'échantillons avec `splitCsv` et créer des meta maps

Commençons par lire les données d'échantillons avec `splitCsv` et les organiser selon le modèle de meta map. Dans le fichier `main.nf`, vous verrez que nous avons déjà commencé le workflow.

```groovy title="main.nf" linenums="1" hl_lines="2"
workflow {
    ch_samplesheet = channel.fromPath("./data/samplesheet.csv")
}
```

!!! note

    Tout au long de ce tutoriel, nous utiliserons le préfixe `ch_` pour toutes les variables de canaux afin d'indiquer clairement qu'il s'agit de canaux Nextflow.

Si vous avez complété la quête annexe [Metadata in workflows](./metadata.md), vous reconnaîtrez ce modèle. Nous utiliserons `splitCsv` pour lire le CSV et structurer immédiatement les données avec une meta map pour séparer les métadonnées des chemins de fichiers.

!!! info

    Nous rencontrerons deux concepts différents appelés `map` dans cette formation :

    - **Structure de données** : La map Groovy (équivalente aux dictionnaires/hashes dans d'autres langages) qui stocke des paires clé-valeur
    - **Opérateur de canal** : L'opérateur `.map()` qui transforme les éléments dans un canal

    Nous clarifierons de quel type nous parlons dans le contexte, mais cette distinction est importante à comprendre lors du travail avec Nextflow.

Appliquez ces modifications à `main.nf` :

=== "Après"

    ```groovy title="main.nf" linenums="2" hl_lines="2-6"
        ch_samples = channel.fromPath("./data/samplesheet.csv")
            .splitCsv(header: true)
            .map{ row ->
              [[id:row.id, repeat:row.repeat, type:row.type], row.bam]
            }
            .view()
    ```

=== "Avant"

    ```groovy title="main.nf" linenums="2" hl_lines="1"
        ch_samplesheet = channel.fromPath("./data/samplesheet.csv")
    ```

Cela combine l'opération `splitCsv` (lecture du CSV avec en-têtes) et l'opération `map` (structuration des données en tuples `[meta, file]`) en une seule étape. Appliquez cette modification et exécutez le pipeline :

```bash
nextflow run main.nf
```

??? success "Sortie de la commande"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `main.nf` [deadly_mercator] DSL2 - revision: bd6b0224e9

    [[id:patientA, repeat:1, type:normal], patientA_rep1_normal.bam]
    [[id:patientA, repeat:1, type:tumor], patientA_rep1_tumor.bam]
    [[id:patientA, repeat:2, type:normal], patientA_rep2_normal.bam]
    [[id:patientA, repeat:2, type:tumor], patientA_rep2_tumor.bam]
    [[id:patientB, repeat:1, type:normal], patientB_rep1_normal.bam]
    [[id:patientB, repeat:1, type:tumor], patientB_rep1_tumor.bam]
    [[id:patientC, repeat:1, type:normal], patientC_rep1_normal.bam]
    [[id:patientC, repeat:1, type:tumor], patientC_rep1_tumor.bam]
    ```

Nous avons maintenant un canal où chaque élément est un tuple `[meta, file]` - métadonnées séparées des chemins de fichiers. Cette structure nous permet de diviser et regrouper notre charge de travail en fonction des champs de métadonnées.

---

## 2. Filtrer et transformer les données

### 2.1. Filtrer les données avec `filter`

Nous pouvons utiliser l'[opérateur `filter`](https://www.nextflow.io/docs/latest/operator.html#filter) pour filtrer les données en fonction d'une condition. Disons que nous voulons traiter uniquement les échantillons normaux. Nous pouvons le faire en filtrant les données en fonction du champ `type`. Insérons cela avant l'opérateur `view`.

=== "Après"

    ```groovy title="main.nf" linenums="2" hl_lines="6"
        ch_samples = channel.fromPath("./data/samplesheet.csv")
            .splitCsv(header: true)
            .map{ row ->
              [[id:row.id, repeat:row.repeat, type:row.type], row.bam]
            }
            .filter { meta, file -> meta.type == 'normal' }
            .view()
    ```

=== "Avant"

    ```groovy title="main.nf" linenums="2"
        ch_samples = channel.fromPath("./data/samplesheet.csv")
            .splitCsv(header: true)
            .map{ row ->
              [[id:row.id, repeat:row.repeat, type:row.type], row.bam]
            }
            .view()
    ```

Exécutez à nouveau le workflow pour voir le résultat filtré :

```bash
nextflow run main.nf
```

??? success "Sortie de la commande"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `main.nf` [admiring_brown] DSL2 - revision: 194d61704d

    [[id:patientA, repeat:1, type:normal], patientA_rep1_normal.bam]
    [[id:patientA, repeat:2, type:normal], patientA_rep2_normal.bam]
    [[id:patientB, repeat:1, type:normal], patientB_rep1_normal.bam]
    [[id:patientC, repeat:1, type:normal], patientC_rep1_normal.bam]
    ```

Nous avons réussi à filtrer les données pour inclure uniquement les échantillons normaux. Récapitulons comment cela fonctionne.

L'opérateur `filter` prend une closure qui est appliquée à chaque élément du canal. Si la closure renvoie `true`, l'élément est inclus ; si elle renvoie `false`, l'élément est exclu.

Dans notre cas, nous voulons conserver uniquement les échantillons où `meta.type == 'normal'`. La closure utilise le tuple `meta,file` pour faire référence à chaque échantillon, accède au type d'échantillon avec `meta.type`, et vérifie s'il est égal à `'normal'`.

Ceci est accompli avec la closure unique que nous avons introduite ci-dessus :

```groovy title="main.nf" linenums="7"
    .filter { meta, file -> meta.type == 'normal' }
```

### 2.2. Créer des canaux filtrés séparés

Actuellement, nous appliquons le filtre au canal créé directement à partir du CSV, mais nous voulons le filtrer de plusieurs façons, réécrivons donc la logique pour créer un canal filtré séparé pour les échantillons normaux :

=== "Après"

    ```groovy title="main.nf" linenums="2" hl_lines="6 8"
        ch_samples = channel.fromPath("./data/samplesheet.csv")
            .splitCsv(header: true)
            .map{ row ->
                [[id:row.id, repeat:row.repeat, type:row.type], row.bam]
            }
        ch_normal_samples = ch_samples
            .filter { meta, file -> meta.type == 'normal' }
        ch_normal_samples
            .view()
    ```

=== "Avant"

    ```groovy title="main.nf" linenums="2"
        ch_samples = channel.fromPath("./data/samplesheet.csv")
            .splitCsv(header: true)
            .map{ row ->
              [[id:row.id, repeat:row.repeat, type:row.type], row.bam]
            }
            .filter { meta, file -> meta.type == 'normal' }
            .view()
    ```

Exécutez le pipeline pour voir les résultats :

```bash
nextflow run main.nf
```

??? success "Sortie de la commande"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `main.nf` [trusting_poisson] DSL2 - revision: 639186ee74

    [[id:patientA, repeat:1, type:normal], patientA_rep1_normal.bam]
    [[id:patientA, repeat:2, type:normal], patientA_rep2_normal.bam]
    [[id:patientB, repeat:1, type:normal], patientB_rep1_normal.bam]
    [[id:patientC, repeat:1, type:normal], patientC_rep1_normal.bam]
    ```

Nous avons réussi à filtrer les données et créé un canal séparé pour les échantillons normaux.

Créons également un canal filtré pour les échantillons tumoraux :

=== "Après"

    ```groovy title="main.nf" linenums="7" hl_lines="3-8"
        ch_normal_samples = ch_samples
            .filter { meta, file -> meta.type == 'normal' }
        ch_tumor_samples = ch_samples
            .filter { meta, file -> meta.type == 'tumor' }
        ch_normal_samples
            .view{'Normal sample: ' + it}
        ch_tumor_samples
            .view{'Tumor sample: ' + it}
    ```

=== "Avant"

    ```groovy title="main.nf" linenums="7" hl_lines="3 4"
        ch_normal_samples = ch_samples
            .filter { meta, file -> meta.type == 'normal' }
        ch_normal_samples
            .view()
    ```

```bash
nextflow run main.nf
```

??? success "Sortie de la commande"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `main.nf` [maniac_boltzmann] DSL2 - revision: 3636b6576b

    Tumor sample: [[id:patientA, repeat:1, type:tumor], patientA_rep1_tumor.bam]
    Tumor sample: [[id:patientA, repeat:2, type:tumor], patientA_rep2_tumor.bam]
    Normal sample: [[id:patientA, repeat:1, type:normal], patientA_rep1_normal.bam]
    Normal sample: [[id:patientA, repeat:2, type:normal], patientA_rep2_normal.bam]
    Normal sample: [[id:patientB, repeat:1, type:normal], patientB_rep1_normal.bam]
    Normal sample: [[id:patientC, repeat:1, type:normal], patientC_rep1_normal.bam]
    Tumor sample: [[id:patientB, repeat:1, type:tumor], patientB_rep1_tumor.bam]
    Tumor sample: [[id:patientC, repeat:1, type:tumor], patientC_rep1_tumor.bam]
    ```

Nous avons séparé les échantillons normaux et tumoraux dans deux canaux différents, et utilisé une closure fournie à `view()` pour les étiqueter différemment dans la sortie : `ch_tumor_samples.view{'Tumor sample: ' + it}`.

### Point clé

Dans cette section, vous avez appris :

- **Filtrer les données** : Comment filtrer les données avec `filter`
- **Diviser les données** : Comment diviser les données en différents canaux en fonction d'une condition
- **Visualiser les données** : Comment utiliser `view` pour afficher les données et étiqueter la sortie de différents canaux

Nous avons maintenant séparé les échantillons normaux et tumoraux dans deux canaux différents. Ensuite, nous joindrons les échantillons normaux et tumoraux sur le champ `id`.

---

## 3. Joindre des canaux par identifiants

Dans la section précédente, nous avons séparé les échantillons normaux et tumoraux dans deux canaux différents. Ceux-ci pourraient être traités indépendamment en utilisant des processes ou workflows spécifiques en fonction de leur type. Mais que se passe-t-il lorsque nous voulons comparer les échantillons normaux et tumoraux du même patient ? À ce stade, nous devons les joindre à nouveau en nous assurant de faire correspondre les échantillons en fonction de leur champ `id`.

Nextflow inclut de nombreuses méthodes pour combiner des canaux, mais dans ce cas, l'opérateur le plus approprié est [`join`](https://www.nextflow.io/docs/latest/operator.html#join). Si vous êtes familier·ère avec SQL, il agit comme l'opération `JOIN`, où nous spécifions la clé sur laquelle joindre et le type de jointure à effectuer.

### 3.1. Utiliser `map` et `join` pour combiner en fonction de l'identifiant patient

Si nous consultons la documentation de [`join`](https://www.nextflow.io/docs/latest/operator.html#join), nous pouvons voir que par défaut, elle joint deux canaux en fonction du premier élément de chaque tuple.

#### 3.1.1. Vérifier la structure des données

Si vous n'avez plus la sortie console disponible, exécutons le pipeline pour vérifier notre structure de données et voir comment nous devons la modifier pour joindre sur le champ `id`.

```bash
nextflow run main.nf
```

??? success "Sortie de la commande"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `main.nf` [maniac_boltzmann] DSL2 - revision: 3636b6576b

    Tumor sample: [[id:patientA, repeat:1, type:tumor], patientA_rep1_tumor.bam]
    Tumor sample: [[id:patientA, repeat:2, type:tumor], patientA_rep2_tumor.bam]
    Normal sample: [[id:patientA, repeat:1, type:normal], patientA_rep1_normal.bam]
    Normal sample: [[id:patientA, repeat:2, type:normal], patientA_rep2_normal.bam]
    Normal sample: [[id:patientB, repeat:1, type:normal], patientB_rep1_normal.bam]
    Normal sample: [[id:patientC, repeat:1, type:normal], patientC_rep1_normal.bam]
    Tumor sample: [[id:patientB, repeat:1, type:tumor], patientB_rep1_tumor.bam]
    Tumor sample: [[id:patientC, repeat:1, type:tumor], patientC_rep1_tumor.bam]
    ```

Nous pouvons voir que le champ `id` est le premier élément dans chaque meta map. Pour que `join` fonctionne, nous devons isoler le champ `id` dans chaque tuple. Après cela, nous pouvons simplement utiliser l'opérateur `join` pour combiner les deux canaux.

#### 3.1.2. Isoler le champ `id`

Pour isoler le champ `id`, nous pouvons utiliser l'[opérateur `map`](https://www.nextflow.io/docs/latest/operator.html#map) pour créer un nouveau tuple avec le champ `id` comme premier élément.

=== "Après"

    ```groovy title="main.nf" linenums="7" hl_lines="3 6"
        ch_normal_samples = ch_samples
            .filter { meta, file -> meta.type == 'normal' }
            .map { meta, file -> [meta.id, meta, file] }
        ch_tumor_samples = ch_samples
            .filter { meta, file -> meta.type == 'tumor' }
            .map { meta, file -> [meta.id, meta, file] }
        ch_normal_samples
            .view{'Normal sample: ' + it}
        ch_tumor_samples
            .view{'Tumor sample: ' + it}
    ```

=== "Avant"

    ```groovy title="main.nf" linenums="7"
        ch_normal_samples = ch_samples
            .filter { meta, file -> meta.type == 'normal' }
        ch_tumor_samples = ch_samples
            .filter { meta, file -> meta.type == 'tumor' }
        ch_normal_samples
            .view{'Normal sample: ' + it}
        ch_tumor_samples
            .view{'Tumor sample: ' + it}
    ```

```bash
nextflow run main.nf
```

??? success "Sortie de la commande"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `main.nf` [mad_lagrange] DSL2 - revision: 9940b3f23d

    Tumor sample: [patientA, [id:patientA, repeat:1, type:tumor], patientA_rep1_tumor.bam]
    Tumor sample: [patientA, [id:patientA, repeat:2, type:tumor], patientA_rep2_tumor.bam]
    Normal sample: [patientA, [id:patientA, repeat:1, type:normal], patientA_rep1_normal.bam]
    Normal sample: [patientA, [id:patientA, repeat:2, type:normal], patientA_rep2_normal.bam]
    Tumor sample: [patientB, [id:patientB, repeat:1, type:tumor], patientB_rep1_tumor.bam]
    Tumor sample: [patientC, [id:patientC, repeat:1, type:tumor], patientC_rep1_tumor.bam]
    Normal sample: [patientB, [id:patientB, repeat:1, type:normal], patientB_rep1_normal.bam]
    Normal sample: [patientC, [id:patientC, repeat:1, type:normal], patientC_rep1_normal.bam]
    ```

C'est peut-être subtil, mais vous devriez pouvoir voir que le premier élément de chaque tuple est le champ `id`.

#### 3.1.3. Combiner les deux canaux

Maintenant, nous pouvons utiliser l'opérateur `join` pour combiner les deux canaux en fonction du champ `id`.

Une fois de plus, nous utiliserons `view` pour afficher les sorties jointes.

=== "Après"

    ```groovy title="main.nf" linenums="7" hl_lines="7-9"
        ch_normal_samples = ch_samples
            .filter { meta, file -> meta.type == 'normal' }
            .map { meta, file -> [meta.id, meta, file] }
        ch_tumor_samples = ch_samples
            .filter { meta, file -> meta.type == 'tumor' }
            .map { meta, file -> [meta.id, meta, file] }
        ch_joined_samples = ch_normal_samples
            .join(ch_tumor_samples)
        ch_joined_samples.view()
    ```

=== "Avant"

    ```groovy title="main.nf" linenums="7" hl_lines="7-10"
        ch_normal_samples = ch_samples
            .filter { meta, file -> meta.type == 'normal' }
            .map { meta, file -> [meta.id, meta, file] }
        ch_tumor_samples = ch_samples
            .filter { meta, file -> meta.type == 'tumor' }
            .map { meta, file -> [meta.id, meta, file] }
        ch_normal_samples
            .view{'Normal sample: ' + it}
        ch_tumor_samples
            .view{'Tumor sample: ' + it}
    ```

```bash
nextflow run main.nf
```

??? success "Sortie de la commande"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `main.nf` [soggy_wiles] DSL2 - revision: 3bc1979889

    [patientA, [id:patientA, repeat:1, type:normal], patientA_rep1_normal.bam, [id:patientA, repeat:1, type:tumor], patientA_rep1_tumor.bam]
    [patientA, [id:patientA, repeat:2, type:normal], patientA_rep2_normal.bam, [id:patientA, repeat:2, type:tumor], patientA_rep2_tumor.bam]
    [patientB, [id:patientB, repeat:1, type:normal], patientB_rep1_normal.bam, [id:patientB, repeat:1, type:tumor], patientB_rep1_tumor.bam]
    [patientC, [id:patientC, repeat:1, type:normal], patientC_rep1_normal.bam, [id:patientC, repeat:1, type:tumor], patientC_rep1_tumor.bam]
    ```

C'est un peu difficile à dire parce que c'est très large, mais vous devriez pouvoir voir que les échantillons ont été joints par le champ `id`. Chaque tuple a maintenant le format :

- `id` : L'identifiant de l'échantillon
- `normal_meta_map` : Les métadonnées de l'échantillon normal incluant le type, le réplicat et le chemin vers le fichier bam
- `normal_sample_file` : Le fichier de l'échantillon normal
- `tumor_meta_map` : Les métadonnées de l'échantillon tumoral incluant le type, le réplicat et le chemin vers le fichier bam
- `tumor_sample` : L'échantillon tumoral incluant le type, le réplicat et le chemin vers le fichier bam

!!! warning

    L'opérateur `join` éliminera tous les tuples non appariés. Dans cet exemple, nous nous sommes assurés que tous les échantillons étaient appariés pour tumeur et normal, mais si ce n'est pas le cas, vous devez utiliser le paramètre `remainder: true` pour conserver les tuples non appariés. Consultez la [documentation](https://www.nextflow.io/docs/latest/operator.html#join) pour plus de détails.

Maintenant vous savez comment utiliser `map` pour isoler un champ dans un tuple, et comment utiliser `join` pour combiner des tuples en fonction du premier champ.
Avec ces connaissances, nous pouvons combiner avec succès des canaux en fonction d'un champ partagé.

Ensuite, nous examinerons la situation où vous voulez joindre sur plusieurs champs.

### 3.2. Joindre sur plusieurs champs

Nous avons 2 réplicats pour l'échantillon A, mais seulement 1 pour les échantillons B et C. Dans ce cas, nous avons pu les joindre efficacement en utilisant le champ `id`, mais que se passerait-il s'ils n'étaient pas synchronisés ? Nous pourrions mélanger les échantillons normaux et tumoraux de différents réplicats !

Pour éviter cela, nous pouvons joindre sur plusieurs champs. Il existe en fait plusieurs façons d'y parvenir, mais nous allons nous concentrer sur la création d'une nouvelle clé de jointure qui inclut à la fois l'`id` de l'échantillon et le numéro de `replicate`.

Commençons par créer une nouvelle clé de jointure. Nous pouvons le faire de la même manière qu'avant en utilisant l'[opérateur `map`](https://www.nextflow.io/docs/latest/operator.html#map) pour créer un nouveau tuple avec les champs `id` et `repeat` comme premier élément.

=== "Après"

    ```groovy title="main.nf" linenums="7" hl_lines="3 6"
        ch_normal_samples = ch_samples
            .filter { meta, file -> meta.type == 'normal' }
            .map { meta, file -> [[meta.id, meta.repeat], meta, file] }
        ch_tumor_samples = ch_samples
            .filter { meta, file -> meta.type == 'tumor' }
            .map { meta, file -> [[meta.id, meta.repeat], meta, file] }
    ```

=== "Avant"

    ```groovy title="main.nf" linenums="7" hl_lines="3 6"
        ch_normal_samples = ch_samples
            .filter { meta, file -> meta.type == 'normal' }
            .map { meta, file -> [meta.id, meta, file] }
        ch_tumor_samples = ch_samples
            .filter { meta, file -> meta.type == 'tumor' }
            .map { meta, file -> [meta.id, meta, file] }
    ```

Maintenant, nous devrions voir que la jointure s'effectue en utilisant à la fois les champs `id` et `repeat`. Exécutez le workflow :

```bash
nextflow run main.nf
```

??? success "Sortie de la commande"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `main.nf` [prickly_wing] DSL2 - revision: 3bebf22dee

    [[patientA, 1], [id:patientA, repeat:1, type:normal], patientA_rep1_normal.bam, [id:patientA, repeat:1, type:tumor], patientA_rep1_tumor.bam]
    [[patientA, 2], [id:patientA, repeat:2, type:normal], patientA_rep2_normal.bam, [id:patientA, repeat:2, type:tumor], patientA_rep2_tumor.bam]
    [[patientB, 1], [id:patientB, repeat:1, type:normal], patientB_rep1_normal.bam, [id:patientB, repeat:1, type:tumor], patientB_rep1_tumor.bam]
    [[patientC, 1], [id:patientC, repeat:1, type:normal], patientC_rep1_normal.bam, [id:patientC, repeat:1, type:tumor], patientC_rep1_tumor.bam]
    ```

Notez comment nous avons un tuple de deux éléments (champs `id` et `repeat`) comme premier élément de chaque résultat joint. Cela démontre comment des éléments complexes peuvent être utilisés comme clé de jointure, permettant un appariement assez complexe entre les échantillons des mêmes conditions.

Si vous souhaitez explorer d'autres façons de joindre sur différentes clés, consultez la [documentation de l'opérateur join](https://www.nextflow.io/docs/latest/operator.html#join) pour des options et exemples supplémentaires.

### 3.3. Utiliser `subMap` pour créer une nouvelle clé de jointure

L'approche précédente perd les noms de champs de notre clé de jointure - les champs `id` et `repeat` deviennent simplement une liste de valeurs. Pour conserver les noms de champs pour un accès ultérieur, nous pouvons utiliser la [méthode `subMap`](<https://docs.groovy-lang.org/latest/html/groovy-jdk/java/util/Map.html#subMap(java.util.Collection)>).

La méthode `subMap` extrait uniquement les paires clé-valeur spécifiées d'une map. Ici, nous extrairons uniquement les champs `id` et `repeat` pour créer notre clé de jointure.

=== "Après"

    ```groovy title="main.nf" linenums="7" hl_lines="3 6"
        ch_normal_samples = ch_samples
            .filter { meta, file -> meta.type == 'normal' }
            .map { meta, file -> [meta.subMap(['id', 'repeat']), meta, file] }
        ch_tumor_samples = ch_samples
            .filter { meta, file -> meta.type == 'tumor' }
            .map { meta, file -> [meta.subMap(['id', 'repeat']), meta, file] }
    ```

=== "Avant"

    ```groovy title="main.nf" linenums="7" hl_lines="3 6"
        ch_normal_samples = ch_samples
            .filter { meta, file -> meta.type == 'normal' }
            .map { meta, file -> [[meta.id, meta.repeat], meta, file] }
        ch_tumor_samples = ch_samples
            .filter { meta, file -> meta.type == 'tumor' }
            .map { meta, file -> [[meta.id, meta.repeat], meta, file] }
    ```

```bash
nextflow run main.nf
```

??? success "Sortie de la commande"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `main.nf` [reverent_wing] DSL2 - revision: 847016c3b7

    [[id:patientA, repeat:1], [id:patientA, repeat:1, type:normal], patientA_rep1_normal.bam, [id:patientA, repeat:1, type:tumor], patientA_rep1_tumor.bam]
    [[id:patientA, repeat:2], [id:patientA, repeat:2, type:normal], patientA_rep2_normal.bam, [id:patientA, repeat:2, type:tumor], patientA_rep2_tumor.bam]
    [[id:patientB, repeat:1], [id:patientB, repeat:1, type:normal], patientB_rep1_normal.bam, [id:patientB, repeat:1, type:tumor], patientB_rep1_tumor.bam]
    [[id:patientC, repeat:1], [id:patientC, repeat:1, type:normal], patientC_rep1_normal.bam, [id:patientC, repeat:1, type:tumor], patientC_rep1_tumor.bam]
    ```

Maintenant, nous avons une nouvelle clé de jointure qui inclut non seulement les champs `id` et `repeat` mais conserve également les noms de champs afin que nous puissions y accéder ultérieurement par nom, par exemple `meta.id` et `meta.repeat`.

### 3.4. Utiliser une closure nommée dans map

Pour éviter la duplication et réduire les erreurs, nous pouvons utiliser une closure nommée. Une closure nommée nous permet de créer une fonction réutilisable que nous pouvons appeler à plusieurs endroits.

Pour ce faire, nous définissons d'abord la closure comme une nouvelle variable :

=== "Après"

    ```groovy title="main.nf" linenums="2" hl_lines="7"
        ch_samples = channel.fromPath("./data/samplesheet.csv")
            .splitCsv(header: true)
            .map{ row ->
                [[id:row.id, repeat:row.repeat, type:row.type], row.bam]
            }

        getSampleIdAndReplicate = { meta, bam -> [ meta.subMap(['id', 'repeat']), meta, file(bam) ] }

        ch_normal_samples = ch_samples
            .filter { meta, file -> meta.type == 'normal' }
    ```

=== "Avant"

    ```groovy title="main.nf" linenums="2"
        ch_samples = channel.fromPath("./data/samplesheet.csv")
            .splitCsv(header: true)
            .map{ row ->
                [[id:row.id, repeat:row.repeat, type:row.type], row.bam]
            }
        ch_normal_samples = ch_samples
            .filter { meta, file -> meta.type == 'normal' }
    ```

Nous avons défini la transformation map comme une variable nommée que nous pouvons réutiliser.

Notez que nous convertissons également le chemin du fichier en objet Path en utilisant `file()` afin que tout process recevant ce canal puisse gérer le fichier correctement (pour plus d'informations, voir [Working with files](./working_with_files.md)).

Implémentons la closure dans notre workflow :

=== "Après"

    ```groovy title="main.nf" linenums="10" hl_lines="3 6"
        ch_normal_samples = ch_samples
            .filter { meta, file -> meta.type == 'normal' }
             .map ( getSampleIdAndReplicate )
        ch_tumor_samples = ch_samples
            .filter { meta, file -> meta.type == 'tumor' }
             .map ( getSampleIdAndReplicate )

    ```

=== "Avant"

    ```groovy title="main.nf" linenums="10" hl_lines="3 6"
        ch_normal_samples = ch_samples
            .filter { meta, file -> meta.type == 'normal' }
            .map { meta, file -> [meta.subMap(['id', 'repeat']), meta, file] }
        ch_tumor_samples = ch_samples
            .filter { meta, file -> meta.type == 'tumor' }
            .map { meta, file -> [meta.subMap(['id', 'repeat']), meta, file] }
    ```

!!! note

    L'opérateur `map` est passé de l'utilisation de `{ }` à l'utilisation de `( )` pour passer la closure comme argument. C'est parce que l'opérateur `map` attend une closure comme argument et `{ }` est utilisé pour définir une closure anonyme. Lors de l'appel d'une closure nommée, utilisez la syntaxe `( )`.

Exécutez le workflow une fois de plus pour vérifier que tout fonctionne toujours :

```bash
nextflow run main.nf
```

??? success "Sortie de la commande"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `main.nf` [angry_meninsky] DSL2 - revision: 2edc226b1d

    [[id:patientA, repeat:1], [id:patientA, repeat:1, type:normal], patientA_rep1_normal.bam, [id:patientA, repeat:1, type:tumor], patientA_rep1_tumor.bam]
    [[id:patientA, repeat:2], [id:patientA, repeat:2, type:normal], patientA_rep2_normal.bam, [id:patientA, repeat:2, type:tumor], patientA_rep2_tumor.bam]
    [[id:patientB, repeat:1], [id:patientB, repeat:1, type:normal], patientB_rep1_normal.bam, [id:patientB, repeat:1, type:tumor], patientB_rep1_tumor.bam]
    [[id:patientC, repeat:1], [id:patientC, repeat:1, type:normal], patientC_rep1_normal.bam, [id:patientC, repeat:1, type:tumor], patientC_rep1_tumor.bam]
    ```

L'utilisation d'une closure nommée nous permet de réutiliser la même transformation à plusieurs endroits, réduisant le risque d'erreurs et rendant le code plus lisible et maintenable.

### 3.5. Réduire la duplication des données

Nous avons beaucoup de données dupliquées dans notre workflow. Chaque élément des échantillons joints répète les champs `id` et `repeat`. Puisque ces informations sont déjà disponibles dans la clé de regroupement, nous pouvons éviter cette redondance. Pour rappel, notre structure de données actuelle ressemble à ceci :

```groovy
[
  [
    "id": "sampleC",
    "repeat": "1",
  ],
  [
    "id": "sampleC",
    "repeat": "1",
    "type": "normal",
  ],
  "sampleC_rep1_normal.bam"
  [
    "id": "sampleC",
    "repeat": "1",
    "type": "tumor",
  ],
  "sampleC_rep1_tumor.bam"
]
```

Puisque les champs `id` et `repeat` sont disponibles dans la clé de regroupement, supprimons-les du reste de chaque élément de canal pour éviter la duplication. Nous pouvons le faire en utilisant la méthode `subMap` pour créer une nouvelle map avec uniquement le champ `type`. Cette approche nous permet de maintenir toutes les informations nécessaires tout en éliminant la redondance dans notre structure de données.

=== "Après"

    ```groovy title="main.nf" linenums="8" hl_lines="1"
        getSampleIdAndReplicate = { meta, bam -> [ meta.subMap(['id', 'repeat']), meta.subMap(['type']), file(bam) ] }
    ```

=== "Avant"

    ```groovy title="main.nf" linenums="8" hl_lines="1"
        getSampleIdAndReplicate = { meta, bam -> [ meta.subMap(['id', 'repeat']), meta, file(bam) ] }
    ```

Maintenant, la closure renvoie un tuple où le premier élément contient les champs `id` et `repeat`, et le deuxième élément contient uniquement le champ `type`. Cela élimine la redondance en stockant les informations `id` et `repeat` une seule fois dans la clé de regroupement, tout en maintenant toutes les informations nécessaires.

Exécutez le workflow pour voir à quoi cela ressemble :

```bash
nextflow run main.nf
```

??? success "Sortie de la commande"

    ```console
    [[id:patientA, repeat:1], [type:normal], /workspaces/training/side-quests/splitting_and_grouping/patientA_rep1_normal.bam, [type:tumor], /workspaces/training/side-quests/splitting_and_grouping/patientA_rep1_tumor.bam]
    [[id:patientA, repeat:2], [type:normal], /workspaces/training/side-quests/splitting_and_grouping/patientA_rep2_normal.bam, [type:tumor], /workspaces/training/side-quests/splitting_and_grouping/patientA_rep2_tumor.bam]
    [[id:patientB, repeat:1], [type:normal], /workspaces/training/side-quests/splitting_and_grouping/patientB_rep1_normal.bam, [type:tumor], /workspaces/training/side-quests/splitting_and_grouping/patientB_rep1_tumor.bam]
    [[id:patientC, repeat:1], [type:normal], /workspaces/training/side-quests/splitting_and_grouping/patientC_rep1_normal.bam, [type:tumor], /workspaces/training/side-quests/splitting_and_grouping/patientC_rep1_tumor.bam]
    ```

Nous pouvons voir que nous ne déclarons les champs `id` et `repeat` qu'une seule fois dans la clé de regroupement et que nous avons le champ `type` dans les données d'échantillon. Nous n'avons perdu aucune information mais nous avons réussi à rendre le contenu de notre canal plus succinct.

### 3.6. Supprimer les informations redondantes

Nous avons supprimé les informations dupliquées ci-dessus, mais nous avons encore d'autres informations redondantes dans nos canaux.

Au début, nous avons séparé les échantillons normaux et tumoraux en utilisant `filter`, puis nous les avons joints en fonction des clés `id` et `repeat`. L'opérateur `join` préserve l'ordre dans lequel les tuples sont fusionnés, donc dans notre cas, avec les échantillons normaux du côté gauche et les échantillons tumoraux du côté droit, le canal résultant maintient cette structure : `id, <éléments normaux>, <éléments tumoraux>`.

Puisque nous connaissons la position de chaque élément dans notre canal, nous pouvons simplifier davantage la structure en supprimant les métadonnées `[type:normal]` et `[type:tumor]`.

=== "Après"

    ```groovy title="main.nf" linenums="8" hl_lines="1"
        getSampleIdAndReplicate = { meta, file -> [ meta.subMap(['id', 'repeat']), file ] }
    ```

=== "Avant"

    ```groovy title="main.nf" linenums="8" hl_lines="1"
        getSampleIdAndReplicate = { meta, file -> [ meta.subMap(['id', 'repeat']), meta.subMap(['type']), file ] }
    ```

Exécutez à nouveau pour voir le résultat :

```bash
nextflow run main.nf
```

??? success "Sortie de la commande"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `main.nf` [confident_leavitt] DSL2 - revision: a2303895bd

    [[id:patientA, repeat:1], patientA_rep1_normal.bam, patientA_rep1_tumor.bam]
    [[id:patientA, repeat:2], patientA_rep2_normal.bam, patientA_rep2_tumor.bam]
    [[id:patientB, repeat:1], patientB_rep1_normal.bam, patientB_rep1_tumor.bam]
    [[id:patientC, repeat:1], patientC_rep1_normal.bam, patientC_rep1_tumor.bam]
    ```

### Point clé

Dans cette section, vous avez appris :

- **Manipuler les tuples** : Comment utiliser `map` pour isoler un champ dans un tuple
- **Joindre les tuples** : Comment utiliser `join` pour combiner des tuples en fonction du premier champ
- **Créer des clés de jointure** : Comment utiliser `subMap` pour créer une nouvelle clé de jointure
- **Closures nommées** : Comment utiliser une closure nommée dans map
- **Jointure sur plusieurs champs** : Comment joindre sur plusieurs champs pour un appariement plus précis
- **Optimisation de la structure de données** : Comment rationaliser la structure du canal en supprimant les informations redondantes

Vous disposez maintenant d'un workflow capable de diviser une feuille d'échantillons, filtrer les échantillons normaux et tumoraux, les joindre par identifiant d'échantillon et numéro de réplicat, puis afficher les résultats.

Il s'agit d'un modèle courant dans les workflows de bioinformatique où vous devez apparier des échantillons ou d'autres types de données après un traitement indépendant, c'est donc une compétence utile. Ensuite, nous examinerons la répétition d'un échantillon plusieurs fois.

## 4. Distribuer les échantillons sur des intervalles

Un modèle clé dans les workflows de bioinformatique est la distribution de l'analyse sur des régions génomiques. Par exemple, l'appel de variants peut être parallélisé en divisant le génome en intervalles (comme des chromosomes ou des régions plus petites). Cette stratégie de parallélisation améliore considérablement l'efficacité du pipeline en distribuant la charge de calcul sur plusieurs cœurs ou nœuds, réduisant ainsi le temps d'exécution global.

Dans la section suivante, nous démontrerons comment distribuer nos données d'échantillons sur plusieurs intervalles génomiques. Nous associerons chaque échantillon à chaque intervalle, permettant le traitement parallèle de différentes régions génomiques. Cela multipliera la taille de notre ensemble de données par le nombre d'intervalles, créant plusieurs unités d'analyse indépendantes qui pourront être rassemblées ultérieurement.

### 4.1. Distribuer les échantillons sur des intervalles en utilisant `combine`

Commençons par créer un canal d'intervalles. Pour simplifier les choses, nous utiliserons simplement 3 intervalles que nous définirons manuellement. Dans un vrai workflow, vous pourriez les lire à partir d'un fichier d'entrée ou même créer un canal avec de nombreux fichiers d'intervalles.

=== "Après"

    ```groovy title="main.nf" linenums="17" hl_lines="2"
            .join(ch_tumor_samples)
        ch_intervals = channel.of('chr1', 'chr2', 'chr3')
    ```

=== "Avant"

    ```groovy title="main.nf" linenums="17" hl_lines="2"
            .join(ch_tumor_samples)
        ch_joined_samples.view()
    ```

Maintenant, souvenez-vous, nous voulons répéter chaque échantillon pour chaque intervalle. Ceci est parfois appelé le produit cartésien des échantillons et des intervalles. Nous pouvons y parvenir en utilisant l'[opérateur `combine`](https://www.nextflow.io/docs/latest/operator.html#combine). Cela prendra chaque élément du canal 1 et le répétera pour chaque élément du canal 2. Ajoutons un opérateur combine à notre workflow :

=== "Après"

    ```groovy title="main.nf" linenums="18" hl_lines="3-5"
        ch_intervals = channel.of('chr1', 'chr2', 'chr3')

        ch_combined_samples = ch_joined_samples
            .combine(ch_intervals)
            .view()
    ```

=== "Avant"

    ```groovy title="main.nf" linenums="18"
        ch_intervals = channel.of('chr1', 'chr2', 'chr3')
    ```

Maintenant, exécutons-le et voyons ce qui se passe :

```bash
nextflow run main.nf
```

??? success "Sortie de la commande"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `main.nf` [mighty_tesla] DSL2 - revision: ae013ab70b

    [[id:patientA, repeat:1], patientA_rep1_normal.bam, patientA_rep1_tumor.bam, chr1]
    [[id:patientA, repeat:1], patientA_rep1_normal.bam, patientA_rep1_tumor.bam, chr2]
    [[id:patientA, repeat:1], patientA_rep1_normal.bam, patientA_rep1_tumor.bam, chr3]
    [[id:patientA, repeat:2], patientA_rep2_normal.bam, patientA_rep2_tumor.bam, chr1]
    [[id:patientA, repeat:2], patientA_rep2_normal.bam, patientA_rep2_tumor.bam, chr2]
    [[id:patientA, repeat:2], patientA_rep2_normal.bam, patientA_rep2_tumor.bam, chr3]
    [[id:patientB, repeat:1], patientB_rep1_normal.bam, patientB_rep1_tumor.bam, chr1]
    [[id:patientB, repeat:1], patientB_rep1_normal.bam, patientB_rep1_tumor.bam, chr2]
    [[id:patientB, repeat:1], patientB_rep1_normal.bam, patientB_rep1_tumor.bam, chr3]
    [[id:patientC, repeat:1], patientC_rep1_normal.bam, patientC_rep1_tumor.bam, chr1]
    [[id:patientC, repeat:1], patientC_rep1_normal.bam, patientC_rep1_tumor.bam, chr2]
    [[id:patientC, repeat:1], patientC_rep1_normal.bam, patientC_rep1_tumor.bam, chr3]
    ```

Succès ! Nous avons répété chaque échantillon pour chaque intervalle de notre liste de 3 intervalles. Nous avons effectivement triplé le nombre d'éléments dans notre canal.

C'est un peu difficile à lire cependant, donc dans la section suivante, nous allons le rendre plus clair.

### 4.2. Organiser le canal

Nous pouvons utiliser l'opérateur `map` pour ranger et refactoriser nos données d'échantillons afin qu'elles soient plus faciles à comprendre. Déplaçons la chaîne d'intervalles vers la map de jointure au premier élément.

=== "Après"

    ```groovy title="main.nf" linenums="20" hl_lines="3-9"
        ch_combined_samples = ch_joined_samples
            .combine(ch_intervals)
            .map { grouping_key, normal, tumor, interval ->
                [
                    grouping_key + [interval: interval],
                    normal,
                    tumor
                ]
            }
            .view()
    ```

=== "Avant"

    ```groovy title="main.nf" linenums="20"
        ch_combined_samples = ch_joined_samples
            .combine(ch_intervals)
            .view()
    ```

Décomposons ce que fait cette opération map étape par étape.

Tout d'abord, nous utilisons des paramètres nommés pour rendre le code plus lisible. En utilisant les noms `grouping_key`, `normal`, `tumor` et `interval`, nous pouvons faire référence aux éléments du tuple par nom au lieu de par index :

```groovy
        .map { grouping_key, normal, tumor, interval ->
```

Ensuite, nous combinons la `grouping_key` avec le champ `interval`. La `grouping_key` est une map contenant les champs `id` et `repeat`. Nous créons une nouvelle map avec l'`interval` et les fusionnons en utilisant l'addition de maps de Groovy (`+`) :

```groovy
                grouping_key + [interval: interval],
```

Enfin, nous retournons cela comme un tuple avec trois éléments : la map de métadonnées combinée, le fichier d'échantillon normal et le fichier d'échantillon tumoral :

```groovy
            [
                grouping_key + [interval: interval],
                normal,
                tumor
            ]
```

Exécutons-le à nouveau et vérifions le contenu du canal :

```bash
nextflow run main.nf
```

??? success "Sortie de la commande"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `main.nf` [sad_hawking] DSL2 - revision: 1f6f6250cd

    [[id:patientA, repeat:1, interval:chr1], patientA_rep1_normal.bam, patientA_rep1_tumor.bam]
    [[id:patientA, repeat:1, interval:chr2], patientA_rep1_normal.bam, patientA_rep1_tumor.bam]
    [[id:patientA, repeat:1, interval:chr3], patientA_rep1_normal.bam, patientA_rep1_tumor.bam]
    [[id:patientA, repeat:2, interval:chr1], patientA_rep2_normal.bam, patientA_rep2_tumor.bam]
    [[id:patientA, repeat:2, interval:chr2], patientA_rep2_normal.bam, patientA_rep2_tumor.bam]
    [[id:patientA, repeat:2, interval:chr3], patientA_rep2_normal.bam, patientA_rep2_tumor.bam]
    [[id:patientB, repeat:1, interval:chr1], patientB_rep1_normal.bam, patientB_rep1_tumor.bam]
    [[id:patientB, repeat:1, interval:chr2], patientB_rep1_normal.bam, patientB_rep1_tumor.bam]
    [[id:patientB, repeat:1, interval:chr3], patientB_rep1_normal.bam, patientB_rep1_tumor.bam]
    [[id:patientC, repeat:1, interval:chr1], patientC_rep1_normal.bam, patientC_rep1_tumor.bam]
    [[id:patientC, repeat:1, interval:chr2], patientC_rep1_normal.bam, patientC_rep1_tumor.bam]
    [[id:patientC, repeat:1, interval:chr3], patientC_rep1_normal.bam, patientC_rep1_tumor.bam]
    ```

Utiliser `map` pour contraindre vos données dans la structure correcte peut être délicat, mais c'est crucial pour une manipulation efficace des données.

Nous avons maintenant chaque échantillon répété sur tous les intervalles génomiques, créant plusieurs unités d'analyse indépendantes qui peuvent être traitées en parallèle. Mais que se passe-t-il si nous voulons rassembler des échantillons connexes ? Dans la section suivante, nous apprendrons comment regrouper des échantillons qui partagent des attributs communs.

### Point clé

Dans cette section, vous avez appris :

- **Distribuer les échantillons sur des intervalles** : Comment utiliser `combine` pour répéter les échantillons sur des intervalles
- **Créer des produits cartésiens** : Comment générer toutes les combinaisons d'échantillons et d'intervalles
- **Organiser la structure du canal** : Comment utiliser `map` pour restructurer les données pour une meilleure lisibilité
- **Préparation au traitement parallèle** : Comment préparer les données pour une analyse distribuée

## 5. Agréger des échantillons en utilisant `groupTuple`

Dans les sections précédentes, nous avons appris comment diviser les données d'un fichier d'entrée et filtrer par des champs spécifiques (dans notre cas, les échantillons normaux et tumoraux). Mais cela ne couvre qu'un seul type de jointure. Et si nous voulons regrouper les échantillons par un attribut spécifique ? Par exemple, au lieu de joindre des paires appariées normal-tumeur, nous pourrions vouloir traiter tous les échantillons de "sampleA" ensemble indépendamment de leur type. Ce modèle est courant dans les workflows de bioin
