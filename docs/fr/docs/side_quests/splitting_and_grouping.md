# Séparation et Regroupement

Nextflow fournit des outils puissants pour travailler avec les données de manière flexible. Une capacité clé est la séparation des données en différents flux, puis le regroupement d'éléments liés. Cela est particulièrement utile dans les workflows de bioinformatique où vous devez traiter différents types d'échantillons séparément avant de combiner les résultats pour l'analyse.

Pensez-y comme au tri du courrier : vous séparez les lettres par destination, traitez chaque pile différemment, puis recombinez les éléments allant à la même personne. Nextflow utilise des opérateurs spéciaux pour accomplir cela avec des données scientifiques. Cette approche est également communément connue sous le nom de modèle **scatter/gather** dans le calcul distribué et les workflows de bioinformatique.

Le système de canaux de Nextflow est au cœur de cette flexibilité. Les canaux connectent différentes parties de votre workflow, permettant aux données de circuler à travers votre analyse. Vous pouvez créer plusieurs canaux à partir d'une seule source de données, traiter chaque canal différemment, puis fusionner les canaux ensemble lorsque nécessaire. Cette approche vous permet de concevoir des workflows qui reflètent naturellement les chemins de ramification et de convergence d'analyses bioinformatiques complexes.

### Objectifs d'apprentissage

Dans cette quête secondaire, vous apprendrez à séparer et regrouper des données en utilisant les opérateurs de canaux de Nextflow.
Nous commencerons avec un fichier CSV contenant des informations sur les échantillons et les fichiers de données associés, puis manipulerons et réorganiserons ces données.

À la fin de cette quête secondaire, vous serez capable de séparer et combiner efficacement des flux de données, en utilisant les techniques suivantes :

- Lire des données à partir de fichiers en utilisant `splitCsv`
- Filtrer et transformer des données avec `filter` et `map`
- Combiner des données liées en utilisant `join` et `groupTuple`
- Créer des combinaisons de données avec `combine` pour le traitement parallèle
- Optimiser la structure des données en utilisant `subMap` et des stratégies de déduplication
- Construire des fonctions réutilisables avec des closures nommées pour vous aider à manipuler les structures de canaux

Ces compétences vous aideront à construire des workflows capables de gérer plusieurs fichiers d'entrée et différents types de données efficacement, tout en maintenant une structure de code propre et maintenable.

### Prérequis

Avant d'entreprendre cette quête secondaire, vous devriez :

- Avoir complété le tutoriel [Hello Nextflow](../hello_nextflow/README.md) ou un cours équivalent pour débutant·es.
- Être à l'aise avec les concepts et mécanismes de base de Nextflow (processus, canaux, opérateurs, travail avec des fichiers, métadonnées)

**Optionnel :** Nous recommandons de compléter d'abord la quête secondaire [Métadonnées dans les workflows](./metadata.md).
Celle-ci couvre les fondamentaux de la lecture de fichiers CSV avec `splitCsv` et de la création de maps de métadonnées, que nous utiliserons intensivement ici.

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

#### Examiner le matériel

Vous trouverez un fichier de workflow principal et un répertoire `data` contenant une feuille d'échantillons nommée `samplesheet.csv`.

```console title="Directory contents"
.
├── data
│   └── samplesheet.csv
└── main.nf
```

La feuille d'échantillons contient des informations sur des échantillons de différent·es patient·es, incluant l'ID du/de la patient·e, le numéro de répétition de l'échantillon, le type (normal ou tumoral), et les chemins vers des fichiers de données hypothétiques (qui n'existent pas réellement, mais nous ferons comme s'ils existaient).

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

Cette feuille d'échantillons liste huit échantillons de trois patient·es (A, B, C).

Pour chaque patient·e, nous avons des échantillons de type `tumor` (provenant généralement de biopsies tumorales) ou `normal` (prélevés sur des tissus sains ou du sang).
Si vous n'êtes pas familier·ère avec l'analyse du cancer, sachez simplement que cela correspond à un modèle expérimental qui utilise des échantillons appariés tumeur/normal pour effectuer des analyses contrastives.

Pour le/la patient·e A spécifiquement, nous avons deux ensembles de réplicats techniques (répétitions).

!!! note

    Ne vous inquiétez pas si vous n'êtes pas familier·ère avec ce plan expérimental, ce n'est pas critique pour comprendre ce tutoriel.

#### Examiner l'assignation

Votre défi est d'écrire un workflow Nextflow qui va :

1. **Lire** les données d'échantillons à partir d'un fichier CSV et les structurer avec des maps de métadonnées
2. **Séparer** les échantillons en différents canaux selon le type (normal vs tumoral)
3. **Joindre** les paires tumeur/normal appariées par ID de patient·e et numéro de réplicat
4. **Distribuer** les échantillons à travers des intervalles génomiques pour le traitement parallèle
5. **Regrouper** les échantillons liés ensemble pour l'analyse en aval

Cela représente un modèle bioinformatique courant où vous devez séparer les données pour un traitement indépendant, puis recombiner les éléments liés pour une analyse comparative.

#### Liste de vérification de préparation

Pensez-vous être prêt·e à vous lancer ?

- [ ] Je comprends l'objectif de ce cours et ses prérequis
- [ ] Mon codespace est opérationnel
- [ ] J'ai défini mon répertoire de travail de manière appropriée
- [ ] Je comprends l'assignation

Si vous pouvez cocher toutes les cases, vous êtes prêt·e à commencer.

---

## 1. Lire les données d'échantillons

### 1.1. Lire les données d'échantillons avec `splitCsv` et créer des maps de métadonnées

Commençons par lire les données d'échantillons avec `splitCsv` et les organiser selon le modèle de map de métadonnées. Dans le `main.nf`, vous verrez que nous avons déjà commencé le workflow.

```groovy title="main.nf" linenums="1" hl_lines="2"
workflow {
    ch_samplesheet = channel.fromPath("./data/samplesheet.csv")
}
```

!!! note

    Tout au long de ce tutoriel, nous utiliserons le préfixe `ch_` pour toutes les variables de canaux afin d'indiquer clairement qu'il s'agit de canaux Nextflow.

Si vous avez complété la quête secondaire [Métadonnées dans les workflows](./metadata.md), vous reconnaîtrez ce modèle. Nous utiliserons `splitCsv` pour lire le CSV et structurer immédiatement les données avec une map de métadonnées pour séparer les métadonnées des chemins de fichiers.

!!! info

    Nous rencontrerons deux concepts différents appelés `map` dans cette formation :

    - **Structure de données** : La map Groovy (équivalent aux dictionnaires/hashes dans d'autres langages) qui stocke des paires clé-valeur
    - **Opérateur de canal** : L'opérateur `.map()` qui transforme les éléments dans un canal

    Nous clarifierons de laquelle nous parlons selon le contexte, mais cette distinction est importante à comprendre lorsque vous travaillez avec Nextflow.

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

Cela combine l'opération `splitCsv` (lecture du CSV avec en-têtes) et l'opération `map` (structuration des données sous forme de tuples `[meta, file]`) en une seule étape. Appliquez ce changement et exécutez le pipeline :

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

Nous avons maintenant un canal où chaque élément est un tuple `[meta, file]` - métadonnées séparées des chemins de fichiers. Cette structure nous permet de séparer et regrouper notre charge de travail en fonction des champs de métadonnées.

---

## 2. Filtrer et transformer les données

### 2.1. Filtrer les données avec `filter`

Nous pouvons utiliser l'[opérateur `filter`](https://www.nextflow.io/docs/latest/operator.html#filter) pour filtrer les données selon une condition. Disons que nous voulons traiter uniquement les échantillons normaux. Nous pouvons le faire en filtrant les données selon le champ `type`. Insérons ceci avant l'opérateur `view`.

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

Nous avons réussi à filtrer les données pour n'inclure que les échantillons normaux. Récapitulons comment cela fonctionne.

L'opérateur `filter` prend une closure qui est appliquée à chaque élément du canal. Si la closure retourne `true`, l'élément est inclus ; si elle retourne `false`, l'élément est exclu.

Dans notre cas, nous voulons garder uniquement les échantillons où `meta.type == 'normal'`. La closure utilise le tuple `meta,file` pour référer à chaque échantillon, accède au type d'échantillon avec `meta.type`, et vérifie s'il est égal à `'normal'`.

Ceci est accompli avec la closure unique que nous avons introduite ci-dessus :

```groovy title="main.nf" linenums="7"
    .filter { meta, file -> meta.type == 'normal' }
```

### 2.2. Créer des canaux filtrés séparés

Actuellement, nous appliquons le filtre au canal créé directement à partir du CSV, mais nous voulons filtrer cela de plusieurs façons, donc réécrivons la logique pour créer un canal filtré séparé pour les échantillons normaux :

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

Nous avons séparé les échantillons normaux et tumoraux en deux canaux différents, et utilisé une closure fournie à `view()` pour les étiqueter différemment dans la sortie : `ch_tumor_samples.view{'Tumor sample: ' + it}`.

### À retenir

Dans cette section, vous avez appris :

- **Filtrer les données** : Comment filtrer les données avec `filter`
- **Séparer les données** : Comment séparer les données en différents canaux selon une condition
- **Visualiser les données** : Comment utiliser `view` pour afficher les données et étiqueter la sortie de différents canaux

Nous avons maintenant séparé les échantillons normaux et tumoraux en deux canaux différents. Ensuite, nous joindrons les échantillons normaux et tumoraux sur le champ `id`.

---

## 3. Joindre les canaux par identifiants

Dans la section précédente, nous avons séparé les échantillons normaux et tumoraux en deux canaux différents. Ceux-ci pourraient être traités indépendamment en utilisant des processus ou workflows spécifiques selon leur type. Mais que se passe-t-il lorsque nous voulons comparer les échantillons normaux et tumoraux du/de la même patient·e ? À ce stade, nous devons les joindre à nouveau en nous assurant de faire correspondre les échantillons selon leur champ `id`.

Nextflow inclut de nombreuses méthodes pour combiner des canaux, mais dans ce cas l'opérateur le plus approprié est [`join`](https://www.nextflow.io/docs/latest/operator.html#join). Si vous êtes familier·ère avec SQL, il agit comme l'opération `JOIN`, où nous spécifions la clé sur laquelle joindre et le type de jointure à effectuer.

### 3.1. Utiliser `map` et `join` pour combiner selon l'ID de patient·e

Si nous vérifions la documentation de [`join`](https://www.nextflow.io/docs/latest/operator.html#join), nous pouvons voir que par défaut il joint deux canaux selon le premier élément de chaque tuple.

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

Nous pouvons voir que le champ `id` est le premier élément de chaque map de métadonnées. Pour que `join` fonctionne, nous devrions isoler le champ `id` dans chaque tuple. Après cela, nous pouvons simplement utiliser l'opérateur `join` pour combiner les deux canaux.

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

Maintenant nous pouvons utiliser l'opérateur `join` pour combiner les deux canaux selon le champ `id`.

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

C'est un peu difficile à dire car c'est très large, mais vous devriez pouvoir voir que les échantillons ont été joints par le champ `id`. Chaque tuple a maintenant le format :

- `id` : L'ID de l'échantillon
- `normal_meta_map` : Les métadonnées de l'échantillon normal incluant le type, le réplicat et le chemin vers le fichier bam
- `normal_sample_file` : Le fichier de l'échantillon normal
- `tumor_meta_map` : Les métadonnées de l'échantillon tumoral incluant le type, le réplicat et le chemin vers le fichier bam
- `tumor_sample` : L'échantillon tumoral incluant le type, le réplicat et le chemin vers le fichier bam

!!! warning

    L'opérateur `join` éliminera tous les tuples non appariés. Dans cet exemple, nous nous sommes assuré·es que tous les échantillons étaient appariés pour tumeur et normal mais si ce n'est pas le cas vous devez utiliser le paramètre `remainder: true` pour conserver les tuples non appariés. Consultez la [documentation](https://www.nextflow.io/docs/latest/operator.html#join) pour plus de détails.

Vous savez maintenant comment utiliser `map` pour isoler un champ dans un tuple, et comment utiliser `join` pour combiner des tuples selon le premier champ.
Avec cette connaissance, nous pouvons combiner avec succès des canaux selon un champ partagé.

Ensuite, nous considérerons la situation où vous voulez joindre sur plusieurs champs.

### 3.2. Joindre sur plusieurs champs

Nous avons 2 réplicats pour l'échantillonA, mais seulement 1 pour l'échantillonB et l'échantillonC. Dans ce cas, nous avons pu les joindre efficacement en utilisant le champ `id`, mais que se passerait-il s'ils n'étaient pas synchronisés ? Nous pourrions mélanger les échantillons normaux et tumoraux de différents réplicats !

Pour éviter cela, nous pouvons joindre sur plusieurs champs. Il existe en fait plusieurs façons d'y parvenir mais nous allons nous concentrer sur la création d'une nouvelle clé de jointure qui inclut à la fois l'`id` de l'échantillon et le numéro de `replicate`.

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

Maintenant nous devrions voir la jointure se produire mais en utilisant à la fois les champs `id` et `repeat`. Exécutez le workflow :

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

Notez comment nous avons un tuple de deux éléments (champs `id` et `repeat`) comme premier élément de chaque résultat joint. Cela démontre comment des éléments complexes peuvent être utilisés comme clé de jointure, permettant des correspondances assez complexes entre échantillons des mêmes conditions.

Si vous voulez explorer plus de façons de joindre sur différentes clés, consultez la [documentation de l'opérateur join](https://www.nextflow.io/docs/latest/operator.html#join) pour des options et exemples supplémentaires.

### 3.3. Utiliser `subMap` pour créer une nouvelle clé de jointure

L'approche précédente perd les noms de champs de notre clé de jointure - les champs `id` et `repeat` deviennent juste une liste de valeurs. Pour conserver les noms de champs pour un accès ultérieur, nous pouvons utiliser la [méthode `subMap`](<https://docs.groovy-lang.org/latest/html/groovy-jdk/java/util/Map.html#subMap(java.util.Collection)>).

La méthode `subMap` extrait uniquement les paires clé-valeur spécifiées d'une map. Ici nous extrairons juste les champs `id` et `repeat` pour créer notre clé de jointure.

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

Maintenant nous avons une nouvelle clé de jointure qui non seulement inclut les champs `id` et `repeat` mais conserve également les noms de champs afin que nous puissions y accéder plus tard par nom, par exemple `meta.id` et `meta.repeat`.

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

Notez que nous convertissons également le chemin de fichier en objet Path en utilisant `file()` afin que tout processus recevant ce canal puisse gérer le fichier correctement (pour plus d'informations voir [Travailler avec des fichiers](./working_with_files.md)).

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

### 3.5. Réduire la duplication de données

Nous avons beaucoup de données dupliquées dans notre workflow. Chaque élément des échantillons joints répète les champs `id` et `repeat`. Puisque cette information est déjà disponible dans la clé de regroupement, nous pouvons éviter cette redondance. Pour rappel, notre structure de données actuelle ressemble à ceci :

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

Maintenant la closure retourne un tuple où le premier élément contient les champs `id` et `repeat`, et le deuxième élément contient uniquement le champ `type`. Cela élimine la redondance en stockant les informations `id` et `repeat` une fois dans la clé de regroupement, tout en maintenant toutes les informations nécessaires.

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

Nous pouvons voir que nous ne déclarons les champs `id` et `repeat` qu'une seule fois dans la clé de regroupement et nous avons le champ `type` dans les données d'échantillon. Nous n'avons perdu aucune information mais nous avons réussi à rendre le contenu de notre canal plus succinct.

### 3.6. Supprimer les informations redondantes

Nous avons supprimé les informations dupliquées ci-dessus, mais nous avons encore d'autres informations redondantes dans nos canaux.

Au début, nous avons séparé les échantillons normaux et tumoraux en utilisant `filter`, puis les avons joints selon les clés `id` et `repeat`. L'opérateur `join` préserve l'ordre dans lequel les tuples sont fusionnés, donc dans notre cas, avec les échantillons normaux du côté gauche et les échantillons tumoraux du côté droit, le canal résultant maintient cette structure : `id, <éléments normaux>, <éléments tumoraux>`.

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

### À retenir

Dans cette section, vous avez appris :

- **Manipuler les tuples** : Comment utiliser `map` pour isoler un champ dans un tuple
- **Joindre les tuples** : Comment utiliser `join` pour combiner des tuples selon le premier champ
- **Créer des clés de jointure** : Comment utiliser `subMap` pour créer une nouvelle clé de jointure
- **Closures nommées** : Comment utiliser une closure nommée dans map
- **Jointure sur plusieurs champs** : Comment joindre sur plusieurs champs pour une correspondance plus précise
- **Optimisation de la structure de données** : Comment rationaliser la structure de canal en supprimant les informations redondantes

Vous avez maintenant un workflow qui peut séparer une feuille d'échantillons, filtrer les échantillons normaux et tumoraux, les joindre ensemble par ID d'échantillon et numéro de réplicat, puis afficher les résultats.

C'est un modèle courant dans les workflows de bioinformatique où vous devez faire correspondre des échantillons ou d'autres types de données après traitement indépendant, c'est donc une compétence utile. Ensuite, nous examinerons la répétition d'un échantillon plusieurs fois.

## 4. Répartir les échantillons sur des intervalles

Un modèle clé dans les workflows de bioinformatique est la distribution de l'analyse à travers des régions génomiques. Par exemple, l'appel de variants peut être parallélisé en divisant le génome en intervalles (comme des chromosomes ou des régions plus petites). Cette stratégie de parallélisation améliore significativement l'efficacité du pipeline en distribuant la charge de calcul sur plusieurs cœurs ou nœuds, réduisant le temps d'exécution global.

Dans la section suivante, nous démontrerons comment distribuer nos données d'échantillons à travers plusieurs intervalles génomiques. Nous apparierons chaque échantillon avec chaque intervalle, permettant le traitement parallèle de différentes régions génomiques. Cela multipliera la taille de notre ensemble de données par le nombre d'intervalles, créant plusieurs unités d'analyse indépendantes qui peuvent être rassemblées plus tard.

### 4.1. Répartir les échantillons sur des intervalles en utilisant `combine`

Commençons par créer un canal d'intervalles. Pour simplifier, nous utiliserons simplement 3 intervalles que nous définirons manuellement. Dans un workflow réel, vous pourriez les lire à partir d'une entrée de fichier ou même créer un canal avec beaucoup de fichiers d'intervalles.

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

Maintenant rappelez-vous, nous voulons répéter chaque échantillon pour chaque intervalle. C'est parfois appelé le produit cartésien des échantillons et intervalles. Nous pouvons y parvenir en utilisant l'[opérateur `combine`](https://www.nextflow.io/docs/latest/operator.html#combine). Cela prendra chaque élément du canal 1 et le répétera pour chaque élément du canal 2. Ajoutons un opérateur combine à notre workflow :

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

Maintenant exécutons-le et voyons ce qui se passe :

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

Succès ! Nous avons répété chaque échantillon pour chaque intervalle unique dans notre liste de 3 intervalles. Nous avons effectivement triplé le nombre d'éléments dans notre canal.

C'est un peu difficile à lire cependant, donc dans la section suivante nous allons le ranger.

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

D'abord, nous utilisons des paramètres nommés pour rendre le code plus lisible. En utilisant les noms `grouping_key`, `normal`, `tumor` et `interval`, nous pouvons référer aux éléments du tuple par nom au lieu de par index :

```groovy
        .map { grouping_key, normal, tumor, interval ->
```

Ensuite, nous combinons la `grouping_key` avec le champ `interval`. La `grouping_key` est une map contenant les champs `id` et `repeat`. Nous créons une nouvelle map avec l'`interval` et les fusionnons en utilisant l'addition de maps de Groovy (`+`) :

```groovy
                grouping_key + [interval: interval],
```

Enfin, nous retournons ceci comme un tuple avec trois éléments : la map de métadonnées combinée, le fichier d'échantillon normal, et le fichier d'échantillon tumoral :

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

Nous avons maintenant chaque échantillon répété à travers tous les intervalles génomiques, créant plusieurs unités d'analyse indépendantes qui peuvent être traitées en parallèle. Mais que faire si nous voulons rassembler des échantillons liés ? Dans la section suivante, nous apprendrons comment regrouper des échantillons qui partagent des attributs communs.

### À retenir

Dans cette section, vous avez appris :

- **Répartir les échantillons sur des intervalles** : Comment utiliser `combine` pour répéter des échantillons sur des intervalles
- **Créer des produits cartésiens** : Comment générer toutes les combinaisons d'échantillons et d'intervalles
- **Organiser la structure de canal** : Comment utiliser `map` pour restructurer les données pour une meilleure lisibilité
- **Préparation au traitement parallèle** : Comment configurer les données pour une analyse distribuée

## 5. Agréger les échantillons en utilisant `groupTuple`

Dans les sections précédentes, nous avons appris comment séparer les données d'un fichier d'entrée et filtrer par champs spécifiques (dans notre cas les échantillons normaux et tumoraux). Mais cela ne couvre qu'un seul type de jointure. Que faire si nous voulons regrouper des échantillons par un attribut spécifique ? Par exemple, au lieu de joindre des paires tumeur-normal appariées, nous pourrions vouloir traiter tous les échantillons de "sampleA" ensemble indépendamment de leur type. Ce modèle est courant dans les workflows de bioinformatique où vous pouvez vouloir traiter des échantillons liés séparément pour des raisons d'efficacité avant de comparer ou combiner les résultats à la fin.

Nextflow inclut des méthodes intégrées pour faire cela, la principale que nous examinerons est `groupTuple`.

Commençons par regrouper tous nos échantillons qui ont les mêmes champs `id` et `interval`, ce qui serait typique d'une analyse où nous voulons regrouper des réplicats techniques mais garder des échantillons significativement différents séparés.

Pour ce faire, nous devons séparer nos variables de regroupement afin de pouvoir les utiliser de manière isolée.

La première étape est similaire à ce que nous avons fait dans la section précédente. Nous devons isoler notre variable de regroupement comme premier élément du tuple. Rappelez-vous, notre premier élément est actuellement une map des champs `id`, `repeat` et `interval` :

```groovy title="main.nf" linenums="1"
{
  "id": "sampleA",
  "repeat": "1",
  "interval": "chr1"
}
```

Nous pouvons réutiliser la méthode `subMap` d'avant pour isoler nos champs `id` et `interval` de la map. Comme avant, nous utiliserons l'opérateur `map` pour appliquer la méthode `subMap` au premier élément du tuple pour chaque échantillon.

=== "Après"

    ```groovy title="main.nf" linenums="20" hl_lines="11-19"
        ch_combined_samples = ch_joined_samples
            .combine(ch_intervals)
            .map { grouping_key, normal, tumor, interval ->
                [
                    grouping_key + [interval: interval],
                    normal,
                    tumor
                ]
            }

        ch_grouped_samples = ch_combined_samples
            .map { grouping_key, normal, tumor ->
                [
                    grouping_key.subMap('id', 'interval'),
                    normal,
                    tumor
                ]
              }
              .view()
    ```

=== "Avant"

    ```groovy title="main.nf" linenums="20" hl_lines="10"
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

Exécutons-le à nouveau et vérifions le contenu du canal :

```bash
nextflow run main.nf
```

??? success "Sortie de la commande"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `main.nf` [hopeful_brenner] DSL2 - revision: 7f4f7fea76

    [[id:patientA, interval:chr1], patientA_rep1_normal.bam, patientA_rep1_tumor.bam]
    [[id:patientA, interval:chr2], patientA_rep1_normal.bam, patientA_rep1_tumor.bam]
    [[id:patientA, interval:chr3], patientA_rep1_normal.bam, patientA_rep1_tumor.bam]
    [[id:patientA, interval:chr1], patientA_rep2_normal.bam, patientA_rep2_tumor.bam]
    [[id:patientA, interval:chr2], patientA_rep2_normal.bam, patientA_rep2_tumor.bam]
    [[id:patientA, interval:chr3], patientA_rep2_normal.bam, patientA_rep2_tumor.bam]
    [[id:patientB, interval:chr1], patientB_rep1_normal.bam, patientB_rep1_tumor.bam]
    [[id:patientB, interval:chr2], patientB_rep1_normal.bam, patientB_rep1_tumor.bam]
    [[id:patientB, interval:chr3], patientB_rep1_normal.bam, patientB_rep1_tumor.bam]
    [[id:patientC, interval:chr1], patientC_rep1_normal.bam, patientC_rep1_tumor.bam]
    [[id:patientC, interval:chr2], patientC_rep1_normal.bam, patientC_rep1_tumor.bam]
    [[id:patientC, interval:chr3], patientC_rep1_normal.bam, patientC_rep1_tumor.bam]
    ```

Nous pouvons voir que nous avons réussi à isoler les champs `id` et `interval`, mais pas encore regroupé les échantillons.

!!! note

    Nous écartons le champ `replicate` ici. C'est parce que nous n'en avons pas besoin pour le traitement en aval ultérieur. Après avoir complété ce tutoriel, voyez si vous pouvez l'inclure sans affecter le regroupement ultérieur !

Regroupons maintenant les échantillons par ce nouvel élément de regroupement, en utilisant l'[opérateur `groupTuple`](https://www.nextflow.io/docs/latest/operator.html#grouptuple).

=== "Après"

    ```groovy title="main.nf" linenums="30" hl_lines="9"
        ch_grouped_samples = ch_combined_samples
            .map { grouping_key, normal, tumor ->
                [
                    grouping_key.subMap('id', 'interval'),
                    normal,
                    tumor
                ]
              }
              .groupTuple()
              .view()
    ```

=== "Avant"

    ```groovy title="main.nf" linenums="30"
        ch_grouped_samples = ch_combined_samples
            .map { grouping_key, normal, tumor ->
                [
                    grouping_key.subMap('id', 'interval'),
                    normal,
                    tumor
                ]
              }
              .view()
    ```

C'est tout ce qu'il y a à faire ! Nous avons juste ajouté une seule ligne de code. Voyons ce qui se passe lorsque nous l'exécutons :

```bash
nextflow run main.nf
```

??? success "Sortie de la commande"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `main.nf` [friendly_jang] DSL2 - revision: a1bee1c55d

    [[id:patientA, interval:chr1], [patientA_rep1_normal.bam, patientA_rep2_normal.bam], [patientA_rep1_tumor.bam, patientA_rep2_tumor.bam]]
    [[id:patientA, interval:chr2], [patientA_rep1_normal.bam, patientA_rep2_normal.bam], [patientA_rep1_tumor.bam, patientA_rep2_tumor.bam]]
    [[id:patientA, interval:chr3], [patientA_rep1_normal.bam, patientA_rep2_normal.bam], [patientA_rep1_tumor.bam, patientA_rep2_tumor.bam]]
    [[id:patientB, interval:chr1], [patientB_rep1_normal.bam], [patientB_rep1_tumor.bam]]
    [[id:patientB, interval:chr2], [patientB_rep1_normal.bam], [patientB_rep1_tumor.bam]]
    [[id:patientB, interval:chr3], [patientB_rep1_normal.bam], [patientB_rep1_tumor.bam]]
    [[id:patientC, interval:chr1], [patientC_rep1_normal.bam], [patientC_rep1_tumor.bam]]
    [[id:patientC, interval:chr2], [patientC_rep1_normal.bam], [patientC_rep1_tumor.bam]]
    [[id:patientC, interval:chr3], [patientC_rep1_normal.bam], [patientC_rep1_tumor.bam]]
    ```

Notez que notre structure de données a changé et qu'au sein de chaque élément de canal les fichiers sont maintenant contenus dans des tuples comme `[patientA_rep1_normal.bam, patientA_rep2_normal.bam]`. C'est parce que lorsque nous utilisons `groupTuple`, Nextflow combine les fichiers uniques pour chaque échantillon d'un groupe. C'est important à retenir lors de la gestion des données en aval.

!!! note

    [`transpose`](https://www.nextflow.io/docs/latest/reference/operator.html#transpose) est l'opposé de groupTuple. Il dépaquette les éléments dans un canal et les aplatit. Essayez d'ajouter `transpose` et d'annuler le regroupement que nous avons effectué ci-dessus !

### À retenir

Dans cette section, vous avez appris :

- **Regrouper des échantillons liés** : Comment utiliser `groupTuple` pour agréger des échantillons par attributs communs
- **Isoler les clés de regroupement** : Comment utiliser `subMap` pour extraire des champs spécifiques pour le regroupement
- **Gérer les structures de données regroupées** : Comment travailler avec la structure imbriquée créée par `groupTuple`
- **Gestion des réplicats techniques** : Comment regrouper des échantillons qui partagent les mêmes conditions expérimentales

---

## Résumé

Dans cette quête secondaire, vous avez appris comment séparer et regrouper des données en utilisant des canaux.

En modifiant les données au fur et à mesure qu'elles circulent dans le pipeline, vous pouvez construire un pipeline évolutif sans utiliser de boucles ou d'instructions while, offrant plusieurs avantages par rapport aux approches plus traditionnelles :

- Nous pouvons évoluer vers autant ou aussi peu d'entrées que nous le souhaitons sans code supplémentaire
- Nous nous concentrons sur la gestion du flux de données à travers le pipeline, au lieu de l'itération
- Nous pouvons être aussi complexes ou simples que nécessaire
- Le pipeline devient plus déclaratif, se concentrant sur ce qui devrait se passer plutôt que sur comment cela devrait se passer
- Nextflow optimisera l'exécution pour nous en exécutant des opérations indépendantes en parallèle

Maîtriser ces opérations de canaux vous permettra de construire des pipelines flexibles et évolutifs qui gèrent des relations de données complexes sans recourir à des boucles ou à la programmation itérative, permettant à Nextflow d'optimiser l'exécution et de paralléliser automatiquement les opérations indépendantes.

### Modèles clés

1.  **Créer des données d'entrée structurées :** Partir d'un fichier CSV avec des maps de métadonnées (en s'appuyant sur les modèles de [Métadonnées dans les workflows](./metadata.md))

    ```groovy
    ch_samples = channel.fromPath("./data/samplesheet.csv")
        .splitCsv(header: true)
        .map{ row ->
          [[id:row.id, repeat:row.repeat, type:row.type], row.bam]
        }
    ```

2.  **Séparer les données en canaux séparés :** Nous avons utilisé `filter` pour diviser les données en flux indépendants selon le champ `type`

    ```groovy
    channel.filter { it.type == 'tumor' }
    ```

3.  **Joindre des échantillons appariés :** Nous avons utilisé `join` pour recombiner des échantillons liés selon les champs `id` et `repeat`

    - Joindre deux canaux par clé (premier élément du tuple)

    ```groovy
    tumor_ch.join(normal_ch)
    ```

    - Extraire la clé de jointure et joindre par cette valeur

    ```groovy
    tumor_ch.map { meta, file -> [meta.id, meta, file] }
        .join(
          normal_ch.map { meta, file -> [meta.id, meta, file] }
        )
    ```

    - Joindre sur plusieurs champs en utilisant subMap

    ```groovy
    tumor_ch.map { meta, file -> [meta.subMap(['id', 'repeat']), meta, file] }
        .join(
          normal_ch.map { meta, file -> [meta.subMap(['id', 'repeat']), meta, file] }
        )
    ```

4.  **Distribuer à travers des intervalles :** Nous avons utilisé `combine` pour créer des produits cartésiens d'échantillons avec des intervalles génomiques pour le traitement parallèle.

    ```groovy
    samples_ch.combine(intervals_ch)
    ```

5.  **Agréger par clés de regroupement :** Nous avons utilisé `groupTuple` pour regrouper par le premier élément de chaque tuple, collectant ainsi des échantillons partageant les champs `id` et `interval` et fusionnant les réplicats techniques.

    ```groovy
    channel.groupTuple()
    ```

6.  **Optimiser la structure de données :** Nous avons utilisé `subMap` pour extraire des champs spécifiques et créé une closure nommée pour rendre les transformations réutilisables.

    - Extraire des champs spécifiques d'une map

    ```groovy
    meta.subMap(['id', 'repeat'])
    ```

    - Utiliser une closure nommée pour des transformations réutilisables

    ```groovy
    getSampleIdAndReplicate = { meta, file -> [meta.subMap(['id', 'repeat']), file] }
    channel.map(getSampleIdAndReplicate)
    ```

### Ressources supplémentaires

- [filter](https://www.nextflow.io/docs/latest/operator.html#filter)
- [map](https://www.nextflow.io/docs/latest/operator.html#map)
- [join](https://www.nextflow.io/docs/latest/operator.html#join)
- [groupTuple](https://www.nextflow.io/docs/latest/operator.html#grouptuple)
- [combine](https://www.nextflow.io/docs/latest/operator.html#combine)

---

## Et ensuite ?

Retournez au [menu des Quêtes secondaires](./index.md) ou cliquez sur le bouton en bas à droite de la page pour passer au sujet suivant de la liste.
