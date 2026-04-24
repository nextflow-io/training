# Métadonnées et meta maps

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Traduction assistée par IA - [en savoir plus et suggérer des améliorations](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

Dans toute analyse scientifique, nous travaillons rarement avec de simples fichiers de données brutes.
Chaque fichier est accompagné d'informations supplémentaires : ce qu'il est, d'où il vient, et ce qui le rend particulier.
Ces informations supplémentaires, c'est ce que nous appelons les métadonnées.

Les métadonnées sont des données qui décrivent d'autres données.
Elles permettent de suivre des détails importants sur les fichiers et les conditions expérimentales, et aident à adapter les analyses aux caractéristiques uniques de chaque jeu de données.

Pensez-y comme à un catalogue de bibliothèque : tandis que les livres contiennent le contenu réel (données brutes), les fiches de catalogue fournissent des informations essentielles sur chaque livre — quand il a été publié, qui l'a écrit, où le trouver (métadonnées).
Dans les pipelines Nextflow, les métadonnées peuvent être utilisées pour :

- Suivre les informations propres à chaque fichier tout au long du workflow
- Configurer les processus en fonction des caractéristiques des fichiers
- Regrouper des fichiers liés pour une analyse conjointe

### Objectifs d'apprentissage

Dans cette quête secondaire, nous allons explorer comment gérer les métadonnées dans les workflows.
En partant d'une feuille de données simple (souvent appelée samplesheet en bioinformatique) contenant des informations de base sur les fichiers, vous apprendrez à :

- Lire et analyser les métadonnées de fichiers à partir de fichiers CSV
- Créer et manipuler des meta maps
- Ajouter de nouveaux champs de métadonnées pendant l'exécution du workflow
- Utiliser les métadonnées pour personnaliser le comportement des processus

Ces compétences vous aideront à construire des pipelines plus robustes et flexibles, capables de gérer des relations complexes entre fichiers et des exigences de traitement variées.

### Prérequis

Avant de vous lancer dans cette quête secondaire, vous devriez :

- Avoir suivi le tutoriel [Hello Nextflow](../hello_nextflow/README.md) ou un cours équivalent pour débutant·es.
- Être à l'aise avec les concepts et mécanismes de base de Nextflow (processus, canaux, opérateurs)

---

## 0. Premiers pas

#### Ouvrir le codespace de formation

Si vous ne l'avez pas encore fait, assurez-vous d'ouvrir l'environnement de formation comme décrit dans la [Configuration de l'environnement](../envsetup/index.md).

[![Open in GitHub Codespaces](https://github.com/codespaces/badge.svg)](https://codespaces.new/nextflow-io/training?quickstart=1&ref=master)

#### Se déplacer dans le répertoire du projet

Déplaçons-nous dans le répertoire où se trouvent les fichiers de ce tutoriel.

```bash
cd side-quests/metadata
```

Vous pouvez configurer VSCode pour qu'il se concentre sur ce répertoire :

```bash
code .
```

#### Examiner les fichiers

Vous trouverez un fichier de workflow principal et un répertoire `data` contenant une feuille de données et quelques fichiers de données.

??? abstract "Contenu du répertoire"

    ```console
    .
    ├── data
    │   ├── bonjour.txt
    │   ├── ciao.txt
    │   ├── guten_tag.txt
    │   ├── hallo.txt
    │   ├── hello.txt
    │   ├── hola.txt
    │   ├── salut.txt
    │   └── datasheet.csv
    ├── main.nf
    └── nextflow.config
    ```

Le workflow dans le fichier `main.nf` est une ébauche que vous allez progressivement développer en un workflow pleinement fonctionnel.

La feuille de données liste les chemins vers les fichiers de données et quelques métadonnées associées, organisées en 3 colonnes :

- `id` : explicite, un identifiant attribué au fichier
- `character` : un nom de personnage, que nous utiliserons plus tard pour dessiner différentes créatures
- `data` : chemins vers des fichiers `.txt` contenant des salutations dans différentes langues

```console title="datasheet.csv"
id,character,recording
sampleA,squirrel,/workspaces/training/side-quests/metadata/data/bonjour.txt
sampleB,tux,/workspaces/training/side-quests/metadata/data/guten_tag.txt
sampleC,sheep,/workspaces/training/side-quests/metadata/data/hallo.txt
sampleD,turkey,/workspaces/training/side-quests/metadata/data/hello.txt
sampleE,stegosaurus,/workspaces/training/side-quests/metadata/data/hola.txt
sampleF,moose,/workspaces/training/side-quests/metadata/data/salut.txt
sampleG,turtle,/workspaces/training/side-quests/metadata/data/ciao.txt
```

Chaque fichier de données contient du texte de salutation dans l'une des cinq langues (fr : français, de : allemand, es : espagnol, it : italien, en : anglais).

Nous vous fournirons également un outil d'analyse linguistique conteneurisé appelé `langid`.

#### Examiner l'exercice

Votre défi est d'écrire un workflow Nextflow qui va :

1. **Identifier** automatiquement la langue dans chaque fichier
2. **Regrouper** les fichiers par famille linguistique (langues germaniques vs langues romanes)
3. **Personnaliser** le traitement de chaque fichier en fonction de sa langue et de ses métadonnées
4. **Organiser** les sorties par groupe linguistique

Cela représente un schéma de workflow typique où les métadonnées propres à chaque fichier guident les décisions de traitement ; exactement le type de problème que les meta maps résolvent élégamment.

#### Liste de vérification

Vous pensez être prêt·e à vous lancer ?

- [ ] Je comprends l'objectif de ce cours et ses prérequis
- [ ] Mon codespace est opérationnel
- [ ] J'ai défini mon répertoire de travail de manière appropriée
- [ ] Je comprends l'exercice

Si vous pouvez cocher toutes les cases, vous êtes prêt·e à commencer.

---

## 1. Charger les métadonnées depuis une feuille de données

Ouvrez le fichier de workflow `main.nf` pour examiner l'ébauche de workflow que nous vous fournissons comme point de départ.

```groovy title="main.nf" linenums="1"
#!/usr/bin/env nextflow

workflow  {

    ch_datasheet = channel.fromPath("./data/datasheet.csv")

}
```

Vous pouvez voir que nous avons mis en place une factory de canal basique pour charger la feuille de données exemple en tant que fichier, mais cela ne lira pas encore le contenu du fichier.
Commençons par ajouter cela.

### 1.1. Lire le contenu avec `splitCsv`

Nous devons choisir un opérateur qui analysera le contenu du fichier de manière appropriée avec un minimum d'effort de notre part.
Puisque notre feuille de données est au format CSV, c'est le rôle de l'opérateur [`splitCsv`](https://www.nextflow.io/docs/latest/reference/operator.html#splitcsv), qui charge chaque ligne du fichier comme un élément dans le canal.

Effectuez les modifications suivantes pour ajouter une opération `splitCsv()` au code de construction du canal, ainsi qu'une opération `view()` pour vérifier que le contenu du fichier est correctement chargé dans le canal.

=== "Après"

    ```groovy title="main.nf" linenums="3" hl_lines="4-5"
    workflow  {

        ch_datasheet = channel.fromPath("./data/datasheet.csv")
            .splitCsv(header: true)
            .view()

    }
    ```

=== "Avant"

    ```groovy title="main.nf" linenums="3"
    workflow  {

        ch_datasheet = channel.fromPath("./data/datasheet.csv")

    }
    ```

Notez que nous utilisons l'option `header: true` pour indiquer à Nextflow de lire la première ligne du fichier CSV comme ligne d'en-tête.

Voyons ce que cela produit, voulez-vous ?
Exécutez le workflow :

```bash
nextflow run main.nf
```

??? success "Sortie de la commande"

    ```console
     N E X T F L O W   ~  version 25.10.2

    Launching `main.nf` [exotic_albattani] DSL2 - revision: c0d03cec83

    [id:sampleA, character:squirrel, recording:/workspaces/training/side-quests/metadata/data/bonjour.txt]
    [id:sampleB, character:tux, recording:/workspaces/training/side-quests/metadata/data/guten_tag.txt]
    [id:sampleC, character:sheep, recording:/workspaces/training/side-quests/metadata/data/hallo.txt]
    [id:sampleD, character:turkey, recording:/workspaces/training/side-quests/metadata/data/hello.txt]
    [id:sampleE, character:stegosaurus, recording:/workspaces/training/side-quests/metadata/data/hola.txt]
    [id:sampleF, character:moose, recording:/workspaces/training/side-quests/metadata/data/salut.txt]
    [id:sampleG, character:turtle, recording:/workspaces/training/side-quests/metadata/data/ciao.txt]
    ```

Nous pouvons voir que l'opérateur a construit une map de paires clé-valeur pour chaque ligne du fichier CSV, avec les en-têtes de colonnes comme clés pour les valeurs correspondantes.

Chaque entrée de la map correspond à une colonne de notre feuille de données :

- `id`
- `character`
- `recording`

C'est parfait ! Cela facilite l'accès à des champs spécifiques de chaque fichier.
Par exemple, nous pourrions accéder à l'identifiant du fichier avec `id` ou au chemin du fichier txt avec `recording`.

??? info "(Optionnel) En savoir plus sur les maps"

    En Groovy, le langage de programmation sur lequel Nextflow est construit, une map est une structure de données clé-valeur similaire aux dictionnaires en Python, aux objets en JavaScript, ou aux hashes en Ruby.

    Voici un script exécutable qui montre comment définir une map et accéder à son contenu en pratique :

    ```groovy title="examples/map_demo.nf"
    #!/usr/bin/env nextflow

    // Créer une map simple
    def my_map = [id:'sampleA', character:'squirrel']

    // Afficher la map entière
    println "map: ${my_map}"

    // Accéder aux valeurs individuelles avec la notation pointée
    println "id: ${my_map.id}"
    println "character: ${my_map.character}"
    ```

    Même s'il n'a pas de bloc `workflow` à proprement parler, Nextflow peut l'exécuter comme s'il s'agissait d'un workflow :

    ```bash
    nextflow run examples/map_demo.nf
    ```

    Et voici ce que vous pouvez vous attendre à voir dans la sortie :

    ```console title="Output"
     N E X T F L O W   ~  version 25.10.2

    Launching `map_demo.nf` [cheesy_plateau] DSL2 - revision: fae5b8496e

    map: [id:sampleA, character:squirrel]
    id: sampleA
    character: squirrel
    ```

### 1.2. Sélectionner des champs spécifiques avec `map`

Supposons que nous voulions accéder à la colonne `character` de la feuille de données et l'afficher.
Nous pouvons utiliser l'opérateur Nextflow `map` pour itérer sur chaque élément de notre canal et sélectionner spécifiquement l'entrée `character` de l'objet map.

Effectuez les modifications suivantes dans le workflow :

=== "Après"

    ```groovy title="main.nf" linenums="3" hl_lines="5-7"
    workflow  {

        ch_datasheet = channel.fromPath("./data/datasheet.csv")
            .splitCsv(header: true)
            .map{ row ->
                row.character
            }
            .view()

    }
    ```

=== "Avant"

    ```groovy title="main.nf" linenums="3"
    workflow  {

        ch_datasheet = channel.fromPath("./data/datasheet.csv")
            .splitCsv(header: true)
            .view()

    }
    ```

Exécutez à nouveau le workflow :

```bash
nextflow run main.nf
```

??? success "Sortie de la commande"

    ```console
     N E X T F L O W   ~  version 25.10.2

    Launching `main.nf` [exotic_albattani] DSL2 - revision: c0d03cec83

    squirrel
    tux
    sheep
    turkey
    stegosaurus
    moose
    turtle
    ```

Parfait ! Nous avons tiré parti de la structure map dérivée de notre feuille de données pour accéder aux valeurs des colonnes individuelles pour chaque ligne.

Maintenant que nous avons réussi à lire la feuille de données et que nous avons accès aux données de chaque ligne, nous pouvons commencer à implémenter la logique de notre pipeline.

### 1.3. Organiser les métadonnées dans une 'meta map'

Dans l'état actuel du workflow, les fichiers d'entrée (sous la clé `recording`) et les métadonnées associées (`id`, `character`) sont tous sur le même plan, comme s'ils étaient tous dans un grand sac.
La conséquence pratique est que chaque processus qui consomme ce canal devrait être configuré en tenant compte de cette structure :

```groovy
    input:
    tuple val(id), val(character), file(recording)
```

C'est acceptable tant que le nombre de colonnes dans la feuille de données ne change pas.
Cependant, si vous ajoutez ne serait-ce qu'une seule colonne à la feuille de données, la forme du canal ne correspondra plus à ce que le processus attend, et le workflow produira des erreurs.
Cela rend également le processus difficile à partager avec d'autres personnes qui pourraient avoir des données d'entrée légèrement différentes, et vous pourriez finir par devoir coder en dur des variables dans le processus qui ne sont pas nécessaires dans le bloc script.

Pour éviter ce problème, nous devons trouver un moyen de maintenir la structure du canal cohérente, quel que soit le nombre de colonnes que contient la feuille de données.

Nous pouvons le faire en regroupant toutes les métadonnées dans un élément au sein du tuple, que nous appellerons la meta map (ou simplement 'meta map').

Effectuez les modifications suivantes dans l'opération `map` :

=== "Après"

    ```groovy title="main.nf" linenums="5" hl_lines="4"
        ch_datasheet = channel.fromPath("./data/datasheet.csv")
            .splitCsv(header: true)
            .map { row ->
                    [ [id: row.id, character: row.character], row.recording ]
            }
            .view()
    ```

=== "Avant"

    ```groovy title="main.nf" linenums="5" hl_lines="4"
        ch_datasheet = channel.fromPath("./data/datasheet.csv")
            .splitCsv(header: true)
            .map{ row ->
                row.character
            }
            .view()
    ```

Nous avons restructuré les éléments de notre canal en un tuple composé de deux éléments : la meta map et l'objet fichier correspondant.

Exécutons le workflow :

```bash
nextflow run main.nf
```

??? success "Sortie de la commande"

    ```console title="View meta map"
     N E X T F L O W   ~  version 25.10.2

    Launching `main.nf` [lethal_booth] DSL2 - revision: 0d8f844c07

    [[id:sampleA, character:squirrel], /workspaces/training/side-quests/metadata/data/bonjour.txt]
    [[id:sampleB, character:tux], /workspaces/training/side-quests/metadata/data/guten_tag.txt]
    [[id:sampleC, character:sheep], /workspaces/training/side-quests/metadata/data/hallo.txt]
    [[id:sampleD, character:turkey], /workspaces/training/side-quests/metadata/data/hello.txt]
    [[id:sampleE, character:stegosaurus], /workspaces/training/side-quests/metadata/data/hola.txt]
    [[id:sampleF, character:moose], /workspaces/training/side-quests/metadata/data/salut.txt]
    [[id:sampleG, character:turtle], /workspaces/training/side-quests/metadata/data/ciao.txt]
    ```

Désormais, chaque élément du canal contient d'abord la meta map et ensuite l'objet fichier correspondant :

```console title="Example output structure"
[
  [id:sampleA, character:squirrel],
  /workspaces/training/side-quests/metadata/data/bonjour.txt
]
```

Par conséquent, l'ajout de colonnes dans la feuille de données rendra davantage de métadonnées disponibles dans la map `meta`, mais ne changera pas la forme du canal.
Cela nous permet d'écrire des processus qui consomment le canal sans avoir à coder en dur les éléments de métadonnées dans la spécification des entrées :

```groovy title="Syntax example"
    input:
    tuple val(meta), file(recording)
```

Il s'agit d'une convention largement utilisée pour organiser les métadonnées dans les workflows Nextflow.

### À retenir

Dans cette section, vous avez appris :

- **Pourquoi les métadonnées sont importantes :** Conserver les métadonnées avec vos données préserve les informations importantes sur les fichiers tout au long du workflow.
- **Comment lire une feuille de données :** Utiliser `splitCsv` pour lire des fichiers CSV avec des informations d'en-tête et transformer les lignes en données structurées
- **Comment créer une meta map :** Séparer les métadonnées des données de fichiers en utilisant la structure tuple `[ [id:valeur, ...], fichier ]`

---

## 2. Manipuler les métadonnées

Maintenant que nos métadonnées sont chargées, faisons quelque chose avec !

Nous allons utiliser un outil appelé [`langid`](https://github.com/saffsd/langid.py) pour identifier la langue contenue dans le fichier d'enregistrement de chaque créature.
L'outil est pré-entraîné sur un ensemble de langues, et à partir d'un extrait de texte, il produira une prédiction de langue et un score de probabilité associé, tous deux vers `stdout`.

### 2.1. Importer le processus et examiner le code

Nous vous fournissons un module de processus pré-écrit appelé `IDENTIFY_LANGUAGE` qui encapsule l'outil `langid`, vous n'avez donc qu'à ajouter une instruction include avant le bloc workflow.

Effectuez la modification suivante dans le workflow :

=== "Après"

    ```groovy title="main.nf" linenums="1" hl_lines="3"
    #!/usr/bin/env nextflow

    include { IDENTIFY_LANGUAGE } from './modules/langid.nf'

    workflow {
    ```

=== "Avant"

    ```groovy title="main.nf" linenums="1"
    #!/usr/bin/env nextflow

    workflow {
    ```

Vous pouvez ouvrir le fichier module pour examiner son code :

```groovy title="modules/langid.nf" linenums="1" hl_lines="9 12"
#!/usr/bin/env nextflow

// Utiliser langid pour prédire la langue de chaque fichier d'entrée
process IDENTIFY_LANGUAGE {

    container 'community.wave.seqera.io/library/pip_langid:b2269f456a5629ff'

    input:
    tuple val(meta), path(file)

    output:
    tuple val(meta), path(file), stdout

    script:
    """
    langid < ${file} -l en,de,fr,es,it | sed -E "s/.*\\('([a-z]+)'.*/\\1/" | tr -d '\\n'
    """
}
```

Comme vous pouvez le voir, la définition des entrées utilise la même structure `tuple val(meta), path(file)` que nous venons d'appliquer à notre canal d'entrée.

La définition des sorties est structurée comme un tuple avec une structure similaire à celle des entrées, sauf qu'elle contient également `stdout` comme troisième élément.
Ce schéma `tuple val(meta), path(file), <sortie>` maintient les métadonnées associées à la fois aux données d'entrée et aux sorties au fur et à mesure qu'elles circulent dans le pipeline.

Notez que nous utilisons ici le qualificateur de sortie [`stdout`](https://www.nextflow.io/docs/latest/process.html#outputs) de Nextflow, car l'outil affiche sa sortie directement sur la console plutôt que d'écrire un fichier ; et nous utilisons `sed` en ligne de commande pour supprimer le score de probabilité, nettoyer la chaîne en supprimant les caractères de nouvelle ligne, et ne retourner que la prédiction de langue.

### 2.2. Ajouter un appel à `IDENTIFY_LANGUAGE`

Maintenant que le processus est disponible pour le workflow, nous pouvons ajouter un appel au processus `IDENTIFY_LANGUAGE` pour l'exécuter sur le canal de données.

Effectuez les modifications suivantes dans le workflow :

=== "Après"

    ```groovy title="main.nf" linenums="7" hl_lines="7-9"
        ch_datasheet = channel.fromPath("./data/datasheet.csv")
            .splitCsv(header: true)
            .map { row ->
                [[id: row.id, character: row.character], row.recording]
            }

        // Exécuter langid pour identifier la langue de chaque salutation
        IDENTIFY_LANGUAGE(ch_datasheet)
        IDENTIFY_LANGUAGE.out.view()
    ```

=== "Avant"

    ```groovy title="main.nf" linenums="7" hl_lines="6"
        ch_datasheet = channel.fromPath("./data/datasheet.csv")
            .splitCsv(header: true)
            .map { row ->
                [[id: row.id, character: row.character], row.recording]
            }
            .view()
    ```

Notez que nous avons supprimé l'opération `.view()` originale dans la construction du canal.

Nous pouvons maintenant exécuter le workflow :

```bash
nextflow run main.nf
```

??? success "Sortie de la commande"

    ```console
     N E X T F L O W   ~  version 25.10.2

    Launching `main.nf` [voluminous_mcnulty] DSL2 - revision: f9bcfebabb

    executor >  local (7)
    [4e/f722fe] IDENTIFY_LANGUAGE (7) [100%] 7 of 7 ✔
    [[id:sampleA, character:squirrel], /workspaces/training/side-quests/metadata/work/eb/f7148ebdd898fbe1136bec6a714acb/bonjour.txt, fr]
    [[id:sampleB, character:tux], /workspaces/training/side-quests/metadata/work/16/71d72410952c22cd0086d9bca03680/guten_tag.txt, de]
    [[id:sampleD, character:turkey], /workspaces/training/side-quests/metadata/work/c4/b7562adddc1cc0b7d414ec45d436eb/hello.txt, en]
    [[id:sampleC, character:sheep], /workspaces/training/side-quests/metadata/work/ea/04f5d979429e4455e14b9242fb3b45/hallo.txt, de]
    [[id:sampleF, character:moose], /workspaces/training/side-quests/metadata/work/5a/6c2b84bf8fadb98e28e216426be079/salut.txt, fr]
    [[id:sampleE, character:stegosaurus], /workspaces/training/side-quests/metadata/work/af/ee7c69bcab891c40d0529305f6b9e7/hola.txt, es]
    [[id:sampleG, character:turtle], /workspaces/training/side-quests/metadata/work/4e/f722fe47271ba7ebcd69afa42964ca/ciao.txt, it]
    ```

Excellent ! Nous avons maintenant une prédiction de la langue parlée par chaque personnage.

Et comme mentionné précédemment, nous avons également inclus le fichier d'entrée et la meta map dans la sortie, ce qui signifie que les deux restent associés aux nouvelles informations que nous venons de produire.
Cela s'avérera utile à l'étape suivante.

!!! note "Note"

    Plus généralement, ce schéma consistant à maintenir la meta map associée aux résultats facilite l'association de résultats liés qui partagent les mêmes identifiants.

    Comme vous l'avez déjà appris, vous ne pouvez pas vous fier à l'ordre des éléments dans les canaux pour faire correspondre les résultats entre eux.
    Vous devez plutôt utiliser des clés pour associer correctement les données, et les meta maps fournissent une structure idéale à cet effet.

    Nous explorons ce cas d'utilisation en détail dans la quête secondaire [Splitting & Grouping](../splitting_and_grouping/).

### 2.3. Enrichir les métadonnées avec les sorties des processus

Étant donné que les résultats que nous venons de produire constituent en eux-mêmes une forme de métadonnées sur le contenu des fichiers, il serait utile de les ajouter à notre meta map.

Cependant, nous ne voulons pas modifier la meta map existante en place.
D'un point de vue technique, il est _possible_ de le faire, mais c'est risqué.

Nous allons donc plutôt créer une nouvelle meta map contenant le contenu de la meta map existante plus une nouvelle paire clé-valeur `lang: lang_id` contenant les nouvelles informations, en utilisant l'opérateur `+` (une fonctionnalité de Groovy).
Et nous combinerons cela avec une opération [`map`](https://www.nextflow.io/docs/latest/operator.html#map) pour remplacer l'ancienne map par la nouvelle.

Voici les modifications à apporter au workflow :

=== "Après"

    ```groovy title="main.nf" linenums="13" hl_lines="3-7"
        // Exécuter langid pour identifier la langue de chaque salutation
        IDENTIFY_LANGUAGE(ch_datasheet)
        IDENTIFY_LANGUAGE.out
            .map { meta, file, lang_id ->
                [meta + [lang: lang_id], file]
            }
            .view()
    ```

=== "Avant"

    ```groovy title="main.nf" linenums="13" hl_lines="3"
        // Exécuter langid pour identifier la langue de chaque salutation
        IDENTIFY_LANGUAGE(ch_datasheet)
        IDENTIFY_LANGUAGE.out.view()
    ```

Si vous n'êtes pas encore familier·ère avec l'opérateur `+`, ou si cela vous semble confus, prenez quelques minutes pour parcourir l'explication détaillée ci-dessous.

??? info "Création de la nouvelle meta map avec l'opérateur `+`"

    **Premièrement, vous devez savoir que nous pouvons fusionner le contenu de deux maps en utilisant l'opérateur Groovy `+`.**

    Supposons que nous ayons les maps suivantes :

    ```groovy
    map1 = [id: 'sampleA', character: 'squirrel']
    map2 = [lang: 'fr']
    ```

    Nous pouvons les fusionner ainsi :

    ```groovy
    new_map = map1 + map2
    ```

    Le contenu de `new_map` sera :

    ```groovy
    [id: 'sampleA', character: 'squirrel', lang: 'fr']
    ```

    Parfait !

    **Mais que faire si vous devez ajouter un champ qui ne fait pas déjà partie d'une map ?**

    Supposons que vous repartiez de `map1`, mais que la prédiction de langue ne soit pas dans sa propre map (il n'y a pas de `map2`).
    Elle est plutôt stockée dans une variable appelée `lang_id`, et vous savez que vous voulez stocker sa valeur (`'fr'`) avec la clé `lang`.

    Vous pouvez en fait faire ce qui suit :

    ```groovy
    new_map = [map1 + [lang: lang_id]]
    ```

    Ici, `[lang: new_info]` crée une nouvelle map sans nom à la volée, et `map1 + ` fusionne `map1` avec la nouvelle map sans nom, produisant le même contenu de `new_map` qu'auparavant.

    Élégant, non ?

    **Transposons maintenant cela dans le contexte d'une opération Nextflow `channel.map()`.**

    Le code devient :

    ```groovy
    .map { map1, lang_id ->
        [map1 + [lang: lang_id]]
    }
    ```

    Cela fait ce qui suit :

    - `map1, lang_id ->` prend les deux éléments du tuple
    - `[map1 + [lang: lang_id]]` crée la nouvelle map comme détaillé ci-dessus

    La sortie est une seule map sans nom avec le même contenu que `new_map` dans notre exemple ci-dessus.
    Nous avons donc effectivement transformé :

    ```groovy
    [id: 'sampleA', character: 'squirrel'], 'fr'
    ```

    en :

    ```groovy
    [id: 'sampleA', character: 'squirrel', lang: 'fr']
    ```

    Vous pouvez voir que si nous remplaçons `map1` par `meta`, c'est essentiellement tout ce dont nous avons besoin pour ajouter la prédiction de langue à notre meta map dans notre workflow.

    Sauf pour une chose !

    Dans le cas de notre workflow, **nous devons également tenir compte de la présence de l'objet `file` dans le tuple**, qui est composé de `meta, file, lang_id`.

    Le code devient donc :

    ```groovy
    .map { meta, file, lang_id ->
        [meta + [lang: lang_id], file]
    }
    ```

    Si vous avez du mal à comprendre pourquoi le `file` semble se déplacer dans l'opération `map`, imaginez qu'au lieu de `[meta + [lang: lang_id], file]`, cette ligne se lise `[new_map, file]`.
    Cela devrait clarifier que nous laissons simplement le `file` à sa place d'origine en deuxième position dans le tuple. Nous avons juste pris la valeur `new_info` et l'avons intégrée dans la map qui est en première position.

    **Et cela nous ramène à la structure de canal `tuple val(meta), path(file)` !**

Une fois que vous êtes certain·e de comprendre ce que fait ce code, exécutez le workflow pour voir si cela a fonctionné :

```bash
nextflow run main.nf -resume
```

??? success "Sortie de la commande"

    ```console
     N E X T F L O W   ~  version 25.10.2

    Launching `main.nf` [cheeky_fermat] DSL2 - revision: d096281ee4

    [4e/f722fe] IDENTIFY_LANGUAGE (7) [100%] 7 of 7, cached: 7 ✔
    [[id:sampleA, character:squirrel, lang:fr], /workspaces/training/side-quests/metadata/work/eb/f7148ebdd898fbe1136bec6a714acb/bonjour.txt]
    [[id:sampleB, character:tux, lang:de], /workspaces/training/side-quests/metadata/work/16/71d72410952c22cd0086d9bca03680/guten_tag.txt]
    [[id:sampleC, character:sheep, lang:de], /workspaces/training/side-quests/metadata/work/ea/04f5d979429e4455e14b9242fb3b45/hallo.txt]
    [[id:sampleD, character:turkey, lang:en], /workspaces/training/side-quests/metadata/work/c4/b7562adddc1cc0b7d414ec45d436eb/hello.txt]
    [[id:sampleF, character:moose, lang:fr], /workspaces/training/side-quests/metadata/work/5a/6c2b84bf8fadb98e28e216426be079/salut.txt]
    [[id:sampleE, character:stegosaurus, lang:es], /workspaces/training/side-quests/metadata/work/af/ee7c69bcab891c40d0529305f6b9e7/hola.txt]
    [[id:sampleG, character:turtle, lang:it], /workspaces/training/side-quests/metadata/work/4e/f722fe47271ba7ebcd69afa42964ca/ciao.txt]
    ```

Oui, c'est bien ça !
Nous avons soigneusement réorganisé la sortie du processus de `meta, file, lang_id` de sorte que `lang_id` est maintenant l'une des clés dans la meta map, et les tuples du canal correspondent à nouveau au modèle `meta, file`.

<!-- TODO (future) Should we also show how to remove a key using subMap?! Or note where to find that. -->

### 2.4. Assigner un groupe linguistique avec des conditions

Maintenant que nous avons nos prédictions de langue, utilisons ces informations pour définir de nouveaux regroupements.

Dans nos données d'exemple, les langues utilisées par nos personnages peuvent être regroupées en langues germaniques (anglais, allemand) et langues romanes (français, espagnol, italien).
Il pourrait être utile d'avoir cette classification facilement disponible plus tard dans le pipeline, alors ajoutons ces informations dans la meta map.

Et bonne nouvelle, c'est encore un cas qui se prête parfaitement à l'utilisation de l'opérateur `map` !

Plus précisément, nous allons définir une variable appelée `lang_group`, utiliser une logique conditionnelle simple pour déterminer quelle valeur assigner à `lang_group` pour chaque donnée.

La syntaxe générale va ressembler à ceci :

```groovy
.map { meta, file ->

    // la logique conditionnelle définissant lang_group va ici

    [meta + [lang_group: lang_group], file]
}
```

Vous pouvez voir que c'est très similaire à l'opération de fusion de maps à la volée que nous avons utilisée à l'étape précédente.
Nous avons juste besoin d'écrire les instructions conditionnelles.

Voici la logique conditionnelle que nous voulons appliquer :

- Définir une variable appelée `lang_group` avec la valeur par défaut `'unknown'`.
- Si `lang` est soit l'allemand (`'de'`) soit l'anglais (`'en'`), changer `lang_group` en `germanic`.
- Sinon si `lang` est inclus dans une liste contenant le français (`'fr'`), l'espagnol (`'es'`) et l'italien (`'it'`), changer `lang_group` en `romance`.

Essayez de l'écrire vous-même si vous savez déjà comment écrire des instructions conditionnelles en Nextflow.

!!! tip "Astuce"

    Vous pouvez accéder à la valeur de `lang` dans l'opération map avec `meta.lang`.

Vous devriez finir par effectuer les modifications suivantes dans le workflow :

=== "Après"

    ```groovy title="main.nf" linenums="13" hl_lines="7-19"
        // Exécuter langid pour identifier la langue de chaque salutation
        IDENTIFY_LANGUAGE(ch_datasheet)
        IDENTIFY_LANGUAGE.out
            .map { meta, file, lang_id ->
                [meta + [lang: lang_id], file]
            }
            .map { meta, file ->

                def lang_group = "unknown"
                if (meta.lang.equals("de") || meta.lang.equals('en')) {
                    lang_group = "germanic"
                }
                else if (meta.lang in ["fr", "es", "it"]) {
                    lang_group = "romance"
                }

                [meta + [lang_group: lang_group], file]
            }
            .view()
    ```

=== "Avant"

    ```groovy title="main.nf" linenums="13" hl_lines="7"
        // Exécuter langid pour identifier la langue de chaque salutation
        IDENTIFY_LANGUAGE(ch_datasheet)
        IDENTIFY_LANGUAGE.out
            .map { meta, file, lang_id ->
                [meta + [lang: lang_id], file]
            }
            .view()
    ```

Voici les points clés :

- Nous utilisons `def lang_group = "unknown"` pour créer la variable `lang_group` avec la valeur par défaut `unknown`.
- Nous utilisons une structure `if {} else if {}` pour la logique conditionnelle, avec des tests `.equals()` alternatifs pour les deux langues germaniques, et un test d'existence dans une liste pour les trois langues romanes.
- Nous utilisons l'opération de fusion `meta + [lang_group:lang_group]` comme précédemment pour générer la meta map mise à jour.

<!-- TODO (future) Add note/links to relevant docs in additional resources section -->

Une fois que tout cela est clair, exécutez à nouveau le workflow pour voir le résultat :

```bash
nextflow run main.nf -resume
```

??? success "Sortie de la commande"

    ```console
     N E X T F L O W   ~  version 25.10.2

    Launching `main.nf` [wise_almeida] DSL2 - revision: 46778c3cd0

    [da/652cc6] IDENTIFY_LANGUAGE (7) [100%] 7 of 7, cached: 7 ✔
    [[id:sampleA, character:squirrel, lang:fr, lang_group:romance], /workspaces/training/side-quests/metadata/data/bonjour.txt]
    [[id:sampleB, character:tux, lang:de, lang_group:germanic], /workspaces/training/side-quests/metadata/data/guten_tag.txt]
    [[id:sampleC, character:sheep, lang:de, lang_group:germanic], /workspaces/training/side-quests/metadata/data/hallo.txt]
    [[id:sampleD, character:turkey, lang:en, lang_group:germanic], /workspaces/training/side-quests/metadata/data/hello.txt]
    [[id:sampleE, character:stegosaurus, lang:es, lang_group:romance], /workspaces/training/side-quests/metadata/data/hola.txt]
    [[id:sampleF, character:moose, lang:fr, lang_group:romance], /workspaces/training/side-quests/metadata/data/salut.txt]
    [[id:sampleG, character:turtle, lang:it, lang_group:romance], /workspaces/training/side-quests/metadata/data/ciao.txt]
    ```

Comme vous pouvez le voir, les éléments du canal maintiennent leur structure `[meta, file]`, mais la meta map inclut maintenant cette nouvelle classification.

### À retenir

Dans cette section, vous avez appris à :

- **Appliquer les métadonnées d'entrée aux canaux de sortie :** Copier les métadonnées de cette façon nous permet d'associer les résultats ultérieurement en fonction du contenu des métadonnées.
- **Créer des clés personnalisées :** Vous avez créé deux nouvelles clés dans votre meta map, en les fusionnant avec `meta + [nouvelle_clé:valeur]` dans la meta map existante. L'une basée sur une valeur calculée par un processus, et l'autre basée sur une condition définie dans l'opérateur `map`.

Ces techniques vous permettent d'associer des métadonnées nouvelles et existantes aux fichiers au fur et à mesure de votre progression dans le pipeline.
Même si vous n'utilisez pas les métadonnées dans un processus, maintenir la meta map associée aux données de cette façon facilite le regroupement de toutes les informations pertinentes.

---

## 3. Utiliser les informations de la meta map dans un processus

Maintenant que vous savez comment créer et mettre à jour la meta map, nous pouvons passer à la partie vraiment intéressante : utiliser réellement les métadonnées dans un processus.

Plus précisément, nous allons ajouter une deuxième étape à notre workflow pour dessiner chaque animal en art ASCII et lui faire dire le texte enregistré dans une bulle de dialogue.
Nous allons le faire en utilisant un outil appelé [`cowpy`](https://github.com/jeffbuttars/cowpy).

??? info "Que fait `cowpy` ?"

    `cowpy` est un outil en ligne de commande qui génère de l'art ASCII pour afficher des textes arbitraires de manière amusante.
    C'est une implémentation Python du classique outil [cowsay](https://en.wikipedia.org/wiki/Cowsay) de Tony Monroe.

    ```console
    cowpy "Hello Nextflow"
    ```

    ```console
    ______________________________________________________
    < Hello Nextflow >
    ------------------------------------------------------
        \   ^__^
          \  (oo)\_______
            (__)\       )\/\
              ||----w |
              ||     ||
    ```

    Optionnellement, vous pouvez sélectionner un personnage (ou 'cowacter') à utiliser à la place de la vache par défaut.

    ```console
    cowpy "Hello Nextflow" -c tux
    ```

    ```console
    __________________
    < Hello Nextflow >
    ------------------
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

Si vous avez suivi le cours Hello Nextflow, vous avez déjà vu cet outil en action.
Sinon, ne vous inquiétez pas ; nous couvrirons tout ce que vous devez savoir au fur et à mesure.

### 3.1. Importer le processus et examiner le code

Nous vous fournissons un module de processus pré-écrit appelé `COWPY` qui encapsule l'outil `cowpy`, vous n'avez donc qu'à ajouter une instruction include avant le bloc workflow.

Effectuez la modification suivante dans le workflow :

=== "Après"

    ```groovy title="main.nf" linenums="1" hl_lines="4"
    #!/usr/bin/env nextflow

    include { IDENTIFY_LANGUAGE } from './modules/langid.nf'
    include { COWPY } from './modules/cowpy.nf'

    workflow {
    ```

=== "Avant"

    ```groovy title="main.nf" linenums="1"
    #!/usr/bin/env nextflow

    include { IDENTIFY_LANGUAGE } from './modules/langid.nf'

    workflow {
    ```

Vous pouvez ouvrir le fichier module pour examiner son code :

```groovy title="modules/cowpy.nf" linenums="1"
#!/usr/bin/env nextflow

// Générer de l'art ASCII avec cowpy
process COWPY {

    publishDir "results/", mode: 'copy'

    container 'community.wave.seqera.io/library/cowpy:1.1.5--3db457ae1977a273'

    input:
    path input_file
    val character

    output:
    path "cowpy-${input_file}"

    script:
    """
    cat ${input_file} | cowpy -c ${character} > cowpy-${input_file}
    """
}
```

Comme vous pouvez le voir, ce processus est actuellement conçu pour prendre un fichier d'entrée (contenant le texte à afficher) et une valeur spécifiant le personnage qui doit être dessiné en art ASCII, généralement fournie au niveau du workflow par un paramètre en ligne de commande.

### 3.2. Passer un champ de la meta map comme entrée

Lorsque nous avons utilisé l'outil `cowpy` dans le cours Hello Nextflow, nous avons utilisé un paramètre en ligne de commande pour déterminer quel personnage utiliser pour dessiner l'image finale.
Cela avait du sens, car nous ne générions qu'une seule image par exécution du pipeline.

Cependant, dans ce tutoriel, nous voulons générer une image appropriée pour chaque sujet que nous traitons, donc utiliser un paramètre en ligne de commande serait trop limitant.

Bonne nouvelle : nous avons une colonne `character` dans notre feuille de données et donc dans notre meta map.
Utilisons-la pour définir le personnage que le processus devrait utiliser pour chaque entrée.

À cette fin, nous devrons faire trois choses :

1. Donner un nom au canal de sortie du processus précédent afin de pouvoir opérer dessus plus facilement.
2. Déterminer comment accéder aux informations d'intérêt
3. Ajouter un appel au deuxième processus et lui fournir les informations de manière appropriée.

Commençons.

#### 3.2.1. Nommer le canal de sortie précédent

Nous avons appliqué les manipulations précédentes directement sur le canal de sortie du premier processus, `IDENTIFY_LANGUAGE.out`.
Afin d'alimenter le contenu de ce canal vers le processus suivant (et de le faire d'une manière claire et facile à lire), nous voulons lui donner son propre nom, `ch_languages`.

Nous pouvons le faire en utilisant l'opérateur [`set`](https://www.nextflow.io/docs/latest/reference/operator.html#set).

Dans le workflow principal, remplacez l'opérateur `.view()` par `.set { ch_languages }`, et ajoutez une ligne pour vérifier que nous pouvons faire référence au canal par son nom.

=== "Après"

    ```groovy title="main.nf" linenums="14" hl_lines="19 21 22"
        // Exécuter langid pour identifier la langue de chaque salutation
        IDENTIFY_LANGUAGE(ch_datasheet)
        IDENTIFY_LANGUAGE.out
            .map { meta, file, lang_id ->
                [meta + [lang: lang_id], file]
            }
            .map { meta, file ->

                def lang_group = "unknown"
                if (meta.lang.equals("de") || meta.lang.equals('en')) {
                    lang_group = "germanic"
                }
                else if (meta.lang in ["fr", "es", "it"]) {
                    lang_group = "romance"
                }

                [meta + [lang_group: lang_group], file]
            }
            .set { ch_languages }

        // Temporaire : inspecter ch_languages
        ch_languages.view()
    ```

=== "Avant"

    ```groovy title="main.nf" linenums="14" hl_lines="19"
        // Exécuter langid pour identifier la langue de chaque salutation
        IDENTIFY_LANGUAGE(ch_datasheet)
        IDENTIFY_LANGUAGE.out
            .map { meta, file, lang_id ->
                [meta + [lang: lang_id], file]
            }
            .map { meta, file ->

                def lang_group = "unknown"
                if (meta.lang.equals("de") || meta.lang.equals('en')) {
                    lang_group = "germanic"
                }
                else if (meta.lang in ["fr", "es", "it"]) {
                    lang_group = "romance"
                }

                [meta + [lang_group: lang_group], file]
            }
            .view()
    ```

Exécutons ceci :

```bash
nextflow run main.nf
```

??? success "Sortie de la commande"

    ```console
     N E X T F L O W   ~  version 25.10.2

    Launching `./main.nf` [friendly_austin] DSL2 - revision: 3dbe460fd6

    [36/cca6a7] IDENTIFY_LANGUAGE (7) | 7 of 7 ✔
    [[id:sampleB, character:tux, lang:de, lang_group:germanic], /workspaces/training/side-quests/metadata/work/e2/6db2402d83cf72081bcd2d11784714/guten_tag.txt]
    [[id:sampleA, character:squirrel, lang:fr, lang_group:romance], /workspaces/training/side-quests/metadata/work/6c/114c818317d169457d6e7336d5d55b/bonjour.txt]
    [[id:sampleC, character:sheep, lang:de, lang_group:germanic], /workspaces/training/side-quests/metadata/work/55/68c69c5efb527f3604ddb3daab8057/hallo.txt]
    [[id:sampleD, character:turkey, lang:en, lang_group:germanic], /workspaces/training/side-quests/metadata/work/2a/4752055ccb5d1370b0ef9da41d3993/hello.txt]
    [[id:sampleE, character:stegosaurus, lang:es, lang_group:romance], /workspaces/training/side-quests/metadata/work/f4/fcd3186dc666d5d239ffa6c37d125d/hola.txt]
    [[id:sampleF, character:moose, lang:fr, lang_group:romance], /workspaces/training/side-quests/metadata/work/c3/3b2627f733f278a7088332a5806108/salut.txt]
    [[id:sampleG, character:turtle, lang:it, lang_group:romance], /workspaces/training/side-quests/metadata/work/36/cca6a7dbfa26ac24f9329787a32e9d/ciao.txt]
    ```

Cela confirme que nous pouvons maintenant faire référence au canal par son nom.

#### 3.2.2. Accéder aux métadonnées du fichier et du personnage

En examinant le code du module, nous savons que le processus `COWPY` s'attend à recevoir un fichier texte et une valeur `character`.
Pour écrire l'appel au processus `COWPY`, nous avons juste besoin de savoir comment extraire l'objet fichier et les métadonnées correspondants de chaque élément du canal.

Comme c'est souvent le cas, la façon la plus simple de le faire est d'utiliser une opération `map`.

Notre canal contient des tuples structurés comme `[meta, file]`, donc nous pouvons accéder directement à l'objet `file`, et nous pouvons accéder à la valeur `character` stockée dans la meta map en y faisant référence comme `meta.character`.

Dans le workflow principal, effectuez les modifications de code suivantes :

=== "Après"

    ```groovy title="main.nf" linenums="34"
        // Temporaire : accéder au fichier et au personnage
        ch_languages.map { meta, file -> file }.view { file -> "File: " + file }
        ch_languages.map { meta, file -> meta.character }.view { character -> "Character: " + character }
    ```

=== "Avant"

    ```groovy title="main.nf" linenums="34"
        // Temporaire : inspecter ch_languages
        ch_languages.view()
    ```

Notez que nous utilisons des closures (comme `{ file -> "File: " + file }`) pour rendre la sortie des opérations `.view` plus lisible.

Exécutons ceci :

```bash
nextflow run main.nf -resume
```

??? success "Sortie de la commande"

    ```console
     N E X T F L O W   ~  version 25.10.2

    Launching `./main.nf` [cheesy_cantor] DSL2 - revision: 15af9c1ec7

    [43/05df08] IDENTIFY_LANGUAGE (7) [100%] 7 of 7, cached: 7 ✔
    Character: squirrel
    File: /workspaces/training/side-quests/metadata/work/8d/4b9498bbccb7a74f04e41877cdc3e5/bonjour.txt
    File: /workspaces/training/side-quests/metadata/work/d3/604274985406e40d79021dea658e60/guten_tag.txt
    Character: tux
    Character: turkey
    File: /workspaces/training/side-quests/metadata/work/d4/fafcc9415b61d2b0fea872e6a05e8a/hello.txt
    File: /workspaces/training/side-quests/metadata/work/02/468ac9efb27f636715e8144b37e9a7/hallo.txt
    Character: sheep
    Character: moose
    Character: stegosaurus
    File: /workspaces/training/side-quests/metadata/work/d4/61a7e1188b4f2742bc72004e226eca/salut.txt
    File: /workspaces/training/side-quests/metadata/work/ae/68364be238c11149c588bf6fc858b1/hola.txt
    File: /workspaces/training/side-quests/metadata/work/43/05df081af5d879ab52e5828fa0357e/ciao.txt
    Character: turtle
    ```

_Les chemins de fichiers et les valeurs de personnage peuvent apparaître dans un ordre différent dans votre sortie._

Cela confirme que nous sommes capables d'accéder au fichier et au personnage pour chaque élément du canal.

#### 3.2.3. Appeler le processus `COWPY`

Maintenant, assemblons tout cela et appelons réellement le processus `COWPY` sur le canal `ch_languages`.

Dans le workflow principal, effectuez les modifications de code suivantes :

=== "Après"

    ```groovy title="main.nf" linenums="34"
        // Exécuter cowpy pour générer de l'art ASCII
        COWPY(
            ch_languages.map { meta, file -> file },
            ch_languages.map { meta, file -> meta.character }
        )
    ```

=== "Avant"

    ```groovy title="main.nf" linenums="34"
        // Temporaire : accéder au fichier et au personnage
        ch_languages.map { meta, file -> [file, meta.character] }
            .view()
    ```

Vous voyez que nous copions simplement les deux opérations map (sans les instructions `.view()`) comme entrées de l'appel au processus.
Assurez-vous juste de ne pas oublier la virgule entre elles !

C'est un peu maladroit, mais nous verrons comment améliorer cela dans la section suivante.

Exécutons ceci :

```bash
nextflow run main.nf -resume
```

??? success "Sortie de la commande"

    ```console
     N E X T F L O W   ~  version 25.10.2

    Launching `main.nf` [suspicious_crick] DSL2 - revision: 25541014c5

    executor >  local (7)
    [43/05df08] IDENTIFY_LANGUAGE (7) [100%] 7 of 7, cached: 7 ✔
    [e7/317c18] COWPY (6)             [100%] 7 of 7 ✔
    ```

Si vous regardez dans le répertoire des résultats, vous devriez voir les fichiers individuels contenant l'art ASCII de chaque salutation prononcée par le personnage correspondant.

??? abstract "Contenu du répertoire et exemple de fichier"

    ```console
    results/
    ├── cowpy-bonjour.txt
    ├── cowpy-ciao.txt
    ├── cowpy-guten_tag.txt
    ├── cowpy-hallo.txt
    ├── cowpy-hello.txt
    ├── cowpy-hola.txt
    └── cowpy-salut.txt
    ```

    ```text title="results/cowpy-bonjour.txt"
     _________________
    / Bonjour         \
    \ Salut, à demain /
    -----------------
      \
        \
                      _ _
          | \__/|  .~    ~.
          /oo `./      .'
          {o__,   \    {
            / .  . )    \
            `-` '-' \    }
          .(   _(   )_.'
          '---.~_ _ _|
    ```

Cela montre que nous avons pu utiliser les informations de la meta map pour paramétrer la commande dans la deuxième étape du pipeline.

Cependant, comme mentionné ci-dessus, une partie du code impliqué était un peu maladroite, car nous devions décompresser les métadonnées tout en étant encore dans le contexte du corps du workflow.
Cette approche fonctionne bien pour utiliser un petit nombre de champs de la meta map, mais ne passerait pas à l'échelle si nous voulions en utiliser beaucoup plus.

Il existe un autre opérateur appelé `multiMap()` qui nous permet de rationaliser un peu cela, mais même ainsi ce n'est pas idéal.

??? info "(Optionnel) Version alternative avec `multiMap()`"

    Au cas où vous vous poseriez la question, nous ne pouvions pas simplement écrire une seule opération `map()` qui produit à la fois le `file` et le `character`, car cela les retournerait sous forme de tuple.
    Nous avons dû écrire deux opérations `map()` séparées afin d'alimenter les éléments `file` et `character` au processus séparément.

    Techniquement, il existe une autre façon de le faire via une seule opération de mapping, en utilisant l'opérateur `multiMap()`, qui est capable d'émettre plusieurs canaux.
    Par exemple, vous pourriez remplacer l'appel à `COWPY` ci-dessus par le code suivant :

    === "Après"

        ```groovy title="main.nf" linenums="34"
            // Exécuter cowpy pour générer de l'art ASCII
            COWPY(
                ch_languages.multiMap { meta, file ->
                    file: file
                    character: meta.character
                }
            )
        ```

    === "Avant"

        ```groovy title="main.nf" linenums="34"
            // Exécuter cowpy pour générer de l'art ASCII
            COWPY(
                ch_languages.map { meta, file -> file },
                ch_languages.map { meta, file -> meta.character }
            )
        ```

    Cela produit exactement le même résultat.

Dans les deux cas, il est maladroit de devoir effectuer une partie du décompactage au niveau du workflow.

Il serait préférable de pouvoir passer l'intégralité de la meta map au processus et de sélectionner ce dont nous avons besoin une fois là-bas.

### 3.3. Passer et utiliser l'intégralité de la meta map

L'intérêt de la meta map est après tout de transmettre toutes les métadonnées ensemble comme un ensemble.
La seule raison pour laquelle nous ne pouvions pas le faire ci-dessus est que le processus n'est pas configuré pour accepter une meta map.
Mais puisque nous contrôlons le code du processus, nous pouvons changer cela.

Modifions le processus `COWPY` pour accepter la structure tuple `[meta, file]` que nous avons utilisée dans le premier processus afin de rationaliser le workflow.

À cette fin, nous devrons faire trois choses :

1. Modifier les définitions d'entrée du module `COWPY`
2. Mettre à jour la commande du processus pour utiliser la meta map
3. Mettre à jour l'appel au processus dans le corps du workflow

Prêt·e ? C'est parti !

#### 3.3.1. Modifier l'entrée du module `COWPY`

Effectuez les modifications suivantes dans le fichier module `cowpy.nf` :

=== "Après"

    ```groovy title="cowpy.nf" linenums="10" hl_lines="2"
    input:
    tuple val(meta), path(input_file)
    ```

=== "Avant"

    ```groovy title="cowpy.nf" linenums="10" hl_lines="2-3"
    input:
    path(input_file)
    val character
    ```

Cela nous permet d'utiliser la structure tuple `[meta, file]` que nous avons vue précédemment dans le tutoriel.

Notez que nous n'avons pas mis à jour la définition des sorties du processus pour inclure la meta map, afin de garder le tutoriel concis, mais n'hésitez pas à le faire vous-même comme exercice en suivant le modèle du processus `IDENTIFY_LANGUAGE`.

#### 3.3.2. Mettre à jour la commande pour utiliser le champ de la meta map

L'intégralité de la meta map est maintenant disponible à l'intérieur du processus, nous pouvons donc faire référence aux informations qu'elle contient directement depuis le bloc de commande.

Effectuez les modifications suivantes dans le fichier module `cowpy.nf` :

=== "Après"

    ```groovy title="cowpy.nf" linenums="16" hl_lines="3"
    script:
    """
    cat ${input_file} | cowpy -c ${meta.character} > cowpy-${input_file}
    """
    ```

=== "Avant"

    ```groovy title="cowpy.nf" linenums="16" hl_lines="3"
    script:
    """
    cat ${input_file} | cowpy -c ${character} > cowpy-${input_file}
    """
    ```

Nous avons remplacé la référence à la valeur `character` précédemment passée comme entrée autonome par la valeur contenue dans la meta map, à laquelle nous faisons référence en utilisant `meta.character`.

Mettons maintenant à jour l'appel au processus en conséquence.

#### 3.3.3. Mettre à jour l'appel au processus et l'exécuter

Le processus s'attend maintenant à ce que son entrée utilise la structure tuple `[meta, file]`, ce qui correspond à ce que le processus précédent produit, nous pouvons donc simplement passer l'intégralité du canal `ch_languages` au processus `COWPY`.

Effectuez les modifications suivantes dans le workflow principal :

=== "Après"

    ```groovy title="main.nf" linenums="34" hl_lines="2"
    // Exécuter cowpy pour générer de l'art ASCII
    COWPY(ch_languages)
    ```

=== "Avant"

    ```groovy title="main.nf" linenums="34" hl_lines="3-4"
    // Exécuter cowpy pour générer de l'art ASCII
    COWPY(
        ch_languages.map { meta, file -> file },
        ch_languages.map { meta, file -> meta.character }
    )
    ```

Cela simplifie considérablement l'appel !

Supprimons les résultats de l'exécution précédente et exécutons-le :

```bash
rm -r results
nextflow run main.nf
```

??? success "Sortie de la commande"

    ```console
     N E X T F L O W   ~  version 25.10.2

    Launching `main.nf` [wise_sammet] DSL2 - revision: 99797b1e92

    executor >  local (14)
    [5d/dffd4e] process > IDENTIFY_LANGUAGE (7) [100%] 7 of 7 ✔
    [25/9243df] process > COWPY (7)             [100%] 7 of 7 ✔
    ```

Si vous regardez dans le répertoire des résultats, vous devriez voir les mêmes sorties qu'auparavant, _c'est-à-dire_ des fichiers individuels contenant l'art ASCII de chaque salutation prononcée par le personnage correspondant.

??? abstract "Contenu du répertoire"

    ```console
    ./results/
    ├── cowpy-bonjour.txt
    ├── cowpy-ciao.txt
    ├── cowpy-guten_tag.txt
    ├── cowpy-hallo.txt
    ├── cowpy-hello.txt
    ├── cowpy-hola.txt
    └── cowpy-salut.txt
    ```

Cela produit donc les mêmes résultats qu'auparavant avec un code plus simple.

Bien sûr, cela suppose que vous êtes en mesure de modifier le code du processus.
Dans certains cas, vous devrez peut-être vous appuyer sur des processus existants que vous n'êtes pas libre de modifier, ce qui limite vos options.
La bonne nouvelle, si vous prévoyez d'utiliser des modules du projet [nf-core](https://nf-co.re/), est que les modules nf-core sont tous configurés pour utiliser la structure tuple `[meta, file]` comme standard.

### 3.4. Résoudre les problèmes d'entrées requises manquantes

La valeur `character` est requise pour que le processus `COWPY` s'exécute avec succès.
Si nous ne définissons pas de valeur par défaut dans un fichier de configuration, nous DEVONS fournir une valeur dans la feuille de données.

**Que se passe-t-il si nous ne le faisons pas ?**
Cela dépend de ce que contient la feuille de données d'entrée et de la version du workflow que nous exécutons.

#### 3.4.1. La colonne character existe mais est vide

Supposons que nous supprimions la valeur du personnage pour l'une des entrées de notre feuille de données pour simuler une erreur de collecte de données :

```csv title="datasheet.csv" linenums="1" hl_lines="2"
id,character,recording
sampleA,,/workspaces/training/side-quests/metadata/data/bonjour.txt
sampleB,tux,/workspaces/training/side-quests/metadata/data/guten_tag.txt
sampleC,sheep,/workspaces/training/side-quests/metadata/data/hallo.txt
sampleD,turkey,/workspaces/training/side-quests/metadata/data/hello.txt
sampleE,stegosaurus,/workspaces/training/side-quests/metadata/data/hola.txt
sampleF,moose,/workspaces/training/side-quests/metadata/data/salut.txt
sampleG,turtle,/workspaces/training/side-quests/metadata/data/ciao.txt
```

Pour l'une ou l'autre version du workflow que nous avons utilisée ci-dessus, la clé `character` sera créée pour toutes les entrées lors de la lecture de la feuille de données, mais pour `sampleA` la valeur sera une chaîne vide.

Cela provoquera une erreur.

??? failure "Sortie de la commande"

    ```console hl_lines="8 11 16 28"
     N E X T F L O W   ~  version 25.10.2

    Launching `main.nf` [marvelous_hirsch] DSL2 - revision: 0dfeee3cc1

    executor >  local (9)
    [c1/c5dd4f] process > IDENTIFY_LANGUAGE (7) [ 85%] 6 of 7
    [d3/b7c415] process > COWPY (2)             [  0%] 0 of 6
    ERROR ~ Error executing process > 'COWPY (1)'

    Caused by:
      Process `COWPY (1)` terminated with an error exit status (2)


    Command executed:

      cat bonjour.txt | cowpy -c  > cowpy-bonjour.txt

    Command exit status:
      2

    Command output:
      (empty)

    Command error:
      usage: cowpy [-h] [-l] [-L] [-t] [-u] [-e EYES] [-c COWACTER] [-E] [-r] [-x]
                  [-C]
                  [msg ...]
      cowpy: error: argument -c/--cowacter: expected one argument

    Work dir:
      /workspaces/training/side-quests/metadata/work/ca/9d49796612a54dec5ed466063c809b

    Container:
      community.wave.seqera.io/library/cowpy:1.1.5--3db457ae1977a273

    Tip: you can try to figure out what's wrong by changing to the process work dir and showing the script file named `.command.sh`

    -- Check '.nextflow.log' file for details
    ```

Lorsque Nextflow exécute la ligne de commande `cowpy` pour cet échantillon, `${meta.character}` est remplacé par une chaîne vide dans la commande `cowpy`, donc l'outil `cowpy` génère une erreur indiquant qu'aucune valeur n'a été fournie pour l'argument `-c`.

#### 3.4.2. La colonne character n'existe pas dans la feuille de données

Supposons maintenant que nous supprimions entièrement la colonne `character` de notre feuille de données :

```csv title="datasheet.csv" linenums="1"
id,recording
sampleA,/workspaces/training/side-quests/metadata/data/bonjour.txt
sampleB,/workspaces/training/side-quests/metadata/data/guten_tag.txt
sampleC,/workspaces/training/side-quests/metadata/data/hallo.txt
sampleD,/workspaces/training/side-quests/metadata/data/hello.txt
sampleE,/workspaces/training/side-quests/metadata/data/hola.txt
sampleF,/workspaces/training/side-quests/metadata/data/salut.txt
sampleG,/workspaces/training/side-quests/metadata/data/ciao.txt
```

Dans ce cas, la clé `character` ne sera pas créée du tout lors de la lecture de la feuille de données.

##### 3.4.2.1. Valeur accédée au niveau du workflow

Si nous utilisons la version du code que nous avons écrite dans la section 3.2, Nextflow tentera d'accéder à la clé `character` dans la meta map AVANT d'appeler le processus `COWPY`.

Il ne trouvera aucun élément correspondant à l'instruction, donc il n'exécutera pas `COWPY` du tout.

??? success "Sortie de la commande"

    ```console hl_lines="7"
     N E X T F L O W   ~  version 25.10.2

    Launching `main.nf` [desperate_montalcini] DSL2 - revision: 0dfeee3cc1

    executor >  local (7)
    [1a/df2544] process > IDENTIFY_LANGUAGE (7) [100%] 7 of 7 ✔
    [-        ] process > COWPY                 -
    ```

Du point de vue de Nextflow, ce workflow s'est exécuté avec succès !
Cependant, aucune des sorties souhaitées ne sera produite.

##### 3.4.2.2. Valeur accédée au niveau du processus

Si nous utilisons la version de la section 3.3, Nextflow passera l'intégralité de la meta map au processus `COWPY` et tentera d'exécuter la commande.

Cela provoquera une erreur, mais différente de celle du premier cas.

??? failure "Sortie de la commande"

    ```console hl_lines="8 11 16"
     N E X T F L O W   ~  version 25.10.2

    Launching `main.nf` [jovial_bohr] DSL2 - revision: eaaf375827

    executor >  local (9)
    [0d/ada9db] process > IDENTIFY_LANGUAGE (5) [ 85%] 6 of 7
    [06/28065f] process > COWPY (2)             [  0%] 0 of 6
    ERROR ~ Error executing process > 'COWPY (2)'

    Caused by:
      Process `COWPY (2)` terminated with an error exit status (1)


    Command executed:

      cat guten_tag.txt | cowpy -c null > cowpy-guten_tag.txt

    Command exit status:
      1

    Command output:
      (empty)

    Command error:
      Traceback (most recent call last):
        File "/opt/conda/bin/cowpy", line 10, in <module>
          sys.exit(main())
                  ~~~~^^
        File "/opt/conda/lib/python3.13/site-packages/cowpy/cow.py", line 1215, in main
          print(cow(eyes=args.eyes,
                ~~~^^^^^^^^^^^^^^^^
                tongue=args.tongue,
                ^^^^^^^^^^^^^^^^^^^
                thoughts=args.thoughts
                ^^^^^^^^^^^^^^^^^^^^^^
                    ).milk(msg)
                    ^
      TypeError: 'str' object is not callable

    Work dir:
      /workspaces/training/side-quests/metadata/work/06/28065f7d9fd7d22bba084aa941b6d6

    Container:
      community.wave.seqera.io/library/cowpy:1.1.5--3db457ae1977a273

    Tip: you can replicate the issue by changing to the process work dir and entering the command `bash .command.run`

    -- Check '.nextflow.log' file for details
    ```

Cela se produit parce que `meta.character` n'existe pas, donc notre tentative d'y accéder retourne `null`. Par conséquent, Nextflow insère littéralement `null` dans la ligne de commande, ce qui n'est bien sûr pas reconnu par l'outil `cowpy`.

#### 3.4.3. Solutions

En dehors de la fourniture d'une valeur par défaut dans la configuration du workflow, il y a deux choses que nous pouvons faire pour gérer cela de manière plus robuste :

1. Implémenter une validation des entrées dans votre workflow pour s'assurer que la feuille de données contient toutes les informations requises. Vous pouvez trouver une [introduction à la validation des entrées](../hello_nf-core/05_input_validation.md) dans le cours de formation Hello nf-core. <!-- TODO (future) pending a proper Validation side quest -->

2. Si vous voulez vous assurer que toute personne utilisant votre module de processus peut immédiatement identifier les entrées requises, vous pouvez également faire de la propriété de métadonnées requise une entrée explicite.

Voici un exemple de comment cela fonctionnerait.

Premièrement, au niveau du processus, mettez à jour la définition des entrées comme suit :

=== "Après"

    ```groovy title="cowpy.nf" linenums="12" hl_lines="2"
        input:
        tuple val(meta), val(character), path(input_file)
    ```

=== "Avant"

    ```groovy title="cowpy.nf" linenums="12" hl_lines="2"
        input:
        tuple val(meta), path(input_file)
    ```

Ensuite, au niveau du workflow, utilisez une opération de mapping pour extraire la propriété `character` des métadonnées et en faire un composant explicite du tuple d'entrée :

=== "Après"

    ```groovy title="main.nf" linenums="37" hl_lines="1"
        COWPY(ch_languages.map{meta, file -> [meta, meta.character, file]})
    ```

=== "Avant"

    ```groovy title="main.nf" linenums="37" hl_lines="1"
        COWPY(ch_languages)
    ```

Cette approche a l'avantage de montrer explicitement que `character` est requis, et rend le processus plus facile à redéployer dans d'autres contextes.

Cela met en évidence un principe de conception important :

**Utilisez la meta map pour les informations optionnelles et descriptives, mais extrayez les valeurs requises comme entrées explicites.**

La meta map est excellente pour maintenir des structures de canaux propres et éviter des structures de canaux arbitraires, mais pour les éléments obligatoires directement référencés dans un processus, les extraire comme entrées explicites crée un code plus robuste et maintenable.

### À retenir

Dans cette section, vous avez appris à utiliser les métadonnées pour personnaliser l'exécution d'un processus, en y accédant soit au niveau du workflow, soit au niveau du processus.

---

## Exercice supplémentaire

Si vous souhaitez vous entraîner à utiliser les informations de la meta map depuis l'intérieur d'un processus, essayez d'utiliser d'autres informations de la meta map telles que `lang` et `lang_group` pour personnaliser la façon dont les sorties sont nommées et/ou organisées.

Par exemple, essayez de modifier le code pour produire ce résultat :

```console title="Results directory contents"
results/
├── germanic
│   ├── de-guten_tag.txt
│   ├── de-hallo.txt
│   └── en-hello.txt
└── romance
    ├── es-hola.txt
    ├── fr-bonjour.txt
    ├── fr-salut.txt
    └── it-ciao.txt
```

<!-- TODO (future) Provide worked out solution -->
<!-- the renaming should use the meta inside the process -->
<!-- the output org should use the meta in the workflow outputs -->

---

## Résumé

Dans cette quête secondaire, vous avez exploré comment travailler efficacement avec les métadonnées dans les workflows Nextflow.

Ce schéma consistant à maintenir les métadonnées explicites et attachées aux données est une bonne pratique fondamentale dans Nextflow, offrant plusieurs avantages par rapport au codage en dur des informations sur les fichiers :

- Les métadonnées des fichiers restent associées aux fichiers tout au long du workflow
- Le comportement des processus peut être personnalisé par fichier
- L'organisation des sorties peut refléter les métadonnées des fichiers
- Les informations sur les fichiers peuvent être enrichies pendant l'exécution du pipeline

Appliquer ce schéma dans votre propre travail vous permettra de construire des workflows bioinformatiques robustes et maintenables.

### Schémas clés

1.  **Lecture et structuration des métadonnées :** Lire des fichiers CSV et créer des meta maps organisées qui restent associées à vos fichiers de données.

    ```groovy
    channel.fromPath('datasheet.csv')
      .splitCsv(header: true)
      .map { row ->
          [ [id:row.id, character:row.character], row.recording ]
      }
    ```

2.  **Enrichissement des métadonnées pendant le workflow :** Ajouter de nouvelles informations à vos métadonnées au fur et à mesure de la progression de votre pipeline en ajoutant des sorties de processus et en dérivant des valeurs via une logique conditionnelle.

    - Ajout de nouvelles clés basées sur la sortie d'un processus

    ```groovy
    .map { meta, file, lang ->
      [ meta + [lang:lang], file ]
    }
    ```

    - Ajout de nouvelles clés en utilisant une clause conditionnelle

    ```groovy
    .map{ meta, file ->
        if ( meta.lang.equals("de") || meta.lang.equals('en') ){
            lang_group = "germanic"
        } else if ( meta.lang in ["fr", "es", "it"] ) {
            lang_group = "romance"
        } else {
            lang_group = "unknown"
        }
    }
    ```

3.  **Personnalisation du comportement des processus :** Utiliser les métadonnées à l'intérieur du processus.

    ```groovy
    cat $input_file | cowpy -c ${meta.character} > cowpy-${input_file}
    ```

### Ressources supplémentaires

- [map](https://www.nextflow.io/docs/latest/operator.html#map)
- [stdout](https://www.nextflow.io/docs/latest/process.html#outputs)

---

## Et ensuite ?

Retournez au [menu des Quêtes secondaires](../) ou cliquez sur le bouton en bas à droite de la page pour passer au sujet suivant dans la liste.
