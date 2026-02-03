# Métadonnées et meta maps

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Traduction assistée par IA - [en savoir plus et suggérer des améliorations](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

Dans toute analyse scientifique, nous travaillons rarement uniquement avec les fichiers de données brutes.
Chaque fichier s'accompagne de ses propres informations supplémentaires : ce qu'il est, d'où il vient et ce qui le rend spécial.
Ces informations supplémentaires sont ce que nous appelons les métadonnées.

Les métadonnées sont des données décrivant d'autres données.
Les métadonnées permettent de suivre les détails importants concernant les fichiers et les conditions expérimentales, et aident à adapter les analyses aux caractéristiques uniques de chaque ensemble de données.

Imaginez cela comme un catalogue de bibliothèque : alors que les livres contiennent le contenu réel (données brutes), les fiches du catalogue fournissent des informations essentielles sur chaque livre—quand il a été publié, qui l'a écrit, où le trouver (métadonnées).
Dans les pipelines Nextflow, les métadonnées peuvent être utilisées pour :

- Suivre les informations spécifiques aux fichiers tout au long du workflow
- Configurer les processus en fonction des caractéristiques des fichiers
- Regrouper les fichiers liés pour une analyse conjointe

### Objectifs d'apprentissage

Dans cette quête secondaire, nous explorerons comment gérer les métadonnées dans les workflows.
En partant d'une simple feuille de données (souvent appelée samplesheet en bioinformatique) contenant des informations de base sur les fichiers, vous apprendrez à :

- Lire et analyser les métadonnées de fichiers à partir de fichiers CSV
- Créer et manipuler des maps de métadonnées
- Ajouter de nouveaux champs de métadonnées pendant l'exécution du workflow
- Utiliser les métadonnées pour personnaliser le comportement des processus

Ces compétences vous aideront à construire des pipelines plus robustes et flexibles qui peuvent gérer des relations de fichiers et des exigences de traitement complexes.

### Prérequis

Avant d'entreprendre cette quête secondaire, vous devriez :

- Avoir complété le tutoriel [Hello Nextflow](../hello_nextflow/README.md) ou un cours équivalent pour débutants.
- Être à l'aise avec l'utilisation des concepts et mécanismes de base de Nextflow (processus, canaux, opérateurs)

---

## 0. Pour commencer

#### Ouvrir le codespace de formation

Si vous ne l'avez pas encore fait, assurez-vous d'ouvrir l'environnement de formation comme décrit dans la [Configuration de l'environnement](../envsetup/index.md).

[![Open in GitHub Codespaces](https://github.com/codespaces/badge.svg)](https://codespaces.new/nextflow-io/training?quickstart=1&ref=master)

#### Se déplacer dans le répertoire du projet

Déplaçons-nous dans le répertoire où se trouvent les fichiers pour ce tutoriel.

```bash
cd side-quests/metadata
```

Vous pouvez configurer VSCode pour se concentrer sur ce répertoire :

```bash
code .
```

#### Examiner les ressources

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

Le workflow dans le fichier `main.nf` est une ébauche que vous développerez progressivement en un workflow pleinement fonctionnel.

La feuille de données répertorie les chemins vers les fichiers de données et certaines métadonnées associées, organisées en 3 colonnes :

- `id` : auto-explicatif, un identifiant attribué au fichier
- `character` : un nom de personnage, que nous utiliserons plus tard pour dessiner différentes créatures
- `data` : chemins vers les fichiers `.txt` qui contiennent des salutations dans différentes langues

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

Nous vous fournirons également un outil d'analyse de langue conteneurisé appelé `langid`.

#### Examiner l'assignation

Votre défi est d'écrire un workflow Nextflow qui va :

1. **Identifier** automatiquement la langue dans chaque fichier
2. **Regrouper** les fichiers par famille de langues (langues germaniques vs langues romanes)
3. **Personnaliser** le traitement de chaque fichier en fonction de sa langue et de ses métadonnées
4. **Organiser** les sorties par groupe de langues

Cela représente un modèle de workflow typique où les métadonnées spécifiques aux fichiers guident les décisions de traitement ; exactement le type de problème que les meta maps résolvent élégamment.

#### Liste de vérification de préparation

Pensez-vous être prêt·e à vous lancer ?

- [ ] Je comprends l'objectif de ce cours et ses prérequis
- [ ] Mon codespace est opérationnel
- [ ] J'ai défini mon répertoire de travail de manière appropriée
- [ ] Je comprends l'assignation

Si vous pouvez cocher toutes les cases, vous êtes prêt·e à commencer.

---

## 1. Charger les métadonnées depuis une feuille de données

Ouvrez le fichier de workflow `main.nf` pour examiner l'ébauche de workflow que nous vous donnons comme point de départ.

```groovy title="main.nf" linenums="1"
#!/usr/bin/env nextflow

workflow  {

    ch_datasheet = channel.fromPath("./data/datasheet.csv")

}
```

Vous pouvez voir que nous avons configuré une factory de canal de base pour charger l'exemple de feuille de données en tant que fichier, mais cela ne lira pas encore le contenu du fichier.
Commençons par ajouter cela.

### 1.1. Lire le contenu avec `splitCsv`

Nous devons choisir un opérateur qui analysera le contenu du fichier de manière appropriée avec un effort minimal de notre part.
Puisque notre feuille de données est au format CSV, c'est un travail pour l'opérateur [`splitCsv`](https://www.nextflow.io/docs/latest/reference/operator.html#splitcsv), qui charge chaque ligne du fichier comme un élément dans le canal.

Effectuez les modifications suivantes pour ajouter une opération `splitCsv()` au code de construction du canal, plus une opération `view()` pour vérifier que le contenu du fichier est correctement chargé dans le canal.

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

Voyons ce qui en ressort, d'accord ?
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

Chaque entrée de map correspond à une colonne dans notre feuille de données :

- `id`
- `character`
- `recording`

C'est excellent ! Cela facilite l'accès à des champs spécifiques de chaque fichier.
Par exemple, nous pourrions accéder à l'identifiant du fichier avec `id` ou au chemin du fichier txt avec `recording`.

??? info "(Optionnel) En savoir plus sur les maps"

    Dans Groovy, le langage de programmation sur lequel Nextflow est construit, une map est une structure de données clé-valeur similaire aux dictionnaires en Python, aux objets en JavaScript ou aux hashes en Ruby.

    Voici un script exécutable qui montre comment vous pouvez définir une map et accéder à son contenu en pratique :

    ```groovy title="examples/map_demo.nf"
    #!/usr/bin/env nextflow

    // Créer une map simple
    def my_map = [id:'sampleA', character:'squirrel']

    // Afficher toute la map
    println "map: ${my_map}"

    // Accéder aux valeurs individuelles en utilisant la notation par point
    println "id: ${my_map.id}"
    println "character: ${my_map.character}"
    ```

    Même s'il n'a pas de bloc `workflow` approprié, Nextflow peut exécuter cela comme s'il s'agissait d'un workflow :

    ```bash
    nextflow run examples/map_demo.nf
    ```

    Et voici ce que vous pouvez vous attendre à voir dans la sortie :

    ```console title="Sortie"
     N E X T F L O W   ~  version 25.10.2

    Launching `map_demo.nf` [cheesy_plateau] DSL2 - revision: fae5b8496e

    map: [id:sampleA, character:squirrel]
    id: sampleA
    character: squirrel
    ```

### 1.2. Sélectionner des champs spécifiques avec `map`

Disons que nous voulons accéder à la colonne `character` de la feuille de données et l'afficher.
Nous pouvons utiliser l'opérateur Nextflow `map` pour itérer sur chaque élément dans notre canal et sélectionner spécifiquement l'entrée `character` de l'objet map.

Effectuez les modifications suivantes au workflow :

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

Maintenant, exécutez à nouveau le workflow :

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

Succès ! Nous avons tiré parti de la structure de map dérivée de notre feuille de données pour accéder aux valeurs de colonnes individuelles pour chaque ligne.

Maintenant que nous avons lu avec succès la feuille de données et que nous avons accès aux données de chaque ligne, nous pouvons commencer à implémenter la logique de notre pipeline.

### 1.3. Organiser les métadonnées dans une 'meta map'

Dans l'état actuel du workflow, les fichiers d'entrée (sous la clé `recording`) et les métadonnées associées (`id`, `character`) sont tous sur le même pied, comme s'ils étaient tous dans un grand sac.
La conséquence pratique est que chaque processus qui consomme ce canal devrait être configuré avec cette structure à l'esprit :

```groovy
    input:
    tuple val(id), val(character), file(recording)
```

C'est bien tant que le nombre de colonnes dans la feuille de données ne change pas.
Cependant, si vous ajoutez ne serait-ce qu'une seule colonne à la feuille de données, la forme du canal ne correspondra plus à ce que le processus attend, et le workflow produira des erreurs.
Cela rend également le processus difficile à partager avec d'autres qui pourraient avoir des données d'entrée légèrement différentes, et vous pourriez finir par devoir coder en dur des variables dans le processus qui ne sont pas nécessaires au bloc script.

Pour éviter ce problème, nous devons trouver un moyen de maintenir la structure du canal cohérente quel que soit le nombre de colonnes que contient la feuille de données.

Nous pouvons le faire en collectant toutes les métadonnées dans un élément au sein du tuple, que nous appellerons la map de métadonnées, ou plus simplement 'meta map'.

Effectuez les modifications suivantes à l'opération `map` :

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

Nous avons restructuré nos éléments de canal en un tuple composé de deux éléments, la meta map et l'objet fichier correspondant.

Exécutons le workflow :

```bash
nextflow run main.nf
```

??? success "Sortie de la commande"

    ```console title="Afficher la meta map"
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

Maintenant, chaque élément dans le canal contient la meta map en premier et l'objet fichier correspondant en second :

```console title="Exemple de structure de sortie"
[
  [id:sampleA, character:squirrel],
  /workspaces/training/side-quests/metadata/data/bonjour.txt
]
```

En conséquence, l'ajout de plus de colonnes dans la feuille de données rendra plus de métadonnées disponibles dans la map `meta`, mais ne changera pas la forme du canal.
Cela nous permet d'écrire des processus qui consomment le canal sans avoir à coder en dur les éléments de métadonnées dans la spécification d'entrée :

```groovy title="Exemple de syntaxe"
    input:
    tuple val(meta), file(recording)
```

C'est une convention largement utilisée pour organiser les métadonnées dans les workflows Nextflow.

### À retenir

Dans cette section, vous avez appris :

- **Pourquoi les métadonnées sont importantes :** Garder les métadonnées avec vos données préserve les informations importantes sur les fichiers tout au long du workflow.
- **Comment lire des feuilles de données :** Utiliser `splitCsv` pour lire les fichiers CSV avec des informations d'en-tête et transformer les lignes en données structurées
- **Comment créer une meta map :** Séparer les métadonnées des données de fichiers en utilisant la structure de tuple `[ [id:value, ...], file ]`

---

## 2. Manipuler les métadonnées

Maintenant que nous avons chargé nos métadonnées, faisons quelque chose avec !

Nous allons utiliser un outil appelé [`langid`](https://github.com/saffsd/langid.py) pour identifier la langue contenue dans le fichier d'enregistrement de chaque créature.
L'outil est pré-entraîné sur un ensemble de langues, et étant donné un extrait de texte, il produira une prédiction de langue et un score de probabilité associé, tous deux vers `stdout`.

### 2.1. Importer le processus et examiner le code

Nous vous fournissons un module de processus pré-écrit appelé `IDENTIFY_LANGUAGE` qui encapsule l'outil `langid`, vous devez donc simplement ajouter une instruction include avant le bloc workflow.

Effectuez la modification suivante au workflow :

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

Vous pouvez ouvrir le fichier de module pour examiner son code :

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

Comme vous pouvez le voir, la définition d'entrée utilise la même structure `tuple val(meta), path(file)` que nous venons d'appliquer à notre canal d'entrée.

La définition de sortie est structurée comme un tuple avec une structure similaire à celle de l'entrée, sauf qu'elle contient également `stdout` comme troisième élément.
Ce modèle `tuple val(meta), path(file), <sortie>` maintient les métadonnées associées à la fois aux données d'entrée et aux sorties alors qu'elles circulent dans le pipeline.

Notez que nous utilisons le qualificateur de sortie [`stdout`](https://www.nextflow.io/docs/latest/process.html#outputs) de Nextflow ici parce que l'outil affiche sa sortie directement sur la console plutôt que d'écrire un fichier ; et nous utilisons `sed` dans la ligne de commande pour supprimer le score de probabilité, nettoyer la chaîne en supprimant les caractères de nouvelle ligne, et ne retourner que la prédiction de langue.

### 2.2. Ajouter un appel à `IDENTIFY_LANGUAGE`

Maintenant que le processus est disponible pour le workflow, nous pouvons ajouter un appel au processus `IDENTIFY_LANGUAGE` pour l'exécuter sur le canal de données.

Effectuez les modifications suivantes au workflow :

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

Et comme noté précédemment, nous avons également inclus le fichier d'entrée et la meta map dans la sortie, ce qui signifie que les deux restent associés aux nouvelles informations que nous venons de produire.
Cela s'avérera utile dans l'étape suivante.

!!! note

    Plus généralement, ce modèle consistant à garder la meta map associée aux résultats facilite l'association de résultats liés qui partagent les mêmes identifiants.

    Comme vous l'aurez déjà appris, vous ne pouvez pas compter sur l'ordre des éléments dans les canaux pour faire correspondre les résultats entre eux.
    Au lieu de cela, vous devez utiliser des clés pour associer correctement les données, et les meta maps fournissent une structure idéale à cette fin.

    Nous explorons ce cas d'usage en détail dans la quête secondaire [Splitting & Grouping](./splitting_and_grouping.md).

### 2.3. Augmenter les métadonnées avec les sorties de processus

Étant donné que les résultats que nous venons de produire sont en eux-mêmes une forme de métadonnées sur le contenu des fichiers, il serait utile de les ajouter à notre meta map.

Cependant, nous ne voulons pas modifier la meta map existante sur place.
D'un point de vue technique, il est _possible_ de le faire, mais c'est risqué.

Donc à la place, nous allons créer une nouvelle meta map contenant le contenu de la meta map existante plus une nouvelle paire clé-valeur `lang: lang_id` contenant les nouvelles informations, en utilisant l'opérateur `+` (une fonctionnalité Groovy).
Et nous combinerons cela avec une opération [`map`](https://www.nextflow.io/docs/latest/operator.html#map) pour remplacer l'ancienne map par la nouvelle.

Voici les modifications que vous devez apporter au workflow :

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

Si vous n'êtes pas encore familier avec l'opérateur `+`, ou si cela semble déroutant, prenez quelques minutes pour parcourir l'explication détaillée ci-dessous.

??? info "Création de la nouvelle meta map en utilisant l'opérateur `+`"

    **D'abord, vous devez savoir que nous pouvons fusionner le contenu de deux maps en utilisant l'opérateur Groovy `+`.**

    Disons que nous avons les maps suivantes :

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

    Disons que vous repartez de `map1`, mais la prédiction de langue n'est pas dans sa propre map (il n'y a pas de `map2`).
    Au lieu de cela, elle est contenue dans une variable appelée `lang_id`, et vous savez que vous voulez stocker sa valeur (`'fr'`) avec la clé `lang`.

    Vous pouvez en fait faire ce qui suit :

    ```groovy
    new_map = [map1 + [lang: lang_id]]
    ```

    Ici, `[lang: new_info]` crée une nouvelle map sans nom à la volée, et `map1 + ` fusionne `map1` avec la nouvelle map sans nom, produisant le même contenu `new_map` qu'auparavant.

    Élégant, n'est-ce pas ?

    **Maintenant transposons cela dans le contexte d'une opération Nextflow `channel.map()`.**

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
    Donc nous avons effectivement transformé :

    ```groovy
    [id: 'sampleA', character: 'squirrel'], 'it'
    ```

    en :

    ```groovy
    [id: 'sampleA', character: 'squirrel', lang: 'fr']
    ```

    Espérons que vous pouvez voir que si nous changeons `map1` en `meta`, c'est essentiellement tout ce dont nous avons besoin pour ajouter la prédiction de langue à notre meta map dans notre workflow.

    Sauf pour une chose !

    Dans le cas de notre workflow, **nous devons également tenir compte de la présence de l'objet `file` dans le tuple**, qui est composé de `meta, file, lang_id`.

    Donc le code ici deviendrait :

    ```groovy
    .map { meta, file, lang_id ->
        [meta + [lang: lang_id], file]
    }
    ```

    Si vous avez du mal à comprendre pourquoi le `file` semble se déplacer dans l'opération `map`, imaginez qu'au lieu de `[meta + [lang: lang_id], file]`, cette ligne lise `[new_map, file]`.
    Cela devrait rendre plus clair que nous laissons simplement le `file` à sa place d'origine en deuxième position dans le tuple. Nous avons juste pris la valeur `new_info` et l'avons intégrée dans la map qui est en première position.

    **Et cela nous ramène à la structure de canal `tuple val(meta), path(file)` !**

Une fois que vous êtes sûr de comprendre ce que fait ce code, exécutez le workflow pour voir si cela a fonctionné :

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

Oui, ça fonctionne !
Nous avons soigneusement réorganisé la sortie du processus de `meta, file, lang_id` pour que `lang_id` soit maintenant l'une des clés dans la meta map, et les tuples du canal correspondent à nouveau au modèle `meta, file`.

### 2.4. Attribuer un groupe de langue en utilisant des conditions

Maintenant que nous avons nos prédictions de langue, utilisons les informations pour attribuer de nouveaux regroupements.

Dans nos données d'exemple, les langues utilisées par nos personnages peuvent être regroupées en langues germaniques (anglais, allemand) et langues romanes (français, espagnol, italien).
Il pourrait être utile d'avoir cette classification facilement disponible quelque part plus tard dans le pipeline, alors ajoutons cette information dans la meta map.

Et, bonne nouvelle, c'est encore un autre cas qui se prête parfaitement à l'utilisation de l'opérateur `map` !

Spécifiquement, nous allons définir une variable appelée `lang_group`, utiliser une logique conditionnelle simple pour déterminer quelle valeur attribuer au `lang_group` pour chaque donnée.

La syntaxe générale va ressembler à ceci :

```groovy
.map { meta, file ->

    // la logique conditionnelle définissant lang_group va ici

    [meta + [lang_group: lang_group], file]
}
```

Vous pouvez voir que c'est très similaire à l'opération de fusion de map à la volée que nous avons utilisée à l'étape précédente.
Nous devons juste écrire les instructions conditionnelles.

Voici la logique conditionnelle que nous voulons appliquer :

- Définir une variable appelée `lang_group` avec la valeur par défaut `'unknown'`.
- Si `lang` est soit allemand (`'de'`) soit anglais (`'en'`), changer `lang_group` en `germanic`.
- Sinon si `lang` est inclus dans une liste contenant français (`'fr'`), espagnol (`'es'`) et italien (`'it'`), changer `lang_group` en `romance`.

Essayez de l'écrire vous-même si vous savez déjà comment écrire des instructions conditionnelles dans Nextflow.

!!! tip

    Vous pouvez accéder à la valeur de `lang` dans l'opération map avec `meta.lang`.

Vous devriez finir par apporter les modifications suivantes au workflow :

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

- Nous utilisons `def lang_group = "unknown"` pour créer la variable `lang_group` avec une valeur par défaut définie sur `unknown`.
- Nous utilisons une structure `if {} else if {}` pour la logique conditionnelle, avec des tests `.equals()` alternatifs pour les deux langues germaniques, et un test d'existence dans une liste pour les trois langues romanes.
- Nous utilisons l'opération de fusion `meta + [lang_group:lang_group]` comme précédemment pour générer la meta map mise à jour.

Une fois que tout cela a du sens, exécutez à nouveau le workflow pour voir le résultat :

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

Dans cette section, vous avez appris comment :

- **Appliquer les métadonnées d'entrée aux canaux de sortie** : Copier les métadonnées de cette manière nous permet d'associer les résultats ultérieurement en fonction du contenu des métadonnées.
- **Créer des clés personnalisées** : Vous avez créé deux nouvelles clés dans votre meta map, les fusionnant avec `meta + [new_key:value]` dans la meta map existante. L'une basée sur une valeur calculée à partir d'un processus, et une basée sur une condition que vous avez définie dans l'opérateur `map`.

Celles-ci vous permettent d'associer des métadonnées nouvelles et existantes avec des fichiers au fur et à mesure que vous progressez dans votre pipeline.
Même si vous n'utilisez pas les métadonnées dans le cadre d'un processus, garder la meta map associée aux données comme ceci facilite le maintien de toutes les informations pertinentes ensemble.

---

## 3. Utiliser les informations de la meta map dans un processus

Maintenant que vous savez comment créer et mettre à jour la meta map, nous pouvons passer à la partie vraiment amusante : utiliser réellement les métadonnées dans un processus.

Plus spécifiquement, nous allons ajouter une deuxième étape à notre workflow pour dessiner chaque animal en art ASCII et le faire dire le texte enregistré dans une bulle de dialogue.
Nous allons faire cela en utilisant un outil appelé [`cowpy`](https://github.com/jeffbuttars/cowpy).

??? info "Que fait `cowpy` ?"

    `cowpy` est un outil en ligne de commande qui génère de l'art ASCII pour afficher des entrées de texte arbitraires de manière amusante.
    C'est une implémentation python du classique outil [cowsay](https://en.wikipedia.org/wiki/Cowsay) par Tony Monroe.

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

    Optionnellement, vous pouvez sélectionner un personnage (ou 'cowacter') à utiliser au lieu de la vache par défaut.

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

Nous vous fournissons un module de processus pré-écrit appelé `COWPY` qui encapsule l'outil `cowpy`, vous devez donc simplement ajouter une instruction include avant le bloc workflow.

Effectuez la modification suivante au workflow :

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

Vous pouvez ouvrir le fichier de module pour examiner son code :

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

Comme vous pouvez le voir, ce processus est actuellement conçu pour prendre un fichier d'entrée (contenant le texte à afficher) et une valeur spécifiant le personnage qui doit être dessiné en art ASCII, généralement fournie au niveau du workflow par un paramètre de ligne de commande.

### 3.2. Passer un champ de meta map comme entrée

Lorsque nous avons utilisé l'outil `cowpy` dans le cours Hello Nextflow, nous avons utilisé un paramètre de ligne de commande pour déterminer quel personnage utiliser pour dessiner l'image finale.
Cela avait du sens, car nous ne générions qu'une seule image par exécution du pipeline.

Cependant, dans ce tutoriel, nous voulons générer une image appropriée pour chaque sujet que nous traitons, donc utiliser un paramètre de ligne de commande serait trop limitant.

Bonne nouvelle : nous avons une colonne `character` dans notre feuille de données et donc, dans notre meta map.
Utilisons cela pour définir le personnage que le processus doit utiliser pour chaque entrée.

À cette fin, nous devrons faire trois choses :

1. Donner un nom au canal de sortie provenant du processus précédent afin de pouvoir l'utiliser plus commodément.
2. Déterminer comment accéder aux informations qui nous intéressent
3. Ajouter un appel au deuxième processus et fournir les informations de manière appropriée.

Commençons.

#### 3.2.1. Nommer le canal de sortie précédent

Nous avons appliqué les manipulations précédentes directement sur le canal de sortie du premier processus, `IDENTIFY_LANGUAGE.out`.
Afin de fournir le contenu de ce canal au processus suivant (et de le faire d'une manière claire et facile à lire) nous voulons lui donner son propre nom, `ch_languages`.

Nous pouvons le faire en utilisant l'opérateur [`set`](https://www.nextflow.io/docs/latest/reference/operator.html#set).

Dans le workflow principal, remplacez l'opérateur `.view()` par `.set { ch_languages }`, et ajoutez une ligne testant que nous pouvons référer au canal par nom.

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

        // Temporaire : jeter un œil dans ch_languages
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

Cela confirme que nous pouvons maintenant référer au canal par nom.

#### 3.2.2. Accéder au fichier et aux métadonnées de personnage

Nous savons en regardant le code du module que le processus `COWPY` s'attend à recevoir un fichier texte et une valeur `character`.
Pour écrire l'appel au processus `COWPY`, nous devons simplement savoir comment extraire l'objet fichier correspondant et les métadonnées de chaque élément du canal.

Comme c'est souvent le cas, la façon la plus simple de le faire est d'utiliser une opération `map`.

Notre canal contient des tuples structurés comme `[meta, file]`, donc nous pouvons accéder à l'objet `file` directement, et nous pouvons accéder à la valeur `character` stockée à l'intérieur de la meta map en y faisant référence comme `meta.character`.

Dans le workflow principal, effectuez les modifications de code suivantes :

=== "Après"

    ```groovy title="main.nf" linenums="34"
        // Temporaire : accéder au fichier et au personnage
        ch_languages.map { meta, file -> file }.view { file -> "Fichier : " + file }
        ch_languages.map { meta, file -> meta.character }.view { character -> "Personnage : " + character }
    ```

=== "Avant"

    ```groovy title="main.nf" linenums="34"
        // Temporaire : jeter un œil dans ch_languages
        ch_languages.view()
    ```

Notez que nous utilisons des closures (telles que `{ file -> "Fichier : " + file }`) pour rendre la sortie des opérations `.view` plus lisible.

Exécutons ceci :

```bash
nextflow run main.nf -resume
```

??? success "Sortie de la commande"

    ```console
     N E X T F L O W   ~  version 25.10.2

    Launching `./main.nf` [cheesy_cantor] DSL2 - revision: 15af9c1ec7

    [43/05df08] IDENTIFY_LANGUAGE (7) [100%] 7 of 7, cached: 7 ✔
    Personnage : squirrel
    Fichier : /workspaces/training/side-quests/metadata/work/8d/4b9498bbccb7a74f04e41877cdc3e5/bonjour.txt
    Fichier : /workspaces/training/side-quests/metadata/work/d3/604274985406e40d79021dea658e60/guten_tag.txt
    Personnage : tux
    Personnage : turkey
    Fichier : /workspaces/training/side-quests/metadata/work/d4/fafcc9415b61d2b0fea872e6a05e8a/hello.txt
    Fichier : /workspaces/training/side-quests/metadata/work/02/468ac9efb27f636715e8144b37e9a7/hallo.txt
    Personnage : sheep
    Personnage : moose
    Personnage : stegosaurus
    Fichier : /workspaces/training/side-quests/metadata/work/d4/61a7e1188b4f2742bc72004e226eca/salut.txt
    Fichier : /workspaces/training/side-quests/metadata/work/ae/68364be238c11149c588bf6fc858b1/hola.txt
    Fichier : /workspaces/training/side-quests/metadata/work/43/05df081af5d879ab52e5828fa0357e/ciao.txt
    Personnage : turtle
    ```

_Les chemins de fichiers et les valeurs de personnages peuvent sortir dans un ordre différent dans votre sortie._

Cela confirme que nous sommes capables d'accéder au fichier et au personnage pour chaque élément du canal.

#### 3.2.3. Appeler le processus `COWPY`

Maintenant, assemblons tout et appelons réellement le processus `COWPY` sur le canal `ch_languages`.

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

Vous voyez que nous copions simplement les deux opérations map (moins les instructions `.view()`) comme entrées de l'appel de processus.
Assurez-vous simplement de ne pas oublier la virgule entre elles !

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

Si vous regardez dans le répertoire results, vous devriez voir les fichiers individuels contenant l'art ASCII de chaque salutation prononcée par le personnage correspondant.

??? abstract "Répertoire et exemple de contenu de fichier"

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

Cela montre que nous avons pu utiliser les informations dans la meta map pour paramétrer la commande dans la deuxième étape du pipeline.

Cependant, comme noté ci-dessus, une partie du code impliqué était un peu maladroite, puisque nous avons dû déballer les métadonnées tout en étant encore dans le contexte du corps du workflow.
Cette approche fonctionne bien pour utiliser un petit nombre de champs de la meta map, mais évoluerait mal si nous voulions en utiliser beaucoup plus.

Il existe un autre opérateur appelé `multiMap()` qui nous permet de rationaliser cela un peu, mais même alors ce n'est pas idéal.

??? info "(Optionnel) Version alternative avec `multiMap()`"

    Au cas où vous vous poseriez la question, nous ne pouvions pas simplement écrire une seule opération `map()` qui sort à la fois le `file` et le `character`, parce que cela les retournerait comme un tuple.
    Nous avons dû écrire deux opérations `map()` séparées afin de fournir les éléments `file` et `character` au processus séparément.

    Techniquement, il existe une autre façon de le faire à travers une seule opération de mapping, en utilisant l'opérateur `multiMap()`, qui est capable d'émettre plusieurs canaux.
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

Dans les deux cas, il est gênant de devoir faire un peu de déballage au niveau du workflow.

Il serait préférable de pouvoir fournir la meta map entière au processus et choisir ce dont nous avons besoin une fois là-bas.

### 3.3. Passer et utiliser la meta map entière

Le but de la meta map est après tout de passer toutes les métadonnées ensemble comme un ensemble.
La seule raison pour laquelle nous ne pouvions pas le faire ci-dessus est que le processus n'est pas configuré pour accepter une meta map.
Mais puisque nous contrôlons le code du processus, nous pouvons changer cela.

Modifions le processus `COWPY` pour accepter la structure de tuple `[meta, file]` que nous avons utilisée dans le premier processus afin de pouvoir rationaliser le workflow.

À cette fin, nous devrons faire trois choses :

1. Modifier les définitions d'entrée du module de processus `COWPY`
2. Mettre à jour la commande du processus pour utiliser la meta map
3. Mettre à jour l'appel de processus dans le corps du workflow

Prêt ? Allons-y !

#### 3.3.1. Modifier l'entrée du module `COWPY`
