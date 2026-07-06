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
- Comprendre pourquoi l'interface « meta map + fichier de données » est une convention largement utilisée
- Ajouter de nouveaux champs de métadonnées pendant l'exécution du workflow
- Utiliser les métadonnées pour personnaliser le comportement des processus et organiser les sorties

Ces compétences vous aideront à construire des pipelines plus robustes et flexibles, capables de gérer des relations complexes entre fichiers et des exigences de traitement variées.

### Prérequis

Avant de vous lancer dans cette quête secondaire, vous devriez :

- Avoir suivi le tutoriel [Hello Nextflow](../../hello_nextflow/index.md) ou un cours équivalent pour débutant·es.
- Être à l'aise avec les concepts et mécanismes de base de Nextflow (processus, canaux, opérateurs)

---

## 0. Premiers pas

#### Ouvrir le codespace de formation

Si vous ne l'avez pas encore fait, assurez-vous d'ouvrir l'environnement de formation comme décrit dans la [Configuration de l'environnement](../../envsetup/index.md).

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

L'éditeur s'ouvre avec le répertoire du projet en focus.

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

Nous allons utiliser un outil appelé [`COWPY`](https://github.com/jeffbuttars/cowpy) pour générer de l'art ASCII de chaque personnage prononçant sa salutation enregistrée.

??? info "Que fait `COWPY` ?"

    `COWPY` est un outil en ligne de commande qui génère de l'art ASCII pour afficher des textes arbitraires de manière amusante.
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

De plus, nous allons utiliser un outil d'analyse linguistique appelé `langid` pour identifier la langue parlée par chaque personnage et organiser les sorties du pipeline en conséquence.

#### Examiner l'exercice

Votre défi est d'écrire un workflow Nextflow qui va :

1. **Générer de l'art ASCII** de chaque personnage
2. **Organiser** les sorties par famille linguistique (langues germaniques vs langues romanes)

Cela représente un schéma de workflow typique où les métadonnées propres à chaque fichier guident les décisions de traitement ; exactement le type de problème que les meta maps résolvent élégamment.

#### Liste de vérification

Vous pensez être prêt·e à vous lancer ?

- [ ] Je comprends l'objectif de ce cours et ses prérequis
- [ ] Mon codespace est opérationnel
- [ ] J'ai défini mon répertoire de travail de manière appropriée
- [ ] Je comprends l'exercice

Si vous pouvez cocher toutes les cases, vous êtes prêt·e à commencer.

---

## 1. Options de base pour charger et utiliser les métadonnées

Ouvrez le fichier de workflow `main.nf` pour examiner l'ébauche de workflow que nous vous fournissons comme point de départ.

```groovy title="main.nf" linenums="1"
#!/usr/bin/env nextflow

workflow  {
    main:
    ch_datasheet = channel.fromPath("./data/datasheet.csv")
        .splitCsv(header: true)
        .view()

    publish:
    cowpy_art = channel.empty()
}

output {
    cowpy_art {
    }
}
```

L'opérateur [`splitCsv`](https://www.nextflow.io/docs/latest/reference/operator.html#splitcsv) lit chaque ligne du fichier comme un élément dans le canal.
C'est la même approche que nous utilisons pour charger des données CSV dans Hello Nextflow, notre cours pour débutant·es.
Consultez [cette section](../../hello_nextflow/02_hello_channels.md#4-read-input-values-from-a-csv-file) si vous avez besoin d'un rappel sur son fonctionnement.

Avec `header: true`, la première ligne est traitée comme en-têtes de colonnes, de sorte que chaque élément devient une map de paires clé-valeur indexées par nom de colonne.

Notez que puisque nous n'exécutons pas encore de processus sur les données, les blocs `publish` et `output` ne sont que des ébauches.

### 1.1. Exécuter le workflow

Exécutez le workflow pour voir comment le contenu du canal est structuré une fois tout chargé :

```bash
nextflow run main.nf
```

??? success "Sortie de la commande"

    ```console
     N E X T F L O W   ~  version 25.10.4

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

Cela facilite l'accès à des champs spécifiques de chaque ligne.
Par exemple, nous pourrions accéder à l'identifiant du fichier avec `id` ou au chemin du fichier txt avec `recording`.

??? info "(Optionnel) En savoir plus sur les maps Groovy"

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
     N E X T F L O W   ~  version 25.10.4

    Launching `map_demo.nf` [cheesy_plateau] DSL2 - revision: fae5b8496e

    map: [id:sampleA, character:squirrel]
    id: sampleA
    character: squirrel
    ```

### 1.2. Sélectionner un champ spécifique avec `map`

Nous allons utiliser l'opérateur `map` pour itérer sur chaque élément d'un canal et extraire uniquement le champ `character`, auquel nous pouvons accéder par son nom en utilisant la notation pointée.

#### 1.2.1. Ajouter l'opération map

Pour accéder à la colonne `character`, ajoutez l'opération `map` avant l'opération `.view()` comme suit :

=== "Après"

    ```groovy title="main.nf" linenums="5" hl_lines="3-5"
        ch_datasheet = channel.fromPath("./data/datasheet.csv")
            .splitCsv(header: true)
            .map { row ->
                row.character
            }
            .view()
    ```

=== "Avant"

    ```groovy title="main.nf" linenums="5"
        ch_datasheet = channel.fromPath("./data/datasheet.csv")
            .splitCsv(header: true)
            .view()
    ```

Cette façon d'accéder à un champ spécifique est expliquée plus en détail dans [cette section](../../hello_nextflow/02_hello_channels.md#43-use-the-map-operator-to-extract-the-greetings) de Hello Nextflow, si vous avez besoin d'un rappel sur son fonctionnement.

#### 1.2.2. Exécuter le workflow

Exécutez le workflow pour vérifier que vous pouvez afficher les noms de personnages extraits.

```bash
nextflow run main.nf
```

??? success "Sortie de la commande"

    ```console
     N E X T F L O W   ~  version 25.10.4

    Launching `main.nf` [exotic_albattani] DSL2 - revision: c0d03cec83

    squirrel
    tux
    sheep
    turkey
    stegosaurus
    moose
    turtle
    ```

Cela montre que nous sommes capables d'accéder aux valeurs de la colonne `character` pour chaque ligne.

Maintenant, faisons quelque chose avec ces données : utilisons les champs `character` et `recording` ensemble pour générer de l'art ASCII avec `COWPY`.

### 1.3. Émettre des sous-canaux avec `multiMap`

Nous vous fournissons un module de processus `COWPY` pré-écrit, vous devez donc d'abord examiner les exigences d'entrée du processus.

Vous pouvez ouvrir le fichier pour voir à quoi ressemble le processus :

```groovy title="modules/cowpy.nf" linenums="1" hl_lines="7 8"
// Générer de l'art ASCII avec cowpy
process COWPY {

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

Comme vous pouvez le voir, le processus prend deux entrées séparées : un fichier d'enregistrement et un nom de personnage.
Il est important de noter que nous avons des valeurs pour les deux, mais elles sont actuellement regroupées à l'intérieur de chaque élément du canal.

Une façon d'extraire plusieurs champs dans des canaux séparés est l'opérateur [`multiMap`](https://www.nextflow.io/docs/latest/reference/operator.html#multimap), qui divise un canal en plusieurs sous-canaux nommés en une seule opération.

#### 1.3.1. Ajouter l'opération multiMap

Remplacez l'opération `map` par `multiMap` :

=== "Après"

    ```groovy title="main.nf" linenums="5" hl_lines="3-6"
        ch_datasheet = channel.fromPath("./data/datasheet.csv")
            .splitCsv(header: true)
            .multiMap { row ->
                file: row.recording
                character: row.character
            }
    ```

=== "Avant"

    ```groovy title="main.nf" linenums="5" hl_lines="3-6"
        ch_datasheet = channel.fromPath("./data/datasheet.csv")
            .splitCsv(header: true)
            .map { row ->
                row.character
            }
            .view()
    ```

Le bloc `multiMap` définit deux sous-canaux nommés (`file` et `character`) à partir de chaque ligne, auxquels nous pouvons accéder en tant que `ch_datasheet.file` et `ch_datasheet.character`.

#### 1.3.2. Appeler COWPY sur les sous-canaux

Maintenant, incluez le processus `COWPY` et passez-lui chaque sous-canal comme argument séparé :

=== "Après"

    ```groovy title="main.nf" linenums="1" hl_lines="3 14"
    #!/usr/bin/env nextflow

    include { COWPY } from './modules/cowpy.nf'

    workflow {
        main:
        ch_datasheet = channel.fromPath("./data/datasheet.csv")
            .splitCsv(header: true)
            .multiMap { row ->
                file: row.recording
                character: row.character
            }

        COWPY(ch_datasheet.file, ch_datasheet.character)

        publish:
        cowpy_art = channel.empty()
    }

    output {
        cowpy_art {
        }
    }
    ```

=== "Avant"

    ```groovy title="main.nf" linenums="1"
    #!/usr/bin/env nextflow

    workflow {
        main:
        ch_datasheet = channel.fromPath("./data/datasheet.csv")
            .splitCsv(header: true)
            .multiMap { row ->
                file: row.recording
                character: row.character
            }

        publish:
        cowpy_art = channel.empty()
    }

    output {
        cowpy_art {
        }
    }
    ```

Cela nous permet de passer les deux champs séparément comme `COWPY` l'exige.

#### 1.3.3. Configurer la publication des sorties

Enfin, ajoutez la sortie de `COWPY` au bloc `publish:` :

=== "Après"

    ```groovy title="main.nf" linenums="14" hl_lines="4"
        COWPY(ch_datasheet.file, ch_datasheet.character)

        publish:
        cowpy_art = COWPY.out
    ```

=== "Avant"

    ```groovy title="main.nf" linenums="14" hl_lines="4"
        COWPY(ch_datasheet.file, ch_datasheet.character)

        publish:
        cowpy_art = channel.empty()
    ```

Cela nous permettra de visualiser facilement les sorties produites par le workflow.

#### 1.3.4. Exécuter le workflow

Exécutez le workflow pour vérifier que `COWPY` s'exécute sur les entrées que nous avons fournies :

```bash
nextflow run main.nf
```

??? success "Sortie de la commande"

    ```console
     N E X T F L O W   ~  version 25.10.2

    Launching `main.nf` [clever_dijkstra] DSL2 - revision: a1b2c3d4e5

    executor >  local (7)
    [3a/f1c290] COWPY (7) [100%] 7 of 7 ✔
    ```

Comme vous pouvez le voir, `COWPY` s'est exécuté sur chaque fichier en utilisant le personnage correct pour chacun.

??? abstract "Contenu du répertoire de résultats"

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

??? example "Contenu de results/cowpy-guten_tag.txt"

    ```console
    $ cat results/cowpy-guten_tag.txt
     _____________________________
    / Guten Tag, wie geht es dir? \
    \ Auf Wiedersehen, bis morgen /
     -----------------------------
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

Cette approche fonctionne, mais présente une limitation : nous avons dû diviser le canal en deux sous-canaux séparés.
Si nous voulions passer davantage de champs au processus, nous devrions les extraire dans encore plus de sous-canaux.
Cela pourrait devenir fastidieux et difficile à maintenir.

Bonne nouvelle : il existe une façon plus simple de faire cela.

### 1.4. Regrouper tout comme une seule entrée du processus

Plutôt que de diviser les champs en canaux séparés, nous pouvons mettre à jour le processus pour recevoir toutes les entrées sous forme d'un seul tuple, ce qui simplifie l'appel au processus.

#### 1.4.1. Mettre à jour le processus COWPY

Mettez à jour `COWPY` pour accepter un tuple correspondant aux trois éléments de chaque ligne :

=== "Après"

    ```groovy title="modules/cowpy.nf" linenums="1" hl_lines="7 10 14"
    // Générer de l'art ASCII avec cowpy
    process COWPY {

        container 'community.wave.seqera.io/library/cowpy:1.1.5--3db457ae1977a273'

        input:
        tuple val(id), val(character), path(recording)

        output:
        path "cowpy-${recording}"

        script:
        """
        cat ${recording} | cowpy -c ${character} > cowpy-${recording}
        """
    }
    ```

=== "Avant"

    ```groovy title="modules/cowpy.nf" linenums="1" hl_lines="7-8 11 15"
    // Générer de l'art ASCII avec cowpy
    process COWPY {

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

Le processus ne prend maintenant qu'une seule entrée contenant tout ce que nous pourrions vouloir lui fournir.

#### 1.4.2. Utiliser `map()` pour créer le tuple d'entrée

Nous devons toujours utiliser une opération de mapping pour énumérer les éléments que nous voulons passer dans le tuple au processus :

=== "Après"

    ```groovy title="main.nf" linenums="5" hl_lines="3-5"
        ch_datasheet = channel.fromPath("./data/datasheet.csv")
            .splitCsv(header: true)
            .map { row ->
                [row.id, row.character, row.recording]
            }
    ```

=== "Avant"

    ```groovy title="main.nf" linenums="5" hl_lines="3-6"
        ch_datasheet = channel.fromPath("./data/datasheet.csv")
            .splitCsv(header: true)
            .multiMap { row ->
                file: row.recording
                character: row.character
            }
    ```

Vous vous demandez peut-être pourquoi nous ne pouvons pas simplement passer la map Groovy entière provenant de `splitCsv` telle quelle.
C'est parce que nous devons indiquer explicitement à Nextflow que le fichier d'enregistrement doit être traité comme un chemin (c'est-à-dire qu'il doit être correctement mis en scène).
Cela se produit au niveau de l'interface d'entrée de `COWPY`, où l'élément `recording` est explicitement désigné comme un `path`.

#### 1.4.3. Mettre à jour l'appel au processus

Enfin, remplaçons les deux entrées séparées dans l'appel au processus par le tuple unique que nous venons de créer :

=== "Après"

    ```groovy title="main.nf" linenums="5" hl_lines="7"
        ch_datasheet = channel.fromPath("./data/datasheet.csv")
            .splitCsv(header: true)
            .map { row ->
                [row.id, row.character, row.recording]
            }

        COWPY(ch_datasheet)
    ```

=== "Avant"

    ```groovy title="main.nf" linenums="5" hl_lines="7"
        ch_datasheet = channel.fromPath("./data/datasheet.csv")
            .splitCsv(header: true)
            .map { row ->
                [row.id, row.character, row.recording]
            }

        COWPY(ch_datasheet.file, ch_datasheet.character)
    ```

Cela simplifie légèrement l'appel au processus.

#### 1.4.4. Exécuter le workflow

Exécutez le workflow pour vérifier que `COWPY` peut toujours traiter les données correctement :

```bash
nextflow run main.nf
```

??? success "Sortie de la commande"

    ```console
     N E X T F L O W   ~  version 25.10.2

    Launching `main.nf` [pedantic_lovelace] DSL2 - revision: b2c3d4e5f6

    executor >  local (7)
    [5e/2a1b34] COWPY (7) [100%] 7 of 7 ✔
    ```

La sortie est la même que les sept fichiers `cowpy-*.txt` qu'auparavant, maintenant produits avec un appel plus simple à `COWPY`.

??? abstract "Contenu du répertoire de résultats"

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

??? example "Contenu de results/cowpy-guten_tag.txt"

    ```console
    $ cat results/cowpy-guten_tag.txt
     _____________________________
    / Guten Tag, wie geht es dir? \
    \ Auf Wiedersehen, bis morgen /
     -----------------------------
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

C'est une légère amélioration par rapport à l'approche `multiMap`.
Mais nous avons toujours dû décompresser la map Groovy originale pour créer le tuple d'entrée, et il existe un couplage fort entre le processus et la feuille de données : la définition des entrées de `COWPY` référence maintenant directement les noms de colonnes `id`, `character` et `recording`.

```groovy
input:
tuple val(id), val(character), path(recording)
```

Si un·e collaborateur·trice utilise une feuille de données structurée différemment — avec des colonnes supplémentaires, ou des colonnes dans un ordre différent — ce processus ne fonctionnera pas sans modification.
Cela rend le processus fragile, car sa structure d'entrée est liée à la composition exacte de la feuille de données.

Pour résoudre ce problème, nous avons besoin d'un moyen de passer toutes les métadonnées sous forme de bundle sans coder en dur sa structure exacte dans l'interface du processus.

### 1.5. Utiliser une interface meta map + fichier

La solution consiste à séparer deux préoccupations distinctes dans le canal : les **métadonnées sur un échantillon**, et le **fichier de données** lui-même.
En regroupant toutes les métadonnées dans une seule map — la « meta map » — nous obtenons un tuple à deux éléments cohérent, quel que soit le nombre de colonnes de métadonnées que contient la feuille de données :

```groovy title="Syntax example"
input:
tuple val(meta), path(file)
```

L'ajout ou la suppression de colonnes dans la feuille de données modifie le contenu de `meta`, mais la forme du tuple `[meta, file]` reste constante.
Les processus qui acceptent cette structure n'ont pas besoin de savoir ni de se soucier du nombre de champs de métadonnées existants.

#### 1.5.1. Réorganiser le contenu du tuple en meta map

Restructurons l'opération `map` pour produire un tuple `[meta, file]` :

=== "Après"

    ```groovy title="main.nf" linenums="5" hl_lines="4 6 8 11"
        ch_datasheet = channel.fromPath("./data/datasheet.csv")
            .splitCsv(header: true)
            .map { row ->
                [[id: row.id, character: row.character], row.recording]
            }
            .view()

        // COWPY(ch_datasheet)  // Sera mis à jour à l'étape suivante

        publish:
        cowpy_art = channel.empty() // COWPY.out
    ```

=== "Avant"

    ```groovy title="main.nf" linenums="5" hl_lines="4 7"
        ch_datasheet = channel.fromPath("./data/datasheet.csv")
            .splitCsv(header: true)
            .map { row ->
                [row.id, row.character, row.recording]
            }

        COWPY(ch_datasheet)

        publish:
        cowpy_art = COWPY.out
    ```

Vous remarquerez que nous avons également ajouté une instruction `view()`, mis en commentaire l'appel à `COWPY` et remplacé `COWPY.out` par `channel.empty()` car la définition des entrées du processus ne correspond pas encore à la nouvelle structure.

#### 1.5.2. Exécuter le workflow pour inspecter le contenu réorganisé

Exécutez le workflow pour voir la nouvelle forme du canal :

```bash
nextflow run main.nf
```

??? success "Sortie de la commande"

    ```console title="View meta map"
     N E X T F L O W   ~  version 25.10.4

    Launching `main.nf` [lethal_booth] DSL2 - revision: 0d8f844c07

    [[id:sampleA, character:squirrel], /workspaces/training/side-quests/metadata/data/bonjour.txt]
    [[id:sampleB, character:tux], /workspaces/training/side-quests/metadata/data/guten_tag.txt]
    [[id:sampleC, character:sheep], /workspaces/training/side-quests/metadata/data/hallo.txt]
    [[id:sampleD, character:turkey], /workspaces/training/side-quests/metadata/data/hello.txt]
    [[id:sampleE, character:stegosaurus], /workspaces/training/side-quests/metadata/data/hola.txt]
    [[id:sampleF, character:moose], /workspaces/training/side-quests/metadata/data/salut.txt]
    [[id:sampleG, character:turtle], /workspaces/training/side-quests/metadata/data/ciao.txt]
    ```

Chaque élément du canal est maintenant un tuple à deux éléments : la meta map en premier, le fichier en second.

```console title="Example element structure"
[
  [id:sampleA, character:squirrel],
  /workspaces/training/side-quests/metadata/data/bonjour.txt
]
```

Si nous ajoutons ultérieurement une colonne `language` à la feuille de données, elle sera disponible en tant que `meta.language` sans nécessiter de modifications de la définition des entrées du processus.

#### 1.5.3. Mettre à jour le processus `COWPY` pour utiliser la meta map

Mettez à jour `COWPY` pour accepter la structure tuple `[meta, file]` :

=== "Après"

    ```groovy title="modules/cowpy.nf" linenums="1" hl_lines="7 10 14"
    // Générer de l'art ASCII avec cowpy
    process COWPY {

        container 'community.wave.seqera.io/library/cowpy:1.1.5--3db457ae1977a273'

        input:
        tuple val(meta), path(input_file)

        output:
        path "cowpy-${input_file}"

        script:
        """
        cat ${input_file} | cowpy -c ${meta.character} > cowpy-${input_file}
        """
    }
    ```

=== "Avant"

    ```groovy title="modules/cowpy.nf" linenums="1" hl_lines="7 10 14"
    // Générer de l'art ASCII avec cowpy
    process COWPY {

        container 'community.wave.seqera.io/library/cowpy:1.1.5--3db457ae1977a273'

        input:
        tuple val(id), val(character), path(recording)

        output:
        path "cowpy-${recording}"

        script:
        """
        cat ${recording} | cowpy -c ${character} > cowpy-${recording}
        """
    }
    ```

Dans le bloc script, `meta.character` accède au champ `character` depuis la meta map.
Tout champ de la meta map est accessible de la même façon.

#### 1.5.4. Mettre à jour l'appel au processus

Restaurez l'appel à `COWPY` et connectez sa sortie pour la publication :

=== "Après"

    ```groovy title="main.nf" linenums="5" hl_lines="7 10"
        ch_datasheet = channel.fromPath("./data/datasheet.csv")
            .splitCsv(header: true)
            .map { row ->
                [[id: row.id, character: row.character], row.recording]
            }

        COWPY(ch_datasheet)

        publish:
        cowpy_art = COWPY.out
    ```

=== "Avant"

    ```groovy title="main.nf" linenums="5" hl_lines="6 8 11"
        ch_datasheet = channel.fromPath("./data/datasheet.csv")
            .splitCsv(header: true)
            .map { row ->
                [[id: row.id, character: row.character], row.recording]
            }
            .view()

        // COWPY(ch_datasheet)  // Sera mis à jour à l'étape suivante

        publish:
        cowpy_art = channel.empty() // COWPY.out
    ```

Nous avons également restauré la publication des sorties.

#### 1.5.5. Exécuter le workflow

Exécutez le workflow pour vérifier que tout fonctionne :

```bash
nextflow run main.nf
```

??? success "Sortie de la commande"

    ```console
     N E X T F L O W   ~  version 25.10.2

    Launching `main.nf` [wise_sammet] DSL2 - revision: 99797b1e92

    executor >  local (7)
    [5d/dffd4e] COWPY (7) [100%] 7 of 7 ✔
    ```

Le répertoire des résultats contient maintenant les fichiers d'art ASCII.

??? abstract "Contenu du répertoire"

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

??? example "Contenu de results/cowpy-guten_tag.txt"

    ```console
    $ cat results/cowpy-guten_tag.txt
     _____________________________
    / Guten Tag, wie geht es dir? \
    \ Auf Wiedersehen, bis morgen /
     -----------------------------
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

Le processus reçoit maintenant toutes les métadonnées sous forme de bundle via `meta`, utilise ce dont il a besoin (`meta.character`), et ignore le reste.

Il s'agit de l'interface standard utilisée par tous les modules [nf-core](https://nf-co.re/).
Le schéma `tuple val(meta), path(file)` apparaît de manière cohérente dans toute la bibliothèque de modules nf-core, c'est pourquoi les workflows qui adoptent cette convention peuvent intégrer des modules nf-core avec un minimum de friction.

### À retenir

Dans cette section, vous avez appris :

- **Comment lire des feuilles de données :** Utiliser `splitCsv` pour analyser des fichiers CSV avec des informations d'en-tête
- **Pourquoi la convention meta map existe :** Séparer les métadonnées des fichiers de données en tuples `[meta, file]` maintient la structure du canal stable à mesure que la feuille de données évolue
- **Comment utiliser les champs de la meta map dans un processus :** Tout champ de la meta map est accessible via la notation pointée dans le bloc script

---

## 2. Manipulations supplémentaires des métadonnées

Maintenant que l'interface meta map est en place, nous pouvons l'enrichir au fur et à mesure que les données circulent dans le pipeline.

Nous allons utiliser un outil appelé [`langid`](https://github.com/saffsd/langid.py) pour identifier la langue dans chaque fichier d'enregistrement.
À partir d'un extrait de texte, il produit une prédiction de langue et un score de probabilité vers `stdout`.

### 2.1. Ajouter une étape d'identification de la langue

Nous vous fournissons un module de processus pré-écrit appelé `IDENTIFY_LANGUAGE` qui encapsule l'outil `langid`.

Ouvrez le fichier module pour examiner son code :

```groovy title="modules/langid.nf" linenums="1" hl_lines="7 10"
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

La définition des entrées utilise la même structure `tuple val(meta), path(file)` que nous venons de construire dans la section 1, donc `ch_datasheet` peut alimenter directement ce processus sans aucune adaptation.

La sortie ajoute `stdout` comme troisième élément : cela capture la prédiction de langue que `langid` affiche sur la console.
La commande `sed` supprime le score de probabilité et le saut de ligne final, ne laissant que le code de langue à deux lettres.

#### 2.1.1. Ajouter un appel à `IDENTIFY_LANGUAGE`

Incluez le module de processus `IDENTIFY_LANGUAGE` et appelez-le sur le canal de la feuille de données :

=== "Après"

    ```groovy title="main.nf" linenums="1" hl_lines="4 14-16"
    #!/usr/bin/env nextflow

    include { COWPY } from './modules/cowpy.nf'
    include { IDENTIFY_LANGUAGE } from './modules/langid.nf'

    workflow {
        main:
        ch_datasheet = channel.fromPath("./data/datasheet.csv")
            .splitCsv(header: true)
            .map { row ->
                [[id: row.id, character: row.character], row.recording]
            }

        // Exécuter langid pour identifier la langue de chaque salutation
        IDENTIFY_LANGUAGE(ch_datasheet)
        IDENTIFY_LANGUAGE.out.view()

        COWPY(ch_datasheet)
    ```

=== "Avant"

    ```groovy title="main.nf" linenums="1"
    #!/usr/bin/env nextflow

    include { COWPY } from './modules/cowpy.nf'

    workflow {
        main:
        ch_datasheet = channel.fromPath("./data/datasheet.csv")
            .splitCsv(header: true)
            .map { row ->
                [[id: row.id, character: row.character], row.recording]
            }

        COWPY(ch_datasheet)
    ```

La sortie principale de ce processus est simplement une chaîne de caractères, il n'y a donc pas de fichiers de sortie à publier.
À la place, nous utilisons `IDENTIFY_LANGUAGE.out.view()` pour visualiser les résultats de l'opération.

#### 2.1.2. Exécuter le workflow

Exécutez le workflow pour produire l'identification de la langue, en utilisant `-resume` pour éviter de réexécuter les tâches `COWPY` :

```bash
nextflow run main.nf -resume
```

??? success "Sortie de la commande"

    ```console
     N E X T F L O W   ~  version 25.10.4

    Launching `main.nf` [voluminous_mcnulty] DSL2 - revision: f9bcfebabb

    executor >  local (14)
    [5d/dffd4e] COWPY (7)             [100%] 7 of 7, cached: 7 ✔
    [4e/f722fe] IDENTIFY_LANGUAGE (7) [100%] 7 of 7 ✔
    [[id:sampleA, character:squirrel], /workspaces/training/side-quests/metadata/work/eb/f7148ebdd898fbe1136bec6a714acb/bonjour.txt, fr]
    [[id:sampleB, character:tux], /workspaces/training/side-quests/metadata/work/16/71d72410952c22cd0086d9bca03680/guten_tag.txt, de]
    [[id:sampleD, character:turkey], /workspaces/training/side-quests/metadata/work/c4/b7562adddc1cc0b7d414ec45d436eb/hello.txt, en]
    [[id:sampleC, character:sheep], /workspaces/training/side-quests/metadata/work/ea/04f5d979429e4455e14b9242fb3b45/hallo.txt, de]
    [[id:sampleF, character:moose], /workspaces/training/side-quests/metadata/work/5a/6c2b84bf8fadb98e28e216426be079/salut.txt, fr]
    [[id:sampleE, character:stegosaurus], /workspaces/training/side-quests/metadata/work/af/ee7c69bcab891c40d0529305f6b9e7/hola.txt, es]
    [[id:sampleG, character:turtle], /workspaces/training/side-quests/metadata/work/4e/f722fe47271ba7ebcd69afa42964ca/ciao.txt, it]
    ```

Nous avons maintenant une prédiction de langue pour chaque fichier du jeu de données.

Notez que le tuple de sortie est composé de `[meta, file, lang_id]`, ce qui signifie que la meta map et le fichier sont transmis aux côtés du nouveau résultat.

!!! note "Note"

    Ce schéma consistant à maintenir la meta map associée aux résultats facilite l'association de résultats entre canaux ultérieurement.
    Vous ne pouvez pas vous fier à l'ordre des éléments dans les canaux pour associer correctement les données.
    Vous devez plutôt utiliser des clés.
    Les meta maps fournissent une structure idéale à cet effet.

    Ce cas d'utilisation est exploré en détail dans la quête secondaire [Splitting & Grouping](../splitting_and_grouping/index.md).

### 2.2. Enrichir les métadonnées avec les sorties des processus

La prédiction de langue est elle-même une métadonnée sur les données contenues dans le fichier.
Plutôt que de la conserver comme élément séparé, intégrons-la dans la meta map.

#### 2.2.1. Créer une nouvelle meta map enrichie

Nous pouvons créer une nouvelle meta map pour remplacer l'originale en utilisant l'opérateur Groovy `+` :

=== "Après"

    ```groovy title="main.nf" linenums="14" hl_lines="3-7"
        // Exécuter langid pour identifier la langue de chaque salutation
        IDENTIFY_LANGUAGE(ch_datasheet)
        IDENTIFY_LANGUAGE.out
            .map { meta, file, lang_id ->
                [meta + [lang: lang_id], file]
            }
            .view()
    ```

=== "Avant"

    ```groovy title="main.nf" linenums="14" hl_lines="3"
        // Exécuter langid pour identifier la langue de chaque salutation
        IDENTIFY_LANGUAGE(ch_datasheet)
        IDENTIFY_LANGUAGE.out.view()
    ```

Le cœur de cette opération est `#!groovy meta + [lang: lang_id]`.

Ce code crée essentiellement une map temporaire avec une seule paire clé-valeur contenant le code de langue (`[lang: lang_id]`), puis utilise l'opérateur Groovy `+` pour la combiner avec la map `meta` originale contenant les métadonnées préexistantes, produisant une nouvelle meta map enrichie.

Pour une explication plus détaillée, consultez l'encadré ci-dessous.

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
    new_map = map1 + [lang: lang_id]
    ```

    Ici, `[lang: lang_id]` crée une nouvelle map sans nom à la volée, et `map1 + ` fusionne `map1` avec la nouvelle map sans nom, produisant le même contenu de `new_map` qu'auparavant.

    Élégant, non ?

    **Transposons maintenant cela dans le contexte d'une opération Nextflow `channel.map()`.**

    Le code devient :

    ```groovy
    .map { map1, lang_id ->
        map1 + [lang: lang_id]
    }
    ```

    Cela fait ce qui suit :

    - `#!groovy map1, lang_id ->` prend les deux éléments du tuple
    - `#!groovy map1 + [lang: lang_id]` crée la nouvelle map comme détaillé ci-dessus

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

    Si vous avez du mal à comprendre pourquoi le `file` semble se déplacer dans l'opération `map`, imaginez qu'au lieu de `#!groovy [meta + [lang: lang_id], file]`, cette ligne se lise `[new_map, file]`.
    Cela devrait clarifier que nous laissons simplement le `file` à sa place d'origine en deuxième position dans le tuple. Nous avons juste pris la valeur `new_info` et l'avons intégrée dans la map qui est en première position.

    **Et cela nous ramène à la structure de canal `tuple val(meta), path(file)` !**

#### 2.2.2. Exécuter le workflow

Une fois que vous êtes certain·e de comprendre ce que fait le code, exécutez le workflow pour voir si cela fonctionne :

```bash
nextflow run main.nf -resume
```

??? success "Sortie de la commande"

    ```console
     N E X T F L O W   ~  version 25.10.4

    Launching `main.nf` [cheeky_fermat] DSL2 - revision: d096281ee4

    [5d/dffd4e] COWPY (7)             [100%] 7 of 7, cached: 7 ✔
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

!!! tip "Supprimer des clés d'une meta map"

    Vous pouvez supprimer une clé d'une meta map en utilisant la méthode Groovy [`subMap`](https://docs.groovy-lang.org/latest/html/groovy-jdk/java/util/Map.html#subMap(java.util.Collection)), qui retourne une nouvelle map contenant uniquement les clés que vous spécifiez :

    ```groovy
    meta.subMap(['id', 'character'])  // retourne une map avec uniquement 'id' et 'character'
    ```

    Cela est utile lorsqu'un processus ou module en aval n'a pas besoin de tous les champs qui se sont accumulés dans la meta map.

### 2.3. Assigner un groupe linguistique avec des conditions

Avec la prédiction de langue dans la meta map, nous pouvons en dériver d'autres métadonnées.
Les langues de notre jeu de données appartiennent à deux familles : germaniques (anglais, allemand) et romanes (français, espagnol, italien).
L'ajout d'un champ `lang_group` rendra cette classification disponible en aval.

#### 2.3.1. Ajouter une opération `map` avec la logique conditionnelle

Nous allons utiliser une deuxième opération `map` avec une logique conditionnelle pour assigner la famille linguistique :

```groovy
.map { meta, file ->

    // la logique conditionnelle définissant lang_group va ici

    [meta + [lang_group: lang_group], file]
}
```

Voici la logique à appliquer :

- Commencer avec `lang_group = 'unknown'` comme valeur par défaut.
- Si `meta.lang` est `'de'` ou `'en'`, définir `lang_group` à `'germanic'`.
- Sinon si `meta.lang` est dans `['fr', 'es', 'it']`, définir `lang_group` à `'romance'`.

!!! tip "Astuce"

    Vous pouvez accéder à la valeur de `lang` dans l'opération map avec `meta.lang`.

Effectuez les modifications suivantes dans le workflow :

=== "Après"

    ```groovy title="main.nf" linenums="14" hl_lines="7-19 21"
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

        ch_languages.view()
    ```

=== "Avant"

    ```groovy title="main.nf" linenums="14" hl_lines="7"
        // Exécuter langid pour identifier la langue de chaque salutation
        IDENTIFY_LANGUAGE(ch_datasheet)
        IDENTIFY_LANGUAGE.out
            .map { meta, file, lang_id ->
                [meta + [lang: lang_id], file]
            }
            .view()
    ```

Points clés :

- `def lang_group = "unknown"` initialise la variable avec une valeur par défaut sûre.
- La structure `if / else if` gère les deux familles linguistiques ; tout le reste reste à `'unknown'`.
- `#!groovy .set { ch_languages }` donne un nom au canal résultant pour l'utiliser à l'étape suivante.

<!-- TODO (future) Add note/links to relevant docs in additional resources section -->

#### 2.3.2. Exécuter le workflow

Exécutez le workflow pour vérifier que cela fonctionne :

```bash
nextflow run main.nf -resume
```

??? success "Sortie de la commande"

    ```console
     N E X T F L O W   ~  version 25.10.4

    Launching `main.nf` [wise_almeida] DSL2 - revision: 46778c3cd0

    [5d/dffd4e] COWPY (7)             [100%] 7 of 7, cached: 7 ✔
    [da/652cc6] IDENTIFY_LANGUAGE (7) [100%] 7 of 7, cached: 7 ✔
    [[id:sampleA, character:squirrel, lang:fr, lang_group:romance], /workspaces/training/side-quests/metadata/data/bonjour.txt]
    [[id:sampleB, character:tux, lang:de, lang_group:germanic], /workspaces/training/side-quests/metadata/data/guten_tag.txt]
    [[id:sampleC, character:sheep, lang:de, lang_group:germanic], /workspaces/training/side-quests/metadata/data/hallo.txt]
    [[id:sampleD, character:turkey, lang:en, lang_group:germanic], /workspaces/training/side-quests/metadata/data/hello.txt]
    [[id:sampleE, character:stegosaurus, lang:es, lang_group:romance], /workspaces/training/side-quests/metadata/data/hola.txt]
    [[id:sampleF, character:moose, lang:fr, lang_group:romance], /workspaces/training/side-quests/metadata/data/salut.txt]
    [[id:sampleG, character:turtle, lang:it, lang_group:romance], /workspaces/training/side-quests/metadata/data/ciao.txt]
    ```

La meta map contient maintenant quatre champs : `id`, `character`, `lang` et `lang_group`.
La structure du canal est toujours `[meta, file]`.

### 2.4. Utiliser les métadonnées pour nommer et organiser les sorties

Avec `lang` et `lang_group` maintenant disponibles dans la meta map, nous pouvons les utiliser pour ajouter un code de langue aux noms des fichiers de sortie et les organiser dans des sous-répertoires par famille linguistique.

Cela nécessite trois modifications : mettre à jour le processus `COWPY` pour renommer sa sortie et inclure `meta` dans ce qu'il émet, mettre à jour l'appel à `COWPY` pour s'exécuter sur `ch_languages`, et mettre à jour le bloc output pour spécifier le chemin du sous-répertoire.

#### 2.4.1. Mettre à jour le processus `COWPY`

Renommez le fichier de sortie en utilisant le code de langue de la meta map, et ajoutez `meta` à la sortie afin que le bloc output puisse accéder à `lang_group` pour le routage vers les sous-répertoires :

=== "Après"

    ```groovy title="modules/cowpy.nf" linenums="9" hl_lines="2 6"
        output:
        tuple val(meta), path("${meta.lang}-${input_file}")

        script:
        """
        cat ${input_file} | cowpy -c ${meta.character} > ${meta.lang}-${input_file}
        """
    ```

=== "Avant"

    ```groovy title="modules/cowpy.nf" linenums="9" hl_lines="2 6"
        output:
        path "cowpy-${input_file}"

        script:
        """
        cat ${input_file} | cowpy -c ${meta.character} > cowpy-${input_file}
        """
    ```

Cela montre comment nous pouvons tirer parti d'autres champs de métadonnées pour personnaliser le comportement d'un processus, sans avoir à modifier du tout la définition des entrées.

#### 2.4.2. Mettre à jour l'appel à `COWPY` pour s'exécuter sur `ch_languages`

Remplacez `COWPY(ch_datasheet)` par `COWPY(ch_languages)` :

=== "Après"

    ```groovy title="main.nf" linenums="32" hl_lines="3"
        .set { ch_languages }

        COWPY(ch_languages)

        publish:
        cowpy_art = COWPY.out
    }
    ```

=== "Avant"

    ```groovy title="main.nf" linenums="32" hl_lines="3 5"
        .set { ch_languages }

        ch_languages.view()

        COWPY(ch_datasheet)

        publish:
        cowpy_art = COWPY.out
    }
    ```

Nous supprimons également la ligne `ch_languages.view()` puisque nous n'avons plus besoin d'inspecter le contenu du canal.

#### 2.4.3. Mettre à jour le bloc output

Ajoutez une closure `path` au bloc `output {}` pour router chaque fichier dans son sous-répertoire de groupe linguistique :

=== "Après"

    ```groovy title="main.nf" linenums="40" hl_lines="3"
    output {
        cowpy_art {
            path { meta, file -> meta.lang_group }
        }
    }
    ```

=== "Avant"

    ```groovy title="main.nf" linenums="40" hl_lines="2 3"
    output {
        cowpy_art {
        }
    }
    ```

Cela montre comment nous pouvons utiliser les métadonnées pour organiser les sorties avec une grande flexibilité.

#### 2.4.4. Exécuter le pipeline complet

Supprimez les résultats précédents et exécutez le pipeline complet :

```bash
rm -r results
nextflow run main.nf
```

??? success "Sortie de la commande"

    ```console
     N E X T F L O W   ~  version 25.10.4

    Launching `main.nf` [suspicious_crick] DSL2 - revision: 25541014c5

    executor >  local (14)
    [5d/dffd4e] IDENTIFY_LANGUAGE (7) [100%] 7 of 7 ✔
    [e7/317c18] COWPY (7)             [100%] 7 of 7 ✔
    ```

Le répertoire des résultats est maintenant organisé par famille linguistique, avec chaque fichier nommé d'après sa langue détectée :

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

La closure `path` dans le bloc `output {}` reçoit chaque tuple `[meta, file]` et retourne `meta.lang_group` comme nom de sous-répertoire.
Le nom du fichier lui-même provient de ce que le processus produit (`#!groovy "${meta.lang}-${input_file}"`).
Les deux éléments de métadonnées (code de langue et groupe linguistique) proviennent de la meta map enrichie construite dans cette section.

### À retenir

Dans cette section, vous avez appris :

- **Comment enrichir la meta map avec les sorties des processus :** Ajouter de nouvelles clés avec `#!groovy meta + [clé: valeur]` maintient la structure du canal `[meta, file]` intacte tout en enrichissant les métadonnées.
- **Comment dériver des métadonnées à partir de métadonnées :** Une logique conditionnelle dans une opération `map` peut calculer de nouveaux champs à partir de champs existants.
- **Comment utiliser les métadonnées pour organiser les sorties :** La closure `path` dans le bloc `output {}` peut lire depuis la meta map pour router les fichiers dans des sous-répertoires.

---

## 3. Considérations de robustesse

Lorsque les valeurs de métadonnées pilotent le comportement des processus, des données manquantes ou incomplètes peuvent causer des problèmes difficiles à diagnostiquer.
Voici ce à quoi vous pouvez vous attendre et comment y faire face.

### 3.1. Que se passe-t-il lorsqu'un champ de métadonnées requis est manquant

La valeur `character` est requise pour que le processus `COWPY` produise un résultat valide.
Le mode d'échec dépend du fait que la colonne existe dans la feuille de données mais est vide, ou qu'elle est totalement absente.

#### 3.1.1. La colonne existe mais une valeur est vide

Supposons qu'une entrée de la feuille de données ait un champ `character` vide :

```csv title="datasheet.csv" linenums="1" hl_lines="2"
id,character,recording
sampleA,,/workspaces/training/side-quests/metadata/data/bonjour.txt
sampleB,tux,/workspaces/training/side-quests/metadata/data/guten_tag.txt
...
```

La clé `character` est créée pour toutes les entrées lors de l'analyse de la feuille de données, mais `meta.character` pour `sampleA` sera une chaîne vide.
Lorsque Nextflow substitue `#!groovy ${meta.character}` dans la commande, l'outil `COWPY` reçoit un argument vide pour `-c` et échoue :

??? failure "Sortie de la commande"

    ```console hl_lines="8 11 16 28"
     N E X T F L O W   ~  version 25.10.4

    Launching `main.nf` [marvelous_hirsch] DSL2 - revision: 0dfeee3cc1

    executor >  local (9)
    [c1/c5dd4f] process > IDENTIFY_LANGUAGE (7) [ 85%] 6 of 7
    [d3/b7c415] process > COWPY (2)             [  0%] 0 of 6
    ERROR ~ Error executing process > 'COWPY (1)'

    Caused by:
      Process `COWPY (1)` terminated with an error exit status (2)


    Command executed:

      cat bonjour.txt | cowpy -c  > fr-bonjour.txt

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

Le message d'erreur (`expected one argument`) pointe vers le drapeau `-c` vide.
L'inspection du fichier `.command.sh` dans le répertoire de travail confirme que la commande a été exécutée avec une valeur vide.

#### 3.1.2. La colonne n'existe pas dans la feuille de données

Si la colonne `character` est totalement absente :

```csv title="datasheet.csv" linenums="1"
id,recording
sampleA,/workspaces/training/side-quests/metadata/data/bonjour.txt
sampleB,/workspaces/training/side-quests/metadata/data/guten_tag.txt
...
```

La clé `character` n'est jamais créée dans la meta map.
Lorsque le script du processus évalue `#!groovy ${meta.character}`, la clé manquante retourne `null`, et Nextflow substitue littéralement la chaîne `null` dans la commande :

??? failure "Sortie de la commande"

    ```console hl_lines="8 11 16"
     N E X T F L O W   ~  version 25.10.4

    Launching `main.nf` [jovial_bohr] DSL2 - revision: eaaf375827

    executor >  local (9)
    [0d/ada9db] process > IDENTIFY_LANGUAGE (5) [ 85%] 6 of 7
    [06/28065f] process > COWPY (2)             [  0%] 0 of 6
    ERROR ~ Error executing process > 'COWPY (2)'

    Caused by:
      Process `COWPY (2)` terminated with an error exit status (1)


    Command executed:

      cat guten_tag.txt | cowpy -c null > de-guten_tag.txt

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

Le `cowpy -c null` dans la commande exécutée est l'indice diagnostique.

### 3.2. Stratégies pour gérer les métadonnées manquantes

Il existe deux approches complémentaires pour rendre les workflows plus robustes face aux métadonnées manquantes.

**1. Validation des entrées**

La solution la plus fiable est de valider la feuille de données avant tout traitement, afin que les problèmes soient détectés tôt avec un message d'erreur clair plutôt que de se manifester comme un échec cryptique de processus en cours d'exécution.
La formation [Hello nf-core](../../hello_nf-core/05_input_validation.md) explique comment ajouter une validation des entrées en utilisant le plugin nf-schema. <!-- TODO (future) pending a proper Validation side quest -->

**2. Entrées de processus explicites pour les valeurs requises**

Si vous souhaitez que l'interface du processus elle-même indique qu'une valeur particulière est obligatoire, envisagez de l'extraire de la meta map comme entrée explicite :

=== "Définition du processus"

    ```groovy title="modules/cowpy.nf" linenums="6"
    input:
    tuple val(meta), val(character), path(input_file)
    ```

=== "Appel dans le workflow"

    ```groovy title="main.nf"
    COWPY(ch_languages.map { meta, file -> [meta, meta.character, file] })
    ```

Cette approche fait de `character` une partie visible et requise du contrat du processus.
Toute personne lisant le module peut immédiatement voir qu'une valeur de personnage doit être fournie.
Si le champ est absent, le workflow échoue clairement au niveau du canal avant même que le processus ne s'exécute.

Cela met en évidence un principe de conception utile :

**Utilisez la meta map pour les informations optionnelles ou descriptives ; extrayez les valeurs requises comme entrées explicites.**

La meta map maintient des structures de canaux propres et stables, mais pour les valeurs qui sont véritablement requises par un processus, les exposer comme entrées nommées améliore la clarté et facilite l'utilisation correcte du module dans d'autres contextes.

### À retenir

Dans cette section, vous avez vu :

- **Comment se manifestent les métadonnées manquantes :** Un champ vide produit un argument vide ; un champ absent produit `null` substitué littéralement dans la commande.
- **Deux stratégies complémentaires :** La validation des entrées pour détecter les problèmes tôt, et les entrées de processus explicites pour communiquer clairement les exigences.

---

## Résumé

Dans cette quête secondaire, vous avez exploré comment travailler efficacement avec les métadonnées dans les workflows Nextflow.

Le schéma de tuple « meta map + fichier de données » est une convention fondamentale dans Nextflow, offrant plusieurs avantages par rapport au passage des métadonnées comme valeurs individuelles :

- La structure du canal reste stable à mesure que la feuille de données évolue
- Le comportement des processus peut être personnalisé par échantillon sans coder en dur les noms de champs
- Les métadonnées sont disponibles tout au long du pipeline pour nommer, regrouper et organiser les sorties
- Les modules écrits pour cette interface sont interchangeables, y compris les modules nf-core

### Schémas clés

1.  **Lecture et structuration des métadonnées :** Analyser une feuille de données CSV et créer une meta map.

    ```groovy
    channel.fromPath('datasheet.csv')
        .splitCsv(header: true)
        .map { row ->
            [ [id: row.id, character: row.character], row.recording ]
        }
    ```

2.  **Enrichissement des métadonnées pendant le workflow :** Ajouter de nouvelles clés à partir des sorties de processus ou d'une logique dérivée.

    ```groovy
    // À partir d'une sortie de processus
    .map { meta, file, lang ->
        [ meta + [lang: lang], file ]
    }

    // À partir d'une logique conditionnelle
    .map { meta, file ->
        def lang_group = "unknown"
        if (meta.lang in ["de", "en"]) { lang_group = "germanic" }
        else if (meta.lang in ["fr", "es", "it"]) { lang_group = "romance" }
        [ meta + [lang_group: lang_group], file ]
    }
    ```

3.  **Utilisation des métadonnées dans un processus :** Accéder à n'importe quel champ via la notation pointée dans le bloc script.

    ```groovy
    cat ${input_file} | cowpy -c ${meta.character} > ${meta.lang}-${input_file}
    ```

4.  **Organisation des sorties par valeur de métadonnée :** Utiliser la closure `path` dans le bloc `output {}`.

    ```groovy
    output {
        cowpy_art {
            path { meta, file -> meta.lang_group }
        }
    }
    ```

### Ressources supplémentaires

- [opérateur map](https://www.nextflow.io/docs/latest/operator.html#map)
- [opérateur multiMap](https://www.nextflow.io/docs/latest/reference/operator.html#multimap)
- [qualificateur de sortie stdout](https://www.nextflow.io/docs/latest/process.html#outputs)

---

## Et ensuite ?

Retournez au [menu des Quêtes secondaires](../index.md) ou cliquez sur le bouton en bas à droite de la page pour passer au sujet suivant dans la liste.
