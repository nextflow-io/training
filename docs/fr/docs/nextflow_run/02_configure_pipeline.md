# Partie 2 : Configurer le pipeline

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Traduction assistée par IA - [en savoir plus et suggérer des améliorations](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

Dans la [Partie 1](./01_run_nextflow.md), vous avez exécuté un pipeline complet à plusieurs étapes qui traite plusieurs entrées en parallèle à l'aide de conteneurs.
Nous allons maintenant examiner comment configurer le comportement du pipeline à l'aide de `nextflow.config` : d'abord en analysant le fichier de configuration que nous vous avons déjà fourni, puis en explorant d'autres façons de fournir une configuration, et enfin en contrôlant comment et où les sorties sont publiées.

---

## 1. Examiner le fichier de configuration principal

Nextflow récupère automatiquement `nextflow.config` depuis le répertoire de travail et applique ses paramètres à chaque exécution.

Nous vous fournissons un fichier de configuration qui couvre quatre domaines : le packaging logiciel, les paramètres de processus, les paramètres du pipeline et les profils d'exécution.

??? full-code "nextflow.config"

    ```groovy title="nextflow.config" linenums="1"
    /*
     * Packaging logiciel
     */
    docker.enabled = true

    /*
     * Paramètres de processus
     */
    process {
        cpus = 1
        memory = 1.GB
    }

    /*
     * Paramètres du pipeline
     */
    params {
        input = 'data/greetings.csv'
        batch = 'batch'
        character = 'turkey'
    }

    /*
     * Profils
     */
    profiles {
        test {
            params.input = 'data/greetings.csv'
            params.batch = 'test'
            params.character = 'tux'
        }
        conda {
            docker.enabled = false
            conda.enabled = true
        }
    }
    ```

Passons en revue chacun d'eux, puis mettons les profils en pratique en exécutant le pipeline avec l'un d'eux.

!!! note "Note"

    Cette configuration couvre l'exécution locale sur une seule machine.
    Nextflow prend également en charge les ordonnanceurs HPC (SLURM, PBS, LSF) et les executors cloud (AWS Batch, Google Cloud Batch, Azure Batch), tous configurés via le même mécanisme `nextflow.config`.
    Consultez la [Partie 1 : S'adapter à votre environnement de calcul](../execution_config/01_packaging_and_execution.md) du cours [Execution Config](../execution_config/index.md) pour une présentation complète de ces options.

### 1.1. Packaging logiciel

Le packaging logiciel désigne la façon dont Nextflow fournit les outils dont vos processus ont besoin, qu'il s'agisse d'une image de conteneur, d'un environnement Conda ou d'autre chose.

```groovy title="nextflow.config" linenums="1"
/*
 * Packaging logiciel
 */
docker.enabled = true
```

Cette ligne active Docker pour chaque processus.
Tout processus qui déclare une directive `container` s'exécute à l'intérieur de l'image spécifiée.

### 1.2. Paramètres de processus

Rappelons qu'un processus est une étape unique de votre pipeline, comme `sayHello` ou `cowpy`.
Nextflow vous permet de configurer un certain nombre de choses sur la façon dont chacun s'exécute réellement : la quantité de CPU et de mémoire allouée, le conteneur ou l'environnement Conda utilisé, et bien plus encore.

```groovy title="nextflow.config" linenums="6"
/*
 * Paramètres de processus
 */
process {
    cpus = 1
    memory = 1.GB
}
```

Cela limite chaque processus à un seul CPU et 1 Go de mémoire.

Nextflow vous permet également de définir des valeurs différentes pour des processus nommés individuellement ou des groupes de processus ; vous apprendrez comment dans la [Partie 2 : Gérer les ressources de calcul et les échecs](../execution_config/02_resources_and_retries.md#12-set-resource-allocations-for-a-specific-process) du cours [Execution Config](../execution_config/index.md).

### 1.3. Paramètres du pipeline

Les paramètres sont les entrées en ligne de commande du pipeline, les mêmes options `--input`, `--batch` et `--character` que vous avez déjà définies directement en ligne de commande.
Définir leurs valeurs par défaut ici signifie que vous n'avez pas à les saisir à chaque fois, bien que, comme vous le verrez plus loin dans cette partie, il existe d'autres façons de les fournir.

```groovy title="nextflow.config" linenums="14"
/*
 * Paramètres du pipeline
 */
params {
    input = 'data/greetings.csv'
    batch = 'batch'
    character = 'turkey'
}
```

Ces valeurs par défaut s'appliquent chaque fois qu'un paramètre n'est pas fourni en ligne de commande, de sorte qu'exécuter `nextflow run main.nf` sans aucune option fonctionne quand même.

### 1.4. Profils

Les profils vous permettent de regrouper un ensemble de paramètres sous un seul nom, afin de pouvoir basculer entre des configurations complètes avec une seule option au lieu de modifier les valeurs manuellement à chaque fois.

```groovy title="nextflow.config" linenums="23"
/*
 * Profils
 */
profiles {
    test {
        params.input = 'data/greetings.csv'
        params.batch = 'test'
        params.character = 'tux'
    }
    conda {
        docker.enabled = false
        conda.enabled = true
    }
}
```

Le profil `test` remplace trois paramètres pour exécuter le pipeline avec un ensemble d'entrées petit et bien défini ; chaque pipeline nf-core en inclut un pour une validation rapide, et c'est une convention qui vaut la peine d'être suivie dans vos propres pipelines.

Le profil `conda` remplace le packaging logiciel Docker par Conda.

Vous activez un profil en passant `-profile <name>` en ligne de commande.

Mettons le profil `test` en pratique.

```bash
nextflow run main.nf -profile test
```

??? success "Sortie de la commande"

    ```console hl_lines="6 7 8 9"
    N E X T F L O W   ~  version 26.04.4

    Launching `main.nf` [reverent_heisenberg] revision: ce74f81996

    executor >  local (8)
    [3d/8a12c7] sayHello (3)       | 3 of 3 ✔
    [e7/e0934f] convertToUpper (2) | 3 of 3 ✔
    [1d/616569] collectGreetings   | 1 of 1 ✔
    [44/7d46cf] cowpy              | 1 of 1 ✔

    Outputs:

      /workspaces/training/nextflow-run/results

      first_output:
        - test/intermediates/Bonjour-output.txt
        - test/intermediates/Hello-output.txt
        - test/intermediates/Hola-output.txt

      uppercased:
        - test/intermediates/UPPER-Bonjour-output.txt
        - test/intermediates/UPPER-Hello-output.txt
        - test/intermediates/UPPER-Hola-output.txt

      collected: test/intermediates/COLLECTED-test-output.txt

      batch_report: test/test-report.txt

      cowpy_art: test/cowpy-COLLECTED-test-output.txt
    ```

Le pipeline s'exécute avec `batch = 'test'` et `character = 'tux'`.
Vérifiez `results/test/` : le nom du batch fait maintenant partie du chemin du répertoire lui-même, et l'art ASCII représente le pingouin tux au lieu d'une dinde.

!!! note "Note"

    Vous pouvez activer plusieurs profils à la fois, et utiliser `nextflow config -profile <name>,<name>` pour voir le résultat entièrement résolu avant d'exécuter quoi que ce soit.
    La combinaison de profils et la façon dont Nextflow résout les conflits entre eux sont traitées en détail dans la [Partie 3 : Utiliser les profils pour changer de configuration](../execution_config/03_profiles.md) du cours [Execution Config](../execution_config/index.md).

### À retenir

Vous savez ce que font les éléments les plus courants d'un fichier `nextflow.config` et comment activer un profil.

### Et ensuite ?

Découvrez d'autres façons de fournir des valeurs de configuration sans modifier le fichier `nextflow.config` principal, utiles pour configurer des exécutions individuelles et pour partager un ensemble précis de paramètres avec quelqu'un d'autre.

---

## 2. Fournir une configuration via des fichiers supplémentaires

Définir des valeurs par défaut dans `nextflow.config` fonctionne bien pour les valeurs qui changent rarement.
Nextflow vous offre également deux mécanismes plus ciblés : un fichier de configuration spécifique à une exécution pour adapter l'exécution à un environnement particulier, et un fichier de paramètres pour partager un ensemble précis de valeurs d'entrée avec un·e collaborateur·trice.

### 2.1. Utiliser un fichier de configuration spécifique à une exécution

Supposons que vous déplaciez le pipeline vers une machine qui n'a pas Docker et que vous souhaitiez donner à chaque processus plus de ressources.
Créez un nouveau fichier de configuration avec uniquement les remplacements dont vous avez besoin :

```groovy title="custom.config" linenums="1"
process {
    cpus = 2
    memory = 2.GB
}

docker.enabled = false
conda.enabled = true
```

Passez-le avec votre pipeline principal via `-c` :

```bash
nextflow run main.nf -c custom.config
```

??? success "Sortie de la commande"

    ```console hl_lines="6 7 8 9"
    N E X T F L O W   ~  version 26.04.4

    Launching `main.nf` [exotic_cray] revision: ce74f81996

    executor >  local (8)
    [77/e02315] sayHello (1)       | 3 of 3 ✔
    [a6/ccf44b] convertToUpper (3) | 3 of 3 ✔
    [56/fd1296] collectGreetings   | 1 of 1 ✔
    [2c/205a94] cowpy              | 1 of 1 ✔

    Outputs:

      /workspaces/training/nextflow-run/results

      first_output:
        - batch/intermediates/Bonjour-output.txt
        - batch/intermediates/Hello-output.txt
        - batch/intermediates/Hola-output.txt

      uppercased:
        - batch/intermediates/UPPER-Bonjour-output.txt
        - batch/intermediates/UPPER-Hello-output.txt
        - batch/intermediates/UPPER-Hola-output.txt

      collected: batch/intermediates/COLLECTED-batch-output.txt

      batch_report: batch/batch-report.txt

      cowpy_art: batch/cowpy-COLLECTED-batch-output.txt
    ```

Nextflow fusionne `custom.config` par-dessus le `nextflow.config` du pipeline, de sorte que chaque processus obtient maintenant 2 CPU et 2 Go de mémoire au lieu des valeurs par défaut, et s'exécute via Conda plutôt que Docker.
`cowpy` est le seul processus avec un package Conda déclaré aux côtés de son conteneur, c'est donc celui pour lequel vous verrez Nextflow construire réellement un environnement :

```console
Creating env using conda: conda-forge::cowpy==1.1.5 [cache /path/to/work/conda/env-898314d566668b6587ad714ae06b8520]
```

Un petit fichier qui ne remplace que l'allocation des ressources et le packaging, sans toucher aux paramètres du pipeline, correspond exactement au modèle attendu par les pipelines nf-core pour les configurations institutionnelles.
Parcourez le dépôt [nf-core/configs](https://github.com/nf-core/configs) pour des exemples concrets.

Cela vous offre un moyen jetable d'adapter un pipeline à un nouvel environnement sans toucher à votre configuration habituelle.

### 2.2. Utiliser un fichier de paramètres

Supposons maintenant que vous ayez besoin de partager un ensemble précis de paramètres d'exécution avec un·e collaborateur·trice, ou de les enregistrer pour une publication.

Nextflow vous permet de fournir des [fichiers de paramètres](https://nextflow.io/docs/latest/config.html#parameter-file) au format YAML ou JSON, qui constituent un moyen plus simple de distribuer un ensemble de valeurs exact et reproductible.

Un fichier de paramètres appelé `test-params.yaml` est déjà fourni dans votre répertoire de travail :

```yaml title="test-params.yaml" linenums="1"
input: "data/greetings.csv"
batch: "yaml"
character: "stegosaurus"
```

La syntaxe utilise des deux-points (`:`) au lieu des signes égal (`=`) utilisés dans `nextflow.config`, car ce fichier est du YAML pur plutôt que du Groovy.

!!! info "Info"

    Une version JSON, `test-params.json`, est également fournie. N'hésitez pas à l'essayer par vous-même ; la syntaxe pour la passer est identique.

Passez le fichier avec `-params-file` :

```bash
nextflow run main.nf -params-file test-params.yaml
```

??? success "Sortie de la commande"

    ```console
    N E X T F L O W   ~  version 26.04.4

    Launching `main.nf` [sharp_faraday] revision: ce74f81996

    executor >  local (8)
    [1c/9ff63e] sayHello (1)       | 3 of 3 ✔
    [3b/bb5691] convertToUpper (2) | 3 of 3 ✔
    [cd/2c1f6e] collectGreetings   | 1 of 1 ✔
    [89/c333bc] cowpy              | 1 of 1 ✔

    Outputs:

      /workspaces/training/nextflow-run/results

      first_output:
        - yaml/intermediates/Bonjour-output.txt
        - yaml/intermediates/Hello-output.txt
        - yaml/intermediates/Hola-output.txt

      uppercased:
        - yaml/intermediates/UPPER-Bonjour-output.txt
        - yaml/intermediates/UPPER-Hello-output.txt
        - yaml/intermediates/UPPER-Hola-output.txt

      collected: yaml/intermediates/COLLECTED-yaml-output.txt

      batch_report: yaml/yaml-report.txt

      cowpy_art: yaml/cowpy-COLLECTED-yaml-output.txt
    ```

??? abstract "Contenu du fichier"

    ```console title="results/yaml/cowpy-COLLECTED-yaml-output.txt"
     _________
    / BONJOUR \
    | HOLA    |
    \ HELLO   /
     ---------
    \                             .       .
     \                           / `.   .' "
      \                  .---.  <    > <    >  .---.
       \                 |    \  \ - ~ ~ - /  /    |
             _____          ..-~             ~-..-~
            |     |   \~~~\.'                    `./~~~/
           ---------   \__/                        \__/
          .'  O    \     /               /       \  "
         (_____,    `._.'               |         }  \/~~~/
          `----.          /       }     |        /    \__/
                `-.      |       /      |       /      `. ,~~|
                    ~-.__|      /_ - ~ ^|      /- _      `..-'
                        |     /        |     /     ~-.     `-. _  _  _
                        |_____|        |_____|         ~ - . _ _ _ _ _>
    ```

Un fichier de paramètres est particulièrement utile lorsqu'un pipeline possède plus d'une poignée de paramètres : il vous permet de les fournir tous en une seule fois, sans une ligne de commande interminable ni aucune modification du script du workflow, et il est facile à distribuer avec vos résultats.

### À retenir

Vous connaissez deux autres façons de fournir une configuration : un fichier de configuration spécifique à une exécution pour adapter l'exécution à un nouvel environnement, et un fichier de paramètres pour partager des valeurs d'entrée exactes et reproductibles.

### Et ensuite ?

Apprenez à contrôler comment et où les sorties de votre pipeline sont publiées.

---

## 3. Gérer les sorties du pipeline

L'auteur·e d'un pipeline décide de la façon dont les sorties sont organisées dans le code, mais vous n'avez pas besoin de toucher à ce code pour contrôler où elles se retrouvent ou comment elles y parviennent.
Nextflow vous offre des moyens de le faire au niveau de la configuration : définir un répertoire de sortie de base et choisir si les fichiers sont copiés ou liés symboliquement.

### 3.1. Personnaliser le répertoire de sortie

Par défaut, Nextflow publie les sorties dans `results/`.
Pointez vers un autre emplacement avec `-output-dir` (ou sa forme abrégée, `-o`) :

```bash
nextflow run main.nf -output-dir outputs
```

??? success "Sortie de la commande"

    ```console hl_lines="6 7 8 9"
    N E X T F L O W   ~  version 26.04.4

    Launching `main.nf` [serene_kimura] revision: ce74f81996

    executor >  local (8)
    [31/df5c15] sayHello (3)       | 3 of 3 ✔
    [08/8bb2e5] convertToUpper (1) | 3 of 3 ✔
    [e5/5814da] collectGreetings   | 1 of 1 ✔
    [cf/8ab8c6] cowpy              | 1 of 1 ✔

    Outputs:

      /workspaces/training/nextflow-run/outputs

      first_output:
        - batch/intermediates/Bonjour-output.txt
        - batch/intermediates/Hello-output.txt
        - batch/intermediates/Hola-output.txt

      uppercased:
        - batch/intermediates/UPPER-Bonjour-output.txt
        - batch/intermediates/UPPER-Hello-output.txt
        - batch/intermediates/UPPER-Hola-output.txt

      collected: batch/intermediates/COLLECTED-batch-output.txt

      batch_report: batch/batch-report.txt

      cowpy_art: batch/cowpy-COLLECTED-batch-output.txt
    ```

??? abstract "Contenu du répertoire"

    ```console
    outputs/batch
    ├── batch-report.txt
    ├── cowpy-COLLECTED-batch-output.txt
    └── intermediates
        ├── Bonjour-output.txt
        ├── COLLECTED-batch-output.txt
        ├── Hello-output.txt
        ├── Hola-output.txt
        ├── UPPER-Bonjour-output.txt
        ├── UPPER-Hello-output.txt
        └── UPPER-Hola-output.txt
    ```

Les sorties se trouvent maintenant dans `outputs/batch/` au lieu du répertoire par défaut `results/batch/`.
Le code du pipeline lui-même décide toujours de la structure à l'intérieur de ce répertoire de base, comme les sous-répertoires `batch/` et `intermediates/` ; `-output-dir` contrôle uniquement l'endroit où cette structure commence.

`-output-dir` n'est en réalité qu'un raccourci en ligne de commande pour l'option de configuration `outputDir`, elle peut donc être placée n'importe où où une configuration peut l'être : directement dans `nextflow.config`, à l'intérieur d'un profil, ou dans un fichier de superposition `-c` comme celui que vous avez utilisé plus tôt dans cette partie.
Par exemple, cet extrait montre le même paramètre placé directement dans `nextflow.config` au lieu d'être passé en ligne de commande :

```groovy title="nextflow.config"
outputDir = 'outputs'
```

Consultez [Configuration file](https://nextflow.io/docs/latest/config.html) dans la référence Nextflow pour la liste complète des endroits où une option de configuration comme celle-ci peut être placée.

### 3.2. Choisir comment les sorties sont publiées

Par défaut, Nextflow publie les sorties sous forme de liens symboliques pointant vers les emplacements des sorties dans `work/`, et non de vraies copies :

```console
$ ls -l results/batch/intermediates/Hello-output.txt
lrwxr-xr-x  ...  Hello-output.txt -> /workspaces/training/nextflow-run/work/b7/b8c4d1.../Hello-output.txt
```

Les auteur·es de pipelines peuvent définir le « mode de publication » sur `'copy'` ou `'move'` pour chaque processus individuel dans le code du workflow.
Ils·elles le font généralement pour les sorties finales du pipeline, tout en laissant le comportement par défaut `'symlink'` pour les fichiers intermédiaires qui peuvent être supprimés une fois que le pipeline complet a été exécuté.

Cela évite de dupliquer les données sur le disque, mais cela signifie que vous ne pouvez pas supprimer les répertoires de tâches dans `work/` sans rompre le lien, perdant ainsi la possibilité d'utiliser `-resume`.
Si vous souhaitez que tous les fichiers de sortie soient correctement copiés, définissez [`workflow.output.mode`](https://nextflow.io/docs/latest/reference/config.html#workflow) sur `'copy'` dans votre configuration de pipeline. (Contrairement à `-output-dir`, il n'y a pas d'option en ligne de commande pour cela ; c'est uniquement via la configuration.)

Essayez de le définir dans `nextflow.config` :

=== "Après"

    ```groovy title="nextflow.config" hl_lines="7"
    params {
        input = 'data/greetings.csv'
        batch = 'batch'
        character = 'turkey'
    }

    workflow.output.mode = 'copy'
    ```

=== "Avant"

    ```groovy title="nextflow.config"
    params {
        input = 'data/greetings.csv'
        batch = 'batch'
        character = 'turkey'
    }
    ```

Puis exécutez le pipeline en changeant le nom du batch afin de pouvoir voir la différence dans les sorties :

```bash
nextflow run main.nf --batch withmode
```

??? success "Sortie de la commande"

    ```console hl_lines="6 7 8 9"
    N E X T F L O W   ~  version 26.04.4

    Launching `main.nf` [angry_noether] revision: ce74f81996

    executor >  local (8)
    [41/2b478d] sayHello (3)       | 3 of 3 ✔
    [bf/dd2840] convertToUpper (2) | 3 of 3 ✔
    [ea/364e97] collectGreetings   | 1 of 1 ✔
    [10/76fe7b] cowpy              | 1 of 1 ✔

    Outputs:

      /workspaces/training/nextflow-run/results

      first_output:
        - withmode/intermediates/Bonjour-output.txt
        - withmode/intermediates/Hello-output.txt
        - withmode/intermediates/Hola-output.txt

      uppercased:
        - withmode/intermediates/UPPER-Bonjour-output.txt
        - withmode/intermediates/UPPER-Hello-output.txt
        - withmode/intermediates/UPPER-Hola-output.txt

      collected: withmode/intermediates/COLLECTED-withmode-output.txt

      batch_report: withmode/withmode-report.txt

      cowpy_art: withmode/cowpy-COLLECTED-withmode-output.txt
    ```

Examinez l'un des fichiers de sortie comme précédemment :

```console
$ ls -l results/withmode/intermediates/Hello-output.txt
-rw-r--r--  ...  Hello-output.txt
```

C'est maintenant un vrai fichier indépendant qui restera disponible même si `work/` est nettoyé.

!!! warning "Avertissement"

    Le paramètre `workflow.output.mode` ne fait que définir une valeur par défaut pour les sorties qui n'ont pas déjà un mode défini dans le code du pipeline.
    Il ne peut pas remplacer un mode codé en dur par l'auteur·e, quelle que soit la valeur que vous lui attribuez.

### À retenir

Vous savez comment personnaliser le répertoire de sortie de base et choisir entre des sorties copiées et des liens symboliques, le tout sans toucher au code du pipeline.

### Et ensuite ?

Passez à la [Partie 3](./03_manage_executions.md), où vous apprendrez à inspecter l'historique des exécutions passées, à générer des rapports d'exécution et à nettoyer les anciens répertoires de travail.

---

## Résumé

Dans cette partie, vous avez appris à :

- Configurer le comportement du pipeline à l'aide de `nextflow.config` et des profils
- Fournir une configuration via un fichier de configuration spécifique à une exécution ou un fichier de paramètres
- Personnaliser le répertoire de sortie et choisir entre des sorties copiées et des liens symboliques
