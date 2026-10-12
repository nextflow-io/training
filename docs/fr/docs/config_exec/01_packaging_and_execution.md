# Partie 1 : Adapter à votre environnement de calcul

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Traduction assistée par IA - [en savoir plus et suggérer des améliorations](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

Dans [Nextflow Run](../nextflow_run/index.md), vous avez configuré les entrées, les paramètres et les sorties d'un pipeline.
Ce cours couvre l'autre aspect de la question : adapter l'exécution d'un pipeline à n'importe quel environnement de calcul sur lequel il est amené à s'exécuter, sans modifier le code du workflow.

!!! example "Scénario"

    Vous avez développé et testé votre pipeline sur votre ordinateur portable en utilisant Docker.
    Vous devez maintenant le transmettre : un·e collaborateur·trice n'a que Conda de configuré, et le cluster HPC de votre institution attend que les jobs passent par son propre ordonnanceur avec ses propres limites de ressources.
    Rien de tout cela ne devrait nécessiter de réécrire le pipeline lui-même.

Le même code de pipeline peut s'exécuter dans tous ces environnements, car rien de tout cela n'est intégré dans le workflow.
L'empaquetage logiciel, la plateforme d'exécution et l'allocation des ressources sont tous contrôlés via la configuration, superposée au code — et c'est ce que couvre ce cours : comment adapter le même pipeline à un nouvel environnement en modifiant la configuration, pas le code.

---

## 1. Sélectionner une technologie d'empaquetage logiciel

Dans [Nextflow Run](../nextflow_run/index.md), vous avez vu un profil `conda` déjà configuré dans `nextflow.config` comme alternative à Docker.
Ici, vous allez construire ce même mécanisme de bascule vous-même, et voir ce qu'il faut pour rendre un processus réellement utilisable avec Conda.

### 1.1. Désactiver Docker et activer Conda

Passez `docker.enabled` à `false` et ajoutez une directive activant Conda.

=== "Après"

    ```groovy title="nextflow.config" linenums="1" hl_lines="1-2"
    docker.enabled = false
    conda.enabled = true
    ```

=== "Avant"

    ```groovy title="nextflow.config" linenums="1"
    docker.enabled = true
    ```

Cela permet à Nextflow de créer et d'utiliser des environnements Conda pour tout processus disposant d'un paquet Conda spécifié.
Le processus `cowpy` n'en a pas encore, alors ajoutons-en un, entièrement depuis la configuration.

### 1.2. Ajouter un paquet Conda via la configuration

Une directive `conda` peut être définie dans la définition du processus elle-même, de la même façon que `container` l'est déjà dans `modules/cowpy.nf`, mais ce n'est pas obligatoire : `withName` vous permet de la définir depuis la configuration, limitée au seul processus `cowpy`.

=== "Après"

    ```groovy title="nextflow.config" linenums="6" hl_lines="3-5"
    process {
        memory = 1.GB
        withName: 'cowpy' {
            conda = 'conda-forge::cowpy==1.1.5'
        }
    }
    ```

=== "Avant"

    ```groovy title="nextflow.config" linenums="6"
    process {
        memory = 1.GB
    }
    ```

Cela ne remplace pas la directive `container` déjà présente dans le code du pipeline, cela ajoute une alternative à côté d'elle, sans toucher du tout à ce code.

!!! tip "Astuce"

    La recherche [Seqera Containers](https://seqera.io/containers/) est un moyen pratique de rechercher l'URI du paquet Conda pour un outil donné, même si vous ne prévoyez pas de construire un conteneur à partir de celui-ci.

### 1.3. Exécuter le workflow pour vérifier qu'il peut utiliser Conda

```bash
nextflow run main.nf --batch conda
```

??? success "Sortie de la commande"

    ```console
    N E X T F L O W   ~  version 26.04.4

    Launching `main.nf` [romantic_pike] revision: c3c85dec78

    executor >  local (8)
    [6d/d48030] sayHello (2)       | 3 of 3 ✔
    [f5/7a9d76] convertToUpper (1) | 3 of 3 ✔
    [1c/79b693] collectGreetings   | 1 of 1 ✔
    Creating env using conda: conda-forge::cowpy==1.1.5 [cache /workspaces/training/config-exec/work/conda/env-898314d566668b6587ad714ae06b8520]
    [bb/64b67c] cowpy              | 1 of 1 ✔

    Outputs:

      /workspaces/training/config-exec/results

      first_output:
        - full_pipeline/intermediates/Hola-output.txt
        - full_pipeline/intermediates/Bonjour-output.txt
        - full_pipeline/intermediates/Hello-output.txt

      uppercased:
        - full_pipeline/intermediates/UPPER-Hello-output.txt
        - full_pipeline/intermediates/UPPER-Bonjour-output.txt
        - full_pipeline/intermediates/UPPER-Hola-output.txt

      collected: full_pipeline/intermediates/COLLECTED-conda-output.txt

      batch_report: full_pipeline/conda-report.txt

      cowpy_art: full_pipeline/cowpy-COLLECTED-conda-output.txt
    ```

Cela produit la même sortie qu'avec Docker, même si les mécanismes sont différents en coulisses : Nextflow récupère le paquet Conda et construit un environnement à partir de celui-ci, au lieu de télécharger une image de conteneur.

!!! info "Info"

    La construction d'un nouvel environnement Conda peut prendre un peu plus de temps que le téléchargement d'un conteneur la première fois, mais le paquet utilisé ici est petit, donc cela devrait être rapide.

Repassez maintenant à Docker pour le reste de ce cours.

```groovy title="nextflow.config" linenums="1"
docker.enabled = true
```

??? tip "Mélanger Docker et Conda"

    Ces paramètres étant limités par processus, vous pouvez les combiner : certains processus utilisent Docker, d'autres utilisent Conda, selon ce qui est disponible pour chaque outil.
    Si une directive `container` (dans le code du pipeline) et une directive `conda` (ici, depuis la configuration) sont toutes deux définies pour le même processus et que les deux systèmes d'empaquetage sont activés, Nextflow donne la priorité aux conteneurs.

### À retenir

Vous savez comment configurer la technologie d'empaquetage logiciel qu'un processus doit utiliser, et comment basculer entre Docker et Conda.

### Et ensuite ?

Apprenez comment modifier la plateforme d'exécution que Nextflow utilise pour exécuter réellement vos tâches.

---

## 2. Sélectionner une plateforme d'exécution

Chaque pipeline que vous avez exécuté jusqu'à présent a utilisé l'executor local : chaque tâche s'exécute sur la même machine que Nextflow lui-même.
Nextflow vérifie les CPU et la mémoire disponibles, et retient les tâches jusqu'à ce que suffisamment de ressources se libèrent.

L'executor local est pratique, mais il ne passe pas à l'échelle au-delà d'une seule machine.
Nextflow prend en charge [de nombreux autres backends d'exécution](https://nextflow.io/docs/latest/executor.html), notamment les ordonnanceurs HPC (Slurm, LSF, SGE, PBS, et d'autres) et les plateformes cloud (AWS Batch, Google Cloud Batch, Azure Batch, Kubernetes, et plus encore).

### 2.1. Cibler un backend différent

L'executor est défini par une directive de processus appelée `executor`.
Par défaut, il est `local`, donc ce qui suit est implicite :

```groovy title="Built-in configuration"
process {
    executor = 'local'
}
```

Pour cibler un backend différent, définissez la directive sur l'executor souhaité.

```groovy title="nextflow.config"
process {
    executor = 'slurm'
}
```

!!! warning "Avertissement"

    L'environnement de formation n'est pas connecté à un cluster HPC, donc vous ne pouvez pas exécuter cela ici.

### 2.2. La syntaxe spécifique au backend est abstraite

La plupart des plateformes HPC exigent que les soumissions de jobs spécifient des demandes de ressources, telles que les CPU, la mémoire et le nom de la file d'attente, en utilisant leur propre syntaxe.
La même demande de 8 CPU et 4 Go de RAM sur une file d'attente appelée `my-science-work` est complètement différente selon l'ordonnanceur.

??? abstract "Exemples"

    ```bash title="Config for SLURM / submit using sbatch"
    #SBATCH -o /path/to/my/task/directory/my-task-1.log
    #SBATCH --no-requeue
    #SBATCH -c 8
    #SBATCH --mem 4096M
    #SBATCH -p my-science-work
    ```

    ```bash title="Config for PBS / submit using qsub"
    #PBS -o /path/to/my/task/directory/my-task-1.log
    #PBS -j oe
    #PBS -q my-science-work
    #PBS -l nodes=1:ppn=8
    #PBS -l mem=4gb
    ```

    ```bash title="Config for SGE / submit using qsub"
    #$ -o /path/to/my/task/directory/my-task-1.log
    #$ -j y
    #$ -terse
    #$ -notify
    #$ -q my-science-work
    #$ -l slots=8
    #$ -l h_rss=4096M,mem_free=4096M
    ```

Nextflow abstrait tout cela : vous spécifiez des propriétés standardisées telles que `cpus`, `memory` et `queue` une seule fois (voir les [directives de processus](https://nextflow.io/docs/latest/reference/process.html#process-directives) pour la liste complète), et Nextflow les traduit en scripts spécifiques au backend approprié au moment de l'exécution.

### 2.3. Voir ce que Nextflow exécute réellement

Cette traduction n'est pas seulement une commodité de fichier de configuration : elle repose sur quelque chose de concret que vous pouvez inspecter dès maintenant, même avec l'executor local.
Dans [Nextflow Run, section 1.3](../nextflow_run/01_run_nextflow.md#13-explore-the-work-directory), vous avez regardé à l'intérieur d'un répertoire de tâche sous `work/` et trouvé `.command.sh`, la commande exacte exécutée par Nextflow.
Ce même répertoire contient également un fichier que vous n'avez pas encore examiné : `.command.run`.

```bash
cat work/0a/0df4a1*/.command.run
```

??? success "Résultat de la commande (extrait)"

    ```console
    #!/bin/bash
    ### ---
    ### name: 'convertToUpper (3)'
    ### container: 'null'
    ### outputs:
    ### - 'UPPER-Bonjour-output.txt'
    ### ...
    set -e
    set -u
    ...
    nxf_launch() {
        /bin/bash -ue /workspaces/training/nextflow-run/work/0a/0df4a1028c2001758b1841cff92fc7/.command.sh
    }
    ...
    ```

`.command.run` est le vrai script que Nextflow transmet pour exécution.
Il enveloppe `.command.sh` avec tout ce qui est nécessaire pour l'exécuter réellement : configuration de l'environnement, staging des entrées/sorties, et rapport du résultat à Nextflow.
Avec l'executor `local`, Nextflow exécute simplement ce script sur la même machine.

C'est exactement ce qui change lorsque vous définissez un `executor` différent.
Pour un ordonnanceur HPC tel que Slurm ou PBS, Nextflow génère le même type de script d'enveloppe, ajoute l'en-tête spécifique à l'ordonnanceur que vous avez vu dans la [section 2.2](#22-backend-specific-syntax-is-abstracted-away) (traduit depuis vos paramètres `cpus`, `memory` et `queue`), et transmet le résultat à la commande de soumission propre à cet ordonnanceur, par exemple `sbatch` pour Slurm.
À partir de là, Nextflow interroge l'ordonnanceur pour connaître l'état du job au lieu de surveiller directement un processus local.
Les backends cloud batch fonctionnent un peu différemment, car ils sont pilotés par des appels API plutôt que par une commande de soumission, mais la même idée sous-jacente s'applique : le même script de tâche s'exécute, seuls la façon dont il est lancé et suivi changent.

### À retenir

Vous savez comment modifier l'executor pour cibler différentes infrastructures de calcul, que Nextflow abstrait la syntaxe de soumission spécifique au backend, et ce qui se passe réellement en coulisses lorsqu'une tâche s'exécute sur un backend différent.

### Et ensuite ?

Passez à la [Partie 2](./02_resources_and_retries.md), où vous apprendrez à profiler et allouer des ressources de calcul, et à gérer les échecs de tâches avec des nouvelles tentatives.

---

## Résumé

Dans cette partie, vous avez appris à :

- Basculer la technologie d'empaquetage logiciel entre Docker et Conda
- Ajouter une directive `conda` à une définition de processus
- Modifier la plateforme d'exécution avec la directive `executor`
- Inspecter ce que Nextflow génère et exécute réellement pour une tâche, et comment cela change selon les executors
