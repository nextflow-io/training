# Partie 2 : Gérer les ressources de calcul et les échecs

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Traduction assistée par IA - [en savoir plus et suggérer des améliorations](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

Dans la [Partie 1](./01_packaging_and_execution.md), vous avez adapté où et comment les tâches d'un pipeline s'exécutent.
Ici, vous allez adapter la quantité de ressources de calcul allouée à chaque tâche, et ce qui se passe lorsqu'une tâche échoue malgré votre meilleure estimation d'allocation.

---

## 1. Contrôler les allocations de ressources de calcul

Par défaut, Nextflow alloue un seul CPU à chaque processus via la directive `cpus`, et n'impose pas de limite mémoire sauf si vous en définissez une :

```groovy title="Built-in configuration"
process {
    cpus = 1
}
```

Vous savez déjà grâce à [Nextflow Run](../nextflow_run/index.md) que la configuration de ce pipeline définit `memory` à 1 Go pour tous les processus.
Mais comment savoir quelles valeurs utiliser réellement pour vos propres pipelines ?

### 1.1. Générer un rapport d'utilisation des ressources

Vous avez déjà généré un rapport d'exécution avec `-with-report` dans [Nextflow Run](../nextflow_run/02_configure_pipeline.md).
Ce même rapport vous permet de savoir de combien de CPU et de mémoire vos processus ont réellement besoin : exécutez le workflow avec des allocations par défaut, enregistrez l'utilisation réelle, puis ajustez en conséquence.

```bash
nextflow run main.nf -with-report report-config-1.html
```

Le rapport est un fichier HTML que vous pouvez ouvrir dans un navigateur.
Il détaille le temps d'exécution et l'utilisation des ressources par processus, notamment le pourcentage des ressources allouées qui a été réellement utilisé.
Voici ce qu'il affiche pour `cowpy` avec les valeurs par défaut actuelles (1 CPU, 1 Go de mémoire) :

| Métrique              | Valeur |
| --------------------- | ------ |
| Utilisation CPU       | 116%   |
| Mémoire maximale utilisée | 6,4 Mo |
| Mémoire allouée       | 1 Go   |

`cowpy` utilise bien moins de 1% de son allocation de 1 Go ; le `%cpu` supérieur à 100% signifie simplement qu'il utilise brièvement plus d'un CPU de traitement à l'intérieur du conteneur, par courtes rafales.

Consultez [Reports](https://nextflow.io/docs/latest/reports.html) pour la liste complète des fonctionnalités disponibles.

### 1.2. Définir les allocations de ressources pour un processus spécifique

Le rapport ci-dessus montre que `cowpy` est confortablement dans les limites de son allocation actuelle, mais supposons que vous souhaitiez lui donner plus de marge, par exemple parce que vous attendez des entrées plus volumineuses en production.
Vous pouvez remplacer les valeurs par défaut pour un seul processus avec `withName`.

=== "Après"

    ```groovy title="nextflow.config" linenums="6" hl_lines="5-6"
    process {
        memory = 1.GB
        withName: 'cowpy' {
            conda = 'conda-forge::cowpy==1.1.5'
            memory = 2.GB
            cpus = 2
        }
    }
    ```

=== "Avant"

    ```groovy title="nextflow.config" linenums="6"
    process {
        memory = 1.GB
        withName: 'cowpy' {
            conda = 'conda-forge::cowpy==1.1.5'
        }
    }
    ```

Avec cette configuration, chaque processus demande 1 Go de mémoire et un seul CPU, sauf `cowpy`, qui demande 2 Go et 2 CPUs (en plus du paramètre `conda` de la [Partie 1](./01_packaging_and_execution.md)).

!!! info "Info"

    Si votre machine dispose de peu de CPUs et que vous en allouez un grand nombre par processus, les appels de tâches peuvent se mettre en file d'attente les uns derrière les autres, car Nextflow ne demandera pas plus de CPUs que ce qui est disponible.

Exécutez-le à nouveau avec un nom de rapport différent, afin de pouvoir comparer avant et après.

```bash
nextflow run main.nf -with-report report-config-2.html
```

??? success "Sortie de la commande"

    ```console
    N E X T F L O W   ~  version 26.04.4

    Launching `main.nf` [voluminous_venter] revision: c3c85dec78

    executor >  local (8)
    [a1/0e96d4] sayHello (1)       | 3 of 3 ✔
    [a3/7173a3] convertToUpper (2) | 3 of 3 ✔
    [4f/a8ae3d] collectGreetings   | 1 of 1 ✔
    [91/3724f8] cowpy              | 1 of 1 ✔

    Outputs:

      /workspaces/training/execution-config/results

      first_output:
        - full_pipeline/intermediates/Bonjour-output.txt
        - full_pipeline/intermediates/Hola-output.txt
        - full_pipeline/intermediates/Hello-output.txt

      uppercased:
        - full_pipeline/intermediates/UPPER-Bonjour-output.txt
        - full_pipeline/intermediates/UPPER-Hello-output.txt
        - full_pipeline/intermediates/UPPER-Hola-output.txt

      collected: full_pipeline/intermediates/COLLECTED-batch-output.txt

      batch_report: full_pipeline/batch-report.txt

      cowpy_art: full_pipeline/cowpy-COLLECTED-batch-output.txt
    ```

Comparaison des deux rapports pour `cowpy` :

| Métrique              | Avant (1 CPU, 1 Go) | Après (2 CPUs, 2 Go) |
| --------------------- | ------------------- | -------------------- |
| Mémoire maximale utilisée | 6,4 Mo          | 6,4 Mo               |
| Utilisation CPU       | 116%                | 118%                 |

Doubler l'allocation n'a pas du tout modifié l'utilisation réelle, ce qui indique que le 1 Go / 1 CPU d'origine était déjà généreux pour cette charge de travail de démonstration.
Sur un pipeline réel traitant des données non triviales, vous vous attendriez à ce que les chiffres diffèrent significativement entre les processus, ce qui explique précisément pourquoi vous profilez avant de décider des allocations, plutôt que de les deviner.

### 1.3. Ajouter des limites de ressources

Selon votre infrastructure de calcul, il peut exister des contraintes strictes sur ce que vous pouvez demander, par exemple un plafond à l'échelle du cluster.
La directive `resourceLimits` vous permet de définir ces limites :

```groovy title="Syntax example"
process {
    resourceLimits = [
        memory: 750.GB,
        cpus: 200,
        time: 30.d
    ]
}
```

Nextflow traduit ces limites en ce que l'exécuteur cible attend.
Si un processus demande plus que la limite, la demande est plafonnée plutôt que rejetée.

!!! warning "Avertissement"

    Ceci ne peut pas être exécuté dans l'environnement de formation, car cela nécessite une infrastructure HPC pour avoir un effet.

??? info "Configurations de référence institutionnelles"

    Le projet nf-core maintient une [collection de fichiers de configuration](https://nf-co.re/configs/) partagés par des institutions du monde entier, couvrant un large éventail d'exécuteurs HPC et cloud.
    Ils constituent un bon point de départ, que votre propre institution en fasse partie ou non.

### À retenir

Vous savez comment générer un rapport de profilage pour évaluer l'utilisation des ressources, remplacer les allocations de ressources pour un processus spécifique, et plafonner les allocations avec `resourceLimits`.

### Et ensuite ?

Apprenez à faire en sorte qu'un pipeline se rétablisse automatiquement lorsqu'une tâche échoue, que votre estimation d'allocation de ressources soit correcte ou non.

---

## 2. Gérer les échecs de tâches avec les nouvelles tentatives

Le profilage vous indique ce dont un processus a besoin la plupart du temps, mais les charges de travail réelles varient : une allocation confortable pour la plupart des entrées peut encore être trop juste pour une entrée inhabituellement volumineuse, et les estimations peuvent simplement être erronées.
Plutôt que de laisser une seule tâche en échec interrompre toute l'exécution, Nextflow peut relancer automatiquement une tâche en échec, en lui allouant optionnellement plus de ressources à chaque tentative.

### 2.1. Relancer automatiquement une tâche en échec

Pour voir cela en action, définissez délibérément l'allocation mémoire de `cowpy` en dessous de ce dont il a réellement besoin : rappelons depuis la [section 1.1](#11-generate-a-resource-utilization-report) qu'il atteint un pic d'environ 6,4 Mo, donc 6 Mo devrait être juste insuffisant.

=== "Après"

    ```groovy title="nextflow.config" linenums="6" hl_lines="5-7"
    process {
        memory = 1.GB
        withName: 'cowpy' {
            conda = 'conda-forge::cowpy==1.1.5'
            memory = 6.MB
            errorStrategy = 'retry'
            maxRetries = 2
        }
    }
    ```

=== "Avant"

    ```groovy title="nextflow.config" linenums="6"
    process {
        memory = 1.GB
        withName: 'cowpy' {
            conda = 'conda-forge::cowpy==1.1.5'
            memory = 2.GB
            cpus = 2
        }
    }
    ```

`errorStrategy` indique à Nextflow quoi faire lorsqu'une tâche échoue : `'retry'` resoumet la tâche au lieu d'arrêter tout le pipeline.
`maxRetries` limite le nombre de tentatives supplémentaires avant que Nextflow abandonne.

```bash
nextflow run main.nf
```

??? failure "Sortie de la commande (abrégée)"

    ```console
    N E X T F L O W   ~  version 26.04.4

    Launching `main.nf` [desperate_brazil] revision: c3c85dec78

    executor >  local (10)
    [67/fe1f49] sayHello (1)       | 3 of 3 ✔
    [8a/f13335] convertToUpper (1) | 3 of 3 ✔
    [39/1b24ed] collectGreetings   | 1 of 1 ✔
    [7a/d5eb6f] cowpy              | 0 of 1, retries: 2 ✘
    [9d/b79eb3] NOTE: Process `cowpy` terminated with an error exit status (137) -- Execution is retried (1)
    [6d/1d9d84] NOTE: Process `cowpy` terminated with an error exit status (137) -- Execution is retried (2)
    ERROR ~ Error executing process > 'cowpy'

    Caused by:
      Process `cowpy` terminated with an error exit status (137)

    Command executed:
      cat COLLECTED-batch-output.txt | cowpy -c "turkey" > cowpy-COLLECTED-batch-output.txt

    Command exit status:
      137

    Command output:
      (empty)

    Command error:
      /usr/local/bin/_activate_current_env.sh: line 35:    14 Killed                  micromamba activate "${ENV_NAME:-base}"

    Work dir:
      /workspaces/training/execution-config/work/7a/d5eb6feeac0eed18d95d3da7a7aeb4

    Tip: when you have fixed the problem you can continue the execution adding the option `-resume` to the run command line

    -- Check '.nextflow.log' file for details
    ```

Le code de sortie 137 est le signal standard pour un arrêt par manque de mémoire : le conteneur ne disposait pas de suffisamment de mémoire pour exécuter `cowpy`. Nextflow a relancé la tâche deux fois, soit trois tentatives au total, conformément à `maxRetries = 2`.
Comme l'allocation mémoire n'a jamais changé entre les tentatives, chaque tentative s'est heurtée au même obstacle ; une fois les nouvelles tentatives épuisées, Nextflow signale l'échec en détail et arrête le pipeline en quittant avec un statut non nul.

Relancer une tâche seul ne résout rien si la cause sous-jacente ne change pas entre les tentatives.

### 2.2. Augmenter les ressources à chaque nouvelle tentative

Dans une directive de processus, `task.attempt` contient le numéro de la tentative en cours, en commençant à 1.
Vous pouvez l'utiliser dans une closure pour augmenter une allocation de ressources à chaque nouvelle tentative.

=== "Après"

    ```groovy title="nextflow.config" linenums="6" hl_lines="5 7"
    process {
        memory = 1.GB
        withName: 'cowpy' {
            conda = 'conda-forge::cowpy==1.1.5'
            memory = { 6.MB * task.attempt }
            errorStrategy = 'retry'
            maxRetries = 3
        }
    }
    ```

=== "Avant"

    ```groovy title="nextflow.config" linenums="6"
    process {
        memory = 1.GB
        withName: 'cowpy' {
            conda = 'conda-forge::cowpy==1.1.5'
            memory = 6.MB
            errorStrategy = 'retry'
            maxRetries = 2
        }
    }
    ```

Exécutez à nouveau le workflow :

```bash
nextflow run main.nf
```

??? success "Sortie de la commande (abrégée)"

    ```console
    N E X T F L O W   ~  version 26.04.4

    Launching `main.nf` [grave_joliot] revision: c3c85dec78

    executor >  local (9)
    [0f/211b8a] sayHello (2)       | 3 of 3 ✔
    [26/301ea2] convertToUpper (3) | 3 of 3 ✔
    [22/5a895b] collectGreetings   | 1 of 1 ✔
    [e1/beee86] cowpy              | 1 of 1, retries: 1 ✔
    [b6/7aed6a] NOTE: Process `cowpy` terminated with an error exit status (137) -- Execution is retried (1)

    Outputs:

      /workspaces/training/execution-config/results

      first_output:
        - full_pipeline/intermediates/Hola-output.txt
        - full_pipeline/intermediates/Hello-output.txt
        - full_pipeline/intermediates/Bonjour-output.txt

      uppercased:
        - full_pipeline/intermediates/UPPER-Hola-output.txt
        - full_pipeline/intermediates/UPPER-Hello-output.txt
        - full_pipeline/intermediates/UPPER-Bonjour-output.txt

      collected: full_pipeline/intermediates/COLLECTED-batch-output.txt

      batch_report: full_pipeline/batch-report.txt

      cowpy_art: full_pipeline/cowpy-COLLECTED-batch-output.txt
    ```

La première tentative échoue toujours à 6 Mo, mais la nouvelle tentative s'exécute avec 12 Mo (`6.MB * 2`) et réussit, et le pipeline se termine avec toutes les sorties publiées.

!!! warning "Avertissement"

    La sortie console inclut toujours une ligne `NOTE:` signalant la première tentative en échec, même si le pipeline dans son ensemble a réussi : Nextflow enregistre chaque nouvelle tentative individuellement, mais un échec suivi d'une nouvelle tentative n'affecte pas le résultat global.
    Vérifiez le résumé `Outputs:`, ou le statut de sortie de la commande, pour confirmer si l'exécution a réellement réussi.

Consultez [Dynamic computing resources](https://nextflow.io/docs/latest/process.html#dynamic-task-resources) dans la documentation Nextflow pour des modèles de nouvelles tentatives plus avancés, notamment la mise à l'échelle basée sur l'erreur spécifique survenue.

### À retenir

Vous savez comment faire en sorte qu'un pipeline relance automatiquement les tâches en échec, et comment augmenter les allocations de ressources à chaque nouvelle tentative en utilisant `task.attempt`.

### Et ensuite ?

Passez à la [Partie 3](./03_profiles.md), où vous apprendrez à regrouper ce type de configuration dans des profils commutables.

---

## Résumé

Dans cette partie, vous avez appris à :

- Générer un rapport de profilage des ressources et définir des allocations de ressources par processus
- Plafonner les demandes de ressources avec `resourceLimits`
- Relancer automatiquement une tâche en échec avec `errorStrategy` et `maxRetries`
- Augmenter une allocation de ressources à chaque nouvelle tentative en utilisant `task.attempt`
