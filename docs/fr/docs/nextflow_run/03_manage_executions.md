# Partie 3 : Gérer les exécutions de workflow

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Traduction assistée par IA - [en savoir plus et suggérer des améliorations](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

Au fur et à mesure que vous exécutez et réexécutez des pipelines, vous accumulez un historique d'exécutions et d'anciens répertoires `work/`.
Dans la [Partie 1](./01_run_nextflow.md#23-use-resume-to-skip-completed-work), vous avez déjà utilisé `-resume` pour ignorer le travail déjà effectué.
Ici, vous apprendrez à générer des rapports sur une exécution, à inspecter l'historique des exécutions passées avec [`nextflow log`](https://nextflow.io/docs/latest/reference/cli.html#log), et à supprimer les anciens répertoires de travail dont vous n'avez plus besoin avec [`nextflow clean`](https://nextflow.io/docs/latest/reference/cli.html#clean).

---

## 1. Générer des rapports de pipeline

Nextflow peut générer plusieurs types de rapports sur une exécution, chacun ajouté avec son propre flag `-with-*` : un rapport d'exécution (`-with-report`), une chronologie d'exécution (`-with-timeline`), un fichier de trace des tâches (`-with-trace`), et un diagramme de workflow (`-with-dag`).
Nous allons générer les deux premiers ici ; consultez [Execution reports](https://nextflow.io/docs/latest/reports.html) dans la référence Nextflow pour les autres.

### 1.1. Générer un rapport d'exécution

Ajoutez `-with-report` à n'importe quelle commande `nextflow run` pour générer un rapport HTML une fois le pipeline terminé :

```bash
nextflow run main.nf --input data/greetings.csv --character turkey -with-report
```

??? success "Sortie de la commande"

    ```console hl_lines="6 7 8 9"
    N E X T F L O W   ~  version 26.04.4

    Launching `main.nf` [intergalactic_dalembert] revision: ce74f81996

    executor >  local (8)
    [34/23f10f] sayHello (2)       | 3 of 3 ✔
    [af/cab69d] convertToUpper (3) | 3 of 3 ✔
    [9e/d73afb] collectGreetings   | 1 of 1 ✔
    [3c/392db0] cowpy              | 1 of 1 ✔

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

Nextflow écrit le rapport dans un fichier nommé `report-<timestamp>.html` dans le répertoire de travail.
Ouvrez-le dans un navigateur pour voir un résumé de l'exécution, un tableau de chaque tâche avec son statut et sa durée d'exécution, ainsi que des graphiques d'utilisation des ressources ventilés par processus.

L'onglet **Tasks** liste chaque tâche exécutée par le pipeline, avec le nom de son processus, son statut et son utilisation des ressources :

![Tableau des tâches du rapport d'exécution](img/execution_report_tasks.png)

Le rapport est particulièrement utile lorsqu'un pipeline prend plus de temps que prévu ou qu'une tâche échoue : le tableau des tâches indique exactement où le temps a été passé et quelles tâches ont réussi ou échoué.

### 1.2. Générer une chronologie d'exécution

Ajoutez `-with-timeline` à une exécution pour obtenir une vue de type diagramme de Gantt indiquant quand chaque tâche s'est exécutée :

```bash
nextflow run main.nf --input data/greetings.csv --character turkey -with-timeline
```

??? success "Sortie de la commande"

    ```console hl_lines="6 7 8 9"
    N E X T F L O W   ~  version 26.04.4

    Launching `main.nf` [jolly_noyce] revision: ce74f81996

    executor >  local (8)
    [ad/e92ef3] sayHello (3)       | 3 of 3 ✔
    [2a/df8a8d] convertToUpper (2) | 3 of 3 ✔
    [be/7fb72a] collectGreetings   | 1 of 1 ✔
    [63/dc9bd6] cowpy              | 1 of 1 ✔

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

Nextflow écrit la chronologie dans un fichier nommé `timeline-<timestamp>.html`.
Ouvrez-le dans un navigateur pour voir une barre pour chaque tâche, positionnée et dimensionnée selon le moment où elle s'est exécutée et sa durée :

![Chronologie d'exécution](img/execution_timeline.png)

La chronologie rend visible en un coup d'œil la forme d'expansion puis de convergence de la [Partie 1](./01_run_nextflow.md#31-run-the-workflow) : les trois tâches `sayHello` s'exécutent en parallèle, puis les trois tâches `convertToUpper`, puis `collectGreetings` et `cowpy` s'exécutent l'une après l'autre car chacune dépend de tout ce qui précède.

### À retenir

Vous savez comment générer un rapport d'exécution HTML avec `-with-report` et une chronologie d'exécution avec `-with-timeline`, et où chercher les autres types de rapports pris en charge par Nextflow.

### Et ensuite ?

Apprenez à inspecter l'historique des exécutions passées.

---

## 2. Inspecter le journal des exécutions passées

Que vous développiez un pipeline ou que vous l'exécutiez en production, vous aurez à un moment ou un autre besoin de consulter des informations sur des exécutions passées.

### 2.1. Le fichier d'historique

Chaque fois que vous lancez un workflow Nextflow, une ligne est écrite dans un fichier journal appelé `history`, dans un répertoire caché appelé `.nextflow` dans le répertoire de travail courant.

??? abstract "Contenu du fichier"

    ```txt title=".nextflow/history" linenums="1"
    2026-09-13 01:47:58	2.8s	determined_ramanujan	OK	ce74f81996b8ce999853fdbb25ba7969	03078950-5011-414f-992b-8fe1bd219ac2	nextflow run main.nf --input data/greetings.csv --character turkey
    2026-09-13 01:48:03	3s	prickly_cuvier	OK	ce74f81996b8ce999853fdbb25ba7969	5a9e8bf7-a7db-4c2c-b1d8-3bd4c7266528	nextflow run main.nf --input data/greetings.csv --character tux
    2026-09-13 01:48:09	2.3s	elegant_panini	OK	ce74f81996b8ce999853fdbb25ba7969	517ccc29-9db7-4d8d-8ea0-c8db6e9ac653	nextflow run main.nf --input data/greetings.csv --character stegosaurus
    2026-09-13 01:48:14	1.5s	deadly_lamport	OK	ce74f81996b8ce999853fdbb25ba7969	517ccc29-9db7-4d8d-8ea0-c8db6e9ac653	nextflow run main.nf --input data/greetings.csv --character stegosaurus -resume
    ```

Chaque ligne vous donne l'horodatage, la durée, le nom de l'exécution, le statut, l'ID de révision, l'ID de session et la ligne de commande complète d'une exécution lancée depuis ce répertoire.

Regardez les deux dernières lignes : ce sont deux invocations distinctes (l'une normale, l'autre avec `-resume`) de la même commande, et elles partagent le même ID de session.
L'ID de session ne change que lorsque vous lancez une nouvelle exécution véritablement nouvelle ; l'utilisation de `-resume` le conserve, ce qui permet à Nextflow de savoir quel cache réutiliser.

### 2.2. Utiliser `nextflow log` pour une vue plus lisible

Lire le fichier d'historique brut fonctionne, mais `nextflow log` formate les mêmes informations avec un en-tête :

```bash
nextflow log
```

??? success "Sortie de la commande"

    ```console linenums="1"
    TIMESTAMP          	DURATION	RUN NAME            	STATUS	REVISION ID	SESSION ID                          	COMMAND
    2026-09-13 01:47:58	2.8s    	determined_ramanujan	OK    	ce74f81996 	03078950-5011-414f-992b-8fe1bd219ac2	nextflow run main.nf --input data/greetings.csv --character turkey
    2026-09-13 01:48:03	3s      	prickly_cuvier      	OK    	ce74f81996 	5a9e8bf7-a7db-4c2c-b1d8-3bd4c7266528	nextflow run main.nf --input data/greetings.csv --character tux
    2026-09-13 01:48:09	2.3s    	elegant_panini      	OK    	ce74f81996 	517ccc29-9db7-4d8d-8ea0-c8db6e9ac653	nextflow run main.nf --input data/greetings.csv --character stegosaurus
    2026-09-13 01:48:14	1.5s    	deadly_lamport      	OK    	ce74f81996 	517ccc29-9db7-4d8d-8ea0-c8db6e9ac653	nextflow run main.nf --input data/greetings.csv --character stegosaurus -resume
    ```

Nextflow regroupe les informations de cache qu'il utilise pour `-resume` sous `.nextflow/cache`, indexées par ID de session.
C'est pourquoi rechercher le bon nom d'exécution ou l'ID de session ici est la première étape chaque fois que vous devez investiguer ou nettoyer une exécution passée.

### À retenir

Vous savez où Nextflow enregistre l'historique des exécutions passées, et comment l'inspecter avec `nextflow log`.

### Et ensuite ?

Apprenez à supprimer les anciens répertoires de travail dont vous n'avez plus besoin.

---

## 3. Supprimer les anciens répertoires de travail

Chaque exécution laisse ses répertoires de tâches dans `work/`, même après que vous avez copié les sorties qui vous intéressent dans `results/`.
Exécutez suffisamment de pipelines pendant le développement et ces sous-répertoires s'accumulent, c'est pourquoi Nextflow fournit `nextflow clean` pour supprimer ceux dont vous n'avez plus besoin.

### 3.1. Définir les critères de suppression

`nextflow clean` propose plusieurs façons de sélectionner ce qu'il faut supprimer ; consultez la [documentation de référence](https://www.nextflow.io/docs/latest/reference/cli.html#clean) pour la liste complète.
Ici, vous allez supprimer tout ce qui provient des exécutions antérieures à une exécution donnée, en utilisant son nom d'exécution.

Recherchez l'exécution la plus récente que vous souhaitez conserver avec `nextflow log` ; dans l'[exemple de la section 2.2](#22-use-nextflow-log-for-a-friendlier-view), il s'agit de `elegant_panini`, la dernière exécution normale avant celle avec `-resume`.
Le nom d'exécution est la chaîne en deux parties générée automatiquement, visible dans la ligne de console `Launching (...)` ou dans la colonne `RUN NAME` de `nextflow log`.

### 3.2. Effectuer une simulation

Ajoutez d'abord `-n` pour vérifier ce qu'une commande donnée supprimerait sans réellement rien supprimer :

```bash
nextflow clean -before elegant_panini -n
```

??? success "Sortie de la commande"

    ```console
    Would remove /workspaces/training/nextflow-run/work/e5/9ec3fe24deb5d1ffa478dff5e5e1d4
    Would remove /workspaces/training/nextflow-run/work/75/36b1ece86b213c331852d0a3dc659b
    Would remove /workspaces/training/nextflow-run/work/be/8157808f28699259b5ff37e4bf9f72
    Would remove /workspaces/training/nextflow-run/work/f3/d828b7ba6b315eb08e0f96168119f5
    Would remove /workspaces/training/nextflow-run/work/94/7a2d7ef67d81b6e66430f300816698
    Would remove /workspaces/training/nextflow-run/work/eb/e37c5c8638c45916de3d9794d6d22f
    Would remove /workspaces/training/nextflow-run/work/3f/f56a705a6370f6c615aff22d3f240a
    Would remove /workspaces/training/nextflow-run/work/25/4922fa35ba2c6087777ddce20e9f2f
    Would remove /workspaces/training/nextflow-run/work/a4/2de7bcc2f8fc974261f7280c1480f7
    Would remove /workspaces/training/nextflow-run/work/48/8c49bded446d4fcb78c8125318c5e3
    Would remove /workspaces/training/nextflow-run/work/8e/c9834a9bdd7c21de2991d959201a9d
    Would remove /workspaces/training/nextflow-run/work/d5/a9055dc8c817c6c1f436494685955c
    Would remove /workspaces/training/nextflow-run/work/27/9473df9e55e35c8869d9571b480926
    Would remove /workspaces/training/nextflow-run/work/d3/ebc066e1d844232b9f5ee0a7d87464
    Would remove /workspaces/training/nextflow-run/work/91/17ef1815ed06b8f177d0135d8ef0dc
    Would remove /workspaces/training/nextflow-run/work/12/0443f9fcff27a552d4bc73aaedf24a
    ```

Cela représente 16 répertoires de tâches : les 8 tâches de l'exécution `turkey` plus les 8 de l'exécution `tux`, exactement le nombre attendu pour deux exécutions complètes de ce pipeline à quatre processus.
L'exécution `elegant_panini` elle-même, ainsi que les tâches en cache réutilisées par l'exécution avec `-resume`, sont laissées intactes.

Votre sortie listera des noms de répertoires différents, et le nombre de lignes dépend du nombre d'exécutions que vous avez effectuées ; si vous ne voyez aucune ligne, soit le nom d'exécution ne correspond à aucun de ceux de votre journal, soit il n'y a rien à supprimer avant lui.

### 3.3. Procéder à la suppression

Une fois que la simulation vous convient, relancez la même commande avec `-f` à la place de `-n` :

```bash
nextflow clean -before elegant_panini -f
```

??? success "Sortie de la commande"

    ```console
    Removed /workspaces/training/nextflow-run/work/e5/9ec3fe24deb5d1ffa478dff5e5e1d4
    Removed /workspaces/training/nextflow-run/work/75/36b1ece86b213c331852d0a3dc659b
    Removed /workspaces/training/nextflow-run/work/be/8157808f28699259b5ff37e4bf9f72
    Removed /workspaces/training/nextflow-run/work/f3/d828b7ba6b315eb08e0f96168119f5
    Removed /workspaces/training/nextflow-run/work/94/7a2d7ef67d81b6e66430f300816698
    Removed /workspaces/training/nextflow-run/work/eb/e37c5c8638c45916de3d9794d6d22f
    Removed /workspaces/training/nextflow-run/work/3f/f56a705a6370f6c615aff22d3f240a
    Removed /workspaces/training/nextflow-run/work/25/4922fa35ba2c6087777ddce20e9f2f
    Removed /workspaces/training/nextflow-run/work/a4/2de7bcc2f8fc974261f7280c1480f7
    Removed /workspaces/training/nextflow-run/work/48/8c49bded446d4fcb78c8125318c5e3
    Removed /workspaces/training/nextflow-run/work/8e/c9834a9bdd7c21de2991d959201a9d
    Removed /workspaces/training/nextflow-run/work/d5/a9055dc8c817c6c1f436494685955c
    Removed /workspaces/training/nextflow-run/work/27/9473df9e55e35c8869d9571b480926
    Removed /workspaces/training/nextflow-run/work/d3/ebc066e1d844232b9f5ee0a7d87464
    Removed /workspaces/training/nextflow-run/work/91/17ef1815ed06b8f177d0135d8ef0dc
    Removed /workspaces/training/nextflow-run/work/12/0443f9fcff27a552d4bc73aaedf24a
    ```

`nextflow clean` vide les répertoires de tâches mais laisse en place les répertoires parents à deux caractères (comme `e5/`).

!!! warning "Avertissement"

    La suppression des répertoires de travail des exécutions passées les retire du cache de Nextflow et supprime toutes les sorties stockées uniquement à cet endroit.
    Cela compromet la capacité de Nextflow à reprendre l'exécution sans réexécuter les processus correspondants ; ne nettoyez donc que les exécutions dont vous êtes certain·e de ne plus avoir besoin de reprendre.
    C'est aussi pourquoi il est préférable de publier tout ce qui vous importe dans `results/` avec `mode 'copy'` plutôt que de vous fier au répertoire `work/` ou à un mode de publication `symlink`.

### À retenir

Vous savez comment supprimer les anciens répertoires de travail avec `nextflow clean`, et pourquoi cette opération vous prive de la possibilité de reprendre l'exécution à partir de ces exécutions.

### Et ensuite ?

Apprenez à exécuter des pipelines directement depuis des dépôts distants tels que GitHub dans la [Partie 4](./04_remote_repositories.md).

---

## Résumé

Dans cette partie, vous avez appris à :

- Générer un rapport d'exécution HTML avec `-with-report` et une chronologie d'exécution avec `-with-timeline`
- Inspecter l'historique des exécutions passées avec `nextflow log`
- Supprimer les anciens répertoires de travail avec `nextflow clean`, et comprendre le compromis sur la reprise d'exécution que cela implique
