# Partie 4 : Exécuter des pipelines distants

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Traduction assistée par IA - [en savoir plus et suggérer des améliorations](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

Jusqu'à présent, vous avez exécuté des scripts de workflow stockés localement.
En pratique, vous souhaiterez souvent exécuter des pipelines publiés dans des dépôts distants, comme GitHub, sans les télécharger vous-même.

Nextflow rend cela simple : vous pouvez exécuter n'importe quel pipeline directement depuis l'URL d'un dépôt Git.

---

## 1. Exécuter un pipeline depuis GitHub

La syntaxe de base pour exécuter un pipeline distant est `nextflow run <repository>`, où `<repository>` peut être un chemin de dépôt GitHub comme `nextflow-io/hello`, une URL complète, ou un chemin vers GitLab, Bitbucket, ou un autre service d'hébergement Git.

### 1.1. Lancer le pipeline

Exécutez le pipeline de démonstration officiel "hello" de Nextflow.
Il s'agit d'un pipeline différent, bien plus simple que celui que vous avez utilisé dans ce cours : il est antérieur au pipeline "Hello" utilisé tout au long de cette formation, et se contente d'afficher un message de bienvenue pour quelques langues codées en dur, donc ne vous attendez pas à l'entrée CSV ou à l'art ASCII auxquels vous êtes habitué·e.

```bash
nextflow run nextflow-io/hello
```

??? success "Sortie de la commande"

    ```console
    N E X T F L O W   ~  version 26.04.4

    Pulling nextflow-io/hello ...
     downloaded from https://github.com/nextflow-io/hello.git
    Launching `https://github.com/nextflow-io/hello` [sleepy_swanson] revision: 3c2cdc9823 [master]

    executor >  local (4)
    [ba/08236d] sayHello (4) | 4 of 4 ✔
    Ciao world!

    Hello world!

    Bonjour world!

    Hola world!
    ```

### 1.2. Trouver où le pipeline est mis en cache

La première fois que vous exécutez un pipeline distant, Nextflow le télécharge et le met en cache localement.
Les exécutions suivantes réutilisent la version mise en cache, sauf si vous demandez explicitement une mise à jour.

Par défaut, Nextflow enregistre les pipelines téléchargés dans `$NXF_HOME/assets`.
Pour savoir où un pipeline spécifique a été enregistré, et quelles révisions sont disponibles, interrogez Nextflow directement :

```bash
nextflow info nextflow-io/hello
```

??? success "Sortie de la commande"

    ```console
     project name: nextflow-io/hello
     repository  : https://github.com/nextflow-io/hello
     local path  : /workspaces/.nextflow/assets/.repos/nextflow-io/hello
     main script : main.nf
     revisions   :
     > master (default)
       mybranch
       testing
       v1.1 [t]
       v1.2 [t]
       v1.3 [t]
    ```

    Nextflow marque chaque révision que vous avez déjà extraite localement avec `>` ; les autres sont disponibles mais n'ont pas encore été récupérées dans une copie de travail.

Vous pouvez également lister tous les pipelines que vous avez téléchargés jusqu'à présent avec `nextflow list` :

```bash
nextflow list
```

??? success "Sortie de la commande"

    ```console
    nextflow-io/hello
    ```

Le cours [Use nf-core](../nfcore_use/01_run_demo.md#12-retrieve-the-pipeline-code) couvre ce mécanisme de mise en cache plus en détail, notamment comment parcourir le code source d'un pipeline téléchargé.

### À retenir

Vous savez comment exécuter un pipeline directement depuis un dépôt GitHub sans le télécharger vous-même, et où le trouver localement par la suite.

### Et ensuite ?

Apprenez à épingler une version spécifique d'un pipeline distant pour garantir la reproductibilité.

---

## 2. Spécifier une version pour la reproductibilité

Par défaut, Nextflow exécute la dernière révision de la branche par défaut.
Vous pouvez épingler une version particulière (tag), une branche, ou un commit en utilisant le flag `-r`.

### 2.1. Épingler une révision spécifique

```bash
nextflow run nextflow-io/hello -r v1.3
```

??? success "Sortie de la commande"

    ```console
    N E X T F L O W   ~  version 26.04.4

    Pulling nextflow-io/hello:v1.3 ...
     downloaded from https://github.com/nextflow-io/hello.git
    Launching `https://github.com/nextflow-io/hello` [sick_carson] revision: 2ce0b0e294 [v1.3]

    executor >  local (4)
    [61/e11f77] sayHello (4) | 4 of 4 ✔
    Ciao world!

    Bonjour world!

    Hello world!

    Hola world!
    ```

Nextflow récupère cette révision la première fois que vous la demandez, d'où les lignes `Pulling` et `downloaded from` ; demander la même révision ultérieurement passe directement à `Launching`.
Épingler une révision exacte est essentiel pour la reproductibilité.
Cela garantit que vous et vos collaborateur·trices exécutez exactement le même code de pipeline, indépendamment de ce qui a changé dans le dépôt depuis lors.

### 2.2. Les révisions s'appliquent uniquement à l'invocation en cours

Épingler une révision avec `-r` n'affecte que l'exécution pour laquelle vous la spécifiez : cela ne modifie pas ce qu'une exécution ultérieure avec un simple `nextflow run` utilisera.
Essayez d'exécuter à nouveau le pipeline sans `-r` :

```bash
nextflow run nextflow-io/hello
```

??? success "Sortie de la commande"

    ```console
    N E X T F L O W   ~  version 26.04.4

    Launching `https://github.com/nextflow-io/hello` [hungry_maxwell] revision: 3c2cdc9823 [master]

    executor >  local (4)
    [ba/08236d] sayHello (4) | 4 of 4 ✔
    Ciao world!

    Hello world!

    Bonjour world!

    Hola world!
    ```

Même si l'exécution précédente avait explicitement épinglé `v1.3`, cette exécution revient directement à la branche par défaut (`master`).
Nextflow conserve une copie de travail locale distincte pour chaque révision que vous avez utilisée, ce que montrent les marqueurs `>` dans `nextflow info`, mais il ne mémorise jamais laquelle vous avez exécutée en dernier.
Vous pouvez trouver le nom de la branche par défaut d'un pipeline en exécutant `nextflow info <pipeline>` ; c'est celle marquée `(default)`.
La reproductibilité vous incombe entièrement : passez toujours `-r` explicitement lorsque cela est important, plutôt que de supposer qu'une révision épinglée lors d'une exécution précédente s'applique encore.

### À retenir

Vous savez comment épingler un pipeline distant à une version, une branche ou un commit spécifique pour une exécution reproductible, et que l'épinglage ne s'applique qu'à cette seule invocation, pas aux exécutions ultérieures.

### Et ensuite ?

Vous avez couvert les fondamentaux de l'exécution et de la gestion des pipelines Nextflow.
Consultez le [Résumé du cours](next_steps.md) pour savoir comment continuer.

---

## Résumé

Dans cette partie, vous avez appris à :

- Exécuter un pipeline directement depuis un dépôt GitHub sans le télécharger
- Épingler un pipeline distant à une révision spécifique pour garantir la reproductibilité, et comprendre que l'épinglage ne s'applique qu'à cette seule invocation
