# Partie 3 : Utiliser les profils pour changer de configuration

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Traduction assistée par IA - [en savoir plus et suggérer des améliorations](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

Tout au long de la [Partie 1](./01_packaging_and_execution.md) et de la [Partie 2](./02_resources_and_retries.md), vous avez accumulé plusieurs options de configuration : l'empaquetage logiciel, la plateforme d'exécution et les allocations de ressources.
En pratique, vous souhaiterez souvent basculer entre des ensembles complets de ces options selon l'endroit où vous exécutez votre pipeline, par exemple un ordinateur portable pour le développement et un cluster HPC pour la production.

Nextflow vous permet de définir autant de [profils](https://nextflow.io/docs/latest/config.html#profiles) que vous le souhaitez, décrivant différentes configurations, et d'en sélectionner un (ou plusieurs) au moment de l'exécution avec un seul indicateur.

Vous en avez déjà utilisé un : le profil `test` de [Nextflow Run](../nextflow_run/index.md) remplace les paramètres d'entrée par un ensemble petit et bien défini.
Vous allez maintenant créer vos propres profils d'infrastructure et les combiner avec celui-ci.

---

## 1. Créer des profils pour différents environnements

### 1.1. Configurer les profils

Ajoutez deux profils à `nextflow.config` : un pour exécuter le pipeline sur un ordinateur portable classique avec Docker, et un pour un cluster HPC universitaire avec un ordonnanceur Slurm et Conda.

=== "Après"

    ```groovy title="nextflow.config" linenums="35" hl_lines="10-19"
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
        my_laptop {
            process.executor = 'local'
            docker.enabled = true
        }
        univ_hpc {
            process.executor = 'slurm'
            conda.enabled = true
            process.resourceLimits = [
                memory: 750.GB,
                cpus: 200,
                time: 30.d
            ]
        }
    }
    ```

=== "Avant"

    ```groovy title="nextflow.config" linenums="35"
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

Le profil `univ_hpc` définit également des limites de ressources, car cela est généralement requis sur une infrastructure HPC partagée.

### 1.2. Exécuter le workflow avec un profil

Sélectionnez un profil au moment de l'exécution avec `-profile`.

```bash
nextflow run main.nf -profile my_laptop
```

??? success "Sortie de la commande"

    ```console
    N E X T F L O W   ~  version 26.04.4

    Launching `main.nf` [nice_kimura] revision: c3c85dec78

    executor >  local (8)
    [df/b86b2f] sayHello (2)       | 3 of 3 ✔
    [b2/e2cb9c] convertToUpper (3) | 3 of 3 ✔
    [af/2ecff7] collectGreetings   | 1 of 1 ✔
    [c6/713c96] cowpy              | 1 of 1 ✔

    Outputs:

      /workspaces/training/execution-config/results

      first_output:
        - full_pipeline/intermediates/Bonjour-output.txt
        - full_pipeline/intermediates/Hello-output.txt
        - full_pipeline/intermediates/Hola-output.txt

      uppercased:
        - full_pipeline/intermediates/UPPER-Bonjour-output.txt
        - full_pipeline/intermediates/UPPER-Hello-output.txt
        - full_pipeline/intermediates/UPPER-Hola-output.txt

      collected: full_pipeline/intermediates/COLLECTED-batch-output.txt

      batch_report: full_pipeline/batch-report.txt

      cowpy_art: full_pipeline/cowpy-COLLECTED-batch-output.txt
    ```

!!! warning "Avertissement"

    Le profil `univ_hpc` ne fonctionnera pas dans l'environnement de formation, car aucun ordonnanceur Slurm n'est disponible.

Si vous trouvez d'autres paramètres qui vont toujours ensemble, ajoutez-les au profil correspondant.
Vous pouvez également créer des profils supplémentaires pour regrouper toute autre combinaison dont vous avez besoin.

### 1.3. Exécuter avec plusieurs profils

Les profils ne sont pas mutuellement exclusifs.
Vous pouvez en activer plusieurs à la fois avec `-profile <profil1>,<profil2>`.
Combinez `my_laptop` avec le profil `test` que vous connaissez déjà depuis Nextflow Run.

```bash
nextflow run main.nf -profile my_laptop,test
```

??? success "Sortie de la commande"

    ```console
    N E X T F L O W   ~  version 26.04.4

    Launching `main.nf` [romantic_volhard] revision: c3c85dec78

    executor >  local (8)
    [22/81de4f] sayHello (3)       | 3 of 3 ✔
    [93/e16b04] convertToUpper (1) | 3 of 3 ✔
    [ca/3e6a91] collectGreetings   | 1 of 1 ✔
    [16/0cc7e3] cowpy              | 1 of 1 ✔

    Outputs:

      /workspaces/training/execution-config/results

      first_output:
        - full_pipeline/intermediates/Hola-output.txt
        - full_pipeline/intermediates/Hello-output.txt
        - full_pipeline/intermediates/Bonjour-output.txt

      uppercased:
        - full_pipeline/intermediates/UPPER-Hello-output.txt
        - full_pipeline/intermediates/UPPER-Bonjour-output.txt
        - full_pipeline/intermediates/UPPER-Hola-output.txt

      collected: full_pipeline/intermediates/COLLECTED-test-output.txt

      batch_report: full_pipeline/test-report.txt

      cowpy_art: full_pipeline/cowpy-COLLECTED-test-output.txt
    ```

Les noms de fichiers individuels reprennent correctement `batch = 'test'` du profil `test` (`COLLECTED-test-output.txt`, etc.).

Si vous combinez des profils qui définissent la même option, Nextflow résout le conflit en utilisant la dernière valeur lue, c'est-à-dire celle qui apparaît le plus tard dans le fichier.
Si les paramètres en conflit proviennent de sources de configuration entièrement différentes, l'[ordre de priorité](https://www.nextflow.io/docs/latest/config.html) standard s'applique.

### À retenir

Vous savez comment définir des profils regroupant une configuration spécifique à l'infrastructure, en sélectionner un au moment de l'exécution avec `-profile`, combiner plusieurs profils en une seule exécution, et comment Nextflow résout les conflits lorsque plusieurs profils définissent la même option.

### Et ensuite ?

Apprenez à inspecter la configuration entièrement résolue avant d'exécuter quoi que ce soit.

---

## 2. Inspecter la configuration résolue

Vous avez déjà utilisé `nextflow config -profile test` dans [Nextflow Run](../nextflow_run/02_configure_pipeline.md) pour vérifier ce à quoi un seul profil se résout.
Cette commande devient particulièrement utile lorsque vous combinez plusieurs profils : comme vous venez de le voir, lorsque deux profils définissent la même option, il peut être difficile de déterminer manuellement quelle valeur l'emporte réellement.
La commande `nextflow config` résout tout cela pour vous, sans exécuter le pipeline.

### 2.1. Résoudre la configuration par défaut

```bash
nextflow config
```

??? success "Sortie de la commande"

    ```groovy
    params {
       input = 'data/greetings.csv'
       batch = 'batch'
       character = 'turkey'
    }

    docker {
       enabled = true
    }

    process {
       memory = '1 GB'
       withName:cowpy {
          conda = 'conda-forge::cowpy==1.1.5'
          memory = '2 GB'
          cpus = 2
       }
    }
    ```

C'est exactement ce qui s'appliquerait si vous exécutiez le pipeline sans indicateurs supplémentaires.

### 2.2. Résoudre la configuration avec des profils activés

Ajoutez les mêmes profils que vous utiliseriez pour une exécution réelle.

```bash
nextflow config -profile my_laptop,test
```

??? success "Sortie de la commande"

    ```groovy
    params {
       input = 'data/greetings.csv'
       batch = 'test'
       character = 'tux'
    }

    docker {
       enabled = true
    }

    process {
       memory = '1 GB'
       withName:cowpy {
          conda = 'conda-forge::cowpy==1.1.5'
          memory = '2 GB'
          cpus = 2
       }
       executor = 'local'
    }
    ```

La comparaison des deux résultats confirme ce qui a changé : `params.batch`, `params.character` et `process.executor` reflètent tous les profils `my_laptop,test`.
Cela devient particulièrement précieux pour les pipelines comportant de nombreuses couches de configuration, où déterminer manuellement les paramètres résolus serait fastidieux et source d'erreurs.

### À retenir

Vous savez comment utiliser `nextflow config` pour inspecter la configuration entièrement résolue pour n'importe quelle combinaison de profils, avant d'exécuter quoi que ce soit.

### Et ensuite ?

Vous avez couvert l'essentiel de la configuration des pipelines Nextflow.
Consultez le [Résumé du cours](next_steps.md) pour savoir comment continuer.

---

## Résumé

Dans cette partie, vous avez appris à :

- Définir des profils regroupant une configuration spécifique à l'infrastructure
- Combiner plusieurs profils en une seule exécution, et comprendre comment les conflits entre eux sont résolus
- Utiliser `nextflow config` pour inspecter la configuration entièrement résolue
