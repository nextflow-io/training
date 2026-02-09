# Débogage des Workflows

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Traduction assistée par IA - [en savoir plus et suggérer des améliorations](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

Le débogage est une compétence essentielle qui peut vous faire gagner des heures de frustration et vous aider à devenir un·e développeur·se Nextflow plus efficace. Tout au long de votre carrière, en particulier lorsque vous débutez, vous rencontrerez des bugs lors de la construction et de la maintenance de vos workflows. L'apprentissage d'approches de débogage systématiques vous aidera à identifier et à résoudre les problèmes rapidement.

### Objectifs d'apprentissage

Dans cette quête annexe, nous explorerons les **techniques de débogage systématiques** pour les workflows Nextflow :

- **Débogage des erreurs de syntaxe** : Utilisation efficace des fonctionnalités de l'IDE et des messages d'erreur de Nextflow
- **Débogage des canaux** : Diagnostic des problèmes de flux de données et de structure des canaux
- **Débogage des processus** : Investigation des échecs d'exécution et des problèmes de ressources
- **Outils de débogage intégrés** : Exploitation du mode preview, de l'exécution stub et des répertoires de travail de Nextflow
- **Approches systématiques** : Une méthodologie en quatre phases pour un débogage efficace

À la fin, vous disposerez d'une méthodologie de débogage robuste qui transforme les messages d'erreur frustrants en feuilles de route claires vers les solutions.

### Prérequis

Avant d'entreprendre cette quête annexe, vous devez :

- Avoir terminé le tutoriel [Hello Nextflow](../hello_nextflow/README.md) ou un cours équivalent pour débutants.
- Être à l'aise avec les concepts et mécanismes de base de Nextflow (processus, canaux, opérateurs)

**Optionnel :** Nous recommandons de terminer d'abord la quête annexe [Fonctionnalités IDE pour le Développement Nextflow](./ide_features.md).
Celle-ci couvre de manière exhaustive les fonctionnalités de l'IDE qui supportent le débogage (coloration syntaxique, détection d'erreurs, etc.), que nous utiliserons intensivement ici.

---

## 0. Commencer

#### Ouvrir l'espace de code de formation

Si vous ne l'avez pas encore fait, assurez-vous d'ouvrir l'environnement de formation comme décrit dans la [Configuration de l'Environnement](../envsetup/index.md).

[![Open in GitHub Codespaces](https://github.com/codespaces/badge.svg)](https://codespaces.new/nextflow-io/training?quickstart=1&ref=master)

#### Se déplacer dans le répertoire du projet

Déplaçons-nous dans le répertoire où se trouvent les fichiers pour ce tutoriel.

```bash
cd side-quests/debugging
```

Vous pouvez configurer VSCode pour se concentrer sur ce répertoire :

```bash
code .
```

#### Examiner le matériel

Vous trouverez un ensemble d'exemples de workflows avec divers types de bugs que nous utiliserons pour la pratique :

??? abstract "Contenu du répertoire"

    ```console
    .
    ├── bad_bash_var.nf
    ├── bad_channel_shape.nf
    ├── bad_channel_shape_viewed_debug.nf
    ├── bad_channel_shape_viewed.nf
    ├── bad_number_inputs.nf
    ├── badpractice_syntax.nf
    ├── bad_resources.nf
    ├── bad_syntax.nf
    ├── buggy_workflow.nf
    ├── data
    │   ├── sample_001.fastq.gz
    │   ├── sample_002.fastq.gz
    │   ├── sample_003.fastq.gz
    │   ├── sample_004.fastq.gz
    │   ├── sample_005.fastq.gz
    │   └── sample_data.csv
    ├── exhausted.nf
    ├── invalid_process.nf
    ├── missing_output.nf
    ├── missing_software.nf
    ├── missing_software_with_stub.nf
    ├── nextflow.config
    └── no_such_var.nf
    ```

Ces fichiers représentent des scénarios de débogage courants que vous rencontrerez dans le développement en conditions réelles.

#### Examiner l'exercice

Votre défi est d'exécuter chaque workflow, d'identifier la ou les erreurs, et de les corriger.

Pour chaque workflow bugué :

1. **Exécuter le workflow** et observer l'erreur
2. **Analyser le message d'erreur** : que vous dit Nextflow ?
3. **Localiser le problème** dans le code en utilisant les indices fournis
4. **Corriger le bug** et vérifier que votre solution fonctionne
5. **Réinitialiser le fichier** avant de passer à la section suivante (utilisez `git checkout <filename>`)

Les exercices progressent des erreurs de syntaxe simples vers des problèmes d'exécution plus subtils.
Les solutions sont discutées en ligne, mais essayez de résoudre chaque problème vous-même avant de lire la suite.

#### Liste de vérification de préparation

Pensez-vous être prêt·e à vous lancer ?

- [ ] Je comprends l'objectif de ce cours et ses prérequis
- [ ] Mon espace de code est opérationnel
- [ ] J'ai défini mon répertoire de travail de manière appropriée
- [ ] Je comprends l'exercice

Si vous pouvez cocher toutes les cases, vous êtes prêt·e à commencer.

---

## 1. Erreurs de Syntaxe

Les erreurs de syntaxe sont le type d'erreur le plus courant que vous rencontrerez lors de l'écriture de code Nextflow. Elles se produisent lorsque le code ne se conforme pas aux règles de syntaxe attendues du DSL Nextflow. Ces erreurs empêchent votre workflow de s'exécuter du tout, il est donc important d'apprendre à les identifier et à les corriger rapidement.

### 1.1. Accolades manquantes

L'une des erreurs de syntaxe les plus courantes, et parfois l'une des plus complexes à déboguer, est **les crochets manquants ou mal appariés**.

Commençons par un exemple pratique.

#### Exécuter le pipeline

```bash
nextflow run bad_syntax.nf
```

??? failure "Sortie de la commande"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `bad_syntax.nf` [stupefied_bhabha] DSL2 - revision: ca6327fad2

    Error bad_syntax.nf:24:1: Unexpected input: '<EOF>'

    ERROR ~ Script compilation failed

     -- Check '.nextflow.log' file for details
    ```

**Éléments clés des messages d'erreur de syntaxe :**

- **Fichier et emplacement** : Indique quel fichier et quelle ligne/colonne contiennent l'erreur (`bad_syntax.nf:24:1`)
- **Description de l'erreur** : Explique ce que l'analyseur a trouvé qu'il n'attendait pas (`Unexpected input: '<EOF>'`)
- **Indicateur EOF** : Le message `<EOF>` (End Of File) indique que l'analyseur a atteint la fin du fichier en attendant toujours plus de contenu - un signe classique d'accolades non fermées

#### Vérifier le code

Maintenant, examinons `bad_syntax.nf` pour comprendre ce qui cause l'erreur :

```groovy title="bad_syntax.nf" hl_lines="14" linenums="1"
#!/usr/bin/env nextflow

process PROCESS_FILES {
    input:
    val sample_name

    output:
    path "${sample_name}_output.txt"

    script:
    """
    echo "Processing ${sample_name}" > ${sample_name}_output.txt
    """
// Accolade fermante manquante pour le processus

workflow {

    // Crée un canal d'entrée
    input_ch = channel.of('sample1', 'sample2', 'sample3')

    // Appelle le processus avec le canal d'entrée
    PROCESS_FILES(input_ch)
}
```

Pour les besoins de cet exemple, nous avons laissé un commentaire pour vous montrer où se trouve l'erreur. L'extension VSCode Nextflow devrait également vous donner quelques indices sur ce qui pourrait ne pas aller, en mettant l'accolade mal appariée en rouge et en soulignant la fin prématurée du fichier :

![Bad syntax](img/bad_syntax.png)

**Stratégie de débogage pour les erreurs d'accolades :**

1. Utilisez la correspondance d'accolades de VS Code (placez le curseur à côté d'une accolade)
2. Vérifiez le panneau Problèmes pour les messages liés aux accolades
3. Assurez-vous que chaque `{` d'ouverture a un `}` de fermeture correspondant

#### Corriger le code

Remplacez le commentaire par l'accolade fermante manquante :

=== "Après"

    ```groovy title="bad_syntax.nf" hl_lines="14" linenums="1"
    #!/usr/bin/env nextflow

    process PROCESS_FILES {
        input:
        val sample_name

        output:
        path "${sample_name}_output.txt"

        script:
        """
        echo "Processing ${sample_name}" > ${sample_name}_output.txt
        """
    }  // Ajoute l'accolade fermante manquante

    workflow {

        // Crée un canal d'entrée
        input_ch = channel.of('sample1', 'sample2', 'sample3')

        // Appelle le processus avec le canal d'entrée
        PROCESS_FILES(input_ch)
    }
    ```

=== "Avant"

    ```groovy title="bad_syntax.nf" hl_lines="14" linenums="1"
    #!/usr/bin/env nextflow

    process PROCESS_FILES {
        input:
        val sample_name

        output:
        path "${sample_name}_output.txt"

        script:
        """
        echo "Processing ${sample_name}" > ${sample_name}_output.txt
        """
    // Accolade fermante manquante pour le processus

    workflow {

        // Crée un canal d'entrée
        input_ch = channel.of('sample1', 'sample2', 'sample3')

        // Appelle le processus avec le canal d'entrée
        PROCESS_FILES(input_ch)
    }
    ```

#### Exécuter le pipeline

Maintenant, exécutez à nouveau le workflow pour confirmer qu'il fonctionne :

```bash
nextflow run bad_syntax.nf
```

??? success "Sortie de la commande"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `bad_syntax.nf` [insane_faggin] DSL2 - revision: 961938ee2b

    executor >  local (3)
    [48/cd7f54] PROCESS_FILES (1) | 3 of 3 ✔
    ```

### 1.2. Utilisation de mots-clés ou directives de processus incorrects

Une autre erreur de syntaxe courante est une **définition de processus invalide**. Cela peut se produire si vous oubliez de définir des blocs requis ou si vous utilisez des directives incorrectes dans la définition du processus.

#### Exécuter le pipeline

```bash
nextflow run invalid_process.nf
```

??? failure "Sortie de la commande"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `invalid_process.nf` [nasty_jepsen] DSL2 - revision: da9758d614

    Error invalid_process.nf:3:1: Invalid process definition -- check for missing or out-of-order section labels
    │   3 | process PROCESS_FILES {
    │     | ^^^^^^^^^^^^^^^^^^^^^^^
    │   4 |     inputs:
    │   5 |     val sample_name
    │   6 |
    ╰   7 |     output:

    ERROR ~ Script compilation failed

     -- Check '.nextflow.log' file for details
    ```

#### Vérifier le code

L'erreur indique une "Définition de processus invalide" et montre le contexte autour du problème. En regardant les lignes 3-7, nous pouvons voir `inputs:` à la ligne 4, qui est le problème. Examinons `invalid_process.nf` :

```groovy title="invalid_process.nf" hl_lines="4" linenums="1"
#!/usr/bin/env nextflow

process PROCESS_FILES {
    inputs:  // ERREUR : Devrait être 'input' et non 'inputs'
    val sample_name

    output:
    path "${sample_name}_output.txt"

    script:
    """
    echo "Processing ${sample_name}" > ${sample_name}_output.txt
    """
}

workflow {

    // Crée un canal d'entrée
    input_ch = channel.of('sample1', 'sample2', 'sample3')

    // Appelle le processus avec le canal d'entrée
    PROCESS_FILES(input_ch)
}
```

En regardant la ligne 4 dans le contexte de l'erreur, nous pouvons repérer le problème : nous utilisons `inputs` au lieu de la directive correcte `input`. L'extension VSCode Nextflow signalera également ceci :

![Invalid process message](img/invalid_process_message.png)

#### Corriger le code

Remplacez le mot-clé incorrect par le bon en consultant [la documentation](https://www.nextflow.io/docs/latest/process.html#) :

=== "Après"

    ```groovy title="invalid_process.nf" hl_lines="4" linenums="1"
    #!/usr/bin/env nextflow

    process PROCESS_FILES {
        input:  // Corrigé : Changé 'inputs' en 'input'
        val sample_name

        output:
        path "${sample_name}_output.txt"

        script:
        """
        echo "Processing ${sample_name}" > ${sample_name}_output.txt
        """
    }

    workflow {

        // Crée un canal d'entrée
        input_ch = channel.of('sample1', 'sample2', 'sample3')

        // Appelle le processus avec le canal d'entrée
        PROCESS_FILES(input_ch)
    }
    ```

=== "Avant"

    ```groovy title="invalid_process.nf" hl_lines="4" linenums="1"
    #!/usr/bin/env nextflow

    process PROCESS_FILES {
        inputs:  // ERREUR : Devrait être 'input' et non 'inputs'
        val sample_name

        output:
        path "${sample_name}_output.txt"

        script:
        """
        echo "Processing ${sample_name}" > ${sample_name}_output.txt
        """
    }

    workflow {

        // Crée un canal d'entrée
        input_ch = channel.of('sample1', 'sample2', 'sample3')

        // Appelle le processus avec le canal d'entrée
        PROCESS_FILES(input_ch)
    }
    ```

#### Exécuter le pipeline

Maintenant, exécutez à nouveau le workflow pour confirmer qu'il fonctionne :

```bash
nextflow run invalid_process.nf
```

??? success "Sortie de la commande"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `invalid_process.nf` [silly_fermi] DSL2 - revision: 961938ee2b

    executor >  local (3)
    [b7/76cd9d] PROCESS_FILES (2) | 3 of 3 ✔
    ```

### 1.3. Utilisation de mauvais noms de variables

Les noms de variables que vous utilisez dans vos blocs de script doivent être valides, dérivés soit des entrées, soit du code groovy inséré avant le script. Mais lorsque vous gérez la complexité au début du développement de pipeline, il est facile de faire des erreurs dans la dénomination des variables, et Nextflow vous le fera savoir rapidement.

#### Exécuter le pipeline

```bash
nextflow run no_such_var.nf
```

??? failure "Sortie de la commande"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `no_such_var.nf` [gloomy_meninsky] DSL2 - revision: 0c4d3bc28c

    Error no_such_var.nf:17:39: `undefined_var` is not defined
    │  17 |     echo "Using undefined variable: ${undefined_var}" >> ${output_pref
    ╰     |                                       ^^^^^^^^^^^^^

    ERROR ~ Script compilation failed

     -- Check '.nextflow.log' file for details
    ```

L'erreur est détectée au moment de la compilation et pointe directement vers la variable non définie à la ligne 17, avec un caret indiquant exactement où se trouve le problème.

#### Vérifier le code

Examinons `no_such_var.nf` :

```groovy title="no_such_var.nf" hl_lines="17" linenums="1"
#!/usr/bin/env nextflow

process PROCESS_FILES {
    input:
    val sample_name

    output:
    path "${sample_name}_processed.txt"

    script:
    // Définit les variables dans le code Groovy avant le script
    def output_prefix = "${sample_name}_processed"
    def timestamp = new Date().format("yyyy-MM-dd")

    """
    echo "Processing ${sample_name} on ${timestamp}" > ${output_prefix}.txt
    echo "Using undefined variable: ${undefined_var}" >> ${output_prefix}.txt  // ERREUR : undefined_var non définie
    """
}

workflow {
    input_ch = channel.of('sample1', 'sample2', 'sample3')
    PROCESS_FILES(input_ch)
}
```

Le message d'erreur indique que la variable n'est pas reconnue dans le template de script, et vous devriez pouvoir voir `${undefined_var}` utilisé dans le bloc script, mais non défini ailleurs.

#### Corriger le code

Si vous obtenez une erreur 'No such variable', vous pouvez la corriger soit en définissant la variable (en corrigeant les noms de variables d'entrée ou en modifiant le code groovy avant le script), soit en la supprimant du bloc script si elle n'est pas nécessaire :

=== "Après"

    ```groovy title="no_such_var.nf" hl_lines="15-17" linenums="1"
    #!/usr/bin/env nextflow

    process PROCESS_FILES {
        input:
        val sample_name

        output:
        path "${sample_name}_output.txt"

        script:
        // Définit les variables dans le code Groovy avant le script
        def output_prefix = "${sample_name}_processed"
        def timestamp = new Date().format("yyyy-MM-dd")

        """
        echo "Processing ${sample_name} on ${timestamp}" > ${output_prefix}.txt
        """  // Supprime la ligne avec undefined_var
    }

    workflow {
        input_ch = channel.of('sample1', 'sample2', 'sample3')
        PROCESS_FILES(input_ch)
    }
    ```

=== "Avant"

    ```groovy title="no_such_var.nf" hl_lines="17" linenums="1"
    #!/usr/bin/env nextflow

    process PROCESS_FILES {
        input:
        val sample_name

        output:
        path "${sample_name}_output.txt"

        script:
        // Définit les variables dans le code Groovy avant le script
        def output_prefix = "${sample_name}_processed"
        def timestamp = new Date().format("yyyy-MM-dd")

        """
        echo "Processing ${sample_name} on ${timestamp}" > ${output_prefix}.txt
        echo "Using undefined variable: ${undefined_var}" >> ${output_prefix}.txt  // ERREUR : undefined_var non définie
        """
    }

    workflow {
        input_ch = channel.of('sample1', 'sample2', 'sample3')
        PROCESS_FILES(input_ch)
    }
    ```

#### Exécuter le pipeline

Maintenant, exécutez à nouveau le workflow pour confirmer qu'il fonctionne :

```bash
nextflow run no_such_var.nf
```

??? success "Sortie de la commande"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `no_such_var.nf` [suspicious_venter] DSL2 - revision: 6ba490f7c5

    executor >  local (3)
    [21/237300] PROCESS_FILES (2) | 3 of 3 ✔
    ```

### 1.4. Mauvaise utilisation des variables Bash

En débutant avec Nextflow, il peut être difficile de comprendre la différence entre les variables Nextflow (Groovy) et Bash. Cela peut générer une autre forme d'erreur de variable incorrecte qui apparaît lors de la tentative d'utilisation de variables dans le contenu Bash du bloc script.

#### Exécuter le pipeline

```bash
nextflow run bad_bash_var.nf
```

??? failure "Sortie de la commande"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `bad_bash_var.nf` [infallible_mandelbrot] DSL2 - revision: 0853c11080

    Error bad_bash_var.nf:13:42: `prefix` is not defined
    │  13 |     echo "Processing ${sample_name}" > ${prefix}.txt
    ╰     |                                          ^^^^^^

    ERROR ~ Script compilation failed

     -- Check '.nextflow.log' file for details
    ```

#### Vérifier le code

L'erreur pointe vers la ligne 13 où `${prefix}` est utilisé. Examinons `bad_bash_var.nf` pour voir ce qui cause le problème :

```groovy title="bad_bash_var.nf" hl_lines="13" linenums="1"
#!/usr/bin/env nextflow

process PROCESS_FILES {
    input:
    val sample_name

    output:
    path "${sample_name}_output.txt"

    script:
    """
    prefix="${sample_name}_output"
    echo "Processing ${sample_name}" > ${prefix}.txt  # ERROR: ${prefix} is Groovy syntax, not Bash
    """
}
```

Dans cet exemple, nous définissons la variable `prefix` en Bash, mais dans un processus Nextflow, la syntaxe `$` que nous avons utilisée pour y faire référence (`${prefix}`) est interprétée comme une variable Groovy, pas Bash. La variable n'existe pas dans le contexte Groovy, donc nous obtenons une erreur 'no such variable'.

#### Corriger le code

Si vous voulez utiliser une variable Bash, vous devez échapper le signe dollar comme ceci :

=== "Après"

    ```groovy title="bad_bash_var.nf" hl_lines="13" linenums="1"
    #!/usr/bin/env nextflow

    process PROCESS_FILES {
        input:
        val sample_name

        output:
        path "${sample_name}_output.txt"

        script:
        """
        prefix="${sample_name}_output"
        echo "Processing ${sample_name}" > \${prefix}.txt  # Fixed: Escaped the dollar sign
        """
    }

    workflow {
        input_ch = channel.of('sample1', 'sample2', 'sample3')
        PROCESS_FILES(input_ch)
    }
    ```

=== "Avant"

    ```groovy title="bad_bash_var.nf" hl_lines="13" linenums="1"
    #!/usr/bin/env nextflow

    process PROCESS_FILES {
        input:
        val sample_name

        output:
        path "${sample_name}_output.txt"

        script:
        """
        prefix="${sample_name}_output"
        echo "Processing ${sample_name}" > ${prefix}.txt  # ERROR: ${prefix} is Groovy syntax, not Bash
        """
    }
    ```

Cela indique à Nextflow d'interpréter ceci comme une variable Bash.

#### Exécuter le pipeline

Maintenant, exécutez à nouveau le workflow pour confirmer qu'il fonctionne :

```bash
nextflow run bad_bash_var.nf
```

??? success "Sortie de la commande"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `bad_bash_var.nf` [naughty_franklin] DSL2 - revision: 58c1c83709

    executor >  local (3)
    [4e/560285] PROCESS_FILES (2) | 3 of 3 ✔
    ```

!!! tip "Variables Groovy vs Variables Bash"

    Pour les manipulations de variables simples comme la concaténation de chaînes ou les opérations de préfixe/suffixe, il est généralement plus lisible d'utiliser des variables Groovy dans la section script plutôt que des variables Bash dans le bloc script :

    ```groovy linenums="1"
    script:
    def output_prefix = "${sample_name}_processed"
    def output_file = "${output_prefix}.txt"
    """
    echo "Processing ${sample_name}" > ${output_file}
    """
    ```

    Cette approche évite le besoin d'échapper les signes dollar et rend le code plus facile à lire et à maintenir.

### 1.5. Instructions en Dehors du Bloc Workflow

L'extension VSCode Nextflow met en évidence les problèmes de structure de code qui causeront des erreurs. Un exemple courant est la définition de canaux en dehors du bloc `workflow {}` - ceci est maintenant imposé comme une erreur de syntaxe.

#### Exécuter le pipeline

```bash
nextflow run badpractice_syntax.nf
```

??? failure "Sortie de la commande"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `badpractice_syntax.nf` [intergalactic_colden] DSL2 - revision: 5e4b291bde

    Error badpractice_syntax.nf:3:1: Statements cannot be mixed with script declarations -- move statements into a process or workflow
    │   3 | input_ch = channel.of('sample1', 'sample2', 'sample3')
    ╰     | ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

    ERROR ~ Script compilation failed

     -- Check '.nextflow.log' file for details
    ```

Le message d'erreur indique clairement le problème : les instructions (comme les définitions de canaux) ne peuvent pas être mélangées avec les déclarations de script en dehors d'un bloc workflow ou process.

#### Vérifier le code

Examinons `badpractice_syntax.nf` pour voir ce qui cause l'erreur :

```groovy title="badpractice_syntax.nf" hl_lines="3" linenums="1"
#!/usr/bin/env nextflow

input_ch = channel.of('sample1', 'sample2', 'sample3')  // ERREUR : Canal défini en dehors du workflow

process PROCESS_FILES {
    input:
    val sample_name

    output:
    path "${sample_name}_processed.txt"

    script:
    // Définit les variables dans le code Groovy avant le script
    def output_prefix = "${sample_name}_processed"
    def timestamp = new Date().format("yyyy-MM-dd")

    """
    echo "Processing ${sample_name} on ${timestamp}" > ${output_prefix}.txt
    """
}

workflow {
    PROCESS_FILES(input_ch)
}
```

L'extension VSCode soulignera également la variable `input_ch` comme étant définie en dehors du bloc workflow :

![Non-lethal syntax error](img/nonlethal.png)

#### Corriger le code

Déplacez la définition du canal à l'intérieur du bloc workflow :

=== "Après"

    ```groovy title="badpractice_syntax.nf" hl_lines="21" linenums="1"
    #!/usr/bin/env nextflow

    process PROCESS_FILES {
        input:
        val sample_name

        output:
        path "${sample_name}_processed.txt"

        script:
        // Définit les variables dans le code Groovy avant le script
        def output_prefix = "${sample_name}_processed"
        def timestamp = new Date().format("yyyy-MM-dd")

        """
        echo "Processing ${sample_name} on ${timestamp}" > ${output_prefix}.txt
        """
    }

    workflow {
        input_ch = channel.of('sample1', 'sample2', 'sample3')  // Déplacé à l'intérieur du bloc workflow
        PROCESS_FILES(input_ch)
    }
    ```

=== "Avant"

    ```groovy title="badpractice_syntax.nf" hl_lines="3" linenums="1"
    #!/usr/bin/env nextflow

    input_ch = channel.of('sample1', 'sample2', 'sample3')  // ERREUR : Canal défini en dehors du workflow

    process PROCESS_FILES {
        input:
        val sample_name

        output:
        path "${sample_name}_processed.txt"

        script:
        // Définit les variables dans le code Groovy avant le script
        def output_prefix = "${sample_name}_processed"
        def timestamp = new Date().format("yyyy-MM-dd")

        """
        echo "Processing ${sample_name} on ${timestamp}" > ${output_prefix}.txt
        """
    }

    workflow {
        PROCESS_FILES(input_ch)
    }
    ```

#### Exécuter le pipeline

Exécutez à nouveau le workflow pour confirmer que la correction fonctionne :

```bash
nextflow run badpractice_syntax.nf
```

??? success "Sortie de la commande"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `badpractice_syntax.nf` [naughty_ochoa] DSL2 - revision: 5e4b291bde

    executor >  local (3)
    [6a/84a608] PROCESS_FILES (2) | 3 of 3 ✔
    ```

Gardez vos canaux d'entrée définis dans le bloc workflow, et en général suivez toutes les autres recommandations faites par l'extension.

### À retenir

Vous pouvez identifier et corriger systématiquement les erreurs de syntaxe en utilisant les messages d'erreur de Nextflow et les indicateurs visuels de l'IDE. Les erreurs de syntaxe courantes incluent les accolades manquantes, les mots-clés de processus incorrects, les variables non définies et l'utilisation inappropriée des variables Bash vs Nextflow. L'extension VSCode aide à détecter bon nombre de ces erreurs avant l'exécution. Avec ces compétences de débogage de syntaxe dans votre boîte à outils, vous serez capable de résoudre rapidement les erreurs de syntaxe Nextflow les plus courantes et de passer à la résolution de problèmes d'exécution plus complexes.

### Et ensuite ?

Apprenez à déboguer des erreurs de structure de canal plus complexes qui se produisent même lorsque la syntaxe est correcte.

---

## 2. Erreurs de Structure des Canaux

Les erreurs de structure des canaux sont plus subtiles que les erreurs de syntaxe car le code est syntaxiquement correct, mais les formes de données ne correspondent pas à ce que les processus attendent. Nextflow tentera d'exécuter le pipeline, mais pourrait constater que le nombre d'entrées ne correspond pas à ce qu'il attend et échouer. Ces erreurs n'apparaissent généralement qu'à l'exécution et nécessitent une compréhension des données circulant dans votre workflow.

!!! tip "Débogage des Canaux avec `.view()`"

    Tout au long de cette section, rappelez-vous que vous pouvez utiliser l'opérateur `.view()` pour inspecter le contenu des canaux à n'importe quel point de votre workflow. C'est l'un des outils de débogage les plus puissants pour comprendre les problèmes de structure des canaux. Nous explorerons cette technique en détail dans la section 2.4, mais n'hésitez pas à l'utiliser pendant que vous travaillez sur les exemples.

    ```groovy
    my_channel.view()  // Affiche ce qui circule dans le canal
    ```

### 2.1. Mauvais Nombre de Canaux d'Entrée

Cette erreur se produit lorsque vous passez un nombre différent de canaux que ce qu'un processus attend.

#### Exécuter le pipeline

```bash
nextflow run bad_number_inputs.nf
```

??? failure "Sortie de la commande"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `bad_number_inputs.nf` [happy_swartz] DSL2 - revision: d83e58dcd3

    Error bad_number_inputs.nf:23:5: Incorrect number of call arguments, expected 1 but received 2
    │  23 |     PROCESS_FILES(samples_ch, files_ch)
    ╰     |     ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

    ERROR ~ Script compilation failed

     -- Check '.nextflow.log' file for details
    ```

#### Vérifier le code

Le message d'erreur indique clairement que l'appel attendait 1 argument mais en a reçu 2, et pointe vers la ligne 23. Examinons `bad_number_inputs.nf` :

```groovy title="bad_number_inputs.nf" hl_lines="5 23" linenums="1"
#!/usr/bin/env nextflow

process PROCESS_FILES {
    input:
        val sample_name  // Le processus n'attend qu'une seule entrée

    output:
        path "${sample_name}_output.txt"

    script:
    """
    echo "Processing ${sample_name}" > ${sample_name}_output.txt
    """
}

workflow {

    // Crée deux canaux séparés
    samples_ch = channel.of('sample1', 'sample2', 'sample3')
    files_ch = channel.of('file1.txt', 'file2.txt', 'file3.txt')

    // ERREUR : Passe 2 canaux mais le processus n'en attend qu'un seul
    PROCESS_FILES(samples_ch, files_ch)
}
```

Vous devriez voir l'appel `PROCESS_FILES` mal apparié, fournissant plusieurs canaux d'entrée alors que le processus n'en définit qu'un. L'extension VSCode soulignera également l'appel du processus en rouge, et fournira un message de diagnostic au survol :

![Incorrect number of args message](img/incorrect_num_args.png)

#### Corriger le code

Pour cet exemple spécifique, le processus attend un seul canal et ne nécessite pas le second canal, nous pouvons donc le corriger en passant uniquement le canal `samples_ch` :

=== "Après"

    ```groovy title="bad_number_inputs.nf" hl_lines="23" linenums="1"
    #!/usr/bin/env nextflow

    process PROCESS_FILES {
        input:
            val sample_name  // Le processus n'attend qu'une seule entrée

        output:
            path "${sample_name}_output.txt"

        script:
        """
        echo "Processing ${sample_name}" > ${sample_name}_output.txt
        """
    }

    workflow {

        // Crée deux canaux séparés
        samples_ch = channel.of('sample1', 'sample2', 'sample3')
        files_ch = channel.of('file1.txt', 'file2.txt', 'file3.txt')

        // Corrigé : Passe uniquement le canal que le processus attend
        PROCESS_FILES(samples_ch)
    }
    ```

=== "Avant"

    ```groovy title="bad_number_inputs.nf" hl_lines="5 23" linenums="1"
    #!/usr/bin/env nextflow

    process PROCESS_FILES {
        input:
            val sample_name  // Le processus n'attend qu'une seule entrée

        output:
            path "${sample_name}_output.txt"

        script:
        """
        echo "Processing ${sample_name}" > ${sample_name}_output.txt
        """
    }

    workflow {

        // Crée deux canaux séparés
        samples_ch = channel.of('sample1', 'sample2', 'sample3')
        files_ch = channel.of('file1.txt', 'file2.txt', 'file3.txt')

        // ERREUR : Passe 2 canaux mais le processus n'en attend qu'un seul
        PROCESS_FILES(samples_ch, files_ch)
    }
    ```

#### Exécuter le pipeline

```bash
nextflow run bad_number_inputs.nf
```

??? success "Sortie de la commande"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `bad_number_inputs.nf` [big_euler] DSL2 - revision: e302bd87be

    executor >  local (3)
    [48/497f7b] PROCESS_FILES (3) | 3 of 3 ✔
    ```

Plus couramment que dans cet exemple, vous pourriez ajouter des entrées supplémentaires à un processus et oublier de mettre à jour l'appel du workflow en conséquence, ce qui peut conduire à ce type d'erreur. Heureusement, c'est l'une des erreurs les plus faciles à comprendre et à corriger, car le message d'erreur est assez clair sur le décalage.

### 2.2. Épuisement des Canaux (Le Processus S'Exécute Moins de Fois Que Prévu)

Certaines erreurs de structure de canal sont beaucoup plus subtiles et ne produisent aucune erreur du tout. Probablement la plus courante de celles-ci reflète un défi auquel les nouveaux·elles utilisateur·trices de Nextflow sont confronté·es en comprenant que les canaux de file peuvent être épuisés et manquer d'éléments, ce qui signifie que le workflow se termine prématurément.

#### Exécuter le pipeline

```bash
nextflow run exhausted.nf
```

??? success "Sortie de la commande"

```console title="Exhausted channel output"
 N E X T F L O W   ~  version 25.10.2

Launching `exhausted.nf` [extravagant_gauss] DSL2 - revision: 08cff7ba2a

executor >  local (1)
[bd/f61fff] PROCESS_FILES (1) [100%] 1 of 1 ✔
```

Ce workflow se termine sans erreur, mais il ne traite qu'un seul échantillon !

#### Vérifier le code

Examinons `exhausted.nf` pour voir si c'est correct :

```groovy title="exhausted.nf" hl_lines="23 24" linenums="1"
#!/usr/bin/env nextflow

process PROCESS_FILES {
    input:
    val reference
    val sample_name

    output:
    path "${output_prefix}.txt"

    script:
    // Définit les variables dans le code Groovy avant le script
    output_prefix = "${reference}_${sample_name}"
    def timestamp = new Date().format("yyyy-MM-dd")

    """
    echo "Processing ${sample_name} on ${timestamp}" > ${output_prefix}.txt
    """
}

workflow {

    reference_ch = channel.of('baseline_reference')
    input_ch = channel.of('sample1', 'sample2', 'sample3')

    PROCESS_FILES(reference_ch, input_ch)
}
```

Le processus ne s'exécute qu'une seule fois au lieu de trois fois car le canal `reference_ch` est un canal de file qui s'épuise après la première exécution du processus. Lorsqu'un canal est épuisé, l'ensemble du processus s'arrête, même si d'autres canaux ont encore des éléments.

C'est un modèle courant où vous avez un seul fichier de référence qui doit être réutilisé sur plusieurs échantillons. La solution consiste à convertir le canal de référence en canal de valeur qui peut être réutilisé indéfiniment.

#### Corriger le code

Il existe plusieurs façons de résoudre ce problème selon le nombre de fichiers affectés.

**Option 1** : Vous avez un seul fichier de référence que vous réutilisez beaucoup. Vous pouvez simplement créer un type de canal de valeur, qui peut être utilisé encore et encore. Il y a trois façons de faire cela :

**1a** Utilisez `channel.value()` :

```groovy title="exhausted.nf (corrigé - Option 1a)" hl_lines="2" linenums="21"
workflow {
    reference_ch = channel.value('baseline_reference')  // Le canal de valeur peut être réutilisé
    input_ch = channel.of('sample1', 'sample2', 'sample3')

    PROCESS_FILES(reference_ch, input_ch)
}
```

**1b** Utilisez l'[opérateur](https://www.nextflow.io/docs/latest/reference/operator.html#first) `first()` :

```groovy title="exhausted.nf (corrigé - Option 1b)" hl_lines="2" linenums="21"
workflow {
    reference_ch = channel.of('baseline_reference').first()  // Convertit en canal de valeur
    input_ch = channel.of('sample1', 'sample2', 'sample3')

    PROCESS_FILES(reference_ch, input_ch)
}
```

**1c.** Utilisez l'[opérateur](https://www.nextflow.io/docs/latest/reference/operator.html#collect) `collect()` :

```groovy title="exhausted.nf (corrigé - Option 1c)" hl_lines="2" linenums="21"
workflow {
    reference_ch = channel.of('baseline_reference').collect()  // Convertit en canal de valeur
    input_ch = channel.of('sample1', 'sample2', 'sample3')

    PROCESS_FILES(reference_ch, input_ch)
}
```

**Option 2** : Dans des scénarios plus complexes, peut-être où vous avez plusieurs fichiers de référence pour tous les échantillons dans le canal d'échantillons, vous pouvez utiliser l'opérateur `combine` pour créer un nouveau canal qui combine les deux canaux en tuples :

```groovy title="exhausted.nf (corrigé - Option 2)" hl_lines="4" linenums="21"
workflow {
    reference_ch = channel.of('baseline_reference','other_reference')
    input_ch = channel.of('sample1', 'sample2', 'sample3')
    combined_ch = reference_ch.combine(input_ch)  // Crée un produit cartésien

    PROCESS_FILES(combined_ch)
}
```

L'opérateur `.combine()` génère un produit cartésien des deux canaux, donc chaque élément dans `reference_ch` sera apparié avec chaque élément dans `input_ch`. Cela permet au processus de s'exécuter pour chaque échantillon tout en utilisant la référence.

Cela nécessite d'ajuster l'entrée du processus. Dans notre exemple, le début de la définition du processus devrait être ajusté comme suit :

```groovy title="exhausted.nf (corrigé - Option 2)" hl_lines="5" linenums="1"
#!/usr/bin/env nextflow

process PROCESS_FILES {
    input:
        tuple val(reference), val(sample_name)
```

Cette approche peut ne pas convenir dans toutes les situations.

#### Exécuter le pipeline

Essayez l'une des corrections ci-dessus et exécutez à nouveau le workflow :

```bash
nextflow run exhausted.nf
```

??? success "Sortie de la commande"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `exhausted.nf` [maniac_leavitt] DSL2 - revision: f372a56a7d

    executor >  local (3)
    [80/0779e9] PROCESS_FILES (3) | 3 of 3 ✔
    ```

Vous devriez maintenant voir tous les trois échantillons être traités au lieu d'un seul.

### 2.3. Mauvaise Structure du Contenu des Canaux

Lorsque les workflows atteignent un certain niveau de complexité, il peut être un peu difficile de suivre les structures internes de chaque canal, et les gens génèrent couramment des décalages entre ce que le processus attend et ce que le canal contient réellement. C'est plus subtil que le problème dont nous avons discuté plus tôt, où le nombre de canaux était incorrect. Dans ce cas, vous pouvez avoir le bon nombre de canaux d'entrée, mais la structure interne d'un ou plusieurs de ces canaux ne correspond pas à ce que le processus attend.

#### Exécuter le pipeline

```bash
nextflow run bad_channel_shape.nf
```

??? failure "Sortie de la commande"

    ```console
    Launching `bad_channel_shape.nf` [hopeful_pare] DSL2 - revision: ffd66071a1

    executor >  local (3)
    executor >  local (3)
    [3f/c2dcb3] PROCESS_FILES (3) [  0%] 0 of 3 ✘
    ERROR ~ Error executing process > 'PROCESS_FILES (1)'

    Caused by:
      Missing output file(s) `[sample1, file1.txt]_output.txt` expected by process `PROCESS_FILES (1)`


    Command executed:

      echo "Processing [sample1, file1.txt]" > [sample1, file1.txt]_output.txt

    Command exit status:
      0

    Command output:
      (empty)

    Work dir:
      /workspaces/training/side-quests/debugging/work/d6/1fb69d1d93300bbc9d42f1875b981e

    Tip: when you have fixed the problem you can continue the execution adding the option `-resume` to the run command line

    -- Check '.nextflow.log' file for details
    ```

#### Vérifier le code

Les crochets dans le message d'erreur fournissent l'indice ici - le processus traite le tuple comme une valeur unique, ce qui n'est pas ce que nous voulons. Examinons `bad_channel_shape.nf` :

```groovy title="bad_channel_shape.nf" hl_lines="5 20-22" linenums="1"
#!/usr/bin/env nextflow

process PROCESS_FILES {
    input:
        val sample_name  // Attend une valeur unique, reçoit un tuple

    output:
        path "${sample_name}_output.txt"

    script:
    """
    echo "Processing ${sample_name}" > ${sample_name}_output.txt
    """
}

workflow {

    // Le canal émet des tuples, mais le processus attend des valeurs uniques
    input_ch = channel.of(
      ['sample1', 'file1.txt'],
      ['sample2', 'file2.txt'],
      ['sample3', 'file3.txt']
    )
    PROCESS_FILES(input_ch)
}
```

Vous pouvez voir que nous générons un canal composé de tuples : `['sample1', 'file1.txt']`, mais le processus attend une valeur unique, `val sample_name`. La commande exécutée montre que le processus essaie de créer un fichier nommé `[sample3, file3.txt]_output.txt`, ce qui n'est pas la sortie prévue.

#### Corriger le code

Pour corriger cela, si le processus nécessite les deux entrées, nous pourrions ajuster le processus pour accepter un tuple :

=== "Option 1: Accepter un tuple dans le processus"

    === "Après"

        ```groovy title="bad_channel_shape.nf" hl_lines="5"  linenums="1"
        #!/usr/bin/env nextflow

        process PROCESS_FILES {
            input:
                tuple val(sample_name), val(file_name)  // Corrigé : Accepte le tuple

            output:
                path "${sample_name}_output.txt"

            script:
            """
            echo "Processing ${sample_name}" > ${sample_name}_output.txt
            """
        }

        workflow {

            // Le canal émet des tuples, mais le processus attend des valeurs uniques
            input_ch = channel.of(
              ['sample1', 'file1.txt'],
              ['sample2', 'file2.txt'],
              ['sample3', 'file3.txt']
            )
            PROCESS_FILES(input_ch)
        }
        ```

    === "Avant"

        ```groovy title="bad_channel_shape.nf" hl_lines="5" linenums="1"
        #!/usr/bin/env nextflow

        process PROCESS_FILES {
            input:
                val sample_name  // Attend une valeur unique, reçoit un tuple

            output:
                path "${sample_name}_output.txt"

            script:
            """
            echo "Processing ${sample_name}" > ${sample_name}_output.txt
            """
        }

        workflow {

            // Le canal émet des tuples, mais le processus attend des valeurs uniques
            input_ch = channel.of(
              ['sample1', 'file1.txt'],
              ['sample2', 'file2.txt'],
              ['sample3', 'file3.txt']
            )
            PROCESS_FILES(input_ch)
        }
        ```

=== "Option 2: Extraire le premier élément"

    === "Après"

        ```groovy title="bad_channel_shape.nf" hl_lines="9" linenums="16"
        workflow {

            // Le canal émet des tuples, mais le processus attend des valeurs uniques
            input_ch = channel.of(
              ['sample1', 'file1.txt'],
              ['sample2', 'file2.txt'],
              ['sample3', 'file3.txt']
            )
            PROCESS_FILES(input_ch.map { it[0] })  // Corrigé : Extrait le premier élément
        }
        ```

    === "Avant"

        ```groovy title="bad_channel_shape.nf" hl_lines="9" linenums="16"
        workflow {

            // Le canal émet des tuples, mais le processus attend des valeurs uniques
            input_ch = channel.of(
              ['sample1', 'file1.txt'],
              ['sample2', 'file2.txt'],
              ['sample3', 'file3.txt']
            )
            PROCESS_FILES(input_ch)
        }
        ```

#### Exécuter le pipeline

Choisissez l'une des solutions et réexécutez le workflow :

```bash
nextflow run bad_channel_shape.nf
```

??? success "Sortie de la commande"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `bad_channel_shape.nf` [clever_thompson] DSL2 - revision: 8cbcae3746

    executor >  local (3)
    [bb/80a958] PROCESS_FILES (2) | 3 of 3 ✔
    ```

### 2.4. Techniques de Débogage des Canaux

#### Utilisation de `.view()` pour l'Inspection des Canaux

L'outil de débogage le plus puissant pour les canaux est l'opérateur `.view()`. Avec `.view()`, vous pouvez comprendre la forme de vos canaux à toutes les étapes pour aider au débogage.

#### Exécuter le pipeline

Exécutez `bad_channel_shape_viewed.nf` pour voir cela en action :

```bash
nextflow run bad_channel_shape_viewed.nf
```

??? success "Sortie de la commande"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `bad_channel_shape_viewed.nf` [maniac_poisson] DSL2 - revision: b4f24dc9da

    executor >  local (3)
    [c0/db76b3] PROCESS_FILES (3) [100%] 3 of 3 ✔
    Channel content: [sample1, file1.txt]
    Channel content: [sample2, file2.txt]
    Channel content: [sample3, file3.txt]
    After mapping: sample1
    After mapping: sample2
    After mapping: sample3
    ```

#### Vérifier le code

Examinons `bad_channel_shape_viewed.nf` pour voir comment `.view()` est utilisé :

```groovy title="bad_channel_shape_viewed.nf" linenums="16" hl_lines="9 11"
workflow {

    // Le canal émet des tuples, mais le processus attend des valeurs uniques
    input_ch = channel.of(
      ['sample1', 'file1.txt'],
      ['sample2', 'file2.txt'],
      ['sample3', 'file3.txt']
    )
    .view { "Channel content: $it" }  // Débogage : Affiche le contenu original du canal
    .map { tuple -> tuple[0] }        // Transformation : Extrait le premier élément
    .view { "After mapping: $it" }    // Débogage : Affiche le contenu transformé du canal

    PROCESS_FILES(input_ch)
}
```

#### Corriger le code

Pour vous éviter d'utiliser excessivement les opérations `.view()` à l'avenir pour comprendre le contenu des canaux, il est conseillé d'ajouter quelques commentaires pour aider :

```groovy title="bad_channel_shape_viewed.nf (avec commentaires)" linenums="16" hl_lines="8 9"
workflow {

    // Le canal émet des tuples, mais le processus attend des valeurs uniques
    input_ch = channel.of(
            ['sample1', 'file1.txt'],
            ['sample2', 'file2.txt'],
            ['sample3', 'file3.txt'],
        ) // [sample_name, file_name]
        .map { tuple -> tuple[0] } // sample_name

    PROCESS_FILES(input_ch)
}
```

Cela deviendra plus important à mesure que vos workflows gagneront en complexité et que la structure des canaux deviendra plus opaque.

#### Exécuter le pipeline

```bash
nextflow run bad_channel_shape_viewed.nf
```

??? success "Sortie de la commande"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `bad_channel_shape_viewed.nf` [marvelous_koch] DSL2 - revision: 03e79cdbad

    executor >  local (3)
    [ff/d67cec] PROCESS_FILES (2) | 3 of 3 ✔
    Channel content: [sample1, file1.txt]
    Channel content: [sample2, file2.txt]
    Channel content: [sample3, file3.txt]
    After mapping: sample1
    After mapping: sample2
    After mapping: sample3
    ```

### À retenir

De nombreuses erreurs de structure de canal peuvent être créées avec une syntaxe Nextflow valide. Vous pouvez déboguer les erreurs de structure de canal en comprenant le flux de données, en utilisant des opérateurs `.view()` pour l'inspection, et en reconnaissant les modèles d'erreur comme les crochets indiquant des structures de tuple inattendues.

### Et ensuite ?

Apprenez les erreurs créées par les définitions de processus.

---

## 3. Erreurs de Structure des Processus

La plupart des erreurs que vous rencontrerez liées aux processus seront liées à des erreurs que vous avez faites dans la formation de la commande, ou à des problèmes liés au logiciel sous-jacent. Cela dit, de manière similaire aux problèmes de canal ci-dessus, vous pouvez faire des erreurs dans la définition du processus qui ne sont pas des erreurs de syntaxe, mais qui causeront des erreurs à l'exécution.

### 3.1. Fichiers de Sortie Manquants

Une erreur courante lors de l'écriture de processus est de faire quelque chose qui génère un décalage entre ce que le processus attend et ce qui est généré.

#### Exécuter le pipeline

```bash
nextflow run missing_output.nf
```

??? failure "Sortie de la commande"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `missing_output.nf` [zen_stone] DSL2 - revision: 37ff61f926

    executor >  local (3)
    executor >  local (3)
    [fd/2642e9] process > PROCESS_FILES (2) [ 66%] 2 of 3, failed: 2
    ERROR ~ Error executing process > 'PROCESS_FILES (3)'

    Caused by:
      Missing output file(s) `sample3.txt` expected by process `PROCESS_FILES (3)`


    Command executed:

      echo "Processing sample3" > sample3_output.txt

    Command exit status:
      0

    Command output:
      (empty)

    Work dir:
      /workspaces/training/side-quests/debugging/work/02/9604d49fb8200a74d737c72a6c98ed

    Tip: when you have fixed the problem you can continue the execution adding the option `-resume` to the run command line

    -- Check '.nextflow.log' file for details
    ```

#### Vérifier le code

Le message d'erreur indique que le processus s'attendait à produire un fichier de sortie nommé `sample3.txt`, mais le script crée en réalité `sample3_output.txt`. Examinons la définition du processus dans `missing_output.nf` :

```groovy title="missing_output.nf" linenums="3" hl_lines="6 10"
process PROCESS_FILES {
    input:
    val sample_name

    output:
    path "${sample_name}.txt"  // Attend : sample3.txt

    script:
    """
    echo "Processing ${sample_name}" > ${sample_name}_output.txt  // Crée : sample3_output.txt
    """
}
```

Vous devriez voir qu'il y a un décalage entre le nom du fichier de sortie dans le bloc `output:`, et celui utilisé dans le script. Ce décalage cause l'échec du processus. Si vous rencontrez ce type d'erreur, retournez vérifier que les sorties correspondent entre votre définition de processus et votre bloc de sortie.

Si le problème n'est toujours pas clair, vérifiez le répertoire de travail lui-même pour identifier les fichiers de sortie réels créés :

```bash
❯ ls -h work/02/9604d49fb8200a74d737c72a6c98ed
sample3_output.txt
```

Pour cet exemple, cela nous mettrait en évidence qu'un suffixe `_output` est incorporé dans le nom du fichier de sortie, contrairement à notre définition `output:`.

#### Corriger le code

Corrigez le décalage en rendant le nom de fichier de sortie cohérent :

=== "Après"

    ```groovy title="missing_output.nf" hl_lines="6 10" linenums="3"
    process PROCESS_FILES {
        input:
        val sample_name

        output:
        path "${sample_name}_output.txt"  // Corrigé : Correspond à la sortie du script

        script:
        """
        echo "Processing ${sample_name}" > ${sample_name}_output.txt
        """
    }
    ```

=== "Avant"

    ```groovy title="missing_output.nf" hl_lines="6 10" linenums="3"
    process PROCESS_FILES {
        input:
        val sample_name

        output:
        path "${sample_name}.txt"  // Attend : sample3.txt

        script:
        """
        echo "Processing ${sample_name}" > ${sample_name}_output.txt  // Crée : sample3_output.txt
        """
    }
    ```

#### Exécuter le pipeline

```bash
nextflow run missing_output.nf
```

??? success "Sortie de la commande"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `missing_output.nf` [elated_hamilton] DSL2 - revision: 961938ee2b

    executor >  local (3)
    [16/1c437c] PROCESS_FILES (3) | 3 of 3 ✔
    ```

### 3.2. Logiciel manquant

Une autre classe d'erreurs se produit en raison d'erreurs dans la fourniture de logiciels. `missing_software.nf` est un workflow syntaxiquement valide, mais il dépend d'un logiciel externe pour fournir la commande `cowpy` qu'il utilise.

#### Exécuter le pipeline

```bash
nextflow run missing_software.nf
```

??? failure "Sortie de la commande"

    ```console hl_lines="12 18"
    ERROR ~ Error executing process > 'PROCESS_FILES (3)'

    Caused by:
      Process `PROCESS_FILES (3)` terminated with an error exit status (127)


    Command executed:

      cowpy sample3 > sample3_output.txt

    Command exit status:
      127

    Command output:
      (empty)

    Command error:
      .command.sh: line 2: cowpy: command not found

    Work dir:
      /workspaces/training/side-quests/debugging/work/82/42a5bfb60c9c6ee63ebdbc2d51aa6e

    Tip: you can try to figure out what's wrong by changing to the process work directory and showing the script file named `.command.sh`

    -- Check '.nextflow.log' file for details
    ```

Le processus n'a pas accès à la commande que nous spécifions. Parfois, c'est parce qu'un script est présent dans le répertoire `bin` du workflow, mais n'a pas été rendu exécutable. D'autres fois, c'est parce que le logiciel n'est pas installé dans le conteneur ou l'environnement où le workflow s'exécute.

#### Vérifier le code

Attention à ce code de sortie `127` - il vous dit exactement le problème. Examinons `missing_software.nf` :

```groovy title="missing_software.nf" linenums="3" hl_lines="3"
process PROCESS_FILES {

    container 'community.wave.seqera.io/library/cowpy:1.1.5--3db457ae1977a273'

    input:
    val sample_name

    output:
    path "${sample_name}_output.txt"

    script:
    """
    cowpy ${sample_name} > ${sample_name}_output.txt
    """
}
```

#### Corriger le code

Nous avons été un peu malhonnêtes ici, et il n'y a en fait rien de mal avec le code. Nous devons simplement spécifier la configuration nécessaire pour exécuter le processus de manière à ce qu'il ait accès à la commande en question. Dans ce cas, le processus a une définition de conteneur, donc tout ce que nous devons faire est d'exécuter le workflow avec Docker activé.

#### Exécuter le pipeline

Nous avons configuré un profil Docker pour vous dans `nextflow.config`, vous pouvez donc exécuter le workflow avec :

```bash
nextflow run missing_software.nf -profile docker
```

??? success "Sortie de la commande"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `missing_software.nf` [awesome_stonebraker] DSL2 - revision: 0296d12839

    executor >  local (3)
    [38/ab20d1] PROCESS_FILES (1) | 3 of 3 ✔
    ```

!!! note

    Pour en savoir plus sur l'utilisation des conteneurs par Nextflow, voir [Hello Nextflow](../hello_nextflow/05_hello_containers.md)

### 3.3. Mauvaise configuration des ressources

En production, vous configurerez les ressources de vos processus. Par exemple, `memory` définit la quantité maximale de mémoire disponible pour votre processus, et si le processus dépasse cette limite, votre ordonnanceur tuera typiquement le processus et renverra un code de sortie `137`. Nous ne pouvons pas le démontrer ici car nous utilisons l'exécuteur `local`, mais nous pouvons montrer quelque chose de similaire avec `time`.

#### Exécuter le pipeline

`bad_resources.nf` a une configuration de processus avec une limite de temps irréaliste de 1 milliseconde :

```bash
nextflow run bad_resources.nf -profile docker
```

??? failure "Sortie de la commande"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `bad_resources.nf` [disturbed_elion] DSL2 - revision: 27d2066e86

    executor >  local (3)
    [c0/ded8e1] PROCESS_FILES (3) | 0 of 3 ✘
    ERROR ~ Error executing process > 'PROCESS_FILES (2)'

    Caused by:
      Process exceeded running time limit (1ms)

    Command executed:

      cowpy sample2 > sample2_output.txt

    Command exit status:
      -

    Command output:
      (empty)

    Work dir:
      /workspaces/training/side-quests/debugging/work/53/f0a4cc56d6b3dc2a6754ff326f1349

    Container:
      community.wave.seqera.io/library/cowpy:1.1.5--3db457ae1977a273

    Tip: you can replicate the issue by changing to the process work dir and entering the command `bash .command.run`

     -- Check '.nextflow.log' file for details
    ```

#### Vérifier le code

Examinons `bad_resources.nf` :

```groovy title="bad_resources.nf" linenums="3" hl_lines="3"
process PROCESS_FILES {

    time '1 ms'  // ERREUR : Limite de temps irréaliste

    input:
    val sample_name

    output:
    path "${sample_name}_output.txt"

    script:
    """
    sleep 1  // Prend 1 seconde, mais la limite de temps est de 1ms
    cowpy ${sample_name} > ${sample_name}_output.txt
    """
}
```

Nous savons que le processus prendra plus d'une seconde (nous avons ajouté un sleep pour nous en assurer), mais le processus est configuré pour expirer après 1 milliseconde. Quelqu'un a été un peu irréaliste avec sa configuration !

#### Corriger le code

Augmentez la limite de temps à une valeur réaliste :

=== "Après"

    ```groovy title="bad_resources.nf" hl_lines="3" linenums="3"
    process PROCESS_FILES {

        time '100 s'  // Corrigé : Limite de temps réaliste

        input:
        val sample_name

        output:
        path "${sample_name}_output.txt"

        script:
        """
        sleep 1
        cowpy ${sample_name} > ${sample_name}_output.txt
        """
    }
    ```

=== "Avant"

    ```groovy title="bad_resources.nf" hl_lines="3" linenums="3"
    process PROCESS_FILES {

        time '1 ms'  // ERREUR : Limite de temps irréaliste

        input:
        val sample_name

        output:
        path "${sample_name}_output.txt"

        script:
        """
        sleep 1  // Prend 1 seconde, mais la limite de temps est de 1ms
        cowpy ${sample_name} > ${sample_name}_output.txt
        """
    }
    ```

#### Exécuter le pipeline

```bash
nextflow run bad_resources.nf -profile docker
```

??? success "Sortie de la commande"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `bad_resources.nf` [friendly_mcclintock] DSL2 - revision: 381567d2c1

    executor >  local (3)
    [c2/9b4c41] PROCESS_FILES (3) | 3 of 3 ✔
    ```

Si vous prenez soin de lire vos messages d'erreur, des échecs comme celui-ci ne devraient pas vous déconcerter trop longtemps. Mais assurez-vous de comprendre les besoins en ressources des commandes que vous exécutez afin de configurer vos directives de ressources de manière appropriée.

### 3.4. Techniques de Débogage des Processus

Lorsque des processus échouent ou se comportent de manière inattendue, vous avez besoin de techniques systématiques pour investiguer ce qui s'est mal passé. Le répertoire de travail contient toutes les informations dont vous avez besoin pour déboguer l'exécution du processus.

#### Utilisation de l'Inspection du Répertoire de Travail

L'outil de débogage le plus puissant pour les processus est l'examen du répertoire de travail. Lorsqu'un processus échoue, Nextflow crée un répertoire de travail pour cette exécution spécifique du processus contenant tous les fichiers nécessaires pour comprendre ce qui s'est passé.

#### Exécuter le pipeline

Utilisons l'exemple `missing_output.nf` de tout à l'heure pour démontrer l'inspection du répertoire de travail (régénérez un décalage de nommage de sortie si nécessaire) :

```bash
nextflow run missing_output.nf
```

??? failure "Sortie de la commande"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `missing_output.nf` [irreverent_payne] DSL2 - revision: 3d5117f7e2

    executor >  local (3)
    [5d/d544a4] PROCESS_FILES (2) | 0 of 3 ✘
    ERROR ~ Error executing process > 'PROCESS_FILES (1)'

    Caused by:
      Missing output file(s) `sample1.txt` expected by process `PROCESS_FILES (1)`

    Command executed:

      echo "Processing sample1" > sample1_output.txt

    Command exit status:
      0

    Command output:
      (empty)

    Work dir:
      /workspaces/training/side-quests/debugging/work/1e/2011154d0b0f001cd383d7364b5244

    Tip: you can replicate the issue by changing to the process work dir and entering the command `bash .command.run`

     -- Check '.nextflow.log' file for details
    ```

#### Vérifier le répertoire de travail

Lorsque vous obtenez cette erreur, le répertoire de travail contient toutes les informations de débogage. Trouvez le chemin du répertoire de travail dans le message d'erreur et examinez son contenu :

```bash
# Trouvez le répertoire de travail dans le message d'erreur
ls work/02/9604d49fb8200a74d737c72a6c98ed/
```

Vous pouvez ensuite examiner les fichiers clés :

##### Vérifier le Script de Commande

Le fichier `.command.sh` montre exactement quelle commande a été exécutée :

```bash
# Voir la commande exécutée
cat work/02/9604d49fb8200a74d737c72a6c98ed/.command.sh
```

Cela révèle :

- **Substitution de variables** : Si les variables Nextflow ont été correctement développées
- **Chemins de fichiers** : Si les fichiers d'entrée ont été correctement localisés
- **Structure de la commande** : Si la syntaxe du script est correcte

Problèmes courants à rechercher :

- **Guillemets manquants** : Les variables contenant des espaces nécessitent des guillemets appropriés
- **Mauvais chemins de fichiers** : Fichiers d'entrée qui n'existent pas ou sont au mauvais endroit
- **Noms de variables incorrects** : Fautes de frappe dans les références de variables
- **Configuration d'environnement manquante** : Commandes qui dépendent d'environnements spécifiques

##### Vérifier la Sortie d'Erreur

Le fichier `.command.err` contient les messages d'erreur réels :

```bash
# Voir la sortie d'erreur
cat work/02/9604d49fb8200a74d737c72a6c98ed/.command.err
```

Ce fichier montrera :

- **Codes de sortie** : 127 (commande introuvable), 137 (tué), etc.
- **Erreurs de permission** : Problèmes d'accès aux fichiers
- **Erreurs logicielles** : Messages d'erreur spécifiques à l'application
- **Erreurs de ressources** : Limite de mémoire/temps dépassée

##### Vérifier la Sortie Standard

Le fichier `.command.out` montre ce que votre commande a produit :

```bash
# Voir la sortie standard
cat work/02/9604d49fb8200a74d737c72a6c98ed/.command.out
```

Cela aide à vérifier :

- **Sortie attendue** : Si la commande a produit les bons résultats
- **Exécution partielle** : Si la commande a démarré mais a échoué en cours de route
- **Informations de débogage** : Toute sortie de diagnostic de votre script

##### Vérifier le Code de Sortie

Le fichier `.exitcode` contient le code de sortie du processus :

```bash
# Voir le code de sortie
cat work/*/*/.exitcode
```

Codes de sortie courants et leurs significations :

- **Code de sortie 127** : Commande introuvable - vérifiez l'installation du logiciel
- **Code de sortie 137** : Processus tué - vérifiez les limites de mémoire/temps

##### Vérifier l'Existence des Fichiers

Lorsque les processus échouent en raison de fichiers de sortie manquants, vérifiez quels fichiers ont réellement été créés :

```bash
# Lister tous les fichiers dans le répertoire de travail
ls -la work/02/9604d49fb8200a74d737c72a6c98ed/
```

Cela aide à identifier :

- **Décalages de noms de fichiers** : Fichiers de sortie avec des noms différents de ceux attendus
- **Problèmes de permissions** : Fichiers qui n'ont pas pu être créés
- **Problèmes de chemin** : Fichiers créés dans de mauvais répertoires

Dans notre exemple précédent, cela nous a confirmé que tandis que notre `sample3.txt` attendu n'était pas présent, `sample3_output.txt` l'était :

```bash
❯ ls -h work/02/9604d49fb8200a74d737c72a6c98ed
sample3_output.txt
```

### À retenir

Le débogage des processus nécessite l'examen des répertoires de travail pour comprendre ce qui s'est mal passé. Les fichiers clés incluent `.command.sh` (le script exécuté), `.command.err` (les messages d'erreur) et `.command.out` (la sortie standard). Les codes de sortie comme 127 (commande introuvable) et 137 (processus tué) fournissent des indices diagnostiques immédiats sur le type de défaillance.

### Et ensuite ?

Découvrez les outils de débogage intégrés de Nextflow et les approches systématiques de dépannage.

---

## 4. Outils de Débogage Intégrés et Techniques Avancées

Nextflow fournit plusieurs outils intégrés puissants pour le débogage et l'analyse de l'exécution des workflows. Ces outils vous aident à comprendre ce qui s'est mal passé, où cela s'est mal passé et comment le corriger efficacement.

### 4.1. Sortie de Processus en Temps Réel

Parfois, vous devez voir ce qui se passe à l'intérieur des processus en cours d'exécution. Vous pouvez activer la sortie de processus en temps réel, qui vous montre exactement ce que chaque tâche fait pendant son exécution.

#### Exécuter le pipeline

`bad_channel_shape_viewed.nf` de nos exemples précédents a affiché le contenu des canaux en utilisant `.view()`, mais nous pouvons aussi utiliser la directive `debug` pour afficher des variables depuis l'intérieur du processus lui-même, ce que nous démontrons dans `bad_channel_shape_viewed_debug.nf`. Exécutez le workflow :

```bash
nextflow run bad_channel_shape_viewed_debug.nf
```

??? success "Sortie de la commande"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `bad_channel_shape_viewed_debug.nf` [agitated_crick] DSL2 - revision: ea3676d9ec

    executor >  local (3)
    [c6/2dac51] process > PROCESS_FILES (3) [100%] 3 of 3 ✔
    Channel content: [sample1, file1.txt]
    Channel content: [sample2, file2.txt]
    Channel content: [sample3, file3.txt]
    After mapping: sample1
    After mapping: sample2
    After mapping: sample3
    Sample name inside process is sample2

    Sample name inside process is sample1

    Sample name inside process is sample3
    ```

#### Vérifier le code

Examinons `bad_channel_shape_viewed_debug.nf` pour voir comment la directive `debug` fonctionne :

```groovy title="bad_channel_shape_viewed_debug.nf" linenums="3" hl_lines="2"
process PROCESS_FILES {
    debug true  // Enable real-time output

    input:
    val sample_name

    output:
    path "${sample_name}_output.txt"

    script:
    """
    echo "Sample name inside process is ${sample_name}"
    echo "Processing ${sample_name}" > ${sample_name}_output.txt
    """
}
```

La directive `debug` peut être un moyen rapide et pratique pour comprendre l'environnement d'un processus.

### 4.2. Mode Preview

Parfois, vous voulez détecter des problèmes avant que des processus ne s'exécutent. Nextflow fournit un flag pour ce type de débogage proactif : `-preview`.

#### Exécuter le pipeline

Le mode preview vous permet de tester la logique du workflow sans exécuter de commandes. Cela peut être très utile pour vérifier rapidement la structure de votre workflow et s'assurer que les processus sont correctement connectés sans exécuter de commandes réelles.

!!! note

    Si vous avez corrigé `bad_syntax.nf` plus tôt, réintroduisez l'erreur de syntaxe en supprimant l'accolade fermante après le bloc script avant d'exécuter cette commande.

Exécutez cette commande :

```bash
nextflow run bad_syntax.nf -preview
```

??? failure "Sortie de la commande"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `bad_syntax.nf` [magical_mercator] DSL2 - revision: 550b9a8873

    Error bad_syntax.nf:24:1: Unexpected input: '<EOF>'

    ERROR ~ Script compilation failed

     -- Check '.nextflow.log' file for details
    ```

Le mode preview est particulièrement utile pour détecter les erreurs de syntaxe tôt sans exécuter de processus. Il valide la structure du workflow et les connexions des processus avant l'exécution.

### 4.3. Exécution Stub pour les Tests de Logique

Parfois, les erreurs sont difficiles à déboguer parce que les commandes prennent trop de temps, nécessitent un logiciel spécial ou échouent pour des raisons complexes. L'exécution stub vous permet de tester la logique du workflow sans exécuter les commandes réelles.

#### Exécuter le pipeline

Lorsque vous développez un processus Nextflow, vous pouvez utiliser la directive `stub` pour définir des commandes 'fictives' qui génèrent des sorties de la forme correcte sans exécuter la commande réelle. Cette approche est particulièrement précieuse lorsque vous voulez vérifier que la logique de votre workflow est correcte avant de traiter les complexités du logiciel réel.

Par exemple, vous souvenez-vous de notre `missing_software.nf` de tout à l'heure ? Celui où nous avions un logiciel manquant qui empêchait le workflow de s'exécuter jusqu'à ce que nous ajoutions `-profile docker` ? `missing_software_with_stub.nf` est un workflow très similaire. Si nous l'exécutons de la même manière, nous générerons la même erreur :

```bash
nextflow run missing_software_with_stub.nf
```

??? failure "Sortie de la commande"

    ```console hl_lines="12 18"
    ERROR ~ Error executing process > 'PROCESS_FILES (3)'

    Caused by:
      Process `PROCESS_FILES (3)` terminated with an error exit status (127)


    Command executed:

      cowpy sample3 > sample3_output.txt

    Command exit status:
      127

    Command output:
      (empty)

    Command error:
      .command.sh: line 2: cowpy: command not found

    Work dir:
      /workspaces/training/side-quests/debugging/work/82/42a5bfb60c9c6ee63ebdbc2d51aa6e

    Tip: you can try to figure out what's wrong by changing to the process work directory and showing the script file named `.command.sh`

    -- Check '.nextflow.log' file for details
    ```

Cependant, ce workflow ne produira pas d'erreurs si nous l'exécutons avec `-stub-run`, même sans le profil `docker` :

```bash
nextflow run missing_software_with_stub.nf -stub-run
```

??? success "Sortie de la commande"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `missing_software_with_stub.nf` [astonishing_shockley] DSL2 - revision: f1f4f05d7d

    executor >  local (3)
    [b5/2517a3] PROCESS_FILES (3) | 3 of 3 ✔
    ```

#### Vérifier le code

Examinons `missing_software_with_stub.nf` :

```groovy title="missing_software.nf (with stub)" hl_lines="16-19" linenums="3"
process PROCESS_FILES {

    container 'community.wave.seqera.io/library/cowpy:1.1.5--3db457ae1977a273'

    input:
    val sample_name

    output:
    path "${sample_name}_output.txt"

    script:
    """
    cowpy ${sample_name} > ${sample_name}_output.txt
    """

    stub:
    """
    touch ${sample_name}_output.txt
    """
}
```

Par rapport à `missing_software.nf`, ce processus a une directive `stub:` spécifiant une commande à utiliser à la place de celle spécifiée dans `script:`, dans le cas où Nextflow est exécuté en mode stub.

La commande `touch` que nous utilisons ici ne dépend d'aucun logiciel ou entrées appropriées, et s'exécutera dans toutes les situations, nous permettant de déboguer la logique du workflow sans nous soucier des détails internes du processus.

**L'exécution stub aide à déboguer :**

- La structure des canaux et le flux de données
- Les connexions et dépendances des processus
- La propagation des paramètres
- La logique du workflow sans dépendances logicielles

### 4.4. Approche Systématique de Débogage

Maintenant que vous avez appris des techniques de débogage individuelles - des fichiers de trace et répertoires de travail au mode preview, l'exécution stub et la surveillance des ressources - combinons-les en une méthodologie systématique. Avoir une approche structurée vous empêche d'être submergé par des erreurs complexes et garantit que vous ne manquerez pas d'indices importants.

Cette méthodologie combine tous les outils que nous avons couverts en un workflow efficace :

**Méthode de Débogage en Quatre Phases :**

**Phase 1 : Résolution des Erreurs de Syntaxe (5 minutes)**

1. Vérifiez les soulignements rouges dans VSCode ou votre IDE
2. Exécutez `nextflow run workflow.nf -preview` pour identifier les problèmes de syntaxe
3. Corrigez toutes les erreurs de syntaxe (accolades manquantes, virgules finales, etc.)
4. Assurez-vous que le workflow s'analyse avec succès avant de continuer

**Phase 2 : Évaluation Rapide (5 minutes)**

1. Lisez attentivement les messages d'erreur d'exécution
2. Vérifiez s'il s'agit d'une erreur d'exécution, de logique ou de ressources
3. Utilisez le mode preview pour tester la logique de base du workflow

**Phase 3 : Investigation Détaillée (15-30 minutes)**

1. Trouvez le répertoire de travail de la tâche échouée
2. Examinez les fichiers de log
3. Ajoutez des opérateurs `.view()` pour inspecter les canaux
4. Utilisez `-stub-run` pour tester la logique du workflow sans exécution

**Phase 4 : Corriger et Valider (15 minutes)**

1. Faites des corrections minimales et ciblées
2. Testez avec resume : `nextflow run workflow.nf -resume`
3. Vérifiez l'exécution complète du workflow

!!! tip "Utilisation de Resume pour un Débogage Efficace"

    Une fois que vous avez identifié un problème, vous avez besoin d'un moyen efficace de tester vos corrections sans perdre de temps à réexécuter les parties réussies de votre workflow. La fonctionnalité `-resume` de Nextflow est inestimable pour le débogage.

    Vous aurez rencontré `-resume` si vous avez travaillé sur [Hello Nextflow](../hello_nextflow/), et il est important de bien l'utiliser lors du débogage pour vous éviter d'attendre pendant que les processus avant votre processus problématique s'exécutent.

    **Stratégie de débogage avec resume :**

    1. Exécutez le workflow jusqu'à l'échec
    2. Examinez le répertoire de travail de la tâche échouée
    3. Corrigez le problème spécifique
    4. Reprenez pour tester uniquement la correction
    5. Répétez jusqu'à ce que le workflow se termine

#### Profil de Configuration de Débogage

Pour rendre cette approche systématique encore plus efficace, vous pouvez créer une configuration de débogage dédiée qui active automatiquement tous les outils dont vous avez besoin :

```groovy title="nextflow.config (debug profile)" linenums="1"
profiles {
    debug {
        process {
            debug = true
            cleanup = false

            // Conservative resources for debugging
            maxForks = 1
            memory = '2.GB'
            cpus = 1
        }
    }
}
```

Ensuite, vous pouvez exécuter le pipeline avec ce profil activé :

```bash
nextflow run workflow.nf -profile debug
```

Ce profil active la sortie en temps réel, préserve les répertoires de travail et limite la parallélisation pour un débogage plus facile.

### 4.5. Exercice Pratique de Débogage

Il est maintenant temps de mettre en pratique l'approche systématique de débogage. Le workflow `buggy_workflow.nf` contient plusieurs erreurs courantes qui représentent les types de problèmes que vous rencontrerez dans le développement réel.

!!! exercise

    Utilisez l'approche systématique de débogage pour identifier et corriger toutes les erreurs dans `buggy_workflow.nf`. Ce workflow tente de traiter des données d'échantillons à partir d'un fichier CSV mais contient plusieurs bugs intentionnels représentant des scénarios de débogage courants.

    Commencez par exécuter le workflow pour voir la première erreur :

    ```bash
    nextflow run buggy_workflow.nf
    ```

    ??? failure "Sortie de la commande"

        ```console
        N E X T F L O W   ~  version 25.10.2

        Launching `buggy_workflow.nf` [wise_ramanujan] DSL2 - revision: d51a8e83fd

        ERROR ~ Range [11, 12) out of bounds for length 11

         -- Check '.nextflow.log' file for details
        ```

        Cette erreur cryptique indique un problème d'analyse autour des lignes 11-12 dans le bloc `params{}`. L'analyseur v2 détecte les problèmes structurels tôt.

    Appliquez la méthode de débogage en quatre phases que vous avez apprise :

    **Phase 1 : Résolution des Erreurs de Syntaxe**
    - Vérifiez les soulignements rouges dans VSCode ou votre IDE
    - Exécutez `nextflow run workflow.nf -preview` pour identifier les problèmes de syntaxe
    - Corrigez toutes les erreurs de syntaxe (accolades manquantes, virgules finales, etc.)
    - Assurez-vous que le workflow s'analyse avec succès avant de continuer

    **Phase 2 : Évaluation Rapide**
    - Lisez attentivement les messages d'erreur d'exécution
    - Identifiez si les erreurs sont liées à l'exécution, à la logique ou aux ressources
    - Utilisez le mode `-preview` pour tester la logique de base du workflow

    **Phase 3 : Investigation Détaillée**
    - Examinez les répertoires de travail des tâches échouées
    - Ajoutez des opérateurs `.view()` pour inspecter les canaux
    - Vérifiez les fichiers de log dans les répertoires de travail
    - Utilisez `-stub-run` pour tester la logique du workflow sans exécution

    **Phase 4 : Corriger et Valider**
    - Faites des corrections ciblées
    - Utilisez `-resume` pour tester les corrections efficacement
    - Vérifiez l'exécution complète du workflow

    **Outils de Débogage à Votre Disposition :**
    ```bash
    # Mode preview pour la vérification de syntaxe
    nextflow run buggy_workflow.nf -preview

    # Profil debug pour une sortie détaillée
    nextflow run buggy_workflow.nf -profile debug

    # Exécution stub pour les tests de logique
    nextflow run buggy_workflow.nf -stub-run

    # Resume après corrections
    nextflow run buggy_workflow.nf -resume
    ```

    ??? solution
        Le `buggy_workflow.nf` contient 9 ou 10 erreurs distinctes (selon la façon dont vous comptez) couvrant toutes les catégories principales de débogage. Voici une analyse systématique de chaque erreur et comment la corriger

        Commençons par les erreurs de syntaxe :

        **Erreur 1 : Erreur de Syntaxe - Virgule Finale**
        ```groovy linenums="21"
        output:
            path "${sample_id}_result.txt",  // ERREUR : Virgule finale
        ```
        **Correction :** Supprimez la virgule finale
        ```groovy linenums="21"
        output:
            path "${sample_id}_result.txt"
        ```

        **Erreur 2 : Erreur de Syntaxe - Accolade Fermante Manquante**
        ```groovy linenums="24"
        script:
        """
        echo "Processing: ${sample}"
        cat ${input_file} > ${sample}_result.txt
        """
        // ERREUR : Accolade fermante manquante pour le processus processFiles
        ```
        **Correction :** Ajoutez l'accolade fermante manquante
        ```groovy linenums="29"
        """
        echo "Processing: ${sample_id}"
        cat ${input_file} > ${sample_id}_result.txt
        """
        }  // Ajout de l'accolade fermante manquante
        ```

        **Erreur 3 : Erreur de Nom de Variable**
        ```groovy linenums="26"
        echo "Processing: ${sample}"     // ERREUR : devrait être sample_id
        cat ${input_file} > ${sample}_result.txt  // ERREUR : devrait être sample_id
        ```
        **Correction :** Utilisez le bon nom de variable d'entrée
        ```groovy linenums="26"
        echo "Processing: ${sample_id}"
        cat ${input_file} > ${sample_id}_result.txt
        ```

        **Erreur 4 : Erreur de Variable Non Définie**
        ```groovy linenums="87"
        heavy_ch = heavyProcess(sample_ids)  // ERREUR : sample_ids non défini
        ```
        **Correction :** Utilisez le bon canal et extrayez les IDs d'échantillons
        ```groovy linenums="87"
        heavy_ch = heavyProcess(input_ch)
        ```

        À ce stade, le workflow s'exécutera, mais nous aurons toujours des erreurs (ex. `Path value cannot be null` dans `processFiles`), causées par une mauvaise structure de canal.

        **Erreur 5 : Erreur de Structure de Canal - Mauvaise Sortie Map**
        ```groovy linenums="83"
        .map { row -> row.sample_id }  // ERREUR : processFiles attend un tuple
        ```
        **Correction :** Retournez la structure de tuple que processFiles attend
        ```groovy linenums="83"
        .map { row -> [row.sample_id, file(row.fastq_path)] }
        ```

        Mais cela cassera notre correction pour l'exécution de `heavyProcess()` ci-dessus, nous devrons donc utiliser un map pour ne passer que les IDs d'échantillons à ce processus :

        **Erreur 6 : Mauvaise structure de canal pour heavyProcess**
        ```groovy linenums="87"
        heavy_ch = heavyProcess(input_ch)  // ERREUR : input_ch a maintenant 2 éléments par émission - heavyProcess n'en a besoin que d'un (le premier)
        ```
        **Correction :** Utilisez le bon canal et extrayez les IDs d'échantillons
        ```groovy linenums="87"
        heavy_ch = heavyProcess(input_ch.map{it[0]})
        ```

        Maintenant nous avançons un peu plus mais recevons une erreur à propos de `No such variable: i`, parce que nous n'avons pas échappé une variable Bash.

        **Erreur 7 : Erreur d'Échappement de Variable Bash**
        ```groovy linenums="48"
        echo "Heavy computation $i for ${sample_id}"  // ERREUR : $i non échappé
        ```
        **Correction :** Échappez la variable bash
        ```groovy linenums="48"
        echo "Heavy computation \${i} for ${sample_id}"
        ```

        Maintenant nous obtenons `Process exceeded running time limit (1ms)`, donc nous corrigeons la limite de temps d'exécution pour le processus concerné :

        **Erreur 8 : Erreur de Configuration de Ressources**
        ```groovy linenums="36"
        time '1 ms'  // ERREUR : Limite de temps irréaliste
        ```
        **Correction :** Augmentez à une limite de temps réaliste
        ```groovy linenums="36"
        time '100 s'
        ```

        Ensuite, nous avons une erreur `Missing output file(s)` à résoudre :

        **Erreur 9 : Décalage de Nom de Fichier de Sortie**
        ```groovy linenums="49"
        done > ${sample_id}.txt  // ERREUR : Mauvais nom de fichier, devrait correspondre à la déclaration de sortie
        ```
        **Correction :** Faites correspondre la déclaration de sortie
        ```groovy linenums="49"
        done > ${sample_id}_heavy.txt
        ```

        Les deux premiers processus se sont exécutés, mais pas le troisième.

        **Erreur 10 : Décalage de Nom de Fichier de Sortie**
        ```groovy linenums="88"
        file_ch = channel.fromPath("*.txt") // Erreur : tentative de prendre l'entrée du pwd plutôt que d'un processus
        handleFiles(file_ch)
        ```
        **Correction :** Prenez la sortie du processus précédent
        ```groovy linenums="88"
        handleFiles(heavyProcess.out)
        ```

        Avec cela, l'ensemble du workflow devrait s'exécuter.

        **Workflow Corrigé Complet :**
        ```groovy linenums="1"
        #!/usr/bin/env nextflow

        /*
        * Buggy workflow for debugging exercises
        * This workflow contains several intentional bugs for learning purposes
        */

        params{
            // Parameters with missing validation
            input: Path = 'data/sample_data.csv'
            output: String = 'results'
        }

        /*
        * Process with input/output mismatch
        */
        process processFiles {
            publishDir "${params.output}/processed", mode: 'copy'

            input:
                tuple val(sample_id), path(input_file)

            output:
                path "${sample_id}_result.txt"

            script:
            """
            echo "Processing: ${sample_id}"
            cat ${input_file} > ${sample_id}_result.txt
            """
        }

        /*
        * Process with resource issues
        */
        process heavyProcess {
            publishDir "${params.output}/heavy", mode: 'copy'

            time '100 s'

            input:
                val sample_id

            output:
                path "${sample_id}_heavy.txt"

            script:
            """
            # Simulate heavy computation
            for i in {1..1000000}; do
                echo "Heavy computation \$i for ${sample_id}"
            done > ${sample_id}_heavy.txt
            """
        }

        /*
        * Process with file handling issues
        */
        process handleFiles {
            publishDir "${params.output}/files", mode: 'copy'

            input:
                path input_file

            output:
                path "processed_${input_file}"

            script:
            """
            if [ -f "${input_file}" ]; then
                cp ${input_file} processed_${input_file}
            fi
            """
        }

        /*
        * Main workflow with channel issues
        */
        workflow {

            // Channel with incorrect usage
            input_ch = channel
                .fromPath(params.input)
                .splitCsv(header: true)
                .map { row -> [row.sample_id, file(row.fastq_path)] }

            processed_ch = processFiles(input_ch)

            heavy_ch = heavyProcess(input_ch.map{it[0]})

            handleFiles(heavyProcess.out)
        }
        ```

**Catégories d'Erreurs Couvertes :**

- **Erreurs de syntaxe** : Accolades manquantes, virgules finales, variables non définies
- **Erreurs de structure de canal** : Mauvaises formes de données, canaux non définis
- **Erreurs de processus** : Décalages de fichiers de sortie, échappement de variables
- **Erreurs de ressources** : Limites de temps irréalistes

**Leçons Clés de Débogage :**

1. **Lisez attentivement les messages d'erreur** - ils pointent souvent directement vers le problème
2. **Utilisez des approches systématiques** - corrigez une erreur à la fois et testez avec `-resume`
3. **Comprenez le flux de données** - les erreurs de structure de canal sont souvent les plus subtiles
4. **Vérifiez les répertoires de travail** - quand les processus échouent, les logs vous disent exactement ce qui s'est mal passé

---

## Résumé

Dans cette quête annexe, vous avez appris un ensemble de techniques systématiques pour déboguer les workflows Nextflow.
L'application de ces techniques dans votre propre travail vous permettra de passer moins de temps à lutter contre votre ordinateur, de résoudre les problèmes plus rapidement et de vous protéger contre les problèmes futurs.

### Modèles clés

**1. Comment identifier et corriger les erreurs de syntaxe** :

- Interprétation des messages d'erreur de Nextflow et localisation des problèmes
- Erreurs de syntaxe courantes : accolades manquantes, mots-clés incorrects, variables non définies
- Distinction entre les variables Nextflow (Groovy) et Bash
- Utilisation des fonctionnalités de l'extension VS Code pour la détection précoce des erreurs

```groovy
// Accolade manquante - cherchez les soulignements rouges dans l'IDE
process FOO {
    script:
    """
    echo "hello"
    """
// } <-- manquante !

// Mauvais mot-clé
inputs:  // Devrait être 'input:'

// Variable non définie - échappez avec un backslash pour les variables Bash
echo "${undefined_var}"      // Variable Nextflow (erreur si non définie)
echo "\${bash_var}"          // Variable Bash (échappée)
```

**2. Comment déboguer les problèmes de structure de canal** :

- Compréhension de la cardinalité des canaux et des problèmes d'épuisement
- Débogage des décalages de structure de contenu des canaux
- Utilisation des opérateurs `.view()` pour l'inspection des canaux
- Reconnaissance des modèles d'erreur comme les crochets dans la sortie

```groovy
// Inspecter le contenu du canal
my_channel.view { "Content: $it" }

// Convertir un canal de file en canal de valeur (empêche l'épuisement)
reference_ch = channel.value('ref.fa')
// ou
reference_ch = channel.of('ref.fa').first()
```

**3. Comment résoudre les problèmes d'exécution des processus** :

- Diagnostic des erreurs de fichiers de sortie manquants
- Compréhension des codes de sortie (127 pour logiciel manquant, 137 pour problèmes de mémoire)
- Investigation des répertoires de travail et fichiers de commande
- Configuration appropriée des ressources

```bash
# Vérifier ce qui a été réellement exécuté
cat work/ab/cdef12/.command.sh

# Vérifier la sortie d'erreur
cat work/ab/cdef12/.command.err

# Code de sortie 127 = commande introuvable
# Code de sortie 137 = tué (limite de mémoire/temps)
```

**4. Comment utiliser les outils de débogage intégrés de Nextflow** :

- Exploitation du mode preview et du débogage en temps réel
- Implémentation de l'exécution stub pour les tests de logique
- Application de resume pour des cycles de débogage efficaces
- Suivi d'une méthodologie de débogage systématique en quatre phases

!!! tip "Référence Rapide de Débogage"

    **Erreurs de syntaxe ?** → Vérifiez les avertissements VSCode, exécutez `nextflow run workflow.nf -preview`

    **Problèmes de canal ?** → Utilisez `.view()` pour inspecter le contenu : `my_channel.view()`

    **Échecs de processus ?** → Vérifiez les fichiers du répertoire de travail :

    - `.command.sh` - le script exécuté
    - `.command.err` - les messages d'erreur
    - `.exitcode` - le statut de sortie (127 = commande introuvable, 137 = tué)

    **Comportement mystérieux ?** → Exécutez avec `-stub-run` pour tester la logique du workflow

    **Corrections effectuées ?** → Utilisez `-resume` pour gagner du temps lors des tests : `nextflow run workflow.nf -resume`

---

### Ressources supplémentaires

- [Guide de dépannage Nextflow](https://www.nextflow.io/docs/latest/troubleshooting.html) : Documentation officielle de dépannage
- [Comprendre les canaux Nextflow](https://www.nextflow.io/docs/latest/channel.html) : Exploration approfondie des types et comportements de canaux
- [Référence des directives de processus](https://www.nextflow.io/docs/latest/process.html#directives) : Toutes les options de configuration de processus disponibles
- [nf-test](https://www.nf-test.com/) : Framework de test pour les pipelines Nextflow
- [Communauté Slack Nextflow](https://www.nextflow.io/slack-invite.html) : Obtenir de l'aide de la communauté

Pour les workflows de production, considérez :

- La mise en place de [Seqera Platform](https://seqera.io/platform/) pour la surveillance et le débogage à grande échelle
- L'utilisation de [Wave containers](https://seqera.io/wave/) pour des environnements logiciels reproductibles

**N'oubliez pas :** Le débogage efficace est une compétence qui s'améliore avec la pratique. La méthodologie systématique et la boîte à outils complète que vous avez acquises ici vous serviront bien tout au long de votre parcours de développement Nextflow.

---

## Et ensuite ?

Retournez au [menu des Quêtes Annexes](./index.md) ou cliquez sur le bouton en bas à droite de la page pour passer au sujet suivant de la liste.
