# Tests avec nf-test

Pouvoir tester systématiquement que chaque partie de votre workflow fait ce qu'elle est censée faire est essentiel pour la reproductibilité et la maintenance à long terme, et peut être d'une grande aide durant le processus de développement.

Prenons un moment pour parler de l'importance des tests. Si vous développez un workflow, l'une des premières choses que vous ferez sera de prendre des données de test dont vous savez qu'elles sont valides et devraient produire un résultat. Vous ajoutez le premier processus au pipeline et le connectez à vos entrées pour le faire fonctionner. Ensuite, pour vérifier que tout fonctionne, vous l'exécutez sur les données de test. En supposant que cela fonctionne, vous passez au processus suivant et exécutez à nouveau les données de test. Vous répétez ce processus jusqu'à obtenir un pipeline qui vous satisfait.

Ensuite, vous ajoutez peut-être un simple paramètre vrai ou faux comme `--skip_process`. Maintenant vous devez exécuter le pipeline deux fois, une fois avec chaque paramètre pour vous assurer qu'il fonctionne comme prévu. Mais attendez, comment vérifier si `--skip_process` saute réellement le processus ? Nous devons fouiller dans les sorties ou vérifier les fichiers de log ! C'est pénible et sujet aux erreurs.

Au fur et à mesure que vous développez votre pipeline, il deviendra rapidement si complexe que tester manuellement chaque itération sera lent et sujet aux erreurs. De plus, si vous trouvez une erreur, il sera très difficile de déterminer exactement d'où provient l'erreur dans votre pipeline. C'est là que les tests interviennent.

Les tests permettent de vérifier systématiquement que chaque partie de votre pipeline fonctionne comme prévu. Les avantages pour un développeur de tests bien écrits sont énormes :

- **Confiance** : Parce que les tests couvrent l'ensemble du pipeline, vous pouvez être sûr que modifier quelque chose n'affecte rien d'autre
- **Fiabilité** : Lorsque plusieurs développeurs travaillent sur le pipeline, ils savent que les autres développeurs n'ont pas cassé le pipeline et chaque composant.
- **Transparence** : Les tests montrent où un pipeline échoue et facilitent le suivi du problème. Ils fonctionnent également comme une forme de documentation, montrant comment exécuter un processus ou un workflow.
- **Rapidité** : Parce que les tests sont automatisés, ils peuvent être exécutés très rapidement et de manière répétée. Vous pouvez itérer rapidement avec moins de crainte d'introduire de nouveaux bugs.

Il existe de nombreux types de tests différents que nous pouvons écrire :

1. **Tests au niveau module** : Pour les processus individuels
2. **Tests au niveau workflow** : Pour un seul workflow
3. **Tests au niveau pipeline** : Pour le pipeline dans son ensemble
4. **Tests de performance** : Pour la vitesse et l'efficacité du pipeline
5. **Tests de stress** : Évaluation de la performance du pipeline dans des conditions extrêmes pour déterminer ses limites

Tester des processus individuels est analogue aux tests unitaires dans d'autres langages. Tester le workflow ou l'ensemble du pipeline est analogue à ce qu'on appelle les tests d'intégration dans d'autres langages, où nous testons les interactions des composants.

[**nf-test**](https://www.nf-test.com/) est un outil qui vous permet d'écrire des tests au niveau module, workflow et pipeline. En bref, il vous permet de vérifier systématiquement que chaque partie individuelle du pipeline fonctionne comme prévu, _de manière isolée_.

### Objectifs d'apprentissage

Dans cette quête secondaire, vous apprendrez à utiliser nf-test pour écrire un test au niveau workflow pour le pipeline ainsi que des tests au niveau module pour les trois processus qu'il appelle.

À la fin de cette quête secondaire, vous serez capable d'utiliser efficacement les techniques suivantes :

- Initialiser nf-test dans votre projet
- Générer des tests au niveau module et workflow
- Ajouter des types courants d'assertions
- Comprendre quand utiliser les snapshots vs. les assertions de contenu
- Exécuter des tests pour un projet entier

Ces compétences vous aideront à mettre en œuvre une stratégie de test complète dans vos projets de pipeline, garantissant qu'ils sont plus robustes et maintenables.

### Prérequis

Avant de vous lancer dans cette quête secondaire, vous devriez :

- Avoir terminé le tutoriel [Hello Nextflow](../hello_nextflow/README.md) ou un cours équivalent pour débutants.
- Être à l'aise avec les concepts et mécanismes de base de Nextflow (processus, canaux, opérateurs, manipulation de fichiers, métadonnées)

---

## 0. Mise en route

#### Ouvrir l'environnement de formation dans Codespaces

Si vous ne l'avez pas encore fait, assurez-vous d'ouvrir l'environnement de formation comme décrit dans la [Configuration de l'environnement](../envsetup/index.md).

[![Open in GitHub Codespaces](https://github.com/codespaces/badge.svg)](https://codespaces.new/nextflow-io/training?quickstart=1&ref=master)

#### Aller dans le répertoire du projet

Déplaçons-nous dans le répertoire où se trouvent les fichiers pour ce tutoriel.

```bash
cd side-quests/nf-test
```

Vous pouvez configurer VSCode pour se concentrer sur ce répertoire :

```bash
code .
```

#### Examiner le matériel

Vous trouverez un fichier de workflow principal et un fichier CSV appelé `greetings.csv` qui contient l'entrée du pipeline.

```console title="Contenu du répertoire"
.
├── greetings.csv
└── main.nf
```

Pour une description détaillée des fichiers, consultez la [mise en route de Hello Nextflow](../hello_nextflow/00_orientation.md).

Le workflow que nous allons tester est un sous-ensemble du workflow Hello construit dans [Hello Workflow](../hello_nextflow/03_hello_workflow.md).

??? example "Que fait le workflow Hello Nextflow ?"

    Si vous n'avez pas suivi la formation [Hello Nextflow](../hello_nextflow/index.md), voici un aperçu rapide de ce que fait ce workflow simple.

    Le workflow prend un fichier CSV contenant des salutations, exécute quatre étapes de transformation consécutives sur elles, et produit un seul fichier texte contenant une image ASCII d'un personnage amusant disant les salutations.

    Les quatre étapes sont implémentées comme des processus Nextflow (`sayHello`, `convertToUpper`, `collectGreetings`, et `cowpy`) stockés dans des fichiers de module séparés.

    1. **`sayHello`:** Écrit chaque salutation dans son propre fichier de sortie (par exemple, "Hello-output.txt")
    2. **`convertToUpper`:** Convertit chaque salutation en majuscules (par exemple, "HELLO")
    3. **`collectGreetings`:** Collecte toutes les salutations en majuscules dans un seul fichier batch
    4. **`cowpy`:** Génère de l'art ASCII en utilisant l'outil `cowpy`

    Les résultats sont publiés dans un répertoire appelé `results/`, et la sortie finale du pipeline (lorsqu'il est exécuté avec les paramètres par défaut) est un fichier texte contenant de l'art ASCII d'un personnage disant les salutations en majuscules.

    Dans cette quête secondaire, nous utilisons une forme intermédiaire du workflow Hello qui ne contient que les deux premiers processus.

Le sous-ensemble avec lequel nous allons travailler est composé de deux processus : `sayHello` et `convertToUpper`.
Vous pouvez voir le code complet du workflow ci-dessous.

??? example "Code du workflow"

    ```groovy title="main.nf"
    /*
    * Pipeline parameters
    */
    params.input_file = "greetings.csv"

    /*
    * Use echo to print 'Hello World!' to standard out
    */
    process sayHello {

        publishDir 'results', mode: 'copy'

        input:
            val greeting

        output:
            path "${greeting}-output.txt"

        script:
        """
        echo '$greeting' > '$greeting-output.txt'
        """
    }

    /*
    * Use a text replace utility to convert the greeting to uppercase
    */
    process convertToUpper {

        publishDir 'results', mode: 'copy'

        input:
            path input_file

        output:
            path "UPPER-${input_file}"

        script:
        """
        cat '$input_file' | tr '[a-z]' '[A-Z]' > UPPER-${input_file}
        """
    }

    workflow {

        // create a channel for inputs from a CSV file
        greeting_ch = channel.fromPath(params.input_file).splitCsv().flatten()

        // emit a greeting
        sayHello(greeting_ch)

        // convert the greeting to uppercase
        convertToUpper(sayHello.out)
    }
    ```

#### Exécuter le workflow

Exécutons le workflow pour nous assurer qu'il fonctionne comme prévu.

```bash
nextflow run main.nf
```

```console title="Résultat de l'exécution du workflow"
 N E X T F L O W   ~  version 24.10.2

Launching `main.nf` [soggy_linnaeus] DSL2 - revision: bbf79d5c31

executor >  local (6)
[f7/c3be66] sayHello (3)       | 3 of 3 ✔
[cd/e15303] convertToUpper (3) | 3 of 3 ✔
```

FÉLICITATIONS ! Vous venez d'exécuter un test !

"Attendez, quoi ? Je viens juste d'exécuter le workflow et ça a fonctionné ! Comment est-ce un test ?"

Bonne question !

Décomposons ce qui vient de se passer.

Vous avez exécuté le workflow avec les paramètres par défaut, vous avez confirmé qu'il fonctionnait et vous êtes satisfait des résultats. C'est l'essence même des tests. Si vous avez suivi le cours de formation Hello Nextflow, vous aurez remarqué que nous avons toujours commencé chaque section en exécutant le workflow que nous utilisions comme point de départ, pour confirmer que tout était correctement configuré.

Tester un logiciel fait essentiellement ce processus pour nous.

#### Examiner l'exercice

Votre défi est d'ajouter des tests standardisés à ce workflow en utilisant nf-test, afin de faciliter la vérification que chaque partie continue de fonctionner comme prévu en cas de modifications ultérieures.

#### Liste de vérification de préparation

Pensez-vous être prêt à vous lancer ?

- [ ] Je comprends l'objectif de ce cours et ses prérequis
- [ ] Mon codespace est opérationnel
- [ ] J'ai défini mon répertoire de travail de manière appropriée
- [ ] J'ai exécuté le workflow avec succès
- [ ] Je comprends l'exercice

Si vous pouvez cocher toutes les cases, vous êtes prêt à commencer.

---

## 1. Initialiser `nf-test`

Le package `nf-test` fournit une commande d'initialisation qui configure quelques éléments afin que nous puissions commencer à développer des tests pour notre projet.

```bash
nf-test init
```

Cela devrait produire la sortie suivante :

```bash
🚀 nf-test 0.9.3
https://www.nf-test.com
(c) 2021 - 2024 Lukas Forer and Sebastian Schoenherr

Project configured. Configuration is stored in nf-test.config
```

Cela crée également un répertoire `tests` contenant une ébauche de fichier de configuration.

### 1.1. Générer un nf-test

`nf-test` est fourni avec un ensemble d'outils pour construire des fichiers nf-test, nous épargnant la majeure partie du travail. Ceux-ci se trouvent sous la sous-commande `generate`. Générons un test pour le pipeline :

```bash
nf-test generate pipeline main.nf
```

```console title="Sortie"
> nf-test generate pipeline main.nf

Load source file '/workspaces/training/side-quests/nf-test/main.nf'
Wrote pipeline test file '/workspaces/training/side-quests/nf-test/tests/main.nf.test

SUCCESS: Generated 1 test files.
```

Cela créera un fichier `main.nf.test` dans le répertoire `tests`. C'est notre fichier de test au niveau pipeline. Si vous exécutez `tree tests/`, vous devriez voir quelque chose comme ceci :

```console title="Contenu du répertoire tests"
tests/
├── main.nf.test
└── nextflow.config
```

Le fichier `main.nf.test` est notre fichier de test au niveau pipeline. Ouvrons-le et examinons son contenu.

```groovy title="tests/main.nf.test"
nextflow_pipeline {

    name "Test Workflow main.nf"
    script "main.nf"

    test("Should run without failures") {

        when {
            params {
                // define parameters here. Example:
                // outdir = "tests/results"
            }
        }

        then {
            assert workflow.success
        }

    }

}
```

Prenons un moment pour comprendre la structure du fichier de test.

Le bloc `nextflow_pipeline` est le point d'entrée pour tous les tests au niveau pipeline. Il contient :

- `name` : Le nom du test.
- `script` : Le chemin vers le script du pipeline.

Le bloc `test` est le test proprement dit. Il contient :

- `when` : Les conditions dans lesquelles le test doit être exécuté. Cela inclut les paramètres qui seront utilisés pour exécuter le pipeline.
- `then` : Les assertions qui doivent être faites. Cela inclut les résultats attendus du pipeline.

En français clair, la logique du test se lit comme suit :
"**Quand** ces _paramètres_ sont fournis à ce _pipeline_, **alors** nous nous attendons à voir ces résultats."

Ce n'est pas un test fonctionnel, nous démontrerons comment le transformer en un dans la section suivante.

### Une note sur les noms de tests

Dans l'exemple ci-dessus, nous avons utilisé le nom par défaut "Should run without failures" qui est approprié pour un test de base qui vérifie simplement si le pipeline s'exécute avec succès. Cependant, au fur et à mesure que nous ajoutons des cas de test plus spécifiques, nous devrions utiliser des noms plus descriptifs qui indiquent ce que nous testons réellement. Par exemple :

- "Should convert input to uppercase" - lors du test de fonctionnalités spécifiques
- "Should handle empty input gracefully" - lors du test de cas limites
- "Should respect max memory parameter" - lors du test de contraintes de ressources
- "Should create expected output files" - lors du test de génération de fichiers

Les bons noms de tests devraient :

1. Commencer par "Should" pour clarifier quel est le comportement attendu
2. Décrire la fonctionnalité ou le scénario spécifique testé
3. Être suffisamment clairs pour que si le test échoue, vous sachiez quelle fonctionnalité est cassée

À mesure que nous ajouterons plus d'assertions et de cas de test spécifiques plus tard, nous utiliserons ces noms plus descriptifs pour clarifier ce que chaque test vérifie.

### 1.2. Exécuter le test

Exécutons le test pour voir ce qui se passe.

```bash
nf-test test tests/main.nf.test
```

```console title="Échec du test pipeline nf-test"
> nf-test test tests/main.nf.test

🚀 nf-test 0.9.3
https://www.nf-test.com
(c) 2021 - 2024 Lukas Forer and Sebastian Schoenherr


Test Workflow main.nf

  Test [693ba951] 'Should run without failures' FAILED (4.652s)

  Assertion failed:

  assert workflow.success
         |        |
         workflow false

  Nextflow stdout:

  ERROR ~ No such file or directory: /workspaces/training/side-quests/nf-test/.nf-test/tests/693ba951a20fec36a5a9292ed1cc8a9f/greetings.csv

   -- Check '/workspaces/training/side-quests/nf-test/.nf-test/tests/693ba951a20fec36a5a9292ed1cc8a9f/meta/nextflow.log' file for details
  Nextflow stderr:

FAILURE: Executed 1 tests in 4.679s (1 failed)
```

Le test échoue ! Que s'est-il passé ?

1. nf-test a essayé d'exécuter le pipeline tel quel, en utilisant les paramètres du bloc `when` :

```groovy title="tests/main.nf.test"
when {
    params {
        // define parameters here. Example:
        // outdir = "tests/results"
    }
}
```

2. nf-test a vérifié le statut du pipeline et l'a comparé au bloc `when` :

```groovy title="tests/main.nf.test"
then {
    assert workflow.success
}
```

Notez comment nf-test a signalé l'échec du pipeline et fourni le message d'erreur de Nextflow :

```console title="Erreur"
ERROR ~ No such file or directory: /workspaces/training/side-quests/nf-test/.nf-test/tests/693ba951a20fec36a5a9292ed1cc8a9f/greetings.csv
```

Alors quel était le problème ? Rappelez-vous que le pipeline a un fichier greetings.csv dans le répertoire du projet. Lorsque nf-test exécute le pipeline, il cherchera ce fichier, mais ne peut pas le trouver. Le fichier est bien là, que se passe-t-il ? Eh bien, si nous regardons le chemin, nous pouvons voir que le test se déroule dans le chemin `./nf-test/tests/longHashString/`. Tout comme Nextflow, nf-test crée un nouveau répertoire pour chaque test afin de garder tout isolé. Le fichier de données ne se trouve pas là, nous devons donc corriger le chemin vers le fichier dans le test original.

Revenons au fichier de test et modifions le chemin vers le fichier dans le bloc `when`.

Vous vous demandez peut-être comment nous allons pointer vers la racine du pipeline dans le test. Puisque c'est une situation courante, nf-test dispose d'une gamme de variables globales que nous pouvons utiliser pour nous faciliter la vie. Vous pouvez trouver la liste complète [ici](https://www.nf-test.com/docs/testcases/global_variables/) mais en attendant, nous utiliserons la variable `projectDir`, qui désigne la racine du projet de pipeline.

_Avant :_

```groovy title="tests/main.nf.test" linenums="1" hl_lines="3 4"
when {
    params {
        // define parameters here. Example:
        // outdir = "tests/results"
    }
}
```

_Après :_

```groovy title="tests/main.nf.test" linenums="1" hl_lines="3"
when {
    params {
        input_file = "${projectDir}/greetings.csv"
    }
}
```

Exécutons à nouveau le test pour voir s'il fonctionne.

```bash title="Réussite du test pipeline nf-test"
nf-test test tests/main.nf.test
```

```console title="Le pipeline réussit"
> nf-test test tests/main.nf.test

🚀 nf-test 0.9.3
https://www.nf-test.com
(c) 2021 - 2024 Lukas Forer and Sebastian Schoenherr


Test Workflow main.nf

  Test [1d4aaf12] 'Should run without failures' PASSED (1.619s)


SUCCESS: Executed 1 tests in 1.626s
```

Succès ! Le pipeline s'exécute avec succès et le test réussit. Exécutez-le autant de fois que vous le souhaitez et vous obtiendrez toujours le même résultat !

Par défaut, la sortie de Nextflow est masquée, mais pour vous convaincre que nf-test exécute définitivement le workflow, vous pouvez utiliser le flag `--verbose` :

```bash
nf-test test tests/main.nf.test --verbose
```

```console title="Le pipeline exécute tous les processus"
> nf-test test tests/main.nf.test

🚀 nf-test 0.9.3
https://www.nf-test.com
(c) 2021 - 2024 Lukas Forer and Sebastian Schoenherr


Test Workflow main.nf

  Test [693ba951] 'Should run without failures'
    > Nextflow 24.10.4 is available - Please consider updating your version to it
    > N E X T F L O W  ~  version 24.10.0
    > Launching `/workspaces/training/side-quests/nf-test/main.nf` [zen_ampere] DSL2 - revision: bbf79d5c31
    > [2b/61e453] Submitted process > sayHello (2)
    > [31/4e1606] Submitted process > sayHello (1)
    > [bb/5209ee] Submitted process > sayHello (3)
    > [83/83db6f] Submitted process > convertToUpper (2)
    > [9b/3428b1] Submitted process > convertToUpper (1)
    > [ca/0ba51b] Submitted process > convertToUpper (3)
    PASSED (5.206s)


SUCCESS: Executed 1 tests in 5.239s
```

### 1.3. Ajouter des assertions

Une vérification simple consiste à s'assurer que notre pipeline exécute tous les processus attendus et n'en saute aucun silencieusement. Rappelez-vous que notre pipeline exécute 6 processus, un appelé `sayHello` et un appelé `convertToUpper` pour chacune des 3 salutations.

Ajoutons une assertion à notre test pour vérifier que le pipeline exécute le nombre attendu de processus. Nous mettrons également à jour le nom de notre test pour mieux refléter ce que nous testons.

**Avant :**

```groovy title="tests/main.nf.test" linenums="1" hl_lines="1"
    test("Should run without failures") {

        when {
            params {
                input_file = "${projectDir}/greetings.csv"
            }
        }

        then {
            assert workflow.success
        }

    }
```

**Après :**

```groovy title="tests/main.nf.test" linenums="1" hl_lines="1 11"
    test("Should run successfully with correct number of processes") {

        when {
            params {
                input_file = "${projectDir}/greetings.csv"
            }
        }

        then {
            assert workflow.success
            assert workflow.trace.tasks().size() == 6
        }

    }
```

Le nom du test reflète maintenant mieux ce que nous vérifions réellement - non seulement que le pipeline s'exécute sans échec, mais qu'il exécute le nombre attendu de processus.

Exécutons à nouveau le test pour voir s'il fonctionne.

```bash title="Réussite du test pipeline nf-test"
nf-test test tests/main.nf.test
```

```console title="Le pipeline réussit avec les assertions"
🚀 nf-test 0.9.3
https://www.nf-test.com
(c) 2021 - 2024 Lukas Forer and Sebastian Schoenherr


Test Workflow main.nf

  Test [1d4aaf12] 'Should run successfully with correct number of processes' PASSED (1.567s)


SUCCESS: Executed 1 tests in 1.588s
```

Succès ! Le pipeline s'exécute avec succès et le test réussit. Nous avons maintenant commencé à tester les détails du pipeline, ainsi que le statut global.

### 1.4. Tester la sortie

Ajoutons une assertion à notre test pour vérifier que le fichier de sortie a été créé. Nous l'ajouterons comme un test séparé, avec un nom informatif, pour faciliter l'interprétation des résultats.

**Avant :**

```groovy title="tests/main.nf.test" linenums="1" hl_lines="14"
    test("Should run successfully with correct number of processes") {

        when {
            params {
                input_file = "${projectDir}/greetings.csv"
            }
        }

        then {
            assert workflow.success
            assert workflow.trace.tasks().size() == 6
        }

    }
```

**Après :**

```groovy title="tests/main.nf.test" linenums="1" hl_lines="14-33"
    test("Should run successfully with correct number of processes") {

        when {
            params {
                input_file = "${projectDir}/greetings.csv"
            }
        }

        then {
            assert workflow.success
            assert workflow.trace.tasks().size() == 6
        }

    }

    test("Should produce correct output files") {

        when {
            params {
                input_file = "${projectDir}/greetings.csv"
            }
        }

        then {
            assert file("$launchDir/results/Bonjour-output.txt").exists()
            assert file("$launchDir/results/Hello-output.txt").exists()
            assert file("$launchDir/results/Holà-output.txt").exists()
            assert file("$launchDir/results/UPPER-Bonjour-output.txt").exists()
            assert file("$launchDir/results/UPPER-Hello-output.txt").exists()
            assert file("$launchDir/results/UPPER-Holà-output.txt").exists()
        }

    }
```

Exécutez à nouveau le test pour voir s'il fonctionne.

```bash title="Réussite du test pipeline nf-test"
nf-test test tests/main.nf.test
```

```console title="Le pipeline réussit avec les assertions de fichiers"
> nf-test test tests/main.nf.test

🚀 nf-test 0.9.3
https://www.nf-test.com
(c) 2021 - 2024 Lukas Forer and Sebastian Schoenherr


Test Workflow main.nf

  Test [f0e08a68] 'Should run successfully with correct number of processes' PASSED (8.144s)
  Test [d7e32a32] 'Should produce correct output files' PASSED (6.994s)


SUCCESS: Executed 2 tests in 15.165s
```

Succès ! Les tests réussissent parce que le pipeline s'est terminé avec succès, le nombre correct de processus a été exécuté et les fichiers de sortie ont été créés. Cela devrait également vous montrer à quel point il est utile de fournir ces noms informatifs pour vos tests.

Ce n'est que la surface, nous pouvons continuer à écrire des assertions pour vérifier les détails du pipeline, mais pour l'instant passons au test des composants internes du pipeline.

### À retenir

Vous savez comment écrire un nf-test pour un pipeline.

### Et ensuite ?

Apprenez à tester un processus Nextflow.

---

## 2. Tester un processus Nextflow

Nous n'avons pas besoin d'écrire des tests pour chaque partie du pipeline, mais plus nous avons de tests, plus nous pouvons être complets sur le pipeline et plus nous pouvons être confiants qu'il fonctionne comme prévu. Dans cette section, nous allons tester les deux processus du pipeline en tant qu'unités individuelles.

### 2.1. Tester le processus `sayHello`

Commençons par le processus `sayHello`.

Utilisons à nouveau la commande `nf-test generate` pour générer des tests pour le processus.

```bash
nf-test generate process main.nf
```

```console title="Sortie"
> nf-test generate process main.nf

Load source file '/workspaces/training/side-quests/nf-test/main.nf'
Wrote process test file '/workspaces/training/side-quests/nf-test/tests/main.sayhello.nf.test
Wrote process test file '/workspaces/training/side-quests/nf-test/tests/main.converttoupper.nf.test

SUCCESS: Generated 2 test files.
```

Concentrons-nous pour l'instant sur le processus `sayhello` dans le fichier `main.sayhello.nf.test`.

Ouvrons le fichier et examinons son contenu.

```groovy title="tests/main.sayhello.nf.test"
nextflow_process {

    name "Test Process sayHello"
    script "main.nf"
    process "sayHello"

    test("Should run without failures") {

        when {
            params {
                // define parameters here. Example:
                // outdir = "tests/results"
            }
            process {
                """
                // define inputs of the process here. Example:
                // input[0] = file("test-file.txt")
                """
            }
        }

        then {
            assert process.success
            assert snapshot(process.out).match()
        }

    }

}
```

Comme précédemment, nous commençons par les détails du test, suivis des blocs `when` et `then`. Cependant, nous avons également un bloc `process` supplémentaire qui nous permet de définir les entrées du processus.

Exécutons le test pour voir s'il fonctionne.

```bash title="Réussite du test pipeline nf-test"
nf-test test tests/main.sayhello.nf.test
```

```console title="Le test du processus échoue"
> nf-test test tests/main.sayhello.nf.test

🚀 nf-test 0.9.3
https://www.nf-test.com
(c) 2021 - 2024 Lukas Forer and Sebastian Schoenherr


Test Process sayHello

  Test [1eaad118] 'Should run without failures' FAILED (4.876s)

  Assertion failed:

  assert process.success
         |       |
         |       false
         sayHello

  Nextflow stdout:

  Process `sayHello` declares 1 input but was called with 0 arguments
  Nextflow stderr:

FAILURE: Executed 1 tests in 4.884s (1 failed)
```

Le test échoue parce que le processus `sayHello` déclare 1 entrée mais a été appelé avec 0 arguments. Corrigeons cela en ajoutant une entrée au processus. Rappelez-vous de [Hello Workflow](../hello_nextflow/03_hello_workflow.md) (et de la section de mise en route ci-dessus) que notre processus `sayHello` prend une seule entrée de valeur, que nous devrons fournir. Nous devrions également corriger le nom du test pour mieux refléter ce que nous testons.

**Avant :**

```groovy title="tests/main.sayhello.nf.test" linenums="1" hl_lines="1 10 11"
    test("Should run without failures") {

        when {
            params {
                // define parameters here. Example:
                // outdir = "tests/results"
            }
            process {
                """
                // define inputs of the process here. Example:
                // input[0] = file("test-file.txt")
                """
            }
        }

        then {
            assert process.success
            assert snapshot(process.out).match()
        }

    }
```

**Après :**

```groovy title="tests/main.sayhello.nf.test" linenums="1" hl_lines="1 10"
    test("Should run without failures and produce correct output") {

        when {
            params {
                // define parameters here. Example:
                // outdir = "tests/results"
            }
            process {
                """
                input[0] = "hello"
                """
            }
        }

        then {
            assert process.success
            assert snapshot(process.out).match()
        }

    }
```

Exécutons à nouveau le test pour voir s'il fonctionne.

```console title="Réussite du test pipeline nf-test"
> nf-test test tests/main.sayhello.nf.test

🚀 nf-test 0.9.3
https://www.nf-test.com
(c) 2021 - 2024 Lukas Forer and Sebastian Schoenherr


Test Process sayHello

  Test [f91a1bcd] 'Should run without failures and produce correct output' PASSED (1.604s)
  Snapshots:
    1 created [Should run without failures and produce correct output]


Snapshot Summary:
  1 created

SUCCESS: Executed 1 tests in 1.611s
```

Succès ! Le test réussit parce que le processus `sayHello` s'est exécuté avec succès et la sortie a été créée.

### 2.2. Examiner le snapshot créé par le test

Si nous regardons le fichier `tests/main.sayhello.nf.test`, nous pouvons voir qu'il utilise une méthode `snapshot()` dans le bloc d'assertion :

```groovy title="tests/main.sayhello.nf.test"
assert snapshot(process.out).match()
```

Cela indique à nf-test de créer un snapshot de la sortie du processus `sayHello`. Jetons un coup d'œil au contenu du fichier de snapshot.

```console title="Contenu du fichier snapshot"
code tests/main.sayhello.nf.test.snap
```

Nous ne l'imprimerons pas ici, mais vous devriez voir un fichier JSON contenant des détails sur le processus et les sorties du processus. En particulier, nous pouvons voir une ligne qui ressemble à ceci :

```json title="Contenu du fichier snapshot"
"0": [
    "hello-output.txt:md5,b1946ac92492d2347c6235b4d2611184"
]
```

Cela représente les sorties créées par le processus `sayHello`, que nous testons explicitement. Si nous réexécutons le test, le programme vérifiera que la nouvelle sortie correspond à la sortie qui a été enregistrée à l'origine. C'est une façon rapide et simple de tester que les sorties du processus ne changent pas, c'est pourquoi nf-test le fournit par défaut.

!!!warning

    Cela signifie que nous devons être sûrs que la sortie que nous enregistrons lors de l'exécution d'origine est correcte !

Si, au cours du développement futur, quelque chose dans le code change et provoque une sortie différente, le test échouera et nous devrons déterminer si le changement est attendu ou non.

- S'il s'avère que quelque chose dans le code s'est cassé, nous devrons le corriger, en s'attendant à ce que le code corrigé passe le test.
- Si c'est un changement attendu (par exemple, l'outil a été amélioré et les résultats sont meilleurs) alors nous devrons mettre à jour le snapshot pour accepter la nouvelle sortie comme référence à faire correspondre. nf-test a un paramètre `--update-snapshot` à cet effet.

Nous pouvons réexécuter le test et voir que le test devrait réussir :

```console title="Réussite du test du processus nf-test avec snapshot"
> nf-test test tests/main.sayhello.nf.test

🚀 nf-test 0.9.3
https://www.nf-test.com
(c) 2021 - 2024 Lukas Forer and Sebastian Schoenherr


Test Process sayHello

  Test [f91a1bcd] 'Should run without failures and produce correct output' PASSED (1.675s)


SUCCESS: Executed 1 tests in 1.685s
```

Succès ! Le test réussit parce que le processus `sayHello` s'est exécuté avec succès et la sortie correspondait au snapshot.

### 2.3. Alternative aux snapshots : Assertions de contenu directes

Bien que les snapshots soient excellents pour détecter tout changement dans la sortie, parfois vous voulez vérifier un contenu spécifique sans être aussi strict sur la correspondance complète du fichier. Par exemple :

- Lorsque des parties de la sortie peuvent changer (horodatages, identifiants aléatoires, etc.) mais que certains contenus clés doivent être présents
- Lorsque vous voulez vérifier des motifs ou valeurs spécifiques dans la sortie
- Lorsque vous voulez rendre le test plus explicite sur ce qui constitue un succès

Voici comment nous pourrions modifier notre test pour vérifier un contenu spécifique :

**Avant :**

```groovy title="tests/main.sayhello.nf.test" linenums="1" hl_lines="1 5 6 17"
    test("Should run without failures and produce correct output") {

        when {
            params {
                // define parameters here. Example:
                // outdir = "tests/results"
            }
            process {
                """
                input[0] = "hello"
                """
            }
        }

        then {
            assert process.success
            assert snapshot(process.out).match()
        }

    }
```

**Après :**

```groovy title="tests/main.sayhello.nf.test" linenums="1" hl_lines="1 5 16 17"
    test("Should run without failures and contain expected greeting") {

        when {
            params {
                // define parameters here
            }
            process {
                """
                input[0] = "hello"
                """
            }
        }

        then {
            assert process.success
            assert path(process.out[0][0]).readLines().contains('hello')
            assert !path(process.out[0][0]).readLines().contains('HELLO')
        }

    }
```

Notez que nf-test voit les sorties du processus comme une liste de listes, donc `process.out[0][0]` récupère la première partie du premier élément du canal (ou 'émission') de ce processus.

Cette approche :

- Rend clair exactement ce que nous attendons dans la sortie
- Est plus résiliente aux changements non pertinents dans la sortie
- Fournit de meilleurs messages d'erreur lorsque les tests échouent
- Permet des validations plus complexes (motifs regex, comparaisons numériques, etc.)

Exécutons le test pour voir s'il fonctionne.

```bash title="Réussite du test pipeline nf-test"
nf-test test tests/main.sayhello.nf.test
```

```console title="Le test du processus échoue"
> nf-test test tests/main.sayhello.nf.test

🚀 nf-test 0.9.3
https://www.nf-test.com
(c) 2021 - 2024 Lukas Forer and Sebastian Schoenherr


Test Process sayHello

  Test [58df4e4b] 'Should run without failures and contain expected greeting' PASSED (7.196s)


SUCCESS: Executed 1 tests in 7.208s
```

### 2.4. Tester le processus `convertToUpper`

Ouvrons le fichier `tests/main.converttoupper.nf.test` et examinons son contenu :

```groovy title="tests/main.converttoupper.nf.test"
nextflow_process {

    name "Test Process convertToUpper"
    script "main.nf"
    process "convertToUpper"

    test("Should run without failures") {

        when {
            params {
                // define parameters here. Example:
                // outdir = "tests/results"
            }
            process {
                """
                // define inputs of the process here. Example:
                // input[0] = file("test-file.txt")
                """
            }
        }

        then {
            assert process.success
            assert snapshot(process.out).match()
        }

    }

}
```

C'est un test similaire au processus `sayHello`, mais il teste le processus `convertToUpper`. Nous savons que celui-ci échouera car tout comme avec `sayHello`, le processus `convertToUpper` prend une seule entrée de chemin, mais nous n'en avons pas spécifié.

Nous devons maintenant fournir un seul fichier d'entrée au processus convertToUpper, qui contient du texte que nous voulons convertir en majuscules. Il existe de nombreuses façons de faire cela :

- Nous pourrions créer un fichier dédié pour tester
- Nous pourrions réutiliser le fichier existant data/greetings.csv
- Nous pourrions le créer à la volée dans le test

Pour l'instant, réutilisons le fichier existant data/greetings.csv en utilisant l'exemple que nous avons utilisé avec le test au niveau pipeline. Comme précédemment, nous pouvons nommer le test pour mieux refléter ce que nous testons, mais cette fois laissons-le "snapshoter" le contenu plutôt que de vérifier des chaînes spécifiques (comme nous l'avons fait dans l'autre processus).

**Avant :**

```groovy title="tests/main.converttoupper.nf.test" linenums="1" hl_lines="1 10 11"
    test("Should run without failures") {

        when {
            params {
                // define parameters here. Example:
                // outdir = "tests/results"
            }
            process {
                """
                // define inputs of the process here. Example:
                // input[0] = file("test-file.txt")
                """
            }
        }

        then {
            assert process.success
            assert snapshot(process.out).match()
        }

    }
```

**Après :**

```groovy title="tests/main.converttoupper.nf.test" linenums="1" hl_lines="1 10"
    test("Should run without failures and produce correct output") {

        when {
            params {
                // define parameters here. Example:
                // outdir = "tests/results"
            }
            process {
                """
                input[0] = "${projectDir}/greetings.csv"
                """
            }
        }

        then {
            assert process.success
            assert snapshot(process.out).match()
        }

    }
```

Et exécutons le test !

```bash title="Réussite du test pipeline nf-test"
nf-test test tests/main.converttoupper.nf.test
```

```console title="Réussite du test du processus convertToUpper nf-test"
> nf-test test tests/main.converttoupper.nf.test

🚀 nf-test 0.9.3
https://www.nf-test.com
(c) 2021 - 2024 Lukas Forer and Sebastian Schoenherr


Test Process convertToUpper

  Test [c59b6044] 'Should run without failures and produce correct output' PASSED (1.755s)
  Snapshots:
    1 created [Should run without failures and produce correct output]


Snapshot Summary:
  1 created

SUCCESS: Executed 1 tests in 1.764s
```

Notez que nous avons créé un fichier snapshot pour le processus `convertToUpper` à `tests/main.converttoupper.nf.test.snap`. Si nous exécutons le test à nouveau, nous devrions voir que nf-test réussit à nouveau.

```bash title="Réussite du test du processus convertToUpper nf-test"
nf-test test tests/main.converttoupper.nf.test
```

```console title="Réussite du test du processus convertToUpper nf-test"
> nf-test test tests/main.converttoupper.nf.test

🚀 nf-test 0.9.3
https://www.nf-test.com
(c) 2021 - 2024 Lukas Forer and Sebastian Schoenherr


Test Process convertToUpper

  Test [c59b6044] 'Should run without failures and produce correct output' PASSED (1.798s)


SUCCESS: Executed 1 tests in 1.811s
```

### À retenir

Vous savez comment écrire des tests pour un processus Nextflow et les exécuter.

### Et ensuite ?

Apprenez à exécuter des tests pour tout en une seule fois !

## 3. Exécuter des tests pour l'ensemble du dépôt

Exécuter nf-test sur chaque composant est correct, mais laborieux et sujet aux erreurs. Ne peut-on pas simplement tout tester en une seule fois ?

Oui, nous le pouvons !

Exécutons nf-test sur l'ensemble du dépôt.

### 3.1. Exécuter nf-test sur l'ensemble du dépôt

Nous pouvons exécuter nf-test sur l'ensemble du dépôt en exécutant la commande `nf-test test`.

```bash
nf-test test .
```

Notez que nous utilisons simplement le `.` pour exécuter tout depuis notre répertoire actuel. Cela inclura tous les tests !

```console title="Réussite du test du dépôt nf-test"
> nf-test test .

🚀 nf-test 0.9.3
https://www.nf-test.com
(c) 2021 - 2024 Lukas Forer and Sebastian Schoenherr


Test Process convertToUpper

  Test [3d26d9af] 'Should run without failures and produce correct output' PASSED (4.155s)

Test Workflow main.nf

  Test [f183df37] 'Should run successfully with correct number of processes' PASSED (3.33s)
  Test [d7e32a32] 'Should produce correct output files' PASSED (3.102s)

Test Process sayHello

  Test [58df4e4b] 'Should run without failures and contain expected greeting' PASSED (2.614s)


SUCCESS: Executed 4 tests in 13.481s
```

Regardez ça ! Nous avons exécuté 4 tests, 1 pour chaque processus et 2 pour l'ensemble du pipeline avec une seule commande. Imaginez à quel point c'est puissant sur une grande base de code !

---

## Résumé

Dans cette quête secondaire, vous avez appris à exploiter les fonctionnalités de nf-test pour créer et exécuter des tests pour des processus individuels ainsi que des tests de bout en bout pour l'ensemble du pipeline.
Vous êtes maintenant conscient des deux principales approches de validation de sortie, les snapshots et les assertions de contenu directes, et savez quand utiliser l'une ou l'autre.
Vous savez également comment exécuter des tests un par un ou pour un projet entier.

L'application de ces techniques dans votre propre travail vous permettra de vous assurer que :

- Votre code fonctionne comme prévu
- Les modifications ne cassent pas les fonctionnalités existantes
- D'autres développeurs peuvent contribuer en toute confiance
- Les problèmes peuvent être identifiés et corrigés rapidement
- Le contenu de sortie correspond aux attentes

### Motifs clés

1. Tests au niveau pipeline :
   - Test de succès de base
   - Vérification du nombre de processus
   - Vérifications d'existence de fichiers de sortie
2. Tests au niveau processus
3. Deux approches de validation de sortie :
   - Utilisation de snapshots pour une vérification complète de la sortie
   - Utilisation d'assertions de contenu directes pour des vérifications de contenu spécifique
4. Exécution de tous les tests dans un dépôt avec une seule commande

### Ressources supplémentaires

Consultez la [documentation nf-test](https://www.nf-test.com/) pour plus de fonctionnalités de test avancées et de meilleures pratiques. Vous pourriez vouloir :

- Ajouter des assertions plus complètes à vos tests
- Écrire des tests pour les cas limites et les conditions d'erreur
- Configurer l'intégration continue pour exécuter les tests automatiquement
- En savoir plus sur d'autres types de tests comme les tests de workflow et de module
- Explorer des techniques de validation de contenu plus avancées

**Rappelez-vous :** Les tests sont une documentation vivante de la façon dont votre code devrait se comporter. Plus vous écrivez de tests, et plus vos assertions sont spécifiques, plus vous pouvez avoir confiance en la fiabilité de votre pipeline.

---

## Et ensuite ?

Retournez au [menu des Quêtes Secondaires](./index.md) ou cliquez sur le bouton en bas à droite de la page pour passer au prochain sujet de la liste.
