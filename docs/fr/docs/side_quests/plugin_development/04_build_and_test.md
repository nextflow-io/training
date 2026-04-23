# Partie 4 : Tests

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Traduction assistée par IA - [en savoir plus et suggérer des améliorations](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

Les plugins sont des logiciels autonomes auxquels les développeur·ses de pipelines doivent pouvoir faire confiance.
Tester chaque fonctionnalité indépendamment, en dehors d'un pipeline, garantit que le plugin fonctionne correctement avant que quiconque ne l'intègre dans un workflow.
Dans cette section, vous allez écrire et exécuter des tests à l'aide du framework de test Spock.

!!! tip "Astuce"

    Vous commencez à partir de cette partie ? Copiez la solution de la Partie 3 pour l'utiliser comme point de départ :

    ```bash
    cp -r solutions/3-custom-functions/* .
    ```

    Puis placez-vous dans le répertoire du plugin :

    ```bash
    cd nf-greeting
    ```

Assurez-vous d'être dans le répertoire du plugin :

```bash
cd nf-greeting
```

---

## 1. Pourquoi tester ?

Une compilation réussie signifie que le code se compile, mais ne vérifie pas qu'il fonctionne comme prévu.
Les tests unitaires sont de petits morceaux de code qui vérifient automatiquement si vos fonctions produisent la bonne sortie pour une entrée donnée.
Par exemple, un test peut vérifier que `#!groovy reverseGreeting("Hello")` retourne `"olleH"`.

Les tests sont utiles car :

- Ils détectent les bugs avant les utilisateur·trices
- Ils vous donnent la confiance nécessaire pour effectuer des modifications sans tout casser
- Ils servent de documentation montrant comment les fonctions doivent être utilisées

---

## 2. Comprendre les tests Spock

Le template de plugin utilise [Spock](https://spockframework.org/), un framework de test pour Groovy.
Spock est déjà configuré dans le projet (via `build.gradle`), vous n'avez donc rien à ajouter.

Si vous avez déjà utilisé des outils de test (comme `pytest` en Python ou `testthat` en R), Spock remplit le même rôle : vous écrivez de petites fonctions qui appellent votre code avec des entrées connues et vérifient les sorties.
La différence est que Spock utilise des blocs étiquetés (`given:`, `expect:`, `when:`, `then:`) qui ressemblent à un processus ou un workflow Nextflow.

Voici la structure de base :

```groovy
def 'should reverse a greeting'() {   // (1)!
    given:                             // (2)!
    def ext = new GreetingExtension()

    expect:                            // (3)!
    ext.reverseGreeting('Hello') == 'olleH'
}
```

1. **Nom du test entre guillemets** : Décrit ce que le test vérifie. Utilisez un langage simple.
2. **Bloc `given:`** : Préparez ce dont vous avez besoin pour le test (créer des objets, préparer des données)
3. **Bloc `expect:`** : Les vérifications proprement dites. Chaque ligne doit être `true` pour que le test réussisse

Cette structure rend les tests lisibles : "Étant donné un objet extension, on s'attend à ce que `reverseGreeting('Hello')` soit égal à `'olleH'`."

---

## 3. Écrire les tests

Écrivez des tests pour les deux fonctions que vous avez créées dans la Partie 3 : `reverseGreeting` et `decorateGreeting`.

### 3.1. Créer la classe de test

```bash
touch src/test/groovy/training/plugin/GreetingExtensionTest.groovy
```

Ouvrez-le dans votre éditeur et ajoutez le squelette de classe de test vide :

```groovy title="src/test/groovy/training/plugin/GreetingExtensionTest.groovy" linenums="1"
package training.plugin

import spock.lang.Specification

/**
 * Tests pour les fonctions de l'extension de salutation
 */
class GreetingExtensionTest extends Specification {  // (1)!

}
```

1. Toutes les classes de test Spock étendent `Specification`. C'est le point de départ de tout fichier de test Spock.

### 3.2. Tester reverseGreeting

Ajoutez une méthode de test à l'intérieur du corps de la classe.
Le bloc `given:` crée une instance de `GreetingExtension`, et le bloc `expect:` vérifie que `reverseGreeting` inverse correctement deux entrées différentes.
Cela teste la fonction directement, sans exécuter de pipeline.

=== "Après"

    ```groovy title="GreetingExtensionTest.groovy" linenums="1" hl_lines="10-17"
    package training.plugin

    import spock.lang.Specification

    /**
     * Tests pour les fonctions de l'extension de salutation
     */
    class GreetingExtensionTest extends Specification {

        def 'should reverse a greeting'() {
            given:
            def ext = new GreetingExtension()            // (1)!

            expect:
            ext.reverseGreeting('Hello') == 'olleH'     // (2)!
            ext.reverseGreeting('Bonjour') == 'ruojnoB'
        }
    }
    ```

    1. Créez une instance de votre extension pour la tester directement, sans exécuter de pipeline
    2. Chaque ligne dans `expect:` est une assertion ; le test réussit uniquement si toutes sont `true`

=== "Avant"

    ```groovy title="GreetingExtensionTest.groovy" linenums="1"
    package training.plugin

    import spock.lang.Specification

    /**
     * Tests pour les fonctions de l'extension de salutation
     */
    class GreetingExtensionTest extends Specification {

    }
    ```

### 3.3. Tester decorateGreeting

Ajoutez une deuxième méthode de test après la première.
Celle-ci vérifie que `decorateGreeting` encadre la chaîne d'entrée avec `***` de chaque côté.

=== "Après"

    ```groovy title="GreetingExtensionTest.groovy" linenums="1" hl_lines="18-25"
    package training.plugin

    import spock.lang.Specification

    /**
     * Tests pour les fonctions de l'extension de salutation
     */
    class GreetingExtensionTest extends Specification {

        def 'should reverse a greeting'() {
            given:
            def ext = new GreetingExtension()

            expect:
            ext.reverseGreeting('Hello') == 'olleH'
            ext.reverseGreeting('Bonjour') == 'ruojnoB'
        }

        def 'should decorate a greeting'() {
            given:
            def ext = new GreetingExtension()

            expect:
            ext.decorateGreeting('Hello') == '*** Hello ***'
        }
    }
    ```

=== "Avant"

    ```groovy title="GreetingExtensionTest.groovy" linenums="1" hl_lines="18"
    package training.plugin

    import spock.lang.Specification

    /**
     * Tests pour les fonctions de l'extension de salutation
     */
    class GreetingExtensionTest extends Specification {

        def 'should reverse a greeting'() {
            given:
            def ext = new GreetingExtension()

            expect:
            ext.reverseGreeting('Hello') == 'olleH'
            ext.reverseGreeting('Bonjour') == 'ruojnoB'
        }
    }
    ```

---

## 4. Exécuter les tests

```bash
make test
```

??? example "Sortie des tests"

    ```console
    BUILD SUCCESSFUL in 5s
    6 actionable tasks: 6 executed
    ```

    **Où sont les résultats des tests ?** Gradle masque la sortie détaillée lorsque tous les tests réussissent.
    "BUILD SUCCESSFUL" signifie que tout a fonctionné.
    Si un test échoue, vous verrez des messages d'erreur détaillés.

??? exercise "Ajouter un test de cas limite"

    Ajoutez un test qui vérifie que `reverseGreeting` gère une chaîne vide.
    Que doit retourner `reverseGreeting('')` ?
    Ajoutez le test, exécutez `make test` et vérifiez qu'il réussit.

    ??? solution "Solution"

        Ajoutez cette méthode de test à `GreetingExtensionTest.groovy` :

        ```groovy
        def 'should handle empty string'() {
            given:
            def ext = new GreetingExtension()

            expect:
            ext.reverseGreeting('') == ''
        }
        ```

        Une chaîne vide inversée est toujours une chaîne vide.

---

## 5. Consulter le rapport de tests

Gradle génère un rapport de tests HTML avec les résultats détaillés de chaque test.
Démarrez un serveur web dans le répertoire du rapport :

```bash
pushd build/reports/tests/test
python -m http.server
```

VS Code vous invitera à ouvrir l'application dans votre navigateur.
Naviguez jusqu'à votre classe de test pour voir les résultats individuels de chaque test :

![Rapport de tests indiquant que tous les tests ont réussi](./img/test_report.png)

Le rapport affiche chaque méthode de test et indique si elle a réussi ou échoué.

Appuyez sur ++ctrl+c++ pour arrêter le serveur, puis revenez au répertoire précédent :

```bash
popd
```

Retournez dans le répertoire principal du projet :

```bash
cd ..
```

---

## À retenir

Vous avez appris que :

- Les tests Spock utilisent une structure lisible `given:`/`expect:`
- Utilisez `make test` pour exécuter les tests et `build/reports/tests/test/` pour le rapport HTML
- Les tests vérifient le comportement et servent de documentation sur la façon dont les fonctions doivent être utilisées

---

## Et ensuite ?

Jusqu'à présent, votre plugin ajoute des fonctions personnalisées que les pipelines peuvent appeler.
Les plugins peuvent également réagir aux événements du workflow (une tâche qui se termine, un fichier publié, le pipeline qui se termine) à l'aide d'observateurs de trace.
Dans la prochaine section, vous allez créer un observateur qui compte les tâches terminées et affiche un résumé lorsque le pipeline se termine.

[Continuer vers la Partie 5 :material-arrow-right:](05_observers.md){ .md-button .md-button--primary }
