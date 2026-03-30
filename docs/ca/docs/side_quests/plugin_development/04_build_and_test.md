# Part 4: Proves

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Traducció assistida per IA - [més informació i suggeriments](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

Els plugins són programari independent en el qual els desenvolupadors de pipelines han de poder confiar.
Provar cada funcionalitat de manera independent, fora d'un pipeline, garanteix que el plugin funcioni correctament abans que ningú l'integri en un workflow.
En aquesta secció, escriureu i executareu proves utilitzant el framework de proves Spock.

!!! tip "Comenceu des d'aquí?"

    Si us incorporeu en aquesta part, copieu la solució de la Part 3 per utilitzar-la com a punt de partida:

    ```bash
    cp -r solutions/3-custom-functions/* .
    ```

    Després canvieu al directori del plugin:

    ```bash
    cd nf-greeting
    ```

Assegureu-vos que sou al directori del plugin:

```bash
cd nf-greeting
```

---

## 1. Per què fer proves?

Una compilació correcta significa que el codi compila, però no comprova que funcioni com s'espera.
Les proves unitàries són petits fragments de codi que comproven automàticament si les vostres funcions produeixen la sortida correcta per a una entrada donada.
Per exemple, una prova podria comprovar que `#!groovy reverseGreeting("Hello")` retorna `"olleH"`.

Les proves són valuoses perquè:

- Detecten errors abans que els usuaris ho facin
- Us donen confiança per fer canvis sense trencar res
- Serveixen com a documentació que mostra com s'han d'utilitzar les funcions

---

## 2. Comprendre les proves Spock

La plantilla del plugin utilitza [Spock](https://spockframework.org/), un framework de proves per a Groovy.
Spock ja està configurat al projecte (mitjançant `build.gradle`), de manera que no cal afegir res.

Si heu utilitzat eines de proves anteriorment (com `pytest` a Python o `testthat` a R), Spock compleix el mateix paper: escriviu petites funcions que criden el vostre codi amb entrades conegudes i comproveu les sortides.
La diferència és que Spock utilitza blocs etiquetats (`given:`, `expect:`, `when:`, `then:`) que s'assemblen a un procés o workflow de Nextflow.

Aquí teniu l'estructura bàsica:

```groovy
def 'should reverse a greeting'() {   // (1)!
    given:                             // (2)!
    def ext = new GreetingExtension()

    expect:                            // (3)!
    ext.reverseGreeting('Hello') == 'olleH'
}
```

1. **Nom de la prova entre cometes**: Descriu què comprova la prova. Utilitzeu anglès senzill.
2. **Bloc `given:`**: Configureu el que necessiteu per a la prova (creeu objectes, prepareu dades)
3. **Bloc `expect:`**: Les comprovacions reals. Cada línia ha de ser `true` perquè la prova passi

Aquesta estructura fa que les proves siguin llegibles: "Donat un objecte d'extensió, s'espera que `reverseGreeting('Hello')` sigui igual a `'olleH'`."

---

## 3. Escriure les proves

Escriviu proves per a les dues funcions que heu creat a la Part 3: `reverseGreeting` i `decorateGreeting`.

### 3.1. Crear la classe de proves

```bash
touch src/test/groovy/training/plugin/GreetingExtensionTest.groovy
```

Obriu-lo a l'editor i afegiu l'esquelet buit de la classe de proves:

```groovy title="src/test/groovy/training/plugin/GreetingExtensionTest.groovy" linenums="1"
package training.plugin

import spock.lang.Specification

/**
 * Proves per a les funcions d'extensió de salutació
 */
class GreetingExtensionTest extends Specification {  // (1)!

}
```

1. Totes les classes de proves Spock estenen `Specification`. Aquest és el punt de partida per a qualsevol fitxer de proves Spock.

### 3.2. Provar reverseGreeting

Afegiu un mètode de prova dins del cos de la classe.
El bloc `given:` crea una instància de `GreetingExtension`, i el bloc `expect:` comprova que `reverseGreeting` inverteix correctament dues entrades diferents.
Això prova la funció directament, sense executar un pipeline.

=== "Després"

    ```groovy title="GreetingExtensionTest.groovy" linenums="1" hl_lines="10-17"
    package training.plugin

    import spock.lang.Specification

    /**
     * Proves per a les funcions d'extensió de salutació
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

    1. Creeu una instància de la vostra extensió per provar-la directament, sense executar un pipeline
    2. Cada línia a `expect:` és una asserció; la prova passa només si totes són `true`

=== "Abans"

    ```groovy title="GreetingExtensionTest.groovy" linenums="1"
    package training.plugin

    import spock.lang.Specification

    /**
     * Proves per a les funcions d'extensió de salutació
     */
    class GreetingExtensionTest extends Specification {

    }
    ```

### 3.3. Provar decorateGreeting

Afegiu un segon mètode de prova després del primer.
Aquest verifica que `decorateGreeting` envolta la cadena d'entrada amb `***` a cada costat.

=== "Després"

    ```groovy title="GreetingExtensionTest.groovy" linenums="1" hl_lines="18-25"
    package training.plugin

    import spock.lang.Specification

    /**
     * Proves per a les funcions d'extensió de salutació
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

=== "Abans"

    ```groovy title="GreetingExtensionTest.groovy" linenums="1" hl_lines="18"
    package training.plugin

    import spock.lang.Specification

    /**
     * Proves per a les funcions d'extensió de salutació
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

## 4. Executar les proves

```bash
make test
```

??? example "Sortida de les proves"

    ```console
    BUILD SUCCESSFUL in 5s
    6 actionable tasks: 6 executed
    ```

    **On són els resultats de les proves?** Gradle amaga la sortida detallada quan totes les proves passen.
    "BUILD SUCCESSFUL" significa que tot ha funcionat.
    Si alguna prova falla, veureu missatges d'error detallats.

??? exercise "Afegiu una prova de cas límit"

    Afegiu una prova que comprovi que `reverseGreeting` gestiona una cadena buida.
    Què hauria de retornar `reverseGreeting('')`?
    Afegiu la prova, executeu `make test` i verifiqueu que passa.

    ??? solution "Solució"

        Afegiu aquest mètode de prova a `GreetingExtensionTest.groovy`:

        ```groovy
        def 'should handle empty string'() {
            given:
            def ext = new GreetingExtension()

            expect:
            ext.reverseGreeting('') == ''
        }
        ```

        Una cadena buida invertida continua sent una cadena buida.

---

## 5. Veure l'informe de proves

Gradle genera un informe de proves HTML amb resultats detallats per a cada prova.
Inicieu un servidor web al directori de l'informe:

```bash
pushd build/reports/tests/test
python -m http.server
```

VS Code us demanarà que obriu l'aplicació al navegador.
Navegueu fins a la vostra classe de proves per veure els resultats individuals de cada prova:

![Informe de proves que mostra que totes les proves han passat](./img/test_report.png)

L'informe mostra cada mètode de prova i si ha passat o ha fallat.

Premeu ++ctrl+c++ per aturar el servidor i torneu al directori anterior:

```bash
popd
```

Torneu al directori principal del projecte:

```bash
cd ..
```

---

## Conclusió

Heu après que:

- Les proves Spock utilitzen una estructura llegible `given:`/`expect:`
- Utilitzeu `make test` per executar les proves i `build/reports/tests/test/` per a l'informe HTML
- Les proves verifiquen el comportament i serveixen com a documentació sobre com s'han d'utilitzar les funcions

---

## Què segueix?

Fins ara, el vostre plugin afegeix funcions personalitzades que els pipelines poden cridar.
Els plugins també poden reaccionar a esdeveniments del workflow (una tasca completada, un fitxer publicat, el pipeline finalitzat) mitjançant observadors de traça.
A la següent secció, construireu un observador que compta les tasques completades i imprimeix un resum quan el pipeline finalitza.

[Continueu a la Part 5 :material-arrow-right:](05_observers.md){ .md-button .md-button--primary }
