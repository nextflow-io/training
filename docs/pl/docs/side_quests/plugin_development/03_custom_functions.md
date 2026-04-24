# Część 3: Własne funkcje

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Tłumaczenie wspomagane przez AI - [dowiedz się więcej i zasugeruj ulepszenia](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

Do końca tej sekcji będziesz mieć własne funkcje w swoim pluginie, zbudowane i zainstalowane lokalnie, działające w prawdziwym workflow'u.

!!! tip "Zaczynasz od tej części?"

    Jeśli dołączasz w tym miejscu, skopiuj rozwiązanie z Części 2 jako punkt startowy:

    ```bash
    cp -r solutions/2-create-project/* .
    ```

---

## 1. Sprawdź, co wygenerował szablon

Zanim napiszesz własne funkcje, przyjrzyj się przykładowej funkcji utworzonej przez szablon, żeby zrozumieć wzorzec.

Przejdź do katalogu pluginu:

```bash
cd nf-greeting
```

Szablon utworzył plik o nazwie `GreetingExtension.groovy`, w którym definiuje się funkcje pluginu.
Otwórz go, żeby zobaczyć punkt startowy:

```bash
cat src/main/groovy/training/plugin/GreetingExtension.groovy
```

```groovy title="Output" hl_lines="29 40-43"
/*
 * Copyright 2025, Seqera Labs
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */

package training.plugin

import groovy.transform.CompileStatic
import nextflow.Session
import nextflow.plugin.extension.Function
import nextflow.plugin.extension.PluginExtensionPoint

/**
 * Implementuje własną funkcję, którą można zaimportować
 * w skryptach Nextflow.
 */
@CompileStatic
class GreetingExtension extends PluginExtensionPoint { // (1)!

    @Override
    protected void init(Session session) {             // (2)!
    }

    /**
     * Przywitaj podany cel.
     *
     * @param target
     */
    @Function                                          // (3)!
    void sayHello(String target) {
        println "Hello, ${target}!"
    }

}
```

1. Klasa, na której opiera się Twoje rozszerzenie. Nextflow wymaga jej, żeby rozpoznać Twoje funkcje.
2. Wywoływana przy ładowaniu pluginu; używaj do inicjalizacji
3. Sprawia, że metodę można wywoływać z workflow'ów za pomocą `include`

Szablon zawiera przykładową funkcję `sayHello`.
Adnotacja `@Function` sprawia, że metodę można wywoływać z workflow'ów Nextflow.
Bez niej metoda istnieje wyłącznie wewnątrz kodu pluginu.

W Groovy (i Javie) metody deklarują typ zwracanej wartości oraz typy swoich parametrów.
Na przykład `String reverseGreeting(String greeting)` deklaruje metodę przyjmującą parametr typu `String` i zwracającą `String`.
Słowo kluczowe `void` oznacza, że metoda nic nie zwraca, jak w przypadku `sayHello` powyżej.
Różni się to od Pythona czy R, gdzie typy nie muszą być deklarowane wprost.

---

## 2. Zastąp sayHello przez reverseGreeting

Funkcja `sayHello` z szablonu jest tylko przykładem.
Zastąp ją własną funkcją, żeby przejść pełny cykl: pisanie, budowanie i używanie funkcji pluginu.

Edytuj `src/main/groovy/training/plugin/GreetingExtension.groovy`, zastępując metodę `sayHello`:

=== "Po"

    ```groovy title="GreetingExtension.groovy" linenums="24" hl_lines="8-14"
    @CompileStatic
    class GreetingExtension extends PluginExtensionPoint {

        @Override
        protected void init(Session session) {
        }

        /**
         * Odwróć string z pozdrowieniem
         */
        @Function                                // (1)!
        String reverseGreeting(String greeting) { // (2)!
            return greeting.reverse()             // (3)!
        }

    }
    ```

    1. Sprawia, że metodę można wywoływać z workflow'ów Nextflow
    2. Przyjmuje String, zwraca String
    3. Wbudowana metoda Groovy do odwracania stringów

=== "Przed"

    ```groovy title="GreetingExtension.groovy" linenums="24" hl_lines="12-20"
    /**
     * Implementuje własną funkcję, którą można zaimportować
     * w skryptach Nextflow.
     */
    @CompileStatic
    class GreetingExtension extends PluginExtensionPoint {

        @Override
        protected void init(Session session) {
        }

        /**
         * Przywitaj podany cel.
         *
         * @param target
         */
        @Function
        void sayHello(String target) {
            println "Hello, ${target}!"
        }

    }
    ```

Kluczowe elementy tej funkcji:

- **`@Function`**: Sprawia, że metodę można wywoływać z workflow'ów Nextflow
- **`String reverseGreeting(String greeting)`**: Przyjmuje String, zwraca String
- **`greeting.reverse()`**: Wbudowana metoda Groovy do odwracania stringów

!!! tip "Metody publiczne i prywatne"

    Metody bez adnotacji `@Function` nie są udostępniane workflow'om Nextflow.
    Możesz dodawać metody pomocnicze do swojej klasy bez obawy, że trafią do przestrzeni nazw workflow'u.

---

## 3. Zbuduj i zainstaluj plugin

Zbuduj i zainstaluj plugin:

```bash
make install
```

!!! tip "Jeśli budowanie się nie powiedzie"

    Uważnie przeczytaj komunikat błędu — zazwyczaj zawiera numer linii i opis problemu.
    Najczęstsze przyczyny to błędy składni (brakujący nawias lub cudzysłów), błędnie napisane nazwy klas i niezgodność typów.
    Jeśli utkniesz, porównaj swój kod znak po znaku z przykładami.

---

## 4. Użyj swojej funkcji w workflow'u

Plugin jest zbudowany i zainstalowany.
Następny krok to użycie `reverseGreeting` w workflow'u, żeby sprawdzić, czy działa od początku do końca.

Wróć do katalogu pipeline'u:

```bash
cd ..
```

Edytuj `greet.nf`, żeby zaimportować i użyć `reverseGreeting`:

=== "Po"

    ```groovy title="greet.nf" hl_lines="4 23-25" linenums="1"
    #!/usr/bin/env nextflow

    include { samplesheetToList } from 'plugin/nf-schema'
    include { reverseGreeting } from 'plugin/nf-greeting'

    params.input = 'greetings.csv'

    process SAY_HELLO {
        input:
            val greeting
        output:
            stdout
        script:
        """
        echo '$greeting'
        """
    }

    workflow {
        greeting_ch = Channel.fromList(samplesheetToList(params.input, 'greetings_schema.json'))
            .map { row -> row[0] }

        greeting_ch
            .map { greeting -> reverseGreeting(greeting) }
            .view { reversed -> "Reversed: $reversed" }

        SAY_HELLO(greeting_ch)
        SAY_HELLO.out.view { result -> "Output: ${result.trim()}" }
    }
    ```

=== "Przed"

    ```groovy title="greet.nf" linenums="1"
    #!/usr/bin/env nextflow

    include { samplesheetToList } from 'plugin/nf-schema'

    params.input = 'greetings.csv'

    process SAY_HELLO {
        input:
            val greeting
        output:
            stdout
        script:
        """
        echo '$greeting'
        """
    }

    workflow {
        greeting_ch = Channel.fromList(samplesheetToList(params.input, 'greetings_schema.json'))
            .map { row -> row[0] }
        SAY_HELLO(greeting_ch)
        SAY_HELLO.out.view { result -> "Output: ${result.trim()}" }
    }
    ```

Uruchom pipeline:

```bash
nextflow run greet.nf
```

??? example "Output"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `greet.nf` [elated_marconi] DSL2 - revision: cd8d52c97c

    Pipeline is starting! 🚀
    executor >  local (5)
    [fe/109754] process > SAY_HELLO (5) [100%] 5 of 5 ✔
    Reversed: olleH
    Reversed: ruojnoB
    Reversed: àloH
    Reversed: oaiC
    Reversed: ollaH
    Output: Hello
    Output: Bonjour
    Output: Holà
    Output: Ciao
    Output: Hallo
    Pipeline complete! 👋
    ```

Twoja pierwsza własna funkcja pluginu działa w prawdziwym workflow'u.
Ten sam wzorzec `include { ... } from 'plugin/...'`, którego używałeś z nf-hello i nf-schema w Części 1, działa również z Twoim własnym pluginem.

---

## 5. Dodaj decorateGreeting

Plugin może udostępniać wiele funkcji.
Dodaj drugą, która opakowuje pozdrowienie ozdobnymi znacznikami — w Części 6 uczynisz ją konfigurowalną.

Edytuj `GreetingExtension.groovy`, dodając `decorateGreeting` po `reverseGreeting`, przed zamykającym nawiasem klamrowym klasy:

=== "Po"

    ```groovy title="GreetingExtension.groovy" linenums="24" hl_lines="16-22"
    @CompileStatic
    class GreetingExtension extends PluginExtensionPoint {

        @Override
        protected void init(Session session) {
        }

        /**
         * Odwróć string z pozdrowieniem
         */
        @Function
        String reverseGreeting(String greeting) {
            return greeting.reverse()
        }

        /**
         * Udekoruj pozdrowienie świątecznymi znacznikami
         */
        @Function
        String decorateGreeting(String greeting) {
            return "*** ${greeting} ***"             // (1)!
        }

    }
    ```

    1. Interpolacja stringów w Groovy: `#!groovy ${...}` wstawia wartość zmiennej do stringa

=== "Przed"

    ```groovy title="GreetingExtension.groovy" linenums="24"
    @CompileStatic
    class GreetingExtension extends PluginExtensionPoint {

        @Override
        protected void init(Session session) {
        }

        /**
         * Odwróć string z pozdrowieniem
         */
        @Function
        String reverseGreeting(String greeting) {
            return greeting.reverse()
        }

    }
    ```

Ta funkcja używa interpolacji stringów w Groovy (`"*** ${greeting} ***"`), żeby osadzić zmienną z pozdrowieniem wewnątrz stringa.

Zbuduj, zainstaluj i zaktualizuj workflow:

```bash
cd nf-greeting && make install && cd ..
```

Zaktualizuj `greet.nf`, żeby zaimportować i użyć również `decorateGreeting`:

=== "Po"

    ```groovy title="greet.nf" hl_lines="4-6 14 16-17 19 33" linenums="1"
    #!/usr/bin/env nextflow

    include { samplesheetToList } from 'plugin/nf-schema'
    // Zaimportuj własne funkcje z naszego pluginu
    include { reverseGreeting } from 'plugin/nf-greeting'
    include { decorateGreeting } from 'plugin/nf-greeting'  // (1)!

    params.input = 'greetings.csv'

    process SAY_HELLO {
        input:
            val greeting
        output:
            path 'greeting.txt'
        script:
        // Użyj własnej funkcji pluginu, żeby udekorować pozdrowienie
        def decorated = decorateGreeting(greeting)  // (2)!
        """
        echo '$decorated' > greeting.txt
        """
    }

    workflow {
        greeting_ch = Channel.fromList(samplesheetToList(params.input, 'greetings_schema.json'))
            .map { row -> row[0] }

        // Zademonstruj użycie funkcji reverseGreeting
        greeting_ch
            .map { greeting -> reverseGreeting(greeting) }
            .view { reversed -> "Reversed: $reversed" }

        SAY_HELLO(greeting_ch)
        SAY_HELLO.out.view { file -> "Decorated: ${file.text.trim()}" }
    }
    ```

    1. Wiele funkcji z tego samego pluginu wymaga osobnych instrukcji `include`
    2. Funkcje pluginu działają również wewnątrz bloków `script:` procesu

=== "Przed"

    ```groovy title="greet.nf" linenums="1" hl_lines="4 12 15 28"
    #!/usr/bin/env nextflow

    include { samplesheetToList } from 'plugin/nf-schema'
    include { reverseGreeting } from 'plugin/nf-greeting'

    params.input = 'greetings.csv'

    process SAY_HELLO {
        input:
            val greeting
        output:
            stdout
        script:
        """
        echo '$greeting'
        """
    }

    workflow {
        greeting_ch = Channel.fromList(samplesheetToList(params.input, 'greetings_schema.json'))
            .map { row -> row[0] }

        greeting_ch
            .map { greeting -> reverseGreeting(greeting) }
            .view { reversed -> "Reversed: $reversed" }

        SAY_HELLO(greeting_ch)
        SAY_HELLO.out.view { result -> "Output: ${result.trim()}" }
    }
    ```

```bash
nextflow run greet.nf
```

??? example "Output"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `greet.nf` [elated_marconi] DSL2 - revision: cd8d52c97c

    Pipeline is starting! 🚀
    executor >  local (5)
    [fe/109754] process > SAY_HELLO (5) [100%] 5 of 5 ✔
    Reversed: olleH
    Reversed: ruojnoB
    Reversed: àloH
    Reversed: oaiC
    Reversed: ollaH
    Decorated: *** Hello ***
    Decorated: *** Bonjour ***
    Decorated: *** Holà ***
    Decorated: *** Ciao ***
    Decorated: *** Hallo ***
    Pipeline complete! 👋
    ```

Funkcje pluginu działają zarówno w skryptach procesów (jak `decorateGreeting` wewnątrz `SAY_HELLO`), jak i w operacjach workflow'u (jak `reverseGreeting` w operatorze `map`).

---

## Podsumowanie

Nauczyłeś się, że:

- Funkcje definiuje się za pomocą adnotacji `@Function` w podklasach `PluginExtensionPoint`
- Funkcje pluginu importowane przez `include` działają identycznie niezależnie od tego, czy pochodzą z własnego pluginu, czy z istniejącego
- Funkcje pluginu działają zarówno w skryptach procesów, jak i w operacjach workflow'u

---

## Co dalej?

Twoje funkcje działają, ale jak dotąd sprawdzałeś to wyłącznie uruchamiając cały pipeline i sprawdzając wyniki wzrokiem.
Takie podejście nie skaluje się: w miarę dodawania kolejnych funkcji potrzebujesz szybszego sposobu na weryfikację, że każda z nich działa poprawnie — zwłaszcza po wprowadzeniu zmian.
Następna sekcja przedstawia testy jednostkowe, które pozwalają automatycznie weryfikować poszczególne funkcje bez uruchamiania pipeline'u.

[Przejdź do Części 4 :material-arrow-right:](04_build_and_test.md){ .md-button .md-button--primary }
