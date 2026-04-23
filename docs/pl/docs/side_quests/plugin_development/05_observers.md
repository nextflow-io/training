# Część 5: Obserwatory śladów

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Tłumaczenie wspomagane przez AI - [dowiedz się więcej i zasugeruj ulepszenia](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

Obserwatory śladów pozwalają wtyczce reagować na zdarzenia workflow'u, takie jak zakończenie zadania, opublikowanie pliku czy zakończenie pipeline'u.
Umożliwia to tworzenie niestandardowych raportów, powiadomień Slack, zbieranie metryk lub integrację z zewnętrznymi systemami monitorowania.
W tej części zbudujesz obserwator, który zlicza ukończone zadania i wyświetla podsumowanie.

!!! tip "Zaczynasz od tej części?"

    Jeśli dołączasz w tym miejscu, skopiuj rozwiązanie z Części 4 jako punkt startowy:

    ```bash
    cp -r solutions/4-build-and-test/* .
    ```

---

## 1. Analiza istniejącego obserwatora śladów

Komunikat „Pipeline is starting!" wyświetlony podczas uruchamiania pipeline'u pochodzi z klasy `GreetingObserver` w Twojej wtyczce.

Przyjrzyj się kodowi obserwatora:

```bash
cat nf-greeting/src/main/groovy/training/plugin/GreetingObserver.groovy
```

```groovy title="Output" hl_lines="30 32-34 37-39"
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
import groovy.util.logging.Slf4j
import nextflow.Session
import nextflow.trace.TraceObserver

/**
 * Implementuje obserwator umożliwiający wykonanie niestandardowej
 * logiki przy zdarzeniach wykonania Nextflow'a.
 */
@Slf4j
@CompileStatic
class GreetingObserver implements TraceObserver {    // (1)!

    @Override
    void onFlowCreate(Session session) {            // (2)!
        println "Pipeline is starting! 🚀"
    }

    @Override
    void onFlowComplete() {                         // (3)!
        println "Pipeline complete! 👋"
    }
}
```

1. Interfejs do podpinania się pod zdarzenia cyklu życia workflow'u
2. Wywoływana przy starcie workflow'u; otrzymuje sesję umożliwiającą dostęp do konfiguracji
3. Wywoływana po pomyślnym zakończeniu workflow'u

Warto zwrócić uwagę na dwie rzeczy:

1. **`class GreetingObserver implements TraceObserver`**: `TraceObserver` to interfejs zdefiniowany przez Nextflow'a. Jeśli Twoja klasa implementuje ten interfejs, Nextflow może się do niej podpiąć i wywoływać Twoje metody w odpowiedzi na zdarzenia.
2. **`@Override`**: Interfejs `TraceObserver` definiuje metody takie jak `onFlowCreate` i `onFlowComplete`. Gdy piszesz metody o tych nazwach i dodajesz adnotację `@Override`, Nextflow wywołuje je we właściwym momencie. Metody, których nie nadpisujesz, są ignorowane.

Pełny zestaw zdarzeń cyklu życia dostępnych w chwili pisania tego materiału:

| Metoda              | Kiedy jest wywoływana            |
| ------------------- | -------------------------------- |
| `onFlowCreate`      | Start workflow'u                 |
| `onFlowComplete`    | Zakończenie workflow'u           |
| `onProcessStart`    | Rozpoczęcie wykonania zadania    |
| `onProcessComplete` | Zakończenie zadania              |
| `onProcessCached`   | Ponowne użycie zadania z cache'u |
| `onFilePublish`     | Opublikowanie pliku              |

Pełna lista dostępna jest w [interfejsie TraceObserver](https://github.com/nextflow-io/nextflow/blob/master/modules/nextflow/src/main/groovy/nextflow/trace/TraceObserver.groovy) w kodzie źródłowym Nextflow'a.

---

## 2. Dodanie obserwatora zliczającego zadania

Celem jest zbudowanie obserwatora, który zlicza ukończone zadania i wyświetla podsumowanie na końcu.
Dodanie nowego obserwatora do wtyczki wymaga dwóch rzeczy: napisania klasy obserwatora i zarejestrowania go w fabryce, aby Nextflow go załadował.

### 2.1. Tworzenie minimalnego obserwatora

Utwórz nowy plik:

```bash
touch nf-greeting/src/main/groovy/training/plugin/TaskCounterObserver.groovy
```

Zacznij od najprostszego możliwego obserwatora, który wyświetla komunikat po zakończeniu każdego zadania:

```groovy title="nf-greeting/src/main/groovy/training/plugin/TaskCounterObserver.groovy" linenums="1"
package training.plugin

import groovy.transform.CompileStatic
import nextflow.processor.TaskHandler       // (1)!
import nextflow.trace.TraceObserver
import nextflow.trace.TraceRecord

/**
 * Obserwator reagujący na zakończenie zadania
 */
@CompileStatic
class TaskCounterObserver implements TraceObserver {  // (2)!

    @Override
    void onProcessComplete(TaskHandler handler, TraceRecord trace) {  // (3)!
        println "✓ Task completed!"
    }
}
```

1. Zaimportuj wymagane klasy: `TraceObserver`, `TaskHandler` i `TraceRecord`
2. Utwórz klasę implementującą `TraceObserver`
3. Nadpisz `onProcessComplete`, aby wykonać kod po zakończeniu zadania

To jest niezbędne minimum:

- Zaimportowanie wymaganych klas (`TraceObserver`, `TaskHandler`, `TraceRecord`)
- Utworzenie klasy implementującej `TraceObserver`
- Nadpisanie `onProcessComplete`, aby coś się działo po zakończeniu zadania

### 2.2. Rejestracja obserwatora

`GreetingFactory` tworzy obserwatory.
Przyjrzyj się jej:

```bash
cat nf-greeting/src/main/groovy/training/plugin/GreetingFactory.groovy
```

```groovy title="Output" hl_lines="25 27-29"
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
import nextflow.trace.TraceObserver
import nextflow.trace.TraceObserverFactory

@CompileStatic
class GreetingFactory implements TraceObserverFactory {

    @Override
    Collection<TraceObserver> create(Session session) {
        return List.<TraceObserver>of(new GreetingObserver())
    }

}
```

Zmodyfikuj `GreetingFactory.groovy`, aby dodać nowy obserwator:

=== "Po"

    ```groovy title="GreetingFactory.groovy" linenums="31" hl_lines="3-6"
    @Override
    Collection<TraceObserver> create(Session session) {
        return [
            new GreetingObserver(),
            new TaskCounterObserver()
        ]
    }
    ```

=== "Przed"

    ```groovy title="GreetingFactory.groovy" linenums="31" hl_lines="3"
    @Override
    Collection<TraceObserver> create(Session session) {
        return List.<TraceObserver>of(new GreetingObserver())
    }
    ```

!!! note "Uwaga"

    Zastąpiliśmy javowy zapis `List.<TraceObserver>of(...)` prostszym literałem listy Groovy `[...]`.
    Oba zwracają `Collection`, ale składnia Groovy jest bardziej czytelna przy dodawaniu wielu elementów.

### 2.3. Budowanie, instalacja i testowanie

```bash
cd nf-greeting && make install && cd ..
nextflow run greet.nf -ansi-log false
```

!!! tip "Wskazówka"

    Domyślnie wyświetlanie postępu ANSI w Nextflow'ie nadpisuje poprzednie linie, pokazując przejrzysty, aktualizowany widok postępu.
    Oznacza to, że widoczna byłaby tylko *ostateczna* liczba zadań, a nie komunikaty pośrednie.

    Użycie `-ansi-log false` wyłącza to zachowanie i pokazuje wszystkie wyjścia sekwencyjnie, co jest niezbędne przy testowaniu obserwatorów wyświetlających komunikaty podczas wykonania.

Powinieneś zobaczyć „✓ Task completed!" wyświetlone pięć razy (raz na zadanie), przeplatane z pozostałymi wyjściami pipeline'u:

```console title="Output (partial)"
...
[9b/df7630] Submitted process > SAY_HELLO (4)
Decorated: *** Hello ***
✓ Task completed!
✓ Task completed!
Decorated: *** Holà ***
✓ Task completed!
...
Pipeline complete! 👋
```

Obserwator działa.
Za każdym razem, gdy zadanie się kończy, Nextflow wywołuje `onProcessComplete`, a nasza implementacja wyświetla komunikat.

??? exercise "Ćwiczenie"

    Spróbuj zmienić komunikat w `onProcessComplete` na własny, przebuduj i uruchom ponownie.
    Potwierdza to, że pełny cykl edycja–budowanie–uruchomienie działa dla obserwatorów.

### 2.4. Dodanie logiki zliczania

Minimalny obserwator potwierdza, że hook działa, ale niczego nie śledzi.

Klasa może przechowywać zmienne (zwane polami lub zmiennymi instancji), które utrzymują się przez cały czas życia obiektu.
Oznacza to, że obserwator może gromadzić stan w trakcie wielu zdarzeń podczas uruchomienia pipeline'u.

Kolejna wersja dodaje zmienną licznika (`taskCount`), która zaczyna od zera.
Za każdym razem, gdy zadanie się kończy, licznik rośnie o jeden.
Po zakończeniu całego workflow'u obserwator wyświetla końcowy wynik.

Zaktualizuj `TaskCounterObserver.groovy` zgodnie z zaznaczonymi zmianami:

```groovy title="nf-greeting/src/main/groovy/training/plugin/TaskCounterObserver.groovy" linenums="1" hl_lines="14 18-19 22-24"
package training.plugin

import groovy.transform.CompileStatic
import nextflow.processor.TaskHandler
import nextflow.trace.TraceObserver
import nextflow.trace.TraceRecord

/**
 * Obserwator zliczający ukończone zadania
 */
@CompileStatic
class TaskCounterObserver implements TraceObserver {

    private int taskCount = 0                // (1)!

    @Override
    void onProcessComplete(TaskHandler handler, TraceRecord trace) {
        taskCount++                          // (2)!
        println "📊 Tasks completed so far: ${taskCount}"
    }

    @Override
    void onFlowComplete() {                  // (3)!
        println "📈 Final task count: ${taskCount}"
    }
}
```

1. `taskCount` to zmienna należąca do obiektu obserwatora. Zachowuje swoją wartość między wywołaniami metod, dzięki czemu może gromadzić licznik przez cały czas trwania workflow'u. `private` oznacza, że tylko ta klasa ma do niej dostęp.
2. `taskCount++` dodaje jeden do licznika. Ta linia wykonuje się za każdym razem, gdy zadanie się kończy, więc licznik rośnie wraz z postępem workflow'u.
3. `onFlowComplete` to drugi hook cyklu życia. Uruchamia się raz po zakończeniu workflow'u, co czyni go dobrym miejscem na wyświetlenie podsumowania.

Podsumowując:

- `taskCount` utrzymuje się między wywołaniami metod, gromadząc licznik przez całe uruchomienie
- `onProcessComplete` zwiększa licznik i wyświetla bieżący wynik za każdym razem, gdy zadanie się kończy
- `onFlowComplete` uruchamia się raz na końcu, wyświetlając końcowy wynik

Przebuduj i przetestuj:

```bash
cd nf-greeting && make install && cd ..
nextflow run greet.nf -ansi-log false
```

??? example "Przykład"

    ```console
    N E X T F L O W  ~  version 25.10.2
    Launching `greet.nf` [pensive_engelbart] DSL2 - revision: 85fefd90d0
    Pipeline is starting! 🚀
    Reversed: olleH
    Reversed: ruojnoB
    Reversed: àloH
    Reversed: oaiC
    Reversed: ollaH
    [be/bd8e72] Submitted process > SAY_HELLO (2)
    [5b/d24c2b] Submitted process > SAY_HELLO (1)
    [14/1f9dbe] Submitted process > SAY_HELLO (3)
    Decorated: *** Bonjour ***
    Decorated: *** Hello ***
    [85/a6b3ad] Submitted process > SAY_HELLO (4)
    📊 Tasks completed so far: 1
    📊 Tasks completed so far: 2
    Decorated: *** Holà ***
    📊 Tasks completed so far: 3
    Decorated: *** Ciao ***
    [3c/be6686] Submitted process > SAY_HELLO (5)
    📊 Tasks completed so far: 4
    Decorated: *** Hallo ***
    📊 Tasks completed so far: 5
    Pipeline complete! 👋
    📈 Final task count: 5
    ```

    Komunikaty licznika przeplatają się z przesyłaniem zadań, ponieważ obserwatory uruchamiają się w miarę kończenia zadań.

---

## 3. Śledzenie opublikowanych plików

Obserwator może również reagować na publikowanie plików.
Metoda `onFilePublish` otrzymuje ścieżki docelową i źródłową, których można użyć do logowania, walidacji lub przetwarzania opublikowanych wyjść.

### 3.1. Dodanie katalogu publikacji

Najpierw zaktualizuj `greet.nf`, aby proces `SAY_HELLO` publikował swoje pliki wyjściowe:

=== "Po"

    ```groovy title="greet.nf" linenums="10" hl_lines="2"
    process SAY_HELLO {
        publishDir 'results'
        input:
            val greeting
        output:
            path 'greeting.txt'
        script:
        // Użyj naszej niestandardowej funkcji wtyczki do ozdobienia powitania
        def decorated = decorateGreeting(greeting)
        """
        echo '$decorated' > greeting.txt
        """
    }
    ```

=== "Przed"

    ```groovy title="greet.nf" linenums="10"
    process SAY_HELLO {
        input:
            val greeting
        output:
            path 'greeting.txt'
        script:
        // Użyj naszej niestandardowej funkcji wtyczki do ozdobienia powitania
        def decorated = decorateGreeting(greeting)
        """
        echo '$decorated' > greeting.txt
        """
    }
    ```

### 3.2. Dodanie metody onFilePublish

Dodaj metodę `onFilePublish` i wymagany import do `TaskCounterObserver.groovy`:

```groovy title="nf-greeting/src/main/groovy/training/plugin/TaskCounterObserver.groovy" linenums="1" hl_lines="5 23-26"
package training.plugin

import groovy.transform.CompileStatic
import nextflow.processor.TaskHandler
import java.nio.file.Path
import nextflow.trace.TraceObserver
import nextflow.trace.TraceRecord

/**
 * Obserwator zliczający ukończone zadania
 */
@CompileStatic
class TaskCounterObserver implements TraceObserver {

    private int taskCount = 0

    @Override
    void onProcessComplete(TaskHandler handler, TraceRecord trace) {
        taskCount++
        println "📊 Tasks completed so far: ${taskCount}"
    }

    @Override
    void onFilePublish(Path destination, Path source) {
        println "📁 Published: ${destination.fileName}"
    }

    @Override
    void onFlowComplete() {
        println "📈 Final task count: ${taskCount}"
    }
}
```

### 3.3. Budowanie i testowanie

```bash
cd nf-greeting && make install && cd ..
nextflow run greet.nf -ansi-log false
```

Powinieneś zobaczyć komunikaty „Published:" dla każdego pliku wyjściowego obok wyjść licznika zadań:

```console title="Output (partial)"
...
📊 Tasks completed so far: 1
📁 Published: greeting.txt
📊 Tasks completed so far: 2
📁 Published: greeting.txt
...
📈 Final task count: 5
Pipeline complete! 👋
```

Metoda `onFilePublish` uruchamia się za każdym razem, gdy Nextflow publikuje plik do katalogu `results`.
Ten wzorzec jest przydatny do budowania dzienników audytu, wyzwalania dalszych działań lub walidacji wyjść w miarę ich powstawania.

---

## Podsumowanie

Nauczyłeś się, że:

- Obserwatory śladów podpinają się pod zdarzenia cyklu życia workflow'u, takie jak `onFlowCreate`, `onProcessComplete`, `onFilePublish` i `onFlowComplete`
- Obserwatory tworzy się przez implementację `TraceObserver` i rejestrację w fabryce
- Obserwatory mogą przechowywać zmienne instancji, aby gromadzić stan w trakcie zdarzeń
- Obserwatory są przydatne do niestandardowego logowania, zbierania metryk, powiadomień i raportowania

---

## Co dalej?

Licznik zadań działa, ale jest zawsze włączony.
W prawdziwej wtyczce użytkownicy powinni móc włączać lub wyłączać funkcje albo dostosowywać zachowanie z poziomu `nextflow.config`, bez edytowania kodu źródłowego wtyczki.
Następna część pokazuje, jak uczynić obserwator konfigurowalnym i jak udostępnić gotową wtyczkę innym.

[Przejdź do Części 6 :material-arrow-right:](06_configuration.md){ .md-button .md-button--primary }
