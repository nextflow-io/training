# Część 1: Podstawy wtyczek

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Tłumaczenie wspomagane przez AI - [dowiedz się więcej i zasugeruj ulepszenia](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

W tej sekcji dowiesz się, jak wtyczki rozszerzają Nextflow'a, a następnie wypróbujesz trzy różne wtyczki w działaniu.

---

## 1. Jak działają wtyczki

Wtyczki rozszerzają Nextflow'a poprzez kilka typów rozszerzeń:

| Typ rozszerzenia    | Co robi                                                    | Przykład                        |
| ------------------- | ---------------------------------------------------------- | ------------------------------- |
| Funkcje             | Dodają własne funkcje wywoływalne z workflow'ów            | `samplesheetToList()`           |
| Monitory workflow'u | Reagują na zdarzenia, takie jak zakończenie zadania        | Własne logowanie, alerty Slack  |
| Executory           | Dodają backendy do wykonywania zadań                       | AWS Batch, Kubernetes           |
| Systemy plików      | Dodają backendy do przechowywania danych                   | S3, Azure Blob                  |

Funkcje i monitory workflow'u (zwane „trace observers" w API Nextflow'a) to najczęstsze typy rozszerzeń dla autorów wtyczek.
Executory i systemy plików są zazwyczaj tworzone przez dostawców platform.

Kolejne ćwiczenia pokazują wtyczki funkcyjne oraz wtyczkę obserwatora, dzięki czemu zobaczysz oba typy w działaniu.

---

## 2. Korzystanie z wtyczek funkcyjnych

Wtyczki funkcyjne dodają wywoływalne funkcje, które importujesz do swoich workflow'ów.
Wypróbujesz dwie: nf-hello (prosty przykład) i nf-schema (szeroko stosowana wtyczka produkcyjna).
Oba ćwiczenia modyfikują ten sam pipeline `hello.nf`, dzięki czemu zobaczysz, jak wtyczki wzbogacają istniejący workflow.

### 2.1. nf-hello: zastąpienie ręcznie napisanego kodu

Wtyczka [nf-hello](https://github.com/nextflow-io/nf-hello) udostępnia funkcję `randomString`, która generuje losowe ciągi znaków.
Pipeline definiuje już własną wbudowaną wersję tej funkcji — zastąpisz ją wersją z wtyczki.

#### 2.1.1. Sprawdź punkt startowy

Przyjrzyj się pipeline'owi:

```bash
cat hello.nf
```

```groovy title="Output"
#!/usr/bin/env nextflow

params.input = 'greetings.csv'

/**
 * Generuje losowy ciąg znaków alfanumerycznych
 */
def randomString(int length) {
    def chars = ('a'..'z') + ('A'..'Z') + ('0'..'9')
    def random = new Random()
    return (1..length).collect { chars[random.nextInt(chars.size())] }.join()
}

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
    greeting_ch = channel.fromPath(params.input)
        .splitCsv(header: true)
        .map { row -> "${row.greeting}_${randomString(8)}" }
    SAY_HELLO(greeting_ch)
    SAY_HELLO.out.view { result -> "Output: ${result.trim()}" }
}
```

Pipeline definiuje własną funkcję `randomString` bezpośrednio w kodzie, a następnie używa jej do dołączenia losowego identyfikatora do każdego powitania.

Uruchom go:

```bash
nextflow run hello.nf
```

```console title="Output"
Output: Hello_aBcDeFgH
Output: Bonjour_xYzWvUtS
Output: Holà_qRsPdMnK
Output: Ciao_jLhGfEcB
Output: Hallo_tNwOiAuR
```

Kolejność wyników i losowe ciągi znaków będą inne, a przy kolejnym uruchomieniu skryptu otrzymasz inny zestaw losowych powitań.

#### 2.1.2. Skonfiguruj wtyczkę

Zastąp wbudowaną funkcję wersją z wtyczki. Dodaj poniższy fragment do pliku `nextflow.config`:

```groovy title="nextflow.config"
// Konfiguracja dla ćwiczeń z tworzenia wtyczek
plugins {
    id 'nf-hello@0.5.0'
}
```

Wtyczki deklaruje się w `nextflow.config` przy użyciu bloku `plugins {}`.
Nextflow automatycznie pobiera je z [Nextflow Plugin Registry](https://registry.nextflow.io/) — centralnego repozytorium wtyczek społecznościowych i oficjalnych.

#### 2.1.3. Użyj funkcji z wtyczki

Zastąp wbudowaną funkcję `randomString` wersją z wtyczki:

=== "Po"

    ```groovy title="hello.nf" hl_lines="3"
    #!/usr/bin/env nextflow

    include { randomString } from 'plugin/nf-hello'

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
        greeting_ch = channel.fromPath(params.input)
            .splitCsv(header: true)
            .map { row -> "${row.greeting}_${randomString(8)}" }
        SAY_HELLO(greeting_ch)
        SAY_HELLO.out.view { result -> "Output: ${result.trim()}" }
    }
    ```

=== "Przed"

    ```groovy title="hello.nf" hl_lines="5-12"
    #!/usr/bin/env nextflow

    params.input = 'greetings.csv'

    /**
     * Generuje losowy ciąg znaków alfanumerycznych
     */
    def randomString(int length) {
        def chars = ('a'..'z') + ('A'..'Z') + ('0'..'9')
        def random = new Random()
        return (1..length).collect { chars[random.nextInt(chars.size())] }.join()
    }

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
        greeting_ch = channel.fromPath(params.input)
            .splitCsv(header: true)
            .map { row -> "${row.greeting}_${randomString(8)}" }
        SAY_HELLO(greeting_ch)
        SAY_HELLO.out.view { result -> "Output: ${result.trim()}" }
    }
    ```

Instrukcja `include` importuje `randomString` z biblioteki, która jest sprawdzona, przetestowana i utrzymywana przez szersze grono współtwórców mogących wykrywać i naprawiać błędy.
Zamiast każdego pipeline'u utrzymującego własną kopię funkcji, każdy pipeline korzystający z wtyczki otrzymuje tę samą zweryfikowaną implementację.
Zmniejsza to ilość zduplikowanego kodu i związane z nim koszty utrzymania.
Składnia `#!groovy include { function } from 'plugin/plugin-id'` to ta sama instrukcja `include`, której używa się dla modułów Nextflow'a, z przedrostkiem `plugin/`.
Możesz przejrzeć [kod źródłowy `randomString`](https://github.com/nextflow-io/nf-hello/blob/e67bddebfa589c7ae51f41bf780c92068dc09e93/plugins/nf-hello/src/main/nextflow/hello/HelloExtension.groovy#L110) w repozytorium nf-hello na GitHub.

#### 2.1.4. Uruchom pipeline

```bash
nextflow run hello.nf
```

```console title="Output"
Pipeline is starting! 🚀
Output: Hello_yqvtclcc
Output: Bonjour_vwwpyzcs
Output: Holà_wrghmgab
Output: Ciao_noniajuy
Output: Hallo_tvrtuxtp
Pipeline complete! 👋
```

(Twoje losowe ciągi znaków będą inne.)

Wyniki nadal mają losowe sufiksy, ale teraz `randomString` pochodzi z wtyczki nf-hello zamiast z kodu wbudowanego.
Komunikaty „Pipeline is starting!" i „Pipeline complete!" są nowe — pochodzą z komponentu obserwatora wtyczki, który poznasz w Części 5.

Nextflow pobiera wtyczki automatycznie przy pierwszym użyciu, więc każdy pipeline deklarujący `nf-hello@0.5.0` otrzymuje dokładnie tę samą przetestowaną funkcję `randomString` bez kopiowania kodu między projektami.

Poznałeś już trzy kroki korzystania z wtyczki funkcyjnej: zadeklaruj ją w `nextflow.config`, zaimportuj funkcję przy użyciu `include` i wywołaj ją w swoim workflow'ie.
Kolejne ćwiczenie stosuje te same kroki do wtyczki produkcyjnej.

### 2.2. nf-schema: walidowane parsowanie CSV

Wtyczka [nf-schema](https://github.com/nextflow-io/nf-schema) jest jedną z najszerzej stosowanych wtyczek Nextflow'a.
Udostępnia funkcję `samplesheetToList`, która parsuje pliki CSV/TSV przy użyciu schematu JSON definiującego oczekiwane kolumny i typy danych.

Pipeline aktualnie odczytuje `greetings.csv` przy użyciu `splitCsv` i ręcznego `map`, ale nf-schema może zastąpić to walidowanym parsowaniem opartym na schemacie.
Plik schematu JSON (`greetings_schema.json`) jest już dostępny w katalogu ćwiczenia.

??? info "Czym jest schemat?"

    Schemat to formalna definicja tego, jak powinny wyglądać poprawne dane.
    Określa takie rzeczy jak oczekiwane kolumny, typ każdej wartości (string, liczba itp.) oraz które pola są wymagane.

    Można go traktować jak kontrakt: jeśli dane wejściowe nie pasują do schematu, narzędzie może wykryć problem wcześnie, zamiast pozwolić mu powodować mylące błędy w dalszej części pipeline'u.

#### 2.2.1. Przyjrzyj się schematowi

```bash
cat greetings_schema.json
```

```json title="Output"
{
  "$schema": "https://json-schema.org/draft/2020-12/schema",
  "type": "array",
  "items": {
    "type": "object",
    "properties": {
      "greeting": {
        "type": "string",
        "description": "The greeting text"
      },
      "language": {
        "type": "string",
        "description": "The language of the greeting"
      }
    },
    "required": ["greeting"]
  }
}
```

Schemat definiuje dwie kolumny (`greeting` i `language`) i oznacza `greeting` jako wymaganą.
Jeśli ktoś poda plik CSV bez kolumny `greeting`, nf-schema wykryje błąd przed uruchomieniem pipeline'u.

#### 2.2.2. Dodaj nf-schema do konfiguracji

Zaktualizuj `nextflow.config`, aby uwzględnić obie wtyczki:

=== "Po"

    ```groovy title="nextflow.config" hl_lines="3"
    plugins {
        id 'nf-hello@0.5.0'
        id 'nf-schema@2.6.1'
    }
    ```

=== "Przed"

    ```groovy title="nextflow.config"
    plugins {
        id 'nf-hello@0.5.0'
    }
    ```

#### 2.2.3. Zaktualizuj hello.nf, aby używał samplesheetToList

Zastąp wejście `splitCsv` funkcją `samplesheetToList`:

=== "Po"

    ```groovy title="hello.nf" hl_lines="4 20 21 22"
    #!/usr/bin/env nextflow

    include { randomString } from 'plugin/nf-hello'
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
        def samplesheet_list = samplesheetToList(params.input, 'greetings_schema.json')
        greeting_ch = Channel.fromList(samplesheet_list)
            .map { row -> "${row[0]}_${randomString(8)}" }
        SAY_HELLO(greeting_ch)
        SAY_HELLO.out.view { result -> "Output: ${result.trim()}" }
    }
    ```

=== "Przed"

    ```groovy title="hello.nf" hl_lines="19 20 21"
    #!/usr/bin/env nextflow

    include { randomString } from 'plugin/nf-hello'

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
        greeting_ch = channel.fromPath(params.input)
            .splitCsv(header: true)
            .map { row -> "${row.greeting}_${randomString(8)}" }
        SAY_HELLO(greeting_ch)
        SAY_HELLO.out.view { result -> "Output: ${result.trim()}" }
    }
    ```

Własny kod parsowania z `splitCsv` i `map` zostaje zastąpiony funkcją `samplesheetToList` — sprawdzoną i przetestowaną, która dodatkowo waliduje samplesheet względem schematu przed uruchomieniem pipeline'u.
Zmniejsza to koszty utrzymania ręcznie napisanej logiki parsowania, jednocześnie poprawiając doświadczenie użytkowników pipeline'u, którzy otrzymują czytelne komunikaty o błędach, gdy ich dane wejściowe nie pasują do oczekiwanego formatu.
Każdy wiersz staje się listą wartości w kolejności kolumn, więc `row[0]` to powitanie, a `row[1]` to język.

#### 2.2.4. Uruchom pipeline

```bash
nextflow run hello.nf
```

```console title="Output"
Pipeline is starting! 🚀
Output: Hello_diozjdwm
Output: Bonjour_speathmm
Output: Holà_dllxnzap
Output: Ciao_wzueddzc
Output: Hallo_hsxwrjbh
Pipeline complete! 👋
```

(Twoje losowe ciągi znaków będą inne.)

Wyniki są takie same, ale teraz schemat waliduje strukturę CSV przed uruchomieniem pipeline'u.
W prawdziwych pipeline'ach ze złożonymi sample sheetami i wieloma kolumnami tego rodzaju walidacja zapobiega błędom, które ręczne `splitCsv` + `map` mogłoby przeoczyć.

#### 2.2.5. Walidacja w działaniu

Aby zobaczyć, co wykrywa walidacja schematu, spróbuj wprowadzić błędy do `greetings.csv`.

Zmień nazwę wymaganej kolumny `greeting` na `message`:

```csv title="greetings.csv" hl_lines="1"
message,language
Hello,English
Bonjour,French
Holà,Spanish
Ciao,Italian
Hallo,German
```

Uruchom pipeline:

```bash
nextflow run hello.nf
```

```console title="Output"
ERROR ~ Validation of samplesheet failed!

The following errors have been detected in greetings.csv:

-> Entry 1: Missing required field(s): greeting
-> Entry 2: Missing required field(s): greeting
-> Entry 3: Missing required field(s): greeting
-> Entry 4: Missing required field(s): greeting
-> Entry 5: Missing required field(s): greeting
```

Pipeline odmawia uruchomienia, ponieważ schemat wymaga kolumny `greeting`, której nie może znaleźć.

Teraz przywróć wymaganą kolumnę, ale zmień nazwę opcjonalnej kolumny `language` na `lang`:

```csv title="greetings.csv" hl_lines="1"
greeting,lang
Hello,English
Bonjour,French
Holà,Spanish
Ciao,Italian
Hallo,German
```

```bash
nextflow run hello.nf
```

Tym razem pipeline uruchamia się, ale wyświetla ostrzeżenie:

```console title="Output (partial)"
WARN: Found the following unidentified headers in greetings.csv:
	- lang
```

Wymagane kolumny powodują błędy krytyczne; opcjonalne — ostrzeżenia.
To właśnie taki wczesny feedback oszczędza czas debugowania w prawdziwych pipeline'ach z dziesiątkami kolumn.

#### 2.2.6. Skonfiguruj zachowanie walidacji

Ostrzeżenie dotyczące `lang` jest przydatne, ale możesz kontrolować jego poziom ważności poprzez konfigurację.
Wtyczki mogą zawierać własne zakresy konfiguracji, które sterują ich zachowaniem.
Wtyczka nf-schema zawiera zakres konfiguracji `validation`; modyfikując ustawienia w tym miejscu, możesz zmienić sposób działania nf-schema.

Dodaj blok `validation` do `nextflow.config`, aby nierozpoznane nagłówki powodowały błąd zamiast ostrzeżenia:

=== "Po"

    ```groovy title="nextflow.config" hl_lines="6-10"
    plugins {
        id 'nf-hello@0.5.0'
        id 'nf-schema@2.6.1'
    }

    validation {
        logging {
            unrecognisedHeaders = "error"
        }
    }
    ```

=== "Przed"

    ```groovy title="nextflow.config"
    plugins {
        id 'nf-hello@0.5.0'
        id 'nf-schema@2.6.1'
    }
    ```

Uruchom pipeline ponownie, pozostawiając kolumnę `lang` bez zmian:

```bash
nextflow run hello.nf
```

```console title="Output (partial)"
Found the following unidentified headers in greetings.csv:
	- lang
 -- Check script 'hello.nf' at line: 20 or see '.nextflow.log' file for more details
```

Pipeline teraz kończy się błędem zamiast ostrzeżeniem.
Kod pipeline'u się nie zmienił — zmieniła się tylko konfiguracja.

Przywróć `greetings.csv` do pierwotnego stanu i usuń blok `validation` przed kontynuowaniem:

```csv title="greetings.csv"
greeting,language
Hello,English
Bonjour,French
Holà,Spanish
Ciao,Italian
Hallo,German
```

```groovy title="nextflow.config"
plugins {
    id 'nf-hello@0.5.0'
    id 'nf-schema@2.6.1'
}
```

Zarówno nf-hello, jak i nf-schema to wtyczki funkcyjne: udostępniają funkcje, które importujesz przy użyciu `include` i wywołujesz w kodzie workflow'u.
Kolejne ćwiczenie pokazuje inny typ wtyczki, który działa bez żadnych instrukcji `include`.

---

## 3. Korzystanie z wtyczki obserwatora: nf-co2footprint

Nie wszystkie wtyczki udostępniają funkcje do importowania.
Wtyczka [nf-co2footprint](https://github.com/nextflow-io/nf-co2footprint) używa **trace observera** do monitorowania zużycia zasobów przez Twój pipeline i szacowania jego śladu węglowego.
Nie musisz zmieniać żadnego kodu pipeline'u — wystarczy dodać ją do konfiguracji.

### 3.1. Dodaj nf-co2footprint do konfiguracji

Zaktualizuj `nextflow.config`:

=== "Po"

    ```groovy title="nextflow.config" hl_lines="4"
    plugins {
        id 'nf-hello@0.5.0'
        id 'nf-schema@2.6.1'
        id 'nf-co2footprint@1.2.0'
    }
    ```

=== "Przed"

    ```groovy title="nextflow.config"
    plugins {
        id 'nf-hello@0.5.0'
        id 'nf-schema@2.6.1'
    }
    ```

### 3.2. Uruchom pipeline

```bash
nextflow run hello.nf
```

Wtyczka generuje kilka komunikatów INFO i WARN podczas wykonywania.
Są one normalne dla małego przykładu uruchamianego na lokalnej maszynie:

```console title="Output (partial)"
nf-co2footprint plugin  ~  version 1.2.0
WARN - [nf-co2footprint] Target zone null not found. Attempting to retrieve carbon intensity for fallback zone GLOBAL.
INFO - [nf-co2footprint] Using fallback carbon intensity from GLOBAL from CI table: 480.0 gCO₂eq/kWh.
WARN - [nf-co2footprint] Executor 'null' not mapped.
WARN - [nf-co2footprint] Fallback to: `machineType = null`, `pue = 1.0`. ...
...
WARN - [nf-co2footprint] No CPU model detected. Using default CPU power draw value (11.41 W).
WARN - [nf-co2footprint] 🔁 Requested memory is null for task 2. Using maximum consumed memory/`peak_rss` (0 GB) for CO₂e footprint computation.
```

Ostrzeżenia dotyczące strefy, executora, modelu procesora i pamięci pojawiają się, ponieważ wtyczka nie może wykryć pełnych szczegółów sprzętowych lokalnego środowiska szkoleniowego.
W środowisku produkcyjnym (np. klastrze HPC lub chmurze) wartości te byłyby dostępne, a szacunki dokładniejsze.

Na końcu poszukaj wiersza podobnego do:

```console title="Output (partial)"
🌱 The workflow run used 126.76 uWh of electricity, resulting in the release of 60.84 ug of CO₂ equivalents into the atmosphere.
```

(Twoje liczby będą inne.)

### 3.3. Przejrzyj raport

Wtyczka generuje pliki wyjściowe w Twoim katalogu roboczym:

```bash
ls co2footprint_*
```

```console title="Output"
co2footprint_report_<timestamp>.html
co2footprint_summary_<timestamp>.txt
co2footprint_trace_<timestamp>.txt
```

Przyjrzyj się podsumowaniu:

```bash
cat co2footprint_summary_*.txt
```

```console title="Output"
Total CO₂e footprint measures of this workflow run (including cached tasks):
  CO₂e emissions: 60.84 ug
  Energy consumption: 126.76 uWh
  CO₂e emissions (market): -

Which equals:
  - 3.48E-7 km travelled by car
  - It takes one tree 0.17s to sequester the equivalent amount of CO₂ from the atmosphere
  - 1.22E-7 % of a flight from Paris to London
```

(Twoje liczby będą inne.)

Pierwsza sekcja pokazuje surowe dane dotyczące zużycia energii i emisji.
Sekcja „Which equals" umieszcza te liczby w perspektywie, przeliczając je na znane odpowiedniki.
Podsumowanie zawiera również sekcję z opcjami konfiguracyjnymi wtyczki oraz cytowanie artykułu naukowego [Green Algorithms](https://doi.org/10.1002/advs.202100707), na którym opiera się metoda obliczeń.

### 3.4. Skonfiguruj wtyczkę

Ostrzeżenie „Target zone null" z sekcji 3.2 pojawiło się, ponieważ wtyczka nie miała skonfigurowanej lokalizacji.
Wtyczka nf-co2footprint definiuje zakres konfiguracji `co2footprint`, w którym możesz ustawić swoją lokalizację geograficzną.

Dodaj blok `co2footprint` do `nextflow.config`:

=== "Po"

    ```groovy title="nextflow.config" hl_lines="7-9"
    plugins {
        id 'nf-hello@0.5.0'
        id 'nf-schema@2.6.1'
        id 'nf-co2footprint@1.2.0'
    }

    co2footprint {
        location = 'GB'
    }
    ```

=== "Przed"

    ```groovy title="nextflow.config"
    plugins {
        id 'nf-hello@0.5.0'
        id 'nf-schema@2.6.1'
        id 'nf-co2footprint@1.2.0'
    }
    ```

!!! tip "Wskazówka"

    Możesz użyć kodu swojego kraju (np. `'US'`, `'DE'`, `'FR'`).

Uruchom pipeline:

```bash
nextflow run hello.nf
```

```console title="Output (partial)"
INFO - [nf-co2footprint] Using fallback carbon intensity from GB from CI table: 163.92 gCO₂eq/kWh.
```

Ostrzeżenie o strefie zniknęło.
Wtyczka używa teraz intensywności emisji dwutlenku węgla właściwej dla Wielkiej Brytanii (163,92 gCO₂eq/kWh) zamiast globalnej wartości zastępczej (480,0 gCO₂eq/kWh).

!!! note "Uwaga"

    Możesz również zobaczyć komunikat `WARN: Unrecognized config option 'co2footprint.location'`.
    Jest on kosmetyczny i można go bezpiecznie zignorować — wtyczka nadal poprawnie odczytuje tę wartość.

W Części 6 stworzysz własny zakres konfiguracji dla swojej wtyczki.

Ta wtyczka działa wyłącznie poprzez mechanizm obserwatora, podłączając się do zdarzeń cyklu życia workflow'u w celu zbierania metryk zasobów i generowania raportu po zakończeniu pipeline'u.

Wypróbowałeś już wtyczki funkcyjne (importowane przy użyciu `include`) oraz wtyczkę obserwatora (aktywowaną wyłącznie przez konfigurację).
To dwa najczęstsze typy rozszerzeń, ale jak pokazuje tabela w sekcji 1, wtyczki mogą również dodawać executory i systemy plików.

---

## 4. Odkrywanie wtyczek

[Nextflow Plugin Registry](https://registry.nextflow.io/) to centralne miejsce do wyszukiwania dostępnych wtyczek.

![Strona wtyczki nf-hello na registry.nextflow.io](img/plugin-registry-nf-hello.png)

Każda strona wtyczki pokazuje jej opis, dostępne wersje, instrukcje instalacji oraz linki do dokumentacji.

---

## 5. Przygotowanie do tworzenia wtyczek

Kolejne sekcje (Części 2–6) używają osobnego pliku pipeline'u, `greet.nf`, który korzysta z nf-schema, ale nie z nf-hello ani nf-co2footprint.

Zaktualizuj `nextflow.config`, aby zachować tylko nf-schema:

```groovy title="nextflow.config"
// Konfiguracja dla ćwiczeń z tworzenia wtyczek
plugins {
    id 'nf-schema@2.6.1'
}
```

Usuń pliki wyjściowe co2footprint:

```bash
rm -f co2footprint_*
```

Plik `hello.nf` zachowuje Twoją pracę z Części 1 jako punkt odniesienia; od tej pory będziesz pracować z `greet.nf`.

---

## Podsumowanie

Użyłeś trzech różnych wtyczek:

- **nf-hello**: Wtyczka funkcyjna udostępniająca `randomString`, importowana przy użyciu `include`
- **nf-schema**: Wtyczka funkcyjna udostępniająca `samplesheetToList` do walidowanego parsowania CSV
- **nf-co2footprint**: Wtyczka obserwatora monitorująca zużycie zasobów automatycznie, bez potrzeby użycia `include`

Kluczowe wzorce:

- Wtyczki deklaruje się w `nextflow.config` przy użyciu `#!groovy plugins { id 'plugin-name@version' }`
- Wtyczki funkcyjne wymagają `#!groovy include { function } from 'plugin/plugin-id'`
- Wtyczki obserwatora działają automatycznie po zadeklarowaniu w konfiguracji
- Wtyczki mogą definiować zakresy konfiguracji (np. `#!groovy validation {}`, `#!groovy co2footprint {}`) w celu dostosowania zachowania
- [Nextflow Plugin Registry](https://registry.nextflow.io/) zawiera listę dostępnych wtyczek

---

## Co dalej?

Kolejne sekcje pokazują, jak zbudować własną wtyczkę.
Jeśli nie interesuje Cię tworzenie wtyczek, możesz zakończyć tutaj lub przejść od razu do [Podsumowania](summary.md).

[Przejdź do Części 2 :material-arrow-right:](02_create_project.md){ .md-button .md-button--primary }
