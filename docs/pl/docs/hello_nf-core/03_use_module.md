# Część 3: Użyj modułu nf-core

W tej trzeciej części szkolenia Hello nf-core pokażemy Ci, jak znaleźć, zainstalować i użyć istniejącego modułu nf-core w Twoim pipeline'ie.

Jedną z wielkich zalet pracy z nf-core jest możliwość wykorzystania gotowych, przetestowanych modułów z repozytorium [nf-core/modules](https://github.com/nf-core/modules).
Zamiast pisać każdy proces od podstaw, możesz zainstalować i użyć modułów utrzymywanych przez społeczność, które przestrzegają najlepszych praktyk.

Aby zademonstrować, jak to działa, zastąpimy niestandardowy moduł `collectGreetings` modułem `cat/cat` z nf-core/modules w pipeline'ie `core-hello`.

??? info "Jak zacząć od tej sekcji"

    Ta sekcja szkolenia zakłada, że ukończyłeś [Część 2: Przepisz Hello dla nf-core](./02_rewrite_hello.md) i masz działający pipeline `core-hello`.

    Jeśli nie ukończyłeś Części 2 lub chcesz zacząć od nowa w tej części, możesz użyć rozwiązania `core-hello-part2` jako punktu wyjścia.
    Uruchom to polecenie z poziomu katalogu `hello-nf-core/`:

    ```bash
    cp -r solutions/core-hello-part2 core-hello
    cd core-hello
    ```

    To da Ci w pełni funkcjonalny pipeline nf-core gotowy do dodawania modułów.
    Możesz sprawdzić, czy działa poprawnie, uruchamiając następujące polecenie:

    ```bash
    nextflow run . --outdir core-hello-results -profile test,docker --validate_params false
    ```

---

## 1. Znajdź i zainstaluj odpowiedni moduł nf-core

Najpierw nauczmy się, jak znaleźć istniejący moduł nf-core i zainstalować go w naszym pipeline'ie.

Naszym celem będzie zastąpienie procesu `collectGreetings`, który używa polecenia Unix `cat` do łączenia wielu plików z powitaniami w jeden.
Łączenie plików to bardzo powszechna operacja, więc jest całkiem prawdopodobne, że może już istnieć moduł w nf-core zaprojektowany do tego celu.

Zanurzmy się w to.

### 1.1. Przeglądaj dostępne moduły na stronie nf-core

Projekt nf-core utrzymuje scentralizowany katalog modułów pod adresem [https://nf-co.re/modules](https://nf-co.re/modules).

Przejdź do strony modułów w przeglądarce internetowej i użyj paska wyszukiwania, aby wyszukać 'concatenate'.

![wyniki wyszukiwania modułów](./img/module-search-results.png)

Jak widzisz, jest całkiem sporo wyników, wiele z nich to moduły zaprojektowane do łączenia bardzo specyficznych typów plików.
Wśród nich powinieneś zobaczyć jeden o nazwie `cat_cat`, który jest ogólnego przeznaczenia.

!!! note "Konwencja nazewnictwa modułów"

    Podkreślenie (`_`) jest używane jako zamiennik znaku ukośnika (`/`) w nazwach modułów.

    Moduły nf-core przestrzegają konwencji nazewnictwa `oprogramowanie/polecenie`, gdy narzędzie udostępnia wiele poleceń, jak `samtools/view` (pakiet samtools, polecenie view) lub `gatk/haplotypecaller` (pakiet GATK, polecenie HaplotypeCaller).
    Dla narzędzi, które udostępniają tylko jedno główne polecenie, moduły używają pojedynczego poziomu, jak `fastqc` lub `multiqc`.

Kliknij na pole modułu `cat_cat`, aby wyświetlić dokumentację modułu.

Strona modułu pokazuje:

- Krótki opis: "A module for concatenation of gzipped or uncompressed files"
- Polecenie instalacji: `nf-core modules install cat/cat`
- Strukturę kanałów wejściowych i wyjściowych
- Dostępne parametry

### 1.2. Wyświetl dostępne moduły z linii poleceń

Alternatywnie możesz również wyszukiwać moduły bezpośrednio z linii poleceń, używając narzędzi nf-core.

```bash
nf-core modules list remote
```

To wyświetli listę wszystkich dostępnych modułów w repozytorium nf-core/modules, choć jest to trochę mniej wygodne, jeśli nie znasz już nazwy modułu, którego szukasz.
Jednak jeśli znasz nazwę, możesz przekierować listę do `grep`, aby znaleźć konkretne moduły:

```bash
nf-core modules list remote | grep 'cat/cat'
```

??? success "Wyjście polecenia"

    ```console
    │ cat/cat
    ```

Pamiętaj tylko, że podejście z `grep` wyciągnie tylko wyniki z wyszukiwanym terminem w nazwie, co nie zadziałałoby dla `cat_cat`.

### 1.3. Uzyskaj szczegółowe informacje o module

Aby zobaczyć szczegółowe informacje o konkretnym module z linii poleceń, użyj polecenia `info`:

```bash
nf-core modules info cat/cat
```

To wyświetla dokumentację modułu, w tym jego wejścia, wyjścia i podstawowe informacje o użyciu.

??? success "Wyjście polecenia"

    ```console

                                              ,--./,-.
              ___     __   __   __   ___     /,-._.--~\
        |\ | |__  __ /  ` /  \ |__) |__         }  {
        | \| |       \__, \__/ |  \ |___     \`-._,-`-,
                                              `._,._,'

        nf-core/tools version 3.4.1 - https://nf-co.re


    ╭─ Module: cat/cat  ─────────────────────────────────────────────────╮
    │ 🌐 Repository: https://github.com/nf-core/modules.git              │
    │ 🔧 Tools: cat                                                      │
    │ 📖 Description: A module for concatenation of gzipped or           │
    │ uncompressed files                                                 │
    ╰────────────────────────────────────────────────────────────────────╯
                      ╷                                          ╷
    📥 Inputs        │Description                               │Pattern
    ╺━━━━━━━━━━━━━━━━━┿━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┿━━━━━━━╸
    input[0]         │                                          │
    ╶─────────────────┼──────────────────────────────────────────┼───────╴
      meta  (map)     │Groovy Map containing sample information  │
                      │e.g. [ id:'test', single_end:false ]      │
    ╶─────────────────┼──────────────────────────────────────────┼───────╴
      files_in  (file)│List of compressed / uncompressed files   │      *
                      ╵                                          ╵
                          ╷                                 ╷
    📥 Outputs           │Description                      │     Pattern
    ╺━━━━━━━━━━━━━━━━━━━━━┿━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┿━━━━━━━━━━━━╸
    file_out             │                                 │
    ╶─────────────────────┼─────────────────────────────────┼────────────╴
      meta  (map)         │Groovy Map containing sample     │
                          │information                      │
    ╶─────────────────────┼─────────────────────────────────┼────────────╴
      ${prefix}  (file)   │Concatenated file. Will be       │ ${file_out}
                          │gzipped if file_out ends with    │
                          │".gz"                            │
    ╶─────────────────────┼─────────────────────────────────┼────────────╴
    versions             │                                 │
    ╶─────────────────────┼─────────────────────────────────┼────────────╴
      versions.yml  (file)│File containing software versions│versions.yml
                          ╵                                 ╵

    💻  Installation command: nf-core modules install cat/cat

    ```

To są dokładnie te same informacje, które możesz znaleźć na stronie internetowej.

### 1.4. Zainstaluj moduł cat/cat

Teraz, gdy znaleźliśmy moduł, którego chcemy, musimy dodać go do kodu źródłowego naszego pipeline'a.

Dobra wiadomość jest taka, że projekt nf-core zawiera narzędzia, które ułatwiają tę część.
W szczególności polecenie `nf-core modules install` umożliwia zautomatyzowanie pobierania kodu i udostępnienia go Twojemu projektowi w jednym kroku.

Przejdź do katalogu Twojego pipeline'a i uruchom polecenie instalacji:

```bash
cd core-hello
nf-core modules install cat/cat
```

Narzędzie może najpierw poprosić Cię o określenie typu repozytorium.
(Jeśli nie, przejdź do "Na koniec narzędzie przystąpi do instalacji modułu.")

??? success "Wyjście polecenia"

    ```console

                                          ,--./,-.
          ___     __   __   __   ___     /,-._.--~\
    |\ | |__  __ /  ` /  \ |__) |__         }  {
    | \| |       \__, \__/ |  \ |___     \`-._,-`-,
                                          `._,._,'

    nf-core/tools version 3.4.1 - https://nf-co.re


    WARNING  'repository_type' not defined in .nf-core.yml
    ? Is this repository a pipeline or a modules repository? (Use arrow keys)
    » Pipeline
      Modules repository
    ```

Jeśli tak, naciśnij enter, aby zaakceptować domyślną odpowiedź (`Pipeline`) i kontynuować.

Narzędzie następnie zaoferuje zmianę konfiguracji Twojego projektu, aby uniknąć tego monitu w przyszłości.

??? success "Wyjście polecenia"

    ```console
        INFO     To avoid this prompt in the future, add the 'repository_type' key to your .nf-core.yml file.
        ? Would you like me to add this config now? [y/n] (y):
    ```

Równie dobrze możemy skorzystać z tego wygodnego narzędzia!
Naciśnij enter, aby zaakceptować domyślną odpowiedź (tak).

Na koniec narzędzie przystąpi do instalacji modułu.

??? success "Wyjście polecenia"

    ```console
    INFO Config added to '.nf-core.yml'
    INFO Reinstalling modules found in 'modules.json' but missing from directory:
    INFO Installing 'cat/cat'
    INFO Use the following statement to include this module:

        include { CAT_CAT } from '../modules/nf-core/cat/cat/main'
    ```

Polecenie automatycznie:

- Pobiera pliki modułu do `modules/nf-core/cat/cat/`
- Aktualizuje `modules.json`, aby śledzić zainstalowany moduł
- Dostarcza Ci poprawną instrukcję `include` do użycia w Twoim workflow'ie

!!! tip

    Zawsze upewnij się, że Twój bieżący katalog roboczy to katalog główny Twojego projektu pipeline'u przed uruchomieniem polecenia instalacji modułu.

Sprawdźmy, czy moduł został poprawnie zainstalowany:

```bash
tree -L 4 modules
```

??? abstract "Zawartość katalogu"

    ```console
    modules
    ├── local
    │   ├── collectGreetings.nf
    │   ├── convertToUpper.nf
    │   ├── cowpy.nf
    │   └── sayHello.nf
    └── nf-core
        └── cat
            └── cat
                ├── environment.yml
                ├── main.nf
                ├── meta.yml
                └── tests

    5 directories, 7 files
    ```

Możesz również zweryfikować instalację, prosząc narzędzie nf-core o wyświetlenie listy lokalnie zainstalowanych modułów:

```bash
nf-core modules list local
```

??? success "Wyjście polecenia"

    ```console
    INFO     Repository type: pipeline
    INFO     Modules installed in '.':

    ┏━━━━━━━━━━━━━┳━━━━━━━━━━━━━━━━━┳━━━━━━━━━━━━━┳━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┳━━━━━━━━━━━━┓
    ┃ Module Name ┃ Repository      ┃ Version SHA ┃ Message                                ┃ Date       ┃
    ┡━━━━━━━━━━━━━╇━━━━━━━━━━━━━━━━━╇━━━━━━━━━━━━━╇━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━╇━━━━━━━━━━━━┩
    │ cat/cat     │ nf-core/modules │ 41dfa3f     │ update meta.yml of all modules (#8747) │ 2025-07-07 │
    └─────────────┴─────────────────┴─────────────┴────────────────────────────────────────┴────────────┘
    ```

To potwierdza, że moduł `cat/cat` jest teraz częścią kodu źródłowego Twojego projektu.

Jednak aby faktycznie użyć nowego modułu, musimy zaimportować go do naszego pipeline'a.

### 1.5. Zaktualizuj importy modułów

Zastąpmy instrukcję `include` dla modułu `collectGreetings` tą dla `CAT_CAT` w sekcji importów pliku workflow'a `workflows/hello.nf`.

Przypominając, narzędzie instalacji modułu dało nam dokładną instrukcję do użycia:

```groovy title="Import statement produced by install command"
include { CAT_CAT } from '../modules/nf-core/cat/cat/main'`
```

Zauważ, że konwencja nf-core polega na używaniu wielkich liter dla nazw modułów podczas ich importowania.

Otwórz [core-hello/workflows/hello.nf](core-hello/workflows/hello.nf) i dokonaj następującej zamiany:

=== "Po"

    ```groovy title="core-hello/workflows/hello.nf" linenums="1" hl_lines="10"
    /*
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    */
    include { paramsSummaryMap       } from 'plugin/nf-schema'
    include { softwareVersionsToYAML } from '../subworkflows/nf-core/utils_nfcore_pipeline'
    include { sayHello               } from '../modules/local/sayHello.nf'
    include { convertToUpper         } from '../modules/local/convertToUpper.nf'
    include { CAT_CAT                } from '../modules/nf-core/cat/cat/main'
    include { cowpy                  } from '../modules/local/cowpy.nf'
    ```

=== "Przed"

    ```groovy title="core-hello/workflows/hello.nf" linenums="1" hl_lines="10"
    /*
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    */
    include { paramsSummaryMap       } from 'plugin/nf-schema'
    include { softwareVersionsToYAML } from '../subworkflows/nf-core/utils_nfcore_pipeline'
    include { sayHello               } from '../modules/local/sayHello.nf'
    include { convertToUpper         } from '../modules/local/convertToUpper.nf'
    include { collectGreetings       } from '../modules/local/collectGreetings.nf'
    include { cowpy                  } from '../modules/local/cowpy.nf'
    ```

Zauważ, jak ścieżka dla modułu nf-core różni się od modułów lokalnych:

- **Moduł nf-core**: `'../modules/nf-core/cat/cat/main'` (odwołuje się do `main.nf`)
- **Moduł lokalny**: `'../modules/local/collectGreetings.nf'` (odwołanie do pojedynczego pliku)

Moduł jest teraz dostępny dla workflow'a, więc wszystko, co musimy zrobić, to zamienić wywołanie `collectGreetings` na użycie `CAT_CAT`. Prawda?

Nie tak szybko.

W tym momencie możesz być kuszony, aby od razu zacząć edytować kod, ale warto poświęcić chwilę na dokładne zbadanie, czego oczekuje nowy moduł i co produkuje.

Zajmiemy się tym jako osobną sekcją, ponieważ obejmuje to nowy mechanizm, którego jeszcze nie omówiliśmy: mapy metadanych.

!!! note

    Opcjonalnie możesz usunąć plik `collectGreetings.nf`:

    ```bash
    rm modules/local/collectGreetings.nf
    ```

    Możesz jednak chcieć go zachować jako punkt odniesienia do zrozumienia różnic między modułami lokalnymi a modułami nf-core.

### Podsumowanie

Wiesz, jak znaleźć moduł nf-core i udostępnić go swojemu projektowi.

### Co dalej?

Oceń, czego wymaga nowy moduł i zidentyfikuj wszelkie ważne zmiany potrzebne do zintegrowania go z pipeline'em.

---

## 2. Oceń wymagania nowego modułu

W szczególności musimy zbadać **interfejs** modułu, tj. jego definicje wejścia i wyjścia, i porównać go z interfejsem modułu, który chcemy zastąpić.
To pozwoli nam określić, czy możemy po prostu potraktować nowy moduł jako bezpośredni zamiennik, czy też będziemy musieli dostosować część okablowania.

Idealnie jest to coś, co powinieneś zrobić _zanim_ w ogóle zainstalujesz moduł, ale hej, lepiej późno niż wcale.
(Tak przy okazji, istnieje polecenie `uninstall`, aby pozbyć się modułów, których już nie chcesz.)

!!! note

    Proces CAT_CAT zawiera dość sprytną obsługę różnych typów kompresji, rozszerzeń plików i tak dalej, które nie są ściśle istotne dla tego, co próbujemy Ci tutaj pokazać, więc zignorujemy większość z tego i skupimy się tylko na częściach, które są ważne.

### 2.1. Porównaj interfejsy obu modułów

Przypominając, tak wygląda interfejs naszego modułu `collectGreetings`:

```groovy title="modules/local/collectGreetings.nf (excerpt)" linenums="1" hl_lines="6-7 10"
process collectGreetings {

    publishDir 'results', mode: 'copy'

    input:
        path input_files
        val batch_name

    output:
        path "COLLECTED-${batch_name}-output.txt" , emit: outfile
```

Moduł `collectGreetings` przyjmuje dwa wejścia:

- `input_files` zawiera jeden lub więcej plików wejściowych do przetworzenia;
- `batch_name` to wartość, której używamy do przypisania nazwy specyficznej dla uruchomienia do pliku wyjściowego, co jest formą metadanych.

Po zakończeniu `collectGreetings` wyemituje pojedynczą ścieżkę pliku z tagiem `outfile`.

Dla porównania interfejs modułu `cat/cat` jest bardziej złożony:

```groovy title="modules/nf-core/cat/cat/main.nf (excerpt)" linenums="1" hl_lines="11 14"
process CAT_CAT {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/pigz:2.3.4' :
        'biocontainers/pigz:2.3.4' }"

    input:
    tuple val(meta), path(files_in)

    output:
    tuple val(meta), path("${prefix}"), emit: file_out
    path "versions.yml"               , emit: versions
```

Moduł CAT_CAT przyjmuje jedno wejście, ale to wejście jest krotką zawierającą dwie rzeczy:

- `meta` to struktura zawierająca metadane, zwana metamapą;
- `files_in` zawiera jeden lub więcej plików wejściowych do przetworzenia, odpowiednik `input_files` z `collectGreetings`.

Po zakończeniu CAT_CAT dostarcza swoje wyjścia w dwóch częściach:

- Kolejna krotka zawierająca metamapę i połączony plik wyjściowy, wyemitowana z tagiem `file_out`;
- Plik `versions.yml`, który przechwytuje informacje o użytej wersji oprogramowania, wyemitowany z tagiem `versions`.

Zauważ również, że domyślnie plik wyjściowy będzie nazwany na podstawie identyfikatora, który jest częścią metadanych (kod nie pokazany tutaj).

Może się to wydawać dużą ilością informacji do śledzenia, patrząc tylko na kod, więc oto diagram, który pomoże Ci zwizualizować, jak wszystko do siebie pasuje.

<figure class="excalidraw">
--8<-- "docs/en/docs/hello_nf-core/img/module_comparison.svg"
</figure>

Widzisz, że oba moduły mają podobne wymagania wejściowe pod względem zawartości (zestaw plików wejściowych plus pewne metadane), ale bardzo różne oczekiwania co do tego, jak ta zawartość jest zapakowana.
Ignorując na razie plik wersji, ich główne wyjście jest również równoważne (połączony plik), z tym że CAT_CAT również emituje metamapę w połączeniu z plikiem wyjściowym.

Różnice w pakowaniu będą dość łatwe do rozwiązania, jak zobaczysz za chwilę.
Jednak aby zrozumieć część z metamapą, musimy przedstawić Ci dodatkowy kontekst.

### 2.2. Zrozumienie metamap

Właśnie powiedzieliśmy Ci, że moduł CAT_CAT oczekuje mapy metadanych jako części swojej krotki wejściowej.
Poświęćmy kilka minut na bliższe przyjrzenie się temu, czym to jest.

**Mapa metadanych**, często nazywana w skrócie **metamapą**, to mapa w stylu Groovy zawierająca informacje o jednostkach danych.
W kontekście pipeline'ów Nextflow jednostki danych mogą być czymkolwiek chcesz: pojedynczymi próbkami, partiami próbek lub całymi zestawami danych.

Zgodnie z konwencją metamapa nf-core jest nazywana `meta` i zawiera wymagane pole `id`, które jest używane do nazywania wyjść i śledzenia jednostek danych.

Na przykład typowa mapa metadanych może wyglądać tak:

```groovy title="Example of sample-level metamap"
[id: 'sample1', single_end: false, strandedness: 'forward']
```

Lub w przypadku, gdy metadane są dołączone na poziomie partii:

```groovy title="Example of batch-level metamap"
[id: 'batch1', date: '25.10.01']
```

Teraz umieśćmy to w kontekście procesu `CAT_CAT`, który oczekuje, że pliki wejściowe będą zapakowane w krotkę z metamapą, i również wyemituje metamapę jako część krotki wyjściowej.

```groovy title="modules/nf-core/cat/cat/main.nf (excerpt)" linenums="1" hl_lines="2 5"
input:
tuple val(meta), path(files_in)

output:
tuple val(meta), path("${prefix}"), emit: file_out
```

W rezultacie każda jednostka danych przemieszcza się przez pipeline z dołączonymi odpowiednimi metadanymi.
Kolejne procesy mogą następnie również łatwo uzyskać dostęp do tych metadanych.

Pamiętasz, jak powiedzieliśmy Ci, że plik wyemitowany przez `CAT_CAT` będzie nazwany na podstawie identyfikatora, który jest częścią metadanych?
To jest odpowiedni kod:

```groovy title="modules/nf-core/cat/cat/main.nf (excerpt)" linenums="35"
prefix   = task.ext.prefix ?: "${meta.id}${getFileSuffix(file_list[0])}"
```

To tłumaczy się mniej więcej tak: jeśli `prefix` jest dostarczony przez system parametrów zewnętrznych zadania (`task.ext`), użyj go do nazwania pliku wyjściowego; w przeciwnym razie utwórz go używając `${meta.id}`, który odpowiada polu `id` w metamapie.

Możesz sobie wyobrazić kanał wejściowy wchodzący do tego modułu z zawartością taką jak ta:

```groovy title="Example input channel contents"
ch_input = [[[id: 'batch1', date: '25.10.01'], ['file1A.txt', 'file1B.txt']],
            [[id: 'batch2', date: '25.10.26'], ['file2A.txt', 'file2B.txt']],
            [[id: 'batch3', date: '25.11.14'], ['file3A.txt', 'file3B.txt']]]
```

Następnie zawartość kanału wyjściowego wychodzącego wygląda tak:

```groovy title="Example output channel contents"
ch_input = [[[id: 'batch1', date: '25.10.01'], 'batch1.txt'],
            [[id: 'batch2', date: '25.10.26'], 'batch2.txt'],
            [[id: 'batch3', date: '25.11.14'], 'batch3.txt']]
```

Jak wspomniano wcześniej, konfiguracja wejściowa `tuple val(meta), path(files_in)` to standardowy wzorzec używany we wszystkich modułach nf-core.

Mamy nadzieję, że zaczynasz widzieć, jak przydatne to może być.
Nie tylko pozwala Ci to nazywać wyjścia na podstawie metadanych, ale możesz również robić takie rzeczy, jak używać ich do stosowania różnych wartości parametrów, a w połączeniu z określonymi operatorami możesz nawet grupować, sortować lub filtrować dane przepływające przez pipeline.

!!! note "Dowiedz się więcej o metadanych"

    Aby uzyskać kompleksowe wprowadzenie do pracy z metadanymi w workflow'ach Nextflow, w tym jak odczytywać metadane z arkuszy próbek i używać ich do dostosowywania przetwarzania, zobacz [Metadane w workflow'ach](../side_quests/metadata).

### 2.3. Podsumuj zmiany do wykonania

Na podstawie tego, co przejrzeliśmy, oto główne zmiany, które musimy wprowadzić w naszym pipeline'ie, aby wykorzystać moduł `cat/cat`:

- Utworzyć metamapę zawierającą nazwę partii;
- Zapakować metamapę w krotkę z zestawem plików wejściowych do połączenia (wychodzących z `convertToUpper`);
- Zamienić wywołanie z `collectGreetings()` na `CAT_CAT`;
- Wyodrębnić plik wyjściowy z krotki wyprodukowanej przez proces `CAT_CAT` przed przekazaniem go do `cowpy`.

To powinno załatwić sprawę! Teraz, gdy mamy plan, jesteśmy gotowi do działania.

### Podsumowanie

Wiesz, jak ocenić interfejs wejściowy i wyjściowy nowego modułu, aby zidentyfikować jego wymagania, i nauczyłeś się, jak metamapy są używane przez pipeline'y nf-core do utrzymywania metadanych ściśle powiązanych z danymi przepływającymi przez pipeline.

### Co dalej?

Zintegruj nowy moduł z workflow'em.

---

## 3. Zintegruj CAT_CAT z workflow'em `hello.nf`

Teraz, gdy wiesz wszystko o metamapach (lub wystarczająco dużo dla celów tego szkolenia), nadszedł czas, aby faktycznie zaimplementować zmiany, które nakreśliliśmy powyżej.

Dla jasności podzielimy to i omówimy każdy krok osobno.

!!! note

    Wszystkie zmiany pokazane poniżej są wprowadzane do logiki workflow'a w bloku `main` w pliku workflow'a `core-hello/workflows/hello.nf`.

### 3.1. Utwórz mapę metadanych

Najpierw musimy utworzyć mapę metadanych dla `CAT_CAT`, pamiętając, że moduły nf-core wymagają, aby metamapa zawierała przynajmniej pole `id`.

Ponieważ nie potrzebujemy żadnych innych metadanych, możemy to uprościć i użyć czegoś takiego:

```groovy title="Syntax example"
def cat_meta = [id: 'test']
```

Z tym że nie chcemy na sztywno kodować wartości `id`; chcemy użyć wartości parametru `params.batch`.
Więc kod staje się:

```groovy title="Syntax example"
def cat_meta = [id: params.batch]
```

Tak, to dosłownie tak proste, aby utworzyć podstawową metamapę.

Dodajmy te linie po wywołaniu `convertToUpper`, usuwając wywołanie `collectGreetings`:

=== "Po"

    ```groovy title="core-hello/workflows/hello.nf" linenums="26" hl_lines="7-8"
        // emit a greeting
        sayHello(ch_samplesheet)

        // convert the greeting to uppercase
        convertToUpper(sayHello.out)

        // create metadata map with batch name as the ID
        def cat_meta = [ id: params.batch ]

        // generate ASCII art of the greetings with cowpy
        cowpy(collectGreetings.out.outfile, params.character)
    ```

=== "Przed"

    ```groovy title="core-hello/workflows/hello.nf" linenums="26" hl_lines="7-8"
        // emit a greeting
        sayHello(ch_samplesheet)

        // convert the greeting to uppercase
        convertToUpper(sayHello.out)

        // collect all the greetings into one file
        collectGreetings(convertToUpper.out.collect(), params.batch)

        // generate ASCII art of the greetings with cowpy
        cowpy(collectGreetings.out.outfile, params.character)
    ```

To tworzy prostą mapę metadanych, gdzie `id` jest ustawione na naszą nazwę partii (która będzie `test` przy użyciu profilu testowego).

### 3.2. Utwórz kanał z krotkami metadanych

Następnie przekształć kanał plików w kanał krotek zawierających metadane i pliki:

=== "Po"

    ```groovy title="core-hello/workflows/hello.nf" linenums="26" hl_lines="10-11"
        // emit a greeting
        sayHello(ch_samplesheet)

        // convert the greeting to uppercase
        convertToUpper(sayHello.out)

        // create metadata map with batch name as the ID
        def cat_meta = [ id: params.batch ]

        // create a channel with metadata and files in tuple format
        ch_for_cat = convertToUpper.out.collect().map { files -> tuple(cat_meta, files) }

        // generate ASCII art of the greetings with cowpy
        cowpy(collectGreetings.out.outfile, params.character)
    ```

=== "Przed"

    ```groovy title="core-hello/workflows/hello.nf" linenums="26"
        // emit a greeting
        sayHello(ch_samplesheet)

        // convert the greeting to uppercase
        convertToUpper(sayHello.out)

        // create metadata map with batch name as the ID
        def cat_meta = [ id: params.batch ]

        // generate ASCII art of the greetings with cowpy
        cowpy(collectGreetings.out.outfile, params.character)
    ```

Dodana linia osiąga dwie rzeczy:

- `.collect()` zbiera wszystkie pliki z wyjścia `convertToUpper` w jedną listę
- `.map { files -> tuple(cat_meta, files) }` tworzy krotkę `[metadane, pliki]` w formacie oczekiwanym przez `CAT_CAT`

To wszystko, co musimy zrobić, aby skonfigurować krotkę wejściową dla `CAT_CAT`.

### 3.3. Wywołaj moduł CAT_CAT

Teraz wywołaj `CAT_CAT` na nowo utworzonym kanale:

=== "Po"

    ```groovy title="core-hello/workflows/hello.nf" linenums="26" hl_lines="13-14"
        // emit a greeting
        sayHello(ch_samplesheet)

        // convert the greeting to uppercase
        convertToUpper(sayHello.out)

        // create metadata map with batch name as the ID
        def cat_meta = [ id: params.batch ]

        // create a channel with metadata and files in tuple format
        ch_for_cat = convertToUpper.out.collect().map { files -> tuple(cat_meta, files) }

        // concatenate files using the nf-core cat/cat module
        CAT_CAT(ch_for_cat)

        // generate ASCII art of the greetings with cowpy
        cowpy(collectGreetings.out.outfile, params.character)
    ```

=== "Przed"

    ```groovy title="core-hello/workflows/hello.nf" linenums="26"
        // emit a greeting
        sayHello(ch_samplesheet)

        // convert the greeting to uppercase
        convertToUpper(sayHello.out)

        // create metadata map with batch name as the ID
        def cat_meta = [ id: params.batch ]

        // create a channel with metadata and files in tuple format
        ch_for_cat = convertToUpper.out.collect().map { files -> tuple(cat_meta, files) }

        // generate ASCII art of the greetings with cowpy
        cowpy(collectGreetings.out.outfile, params.character)
    ```

To kończy najtrudniejszą część tej zamiany, ale jeszcze nie skończyliśmy: nadal musimy zaktualizować sposób, w jaki przekazujemy połączone wyjście do procesu `cowpy`.

### 3.4. Wyodrębnij plik wyjściowy z krotki dla `cowpy`

Wcześniej proces `collectGreetings` po prostu produkował plik, który mogliśmy przekazać bezpośrednio do `cowpy`.
Jednak proces `CAT_CAT` produkuje krotkę, która zawiera metamapę oprócz pliku wyjściowego.

Ponieważ `cowpy` nie akceptuje jeszcze krotek metadanych (naprawimy to w następnej części szkolenia), musimy wyodrębnić plik wyjściowy z krotki wyprodukowanej przez `CAT_CAT` przed przekazaniem go do `cowpy`:

=== "Po"

    ```groovy title="core-hello/workflows/hello.nf" linenums="26" hl_lines="16-17 20"
        // emit a greeting
        sayHello(ch_samplesheet)

        // convert the greeting to uppercase
        convertToUpper(sayHello.out)

        // create metadata map with batch name as the ID
        def cat_meta = [ id: params.batch ]

        // create a channel with metadata and files in tuple format
        ch_for_cat = convertToUpper.out.collect().map { files -> tuple(cat_meta, files) }

        // concatenate the greetings
        CAT_CAT(ch_for_cat)

        // extract the file from the tuple since cowpy doesn't use metadata yet
        ch_for_cowpy = CAT_CAT.out.file_out.map{ meta, file -> file }

        // generate ASCII art of the greetings with cowpy
        cowpy(ch_for_cowpy, params.character)
    ```

=== "Przed"

    ```groovy title="core-hello/workflows/hello.nf" linenums="26" hl_lines="17"
        // emit a greeting
        sayHello(ch_samplesheet)

        // convert the greeting to uppercase
        convertToUpper(sayHello.out)

        // create metadata map with batch name as the ID
        def cat_meta = [ id: params.batch ]

        // create a channel with metadata and files in tuple format
        ch_for_cat = convertToUpper.out.collect().map { files -> tuple(cat_meta, files) }

        // concatenate the greetings
        CAT_CAT(ch_for_cat)

        // generate ASCII art of the greetings with cowpy
        cowpy(collectGreetings.out.outfile, params.character)
    ```

Operacja `.map{ meta, file -> file }` wyodrębnia plik z krotki `[metadane, plik]` wyprodukowanej przez `CAT_CAT` do nowego kanału `ch_for_cowpy`.

Następnie to tylko kwestia przekazania `ch_for_cowpy` do `cowpy` zamiast `collectGreetings.out.outfile` w tej ostatniej linii.

!!! note

    W następnej części szkolenia zaktualizujemy `cowpy`, aby działał bezpośrednio z krotkami metadanych, więc ten krok wyodrębniania nie będzie już konieczny.

### 3.5. Przetestuj workflow

Przetestujmy, czy workflow działa z nowo zintegrowanym modułem `cat/cat`:

```bash
nextflow run . --outdir core-hello-results -profile test,docker --validate_params false
```

To powinno działać dość szybko.

??? success "Wyjście polecenia"

    ```console
    N E X T F L O W ~ version 25.04.3

        Launching `./main.nf` [evil_pike] DSL2 - revision: b9e9b3b8de

        Input/output options
          input                     : /workspaces/training/hello-nf-core/core-hello/assets/greetings.csv
          outdir                    : core-hello-results

        Institutional config options
          config_profile_name       : Test profile
          config_profile_description: Minimal test dataset to check pipeline function

        Generic options
          validate_params           : false
          trace_report_suffix       : 2025-10-30_18-50-58

        Core Nextflow options
          runName                   : evil_pike
          containerEngine           : docker
          launchDir                 : /workspaces/training/hello-nf-core/core-hello
          workDir                   : /workspaces/training/hello-nf-core/core-hello/work
          projectDir                : /workspaces/training/hello-nf-core/core-hello
          userName                  : root
          profile                   : test,docker
          configFiles               : /workspaces/training/hello-nf-core/core-hello/nextflow.config

        !! Only displaying parameters that differ from the pipeline defaults !!
        ------------------------------------------------------
        executor >  local (8)
        [b3/f005fd] CORE_HELLO:HELLO:sayHello (3)       [100%] 3 of 3 ✔
        [08/f923d0] CORE_HELLO:HELLO:convertToUpper (3) [100%] 3 of 3 ✔
        [34/3729a9] CORE_HELLO:HELLO:CAT_CAT (test)     [100%] 1 of 1 ✔
        [24/df918a] CORE_HELLO:HELLO:cowpy              [100%] 1 of 1 ✔
        -[core/hello] Pipeline completed successfully-
    ```

Zauważ, że `CAT_CAT` pojawia się teraz na liście wykonywanych procesów zamiast `collectGreetings`.

I to wszystko! Używamy teraz solidnego, utrzymywanego przez społeczność modułu zamiast niestandardowego kodu na poziomie prototypu dla tego kroku w pipeline'ie.

### Podsumowanie

Teraz wiesz, jak:

- Znaleźć i zainstalować moduły nf-core
- Ocenić wymagania modułu nf-core
- Utworzyć prostą mapę metadanych do użycia z modułem nf-core
- Zintegrować moduł nf-core z Twoim workflow'em

### Co dalej?

Naucz się dostosowywać swoje moduły lokalne, aby przestrzegały konwencji nf-core.
Pokażemy Ci również, jak tworzyć nowe moduły nf-core z szablonu przy użyciu narzędzi nf-core.
