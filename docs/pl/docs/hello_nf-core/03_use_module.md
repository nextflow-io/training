# Część 3: Użycie modułu nf-core

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Tłumaczenie wspomagane przez AI - [dowiedz się więcej i zasugeruj ulepszenia](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

W trzeciej części kursu szkoleniowego Hello nf-core pokażemy, jak znaleźć, zainstalować i użyć istniejącego modułu nf-core w Swoim pipeline'ie.

Jedną z głównych korzyści pracy z nf-core jest możliwość wykorzystania wcześniej przygotowanych, przetestowanych modułów z repozytorium [nf-core/modules](https://github.com/nf-core/modules).
Zamiast pisać każdy proces od podstaw, możesz zainstalować i używać gotowych komponentów utrzymywanych przez społeczność, które przestrzegają najlepszych praktyk.

Aby pokazać, jak to działa, zastąpimy niestandardowy moduł `collectGreetings` modułem `cat/cat` z nf-core/modules w pipeline'ie `core-hello`.

??? info "Jak zacząć od tej sekcji"

    Ta sekcja kursu zakłada, że ukończyłeś [Część 2: Przepisanie Hello dla nf-core](./02_rewrite_hello.md) i masz działający pipeline `core-hello`.

    Jeśli nie ukończyłeś Części 2 lub chcesz zacząć od nowa dla tej części, możesz użyć rozwiązania `core-hello-part2` jako punktu wyjścia.
    Uruchom to polecenie z poziomu katalogu `hello-nf-core/`:

    ```bash
    cp -r solutions/core-hello-part2 core-hello
    cd core-hello
    ```

    Otrzymasz w ten sposób w pełni funkcjonalny pipeline nf-core gotowy do dodawania modułów.
    Możesz sprawdzić, czy działa poprawnie, uruchamiając następujące polecenie:

    ```bash
    nextflow run . --outdir core-hello-results -profile test,docker --validate_params false
    ```

---

## 1. Znajdź i zainstaluj odpowiedni moduł nf-core

Najpierw nauczmy się, jak znaleźć istniejący moduł nf-core i zainstalować go w naszym pipeline'ie.

Będziemy dążyć do zastąpienia procesu `collectGreetings`, który używa polecenia Unix `cat` do łączenia wielu plików z powitaniami w jeden.
Ta operacja jest bardzo powszechna, więc prawdopodobnie istnieje już moduł w nf-core zaprojektowany do tego celu.

Zagłębmy się w to.

### 1.1. Przeglądanie dostępnych modułów na stronie nf-core

Projekt nf-core utrzymuje scentralizowany katalog modułów pod adresem [https://nf-co.re/modules](https://nf-co.re/modules).

Przejdź do strony modułów w Swojej przeglądarce internetowej i użyj paska wyszukiwania, aby wyszukać 'concatenate'.

![wyniki wyszukiwania modułów](./img/module-search-results.png)

Jak widać, jest sporo wyników, wiele z nich to moduły zaprojektowane do łączenia bardzo specyficznych typów plików.
Wśród nich powinieneś zobaczyć jeden o nazwie `cat_cat`, który jest ogólnego przeznaczenia.

!!! note "Konwencja nazewnictwa modułów"

    Podkreślenie (`_`) jest używane jako zastępnik znaku ukośnika (`/`) w nazwach modułów.

    Moduły nf-core przestrzegają konwencji nazewnictwa `software/command`, gdy narzędzie dostarcza wiele poleceń, jak `samtools/view` (pakiet samtools, polecenie view) lub `gatk/haplotypecaller` (pakiet GATK, polecenie HaplotypeCaller).
    Dla narzędzi, które dostarczają tylko jedno główne polecenie, moduły używają jednego poziomu, jak `fastqc` lub `multiqc`.

Kliknij na pole modułu `cat_cat`, aby wyświetlić dokumentację modułu.

Strona modułu pokazuje:

- Krótki opis: "A module for concatenation of gzipped or uncompressed files"
- Polecenie instalacji: `nf-core modules install cat/cat`
- Strukturę kanałów wejściowych i wyjściowych
- Dostępne parametry

### 1.2. Wyświetlanie dostępnych modułów z wiersza poleceń

Alternatywnie, możesz również wyszukiwać moduły bezpośrednio z wiersza poleceń używając narzędzi nf-core.

```bash
nf-core modules list remote
```

To wyświetli listę wszystkich dostępnych modułów w repozytorium nf-core/modules, choć jest to nieco mniej wygodne, jeśli nie znasz już nazwy modułu, którego szukasz.
Jednak jeśli znasz, możesz przekierować listę do `grep`, aby znaleźć konkretne moduły:

```bash
nf-core modules list remote | grep 'cat/cat'
```

??? success "Wynik polecenia"

    ```console
    │ cat/cat
    ```

Pamiętaj tylko, że podejście z `grep` wyciągnie tylko wyniki z wyszukiwanym terminem w nazwie, co nie zadziałałoby dla `cat_cat`.

### 1.3. Uzyskanie szczegółowych informacji o module

Aby zobaczyć szczegółowe informacje o konkretnym module z wiersza poleceń, użyj polecenia `info`:

```bash
nf-core modules info cat/cat
```

To wyświetla dokumentację modułu, w tym jego wejścia, wyjścia i podstawowe informacje o użyciu.

??? success "Wynik polecenia"

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

### 1.4. Instalacja modułu cat/cat

Teraz, gdy znaleźliśmy moduł, którego chcemy, musimy dodać go do kodu źródłowego naszego pipeline.

Dobra wiadomość jest taka, że projekt nf-core zawiera narzędzia, które ułatwiają tę część.
Konkretnie, polecenie `nf-core modules install` umożliwia zautomatyzowanie pobierania kodu i udostępnienia go projektowi w jednym kroku.

Przejdź do katalogu Swojego pipeline i uruchom polecenie instalacji:

```bash
cd core-hello
nf-core modules install cat/cat
```

Narzędzie może najpierw poprosić Cię o określenie typu repozytorium.
(Jeśli nie, przejdź do "Na końcu narzędzie przystąpi do instalacji modułu.")

??? success "Wynik polecenia"

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

??? success "Wynik polecenia"

    ```console
        INFO     To avoid this prompt in the future, add the 'repository_type' key to your .nf-core.yml file.
        ? Would you like me to add this config now? [y/n] (y):
    ```

Warto skorzystać z tego wygodnego narzędzia!
Naciśnij enter, aby zaakceptować domyślną odpowiedź (tak).

Na końcu narzędzie przystąpi do instalacji modułu.

??? success "Wynik polecenia"

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
- Dostarcza Ci prawidłową instrukcję `include` do użycia w Swoim workflow'ie

!!! tip

    Zawsze upewnij się, że Twój bieżący katalog roboczy to katalog główny projektu pipeline przed uruchomieniem polecenia instalacji modułu.

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

Możesz również zweryfikować instalację, prosząc narzędzie nf-core o wyświetlenie lokalnie zainstalowanych modułów:

```bash
nf-core modules list local
```

??? success "Wynik polecenia"

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

Jednak aby faktycznie użyć nowego modułu, musimy go zaimportować do naszego pipeline.

### 1.5. Aktualizacja importów modułów

Zastąpmy instrukcję `include` dla modułu `collectGreetings` instrukcją dla `CAT_CAT` w sekcji importów pliku workflow `workflows/hello.nf`.

Przypominając, narzędzie instalacji modułu podało nam dokładną instrukcję do użycia:

```groovy title="Instrukcja importu wygenerowana przez polecenie instalacji"
include { CAT_CAT } from '../modules/nf-core/cat/cat/main'`
```

Zauważ, że konwencja nf-core polega na użyciu wielkich liter dla nazw modułów podczas ich importowania.

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

Zauważ, jak ścieżka dla modułu nf-core różni się od komponentów lokalnych:

- **Komponent nf-core**: `'../modules/nf-core/cat/cat/main'` (odniesienie do `main.nf`)
- **Komponent lokalny**: `'../modules/local/collectGreetings.nf'` (odniesienie do pojedynczego pliku)

CAT_CAT jest teraz dostępny dla workflow, więc wszystko, co musimy zrobić, to zamienić wywołanie `collectGreetings` na jego użycie. Prawda?

Nie tak szybko.

W tym momencie możesz być kuszony, aby zacząć edytować kod, ale warto poświęcić chwilę, aby dokładnie sprawdzić, czego oczekuje nowy moduł i co produkuje.

Zajmiemy się tym jako osobną sekcją, ponieważ obejmuje to nowy mechanizm, którego jeszcze nie omówiliśmy: mapy metadanych.

!!! note

    Opcjonalnie możesz usunąć plik `collectGreetings.nf`:

    ```bash
    rm modules/local/collectGreetings.nf
    ```

    Możesz jednak chcieć go zachować jako punkt odniesienia do zrozumienia różnic między modułami lokalnymi a modułami nf-core.

### Podsumowanie

Wiesz, jak znaleźć moduł nf-core i udostępnić go Swojemu projektowi.

### Co dalej?

Oceń, czego wymaga nowy moduł i zidentyfikuj wszelkie ważne zmiany potrzebne do zintegrowania go z pipeline'em.

---

## 2. Ocena wymagań nowego modułu

Konkretnie, musimy zbadać **interfejs** modułu, tj. jego definicje wejść i wyjść, i porównać go z interfejsem modułu, który chcemy zastąpić.
To pozwoli nam określić, czy możemy po prostu traktować nowy moduł jako zamiennik typu "drop-in", czy też będziemy musieli dostosować część połączeń.

Najlepiej byłoby zrobić to _przed_ zainstalowaniem modułu, ale hej, lepiej późno niż wcale.
(Na marginesie, istnieje polecenie `uninstall`, aby pozbyć się modułów, których nie chcesz już używać.)

!!! note

    Proces CAT_CAT zawiera dość sprytne obsługiwanie różnych typów kompresji, rozszerzeń plików itp., które nie są ściśle istotne dla tego, co próbujemy Ci tutaj pokazać, więc zignorujemy większość z tego i skupimy się tylko na częściach, które są ważne.

### 2.1. Porównanie interfejsów dwóch modułów

Przypominając, tak wygląda interfejs naszego modułu `collectGreetings`:

```groovy title="modules/local/collectGreetings.nf (fragment)" linenums="1" hl_lines="6-7 10"
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

Po zakończeniu `collectGreetings` wyprowadza pojedynczą ścieżkę pliku, emitowaną z tagiem `outfile`.

W porównaniu, interfejs modułu `cat/cat` jest bardziej złożony:

```groovy title="modules/nf-core/cat/cat/main.nf (fragment)" linenums="1" hl_lines="11 14"
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

CAT_CAT przyjmuje pojedyncze wejście, ale jest to krotka składająca się z dwóch elementów:

- `meta` to struktura przechowująca metadane, nazywana metamapą;
- `files_in` to jeden lub więcej plików do przetworzenia, równoważne `input_files` z `collectGreetings`.

Po zakończeniu CAT_CAT dostarcza swoje wyjścia w dwóch częściach:

- Kolejna krotka zawierająca metamapę i połączony plik wyjściowy, emitowana z tagiem `file_out`;
- Plik `versions.yml`, który przechwytuje informacje o wersji oprogramowania, które zostało użyte, emitowany z tagiem `versions`.

Zauważ również, że domyślnie plik wyjściowy będzie nazwany na podstawie identyfikatora będącego częścią metadanych (kod nie pokazany tutaj).

Może to wydawać się dużo do śledzenia, patrząc tylko na kod, więc oto diagram, który pomoże Ci zwizualizować, jak wszystko do siebie pasuje.

<figure class="excalidraw">
--8<-- "docs/en/docs/hello_nf-core/img/module_comparison.svg"
</figure>

Widać, że oba moduły mają podobne wymagania wejściowe pod względem zawartości (zestaw plików wejściowych plus niektóre metadane), ale bardzo różne oczekiwania co do sposobu pakowania tej zawartości.
Ignorując na razie plik wersji, ich główne wyjście jest również równoważne (połączony plik), z wyjątkiem tego, że CAT_CAT emituje również metamapę w połączeniu z plikiem wyjściowym.

Różnice w pakowaniu będą dość łatwe do obsłużenia, jak zobaczysz za chwilę.
Jednak aby zrozumieć część z metamapą, musimy przedstawić Ci dodatkowy kontekst.

### 2.2. Zrozumienie metamap

Właśnie powiedzieliśmy, że moduł CAT_CAT oczekuje mapy metadanych jako części Swojej krotki wejściowej.
Poświęćmy kilka minut na bliższe przyjrzenie się temu, czym to jest.

**Mapa metadanych**, często nazywana w skrócie **metamapą**, to mapa w stylu Groovy zawierająca informacje o jednostkach danych.
W kontekście pipeline Nextflow jednostki danych mogą być czymkolwiek: pojedynczymi próbkami, partiami próbek lub całymi zbiorami danych.

Zgodnie z konwencją, metamapa nf-core jest nazywana `meta` i zawiera wymagane pole `id`, które jest używane do nazywania wyjść i śledzenia jednostek danych.

Na przykład, typowa mapa metadanych może wyglądać tak:

```groovy title="Przykład metamapy na poziomie próbki"
[id: 'sample1', single_end: false, strandedness: 'forward']
```

Lub w przypadku, gdy metadane są dołączone na poziomie partii:

```groovy title="Przykład metamapy na poziomie partii"
[id: 'batch1', date: '25.10.01']
```

Teraz umieśćmy to w kontekście procesu `CAT_CAT`, który oczekuje, że pliki wejściowe będą zapakowane w krotkę z metamapą, i również wyprowadza metamapę jako część krotki wyjściowej.

```groovy title="modules/nf-core/cat/cat/main.nf (fragment)" linenums="1" hl_lines="2 5"
input:
tuple val(meta), path(files_in)

output:
tuple val(meta), path("${prefix}"), emit: file_out
```

W rezultacie każda jednostka danych przemieszcza się przez pipeline z dołączonymi odpowiednimi metadanymi.
Kolejne procesy mogą następnie również łatwo uzyskać dostęp do tych metadanych.

Pamiętasz, jak mówiliśmy, że plik wyprowadzany przez `CAT_CAT` będzie nazwany na podstawie identyfikatora będącego częścią metadanych?
To jest odpowiedni kod:

```groovy title="modules/nf-core/cat/cat/main.nf (fragment)" linenums="35"
prefix   = task.ext.prefix ?: "${meta.id}${getFileSuffix(file_list[0])}"
```

Oznacza to mniej więcej następująco: jeśli `prefix` jest dostarczony przez system parametrów zewnętrznych zadania (`task.ext`), użyj tego do nazwania pliku wyjściowego; w przeciwnym razie utwórz go używając `${meta.id}`, który odpowiada polu `id` w metamapie.

Możesz sobie wyobrazić kanał wejściowy wchodzący do tego modułu z zawartością taką jak ta:

```groovy title="Przykład zawartości kanału wejściowego"
ch_input = [[[id: 'batch1', date: '25.10.01'], ['file1A.txt', 'file1B.txt']],
            [[id: 'batch2', date: '25.10.26'], ['file2A.txt', 'file2B.txt']],
            [[id: 'batch3', date: '25.11.14'], ['file3A.txt', 'file3B.txt']]]
```

Następnie zawartość kanału wyjściowego wychodzącego wyglądałaby tak:

```groovy title="Przykład zawartości kanału wyjściowego"
ch_input = [[[id: 'batch1', date: '25.10.01'], 'batch1.txt'],
            [[id: 'batch2', date: '25.10.26'], 'batch2.txt'],
            [[id: 'batch3', date: '25.11.14'], 'batch3.txt']]
```

Jak wspomniano wcześniej, konfiguracja wejściowa `tuple val(meta), path(files_in)` jest standardowym wzorcem używanym we wszystkich modułach nf-core.

Miejmy nadzieję, że zaczynasz widzieć, jak przydatne może to być.
Nie tylko pozwala to nazwać wyjścia na podstawie metadanych, ale możesz również robić takie rzeczy, jak używać ich do stosowania różnych wartości parametrów, a w połączeniu z określonymi operatorami możesz nawet grupować, sortować lub filtrować dane, gdy przepływają przez pipeline.

!!! note "Dowiedz się więcej o metadanych"

    Aby uzyskać kompleksowe wprowadzenie do pracy z metadanymi w workflow'ach Nextflow, w tym jak odczytywać metadane z arkuszy próbek i używać ich do dostosowywania przetwarzania, zobacz side quest [Metadane w workflow'ach](../side_quests/metadata).

### 2.3. Podsumowanie zmian do wprowadzenia

Na podstawie tego, co przejrzeliśmy, oto główne zmiany, które musimy wprowadzić w naszym pipeline'ie, aby wykorzystać moduł `cat/cat`:

- Utworzyć metamapę zawierającą nazwę partii;
- Zapakować metamapę w krotkę z zestawem plików wejściowych do połączenia (wychodzących z `convertToUpper`);
- Zmienić wywołanie z `collectGreetings()` na `CAT_CAT`;
- Wyodrębnić plik wyjściowy z krotki wytworzonej przez proces `CAT_CAT` przed przekazaniem go do `cowpy`.

To powinno załatwić sprawę! Teraz, gdy mamy plan, jesteśmy gotowi do działania.

### Podsumowanie

Wiesz, jak ocenić interfejs wejściowy i wyjściowy nowego modułu, aby zidentyfikować jego wymagania. Nauczyłeś się również, jak metamapy są używane przez nf-core do utrzymywania informacji kontekstowych ściśle powiązanych z danymi podczas ich przepływu przez workflow.

### Co dalej?

Zintegrować nowy moduł z workflow'iem.

---

## 3. Integracja CAT_CAT z workflow'iem `hello.nf`

Teraz, gdy wiesz wszystko o metamapach (lub wystarczająco dużo dla celów tego kursu), nadszedł czas, aby faktycznie zaimplementować zmiany, które opisaliśmy powyżej.

Dla jasności podzielimy to i omówimy każdy krok osobno.

!!! note

    Wszystkie zmiany pokazane poniżej są dokonywane w logice workflow'u w bloku `main` w pliku workflow'u `core-hello/workflows/hello.nf`.

### 3.1. Utworzenie mapy metadanych

Najpierw musimy utworzyć mapę metadanych dla `CAT_CAT`, pamiętając, że moduły nf-core wymagają, aby metamapa zawierała co najmniej pole `id`.

Ponieważ nie potrzebujemy innych metadanych, możemy to uprościć i użyć czegoś takiego:

```groovy title="Przykład składni"
def cat_meta = [id: 'test']
```

Z tym że nie chcemy sztywno kodować wartości `id`; chcemy użyć wartości parametru `params.batch`.
Więc kod staje się:

```groovy title="Przykład składni"
def cat_meta = [id: params.batch]
```

Tak, jest to dosłownie tak proste, aby utworzyć podstawową metamapę.

Dodajmy te linie po wywołaniu `convertToUpper`, usuwając wywołanie `collectGreetings`:

=== "Po"

    ```groovy title="core-hello/workflows/hello.nf" linenums="26" hl_lines="7-8"
        // wyemituj pozdrowienie
        sayHello(ch_samplesheet)

        // przekształć pozdrowienie na wielkie litery
        convertToUpper(sayHello.out)

        // utwórz mapę metadanych z nazwą partii jako ID
        def cat_meta = [ id: params.batch ]

        // wygeneruj grafikę ASCII pozdrowień za pomocą cowpy
        cowpy(collectGreetings.out.outfile, params.character)
    ```

=== "Przed"

    ```groovy title="core-hello/workflows/hello.nf" linenums="26" hl_lines="7-8"
        // wyemituj pozdrowienie
        sayHello(ch_samplesheet)

        // przekształć pozdrowienie na wielkie litery
        convertToUpper(sayHello.out)

        // zbierz wszystkie powitania w jeden plik
        collectGreetings(convertToUpper.out.collect(), params.batch)

        // wygeneruj grafikę ASCII pozdrowień za pomocą cowpy
        cowpy(collectGreetings.out.outfile, params.character)
    ```

To tworzy prostą mapę metadanych, gdzie `id` jest ustawione na naszą nazwę partii (która będzie `test` przy użyciu profilu testowego).

### 3.2. Utworzenie kanału z krotkami metadanych

Następnie przekształć kanał plików w kanał krotek zawierających metadane i pliki:

=== "Po"

    ```groovy title="core-hello/workflows/hello.nf" linenums="26" hl_lines="10-11"
        // wyemituj pozdrowienie
        sayHello(ch_samplesheet)

        // przekształć pozdrowienie na wielkie litery
        convertToUpper(sayHello.out)

        // utwórz mapę metadanych z nazwą partii jako ID
        def cat_meta = [ id: params.batch ]

        // utwórz kanał z metadanymi i plikami w formacie krotki
        ch_for_cat = convertToUpper.out.collect().map { files -> tuple(cat_meta, files) }

        // wygeneruj grafikę ASCII pozdrowień za pomocą cowpy
        cowpy(collectGreetings.out.outfile, params.character)
    ```

=== "Przed"

    ```groovy title="core-hello/workflows/hello.nf" linenums="26"
        // wyemituj pozdrowienie
        sayHello(ch_samplesheet)

        // przekształć pozdrowienie na wielkie litery
        convertToUpper(sayHello.out)

        // utwórz mapę metadanych z nazwą partii jako ID
        def cat_meta = [ id: params.batch ]

        // wygeneruj grafikę ASCII pozdrowień za pomocą cowpy
        cowpy(collectGreetings.out.outfile, params.character)
    ```

Linia, którą dodaliśmy, osiąga dwie rzeczy:

- `.collect()` zbiera wszystkie elementy z wyjścia `convertToUpper` w jedną listę
- `.map { files -> tuple(cat_meta, files) }` tworzy krotkę `[metadata, files]` w formacie oczekiwanym przez `CAT_CAT`

To wszystko, co musimy zrobić, aby przygotować dane wejściowe dla `CAT_CAT`.

### 3.3. Wywołanie modułu CAT_CAT

Teraz wywołaj `CAT_CAT` na nowo utworzonym kanale:

=== "Po"

    ```groovy title="core-hello/workflows/hello.nf" linenums="26" hl_lines="13-14"
        // wyemituj pozdrowienie
        sayHello(ch_samplesheet)

        // przekształć pozdrowienie na wielkie litery
        convertToUpper(sayHello.out)

        // utwórz mapę metadanych z nazwą partii jako ID
        def cat_meta = [ id: params.batch ]

        // utwórz kanał z metadanymi i plikami w formacie krotki
        ch_for_cat = convertToUpper.out.collect().map { files -> tuple(cat_meta, files) }

        // połącz pliki używając modułu nf-core cat/cat
        CAT_CAT(ch_for_cat)

        // wygeneruj grafikę ASCII pozdrowień za pomocą cowpy
        cowpy(collectGreetings.out.outfile, params.character)
    ```

=== "Przed"

    ```groovy title="core-hello/workflows/hello.nf" linenums="26"
        // wyemituj pozdrowienie
        sayHello(ch_samplesheet)

        // przekształć pozdrowienie na wielkie litery
        convertToUpper(sayHello.out)

        // utwórz mapę metadanych z nazwą partii jako ID
        def cat_meta = [ id: params.batch ]

        // utwórz kanał z metadanymi i plikami w formacie krotki
        ch_for_cat = convertToUpper.out.collect().map { files -> tuple(cat_meta, files) }

        // wygeneruj grafikę ASCII pozdrowień za pomocą cowpy
        cowpy(collectGreetings.out.outfile, params.character)
    ```

To kończy najtrudniejszą część tej zamiany, ale jeszcze nie skończyliśmy: nadal musimy zaktualizować sposób, w jaki przekazujemy połączone wyjście do procesu `cowpy`.

### 3.4. Wyodrębnienie pliku wyjściowego z krotki dla `cowpy`

Wcześniej proces `collectGreetings` po prostu produkował plik, który mogliśmy przekazać bezpośrednio do `cowpy`.
Jednak proces `CAT_CAT` produkuje krotkę, która zawiera metamapę oprócz pliku wyjściowego.

Ponieważ `cowpy` nie akceptuje jeszcze krotek metadanych (naprawimy to w następnej części kursu), musimy wyodrębnić plik wyjściowy z krotki wytworzonej przez `CAT_CAT` przed przekazaniem go do `cowpy`:

=== "Po"

    ```groovy title="core-hello/workflows/hello.nf" linenums="26" hl_lines="16-17 20"
        // wyemituj pozdrowienie
        sayHello(ch_samplesheet)

        // przekształć pozdrowienie na wielkie litery
        convertToUpper(sayHello.out)

        // utwórz mapę metadanych z nazwą partii jako ID
        def cat_meta = [ id: params.batch ]

        // utwórz kanał z metadanymi i plikami w formacie krotki
        ch_for_cat = convertToUpper.out.collect().map { files -> tuple(cat_meta, files) }

        // połącz powitania
        CAT_CAT(ch_for_cat)

        // wyodrębnij plik z krotki, ponieważ cowpy nie używa jeszcze metadanych
        ch_for_cowpy = CAT_CAT.out.file_out.map{ meta, file -> file }

        // wygeneruj grafikę ASCII z powitaniami za pomocą cowpy
        cowpy(ch_for_cowpy, params.character)
    ```

=== "Przed"

    ```groovy title="core-hello/workflows/hello.nf" linenums="26" hl_lines="17"
        // wyemituj pozdrowienie
        sayHello(ch_samplesheet)

        // przekształć pozdrowienie na wielkie litery
        convertToUpper(sayHello.out)

        // utwórz mapę metadanych z nazwą partii jako ID
        def cat_meta = [ id: params.batch ]

        // utwórz kanał z metadanymi i plikami w formacie krotki
        ch_for_cat = convertToUpper.out.collect().map { files -> tuple(cat_meta, files) }

        // połącz powitania
        CAT_CAT(ch_for_cat)

        // wygeneruj grafikę ASCII z powitaniami za pomocą cowpy
        cowpy(collectGreetings.out.outfile, params.character)
    ```

Operacja `.map{ meta, file -> file }` wyodrębnia plik z krotki `[metadata, file]` wytworzonej przez `CAT_CAT` do nowego kanału, `ch_for_cowpy`.

Następnie wystarczy przekazać `ch_for_cowpy` do `cowpy` zamiast `collectGreetings.out.outfile` w tej ostatniej linii.

!!! note

    W następnej części kursu zaktualizujemy `cowpy`, aby pracował bezpośrednio z krotkami metadanych, więc ten krok ekstrakcji nie będzie już potrzebny.

### 3.5. Testowanie workflow

Przetestujmy, czy workflow działa z nowo zintegrowanym modułem `cat/cat`:

```bash
nextflow run . --outdir core-hello-results -profile test,docker --validate_params false
```

To powinno działać dość szybko.

??? success "Wynik polecenia"

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

I to wszystko! Teraz używamy solidnego, utrzymywanego przez społeczność modułu zamiast niestandardowego kodu na poziomie prototypu dla tego kroku w pipeline'ie.

### Podsumowanie

Teraz wiesz, jak:

- Znaleźć i zainstalować moduły nf-core
- Ocenić wymagania modułu nf-core
- Utworzyć prostą mapę metadanych do użycia z modułem nf-core
- Zintegrować moduł nf-core ze Swoim workflow'em

### Co dalej?

Naucz się dostosowywać Swoje moduły lokalne, aby przestrzegały konwencji nf-core.
Pokażemy Ci również, jak tworzyć nowe moduły nf-core z szablonu za pomocą narzędzi nf-core.
