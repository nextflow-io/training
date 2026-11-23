# Część 2: Konfiguracja pipeline'u

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Tłumaczenie wspomagane przez AI - [dowiedz się więcej i zasugeruj ulepszenia](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

W [Części 1](./01_run_nextflow.md) uruchomiłeś/uruchomiłaś kompletny, wieloetapowy pipeline, który przetwarza wiele danych wejściowych równolegle z użyciem kontenerów.
Teraz przyjrzymy się, jak konfigurować zachowanie pipeline'u przy użyciu `nextflow.config`: najpierw przeanalizujemy plik konfiguracyjny, który już otrzymałeś/otrzymałaś, następnie zbadamy kilka innych sposobów dostarczania konfiguracji, a na końcu zajmiemy się kontrolowaniem tego, jak i gdzie publikowane są wyniki.

---

## 1. Analiza głównego pliku konfiguracyjnego

Nextflow automatycznie wczytuje `nextflow.config` z bieżącego katalogu roboczego i stosuje jego ustawienia do każdego uruchomienia.

Dostarczamy Ci plik konfiguracyjny obejmujący cztery obszary: pakowanie oprogramowania, ustawienia procesów, parametry pipeline'u oraz profile wykonania.

??? full-code "nextflow.config"

    ```groovy title="nextflow.config" linenums="1"
    /*
     * Pakowanie oprogramowania
     */
    docker.enabled = true

    /*
     * Ustawienia procesów
     */
    process {
        cpus = 1
        memory = 1.GB
    }

    /*
     * Parametry pipeline'u
     */
    params {
        input = 'data/greetings.csv'
        batch = 'batch'
        character = 'turkey'
    }

    /*
     * Profile
     */
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

Omówimy każdy z tych obszarów, a następnie wykorzystamy profile w praktyce, uruchamiając pipeline z jednym z nich.

!!! note "Uwaga"

    Ta konfiguracja dotyczy lokalnego wykonania na pojedynczej maszynie.
    Nextflow obsługuje również harmonogramy HPC (SLURM, PBS, LSF) oraz executory chmurowe (AWS Batch, Google Cloud Batch, Azure Batch) — wszystkie konfigurowane za pomocą tego samego mechanizmu `nextflow.config`.
    Pełny opis tych opcji znajdziesz w [Części 1: Dostosowanie do środowiska obliczeniowego](../config_exec/01_packaging_and_execution.md) kursu [Configure Execution](../config_exec/index.md).

### 1.1. Pakowanie oprogramowania

Pakowanie oprogramowania to sposób, w jaki Nextflow dostarcza narzędzia potrzebne Twoim procesom — czy to obraz kontenera, środowisko Conda, czy coś innego.

```groovy title="nextflow.config" linenums="1"
/*
 * Pakowanie oprogramowania
 */
docker.enabled = true
```

Ta linia włącza Docker dla każdego procesu.
Każdy proces, który deklaruje dyrektywę `container`, uruchamia się wewnątrz wskazanego obrazu.

### 1.2. Ustawienia procesów

Pamiętaj, że proces to pojedynczy krok w Twoim pipeline'ie, taki jak `sayHello` czy `cowpy`.
Nextflow pozwala skonfigurować wiele aspektów działania każdego z nich: ile CPU i pamięci otrzymuje, jakiego kontenera lub środowiska Conda używa i wiele więcej.

```groovy title="nextflow.config" linenums="6"
/*
 * Ustawienia procesów
 */
process {
    cpus = 1
    memory = 1.GB
}
```

To ogranicza każdy proces do jednego CPU i 1 GB pamięci.

Nextflow pozwala też ustawiać różne wartości dla poszczególnych nazwanych procesów lub grup procesów — dowiesz się, jak to zrobić, w [Części 2: Zarządzanie zasobami obliczeniowymi i błędami](../config_exec/02_resources_and_retries.md#12-set-resource-allocations-for-a-specific-process) kursu [Configure Execution](../config_exec/index.md).

### 1.3. Parametry pipeline'u

Parametry to dane wejściowe pipeline'u przekazywane z wiersza poleceń — te same flagi `--input`, `--batch` i `--character`, które już ustawiałeś/ustawiałaś bezpośrednio w wierszu poleceń.
Zdefiniowanie ich wartości domyślnych tutaj oznacza, że nie musisz ich wpisywać za każdym razem. Jak zobaczysz w dalszej części, istnieje też kilka innych sposobów ich dostarczania.

```groovy title="nextflow.config" linenums="14"
/*
 * Parametry pipeline'u
 */
params {
    input = 'data/greetings.csv'
    batch = 'batch'
    character = 'turkey'
}
```

Te wartości domyślne są stosowane zawsze, gdy parametr nie zostanie podany w wierszu poleceń — dzięki temu uruchomienie `nextflow run main.nf` bez żadnych flag nadal działa.

### 1.4. Profile

Profile pozwalają zgrupować zestaw ustawień pod jedną nazwą, dzięki czemu możesz przełączać się między całymi konfiguracjami za pomocą jednej flagi, zamiast ręcznie zmieniać wartości za każdym razem.

```groovy title="nextflow.config" linenums="23"
/*
 * Profile
 */
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

Profil `test` nadpisuje trzy parametry, aby uruchomić pipeline z małym, dobrze zdefiniowanym zestawem danych wejściowych. Każdy pipeline nf-core zawiera taki profil do szybkiej walidacji — to konwencja warta stosowania również we własnych pipeline'ach.

Profil `conda` przełącza pakowanie oprogramowania z Docker na Conda.

Profil aktywuje się, przekazując `-profile <nazwa>` w wierszu poleceń.

Wypróbujmy profil `test` w praktyce.

```bash
nextflow run main.nf -profile test
```

??? success "Wyjście polecenia"

    ```console hl_lines="6 7 8 9"
    N E X T F L O W   ~  version 26.04.4

    Launching `main.nf` [reverent_heisenberg] revision: ce74f81996

    executor >  local (8)
    [3d/8a12c7] sayHello (3)       | 3 of 3 ✔
    [e7/e0934f] convertToUpper (2) | 3 of 3 ✔
    [1d/616569] collectGreetings   | 1 of 1 ✔
    [44/7d46cf] cowpy              | 1 of 1 ✔

    Outputs:

      /workspaces/training/nextflow-run/results

      first_output:
        - test/intermediates/Bonjour-output.txt
        - test/intermediates/Hello-output.txt
        - test/intermediates/Hola-output.txt

      uppercased:
        - test/intermediates/UPPER-Bonjour-output.txt
        - test/intermediates/UPPER-Hello-output.txt
        - test/intermediates/UPPER-Hola-output.txt

      collected: test/intermediates/COLLECTED-test-output.txt

      batch_report: test/test-report.txt

      cowpy_art: test/cowpy-COLLECTED-test-output.txt
    ```

Pipeline uruchamia się z `batch = 'test'` i `character = 'tux'`.
Sprawdź katalog `results/test/`: nazwa partii jest teraz częścią ścieżki katalogu, a grafika ASCII przedstawia pingwina tux zamiast indyka.

!!! note "Uwaga"

    Możesz aktywować kilka profili jednocześnie i użyć `nextflow config -profile <nazwa>,<nazwa>`, aby zobaczyć w pełni rozwiązany wynik przed uruchomieniem czegokolwiek.
    Łączenie profili i sposób, w jaki Nextflow rozwiązuje konflikty między nimi, jest szczegółowo omówione w [Części 3: Używanie profili do przełączania konfiguracji](../config_exec/03_profiles.md) kursu [Configure Execution](../config_exec/index.md).

### Podsumowanie

Wiesz już, do czego służą najczęstsze elementy pliku `nextflow.config` i jak aktywować profil.

### Co dalej?

Poznaj kilka innych sposobów dostarczania wartości konfiguracyjnych bez modyfikowania głównego pliku `nextflow.config` — przydatnych do konfigurowania pojedynczych uruchomień oraz do udostępniania dokładnego zestawu ustawień innym osobom.

---

## 2. Dostarczanie konfiguracji za pomocą plików uzupełniających

Ustawianie wartości domyślnych w `nextflow.config` sprawdza się dobrze dla wartości, które rzadko się zmieniają.
Nextflow oferuje też dwa bardziej ukierunkowane mechanizmy: plik konfiguracyjny specyficzny dla danego uruchomienia, służący do dostosowania wykonania do konkretnego środowiska, oraz plik parametrów do udostępniania dokładnego zestawu wartości wejściowych współpracownikom.

### 2.1. Użycie pliku konfiguracyjnego specyficznego dla uruchomienia

Załóżmy, że przenosisz pipeline na maszynę bez Docker i chcesz dać każdemu procesowi więcej zasobów.
Utwórz nowy plik konfiguracyjny zawierający tylko potrzebne nadpisania:

```groovy title="custom.config" linenums="1"
process {
    cpus = 2
    memory = 2.GB
}

docker.enabled = false
conda.enabled = true
```

Przekaż go razem z głównym pipeline'em za pomocą flagi `-c`:

```bash
nextflow run main.nf -c custom.config
```

??? success "Wyjście polecenia"

    ```console hl_lines="6 7 8 9"
    N E X T F L O W   ~  version 26.04.4

    Launching `main.nf` [exotic_cray] revision: ce74f81996

    executor >  local (8)
    [77/e02315] sayHello (1)       | 3 of 3 ✔
    [a6/ccf44b] convertToUpper (3) | 3 of 3 ✔
    [56/fd1296] collectGreetings   | 1 of 1 ✔
    [2c/205a94] cowpy              | 1 of 1 ✔

    Outputs:

      /workspaces/training/nextflow-run/results

      first_output:
        - batch/intermediates/Bonjour-output.txt
        - batch/intermediates/Hello-output.txt
        - batch/intermediates/Hola-output.txt

      uppercased:
        - batch/intermediates/UPPER-Bonjour-output.txt
        - batch/intermediates/UPPER-Hello-output.txt
        - batch/intermediates/UPPER-Hola-output.txt

      collected: batch/intermediates/COLLECTED-batch-output.txt

      batch_report: batch/batch-report.txt

      cowpy_art: batch/cowpy-COLLECTED-batch-output.txt
    ```

Nextflow scala `custom.config` z własnym `nextflow.config` pipeline'u, więc każdy proces otrzymuje teraz 2 CPU i 2 GB pamięci zamiast wartości domyślnych i działa przez Conda zamiast Docker.
`cowpy` to jedyny proces z zadeklarowanym pakietem Conda obok kontenera, więc to właśnie dla niego Nextflow faktycznie zbuduje środowisko:

```console
Creating env using conda: conda-forge::cowpy==1.1.5 [cache /path/to/work/conda/env-898314d566668b6587ad714ae06b8520]
```

Mały plik nadpisujący tylko alokację zasobów i pakowanie, bez dotykania parametrów pipeline'u, to dokładnie wzorzec, którego pipeline'y nf-core oczekują od konfiguracji instytucjonalnych.
Przejrzyj repozytorium [nf-core/configs](https://github.com/nf-core/configs), aby zobaczyć przykłady z rzeczywistych wdrożeń.

To wygodny sposób na dostosowanie pipeline'u do nowego środowiska bez ingerowania w normalną konfigurację.

### 2.2. Użycie pliku parametrów

Załóżmy, że musisz udostępnić dokładny zestaw parametrów uruchomienia współpracownikowi lub zapisać je na potrzeby publikacji.

Nextflow pozwala dostarczać [pliki parametrów](https://nextflow.io/docs/latest/config.html#parameter-file) w formacie YAML lub JSON — to prostszy sposób na dystrybucję dokładnego, odtwarzalnego zestawu wartości.

W Twoim bieżącym katalogu roboczym znajduje się już plik parametrów `test-params.yaml`:

```yaml title="test-params.yaml" linenums="1"
input: "data/greetings.csv"
batch: "yaml"
character: "stegosaurus"
```

Składnia używa dwukropków (`:`) zamiast znaków równości (`=`) stosowanych w `nextflow.config`, ponieważ ten plik jest zwykłym YAML, a nie Groovy.

!!! info "Info"

    Dostępna jest również wersja JSON: `test-params.json`. Możesz ją wypróbować samodzielnie — składnia przekazywania jest identyczna.

Przekaż plik za pomocą `-params-file`:

```bash
nextflow run main.nf -params-file test-params.yaml
```

??? success "Wyjście polecenia"

    ```console
    N E X T F L O W   ~  version 26.04.4

    Launching `main.nf` [sharp_faraday] revision: ce74f81996

    executor >  local (8)
    [1c/9ff63e] sayHello (1)       | 3 of 3 ✔
    [3b/bb5691] convertToUpper (2) | 3 of 3 ✔
    [cd/2c1f6e] collectGreetings   | 1 of 1 ✔
    [89/c333bc] cowpy              | 1 of 1 ✔

    Outputs:

      /workspaces/training/nextflow-run/results

      first_output:
        - yaml/intermediates/Bonjour-output.txt
        - yaml/intermediates/Hello-output.txt
        - yaml/intermediates/Hola-output.txt

      uppercased:
        - yaml/intermediates/UPPER-Bonjour-output.txt
        - yaml/intermediates/UPPER-Hello-output.txt
        - yaml/intermediates/UPPER-Hola-output.txt

      collected: yaml/intermediates/COLLECTED-yaml-output.txt

      batch_report: yaml/yaml-report.txt

      cowpy_art: yaml/cowpy-COLLECTED-yaml-output.txt
    ```

??? abstract "Zawartość pliku"

    ```console title="results/yaml/cowpy-COLLECTED-yaml-output.txt"
     _________
    / BONJOUR \
    | HOLA    |
    \ HELLO   /
     ---------
    \                             .       .
     \                           / `.   .' "
      \                  .---.  <    > <    >  .---.
       \                 |    \  \ - ~ ~ - /  /    |
             _____          ..-~             ~-..-~
            |     |   \~~~\.'                    `./~~~/
           ---------   \__/                        \__/
          .'  O    \     /               /       \  "
         (_____,    `._.'               |         }  \/~~~/
          `----.          /       }     |        /    \__/
                `-.      |       /      |       /      `. ,~~|
                    ~-.__|      /_ - ~ ^|      /- _      `..-'
                        |     /        |     /     ~-.     `-. _  _  _
                        |_____|        |_____|         ~ - . _ _ _ _ _>
    ```

Plik parametrów jest szczególnie wartościowy, gdy pipeline ma więcej niż kilka parametrów: pozwala dostarczyć je wszystkie naraz, bez rozbudowanego wiersza poleceń i bez żadnych zmian w skrypcie workflow'u, a jego dystrybucja razem z wynikami jest prosta.

### Podsumowanie

Znasz już dwa kolejne sposoby dostarczania konfiguracji: plik konfiguracyjny specyficzny dla uruchomienia, służący do dostosowania wykonania do nowego środowiska, oraz plik parametrów do udostępniania dokładnych, odtwarzalnych wartości wejściowych.

### Co dalej?

Dowiedz się, jak kontrolować sposób i miejsce publikowania wyników pipeline'u.

---

## 3. Zarządzanie wynikami pipeline'u

Autor pipeline'u decyduje w kodzie o organizacji wyników, ale nie musisz dotykać tego kodu, aby kontrolować, gdzie trafiają i w jaki sposób.
Nextflow oferuje do tego mechanizmy na poziomie konfiguracji: ustawienie bazowego katalogu wyjściowego oraz wybór między kopiowaniem a dowiązaniami symbolicznymi plików.

### 3.1. Dostosowanie katalogu wyjściowego

Domyślnie Nextflow publikuje wyniki w katalogu `results/`.
Możesz wskazać inną lokalizację za pomocą `-output-dir` (lub skróconej formy `-o`):

```bash
nextflow run main.nf -output-dir outputs
```

??? success "Wyjście polecenia"

    ```console hl_lines="6 7 8 9"
    N E X T F L O W   ~  version 26.04.4

    Launching `main.nf` [serene_kimura] revision: ce74f81996

    executor >  local (8)
    [31/df5c15] sayHello (3)       | 3 of 3 ✔
    [08/8bb2e5] convertToUpper (1) | 3 of 3 ✔
    [e5/5814da] collectGreetings   | 1 of 1 ✔
    [cf/8ab8c6] cowpy              | 1 of 1 ✔

    Outputs:

      /workspaces/training/nextflow-run/outputs

      first_output:
        - batch/intermediates/Bonjour-output.txt
        - batch/intermediates/Hello-output.txt
        - batch/intermediates/Hola-output.txt

      uppercased:
        - batch/intermediates/UPPER-Bonjour-output.txt
        - batch/intermediates/UPPER-Hello-output.txt
        - batch/intermediates/UPPER-Hola-output.txt

      collected: batch/intermediates/COLLECTED-batch-output.txt

      batch_report: batch/batch-report.txt

      cowpy_art: batch/cowpy-COLLECTED-batch-output.txt
    ```

??? abstract "Zawartość katalogu"

    ```console
    outputs/batch
    ├── batch-report.txt
    ├── cowpy-COLLECTED-batch-output.txt
    └── intermediates
        ├── Bonjour-output.txt
        ├── COLLECTED-batch-output.txt
        ├── Hello-output.txt
        ├── Hola-output.txt
        ├── UPPER-Bonjour-output.txt
        ├── UPPER-Hello-output.txt
        └── UPPER-Hola-output.txt
    ```

Wyniki trafiają teraz do `outputs/batch/` zamiast domyślnego `results/batch/`.
Własny kod pipeline'u nadal decyduje o strukturze wewnątrz tego bazowego katalogu — takich jak podkatalogi `batch/` i `intermediates/`. Flaga `-output-dir` kontroluje jedynie punkt startowy tej struktury.

`-output-dir` to tak naprawdę skrót wiersza poleceń dla opcji konfiguracyjnej `outputDir`, więc można jej użyć wszędzie tam, gdzie można umieścić konfigurację: bezpośrednio w `nextflow.config`, wewnątrz profilu lub w pliku nakładki `-c`, takim jak ten użyty wcześniej w tej części.
Na przykład poniższy fragment pokazuje to samo ustawienie umieszczone bezpośrednio w `nextflow.config` zamiast przekazywane w wierszu poleceń:

```groovy title="nextflow.config"
outputDir = 'outputs'
```

Pełną listę miejsc, w których może znajdować się taka opcja konfiguracyjna, znajdziesz w sekcji [Configuration file](https://nextflow.io/docs/latest/config.html) dokumentacji referencyjnej Nextflow'a.

### 3.2. Wybór sposobu publikowania wyników

Domyślnie Nextflow publikuje wyniki jako dowiązania symboliczne wskazujące na lokalizacje plików w katalogu `work/`, a nie jako prawdziwe kopie:

```console
$ ls -l results/batch/intermediates/Hello-output.txt
lrwxr-xr-x  ...  Hello-output.txt -> /workspaces/training/nextflow-run/work/b7/b8c4d1.../Hello-output.txt
```

Autorzy pipeline'u mogą ustawić tryb publikowania (`'copy'` lub `'move'`) dla każdego procesu z osobna w kodzie workflow'u.
Zazwyczaj robią to dla końcowych wyników pipeline'u, pozostawiając domyślne zachowanie `'symlink'` dla plików pośrednich, które można usunąć po zakończeniu całego pipeline'u.

Takie podejście pozwala uniknąć duplikowania danych na dysku, ale oznacza, że nie można usunąć katalogów zadań w `work/` bez zerwania dowiązań i utraty możliwości korzystania z `-resume`.
Jeśli chcesz, aby wszystkie pliki wyjściowe były właściwie skopiowane, ustaw [`workflow.output.mode`](https://nextflow.io/docs/latest/reference/config.html#workflow) na `'copy'` w konfiguracji pipeline'u. (W przeciwieństwie do `-output-dir`, nie ma dla tego flagi wiersza poleceń — to wyłącznie opcja konfiguracyjna.)

Spróbuj ustawić to w `nextflow.config`:

=== "Po"

    ```groovy title="nextflow.config" hl_lines="7"
    params {
        input = 'data/greetings.csv'
        batch = 'batch'
        character = 'turkey'
    }

    workflow.output.mode = 'copy'
    ```

=== "Przed"

    ```groovy title="nextflow.config"
    params {
        input = 'data/greetings.csv'
        batch = 'batch'
        character = 'turkey'
    }
    ```

Następnie uruchom pipeline, zmieniając nazwę partii, aby zobaczyć różnicę w wynikach:

```bash
nextflow run main.nf --batch withmode
```

??? success "Wyjście polecenia"

    ```console hl_lines="6 7 8 9"
    N E X T F L O W   ~  version 26.04.4

    Launching `main.nf` [angry_noether] revision: ce74f81996

    executor >  local (8)
    [41/2b478d] sayHello (3)       | 3 of 3 ✔
    [bf/dd2840] convertToUpper (2) | 3 of 3 ✔
    [ea/364e97] collectGreetings   | 1 of 1 ✔
    [10/76fe7b] cowpy              | 1 of 1 ✔

    Outputs:

      /workspaces/training/nextflow-run/results

      first_output:
        - withmode/intermediates/Bonjour-output.txt
        - withmode/intermediates/Hello-output.txt
        - withmode/intermediates/Hola-output.txt

      uppercased:
        - withmode/intermediates/UPPER-Bonjour-output.txt
        - withmode/intermediates/UPPER-Hello-output.txt
        - withmode/intermediates/UPPER-Hola-output.txt

      collected: withmode/intermediates/COLLECTED-withmode-output.txt

      batch_report: withmode/withmode-report.txt

      cowpy_art: withmode/cowpy-COLLECTED-withmode-output.txt
    ```

Sprawdź jeden z plików wyjściowych tak jak poprzednio:

```console
$ ls -l results/withmode/intermediates/Hello-output.txt
-rw-r--r--  ...  Hello-output.txt
```

Teraz to prawdziwy, niezależny plik, który pozostanie dostępny nawet po wyczyszczeniu katalogu `work/`.

!!! warning "Ostrzeżenie"

    Ustawienie `workflow.output.mode` wypełnia jedynie wartość domyślną dla wyników, które nie mają jeszcze ustawionego trybu w kodzie pipeline'u.
    Nie może nadpisać trybu zakodowanego na stałe przez autora, niezależnie od tego, co ustawisz.

### Podsumowanie

Wiesz już, jak dostosować bazowy katalog wyjściowy i wybrać między kopiowaniem a dowiązaniami symbolicznymi — bez dotykania kodu pipeline'u.

### Co dalej?

Przejdź do [Części 3](./03_manage_executions.md), gdzie nauczysz się przeglądać historię poprzednich uruchomień, generować raporty wykonania i czyścić stare katalogi robocze.

---

## Podsumowanie

W tej części nauczyłeś/nauczyłaś się:

- Konfigurować zachowanie pipeline'u przy użyciu `nextflow.config` i profili
- Dostarczać konfigurację za pomocą pliku konfiguracyjnego specyficznego dla uruchomienia lub pliku parametrów
- Dostosowywać katalog wyjściowy i wybierać między kopiowaniem a dowiązaniami symbolicznymi wyników
