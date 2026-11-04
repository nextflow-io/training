# Część 2: Konfiguracja uruchomienia pipeline'u

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Tłumaczenie wspomagane przez AI - [dowiedz się więcej i zasugeruj ulepszenia](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

W [Części 1](./01_run_demo.md) znalazłeś i uruchomiłeś pipeline nf-core/demo przy użyciu profilu testowego.
Teraz przyjrzymy się, jak konfigurować uruchomienie pipeline'u: ustawianiu parametrów, walidacji oraz dostosowywaniu alokacji zasobów i argumentów narzędzi.

Jak wyjaśniono w [Hello Config](../hello_nextflow/06_hello_config.md), chcemy mieć możliwość zmiany danych wejściowych i sposobu działania pipeline'u bez modyfikowania jego kodu.
W tym celu Nextflow obsługuje wiele sposobów kontrolowania konfiguracji pipeline'u, co może być nieco przytłaczające.

Projekt nf-core określa konwencje organizowania elementów konfiguracji, rozróżniając dwa rodzaje konfiguracji na najwyższym poziomie: **parametry pipeline'u** oraz **konfigurację** w ścisłym sensie.

- **Parametry pipeline'u** (ustawiane przez system `params`) obejmują zazwyczaj takie elementy jak pliki wejściowe, flagi sterujące zachowaniem narzędzi i parametry analizy.
- **Konfiguracja** w ścisłym sensie odnosi się do logistyki uruchamiania pipeline'u, czyli executora, alokacji zasobów obliczeniowych i tym podobnych.

<figure class="excalidraw">
    --8<-- "docs/en/docs/nfcore_use/img/params_vs_config.excalidraw.svg"
</figure>

Zacznijmy od parametrów pipeline'u, a następnie przyjrzymy się konfiguracji w ścisłym sensie.

---

## 1. Parametry pipeline'u

Dla wszystkich pipeline'ów nf-core możesz uzyskać pełną listę parametrów bezpośrednio z wiersza poleceń, używając flagi `--help`, która sama w sobie jest parametrem pipeline'u.

### 1.1. Pobieranie listy parametrów za pomocą `--help`

Uruchom polecenie pomocy dla pipeline'u demo:

```bash
nextflow run nf-core/demo --help
```

??? success "Wyjście polecenia"

    ```console
     N E X T F L O W   ~  version 26.04.4

    Launching `https://github.com/nf-core/demo` [adoring_meucci] revision: 32893afef8 [master]


    ------------------------------------------------------
                                            ,--./,-.
            ___     __   __   __   ___     /,-._.--~'
      |\ | |__  __ /  ` /  \ |__) |__         }  {
      | \| |       \__, \__/ |  \ |___     \`-._,-`-,
                                            `._,._,'
      nf-core/demo 1.2.0
    ------------------------------------------------------
    Typical pipeline command:

      nextflow run nf-core/demo -profile <docker/singularity/.../institute> --input samplesheet.csv --outdir <OUTDIR>


    Input/output options
      --input                       [string] Path to a metadata file containing information about the samples in the experiment.
      --outdir                      [string] The output directory where the results will be saved. You have to use absolute paths to storage on Cloud infrastructure.
      --email                       [string] Email address for completion summary.
      --multiqc_title               [string] MultiQC report title. Printed as page header, used for filename if not otherwise specified.

    Reference genome options
      --genome                      [string] Name of iGenomes reference.
      --fasta                       [string] Path to FASTA genome file.

    Process skipping options
      --skip_trim                   [boolean] Skip trimming fastq files with seqtk

    Generic options
      --multiqc_methods_description [string]          Custom MultiQC yaml file containing HTML including a methods description.
      --help                        [boolean, string] Display the help message.
      --help_full                   [boolean]         Display the full detailed help message.
      --show_hidden                 [boolean]         Display hidden parameters in the help message (only works when --help or --help_full are provided).
     !! Hiding 19 param(s), use the `--showHidden` parameter to show them !!
    ------------------------------------------------------

    * The pipeline
        https://doi.org/10.5281/zenodo.12192442

    * The nf-core framework
        https://doi.org/10.1038/s41587-020-0439-x

    * Software dependencies
        https://github.com/nf-core/demo/blob/master/CITATIONS.md
    ```

Jak widać, wynik grupuje parametry w kategorie (opcje wejścia/wyjścia, opcje genomu referencyjnego itd.) wraz z typami i opisami każdego z nich.

Kategoryzacja ta jest określana przez plik schematu, który omówimy poniżej.
W zwykłych pipeline'ach Nextflow `--help` działa tylko wtedy, gdy deweloper zaimplementował tę funkcję ręcznie.

!!! tip "Wskazówka"

    Użyj `--help --show_hidden`, aby zobaczyć dodatkowe parametry ukryte domyślnie, takie jak `--publish_dir_mode` czy `--monochrome_logs`.

### 1.2. Ustawianie wartości parametrów

Jak omówiono w [Hello Config](../hello_nextflow/06_hello_config.md), wartości parametrów możesz ustawiać w wierszu poleceń za pomocą `--nazwa_parametru` lub zebrać zestaw parametrów w pliku YAML i przekazać go za pomocą `-params-file`.
Oba podejścia działają tak samo w pipeline'ach nf-core.

Na przykład, aby pominąć krok przycinania, chcemy ustawić parametr boolean `skip_trim` na `true`.
W Twoim katalogu roboczym znajduje się plik parametrów `my_params.yml` z już ustawioną tą wartością:

```yaml title="my_params.yml"
skip_trim: true
```

Przekaż go za pomocą `-params-file`:

```bash
nextflow run nf-core/demo -profile test,docker --outdir demo-results-notrim -params-file my_params.yml
```

??? success "Wyjście polecenia"

    ```console
     N E X T F L O W   ~  version 26.04.4

    Launching `https://github.com/nf-core/demo` [focused_heisenberg] revision: 32893afef8 [master]


    ------------------------------------------------------
                                            ,--./,-.
            ___     __   __   __   ___     /,-._.--~'
      |\ | |__  __ /  ` /  \ |__) |__         }  {
      | \| |       \__, \__/ |  \ |___     \`-._,-`-,
                                            `._,._,'
      nf-core/demo 1.2.0
    ------------------------------------------------------

    Input/output options
      input                     : https://raw.githubusercontent.com/nf-core/test-datasets/viralrecon/samplesheet/samplesheet_test_illumina_amplicon.csv
      outdir                    : demo-results-notrim

    Process skipping options
      skip_trim                 : true

    Institutional config options
      config_profile_name       : Test profile
      config_profile_description: Minimal test dataset to check pipeline function

    Generic options
      trace_report_suffix       : 2026-07-03_22-08-47

    Core Nextflow options
      revision                  : master
      runName                   : focused_heisenberg
      containerEngine           : docker
      launchDir                 : /workspaces/training/nfcore-use
      workDir                   : /workspaces/training/nfcore-use/work
      projectDir                : /workspaces/.nextflow/assets/.repos/nf-core/demo/clones/32893afef8076a03a2767a020b3f0cab2e0b40b2
      userName                  : root
      profile                   : test,docker
      configFiles               : /workspaces/.nextflow/assets/.repos/nf-core/demo/clones/32893afef8076a03a2767a020b3f0cab2e0b40b2/nextflow.config

    !! Only displaying parameters that differ from the pipeline defaults !!
    ------------------------------------------------------

    * The pipeline
        https://doi.org/10.5281/zenodo.12192442

    * The nf-core framework
        https://doi.org/10.1038/s41587-020-0439-x

    * Software dependencies
        https://github.com/nf-core/demo/blob/master/CITATIONS.md

    executor >  local (5)
    [7a/f3599e] NFCORE_DEMO:DEMO:FASTQC (SAMPLE3_SE) | 3 of 3 ✔
    [b0/2f0bdc] NFCORE_DEMO:DEMO:COWPY               | 1 of 1 ✔
    [c3/3c2278] NFCORE_DEMO:DEMO:MULTIQC (demo)      | 1 of 1 ✔
    -[nf-core/demo] Pipeline completed successfully-
    ```

Proces `SEQTK_TRIM` nie pojawia się już w wynikach.

!!! warning "Ostrzeżenie"

    **Ustawianie parametrów boolean w wierszu poleceń**

    Począwszy od Nextflow w wersji 26.04, wszystkie wartości podawane w wierszu poleceń są traktowane jako ciągi znaków (string).
    W przypadku parametru boolean takiego jak `skip_trim`, przekazanie go jako samodzielnej flagi (`--skip_trim`) lub jako `--skip_trim true` jest interpretowane jako **string** `"true"`, co powoduje błąd walidacji schematu:

    ```console
    * --skip_trim (true): Value is [string] but should be [boolean]
    ```

    Aby ustawić parametr boolean na prawdziwą wartość `true`/`false`, użyj `-params-file` jak pokazano powyżej lub ustaw go w pliku konfiguracyjnym.
    Parametry typu string, integer i ścieżki do plików nie są tym dotknięte i nadal można je ustawiać bezpośrednio w wierszu poleceń.
    W tym kursie ten wzorzec jest stosowany konsekwentnie dla parametrów boolean.

    **Używanie własnych plików konfiguracyjnych**

    Choć technicznie możliwe jest ustawianie parametrów pipeline'u w niestandardowym pliku konfiguracyjnym przekazywanym za pomocą `-c`, może to nie nadpisywać wartości domyślnych już ustawionych w pliku `nextflow.config` pipeline'u, w zależności od reguł pierwszeństwa konfiguracji Nextflow.
    Użycie `--nazwa_parametru` w wierszu poleceń lub `-params-file` jest bardziej niezawodne, ponieważ te metody zawsze mają pierwszeństwo.

    Jako ogólna zasada: jeśli parametr pojawia się w wynikach `--help`, ustaw go przez wiersz poleceń lub plik parametrów, a nie przez plik konfiguracyjny.

### 1.3. Walidacja parametrów

Ciekawostka: polecenie `--help` działa dla wszystkich pipeline'ów nf-core, ponieważ projekt nf-core wymaga od deweloperów formalnego zdefiniowania wszystkich parametrów pipeline'u w pliku schematu JSON (`nextflow_schema.json`).
Schemat ten rejestruje typ, opis, wartość domyślną i grupowanie każdego parametru.

Poza obsługą wyników `--help`, plik schematu umożliwia również automatyczną walidację podczas uruchamiania.
Oznacza to, że Nextflow może sprawdzić, czy każdy przekazany parametr istnieje i ma odpowiednią wartość (właściwego typu, w dozwolonym zakresie wartości itd.).

Omawiamy to szczegółowo w [sekcji dotyczącej walidacji danych wejściowych](../nfcore_build/04_input_validation.md), ale możesz już zobaczyć to w działaniu, podając pipeline'owi demo nieprawidłowe dane wejściowe.

#### 1.3.1. Nierozpoznane parametry

Spróbuj przekazać parametr, który nie istnieje:

```bash
nextflow run nf-core/demo -profile test,docker --outdir demo-results --foobar "invalid"
```

Wynik w konsoli zawiera ostrzeżenie:

```console
WARN: The following invalid input values have been detected:

* --foobar: invalid
```

Pipeline nadal działa, ale ostrzeżenie natychmiast informuje Cię, że `--foobar` nie jest rozpoznanym parametrem.
Ma to na celu zwrócenie Twojej uwagi na niekrytyczne literówki, takie jak użycie `--outDir` zamiast `--outdir`, co może pomóc uniknąć marnowania czasu i zasobów obliczeniowych.

#### 1.3.2. Nieprawidłowe wartości parametrów

Walidacja sprawdza również **wartości** parametrów.
Parametr `--skip_trim` jest flagą boolean, więc przekazanie wartości tekstowej powoduje natychmiastowe zatrzymanie pipeline'u:

```bash
nextflow run nf-core/demo -profile test,docker --outdir demo-results --skip_trim yes
```

```console
ERROR ~ Validation of pipeline parameters failed!

 -- Check '.nextflow.log' file for details
The following invalid input values have been detected:

* --skip_trim (yes): Value is [string] but should be [boolean]

 -- Check script '/workspaces/.nextflow/assets/.repos/nf-core/demo/clones/32893afef8076a03a2767a020b3f0cab2e0b40b2/subworkflows/nf-core/utils_nfschema_plugin/main.nf' at line: 68 or see '.nextflow.log' file for more details
```

Pipeline zatrzymuje się przed uruchomieniem jakichkolwiek procesów, oszczędzając Ci nieudanego lub nieprawidłowego wykonania.
Jak wspomniano w [1.2](#12-set-parameter-values), parametry boolean powinny być ustawiane na prawdziwą wartość `true`/`false` w pliku parametrów, a nie przekazywane w wierszu poleceń, ponieważ wartości z wiersza poleceń są traktowane jako ciągi znaków.

### 1.4. Walidacja danych wejściowych

Ta sama logika walidacji może być również używana do sprawdzania poprawności plików wejściowych.
Na przykład, jeśli pipeline oczekuje samplesheet jako głównego wejścia danych (co dotyczy wielu, jeśli nie większości pipeline'ów nf-core), deweloper może dostarczyć schemat wejściowy (odrębny od schematu parametrów) opisujący, jak powinien być zorganizowany plik wejściowy.

Następnie, w czasie wykonywania, Nextflow może sprawdzić, czy dostarczony plik wejściowy jest prawidłowy.

Omawiamy to również szczegółowo w [sekcji dotyczącej walidacji danych wejściowych](../nfcore_build/04_input_validation.md), ale możesz już zobaczyć to w działaniu, podając pipeline'owi demo nieprawidłowy samplesheet.

Pipeline `nf-core/demo` oczekuje pliku CSV z kolumnami `sample`, `fastq_1` i `fastq_2`.
Jest to zdefiniowane w pliku schematu (`assets/schema_input.json`), który określa oczekiwaną strukturę, typy kolumn i ograniczenia.

??? abstract "Plik schematu dla danych wejściowych"

    ```json title="assets/schema_input.json"
    {
        "$schema": "https://json-schema.org/draft/2020-12/schema",
        "$id": "https://raw.githubusercontent.com/nf-core/demo/master/assets/schema_input.json",
        "title": "nf-core/demo pipeline - params.input schema",
        "description": "Schema for the file provided with params.input",
        "type": "array",
        "items": {
            "type": "object",
            "properties": {
                "sample": {
                    "type": "string",
                    "pattern": "^\\S+$",
                    "errorMessage": "Sample name must be provided and cannot contain spaces",
                    "meta": ["id"]
                },
                "fastq_1": {
                    "type": "string",
                    "format": "file-path",
                    "exists": true,
                    "pattern": "^([\\S\\s]*\\/)?[^\\s\\/]+\\.f(ast)?q\\.gz$",
                    "errorMessage": "FastQ file for reads 1 must be provided, cannot contain spaces and must have extension '.fq.gz' or '.fastq.gz'"
                },
                "fastq_2": {
                    "type": "string",
                    "format": "file-path",
                    "exists": true,
                    "pattern": "^([\\S\\s]*\\/)?[^\\s\\/]+\\.f(ast)?q\\.gz$",
                    "errorMessage": "FastQ file for reads 2 cannot contain spaces and must have extension '.fq.gz' or '.fastq.gz'"
                }
            },
            "required": ["sample", "fastq_1"]
        }
    }
    ```

Schemat określa, że `sample` i `fastq_1` są wymagane, podczas gdy `fastq_2` jest opcjonalne (obsługując zarówno dane paired-end, jak i single-end).
Ścieżki do plików są walidowane pod kątem istnienia i wzorca rozszerzenia.

Aby to zademonstrować, w Twoim katalogu roboczym znajduje się nieprawidłowy samplesheet o nazwie `malformed_samplesheet.csv`:

```csv title="malformed_samplesheet.csv"
sample,fastq_2
SAMPLE1,/not/a/real/file.fastq.gz
```

W tym samplesheet brakuje wymaganej kolumny `fastq_1`, a ścieżka do pliku w `fastq_2` nie istnieje.

Uruchom pipeline demo, używając `malformed_samplesheet.csv` jako danych wejściowych:

```bash
nextflow run nf-core/demo -profile test,docker --outdir demo-results --input malformed_samplesheet.csv
```

```console
ERROR ~ Validation of pipeline parameters failed!

 -- Check '.nextflow.log' file for details
The following invalid input values have been detected:

* --input (malformed_samplesheet.csv): Validation of file failed:
    -> Entry 1: Error for field 'fastq_2' (/not/a/real/file.fastq.gz): the file or directory '/not/a/real/file.fastq.gz' does not exist (FastQ file for reads 2 cannot contain spaces and must have extension '.fq.gz' or '.fastq.gz')
    -> Entry 1: Missing required field(s): fastq_1

 -- Check script '/workspaces/.nextflow/assets/.repos/nf-core/demo/clones/32893afef8076a03a2767a020b3f0cab2e0b40b2/subworkflows/nf-core/utils_nfschema_plugin/main.nf' at line: 68 or see '.nextflow.log' file for more details
```

Jak widać, pipeline zatrzymuje się natychmiast i zgłasza **wszystkie** błędy walidacji naraz.
nf-schema nie zatrzymuje się na pierwszym błędzie — zbiera wszystkie problemy i wyświetla je razem, dzięki czemu możesz naprawić wszystko za jednym razem, zamiast odkrywać kolejne problemy jeden po drugim.

Każdy błąd wskazuje dokładny wpis i pole, które spowodowało problem, więc możesz poprawić samplesheet i ponownie uruchomić pipeline z pewnością, że nie zakończy się niepowodzeniem w późniejszym momencie, gdy Nextflow będzie próbował uzyskać dostęp do ścieżki pliku.

Dla deweloperów wszystko to jest omówione szczegółowo w [Części 4 kursu Build with nf-core](../nfcore_build/04_input_validation.md).

### Podsumowanie

Wiesz już, jak uzyskać pełną listę parametrów pipeline'u za pomocą `--help`, ustawiać je przez wiersz poleceń lub plik parametrów, oraz jak Nextflow waliduje zarówno wartości parametrów, jak i pliki wejściowe względem schematów pipeline'u.

### Co dalej?

Dowiedz się o drugim rodzaju konfiguracji: jak pipeline jest uruchamiany, obejmując alokację zasobów i argumenty narzędzi.

---

## 2. Konfiguracja

Konfiguracja w ścisłym sensie kontroluje **sposób** uruchamiania pipeline'u: alokację zasobów, argumenty specyficzne dla narzędzi, miejsce wykonywania zadań oraz używany system pakowania oprogramowania.

Pipeline'y nf-core zawierają domyślną konfigurację w `nextflow.config` i katalogu `conf/`.
Przed nadpisaniem czegokolwiek warto wiedzieć, gdzie znajdują się wartości domyślne.

### 2.1. Przeglądanie plików konfiguracyjnych

W [Części 1](./01_run_demo.md) widziałeś już, że kod źródłowy pipeline'u znajduje się w `$NXF_HOME/assets`.
Używając dowiązania symbolicznego `pipelines` utworzonego w [Części 1](./01_run_demo.md), wylistuj pliki konfiguracyjne, aby zobaczyć, co jest dostępne:

```bash
ls pipelines/nf-core/demo/conf/
```

```console
base.config
containers_conda_lock_files_amd64.config
containers_conda_lock_files_arm64.config
containers_docker_amd64.config
containers_docker_arm64.config
containers_singularity_https_amd64.config
containers_singularity_https_arm64.config
containers_singularity_oras_amd64.config
containers_singularity_oras_arm64.config
igenomes.config
igenomes_ignored.config
modules.config
test.config
test_full.config
```

<figure class="excalidraw">
--8<-- "docs/en/docs/nfcore_use/img/nfcore_config_files.excalidraw.svg"
</figure>

Najważniejsze pliki konfiguracyjne to:

- **`conf/base.config`**: Definiuje etykiety zasobów (`process_low`, `process_medium`, `process_high`), które przypisują procesory, pamięć i czas do procesów. Gdy widzisz, że proces używa więcej zasobów niż oczekiwano, to właśnie stąd pochodzą te wartości domyślne.
- **`conf/modules.config`**: Ustawia argumenty narzędzi dla poszczególnych procesów (`ext.args`) oraz ustawienia publikowania wyników (`publishDir`). Otwórz ten plik, aby zobaczyć, jakie argumenty każde narzędzie otrzymuje domyślnie.
- **`conf/test.config`**: Profil testowy użyty w [Części 1](./01_run_demo.md), który ogranicza zasoby za pomocą `resourceLimits` i ustawia testowy samplesheet. Aktywowany za pomocą `-profile test`.
  Istnieje również `conf/test_full.config` do uruchamiania z pełnowymiarowym zestawem danych testowych, przydatny do benchmarkingu.

Centralny plik `nextflow.config` ładuje wszystkie powyższe i ustawia odpowiednie wartości domyślne dla wszystkiego.

Jeśli chcesz zmodyfikować którekolwiek z ustawień określonych w tych plikach, nie modyfikuj żadnego z nich bezpośrednio.
Zamiast tego utwórz własny plik konfiguracyjny i przekaż go za pomocą `-c`.
Podane przez Ciebie wartości nadpiszą wartości domyślne ustawione w tamtych plikach.

Wypróbujmy to w praktyce.

### 2.2. Dostosowywanie zasobów procesów i argumentów narzędzi

Moduły nf-core obsługują dwa typowe rodzaje nadpisywania konfiguracji: **alokację zasobów** (procesory, pamięć, czas) oraz **argumenty narzędzi** za pomocą `ext.args`.

Wiele narzędzi wiersza poleceń ma argumenty, które nie są wystarczająco powszechnie używane, aby być udostępnione jako parametry pipeline'u.
Konwencja `ext.args` pozwala przekazywać te argumenty do bazowego narzędzia przez plik konfiguracyjny.

Plik `custom.config` dostarczony w Twoim katalogu roboczym demonstruje oba rodzaje nadpisywania:

```groovy title="custom.config" linenums="1"
process {
    withName: 'FASTQC' {
        cpus = 2
        memory = 4.GB
    }
    withName: 'SEQTK_TRIM' {
        ext.args = '-b 5'
    }
}
```

Pierwszy blok nadpisuje alokację zasobów dla `FASTQC`.
Domyślnie `FASTQC` używa etykiety `process_medium` z `base.config`, która przydziela 6 procesorów i 36 GB pamięci; tutaj ograniczamy to do 2 procesorów i 4 GB.

Drugi blok przekazuje dodatkowy argument do `SEQTK_TRIM` za pomocą `ext.args`.
Flaga `-b 5` mówi `seqtk trimfq`, aby przyciął 5 zasad z początku każdego odczytu, oprócz przycinania jakościowego.

Uruchom pipeline z tą konfiguracją:

```bash
nextflow run nf-core/demo -profile test,docker --outdir demo-results-custom -c custom.config
```

??? success "Wyjście polecenia"

    ```console
    executor >  local (8)
    [95/b32876] NFCORE_DEMO:DEMO:FASTQC (SAMPLE1_PE)     | 3 of 3 ✔
    [17/428668] NFCORE_DEMO:DEMO:SEQTK_TRIM (SAMPLE1_PE) | 3 of 3 ✔
    [cf/85991a] NFCORE_DEMO:DEMO:COWPY                   | 1 of 1 ✔
    [3c/94a7a0] NFCORE_DEMO:DEMO:MULTIQC (demo)          | 1 of 1 ✔
    -[nf-core/demo] Pipeline completed successfully-
    ```

Flaga `-c` dodaje Twoją konfigurację na wierzch wbudowanej konfiguracji pipeline'u.

Aby sprawdzić, czy nadpisanie `ext.args` zadziałało, znajdź hash katalogu roboczego procesu `SEQTK_TRIM` z wyników uruchomienia (np. `work/17/428668...`) i sprawdź plik `.command.sh` wewnątrz niego:

```bash
cat work/17/428668/.command.sh
```

??? success "Wyjście polecenia"

    ```console
    #!/usr/bin/env bash -e -u -o pipefail
    printf "%s\n" sample1_R1.fastq.gz sample1_R2.fastq.gz | while read f;
    do
        seqtk \
            trimfq \
            -b 5 \
            $f \
            | gzip --no-name > SAMPLE1_PE_$(basename $f)
    done
    ...
    ```

Powinieneś zobaczyć `-b 5` w poleceniu `seqtk trimfq`.

Ważna rzecz do zapamiętania dotycząca `ext.args`: jeśli moduł ma już ustawioną wartość domyślną, Twoja wartość **całkowicie ją zastąpi**, a nie dołączy do niej.
Na przykład `FASTQC` ma domyślnie ustawione `ext.args = '--quiet'` w `conf/modules.config`:

```groovy title="conf/modules.config" linenums="21" hl_lines="2"
    withName: FASTQC {
        ext.args   = '--quiet'
        publishDir = [
            path: { "${params.outdir}/fastqc/${meta.id}" },
            mode: params.publish_dir_mode,
            pattern: "*.{html,json}",
        ]
    }
```

Jeśli ustawisz `ext.args = '--kmers 8'` dla `FASTQC`, flaga `--quiet` nie będzie już stosowana.
Aby zachować obie, ustaw `ext.args = '--quiet --kmers 8'`.

Zawsze sprawdzaj domyślną konfigurację modułu przed nadpisaniem `ext.args`.

### Podsumowanie

Wiesz już, gdzie znajdują się domyślne ustawienia konfiguracji pipeline'ów nf-core oraz jak nadpisywać alokacje zasobów i argumenty narzędzi za pomocą niestandardowego pliku konfiguracyjnego.

### Co dalej?

Przejdź do [Części 3](./03_run_production_pipeline.md), gdzie zastosujesz zdobytą wiedzę w prawdziwym pipeline'ie produkcyjnym.

---

## Podsumowanie

W tej części nauczyłeś się:

- Uzyskiwać pomoc, ustawiać parametry oraz rozumieć walidację parametrów i danych wejściowych
- Dostosowywać alokację zasobów i argumenty narzędzi za pomocą plików konfiguracyjnych
