# Część 3: Uruchamianie produkcyjnego pipeline'u

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Tłumaczenie wspomagane przez AI - [dowiedz się więcej i zasugeruj ulepszenia](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

W [Części 2](./02_configure_execution.md) nauczyłeś się ustawiać parametry i dostosowywać konfigurację dla nf-core/demo.
Teraz zastosujemy zdobytą wiedzę do prawdziwego produkcyjnego pipeline'u — nf-core/rnaseq.

---

## 1. Pobieranie i uruchamianie nf-core/rnaseq

Do tej pory korzystaliśmy z `nf-core/demo` — minimalnego pipeline'u zaprojektowanego na potrzeby szkolenia.
Teraz pobierzemy prawdziwy produkcyjny pipeline i uruchomimy go z profilem testowym.

Pipeline `nf-core/rnaseq` wykonuje podstawowe kroki analizy masowego sekwencjonowania RNA: kontrolę jakości, przycinanie adapterów, dopasowanie odczytów oraz kwantyfikację na poziomie genów.
Jest to prawdopodobnie najszerzej stosowany pipeline nf-core.

### 1.1. Pobieranie pipeline'u

Uruchom poniższe polecenie, aby go pobrać.

```bash
nextflow pull nf-core/rnaseq
```

??? success "Wyjście polecenia"

    ```console
    Checking nf-core/rnaseq ...
     downloaded from https://github.com/nf-core/rnaseq.git - revision: 1f03b53ef7 [master]
    ```

Pipeline jest teraz zapisany w lokalnej pamięci podręcznej i gotowy do uruchomienia.

### 1.2. Uruchamianie profilu testowego

Uruchom go z profilem testowym i Docker:

```bash
nextflow run nf-core/rnaseq -profile test,docker --outdir rnaseq-results
```

??? failure "Wyjście polecenia"

    ```console
     N E X T F L O W   ~  version 26.04.4

    Launching `https://github.com/nf-core/rnaseq` [suspicious_dijkstra] revision: 1f03b53ef7 [master]

    ------------------------------------------------------
                                            ,--./,-.
            ___     __   __   __   ___     /,-._.--~'
      |\ | |__  __ /  ` /  \ |__) |__         }  {
      | \| |       \__, \__/ |  \ |___     \`-._,-`-,
                                            `._,._,'
      nf-core/rnaseq 3.26.0
    ------------------------------------------------------
    ...
    [-        ] NFCORE_RNASEQ:RNASEQ:FASTQ_QC_TRIM_FILTER_SETSTRANDEDNESS:FQ_LINT  -
    Plus 47 more processes waiting for tasks…

    Execution cancelled -- Finishing pending tasks before exit
    -[nf-core/rnaseq] Pipeline completed with errors-
    ERROR ~ Error executing process > 'NFCORE_RNASEQ:RNASEQ:FASTQ_QC_TRIM_FILTER_SETSTRANDEDNESS:FQ_LINT (WT_REP2)'

    Caused by:
      Process requirement exceeds available memory -- req: 12 GB; avail: 7.7 GB


    Command executed:

      fq lint \
          --disable-validator P001 \
          SRR6357072_1.fastq.gz SRR6357072_2.fastq.gz > WT_REP2.fq_lint.txt

    Command exit status:
      -

    Command output:
      (empty)

    Work dir:
      /workspaces/training/nfcore-use/work/xx/xxxxxxxxxxxxxxxxxxxxxx

    Container:
      quay.io/biocontainers/fq:0.12.0--h9ee0642_0

    Tip: you can replicate the issue by changing to the process work dir and entering the command `bash .command.run`

     -- Check '.nextflow.log' file for details
    ERROR ~ Pipeline failed. Please refer to troubleshooting docs for common issues: https://nf-co.re/docs/running/troubleshooting
    ```

Kluczowy wiersz tego błędu to:

```console
Process requirement exceeds available memory -- req: 12 GB; avail: 7.7 GB
```

Domyślna maszyna Codespaces ma 8 GB RAM — tyle samo, ile wynosi typowy limit Docker Desktop.
Pipeline żąda 12 GB dla procesu `FQ_LINT`, co przekracza możliwości maszyny.

Te 12 GB pochodzi z etykiety zasobów `process_low` zdefiniowanej w `conf/base.config`:

```groovy title="conf/base.config"
withLabel:process_low {
    cpus   = { 2     * task.attempt }
    memory = { 12.GB * task.attempt }
    time   = { 4.h   * task.attempt }
}
```

Jednym z rozwiązań byłoby użycie maszyny o większej ilości zasobów, ale na potrzeby testów chcemy móc uruchamiać pipeline na dowolnym dostępnym sprzęcie.
Lepszym podejściem jest nadpisanie domyślnych zasobów w niestandardowym pliku konfiguracyjnym.

### 1.3. Ponowne uruchomienie z niestandardową konfiguracją

Udostępniamy Ci niestandardowy plik konfiguracyjny, który nadpisuje domyślne zasoby przypisane do etykiet.

??? full-code "laptop.config"

    ```groovy title="laptop.config"
    process {
        withLabel: 'process_low' {
            cpus   = 2
            memory = 6.GB
        }
        withLabel: 'process_medium' {
            cpus   = 4
            memory = 6.GB
        }
        withLabel: 'process_high' {
            cpus   = 6
            memory = 6.GB
        }
        withLabel: 'process_high_memory' {
            memory = 6.GB
        }
    }
    ```

W [Części 2](./02_configure_execution.md) poznałeś `withName:`, które pozwala wskazać pojedynczy proces po nazwie.
Tutaj używamy `withLabel:`, aby jednocześnie objąć wszystkie procesy współdzielące daną etykietę.

Plik ten jest już obecny w Twoim katalogu roboczym.
Przekaż go za pomocą `-c`, aby zastosować nadpisania:

```bash
nextflow run nf-core/rnaseq -profile test,docker -c laptop.config --outdir rnaseq-results
```

??? success "Wynik polecenia (uruchamianie pipeline'u)"

    ```console
     N E X T F L O W   ~  version 26.04.4

    Launching `https://github.com/nf-core/rnaseq` [romantic_faraday] revision: 1f03b53ef7 [master]

    ------------------------------------------------------
                                            ,--./,-.
            ___     __   __   __   ___     /,-._.--~'
      |\ | |__  __ /  ` /  \ |__) |__         }  {
      | \| |       \__, \__/ |  \ |___     \`-._,-`-,
                                            `._,._,'
      nf-core/rnaseq 3.26.0
    ------------------------------------------------------
    ...
    executor >  local (7)
    [xx/xxxxxx] NFCORE_RNASEQ:RNASEQ:FQ_LINT (RAP1_IAA_30M_REP1)   | 3 of 5
    [xx/xxxxxx] NFCORE_RNASEQ:RNASEQ:FASTQC (RAP1_IAA_30M_REP1)    | 2 of 5
    ...
    ```

Pipeline działa, a Ty możesz obserwować kolejno kończące się zadania.
Na tym minimalnym zbiorze testowym ukończy się w ciągu 15–20 minut, wykonując łącznie ponad 200 zadań.

Prawdziwe eksperymenty RNA-seq obejmują zazwyczaj dziesiątki próbek i trwają godziny lub dni.
Nextflow obsługuje harmonogramery HPC (SLURM, PBS, LSF) oraz platformy chmurowe (AWS, Google Cloud, Azure), które mogą znacznie skrócić czas wykonania dzięki rozdzieleniu pracy na wiele węzłów.
Konfiguracja tych środowisk wiąże się jednak ze znaczną złożonością.

Platforma Seqera (stworzona przez twórców Nextflow) udostępnia interfejs webowy do uruchamiania pipeline'ów Nextflow na infrastrukturze HPC lub chmurowej — własnej lub zarządzanej — z możliwościami zarządzania obliczeniami i danymi, które upraszczają uruchamianie pipeline'ów na dużą skalę.

!!! tip "Wskazówka"

    Naukowcy akademiccy mogą korzystać z Seqera Platform bezpłatnie w ramach [programu akademickiego Seqera](https://seqera.io/academic-program/).

### Podsumowanie

Pobrałeś `nf-core/rnaseq`, zobaczyłeś, jak działają etykiety zasobów nf-core, i nauczyłeś się je nadpisywać za pomocą niestandardowego pliku konfiguracyjnego.
Co ważniejsze, przekonałeś się, że lokalne uruchamianie to punkt wyjścia, a nie cel dla analiz w prawdziwej skali.

### Co dalej?

Poznałeś podstawy uruchamiania pipeline'ów nf-core.
Zajrzyj do sekcji [Następne kroki](next_steps.md), aby dowiedzieć się, co robić dalej.

---

## Podsumowanie

W tej części nauczyłeś się:

- Pobierać i uruchamiać produkcyjny pipeline (nf-core/rnaseq) oraz nadpisywać jego domyślne etykiety zasobów
