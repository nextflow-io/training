````groovy
        echo "Processing ${reads}" > ${meta.id}_report.txt
        echo "Sample: ${meta.id}" >> ${meta.id}_report.txt
        echo "Processed by: ${USER}" >> ${meta.id}_report.txt
        echo "Hostname: $(hostname)" >> ${meta.id}_report.txt
        echo "Date: $(date)" >> ${meta.id}_report.txt
        """
    ```

=== "Przed"

    ```groovy title="modules/generate_report.nf" linenums="10"
        script:
        """
        echo "Processing ${reads}" > ${meta.id}_report.txt
        echo "Sample: ${meta.id}" >> ${meta.id}_report.txt
        """
    ```

Jeśli uruchomisz to, zauważysz błąd - Nextflow próbuje zinterpretować `${USER}` jako zmienną Nextflow, która nie istnieje.

??? failure "Wyjście polecenia"

    ```console
    Error modules/generate_report.nf:15:27: `USER` is not defined
    │  15 |     echo "Processed by: ${USER}" >> ${meta.id}_report.txt
    ╰     |                           ^^^^

    ERROR ~ Script compilation failed
    ```

Musimy to escape'ować, aby Bash mógł to obsłużyć.

Napraw to, escape'ując zmienne powłoki i podstawienia poleceń ukośnikiem odwrotnym (`\`):

=== "Po"

    ```groovy title="modules/generate_report.nf" linenums="10" hl_lines="5-7"
        script:
        """
        echo "Processing ${reads}" > ${meta.id}_report.txt
        echo "Sample: ${meta.id}" >> ${meta.id}_report.txt
        echo "Processed by: \${USER}" >> ${meta.id}_report.txt
        echo "Hostname: \$(hostname)" >> ${meta.id}_report.txt
        echo "Date: \$(date)" >> ${meta.id}_report.txt
        """
    ```

=== "Przed"

    ```groovy title="modules/generate_report.nf" linenums="10"
        script:
        """
        echo "Processing ${reads}" > ${meta.id}_report.txt
        echo "Sample: ${meta.id}" >> ${meta.id}_report.txt
        echo "Processed by: ${USER}" >> ${meta.id}_report.txt
        echo "Hostname: $(hostname)" >> ${meta.id}_report.txt
        echo "Date: $(date)" >> ${meta.id}_report.txt
        """
    ```

Teraz działa! Ukośnik odwrotny (`\`) mówi Nextflow "nie interpretuj tego, przekaż to do Bash."

### Podsumowanie

W tej sekcji nauczyłeś się technik **przetwarzania ciągów znaków**:

- **Wyrażenia regularne do parsowania plików**: Używanie operatora `=~` i wzorców regex (`~/wzorzec/`) do wyodrębniania metadanych ze złożonych konwencji nazewnictwa plików
- **Dynamiczne generowanie skryptów**: Używanie logiki warunkowej (if/else, operatory ternarne) do generowania różnych ciągów skryptowych na podstawie charakterystyki wejścia
- **Interpolacja zmiennych**: Zrozumienie, kiedy Nextflow interpretuje ciągi znaków, a kiedy robi to powłoka
  - `${var}` - Zmienne Nextflow (interpolowane przez Nextflow w czasie kompilacji workflow)
  - `\${var}` - Zmienne środowiskowe powłoki (escape'owane, przekazywane do bash w czasie wykonania)
  - `\$(cmd)` - Podstawienie polecenia powłoki (escape'owane, wykonywane przez bash w czasie wykonania)

Te wzorce przetwarzania i generowania ciągów znaków są niezbędne do obsługi różnorodnych formatów plików i konwencji nazewnictwa, z którymi spotkasz się w rzeczywistych workflow bioinformatycznych.

---

## 3. Tworzenie Funkcji Wielokrotnego Użytku

Złożona logika workflow umieszczona inline w operatorach kanałów lub definicjach procesów zmniejsza czytelność i łatwość utrzymania. **Funkcje** pozwalają wyodrębnić tę logikę do nazwanych, wielokrotnie używalnych komponentów.

Nasza operacja map stała się długa i złożona. Wyodrębnijmy ją do funkcji wielokrotnego użytku używając słowa kluczowego `def`.

Aby zilustrować, jak to wygląda z naszym istniejącym workflow, wprowadź poniższą modyfikację, używając `def` do zdefiniowania funkcji wielokrotnego użytku o nazwie `separateMetadata`:

=== "Po"

    ```groovy title="main.nf" linenums="1" hl_lines="4-24 29"
    include { FASTP } from './modules/fastp.nf'
    include { GENERATE_REPORT } from './modules/generate_report.nf'

    def separateMetadata(row) {
        def sample_meta = [
            id: row.sample_id.toLowerCase(),
            organism: row.organism,
            tissue: row.tissue_type.replaceAll('_', ' ').toLowerCase(),
            depth: row.sequencing_depth.toInteger(),
            quality: row.quality_score.toDouble()
        ]
        def fastq_path = file(row.file_path)

        def m = (fastq_path.name =~ /^(.+)_S(\d+)_L(\d{3})_(R[12])_(\d{3})\.fastq(?:\.gz)?$/)
        def file_meta = m ? [
            sample_num: m[0][2].toInteger(),
            lane: m[0][3],
            read: m[0][4],
            chunk: m[0][5]
        ] : [:]

        def priority = sample_meta.quality > 40 ? 'high' : 'normal'
        return tuple(sample_meta + file_meta + [priority: priority], fastq_path)
    }

    workflow {
        ch_samples = channel.fromPath("./data/samples.csv")
            .splitCsv(header: true)
            .map{ row -> separateMetadata(row) }

        ch_fastp = FASTP(ch_samples)
        GENERATE_REPORT(ch_samples)
    }
    ```

=== "Przed"

    ```groovy title="main.nf" linenums="1" hl_lines="7-27"
    include { FASTP } from './modules/fastp.nf'
    include { GENERATE_REPORT } from './modules/generate_report.nf'

    workflow {
        ch_samples = channel.fromPath("./data/samples.csv")
            .splitCsv(header: true)
            .map { row ->
                def sample_meta = [
                    id: row.sample_id.toLowerCase(),
                    organism: row.organism,
                    tissue: row.tissue_type.replaceAll('_', ' ').toLowerCase(),
                    depth: row.sequencing_depth.toInteger(),
                    quality: row.quality_score.toDouble()
                ]
                def fastq_path = file(row.file_path)

                def m = (fastq_path.name =~ /^(.+)_S(\d+)_L(\d{3})_(R[12])_(\d{3})\.fastq(?:\.gz)?$/)
                def file_meta = m ? [
                    sample_num: m[0][2].toInteger(),
                    lane: m[0][3],
                    read: m[0][4],
                    chunk: m[0][5]
                ] : [:]

                def priority = sample_meta.quality > 40 ? 'high' : 'normal'
                return tuple(sample_meta + file_meta + [priority: priority], fastq_path)
            }

        ch_fastp = FASTP(ch_samples)
        GENERATE_REPORT(ch_samples)
    }
    ```

Wyodrębniając tę logikę do funkcji, zredukowaliśmy rzeczywistą logikę workflow do czegoś znacznie czystszego:

```groovy title="minimalny workflow"
    ch_samples = channel.fromPath("./data/samples.csv")
        .splitCsv(header: true)
        .map{ row -> separateMetadata(row) }

    ch_fastp = FASTP(ch_samples)
    GENERATE_REPORT(ch_samples)
````

To sprawia, że logika workflow jest znacznie łatwiejsza do przeczytania i zrozumienia na pierwszy rzut oka. Funkcja `separateMetadata` enkapsuluje całą złożoną logikę parsowania i wzbogacania metadanych, czyniąc ją wielokrotnie używalną i testowalną.

Uruchom workflow, aby upewnić się, że wciąż działa:

```bash
nextflow run main.nf
```

??? success "Wyjście polecenia"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `main.nf` [admiring_panini] DSL2 - revision: 8cc832e32f

    executor >  local (6)
    [8c/2e3f91] process > FASTP (3)           [100%] 3 of 3 ✔
    [7a/1b4c92] process > GENERATE_REPORT (3) [100%] 3 of 3 ✔
    ```

Wyjście powinno pokazywać, że oba procesy zakończyły się pomyślnie. Workflow jest teraz znacznie czystszy i łatwiejszy w utrzymaniu, z całą złożoną logiką przetwarzania metadanych enkapsulowaną w funkcji `separateMetadata`.

### Podsumowanie

W tej sekcji nauczyłeś się **tworzenia funkcji**:

- **Definiowanie funkcji za pomocą `def`**: Słowo kluczowe do tworzenia nazwanych funkcji (jak `def` w Python lub `function` w JavaScript)
- **Zakres funkcji**: Funkcje zdefiniowane na poziomie skryptu są dostępne w całym workflow Nextflow
- **Wartości zwracane**: Funkcje automatycznie zwracają ostatnie wyrażenie lub używają jawnego `return`
- **Czystszy kod**: Wyodrębnianie złożonej logiki do funkcji to fundamentalna praktyka inżynierii oprogramowania w każdym języku

Następnie zbadamy, jak używać closures w dyrektywach procesów do dynamicznej alokacji zasobów.

---

## 4. Dynamiczne Dyrektywy Zasobów z Closures

Do tej pory używaliśmy skryptowania w bloku `script` procesów. Ale **closures** (przedstawione w Sekcji 1.1) są również niezwykle użyteczne w dyrektywach procesów, szczególnie do dynamicznej alokacji zasobów. Dodajmy dyrektywy zasobów do naszego procesu FASTP, które dostosowują się w oparciu o charakterystykę próbki.

### 4.1. Alokacja zasobów specyficzna dla próbki

Obecnie nasz proces FASTP używa domyślnych zasobów. Uczyńmy go inteligentniejszym, alokując więcej CPU dla próbek o wysokiej głębokości. Edytuj `modules/fastp.nf`, aby uwzględnić dynamiczną dyrektywę `cpus` i statyczną dyrektywę `memory`:

=== "Po"

    ```groovy title="modules/fastp.nf" linenums="1" hl_lines="4-5"
    process FASTP {
        container 'community.wave.seqera.io/library/fastp:0.24.0--62c97b06e8447690'

        cpus { meta.depth > 40000000 ? 2 : 1 }
        memory 2.GB

        input:
        tuple val(meta), path(reads)
    ```

=== "Przed"

    ```groovy title="modules/fastp.nf" linenums="1"
    process FASTP {
        container 'community.wave.seqera.io/library/fastp:0.24.0--62c97b06e8447690'

        input:
        tuple val(meta), path(reads)
    ```

Closure `{ meta.depth > 40000000 ? 2 : 1 }` używa **operatora ternarnego** (omówionego w Sekcji 1.1) i jest ewaluowane dla każdego zadania, umożliwiając alokację zasobów per próbka. Próbki o wysokiej głębokości (>40M odczytów) otrzymują 2 CPU, podczas gdy inne otrzymują 1 CPU.

!!! note "Dostęp do Zmiennych Wejściowych w Dyrektywach"

    Closure może uzyskać dostęp do dowolnych zmiennych wejściowych (jak tutaj `meta`), ponieważ Nextflow ewaluuje te closures w kontekście każdego wykonania zadania.

Uruchom workflow ponownie z opcją `-ansi-log false`, aby łatwiej zobaczyć hasze zadań.

```bash
nextflow run main.nf -ansi-log false
```

??? success "Wyjście polecenia"

    ```console
    N E X T F L O W  ~  version 25.10.2
    Launching `main.nf` [fervent_albattani] DSL2 - revision: fa8f249759
    [bd/ff3d41] Submitted process > FASTP (2)
    [a4/a3aab2] Submitted process > FASTP (1)
    [48/6db0c9] Submitted process > FASTP (3)
    [ec/83439d] Submitted process > GENERATE_REPORT (3)
    [bd/15d7cc] Submitted process > GENERATE_REPORT (2)
    [42/699357] Submitted process > GENERATE_REPORT (1)
    ```

Możesz sprawdzić dokładne polecenie `docker`, które zostało uruchomione, aby zobaczyć alokację CPU dla danego zadania:

```console title="Sprawdź polecenie docker"
cat work/48/6db0c9e9d8aa65e4bb4936cd3bd59e/.command.run | grep "docker run"
```

Powinieneś zobaczyć coś takiego:

```bash title="polecenie docker"
    docker run -i --cpu-shares 4096 --memory 2048m -e "NXF_TASK_WORKDIR" -v /workspaces/training/side-quests/essential_scripting_patterns:/workspaces/training/side-quests/essential_scripting_patterns -w "$NXF_TASK_WORKDIR" --name $NXF_BOXID community.wave.seqera.io/library/fastp:0.24.0--62c97b06e8447690 /bin/bash -ue /workspaces/training/side-quests/essential_scripting_patterns/work/48/6db0c9e9d8aa65e4bb4936cd3bd59e/.command.sh
```

W tym przykładzie wybraliśmy przykład, który zażądał 2 CPU (`--cpu-shares 2048`), ponieważ była to próbka o wysokiej głębokości, ale powinieneś zobaczyć różne alokacje CPU w zależności od głębokości próbki. Spróbuj tego również dla innych zadań.

### 4.2. Strategie ponownych prób

Innym potężnym wzorcem jest użycie `task.attempt` do strategii ponownych prób. Aby pokazać, dlaczego jest to użyteczne, zaczniemy od zmniejszenia alokacji pamięci dla FASTP do mniej niż potrzebuje. Zmień dyrektywę `memory` w `modules/fastp.nf` na `1.GB`:

=== "Po"

    ```groovy title="modules/fastp.nf" linenums="1" hl_lines="5"
    process FASTP {
        container 'community.wave.seqera.io/library/fastp:0.24.0--62c97b06e8447690'

        cpus { meta.depth > 40000000 ? 4 : 2 }
        memory 1.GB

        input:
        tuple val(meta), path(reads)
    ```

=== "Przed"

    ```groovy title="modules/fastp.nf" linenums="1" hl_lines="5"
    process FASTP {
        container 'community.wave.seqera.io/library/fastp:0.24.0--62c97b06e8447690'

        cpus { meta.depth > 40000000 ? 4 : 2 }
        memory 2.GB

        input:
        tuple val(meta), path(reads)
    ```

... i uruchom workflow ponownie:

```bash
nextflow run main.nf
```

??? failure "Wyjście polecenia"

    ```console hl_lines="2 11"
    Command exit status:
      137

    Command output:
      (empty)

    Command error:
      Detecting adapter sequence for read1...
      No adapter detected for read1

      .command.sh: line 7:   101 Killed                  fastp --in1 SAMPLE_002_S2_L001_R1_001.fastq --out1 sample_002_trimmed.fastq.gz --json sample_002.fastp.json --html sample_002.fastp.html --thread 2
    ```

To wskazuje, że proces został zabity za przekroczenie limitów pamięci.

To bardzo powszechny scenariusz w rzeczywistych workflow - czasami po prostu nie wiesz, ile pamięci będzie potrzebować zadanie, dopóki go nie uruchomisz.

Aby uczynić nasz workflow bardziej odpornym, możemy zaimplementować strategię ponownych prób, która zwiększa alokację pamięci przy każdej próbie, ponownie używając closure Groovy. Zmodyfikuj dyrektywę `memory`, aby mnożyć bazową pamięć przez `task.attempt`, i dodaj dyrektywy `errorStrategy 'retry'` oraz `maxRetries 2`:

=== "Po"

    ```groovy title="modules/fastp.nf" linenums="1" hl_lines="5-7"
    process FASTP {
        container 'community.wave.seqera.io/library/fastp:0.24.0--62c97b06e8447690'

        cpus { meta.depth > 40000000 ? 4 : 2 }
        memory { 1.GB * task.attempt }
        errorStrategy 'retry'
        maxRetries 2

        input:
        tuple val(meta), path(reads)
    ```

=== "Przed"

    ```groovy title="modules/fastp.nf" linenums="1" hl_lines="5"
    process FASTP {
        container 'community.wave.seqera.io/library/fastp:0.24.0--62c97b06e8447690'

        cpus { meta.depth > 40000000 ? 4 : 2 }
        memory 2.GB

        input:
        tuple val(meta), path(reads)
    ```

Teraz, jeśli proces zawiedzie z powodu niewystarczającej pamięci, Nextflow spróbuje ponownie z większą ilością pamięci:

- Pierwsza próba: 1 GB (task.attempt = 1)
- Druga próba: 2.GB (task.attempt = 2)

... i tak dalej, do limitu `maxRetries`.

### Podsumowanie

Dynamiczne dyrektywy z closures pozwalają Ci:

- Alokować zasoby na podstawie charakterystyki wejścia
- Implementować automatyczne strategie ponownych prób ze zwiększającymi się zasobami
- Łączyć wiele czynników (metadane, numer próby, priorytety)
- Używać logiki warunkowej do złożonych obliczeń zasobów

To sprawia, że Twoje workflow są zarówno bardziej wydajne (nie przealokują zasobów), jak i bardziej odporne (automatyczne ponowne próby z większymi zasobami).

---

## 5. Logika Warunkowa i Kontrola Procesów

Wcześniej używaliśmy `.map()` ze skryptowaniem do transformacji danych kanału. Teraz użyjemy logiki warunkowej do kontrolowania, które procesy wykonują się na podstawie danych — niezbędne dla elastycznych workflow dostosowujących się do różnych typów próbek.

[Operatory przepływu danych](https://www.nextflow.io/docs/latest/reference/operator.html) Nextflow przyjmują closures ewaluowane w czasie wykonania, umożliwiając logice warunkowej kierowanie decyzjami workflow na podstawie zawartości kanału.

### 5.1. Kierowanie za pomocą `.branch()`

Na przykład, wyobraźmy sobie, że nasze próbki sekwencjonowania muszą być przycinane za pomocą FASTP tylko jeśli są próbkami ludzkimi z pokryciem powyżej pewnego progu. Próbki mysie lub próbki o niskim pokryciu powinny być uruchamiane z Trimgalore (to jest wymyślony przykład, ale ilustruje punkt).

Dostarczyliśmy prosty proces Trimgalore w `modules/trimgalore.nf`, spójrz, jeśli chcesz, ale szczegóły nie są ważne dla tego ćwiczenia. Kluczowym punktem jest to, że chcemy kierować próbki na podstawie ich metadanych.

Dołącz nowy moduł z `modules/trimgalore.nf`:

=== "Po"

    ```groovy title="main.nf" linenums="1" hl_lines="2"
    include { FASTP } from './modules/fastp.nf'
    include { TRIMGALORE } from './modules/trimgalore.nf'
    ```

=== "Przed"

    ```groovy title="main.nf" linenums="1"
    include { FASTP } from './modules/fastp.nf'
    ```

... a następnie zmodyfikuj swój workflow `main.nf`, aby rozgałęziać próbki na podstawie ich metadanych i kierować je przez odpowiedni proces przycinania, tak:

=== "Po"

    ```groovy title="main.nf" linenums="28" hl_lines="5-12"
        ch_samples = channel.fromPath("./data/samples.csv")
            .splitCsv(header: true)
            .map { row -> separateMetadata(row) }

        trim_branches = ch_samples
            .branch { meta, reads ->
                fastp: meta.organism == 'human' && meta.depth >= 30000000
                trimgalore: true
            }

        ch_fastp = FASTP(trim_branches.fastp)
        ch_trimgalore = TRIMGALORE(trim_branches.trimgalore)
        GENERATE_REPORT(ch_samples)
    ```

=== "Przed"

    ```groovy title="main.nf" linenums="28" hl_lines="5"
        ch_samples = channel.fromPath("./data/samples.csv")
            .splitCsv(header: true)
            .map { row -> separateMetadata(row) }

        ch_fastp = FASTP(ch_samples)
        GENERATE_REPORT(ch_samples)
    ```

Uruchom ten zmodyfikowany workflow:

```bash
nextflow run main.nf
```

??? success "Wyjście polecenia"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `main.nf` [adoring_galileo] DSL2 - revision: c9e83aaef1

    executor >  local (6)
    [1d/0747ac] process > FASTP (2)           [100%] 2 of 2 ✔
    [cc/c44caf] process > TRIMGALORE (1)      [100%] 1 of 1 ✔
    [34/bd5a9f] process > GENERATE_REPORT (1) [100%] 3 of 3 ✔
    ```

Tutaj użyliśmy małych, ale potężnych wyrażeń warunkowych wewnątrz operatora `.branch{}` do kierowania próbek na podstawie ich metadanych. Próbki ludzkie z wysokim pokryciem przechodzą przez `FASTP`, podczas gdy wszystkie inne próbki przechodzą przez `TRIMGALORE`.

### 5.2. Używanie `.filter()` z Prawdziwością

Innym potężnym wzorcem do kontrolowania wykonywania workflow jest operator `.filter()`, który używa closure do określenia, które elementy powinny kontynuować przez pipeline. Wewnątrz closure filtra napiszesz **wyrażenia boolowskie**, które decydują, które elementy przechodzą.

Nextflow (jak wiele dynamicznych języków) ma koncepcję **"prawdziwości" (truthiness)**, która określa, jakie wartości są ewaluowane jako `true` lub `false` w kontekstach boolowskich:

- **Truthy (prawdziwe)**: Wartości nie-null, niepuste ciągi znaków, niezerowe liczby, niepuste kolekcje
- **Falsy (fałszywe)**: `null`, puste ciągi `""`, zero `0`, puste kolekcje `[]` lub `[:]`, `false`

To oznacza, że samo `meta.id` (bez jawnego `!= null`) sprawdza, czy ID istnieje i nie jest puste. Użyjmy tego do filtrowania próbek, które nie spełniają naszych wymagań jakościowych.

Dodaj następujący fragment przed operacją branch:

=== "Po"

    ```groovy title="main.nf" linenums="28" hl_lines="5-11"
        ch_samples = channel.fromPath("./data/samples.csv")
            .splitCsv(header: true)
            .map { row -> separateMetadata(row) }

        // Odfiltruj nieprawidłowe lub niskiej jakości próbki
        ch_valid_samples = ch_samples
            .filter { meta, reads ->
                meta.id && meta.organism && meta.depth >= 25000000
            }

        trim_branches = ch_valid_samples
            .branch { meta, reads ->
                fastp: meta.organism == 'human' && meta.depth >= 30000000
                trimgalore: true
            }
    ```

=== "Przed"

    ```groovy title="main.nf" linenums="28" hl_lines="5"
        ch_samples = channel.fromPath("./data/samples.csv")
            .splitCsv(header: true)
            .map { row -> separateMetadata(row) }

        trim_branches = ch_samples
            .branch { meta, reads ->
                fastp: meta.organism == 'human' && meta.depth >= 30000000
                trimgalore: true
            }
    ```

Uruchom workflow ponownie:

```bash
nextflow run main.nf
```

??? success "Wyjście polecenia"

    ```console
    N E X T F L O W  ~  version 25.10.2
    Launching `main.nf` [lonely_williams] DSL2 - revision: d0b3f121ec
    [94/b48eac] Submitted process > FASTP (2)
    [2c/d2b28f] Submitted process > GENERATE_REPORT (2)
    [65/2e3be4] Submitted process > GENERATE_REPORT (1)
    [94/b48eac] NOTE: Process `FASTP (2)` terminated with an error exit status (137) -- Execution is retried (1)
    [3e/0d8664] Submitted process > TRIMGALORE (1)
    [6a/9137b0] Submitted process > FASTP (1)
    [6a/9137b0] NOTE: Process `FASTP (1)` terminated with an error exit status (137) -- Execution is retried (1)
    [83/577ac0] Submitted process > GENERATE_REPORT (3)
    [a2/5117de] Re-submitted process > FASTP (1)
    [1f/a1a4ca] Re-submitted process > FASTP (2)
    ```

Ponieważ wybraliśmy filtr, który wyklucza niektóre próbki, wykonano mniej zadań.

Wyrażenie filtru `meta.id && meta.organism && meta.depth >= 25000000` łączy prawdziwość z jawnymi porównaniami:

- `meta.id && meta.organism` sprawdza, że oba pola istnieją i nie są puste (używając prawdziwości)
- `meta.depth >= 25000000` zapewnia wystarczającą głębokość sekwencjonowania za pomocą jawnego porównania

!!! note "Prawdziwość w Praktyce"

    Wyrażenie `meta.id && meta.organism` jest bardziej zwięzłe niż pisanie:
    ```groovy
    meta.id != null && meta.id != '' && meta.organism != null && meta.organism != ''
    ```

    To sprawia, że logika filtrowania jest znacznie czystsza i łatwiejsza do przeczytania.

### Podsumowanie

W tej sekcji nauczyłeś się używać logiki warunkowej do kontrolowania wykonywania workflow, używając interfejsów closure operatorów Nextflow takich jak `.branch{}` i `.filter{}`, wykorzystując prawdziwość do pisania zwięzłych wyrażeń warunkowych.

Nasz pipeline teraz inteligentnie kieruje próbki przez odpowiednie procesy, ale workflow produkcyjne muszą obsługiwać nieprawidłowe dane z gracją. Uczyńmy nasz workflow odpornym na brakujące lub null wartości.

---

## 6. Bezpieczna Nawigacja i Operatory Elvis

Nasza funkcja `separateMetadata` obecnie zakłada, że wszystkie pola CSV są obecne i prawidłowe. Ale co się dzieje z niekompletnymi danymi? Sprawdźmy.

### 6.1. Problem: Dostęp do Właściwości, Które Nie Istnieją

Powiedzmy, że chcemy dodać wsparcie dla opcjonalnych informacji o uruchomieniu sekwencjonowania. W niektórych laboratoriach próbki mogą mieć dodatkowe pole dla ID uruchomienia sekwencjonowania lub numeru partii, ale nasz obecny CSV nie ma tej kolumny. Spróbujmy uzyskać do niej dostęp mimo to.

Zmodyfikuj funkcję `separateMetadata`, aby uwzględnić pole run_id:

=== "Po"

    ```groovy title="main.nf" linenums="5" hl_lines="9"
    def separateMetadata(row) {
        def sample_meta = [
            id: row.sample_id.toLowerCase(),
            organism: row.organism,
            tissue: row.tissue_type.replaceAll('_', ' ').toLowerCase(),
            depth: row.sequencing_depth.toInteger(),
            quality: row.quality_score.toDouble()
        ]
        def run_id = row.run_id.toUpperCase()
    ```

=== "Przed"

    ```groovy title="main.nf" linenums="5"
    def separateMetadata(row) {
        def sample_meta = [
            id: row.sample_id.toLowerCase(),
            organism: row.organism,
            tissue: row.tissue_type.replaceAll('_', ' ').toLowerCase(),
            depth: row.sequencing_depth.toInteger(),
            quality: row.quality_score.toDouble()
        ]
    ```

Teraz uruchom workflow:

```bash
nextflow run main.nf
```

??? failure "Wyjście polecenia"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `main.nf` [trusting_torvalds] DSL2 - revision: b56fbfbce2

    ERROR ~ Cannot invoke method toUpperCase() on null object

    -- Check script 'main.nf' at line: 13 or see '.nextflow.log' file for more details
    ```

To powoduje awarię z NullPointerException.

Problem polega na tym, że `row.run_id` zwraca `null`, ponieważ kolumna `run_id` nie istnieje w naszym CSV. Kiedy próbujemy wywołać `.toUpperCase()` na `null`, następuje awaria. Tu operator bezpiecznej nawigacji ratuje sytuację.

### 6.2. Operator Bezpiecznej Nawigacji (`?.`)

Operator bezpiecznej nawigacji (`?.`) zwraca `null` zamiast rzucać wyjątek, gdy jest wywoływany na wartości `null`. Jeśli obiekt przed `?.` jest `null`, całe wyrażenie jest ewaluowane jako `null` bez wykonywania metody.

Zaktualizuj funkcję, aby używać bezpiecznej nawigacji:

=== "Po"

    ```groovy title="main.nf" linenums="4" hl_lines="9"
    def separateMetadata(row) {
        def sample_meta = [
            id: row.sample_id.toLowerCase(),
            organism: row.organism,
            tissue: row.tissue_type.replaceAll('_', ' ').toLowerCase(),
            depth: row.sequencing_depth.toInteger(),
            quality: row.quality_score.toDouble()
        ]
        def run_id = row.run_id?.toUpperCase()
    ```

=== "Przed"

    ```groovy title="main.nf" linenums="4" hl_lines="9"
    def separateMetadata(row) {
        def sample_meta = [
            id: row.sample_id.toLowerCase(),
            organism: row.organism,
            tissue: row.tissue_type.replaceAll('_', ' ').toLowerCase(),
            depth: row.sequencing_depth.toInteger(),
            quality: row.quality_score.toDouble()
        ]
        def run_id = row.run_id.toUpperCase()
    ```

Uruchom ponownie:

```bash
nextflow run main.nf
```

??? success "Wyjście polecenia"

    ```console
    <!-- TODO: wyjście -->
    ```

Brak awarii! Workflow teraz obsługuje brakujące pole z gracją. Kiedy `row.run_id` jest `null`, operator `?.` zapobiega wywołaniu `.toUpperCase()`, a `run_id` staje się `null` zamiast powodować wyjątek.

### 6.3. Operator Elvis (`?:`) dla Wartości Domyślnych

Operator Elvis (`?:`) dostarcza wartości domyślne, gdy lewa strona jest "falsy" (jak wyjaśniono wcześniej). Jest nazwany na cześć Elvisa Presleya, ponieważ `?:` wygląda jak jego słynne włosy i oczy, gdy oglądane z boku!

Teraz, gdy używamy bezpiecznej nawigacji, `run_id` będzie `null` dla próbek bez tego pola. Użyjmy operatora Elvis, aby dostarczyć wartość domyślną i dodajmy ją do naszej mapy `sample_meta`:

=== "Po"

    ```groovy title="main.nf" linenums="5" hl_lines="9-10"
    def separateMetadata(row) {
        def sample_meta = [
            id: row.sample_id.toLowerCase(),
            organism: row.organism,
            tissue: row.tissue_type.replaceAll('_', ' ').toLowerCase(),
            depth: row.sequencing_depth.toInteger(),
            quality: row.quality_score.toDouble()
        ]
        def run_id = row.run_id?.toUpperCase() ?: 'UNSPECIFIED'
        sample_meta.run = run_id
    ```

=== "Przed"

    ```groovy title="main.nf" linenums="5" hl_lines="9"
    def separateMetadata(row) {
        def sample_meta = [
            id: row.sample_id.toLowerCase(),
            organism: row.organism,
            tissue: row.tissue_type.replaceAll('_', ' ').toLowerCase(),
            depth: row.sequencing_depth.toInteger(),
            quality: row.quality_score.toDouble()
        ]
        def run_id = row.run_id?.toUpperCase()
    ```

Dodaj również operator `view()` w workflow, aby zobaczyć wyniki:

=== "Po"

    ```groovy title="main.nf" linenums="30" hl_lines="4"
        ch_samples = channel.fromPath("./data/samples.csv")
            .splitCsv(header: true)
            .map{ row -> separateMetadata(row) }
            .view()
    ```

=== "Przed"

    ```groovy title="main.nf" linenums="30"
        ch_samples = channel.fromPath("./data/samples.csv")
            .splitCsv(header: true)
            .map{ row -> separateMetadata(row) }
    ```

i uruchom workflow:

```bash
nextflow run main.nf
```

??? success "Wyjście polecenia"

    ```console
    [[id:sample_001, organism:human, tissue:liver, depth:30000000, quality:38.5, run:UNSPECIFIED, sample_num:1, lane:001, read:R1, chunk:001, priority:normal], /workspaces/training/side-quests/essential_scripting_patterns/data/sequences/SAMPLE_001_S1_L001_R1_001.fastq]
    [[id:sample_002, organism:mouse, tissue:brain, depth:25000000, quality:35.2, run:UNSPECIFIED, sample_num:2, lane:001, read:R1, chunk:001, priority:normal], /workspaces/training/side-quests/essential_scripting_patterns/data/sequences/SAMPLE_002_S2_L001_R1_001.fastq]
    [[id:sample_003, organism:human, tissue:kidney, depth:45000000, quality:42.1, run:UNSPECIFIED, sample_num:3, lane:001, read:R1, chunk:001, priority:high], /workspaces/training/side-quests/essential_scripting_patterns/data/sequences/SAMPLE_003_S3_L001_R1_001.fastq]
    ```

Idealne! Teraz wszystkie próbki mają pole `run` z ich rzeczywistym ID uruchomienia (wielkimi literami) lub wartością domyślną 'UNSPECIFIED'. Kombinacja `?.` i `?:` zapewnia zarówno bezpieczeństwo (brak awarii), jak i sensowne wartości domyślne.

Usuń teraz operator `.view()`, teraz gdy potwierdziliśmy, że działa.

!!! tip "Łączenie Bezpiecznej Nawigacji i Elvis"

    Wzorzec `value?.method() ?: 'default'` jest powszechny w workflow produkcyjnych:

    - `value?.method()` - Bezpiecznie wywołuje metodę, zwraca `null` jeśli `value` jest `null`
    - `?: 'default'` - Dostarcza wartość zapasową, jeśli wynik jest `null`

    Ten wzorzec obsługuje brakujące/niekompletne dane z gracją.

Używaj tych operatorów konsekwentnie w funkcjach, closures operatorów (`.map{}`, `.filter{}`), skryptach procesów i plikach konfiguracyjnych. Zapobiegają awariom podczas obsługi rzeczywistych danych.

### Podsumowanie

- **Bezpieczna nawigacja (`?.`)**: Zapobiega awariom na wartościach null - zwraca null zamiast rzucać wyjątek
- **Operator Elvis (`?:`)**: Dostarcza wartości domyślne - `value ?: 'default'`
- **Łączenie**: `value?.method() ?: 'default'` to powszechny wzorzec

Te operatory sprawiają, że workflow są odporne na niekompletne dane - niezbędne dla pracy w rzeczywistości.

---

## 7. Walidacja z `error()` i `log.warn`

Czasami musisz natychmiast zatrzymać workflow, jeśli parametry wejściowe są nieprawidłowe. W Nextflow możesz używać wbudowanych funkcji takich jak `error()` i `log.warn`, jak również standardowych konstrukcji programistycznych takich jak instrukcje `if` i logika boolowska, aby implementować logikę walidacji. Dodajmy walidację do naszego workflow.

Utwórz funkcję walidacji przed blokiem workflow, wywołaj ją z workflow i zmień tworzenie kanału, aby używać parametru dla ścieżki pliku CSV. Jeśli parametr brakuje lub plik nie istnieje, wywołaj `error()`, aby zatrzymać wykonanie z jasnym komunikatem.

=== "Po"

    ```groovy title="main.nf" linenums="1" hl_lines="5-20 23-24"
    include { FASTP } from './modules/fastp.nf'
    include { TRIMGALORE } from './modules/trimgalore.nf'
    include { GENERATE_REPORT } from './modules/generate_report.nf'

    def validateInputs() {
        // Sprawdź, czy parametr wejściowy jest podany
        if (!params.input) {
            error("Ścieżka pliku CSV wejściowego nie została podana. Proszę określić --input <file.csv>")
        }

        // Sprawdź, czy plik CSV istnieje
        if (!file(params.input).exists()) {
            error("Plik CSV wejściowy nie został znaleziony: ${params.input}")
        }
    }
    ...
    workflow {
        validateInputs()
        ch_samples = channel.fromPath(params.input)
    ```

=== "Przed"

    ```groovy title="main.nf" linenums="1"
    include { FASTP } from './modules/fastp.nf'
    include { TRIMGALORE } from './modules/trimgalore.nf'
    include { GENERATE_REPORT } from './modules/generate_report.nf'

    ...
    workflow {
        ch_samples = channel.fromPath("./data/samples.csv")
    ```

Teraz spróbuj uruchomić bez pliku CSV:

```bash
nextflow run main.nf
```

??? failure "Wyjście polecenia"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `main.nf` [confident_coulomb] DSL2 - revision: 07059399ed

    WARN: Access to undefined parameter `input` -- Initialise it to a default value eg. `params.input = some_value`
    Ścieżka pliku CSV wejściowego nie została podana. Proszę określić --input <file.csv>
    ```

Workflow zatrzymuje się natychmiast z jasnym komunikatem błędu zamiast zawieść w tajemniczy sposób później

Teraz uruchom z nieistniejącym plikiem:

```bash
nextflow run main.nf --input ./data/nonexistent.csv
```

??? failure "Wyjście polecenia"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `main.nf` [cranky_gates] DSL2 - revision: 26839ae3eb

    Plik CSV wejściowy nie został znaleziony: ./data/nonexistent.csv
    ```

Na koniec uruchom z prawidłowym plikiem:

```bash
nextflow run main.nf --input ./data/samples.csv
```

??? success "Wyjście polecenia"

    ```console
    <!-- TODO: wyjście -->
    ```

Tym razem uruchamia się pomyślnie.

Możesz również dodać walidację wewnątrz funkcji `separateMetadata`. Użyjmy niekrytycznego `log.warn` do wydawania ostrzeżeń dla próbek o niskiej głębokości sekwencjonowania, ale nadal pozwólmy workflow kontynuować:

=== "Po"

    ```groovy title="main.nf" linenums="1" hl_lines="3-6"
        def priority = sample_meta.quality > 40 ? 'high' : 'normal'

        // Waliduj, czy dane mają sens
        if (sample_meta.depth < 30000000) {
            log.warn "Niska głębokość sekwencjonowania dla ${sample_meta.id}: ${sample_meta.depth}"
        }

        return tuple(sample_meta + file_meta + [priority: priority], fastq_path)
    }
    ```

=== "Przed"

    ```groovy title="main.nf" linenums="1"
        def priority = sample_meta.quality > 40 ? 'high' : 'normal'

        return tuple(sample_meta + file_meta + [priority: priority], fastq_path)
    }
    ```

Uruchom workflow ponownie z oryginalnym CSV:

```bash
nextflow run main.nf --input ./data/samples.csv
```

??? warning "Wyjście polecenia"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `main.nf` [awesome_goldwasser] DSL2 - revision: a31662a7c1

    executor >  local (5)
    [ce/df5eeb] process > FASTP (2)           [100%] 2 of 2 ✔
    [-        ] process > TRIMGALORE          -
    [d1/7d2b4b] process > GENERATE_REPORT (3) [100%] 3 of 3 ✔
    WARN: Niska głębokość sekwencjonowania dla sample_002: 25000000
    ```

Widzimy ostrzeżenie o niskiej głębokości sekwencjonowania dla jednej z próbek.

### Podsumowanie

- **`error()`**: Natychmiast zatrzymuje workflow z jasnym komunikatem
- **`log.warn`**: Wydaje ostrzeżenia bez zatrzymywania workflow
- **Wczesna walidacja**: Sprawdzaj wejścia przed przetwarzaniem, aby szybko zawodzić z pomocnymi błędami
- **Funkcje walidacyjne**: Twórz logikę walidacji wielokrotnego użytku, którą można wywoływać na początku workflow

Odpowiednia walidacja sprawia, że workflow są bardziej odporne i przyjazne dla użytkownika, wyłapując problemy wcześnie z jasnymi komunikatami błędów.

---

## 8. Handlery Zdarzeń Workflow

Do tej pory pisaliśmy kod w naszych skryptach workflow i definicjach procesów. Ale jest jeszcze jedna ważna funkcja, o której powinieneś wiedzieć: handlery zdarzeń workflow.

Handlery zdarzeń to closures, które uruchamiają się w określonych punktach cyklu życia Twojego workflow. Są idealne do dodawania logowania, powiadomień lub operacji porządkowych. Te handlery powinny być zdefiniowane w Twoim skrypcie workflow obok definicji workflow.

### 8.1. Handler `onComplete`

Najczęściej używanym handlerem zdarzeń jest `onComplete`, który uruchamia się, gdy Twój workflow się kończy (czy zakończył się sukcesem, czy niepowodzeniem). Dodajmy jeden, aby podsumować wyniki naszego pipeline.

Dodaj handler zdarzeń do Twojego pliku `main.nf`, wewnątrz definicji workflow:

=== "Po"

    ```groovy title="main.nf" linenums="66" hl_lines="5-16"
        ch_fastp = FASTP(trim_branches.fastp)
        ch_trimgalore = TRIMGALORE(trim_branches.trimgalore)
        GENERATE_REPORT(ch_samples)

        workflow.onComplete = {
            println ""
            println "Podsumowanie wykonania pipeline:"
            println "=========================="
            println "Zakończono: ${workflow.complete}"
            println "Czas trwania: ${workflow.duration}"
            println "Sukces     : ${workflow.success}"
            println "katalog roboczy: ${workflow.workDir}"
            println "status wyjścia : ${workflow.exitStatus}"
            println ""
        }
    }
    ```

=== "Przed"

    ```groovy title="main.nf" linenums="66" hl_lines="4"
        ch_fastp = FASTP(trim_branches.fastp)
        ch_trimgalore = TRIMGALORE(trim_branches.trimgalore)
        GENERATE_REPORT(ch_samples)
    }
    ```

To closure uruchamia się, gdy workflow się kończy. Wewnątrz masz dostęp do obiektu `workflow`, który dostarcza użytecznych właściwości o wykonaniu.

Uruchom swój workflow, a zobaczysz to podsumowanie pojawiające się na końcu!

```bash
nextflow run main.nf --input ./data/samples.csv -ansi-log false
```

??? success "Wyjście polecenia"

    ```console
    N E X T F L O W  ~  version 25.10.2
    Launching `main.nf` [marvelous_boltzmann] DSL2 - revision: a31662a7c1
    WARN: Niska głębokość sekwencjonowania dla sample_002: 25000000
    [9b/d48e40] Submitted process > FASTP (2)
    [6a/73867a] Submitted process > GENERATE_REPORT (2)
    [79/ad0ac5] Submitted process > GENERATE_REPORT (1)
    [f3/bda6cb] Submitted process > FASTP (1)
    [34/d5b52f] Submitted process > GENERATE_REPORT (3)

    Podsumowanie wykonania pipeline:
    ==========================
    Zakończono: 2025-10-10T12:14:24.885384+01:00
    Czas trwania: 2.9s
    Sukces     : true
    katalog roboczy: /workspaces/training/side-quests/essential_scripting_patterns/work
    status wyjścia : 0
    ```

Uczyńmy to bardziej użytecznym, dodając logikę warunkową:

=== "Po"

    ```groovy title="main.nf" linenums="66" hl_lines="5-22"
        ch_fastp = FASTP(trim_branches.fastp)
        ch_trimgalore = TRIMGALORE(trim_branches.trimgalore)
        GENERATE_REPORT(ch_samples)

        workflow.onComplete = {
            println ""
            println "Podsumowanie wykonania pipeline:"
            println "=========================="
            println "Zakończono: ${workflow.complete}"
            println "Czas trwania: ${workflow.duration}"
            println "Sukces     : ${workflow.success}"
            println "katalog roboczy: ${workflow.workDir}"
            println "status wyjścia : ${workflow.exitStatus}"
            println ""

            if (workflow.success) {
                println "✅ Pipeline zakończył się pomyślnie!"
            } else {
                println "❌ Pipeline zawiódł!"
                println "Błąd: ${workflow.errorMessage}"
            }
        }
    }
    ```

=== "Przed"

    ```groovy title="main.nf" linenums="66" hl_lines="5-16"
        ch_fastp = FASTP(trim_branches.fastp)
        ch_trimgalore = TRIMGALORE(trim_branches.trimgalore)
        GENERATE_REPORT(ch_samples)

        workflow.onComplete = {
            println ""
            println "Podsumowanie wykonania pipeline:"
            println "=========================="
            println "Zakończono: ${workflow.complete}"
            println "Czas trwania: ${workflow.duration}"
            println "Sukces     : ${workflow.success}"
            println "katalog roboczy: ${workflow.workDir}"
            println "status wyjścia : ${workflow.exitStatus}"
            println ""
        }
    }
    ```

Teraz otrzymujemy jeszcze bardziej informacyjne podsumowanie, zawierające komunikat sukcesu/niepowodzenia i katalog wyjściowy, jeśli został określony:

<!-- TODO: dodaj polecenie uruchomienia -->

??? success "Wyjście polecenia"

    ```console
    N E X T F L O W  ~  version 25.10.2
    Launching `main.nf` [boring_linnaeus] DSL2 - revision: a31662a7c1
    WARN: Niska głębokość sekwencjonowania dla sample_002: 25000000
    [e5/242efc] Submitted process > FASTP (2)
    [3b/74047c] Submitted process > GENERATE_REPORT (3)
    [8a/7a57e6] Submitted process > GENERATE_REPORT (1)
    [a8/b1a31f] Submitted process > GENERATE_REPORT (2)
    [40/648429] Submitted process > FASTP (1)

    Podsumowanie wykonania pipeline:
    ==========================
    Zakończono: 2025-10-10T12:16:00.522569+01:00
    Czas trwania: 3.6s
    Sukces     : true
    katalog roboczy: /workspaces/training/side-quests/essential_scripting_patterns/work
    status wyjścia : 0

    ✅ Pipeline zakończył się pomyślnie!
    ```

Możesz również zapisać podsumowanie do pliku używając operacji na plikach:

```groovy title="main.nf - Zapisywanie podsumowania do pliku"
workflow {
    // ... kod Twojego workflow ...

    workflow.onComplete = {
        def summary = """
        Podsumowanie Wykonania Pipeline
        ===========================
        Zakończono: ${workflow.complete}
        Czas trwania: ${workflow.duration}
        Sukces  : ${workflow.success}
        Polecenie: ${workflow.commandLine}
        """

        println summary

        // Zapisz do pliku logów
        def log_file = file("${workflow.launchDir}/pipeline_summary.txt")
        log_file.text = summary
    }
}
```

### 8.2. Handler `onError`

Poza `onComplete`, jest jeszcze jeden handler zdarzeń, którego możesz użyć: `onError`, który uruchamia się tylko, jeśli workflow zawiedzie:

```groovy title="main.nf - handler onError"
workflow {
    // ... kod Twojego workflow ...

    workflow.onError = {
        println "="* 50
        println "Wykonanie pipeline zawiodło!"
        println "Komunikat błędu: ${workflow.errorMessage}"
        println "="* 50

        // Zapisz szczegółowy log błędów
        def error_file = file("${workflow.launchDir}/error.log")
        error_file.text = """
        Raport Błędu Workflow
        =====================
        Czas: ${new Date()}
        Błąd: ${workflow.errorMessage}
        Raport błędu: ${workflow.errorReport ?: 'Brak szczegółowego raportu'}
        """

        println "Szczegóły błędu zapisane do: ${error_file}"
    }
}
```

Możesz używać wielu handlerów razem w swoim skrypcie workflow:

```groovy title="main.nf - Połączone handlery"
workflow {
    // ... kod Twojego workflow ...

    workflow.onError = {
        println "Workflow zawiódł: ${workflow.errorMessage}"
    }

    workflow.onComplete = {
        def duration_mins = workflow.duration.toMinutes().round(2)
        def status = workflow.success ? "SUKCES ✅" : "NIEPOWODZENIE ❌"

        println """
        Pipeline zakończony: ${status}
        Czas trwania: ${duration_mins} minut
        """
    }
}
```

### Podsumowanie

W tej sekcji nauczyłeś się:

- **Closures handlerów zdarzeń**: Closures w Twoim skrypcie workflow, które uruchamiają się w różnych punktach cyklu życia
- **Handler `onComplete`**: Do podsumowań wykonania i raportowania wyników
- **Handler `onError`**: Do obsługi błędów i logowania niepowodzeń
- **Właściwości obiektu workflow**: Dostęp do `workflow.success`, `workflow.duration`, `workflow.errorMessage`, itp.

Handlery zdarzeń pokazują, jak możesz używać pełnej mocy języka Nextflow w swoich skryptach workflow, aby dodać zaawansowane możliwości logowania i powiadamiania.

---

## Podsumowanie

Gratulacje, udało Ci się!

W trakcie tego side questa zbudowałeś kompleksowy pipeline przetwarzania próbek, który ewoluował od podstawowej obsługi metadanych do zaawansowanego, gotowego do produkcji workflow.
Każda sekcja budowała na poprzedniej, demonstrując, jak konstrukcje programistyczne przekształcają proste workflow w potężne systemy przetwarzania danych, z następującymi korzyściami:

- **Czystszy kod**: Zrozumienie przepływu danych vs skryptowania pomaga pisać bardziej zorganizowane workflow
- **Odporna obsługa**: Bezpieczna nawigacja i operatory Elvis sprawiają, że workflow są odporne na brakujące dane
- **Elastyczne przetwarzanie**: Logika warunkowa pozwala Twoim workflow odpowiednio przetwarzać różne typy próbek
- **Adaptacyjne zasoby**: Dynamiczne dyrektywy optymalizują użycie zasobów na podstawie charakterystyki wejścia

Ta progresja odzwierciedla rzeczywistą ewolucję pipeline'ów bioinformatycznych, od prototypów badawczych obsługujących kilka próbek do systemów produkcyjnych przetwarzających tysiące próbek w laboratoriach i instytucjach.
Każde wyzwanie, które rozwiązałeś, i wzorzec, którego się nauczyłeś, odzwierciedla faktyczne problemy, z którymi borykają się deweloperzy przy skalowaniu workflow Nextflow.

Zastosowanie tych wzorców w Twojej własnej pracy umożliwi Ci budowanie odpornych, gotowych do produkcji workflow.

### Kluczowe wzorce

1.  **Przepływ Danych vs Skryptowanie:** Nauczyłeś się rozróżniać między operacjami przepływu danych (orkiestracja kanałów) a skryptowaniem (kod, który manipuluje danymi), włączając w to kluczowe różnice między operacjami na różnych typach jak `collect` na Channel vs List.

    - Przepływ danych: orkiestracja kanałów

    ```groovy
    channel.fromPath('*.fastq').splitCsv(header: true)
    ```

    - Skryptowanie: przetwarzanie danych na kolekcjach

    ```groovy
    sample_data.collect { it.toUpperCase() }
    ```

2.  **Zaawansowane Przetwarzanie Ciągów Znaków**: Opanowałeś wyrażenia regularne do parsowania nazw plików, dynamiczne generowanie skryptów w procesach oraz interpolację zmiennych (Nextflow vs Bash vs Shell).

    - Dopasowywanie wzorców

    ```groovy
    filename =~ ~/^(\w+)_(\w+)_(\d+)\.fastq$/
    ```

    - Funkcja z warunkowym zwrotem

    ```groovy
    def parseSample(filename) {
        def matcher = filename =~ pattern
        return matcher ? [valid: true, data: matcher[0]] : [valid: false]
    }
    ```

    - Kolekcja plików do argumentów polecenia (w bloku skryptowym procesu)

    ```groovy
    script:
    def file_args = input_files.collect { file -> "--input ${file}" }.join(' ')
    """
    analysis_tool ${file_args} --output results.txt
    """
    ```

3.  **Tworzenie Funkcji Wielokrotnego Użytku**: Nauczyłeś się wyodrębniać złożoną logikę do nazwanych funkcji, które można wywoływać z operatorów kanałów, czyniąc workflow bardziej czytelnymi i łatwymi w utrzymaniu.

    - Zdefiniuj nazwaną funkcję

    ```groovy
    def separateMetadata(row) {
        def sample_meta = [ /* kod ukryty dla zwięzłości */ ]
        def fastq_path = file(row.file_path)
        def m = (fastq_path.name =~ /^(.+)_S(\d+)_L(\d{3})_(R[12])_(\d{3})\.fastq(?:\.gz)?$/)
        def file_meta = m ? [ /* kod ukryty dla zwięzłości */ ] : [:]
        def priority = sample_meta.quality > 40 ? 'high' : 'normal'

        return tuple(sample_meta + file_meta + [priority: priority], fastq_path)
    }
    ```

    - Wywołaj nazwaną funkcję w workflow

    ```groovy
    workflow {
        ch_samples = channel.fromPath("./data/samples.csv")
            .splitCsv(header: true)
            .map{ row -> separateMetadata(row) }

        ch_fastp = FASTP(ch_samples)
    }
    ```

4.  **Dynamiczne Dyrektywy Zasobów z Closures**: Zbadałeś używanie closures w dyrektywach procesów do adaptacyjnej alokacji zasobów na podstawie charakterystyki wejścia.

    - Nazwane closures i kompozycja

    ```groovy
    def enrichData = normalizeId >> addQualityCategory >> addFlags
    def processor = generalFunction.curry(fixedParam)
    ```

    - Closures z dostępem do zakresu

    ```groovy
    def collectStats = { data -> stats.count++; return data }
    ```

5.  **Logika Warunkowa i Kontrola Procesów**: Dodałeś inteligentne kierowanie używając operatorów `.branch()` i `.filter()`, wykorzystując prawdziwość do zwięzłych wyrażeń warunkowych.

    - Użyj `.branch()` do kierowania danych przez różne gałęzie workflow

    ```groovy
    trim_branches = ch_samples
    .branch { meta, reads ->
        fastp: meta.organism == 'human' && meta.depth >= 30000000
        trimgalore: true
    }

    ch_fastp = FASTP(trim_branches.fastp)
    ch_trimgalore = TRIMGALORE(trim_branches.trimgalore)
    ```

    - Ewaluacja boolowska z Groovy Truth

    ```groovy
    if (sample.files) println "Ma pliki"
    ```

    - Użyj `filter()` do podzbioru danych z 'prawdziwością'

    ```groovy
    ch_valid_samples = ch_samples
        .filter { meta, reads ->
            meta.id && meta.organism && meta.depth >= 25000000
        }
    ```

6.  **Bezpieczna Nawigacja i Operatory Elvis**: Uczyniłeś pipeline odpornym na brakujące dane używając `?.` do bezpiecznego dostępu do właściwości względem null i `?:` do dostarczania wartości domyślnych.

    ```groovy
    def id = data?.sample?.id ?: 'unknown'
    ```

7.  **Walidacja z error() i log.warn**: Nauczyłeś się walidować wejścia wcześnie i szybko zawodzić z jasnymi komunikatami błędów.

    ```groovy
    try {
        def errors = validateSample(sample)
        if (errors) throw new RuntimeException("Nieprawidłowy: ${errors.join(', ')}")
    } catch (Exception e) {
        println "Błąd: ${e.message}"
    }
    ```

8.  **Handlery Zdarzeń Konfiguracji**: Nauczyłeś się używać handlerów zdarzeń workflow (`onComplete` i `onError`) do logowania, powiadomień i zarządzania cyklem życia.

    - Używanie `onComplete` do logowania i powiadamiania

    ```groovy
    workflow.onComplete = {
        println "Sukces     : ${workflow.success}"
        println "status wyjścia : ${workflow.exitStatus}"

        if (workflow.success) {
            println "✅ Pipeline zakończył się pomyślnie!"
        } else {
            println "❌ Pipeline zawiódł!"
            println "Błąd: ${workflow.errorMessage}"
        }
    }
    ```

    - Używanie `onError` do podejmowania działań specyficznie w przypadku niepowodzenia

    ```groovy
    workflow.onError = {
        // Zapisz szczegółowy log błędów
        def error_file = file("${workflow.launchDir}/error.log")
        error_file.text = """
        Czas: ${new Date()}
        Błąd: ${workflow.errorMessage}
        Raport błędu: ${workflow.errorReport ?: 'Brak szczegółowego raportu'}
        """

        println "Szczegóły błędu zapisane do: ${error_file}"
    }
    ```

### Dodatkowe zasoby

- [Dokumentacja Języka Nextflow](https://nextflow.io/docs/latest/reference/syntax.html)
- [Operatory Nextflow](https://www.nextflow.io/docs/latest/operator.html)
- [Składnia Skryptów Nextflow](https://www.nextflow.io/docs/latest/script.html)
- [Biblioteka Standardowa Nextflow](https://nextflow.io/docs/latest/reference/stdlib.html)

Pamiętaj, aby sprawdzić te zasoby, gdy będziesz musiał zbadać bardziej zaawansowane funkcje.

Będziesz korzystać z ćwiczeń i rozwijania swoich umiejętności, aby:

- Pisać czystsze workflow z odpowiednim oddzieleniem między przepływem danych a skryptowaniem
- Opanować interpolację zmiennych, aby uniknąć powszechnych pułapek ze zmiennymi Nextflow, Bash i powłoki
- Używać dynamicznych dyrektyw zasobów dla wydajnych, adaptacyjnych workflow
- Przekształcać kolekcje plików w prawidłowo sformatowane argumenty linii poleceń
- Obsługiwać różne konwencje nazewnictwa plików i formaty wejściowe z gracją, używając regex i przetwarzania ciągów znaków
- Budować wielokrotnie używalny, łatwy w utrzymaniu kod używając zaawansowanych wzorców closure i programowania funkcyjnego
- Przetwarzać i organizować złożone zbiory danych używając operacji na kolekcjach
- Dodawać walidację, obsługę błędów i logowanie, aby uczynić swoje workflow gotowymi do produkcji
- Implementować zarządzanie cyklem życia workflow za pomocą handlerów zdarzeń

---

## Co dalej?

Wróć do [menu Side Quests](./index.md) lub kliknij przycisk w prawym dolnym rogu strony, aby przejść do następnego tematu na liście.
