# Część 1: Uruchomienie demonstracyjnego pipeline'a

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Tłumaczenie wspomagane przez AI - [dowiedz się więcej i zasugeruj ulepszenia](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

W tej pierwszej części szkolenia Hello nf-core pokażemy Ci, jak znaleźć i wypróbować pipeline nf-core, zrozumieć, jak zorganizowany jest kod oraz rozpoznać, czym różni się on od zwykłego kodu Nextflow przedstawionego w [Hello Nextflow](../hello_nextflow/index.md).

Będziemy używać pipeline'a o nazwie nf-core/demo, który jest utrzymywany przez projekt nf-core jako część inwentarza pipeline'ów służących do demonstracji struktury kodu i operacji narzędzi.

Upewnij się, że Twój katalog roboczy jest ustawiony na `hello-nf-core/`, jak wskazano na stronie [Pierwsze kroki](./00_orientation.md).

---

## 1. Znajdź i pobierz pipeline nf-core/demo

Zacznijmy od zlokalizowania pipeline'a nf-core/demo na stronie projektu pod adresem [nf-co.re](https://nf-co.re), która centralizuje wszystkie informacje, takie jak: ogólna dokumentacja i artykuły pomocnicze, dokumentacja dla każdego z pipeline'ów, wpisy na blogu, ogłoszenia wydarzeń i tak dalej.

### 1.1. Znajdź pipeline na stronie internetowej

W przeglądarce internetowej przejdź do [https://nf-co.re/pipelines/](https://nf-co.re/pipelines/) i wpisz `demo` w pasku wyszukiwania.

![wyniki wyszukiwania](./img/search-results.png)

Kliknij nazwę pipeline'a, `demo`, aby uzyskać dostęp do strony dokumentacji pipeline'a.

Każdy wydany pipeline ma dedykowaną stronę zawierającą następujące sekcje dokumentacji:

- **Introduction:** Wprowadzenie i przegląd pipeline'a
- **Usage:** Opisy sposobu wykonywania pipeline'a
- **Parameters:** Zgrupowane parametry pipeline'a z opisami
- **Output:** Opisy i przykłady oczekiwanych plików wyjściowych
- **Results:** Przykładowe pliki wyjściowe wygenerowane z pełnego zestawu danych testowych
- **Releases & Statistics:** Historia wersji pipeline'a i statystyki

Gdy rozważasz przyjęcie nowego pipeline'a, powinieneś najpierw uważnie przeczytać jego dokumentację, aby zrozumieć, co robi i jak powinien być skonfigurowany, zanim spróbujesz go uruchomić.

Spójrz teraz i zobacz, czy możesz dowiedzieć się:

- Które narzędzia uruchomi pipeline (Sprawdź zakładkę: `Introduction`)
- Jakie wejścia i parametry pipeline akceptuje lub wymaga (Sprawdź zakładkę: `Parameters`)
- Jakie są wyjścia produkowane przez pipeline (Sprawdź zakładkę: `Output`)

#### 1.1.1. Przegląd pipeline'a

Zakładka `Introduction` dostarcza przegląd pipeline'a, w tym wizualną reprezentację (zwaną mapą metra) i listę narzędzi uruchamianych jako część pipeline'a.

![mapa metra pipeline'a](./img/nf-core-demo-subway-cropped.png)

1. Read QC (FASTQC)
2. Adapter and quality trimming (SEQTK_TRIM)
3. Present QC for raw reads (MULTIQC)

#### 1.1.2. Przykładowa linia poleceń

Dokumentacja dostarcza również przykładowy plik wejściowy (omówiony dokładniej poniżej) i przykładową linię poleceń.

```bash
nextflow run nf-core/demo \
  -profile <docker/singularity/.../institute> \
  --input samplesheet.csv \
  --outdir <OUTDIR>
```

Zauważysz, że przykładowe polecenie NIE określa pliku workflow'u, tylko referencję do repozytorium pipeline'a, `nf-core/demo`.

Gdy jest wywoływany w ten sposób, Nextflow zakłada, że kod jest zorganizowany w określony sposób.
Pobierzmy kod, abyśmy mogli zbadać tę strukturę.

### 1.2. Pobierz kod pipeline'a

Gdy ustaliliśmy już, że pipeline wydaje się być odpowiedni dla naszych celów, wypróbujmy go.
Na szczęście Nextflow ułatwia pobieranie pipeline'ów z poprawnie sformatowanych repozytoriów bez konieczności ręcznego pobierania czegokolwiek.

Wróćmy do terminala i uruchommy następujące polecenie:

```bash
nextflow pull nf-core/demo
```

??? success "Wyjście polecenia"

    ```console
    Checking nf-core/demo ...
    downloaded from https://github.com/nf-core/demo.git - revision: 04060b4644 [master]
    ```

Nextflow wykonuje `pull` kodu pipeline'a, co oznacza, że pobiera całe repozytorium na Twój dysk lokalny.

Wyjaśnijmy, że możesz to zrobić z każdym workflow'em Nextflow, który jest odpowiednio skonfigurowany w GitHub, nie tylko z pipeline'ami nf-core.
Jednak nf-core to największa otwarta kolekcja workflow'ów Nextflow.

Możesz uzyskać od Nextflow listę pipeline'ów, które pobrałeś w ten sposób:

```bash
nextflow list
```

??? success "Wyjście polecenia"

    ```console
    nf-core/demo
    ```

Zauważysz, że pliki nie znajdują się w Twoim bieżącym katalogu roboczym.
Domyślnie Nextflow zapisuje je w `$NXF_HOME/assets`.

```bash
tree -L 2 $NXF_HOME/assets/
```

```console title="Zawartość katalogu"
/workspaces/.nextflow/assets/
└── nf-core
    └── demo

2 directories, 0 files
```

!!! note

    Pełna ścieżka może się różnić w Twoim systemie, jeśli nie używasz naszego środowiska szkoleniowego.

Nextflow celowo trzyma pobrany kod źródłowy 'z dala od drogi' w oparciu o zasadę, że te pipeline'y powinny być używane bardziej jak biblioteki niż kod, z którym bezpośrednio współdziałasz.

Jednak dla celów tego szkolenia chcemy móc zagłębiać się i zobaczyć, co tam jest.
Więc aby to ułatwić, stwórzmy dowiązanie symboliczne do tej lokalizacji z naszego bieżącego katalogu roboczego.

```bash
ln -s $NXF_HOME/assets pipelines
```

To tworzy skrót, który ułatwia eksplorację kodu, który właśnie pobraliśmy.

```bash
tree -L 2 pipelines
```

```console title="Zawartość katalogu"
pipelines
└── nf-core
    └── demo

2 directories, 0 files
```

Teraz możemy łatwiej zajrzeć do kodu źródłowego w razie potrzeby.

Ale najpierw spróbujmy uruchomić nasz pierwszy pipeline nf-core!

### Podsumowanie

Wiesz, jak znaleźć pipeline za pośrednictwem strony nf-core i pobrać lokalną kopię kodu źródłowego.

### Co dalej?

Dowiedz się, jak wypróbować pipeline nf-core przy minimalnym wysiłku.

---

## 2. Wypróbuj pipeline z jego profilem testowym

Wygodnie, każdy pipeline nf-core jest dostarczany z profilem testowym.
Jest to minimalny zestaw ustawień konfiguracyjnych dla pipeline'a do uruchomienia z użyciem małego zestawu danych testowych hostowanego w repozytorium [nf-core/test-datasets](https://github.com/nf-core/test-datasets).
To świetny sposób, aby szybko wypróbować pipeline na małą skalę.

!!! note

    System profili konfiguracyjnych Nextflow pozwala łatwo przełączać się między różnymi silnikami kontenerów lub środowiskami wykonawczymi.
    Aby uzyskać więcej szczegółów, zobacz [Hello Nextflow Część 6: Konfiguracja](../hello_nextflow/06_hello_config.md).

### 2.1. Zbadaj profil testowy

Dobrą praktyką jest sprawdzenie, co określa profil testowy pipeline'a przed jego uruchomieniem.
Profil `test` dla `nf-core/demo` znajduje się w pliku konfiguracyjnym `conf/test.config` i jest pokazany poniżej.

```groovy title="conf/test.config" linenums="1" hl_lines="8 26"
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Plik konfiguracyjny Nextflow do uruchamiania minimalnych testów
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Definiuje pliki wejściowe i wszystko, co jest wymagane do uruchomienia szybkiego i prostego testu pipeline'a.

    Użyj w następujący sposób:
        nextflow run nf-core/demo -profile test,<docker/singularity> --outdir <OUTDIR>

----------------------------------------------------------------------------------------
*/

process {
    resourceLimits = [
        cpus: 4,
        memory: '4.GB',
        time: '1.h'
    ]
}

params {
    config_profile_name        = 'Test profile'
    config_profile_description = 'Minimal test dataset to check pipeline function'

    // Dane wejściowe
    input  = 'https://raw.githubusercontent.com/nf-core/test-datasets/viralrecon/samplesheet/samplesheet_test_illumina_amplicon.csv'

}
```

Od razu zauważysz, że blok komentarza na górze zawiera przykład użycia pokazujący, jak uruchomić pipeline z tym profilem testowym.

```groovy title="conf/test.config" linenums="7"
Użyj w następujący sposób:
        nextflow run nf-core/demo -profile test,<docker/singularity> --outdir <OUTDIR>
```

Jedyne rzeczy, które musimy dostarczyć, to to, co jest pokazane między nawiasami kątowymi w przykładowym poleceniu: `<docker/singularity>` i `<OUTDIR>`.

Przypominając, `<docker/singularity>` odnosi się do wyboru systemu kontenerów. Wszystkie pipeline'y nf-core są zaprojektowane tak, aby były użyteczne z kontenerami (Docker, Singularity, itp.) w celu zapewnienia powtarzalności i eliminacji problemów z instalacją oprogramowania.
Więc będziemy musieli określić, czy chcemy użyć Docker czy Singularity do testowania pipeline'a.

Część `--outdir <OUTDIR>` odnosi się do katalogu, w którym Nextflow zapisze wyjścia pipeline'a.
Musimy podać dla niego nazwę, którą możemy po prostu wymyślić.
Jeśli nie istnieje już, Nextflow utworzy go dla nas w czasie wykonywania.

Przechodząc do sekcji po bloku komentarza, profil testowy pokazuje nam, co zostało wstępnie skonfigurowane do testowania: przede wszystkim parametr `input` jest już ustawiony tak, aby wskazywał na zestaw danych testowych, więc nie musimy dostarczać własnych danych.
Jeśli podążysz za linkiem do wstępnie skonfigurowanego wejścia, zobaczysz, że jest to plik csv zawierający identyfikatory próbek i ścieżki plików dla kilku próbek eksperymentalnych.

```csv title="samplesheet_test_illumina_amplicon.csv"
sample,fastq_1,fastq_2
SAMPLE1_PE,https://raw.githubusercontent.com/nf-core/test-datasets/viralrecon/illumina/amplicon/sample1_R1.fastq.gz,https://raw.githubusercontent.com/nf-core/test-datasets/viralrecon/illumina/amplicon/sample1_R2.fastq.gz
SAMPLE2_PE,https://raw.githubusercontent.com/nf-core/test-datasets/viralrecon/illumina/amplicon/sample2_R1.fastq.gz,https://raw.githubusercontent.com/nf-core/test-datasets/viralrecon/illumina/amplicon/sample2_R2.fastq.gz
SAMPLE3_SE,https://raw.githubusercontent.com/nf-core/test-datasets/viralrecon/illumina/amplicon/sample1_R1.fastq.gz,
SAMPLE3_SE,https://raw.githubusercontent.com/nf-core/test-datasets/viralrecon/illumina/amplicon/sample2_R1.fastq.gz,
```

Nazywa się to samplesheet i jest to najczęstsza forma wejścia do pipeline'ów nf-core.

!!! note

    Nie martw się, jeśli nie jesteś zaznajomiony z formatami i typami danych, nie jest to ważne dla tego, co następuje.

Więc potwierdza to, że mamy wszystko, czego potrzebujemy, aby wypróbować pipeline.

### 2.2. Uruchom pipeline

Zdecydujmy się użyć Docker dla systemu kontenerów i `demo-results` jako katalogu wyjściowego, i jesteśmy gotowi do uruchomienia polecenia testowego:

```bash
nextflow run nf-core/demo -profile docker,test --outdir demo-results
```

??? success "Wyjście polecenia"

    ```console
     N E X T F L O W   ~  version 25.04.3

    Launching `https://github.com/nf-core/demo` [magical_pauling] DSL2 - revision: db7f526ce1 [master]


    ------------------------------------------------------
                                            ,--./,-.
            ___     __   __   __   ___     /,-._.--~'
      |\ | |__  __ /  ` /  \ |__) |__         }  {
      | \| |       \__, \__/ |  \ |___     \`-._,-`-,
                                            `._,._,'
      nf-core/demo 1.0.2
    ------------------------------------------------------
    Input/output options
      input                     : https://raw.githubusercontent.com/nf-core/test-datasets/viralrecon/samplesheet/samplesheet_test_illumina_amplicon.csv
      outdir                    : demo-results

    Institutional config options
      config_profile_name       : Test profile
      config_profile_description: Minimal test dataset to check pipeline function

    Generic options
      trace_report_suffix       : 2025-11-21_04-57-41

    Core Nextflow options
      revision                  : master
      runName                   : magical_pauling
      containerEngine           : docker
      launchDir                 : /workspaces/training/hello-nf-core
      workDir                   : /workspaces/training/hello-nf-core/work
      projectDir                : /workspaces/.nextflow/assets/nf-core/demo
      userName                  : root
      profile                   : docker,test
      configFiles               : /workspaces/.nextflow/assets/nf-core/demo/nextflow.config

    !! Only displaying parameters that differ from the pipeline defaults !!
    ------------------------------------------------------
    * The pipeline
        https://doi.org/10.5281/zenodo.12192442

    * The nf-core framework
        https://doi.org/10.1038/s41587-020-0439-x

    * Software dependencies
        https://github.com/nf-core/demo/blob/master/CITATIONS.md


    executor >  local (7)
    [ff/a6976b] NFCORE_DEMO:DEMO:FASTQC (SAMPLE3_SE)     | 3 of 3 ✔
    [39/731ab7] NFCORE_DEMO:DEMO:SEQTK_TRIM (SAMPLE3_SE) | 3 of 3 ✔
    [7c/78d96e] NFCORE_DEMO:DEMO:MULTIQC                 | 1 of 1 ✔
    -[nf-core/demo] Pipeline completed successfully-
    ```

Jeśli Twoje wyjście pasuje do tego, gratulacje! Właśnie uruchomiłeś Swój pierwszy pipeline nf-core.

Zauważysz, że jest znacznie więcej wyjścia konsoli niż podczas uruchamiania podstawowego pipeline'a Nextflow.
Jest nagłówek, który zawiera podsumowanie wersji pipeline'a, wejść i wyjść oraz kilku elementów konfiguracji.

!!! note

    Twoje wyjście pokaże różne znaczniki czasu, nazwy wykonań i ścieżki plików, ale ogólna struktura i wykonanie procesów powinny być podobne.

Przechodząc do wyjścia wykonania, spójrzmy na linie, które mówią nam, jakie procesy zostały uruchomione:

```console
    [ff/a6976b] NFCORE_DEMO:DEMO:FASTQC (SAMPLE3_SE)     | 3 of 3 ✔
    [39/731ab7] NFCORE_DEMO:DEMO:SEQTK_TRIM (SAMPLE3_SE) | 3 of 3 ✔
    [7c/78d96e] NFCORE_DEMO:DEMO:MULTIQC                 | 1 of 1 ✔
```

To mówi nam, że zostały uruchomione trzy procesy, odpowiadające trzem narzędziom pokazanym na stronie dokumentacji pipeline'a na stronie nf-core: FASTQC, SEQTK_TRIM i MULTIQC.

Pełne nazwy procesów, jak pokazano tutaj, takie jak `NFCORE_DEMO:DEMO:MULTIQC`, są dłuższe niż to, co mogłeś zobaczyć we wstępnym materiale Hello Nextflow.
Zawierają one nazwy ich workflow'ów nadrzędnych i odzwierciedlają modularność kodu pipeline'a.
Zajmiemy się tym bardziej szczegółowo za chwilę.

### 2.3. Zbadaj wyjścia pipeline'a

Na koniec spójrzmy na katalog `demo-results` wyprodukowany przez pipeline.

```bash
tree -L 2 demo-results
```

??? abstract "Zawartość katalogu"

    ```console
    demo-results
    ├── fastqc
    │   ├── SAMPLE1_PE
    │   ├── SAMPLE2_PE
    │   └── SAMPLE3_SE
    ├── fq
    │   ├── SAMPLE1_PE
    │   ├── SAMPLE2_PE
    │   └── SAMPLE3_SE
    ├── multiqc
    │   ├── multiqc_data
    │   ├── multiqc_plots
    │   └── multiqc_report.html
    └── pipeline_info
        ├── execution_report_2025-11-21_04-57-41.html
        ├── execution_timeline_2025-11-21_04-57-41.html
        ├── execution_trace_2025-11-21_04-57-41.txt
        ├── nf_core_demo_software_mqc_versions.yml
        ├── params_2025-11-21_04-57-46.json
        └── pipeline_dag_2025-11-21_04-57-41.html
    ```

To może wydawać się dużo.
Aby dowiedzieć się więcej o wyjściach pipeline'a `nf-core/demo`, sprawdź jego [stronę dokumentacji](https://nf-co.re/demo/1.0.2/docs/output/).

Na tym etapie ważne jest zaobserwowanie, że wyniki są zorganizowane według modułu, a dodatkowo istnieje katalog o nazwie `pipeline_info` zawierający różne raporty z datami dotyczące wykonania pipeline'a.

Na przykład plik `execution_timeline_*` pokazuje, jakie procesy zostały uruchomione, w jakiej kolejności i jak długo trwało ich uruchomienie:

![raport osi czasu wykonania](./img/execution_timeline.png)

!!! note

    Tutaj zadania nie zostały uruchomione równolegle, ponieważ działamy na minimalistycznej maszynie w Github Codespaces.
    Aby zobaczyć ich równoległe uruchomienie, spróbuj zwiększyć alokację CPU Swojego codespace i limity zasobów w konfiguracji testowej.

Te raporty są generowane automatycznie dla wszystkich pipeline'ów nf-core.

### Podsumowanie

Wiesz, jak uruchomić pipeline nf-core używając jego wbudowanego profilu testowego i gdzie znaleźć jego wyjścia.

### Co dalej?

Dowiedz się, jak kod pipeline'a jest zorganizowany.

---

## 3. Zbadaj strukturę kodu pipeline'a

Teraz, gdy pomyślnie uruchomiliśmy pipeline jako użytkownicy, zmieńmy naszą perspektywę, aby przyjrzeć się wewnętrznej organizacji pipeline'ów nf-core.

Projekt nf-core egzekwuje silne wytyczne dotyczące tego, jak powinny być strukturyzowane pipeline'y oraz jak powinien być zorganizowany, skonfigurowany i udokumentowany kod.
Zrozumienie tego, jak to wszystko jest zorganizowane, jest pierwszym krokiem w kierunku tworzenia własnych pipeline'ów kompatybilnych z nf-core, co podejmiemy w Części 2 tego kursu.

Spójrzmy, jak kod pipeline'a jest zorganizowany w repozytorium `nf-core/demo`, używając dowiązania symbolicznego `pipelines`, które utworzyliśmy wcześniej.

Możesz użyć `tree` lub użyć eksploratora plików, aby znaleźć i otworzyć katalog `nf-core/demo`.

```bash
tree -L 1 pipelines/nf-core/demo
```

??? abstract "Zawartość katalogu"

    ```console
    pipelines/nf-core/demo
    ├── assets
    ├── CHANGELOG.md
    ├── CITATIONS.md
    ├── CODE_OF_CONDUCT.md
    ├── conf
    ├── docs
    ├── LICENSE
    ├── main.nf
    ├── modules
    ├── modules.json
    ├── nextflow.config
    ├── nextflow_schema.json
    ├── nf-test.config
    ├── README.md
    ├── ro-crate-metadata.json
    ├── subworkflows
    ├── tests
    ├── tower.yml
    └── workflows
    ```

Dzieje się tam dużo, więc zajmiemy się tym krok po kroku.

Po pierwsze, zauważmy, że na najwyższym poziomie możesz znaleźć plik README z informacjami podsumowującymi, a także pliki akcesoriów, które podsumowują informacje o projekcie, takie jak licencjonowanie, wytyczne dotyczące wkładu, cytowania i kodeks postępowania.
Szczegółowa dokumentacja pipeline'a znajduje się w katalogu `docs`.
Cała ta zawartość jest używana do programowego generowania stron internetowych na stronie nf-core, więc są one zawsze aktualne z kodem.

Teraz, co do reszty, podzielimy naszą eksplorację na trzy etapy:

1. Komponenty kodu pipeline'a (`main.nf`, `workflows`, `subworkflows`, `modules`)
2. Konfiguracja pipeline'a
3. Wejścia i walidacja

Zacznijmy od komponentów kodu pipeline'a.
Skoncentrujemy się na hierarchii plików i organizacji strukturalnej, zamiast zagłębiać się w kod w poszczególnych plikach.

### 3.1. Komponenty kodu pipeline'a

Standardowa organizacja kodu pipeline'a nf-core podąża za modularną strukturą zaprojektowaną tak, aby maksymalizować ponowne użycie kodu, jak wprowadzono w [Hello Modules](../hello_nextflow/04_hello_modules.md), Części 4 kursu [Hello Nextflow](../hello_nextflow/index.md), chociaż w prawdziwie nf-core'owym stylu jest to zaimplementowane z odrobiną dodatkowej złożoności.
Konkretnie, pipeline'y nf-core obficie wykorzystują subworkflow'y, tj. skrypty importowane przez nadrzędny workflow.

To może brzmieć trochę abstrakcyjnie, więc spójrzmy, jak jest to używane w praktyce w pipeline'ie `nf-core/demo`.

!!! note

    Nie przejdziemy przez faktyczny kod opisujący _sposób_, w jaki te komponenty modułowe są połączone, ponieważ istnieje pewna dodatkowa złożoność związana z użyciem subworkflow'ów, która może być myląca. Zrozumienie tego nie jest konieczne na tym etapie szkolenia.
    Na razie skoncentrujemy się na ogólnej organizacji i logice.

#### 3.1.1. Ogólny przegląd

Oto jak wyglądają relacje między odpowiednimi komponentami kodu dla pipeline'a `nf-core/demo`:

<figure class="excalidraw">
    --8<-- "docs/en/docs/hello_nf-core/img/nf-core_demo_code_organization.svg"
</figure>

Istnieje tak zwany skrypt _punktu wejścia_ o nazwie `main.nf`, który działa jako wrapper dla dwóch rodzajów zagnieżdżonych workflow'ów: workflow'u zawierającego rzeczywistą logikę analizy, zlokalizowanego w `workflows/` i nazwanego `demo.nf`, oraz zestawu workflow'ów pomocniczych zlokalizowanych w `subworkflows/`.
Workflow `demo.nf` wywołuje **moduły** zlokalizowane w `modules/`; zawierają one **procesy**, które wykonają rzeczywiste kroki analizy.

!!! note

    Subworkflow'y nie są ograniczone do funkcji pomocniczych i mogą wykorzystywać moduły procesów.

    Pipeline `nf-core/demo` pokazany tutaj jest po prostszej stronie spektrum, ale inne pipeline'y nf-core (takie jak `nf-core/rnaseq`) wykorzystują subworkflow'y zaangażowane w rzeczywistą analizę.

Teraz przejrzyjmy te komponenty po kolei.

#### 3.1.2. Skrypt punktu wejścia: `main.nf`

Skrypt `main.nf` jest punktem wejścia, od którego Nextflow rozpoczyna, gdy wykonujemy `nextflow run nf-core/demo`.
Oznacza to, że gdy uruchamiasz `nextflow run nf-core/demo` aby uruchomić pipeline, Nextflow automatycznie znajduje i wykonuje skrypt `main.nf`.
Działa to dla każdego pipeline'a Nextflow, który podąża za tą konwencjonalną nazwą i strukturą, nie tylko dla pipeline'ów nf-core.

Użycie skryptu punktu wejścia ułatwia uruchamianie ustandaryzowanych 'pomocniczych' subworkflow'ów przed i po uruchomieniu właściwego skryptu analizy.
Przejdziemy przez nie po tym, jak przejrzymy rzeczywisty workflow analizy i jego moduły.

#### 3.1.3. Skrypt analizy: `workflows/demo.nf`

Workflow `workflows/demo.nf` to miejsce, w którym przechowywana jest centralna logika pipeline'a.
Jest on ustrukturyzowany podobnie jak normalny workflow Nextflow, z wyjątkiem tego, że jest zaprojektowany tak, aby być wywoływanym z workflow'u nadrzędnego, co wymaga kilku dodatkowych funkcji.
Omówimy odpowiednie różnice w następnej części tego kursu, gdy zajmiemy się konwersją prostego pipeline'a Hello z Hello Nextflow do formy kompatybilnej z nf-core.

Workflow `demo.nf` wywołuje **moduły** zlokalizowane w `modules/`, które przejrzymy w następnej kolejności.

!!! note

    Niektóre workflow'y analizy nf-core wyświetlają dodatkowe poziomy zagnieżdżenia poprzez wywoływanie subworkflow'ów niższego poziomu.
    Jest to głównie używane do opakowywania dwóch lub więcej modułów, które są powszechnie używane razem, w łatwe do ponownego użycia segmenty pipeline'a.
    Możesz zobaczyć kilka przykładów, przeglądając dostępne [subworkflow'y nf-core](https://nf-co.re/subworkflows/) na stronie nf-core.

    Gdy skrypt analizy używa subworkflow'ów, są one przechowywane w katalogu `subworkflows/`.

#### 3.1.4. Moduły

Moduły to miejsce, w którym znajduje się kod procesu, jak opisano w [Części 4 kursu szkoleniowego Hello Nextflow](../hello_nextflow/04_hello_modules.md).

W projekcie nf-core moduły są organizowane przy użyciu wielopoziomowej zagnieżdżonej struktury, która odzwierciedla zarówno ich pochodzenie, jak i ich zawartość.
Na najwyższym poziomie moduły są różnicowane jako `nf-core` lub `local` (nie będące częścią projektu nf-core), a następnie dalej umieszczane w katalogu nazwanym według narzędzia(i), które opakowują.
Jeśli narzędzie należy do zestawu narzędzi (tj. pakietu zawierającego wiele narzędzi), istnieje pośredni poziom katalogu nazwany według zestawu narzędzi.

Możesz zobaczyć to zastosowane w praktyce do modułów pipeline'a `nf-core/demo`:

```bash
tree -L 3 pipelines/nf-core/demo/modules
```

??? abstract "Zawartość katalogu"

    ```console
    pipelines/nf-core/demo/modules
    └── nf-core
        ├── fastqc
        │   ├── environment.yml
        │   ├── main.nf
        │   ├── meta.yml
        │   └── tests
        ├── multiqc
        │   ├── environment.yml
        │   ├── main.nf
        │   ├── meta.yml
        │   └── tests
        └── seqtk
            └── trim

    7 directories, 6 files
    ```

Tutaj widzisz, że moduły `fastqc` i `multiqc` znajdują się na najwyższym poziomie w modułach `nf-core`, podczas gdy moduł `trim` znajduje się pod zestawem narzędzi, do którego należy, `seqtk`.
W tym przypadku nie ma żadnych modułów `local`.

Plik kodu modułu opisujący proces zawsze nazywa się `main.nf` i jest accompanied by testami i plikami `.yml`, które na razie zignorujem.

Razem wzięte, workflow punktu wejścia, workflow analizy i moduły są wystarczające do uruchomienia 'interesujących' części pipeline'a.
Jednak wiemy, że są tam również subworkflow'y pomocnicze, więc spójrzmy na nie teraz.

#### 3.1.5. Subworkflows pomocnicze

Podobnie jak moduły, subworkflow'y są różnicowane na katalogi `local` i `nf-core`, a każdy subworkflow ma swoją własną zagnieżdżoną strukturę katalogów ze swoim własnym skryptem `main.nf`, testami i plikiem `.yml`.

```bash
tree -L 3 pipelines/nf-core/demo/subworkflows
```

??? abstract "Zawartość katalogu"

    ```console
    pipelines/nf-core/demo/subworkflows
    ├── local
    │   └── utils_nfcore_demo_pipeline
    │       └── main.nf
    └── nf-core
        ├── utils_nextflow_pipeline
        │   ├── main.nf
        │   ├── meta.yml
        │   └── tests
        ├── utils_nfcore_pipeline
        │   ├── main.nf
        │   ├── meta.yml
        │   └── tests
        └── utils_nfschema_plugin
            ├── main.nf
            ├── meta.yml
            └── tests

    9 directories, 7 files
    ```

Jak zauważono powyżej, pipeline `nf-core/demo` nie zawiera żadnych subworkflow'ów specyficznych dla analizy, więc wszystkie subworkflow'y, które tutaj widzimy, są tak zwanymi workflow'ami 'pomocniczymi' lub 'użytkowymi', jak wskazuje prefiks `utils_` w ich nazwach.
Te subworkflow'y to te, które produkują wymyślny nagłówek nf-core w wyjściu konsoli, między innymi funkcjami akcesoriów.

!!! tip

    Poza ich wzorcem nazewnictwa, inną wskazówką, że te subworkflow'y nie wykonują żadnej naprawdę związanej z analizą funkcji, jest to, że nie wywołują żadnych procesów w ogóle.

To kończy zestawienie podstawowych komponentów kodu, które stanowią pipeline `nf-core/demo`.
Teraz spójrzmy na pozostałe elementy, o których powinieneś wiedzieć trochę przed zagłębieniem się w rozwój: konfigurację pipeline'a i walidację wejścia.

### 3.2. Konfiguracja pipeline'a

Nauczyłeś się wcześniej, że Nextflow oferuje wiele opcji konfigurowania wykonania pipeline'a, czy to w zakresie wejść i parametrów, zasobów obliczeniowych i innych aspektów orkiestracji.
Projekt nf-core stosuje wysoce ustandaryzowane wytyczne dla konfiguracji pipeline'ów, które mają na celu budowanie na elastycznych możliwościach dostosowywania Nextflow w sposób zapewniający większą spójność i łatwość konserwacji między pipeline'ami.

Centralny plik konfiguracyjny `nextflow.config` jest używany do ustawiania domyślnych wartości dla parametrów i innych opcji konfiguracyjnych.
Większość tych opcji konfiguracyjnych jest stosowana domyślnie, podczas gdy inne (np. profile zależności oprogramowania) są włączone jako opcjonalne profile.

Istnieje kilka dodatkowych plików konfiguracyjnych, które są przechowywane w katalogu `conf` i które mogą być dodane do konfiguracji domyślnie lub opcjonalnie jako profile:

- `base.config`: Plik konfiguracyjny 'czystej karty', odpowiedni do ogólnego użytku w większości środowisk obliczeniowych o wysokiej wydajności. Definiuje szerokie przedziały użycia zasobów, na przykład, które są wygodne do zastosowania do modułów.
- `modules.config`: Dodatkowe dyrektywy i argumenty modułów.
- `test.config`: Profil do uruchomienia pipeline'a z minimalnymi danymi testowymi, który użyliśmy, gdy uruchomiliśmy pipeline demo.
- `test_full.config`: Profil do uruchomienia pipeline'a z pełnowymiarowym zestawem danych testowych.

Dotniemy kilku z tych plików później w kursie.

### 3.3. Wejścia i walidacja

Jak zauważyliśmy wcześniej, gdy badaliśmy profil testowy pipeline'a `nf-core/demo`, jest on zaprojektowany tak, aby przyjmować jako wejście samplesheet zawierający ścieżki plików i identyfikatory próbek.
Ścieżki plików były powiązane z rzeczywistymi danymi zlokalizowanymi w repozytorium `nf-core/test-datasets`.

Przykładowy samplesheet jest również dostarczany w katalogu `assets`, chociaż ścieżki w tym nie są rzeczywiste.

```csv title="assets/samplesheet.csv" linenums="1"
sample,fastq_1,fastq_2
SAMPLE_PAIRED_END,/path/to/fastq/files/AEG588A1_S1_L002_R1_001.fastq.gz,/path/to/fastq/files/AEG588A1_S1_L002_R2_001.fastq.gz
SAMPLE_SINGLE_END,/path/to/fastq/files/AEG588A4_S4_L003_R1_001.fastq.gz,

```

Ten konkretny samplesheet jest dość prosty, ale niektóre pipeline'y działają na samplesheets, które są bardziej złożone, z o wiele większą ilością metadanych związanych z podstawowymi wejściami.

Niestety, ponieważ te pliki mogą być trudne do sprawdzenia gołym okiem, nieprawidłowe formatowanie danych wejściowych jest bardzo powszechnym źródłem awarii pipeline'a.
Powiązanym problemem jest sytuacja, gdy parametry są nieprawidłowo podane.

Rozwiązaniem tych problemów jest uruchomienie automatycznych kontroli walidacyjnych na wszystkich plikach wejściowych, aby upewnić się, że zawierają oczekiwane typy informacji, poprawnie sformatowane, oraz na parametrach, aby upewnić się, że są oczekiwanego typu.
Nazywa się to walidacją wejścia i idealnie powinno być wykonane _przed_ próbą uruchomienia pipeline'a, zamiast czekać, aż pipeline się nie powiedzie, aby dowiedzieć się, że był problem z wejściami.

Podobnie jak w przypadku konfiguracji, projekt nf-core jest bardzo zopiniowany na temat walidacji wejścia i zaleca użycie [wtyczki nf-schema](https://nextflow-io.github.io/nf-schema/latest/), wtyczki Nextflow, która zapewnia kompleksowe możliwości walidacji dla pipeline'ów Nextflow.

Omówimy ten temat bardziej szczegółowo w Części 5 tego kursu.
Na razie po prostu bądź świadomy, że istnieją dwa pliki JSON dostarczone w tym celu: `nextflow_schema.json` i `assets/schema_input.json`.

`nextflow_schema.json` to plik używany do przechowywania informacji o parametrach pipeline'a, w tym typ, opis i tekst pomocy w formacie czytelnym maszynowo.
Jest to używane do różnych celów, w tym automatycznej walidacji parametrów, generowania tekstu pomocy i interaktywnego renderowania formularza parametrów w interfejsach użytkownika.

`schema_input.json` to plik używany do definiowania struktury samplesheet wejściowego.
Każda kolumna może mieć typ, wzorzec, opis i tekst pomocy w formacie czytelnym maszynowo.
Schemat jest używany do różnych celów, w tym automatycznej walidacji i dostarczania pomocnych komunikatów o błędach.

### Podsumowanie

Wiesz, jakie są główne komponenty pipeline'a nf-core i jak kod jest zorganizowany; gdzie znajdują się główne elementy konfiguracji; i jesteś świadomy, do czego służy walidacja wejścia.

### Co dalej?

Zrób sobie przerwę! To było dużo. Gdy poczujesz się odświeżony i gotowy, przejdź do następnej sekcji, aby zastosować to, czego się nauczyłeś, do napisania pipeline'a kompatybilnego z nf-core.

!!! tip

    Jeśli chciałbyś dowiedzieć się, jak komponować workflow z subworkflow'ów przed przejściem do następnej części, sprawdź [Workflows of Workflows](../side_quests/workflows_of_workflows.md) Side Quest.
