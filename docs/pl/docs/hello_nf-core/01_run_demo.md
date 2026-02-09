# Część 1: Uruchom demonstracyjny pipeline

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Tłumaczenie wspomagane przez AI - [dowiedz się więcej i zaproponuj poprawki](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

W tej pierwszej części szkolenia Hello nf-core pokażemy Ci, jak znaleźć i wypróbować pipeline nf-core, zrozumieć organizację kodu oraz rozpoznać, czym różni się on od zwykłego kodu Nextflow'a przedstawionego w [Hello Nextflow](../hello_nextflow/index.md).

Użyjemy pipeline'u o nazwie nf-core/demo, który jest utrzymywany przez projekt nf-core jako część jego kolekcji pipeline'ów służących do demonstracji struktury kodu i działania narzędzi.

Upewnij się, że Twój katalog roboczy jest ustawiony na `hello-nf-core/`, zgodnie z instrukcjami na stronie [Pierwsze kroki](./00_orientation.md).

---

## 1. Znajdź i pobierz pipeline nf-core/demo

Zacznijmy od zlokalizowania pipeline'u nf-core/demo na stronie projektu [nf-co.re](https://nf-co.re), która centralizuje wszystkie informacje, takie jak: ogólna dokumentacja i artykuły pomocnicze, dokumentacja dla każdego z pipeline'ów, wpisy na blogu, ogłoszenia o wydarzeniach i tak dalej.

### 1.1. Znajdź pipeline na stronie internetowej

W przeglądarce internetowej przejdź do [https://nf-co.re/pipelines/](https://nf-co.re/pipelines/) i wpisz `demo` w pasku wyszukiwania.

![wyniki wyszukiwania](./img/search-results.png)

Kliknij nazwę pipeline'u, `demo`, aby uzyskać dostęp do strony dokumentacji pipeline'u.

Każdy wydany pipeline ma dedykowaną stronę zawierającą następujące sekcje dokumentacji:

- **Introduction:** Wprowadzenie i przegląd pipeline'u
- **Usage:** Opisy sposobu wykonywania pipeline'u
- **Parameters:** Zgrupowane parametry pipeline'u z opisami
- **Output:** Opisy i przykłady oczekiwanych plików wyjściowych
- **Results:** Przykładowe pliki wyjściowe wygenerowane z pełnego zestawu danych testowych
- **Releases & Statistics:** Historia wersji pipeline'u i statystyki

Zawsze, gdy rozważasz przyjęcie nowego pipeline'u, powinieneś najpierw uważnie przeczytać dokumentację pipeline'u, aby zrozumieć, co robi i jak powinien być skonfigurowany, zanim spróbujesz go uruchomić.

Spójrz teraz i sprawdź, czy możesz dowiedzieć się:

- Które narzędzia uruchomi pipeline (Sprawdź zakładkę: `Introduction`)
- Jakie dane wejściowe i parametry akceptuje lub wymaga pipeline (Sprawdź zakładkę: `Parameters`)
- Jakie są wyjścia produkowane przez pipeline (Sprawdź zakładkę: `Output`)

#### 1.1.1. Przegląd pipeline'u

Zakładka `Introduction` zawiera przegląd pipeline'u, w tym wizualną reprezentację (nazywaną mapą metra) oraz listę narzędzi uruchamianych w ramach pipeline'u.

![mapa metra pipeline'u](./img/nf-core-demo-subway-cropped.png)

1. Read QC (FASTQC)
2. Adapter and quality trimming (SEQTK_TRIM)
3. Present QC for raw reads (MULTIQC)

#### 1.1.2. Przykładowa linia poleceń

Dokumentacja zawiera również przykładowy plik wejściowy (omówiony bardziej szczegółowo poniżej) oraz przykładową linię poleceń.

```bash
nextflow run nf-core/demo \
  -profile <docker/singularity/.../institute> \
  --input samplesheet.csv \
  --outdir <OUTDIR>
```

Zauważysz, że przykładowe polecenie NIE określa pliku workflow'a, tylko odwołanie do repozytorium pipeline'u, `nf-core/demo`.

Wywołany w ten sposób, Nextflow założy, że kod jest zorganizowany w określony sposób.
Pobierzmy kod, abyśmy mogli zbadać tę strukturę.

### 1.2. Pobierz kod pipeline'u

Gdy już ustalimy, że pipeline wydaje się odpowiedni do naszych celów, wypróbujmy go.
Na szczęście Nextflow ułatwia pobieranie pipeline'ów z prawidłowo sformatowanych repozytoriów bez konieczności ręcznego pobierania czegokolwiek.

Wróćmy do terminala i uruchommy następujące polecenie:

```bash
nextflow pull nf-core/demo
```

??? success "Wyjście polecenia"

    ```console
    Checking nf-core/demo ...
    downloaded from https://github.com/nf-core/demo.git - revision: 04060b4644 [master]
    ```

Nextflow wykonuje `pull` kodu pipeline'u, co oznacza, że pobiera pełne repozytorium na Twój dysk lokalny.

Dla jasności, możesz to zrobić z dowolnym pipeline'em Nextflow'a, który jest odpowiednio skonfigurowany w GitHub, nie tylko z pipeline'ami nf-core.
Jednak nf-core jest największą otwartą kolekcją pipeline'ów Nextflow'a.

Możesz poprosić Nextflow'a o podanie listy pipeline'ów, które pobrałeś w ten sposób:

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

!!! note "Uwaga"

    Pełna ścieżka może się różnić w Twoim systemie, jeśli nie korzystasz z naszego środowiska szkoleniowego.

Nextflow celowo przechowuje pobrany kod źródłowy 'z dala od oczu', kierując się zasadą, że te pipeline'y powinny być używane bardziej jak biblioteki niż kod, z którym bezpośrednio wchodziłbyś w interakcję.

Jednak dla celów tego szkolenia chcemy móc zajrzeć do środka i zobaczyć, co tam jest.
Aby to ułatwić, stwórzmy dowiązanie symboliczne do tej lokalizacji z naszego bieżącego katalogu roboczego.

```bash
ln -s $NXF_HOME/assets pipelines
```

To tworzy skrót, który ułatwia eksplorację pobranego kodu.

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

Ale najpierw uruchommy nasz pierwszy pipeline nf-core!

### Podsumowanie

Teraz wiesz, jak znaleźć pipeline za pośrednictwem strony nf-core i pobrać lokalną kopię kodu źródłowego.

### Co dalej?

Dowiedz się, jak wypróbować pipeline nf-core przy minimalnym wysiłku.

---

## 2. Wypróbuj pipeline z jego profilem testowym

Wygodnie, każdy pipeline nf-core jest dostarczany z profilem testowym.
Jest to minimalny zestaw ustawień konfiguracyjnych dla pipeline'u do uruchomienia przy użyciu małego zestawu danych testowych hostowanego w repozytorium [nf-core/test-datasets](https://github.com/nf-core/test-datasets).
To świetny sposób na szybkie wypróbowanie pipeline'u na małą skalę.

!!! note "Uwaga"

    System profili konfiguracyjnych Nextflow'a pozwala łatwo przełączać się między różnymi silnikami kontenerów lub środowiskami wykonawczymi.
    Więcej szczegółów znajdziesz w [Hello Nextflow Część 6: Konfiguracja](../hello_nextflow/06_hello_config.md).

### 2.1. Zbadaj profil testowy

Dobrą praktyką jest sprawdzenie, co określa profil testowy pipeline'u przed jego uruchomieniem.
Profil `test` dla `nf-core/demo` znajduje się w pliku konfiguracyjnym `conf/test.config` i jest pokazany poniżej.

```groovy title="conf/test.config" linenums="1" hl_lines="8 26"
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Nextflow config file for running minimal tests
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Defines input files and everything required to run a fast and simple pipeline test.

    Use as follows:
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

    // Input data
    input  = 'https://raw.githubusercontent.com/nf-core/test-datasets/viralrecon/samplesheet/samplesheet_test_illumina_amplicon.csv'

}
```

Od razu zauważysz, że blok komentarza na górze zawiera przykład użycia pokazujący, jak uruchomić pipeline z tym profilem testowym.

```groovy title="conf/test.config" linenums="7"
Use as follows:
        nextflow run nf-core/demo -profile test,<docker/singularity> --outdir <OUTDIR>
```

Jedyne rzeczy, które musimy podać, to to, co jest pokazane między nawiasami ostrymi w przykładowym poleceniu: `<docker/singularity>` i `<OUTDIR>`.

Przypominając, `<docker/singularity>` odnosi się do wyboru systemu kontenerów. Wszystkie pipeline'y nf-core są zaprojektowane tak, aby można było ich używać z kontenerami (Docker, Singularity itp.) w celu zapewnienia powtarzalności i wyeliminowania problemów z instalacją oprogramowania.
Musimy więc określić, czy chcemy użyć Dockera czy Singularity do przetestowania pipeline'u.

Część `--outdir <OUTDIR>` odnosi się do katalogu, w którym Nextflow zapisze wyjścia pipeline'u.
Musimy podać dla niego nazwę, którą możemy po prostu wymyślić.
Jeśli jeszcze nie istnieje, Nextflow utworzy go dla nas w czasie wykonywania.

Przechodząc do sekcji po bloku komentarza, profil testowy pokazuje nam, co zostało wstępnie skonfigurowane do testowania: przede wszystkim parametr `input` jest już ustawiony tak, aby wskazywał na zestaw danych testowych, więc nie musimy dostarczać własnych danych.
Jeśli podążysz za linkiem do wstępnie skonfigurowanego wejścia, zobaczysz, że jest to plik csv zawierający identyfikatory próbek i ścieżki plików dla kilku próbek eksperymentalnych.

```csv title="samplesheet_test_illumina_amplicon.csv"
sample,fastq_1,fastq_2
SAMPLE1_PE,https://raw.githubusercontent.com/nf-core/test-datasets/viralrecon/illumina/amplicon/sample1_R1.fastq.gz,https://raw.githubusercontent.com/nf-core/test-datasets/viralrecon/illumina/amplicon/sample1_R2.fastq.gz
SAMPLE2_PE,https://raw.githubusercontent.com/nf-core/test-datasets/viralrecon/illumina/amplicon/sample2_R1.fastq.gz,https://raw.githubusercontent.com/nf-core/test-datasets/viralrecon/illumina/amplicon/sample2_R2.fastq.gz
SAMPLE3_SE,https://raw.githubusercontent.com/nf-core/test-datasets/viralrecon/illumina/amplicon/sample1_R1.fastq.gz,
SAMPLE3_SE,https://raw.githubusercontent.com/nf-core/test-datasets/viralrecon/illumina/amplicon/sample2_R1.fastq.gz,
```

Nazywa się to arkuszem próbek (samplesheet) i jest to najczęstsza forma wejścia do pipeline'ów nf-core.

!!! note "Uwaga"

    Nie martw się, jeśli nie znasz formatów i typów danych, nie jest to ważne dla tego, co następuje.

To potwierdza, że mamy wszystko, czego potrzebujemy, aby wypróbować pipeline.

### 2.2. Uruchom pipeline

Zdecydujmy się użyć Dockera jako systemu kontenerów i `demo-results` jako katalogu wyjściowego, i jesteśmy gotowi do uruchomienia polecenia testowego:

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

Jeśli Twoje wyjście pasuje do tego, gratulacje! Właśnie uruchomiłeś swój pierwszy pipeline nf-core.

Zauważysz, że jest znacznie więcej wyjścia w konsoli niż podczas uruchamiania podstawowego pipeline'u Nextflow'a.
Jest nagłówek zawierający podsumowanie wersji pipeline'u, danych wejściowych i wyjściowych oraz kilku elementów konfiguracji.

!!! note "Uwaga"

    Twoje wyjście pokaże różne znaczniki czasu, nazwy wykonań i ścieżki plików, ale ogólna struktura i wykonanie procesów powinny być podobne.

Przechodząc do wyjścia wykonania, spójrzmy na linie, które mówią nam, jakie procesy zostały uruchomione:

```console
    [ff/a6976b] NFCORE_DEMO:DEMO:FASTQC (SAMPLE3_SE)     | 3 of 3 ✔
    [39/731ab7] NFCORE_DEMO:DEMO:SEQTK_TRIM (SAMPLE3_SE) | 3 of 3 ✔
    [7c/78d96e] NFCORE_DEMO:DEMO:MULTIQC                 | 1 of 1 ✔
```

To mówi nam, że uruchomiono trzy procesy, odpowiadające trzem narzędziom pokazanym na stronie dokumentacji pipeline'u na stronie nf-core: FASTQC, SEQTK_TRIM i MULTIQC.

Pełne nazwy procesów, jak pokazano tutaj, takie jak `NFCORE_DEMO:DEMO:MULTIQC`, są dłuższe niż to, co mogłeś zobaczyć we wprowadzającym materiale Hello Nextflow.
Zawierają one nazwy ich nadrzędnych workflow'ów i odzwierciedlają modularność kodu pipeline'u.
Wejdziemy w więcej szczegółów na ten temat za chwilę.

### 2.3. Zbadaj wyjścia pipeline'u

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
Aby dowiedzieć się więcej o wyjściach pipeline'u `nf-core/demo`, sprawdź jego [stronę dokumentacji](https://nf-co.re/demo/1.0.2/docs/output/).

Na tym etapie ważne jest zaobserwowanie, że wyniki są zorganizowane według modułów, a dodatkowo istnieje katalog o nazwie `pipeline_info` zawierający różne raporty z oznaczeniami czasu dotyczące wykonania pipeline'u.

Na przykład plik `execution_timeline_*` pokazuje, jakie procesy zostały uruchomione, w jakiej kolejności i jak długo trwało ich uruchomienie:

![raport osi czasu wykonania](./img/execution_timeline.png)

!!! note "Uwaga"

    Tutaj zadania nie były uruchamiane równolegle, ponieważ działamy na minimalistycznej maszynie w Github Codespaces.
    Aby zobaczyć ich równoległe uruchomienie, spróbuj zwiększyć alokację CPU Twojego codespace'a i limity zasobów w konfiguracji testowej.

Te raporty są generowane automatycznie dla wszystkich pipeline'ów nf-core.

### Podsumowanie

Wiesz, jak uruchomić pipeline nf-core używając jego wbudowanego profilu testowego i gdzie znaleźć jego wyjścia.

### Co dalej?

Dowiedz się, jak zorganizowany jest kod pipeline'u.

---

## 3. Zbadaj strukturę kodu pipeline'u

Teraz, gdy pomyślnie uruchomiliśmy pipeline jako użytkownicy, zmieńmy perspektywę, aby przyjrzeć się, jak pipeline'y nf-core są zorganizowane wewnętrznie.

Projekt nf-core wymusza silne wytyczne dotyczące struktury pipeline'ów oraz sposobu organizacji, konfiguracji i dokumentowania kodu.
Zrozumienie, jak to wszystko jest zorganizowane, jest pierwszym krokiem w kierunku tworzenia własnych pipeline'ów kompatybilnych z nf-core, co podejmiemy w Części 2 tego kursu.

Spójrzmy, jak kod pipeline'u jest zorganizowany w repozytorium `nf-core/demo`, używając dowiązania symbolicznego `pipelines`, które utworzyliśmy wcześniej.

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

Dzieje się tam wiele, więc podejmiemy to krok po kroku.

Po pierwsze, zauważmy, że na najwyższym poziomie możesz znaleźć plik README z informacjami podsumowującymi, a także pliki pomocnicze, które podsumowują informacje o projekcie, takie jak licencjonowanie, wytyczne dotyczące wkładu, cytowanie i kodeks postępowania.
Szczegółowa dokumentacja pipeline'u znajduje się w katalogu `docs`.
Cała ta zawartość jest używana do programowego generowania stron internetowych na stronie nf-core, więc są zawsze aktualne z kodem.

Teraz, dla reszty, podzielimy naszą eksplorację na trzy etapy:

1. Komponenty kodu pipeline'u (`main.nf`, `workflows`, `subworkflows`, `modules`)
2. Konfiguracja pipeline'u
3. Dane wejściowe i walidacja

Zacznijmy od komponentów kodu pipeline'u.
Skupimy się na hierarchii plików i organizacji strukturalnej, zamiast zagłębiać się w kod w poszczególnych plikach.

### 3.1. Komponenty kodu pipeline'u

Standardowa organizacja kodu pipeline'u nf-core podąża za modularną strukturą, która jest zaprojektowana tak, aby maksymalizować ponowne wykorzystanie kodu, jak wprowadzono w [Hello Modules](../hello_nextflow/04_hello_modules.md), Część 4 kursu [Hello Nextflow](../hello_nextflow/index.md), chociaż w prawdziwym stylu nf-core jest to zaimplementowane z odrobiną dodatkowej złożoności.
W szczególności pipeline'y nf-core obficie wykorzystują subworkflow'y, tj. skrypty workflow'ów importowane przez nadrzędny workflow.

To może brzmieć nieco abstrakcyjnie, więc spójrzmy, jak jest to używane w praktyce w pipeline'ie `nf-core/demo`.

!!! note "Uwaga"

    Nie będziemy omawiać rzeczywistego kodu dotyczącego _sposobu_ łączenia tych modularnych komponentów, ponieważ istnieje pewna dodatkowa złożoność związana z użyciem subworkflow'ów, która może być myląca, a zrozumienie tego nie jest konieczne na tym etapie szkolenia.
    Na razie skupimy się na ogólnej organizacji i logice.

#### 3.1.1. Ogólny przegląd

Oto jak wyglądają relacje między odpowiednimi komponentami kodu dla pipeline'u `nf-core/demo`:

<figure class="excalidraw">
    --8<-- "docs/en/docs/hello_nf-core/img/nf-core_demo_code_organization.svg"
</figure>

Istnieje tak zwany skrypt _punktu wejścia_ o nazwie `main.nf`, który działa jako opakowanie dla dwóch rodzajów zagnieżdżonych workflow'ów: workflow'a zawierającego rzeczywistą logikę analizy, znajdującego się w `workflows/` i nazwanego `demo.nf`, oraz zestawu workflow'ów pomocniczych znajdujących się w `subworkflows/`.
Workflow `demo.nf` wywołuje **moduły** znajdujące się w `modules/`; zawierają one **procesy**, które będą wykonywać rzeczywiste kroki analizy.

!!! note "Uwaga"

    Subworkflow'y nie są ograniczone do funkcji pomocniczych i mogą wykorzystywać moduły procesów.

    Pipeline `nf-core/demo` pokazany tutaj jest po prostszej stronie spektrum, ale inne pipeline'y nf-core (takie jak `nf-core/rnaseq`) wykorzystują subworkflow'y zaangażowane w rzeczywistą analizę.

Teraz przejrzyjmy te komponenty po kolei.

#### 3.1.2. Skrypt punktu wejścia: `main.nf`

Skrypt `main.nf` jest punktem wejścia, od którego Nextflow zaczyna, gdy wykonujemy `nextflow run nf-core/demo`.
Oznacza to, że gdy uruchamiasz `nextflow run nf-core/demo`, aby uruchomić pipeline, Nextflow automatycznie znajduje i wykonuje skrypt `main.nf`.
Działa to dla każdego pipeline'u Nextflow'a, który podąża za tą konwencjonalną nazwą i strukturą, nie tylko dla pipeline'ów nf-core.

Użycie skryptu punktu wejścia ułatwia uruchamianie standardowych subworkflow'ów 'pomocniczych' przed i po uruchomieniu rzeczywistego skryptu analizy.
Przejdziemy przez nie po przejrzeniu rzeczywistego workflow'a analizy i jego modułów.

#### 3.1.3. Skrypt analizy: `workflows/demo.nf`

Workflow `workflows/demo.nf` to miejsce, w którym przechowywana jest centralna logika pipeline'u.
Jest zbudowany podobnie jak normalny workflow Nextflow'a, z wyjątkiem tego, że jest zaprojektowany tak, aby być wywoływanym z nadrzędnego workflow'a, co wymaga kilku dodatkowych funkcji.
Omówimy odpowiednie różnice w następnej części tego kursu, gdy podejmiemy konwersję prostego pipeline'u Hello z Hello Nextflow do formy kompatybilnej z nf-core.

Workflow `demo.nf` wywołuje **moduły** znajdujące się w `modules/`, które przejrzymy dalej.

!!! note "Uwaga"

    Niektóre workflow'y analizy nf-core wyświetlają dodatkowe poziomy zagnieżdżenia, wywołując subworkflow'y niższego poziomu.
    Jest to głównie używane do opakowywania dwóch lub więcej modułów, które są powszechnie używane razem, w łatwo wielokrotnie używane segmenty pipeline'u.
    Możesz zobaczyć kilka przykładów, przeglądając dostępne [subworkflow'y nf-core](https://nf-co.re/subworkflows/) na stronie nf-core.

    Gdy skrypt analizy używa subworkflow'ów, są one przechowywane w katalogu `subworkflows/`.

#### 3.1.4. Moduły

Moduły to miejsce, w którym znajduje się kod procesu, jak opisano w [Części 4 kursu szkoleniowego Hello Nextflow](../hello_nextflow/04_hello_modules.md).

W projekcie nf-core moduły są organizowane przy użyciu wielopoziomowej zagnieżdżonej struktury, która odzwierciedla zarówno ich pochodzenie, jak i zawartość.
Na najwyższym poziomie moduły są rozróżniane jako `nf-core` lub `local` (nie będące częścią projektu nf-core), a następnie dalej umieszczane w katalogu nazwanym po narzędziu (narzędziach), które opakowują.
Jeśli narzędzie należy do zestawu narzędzi (tj. pakietu zawierającego wiele narzędzi), istnieje pośredni poziom katalogu nazwany po zestawie narzędzi.

Możesz zobaczyć to zastosowane w praktyce do modułów pipeline'u `nf-core/demo`:

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
W tym przypadku nie ma modułów `local`.

Plik kodu modułu opisujący proces zawsze nazywa się `main.nf` i towarzyszy mu testy oraz pliki `.yml`, które na razie zignorujemy.

Razem wzięte, workflow punktu wejścia, workflow analizy i moduły są wystarczające do uruchomienia 'interesujących' części pipeline'u.
Jednak wiemy, że są tam również subworkflow'y pomocnicze, więc spójrzmy na nie teraz.

#### 3.1.5. Subworkflow'y pomocnicze

Podobnie jak moduły, subworkflow'y są rozróżniane na katalogi `local` i `nf-core`, a każdy subworkflow ma własną zagnieżdżoną strukturę katalogów z własnym skryptem `main.nf`, testami i plikiem `.yml`.

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

Jak wspomniano powyżej, pipeline `nf-core/demo` nie zawiera żadnych subworkflow'ów specyficznych dla analizy, więc wszystkie subworkflow'y, które tu widzimy, to tak zwane workflow'y 'pomocnicze' lub 'użytkowe', jak oznacza prefiks `utils_` w ich nazwach.
Te subworkflow'y to te, które produkują fantazyjny nagłówek nf-core w wyjściu konsoli, między innymi funkcjami pomocniczymi.

!!! tip "Wskazówka"

    Oprócz wzorca nazewnictwa, innym wskazaniem, że te subworkflow'y nie wykonują żadnej prawdziwie związanej z analizą funkcji, jest to, że w ogóle nie wywołują żadnych procesów.

To kończy przegląd podstawowych komponentów kodu, które stanowią pipeline `nf-core/demo`.
Teraz spójrzmy na pozostałe elementy, o których powinieneś wiedzieć trochę przed zagłębieniem się w rozwój: konfigurację pipeline'u i walidację danych wejściowych.

### 3.2. Konfiguracja pipeline'u

Nauczyłeś się wcześniej, że Nextflow oferuje wiele opcji konfigurowania wykonania pipeline'u, czy to pod względem danych wejściowych i parametrów, zasobów obliczeniowych i innych aspektów orkiestracji.
Projekt nf-core stosuje wysoce ustandaryzowane wytyczne dotyczące konfiguracji pipeline'u, które mają na celu budowanie na elastycznych opcjach dostosowywania Nextflow'a w sposób, który zapewnia większą spójność i łatwość utrzymania w pipeline'ach.

Centralny plik konfiguracyjny `nextflow.config` jest używany do ustawiania wartości domyślnych dla parametrów i innych opcji konfiguracyjnych.
Większość tych opcji konfiguracyjnych jest stosowana domyślnie, podczas gdy inne (np. profile zależności oprogramowania) są dołączone jako opcjonalne profile.

Istnieje kilka dodatkowych plików konfiguracyjnych, które są przechowywane w folderze `conf` i które mogą być dodane do konfiguracji domyślnie lub opcjonalnie jako profile:

- `base.config`: Plik konfiguracyjny 'czystej tablicy', odpowiedni do ogólnego użytku w większości środowisk obliczeniowych o wysokiej wydajności. Definiuje to szerokie przedziały użycia zasobów, na przykład, które są wygodne do zastosowania do modułów.
- `modules.config`: Dodatkowe dyrektywy i argumenty modułów.
- `test.config`: Profil do uruchomienia pipeline'u z minimalnymi danymi testowymi, którego użyliśmy, gdy uruchomiliśmy pipeline demo.
- `test_full.config`: Profil do uruchomienia pipeline'u z pełnowymiarowym zestawem danych testowych.

Dotniemy kilku z tych plików później w kursie.

### 3.3. Dane wejściowe i walidacja

Jak zauważyliśmy wcześniej, gdy badaliśmy profil testowy pipeline'u `nf-core/demo`, jest on zaprojektowany tak, aby przyjmować jako wejście arkusz próbek zawierający ścieżki plików i identyfikatory próbek.
Ścieżki plików wskazywały na rzeczywiste dane znajdujące się w repozytorium `nf-core/test-datasets`.

Przykładowy arkusz próbek jest również dostarczany w katalogu `assets`, chociaż ścieżki w tym nie są prawdziwe.

```csv title="assets/samplesheet.csv" linenums="1"
sample,fastq_1,fastq_2
SAMPLE_PAIRED_END,/path/to/fastq/files/AEG588A1_S1_L002_R1_001.fastq.gz,/path/to/fastq/files/AEG588A1_S1_L002_R2_001.fastq.gz
SAMPLE_SINGLE_END,/path/to/fastq/files/AEG588A4_S4_L003_R1_001.fastq.gz,

```

Ten konkretny arkusz próbek jest dość prosty, ale niektóre pipeline'y działają na arkuszach próbek, które są bardziej złożone, z dużo większą ilością metadanych powiązanych z podstawowymi danymi wejściowymi.

Niestety, ponieważ te pliki mogą być trudne do sprawdzenia wzrokiem, nieprawidłowe formatowanie danych wejściowych jest bardzo częstym źródłem awarii pipeline'u.
Powiązanym problemem jest sytuacja, gdy parametry są podawane nieprawidłowo.

Rozwiązaniem tych problemów jest uruchomienie automatycznych kontroli walidacyjnych na wszystkich plikach wejściowych, aby upewnić się, że zawierają oczekiwane typy informacji, prawidłowo sformatowane, oraz na parametrach, aby upewnić się, że są oczekiwanego typu.
Nazywa się to walidacją danych wejściowych i idealnie powinno być wykonywane _przed_ próbą uruchomienia pipeline'u, zamiast czekać, aż pipeline zawiedzie, aby dowiedzieć się, że był problem z danymi wejściowymi.

Podobnie jak w przypadku konfiguracji, projekt nf-core ma bardzo zdecydowane zdanie na temat walidacji danych wejściowych i zaleca użycie [wtyczki nf-schema](https://nextflow-io.github.io/nf-schema/latest/), wtyczki Nextflow'a, która zapewnia kompleksowe możliwości walidacji dla pipeline'ów Nextflow'a.

Omówimy ten temat bardziej szczegółowo w Części 5 tego kursu.
Na razie po prostu bądź świadomy, że istnieją dwa pliki JSON dostarczone w tym celu, `nextflow_schema.json` i `assets/schema_input.json`.

`nextflow_schema.json` to plik używany do przechowywania informacji o parametrach pipeline'u, w tym typie, opisie i tekście pomocy w formacie czytelnym maszynowo.
Jest to używane do różnych celów, w tym automatycznej walidacji parametrów, generowania tekstu pomocy i renderowania interaktywnych formularzy parametrów w interfejsach użytkownika.

`schema_input.json` to plik używany do definiowania struktury arkusza próbek wejściowych.
Każda kolumna może mieć typ, wzorzec, opis i tekst pomocy w formacie czytelnym maszynowo.
Schemat jest używany do różnych celów, w tym automatycznej walidacji i dostarczania pomocnych komunikatów o błędach.

### Podsumowanie

Wiesz, jakie są główne komponenty pipeline'u nf-core i jak zorganizowany jest kod; gdzie znajdują się główne elementy konfiguracji; i jesteś świadomy, do czego służy walidacja danych wejściowych.

### Co dalej?

Zrób sobie przerwę! To było dużo. Gdy poczujesz się odświeżony i gotowy, przejdź do następnej sekcji, aby zastosować to, czego się nauczyłeś, do napisania pipeline'u kompatybilnego z nf-core.

!!! tip "Wskazówka"

    Jeśli chciałbyś dowiedzieć się, jak komponować workflow'y z subworkflow'ami przed przejściem do następnej części, sprawdź [Workflows of Workflows](../side_quests/workflows_of_workflows.md) Side Quest.
