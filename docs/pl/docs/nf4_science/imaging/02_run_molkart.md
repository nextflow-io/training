# Część 2: Uruchomienie nf-core/molkart

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Tłumaczenie wspomagane przez AI - [dowiedz się więcej i zasugeruj ulepszenia](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

W Części 1 uruchomiliśmy prosty workflow Hello World, aby zrozumieć podstawy wykonywania Nextflow'a.
Teraz uruchomimy rzeczywisty pipeline do bioobrażowania: **nf-core/molkart**.

Ten pipeline przetwarza dane transkryptomiki przestrzennej Molecular Cartography z Resolve Bioscience.
Jednak wzorce Nextflow'a, których się tutaj nauczysz, mają zastosowanie do każdego pipeline'u nf-core lub produkcyjnego workflow'u.

## 1. Zrozumienie pipeline'ów nf-core

Zanim uruchomimy pipeline, zrozummy czym jest nf-core i dlaczego ma znaczenie przy uruchamianiu workflow'ów.

### 1.1. Czym jest nf-core?

[nf-core](https://nf-co.re/) to wspierana przez społeczność kolekcja wysokiej jakości pipeline'ów Nextflow'a.
Wszystkie pipeline'y nf-core mają tę samą strukturę i konwencje, co oznacza, że nauczywszy się uruchamiać jeden, możesz uruchomić każdy z nich.

Kluczowe cechy pipeline'ów nf-core:

- **Znormalizowana struktura**: Wszystkie pipeline'y mają spójne nazwy parametrów i wzorce użycia
- **Wbudowane dane testowe**: Każdy pipeline zawiera profile testowe do szybkiej walidacji
- **Kompleksowa dokumentacja**: Szczegółowe instrukcje użycia i opisy parametrów
- **Kontrola jakości**: Automatyczne raporty QC przy użyciu MultiQC
- **Wsparcie dla kontenerów**: Gotowe kontenery zapewniające powtarzalność

!!! tip "Chcesz dowiedzieć się więcej o nf-core?"

    Aby zapoznać się ze szczegółowym wprowadzeniem do tworzenia pipeline'ów nf-core, sprawdź kurs szkoleniowy [Hello nf-core](../../hello_nf-core/index.md).
    Obejmuje on tworzenie i dostosowywanie pipeline'ów nf-core od podstaw.

### 1.2. Pipeline molkart

![Pipeline nf-core/molkart](img/molkart.png)

Pipeline [nf-core/molkart](https://nf-co.re/molkart) przetwarza dane obrazowania transkryptomiki przestrzennej przez kilka etapów:

1. **Przetwarzanie wstępne obrazu**: Wypełnianie wzoru siatki i opcjonalne wzmocnienie kontrastu
2. **Segmentacja komórek**: Wiele opcji algorytmów (Cellpose, Mesmer, ilastik, Stardist)
3. **Przypisanie punktów**: Przypisywanie punktów transkryptów do zsegmentowanych komórek
4. **Kontrola jakości**: Generowanie kompleksowych raportów QC

Kluczowe wyniki to:

- Tabele liczby transkryptów według komórek
- Maski segmentacji
- Raport kontroli jakości MultiQC

---

## 2. Uruchomienie molkart z danymi testowymi

Zanim zaczniemy, sklonujmy repozytorium molkart lokalnie, abyśmy mogli sprawdzić jego kod:

```bash
cd /workspaces/training/nf4-science/imaging
git clone --branch 1.2.0 --depth 1 https://github.com/nf-core/molkart
```

To tworzy katalog `molkart/` zawierający kompletny kod źródłowy pipeline'u.

!!! note "Dlaczego klonujemy lokalnie?"

    Zazwyczaj uruchamiałbyś pipeline'y nf-core bezpośrednio z GitHub'a używając `nextflow run nf-core/molkart -r 1.2.0`.
    Nextflow automatycznie pobiera żądaną wersję pipeline'u do `$HOME/.nextflow/assets/nf-core/molkart` i uruchamia go stamtąd.
    Jednak dla celów tego szkolenia klonujemy pipeline do innego katalogu lokalnego, abyśmy mogli łatwiej sprawdzić kod.

### 2.1. Zrozumienie wymagań kontenerów

Zanim uruchomimy pełny pipeline, nauczmy się dlaczego kontenery są niezbędne dla pipeline'ów nf-core.

Parametry pipeline'u podamy za pomocą pliku parametrów.
Plik parametrów to plik YAML zawierający listę parametrów i ich wartości, który zachowuje typy wartości (takie jak liczby całkowite) i skraca polecenie uruchomienia.

W katalogu roboczym znajduje się już gotowy plik `params.yaml`:

```yaml title="params.yaml"
input: "data/samplesheet.csv"
outdir: "results"
mindagap_tilesize: 90
mindagap_boxsize: 7
mindagap_loopnum: 100
clahe_pyramid_tile: 368
segmentation_method: "cellpose"
```

Znaczenie tych parametrów:

- `input`: Ścieżka do arkusza próbek zawierającego metadane próbek
- `mindagap_tilesize`, `mindagap_boxsize`, `mindagap_loopnum`: Parametry dla wypełniania wzoru siatki
- `clahe_pyramid_tile`: Rozmiar kernela dla wzmocnienia kontrastu
- `segmentation_method`: Który algorytm(y) użyć do segmentacji komórek
- `outdir`: Gdzie zapisać wyniki

Spróbujmy uruchomić pipeline z tymi parametrami:

```bash
nextflow run ./molkart -params-file params.yaml
```

!!! Warning "To polecenie zakończy się niepowodzeniem - to zamierzone!"

    Celowo uruchamiamy to bez kontenerów, aby pokazać dlaczego są potrzebne.

Po kilku chwilach zobaczysz błąd taki jak ten:

??? failure "Wyjście polecenia"

    ```console
    ERROR ~ Error executing process > 'NFCORE_MOLKART:MOLKART:MINDAGAP_DUPLICATEFINDER (mem_only)'

    Caused by:
      Process `NFCORE_MOLKART:MOLKART:MINDAGAP_DUPLICATEFINDER (mem_only)` terminated with an error exit status (127)

    Command executed:

      duplicate_finder.py \
          spots.txt \
          90

    Command exit status:
      127

    Command error:
      .command.sh: line 3: duplicate_finder.py: command not found
    ```

**Co się dzieje?**

Błąd `command not found` (status wyjścia 127) oznacza, że Nextflow próbował uruchomić `duplicate_finder.py`, ale nie mógł go znaleźć w Twoim systemie.
Dzieje się tak, ponieważ:

1. Pipeline oczekuje zainstalowanego specjalistycznego oprogramowania bioinformatycznego
2. Te narzędzia (takie jak `duplicate_finder.py`, `apply_clahe.dask.py`, itp.) nie są częścią standardowych dystrybucji Linuksa
3. Bez kontenerów Nextflow próbuje uruchomić polecenia bezpośrednio na Twojej lokalnej maszynie

**Skąd mają pochodzić te narzędzia?**

Sprawdźmy jeden z modułów procesu, aby zobaczyć jak deklaruje swoje wymagania dotyczące oprogramowania.

Otwórz moduł przetwarzania wstępnego CLAHE:

```bash
code molkart/modules/local/clahe/main.nf
```

Spójrz na linię 5 - zobaczysz:

```groovy
container 'ghcr.io/schapirolabor/molkart-local:v0.0.4'
```

Ta linia informuje Nextflow'a: "Aby uruchomić ten proces, użyj obrazu Docker'a `ghcr.io/schapirolabor/molkart-local:v0.0.4`, który zawiera wszystkie wymagane oprogramowanie."

Każdy proces deklaruje, który obraz kontenera dostarcza wymagane narzędzia.
Jednak Nextflow używa tych kontenerów tylko wtedy, gdy mu powiesz!

**Rozwiązanie: Włącz Docker w konfiguracji**

### 2.2. Konfiguracja Docker'a i uruchomienie pipeline'u

Aby włączyć Docker'a, musimy zmienić `docker.enabled` z `false` na `true` w pliku `nextflow.config`.

Otwórz plik konfiguracyjny:

```bash
code nextflow.config
```

Zmień `docker.enabled = false` na `docker.enabled = true`:

```groovy
docker.enabled = true
process {
    resourceLimits = [
        cpus: 2,
        memory: '7.GB',
    ]
}
```

Teraz uruchom pipeline ponownie, tym razem z wszystkimi trzema metodami segmentacji, abyśmy mogli je później porównać:

```bash
nextflow run ./molkart -params-file params.yaml --segmentation_method "mesmer,cellpose,stardist"
```

Tym razem Nextflow:

1. Odczyta ustawienie `docker.enabled = true` z konfiguracji
2. Pobierze wymagane obrazy Docker'a (tylko za pierwszym razem)
3. Uruchomi każdy proces wewnątrz określonego kontenera
4. Wykona się pomyślnie, ponieważ wszystkie narzędzia są dostępne wewnątrz kontenerów

!!! Tip "Dlaczego kontenery są ważne"

    Większość pipeline'ów nf-core **wymaga** konteneryzacji (Docker, Singularity, Podman, itp.), ponieważ:

    - Używają specjalistycznego oprogramowania bioinformatycznego niedostępnego w standardowych środowiskach
    - Kontenery zapewniają powtarzalność - dokładnie te same wersje oprogramowania działają wszędzie
    - Nie musisz ręcznie instalować dziesiątek narzędzi i ich zależności

    Aby uzyskać więcej szczegółów o kontenerach w Nextflow'ie, zobacz [Hello Containers](../../hello_nextflow/05_hello_containers.md) ze szkolenia Hello Nextflow.

### 2.3. Monitorowanie wykonania

Podczas działania pipeline'u zobaczysz wyjście podobne do tego:

??? success "Wyjście polecenia"

    ```console
     N E X T F L O W   ~  version 26.04.4

    Launching `./molkart/main.nf` [exotic_banach] revision: e4152308ec

    WARN: Unrecognized config option 'validation.help.enabled'
    WARN: Unrecognized config option 'validation.defaultIgnoreParams'
    WARN: Unrecognized config option 'validation.monochromeLogs'

    ------------------------------------------------------
                                            ,--./,-.
            ___     __   __   __   ___     /,-._.--~'
      |\ | |__  __ /  ` /  \ |__) |__         }  {
      | \| |       \__, \__/ |  \ |___     \`-._,-`-,
                                            `._,._,'
      nf-core/molkart 1.2.0
    ------------------------------------------------------
    Segmentation methods and options
      segmentation_method: mesmer,cellpose,stardist

    Image preprocessing
      mindagap_boxsize   : 7
      mindagap_loopnum   : 100
      clahe_kernel       : 25
      mindagap_tilesize  : 90
      clahe_pyramid_tile : 368

    Input/output options
      input              : data/samplesheet.csv
      outdir             : results

    Generic options
      trace_report_suffix: 2026-06-23_16-25-30

    Core Nextflow options
      runName            : exotic_banach
      containerEngine    : docker
      launchDir          : /workspaces/training/nf4-science/imaging
      workDir            : /workspaces/training/nf4-science/imaging/work
      projectDir         : /workspaces/training/nf4-science/imaging/molkart
      userName           : root
      profile            : standard
      configFiles        : /workspaces/training/nf4-science/imaging/molkart/nextflow.config, /workspaces/training/nf4-science/imaging/nextflow.config

    !! Only displaying parameters that differ from the pipeline defaults !!
    ------------------------------------------------------
    * The pipeline
        https://doi.org/10.5281/zenodo.10650748

    * The nf-core framework
        https://doi.org/10.1038/s41587-020-0439-x

    * Software dependencies
        https://github.com/nf-core/molkart/blob/master/CITATIONS.md

    executor >  local (22)
    [b4/e57ff1] NFC…NDAGAP_MINDAGAP (mem_only) | 2 of 2 ✔
    [2e/cf8910] NFC…T:MOLKART:CLAHE (mem_only) | 2 of 2 ✔
    [94/b1ef70] NFC…RT:CREATE_STACK (mem_only) | 1 of 1 ✔
    [58/4b6426] NFC…DUPLICATEFINDER (mem_only) | 1 of 1 ✔
    [c2/ef7002] NFC…DEEPCELL_MESMER (mem_only) | 1 of 1 ✔
    [4f/ddd328] NFC…OLKART:STARDIST (mem_only) | 1 of 1 ✔
    [8c/dd0362] NFC…OLKART:CELLPOSE (mem_only) | 1 of 1 ✔
    [08/2ddcb5] NFC…KART:MASKFILTER (mem_only) | 3 of 3 ✔
    [56/e30de0] NFC…LKART:SPOT2CELL (mem_only) | 3 of 3 ✔
    [b0/5aa635] NFC…:CREATE_ANNDATA (mem_only) | 3 of 3 ✔
    [5e/39d5d0] NFC…LKART:MOLKARTQC (mem_only) | 3 of 3 ✔
    [8e/8bd365] NFCORE_MOLKART:MOLKART:MULTIQC | 1 of 1 ✔
    -[nf-core/molkart] Pipeline completed successfully-
    Completed at: 23-Jun-2026 16:31:40
    Duration    : 6m 9s
    CPU hours   : (a few seconds)
    Succeeded   : 22
    ```

Zauważ, jak to wyjście jest bardziej szczegółowe niż nasz przykład Hello World dzięki konwencjom nf-core, których przestrzega pipeline:

- Pipeline pokazuje swoją wersję i logo
- Wyświetlane są parametry konfiguracji
- Wiele procesów działa równolegle (wskazane przez wiele linii procesów)
- Nazwy procesów zawierają pełną ścieżkę modułu (np. `NFCORE_MOLKART:MOLKART:MINDAGAP_MINDAGAP`)

### 2.4. Zrozumienie wykonania procesu

Linia executor `executor > local (22)` informuje Cię:

- **executor**: Które środowisko obliczeniowe jest używane (`local` = Twoja maszyna)
- **(22)**: Całkowita liczba uruchomionych zadań

Każda linia procesu pokazuje:

- **Hash** (`[b4/e57ff1]`): Identyfikator katalogu roboczego (jak wcześniej)
- **Nazwa procesu**: Pełna ścieżka modułu i nazwa procesu
- **Identyfikator wejścia**: Nazwa próbki w nawiasach
- **Postęp**: Liczba zadań i status ukończenia (np. `1 of 1 ✔`)

### Podsumowanie

Wiesz jak uruchomić pipeline nf-core z danymi testowymi i interpretować wyjście jego wykonania.

### Co dalej?

Naucz się gdzie znaleźć wyniki i jak je interpretować.

---

## 3. Znajdowanie i badanie wyjść

Kiedy pipeline zakończy się pomyślnie, zobaczysz komunikat o zakończeniu i podsumowanie wykonania.

### 3.1. Zlokalizowanie katalogu wyników

Domyślnie pipeline'y nf-core zapisują wyniki do katalogu określonego przez parametr `outdir`, który ustawiliśmy na `results/`.

Wyświetl zawartość:

```bash
tree results/
```

Powinieneś zobaczyć kilka podkatalogów:

```console title="results/"
results/
├── anndata/
├── clahe/
├── mindagap/
├── molkartqc/
├── multiqc/
├── pipeline_info/
├── segmentation/
├── spot2cell/
└── stack/
```

Każdy podkatalog zawiera wyniki z określonego etapu pipeline'u:

- **mindagap/**: Obrazy z wypełnioną siatką z kroku przetwarzania wstępnego MindaGap
- **clahe/**: Obrazy ze wzmocnionym kontrastem z przetwarzania wstępnego CLAHE
- **stack/**: Stosy obrazów wielokanałowych utworzone do segmentacji
- **segmentation/**: Wyniki segmentacji z różnych algorytmów (cellpose/, mesmer/, stardist/, filtered_masks/)
- **spot2cell/**: Tabele liczby transkryptów według komórek
- **anndata/**: Obiekty AnnData zawierające macierze transkryptów według komórek i współrzędne przestrzenne
- **molkartqc/**: Metryki kontroli jakości dla przypisania punktów
- **multiqc/**: Kompleksowy raport kontroli jakości
- **pipeline_info/**: Raporty wykonania i logi

### 3.2. Badanie raportu MultiQC

Raport MultiQC to kompleksowy plik HTML, który agreguje metryki jakości ze wszystkich kroków pipeline'u.

Otwórz raport w przeglądarce plików, a następnie kliknij przycisk "Show Preview", aby zobaczyć go renderowanego bezpośrednio w VS Code.

Raport zawiera:

- Ogólne statystyki dla wszystkich próbek
- Metryki przetwarzania wstępnego
- Metryki jakości segmentacji
- Liczbę wykrytych komórek i punktów

!!! Tip "Wskazówka"

    Raporty MultiQC są zwykle dołączane do wszystkich pipeline'ów nf-core.
    Zawsze zapewniają ogólny przegląd wykonania pipeline'u i jakości danych.

### 3.3. Badanie tabel transkryptów według komórek

Najważniejszym wynikiem naukowym jest tabela liczby transkryptów według komórek.
Mówi ona, ile transkryptów każdego typu zostało wykrytych w każdej komórce.

Przejdź do katalogu spot2cell:

```bash
ls results/spot2cell/
```

Znajdziesz pliki takie jak:

- `cellxgene_mem_only_cellpose.csv`: Tabela transkryptów według komórek używająca segmentacji Cellpose
- `cellxgene_mem_only_mesmer.csv`: Tabela transkryptów według komórek używająca segmentacji Mesmer
- `cellxgene_mem_only_stardist.csv`: Tabela transkryptów według komórek używająca segmentacji Stardist

Uruchomiliśmy tylko jedną próbkę w tym zestawie danych testowych, ale w prawdziwym eksperymencie mielibyśmy te tabele dla każdej próbki.
Zauważ, jak Nextflow potrafi przetwarzać wiele metod segmentacji równolegle, ułatwiając porównywanie wyników.

### 3.4. Przeglądanie raportów wykonania

Nextflow automatycznie generuje kilka raportów wykonania.

Sprawdź katalog pipeline_info:

```bash
ls results/pipeline_info/
```

Kluczowe pliki:

- **execution_report.html**: Oś czasu i wizualizacja użycia zasobów
- **execution_timeline.html**: Wykres Gantta wykonania procesu
- **execution_trace.txt**: Szczegółowe metryki wykonania zadań
- **pipeline_dag.html**: Skierowany graf acykliczny pokazujący strukturę workflow'u

Otwórz raport wykonania, aby zobaczyć użycie zasobów:

```bash
code results/pipeline_info/execution_report.html
```

Pokazuje on:

- Jak długo trwał każdy proces
- Użycie CPU i pamięci
- Które zadania były w pamięci podręcznej, a które wykonane

!!! Tip "Wskazówka"

    Te raporty są niezwykle przydatne do optymalizacji alokacji zasobów i rozwiązywania problemów z wydajnością.

### Podsumowanie

Wiesz jak zlokalizować wyniki pipeline'u, badać raporty kontroli jakości i uzyskać dostęp do metryk wykonania.

### Co dalej?

Naucz się o katalogu roboczym i jak Nextflow zarządza plikami pośrednimi.

---

## 4. Eksploracja katalogu roboczego

Tak jak w naszym przykładzie Hello World, cała rzeczywista praca odbywa się w katalogu `work/`.

### 4.1. Zrozumienie struktury katalogu roboczego

Katalog roboczy zawiera podkatalog dla każdego zadania, które zostało wykonane.
Dla tego pipeline'u z 22 zadaniami będzie 22 podkatalogów roboczych.

Wyświetl listę katalogu roboczego:

```bash
ls -d work/*/*/ | head -5
```

To pokazuje pierwsze 5 katalogów zadań.

### 4.2. Inspekcja katalogu zadania

Wybierz jeden z hashy procesu segmentacji z wyjścia konsoli (np. `[3m/4n5o6p]`) i zajrzyj do środka:

```bash
ls -la work/3m/4n5o6p*/
```

Zobaczysz:

- **Pliki .command.\***: Skrypty wykonania Nextflow'a i logi (jak wcześniej)
- **Przygotowane pliki wejściowe**: Dowiązania symboliczne do rzeczywistych plików wejściowych
- **Pliki wyjściowe**: Maski segmentacji, wyniki pośrednie, itp.

Kluczowa różnica od Hello World:

- Rzeczywiste pipeline'y przygotowują duże pliki wejściowe (obrazy, dane referencyjne)
- Pliki wyjściowe mogą być dość duże (maski segmentacji, przetworzone obrazy)
- Wiele plików wejściowych i wyjściowych na zadanie

!!! Tip "Wskazówka"

    Jeśli proces się nie powiedzie, możesz przejść do jego katalogu roboczego, sprawdzić `.command.err` w poszukiwaniu komunikatów o błędach, a nawet ponownie uruchomić `.command.sh` ręcznie, aby debugować problem.

### 4.3. Czyszczenie katalogu roboczego

Katalog roboczy może stać się dość duży po wielu uruchomieniach pipeline'u.
Jak nauczyliśmy się w Części 1, możesz użyć `nextflow clean`, aby usunąć katalogi robocze ze starych uruchomień.

Jednak dla pipeline'ów nf-core z dużymi plikami pośrednimi szczególnie ważne jest regularne czyszczenie.

### Podsumowanie

Rozumiesz jak pipeline'y nf-core organizują swoje katalogi robocze i jak sprawdzać poszczególne zadania w celu debugowania.

### Co dalej?

Naucz się o pamięci podręcznej Nextflow'a i jak wznowić nieudane uruchomienia pipeline'u.

---

## 5. Wznowienie uruchomienia pipeline'u

Jedną z najpotężniejszych funkcji Nextflow'a jest możliwość wznowienia pipeline'u od momentu niepowodzenia.

### 5.1. Mechanizm pamięci podręcznej

Kiedy uruchamiasz pipeline z `-resume`, Nextflow:

1. Sprawdza pamięć podręczną dla każdego zadania
2. Jeśli dane wejściowe, kod i parametry są identyczne, ponownie używa wyniku z pamięci podręcznej
3. Ponownie uruchamia tylko zadania, które się zmieniły lub zakończyły się niepowodzeniem

Jest to niezbędne dla długo działających pipeline'ów, gdzie niepowodzenia mogą wystąpić późno w wykonaniu.

### 5.2. Wypróbuj resume z molkart

Uruchom to samo polecenie ponownie, ale dodaj `-resume`:

```bash
nextflow run ./molkart -params-file params.yaml --segmentation_method "mesmer,cellpose,stardist" -resume
```

Powinieneś zobaczyć wyjście takie jak:

```console
executor >  local (1)
[43/e03702] NFC…NDAGAP_MINDAGAP (mem_only) | 2 of 2, cached: 2 ✔
[2e/cf8910] NFC…T:MOLKART:CLAHE (mem_only) | 2 of 2, cached: 2 ✔
[94/b1ef70] NFC…RT:CREATE_STACK (mem_only) | 1 of 1, cached: 1 ✔
[58/4b6426] NFC…DUPLICATEFINDER (mem_only) | 1 of 1, cached: 1 ✔
[c2/ef7002] NFC…DEEPCELL_MESMER (mem_only) | 1 of 1, cached: 1 ✔
[4f/ddd328] NFC…OLKART:STARDIST (mem_only) | 1 of 1, cached: 1 ✔
[8c/dd0362] NFC…OLKART:CELLPOSE (mem_only) | 1 of 1, cached: 1 ✔
[08/2ddcb5] NFC…KART:MASKFILTER (mem_only) | 3 of 3, cached: 3 ✔
[56/e30de0] NFC…LKART:SPOT2CELL (mem_only) | 3 of 3, cached: 3 ✔
[b0/5aa635] NFC…:CREATE_ANNDATA (mem_only) | 3 of 3, cached: 3 ✔
[5e/39d5d0] NFC…LKART:MOLKARTQC (mem_only) | 3 of 3, cached: 3 ✔
[73/239f45] NFCORE_MOLKART:MOLKART:MULTIQC | 1 of 1 ✔
-[nf-core/molkart] Pipeline completed successfully-
```

Zauważ adnotację `cached: N` przy każdym procesie przetwarzania wstępnego i segmentacji - te zadania zostały ponownie użyte zamiast wykonane od nowa.

### 5.3. Kiedy resume jest przydatne

Resume jest szczególnie wartościowe, gdy:

- Pipeline kończy się niepowodzeniem z powodu limitów zasobów (brak pamięci, przekroczenie limitu czasu)
- Musisz zmodyfikować procesy końcowe bez ponownego uruchamiania kroków początkowych
- Twoje połączenie sieciowe zostanie przerwane podczas pobierania danych
- Chcesz dodać dodatkowe wyniki bez powtarzania obliczeń

!!! Warning "Ostrzeżenie"

    Resume działa tylko wtedy, gdy nie zmieniłeś danych wejściowych, kodu pipeline'u ani parametrów.
    Jeśli zmienisz którekolwiek z nich, Nextflow poprawnie ponownie uruchomi dotknięte zadania.

### Podsumowanie

Wiesz już, jak używać `-resume`, aby efektywnie ponownie uruchamiać pipeline'y bez powtarzania udanych zadań.

### Co dalej?

Teraz, gdy możesz uruchomić nf-core/molkart z danymi testowymi, jesteś gotowy nauczyć się, jak skonfigurować go dla własnych zestawów danych.
