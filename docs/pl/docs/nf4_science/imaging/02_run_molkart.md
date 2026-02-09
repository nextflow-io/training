# Część 2: Uruchomienie nf-core/molkart

W Części 1 uruchomiliśmy prosty workflow Hello World, aby zrozumieć podstawy wykonywania Nextflow'a.
Teraz uruchomimy rzeczywisty pipeline bioinformatyczny: **nf-core/molkart**.

Ten pipeline przetwarza dane transkryptomiki przestrzennej Molecular Cartography od Resolve Bioscience.
Wzorce Nextflow'a, których się tutaj nauczysz, mają jednak zastosowanie do każdego pipeline'u nf-core i każdego workflow'a produkcyjnego.

## 1. Zrozumienie pipeline'ów nf-core

Zanim uruchomimy pipeline'a, zrozummy czym jest nf-core i dlaczego ma znaczenie dla uruchamiania workflow'ów.

### 1.1. Czym jest nf-core?

[nf-core](https://nf-co.re/) to wspierana przez społeczność kolekcja wysokiej jakości pipeline'ów Nextflow'a.
Wszystkie pipeline'y nf-core mają tę samą strukturę i konwencje, co oznacza, że gdy nauczysz się uruchamiać jeden z nich, będziesz potrafić uruchomić każdy.

Kluczowe cechy pipeline'ów nf-core:

- **Ustandaryzowana struktura**: Wszystkie pipeline'y mają spójne nazwy parametrów i wzorce użycia
- **Wbudowane dane testowe**: Każdy pipeline zawiera profile testowe do szybkiej walidacji
- **Kompleksowa dokumentacja**: Szczegółowe instrukcje użycia i opisy parametrów
- **Kontrola jakości**: Automatyczne raporty QC przy użyciu MultiQC
- **Wsparcie dla kontenerów**: Gotowe kontenery zapewniające powtarzalność

!!! tip "Chcesz dowiedzieć się więcej o nf-core?"

    Aby poznać dogłębne wprowadzenie do tworzenia pipeline'ów nf-core, sprawdź kurs szkoleniowy [Hello nf-core](../../hello_nf-core/index.md).
    Obejmuje on tworzenie i dostosowywanie pipeline'ów nf-core od podstaw.

### 1.2. Pipeline molkart

![Pipeline nf-core/molkart](img/molkart.png)

Pipeline [nf-core/molkart](https://nf-co.re/molkart) przetwarza dane obrazowania transkryptomiki przestrzennej przez kilka etapów:

1. **Wstępne przetwarzanie obrazu**: Wypełnianie wzoru siatki i opcjonalne zwiększanie kontrastu
2. **Segmentacja komórek**: Wiele opcji algorytmów (Cellpose, Mesmer, ilastik, Stardist)
3. **Przypisywanie punktów**: Przypisanie punktów transkryptów do segmentowanych komórek
4. **Kontrola jakości**: Generowanie kompleksowych raportów QC

Kluczowe wyniki to:

- Tabele liczności komórka-transkrypt
- Maski segmentacji
- Raport kontroli jakości MultiQC

---

## 2. Uruchomienie molkart z danymi testowymi

Zanim zaczniemy, sklonujmy repozytorium molkart lokalnie, abyśmy mogli przejrzeć jego kod:

```bash
cd /workspaces/training/nf4-science/imaging
git clone --branch 1.2.0 --depth 1 https://github.com/nf-core/molkart
```

To tworzy katalog `molkart/` zawierający kompletny kod źródłowy pipeline'u.

!!! note "Dlaczego klonujemy lokalnie?"

    Zazwyczaj uruchamiałbyś pipeline'y nf-core bezpośrednio z GitHuba używając `nextflow run nf-core/molkart -r 1.2.0`.
    Nextflow automatycznie pobiera żądaną wersję pipeline'u do `$HOME/.nextflow/assets/nf-core/molkart` i uruchamia go stamtąd.
    Jednak na potrzeby tego szkolenia klonujemy pipeline'a do innego lokalnego katalogu, abyśmy mogli łatwiej przejrzeć kod.

### 2.1. Zrozumienie wymagań dotyczących kontenerów

Zanim uruchomimy pełny pipeline, dowiedzmy się, dlaczego kontenery są niezbędne dla pipeline'ów nf-core.

Spróbujmy uruchomić pipeline'a używając zestawu danych testowych i parametrów z konfiguracji testowej molkart:

```bash
nextflow run ./molkart \
  --input 'data/samplesheet.csv' \
  --mindagap_tilesize 90 \
  --mindagap_boxsize 7 \
  --mindagap_loopnum 100 \
  --clahe_pyramid_tile 368 \
  --segmentation_method "mesmer,cellpose,stardist" \
  --outdir results
```

Rozbijmy te parametry:

- `--input`: Ścieżka do arkusza próbek zawierającego metadane próbek
- `--mindagap_tilesize`, `--mindagap_boxsize`, `--mindagap_loopnum`: Parametry wypełniania wzoru siatki
- `--clahe_pyramid_tile`: Rozmiar jądra do zwiększania kontrastu
- `--segmentation_method`: Który algorytm (lub algorytmy) użyć do segmentacji komórek
- `--outdir`: Gdzie zapisać wyniki

!!! Warning "To polecenie zakończy się niepowodzeniem - to zamierzone!"

    Celowo uruchamiamy to bez kontenerów, aby pokazać, dlaczego są potrzebne.

Po kilku chwilach zobaczysz błąd podobny do tego:

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
2. Te narzędzia (takie jak `duplicate_finder.py`, `apply_clahe.dask.py` itp.) nie są częścią standardowych dystrybucji Linuksa
3. Bez kontenerów Nextflow próbuje uruchamiać polecenia bezpośrednio na Twojej lokalnej maszynie

**Skąd mają pochodzić te narzędzia?**

Przyjrzyjmy się jednemu z modułów procesów, aby zobaczyć, jak deklaruje swoje wymagania dotyczące oprogramowania.

Otwórz moduł wstępnego przetwarzania CLAHE:

```bash
code molkart/modules/local/clahe/main.nf
```

Spójrz na linię 5 - zobaczysz:

```groovy
container 'ghcr.io/schapirolabor/molkart-local:v0.0.4'
```

Ta linia mówi Nextflow'owi: "Aby uruchomić ten proces, użyj obrazu Dockera `ghcr.io/schapirolabor/molkart-local:v0.0.4`, który zawiera wszystkie wymagane oprogramowanie."

Każdy proces deklaruje, który obraz kontenera dostarcza jego wymaganych narzędzi.
Jednak Nextflow używa tych kontenerów tylko wtedy, gdy mu to powiesz!

**Rozwiązanie: Włącz Dockera w konfiguracji**

### 2.2. Konfiguracja Dockera i uruchomienie pipeline'u

Aby włączyć Dockera, musimy zmienić `docker.enabled` z `false` na `true` w pliku `nextflow.config`.

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

Teraz uruchom pipeline'a ponownie tym samym poleceniem:

```bash
nextflow run ./molkart \
  --input 'data/samplesheet.csv' \
  --mindagap_tilesize 90 \
  --mindagap_boxsize 7 \
  --mindagap_loopnum 100 \
  --clahe_pyramid_tile 368 \
  --segmentation_method "cellpose,mesmer,stardist" \
  --outdir results
```

Tym razem Nextflow:

1. Odczyta ustawienie `docker.enabled = true` z konfiguracji
2. Pobierze wymagane obrazy Dockera (tylko za pierwszym razem)
3. Uruchomi każdy proces wewnątrz określonego kontenera
4. Wykona się pomyślnie, ponieważ wszystkie narzędzia są dostępne wewnątrz kontenerów

!!! Tip "Dlaczego kontenery mają znaczenie"

    Większość pipeline'ów nf-core **wymaga** konteneryzacji (Docker, Singularity, Podman itp.), ponieważ:

    - Używają specjalistycznego oprogramowania bioinformatycznego niedostępnego w standardowych środowiskach
    - Kontenery zapewniają powtarzalność - dokładnie te same wersje oprogramowania działają wszędzie
    - Nie musisz ręcznie instalować dziesiątek narzędzi i ich zależności

    Więcej szczegółów o kontenerach w Nextflow znajdziesz w [Hello Containers](../../hello_nextflow/05_hello_containers.md) ze szkolenia Hello Nextflow.

### 2.3. Monitorowanie wykonania

Podczas działania pipeline'u zobaczysz wyjście podobne do tego:

??? success "Wyjście polecenia"

    ```console
    Nextflow 25.04.8 is available - Please consider updating your version to it

    N E X T F L O W   ~  version 25.04.3

    Launching `https://github.com/nf-core/molkart` [soggy_kalam] DSL2 - revision: 5e54b29cb3 [dev]


    ------------------------------------------------------
                                            ,--./,-.
            ___     __   __   __   ___     /,-._.--~'
      |\ | |__  __ /  ` /  \ |__) |__         }  {
      | \| |       \__, \__/ |  \ |___     \`-._,-`-,
                                            `._,._,'
      nf-core/molkart 1.2.0dev
    ------------------------------------------------------
    Segmentation methods and options
      segmentation_method       : mesmer,cellpose,stardist

    Image preprocessing
      mindagap_boxsize          : 7
      mindagap_loopnum          : 100
      clahe_kernel              : 25
      mindagap_tilesize         : 90
      clahe_pyramid_tile        : 368

    Input/output options
      input                     : https://raw.githubusercontent.com/nf-core/test-datasets/molkart/test_data/samplesheets/samplesheet_membrane.csv
      outdir                    : results

    Institutional config options
      config_profile_name       : Test profile
      config_profile_description: Minimal test dataset to check pipeline function

    Generic options
      trace_report_suffix       : 2025-10-18_22-22-21

    Core Nextflow options
      revision                  : dev
      runName                   : soggy_kalam
      containerEngine           : docker
      launchDir                 : /workspaces/training/nf4-science/imaging
      workDir                   : /workspaces/training/nf4-science/imaging/work
      projectDir                : /workspaces/.nextflow/assets/nf-core/molkart
      userName                  : root
      profile                   : docker,test
      configFiles               :

    !! Only displaying parameters that differ from the pipeline defaults !!
    ------------------------------------------------------
    * The pipeline
        https://doi.org/10.5281/zenodo.10650748

    * The nf-core framework
        https://doi.org/10.1038/s41587-020-0439-x

    * Software dependencies
        https://github.com/nf-core/molkart/blob/master/CITATIONS.md

    executor >  local (22)
    [c1/da5009] NFCORE_MOLKART:MOLKART:MINDAGAP_MINDAGAP (mem_only)        [100%] 2 of 2 ✔
    [73/8f5e8a] NFCORE_MOLKART:MOLKART:CLAHE (mem_only)                    [100%] 2 of 2 ✔
    [ec/8f84d5] NFCORE_MOLKART:MOLKART:CREATE_STACK (mem_only)             [100%] 1 of 1 ✔
    [a2/99349b] NFCORE_MOLKART:MOLKART:MINDAGAP_DUPLICATEFINDER (mem_only) [100%] 1 of 1 ✔
    [95/c9b4b1] NFCORE_MOLKART:MOLKART:DEEPCELL_MESMER (mem_only)          [100%] 1 of 1 ✔
    [d4/1ebd1e] NFCORE_MOLKART:MOLKART:STARDIST (mem_only)                 [100%] 1 of 1 ✔
    [3e/3c0736] NFCORE_MOLKART:MOLKART:CELLPOSE (mem_only)                 [100%] 1 of 1 ✔
    [a0/415c6a] NFCORE_MOLKART:MOLKART:MASKFILTER (mem_only)               [100%] 3 of 3 ✔
    [14/a830c9] NFCORE_MOLKART:MOLKART:SPOT2CELL (mem_only)                [100%] 3 of 3 ✔
    [b5/391836] NFCORE_MOLKART:MOLKART:CREATE_ANNDATA (mem_only)           [100%] 3 of 3 ✔
    [77/aed558] NFCORE_MOLKART:MOLKART:MOLKARTQC (mem_only)                [100%] 3 of 3 ✔
    [e6/b81475] NFCORE_MOLKART:MOLKART:MULTIQC                             [100%] 1 of 1 ✔
    -[nf-core/molkart] Pipeline completed successfully-
    Completed at: 19-Oct-2025 22:23:01
    Duration    : 2m 52s
    CPU hours   : 0.1
    Succeeded   : 22
    ```

Zauważ, że to wyjście jest bardziej szczegółowe niż nasz przykład Hello World, ze względu na konwencje nf-core, których przestrzega pipeline:

- Pipeline pokazuje swoją wersję i logo
- Wyświetlane są parametry konfiguracji
- Wiele procesów działa równolegle (wskazane przez wiele linii procesów)
- Nazwy procesów zawierają pełną ścieżkę modułu (np. `NFCORE_MOLKART:MOLKART:MINDAGAP_MINDAGAP`)

### 2.4. Zrozumienie wykonywania procesów

Linia executora `executor > local (22)` mówi Ci:

- **executor**: Które środowisko obliczeniowe jest używane (`local` = Twoja maszyna)
- **(22)**: Całkowita liczba uruchomionych zadań

Każda linia procesu pokazuje:

- **Hash** (`[1a/2b3c4d]`): Identyfikator katalogu roboczego (jak wcześniej)
- **Nazwa procesu**: Pełna ścieżka modułu i nazwa procesu
- **Identyfikator wejścia**: Nazwa próbki w nawiasach
- **Postęp**: Procent ukończenia i liczba (np. `1 of 1 ✔`)

### Podsumowanie

Wiesz, jak uruchomić pipeline nf-core z danymi testowymi i interpretować jego wyjście wykonania.

### Co dalej?

Dowiedz się, gdzie znaleźć wyniki i jak je interpretować.

---

## 3. Znajdowanie i badanie wyników

Gdy pipeline zakończy się pomyślnie, zobaczysz komunikat o zakończeniu i podsumowanie wykonania.

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

- **mindagap/**: Obrazy z wypełnioną siatką z kroku wstępnego przetwarzania MindaGap
- **clahe/**: Obrazy ze zwiększonym kontrastem z wstępnego przetwarzania CLAHE
- **stack/**: Wielokanałowe stosy obrazów utworzone do segmentacji
- **segmentation/**: Wyniki segmentacji z różnych algorytmów (cellpose/, mesmer/, stardist/, filtered_masks/)
- **spot2cell/**: Tabele liczności komórka-transkrypt
- **anndata/**: Obiekty AnnData zawierające macierze komórka-transkrypt i współrzędne przestrzenne
- **molkartqc/**: Metryki kontroli jakości przypisywania punktów
- **multiqc/**: Kompleksowy raport kontroli jakości
- **pipeline_info/**: Raporty wykonania i logi

### 3.2. Przejrzenie raportu MultiQC

Raport MultiQC to kompleksowy plik HTML agregujący metryki jakości ze wszystkich kroków pipeline'u.

Otwórz raport w przeglądarce plików, a następnie kliknij przycisk "Show Preview", aby zobaczyć go wyrenderowanego bezpośrednio w VS Code.

Raport zawiera:

- Ogólne statystyki dla wszystkich próbek
- Metryki wstępnego przetwarzania
- Metryki jakości segmentacji
- Liczbę wykrytych komórek i punktów

!!! Tip

    Raporty MultiQC są zazwyczaj dołączane do wszystkich pipeline'ów nf-core.
    Zawsze zapewniają ogólny przegląd wykonania pipeline'u i jakości danych.

### 3.3. Przejrzenie tabel komórka-transkrypt

Najważniejszym wynikiem naukowym jest tabela liczności komórka-transkrypt.
Mówi Ci ona, ile każdego transkryptu wykryto w każdej komórce.

Przejdź do katalogu spot2cell:

```bash
ls results/spot2cell/
```

Znajdziesz pliki takie jak:

- `cellxgene_mem_only_cellpose.csv`: Tabela komórka-transkrypt używająca segmentacji Cellpose
- `cellxgene_mem_only_mesmer.csv`: Tabela komórka-transkrypt używająca segmentacji Mesmer
- `cellxgene_mem_only_stardist.csv`: Tabela komórka-transkrypt używająca segmentacji Stardist

Uruchomiliśmy tylko 1 próbkę w tym zestawie danych testowych, ale w prawdziwym eksperymencie mielibyśmy te tabele dla każdej próbki.
Zauważ, jak Nextflow jest w stanie przetwarzać wiele metod segmentacji równolegle, ułatwiając porównywanie wyników.

### 3.4. Przeglądanie raportów wykonania

Nextflow automatycznie generuje kilka raportów wykonania.

Sprawdź katalog pipeline_info:

```bash
ls results/pipeline_info/
```

Kluczowe pliki:

- **execution_report.html**: Wizualizacja osi czasu i użycia zasobów
- **execution_timeline.html**: Wykres Gantta wykonania procesów
- **execution_trace.txt**: Szczegółowe metryki wykonania zadań
- **pipeline_dag.html**: Skierowany graf acykliczny pokazujący strukturę workflow'a

Otwórz raport wykonania, aby zobaczyć użycie zasobów:

```bash
code results/pipeline_info/execution_report.html
```

To pokazuje:

- Jak długo trwał każdy proces
- Użycie CPU i pamięci
- Które zadania były w pamięci podręcznej, a które wykonane

!!! Tip

    Te raporty są niezwykle przydatne do optymalizacji alokacji zasobów i rozwiązywania problemów z wydajnością.

### Podsumowanie

Wiesz, jak zlokalizować wyniki pipeline'u, przejrzeć raporty kontroli jakości i uzyskać dostęp do metryk wykonania.

### Co dalej?

Dowiedz się o katalogu roboczym i jak Nextflow zarządza plikami pośrednimi.

---

## 4. Eksploracja katalogu roboczego

Tak jak w naszym przykładzie Hello World, cała rzeczywista praca odbywa się w katalogu `work/`.

### 4.1. Zrozumienie struktury katalogu roboczego

Katalog roboczy zawiera podkatalog dla każdego zadania, które zostało wykonane.
Dla tego pipeline'u z 12 zadaniami będzie 12 podkatalogów roboczych.

Wyświetl katalog roboczy:

```bash
ls -d work/*/*/ | head -5
```

To pokazuje pierwsze 5 katalogów zadań.

### 4.2. Przejrzenie katalogu zadania

Wybierz jeden z haszy procesu segmentacji z wyjścia konsoli (np. `[3m/4n5o6p]`) i zajrzyj do środka:

```bash
ls -la work/3m/4n5o6p*/
```

Zobaczysz:

- **Pliki .command.\***: Skrypty wykonania Nextflow'a i logi (jak wcześniej)
- **Przygotowane pliki wejściowe**: Dowiązania symboliczne do rzeczywistych plików wejściowych
- **Pliki wyjściowe**: Maski segmentacji, wyniki pośrednie itp.

Kluczowa różnica w stosunku do Hello World:

- Prawdziwe pipeline'y przygotowują duże pliki wejściowe (obrazy, dane referencyjne)
- Pliki wyjściowe mogą być dość duże (maski segmentacji, przetworzone obrazy)
- Wiele plików wejściowych i wyjściowych na zadanie

!!! Tip

    Jeśli proces się nie powiedzie, możesz przejść do jego katalogu roboczego, przejrzeć `.command.err` w poszukiwaniu komunikatów o błędach, a nawet ręcznie ponownie uruchomić `.command.sh`, aby debugować problem.

### 4.3. Czyszczenie katalogu roboczego

Katalog roboczy może stać się dość duży po wielu uruchomieniach pipeline'u.
Jak nauczyliśmy się w Części 1, możesz użyć `nextflow clean`, aby usunąć katalogi robocze ze starych uruchomień.

Jednak dla pipeline'ów nf-core z dużymi plikami pośrednimi szczególnie ważne jest regularne czyszczenie.

### Podsumowanie

Rozumiesz, jak pipeline'y nf-core organizują swoje katalogi robocze i jak przeglądać poszczególne zadania w celu debugowania.

### Co dalej?

Dowiedz się o pamięci podręcznej Nextflow'a i jak wznowić nieudane uruchomienia pipeline'u.

---

## 5. Wznowienie uruchomienia pipeline'u

Jedną z najpotężniejszych funkcji Nextflow'a jest możliwość wznowienia pipeline'u od punktu awarii.

### 5.1. Mechanizm pamięci podręcznej

Gdy uruchamiasz pipeline'a z `-resume`, Nextflow:

1. Sprawdza pamięć podręczną dla każdego zadania
2. Jeśli dane wejściowe, kod i parametry są identyczne, ponownie używa zbuforowanego wyniku
3. Ponownie uruchamia tylko zadania, które się zmieniły lub nie powiodły

Jest to niezbędne dla długo działających pipeline'ów, gdzie awarie mogą wystąpić późno w wykonaniu.

### 5.2. Wypróbowanie resume z molkart

Uruchom to samo polecenie ponownie, ale dodaj `-resume`:

```bash
nextflow run ./molkart \
  --input 'data/samplesheet.csv' \
  --mindagap_tilesize 90 \
  --mindagap_boxsize 7 \
  --mindagap_loopnum 100 \
  --clahe_pyramid_tile 368 \
  --segmentation_method "cellpose" \
  --outdir results \
  -resume
```

Powinieneś zobaczyć wyjście podobne do:

```console
executor >  local (0)
[1a/2b3c4d] NFCORE_MOLKART:MOLKART:MINDAGAP_MINDAGAP (mem_only)        [100%] 2 of 2, cached: 2 ✔
[5e/6f7g8h] NFCORE_MOLKART:MOLKART:CLAHE (mem_only)                    [100%] 2 of 2, cached: 2 ✔
[7f/8g9h0i] NFCORE_MOLKART:MOLKART:CREATE_STACK (mem_only)             [100%] 1 of 1, cached: 1 ✔
[9h/0i1j2k] NFCORE_MOLKART:MOLKART:MINDAGAP_DUPLICATEFINDER (mem_only) [100%] 1 of 1, cached: 1 ✔
[2k/3l4m5n] NFCORE_MOLKART:MOLKART:CELLPOSE (mem_only)                 [100%] 1 of 1, cached: 1 ✔
...
```

Zauważ `cached: 2` lub `cached: 1` dla każdego procesu - nic nie zostało ponownie wykonane!

### 5.3. Kiedy resume jest przydatne

Resume jest szczególnie wartościowe, gdy:

- Pipeline nie powiedzie się z powodu limitów zasobów (brak pamięci, przekroczony limit czasu)
- Musisz zmodyfikować procesy końcowe bez ponownego uruchamiania kroków początkowych
- Twoje połączenie sieciowe zostanie przerwane podczas pobierania danych
- Chcesz dodać dodatkowe wyniki bez powtarzania obliczeń

!!! Warning

    Resume działa tylko wtedy, gdy nie zmieniłeś danych wejściowych, kodu pipeline'u ani parametrów.
    Jeśli zmienisz którykolwiek z nich, Nextflow poprawnie ponownie uruchomi dotknięte zadania.

### Podsumowanie

Wiesz, jak używać `-resume`, aby efektywnie ponownie uruchamiać pipeline'y bez powtarzania udanych zadań.

### Co dalej?

Teraz, gdy potrafisz uruchomić nf-core/molkart z danymi testowymi, jesteś gotowy, aby nauczyć się, jak skonfigurować go dla Twoich własnych zestawów danych.
