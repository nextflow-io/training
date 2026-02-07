# Część 4: Konfiguracja

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Tłumaczenie wspomagane przez AI - [dowiedz się więcej i zasugeruj ulepszenia](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

W częściach 1-3 nauczyliśmy się jak uruchamiać Nextflow'a, jak uruchomić pipeline nf-core oraz jak zarządzać wejściami za pomocą plików parametrów i arkuszy próbek.
Teraz zbadamy jak konfigurować pipeline'y dla różnych środowisk obliczeniowych używając **plików konfiguracyjnych** i **profili**.

## Cele nauczania

Pod koniec tej części będziesz potrafić:

- Zrozumieć jak Nextflow rozwiązuje konfigurację z wielu źródeł
- Używać wbudowanych profili nf-core dla kontenerów i testowania
- Tworzyć własne profile dla różnych środowisk obliczeniowych
- Dostosowywać żądania zasobów używając etykiet procesów
- Zarządzać limitami zasobów w ograniczonych środowiskach
- Sprawdzać rozwiązaną konfigurację za pomocą `nextflow config`

---

## 1. Zrozumienie konfiguracji Nextflow'a

### 1.1. Czym jest plik konfiguracyjny?

Nextflow używa plików konfiguracyjnych do oddzielenia **logiki workflow'u** (co robić) od **ustawień wykonania** (jak i gdzie to robić).

Pliki konfiguracyjne kontrolują:

- Silniki kontenerów (Docker, Singularity, Conda)
- Zasoby obliczeniowe (CPU, pamięć, czas)
- Platformy wykonania (lokalnie, HPC, chmura)
- Parametry pipeline'u

### 1.2. Kolejność pierwszeństwa konfiguracji

Nextflow wczytuje konfigurację z wielu źródeł, przy czym późniejsze źródła nadpisują wcześniejsze:

1. **Konfiguracja pipeline'u**: `nextflow.config` w repozytorium pipeline'u
2. **Konfiguracja katalogu**: `nextflow.config` w bieżącym katalogu roboczym
3. **Konfiguracja użytkownika**: `~/.nextflow/config`
4. **Linia poleceń**: Parametry i opcje przekazane bezpośrednio

To warstwowe podejście pozwala zachować wartości domyślne w pipeline'ie, nadpisać je ustawieniami użytkownika i dokonać szybkich dostosowań w linii poleceń.

### 1.3. Nasza obecna konfiguracja

Przyjrzyjmy się konfiguracji, której używaliśmy:

```groovy title="nextflow.config"
docker.enabled = true
process {
    resourceLimits = [
        cpus: 2,
        memory: '7.GB',
    ]
}

```

Zakomentujmy lub zmieńmy z powrotem linię `docker.enabled = true` z Części 2 i dowiedzmy się, jak możemy osiągnąć ten sam rezultat używając profilu w molkart.

---

## 2. Używanie profili

### 2.1. Czym są profile?

Profile to nazwane zestawy konfiguracji, które można aktywować za pomocą flagi `-profile` przy użyciu polecenia `nextflow run`.
Ułatwiają przełączanie między różnymi scenariuszami obliczeniowymi bez edycji plików konfiguracyjnych.

Wszystkie pipeline'y nf-core zawierają szereg domyślnych profili, z których możemy skorzystać.

### 2.2. Sprawdzanie wbudowanych profili

Sprawdźmy je w pliku `molkart/nextflow.config` powiązanym z bazą kodu pipeline'u:

```bash
code molkart/nextflow.config
```

Wyszukaj blok `profiles`:

```groovy title="molkart/nextflow.config (fragment)"
profiles {
    docker {
        docker.enabled          = true
        singularity.enabled     = false
        conda.enabled           = false
    }
    singularity {
        singularity.enabled     = true
        singularity.autoMounts  = true
        docker.enabled          = false
    }
    conda {
        conda.enabled           = true
        docker.enabled          = false
        conda.channels          = ['conda-forge', 'bioconda']
    }
}
```

Popularne profile kontenerów:

- `docker`: Używaj kontenerów Docker (najczęściej dla lokalnego rozwoju)
- `singularity`: Używaj Singularity/Apptainer (powszechne na HPC)
- `conda`: Używaj środowisk Conda
- `apptainer`: Używaj kontenerów Apptainer

### 2.3. Ponowne uruchomienie z profilami zamiast nextflow.config

Teraz, gdy wyłączyliśmy konfigurację Docker w naszym lokalnym pliku `nextflow.config` i rozumiemy profile, uruchommy ponownie pipeline używając flagi `-profile`.

Wcześniej w Części 3 utworzyliśmy plik `params.yaml` z naszymi własnymi parametrami.
Możemy teraz połączyć to z wbudowanym profilem Docker:

```bash
nextflow run ./molkart \
  -profile docker \
  -params-file params.yaml \
  -resume
```

Rozłóżmy to, co robi każda flaga:

- `-profile docker`: Aktywuje profil Docker z `nextflow.config` molkart, który ustawia `docker.enabled = true`
- `-params-file params.yaml`: Wczytuje wszystkie parametry pipeline'u z naszego pliku YAML
- `-resume`: Ponownie używa wyników z cache'u z poprzednich uruchomień

Ponieważ używamy `-resume`, Nextflow sprawdzi, czy coś się zmieniło od ostatniego uruchomienia.
Jeśli parametry, wejścia i kod są takie same, wszystkie zadania będą pobrane z cache'u i pipeline zakończy się niemal natychmiast.

```console title="Wyjście (fragment)"
executor >  local (12)
...
[1a/2b3c4d] NFCORE_MOLKART:MOLKART:MINDAGAP_MINDAGAP (mem_only)   [100%] 2 of 2, cached: 2 ✔
[5e/6f7g8h] NFCORE_MOLKART:MOLKART:CLAHE (mem_only)               [100%] 2 of 2, cached: 2 ✔
...
-[nf-core/molkart] Pipeline completed successfully-
```

Zauważ, że wszystkie procesy pokazują `cached: 2` lub `cached: 1` - nic nie zostało ponownie wykonane!

### 2.4. Profile testowe

Profile testowe zapewniają szybkie sposoby określania domyślnych parametrów wejściowych i plików danych, aby umożliwić weryfikację działania pipeline'u.
Pipeline'y nf-core zawsze zawierają co najmniej dwa profile testowe:

- `test`: Mały zestaw danych z szybkimi parametrami do szybkiego testowania
- `test_full`: Bardziej kompleksowy test z większymi danymi

Przyjrzyjmy się bliżej profilowi `test` w molkart, który jest dołączany za pomocą dyrektywy `includeConfig`:

```groovy title="molkart/nextflow.config (fragment)"
profiles {
  ...
    test      { includeConfig 'conf/test.config'      }
}
```

Oznacza to, że za każdym razem gdy uruchamiamy pipeline z `-profile test`, Nextflow wczyta konfigurację z `conf/test.config`.

```groovy title="molkart/conf/test.config (fragment)"
params {
    config_profile_name        = 'Test profile'
    config_profile_description = 'Minimal test dataset to check pipeline function'

    input = 'https://raw.githubusercontent.com/nf-core/test-datasets/molkart/test_data/samplesheets/samplesheet_membrane.csv'
    mindagap_tilesize = 90
    mindagap_boxsize = 7
    mindagap_loopnum = 100
    clahe_pyramid_tile = 368
    segmentation_method = "mesmer,cellpose,stardist"
}

process {
    resourceLimits = [
        cpus: 4,
        memory: '15.GB',
        time: '1.h'
    ]
}
```

Zauważ, że ten profil zawiera te same parametry, których użyliśmy wcześniej w naszym pliku `params.yaml`.

Możesz aktywować wiele profili oddzielając je przecinkami.
Użyjmy tego do przetestowania naszego pipeline'u bez potrzeby używania pliku params:

```bash
nextflow run ./molkart -profile docker,test --outdir results -resume
```

To łączy:

- `docker`: Włącz kontenery Docker
- `test`: Użyj zestawu danych testowych i parametrów

Profile są stosowane od lewej do prawej, więc późniejsze profile nadpisują wcześniejsze, jeśli ustawiają te same wartości.

### Podsumowanie

Pipeline'y nf-core zawierają wbudowane profile dla kontenerów, testowania i specjalnych środowisk.
Możesz łączyć wiele profili, aby zbudować potrzebną konfigurację.

### Co dalej?

Naucz się jak tworzyć własne profile dla różnych środowisk obliczeniowych.

---

## 3. Tworzenie własnych profili

### 3.1. Tworzenie profili do przełączania między lokalnym rozwojem a wykonaniem na HPC

Stwórzmy własne profile dla dwóch scenariuszy:

1. Lokalny rozwój z Docker
2. Uniwersytecki HPC ze schedulerem Slurm i Singularity

Dodaj następujący kod do Swojego `nextflow.config`:

```groovy title="nextflow.config"
profiles {
    local_dev {
        docker.enabled = true
        process.executor = 'local'
    }

    hpc_cluster {
        singularity.enabled = true
        process.executor = 'slurm'
        process.queue = 'standard_queue'
        singularity.cacheDir = '/shared/containers'
    }
}
```

Teraz możesz łatwo przełączać się między środowiskami:

```bash
# Dla lokalnego rozwoju
nextflow run ./molkart -profile local_dev --input data/samplesheet.csv --outdir results

# Dla HPC (gdy jest dostępny)
nextflow run ./molkart -profile hpc_cluster --input data/samplesheet.csv --outdir results
```

!!! note "Uwaga"

    Nie możemy przetestować profilu HPC w tym środowisku szkoleniowym, ponieważ nie mamy dostępu do schedulera Slurm.
    Ale to pokazuje, jak skonfigurowałbyś go do rzeczywistego użycia.

### 3.2. Użycie `nextflow config` do sprawdzenia konfiguracji

Polecenie `nextflow config` pokazuje w pełni rozwiązaną konfigurację bez uruchamiania pipeline'u.

Wyświetl domyślną konfigurację:

```bash
nextflow config ./molkart
```

Wyświetl konfigurację z konkretnym profilem:

```bash
nextflow config -profile local_dev ./molkart
```

Jest to niezwykle przydatne do:

- Debugowania problemów z konfiguracją
- Zrozumienia, jakie wartości będą faktycznie użyte
- Sprawdzania, jak wiele profili wchodzi w interakcje

### Podsumowanie

Własne profile pozwalają na przełączanie między różnymi środowiskami obliczeniowymi za pomocą pojedynczej flagi linii poleceń.
Użyj `nextflow config` do sprawdzenia rozwiązanej konfiguracji przed uruchomieniem.

### Co dalej?

Naucz się jak dostosować żądania zasobów dla poszczególnych procesów używając systemu etykiet procesów nf-core.

---

## 4. Dostosowywanie żądań zasobów

### 4.1. Zrozumienie etykiet procesów w pipeline'ach nf-core

Dla uproszczenia, pipeline'y nf-core używają [**etykiet procesów**](https://www.nextflow.io/docs/latest/reference/process.html#process-label) do standaryzacji alokacji zasobów we wszystkich pipeline'ach.
Każdy proces jest oznaczony etykietą taką jak `process_low`, `process_medium` lub `process_high` w celu opisania odpowiednio niskich, średnich lub wysokich wymagań zasobów obliczeniowych.
Te etykiety są konwertowane na konkretne żądania zasobów w jednym z plików konfiguracyjnych znajdujących się w katalogu `conf/` pipeline'u.

```groovy title="molkart/conf/base.config (fragment)"
process {
    cpus   = { 1      * task.attempt }
    memory = { 6.GB   * task.attempt }
    time   = { 4.h    * task.attempt }

    errorStrategy = { task.exitStatus in ((130..145) + 104) ? 'retry' : 'finish' }
    maxRetries    = 1
    maxErrors     = '-1'

    withLabel:process_single {
        cpus   = { 1                   }
        memory = { 6.GB * task.attempt }
        time   = { 4.h  * task.attempt }
    }
    withLabel:process_low {
        cpus   = { 2     * task.attempt }
        memory = { 12.GB * task.attempt }
        time   = { 4.h   * task.attempt }
    }
    withLabel:process_medium {
        cpus   = { 6     * task.attempt }
        memory = { 36.GB * task.attempt }
        time   = { 8.h   * task.attempt }
    }
    withLabel:process_high {
        cpus   = { 12    * task.attempt }
        memory = { 72.GB * task.attempt }
        time   = { 16.h  * task.attempt }
    }
}
```

Zauważ mnożnik `task.attempt` - pozwala to kolejnym ponownym próbom zadania zażądać więcej zasobów, jeśli pipeline jest ustawiony z `process.maxRetries > 1`.

### 4.2. Nadpisywanie zasobów dla konkretnych procesów

Aby uzyskać szczegółową kontrolę, wskaż poszczególne procesy po nazwie:

```groovy title="nextflow.config"
process {
    withName: 'NFCORE_MOLKART:MOLKART:CELLPOSE' {
        cpus   = 16
        memory = 32.GB
    }
}
```

Jeśli spróbujemy uruchomić ten pipeline z powyższym nadpisaniem, proces `CELLPOSE` zażąda 16 CPU i 32 GB pamięci zamiast wartości domyślnej zdefiniowanej przez jego etykietę.
Spowoduje to niepowodzenie pipeline'u w naszym obecnym środowisku, ponieważ nie mamy dostępnej takiej ilości pamięci RAM.
W następnej sekcji nauczymy się, jak zapobiegać tego typu awariom.

!!! tip "Wskazówka"

    Aby znaleźć nazwy procesów, sprawdź wyjście wykonania pipeline'u lub sprawdź `.nextflow.log`.
    Nazwy procesów mają wzór `WORKFLOW:SUBWORKFLOW:PROCESS`.

### Podsumowanie

Pipeline'y nf-core używają etykiet procesów do standaryzacji alokacji zasobów.
Możesz nadpisać zasoby według etykiety (wpływa na wiele procesów) lub według nazwy (wpływa na jeden konkretny proces).

### Co dalej?

Naucz się jak zarządzać limitami zasobów w ograniczonych środowiskach takich jak GitHub Codespaces.

---

## 5. Zarządzanie zasobami w ograniczonych środowiskach

### 5.1. Problem limitów zasobów

Gdybyśmy spróbowali uruchomić molkart z procesem żądającym 16 CPU i 32 GB pamięci (jak pokazano w sekcji 4.2), nie powiódłby się w naszym obecnym środowisku, ponieważ nie mamy dostępnych tylu zasobów.
W środowisku klastra z większymi węzłami takie żądania zostałyby przesłane do schedulera.

W ograniczonych środowiskach takich jak GitHub Codespaces, bez limitów, Nextflow odmówiłby uruchomienia procesów przekraczających dostępne zasoby.

### 5.2. Ustawianie limitów zasobów

Dyrektywa `resourceLimits` ogranicza żądania zasobów do określonych wartości:

```groovy title="nextflow.config"
process {
    resourceLimits = [ cpus: 2, memory: 7.GB ]
}
```

Mówi to Nextflow'owi: "Jeśli jakikolwiek proces zażąda więcej niż 2 CPU lub 7 GB pamięci, ogranicz to do tych limitów."

### 5.3. Dodawanie limitów zasobów do własnych profili

Zaktualizuj Swoje własne profile, aby zawierały odpowiednie limity:

```groovy title="nextflow.config"
profiles {
    local_dev {
        docker.enabled = true
        process.executor = 'local'
        process.resourceLimits = [
            cpus: 2,
            memory: 7.GB
        ]
    }

    hpc_cluster {
        singularity.enabled = true
        process.executor = 'slurm'
        process.queue = 'batch'
        process.resourceLimits = [
            cpus: 32,
            memory: 128.GB,
            time: 24.h
        ]
    }
}
```

!!! warning "Ostrzeżenie"

    Ustawienie limitów zasobów zbyt nisko może spowodować niepowodzenie procesów lub spowolnić ich działanie.
    Pipeline może potrzebować użyć mniej wymagających pamięciowo algorytmów lub przetwarzać dane w mniejszych fragmentach.

### Podsumowanie

Użyj `resourceLimits`, aby uruchamiać pipeline'y w środowiskach z ograniczonymi zasobami poprzez ograniczenie żądań zasobów procesów.
Różne profile mogą mieć różne limity odpowiednie dla Swojego środowiska.

### Co dalej?

Ukończyłeś podstawowe szkolenie Nextflow for Bioimaging!

---

## Podsumowanie

Teraz rozumiesz, jak konfigurować pipeline'y Nextflow'a dla różnych środowisk obliczeniowych.

Kluczowe umiejętności, których się nauczyłeś:

- **Kolejność pierwszeństwa konfiguracji**: Jak Nextflow rozwiązuje ustawienia z wielu źródeł
- **Profile nf-core**: Używanie wbudowanych profili dla kontenerów, testowania i narzędzi
- **Własne profile**: Tworzenie własnych profili dla różnych środowisk
- **Etykiety procesów**: Zrozumienie i nadpisywanie żądań zasobów według etykiety
- **Limity zasobów**: Zarządzanie ograniczonymi środowiskami za pomocą `resourceLimits`
- **Sprawdzanie konfiguracji**: Używanie `nextflow config` do debugowania i weryfikacji ustawień

Te umiejętności konfiguracyjne można przenieść na każdy pipeline Nextflow'a i pomogą Ci efektywnie uruchamiać workflow'e na lokalnych maszynach, klastrach HPC i platformach chmurowych.

### Co dalej?

Gratulacje ukończenia kursu Nextflow for Bioimaging!

Następne kroki:

- Wypełnij ankietę kursu, aby przekazać informację zwrotną
- Sprawdź [Hello Nextflow](../hello_nextflow/index.md), aby dowiedzieć się więcej o tworzeniu workflow'ów
- Poznaj [Hello nf-core](../hello_nf-core/index.md), aby zgłębić narzędzia nf-core
- Przeglądaj inne kursy w [kolekcjach szkoleniowych](../training_collections/index.md)
