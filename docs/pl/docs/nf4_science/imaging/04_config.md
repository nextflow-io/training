# Część 4: Konfiguracja

W częściach 1-3 nauczyliśmy się uruchamiać Nextflow, uruchamiać pipeline nf-core oraz zarządzać danymi wejściowymi przy użyciu plików parametrów i arkuszy próbek.
Teraz zbadamy, jak konfigurować pipeline'y dla różnych środowisk obliczeniowych przy użyciu **plików konfiguracyjnych** i **profili**.

## Cele nauczania

Pod koniec tej części będziesz w stanie:

- Zrozumieć, jak Nextflow rozwiązuje konfigurację z wielu źródeł
- Używać wbudowanych profili nf-core dla kontenerów i testowania
- Tworzyć niestandardowe profile dla różnych środowisk obliczeniowych
- Dostosowywać żądania zasobów przy użyciu etykiet procesów
- Zarządzać limitami zasobów w ograniczonych środowiskach
- Sprawdzać rozwiązaną konfigurację za pomocą `nextflow config`

---

## 1. Zrozumienie konfiguracji Nextflow

### 1.1. Czym jest plik konfiguracyjny?

Nextflow używa plików konfiguracyjnych do oddzielenia **logiki workflow'a** (co robić) od **ustawień wykonania** (jak i gdzie to robić).

Pliki konfiguracyjne kontrolują:

- Silniki kontenerów (Docker, Singularity, Conda)
- Zasoby obliczeniowe (procesory, pamięć, czas)
- Platformy wykonawcze (lokalne, HPC, chmura)
- Parametry pipeline'u

### 1.2. Pierwszeństwo konfiguracji

Nextflow wczytuje konfigurację z wielu źródeł, przy czym późniejsze źródła nadpisują wcześniejsze:

1. **Konfiguracja pipeline'u**: `nextflow.config` w repozytorium pipeline'u
2. **Konfiguracja katalogu**: `nextflow.config` w Twoim bieżącym katalogu roboczym
3. **Konfiguracja użytkownika**: `~/.nextflow/config`
4. **Wiersz poleceń**: Parametry i opcje przekazane bezpośrednio

To warstwowe podejście pozwala zachować wartości domyślne w pipeline'ie, nadpisać je ustawieniami specyficznymi dla użytkownika i dokonać szybkich zmian w wierszu poleceń.

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

Zakomentujmy lub zmieńmy z powrotem linię `docker.enabled = true` z części 2 i zastanówmy się, jak możemy osiągnąć ten sam rezultat używając profilu w molkart.

---

## 2. Używanie profili

### 2.1. Czym są profile?

Profile to nazwane zestawy konfiguracji, które można aktywować za pomocą flagi `-profile` w poleceniu `nextflow run`.
Ułatwiają przełączanie między różnymi scenariuszami obliczeniowymi bez edytowania plików konfiguracyjnych.

Wszystkie pipeline'y nf-core zawierają szereg domyślnych profili, z których możemy skorzystać.

### 2.2. Sprawdzanie wbudowanych profili

Sprawdźmy je w pliku `molkart/nextflow.config` powiązanym z kodem pipeline'u:

```bash
code molkart/nextflow.config
```

Wyszukaj blok `profiles`:

```groovy title="molkart/nextflow.config (excerpt)"
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

- `docker`: Używa kontenerów Docker (najczęściej do lokalnego rozwoju)
- `singularity`: Używa Singularity/Apptainer (często na HPC)
- `conda`: Używa środowisk Conda
- `apptainer`: Używa kontenerów Apptainer

### 2.3. Ponowne uruchomienie z profilami zamiast nextflow.config

Teraz, gdy wyłączyliśmy konfigurację dockera w naszym lokalnym pliku `nextflow.config` i rozumiemy profile, uruchommy ponownie pipeline używając flagi `-profile`.

Wcześniej w części 3 stworzyliśmy plik `params.yaml` z naszymi niestandardowymi parametrami.
Możemy teraz połączyć go z wbudowanym profilem Docker:

```bash
nextflow run ./molkart \
  -profile docker \
  -params-file params.yaml \
  -resume
```

Przeanalizujmy, co robi każda flaga:

- `-profile docker`: Aktywuje profil Docker z `nextflow.config` molkart, który ustawia `docker.enabled = true`
- `-params-file params.yaml`: Wczytuje wszystkie parametry pipeline'u z naszego pliku YAML
- `-resume`: Ponownie używa wyników z pamięci podręcznej z poprzednich uruchomień

Ponieważ używamy `-resume`, Nextflow sprawdzi, czy coś się zmieniło od ostatniego uruchomienia.
Jeśli parametry, dane wejściowe i kod są takie same, wszystkie zadania zostaną pobrane z pamięci podręcznej, a pipeline zakończy się niemal natychmiast.

```console title="Output (excerpt)"
executor >  local (12)
...
[1a/2b3c4d] NFCORE_MOLKART:MOLKART:MINDAGAP_MINDAGAP (mem_only)   [100%] 2 of 2, cached: 2 ✔
[5e/6f7g8h] NFCORE_MOLKART:MOLKART:CLAHE (mem_only)               [100%] 2 of 2, cached: 2 ✔
...
-[nf-core/molkart] Pipeline completed successfully-
```

Zauważ, że wszystkie procesy pokazują `cached: 2` lub `cached: 1` - nic nie zostało ponownie wykonane!

### 2.4. Profile testowe

Profile testowe zapewniają szybkie sposoby określenia domyślnych parametrów wejściowych i plików danych, aby umożliwić Ci weryfikację działania pipeline'u.
Pipeline'y nf-core zawsze zawierają co najmniej dwa profile testowe:

- `test`: Mały zestaw danych z szybkimi parametrami do szybkiego testowania
- `test_full`: Bardziej kompleksowy test z większymi danymi

Przyjrzyjmy się bliżej profilowi `test` w molkart, który jest dołączany za pomocą dyrektywy `includeConfig`:

```groovy title="molkart/nextflow.config (excerpt)"
profiles {
  ...
    test      { includeConfig 'conf/test.config'      }
}
```

Oznacza to, że za każdym razem, gdy uruchamiamy pipeline z `-profile test`, Nextflow wczyta konfigurację z `conf/test.config`.

```groovy title="molkart/conf/test.config (excerpt)"
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

Możesz aktywować wiele profili, oddzielając je przecinkami.
Użyjmy tego do przetestowania naszego pipeline'u bez potrzeby używania pliku parametrów:

```bash
nextflow run ./molkart -profile docker,test --outdir results -resume
```

To łączy:

- `docker`: Włącza kontenery Docker
- `test`: Używa zestawu danych testowych i parametrów

Profile są stosowane od lewej do prawej, więc późniejsze profile nadpisują wcześniejsze, jeśli ustawiają te same wartości.

### Podsumowanie

Pipeline'y nf-core zawierają wbudowane profile dla kontenerów, testowania i specjalnych środowisk.
Możesz łączyć wiele profili, aby zbudować potrzebną Ci konfigurację.

### Co dalej?

Naucz się tworzyć własne niestandardowe profile dla różnych środowisk obliczeniowych.

---

## 3. Tworzenie niestandardowych profili

### 3.1. Tworzenie profili do przełączania między lokalnym rozwojem a wykonaniem na HPC

Stwórzmy niestandardowe profile dla dwóch scenariuszy:

1. Lokalny rozwój z Dockerem
2. Uniwersytecki HPC z harmonogramem Slurm i Singularity

Dodaj następujący kod do swojego `nextflow.config`:

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
# Do lokalnego rozwoju
nextflow run ./molkart -profile local_dev --input data/samplesheet.csv --outdir results

# Do HPC (gdy dostępne)
nextflow run ./molkart -profile hpc_cluster --input data/samplesheet.csv --outdir results
```

!!! note "Uwaga"

    Nie możemy przetestować profilu HPC w tym środowisku szkoleniowym, ponieważ nie mamy dostępu do harmonogramu Slurm.
    Ale to pokazuje, jak skonfigurowałbyś go do rzeczywistego użycia.

### 3.2. Używanie `nextflow config` do sprawdzania konfiguracji

Polecenie `nextflow config` pokazuje w pełni rozwiązaną konfigurację bez uruchamiania pipeline'u.

Wyświetl domyślną konfigurację:

```bash
nextflow config ./molkart
```

Wyświetl konfigurację z określonym profilem:

```bash
nextflow config -profile local_dev ./molkart
```

Jest to niezwykle przydatne do:

- Debugowania problemów z konfiguracją
- Zrozumienia, jakie wartości faktycznie zostaną użyte
- Sprawdzania, jak wiele profili ze sobą współdziała

### Podsumowanie

Niestandardowe profile pozwalają przełączać się między różnymi środowiskami obliczeniowymi za pomocą jednej flagi wiersza poleceń.
Używaj `nextflow config` do sprawdzania rozwiązanej konfiguracji przed uruchomieniem.

### Co dalej?

Naucz się dostosowywać żądania zasobów dla poszczególnych procesów przy użyciu systemu etykiet procesów nf-core.

---

## 4. Dostosowywanie żądań zasobów

### 4.1. Zrozumienie etykiet procesów w pipeline'ach nf-core

Dla uproszczenia pipeline'y nf-core używają [**etykiet procesów**](https://www.nextflow.io/docs/latest/reference/process.html#process-label) do standaryzacji alokacji zasobów we wszystkich pipeline'ach.
Każdy proces jest oznaczony etykietą taką jak `process_low`, `process_medium` lub `process_high`, aby opisać odpowiednio niskie, średnie lub wysokie wymagania dotyczące zasobów obliczeniowych.
Te etykiety są konwertowane na konkretne żądania zasobów w jednym z plików konfiguracyjnych znajdujących się w katalogu `conf/` pipeline'u.

```groovy title="molkart/conf/base.config (excerpt)"
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

Zauważ mnożnik `task.attempt` - pozwala to kolejnym ponownym próbom zadania żądać więcej zasobów, jeśli pipeline jest ustawiony z `process.maxRetries > 1`.

### 4.2. Nadpisywanie zasobów dla konkretnych procesów

Aby uzyskać szczegółową kontrolę, kieruj się do poszczególnych procesów po nazwie:

```groovy title="nextflow.config"
process {
    withName: 'NFCORE_MOLKART:MOLKART:CELLPOSE' {
        cpus   = 16
        memory = 32.GB
    }
}
```

Jeśli spróbujemy uruchomić ten pipeline z powyższym nadpisaniem, proces `CELLPOSE` zażąda 16 procesorów i 32 GB pamięci zamiast wartości domyślnej zdefiniowanej przez jego etykietę.
Spowoduje to niepowodzenie pipeline'u w naszym obecnym środowisku, ponieważ nie mamy dostępnej takiej ilości pamięci RAM.
W następnej sekcji nauczymy się, jak zapobiegać tego typu awariom.

!!! tip "Wskazówka"

    Aby znaleźć nazwy procesów, sprawdź dane wyjściowe wykonania pipeline'u lub sprawdź `.nextflow.log`.
    Nazwy procesów mają wzór `WORKFLOW:SUBWORKFLOW:PROCESS`.

### Podsumowanie

Pipeline'y nf-core używają etykiet procesów do standaryzacji alokacji zasobów.
Możesz nadpisać zasoby według etykiety (wpływa na wiele procesów) lub według nazwy (wpływa na jeden konkretny proces).

### Co dalej?

Naucz się zarządzać limitami zasobów w ograniczonych środowiskach, takich jak GitHub Codespaces.

---

## 5. Zarządzanie zasobami w ograniczonych środowiskach

### 5.1. Problem limitów zasobów

Gdybyśmy spróbowali uruchomić molkart z procesem żądającym 16 procesorów i 32 GB pamięci (jak pokazano w sekcji 4.2), zakończyłoby się to niepowodzeniem w naszym obecnym środowisku, ponieważ nie mamy dostępnych tylu zasobów.
W środowisku klastra z większymi węzłami takie żądania zostałyby przesłane do harmonogramu.

W ograniczonych środowiskach, takich jak GitHub Codespaces, bez limitów Nextflow odmówiłby uruchomienia procesów przekraczających dostępne zasoby.

### 5.2. Ustawianie limitów zasobów

Dyrektywa `resourceLimits` ogranicza żądania zasobów do określonych wartości:

```groovy title="nextflow.config"
process {
    resourceLimits = [ cpus: 2, memory: 7.GB ]
}
```

To mówi Nextflow'owi: "Jeśli jakikolwiek proces zażąda więcej niż 2 procesory lub 7 GB pamięci, ogranicz to do tych limitów."

### 5.3. Dodawanie limitów zasobów do niestandardowych profili

Zaktualizuj swoje niestandardowe profile, aby zawierały odpowiednie limity:

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

    Ustawienie zbyt niskich limitów zasobów może spowodować niepowodzenie procesów lub ich wolne działanie.
    Pipeline może potrzebować użyć mniej pamięciochłonnych algorytmów lub przetwarzać dane w mniejszych porcjach.

### Podsumowanie

Używaj `resourceLimits`, aby uruchamiać pipeline'y w środowiskach o ograniczonych zasobach, ograniczając żądania zasobów procesów.
Różne profile mogą mieć różne limity odpowiednie dla ich środowiska.

### Co dalej?

Ukończyłeś podstawowe szkolenie Nextflow dla Bioimaging!

---

## Podsumowanie

Teraz rozumiesz, jak konfigurować pipeline'y Nextflow dla różnych środowisk obliczeniowych.

Kluczowe umiejętności, których się nauczyłeś:

- **Pierwszeństwo konfiguracji**: Jak Nextflow rozwiązuje ustawienia z wielu źródeł
- **Profile nf-core**: Używanie wbudowanych profili dla kontenerów, testowania i narzędzi
- **Niestandardowe profile**: Tworzenie własnych profili dla różnych środowisk
- **Etykiety procesów**: Zrozumienie i nadpisywanie żądań zasobów według etykiety
- **Limity zasobów**: Zarządzanie ograniczonymi środowiskami za pomocą `resourceLimits`
- **Sprawdzanie konfiguracji**: Używanie `nextflow config` do debugowania i weryfikacji ustawień

Te umiejętności konfiguracyjne są przydatne w każdym pipeline'ie Nextflow i pomogą Ci efektywnie uruchamiać workflow'y na komputerach lokalnych, klastrach HPC i platformach chmurowych.

### Co dalej?

Gratulacje ukończenia kursu Nextflow dla Bioimaging!

Następne kroki:

- Wypełnij ankietę kursu, aby przekazać opinię
- Sprawdź [Hello Nextflow](../hello_nextflow/index.md), aby dowiedzieć się więcej o tworzeniu workflow'ów
- Poznaj [Hello nf-core](../hello_nf-core/index.md), aby zgłębić narzędzia nf-core
- Przeglądaj inne kursy w [kolekcjach szkoleniowych](../training_collections/index.md)
