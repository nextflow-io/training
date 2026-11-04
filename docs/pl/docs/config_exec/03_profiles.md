# Część 3: Używanie profili do przełączania konfiguracji

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Tłumaczenie wspomagane przez AI - [dowiedz się więcej i zasugeruj ulepszenia](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

W [Części 1](./01_packaging_and_execution.md) i [Części 2](./02_resources_and_retries.md) zebrałeś kilka opcji konfiguracyjnych: pakowanie oprogramowania, platformę wykonawczą i przydziały zasobów.
W praktyce często będziesz chciał przełączać się między całymi zestawami tych opcji w zależności od środowiska uruchomieniowego — na przykład laptopa do programowania i klastra HPC do produkcji.

Nextflow pozwala Ci skonfigurować dowolną liczbę [profili](https://nextflow.io/docs/latest/config.html#profiles) opisujących różne konfiguracje i wybrać jeden (lub kilka) w czasie wykonania za pomocą jednej flagi.

Jeden z nich już znasz: profil `test` z [Nextflow Run](../nextflow_run/index.md) nadpisuje parametry wejściowe na mały, dobrze zdefiniowany zestaw.
Teraz stworzysz własne profile infrastrukturalne i połączysz je z nim.

---

## 1. Tworzenie profili dla różnych środowisk

### 1.1. Konfiguracja profili

Dodaj dwa profile do `nextflow.config`: jeden do uruchamiania na zwykłym laptopie z Docker, a drugi dla uczelnianego klastra HPC z harmonogramem Slurm i Conda.

=== "Po"

    ```groovy title="nextflow.config" linenums="35" hl_lines="10-19"
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
        my_laptop {
            process.executor = 'local'
            docker.enabled = true
        }
        univ_hpc {
            process.executor = 'slurm'
            conda.enabled = true
            process.resourceLimits = [
                memory: 750.GB,
                cpus: 200,
                time: 30.d
            ]
        }
    }
    ```

=== "Przed"

    ```groovy title="nextflow.config" linenums="35"
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

Profil `univ_hpc` ustawia również limity zasobów, ponieważ jest to zazwyczaj wymagane we współdzielonej infrastrukturze HPC.

### 1.2. Uruchamianie workflow'u z profilem

Wybierz profil w czasie wykonania za pomocą `-profile`.

```bash
nextflow run main.nf -profile my_laptop
```

??? success "Wyjście polecenia"

    ```console
    N E X T F L O W   ~  version 26.04.4

    Launching `main.nf` [nice_kimura] revision: c3c85dec78

    executor >  local (8)
    [df/b86b2f] sayHello (2)       | 3 of 3 ✔
    [b2/e2cb9c] convertToUpper (3) | 3 of 3 ✔
    [af/2ecff7] collectGreetings   | 1 of 1 ✔
    [c6/713c96] cowpy              | 1 of 1 ✔

    Outputs:

      /workspaces/training/config-exec/results

      first_output:
        - full_pipeline/intermediates/Bonjour-output.txt
        - full_pipeline/intermediates/Hello-output.txt
        - full_pipeline/intermediates/Hola-output.txt

      uppercased:
        - full_pipeline/intermediates/UPPER-Bonjour-output.txt
        - full_pipeline/intermediates/UPPER-Hello-output.txt
        - full_pipeline/intermediates/UPPER-Hola-output.txt

      collected: full_pipeline/intermediates/COLLECTED-batch-output.txt

      batch_report: full_pipeline/batch-report.txt

      cowpy_art: full_pipeline/cowpy-COLLECTED-batch-output.txt
    ```

!!! warning "Ostrzeżenie"

    Profil `univ_hpc` nie zadziała w środowisku szkoleniowym, ponieważ nie ma dostępnego harmonogramu Slurm.

Jeśli znajdziesz inne ustawienia, które zawsze powinny być razem, dodaj je do odpowiedniego profilu.
Możesz też tworzyć dodatkowe profile grupujące dowolne inne kombinacje, których potrzebujesz.

### 1.3. Uruchamianie z wieloma profilami

Profile nie wykluczają się wzajemnie.
Możesz aktywować kilka naraz za pomocą `-profile <profil1>,<profil2>`.
Połącz `my_laptop` z profilem `test`, który już znasz z Nextflow Run.

```bash
nextflow run main.nf -profile my_laptop,test
```

??? success "Wyjście polecenia"

    ```console
    N E X T F L O W   ~  version 26.04.4

    Launching `main.nf` [romantic_volhard] revision: c3c85dec78

    executor >  local (8)
    [22/81de4f] sayHello (3)       | 3 of 3 ✔
    [93/e16b04] convertToUpper (1) | 3 of 3 ✔
    [ca/3e6a91] collectGreetings   | 1 of 1 ✔
    [16/0cc7e3] cowpy              | 1 of 1 ✔

    Outputs:

      /workspaces/training/config-exec/results

      first_output:
        - full_pipeline/intermediates/Hola-output.txt
        - full_pipeline/intermediates/Hello-output.txt
        - full_pipeline/intermediates/Bonjour-output.txt

      uppercased:
        - full_pipeline/intermediates/UPPER-Hello-output.txt
        - full_pipeline/intermediates/UPPER-Bonjour-output.txt
        - full_pipeline/intermediates/UPPER-Hola-output.txt

      collected: full_pipeline/intermediates/COLLECTED-test-output.txt

      batch_report: full_pipeline/test-report.txt

      cowpy_art: full_pipeline/cowpy-COLLECTED-test-output.txt
    ```

Nazwy poszczególnych plików poprawnie uwzględniają `batch = 'test'` z profilu `test` (`COLLECTED-test-output.txt` i tak dalej).

Jeśli połączysz profile ustawiające tę samą opcję, Nextflow rozwiązuje konflikt, używając wartości odczytanej jako ostatnia — czyli tej, która pojawia się później w pliku.
Gdy sprzeczne ustawienia pochodzą z zupełnie różnych źródeł konfiguracji, obowiązuje standardowa [kolejność pierwszeństwa](https://www.nextflow.io/docs/latest/config.html).

### Podsumowanie

Wiesz już, jak definiować profile grupujące konfigurację specyficzną dla infrastruktury, wybierać jeden z nich w czasie wykonania za pomocą `-profile`, łączyć wiele profili w jednym uruchomieniu oraz jak Nextflow rozwiązuje konflikty, gdy więcej niż jeden profil ustawia tę samą opcję.

### Co dalej?

Dowiedz się, jak sprawdzić w pełni rozwiązaną konfigurację przed uruchomieniem czegokolwiek.

---

## 2. Sprawdzanie rozwiązanej konfiguracji

Polecenia `nextflow config -profile test` użyłeś już w [Nextflow Run](../nextflow_run/02_configure_pipeline.md), aby sprawdzić, do czego sprowadza się pojedynczy profil.
Staje się ono szczególnie przydatne, gdy łączysz wiele profili: jak właśnie widziałeś, gdy dwa profile ustawiają tę samą opcję, ręczne ustalenie, która wartość faktycznie wygrywa, może być kłopotliwe.
Polecenie `nextflow config` rozwiązuje to wszystko za Ciebie, bez uruchamiania pipeline'u.

### 2.1. Rozwiązywanie domyślnej konfiguracji

```bash
nextflow config
```

??? success "Wyjście polecenia"

    ```groovy
    params {
       input = 'data/greetings.csv'
       batch = 'batch'
       character = 'turkey'
    }

    docker {
       enabled = true
    }

    process {
       memory = '1 GB'
       withName:cowpy {
          conda = 'conda-forge::cowpy==1.1.5'
          memory = '2 GB'
          cpus = 2
       }
    }
    ```

To dokładnie to, co zostałoby zastosowane, gdybyś uruchomił pipeline bez żadnych dodatkowych flag.

### 2.2. Rozwiązywanie konfiguracji z aktywnymi profilami

Dodaj te same profile, których użyłbyś do rzeczywistego uruchomienia.

```bash
nextflow config -profile my_laptop,test
```

??? success "Wyjście polecenia"

    ```groovy
    params {
       input = 'data/greetings.csv'
       batch = 'test'
       character = 'tux'
    }

    docker {
       enabled = true
    }

    process {
       memory = '1 GB'
       withName:cowpy {
          conda = 'conda-forge::cowpy==1.1.5'
          memory = '2 GB'
          cpus = 2
       }
       executor = 'local'
    }
    ```

Porównanie obu wyników potwierdza, co się zmieniło: `params.batch`, `params.character` i `process.executor` odzwierciedlają profile `my_laptop,test`.
Jest to szczególnie cenne w przypadku pipeline'ów z wieloma warstwami konfiguracji, gdzie ręczne ustalanie rozwiązanych ustawień byłoby żmudne i podatne na błędy.

### Podsumowanie

Wiesz już, jak używać `nextflow config` do sprawdzania w pełni rozwiązanej konfiguracji dla dowolnej kombinacji profili — przed uruchomieniem czegokolwiek.

### Co dalej?

Poznałeś podstawy konfigurowania pipeline'ów Nextflow'a.
Zajrzyj do [podsumowania kursu](next_steps.md), aby dowiedzieć się, co robić dalej.

---

## Podsumowanie

W tej części nauczyłeś się:

- Definiować profile grupujące konfigurację specyficzną dla infrastruktury
- Łączyć wiele profili w jednym uruchomieniu i rozumieć, jak rozwiązywane są konflikty między nimi
- Używać `nextflow config` do sprawdzania w pełni rozwiązanej konfiguracji
