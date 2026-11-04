# Część 2: Zarządzanie zasobami obliczeniowymi i błędami

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Tłumaczenie wspomagane przez AI - [dowiedz się więcej i zasugeruj ulepszenia](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

W [Części 1](./01_packaging_and_execution.md) dostosowałeś miejsce i sposób uruchamiania zadań pipeline'u.
Tutaj dowiesz się, ile zasobów obliczeniowych przydzielić każdemu zadaniu i co się dzieje, gdy zadanie kończy się błędem mimo Twoich najlepszych szacunków dotyczących alokacji.

---

## 1. Kontrola alokacji zasobów obliczeniowych

Domyślnie Nextflow przydziela jeden CPU każdemu procesowi za pomocą dyrektywy `cpus` i nie narzuca limitu pamięci, chyba że go ustawisz:

```groovy title="Built-in configuration"
process {
    cpus = 1
}
```

Z modułu [Nextflow Run](../nextflow_run/index.md) wiesz już, że konfiguracja tego pipeline'u ustawia `memory` na 1 GB dla wszystkich procesów.
Skąd jednak wiadomo, jakich wartości faktycznie użyć we własnych pipeline'ach?

### 1.1. Generowanie raportu wykorzystania zasobów

W module [Nextflow Run](../nextflow_run/02_configure_pipeline.md) wygenerowałeś już raport wykonania za pomocą `-with-report`.
Ten sam raport pozwala sprawdzić, ile CPU i pamięci faktycznie potrzebują Twoje procesy: uruchom workflow z domyślnymi alokacjami, zarejestruj rzeczywiste użycie, a następnie dostosuj wartości.

```bash
nextflow run main.nf -with-report report-config-1.html
```

Raport to plik HTML, który możesz otworzyć w przeglądarce.
Zawiera zestawienie czasu wykonania i wykorzystania zasobów dla każdego procesu, w tym jaki procent przydzielonych zasobów został faktycznie użyty.
Oto co pokazuje dla `cowpy` przy obecnych ustawieniach domyślnych (1 CPU, 1 GB pamięci):

| Metryka                  | Wartość |
| ------------------------ | ------- |
| Użycie CPU               | 116%    |
| Szczytowe użycie pamięci | 6,4 MB  |
| Przydzielona pamięć      | 1 GB    |

`cowpy` używa znacznie poniżej 1% swojej alokacji 1 GB; wartość `%cpu` powyżej 100% oznacza jedynie, że przez krótkie chwile proces korzysta z więcej niż jednego CPU wewnątrz kontenera.

Pełną listę dostępnych funkcji znajdziesz w dokumentacji [Reports](https://nextflow.io/docs/latest/reports.html).

### 1.2. Ustawianie alokacji zasobów dla konkretnego procesu

Powyższy raport pokazuje, że `cowpy` wygodnie mieści się w swojej bieżącej alokacji, ale załóżmy, że chcesz dać mu więcej zapasu — na przykład dlatego, że spodziewasz się większych danych wejściowych w środowisku produkcyjnym.
Możesz nadpisać wartości domyślne dla pojedynczego procesu za pomocą `withName`.

=== "Po"

    ```groovy title="nextflow.config" linenums="6" hl_lines="5-6"
    process {
        memory = 1.GB
        withName: 'cowpy' {
            conda = 'conda-forge::cowpy==1.1.5'
            memory = 2.GB
            cpus = 2
        }
    }
    ```

=== "Przed"

    ```groovy title="nextflow.config" linenums="6"
    process {
        memory = 1.GB
        withName: 'cowpy' {
            conda = 'conda-forge::cowpy==1.1.5'
        }
    }
    ```

Po wprowadzeniu tej zmiany każdy proces żąda 1 GB pamięci i jednego CPU, z wyjątkiem `cowpy`, który żąda 2 GB i 2 CPU (oprócz ustawienia `conda` z [Części 1](./01_packaging_and_execution.md)).

!!! info "Info"

    Jeśli Twój komputer ma niewiele CPU, a Ty przydzielasz ich dużo na proces, wywołania zadań mogą ustawiać się w kolejce, ponieważ Nextflow nie zażąda więcej CPU niż jest dostępnych.

Uruchom pipeline ponownie z inną nazwą pliku raportu, żeby móc porównać wyniki przed i po zmianie.

```bash
nextflow run main.nf -with-report report-config-2.html
```

??? success "Wyjście polecenia"

    ```console
    N E X T F L O W   ~  version 26.04.4

    Launching `main.nf` [voluminous_venter] revision: c3c85dec78

    executor >  local (8)
    [a1/0e96d4] sayHello (1)       | 3 of 3 ✔
    [a3/7173a3] convertToUpper (2) | 3 of 3 ✔
    [4f/a8ae3d] collectGreetings   | 1 of 1 ✔
    [91/3724f8] cowpy              | 1 of 1 ✔

    Outputs:

      /workspaces/training/config-exec/results

      first_output:
        - full_pipeline/intermediates/Bonjour-output.txt
        - full_pipeline/intermediates/Hola-output.txt
        - full_pipeline/intermediates/Hello-output.txt

      uppercased:
        - full_pipeline/intermediates/UPPER-Bonjour-output.txt
        - full_pipeline/intermediates/UPPER-Hello-output.txt
        - full_pipeline/intermediates/UPPER-Hola-output.txt

      collected: full_pipeline/intermediates/COLLECTED-batch-output.txt

      batch_report: full_pipeline/batch-report.txt

      cowpy_art: full_pipeline/cowpy-COLLECTED-batch-output.txt
    ```

Porównanie raportów dla `cowpy`:

| Metryka                  | Przed (1 CPU, 1 GB) | Po (2 CPU, 2 GB) |
| ------------------------ | ------------------- | ---------------- |
| Szczytowe użycie pamięci | 6,4 MB              | 6,4 MB           |
| Użycie CPU               | 116%                | 118%             |

Podwojenie alokacji nie zmieniło rzeczywistego użycia zasobów, co oznacza, że oryginalne 1 GB / 1 CPU było już hojne dla tego przykładowego obciążenia.
W prawdziwym pipeline'ie przetwarzającym nietrywalne dane liczby będą się znacząco różnić między procesami — właśnie dlatego warto profilować przed podjęciem decyzji o alokacji, zamiast zgadywać.

### 1.3. Dodawanie limitów zasobów

W zależności od infrastruktury obliczeniowej mogą obowiązywać twarde ograniczenia dotyczące tego, czego można żądać — na przykład limit obowiązujący w całym klastrze.
Dyrektywa `resourceLimits` pozwala ustawić takie limity:

```groovy title="Syntax example"
process {
    resourceLimits = [
        memory: 750.GB,
        cpus: 200,
        time: 30.d
    ]
}
```

Nextflow tłumaczy te wartości na format oczekiwany przez docelowy executor.
Jeśli proces żąda więcej niż wynosi limit, żądanie zostaje obcięte, a nie odrzucone.

!!! warning "Ostrzeżenie"

    Tej funkcji nie można przetestować w środowisku szkoleniowym, ponieważ do działania wymaga infrastruktury HPC.

??? info "Instytucjonalne konfiguracje referencyjne"

    Projekt nf-core utrzymuje [kolekcję plików konfiguracyjnych](https://nf-co.re/configs/) udostępnianych przez instytucje z całego świata, obejmujących szeroki zakres executorów HPC i chmurowych.
    Stanowią one dobry punkt wyjścia niezależnie od tego, czy Twoja instytucja znajduje się wśród nich.

### Podsumowanie

Wiesz już, jak generować raport profilowania w celu oceny wykorzystania zasobów, nadpisywać alokacje zasobów dla konkretnego procesu oraz ograniczać alokacje za pomocą `resourceLimits`.

### Co dalej?

Dowiedz się, jak sprawić, żeby pipeline automatycznie odtwarzał się po błędzie zadania — niezależnie od tego, czy Twoje szacunki dotyczące zasobów były trafne.

---

## 2. Obsługa błędów zadań za pomocą ponownych prób

Profilowanie mówi Ci, czego proces potrzebuje przez większość czasu, ale rzeczywiste obciążenia są zmienne: alokacja wystarczająca dla większości danych wejściowych może okazać się zbyt mała dla wyjątkowo dużych, a szacunki bywają po prostu błędne.
Zamiast pozwolić, żeby jedno nieudane zadanie zatrzymało całe uruchomienie, Nextflow może automatycznie ponowić próbę wykonania zadania, opcjonalnie przydzielając mu więcej zasobów przy każdej kolejnej próbie.

### 2.1. Automatyczne ponawianie nieudanego zadania

Żeby zobaczyć to w działaniu, celowo ustaw alokację pamięci dla `cowpy` poniżej faktycznych potrzeb: z [punktu 1.1](#11-generate-a-resource-utilization-report) wiesz, że szczyt użycia wynosi około 6,4 MB, więc 6 MB powinno być tuż poniżej wymaganego minimum.

=== "Po"

    ```groovy title="nextflow.config" linenums="6" hl_lines="5-7"
    process {
        memory = 1.GB
        withName: 'cowpy' {
            conda = 'conda-forge::cowpy==1.1.5'
            memory = 6.MB
            errorStrategy = 'retry'
            maxRetries = 2
        }
    }
    ```

=== "Przed"

    ```groovy title="nextflow.config" linenums="6"
    process {
        memory = 1.GB
        withName: 'cowpy' {
            conda = 'conda-forge::cowpy==1.1.5'
            memory = 2.GB
            cpus = 2
        }
    }
    ```

`errorStrategy` mówi Nextflow'owi, co zrobić, gdy zadanie zakończy się błędem: wartość `'retry'` powoduje ponowne przesłanie zadania zamiast zatrzymania całego pipeline'u.
`maxRetries` ogranicza liczbę dodatkowych prób przed ostatecznym poddaniem się przez Nextflow'a.

```bash
nextflow run main.nf
```

??? failure "Wynik polecenia (skrócony)"

    ```console
    N E X T F L O W   ~  version 26.04.4

    Launching `main.nf` [desperate_brazil] revision: c3c85dec78

    executor >  local (10)
    [67/fe1f49] sayHello (1)       | 3 of 3 ✔
    [8a/f13335] convertToUpper (1) | 3 of 3 ✔
    [39/1b24ed] collectGreetings   | 1 of 1 ✔
    [7a/d5eb6f] cowpy              | 0 of 1, retries: 2 ✘
    [9d/b79eb3] NOTE: Process `cowpy` terminated with an error exit status (137) -- Execution is retried (1)
    [6d/1d9d84] NOTE: Process `cowpy` terminated with an error exit status (137) -- Execution is retried (2)
    ERROR ~ Error executing process > 'cowpy'

    Caused by:
      Process `cowpy` terminated with an error exit status (137)

    Command executed:
      cat COLLECTED-batch-output.txt | cowpy -c "turkey" > cowpy-COLLECTED-batch-output.txt

    Command exit status:
      137

    Command output:
      (empty)

    Command error:
      /usr/local/bin/_activate_current_env.sh: line 35:    14 Killed                  micromamba activate "${ENV_NAME:-base}"

    Work dir:
      /workspaces/training/config-exec/work/7a/d5eb6feeac0eed18d95d3da7a7aeb4

    Tip: when you have fixed the problem you can continue the execution adding the option `-resume` to the run command line

    -- Check '.nextflow.log' file for details
    ```

Kod wyjścia 137 to standardowy sygnał zabicia procesu z powodu braku pamięci: kontener nie miał wystarczająco dużo pamięci, żeby w ogóle uruchomić `cowpy`.
Nextflow ponowił próbę dwukrotnie — łącznie trzy próby, zgodnie z `maxRetries = 2`.
Ponieważ alokacja pamięci nie zmieniała się między próbami, każda z nich napotykała tę samą przeszkodę. Po wyczerpaniu prób Nextflow zgłasza błąd i zatrzymuje pipeline, kończąc działanie z niezerowym kodem wyjścia.

Samo ponawianie prób nic nie naprawia, jeśli przyczyna błędu nie zmienia się między kolejnymi próbami.

### 2.2. Zwiększanie zasobów przy każdej ponownej próbie

Wewnątrz dyrektywy procesu `task.attempt` przechowuje numer bieżącej próby, zaczynając od 1.
Możesz użyć go w domknięciu, żeby skalować alokację zasobów w górę przy każdej ponownej próbie.

=== "Po"

    ```groovy title="nextflow.config" linenums="6" hl_lines="5 7"
    process {
        memory = 1.GB
        withName: 'cowpy' {
            conda = 'conda-forge::cowpy==1.1.5'
            memory = { 6.MB * task.attempt }
            errorStrategy = 'retry'
            maxRetries = 3
        }
    }
    ```

=== "Przed"

    ```groovy title="nextflow.config" linenums="6"
    process {
        memory = 1.GB
        withName: 'cowpy' {
            conda = 'conda-forge::cowpy==1.1.5'
            memory = 6.MB
            errorStrategy = 'retry'
            maxRetries = 2
        }
    }
    ```

Uruchom workflow ponownie:

```bash
nextflow run main.nf
```

??? success "Wynik polecenia (skrócony)"

    ```console
    N E X T F L O W   ~  version 26.04.4

    Launching `main.nf` [grave_joliot] revision: c3c85dec78

    executor >  local (9)
    [0f/211b8a] sayHello (2)       | 3 of 3 ✔
    [26/301ea2] convertToUpper (3) | 3 of 3 ✔
    [22/5a895b] collectGreetings   | 1 of 1 ✔
    [e1/beee86] cowpy              | 1 of 1, retries: 1 ✔
    [b6/7aed6a] NOTE: Process `cowpy` terminated with an error exit status (137) -- Execution is retried (1)

    Outputs:

      /workspaces/training/config-exec/results

      first_output:
        - full_pipeline/intermediates/Hola-output.txt
        - full_pipeline/intermediates/Hello-output.txt
        - full_pipeline/intermediates/Bonjour-output.txt

      uppercased:
        - full_pipeline/intermediates/UPPER-Hola-output.txt
        - full_pipeline/intermediates/UPPER-Hello-output.txt
        - full_pipeline/intermediates/UPPER-Bonjour-output.txt

      collected: full_pipeline/intermediates/COLLECTED-batch-output.txt

      batch_report: full_pipeline/batch-report.txt

      cowpy_art: full_pipeline/cowpy-COLLECTED-batch-output.txt
    ```

Pierwsza próba nadal kończy się niepowodzeniem przy 6 MB, ale ponowna próba uruchamia się z 12 MB (`6.MB * 2`) i kończy się sukcesem — pipeline kończy działanie z opublikowanymi wszystkimi wynikami.

!!! warning "Ostrzeżenie"

    Wynik w konsoli nadal zawiera linię `NOTE:` informującą o nieudanej pierwszej próbie, mimo że pipeline jako całość zakończył się sukcesem: Nextflow rejestruje każdą ponowną próbę osobno, ale nieudana próba nie wpływa na ostateczny wynik.
    Sprawdź podsumowanie `Outputs:` lub kod wyjścia polecenia, żeby potwierdzić, czy uruchomienie faktycznie się powiodło.

Więcej zaawansowanych wzorców ponawiania prób, w tym skalowania w zależności od konkretnego błędu, znajdziesz w dokumentacji Nextflow'a: [Dynamic computing resources](https://nextflow.io/docs/latest/process.html#dynamic-task-resources).

### Podsumowanie

Wiesz już, jak sprawić, żeby pipeline automatycznie ponawiał nieudane zadania, oraz jak skalować alokacje zasobów przy każdej ponownej próbie za pomocą `task.attempt`.

### Co dalej?

Przejdź do [Części 3](./03_profiles.md), gdzie dowiesz się, jak grupować taką konfigurację w przełączalne profile.

---

## Podsumowanie

W tej części nauczyłeś się:

- Generować raport profilowania zasobów i ustawiać alokacje zasobów dla poszczególnych procesów
- Ograniczać żądania zasobów za pomocą `resourceLimits`
- Automatycznie ponawiać nieudane zadania za pomocą `errorStrategy` i `maxRetries`
- Skalować alokację zasobów w górę przy każdej ponownej próbie za pomocą `task.attempt`
