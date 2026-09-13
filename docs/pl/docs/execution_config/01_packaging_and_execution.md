# Część 1: Dostosowanie do środowiska obliczeniowego

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Tłumaczenie wspomagane przez AI - [dowiedz się więcej i zasugeruj ulepszenia](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

W module [Nextflow Run](../nextflow_run/index.md) skonfigurowałeś wejścia, parametry i wyjścia pipeline'u.
Ten kurs obejmuje drugą połowę zagadnienia: dostosowanie wykonania pipeline'u do dowolnego środowiska obliczeniowego, bez modyfikowania kodu workflow'u.

!!! example "Scenariusz"

    Opracowałeś i przetestowałeś swój pipeline na laptopie przy użyciu Docker'a.
    Teraz musisz go przekazać dalej: współpracownik ma skonfigurowaną tylko Condę, a klaster HPC Twojej instytucji wymaga, żeby zadania przechodziły przez własny scheduler z własnymi limitami zasobów.
    Żadna z tych zmian nie powinna wymagać przepisywania samego pipeline'u.

Ten sam kod pipeline'u może działać we wszystkich tych miejscach, ponieważ żadna z tych kwestii nie jest wbudowana w workflow.
Pakowanie oprogramowania, platforma wykonawcza i alokacja zasobów są kontrolowane przez konfigurację nałożoną na kod — i właśnie to obejmuje ten kurs: jak dostosować ten sam pipeline do nowego środowiska, zmieniając konfigurację, a nie kod.

---

## 1. Wybór technologii pakowania oprogramowania

W module [Nextflow Run](../nextflow_run/index.md) widziałeś profil `conda` skonfigurowany w `nextflow.config` jako alternatywę dla Docker'a.
Tutaj samodzielnie zbudujesz ten sam przełącznik i zobaczysz, co jest potrzebne, żeby proces był faktycznie użyteczny z Condą.

### 1.1. Wyłączenie Docker'a i włączenie Condy

Zmień `docker.enabled` na `false` i dodaj dyrektywę włączającą Condę.

=== "Po"

    ```groovy title="nextflow.config" linenums="1" hl_lines="1-2"
    docker.enabled = false
    conda.enabled = true
    ```

=== "Przed"

    ```groovy title="nextflow.config" linenums="1"
    docker.enabled = true
    ```

Dzięki temu Nextflow może tworzyć i używać środowisk Conda dla każdego procesu, który ma określony pakiet Conda.
Proces `cowpy` jeszcze go nie ma, więc dodajmy go — wyłącznie z poziomu konfiguracji.

### 1.2. Dodanie pakietu Conda przez konfigurację

Dyrektywę `conda` można ustawić bezpośrednio w definicji procesu, tak jak dyrektywa `container` w `modules/cowpy.nf`, ale nie jest to konieczne: `withName` pozwala ustawić ją z poziomu konfiguracji, ograniczając zakres wyłącznie do procesu `cowpy`.

=== "Po"

    ```groovy title="nextflow.config" linenums="6" hl_lines="3-5"
    process {
        memory = 1.GB
        withName: 'cowpy' {
            conda = 'conda-forge::cowpy==1.1.5'
        }
    }
    ```

=== "Przed"

    ```groovy title="nextflow.config" linenums="6"
    process {
        memory = 1.GB
    }
    ```

To nie zastępuje dyrektywy `container` znajdującej się już w kodzie pipeline'u — dodaje alternatywę obok niej, bez żadnej ingerencji w ten kod.

!!! tip "Wskazówka"

    Wyszukiwarka [Seqera Containers](https://seqera.io/containers/) to wygodny sposób na znalezienie URI pakietu Conda dla danego narzędzia, nawet jeśli nie planujesz budować z niego kontenera.

### 1.3. Uruchomienie workflow'u w celu weryfikacji działania z Condą

```bash
nextflow run main.nf --batch conda
```

??? success "Wyjście polecenia"

    ```console
    N E X T F L O W   ~  version 26.04.4

    Launching `main.nf` [romantic_pike] revision: c3c85dec78

    executor >  local (8)
    [6d/d48030] sayHello (2)       | 3 of 3 ✔
    [f5/7a9d76] convertToUpper (1) | 3 of 3 ✔
    [1c/79b693] collectGreetings   | 1 of 1 ✔
    Creating env using conda: conda-forge::cowpy==1.1.5 [cache /workspaces/training/execution-config/work/conda/env-898314d566668b6587ad714ae06b8520]
    [bb/64b67c] cowpy              | 1 of 1 ✔

    Outputs:

      /workspaces/training/execution-config/results

      first_output:
        - full_pipeline/intermediates/Hola-output.txt
        - full_pipeline/intermediates/Bonjour-output.txt
        - full_pipeline/intermediates/Hello-output.txt

      uppercased:
        - full_pipeline/intermediates/UPPER-Hello-output.txt
        - full_pipeline/intermediates/UPPER-Bonjour-output.txt
        - full_pipeline/intermediates/UPPER-Hola-output.txt

      collected: full_pipeline/intermediates/COLLECTED-conda-output.txt

      batch_report: full_pipeline/conda-report.txt

      cowpy_art: full_pipeline/cowpy-COLLECTED-conda-output.txt
    ```

Wynik jest taki sam jak przy uruchomieniu z Docker'em, choć mechanizm działania jest inny: Nextflow pobiera pakiet Conda i buduje z niego środowisko, zamiast pobierać obraz kontenera.

!!! info "Info"

    Budowanie nowego środowiska Conda może trwać nieco dłużej niż pierwsze pobranie kontenera, ale pakiet użyty tutaj jest mały, więc powinno być szybko.

Teraz wróć do Docker'a na potrzeby pozostałej części kursu.

```groovy title="nextflow.config" linenums="1"
docker.enabled = true
```

??? tip "Mieszanie Docker'a i Condy"

    Ponieważ te ustawienia mają zakres per proces, można je mieszać: część procesów używa Docker'a, inne Condy — w zależności od tego, co jest dostępne dla danego narzędzia.
    Jeśli dla tego samego procesu ustawiona jest zarówno dyrektywa `container` (w kodzie pipeline'u), jak i dyrektywa `conda` (tutaj, z konfiguracji), a oba systemy pakowania są włączone, Nextflow priorytetyzuje kontenery.

### Podsumowanie

Wiesz już, jak skonfigurować technologię pakowania oprogramowania dla procesu i jak przełączać się między Docker'em a Condą.

### Co dalej?

Dowiedz się, jak zmienić platformę wykonawczą, której Nextflow używa do uruchamiania zadań.

---

## 2. Wybór platformy wykonawczej

Każdy pipeline, który do tej pory uruchamiałeś, korzystał z lokalnego executora: każde zadanie działa na tej samej maszynie co sam Nextflow.
Nextflow sprawdza dostępne procesory i pamięć, wstrzymując zadania do momentu zwolnienia wystarczających zasobów.

Lokalny executor jest wygodny, ale nie skaluje się poza pojedynczą maszynę.
Nextflow obsługuje [wiele innych backendów wykonawczych](https://nextflow.io/docs/latest/executor.html), w tym schedulery HPC (Slurm, LSF, SGE, PBS i inne) oraz platformy chmurowe (AWS Batch, Google Cloud Batch, Azure Batch, Kubernetes i inne).

### 2.1. Wybór innego backendu

Executor ustawia się za pomocą dyrektywy procesu o nazwie `executor`.
Domyślnie jest to `local`, więc poniższy zapis jest domyślnie przyjmowany:

```groovy title="Built-in configuration"
process {
    executor = 'local'
}
```

Aby wybrać inny backend, ustaw dyrektywę na żądany executor.

```groovy title="nextflow.config"
process {
    executor = 'slurm'
}
```

!!! warning "Ostrzeżenie"

    Środowisko szkoleniowe nie jest podłączone do klastra HPC, więc tego nie można tutaj uruchomić.

### 2.2. Abstrakcja składni specyficznej dla backendu

Większość platform HPC wymaga, żeby zgłoszenia zadań określały zapotrzebowanie na zasoby — takie jak procesory, pamięć i nazwa kolejki — używając własnej składni.
To samo zapotrzebowanie na 8 procesorów i 4 GB RAM w kolejce o nazwie `my-science-work` wygląda zupełnie inaczej w zależności od schedulera.

??? abstract "Przykłady"

    ```bash title="Config for SLURM / submit using sbatch"
    #SBATCH -o /path/to/my/task/directory/my-task-1.log
    #SBATCH --no-requeue
    #SBATCH -c 8
    #SBATCH --mem 4096M
    #SBATCH -p my-science-work
    ```

    ```bash title="Config for PBS / submit using qsub"
    #PBS -o /path/to/my/task/directory/my-task-1.log
    #PBS -j oe
    #PBS -q my-science-work
    #PBS -l nodes=1:ppn=8
    #PBS -l mem=4gb
    ```

    ```bash title="Config for SGE / submit using qsub"
    #$ -o /path/to/my/task/directory/my-task-1.log
    #$ -j y
    #$ -terse
    #$ -notify
    #$ -q my-science-work
    #$ -l slots=8
    #$ -l h_rss=4096M,mem_free=4096M
    ```

Nextflow abstrahuje to wszystko: raz określasz standardowe właściwości, takie jak `cpus`, `memory` i `queue` (pełna lista w dokumentacji [dyrektyw procesu](https://nextflow.io/docs/latest/reference/process.html#process-directives)), a Nextflow tłumaczy je na odpowiednie skrypty specyficzne dla danego backendu w czasie wykonania.

### 2.3. Podgląd tego, co Nextflow faktycznie uruchamia

To tłumaczenie to nie tylko wygoda na poziomie pliku konfiguracyjnego — stoi za nim coś konkretnego, co możesz sprawdzić już teraz, nawet przy lokalnym executorze.
W module [Nextflow Run, sekcja 1.3](../nextflow_run/01_run_nextflow.md#13-explore-the-work-directory) zajrzałeś do katalogu zadania w `work/` i znalazłeś `.command.sh` — dokładne polecenie uruchomione przez Nextflow.
Ten sam katalog zawiera też plik, którego jeszcze nie oglądałeś: `.command.run`.

```bash
cat work/0a/0df4a1*/.command.run
```

??? success "Wynik polecenia (fragment)"

    ```console
    #!/bin/bash
    ### ---
    ### name: 'convertToUpper (3)'
    ### container: 'null'
    ### outputs:
    ### - 'UPPER-Bonjour-output.txt'
    ### ...
    set -e
    set -u
    ...
    nxf_launch() {
        /bin/bash -ue /workspaces/training/nextflow-run/work/0a/0df4a1028c2001758b1841cff92fc7/.command.sh
    }
    ...
    ```

`.command.run` to właściwy skrypt, który Nextflow przekazuje do wykonania.
Opakowuje `.command.sh` we wszystko, co jest potrzebne do jego uruchomienia: konfigurację środowiska, staging wejść i wyjść oraz raportowanie wyniku z powrotem do Nextflow.
Przy executorze `local` Nextflow po prostu uruchamia ten skrypt na tej samej maszynie.

Dokładnie to zmienia się, gdy ustawisz inny `executor`.
Dla schedulera HPC, takiego jak Slurm lub PBS, Nextflow generuje ten sam rodzaj skryptu opakowującego, dodaje nagłówek specyficzny dla schedulera, który widziałeś w [sekcji 2.2](#22-backend-specific-syntax-is-abstracted-away) (przetłumaczony z ustawień `cpus`, `memory` i `queue`), a następnie przekazuje wynik do własnego polecenia zgłoszeniowego schedulera — na przykład `sbatch` dla Slurma.
Od tego momentu Nextflow odpytuje scheduler o status zadania, zamiast bezpośrednio obserwować lokalny proces.
Backendowe platformy chmurowe działają nieco inaczej, bo są sterowane wywołaniami API zamiast poleceniem zgłoszeniowym, ale ta sama podstawowa idea obowiązuje: ten sam skrypt zadania jest uruchamiany, zmienia się tylko sposób jego uruchomienia i śledzenia.

### Podsumowanie

Wiesz już, jak zmienić executor, żeby korzystać z innej infrastruktury obliczeniowej, że Nextflow abstrahuje składnię zgłoszeń specyficzną dla backendu, oraz co faktycznie dzieje się za kulisami, gdy zadanie jest uruchamiane na innym backendzie.

### Co dalej?

Przejdź do [Części 2](./02_resources_and_retries.md), gdzie dowiesz się, jak profilować i alokować zasoby obliczeniowe oraz obsługiwać błędy zadań przy użyciu mechanizmu ponawiania.

---

## Podsumowanie

W tej części nauczyłeś się:

- Przełączać technologię pakowania oprogramowania między Docker'em a Condą
- Dodawać dyrektywę `conda` do definicji procesu
- Zmieniać platformę wykonawczą za pomocą dyrektywy `executor`
- Sprawdzać, co Nextflow faktycznie generuje i uruchamia dla zadania, oraz jak to się zmienia w zależności od executora
