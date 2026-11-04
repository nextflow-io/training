# Część 2: Uruchamianie pipeline'ów z wiersza poleceń

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Tłumaczenie wspomagane przez AI - [dowiedz się więcej i zasugeruj ulepszenia](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

W [Części 1](./01_run_with_seqera.md) uruchomiłeś/uruchomiłaś nf-core/rnaseq z interfejsu webowego Seqera.
Teraz zrobimy to samo z wiersza poleceń przy użyciu `tw` CLI i dodamy nowy pipeline do Twojego workspace'u.

---

## 1. Uruchamianie pipeline'ów z wiersza poleceń

W widoku uruchomienia kliknij zakładkę **Command line**.
Zobaczysz dokładne polecenie `nextflow run`, które Platform skonstruował i wysłał w Twoim imieniu — tego samego rodzaju polecenie, które uruchamiałeś/uruchamiałaś ręcznie w kursie Use nf-core.

Platform nie zastępuje Nextflow'a — orkiestruje go.
Wszystko, co możesz zrobić przez interfejs webowy, możesz też wykonać z terminala przy użyciu `tw` CLI, czyli narzędzia wiersza poleceń do interakcji z API Platform.
Jest to przydatne przy automatyzowaniu uruchomień ze skryptów lub pipeline'ów CI/CD.

Zrobimy to teraz z tego samego codespace'u, którego używałeś/używałaś we wcześniejszych kursach.

### 1.1. Instalacja tw CLI

Uruchom następujące polecenia w terminalu Codespace, aby pobrać i zainstalować plik binarny `tw`:

```bash
curl -fsSL https://github.com/seqeralabs/tower-cli/releases/latest/download/tw-linux-x86_64 -o tw
chmod +x tw
sudo mv tw /usr/local/bin/
```

Zweryfikuj instalację:

```bash
tw --version
```

??? success "Wyjście polecenia"

    ```console
    Tower CLI version 0.40.0 (build 26db579)
    ```

`tw` CLI jest zainstalowane i gotowe do konfiguracji.

### 1.2. Uzyskanie tokenu dostępu

`tw` CLI uwierzytelnia się w Seqera przy użyciu osobistego tokenu dostępu.

1. W interfejsie webowym Seqera kliknij swój awatar w prawym górnym rogu i wybierz **Your tokens**.
2. Kliknij **Add token**, nadaj mu nazwę (np. `training`) i kliknij **Add**.
3. Skopiuj wartość tokenu — zostanie wyświetlona tylko raz.
   Jeśli nie zapiszesz go od razu, będziesz musiał/musiała wygenerować nowy.

### 1.3. Konfiguracja CLI

Dla wygody skonfigurujemy plik konfiguracyjny zawierający token dostępu, który właśnie wygenerowałeś/wygenerowałaś, oraz identyfikator workspace'u.

Otwórz plik `.seqera_config` w tym katalogu w edytorze i ustaw dwie zmienne:

- **`TOWER_ACCESS_TOKEN`**: token wygenerowany w sekcji 1.2
- **`TOWER_WORKSPACE_ID`**: numeryczny identyfikator Twojego workspace'u (kolumna `ID` w wynikach `tw workspaces list`, które uruchomisz w sekcji 1.4)

Po uzupełnieniu wartości załaduj konfigurację:

```bash
source .seqera_config
```

Zweryfikuj połączenie:

```bash
tw info
```

??? success "Wyjście polecenia"

    ```console
        Details
    -------------------------+-----------------------------
     Tower API endpoint      | https://api.cloud.seqera.io
     Tower API version       | 1.150.0
     Tower version           | 26.1.0-cycle54
     CLI version             | 0.30.0 (fde9dec)
     CLI minimum API version | 1.148.0
     Authenticated user      | <your-name>

    System health status
    ---------------------------------------+----
     Remote API server connection check    | OK
     Tower API version check               | OK
     Authentication API credential's token | OK
    ```

`tw` CLI jest teraz uwierzytelnione i połączone z Twoim kontem Seqera.
Uruchamiaj `source .seqera_config` na początku każdej sesji Codespace, aby ponownie załadować konfigurację.

!!! tip "Wskazówka"

    Jeśli Twój workspace nie ma ustawionego podstawowego środowiska obliczeniowego, możesz dodać `export TOWER_COMPUTE_ENV=<compute-env-name>` do pliku konfiguracyjnego, aby ustawić domyślne.
    Każdą wartość konfiguracyjną można nadpisać w wierszu poleceń, przekazując flagę jawnie (np. `--compute-env other-env`).
    Pełna lista opcji i zmiennych środowiskowych znajduje się w [dokumentacji tw CLI](https://docs.seqera.io/platform/latest/cli/reference).

### 1.4. Eksploracja workspace'u z CLI

Wyświetl listę workspace'ów, do których masz dostęp:

```bash
tw workspaces list
```

??? success "Wyjście polecenia"

    ```console
    Available workspaces:
    ID              | Name            | Full Name                | Visibility
    --------------- | --------------- | ------------------------ | ----------
    <workspace-id>  | my-workspace    | my-org/my-workspace      | PRIVATE
    ```

Wyświetl uruchomienia w swoim workspace'ie, w tym uruchomienie nf-core/rnaseq, które właśnie zainicjowałeś/zainicjowałaś:

```bash
tw runs list
```

??? success "Wyjście polecenia"

    ```console
    Pipeline runs for my-org/my-workspace workspace:
    ID        | Status   | Name          | Pipeline               | Run name
    --------- | -------- | ------------- | ---------------------- | --------
    <run-id>  | RUNNING  | ...           | nf-core/rnaseq         | happy_curie
    ```

To samo uruchomienie, które monitorujesz w interfejsie webowym, jest widoczne tutaj.

!!! note "Uwaga"

    Ponieważ `TOWER_WORKSPACE_ID` jest ustawione w `.seqera_config`, możesz pominąć `--workspace` we wszystkich poleceniach `tw`.
    Bez konfiguracji należałoby przekazać je jawnie:

    ```bash
    tw runs list --workspace <org>/<workspace>
    ```

Wszystko widoczne w interfejsie webowym jest dostępne z CLI.

### 1.5. Uruchamianie nf-core/rnaseq z CLI

Pipeline dodany do Twojego workspace'u w [Części 1](./01_run_with_seqera.md) jest dostępny po nazwie w CLI.
Uruchom go z profilem `test`:

```bash
tw launch nf-core-rnaseq --profile test
```

??? success "Wyjście polecenia"

    ```console
    Launching pipeline nf-core-rnaseq
    Run name: focused_einstein
    https://cloud.seqera.io/orgs/my-org/workspaces/my-workspace/watch/<run-id>
    ```

Otwórz link w przeglądarce i potwierdź, że uruchomienie pojawia się w panelu **Runs**.

Gdy zobaczysz, że działa, potwierdzisz, że CLI i interfejs webowy to dwa widoki tego samego workspace'u.

!!! note "Uwaga"

    Możesz też przekazać pełny adres URL GitHub bezpośrednio do `tw launch` bez wcześniejszego dodawania pipeline'u do workspace'u.
    Jednak jawne dodanie pipeline'u przed uruchomieniem jest generalnie lepszym rozwiązaniem: zapisuje konfigurację na potrzeby przyszłych uruchomień, udostępnia go po nazwie i sprawia, że jest widoczny dla wszystkich członków workspace'u w Launchpadzie.

    Pipeline można dodać do workspace'u bezpośrednio z wiersza poleceń przy użyciu `tw`.
    Następna sekcja pokazuje, jak to zrobić z pipeline'em nf-core/demo.

### Podsumowanie

Wiesz już, jak uwierzytelnić `tw` CLI, przeglądać swój workspace i uruchamiać zapisany pipeline z terminala.

### Co dalej?

Dodaj nowy pipeline do swojego workspace'u z wiersza poleceń i uruchom go.

---

## 2. Dodawanie nowego pipeline'u i jego uruchamianie

Każdy pipeline Nextflow'a dostępny na GitHub można dodać do workspace'u za pomocą `tw pipelines add`, o ile posiada punkt wejścia `main.nf` i plik `nextflow.config` w swoim katalogu głównym.
nf-core/demo to dobry przykład do ćwiczeń: uruchamiałeś/uruchamiałaś go już w kursie Use nf-core, więc wiesz, co robi i czego się spodziewać.

### 2.1. Dodawanie nf-core/demo do workspace'u

Uruchom następujące polecenie, aby zarejestrować pipeline w swoim workspace'ie:

```bash
tw pipelines add \
    --name nf-core-demo \
    https://github.com/nf-core/demo
```

??? success "Wyjście polecenia"

    ```console
    New pipeline 'nf-core-demo' added at my-org/my-workspace workspace
    ```

Pipeline jest teraz zarejestrowany i pojawi się w Launchpadzie.

### 2.2. Weryfikacja obecności w Launchpadzie

Wyświetl listę pipeline'ów w swoim workspace'ie, aby potwierdzić dodanie:

```bash
tw pipelines list
```

??? success "Wyjście polecenia"

    ```console
    Pipelines at my-org/my-workspace:
    ID   | Name            | Repository
    ---- | --------------- | -----------------------------------------
    ...  | nf-core-demo    | https://github.com/nf-core/demo
    ...  | nf-core-rnaseq  | ...
    ```

Otwórz swój workspace w przeglądarce i kliknij **Launchpad**, aby potwierdzić, że nf-core/demo pojawia się obok nf-core/rnaseq.

!!! tip "Wskazówka"

    Pipeline'y możesz też dodawać przez interfejs webowy: na lewym pasku bocznym kliknij **Launchpad**, następnie **Add pipeline** i wypełnij formularz.

Kliknij przycisk **Launch** przy pozycji nf-core/demo, aby otworzyć formularz uruchomienia.
Zobaczysz, że parametry `input` i `outdir` są podświetlone na czerwono — są to pola wymagane bez wartości domyślnych, ponieważ `tw pipelines add` rejestruje jedynie źródło pipeline'u bez wstępnej konfiguracji parametrów.
Kolejne dwie sekcje pokazują, jak podać te wartości: najpierw przez formularz webowy, a następnie z wiersza poleceń.

### 2.3. Uruchamianie nf-core/demo z interfejsu webowego

Mając otwarty formularz uruchomienia, wypełnij dwa wymagane parametry.

Dla `input` podaj adres URL testowego samplesheet z profilu testowego nf-core/demo.
Znajdziesz go w pliku `conf/test.config` w repozytorium pipeline'u, który przeglądałeś/przeglądałaś w kursie Use nf-core:

```
https://raw.githubusercontent.com/nf-core/test-datasets/viralrecon/samplesheet/samplesheet_test_illumina_amplicon.csv
```

Dla `outdir` podaj ścieżkę do magazynu w chmurze, gdzie pipeline może zapisać wyniki.
Użyj bucketu skonfigurowanego dla Twojego workspace'u z podkatalogiem, aby zachować porządek wśród uruchomień:

```
s3://my-bucket/demo-results
```

Po wypełnieniu obu pól kliknij niebieski przycisk **Launch**.

Uruchomienie pojawi się w panelu **Runs** i powinno zakończyć się w ciągu kilku minut na zbiorze testowym.
Kliknij w uruchomienie, aby przejrzeć tabelę zadań i ewentualne raporty wykonania.

### 2.4. Uruchamianie nf-core/demo z CLI

W odróżnieniu od `nextflow run`, polecenie `tw launch` nie przyjmuje indywidualnych flag parametrów, takich jak `--input` czy `--outdir`.
Parametry należy podać w pliku w formacie YAML lub JSON, przekazanym za pomocą `--params-file`.
Takie podejście sprzyja odtwarzalności: zapisany plik parametrów dokumentuje dokładnie, jakie wartości zostały użyte w danym uruchomieniu, co ułatwia jego powtórzenie lub udostępnienie konfiguracji.

Utwórz plik parametrów w swoim katalogu roboczym:

```bash
touch params.yaml
```

Otwórz go w edytorze i dodaj ścieżkę wyjściową:

```yaml title="params.yaml"
outdir: "s3://my-bucket/demo-results-cli"
```

Teraz możesz uruchomić pipeline, używając profilu `test` (który dostarcza samplesheet dla `input`) i pliku parametrów (który dostarcza `outdir`):

```bash
tw launch nf-core-demo --profile test --params-file params.yaml
```

??? success "Wyjście polecenia"

    ```console
    Launching pipeline nf-core-demo
    Run name: focused_feynman
    https://cloud.seqera.io/orgs/my-org/workspaces/my-workspace/watch/<run-id>
    ```

Otwórz link, aby potwierdzić, że uruchomienie pojawia się w panelu **Runs**.

!!! tip "Wskazówka"

    Plik parametrów możesz dołączyć już na etapie początkowej konfiguracji, jeśli chcesz ustawić pewne wartości domyślne oraz dodatkowe właściwości odpowiadające temu, co zrobiliśmy wcześniej przez formularz webowy:

    ```bash
    tw pipelines add \
      --name nf-core-demo-gg \
      --description "Demo pipeline with defaults for testing" \
      --params-file params.yaml \
      --profile test \
      https://github.com/nf-core/demo
    ```

### Podsumowanie

Wiesz już, jak dodać dowolny pipeline Nextflow'a hostowany na GitHub do swojego workspace'u i uruchomić go — zarówno z interfejsu webowego przez ręczne wypełnienie parametrów, jak i z `tw` CLI przez połączenie profilu z plikiem parametrów.

---

## Podsumowanie

W tej części nauczyłeś/nauczyłaś się:

- Uwierzytelniać `tw` CLI i uruchamiać zapisany pipeline z terminala
- Dodawać nowy pipeline z GitHub przy użyciu CLI i weryfikować jego obecność w Launchpadzie
- Uruchamiać pipeline z interfejsu webowego Seqera przez ręczne wypełnienie wymaganych parametrów
- Uruchamiać pipeline z CLI przy użyciu profilu Nextflow'a i pliku parametrów
