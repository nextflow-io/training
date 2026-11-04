# Część 3: Zarządzanie uruchomieniami workflow'u

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Tłumaczenie wspomagane przez AI - [dowiedz się więcej i zasugeruj ulepszenia](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

Uruchamiając pipeline'y wielokrotnie, gromadzisz historię uruchomień i stare katalogi `work/`.
W [Części 1](./01_run_nextflow.md#23-use-resume-to-skip-completed-work) korzystałeś już z flagi `-resume`, aby pominąć pracę wykonaną wcześniej.
Tutaj nauczysz się generować raporty z uruchomień, przeglądać historię poprzednich uruchomień za pomocą [`nextflow log`](https://nextflow.io/docs/latest/reference/cli.html#log) oraz usuwać stare katalogi robocze, których już nie potrzebujesz, przy użyciu [`nextflow clean`](https://nextflow.io/docs/latest/reference/cli.html#clean).

---

## 1. Generowanie raportów pipeline'u

Nextflow może generować kilka rodzajów raportów z uruchomień — każdy dodaje się osobną flagą `-with-*`: raport wykonania (`-with-report`), oś czasu wykonania (`-with-timeline`), plik śledzenia zadań (`-with-trace`) oraz diagram workflow'u (`-with-dag`).
Tutaj wygenerujemy pierwsze dwa; pozostałe znajdziesz w sekcji [Execution reports](https://nextflow.io/docs/latest/reports.html) w dokumentacji Nextflow'a.

### 1.1. Generowanie raportu wykonania

Dodaj `-with-report` do dowolnego polecenia `nextflow run`, aby po zakończeniu pipeline'u wygenerować raport HTML:

```bash
nextflow run main.nf --input data/greetings.csv --character turkey -with-report
```

??? success "Wyjście polecenia"

    ```console hl_lines="6 7 8 9"
    N E X T F L O W   ~  version 26.04.4

    Launching `main.nf` [intergalactic_dalembert] revision: ce74f81996

    executor >  local (8)
    [34/23f10f] sayHello (2)       | 3 of 3 ✔
    [af/cab69d] convertToUpper (3) | 3 of 3 ✔
    [9e/d73afb] collectGreetings   | 1 of 1 ✔
    [3c/392db0] cowpy              | 1 of 1 ✔

    Outputs:

      /workspaces/training/nextflow-run/results

      first_output:
        - batch/intermediates/Bonjour-output.txt
        - batch/intermediates/Hello-output.txt
        - batch/intermediates/Hola-output.txt

      uppercased:
        - batch/intermediates/UPPER-Bonjour-output.txt
        - batch/intermediates/UPPER-Hello-output.txt
        - batch/intermediates/UPPER-Hola-output.txt

      collected: batch/intermediates/COLLECTED-batch-output.txt

      batch_report: batch/batch-report.txt

      cowpy_art: batch/cowpy-COLLECTED-batch-output.txt
    ```

Nextflow zapisuje raport do pliku o nazwie `report-<timestamp>.html` w bieżącym katalogu roboczym.
Otwórz go w przeglądarce, aby zobaczyć podsumowanie wykonania, tabelę wszystkich zadań z ich statusem i czasem działania oraz wykresy użycia zasobów z podziałem na procesy.

Zakładka **Tasks** wyświetla każde zadanie uruchomione przez pipeline wraz z nazwą procesu, statusem i zużyciem zasobów:

![Tabela zadań w raporcie wykonania](img/execution_report_tasks.png)

Raport jest szczególnie przydatny, gdy pipeline działa dłużej niż oczekiwano lub gdy jakieś zadanie kończy się błędem — tabela zadań pokazuje dokładnie, gdzie czas był spędzony i które zadania zakończyły się sukcesem, a które niepowodzeniem.

### 1.2. Generowanie osi czasu wykonania

Dodaj `-with-timeline` do uruchomienia, aby uzyskać widok w stylu wykresu Gantta pokazujący, kiedy każde zadanie było wykonywane:

```bash
nextflow run main.nf --input data/greetings.csv --character turkey -with-timeline
```

??? success "Wyjście polecenia"

    ```console hl_lines="6 7 8 9"
    N E X T F L O W   ~  version 26.04.4

    Launching `main.nf` [jolly_noyce] revision: ce74f81996

    executor >  local (8)
    [ad/e92ef3] sayHello (3)       | 3 of 3 ✔
    [2a/df8a8d] convertToUpper (2) | 3 of 3 ✔
    [be/7fb72a] collectGreetings   | 1 of 1 ✔
    [63/dc9bd6] cowpy              | 1 of 1 ✔

    Outputs:

      /workspaces/training/nextflow-run/results

      first_output:
        - batch/intermediates/Bonjour-output.txt
        - batch/intermediates/Hello-output.txt
        - batch/intermediates/Hola-output.txt

      uppercased:
        - batch/intermediates/UPPER-Bonjour-output.txt
        - batch/intermediates/UPPER-Hello-output.txt
        - batch/intermediates/UPPER-Hola-output.txt

      collected: batch/intermediates/COLLECTED-batch-output.txt

      batch_report: batch/batch-report.txt

      cowpy_art: batch/cowpy-COLLECTED-batch-output.txt
    ```

Nextflow zapisuje oś czasu do pliku o nazwie `timeline-<timestamp>.html`.
Otwórz go w przeglądarce, aby zobaczyć pasek dla każdego zadania, rozmieszczony i wyskalowany według czasu rozpoczęcia i czasu trwania:

![Oś czasu wykonania](img/execution_timeline.png)

Oś czasu uwidacznia charakterystyczny kształt „rozejście–zbieżność" z [Części 1](./01_run_nextflow.md#31-run-the-workflow): trzy zadania `sayHello` działają równolegle, następnie trzy zadania `convertToUpper`, po czym `collectGreetings` i `cowpy` wykonują się jedno po drugim, gdyż każde z nich zależy od wszystkich poprzednich.

### Podsumowanie

Wiesz już, jak generować raport wykonania HTML za pomocą `-with-report` oraz oś czasu wykonania za pomocą `-with-timeline`, a także gdzie szukać pozostałych typów raportów obsługiwanych przez Nextflow'a.

### Co dalej?

Dowiedz się, jak przeglądać historię poprzednich uruchomień.

---

## 2. Przeglądanie historii poprzednich uruchomień

Niezależnie od tego, czy rozwijasz pipeline, czy uruchamiasz go produkcyjnie, w pewnym momencie będziesz potrzebować informacji o poprzednich uruchomieniach.

### 2.1. Plik historii

Za każdym razem, gdy uruchamiasz workflow Nextflow'a, do pliku dziennika o nazwie `history` — znajdującego się w ukrytym katalogu `.nextflow` w bieżącym katalogu roboczym — dopisywana jest nowa linia.

??? abstract "Zawartość pliku"

    ```txt title=".nextflow/history" linenums="1"
    2026-09-13 01:47:58	2.8s	determined_ramanujan	OK	ce74f81996b8ce999853fdbb25ba7969	03078950-5011-414f-992b-8fe1bd219ac2	nextflow run main.nf --input data/greetings.csv --character turkey
    2026-09-13 01:48:03	3s	prickly_cuvier	OK	ce74f81996b8ce999853fdbb25ba7969	5a9e8bf7-a7db-4c2c-b1d8-3bd4c7266528	nextflow run main.nf --input data/greetings.csv --character tux
    2026-09-13 01:48:09	2.3s	elegant_panini	OK	ce74f81996b8ce999853fdbb25ba7969	517ccc29-9db7-4d8d-8ea0-c8db6e9ac653	nextflow run main.nf --input data/greetings.csv --character stegosaurus
    2026-09-13 01:48:14	1.5s	deadly_lamport	OK	ce74f81996b8ce999853fdbb25ba7969	517ccc29-9db7-4d8d-8ea0-c8db6e9ac653	nextflow run main.nf --input data/greetings.csv --character stegosaurus -resume
    ```

Każda linia zawiera znacznik czasu, czas trwania, nazwę uruchomienia, status, identyfikator rewizji, identyfikator sesji oraz pełne polecenie użyte do uruchomienia z tego katalogu.

Przyjrzyj się dwóm ostatnim liniom: to dwa osobne wywołania (jedno zwykłe, drugie z `-resume`) dokładnie tego samego polecenia, które mają ten sam identyfikator sesji.
Identyfikator sesji zmienia się tylko wtedy, gdy uruchamiasz naprawdę nowe uruchomienie; użycie `-resume` go zachowuje — właśnie po nim Nextflow wie, której pamięci podręcznej użyć.

### 2.2. Użycie `nextflow log` dla wygodniejszego podglądu

Czytanie surowego pliku historii działa, ale `nextflow log` formatuje te same informacje z nagłówkiem:

```bash
nextflow log
```

??? success "Wyjście polecenia"

    ```console linenums="1"
    TIMESTAMP          	DURATION	RUN NAME            	STATUS	REVISION ID	SESSION ID                          	COMMAND
    2026-09-13 01:47:58	2.8s    	determined_ramanujan	OK    	ce74f81996 	03078950-5011-414f-992b-8fe1bd219ac2	nextflow run main.nf --input data/greetings.csv --character turkey
    2026-09-13 01:48:03	3s      	prickly_cuvier      	OK    	ce74f81996 	5a9e8bf7-a7db-4c2c-b1d8-3bd4c7266528	nextflow run main.nf --input data/greetings.csv --character tux
    2026-09-13 01:48:09	2.3s    	elegant_panini      	OK    	ce74f81996 	517ccc29-9db7-4d8d-8ea0-c8db6e9ac653	nextflow run main.nf --input data/greetings.csv --character stegosaurus
    2026-09-13 01:48:14	1.5s    	deadly_lamport      	OK    	ce74f81996 	517ccc29-9db7-4d8d-8ea0-c8db6e9ac653	nextflow run main.nf --input data/greetings.csv --character stegosaurus -resume
    ```

Nextflow przechowuje informacje o pamięci podręcznej używanej przez `-resume` w katalogu `.nextflow/cache`, z kluczem opartym na identyfikatorze sesji.
Dlatego odszukanie właściwej nazwy uruchomienia lub identyfikatora sesji jest pierwszym krokiem, gdy chcesz zbadać lub wyczyścić poprzednie uruchomienie.

### Podsumowanie

Wiesz już, gdzie Nextflow zapisuje historię poprzednich uruchomień i jak ją przeglądać za pomocą `nextflow log`.

### Co dalej?

Dowiedz się, jak usuwać stare katalogi robocze, których już nie potrzebujesz.

---

## 3. Usuwanie starych katalogów roboczych

Każde uruchomienie pozostawia katalogi zadań w `work/`, nawet po skopiowaniu interesujących Cię wyników do `results/`.
Przy intensywnym rozwijaniu pipeline'ów te podkatalogi szybko się kumulują, dlatego Nextflow udostępnia polecenie `nextflow clean` do usuwania tych, których już nie potrzebujesz.

### 3.1. Określenie kryteriów usuwania

`nextflow clean` obsługuje kilka sposobów wyboru tego, co usunąć; pełną listę znajdziesz w [dokumentacji referencyjnej](https://www.nextflow.io/docs/latest/reference/cli.html#clean).
Tutaj usuniemy wszystko z uruchomień poprzedzających wybrane uruchomienie, korzystając z jego nazwy.

Odszukaj nazwę najnowszego uruchomienia, które chcesz zachować, używając `nextflow log`; w [przykładzie z sekcji 2.2](#22-use-nextflow-log-for-a-friendlier-view) jest to `elegant_panini` — ostatnie zwykłe uruchomienie przed tym z `-resume`.
Nazwa uruchomienia to generowany automatycznie dwuczłonowy ciąg znaków widoczny w linii `Launching (...)` w konsoli lub w kolumnie `RUN NAME` w wynikach `nextflow log`.

### 3.2. Próbne uruchomienie

Najpierw dodaj `-n`, aby sprawdzić, co dane polecenie by usunęło, bez faktycznego usuwania czegokolwiek:

```bash
nextflow clean -before elegant_panini -n
```

??? success "Wyjście polecenia"

    ```console
    Would remove /workspaces/training/nextflow-run/work/e5/9ec3fe24deb5d1ffa478dff5e5e1d4
    Would remove /workspaces/training/nextflow-run/work/75/36b1ece86b213c331852d0a3dc659b
    Would remove /workspaces/training/nextflow-run/work/be/8157808f28699259b5ff37e4bf9f72
    Would remove /workspaces/training/nextflow-run/work/f3/d828b7ba6b315eb08e0f96168119f5
    Would remove /workspaces/training/nextflow-run/work/94/7a2d7ef67d81b6e66430f300816698
    Would remove /workspaces/training/nextflow-run/work/eb/e37c5c8638c45916de3d9794d6d22f
    Would remove /workspaces/training/nextflow-run/work/3f/f56a705a6370f6c615aff22d3f240a
    Would remove /workspaces/training/nextflow-run/work/25/4922fa35ba2c6087777ddce20e9f2f
    Would remove /workspaces/training/nextflow-run/work/a4/2de7bcc2f8fc974261f7280c1480f7
    Would remove /workspaces/training/nextflow-run/work/48/8c49bded446d4fcb78c8125318c5e3
    Would remove /workspaces/training/nextflow-run/work/8e/c9834a9bdd7c21de2991d959201a9d
    Would remove /workspaces/training/nextflow-run/work/d5/a9055dc8c817c6c1f436494685955c
    Would remove /workspaces/training/nextflow-run/work/27/9473df9e55e35c8869d9571b480926
    Would remove /workspaces/training/nextflow-run/work/d3/ebc066e1d844232b9f5ee0a7d87464
    Would remove /workspaces/training/nextflow-run/work/91/17ef1815ed06b8f177d0135d8ef0dc
    Would remove /workspaces/training/nextflow-run/work/12/0443f9fcff27a552d4bc73aaedf24a
    ```

To 16 katalogów zadań: 8 z uruchomienia `turkey` i 8 z uruchomienia `tux` — dokładnie tyle, ile można się spodziewać po dwóch pełnych uruchomieniach tego czteroprocesowego pipeline'u.
Samo uruchomienie `elegant_panini` oraz zadania z pamięci podręcznej ponownie użyte przez uruchomienie z `-resume` pozostają nienaruszone.

Twój wynik będzie zawierał inne nazwy katalogów, a liczba linii zależy od tego, ile uruchomień wykonałeś. Jeśli nie widzisz żadnych linii, albo nazwa uruchomienia nie pasuje do żadnej w Twoim dzienniku, albo nie ma nic do usunięcia przed nią.

### 3.3. Właściwe usuwanie

Gdy próbne uruchomienie wygląda poprawnie, uruchom to samo polecenie z `-f` zamiast `-n`:

```bash
nextflow clean -before elegant_panini -f
```

??? success "Wyjście polecenia"

    ```console
    Removed /workspaces/training/nextflow-run/work/e5/9ec3fe24deb5d1ffa478dff5e5e1d4
    Removed /workspaces/training/nextflow-run/work/75/36b1ece86b213c331852d0a3dc659b
    Removed /workspaces/training/nextflow-run/work/be/8157808f28699259b5ff37e4bf9f72
    Removed /workspaces/training/nextflow-run/work/f3/d828b7ba6b315eb08e0f96168119f5
    Removed /workspaces/training/nextflow-run/work/94/7a2d7ef67d81b6e66430f300816698
    Removed /workspaces/training/nextflow-run/work/eb/e37c5c8638c45916de3d9794d6d22f
    Removed /workspaces/training/nextflow-run/work/3f/f56a705a6370f6c615aff22d3f240a
    Removed /workspaces/training/nextflow-run/work/25/4922fa35ba2c6087777ddce20e9f2f
    Removed /workspaces/training/nextflow-run/work/a4/2de7bcc2f8fc974261f7280c1480f7
    Removed /workspaces/training/nextflow-run/work/48/8c49bded446d4fcb78c8125318c5e3
    Removed /workspaces/training/nextflow-run/work/8e/c9834a9bdd7c21de2991d959201a9d
    Removed /workspaces/training/nextflow-run/work/d5/a9055dc8c817c6c1f436494685955c
    Removed /workspaces/training/nextflow-run/work/27/9473df9e55e35c8869d9571b480926
    Removed /workspaces/training/nextflow-run/work/d3/ebc066e1d844232b9f5ee0a7d87464
    Removed /workspaces/training/nextflow-run/work/91/17ef1815ed06b8f177d0135d8ef0dc
    Removed /workspaces/training/nextflow-run/work/12/0443f9fcff27a552d4bc73aaedf24a
    ```

`nextflow clean` opróżnia katalogi zadań, ale pozostawia dwuznakowe katalogi nadrzędne (takie jak `e5/`) na miejscu.

!!! warning "Ostrzeżenie"

    Usunięcie katalogów roboczych z poprzednich uruchomień usuwa je z pamięci podręcznej Nextflow'a i kasuje wszelkie wyniki przechowywane wyłącznie tam.
    Uniemożliwia to wznowienie wykonania bez ponownego uruchamiania odpowiednich procesów, dlatego usuwaj tylko te uruchomienia, z których na pewno nie będziesz potrzebować wznawiać.
    To również powód, dla którego warto publikować wszystko, na czym Ci zależy, do `results/` z opcją `mode 'copy'`, zamiast polegać na katalogu `work/` lub trybie publikowania `symlink`.

### Podsumowanie

Wiesz już, jak usuwać stare katalogi robocze za pomocą `nextflow clean`, i rozumiesz, że wiąże się to z utratą możliwości wznawiania tych uruchomień.

### Co dalej?

Dowiedz się, jak uruchamiać pipeline'y bezpośrednio ze zdalnych repozytoriów, takich jak GitHub, w [Części 4](./04_remote_repositories.md).

---

## Podsumowanie

W tej części nauczyłeś się:

- Generować raport wykonania HTML za pomocą `-with-report` oraz oś czasu wykonania za pomocą `-with-timeline`
- Przeglądać historię poprzednich uruchomień za pomocą `nextflow log`
- Usuwać stare katalogi robocze za pomocą `nextflow clean` i rozumieć związany z tym kompromis dotyczący wznawiania
