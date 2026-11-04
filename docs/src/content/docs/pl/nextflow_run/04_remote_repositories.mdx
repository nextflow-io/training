# Część 4: Uruchamianie zdalnych pipeline'ów

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Tłumaczenie wspomagane przez AI - [dowiedz się więcej i zasugeruj ulepszenia](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

Do tej pory uruchamiałeś skrypty workflow'ów przechowywane lokalnie.
W praktyce często będziesz chciał uruchamiać pipeline'y opublikowane w zdalnych repozytoriach, takich jak GitHub, bez konieczności ich samodzielnego pobierania.

Nextflow sprawia, że jest to proste: możesz uruchomić dowolny pipeline bezpośrednio z adresu URL repozytorium Git.

---

## 1. Uruchamianie pipeline'u z GitHub

Podstawowa składnia uruchamiania zdalnego pipeline'u to `nextflow run <repository>`, gdzie `<repository>` może być ścieżką do repozytorium GitHub, np. `nextflow-io/hello`, pełnym adresem URL lub ścieżką do GitLab, Bitbucket albo innej usługi hostingowej Git.

### 1.1. Uruchomienie pipeline'u

Uruchom oficjalny demonstracyjny pipeline Nextflow'a „hello".
Jest to inny, znacznie prostszy pipeline niż ten, z którym pracowałeś w tym kursie: powstał wcześniej niż pipeline „Hello" używany w tym szkoleniu i po prostu wypisuje powitanie w kilku zakodowanych na stałe językach — nie spodziewaj się więc wejścia CSV ani ASCII art, do których jesteś przyzwyczajony.

```bash
nextflow run nextflow-io/hello
```

??? success "Wyjście polecenia"

    ```console
    N E X T F L O W   ~  version 26.04.4

    Pulling nextflow-io/hello ...
     downloaded from https://github.com/nextflow-io/hello.git
    Launching `https://github.com/nextflow-io/hello` [sleepy_swanson] revision: 3c2cdc9823 [master]

    executor >  local (4)
    [ba/08236d] sayHello (4) | 4 of 4 ✔
    Ciao world!

    Hello world!

    Bonjour world!

    Hola world!
    ```

### 1.2. Znajdowanie miejsca, w którym pipeline jest zapisany w cache

Przy pierwszym uruchomieniu zdalnego pipeline'u Nextflow pobiera go i zapisuje lokalnie w cache.
Kolejne uruchomienia korzystają z wersji zapisanej w cache, chyba że jawnie zażądasz aktualizacji.

Domyślnie Nextflow zapisuje pobrane pipeline'y w katalogu `$NXF_HOME/assets`.
Aby dowiedzieć się, gdzie trafił konkretny pipeline i jakie rewizje są dostępne, zapytaj Nextflow'a bezpośrednio:

```bash
nextflow info nextflow-io/hello
```

??? success "Wyjście polecenia"

    ```console
     project name: nextflow-io/hello
     repository  : https://github.com/nextflow-io/hello
     local path  : /workspaces/.nextflow/assets/.repos/nextflow-io/hello
     main script : main.nf
     revisions   :
     > master (default)
       mybranch
       testing
       v1.1 [t]
       v1.2 [t]
       v1.3 [t]
    ```

    Nextflow oznacza każdą rewizję, którą już pobrałeś lokalnie, symbolem `>`; pozostałe są dostępne, ale nie zostały jeszcze pobrane do kopii roboczej.

Możesz też wylistować wszystkie dotychczas pobrane pipeline'y za pomocą `nextflow list`:

```bash
nextflow list
```

??? success "Wyjście polecenia"

    ```console
    nextflow-io/hello
    ```

Kurs [Use nf-core](../nfcore_use/01_run_demo.md#12-retrieve-the-pipeline-code) omawia ten mechanizm cache'owania bardziej szczegółowo, w tym sposób przeglądania kodu źródłowego pobranego pipeline'u.

### Podsumowanie

Wiesz już, jak uruchomić pipeline bezpośrednio z repozytorium GitHub bez jego samodzielnego pobierania, oraz gdzie znaleźć go lokalnie po uruchomieniu.

### Co dalej?

Dowiedz się, jak przypiąć konkretną wersję zdalnego pipeline'u w celu zapewnienia odtwarzalności.

---

## 2. Określanie wersji dla odtwarzalności

Domyślnie Nextflow uruchamia najnowszą rewizję z domyślnej gałęzi.
Możesz przypiąć konkretną wersję (tag), gałąź lub commit za pomocą flagi `-r`.

### 2.1. Przypinanie konkretnej rewizji

```bash
nextflow run nextflow-io/hello -r v1.3
```

??? success "Wyjście polecenia"

    ```console
    N E X T F L O W   ~  version 26.04.4

    Pulling nextflow-io/hello:v1.3 ...
     downloaded from https://github.com/nextflow-io/hello.git
    Launching `https://github.com/nextflow-io/hello` [sick_carson] revision: 2ce0b0e294 [v1.3]

    executor >  local (4)
    [61/e11f77] sayHello (4) | 4 of 4 ✔
    Ciao world!

    Bonjour world!

    Hello world!

    Hola world!
    ```

Nextflow pobiera daną rewizję przy pierwszym żądaniu — stąd linie `Pulling` i `downloaded from`; kolejne żądania tej samej rewizji pomijają ten krok i przechodzą od razu do `Launching`.
Przypinanie dokładnej rewizji jest kluczowe dla odtwarzalności.
Gwarantuje, że Ty i Twoi współpracownicy uruchamiacie dokładnie ten sam kod pipeline'u, niezależnie od tego, co zmieniło się w repozytorium od tamtej pory.

### 2.2. Rewizje obowiązują tylko dla danego wywołania

Przypięcie rewizji za pomocą `-r` wpływa wyłącznie na uruchomienie, w którym ją podajesz — nie zmienia tego, czego użyje późniejsze, zwykłe `nextflow run`.
Spróbuj uruchomić pipeline ponownie bez `-r`:

```bash
nextflow run nextflow-io/hello
```

??? success "Wyjście polecenia"

    ```console
    N E X T F L O W   ~  version 26.04.4

    Launching `https://github.com/nextflow-io/hello` [hungry_maxwell] revision: 3c2cdc9823 [master]

    executor >  local (4)
    [ba/08236d] sayHello (4) | 4 of 4 ✔
    Ciao world!

    Hello world!

    Bonjour world!

    Hola world!
    ```

Mimo że poprzednie uruchomienie jawnie przypinało `v1.3`, to uruchomienie wraca do domyślnej gałęzi (`master`).
Nextflow przechowuje osobną lokalną kopię roboczą dla każdej używanej rewizji — to właśnie pokazują znaczniki `>` w wynikach `nextflow info` — ale nigdy nie zapamiętuje, którą uruchomiłeś ostatnio.
Nazwę domyślnej gałęzi pipeline'u możesz sprawdzić poleceniem `nextflow info <pipeline>` — jest to ta oznaczona jako `(default)`.
Odtwarzalność leży całkowicie po Twojej stronie: zawsze podawaj `-r` jawnie, gdy ma to znaczenie, zamiast zakładać, że rewizja przypięta we wcześniejszym uruchomieniu nadal obowiązuje.

### Podsumowanie

Wiesz już, jak przypiąć zdalny pipeline do konkretnej wersji, gałęzi lub commitu w celu odtwarzalnego wykonania, oraz że przypięcie obowiązuje tylko dla danego wywołania, a nie dla kolejnych uruchomień.

### Co dalej?

Poznałeś już podstawy uruchamiania pipeline'ów Nextflow'a i zarządzania nimi.
Zajrzyj do [Podsumowania kursu](next_steps.md), aby dowiedzieć się, co robić dalej.

---

## Podsumowanie

W tej części nauczyłeś się:

- Uruchamiać pipeline bezpośrednio z repozytorium GitHub bez jego pobierania
- Przypinać zdalny pipeline do konkretnej rewizji w celu zapewnienia odtwarzalności oraz rozumieć, że przypięcie obowiązuje tylko dla danego wywołania
