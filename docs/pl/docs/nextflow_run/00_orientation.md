# Rozpoczęcie pracy

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Tłumaczenie wspomagane przez AI - [dowiedz się więcej i zasugeruj ulepszenia](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

## Uruchom środowisko szkoleniowe

Aby skorzystać z gotowego środowiska, które udostępniamy na GitHub Codespaces, kliknij przycisk "Open in GitHub Codespaces" poniżej. Inne opcje znajdziesz w sekcji [Opcje środowiska](../envsetup/index.md).

Zalecamy otwarcie środowiska szkoleniowego w nowej karcie lub oknie przeglądarki (użyj prawego przycisku myszy, ctrl-click lub cmd-click w zależności od sprzętu), aby móc czytać dalej podczas ładowania środowiska.
Będziesz musiał mieć te instrukcje otwarte równolegle, aby przejść przez kurs.

[![Open in GitHub Codespaces](https://github.com/codespaces/badge.svg)](https://codespaces.new/nextflow-io/training?quickstart=1&ref=master)

### Podstawy środowiska

To środowisko szkoleniowe zawiera całe oprogramowanie, kod i dane potrzebne do pracy z kursem, więc nie musisz niczego instalować samodzielnie.

Codespace jest skonfigurowany z interfejsem VSCode, który zawiera eksplorator systemu plików, edytor kodu i powłokę terminala.
Wszystkie instrukcje podane podczas kursu (np. 'otwórz plik', 'edytuj kod' lub 'uruchom to polecenie') odnoszą się do tych trzech części interfejsu VSCode, chyba że określono inaczej.

Jeśli pracujesz nad tym kursem samodzielnie, zapoznaj się z [podstawami środowiska](../envsetup/01_setup.md), aby uzyskać więcej szczegółów.

### Wymagania dotyczące wersji

To szkolenie wymaga Nextflow'a 25.10.2 lub nowszego z włączonym parserem składni v2 (domyślnym od wersji 25.10+).
Jeśli używasz lokalnego lub niestandardowego środowiska, upewnij się, że używasz prawidłowych ustawień, jak opisano [tutaj](../info/nxf_versions.md).

## Przygotuj się do pracy

Po uruchomieniu codespace'a musisz wykonać dwie czynności przed rozpoczęciem szkolenia: ustawić katalog roboczy i przejrzeć dostarczone materiały.

### Ustaw katalog roboczy

Domyślnie codespace otwiera się w głównym folderze wszystkich kursów szkoleniowych.
Na potrzeby tego kursu przejdź do katalogu `nextflow-run/`:

```bash
cd nextflow-run/
```

Następnie ustaw VSCode tak, aby koncentrował się na tym katalogu, dzięki czemu tylko odpowiednie pliki będą wyświetlane na pasku bocznym eksploratora plików:

```bash
code .
```

!!! tip "Wskazówka"

    Jeśli z jakiegokolwiek powodu wyjdziesz z tego katalogu (np. codespace przejdzie w stan uśpienia), zawsze możesz użyć pełnej ścieżki, aby do niego wrócić, zakładając, że pracujesz w środowisku szkoleniowym Github Codespaces:

    ```bash
    cd /workspaces/training/nextflow-run
    ```

### Przejrzyj dostarczone materiały

Możesz przeglądać materiały kursu, używając eksploratora plików po lewej stronie lub polecenia `tree`.
Uruchom poniższe polecenie w terminalu, aby zobaczyć pełną strukturę:

```bash
tree . -L 2
```

??? abstract "Zawartość katalogu"

    ```console
    .
    ├── 1-hello.nf
    ├── 2-inputs.nf
    ├── data
    │   ├── greetings-extended.csv
    │   └── greetings.csv
    ├── main.nf
    ├── modules
    │   ├── collectGreetings.nf
    │   ├── convertToUpper.nf
    │   ├── cowpy.nf
    │   └── sayHello.nf
    └── nextflow.config
    ```

**Pliki `.nf`** to skrypty workflow'ów o rosnącym stopniu złożoności, używane w tej kolejności w trakcie kursu.

**Katalog `data/`** zawiera pliki CSV, których użyjemy od sekcji 2.

**Katalog `modules/`** zawiera definicje procesów używane przez `main.nf`.

**Plik `nextflow.config`** to plik konfiguracyjny, który ustawia minimalne właściwości środowiska. Na razie możesz go zignorować; omówimy go w sekcji 4.

## Lista kontrolna gotowości

Myślisz, że jesteś gotowy, aby zacząć?

- [ ] Rozumiem cel tego kursu i jego wymagania wstępne
- [ ] Moje środowisko jest uruchomione i działa
- [ ] Ustawiłem odpowiednio Swój katalog roboczy

Jeśli możesz zaznaczyć wszystkie pola, jesteś gotowy do działania.

**Aby przejść do [Części 1: Uruchamianie Nextflow'a](./01_run_nextflow.md), kliknij strzałkę w prawym dolnym rogu tej strony.**
