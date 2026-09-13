# Pierwsze kroki

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Tłumaczenie wspomagane przez AI - [dowiedz się więcej i zasugeruj ulepszenia](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

## Uruchomienie środowiska szkoleniowego

Aby skorzystać z gotowego środowiska dostępnego na GitHub Codespaces, kliknij przycisk „Open in GitHub Codespaces" poniżej. Inne opcje znajdziesz w sekcji [Opcje środowiska](../envsetup/index.md).

Zalecamy otwarcie środowiska szkoleniowego w nowej karcie lub oknie przeglądarki (użyj prawego przycisku myszy, ctrl+klik lub cmd+klik, w zależności od Twojego sprzętu), żebyś mógł czytać dalej, podczas gdy środowisko się ładuje.
Podczas pracy z kursem trzymaj te instrukcje otwarte równolegle.

[![Open in GitHub Codespaces](https://github.com/codespaces/badge.svg)](https://codespaces.new/nextflow-io/training?quickstart=1&ref=master)

### Podstawy środowiska

Środowisko szkoleniowe zawiera całe oprogramowanie, kod i dane potrzebne do ukończenia kursu — nie musisz niczego instalować samodzielnie.

Codespace jest skonfigurowany z interfejsem VSCode, który obejmuje eksplorator plików, edytor kodu oraz terminal powłoki.
Wszystkie instrukcje podawane w trakcie kursu (np. „otwórz plik", „edytuj kod" lub „uruchom to polecenie") odnoszą się do tych trzech elementów interfejsu VSCode, chyba że zaznaczono inaczej.

Jeśli pracujesz z tym kursem samodzielnie, zapoznaj się z [podstawami środowiska](../envsetup/01_setup.md), gdzie znajdziesz więcej szczegółów.

### Wymagania dotyczące wersji

Ten kurs wymaga Nextflow'a w wersji 25.10.2 lub nowszej, z włączonym parserem składni v2 (domyślnym od wersji 25.10+).
Jeśli korzystasz z lokalnego lub niestandardowego środowiska, upewnij się, że używasz właściwych ustawień opisanych [tutaj](../info/nxf_versions.md).

## Przygotowanie do pracy

Gdy Twój codespace jest już uruchomiony, przed rozpoczęciem pracy należy zrobić dwie rzeczy: ustawić katalog roboczy i zapoznać się z dostarczonymi materiałami.

### Ustawienie katalogu roboczego

Domyślnie codespace otwiera się w katalogu głównym wszystkich kursów szkoleniowych.
Na potrzeby tego kursu przejdź do katalogu `execution-config/`:

```bash
cd execution-config/
```

Następnie skieruj VSCode na ten katalog, żeby w pasku bocznym eksploratora plików wyświetlały się tylko odpowiednie pliki:

```bash
code .
```

!!! tip "Wskazówka"

    Jeśli z jakiegoś powodu opuścisz ten katalog (np. Twój codespace przejdzie w tryb uśpienia), zawsze możesz użyć pełnej ścieżki, aby do niego wrócić — zakładając, że pracujesz w środowisku szkoleniowym GitHub Codespaces:

    ```bash
    cd /workspaces/training/execution-config
    ```

### Zapoznanie się z dostarczonymi materiałami

Materiały kursu możesz przeglądać za pomocą eksploratora plików po lewej stronie lub polecenia `tree`.
Uruchom poniższe polecenie w terminalu, aby zobaczyć pełną strukturę:

```bash
tree . -L 2
```

??? abstract "Zawartość katalogu"

    ```console
    .
    ├── data
    │   └── greetings.csv
    ├── main.nf
    ├── modules
    │   ├── collectGreetings.nf
    │   ├── convertToUpper.nf
    │   ├── cowpy.nf
    │   └── sayHello.nf
    └── nextflow.config
    ```

Pliki **`main.nf`** i **`modules/`** to ten sam wieloetapowy pipeline co w [Nextflow Run](../nextflow_run/index.md), a plik **`nextflow.config`** to ta sama konfiguracja, którą już tam widziałeś.
W trakcie ćwiczeń rozszerzysz oba.

Katalog **`data/`** zawiera plik wejściowy CSV, z którego korzysta pipeline.

## Lista kontrolna gotowości

Myślisz, że jesteś gotowy, żeby zacząć?

- [ ] Rozumiem cel tego kursu i jego wymagania wstępne
- [ ] Moje środowisko jest uruchomione i działa
- [ ] Ustawiłem odpowiedni katalog roboczy

Jeśli możesz zaznaczyć wszystkie pola, możesz śmiało ruszać dalej.

**Aby przejść do [Części 1: Dostosowanie do środowiska obliczeniowego](./01_packaging_and_execution.md), kliknij strzałkę w prawym dolnym rogu tej strony.**
