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

To szkolenie jest przeznaczone dla Nextflow 25.10.2 lub nowszego **z WŁĄCZONYM parserem składni v2**.
Jeśli używasz lokalnego lub niestandardowego środowiska, upewnij się, że używasz prawidłowych ustawień, jak opisano [tutaj](../info/nxf_versions.md).

## Przygotuj się do pracy

Po uruchomieniu codespace'a musisz wykonać dwie czynności przed rozpoczęciem szkolenia: ustawić katalog roboczy dla tego konkretnego kursu i przejrzeć dostarczone materiały.

### Ustaw katalog roboczy

Domyślnie codespace otwiera się z katalogiem roboczym ustawionym w głównym folderze wszystkich kursów szkoleniowych, ale w tym kursie będziemy pracować w katalogu `nextflow-run/`.

Zmień katalog teraz, uruchamiając to polecenie w terminalu:

```bash
cd nextflow-run/
```

Możesz ustawić VSCode tak, aby koncentrował się na tym katalogu, dzięki czemu tylko odpowiednie pliki będą wyświetlane na pasku bocznym eksploratora plików:

```bash
code .
```

!!! tip "Wskazówka"

    Jeśli z jakiegokolwiek powodu wyjdziesz z tego katalogu (np. codespace przejdzie w stan uśpienia), zawsze możesz użyć pełnej ścieżki, aby do niego wrócić, zakładając, że pracujesz w środowisku szkoleniowym Github Codespaces:

    ```bash
    cd /workspaces/training/nextflow-run
    ```

Teraz przyjrzyjmy się zawartości.

### Przejrzyj dostarczone materiały

Możesz przeglądać zawartość tego katalogu, używając eksploratora plików po lewej stronie obszaru roboczego szkolenia.
Alternatywnie możesz użyć polecenia `tree`.

W całym kursie używamy wyjścia `tree` do przedstawienia struktury i zawartości katalogów w czytelnej formie, czasami z drobnymi modyfikacjami dla przejrzystości.

Tutaj generujemy spis treści do drugiego poziomu:

```bash
tree . -L 2
```

??? abstract "Zawartość katalogu"

    ```console
    .
    ├── 1-hello.nf
    ├── 2a-inputs.nf
    ├── 2b-multistep.nf
    ├── 2c-modules.nf
    ├── 2d-container.nf
    ├── 3-main.nf
    ├── data
    │   └── greetings.csv
    ├── modules
    │   ├── collectGreetings.nf
    │   ├── convertToUpper.nf
    │   ├── cowpy.nf
    │   └── sayHello.nf
    ├── nextflow.config
    ├── solutions
    │   ├── 3-main.nf
    │   ├── modules
    │   └── nextflow.config
    ├── test-params.json
    └── test-params.yaml
    ```

Kliknij na kolorowe pole, aby rozwinąć sekcję i zobaczyć jej zawartość.
Używamy zwijanych sekcji takich jak ta do wyświetlania oczekiwanego wyjścia poleceń oraz zawartości katalogów i plików w zwięzły sposób.

- **Pliki `.nf`** to skrypty workflow'ów ponumerowane według części kursu, w której są używane.

- **Plik `nextflow.config`** to plik konfiguracyjny, który ustawia minimalne właściwości środowiska.
  Na razie możesz go zignorować.

- **Plik `greetings.csv`** w katalogu `data/` zawiera dane wejściowe, których użyjemy w większości kursu. Jest opisany w Części 2 (Uruchamianie pipeline'ów), kiedy wprowadzamy go po raz pierwszy.

- **Pliki `test-params.*`** to pliki konfiguracyjne, których użyjemy w Części 3 (Konfiguracja). Na razie możesz je zignorować.

- **Katalog `solutions`** zawiera końcowy stan workflow'u i jego plików pomocniczych (config i moduły), które są wynikiem ukończenia kursu.
  Są przeznaczone do użycia jako odniesienie do sprawdzenia Twojej pracy i rozwiązywania problemów.

## Lista kontrolna gotowości

Myślisz, że jesteś gotowy, aby zacząć?

- [ ] Rozumiem cel tego kursu i jego wymagania wstępne
- [ ] Moje środowisko jest uruchomione i działa
- [ ] Ustawiłem odpowiednio Swój katalog roboczy

Jeśli możesz zaznaczyć wszystkie pola, jesteś gotowy do działania.

**Aby przejść do [Części 1: Uruchamianie podstawowych operacji](./01_basics.md), kliknij strzałkę w prawym dolnym rogu tej strony.**
