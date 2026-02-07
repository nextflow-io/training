# Rozpoczęcie pracy

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Tłumaczenie wspomagane przez AI - [dowiedz się więcej i zasugeruj ulepszenia](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

<!--
<div class="video-wrapper">
  <iframe width="560" height="315" src="https://www.youtube.com/embed/G3CV-FcV-rc?si=nyLvwhrSB2m1NPc5&amp;list=PLPZ8WHdZGxmXiHf8B26oB_fTfoKQdhlik" title="YouTube video player" frameborder="0" allow="accelerometer; autoplay; clipboard-write; encrypted-media; gyroscope; picture-in-picture; web-share" referrerpolicy="strict-origin-when-cross-origin" allowfullscreen></iframe>
</div>

/// caption
:fontawesome-brands-youtube:{ .youtube } Zobacz [całą playlistę](https://www.youtube.com/playlist?list=PLPZ8WHdZGxmXiHf8B26oB_fTfoKQdhlik) na kanale YouTube Nextflow.

:green_book: Transkrypcja wideo jest dostępna [tutaj](./transcripts/00_orientation.md).
///

!!! tip "Wskazówka"

    Filmy na YouTube mają specjalne funkcje!

    - :fontawesome-solid-closed-captioning: Wysokiej jakości (ręcznie przygotowane) napisy. Włącz je ikoną :material-subtitles:
    - :material-bookmark: Rozdziały wideo na osi czasu odpowiadające nagłówkom strony.

-->

## Uruchom środowisko szkoleniowe

Aby użyć prekonfigurowanego środowiska, które udostępniamy na GitHub Codespaces, kliknij przycisk "Open in GitHub Codespaces" poniżej. Inne opcje znajdziesz w [Opcjach środowiska](../envsetup/index.md).

Zalecamy otwarcie środowiska szkoleniowego w nowej karcie lub oknie przeglądarki (użyj prawego przycisku myszy, ctrl-click lub cmd-click w zależności od sprzętu), abyś mógł czytać dalej podczas ładowania środowiska.
Będziesz musiał trzymać te instrukcje otwarte równolegle, aby pracować z kursem.

[![Open in GitHub Codespaces](https://github.com/codespaces/badge.svg)](https://codespaces.new/nextflow-io/training?quickstart=1&ref=master)

### Podstawy środowiska

To środowisko szkoleniowe zawiera całe oprogramowanie, kod i dane niezbędne do pracy z kursem szkoleniowym, więc nie musisz niczego instalować samodzielnie.

Codespace jest skonfigurowany z interfejsem VSCode, który obejmuje eksplorator systemu plików, edytor kodu i powłokę terminala.
Wszystkie instrukcje podane podczas kursu (np. 'otwórz plik', 'edytuj kod' lub 'uruchom to polecenie') odnoszą się do tych trzech części interfejsu VSCode, chyba że zaznaczono inaczej.

Jeśli pracujesz z tym kursem samodzielnie, zapoznaj się z [podstawami środowiska](../envsetup/01_setup.md), aby uzyskać więcej szczegółów.

### Wymagania wersji

To szkolenie jest zaprojektowane dla Nextflow'a 25.10.2 lub nowszego **z WŁĄCZONYM parserem składni v2**.
Jeśli używasz lokalnego lub niestandardowego środowiska, upewnij się, że używasz prawidłowych ustawień, jak udokumentowano [tutaj](../info/nxf_versions.md).

## Przygotuj się do pracy

Po uruchomieniu Twojego codespace musisz wykonać dwie czynności zanim zanurzysz się w szkoleniu: ustawić katalog roboczy dla tego konkretnego kursu i przejrzeć dostarczone materiały.

### Ustaw katalog roboczy

Domyślnie codespace otwiera się z główną lokalizacją wszystkich kursów szkoleniowych jako katalogiem roboczym, ale tym razem będziemy pracować w katalogu `hello-nextflow/`.

Zmień teraz katalog, uruchamiając to polecenie w terminalu:

```bash
cd hello-nextflow/
```

Możesz ustawić VSCode, aby skupił się na tym katalogu, tak aby w pasku bocznym eksploratora plików wyświetlały się tylko odpowiednie pliki:

```bash
code .
```

!!! tip "Wskazówka"

    Jeśli z jakiegokolwiek powodu wyjdziesz z tego katalogu (np. Twój codespace przejdzie w stan uśpienia), zawsze możesz użyć pełnej ścieżki, aby do niego wrócić, zakładając, że uruchamiasz to w środowisku szkoleniowym GitHub Codespaces:

    ```bash
    cd /workspaces/training/hello-nextflow
    ```

Teraz przyjrzyjmy się zawartości.

### Przeglądaj dostarczone materiały

Możesz przeglądać zawartość tego katalogu za pomocą eksploratora plików po lewej stronie przestrzeni roboczej szkolenia.
Alternatywnie możesz użyć polecenia `tree`.

W trakcie kursu używamy wyjścia `tree` do reprezentowania struktury katalogów i zawartości w czytelnej formie, czasami z drobnymi modyfikacjami dla przejrzystości.

Tutaj generujemy spis treści do drugiego poziomu w głąb:

```bash
tree . -L 2
```

??? abstract "Zawartość katalogu"

    ```console
    .
    ├── data
    │   └── greetings.csv
    ├── hello-channels.nf
    ├── hello-config.nf
    ├── hello-containers.nf
    ├── hello-modules.nf
    ├── hello-workflow.nf
    ├── hello-world.nf
    ├── nextflow.config
    ├── solutions
    │   ├── 1-hello-world
    │   ├── 2-hello-channels
    │   ├── 3-hello-workflow
    │   ├── 4-hello-modules
    │   ├── 5-hello-containers
    │   └── 6-hello-config
    ├── test-params.json
    └── test-params.yaml
    ```

Kliknij na kolorowe pole, aby rozwinąć sekcję i zobaczyć jej zawartość.
Używamy takich rozwijanych sekcji, aby zwięźle włączać oczekiwane wyjście poleceń.

- **Pliki `.nf`** to skrypty workflow'ów nazwane na podstawie części kursu, w której są używane.

- **Plik `nextflow.config`** to plik konfiguracyjny, który ustawia minimalne właściwości środowiska.
  Na razie możesz go zignorować.

- **Plik `greetings.csv`** w `data/` zawiera dane wejściowe, których użyjemy w większości kursu. Jest opisany w Części 2 (Channels), kiedy wprowadzamy go po raz pierwszy.

- **Pliki `test-params.*`** to pliki konfiguracyjne, których użyjemy w Części 6 (Configuration). Na razie możesz je zignorować.

- **Katalog `solutions`** zawiera ukończone skrypty workflow'ów, które są wynikiem każdego etapu kursu.
  Służą jako odniesienie do weryfikacji Twojej pracy i rozwiązywania problemów.

## Lista kontrolna gotowości

Myślisz, że jesteś gotowy do zanurzenia się?

- [ ] Rozumiem cel tego kursu i jego wymagania wstępne
- [ ] Moje środowisko jest uruchomione
- [ ] Ustawiłem odpowiednio Swój katalog roboczy

Jeśli możesz zaznaczyć wszystkie pola, możesz zaczynać.

**Aby przejść do [Części 1: Hello World](./01_hello_world.md), kliknij strzałkę w prawym dolnym rogu tej strony.**
