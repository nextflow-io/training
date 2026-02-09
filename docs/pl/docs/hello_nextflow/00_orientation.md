# Pierwsze kroki

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Tłumaczenie wspomagane przez AI - [dowiedz się więcej i zaproponuj poprawki](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

<div class="video-wrapper">
  <iframe width="560" height="315" src="https://www.youtube.com/embed/PIjOdFaYwWA?si=gZxlXgkVxuLEzOsC" title="YouTube video player" frameborder="0" allow="accelerometer; autoplay; clipboard-write; encrypted-media; gyroscope; picture-in-picture; web-share" referrerpolicy="strict-origin-when-cross-origin" allowfullscreen></iframe>
</div>

/// caption
:fontawesome-brands-youtube:{ .youtube } Zobacz [całą playlistę](https://youtube.com/playlist?list=PLPZ8WHdZGxmWKozQuzr27jyMGqp9kElVK&si=eF7cLR62goy-lc6n) na kanale YouTube Nextflow.

:green_book: Transkrypcja wideo jest dostępna [tutaj](./transcripts/00_orientation.md).
///

!!! tip "Wskazówka"

    Filmy na YouTube mają kilka super mocy!

    - :fontawesome-solid-closed-captioning: Wysokiej jakości (ręcznie opracowane) napisy. Włącz je za pomocą ikony :material-subtitles:
    - :material-bookmark: Rozdziały wideo na osi czasu odpowiadające nagłówkom stron.

## Uruchom środowisko szkoleniowe

Aby skorzystać z gotowego środowiska udostępnionego przez nas na GitHub Codespaces, kliknij przycisk „Open in GitHub Codespaces" poniżej. Inne opcje znajdziesz w sekcji [Opcje środowiska](../envsetup/index.md).

Zalecamy otwarcie środowiska szkoleniowego w nowej karcie lub oknie przeglądarki (użyj prawego przycisku myszy, ctrl+klik lub cmd+klik w zależności od Twojego sprzętu), abyś mógł czytać instrukcje podczas ładowania środowiska.
Musisz mieć te instrukcje otwarte równolegle, aby przejść przez kurs.

[![Open in GitHub Codespaces](https://github.com/codespaces/badge.svg)](https://codespaces.new/nextflow-io/training?quickstart=1&ref=master)

### Podstawy środowiska

To środowisko szkoleniowe zawiera wszystkie niezbędne oprogramowanie, kod i dane potrzebne do pracy z kursem, więc nie musisz niczego instalować samodzielnie.

Codespace jest skonfigurowany z interfejsem VSCode, który obejmuje eksplorator systemu plików, edytor kodu i powłokę terminala.
Wszystkie instrukcje podawane podczas kursu (np. „otwórz plik", „edytuj kod" lub „uruchom to polecenie") odnoszą się do tych trzech części interfejsu VSCode, chyba że określono inaczej.

Jeśli pracujesz nad tym kursem samodzielnie, zapoznaj się z [podstawami środowiska](../envsetup/01_setup.md), aby uzyskać więcej szczegółów.

### Wymagania dotyczące wersji

To szkolenie jest przeznaczone dla Nextflow 25.10.2 lub nowszego **z WŁĄCZONYM parserem składni v2**.
Jeśli korzystasz ze środowiska lokalnego lub niestandardowego, upewnij się, że używasz prawidłowych ustawień zgodnie z dokumentacją [tutaj](../info/nxf_versions.md).

## Przygotuj się do pracy

Gdy Twój codespace już działa, musisz wykonać dwie czynności przed rozpoczęciem szkolenia: ustawić katalog roboczy dla tego konkretnego kursu i przyjrzeć się udostępnionym materiałom.

### Ustaw katalog roboczy

Domyślnie codespace otwiera się z katalogiem roboczym ustawionym w katalogu głównym wszystkich kursów szkoleniowych, ale w tym kursie będziemy pracować w katalogu `hello-nextflow/`.

Zmień katalog teraz, uruchamiając to polecenie w terminalu:

```bash
cd hello-nextflow/
```

Możesz ustawić VSCode tak, aby skupił się na tym katalogu, dzięki czemu tylko odpowiednie pliki będą widoczne w pasku bocznym eksploratora plików:

```bash
code .
```

!!! tip "Wskazówka"

    Jeśli z jakiegokolwiek powodu wyjdziesz z tego katalogu (np. Twój codespace przejdzie w stan uśpienia), zawsze możesz użyć pełnej ścieżki, aby do niego wrócić, zakładając, że pracujesz w środowisku szkoleniowym Github Codespaces:

    ```bash
    cd /workspaces/training/hello-nextflow
    ```

Teraz przyjrzyjmy się zawartości.

### Przejrzyj udostępnione materiały

Możesz przeglądać zawartość tego katalogu za pomocą eksploratora plików po lewej stronie obszaru roboczego szkolenia.
Alternatywnie możesz użyć polecenia `tree`.

W trakcie kursu używamy wyniku polecenia `tree` do przedstawienia struktury katalogów i zawartości w czytelnej formie, czasami z niewielkimi modyfikacjami dla przejrzystości.

Tutaj generujemy spis treści do drugiego poziomu w dół:

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
Używamy takich zwijanych sekcji, aby w zwięzły sposób uwzględnić oczekiwane wyniki poleceń.

- **Pliki `.nf`** to skrypty workflow'ów, których nazwy odpowiadają części kursu, w której są używane.

- **Plik `nextflow.config`** to plik konfiguracyjny, który ustawia minimalne właściwości środowiska.
  Na razie możesz go zignorować.

- **Plik `greetings.csv`** w katalogu `data/` zawiera dane wejściowe, których będziemy używać w większości kursu. Jest opisany w Części 2 (Kanały), gdy wprowadzamy go po raz pierwszy.

- **Pliki `test-params.*`** to pliki konfiguracyjne, których użyjemy w Części 6 (Konfiguracja). Na razie możesz je zignorować.

- **Katalog `solutions`** zawiera ukończone skrypty workflow'ów będące wynikiem każdego kroku kursu.
  Są przeznaczone do użycia jako punkt odniesienia do sprawdzenia Twojej pracy i rozwiązywania ewentualnych problemów.

## Lista kontrolna gotowości

Myślisz, że jesteś gotowy, aby zacząć?

- [ ] Rozumiem cel tego kursu i jego wymagania wstępne
- [ ] Moje środowisko działa
- [ ] Ustawiłem odpowiednio mój katalog roboczy

Jeśli możesz zaznaczyć wszystkie pola, możesz zaczynać.

**Aby przejść do [Części 1: Hello World](./01_hello_world.md), kliknij strzałkę w prawym dolnym rogu tej strony.**
