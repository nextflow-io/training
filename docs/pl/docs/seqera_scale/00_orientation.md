# Pierwsze kroki

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Tłumaczenie wspomagane przez AI - [dowiedz się więcej i zasugeruj ulepszenia](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

## Uruchomienie środowiska szkoleniowego

Aby skorzystać z gotowego środowiska dostępnego na GitHub Codespaces, kliknij przycisk „Open in GitHub Codespaces" poniżej. Inne opcje znajdziesz w sekcji [Opcje środowiska](../envsetup/index.md).

Zalecamy otwarcie środowiska szkoleniowego w nowej karcie lub oknie przeglądarki (użyj prawego przycisku myszy, ctrl+klik lub cmd+klik, w zależności od Twojego sprzętu), żebyś mógł/mogła czytać dalej, podczas gdy środowisko się ładuje.
Instrukcje muszą pozostać otwarte równolegle przez cały czas pracy z kursem.

[![Open in GitHub Codespaces](https://github.com/codespaces/badge.svg)](https://codespaces.new/nextflow-io/training?quickstart=1&ref=master)

### Podstawy środowiska

Środowisko szkoleniowe zawiera całe oprogramowanie, kod i dane niezbędne do ukończenia kursu — nie musisz niczego instalować samodzielnie.

Codespace jest skonfigurowany z interfejsem VSCode, który obejmuje eksplorator plików, edytor kodu oraz terminal powłoki.
Wszystkie instrukcje podawane w trakcie kursu (np. „otwórz plik", „edytuj kod" lub „uruchom to polecenie") odnoszą się do tych trzech elementów interfejsu VSCode, chyba że zaznaczono inaczej.

Jeśli pracujesz z tym kursem samodzielnie, zapoznaj się z [podstawami środowiska](../envsetup/01_setup.md), gdzie znajdziesz więcej szczegółów.

## Przygotowanie do pracy

Gdy Twój codespace jest już uruchomiony, przed rozpoczęciem pracy należy wykonać dwie rzeczy: ustawić katalog roboczy i zapoznać się z dostarczonymi materiałami.

### Ustawienie katalogu roboczego

Domyślnie codespace otwiera się w katalogu głównym wszystkich kursów szkoleniowych.
Na potrzeby tego kursu przejdź do katalogu `seqera-scale/`:

```bash
cd seqera-scale/
```

Następnie skonfiguruj VSCode tak, aby skupiał się na tym katalogu — dzięki temu na pasku bocznym eksploratora plików będą widoczne tylko odpowiednie pliki:

```bash
code .
```

!!! tip "Wskazówka"

    Jeśli z jakiegoś powodu opuścisz ten katalog (np. Twój codespace przejdzie w tryb uśpienia), zawsze możesz wrócić do niego, używając pełnej ścieżki — zakładając, że pracujesz w środowisku szkoleniowym GitHub Codespaces:

    ```bash
    cd /workspaces/training/seqera-scale
    ```

### Zapoznanie się z dostarczonymi materiałami

Materiały kursu możesz przeglądać za pomocą eksploratora plików po lewej stronie lub polecenia `tree`.
Uruchom poniższe polecenie w terminalu, aby zobaczyć pełną strukturę:

```bash
tree -a .
```

??? abstract "Zawartość katalogu"

    ```console
    .
    └── .seqera_config
    ```

Plik **`.seqera_config`** to szablon, który wypełnisz w sekcji 3, aby skonfigurować CLI `tw` z Twoim tokenem dostępu i przestrzenią roboczą Seqera.

## Lista kontrolna gotowości

Gotowy/Gotowa do działania?

- [ ] Rozumiem cel tego kursu i jego wymagania wstępne
- [ ] Moje środowisko jest uruchomione i działa
- [ ] Ustawiłem/Ustawiłam odpowiedni katalog roboczy

Jeśli możesz zaznaczyć wszystkie pola, możesz zaczynać.

**Aby przejść do [Części 1: Uruchamianie pipeline'ów z interfejsu webowego](./01_run_with_seqera.md), kliknij strzałkę w prawym dolnym rogu tej strony.**
