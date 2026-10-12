# Pierwsze kroki

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Tłumaczenie wspomagane przez AI - [dowiedz się więcej i zasugeruj ulepszenia](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

## Uruchomienie środowiska szkoleniowego

Aby skorzystać z gotowego środowiska dostępnego na GitHub Codespaces, kliknij poniższy przycisk „Open in GitHub Codespaces". Inne opcje znajdziesz w sekcji [Konfiguracja środowiska](../envsetup/index.md).

Zalecamy otwarcie środowiska szkoleniowego w nowej karcie lub oknie przeglądarki (użyj prawego przycisku myszy, ctrl+klik lub cmd+klik, w zależności od Twojego sprzętu), żebyś mógł/mogła czytać dalej, podczas gdy środowisko się ładuje.
Instrukcje muszą pozostać otwarte równolegle przez cały czas pracy z kursem.

[![Open in GitHub Codespaces](https://github.com/codespaces/badge.svg)](https://codespaces.new/nextflow-io/training?quickstart=1&ref=master)

### Podstawy środowiska

Środowisko szkoleniowe zawiera całe oprogramowanie, kod i dane niezbędne do ukończenia kursu — nie musisz niczego instalować samodzielnie.

Codespace jest skonfigurowany z interfejsem VSCode, który obejmuje eksplorator plików, edytor kodu oraz terminal z powłoką.
Wszystkie instrukcje podawane w trakcie kursu (np. „otwórz plik", „edytuj kod" lub „uruchom to polecenie") odnoszą się do tych trzech elementów interfejsu VSCode, chyba że zaznaczono inaczej.

Jeśli pracujesz z tym kursem samodzielnie, zapoznaj się z [podstawami środowiska](../envsetup/01_setup.md), gdzie znajdziesz więcej szczegółów.

### Wymagania dotyczące wersji

To szkolenie działa z Nextflow'em w wersji 25.10.2 lub nowszej **z parserem składni v2**, który jest domyślny od Nextflow'a 26.04.
W naszym środowisku szkoleniowym nie musisz nic robić: działa tam Nextflow 26.04.4 z parserem v2. Jeśli korzystasz z lokalnego lub własnego środowiska, zapoznaj się z [informacjami o wersjach](../info/nxf_versions.md).

!!! warning "nf-core/demo wymaga Nextflow'a 25.10.4 lub nowszego"

    Pipeline `nf-core/demo` używany w Części 1 wymusza własną minimalną wersję Nextflow'a (`>=25.10.4`), która jest bardziej restrykcyjna niż ogólne wymaganie szkolenia wynoszące 25.10.2.
    Nasze środowisko szkoleniowe już spełnia ten wymóg; jeśli korzystasz z lokalnego lub własnego środowiska, upewnij się, że masz Nextflow'a w wersji 25.10.4 lub nowszej.

To szkolenie wymaga dodatkowo **nf-core tools w wersji 4.0.2**.
Korzystanie z innej wersji narzędzi nf-core może utrudnić śledzenie materiału.

Zainstalowaną wersję możesz sprawdzić poleceniem `nf-core --version`.

!!! warning "Zgodność z parserem v2"

    Wiele pipeline'ów nf-core nie obsługuje jeszcze parsera składni v2.
    Jeśli uruchomisz pipeline nf-core inny niż te używane w tym kursie i napotkasz błędy, może być konieczne przełączenie się na parser v1 poprzez ustawienie `export NXF_SYNTAX_PARSER=v1`.
    Szczegóły znajdziesz w [informacjach o wersjach](../info/nxf_versions.md).

## Przygotowanie do pracy

Gdy Twój codespace jest już uruchomiony, przed rozpoczęciem szkolenia musisz wykonać dwie czynności: ustawić katalog roboczy dla tego konkretnego kursu oraz zapoznać się z dostarczonymi materiałami.

### Ustawienie katalogu roboczego

Domyślnie codespace otwiera się z katalogiem roboczym ustawionym w katalogu głównym wszystkich kursów szkoleniowych, jednak w tym kursie będziemy pracować w katalogu `nfcore-use/`.

Zmień katalog, uruchamiając w terminalu następujące polecenie:

```bash
cd nfcore-use/
```

!!! tip "Wskazówka"

    Jeśli z jakiegoś powodu opuścisz ten katalog (np. Twój codespace przejdzie w tryb uśpienia), zawsze możesz użyć pełnej ścieżki, aby do niego wrócić — zakładając, że pracujesz w środowisku szkoleniowym GitHub Codespaces:

    ```bash
    cd /workspaces/training/nfcore-use
    ```

Następnie zapoznaj się z zawartością tego katalogu.

### Zapoznanie się z dostarczonymi materiałami

Zawartość katalogu możesz przeglądać za pomocą eksploratora plików po lewej stronie obszaru roboczego.
Alternatywnie możesz użyć polecenia `tree`.

```bash
tree .
```

??? abstract "Zawartość katalogu"

    ```console
    .
    ├── custom.config
    ├── laptop.config
    ├── malformed_samplesheet.csv
    └── my_params.yml
    ```

- **Plik `laptop.config`** to plik konfiguracyjny, którego użyjemy w sekcji 4, aby ograniczyć zużycie zasobów podczas lokalnego uruchamiania pipeline'u w skali produkcyjnej.
  Na razie możesz go zignorować.
- **Pliki `my_params.yml`, `malformed_samplesheet.csv` i `custom.config`** są używane w Części 2 do zademonstrowania ustawiania parametrów z pliku, walidacji danych wejściowych oraz nadpisywania konfiguracji na poziomie procesu.
  Możesz je również zignorować do tego momentu.

## Lista kontrolna gotowości

Myślisz, że jesteś gotowy/gotowa do działania?

- [ ] Rozumiem cel tego kursu i jego wymagania wstępne
- [ ] Moje środowisko jest uruchomione i działa
- [ ] Korzystam z nf-core tools w wersji 4.0.2 (sprawdź poleceniem `nf-core --version`)
- [ ] Ustawiłem/ustawiłam odpowiedni katalog roboczy

Jeśli możesz zaznaczyć wszystkie pola, możesz zaczynać.

**Aby przejść do Części 1, kliknij strzałkę w prawym dolnym rogu tej strony.**
