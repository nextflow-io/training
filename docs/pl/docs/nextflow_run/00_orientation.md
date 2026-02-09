# Pierwsze kroki

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Tłumaczenie wspomagane przez AI - [dowiedz się więcej i zaproponuj poprawki](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

## Uruchom środowisko szkoleniowe

Aby skorzystać z gotowego środowiska udostępnionego przez nas na GitHub Codespaces, kliknij przycisk „Open in GitHub Codespaces" poniżej. Inne opcje znajdziesz w sekcji [Opcje środowiska](../envsetup/index.md).

Zalecamy otwarcie środowiska szkoleniowego w nowej karcie lub oknie przeglądarki (użyj prawego przycisku myszy, ctrl+klik lub cmd+klik w zależności od Twojego sprzętu), abyś mógł czytać instrukcje podczas ładowania środowiska.
Musisz mieć te instrukcje otwarte równolegle, aby przejść przez kurs.

[![Open in GitHub Codespaces](https://github.com/codespaces/badge.svg)](https://codespaces.new/nextflow-io/training?quickstart=1&ref=master)

### Podstawy środowiska

To środowisko szkoleniowe zawiera wszystkie niezbędne oprogramowanie, kod i dane potrzebne do przejścia przez kurs, więc nie musisz niczego instalować samodzielnie.

Codespace jest skonfigurowany z interfejsem VSCode, który obejmuje eksplorator systemu plików, edytor kodu i powłokę terminala.
Wszystkie instrukcje podane w trakcie kursu (np. „otwórz plik", „edytuj kod" lub „uruchom to polecenie") odnoszą się do tych trzech części interfejsu VSCode, chyba że zaznaczono inaczej.

Jeśli przechodzisz przez ten kurs samodzielnie, zapoznaj się z [podstawami środowiska](../envsetup/01_setup.md), aby uzyskać więcej szczegółów.

### Wymagania dotyczące wersji

To szkolenie jest przeznaczone dla Nextflow'a w wersji 25.10.2 lub nowszej **z WŁĄCZONYM parserem składni v2**.
Jeśli korzystasz z lokalnego lub niestandardowego środowiska, upewnij się, że używasz prawidłowych ustawień zgodnie z dokumentacją [tutaj](../info/nxf_versions.md).

## Przygotuj się do pracy

Gdy Twój codespace już działa, musisz wykonać dwie czynności przed rozpoczęciem szkolenia: ustawić katalog roboczy dla tego konkretnego kursu i przyjrzeć się udostępnionym materiałom.

### Ustaw katalog roboczy

Domyślnie codespace otwiera się z katalogiem roboczym ustawionym w katalogu głównym wszystkich kursów szkoleniowych, ale w tym kursie będziemy pracować w katalogu `nextflow-run/`.

Zmień katalog teraz, uruchamiając to polecenie w terminalu:

```bash
cd nextflow-run/
```

Możesz ustawić VSCode tak, aby skupił się na tym katalogu, dzięki czemu tylko odpowiednie pliki będą widoczne na pasku bocznym eksploratora plików:

```bash
code .
```

!!! tip "Wskazówka"

    Jeśli z jakiegokolwiek powodu wyjdziesz z tego katalogu (np. Twój codespace przejdzie w stan uśpienia), zawsze możesz użyć pełnej ścieżki, aby do niego wrócić, zakładając, że pracujesz w środowisku szkoleniowym Github Codespaces:

    ```bash
    cd /workspaces/training/nextflow-run
    ```

Teraz przyjrzyjmy się zawartości.

### Poznaj udostępnione materiały

Możesz przeglądać zawartość tego katalogu, używając eksploratora plików po lewej stronie obszaru roboczego szkolenia.
Alternatywnie możesz użyć polecenia `tree`.

W trakcie kursu używamy wyjścia polecenia `tree` do przedstawienia struktury katalogów i ich zawartości w czytelnej formie, czasami z drobnymi modyfikacjami dla przejrzystości.

Tutaj generujemy spis treści do drugiego poziomu w dół:

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
Używamy takich zwijanych sekcji, aby wyświetlać oczekiwane wyjście poleceń, a także zawartość katalogów i plików w zwięzły sposób.

- **Pliki `.nf`** to skrypty workflow'ów, które są ponumerowane zgodnie z częścią kursu, w której są używane.

- **Plik `nextflow.config`** to plik konfiguracyjny, który ustawia minimalne właściwości środowiska.
  Możesz go na razie zignorować.

- **Plik `greetings.csv`** w katalogu `data/` zawiera dane wejściowe, których będziemy używać w większości kursu. Jest opisany w Części 2 (Uruchamianie pipeline'ów), gdy wprowadzamy go po raz pierwszy.

- **Pliki `test-params.*`** to pliki konfiguracyjne, których użyjemy w Części 3 (Konfiguracja). Możesz je na razie zignorować.

- **Katalog `solutions`** zawiera końcowy stan workflow'a i jego plików pomocniczych (konfiguracji i modułów), które powstają po ukończeniu kursu.
  Są przeznaczone do użycia jako punkt odniesienia do sprawdzenia Twojej pracy i rozwiązywania ewentualnych problemów.

## Lista kontrolna gotowości

Myślisz, że jesteś gotowy, aby zacząć?

- [ ] Rozumiem cel tego kursu i jego wymagania wstępne
- [ ] Moje środowisko działa
- [ ] Ustawiłem odpowiednio mój katalog roboczy

Jeśli możesz zaznaczyć wszystkie pola, możesz zaczynać.

**Aby przejść do [Części 1: Uruchamianie podstawowych operacji](./01_basics.md), kliknij strzałkę w prawym dolnym rogu tej strony.**
