# Pierwsze kroki

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Tłumaczenie wspomagane przez AI - [dowiedz się więcej i zasugeruj ulepszenia](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

## Uruchom środowisko szkoleniowe

Aby skorzystać z gotowego środowiska, które udostępniamy na GitHub Codespaces, kliknij poniższy przycisk „Open in GitHub Codespaces". Inne opcje znajdziesz w sekcji [Opcje środowiska](../envsetup/index.md).

Zalecamy otwarcie środowiska szkoleniowego w nowej karcie lub oknie przeglądarki (użyj prawego przycisku myszy, ctrl+klik lub cmd+klik w zależności od Twojego sprzętu), abyś mógł/mogła czytać instrukcje podczas ładowania środowiska.
Musisz mieć te instrukcje otwarte równolegle, aby przejść przez kurs.

[![Open in GitHub Codespaces](https://github.com/codespaces/badge.svg)](https://codespaces.new/nextflow-io/training?quickstart=1&ref=master)

### Podstawy środowiska

To środowisko szkoleniowe zawiera wszystkie niezbędne oprogramowanie, kod i dane potrzebne do przejścia przez kurs, więc nie musisz niczego instalować samodzielnie.

Codespace jest skonfigurowany z interfejsem VSCode, który obejmuje eksplorator systemu plików, edytor kodu i terminal powłoki.
Wszystkie instrukcje podawane podczas kursu (np. „otwórz plik", „edytuj kod" lub „uruchom to polecenie") odnoszą się do tych trzech części interfejsu VSCode, chyba że zaznaczono inaczej.

Jeśli przechodzisz przez ten kurs samodzielnie, zapoznaj się z [podstawami środowiska](../envsetup/01_setup.md), aby uzyskać więcej szczegółów.

### Wymagania dotyczące wersji

To szkolenie jest przeznaczone dla **Nextflow 25.10.2** lub nowszego **z WYŁĄCZONYM parserem składni v2**.

#### Jeśli korzystasz z naszego środowiska szkoleniowego:

MUSISZ uruchomić następujące polecenie przed dalszym postępowaniem:

```bash
export NXF_SYNTAX_PARSER=v1
```

#### Jeśli korzystasz ze środowiska lokalnego lub niestandardowego:

Upewnij się, że używasz prawidłowych ustawień zgodnie z dokumentacją [tutaj](../info/nxf_versions.md).

Szkolenie wymaga dodatkowo **nf-core tools 3.4.1**.
Jeśli używasz innej wersji narzędzi nf-core, możesz mieć trudności z śledzeniem kursu.

Możesz sprawdzić, jaka wersja jest zainstalowana w Twoim środowisku, używając polecenia `nf-core --version`.

## Przygotuj się do pracy

Gdy Twój codespace już działa, musisz wykonać dwie czynności przed rozpoczęciem szkolenia: ustawić katalog roboczy dla tego konkretnego kursu i przyjrzeć się udostępnionym materiałom.

### Ustaw katalog roboczy

Domyślnie codespace otwiera się z katalogiem roboczym ustawionym w katalogu głównym wszystkich kursów szkoleniowych, ale w tym kursie będziemy pracować w katalogu `hello-nf-core/`.

Zmień katalog teraz, uruchamiając to polecenie w terminalu:

```bash
cd hello-nf-core/
```

!!! tip "Wskazówka"

    Jeśli z jakiegokolwiek powodu wyjdziesz z tego katalogu (np. Twój codespace przejdzie w stan uśpienia), zawsze możesz użyć pełnej ścieżki, aby do niego wrócić, zakładając, że pracujesz w środowisku szkoleniowym Github Codespaces:

    ```bash
    cd /workspaces/training/hello-nf-core
    ```

Teraz przyjrzyjmy się zawartości tego katalogu.

### Przejrzyj udostępnione materiały

Możesz przeglądać zawartość tego katalogu, używając eksploratora plików po lewej stronie obszaru roboczego szkolenia.
Alternatywnie możesz użyć polecenia `tree`.

W trakcie kursu używamy wyniku polecenia `tree` do przedstawienia struktury katalogów i ich zawartości w czytelnej formie, czasami z drobnymi modyfikacjami dla przejrzystości.

Tutaj generujemy spis treści do drugiego poziomu w dół:

```bash
tree . -L 2
```

??? abstract "Zawartość katalogu"

    ```console
    .
    ├── greetings.csv
    ├── original-hello
    │   ├── hello.nf
    │   ├── modules
    │   └── nextflow.config
    └── solutions
        ├── composable-hello
        ├── core-hello-part2
        ├── core-hello-part3
        ├── core-hello-part4
        ├── core-hello-part5
        └── core-hello-start
    ```

Kliknij na kolorowe pole, aby rozwinąć sekcję i zobaczyć jej zawartość.
Używamy takich zwijanych sekcji, aby w zwięzły sposób uwzględnić oczekiwane wyniki poleceń.

- **Plik `greetings.csv`** to plik CSV zawierający minimalne dane kolumnowe, których używamy do celów testowych.

- **Katalog `original-hello`** zawiera kopię kodu źródłowego powstałego w wyniku przejścia przez kompletną serię szkoleń Hello Nextflow (z włączonym Dockerem).

- **Katalog `solutions`** zawiera ukończone skrypty workflow'ów, które powstają w wyniku każdego kroku kursu.
  Są one przeznaczone do użycia jako punkt odniesienia do sprawdzenia Twojej pracy i rozwiązywania ewentualnych problemów.

## Lista kontrolna gotowości

Myślisz, że jesteś gotowy/gotowa do rozpoczęcia?

- [ ] Rozumiem cel tego kursu i jego wymagania wstępne
- [ ] Moje środowisko jest uruchomione i działa
- [ ] Upewniłem/upewniłam się, że parser składni jest ustawiony na **v1**
- [ ] Ustawiłem/ustawiłam odpowiednio mój katalog roboczy

Jeśli możesz zaznaczyć wszystkie pola, możesz zaczynać.

**Aby przejść do Części 1, kliknij strzałkę w prawym dolnym rogu tej strony.**
