# Pierwsze kroki

## Uruchomienie środowiska szkoleniowego

Aby skorzystać z przygotowanego przez nas środowiska w GitHub Codespaces, kliknij przycisk "Open in GitHub Codespaces" poniżej. Aby poznać inne opcje, zobacz [Opcje środowiska](../envsetup/index.md).

Zalecamy otwarcie środowiska szkoleniowego w nowej karcie lub oknie przeglądarki (użyj kliknięcia prawym przyciskiem myszy, ctrl+kliknięcie lub cmd+kliknięcie w zależności od sprzętu), aby móc czytać instrukcje podczas ładowania środowiska.
Musisz zachować te instrukcje otwarte równolegle, aby przejść przez kurs.

[![Open in GitHub Codespaces](https://github.com/codespaces/badge.svg)](https://codespaces.new/nextflow-io/training?quickstart=1&ref=master)

### Podstawy środowiska

To środowisko szkoleniowe zawiera wszystkie oprogramowanie, kod i dane niezbędne do przejścia przez kurs szkoleniowy, więc nie musisz niczego samodzielnie instalować.

Codespace jest skonfigurowany z interfejsem VSCode, który zawiera eksplorator systemu plików, edytor kodu i powłokę terminala.
Wszystkie instrukcje podane podczas kursu (np. „otwórz plik", „edytuj kod" lub „uruchom to polecenie") odnoszą się do tych trzech części interfejsu VSCode, chyba że określono inaczej.

Jeśli przechodzisz przez ten kurs samodzielnie, zapoznaj się z [podstawami środowiska](../envsetup/01_setup.md), aby uzyskać więcej szczegółów.

### Wymagania dotyczące wersji

To szkolenie jest przeznaczone dla **Nextflow 25.10.2** lub nowszego **z WYŁĄCZONYM parserem składni v2**.

#### Jeśli korzystasz z naszego środowiska szkoleniowego:

MUSISZ uruchomić następujące polecenie przed kontynuowaniem:

```bash
export NXF_SYNTAX_PARSER=v1
```

#### Jeśli korzystasz ze środowiska lokalnego lub niestandardowego:

Upewnij się, że używasz prawidłowych ustawień zgodnie z dokumentacją [tutaj](../info/nxf_versions.md).

Szkolenie dodatkowo wymaga **nf-core tools 3.4.1**.
Jeśli używasz innej wersji narzędzi nf-core, możesz mieć trudności z podążaniem za kursem.

Możesz sprawdzić, jaka wersja jest zainstalowana w Twoim środowisku za pomocą polecenia `nf-core --version`.

## Przygotowanie do pracy

Gdy Twój codespace już działa, musisz wykonać dwie rzeczy przed rozpoczęciem szkolenia: ustawić katalog roboczy dla tego konkretnego kursu i przyjrzeć się dostarczonym materiałom.

### Ustawienie katalogu roboczego

Domyślnie codespace otwiera się z katalogiem roboczym ustawionym na katalog główny wszystkich kursów szkoleniowych, ale w tym kursie będziemy pracować w katalogu `hello-nf-core/`.

Zmień teraz katalog, uruchamiając to polecenie w terminalu:

```bash
cd hello-nf-core/
```

!!! tip "Wskazówka"

    Jeśli z jakiegokolwiek powodu wyjdziesz z tego katalogu (np. Twój codespace przejdzie w stan uśpienia), zawsze możesz użyć pełnej ścieżki, aby do niego powrócić, zakładając, że pracujesz w środowisku szkoleniowym Github Codespaces:

    ```bash
    cd /workspaces/training/hello-nf-core
    ```

Teraz przyjrzyjmy się zawartości tego katalogu.

### Eksploracja dostarczonych materiałów

Możesz przeglądać zawartość tego katalogu używając eksploratora plików po lewej stronie obszaru roboczego szkolenia.
Alternatywnie możesz użyć polecenia `tree`.

W trakcie kursu używamy wyjścia `tree` do przedstawienia struktury katalogów i ich zawartości w czytelnej formie, czasami z drobnymi modyfikacjami dla przejrzystości.

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

Kliknij na kolorowe pole, aby rozwinąć sekcję i wyświetlić jej zawartość.
Używamy zwijanych sekcji w ten sposób, aby w zwięzły sposób uwzględnić oczekiwane wyjście poleceń.

- **Plik `greetings.csv`** to plik CSV zawierający minimalne dane kolumnowe, których używamy do celów testowych.

- **Katalog `original-hello`** zawiera kopię kodu źródłowego powstałego w wyniku przejścia przez kompletną serię szkoleń Hello Nextflow (z włączonym Docker).

- **Katalog `solutions`** zawiera ukończone skrypty workflow'ów powstałe w wyniku każdego kroku kursu.
  Są one przeznaczone do użycia jako punkt odniesienia do sprawdzenia Twojej pracy i rozwiązywania problemów.

## Lista gotowości

Myślisz, że jesteś gotowy, aby zacząć?

- [ ] Rozumiem cel tego kursu i jego wymagania wstępne
- [ ] Moje środowisko jest uruchomione i działa
- [ ] Upewniłem się, że parser składni jest ustawiony na **v1**
- [ ] Ustawiłem odpowiednio mój katalog roboczy

Jeśli możesz zaznaczyć wszystkie pola, możesz zaczynać.

**Aby przejść do Części 1, kliknij strzałkę w prawym dolnym rogu tej strony.**
