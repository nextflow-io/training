# Pierwsze kroki

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Tłumaczenie wspomagane przez AI - [dowiedz się więcej i zasugeruj ulepszenia](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

## Uruchom środowisko szkoleniowe

Aby skorzystać z gotowego środowiska, które udostępniamy na GitHub Codespaces, kliknij przycisk „Open in GitHub Codespaces" poniżej. Inne opcje znajdziesz w sekcji [Opcje środowiska](../../envsetup/index.md).

Zalecamy otwarcie środowiska szkoleniowego w nowej karcie lub oknie przeglądarki (użyj kliknięcia prawym przyciskiem myszy, ctrl+kliknięcie lub cmd+kliknięcie w zależności od Twojego sprzętu), abyś mógł/mogła czytać instrukcje podczas ładowania środowiska.
Musisz mieć te instrukcje otwarte równolegle, aby przejść przez kurs.

[![Open in GitHub Codespaces](https://github.com/codespaces/badge.svg)](https://codespaces.new/nextflow-io/training?quickstart=1&ref=master)

### Podstawy środowiska

To środowisko szkoleniowe zawiera wszystkie niezbędne oprogramowanie, kod i dane potrzebne do przejścia przez kurs, więc nie musisz niczego instalować samodzielnie.

Codespace jest skonfigurowany z interfejsem VSCode, który obejmuje eksplorator systemu plików, edytor kodu i terminal powłoki.
Wszystkie instrukcje podawane podczas kursu (np. „otwórz plik", „edytuj kod" lub „uruchom to polecenie") odnoszą się do tych trzech części interfejsu VSCode, chyba że zaznaczono inaczej.

Jeśli przechodzisz przez ten kurs samodzielnie, zapoznaj się z [podstawami środowiska](../../envsetup/01_setup.md), aby uzyskać więcej szczegółów.

### Wymagania dotyczące wersji

To szkolenie jest przeznaczone dla Nextflow'a w wersji 25.10.2 lub nowszej **z WŁĄCZONYM parserem składni v2**.
Jeśli korzystasz z lokalnego lub niestandardowego środowiska, upewnij się, że używasz prawidłowych ustawień zgodnie z dokumentacją [tutaj](../../info/nxf_versions.md).

## Przygotuj się do pracy

Gdy Twój codespace już działa, musisz wykonać dwie czynności przed rozpoczęciem szkolenia: ustawić katalog roboczy dla tego konkretnego kursu i przyjrzeć się udostępnionym materiałom.

### Ustaw katalog roboczy

Domyślnie codespace otwiera się z katalogiem roboczym ustawionym w katalogu głównym wszystkich kursów szkoleniowych, ale w tym kursie będziemy pracować w katalogu `nf4-science/genomics/`.

Zmień teraz katalog, uruchamiając to polecenie w terminalu:

```bash
cd nf4-science/genomics/
```

Możesz ustawić VSCode tak, aby skupił się na tym katalogu, dzięki czemu tylko odpowiednie pliki będą widoczne w pasku bocznym eksploratora plików:

```bash
code .
```

!!! tip "Wskazówka"

    Jeśli z jakiegokolwiek powodu wyjdziesz z tego katalogu (np. Twój codespace przejdzie w stan uśpienia), zawsze możesz użyć pełnej ścieżki, aby do niego wrócić, zakładając, że pracujesz w środowisku szkoleniowym Github Codespaces:

    ```bash
    cd /workspaces/training/nf4-science/genomics
    ```

Teraz przyjrzyjmy się zawartości.

### Przejrzyj udostępnione materiały

Możesz przeglądać zawartość tego katalogu, korzystając z eksploratora plików po lewej stronie obszaru roboczego szkolenia.
Alternatywnie możesz użyć polecenia `tree`.

W trakcie kursu używamy wyniku polecenia `tree` do przedstawienia struktury katalogów i ich zawartości w czytelnej formie, czasami z niewielkimi modyfikacjami dla przejrzystości.

Tutaj generujemy spis treści do drugiego poziomu w dół:

```bash
tree . -L 2
```

??? abstract "Zawartość katalogu"

    ```console
    .
    ├── data
    │   ├── bam
    │   ├── ref
    │   ├── sample_bams.txt
    │   └── samplesheet.csv
    ├── genomics.nf
    ├── modules
    │   ├── gatk_haplotypecaller.nf
    │   └── samtools_index.nf
    ├── nextflow.config
    └── solutions
        ├── modules
        ├── nf-test.config
        ├── part2
        └── tests

    8 directories, 8 files
    ```

Kliknij na kolorowe pole, aby rozwinąć sekcję i zobaczyć jej zawartość.
Używamy takich zwijanych sekcji, aby wyświetlać oczekiwane wyniki poleceń, a także zawartość katalogów i plików w zwięzły sposób.

- **Plik `genomics.nf`** to skrypt workflow'a, który będziesz rozbudowywać w trakcie kursu.

- **Katalog `modules`** zawiera szkieletowe pliki modułów, które wypełnisz podczas kursu.

- **Plik `nextflow.config`** to plik konfiguracyjny, który ustawia minimalne właściwości środowiska.
  Na razie możesz go zignorować.

- **Katalog `data`** zawiera dane wejściowe i powiązane zasoby, opisane później w kursie.

- **Katalog `solutions`** zawiera ukończone pliki modułów oraz rozwiązanie Części 2, które może służyć jako punkt wyjścia do Części 3.
  Są one przeznaczone do użycia jako punkt odniesienia do sprawdzenia Twojej pracy i rozwiązywania ewentualnych problemów.

## Lista kontrolna gotowości

Myślisz, że jesteś gotowy/gotowa, aby zacząć?

- [ ] Rozumiem cel tego kursu i jego wymagania wstępne
- [ ] Moje środowisko działa
- [ ] Ustawiłem/ustawiłam odpowiednio mój katalog roboczy

Jeśli możesz zaznaczyć wszystkie pola, możesz zaczynać.

**Aby przejść do [Części 1: Przegląd metody i testowanie ręczne](./01_method.md), kliknij strzałkę w prawym dolnym rogu tej strony.**
