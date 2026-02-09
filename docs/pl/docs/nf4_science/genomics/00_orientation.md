# Pierwsze kroki

## Uruchom środowisko szkoleniowe

Aby skorzystać z gotowego środowiska, które udostępniamy w GitHub Codespaces, kliknij przycisk „Otwórz w GitHub Codespaces" poniżej. Inne opcje znajdziesz w sekcji [Opcje środowiska](../../envsetup/index.md).

Zalecamy otwarcie środowiska szkoleniowego w nowej karcie lub oknie przeglądarki (użyj prawego przycisku myszy, ctrl+klik lub cmd+klik w zależności od Twojego sprzętu), abyś mógł/mogła czytać instrukcje podczas ładowania środowiska.
Będziesz potrzebować tych instrukcji otwartych równolegle, aby przejść przez kurs.

[![Otwórz w GitHub Codespaces](https://github.com/codespaces/badge.svg)](https://codespaces.new/nextflow-io/training?quickstart=1&ref=master)

### Podstawy środowiska

To środowisko szkoleniowe zawiera całe oprogramowanie, kod i dane niezbędne do pracy z tym kursem, więc nie musisz niczego instalować samodzielnie.

Codespace jest skonfigurowany z interfejsem VSCode, który obejmuje eksplorator systemu plików, edytor kodu i terminal powłoki.
Wszystkie instrukcje podane podczas kursu (np. „otwórz plik", „edytuj kod" lub „uruchom to polecenie") odnoszą się do tych trzech części interfejsu VSCode, chyba że zaznaczono inaczej.

Jeśli przechodzisz przez ten kurs samodzielnie, zapoznaj się z [podstawami środowiska](../../envsetup/01_setup.md), aby uzyskać więcej szczegółów.

### Wymagania wersji

To szkolenie jest zaprojektowane dla Nextflow w wersji 25.10.2 lub nowszej **z WŁĄCZONYM parserem składni v2**.
Jeśli używasz lokalnego lub niestandardowego środowiska, upewnij się, że korzystasz z odpowiednich ustawień zgodnie z dokumentacją [tutaj](../../info/nxf_versions.md).

## Przygotuj się do pracy

Gdy Twój codespace już działa, musisz zrobić dwie rzeczy przed przystąpieniem do szkolenia: ustawić katalog roboczy dla tego konkretnego kursu i przyjrzeć się dostarczonym materiałom.

### Ustaw katalog roboczy

Domyślnie codespace otwiera się z katalogiem roboczym ustawionym w głównym katalogu wszystkich kursów szkoleniowych, ale w tym kursie będziemy pracować w katalogu `nf4-science/genomics/`.

Przejdź teraz do tego katalogu, uruchamiając to polecenie w terminalu:

```bash
cd nf4-science/genomics/
```

Możesz ustawić VSCode tak, aby skupił się na tym katalogu, dzięki czemu tylko odpowiednie pliki będą widoczne w pasku bocznym eksploratora plików:

```bash
code .
```

!!! tip "Wskazówka"

    Jeśli z jakiegokolwiek powodu wyjdziesz z tego katalogu (np. Twój codespace przejdzie w stan uśpienia), zawsze możesz użyć pełnej ścieżki, aby do niego wrócić, zakładając, że pracujesz w środowisku szkoleniowym GitHub Codespaces:

    ```bash
    cd /workspaces/training/nf4-science/genomics
    ```

Teraz spójrzmy na zawartość.

### Poznaj dostarczone materiały

Możesz przeglądać zawartość tego katalogu, używając eksploratora plików po lewej stronie przestrzeni roboczej szkolenia.
Alternatywnie możesz użyć polecenia `tree`.

Przez cały kurs używamy wyjścia z `tree` do przedstawienia struktury i zawartości katalogów w czytelnej formie, czasami z niewielkimi modyfikacjami dla przejrzystości.

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
Używamy rozwijanych sekcji takich jak ta, aby zwięźle wyświetlić oczekiwane wyjście poleceń oraz zawartość katalogów i plików.

- **Plik `genomics.nf`** to skrypt workflow'u, który będziesz rozbudowywać w trakcie kursu.

- **Katalog `modules`** zawiera szkieletowe pliki modułów, które wypełnisz podczas kursu.

- **Plik `nextflow.config`** to plik konfiguracyjny, który ustawia minimalne właściwości środowiska.
  Na razie możesz go zignorować.

- **Katalog `data`** zawiera dane wejściowe i powiązane zasoby, opisane później w kursie.

- **Katalog `solutions`** zawiera ukończone pliki modułów oraz rozwiązanie części 2, które może służyć jako punkt wyjścia dla części 3.
  Są przeznaczone do użycia jako punkt odniesienia do sprawdzenia Twojej pracy i rozwiązywania problemów.

## Lista kontrolna gotowości

Myślisz, że jesteś gotowy/gotowa, aby zacząć?

- [ ] Rozumiem cel tego kursu i jego wymagania wstępne
- [ ] Moje środowisko jest uruchomione i działa
- [ ] Ustawiłem/ustawiłam odpowiednio mój katalog roboczy

Jeśli możesz zaznaczyć wszystkie pola, możesz zaczynać.

**Aby przejść do [Części 1: Przegląd metody i ręczne testowanie](./01_method.md), kliknij strzałkę w prawym dolnym rogu tej strony.**
