# Pierwsze kroki

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Tłumaczenie wspomagane przez AI - [dowiedz się więcej i zasugeruj ulepszenia](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

## Uruchom środowisko szkoleniowe

Aby skorzystać z gotowego środowiska udostępnianego przez nas na GitHub Codespaces, kliknij przycisk „Open in GitHub Codespaces" poniżej. Inne opcje znajdziesz w sekcji [Opcje środowiska](../../envsetup/index.md).

Zalecamy otwarcie środowiska szkoleniowego w nowej karcie lub oknie przeglądarki (użyj kliknięcia prawym przyciskiem myszy, ctrl+klik lub cmd+klik w zależności od Twojego sprzętu), abyś mógł czytać instrukcje podczas ładowania środowiska.
Musisz mieć te instrukcje otwarte równolegle, aby przejść przez kurs.

[![Open in GitHub Codespaces](https://github.com/codespaces/badge.svg)](https://codespaces.new/nextflow-io/training?quickstart=1&ref=master)

### Podstawy środowiska

To środowisko szkoleniowe zawiera wszystkie oprogramowanie, kod i dane niezbędne do pracy z kursem, więc nie musisz niczego instalować samodzielnie.

Codespace jest skonfigurowany z interfejsem VSCode, który obejmuje eksplorator systemu plików, edytor kodu i powłokę terminala.
Wszystkie instrukcje podawane w trakcie kursu (np. „otwórz plik", „edytuj kod" lub „uruchom to polecenie") odnoszą się do tych trzech części interfejsu VSCode, chyba że zaznaczono inaczej.

Jeśli pracujesz nad tym kursem samodzielnie, zapoznaj się z [podstawami środowiska](../../envsetup/01_setup.md), aby uzyskać więcej szczegółów.

### Wymagania dotyczące wersji

To szkolenie jest zaprojektowane dla Nextflow'a 25.10.2 lub nowszego **z WŁĄCZONYM parserem składni v2**.
Jeśli używasz lokalnego lub niestandardowego środowiska, upewnij się, że używasz prawidłowych ustawień zgodnie z dokumentacją [tutaj](../../info/nxf_versions.md).

## Przygotuj się do pracy

Gdy Twój codespace już działa, musisz wykonać dwie czynności przed rozpoczęciem szkolenia: ustawić katalog roboczy dla tego konkretnego kursu i przyjrzeć się dostarczonym materiałom.

### Ustaw katalog roboczy

Domyślnie codespace otwiera się z katalogiem roboczym ustawionym w katalogu głównym wszystkich kursów szkoleniowych, ale w tym kursie będziemy pracować w katalogu `nf4-science/rnaseq/`.

Zmień katalog teraz, uruchamiając to polecenie w terminalu:

```bash
cd nf4-science/rnaseq/
```

Możesz ustawić VSCode tak, aby skupił się na tym katalogu, dzięki czemu tylko odpowiednie pliki będą widoczne w pasku bocznym eksploratora plików:

```bash
code .
```

!!! tip "Wskazówka"

    Jeśli z jakiegokolwiek powodu opuścisz ten katalog (np. Twój codespace przejdzie w stan uśpienia), zawsze możesz użyć pełnej ścieżki, aby do niego wrócić, zakładając, że pracujesz w środowisku szkoleniowym Github Codespaces:

    ```bash
    cd /workspaces/training/nf4-science/rnaseq
    ```

Teraz przyjrzyjmy się zawartości.

### Poznaj dostarczone materiały

Możesz eksplorować zawartość tego katalogu, używając eksploratora plików po lewej stronie obszaru roboczego szkolenia.
Alternatywnie możesz użyć polecenia `tree`.

W trakcie kursu używamy wyjścia `tree` do przedstawienia struktury katalogów i zawartości w czytelnej formie, czasami z drobnymi modyfikacjami dla przejrzystości.

Tutaj generujemy spis treści do trzeciego poziomu:

```bash
tree . -L 3
```

??? abstract "Zawartość katalogu"

    ```console
    .
    ├── data
    │   ├── genome.fa
    │   ├── paired-end.csv
    │   ├── reads
    │   │   ├── ENCSR000COQ1_1.fastq.gz
    │   │   ├── ENCSR000COQ1_2.fastq.gz
    │   │   ├── ENCSR000COQ2_1.fastq.gz
    │   │   ├── ENCSR000COQ2_2.fastq.gz
    │   │   ├── ENCSR000COR1_1.fastq.gz
    │   │   ├── ENCSR000COR1_2.fastq.gz
    │   │   ├── ENCSR000COR2_1.fastq.gz
    │   │   ├── ENCSR000COR2_2.fastq.gz
    │   │   ├── ENCSR000CPO1_1.fastq.gz
    │   │   ├── ENCSR000CPO1_2.fastq.gz
    │   │   ├── ENCSR000CPO2_1.fastq.gz
    │   │   └── ENCSR000CPO2_2.fastq.gz
    │   └── single-end.csv
    ├── nextflow.config
    ├── rnaseq.nf
    └── solutions
        ├── modules
        │   ├── fastqc.nf
        │   ├── fastqc_pe.nf
        │   ├── hisat2_align.nf
        │   ├── hisat2_align_pe.nf
        │   ├── multiqc.nf
        │   ├── trim_galore.nf
        │   └── trim_galore_pe.nf
        ├── rnaseq-2.1.nf
        ├── rnaseq-2.2.nf
        ├── rnaseq-2.3.nf
        ├── rnaseq-3.1.nf
        ├── rnaseq-3.2.nf
        └── rnaseq_pe-3.3.nf
    ```

Kliknij na kolorowe pole, aby rozwinąć sekcję i zobaczyć jej zawartość.
Używamy takich zwijanych sekcji, aby wyświetlać oczekiwane wyjście poleceń, a także zawartość katalogów i plików w zwięzły sposób.

- **Plik `rnaseq.nf`** to zarys skryptu workflow'u, który będziesz rozwijać w trakcie kursu.

- **Katalog `modules`** zawiera zarysy modułów procesów, które wypełnisz podczas kursu.

- **Plik `nextflow.config`** to plik konfiguracyjny, który ustawia minimalne właściwości środowiska.
  Na razie możesz go zignorować.

- **Katalog `data`** zawiera dane wejściowe i powiązane zasoby, opisane później w kursie.

- **Katalog `solutions`** zawiera ukończone skrypty workflow'ów i moduły, które są wynikiem każdego etapu kursu.
  Mają one służyć jako punkt odniesienia do sprawdzenia Twojej pracy i rozwiązywania ewentualnych problemów.
  Rozwiązanie z Części 2 może być użyte jako punkt wyjścia dla Części 3.

## Lista kontrolna gotowości

Myślisz, że jesteś gotowy, aby zacząć?

- [ ] Rozumiem cel tego kursu i jego wymagania wstępne
- [ ] Moje środowisko działa
- [ ] Ustawiłem odpowiednio mój katalog roboczy

Jeśli możesz zaznaczyć wszystkie pola, możesz zaczynać.

**Aby przejść do [Części 1: Przegląd metody](./01_method.md), kliknij strzałkę w prawym dolnym rogu tej strony.**
