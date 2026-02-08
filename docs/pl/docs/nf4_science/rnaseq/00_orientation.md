# Orientacja

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Tłumaczenie wspomagane przez AI - [dowiedz się więcej i zasugeruj ulepszenia](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

Środowisko szkoleniowe zawiera wszystkie oprogramowanie, kod i dane niezbędne do pracy z tym kursem, więc nie musisz niczego instalować samodzielnie.
Potrzebujesz jednak (darmowego) konta, aby się zalogować, i powinieneś poświęcić kilka minut na zapoznanie się z interfejsem.

Jeśli jeszcze tego nie zrobiłeś, ukończ mini-kurs [Konfiguracja środowiska](../../envsetup/) przed dalszą pracą.

## Dostarczone materiały

W trakcie tego kursu będziemy pracować w katalogu `nf4-science/rnaseq/`, do którego należy przejść po otwarciu obszaru roboczego szkolenia.
Ten katalog zawiera wszystkie pliki kodu, dane testowe i pliki pomocnicze, których będziesz potrzebować.

Możesz swobodnie eksplorować zawartość tego katalogu; najłatwiejszym sposobem jest użycie eksploratora plików po lewej stronie obszaru roboczego szkolenia w interfejsie VSCode.
Alternatywnie możesz użyć polecenia `tree`.
W trakcie kursu używamy wyjścia `tree` do przedstawienia struktury katalogów i zawartości w czytelnej formie, czasami z drobnymi modyfikacjami dla przejrzystości.

Tutaj generujemy spis treści do drugiego poziomu:

```bash
tree . -L 3
```

??? success "Zawartość katalogu"

    ```console
    rnaseq
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

!!!note "Uwaga"

    Nie martw się, jeśli wydaje się to dużo; przejdziemy przez odpowiednie elementy na każdym etapie kursu.
    To ma jedynie na celu dać Ci przegląd.

**Oto podsumowanie tego, co powinieneś wiedzieć, aby zacząć:**

- **Plik `rnaseq.nf`** to zarys skryptu workflow'u, który będziemy rozwijać.

- **Plik `nextflow.config`** to plik konfiguracyjny, który ustawia minimalne właściwości środowiska. Na razie możesz go zignorować.

- **Katalog `data`** zawiera dane wejściowe i powiązane zasoby:

  - _Genom referencyjny_ o nazwie `genome.fa` składający się z małego regionu ludzkiego chromosomu 20 (z hg19/b37).
  - _Dane RNAseq_, które zostały ograniczone do małego regionu, aby zmniejszyć rozmiar plików, w katalogu `reads/`.
  - _Pliki CSV_ zawierające identyfikatory i ścieżki przykładowych plików danych, do przetwarzania wsadowego.

- **Katalog `solutions`** zawiera ukończone skrypty workflow'ów i moduły, które są wynikiem każdego etapu kursu.
  Mają one służyć jako punkt odniesienia do sprawdzenia Twojej pracy i rozwiązywania ewentualnych problemów.
  Numer w nazwie pliku odpowiada etapowi odpowiedniej części kursu.

!!!tip "Wskazówka"

    Jeśli z jakiegokolwiek powodu opuścisz ten katalog, zawsze możesz uruchomić to polecenie, aby do niego wrócić:

    ```bash
    cd /workspaces/training/nf4-science/rnaseq
    ```

Teraz, aby rozpocząć kurs, kliknij strzałkę w prawym dolnym rogu tej strony.
