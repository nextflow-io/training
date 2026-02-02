# Orientacja

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Tłumaczenie wspomagane przez AI - [dowiedz się więcej i zasugeruj ulepszenia](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

Środowisko szkoleniowe zawiera wszystkie oprogramowanie, kod i dane niezbędne do pracy z tym kursem szkoleniowym, więc nie musisz niczego instalować samodzielnie.
Jednakże potrzebujesz (bezpłatnego) konta, aby się zalogować, i powinieneś/powinnaś poświęcić kilka minut na zapoznanie się z interfejsem.

Jeśli jeszcze tego nie zrobiłeś/zrobiłaś, skorzystaj z [tego linku](../../../envsetup/) przed dalszym przejściem.

## Dostarczone materiały

Przez cały kurs szkoleniowy będziemy pracować w katalogu `nf4-science/genomics/`, do którego musisz przejść po otwarciu przestrzeni roboczej szkolenia.
Ten katalog zawiera wszystkie pliki kodu, dane testowe i pliki pomocnicze, których będziesz potrzebować.

Możesz swobodnie eksplorować zawartość tego katalogu; najłatwiejszym sposobem jest użycie eksploratora plików po lewej stronie przestrzeni roboczej szkolenia w interfejsie VSCode.
Alternatywnie możesz użyć polecenia `tree`.
W trakcie kursu używamy wyjścia z `tree` do przedstawienia struktury i zawartości katalogów w czytelnej formie, czasami z niewielkimi modyfikacjami dla przejrzystości.

Tutaj generujemy spis treści do drugiego poziomu w dół:

```bash
tree . -L 2
```

Jeśli uruchomisz to wewnątrz `nf4-science/genomics`, powinieneś/powinnaś zobaczyć następujące wyjście:

```console title="Zawartość katalogu"

.
├── data
│   ├── bam
│   ├── ref
│   ├── sample_bams.txt
│   └── samplesheet.csv
├── genomics-1.nf
├── genomics-2.nf
├── genomics-3.nf
├── genomics-4.nf
├── nextflow.config
└── solutions
    ├── modules
    ├── nf-test.config
    └── tests

6 directories, 8 files

```

!!!note "Uwaga"

    Nie martw się, jeśli wydaje się to dużo; przejdziemy przez odpowiednie elementy na każdym etapie kursu.
    Ma to na celu jedynie dać Ci przegląd.

**Oto podsumowanie tego, co powinieneś wiedzieć, aby rozpocząć:**

- **Pliki `.nf`** to skrypty workflow'u, które są nazwane na podstawie tego, w jakiej części kursu są używane.

- **Plik `nextflow.config`** to plik konfiguracyjny, który ustawia minimalne właściwości środowiska.
  Na razie możesz go zignorować.

- **Katalog `data`** zawiera dane wejściowe i powiązane zasoby, opisane później w kursie.

- **Katalog `solutions`** zawiera pliki modułów i konfiguracje testów, które są wynikiem części 3 i 4 kursu.
  Są przeznaczone do użycia jako punkt odniesienia do sprawdzenia Twojej pracy i rozwiązywania problemów.

!!!tip "Wskazówka"

    Jeśli z jakiegokolwiek powodu wyjdziesz z tego katalogu, zawsze możesz uruchomić to polecenie, aby do niego wrócić:

    ```bash
    cd /workspaces/training/nf4-science/genomics
    ```

Teraz, aby rozpocząć kurs, kliknij strzałkę w prawym dolnym rogu tej strony.
