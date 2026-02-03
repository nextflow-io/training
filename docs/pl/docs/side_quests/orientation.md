# Orientacja

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Tłumaczenie wspomagane przez AI - [dowiedz się więcej i zasugeruj ulepszenia](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

Środowisko GitHub Codespaces zawiera wszystkie oprogramowanie, kod i dane niezbędne do przepracowania tego kursu szkoleniowego, więc nie musisz niczego instalować samodzielnie.
Jednak potrzebujesz (bezpłatnego) konta, aby się zalogować, i powinieneś poświęcić kilka minut na zapoznanie się z interfejsem.

Jeśli jeszcze tego nie zrobiłeś, proszę kliknij [ten link](../../envsetup/) przed kontynuowaniem.

## Dostarczone materiały

Podczas tego kursu szkoleniowego będziemy pracować w katalogu `side-quests/`.
Ten katalog zawiera wszystkie pliki kodu, dane testowe i pliki pomocnicze, które będą potrzebne.

Zachęcamy do zapoznania się z zawartością tego katalogu; najłatwiejszym sposobem jest użycie eksploratora plików po lewej stronie obszaru roboczego GitHub Codespaces.
Alternatywnie możesz użyć polecenia `tree`.
Podczas kursu używamy wyjścia `tree` do przedstawienia struktury i zawartości katalogów w czytelnej formie, czasami z niewielkimi modyfikacjami dla przejrzystości.

Tutaj generujemy spis treści do drugiego poziomu w dół:

```bash
tree . -L 2
```

Jeśli uruchomisz to wewnątrz `side-quests`, powinieneś zobaczyć następujące wyjście:

```console title="Zawartość katalogu"
.
├── metadata
├── nf-core
├── nf-test
├── solutions
├── splitting_and_grouping
└── workflows_of_workflows
```

**Oto podsumowanie tego, co powinieneś wiedzieć na początek:**

- **Każdy katalog odpowiada indywidualnemu zadaniu pobocznemu.**
  Ich zawartość jest szczegółowo opisana na stronie odpowiedniego zadania pobocznego.

- **Katalog `solutions`** zawiera ukończone skrypty workflow i/lub modułów, które powstają w wyniku wykonania różnych kroków każdego zadania pobocznego.
  Mają one służyć jako punkt odniesienia do sprawdzenia Twojej pracy i rozwiązywania ewentualnych problemów.

!!!tip "Wskazówka"

    Jeśli z jakiegokolwiek powodu opuścisz ten katalog, zawsze możesz uruchomić to polecenie, aby do niego wrócić:

    ```bash
    cd /workspaces/training/side-quests
    ```

Teraz, aby rozpocząć kurs, kliknij strzałkę w prawym dolnym rogu tej strony.
