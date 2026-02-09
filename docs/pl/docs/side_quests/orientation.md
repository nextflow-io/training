# Orientacja

Środowisko GitHub Codespaces zawiera wszystkie niezbędne oprogramowanie, kod i dane potrzebne do przejścia przez ten kurs szkoleniowy, więc nie musisz niczego instalować samodzielnie.
Potrzebujesz jednak (darmowego) konta, aby się zalogować, i powinieneś poświęcić kilka minut na zapoznanie się z interfejsem.

Jeśli jeszcze tego nie zrobiłeś, przejdź pod [ten link](../../envsetup/) przed dalszą lekturą.

## Dostarczone materiały

W trakcie tego kursu szkoleniowego będziemy pracować w katalogu `side-quests/`.
Katalog ten zawiera wszystkie pliki kodu, dane testowe i pliki pomocnicze, których będziesz potrzebować.

Zachęcamy do zapoznania się z zawartością tego katalogu; najłatwiejszym sposobem jest użycie eksploratora plików po lewej stronie obszaru roboczego GitHub Codespaces.
Alternatywnie możesz użyć polecenia `tree`.
W trakcie kursu używamy wyjścia polecenia `tree` do przedstawienia struktury i zawartości katalogów w czytelnej formie, czasami z niewielkimi modyfikacjami dla przejrzystości.

Tutaj generujemy spis treści do drugiego poziomu w dół:

```bash
tree . -L 2
```

Jeśli uruchomisz to polecenie wewnątrz `side-quests`, powinieneś zobaczyć następujące wyjście:

```console title="Directory contents"
.
├── metadata
├── nf-core
├── nf-test
├── solutions
├── splitting_and_grouping
└── workflows_of_workflows
```

**Oto podsumowanie tego, co powinieneś wiedzieć na początek:**

- **Każdy katalog odpowiada pojedynczemu zadaniu dodatkowemu.**
  Ich zawartość jest szczegółowo opisana na stronie odpowiadającego zadania dodatkowego.

- **Katalog `solutions`** zawiera ukończone skrypty workflow'ów i/lub modułów, które powstają w wyniku przejścia przez różne etapy każdego zadania dodatkowego.
  Są one przeznaczone do użycia jako punkt odniesienia do sprawdzenia Twojej pracy i rozwiązywania ewentualnych problemów.

!!!tip "Wskazówka"

    Jeśli z jakiegokolwiek powodu wyjdziesz z tego katalogu, zawsze możesz uruchomić to polecenie, aby do niego wrócić:

    ```bash
    cd /workspaces/training/side-quests
    ```

Teraz, aby rozpocząć kurs, kliknij strzałkę w prawym dolnym rogu tej strony.
