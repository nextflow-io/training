# Orientacja

Ta orientacja zakłada, że otworzyłeś już środowisko szkoleniowe, klikając przycisk "Open in GitHub Codespaces".
Jeśli jeszcze tego nie zrobiłeś, zrób to teraz, najlepiej w drugim oknie lub karcie przeglądarki, aby móc wrócić do tych instrukcji.

[![Open in GitHub Codespaces](https://github.com/codespaces/badge.svg)](https://github.com/codespaces/new?hide_repo_select=true&ref=master&repo=290790519&skip_quickstart=true&machine=premiumLinux&devcontainer_path=.devcontainer%2Fdevcontainer.json)

!!!warning "Wymagany rozmiar maszyny"

    Upewnij się, że wybierasz **maszynę 8-rdzeniową** podczas tworzenia Codespace dla tego szkolenia. Workflow'y bioimagingu wymagają dodatkowych zasobów obliczeniowych.

## GitHub Codespaces

Środowisko GitHub Codespaces zawiera wszystkie oprogramowanie, kod i dane niezbędne do pracy z tym szkoleniem, więc nie musisz niczego instalować samodzielnie.
Potrzebujesz jednak (darmowego) konta GitHub, aby się zalogować, a jeśli nie znasz interfejsu, poświęć kilka minut na zapoznanie się z nim, wykonując mini-kurs [Orientacja w GitHub Codespaces](../../envsetup/index.md).

## Wstępne pobieranie obrazów Docker

Po otwarciu Codespace pobierzmy wszystkie obrazy Docker, które będą nam potrzebne podczas tego szkolenia.
Zaoszczędzi to czas później i zapewni płynne wykonywanie workflow'ów.

Otwórz nową kartę terminala i uruchom następujące polecenie:

```bash
nextflow run nf-core/molkart -profile docker,test -stub -resume --outdir results
```

To polecenie pobierze wszystkie niezbędne obrazy Docker w tle.
Możesz kontynuować resztę orientacji, podczas gdy to się wykonuje.

!!!tip

    Flaga `-stub` pozwala pipeline'owi działać szybko bez przetwarzania prawdziwych danych, co jest idealne do pobierania obrazów. Możesz monitorować postęp w karcie terminala.

## Katalog roboczy

Podczas tego szkolenia będziemy pracować w katalogu `nf4-science/imaging/`.

Zmień teraz katalog, uruchamiając to polecenie w terminalu:

```bash
cd nf4-science/imaging/
```

!!!tip

    Jeśli z jakiegokolwiek powodu wyjdziesz z tego katalogu, zawsze możesz użyć pełnej ścieżki, aby do niego wrócić, zakładając, że pracujesz w środowisku szkoleniowym GitHub Codespaces:

    ```bash
    cd /workspaces/training/nf4-science/imaging
    ```

**Teraz, aby rozpocząć kurs, kliknij strzałkę w prawym dolnym rogu tej strony.**
