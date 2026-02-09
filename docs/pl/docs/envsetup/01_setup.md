# GitHub Codespaces

GitHub Codespaces to platforma internetowa, która pozwala nam udostępnić wstępnie skonfigurowane środowisko szkoleniowe, oparte na maszynach wirtualnych w chmurze.
Platforma jest zarządzana przez GitHub (który należy do Microsoftu) i jest dostępna za darmo (z limitami użytkowania) dla każdego, kto posiada konto GitHub.

!!! warning "Ostrzeżenie"

    Konta powiązane z organizacjami mogą podlegać pewnym dodatkowym ograniczeniom.
    Jeśli to Twój przypadek, możesz potrzebować użyć niezależnego konta osobistego lub zamiast tego skorzystać z instalacji lokalnej.

## Tworzenie konta GitHub

Możesz utworzyć darmowe konto GitHub na [stronie głównej GitHub](https://github.com/).

## Uruchamianie GitHub Codespace

Po zalogowaniu się do GitHub, otwórz ten link w przeglądarce, aby otworzyć środowisko szkoleniowe Nextflow: <https://codespaces.new/nextflow-io/training?quickstart=1&ref=master>

Alternatywnie możesz kliknąć przycisk pokazany poniżej, który jest powtarzany w każdym kursie szkoleniowym (zazwyczaj na stronie Orientacji).

[![Open in GitHub Codespaces](https://github.com/codespaces/badge.svg)](https://codespaces.new/nextflow-io/training?quickstart=1&ref=master)

Powinna pojawić się strona, na której możesz utworzyć nowy GitHub Codespace:

![Create a GitHub Codespace](img/codespaces_create.png)

### Konfiguracja

Do ogólnego użytku nie powinieneś potrzebować niczego konfigurować.
O ile nie określono inaczej w kursie, który rozpoczynasz, możesz po prostu kliknąć główny przycisk, aby kontynuować.

Możliwe jest jednak dostosowanie środowiska poprzez kliknięcie przycisku "Change options".

??? info "Opcje konfiguracji"

    Jeśli klikniesz przycisk "Change options", otrzymasz możliwość dostosowania następujących elementów:

    #### Branch

    Pozwala to wybrać inną wersję materiałów szkoleniowych.
    Gałąź `master` zazwyczaj zawiera poprawki błędów oraz materiały, które zostały niedawno opracowane i zatwierdzone, ale nie zostały jeszcze opublikowane na stronie internetowej.
    Inne gałęzie zawierają prace w toku, które mogą nie być w pełni funkcjonalne.

    #### Machine type

    Pozwala to dostosować maszynę wirtualną, której będziesz używać do pracy ze szkoleniem.

    Użycie maszyny z większą liczbą rdzeni pozwala lepiej wykorzystać zdolność Nextflow'a do paralelizacji wykonywania workflow'ów.
    Jednak szybciej wyczerpie to Twój darmowy limit, dlatego nie zalecamy zmiany tego ustawienia, chyba że jest to zalecane w instrukcjach kursu, który planujesz odbyć.

    Zobacz "Limity GitHub Codespaces" poniżej, aby uzyskać więcej informacji o limitach.

### Czas uruchamiania

Otwieranie nowego środowiska GitHub Codespaces po raz pierwszy może zająć kilka minut, ponieważ system musi skonfigurować Twoją maszynę wirtualną, więc nie martw się, jeśli pojawi się czas oczekiwania.
Nie powinno to jednak zająć więcej niż pięć minut.

## Nawigacja po interfejsie szkoleniowym

Po załadowaniu GitHub Codespaces powinieneś zobaczyć coś podobnego do poniższego (które może otworzyć się w trybie jasnym w zależności od preferencji Twojego konta):

![GitHub Codespaces welcome](img/codespaces_welcome.png)

To jest interfejs IDE VSCode, popularnej aplikacji do tworzenia kodu, którą zalecamy używać do rozwoju w Nextflow.

- **Główny edytor** to miejsce, w którym będą otwierane kod Nextflow'a i inne pliki tekstowe. To tutaj będziesz edytować kod. Po otwarciu codespace wyświetli się podgląd pliku `README.md`.
- **Terminal** poniżej głównego edytora pozwala uruchamiać polecenia. To tutaj będziesz uruchamiać wszystkie polecenia podane w instrukcjach kursu.
- **Panel boczny** pozwala dostosować środowisko i wykonywać podstawowe zadania (kopiowanie, wklejanie, otwieranie plików, wyszukiwanie, git itp.). Domyślnie jest otwarty na eksploratorze plików, który pozwala przeglądać zawartość repozytorium. Kliknięcie pliku w eksploratorze otworzy go w głównym oknie edytora.

Możesz dostosować względne proporcje paneli okna według własnych preferencji.

<!-- TODO (future) Link to development best practices side quest? -->

## Inne uwagi dotyczące używania GitHub Codespaces

### Wznawianie sesji

Po utworzeniu środowiska możesz łatwo je wznowić lub zrestartować i kontynuować od miejsca, w którym skończyłeś.
Twoje środowisko wygaśnie po 30 minutach bezczynności i zapisze Twoje zmiany przez maksymalnie 2 tygodnie.

Możesz ponownie otworzyć środowisko z <https://github.com/codespaces/>.
Poprzednie środowiska będą wymienione.
Kliknij sesję, aby ją wznowić.

![List GitHub Codespace sessions](img/codespaces_list.png)

Jeśli zapisałeś URL swojego poprzedniego środowiska GitHub Codespaces, możesz po prostu otworzyć go w przeglądarce.
Alternatywnie kliknij ten sam przycisk, którego użyłeś do jego utworzenia:

[![Open in GitHub Codespaces](https://github.com/codespaces/badge.svg)](https://codespaces.new/nextflow-io/training?quickstart=1&ref=master)

Powinieneś zobaczyć poprzednią sesję, domyślną opcją jest jej wznowienie:

![Resume a GitHub Codespace](img/codespaces_resume.png)

### Zapisywanie plików na lokalnym komputerze

Aby zapisać dowolny plik z panelu eksploratora, kliknij plik prawym przyciskiem myszy i wybierz `Download`.

### Zarządzanie limitami GitHub Codespaces

GitHub Codespaces daje Ci do 15 GB-miesięcy przestrzeni dyskowej miesięcznie i 120 godzin rdzeniowych miesięcznie.
Jest to równoważne około 60 godzinom działania domyślnego środowiska przy użyciu standardowego obszaru roboczego (2 rdzenie, 8 GB RAM i 32 GB przestrzeni dyskowej).

Możesz je tworzyć z większymi zasobami (zobacz wyjaśnienie powyżej), ale szybciej wyczerpie to Twoje darmowe użytkowanie i będziesz mieć mniej godzin dostępu do tej przestrzeni.
Na przykład, jeśli wybierzesz maszynę 4-rdzeniową zamiast domyślnej 2-rdzeniowej, Twój limit wyczerpie się w połowie czasu.

Opcjonalnie możesz wykupić dostęp do większych zasobów.

Więcej informacji znajdziesz w dokumentacji GitHub:
[About billing for GitHub Codespaces](https://docs.github.com/en/billing/managing-billing-for-your-products/managing-billing-for-github-codespaces/about-billing-for-github-codespaces)
