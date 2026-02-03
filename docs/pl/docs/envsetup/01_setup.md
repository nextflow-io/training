# GitHub Codespaces

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Tłumaczenie wspomagane przez AI - [dowiedz się więcej i zasugeruj ulepszenia](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

GitHub Codespaces to platforma internetowa, która pozwala nam dostarczyć prekonfigurowane środowisko do szkoleń, wspierane przez maszyny wirtualne w chmurze.
Platforma jest obsługiwana przez GitHub (który należy do Microsoft) i jest dostępna bezpłatnie (z limitami użycia) dla każdego z kontem GitHub.

!!! warning "Ostrzeżenie"

    Konta powiązane z organizacjami mogą podlegać pewnym dodatkowym ograniczeniom.
    W takim przypadku może być konieczne użycie niezależnego konta osobistego lub lokalna instalacja.

## Tworzenie konta GitHub

Możesz utworzyć darmowe konto GitHub na [stronie głównej GitHub](https://github.com/).

## Uruchamianie GitHub Codespace

Po zalogowaniu się do GitHub otwórz ten link w przeglądarce, aby otworzyć środowisko szkoleniowe Nextflow: <https://codespaces.new/nextflow-io/training?quickstart=1&ref=master>

Alternatywnie możesz kliknąć przycisk pokazany poniżej, który jest powtórzony w każdym kursie szkoleniowym (zazwyczaj na stronie Orientacji).

[![Open in GitHub Codespaces](https://github.com/codespaces/badge.svg)](https://codespaces.new/nextflow-io/training?quickstart=1&ref=master)

Powinieneś zobaczyć stronę, na której możesz utworzyć nowy GitHub Codespace:

![Create a GitHub Codespace](img/codespaces_create.png)

### Konfiguracja

Do ogólnego użytku nie powinieneś potrzebować niczego konfigurować.
O ile nie zaznaczono inaczej w kursie, który rozpoczynasz, możesz po prostu kliknąć główny przycisk, aby kontynuować.

Jednak możliwe jest dostosowanie środowiska, klikając przycisk "Change options".

??? info "Opcje konfiguracji"

    Jeśli klikniesz przycisk "Change options", będziesz mieć możliwość dostosowania następujących elementów:

    #### Branch

    Pozwala wybrać inną wersję materiałów szkoleniowych.
    Branch `master` zazwyczaj zawiera poprawki błędów i materiały, które zostały niedawno opracowane i zatwierdzone, ale nie zostały jeszcze opublikowane na stronie.
    Inne branche zawierają prace w toku, które mogą nie być w pełni funkcjonalne.

    #### Machine type

    Pozwala dostosować maszynę wirtualną, której będziesz używać do pracy ze szkoleniem.

    Używanie maszyny z większą liczbą rdzeni pozwala lepiej wykorzystać zdolność Nextflow do równoległego wykonywania workflow'ów. Jednak zużyje to szybciej Twój darmowy limit. Nie zalecamy zmiany tego ustawienia, chyba że wymaga tego kurs, który zamierzasz realizować.

    Zobacz 'Limity GitHub Codespaces' poniżej, aby uzyskać więcej szczegółów o limitach.

### Czas uruchamiania

Otwieranie nowego środowiska GitHub Codespaces po raz pierwszy może zająć kilka minut, ponieważ system musi skonfigurować Twoją maszynę wirtualną, więc nie martw się, jeśli jest czas oczekiwania.
Jednak nie powinno to trwać dłużej niż pięć minut.

## Nawigacja w interfejsie szkoleniowym

Po załadowaniu GitHub Codespaces powinieneś zobaczyć coś podobnego do poniższego (może otworzyć się w trybie jasnym w zależności od preferencji Twojego konta):

![GitHub Codespaces welcome](img/codespaces_welcome.png)

To jest interfejs VSCode IDE, popularnej aplikacji do tworzenia kodu, którą polecamy do pracy z Nextflow.

- **Główny edytor** to miejsce, gdzie będą otwierane kod Nextflow i inne pliki tekstowe. Tutaj będziesz edytować kod. Po otwarciu codespace wyświetli podgląd pliku `README.md`.
- **Terminal** poniżej głównego edytora pozwala uruchamiać polecenia. Tutaj będziesz uruchamiać wszystkie polecenia podane w instrukcjach kursu.
- **Pasek boczny** pozwala dostosować środowisko i wykonywać podstawowe zadania (kopiowanie, wklejanie, otwieranie plików, wyszukiwanie, git itp.). Domyślnie jest otwarty na eksploratorze plików, który pozwala przeglądać zawartość repozytorium. Kliknięcie pliku w eksploratorze otworzy go w głównym oknie edytora.

Możesz dostosować względne proporcje paneli okna według własnych upodobań.

<!-- TODO (future) Link to development best practices side quest? -->

## Inne uwagi dotyczące korzystania z GitHub Codespaces

### Wznawianie sesji

Po utworzeniu codespace możesz je łatwo wznowić lub uruchomić ponownie i kontynuować od miejsca, w którym skończyłeś.
Sesja automatycznie zakończy się po 30 minutach bezczynności i zachowa Twoje zmiany przez maksymalnie 2 tygodnie.

Możesz ponownie otworzyć środowisko z <https://github.com/codespaces/>.
Poprzednie środowiska będą wyświetlone na liście.
Kliknij sesję, aby ją wznowić.

![List GitHub Codespace sessions](img/codespaces_list.png)

Jeśli zapisałeś URL Swojego poprzedniego środowiska GitHub Codespaces, możesz po prostu otworzyć go w przeglądarce.
Alternatywnie kliknij ten sam przycisk, którego użyłeś do jego utworzenia:

[![Open in GitHub Codespaces](https://github.com/codespaces/badge.svg)](https://codespaces.new/nextflow-io/training?quickstart=1&ref=master)

Powinieneś zobaczyć poprzednią sesję, domyślną opcją jest jej wznowienie:

![Resume a GitHub Codespace](img/codespaces_resume.png)

### Zapisywanie plików na lokalną maszynę

Aby zapisać dowolny plik z panelu eksploratora, kliknij prawym przyciskiem myszy na pliku i wybierz `Download`.

### Zarządzanie limitami GitHub Codespaces

GitHub Codespaces daje Ci do 15 GB-miesięcy przestrzeni dyskowej miesięcznie i 120 godzin rdzeniowych miesięcznie.
Jest to równoważne około 60 godzinom domyślnego czasu działania środowiska przy użyciu standardowej przestrzeni roboczej (2 rdzenie, 8 GB RAM i 32 GB przestrzeni dyskowej).

Możesz tworzyć je z większymi zasobami (patrz wyjaśnienie powyżej), ale spowoduje to szybsze zużycie darmowego limitu i będziesz mieć mniej godzin dostępu do tej przestrzeni.
Na przykład, jeśli wybierzesz maszynę 4-rdzeniową zamiast domyślnej 2-rdzeniowej, Twój limit wyczerpie się o połowę szybciej.

Opcjonalnie możesz zakupić dostęp do większych zasobów.

Więcej informacji znajdziesz w dokumentacji GitHub:
[About billing for GitHub Codespaces](https://docs.github.com/en/billing/managing-billing-for-your-products/managing-billing-for-github-codespaces/about-billing-for-github-codespaces)
