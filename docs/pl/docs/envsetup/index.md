# Opcje środowiska

Naszym celem jest zapewnienie spójnego i dokładnie przetestowanego środowiska, które pozwoli Ci skupić się na nauce Nextflow'a bez konieczności poświęcania czasu i wysiłku na zarządzanie oprogramowaniem.
W tym celu opracowaliśmy skonteneryzowane środowisko zawierające wszystkie niezbędne oprogramowanie, pliki kodu i przykładowe dane potrzebne do pracy ze wszystkimi naszymi kursami.

To skonteneryzowane środowisko można uruchomić od razu w Github Codespaces lub lokalnie w VS Code z rozszerzeniem Devcontainers.

<div class="grid cards" markdown>

-   :material-cloud-outline:{ .lg .middle } __Github Codespaces__

    ---

    GitHub Codespaces to usługa internetowa, która pozwala nam zapewnić gotowe środowisko szkoleniowe ze wszystkimi narzędziami i danymi, oparte na maszynach wirtualnych w chmurze. Jest dostępne bezpłatnie dla każdego posiadacza konta Github.

    [Użyj Github Codespaces:material-arrow-right:](01_setup.md){ .md-button .md-button--primary .mt-1 }

-   :material-laptop:{ .lg .middle } __Lokalne Devcontainers__

    ---

    VS Code z Devcontainers zapewnia lokalnie uruchamiane skonteneryzowane środowisko programistyczne ze wszystkimi wstępnie skonfigurowanymi narzędziami szkoleniowymi. Oferuje to samo gotowe środowisko co Codespaces, ale działa w całości na Twoim lokalnym sprzęcie.

    [Użyj Devcontainers lokalnie :material-arrow-right:](03_devcontainer.md){ .md-button .md-button--primary .mt-1 }

</div>

## Instrukcje instalacji ręcznej

Jeśli żadna z powyższych opcji nie odpowiada Twoim potrzebom, możesz odtworzyć to środowisko na własnym systemie lokalnym, instalując zależności oprogramowania ręcznie i klonując repozytorium szkoleniowe.

[Instalacja ręczna :material-arrow-right:](02_local.md){ .md-button .md-button--primary .mt-1 }

---

!!! info "Wycofanie Gitpod"

    Szkolenia Nextflow do lutego 2025 roku korzystały z [Gitpod](https://gitpod.io).
    Jednak twórcy Gitpod zdecydowali się wycofać darmową funkcjonalność na rzecz systemu [Gitpod Flex](https://www.gitpod.io/blog/introducing-gitpod-flex).
    Z tego powodu przeszliśmy na GitHub Codespaces, które również oferują środowisko programistyczne uruchamiane jednym kliknięciem, bez wcześniejszej konfiguracji.

    W zależności od tego, kiedy zarejestrowałeś się w Gitpod i kiedy dokładnie wycofają usługę, możesz nadal być w stanie uruchomić szkolenie w ich starym cloudowym IDE, choć nie możemy zagwarantować niezawodnego dostępu w przyszłości:
    [Otwórz w Gitpod](https://gitpod.io/#https://github.com/nextflow-io/training).
