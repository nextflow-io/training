# Instalacja ręczna

Możliwe jest ręczne zainstalowanie wszystkiego, czego potrzebujesz do uruchomienia szkolenia we własnym lokalnym środowisku.

Poniżej udokumentowaliśmy, jak to zrobić na standardowych systemach zgodnych z POSIX (zakładając komputer osobisty, taki jak laptop).
Pamiętaj, że niektóre szczegóły mogą się różnić w zależności od Twojego konkretnego systemu.

!!! tip "Wskazówka"

    Zanim przejdziesz dalej, czy rozważyłeś/rozważyłaś użycie [podejścia Devcontainers](03_devcontainer.md)?
    Zapewnia ono wszystkie niezbędne narzędzia i zależności bez konieczności ręcznej instalacji.

## Ogólne wymagania dotyczące oprogramowania

Nextflow'a można używać na dowolnym systemie zgodnym z POSIX (Linux, macOS, Windows Subsystem for Linux itp.) z zainstalowaną Javą.
Nasze kursy szkoleniowe mają kilka dodatkowych wymagań.

W sumie będziesz potrzebować następującego zainstalowanego oprogramowania:

- Bash lub równoważna powłoka
- [Java 11 (lub nowsza, do wersji 21)](https://www.oracle.com/technetwork/java/javase/downloads/index.html)
- [Git](https://git-scm.com/)
- [Docker](https://docs.docker.com/get-docker/)
- [Conda](https://conda.io/) 4.5 (lub nowsza)
- [VSCode](https://code.visualstudio.com) z [rozszerzeniem Nextflow](https://www.nextflow.io/docs/latest/developer-env.html#devenv-nextflow)

Aplikacja VSCode jest technicznie opcjonalna, ale zdecydowanie zalecamy jej użycie zarówno do pracy z kursami, jak i do ogólnej pracy deweloperskiej z Nextflow'em.

Dokumentacja Nextflow'a zawiera instrukcje instalacji tych zależności w sekcji [Environment setup](https://www.nextflow.io/docs/latest/developer-env.html).

## Nextflow i narzędzia nf-core

Będziesz potrzebować zainstalować sam Nextflow, a także narzędzia nf-core, zgodnie ze szczegółami w artykułach linkowanych poniżej:

- [Instalacja Nextflow'a](https://www.nextflow.io/docs/latest/install.html)
- [Narzędzia nf-core](https://nf-co.re/docs/nf-core-tools/installation)

Zalecamy użycie opcji samoinstalacji dla Nextflow'a i opcji PyPI dla narzędzi nf-core.

!!! warning "Ostrzeżenie"

    <!-- Any update to this content needs to be copied to the home page -->
    **Od stycznia 2026 roku wszystkie nasze kursy szkoleniowe Nextflow wymagają wersji Nextflow 25.10.2 lub nowszej, z aktywowaną ścisłą składnią v2, chyba że zaznaczono inaczej.**

    Więcej informacji o wymaganiach dotyczących wersji i ścisłej składni v2 znajdziesz w przewodniku [Wersje Nextflow'a](../info/nxf_versions.md).

    Starsze wersje materiałów szkoleniowych odpowiadające wcześniejszej składni są dostępne za pomocą selektora wersji w pasku menu tej strony internetowej.

## Materiały szkoleniowe

Najłatwiejszym sposobem pobrania materiałów szkoleniowych jest sklonowanie całego repozytorium za pomocą tego polecenia:

```bash
git clone https://github.com/nextflow-io/training.git
```

Każdy kurs ma swój własny katalog.
Aby przejść przez kurs, otwórz okno terminala (najlepiej z poziomu aplikacji VSCode) i przejdź do odpowiedniego katalogu poleceniem `cd`.

Następnie możesz postępować zgodnie z instrukcjami kursu dostępnymi na stronie internetowej.
