# Ręczna instalacja

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Tłumaczenie wspomagane przez AI - [dowiedz się więcej i zasugeruj ulepszenia](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

Możesz ręcznie zainstalować wszystko, co jest potrzebne do uruchomienia szkoleń we własnym lokalnym środowisku.

Udokumentowaliśmy tutaj, jak to zrobić na standardowych systemach zgodnych z POSIX (zakładając komputer osobisty, taki jak laptop).
Pamiętaj, że niektóre szczegóły mogą się różnić w zależności od Twojego konkretnego systemu.

!!! tip "Wskazówka"

    Zanim przejdziesz dalej, zastanów się, czy nie lepiej byłoby skorzystać z [podejścia opartego na Devcontainers](03_devcontainer.md)?
    Zapewnia ono dostęp do wszystkich niezbędnych narzędzi i zależności bez konieczności ręcznej instalacji.

## Ogólne wymagania programowe

Nextflow'a można używać na dowolnym systemie zgodnym z POSIX (Linux, macOS, Windows Subsystem for Linux itp.) z zainstalowaną Javą.
Nasze kursy szkoleniowe mają kilka dodatkowych wymagań.

W sumie będziesz potrzebować następującego oprogramowania:

- Bash lub równoważna powłoka
- [Java 11 (lub nowsza, do 21)](https://www.oracle.com/technetwork/java/javase/downloads/index.html)
- [Git](https://git-scm.com/)
- [Docker](https://docs.docker.com/get-docker/)
- [Conda](https://conda.io/) 4.5 (lub nowsza)
- [VSCode](https://code.visualstudio.com) z [rozszerzeniem Nextflow](https://www.nextflow.io/docs/latest/developer-env.html#devenv-nextflow)

Aplikacja VSCode jest technicznie opcjonalna, ale zdecydowanie zalecamy jej użycie zarówno do pracy z kursami, jak i do ogólnego rozwijania projektów w Nextflow'ie.

Podręcznik dokumentacji Nextflow'a zawiera instrukcje instalacji tych zależności w sekcji [Environment setup](https://www.nextflow.io/docs/latest/developer-env.html).

## Nextflow i narzędzia nf-core

Będziesz musiał zainstalować sam Nextflow oraz narzędzia nf-core, jak opisano w artykułach podlinkowanych poniżej:

- [Instalacja Nextflow'a](https://www.nextflow.io/docs/latest/install.html)
- [Narzędzia nf-core](https://nf-co.re/docs/nf-core-tools/installation)

Zalecamy użycie opcji samoinstalacji dla Nextflow'a i opcji PyPI dla narzędzi nf-core.

!!! warning "Kompatybilność wersji"

    <!-- Any update to this content needs to be copied to the home page -->
    **Od stycznia 2026 roku wszystkie nasze kursy szkoleniowe Nextflow wymagają wersji Nextflow'a 25.10.2 lub nowszej, z włączoną ścisłą składnią v2, chyba że zaznaczono inaczej.**

    Więcej informacji o wymaganiach wersji i ścisłej składni v2 znajdziesz w przewodniku [Wersje Nextflow'a](../info/nxf_versions.md).

    Starsze wersje materiałów szkoleniowych odpowiadające wcześniejszej składni są dostępne poprzez selektor wersji w pasku menu tej strony.

## Materiały szkoleniowe

Najprostszym sposobem pobrania materiałów szkoleniowych jest sklonowanie całego repozytorium za pomocą tego polecenia:

```bash
git clone https://github.com/nextflow-io/training.git
```

Każdy kurs ma swój własny katalog.
Aby rozpocząć pracę z kursem, otwórz okno terminala (najlepiej z poziomu aplikacji VSCode) i przejdź za pomocą `cd` do odpowiedniego katalogu.

Następnie możesz postępować zgodnie z instrukcjami kursu dostarczonymi na stronie internetowej.
