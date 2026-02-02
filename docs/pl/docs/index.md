---
title: Strona główna
description: Witamy w portalu szkoleniowym społeczności Nextflow!
hide:
  - toc
  - footer
---

# Szkolenie Nextflow

<div class="grid cards" markdown>

-   :material-book-open-variant:{ .lg .middle } __Kursy do samodzielnej nauki__

    ---

    **Witamy w portalu szkoleniowym społeczności Nextflow!**

    Poniższe kursy szkoleniowe są zaprojektowane jako materiały do samodzielnej nauki.
    Możesz przerabiać je we własnym tempie w dowolnym czasie, korzystając ze środowiska internetowego, które udostępniamy przez GitHub Codespaces, lub we własnym środowisku.

    [Przeglądaj kursy :material-arrow-right:](#katalog-kursow-szkoleniowych-nextflow){ .md-button .md-button--primary .mt-1 }

-   :material-information-outline:{ .lg .middle } __Dodatkowe informacje__

    ---

    ??? warning "Kompatybilność wersji"

        <!-- Any update to this content needs to be copied to the local installation page -->
        **Od stycznia 2026 roku wszystkie nasze kursy szkoleniowe Nextflow wymagają wersji Nextflow 25.10.2 lub nowszej, z włączoną ścisłą składnią, chyba że zaznaczono inaczej.**

        Więcej informacji o wymaganiach wersji i ścisłej składni znajdziesz w [przewodniku migracji dokumentacji Nextflow](https://nextflow.io/docs/latest/strict-syntax.html).

        Starsze wersje materiałów szkoleniowych odpowiadające wcześniejszej składni są dostępne poprzez selektor wersji w pasku menu tej strony.

    ??? terminal "Opcje środowiska"

        Udostępniamy internetowe środowisko szkoleniowe, w którym wszystko, czego potrzebujesz do szkolenia, jest preinstalowane, dostępne przez GitHub Codespaces (wymaga darmowego konta GitHub).

        [![Open in GitHub Codespaces](https://github.com/codespaces/badge.svg)](https://codespaces.new/nextflow-io/training?quickstart=1&ref=master)

        Jeśli to nie odpowiada Twoim potrzebom, zapoznaj się z innymi [opcjami środowiska](./envsetup/index.md).

    ??? learning "Wydarzenia szkoleniowe"

        Jeśli wolisz uczestniczyć w szkoleniu Nextflow w ramach zorganizowanego wydarzenia, istnieje wiele możliwości. Polecamy sprawdzić następujące opcje:

        - **[Tygodnie szkoleniowe]()** organizowane kwartalnie przez zespół społeczności
        - **[Wydarzenia Seqera](https://seqera.io/events/)** obejmują wydarzenia szkoleniowe organizowane przez Seqera (szukaj 'Seqera Sessions' i 'Nextflow Summit')
        - **[Ambasadorzy Nextflow]()** organizują wydarzenia dla swojej lokalnej społeczności
        - **[Wydarzenia nf-core](https://nf-co.re/events)** obejmują hackathony społecznościowe

    ??? people "Informacje dla prowadzących szkolenia"

        Jeśli jesteś instruktorem prowadzącym własne szkolenia, możesz korzystać z naszych materiałów bezpośrednio z portalu szkoleniowego, pod warunkiem podania odpowiedniego źródła. Szczegóły znajdziesz w sekcji 'Licencja i polityka współpracy' poniżej.

        Ponadto chętnie dowiemy się, jak moglibyśmy lepiej wspierać Twoje działania szkoleniowe! Skontaktuj się z nami pod adresem [community@seqera.io](mailto:community@seqera.io) lub na forum społecznościowym (zobacz stronę [Pomoc](help.md)).

    ??? licensing "Licencja open-source i polityka współpracy"

        [![Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International (CC BY-NC-SA 4.0)](assets/img/cc_by-nc-sa.svg){ align=right }](https://creativecommons.org/licenses/by-nc-sa/4.0/)

        Te materiały szkoleniowe są tworzone i utrzymywane przez [Seqera](https://seqera.io) i udostępniane na licencji open-source ([CC BY-NC-SA](https://creativecommons.org/licenses/by-nc-sa/4.0/)) dla dobra społeczności. Jeśli chcesz wykorzystać te materiały w sposób wykraczający poza zakres licencji (zwróć uwagę na ograniczenia dotyczące użytku komercyjnego i redystrybucji),         skontaktuj się z nami pod adresem [community@seqera.io](mailto:community@seqera.io), aby omówić Swoją prośbę.

        Z przyjemnością przyjmujemy ulepszenia, poprawki i zgłoszenia błędów od społeczności. Każda strona ma ikonę :material-file-edit-outline: w prawym górnym rogu, która prowadzi do repozytorium kodu, gdzie możesz zgłaszać problemy lub proponować zmiany w materiałach szkoleniowych poprzez pull request. Więcej szczegółów znajdziesz w pliku `README.md` w repozytorium.

</div>

!!! note "Tłumaczenie wspomagane przez AI"

    To tłumaczenie zostało utworzone przy użyciu sztucznej inteligencji i zweryfikowane przez ludzkich tłumaczy.
    Zachęcamy do przekazywania opinii i sugestii ulepszeń.
    Więcej informacji znajdziesz w naszym [przewodniku tłumaczenia](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md).

## Katalog kursów szkoleniowych Nextflow

<div class="grid cards" markdown>

-   :material-walk:{ .lg .middle } __Ścieżka wprowadzająca__

    ---

    ### :material-compass:{.nextflow-primary} Nextflow dla początkujących {.mt-1}

    Kursy niezależne od dziedziny, przeznaczone dla osób, które są całkowicie nowe w Nextflow. Każdy z nich składa się z serii modułów szkoleniowych. Te moduły pomagają uczącym się stopniowo rozwijać swoje umiejętności.

    ??? courses "**Hello Nextflow:** Naucz się tworzyć własne pipeline'y"

        Ten kurs obejmuje podstawowe komponenty języka Nextflow w wystarczającym szczegółach, aby umożliwić tworzenie prostych, ale w pełni funkcjonalnych pipeline'ów, a także kluczowe elementy projektowania, tworzenia i konfigurowania pipeline'ów.

        [Rozpocznij szkolenie Hello Nextflow :material-arrow-right:](hello_nextflow/index.md){ .md-button .md-button--secondary }

    ??? courses "**Nextflow Run:** Naucz się uruchamiać istniejące pipeline'y"

        Zwięzłe wprowadzenie do uruchamiania i konfigurowania pipeline'ów Nextflow, oparte na kursie Hello Nextflow dla programistów, ale z mniejszym naciskiem na kod. Obejmuje wykonywanie, wyjścia, podstawową strukturę kodu i konfigurację dla różnych środowisk obliczeniowych.

        [Rozpocznij szkolenie Nextflow Run :material-arrow-right:](nextflow_run/index.md){ .md-button .md-button--secondary }

    ---

    ### :material-microscope:{.nextflow-primary} Nextflow dla nauki {.mt-1}

    Naucz się stosować koncepcje i komponenty przedstawione w 'Hello Nextflow' do konkretnych zastosowań naukowych.

    ??? courses "**Nextflow dla genomiki** (wykrywanie wariantów)"

        Dla badaczy, którzy chcą nauczyć się tworzyć własne pipeline'y genomiczne. Poprzez przypadek wykrywania wariantów pokazuje, jak budować prosty, ale funkcjonalny pipeline genomiczny.

        [Rozpocznij szkolenie Nextflow dla genomiki :material-arrow-right:](nf4_science/genomics/){ .md-button .md-button--secondary }

    ??? courses "**Nextflow dla RNAseq** (bulk RNAseq)"

        Dla badaczy, którzy chcą nauczyć się tworzyć własne pipeline'y RNAseq. Poprzez przypadek przetwarzania bulk RNAseq demonstruje, jak budować prosty, ale funkcjonalny pipeline analityczny.

        [Rozpocznij szkolenie Nextflow dla RNAseq :material-arrow-right:](nf4_science/rnaseq/){ .md-button .md-button--secondary }

    ??? courses "**Nextflow dla obrazowania** (spatial omics)"

        Dla badaczy zajmujących się obrazowaniem i spatial omics, którzy chcą nauczyć się uruchamiać i dostosowywać workflow'y analityczne. Oparty na nf-core/molkart, kurs demonstruje konfigurację, wykonywanie oraz zarządzanie wejściami dla workflow'ów Nextflow.

        [Rozpocznij szkolenie Nextflow dla obrazowania :material-arrow-right:](nf4_science/imaging/){ .md-button .md-button--secondary }

-   :material-run:{ .lg .middle } __Ścieżka zaawansowana__

    ---

    ### :material-bridge:{.nextflow-primary} Od Nextflow do nf-core {.mt-1}

    Naucz się wykorzystywać kod i najlepsze praktyki z projektu społeczności [nf-core](https://nf-co.re/).

    Te kursy pomagają przejść od podstaw Nextflow do najlepszych praktyk nf-core. Zrozum, jak i dlaczego społeczność buduje swoje workflow'y, oraz jak możesz współtworzyć i ponownie stosować te techniki.

    ??? courses "**Hello nf-core:** Rozpocznij pracę z nf-core"

        Dla programistów, którzy chcą nauczyć się uruchamiać i tworzyć workflow'y zgodne z [nf-core](https://nf-co.re/). Kurs szczegółowo omawia strukturę tych workflow'ów, umożliwiając budowanie prostych, ale funkcjonalnych rozwiązań. Uczestnicy poznają szablon, najlepsze praktyki deweloperskie oraz istniejące moduły społeczności.

        [Rozpocznij szkolenie Hello nf-core :material-arrow-right:](hello_nf-core/index.md){ .md-button .md-button--secondary }

    ---

    ### :material-rocket-launch:{.nextflow-primary} Zaawansowane szkolenie Nextflow {.mt-1}

    Poznaj zaawansowane koncepcje i mechanizmy tworzenia i wdrażania pipeline'ów Nextflow dla rzeczywistych przypadków użycia.

    ??? courses "**Side Quests:** Głębokie zanurzenia w samodzielne tematy"

        Samodzielne mini-kursy przeznaczone dla programistów Nextflow, którzy chcą poszerzyć swój zakres i/lub pogłębić umiejętności w określonych tematach. Są przedstawione liniowo, ale można je realizować w dowolnej kolejności (zobacz zależności w przeglądzie każdego mini-kursu).

        [Przeglądaj Side Quests :material-arrow-right:](side_quests/){ .md-button .md-button--secondary }

    ??? courses "**Kolekcje szkoleniowe:** Polecane ścieżki nauki przez Side Quests"

        Kolekcje szkoleniowe łączą wiele Side Quests, aby zapewnić kompleksowe doświadczenie edukacyjne wokół określonego tematu lub przypadku użycia.

        [Przeglądaj kolekcje szkoleniowe :material-arrow-right:](training_collections/){ .md-button .md-button--secondary }

</div>

!!! info "Szukasz zarchiwizowanych materiałów szkoleniowych?"

    Starsze materiały szkoleniowe (Fundamentals Training, Advanced Training i inne eksperymentalne kursy) zostały usunięte z portalu szkoleniowego, ponieważ są niekompatybilne ze ścisłą składnią Nextflow 3.0.
    Jeśli potrzebujesz dostępu do tych materiałów, są one dostępne w [historii git](https://github.com/nextflow-io/training) sprzed stycznia 2026 roku.

---

<div markdown class="homepage_logos">

![Seqera](assets/img/seqera_logo.png#only-light)

![Seqera](assets/img/seqera_logo_dark.png#only-dark)

</div>
