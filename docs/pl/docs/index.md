---
title: Strona główna
description: Witaj w portalu szkoleń społeczności Nextflow!
hide:
  - toc
  - footer
---

# Szkolenia Nextflow

<div class="grid cards" markdown>

-   :material-book-open-variant:{ .lg .middle } __Kursy samoobsługowe__

    ---

    **Witaj w portalu szkoleń społeczności Nextflow!**

    Kursy szkoleniowe wymienione poniżej zostały zaprojektowane jako materiały samoobsługowe.
    Możesz przejść przez nie samodzielnie w dowolnym momencie, korzystając z internetowego środowiska udostępnianego przez nas za pośrednictwem Github Codespaces lub we własnym środowisku.

    [Przeglądaj kursy :material-arrow-right:](#katalog-kursow-szkoleniowych-nextflow){ .md-button .md-button--primary .mt-1 }

-   :material-information-outline:{ .lg .middle } __Dodatkowe informacje__

    ---

    ??? warning "Zgodność wersji"

        <!-- Any update to this content needs to be copied to the local installation page -->
        **Od stycznia 2026 roku wszystkie nasze kursy szkoleniowe Nextflow wymagają wersji Nextflow 25.10.2 lub nowszej, z aktywowaną ścisłą składnią, chyba że zaznaczono inaczej.**

        Więcej informacji o wymaganiach dotyczących wersji i ścisłej składni znajdziesz w [przewodniku migracji w dokumentacji Nextflow](https://nextflow.io/docs/latest/strict-syntax.html).

        Starsze wersje materiałów szkoleniowych odpowiadające wcześniejszej składni są dostępne za pośrednictwem selektora wersji w pasku menu tej strony.

    ??? terminal "Opcje środowiska"

        Udostępniamy internetowe środowisko szkoleniowe, w którym wszystko, czego potrzebujesz do odbycia szkolenia, jest preinstalowane. Jest ono dostępne za pośrednictwem Github Codespaces (wymaga bezpłatnego konta GitHub).

        [![Open in GitHub Codespaces](https://github.com/codespaces/badge.svg)](https://codespaces.new/nextflow-io/training?quickstart=1&ref=master)

        Jeśli to nie odpowiada Twoim potrzebom, zapoznaj się z innymi [opcjami środowiska](./envsetup/index.md).

    ??? learning "Wydarzenia szkoleniowe"

        Jeśli wolisz wziąć udział w szkoleniu Nextflow w ramach zorganizowanego wydarzenia, istnieje wiele możliwości. Polecamy sprawdzić następujące opcje:

        - **[Training Weeks]()** organizowane kwartalnie przez zespół społeczności
        - **[Seqera Events](https://seqera.io/events/)** obejmują stacjonarne wydarzenia szkoleniowe organizowane przez Seqera (szukaj 'Seqera Sessions' i 'Nextflow Summit')
        - **[Nextflow Ambassadors]()** organizują wydarzenia dla swoich lokalnych społeczności
        - **[nf-core events](https://nf-co.re/events)** obejmują hackathony społeczności

    ??? people "Informacje dla trenerów"

        Jeśli jesteś instruktorem prowadzącym własne szkolenia, możesz swobodnie korzystać z naszych materiałów bezpośrednio z portalu szkoleniowego, o ile odpowiednio je przypisujesz. Zobacz szczegóły w sekcji 'Autorstwo i wkład' poniżej.

        Ponadto chętnie usłyszymy od Ciebie, jak moglibyśmy lepiej wspierać Twoje działania szkoleniowe! Skontaktuj się z nami pod adresem [community@seqera.io](mailto:community@seqera.io) lub na forum społeczności (zobacz stronę [Pomoc](help.md)).

    ??? licensing "Licencja open-source i polityka wkładu"

        [![Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International (CC BY-NC-SA 4.0)](assets/img/cc_by-nc-sa.svg){ align=right }](https://creativecommons.org/licenses/by-nc-sa/4.0/)

        Te materiały szkoleniowe są rozwijane i utrzymywane przez [Seqera](https://seqera.io) i udostępniane na licencji open-source ([CC BY-NC-SA](https://creativecommons.org/licenses/by-nc-sa/4.0/)) z korzyścią dla społeczności. Jeśli chcesz wykorzystać te materiały w sposób wykraczający poza zakres licencji (zwróć uwagę na ograniczenia dotyczące użytku komercyjnego i redystrybucji), skontaktuj się z nami pod adresem [community@seqera.io](mailto:community@seqera.io), aby omówić Twoją prośbę.

        Z zadowoleniem przyjmujemy ulepszenia, poprawki i zgłoszenia błędów od społeczności. Każda strona ma ikonę :material-file-edit-outline: w prawym górnym rogu, która prowadzi do repozytorium kodu, gdzie możesz zgłaszać problemy lub proponować zmiany w materiałach źródłowych szkolenia za pomocą pull requesta. Zobacz `README.md` w repozytorium, aby uzyskać więcej szczegółów.

</div>

!!! note "Tłumaczenie wspomagane przez AI"

    To tłumaczenie zostało stworzone przy użyciu sztucznej inteligencji i zweryfikowane przez ludzkich tłumaczy.
    Zachęcamy do przesyłania opinii i sugestii ulepszeń.
    Więcej informacji znajdziesz w naszym [przewodniku tłumaczenia](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md).

## Katalog kursów szkoleniowych Nextflow

<div class="grid cards" markdown>

-   :material-walk:{ .lg .middle } __Ścieżka wprowadzająca__

    ---

    ### :material-compass:{.nextflow-primary} Nextflow dla początkujących {.mt-1}

    Kursy niezależne od dziedziny, przeznaczone dla osób całkowicie nowych w Nextflow. Każdy kurs składa się z serii modułów szkoleniowych zaprojektowanych tak, aby pomóc uczącym się stopniowo rozwijać swoje umiejętności.

    ??? courses "**Hello Nextflow:** Naucz się tworzyć własne pipeline'y"

        Ten kurs obejmuje podstawowe komponenty języka Nextflow w wystarczającym stopniu szczegółowości, aby umożliwić tworzenie prostych, ale w pełni funkcjonalnych pipeline'ów, a także kluczowe elementy projektowania, rozwoju i konfiguracji pipeline'ów.

        [Rozpocznij szkolenie Hello Nextflow :material-arrow-right:](hello_nextflow/index.md){ .md-button .md-button--secondary }

    ??? courses "**Nextflow Run:** Naucz się uruchamiać istniejące pipeline'y"

        Zwięzłe wprowadzenie do uruchamiania i konfigurowania pipeline'ów Nextflow, oparte na kursie deweloperskim Hello Nextflow, ale z mniejszym naciskiem na kod. Obejmuje wykonywanie, wyjścia, podstawową strukturę kodu i konfigurację dla różnych środowisk obliczeniowych.

        [Rozpocznij szkolenie Nextflow Run :material-arrow-right:](nextflow_run/index.md){ .md-button .md-button--secondary }

    ---

    ### :material-microscope:{.nextflow-primary} Nextflow dla nauki {.mt-1}

    Naucz się stosować koncepcje i komponenty przedstawione w 'Hello Nextflow' do konkretnych przypadków użycia naukowego.

    ??? courses "**Nextflow for Genomics** (wykrywanie wariantów)"

        Dla badaczy, którzy chcą nauczyć się tworzyć własne pipeline'y genomiczne. Kurs wykorzystuje przypadek użycia wykrywania wariantów, aby zademonstrować, jak opracować prosty, ale funkcjonalny pipeline genomiczny.

        [Rozpocznij szkolenie Nextflow for Genomics :material-arrow-right:](nf4_science/genomics/){ .md-button .md-button--secondary }

    ??? courses "**Nextflow for RNAseq** (bulk RNAseq)"

        Dla badaczy, którzy chcą nauczyć się tworzyć własne pipeline'y RNAseq. Kurs wykorzystuje przypadek użycia przetwarzania bulk RNAseq, aby zademonstrować, jak opracować prosty, ale funkcjonalny pipeline RNAseq.

        [Rozpocznij szkolenie Nextflow for RNAseq :material-arrow-right:](nf4_science/rnaseq/){ .md-button .md-button--secondary }

    ??? courses "**Nextflow for Imaging** (omika przestrzenna)"

        Dla badaczy zajmujących się obrazowaniem i omiką przestrzenną, którzy chcą nauczyć się uruchamiać i dostosowywać pipeline'y analityczne. Kurs wykorzystuje pipeline nf-core/molkart, aby zapewnić biologicznie istotny przykład demonstrujący, jak uruchamiać, konfigurować i zarządzać danymi wejściowymi dla workflow'ów Nextflow.

        [Rozpocznij szkolenie Nextflow for Imaging :material-arrow-right:](nf4_science/imaging/){ .md-button .md-button--secondary }

-   :material-run:{ .lg .middle } __Ścieżka zaawansowana__

    ---

    ### :material-bridge:{.nextflow-primary} Od Nextflow do nf-core {.mt-1}

    Naucz się wykorzystywać kod i najlepsze praktyki z projektu społecznościowego [nf-core](https://nf-co.re/).

    Te kursy pomogą Ci przejść od podstaw Nextflow do najlepszych praktyk nf-core.
    Zrozum, jak i dlaczego społeczność nf-core buduje pipeline'y oraz jak możesz wnosić wkład i ponownie wykorzystywać te techniki.

    ??? courses "**Hello nf-core:** Rozpocznij pracę z nf-core"

        Dla deweloperów, którzy chcą nauczyć się uruchamiać i tworzyć pipeline'y zgodne z [nf-core](https://nf-co.re/). Kurs obejmuje strukturę pipeline'ów nf-core w wystarczającym stopniu szczegółowości, aby umożliwić tworzenie prostych, ale w pełni funkcjonalnych pipeline'ów zgodnych z szablonem nf-core i najlepszymi praktykami rozwoju, a także korzystanie z istniejących modułów nf-core.

        [Rozpocznij szkolenie Hello nf-core :material-arrow-right:](hello_nf-core/index.md){ .md-button .md-button--secondary }

    ---

    ### :material-rocket-launch:{.nextflow-primary} Zaawansowane szkolenie Nextflow {.mt-1}

    Poznaj zaawansowane koncepcje i mechanizmy tworzenia i wdrażania pipeline'ów Nextflow w celu rozwiązywania rzeczywistych przypadków użycia.

    ??? courses "**Side Quests:** Głębokie zanurzenia w samodzielne tematy"

        Samodzielne mini-kursy przeznaczone dla deweloperów Nextflow, którzy chcą poszerzyć swój zakres i/lub pogłębić umiejętności w określonych tematach. Są prezentowane liniowo, ale można je realizować w dowolnej kolejności (zobacz zależności w przeglądzie każdego mini-kursu).

        [Przeglądaj Side Quests :material-arrow-right:](side_quests/){ .md-button .md-button--secondary }

    ??? courses "**Training Collections:** Zalecane ścieżki uczenia się przez Side Quests"

        Training Collections łączą wiele Side Quests, aby zapewnić kompleksowe doświadczenie edukacyjne wokół określonego tematu lub przypadku użycia.

        [Przeglądaj Training Collections :material-arrow-right:](training_collections/){ .md-button .md-button--secondary }

</div>

!!! info "Szukasz zarchiwizowanych materiałów szkoleniowych?"

    Starsze materiały szkoleniowe (Fundamentals Training, Advanced Training i inne eksperymentalne kursy) zostały usunięte z portalu szkoleniowego, ponieważ są niezgodne ze ścisłą składnią Nextflow 3.0.
    Jeśli potrzebujesz dostępu do tych materiałów, są one dostępne w [historii git](https://github.com/nextflow-io/training) sprzed stycznia 2026 roku.

---

<div markdown class="homepage_logos">

![Seqera](assets/img/seqera_logo.png#only-light)

![Seqera](assets/img/seqera_logo_dark.png#only-dark)

</div>
