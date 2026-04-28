---
title: Strona główna
description: Witaj w portalu szkoleniowym społeczności Nextflow!
hide:
  - toc
  - footer
---

# Szkolenie Nextflow

<div class="grid cards" markdown>

-   :material-book-open-variant:{ .lg .middle } __Kursy do samodzielnej nauki__

    ---

    **Witaj w portalu szkoleniowym społeczności Nextflow!**

    Kursy wymienione poniżej zostały zaprojektowane jako materiały do samodzielnej nauki.
    Możesz je przerabiać we własnym tempie, w dowolnym momencie — zarówno w środowisku webowym udostępnianym przez nas za pośrednictwem Github Codespaces, jak i we własnym środowisku.

    [Przeglądaj kursy :material-arrow-right:](#catalog-of-nextflow-training-courses){ .md-button .md-button--primary .mt-1 }

-   :material-information-outline:{ .lg .middle } __Dodatkowe informacje__

    ---

    ??? warning "Zgodność wersji"

        <!-- Wszelkie zmiany tej treści należy skopiować na stronę instalacji lokalnej -->
        **Od stycznia 2026 roku wszystkie nasze kursy szkoleniowe Nextflow wymagają Nextflow w wersji 25.10.2 lub nowszej z aktywną ścisłą składnią, o ile nie zaznaczono inaczej.**

        Więcej informacji na temat wymagań dotyczących wersji i ścisłej składni znajdziesz w [przewodniku migracji dokumentacji Nextflow](https://nextflow.io/docs/latest/strict-syntax.html).

        Starsze wersje materiałów szkoleniowych odpowiadające poprzedniej składni są dostępne za pomocą selektora wersji na pasku menu tej strony.

    ??? terminal "Opcje środowiska"

        Udostępniamy webowe środowisko szkoleniowe, w którym wszystko, czego potrzebujesz do nauki, jest już zainstalowane — dostępne przez Github Codespaces (wymaga bezpłatnego konta GitHub).

        [![Otwórz w GitHub Codespaces](https://github.com/codespaces/badge.svg)](https://codespaces.new/nextflow-io/training?quickstart=1&ref=master)

        Jeśli to rozwiązanie Ci nie odpowiada, zapoznaj się z innymi [opcjami środowiska](./envsetup/index.md).

    ??? learning "Wydarzenia szkoleniowe"

        Jeśli wolisz uczyć się Nextflow w ramach zorganizowanego wydarzenia, masz do wyboru wiele możliwości. Polecamy sprawdzić następujące opcje:

        - **[Training Weeks]()** — organizowane kwartalnie przez zespół społeczności
        - **[Seqera Events](https://seqera.io/events/)** — stacjonarne wydarzenia szkoleniowe organizowane przez Seqera (szukaj „Seqera Sessions" i „Nextflow Summit")
        - **[Nextflow Ambassadors]()** — organizują wydarzenia dla swoich lokalnych społeczności
        - **[nf-core events](https://nf-co.re/events)** — hackathony społecznościowe

    ??? people "Informacje dla prowadzących szkolenia"

        Jeśli jesteś instruktorem prowadzącym własne szkolenia, możesz korzystać z naszych materiałów bezpośrednio z portalu szkoleniowego, pod warunkiem odpowiedniego podania źródła. Szczegóły znajdziesz w sekcji „Licencja i wkład" poniżej.

        Chętnie dowiemy się też, jak możemy lepiej wspierać Twoje działania szkoleniowe! Skontaktuj się z nami pod adresem [community@seqera.io](mailto:community@seqera.io) lub na forum społeczności (zob. strona [Pomoc](help.md)).

    ??? licensing "Licencja open-source i polityka wkładu"

        [![Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International (CC BY-NC-SA 4.0)](assets/img/cc_by-nc-sa.svg){ align=right }](https://creativecommons.org/licenses/by-nc-sa/4.0/)

        Materiały szkoleniowe są opracowywane i utrzymywane przez [Seqera](https://seqera.io) i udostępniane na licencji open-source ([CC BY-NC-SA](https://creativecommons.org/licenses/by-nc-sa/4.0/)) z korzyścią dla społeczności. Jeśli chcesz wykorzystać te materiały w sposób wykraczający poza zakres licencji (zwróć uwagę na ograniczenia dotyczące użytku komercyjnego i redystrybucji), skontaktuj się z nami pod adresem [community@seqera.io](mailto:community@seqera.io), aby omówić swoją prośbę.

        Chętnie przyjmujemy ulepszenia, poprawki i zgłoszenia błędów od społeczności. Na każdej stronie w prawym górnym rogu znajduje się ikona :material-file-edit-outline:, która prowadzi do repozytorium kodu — tam możesz zgłaszać problemy lub proponować zmiany w materiałach szkoleniowych za pomocą pull request. Więcej szczegółów znajdziesz w pliku `README.md` w repozytorium.

</div>

!!! note "Tłumaczenie wspomagane przez AI"

    To tłumaczenie zostało stworzone przy użyciu sztucznej inteligencji i zweryfikowane przez ludzkich tłumaczy.
    Zachęcamy do przesyłania opinii i sugestii ulepszeń.
    Więcej informacji znajdziesz w naszym [przewodniku tłumaczenia](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md).

## Katalog kursów szkoleniowych Nextflow

<div class="grid cards" markdown>

-   :material-walk:{ .lg .middle } __Ścieżka dla początkujących__

    ---

    ### :material-compass:{.nextflow-primary} Nextflow dla nowych użytkowników {.mt-1}

    Kursy niezależne od dziedziny, przeznaczone dla osób zupełnie nowych w Nextflow. Każdy kurs składa się z serii modułów szkoleniowych, które pomagają stopniowo rozwijać umiejętności.

    ??? courses "**Hello Nextflow:** Naucz się tworzyć własne pipeline'y"

        Kurs obejmuje podstawowe komponenty języka Nextflow w stopniu wystarczającym do tworzenia prostych, ale w pełni funkcjonalnych pipeline'ów, a także kluczowe elementy projektowania, rozwijania i konfigurowania pipeline'ów.

        [Rozpocznij szkolenie Hello Nextflow :material-arrow-right:](hello_nextflow/index.md){ .md-button .md-button--secondary }

    ??? courses "**Nextflow Run:** Naucz się uruchamiać istniejące pipeline'y"

        Zwięzłe wprowadzenie do uruchamiania i konfigurowania pipeline'ów Nextflow, oparte na kursie deweloperskim Hello Nextflow, ale z mniejszym naciskiem na kod. Obejmuje uruchamianie, wyniki, podstawową strukturę kodu oraz konfigurację dla różnych środowisk obliczeniowych.

        [Rozpocznij szkolenie Nextflow Run :material-arrow-right:](nextflow_run/index.md){ .md-button .md-button--secondary }

    ---

    ### :material-microscope:{.nextflow-primary} Nextflow dla nauki {.mt-1}

    Naucz się stosować koncepcje i komponenty przedstawione w „Hello Nextflow" do konkretnych zastosowań naukowych.

    ??? courses "**Nextflow for Genomics** (wykrywanie wariantów)"

        Dla badaczy, którzy chcą nauczyć się tworzyć własne pipeline'y genomiczne. Kurs wykorzystuje przypadek użycia wykrywania wariantów, aby pokazać, jak zbudować prosty, ale funkcjonalny pipeline genomiczny.

        [Rozpocznij szkolenie Nextflow for Genomics :material-arrow-right:](nf4_science/genomics/index.md){ .md-button .md-button--secondary }

    ??? courses "**Nextflow for RNAseq** (bulk RNAseq)"

        Dla badaczy, którzy chcą nauczyć się tworzyć własne pipeline'y RNAseq. Kurs wykorzystuje przypadek użycia przetwarzania bulk RNAseq, aby pokazać, jak zbudować prosty, ale funkcjonalny pipeline RNAseq.

        [Rozpocznij szkolenie Nextflow for RNAseq :material-arrow-right:](nf4_science/rnaseq/index.md){ .md-button .md-button--secondary }

    ??? courses "**Nextflow for Imaging** (spatial omics)"

        Dla badaczy zajmujących się obrazowaniem i spatial omics, którzy chcą nauczyć się uruchamiać i dostosowywać pipeline'y analityczne. Kurs wykorzystuje pipeline nf-core/molkart, aby na biologicznie istotnym przykładzie pokazać, jak uruchamiać Nextflow'a, konfigurować go i zarządzać danymi wejściowymi.

        [Rozpocznij szkolenie Nextflow for Imaging :material-arrow-right:](nf4_science/imaging/index.md){ .md-button .md-button--secondary }

-   :material-run:{ .lg .middle } __Ścieżka zaawansowana__

    ---

    ### :material-bridge:{.nextflow-primary} Od Nextflow do nf-core {.mt-1}

    Naucz się korzystać z kodu i najlepszych praktyk projektu społecznościowego [nf-core](https://nf-co.re/).

    Te kursy przeprowadzą Cię od podstaw Nextflow do najlepszych praktyk nf-core.
    Dowiedz się, jak i dlaczego społeczność nf-core buduje pipeline'y, oraz jak możesz wnosić wkład i ponownie wykorzystywać te techniki.

    ??? courses "**Hello nf-core:** Pierwsze kroki z nf-core"

        Dla deweloperów, którzy chcą nauczyć się uruchamiać i rozwijać pipeline'y zgodne z [nf-core](https://nf-co.re/). Kurs omawia strukturę pipeline'ów nf-core w stopniu wystarczającym do tworzenia prostych, ale w pełni funkcjonalnych pipeline'ów zgodnych z szablonem i najlepszymi praktykami nf-core, a także do korzystania z istniejących modułów nf-core.

        [Rozpocznij szkolenie Hello nf-core :material-arrow-right:](hello_nf-core/index.md){ .md-button .md-button--secondary }

    ---

    ### :material-rocket-launch:{.nextflow-primary} Zaawansowane szkolenie Nextflow {.mt-1}

    Poznaj zaawansowane koncepcje i mechanizmy tworzenia oraz wdrażania pipeline'ów Nextflow w rzeczywistych zastosowaniach.

    ??? courses "**Side Quests:** Dogłębne omówienie wybranych tematów"

        Samodzielne mini-kursy przeznaczone dla deweloperów Nextflow, którzy chcą poszerzyć swoje horyzonty lub pogłębić wiedzę w konkretnych obszarach. Są ułożone liniowo, ale można je realizować w dowolnej kolejności (zob. zależności w przeglądzie każdego mini-kursu).

        [Przeglądaj Side Quests :material-arrow-right:](side_quests/index.md){ .md-button .md-button--secondary }

    ??? courses "**Training Collections:** Zalecane ścieżki nauki przez Side Quests"

        Training Collections łączą wiele Side Quests, aby zapewnić kompleksowe doświadczenie edukacyjne wokół określonego tematu lub przypadku użycia.

        [Przeglądaj Training Collections :material-arrow-right:](training_collections/index.md){ .md-button .md-button--secondary }

</div>

!!! info "Szukasz archiwalnych materiałów szkoleniowych?"

    Starsze materiały szkoleniowe (Fundamentals Training, Advanced Training i inne kursy eksperymentalne) zostały usunięte z portalu szkoleniowego, ponieważ są niezgodne ze ścisłą składnią Nextflow 3.0.
    Jeśli potrzebujesz dostępu do tych materiałów, są one dostępne w [historii git](https://github.com/nextflow-io/training) sprzed stycznia 2026 roku.

---

<div markdown class="homepage_logos">

![Seqera](assets/img/seqera_logo.png#only-light)

![Seqera](assets/img/seqera_logo_dark.png#only-dark)

</div>
