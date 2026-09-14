---
title: Strona główna
description: Witaj w portalu szkoleń społeczności Nextflow!
hide:
  - toc
  - footer
---

# Szkolenie Nextflow

<div class="grid cards" markdown>

-   :material-book-open-variant:{ .lg .middle } __Kursy do samodzielnej nauki__

    ---

    **Witaj w portalu szkoleń społeczności Nextflow!**

    Przerabiaj poniższe kursy we własnym tempie — w naszym środowisku webowym lub swoim własnym.
    Każdy kurs jest praktyczny i składa się z ćwiczeń zorientowanych na cel, które możesz ukończyć samodzielnie.

    [Przeglądaj kursy :material-arrow-down:](#catalog-of-nextflow-training-courses){ .md-button .md-button--primary .mt-1 }

-   :material-account-group-outline:{ .lg .middle } __Wydarzenia szkoleniowe__

    ---

    **Szukasz czegoś więcej niż samodzielnej nauki?**

    Znajdź ustrukturyzowane wydarzenia szkoleniowe, wskazówki dotyczące prowadzenia własnych szkoleń oraz naszą licencję open-source i politykę wkładu.

    [Zobacz wydarzenia szkoleniowe :material-arrow-right:](training_events.md){ .md-button .md-button--secondary .mt-1 }

</div>

!!! note "Tłumaczenie wspomagane przez AI"

    To tłumaczenie zostało stworzone przy użyciu sztucznej inteligencji i zweryfikowane przez ludzkich tłumaczy.
    Zachęcamy do przesyłania opinii i sugestii ulepszeń.
    Więcej informacji znajdziesz w naszym [przewodniku tłumaczenia](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md).

## Catalog of Nextflow training courses

<div class="grid cards" markdown>

-   :material-account:{ .lg .middle } __Dla użytkowników__

    ---

    ### :material-play-circle:{.nextflow-primary} Uruchamianie pipeline'ów {.mt-1}

    Naucz się uruchamiać istniejące pipeline'y bez pisania kodu.

    ??? courses "**Nextflow Run:** Uruchamianie pipeline'ów z Nextflow"

        Szybkie wprowadzenie do uruchamiania pipeline'ów Nextflow, które nie wymaga znajomości kodu. Obejmuje uruchamianie pipeline'ów, pobieranie wyników, korzystanie z kontenerów oraz podstawową konfigurację wykonania.

        [Zobacz szkolenie :material-arrow-right:](nextflow_run/index.md){ .md-button .md-button--secondary }

    ??? courses "**Use nf-core:** Znajdowanie i uruchamianie pipeline'ów tworzonych przez społeczność"

        Szybkie wprowadzenie do wyszukiwania, uruchamiania i konfigurowania pipeline'ów z projektu społecznościowego nf-core — zaczynając od minimalnego demonstracyjnego pipeline'u, a kończąc na pipeline'ie do analizy produkcyjnej.

        [Zobacz szkolenie :material-arrow-right:](nfcore_use/index.md){ .md-button .md-button--secondary }

    ??? courses "**Scale with Seqera:** Uruchamianie i monitorowanie pipeline'ów na dużą skalę"

        Praktyczne wprowadzenie do uruchamiania i monitorowania pipeline'ów Nextflow za pomocą Seqera Platform — zarówno przez interfejs webowy, jak i wiersz poleceń.

        [Zobacz szkolenie :material-arrow-right:](seqera_scale/index.md){ .md-button .md-button--secondary }

    ---

    ### :material-tune:{.nextflow-primary} Zarządzanie wykonaniem {.mt-1}

    Naucz się efektywnie zarządzać wykonaniem pipeline'ów.

    ??? courses "**Configure Execution:** Konfigurowanie zasobów, ponownych prób i profili wykonania"

        Praktyczne wprowadzenie do konfigurowania wykonania pipeline'ów Nextflow: dostosowywanie do różnych środowisk obliczeniowych, kontrolowanie przydziału zasobów i ponownych prób oraz przełączanie między predefiniowanymi profilami konfiguracyjnymi.

        [Zobacz szkolenie :material-arrow-right:](config_exec/index.md){ .md-button .md-button--secondary }

    !!! info compact "Więcej tematów już wkrótce"

        W tej sekcji planowane są zagadnienia dotyczące strojenia wydajności, wykonania na HPC/w chmurze i inne.
        Zagłosuj na to, co chcesz zobaczyć jako następne, w naszej [krótkiej ankiecie](https://seqera.typeform.com/to/JCs91e8v).

-   :material-code-tags:{ .lg .middle } __Dla deweloperów__

    ---

    ### :material-wrench:{.nextflow-primary} Pisanie pipeline'ów {.mt-1}

    Naucz się tworzyć własne pipeline'y Nextflow.

    ??? courses "**Hello Nextflow:** Tworzenie własnych pipeline'ów od podstaw"

        Kurs obejmuje podstawowe komponenty języka Nextflow w stopniu wystarczającym do tworzenia prostych, ale w pełni funkcjonalnych pipeline'ów, a także kluczowe elementy projektowania, rozwijania i konfigurowania pipeline'ów.

        [Zobacz szkolenie :material-arrow-right:](hello_nextflow/index.md){ .md-button .md-button--secondary }

    ??? courses "**Build with nf-core:** Korzystanie z narzędzi i zasad nf-core"

        Dla deweloperów Nextflow, którzy chcą nauczyć się tworzyć pipeline'y zgodne z [nf-core](https://nf-co.re/).
        Kurs omawia strukturę pipeline'ów nf-core w stopniu wystarczającym do tworzenia prostych, ale w pełni funkcjonalnych pipeline'ów, które korzystają z szablonu nf-core i najlepszych praktyk deweloperskich, a także z istniejących modułów nf-core.

        [Zobacz szkolenie :material-arrow-right:](nfcore_build/index.md){ .md-button .md-button--secondary }

    ??? catalog "**Side Quests:** Zagłęb się w zaawansowane tematy Nextflow"

        Zbiór samodzielnych mini-kursów przeznaczonych dla deweloperów Nextflow, którzy chcą poszerzyć swoje horyzonty lub pogłębić umiejętności w konkretnych obszarach.
        Są ułożone liniowo, ale można je przerabiać w dowolnej kolejności (zależności znajdziesz w opisie każdego mini-kursu).

        [Przeglądaj Side Quests :material-arrow-right:](side_quests/index.md){ .md-button .md-button--secondary }

    ---

    ### :material-microscope:{.nextflow-primary} Nextflow dla nauki {.mt-1}

    Naucz się tworzyć pipeline'y Nextflow do konkretnych zastosowań naukowych.

    ??? courses "**Genomics:** Tworzenie pipeline'u do wywoływania wariantów"

        Kurs dla badaczy, którzy chcą nauczyć się tworzyć własne pipeline'y genomiczne — na przykładzie wywoływania wariantów, ilustrującego kluczowe wzorce deweloperskie w Nextflow.

        [Zobacz szkolenie :material-arrow-right:](nf4_science/genomics/index.md){ .md-button .md-button--secondary }

    ??? courses "**RNAseq:** Tworzenie pipeline'u do przetwarzania bulk RNAseq"

        Kurs dla badaczy, którzy chcą nauczyć się tworzyć własne pipeline'y RNAseq — na przykładzie przetwarzania bulk RNAseq, ilustrującego kluczowe wzorce deweloperskie w Nextflow.

        [Zobacz szkolenie :material-arrow-right:](nf4_science/rnaseq/index.md){ .md-button .md-button--secondary }

    ??? courses "**Bioimaging:** Uruchamianie i konfigurowanie pipeline'ów do obrazowania"

        Kurs dla badaczy, którzy chcą nauczyć się uruchamiać i konfigurować pipeline'y do bioimagingu — na przykładzie nf-core/molkart, ilustrującego kluczowe wzorce użytkowania Nextflow.

        [Zobacz szkolenie :material-arrow-right:](nf4_science/imaging/index.md){ .md-button .md-button--secondary }

</div>

## Konfiguracja i pomoc

<div class="grid cards mb-4" markdown>

-   :material-cog-outline:{ .lg .middle } __Środowisko szkoleniowe__

    ---

    Opcje konfiguracji środowiska do szkoleń z Nextflow.

    [Zobacz środowiska szkoleniowe :material-arrow-right:](envsetup/index.md){ .md-button .md-button--secondary }

-   :material-tag-outline:{ .lg .middle } __Wersje Nextflow__

    ---

    Informacje o wersjach składni Nextflow i zarządzaniu ich ewolucją.

    [Sprawdź wymagania dotyczące wersji :material-arrow-right:](info/nxf_versions.md){ .md-button .md-button--secondary }

-   :material-file-code-outline:{ .lg .middle } __Pipeline Hello__

    ---

    Omówienie tego, co robi pipeline Hello i jak jest zbudowany.

    [Przeczytaj omówienie :material-arrow-right:](info/hello_pipeline.md){ .md-button .md-button--secondary }

-   :material-lifebuoy:{ .lg .middle } __Uzyskiwanie pomocy__

    ---

    Przydatne zasoby, gdy napotkasz problem podczas szkolenia z Nextflow.

    [Znajdź pomoc :material-arrow-right:](help.md){ .md-button .md-button--secondary }

</div>
