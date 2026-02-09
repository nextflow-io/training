---
title: Hello nf-core
hide:
  - toc
page_type: index_page
index_type: course
additional_information:
  technical_requirements: true
  learning_objectives:
    - Pobierać, uruchamiać i zarządzać wykonywaniem pipeline'ów nf-core
    - Opisywać strukturę kodu i organizację projektów pipeline'ów nf-core
    - Tworzyć podstawowe pipeline'y kompatybilne z nf-core na podstawie szablonu
    - Ulepszać zwykłe workflow'y Nextflow'a, aby spełniały standardy nf-core
    - Dodawać moduły nf-core do pipeline'ów kompatybilnych z nf-core
    - Wnosić własne moduły do nf-core
    - Walidować dane wejściowe i parametry przy użyciu narzędzi nf-core
  audience_prerequisites:
    - "**Odbiorcy:** Ten kurs jest przeznaczony dla osób, które znają już podstawy Nextflow'a i chcą nauczyć się korzystać z zasobów i najlepszych praktyk nf-core."
    - "**Umiejętności:** Zakładamy znajomość wiersza poleceń, podstawowych koncepcji skryptowania oraz popularnych formatów plików."
    - "**Kursy:** Musisz ukończyć kurs [Hello Nextflow](../hello_nextflow/index.md) lub jego odpowiednik."
    - "**Dziedzina:** Wszystkie ćwiczenia są niezależne od dziedziny naukowej, więc nie jest wymagana wcześniejsza wiedza specjalistyczna."
---

# Hello nf-core

**Hello nf-core to praktyczne wprowadzenie do korzystania z zasobów i najlepszych praktyk nf-core.**

![logo nf-core](./img/nf-core-logo.png#only-light)
![logo nf-core](./img/nf-core-logo-darkbg.png#only-dark)

Pracując nad praktycznymi przykładami i kierowanymi ćwiczeniami, nauczysz się używać i rozwijać moduły oraz pipeline'y kompatybilne z nf-core, a także efektywnie wykorzystywać narzędzia nf-core.

Wyniesiesz umiejętności i pewność siebie potrzebne do rozpoczęcia tworzenia pipeline'ów zgodnie z najlepszymi praktykami nf-core.

<!-- additional_information -->

## Przegląd kursu

Ten kurs został zaprojektowany jako praktyczny, z ćwiczeniami ukierunkowanymi na cel i ustrukturyzowanymi tak, aby stopniowo wprowadzać nowe informacje.

Zostaniesz zapoznany z [**nf-core**](https://nf-co.re/) – wspólnym wysiłkiem społeczności mającym na celu rozwijanie i utrzymywanie wyselekcjonowanego zestawu pipeline'ów naukowych zbudowanych przy użyciu Nextflow'a, a także odpowiednich narzędzi i wytycznych promujących otwarte tworzenie, testowanie i wzajemną ocenę ([Nat Biotechnol 38, 276–278 (2020)](https://www.nature.com/articles/s41587-020-0439-x), [Genome Biol 26, 228 (2025)](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-025-03673-9)).

Pipeline'y tworzone przez społeczność nf-core są zaprojektowane jako modułowe, skalowalne i przenośne, co pozwala badaczom łatwo je dostosowywać i uruchamiać przy użyciu własnych danych i zasobów obliczeniowych.
Wytyczne dotyczące najlepszych praktyk egzekwowane przez projekt dodatkowo zapewniają, że pipeline'y są solidne, dobrze udokumentowane i zwalidowane na rzeczywistych zestawach danych.
Pomaga to zwiększyć niezawodność i powtarzalność analiz naukowych, a ostatecznie umożliwia badaczom przyspieszenie ich odkryć naukowych.

W tym kursie nie omówimy wszystkiego, co można wiedzieć o pipeline'ach nf-core, ponieważ nf-core obejmuje wiele funkcji i konwencji wypracowanych przez społeczność przez lata.
Zamiast tego skupimy się na podstawowych koncepcjach, które pomogą Ci zacząć i zrozumieć, jak działa nf-core.

### Plan lekcji

Podzieliliśmy kurs na pięć części, z których każda skupia się na konkretnych aspektach korzystania z zasobów nf-core.

| Rozdział kursu                                             | Podsumowanie                                                                                                                                                                      | Szacowany czas |
| ---------------------------------------------------------- | --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- | -------------- |
| [Część 1: Uruchom demonstracyjny pipeline](./01_run_demo.md)            | Uruchom istniejący pipeline nf-core i zbadaj jego strukturę kodu, aby zorientować się, co odróżnia te pipeline'y od podstawowych workflow'ów Nextflow'a               | 30 min            |
| [Część 2: Przepisz Hello dla nf-core](./02_rewrite_hello.md) | Dostosuj istniejący workflow do szkieletu szablonu nf-core, zaczynając od prostego workflow'a stworzonego w kursie [Hello Nextflow](../hello_nextflow/index.md) | 60 min            |
| [Część 3: Użyj modułu nf-core](./03_use_module.md)        | Poznaj bibliotekę modułów społeczności i naucz się integrować gotowe, przetestowane moduły opakowujące popularne narzędzia bioinformatyczne                                       | 30 min            |
| [Część 4: Stwórz moduł nf-core](./04_make_module.md)      | Stwórz własny moduł w stylu nf-core, używając specyficznej struktury, konwencji nazewnictwa i wymagań dotyczących metadanych ustanowionych przez nf-core                            | 30 min            |
| [Część 5: Dodaj walidację danych wejściowych](./05_input_validation.md)   | Zaimplementuj walidację danych wejściowych zarówno dla parametrów wiersza poleceń, jak i plików danych wejściowych przy użyciu nf-schema                                                                   | 30 min            |

Pod koniec tego kursu będziesz w stanie wykorzystać ogromne bogactwo zasobów oferowanych przez projekt nf-core.

Gotowy, aby rozpocząć kurs?

[Zacznij naukę :material-arrow-right:](00_orientation.md){ .md-button .md-button--primary }
