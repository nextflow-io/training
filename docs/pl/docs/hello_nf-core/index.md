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
    - Tworzyć podstawowy pipeline kompatybilny z nf-core na podstawie szablonu
    - Uaktualniać prosty workflow Nextflow do standardów nf-core
    - Dodawać moduły nf-core do pipeline'u kompatybilnego z nf-core
    - Współtworzyć własne moduły dla nf-core
    - Walidować dane wejściowe i parametry przy użyciu narzędzi nf-core
  audience_prerequisites:
    - "**Odbiorcy:** Ten kurs jest przeznaczony dla osób, które są już zaznajomione z podstawami Nextflow i chcą nauczyć się korzystać z zasobów i najlepszych praktyk nf-core."
    - "**Umiejętności:** Zakłada się znajomość wiersza poleceń, podstawowych koncepcji skryptowania oraz popularnych formatów plików."
    - "**Kursy:** Wymagane jest ukończenie kursu [Hello Nextflow](../hello_nextflow/index.md) lub równoważnego."
    - "**Dziedzina:** Wszystkie ćwiczenia są niezależne od dziedziny, więc nie jest wymagana wcześniejsza wiedza naukowa."
---

# Hello nf-core

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Tłumaczenie wspomagane przez AI - [dowiedz się więcej i zasugeruj ulepszenia](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

**Hello nf-core to praktyczne wprowadzenie do korzystania z zasobów i najlepszych praktyk nf-core.**

![nf-core logo](./img/nf-core-logo.png)

Pracując nad praktycznymi przykładami i prowadzonymi ćwiczeniami, nauczysz się używać i rozwijać moduły i pipeline'y kompatybilne z nf-core oraz efektywnie wykorzystywać narzędzia nf-core.

Wyniesiesz umiejętności i pewność siebie, aby zacząć rozwijać pipeline'y zgodnie z najlepszymi praktykami nf-core.

<!-- additional_information -->

## Przegląd kursu

Ten kurs został zaprojektowany jako praktyczny, z ćwiczeniami zorientowanymi na cel, ustrukturyzowanymi tak, aby stopniowo wprowadzać informacje.

Zostaniesz zapoznany z [**nf-core**](https://nf-co.re/), inicjatywą społeczności mającą na celu rozwój i utrzymanie wyselekcjonowanego zestawu naukowych pipeline'ów zbudowanych przy użyciu Nextflow, a także odpowiednich narzędzi i wytycznych promujących otwartą współpracę, testowanie i przeglądy kodu ([Nat Biotechnol 38, 276–278 (2020)](https://www.nature.com/articles/s41587-020-0439-x), [Genome Biol 26, 228 (2025)](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-025-03673-9)).

Pipeline'y opracowane przez społeczność nf-core zostały zaprojektowane jako modułowe, skalowalne i przenośne, pozwalając badaczom na łatwe dostosowanie i wykonanie ich przy użyciu własnych danych i zasobów obliczeniowych.
Wytyczne dotyczące najlepszych praktyk egzekwowane przez projekt dodatkowo zapewniają, że te rozwiązania są niezawodne, dobrze udokumentowane i zwalidowane na rzeczywistych zestawach danych.
Pomaga to zwiększyć niezawodność i odtwarzalność analiz naukowych i ostatecznie umożliwia badaczom przyspieszenie ich odkryć naukowych.

Nie omówimy wszystkiego, co można wiedzieć o pipeline'ach nf-core w tym kursie, ponieważ nf-core obejmuje wiele funkcji i konwencji opracowanych przez społeczność przez lata.
Zamiast tego skupimy się na podstawowych koncepcjach, które pomogą Ci zacząć i zrozumieć, jak działa nf-core.

### Plan lekcji

Podzieliliśmy to na pięć części, z których każda skupi się na określonych aspektach korzystania z zasobów nf-core.

| Rozdział kursu                                                  | Podsumowanie                                                                                                                                                               | Szacowany czas |
| --------------------------------------------------------------- | -------------------------------------------------------------------------------------------------------------------------------------------------------------------------- | -------------- |
| [Część 1: Uruchomienie demo pipeline'a](./01_run_demo.md)       | Uruchomienie istniejącego pipeline'a nf-core i zbadanie jego struktury kodu, aby uzyskać wyobrażenie o tym, co wyróżnia te pipeline'y od podstawowych workflow'ów Nextflow | 30 min         |
| [Część 2: Przepisanie Hello dla nf-core](./02_rewrite_hello.md) | Dostosowanie istniejącego workflow'u do szkieletu szablonu nf-core, zaczynając od prostego workflow'u utworzonego w kursie [Hello Nextflow](../hello_nextflow/index.md)    | 60 min         |
| [Część 3: Użycie modułu nf-core](./03_use_module.md)            | Eksploracja biblioteki modułów społeczności i nauka integrowania gotowych, przetestowanych modułów opakowujących popularne narzędzia bioinformatyczne                      | 30 min         |
| [Część 4: Stworzenie modułu nf-core](./04_make_module.md)       | Tworzenie własnego modułu w stylu nf-core przy użyciu określonej struktury, konwencji nazewnictwa i wymagań dotyczących metadanych ustanowionych przez nf-core             | 30 min         |
| [Część 5: Dodanie walidacji wejścia](./05_input_validation.md)  | Implementacja walidacji wejścia zarówno dla parametrów wiersza poleceń, jak i plików danych wejściowych przy użyciu nf-schema                                              | 30 min         |

Pod koniec tego kursu będziesz w stanie wykorzystać ogromne bogactwo zasobów oferowanych przez projekt nf-core.

Gotowy rozpocząć kurs?

[Rozpocznij naukę :material-arrow-right:](00_orientation.md){ .md-button .md-button--primary }
