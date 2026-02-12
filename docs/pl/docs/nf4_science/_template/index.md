---
title: Nextflow dla {DOMAIN}
hide:
  - toc
page_type: index_page
index_type: course
additional_information:
  technical_requirements: true
  learning_objectives:
    - Napisać liniowy workflow'a, aby zastosować {METHOD} do pojedynczej próbki
    - Odpowiednio obsługiwać pliki pomocnicze, takie jak {ACCESSORY_FILES}
    - Wykorzystać paradygmat przepływu danych Nextflow'a do zrównoleglenia przetwarzania próbek
    - Zaimplementować agregację wielu próbek przy użyciu odpowiednich operatorów kanałów
  audience_prerequisites:
    - "**Odbiorcy:** Ten kurs jest przeznaczony dla badaczy z dziedziny {DOMAIN} i pokrewnych obszarów, którzy chcą tworzyć lub dostosowywać pipeline'y do analizy danych."
    - "**Umiejętności:** Zakładamy pewną znajomość wiersza poleceń, podstawowych koncepcji skryptowania oraz popularnych formatów plików w {DOMAIN}."
    - "**Wymagania wstępne:** Podstawowe koncepcje i narzędzia Nextflow'a omówione w kursie [Hello Nextflow](../../hello_nextflow/)."
---

# Nextflow dla {DOMAIN}

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Tłumaczenie wspomagane przez AI - [dowiedz się więcej i zasugeruj ulepszenia](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

**Praktyczny kurs stosowania Nextflow'a do rzeczywistego przypadku użycia w {DOMAIN}: {METHOD_SHORT_DESCRIPTION}.**

Ten kurs opiera się na szkoleniu dla początkujących [Hello Nextflow](../../hello_nextflow/) i pokazuje, jak używać Nextflow'a w specyficznym kontekście dziedziny {DOMAIN}.
Zaimplementujesz pipeline {METHOD} przy użyciu narzędzi [{TOOL_A}]({TOOL_A_URL}) i [{TOOL_B}]({TOOL_B_URL}).

<!-- additional_information -->

## Przegląd kursu

Ten kurs ma charakter praktyczny, z ćwiczeniami zorientowanymi na cel i ustrukturyzowanymi tak, aby stopniowo wprowadzać nowe informacje.

Zaczniesz od ręcznego uruchamiania narzędzi analitycznych w terminalu, aby zrozumieć metodologię, a następnie stopniowo zbudujesz pipeline Nextflow'a, który automatyzuje i skaluje analizę.

### Plan lekcji

Podzieliliśmy to na trzy części, z których każda koncentruje się na konkretnych aspektach stosowania Nextflow'a do przypadku użycia w {DOMAIN}.

| Rozdział kursu                                            | Podsumowanie                                                                                                  | Szacowany czas |
| --------------------------------------------------------- | ------------------------------------------------------------------------------------------------------------- | -------------- |
| [Część 1: Przegląd metody](./01_method.md)                | Zrozumienie metodologii {METHOD} i ręczne uruchamianie narzędzi                                               | 30 min         |
| [Część 2: Przetwarzanie pojedynczej próbki](./02_single_sample.md) | Budowanie pipeline'u, który {PART2_SUMMARY}, a następnie skalowanie do wielu próbek                           | 60 min         |
| [Część 3: Agregacja wielu próbek](./03_multi_sample.md)  | Dodawanie {AGGREGATION_SUMMARY} dla wielu próbek przy użyciu operatorów kanałów do agregacji wyników próbek | 45 min         |

Pod koniec tego kursu będziesz w stanie zastosować podstawowe koncepcje i narzędzia Nextflow'a do typowego przypadku użycia w {DOMAIN}.

Gotowy, aby rozpocząć kurs?

[Zacznij :material-arrow-right:](00_orientation.md){ .md-button .md-button--primary }

<!-- Clearfix for float -->
<div style="content: ''; clear: both; display: table;"></div>
