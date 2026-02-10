---
title: Nextflow dla genomiki
hide:
  - toc
page_type: index_page
index_type: course
additional_information:
  technical_requirements: true
  learning_objectives:
    - Napisać liniowy workflow do wykrywania wariantów w pojedynczej próbce
    - Odpowiednio obsługiwać pliki pomocnicze, takie jak pliki indeksów i zasoby genomu referencyjnego
    - Wykorzystać paradygmat przepływu danych Nextflow'a do zrównoleglenia wykrywania wariantów w próbkach
    - Zaimplementować wspólne wywoływanie wariantów dla wielu próbek przy użyciu odpowiednich operatorów kanałów
  audience_prerequisites:
    - "**Odbiorcy:** Ten kurs jest przeznaczony dla badaczy w dziedzinie genomiki i pokrewnych dziedzinach, którzy chcą tworzyć lub dostosowywać pipeline'y do analizy danych."
    - "**Umiejętności:** Zakładamy pewną znajomość wiersza poleceń, podstawowych koncepcji skryptowania oraz popularnych formatów plików genomicznych."
    - "**Wymagania wstępne:** Podstawowe koncepcje i narzędzia Nextflow'a omówione w kursie [Hello Nextflow](../../hello_nextflow/)."
---

# Nextflow dla genomiki

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Tłumaczenie wspomagane przez AI - [dowiedz się więcej i zasugeruj ulepszenia](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

**Praktyczny kurs stosowania Nextflow'a w rzeczywistym przypadku użycia z genomiki: wykrywanie wariantów za pomocą GATK.**

Ten kurs opiera się na szkoleniu dla początkujących [Hello Nextflow](../../hello_nextflow/) i pokazuje, jak używać Nextflow'a w konkretnym kontekście dziedziny genomiki.
Zaimplementujesz pipeline do wykrywania wariantów za pomocą [GATK](https://gatk.broadinstitute.org/) (Genome Analysis Toolkit), szeroko stosowanego pakietu oprogramowania do analizy danych z sekwencjonowania o wysokiej przepustowości.

<!-- additional_information -->

## Przegląd kursu

Ten kurs jest praktyczny, z ćwiczeniami zorientowanymi na cel i ustrukturyzowanymi tak, aby stopniowo wprowadzać informacje.

Zaczniesz od ręcznego uruchamiania narzędzi do wykrywania wariantów w terminalu, aby zrozumieć metodologię, a następnie stopniowo zbudujesz pipeline Nextflow'a, który automatyzuje i skaluje analizę.

### Plan lekcji

Podzieliliśmy to na trzy części, z których każda koncentruje się na konkretnych aspektach stosowania Nextflow'a w przypadku użycia z genomiki.

| Rozdział kursu                                                                 | Podsumowanie                                                                                                     | Szacowany czas |
| ------------------------------------------------------------------------------ | ---------------------------------------------------------------------------------------------------------------- | -------------- |
| [Część 1: Przegląd metody](./01_method.md)                                     | Zrozumienie metodologii wykrywania wariantów i ręczne uruchamianie narzędzi                                      | 30 min         |
| [Część 2: Wykrywanie wariantów w próbkach](./02_per_sample_variant_calling.md) | Budowanie pipeline'u, który indeksuje pliki BAM i wykrywa warianty, a następnie skalowanie do wielu próbek       | 60 min         |
| [Część 3: Wspólne wywoływanie w kohorcie](./03_joint_calling.md)               | Dodawanie wielopróbkowego wspólnego genotypowania przy użyciu operatorów kanałów do agregowania wyników z próbek | 45 min         |

Pod koniec tego kursu będziesz w stanie zastosować podstawowe koncepcje i narzędzia Nextflow'a do typowego przypadku użycia z genomiki.

Gotowy, aby rozpocząć kurs?

[Zacznij :material-arrow-right:](00_orientation.md){ .md-button .md-button--primary }

<!-- Clearfix for float -->
<div style="content: ''; clear: both; display: table;"></div>
