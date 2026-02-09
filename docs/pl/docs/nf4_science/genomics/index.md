---
title: Nextflow dla Genomiki
hide:
  - toc
page_type: index_page
index_type: course
additional_information:
  technical_requirements: true
  learning_objectives:
    - Napisać liniowy workflow do wykrywania wariantów dla pojedynczej próbki
    - Odpowiednio obsługiwać pliki pomocnicze, takie jak pliki indeksów i zasoby genomu referencyjnego
    - Wykorzystać paradygmat przepływu danych Nextflow do zrównoleglenia wykrywania wariantów dla poszczególnych próbek
    - Zaimplementować wielopróbkowe wspólne genotypowanie, używając odpowiednich operatorów kanałów
  audience_prerequisites:
    - "**Odbiorcy:** Ten kurs jest przeznaczony dla badaczy w dziedzinie genomiki i pokrewnych obszarach, którzy chcą tworzyć lub dostosowywać pipeline'y analizy danych."
    - "**Umiejętności:** Zakładamy pewną znajomość wiersza poleceń, podstawowych koncepcji skryptowania oraz popularnych formatów plików genomicznych."
    - "**Wymagania wstępne:** Podstawowe koncepcje i narzędzia Nextflow omówione w [Hello Nextflow](../../hello_nextflow/)."
---

# Nextflow dla Genomiki

**Praktyczny kurs stosowania Nextflow do rzeczywistego przypadku użycia w genomice: wykrywania wariantów za pomocą GATK.**

Ten kurs opiera się na [Hello Nextflow](../../hello_nextflow/) – szkoleniu dla początkujących i demonstruje, jak używać Nextflow w specyficznym kontekście dziedziny genomiki.
Zaimplementujesz pipeline wykrywania wariantów z [GATK](https://gatk.broadinstitute.org/) (Genome Analysis Toolkit), szeroko stosowanym pakietem oprogramowania do analizy danych z sekwencjonowania wysokoprzepustowego.

<!-- additional_information -->

## Przegląd kursu

Kurs ma charakter praktyczny, z ćwiczeniami ukierunkowanymi na cel, zorganizowanymi tak, aby stopniowo wprowadzać nowe informacje.

Zaczniesz od ręcznego uruchamiania narzędzi do wykrywania wariantów w terminalu, aby zrozumieć metodologię, a następnie krok po kroku zbudujesz pipeline Nextflow automatyzujący i skalujący analizę.

### Plan lekcji

Podzieliliśmy kurs na trzy części, z których każda koncentruje się na konkretnych aspektach stosowania Nextflow do przypadku użycia w genomice.

| Rozdział kursu                                                                              | Podsumowanie                                                                                                    | Szacowany czas trwania |
| ------------------------------------------------------------------------------------------- | --------------------------------------------------------------------------------------------------------------- | ---------------------- |
| [Część 1: Przegląd metody](./01_method.md)                                                  | Zrozumienie metodologii wykrywania wariantów i ręczne uruchamianie narzędzi                                     | 30 min                 |
| [Część 2: Wykrywanie wariantów dla pojedynczych próbek](./02_per_sample_variant_calling.md) | Budowa pipeline'u indeksującego pliki BAM i wykrywającego warianty, a następnie jego skalowanie na wiele próbek | 60 min                 |
| [Część 3: Wspólne genotypowanie kohort](./03_joint_calling.md)                              | Dodanie wielopróbkowego wspólnego genotypowania przy użyciu operatorów kanałów do agregacji wyników             | 45 min                 |

Pod koniec tego kursu będziesz w stanie zastosować podstawowe koncepcje i narzędzia Nextflow do typowego przypadku użycia w genomice.

Gotowy, aby rozpocząć kurs?

[Rozpocznij :material-arrow-right:](00_orientation.md){ .md-button .md-button--primary }

<!-- Clearfix for float -->
<div style="content: ''; clear: both; display: table;"></div>
