---
title: RNAseq
hide:
  - toc
page_type: index_page
index_type: course
additional_information:
  technical_requirements: true
  learning_objectives:
    - Napisać liniowy workflow do zastosowania podstawowych metod przetwarzania i QC dla RNAseq
    - Odpowiednio obsługiwać pliki specyficzne dla dziedziny, takie jak FASTQ i zasoby genomu referencyjnego
    - Obsługiwać dane sekwencjonowania single-end i paired-end
    - Wykorzystać paradygmat przepływu danych Nextflow'a do zrównoleglenia przetwarzania RNAseq na poziomie próbek
    - Agregować raporty QC z wielu etapów i próbek przy użyciu odpowiednich operatorów kanałów
  audience_prerequisites:
    - "**Odbiorcy:** Ten kurs jest przeznaczony dla badaczy w dziedzinie transkryptomiki i pokrewnych obszarów, którzy chcą tworzyć lub dostosowywać pipeline'y analizy danych."
    - "**Umiejętności:** Zakładamy pewną znajomość linii poleceń, podstawowych koncepcji skryptowania oraz powszechnych formatów plików RNAseq."
    - "**Wymagania wstępne:** Podstawowe koncepcje i narzędzia Nextflow omówione w [Hello Nextflow](../../hello_nextflow/)."
---

# Nextflow dla RNAseq

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Tłumaczenie wspomagane przez AI - [dowiedz się więcej i zasugeruj ulepszenia](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

**Praktyczny kurs stosowania Nextflow'a w rzeczywistym przypadku użycia z dziedziny transkryptomiki: przetwarzanie masowego RNAseq przy użyciu Trim Galore, HISAT2 i FastQC.**

Ten kurs opiera się na szkoleniu dla początkujących [Hello Nextflow](../../hello_nextflow/) i demonstruje, jak używać Nextflow'a w konkretnym kontekście analizy masowego RNAseq.
Zaimplementujesz pipeline przetwarzania, który przycina sekwencje adapterów, dopasowuje odczyty do genomu referencyjnego i przeprowadza kontrolę jakości (QC) na kilku etapach.

<!-- additional_information -->

## Przegląd kursu

Ten kurs jest praktyczny, z ćwiczeniami zorientowanymi na cel, ustrukturyzowanymi tak, aby stopniowo wprowadzać informacje.

Zaczniesz od ręcznego uruchamiania narzędzi przetwarzania w terminalu, aby zrozumieć metodologię, a następnie stopniowo zbudujesz pipeline Nextflow'a, który automatyzuje i skaluje analizę.

### Plan lekcji

Podzieliliśmy to na trzy części, z których każda koncentruje się na konkretnych aspektach stosowania Nextflow'a w przypadku użycia RNAseq.

| Rozdział kursu                                                             | Podsumowanie                                                                                                               | Szacowany czas |
| -------------------------------------------------------------------------- | -------------------------------------------------------------------------------------------------------------------------- | -------------- |
| [Część 1: Przegląd metody](./01_method.md)                                 | Zrozumienie metodologii przetwarzania RNAseq i ręczne uruchamianie narzędzi                                                | 30 min         |
| [Część 2: Implementacja dla pojedynczej próbki](./02_single-sample.md)     | Budowanie pipeline'u, który przycina, dopasowuje i przeprowadza QC pojedynczej próbki, a następnie skaluje do wielu próbek | 60 min         |
| [Część 3: Implementacja dla wielu próbek paired-end](./03_multi-sample.md) | Rozszerzenie pipeline'u o obsługę danych paired-end i agregację raportów QC z wielu próbek                                 | 45 min         |

Pod koniec tego kursu będziesz w stanie zastosować podstawowe koncepcje i narzędzia Nextflow do typowego przypadku użycia RNAseq.

Gotowy, aby rozpocząć kurs?

[Zacznij :material-arrow-right:](00_orientation.md){ .md-button .md-button--primary }

<!-- Clearfix for float -->
<div style="content: ''; clear: both; display: table;"></div>
