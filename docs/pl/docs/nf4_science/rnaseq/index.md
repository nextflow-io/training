---
title: Nextflow dla RNAseq
hide:
  - toc
---

# Nextflow dla RNAseq

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Tłumaczenie wspomagane przez AI - [dowiedz się więcej i zasugeruj ulepszenia](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

Ten kurs szkoleniowy jest przeznaczony dla badaczy w dziedzinie transkryptomiki i pokrewnych obszarów, którzy są zainteresowani tworzeniem lub dostosowywaniem pipeline'ów analizy danych.
Opiera się on na szkoleniu dla początkujących [Hello Nextflow](../../hello_nextflow/) i demonstruje, jak używać Nextflow'a w konkretnym kontekście analizy masowego RNAseq.

W szczególności, ten kurs pokazuje, jak zaimplementować prosty pipeline przetwarzania masowego RNAseq do przycinania sekwencji adapterów, dopasowywania odczytów do genomu referencyjnego i przeprowadzania kontroli jakości (QC) na kilku etapach.

Zacznijmy! Kliknij przycisk "Open in GitHub Codespaces" poniżej, aby uruchomić środowisko szkoleniowe (najlepiej w osobnej karcie), a następnie kontynuuj czytanie podczas ładowania.

[![Open in GitHub Codespaces](https://github.com/codespaces/badge.svg)](https://codespaces.new/nextflow-io/training?quickstart=1&ref=master)

## Cele szkoleniowe

Przechodząc przez ten kurs, nauczysz się, jak zastosować podstawowe koncepcje i narzędzia Nextflow do typowego przypadku użycia RNAseq.

Pod koniec tego warsztatu będziesz w stanie:

- Napisać liniowy workflow do zastosowania podstawowych metod przetwarzania i QC dla RNAseq
- Odpowiednio obsługiwać pliki specyficzne dla dziedziny, takie jak FASTQ i zasoby genomu referencyjnego
- Obsługiwać dane sekwencjonowania single-end i paired-end
- Wykorzystać paradygmat przepływu danych Nextflow'a do zrównoleglenia przetwarzania RNAseq na poziomie próbek
- Agregować raporty QC z wielu etapów i próbek przy użyciu odpowiednich operatorów kanałów

<!-- TODO
- Configure pipeline execution and manage and optimize resource allocations
- Implement per-step and end-to-end pipeline tests that handle RNAseq-specific idiosyncrasies appropriately
-->
<!-- TODO for future expansion: add metadata/samplesheet handling -->

## Wymagania wstępne

Kurs zakłada minimalną znajomość następujących zagadnień:

- Narzędzia i formaty plików powszechnie używane w tej dziedzinie naukowej
- Doświadczenie z linią poleceń
- Podstawowe koncepcje i narzędzia Nextflow omówione w szkoleniu dla początkujących [Hello Nextflow](../../hello_nextflow/)

Aby zapoznać się z wymaganiami technicznymi i konfiguracją środowiska, zobacz mini-kurs [Environment Setup](../../envsetup/).
