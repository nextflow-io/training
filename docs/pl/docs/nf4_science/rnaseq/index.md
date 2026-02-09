# Nextflow dla RNAseq

Ten kurs szkoleniowy jest przeznaczony dla badaczy w dziedzinie transkryptomiki i pokrewnych dziedzin, którzy są zainteresowani tworzeniem lub dostosowywaniem pipeline'ów do analizy danych.
Opiera się na szkoleniu dla początkujących [Hello Nextflow](../../hello_nextflow/) i pokazuje, jak używać Nextflow'a w konkretnym kontekście analizy masowego RNAseq.

W szczególności, ten kurs pokazuje, jak zaimplementować prosty pipeline do przetwarzania masowego RNAseq, który przycina sekwencje adapterów, dopasowuje odczyty do genomu referencyjnego i przeprowadza kontrolę jakości (QC) na kilku etapach.

Zaczynajmy! Kliknij poniższy przycisk „Open in GitHub Codespaces", aby uruchomić środowisko szkoleniowe (najlepiej w osobnej karcie), a następnie czytaj dalej, podczas gdy się ładuje.

[![Open in GitHub Codespaces](https://github.com/codespaces/badge.svg)](https://codespaces.new/nextflow-io/training?quickstart=1&ref=master)

## Cele szkolenia

Przechodząc przez ten kurs, nauczysz się, jak zastosować podstawowe koncepcje i narzędzia Nextflow do typowego przypadku użycia RNAseq.

Pod koniec tego warsztatu będziesz w stanie:

- Napisać liniowy workflow, aby zastosować podstawowe metody przetwarzania i kontroli jakości RNAseq
- Odpowiednio obsługiwać pliki specyficzne dla dziedziny, takie jak FASTQ i zasoby genomu referencyjnego
- Obsługiwać dane sekwencjonowania single-end i paired-end
- Wykorzystać paradygmat przepływu danych Nextflow'a do zrównoleglenia przetwarzania RNAseq dla poszczególnych próbek
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

Informacje dotyczące wymagań technicznych i konfiguracji środowiska znajdziesz w mini-kursie [Konfiguracja środowiska](../../envsetup/).
