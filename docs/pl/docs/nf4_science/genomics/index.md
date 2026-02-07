---
title: Nextflow dla Genomiki
hide:
  - toc
---

# Nextflow dla Genomiki

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Tłumaczenie wspomagane przez AI - [dowiedz się więcej i zasugeruj ulepszenia](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

Ten kurs szkoleniowy jest przeznaczony dla badaczy w dziedzinie genomiki i pokrewnych obszarach, którzy są zainteresowani tworzeniem lub dostosowywaniem pipeline'ów analizy danych.
Opiera się na [Hello Nextflow](../../hello_nextflow/) – szkoleniu dla początkujących i demonstruje, jak używać Nextflow w specyficznym kontekście dziedziny genomiki.

W szczególności, ten kurs pokazuje, jak zaimplementować prosty pipeline wykrywania wariantów z [GATK](https://gatk.broadinstitute.org/) (Genome Analysis Toolkit), szeroko stosowanym pakietem oprogramowania do analizy danych z sekwencjonowania wysokoprzepustowego.

Zaczynajmy! Kliknij przycisk "Open in GitHub Codespaces" poniżej, aby uruchomić środowisko szkoleniowe (najlepiej w oddzielnej karcie), a następnie czytaj dalej, podczas gdy się ładuje.

[![Open in GitHub Codespaces](https://github.com/codespaces/badge.svg)](https://codespaces.new/nextflow-io/training?quickstart=1&ref=master)

## Cele szkoleniowe

Przechodząc przez ten kurs, nauczysz się, jak stosować podstawowe koncepcje i narzędzia Nextflow do typowego przypadku użycia w genomice.

Pod koniec tego warsztatu będziesz w stanie:

- Napisać liniowy workflow do wykrywania wariantów dla pojedynczej próbki
- Odpowiednio obsługiwać pliki pomocnicze, takie jak pliki indeksów i zasoby genomu referencyjnego
- Wykorzystać paradygmat przepływu danych Nextflow do zrównoleglenia wykrywania wariantów dla poszczególnych próbek
- Zaimplementować wykrywanie wariantów dla wielu próbek, używając odpowiednich operatorów kanałów
- Zaimplementować testy pipeline'u dla poszczególnych kroków i całościowe, które odpowiednio obsługują specyficzne dla genomiki osobliwości

<!-- TODO for future expansion: add metadata/samplesheet handling -->

## Wymagania wstępne

Kurs zakłada minimalną znajomość następujących zagadnień:

- Narzędzia i formaty plików powszechnie używane w tej dziedzinie naukowej
- Doświadczenie z wierszem poleceń
- Podstawowe koncepcje i narzędzia Nextflow omówione w [Hello Nextflow](../../hello_nextflow/) – szkoleniu dla początkujących

Aby zapoznać się z wymaganiami technicznymi i konfiguracją środowiska, zobacz mini-kurs [Konfiguracja Środowiska](../../envsetup/).
