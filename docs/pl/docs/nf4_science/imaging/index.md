---
title: Nextflow run dla Obrazowania
hide:
  - toc
---

# Nextflow run dla Obrazowania

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Tłumaczenie wspomagane przez AI - [dowiedz się więcej i zasugeruj ulepszenia](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

Ten kurs szkoleniowy jest przeznaczony dla naukowców zajmujących się obrazowaniem i biologią przestrzenną, którzy są zainteresowani uruchamianiem i dostosowywaniem pipeline'ów analizy danych.
Kurs uczy podstawowych koncepcji Nextflow związanych z uruchamianiem, organizowaniem i konfigurowaniem workflow'ów przy użyciu [nf-core/molkart](https://nf-co.re/molkart), pipeline'u do przetwarzania danych transkryptomiki przestrzennej Molecular Cartography.
Umiejętności, których się tutaj nauczysz, można przenieść do dowolnego pipeline'u Nextflow lub nf-core.

Zaczynajmy! Kliknij przycisk "Open in GitHub Codespaces" poniżej, aby uruchomić środowisko szkoleniowe (najlepiej w osobnej karcie), a następnie czytaj dalej podczas jego ładowania.

[![Open in GitHub Codespaces](https://github.com/codespaces/badge.svg)](https://codespaces.new/nextflow-io/training?quickstart=1&ref=master)

## Cele szkolenia

Pracując nad tym kursem, nauczysz się, jak stosować podstawowe koncepcje i narzędzia Nextflow do uruchamiania pipeline'ów analizy obrazowania.

Pod koniec tego warsztatu będziesz w stanie:

- Uruchomić workflow Nextflow lokalnie i monitorować jego wykonywanie
- Znaleźć i zinterpretować wyjścia (wyniki) oraz pliki dziennika wygenerowane przez Nextflow
- Uruchomić pipeline nf-core z danymi testowymi oraz własnymi danymi wejściowymi
- Skonfigurować wykonywanie pipeline'u przy użyciu profili i plików parametrów
- Zarządzać danymi wejściowymi za pomocą arkuszy próbek i parametrów wiersza poleceń

## Odbiorcy i wymagania wstępne

Ten kurs zakłada podstawową znajomość następujących zagadnień:

- Doświadczenie z wierszem poleceń
- Podstawowa znajomość formatów plików obrazowych (obrazy TIFF, dane tabelaryczne)

Wymagania techniczne i konfigurację środowiska można znaleźć w mini-kursie [Environment Setup](../../envsetup/).
