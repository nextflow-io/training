```markdown
---
title: Nextflow run dla Obrazowania
hide:
  - toc
---

# Nextflow run dla Obrazowania

Ten kurs szkoleniowy jest przeznaczony dla naukowców zajmujących się obrazowaniem i biologią przestrzenną, którzy są zainteresowani uruchamianiem i dostosowywaniem potoków analizy danych.
Kurs uczy podstawowych koncepcji Nextflow związanych z uruchamianiem, organizowaniem i konfigurowaniem workflow przy użyciu [nf-core/molkart](https://nf-co.re/molkart), pipeline do przetwarzania danych transkryptomiki przestrzennej Molecular Cartography.
Umiejętności, których się tutaj nauczysz, można przenieść do dowolnego pipeline Nextflow lub nf-core.

Zaczynajmy! Kliknij przycisk "Open in GitHub Codespaces" poniżej, aby uruchomić środowisko szkoleniowe (najlepiej w osobnej karcie), a następnie czytaj dalej podczas jego ładowania.

[![Open in GitHub Codespaces](https://github.com/codespaces/badge.svg)](https://codespaces.new/nextflow-io/training?quickstart=1&ref=master)

## Cele szkolenia

Pracując nad tym kursem, nauczysz się, jak stosować podstawowe koncepcje i narzędzia Nextflow do uruchamiania potoków analizy obrazowania.

Pod koniec tego warsztatu będziesz w stanie:

- Uruchomić workflow Nextflow lokalnie i monitorować jego wykonywanie
- Znaleźć i zinterpretować wyjścia (wyniki) i pliki dziennika wygenerowane przez Nextflow
- Uruchomić pipeline nf-core z danymi testowymi i własnymi danymi wejściowymi
- Skonfigurować wykonywanie pipeline przy użyciu profili i plików parametrów
- Zarządzać danymi wejściowymi za pomocą arkuszy próbek i parametrów wiersza poleceń

## Odbiorcy i wymagania wstępne

Ten kurs zakłada podstawową znajomość następujących zagadnień:

- Doświadczenie z wierszem poleceń
- Podstawowa znajomość formatów plików obrazowych (obrazy TIFF, dane tabelaryczne)

Wymagania techniczne i konfigurację środowiska można znaleźć w mini-kursie [Environment Setup](../../envsetup/).
```
