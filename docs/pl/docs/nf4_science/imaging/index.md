# Nextflow run for Imaging

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Tłumaczenie wspomagane przez AI - [dowiedz się więcej i zaproponuj poprawki](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

Ten kurs szkoleniowy jest przeznaczony dla badaczy zajmujących się obrazowaniem i biologią przestrzenną, którzy chcą uruchamiać i dostosowywać pipeline'y analizy danych.
Uczy podstawowych konceptów Nextflow związanych z uruchamianiem, organizowaniem i konfigurowaniem workflow'ów przy użyciu [nf-core/molkart](https://nf-co.re/molkart), pipeline'u do przetwarzania danych transkryptomiki przestrzennej Molecular Cartography.
Umiejętności, których się tutaj nauczysz, są przydatne w pracy z dowolnym pipeline'em Nextflow lub nf-core.

Zaczynajmy! Kliknij poniższy przycisk "Open in GitHub Codespaces", aby uruchomić środowisko szkoleniowe (najlepiej w osobnej karcie), a następnie czytaj dalej, podczas gdy się ładuje.

[![Open in GitHub Codespaces](https://github.com/codespaces/badge.svg)](https://codespaces.new/nextflow-io/training?quickstart=1&ref=master)

## Cele szkolenia

Przechodząc przez ten kurs, nauczysz się stosować podstawowe koncepty i narzędzia Nextflow do uruchamiania pipeline'ów analizy obrazów.

Pod koniec tego warsztatu będziesz w stanie:

- Uruchomić workflow Nextflow'a lokalnie i monitorować jego wykonanie
- Znaleźć i zinterpretować wyniki oraz pliki logów wygenerowane przez Nextflow'a
- Uruchomić pipeline nf-core z danymi testowymi i własnymi danymi wejściowymi
- Skonfigurować wykonanie pipeline'u przy użyciu profili i plików parametrów
- Zarządzać danymi wejściowymi przy użyciu arkuszy próbek i parametrów wiersza poleceń

## Odbiorcy i wymagania wstępne

Ten kurs zakłada minimalne zaznajomienie z następującymi zagadnieniami:

- Doświadczenie z wierszem poleceń
- Podstawowa znajomość formatów plików obrazowych (obrazy TIFF, dane tabelaryczne)

Wymagania techniczne i konfigurację środowiska znajdziesz w mini-kursie [Konfiguracja środowiska](../../envsetup/).
