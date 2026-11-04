---
title: Skalowanie z Seqera
hide:
  - toc
page_type: index_page
index_type: course
additional_information:
  technical_requirements: true
  learning_objectives:
    - Zarejestruj się na Seqera Platform i poznaj Community Showcase
    - Dodaj pipeline do przestrzeni roboczej i uruchom go z poziomu interfejsu webowego
    - Uwierzytelnij się i uruchamiaj pipeline'y z wiersza poleceń przy użyciu `tw` CLI
    - Zarejestruj pipeline hostowany na GitHub i uruchom go obydwoma sposobami
  audience_prerequisites:
    - "**Odbiorcy:** Ten kurs jest przeznaczony dla osób, które chcą uruchamiać pipeline'y Nextflow na dużą skalę przy użyciu Seqera Platform."
    - "**Umiejętności:** Zakładamy znajomość uruchamiania pipeline'ów nf-core z wiersza poleceń."
    - "**Kursy:** Wymagane ukończenie [Nextflow Run](../nextflow_run/index.md) i [Use nf-core](../nfcore_use/index.md) lub równoważne doświadczenie w uruchamianiu lokalnych pipeline'ów oraz pipeline'u `nf-core/rnaseq`."
---

# Skalowanie z Seqera

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Tłumaczenie wspomagane przez AI - [dowiedz się więcej i zasugeruj ulepszenia](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

**Skalowanie z Seqera to praktyczne wprowadzenie do uruchamiania i monitorowania pipeline'ów Nextflow przy użyciu Seqera Platform.**

Pracując na konkretnych przykładach, skonfigurujesz dostęp do Seqera Platform, uruchomisz pipeline produkcyjny zarówno z poziomu interfejsu webowego, jak i wiersza poleceń, a także dodasz nowy pipeline do swojej przestrzeni roboczej.

Po ukończeniu kursu będziesz mieć umiejętności i pewność siebie potrzebne do samodzielnego uruchamiania i monitorowania pipeline'ów na Seqera Platform.

<!-- additional_information -->

## Przegląd kursu

Kurs ma charakter praktyczny i opiera się na pipeline'ach, które uruchamiałeś/uruchamiałaś w ramach [Use nf-core](../nfcore_use/index.md).

Zaczniesz od rejestracji na Seqera Platform i uruchomienia `nf-core/rnaseq` — pipeline'u produkcyjnego — z poziomu interfejsu webowego.
Następnie przejdziesz do narzędzia wiersza poleceń `tw`, aby zrobić to samo z terminala, a na koniec zarejestrujesz nowy pipeline, `nf-core/demo`, i uruchomisz go obydwoma sposobami.

### Plan kursu

| Rozdział kursu                                                                     | Opis                                                                                                    | Szacowany czas |
| ---------------------------------------------------------------------------------- | ------------------------------------------------------------------------------------------------------- | -------------- |
| [Część 1: Uruchamianie pipeline'ów z interfejsu webowego](./01_run_with_seqera.md) | Skonfiguruj dostęp do Seqera Platform i uruchom pipeline produkcyjny z poziomu interfejsu webowego      | 20 min         |
| [Część 2: Uruchamianie pipeline'ów z wiersza poleceń](./02_launch_from_cli.md)     | Uwierzytelnij `tw` CLI, uruchom zapisany pipeline i zarejestruj nowy pipeline z poziomu wiersza poleceń | 25 min         |

Po ukończeniu kursu będziesz swobodnie uruchamiać i monitorować pipeline'y Nextflow na Seqera Platform — niezależnie od tego, czy wolisz pracować z interfejsu webowego, czy z wiersza poleceń.

Gotowy/Gotowa, żeby zacząć?

[Rozpocznij naukę :material-arrow-right:](00_orientation.md){ .md-button .md-button--primary }
