---
title: Execution Config
hide:
  - toc
page_type: index_page
index_type: course
additional_information:
  technical_requirements: true
  learning_objectives:
    - Przełączanie technologii pakowania oprogramowania między Docker a Conda
    - Wybór platformy wykonawczej i zrozumienie, jak Nextflow dostosowuje wykonanie zadań do jej specyfiki
    - Kontrolowanie przydziału zasobów obliczeniowych i automatyczne ponawianie zadań zakończonych niepowodzeniem
    - Definiowanie i łączenie profili w celu przełączania między gotowymi konfiguracjami
  audience_prerequisites:
    - "**Odbiorcy:** Ten kurs jest przeznaczony dla osób, które już wiedzą, jak uruchamiać lokalne pipeline'y Nextflow i chcą głębiej poznać konfigurację wykonania."
    - "**Umiejętności:** Zakładana jest podstawowa znajomość wiersza poleceń."
    - "**Kursy:** Wymagane ukończenie kursu [Nextflow Run](../nextflow_run/index.md) lub równoważna znajomość uruchamiania lokalnego pipeline'u przy użyciu `nextflow run`."
---

# Configure Execution

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Tłumaczenie wspomagane przez AI - [dowiedz się więcej i zasugeruj ulepszenia](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>


**Configure Execution to praktyczne wprowadzenie do dostosowywania wykonania pipeline'ów Nextflow do różnych środowisk obliczeniowych.**

Pracując nad ćwiczeniami zorientowanymi na konkretne cele, nauczysz się przełączać technologie pakowania oprogramowania, wybierać platformę wykonawczą, kontrolować przydział zasobów obliczeniowych i ponawianie zadań, a także grupować konfigurację w przełączalne profile.

Po ukończeniu kursu będziesz mieć umiejętności i pewność siebie, by konfigurować wykonanie pipeline'ów Nextflow jak profesjonalista.

<!-- additional_information -->

## Przegląd kursu

Ten kurs ma charakter praktyczny i opiera się na umiejętnościach zdobytych w ramach kursu [Nextflow Run](../nextflow_run/index.md).

Weźmiesz ten sam wieloetapowy pipeline z tamtego kursu i będziesz stopniowo dostosowywać jego konfigurację do różnych środowisk obliczeniowych, a następnie zgrupujesz wszystko w profile, między którymi możesz przełączać się w czasie wykonania.

### Plan zajęć

| Rozdział kursu                                                                          | Opis                                                                                         | Szacowany czas |
| --------------------------------------------------------------------------------------- | -------------------------------------------------------------------------------------------- | -------------- |
| [Część 1: Dostosowanie do środowiska obliczeniowego](./01_packaging_and_execution.md)   | Przełączanie technologii pakowania oprogramowania i wybór platformy wykonawczej              | 20 min         |
| [Część 2: Zarządzanie zasobami obliczeniowymi i błędami](./02_resources_and_retries.md) | Kontrolowanie przydziału zasobów i automatyczne ponawianie zadań zakończonych niepowodzeniem | 15 min         |
| [Część 3: Używanie profili do przełączania konfiguracji](./03_profiles.md)              | Definiowanie i łączenie profili oraz inspekcja w pełni rozwiązanej konfiguracji              | 15 min         |

Po ukończeniu tego kursu będziesz swobodnie konfigurować pipeline'y Nextflow dla różnych środowisk obliczeniowych i przełączać się między nimi bez zbędnych komplikacji.

Gotowy, żeby zacząć?

[Rozpocznij naukę :material-arrow-right:](00_orientation.md){ .md-button .md-button--primary }
