---
title: Nextflow Run
hide:
  - toc
page_type: index_page
index_type: course
additional_information:
  technical_requirements: true
  learning_objectives:
    - Uruchamianie i zarządzanie wykonywaniem workflow'ów Nextflow
    - Znajdowanie i interpretowanie wyników oraz plików dziennika
    - Rozpoznawanie podstawowych komponentów Nextflow w prostym wieloetapowym workflow'ie
    - Konfigurowanie wykonywania pipeline'u na popularnych platformach obliczeniowych, w tym HPC i chmurze
    - Podsumowanie najlepszych praktyk dotyczących odtwarzalności, przenośności i ponownego wykorzystania kodu, które czynią pipeline'y zgodnymi z zasadami FAIR, w tym modułowość kodu i kontenery oprogramowania
  audience_prerequisites:
    - "**Odbiorcy:** Ten kurs jest przeznaczony dla osób, które są całkowicie nowe w Nextflow i chcą uruchamiać istniejące pipeline'y."
    - "**Umiejętności:** Zakładana jest pewna znajomość wiersza poleceń, podstawowych koncepcji skryptowania i popularnych formatów plików."
    - "**Dziedzina:** Wszystkie ćwiczenia są niezależne od dziedziny, więc nie jest wymagana wcześniejsza wiedza naukowa."
---

# Nextflow Run

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Tłumaczenie wspomagane przez AI - [dowiedz się więcej i zasugeruj ulepszenia](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

**Nextflow Run to praktyczne wprowadzenie do uruchamiania odtwarzalnych i skalowalnych workflow'ów analizy danych.**

Pracując z praktycznymi przykładami i ćwiczeniami z przewodnikiem, poznasz podstawy używania Nextflow'a, w tym jak wykonywać pipeline'y, zarządzać plikami i zależnościami oprogramowania, bez wysiłku równoleglić wykonywanie i uruchamiać workflow'y w różnych środowiskach obliczeniowych.

Zdobędziesz umiejętności i pewność siebie, aby rozpocząć uruchamianie workflow'ów z Nextflow'em.

<!-- additional_information -->

## Przegląd kursu

### Co będziesz robić

Ten kurs jest praktyczny, z ćwiczeniami zorientowanymi na cele, zaprojektowanymi tak, aby stopniowo wprowadzać informacje.

Wykonasz kilka wariantów pipeline'u Nextflow, który przetwarza tekst.
Zaczniesz od prostej wersji składającej się z jednego kroku, a ostatecznie przejdziesz do wieloetapowej implementacji pobierającej plik CSV z tabelarycznymi danymi wejściowymi, wykonującej kilka transformacji i generującej pojedynczy plik tekstowy zawierający obraz ASCII postaci mówiącej przekształcone słowa.

Ten kurs koncentruje się na uruchamianiu pipeline'ów (nazwany od podstawowego polecenia `nextflow run`).
Jeśli szukasz wprowadzenia do tworzenia pipeline'ów Nextflow, zobacz [Hello Nextflow](../hello_nextflow/index.md).

### Plan lekcji

Podzieliliśmy to na trzy części, z których każda koncentruje się na konkretnych aspektach uruchamiania i zarządzania pipeline'ami napisanymi w Nextflow'ie.

| Rozdział kursu                                                    | Podsumowanie                                                                                                                                 | Szacowany czas |
| ----------------------------------------------------------------- | -------------------------------------------------------------------------------------------------------------------------------------------- | -------------- |
| [Część 1: Podstawowe operacje](./01_basics.md)                    | Uruchamianie i zarządzanie wykonywaniem prostego workflow'u                                                                                  | 30 min         |
| [Część 2: Uruchamianie prawdziwych pipeline'ów](./02_pipeline.md) | Przetwarzanie złożonych danych wejściowych, uruchamianie wieloetapowych workflow'ów, używanie kontenerów i bezproblemowe równoleglenie pracy | 60 min         |
| [Część 3: Konfiguracja uruchamiania](./03_config.md)              | Dostosowywanie zachowania pipeline'u i optymalizacja użycia w różnych środowiskach obliczeniowych                                            | 60 min         |

Pod koniec tego kursu będziesz dobrze przygotowany do podjęcia kolejnych kroków na swojej drodze ku uruchamianiu odtwarzalnych workflow'ów dla potrzeb obliczeń naukowych.

Gotowy do rozpoczęcia kursu?

[Rozpocznij naukę :material-arrow-right:](00_orientation.md){ .md-button .md-button--primary }
