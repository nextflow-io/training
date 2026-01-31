---
title: Nextflow Run
hide:
  - toc
page_type: index_page
index_type: course
additional_information:
  technical_requirements: true
  learning_objectives:
    - Uruchamianie i zarządzanie wykonywaniem workflow Nextflow
    - Znajdowanie i interpretowanie wyników oraz plików dziennika
    - Rozpoznawanie podstawowych komponentów Nextflow w prostym wieloetapowym workflow
    - Konfigurowanie wykonywania pipeline na popularnych platformach obliczeniowych, w tym HPC i chmurze
    - Podsumowanie najlepszych praktyk dotyczących odtwarzalności, przenośności i ponownego wykorzystania kodu, które czynią pipeline zgodnymi z zasadami FAIR, w tym modułowość kodu i kontenery oprogramowania
  audience_prerequisites:
    - "**Odbiorcy:** Ten kurs jest przeznaczony dla osób, które są całkowicie nowe w Nextflow i chcą uruchamiać istniejące pipeline."
    - "**Umiejętności:** Zakładana jest pewna znajomość wiersza poleceń, podstawowych koncepcji skryptowania i popularnych formatów plików."
    - "**Dziedzina:** Wszystkie ćwiczenia są niezależne od domeny, więc nie jest wymagana wcześniejsza wiedza naukowa."
---

# Nextflow Run

**Nextflow Run to praktyczne wprowadzenie do uruchamiania odtwarzalnych i skalowalnych workflow analizy danych.**

Pracując z praktycznymi przykładami i ćwiczeniami z przewodnikiem, poznasz podstawy używania Nextflow, w tym jak wykonywać pipeline, zarządzać plikami i zależnościami oprogramowania, bez wysiłku równoleglić wykonywanie i uruchamiać workflow w różnych środowiskach obliczeniowych.

Zdobędziesz umiejętności i pewność siebie, aby rozpocząć uruchamianie workflow z Nextflow.

<!-- additional_information -->

## Przegląd kursu

### Co będziesz robić

Ten kurs jest praktyczny, z ćwiczeniami zorientowanymi na cele, zaprojektowanymi tak, aby stopniowo wprowadzać informacje.

Wykonasz kilka wersji pipeline Nextflow, który przetwarza tekstowe dane wejściowe.
Zaczniesz od prostej wersji składającej się z jednego kroku, a ostatecznie przejdziesz do wieloetapowej wersji, która pobiera plik CSV z tabelarycznymi danymi tekstowymi, wykonuje kilka kroków transformacji i generuje pojedynczy plik tekstowy zawierający obraz ASCII postaci mówiącej przekształcony tekst.

Ten kurs koncentruje się na uruchamianiu pipeline (nazwany od podstawowego polecenia `nextflow run`).
Jeśli szukasz wprowadzenia do tworzenia pipeline Nextflow, zobacz [Hello Nextflow](../hello_nextflow/index.md).

### Plan lekcji

Podzieliliśmy to na trzy części, z których każda koncentruje się na konkretnych aspektach uruchamiania i zarządzania pipeline napisanymi w Nextflow.

| Rozdział kursu                                                 | Podsumowanie                                                                                                                                    | Szacowany czas |
| -------------------------------------------------------------- | ----------------------------------------------------------------------------------------------------------------------------------------------- | -------------- |
| [Część 1: Podstawowe operacje](./01_basics.md)                 | Uruchamianie i zarządzanie wykonywaniem prostego workflow                                                                                       | 30 min         |
| [Część 2: Uruchamianie prawdziwych pipeline](./02_pipeline.md) | Przetwarzanie złożonych danych wejściowych, uruchamianie wieloetapowych workflow, używanie kontenerów i bezproblemowe równoleglenie wykonywania | 60 min         |
| [Część 3: Konfiguracja uruchamiania](./03_config.md)           | Dostosowywanie zachowania pipeline i optymalizacja użycia w różnych środowiskach obliczeniowych                                                 | 60 min         |

Pod koniec tego kursu będziesz dobrze przygotowany do podjęcia kolejnych kroków w swojej podróży ku uruchamianiu odtwarzalnych workflow dla potrzeb obliczeń naukowych.

Gotowy do rozpoczęcia kursu?

[Rozpocznij naukę :material-arrow-right:](00_orientation.md){ .md-button .md-button--primary }
