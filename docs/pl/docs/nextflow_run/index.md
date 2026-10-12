---
title: Nextflow Run
hide:
  - toc
page_type: index_page
index_type: course
additional_information:
  technical_requirements: true
  learning_objectives:
    - Uruchamianie pipeline'ów Nextflow z wiersza poleceń i zarządzanie nimi
    - Zrozumienie, w jaki sposób kanały i operatory umożliwiają wydajne workflow'y z wieloma wejściami i wieloma krokami
    - Używanie kontenerów do zarządzania zależnościami oprogramowania i zapewnienia odtwarzalności
    - Konfigurowanie wykonywania pipeline'u i jego wyników
    - Generowanie raportów wykonania, przeglądanie historii poprzednich uruchomień i czyszczenie starych katalogów roboczych
    - Uruchamianie pipeline'ów bezpośrednio ze zdalnych repozytoriów, takich jak GitHub
  audience_prerequisites:
    - "**Odbiorcy:** Ten kurs jest przeznaczony dla osób, które są całkowicie nowe w Nextflow i chcą uruchamiać istniejące pipeline'y."
    - "**Umiejętności:** Zakładana jest pewna znajomość wiersza poleceń, podstawowych koncepcji skryptowania i popularnych formatów plików."
    - "**Dziedzina:** Wszystkie ćwiczenia są niezależne od dziedziny, więc nie jest wymagana wcześniejsza wiedza naukowa."
---

# Nextflow Run

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Tłumaczenie wspomagane przez AI - [dowiedz się więcej i zasugeruj ulepszenia](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

**Nextflow Run to praktyczne wprowadzenie do uruchamiania odtwarzalnych i skalowalnych workflow'ów analizy danych.**

Pracując przez serię ćwiczeń zorientowanych na cele, poznasz podstawy uruchamiania pipeline'ów Nextflow i zarządzania nimi, zrozumiesz, jak kanały i operatory umożliwiają równoległe przetwarzanie wielu wejść, oraz nauczysz się używać kontenerów do zarządzania zależnościami oprogramowania.

Zdobędziesz umiejętności i pewność siebie, aby rozpocząć uruchamianie workflow'ów z Nextflow'em.

<!-- additional_information -->

## Przegląd kursu

Ten kurs jest praktyczny, z ćwiczeniami zorientowanymi na cele, zaprojektowanymi tak, aby stopniowo wprowadzać informacje.

Wykonasz kilka wariantów pipeline'u Nextflow przetwarzającego tekstowe dane wejściowe — zaczniesz od prostej wersji składającej się z jednego kroku, a następnie przejdziesz do wieloetapowej implementacji pobierającej plik CSV z danymi wejściowymi, wykonującej kilka kroków transformacji i generującej pojedynczy plik tekstowy zawierający grafikę ASCII wygenerowaną przez narzędzie uruchomione w kontenerze.

Ten kurs koncentruje się na uruchamianiu pipeline'ów (nazwany od podstawowego polecenia `nextflow run`).
Jeśli szukasz wprowadzenia do tworzenia pipeline'ów Nextflow, zobacz [Hello Nextflow](../hello_nextflow/index.md).

!!! note "Uwaga"

    Szukasz poprzedniej wersji tego kursu? Została zastąpiona przez wersję na tej stronie, ale nadal można ją przeglądać w [wydaniu 3.6.1](https://training.nextflow.io/3.6.1/nextflow_run/) serwisu szkoleniowego.

### Plan lekcji

| Rozdział kursu                                                              | Podsumowanie                                                                                                           | Szacowany czas |
| --------------------------------------------------------------------------- | ---------------------------------------------------------------------------------------------------------------------- | -------------- |
| [Część 1: Uruchamianie Nextflow'a](./01_run_nextflow.md)                    | Uruchamianie pipeline'ów Nextflow i zarządzanie nimi oraz zrozumienie podstawowych mechanizmów workflow'u              | 25 min         |
| [Część 2: Konfiguracja pipeline'u](./02_configure_pipeline.md)              | Konfigurowanie wykonywania pipeline'u i jego wyników przy użyciu `nextflow.config`                                     | 20 min         |
| [Część 3: Zarządzanie uruchomieniami workflow'u](./03_manage_executions.md) | Generowanie raportów wykonania, przeglądanie historii poprzednich uruchomień i czyszczenie starych katalogów roboczych | 10 min         |
| [Część 4: Uruchamianie zdalnych pipeline'ów](./04_remote_repositories.md)   | Uruchamianie pipeline'u bezpośrednio z GitHub i przypinanie go do konkretnej wersji                                    | 10 min         |

Pod koniec tego kursu będziesz dobrze przygotowany do podjęcia kolejnych kroków na swojej drodze ku uruchamianiu odtwarzalnych workflow'ów dla potrzeb obliczeń naukowych.

Gotowy do rozpoczęcia kursu?

[Rozpocznij naukę :material-arrow-right:](00_orientation.md){ .md-button .md-button--primary }
