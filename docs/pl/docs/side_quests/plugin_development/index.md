---
title: Tworzenie wtyczek
hide:
  - toc
---

# Tworzenie wtyczek

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Tłumaczenie wspomagane przez AI - [dowiedz się więcej i zasugeruj ulepszenia](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

System wtyczek Nextflow'a pozwala rozszerzać język o własne funkcje, hooki monitorujące, backendy wykonawcze i wiele więcej.
Wtyczki umożliwiają społeczności dodawanie funkcji do Nextflow'a bez modyfikowania jego rdzenia, co czyni je idealnym rozwiązaniem do współdzielenia wielokrotnie używanej funkcjonalności między pipeline'ami.

Podczas tego szkolenia nauczysz się korzystać z istniejących wtyczek oraz opcjonalnie tworzyć własne.

## Odbiorcy i wymagania wstępne

Część 1 dotyczy korzystania z istniejących wtyczek i jest przydatna dla wszystkich użytkowników Nextflow'a.
Części 2–6 obejmują budowanie własnych wtyczek i wymagają znajomości kodu Groovy oraz narzędzi do budowania.
Wcześniejsza znajomość Javy ani Groovy nie jest wymagana.

**Wymagania wstępne**

- Konto GitHub LUB lokalna instalacja opisana [tutaj](../../envsetup/02_local).
- Ukończony kurs [Hello Nextflow](../../hello_nextflow/index.md) lub równoważny.
- Java 21 lub nowsza (dostępna w środowisku szkoleniowym; potrzebna tylko do Części 2–6).

**Katalog roboczy:** `side-quests/plugin_development`

## Cele szkolenia

Po ukończeniu tego szkolenia będziesz potrafić:

**Korzystanie z wtyczek (Część 1):**

- Instalować i konfigurować istniejące wtyczki w swoich workflow'ach
- Importować i używać funkcji dostarczanych przez wtyczki

**Tworzenie wtyczek (Części 2–6):**

- Tworzyć nowy projekt wtyczki przy użyciu wbudowanego generatora projektów Nextflow'a
- Implementować własne funkcje wywoływalne z poziomu workflow'ów
- Budować, testować i instalować wtyczkę lokalnie
- Monitorować zdarzenia workflow'u (np. zakończenie zadania, start/koniec pipeline'u) na potrzeby własnego logowania lub powiadomień
- Dodawać opcje konfiguracyjne, aby wtyczki były konfigurowalne
- Dystrybuować swoją wtyczkę

## Plan lekcji

#### Część 1: Podstawy wtyczek

Korzystanie z istniejących wtyczek w workflow'u Nextflow'a i konfigurowanie ich zachowania.

#### Część 2: Tworzenie projektu wtyczki

Generowanie nowego projektu wtyczki i analiza jego struktury.

#### Część 3: Własne funkcje

Implementowanie własnych funkcji, budowanie wtyczki i uruchamianie jej w workflow'u.

#### Część 4: Testowanie

Pisanie i uruchamianie testów jednostkowych przy użyciu frameworka Spock.

#### Część 5: Monitorowanie workflow'u

Reagowanie na zdarzenia, takie jak zakończenie zadania, w celu zbudowania licznika zadań.

#### Część 6: Konfiguracja i dystrybucja

Odczytywanie ustawień z `nextflow.config`, aby uczynić wtyczkę konfigurowalną, a następnie nauka jej udostępniania.

Gotowy, żeby rozpocząć kurs?

[Zacznij naukę :material-arrow-right:](01_plugin_basics.md){ .md-button .md-button--primary }
