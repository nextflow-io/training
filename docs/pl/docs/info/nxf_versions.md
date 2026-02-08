---
title: Wersje Nextflow
description: Zrozumienie i zarządzanie ewolucją wersji składni Nextflow
hide:
  - toc
  - footer
---

## Aktualnie obsługiwana wersja składni Nextflow i wymagania

Od wersji 3.0 portalu szkoleniowego wszystkie nasze kursy opierają się na wersji 25.10.2 Nextflow, chyba że na stronie indeksu kursu podano inaczej (z wyjątkiem przestarzałych lub zarchiwizowanych materiałów, które mogą nie zawierać informacji o wersji).

Ponieważ kursy używają teraz typowanych wejść na poziomie workflow oraz dyrektyw wyjściowych na poziomie workflow, wymagają użycia parsera składni V2.
Jeśli zamierzasz korzystać ze środowiska, które udostępniamy poprzez [Github Codespaces](../envsetup/01_setup.md) lub [lokalne devcontainers](../envsetup/03_devcontainer.md), nie musisz nic robić, chyba że w instrukcjach kursu zaznaczono inaczej.
Jednak jeśli wolisz pracować przez szkolenia we własnym środowisku ([Instalacja ręczna](../envsetup/02_local.md)), musisz upewnić się, że używasz Nextflow w wersji 25.10.2 lub nowszej z włączonym parserem składni v2.

## Starsze wersje materiałów szkoleniowych

Nasze materiały szkoleniowe są wersjonowane od lutego 2025 roku.

Możesz uzyskać dostęp do starszych wersji materiałów szkoleniowych, które działają z wersjami Nextflow **przed 25.10.2** poprzez menu rozwijane na górze każdej strony, które pokazuje numerowaną wersję materiałów szkoleniowych.
Gdy wybierzesz starszą wersję materiałów szkoleniowych, linki do środowiska szkoleniowego automatycznie określą odpowiednią wersję środowiska.

## Inne informacje o wersjach składni Nextflow

Nextflow ma dwie odrębne koncepcje wersjonowania, które są czasami mylone: **warianty DSL** i **warianty parsera składni**.

**DSL1 vs DSL2** odnosi się do zasadniczo różnych sposobów pisania pipeline'ów Nextflow.
DSL1 był oryginalną składnią, gdzie procesy były niejawnie połączone przez kanały.
DSL2, wprowadzony w Nextflow 20.07, dodał funkcje modularności: możliwość importowania procesów i workflow z innych plików, jawne bloki `workflow` oraz nazwane wyjścia procesów.
DSL1 został oznaczony jako przestarzały w Nextflow 22.03 i usunięty w 22.12.
Cały nowoczesny kod Nextflow używa DSL2.

**Parser składni v1 vs v2** odnosi się do różnych parserów, które oba działają z kodem DSL2.
Parser v1 to oryginalny, bardziej tolerancyjny parser.
Parser v2 jest bardziej rygorystyczny i umożliwia nowe funkcje języka, takie jak statyczne typowanie (typowane wejścia i wyjścia) oraz dyrektywy wyjściowe na poziomie workflow.
Parser v2 zapewnia również lepsze komunikaty o błędach i wykrywa więcej błędów podczas parsowania, a nie w czasie wykonywania.
Parser v2 stanie się domyślny w Nextflow 26.04.

Podsumowując: DSL2 to język, który piszesz. Wariant parsera składni określa, jak ściśle ten język jest interpretowany i jakie zaawansowane funkcje są dostępne.

### Sprawdzanie i ustawianie wersji Nextflow

Możesz sprawdzić, jaka wersja Nextflow jest zainstalowana w Twoim systemie za pomocą polecenia `nextflow --version`.

Więcej informacji o tym, jak zaktualizować wersję Nextflow, znajdziesz w dokumentacji referencyjnej na temat [Aktualizacji Nextflow](https://www.nextflow.io/docs/latest/updating-nextflow.html).

### Włączanie parsera składni v2

Aby **włączyć** parser składni v2 dla bieżącej sesji, uruchom następujące polecenie w terminalu:

```bash
export NXF_SYNTAX_PARSER=v2
```

Aby ustawić to na stałe (do czasu, gdy v2 stanie się domyślny w Nextflow 26.04), dodaj polecenie export do profilu Twojej powłoki (`~/.bashrc`, `~/.zshrc`, itp.):

```bash
echo 'export NXF_SYNTAX_PARSER=v2' >> ~/.bashrc
source ~/.bashrc
```

Zauważ, że zmienna środowiskowa `NXF_SYNTAX_PARSER=v2` jest tymczasowym wymogiem.
Od Nextflow 26.04 parser v2 stanie się domyślny i to ustawienie nie będzie już potrzebne.

### Wyłączanie parsera składni v2

Aby **wyłączyć** parser składni v2 dla bieżącej sesji, uruchom następujące polecenie w terminalu:

```bash
export NXF_SYNTAX_PARSER=v1
```

<!-- Will it be possible to disable it in versions after 26.04? -->

### Migracja istniejącego kodu

Wskazówki dotyczące migracji istniejącego kodu do zgodności z nowszymi wersjami Nextflow znajdziesz w [Notatkach migracyjnych](https://www.nextflow.io/docs/latest/migrations/index.html) w dokumentacji referencyjnej.

Te dwa artykuły są szczególnie pomocne przy migracji do najnowszej wersji:

- [Migracja do wyjść workflow](https://www.nextflow.io/docs/latest/tutorials/workflow-outputs.html)
- [Migracja do typów statycznych](https://www.nextflow.io/docs/latest/tutorials/static-types.html)

Obie te funkcje są omówione w ramach szkolenia dla początkujących, począwszy od wersji 3.0 materiałów szkoleniowych.

W zależności od generacji kodu Nextflow, który zamierzasz migrować, możesz wykonać większość pracy za pomocą lintera Nextflow, używając polecenia `nextflow lint -format`.
Zobacz dokumentację CLI dla [`lint`](https://www.nextflow.io/docs/latest/reference/cli.html#lint), aby uzyskać więcej szczegółów.

Mamy nadzieję, że będzie to pomocne.
Jeśli potrzebujesz pomocy, skontaktuj się na Slacku lub na forum.

---

<div markdown class="homepage_logos">

![Seqera](../assets/img/seqera_logo.png#only-light)

![Seqera](../assets/img/seqera_logo_dark.png#only-dark)

</div>
