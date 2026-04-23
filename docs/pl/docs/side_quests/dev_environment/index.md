# Środowisko programistyczne

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Tłumaczenie wspomagane przez AI - [dowiedz się więcej i zasugeruj ulepszenia](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

Nowoczesne zintegrowane środowiska programistyczne (IDE) mogą diametralnie zmienić Twoje doświadczenia z rozwijaniem Nextflow'a. Ten dodatkowy moduł skupia się na wykorzystaniu VS Code i jego rozszerzenia Nextflow do szybszego pisania kodu, wczesnego wykrywania błędów oraz sprawnego poruszania się po złożonych workflow'ach.

!!! note "To nie jest tradycyjny samouczek"

    W odróżnieniu od innych modułów szkoleniowych, ten przewodnik jest zorganizowany jako zbiór szybkich wskazówek, porad i praktycznych przykładów, a nie jako samouczek krok po kroku. Każdą sekcję można eksplorować niezależnie, w zależności od Twoich zainteresowań i aktualnych potrzeb. Możesz swobodnie przeskakiwać między sekcjami i skupiać się na funkcjach, które będą najbardziej przydatne w Twoim procesie tworzenia workflow'ów.

## Co powinieneś wiedzieć wcześniej

Ten przewodnik zakłada, że ukończyłeś kurs [Hello Nextflow](../hello_nextflow/) i swobodnie poruszasz się po podstawowych koncepcjach Nextflow, w tym:

- **Podstawowa struktura workflow'u**: Rozumienie procesów, workflow'ów i sposobu ich łączenia
- **Operacje na kanałach**: Tworzenie kanałów, przekazywanie danych między procesami i używanie podstawowych operatorów
- **Moduły i organizacja**: Tworzenie wielokrotnego użytku modułów i używanie instrukcji include
- **Podstawy konfiguracji**: Używanie `nextflow.config` do parametrów, dyrektyw procesów i profili

## Czego się tu nauczysz

Ten przewodnik skupia się na **funkcjach produktywności IDE**, które uczynią Cię bardziej efektywnym programistą Nextflow:

- **Zaawansowane podświetlanie składni**: Rozumienie tego, co VS Code pokazuje Ci o strukturze kodu
- **Inteligentne autouzupełnianie**: Wykorzystanie kontekstowych podpowiedzi do szybszego pisania kodu
- **Wykrywanie błędów i diagnostyka**: Wychwytywanie błędów składniowych przed uruchomieniem workflow'u
- **Nawigacja po kodzie**: Szybkie przemieszczanie się między procesami, modułami i definicjami
- **Formatowanie i organizacja**: Utrzymywanie spójnego, czytelnego stylu kodu
- **Programowanie wspomagane przez AI** (opcjonalnie): Używanie nowoczesnych narzędzi AI zintegrowanych z IDE

!!! info "Dlaczego funkcje IDE dopiero teraz?"

    Prawdopodobnie używałeś już VS Code podczas kursu [Hello Nextflow](../hello_nextflow/), ale skupialiśmy się wtedy na nauce podstaw Nextflow, a nie na funkcjach IDE. Teraz, gdy swobodnie poruszasz się po podstawowych koncepcjach Nextflow, takich jak procesy, workflow'y, kanały i moduły, jesteś gotowy, aby w pełni wykorzystać zaawansowane funkcje IDE, które uczynią Cię bardziej efektywnym programistą.

    Potraktuj to jako „przejście na wyższy poziom" swojego środowiska programistycznego — ten sam edytor, którego używałeś, ma o wiele potężniejsze możliwości, które stają się naprawdę wartościowe, gdy rozumiesz już, w czym Ci pomagają.

---

## 0. Konfiguracja i rozgrzewka

Skonfigurujmy przestrzeń roboczą przeznaczoną do eksploracji funkcji IDE:

```bash title="Navigate to the IDE features directory"
cd side-quests/ide_features
```

Otwórz ten katalog w VS Code:

```bash title="Open VS Code in current directory"
code .
```

Katalog `ide_features` zawiera przykładowe workflow'y demonstrujące różne funkcje IDE:

```bash title="Show directory structure"
tree .
```

```console title="Project structure"
tree .
.
├── basic_workflow.nf
├── complex_workflow.nf
├── data
│   ├── sample_001.fastq.gz
│   ├── sample_002.fastq.gz
│   ├── sample_003.fastq.gz
│   ├── sample_004.fastq.gz
│   ├── sample_005.fastq.gz
│   └── sample_data.csv
├── modules
│   ├── fastqc.nf
│   ├── star.nf
│   └── utils.nf
└── nextflow.config

3 directories, 12 files
```

!!! note "O plikach przykładowych"

    - `basic_workflow.nf` to działający podstawowy workflow, który możesz uruchomić i modyfikować
    - `complex_workflow.nf` służy wyłącznie do celów ilustracyjnych i demonstrowania funkcji nawigacji — może nie uruchomić się poprawnie, ale pokazuje realistyczną strukturę wieloplikowego workflow'u

### Skróty klawiszowe

Niektóre funkcje opisane w tym przewodniku korzystają z opcjonalnych skrótów klawiszowych. Jeśli korzystasz z tego materiału przez GitHub Codespaces w przeglądarce, skróty mogą nie działać zgodnie z oczekiwaniami, ponieważ są używane do innych celów w Twoim systemie.

Jeśli uruchamiasz VS Code lokalnie — tak jak prawdopodobnie będziesz robić podczas faktycznego pisania workflow'ów — skróty będą działać zgodnie z opisem.

Jeśli używasz Maca, niektóre (nie wszystkie) skróty klawiszowe będą używać „cmd" zamiast „ctrl". Będziemy to zaznaczać w tekście jako `Ctrl/Cmd`.

### 0.1. Instalacja rozszerzenia Nextflow

!!! note "Używasz już devcontainerów?"

    Jeśli pracujesz w **GitHub Codespaces** lub używasz **lokalnego devcontainera**, rozszerzenie Nextflow jest prawdopodobnie już zainstalowane i skonfigurowane. Możesz pominąć poniższe kroki ręcznej instalacji i przejść bezpośrednio do eksploracji funkcji rozszerzenia.

Aby zainstalować rozszerzenie ręcznie:

1. Otwórz VS Code
2. Przejdź do widoku rozszerzeń, klikając ikonę rozszerzeń po lewej stronie: ![ikona rozszerzeń](img/extensions_icon.png) (skrót `Ctrl/Cmd+Shift+X`, jeśli uruchamiasz VS Code lokalnie)
3. Wyszukaj „Nextflow"
4. Zainstaluj oficjalne rozszerzenie Nextflow

![Instalacja rozszerzenia Nextflow](img/install_extension.png)

### 0.2. Układ przestrzeni roboczej

Ponieważ używałeś VS Code przez cały kurs Hello Nextflow, znasz już podstawy. Oto jak efektywnie zorganizować przestrzeń roboczą na potrzeby tej sesji:

- **Obszar edytora**: Do przeglądania i edytowania plików. Możesz podzielić go na wiele paneli, aby porównywać pliki obok siebie.
- **Eksplorator plików** (kliknij ![ikona eksploratora plików](img/files_icon.png)) (`Ctrl/Cmd+Shift+E`): Lokalne pliki i foldery w Twoim systemie. Trzymaj go otwartego po lewej stronie, aby nawigować między plikami.
- **Zintegrowany terminal** (`Ctrl+Shift+` backtick dla Windows i MacOS): Terminal do interakcji z komputerem na dole ekranu. Używaj go do uruchamiania Nextflow'a i innych poleceń.
- **Panel problemów** (`Ctrl+Shift+M`): VS Code wyświetla tu wykryte błędy i problemy. Przydatny do szybkiego przeglądu problemów.

Możesz przeciągać panele lub je ukrywać (`Ctrl/Cmd+B` przełącza pasek boczny), aby dostosować układ podczas pracy z przykładami.

### Podsumowanie

Masz skonfigurowany VS Code z rozszerzeniem Nextflow i rozumiesz układ przestrzeni roboczej dla efektywnego programowania.

### Co dalej?

Dowiedz się, jak podświetlanie składni pomaga zrozumieć strukturę kodu Nextflow na pierwszy rzut oka.

---

## 1. Podświetlanie składni i struktura kodu

Teraz, gdy Twoja przestrzeń robocza jest skonfigurowana, przyjrzyjmy się, jak podświetlanie składni w VS Code pomaga czytać i pisać kod Nextflow bardziej efektywnie.

### 1.1. Elementy składni Nextflow

Otwórz `basic_workflow.nf`, aby zobaczyć podświetlanie składni w akcji:

![Prezentacja składni](img/syntax_showcase.png)

Zwróć uwagę, jak VS Code wyróżnia:

- **Słowa kluczowe** (`process`, `workflow`, `input`, `output`, `script`) w wyraźnych kolorach
- **Literały łańcuchowe** i **parametry** z różnym stylem
- **Komentarze** w stonowanym kolorze
- **Zmienne** i **wywołania funkcji** z odpowiednim wyróżnieniem
- **Bloki kodu** z właściwymi prowadnicami wcięć

!!! note "Kolory zależne od motywu"

    Konkretne kolory, które widzisz, zależą od Twojego motywu VS Code (tryb ciemny/jasny), ustawień kolorów i wszelkich wprowadzonych przez Ciebie dostosowań. Ważne jest to, że różne elementy składni są wizualnie odróżnione od siebie, co ułatwia zrozumienie struktury kodu niezależnie od wybranego schematu kolorów.

### 1.2. Rozumienie struktury kodu

Podświetlanie składni pomaga szybko zidentyfikować:

- **Granice procesów**: Wyraźne rozróżnienie między różnymi procesami
- **Bloki wejścia/wyjścia**: Łatwe do zauważenia definicje przepływu danych
- **Bloki skryptów**: Faktyczne wykonywane polecenia
- **Operacje na kanałach**: Kroki transformacji danych
- **Dyrektywy konfiguracyjne**: Ustawienia specyficzne dla procesu

Ta wizualna organizacja staje się nieoceniona podczas pracy ze złożonymi workflow'ami zawierającymi wiele procesów i skomplikowane przepływy danych.

### Podsumowanie

Rozumiesz, jak podświetlanie składni w VS Code pomaga czytać strukturę kodu Nextflow i identyfikować różne elementy języka, co przyspiesza programowanie.

### Co dalej?

Dowiedz się, jak inteligentne autouzupełnianie przyspiesza pisanie kodu dzięki kontekstowym podpowiedziom.

---

## 2. Inteligentne autouzupełnianie

Funkcje autouzupełniania w VS Code pomagają pisać kod szybciej i z mniejszą liczbą błędów, sugerując odpowiednie opcje na podstawie kontekstu.

### 2.1. Kontekstowe podpowiedzi

Opcje autouzupełniania różnią się w zależności od miejsca w kodzie:

#### Operacje na kanałach

Otwórz ponownie `basic_workflow.nf` i spróbuj wpisać `channel.` w bloku workflow:

![Autouzupełnianie kanałów](img/autocomplete_channel.png)

Zobaczysz podpowiedzi dla:

- `fromPath()` — Utwórz kanał ze ścieżek plików
- `fromFilePairs()` — Utwórz kanał z par plików
- `of()` — Utwórz kanał z wartości
- `fromSRA()` — Utwórz kanał z akcesjów SRA
- I wiele więcej...

Pomaga to szybko znaleźć odpowiednią fabrykę kanałów bez konieczności pamiętania dokładnych nazw metod.

Możesz też odkrywać operatory dostępne dla kanałów. Na przykład wpisz `FASTQC.out.html.`, aby zobaczyć dostępne operacje:

![Autouzupełnianie operatorów kanałów](img/autocomplete_operators.png)

#### Dyrektywy procesów

Wewnątrz bloku skryptu procesu wpisz `task.`, aby zobaczyć dostępne właściwości środowiska uruchomieniowego:

![Autouzupełnianie właściwości zadania](img/autocomplete_task.png)

#### Konfiguracja

Otwórz `nextflow.config` i wpisz `process.` w dowolnym miejscu, aby zobaczyć dostępne dyrektywy procesów:

![Autouzupełnianie konfiguracji](img/autocomplete_config.png)

Zobaczysz podpowiedzi dla:

- `executor`
- `memory`
- `cpus`

Oszczędza to czas podczas konfigurowania procesów i działa w różnych zakresach konfiguracji. Na przykład spróbuj wpisać `docker.`, aby zobaczyć opcje konfiguracji specyficzne dla Docker'a.

### Podsumowanie

Możesz używać inteligentnego autouzupełniania VS Code do odkrywania dostępnych operacji na kanałach, dyrektyw procesów i opcji konfiguracji bez konieczności zapamiętywania składni.

### Co dalej?

Dowiedz się, jak wykrywanie błędów w czasie rzeczywistym pomaga wychwytywać problemy przed uruchomieniem workflow'u, już na etapie czytania kodu.

## 3. Wykrywanie błędów i diagnostyka

Wykrywanie błędów w czasie rzeczywistym w VS Code pomaga wychwytywać problemy przed uruchomieniem workflow'u.

### 3.1. Wykrywanie błędów składniowych

Stwórzmy celowy błąd, aby zobaczyć wykrywanie w akcji. Otwórz `basic_workflow.nf` i zmień nazwę procesu z `FASTQC` na `FASTQ` (lub inną nieprawidłową nazwę). VS Code natychmiast wyróżni błąd w bloku workflow czerwonym falistym podkreśleniem:

![Podkreślenie błędu](img/error_underline.png)

### 3.2. Panel problemów

Poza indywidualnym wyróżnianiem błędów, VS Code udostępnia scentralizowany panel problemów, który agreguje wszystkie błędy, ostrzeżenia i komunikaty informacyjne z całej przestrzeni roboczej. Otwórz go za pomocą `Ctrl/Cmd+Shift+M` i użyj ikony filtra, aby wyświetlić tylko błędy dotyczące bieżącego pliku:

![Filtrowanie panelu problemów](img/active_file.png)

Kliknij dowolny problem, aby przejść bezpośrednio do problematycznej linii.

![Panel problemów](img/problems_panel.png)

Napraw błąd, zmieniając nazwę procesu z powrotem na `FASTQC`.

### 3.3. Typowe wzorce błędów

Typowe błędy w składni Nextflow obejmują:

- **Brakujące nawiasy klamrowe**: Niedopasowane `{` lub `}`
- **Niekompletne bloki**: Brakujące wymagane sekcje w procesach
- **Nieprawidłowa składnia**: Zniekształcony DSL Nextflow
- **Literówki w słowach kluczowych**: Błędnie napisane dyrektywy procesów
- **Niezgodności kanałów**: Niezgodności typów

Serwer języka Nextflow wyróżnia te problemy w panelu problemów. Możesz sprawdzać je wcześnie, aby unikać błędów składniowych podczas uruchamiania pipeline'u.

### Podsumowanie

Możesz używać wykrywania błędów i panelu problemów w VS Code do wychwytywania błędów składniowych przed uruchomieniem workflow'u, oszczędzając czas i unikając frustracji.

### Co dalej?

Dowiedz się, jak efektywnie nawigować między procesami, modułami i definicjami w złożonych workflow'ach.

---

## 4. Nawigacja po kodzie i zarządzanie symbolami

Efektywna nawigacja jest kluczowa podczas pracy ze złożonymi workflow'ami obejmującymi wiele plików. Aby to zrozumieć, zastąp definicję procesu w `basic_workflow.nf` importem modułu, który dla Ciebie przygotowaliśmy:

=== "Po"

    ```groovy title="basic_workflow.nf" linenums="3"
    include { FASTQC } from './modules/fastqc.nf'
    ```

=== "Przed"

    ```groovy title="basic_workflow.nf" linenums="3"
    process FASTQC {
        tag "${sample_id}"
        publishDir "${params.output_dir}/fastqc", mode: 'copy'

        input:
        tuple val(sample_id), path(reads)

        output:
        tuple val(sample_id), path("*.html"), emit: html
        tuple val(sample_id), path("*.zip"), emit: zip

        script:
        def args = task.ext.args ?: ''
        """
        fastqc \\
            ${args} \\
            --threads ${task.cpus} \\
            ${reads}
        """
    }
    ```

### 4.1. Przejdź do definicji

Jeśli najedziesz myszką na nazwę procesu, np. `FASTQC`, zobaczysz wyskakujące okienko z interfejsem modułu (wejściami i wyjściami):

![Przejdź do definicji](img/syntax.png)

Ta funkcja jest szczególnie wartościowa podczas tworzenia workflow'ów, ponieważ pozwala zrozumieć interfejs modułu bez bezpośredniego otwierania pliku modułu.

Możesz szybko przejść do dowolnej definicji procesu, modułu lub zmiennej używając **Ctrl/Cmd+kliknięcie**. Najedź myszką na link do pliku modułu na górze skryptu i podążaj za linkiem zgodnie z sugestią:

![Podążaj za linkiem](img/follow_link.png)

To samo działa dla nazw procesów. Wróć do `basic_workflow.nf` i wypróbuj to na nazwie procesu `FASTQC` w bloku workflow. Przeniesie Cię bezpośrednio do nazwy procesu (która w tym przykładzie jest taka sama jak plik modułu, ale może znajdować się w połowie znacznie większego pliku).

Aby wrócić do poprzedniego miejsca, użyj **Alt+←** (lub **Ctrl+-** na Macu). To potężny sposób na eksplorację kodu bez gubienia miejsca, w którym się znajdowałeś.

Teraz przyjrzyjmy się nawigacji w bardziej złożonym workflow'u, używając `complex_workflow.nf` (wspomnianego wcześniej pliku wyłącznie ilustracyjnego). Ten workflow zawiera wiele procesów zdefiniowanych w osobnych plikach modułów, a także kilka zdefiniowanych bezpośrednio. Choć złożone struktury wieloplikowe mogą być trudne do ręcznego przeglądania, możliwość przeskakiwania do definicji znacznie ułatwia eksplorację.

1. Otwórz `complex_workflow.nf`
2. Przejdź do definicji modułów
3. Użyj **Alt+←** (lub **Ctrl+-**), aby wrócić
4. Przejdź do nazwy procesu `FASTQC` w bloku workflow. Przeniesie Cię bezpośrednio do nazwy procesu (która w tym przykładzie jest taka sama jak plik modułu, ale może znajdować się w połowie znacznie większego pliku).
5. Wróć ponownie
6. Przejdź do procesu `TRIM_GALORE` w bloku workflow. Jest on zdefiniowany bezpośrednio, więc nie przeniesie Cię do osobnego pliku, ale nadal pokaże Ci definicję procesu, a Ty będziesz mógł wrócić do poprzedniego miejsca.

### 4.2. Nawigacja po symbolach

Mając nadal otwarty `complex_workflow.nf`, możesz uzyskać przegląd wszystkich symboli w pliku, wpisując `@` w pasku wyszukiwania na górze VS Code (skrót klawiszowy to `Ctrl/Cmd+Shift+O`, ale może nie działać w Codespaces). Otwiera to panel nawigacji po symbolach, który wyświetla wszystkie symbole w bieżącym pliku:

![Nawigacja po symbolach](img/symbols.png)

Pokazuje to:

- Wszystkie definicje procesów
- Definicje workflow'ów (w tym pliku zdefiniowane są dwa workflow'y)
- Definicje funkcji

Zacznij pisać, aby filtrować wyniki.

### 4.3. Znajdź wszystkie odwołania

Wiedza o tym, gdzie proces lub zmienna jest używana w całej bazie kodu, może być bardzo pomocna. Na przykład, jeśli chcesz znaleźć wszystkie odwołania do procesu `FASTQC`, zacznij od przejścia do jego definicji. Możesz to zrobić, otwierając bezpośrednio `modules/fastqc.nf` lub używając funkcji szybkiej nawigacji VS Code z `Ctrl/Cmd+kliknięcie`, jak robiliśmy to wcześniej. Po przejściu do definicji procesu kliknij prawym przyciskiem myszy na nazwie procesu `FASTQC` i wybierz „Find All References" z menu kontekstowego, aby zobaczyć wszystkie miejsca, w których jest używany.

![Znajdź odwołania](img/references.png)

Ta funkcja wyświetla wszystkie miejsca, w których `FASTQC` jest przywoływany w Twojej przestrzeni roboczej, w tym jego użycie w dwóch odrębnych workflow'ach. Ta informacja jest kluczowa przy ocenie potencjalnego wpływu modyfikacji procesu `FASTQC`.

### 4.4. Panel konspektu

Panel konspektu, znajdujący się na pasku bocznym eksploratora (kliknij ![ikona eksploratora](img/files_icon.png)), zapewnia wygodny przegląd wszystkich symboli w bieżącym pliku. Funkcja ta pozwala szybko nawigować i zarządzać strukturą kodu, wyświetlając funkcje, zmienne i inne kluczowe elementy w widoku hierarchicznym.

![Panel konspektu](img/outline.png)

Używaj panelu konspektu do szybkiego przechodzenia do różnych części kodu bez korzystania z przeglądarki plików.

### 4.5. Wizualizacja DAG

Rozszerzenie Nextflow dla VS Code może wizualizować Twój workflow jako skierowany graf acykliczny (DAG). Pomaga to zrozumieć przepływ danych i zależności między procesami. Otwórz `complex_workflow.nf` i kliknij przycisk „Preview DAG" nad `workflow {` (drugi blok `workflow` w tym pliku):

![Podgląd DAG](img/dag_preview.png)

To jest tylko workflow „wejściowy", ale możesz też podejrzeć DAG dla wewnętrznych workflow'ów, klikając przycisk „Preview DAG" nad workflow'em `RNASEQ_PIPELINE {` wyżej w pliku:

![Podgląd DAG wewnętrznego workflow'u](img/dag_preview_inner.png)

W tym workflow'ie możesz używać węzłów w DAG do nawigacji do odpowiednich definicji procesów w kodzie. Kliknij węzeł, a przeniesie Cię do odpowiedniej definicji procesu w edytorze. Szczególnie gdy workflow rozrośnie się do dużych rozmiarów, może to naprawdę pomóc w nawigacji po kodzie i zrozumieniu, jak procesy są ze sobą połączone.

### Podsumowanie

Możesz efektywnie nawigować po złożonych workflow'ach, używając przejścia do definicji, wyszukiwania symboli, znajdowania odwołań i wizualizacji DAG, aby zrozumieć strukturę kodu i zależności.

### Co dalej?

Dowiedz się, jak efektywnie pracować z wieloma powiązanymi plikami w większych projektach Nextflow.

## 5. Praca z wieloma plikami

Prawdziwe programowanie w Nextflow wiąże się z pracą z wieloma powiązanymi plikami. Przyjrzyjmy się, jak VS Code pomaga efektywnie zarządzać złożonymi projektami.

### 5.1. Szybka nawigacja między plikami

Mając otwarty `complex_workflow.nf`, zauważysz, że importuje on kilka modułów. Przećwiczmy szybką nawigację między nimi.

Naciśnij **Ctrl+P** (lub **Cmd+P**) i zacznij wpisywać „fast":

VS Code pokaże Ci pasujące pliki. Wybierz `modules/fastqc.nf`, aby natychmiast tam przejść. Jest to znacznie szybsze niż klikanie przez eksplorator plików, gdy wiesz mniej więcej, jakiego pliku szukasz.

Wypróbuj to z innymi wzorcami:

- Wpisz „star", aby znaleźć plik modułu dopasowania STAR (`star.nf`)
- Wpisz „utils", aby znaleźć plik funkcji narzędziowych (`utils.nf`)
- Wpisz „config", aby przejść do plików konfiguracyjnych (`nextflow.config`)

### 5.2. Podzielony edytor do pracy z wieloma plikami

Podczas pracy z modułami często trzeba jednocześnie widzieć zarówno główny workflow, jak i definicje modułów. Skonfigurujmy to:

1. Otwórz `complex_workflow.nf`
2. Otwórz `modules/fastqc.nf` w nowej karcie
3. Kliknij prawym przyciskiem myszy na karcie `modules/fastqc.nf` i wybierz „Split Right"
4. Teraz możesz widzieć oba pliki obok siebie

![Podzielony edytor](img/split_editor.png)

Jest to nieocenione, gdy:

- Sprawdzasz interfejsy modułów podczas pisania wywołań workflow'u, a podgląd nie wystarcza
- Porównujesz podobne procesy w różnych modułach
- Debugujesz przepływ danych między workflow'em a modułami

### 5.3. Wyszukiwanie w całym projekcie

Czasami trzeba znaleźć, gdzie konkretne wzorce są używane w całym projekcie. Naciśnij `Ctrl/Cmd+Shift+F`, aby otworzyć panel wyszukiwania.

Spróbuj wyszukać `publishDir` w całej przestrzeni roboczej:

![Wyszukiwanie w projekcie](img/project_search.png)

Pokazuje to każdy plik, który używa katalogów publikowania, pomagając Ci:

- Zrozumieć wzorce organizacji wyjść
- Znaleźć przykłady konkretnych dyrektyw
- Zapewnić spójność między modułami

### Podsumowanie

Możesz zarządzać złożonymi projektami wieloplikowymi, używając szybkiej nawigacji między plikami, podzielonych edytorów i wyszukiwania w całym projekcie, aby efektywnie pracować z workflow'ami i modułami.

### Co dalej?

Dowiedz się, jak funkcje formatowania i utrzymania kodu utrzymują Twoje workflow'y zorganizowane i czytelne.

---

## 6. Formatowanie i utrzymanie kodu

Właściwe formatowanie kodu jest niezbędne nie tylko ze względów estetycznych, ale także dla poprawy czytelności, zrozumienia i łatwości aktualizacji złożonych workflow'ów.

### 6.1. Automatyczne formatowanie w akcji

Otwórz `basic_workflow.nf` i celowo zepsuj formatowanie:

- Usuń część wcięć: Zaznacz cały dokument i naciśnij `shift+tab` wiele razy, aby usunąć jak najwięcej wcięć.
- Dodaj dodatkowe spacje w losowych miejscach: w instrukcji `channel.fromPath` dodaj 30 spacji po `(`.
- Złam niektóre linie w niezręczny sposób: Dodaj nową linię między operatorem `.view {` a łańcuchem `Processing sample:`, ale nie dodawaj odpowiadającego nowego wiersza przed zamykającym nawiasem `}`.

Teraz naciśnij `Shift+Alt+F` (lub `Shift+Option+F` na MacOS), aby automatycznie sformatować:

VS Code natychmiast:

- Naprawia wcięcia, aby wyraźnie pokazać strukturę procesów
- Wyrównuje podobne elementy spójnie
- Usuwa zbędne białe znaki
- Utrzymuje czytelne podziały wierszy

Pamiętaj, że automatyczne formatowanie może nie rozwiązać każdego problemu ze stylem kodu. Serwer języka Nextflow stara się utrzymywać kod w porządku, ale szanuje też Twoje osobiste preferencje w pewnych obszarach. Na przykład, jeśli usuniesz wcięcia wewnątrz bloku `script` procesu, formater pozostawi je bez zmian, ponieważ możesz celowo preferować taki styl.

Obecnie nie ma ścisłego wymuszania stylu dla Nextflow, więc serwer języka oferuje pewną elastyczność. Będzie jednak konsekwentnie stosować reguły formatowania wokół definicji metod i funkcji, aby zachować przejrzystość.

### 6.2. Funkcje organizacji kodu

#### Szybkie komentowanie

Zaznacz blok kodu w swoim workflow'ie i naciśnij **Ctrl+/** (lub **Cmd+/**), aby go zakomentować:

```groovy
// workflow {
//     ch_input = channel.fromPath(params.input)
//         .splitCsv(header: true)
//         .map { row -> [row.sample_id, file(row.fastq_path)] }
//
//     FASTQC(ch_input)
// }
```

Jest to idealne do:

- Tymczasowego wyłączania części workflow'ów podczas programowania
- Dodawania wyjaśniających komentarzy do złożonych operacji na kanałach
- Dokumentowania sekcji workflow'u

Naciśnij ponownie **Ctrl+/** (lub **Cmd+/**), aby odkomentować kod.

#### Zwijanie kodu dla przeglądu

W `complex_workflow.nf` zwróć uwagę na małe strzałki obok definicji procesów. Kliknij je, aby zwinąć (zredukować) procesy:

![Zwijanie kodu](img/code_folding.png)

Daje to ogólny przegląd struktury workflow'u bez zagłębiania się w szczegóły implementacji.

#### Dopasowywanie nawiasów klamrowych

Umieść kursor obok dowolnego nawiasu `{` lub `}`, a VS Code wyróżni pasujący nawias. Użyj **Ctrl+Shift+\\** (lub **Cmd+Shift+\\**), aby przeskakiwać między pasującymi nawiasami.

Jest to kluczowe dla:

- Rozumienia granic procesów
- Znajdowania brakujących lub nadmiarowych nawiasów
- Nawigacji po zagnieżdżonych strukturach workflow'u

#### Wieloliniowe zaznaczanie i edycja

Do jednoczesnej edycji wielu linii VS Code oferuje zaawansowane możliwości wielokursorowe:

- **Wieloliniowe zaznaczanie**: Przytrzymaj **Ctrl+Alt** (lub **Cmd+Option** na MacOS) i używaj klawiszy strzałek, aby zaznaczyć wiele linii
- **Wieloliniowe wcięcia**: Zaznacz wiele linii i użyj **Tab**, aby dodać wcięcie, lub **Shift+Tab**, aby je usunąć z całych bloków

Jest to szczególnie przydatne do:

- Spójnego wcięcia całych bloków procesów
- Dodawania komentarzy do wielu linii jednocześnie
- Edytowania podobnych definicji parametrów w wielu procesach

### Podsumowanie

Możesz utrzymywać czysty, czytelny kod, używając automatycznego formatowania, funkcji komentowania, zwijania kodu, dopasowywania nawiasów i wieloliniowej edycji, aby efektywnie organizować złożone workflow'y.

### Co dalej?

Dowiedz się, jak VS Code integruje się z Twoim szerszym procesem programowania, wykraczając poza samą edycję kodu.

---

## 7. Integracja z procesem programowania

VS Code dobrze integruje się z Twoim procesem programowania, wykraczając poza samą edycję kodu.

### 7.1. Integracja z kontrolą wersji

!!! note "Codespaces i integracja z Git"

    Jeśli pracujesz w **GitHub Codespaces**, niektóre funkcje integracji z Git mogą nie działać zgodnie z oczekiwaniami, szczególnie skróty klawiszowe dla kontroli źródła. Możliwe też, że podczas początkowej konfiguracji odmówiłeś otwarcia katalogu jako repozytorium Git — to jest w porządku na potrzeby szkolenia.

Jeśli Twój projekt jest repozytorium git (jak w tym przypadku), VS Code pokazuje:

- Zmodyfikowane pliki z kolorowymi wskaźnikami
- Status Git na pasku stanu
- Widoki różnic inline
- Możliwości commitowania i pushowania

Otwórz panel kontroli źródła, używając przycisku kontroli źródła (![ikona kontroli źródła](img/source_control_icon.png)) (`Ctrl+Shift+G` lub `Cmd+Shift+G`, jeśli pracujesz z VS Code lokalnie), aby zobaczyć zmiany git i tworzyć commity bezpośrednio w edytorze.

![Panel kontroli źródła](img/source_control.png)

### 7.2. Uruchamianie i inspekcja workflow'ów

Uruchommy workflow, a następnie sprawdźmy wyniki. W zintegrowanym terminalu (`Ctrl+Shift+` backtick w Windows i MacOS) uruchom podstawowy workflow:

```bash title="Run the basic workflow"
nextflow run basic_workflow.nf --input data/sample_data.csv --output_dir results
```

Podczas działania workflow'u zobaczysz dane wyjściowe w czasie rzeczywistym w terminalu. Po zakończeniu możesz używać VS Code do inspekcji wyników bez opuszczania edytora:

1. **Przejdź do katalogów roboczych**: Użyj eksploratora plików lub terminala, aby przeglądać `.nextflow/work`
2. **Otwieraj pliki dziennika**: Kliknij na ścieżki plików dziennika w danych wyjściowych terminala, aby otworzyć je bezpośrednio w VS Code
3. **Sprawdzaj wyjścia**: Przeglądaj opublikowane katalogi wyników w eksploratorze plików
4. **Przeglądaj raporty wykonania**: Otwieraj raporty HTML bezpośrednio w VS Code lub przeglądarce

Dzięki temu wszystko jest w jednym miejscu, zamiast przełączać się między wieloma aplikacjami.

### Podsumowanie

Możesz zintegrować VS Code z kontrolą wersji i wykonywaniem workflow'ów, aby zarządzać całym procesem programowania z jednego interfejsu.

### Co dalej?

Zobacz, jak wszystkie te funkcje IDE współpracują ze sobą w codziennym procesie programowania.

---

## 8. Podsumowanie i szybkie notatki

Oto kilka szybkich notatek dotyczących każdej z omówionych powyżej funkcji IDE:

### 8.1. Rozpoczynanie nowej funkcji

1. **Szybkie otwieranie pliku** (`Ctrl+P` lub `Cmd+P`) w celu znalezienia odpowiednich istniejących modułów
2. **Podzielony edytor** do przeglądania podobnych procesów obok siebie
3. **Nawigacja po symbolach** (`Ctrl+Shift+O` lub `Cmd+Shift+O`) w celu zrozumienia struktury pliku
4. **Autouzupełnianie** do szybkiego pisania nowego kodu

### 8.2. Debugowanie problemów

1. **Panel problemów** (`Ctrl+Shift+M` lub `Cmd+Shift+M`) do jednoczesnego przeglądania wszystkich błędów
2. **Przejdź do definicji** (`Ctrl+kliknięcie` lub `Cmd+kliknięcie`) w celu zrozumienia interfejsów procesów
3. **Znajdź wszystkie odwołania**, aby zobaczyć, jak procesy są używane
4. **Wyszukiwanie w całym projekcie** w celu znalezienia podobnych wzorców lub problemów

### 8.3. Refaktoryzacja i ulepszanie

1. **Wyszukiwanie w całym projekcie** (`Ctrl+Shift+F` lub `Cmd+Shift+F`) w celu znalezienia wzorców
2. **Automatyczne formatowanie** (`Shift+Alt+F` lub `Shift+Option+F`) dla zachowania spójności
3. **Zwijanie kodu** w celu skupienia się na strukturze
4. **Integracja z Git** do śledzenia zmian

---

## Podsumowanie

Odbyłeś błyskawiczną wycieczkę po funkcjach IDE VS Code dla programowania w Nextflow. Narzędzia te znacząco zwiększą Twoją produktywność poprzez:

- **Redukcję błędów** dzięki sprawdzaniu składni w czasie rzeczywistym
- **Przyspieszenie programowania** dzięki inteligentnemu autouzupełnianiu
- **Poprawę nawigacji** w złożonych workflow'ach wieloplikowych
- **Utrzymanie jakości** dzięki spójnemu formatowaniu
- **Pogłębienie rozumienia** dzięki zaawansowanemu podświetlaniu i wizualizacji struktury

Nie oczekujemy, że zapamiętasz wszystko, ale teraz wiesz, że te funkcje istnieją i będziesz mógł je znaleźć, gdy będą Ci potrzebne. W miarę dalszego rozwijania workflow'ów Nextflow, te funkcje IDE staną się drugą naturą, pozwalając Ci skupić się na pisaniu wysokiej jakości kodu, zamiast zmagać się ze składnią i strukturą.

### Co dalej?

Zastosuj te umiejętności IDE podczas pracy z innymi modułami szkoleniowymi, na przykład:

- **[nf-test](nf-test.md)**: Twórz kompleksowe zestawy testów dla swoich workflow'ów
- **[Hello nf-core](../../hello_nf-core/)**: Buduj pipeline'y produkcyjnej jakości zgodne ze standardami społeczności

Prawdziwa moc tych funkcji IDE ujawnia się podczas pracy nad większymi, bardziej złożonymi projektami. Zacznij stopniowo włączać je do swojego procesu pracy — po kilku sesjach staną się drugą naturą i zmienią Twoje podejście do programowania w Nextflow.

Od wychwytywania błędów, zanim Cię spowolnią, po sprawne poruszanie się po złożonych bazach kodu — te narzędzia uczynią Cię pewniejszym i bardziej efektywnym programistą.

Miłego kodowania!
