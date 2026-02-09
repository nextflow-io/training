# Środowisko programistyczne

Nowoczesne zintegrowane środowiska programistyczne (IDE) mogą dramatycznie zmienić Twoje doświadczenie z tworzeniem kodu w Nextflow. Ten dodatkowy quest skupia się konkretnie na wykorzystaniu VS Code i jego rozszerzenia Nextflow do szybszego pisania kodu, wczesnego wykrywania błędów i efektywnej nawigacji po złożonych workflow'ach.

!!! note "To nie jest tradycyjny samouczek"

    W przeciwieństwie do innych modułów szkoleniowych, ten przewodnik jest zorganizowany jako zbiór szybkich wskazówek, porad i praktycznych przykładów, a nie jako samouczek krok po kroku. Każda sekcja może być eksplorowana niezależnie, w zależności od Twoich zainteresowań i bieżących potrzeb programistycznych. Możesz swobodnie przeskakiwać między sekcjami i skupić się na funkcjach, które będą najbardziej przydatne w Twoim procesie tworzenia workflow'ów.

## Co powinieneś już wiedzieć

Ten przewodnik zakłada, że ukończyłeś kurs szkoleniowy [Hello Nextflow](../hello_nextflow/) i czujesz się komfortowo z podstawowymi koncepcjami Nextflow, w tym:

- **Podstawowa struktura workflow'a**: Rozumienie procesów, workflow'ów i sposobu ich łączenia
- **Operacje na kanałach**: Tworzenie kanałów, przekazywanie danych między procesami i używanie podstawowych operatorów
- **Moduły i organizacja**: Tworzenie modułów wielokrotnego użytku i używanie instrukcji include
- **Podstawy konfiguracji**: Używanie `nextflow.config` do parametrów, dyrektyw procesów i profili

## Czego się tutaj nauczysz

Ten przewodnik skupia się na **funkcjach produktywności IDE**, które uczynią Cię bardziej efektywnym programistą Nextflow:

- **Zaawansowane podświetlanie składni**: Zrozumienie, co VS Code pokazuje Ci o strukturze Twojego kodu
- **Inteligentne autouzupełnianie**: Wykorzystanie kontekstowych sugestii do szybszego pisania kodu
- **Wykrywanie błędów i diagnostyka**: Wychwytywanie błędów składniowych przed uruchomieniem workflow'a
- **Nawigacja po kodzie**: Szybkie przemieszczanie się między procesami, modułami i definicjami
- **Formatowanie i organizacja**: Utrzymywanie spójnego, czytelnego stylu kodu
- **Programowanie wspomagane AI** (opcjonalnie): Używanie nowoczesnych narzędzi AI zintegrowanych z Twoim IDE

!!! info "Dlaczego funkcje IDE teraz?"

    Prawdopodobnie już używałeś VS Code podczas kursu [Hello Nextflow](../hello_nextflow/), ale skupiliśmy się na nauce podstaw Nextflow, a nie funkcji IDE. Teraz, gdy czujesz się komfortowo z podstawowymi koncepcjami Nextflow, takimi jak procesy, workflow'y, kanały i moduły, jesteś gotowy wykorzystać zaawansowane funkcje IDE, które uczynią Cię bardziej efektywnym programistą.

    Pomyśl o tym jak o "awansowaniu" Twojego środowiska programistycznego - ten sam edytor, którego używałeś, ma znacznie potężniejsze możliwości, które stają się naprawdę wartościowe, gdy rozumiesz, w czym Ci pomagają.

---

## 0. Konfiguracja i rozgrzewka

Skonfigurujmy przestrzeń roboczą specjalnie do eksploracji funkcji IDE:

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
    - `complex_workflow.nf` jest przeznaczony wyłącznie do ilustracji, aby zademonstrować funkcje nawigacji - może nie działać poprawnie, ale pokazuje realistyczną strukturę wieloplikowego workflow'a

### Skróty klawiszowe

Niektóre funkcje w tym przewodniku będą używać opcjonalnych skrótów klawiszowych. Możliwe, że korzystasz z tego materiału przez GitHub Codespaces w przeglądarce, a w takim przypadku czasami skróty nie będą działać zgodnie z oczekiwaniami, ponieważ są używane do innych rzeczy w Twoim systemie.

Jeśli uruchamiasz VS Code lokalnie, tak jak prawdopodobnie będziesz robić podczas faktycznego pisania workflow'ów, skróty będą działać zgodnie z opisem.

Jeśli używasz Maca, niektóre (nie wszystkie) skróty klawiszowe będą używać "cmd" zamiast "ctrl", co oznaczymy w tekście jako `Ctrl/Cmd`.

### 0.1. Instalacja rozszerzenia Nextflow

!!! note "Już używasz Devcontainerów?"

    Jeśli pracujesz w **GitHub Codespaces** lub używasz **lokalnego devcontainera**, rozszerzenie Nextflow jest prawdopodobnie już zainstalowane i skonfigurowane. Możesz pominąć poniższe kroki ręcznej instalacji i przejść bezpośrednio do eksploracji funkcji rozszerzenia.

Aby zainstalować rozszerzenie ręcznie:

1. Otwórz VS Code
2. Przejdź do widoku rozszerzeń, klikając ikonę rozszerzeń po lewej stronie: ![ikona rozszerzeń](img/extensions_icon.png) (skrót `Ctrl/Cmd+Shift+X`, jeśli uruchamiasz VSCode lokalnie)
3. Wyszukaj "Nextflow"
4. Zainstaluj oficjalne rozszerzenie Nextflow

![Instalacja rozszerzenia Nextflow](img/install_extension.png)

### 0.2. Layout przestrzeni roboczej

Ponieważ używałeś VS Code podczas całego kursu Hello Nextflow, znasz już podstawy. Oto jak efektywnie zorganizować przestrzeń roboczą na tę sesję:

- **Obszar edytora**: Do przeglądania i edycji plików. Możesz podzielić go na wiele paneli, aby porównywać pliki obok siebie.
- **Eksplorator plików** kliknij (![ikona eksploratora plików](img/files_icon.png)) (`Ctrl/Cmd+Shift+E`): Lokalne pliki i foldery w Twoim systemie. Trzymaj to otwarte po lewej stronie, aby nawigować między plikami
- **Zintegrowany terminal** (`Ctrl+Shift+` backtick zarówno dla Windows, jak i MacOS): Terminal do interakcji z komputerem na dole. Użyj go do uruchamiania Nextflow'a lub innych poleceń.
- **Panel problemów** (`Ctrl+Shift+M`): VS Code pokaże tutaj wszelkie wykryte błędy i problemy. Jest to przydatne do szybkiego podświetlania problemów.

Możesz przeciągać panele lub je ukrywać (`Ctrl/Cmd+B`, aby przełączyć pasek boczny), aby dostosować swój layout podczas pracy z przykładami.

### Podsumowanie

Masz skonfigurowany VS Code z rozszerzeniem Nextflow i rozumiesz layout przestrzeni roboczej dla efektywnego programowania.

### Co dalej?

Dowiedz się, jak podświetlanie składni pomaga zrozumieć strukturę kodu Nextflow na pierwszy rzut oka.

---

## 1. Podświetlanie składni i struktura kodu

Teraz, gdy Twoja przestrzeń robocza jest skonfigurowana, zbadajmy, jak podświetlanie składni VS Code pomaga Ci efektywniej czytać i pisać kod Nextflow.

### 1.1. Elementy składni Nextflow

Otwórz `basic_workflow.nf`, aby zobaczyć podświetlanie składni w akcji:

![Prezentacja składni](img/syntax_showcase.png)

Zauważ, jak VS Code podświetla:

- **Słowa kluczowe** (`process`, `workflow`, `input`, `output`, `script`) w wyraźnych kolorach
- **Literały napisowe** i **parametry** z różnym stylem
- **Komentarze** w stonowanym kolorze
- **Zmienne** i **wywołania funkcji** z odpowiednim wyróżnieniem
- **Bloki kodu** z odpowiednimi prowadnicami wcięć

!!! note "Kolory zależne od motywu"

    Konkretne kolory, które widzisz, będą zależeć od Twojego motywu VS Code (tryb ciemny/jasny), ustawień kolorów i wszelkich dostosowań, które wykonałeś. Ważne jest, że różne elementy składni są wizualnie odróżnione od siebie, co ułatwia zrozumienie struktury kodu niezależnie od wybranego schematu kolorów.

### 1.2. Zrozumienie struktury kodu

Podświetlanie składni pomaga szybko zidentyfikować:

- **Granice procesów**: Wyraźne rozróżnienie między różnymi procesami
- **Bloki wejścia/wyjścia**: Łatwe do zauważenia definicje przepływu danych
- **Bloki skryptów**: Faktyczne polecenia, które są wykonywane
- **Operacje na kanałach**: Kroki transformacji danych
- **Dyrektywy konfiguracyjne**: Ustawienia specyficzne dla procesów

Ta wizualna organizacja staje się nieoceniona podczas pracy ze złożonymi workflow'ami zawierającymi wiele procesów i skomplikowane przepływy danych.

### Podsumowanie

Rozumiesz, jak podświetlanie składni VS Code pomaga czytać strukturę kodu Nextflow i identyfikować różne elementy języka dla szybszego programowania.

### Co dalej?

Dowiedz się, jak inteligentne autouzupełnianie przyspiesza pisanie kodu dzięki kontekstowym sugestiom.

---

## 2. Inteligentne autouzupełnianie

Funkcje autouzupełniania VS Code pomagają pisać kod szybciej i z mniejszą liczbą błędów, sugerując odpowiednie opcje w zależności od kontekstu.

### 2.1. Sugestie kontekstowe

Opcje autouzupełniania różnią się w zależności od tego, gdzie jesteś w swoim kodzie:

#### Operacje na kanałach

Otwórz ponownie `basic_workflow.nf` i spróbuj wpisać `channel.` w bloku workflow:

![Autouzupełnianie kanałów](img/autocomplete_channel.png)

Zobaczysz sugestie dla:

- `fromPath()` - Utwórz kanał ze ścieżek plików
- `fromFilePairs()` - Utwórz kanał z par plików
- `of()` - Utwórz kanał z wartości
- `fromSRA()` - Utwórz kanał z akcesji SRA
- I wiele więcej...

To pomaga szybko znaleźć odpowiednią fabrykę kanałów bez konieczności zapamiętywania dokładnych nazw metod.

Możesz również odkryć operatory dostępne do zastosowania na kanałach. Na przykład wpisz `FASTQC.out.html.`, aby zobaczyć dostępne operacje:

![Autouzupełnianie operatorów kanałów](img/autocomplete_operators.png)

#### Dyrektywy procesów

Wewnątrz bloku skryptu procesu wpisz `task.`, aby zobaczyć dostępne właściwości runtime:

![Autouzupełnianie właściwości zadań](img/autocomplete_task.png)

#### Konfiguracja

Otwórz nextflow.config i wpisz `process.` gdziekolwiek, aby zobaczyć dostępne dyrektywy procesów:

![Autouzupełnianie konfiguracji](img/autocomplete_config.png)

Zobaczysz sugestie dla:

- `executor`
- `memory`
- `cpus`

To oszczędza czas podczas konfigurowania procesów i działa w różnych zakresach konfiguracji. Na przykład spróbuj wpisać `docker.`, aby zobaczyć opcje konfiguracji specyficzne dla Dockera.

### Podsumowanie

Możesz używać inteligentnego autouzupełniania VS Code do odkrywania dostępnych operacji na kanałach, dyrektyw procesów i opcji konfiguracji bez zapamiętywania składni.

### Co dalej?

Dowiedz się, jak wykrywanie błędów w czasie rzeczywistym pomaga wychwytywać problemy przed uruchomieniem workflow'a, po prostu czytając kod.

## 3. Wykrywanie błędów i diagnostyka

Wykrywanie błędów w czasie rzeczywistym w VS Code pomaga wychwytywać problemy przed uruchomieniem workflow'a.

### 3.1. Wykrywanie błędów składniowych

Stwórzmy celowy błąd, aby zobaczyć wykrywanie w akcji. Otwórz `basic_workflow.nf` i zmień nazwę procesu z `FASTQC` na `FASTQ` (lub dowolną inną nieprawidłową nazwę). VS Code natychmiast podświetli błąd w bloku workflow czerwoną falistą linią:

![Podkreślenie błędu](img/error_underline.png)

### 3.2. Panel problemów

Poza indywidualnym podświetlaniem błędów, VS Code zapewnia scentralizowany panel problemów, który agreguje wszystkie błędy, ostrzeżenia i komunikaty informacyjne w całej przestrzeni roboczej. Otwórz go za pomocą `Ctrl/Cmd+Shift+M` i użyj ikony filtra, aby pokazać tylko błędy istotne dla bieżącego pliku:

![Filtrowanie panelu problemów](img/active_file.png)

Kliknij dowolny problem, aby przejść bezpośrednio do problematycznej linii

![Panel problemów](img/problems_panel.png)

Napraw błąd, zmieniając nazwę procesu z powrotem na `FASTQC`.

### 3.3. Typowe wzorce błędów

Typowe błędy w składni Nextflow obejmują:

- **Brakujące nawiasy**: Niedopasowane `{` lub `}`
- **Niekompletne bloki**: Brakujące wymagane sekcje w procesach
- **Nieprawidłowa składnia**: Źle sformułowany DSL Nextflow
- **Literówki w słowach kluczowych**: Błędnie napisane dyrektywy procesów
- **Niedopasowania kanałów**: Niezgodności typów

Serwer języka Nextflow podświetla te problemy w panelu problemów. Możesz sprawdzić je wcześnie, aby uniknąć błędów składniowych podczas uruchamiania pipeline'u.

### Podsumowanie

Możesz używać wykrywania błędów VS Code i panelu problemów do wychwytywania błędów składniowych i problemów przed uruchomieniem workflow'a, oszczędzając czas i zapobiegając frustracji.

### Co dalej?

Dowiedz się, jak efektywnie nawigować między procesami, modułami i definicjami w złożonych workflow'ach.

---

## 4. Nawigacja po kodzie i zarządzanie symbolami

Efektywna nawigacja jest kluczowa podczas pracy ze złożonymi workflow'ami rozciągającymi się na wiele plików. Aby to zrozumieć, zastąp definicję procesu w `basic_workflow.nf` importem dla modułu, który Ci dostarczyliśmy:

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

Jeśli najedziesz myszką na nazwę procesu, taką jak `FASTQC`, zobaczysz popup z interfejsem modułu (wejścia i wyjścia):

![Przejdź do definicji](img/syntax.png)

Ta funkcja jest szczególnie wartościowa podczas tworzenia workflow'ów, ponieważ pozwala zrozumieć interfejs modułu bez bezpośredniego otwierania pliku modułu.

Możesz szybko nawigować do dowolnej definicji procesu, modułu lub zmiennej używając **Ctrl/Cmd-klik**. Najedź myszką na link do pliku modułu na górze skryptu i podążaj za linkiem zgodnie z sugestią:

![Podążaj za linkiem](img/follow_link.png)

To samo działa dla nazw procesów. Wróć do `basic_workflow.nf` i spróbuj tego na nazwie procesu `FASTQC` w bloku workflow. To linkuje Cię bezpośrednio do nazwy procesu (która w tym przykładzie jest taka sama jak plik modułu, ale może być w połowie znacznie większego pliku).

Aby wrócić do miejsca, w którym byłeś, użyj **Alt+←** (lub **Ctrl+-** na Macu). To potężny sposób na eksplorację kodu bez gubienia swojego miejsca.

Teraz zbadajmy nawigację w bardziej złożonym workflow'ie używając `complex_workflow.nf` (pliku wyłącznie ilustracyjnego wspomnianego wcześniej). Ten workflow zawiera wiele procesów zdefiniowanych w oddzielnych plikach modułów, a także kilka wbudowanych. Chociaż złożone struktury wieloplikowe mogą być trudne do ręcznej nawigacji, możliwość przeskakiwania do definicji sprawia, że eksploracja staje się znacznie łatwiejsza.

1. Otwórz `complex_workflow.nf`
2. Nawiguj do definicji modułów
3. Użyj **Alt+←** (lub **Ctrl+-**), aby nawigować z powrotem
4. Nawiguj do nazwy procesu `FASTQC` w bloku workflow. To linkuje Cię bezpośrednio do nazwy procesu (która w tym przykładzie jest taka sama jak plik modułu, ale może być w połowie znacznie większego pliku).
5. Nawiguj z powrotem ponownie
6. Nawiguj do procesu `TRIM_GALORE` w bloku workflow. Jest on zdefiniowany wewnętrznie, więc nie przeniesie Cię do oddzielnego pliku, ale nadal pokaże Ci definicję procesu i nadal możesz nawigować z powrotem do miejsca, w którym byłeś.

### 4.2. Nawigacja po symbolach

Z nadal otwartym `complex_workflow.nf` możesz uzyskać przegląd wszystkich symboli w pliku, wpisując `@` w pasku wyszukiwania na górze VSCode (skrót klawiszowy to `Ctrl/Cmd+Shift+O`, ale może nie działać w Codespaces). To otwiera panel nawigacji po symbolach, który wyświetla wszystkie symbole w bieżącym pliku:

![Nawigacja po symbolach](img/symbols.png)

To pokazuje:

- Wszystkie definicje procesów
- Definicje workflow'ów (w tym pliku są zdefiniowane dwa workflow'y)
- Definicje funkcji

Zacznij pisać, aby filtrować wyniki.

### 4.3. Znajdź wszystkie odniesienia

Zrozumienie, gdzie proces lub zmienna jest używana w całej bazie kodu, może być bardzo pomocne. Na przykład, jeśli chcesz znaleźć wszystkie odniesienia do procesu `FASTQC`, zacznij od nawigacji do jego definicji. Możesz to zrobić, otwierając bezpośrednio `modules/fastqc.nf` lub używając funkcji szybkiej nawigacji VS Code z **Ctrl/Cmd-klik**, jak robiliśmy powyżej. Gdy jesteś przy definicji procesu, kliknij prawym przyciskiem myszy na nazwie procesu `FASTQC` i wybierz "Find All References" z menu kontekstowego, aby zobaczyć wszystkie wystąpienia, w których jest używany.

![Znajdź odniesienia](img/references.png)

Ta funkcja wyświetla wszystkie wystąpienia, w których `FASTQC` jest przywoływany w Twojej przestrzeni roboczej, w tym jego użycie w dwóch odrębnych workflow'ach. Ten wgląd jest kluczowy dla oceny potencjalnego wpływu modyfikacji procesu `FASTQC`.

### 4.4. Panel konspektu

Panel konspektu, znajdujący się w pasku bocznym Eksploratora (kliknij ![Ikona Eksploratora](img/files_icon.png)), zapewnia wygodny przegląd wszystkich symboli w Twoim bieżącym pliku. Ta funkcja pozwala szybko nawigować i zarządzać strukturą kodu, wyświetlając funkcje, zmienne i inne kluczowe elementy w widoku hierarchicznym.

![Panel konspektu](img/outline.png)

Użyj panelu konspektu, aby szybko nawigować do różnych części kodu bez używania przeglądarki plików.

### 4.5. Wizualizacja DAG

Rozszerzenie Nextflow VS Code może wizualizować Twój workflow jako skierowany graf acykliczny (DAG). To pomaga zrozumieć przepływ danych i zależności między procesami. Otwórz `complex_workflow.nf` i kliknij przycisk "Preview DAG" nad `workflow {` (drugi blok `workflow` w tym pliku):

![Podgląd DAG](img/dag_preview.png)

To jest tylko workflow 'wejściowy', ale możesz również podejrzeć DAG dla wewnętrznych workflow'ów, klikając przycisk "Preview DAG" nad workflow `RNASEQ_PIPELINE {` wyżej:

![Podgląd DAG wewnętrznego workflow'a](img/dag_preview_inner.png)

Dla tego workflow'a możesz używać węzłów w DAG do nawigacji do odpowiednich definicji procesów w kodzie. Kliknij węzeł, a przeniesie Cię do odpowiedniej definicji procesu w edytorze. Szczególnie gdy workflow rośnie do dużego rozmiaru, może to naprawdę pomóc w nawigacji po kodzie i zrozumieniu, jak procesy są połączone.

### Podsumowanie

Możesz efektywnie nawigować po złożonych workflow'ach używając przejścia do definicji, wyszukiwania symboli, znajdowania odniesień i wizualizacji DAG, aby zrozumieć strukturę kodu i zależności.

### Co dalej?

Dowiedz się, jak efektywnie pracować z wieloma powiązanymi plikami w większych projektach Nextflow.

## 5. Praca z wieloma plikami

Rzeczywiste programowanie w Nextflow obejmuje pracę z wieloma powiązanymi plikami. Zbadajmy, jak VS Code pomaga efektywnie zarządzać złożonymi projektami.

### 5.1. Szybka nawigacja po plikach

Z otwartym `complex_workflow.nf` zauważysz, że importuje on kilka modułów. Poćwiczmy szybką nawigację między nimi.

Naciśnij **Ctrl+P** (lub **Cmd+P**) i zacznij wpisywać "fast":

VS Code pokaże Ci pasujące pliki. Wybierz `modules/fastqc.nf`, aby tam natychmiast przeskoczyć. To jest znacznie szybsze niż klikanie przez eksplorator plików, gdy wiesz mniej więcej, jakiego pliku szukasz.

Spróbuj tego z innymi wzorcami:

- Wpisz "star", aby znaleźć plik modułu dopasowania STAR (`star.nf`)
- Wpisz "utils", aby znaleźć plik funkcji użytkowych (`utils.nf`)
- Wpisz "config", aby przeskoczyć do plików konfiguracyjnych (`nextflow.config`)

### 5.2. Podzielony edytor dla programowania wieloplikowego

Podczas pracy z modułami często musisz widzieć jednocześnie zarówno główny workflow, jak i definicje modułów. Skonfigurujmy to:

1. Otwórz `complex_workflow.nf`
2. Otwórz `modules/fastqc.nf` w nowej zakładce
3. Kliknij prawym przyciskiem myszy na zakładce `modules/fastqc.nf` i wybierz "Split Right"
4. Teraz możesz widzieć oba pliki obok siebie

![Podzielony edytor](img/split_editor.png)

To jest nieocenione, gdy:

- Sprawdzasz interfejsy modułów podczas pisania wywołań workflow, a podgląd nie wystarcza
- Porównujesz podobne procesy w różnych modułach
- Debugujesz przepływ danych między workflow'em a modułami

### 5.3. Wyszukiwanie w całym projekcie

Czasami musisz znaleźć, gdzie określone wzorce są używane w całym projekcie. Naciśnij `Ctrl/Cmd+Shift+F`, aby otworzyć panel wyszukiwania.

Spróbuj wyszukać `publishDir` w całej przestrzeni roboczej:

![Wyszukiwanie w projekcie](img/project_search.png)

To pokazuje Ci każdy plik, który używa katalogów publikacji, pomagając Ci:

- Zrozumieć wzorce organizacji wyjścia
- Znaleźć przykłady konkretnych dyrektyw
- Zapewnić spójność w modułach

### Podsumowanie

Możesz zarządzać złożonymi projektami wieloplikowymi używając szybkiej nawigacji po plikach, podzielonych edytorów i wyszukiwania w całym projekcie, aby efektywnie pracować z workflow'ami i modułami.

### Co dalej?

Dowiedz się, jak funkcje formatowania kodu i konserwacji utrzymują Twoje workflow'y zorganizowane i czytelne.

---

## 6. Formatowanie kodu i konserwacja

Odpowiednie formatowanie kodu jest istotne nie tylko dla estetyki, ale także dla zwiększenia czytelności, zrozumienia i łatwości aktualizacji złożonych workflow'ów.

### 6.1. Automatyczne formatowanie w akcji

Otwórz `basic_workflow.nf` i celowo zepsuj formatowanie:

- Usuń trochę wcięć: Zaznacz cały dokument i naciśnij `shift+tab` wiele razy, aby usunąć jak najwięcej wcięć.
- Dodaj dodatkowe spacje w losowych miejscach: w instrukcji `channel.fromPath` dodaj 30 spacji po `(`.
- Złam niektóre linie niezręcznie: Dodaj nową linię między operatorem `.view {` a napisem `Processing sample:`, ale nie dodawaj odpowiedniej nowej linii przed zamykającym nawiasem `}`.

Teraz naciśnij `Shift+Alt+F` (lub `Shift+Option+F` na MacOS), aby automatycznie sformatować:

VS Code natychmiast:

- Naprawia wcięcia, aby wyraźnie pokazać strukturę procesu
- Wyrównuje podobne elementy konsekwentnie
- Usuwa niepotrzebne białe znaki
- Utrzymuje czytelne łamanie linii

Zauważ, że automatyczne formatowanie może nie rozwiązać każdego problemu ze stylem kodu. Serwer języka Nextflow stara się utrzymać Twój kod w porządku, ale także szanuje Twoje osobiste preferencje w niektórych obszarach. Na przykład, jeśli usuniesz wcięcia wewnątrz bloku `script` procesu, formater pozostawi je takimi, jakie są, ponieważ możesz celowo preferować ten styl.

Obecnie nie ma ścisłego wymuszania stylu dla Nextflow, więc serwer języka oferuje pewną elastyczność. Jednak będzie konsekwentnie stosować reguły formatowania wokół definicji metod i funkcji, aby utrzymać przejrzystość.

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

To jest idealne do:

- Tymczasowego wyłączania części workflow'ów podczas programowania
- Dodawania wyjaśniających komentarzy do złożonych operacji na kanałach
- Dokumentowania sekcji workflow'a

Użyj **Ctrl+/** (lub **Cmd+/**) ponownie, aby odkomentować kod.

#### Zwijanie kodu dla przeglądu

W `complex_workflow.nf` zauważ małe strzałki obok definicji procesów. Kliknij je, aby zwinąć (zkolapsować) procesy:

![Zwijanie kodu](img/code_folding.png)

To daje Ci przegląd wysokiego poziomu struktury Twojego workflow'a bez gubienia się w szczegółach implementacji.

#### Dopasowywanie nawiasów

Umieść kursor obok dowolnego nawiasu `{` lub `}`, a VS Code podświetli pasujący nawias. Użyj **Ctrl+Shift+\\** (lub **Cmd+Shift+\\**), aby przeskakiwać między pasującymi nawiasami.

To jest kluczowe dla:

- Zrozumienia granic procesów
- Znajdowania brakujących lub dodatkowych nawiasów
- Nawigacji po zagnieżdżonych strukturach workflow'a

#### Wieloliniowe zaznaczanie i edycja

Do edycji wielu linii jednocześnie VS Code oferuje potężne możliwości wielokursorowe:

- **Wieloliniowe zaznaczanie**: Przytrzymaj **Ctrl+Alt** (lub **Cmd+Option** dla MacOS) i użyj klawiszy strzałek, aby zaznaczyć wiele linii
- **Wieloliniowe wcięcia**: Zaznacz wiele linii i użyj **Tab**, aby wciąć lub **Shift+Tab**, aby wyciąć całe bloki

To jest szczególnie przydatne do:

- Konsekwentnego wcięcia całych bloków procesów
- Dodawania komentarzy do wielu linii naraz
- Edycji podobnych definicji parametrów w wielu procesach

### Podsumowanie

Możesz utrzymywać czysty, czytelny kod używając automatycznego formatowania, funkcji komentowania, zwijania kodu, dopasowywania nawiasów i wieloliniowej edycji, aby efektywnie organizować złożone workflow'y.

### Co dalej?

Dowiedz się, jak VS Code integruje się z Twoim szerszym workflow'em programistycznym poza samą edycją kodu.

---

## 7. Integracja workflow'a programistycznego

VS Code dobrze integruje się z Twoim workflow'em programistycznym poza samą edycją kodu.

### 7.1. Integracja kontroli wersji

!!! note "Codespaces i integracja Git"

    Jeśli pracujesz w **GitHub Codespaces**, niektóre funkcje integracji Git mogą nie działać zgodnie z oczekiwaniami, szczególnie skróty klawiszowe dla kontroli źródła. Mogłeś również odmówić otwarcia katalogu jako repozytorium Git podczas początkowej konfiguracji, co jest w porządku dla celów szkoleniowych.

Jeśli Twój projekt jest repozytorium git (jak to jest), VS Code pokazuje:

- Zmodyfikowane pliki z kolorowymi wskaźnikami
- Status Git na pasku stanu
- Widoki różnic wbudowane
- Możliwości commitowania i pushowania

Otwórz panel kontroli źródła używając przycisku kontroli źródła (![Ikona kontroli źródła](img/source_control_icon.png)) (`Ctrl+Shift+G` lub `Cmd+Shift+G`, jeśli pracujesz z VSCode lokalnie), aby zobaczyć zmiany git i zatwierdzać commity bezpośrednio w edytorze.

![Panel kontroli źródła](img/source_control.png)

### 7.2. Uruchamianie i inspekcja workflow'ów

Uruchommy workflow i następnie zbadajmy wyniki. W zintegrowanym terminalu (`Ctrl+Shift+` backtick zarówno w Windows, jak i MacOS) uruchom podstawowy workflow:

```bash title="Run the basic workflow"
nextflow run basic_workflow.nf --input data/sample_data.csv --output_dir results
```

Podczas działania workflow'a zobaczysz wyjście w czasie rzeczywistym w terminalu. Po zakończeniu możesz użyć VS Code do inspekcji wyników bez opuszczania edytora:

1. **Nawiguj do katalogów roboczych**: Użyj eksploratora plików lub terminala, aby przeglądać `.nextflow/work`
2. **Otwórz pliki logów**: Kliknij ścieżki plików logów w wyjściu terminala, aby otworzyć je bezpośrednio w VS Code
3. **Zbadaj wyjścia**: Przeglądaj katalogi opublikowanych wyników w eksploratorze plików
4. **Wyświetl raporty wykonania**: Otwórz raporty HTML bezpośrednio w VS Code lub przeglądarce

To utrzymuje wszystko w jednym miejscu, zamiast przełączać się między wieloma aplikacjami.

### Podsumowanie

Możesz zintegrować VS Code z kontrolą wersji i wykonywaniem workflow'ów, aby zarządzać całym procesem programistycznym z jednego interfejsu.

### Co dalej?

Zobacz, jak wszystkie te funkcje IDE współpracują w Twoim codziennym workflow'ie programistycznym.

---

## 8. Podsumowanie i szybkie notatki

Oto kilka szybkich notatek na temat każdej z funkcji IDE omówionych powyżej:

### 8.1. Rozpoczynanie nowej funkcjonalności

1. **Szybkie otwieranie plików** (`Ctrl+P` lub `Cmd+P`), aby znaleźć odpowiednie istniejące moduły
2. **Podzielony edytor**, aby wyświetlić podobne procesy obok siebie
3. **Nawigacja po symbolach** (`Ctrl+Shift+O` lub `Cmd+Shift+O`), aby zrozumieć strukturę pliku
4. **Autouzupełnianie**, aby szybko pisać nowy kod

### 8.2. Debugowanie problemów

1. **Panel problemów** (`Ctrl+Shift+M` lub `Cmd+Shift+M`), aby zobaczyć wszystkie błędy naraz
2. **Przejdź do definicji** (`Ctrl-klik` lub `Cmd-klik`), aby zrozumieć interfejsy procesów
3. **Znajdź wszystkie odniesienia**, aby zobaczyć, jak procesy są używane
4. **Wyszukiwanie w całym projekcie**, aby znaleźć podobne wzorce lub problemy

### 8.3. Refaktoryzacja i ulepszanie

1. **Wyszukiwanie w całym projekcie** (`Ctrl+Shift+F` lub `Cmd+Shift+F`), aby znaleźć wzorce
2. **Automatyczne formatowanie** (`Shift+Alt+F` lub `Shift+Option+F`), aby utrzymać spójność
3. **Zwijanie kodu**, aby skupić się na strukturze
4. **Integracja Git**, aby śledzić zmiany

---

## Podsumowanie

Właśnie odbyłeś błyskawiczną wycieczkę po funkcjach IDE VS Code dla programowania w Nextflow. Te narzędzia uczynią Cię znacznie bardziej produktywnym poprzez:

- **Redukcję błędów** dzięki sprawdzaniu składni w czasie rzeczywistym
- **Przyspieszenie programowania** dzięki inteligentnemu autouzupełnianiu
- **Poprawę nawigacji** w złożonych wieloplikowych workflow'ach
- **Utrzymanie jakości** dzięki konsekwentnemu formatowaniu
- **Zwiększenie zrozumienia** dzięki zaawansowanemu podświetlaniu i wizualizacji struktury

Nie oczekujemy, że zapamiętasz wszystko, ale teraz wiesz, że te funkcje istnieją i będziesz mógł je znaleźć, gdy będziesz ich potrzebować. W miarę jak będziesz kontynuować programowanie workflow'ów Nextflow, te funkcje IDE staną się drugą naturą, pozwalając Ci skupić się na pisaniu wysokiej jakości kodu, zamiast zmagać się ze składnią i strukturą.

### Co dalej?

Zastosuj te umiejętności IDE podczas pracy z innymi modułami szkoleniowymi, na przykład:

- **[nf-test](nf-test.md)**: Twórz kompleksowe zestawy testów dla swoich workflow'ów
- **[Hello nf-core](../../hello_nf-core/)**: Buduj pipeline'y o jakości produkcyjnej ze standardami społeczności

Prawdziwa moc tych funkcji IDE ujawnia się, gdy pracujesz nad większymi, bardziej złożonymi projektami. Zacznij stopniowo włączać je do swojego workflow'a - w ciągu kilku sesji staną się drugą naturą i zmienią sposób, w jaki podchodzisz do programowania w Nextflow.

Od wychwytywania błędów, zanim Cię spowolnią, po nawigację po złożonych bazach kodu z łatwością - te narzędzia uczynią Cię bardziej pewnym siebie i efektywnym programistą.

Miłego kodowania!
