# Część 1: Uruchamianie podstawowych operacji

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Tłumaczenie wspomagane przez AI - [dowiedz się więcej i zasugeruj ulepszenia](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

W tej pierwszej części kursu szkoleniowego Nextflow dla Bioimagingu użyjemy bardzo prostego, niezależnego od dziedziny przykładu Hello World, aby zademonstrować podstawowe operacje i wskazać odpowiednie komponenty kodu Nextflow.

## 1. Uruchomienie workflow'u

Udostępniamy skrypt workflow'u o nazwie `hello-world.nf`, który przyjmuje dane wejściowe za pomocą argumentu wiersza poleceń o nazwie `--greeting` i tworzy plik tekstowy zawierający to powitanie.
Na razie nie będziemy przeglądać kodu; najpierw zobaczmy, jak wygląda jego uruchomienie.

### 1.1. Uruchomienie workflow'u i monitorowanie wykonania

W terminalu uruchom następujące polecenie:

```bash
nextflow run hello-world.nf --greeting 'Hello World!'
```

Wyjście w konsoli powinno wyglądać mniej więcej tak:

```console title="Wyjście" linenums="1"
 N E X T F L O W   ~  version 25.04.3

Launching `hello-world.nf` [goofy_torvalds] DSL2 - revision: c33d41f479

executor >  local (1)
[a3/7be2fa] sayHello | 1 of 1 ✔
```

Gratulacje, właśnie uruchomiłeś Swój pierwszy workflow Nextflow!

Najważniejszym wyjściem jest tutaj ostatnia linia (linia 6):

```console title="Wyjście" linenums="6"
[a3/7be2fa] sayHello | 1 of 1 ✔
```

To mówi nam, że proces `sayHello` został pomyślnie wykonany raz (`1 of 1 ✔`).

To świetnie, ale możesz się zastanawiać: gdzie jest wyjście?

### 1.2. Znalezienie pliku wyjściowego w katalogu `results`

Ten workflow jest skonfigurowany do publikowania Swojego wyjścia do katalogu o nazwie `results`.
Jeśli spojrzysz na Swój bieżący katalog, zobaczysz, że kiedy uruchomiłeś workflow'a, Nextflow utworzył nowy katalog o nazwie `results`, który zawiera plik o nazwie `output.txt`.

```console title="results/" linenums="1"
results
└── output.txt
```

Otwórz plik; jego zawartość powinna odpowiadać powitaniu, które podałeś w wierszu poleceń.

<details>
  <summary>Zawartość pliku</summary>

```console title="results/output.txt" linenums="1"
Hello World!
```

</details>

Świetnie, nasz workflow zrobił to, co powinien!

Należy jednak pamiętać, że 'opublikowany' wynik jest kopią (lub w niektórych przypadkach dowiązaniem symbolicznym) rzeczywistego wyjścia wygenerowanego przez Nextflow'a podczas wykonywania workflow'u.

Teraz zajrzymy pod maskę, aby zobaczyć, gdzie Nextflow faktycznie wykonał pracę.

!!! warning "Ostrzeżenie"

    Nie wszystkie workflow'y będą skonfigurowane do publikowania wyjść do katalogu results; nazwa katalogu może być też inna.
    Nieco dalej w tej sekcji pokażemy, jak dowiedzieć się, gdzie to zachowanie jest określone.

### 1.3. Znalezienie oryginalnego wyjścia i logów w katalogu `work/`

Kiedy uruchamiasz workflow, Nextflow tworzy odrębny 'katalog zadania' dla każdego pojedynczego wywołania każdego procesu w workflow'ie (czyli każdego kroku w pipeline'ie).
Dla każdego z nich przygotuje niezbędne wejścia, wykona odpowiednie instrukcje i zapisze wyjścia oraz pliki logów w tym samym katalogu, automatycznie nazwanym przy użyciu hasha w celu zapewnienia unikalności.

Wszystkie te katalogi zadań będą znajdować się w katalogu o nazwie `work` w Twoim bieżącym katalogu (gdzie uruchamiasz polecenie).

To może brzmieć zagmatwanie, więc zobaczmy, jak to wygląda w praktyce.

Wracając do wyjścia konsoli dla workflow'u, który uruchomiliśmy wcześniej, mieliśmy tę linię:

```console title="Fragment wyjścia polecenia" linenums="6"
[a3/7be2fa] sayHello | 1 of 1 ✔
```

Widzisz, jak linia zaczyna się od `[a3/7be2fa]`?
To skrócona forma ścieżki katalogu zadania dla tego jednego wywołania procesu i mówi Ci, gdzie znaleźć wyjście wywołania procesu `sayHello` w ścieżce katalogu `work/`.

Możesz znaleźć pełną ścieżkę, wpisując następujące polecenie (zastępując `a3/7be2fa` tym, co widzisz w Swoim własnym terminalu) i naciskając klawisz tab, aby automatycznie uzupełnić ścieżkę, lub dodając gwiazdkę:

```bash
tree work/a3/7be2fa*
```

Powinno to zwrócić pełną ścieżkę katalogu: `work/a3/7be2fa7be2fad5e71e5f49998f795677fd68`

Zobaczmy, co tam jest.

!!! Tip "Wskazówka"

    Jeśli przeglądasz zawartość podkatalogu zadania w eksploratorze plików VSCode, zobaczysz wszystkie pliki od razu.
    Jednak pliki logów są ustawione jako niewidoczne w terminalu, więc jeśli chcesz użyć `ls` lub `tree` do ich przeglądania, będziesz musiał ustawić odpowiednią opcję wyświetlania niewidocznych plików.

    ```bash
    tree -a work
    ```

Dokładne nazwy podkatalogów będą różne w Twoim systemie.

<details>
  <summary>Zawartość katalogu</summary>

```console title="work/"
work
└── a3
    └── 7be2fad5e71e5f49998f795677fd68
        ├── .command.begin
        ├── .command.err
        ├── .command.log
        ├── .command.out
        ├── .command.run
        ├── .command.sh
        ├── .exitcode
        └── output.txt
```

</details>

Powinieneś od razu rozpoznać plik `output.txt`, który jest w rzeczywistości oryginalnym wyjściem procesu `sayHello`, opublikowanym w katalogu `results`.
Jeśli go otworzysz, znajdziesz ponownie powitanie `Hello World!`.

<details>
  <summary>Zawartość pliku output.txt</summary>

```console title="work/a3/7be2fa7be2fad5e71e5f49998f795677fd68/output.txt" linenums="1"
Hello World!
```

</details>

Co więc z tymi wszystkimi innymi plikami?

To są pliki pomocnicze i logi, które Nextflow zapisał jako część wykonania zadania:

- **`.command.begin`**: Plik wartowniczy utworzony natychmiast po uruchomieniu zadania.
- **`.command.err`**: Komunikaty o błędach (`stderr`) wyemitowane przez wywołanie procesu
- **`.command.log`**: Kompletny log wyjścia wyemitowany przez wywołanie procesu
- **`.command.out`**: Regularne wyjście (`stdout`) wywołania procesu
- **`.command.run`**: Pełny skrypt uruchomiony przez Nextflow'a w celu wykonania wywołania procesu
- **`.command.sh`**: Polecenie, które faktycznie zostało uruchomione przez wywołanie procesu
- **`.exitcode`**: Kod wyjścia wynikający z polecenia

Plik `.command.sh` jest szczególnie przydatny, ponieważ pokazuje główne polecenie wykonane przez Nextflow'a, nie włączając całej księgowości i konfiguracji zadania/środowiska.

<details>
  <summary>Zawartość pliku</summary>

```console title="work/a3/7be2fa7be2fad5e71e5f49998f795677fd68/.command.sh" linenums="1"
#!/bin/bash -ue
echo 'Hello World!' > output.txt

```

</details>

!!! Tip "Wskazówka"

    Gdy coś pójdzie nie tak i musisz rozwiązać problem, może być przydatne sprawdzenie skryptu `command.sh`, aby zweryfikować dokładnie, jakie polecenie Nextflow skomponował na podstawie instrukcji workflow'u, interpolacji zmiennych i tak dalej.

### 1.4. Ćwiczenie opcjonalne: ponowne uruchomienie z różnymi powitaniami

Spróbuj ponownie uruchomić workflow kilka razy z różnymi wartościami dla argumentu `--greeting`, a następnie sprawdź zarówno zawartość katalogu `results/`, jak i katalogów zadań.

Zauważ, jak wyjścia i logi izolowanych katalogów zadań są zachowywane, podczas gdy zawartość katalogu `results` jest nadpisywana wyjściem kolejnych wykonań.

### Podsumowanie

Wiesz, jak uruchomić prosty skrypt Nextflow, monitorować jego wykonanie i znaleźć jego wyjścia.

### Co dalej?

Dowiedz się, jak czytać podstawowy skrypt Nextflow i zidentyfikować, jak jego komponenty odnoszą się do jego funkcjonalności.

---

## 2. Przegląd początkowego skryptu workflow'u Hello World

To, co zrobiliśmy, to w zasadzie potraktowanie skryptu workflow'u jak czarnej skrzynki.
Teraz, gdy zobaczyliśmy, co robi, otwórzmy pudełko i zajrzyjmy do środka.

_Celem tutaj nie jest zapamiętanie składni kodu Nextflow, ale zbudowanie podstawowej intuicji dotyczącej głównych komponentów i sposobu ich organizacji._

### 2.1. Przegląd ogólnej struktury kodu

Otwórzmy skrypt `hello-world.nf` w panelu edytora.

<details>
  <summary>Kod</summary>

```groovy title="hello-world.nf" linenums="1"
#!/usr/bin/env nextflow

/*
 * Użyj echo do wypisania pozdrowienia do pliku
 */
process sayHello {

    publishDir 'results', mode: 'copy'

    input:
    val greeting

    output:
    path 'output.txt'

    script:
    """
    echo '$greeting' > output.txt
    """
}

workflow {

    // wyemituj pozdrowienie
    sayHello(params.greeting)
}
```

</details>

Skrypt Nextflow obejmuje dwa główne typy podstawowych komponentów: jeden lub więcej **procesów** oraz sam **workflow**.
Każdy **proces** opisuje, jakie operacje powinien wykonać odpowiedni krok w pipeline'ie, natomiast **workflow** opisuje logikę przepływu danych łączącą poszczególne kroki.

Przyjrzyjmy się najpierw bliżej blokowi **process**, a następnie przyjrzymy się blokowi **workflow**.

### 2.2. Definicja `process`

Pierwszy blok kodu opisuje **proces**.
Definicja procesu zaczyna się od słowa kluczowego `process`, po którym następuje nazwa procesu, a na końcu treść procesu oddzielona klamrami.
Treść procesu musi zawierać blok skryptu, który określa polecenie do uruchomienia; może to być wszystko, co można uruchomić w terminalu wiersza poleceń.

Tutaj mamy **proces** o nazwie `sayHello`, który przyjmuje zmienną **wejściową** o nazwie `greeting` i zapisuje Swoje **wyjście** do pliku o nazwie `output.txt`.

<details>
  <summary>Kod</summary>

```groovy title="hello-world.nf" linenums="3"
/*
 * Użyj echo do wypisania pozdrowienia do pliku
 */
process sayHello {

    publishDir 'results', mode: 'copy'

    input:
    val greeting

    output:
    path 'output.txt'

    script:
    """
    echo '$greeting' > output.txt
    """
}
```

</details>

To bardzo minimalna definicja procesu, która zawiera tylko definicję `input`, definicję `output` oraz `script` do wykonania.

Definicja `input` zawiera kwalifikator `val`, który mówi Nextflow'owi, aby spodziewał się wartości jakiegoś rodzaju (może to być ciąg znaków, liczba, cokolwiek).

Definicja `output` zawiera kwalifikator `path`, który mówi Nextflow'owi, że powinno to być traktowane jako ścieżka (obejmuje zarówno ścieżki katalogów, jak i pliki).

!!! Tip "Wskazówka"

    Definicja wyjścia nie _określa_, jakie wyjście zostanie utworzone.
    Po prostu _deklaruje_, gdzie znaleźć oczekiwane pliki wyjściowe, aby Nextflow mógł je znaleźć po zakończeniu wykonania.

    Jest to konieczne do weryfikacji, że polecenie zostało wykonane pomyślnie oraz do przekazania wyjścia do procesów następnych, jeśli jest to potrzebne.
    Wyjście wytworzone, które nie pasuje do tego, co jest zadeklarowane w bloku wyjścia, nie zostanie przekazane do procesów następnych.

W rzeczywistym pipeline'ie proces zazwyczaj zawiera dodatkowe informacje, takie jak dyrektywy procesu, które przedstawimy za chwilę.

### 2.3. Definicja `workflow`

Drugi blok kodu opisuje sam **workflow**.
Definicja workflow'u zaczyna się od słowa kluczowego `workflow`, po którym następuje opcjonalna nazwa, a następnie treść workflow'u oddzielona klamrami.

Tutaj mamy **workflow**, który składa się z jednego wywołania procesu `sayHello`; przyjmuje on wejście `params.greeting`, które przechowuje wartość podaną przez nas do parametru `--greeting`.

```groovy title="hello-world.nf" linenums="22"
workflow {

    // wyemituj pozdrowienie
    sayHello(params.greeting)
}
```

To bardzo minimalna definicja **workflow'u**.
W rzeczywistym pipeline'ie workflow zazwyczaj zawiera wiele wywołań **procesów** połączonych **kanałami**, i mogą być ustawione domyślne wartości dla zmiennych wejściowych.

Zobaczymy to w akcji, gdy uruchomimy nf-core/molkart w Części 2 kursu.

### 2.4. System `params` parametrów wiersza poleceń

`params.greeting`, który przekazujemy do wywołania procesu `sayHello()`, to elegancki fragment kodu Nextflow i warto poświęcić na niego dodatkową minutę.

Jak wspomniano powyżej, tak przekazujemy wartość parametru wiersza poleceń `--greeting` do wywołania procesu `sayHello()`.
W rzeczywistości samo zadeklarowanie `params.someParameterName` umożliwi nam podanie workflow'owi parametru o nazwie `--someParameterName` z wiersza poleceń.

!!! Tip "Wskazówka"

    Te parametry workflow'u zadeklarowane przy użyciu systemu `params` zawsze przyjmują dwa myślniki (`--`).
    To odróżnia je od parametrów poziomu Nextflow'a, które przyjmują tylko jeden myślnik (`-`).

### Podsumowanie

Wiesz teraz, jak jest zbudowany prosty workflow Nextflow i jak podstawowe komponenty odnoszą się do jego funkcjonalności.

### Co dalej?

Naucz się wygodnie zarządzać wykonaniami workflow'u.

---

## 3. Zarządzanie wykonaniami workflow'u

Wiedza, jak uruchamiać workflow'e i pobierać wyjścia, jest świetna, ale szybko odkryjesz, że jest kilka innych aspektów zarządzania workflow'em, które ułatwią Ci życie.

Tutaj pokażemy Ci, jak wykorzystać funkcję `resume` do ponownego uruchomienia tego samego workflow'u, jak sprawdzić logi wykonania za pomocą `nextflow log` oraz jak usunąć starsze katalogi robocze za pomocą `nextflow clean`.

### 3.1. Ponowne uruchomienie workflow'u z `-resume`

Czasami będziesz chciał ponownie uruchomić pipeline, który już wcześniej uruchomiłeś, bez powtarzania jakiejkolwiek pracy, która została już pomyślnie zakończona.

Nextflow ma opcję o nazwie `-resume`, która pozwala Ci to zrobić.
W szczególności, w tym trybie wszystkie procesy, które zostały już uruchomione z dokładnie tym samym kodem, ustawieniami i wejściami, zostaną pominięte.
Oznacza to, że Nextflow uruchomi tylko procesy, które dodałeś lub zmodyfikowałeś od ostatniego uruchomienia, lub którym przekazujesz nowe ustawienia lub wejścia.

Są dwie kluczowe zalety tego podejścia:

- Jeśli jesteś w trakcie opracowywania pipeline'u, możesz iterować szybciej, ponieważ musisz uruchomić tylko proces(y), nad którymi aktywnie pracujesz, aby przetestować Swoje zmiany.
- Jeśli uruchamiasz pipeline w środowisku produkcyjnym i coś pójdzie nie tak, w wielu przypadkach możesz naprawić problem i ponownie uruchomić pipeline, a on wznowi działanie od punktu awarii, co może zaoszczędzić Ci dużo czasu i mocy obliczeniowej.

Aby z niego skorzystać, po prostu dodaj `-resume` do Swojego polecenia i uruchom je:

```bash
nextflow run hello-world.nf --greeting 'Hello World!' -resume
```

??? success "Wyjście polecenia"

    ```console
    N E X T F L O W   ~  version 25.04.3

    Launching `hello-world.nf` [tiny_noyce] DSL2 - revision: c33d41f479

    [a3/7be2fa] process > sayHello [100%] 1 of 1, cached: 1 ✔
    ```

Zwróć uwagę na fragment `cached:`, który został dodany w linii statusu procesu (linia 5), co oznacza, że Nextflow rozpoznał, iż już wykonał tę pracę i po prostu ponownie wykorzystał wynik z poprzedniego pomyślnego uruchomienia.

Możesz również zobaczyć, że hash podkatalogu roboczego jest taki sam jak w poprzednim uruchomieniu.
Nextflow dosłownie wskazuje Ci na poprzednie wykonanie, mówiąc "Już to zrobiłem tam."

!!! Tip "Wskazówka"

    Gdy ponownie uruchamiasz pipeline z `resume`, Nextflow nie nadpisuje żadnych plików zapisanych do katalogu `publishDir` przez żadne wywołanie procesu, które zostało wcześniej pomyślnie uruchomione.

### 3.2. Przegląd logu wcześniejszych wykonań

Za każdym razem, gdy uruchamiasz workflow Nextflow, linia jest zapisywana do pliku logu o nazwie `history` w ukrytym katalogu o nazwie `.nextflow` w bieżącym katalogu roboczym.

Bardziej wygodnym sposobem dostępu do tych informacji jest użycie polecenia `nextflow log`.

```bash
nextflow log
```

To wyświetli zawartość pliku logu w terminalu, pokazując Ci znacznik czasu, nazwę uruchomienia, status i pełną linię poleceń dla każdego uruchomienia Nextflow, które zostało uruchomione z bieżącego katalogu roboczego.

### 3.3. Usuwanie starszych katalogów roboczych

Podczas procesu rozwoju zazwyczaj uruchamiasz Swoje robocze pipeline'y wiele razy, co może prowadzić do nagromadzenia bardzo wielu plików w wielu podkatalogach.
Ponieważ podkatalogi są nazwane losowo, trudno jest stwierdzić na podstawie ich nazw, które są starsze, a które nowsze.

Nextflow zawiera wygodne podpolecenie `clean`, które może automatycznie usuwać podkatalogi robocze dla wcześniejszych uruchomień, o które już Ci nie zależy, z kilkoma [opcjami](https://www.nextflow.io/docs/latest/reference/cli.html#clean) do kontrolowania, co zostanie usunięte.

Możesz użyć logu Nextflow, aby wyszukać uruchomienie na podstawie jego znacznika czasu i/lub linii poleceń, a następnie użyć `nextflow clean -before <run_name> -f`, aby usunąć katalogi robocze z wcześniejszych uruchomień.

!!! Warning "Ostrzeżenie"

    Usuwanie podkatalogów roboczych z wcześniejszych uruchomień usuwa je z pamięci podręcznej Nextflow'a i usuwa wszystkie wyjścia, które były przechowywane w tych katalogach.
    Oznacza to, że przerywa zdolność Nextflow'a do wznowienia wykonania bez ponownego uruchamiania odpowiednich procesów.

    Jesteś odpowiedzialny za zapisanie wszelkich wyjść, na których Ci zależy lub na których planujesz polegać! Jeśli używasz dyrektywy `publishDir` w tym celu, upewnij się, że używasz trybu `copy`, a nie trybu `symlink`.

### Podsumowanie

Wiesz, jak ponownie uruchomić pipeline bez powtarzania kroków, które zostały już uruchomione w identyczny sposób, sprawdzić log wykonania i użyć polecenia `nextflow clean` do czyszczenia starych katalogów roboczych.

### Co dalej?

Teraz, gdy rozumiesz podstawowe operacje Nextflow, jesteś gotowy, aby uruchomić prawdziwy pipeline bioimagingu z nf-core/molkart.
