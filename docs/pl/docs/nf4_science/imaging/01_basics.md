# Część 1: Podstawowe operacje

W pierwszej części szkolenia Nextflow dla bioimagingu użyjemy bardzo prostego, niezależnego od dziedziny przykładu Hello World, aby zademonstrować podstawowe operacje i wskazać odpowiadające im komponenty kodu Nextflow.

## 1. Uruchom workflow'a

Przygotowaliśmy dla Ciebie skrypt workflow'a o nazwie `hello-world.nf`, który przyjmuje dane wejściowe przez argument wiersza poleceń o nazwie `--greeting` i tworzy plik tekstowy zawierający to powitanie.
Na razie nie będziemy przeglądać kodu; najpierw zobaczmy, jak wygląda jego uruchomienie.

### 1.1. Uruchom workflow'a i monitoruj wykonanie

W terminalu uruchom następujące polecenie:

```bash
nextflow run hello-world.nf --greeting 'Hello World!'
```

Wyjście w konsoli powinno wyglądać mniej więcej tak:

```console title="Output" linenums="1"
 N E X T F L O W   ~  version 25.04.3

Launching `hello-world.nf` [goofy_torvalds] DSL2 - revision: c33d41f479

executor >  local (1)
[a3/7be2fa] sayHello | 1 of 1 ✔
```

Gratulacje, właśnie uruchomiłeś swój pierwszy workflow Nextflow'a!

Najważniejszym wyjściem jest tutaj ostatnia linia (linia 6):

```console title="Output" linenums="6"
[a3/7be2fa] sayHello | 1 of 1 ✔
```

Informuje nas ona, że proces `sayHello` został pomyślnie wykonany raz (`1 of 1 ✔`).

To świetnie, ale możesz się zastanawiać: gdzie jest wyjście?

### 1.2. Znajdź plik wyjściowy w katalogu `results`

Ten workflow jest skonfigurowany tak, aby publikować swoje wyjście do katalogu o nazwie `results`.
Jeśli spojrzysz na swój bieżący katalog, zobaczysz, że gdy uruchomiłeś workflow'a, Nextflow utworzył nowy katalog o nazwie `results`, który zawiera plik o nazwie `output.txt`.

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

Pamiętaj jednak, że „opublikowany" wynik jest kopią (lub w niektórych przypadkach dowiązaniem symbolicznym) rzeczywistego wyjścia wygenerowanego przez Nextflow'a podczas wykonywania workflow'a.

Teraz zajrzymy pod maskę, aby zobaczyć, gdzie Nextflow faktycznie wykonał pracę.

!!! warning "Ostrzeżenie"

    Nie wszystkie workflow'y będą skonfigurowane tak, aby publikować wyjścia do katalogu results, i/lub nazwa katalogu może być inna.
    Nieco dalej w tej sekcji pokażemy Ci, jak dowiedzieć się, gdzie to zachowanie jest określone.

### 1.3. Znajdź oryginalne wyjście i logi w katalogu `work/`

Gdy uruchamiasz workflow'a, Nextflow tworzy odrębny „katalog zadania" dla każdego pojedynczego wywołania każdego procesu w workflow'ie (= każdego kroku w pipeline'ie).
Dla każdego z nich przygotuje niezbędne dane wejściowe, wykona odpowiednie instrukcje i zapisze wyjścia oraz pliki logów w tym jednym katalogu, który jest automatycznie nazywany przy użyciu skrótu, aby był unikalny.

Wszystkie te katalogi zadań będą znajdować się w katalogu o nazwie `work` w Twoim bieżącym katalogu (skąd uruchamiasz polecenie).

To może brzmieć zagmatwanie, więc zobaczmy, jak to wygląda w praktyce.

Wracając do wyjścia konsoli dla workflow'a, który uruchomiliśmy wcześniej, mieliśmy tę linię:

```console title="Excerpt of command output" linenums="6"
[a3/7be2fa] sayHello | 1 of 1 ✔
```

Widzisz, jak linia zaczyna się od `[a3/7be2fa]`?
To jest skrócona forma ścieżki katalogu zadania dla tego jednego wywołania procesu i mówi Ci, gdzie znaleźć wyjście wywołania procesu `sayHello` w ścieżce katalogu `work/`.

Możesz znaleźć pełną ścieżkę, wpisując następujące polecenie (zastępując `a3/7be2fa` tym, co widzisz w swoim własnym terminalu) i naciskając klawisz tab, aby automatycznie uzupełnić ścieżkę, lub dodając gwiazdkę:

```bash
tree work/a3/7be2fa*
```

Powinno to dać pełną ścieżkę katalogu: `work/a3/7be2fa7be2fad5e71e5f49998f795677fd68`

Zobaczmy, co tam jest.

!!! Tip "Wskazówka"

    Jeśli przeglądasz zawartość podkatalogu zadania w eksploratorze plików VSCode, zobaczysz wszystkie pliki od razu.
    Jednak pliki logów są ustawione jako niewidoczne w terminalu, więc jeśli chcesz użyć `ls` lub `tree` do ich wyświetlenia, musisz ustawić odpowiednią opcję wyświetlania niewidocznych plików.

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

Powinieneś od razu rozpoznać plik `output.txt`, który w rzeczywistości jest oryginalnym wyjściem procesu `sayHello`, które zostało opublikowane w katalogu `results`.
Jeśli go otworzysz, znajdziesz ponownie powitanie `Hello World!`.

<details>
  <summary>Zawartość pliku output.txt</summary>

```console title="work/a3/7be2fa7be2fad5e71e5f49998f795677fd68/output.txt" linenums="1"
Hello World!
```

</details>

A co z tymi wszystkimi innymi plikami?

To pliki pomocnicze i logi, które Nextflow zapisał w ramach wykonania zadania:

- **`.command.begin`**: Plik wartowniczy utworzony natychmiast po uruchomieniu zadania.
- **`.command.err`**: Komunikaty o błędach (`stderr`) wyemitowane przez wywołanie procesu
- **`.command.log`**: Kompletne wyjście logu wyemitowane przez wywołanie procesu
- **`.command.out`**: Regularne wyjście (`stdout`) wywołania procesu
- **`.command.run`**: Pełny skrypt uruchomiony przez Nextflow'a w celu wykonania wywołania procesu
- **`.command.sh`**: Polecenie, które faktycznie zostało uruchomione przez wywołanie procesu
- **`.exitcode`**: Kod wyjścia wynikający z polecenia

Plik `.command.sh` jest szczególnie przydatny, ponieważ pokazuje główne polecenie wykonane przez Nextflow'a, nie uwzględniając całej księgowości i konfiguracji zadania/środowiska.

<details>
  <summary>Zawartość pliku</summary>

```console title="work/a3/7be2fa7be2fad5e71e5f49998f795677fd68/.command.sh" linenums="1"
#!/bin/bash -ue
echo 'Hello World!' > output.txt

```

</details>

!!! Tip "Wskazówka"

    Gdy coś pójdzie nie tak i musisz rozwiązać problem, przydatne może być spojrzenie na skrypt `command.sh`, aby sprawdzić dokładnie, jakie polecenie Nextflow skomponował na podstawie instrukcji workflow'a, interpolacji zmiennych i tak dalej.

### 1.4. Ćwiczenie opcjonalne: uruchom ponownie z różnymi powitaniami

Spróbuj uruchomić workflow'a kilka razy z różnymi wartościami dla argumentu `--greeting`, a następnie przyjrzyj się zarówno zawartości katalogu `results/`, jak i katalogów zadań.

Zauważ, jak wyjścia i logi izolowanych katalogów zadań są zachowywane, podczas gdy zawartość katalogu `results` jest nadpisywana przez wyjście kolejnych wykonań.

### Podsumowanie

Wiesz, jak uruchomić prosty skrypt Nextflow, monitorować jego wykonanie i znaleźć jego wyjścia.

### Co dalej?

Naucz się czytać podstawowy skrypt Nextflow i identyfikować, jak jego komponenty odnoszą się do jego funkcjonalności.

---

## 2. Przeanalizuj początkowy skrypt workflow'a Hello World

To, co zrobiliśmy, to w zasadzie potraktowanie skryptu workflow'a jak czarnej skrzynki.
Teraz, gdy zobaczyliśmy, co robi, otwórzmy skrzynkę i zajrzyjmy do środka.

_Celem nie jest zapamiętanie składni kodu Nextflow, ale wyrobienie podstawowej intuicji dotyczącej głównych komponentów i sposobu ich organizacji._

### 2.1. Przeanalizuj ogólną strukturę kodu

Otwórzmy skrypt `hello-world.nf` w panelu edytora.

<details>
  <summary>Kod</summary>

```groovy title="hello-world.nf" linenums="1"
#!/usr/bin/env nextflow

/*
 * Use echo to print a greeting to a file
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

    // emit a greeting
    sayHello(params.greeting)
}
```

</details>

Skrypt Nextflow obejmuje dwa główne typy podstawowych komponentów: jeden lub więcej **procesów** oraz sam **workflow**.
Każdy **proces** opisuje, jakie operacje powinien wykonać odpowiedni krok w pipeline'ie, podczas gdy **workflow** opisuje logikę przepływu danych łączącą różne kroki.

Przyjrzyjmy się najpierw bliżej blokowi **process**, a następnie spojrzymy na blok **workflow**.

### 2.2. Definicja `process`

Pierwszy blok kodu opisuje **proces**.
Definicja procesu zaczyna się od słowa kluczowego `process`, po którym następuje nazwa procesu, a na końcu ciało procesu ograniczone nawiasami klamrowymi.
Ciało procesu musi zawierać blok skryptu, który określa polecenie do uruchomienia, które może być wszystkim, co mógłbyś uruchomić w terminalu wiersza poleceń.

Tutaj mamy **proces** o nazwie `sayHello`, który przyjmuje zmienną **wejściową** o nazwie `greeting` i zapisuje swoje **wyjście** do pliku o nazwie `output.txt`.

<details>
  <summary>Kod</summary>

```groovy title="hello-world.nf" linenums="3"
/*
 * Use echo to print a greeting to a file
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

Definicja `input` zawiera kwalifikator `val`, który mówi Nextflow'owi, aby oczekiwał wartości jakiegoś rodzaju (może to być ciąg znaków, liczba, cokolwiek).

Definicja `output` zawiera kwalifikator `path`, który mówi Nextflow'owi, że powinno to być traktowane jako ścieżka (obejmuje zarówno ścieżki katalogów, jak i pliki).

!!! Tip "Wskazówka"

    Definicja wyjścia nie _określa_, jakie wyjście zostanie utworzone.
    Po prostu _deklaruje_, gdzie znaleźć oczekiwane pliki wyjściowe, aby Nextflow mógł ich szukać po zakończeniu wykonania.

    Jest to niezbędne do weryfikacji, czy polecenie zostało wykonane pomyślnie, oraz do przekazania wyjścia do procesów podrzędnych, jeśli jest to potrzebne.
    Wyjście wygenerowane, które nie pasuje do tego, co jest zadeklarowane w bloku wyjścia, nie zostanie przekazane do procesów podrzędnych.

W rzeczywistym pipeline'ie proces zwykle zawiera dodatkowe informacje, takie jak dyrektywy procesu, które przedstawimy za chwilę.

### 2.3. Definicja `workflow`

Drugi blok kodu opisuje sam **workflow**.
Definicja workflow'a zaczyna się od słowa kluczowego `workflow`, po którym następuje opcjonalna nazwa, a następnie ciało workflow'a ograniczone nawiasami klamrowymi.

Tutaj mamy **workflow**, który składa się z jednego wywołania procesu `sayHello`, które przyjmuje wejście `params.greeting`, przechowujące wartość, którą podaliśmy parametrowi `--greeting`.

```groovy title="hello-world.nf" linenums="22"
workflow {

    // emit a greeting
    sayHello(params.greeting)
}
```

To bardzo minimalna definicja **workflow'a**.
W rzeczywistym pipeline'ie workflow zazwyczaj zawiera wiele wywołań **procesów** połączonych **kanałami**, i mogą być ustawione wartości domyślne dla zmiennych wejściowych.

Zobaczymy to w akcji, gdy uruchomimy nf-core/molkart w części 2 kursu.

### 2.4. System `params` parametrów wiersza poleceń

`params.greeting`, które przekazujemy do wywołania procesu `sayHello()`, to sprytny fragment kodu Nextflow i warto poświęcić mu dodatkową minutę.

Jak wspomniano powyżej, w ten sposób przekazujemy wartość parametru wiersza poleceń `--greeting` do wywołania procesu `sayHello()`.
W rzeczywistości samo zadeklarowanie `params.someParameterName` umożliwi nam podanie workflow'owi parametru o nazwie `--someParameterName` z wiersza poleceń.

!!! Tip "Wskazówka"

    Te parametry workflow'a zadeklarowane przy użyciu systemu `params` zawsze przyjmują dwie myślniki (`--`).
    Odróżnia to je od parametrów na poziomie Nextflow'a, które przyjmują tylko jeden myślnik (`-`).

### Podsumowanie

Teraz wiesz, jak zbudowany jest prosty workflow Nextflow i jak podstawowe komponenty odnoszą się do jego funkcjonalności.

### Co dalej?

Naucz się wygodnie zarządzać wykonaniami workflow'a.

---

## 3. Zarządzaj wykonaniami workflow'a

Wiedza o tym, jak uruchamiać workflow'y i pobierać wyjścia, jest świetna, ale szybko odkryjesz, że jest kilka innych aspektów zarządzania workflow'ami, które ułatwią Ci życie.

Tutaj pokażemy Ci, jak wykorzystać funkcję `resume` na wypadek, gdybyś musiał ponownie uruchomić ten sam workflow, jak sprawdzić logi wykonania za pomocą `nextflow log` oraz jak usunąć starsze katalogi robocze za pomocą `nextflow clean`.

### 3.1. Uruchom ponownie workflow'a z `-resume`

Czasami będziesz chciał ponownie uruchomić pipeline, który już wcześniej uruchomiłeś, bez powtarzania jakiejkolwiek pracy, która została już pomyślnie ukończona.

Nextflow ma opcję o nazwie `-resume`, która pozwala Ci to zrobić.
W szczególności w tym trybie wszystkie procesy, które zostały już uruchomione z dokładnie tym samym kodem, ustawieniami i wejściami, zostaną pominięte.
Oznacza to, że Nextflow uruchomi tylko procesy, które dodałeś lub zmodyfikowałeś od ostatniego uruchomienia, lub którym przekazujesz nowe ustawienia lub wejścia.

Istnieją dwie kluczowe zalety tego podejścia:

- Jeśli jesteś w trakcie tworzenia pipeline'a, możesz iterować szybciej, ponieważ musisz uruchomić tylko proces(y), nad którymi aktywnie pracujesz, aby przetestować swoje zmiany.
- Jeśli uruchamiasz pipeline w produkcji i coś pójdzie nie tak, w wielu przypadkach możesz naprawić problem i ponownie uruchomić pipeline, a on wznowi działanie od punktu awarii, co może zaoszczędzić Ci dużo czasu i mocy obliczeniowej.

Aby z tego skorzystać, po prostu dodaj `-resume` do swojego polecenia i uruchom je:

```bash
nextflow run hello-world.nf --greeting 'Hello World!' -resume
```

??? success "Wyjście polecenia"

    ```console
    N E X T F L O W   ~  version 25.04.3

    Launching `hello-world.nf` [tiny_noyce] DSL2 - revision: c33d41f479

    [a3/7be2fa] process > sayHello [100%] 1 of 1, cached: 1 ✔
    ```

Poszukaj fragmentu `cached:`, który został dodany w linii statusu procesu (linia 5), co oznacza, że Nextflow rozpoznał, że już wykonał tę pracę i po prostu ponownie wykorzystał wynik z poprzedniego pomyślnego uruchomienia.

Możesz również zobaczyć, że skrót podkatalogu roboczego jest taki sam jak w poprzednim uruchomieniu.
Nextflow dosłownie wskazuje Ci poprzednie wykonanie i mówi: „Już to zrobiłem tam".

!!! Tip "Wskazówka"

    Gdy ponownie uruchamiasz pipeline z `resume`, Nextflow nie nadpisuje żadnych plików zapisanych do katalogu `publishDir` przez żadne wywołanie procesu, które zostało wcześniej pomyślnie uruchomione.

### 3.2. Sprawdź log poprzednich wykonań

Za każdym razem, gdy uruchamiasz workflow Nextflow'a, linia jest zapisywana do pliku logu o nazwie `history`, w ukrytym katalogu o nazwie `.nextflow` w bieżącym katalogu roboczym.

Wygodniejszym sposobem dostępu do tych informacji jest użycie polecenia `nextflow log`.

```bash
nextflow log
```

Spowoduje to wyświetlenie zawartości pliku logu w terminalu, pokazując znacznik czasu, nazwę uruchomienia, status i pełny wiersz poleceń dla każdego uruchomienia Nextflow'a, które zostało uruchomione z bieżącego katalogu roboczego.

### 3.3. Usuń starsze katalogi robocze

Podczas procesu tworzenia zazwyczaj uruchomisz swoje wersje robocze pipeline'ów wiele razy, co może prowadzić do gromadzenia się bardzo wielu plików w wielu podkatalogach.
Ponieważ podkatalogi są nazywane losowo, trudno jest określić po ich nazwach, które są starsze, a które nowsze.

Nextflow zawiera wygodne podpolecenie `clean`, które może automatycznie usuwać podkatalogi robocze dla poprzednich uruchomień, które Cię już nie interesują, z kilkoma [opcjami](https://www.nextflow.io/docs/latest/reference/cli.html#clean) kontrolującymi, co zostanie usunięte.

Możesz użyć logu Nextflow'a, aby wyszukać uruchomienie na podstawie jego znacznika czasu i/lub wiersza poleceń, a następnie użyć `nextflow clean -before <run_name> -f`, aby usunąć katalogi robocze z wcześniejszych uruchomień.

!!! Warning "Ostrzeżenie"

    Usunięcie podkatalogów roboczych z poprzednich uruchomień usuwa je z pamięci podręcznej Nextflow'a i usuwa wszelkie wyjścia, które były przechowywane w tych katalogach.
    Oznacza to, że przerywa zdolność Nextflow'a do wznowienia wykonania bez ponownego uruchamiania odpowiednich procesów.

    Jesteś odpowiedzialny za zapisanie wszelkich wyjść, na których Ci zależy lub na których planujesz polegać! Jeśli używasz dyrektywy `publishDir` w tym celu, upewnij się, że używasz trybu `copy`, a nie trybu `symlink`.

### Podsumowanie

Wiesz, jak ponownie uruchomić pipeline bez powtarzania kroków, które zostały już uruchomione w identyczny sposób, sprawdzić log wykonania i użyć polecenia `nextflow clean` do czyszczenia starych katalogów roboczych.

### Co dalej?

Teraz, gdy rozumiesz podstawowe operacje Nextflow'a, jesteś gotowy do uruchomienia prawdziwego pipeline'a bioimagingu z nf-core/molkart.
