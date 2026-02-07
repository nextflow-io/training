# Część 5: Hello Containers

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Tłumaczenie wspomagane przez AI - [dowiedz się więcej i zasugeruj ulepszenia](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

<!--
<div class="video-wrapper">
  <iframe width="560" height="315" src="https://www.youtube.com/embed/5PyOWjKnNmg?si=QinuAnFwFj-Z8CrO&amp;list=PLPZ8WHdZGxmXiHf8B26oB_fTfoKQdhlik" title="YouTube video player" frameborder="0" allow="accelerometer; autoplay; clipboard-write; encrypted-media; gyroscope; picture-in-picture; web-share" referrerpolicy="strict-origin-when-cross-origin" allowfullscreen></iframe>
</div>

/// caption
:fontawesome-brands-youtube:{ .youtube } Obejrzyj [całą playlistę](https://www.youtube.com/playlist?list=PLPZ8WHdZGxmXiHf8B26oB_fTfoKQdhlik) na kanale YouTube Nextflow.

:green_book: Transkrypcja wideo jest dostępna [tutaj](./transcripts/05_hello_containers.md).
///
-->

W Częściach 1-4 tego kursu nauczyłeś się używać podstawowych elementów budulcowych Nextflow'a do składania prostego workflow'a zdolnego do przetwarzania tekstu, równoległego wykonywania wielu wejść i zbierania wyników do dalszego przetwarzania.

Jednak byłeś ograniczony do podstawowych narzędzi UNIX dostępnych w Twoim środowisku.
Rzeczywiste zadania często wymagają różnych narzędzi i pakietów, które nie są domyślnie dostępne.
Zazwyczaj musiałbyś zainstalować te narzędzia, zarządzać ich zależnościami i rozwiązywać konflikty.

To wszystko jest bardzo nużące i irytujące, więc pokażemy Ci, jak używać **kontenerów**, aby rozwiązać ten problem znacznie wygodniej.

**Kontener** to lekka, samodzielna, wykonywalna jednostka oprogramowania utworzona z **obrazu kontenera**, która zawiera wszystko potrzebne do uruchomienia aplikacji, w tym kod, biblioteki systemowe i ustawienia.
Jak można się domyślić, będzie to bardzo pomocne w zwiększeniu powtarzalności Twoich pipeline'ów.

Zauważ, że będziemy tego uczyć używając [Docker](https://www.docker.com/get-started/), ale pamiętaj, że Nextflow obsługuje również [kilka innych technologii kontenerowych](https://nextflow.io/docs/latest/container.html).

??? info "Jak zacząć od tej sekcji"

    Ta sekcja kursu zakłada, że ukończyłeś Części 1-4 kursu [Hello Nextflow](./index.md) i masz kompletny działający pipeline.

    Jeśli zaczynasz kurs od tego miejsca, musisz skopiować katalog `modules` z rozwiązań:

    ```bash
    cp -r solutions/4-hello-modules/modules .
    ```

---

## 0. Rozgrzewka: Uruchom `hello-containers.nf`

Użyjemy skryptu workflow `hello-containers.nf` jako punktu wyjścia.
Jest on równoważny skryptowi utworzonemu podczas pracy nad Częścią 4 tego szkolenia, z tą różnicą, że zmieniliśmy miejsca docelowe wyjść:

```groovy title="hello-containers.nf" linenums="37" hl_lines="3 7 11 15"
output {
    first_output {
        path 'hello_containers'
        mode 'copy'
    }
    uppercased {
        path 'hello_containers'
        mode 'copy'
    }
    collected {
        path 'hello_containers'
        mode 'copy'
    }
    batch_report {
        path 'hello_containers'
        mode 'copy'
    }
}
```

Aby upewnić się, że wszystko działa, uruchom skrypt raz przed wprowadzeniem jakichkolwiek zmian:

```bash
nextflow run hello-containers.nf
```

??? success "Wynik polecenia"

    ```console
     N E X T F L O W   ~  version 25.10.2

    Launching `hello-containers.nf` [nice_escher] DSL2 - revision: d5dfdc9872

    executor > local (7)
    [5a/ec1fa1] sayHello (2) [100%] 3 of 3 ✔
    [30/32b5b8] convertToUpper (3) [100%] 3 of 3 ✔
    [d3/be01bc] collectGreetings [100%] 1 of 1 ✔

    ```

Jak poprzednio, pliki wyjściowe znajdziesz w katalogu określonym w bloku `output` (`results/hello_containers/`).

??? abstract "Zawartość katalogu"

    ```console
    results/hello_containers/
    ├── Bonjour-output.txt
    ├── COLLECTED-batch-output.txt
    ├── Hello-output.txt
    ├── Holà-output.txt
    ├── batch-report.txt
    ├── UPPER-Bonjour-output.txt
    ├── UPPER-Hello-output.txt
    └── UPPER-Holà-output.txt
    ```

Jeśli to zadziałało, jesteś gotowy do nauki używania kontenerów.

---

## 1. Użyj kontenera 'ręcznie'

Chcemy dodać krok do naszego workflow'a, który będzie używał kontenera do wykonania.

Jednak najpierw omówimy podstawowe koncepcje i operacje, aby utrwalić Twoje zrozumienie tego, czym są kontenery, zanim zaczniemy ich używać w Nextflow'ie.

### 1.1. Pobierz obraz kontenera

Aby użyć kontenera, zazwyczaj pobierasz lub _ściągasz_ obraz kontenera z rejestru kontenerów, a następnie uruchamiasz go, aby utworzyć działającą instancję.

Ogólna składnia jest następująca:

```bash title="Składnia"
docker pull '<kontener>'
```

Część `docker pull` to instrukcja dla systemu kontenerowego, aby pobrał obraz kontenera z repozytorium.

Część `'<kontener>'` to adres URI obrazu kontenera.

Jako przykład pobierzmy obraz kontenera zawierający [cowpy](https://github.com/jeffbuttars/cowpy), pythonową implementację narzędzia o nazwie `cowsay`, które generuje grafikę ASCII do wyświetlania dowolnych tekstów wejściowych w zabawny sposób.

```txt title="Przykład"
 ________________________
< Are we having fun yet? >
 ------------------------
    \                                  ___-------___
     \                             _-~~             ~~-_
      \                         _-~                    /~-_
             /^\__/^\         /~  \                   /    \
           /|  O|| O|        /      \_______________/        \
          | |___||__|      /       /                \          \
          |          \    /      /                    \          \
          |   (_______) /______/                        \_________ \
          |         / /         \                      /            \
           \         \^\\         \                  /               \     /
             \         ||           \______________/      _-_       //\__//
               \       ||------_-~~-_ ------------- \ --/~   ~\    || __/
                 ~-----||====/~     |==================|       |/~~~~~
                  (_(__/  ./     /                    \_\      \.
                         (_(___/                         \_____)_)
```

Istnieją różne repozytoria, w których można znaleźć opublikowane kontenery.
Użyliśmy usługi [Seqera Containers](https://seqera.io/containers/) do wygenerowania tego obrazu kontenera Docker z pakietu Conda `cowpy`: `'community.wave.seqera.io/library/cowpy:1.1.5--3db457ae1977a273'`.

Uruchom pełne polecenie pobierania:

```bash
docker pull 'community.wave.seqera.io/library/cowpy:1.1.5--3db457ae1977a273'
```

??? success "Wynik polecenia"

    ```console
    1.1.5--3db457ae1977a273: Pulling from library/cowpy
    dafa2b0c44d2: Pull complete
    dec6b097362e: Pull complete
    f88da01cff0b: Pull complete
    4f4fb700ef54: Pull complete
    92dc97a3ef36: Pull complete
    403f74b0f85e: Pull complete
    10b8c00c10a5: Pull complete
    17dc7ea432cc: Pull complete
    bb36d6c3110d: Pull complete
    0ea1a16bbe82: Pull complete
    030a47592a0a: Pull complete
    c23bdb422167: Pull complete
    e1686ff32a11: Pull complete
    Digest: sha256:1ebc0043e8cafa61203bf42d29fd05bd14e7b4298e5e8cf986504c15f5aa4160
    Status: Downloaded newer image for community.wave.seqera.io/library/cowpy:1.1.5--3db457ae1977a273
    community.wave.seqera.io/library/cowpy:1.1.5--3db457ae1977a273
    ```

Jeśli nigdy wcześniej nie pobierałeś tego obrazu, może to zająć minutę.
Po zakończeniu masz lokalną kopię obrazu kontenera.

### 1.2. Użyj kontenera do uruchomienia `cowpy` jako jednorazowego polecenia

Bardzo częstym sposobem używania kontenerów jest uruchamianie ich bezpośrednio, czyli nieinteraktywnie.
To świetne rozwiązanie do uruchamiania jednorazowych poleceń.

Ogólna składnia jest następująca:

```bash title="Składnia"
docker run --rm '<kontener>' [polecenie narzędzia]
```

Część `docker run --rm '<kontener>'` to instrukcja dla systemu, aby uruchomił instancję kontenera z obrazu i wykonał w niej polecenie.
Flaga `--rm` mówi systemowi, aby usunął instancję kontenera po zakończeniu polecenia.

Składnia `[polecenie narzędzia]` zależy od używanego narzędzia i konfiguracji kontenera.
Zacznijmy po prostu od `cowpy`.

W pełni złożone polecenie wykonania kontenera wygląda tak; uruchom je:

```bash
docker run --rm 'community.wave.seqera.io/library/cowpy:1.1.5--3db457ae1977a273' cowpy
```

??? success "Wynik polecenia"

    ```console
    ______________________________________________________
    < Cowacter, eyes:default, tongue:False, thoughts:False >
    ------------------------------------------------------
        \   ^__^
          \  (oo)\_______
            (__)\       )\/\
              ||----w |
              ||     ||
    ```

System uruchomił kontener, wykonał polecenie `cowpy` z jego parametrami, wysłał wynik do konsoli i na koniec wyłączył instancję kontenera.

### 1.3. Użyj kontenera do uruchomienia `cowpy` interaktywnie

Możesz również uruchomić kontener interaktywnie, co daje Ci wiersz poleceń wewnątrz kontenera i pozwala eksperymentować z poleceniem.

#### 1.3.1. Uruchom kontener

Aby uruchomić interaktywnie, po prostu dodajemy `-it` do polecenia `docker run`.
Opcjonalnie możemy określić powłokę, której chcemy używać wewnątrz kontenera, dodając np. `/bin/bash` do polecenia.

```bash
docker run --rm -it 'community.wave.seqera.io/library/cowpy:1.1.5--3db457ae1977a273' /bin/bash
```

Zauważ, że Twój wiersz poleceń zmienia się na coś w rodzaju `(base) root@b645838b3314:/tmp#`, co wskazuje, że jesteś teraz wewnątrz kontenera.

Możesz to sprawdzić, uruchamiając `ls /`, aby wylistować zawartość katalogu od korzenia systemu plików:

```bash
ls /
```

??? abstract "Wynik polecenia"

    ```console
    bin  boot  dev  etc  home  lib  lib64  media  mnt  opt  proc  root  run  sbin  srv  sys  tmp  usr  var
    ```

Używamy tutaj `ls` zamiast `tree`, ponieważ narzędzie `tree` nie jest dostępne w tym kontenerze.
Widać, że system plików wewnątrz kontenera różni się od systemu plików na Twoim systemie hosta.

Jednym z ograniczeń tego, co właśnie zrobiliśmy, jest to, że kontener jest domyślnie całkowicie odizolowany od systemu hosta.
Oznacza to, że kontener nie może uzyskać dostępu do żadnych plików na systemie hosta, chyba że wyraźnie mu na to pozwolisz.

Pokażemy Ci, jak to zrobić za chwilę.

#### 1.3.2. Uruchom żądane polecenie(a) narzędzia

Teraz, gdy jesteś wewnątrz kontenera, możesz uruchomić polecenie `cowpy` bezpośrednio i podać mu parametry.
Na przykład dokumentacja narzędzia mówi, że możemy zmienić postać ('cowacter') za pomocą `-c`.

```bash
cowpy "Hello Containers" -c tux
```

??? success "Wynik polecenia"

    ```console
    __________________
    < Hello Containers >
    ------------------
      \
        \
            .--.
          |o_o |
          |:_/ |
          //   \ \
        (|     | )
        /'\_   _/`\
        \___)=(___/
    ```

Teraz wynik pokazuje pingwina Linuxa, Tuxa, zamiast domyślnej krowy, ponieważ określiliśmy parametr `-c tux`.

Ponieważ jesteś wewnątrz kontenera, możesz uruchamiać polecenie `cowpy` ile razy chcesz, zmieniając parametry wejściowe, bez konieczności zajmowania się poleceniami Docker.

!!! tip "Wskazówka"

    Użyj flagi '-c', aby wybrać inną postać, w tym:
    `beavis`, `cheese`, `daemon`, `dragonandcow`, `ghostbusters`, `kitty`, `moose`, `milk`, `stegosaurus`, `turkey`, `turtle`, `tux`

To fajne. Jeszcze fajniej byłoby, gdybyśmy mogli podać nasz plik `greetings.csv` jako wejście.
Ale ponieważ nie mamy dostępu do systemu plików, nie możemy.

Naprawmy to.

#### 1.3.3. Wyjdź z kontenera

Aby wyjść z kontenera, możesz wpisać `exit` w wierszu poleceń lub użyć skrótu klawiaturowego ++ctrl+d++.

```bash
exit
```

Twój wiersz poleceń powinien teraz wrócić do stanu sprzed uruchomienia kontenera.

#### 1.3.4. Zamontuj dane w kontenerze

Jak wspomniano wcześniej, kontener jest domyślnie odizolowany od systemu hosta.

Aby pozwolić kontenerowi na dostęp do systemu plików hosta, możesz **zamontować** **wolumen** z systemu hosta do kontenera, używając następującej składni:

```bash title="Składnia"
-v <ścieżka_zewnętrzna>:<ścieżka_wewnętrzna>
```

W naszym przypadku `<ścieżka_zewnętrzna>` będzie bieżącym katalogiem roboczym, więc możemy po prostu użyć kropki (`.`), a `<ścieżka_wewnętrzna>` to alias, który wymyślamy; nazwijmy go `/my_project` (ścieżka wewnętrzna musi być bezwzględna).

Aby zamontować wolumen, zastępujemy ścieżki i dodajemy argument montowania wolumenu do polecenia docker run w następujący sposób:

```bash
docker run --rm -it -v .:/my_project 'community.wave.seqera.io/library/cowpy:1.1.5--3db457ae1977a273' /bin/bash
```

To montuje bieżący katalog roboczy jako wolumen, który będzie dostępny pod `/my_project` wewnątrz kontenera.

Możesz sprawdzić, czy to działa, wylistowując zawartość `/my_project`:

```bash
ls /my_project
```

??? success "Wynik polecenia"

    ```console
    data               hello-config.nf      hello-modules.nf   hello-world.nf  nextflow.config  solutions         work
    hello-channels.nf  hello-containers.nf  hello-workflow.nf  modules         results          test-params.json
    ```

Teraz możesz zobaczyć zawartość katalogu roboczego z wnętrza kontenera, w tym plik `greetings.csv` w katalogu `data/`.

To skutecznie utworzyło tunel przez ścianę kontenera, którego możesz użyć do dostępu do tej części systemu plików.

#### 1.3.5. Użyj zamontowanych danych

Teraz, gdy zamontowaliśmy katalog roboczy w kontenerze, możemy użyć polecenia `cowpy` do wyświetlenia zawartości pliku `greetings.csv`.

W tym celu użyjemy `cat /my_project/data/greetings.csv | `, aby przekierować zawartość pliku CSV do polecenia `cowpy`.

```bash
cat /my_project/data/greetings.csv | cowpy -c turkey
```

??? success "Wynik polecenia"

    ```console title="data/greetings.csv"
     ____________________
    / Hello,English,123  \
    | Bonjour,French,456 |
    \ Holà,Spanish,789   /
    --------------------
      \                                  ,+*^^*+___+++_
      \                           ,*^^^^              )
        \                       _+*                     ^**+_
        \                    +^       _ _++*+_+++_,         )
                  _+^^*+_    (     ,+*^ ^          \+_        )
                {       )  (    ,(    ,_+--+--,      ^)      ^\
                { (\@)    } f   ,(  ,+-^ __*_*_  ^^\_   ^\       )
              {:;-/    (_+*-+^^^^^+*+*<_ _++_)_    )    )      /
              ( /  (    (        ,___    ^*+_+* )   <    <      \
              U _/     )    *--<  ) ^\-----++__)   )    )       )
                (      )  _(^)^^))  )  )\^^^^^))^*+/    /       /
              (      /  (_))_^)) )  )  ))^^^^^))^^^)__/     +^^
            (     ,/    (^))^))  )  ) ))^^^^^^^))^^)       _)
              *+__+*       (_))^)  ) ) ))^^^^^^))^^^^^)____*^
              \             \_)^)_)) ))^^^^^^^^^^))^^^^)
              (_             ^\__^^^^^^^^^^^^))^^^^^^^)
                ^\___            ^\__^^^^^^))^^^^^^^^)\\
                      ^^^^^\uuu/^^\uuu/^^^^\^\^\^\^\^\^\^\
                        ___) >____) >___   ^\_\_\_\_\_\_\)
                        ^^^//\\_^^//\\_^       ^(\_\_\_\)
                          ^^^ ^^ ^^^ ^
    ```

To produkuje pożądaną grafikę ASCII indyka recytującego nasze przykładowe pozdrowienia!
Tyle że tutaj indyk powtarza pełne wiersze zamiast samych pozdrowień.
Wiemy już, że nasz workflow Nextflow'a zrobi to lepiej!

Poeksperymentuj z tym poleceniem.
Gdy skończysz, wyjdź z kontenera jak poprzednio:

```bash
exit
```

Wrócisz do normalnej powłoki.

### Podsumowanie

Wiesz już, jak pobrać kontener i uruchomić go jednorazowo lub interaktywnie. Wiesz również, jak udostępnić Swoje dane z wnętrza kontenera, co pozwala wypróbować każde interesujące Cię narzędzie na prawdziwych danych bez konieczności instalowania oprogramowania na Swoim systemie.

### Co dalej?

Naucz się używać kontenerów do wykonywania procesów Nextflow'a.

---

## 2. Używaj kontenerów w Nextflow'ie

Nextflow ma wbudowaną obsługę uruchamiania procesów wewnątrz kontenerów, co pozwala korzystać z narzędzi, których nie masz zainstalowanych w swoim środowisku obliczeniowym.
Oznacza to, że możesz użyć dowolnego obrazu kontenera do wykonywania procesów, a Nextflow zajmie się pobieraniem obrazu, montowaniem danych i obsługą całego cyklu życia kontenera.

Aby to zademonstrować, dodamy krok `cowpy` do pipeline'u, który rozwijaliśmy, po kroku `collectGreetings`.

<figure class="excalidraw">
--8<-- "docs/en/docs/hello_nextflow/img/hello-pipeline-cowpy.svg"
</figure>

### 2.1. Napisz moduł `cowpy`

Najpierw utwórzmy moduł procesu `cowpy`.

#### 2.1.1. Utwórz plik dla nowego modułu

Utwórz pusty plik dla modułu o nazwie `cowpy.nf`.

```bash
touch modules/cowpy.nf
```

To daje nam miejsce na umieszczenie kodu procesu.

#### 2.1.2. Skopiuj kod procesu `cowpy` do pliku modułu

Możemy wzorować nasz proces `cowpy` na innych procesach, które napisaliśmy wcześniej.

```groovy title="modules/cowpy.nf" linenums="1"
#!/usr/bin/env nextflow

// Wygeneruj grafikę ASCII za pomocą cowpy
process cowpy {

    input:
    path input_file
    val character

    output:
    path "cowpy-${input_file}"

    script:
    """
    cat ${input_file} | cowpy -c "${character}" > cowpy-${input_file}
    """

}
```

Proces oczekuje `input_file` zawierającego pozdrowienia oraz wartości `character`.

Wyjściem będzie nowy plik tekstowy zawierający grafikę ASCII wygenerowaną przez narzędzie `cowpy`.

### 2.2. Dodaj cowpy do workflow'a

Teraz musimy zaimportować moduł i wywołać proces.

#### 2.2.1. Importuj proces `cowpy` do `hello-containers.nf`

Wstaw deklarację importu powyżej bloku workflow i wypełnij ją odpowiednio.

=== "Po"

    ```groovy title="hello-containers.nf" linenums="3" hl_lines="5"
    // Dołącz moduły
    include { sayHello } from './modules/sayHello.nf'
    include { convertToUpper } from './modules/convertToUpper.nf'
    include { collectGreetings } from './modules/collectGreetings.nf'
    include { cowpy } from './modules/cowpy.nf'
    ```

=== "Przed"

    ```groovy title="hello-containers.nf" linenums="3"
    // Dołącz moduły
    include { sayHello } from './modules/sayHello.nf'
    include { convertToUpper } from './modules/convertToUpper.nf'
    include { collectGreetings } from './modules/collectGreetings.nf'
    ```

Teraz moduł `cowpy` jest dostępny do użycia w workflow'ie.

#### 2.2.2. Dodaj wywołanie procesu `cowpy` w workflow'ie

Połączmy proces `cowpy()` z wyjściem procesu `collectGreetings()`, który jak pamiętasz produkuje dwa wyjścia:

- `collectGreetings.out.outfile` zawiera plik wyjściowy <--_to, czego potrzebujemy_
- `collectGreetings.out.report` zawiera plik raportu z liczbą pozdrowień na partię

W bloku workflow wprowadź następującą zmianę w kodzie:

=== "Po"

    ```groovy title="hello-containers.nf" linenums="19" hl_lines="12-13"
        main:
        // utwórz kanał dla danych wejściowych z pliku CSV
        greeting_ch = channel.fromPath(params.input)
                            .splitCsv()
                            .map { line -> line[0] }
        // wyemituj pozdrowienie
        sayHello(greeting_ch)
        // przekształć pozdrowienie na wielkie litery
        convertToUpper(sayHello.out)
        // zbierz wszystkie pozdrowienia do jednego pliku
        collectGreetings(convertToUpper.out.collect(), params.batch)
        // wygeneruj grafikę ASCII pozdrowień za pomocą cowpy
        cowpy(collectGreetings.out.outfile, params.character)
    ```

=== "Przed"

    ```groovy title="hello-containers.nf" linenums="19"
        main:
        // utwórz kanał dla danych wejściowych z pliku CSV
        greeting_ch = channel.fromPath(params.input)
                            .splitCsv()
                            .map { line -> line[0] }
        // wyemituj pozdrowienie
        sayHello(greeting_ch)
        // przekształć pozdrowienie na wielkie litery
        convertToUpper(sayHello.out)
        // zbierz wszystkie pozdrowienia do jednego pliku
        collectGreetings(convertToUpper.out.collect(), params.batch)
    ```

Zauważ, że zadeklarowaliśmy nowy parametr CLI, `params.character`, aby określić, która postać ma wypowiedzieć pozdrowienia.

#### 2.2.3. Dodaj parametr `character` do bloku `params`

To technicznie jest opcjonalne, ale jest to zalecana praktyka i okazja do ustawienia domyślnej wartości dla postaci.

=== "Po"

    ```groovy title="hello-containers.nf" linenums="9" hl_lines="7"
    /*
    * Parametry pipeline'u
    */
    params {
        input: Path = 'data/greetings.csv'
        batch: String = 'batch'
        character: String = 'turkey'
    }
    ```

=== "Przed"

    ```groovy title="hello-containers.nf" linenums="9"
    /*
    * Parametry pipeline'u
    */
    params {
        input: Path = 'data/greetings.csv'
        batch: String = 'batch'
    }
    ```

Teraz możemy być leniwi i pominąć wpisywanie parametru postaci w naszych poleceniach.

#### 2.2.4. Zaktualizuj wyjścia workflow'a

Musimy zaktualizować wyjścia workflow'a, aby publikować wyjście procesu `cowpy`.

##### 2.2.4.1. Zaktualizuj sekcję `publish:`

W bloku `workflow` wprowadź następującą zmianę w kodzie:

=== "Po"

    ```groovy title="hello-containers.nf" linenums="34" hl_lines="6"
        publish:
        first_output = sayHello.out
        uppercased = convertToUpper.out
        collected = collectGreetings.out.outfile
        batch_report = collectGreetings.out.report
        cowpy_art = cowpy.out
    ```

=== "Przed"

    ```groovy title="hello-containers.nf" linenums="34"
        publish:
        first_output = sayHello.out
        uppercased = convertToUpper.out
        collected = collectGreetings.out.outfile
        batch_report = collectGreetings.out.report
    ```

Proces `cowpy` produkuje tylko jedno wyjście, więc możemy się do niego odnieść w zwykły sposób, dodając `.out`.

Ale na razie dokończmy aktualizację wyjść na poziomie workflow'a.

##### 2.2.4.2. Zaktualizuj blok `output`

Musimy dodać końcowe wyjście `cowpy_art` do bloku `output`. Przy okazji edytujmy również miejsca docelowe publikacji, ponieważ teraz nasz pipeline jest kompletny i wiemy, które wyjścia naprawdę nas interesują.

W bloku `output` wprowadź następujące zmiany w kodzie:

=== "Po"

    ```groovy title="hello-containers.nf" linenums="42" hl_lines="3 7 11 15 18-21"
    output {
        first_output {
            path 'hello_containers/intermediates'
            mode 'copy'
        }
        uppercased {
            path 'hello_containers/intermediates'
            mode 'copy'
        }
        collected {
            path 'hello_containers/intermediates'
            mode 'copy'
        }
        batch_report {
            path 'hello_containers'
            mode 'copy'
        }
        cowpy_art {
            path 'hello_containers'
            mode 'copy'
        }
    }
    ```

=== "Przed"

    ```groovy title="hello-containers.nf" linenums="42" hl_lines="3 7 11 15"
    output {
        first_output {
            path 'hello_containers'
            mode 'copy'
        }
        uppercased {
            path 'hello_containers'
            mode 'copy'
        }
        collected {
            path 'hello_containers'
            mode 'copy'
        }
        batch_report {
            path 'hello_containers'
            mode 'copy'
        }
    }
    ```

Teraz opublikowane wyjścia będą nieco lepiej zorganizowane.

#### 2.2.5. Uruchom workflow

Podsumowując, oto co chcemy osiągnąć:

<figure class="excalidraw">
--8<-- "docs/en/docs/hello_nextflow/img/hello_pipeline_complete.svg"
</figure>

Myślisz, że zadziała?

Usuńmy poprzednie opublikowane wyjścia, aby mieć czystą kartę, i uruchommy workflow z flagą `-resume`.

```bash
rm -r hello_containers/
nextflow run hello-containers.nf -resume
```

??? failure "Wynik polecenia (zredagowany dla przejrzystości)"

    ```console hl_lines="10 13 20-21 26-27"
     N E X T F L O W   ~  version 25.10.2

    Launching `hello-containers.nf` [lonely_woese] DSL2 - revision: abf1dccf7f

    executor >  local (1)
    [c9/f5c686] sayHello (3)       [100%] 3 of 3, cached: 3 ✔
    [ef/3135a8] convertToUpper (3) [100%] 3 of 3, cached: 3 ✔
    [7f/f435e3] collectGreetings   [100%] 1 of 1, cached: 1 ✔
    [9b/02e776] cowpy              [  0%] 0 of 1 ✘
    ERROR ~ Error executing process > 'cowpy'

    Caused by:
      Process `cowpy` terminated with an error exit status (127)


    Command executed:

      cat COLLECTED-batch-output.txt | cowpy -c "turkey" > cowpy-COLLECTED-batch-output.txt

    Command exit status:
      127

    Command output:
      (empty)

    Command error:
      .command.sh: line 2: cowpy: command not found

    Work dir:
      /workspaces/training/hello-nextflow/work/9b/02e7761db848f82db3c3e59ff3a9b6

    Tip: when you have fixed the problem you can continue the execution adding the option `-resume` to the run command line

    -- Check '.nextflow.log' file for details
    ERROR ~ Cannot access first() element from an empty List

    -- Check '.nextflow.log' file for details
    ```

O nie, jest błąd!
Kod błędu podany przez `error exit status (127)` oznacza, że żądany plik wykonywalny nie został znaleziony.

To ma sens, ponieważ wywołujemy narzędzie `cowpy`, ale tak naprawdę nie określiliśmy jeszcze kontenera (ups).

### 2.3. Użyj kontenera do uruchomienia procesu `cowpy`

Musimy określić kontener i powiedzieć Nextflow'owi, aby użył go dla procesu `cowpy()`.

#### 2.3.1. Określ kontener dla `cowpy`

Możemy użyć tego samego obrazu, którego używaliśmy bezpośrednio w pierwszej sekcji tego samouczka.

Edytuj moduł `cowpy.nf`, aby dodać dyrektywę `container` do definicji procesu w następujący sposób:

=== "Po"

    ```groovy title="modules/cowpy.nf" linenums="4" hl_lines="3"
    process cowpy {

        container 'community.wave.seqera.io/library/cowpy:1.1.5--3db457ae1977a273'

        input:
        path input_file
        val character

        output:
        path "cowpy-${input_file}"

        script:
        """
        cat ${input_file} | cowpy -c "${character}" > cowpy-${input_file}
        """
    }
    ```

=== "Przed"

    ```groovy title="modules/cowpy.nf" linenums="4"
    process cowpy {

        input:
        path input_file
        val character

        output:
        path "cowpy-${input_file}"

        script:
        """
        cat ${input_file} | cowpy -c "${character}" > cowpy-${input_file}
        """
    }
    ```

To mówi Nextflow'owi, że _jeśli użycie Docker jest włączone_, powinien użyć określonego tutaj obrazu kontenera do wykonania procesu.

#### 2.3.2. Włącz użycie Docker przez plik `nextflow.config`

Zauważ, że wspomnieliśmy _'jeśli użycie Docker jest włączone'_. Domyślnie nie jest, więc musimy włączyć tę opcję w Nextflow'ie.
W tym celu nieco uprzedzamy temat następnej i ostatniej części tego kursu (Część 6), która obejmuje konfigurację.

Jednym z głównych sposobów konfigurowania wykonywania workflow'a, które oferuje Nextflow, jest użycie pliku `nextflow.config`.
Gdy taki plik jest obecny w bieżącym katalogu, Nextflow automatycznie go załaduje i zastosuje zawartą w nim konfigurację.

Dostarczyliśmy plik `nextflow.config` z jedną linią kodu, która jawnie wyłącza Docker: `docker.enabled = false`.

Teraz zmieńmy to na `true`, aby włączyć Docker:

=== "Po"

    ```console title="nextflow.config" linenums="1" hl_lines="1"
    docker.enabled = true
    ```

=== "Przed"

    ```console title="nextflow.config" linenums="1" hl_lines="1"
    docker.enabled = false
    ```

!!! tip "Wskazówka"

    Możliwe jest włączenie wykonywania Docker z wiersza poleceń, dla pojedynczego uruchomienia, używając parametru `-with-docker <kontener>`.
    Jednak pozwala to tylko na określenie jednego kontenera dla całego workflow'a, podczas gdy podejście, które właśnie pokazaliśmy, pozwala na określenie innego kontenera dla każdego procesu.
    To jest lepsze dla modularności, utrzymania kodu i powtarzalności.

#### 2.3.3. Uruchom workflow z włączonym Docker

Uruchom workflow z flagą `-resume`:

```bash
nextflow run hello-containers.nf -resume
```

??? success "Wynik polecenia"

    ```console
     N E X T F L O W   ~  version 25.10.2

    Launching `hello-containers.nf` [drunk_perlman] DSL2 - revision: abf1dccf7f

    executor >  local (1)
    [c9/f5c686] sayHello (3)       [100%] 3 of 3, cached: 3 ✔
    [ef/3135a8] convertToUpper (3) [100%] 3 of 3, cached: 3 ✔
    [7f/f435e3] collectGreetings   [100%] 1 of 1, cached: 1 ✔
    [98/656c6c] cowpy              [100%] 1 of 1 ✔
    ```

Tym razem rzeczywiście działa!
Jak zwykle wyjścia workflow'a znajdziesz w odpowiednim katalogu results, choć tym razem są nieco lepiej zorganizowane, z tylko raportem i końcowym wyjściem na najwyższym poziomie, a wszystkie pliki pośrednie schowane w podkatalogu.

??? abstract "Zawartość katalogu"

    ```console
    results/hello_containers/
    ├── cowpy-COLLECTED-batch-output.txt
    ├── intermediates
    │   ├── Bonjour-output.txt
    │   ├── COLLECTED-batch-output.txt
    │   ├── Hello-output.txt
    │   ├── Holà-output.txt
    │   ├── UPPER-Bonjour-output.txt
    │   ├── UPPER-Hello-output.txt
    │   └── UPPER-Holà-output.txt
    └── batch-report.txt
    ```

Końcowa grafika ASCII znajduje się w katalogu `results/hello_containers/`, pod nazwą `cowpy-COLLECTED-batch-output.txt`.

??? abstract "Zawartość pliku"

    ```console title="results/hello_containers/cowpy-COLLECTED-batch-output.txt"
    _________
    / HOLà    \
    | HELLO   |
    \ BONJOUR /
    ---------
      \                                  ,+*^^*+___+++_
      \                           ,*^^^^              )
        \                       _+*                     ^**+_
        \                    +^       _ _++*+_+++_,         )
                  _+^^*+_    (     ,+*^ ^          \+_        )
                {       )  (    ,(    ,_+--+--,      ^)      ^\
                { (\@)    } f   ,(  ,+-^ __*_*_  ^^\_   ^\       )
              {:;-/    (_+*-+^^^^^+*+*<_ _++_)_    )    )      /
              ( /  (    (        ,___    ^*+_+* )   <    <      \
              U _/     )    *--<  ) ^\-----++__)   )    )       )
                (      )  _(^)^^))  )  )\^^^^^))^*+/    /       /
              (      /  (_))_^)) )  )  ))^^^^^))^^^)__/     +^^
            (     ,/    (^))^))  )  ) ))^^^^^^^))^^)       _)
              *+__+*       (_))^)  ) ) ))^^^^^^))^^^^^)____*^
              \             \_)^)_)) ))^^^^^^^^^^))^^^^)
              (_             ^\__^^^^^^^^^^^^))^^^^^^^)
                ^\___            ^\__^^^^^^))^^^^^^^^)\\
                      ^^^^^\uuu/^^\uuu/^^^^\^\^\^\^\^\^\^\
                        ___) >____) >___   ^\_\_\_\_\_\_\)
                        ^^^//\\_^^//\\_^       ^(\_\_\_\)
                          ^^^ ^^ ^^^ ^
    ```

I oto jest, nasz piękny indyk wypowiadający pozdrowienia zgodnie z oczekiwaniami.

#### 2.3.4. Sprawdź, jak Nextflow uruchomił skonteneryzowane zadanie

Na zakończenie tej sekcji przyjrzyjmy się podkatalogowi work dla jednego z wywołań procesu `cowpy`, aby uzyskać trochę więcej wglądu w to, jak Nextflow pracuje z kontenerami pod maską.

Sprawdź wynik polecenia `nextflow run`, aby znaleźć ścieżkę do podkatalogu work dla procesu `cowpy`.
Patrząc na to, co otrzymaliśmy dla pokazanego powyżej uruchomienia, linia logu konsoli dla procesu `cowpy` zaczyna się od `[98/656c6c]`.
Odpowiada to następującej skróconej ścieżce katalogu: `work/98/656c6c`.

W tym katalogu znajdziesz plik `.command.run`, który zawiera wszystkie polecenia, które Nextflow uruchomił w Twoim imieniu podczas wykonywania pipeline'u.

??? abstract "Zawartość pliku"

    ```console title="work/98/656c6c90cce1667c094d880f4b6dcc/.command.run"
    #!/bin/bash
    ### ---
    ### name: 'cowpy'
    ### container: 'community.wave.seqera.io/library/cowpy:1.1.5--3db457ae1977a273'
    ### outputs:
    ### - 'cowpy-COLLECTED-batch-output.txt'
    ### ...
    set -e
    set -u
    NXF_DEBUG=${NXF_DEBUG:=0}; [[ $NXF_DEBUG > 1 ]] && set -x
    NXF_ENTRY=${1:-nxf_main}


    nxf_sleep() {
      sleep $1 2>/dev/null || sleep 1;
    }

    nxf_date() {
        local ts=$(date +%s%3N);
        if [[ ${#ts} == 10 ]]; then echo ${ts}000
        elif [[ $ts == *%3N ]]; then echo ${ts/\%3N/000}
        elif [[ $ts == *3N ]]; then echo ${ts/3N/000}
        elif [[ ${#ts} == 13 ]]; then echo $ts
        else echo "Unexpected timestamp value: $ts"; exit 1
        fi
    }

    nxf_env() {
        echo '============= task environment ============='
        env | sort | sed "s/\(.*\)AWS\(.*\)=\(.\{6\}\).*/\1AWS\2=\3xxxxxxxxxxxxx/"
        echo '============= task output =================='
    }

    nxf_kill() {
        declare -a children
        while read P PP;do
            children[$PP]+=" $P"
        done < <(ps -e -o pid= -o ppid=)

        kill_all() {
            [[ $1 != $$ ]] && kill $1 2>/dev/null || true
            for i in ${children[$1]:=}; do kill_all $i; done
        }

        kill_all $1
    }

    nxf_mktemp() {
        local base=${1:-/tmp}
        mkdir -p "$base"
        if [[ $(uname) = Darwin ]]; then mktemp -d $base/nxf.XXXXXXXXXX
        else TMPDIR="$base" mktemp -d -t nxf.XXXXXXXXXX
        fi
    }

    nxf_fs_copy() {
      local source=$1
      local target=$2
      local basedir=$(dirname $1)
      mkdir -p $target/$basedir
      cp -fRL $source $target/$basedir
    }

    nxf_fs_move() {
      local source=$1
      local target=$2
      local basedir=$(dirname $1)
      mkdir -p $target/$basedir
      mv -f $source $target/$basedir
    }

    nxf_fs_rsync() {
      rsync -rRl $1 $2
    }

    nxf_fs_rclone() {
      rclone copyto $1 $2/$1
    }

    nxf_fs_fcp() {
      fcp $1 $2/$1
    }

    on_exit() {
        local last_err=$?
        local exit_status=${nxf_main_ret:=0}
        [[ ${exit_status} -eq 0 && ${nxf_unstage_ret:=0} -ne 0 ]] && exit_status=${nxf_unstage_ret:=0}
        [[ ${exit_status} -eq 0 && ${last_err} -ne 0 ]] && exit_status=${last_err}
        printf -- $exit_status > /workspaces/training/hello-nextflow/work/98/656c6c90cce1667c094d880f4b6dcc/.exitcode
        set +u
        docker rm $NXF_BOXID &>/dev/null || true
        exit $exit_status
    }

    on_term() {
        set +e
        docker stop $NXF_BOXID
    }

    nxf_launch() {
        docker run -i --cpu-shares 1024 -e "NXF_TASK_WORKDIR" -v /workspaces/training/hello-nextflow/work:/workspaces/training/hello-nextflow/work -w "$NXF_TASK_WORKDIR" --name $NXF_BOXID community.wave.seqera.io/library/cowpy:1.1.5--3db457ae1977a273 /bin/bash -ue /workspaces/training/hello-nextflow/work/98/656c6c90cce1667c094d880f4b6dcc/.command.sh
    }

    nxf_stage() {
        true
        # stage input files
        rm -f COLLECTED-batch-output.txt
        ln -s /workspaces/training/hello-nextflow/work/7f/f435e3f2cf95979b5f3d7647ae6696/COLLECTED-batch-output.txt COLLECTED-batch-output.txt
    }

    nxf_unstage_outputs() {
        true
    }

    nxf_unstage_controls() {
        true
    }

    nxf_unstage() {
        if [[ ${nxf_main_ret:=0} == 0 ]]; then
            (set -e -o pipefail; (nxf_unstage_outputs | tee -a .command.out) 3>&1 1>&2 2>&3 | tee -a .command.err)
            nxf_unstage_ret=$?
        fi
        nxf_unstage_controls
    }

    nxf_main() {
        trap on_exit EXIT
        trap on_term TERM INT USR2
        trap '' USR1

        [[ "${NXF_CHDIR:-}" ]] && cd "$NXF_CHDIR"
        export NXF_BOXID="nxf-$(dd bs=18 count=1 if=/dev/urandom 2>/dev/null | base64 | tr +/ 0A | tr -d '\r\n')"
        NXF_SCRATCH=''
        [[ $NXF_DEBUG > 0 ]] && nxf_env
        touch /workspaces/training/hello-nextflow/work/98/656c6c90cce1667c094d880f4b6dcc/.command.begin
        set +u
        set -u
        [[ $NXF_SCRATCH ]] && cd $NXF_SCRATCH
        export NXF_TASK_WORKDIR="$PWD"
        nxf_stage

        set +e
        (set -o pipefail; (nxf_launch | tee .command.out) 3>&1 1>&2 2>&3 | tee .command.err) &
        pid=$!
        wait $pid || nxf_main_ret=$?
        nxf_unstage
    }

    $NXF_ENTRY

    ```

Jeśli wyszukasz `nxf_launch` w tym pliku, powinieneś zobaczyć coś takiego:

```console
nxf_launch() {
    docker run -i --cpu-shares 1024 -e "NXF_TASK_WORKDIR" -v /workspaces/training/hello-nextflow/work:/workspaces/training/hello-nextflow/work -w "$NXF_TASK_WORKDIR" --name $NXF_BOXID community.wave.seqera.io/library/cowpy:1.1.5--3db457ae1977a273 /bin/bash -ue /workspaces/training/hello-nextflow/work/98/656c6c90cce1667c094d880f4b6dcc/.command.sh
}
```

Jak widać, Nextflow używa polecenia `docker run` do uruchomienia wywołania procesu.
Montuje również odpowiedni podkatalog work do kontenera, ustawia odpowiednio katalog roboczy wewnątrz kontenera i uruchamia nasz szablonowy skrypt bash w pliku `.command.sh`.

Cała ciężka praca, którą musieliśmy wykonać ręcznie w pierwszej sekcji? Nextflow robi to za nas za kulisami!

```txt
 _______________________
< Hurray for robots...! >
 -----------------------
                                   ,-----.
                                   |     |
                                ,--|     |-.
                         __,----|  |     | |
                       ,;::     |  `_____' |
                       `._______|    i^i   |
                                `----| |---'| .
                           ,-------._| |== ||//
                           |       |_|P`.  /'/
                           `-------' 'Y Y/'/'
                                     .==\ /_\
   ^__^                             /   /'|  `i
   (oo)\_______                   /'   /  |   |
   (__)\       )\/\             /'    /   |   `i
       ||----w |           ___,;`----'.___L_,-'`\__
       ||     ||          i_____;----\.____i""\____\
```

### Podsumowanie

Wiesz już, jak używać kontenerów w Nextflow'ie do uruchamiania procesów.

### Co dalej?

Zrób przerwę!

Gdy będziesz gotowy, przejdź do [**Część 6: Hello Config**](./06_hello_config.md), aby dowiedzieć się, jak konfigurować wykonywanie pipeline'u, aby dopasować go do Twojej infrastruktury, a także zarządzać konfiguracją wejść i parametrów.

To ostatnia część i potem ukończysz ten kurs!

---

## Quiz

<quiz>
Czym jest kontener?
- [ ] Rodzaj maszyny wirtualnej
- [ ] Format kompresji plików
- [x] Lekka, samodzielna, wykonywalna jednostka zawierająca wszystko potrzebne do uruchomienia aplikacji
- [ ] Protokół sieciowy
</quiz>

<quiz>
Jaka jest różnica między obrazem kontenera a instancją kontenera?
- [ ] To to samo
- [x] Obraz to szablon; instancja to działający kontener utworzony z tego obrazu
- [ ] Instancja to szablon; obraz to działający kontener
- [ ] Obrazy są dla Docker; instancje są dla Singularity
</quiz>

<quiz>
Co robi flaga `-v` w poleceniu `docker run`?
- [ ] Włącza szczegółowe wyjście
- [ ] Waliduje kontener
- [x] Montuje wolumen z systemu hosta do kontenera
- [ ] Określa wersję kontenera

Dowiedz się więcej: [1.3.4. Zamontuj dane w kontenerze](#134-zamontuj-dane-w-kontenerze)
</quiz>

<quiz>
Dlaczego trzeba montować wolumeny podczas używania kontenerów?
- [ ] Aby poprawić wydajność kontenera
- [ ] Aby oszczędzić miejsce na dysku
- [x] Ponieważ kontenery są domyślnie odizolowane od systemu plików hosta
- [ ] Aby włączyć sieć

Dowiedz się więcej: [1.3.4. Zamontuj dane w kontenerze](#134-zamontuj-dane-w-kontenerze)
</quiz>

<quiz>
Jak określasz kontener dla procesu Nextflow'a?
- [ ] `docker 'container-uri'`
- [ ] `image 'container-uri'`
- [x] `container 'container-uri'`
- [ ] `use 'container-uri'`

Dowiedz się więcej: [2.3.1. Określ kontener dla cowpy](#231-okresl-kontener-dla-cowpy)
</quiz>

<quiz>
Które ustawienie `nextflow.config` włącza Docker dla Twojego workflow'a?
- [ ] `#!groovy process.docker = true`
- [x] `#!groovy docker.enabled = true`
- [ ] `#!groovy container.engine = 'docker'`
- [ ] `#!groovy docker.activate = true`

Dowiedz się więcej: [2.3.2. Włącz użycie Docker przez plik `nextflow.config`](#232-wlacz-uzycie-docker-przez-plik-nextflowconfig)
</quiz>

<quiz>
Co Nextflow automatycznie obsługuje podczas uruchamiania procesu w kontenerze? (Wybierz wszystkie pasujące)
- [x] Pobieranie obrazu kontenera w razie potrzeby
- [x] Montowanie katalogu work
- [x] Uruchamianie skryptu procesu wewnątrz kontenera
- [x] Czyszczenie instancji kontenera po wykonaniu

Dowiedz się więcej: [2.3.4. Sprawdź, jak Nextflow uruchomił skonteneryzowane zadanie](#234-sprawdz-jak-nextflow-uruchomil-skonteneryzowane-zadanie)
</quiz>
