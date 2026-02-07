# Część 3: Hello Workflow - Transkrypcja

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Tłumaczenie wspomagane przez AI - [dowiedz się więcej i zasugeruj ulepszenia](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

<div class="video-wrapper">
  <iframe width="560" height="315" src="https://www.youtube.com/embed/zJP7cUYPEbA?si=Irl9nAQniDyICp2b&amp;list=PLPZ8WHdZGxmXiHf8B26oB_fTfoKQdhlik" title="YouTube video player" frameborder="0" allow="accelerometer; autoplay; clipboard-write; encrypted-media; gyroscope; picture-in-picture; web-share" referrerpolicy="strict-origin-when-cross-origin" allowfullscreen></iframe>
</div>

!!!note "Ważne uwagi"

    Ta strona zawiera tylko transkrypcję. Aby uzyskać pełne instrukcje krok po kroku, wróć do [materiałów szkoleniowych](../03_hello_workflow.md).

    Numery sekcji pokazane w transkrypcji są podane tylko orientacyjnie i mogą nie zawierać wszystkich numerów sekcji w materiałach.

## Powitanie

Witam, witajcie w trzeciej części kursu szkoleniowego "Hello Nextflow".

Ten rozdział nosi tytuł "Hello Workflow".

W rozdziale drugim zbudowaliśmy prosty workflow składający się z jednego procesu, ale w rzeczywistości pipeline'y są użyteczne, ponieważ mogą łączyć wiele kroków analizy razem.

W tym rozdziale weźmiemy ten początkowy przykład i rozszerzymy go, aby był nieco bardziej realistyczny.

Dodamy kilka dodatkowych kroków i przyjrzymy się, jak używamy kanałów do łączenia tych kroków.

Będziemy się zajmować wieloma zadaniami, które mogą się zwijać do jednego procesu, i przyjrzymy się procesom, które mogą mieć wiele wejść i wiele wyjść.

Dobra, zaczynajmy.

Zacznijmy. Tak jak poprzednio: przejdźmy do training.nextflow.io, Hello Nextflow, rozdział trzeci, Hello Workflow — i otwórzmy nasz obszar roboczy. Uporządkowałem wszystkie pliki robocze z poprzednich rozdziałów i teraz otwieram Hello Workflow.

To jest ten sam plik, nad którym pracowaliśmy do tej pory, więc powinien wyglądać znajomo. Mamy nasz proces `SAY_HELLO`. Mamy nasz `params.greeting` z plikiem CSV greetings i mamy nasz workflow na dole, który ładuje ten plik CSV, tworzy kanał i przekazuje go do naszego procesu.

## 0. Rozgrzewka: Uruchomienie hello-workflow.nf

Jeśli chcesz, możemy to wypróbować i dokładnie sprawdzić, czy działa zgodnie z oczekiwaniami. Otwórzmy terminal i uruchommy `nextflow run hello-workflow.nf`, a następnie naciśnijmy enter.

Dobra, świetnie. Nasze trzy procesy działają. Mamy nasz katalog `results` z naszymi trzema wyjściami: Bonjour, Hello, Holà. Więc zamknijmy te pliki, zamknijmy terminal i wróćmy do skryptu.

## 1. Dodanie drugiego kroku do workflow'u

Dobra. W naszym przykładzie pozostajemy przy podstawach i staramy się być niezależni od dziedziny. Więc nasz drugi proces po prostu będzie manipulował tymi łańcuchami znaków, tymi słowami, w prosty sposób. Użyjemy polecenia Unix translate, aby wziąć te pliki i zamienić je na duże litery. Robimy to za pomocą polecenia `tr`.

## 1.1. Zdefiniowanie polecenia zamiany na wielkie litery i przetestowanie go w terminalu

Możemy to wypróbować w terminalu bash i sprawdzić, czy działa. Więc robisz `echo Hello World`, a następnie przekazujesz to z symbolem pipe do `tr`, i podajemy mu wzorzec rozpoznawania, od a do z, i na co powinien to przetłumaczyć: od A do Z wielkimi literami.

To jest bardzo proste, ponieważ dosłownie używa znaków od A do Z. Więc nie będzie działać na niczym ze znakami diakrytycznymi. Ale do celów przykładu zrozumiecie ideę.

Nacisnę enter i wypisuje w terminalu HELLO WORLD wielkimi literami. I tak jak poprzednio, moglibyśmy przekierować to do pliku, gdybyśmy chcieli. Outfile.

Dobra. Posprzątajmy to.

## 1.1. Napisanie kroku zamiany na wielkie litery jako proces Nextflow'a

Wróćmy do naszego skryptu i napiszmy nowy proces do obsługi tego polecenia bash. Skopiuję poprzedni proces, wkleję go poniżej i nazwę `CONVERT_TO_UPPER` — dla uppercase. Użyję tego samego `publishDir results`, ale wprowadzę kilka zmian tutaj. Zamiast przyjmować `val`, przyjmę `path input_file` i dodam tutaj prefiks `upper`, aby nasze pliki wyjściowe nie nadpisały poprzedniego wyjścia. I użyję nazwy zmiennej z wejścia. A potem zmienię skrypt tutaj na dole, i zamiast tego użyję `cat` na pliku wejściowym i tak jak zrobiliśmy w Bashu `tr a-z A-Z` oraz `upper_${input_file}.txt`. Dobra, kliknijmy zapisz.

## 1.2. Dodanie wywołania nowego procesu w bloku workflow

Teraz jeśli przewinę w dół, musimy faktycznie wywołać ten proces. Samo dodanie procesu do skryptu nie wystarczy. Musimy powiedzieć Nextflow'owi, że musimy uruchomić ten proces i gdzie to zrobić.

Więc tutaj dodam `CONVERT_TO_UPPER` i

dobra, dostajemy błąd mówiący, że oczekuje argumentu. Na pewno, musimy coś przekazać do tego procesu, aby faktycznie miał coś do zrobienia.

## 1.3. Przekazanie wyjścia pierwszego procesu do drugiego procesu

Weźmiemy wyjście z tego procesu. Więc biorę nazwę, `SAY_HELLO`, i gdy robię `.out`.

Dla tak prostego przykładu jak ten, gdzie mamy proces, który ma tylko jedno wyjście i przekazujemy je do nowego procesu z jednym wejściem, to powinno być wszystko, czego potrzebujemy. Więc kliknę zapisz, otworzę terminal i spróbujmy uruchomić to ponownie.

## 1.4. Ponowne uruchomienie workflow'u

Teraz nie uporządkowałem mojego katalogu `work` od ostatniego razu, kiedy uruchomiłem ten workflow. Uruchomię go ponownie i użyję tego jako okazji do pokazania, jak działa częściowe cache'owanie. Więc jeśli dodam `-resume` — miejmy nadzieję, że powinien ponownie użyć wyjść z pierwszego procesu, które były dokładnie takie same jak ostatnim razem. Ale teraz mamy nowy proces tutaj, który nie był wcześniej uruchamiany, który działa od podstaw. I na pewno, możesz zobaczyć, że pierwszy proces użył cache'owanych wyjść, a drugie wyjście uruchomiło trzy z trzech. Możesz również zobaczyć, że mamy teraz oba nasze procesy tutaj: nasz pierwszy proces, `SAY_HELLO`, uruchomiony trzy razy, i nasz drugi proces `CONVERT_TO_UPPER` uruchomiony trzy razy.

Jeśli uruchomię to ponownie, jako przypomnienie, z `-ansi-log false`, powinniśmy zobaczyć, że sześć różnych zadań procesu uruchomiło się, trzy dla każdego z nich. Więc to robi dokładnie to, czego się spodziewaliśmy. Pierwszy proces działa trzy razy, przekazując te wyjścia do drugiego procesu, który następnie działa trzy razy.

Przyjrzyjmy się wnętrzu katalogu `work` i zobaczmy, jak Nextflow obsługuje te pliki wejściowe. Jeśli wezmę ten katalog hash tutaj z drugiego procesu, możemy użyć polecenia `tree` ponownie z `-a` tylko po to, aby spojrzeć na te pliki. Widzisz tutaj, że mamy nasz plik wejściowy, którym jest plik `Bonjour-output.txt`, i to jest faktycznie symlink. To właśnie pokazuje nam ta strzałka i wskazuje na plik w poprzednim katalogu `work`.

To ma sens. Nextflow obsługuje wykonanie każdego zadania we własnym zamkniętym katalogu, więc jest całkowicie samodzielne. Jednakże musi dostarczyć pliki z poprzednich kroków jako wejście. Zamiast sięgać poza katalog `work`, aby pobrać te pliki, Nextflow umieszcza je w katalogu `work`.

Jeśli mamy współdzielony system plików jak tutaj, robi to za pomocą symlinku, tak aby nie używało dodatkowej przestrzeni na pliki. Jeśli używamy pamięci chmurowej z bucket'ami w różnych lokalizacjach, pobrałby te pliki i faktycznie skopiował je do katalogu `work`.

Przyjrzyjmy się plikowi `command.sh`. Jeśli zrobię `code work/command.sh`, możesz zobaczyć na pewno, że uzyskuje dostęp do tego pliku z lokalnego katalogu. Więc wszystko jest bardzo samodzielne i czyste.

Możemy również sprawdzić katalog `results` i upewnić się, że te pliki zostały poprawnie wyprowadzone. I na pewno, w `results` możemy zobaczyć wszystkie pliki wyjściowe z pierwszego procesu i wszystkie pliki wyjściowe z drugiego. I wszystkie są wielkimi literami, jak się spodziewaliśmy.

Tu zaczyna się pokazywać moc Nextflow'a. Za pomocą bardzo minimalnego kodu Nextflow obsługiwał równoległe wykonywanie tych zadań z czystą enkapsulacją w oddzielnych katalogach `work` oraz umieszczanie plików wejściowych i wyjściowych oraz publikowanie plików — wszystko automatycznie, od razu. Więc widzicie, jak, gdy skalujemy złożoność naszych workflow'ów analizy, ta funkcjonalność jest naprawdę, naprawdę cenna.

## 2. Dodanie trzeciego kroku do zebrania wszystkich powitań

Dobra. Te kroki były jeden-do-jednego. Mieliśmy jedno wyjście z pierwszego procesu trafiające do jednego wejścia dla drugiego procesu. Dalej będziemy mówić o tym, jak zebrać te różne wyjścia w jedno zadanie procesu, co znowu jest bardzo powszechną rzeczą do zrobienia. Więc szybko otwórzmy terminal i zróbmy próbny przebieg tego.

## 2.1. Zdefiniowanie polecenia zbierania i przetestowanie go w terminalu

Oszukam i skopiuję przykładowy kod bash z materiału szkoleniowego i po prostu nacisnę enter.

To, co możemy tutaj zobaczyć, to uruchomiliśmy to polecenie `echo` trzy razy do trzech różnych plików wyjściowych, które mogę tutaj zobaczyć. A następnie użyliśmy polecenia `cat`, aby wypisać wyjście każdego z tych trzech różnych plików i przekierować to do jednego zebranego pliku.

I jeśli zrobię `cat COLLECTED-output`, możesz zobaczyć, że zawiera zawartość tych trzech różnych plików, teraz w jednym pliku.

## 2.2. Utworzenie nowego procesu do wykonania kroku zbierania

Zobaczmy, czy możemy odtworzyć to samo w naszym pipeline'ie Nextflow.

Przewińmy w górę i stwórzmy trzeci proces. Skopiuję ten poprzedni, a tym razem nazwę go `COLLECT_GREETINGS`.

W terminalu bash nazwaliśmy to `collected-output.txt`. Więc powiem to samo: `path output` tutaj. I zrobię przekierowanie tutaj, więc jest zapisywane w ten sam sposób.

Dobra. Musimy zmienić to, co dzieje się na początku tego polecenia i musimy pomyśleć o tym, czym jest tutaj plik wejściowy. W rzeczywistości ten proces będzie przyjmował wiele plików wejściowych. Zachowam `path` i zmienię to na nową zmienną o nazwie `input_files`, w liczbie mnogiej.

Następnie ponownie użyję `cat` na nich, tak jak zrobiliśmy w naszym skrypcie bash. I użyję tutaj zmiennej.

Teraz możesz pomyśleć, że to nie zadziała. Widzieliśmy wcześniej błędy, gdzie tablica łańcuchów lub tablica ścieżek została przekazana do procesu i to spowodowało błąd. Ale w rzeczywistości tutaj Nextflow obsłuży to automatycznie dla nas we właściwy sposób. Weźmie kilka różnych plików wejściowych i po prostu wypisze różne ścieżki plików tutaj.

Oczywiście pomaga to, że polecenie `cat` może przyjmować serię nazw plików w ten sposób. Gdybym używał innego polecenia, które wymagało argumentu przed każdą ścieżką pliku lub czegoś, musielibyśmy mieć tutaj trochę więcej kodu i logiki, aby móc obsłużyć iterację tych ścieżek plików. Ale w tym przypadku powinno po prostu zadziałać.

## 2.3. Dodanie kroku zbierania do workflow'u

Dobra, zejdźmy do workflow'u i dodajmy nasz nowy proces: `COLLECT_GREETINGS`. I znowu weźmy wyjście z `CONVERT_TO_UPPER.out`. Zapiszmy to.

Wypróbujmy to. `nextflow run hello-workflow`.

Dobra, workflow się uruchomił, ale coś jest trochę dziwne. Mamy trzy wykonania pierwszego kroku, czego się spodziewamy. Trzy zadania dla drugiego, ale mamy również trzy zadania na końcu, kiedy spodziewaliśmy się mieć tylko jedno zadanie tutaj łączące wszystkie wyjścia.

Jeśli przejdziemy do naszego katalogu `results`, widzimy również, że `collected-output` ma tylko jedną wartość zamiast wszystkich trzech. To dlatego, że ten plik wyjściowy był nadpisywany trzy razy z trzema różnymi wartościami.

To ma sens, ponieważ przekazaliśmy jedno wyjście do jednego wejścia tutaj w ten sam sposób, jak zrobiliśmy to w poprzednim kroku.

## 2.4. Użycie operatora do zebrania powitań w jedno wejście

Więc potrzebujemy operatora tutaj, aby wziąć ten kanał z trzema elementami i zwinąć je do jednego elementu, tak aby ten końcowy proces uruchomił się tylko raz.

Aby to zrobić, użyjemy operatora `collect`. Mogę to zrobić bezpośrednio w workflow'ie. Mogę zrobić `.out` i połączyć łańcuchowo z operatorem tutaj na końcu `.collect`.

Kliknę zapisz. A następnie dla celów tego szkolenia dodam również kilka operatorów `view`, jak robiliśmy wcześniej, abyśmy mogli spojrzeć na ten kanał przed i po użyciu operatora `collect`, więc możemy zrozumieć, co się dzieje.

Wezmę ten kanał, pozbędę się `collect` i dodam `.view("before collect: ")`, a następnie zduplikuję tę linię, dodam operator `collect` i zmienię to na `"after collect: "`.

To jest oddzielne od miejsca, gdzie to wywołujemy, ale to w porządku, ponieważ używamy tych samych wywołań operatorów na tym samym kanale wyjściowym.

Dobra, kliknijmy zapisz i wypróbujmy to w terminalu. Uruchomię `nextflow run hello-workflow`. Ponownie uruchamiam nasz skrypt.

Dobra. To wygląda lepiej. Tak jak poprzednio możemy zobaczyć, że pierwsze dwa procesy uruchamiają się trzy razy, a teraz nasz końcowy proces uruchomił się tylko raz.

Jeśli spojrzymy na to, co zostało wypisane przez operator `view`, tutaj na dole, powiedzieliśmy `before collect`, co jest tym wyjściem tutaj, i to jest wypisane trzy razy. I widzisz, że jest pojedyncza ścieżka dla każdego z nich. A potem `after collect`, możesz zobaczyć, że mamy tę tablicę trzech ścieżek. Więc to jest zgodne z oczekiwaniami.

Dobra, sprawdźmy plik `results` i zobaczmy, czy tym razem jest zgodny z oczekiwaniami. Na pewno, są teraz trzy linie w pliku — to z powodzeniem połączyło te trzy wyjścia w jeden plik wyjściowy. Fantastycznie.

Dobra, posprzątam i przejdźmy do następnego kroku. I usunę te instrukcje `view` tylko po to, aby zachować porządek.

## 3. Przekazanie więcej niż jednego wejścia do procesu w celu unikalnego nazwania końcowego pliku wyjściowego

Dobra. Do tej pory wszystkie nasze procesy przyjmowały tylko jedno wejście. Teraz zamierzamy wykonać ćwiczenie, w którym dodamy więcej niż jedno wejście do procesu, aby zobaczyć, jak to działa. Aby to zrobić, użyjemy tego przykładu `COLLECT_GREETINGS`.

Za każdym razem, gdy uruchamiałem workflow, nadpisywał ten plik w katalogu `results`, co może nie być tym, czego chcemy.

## 3.1. Modyfikacja procesu zbierania, aby akceptował zdefiniowaną przez użytkownika nazwę dla pliku wyjściowego

Więc dla tego przykładu przekażemy dodatkowy parametr, abyśmy mogli dostosować nazwę pliku wyjściowego.

Dodanie drugiego wejścia do procesu jest bardzo proste. Po prostu dodaję drugą linię w bloku `input`. Tym razem to będzie `val`, a nie `path`, ponieważ chcemy przekazać łańcuch i nazwę to `batch_name`.

Mogę teraz użyć tej zmiennej w bloku `script` i powiem `collected-${batch_name}`.

Używam tutaj nawiasów klamrowych wokół nazwy zmiennej. To tylko po to, aby oddzielić ją od reszty łańcucha i prawdopodobnie nie jest to potrzebne w tym przypadku, ale myślę, że ułatwia to czytanie.

Dobra. Na koniec pamiętaj, aby zaktualizować ścieżkę wyjścia, ponieważ teraz nazwa pliku się zmieniła, więc zrobię to samo i umieszczę `batch_name` w wyjściu `path` zgodnie z oczekiwaniami.

## 3.2. Dodanie parametru wiersza poleceń batch

Teraz musimy przekazać nazwę batch skądś i stworzę drugi parametr, aby to zrobić, abyśmy mogli to zrobić w wierszu poleceń, gdy uruchamiamy workflow.

Więc zrobię `params.batch_name`, a domyślnie nazwijmy to `'test_batch'`. Teraz mogę użyć tej specjalnej zmiennej parametrów tam, gdzie wywołujemy proces.

I na pewno VS Code mówi nam, że nie ma wystarczającej liczby argumentów dla tego procesu teraz i że oczekuje drugiego wejścia.

Po prostu robię przecinek i przekazuję naszą nową zmienną i błąd znika.

Zauważ, że kolejność wejść tutaj jest naprawdę ważna. Pierwsze wejście procesu było `path`, a drugie wejście to `name`. Jeśli zmienię kolejność tutaj, muszę również zmienić kolejność, gdy wywołuję proces. W przeciwnym razie Nextflow przekaże niewłaściwy kanał do niewłaściwego wejścia.

## 3.3. Uruchomienie workflow'u

Dobra, wypróbujmy to i zobaczmy, czy działa. Zróbmy `nextflow run hello-workflow`. Dobra, uruchomił się jak poprzednio. Przyjrzyjmy się katalogowi `results`.

Na pewno, nasza nazwa pliku tutaj nazywa się teraz `collected-test_batch-output.txt`. Fantastycznie.

A teraz zobaczmy, czy możemy to nadpisać, uruchamiając ponownie. Tym razem zrobię `--batch_name`, aby dopasować tę nazwę specjalnej zmiennej parametru tutaj. I nazwę to `demo_output`.

Uruchomię workflow ponownie i zobaczymy, czy coś się stanie.

Dobra, teraz mamy `collected-demo_output-output.txt`. I ponieważ ta nazwa pliku różni się od tamtej, nie nadpisała jej. Obie są teraz obecne w katalogu `results`.

## 4. Dodanie wyjścia do kroku zbierania

Dobra, więc tam pokazaliśmy podawanie wielu wejść do procesu, ale co z wieloma wyjściami? W tym przykładzie obliczymy liczbę powitań, które są przetwarzane i wyprowadzamy to jako drugorzędne wyjście dla tego kroku `COLLECT_GREETINGS`.

## 4.1. Modyfikacja procesu, aby liczył i wyprowadzał liczbę powitań

Zrobimy tutaj małą sztuczkę. Procesy Nextflow'a mają ten blok `script` z wieloliniowym łańcuchem, i jest on przekazywany jako wyjście bash do `command.sh`. Ale możemy faktycznie napisać dowolny niestandardowy kod powyżej tego, i zostanie on wykonany jako część zadania, ale nie zostanie włączony w skrypt bash.

Jedną z wbudowanych funkcji w składni Nextflow'a jest funkcja `size`. Więc wezmę wejście `path` i powiem `count_greetings`, tylko po to, aby zdefiniować nazwę zmiennej. Wezmę `input_files` i wywołam `size()` na nich.

Ta funkcja policzy rozmiar tego kanału wejściowego i przypisze go do zmiennej.

Możemy teraz zwrócić tę zmienną jako część bloku `output`. Więc mówimy `val`, ponieważ jest to wartość, a nie plik. I `count_greetings`.

Teraz to samo w sobie jest wystarczające i moglibyśmy teraz uzyskać dostęp do tych różnych wyjść z tego procesu. Jednakże musielibyśmy uzyskać do nich dostęp w sposób pozycyjny. Więc używając klucza indeksu takiego jak zero i jeden.

Aby ułatwić sobie dostęp do wyjść, możemy je nazwać i robimy to za pomocą instrukcji `emit`.

Więc robię przecinek `emit outfile` lub cokolwiek chcę to nazwać. I robię tutaj `emit count`. To jest w zasadzie tylko dekorator, który pomaga nam pisać nieco czystszy kod, abyśmy mogli łatwo odwoływać się do konkretnych wyjść później w bloku `workflow`.

## 4.2. Raportowanie wyjścia na końcu workflow'u

Dobra. Jeśli przewinę w dół do bloku `workflow`, mogę teraz wziąć wyjścia `COLLECT_GREETINGS`, zrobić `COLLECT_GREETINGS.out`, i możemy zobaczyć nasze dwa nazwane wyjścia są sugerowane tutaj przez rozszerzenie VS Code. Bardzo przydatne.

Więc zrobię `.count`, aby uzyskać wartość count, którą właśnie stworzyliśmy, i zrobię `view`, aby zostało wypisane w wierszu poleceń. Więc możemy to zobaczyć, gdy uruchamiamy workflow.

Napiszmy coś w domknięciu tutaj tylko po to, aby było to trochę ładniejsze: `{ num_greetings -> "There were ${num_greetings} greetings" }`.

I tak naprawdę nie zależy nam na drugim wyjściu, ponieważ nie używamy go jako wejścia dla żadnych innych procesów. Ale widzicie, jak moglibyśmy łatwo przekazać to jako wejście do innego procesu, gdybyśmy chcieli, w dalszej części.

## 4.3. Uruchomienie workflow'u

Kliknę zapisz. Przyjrzyjmy się terminalowi i wypróbujmy to.

Dobra, fantastycznie. Proszę bardzo: `There were 3 greetings`. To dokładnie tak.

Dobra, świetna robota. To koniec tego rozdziału. Wszystko zrobione, że dotarliście tak daleko. Teraz zaczynasz budować dość realistyczny workflow, gdzie jesteśmy w stanie obsługiwać wejścia i wyjścia oraz logikę w naszym workflow'ie.

Gdy te pliki workflow stają się dłuższe, zaczynają być nieco nieporęczne. Więc w następnym rozdziale przyjrzymy się, jak możemy modularyzować kod Nextflow'a do oddzielnych plików, aby łatwiej było znaleźć i utrzymać kod w workflow'ie.

Dołącz do nas w następnym filmie dla rozdziału czwartego: Hello Modules.

[Następna transkrypcja wideo :octicons-arrow-right-24:](04_hello_modules.md)
