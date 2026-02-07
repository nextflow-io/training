# Część 3: Hello Workflow - Transkrypcja

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Tłumaczenie wspomagane przez AI - [dowiedz się więcej i zasugeruj ulepszenia](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

<div class="video-wrapper">
  <iframe width="560" height="315" src="https://www.youtube.com/embed/zJP7cUYPEbA?si=Irl9nAQniDyICp2b&amp;list=PLPZ8WHdZGxmXiHf8B26oB_fTfoKQdhlik" title="YouTube video player" frameborder="0" allow="accelerometer; autoplay; clipboard-write; encrypted-media; gyroscope; picture-in-picture; web-share" referrerpolicy="strict-origin-when-cross-origin" allowfullscreen></iframe>
</div>

!!!note "Ważne uwagi"

    Ta strona zawiera tylko transkrypcję. Aby uzyskać pełne instrukcje krok po kroku, wróć do [materiałów szkoleniowych](../03_hello_workflow.md).

    Numery sekcji pokazane w transkrypcji są podane tylko orientacyjnie i mogą nie zawierać wszystkich numerów sekcji z materiałów.

## Powitanie

Cześć, witamy w trzeciej części kursu "Hello Nextflow".

Ten rozdział nosi tytuł "Hello Workflow".

W rozdziale drugim zbudowaliśmy prosty workflow składający się z jednego procesu, ale w rzeczywistości pipeline'y są użyteczne, ponieważ mogą łączyć wiele kroków analizy razem.

W tym rozdziale weźmiemy ten początkowy przykład i rozszerzymy go, aby był nieco bardziej realistyczny.

Dodamy kilka dodatkowych kroków i zobaczymy, jak używamy kanałów do łączenia tych kroków.

Przyjrzymy się wielu zadaniom, które mogą się zwinąć do jednego procesu, oraz procesom, które mogą mieć wiele wejść i wiele wyjść.

Dobrze, zaczynajmy.

Więc zacznijmy. Tak jak poprzednio. Przejdźmy do training.nextflow.io. Hello Nextflow, rozdział trzeci. Hello Workflow. I otwórzmy nasz obszar roboczy. Uporządkowałem wszystkie moje pliki robocze z poprzednich rozdziałów i zamierzam otworzyć Hello Workflow.

To jest ten sam plik, nad którym pracowaliśmy do tej pory, więc powinien wyglądać znajomo. Mamy nasz proces say hello. Mamy nasz params.greeting z plikiem CSV greetings i mamy nasz workflow na dole, który ładuje ten plik CSV, tworzy kanał i przekazuje go do naszego procesu.

## 0. Rozgrzewka: Uruchomienie hello-workflow.nf

Jeśli chcesz, możemy to wypróbować i dokładnie sprawdzić, czy działa zgodnie z oczekiwaniami. Otwórzmy terminal, aby uruchomić nextflow run hello workflow nf i kliknij enter.

Dobrze, świetnie. Nasze trzy procesy działają. Mamy nasz katalog results z naszymi trzema wyjściami. Bonjour. Hello. Holà. Więc zamknijmy te pliki, zamknijmy terminal, wróćmy do skryptu.

## 1. Dodanie drugiego kroku do workflow'a

W naszym przykładzie pozostajemy przy podstawach i staramy się być niezależni od dziedziny. Więc nasz drugi proces po prostu będzie manipulował tymi łańcuchami znaków, tymi słowami, w prosty sposób. Użyjemy polecenia Unix translate, aby wziąć te pliki i zamienić je na duże litery. Robimy to za pomocą polecenia "tr".

## 1.1. Zdefiniowanie polecenia zamiany na wielkie litery i przetestowanie go w terminalu

Możemy to wypróbować po prostu w terminalu bash i sprawdzić, czy działa. Więc robisz echo, Hello World, a następnie przekazujesz to z symbolem pipe do tr i podajemy mu wzorzec rozpoznawania, a do z oraz na co powinien to przetłumaczyć. A do Z wielkimi literami.

To jest bardzo proste, ponieważ dosłownie używa znaków od A do Z. Więc nie zadziała na niczym, co jest akcentowane ani nic podobnego. Ale dla celów przykładu zrozumiecie ideę.

Zamierzam nacisnąć enter i wypisuje w terminalu HELLO WORLD wielkimi literami. I tak jak poprzednio, moglibyśmy przekierować to do pliku, gdybyśmy chcieli. Outfile.

Dobrze. Posprzątajmy to.

## 1.1. Napisanie kroku zamiany na wielkie litery jako proces Nextflow'a

Wróćmy do naszego skryptu i napiszmy nowy proces do obsługi tego polecenia bash. Zamierzam skopiować poprzedni proces, wkleić go poniżej i nazwać convert to upper. Dla uppercase. Zamierzam użyć tego samego publishDir results, ale zamierzam wprowadzić kilka zmian tutaj. Zamiast przyjmować val, zamierzam przyjąć path input file i zamierzam mieć tutaj prefiks upper, aby nasze pliki wyjściowe nie nadpisały wyjścia. I zamierzam użyć nazwy zmiennej z wejścia. A potem zamierzam zmienić skrypt tutaj na dole i zamiast tego zamierzam użyć cat na pliku wejściowym i tak jak zrobiliśmy w Bash TR, a-z, upper input file .txt. Dobrze, kliknijmy zapisz.

## 1.2. Dodanie wywołania nowego procesu w bloku workflow

Teraz jeśli przewinę w dół, musimy faktycznie wywołać ten proces. Samo dodanie procesu do skryptu nie wystarczy. Musimy powiedzieć Nextflow'owi, że musimy uruchomić ten proces i gdzie to zrobić.

Więc zamierzam tutaj, convert to upper i

dobrze, dostajemy błąd mówiący, że oczekuje argumentu. Na pewno, musimy coś przekazać do tego procesu, aby faktycznie miał coś do zrobienia.

## 1.3. Przekazanie wyjścia pierwszego procesu do drugiego procesu

To, co zamierzamy zrobić, to weźmiemy wyjście z tego procesu. Więc biorę nazwę, say hello i gdy robię dot out.

Dla tak prostego przykładu jak ten, gdzie mamy proces, który ma tylko jedno wyjście i przekazujemy je do nowego procesu, więc ma jedno wejście, to powinno być wszystko, czego potrzebujemy. Więc zamierzam kliknąć zapisz, otworzyć terminal i spróbujmy uruchomić to ponownie.

## 1.4. Ponowne uruchomienie workflow'a

Teraz nie uporządkowałem mojego katalogu work od ostatniego razu, kiedy uruchomiłem ten workflow. Zamierzam uruchomić go ponownie i zamierzam użyć tego jako okazji do pokazania, jak działa częściowe cache'owanie. Więc jeśli zrobię pojedynczą kreskę resume. Miejmy nadzieję, że powinien ponownie użyć wyjść z pierwszego procesu, które były dokładnie takie same jak ostatnim razem, kiedy uruchamiałem. Ale teraz mamy nowy proces tutaj, który nie był wcześniej uruchamiany, który działa od podstaw. I na pewno, widać, że pierwszy proces użył cache'owanych wyjść, a drugie wyjście uruchomiło trzy z trzech. Widzisz również, że mamy teraz oba nasze procesy tutaj: nasz pierwszy proces say hello uruchomiony trzy razy, a nasz drugi proces convert to upper uruchomiony trzy razy.

Jeśli uruchomię to ponownie, jako przypomnienie, z -ansi-log false, powinniśmy zobaczyć, że sześć różnych zadań procesu się uruchomiło, trzy dla każdego z nich. Więc to robi dokładnie to, czego się spodziewaliśmy. Pierwszy proces działa trzy razy, przekazując te wyjścia do drugiego procesu, który następnie działa trzy razy.

Więc przyjrzyjmy się wnętrzu katalogu work i zobaczmy, jak Nextflow obsługuje te pliki wejściowe. Jeśli wezmę ten katalog hash tutaj z drugiego procesu, możemy użyć polecenia tree ponownie z -a tylko po to, aby spojrzeć na te pliki. Widzisz tutaj, że mamy nasz plik wejściowy, którym jest plik Bonjour-output.txt i to jest faktycznie symlink. To właśnie pokazuje nam ta strzałka i wskazuje na plik w poprzednim katalogu work.

To ma sens. Nextflow obsługuje wykonanie każdego zadania we własnym zamkniętym katalogu, więc jest całkowicie samodzielne. Jednakże musi dostarczyć pliki z poprzednich kroków jako wejście. Zamiast sięgać poza katalog work, aby pobrać te pliki, Nextflow umieszcza je w katalogu work.

Jeśli mamy współdzielony system plików jak tutaj, robi to za pomocą symlinku, tak aby nie używało dodatkowej przestrzeni na pliki. Jeśli używamy pamięci w chmurze z bucket'ami w różnych lokalizacjach, pobrałby te pliki i faktycznie skopiował je do katalogu work.

Przyjrzyjmy się plikowi command sh. Jeśli zrobię code work, command sh, widać, na pewno, że uzyskuje dostęp do tego pliku z lokalnego katalogu. Więc wszystko jest bardzo samodzielne i czyste.

Możemy również sprawdzić katalog results i upewnić się, że te pliki zostały poprawnie wyprowadzone. I na pewno, w results widzimy wszystkie pliki wyjściowe z pierwszego procesu oraz wszystkie pliki wyjściowe z drugiego. I wszystkie są wielkimi literami, jak się spodziewaliśmy.

To tutaj zaczyna się pokazywać moc Nextflow'a. Za pomocą bardzo minimalnego kodu Nextflow obsługiwał równoległe wykonywanie tych zadań z czystą enkapsulacją w oddzielnych katalogach work oraz umieszczanie plików wejściowych i wyjściowych oraz publikowanie plików — wszystko automatycznie dla nas od razu. Więc widzicie, jak gdy skalujemy złożoność naszych workflow'ów analizy, ta funkcjonalność jest naprawdę bardzo cenna.

## 2. Dodanie trzeciego kroku do zebrania wszystkich powitań

Dobrze. Te kroki były jeden-do-jednego. Mieliśmy jedno wyjście z pierwszego procesu trafiające do jednego wejścia dla drugiego procesu. Dalej będziemy mówić o tym, jak zebrać te różne wyjścia w jedno zadanie procesu, co znowu jest bardzo powszechną rzeczą do zrobienia. Więc szybko otwórzmy terminal i zróbmy próbny przebieg tego.

## 2.1. Zdefiniowanie polecenia zbierania i przetestowanie go w terminalu

Zamierzam oszukać i skopiować przykładowy kod bash z materiału szkoleniowego i po prostu nacisnąć enter.

To, co tutaj widzimy, to uruchomiliśmy to polecenie echo trzy razy do trzech różnych plików wyjściowych, które widzę tutaj. A następnie użyliśmy polecenia cat, aby wypisać wyjście każdego z tych trzech różnych plików i przekierować to do jednego zebranego pliku.

I jeśli zrobię "cat COLLECTED-output", widzisz, że zawiera zawartość tych trzech różnych plików teraz w jednym pliku.

## 2.2. Utworzenie nowego procesu do wykonania kroku zbierania

Więc zobaczmy, czy możemy odtworzyć to samo w naszym pipeline'ie Nextflow'a.

Przewińmy w górę i stwórzmy trzeci proces. Zamierzam skopiować ten poprzedni, a tym razem zamierzam nazwać go Collect Greetings.

W terminalu bash nazwaliśmy to collected output txt. Więc zamierzam powiedzieć to samo path output tutaj. I zamierzam zrobić przekierowanie tutaj, więc jest zapisywane w ten sam sposób.

Dobrze. Musimy zmienić to, co dzieje się na początku tego polecenia i musimy pomyśleć o tym, czym jest tutaj plik wejściowy. W rzeczywistości ten proces będzie przyjmował wiele plików wejściowych. Zamierzam zachować path i zamierzam zmienić to na nową zmienną o nazwie input files, w liczbie mnogiej.

Następnie zamierzam ponownie użyć cat na nich, tak jak zrobiliśmy w naszym skrypcie bash. I zamierzam użyć tutaj zmiennej.

Teraz możesz pomyśleć, że to nie zadziała. Widzieliśmy wcześniej błędy, gdzie tablica łańcuchów lub tablica ścieżek została przekazana do procesu i to spowodowało błąd. Ale w rzeczywistości tutaj Nextflow obsłuży to automatycznie dla nas we właściwy sposób. Weźmie kilka różnych plików wejściowych i po prostu wypisze różne ścieżki plików tutaj.

Oczywiście pomaga to, że polecenie cat może przyjmować serię nazw plików w ten sposób. Gdybym używał innego polecenia, które wymagało argumentu przed każdą ścieżką pliku lub czegoś, musielibyśmy mieć tutaj trochę więcej kodu i logiki, aby móc obsłużyć iterację tych ścieżek plików. Ale w tym przypadku powinno po prostu zadziałać.

## 2.3. Dodanie kroku zbierania do workflow'a

Dobrze, zejdźmy do workflow'a i dodajmy nasz nowy proces. Collect greetings. I znowu weźmy wyjście z convert to upper out. Zapiszmy to.

Wypróbujmy to. nextflow run hello workflow.

Dobrze, workflow się uruchomił, ale coś jest trochę dziwne. Mamy trzy wykonania pierwszego kroku, czego się spodziewamy. Trzy zadania dla drugiego, ale mamy również trzy zadania na końcu, kiedy spodziewaliśmy się mieć tylko jedno zadanie tutaj łączące wszystkie wyjścia.

Jeśli przejdziemy do naszego katalogu results, widzimy również, że collected output ma tylko jedną wartość zamiast wszystkich trzech. To dlatego, że ten plik wyjściowy był nadpisywany trzy razy z trzema różnymi wartościami.

To ma sens, ponieważ przekazaliśmy jedno wyjście do jednego wejścia tutaj w ten sam sposób, jak zrobiliśmy to w poprzednim kroku.

## 2.4. Użycie operatora do zebrania powitań w jedno wejście

Więc potrzebujemy operatora tutaj, aby wziąć ten kanał z trzema elementami i zwinąć je do jednego elementu, tak aby ten końcowy proces uruchomił się tylko raz.

Aby to zrobić, użyjemy operatora collect. Mogę to zrobić bezpośrednio w workflow'ie. Mogę zrobić .out i połączyć łańcuchowo z operatorem tutaj na końcu .collect.

Kliknij zapisz. A następnie dla celów tego szkolenia zamierzam również zrobić kilka operatorów view, jak robiliśmy wcześniej, abyśmy mogli spojrzeć na ten kanał przed i po użyciu operatora collect, więc możemy zrozumieć, co się dzieje.

Zamierzam wziąć ten kanał, pozbyć się collect i dot view greetings, a następnie zamierzam zduplikować tę linię, dodać operator collect. I zmienić to na after.

To jest oddzielne od miejsca, gdzie to wywołujemy, ale jest w porządku, ponieważ używamy tych samych wywołań operatorów na tym samym kanale wyjściowym.

Dobrze, kliknijmy zapisz i wypróbujmy to w terminalu. Zamierzam uruchomić nextflow run. Hello, workflow. Ponownie uruchomić nasz skrypt.

Dobrze. To wygląda lepiej. Tak jak poprzednio widzimy, że pierwsze dwa procesy uruchamiają się trzy razy, a teraz nasz końcowy proces uruchomił się tylko raz.

Jeśli spojrzymy na to, co zostało wypisane przez operator view, tutaj na dole, powiedzieliśmy before collect, co jest tym wyjściem tutaj i jest wypisane trzy razy. I widzisz, że jest pojedyncza ścieżka dla każdego z nich. A potem after collect, widać, że mamy tę tablicę trzech ścieżek. Więc to jest zgodne z oczekiwaniami.

Dobrze, sprawdźmy plik results i zobaczmy, czy tym razem jest zgodny z oczekiwaniami. Na pewno, są teraz trzy linie w pliku — to z powodzeniem połączyło te trzy wyjścia w jeden plik wyjściowy. Fantastycznie.

Dobrze, zamierzam posprzątać i przejdźmy do następnego kroku. I zamierzam usunąć te instrukcje view tylko po to, aby zachować porządek.

## 3. Przekazanie więcej niż jednego wejścia do procesu w celu unikalnego nazwania końcowego pliku wyjściowego

Dobrze. Do tej pory wszystkie nasze procesy przyjmowały tylko jedno wejście. Teraz zamierzamy wykonać ćwiczenie, w którym dodamy więcej niż jedno wejście do procesu, aby zobaczyć, jak to działa. Aby to zrobić, użyjemy tego przykładu collect greetings.

Za każdym razem, gdy uruchamiałem workflow, nadpisywał ten plik w katalogu results, co może nie być tym, czego chcemy.

## 3.1. Modyfikacja procesu zbierającego, aby akceptował zdefiniowaną przez użytkownika nazwę dla pliku wyjściowego

Więc dla tego przykładu zamierzamy przekazać dodatkowy parametr, abyśmy mogli dostosować nazwę pliku wyjściowego.

Dodanie drugiego wejścia do procesu jest bardzo proste. Po prostu dodaję drugą linię w bloku input. Tym razem to będzie value, a nie path, ponieważ chcemy przekazać łańcuch i zamierzam nazwać to batch underscore name.

Mogę teraz użyć tej zmiennej w bloku script i zamierzam powiedzieć collected dash dollar batch name.

Używam tutaj nawiasów klamrowych wokół nazwy zmiennej. To tylko po to, aby oddzielić ją od reszty łańcucha i prawdopodobnie nie jest to potrzebne w tym przypadku, ale myślę, że ułatwia to czytanie.

Dobrze. Na koniec pamiętaj, aby zaktualizować ścieżkę wyjścia, ponieważ teraz nazwa pliku się zmieniła, więc zamierzam zrobić to samo i umieścić batch name w wyjściu path zgodnie z oczekiwaniami.

## 3.2. Dodanie parametru wiersza poleceń batch

Teraz musimy przekazać nazwę batch skądś i zamierzam stworzyć drugi parametr, aby to zrobić, abyśmy mogli to zrobić w wierszu poleceń, gdy uruchamiamy workflow.

Więc zamierzam zrobić params batch name, a domyślnie nazwijmy to test batch. Teraz mogę użyć tej specjalnej zmiennej parametrów tam, gdzie wywołujemy proces.

I na pewno VS Code mówi nam, że nie ma wystarczającej liczby argumentów dla tego procesu teraz i że oczekuje drugiego wejścia.

Po prostu robię przecinek i przekazuję naszą nową zmienną i błąd znika.

Zauważ, że kolejność wejść tutaj jest naprawdę ważna. Pierwsze wejście procesu było path, a drugie wejście to name. Jeśli zmienię kolejność tutaj, muszę również zmienić kolejność, gdy wywołuję proces. W przeciwnym razie Nextflow przekaże niewłaściwy kanał do niewłaściwego wejścia.

## 3.3. Uruchomienie workflow'a

Dobrze, wypróbujmy to i zobaczmy, czy działa. Zróbmy "nextflow run hello-workflow. Dobrze, uruchomił się jak poprzednio. Przyjrzyjmy się katalogowi results.

Na pewno, nasza nazwa pliku tutaj nazywa się teraz "collected test batch output txt". Fantastycznie.

A teraz zobaczmy, czy możemy to nadpisać, uruchamiając ponownie. Tym razem zamierzam zrobić --batch_name, aby dopasować tę nazwę specjalnej zmiennej parametru tutaj. I zamierzam nazwać to demo output.

Uruchom workflow ponownie i zobaczymy, czy coś się stanie.

Dobrze, teraz mamy collected demo output .txt. I ponieważ ta nazwa pliku różni się od tamtej, nie nadpisała jej. Obie są teraz obecne w katalogu results.

## 4. Dodanie wyjścia do kroku zbierającego

Dobrze, więc tam pokazaliśmy podawanie wielu wejść do procesu, ale co z wieloma wyjściami? W tym przykładzie zamierzamy obliczyć liczbę powitań, które są przetwarzane i wyprowadzić to jako drugorzędne wyjście dla tego kroku collect greeting.

## 4.1. Modyfikacja procesu, aby liczył i wyprowadzał liczbę powitań

Zamierzamy zrobić tutaj małą sztuczkę. Procesy Nextflow'a mają ten blok script z wieloliniowym łańcuchem i jest on przekazywany jako wyjście bash do dot command dot sh. Ale możemy faktycznie napisać dowolny niestandardowy kod powyżej tego i zostanie on wykonany jako część zadania, ale nie zostanie włączony w skrypt bash.

Jedną z wbudowanych funkcji w składni Nextflow'a jest funkcja size. Więc zamierzam wziąć wejście path i zamierzam powiedzieć count underscore greetings, tylko po to, aby zdefiniować nazwę zmiennej. Zamierzam wziąć pliki wejściowe i zamierzam wywołać "size" na nich.

Ta funkcja policzy rozmiar tego kanału wejściowego i przypisze go do zmiennej.

Możemy teraz zwrócić tę zmienną jako część bloku output. Więc mówimy val, ponieważ jest to wartość, a nie plik. I count greetings.

Teraz to samo w sobie jest wystarczające i moglibyśmy teraz uzyskać dostęp do tych różnych wyjść z tego procesu. Jednakże musielibyśmy uzyskać do nich dostęp w sposób pozycyjny. Więc używając klucza indeksu takiego jak zero i jeden.

Aby ułatwić sobie dostęp do wyjść, możemy je nazwać i robimy to za pomocą instrukcji emit.

Więc robimy przecinek emit out file lub cokolwiek chcę to nazwać. I robię tutaj emit count. To jest w zasadzie tylko dekorator, który pomaga nam pisać nieco czystszy kod, abyśmy mogli łatwo odwoływać się do konkretnych wyjść później w bloku workflow.

## 4.2. Raportowanie wyjścia na końcu workflow'a

Dobrze. Jeśli przewinę w dół do bloku workflow, mogę teraz wziąć wyjścia collect greetings, zrobić collect greetings, dot out i widzimy, że nasze dwa nazwane wyjścia są sugerowane tutaj przez rozszerzenie VS Code. Bardzo przydatne.

Więc zamierzam zrobić dot count, aby uzyskać wartość count, którą właśnie stworzyliśmy i zamierzam zrobić view, aby zostało wypisane w wierszu poleceń. Więc możemy to zobaczyć, gdy uruchamiamy workflow.

Napiszmy coś w domknięciu tutaj tylko po to, aby było to trochę ładniejsze. num greetings, there were greetings greetings.

I tak naprawdę nie zależy nam na drugim wyjściu, ponieważ nie używamy go jako wejścia dla żadnych innych procesów. Ale widzicie, jak moglibyśmy łatwo przekazać to jako wejście do innego procesu, gdybyśmy chcieli, w dalszej części.

## 4.3. Uruchomienie workflow'a

Zamierzamy kliknąć zapisz. Przyjrzyjmy się terminalowi i wypróbujmy to.

Dobrze, fantastycznie. Proszę bardzo. There are three greetings. To dokładnie tak.

Dobrze, świetna robota. To koniec tego rozdziału. Gratulacje, że dotarliście tak daleko. Teraz zaczynasz budować dość realistyczny workflow, gdzie jesteśmy w stanie obsługiwać wejścia i wyjścia oraz logikę w naszym workflow'ie.

Gdy te pliki workflow'ów stają się dłuższe, zaczynają być nieco nieporęczne. Więc w następnym rozdziale przyjrzymy się, jak możemy modularyzować kod Nextflow'a do oddzielnych plików, aby łatwiej było znaleźć i utrzymać kod w workflow'ie.

Dołącz do nas w następnym filmie dla rozdziału czwartego. Hello Modules.

[Następna transkrypcja wideo :octicons-arrow-right-24:](04_hello_modules.md)
