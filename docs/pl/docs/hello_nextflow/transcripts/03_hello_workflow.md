# Część 3: Hello Workflow - Transkrypcja wideo

<div class="video-wrapper">
  <iframe width="560" height="315" src="https://www.youtube.com/embed/_aO56V3iXGI?si=Irl9nAQniDyICp2b&amp;list=PLPZ8WHdZGxmWKozQuzr27jyMGqp9kElVK&amp;cc_load_policy=1" title="YouTube video player" frameborder="0" allow="accelerometer; autoplay; clipboard-write; encrypted-media; gyroscope; picture-in-picture; web-share" referrerpolicy="strict-origin-when-cross-origin" allowfullscreen></iframe>
</div>

!!!note "Ważne uwagi"

    Ta strona zawiera wyłącznie transkrypcję. Pełną instrukcję krok po kroku znajdziesz w [materiałach szkoleniowych](../03_hello_workflow.md).

    Numery sekcji w transkrypcji służą wyłącznie celom orientacyjnym i mogą nie pokrywać się ze wszystkimi numerami sekcji w materiałach.

## Powitanie i podsumowanie

Cześć i witaj z powrotem w części trzeciej kursu Hello Nextflow. Ta część nosi nazwę Hello Workflow i to właśnie tutaj zaczynamy uzasadniać nazwę "pipeline" czy "workflow".

Weźmiemy nasz prosty skrypt pipeline'u, który do tej pory zawierał tylko jeden proces, i zaczniemy dodawać kolejne procesy, aby zobaczyć, jak Nextflow radzi sobie z orkiestracją i przepływem danych przez pipeline'a.

Wróćmy do naszego środowiska Codespaces. Widzicie, że usunąłem wszystkie katalogi .nextflow\* i katalogi work, żeby utrzymać porządek. Nie martwcie się, jeśli nadal macie te pliki z poprzednich części kursu.

Będziemy pracować na pliku o nazwie hello-workflow.nf. Jak wcześniej, reprezentuje on skrypt, który zbudowaliśmy do tej pory, i daje nam czysty punkt wyjścia. Ponownie, w sekcji output widzimy, że ścieżka to teraz hello_workflow. Więc publikowane pliki powinny trafiać do innego podkatalogu w katalogu results.

Podsumujmy, gdzie jesteśmy. Mamy pojedynczy proces z jednym wejściem (greeting) i jednym wyjściem (greeting file). A następnie prosty skrypt Bash, który wykonuje polecenie echo do pliku.

Mamy pojedyncze wejście workflow'a - blok params - gdzie mówimy, że oczekuje ścieżki, a wartość domyślna to data/greetings.csv, czyli ten plik tutaj.

W samym workflow'ie mamy blok main. Tworzymy kanał, parsujemy CSV na wiersze, a następnie bierzemy pierwszy element każdej tablicy i przekazujemy ten kanał do procesu, który generuje trzy zadania. Publikujemy z workflow'a wyjścia z tego procesu.

Na końcu w bloku output mówimy Nextflow'owi, żeby opublikował te pliki z tego kanału do katalogu o nazwie hello_workflow. Oraz żeby skopiował te pliki, a nie tworzył dowiązań symbolicznych.

## 1. Dodanie drugiego kroku do workflow'a

Dobrze, w tej części dodamy drugi proces do naszego workflow'a. Weźmiemy wyjścia z procesu sayHello i przetworzymy je w drugim kroku, który zamieni wszystkie litery w tych plikach na wielkie - convertToUppercase.

To tylko przykładowy przypadek, kolejne proste przetwarzanie łańcuchów znaków, ale pokazuje, jak możemy zastosować logikę w workflow'ie.

Użyjemy do tego polecenia bash o nazwie "tr", skrót od "translate". To polecenie uniksowe istniejące od zawsze. Jeśli nie znacie go, to nic dziwnego. Sam chyba nigdy nie używałem go przed szkoleniem, ale możecie je szybko wypróbować w terminalu. Jeśli wpiszę "echo 'hello world'", a następnie przekieruję przez pipe do 'tr', a potem w cudzysłowach podam zakres znaków, więc od A do Z, małe litery, a następnie chcę od A do Z wielkie litery. Po prostu mówi: przetłumacz te litery na te litery.

Gdy nacisnę enter, zobaczysz, że wszystko zostało zamienione na wielkie litery. Bardzo fajne, jeśli lubicie krzyczeć na ludzi.

To bardzo prosty styl polecenia bash, którego użyjemy w naszym drugim procesie.

## 1.2. Napisanie kroku zamiany na wielkie litery jako procesu Nextflow

Więc jeśli wrócę do mojego skryptu, trochę oszukam i po prostu skopiuję kod z dokumentacji szkoleniowej. Ale widzicie dokładnie, co się dzieje.

Mamy nowy proces tutaj. Ten nazwaliśmy convertToUpper, ale moglibyśmy nazwać go jak tylko chcemy.

Mamy jedno wejście typu path, tak jak wcześniej. To nie jest kanał wartości, to kanał ścieżki. A następnie pojedyncze wyjście.

W bloku script robimy "cat" na pliku wejściowym. Możemy umieścić to w klamrach, jeśli chcemy, co pobiera tę zmienną. Uruchamiamy to samo polecenie bash w pipe i zapisujemy wyniki do pliku z tą nazwą, która jest przechwytywana przez wyjściową ścieżkę.

Teraz musimy coś zrobić z tym nowym procesem. Więc przejdziemy do workflow'a, gdzie budujemy różną logikę workflow'a, i po pierwszym procesie uruchomimy nasz drugi proces. Więc convertToUpper to nazwa procesu tutaj.

Przyjmuje wejście, więc nie możemy go po prostu wywołać samodzielnie. Chcemy przetworzyć wyjście pierwszego procesu. Więc tak jak zrobiliśmy to z sayHello out, gdzie publikujemy te wyniki, chcemy użyć tych samych wyników tutaj jako wejścia, więc możemy to skopiować i tam wstawić.

Chcemy proces sayHello z ".out", a Nextflow wie, że oznacza to prosty pojedynczy rekord wyjściowy tutaj, którym jest ten plik. Więc to zostanie przekazane jako wejście do drugiego procesu.

## 1.5. Konfiguracja publikowania wyjścia workflow'a

Dobra. I w końcu, żebyśmy faktycznie zapisali wyniki tego drugiego procesu, musimy je również opublikować z workflow'a, a następnie zdefiniować je w bloku output, ta sama składnia co wcześniej. Więc możemy to skopiować i powiedzieć second outputs, lub nazwać to jak chcecie.

Weźmiemy nazwę procesu, który nas interesuje, convertToUpper out, a następnie tutaj w dole w bloku output. Dodamy to i moglibyśmy użyć tych samych atrybutów tutaj. Więc chcemy również te pliki w podkatalogu Hello Workflow i chcemy je również skopiować.

Świetnie. Spróbujmy to uruchomić. Więc jeśli wywołam terminal i wpiszę "nextflow run hello-workflow.nf", zobaczymy, co to da. Sprawdzimy, czy wygląda to inaczej niż w poprzednich częściach.

Więc uruchamia Nextflow. W dokumentacji jest napisane, żeby zrobić to z "-resume", ale usunąłem cały mój katalog work, więc i tak by to nie zrobiło różnicy tutaj. Ale jeśli Wy tego nie zrobiliście, to też zadziała.

I wygląda to prawie dokładnie tak samo. Ale widzicie teraz drugą linię wyjścia tutaj, gdzie możecie zobaczyć nazwę drugiego procesu, który właśnie dodaliśmy. I rzeczywiście, widzicie, że uruchomił się trzy razy pomyślnie.

Świetnie. Gdybym miał poprzednie katalogi work i zrobiłbym to z "-resume", te byłyby zacache'owane tylko dla pierwszego kroku w pipeline'ie. Bo te wyjścia były dokładnie takie same, więc Nextflow wiedziałby, żeby je ponownie użyć.

I więc widzicie, jak możecie używać -resume do iteracyjnego budowania workflow'a, krok po kroku, jeśli trzeba.

Dobra, spójrzmy do katalogu results tutaj na górze i zobaczmy, czy zadziałało. Widzimy, że mamy tutaj więcej plików. Mamy nasze oryginalne pliki jak wcześniej z pierwszego procesu. I rzeczywiście, mamy nasze pliki upper i litery są wszystkie wielkie, więc zadziałało. Naprawdę miło to widzieć.

Ciekawe jest również zajrzenie do tych katalogów work. Jak wcześniej, hash tutaj odpowiada katalogom work. Więc jeśli spojrzę do "ls work", a następnie rozwinę to, zobaczymy różne pliki tutaj.

Widzimy plik wyjściowy z pierwszego procesu, który został tutaj wciągnięty jako wejście. I widzimy nowy plik wyjściowy, który został wygenerowany.

Teraz jeśli zrobię to z "-la", żeby wylistować i pokazać wszystkie pliki, zobaczymy kilka więcej rzeczy. Po pierwsze, zobaczycie, że ten plik jest w rzeczywistości dowiązaniem symbolicznym do pierwszego procesu. To prawie zawsze jest dowiązanie symboliczne, jeśli to możliwe, żeby zaoszczędzić miejsce na dysku. Nie publikujemy tutaj plików i po prostu odwołuje się do tego pliku z pierwszego zadania do drugiego zadania, tak że wszystko jest zamknięte w tym jednym katalogu roboczym, bezpieczne i odizolowane od wszystkiego innego.

I to musi tam być, bo jeśli spojrzymy na plik .command.sh, więc jeśli zrobię "cat work/b8/56\*", zobaczysz, że części plików tutaj są względne, więc robi cat tego pliku wejściowego, który został dowiązany symbolicznie do tego samego katalogu roboczego.

Tak więc tak będzie wyglądał każdy katalog work. Gdy spojrzysz na to w Nextflow, będziesz miał wszystkie pliki wejściowe tam przygotowane do tego katalogu roboczego. A następnie będziesz miał również wszystkie pliki wyjściowe, które zostały utworzone. Więc to świetnie. To wygląda, jak się spodziewaliśmy.

## 2.1. Zdefiniowanie polecenia zbierającego i przetestowanie go w terminalu

Dobra, wróćmy do naszego workflow'a. Jaki jest następny krok, który chcemy wykonać?

Mamy teraz dwa procesy i pobierają one ten jeden plik CSV, parsują go i dzielą. A następnie mamy trzy zadania dla każdego z tych procesów i Nextflow zajmuje się paralelizacją tego wszystkiego, więc wszystko działa równolegle, gdzie to możliwe.

Ten sposób dzielenia pracy, żeby uruchamiać rzeczy równolegle, jest bardzo powszechny. I odwrotnością tego jest zbieranie wszystkiego z powrotem. Więc to zrobimy z naszym ostatnim procesem w workflow'ie - będziemy mieć tutaj trzeci, który pobiera te trzy różne wyjścia i łączy je wszystkie w jeden plik.

Możemy to zrobić całkiem prosto w terminalu, żeby poczuć, jak to będzie wyglądało.

Jeśli przejdę do katalogu results. Więc, "cd results/hello_workflow/", mamy tutaj wszystkie pliki UPPER. Mogę po prostu użyć "cat", którego używamy do wypisania zawartości tego pliku, i możesz podać wiele plików do "cat" i będzie czytał jeden po drugim.

Więc mogę powiedzieć "UPPER-\*", co daje mi tę samą listę trzech nazw plików z rozwinięciem Bash. I mogę powiedzieć combined.txt. Myślę, że w dokumentacji wymienione są dokładne nazwy plików, ale to robi to samo.

Teraz, jeśli użyję "cat combined.txt", zobaczymy, że mamy zawartość wszystkich trzech tych plików.

Więc to w zasadzie wszystko, co ten proces będzie robił - spróbujemy dać mu wszystkie różne pliki wyjściowe z poprzedniego procesu w jednym zadaniu procesu, a następnie "cat" je razem i zapiszemy plik wyjściowy.

## 2.2. Utworzenie nowego procesu do wykonania kroku zbierającego

Dobra, więc dodajmy nasz nowy proces. Wkleję to z materiałów szkoleniowych i widzicie, że zostawiono nam trochę ćwiczenia dla czytelnika tutaj z tymi znakami zapytania. Ale widzicie ogólny zarys procesu to w zasadzie to, co właśnie zrobiliśmy w terminalu, gdzie robimy "cat" wielu plików wejściowych i zapisujemy to do pliku wyjściowego tutaj o nazwie collected, a następnie wyjście oczekuje tej pojedynczej ścieżki ponownie.

Więc potrzebujemy jakiegoś wejścia tutaj i będzie to zbiór ścieżek. Więc ponownie, definiujemy wejściowy kanał ścieżki i nazwijmy go input_files. Teraz, to wcześniej dawało nam pojedynczą ścieżkę tutaj, ale ścieżka może również mieć wiele plików tutaj, nawet jeśli to nadal pojedyncza deklaracja.

Skopiuję to tutaj w dół, bo chcemy zrobić "cat" tych plików. I możecie pomyśleć, że mamy jakieś problemy tutaj z drukowaniem tablicy czy coś w tym stylu, ale Nextflow jest ogólnie dość rozsądny, jeśli chodzi o to. I jeśli dostanie kanał z wieloma plikami w nim w ten sposób, umieści je wszystkie razem z separatorami spacji. Więc to da nam poprawną składnię.

To świetnie. Więc teraz podłączmy nasz nowy proces. Przechodzę do workflow'a. Powiem combine the outputs, nowa nazwa procesu, i tak samo jak wcześniej. Wezmę ten poprzedni proces, convertToUpper i zrobię ".out".

Świetnie. Wypróbujmy to i zobaczmy, czy działa w terminalu. Jeśli po prostu wrócę o kilka katalogów w górę, a następnie uruchomię ponownie polecenie Nextflow, zobaczymy, co się stanie.

Więc workflow został uruchomiony i teraz widzicie, że mamy trzy różne nazwy procesów, co jest świetne. Pierwsze dwa wyglądają tak samo jak wcześniej, a trzeci nowy się uruchamia, co jest dobre.

Jednak jest coś dziwnego. Chcieliśmy połączyć te pliki wyjściowe w jeden plik, a jednak ten proces, widzimy, uruchomił się trzy razy, a nie raz.

Rzeczywiście, jeśli przejdziemy do jednego z tych katalogów work. I zrobimy "cat work/" "collected", wtedy zobaczymy. Jest tutaj tylko jedno słowo, nie trzy.

I więc co się stało, to Nextflow kontynuował tę paralelizację tak jak w poprzednich krokach. I ten proces dał nam kanał z trzema elementami, i te trzy elementy kanału zostały przekazane do naszego procesu downstream, który wygenerował trzy zadania procesu.

W zasadzie próbował zebrać trzy oddzielne razy i za każdym razem miał tylko pojedynczy plik, więc po prostu zrobił cat pojedynczego pliku do wyjścia, i w rzeczywistości możemy to również zobaczyć w pliku .command.sh.

Jeśli zrobię .command.sh, zobaczymy, że ma tylko pojedynczą nazwę pliku tutaj i tylko pojedynczy plik został przygotowany do tego katalogu roboczego.

## 2.3. Dodanie kroku zbierającego do workflow'a

Więc jakoś musimy powiedzieć Nextflow, żeby zebrał wszystkie te wyjścia z poprzedniego procesu i dał je temu procesowi downstream jako pojedynczy element kanału, a nie trzy.

Robimy to za pomocą operatora kanału o nazwie _collect_.

To super przydatny operator, który zobaczycie w pipeline'ach Nextflow cały czas. To jest kanał tutaj, ten kanał wyjściowy, tak samo jak ten, który stworzyliśmy na górze. I więc możemy dołączać do niego operatory kanału, tak jak robiliśmy wcześniej. Możemy po prostu zrobić kropkę, a następnie w tym przypadku collect, nawiasy.

I to wszystko, czego potrzebujemy. To wtedy zmanipuluje ten kanał, zanim zostanie przekazany do tego procesu.

Jeśli chcecie zobaczyć, co się z nim dzieje, możemy również go wyświetlić tutaj. Więc tutaj, to nie jest związane z uruchamianiem tego procesu w ogóle, więc mogłem to umieścić w dowolnym punkcie po uruchomieniu tego procesu. Ale bierzemy ten sam kanał wyjściowy i patrzymy na niego z .view, a następnie patrzymy na niego ponownie z .collect.view.

I gdy to uruchomimy, pokaże nam dwie różne struktury tego kanału, przed i po collect. Więc spróbujmy tego teraz. Dobra, właśnie trochę pomniejszyłem, bo niektóre wyjścia są dość długie, ale jeśli uruchomię pipeline'a, zobaczymy, czy działa.

Mam nadzieję, że trzeci proces uruchomi się tylko raz, bo zbiera wyjścia i rzeczywiście, widzicie collectGreetings jako jeden z jednego. Więc uruchomił się tylko jedno zadanie.

A następnie jeśli spojrzymy na instrukcje view, mamy trzy instrukcje view dla trzech elementów przed, z jedną ścieżką pliku w każdym.

A następnie po tej instrukcji collect, została wyzwolona tylko raz, bo jest pojedynczy element w tym kanale. I teraz mamy tę listę trzech różnych ścieżek plików.

To dokładnie to, na co liczyliśmy. I widzicie, mam nadzieję, że to w zasadzie odwrotność tego operatora `map`, którego użyliśmy, żeby przejść z tablic CSV do oddzielnych elementów kanału. Teraz bierzemy oddzielne elementy kanału i wkładamy je z powrotem do pojedynczej tablicy.

Świetnie, możemy wyczyścić te instrukcje view. Nie potrzebujemy ich już. Możemy przejść do następnego kroku.

Zanim pójdę dalej i zanim zapomnę, dodam nową instrukcję publish tutaj. Third output. Możecie nazwać to czymś bardziej semantycznym i opisowym w Waszym workflow'ie. A następnie dodam to do bloku output ponownie i powiem path 'hello_workflow' mode 'copy'. Po prostu żeby plik wyjściowy wygenerowany przez ten proces został zapisany do naszego katalogu results tutaj na górze.

Tylko żeby szybko sprawdzić, czy działa. Powinno być teraz trochę czystsze, bo nie mamy tych instrukcji view. I zobaczymy, czy otrzymamy nasz nowy plik wyjściowy tutaj na górze. Jedno z jednego zadania uruchomiło się, dostaliśmy nowy plik o nazwie collected, i teraz mamy wszystkie trzy te słowa. Fantastycznie. Co dalej?

## 3. Przekazanie dodatkowych parametrów do procesu

Dobra. Następnie przyjrzymy się obsłudze wielu wejść do pojedynczego procesu. Do tej pory widzicie, że wszystkie nasze procesy pobierają tylko jedną rzecz jako wejście. Wszystkie mają pojedynczą linię pod swoim wejściem.

Zademonstrujemy to pozwalając Nextflow określić inny identyfikator batch, tak że może uruchomicie ten workflow wielokrotnie i możecie nadać mu inny batch ID za każdym razem.

Po prostu dodam drugą linię we wejściu tutaj dla collectGreetings. I nazwę to "val", bo to jest łańcuch znaków. Teraz to jest wartość, nie ścieżka, i nazwę to "batch_name".

Następnie edytuję skrypt tutaj na dole, żeby użyć tej zmiennej, i spróbuję umieścić ją w tym samym miejscu co materiał szkoleniowy. Więc umieszczę ją w środku tej ścieżki pliku COLLECTED-$\{batch_name\}-output.

Jeszcze nie skończyliśmy. Pamiętajcie, że musimy powiedzieć Nextflow, jakie będą nazwy plików wyjściowych. Więc musimy też zrobić to samo tutaj na górze: COLLECTED-$\{batch_name\}-output.txt".

Fantastycznie. Nextflow teraz dostaje drugie wejście zmiennej i interpoluje to do skryptu i wyjścia.

Jeszcze jedna rzecz, teraz musimy znaleźć, gdzie to jest wywoływane, i musimy przekazać drugie wejście do procesu. To jest jak każde inne wejście do funkcji w każdym innym języku.

Tak jak zrobiliśmy wcześniej w szkoleniu, użyję specjalnego "params" tutaj, i nazwiemy to "params.batch", żebyśmy mogli mieć opcję CLI --batch. I teraz widzicie, że nasz proces tutaj ma dwa oddzielne wejścia oddzielone przecinkiem, które są przekazywane.

Naprawdę ważne jest, żeby uzyskać właściwą kolejność, więc kolejność argumentów tutaj dla kanału i następnie parametru musi się zgadzać. Kanał i nazwa batch tutaj. To jest po prostu dopasowanie pozycyjne.

Dobra. Mogę teraz uruchomić ten pipeline od razu z --batch, ale zróbmy najpierw właściwą rzecz i zdefiniujmy to we wejściu tutaj w Params. Więc dodam to do batch i następnie powiemy, że to jest łańcuch znaków i dajmy wartość domyślną. Więc po prostu nazwijmy to batch. Dobra? Teraz spróbujmy uruchomić workflow.

--batch Trio. Myślę, że tak jest napisane w materiale szkoleniowym, ale moglibyśmy użyć dowolnego łańcucha znaków tam. I mam nadzieję, że zobaczymy, że plik wyjściowy results pojawi się tutaj.

I rzeczywiście, COLLECTED-trio-output - zadziałało poprawnie. Zmienił nazwę naszego pliku. I możecie sobie wyobrazić, że teraz to jest przydatne, bo jeśli uruchomię to ponownie z inną nazwą batch, jak replicate_two, to da nam inną nazwę batch tutaj na górze.

I w tym przypadku nie nadpisze plików wyjściowych. Więc to miłe.

## 4. Dodanie wyjścia do kroku zbierającego

Dobra, więc mamy teraz wiele wejść do naszego procesu tutaj. Ale co się stanie, jeśli chcemy utworzyć wiele wyjść? Naszym przykładem tutaj będzie to, że utworzymy raport dla tego procesu, mówiący po prostu, ile plików zostało zebranych.

I zrobimy to za pomocą polecenia echo tutaj. Więc możemy powiedzieć echo. There were, skopiuję to z materiału szkoleniowego, żebyście nie musieli patrzeć, jak to wpisuję.

There were $\{count_greetings\} greetings in this batch, i zapiszemy to do nowego pliku teraz o nazwie $\{batch_name\}, więc ta sama zmienna, możemy jej używać tyle razy, ile chcemy, report.txt.

## 4.1.1. Policzenie liczby zebranych powitań

Musimy jakoś to faktycznie obliczyć. Moglibyśmy zrobić tę logikę w skrypcie Bash, gdybyśmy chcieli, używając logiki Bash. Jednak możemy również po prostu robić skryptowanie bezpośrednio w kodzie Nextflow, o ile znajduje się w bloku script w procesie i nad zacytowaną sekcją.

Cokolwiek tutaj nie zostanie włączone do ostatecznego wyrenderowanego skryptu i zostanie po prostu wykonane przez Nextflow, gdy renderuje zadanie.

Więc tutaj po prostu robimy trochę logiki. Tworzymy nową zmienną o nazwie count_greetings. Bierzemy tutaj kanał plików wejściowych i wywołujemy na nim .size().

Dobra, ta funkcja da mi tutaj liczbę do tej zmiennej, i teraz nasze ostrzeżenie zniknęło, bo ta zmienna jest definiowana.

Dobra, więc tworzymy ten drugi plik w katalogu work, ale musimy powiedzieć Nextflow, żeby oczekiwał go jako opublikowane wyjście tego procesu. Więc robimy to dokładnie tą samą składnią, jakiej użyliśmy dla pierwszego pliku.

Mówimy path, bo znowu, moglibyśmy publikować zmienną tutaj, gdybyśmy chcieli, z "val", ale powiemy "path". A następnie oczekiwaną nazwę pliku. Zauważcie, że nie jest tutaj podświetlona. To dlatego, że użyłem pojedynczych cudzysłowów. Muszę użyć podwójnych cudzysłowów.

## 4.1.2. Wyemitowanie pliku raportu i nazwanie wyjść

Dobra, to świetnie. I moglibyśmy teraz zacząć uzyskiwać dostęp do tych wyjść tutaj na dole tak jak zrobiłem tutaj. Ale teraz to jest tablica różnych obiektów, więc mógłbym zrobić collectGreetings.out[0], żeby dostać pierwszy, lub jeden, żeby dostać drugi, czyli nasz nowy raport.

Ale naprawdę nie lubię tego robić, bo dość łatwo pomylić się w liczeniu indeksów. I siedzicie tam licząc linie często i dodajecie nowe wyjście i nagle wszystko się psuje. Więc

znacznie lepiej jest odwoływać się do wszystkiego po nazwie zamiast tego. I możemy to zrobić za pomocą specjalnego klucza tutaj o nazwie "emit".

Więc możemy to nazwać, jak chcemy. Nazwijmy to emit outfile, i emit reports. Jeśli je zdefiniujecie i możecie to zrobić dla jednego lub wielu, zależy od Was. Teraz mogę przejść tutaj na dół i zamiast tego mogę przejść do dot out dot reports i po prostu wywołać to po nazwie, co jest znacznie łatwiejsze do zrozumienia Waszego kodu, gdy go czytacie, i jest bezpieczniejsze dla zmian w kodzie.

Dodałem tutaj .out.report, ale faktycznie muszę mieć dwa różne wyjścia publikowane. Więc zmienię nazwę na coś bardziej interesującego jak collected i report i czy tak to nazwałem? Nazwałem to out file, przepraszam. Więc ta nazwa emit tutaj outfile i report, bo publikujemy dwa różne kanały wyjściowe i więc musimy odwołać się do nich obu w bloku publish.

Następnie musimy również zdefiniować je w bloku output. Więc zmieniłem nazwę tego na collected, i znowu, dla reports, trochę rozwlekłe tutaj, ale jest to naprawdę przydatne, gdy przychodzisz czytać nowy workflow, żeby zobaczyć wszystkie różne wyjścia tutaj, wszystkie różne kanały wymienione obok siebie, i są sposoby, żeby to uczynić mniej rozwlekłym, których dotknę później.

Dobra, spróbujmy i uruchommy nasz workflow i zobaczmy, co się stanie.

Mam nadzieję, że teraz powinien działać w zasadzie tak samo jak wcześniej. I otrzymamy nowy plik wyjściowy tutaj na górze o nazwie replicate_two, report. I proszę. Został otwarty i mówi, że są trzy powitania w batch, co jest tym, czego oczekiwaliśmy, więc jest idealnie.

Jeśli przejdę do katalogu work tutaj tylko po to, żeby Wam udowodnić, że zostało wykonane w kodzie Nextflow, a nie w skrypcie bash, mogę przejść do cat work/ command.sh, i zobaczycie tutaj, że po prostu wypisuje ten łańcuch znaków bezpośrednio. There were three greetings in this batch, i więc ta zmienna została zinterpolowana przez Nextflow. Została obliczona w bloku script, zanim napisał plik .command.sh. Więc wynikowe obliczenie zmiennej jest w zasadzie zakodowane na stałe do tego, zanim zostanie wykonane w Waszym środowisku obliczeniowym w tym przypadku.

I więc widzicie to oddzielenie między blokiem script tutaj a wszystkim powyżej niego. Mam nadzieję, że to ma sens.

## Podsumowanie i quiz

Dobra, to koniec tej części kursu Hello Nextflow. Więc jak wcześniej, idźcie i sprawdźcie quiz. Zróbcie to na stronie internetowej lub w CLI, przejdźcie przez kilka pytań i sprawdźcie, czy zrozumieliście część materiału, który omówiliśmy. Zobaczcie, czy jest tam coś, co podkreśla coś, czego mogliście nie zrozumieć. Niewiele pytań. Łatwe i miłe do zrobienia. Lub możecie zrobić to na stronie internetowej tutaj na dole również.

I zróbcie sobie małą przerwę, małą przechadzkę i wracajcie i dołączcie do nas w części czwartej kursu Hello Nextflow, gdzie będziemy mówić o modułach. Dziękuję bardzo.
