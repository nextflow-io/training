# Część 3: Hello Workflow - Transkrypcja wideo

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Tłumaczenie wspomagane przez AI - [dowiedz się więcej i zaproponuj poprawki](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

<div class="video-wrapper">
  <iframe width="560" height="315" src="https://www.youtube.com/embed/_aO56V3iXGI?si=Irl9nAQniDyICp2b&amp;list=PLPZ8WHdZGxmWKozQuzr27jyMGqp9kElVK&amp;cc_load_policy=1" title="YouTube video player" frameborder="0" allow="accelerometer; autoplay; clipboard-write; encrypted-media; gyroscope; picture-in-picture; web-share" referrerpolicy="strict-origin-when-cross-origin" allowfullscreen></iframe>
</div>

!!!note "Ważne uwagi"

    Ta strona zawiera tylko transkrypcję. Pełne instrukcje krok po kroku znajdziesz w [materiałach szkoleniowych](../03_hello_workflow.md).

    Numery sekcji pokazane w transkrypcji mają charakter orientacyjny i mogą nie obejmować wszystkich numerów sekcji w materiałach.

## Powitanie i podsumowanie

Cześć i witaj ponownie w trzeciej części Hello Nextflow. Ta część nazywa się Hello Workflow i to właśnie w tej części szkolenia zaczynamy naprawdę uzasadniać nazwę pipeline czy workflow.

Weźmiemy nasz prosty skrypt pipeline'u z jednym procesem i zaczniemy dodawać kolejne procesy, aby zobaczyć, jak Nextflow radzi sobie z orkiestracją i przepływem danych przez pipeline.

Wróćmy do naszego środowiska Codespaces. Widzisz, że usunąłem wszystkie katalogi .nextflow\* oraz katalogi work, żeby utrzymać porządek. Nie martw się, jeśli nadal masz te pliki z poprzednich części szkolenia.

Będziemy pracować z plikiem o nazwie hello-workflow.nf. Jak poprzednio, reprezentuje on skrypt, który zbudowaliśmy do tej pory, i daje nam czysty punkt wyjścia. I znowu, w sekcji wyjściowej widzimy, że ścieżka to teraz hello_workflow. Więc opublikowane pliki powinny trafiać do innego podkatalogu w Twoim folderze results.

Podsumujmy, gdzie jesteśmy. Mamy jeden proces z jednym wejściem greeting i jednym wyjściem greeting file. Następnie prosty skrypt Bash, który po prostu wykonuje polecenie echo do pliku.

Mamy jedno wejście workflow'u, blok params, w którym mówimy, że oczekuje ścieżki, a wartością domyślną jest data/greetings.csv, czyli ten plik tutaj.

Następnie w samym workflow'ie mamy blok main. Tworzymy kanał, parsujemy CSV na wiersze, a następnie bierzemy pierwszy element każdej tablicy i przekazujemy ten kanał do procesu, który generuje trzy zadania. Publikujemy z workflow'u wyjścia z tego procesu.

I wreszcie, w bloku output mówimy Nextflow'owi, aby opublikował te pliki z tego kanału do katalogu o nazwie hello_workflow i skopiował te pliki zamiast tworzyć dowiązania symboliczne.

## 1. Dodaj drugi krok do workflow'u

Dobrze, w tej części dodamy drugi proces do naszego workflow'u. Weźmiemy wyjścia z procesu sayHello i przetworzymy je w drugim kroku, który zamieni wszystkie litery w tych plikach na wielkie - convertToUppercase.

To tylko prosty przykład, znowu prosta obróbka ciągów znaków, ale pokazuje, jak możemy wykorzystać logikę w workflow'ie.

Użyjemy do tego polecenia bash o nazwie "tr", co jest skrótem od translate. To polecenie Unix, które istnieje od zawsze. Jeśli go nie znasz, nie dziwię się. Myślę, że nigdy nie używałem go przed szkoleniem, ale możesz szybko wypróbować je w terminalu. Jeśli napiszę "echo 'hello world'" i przekieruję do 'tr', a następnie w cudzysłowach podam zakres znaków, więc od A do Z małymi literami, a potem od A do Z wielkimi literami - to mówi: przetłumacz te litery na tamte litery.

Kiedy nacisnę enter, zobaczysz, że wszystko jest teraz wielkimi literami. Bardzo przydatne, jeśli lubisz krzyczeć na ludzi.

To bardzo prosty styl polecenia bash, którego użyjemy w naszym drugim procesie.

## 1.2. Napisz krok zamiany na wielkie litery jako proces Nextflow'a

Więc jeśli wrócę do mojego skryptu, trochę oszukam i po prostu skopiuję kod z dokumentacji szkoleniowej. Ale widzisz dokładnie, co się dzieje.

Mamy tutaj nowy proces. Ten nazwaliśmy convertToUpper, ale moglibyśmy nazwać go jak chcemy.

Mamy jedno wejście path, jak poprzednio. To nie jest kanał wartości, to kanał ścieżki. A następnie jedno wyjście.

W bloku script robimy "cat" na pliku wejściowym. Możemy umieścić to w klamrowych nawiasach, jeśli chcemy, co pobiera tę zmienną. I uruchamiamy to samo polecenie bash w potoku i zapisujemy wyniki do pliku o tej nazwie, która jest przechwytywana przez wyjściową ścieżkę.

Teraz musimy coś zrobić z tym nowym procesem. Więc zejdziemy do workflow'u, gdzie budujemy różną logikę workflow'u, i po pierwszym procesie uruchomimy nasz drugi proces. Więc convertToUpper to nazwa procesu tutaj.

Przyjmuje wejście, więc nie możemy po prostu go wywołać samodzielnie. Chcemy przetworzyć wyjście pierwszego procesu. Więc tak jak zrobiliśmy to z sayHello out, gdzie publikujemy te wyniki, chcemy użyć tych samych wyników tutaj jako wejścia, więc możemy je skopiować i umieścić tam.

Chcemy proces sayHello ".out", a Nextflow wie, że oznacza to prosty pojedynczy rekord wyjściowy tutaj, którym jest ten plik. Więc zostanie przekazany jako wejście do drugiego procesu.

## 1.5. Skonfiguruj publikowanie wyjścia workflow'u

Dobrze. I wreszcie, żebyśmy faktycznie zapisali wyniki tego drugiego procesu, musimy również opublikować je z workflow'u, a następnie zdefiniować je w bloku output, ta sama składnia jak poprzednio. Więc możemy to skopiować i powiedzieć second outputs, lub jak chcesz to nazwać.

Weź nazwę procesu, który nas interesuje, convertToUpper out, a następnie tutaj w bloku output dodaj to i moglibyśmy użyć tych samych atrybutów. Więc również chcemy te pliki w podkatalogu Hello Workflow i również chcemy je skopiować.

Świetnie. Spróbujmy to uruchomić. Więc jeśli wywołam terminal i napiszę "nextflow run hello-workflow.nf", zobaczymy, co zrobi. Sprawdzimy, czy wygląda inaczej niż w poprzednich częściach.

Więc uruchamia Nextflow'a. W dokumentacji mówi się, żeby zrobić to z "-resume", ale usunąłem cały mój katalog work, więc nie zrobiłoby to różnicy. Ale jeśli Ty tego nie zrobiłeś, to również zadziała.

I wygląda prawie dokładnie tak samo. Ale widzisz teraz drugą linię wyjścia tutaj, gdzie widzisz nazwę drugiego procesu, który właśnie dodaliśmy. I rzeczywiście, widzisz, że uruchomił się trzy razy pomyślnie.

Wspaniale. Gdybym miał moje poprzednie katalogi work i zrobiłbym to z "-resume", te byłyby zbuforowane tylko dla pierwszego kroku w pipeline'ie. Ponieważ te wyjścia były dokładnie takie same, więc Nextflow wiedziałby, żeby użyć ich ponownie.

I widzisz, jak możesz użyć -resume do iteracyjnego budowania workflow'u, krok po kroku, jeśli potrzebujesz.

Dobrze, spójrzmy do katalogu results tutaj i zobaczmy, czy zadziałało. Widzimy, że mamy tutaj więcej plików. Mamy nasze oryginalne pliki jak poprzednio z pierwszego procesu. I rzeczywiście, mamy nasze pliki upper i litery są wszystkie wielkie, więc zadziałało. Naprawdę miło to widzieć.

Interesujące jest również sprawdzenie wnętrza tych katalogów work. Jak poprzednio, hash tutaj odpowiada katalogom work. Więc jeśli spojrzę do "ls work" i rozwinę to, zobaczymy różne pliki tutaj.

Widzimy plik wyjściowy z pierwszego procesu, który został pobrany tutaj jako wejście. I widzimy nowy plik wyjściowy, który został wygenerowany.

Teraz jeśli zrobię to z "-la", aby wyświetlić wszystkie pliki, zobaczymy kilka więcej rzeczy. Po pierwsze, zobaczysz, że ten plik jest faktycznie dowiązaniem symbolicznym do pierwszego procesu. To jest zasadniczo zawsze dowiązanie symboliczne, jeśli może być, aby zaoszczędzić miejsce na dysku. Nie publikujemy tutaj plików i po prostu odwołuje się do tego pliku z pierwszego zadania do drugiego zadania, tak że wszystko jest zamknięte w tym jednym katalogu roboczym i bezpieczne oraz odizolowane od wszystkiego innego.

I to musi tam być, ponieważ jeśli spojrzymy na plik .command.sh, więc jeśli zrobię "cat work/b8/56\*", zobaczysz, że części pliku tutaj są względne, więc robi cat tego pliku wejściowego, który został dowiązany symbolicznie do tego samego katalogu roboczego.

Tak więc tak będzie wyglądał każdy katalog work. Kiedy spojrzysz na to w Nextflow'ie, będziesz miał wszystkie pliki wejściowe przygotowane w tym katalogu roboczym. A następnie będziesz miał również wszelkie pliki wyjściowe, które zostały utworzone. Więc to świetnie. Wygląda tak, jak oczekiwaliśmy.

## 2.1. Zdefiniuj polecenie zbierania i przetestuj je w terminalu

Dobrze, wróćmy do naszego workflow'u. Jaki jest następny krok, który chcemy zrobić?

Mamy teraz dwa procesy i pobierają one ten jeden plik CSV, parsują go i dzielą. A następnie mamy trzy zadania dla każdego z tych procesów i Nextflow obsługuje paralelizację tego wszystkiego, więc wszystko działa równolegle, gdzie to możliwe.

Ten sposób dzielenia pracy, aby uruchamiać rzeczy równolegle, jest bardzo powszechny. A odwrotnością tego jest zebranie wszystkiego z powrotem. Więc to zrobimy z naszym ostatnim procesem w workflow'ie - będziemy mieli tutaj trzeci, który weźmie te trzy różne wyjścia i połączy je wszystkie w jeden plik.

Możemy to zrobić dość prosto w terminalu, żeby poczuć, jak to będzie wyglądać.

Jeśli przejdę do folderu results. Więc "cd results/hello_workflow/", i mamy tutaj wszystkie pliki UPPER. Mogę po prostu użyć "cat", którego używamy do wydrukowania zawartości tego pliku, i możesz podać wiele plików do "cat" i będzie czytał jeden po drugim.

Więc mogę powiedzieć "UPPER-\*", co daje mi tę samą listę trzech nazw plików z rozwinięciem Bash. I mogę powiedzieć combined.txt. Myślę, że w dokumentacji wymienia dokładne nazwy plików, ale robi to samo.

Teraz, jeśli użyję "cat combined.txt", zobaczymy, że mamy zawartość plików wszystkich trzech tych plików.

Więc to zasadniczo wszystko, co ten proces będzie robił - spróbujemy dać mu wszystkie różne pliki wyjściowe z poprzedniego procesu w jednym zadaniu procesu, a następnie zrobimy "cat" ich razem i zapiszemy plik wyjściowy.

## 2.2. Utwórz nowy proces do wykonania kroku zbierania

Dobrze, więc dodajmy nasz nowy proces. Wkleję to z materiałów szkoleniowych i widzisz, że zostawił nam trochę ćwiczenia dla czytelnika tutaj z tymi znakami zapytania. Ale widzisz ogólny zarys procesu to zasadniczo to, co właśnie zrobiliśmy w terminalu, gdzie robimy "cat" wielu plików wejściowych i zapisujemy to do pliku wyjściowego tutaj o nazwie collected, a następnie wyjście oczekuje tej pojedynczej ścieżki ponownie.

Więc potrzebujemy tutaj jakiegoś wejścia i będzie to zestaw ścieżek. Więc znowu definiujemy kanał wejściowy path i nazwijmy go input_files. Teraz, to poprzednio dawało nam pojedynczą ścieżkę tutaj, ale ścieżka może również mieć wiele plików tutaj, nawet jeśli to nadal pojedyncza deklaracja.

Skopiuję to tutaj, ponieważ chcemy zrobić "cat" tych plików. I możesz pomyśleć, że mamy tutaj jakieś problemy z drukowaniem tablicy lub czymś takim, ale Nextflow jest generalnie dość rozsądny, jeśli chodzi o to. I jeśli dostanie kanał z wieloma plikami w nim w ten sposób, umieści je wszystkie razem z separatorami spacji. Więc da nam to poprawną składnię.

To świetnie. Więc teraz podłączmy nasz nowy proces. Przejdę do workflow'u. Powiem combine the outputs, nowa nazwa procesu, i tak samo jak poprzednio wezmę ten poprzedni proces, convertToUpper i zrobię ".out".

Świetnie. Wypróbujmy to i zobaczmy, czy działa w terminalu. Jeśli po prostu wrócę o kilka katalogów w górę, a następnie ponownie uruchomię polecenie Nextflow, zobaczymy, co się stanie.

Więc workflow został uruchomiony i teraz widzisz, że mamy trzy różne nazwy procesów, co jest świetne. Pierwsze dwa wyglądają tak samo jak poprzednio, a trzeci nowy się uruchamia, co jest dobre.

Jednak jest coś trochę dziwnego tutaj. Chcieliśmy połączyć te pliki wyjściowe w jeden plik, a jednak ten proces, jak widzimy, uruchomił się trzy razy, a nie raz.

Rzeczywiście, jeśli wejdziemy do jednego z tych katalogów work i zrobimy "cat work/" "collected", zobaczymy, że jest tutaj tylko jedno słowo, a nie trzy.

I więc to, co się stało, to że Nextflow kontynuował tę paralelizację tak samo, jak w poprzednich krokach. I ten proces dał nam kanał z trzema elementami, a te trzy elementy kanału zostały przekazane do naszego procesu downstream, który wygenerował trzy zadania procesu.

Zasadniczo próbował zebrać trzy oddzielne razy i za każdym razem miał tylko jeden plik, więc po prostu zrobił cat pojedynczego pliku do wyjścia, i faktycznie możemy to zobaczyć również w pliku .command.sh.

Jeśli zrobię .command.sh, zobaczymy, że ma tylko jedną nazwę pliku tutaj i tylko jeden plik został przygotowany w tym katalogu roboczym.

## 2.3. Dodaj krok zbierania do workflow'u

Więc jakoś musimy powiedzieć Nextflow'owi, aby zebrał wszystkie te wyjścia z poprzedniego procesu i dał je temu procesowi downstream jako pojedynczy element kanału, a nie trzy.

Robimy to za pomocą operatora kanału o nazwie _collect_.

To jest super przydatny operator, który zobaczysz w pipeline'ach Nextflow cały czas. To jest kanał tutaj, ten kanał wyjściowy, taki sam jak ten, który stworzyliśmy na górze. I więc możemy dołączać do niego operatory kanału tak samo, jak robiliśmy to wcześniej. Możemy po prostu zrobić kropkę, a następnie w tym przypadku collect, nawiasy.

I to wszystko, czego potrzebujemy. To będzie wtedy manipulować tym kanałem, zanim zostanie przekazany do tego procesu.

Jeśli chcesz zobaczyć, co się z nim dzieje, możemy również to wyświetlić tutaj. Więc tutaj, to nie jest związane z uruchamianiem tego procesu w ogóle, więc mogłem umieścić to w dowolnym momencie po uruchomieniu tego procesu. Ale bierzemy ten sam kanał wyjściowy i patrzymy na niego z .view, a następnie patrzymy na niego ponownie z .collect.view.

I kiedy to uruchomimy, pokaże nam dwie różne struktury tego kanału, przed i po collect. Więc spróbujmy tego teraz. Dobrze, właśnie trochę pomniejszyłem, ponieważ niektóre wyjścia są dość długie, ale jeśli uruchomię pipeline, zobaczymy, czy działa.

Mam nadzieję, że trzeci proces uruchomi się tylko raz, ponieważ zbiera wyjścia i rzeczywiście, widzisz collectGreetings jako jeden z jednego. Więc uruchomił się tylko jedno zadanie.

A następnie, jeśli spojrzymy na instrukcje view, mamy trzy instrukcje view dla trzech elementów przed, z jedną ścieżką pliku w każdym.

A następnie po tej instrukcji collect, która została uruchomiona tylko raz, ponieważ jest pojedynczy element w tym kanale. I teraz mamy tę listę trzech różnych ścieżek plików.

To dokładnie to, na co liczyliśmy. I widzisz, mam nadzieję, że to jest zasadniczo odwrotność tego operatora "map", którego użyliśmy, aby przejść z tablic CSV do oddzielnych elementów kanału. Teraz bierzemy oddzielne elementy kanału i wkładamy z powrotem do pojedynczej tablicy.

Świetnie, możemy usunąć te instrukcje view. Nie potrzebujemy ich już. Możemy przejść do następnego kroku.

Zanim pójdę dalej i zanim zapomnę, dodam tutaj nową instrukcję publish. Third output. Możesz nazwać to czymś bardziej semantycznym i opisowym w swoim workflow'ie. A następnie dodam to do bloku output ponownie i powiem path 'hello_workflow' mode 'copy'. Tylko po to, aby plik wyjściowy wygenerowany przez ten proces został zapisany w naszym folderze results tutaj.

Tylko żeby szybko sprawdzić, czy to działa. Powinno być teraz trochę czystsze, ponieważ nie mamy tych instrukcji view. I zobaczymy, czy dostaniemy nasz nowy plik wyjściowy tutaj. Jedno z jednego zadanie uruchomiło się, dostaliśmy nowy plik o nazwie collected, i teraz mamy wszystkie trzy te słowa. Fantastycznie. Co dalej?

## 3. Przekaż dodatkowe parametry do procesu

Dobrze. Następnie przyjrzymy się obsłudze wielu wejść do jednego procesu. Do tej pory widzisz, że wszystkie nasze procesy przyjmują tylko jedną rzecz jako wejście. Wszystkie mają pojedynczą linię pod swoim wejściem.

Zademonstrujemy to, pozwalając Nextflow'owi określić inny identyfikator partii, tak że może uruchomisz ten workflow wiele razy i możesz podać mu inny identyfikator partii za każdym razem.

Po prostu dodam drugą linię w wejściu tutaj dla collectGreetings. I nazwę to "val", ponieważ to jest ciąg znaków. Teraz to jest wartość, a nie ścieżka, i nazwę to "batch_name".

Następnie edytuję skrypt tutaj, aby użyć tej zmiennej, i spróbuję umieścić to w tym samym miejscu co materiał szkoleniowy. Więc umieszczę to w środku tej ścieżki pliku COLLECTED-$\{batch_name\}-output.

Jeszcze nie skończyliśmy. Pamiętaj, że musimy powiedzieć Nextflow'owi, jakie będą nazwy plików wyjściowych. Więc musimy również zrobić to samo tutaj: COLLECTED-$\{batch_name\}-output.txt".

Fantastycznie. Nextflow teraz otrzymuje drugie wejście zmiennej i interpoluje to do skryptu i wyjścia.

Ostatnia rzecz, teraz musimy znaleźć, gdzie to jest wywoływane, i musimy przekazać drugie wejście do procesu. To jest tak samo jak każde inne wejście do funkcji w każdym innym języku.

Tak jak zrobiliśmy wcześniej w szkoleniu, użyję tutaj specjalnego "params" i nazwiemy to "params.batch", abyśmy mogli mieć opcję CLI --batch. I teraz widzisz, że nasz proces tutaj ma dwa oddzielne wejścia oddzielone przecinkami, które są przekazywane.

Naprawdę ważne jest, aby uzyskać właściwą kolejność, więc kolejność argumentów tutaj dla kanału, a następnie parametru musi się zgadzać z kanałem i batch name tam. To jest tylko dopasowanie pozycyjne.

Dobrze. Mogę teraz uruchomić ten pipeline od razu z --batch, ale najpierw zróbmy właściwą rzecz i zdefiniujmy to w wejściu tutaj w Params. Więc dodam to do batch, a następnie powiemy, że to jest ciąg znaków i dajmy mu wartość domyślną. Więc nazwijmy to po prostu batch. Dobrze? Teraz spróbujmy uruchomić workflow.

--batch Trio. Myślę, że tak mówi w materiale szkoleniowym, ale moglibyśmy użyć dowolnego ciągu znaków. I mam nadzieję, że zobaczymy ten plik wyjściowy wyników tutaj.

I rzeczywiście, COLLECTED-trio-output - to zadziałało poprawnie. Zmieniło nazwę naszego pliku. I możesz sobie teraz wyobrazić, że to jest przydatne, ponieważ jeśli uruchomię to ponownie z inną nazwą partii, jak replicate_two, to da nam inną nazwę partii tutaj.

I nie nadpisze wtedy plików wyjściowych w tym przypadku. Więc to miłe.

## 4. Dodaj wyjście do kroku zbierającego

Dobrze, więc mamy teraz wiele wejść do naszego procesu tutaj. Ale co się stanie, jeśli chcemy utworzyć wiele wyjść? Naszym przykładem tutaj jest to, że utworzymy raport dla tego procesu, mówiący, ile plików zostało zebranych.

I zrobimy to za pomocą polecenia echo tutaj. Więc możemy powiedzieć echo. Było, skopiuję to z materiału szkoleniowego, żebyś nie musiał patrzeć, jak to piszę.

Było $\{count_greetings\} powitań w tej partii, i zapiszemy to do nowego pliku teraz o nazwie $\{batch_name\}, więc ta sama zmienna, możemy użyć jej tyle razy, ile chcemy, report.txt.

## 4.1.1. Policz liczbę zebranych powitań

Musimy to jakoś obliczyć. Moglibyśmy zrobić tę logikę w skrypcie Bash, gdybyśmy chcieli, używając logiki Bash. Jednak możemy również po prostu wykonać skryptowanie bezpośrednio w kodzie Nextflow'a, o ile jest to w bloku script w procesie i powyżej sekcji w cudzysłowach.

Wszystko tutaj nie będzie zawarte w ostatecznym wyrenderowanym skrypcie i zostanie po prostu wykonane przez Nextflow'a, gdy renderuje zadanie.

Więc tutaj po prostu robimy trochę logiki. Tworzymy nową zmienną o nazwie count_greetings. Bierzemy kanał plików wejściowych tutaj i wywołujemy na nim .size().

Dobrze, ta funkcja da mi tutaj liczbę do tej zmiennej, i teraz nasze ostrzeżenie zniknęło, ponieważ ta zmienna jest definiowana.

Dobrze, więc tworzymy ten drugi plik w katalogu work, ale musimy powiedzieć Nextflow'owi, aby oczekiwał go jako opublikowane wyjście tego procesu. Więc robimy to dokładnie tą samą składnią, co dla pierwszego pliku.

Mówimy path, ponieważ to jest, znowu, moglibyśmy publikować tutaj zmienną, gdybyśmy chcieli z "val", ale powiemy "path". A następnie oczekiwana nazwa pliku. Zauważ, że nie jest tutaj podświetlona. To dlatego, że użyłem pojedynczych cudzysłowów. Muszę użyć podwójnych cudzysłowów.

## 4.1.2. Wyemituj plik raportu i nazwij wyjścia

Dobrze, to świetnie. I moglibyśmy teraz zacząć uzyskiwać dostęp do tych wyjść tutaj tak samo, jak zrobiłem tutaj. Ale teraz to jest tablica różnych obiektów, więc mógłbym zrobić collectGreetings.out[0], aby uzyskać pierwszy, lub jeden, aby uzyskać drugi, który jest naszym nowym raportem.

Ale naprawdę nie lubię tego robić, ponieważ dość łatwo pomylić się w liczeniu indeksów. I siedzisz tam, licząc linie dużo i dodajesz nowe wyjście i nagle wszystko się psuje. Więc

o wiele lepiej jest odwoływać się do wszystkiego po nazwie. I możemy to zrobić za pomocą specjalnego klucza tutaj o nazwie "emit".

Więc możemy nazwać to, jak chcemy. Nazwijmy to emit outfile i emit reports. Jeśli zdefiniujesz te i możesz to zrobić na jednym lub wielu, to zależy od Ciebie. Teraz mogę zejść tutaj i zamiast tego mogę zrobić kropkę out kropka reports i po prostu wywołać to po nazwie, co jest o wiele łatwiejsze do zrozumienia Twojego kodu, gdy go czytasz, i jest bezpieczniejsze dla zmian w kodzie.

Dodałem .out.report tutaj, ale właściwie muszę mieć dwa różne wyjścia publikowane. Więc zmienię nazwę na coś bardziej interesującego, jak collected i report i czy tak to nazwałem? Nazwałem to out file, przepraszam. Więc ta nazwa emit tutaj outfile i report, ponieważ publikujemy dwa różne kanały wyjściowe i więc musimy odwołać się do obu z nich w bloku publish.

Następnie musimy również zdefiniować je w bloku output. Więc zmieniłem nazwę tego collected, i znowu, dla reports, trochę rozwlekłe tutaj, ale jest to naprawdę przydatne, gdy wchodzisz, aby przeczytać nowy workflow, aby zobaczyć wszystkie różne wyjścia tutaj, wszystkie różne kanały wymienione obok siebie, i są sposoby, aby uczynić to mniej rozwlekłym, czego dotknę później.

Dobrze, spróbujmy to i uruchommy nasz workflow i zobaczmy, co się stanie.

Mam nadzieję, że teraz powinien działać zasadniczo tak samo, jak poprzednio. I dostaniemy nowy plik wyjściowy tutaj o nazwie replicate_two, report. I proszę. Otworzył się i mówi, że są trzy powitania w partii, co jest tym, czego oczekiwaliśmy, więc jest idealnie.

Jeśli wejdę do katalogu work tutaj, tylko żeby Ci udowodnić, że zostało wykonane w kodzie Nextflow'a, a nie w skrypcie bash, mogę zrobić cat work/ command.sh, i zobaczysz tutaj, że po prostu wyświetla ten ciąg bezpośrednio. Było trzy powitania w tej partii, i więc ta zmienna została interpolowana przez Nextflow'a. Została obliczona w bloku script, zanim napisał plik .command.sh. Więc wynikowe obliczenie zmiennej jest zasadniczo zakodowane na stałe w tym, zanim zostanie wykonane w Twoim środowisku obliczeniowym w tym przypadku.

I więc widzisz to rozdzielenie między blokiem script tutaj i wszystkim powyżej niego. Mam nadzieję, że to ma sens.

## Podsumowanie i quiz

Dobrze, to koniec tej części Hello Nextflow. Więc jak poprzednio, sprawdź quiz. Zrób to na stronie internetowej lub w CLI, przejdź przez niektóre pytania i po prostu sprawdź, czy zrozumiałeś część materiału, który omówiliśmy. Zobacz, czy jest tam coś, co podkreśla coś, czego mogłeś nie zrozumieć. Niezbyt wiele pytań. Miłe i łatwe do zrobienia. Lub możesz to zrobić na stronie internetowej tutaj również.

I zrób sobie małą przerwę, małą przechadzkę i wróć i dołącz do nas w czwartej części Hello Nextflow, gdzie porozmawiamy o modułach. Dziękuję bardzo.
