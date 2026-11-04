# Część 3: Hello Workflow - Transkrypcja wideo

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Tłumaczenie wspomagane przez AI - [dowiedz się więcej i zasugeruj ulepszenia](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

<div class="video-wrapper">
  <iframe width="560" height="315" src="https://www.youtube.com/embed/_aO56V3iXGI?si=Irl9nAQniDyICp2b&amp;list=PLPZ8WHdZGxmWKozQuzr27jyMGqp9kElVK&amp;cc_load_policy=1" title="YouTube video player" frameborder="0" allow="accelerometer; autoplay; clipboard-write; encrypted-media; gyroscope; picture-in-picture; web-share" referrerpolicy="strict-origin-when-cross-origin" allowfullscreen></iframe>
</div>

!!!note "Ważne uwagi"

    Ta strona zawiera tylko transkrypcję. Pełne instrukcje krok po kroku znajdziesz w [materiałach szkoleniowych](../03_hello_workflow.md).

    Numery sekcji pokazane w transkrypcji służą wyłącznie celom orientacyjnym i mogą nie obejmować wszystkich numerów sekcji w materiałach.

## Powitanie i podsumowanie

Cześć i witaj ponownie w trzeciej części Hello Nextflow. Ta część nazywa się Hello Workflow i to właśnie tutaj zaczynamy naprawdę uzasadniać nazwę pipeline czy workflow.

Weźmiemy nasz prosty skrypt pipeline'u z jednym procesem i zaczniemy dodawać kolejne procesy, aby zobaczyć, jak Nextflow obsługuje tę orkiestrację i przepływ danych przez pipeline.

Wróćmy do naszego Codespaces. Widzisz, że usunąłem wszystkie katalogi .nextflow\* oraz katalogi work, żeby było czysto. Nie martw się, jeśli nadal masz te pliki z poprzednich części szkolenia.

Będziemy pracować z plikiem o nazwie hello-workflow.nf. Jak poprzednio, reprezentuje on skrypt, który zbudowaliśmy do tej pory, i daje nam czysty punkt wyjścia. I znowu, w sekcji output widzimy, że ścieżka to teraz hello_workflow. Więc opublikowane pliki powinny trafiać do innego podkatalogu w Twoim folderze results.

Podsumujmy, gdzie jesteśmy. Mamy jeden proces z jednym wejściem greeting, jednym wyjściem greeting file. A potem prosty skrypt Bash, który po prostu wykonuje polecenie echo do pliku.

Mamy jedno wejście workflow'u, blok params, gdzie mówimy, że oczekuje ścieżki, a domyślnie to data/greetings.csv, czyli ten plik tutaj.

Następnie w samym workflow mamy blok main. Tworzymy kanał. Parsujemy CSV na wiersze, a potem bierzemy pierwszy element każdej tablicy i przekazujemy ten kanał do procesu, który następnie generuje trzy zadania, i publikujemy z workflow'u wyjścia z tego procesu.

I wreszcie, w bloku output mówimy Nextflow'owi, żeby opublikował te pliki z tego kanału do katalogu o nazwie hello_workflow. I żeby skopiował te pliki zamiast tworzyć dowiązania symboliczne.

## 1. Dodaj drugi krok do workflow'u

Dobrze, w tej części dodamy drugi proces do naszego workflow'u. Weźmiemy wyjścia z procesu sayHello i przetworzymy je w drugim kroku, który zamieni wszystkie litery w tych plikach - convertToUppercase.

To tylko głupi przykład, znowu proste przetwarzanie ciągów znaków, ale pokazuje, jak możemy wykorzystać logikę w workflow.

Użyjemy do tego polecenia bash o nazwie "tr", co jest skrótem od translate. To polecenie Unix, które istnieje od zawsze. Jeśli go nie znasz, nie dziwię się. Nie sądzę, żebym kiedykolwiek używał go przed szkoleniem, ale możesz szybko wypróbować je w terminalu. Jeśli zrobię "echo 'hello world'" i potem przekieruję do 'tr', a następnie w cudzysłowach podam zakres znaków, więc od A do Z, małe litery, a potem chcę od A do Z, wielkie litery. I po prostu mówi: przetłumacz te litery na te litery.

I kiedy nacisnę enter, zobaczysz, że wszystko jest teraz wielkimi literami. Bardzo fajne, jeśli lubisz krzyczeć na ludzi.

To bardzo prosty styl polecenia bash, którego użyjemy w naszym drugim procesie.

## 1.2. Napisz krok zamiany na wielkie litery jako proces Nextflow

Więc jeśli wrócę do mojego skryptu, trochę oszukam i po prostu skopiuję kod z dokumentacji szkoleniowej. Ale widzisz dokładnie, co się dzieje.

Mamy tutaj nowy proces. Ten nazwaliśmy convertToUpper, ale moglibyśmy nazwać go jak chcemy.

Mamy jedno wejście path, jak poprzednio. To nie jest kanał wartości, to kanał ścieżki. A potem jedno wyjście.

W bloku script robimy "cat" na pliku wejściowym. I możemy umieścić to w klamrach, jeśli chcemy. co pobiera tę zmienną. I uruchamiamy to samo polecenie bash w potoku i zapisujemy wyniki do pliku o tej nazwie, a to jest przechwytywane przez wyjściową ścieżkę.

Teraz musimy coś zrobić z tym nowym procesem. Więc zejdziemy do workflow, gdzie budujemy różną logikę workflow'u, i po tym pierwszym procesie uruchomimy nasz drugi proces. Więc convertToUpper to nazwa procesu tutaj.

Przyjmuje wejście, więc nie możemy po prostu go wywołać samodzielnie. Chcemy przetworzyć wyjście pierwszego procesu. Więc tak jak zrobiliśmy z tym, sayHello out, gdzie publikujemy te wyniki. Chcemy użyć tych samych wyników tutaj jako wejścia, więc możemy je skopiować i umieścić tam.

Chcemy proces sayHello ".out", a Nextflow wie, że to oznacza prosty pojedynczy rekord wyjściowy tutaj, którym jest ten plik. Więc zostanie on przekazany jako wejście do drugiego procesu.

## 1.5. Skonfiguruj publikowanie wyjścia workflow'u

Dobrze. I wreszcie, żebyśmy faktycznie zapisali wyniki tego drugiego procesu, musimy je również opublikować z workflow'u, a następnie zdefiniować je w bloku output, ta sama składnia jak poprzednio. Więc możemy to skopiować i powiedzieć second outputs, lub jak chcesz to nazwać.

Weź nazwę procesu, który nas interesuje, convertToUpper out, a potem tutaj w bloku output. Dodaj to i moglibyśmy zrobić te same atrybuty tutaj. Więc chcemy również te pliki w podkatalogu Hello Workflow i chcemy je również skopiować.

Świetnie. Spróbujmy to uruchomić. Więc jeśli wywołam terminal i zrobię "nextflow run hello-workflow.nf", zobaczymy, co zrobi. Zobaczymy, czy wygląda inaczej niż poprzednie części.

Więc uruchamia Nextflow. W dokumentacji mówi, żeby zrobić to z "-resume", ale usunąłem cały mój katalog work, więc nie zrobiłoby to tutaj różnicy. Ale jeśli Ty to zrobiłeś, to też zadziała.

I wygląda prawie dokładnie tak samo. Ale widzisz teraz, że jest druga linia wyjścia tutaj, gdzie widzisz nazwę drugiego procesu, który właśnie dodaliśmy. I rzeczywiście, widzisz, że uruchomił się trzy razy pomyślnie.

Wspaniale. Gdybym miał moje poprzednie katalogi work i zrobiłbym to z "-resume", te byłyby zbuforowane tylko pierwszy krok w pipeline. Bo te wyjścia były dokładnie takie same, więc Nextflow wiedziałby, żeby użyć ich ponownie.

I widzisz, jak możesz użyć -resume do iteracyjnego budowania swojego workflow'u, krok po kroku, jeśli potrzebujesz.

Dobrze, spójrzmy do katalogu results tutaj i zobaczmy, czy zadziałało. Widzimy, że mamy tutaj więcej plików. Mamy nasze oryginalne pliki jak poprzednio z pierwszego procesu. I rzeczywiście, mamy nasze pliki upper i litery są wszystkie wielkie, więc zadziałało. Naprawdę miło to widzieć.

Interesujące jest również sprawdzenie wnętrza tych katalogów work. Jak poprzednio, hash tutaj odpowiada katalogom work. Więc jeśli spojrzę do "ls work", a potem rozwinę to, zobaczymy różne pliki tutaj.

Widzimy plik wyjściowy z pierwszego procesu, który został pobrany tutaj jako wejście. I widzimy nowy plik wyjściowy, który został wygenerowany.

Teraz jeśli zrobię to z "-la", żeby wylistować i pokazać wszystkie pliki, zobaczymy kilka więcej rzeczy. Po pierwsze, zobaczysz, że ten plik jest faktycznie dowiązaniem symbolicznym do pierwszego procesu. To jest zasadniczo zawsze dowiązanie symboliczne, jeśli może być, żeby zaoszczędzić miejsce na dysku. Nie publikujemy tutaj plików i po prostu odwołuje się do tego pliku z pierwszego zadania do drugiego zadania, tak że wszystko jest zamknięte w tym jednym katalogu roboczym, bezpieczne i odizolowane od wszystkiego innego.

I to musi tam być, ponieważ jeśli spojrzymy na plik .command.sh, więc jeśli zrobię "cat work/b8/56\*", zobaczysz, że części pliku tutaj są względne, więc robi cat tego pliku wejściowego, który został dowiązany symbolicznie do tego samego katalogu roboczego.

Więc tak będzie wyglądał każdy katalog work. Kiedy spojrzysz na to w Nextflow, będziesz miał wszystkie pliki wejściowe tam przygotowane w tym katalogu roboczym. A potem będziesz miał również wszystkie pliki wyjściowe, które zostały utworzone. Więc to świetnie. Wygląda tak, jak oczekujemy.

## 2.1. Zdefiniuj polecenie zbierania i przetestuj je w terminalu

Dobrze, wróćmy do naszego workflow'u. Jaki jest następny krok, który chcemy zrobić?

Mamy teraz dwa procesy i pobierają one ten jeden plik CSV, parsują go i dzielą. A potem mamy trzy zadania dla każdego z tych procesów i Nextflow obsługuje paralelizację tego wszystkiego, więc wszystko działa obok siebie, gdzie to możliwe.

Ten sposób dzielenia pracy, żeby uruchamiać rzeczy równolegle, jest bardzo powszechny. A odwrotnością tego jest zebranie wszystkiego z powrotem. Więc to zrobimy z naszym ostatnim procesem w workflow - będziemy mieli tutaj trzeci, który weźmie te trzy różne wyjścia i połączy je wszystkie w jeden plik.

Możemy to zrobić całkiem prosto w terminalu, żeby poczuć, jak to będzie wyglądać.

Jeśli przejdę do folderu results. Więc, "cd results/hello_workflow/", i mamy tutaj wszystkie pliki UPPER. Mogę po prostu użyć "cat", którego używamy do wydrukowania zawartości tego pliku, i możesz podać wiele plików do "cat" i przeczyta je jeden po drugim.

Więc mogę powiedzieć "UPPER-\*", co daje mi tę samą listę trzech nazw plików z rozwinięciem Bash. I mogę powiedzieć combined.txt. Myślę, że w dokumentacji wymienia dokładne nazwy plików, ale robi to samo.

Teraz, jeśli użyję "cat combined.txt", zobaczymy, że mamy zawartość pliku wszystkich trzech tych plików.

Więc to zasadniczo wszystko, co ten proces będzie robił - spróbujemy dać mu wszystkie różne pliki wyjściowe z poprzedniego procesu w jednym zadaniu procesu, a potem zrobimy "cat" na nich razem i zapiszemy plik wyjściowy.

## 2.2. Utwórz nowy proces do wykonania kroku zbierania

Dobrze, więc dodajmy nasz nowy proces. Wkleję to z materiałów szkoleniowych i widzisz, że zostawił nam trochę ćwiczenia dla czytelnika tutaj z tymi znakami zapytania. Ale widzisz ogólny zarys procesu to zasadniczo to, co właśnie zrobiliśmy w terminalu, gdzie robimy "cat" na kilku plikach wejściowych i zapisujemy to do pliku wyjściowego tutaj o nazwie collected, a potem wyjście oczekuje tej pojedynczej ścieżki ponownie.

Więc potrzebujemy tutaj jakiegoś wejścia i będzie to zestaw ścieżek. Więc znowu, definiujemy kanał wejściowy path i nazwijmy go input_files. Teraz, to poprzednio dawało nam pojedynczą ścieżkę tutaj, ale ścieżka może również mieć wiele plików tutaj, mimo że to nadal pojedyncza deklaracja.

Skopiuję to tutaj, ponieważ chcemy zrobić "cat" na tych plikach. I możesz pomyśleć, że mamy tutaj jakieś problemy z drukowaniem tablicy lub czymś takim, ale Nextflow jest generalnie całkiem rozsądny, jeśli chodzi o to. I jeśli dostanie kanał z wieloma plikami w nim jak ten, umieści je wszystkie razem z separatorami spacji. Więc to da nam poprawną składnię.

To świetnie. Więc teraz podłączmy nasz nowy proces. Idę do workflow. Powiem combine the outputs, nowa nazwa procesu, i tak samo jak poprzednio. Wezmę ten poprzedni proces, convertToUpper i zrobię ".out".

Świetnie. Wypróbujmy to i zobaczmy, czy działa w terminalu. Jeśli po prostu wrócę o kilka katalogów w górę, a potem ponownie uruchomię polecenie Nextflow, zobaczymy, co się stanie.

Więc workflow został uruchomiony i teraz widzisz, że mamy trzy różne nazwy procesów, co jest świetne. Pierwsze dwa wyglądają tak samo jak poprzednio, a trzeci nowy się uruchamia, co jest dobre.

Jednak jest coś trochę dziwnego tutaj. Chcieliśmy połączyć te pliki wyjściowe w jeden plik, a jednak ten proces, jak widzimy, uruchomił się trzy razy, nie raz.

Rzeczywiście, jeśli wejdziemy do jednego z tych katalogów work. I zrobimy "cat work/" "collected", to zobaczymy. Jest tutaj tylko jedno słowo, nie trzy.

I więc to, co się stało, to że Nextflow kontynuował tę paralelizację tak samo jak w poprzednich krokach. I ten proces dał nam kanał z trzema elementami, a te trzy elementy kanału zostały przekazane do naszego procesu downstream, który wygenerował trzy zadania procesu.

Zasadniczo próbował zebrać trzy oddzielne razy i za każdym razem miał tylko jeden plik, więc po prostu zrobił cat pojedynczego pliku do wyjścia, i faktycznie, możemy to zobaczyć również w pliku .command.sh.

Jeśli zrobię .command.sh, zobaczymy, że ma tylko jedną nazwę pliku tutaj i tylko jeden plik został przygotowany w tym katalogu roboczym.

## 2.3. Dodaj krok zbierania do workflow'u

Więc jakoś musimy powiedzieć Nextflow, żeby zebrał wszystkie te wyjścia z poprzedniego procesu i dał je temu procesowi downstream jako pojedynczy element kanału, a nie trzy.

Robimy to za pomocą operatora kanału o nazwie _collect_.

To jest super przydatny operator, który zobaczysz w pipeline'ach Nextflow cały czas. To jest kanał tutaj, ten kanał wyjściowy, tak samo jak ten, który stworzyliśmy na górze. I więc możemy dołączać do niego operatory kanału tak samo jak poprzednio. Możemy po prostu zrobić kropkę, a potem w tym przypadku, collect, nawiasy.

I to wszystko, czego potrzebujemy. To następnie zmanipuluje ten kanał, zanim zostanie przekazany do tego procesu.

Jeśli chcesz zobaczyć, co się z nim dzieje, możemy również go wyświetlić tutaj. Więc tutaj, to nie jest związane z uruchamianiem tego procesu w ogóle, więc mogłem umieścić to w dowolnym momencie po uruchomieniu tego procesu. Ale bierzemy ten sam kanał wyjściowy i patrzymy na niego z .view, a potem patrzymy na niego ponownie z .collect.view.

I kiedy to uruchomimy, pokaże nam dwie różne struktury tego kanału, przed i po collect. Więc spróbujmy tego teraz. Dobrze, właśnie trochę oddaliłem, ponieważ niektóre wyjścia są dość długie, ale jeśli uruchomię pipeline, zobaczymy, czy działa.

Mam nadzieję, że trzeci proces uruchomi się tylko raz, ponieważ zbiera wyjścia i rzeczywiście, widzisz collectGreetings jako jeden z jednego. Więc uruchomił się tylko jedno zadanie.

A potem jeśli spojrzymy na instrukcje view, mamy trzy instrukcje view dla trzech elementów przed, z jedną ścieżką pliku w każdym.

A potem po tej instrukcji collect, to zostało uruchomione tylko raz, ponieważ jest pojedynczy element w tym kanale. I teraz mamy tę listę trzech różnych ścieżek plików.

To dokładnie to, na co liczyliśmy. I widzisz, mam nadzieję, że to jest zasadniczo odwrotność tego operatora "map", którego użyliśmy, żeby przejść z tablic CSV do oddzielnych elementów kanału. Teraz bierzemy oddzielne elementy kanału i wkładamy z powrotem do pojedynczej tablicy.

Świetnie, możemy wyczyścić te instrukcje view. Nie potrzebujemy ich już. Możemy przejść do następnego kroku.

Zanim pójdę dalej i zanim zapomnę, dodam tutaj nową instrukcję publish. Third output. Możesz nazwać to czymś bardziej semantycznym i opisowym w swoim workflow. A potem dodam to do bloku output ponownie i powiem path 'hello_workflow' mode 'copy'. Tylko po to, żeby plik wyjściowy wygenerowany przez ten proces został zapisany do naszego folderu results tutaj.

Tylko żeby szybko sprawdzić, czy to działa. Powinno być teraz trochę czystsze, bo nie mamy tych instrukcji view. I zobaczymy, czy dostaniemy nasz nowy plik wyjściowy tutaj. Jedno z jednego zadanie uruchomiło się, dostaliśmy nowy plik o nazwie collected, i teraz mamy wszystkie trzy te słowa. Fantastycznie. Co dalej?

## 3. Przekaż dodatkowe parametry do procesu

Dobrze. Następnie przyjrzymy się obsłudze wielu wejść do jednego procesu. Do tej pory widzisz, że wszystkie nasze procesy przyjmują tylko jedną rzecz jako wejście. Wszystkie mają pojedynczą linię pod swoim wejściem.

Zademonstrujemy to, pozwalając Nextflow określić inny identyfikator partii, tak że może uruchomisz ten workflow wiele razy i możesz podać mu inny ID partii za każdym razem.

Po prostu dodam drugą linię we wejściu tutaj dla collectGreetings. I nazwę to "val", bo to jest ciąg znaków. Teraz to jest wartość, nie ścieżka, i nazwę to "batch_name".

Następnie edytuję skrypt tutaj, żeby użyć tej zmiennej, i spróbuję umieścić to w tym samym miejscu co materiał szkoleniowy. Więc umieszczę to w środku tej ścieżki pliku COLLECTED-$\{batch_name\}-output.

Jeszcze nie skończyliśmy. Pamiętaj, że musimy powiedzieć Nextflow, jakie będą nazwy plików wyjściowych. Więc musimy również zrobić to samo tutaj: COLLECTED-$\{batch_name\}-output.txt".

Fantastycznie. Nextflow teraz dostaje drugie wejście zmiennej i interpoluje to do skryptu i wyjścia.

Ostatnia rzecz, teraz musimy znaleźć, gdzie to jest wywoływane, i musimy przekazać drugie wejście do procesu. To jest tak samo jak każde inne wejście do funkcji w każdym innym języku.

Tak jak zrobiliśmy wcześniej w szkoleniu, użyję specjalnego "params" tutaj, i nazwiemy to "params.batch", żebyśmy mogli mieć opcję CLI --batch. I teraz widzisz, że nasz proces tutaj ma dwa oddzielne wejścia oddzielone przecinkami, które są przekazywane.

Naprawdę ważne jest, żeby uzyskać właściwą kolejność, więc kolejność argumentów tutaj dla kanału, a potem parametru musi się zgadzać. Kanał i batch name tam. To jest tylko dopasowanie pozycyjne.

Dobrze. Mogę teraz uruchomić ten pipeline od razu z --batch, ale najpierw zróbmy właściwą rzecz i zdefiniujmy to we wejściu tutaj w Params. Więc dodam to do batch, a potem powiemy, że to jest ciąg znaków i dajmy mu wartość domyślną. Więc nazwijmy to po prostu batch. Dobrze? Teraz spróbujmy uruchomić workflow.

--batch Trio. Myślę, że tak mówi w materiale szkoleniowym, ale moglibyśmy użyć dowolnego ciągu, jakiego chcemy. I mam nadzieję, że zobaczymy ten plik wyjściowy results tutaj.

I rzeczywiście, COLLECTED-trio-output - to zadziałało poprawnie. Zmieniło nazwę naszego pliku. I możesz sobie teraz wyobrazić, że to jest przydatne, bo jeśli uruchomię to ponownie z inną nazwą partii, jak replicate_two, to da nam inną nazwę partii tutaj.

I nie nadpisze wtedy plików wyjściowych w tym przypadku. Więc to miłe.

## 4. Dodaj wyjście do kroku zbierającego

Dobrze, więc mamy teraz wiele wejść do naszego procesu tutaj. Ale co się stanie, jeśli chcemy utworzyć wiele wyjść? Naszym przykładem tutaj jest to, że utworzymy raport dla tego procesu, mówiący po prostu, ile plików zostało zebranych.

I zrobimy to za pomocą polecenia echo tutaj. Więc możemy powiedzieć echo. Było, skopiuję to z materiału szkoleniowego, żebyś nie musiał patrzeć, jak to wpisuję.

Było $\{count_greetings\} powitań w tej partii, i zapiszemy to do nowego pliku teraz o nazwie $\{batch_name\}, więc ta sama zmienna, możemy użyć jej tyle razy, ile chcemy, report.txt.

## 4.1.1. Policz liczbę zebranych powitań

Musimy jakoś to obliczyć. Moglibyśmy zrobić tę logikę w skrypcie Bash, gdybyśmy chcieli, używając logiki Bash. Jednak możemy również po prostu wykonać skryptowanie bezpośrednio w kodzie Nextflow, o ile jest to w bloku script w procesie i powyżej sekcji w cudzysłowach.

Wszystko tutaj nie będzie zawarte w ostatecznym wyrenderowanym skrypcie i zostanie po prostu wykonane przez Nextflow, gdy renderuje zadanie.

Więc tutaj po prostu robimy trochę logiki. Tworzymy nową zmienną o nazwie count_greetings. Bierzemy kanał input files tutaj i wywołujemy na nim .size().

Dobrze, ta funkcja da mi tutaj liczbę do tej zmiennej, i teraz nasze ostrzeżenie zniknęło, ponieważ ta zmienna jest definiowana.

Dobrze, więc tworzymy ten drugi plik w katalogu work, ale musimy powiedzieć Nextflow, żeby oczekiwał go jako opublikowane wyjście tego procesu. Więc robimy to dokładnie tą samą składnią, co dla pierwszego pliku.

Mówimy path, bo to jest, znowu, moglibyśmy publikować tutaj zmienną, gdybyśmy chcieli z "val", ale powiemy "path". A potem oczekiwana nazwa pliku. Zauważ, że nie jest tutaj podświetlona. To dlatego, że użyłem pojedynczych cudzysłowów. Muszę użyć podwójnych cudzysłowów.

## 4.1.2. Wyemituj plik raportu i nazwij wyjścia

Dobrze, to świetnie. I moglibyśmy teraz zacząć uzyskiwać dostęp do tych wyjść tutaj tak samo jak tutaj. Ale teraz to jest tablica różnych obiektów, więc mogłem zrobić collectGreetings.out[0], żeby dostać pierwsze, lub jeden, żeby dostać drugie, którym jest nasz nowy raport.

Ale naprawdę nie lubię tego robić, bo dość łatwo pomylić się w liczeniu indeksów. I siedzisz tam licząc linie dużo i dodajesz nowe wyjście i nagle wszystko się psuje. Więc

o wiele lepiej jest odwoływać się do wszystkiego po nazwie zamiast tego. I możemy to zrobić za pomocą specjalnego klucza tutaj o nazwie "emit".

Więc możemy nazwać to jak chcemy. Nazwijmy to emit outfile, i emit reports. Jeśli zdefiniujesz te i możesz to zrobić na jednym lub wielu, to zależy od Ciebie. Teraz mogę zejść tutaj i zamiast tego mogę zrobić kropka out kropka reports i po prostu wywołać to po nazwie, co jest o wiele łatwiejsze do zrozumienia Twojego kodu, kiedy go czytasz, i jest bezpieczniejsze dla zmian w kodzie.

Dodałem .out.report tutaj, ale faktycznie muszę mieć dwa różne wyjścia publikowane. Więc zmienię nazwę na coś ciekawszego jak collected i report i czy tak to nazwałem? Nazwałem to out file, przepraszam. Więc ta nazwa emit tutaj outfile i report. bo publikujemy dwa różne kanały wyjściowe i więc musimy odwołać się do obu w bloku publish.

Następnie musimy również zdefiniować je w bloku output. Więc zmieniłem nazwę tego collected, i znowu, dla reports, trochę rozwlekłe tutaj, ale jest naprawdę przydatne, kiedy wchodzisz, żeby przeczytać nowy workflow, zobaczyć wszystkie różne wyjścia tutaj, wszystkie różne kanały wymienione obok siebie, i są sposoby, żeby to uczynić mniej rozwlekłym, czego dotknę później.

Dobrze, spróbujmy to uruchomić i zobaczymy, co się stanie.

Mam nadzieję, że teraz powinno działać zasadniczo tak samo jak poprzednio. I dostaniemy nowy plik wyjściowy tutaj o nazwie replicate_two, report. I proszę. Otworzył się i mówi, że są trzy powitania w partii, co jest tym, czego oczekiwaliśmy, więc jest idealnie.

Jeśli wejdę do katalogu work tutaj tylko po to, żeby Ci udowodnić, że zostało wykonane w kodzie Nextflow, a nie w skrypcie bash, mogę zrobić cat work/ command.sh, i zobaczysz tutaj, że po prostu wypisuje ten ciąg bezpośrednio. Było trzy powitania w tej partii, i więc ta zmienna została zinterpolowana przez Nextflow. Została obliczona w bloku script, zanim napisał plik .command.sh. Więc wynikowe obliczenie zmiennej jest zasadniczo zakodowane na stałe w tym, zanim zostanie wykonane w Twoim środowisku obliczeniowym w tym przypadku.

I więc widzisz to rozdzielenie między blokiem script. Tutaj i wszystkim powyżej niego. Mam nadzieję, że to ma sens.

## Podsumowanie i quiz

Dobrze, to koniec tej części Hello Nextflow. Więc jak poprzednio, idź i sprawdź quiz. Zrób to na stronie internetowej lub w CLI, przejdź przez niektóre pytania i po prostu sprawdź, czy zrozumiałeś część materiału, który omówiliśmy. Zobacz, czy jest tam coś, co podkreśla coś, czego mogłeś nie zrozumieć. Niezbyt wiele pytań. Miłe i łatwe do zrobienia. Lub możesz to zrobić na stronie internetowej tutaj również.

I zrób sobie małą przerwę, małą przechadzkę i wróć i dołącz do nas w części czwartej Hello, Nextflow, gdzie porozmawiamy o modułach. Dziękuję bardzo.
