# Część 2: Hello Channels - Transkrypcja

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Tłumaczenie wspomagane przez AI - [dowiedz się więcej i zasugeruj ulepszenia](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

<div class="video-wrapper">
  <iframe width="560" height="315" src="https://www.youtube.com/embed/lJ41WMMm44M?si=xCItHLiOQWqoqBB9&amp;list=PLPZ8WHdZGxmXiHf8B26oB_fTfoKQdhlik" title="YouTube video player" frameborder="0" allow="accelerometer; autoplay; clipboard-write; encrypted-media; gyroscope; picture-in-picture; web-share" referrerpolicy="strict-origin-when-cross-origin" allowfullscreen></iframe>
</div>

!!!note "Ważne informacje"

    Ta strona zawiera jedynie transkrypcję. Aby uzyskać pełne instrukcje krok po kroku, wróć do [materiałów szkoleniowych](../02_hello_channels.md).

    Numery sekcji pokazane w transkrypcji są podane jedynie orientacyjnie i mogą nie obejmować wszystkich numerów sekcji w materiałach.

## Powitanie

Cześć, witamy w drugiej części Hello Nextflow.

Ten rozdział nazywa się Hello Channels. Będziemy mówić o tej fundamentalnej części Nextflow.

Channels to elementy, które łączą różne kroki w Twoim pipeline, sposób, w jaki Twoje dane i logika przepływają przez Twój workflow.

Dobra, zanurzmy się w to.

Zacznijmy od przejścia do training.nextflow.io

Hello Nextflow w pasku bocznym i kliknięcia w część drugą, Hello Channels.

Wszystkie materiały są tutaj napisane, więc możesz podążać we własnym tempie i nadrobić wszystko, co mogłeś przegapić.

Gdy już otworzysz stronę internetową, możesz załadować Codespaces i będziemy kontynuować od miejsca, w którym zakończyliśmy ostatni rozdział.

## 0. Rozgrzewka: Uruchom hello-channels.nf

W tym rozdziale będziemy edytować inny plik. Ten nazywa się Hello Channels, więc możesz go znaleźć w pasku bocznym, kliknij dwukrotnie, aby otworzyć.

Teraz, jeśli właśnie przyszedłeś z rozdziału pierwszego, ten plik będzie Ci bardzo znajomy. Punkt wyjścia tutaj to w zasadzie miejsce, w którym zakończyliśmy rozdział pierwszy, z naszym processem o nazwie sayHello, naszym wejściem, wyjściem, naszym publishDir i naszym params.greeting, oraz naszym prostym workflow.

Zaczynamy od nowego pliku, więc to równe warunki dla wszystkich, ale możesz kontynuować z poprzednim plikiem, jeśli wolisz.

Zauważ, że usunąłem również wszystkie pliki .nextflow\* i katalogi work tutaj, po prostu aby był czysty punkt wyjścia. Nie ma znaczenia, czy to zrobisz, czy nie, to zależy od Ciebie.

Dobra. Zacznijmy od sprawdzenia, czy ten pipeline nadal działa zgodnie z naszymi oczekiwaniami. Wywołam terminal tutaj.

Wpiszę "nextflow run hello-channels.nf" i nacisnę enter.

To uruchomi ten mały workflow, uruchomi nasz krok sayHello, wygeneruje katalog work z tym hashem, i tutaj jest nasz folder results, a tam jest nasz plik wyjściowy, dokładnie tak, jak oczekiwaliśmy od naszego domyślnego params.greeting.

To wspaniałe. Dokładnie to samo co w rozdziale pierwszym, działa zgodnie z oczekiwaniami.

## 1. Dostarczanie zmiennych wejściowych przez channel jawnie

W rozdziale pierwszym faktycznie już używałeś channels, po prostu nie zdawałeś sobie z tego sprawy. Gdy określiliśmy tutaj string, Nextflow automatycznie utworzył dla nas channel wokół tego stringa, po prostu dlatego, że wiedział, że wywołujemy process, więc potrzebowaliśmy input channel.

Pierwszą rzeczą, którą zrobimy, jest uczynienie tego jawnym poprzez faktyczne wypisanie samego channel.

## 1.1. Utwórz input channel

Więc przejdę do workflow tutaj na dole skryptu i powiem greeting_ch. To jest konwencja, której często używamy w kodzie Nextflow, aby mieć podkreślenie ch na końcu nazwy zmiennej, gdy jest to channel, po prostu aby łatwo było zidentyfikować, że to jest channel, ale nie musisz tego robić. Równa się channel of Hello Channels.

To, czego właśnie użyliśmy, nazywa się "Channel Factory" w języku Nextflow. To jest ta rzecz tutaj, ustawiamy tę zmienną na nowy channel, a ta fabryka channel tutaj tworzy dla nas channel w określony sposób.

Istnieje kilka różnych fabryk channel, które ma Nextflow, aby tworzyć channels z różnych typów wejść. Dot of jest najprościejszy i po prostu przyjmuje wszystkie stringi, które mu przekażemy.

Zauważ, że gdy najadę na te słowa w VS Code, rozszerzenie Nextflow pokazuje mi popup wyjaśniający, co robi ta składnia, a na dole tego okna popup jest również tekst "read more".

Jeśli kliknę to, otworzy dokumentację Nextflow w nowej zakładce i zabierze mnie bezpośrednio do dokumentacji dla tej konkretnej rzeczy. W tym przypadku dla channel.of.

## 1.2. Dodaj channel jako wejście do wywołania procesu

Zauważ, że rozszerzenie również daje nam ostrzeżenie mówiące, że utworzyliśmy tutaj nowy channel, ale nic go nie używa.

Więc naprawmy to. Wezmę nową nazwę channel i zastąpię ten params.greeting naszym nowym channel.

Zauważ, że nie używamy już teraz flagi wiersza poleceń --greeting, params.greeting nie jest używany, wracamy do zakodowania na stałe tego stringa. To w porządku. Staram się po prostu utrzymać rzeczy prostymi. Wrócimy później i użyjemy params ponownie.

## 1.3. Uruchom ponownie polecenie workflow

Dobra, sprawdźmy tylko, czy to działa. Wywołam terminal i zauważ ponownie. Nextflow run hello channels. Sprawdź output.txt, i oto jest.

Świetny, trochę nudny przykład, robiący dokładnie to samo co wcześniej, ale teraz przynajmniej logika jest bardziej czytelna. Jesteśmy jawni w pisaniu nowego channel.

Właściwie właśnie napisaliśmy więcej kodu, aby zrobić to samo. Ale to zacznie mieć więcej sensu, gdy staniemy się nieco bardziej skomplikowani w sposobie tworzenia naszych channels.

## 2. Zmodyfikuj workflow, aby działał na wielu wartościach wejściowych

Dobra, uczyńmy to nieco ciekawszym. Bardzo rzadko chcesz uruchomić pipeline Nextflow na pojedynczym wejściu, więc dajmy mu kilka wejść.

## 2.1. Załaduj wiele powitań do input channel

Z dokumentacji tutaj. Skopiuję te różne stringi, trzy z nich. Hello, Bonjour, Olà. Oh, mam nadzieję. Copilot sugeruje kilka innych. Więc zatabulujmy i wprowadźmy je.

Dokumentacja Nextflow tutaj mówi nam, że możemy przekazać wiele wartości do tego operatora, więc powinno działać, ale wypróbujmy to i zobaczmy, co się stanie.

## 2.1.2. Uruchom polecenie i spójrz na wyjście logów

Cóż. Tak i nie. Zobaczmy. Mówi, że pięć z pięciu zadań zostało uruchomionych tutaj, ale pokazuje nam tylko jeden hash, co jest trochę dziwne. To w porządku. Wszystko jest zgodne z oczekiwaniami tutaj. Domyślnie Nextflow używa specjalnego typu wyjścia do terminala zwanego kodami kontrolnymi ANSI, co oznacza, że nadpisuje pewne linie, aby dać ładny skompresowany widok wszystkich różnych procesów, które są uruchamiane.

Ma to o wiele większy sens, gdy masz większe workflows i uruchamiasz setki lub tysiące różnych próbek. Po prostu możesz wygenerować tak wiele wyjścia na terminalu, że niemożliwe jest spojrzenie na nie, podczas gdy ten aktualizujący się widok daje ci postęp w czasie rzeczywistym.

## 2.1.3. Uruchom polecenie ponownie z opcją -ansi-log false

Jeśli chcesz, możesz uruchomić to ponownie, a tym razem użyję dodatkowego argumentu rdzenia Nextflow z pojedynczym myślnikiem mówiącym, "-ansi-log false". To używa poprzedniej wersji wyjścia logów Nextflow. I tutaj możesz zobaczyć wszystkie indywidualne procesy, które zostały uruchomione.

To zależy od Ciebie, czy to zrobisz, czy nie. Wyjście z Nextflow jest dokładnie takie samo w obu przypadkach.

## 2.2. Upewnij się, że nazwy plików wyjściowych będą unikalne

Dobra, spójrzmy więc na pliki wyjściowe, następnie przejdziemy do results. Ale jest tylko pojedynczy plik wyjściowy. Co się stało? Widzieliśmy, że process był uruchamiany wiele razy. Możemy przejść do katalogu work i zobaczyć wszystkie różne hashe, wszystkie zadania zostały wykonane prawidłowo. Ale jeśli pamiętasz w naszym procesie tutaj, zapisujemy wszystko do pliku output.txt, a następnie publikujemy to do tego katalogu.

Więc ten sam plik został utworzony pięć razy, a następnie został nadpisany pięć razy. I po prostu mamy to, które zadanie wykonało się ostatnie.

## 2.2.1. Skonstruuj dynamiczną nazwę pliku wyjściowego

Sposób, w jaki to naprawiamy, to użycie dynamicznej nazwy pliku wyjściowego. Tutaj już mamy zmienną o nazwie greeting w procesie, więc możemy użyć jej w nazwie pliku wyjściowego. Kopiuję to i robię $greeting-output.txt.

Otoczę to w cudzysłowy, po prostu żeby bash nie pomylił się przez jakiekolwiek spacje, które mogą się tu pojawić. A następnie wezmę tę samą nazwę pliku i zaktualizuję tutaj wyjście.

To jest naprawdę ważne, aby wyjście pasowało do tego, ponieważ w przeciwnym razie ten plik nie zostanie znaleziony i Nextflow się zawiesi.

Zamierzam dokonać jeszcze jednej naprawdę ważnej edycji, czyli zamienię te pojedyncze cudzysłowy na podwójne cudzysłowy. Zauważ, że kolor kodu zmienił się, gdy to zrobiłem. Ta zmienna jest rozwijana tylko wtedy, gdy używamy podwójnych cudzysłowów. Jeśli użyję tutaj pojedynczych cudzysłowów, jest używana jako wartość literalna, i otrzymałbym pojedynczy plik o nazwie $greeting-output, co nie jest tym, czego chcę.

## 2.2.2. Uruchom workflow

Więc wróćmy do podwójnych cudzysłowów i spróbujmy.

Po prostu zamierzam uporządkować mój katalog przed rozpoczęciem, więc będzie łatwo zobaczyć nowe pliki. Zamierzam usunąć wszystko o nazwie .nextflow, work i results.

I zamierzam uruchomić to polecenie Nextflow ponownie i zobaczmy, jakie pliki zostały utworzone. Więc uruchamia pięć procesów tam. Jeśli bardzo uważnie patrzyłeś, mogłeś zobaczyć, że ta linia się aktualizowała podczas działania.

I teraz możemy przejść do katalogu results, i rzeczywiście, mamy pięć różnych wyjść, i wszystkie są poprzedzone różnymi powitaniami.

Jeśli otworzę każdy z nich, zobaczymy, że każdy zawiera odpowiednie powitanie. Fantastycznie. To jest to, czego chcemy.

## 3. Użyj operatora do przekształcenia zawartości channel

Dobra, więc teraz wiemy, czym są channels i wiemy, czym są fabryki channel. A co z operatorami? To kolejny termin dla części języka Nextflow, który jest serią funkcji, które pozwalają nam operować na channels, aby zrobić z nimi pewne rzeczy. Nextflow zawiera zestaw operatorów, które pozwalają nam manipulować channels na różne sposoby.

## 3.1. Dostarcz tablicę wartości jako wejście do channel

Przejdźmy przez to na przykładzie. Powiedzmy, że chcemy wziąć te stringi wejściowe, ale zamiast po prostu umieszczać je bezpośrednio w fabryce channel, chcemy zdefiniować je jako tablicę.

## 3.1.1. Ustaw zmienną wejściową

Więc wezmę je i zrobię to jako nową linię powyżej i powiem, greetings, array.

Proszę bardzo. Wezmę tę zmienną array i umieszczę ją w channel.of, i nacisnę save.

## 3.1.3. Uruchom workflow

Teraz zobaczmy, co się stanie. Wracam do mojego terminala. Po prostu zamierzam uporządkować wszystkie te pliki tymczasowe ponownie. I uruchommy workflow.

Niedobrze. Dobra. To się zepsuło. W porządku. Spodziewałem się, że to się zepsuje tym razem. Debugowanie tego, co idzie nie tak, gdy workflow Nextflow zawodzi, jest kluczową częścią bycia deweloperem Nextflow. To będzie się zdarzać często i ważne jest, aby zrozumieć, co mówi komunikat o błędzie i jak sobie z tym radzić.

Komunikaty o błędach Nextflow są w rzeczywistości dość strukturalne. Mówi nam, który process poszedł nie tak. Podaje nam komunikat o błędzie z powodu. Mówi, jakie było polecenie, które próbowało uruchomić w ramach tego konkretnego zadania, jaki był status wyjścia, jakie było wyjście i gdzie był katalog work tego zadania.

Zauważ, że mogę kliknąć to opcją w VS Code, a otwiera się to w pasku bocznym, więc mogę przejść tam bezpośrednio i wyświetlić wszystkie te ukryte pliki, o których mówiliśmy w poprzednim rozdziale, włącznie z plikiem .command.sh. Jak widać, jest to to samo co polecenia, które zostały wykonane tutaj.

Patrząc na ten plik, możemy poczuć, co mogło pójść nie tak tutaj zamiast uruchamiania pojedynczego zadania dla każdego elementu w tablicy, jak to było ostatnim razem, po prostu dostarczyło całą tablicę na raz jako string. Więc musimy rozpakować tę tablicę na indywidualne wartości, zanim przekażemy ją do channel. Wróćmy i zobaczmy, czy możemy to zrobić za pomocą operatora.

## 3.2. Użyj operatora do przekształcenia zawartości channel

W tym przypadku nie zamierzamy zmieniać tablicy przed przekazaniem jej do channel. Zamierzamy dostosować channel tak, aby zachowywał się w sposób, jakiego oczekujemy. Zamierzamy to zrobić, używając operatora flatten, może zrobić dot, zacznij pisać i zobaczymy, że rozszerzenie VS Code zaczyna sugerować wszystkie różne operatory, które mamy dostępne.

## 3.2.1. Dodaj operator flatten()

I zamierzam wybrać flatten. Zauważ, że białe znaki nie mają znaczenia w tym kontekście dla Nextflow. Więc możesz umieścić te operatory w nowej linii, jeśli chcesz. Więc mogę upuścić to tutaj i wciąć, żeby znajdowało się pod ".of" i zobaczysz, że ludzie często łańcuchują wiele operatorów w ten sposób na channel i wcięli to w ten sposób, aby było łatwiej to czytać.

Możesz również zobaczyć, tak jak wcześniej, mogę najechać na to i przeczytać, co robi operator flatten, a także podążyć za linkiem do dokumentacji, jeśli chcę.

Więc ten operator bierze ten channel, który ma w sobie pojedynczą tablicę i rozdziela wartości tablicy.

## 3.2.2. Dodaj view() aby sprawdzić zawartość channel

Możemy zajrzeć do channels za pomocą specjalnego operatora view, i zamierzam dodać kilka z nich tutaj. To jest trochę jak używanie instrukcji print w innych językach. Więc zamierzam zrobić dot view, a następnie zamierzam użyć tych kręconych nawiasów.

To nazywa się closure. To zasadniczo daje dodatkowy kod do operatora view, który wykona na każdym elemencie w channel. W tym przypadku zamierzam powiedzieć greeting before flatten. Greeting.

Definiuję tutaj zmienną, która jest tylko w zakresie tego closure. Więc ta zmienna jest używana tylko tutaj i mogłem nazwać ją, jak chciałem. To naprawdę nie ma znaczenia. Po prostu używam greeting, aby było łatwo czytać.

W niektórych pipeline Nextflow możesz zobaczyć, że ludzie używają specjalnej niejawnej zmiennej o nazwie "$it". Tak jak to. To jest specjalna zmienna w kodzie Nextflow, która jest skrótem, więc nie musisz robić małej definicji zmiennej. Jednak z czasem myślimy, że to nie jest bardzo jasne dla ludzi, którzy są nowi w Nextflow, i teraz zniechęcamy do używania "$it".

Więc zamierzam trzymać się poprzedniego zachowania greeting i używać tego w ten sposób, ponieważ jest to bardziej jawne i jaśniejsze, co się dzieje.

Następnie skopiuję tę linię i zrobię dokładnie to samo ponownie po argumentach flatten. Operator view jest trochę specjalny, ponieważ robi coś na elementach, ale także po prostu kontynuuje przekazywanie ich do następnego operatora, więc możemy połączyć go w środku łańcucha operacji w ten sposób, a on wydrukuje tam status i będzie kontynuował. Więc miejmy nadzieję, że to pokaże nam, jak wygląda channel przed i po operatorze flatten.

## 3.2.3. Uruchom workflow

Wypróbujmy to. Wyczyść. Wyczyść wszystko w przestrzeni roboczej. Uruchom pipeline ponownie.

Dobra, więc możemy zobaczyć, że uruchomił nasze pięć procesów. Ponownie, nie zawiesił się z błędem, więc to zdecydowanie dobrze. I teraz mamy before flatten i rzeczywiście mamy naszą tablicę i mamy after flatten, wydrukowane pięć razy, raz dla każdego elementu tablicy. To dokładnie to, na co liczyliśmy. Więc to naprawdę dobra wiadomość. I to pasuje dokładnie do tego, czego oczekiwalibyśmy od kodu.

Nie potrzebujemy już tych instrukcji debugowania, więc mogę je albo zakomentować, albo usunąć. Zamierzam je usunąć, żeby utrzymać mój kod ładny i czysty. Dobra, świetnie. Ten przykład działa teraz ładnie i możemy zacząć widzieć, jak channels mogą robić nieco bardziej skomplikowaną logikę.

## 4. Użyj operatora do parsowania wartości wejściowych z pliku CSV

Teraz spróbujemy to zrobić, używając pliku z serią wejść zamiast tego. To jest bardzo powszechny sposób pisania pipeline Nextflow przy użyciu arkusza próbek lub CSV z metadanymi.

## 4.1. Zmodyfikuj skrypt, aby oczekiwał pliku CSV jako źródła powitań

Jeśli przejdę do paska bocznego, możesz zobaczyć greetings.csv w repozytorium przykładowym, i to jest bardzo, bardzo prosty plik CSV, który po prostu zawiera trzy linie z trzema różnymi powitaniami. Zobaczmy, czy możemy użyć tego pliku CSV w naszym workflow.

Teraz zamierzam wrócić do używania params, jak robiliśmy to w rozdziale pierwszym, abyśmy mogli mieć wejście wiersza poleceń.

Zamierzam usunąć tę tablicę greetings.

## 4.1.1. Przełącz parametr wejściowy na plik CSV

Zamierzam ustawić params greeting na nazwę pliku, która jest greetings.csv, i zamierzam użyć tej specjalnej zmiennej do wygenerowania channel. Zamierzam umieścić to tam, a błędy znikają. Pamiętaj, że to ustawia tę zmienną domyślnie teraz. Więc jeśli uruchomię pipeline bez żadnych argumentów, użyje greetings.csv, ale mogłem zrobić --greeting, aby nadpisać tę zmienną, gdybym chciał.

## 4.1.2. Przełącz się na fabrykę channel zaprojektowaną do obsługi pliku

Dobra, przekazujemy teraz plik zamiast stringa lub tablicy stringów, więc prawdopodobnie potrzebujemy innej fabryki channel.

Pozbędziemy się "of", którego używaliśmy do tej pory, a zamiast tego użyjemy .fromPath. To robi dokładnie to, jak brzmi. Tworzy channel ze ścieżkami zamiast wartości, używając nazwy pliku string lub glob. Zamierzam również usunąć operator flatten, ponieważ już go nie potrzebujemy, teraz, gdy przekazujemy plik.

## 4.1.3. Uruchom workflow

Zamierzam nacisnąć save, otworzyć terminal, uruchomić workflow i zobaczyć, co się stanie.

Dobra. Znowu się zawiesiło. Nie martw się. Tego też się spodziewałem. Spójrzmy na komunikat o błędzie i zobaczmy, czy możemy dowiedzieć się, co idzie nie tak. Tutaj możemy zobaczyć wykonane polecenie, i trochę jak wcześniej, gdzie mieliśmy wydrukowaną całą tablicę. Teraz mamy ścieżkę pliku echo do polecenia, zamiast przechodzenia przez zawartość pliku.

## 4.2. Użyj operatora splitCsv() do parsowania pliku

Więc aby użyć zawartości pliku zamiast tego, potrzebujemy innego operatora. Operator, którego zamierzamy użyć dla tego, nazywa się splitCsv. Ma sens, ponieważ to jest plik CSV, który ładujemy.

## 4.2.1. Zastosuj splitCsv() do channel

Ok, więc splitCsv. Zamknij nawias. Nie potrzebujemy tutaj żadnych argumentów. I znowu zamierzam użyć kilku operatorów view, aby dać pewien wgląd w to, co się tutaj dzieje.

.view csv after splitCsv. Before split Csv.

## 4.2.2. Uruchom workflow ponownie

Dobra, spróbujmy to uruchomić i zobaczymy, co się stanie.

Dobra, tym razem mamy trochę więcej wyjścia, ale nadal się zawiodło. Możemy spojrzeć na instrukcje view, i tutaj możesz zobaczyć before split CSV, i mamy ścieżkę pliku, jak widzieliśmy w poprzednim komunikacie o błędzie. After split CSV, teraz mamy trzy wartości odpowiadające trzem linii w pliku CSV.

Jednak możesz zobaczyć, że każda z tych wartości jest otoczona nawiasami kwadratowymi. Więc każda z nich była tablicą sama w sobie, i to dało nam ten sam obszar, który mieliśmy wcześniej, gdzie próbuje echo tablicę zamiast tylko pojedynczego stringa.

Jeśli pomyślimy o pliku CSV, ma to trochę sensu. Zazwyczaj plik CSV będzie miał wiersze i kolumny, więc split CSV robi dwuwymiarową tablicę. Pierwszy wymiar tablicy to każdy wiersz, a następnie jest drugi wymiar, który jest każdą kolumną dla każdego wiersza.

Więc tutaj mamy tylko pojedynczą wartość w każdej linii, więc mamy pojedynczą kolumnę, więc mamy tablicę jednooelementową dla każdej linii pliku.

To w porządku. Po prostu potrzebujemy kolejnego operatora, aby zwinąć tę tablicę dla każdej linii parsowanego pliku CSV. Posprzątajmy to. Pozbądźmy się terminala i zobaczmy, co możemy zrobić.

## 4.3. Użyj operatora map() do wyodrębnienia powitań

Teraz moglibyśmy użyć operatora flatten ponownie, którego używaliśmy wcześniej. Widzieliśmy, jak może zwinąć tablicę w serię wartości, co bardzo dobrze by tutaj zadziałało. Ale zamierzam wykorzystać okazję, aby zademonstrować inny operator, który jest bardzo powszechny w workflows, zwany operatorem map.

## 4.3.1. Zastosuj map() do channel

Zamierzam zrobić dot map i zamierzam zrobić item item[0].

Jeśli piszesz wiele innych języków kodu, możesz być już zaznajomiony z operatorem map. Bierze iterowalny, taki jak tablica lub channel, i wykonuje jakąś operację na każdej wartości tego.

Tutaj mówimy, że powinniśmy zdefiniować zmienną o nazwie item w zakresie tego closure, a następnie chcemy zwrócić, tylko pierwszą wartość w tej tablicy. Więc item indeks zero.

To jest skuteczne spłaszczanie tablicy. Możesz zobaczyć, jak moglibyśmy rozszerzyć to, aby było bardziej złożone, chociaż: gdyby nasz plik CSV miał sześć kolumn, ale jesteśmy zainteresowani tylko czwartą kolumną, moglibyśmy uzyskać dostęp do konkretnego indeksu tutaj. Lub wykonać jakikolwiek inny rodzaj operacji na wartości przed przekazaniem jej do przetwarzania downstream.

Więc operator map jest niezwykle elastyczny i bardzo potężny do modyfikowania channels w locie. Wstawmy kolejną instrukcję view, aby zobaczyć, co robi w naszym wykonaniu. Może adjudicat tę linię i przenieść ją w dół. I after map.

## 4.3.2. Uruchom workflow jeszcze raz

Wywołajmy terminal i spróbujmy uruchomić workflow.

Dobra, tym razem nie ma błędów. To dobry znak. Możemy teraz przejść przez wszystkie te różne wyjścia z instrukcji view. Before split CSV, mieliśmy pojedynczą ścieżkę. After split CSV, mieliśmy tablice jednowartościowe, a następnie after map, mamy tylko wartości bez żadnej składni tablicy. Przejdźmy do katalogu results, i tutaj są nasze pliki wyjściowe zachowujące się dokładnie tak, jak chcieliśmy.

Jest mały bonus tutaj. Możesz faktycznie zobaczyć, że operatory view są nieco pomieszane w kolejności, w jakiej wykonały wyjście. To dlatego, że Nextflow wykonuje równoległość tych różnych zadań. Więc po podzieleniu CSV, są trzy elementy w tym channel, i obsługuje przetwarzanie tych trzech elementów równolegle automatycznie. To oznacza, że kolejność wyjść jest stochastyczna i może się różnić. W tym przypadku po prostu zdarzyło się, że niektóre z operatorów view zwróciły po zakończeniu kolejnego kroku, więc przyszło w tej kolejności.

Jeśli uruchomię ten sam workflow ponownie. To rzeczywiście, przyszło w innej kolejności i tym razem mamy split CSV i mapy w kolejności, jakiej byśmy oczekiwali.

Więc po prostu pamiętaj, nie możesz polegać na kolejności wyjść z zadania procesu, ponieważ Nextflow obsługuje tę równoległość dla Ciebie automatycznie. Nextflow robi to dla Ciebie Swoją logiką przepływu danych, i to jest prawdziwa moc Nextflow.

Dobra, to prawdopodobnie jeden z najważniejszych rozdziałów całego szkolenia. Gdy zrozumiesz channels, fabryki channel i operatory, zaczniesz włączać się w siłę Nextflow i to, co czyni go wyjątkowym jako język programowania. Ta funkcjonalność pozwala Nextflow zrównoleglać wszystkie Twoje workflows dla Ciebie i generować niezwykle złożoną logikę workflow z bardzo czystą składnią i modelem przepływu danych push. To może być na początku trochę dziwna koncepcja, ale gdy już przyzwyczaisz się do pisania kodu w ten sposób, szybko poczuje się to naturalne i zanim się zorientujesz, będziesz pisać fantastyczne workflows.

Zrób sobie przerwę, filiżankę herbaty, spacer dookoła i przejdźmy do rozdziału trzeciego, gdzie zaczynamy rozszerzać te koncepcje na bardziej złożone workflows. Do zobaczenia w następnym filmie.

[Następna transkrypcja wideo :octicons-arrow-right-24:](03_hello_workflow.md)
