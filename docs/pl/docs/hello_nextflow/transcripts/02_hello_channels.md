# Część 2: Hello Channels - Transkrypcja wideo

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Tłumaczenie wspomagane przez AI - [dowiedz się więcej i zasugeruj ulepszenia](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

<div class="video-wrapper">
  <iframe width="560" height="315" src="https://www.youtube.com/embed/yDR66fzAMOg?si=xCItHLiOQWqoqBB9&amp;list=PLPZ8WHdZGxmWKozQuzr27jyMGqp9kElVK&amp;cc_load_policy=1" title="YouTube video player" frameborder="0" allow="accelerometer; autoplay; clipboard-write; encrypted-media; gyroscope; picture-in-picture; web-share" referrerpolicy="strict-origin-when-cross-origin" allowfullscreen></iframe>
</div>

!!!note "Ważne uwagi"

    Ta strona zawiera tylko transkrypcję. Pełne instrukcje krok po kroku znajdziesz w [materiałach szkoleniowych](../02_hello_channels.md).

    Numery sekcji pokazane w transkrypcji mają charakter orientacyjny i mogą nie obejmować wszystkich numerów sekcji w materiałach.

## Powitanie

Witaj ponownie w Części 2 szkolenia Hello Nextflow. Ten rozdział nosi tytuł Hello Channels.

Kanały są jak klej w Twoim pipeline'ie Nextflow. To elementy, które łączą wszystkie różne procesy, których Nextflow używa do przekazywania informacji i orkiestracji Twojego workflow'a.

Kanały mają jeszcze jeden aspekt - operatory. To w zasadzie funkcje, których możemy używać na kanałach, aby modyfikować ich zawartość. Przejdźmy do VS Code i zobaczmy, gdzie jesteśmy.

Jestem bardzo przybliżony w tym VS Code, więc żeby zachować porządek, usunąłem wszystkie pliki _.nextflow\*_ oraz katalog _work/_, _results/_ i wszystko z Rozdziału Pierwszego. Po prostu zaczynam od czystej karty. Ale nie martw się zbytnio. Jeśli nie chcesz, możesz zostawić te pliki. Nie spowodują żadnych problemów.

Zaczniemy od pracy nad _hello-channels.nf_ w tym rozdziale, a jeśli go otworzę, powinien wyglądać bardzo podobnie do pliku, nad którym pracowaliśmy wcześniej. Może się zdarzyć, że różne części są w innych miejscach skryptu, ale wszystko powinno być zasadniczo takie samo.

Jedna rzecz, która jest inna, to ścieżka w bloku `output` - teraz to _hello_channels_ dla tej części, co oznacza, że pliki wynikowe będą przechowywane w innym podkatalogu w Twoich wynikach, jeśli nadal tam je masz. Powinno to być ładne i czyste miejsce do rozpoczęcia bez zamieszania z plikami wyjściowymi.

Dobrze, więc szybko przypomnijmy sobie, co robi ten skrypt, gdy uruchamiamy ten workflow. Wykonujemy _"nextflow run hello-channels.nf"_. Możemy dodać _"--input myinput"_, a gdy to uruchomimy, użyje tego parametru, `params.input`, który został przekazany jako zmienna do procesu `sayHello` tutaj na górze, który trafia do `greeting` i zostaje zapisany do `output.txt`. I możemy to zobaczyć w pliku wyników. Świetnie.

## 1. Dostarczanie zmiennych wejściowych przez kanał w sposób jawny

To miłe. Ale jest dość uproszczone. Mamy jedną zmienną w tym parametrze, która trafia do procesu, który uruchamia się raz i nie skaluje się naprawdę. Nie możemy podać mu wielu różnych plików do utworzenia. Nie możemy podać mu wielu różnych powitań. Mamy tylko jedno.

W rzeczywistości Nextflow polega na skalowaniu Twojej analizy. Więc prawdopodobnie chcesz, żeby robił więcej niż jedną rzecz. I robimy to za pomocą _kanałów_.

Kanały to dość unikalny koncept dla wielu osób zaczynających pracę z Nextflow. Pochodzi z koncepcji programowania funkcyjnego i może zająć trochę czasu, zanim to zrozumiesz, ale gdy już kliknie, naprawdę odblokowują moc Nextflow i są kluczowe dla sposobu pisania Twoich workflow'ów.

## 1.1. Tworzenie kanału wejściowego

Zacznijmy od wzięcia tego skryptu i sprawienia, by używał _kanału_ zamiast tylko _parametru_.

Przechodzimy do workflow'a, który zawiera całą naszą logikę workflow'a dotyczącą łączenia rzeczy. I zamierzam tu wejść i utworzyć nowy kanał.

Tworzę nowy kanał.

I nazwę go "_greeting_ch"_. To konwencja, żeby dodawać "_\_ch"_ w ten sposób, po prostu żebyś mógł pamiętać, że ta zmienna jest kanałem. Ale możesz nazwać to, jak chcesz.

A potem powiem `equals`, i zrobię _"channel.of"._

`Channel` to jak przestrzeń nazw dla wszystkiego, co dotyczy kanałów. Małe "c", jeśli wcześniej używałeś Nextflow. A _".of"_ to coś, co nazywa się fabryką kanałów (Channel factory), co jest w zasadzie sposobem na utworzenie kanału.

Jest wiele różnych fabryk kanałów. Jeśli zrobię tylko ".", zobaczysz, że VS Code sugeruje ich mnóstwo, ale _".of"_ jest najprostszy i po prostu przyjmuje dane wejściowe tutaj.

Więc mogę zrobić nawiasy i powiem _"Hello Channels!"_.

Świetnie. Mam kanał. Fantastycznie. Mogę nacisnąć zapisz, mogę uruchomić ponownie, ale nic ciekawego się nie stanie. VS Code dał mi pomarańczową linię ostrzeżenia tutaj i powiedział, że to jest ustawione: utworzyłeś to, ale nigdy tego nie użyłeś do niczego. Ten kanał nie jest konsumowany.

Dobrze, więc jak go używamy? Bardzo prosto. Wezmę to, skopiuję i usunę _`params.input`_ i zamiast tego wstawię _"greeting_ch"_ tutaj. Więc przekażemy ten kanał jako wejście do `sayHello`.

Zauważ, że na razie zakodowałem ten ciąg znaków na stałe. To trochę krok wstecz po naszym ładnym parametrze, którego używaliśmy na końcu ostatniego rozdziału, ale po prostu upraszcza sprawy na początek, żebyś mógł zobaczyć logikę.

Dobrze, przejdę do mojego terminala i uruchomię workflow ponownie. Tym razem bez żadnego _"--input"_, i uruchomi się i użyje tego kanału, który utworzyliśmy i mam nadzieję, że powinniśmy mieć plik tutaj w _results/hello_channels/_ i teraz mówi "Hello Channels!". Fantastycznie. Więc to jest to, czego oczekiwaliśmy od naszego kanału tutaj. Świetnie.

## 1.4. Używanie `view()` do inspekcji zawartości kanału

Jeszcze jedna rzecz do dodania tutaj, tylko krótkie wprowadzenie do innej funkcji, której możemy używać na kanałach, zwanej "_.view"_.

Jest to analogiczne do polecenia _print_ w Pythonie lub innych językach, których możesz być przyzwyczajony, i po prostu zrzuca zawartość tego kanału do terminala, gdy go uruchamiamy.

Więc robię "_.view"_, a potem jeśli uruchomię workflow ponownie, powinien wydrukować do terminala, jaka jest zawartość tego kanału w momencie, gdy go utworzyliśmy.

I rzeczywiście, możesz zobaczyć, że wydrukował się do terminala tutaj. _"Hello Channels!"_.

Zauważ, że możesz łamać te rzeczy w wielu liniach, jeśli chcesz, i w rzeczywistości automatyczny formatter Nextflow będzie próbował to dla Ciebie zrobić. Białe znaki nie są tutaj naprawdę ważne, więc możesz łączyć te rzeczy jedną po drugiej.

## 2. Modyfikacja workflow'a do uruchamiania na wielu wartościach wejściowych

Dobrze, więc nasz kanał ma jedną rzecz, co jest miłe, ale to w zasadzie to samo co wcześniej. Więc uczyńmy to trochę bardziej skomplikowanym. Dodajmy kilka więcej rzeczy do naszego kanału.

Fabryka kanałów "_.of()"_ może przyjąć wiele elementów, więc napiszmy kilka więcej. Zrobimy _Hello, Bonjour, Hej_. A potem możemy uruchomić ten workflow ponownie i zobaczymy, co się stanie.

Powinien się uruchomić ponownie. I wydrukowaliśmy teraz. _"Hello", "Bonjour"_ i _"Hej"_ do terminala za pomocą naszej instrukcji `view`. Fantastycznie.

## 2.1.2. Uruchamianie polecenia i przeglądanie logów

Możesz myśleć, że w tym momencie skończyliśmy. Ale właściwie jest tu pewna pułapka, która nas złapie. Jeśli spojrzymy na nasz plik wyjściowy tutaj. Możesz zobaczyć, że ma _"Hello"_ w środku, ale nie ma żadnych innych wyników. W rzeczywistości to tylko ten jeden.

Jeśli uruchomimy ten workflow wiele razy, możemy nawet zobaczyć, że czasami ma _"Bonjour"_, czasami ma _"Hej"_. To trochę losowe.

Jeśli spojrzymy na terminal, możemy zobaczyć, że uruchomił się trzy razy i możemy zobaczyć różne wyjścia `view`. Ale jeśli przejdę do katalogu roboczego, mogę zrobić _"cat work"_. Wstawić ten hash i rozwinąć to i _output.txt_. Możesz zobaczyć, że ten plik w katalogu roboczym jest inny niż katalog wyników, a ten to _"Hej"._ Więc coś nie działa tutaj prawidłowo.

A kluczem jest to, że uruchomiły się trzy zadania. Wyjście Nextflow próbuje to podsumować w miarę postępu przetwarzania, tak aby nie przejęło całkowicie Twojego całego terminala, a to logowanie ANSI używa kodów ucieczki ANSI, które w zasadzie nadpisały inne zadania. Więc po prostu pokazuje Ci ostatnie, które akurat zostało zaktualizowane.

## 2.1.3. Ponowne uruchomienie polecenia z opcją `-ansi-log false`

Jest kilka rzeczy, które możemy zrobić, aby to lepiej zrozumieć. Możemy zajrzeć do samego katalogu roboczego i zobaczysz wszystkie różne katalogi robocze tam, ale to trochę mylące, bo będą pomieszane z różnymi uruchomieniami Nextflow.

Albo możemy powiedzieć Nextflow, żeby nie używał kodów ucieczki ANSI.

Więc jeśli uruchomię polecenie ponownie, ale tym razem powiem _"-ansi-log false"_, żeby to wyłączyć, mógłbym też użyć zmiennych środowiskowych _`$NO_COLOR`_ lub _"`$NXF_ANSI_LOG=false`"_. Wtedy używa bardziej staromodnego stylu logowania Nextflow bez żadnych z tych kodów ucieczki. Po prostu drukuje bezpośrednio do terminala bez żadnych sprytnych aktualizacji.

I teraz możemy zobaczyć wszystkie trzy z tych procesów, które się uruchomiły. I każdy z nich ma swój własny hash zadania. A jeśli wejdziemy do tych katalogów roboczych, zobaczymy trzy różne powitania, które określiliśmy.

Więc to ma teraz więcej sensu. Mam nadzieję, że rozumiesz, że Nextflow to robił, po prostu był trochę sprytny z tym, co Ci pokazywał w terminalu z tymi katalogami roboczymi.

Jednak to naprawiło jeden problem z katalogami roboczymi, ale nie naprawiło problemu z plikiem wyjściowym. Nadal mamy tylko jeden plik wyjściowy, który mówi _"Hello"_.

## 2.2. Zapewnienie unikalności nazw plików wyjściowych

Teraz, żeby to zrozumieć, musimy wrócić do naszego skryptu workflow'a. Generujemy nasz kanał tutaj, przekazujemy go do naszego procesu, a jeśli spojrzymy na proces, zapisujemy powitanie do pliku o nazwie _"output.txt"_ i przekazujemy ten plik wyjściowy z powrotem do bloku `output` tutaj na dole, publikując go.

Jednak za każdym razem, gdy ten proces uruchamia się trzy razy, te trzy różne zadania. Wszystkie generują plik o nazwie _"output.txt"_, wszystkie te pliki wyjściowe są publikowane do katalogu wyników i wszystkie się nawzajem nadpisują. Więc jakikolwiek plik wynikowy tam otrzymasz, to po prostu ostatni, który został wygenerowany, ale zniszczył wszystkie inne. To nie jest to, czego naprawdę chcemy.

## 2.2.1. Konstruowanie dynamicznej nazwy pliku wyjściowego

Są różne sposoby, aby sobie z tym poradzić, ale najprostszy na razie to po prostu utworzenie różnych unikalnych nazw plików. Więc za każdym razem, gdy zadanie uruchamia się z innym powitaniem, wygeneruje inny plik wyjściowy, który nie będzie już kolidował podczas publikowania. A wtedy otrzymamy trzy unikalne pliki wyjściowe.

Robimy to dokładnie w ten sam sposób. Możemy użyć tej zmiennej w dowolnym miejscu w bloku `script` i możemy jej użyć wiele razy.

Więc mogę wkleić to tutaj, _"`$\{greeting\}_output.txt`"_, a potem muszę też wkleić to tutaj na górze, ponieważ nie tworzymy już pliku o nazwie `output.txt`. Więc jeśli tego nie zaktualizuję, Nextflow się wysypie z błędem mówiącym, że oczekiwał pliku, który nigdy nie został wygenerowany.

Więc muszę zrobić to samo tam i muszę użyć podwójnych cudzysłowów, a nie pojedynczych, żeby ta zmienna była zrozumiana.

Dobrze, wypróbujmy to i zobaczmy, czy zadziałało. Uruchomimy workflow ponownie. Mam nadzieję, że pokaże nam trzy różne zadania w trzech różnych katalogach roboczych. I rzeczywiście, możesz zobaczyć w folderze wyników tutaj po lewej. Mamy teraz trzy różne pliki z trzema różnymi nazwami plików i każdy z różną zawartością, której oczekujemy. Więc pliki nie nadpisują się już nawzajem i wszystko jest tam, jak oczekujemy.

To trochę trywialna konfiguracja, przez którą przeszliśmy tutaj, ale podkreśla niektóre z kluczowych koncepcji, które musisz zrozumieć o tym, jak działa publikowanie plików, i niektóre z rzeczy, w które możesz wpaść jako pułapki. Więc mam nadzieję, że możesz tego uniknąć w swoich własnych workflow'ach.

Warto też zauważyć, że to, co tutaj zrobiliśmy, jest trochę niepraktyczne w rzeczywistych sytuacjach. Wzięliśmy jakieś dane wejściowe i używamy tych danych, ale także nazywamy plik po tych danych, czego zwykle nie można zrobić.

Więc w prawdziwych, bardziej dojrzałych pipeline'ach Nextflow często będziesz przekazywać obiekt meta ze wszystkimi metadanymi powiązanymi z daną próbką. Możesz wtedy tworzyć dynamiczne nazwy plików na podstawie tego, co jest o wiele bardziej praktyczne.

Jeśli jesteś zainteresowany, jak to zrobić zgodnie z najlepszymi praktykami, jest side quest na _training.nextflow.io_, który dotyczy specyficznie metadanych i map meta, więc możesz tam zagłębić się w więcej szczegółów.

## 3. Dostarczanie wielu danych wejściowych przez tablicę

Dobrze. Następnie zbadamy trochę, jak kanały są zbudowane i jak różnią się od innych rodzajów struktur danych w języku kodowania. I pomyślę trochę o tym, jak mógłbym potencjalnie użyć tablicy, która może być znajomą koncepcją, jeśli przyszedłeś z innych języków.

Czy mogę użyć tablicy w kanale? Spróbujmy. Utworzę tablicę, a skopiowałem to z dokumentacji, _"greetings_array"_ i _"Hello", "Bonjour"_ i _"Holà"_. A potem wstawię to tutaj zamiast moich zakodowanych na stałe ciągów znaków. Więc powiem `channel.of` _"greetings_array"_, przekazując tę tablicę do kanału. Spróbujmy.

Wywołuję terminal i uruchamiam pipeline.

Dobrze. Możesz zobaczyć, że instrukcja `view` tutaj wydrukowała naszą tablicę zgodnie z oczekiwaniami, ale potem cały ten czerwony tekst, lub nie będzie czerwony, jeśli nadal masz wyłączony _"-ansi-log"_, ale cały ten czerwony tekst mówi nam, że coś poszło nie tak.

Nie mamy już ładnego zielonego ptaszka. Mamy czerwony krzyżyk, a jeśli tylko zrobię to trochę szersze, żeby było łatwiej czytać, Nextflow mówi nam, co poszło nie tak.

Więc rozłóżmy to na sekcje. Mówi, że błąd został spowodowany przez, a potem powód błędu, którym są brakujące pliki wyjściowe. Więc w zasadzie ten blok `output` powiedział, że ten plik powinien zostać utworzony, a nie został. Następnie mówi, że to jest polecenie, które zostało wykonane. Więc to jest w zasadzie zawartość tego pliku _.command.sh_. Tak to wyglądało po tym, jak wszystkie te zmienne zostały wstawione.

I możesz zobaczyć tutaj, że nasze polecenie `echo` zostało uruchomione tylko raz i użyło całej tablicy, ale w reprezentacji ciągu znaków, co nie jest tym, czego chcieliśmy.

A potem polecenie zakończyło się w ten sposób, i to był katalog roboczy, gdzie możemy pójść i zobaczyć pliki, aby zrozumieć trochę więcej.

Dobrze. Więc co się stało, to Nextflow po prostu przekazał całą tę tablicę jako pojedynczy element kanału do procesu, co oznaczało, że proces uruchomił się tylko raz. Miał jedno zadanie i nie użył danych w strukturze, której oczekiwaliśmy.

## 3.2. Używanie operatora do transformacji zawartości kanału

Więc musimy coś zrobić z tym kanałem najpierw, zanim będzie mógł być użyty. I to przygotowuje scenę do używania operatorów, które są specjalnymi funkcjami, których możemy używać na kanałach do manipulowania zawartością kanałów.

W tym przypadku użyjemy czegoś, co nazywa się _flatten_. Które przekazujemy na końcu kanału tutaj. Więc tworzymy kanał, a potem uruchamiamy _flatten_. I znowu, jeśli najedziemy na to, pokazuje nam dokumentację dla tego polecenia od razu w VS Code, co jest bardzo pomocne. Możesz też znaleźć wszystkie te dokumenty na stronie Nextflow, w dokumentacji.

Mógłbym po prostu uruchomić ten kod teraz i zobaczyć, czy działa, ale to też dobra okazja, żeby wprowadzić, jak robić dynamiczny kod w operatorach i w kodzie Nextflow, które nazywają się domknięciami (closures).

Więc dodam z powrotem polecenie `view` tutaj, zanim uruchomimy _flatten_. A tutaj to ma te kręcone nawiasy, które są dynamicznym domknięciem. I jest po prostu jakiś arbitralny kod w środku, który zostanie wykonany w kontekście operatora `view`.

Tutaj to mówi, weź `greeting`, które jest wejściem operatora `view`, i to jest tutaj. Mógłbym nazwać to, jak chciałem, mógłbym nazwać to _"foo"_ i po prostu muszę się do tego odwoływać jako _"foo"_ później. A potem mówię z tym, zwróć to.

A potem ustawiam zwracanie ciągu znaków, który mówi `before flatten` dla zmiennej. Bardzo proste.

Teraz dodam kolejne dokładnie takie samo, ale powiem `after flatten`.

Więc to robi, ponieważ to działa sekwencyjnie, zobaczysz, jak wygląda kanał, zanim uruchomimy _flatten_, a potem znowu po uruchomieniu _flatten_.

A potem ten kanał `greeting` jest nadal tworzony, więc nadal będzie przekazany do procesu. I mam nadzieję, że teraz workflow się uruchomi. Wypróbujmy to.

Świetnie. Więc przede wszystkim pipeline tym razem się nie wysypał. Mieliśmy trzy procesy, które uruchomiły się prawidłowo i mamy mały ptaszek. A potem możemy zobaczyć, że nasze instrukcje `view` zadziałały.

Mamy `before flatten`, która jest tą tablicą, którą widzieliśmy wcześniej z awarii, a potem mamy trzy razy `after flatten` został wywołany, gdzie mamy _"Hello", "Bonjour"_ i wszystkie te trzy oddzielne elementy w tablicy, które są teraz, jak mieliśmy nadzieję, trzema oddzielnymi elementami w kanale.

I możesz zobaczyć, że operator `view` został uruchomiony trzy razy. I to dlatego, że ten kanał po _flatten_ ma teraz trzy elementy. I więc operator zostaje wywołany trzy razy.

Bardzo szybko, wspomniałbym tylko, że kiedy tworzyłem fabryki kanałów wcześniej, zrobiłem _"."_, a potem zobaczyliśmy, że było wiele różnych sposobów tworzenia kanałów, a jeden z nich nazywa się "_fromList"_. I to jest właściwie specjalnie zaprojektowane, żeby zrobić tę samą operację. Więc mogliśmy po prostu zrobić `fromList greetings_array`, i to zadziała. To nieco czystsza i ładniejsza składnia. Ale dla celów tej demonstracji chcieliśmy to zrobić bardziej krok po kroku, żebyś mógł zobaczyć, jak kanał jest manipulowany i jak różne operatory mogą zmieniać to, co jest w zawartości kanału.

## 4. Odczytywanie wartości wejściowych z pliku CSV

Dobrze, jak możemy to uczynić bardziej realistycznym? Prawdopodobnie nie będziesz chciał tworzyć mnóstwa kodu w swoim pipeline'ie Nextflow z zakodowanymi na stało tablicami. Prawdopodobnie będziesz chciał wziąć dane z zewnątrz, gdy uruchamiasz, a te dane prawie na pewno będą w plikach.

Więc następną rzeczą, którą zrobimy, jest to, że zreplikujemy to, ale zamiast brać dane z pojedynczego parametru CLI lub z zakodowanego na stało ciągu znaków lub tablicy, weźmiemy to z pliku.

Więc pozbądźmy się naszego `greetings_array`. A teraz zmienimy tę fabrykę kanałów ponownie. Właśnie powiedziałem, że jest kilka do wyboru i jest jedna zwana _".fromPath"_. I powiem jej, żeby w tym przypadku wzięła _`params.input`_, co wraca do naszego wejścia, którego używaliśmy wcześniej.

Teraz ten parametr nie jest jeszcze gotowy do użycia. Nadal mówimy, że to ciąg znaków i jest zakodowany na stałe tutaj z domyślną wartością, ale moglibyśmy nadpisać ten ciąg znaków. Teraz chcemy, żeby to był plik zamiast tego. Więc typ jest inny. To nie jest już _String_. To _Path_.

A potem możemy ustawić wartość domyślną, jeśli chcemy, znowu na `Path`. A jeśli spojrzę w eksplorator po lewej, możesz zobaczyć w tym repozytorium, w tym katalogu roboczym, mam katalog o nazwie `data`. Mam tam plik o nazwie _"greetings.csv"._

Więc mogę po prostu ustawić wartość domyślną tutaj na _"data/greetings.csv"_. Teraz, gdy uruchomię ten pipeline ponownie bez żadnych opcji wiersza poleceń, użyje tej wartości domyślnej. Wie, że to ścieżka, więc wie, że powinien to obsługiwać jako ścieżkę, a nie ciąg znaków.

A potem przekaże to do fabryki kanałów z tego _`params.input`_ i utworzy nasz kanał, który następnie będzie użyty w tym procesie zwanym _sayHello_. Wypróbujmy to.

Dobrze. Nie powiodło się. Nie martw się. To było oczekiwane. A jeśli śledzisz materiał szkoleniowy, zobaczysz, że tam też było oczekiwane. Zobaczmy, co się dzieje tutaj.

Próbował uruchomić pipeline. Próbował wykonać proces i dostał dość podobny błąd do tego, który widzieliśmy wcześniej.

Tutaj mówi: próbowaliśmy uruchomić `echo`, ale zamiast wyechować zawartość tego pliku CSV, po prostu wyechował ścieżkę. I możesz zobaczyć, że to pełna ścieżka bezwzględna tutaj do tego pliku CSV.

A potem rzeczywiście, ponieważ próbował zapisać to do tej naprawdę skomplikowanej ścieżki, nie bardzo wiedział, co zrobić. I było to poza zakresem katalogu roboczego procesu.

Wspomniałem na początku, że Nextflow enkapsuluje każde wykonane zadanie w specjalnym katalogu roboczym. A jeśli próbujesz zapisać dane, które są poza tym katalogiem roboczym, Nextflow Cię zatrzyma jako środek ostrożności. I to właśnie się tutaj stało. Próbowaliśmy zapisać do ścieżki bezwzględnej i Nextflow zawiódł i nas powstrzymał.

## 4.2. Używanie operatora `splitCsv()` do parsowania pliku

Dobrze, spójrzmy na ten kanał i zobaczmy, jak wygląda. Możemy zrobić _".view"_, a skopiowałem to ze strony. Więc `.view`, i mamy dynamiczne domknięcie tutaj i mówimy nazwa zmiennej "_csv"_ jako wejście. Więc to zawartość kanału, i mówimy `before splitCsv`, i tak to wygląda.

Jeśli uruchomię to ponownie, nadal się nie powiedzie, ale pokaże nam, co jest w tym kanale. To nie jest szczególnie ekscytujące. To ta zmienna _path_. Więc możesz zobaczyć, że to po prostu ciąg znaków tutaj, ponieważ jest drukowany do terminala, ale to jest obiekt _path_, który zawiera informacje i metadane o tym pliku.

Nie chcemy przekazywać metadanych pliku do wejścia. Chcemy przekazać zawartość tego pliku. Jeśli spojrzymy na plik _greetings.csv_, możesz zobaczyć tutaj, że ma te różne zmienne tutaj. _Hello, Bonjour, Holà_ znowu. I to są rzeczy, które naprawdę chcemy przekazywać do naszego procesu, a nie tylko sam plik jako pojedynczy obiekt.

Więc musimy sparsować ten plik CSV. Musimy go rozpakować, dostać się do zawartości pliku CSV, a potem przekazać zawartość w kanale do procesu.

Jak prawdopodobnie możesz wywnioskować z komunikatu logowania, chcemy użyć _splitCsv_, który jest kolejnym operatorem, kolejnym operatorem kanału. Więc jeśli zrobię "_dot" "s"_, a potem zobaczysz, że jest automatycznie sugerowany. Ups, _splitCsv_ i nawiasy.

A potem po _splitCsv_, wstawię kolejną instrukcję `view`, żebyśmy mogli zobaczyć, jak to wygląda potem. Uruchommy pipeline i zobaczmy, co mamy.

Dobrze. Nadal się nie powiodło, ale w nowy i ekscytujący sposób, co jest postępem.

Tym razem znowu mamy jakiś problem z naszym skryptem, który został wyrenderowany. Teraz. Nie mamy już końcowej ścieżki, ale mamy tablicę zmiennych, która wygląda bardzo podobnie do błędu, który mieliśmy wcześniej, gdy przekazywaliśmy tablicę jako stałe wejście.

Z naszym logowaniem z operatora `view` możemy zobaczyć, że `before splitCsv` była ścieżka. I rzeczywiście, po _splitCsv_ mamy trzy różne wyjścia i każde z tych wyjść wygląda strasznie podobnie do każdego z wierszy z pliku _greetings.csv_, co ma sens.

Więc co się tutaj stało, to Nextflow sparsował ten plik CSV, dał nam trzy obiekty, jedną tablicę dla każdej linii pliku CSV. Więc potem trzy razy przekazaliśmy tablicę zmiennych do kanału zamiast pojedynczej wartości ciągu znaków.

Dobrze, więc ostatnim razem mieliśmy ten problem, użyliśmy _flatten_. Spróbujmy bardzo szybko. Wypróbujmy `flatten` i zobaczmy, co się stanie.

Mogę nazwać te zmienne, jak chcę. Więc nazwę to _myarray_, bo to nie jest już naprawdę CSV. Spróbujmy uruchomić to ponownie i zobaczmy, co się stanie z _flatten_.

Więc tym razem uruchomimy, sparsowaliśmy CSV na trzy obiekty tablicowe, a potem je spłaszczyliśmy. I tym razem przeszło. I pipeline Nextflow się uruchomił. Jednak możesz zobaczyć, że _flatten_ naprawdę się rozpędza i spłaszcza wszystko. I więc otrzymujemy trzy niezależne wpisy tablicowe dla każdego wiersza. I więc uruchomił proces trzy razy dla każdego wiersza CSV. I teraz mamy całą masę plików wyników, i 123, 456, i wszelkiego rodzaju rzeczy, nie tylko tę pierwszą kolumnę CSV, której naprawdę chcieliśmy.

## 4.3. Używanie operatora `map()` do wyodrębnienia powitań

Więc jak dostać się tylko do pierwszej kolumny? Jeśli `flatten` jest tutaj zbyt uproszczony, potrzebujemy bardziej złożonego operatora, gdzie możemy faktycznie dostosować i powiedzieć mu, czego chcemy z CSV.

Aby to zrobić, użyjemy _map_. W zasadzie _map_ po prostu mówi, uruchom jakiś kod, jakąś funkcję na każdym elemencie, który otrzymam i zrób jakąś transformację na nim. A ponieważ jest tak elastyczny, zobaczysz, że pojawia się w kodzie Nextflow cały czas.

Sam w sobie nic nie robi. Więc nie chcemy zwykłych nawiasów, chcemy tutaj domknięcie i musimy mu powiedzieć, co zrobić. Więc powiem _"row"_, bo dostaje wiersze z CSV, więc to logiczna nazwa zmiennej. Jest wejściem. I chcę zwrócić tylko pierwszy element tej tablicy.

Tablice w Nextflow są indeksowane od zera, więc powiemy tylko pierwszy element, który jest wierszem zero. Gdybym chciał drugą kolumnę, mógłbym być jeden lub trzecią kolumnę być dwa, i tak dalej. Możemy zwrócić tutaj, co chcemy, ale zwrócę tylko pierwszą wartość.

I teraz możemy uruchomić pipeline ponownie i zobaczyć, czy robi to, czego oczekujemy.

I rzeczywiście, po _splitCsv_ mamy nasze tablice, a potem po _map_ mamy nasze ładne czyste ciągi znaków, tylko _"Hello", "Bonjour"_ i _"Holà"_. I pipeline teraz robi to, czego chcemy. Fantastycznie.

Więc możemy teraz pozbyć się wszystkich tych poleceń `view`. Nie potrzebujemy ich już.

## Podsumowanie

Skończyliśmy nasze debugowanie i to jest kod, z którym kończymy. Bierzemy nasz parametr CLI zwany _input_, który jest sklasyfikowany jako _Path_. Nextflow znajduje ścieżkę, ładuje ją i rozumie plik CSV. Zwraca wszystkie różne wiersze. A potem mapujemy tylko pierwszy element tego wiersza do kanału, który daje nam zawartość kanału, która jest przekazywana do procesu.

A proces uruchamia się na każdym elemencie w kanale, których jest trzy. I uruchamia proces trzy razy, dając mu trzy zadania. A te wyniki są następnie publikowane z workflow'a, przechwycone przez wyjście procesu. Publikowane z workflow'a i zapisane w bloku `output` do podkatalogu o nazwie _"hello_channels"_.

Całkiem fajnie. Zbliżamy się teraz do czegoś, co bardziej przypomina prawdziwy pipeline Nextflow, który mógłbyś uruchomić dla jakiejś prawdziwej analizy.

## Podsumowanie

Dobrze. Mam nadzieję, że teraz masz wyczucie, czym są kanały i operatory Nextflow i jak operatory działają na kanałach i jak możesz je tworzyć.

Kanały, jak powiedziałem na początku tego wideo, są klejem Nextflow. I możesz zobaczyć tutaj, że możemy wziąć różne wejścia i manipulować nimi i wziąć te dane, a potem przekazać je do logiki workflow'a niżej.

I ten blok `workflow` tutaj jest naprawdę miejscem, gdzie budujesz całą tę paralelizację i całą sprytną logikę, i wyjaśniasz Nextflow, jak zbudować Twój DAG workflow'a i jak orkiestrować Twój pipeline.

Kanały nie są najłatwiejszą koncepcją do zrozumienia. Więc zrób sobie przerwę, pomyśl trochę o tym, może przeczytaj materiał ponownie i naprawdę upewnij się, że masz te koncepcje opanowane, ponieważ to jest kluczowe dla Twojego zrozumienia Nextflow i im lepiej rozumiesz kanały i różne operatory kanałów i różne fabryki kanałów. Tym więcej zabawy będziesz miał pisząc Nextflow i tym potężniejsze będą Twoje pipeline'y.

To nie jest to samo co zwykłe programowanie w Pythonie lub innych językach. Nie używamy tutaj instrukcji _if_, to jest funkcyjne programowanie przepływowe używające kanałów i operatorów. Więc to trochę inne, ale jest też super potężne.

To koniec tego rozdziału. Idź i zrób sobie szybką przerwę, a zobaczę Cię w następnym wideo dla części trzeciej, gdzie przejdziemy przez Hello Workflow i porozmawiamy trochę więcej o workflow'ach.

Tak jak w poprzednim rozdziale, jest kilka pytań quizowych na dole tej strony internetowej tutaj, więc możesz szybko przez nie przejść i upewnić się, że rozumiesz wszystkie różne części materiału, który właśnie zrobiliśmy. A poza tym zobaczę Cię w następnym wideo. Dziękuję bardzo.

Dobrze.
