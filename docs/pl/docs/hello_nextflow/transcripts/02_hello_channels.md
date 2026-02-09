# Część 2: Hello Channels - Transkrypcja wideo

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Tłumaczenie wspomagane przez AI - [dowiedz się więcej i zaproponuj poprawki](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

<div class="video-wrapper">
  <iframe width="560" height="315" src="https://www.youtube.com/embed/yDR66fzAMOg?si=xCItHLiOQWqoqBB9&amp;list=PLPZ8WHdZGxmWKozQuzr27jyMGqp9kElVK&amp;cc_load_policy=1" title="YouTube video player" frameborder="0" allow="accelerometer; autoplay; clipboard-write; encrypted-media; gyroscope; picture-in-picture; web-share" referrerpolicy="strict-origin-when-cross-origin" allowfullscreen></iframe>
</div>

!!! note "Ważne uwagi"

    Ta strona zawiera tylko transkrypcję. Pełne instrukcje krok po kroku znajdziesz w [materiałach szkoleniowych](../02_hello_channels.md).

    Numery sekcji pokazane w transkrypcji mają charakter orientacyjny i mogą nie obejmować wszystkich numerów sekcji w materiałach.

## Powitanie

Witaj ponownie w części 2 szkolenia Hello Nextflow. Ten rozdział nosi tytuł Hello Channels.

Kanały są jak klej w Twoim pipeline'ie Nextflow'a. To elementy, które łączą wszystkie różne procesy, których Nextflow używa do przekazywania informacji i orkiestracji workflow'u.

Kanały mają jeszcze jeden aspekt - operatory. To w zasadzie funkcje, których możemy używać na kanałach, aby modyfikować ich zawartość. Przejdźmy do VS Code i zobaczmy, gdzie jesteśmy.

Bardzo mocno przybliżyłem widok w VS Code, więc żeby zachować porządek, usunąłem wszystkie pliki _.nextflow\*_ oraz katalog _work/_, _results/_ i wszystko z rozdziału pierwszego. Po prostu zaczynam od czystej karty. Ale nie martw się zbytnio - jeśli nie chcesz, możesz zostawić te pliki. Nie spowodują żadnych problemów.

Zaczniemy od pracy nad plikiem _hello-channels.nf_ w tym rozdziale. Gdy go otworzę, powinien wyglądać bardzo podobnie do pliku, nad którym pracowaliśmy wcześniej. Może się zdarzyć, że różne części są w innych miejscach skryptu, ale wszystko powinno być zasadniczo takie samo.

Jedna rzecz jest inna - ścieżka w bloku output to teraz _hello_channels_ dla tej części, co oznacza, że pliki wynikowe będą przechowywane w innym podkatalogu w Twoim katalogu results, jeśli nadal go tam masz. Powinno to być czyste miejsce do rozpoczęcia bez zamieszania z plikami wyjściowymi.

Dobrze, szybko przypomnijmy sobie, co robi ten skrypt. Gdy uruchamiamy ten workflow, wykonujemy _"nextflow run hello-channels.nf"_. Możemy dodać _"--input myinput"_, a gdy to uruchomimy, użyje tego parametru, params.input, który zostanie przekazany jako zmienna do procesu sayHello tutaj na górze, trafi do greeting i zostanie zapisany do output.txt. Możemy to zobaczyć w pliku wynikowym. Świetnie.

## 1. Dostarczanie zmiennych wejściowych przez kanał w sposób jawny

To miłe. Ale jest dość uproszczone. Mamy jedną zmienną w tym parametrze, która trafia do procesu uruchamianego raz i nie skaluje się zbyt dobrze. Nie możemy podać wielu różnych plików do utworzenia. Nie możemy podać wielu różnych powitań. Mamy tylko jedno.

W rzeczywistości Nextflow służy do skalowania Twojej analizy. Prawdopodobnie chcesz, żeby robił więcej niż jedną rzecz. I robimy to za pomocą _kanałów_.

Kanały to dość unikalny koncept dla wielu osób zaczynających pracę z Nextflow'em. Pochodzi z koncepcji programowania funkcyjnego i może zająć trochę czasu, zanim to zrozumiesz, ale gdy już kliknie, naprawdę odblokują moc Nextflow'a i są kluczowe dla pisania workflow'ów.

## 1.1. Tworzenie kanału wejściowego

Zacznijmy od wzięcia tego skryptu i sprawienia, by używał _kanału_ zamiast tylko _parametru_.

Przechodzimy do workflow'u, gdzie znajduje się cała nasza logika łączenia rzeczy. Wejdę tutaj i utworzę nowy kanał.

Utworzę nowy kanał.

I nazwę go "_greeting_ch"_. Konwencją jest dodawanie "_\_ch"_ w ten sposób, żebyś pamiętał, że ta zmienna jest kanałem. Ale możesz nazwać to, jak chcesz.

Następnie napiszę equals i zrobię _"channel.of"._

Channel to przestrzeń nazw dla wszystkiego, co dotyczy kanałów. Małe "c", jeśli wcześniej używałeś Nextflow'a. A _".of"_ to coś, co nazywa się fabryką kanałów (Channel factory), co jest w zasadzie sposobem na utworzenie kanału.

Istnieje wiele różnych fabryk kanałów. Jeśli napiszę tylko ".", zobaczysz, że VS Code sugeruje ich mnóstwo, ale _".of"_ jest najprostszy i po prostu przyjmuje dane wejściowe.

Mogę więc dodać nawiasy i napiszę _"Hello Channels!"_.

Świetnie. Mam kanał. Fantastycznie. Mogę zapisać, mogę uruchomić ponownie, ale nic ciekawego się nie stanie. VS Code dał mi pomarańczową linię ostrzeżenia tutaj i powiedział, że to jest ustawione: utworzyłeś to, ale nigdy tego nie użyłeś do niczego. Ten kanał nie jest konsumowany.

Dobrze, więc jak go użyć? Bardzo prosto. Skopiuję to i usunę _params.input_, a zamiast tego wstawię _"greeting_ch"_. Przekażemy więc ten kanał jako wejście do sayHello.

Zauważ, że na razie zakodowałem ten string na sztywno. To trochę krok wstecz po naszym ładnym parametrze, którego używaliśmy na końcu ostatniego rozdziału, ale po prostu upraszcza sprawy na początek, żebyś mógł zobaczyć logikę.

Dobrze, przejdę do terminala i uruchomię workflow ponownie. Tym razem bez żadnego _"--input"_, uruchomi się i użyje kanału, który utworzyliśmy. Mam nadzieję, że powinniśmy mieć plik tutaj w _results/hello_channels/_ i teraz mówi "Hello Channels!". Fantastycznie. To jest to, czego oczekiwaliśmy od naszego kanału. Świetnie.

## 1.4. Używanie view() do inspekcji zawartości kanału

Jeszcze jedna rzecz do dodania - szybkie wprowadzenie do innej funkcji, której możemy używać na kanałach, zwanej "_.view"_.

Jest to analogiczne do polecenia _print_ w Pythonie lub innych językach, których możesz używać, i po prostu wyrzuca zawartość tego kanału do terminala, gdy go uruchamiamy.

Więc napiszę "_.view"_, a następnie jeśli uruchomię workflow ponownie, powinien wypisać do terminala zawartość tego kanału w momencie, gdy go utworzyliśmy.

I rzeczywiście, możesz zobaczyć, że wypisało się do terminala. _"Hello Channels!"_.

Zauważ, że możesz łamać te rzeczy na wiele linii, jeśli chcesz, i w rzeczywistości automatyczny formatter Nextflow'a będzie próbował to dla Ciebie robić. Białe znaki nie są tutaj naprawdę ważne, więc możesz łączyć te rzeczy jedną po drugiej.

## 2. Modyfikacja workflow'u do pracy z wieloma wartościami wejściowymi

Dobrze, więc nasz kanał ma jedną rzecz, co jest miłe, ale to w zasadzie to samo co wcześniej. Więc skomplikujmy to trochę. Dodajmy kilka więcej rzeczy do naszego kanału.

Fabryka kanałów "_.of()"_ może przyjąć wiele elementów, więc napiszmy kilka więcej. Zrobimy _Hello, Bonjour, Hej_. A potem możemy uruchomić ten workflow ponownie i zobaczymy, co się stanie.

Powinien się uruchomić ponownie. I teraz wypisaliśmy _"Hello", "Bonjour"_ i _"Hej"_ do terminala za pomocą naszej instrukcji view. Fantastycznie.

## 2.1.2. Uruchomienie polecenia i sprawdzenie logów

Możesz myśleć, że w tym momencie skończyliśmy. Ale właściwie jest tutaj pewna pułapka, która nas złapie. Jeśli spojrzymy na nasz plik wyjściowy tutaj, możesz zobaczyć, że ma _"Hello"_, ale nie ma żadnych innych wyników. W rzeczywistości to tylko ten jeden.

Jeśli uruchomimy ten workflow wiele razy, możemy nawet zobaczyć, że czasami ma _"Bonjour"_, czasami _"Hej"_. To trochę losowe.

Jeśli spojrzymy na terminal, widzimy, że uruchomił się trzy razy i widzimy różne wyjścia view. Ale jeśli przejdę do katalogu work, mogę zrobić _"cat work"_. Wstawię ten hash, rozwinę to i _output.txt_. Możesz zobaczyć, że ten plik w katalogu work jest inny niż w katalogu results, a ten to _"Hej"_. Więc coś nie działa tutaj prawidłowo.

Kluczem jest to, że uruchomiły się trzy zadania. Wyjście Nextflow'a próbuje to podsumować w miarę postępu przetwarzania, żeby nie przejęło całkowicie Twojego terminala, a to logowanie ANSI używa kodów escape ANSI i w zasadzie nadpisało inne zadania. Więc pokazuje Ci tylko ostatnie, które akurat zostało zaktualizowane.

## 2.1.3. Ponowne uruchomienie polecenia z opcją -ansi-log false

Jest kilka rzeczy, które możemy zrobić, żeby to lepiej zrozumieć. Możemy zajrzeć do samego katalogu work i zobaczyć wszystkie różne katalogi work tam, ale to trochę mylące, bo będą pomieszane z różnymi uruchomieniami Nextflow'a.

Albo możemy powiedzieć Nextflow'owi, żeby nie używał kodów escape ANSI.

Więc jeśli uruchomię polecenie ponownie, ale tym razem powiem _"-ansi-log false"_, żeby to wyłączyć, mógłbym też użyć zmiennych środowiskowych _$NO_COLOR_ lub _"$NXF_ANSI_LOG=false"_. Wtedy używa bardziej staromodnego stylu logowania Nextflow'a bez żadnych z tych kodów escape. Po prostu wypisuje bezpośrednio do terminala bez żadnych sprytnych aktualizacji.

I teraz możemy zobaczyć wszystkie trzy te procesy, które się uruchomiły. I każdy z nich ma swój własny hash zadania. A jeśli wejdziemy do tych katalogów work, zobaczymy trzy różne powitania, które określiliśmy.

Więc to ma teraz więcej sensu. Mam nadzieję, że rozumiesz, że Nextflow to robił, po prostu był trochę sprytny z tym, co Ci pokazywał w terminalu z tymi katalogami work.

Jednak to naprawiło jeden problem z katalogami work, ale nie naprawiło problemu z plikiem wyjściowym. Nadal mamy tylko jeden plik wyjściowy, który mówi _"Hello"_.

## 2.2. Zapewnienie unikalności nazw plików wyjściowych

Teraz, żeby to zrozumieć, musimy wrócić do naszego skryptu workflow'u. Generujemy nasz kanał tutaj, przekazujemy go do naszego procesu, a jeśli spojrzymy na proces, zapisujemy greeting do pliku o nazwie _"output.txt"_ i przekazujemy ten plik wyjściowy z powrotem do bloku output tutaj na dole, publikując go.

Jednak za każdym razem, gdy ten proces uruchamia się trzy razy, te trzy różne zadania wszystkie generują plik o nazwie _"output.txt"_, wszystkie te pliki wyjściowe są publikowane do katalogu results i wszystkie się nawzajem nadpisują. Więc jakikolwiek plik wynikowy tam otrzymasz, to po prostu ostatni, który został wygenerowany, ale zniszczył wszystkie pozostałe. To nie jest to, czego chcemy.

## 2.2.1. Konstruowanie dynamicznej nazwy pliku wyjściowego

Są różne sposoby, żeby sobie z tym poradzić, ale najprostszy na razie to po prostu utworzenie różnych unikalnych nazw plików. Więc za każdym razem, gdy zadanie uruchamia się z innym powitaniem, wygeneruje inny plik wyjściowy, który nie będzie już kolidował podczas publikowania. I wtedy otrzymamy trzy unikalne pliki wyjściowe.

Robimy to dokładnie w ten sam sposób. Możemy użyć tej zmiennej gdziekolwiek w bloku script i możemy jej użyć wiele razy.

Mogę więc wkleić to tutaj, _"$\{greeting\}\_output.txt"_, a następnie muszę też wkleić to tutaj na górze, ponieważ nie tworzymy już pliku o nazwie _output.txt_. Więc jeśli tego nie zaktualizuję, Nextflow zawali się z błędem mówiącym, że oczekiwał pliku, który nigdy nie został wygenerowany.

Więc muszę zrobić to samo tam i muszę użyć podwójnych cudzysłowów, a nie pojedynczych, żeby ta zmienna była rozumiana.

Dobrze, wypróbujmy to i zobaczmy, czy zadziałało. Uruchomimy workflow ponownie. Mam nadzieję, że pokaże nam trzy różne zadania w trzech różnych katalogach work. I rzeczywiście, możesz zobaczyć w folderze results tutaj po lewej. Mamy teraz trzy różne pliki z trzema różnymi nazwami plików i każdy z różną zawartością, której oczekujemy. Więc pliki nie nadpisują się już nawzajem i wszystko jest tam, jak oczekujemy.

To trochę trywialna konfiguracja, przez którą przeszliśmy, ale podkreśla niektóre kluczowe koncepcje, które musisz zrozumieć o tym, jak działa publikowanie plików, i niektóre pułapki, w które możesz wpaść. Więc mam nadzieję, że możesz tego uniknąć we własnych workflow'ach.

Warto też zauważyć, że to, co zrobiliśmy tutaj, jest trochę niepraktyczne w rzeczywistych sytuacjach. Wzięliśmy jakieś dane wejściowe i używamy tych danych, ale także nazywamy plik po tych danych, czego zwykle nie można robić.

Więc w prawdziwych, bardziej dojrzałych pipeline'ach Nextflow'a często będziesz przekazywał obiekt meta ze wszystkimi metadanymi powiązanymi z daną próbką. Możesz wtedy tworzyć dynamiczne nazwy plików na podstawie tego, co jest o wiele bardziej praktyczne.

Jeśli jesteś zainteresowany, jak to robić zgodnie z najlepszymi praktykami, jest side quest na _training.nextflow.io_, który dotyczy specyficznie metadanych i map meta, więc możesz tam zagłębić się w szczegóły.

## 3. Dostarczanie wielu danych wejściowych przez tablicę

Dobrze. Teraz zbadamy trochę, jak kanały są zbudowane i czym różnią się od innych rodzajów struktur danych w języku programowania. I pomyślę trochę o tym, jak mógłbym potencjalnie użyć tablicy, która może być znajomym konceptem, jeśli przyszedłeś z innych języków.

Czy mogę użyć tablicy w kanale? Spróbujmy. Utworzę tablicę, a skopiowałem to z dokumentacji, _"greetings_array"_ i _"Hello", "Bonjour"_ i _"Holà"_. A potem wstawię to tutaj zamiast moich zakodowanych na sztywno stringów. Więc powiem "channel.of" _"greetings_array"_, przekazując tę tablicę do kanału. Spróbujmy.

Wywołam terminal i uruchomię pipeline.

Dobrze. Możesz zobaczyć, że instrukcja view tutaj wypisała naszą tablicę zgodnie z oczekiwaniami, ale potem cały ten czerwony tekst, lub nie będzie czerwony, jeśli nadal masz wyłączone _"-ansi-log"_, ale cały ten czerwony tekst mówi nam, że coś poszło nie tak.

Nie mamy już ładnego zielonego ptaszka. Mamy czerwony krzyżyk, a jeśli tylko trochę to poszerzę, żeby było łatwiej czytać, Nextflow mówi nam, co poszło nie tak.

Więc rozłóżmy to na sekcje. Mówi, że błąd został spowodowany przez, a następnie powód błędu, którym są brakujące pliki wyjściowe. Więc w zasadzie ten blok output powiedział, że ten plik powinien zostać utworzony, a nie został. Następnie mówi, że to jest polecenie, które zostało wykonane. Więc to jest w zasadzie zawartość tego pliku _.command.sh_. Tak to wyglądało po tym, jak wszystkie te zmienne zostały wstawione.

I możesz zobaczyć tutaj, że nasze polecenie echo faktycznie zostało uruchomione tylko raz i użyło całej tablicy, ale w reprezentacji stringowej, co nie jest tym, czego chcieliśmy.

A potem polecenie zakończyło się w ten sposób, i to był katalog work, gdzie możemy pójść i zobaczyć pliki, żeby lepiej to zrozumieć.

Dobrze. Więc co się stało, to Nextflow po prostu przekazał całą tę tablicę jako pojedynczy element kanału do procesu, co oznaczało, że proces uruchomił się tylko raz. Miał jedno zadanie i nie użył danych w strukturze, której oczekiwaliśmy.

## 3.2. Używanie operatora do transformacji zawartości kanału

Więc musimy coś zrobić z tym kanałem najpierw, zanim będzie mógł być użyty. I to przygotowuje scenę do używania operatorów, które są specjalnymi funkcjami, których możemy używać na kanałach do manipulowania zawartością kanałów.

W tym przypadku użyjemy czegoś, co nazywa się _flatten_. Które przekazujemy na końcu kanału tutaj. Więc tworzymy kanał, a następnie uruchamiamy _flatten_. I znowu, jeśli najedziemy na to, pokazuje nam dokumentację dla tego polecenia od razu w VS Code, co jest bardzo pomocne. Możesz też znaleźć wszystkie te dokumenty na stronie Nextflow'a, w dokumentacji.

Mógłbym po prostu uruchomić ten kod teraz i zobaczyć, czy działa, ale to też dobra okazja, żeby wprowadzić, jak robić dynamiczny kod w operatorach i w kodzie Nextflow'a, które nazywają się closures.

Więc dodam z powrotem polecenie view tutaj, zanim uruchomimy _flatten_. A tutaj to ma te klamrowe nawiasy, które są dynamicznym closure. I jest po prostu jakiś arbitralny kod tutaj, który zostanie wykonany w kontekście operatora view.

Tutaj to mówi: weź greeting, które jest wejściem operatora view, i to jest tutaj. Mógłbym nazwać to, jak chcę, mógłbym nazwać to _"foo"_ i po prostu muszę się do tego odwoływać jako _"foo"_ później. A potem mówię: z tym zwróć to.

A potem ustawiam zwracanie stringa, który mówi before the flatten dla zmiennej. Bardzo proste.

Teraz dodam kolejny dokładnie taki sam, ale powiem after _flatten_.

Więc co to robi, ponieważ to działa sekwencyjnie, zobaczysz, jak wygląda kanał, zanim uruchomimy _flatten_, a potem znowu po uruchomieniu _flatten_.

A potem ten kanał greeting jest nadal utworzony, więc nadal będzie przekazany do procesu. I mam nadzieję, że teraz workflow się uruchomi. Spróbujmy.

Świetnie. Po pierwsze, pipeline tym razem się nie zawiesił. Mieliśmy trzy procesy, które uruchomiły się prawidłowo i mamy mały ptaszek. A potem możemy zobaczyć, że nasze instrukcje view zadziałały.

Mamy before _flatten_, która jest tą tablicą, którą widzieliśmy wcześniej z błędu, a potem mamy trzy razy after _flatten_ został wywołany, gdzie mamy _"Hello", "Bonjour"_ i wszystkie te trzy oddzielne elementy w tablicy, które są teraz, jak mieliśmy nadzieję, trzema oddzielnymi elementami w kanale.

I możesz zobaczyć, że operator _view_ został uruchomiony trzy razy. I to dlatego, że ten kanał po _flatten_ ma teraz trzy elementy. I więc operator zostaje wywołany trzy razy.

Bardzo szybko, wspomniałbym tylko, że gdy wcześniej tworzyłem fabryki kanałów, zrobiłem _"."_, i wtedy zobaczyliśmy, że było wiele różnych sposobów tworzenia kanałów, a jeden z nich nazywa się "_fromList"_. I to jest faktycznie specjalnie zaprojektowane, żeby robić tę samą operację. Więc mogliśmy po prostu zrobić fromList greetings array i to zadziała. To nieco czystsza i ładniejsza składnia. Ale dla celów tej demonstracji chcieliśmy to zrobić bardziej krok po kroku, żebyś mógł zobaczyć, jak kanał jest manipulowany i jak różne operatory mogą zmieniać zawartość kanału.

## 4. Odczytywanie wartości wejściowych z pliku CSV

Dobrze, jak możemy to uczynić bardziej realistycznym? Prawdopodobnie nie będziesz chciał tworzyć mnóstwa kodu w swoim pipeline'ie Nextflow'a z zakodowanymi na sztywno tablicami. Prawdopodobnie będziesz chciał wziąć dane z zewnątrz, gdy uruchamiasz, a te dane prawie na pewno będą w plikach.

Więc następną rzeczą, którą zrobimy, jest to, że zreplikujemy to, ale zamiast brać dane z pojedynczego parametru CLI lub z zakodowanego na sztywno stringa lub tablicy, weźmiemy je z pliku.

Więc pozbądźmy się naszego greetings array. A teraz zmienimy tę fabrykę kanałów ponownie. Właśnie powiedziałem, że jest kilka do wyboru i jest jedna zwana _".fromPath"_. I powiem jej, żeby w tym przypadku wzięła _params.input_, co wraca do naszego inputu, którego używaliśmy wcześniej.

Teraz ten parametr nie jest jeszcze gotowy do użycia. Nadal mówimy, że to string i jest zakodowany na sztywno tutaj z domyślną wartością, ale moglibyśmy nadpisać ten string. Teraz chcemy, żeby to był plik. Więc typ jest inny. To nie jest już _String_. To _Path_.

A potem możemy ustawić domyślną wartość, jeśli chcemy, znowu na Path. A jeśli spojrzę w explorer po lewej, możesz zobaczyć w tym repozytorium, w tym katalogu roboczym, mam katalog o nazwie data. Mam tam plik o nazwie _"greetings.csv"._

Więc mogę po prostu ustawić domyślną wartość tutaj na _"data/greetings.csv"_. Teraz, gdy uruchomię ten pipeline ponownie bez żadnych opcji linii poleceń, użyje tej domyślnej wartości. Wie, że to ścieżka, więc wie, że powinien to obsługiwać jako ścieżkę, a nie string.

A potem przekaże to do fabryki kanałów z tego _params.input_ i utworzy nasz kanał, który następnie będzie użyty w tym procesie zwanym _sayHello_. Spróbujmy.

Dobrze. Nie powiodło się. Nie martw się. To było oczekiwane. A jeśli śledzisz materiały szkoleniowe, zobaczysz, że tam też było oczekiwane. Zobaczmy, co się dzieje.

Próbował uruchomić pipeline. Próbował wykonać proces i dostał dość podobny błąd do tego, który widzieliśmy wcześniej.

Tutaj mówi: próbowaliśmy uruchomić _echo_, ale zamiast wypisać zawartość tego pliku CSV, po prostu wypisał ścieżkę. I możesz zobaczyć, że to pełna ścieżka bezwzględna tutaj do tego pliku CSV.

A potem oczywiście, ponieważ próbował zapisać to do tej naprawdę skomplikowanej ścieżki, nie bardzo wiedział, co zrobić. I było to poza zakresem katalogu roboczego procesu.

Wspomniałem na początku, że Nextflow enkapsuluje każde wykonane zadanie w specjalnym katalogu roboczym. A jeśli próbujesz zapisać dane, które są poza tym katalogiem roboczym, Nextflow Cię zatrzyma jako środek ostrożności. I to właśnie się tutaj stało. Próbowaliśmy zapisać do ścieżki bezwzględnej i Nextflow zawiódł i nas powstrzymał.

## 4.2. Używanie operatora splitCsv() do parsowania pliku

Dobrze, spójrzmy na ten kanał i zobaczmy, jak wygląda. Możemy zrobić _".view"_, a skopiowałem to ze strony. Więc _.view_, i mamy dynamiczne closure tutaj i mówimy nazwa zmiennej "_csv"_ jako wejście. Więc to jest zawartość kanału, i mówimy before splitCsv, i tak to wygląda.

Jeśli uruchomię to ponownie, nadal się nie powiedzie, ale pokaże nam, co jest w tym kanale. To nie jest szczególnie ekscytujące. To jest ta zmienna _path_. Więc możesz zobaczyć, że to po prostu string tutaj, ponieważ jest wypisywany do terminala, ale to jest obiekt _path_, który zawiera informacje i metadane o tym pliku.

Nie chcemy przekazywać metadanych pliku do wejścia. Chcemy przekazać zawartość tego pliku. Jeśli spojrzymy na plik _greetings.csv_, możesz zobaczyć tutaj, że ma te różne zmienne. _Hello, Bonjour, Holà_ znowu. I to są rzeczy, które naprawdę chcemy przekazywać do naszego procesu, a nie tylko sam plik jako pojedynczy obiekt.

Więc musimy sparsować ten plik CSV. Musimy go rozpakować, dostać się do zawartości pliku CSV, a następnie przekazać zawartość w kanale do procesu.

Jak prawdopodobnie możesz wywnioskować z komunikatu logowania, chcemy użyć _splitCsv_, który jest kolejnym operatorem, kolejnym operatorem kanału. Więc jeśli zrobię "_dot" "s"_, a potem zobaczysz, że jest automatycznie sugerowany. Ups, _splitCsv_ i jakieś nawiasy.

A potem po _splitCsv_ wstawię kolejną instrukcję _view_, żebyśmy mogli zobaczyć, jak to wygląda potem. Uruchommy pipeline i zobaczmy, co mamy.

Dobrze. Nadal się nie powiodło, ale w nowy i ekscytujący sposób, co jest postępem.

Tym razem znowu mamy jakiś problem z naszym skryptem, który został wyrenderowany. Teraz nie mamy już końcowej ścieżki, ale mamy tablicę zmiennych, która wygląda bardzo podobnie do błędu, który mieliśmy wcześniej, gdy przekazywaliśmy tablicę jako stałe wejście.

Z naszym logowaniem z operatora view możemy zobaczyć, że before _splitCsv_ była ścieżka. I rzeczywiście, po _splitCsv_ mamy trzy różne wyjścia i każde z tych wyjść wygląda bardzo podobnie do każdego z wierszy z pliku _greetings.csv_, co ma sens.

Więc co się tutaj stało, to Nextflow sparsował ten plik CSV, dał nam trzy obiekty, jedną tablicę dla każdej linii pliku CSV. Więc potem trzy razy przekazaliśmy tablicę zmiennych do kanału zamiast pojedynczej wartości stringa.

Dobrze, więc ostatnim razem, gdy mieliśmy ten problem, użyliśmy _flatten_. Spróbujmy bardzo szybko flatten i zobaczmy, co się stanie.

Mogę nazwać te zmienne, jak chcę. Więc nazwę to _myarray_, bo to już nie jest naprawdę CSV. Spróbujmy uruchomić to ponownie i zobaczmy, co się stanie z _flatten_.

Więc tym razem uruchomimy, sparsowaliśmy CSV na trzy obiekty tablicowe, a potem je spłaszczyliśmy. I tym razem przeszło. I pipeline Nextflow'a się uruchomił. Jednak możesz zobaczyć, że _flatten_ naprawdę się rozszalał i spłaszczył wszystko. I więc dostajemy trzy niezależne wpisy tablicowe dla każdego wiersza. I więc uruchomił proces trzy razy dla każdego wiersza CSV. I teraz mamy całą masę plików wynikowych, i 123, 456, i wszelkiego rodzaju rzeczy, a nie tylko tę pierwszą kolumnę CSV, której naprawdę chcieliśmy.

## 4.3. Używanie operatora map() do wyodrębnienia powitań

Więc jak dostać się tylko do pierwszej kolumny? Jeśli flatten jest tutaj zbyt uproszczony, potrzebujemy bardziej złożonego operatora, gdzie możemy faktycznie dostosować i powiedzieć mu, czego chcemy z CSV.

Żeby to zrobić, użyjemy _map_. W zasadzie _map_ po prostu mówi: uruchom jakiś kod, jakąś funkcję na każdym elemencie, który dostaję i zrób jakąś transformację na nim. A ponieważ jest tak elastyczny, zobaczysz, że pojawia się w kodzie Nextflow'a cały czas.

Sam w sobie nic nie robi. Więc nie chcemy zwykłych nawiasów, chcemy tutaj closure i musimy powiedzieć mu, co robić. Więc powiem _"row"_, bo dostaje wiersze z CSV, więc to logiczna nazwa zmiennej. To wejście. I chcę zwrócić tylko pierwszy element tej tablicy.

Tablice w Nextflow'ie są indeksowane od zera, więc powiemy tylko pierwszy element, który jest row zero. Gdybym chciał drugą kolumnę, mógłbym być jeden lub trzecią kolumnę być dwa, i tak dalej. Możemy zwrócić tutaj, co chcemy, ale zwrócę tylko pierwszą wartość.

I teraz możemy uruchomić pipeline ponownie i zobaczyć, czy robi to, czego oczekujemy.

I rzeczywiście, po _splitCsv_ mamy nasze tablice, a potem po _map_ mamy nasze ładne czyste stringi, tylko _"Hello", "Bonjour"_ i _"Holà"_. I pipeline teraz robi to, czego chcemy. Fantastycznie.

Więc możemy teraz pozbyć się wszystkich tych poleceń view. Nie potrzebujemy ich już.

## Podsumowanie

Skończyliśmy nasze debugowanie i to jest kod, z którym kończymy. Bierzemy nasz parametr CLI zwany _input_, który jest sklasyfikowany jako _Path_. Nextflow znajduje ścieżkę, ładuje ją i rozumie plik CSV. Zwraca wszystkie różne wiersze. A potem mapujemy tylko pierwszy element tego wiersza do kanału, który daje nam zawartość kanału, która jest przekazywana do procesu.

A proces działa na każdym elemencie w kanale, których jest trzy. I uruchamia proces trzy razy, dając mu trzy zadania. A te wyniki są następnie publikowane z workflow'u, przechwytywane przez wyjście procesu. Publikowane z workflow'u i zapisywane w bloku output do podkatalogu o nazwie _"hello_channels"_.

Całkiem fajnie. Zbliżamy się teraz do czegoś, co bardziej przypomina prawdziwy pipeline Nextflow'a, który mógłbyś uruchomić do jakiejś prawdziwej analizy.

## Podsumowanie

Dobrze. Mam nadzieję, że teraz masz wyczucie, czym są kanały i operatory Nextflow'a i jak operatory działają na kanałach i jak możesz je tworzyć.

Kanały, jak powiedziałem na początku tego wideo, są klejem Nextflow'a. I możesz zobaczyć tutaj, że możemy brać różne wejścia i manipulować nimi i brać te dane, a następnie przekazywać je do logiki workflow'u niżej.

A ten blok workflow tutaj jest naprawdę miejscem, gdzie budujesz całą tę paralelizację i całą sprytną logikę, i wyjaśniasz Nextflow'owi, jak zbudować Twój DAG workflow'u i jak orkiestrować Twój pipeline.

Kanały nie są najłatwiejszym konceptem do zrozumienia. Więc zrób sobie przerwę, pomyśl trochę o tym, może przeczytaj materiał ponownie i naprawdę upewnij się, że masz te koncepcje opanowane, ponieważ to jest kluczowe dla Twojego zrozumienia Nextflow'a i im lepiej rozumiesz kanały i różne operatory kanałów i różne fabryki kanałów, tym więcej frajdy będziesz miał pisząc w Nextflow'ie i tym potężniejsze będą Twoje pipeline'y.

To nie jest to samo co zwykłe programowanie w Pythonie lub innych językach. Nie używamy tutaj instrukcji _if_, to jest funkcyjne programowanie przepływowe używające kanałów i operatorów. Więc to trochę inne, ale też super potężne.

To koniec tego rozdziału. Idź i zrób sobie krótką przerwę, a zobaczę Cię w następnym wideo dla części trzeciej, gdzie przejdziemy przez Hello Workflow i porozmawiamy trochę więcej o workflow'ach.

Tak jak w poprzednim rozdziale, jest kilka pytań quizowych na dole tej strony internetowej, więc możesz szybko przez nie przejść i upewnić się, że rozumiesz wszystkie różne części materiału, który właśnie przeszliśmy. A poza tym zobaczę Cię w następnym wideo. Dziękuję bardzo.

Dobrze.
