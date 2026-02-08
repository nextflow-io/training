# Część 2: Hello Channels - Transkrypcja Wideo

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Tłumaczenie wspomagane przez AI - [dowiedz się więcej i zaproponuj poprawki](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

<div class="video-wrapper">
  <iframe width="560" height="315" src="https://www.youtube.com/embed/yDR66fzAMOg?si=xCItHLiOQWqoqBB9&amp;list=PLPZ8WHdZGxmWKozQuzr27jyMGqp9kElVK&amp;cc_load_policy=1" title="YouTube video player" frameborder="0" allow="accelerometer; autoplay; clipboard-write; encrypted-media; gyroscope; picture-in-picture; web-share" referrerpolicy="strict-origin-when-cross-origin" allowfullscreen></iframe>
</div>

!!! note "Ważne uwagi"

    Ta strona zawiera jedynie transkrypcję. Pełne instrukcje krok po kroku znajdziesz w [materiale szkoleniowym](../02_hello_channels.md).

    Numery sekcji pokazane w transkrypcji służą jedynie celom informacyjnym i mogą nie obejmować wszystkich numerów sekcji w materiałach.

## Witamy

Witaj z powrotem w Części 2 Hello Nextflow. Ten rozdział nosi tytuł Hello Channels.

Kanały są jak spoiwo Twojego pipeline'u Nextflow. To fragmenty, które trzymają razem wszystkie różne procesy, których Nextflow używa do przekazywania informacji i orkiestracji Twojego workflow'u.

Kanały mają jeszcze jeden aspekt - operatory. To zasadniczo funkcje, których możemy używać na kanałach, aby modyfikować ich zawartość. Przejdźmy do VS Code i zobaczmy, gdzie się znajdujemy.

Jestem bardzo przybliżony w tym VS Code, więc żeby zachować porządek, usunąłem wszystkie pliki _.nextflow\*_ oraz katalog _work/_, _results/_ i wszystko z Rozdziału Pierwszego. Po prostu zaczynam tu od nowa. Ale nie przejmuj się tym zbytnio. Jeśli nie chcesz, możesz zostawić te pliki. Nie spowodują żadnych problemów.

Zaczniemy od pracy nad _hello-channels.nf_ w tym rozdziale, a jeśli go otworzę, powinien wyglądać bardzo podobnie do pliku, nad którym pracowaliśmy wcześniej. Może się zdarzyć, że różne fragmenty znajdują się w innych miejscach skryptu, ale wszystko powinno być zasadniczo takie samo.

Jedną różnicą jest to, że ścieżka w bloku output to teraz _hello_channels_ dla tej części, co oznacza, że pliki wynikowe będą przechowywane w innym podkatalogu w Twoim katalogu results, jeśli nadal go tam masz. Powinno być więc miłe i czyste miejsce do rozpoczęcia bez dezorientacji związanej z wynikami.

Dobrze, przypomnijmy sobie szybko, co robi ten skrypt, gdy uruchamiamy ten workflow. Wykonujemy _"nextflow run hello-channels.nf"_. Możemy dodać _"--input myinput"_, a gdy to uruchomimy, użyje tego parametru params.input, który został przekazany jako zmienna do procesu sayHello tutaj na górze, która trafia do greeting i zostaje zapisana do output.txt. I możemy to zobaczyć w pliku results. Świetnie.

## 1. Dostarczanie zmiennych wejściowych za pomocą kanału jawnie

To miłe. Ale jest dość uproszczone. Mamy jedną zmienną w tym parametrze, która trafia do procesu uruchamianego raz i nie skaluje się zbyt dobrze. Nie możemy przekazać jej wielu różnych plików do utworzenia. Nie możemy jej dać wielu różnych powitań. Mamy tylko jedno.

W rzeczywistości Nextflow służy do skalowania Twojej analizy. Prawdopodobnie więc chcesz, żeby robił więcej niż jedną rzecz. I robimy to za pomocą _kanałów_.

Kanały to trochę unikalny koncept dla wielu osób zaczynających pracę z Nextflow. Wywodzi się z koncepcji programowania funkcyjnego i może zająć trochę czasu, zanim się go zrozumie, ale gdy już to zrozumiesz, naprawdę odblokują moc Nextflow i są kluczowe dla tego, jak piszesz swoje workflow'y.

## 1.1. Utworzenie kanału wejściowego

Zacznijmy od wzięcia tego skryptu i sprawienia, by używał _kanału_ zamiast tylko _parametru_.

Przechodzimy do workflow'u, gdzie znajduje się cała logika naszego workflow'u dotycząca łączenia rzeczy. I tutaj utworzę nowy kanał.

Utworzę nowy kanał.

I nazwę go "_greeting_ch"_. To konwencja, by dodawać "_\_ch"_ w ten sposób, żebyś mógł pamiętać, że ta zmienna jest kanałem. Ale możesz nazwać to, jak chcesz.

Potem napiszę równa się, i zrobię _"channel.of"._

Channel to jak przestrzeń nazw dla wszystkiego związanego z kanałami. Małe "c", jeśli używałeś wcześniej Nextflow. A _".of"_ to coś, co nazywa się fabryką kanałów (Channel factory), która jest zasadniczo sposobem tworzenia kanału.

Istnieje wiele różnych fabryk kanałów. Jeśli napiszę tylko ".", zobaczysz, że VS Code sugeruje ich mnóstwo, ale _".of"_ jest najprostszy i po prostu przyjmuje wartość wejściową.

Więc mogę dodać nawiasy i napiszę _"Hello Channels!"_.

Świetnie. Mam kanał. Fantastycznie. Mogę nacisnąć zapisz, mogę to uruchomić ponownie, ale nic ciekawego się nie stanie. VS Code pokazał mi pomarańczową linię ostrzeżenia tutaj i powiedział, że to jest ustawione: utworzyłeś to, ale nigdy tego nie użyłeś do niczego. Ten kanał nie jest wykorzystywany.

Dobrze, więc jak go użyjemy? Bardzo prosto. Skopiuję to, usunę _params.input_ i zamiast tego wstawię tutaj _"greeting_ch"_. Więc przekażemy ten kanał jako wejście do sayHello.

Zauważ, że na razie zakodowałem ten string na sztywno. To trochę krok wstecz po naszym ładnym parametrze, którego używaliśmy na końcu ostatniego rozdziału, ale pozwala to na prostotę na początek, żebyś mógł zobaczyć logikę.

Dobrze, przejdę do mojego terminala i uruchomię workflow ponownie. Tym razem bez żadnego _"--input"_, a on się uruchomi i użyje kanału, który utworzyliśmy, i miejmy nadzieję, że powinniśmy mieć plik tutaj w _results/hello_channels/_, który teraz mówi "Hello Channels!". Fantastycznie. Więc tego oczekujemy od naszego kanału. Świetnie.

## 1.4. Użycie view() do inspekcji zawartości kanału

Jeszcze jedna rzecz do dodania tutaj, szybkie wprowadzenie do kolejnej funkcji, której możemy używać na kanałach, zwanej "_.view"_.

To jest analogiczne do polecenia _print_ w Pythonie lub innych językach, do których możesz być przyzwyczajony, i po prostu wyrzuca zawartość tego kanału do terminala, gdy go uruchamiamy.

Więc zrobię "_.view"_, a potem jeśli ponownie uruchomię workflow, powinno wydrukować do terminala, jaka jest zawartość tego kanału w momencie jego utworzenia.

I rzeczywiście, widzisz, że jest to wydrukowane do terminala tutaj. _"Hello Channels!"_.

Zauważ, że możesz dzielić te rzeczy na wiele linii, jeśli chcesz, i w rzeczywistości automatyczny formater Nextflow będzie próbował to dla Ciebie robić. Białe znaki nie są tu naprawdę ważne, więc możesz łączyć te rzeczy jedną po drugiej.

## 2. Modyfikacja workflow'u do działania na wielu wartościach wejściowych

Dobrze, więc nasz kanał ma jedną rzecz, co jest miłe, ale jest to zasadniczo to samo, co było wcześniej. Więc skomplikujmy to trochę. Dodajmy kilka więcej rzeczy do naszego kanału.

Fabryka kanałów "_.of()"_ może przyjąć wiele elementów, więc napiszmy kilka więcej. Zrobimy _Hello, Bonjour, Hej_. I wtedy możemy uruchomić ten workflow ponownie i zobaczymy, co się stanie.

Powinno się uruchomić ponownie. I teraz wydrukowaliśmy _"Hello", "Bonjour"_ i _"Hej"_ do terminala za pomocą naszej instrukcji view. Fantastycznie.

## 2.1.2. Uruchomienie polecenia i sprawdzenie logów wyjściowych

Możesz myśleć, że w tym momencie skończyliśmy. Ale właściwie jest tu pewna pułapka, która nas złapie. Jeśli spojrzymy na nasz plik wyjściowy tutaj, zobaczysz, że ma _"Hello"_, ale nie ma żadnych innych wyników. W rzeczywistości jest tylko ten jeden.

Jeśli uruchomimy ten workflow wielokrotnie, możemy nawet zobaczyć, że czasami ma _"Bonjour"_, czasami _"Hej"_. To trochę losowe.

Jeśli spojrzymy na terminal, możemy zobaczyć, że uruchomił się trzy razy i możemy zobaczyć różne wyniki view. Ale jeśli przejdę do katalogu work, mogę zrobić _"cat work"_. Wstaw ten hash i rozwiń to oraz _output.txt_. Widzisz, że ten plik w katalogu work jest inny niż katalog results, a ten to _"Hej"._ Więc coś tu nie działa całkiem prawidłowo.

Kluczowe jest to, że mamy trzy zadania, które się uruchomiły. Output Nextflow próbuje to podsumować w miarę postępu przetwarzania, aby nie przejął całkowicie Twojego terminala, a to logowanie ANSI używa kodów escape ANSI, co zasadniczo nadpisało inne zadania. Więc pokazuje Ci tylko ostatnie, które zostało zaktualizowane.

## 2.1.3. Ponowne uruchomienie polecenia z opcją -ansi-log false

Jest kilka rzeczy, które możemy zrobić, aby lepiej to zrozumieć. Możemy zajrzeć do samego katalogu work i zobaczyć wszystkie różne katalogi work tam, ale to jest trochę mylące, ponieważ będą one zmieszane z różnymi uruchomieniami Nextflow.

Albo możemy powiedzieć Nextflow, żeby nie używał kodów escape ANSI.

Więc jeśli uruchomię polecenie ponownie, ale tym razem powiem _"-ansi-log false"_, żeby to wyłączyć, mógłbym też użyć zmiennych środowiskowych _$NO_COLOR_ lub _"$NXF_ANSI_LOG=false"_. Wtedy używa bardziej staromodnego stylu logowania Nextflow bez żadnych z tych kodów escape. Po prostu drukuje bezpośrednio do terminala bez żadnych sprytnych aktualizacji.

I teraz możemy zobaczyć wszystkie trzy procesy, które się uruchomiły. I każdy z nich ma swój własny hash zadania. A jeśli wejdziemy do tych katalogów work, zobaczymy trzy różne powitania, które określiliśmy.

Więc to ma teraz więcej sensu. Miejmy nadzieję, że rozumiesz, że Nextflow to robił, po prostu był trochę sprytny z tym, co pokazywał Ci w terminalu z tymi katalogami work.

Jednak to naprawiło jeden problem z katalogami work, ale nie naprawiło problemu z plikiem wyjściowym. Nadal mamy tylko jeden plik wyjściowy, który mówi _"Hello"_.

## 2.2. Upewnienie się, że nazwy plików wyjściowych będą unikalne

Teraz, aby to zrozumieć, musimy wrócić do naszego skryptu workflow. Generujemy tutaj nasz kanał, przekazujemy go do naszego procesu, a jeśli spojrzymy na proces, zapisujemy greeting do pliku o nazwie _"output.txt"_ i przekazujemy ten plik wyjściowy z powrotem do bloku output tutaj na dole, publikując go.

Jednak za każdym razem, gdy ten proces uruchamia się te trzy razy - te trzy różne zadania - wszystkie generują plik o nazwie _"output.txt"_, wszystkie te pliki wyjściowe są publikowane w katalogu results i wszystkie się wzajemnie nadpisują. Więc jakikolwiek plik wynikowy tam dostaniesz, jest po prostu ostatnim, który został wygenerowany, ale zniszczył wszystkie pozostałe. To nie jest to, czego chcemy.

## 2.2.1. Konstruowanie dynamicznej nazwy pliku wyjściowego

Istnieją różne sposoby radzenia sobie z tym, ale najprostszym na razie jest po prostu utworzenie różnych unikalnych nazw plików. Więc za każdym razem, gdy zadanie uruchamia się z innym powitaniem, wygeneruje inny plik wyjściowy, który nie będzie już kolidował podczas publikowania. I wtedy otrzymamy trzy unikalne pliki wyjściowe.

Robimy to dokładnie w ten sam sposób. Możemy użyć tej zmiennej gdziekolwiek w bloku script i możemy jej użyć wiele razy.

Więc mogę wkleić to tutaj, _"$\{greeting\}\_output.txt"_, a następnie muszę również wkleić to tutaj na górze, ponieważ nie tworzymy już pliku o nazwie _output.txt_. Więc jeśli tego nie zaktualizuję, Nextflow ulegnie awarii z błędem mówiącym, że oczekiwał pliku, który nigdy nie został wygenerowany.

Więc muszę zrobić to samo tam i muszę użyć podwójnych cudzysłowów, a nie pojedynczych, żeby ta zmienna była zrozumiana.

Dobrze, wypróbujmy to i zobaczmy, czy zadziałało. Uruchomimy workflow ponownie. Miejmy nadzieję, że pokaże nam trzy różne zadania w trzech różnych katalogach work. I rzeczywiście, widzisz tutaj na górze w folderze results po lewej stronie. Mamy teraz trzy różne pliki z trzema różnymi nazwami plików, każdy z różną zawartością, jakiej oczekujemy. Więc pliki już się nie nadpisują i wszystko jest tam, jak oczekujemy.

To jest trochę trywialna konfiguracja, przez którą przeszliśmy tutaj, ale podkreśla niektóre kluczowe koncepcje, które musisz zrozumieć o tym, jak działa publikowanie plików, i niektóre pułapki, w które możesz wpaść. Więc miejmy nadzieję, że możesz tego uniknąć w swoich własnych workflow'ach.

Warto również zauważyć, że to, co tutaj zrobiliśmy, jest trochę niepraktyczne w rzeczywistych sytuacjach. Wzięliśmy jakieś dane wejściowe i używamy tych danych, ale również nazywamy plik po tych danych, czego zwykle nie można robić.

Więc w rzeczywistych, bardziej dojrzałych pipeline'ach Nextflow często będziesz przekazywał obiekt meta ze wszystkimi metadanymi powiązanymi z daną próbką. Możesz wtedy tworzyć dynamiczne nazwy plików na podstawie tego, co jest o wiele bardziej praktyczne.

Jeśli jesteś zainteresowany, jak to zrobić zgodnie z najlepszymi praktykami, na _training.nextflow.io_ jest quest poboczny, który dotyczy konkretnie metadanych i map meta, więc możesz tam zagłębić się w więcej szczegółów.

## 3. Dostarczanie wielu danych wejściowych za pomocą tablicy

Dobrze. Następnie zbadamy trochę, jak kanały są zbudowane i jak różnią się od innych rodzajów struktur danych w języku programowania. I pomyślę trochę o tym, jak mogę potencjalnie użyć tablicy, która może być znajomą koncepcją, jeśli przyszedłeś z innych języków.

Czy mogę użyć tablicy w kanale? Spróbujmy tego. Utworzę tablicę, a to skopiowałem z dokumentacji, _"greetings_array"_ i _"Hello", "Bonjour"_ i _"Holà"_. A potem wstawię to tutaj zamiast moich zakodowanych na sztywno stringów. Więc powiem "channel.of" _"greetings_array"_, przekazując tę tablicę do kanału. Spróbujmy tego.

Wywołam terminal i uruchomię pipeline.

Dobrze. Widzisz, że instrukcja view tutaj wydrukowała naszą tablicę zgodnie z oczekiwaniami, ale potem cały ten czerwony tekst, lub nie będzie czerwony, jeśli nadal masz _"-ansi-log"_ wyłączony, ale cały ten czerwony tekst mówi nam, że coś poszło nie tak.

Nie mamy już ładnego zielonego znacznika. Mamy czerwony krzyżyk, a jeśli po prostu zrobię to trochę szersze, żeby było łatwiej czytać, Nextflow mówi nam, co poszło nie tak.

Więc rozłóżmy to na części. Mówi, że błąd został spowodowany przez, a następnie przyczynę błędu, którą są brakujące pliki wyjściowe. Więc zasadniczo ten blok output powiedział, że ten plik powinien zostać utworzony, a nie został. Następnie mówi, że to jest polecenie, które zostało wykonane. Więc to jest zasadniczo zawartość tego pliku _.command.sh_. Tak to wyglądało po wstawieniu wszystkich tych zmiennych.

I widzisz tutaj, że nasze polecenie echo zostało faktycznie uruchomione tylko raz i użyło całej tablicy, ale w reprezentacji stringowej, co nie jest tym, czego chcieliśmy.

A potem polecenie zakończyło się w ten sposób, i to był katalog work, gdzie możemy pójść i zobaczyć pliki, aby lepiej to zrozumieć.

Dobrze. Więc co się stało, to Nextflow po prostu przekazał całą tę tablicę jako pojedynczy element kanału do procesu, co oznaczało, że proces uruchomił się tylko raz. Miał jedno zadanie i nie użył danych w strukturze, jakiej oczekiwaliśmy.

## 3.2. Użycie operatora do transformacji zawartości kanału

Więc musimy zrobić coś z tym kanałem najpierw, zanim będzie mógł być użyty. I to przygotowuje grunt pod użycie operatorów, które są specjalnymi funkcjami, których możemy używać na kanałach do manipulowania zawartością kanału.

W tym przypadku użyjemy czegoś, co nazywa się _flatten_. Które przekazujemy na końcu kanału tutaj. Więc tworzymy kanał, a następnie uruchamiamy _flatten_. I znowu, jeśli najedziesz na to myszką, pokaże Ci dokumentację dla tego polecenia od razu w VS Code, co jest bardzo pomocne. Możesz również znaleźć wszystkie te dokumenty na stronie Nextflow, w dokumentacji.

Mógłbym po prostu uruchomić ten kod teraz i zobaczyć, czy działa, ale to również dobra okazja, aby wprowadzić, jak robić dynamiczny kod w operatorach i w kodzie Nextflow, które nazywane są closures.

Więc dodam z powrotem polecenie view tutaj, zanim uruchomimy _flatten_. I tutaj ten ma te kręcone nawiasy, które są dynamicznym closure. I jest tam tylko jakiś dowolny kod, który zostanie wykonany w kontekście operatora view.

Tutaj to mówi weź greeting, które jest wejściem operatora view, i to jest tutaj. Mógłbym nazwać to, jak chcę, mógłbym nazwać to _"foo"_ i po prostu muszę się do tego później odwoływać jako _"foo"_. A potem mówię z tym, zwróć to.

A następnie ustawiamy zwracanie stringa, który mówi przed flatten dla zmiennej. Bardzo proste.

Teraz dodam kolejny dokładnie taki sam, ale powiem po _flatten_.

Więc to co to robi, ponieważ to działa sekwencyjnie, zobaczysz, jak wygląda kanał, zanim uruchomimy _flatten_, a następnie ponownie po uruchomieniu _flatten_.

A potem ten kanał greeting jest nadal utworzony, więc nadal będzie przekazany do procesu. I miejmy nadzieję, że teraz workflow się uruchomi. Spróbujmy tego.

Świetnie. Po pierwsze, pipeline tym razem się nie zawiesił. Mieliśmy trzy procesy, które uruchomiły się prawidłowo i mamy mały znacznik. A następnie możemy zobaczyć, że nasze instrukcje view zadziałały.

Mamy przed _flatten_, która jest tą tablicą, którą widzieliśmy wcześniej z błędu, a następnie mamy trzy razy po _flatten_ została wywołana, gdzie mamy _"Hello", "Bonjour"_ i te wszystkie trzy oddzielne elementy w tablicy, które są teraz, jak mieliśmy nadzieję, trzema oddzielnymi elementami w kanale.

I widzisz, że operator _view_ został uruchomiony trzy razy. I to dlatego, że ten kanał po _flatten_ ma teraz trzy elementy. I więc operator jest wywoływany trzy razy.

Bardzo szybko, wspomnę tylko, że gdy tworzyłem fabryki kanałów wcześniej, zrobiłem _"."_, a następnie zobaczyliśmy, że było wiele różnych sposobów tworzenia kanałów, a jeden z nich nazywa się "_fromList"_. I to jest faktycznie specjalnie zaprojektowane, aby wykonać tę samą operację. Więc mogliśmy po prostu zrobić from list greetings array, i to by zadziałało. To nieco czystsza i ładniejsza składnia. Ale dla celów tej demonstracji chcieliśmy zrobić to bardziej krok po kroku, żebyś mógł zobaczyć, jak kanał jest manipulowany i jak różne operatory mogą zmienić zawartość kanału.

## 4. Odczytywanie wartości wejściowych z pliku CSV

Dobrze, jak możemy to uczynić bardziej realistycznym? Prawdopodobnie nie będziesz chciał tworzyć mnóstwa kodu w swoim pipeline Nextflow z zakodowanymi na sztywno tablicami. Prawdopodobnie będziesz chciał pobrać dane z zewnątrz podczas uruchamiania, a te dane prawie na pewno będą w plikach.

Więc następną rzeczą, którą zrobimy, będzie replikacja tego, ale zamiast brać dane z pojedynczego parametru CLI lub z zakodowanego na sztywno stringa lub tablicy, weźmiemy je z pliku.

Więc pozbądźmy się naszych greetings array. A teraz zmienimy tę fabrykę kanałów ponownie. Właśnie powiedziałem, że jest mnóstwo do wyboru i jest jedna zwana _".fromPath"_. I powiem jej, żeby w tym przypadku wzięła _params.input_, co wraca do naszego inputu, którego używaliśmy wcześniej.

Teraz ten parametr nie jest jeszcze gotowy do użycia. Nadal mówimy, że to jest string i jest zakodowany na sztywno tutaj z domyślną wartością, ale możemy nadpisać ten string. Teraz chcemy, żeby to był plik. Więc typ jest inny. To już nie jest _String_. To _Path_.

A następnie możemy ustawić wartość domyślną, jeśli chcemy, ponownie na Path. A jeśli spojrzę w explorer po lewej stronie, zobaczysz w tym repozytorium, w tym katalogu roboczym, że mam katalog o nazwie data. Mam tam plik o nazwie _"greetings.csv"._

Więc mogę po prostu ustawić wartość domyślną tutaj na _"data/greetings.csv"_. Teraz, gdy uruchomię ten pipeline ponownie bez żadnych opcji wiersza poleceń, użyje tej wartości domyślnej. Wie, że to path, więc wie, że powinno to traktować jako ścieżkę, a nie string.

A następnie przekaże to do fabryki kanałów z tego _params.input_ i utworzy nasz kanał, który następnie będzie użyty w tym procesie o nazwie _sayHello_. Spróbujmy tego.

Dobrze. Nie powiodło się. Nie martw się. To było oczekiwane. A jeśli śledzisz materiał szkoleniowy, zobaczysz, że tam też było oczekiwane. Zobaczmy, co się dzieje.

Próbował uruchomić pipeline. Próbował wykonać proces i dostał dość podobny błąd do tego, który widzieliśmy wcześniej.

Tutaj mówi: próbowaliśmy uruchomić _echo_, ale zamiast echować zawartość tego pliku CSV, po prostu echowało ścieżkę. I widzisz, że to jest pełna ścieżka absolutna tutaj do tego pliku CSV.

A następnie, rzecz jasna, ponieważ próbował zapisać to do tej naprawdę skomplikowanej ścieżki, nie bardzo wiedział, co zrobić. I było to poza zakresem katalogu work procesu.

Wspomniałem na początku, że Nextflow enkapsuluje każde wykonane zadanie w specjalnym katalogu work. A jeśli spróbujesz zapisać dane poza tym katalogiem work, Nextflow Cię zatrzyma jako środek ostrożności. I to właśnie się tutaj stało. Próbowaliśmy zapisać do ścieżki absolutnej, a Nextflow zawiódł i nas powstrzymał.

## 4.2. Użycie operatora splitCsv() do parsowania pliku

Dobrze, przyjrzyjmy się temu kanałowi i zobaczmy, jak wygląda. Możemy zrobić _".view"_, a to skopiowałem ze strony internetowej. Więc _.view_, i mamy dynamiczne closure tutaj i mówimy nazwę zmiennej "_csv"_ jako wejście. Więc to jest zawartość kanału, i mówimy przed splitCsv, i tak to wygląda.

Jeśli uruchomię to ponownie, nadal się nie powiedzie, ale pokaże nam, co jest w tym kanale. Nie jest to szczególnie ekscytujące. To jest ta zmienna _path_. Więc widzisz, że to jest po prostu string tutaj, ponieważ jest drukowany do terminala, ale jest to obiekt _path_, który zawiera informacje i metadane o tym pliku.

Nie chcemy przekazywać metadanych pliku do wejścia. Chcemy przekazać zawartość tego pliku. Jeśli spojrzymy na plik _greetings.csv_, zobaczysz tutaj, że ma te różne zmienne tutaj. _Hello, Bonjour, Holà_ ponownie. I to są rzeczy, które naprawdę chcemy przekazywać do naszego procesu, a nie tylko sam plik jako pojedynczy obiekt.

Więc musimy sparsować ten plik CSV. Musimy go rozpakować, dostać się do zawartości pliku CSV, a następnie przekazać zawartość w kanale do procesu.

Jak prawdopodobnie możesz wywnioskować z komunikatu dziennika, chcemy użyć _splitCsv_, który jest kolejnym operatorem, kolejnym operatorem kanału. Więc jeśli zrobię "_dot" "s"_, a następnie zobaczysz, że jest autosugerowany. Ups, _splitCsv_ i nawiasy.

A następnie po _splitCsv_, wstawię kolejną instrukcję _view_, żebyśmy mogli zobaczyć, jak to wygląda później. Uruchommy pipeline i zobaczmy, co dostaliśmy.

Dobrze. Nadal się nie powiodło, ale w nowy i ekscytujący sposób, co jest postępem.

Tym razem znowu mamy jakiś problem z naszym skryptem, który został wyrenderowany. Teraz nie mamy już końcowej ścieżki, ale mamy tablicę zmiennych, która wygląda bardzo podobnie do błędu, który mieliśmy wcześniej, gdy przekazywaliśmy tablicę jako stałe wejście.

Dzięki naszemu logowaniu z operatora view możemy zobaczyć, że przed _splitCsv_ była ścieżka. I rzeczywiście, po _splitCsv_, mamy trzy różne wyjścia, a każde z tych wyjść wygląda bardzo podobnie do każdego z wierszy z pliku _greetings.csv_, co ma sens.

Więc to, co się tutaj stało, to Nextflow sparsował ten plik CSV, dając nam trzy obiekty, jedną tablicę dla każdej linii pliku CSV. Więc trzy razy przekazaliśmy tablicę zmiennych do kanału zamiast pojedynczej wartości stringa.

Dobrze, więc ostatnim razem, gdy mieliśmy ten problem, użyliśmy _flatten_. Spróbujmy bardzo szybko flatten i zobaczmy, co się stanie.

Mogę nazwać te zmienne, jak chcę. Więc nazwę to _myarray_, ponieważ to już tak naprawdę nie jest CSV. Spróbujmy uruchomić to ponownie i zobaczmy, co się stanie z _flatten_.

Więc tym razem uruchomimy, sparsowaliśmy CSV na trzy obiekty tablicowe, a następnie spłaszczyliśmy to. I tym razem to przeszło. I pipeline Nextflow się uruchomił. Jednak widzisz, że _flatten_ naprawdę się rozkręcił i spłaszczył wszystko. I więc otrzymujemy trzy niezależne wpisy tablicy dla każdego wiersza. I więc uruchomił proces trzy razy dla każdego wiersza CSV. A teraz mamy całą masę plików wynikowych, 123, 456 i wszystkie rodzaje rzeczy, a nie tylko tę pierwszą kolumnę CSV, której naprawdę chcieliśmy.

## 4.3. Użycie operatora map() do wydobycia powitań

Więc jak dostać się tylko do pierwszej kolumny? Jeśli flatten jest tu zbyt uproszczony, potrzebujemy bardziej złożonego operatora, gdzie możemy faktycznie dostosować i powiedzieć mu, czego chcemy z CSV.

Aby to zrobić, użyjemy _map_. Zasadniczo _map_ po prostu mówi, uruchom jakiś kod, jakąś funkcję nad każdym elementem, który dostaję, i wykonaj jakąś transformację na nim. A ponieważ jest tak elastyczny, zobaczysz, że pojawia się w kodzie Nextflow cały czas.

Sam w sobie nic nie robi. Więc nie chcemy zwykłych nawiasów, chcemy tutaj closure i musimy mu powiedzieć, co robić. Więc powiem _"row"_, ponieważ to dostaje wiersze z CSV, więc to logiczna nazwa zmiennej. Jest wejściem. I chcę zwrócić tylko pierwszy element tej tablicy.

Tablice w Nextflow są indeksowane od zera, więc powiemy tylko pierwszy element, którym jest wiersz zero. Gdybym chciał drugą kolumnę, mógłbym być jeden lub trzecia kolumna byłaby dwa, i tak dalej. Możemy zwrócić tutaj, co chcemy, ale zwrócę tylko pierwszą wartość.

A teraz możemy uruchomić pipeline ponownie i zobaczyć, czy robi to, czego oczekujemy.

I rzeczywiście, po _splitCsv_ mamy nasze tablice, a następnie po _map_, mamy nasze ładne, czyste stringi, tylko _"Hello", "Bonjour"_ i _"Holà"_. I pipeline teraz robi to, czego chcemy. Fantastycznie.

Więc możemy teraz pozbyć się wszystkich tych poleceń view. Nie potrzebujemy ich już.

## Podsumowanie

Skończyliśmy naszą debugowanie i to jest kod, z którym kończymy. Bierzemy nasz parametr CLI o nazwie _input_, który jest sklasyfikowany jako _Path_. Nextflow znajduje ścieżkę, ładuje ją i rozumie plik CSV. Zwraca wszystkie różne wiersze. A następnie mapujemy tylko pierwszy element tego wiersza do kanału, który daje nam zawartość kanału, która jest przekazywana do procesu.

A proces działa na każdym elemencie w kanale, których jest trzy. I uruchamia proces trzy razy, dając mu trzy zadania. A te wyniki są następnie publikowane z workflow'u, przechwycone przez wyjście procesu. Publikowane z workflow'u i zapisywane w bloku output do podkatalogu o nazwie _"hello_channels"_.

Całkiem fajne. Zbliżamy się teraz do czegoś, co bardziej przypomina prawdziwy pipeline Nextflow, który możesz uruchomić dla jakiejś prawdziwej analizy.

## Podsumowanie

Dobrze. Miejmy nadzieję, że teraz rozumiesz, czym są kanały i operatory Nextflow i jak operatory działają na kanałach oraz jak możesz je tworzyć.

Kanały, jak powiedziałem na początku tego wideo, są spoiwem Nextflow. I widzisz tutaj, że możemy brać różne dane wejściowe i manipulować nimi oraz pobierać te dane, a następnie przekazywać je do logiki downstream workflow.

A ten blok workflow tutaj jest naprawdę miejscem, gdzie budujesz całą paralelizację i całą sprytną logikę i wyjaśniasz Nextflow, jak zbudować Twój DAG workflow i jak zorkiestrować Twój pipeline.

Kanały nie są najłatwiejszą koncepcją do zrozumienia. Więc zrób sobie przerwę, pomyśl o tym trochę, może przeczytaj materiał ponownie i naprawdę upewnij się, że zrozumiałeś te koncepcje, ponieważ jest to kluczowe dla Twojego zrozumienia Nextflow i im lepiej rozumiesz kanały oraz różne operatory kanałów i różne fabryki kanałów, tym więcej zabawy będziesz mieć pisząc Nextflow i tym potężniejsze będą Twoje pipeline'y.

To nie jest to samo, co regularne programowanie w Pythonie lub innych językach. Nie używamy tutaj instrukcji _if_, to jest funkcjonalne programowanie przepływu przy użyciu kanałów i operatorów. Więc jest trochę inne, ale jest też super potężne.

To koniec tego rozdziału. Idź i zrób sobie krótką przerwę, a zobaczę się z Tobą w następnym wideo dla części trzeciej, gdzie przejdziemy przez Hello Workflow i porozmawiamy trochę więcej o workflow'ach.

Tak jak w poprzednim rozdziale, na dole strony internetowej tutaj jest kilka pytań quizowych, więc możesz przez nie szybko przejść i upewnić się, że rozumiesz wszystkie różne części materiału, który właśnie przeszliśmy. A poza tym zobaczę się z Tobą w następnym wideo. Dziękuję bardzo.

Dobrze.
