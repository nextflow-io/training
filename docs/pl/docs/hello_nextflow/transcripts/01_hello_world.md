# Część 1: Hello World - Transkrypcja wideo

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Tłumaczenie wspomagane przez AI - [dowiedz się więcej i zasugeruj ulepszenia](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

<div class="video-wrapper">
  <iframe width="560" height="315" src="https://www.youtube.com/embed/tOukLxWCHiA?si=F0t9LFYLjAWoyRXj&amp;list=PLPZ8WHdZGxmWKozQuzr27jyMGqp9kElVK&amp;cc_load_policy=1" title="YouTube video player" frameborder="0" allow="accelerometer; autoplay; clipboard-write; encrypted-media; gyroscope; picture-in-picture; web-share" referrerpolicy="strict-origin-when-cross-origin" allowfullscreen></iframe>
</div>

!!!note "Ważne uwagi"

    Ta strona zawiera tylko transkrypcję. Pełne instrukcje krok po kroku znajdziesz w [materiałach szkoleniowych](../01_hello_world.md).

    Numery sekcji pokazane w transkrypcji mają charakter orientacyjny i mogą nie obejmować wszystkich numerów sekcji w materiałach.

## Powitanie

Cześć i witaj ponownie.

Jesteś teraz w części pierwszej szkolenia "Hello Nextflow" o nazwie "Hello World". W tym rozdziale zaczniemy budować podstawową wiedzę o Nextflow'ie.

Mam nadzieję, że masz już skonfigurowane środowisko w Codespaces lub w równoważnym miejscu z uruchomionym VS Code i folder Hello Nextflow w przestrzeni roboczej w Eksploratorze ze wszystkimi tymi plikami.

Zaczniemy od wykonania bardzo prostych rzeczy w terminalu przy użyciu Bash'a, a następnie zobaczymy, czy możemy zrobić to samo w Nextflow'ie, żebyś poczuł, jak wygląda składnia.

## 0. Rozgrzewka

Zacznijmy naprawdę prosto. Po prostu użyjmy "echo", żeby wypisać coś do terminala. "Hello World". Naciskam enter i to trafia do terminala. Hello World. Mam nadzieję, że to nie jest niespodzianką dla nikogo oglądającego to szkolenie.

Dobra, zróbmy coś z tym. Zamiast po prostu wypisywać to do terminala, zapiszmy to do pliku. Nacisnę strzałkę w górę na klawiaturze, co przewija historię Bash'a, więc dostanę moje ostatnie polecenie, i dodam na końcu znak większości, który przekierowuje wyjście z tego polecenia do pliku, i nazwę go output.txt.

Enter ponownie, żeby uruchomić to polecenie, tym razem nic w terminalu, ale widzimy po lewej stronie, że pojawił się nowy plik o nazwie output.txt.

Możemy go wyświetlić w terminalu czymś takim jak cat. Więc cat output.txt i rzeczywiście mówi "Hello World". Możemy też kliknąć go dwukrotnie i otworzy się w edytorze kodu w VS Code.

## 1.1. Przeanalizuj kod

W porządku. Mówiłem, że to proste. Co dalej? Spróbujmy powtórzyć ten proces, ale tym razem zróbmy to w Nextflow'ie.

Jak mówiłem, wszystkie różne rozdziały w tym szkoleniu zaczynają się od skryptu i ten nazywa się Hello World. Więc znajdę Hello World. Podgląd pojawia się po pojedynczym kliknięciu, kliknę dwukrotnie, żeby otworzyć go w edytorze. I szybko pozbędę się terminala.

To bardzo prosty skrypt, tak prosty, jak to tylko możliwe. Ma tylko 22 linie długości i robi w zasadzie to samo. W rzeczywistości część tego powinna wyglądać znajomo. Widzimy nasze polecenie bash przekierowujące do pliku.

Dobra. Co jeszcze? Również w tym pliku możemy zacząć widzieć niektóre podstawowe koncepcje Nextflow'a. Mamy proces na czerwono tutaj i workflow. To specjalne słowa kluczowe i specjalna terminologia w Nextflow'ie.

## 1.1.1. Definicja procesu

Różne procesy w workflow'ie opakowują różne logiczne jednostki Twojego workflow'u. Każdy proces robi jedną rzecz.

Kiedy go uruchamiamy, generuje zadanie lub wiele zadań, które są faktycznymi krokami wykonywanymi w pipeline'ie. Wszystkie procesy są następnie orkiestrowane w bloku workflow, który widzimy na dole, i w tym przypadku po prostu uruchamia ten jeden proces.

Nazwa procesu następuje po tym słowie kluczowym tutaj i może być w zasadzie dowolna. A następnie zawartość procesu znajduje się w tych nawiasach klamrowych.

Jest tylko jeden wymóg dla procesu, a mianowicie to, że musi zawierać jakiś blok script lub exec. Jest to w potrójnych cudzysłowach tutaj i to jest skrypt bash, który zostaje zapisany do katalogu roboczego, gdy uruchamiamy pipeline'a, i to jest rzecz, która faktycznie działa na Twoim komputerze lub serwerze.

To jest zazwyczaj bash, ale możesz też umieścić tutaj inny hash bang na górze i może to być skrypt Pythona lub skrypt R. Nie ma znaczenia. Cokolwiek jest w tym skrypcie, zostanie wykonane.

Jest jeszcze jedna rzecz, którą dodaliśmy do tego procesu tutaj, a mianowicie deklaracja wyjścia. To mówi Nextflow'owi, że ten proces oczekuje pliku wyjściowego o nazwie output.txt. Mówi, że to jest path, więc powinien być obsługiwany jak plik, a nie powiedzmy, gdyby to było val, mówiłoby, że to jest jak zmienna lub wartość.

Zauważ, że to nie tworzy tego pliku. Nie generuje go faktycznie. To robi skrypt tutaj na dole. Po prostu mówi Nextflow'owi, żeby oczekiwał pliku wyjściowego o tej nazwie.

## 1.1.2. Definicja workflow'u

Dobra. A następnie na dole mamy workflow tutaj, i znowu mamy deklarację. Ten nazywa się Main. To jest odpowiednik bloku script dla workflow'u, jeśli chcesz. To jest część workflow'u, która coś robi. I w tym przypadku mówimy, wywołaj proces o nazwie sayHello.

Normalnie, oczywiście, Twój pipeline będzie wyglądał znacznie bardziej skomplikowanie. Prawdopodobnie będziesz miał więcej niż jeden proces i użyjesz kanałów do orkiestracji przepływu danych między nimi. Przejdziemy do tego w następnych częściach tego szkolenia, ale na razie to wystarczy. To jest poprawny pipeline, który powinien działać.

Mogę nawet kliknąć preview DAG tutaj w VS Code. DAG to reprezentacja struktury przepływu danych w pipeline'ie i możemy zobaczyć go wyrenderowanego z boku jako diagram mermaid. W tym przypadku jest bardzo prosty. Jest jedno pudełko, które jest workflow'em i jeden proces, który nazywa się sayHello, ale może to wyglądać ciekawiej w miarę postępów.

## 1.2. Uruchom workflow'a

Dobra, spróbujmy uruchomić tego workflow'a i zobaczmy, co się stanie.

Przywrócę terminal na dole, wyczyszczę wyjście i wpiszę Nextflow Run. A następnie po prostu wpiszę nazwę skryptu, która brzmi hello-world.nf. I nacisnę enter.

Dobra, ma standardowe rzeczy na górze, które mówią nam, że Nextflow został uruchomiony i która wersja działała i jaka była nazwa skryptu i wszystko.

I naprawdę ważną rzeczą, której szukamy tutaj, jest _tutaj_, co jest podsumowaniem różnych zadań, które zostały wykonane.

Jeśli Twoje wygląda tak z małym zielonym ptaszkiem, to gratulacje. Właśnie uruchomiłeś swój pierwszy pipeline. Fantastycznie.

Mówi nam tutaj nazwę procesu, który został uruchomiony, który nazywał się Say Hello, i powiedział nam, że uruchomił się raz i że był udany. To aktualizuje się w miarę postępów, więc gdy uruchamiasz większy pipeline, zobaczysz postęp reprezentowany tutaj. Ale ponieważ to jest tak małe, działa w zasadzie natychmiast.

## 1.2.2. Znajdź wyjście i logi w katalogu work

Teraz, gdy uruchamiasz pipeline Nextflow'a, każdy z tych procesów jest zszyty razem, i każdy proces, jak mówiłem wcześniej, może generować zadania - jedno lub wiele. Więc w tym przypadku mieliśmy pojedyncze zadanie z tego procesu. Po prostu uruchomił się raz i to zostało zrobione pod tym _hashem_ zadania.

Nextflow nie zajmuje się plikami w Twoim katalogu roboczym bezpośrednio, tworzy specjalny folder o nazwie work. I jeśli zrobię "ls", zobaczymy, że pojawił się tutaj: _work_, a w nim są podkatalogi dla każdego pojedynczego zadania, które się uruchamia. I to pasuje do tego hasha. Więc możesz zobaczyć, jeśli przejdę do "ls work/c4", a następnie jest skrócone, ale zaczyna się 203, i to jest katalog roboczy, który został utworzony przez ten proces, gdy uruchomiliśmy pipeline'a. I możesz go zobaczyć również z boku.

Gdy wylistuję te pliki, możesz zobaczyć, że plik output.txt został wygenerowany. Możesz go zobaczyć również tutaj. I jest kilka ukrytych plików, które nie pokazują się przy zwykłym "ls".

Jeśli kliknę na output.txt, rzeczywiście mamy nasze wyjście. Fantastycznie. Więc pipeline zadziałał.

Może się wydawać, że to dość dużo kodu dla uruchomienia tego, co było w zasadzie jednoliniowym skryptem bash, ale będzie to miało więcej sensu, gdy nasze procesy staną się bardziej skomplikowane. A ten katalog work z Nextflow'em i te pliki, które są tworzone, to naprawdę kręgosłup tego, co czyni Nextflow'a tak potężnym.

Każde zadanie, każdy element pipeline'u jest odizolowany od każdego innego zadania. Jest odtwarzalny. Nie kolidują ze sobą i wszystko może działać równolegle. To jest naprawdę miły sposób, gdy się do tego przyzwyczaisz, ze względu na tę izolację, że możesz wejść i zobaczyć dokładnie, co się stało dla pojedynczego zadania i debugować.

Rzućmy szybkie spojrzenie na te inne pliki w katalogu work. Od góry do dołu mamy plik o nazwie _.command.begin_. Jest pusty. To jest po prostu tak zwany plik wartowniczy, utworzony przez Nextflow'a mówiący, dobra, zaczynam zadanie. Nic ciekawego tam.

Następnie jest _.command.error_, _.command.log_ i _.command.out_. To wszystko są wyjścia z polecenia bash lub tego skryptu, który został uruchomiony. To jest standardowy błąd. To jest standardowe wyjście, i to są te dwa połączone w kolejności, w jakiej wyszły. Więc dostajesz logiczną kolejność.

Dobra, te wszystkie były również puste dla tego, więc niezbyt interesujące, ale rzeczy stają się ciekawsze, gdy dojdziesz do _.command.run_.

To jest zazwyczaj bardzo długi skrypt. I to jest to, co Nextflow faktycznie wykonuje. Jeśli tu wejdziesz, zaczniesz widzieć całą wewnętrzną logikę Nextflow'a i zobaczysz, co robi i jak wykonuje Twój proces. To będzie zależeć od tego, gdzie uruchamiasz, czy uruchamiamy lokalnie, czy przesyłamy to jako zadanie do SLURM'a, w którym to przypadku będziemy mieli nagłówki SLURM'a na górze. Wszystkie te różne konfiguracje.

Generalnie nie musisz nigdy zaglądać do tego pliku. Jest autogenerowany przez Nextflow'a i nie ma w nim nic szczególnie unikalnego dla Twojego pipeline'u. Ale to jest naprawdę rdzeń tego, co działa.

Następny jest znacznie ciekawszy. _.command.sh_ to wygenerowany skrypt, który pochodzi z Twojego procesu, i tutaj możesz zobaczyć, że Nextflow dodał nagłówek Bash, a następnie wykonał nasze polecenie, które było w naszym bloku script.

I to wszystko, co robi plik _.command.run_, po prostu uruchamia ten plik _.command.sh_.

To jest naprawdę przydatny, ten, na który zazwyczaj patrzysz najczęściej, gdy próbujesz coś debugować i sprawdzić, czy logika Twojego pipeline'u Nextflow'a robi to, czego oczekujesz.

Na koniec mamy plik o nazwie _.exitcode_, i to po prostu przechwytuje kod wyjścia z zadania, który w tym przypadku był udany. Więc kod wyjścia wynosił zero.

Jeśli coś pójdzie nie tak, zabraknie pamięci lub coś innego i się nie powiedzie, to jest bardzo przydatne do zrozumienia, co poszło nie tak.

## 1.3. Uruchom workflow'a ponownie

Jeszcze jedna rzecz do zrozumienia o katalogach work jest taka, że jeśli będę uruchamiał ten pipeline wielokrotnie, więc jeśli _"nextflow run hello-world.nf"_, zrobi dokładnie to samo, ale tym razem będzie miał nowy identyfikator zadania. Możesz zobaczyć, że ten hash tutaj jest inny, a teraz jeśli spojrzę w work, są dwa katalogi hash. I te są znowu oddzielne od siebie.

Więc za każdym razem, gdy uruchamiasz workflow Nextflow'a, chyba że użyjesz resume, który używa cache'u, do czego dojdziemy później, będzie ponownie uruchamiał te procesy w nowych katalogach work, które są oddzielne od siebie. Nie dostaniesz żadnych kolizji nazw plików, nie będziesz miał żadnych takich problemów. Wszystko jest odizolowane i czyste.

A jeśli wejdziemy do tego katalogu, możesz zobaczyć wszystkie te same pliki i ten sam _output.txt_, który został odtworzony od zera.

## 2. Publikuj wyjścia

Dobra, to świetne dla samego Nextflow'a, podczas gdy uruchamia Twój pipeline, żeby wszystkie rzeczy były oddzielone od siebie i czyste i mogły być zarządzane.

Ale nie jest to super przydatne, jeśli jesteś osobą próbującą eksplorować swoje wyniki. Nie chcesz naprawdę przeszukiwać tysięcy różnych katalogów work próbując znaleźć swoje pliki wynikowe. I nie powinieneś. Katalogi work nie są przeznaczone do bycia końcowym stanem, gdzie Twoje pliki są tworzone.

Robimy to poprzez publikowanie naszych plików.

## 2.1.1. Zadeklaruj wyjście procesu sayHello

Więc jeśli wrócę do naszego skryptu, będziemy pracować w naszym bloku workflow tutaj. Powiemy mu, jakich plików oczekiwać, które pliki nas interesują, a następnie utworzymy nowy blok poniżej o nazwie blok output.

To jest nowa składnia, która pojawiła się z parserem składni i jest domyślna w wersji 26.04 Nextflow'a. Więc jeśli używałeś Nextflow'a trochę wcześniej, to jest jedna z rzeczy, które są nowe.

Więc mamy blok main, a następnie powiem publish i powiem Nextflow'owi, czego oczekiwać od publikowania. Nazwiemy to _first_output_, i nazwiemy to _sayHello.out_.

Przypadkowo zrobiłem literówkę tam, ale to jest dobra okazja, żeby również wskazać niektóre funkcje rozszerzenia Nextflow VS Code. Możesz zobaczyć, że od razu dało mi małą czerwoną falistą linię pod tym mówiącą, że coś jest nie tak. A jeśli nawiodę na to kursorem, powie mi, że ta zmienna nie jest zdefiniowana. Nie wiem, co to jest.

To dość oczywiste w tym przypadku, zrobiłem literówkę. Chciałem wpisać sayHello, i wtedy falista linia znika.

Teraz jest fioletowa. Parser składni Nextflow'a wie, że to jest proces i gdy nawiodę na to kursorem, daje mi zredukowaną reprezentację tego, jak wygląda ten proces. Więc mogę bardzo szybko na pierwszy rzut oka zobaczyć, że nie przyjmuje żadnych wejść i daje nam to wyjście. Więc praca w VS Code z tym rozszerzeniem daje Ci wiele kontekstowych informacji podczas pisania kodu.

Zauważ, że możemy odwoływać się do wyjścia z tego procesu za pomocą składni _.out_. I w tej chwili możemy to nazwać, jak chcemy, to jest po prostu arbitralna nazwa zmiennej.

## 2.1.2. Dodaj blok output: do skryptu

Gdzie to staje się ważne, to gdy robimy nasz nowy blok tutaj, i to jest poniżej bloku workflow teraz, nie jesteśmy już wewnątrz workflow. Nawiasy klamrowe znowu. I to jest miejsce, gdzie po prostu mówimy Nextflow'owi, gdzie umieścić wszystkie pliki, które są tworzone przez workflow.

Teraz wezmę tę nazwę zmiennej, którą utworzyłem tutaj, i umieszczę ją tam i dam nawiasy klamrowe dla tego. I powiem Nextflow'owi, żeby użył path. Ups. Path, w cudzysłowie. I użyję kropki. To po prostu mówi Nextflow'owi, żeby umieścił plik w katalogu głównym katalogu results. Więc nie w żadnych podkatalogach ani nic.

Spróbujmy uruchomić naszego workflow'a ponownie. Jeśli zrobię _"nextflow run hello-world.nf"_, to miejmy nadzieję, że powinno wyglądać w zasadzie dokładnie tak samo. Nic się naprawdę nie zmieniło z Nextflow'em tutaj. Uruchamia te same rzeczy. Po prostu robi je w katalogach work ponownie.

Ale teraz jeśli zrobię _"ls results/"_, zobaczysz, że jest nowy katalog tutaj, który został utworzony o nazwie results, który jest domyślnym katalogiem bazowym dla publikowania workflow'u. A w nim jest plik o nazwie _output.txt_.

Jeśli zrobię _"ls -l results"_, zobaczysz, że to jest faktycznie dowiązanie symboliczne do katalogu work. Więc to nie jest prawdziwy plik, jest połączony z katalogiem work i zebrał wszystkie pliki tam dla nas.

## 2.2. Ustaw niestandardową lokalizację

"Results" to domyślna nazwa dla tej ścieżki. Jeśli uruchomię workflow'a ponownie, i tym razem zrobię _dash_ pojedynczy myślnik, to jest, bo to jest podstawowa opcja Nextflow'a. _"-output-dir **my**results"._ Mógłbym też po prostu zrobić _"-o"_ w skrócie. Wtedy ustawi inny katalog bazowy dla miejsca, gdzie pliki są przechowywane i jeszcze raz, tutaj w _myresults/_, teraz mamy _output.txt_.

To świetne, ale prawdopodobnie nie chcemy wszystkich plików po prostu w katalogu głównym. Chcemy trochę organizacji, więc możemy również utworzyć podkatalog tutaj o nazwie, jakiej chcemy. Powiedzmy _"path 'hello_world'"_, i po prostu uruchomię to ponownie. _"nextflow run hello-world.nf"_. Powinno pójść do katalogu results do podkatalogu i rzeczywiście, teraz pod results tutaj na górze mamy _hello_world/_ i mamy _output.txt_.

Ważna rzecz do zauważenia, stary plik _output.txt_ nadal tam jest. Katalog results nie jest czyszczony, gdy to robisz. Po prostu nowe pliki są tam kopiowane. Nadpiszą pliki, które już tam są, jeśli mają tę samą nazwę pliku, ale nie wyczyszczą starych. Więc musisz być trochę ostrożny, gdy ponownie uruchamiasz pipeline'y. Jeśli nie chcesz, żeby były na wierzchu plików, które już tam są. Upewnij się, że używasz pustego katalogu.

## 2.3. Ustaw tryb publikowania na copy

Dobra, wspomniałem, że te pliki są dowiązaniami symbolicznymi, więc jeśli zrobię _"ls -l results/hello_world/"_, możesz zobaczyć, że to jest dowiązanie symboliczne do katalogu work. To jest generalnie dobra rzecz, jeśli pracujesz na czymś takim jak HPC, i to są naprawdę ogromne pliki i nie chcesz ich duplikować, ponieważ oznacza to, że pliki są przechowywane tylko raz w systemie plików.

Jednak oznacza to, że jeśli usuniesz katalog work: jeśli zrobię _"rm -r work"_ i wyczyszczę wszystkie te pliki pośrednie, które zostały utworzone. Teraz, jeśli spróbuję przeczytać ten plik _"results/hello_world/"_. Będzie wskazywał jako dowiązanie symboliczne na plik, który już nie istnieje i dane zniknęły na zawsze i są nieodwracalne, co może nie być świetne.

Więc generalnie mówimy, że dobrą praktyką jest kopiowanie plików zamiast dowiązań symbolicznych, jeśli możesz, ponieważ jest to bezpieczniejsze. Po prostu bądź świadomy, że użyje dwa razy więcej miejsca na dysku, chyba że usuniesz te katalogi work.

Żeby to zrobić z blokiem output, przejdę do first output tutaj. Ustawiłem path wcześniej i teraz ustawię mode i możesz zobaczyć, gdy piszę, rozszerzenie VS code sugeruje rzeczy, wie, że to jest dyrektywa output tutaj. I powiem copy. Zapisuję.

Uruchommy workflow'a ponownie. Utworzy pliki ponownie, nowy katalog work.

Teraz, jeśli przejdę do _"ls -l results/hello_world/"_ możesz zobaczyć, że to jest prawdziwy plik i nie jest już dowiązaniem symbolicznym, i Nextflow to skopiował. Dobrze wiedzieć. Więc path i mode to rzeczy, które będziesz pisał dość często.

Teraz, oczywiście, to jest bardzo proste. Uczynimy to bardziej złożonym i potężnym w miarę postępów i zobaczysz, jak uczynić te rzeczy dynamicznymi i nie zbyt rozwlekłymi.

## 2.4. Uwaga o dyrektywach publishDir na poziomie procesu

Teraz, powiedziałem, gdy zaczynaliśmy to, że to jest dość nowa forma składni. Jest dostępna tylko w najnowszych wersjach Nextflow'a, gdy to nagrywam, i nazywa się Workflow Outputs.

Jeśli tego używasz, to świetnie. Odblokowuje wiele innych fajnych funkcji w Nextflow'ie, takich jak Nextflow Lineage, aby pomóc śledzić pochodzenie tych plików w miarę ich tworzenia, i wkrótce będzie domyślne w 26.04. A w późniejszym terminie w przyszłości będzie to jedyny sposób pisania Twoich workflow'ów.

Jednak, ponieważ jesteśmy teraz w tej fazie przejściowej, możesz dobrze zobaczyć pipeline'y w naturze, które używają czegoś zwanego publishDir, co jest starym sposobem robienia tego, i to jest zdefiniowane nie na poziomie workflow i output, ale jest zdefiniowane na poziomie procesu.

I ta deklaracja mówi w zasadzie to samo. Mówi, opublikuj pliki wynikowe do katalogu o nazwie results i użyj trybu copy. Więc możesz zobaczyć, że składnia jest bardzo podobna. Ale gdy piszesz nowe pipeline'y teraz, staraj się nie używać tej dyrektywy publishDir, nawet jeśli ją widzisz, w wynikach AI lub w dokumentacji lub innych pipeline'ach, ponieważ to jest stary sposób robienia tego.

W 2026 wszyscy powinniśmy używać workflow outputs.

To wszystko jest udokumentowane, jeśli to robisz i używałeś Nextflow'a wcześniej, możesz przejść do dokumentacji Nextflow tutaj, nextflow.io/docs/. I jeśli przewinę w dół do tutorials, jest tutorial o nazwie _Migrating to Workflow Outputs_.

Jest naprawdę dobry. Przechodzi przez całą składnię, jak jest równoważna starej składni, dlaczego to zmieniliśmy, i ma oś czasu i wszystko. I przechodzi przez wszystkie różne scenariusze z mnóstwem przykładów. Więc możesz łatwo przekonwertować istniejący kod Nextflow'a na nową składnię.

## 3.1. Zmień proces sayHello, żeby oczekiwał zmiennego wejścia

Dobra, więc mamy nasz prosty skrypt, który uruchamia proces, tworzy plik, mówi Nextflow'owi, że to jest wyjście, a następnie mówimy Nextflow'owi, gdzie zapisać ten plik. To dobry początek.

Ale byłoby ciekawiej, gdyby nie wszystko było zakodowane na stałe. Więc następnie pomyślmy o tym, jak powiedzieć Nextflow'owi, że ten proces może przyjąć zmienne wejście, które jest czymś, co możemy kontrolować w czasie uruchomienia, gdy uruchamiamy workflow'a.

Musimy zrobić kilka różnych rzeczy, żeby to się stało.

Po pierwsze, musimy powiedzieć temu procesowi, że może przyjąć zmienną wejściową i wpiszemy _input_ tutaj jako nowy blok deklaracji. I nazwiemy to _"val greeting"_.

Część val jest odpowiednikiem path tutaj na dole. Mówi Nextflow'owi, że to jest zmienna, jak string w tym przypadku. A jeśli nawiedziesz na to kursorem ponownie, powie Ci z rozszerzenia, co to oznacza.

Następnie powiemy Nextflow'owi, co z tym zrobić. Nie wystarczy po prostu powiedzieć, że jest zmienna. Musisz powiedzieć w skrypcie, jak użyć tej zmiennej. I więc pozbędę się tego zakodowanego na stałe stringa tutaj, i umieszczę zmienną.

Szybko zrobię to bez nawiasów klamrowych, żeby pokazać Ci, że to jest dozwolone, i to jest stary sposób robienia tego. Ale teraz z nową składnią naprawdę zalecamy umieszczanie tego w nawiasach klamrowych w ten sposób, i to sprawia, że jest naprawdę jasne, że to jest interpolowane przez Nextflow'a tutaj.

Świetnie. Więc _"input greeting"_ idzie do _$\{greeting\}._ Ostatnia rzecz to musimy powiedzieć Nextflow'owi na poziomie workflow, że ten proces teraz przyjmuje wejście. I żeby to zrobić, w zasadzie damy mu zmienną.

## 3.2. Skonfiguruj parametr wiersza poleceń do przechwytywania wejścia użytkownika

Moglibyśmy zakodować to na stałe ponownie, jak Hello World, i to by działało dobrze, ale oczywiście nie daje nam to naprawdę żadnej przewagi. Chcieliśmy być w stanie skonfigurować to w czasie uruchomienia, więc chcemy być w stanie zrobić to w CLI, gdy uruchamiamy Nextflow'a.

A sposób, w jaki to robimy, to specjalna koncepcja Nextflow'a zwana _params_. Nazwiemy to _params.input_.

To, co to robi, to eksponuje tę zmienną input w CLI i tam używamy podwójnego myślnika, gdy uruchamiamy Nextflow'a.

Mogę to nazwać, jak chcę, mogę to nazwać _hello, greeting_. Nie ma znaczenia. Cokolwiek tam zrobię, będzie eksponowane jako opcja CLI, gdy uruchamiamy pipeline'a. I to jest prawdziwa sztuczka magiczna Nextflow'a, ponieważ oznacza to, że możesz bardzo szybko zbudować swój skrypt workflow'u z tymi parametrami, i w zasadzie budujesz niestandardowe CLI dla swojego pipeline'u, czyniąc go naprawdę łatwym do dostosowania różnych opcji w locie, gdy uruchamiasz.

Więc. Spróbujmy tego. Wracam do naszego terminala. Mamy nasze polecenie _"nextflow run"_ tutaj. A teraz zrobię _"--input"_, co pasuje do _"params.input"_, które widzieliśmy wcześniej. Myślę, że w dokumentacji jest po francusku. Geraldine lubi mówić po francusku. Zrobię to po szwedzku, bo mieszkam w Szwecji. więc powiem, "_Hej Världen_" i nacisnę enter.

Można użyć pojedynczych cudzysłowów lub podwójnych cudzysłowów, to po prostu wpływa na to, jak Bash to interpretuje.

Uruchamia pipeline Nextflow'a dokładnie w ten sam sposób. Możesz zobaczyć katalog roboczy i wszystko jest takie samo. Ale teraz jeśli przejdę do _"results/hello_world/output"_. Możemy zobaczyć nasz ładny szwedzki tutaj zamiast tego.

Więc dynamicznie przekazaliśmy wejście z CLI do parametru. Przekazaliśmy to jako wejście do procesu i proces to zinterpretował i umieścił w bloku script, który następnie dynamicznie zmienił wyjście tego wyniku skryptu. Całkiem fajnie.

Dość złożona logika z bardzo małą składnią tutaj. I możesz mieć nadzieję zobaczyć, jak to teraz zaczyna się skalować. I tak naprawdę budujemy logikę i możliwość dostosowania naszych pipeline'ów do skryptu Nextflow'a.

## 3.4. Użyj wartości domyślnych dla parametrów wiersza poleceń

Dobra, to świetne. Problem jednak teraz jest taki, że za każdym razem, gdy uruchamiam ten pipeline, muszę zrobić dash, input, żeby się uruchomił.

Jeśli spróbuję uruchomić bez tego parametru, teraz Nextflow wyrzuci błąd mówiący, że potrzebował tego parametru i nie został ustawiony. i więc nie wiedział, co zrobić.

To jest fajne nowe coś, przy okazji. W przeszłości Nextflow po prostu uruchomiłby się z pustym stringiem i miałbyś wszelkiego rodzaju dziwne błędy, które byłyby trudne do zrozumienia. Ale w nowym parserze składni Nextflow'a jest trochę bardziej ostrożny i mówi Ci od razu.

Więc nie zawsze chcemy określać każdą pojedynczą opcję. Dobrą praktyką jest określanie rozsądnych wartości domyślnych. Więc jak to robimy w naszym skrypcie?

Zauważysz, że gdy to napisaliśmy, po prostu umieściliśmy _params.input_ bezpośrednio tam, gdzie go używamy. Więc oczywistym rozwiązaniem jest zdefiniowanie wartości domyślnej, i robimy to na górze skryptu tutaj w specjalnym bloku params w workflow'ie. To jest w skrypcie workflow tutaj.

Znowu, trochę nowej składni tutaj, więc zwróć uwagę. To są naprawdę fajne rzeczy. Mamy nazwę parametru, który będzie oczekiwany tutaj.

A następnie po tym znaku dwukropka definiujemy typ zmiennej. Nie musisz tego robić, możesz po prostu zostawić to puste, ale to jest naprawdę miłe. Mówi Nextflow'owi, że oczekujemy stringa i traktuj go jako taki.

Jeśli chcemy liczby zamiast tego, na przykład, moglibyśmy napisać float, i to by powiedziało, że chcemy liczby zmiennoprzecinkowej. A jeśli spróbujemy uruchomić z tym, wtedy wyrzuci błąd. Jeśli damy mu stringa, który nie jest float'em. I również przekaże go jako taki. Jeśli zrobimy string, wtedy wie, że to jest string. I nawet jeśli ma wiodące zera i jest całkowicie numeryczny, nadal przekaże go jako faktyczny string.

Więc to bezpieczeństwo typów jest bardzo nową funkcją Nextflow'a, ale naprawdę potężną, żeby uczynić Twój kod bezpieczniejszym do pisania i uruchamiania.

Następnie po tym mamy znak równości, a następnie wartość domyślną tutaj. Nextflow został napisany w Barcelonie pierwotnie, więc wydaje się odpowiednie, że mamy trochę hiszpańskiego tutaj, _"Holà mundo!"_ jako wartość domyślną.

Dobra, zapiszę ten skrypt, wrócę, uruchomię skrypt ponownie bez _--input_. I tym razem powinien się uruchomić i utworzy nasz nowy plik w _results_. I w tym pliku teraz mówi _"Holà mundo!"_.

To jest tylko wartość domyślna, więc nie oznacza to, że nie możemy nadal robić tego samego co wcześniej. Jeśli wrócę i znajdę mój stary skrypt tutaj, _"Hej Världen"_, ponieważ robię _--input_ w wierszu poleceń, to nadpisze tę wartość domyślną i użyje tego ponownie w pliku output.txt.

Więc to w skrypcie jest tylko wartością domyślną, którą ustawiam.

W miarę jak budujemy nasz workflow, żeby był bardziej złożony i zawierał więcej parametrów, ten blok params na górze skryptu zacznie zbierać je wszystkie w jednym miejscu.

I kończysz z tą całkiem ładną symetrią w swoim skrypcie, gdzie faktycznie masz wszystkie swoje wejścia workflow tutaj i wyjścia workflow na dole. I jest bardzo jasne, jaki jest interfejs Twojego workflow'u do świata zewnętrznego. Więc możesz bardzo szybko podnieść nowy pipeline z nową składnią i zrozumieć, jak go używać.

Jeszcze jedna ostatnia fajna rzecz. Nie musimy ustawiać wartości domyślnej z tym. Jeśli zrobimy params input, ale nie ustawimy wartości domyślnej, wtedy mówi Nextflow'owi, że ten parametr jest wymagany, i znowu, pipeline nie uruchomi się bez niego, ale da Ci bardziej użyteczny komunikat o błędzie zamiast czegoś o tym, że jest null.

Więc mówi, że oczekujemy, że jego input jest wymagany, ale nie został określony w wierszu poleceń. Bardzo miłe.

Dobra, więc miejmy nadzieję, że teraz jest jasne, jak skonfigurować Twój pipeline Nextflow'a ze zmiennymi wejściami i parametrami, jak ustawić wartość domyślną, ustawić typy, może to być Boolean true false flaga lub integer lub różne typy tutaj. Jak przekazać je do Twojego workflow'u, gdzie to przechodzi, a następnie interpoluje do Twojego procesu. A także wiesz, jak dostosować te w wierszu poleceń, gdy uruchamiasz Nextflow'a. To zaczyna wyglądać ciekawiej niż nasze proste polecenie bash.

## 4. Zarządzaj wykonaniami workflow'u

Dobra. Co dalej? W ostatniej części tego rozdziału porozmawiamy trochę o tym, jak zarządzać wszystkimi różnymi wykonaniami workflow'u. Jeśli spojrzysz w moim pasku bocznym tutaj i Eksploratorze pod work, zobaczysz, że uruchomiłem kilka różnych pipeline'ów i te katalogi work stają się dość długie, jest ich wiele.

A druga rzecz jest taka, że jak mówiłem wcześniej, za każdym razem, gdy ponownie uruchamiam ten pipeline, tworzy nowy zestaw katalogów work i ponownie uruchamia wszystkie procesy od zera, co jest dobrą rzeczą. To jest zamierzone zachowanie. Jest odtwarzalne i regeneruje wszystko na świeżo. Ale oczywiście, jeśli uruchamiasz bardzo długo działające procesy, irytujące jest zawsze musieć zaczynać swój pipeline od początku, jeśli się zawiesił w połowie drogi, lub jeśli zmienisz coś na końcu pipeline'u.

## 4.1. Uruchom ponownie workflow'a z -resume

Na szczęście Nextflow jest naprawdę dobry w wiedzy, co zostało wcześniej uruchomione i co jest dostępne, i ponowne użycie tych starych wyników jest bardzo proste. Po prostu dodajemy nową flagę na końcu polecenia _"-resume"_.

Teraz zauważ, że są dwa myślniki na input, bo to jest parametr. Jest tylko jeden myślnik na resume, ponieważ to jest podstawowa opcja Nextflow'a.

To wprawia ludzi w zakłopotanie cały czas, nawet jeśli używasz Nextflow'a od długiego czasu. Więc zawsze pamiętaj jeden lub dwa myślniki. Zależy, czy to jest podstawowa opcja Nextflow'a.

Dobra, więc teraz robię _-resume_ i uruchamiam dokładnie tego samego workflow'a ponownie. I tym razem powinno wyglądać prawie dokładnie tak samo z jedną kluczową różnicą.

W wyjściu tutaj możesz zobaczyć, że wyniki zostały zbuforowane. I w rzeczywistości ten hash zadania tutaj jest dokładnie taki sam jak poprzednie uruchomienie, i po prostu ponownie użył tego katalogu work w całości. Wejścia i wyjścia i skrypt były wszystkie niezmodyfikowane. I więc po prostu bierze ten plik stamtąd i jeśli są kroki downstream w procesie, przekazałby je do następnego kroku w pipeline'ie.

Więc nadal uruchamia cały pipeline od początku do końca, ale używa zbuforowanych wyników dla każdego z tych zadań, gdzie może.

Teraz, gdy robisz _-resume_, po prostu wznawia ostatnie uruchomienie pipeline'u w Twoim katalogu roboczym, cokolwiek to było. Ale możesz faktycznie wznowić z dowolnego poprzedniego uruchomienia, które tam zrobiłeś. I zrobiliśmy już całkiem sporo teraz.

## 4.2. Sprawdź log poprzednich wykonań

Żeby spojrzeć na wszystkie z nich, możemy zrobić _"nextflow log"_ zamiast _"nextflow run"_, i to da nam ładne wyjście pokazujące wszystkie te różne.. Muszę zmniejszyć mój ekran trochę, żebyśmy mogli to zobaczyć, wszystkie te różne uruchomienia, gdy je zrobiliśmy, identyfikator sesji, polecenie i wszystko.

I możemy tu zajrzeć i możemy wziąć nazwę uruchomienia któregokolwiek z nich i następnie wznowić jeden z tych konkretnych. Więc mogę wrócić i mogę wznowić ten o nazwie _hungry_ekeblad_. I po prostu umieszczam to po _resume_.

Jeśli jesteś ciekaw, przy okazji, wszystkie te przymiotniki i nazwy naukowców są w kodzie źródłowym Nextflow'a. To naprawdę dobry sposób, żeby dostać swój pierwszy pull request do Nextflow'a, idąc i znajdując to i dodając swojego ulubionego naukowca.

I tak czy inaczej, więc zrobiłem to i wrócił i spojrzał na zbuforowane wyniki z tego uruchomienia workflow'u, zdał sobie sprawę, że nadal może ich użyć, i zrobił to. Więc dostałem zbuforowane wyniki ponownie.

## 4.3. Usuń starsze katalogi work

To świetne. Co jeśli chcę wyczyścić te katalogi work? Jest ich mnóstwo tutaj. Jest mnóstwo plików. Może wiem na pewno, że chcę wznowić z ostatnich kilku uruchomień pipeline'u, ale nie obchodzą mnie wszystkie te sprzed tego.

Wtedy mogę wybrać jeden tutaj i mogę użyć innego polecenia Nextflow'a, które jest _"nextflow clean"_, i mogę zrobić _"nextflow clean"_, zrobię _"-before"_, i konkretną nazwę uruchomienia, która w tym przypadku była _reverent_pike_ i zrobię _"-n"_, co mówi Nextflow'owi, żeby po prostu zrobił suchy przebieg. Więc po prostu mówi mi, co usunie. Bez faktycznego robienia czegokolwiek, więc usunąłby te katalogi work.

To wygląda rozsądnie. Więc zrobię to samo polecenie ponownie, ale zamiast _"-n"_ zrobię _"-f"_, żeby faktycznie wykonać czyszczenie. I tym razem faktycznie usunął wszystkie te katalogi. A jeśli wejdę i spojrzę na katalogi work, teraz wygląda znacznie lżej. Fantastycznie.

Więc tak czyścisz wszystkie swoje lokalne katalogi work w dość bezpieczny sposób bez całkowitego niszczenia cache'u. Więc nadal możesz wznowić, jeśli chcesz.

Jeśli kiedykolwiek zapomnisz, jakie są te flagi dla każdego polecenia Nextflow'a, możesz zrobić _"nextflow help"_, a następnie nazwę polecenia. Więc jeśli zrobię _"nextflow help clean"_, możesz zobaczyć wszystkie różne opcje: _-after, -before, -but_, wszystkie różne sposoby konfiguracji tego zachowania czyszczenia. Całkiem fajnie.

## Podsumowanie

Dobra, to koniec części pierwszej Hello Nextflow. To dość intensywny początek szkolenia, ale miejmy nadzieję, że teraz masz całkiem dobre zrozumienie tego, jak wygląda skrypt Nextflow'a; z różnymi kluczowymi częściami, procesami, workflow'ami, wyjściami i parametrami. Wiesz, jak je skonfigurować z podstawowymi nadpisaniami z wiersza poleceń, jak zrobić dynamiczny blok input z dynamicznym skryptem i wiesz, jak zarządzać wszystkimi swoimi wykonaniami obciążenia: widzieć, co już uruchomiłeś, wznawiać, czyścić. Jest dużo rzeczy. Przeszedłeś długą drogę. Więc jeśli chcesz zrobić przerwę i mieć szybki spacer i filiżankę herbaty, teraz jest prawdopodobnie dobry czas. Zasłużyłeś na to.

Od tego momentu w zasadzie budujemy na tym fundamencie. Jak możemy uczynić to bardziej złożonym, bardziej potężnym? Jak możemy uczynić to bardziej elastycznym? Robić rzeczy, które chcemy zrobić, naszą analizę na skalę.

## Quiz

Teraz jeśli przewiniesz w dół do części pierwszej, hello world, na stronie internetowej zobaczysz mały quiz i to jest coś nowego, co zrobiliśmy dla tej wersji szkolenia Nextflow'a. I możesz przejść i sprawdzić się, żeby sprawdzić, czy zrozumiałeś cały materiał, który zrobiliśmy w tym rozdziale.

To nie jest wysyłane do nas ani nic, jest po prostu przechowywane w Twojej przeglądarce. Więc nie wiemy, jakie są Twoje odpowiedzi, ale to jest po prostu mała samokontrola, żeby upewnić się, że niczego nie przegapiłeś lub niczego nie źle zrozumiałeś. I możesz spróbować tyle razy, ile chcesz.

Jeśli jesteś jak ja, może chcesz zostać w terminalu w swojej instancji VS Code, w którym to przypadku możesz wpisać polecenie _quiz_, a następnie po prostu powiedzieć mu, w którym rozdziale jesteś. Więc robimy _"Hello World"_, a następnie możesz zrobić dokładnie te same pytania quizowe, które są w przeglądarce internetowej, ale po prostu w swoim terminalu.

Fajnie. Dobra. Mam nadzieję, że Ci się to podoba. Baw się trochę i zobaczymy się w następnym rozdziale za chwilę, żeby porozmawiać o kanałach Nextflow'a.
