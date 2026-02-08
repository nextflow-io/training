# Część 1: Hello World - Transkrypcja wideo

<div class="video-wrapper">
  <iframe width="560" height="315" src="https://www.youtube.com/embed/tOukLxWCHiA?si=F0t9LFYLjAWoyRXj&amp;list=PLPZ8WHdZGxmWKozQuzr27jyMGqp9kElVK&amp;cc_load_policy=1" title="YouTube video player" frameborder="0" allow="accelerometer; autoplay; clipboard-write; encrypted-media; gyroscope; picture-in-picture; web-share" referrerpolicy="strict-origin-when-cross-origin" allowfullscreen></iframe>
</div>

!!!note "Ważne uwagi"

    Ta strona zawiera wyłącznie transkrypcję. Pełne instrukcje krok po kroku znajdziesz w [materiale szkoleniowym](../01_hello_world.md).

    Numery sekcji pokazane w transkrypcji służą jedynie celom informacyjnym i mogą nie obejmować wszystkich numerów sekcji z materiałów.

## Powitanie

Witaj ponownie.

Jesteś teraz w części pierwszej kursu "Hello Nextflow" o nazwie "Hello World". W tym rozdziale zaczniemy budować zrozumienie absolutnych podstaw Nextflow'a.

Miejmy nadzieję, że jesteś już skonfigurowany w Codespaces lub równoważnym środowisku z uruchomionym VS Code i masz folder Hello Nextflow w przestrzeni roboczej w Eksplorerze ze wszystkimi tymi różnymi plikami tutaj.

Zaczniemy od wykonania bardzo podstawowych rzeczy w terminalu przy użyciu Bash, a następnie zobaczymy, czy potrafimy zrobić to samo w Nextflow, żebyś poczuł, jak wygląda składnia.

## 0. Rozgrzewka

Zacznijmy naprawdę od prostych rzeczy. Po prostu zacznijmy od "echo", aby wydrukować coś w terminalu. "Hello World". Naciskam enter i to trafia do terminala. Hello World. Mam nadzieję, że to nie jest zaskoczeniem dla nikogo oglądającego ten kurs.

Dobra, zróbmy z tym coś. Zamiast po prostu drukować to w terminalu, zapiszmy to do pliku. Nacisnę strzałkę w górę na klawiaturze, co przewija historię Bash, więc daje mi ostatnie polecenie, i dodam na jego końcu tam mały znak większości, który przekierowuje wyjście z tego polecenia do pliku, i nazwę go output.txt.

Enter ponownie, aby uruchomić to polecenie, tym razem nic w terminalu, ale możemy zobaczyć po lewej stronie, że pojawił się nowy plik, nazwany output.txt.

Możemy wyświetlić go w terminalu czymś takim jak cat. Więc cat output.txt i rzeczywiście mówi "Hello World". Możemy też go dwukrotnie kliknąć i otworzy się w edytorze kodu w VS Code.

## 1.1. Zbadaj kod

W porządku. Mówiłem Ci, że to proste. Co dalej? Spróbujmy teraz wziąć ten proces i zrobić to ponownie, ale tym razem zróbmy to w Nextflow.

Jak powiedziałem, wszystkie różne rozdziały w tym kursie zaczynają się od skryptu i ten nazywa się Hello World. Więc znajdę Hello World. Pokazuje jego podgląd, gdy kliknę go raz, dwukrotnie kliknę, aby otworzyć go w edytorze tutaj. I szybko pozbędę się terminala.

To jest bardzo prosty skrypt, więc jeden z najprostszych jakie mogą być. Ma tylko 22 linijki długości i robi w zasadzie to samo. W rzeczywistości niektóre z tego powinny wyglądać znajomo. Możemy zobaczyć nasze polecenie bash przekierowujące do pliku tam.

Dobrze. Co jeszcze? Również w tym pliku możemy zacząć widzieć niektóre z podstawowych koncepcji Nextflow'a. Mamy proces na czerwono tutaj i workflow. To są specjalne słowa kluczowe i specjalna terminologia w Nextflow.

## 1.1.1. Definicja procesu

Różne procesy w workflow opakowują różne jednostki logiczne Twojego workflow'u. Każdy proces robi jedną rzecz.

Kiedy go uruchamiamy, generuje zadanie lub wiele zadań, które są faktycznymi krokami wykonywanymi w pipeline. Wszystkie procesy są następnie orkiestrowane w bloku workflow, który widzimy na dole, i w tym przypadku po prostu uruchamia ten jeden proces.

Nazwa procesu następuje po tym słowie kluczowym tutaj i może to być w zasadzie cokolwiek. A zawartość procesu znajduje się w tych klamrach.

Istnieje naprawdę tylko jeden wymóg dla procesu, którym jest to, że zawiera jakiś rodzaj bloku script lub exec. Jest to w potrójnych cudzysłowach tutaj i to jest skrypt bash, który jest zapisywany do katalogu roboczego, gdy uruchamiamy pipeline, i jest to rzecz, która faktycznie działa na Twoim komputerze lub serwerze.

To zazwyczaj bash, ale możesz też umieścić tutaj inny hash bang na górze i może to być skrypt Pythona lub skrypt R. Nie ma to znaczenia. Cokolwiek jest w tym skrypcie, zostanie wykonane.

Jest jeszcze jedna rzecz, którą dodaliśmy do tego procesu tutaj, którą jest deklaracja wyjścia. To mówi Nextflow, że ten proces oczekuje pliku wyjściowego o nazwie output.txt. Mówi, że jest to path, więc powinien być traktowany jak plik, a nie powiedzmy, gdyby to było val, powiedziałoby, że jest jak zmienna lub wartość.

Zauważ, że to nie tworzy tego pliku. Nie generuje go faktycznie. To robi skrypt tutaj na dole. Po prostu mówi Nextflow, aby oczekiwał pliku wyjściowego z tą nazwą.

## 1.1.2. Definicja workflow

Dobrze. A na dole mamy tutaj workflow, i znowu mamy deklarację. Ta nazywa się Main. To jest odpowiednik bloku skryptu dla workflow, jeśli wolisz. To część workflow, która coś robi. I w tym przypadku mówimy, wywołaj proces o nazwie sayHello.

Oczywiście normalnie Twój pipeline będzie wyglądał znacznie bardziej złożenie. Prawdopodobnie będziesz mieć więcej niż jeden proces i użyjesz kanałów do orkiestracji przepływu danych między nimi. Przejdziemy do tego w kolejnych częściach tego kursu, ale na razie to wystarczy. To jest ważny pipeline, który powinien działać.

Mogę nawet kliknąć preview DAG tutaj w VS Code. DAG jest reprezentacją struktury przepływu danych w pipeline, i możemy zobaczyć go wyrenderowanego po stronie jako diagram mermaid. W tym przypadku jest bardzo prosty. Jest jedno pole, które jest workflow i jeden proces, który nazywa się sayHello, ale to może wyglądać ciekawiej w miarę postępów.

## 1.2. Uruchom workflow'a

Dobrze, spróbujmy uruchomić ten workflow i zobaczmy, co się stanie.

Przywołam znowu terminal na dole, wyczyszczę wyjście i wpiszę Nextflow Run. A następnie po prostu wpiszę nazwę skryptu, która to hello-world.nf. I nacisnę enter.

Dobrze, ma standardowe rzeczy na górze, które mówią nam, że Nextflow został uruchomiony, która wersja była uruchomiona, jaka była nazwa skryptu i wszystko inne.

I naprawdę ważną rzeczą, której szukamy tutaj, jest _tutaj_, co jest podsumowaniem różnych zadań, które zostały wykonane.

Jeśli Twoje wygląda tak z małym zielonym ptaszkiem, to brawo. Właśnie uruchomiłeś swój pierwszy pipeline. Fantastycznie.

Mówi nam tutaj nazwę procesu, który się uruchomił, który nazywał się Say Hello, i powiedział nam, że uruchomił się raz i że zakończył się sukcesem. To aktualizuje się w miarę postępu, więc gdy uruchamiasz większy pipeline, zobaczysz postęp reprezentowany tutaj. Ale ponieważ to jest takie małe, uruchamia się zasadniczo natychmiast.

## 1.2.2. Znajdź wyjście i logi w katalogu roboczym

Teraz gdy uruchamiasz pipeline Nextflow, każdy z tych procesów jest zestawiany razem, i każdy proces, jak powiedziałem wcześniej, może generować zadania jedno lub wiele. Więc w tym przypadku mieliśmy pojedyncze zadanie z tego procesu. Po prostu uruchomił się raz i to zostało zrobione pod tym _hashem_ zadania.

Nextflow nie zajmuje się bezpośrednio plikami w Twoim katalogu roboczym, tworzy specjalny folder o nazwie work. I jeśli zrobię "ls", zobaczymy, że pojawił się tutaj: _work_, a wewnątrz są podkatalogi dla każdego pojedynczego zadania, które się uruchamia. I to pasuje do tego hasza. Więc możesz zobaczyć, jeśli przejdę do "ls work/c4", a następnie jest obcięte, ale zaczyna się 203, i to jest katalog roboczy, który został utworzony przez ten proces, gdy uruchomiliśmy pipeline. I możesz go też zobaczyć po stronie.

Gdy wymieniam te pliki, możesz zobaczyć, że plik output.txt został wygenerowany. Możesz go zobaczyć tutaj również. I jest kilka ukrytych plików, które nie pokazują się przy zwykłym "ls".

Jeśli kliknę na output.txt, rzeczywiście mamy nasze wyjście. Fantastycznie. Więc pipeline zadziałał.

Może się wydawać, że to dość dużo kodu dla uruchomienia czegoś, co było zasadniczo jednowierszowym skryptem bash, ale to nabierze więcej sensu, gdy nasze procesy staną się bardziej skomplikowane. A ten katalog roboczy z Nextflow i te pliki, które są tworzone, to naprawdę kręgosłup tego, co czyni Nextflow tak potężnym.

Każde zadanie, każdy element pipeline'u jest odizolowany od każdego innego zadania. Jest powtarzalny. Nie kolidują ze sobą i wszystko może działać równolegle. To jest w rzeczywistości naprawdę przyjemny sposób, gdy się do niego przyzwyczaisz, z powodu tej izolacji, że możesz wejść i zobaczyć dokładnie co się stało dla pojedynczego zadania i debugować.

Rzućmy szybkie spojrzenie na te inne pliki w katalogu roboczym. Od góry do dołu mamy plik o nazwie _.command.begin_. Jest pusty. To po prostu tak zwany plik wartowniczy, stworzony przez Nextflow mówiący, okej, zaczynam zadanie. Nic interesującego tam.

Następnie jest _.command.error_, _.command.log_ i _.command.out_. To są wszystkie wyjścia z polecenia bash lub tego skryptu, który się uruchomił. To jest standard error. To jest standard out, i to są te dwa połączone tak, jak wyszły. Więc dostajesz logiczną kolejność.

Dobrze, te wszystkie były też puste dla tego, więc niezbyt interesujące, ale rzeczy stają się bardziej interesujące, gdy dojdziesz do _.command.run_.

To jest zazwyczaj bardzo długi skrypt. I to jest to, co Nextflow faktycznie wykonuje. Jeśli tu wejdziesz, zaczniesz widzieć całą wewnętrzną logikę Nextflow i zobaczyć, co robi i jak wykonuje Twój proces. To będzie zależeć od tego, gdzie uruchamiasz, czy uruchamiamy lokalnie, czy przesyłając to jako zadanie do SLURM, w którym to przypadku będziemy mieli nagłówki SLURM na górze. Wszystkie te różne konfiguracje.

Ogólnie nie musisz nigdy zaglądać do tego pliku. Jest autogenerowany przez Nextflow i nie ma w nim nic szczególnie unikalnego dla Twojego pipeline'u. Ale to jest naprawdę rdzeń tego, co się uruchamia.

Następny jest znacznie bardziej interesujący. _.command.sh_ to wygenerowany skrypt, który pochodzi z Twojego procesu, i tutaj możesz zobaczyć, że Nextflow dodał nagłówek Bash, a następnie wykonał nasze polecenie, które było w naszym bloku skryptu.

I to wszystko, co robi plik _.command.run_, po prostu uruchamia ten plik _.command.sh_.

To jest naprawdę użyteczny, ten, na który zwykle patrzysz najwięcej, gdy próbujesz coś debugować i sprawdzać, czy logika Twojego pipeline'u Nextflow robi to, czego oczekujesz.

Na koniec mamy plik o nazwie _.exitcode_, i to po prostu przechwytuje kod wyjścia z zadania, który w tym przypadku zakończył się sukcesem. Więc kod wyjścia był zero.

Jeśli coś pójdzie nie tak, skończy Ci się pamięć lub coś innego i się nie powiedzie, to bardzo przydatne do zrozumienia, co poszło nie tak.

## 1.3. Uruchom workflow'a ponownie

Jeszcze jedna rzecz do zrozumienia o katalogach roboczych jest taka, że jeśli będę uruchamiał ten pipeline wielokrotnie, więc jeśli zrobię _"nextflow run hello-world.nf"_, zrobi dokładnie to samo, ale tym razem będzie miał nowy id zadania. Możesz zobaczyć, że ten hasz tutaj jest inny, i teraz jeśli spojrzę w work, są dwa katalogi hash. I te są znowu oddzielone od siebie.

Więc za każdym razem, gdy uruchamiasz workflow Nextflow, chyba że użyjesz resume, który używa cache, dotrzemy do tego później, będzie ponownie uruchamiał te procesy w nowych katalogach roboczych, które są od siebie oddzielone. Nie dostaniesz żadnych kolizji nazw plików, nie będziesz mieć żadnych takich problemów. Wszystko jest odizolowane i czyste.

A jeśli wejdziemy do tego katalogu, możesz zobaczyć wszystkie te same pliki i ten sam _output.txt_, który został odtworzony od zera.

## 2. Publikuj wyjścia

Dobrze, to świetne dla Nextflow'a dla niego samego, gdy uruchamia Twój pipeline, żeby wszystkie rzeczy były oddzielone od siebie i czyste i mogły być zarządzane.

Ale nie jest to super użyteczne, jeśli jesteś osobą próbującą eksplorować swoje wyniki. Nie chcesz naprawdę przekopywać się przez tysiące i tysiące różnych katalogów roboczych próbując znaleźć swoje pliki wynikowe. I nie powinieneś tego robić. Katalogi robocze nie są przeznaczone do bycia ostatecznym stanem, gdzie Twoje pliki są tworzone.

Robimy to publikując nasze pliki.

## 2.1.1. Zadeklaruj wyjście procesu sayHello

Więc jeśli wrócę do naszego skryptu, będziemy pracować w naszym bloku workflow tutaj. Powiemy mu, jakich plików oczekiwać, które pliki nas obchodzą, a następnie stworzymy nowy blok poniżej o nazwie blok wyjścia.

To jest nowa składnia, która pojawiła się z parserem składni i jest domyślna w wersji 26.04 Nextflow. Więc jeśli używałeś Nextflow trochę wcześniej, to jest jedna z rzeczy, które są nowe.

Więc mamy blok main, a następnie powiem publish i powiem Nextflow, czego oczekiwać od publikowania. Nazwiemy to _first_output_, i nazwiemy to _sayHello.out_.

Przypadkowo zrobiłem tam literówkę, ale to dobra okazja, aby też wskazać niektóre funkcje rozszerzenia Nextflow VS Code. Możesz zobaczyć, że od razu dało mi małą falującą czerwoną linię pod tym mówiąc, że coś jest nie tak. A jeśli na to najadę, powie mi, że ta zmienna nie jest zdefiniowana. Nie wiem, co to jest.

To dość oczywiste w tym przypadku, zrobiłem literówkę. Chciałem napisać sayHello, i wtedy falująca linia znika.

Teraz jest fioletowa. Parser składni Nextflow wie, że to jest proces i gdy na to najadę, daje mi zredukowaną reprezentację tego, jak ten proces wygląda. Więc mogę bardzo szybko na pierwszy rzut oka zobaczyć, że nie przyjmuje żadnych wejść i daje nam to wyjście. Więc praca w VS Code z tym rozszerzeniem daje Ci dużo kontekstowych informacji, gdy piszesz kod.

Zauważ, że możemy odnieść się do wyjścia z tego procesu ze składnią _.out_. I w tej chwili możemy nazwać to, jak chcemy, to tylko arbitralna nazwa zmiennej.

## 2.1.2. Dodaj blok output: do skryptu

Gdzie to staje się ważne, to gdy robimy nasz nowy blok tutaj, i to jest teraz poniżej bloku workflow, nie jesteśmy już wewnątrz workflow. Klamry ponownie. I tutaj po prostu mówimy Nextflow, gdzie umieścić wszystkie pliki, które są tworzone przez workflow.

Teraz wezmę tę nazwę zmiennej, którą stworzyłem tutaj, i umieszczę ją tam i dam niektóre klamry dla tego. I powiem Nextflow, aby użył path. Ups. Path, w cudzysłowach. I użyję kropki. To po prostu mówi Nextflow, aby umieścił plik w katalogu głównym katalogu results. Więc nie w żadnych podkatalogach ani nic.

Spróbujmy uruchomić nasz workflow ponownie. Jeśli zrobię _"nextflow run hello-world.nf"_, to miejmy nadzieję, że powinno wyglądać zasadniczo dokładnie tak samo. Nic naprawdę się nie zmieniło z Nextflow tutaj. Uruchamia te same rzeczy. Po prostu robi je w katalogach roboczych ponownie.

Ale teraz jeśli zrobię _"ls results/"_, zobaczysz, że jest nowy katalog tutaj, który został utworzony o nazwie results, który jest domyślnym katalogiem bazowym dla publikowania workflow. I tam jest plik o nazwie _output.txt_.

Jeśli zrobię _"ls -l results"_, zobaczysz, że to jest faktycznie soft link do katalogu roboczego. Więc to nie jest prawdziwy plik, jest połączony z katalogiem roboczym i zebrał dla nas wszystkie pliki tam.

## 2.2. Ustaw niestandardową lokalizację

"Results" to domyślna nazwa dla tej ścieżki. Jeśli uruchomię workflow ponownie, i tym razem zrobię _dash_ pojedynczy łącznik, to dlatego, że to jest podstawowa opcja Nextflow. _"-Output-dir **my** results"._ Mogę też po prostu zrobić _"-o"_ w skrócie. Następnie ustawi inny katalog bazowy, gdzie pliki są przechowywane i jeszcze raz, tutaj na górze w _myresults/_, teraz mamy _output.txt_.

To świetne, ale prawdopodobnie nie chcemy wszystkich plików tylko w katalogu głównym. Chcemy trochę organizacji, więc możemy też stworzyć podkatalog tutaj zwany czymkolwiek chcemy. Powiedzmy _"path 'hello_world'"_, i po prostu uruchomię to ponownie. _"nextflow run hello-world.nf"_. Powinno przejść do katalogu results do podkatalogu i rzeczywiście, teraz pod results tutaj na górze mamy _hello_world/_ i mamy _output.txt_.

Ważna rzecz do zauważenia, stary plik _output.txt_ jest nadal tam. Katalog results nie jest czyszczony, gdy to robisz. Po prostu nowe pliki są tam kopiowane. Nadpiszą pliki, które są już tam, jeśli mają tę samą nazwę pliku, ale nie wyczyszczą starych. Więc musisz być trochę ostrożny, gdy ponownie uruchamiasz pipeline'y. Jeśli nie chcesz, aby były na wierzchu plików, które już tam są. Upewnij się, że używasz pustego, czystego katalogu.

## 2.3. Ustaw tryb publikowania na kopiowanie

Dobrze, wspomniałem, że te pliki są soft linkami, więc jeśli zrobię _"ls -l results/hello_world/"_, możesz zobaczyć, że soft linkuje do katalogu roboczego. To jest na ogół dobra rzecz, jeśli pracujesz na czymś takim jak HPC, i to są naprawdę ogromne pliki i nie chcesz ich duplikować, ponieważ oznacza to, że pliki są przechowywane tylko raz w systemie plików.

Jednakże oznacza to, że jeśli usuniesz katalog roboczy: jeśli zrobię _"rm -r work"_ i wyczyszczę wszystkie te pliki pośrednie, które zostały utworzone. Teraz, jeśli spróbuję odczytać ten plik _"results/hello_world/"_. Będzie wskazywał jako soft link na plik, który już nie istnieje i dane są utracone na zawsze i nie można ich odzyskać, co może nie być świetne.

Więc ogólnie, powiem, że to dobra praktyka, aby kopiować pliki zamiast soft linkować, jeśli możesz, ponieważ jest to bezpieczniejsze. Po prostu bądź świadomy, że użyje to dwa razy więcej miejsca na dysku, chyba że usuniesz te katalogi robocze.

Aby to zrobić z blokiem wyjścia, przejdę do pierwszego wyjścia tutaj. Ustawiłem wcześniej ścieżkę i teraz ustawię tryb i możesz zobaczyć, gdy piszę, rozszerzenie VS code sugeruje rzeczy, wie, że jest to dyrektywa wyjściowa tutaj. I powiem copy. Naciskam zapisz.

Uruchommy ponownie workflow. Będzie tworzyć pliki ponownie, nowy katalog roboczy.

Teraz, jeśli przejdę do _"ls -l results/hello_world/"_ możesz zobaczyć, że to jest prawdziwy plik i nie jest już soft linkiem, i Nextflow skopiował to. Dobrze wiedzieć. Więc path i mode to rzeczy, które będziesz pisał dość często.

Teraz, oczywiście, to jest bardzo proste. Uczynimy to bardziej złożonym i potężnym w miarę postępów i zobaczysz, jak uczynić te rzeczy dynamicznymi i niezbyt rozwlekłymi.

## 2.4. Uwaga o dyrektywach publishDir na poziomie procesu

Teraz, powiedziałem na początku, że to jest dość nowa forma składni. Jest dostępna tylko w najnowszych wersjach Nextflow, gdy to nagrywam, i nazywa się Workflow Outputs.

Jeśli tego użyjesz, to świetnie. Odblokowuje to wiele innych fajnych funkcji w Nextflow, takich jak Nextflow Lineage, aby pomóc śledzić dziedzictwo tych plików, gdy są tworzone, i wkrótce będzie to domyślne w 26.04. A w późniejszym terminie w przyszłości będzie to jedyny sposób na pisanie Twoich workflow'ów.

Jednakże, gdy jesteśmy w tej fazie przejściowej teraz, możesz dobrze zobaczyć pipeline'y w naturze, które używają czegoś o nazwie publishDir, co jest starym sposobem na to, i jest to zdefiniowane nie na poziomie workflow i wyjścia, ale jest to zdefiniowane na poziomie procesu.

I ta deklaracja mówi zasadniczo to samo. Mówi, publikuj pliki wyników do katalogu o nazwie results i użyj trybu kopiowania. Więc możesz zobaczyć, że składnia jest bardzo podobna. Ale gdy piszesz nowe pipeline'y teraz, staraj się nie używać tej dyrektywy publishDir, nawet jeśli ją zobaczysz, w wynikach AI lub w dokumentacji lub innych pipeline'ach, ponieważ to jest stary sposób na to.

W 2026 wszyscy powinniśmy używać workflow outputs.

To wszystko jest udokumentowane, jeśli robisz to i używałeś Nextflow wcześniej, możesz przejść do dokumentacji Nextflow tutaj, nextflow.io/docs/. I jeśli przewinę w dół do tutorials, jest tutorial o nazwie _Migrating to Workflow Outputs_.

Jest naprawdę dobry. Przechodzi przez całą składnię, jak jest równoważna starej składni, dlaczego to zmieniliśmy i mają oś czasu i wszystko. I przechodzi przez wszystkie różne scenariusze z mnóstwem przykładów. Więc możesz łatwo przekonwertować istniejący kod Nextflow na nową składnię.

## 3.1. Zmień proces sayHello, aby oczekiwał zmiennego wejścia

Dobrze, więc mamy nasz prosty skrypt, który uruchamia proces, tworzy plik, mówi Nextflow, że to wyjście, a następnie mówimy Nextflow, gdzie zapisać ten plik. To dobry początek.

Ale byłoby ciekawiej, gdyby nie wszystko było zakodowane na stałe. Więc następnie pomyślmy o tym, jak powiedzieć Nextflow, że ten proces może przyjąć zmienne wejście, którym jest coś, co możemy kontrolować w czasie uruchomienia, gdy uruchamiamy workflow.

Musimy zrobić kilka różnych rzeczy, aby to się stało.

Po pierwsze, musimy powiedzieć temu procesowi, że może przyjąć zmienną wejściową i wpisujemy _input_ tutaj jako nowy blok deklaracji. I nazwiemy to _"val greeting"_.

Bit val jest odpowiednikiem path tutaj na dole. Mówi Nextflow, że to jest zmienna, jak string w tym przypadku. A jeśli na to najedziesz ponownie, powie Ci z rozszerzenia, co to znaczy.

Następnie powiemy Nextflow, co z tym zrobić. Nie wystarczy po prostu powiedzieć, że jest zmienna. Musisz powiedzieć w skrypcie, jak użyć tej zmiennej. Więc pozbędę się tego zakodowanego na stałe stringa tutaj i umieszczę zmienną.

Szybko zrobię to bez klamer, tylko żeby Ci pokazać, że jest to dozwolone, i to jest stary sposób na to. Ale teraz z nową składnią naprawdę zalecamy umieszczenie tego w klamrach tak, i sprawia to, że jest naprawdę jasne, że to jest interpolowane przez Nextflow tutaj.

Świetnie. Więc _"input greeting"_ wchodzi do _$\{greeting\}._ Ostatnia rzecz, którą musimy zrobić, to powiedzieć Nextflow na poziomie workflow, że ten proces teraz przyjmuje wejście. I aby to zrobić, zasadniczo damy mu zmienną.

## 3.2. Skonfiguruj parametr wiersza poleceń, aby przechwycić wejście użytkownika

Moglibyśmy zakodować to ponownie na stałe, jak Hello World, i to by dobrze zadziałało, ale oczywiście nie daje nam to naprawdę żadnej przewagi. Chcieliśmy być w stanie skonfigurować to w czasie uruchomienia, więc chcemy być w stanie zrobić to w CLI, gdy uruchamiasz Nextflow.

I sposób, w jaki to robimy, to specjalna koncepcja Nextflow o nazwie _params_. Nazwiemy to _params.input_.

Co to robi, to eksponuje tę zmienną wejściową w CLI i tam używamy podwójnego myślnika, gdy uruchamiamy Nextflow.

Mogę to nazwać, jak chcę, mogę to nazwać _hello, greeting_. Nie ma to znaczenia. Cokolwiek tam zrobię, będzie eksponowane jako opcja CLI, gdy uruchamiamy pipeline. I to jest prawdziwa sztuczka magiczna przez Nextflow, ponieważ oznacza to, że możesz zbudować swój skrypt workflow bardzo szybko z tymi parametrami i zasadniczo budujesz niestandardowe CLI dla swojego pipeline'u, sprawiając, że jest naprawdę łatwo dostosować różne opcje w locie, gdy uruchamiasz.

Więc. Spróbujmy tego. Wróćmy do naszego terminala. Mamy nasze polecenie _"nextflow run"_ tutaj. A teraz zrobię _"--input"_, co pasuje do _"params.input"_, które widzieliśmy wcześniej. Myślę, że w dokumentach jest po francusku. Geraldine lubi mówić po francusku. Zrobię to po szwedzku, bo mieszkam w Szwecji. więc powiem, "_Hej Världen_" i nacisnę enter.

Można użyć pojedynczych lub podwójnych cudzysłowów, to tylko wpływa na to, jak Bash to interpretuje.

Uruchamia pipeline Nextflow dokładnie w ten sam sposób. Możesz zobaczyć katalog roboczy i wszystko jest takie samo. Ale teraz jeśli przejdę do _"results/hello_world/output"_. Możemy zobaczyć nasz miły szwedzki tutaj zamiast tego.

Więc dynamicznie przekazaliśmy wejście z CLI do parametru. Przekazaliśmy to jako wejście do procesu i proces to zinterpretował i umieścił to w bloku skryptu, który następnie dynamicznie zmienił wyjście tego wyniku skryptu. Całkiem fajnie.

Dość złożona logika z bardzo małą składnią tutaj. I możesz miejmy nadzieję zobaczyć, jak to teraz zaczyna się skalować. I tak naprawdę budujemy logikę i konfigurowalność naszych pipeline'ów do skryptu Nextflow.

## 3.4. Użyj domyślnych wartości dla parametrów wiersza poleceń

Dobrze, to świetne. Problem jednak teraz jest taki, że za każdym razem, gdy uruchamiam ten pipeline, muszę zrobić dash, input, aby się uruchomił.

Jeśli spróbuję uruchomić bez tego parametru, teraz Nextflow rzuci błąd mówiąc, że potrzebował tego parametru i nie został ustawiony. i więc nie wiedział, co zrobić.

To jest fajna nowa rzecz, przy okazji. W przeszłości Nextflow po prostu uruchomiłby się z pustym stringiem i miałbyś wszelkiego rodzaju dziwne błędy, które byłyby trudne do zrozumienia. Ale w nowym parserze składni Nextflow jest trochę bardziej ostrożny i mówi Ci od razu.

Więc nie zawsze chcemy specyfikować każdą pojedynczą opcję. Dobrą praktyką jest określanie sensownych wartości domyślnych. Więc jak to robimy w naszym skrypcie?

Zauważysz, że gdy to napisaliśmy, po prostu umieściliśmy _params.input_ bezpośrednio tam, gdzie tego używamy. Więc oczywistym rozwiązaniem jest zdefiniowanie wartości domyślnej i robimy to na górze skryptu tutaj w specjalnym bloku params w workflow. To jest w skrypcie workflow tutaj.

Znowu, trochę nowej składni tutaj, więc zwróć uwagę. To naprawdę fajne rzeczy. Mamy nazwę parametru, który będzie tu oczekiwany.

A potem po tym dwukropku definiujemy typ zmiennej. Nie musisz tego robić, możesz po prostu zostawić to puste, ale jest naprawdę miłe. Mówi Nextflow, że oczekujemy stringa i traktujmy go jako taki.

Jeśli chcemy liczbę zamiast tego, na przykład, moglibyśmy napisać float i to powiedziałoby, że chcemy liczby zmiennoprzecinkowej. A jeśli spróbujemy uruchomić z tym, to rzuci błąd. Jeśli damy mu string, który nie jest float. A także przekaże go jako taki. Jeśli zrobimy string, wtedy wie, że to string. I nawet jeśli ma wiodące zera i jest całkowicie numeryczny, nadal przekaże go jako faktyczny string.

Więc to bezpieczeństwo typu jest bardzo nową funkcją Nextflow, ale naprawdę potężną, aby Twój kod był bezpieczniejszy do pisania i uruchamiania.

Następnie po tym mamy symbol równości, a następnie wartość domyślną tutaj. Nextflow został napisany w Barcelonie pierwotnie, więc wydaje się odpowiednie, że mamy trochę hiszpańskiego tutaj, _"Holà mundo!"_ jako wartość domyślną.

Dobrze, zapiszę ten skrypt, wrócę, uruchomię skrypt ponownie bez _--input_. I tym razem powinien się uruchomić i utworzy nasz nowy plik w _results_. I w tym pliku teraz mówi _"Holà mundo!"_.

To jest tylko wartość domyślna, więc nie oznacza to, że nadal nie możemy zrobić tego samego co wcześniej. Jeśli wrócę i znajdę mój stary skrypt tutaj, _"Hej Världen"_, ponieważ robię _--input_ w wierszu poleceń, to nadpisze tę wartość domyślną i użyje tego ponownie w pliku output.txt.

Więc to w skrypcie jest tylko wartość domyślna, którą ustawiam.

W miarę jak budujemy nasz workflow, aby był bardziej złożony i zawierał więcej parametrów, ten blok params na górze skryptu zacznie zbierać je wszystkie w jednym miejscu.

I kończysz z całkiem miłą symetrią w swoim skrypcie, gdzie faktycznie masz wszystkie swoje wejścia workflow tutaj i swoje wyjścia workflow na dole. I jest bardzo jasne, jaki jest interfejs Twojego workflow'u do świata zewnętrznego. Więc możesz bardzo szybko podjąć nowy pipeline z nową składnią i zrozumieć, jak go używać.

Jeszcze jedna fajna rzecz. Nie musimy ustawiać wartości domyślnej z tym. Jeśli zrobimy params input, ale nie ustawimy wartości domyślnej, to mówi Nextflow, że ten parametr jest wymagany, i znowu pipeline nie uruchomi się bez niego, ale da Ci bardziej użyteczny komunikat o błędzie zamiast czegoś o tym, że jest null.

Więc mówi, że oczekujemy jego input jest wymagany, ale nie został określony w wierszu poleceń. Bardzo ładnie.

Dobrze, więc miejmy nadzieję, że teraz jest jasne, jak skonfigurować Twój pipeline Nextflow ze zmiennymi wejściami i parametrami, jak ustawić wartość domyślną, ustawić typy, może to być Boolean true false flaga lub integer lub różne typy tutaj. Jak przekazać je do Twojego workflow, gdzie przechodzi, a następnie interpoluje do Twojego procesu. A także wiesz, jak dostosować te w wierszu poleceń, gdy uruchamiasz Nextflow. To zaczyna wyglądać ciekawiej niż nasze proste polecenie bash.

## 4. Zarządzaj wykonaniami workflow

Dobrze. Co dalej? Na końcową część tego rozdziału porozmawiamy trochę o tym, jak zarządzać wszystkimi różnymi wykonaniami workflow. Jeśli spojrzysz w moim pasku bocznym tutaj i Eksplorerze pod work, zobaczysz, że uruchomiłem kilka różnych pipeline'ów i te katalogi robocze stają się dość długie, jest ich dużo.

I druga rzecz jest taka, że jak powiedziałem wcześniej, za każdym razem, gdy ponownie uruchamiam ten pipeline, tworzy nowy zestaw katalogów roboczych i ponownie uruchamia wszystkie procesy od zera, co jest dobrą rzeczą. To jest zamierzone zachowanie. Jest powtarzalne i regeneruje wszystko na świeżo. Ale oczywiście, jeśli uruchamiasz bardzo długo trwające procesy, irytujące jest zawsze konieczność rozpoczynania pipeline'u od początku, jeśli zawiesił się w połowie lub jeśli zmienisz coś na końcu pipeline'u.

## 4.1. Ponownie uruchom workflow z -resume

Na szczęście Nextflow jest naprawdę dobry w wiedzy, co zostało wcześniej uruchomione i co jest dostępne, i aby ponownie użyć tych starych wyników jest bardzo proste. Po prostu dodajemy nową flagę na końcu polecenia _"-resume"_.

Teraz zauważ, że są dwa myślniki na input, bo to parametr. Jest tylko jeden myślnik na resume, ponieważ to jest podstawowa opcja Nextflow.

To podpycha ludzi cały czas, nawet jeśli używałeś Nextflow od dawna. Więc zawsze pamiętaj jeden lub dwa myślniki. Zależy, czy to jest podstawowa opcja Nextflow.

Dobrze, więc teraz robię _-resume_ i uruchamiam dokładnie ten sam workflow ponownie. I tym razem powinno wyglądać dość dokładnie tak samo z jedną kluczową różnicą.

W wyjściu tutaj możesz zobaczyć, że wyniki zostały zbuforowane. I w rzeczywistości ten hasz zadania tutaj jest dokładnie taki sam jak poprzednie uruchomienie i po prostu ponownie użył tego katalogu roboczego w całości. Wejścia i wyjścia i skrypt były wszystkie niezmodyfikowane. I więc po prostu bierze ten plik z tego i jeśli są kroki downstream w procesie, przekazałby je do następnego kroku w pipeline.

Więc nadal uruchamia cały pipeline od początku do końca, ale używa buforowanych wyników dla każdego z tych zadań, gdzie może.

Teraz, gdy robisz _-resume_, po prostu wznawia ostatni uruchomiony pipeline w Twoim katalogu roboczym, cokolwiek to było. Ale możesz faktycznie wznowić z dowolnego poprzedniego uruchomienia, które tam robiłeś. I zrobiliśmy już całkiem sporo teraz.

## 4.2. Zbadaj log przeszłych wykonań

Aby spojrzeć na wszystkie z nich, możemy zrobić _"nextflow log"_ zamiast _"nextflow run"_, i to da nam ładne wyjście pokazujące wszystkie te różne.. Muszę zmniejszyć mój ekran, żebyśmy mogli to zobaczyć, wszystkie te różne uruchomienia, gdy je zrobiliśmy, id sesji, polecenie i wszystko.

I możemy tu spojrzeć i możemy wziąć nazwę uruchomienia któregokolwiek z nich i następnie wznowić jeden z tych konkretnych. Więc mogę wrócić i mogę wznowić ten o nazwie _hungry_ekeblad_. I po prostu umieszczam to po _resume_.

Jeśli jesteś ciekawy, przy okazji, wszystkie te przymiotniki i nazwy naukowców są w kodzie źródłowym Nextflow. To naprawdę dobry sposób, aby dostać swój pierwszy pull request do Nextflow, idąc i znajdując to i dodając swojego ulubionego naukowca.

I tak czy inaczej, więc zrobiłem to i wróciło i spojrzało na zbuforowane wyniki z tego uruchomienia workflow, zdało sobie sprawę, że może je nadal ponownie użyć i zrobiło to. Więc dostałem zbuforowane wyniki ponownie.

## 4.3. Usuń starsze katalogi robocze

To świetne. Co jeśli chcę wyczyścić te katalogi robocze? Jest ich tutaj mnóstwo. Jest mnóstwo plików. Może wiem na pewno, że chcę wznowić z ostatnich kilku uruchomień pipeline, ale nie obchodzą mnie wszystkie te sprzed tego.

Wtedy mogę wybrać jedno tutaj i mogę użyć innego polecenia Nextflow, którym jest _"nextflow clean"_, i mogę zrobić _"nextflow clean"_, zrobię _"-before"_, i konkretną nazwę uruchomienia, która w tym przypadku była _reverent_pike_ i zrobię _"-n"_, co mówi Nextflow, aby po prostu zrobił dry run. Więc po prostu mówi mi, co usunie. Bez faktycznego robienia czegokolwiek, więc usunąłby te katalogi robocze.

To wygląda sensownie. Więc zrobię to samo polecenie ponownie, ale zamiast _"-n"_ zrobię _"-f"_, aby faktycznie wykonać czyszczenie. I tym razem faktycznie usunął wszystkie te katalogi. A jeśli wejdę i spojrzę na katalogi robocze, wygląda teraz znacznie lżej. Fantastycznie.

Więc tak czyścisz wszystkie swoje lokalne katalogi robocze w całkiem bezpieczny sposób bez kompletnego niszczenia cache. Więc nadal możesz wznowić, jeśli chcesz.

Jeśli kiedykolwiek zapomnisz, jakie są te flagi dla każdego polecenia Nextflow, możesz zrobić _"nextflow help"_, a następnie nazwę polecenia. Więc jeśli zrobię _"nextflow help clean"_, możesz zobaczyć wszystkie różne opcje: _-after, -before, -but_, wszystkie różne sposoby konfiguracji tego zachowania czyszczenia. Całkiem fajnie.

## Podsumowanie

Dobrze, to koniec części pierwszej Hello Nextflow. To całkiem intensywny początek kursu, ale miejmy nadzieję, że masz teraz całkiem dobre zrozumienie, jak wygląda skrypt Nextflow; z różnymi kluczowymi częściami, procesami, workflow'ami, wyjściami i parametrami. Wiesz, jak je konfigurować z podstawowymi nadpisaniami z wiersza poleceń, jak zrobić dynamiczny blok wejściowy z dynamicznym skryptem i wiesz, jak zarządzać wszystkimi swoimi wykonaniami obciążenia: widząc, co już uruchomiłeś, wznawiając, czyszcząc. Jest dużo rzeczy. Przeszedłeś długą drogę. Więc jeśli chcesz zrobić przerwę i mieć krótki spacer i herbatę, teraz jest prawdopodobnie dobry czas. Zasłużyłeś na to.

Od teraz zasadniczo budujemy na tym fundamencie. Jak możemy to uczynić bardziej złożonym, bardziej potężnym? Jak możemy to uczynić bardziej elastycznym? Robić rzeczy, które chcemy zrobić naszą analizę na dużą skalę.

## Quiz

Teraz jeśli przewiniesz w dół do części pierwszej, hello world, na stronie internetowej zobaczysz mały quiz i to jest coś nowego, co zrobiliśmy dla tej wersji szkolenia Nextflow. I możesz przejść i sprawdzić się, aby sprawdzić, czy zrozumiałeś wszystkie materiały, które zrobiliśmy w tym rozdziale.

Nie jest to do nas wysyłane ani nic, po prostu jest przechowywane w Twojej przeglądarce. Więc nie wiemy, jakie są Twoje odpowiedzi, ale to tylko mała samodzielna kontrola, aby upewnić się, że nie przeoczyłeś niczego lub czegoś nie zrozumiałeś. I możesz spróbować tyle razy, ile chcesz.

Jeśli jesteś jak ja, może chcesz pozostać w terminalu w swoim instancji VS Code, w którym to przypadku możesz wpisać polecenie _quiz_, a następnie po prostu powiedzieć mu, w którym rozdziale jesteś. Więc robimy _"Hello World"_, a następnie możesz zrobić dokładnie ten sam quiz, pytania, które są w przeglądarce internetowej, ale tylko w swoim terminalu.

Fajnie. Dobrze. Mam nadzieję, że Ci się to podoba. Pobaw się trochę i zobaczymy się w następnym rozdziale za chwilę, aby porozmawiać o wszystkim o kanałach Nextflow.
