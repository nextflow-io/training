# Część 6: Hello Config - Transkrypcja

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Tłumaczenie wspomagane przez AI - [dowiedz się więcej i zasugeruj ulepszenia](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

<div class="video-wrapper">
  <iframe width="560" height="315" src="https://www.youtube.com/embed/IuDO2HeKvXk?si=tnXTi6mRkITY0zW_&amp;list=PLPZ8WHdZGxmXiHf8B26oB_fTfoKQdhlik" title="YouTube video player" frameborder="0" allow="accelerometer; autoplay; clipboard-write; encrypted-media; gyroscope; picture-in-picture; web-share" referrerpolicy="strict-origin-when-cross-origin" allowfullscreen></iframe>
</div>

!!!note "Ważne informacje"

    Ta strona zawiera jedynie transkrypcję. Aby uzyskać pełne instrukcje krok po kroku, wróć do [materiałów szkoleniowych](../06_hello_config.md).

    Numery sekcji pokazane w transkrypcji są podane jedynie w celach orientacyjnych i mogą nie obejmować wszystkich numerów sekcji w materiałach.

## Powitanie

Cześć, witajcie w szóstej części kursu szkoleniowego Hello Nextflow.

Ten rozdział nazywa się Hello Config i jest ostatnią częścią naszego kursu szkoleniowego.

W tym rozdziale będziemy omawiać konfigurację Nextflow'a. Konfiguracja Nextflow'a jest naprawdę potężna. Pozwala nam uruchamiać ten sam pipeline na wielu różnych infrastrukturach obliczeniowych, z różnym dostarczaniem oprogramowania i różnymi opcjami w samym pipeline.

Oznacza to, że możecie wziąć pipeline'y Nextflow'a zbudowane przez inne osoby i uruchomić je na Swoim systemie, nawet jeśli mogły być zbudowane dla zupełnie innej infrastruktury. Ta możliwość konfigurowania Nextflow'a sprawia, że workflow'y są naprawdę przenośne i możliwe do udostępniania.

W tym rozdziale będziemy używać workflow'a, który zbudowaliśmy we wcześniejszych częściach, ale nie będziemy w ogóle edytować jego kodu. Przyjrzymy się tylko naszemu plikowi konfiguracyjnemu Nextflow'a i zobaczymy, jak zmiana konfiguracji zmienia sposób działania Nextflow'a.

Dobrze, zaczynajmy.

Tak jak wcześniej, zacznijmy od przejścia do training.nextflow.io. Po lewej stronie przejdźcie do Hello Nextflow i rozdziału szóstego - Hello config. Teraz przejdę do mojego środowiska GitHub Codespaces i sprawdzę skrypt, którego będziemy używać.

## 0. Rozgrzewka: Sprawdź, czy Docker jest włączony i uruchom workflow Hello Config

Ten nazywa się Hello Config i zaczyna od tego miejsca, w którym zakończyliśmy wcześniej. Wygląda więc dokładnie tak samo - z trzema parametrami: Greetings dla pliku CSV, batch dla nazwy partii wyjściowej i character dla nazwy cowpy. Mamy też nasze cztery importy różnych procesów, a następnie workflow, w którym je łączymy.

Faktycznie zamknę teraz ten plik, ponieważ w tym rozdziale nie będziemy w ogóle dotykać pliku Nextflow'a. Będziemy pracować wyłącznie w pliku konfiguracyjnym. Jeśli spojrzę do pliku Nextflow dot config, który krótko omówiliśmy w poprzednim rozdziale piątym, możemy zobaczyć, że mamy tutaj pojedynczą instrukcję: Docker enabled równa się true, która mówi Nextflow'owi, aby używał Dockera podczas wykonywania tego workflow'a.

Używam pliku Nextflow dot config w katalogu głównym pipeline'u, który jest ładowany automatycznie, gdy uruchamiam Nextflow'a. Ale pamiętajcie, Nextflow może ładować pliki konfiguracyjne z wielu miejsc.

Jeśli sprawdzę w dokumentacji Nextflow'a - przejdę do Configuration - możecie zobaczyć listę tych miejsc oraz priorytet, w jakim są ładowane.

Dobrze. Sprawdźmy, czy nasz workflow wykonuje się tak, jak oczekujemy. Otworzę terminal. Robię Nextflow. Run. Hello, config. I naciskam enter. Powinniśmy zobaczyć te cztery procesy uruchomione, kończące się poleceniem cowpy. Rzeczywiście, to zadziałało prawidłowo. Miałem włączonego Dockera, pobrał obraz Dockera i uruchomił dla mnie cowpy, tak jak na końcu rozdziału piątego.

## 1. Określ, jakiej technologii pakowania oprogramowania użyć

Dobrze. Powiedzmy, że pracuję na HPC i nie mam zainstalowanego Dockera. Najlepszą rzeczą do zrobienia w tym scenariuszu byłoby użycie Singularity lub Apptainer. Gdybym miał to zrobić, wszedłbym do modułu cowpy i zmieniłbym ten kontener, aby używał obrazu singularity, jak pokazałem w poprzednim rozdziale, z oras://, który również można uzyskać z Seqera Containers.

Następnie przeszedłbym do pliku Nextflow dot config, ustawił Docker enabled na false i zrobił singularity enabled równa się true. Lub, jeśli używam Apptainer, apptainer enabled równa się true - i to by zadziałało.

Nextflow obsługuje również inne technologie oprócz kontenerów - coś, co możecie znać, to conda. Tutaj możemy zrobić conda enabled równa się true i ustawić Docker na false. conda nie używa tej samej dyrektywy container. Zamiast tego możemy dodać tutaj nową dyrektywę o nazwie conda. Następnie określamy pakiet conda, którego chcemy użyć. Dobrą praktyką jest bycie tak szczegółowym, jak to możliwe, aby spróbować uczynić pipeline jak najbardziej powtarzalnym. Więc zamierzam określić kanał conda - conda-forge - a następnie cowpy i dokładną wersję, którą była 1.1.5.

Mógłbym też po prostu napisać cowpy, gdybym chciał, ale mogłoby to dać inną wersję cowpy przy różnych wykonaniach pipeline'u.

Fajną rzeczą w tym jest to, że nie dotknąłem w ogóle dyrektywy docker. Ten obraz Dockera jest nadal tam. Po prostu dostarczam teraz dwie alternatywy, a te mogą być włączane lub wyłączane tylko za pomocą pliku konfiguracyjnego.

## 1.3. Uruchom workflow, aby zweryfikować, że może używać Conda

Conda jest teraz włączona, więc spróbujmy.

Świetnie. Działa i możecie zobaczyć, że jest tu komunikat od Nextflow'a mówiący, że Nextflow tworzy dla mnie środowisko conda i używa tej lokalizacji cache.

W tle Nextflow uruchamia dla mnie polecenia "conda create", aby stworzyć nowe izolowane środowisko conda z tylko tymi pakietami, których chcę, a następnie instaluje i pobiera te pakiety conda, aby mogło uruchomić proces.

Możecie zobaczyć, że zajęło to trochę czasu, ponieważ tworzyło środowisko i instalowało oprogramowanie po raz pierwszy. Jednak buforuje to środowisko, więc jeśli uruchomię to samo polecenie Nextflow ponownie, powinno być znacznie szybsze, ponieważ będzie ponownie używać tego samego środowiska conda.

Jedną z fajnych rzeczy w tym jest to, że te dyrektywy mogą być określone na poziomie procesu, a nie tylko całego workflow'a. Więc jeśli chcecie, możecie mieszać i dopasowywać, jaka technologia jest używana dla różnych procesów.

## 2. Przydziel zasoby obliczeniowe za pomocą dyrektyw procesu

Plik konfiguracyjny Nextflow'a może robić znacznie więcej niż tylko pakowanie oprogramowania. Możemy również powiedzieć Nextflow'owi, jak faktycznie uruchamiać kroki w pipeline. Jednym z przykładów jest powiedzenie systemowi hostującemu, jakie zasoby powinny być udostępnione każdemu wykonywanemu zadaniu.

Domyślnie Nextflow nie daje zbyt wiele. Daje pojedynczy procesor CPU i tylko dwa gigabajty pamięci każdemu procesowi.

To jest prawdopodobnie coś, co chcielibyśmy zmienić, aby procesy, które zajmują dużo czasu, mogły mieć więcej zasobów i działać szybciej. Ale może być trudno wiedzieć, co przydzielić procesowi. Nextflow ma kilka fajnych sztuczek w rękawie, aby Wam w tym pomóc.

## 2.1. Uruchom workflow, aby wygenerować raport wykorzystania zasobów

Uruchommy workflow ponownie. Tym razem zamierzam dodać dodatkowy argument, którym jest dash with reports. To jest podstawowa opcja Nextflow'a, więc jest to pojedynczy łącznik. A następnie dowolna nazwa pliku. W tym przypadku zamierzam nazwać to report config one html.

Zamierzam uruchomić workflow ponownie. Będzie działać dokładnie tak samo jak wcześniej, ale da mi dodatkowy raport pomocniczy, który możecie zobaczyć - pojawił się teraz tutaj na pasku bocznym.

Kliknę prawym przyciskiem myszy na ten plik, kliknę download, co pobierze go z GitHub Codespaces do mojego lokalnego systemu, abym mógł łatwo obejrzeć go w przeglądarce internetowej tutaj.

Ten raport może być wygenerowany dla każdego uruchomienia Nextflow'a i zawiera wiele informacji. Zaczyna się na górze od metadanych o tym, jakie polecenie zostało użyte, kiedy workflow został uruchomiony, ile czasu to zajęło, ale gdy przewijamy w dół, otrzymujemy bardziej szczegółowe informacje o zasobach, które zostały użyte przez każdy krok w pipeline.

Ponieważ każdy proces uruchamia się wiele razy dla różnych zadań, mamy wykres pudełkowy pokazujący zmienność zasobów, których użyliśmy dla każdego procesu.

Jeśli przewinę trochę dalej w dół, zobaczę podobne informacje o użytej pamięci i czasie trwania zadania. Także odczyt i zapis dysku.

Możecie sobie wyobrazić, że dla dużego pipeline'u z długo działającymi zadaniami może to być bardzo pouczające o tym, jak dostroić konfigurację zasobów, o które prosimy, aby nie żądać zbyt wiele, ale także aby zapewnić wystarczająco dużo, że działa szybko.

Jeśli będę dalej przewijać raport, zobaczymy również tabelę zadań, która pokazuje nam szczegółowe informacje o każdym pojedynczym zadaniu, które zostało uruchomione w workflow. Obejmuje to informacje takie jak rozwiązany skrypt, który został uruchomiony.

Dobrze, wróćmy do naszego pliku konfiguracyjnego. Zobaczyliśmy, że naprawdę nie potrzebowaliśmy zbyt wiele dla naszego workflow'a, więc powiedzmy Nextflow'owi, że potrzebujemy tylko jednego gigabajta pamięci dla każdego procesu w workflow.

Teraz, gdy definiujemy to w ten sposób na poziomie procesu, jest to stosowane do każdego pojedynczego procesu w pipeline.

## 2.3. Ustaw alokacje zasobów dla pojedynczego procesu

Dla przykładu, udajmy, że cowpy naprawdę wykonuje dużo ciężkiej pracy i potrzebuje więcej zasobów niż inne zadania. Możemy zdefiniować dodatkowy blok konfiguracji tutaj, który stosuje się tylko do tego procesu, używając with name cowpy.

To jest nazywane selektorem konfiguracji i możemy tutaj zdefiniować różne wzorce, aby dopasować różne procesy. Na przykład, mógłbym zrobić cow star. Następnie dostaję nawiasy klamrowe i dajmy mu dwa gigabajty pamięci zamiast jednego i powiedzmy dwa procesory CPU.

Teraz Nextflow będzie dawać każdemu procesowi w workflow jeden gigabajt oprócz tego żądania, które jest bardziej szczegółowe. Więc to je nadpisuje. I tylko dla procesów, które nazywają się cowpy, otrzymają dwa gigabajty pamięci i dwa procesory CPU.

Zauważcie, że Nextflow jest mądry w kwestii wykorzystania zasobów. Więc jeśli zaczniesz ustawiać te liczby na wyższych wartościach, zobaczysz, że Nextflow zaczyna kolejkować przesyłanie zadań jedno po drugim, zamiast uruchamiać je wszystkie równolegle, aby nie nadmiernie żądać dostępnych zasobów.

## 2.4. Uruchom workflow ze zmodyfikowaną konfiguracją

Spróbujmy uruchomić workflow ponownie i zapiszmy nowy raport tym razem.

Dobrze, możemy pobrać ten plik i rzucić okiem.

Tak, nic dziwnego, wygląda to zasadniczo dokładnie tak samo, ponieważ to jest fikcyjny workflow, który nic prawdziwego nie robi. Ale możecie sobie wyobrazić, jak to iteracyjne podejście do definiowania limitów i wykonywania prawdziwych workflow'ów z tego rodzaju raportowaniem pozwala na podejście oparte na dowodach do ustalania odpowiedniej konfiguracji i naprawdę maksymalne wykorzystanie dostępnych dla Was zasobów obliczeniowych.

Możecie zacząć być naprawdę sprytni w tej kwestii. Nextflow ma wbudowaną zdolność do ponawiania prób w przypadku awarii i możecie to wykorzystać w Swoim pliku konfiguracyjnym, używając domknięcia takiego jak to i dynamicznie ustawiając dostępne zasoby. Więc tutaj powiedziałem Nextflow'owi, aby pomnożył te dwa gigabajty przez próbę ponowienia. Więc druga próba otrzyma cztery gigabajty, trzecia próba otrzyma sześć gigabajtów i tak dalej. To trochę wykracza poza zakres tego kursu szkoleniowego, ale jeśli jesteście zainteresowani, sprawdźcie dokumentację Nextflow'a, która ma fajną sekcję o dynamicznej logice ponownych prób.

## 2.5. Dodaj limity zasobów

Teraz, jedną rzeczą, którą możecie zauważyć, jest to, że tego rodzaju rzecz może ułatwić przypadkowe wyjście poza zasoby dostępne w Waszym systemie. Jeśli zażądacie więcej zasobów niż jest dostępnych, Nextflow wyrzuci błąd dotyczący Waszej konfiguracji i zatrzyma uruchomienie. Aby tego uniknąć, możecie użyć czegoś, co nazywa się limitami zasobów.

Pod zakresem process w naszym workflow możemy zdefiniować limity zasobów w ten sposób, co przyjmuje tablicę, i możemy określić maksymalną pamięć, procesory CPU i czas, które są dostępne w tym systemie.

Ustawienie tutaj wysokich wartości nie zwiększa ilości żądanych zasobów. Nadal będziemy używać jednego gigabajta w naszych żądaniach, ale oznacza to, że jeśli którekolwiek z tych żądań osiągnie 750, osiągną ten pułap i nic więcej nie będzie żądane, co oznacza, że Nextflow będzie kontynuować działanie i nie zawiesi się z powodu niedostępnych zasobów.

Więc to jest fajne zabezpieczenie do użycia, szczególnie jeśli używacie dynamicznej logiki z Waszą alokacją zasobów.

Druga sytuacja, w której jest to naprawdę użyteczne, jest wtedy, gdy używacie pipeline'ów, które są publiczne i nie są kontrolowane przez Was. Mogą pochodzić z domyślnymi konfiguracjami, a Nextflow automatycznie przyjmie właściwe podejście do ograniczenia wszelkich żądań zasobów, aby działać na Waszym systemie.

Dobrze, świetnie. Mówiliśmy o oprogramowaniu. Mówiliśmy o alokacji zasobów i opisaliśmy różne zakresy konfiguracji, zarówno dla wszystkich procesów, jak i konkretnych procesów.

## 3. Użyj pliku parametrów do przechowywania parametrów workflow

Dobrze, następnie zwrócimy naszą uwagę na parametry. Możemy zdefiniować parametry w pliku konfiguracyjnym tak samo, jak zrobiliśmy to wcześniej w skrypcie Nextflow'a. Więc params dot greeting równa się hello lub użyj zakresu params i ustaw foo równa się bar.

I to świetnie nadaje się do ustawiania wartości domyślnych dla Waszego workflow'a. Jednak podczas uruchamiania pipeline'ów może być miło określić parametry w pliku JSON lub YAML.

Używanie takiego pliku jest znacznie lepsze niż określanie opcji linii poleceń z dash dash. Ponieważ podczas uruchamiania workflow'a możecie musieć określić wiele parametrów i może być żmudne, aby je wszystkie napisać w jednym CLI - a także podatne na błędy. Ponadto jest mało prawdopodobne, że zapamiętacie wszystkie parametry, których użyliście, więc jeśli zakodujecie to w pliku, łatwiej jest uruchomić workflow ponownie, używając tych samych parametrów w przyszłości.

Mamy tutaj przykładowy plik o nazwie test params i możecie zobaczyć, że określa to trzy parametry, które mamy w naszym workflow z trzema różnymi wartościami. Osobiście uważam, że YAML jest łatwiejszy do pisania niż JSON. Więc tylko po to, aby zademonstrować, że to działa, zamierzam utworzyć nowy plik o nazwie Test yaml i skopiować je, pozbyć się cudzysłowów. I nacisnąć zapisz.

Te pliki JSON i YAML mogą być łatwiejsze do pisania, ponieważ są bardziej znaną składnią. Ale zauważcie, że są one tylko dla parametrów i przyjmują tylko składnię klucz-wartość taką jak ta.

## 3.1. Uruchom workflow używając pliku parametrów

Spróbujmy. Zrobię to samo polecenie co wcześniej. Pozbędę się raportu i zamierzam zrobić dash params file test params yaml.

Nie, to jest podstawowa opcja Nextflow'a, więc jest to pojedynczy łącznik.

Dobrze. Uruchomiło workflow i użyło parametrów z tego pliku YAML zamiast określania ich wszystkich w linii poleceń. Może wydawać się to przesadą tylko dla tego prostego przykładu, ale możecie sobie wyobrazić, że jeśli macie 10 lub 20 różnych parametrów, może być uciążliwe wpisywanie ich ręcznie, a to jest po prostu znacznie łatwiejsze do edycji w edytorze kodu i przechowywania ze względu na powtarzalność.

## 3. Określ, jaki(e) executor(y) powinny być użyte do wykonania pracy

Dobrze. Mówiliśmy o pakowaniu oprogramowania z Docker i conda. Mówiliśmy o wymaganiach zasobów procesu z procesorami CPU i pamięcią. I mówiliśmy trochę o tym, jak określić parametry podczas uruchamiania workflow'ów.

Końcowe części konfiguracji to naprawdę wykonanie - sama podstawowa infrastruktura obliczeniowa - i to jest prawdziwy klejnot w koronie Nextflow'a: to, że możemy uruchamiać te same workflow'y na wielu różnych infrastrukturach obliczeniowych.

Faktycznie zamierzam przełączyć się na pisemny materiał szkoleniowy przez sekundę. W tej części szkolenia możemy zobaczyć kilka różnych przykładów tego, jak różne executors - w tym przypadku harmonogramy HPC - definiują wymagania zasobów potrzebne do przesłania zadania.

Więc dla Slurm macie te nagłówki SBATCH, które definiują dash dash mem i numer CPU. Jeśli używacie PBS, macie różne nagłówki, a jeśli używacie Grid Engine, macie znowu różne nagłówki.

Możecie sobie wyobrazić, że jest jeszcze bardziej różnie, jeśli chcecie uruchomić w chmurze - czy to AWS batch, Google Cloud, Azure czy więcej.

Każda z tych podstawowych infrastruktur obliczeniowych nazywana jest executorem i Nextflow wie, jak rozmawiać ze wszystkimi tymi różnymi executorami, aby przesłać zadania z poprawną składnią.

Dobra wiadomość jest taka, że nie musicie o tym wiedzieć. Wszystko, co musicie zrobić, to powiedzieć Nextflow'owi, którego executora użyć.

## 3.1. Targetowanie innego backendu

Wracamy do naszego pliku konfiguracyjnego i w procesie robimy executor, i zamierzam wpisać local.

Local jest faktycznie domyślny - jeśli nie określicie żadnego innego executora, local jest tym, co będzie używane, a to po prostu oznacza Wasz system hostujący, gdziekolwiek uruchomiliście Nextflow'a.

Mógłbym określić zamiast tego Slurm. I to przesłałoby zadania Slurm, lub mógłbym powiedzieć AWS batch, i to przesłałoby zadania do AWS batch.

W niektórych przypadkach potrzebujecie dodatkowej konfiguracji - na przykład uruchamianie w chmurze będzie wymagać pewnych poświadczeń - ale naprawdę to jest rdzeń tego, i może to być tak proste jak jedna lub dwie linie konfiguracji, aby uruchomić Wasz workflow w zupełnie innym środowisku obliczeniowym.

Mimo że działamy na prostym systemie w Codespaces, nadal mogę trochę się z tym pobawić i udawać, że działamy na Slurm. Jeśli następnie uruchomię workflow ponownie - Nextflow run, hello config - to się nie powiedzie, ponieważ nie będzie mogło przesłać zadań do Slurm. Ale nadal możemy wejść do katalogów work i zobaczyć, co Nextflow zrobił. Więc jeśli przejdziemy do tego katalogu work i spojrzymy na Command Run, możecie zobaczyć na górze tego pliku, że mamy teraz te linie nagłówka sbatch, które próbowały określić zasoby potrzebne dla zadania Slurm.

## 4. Użyj profili do wyboru wstępnie ustawionych konfiguracji

Dobrze, prawie jesteśmy na miejscu. Ostatnia część tego rozdziału to mówienie o profilach konfiguracyjnych. Jeśli uruchamiacie Swój pipeline na kilku różnych systemach, mogłoby być irytujące mieć wszystkie te różne pliki konfiguracyjne Nextflow'a, które musicie określać za każdym razem.

Zamiast tego możecie zakodować grupowania konfiguracji w Waszym pliku konfiguracyjnym Nextflow'a i włączać oraz wyłączać te grupy, używając flagi profile. Zobaczmy, jak to wygląda.

## 4.1. Utwórz profile do przełączania między lokalnym rozwojem a wykonaniem na HPC

Zamierzamy utworzyć dwa profile w naszym przykładzie tutaj - jeden dla mojego laptopa i jeden dla cięższego systemu HPC. Zamierzam trochę oszukać i po prostu skopiować kod z materiału szkoleniowego i umieścić go tutaj.

Mamy nowy zakres o nazwie profiles, a następnie mamy nazwę dla każdego profilu, która może być dowolna. I w tym mamy konfigurację, która wygląda dokładnie tak samo jak konfiguracja najwyższego poziomu, którą już napisaliśmy. Więc znowu mamy zakres process. Zakres Docker.

W profilu o nazwie my laptop mówię, aby uruchomić używając executora local, więc na moim systemie hostującym, i aby używać Dockera.

W profilu university HPC tutaj mówię, aby używać Slurm do przesyłania zadań, aby używać conda zamiast Dockera, i określam różne limity zasobów, które mogą pasować do rozmiaru systemu węzłów na HPC, którego używam.

Domyślnie żadna z tej konfiguracji nie będzie używana, gdy uruchomię Nextflow'a - muszę określić, że chcę użyć jednego z tych profili.

## 4.2. Uruchom workflow z profilem

Zróbmy nextflow run hello config. I zamierzam zrobić dash profile, pojedynczy łącznik, ponieważ to jest podstawowa opcja Nextflow'a. A następnie nazwa, którą mu nadałem, którą jest my laptop. Nextflow powinien teraz użyć bloku konfiguracji, który został określony w tym profilu konfiguracyjnym i zastosować go podczas uruchamiania Nextflow'a. Gdybym chciał użyć drugiego bloku konfiguracji, musiałbym tylko przełączyć tę nazwę profilu. Znacznie łatwiejsze do zapamiętania. Znacznie łatwiejsze w użyciu.

## 4.3. Utwórz profil testowy

Zauważcie, profile mogą mieć dowolny rodzaj konfiguracji, więc nie musi to być związane z Waszym środowiskiem wykonawczym. Na przykład, stwórzmy tutaj nowy profil, który ma zestaw parametrów. Możemy zmienić to na tux i zmienić na my profile, a teraz gdy zrobimy profile test, to określi te parametry, które nadpiszą parametry, które są określone na najwyższym poziomie workflow'a.

Kiedy uruchamiacie Nextflow'a, możecie połączyć wiele profili i będą one stosowane kolejno.

## 4.4. Uruchom workflow lokalnie z profilem testowym

Więc mogę wziąć poprzednie polecenie i zrobić przecinek test. To zastosuje konfigurację my laptop najpierw, a następnie zastosuje konfigurację test. Jeśli jest jakiekolwiek nakładanie się, to profil po prawej nadpisze wszelką konfigurację w poprzednich profilach. Jeśli nacisnę enter, zobaczmy, co się stanie.

Dobrze, mamy tutaj nowy plik wyników. Możecie zobaczyć My Profile, który określiłem jako jedną z opcji. I możemy również zobaczyć cowpy, my profile, i rzeczywiście, jest tux. Więc to zadziałało.

## Podsumowanie

Dobrze! Niesamowite. To wszystko. Dotarliście do końca kursu. Dostajecie trochę konfetti celebracyjnego. Brawo za ukończenie tego rozdziału.

[Transkrypcja następnego filmu :octicons-arrow-right-24:](07_next_steps.md)
