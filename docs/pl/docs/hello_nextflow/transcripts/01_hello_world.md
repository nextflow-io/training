# Część 1: Hello World - Transkrypcja

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Tłumaczenie wspomagane przez AI - [dowiedz się więcej i zasugeruj ulepszenia](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

<div class="video-wrapper">
  <iframe width="560" height="315" src="https://www.youtube.com/embed/8X2hHI-9vms?si=F0t9LFYLjAWoyRXj&amp;list=PLPZ8WHdZGxmXiHf8B26oB_fTfoKQdhlik" title="YouTube video player" frameborder="0" allow="accelerometer; autoplay; clipboard-write; encrypted-media; gyroscope; picture-in-picture; web-share" referrerpolicy="strict-origin-when-cross-origin" allowfullscreen></iframe>
</div>

!!!note "Ważne uwagi"

    Ta strona zawiera tylko transkrypcję. Pełne instrukcje krok po kroku znajdziesz w [materiałach szkoleniowych](../01_hello_world.md).

    Numery sekcji w transkrypcji są podane wyłącznie w celach informacyjnych i mogą nie obejmować wszystkich numerów sekcji z materiałów.

## Powitanie

Cześć, witamy w pierwszym rozdziale Hello Nextflow.

W tej pierwszej części sześcioczęściowego kursu przejdziemy przez podstawy Nextflow'a. Zaczniemy od uruchomienia kilku poleceń w terminalu, a następnie zobaczymy, jak przekształcić te polecenia Bash w skrypt Nextflow'a.

Spróbujemy uruchomić pierwszy pipeline Nextflow'a, zobaczymy, co Nextflow robi, gdzie działa, jakie pliki tworzy i jaki jest ich cel.

Dobra, zaczynajmy.

## training.nextflow.io

Przede wszystkim wejdź na training.nextflow.io. Tak jak wcześniej, wszystkie materiały są tu dostępne. Będę pracował przez nie krok po kroku. Będę pokazywał mój ekran podczas wykonywania kroków szkolenia, ale wszystko, co mówię, znajdujesz w materiałach, więc możesz je śledzić we własnym tempie.

To wideo ma również włączone napisy, więc możesz je włączyć i śledzić dokładnie to, co mówię.

Dobra, przejdźmy do Hello Nextflow. To kurs, który dzisiaj będziemy realizować. Orientację już zrobiliśmy w pierwszym wideo, więc przejdziemy od razu do części pierwszej: Hello World.

Teraz opuszczę materiały i przejdę do mojego środowiska Code Spaces. To jest to, co skonfigurowaliśmy w pierwszym wideo. Mam nadzieję, że Twoje środowisko wygląda podobnie. Używam VS Code, patrzę na materiały szkoleniowe i zmieniłem katalog na hello Nextflow.

## 0. Rozgrzewka: Uruchom Hello World bezpośrednio

Dobra. Zacznijmy od kilku podstaw, które mam nadzieję będą wszystkim znajome. Zacznę od napisania bardzo prostego polecenia w terminalu. Tutaj na dole napiszę 'echo Hello World!"', nacisnę enter i, bez niespodzianek, terminal robi to, czego od niego żądam, i zwraca ten ciąg znaków: Hello world.

Dobra, teraz nacisnę strzałkę w górę, żeby przywołać to polecenie, i trochę je zmodyfikuję. Tym razem przekierujmy to wyjście do pliku. Zapiszę je do output.txt, nacisnę enter — nic tym razem w terminalu, ponieważ wyjście nie trafiło do terminala. Trafiło do tego pliku.

Mogę wtedy odczytać ten plik, robiąc 'cat output.txt', nacisnąć tab, żeby automatycznie rozszerzyć nazwę pliku — i proszę. Plik jest tam.

Widzę również ten plik na pasku bocznym w eksploratorze plików w VS Code. Mogę go dwukrotnie kliknąć i otworzyć tutaj. Jeśli chcesz otworzyć go w VS Code bez klikania, możesz również zrobić "code", a następnie "output.txt" — robi to samo.

Świetnie. To pierwszy krok. Bardzo prosty.

## 1. Przeanalizuj startowy skrypt workflow'u Hello World

Dobra. Teraz zrobimy dokładnie to samo, ale w Nextflow, zamiast bezpośrednio w terminalu.

Użyjemy pierwszego przykładowego skryptu na początek. Ten plik nazywa się Hello World. Mogę zrobić "ls", żeby go zobaczyć w terminalu, a jestem na Macu, więc mogę zrobić command+click, żeby otworzyć ten plik, lub mogłem po prostu dwukrotnie kliknąć na pasku bocznym tutaj.

Jest kilka rzeczy, które możemy zobaczyć w tym pliku. Na samej górze jest instrukcja hash mówiąca, że to jest plik Nextflow'a i tak mógłby być wykonany. Są tu jakieś komentarze, zwykłe komentarze kodu w jasnoszarym kolorze, które nie wpływają na wykonanie — tylko pomagają nam czytać skrypt.

A potem są dwie główne struktury: `process` i `workflow`.

Procesy w Nextflow to kroki pipeline'u. To części, które faktycznie wykonują logikę i przetwarzanie.

Workflow na dole łączy te procesy razem i kontroluje logikę workflow'u — jak wszystko się ze sobą łączy.

Zaczniemy od spojrzenia na `process`. Za chwilę wrócimy do `workflow`.

## 1.2 Definicja procesu

Więc każdy `process` zaczyna się od słowa kluczowego `process`. Ma nazwę, a następnie nawiasy klamrowe — wszystko w tych nawiasach to ten pojedynczy proces.

Proces musi mieć sekcję `script`, a zawarta tutaj jest fragmentu bash w wielowierszowym ciągu znaków. To część kodu faktycznie wykonywana w środowisku obliczeniowym.

Mamy również tutaj instrukcję `output`, która mówi Nextflow'owi, jakie pliki mają zostać utworzone przez skrypt. Zauważ, że `output` tutaj ma słowo kluczowe `path`, które mówi Nextflow'owi, że to jest plik, a nie wartość ani ciąg znaków.

W bloku `script` jest to po prostu zwykła instrukcja bash — dokładnie taka sama jak to, co napisaliśmy w terminalu. Wykonujemy echo hello world do pliku o nazwie output.txt. Ten output.txt jest następnie przechwytywany przez definicję `output`. Definicja `output` tak naprawdę nic nie robi. Tylko mówi Nextflow'owi, czego się spodziewać — gdyby ten plik nie został utworzony, Nextflow zgłosiłby błąd.

Zauważ, że ten przykład nie jest świetny, ponieważ sztywno zakodowaliśmy nazwę pliku tutaj: output.txt i output.txt. Gdyby którykolwiek z nich został zmieniony, spowodowałoby to błąd w workflow.

Jest lepszy sposób, aby to zrobić za pomocą zmiennych — omówimy to za chwilę.

## 1.3 Definicja workflow

Dobra. Przechodząc do `workflow`, widzimy, że mamy komentarz, a następnie uruchamiamy proces o nazwie sayHello. To jest to samo słowo kluczowe, które jest tutaj na górze. To jest mniej więcej najprostszy workflow, jaki może być. Po prostu wywołujemy pojedynczy proces bez zmiennego wejścia, więc nie łączymy go z niczym innym. W późniejszej części tego kursu porozmawiamy o tym, jak uczynić to bardziej zaawansowanym, używając zmiennych wejść i łącząc rzeczy z kanałami.

## 2. Uruchom workflow

Dobra, to wszystko, czego potrzebujemy. Zobaczmy, czy możemy to uruchomić i zobaczyć, co się stanie. Wyczyśzczę terminal, a następnie zrobię "nextflow run" i wywołam nazwę pliku, która to hello-world.nf. To wszystko, czego potrzebujemy, aby uruchomić pipeline Nextflow'a. Ten pipeline nie przyjmuje żadnego wejścia, więc nie potrzebujemy żadnych innych argumentów.

Naciśnijmy enter i zobaczmy, co się stanie.

Dobra. Mam nadzieję, że masz jakieś wyjście, które wygląda tak. Mamy kilka informacji mówiących nam, że Nextflow został uruchomiony i jaka była używana wersja. Mówi nam, który skrypt został uruchomiony i podaje losowo wygenerowaną nazwę dla tego konkretnego wykonania workflow'u. W tym przypadku mój nazywał się "gloomy_crick".

Najważniejszą częścią jest jednak to, że mówi nam, które kroki zostały uruchomione w pipeline. Widzisz, że nasz proces o nazwie sayHello został uruchomiony. Został uruchomiony raz i został ukończony w stu procentach.

Ta część tutaj to hash dla tego konkretnego zadania workflow'u. Każdy proces działa jeden lub więcej razy, a każde z tych wykonań nazywa się zadaniem.

## 2.2. Znajdź wyjście i logi w katalogu work

Każde zadanie ma swój własny izolowany katalog, w którym działa — jest oddzielony od reszty wykonania workflow'u. Ten hash odpowiada strukturze plików w katalogu work. Jeśli zrobię "tree work", widzimy a0, a następnie dłuższą wersję krótkiego hasha, a następnie nasz plik output.txt. Możesz również zobaczyć to na pasku bocznym.

Widzisz na pasku bocznym, że jest tu kilka dodatkowych plików. Powodem, dla którego nie pojawiły się one w terminalu, jest to, że są to ukryte pliki — zaczynają się od kropki. I rzeczywiście, jeśli zrobię "tree -a" dla wszystkich i "work", widzimy je tutaj.

Te pliki z kropką są obecne w każdym pojedynczym katalogu work, który tworzy Nextflow. Każdy z nich ma nieco inne zadanie. Po pierwsze .command.begin zawiera po prostu kilka instrukcji dla Nextflow'a, które konfigurują zadanie przed jego uruchomieniem. .command.run to rzeczywiste instrukcje wykonywane przez sam Nextflow. Następnie .command.sh jest prawdopodobnie najbardziej interesujący. To jest skrypt, który został rozwiązany z naszego bloku `script` procesu.

Jeśli go otworzę, widzisz, że mamy nasze "echo Hello World" do pliku output.txt. To jest dokładnie to samo co nasz proces w tym przypadku, ale jeśli mamy jakiekolwiek zmienne w naszym kodzie Nextflow'a, każde zadanie będzie miało inny .command.sh i możesz zobaczyć, jak te zmienne zostały rozwiązane.

Pozostałe pliki dotyczą tego, jak zadanie zostało wykonane. Więc .command.err, .log i .out to standardowy błąd, standardowe wyjście i oba połączone. A .exitcode mówi Nextflow'owi, jak to zadanie zostało wykonane — z jakim kodem wyjścia, czy zakończyło się sukcesem, czy nie.

Na koniec mamy nasz plik output.txt i oczywiście "Hello World" — to jest to, czego oczekiwaliśmy i to zostało utworzone.

Dobra, świetnie. To było Twoje pierwsze uruchomienie Nextflow'a. Gratulacje. Naprawdę jest tak proste.

Następnie przejdziemy do tego, jak zrobić to trochę wygodniej, żebyśmy nie musieli edytować kodu za każdym razem, gdy chcemy dokonać zmiany w sposobie uruchamiania pipeline'u.

## 3. Zarządzaj wykonaniami workflow

Ta struktura katalogów jest świetna do utrzymywania wszystkich zadań oddzielonych i wszystkiego uporządkowanego, ale oczywiście nie jest zbyt wygodne znalezienie plików wyjściowych. Nie chcesz przeszukiwać mnóstwa zagnieżdżonych katalogów, próbując znaleźć wyniki Twojego pipeline'u.

## 3.1. Publikuj wyjścia

Dobra wiadomość jest taka, że nie musisz. Katalogi work są naprawdę tylko dla Nextflow'a do własnego użytku. Więc użyjemy funkcji Nextflow'a zwanej "publishDir".

Wracamy do workflow'u, idziemy do procesu. Możemy dodać tutaj nową instrukcję zwaną dyrektywą. To jest to, jak Nextflow nazywa te rzeczy na górze procesów, które rozszerzają sposób działania funkcjonalności — ta, której użyjemy, nazywa się publishDir.

Widzisz, że zacząłem tu pisać i rozszerzenie Nextflow dla VS Code zasugerowało mi dyrektywę, więc mogę po prostu nacisnąć enter.

Dobra. Podam katalog o nazwie "results" i powiemy mu, żeby skopiował tam pliki wyjściowe. Więc powiem mode copy. Świetnie. Zapiszę i uruchommy workflow ponownie.

nextflow run hello-world.nf

Działa dokładnie tak samo. Chociaż zauważ, że mamy teraz nieco inny hash. Nextflow użyje innego hasha za każdym razem, gdy uruchomisz workflow. W rezultacie mamy inny zestaw katalogów work. Jeden nazywa się EB, ale widzisz, że wszystkie pliki są takie same. Jednak tym razem nowością jest to, że mamy również katalog o nazwie "results".

W "results" tutaj mamy nasz plik wyjściowy. To jest to, co kazaliśmy zrobić Nextflow'owi. Powiedzieliśmy: zapisz pliki wynikowe w katalogu o nazwie "results" i skopiuj je tam. Teraz jest to znacznie łatwiejsze do znalezienia. Jest tam po prostu obok miejsca, w którym uruchomiliśmy workflow — wszystkie różne pliki mogą być tam zorganizowane, jak chcemy, niezależnie od tego, gdzie lub jak Nextflow uruchomił faktyczne wykonanie.

Zauważ, że publishDir może obsługiwać dowiązania symboliczne, co jest dobre, jeśli pracujesz na współdzielonym systemie plików i chcesz zaoszczędzić na miejscu. A także nie musisz definiować wszystkich plików, które są tworzone przez proces jako `output`.

Nextflow skopiuje tylko rzeczy, które są zdefiniowane w tym bloku `output`. Więc jeśli masz pliki pośrednie tworzone przez krok, które nie są potrzebne w dalszej części tego procesu, po prostu nie definiujesz ich w `output` i nie pojawią się w publishDir. To jest więc sposób na utrzymanie czystości plików wyjściowych z pipeline'u i łatwe usuwanie plików pośrednich po zakończeniu pracy.

Krótka uwaga tutaj. Jest nowa składnia Nextflow'a, która nadchodzi, nazywana definicjami wyjścia workflow, która ostatecznie zastąpi publishDir. Daje nam to sposób na zdefiniowanie wszystkich wyjść z workflow'u na poziomie pipeline'u w dół w bloku `workflow`. Jest to opisane w dokumentacji Nextflow'a, jeśli chcesz to wypróbować. Ale na razie publishDir będzie jeszcze przez jakiś czas, więc nadal mamy to w szkoleniu na rok 2025.

## 3.2. Uruchom ponownie workflow z -resume

Dobra. Wspomniałem, że katalog work ma teraz dwa zestawy wyników z innym hashem z każdego uruchomienia workflow'u. To dobrze. Jednak czasami nie chcemy ponownie obliczać kroków za każdym razem, jeśli nie musimy.

Może budujesz iteracyjnie Twój workflow i dodajesz kroki, a chcesz, aby pierwsze kroki po prostu użyły wersji z cache. Lub może coś poszło nie tak w Twoim systemie obliczeniowym w połowie workflow'u i chcesz, aby kontynuował od miejsca, w którym przerwał, ale pominął kroki, które już wykonał.

Nextflow ma wbudowaną funkcjonalność do tego, zwaną resume. Wypróbujmy to. Więc przede wszystkim, po prostu spojrzę na katalog work, żebyśmy mogli pamiętać, co tam było.

A następnie zrobię "nextflow run hello-world.nf" i dodam jedno polecenie tutaj: "-resume".

Zauważ, pojedynczy myślnik — to naprawdę ważne. Uruchomię to, a wyjście będzie wyglądać zasadniczo dokładnie tak samo, z kilkoma małymi różnicami.

Zauważ tutaj, że mówi "cached" na szaro. To oznacza, że Nextflow nie uruchomił zadania. Tym razem znalazł coś, co pasowało do wymagań, i użył tych wyjść bezpośrednio, zamiast ponownie uruchamiać krok.

I rzeczywiście, jeśli spojrzysz na hash tutaj, widzisz, że odpowiada istniejącemu hashowi, który mieliśmy z poprzedniego uruchomienia.

## 3.3. Usuń starsze katalogi work

Dobra. Ale jeśli rozwijasz iteracyjnie, zbudujesz wiele tych plików workflow'u. To może być problem, jeśli możesz mieć mało miejsca.

Nextflow może nam pomóc wyczyścić te katalogi work za pomocą kilku poleceń pomocniczych. Jeśli zrobię "nextflow log", to da mi listę wszystkich różnych uruchomień workflow'u, które wykonałem w tym katalogu. Mają tutaj nazwy uruchomień. Widzisz ten gloomy crick, który był pierwszym uruchomionym, a następnie te dwa nowe.

Możemy teraz wziąć tę nazwę i użyć ich z poleceniem "nextflow clean". Mogę określić pojedynczą nazwę uruchomienia. Lub jeszcze lepiej, mogę powiedzieć Nextflow'owi, żeby usunął wszystko sprzed pojedynczej nazwy workflow'u za pomocą "-before" i podam "stupefied_shaw". To było moje najnowsze uruchomienie, "-n".

Polecenie "-n" powiedziało Nextflow'owi, aby zrobił to jako próbne uruchomienie bez faktycznego usuwania czegokolwiek na prawdę — mówi nam, które z katalogów hash zostałyby usunięte. Rzeczywiście, to tylko ten z pierwszego wykonania. Oba drugie wykonania używają tego samego katalogu hash.

Uruchomię to ponownie, ale teraz zamiast "-n" dla próbnego uruchomienia, zrobię "-f" dla force — usunął ten katalog hash. Teraz jeśli zrobię "tree work", widzimy, że mamy tylko ten plik wyjściowy.

Świetnie. Więc udało nam się wyczyścić sporo miejsca na dysku.

Kilka rzeczy do zauważenia podczas usuwania katalogów work: jeśli stworzyłeś dowiązanie symboliczne do Twojego katalogu wyników, te źródła dowiązań symbolicznych zostaną teraz usunięte i Twoje wyniki znikną na zawsze. Dlatego użycie trybu copy jest bezpieczniejsze i ogólnie to zalecamy.

Po drugie, funkcjonalność resume Nextflow'a opiera się na tych katalogach work. Więc jeśli je usuniesz i uruchomisz Nextflow ponownie, funkcjonalność resume nie będzie już działać. Więc to do Ciebie należy śledzenie, których rzeczy możesz potrzebować lub nie — usuwaj rzeczy tylko wtedy, gdy masz pewność, że jest to bezpieczne.

Inna rzecz, którą możemy zrobić, to możemy po prostu usunąć cały katalog work, jeśli skończyliśmy uruchamianie workflow'u i mamy pewność, że już go nie potrzebujemy.

Więc mogę zrobić "rm -r work". Wiem, że nie było tam nic ważnego. Mam moje wyniki, na których mi zależy, w katalogu results, gdzie je skopiowaliśmy. I więc było bezpiecznie usunąć katalog work. To do Ciebie należy, które z tych podejść zastosujesz.

## 4. Użyj zmiennego wejścia przekazanego z linii poleceń

Dobra, co dalej? Wspomniałem, że sztywno zakodowaliśmy niektóre wartości w naszym skrypcie workflow'u tutaj — plik output.txt — i że może być lepszy sposób, aby to zrobić.

Zacznijmy od tego. Zrobimy trzy rzeczy. Dodamy nowe wejście do procesu. Powiemy skryptowi procesu, jak użyć tego wejścia, a następnie podłączymy to w workflow, abyśmy mogli używać tego dynamicznie z flagą linii poleceń podczas uruchamiania Nextflow'a.

Więc przede wszystkim, dodajmy blok `input` tutaj. Tak samo jak `output`. To jest nowa sekcja dla procesu i powiem "val greeting".

Zauważ tutaj, mówię "val", co oznacza, że to jest zmienna, a nie `path`.

Mogę wtedy przejść do skryptu i mogę wziąć ten sztywno zakodowany tekst tutaj i zrobić $greeting. To działa tak jak każdy inny język programowania. Definiujemy tutaj zmienną i odwołujemy się do niej w tym bloku `script`. Kiedy Nextflow uruchamia ten proces, zmienna zostanie zinterpretowana. A kiedy zajrzymy do tego pliku .command.sh, zobaczymy faktyczny sztywno zakodowany ciąg znaków tutaj zamiast tego.

## 4.1.3. Ustaw parametr CLI i podaj go jako wejście do wywołania procesu

Dobra, ale gdzie podajemy zmienną? Następnie przechodzimy do sekcji `workflow` i widzisz, że rozszerzenie tutaj mówi, że teraz oczekujemy wejścia i dało mi ostrzeżenie.

Teraz najprostszą rzeczą, którą moglibyśmy zrobić, jest po prostu sztywne zakodowanie tego. Mógłbym napisać "Hello World" i podać to wejście ciągu znaków do procesu. Ale znowu, to tak naprawdę nie rozwiązałoby żadnych problemów. Nadal musielibyśmy wracać i edytować kod pipeline'u za każdym razem, gdy chcielibyśmy coś zmienić — co nie jest dobre.

Dobra wiadomość jest taka, że Nextflow ma wbudowany system obsługi argumentów linii poleceń zwanych parametrami. Więc zamiast tego mogę użyć jednej z tych specjalnych zmiennych zwanych `params` i mogę nazwać to jak chcę, ale powiem greeting, żeby pasowało do logiki workflow'u.

Zapisz i zobaczmy, co możemy z tym zrobić.

Więc jeśli wrócę do terminala — robimy "nextflow run hello-world.nf". Tak jak wcześniej, ale kluczowa różnica polega na tym, że robimy --greeting

Zauważ, że są tutaj dwa myślniki, ponieważ to jest parametr. Kiedy wznawialiśmy workflow wcześniej, to był pojedynczy myślnik. To dlatego, że resume jest podstawową opcją Nextflow'a, a to jest parametr, który jest specyficzny dla naszego pipeline'u.

Nie myl tych dwóch. Łatwo to zrobić. Gdybyś zrobił --resume zamiast tylko jednego myślnika, to byłoby "params.resume", co by nic nie zrobiło. Podobnie, gdybyś zrobił pojedynczy myślnik tutaj, Nextflow nie rozpoznałby tego jako kluczowego argumentu.

Więc to --greeting, które odpowiada parameters greeting.

Mogę teraz śledzić to dowolnym tekstem, jaki chcę. Więc jestem obecnie w Szwecji, więc powiem "Hej världen".

Więc uruchommy to, zobaczmy, co się stanie — moment prawdy.

Dobra, więc widzisz, że proces został uruchomiony ponownie, tak jak wcześniej: sayHello z pojedynczym wykonaniem.

To nadpisze plik, który był w katalogu publishDir "results". Więc bądź ostrożny, gdy ponownie uruchamiasz pliki, ponieważ rzeczy w publishDir zostaną nadpisane.

Mogę teraz zrobić "code results/output.txt" i rzeczywiście, nasze wyjście zostało zaktualizowane i teraz mówi "Hej världen".

## 4.2. Użyj wartości domyślnych dla parametrów linii poleceń

Dobra, to świetnie. Ale problem polega teraz na tym, że nasz workflow polega na tym, że zawsze definiujemy ten parametr — a miło jest mieć rozsądne wartości domyślne, aby rzeczy działały w rozsądny sposób dla Twojego workflow'u, chyba że nadpiszesz domyślne.

Więc sposób, w jaki to robimy, to ustawienie wartości domyślnej dla parametru w naszym skrypcie workflow'u.

Więc jeśli wrócę do mojego pliku hello-world.nf, mogę przejść do skryptu tuż nad `workflow`, wpisać "params.greeting" i zdefiniować to jak każdą inną zmienną. Więc umieśćmy tutaj ciąg znaków i powiedzmy "Holà mundo!"

Teraz ten parametr ma zdefiniowaną wartość domyślną, która będzie używana tutaj, lub nadal możemy nadpisać to w linii poleceń za pomocą --greeting, tak jak robiliśmy wcześniej.

Więc sprawdźmy, czy działa. "nextflow run hello-world.nf"

Tym razem bez argumentów linii poleceń — sprawdźmy, czy zrobił właściwą rzecz.

"code results/output.txt". I proszę. Dostaliśmy naszą wartość domyślną.

Dobra, spróbujmy jeszcze raz, tylko sprawdźmy, że nie mówię nieprawdy. Uruchommy to ponownie, ale zróbmy --greeting i użyjmy przykładu z materiałów szkoleniowych, powiedzmy "Konnichiwa!"

Ponownie uruchamia workflow i rzeczywiście, nasz plik wyjściowy na górze został właśnie zaktualizowany nową wartością, którą podaliśmy w linii poleceń.

Świetnie. To jest naprawdę centralny aspekt pisania dowolnego workflow'u Nextflow'a. Definiowanie rozsądnych wartości domyślnych w kodzie pipeline'u, ale ułatwianie konfiguracji użytkownikowi końcowemu poprzez posiadanie argumentów linii poleceń w terminalu.

Zauważ, że użytkownik końcowy może nadpisać konfigurację w wielu różnych miejscach. Możesz mieć plik konfiguracyjny w Twoim katalogu domowym, który jest stosowany do każdego pojedynczego uruchomienia Nextflow'a, które robisz. Możesz mieć plik konfiguracyjny w katalogu startowym. Możesz mieć plik konfiguracyjny w katalogu pipeline'u. Wszystkie te różne lokalizacje konfiguracji są ładowane w określonej kolejności, która jest opisana w dokumentacji Nextflow'a.

Dobra, to koniec sekcji pierwszej. Mieliśmy nasz pierwszy skrypt workflow'u w Nextflow z `process` i `workflow`. Przyjrzeliśmy się wejściom, wyjściom, skryptom i publikowaniu oraz jak podłączyć parametry i kanał wejściowy do naszego procesu.

Gratulacje, Twój pierwszy krok w kierunku pisania kodu Nextflow'a jest ukończony.

Zrób sobie małą przerwę, a zobaczę Cię za kilka minut w rozdziale drugim.

[Następna transkrypcja wideo :octicons-arrow-right-24:](02_hello_channels.md)
