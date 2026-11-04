# Część 6: Hello Config - Transkrypcja wideo

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Tłumaczenie wspomagane sztuczną inteligencją - [dowiedz się więcej i prześlij sugestie](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

<div class="video-wrapper">
  <iframe width="560" height="315" src="https://www.youtube.com/embed/FcZTiE25TeA?si=tnXTi6mRkITY0zW_&amp;list=PLPZ8WHdZGxmWKozQuzr27jyMGqp9kElVK&amp;cc_load_policy=1" title="YouTube video player" frameborder="0" allow="accelerometer; autoplay; clipboard-write; encrypted-media; gyroscope; picture-in-picture; web-share" referrerpolicy="strict-origin-when-cross-origin" allowfullscreen></iframe>
</div>

!!!note "Ważna informacja"

    Ta strona zawiera wyłącznie transkrypcję. Aby uzyskać pełne instrukcje krok po kroku, wróć do [materiałów kursu](../06_hello_config.md).

    Numery sekcji widoczne w transkrypcji mają charakter orientacyjny i mogą nie obejmować wszystkich numerów sekcji obecnych w materiałach.

## Powitanie

Witaj ponownie w części szóstej szkolenia Hello Nextflow. Ta sekcja dotyczy plików konfiguracyjnych i jest ostatnią częścią tego kursu.

Nextflow szczególnie dobrze radzi sobie z dwoma rzeczami: odtwarzalnością i przenośnością. Pliki konfiguracyjne to miejsce, gdzie widać prawdziwą siłę tej drugiej cechy. Możliwość skonfigurowania pipeline'u Nextflow'a tak, aby działał na różne sposoby i funkcjonował w różnych systemach, bez konieczności edycji podstawowego kodu pipeline'u.

Ta supermoc pozwala na ponowne wykorzystanie pipeline'ów Nextflow'a przez inne osoby w innych miejscach lub w różnych infrastrukturach, do których Ty sam możesz mieć dostęp.

Oznacza to, że możesz rozwijać kod pipeline'u na swoim laptopie, przenieść go do chmury, uruchomić na swoim HPC - to ten sam kod pipeline'u, który działa wszędzie.

W tej sekcji przejdziemy przez kilka tematów. Zaczniemy od tego, jak Nextflow obsługuje pliki konfiguracyjne, skąd je ładuje, jak je pisać i jak je strukturyzować, oraz od rozdzielenia między samym pipeline'em a tym, co powinno trafić do pliku konfiguracyjnego.

Następnie przejdziemy do kilku typowych przypadków użycia, takich jak zmiana miejsca przechowywania plików wyjściowych, a także do sposobu na to, by pipeline działał w różnych infrastrukturach - zarówno przy użyciu różnych rodzajów pakietów oprogramowania, jak i przy przesyłaniu zadań do różnych infrastruktur.

## Hierarchia plików konfiguracyjnych

Dobra, zaczynajmy. Jeśli chodzi o ładowanie plików konfiguracyjnych, Nextflow może pobierać je z wielu różnych miejsc, co jest dobrą rzeczą, a jednocześnie może być nieco ryzykowne, ponieważ czasami może być trudno wiedzieć, skąd pobiera plik konfiguracyjny i w jakiej kolejności ładuje różne elementy.

Dlatego naprawdę polecam kliknięcie tego linku, który prowadzi do dokumentacji Nextflow'a. Na tej stronie konfiguracji wymienione są kluczowe miejsca, z których ładowana jest konfiguracja, i, co ważne, kolejność, w jakiej te elementy są ładowane.

Widzisz, że możesz umieścić plik konfiguracyjny w katalogu domowym Nextflow'a, który zwykle to ".nextflow" w Twoim katalogu domowym. Ten plik będzie zawsze ładowany przy każdym uruchomieniu Nextflow'a w Twoim systemie.

Kolejne miejsce to plik w głównym katalogu Twojego repozytorium lub katalogu pipeline'u o nazwie "nextflow.config".

Następnie kolejny plik o nazwie "nextflow.config", ale tym razem w katalogu, z którego uruchamiasz Nextflow'a: w katalogu uruchomieniowym.

Na koniec możesz podać ścieżki do plików konfiguracyjnych w linii poleceń przy użyciu argumentu "-c", i możesz to zrobić wielokrotnie. Są one stosowane w kolejności, w jakiej je określisz.

Możesz podać pliki konfiguracyjne we wszystkich tych lokalizacjach, jeśli chcesz, i będą one ładowane iteracyjnie, każdy nadpisując poprzedni tylko w tych zakresach konfiguracji, w których kolidują.

To naprawdę potężny system, ponieważ oznacza, że możesz ustawić rozsądne wartości domyślne, a następnie stawać się coraz bardziej szczegółowym, w miarę jak precyzujesz konfigurację.

## 0. Rozgrzewka: Uruchom hello-config.nf

Dobra, zamknijmy to i przejdźmy do naszego Codespaces, żeby zacząć. Jak zwykle posprzątałem tutaj, usunąłem poprzednie katalogi z wynikami, moje katalogi Nextflow i work i tak dalej. Nie martw się, jeśli nadal masz te pliki. To dlatego, że jestem bardzo przybliżony, więc inaczej szybko robi się tu bałagan.

Będziemy pracować z hello-config.nf, ostatnim plikiem w naszym katalogu, i powinien on kontynuować to, na czym skończyliśmy w poprzedniej sekcji.

Mamy więc nasze cztery różne procesy, które są dołączane z plików modułowych. Mamy parametry naszego pipeline'u, blok workflow, w którym wywołujemy różne procesy i łączymy kanały, publikujemy kanały wyjściowe, a następnie blok output na dole, gdzie definiujemy, gdzie te pliki powinny być przechowywane i jak powinny być kopiowane.

Mamy też już plik "nextflow.config" z ostatniego rozdziału, gdzie włączyliśmy Docker'a, i będziemy dzisiaj rozbudowywać ten plik.

Jak poprzednio, zmieniliśmy ścieżkę wyjściową w tym głównym skrypcie na hello config, żeby nie kolidowała z poprzednimi wynikami, które wygenerowałeś.

Dobra, sprawdźmy szybko, czy wszystko nadal działa zgodnie z oczekiwaniami. Otwórzmy terminal i wykonajmy nextflow run hello-config.nf. Nextflow się ładuje. Powinien uruchomić nasze cztery różne procesy. Wygenerować trochę ładnej grafiki ASCII przy użyciu cowpy, a następnie zapisać nasze wyniki do plików wynikowych w tym katalogu.

Mogę szybko zajrzeć tutaj, żeby upewnić się, że te pliki wyglądają zgodnie z oczekiwaniami, i rzeczywiście, tam jest nasz gigantyczny indyk. Świetnie.

## 1.1. Przenieś wartości domyślne do nextflow.config

Teraz pierwszą rzeczą, którą zrobimy, jest przeniesienie niektórych elementów z naszego skryptu do pliku konfiguracyjnego.

A najbardziej zależy nam na parametrach na tym etapie. Chcemy przenieść wartości domyślne do pliku konfiguracyjnego, żeby było jaśniejsze, jakie są wartości domyślne i łatwiej było je nadpisać.

Wezmę ten blok params tutaj ze skryptu i umieszczę go w pliku konfiguracyjnym. Musimy być tutaj trochę ostrożni, bo składnia różni się nieco między plikami konfiguracyjnymi a skryptami. Plik konfiguracyjny nie może przyjmować deklaracji typów, bo tak naprawdę nie definiujemy tych parametrów, tylko się do nich odwołujemy. Więc pozbędę się tych deklaracji.

Ale poza tym jest bardzo podobnie. Mamy blok params, a następnie mamy nasze różne parametry wejściowe, parametr batch, parametr character.

Mogę teraz wrócić do mojego skryptu i nie muszę już definiować tych wartości domyślnych, ponieważ te wartości są teraz w moim pliku Nextflow config.

Jednak zostawiam nazwy parametrów i ich typy, żeby Nextflow znał te informacje i nadal mógł wykonywać kontrolę typów i wszystko inne.

Dobra. Zapiszmy te pliki i szybko sprawdźmy, czy wszystko nadal działa tak samo jak wcześniej. Nie powinno być żadnych zmian. Zachowaliśmy te same wartości. Po prostu przenieśliśmy miejsce, w którym zostały zdefiniowane.

Świetnie.

## 1.2. Użyj pliku konfiguracyjnego specyficznego dla uruchomienia

Teraz, do tej pory uruchamialiśmy Nextflow'a z tego samego katalogu, w którym mamy nasz skrypt pipeline'u. Więc nasz katalog uruchomieniowy i katalog pipeline'u to w zasadzie to samo.

Aby pokazać, jak możemy mieć różne pliki konfiguracyjne z różnymi katalogami uruchomieniowymi, utworzymy teraz nowy podkatalog.

Więc powiem mkdir, i nazwiemy go tux-run.

A potem wykonam cd, zmienię katalog na tux-run. Zauważ, że teraz nasz katalog roboczy nie jest już w tym samym katalogu co skrypty pipeline'u.

Dobra, stwórzmy nowy plik "nextflow.config". Więc touch Nextflow config, i po prostu otwórzmy go w VS Code. Widzisz też na pasku bocznym tutaj, że jesteśmy teraz w tym podkatalogu.

Teraz możemy wziąć ten sam blok params, który mieliśmy w głównym nextflow.config, skopiować go tutaj i teraz możemy zmienić te wartości.

Po pierwsze, data to teraz inna ścieżka względna, ponieważ jesteśmy w podkatalogu, więc musimy to zaktualizować. A następnie zmienimy batch na experiment i zmienimy character z Turkey na tux.

Teraz klikam zapisz, i wypróbujmy to. Tak jak w przypadku data, muszę teraz powiedzieć ../ żeby dostać się do skryptu. Więc to Hello config. I naciskam enter.

Kod pipeline'u w ogóle się nie zmienił, ale teraz będziemy mieć dwa zestawy ładowanej konfiguracji, a plik konfiguracyjny z katalogu uruchomieniowego powinien nadpisać wartości domyślne, które zostały ustawione w pipeline'owym nextflow.config, i powinniśmy otrzymać różne zestawy wyników.

Rzeczywiście, w naszym katalogu tutaj, w tux-run, widzisz, że mamy katalog dot Nextflow i katalog work i to dlatego, że są one zawsze tworzone w Twoim katalogu uruchomieniowym. Więc są to inne katalogi work i results niż te z wcześniejszych uruchomień.

Teraz, jeśli spojrzę w results, możemy zobaczyć nasz plik collected i jest tam nasz mały znak tux. Więc widzisz, że te parametry zostały prawidłowo zinterpretowane.

## 1.3. Użyj pliku parametrów

Dobra. Wcześniej, gdy mówiłem o różnych plikach konfiguracyjnych, które można załadować, pominąłem jedno inne miejsce, z którego możemy pobrać konfigurację.

Możemy ją uzyskać z linii poleceń, jak widzieliśmy z myślnikami i nazwami parametrów, ale możemy też dostarczyć plik YAML lub JSON, składający się tylko z parametrów.

Plik konfiguracyjny może zawierać różne typy zakresów, ale te pliki to tylko parametry, i jest to przyjazny dla użytkownika sposób na dostarczenie wielu parametrów na raz, i być może nieco bardziej odtwarzalny sposób, ponieważ zapisujesz je do pliku, więc łatwo je odzyskać później.

Więc wróćmy do naszego terminala i zanim zapomnimy, upewnijmy się, że wracamy o katalog wyżej, więc nie jestem już w podkatalogu, i zobaczę plik YAML, który tu mamy o nazwie test-params.yaml.

Więc jeśli po prostu zrobię code test-params.yaml, widzisz, że to tylko zwykły plik YAML. Nic specjalnego. Z kluczami będącymi naszymi nazwami parametrów, z formatowaniem YAML, więc dwukropek tutaj, a następnie wartość.

Zauważ, że to nie jest kod Nextflow'a, więc nie możemy tu umieszczać rzeczy takich jak zmienne. To tylko wartości statyczne.

Również dlatego, że JSON w rzeczywistości jest parsowany jako YAML, możemy również mieć plik test-params.json, który wygląda bardzo podobnie. To po prostu inny format danych.

Więc mamy tutaj dwa różne pliki testowe i mamy nieco inne zmienne.

Dobra, więc jak przekazujemy je do Nextflow'a? To bardzo proste. Robimy Nextflow run hello config, jak wcześniej. I zamiast "-c" dla pliku konfiguracyjnego, lub ładowania tych domyślnych nazw plików, robimy -params-file. Pojedynczy myślnik, ponieważ to podstawowa opcja Nextflow'a.

A następnie podajemy ścieżkę do tego pliku. Więc zrobię "-params-file test-params.yaml", i zobaczymy, czy są prawidłowo załadowane.

Dobra. Uruchomił się. Przypomnijmy sobie, co było w tym pliku YAML. Więc batch został ustawiony na YAML, więc tak powinien się nazywać, i powinien mieć stegozaura. Więc przejdźmy i zobaczmy w results. I mamy COLLECTED-yaml. Więc zobaczmy, czy mamy Stegozaura. Fantastycznie, Stegozaur w kapeluszu. Tego nam było trzeba.

Więc to zadziałało naprawdę dobrze, i dokładnie tak samo jest z plikiem JSON. Po prostu zamieniamy rozszerzenie pliku tutaj i Nextflow wie, jak to odczytać.

I w tym przypadku powinniśmy mieć batch o nazwie JSON i powinniśmy mieć żółwia. Zobaczmy. Wspaniale. Jedno z moich ulubionych narzędzi CLI.

## 2.1. Dostosuj katalog wyjściowy przy użyciu -output-dir

Dobra, więc to było głównie myślenie o wejściach do pipeline'u i zmianie parametrów. A co z wyjściami?

Teraz, chociaż zmienialiśmy podkatalogi przy użyciu params, mogłeś zauważyć, że wszystkie nasze pliki nadal trafiają do results.

Możemy zmienić ten katalog bazowy, do którego publikowane są wszystkie pliki, za pomocą flagi linii poleceń o nazwie -output-dir. Więc jeśli zrobię Nextflow run hello config, a następnie zrobię -output-dir, i nazwiemy to "custom-outdir-cli". Nie mogę pisać. Żebyśmy pamiętali, skąd pochodzą te pliki.

To podstawowa opcja Nextflow'a i jest to bardzo nowa opcja. Została dodana dopiero niedawno, i to jedna z rzeczy, które możemy zrobić z nowym parserem języka i wszystkim.

To trochę długo się pisze. Możesz też po prostu nazwać to "-o", jeśli chcesz. Więc jeśli wrócę. Mogę po prostu skrócić to do "-o", co jest nieco prostsze.

Dobra. Uruchamiamy to. Nie zmieniliśmy niczego w naszym pipeline'ie ani nawet w naszej konfiguracji na tym etapie, i powinno, mam nadzieję, zapisać wszystkie nasze wyniki do innego katalogu najwyższego poziomu. I możesz sobie wyobrazić, że możesz ustawić to na zasadniczo dowolną ścieżkę, którą chcesz.

Właśnie się pojawił u góry. Mamy custom-outdir-cli, i wszystkie pliki są tam zorganizowane w dokładnie taki sam sposób, z tymi samymi podkatalogami i nazwami plików. Więc to naprawdę łatwy sposób na zmianę miejsca, gdzie pipeline publikuje swoje wyniki, bez zastanawiania się zbytnio nad tym, jak te wyniki są zorganizowane.

## 2.1.2. Usuń zakodowane ścieżki z bloku output

Jeśli zajrzę do tego katalogu, możemy zobaczyć, że nadal mamy podkatalog o nazwie Hello Config, co wydaje się teraz trochę zbędne.

Więc załadujmy nasz skrypt ponownie i możemy teraz usunąć ten podkatalog z bloku output na dole. Bo tak naprawdę nie potrzebujemy go już. Więc możemy po prostu to teraz zrobić, usunąć to stąd. A jeśli to tylko to, możesz albo całkowicie to usunąć, albo zostawić jako pusty ciąg znaków. Zostawię to jako pusty ciąg znaków na razie, ponieważ wrócimy i umieścimy tam kilka różnych rzeczy w przyszłości. Ale jeśli nie dbasz o podkatalogi, najczystszym rozwiązaniem jest całkowite usunięcie deklaracji path.

Dobra, zapiszmy. Szybko wypróbujmy to ponownie. Właściwie usunę mój katalog "custom-outdir-cli", żebyśmy nie byli zdezorientowani przez istniejące pliki tam. Pamiętaj, że kiedy publikujesz rzeczy, nie usuwa to plików, które były tam wcześniej. Po prostu dodaje nowe. Uruchommy to polecenie ponownie, custom-outdir-cli.

A teraz jeśli zrobisz "ls custom-outdir-cli", nie ma już katalogu o nazwie Hello Config.

## 2.2.1. Ustaw outputDir w pliku konfiguracyjnym

Dobra, flaga linii poleceń tutaj, "-o" lub "-output-dir" jest dobra. Ale co z ustawianiem wartości domyślnych dla tego w konfiguracji? Jak to zrobić?

Otwieram plik "nextflow.config", zamykam wszystko inne i pozbywam się tego. Możemy dodać nową opcję konfiguracyjną tutaj, którą właśnie skopiowałem z materiałów szkoleniowych, i nazywa się outputDir.

Nie znajduje się pod żadnym zakresem. Nie pod params ani niczym. Jest na najwyższym poziomie, i możemy ustawić to na ciąg znaków. Teraz prostą rzeczą do zrobienia jest po prostu zmiana tego na cokolwiek innego niż results jako zakodowany ciąg. Ale ponieważ to jest w pliku konfiguracyjnym Nextflow'a, możemy być tutaj trochę sprytni i również włączyć zmienne.

I widzisz tutaj, że włączyliśmy zmienną params, params.batch, która jest częścią tego ciągu znaków. To oznacza, że możemy ponownie wykorzystywać zmienne, które pochodzą z innych miejsc. I w tym przypadku, jeśli zrobimy --batch, gdy uruchamiamy Nextflow Pipeline, otrzymamy podkatalog w naszej niestandardowej ścieżce w oparciu o to, jaka była nazwa batch.

Dobra, więc wypróbujmy to i po prostu szybko zobaczmy, jak, jak wyglądają wyniki. Więc jeśli zrobię Nextflow run hello config i --batch my_run. Przypomnijmy sobie, jak wyglądała konfiguracja. Więc to custom-outdir-config.

Tree custom-outdir-config. I widzisz, batch nazywał się my_run. A potem mamy ten podkatalog o nazwie my_run. Więc ta dynamiczna ścieżka pliku zadziałała.

I nie tylko to, nie trafił już do domyślnego katalogu results, i nie musiałem określać niczego w linii poleceń, żeby zmienić katalog bazowy. Więc pomyślnie zresetowaliśmy wartość domyślną dla domyślnego outputDir.

## 2.2.2. Podkatalogi z nazwami batch i procesów

Dobra, pójdźmy o krok dalej. To dynamiczna zmienna w pliku konfiguracyjnym. A co ze skryptem? Teraz, do tej pory mieliśmy te ścieżki tutaj i one również mogą być dynamiczne. Więc zamiast po prostu kodować coś na sztywno, możemy umieścić kilka nawiasów klamrowych i umieścić coś dynamicznego.

Na przykład mamy nasze procesy o nazwie sayHello. Moglibyśmy zrobić sayHello.name, który jest atrybutem procesu, co jest trochę nudne, bo to po prostu "sayHello" w tym przypadku. Ale jest zmienne.

Więc to daje Ci pomysł. Więc możemy umieścić to tutaj i powiedzieć convertToUpper.name, collectGreetings.name, collectGreetings.name znowu, i cowpy.

Teraz kiedy uruchomimy, katalog bazowy nadal będzie custom-outdir-config. I będzie w podkatalogu o nazwie params.batch, ale podkatalogi pod tym powinny być zorganizowane według nazwy procesu.

Wypróbujmy to i zobaczmy, czy to działa. Więc usunę poprzedni katalog, żebyśmy się nie pomylili, i użyję dokładnie tego samego polecenia Nextflow Run.

Powinno uruchomić się w ten sam sposób. Mogłem używać dash resume na wszystkich tych, żeby było trochę szybciej i użyć wcześniej obliczonych wyników. Teraz, jeśli zrobię tree custom-outdir-config, widzisz, że to nie jest w results, to w naszym katalogu bazowym z nazwą batch. I widzisz, że wszystkie wyniki są teraz zorganizowane w podkatalogach nazwanych według procesu. Więc mamy dwa różne miejsca, gdzie definiujemy dynamiczne ścieżki wyjściowe tutaj.

Dobra. Ostatnia rzecz, dodajmy z powrotem te foldery pośrednie, które mieliśmy wcześniej, bo były dość przyjemne. Intermediates.

I możemy również pomyśleć trochę o tym params.batch, może jako programista pipeline'u naprawdę lubiłem mieć to w podkatalogu, ale jeśli użytkownicy końcowi pipeline'u ustawiają "-o" lub -output-dir w CLI, całkowicie nadpisuje to całą instrukcję, i tracimy ten podkatalog.

Więc możemy wziąć tę dynamiczną ścieżkę z konfiguracji outputDir, która zostałaby nadpisana, i umieścić ją w output path, która nie jest nadpisywana.

Więc możemy zrobić params.batch slash intermediates slash sayHello.name, i zrobić to wszystko w ciągu znaków w podwójnych cudzysłowach, żeby było interpolowane przez Nextflow'a.

Mogę teraz skopiować, ups. Skopiować to do innych procesów. Pamiętaj, żeby umieścić je wszystkie w cudzysłowach. I usunąć intermediates z tych konkretnych wyjść.

Dobra? Wygląda to teraz nieco bardziej złożenie, ale widzisz, że naprawdę zaczynamy budować ładnie zorganizowaną strukturę katalogów wyjściowych w naszym kodzie.

I co naprawdę fajne, to że ta dodatkowa złożoność w kodzie nie przechodzi do CLI. Więc możemy uruchomić nasze polecenie z -output-dir i jakimikolwiek zmiennymi batch, po prostu myśląc o tym, jak uruchomić pipeline i nie myśląc zbytnio o tym, co jest w kodzie. A nasze pliki wyjściowe zostaną skonstruowane naprawdę ładnie w bardzo dobrze zorganizowany sposób, co jest miłe dla osób używających pipeline'u.

Świetnie. Gdy to piszę, zdaję sobie sprawę, że popełniłem błąd. Zobaczcie, czy ktoś mnie przyłapał. Mamy collectGreetings.name, więc coś poszło nie tak. I tak, rzeczywiście, przypadkowo zapomniałem umieścić to w nawiasach klamrowych.

Więc pamiętaj, bądź ostrożny, gdy piszesz swój kod i upewnij się, że mówisz Nextflow'owi, co jest zmienną, a co tylko ciągiem znaków. Bo zrobi dokładnie to, co mu powiesz. I nic więcej. Jak wszystkie dobre komputery. Dobra, to powinno to naprawić.

## 2.3. Ustaw tryb publikowania na poziomie workflow'u

Jest jedna część tego skryptu, której nadal nie lubię, a jest to fakt, że piszemy mode copy wielokrotnie, a jeśli jest jedna rzecz, której nie lubimy, to powtarzanie się.

Więc możemy to trochę posprzątać, biorąc to i przenosząc do konfiguracji. I w rzeczywistości możemy ustawić to dla całego pipeline'u za jednym razem. Więc nie musimy mówić tego wielokrotnie.

Przechodzimy do naszego pliku konfiguracyjnego i mamy nowy zakres tutaj o nazwie workflow. I możemy użyć nawiasów klamrowych lub możemy użyć notacji kropkowej. Nie ma to znaczenia. Materiały szkoleniowe używają notacji kropkowej. Mogę powiedzieć output i możemy mieszać i dopasowywać, więc mode equals copy. Świetnie.

I teraz możemy wrócić tutaj i je usunąć. Moglibyśmy je zostawić na miejscu. Konfiguracja zasadniczo nadpisuje to, co jest tutaj napisane, ale ponieważ mamy to w konfiguracji na poziomie pipeline'u, i te dwa pliki są dostarczane razem, nie ma powodu, żeby to robić dwukrotnie.

Dobra. Tylko sprawdźmy się, bo najwyraźniej popełniamy błędy. Uruchommy to ponownie i sprawdźmy, czy prawidłowo używamy trybu copy do publikowania plików. Więc uruchomimy skrypt ponownie i tym razem umieściliśmy wyniki w katalogu o nazwie config-output-mode, zobaczmy, jak wyglądają tam pliki.

A potem jeśli zrobię "ls -l", żeby spojrzeć na batch, i możemy spojrzeć na cowpy, na przykład. I powinniśmy zobaczyć, tak, że to jest właściwy plik tutaj, który nie jest dowiązaniem symbolicznym, więc ten atrybut konfiguracji został zastosowany prawidłowo.

## 3. Wybierz technologię pakowania oprogramowania

Dobra. Do tej pory koncentrowaliśmy się na wejściach i wyjściach, plikach, z którymi workflow pracuje. Ale co z infrastrukturą? Powiedziałem na początku, że Nextflow pozwala uruchamiać ten sam pipeline w różnych konfiguracjach obliczeniowych. Więc jak to wygląda?

Żeby to pokazać, przełączymy się z używania Docker'a do uruchamiania cowpy, i zamiast tego użyjemy Condy do tego samego.

Mogę to zrobić bardzo prosto. Jeśli przejdę do code, "nextflow.config". Jeśli pamiętasz na górze, zdefiniowaliśmy docker.enabled wcześniej, w ostatnim rozdziale, żebyśmy mogli użyć kontenera z cowpy.

Powiem Nextflow'owi, żeby nie używał Docker'a. Ustawię to na false. I powiem Conda enabled equals true. Więc powiem Nextflow'owi, proszę użyj Condy.

Teraz samo włączenie Condy nie wystarczy. Tak jak zrobiliśmy z Docker'em, musimy powiedzieć Nextflow'owi, skąd może pobrać potrzebne oprogramowanie.

Więc jeśli przejdziemy do modules tutaj. I otworzymy skrypt cowpy. Możemy zobaczyć, że mamy deklarację container na górze. I container jest używany przez Docker'a, ale także Singularity, Apptainer i wiele innych narzędzi programowych.

Ale nie może być używany dla Condy, więc mamy oddzielną deklarację o nazwie "conda", i moglibyśmy po prostu napisać "cowpy". I to pozostawi rozwiązanie pakietu conda, aby wymyślić najlepszy sposób na rozwiązanie tego, zgodnie z Twoim lokalnym środowiskiem conda.

Lub dobrą praktyką jest zrobić to, co mówią materiały szkoleniowe, czyli zdefiniować konkretny kanał conda z jego notacją podwójnego dwukropka i zdecydowanie zdefiniować konkretną wersję oprogramowania, tak aby każda osoba, która uruchamia pipeline, otrzymała tę samą wersję.

Zauważ, że kontenery są nieco lepsze pod tym względem, ponieważ kiedy instalujesz coś z Condą, nadal będzie rozwiązywać wszystkie zależności dla tego pakietu, i mogą się one zmieniać z czasem. Nazywa się to dryfem zależności.

Więc kontenery jednak blokują cały stos zależności oprogramowania aż do samego dołu, więc możesz być nieco bardziej pewny, że A, to zadziała, i B, będzie odtwarzalne.

Więc jeśli jesteś w stanie używać Docker'a lub Singularity lub Apptainer'a, zdecydowanie polecam.

Teraz to, co jest fajne w tym, to to, że plik modułu, który jest napisany przez programistę pipeline'u, ma teraz zarówno Container jak i Conda, i więc mówimy osobie, która uruchamia ten pipeline, nie obchodzi nas, jakiego rozwiązania pakowania oprogramowania używasz. Będzie działać zarówno z Docker'em jak i z Condą, i to jest miejsce, gdzie pobrać oprogramowanie w obu przypadkach.

Możemy otworzyć terminal i wypróbujmy to. Więc Nextflow run hello config --batch conda. I po raz pierwszy, gdy to uruchomi się z condą, będzie trochę wolno, gdy dojdzie do tego konkretnego procesu, ponieważ musi uruchomić "conda install".

I tworzy specjalne środowisko conda tylko dla tego jednego procesu. Więc nie używa mojego globalnego środowiska conda, które mam w moim terminalu. Tworzy jedno tylko dla tego jednego procesu. To jest dobre, ponieważ unika rzeczy takich jak konflikty zależności między różnymi procesami w Twoim workflow'ie. Jeśli Twoje procesy mają narzędzia, które potrzebują różnych wersji Pythona lub podobnych rzeczy, to jest okej, ponieważ używają różnych środowisk conda.

Nextflow cache'uje te środowiska conda lokalnie, widzisz, że mówi Ci, gdzie jest ta ścieżka, jest w katalogu work tutaj. I więc następnym razem, gdy uruchomię ten skrypt z Condą, będzie znacznie szybciej, bo znajdzie to istniejące środowisko conda i po prostu je ponownie użyje. Ale za pierwszym razem musi pójść i je pobrać, rozwiązać, pobrać wszystkie zależności i wszystko skonfigurować.

Dobra, świetnie, uruchomiło się. Możemy sobie po prostu przypomnieć, do czego pipeline jest obecnie skonfigurowany. Jeśli spojrzymy w plik konfiguracyjny, to było "custom-outdir-config" teraz dla mnie. Zobacz, czy przejdę do tego katalogu bazowego. I zrobiłem --batch conda. Jest nasz podkatalog conda. Więc to zadziałało i jest nasze wyjście cowpy.

Więc pobrało cowpy, zainstalowało go w moim lokalnym systemie używając condy i uruchomiło proces. I co świetne, jako użytkownik końcowy nie musiałem w ogóle myśleć o żadnym zarządzaniu oprogramowaniem. Nextflow po prostu to dla mnie załatwił. Powiedziałem, muszę użyć condy w tym systemie. Programista pipeline'u powiedział, których pakietów potrzebowałem. A Nextflow zrobił resztę. Bardzo potężne.

Zauważ, że możesz faktycznie używać mieszanki różnych technologii. Więc mogę włączyć Docker'a dla konkretnych procesów, i condę dla innych procesów, lub powiedzieć, że niektóre procesy powinny po prostu użyć jakiegokolwiek lokalnego oprogramowania, które zainstalowałem. To dość niezwykłe, ale jest możliwe, i w niektórych przypadkach, na przykład, jeśli używasz określonego oprogramowania, które może być trudne do zapakowania w Docker'ze, masz wyjście awaryjne.

## 4. Wybierz platformę wykonawczą

Więc to jest pakowanie oprogramowania. Druga część przenośności do innych systemów to miejsce, gdzie faktycznie wykonują się zadania. W tej chwili uruchamiam się zasadniczo na moim laptopie lub w tym Codespaces, co jest pojedynczym komputerem. Nic specjalnego. Nextflow jest trochę sprytny w paralelizacji zadań najlepiej jak potrafi, ale to wszystko na jednym systemie.

Teraz, jeśli uruchamiasz się na HPC, prawdopodobnie masz jakiś harmonogram zadań, taki jak SLURM lub PBS lub coś, i wysyłasz zadania do tego harmonogramu, a on rozdziela wszystkie zadania do różnych węzłów obliczeniowych.

Innym sposobem uruchamiania jest chmura. Więc może używasz AWS Batch, lub Azure Cloud, lub Google. I wszystkie te działają w podobnym systemie, gdzie masz harmonogram i wysyłasz zadania i są one przesyłane do różnych miejsc do obliczenia.

Teraz w odległej przeszłości, gdy zaczynałem zajmować się bioinformatyką, oprogramowanie wszystkich do uruchamiania analiz było bardzo związane z ich infrastrukturą obliczeniową, co czyniło replikację prawie niemożliwą.

Ale z tym rozdzieleniem konfiguracji w Nextflow'ie, i ze zdolnością Nextflow'a do interakcji z bardzo wieloma różnymi backendami infrastruktury obliczeniowej, bardzo proste jest wzięcie naszego pipeline'u bez modyfikacji kodu pipeline'u w ogóle i po prostu to przełączenie.

## 4.1. Kierowanie na inny backend

Więc jeśli przejdziemy do naszego pliku "nextflow.config", i możemy teraz umieścić trochę konfiguracji na poziomie procesu. Więc jeśli umieszczę na górze zakres process i mogę ustawić executor, i tutaj jest ustawiony na local, co jest wartością domyślną.

Zauważ, że ponieważ to jest poziom procesu, możemy kierować rzeczy do różnych procesów. I więc możesz faktycznie ustawić executory, aby były specyficzne dla procesu i mieć hybrydowe wykonanie, gdzie niektóre zadania mogą być uruchamiane lokalnie, gdziekolwiek zadanie Nextflow'a jest wykonywane. Niektóre są przesyłane do różnych HPC, a niektóre mogą być przesyłane do chmury. Możesz być tak sprytny, jak chcesz.

Teraz bardzo trudno jest to zademonstrować w środowisku szkoleniowym takim jak to, ponieważ nie mam HPC, do którego mogę przesyłać. Ale mogę zrobić, jeśli wpiszę slurm, możemy trochę oszukać i możesz poczuć to.

I to jest najbardziej interesujące dla osób, które są przyzwyczajone do uruchamiania na SLURM i wiedzą, jak wyglądają nagłówki SLURM. Ale jeśli zrobię Nextflow run, hello config. To się nie powiedzie, ponieważ będzie próbować przesłać zadania do klastra, który nie istnieje. Więc dostaniemy jakiś błąd o tym, że sbatch nie jest dostępny.

Tak, napisane. To jest narzędzie. To narzędzie CLI, którego używasz do przesyłania zadań do klastra slurm. Ale co możemy zrobić, to możemy przejść i zajrzeć do naszego katalogu work tutaj przez command, kliknąć, otworzyć ten katalog i spojrzeć na .command.run. I widzisz na górze pliku .command.run, mamy nasze nagłówki sbatch, mówiące teoretycznemu klastrowi SLURM, jak obsłużyć to przesłanie zadania.

I więc widzisz, że Nextflow jest sprytny, robi wszystkie właściwe rzeczy. Po prostu nie mieliśmy klastra, do którego moglibyśmy przesłać.

## 5. Kontroluj alokację zasobów obliczeniowych

Co jeszcze jest inne między różnymi infrastrukturami obliczeniowymi? Inną rzeczą jest to, ile dostępnych zasobów masz, i w rzeczywistości w wielu środowiskach obliczeniowych jest to wymóg, że musisz określić, ile CPU i ile pamięci zadanie potrzebuje.

Znowu, Nextflow abstrahuje to dla nas, tak że nie jest już specyficzne dla pojedynczego typu środowiska obliczeniowego, i możemy wpisać w zakresie poziomu procesu tutaj. CPUs equals one, memory equals two gigabytes. Nasz pipeline nie jest zbyt wymagający, więc to powinno być w porządku.

Teraz po prostu zgadłem te liczby tutaj, ale skąd wiesz, jaka jest rozsądna ilość zasobów do użycia? To dość trudne zadanie, żeby przejść i przeszukać wszystkie te różne procesy dużego pipeline'u wielu próbek i zrozumieć, jakie było wykorzystanie zasobów.

Więc dobrym podejściem do tego jest ustawienie tych wartości na wysokie liczby na początek, żeby Twój pipeline działał bez błędów, a następnie poprosić Nextflow'a o wygenerowanie raportu o użyciu dla Ciebie.

To jest super łatwe do zrobienia, więc wrócę do terminala. O, muszę pamiętać, żeby ustawić to z powrotem na local, żeby mój pipeline faktycznie działał. I powiem Nextflow run, i użyję flagi linii poleceń -with-report.

I mogę zostawić to puste i da domyślną nazwę pliku, ale dam mu konkretną nazwę pliku, żeby zostało zapisane w konkretnym miejscu.

Naciśnij Enter, i pipeline działa dokładnie jak normalnie, ale kiedy się skończy, wygeneruje ładny raport HTML dla mnie.

Więc na pasku bocznym tutaj mam ten plik HTML. Gdybym uruchamiał to lokalnie, po prostu bym go otworzył. Ponieważ jestem w Codespaces, kliknę prawym przyciskiem myszy i kliknę pobierz, co pobierze go na mój lokalny komputer. I mogę go łatwo otworzyć w przeglądarce internetowej.

Nextflow może wygenerować taki raport dla każdego pipeline'u i ma kilka naprawdę fajnych informacji. Więc dobrą praktyką jest zawsze zapisywać te rzeczy. Mówi nam, kiedy uruchomiliśmy, gdzie uruchomiliśmy, czy to się powiodło czy nie, jakie parametry zostały użyte, jakie było polecenie CLI, takie rzeczy.

I są też te wykresy dotyczące użycia zasobów. Więc mówi nam, jaki procent wywołań CPU został użyty dla każdego procesu jako wykres pudełkowy tutaj, ponieważ jest wiele zadań dla każdego procesu, więc możemy zobaczyć dystrybucję.

Widzisz nasze procesy tutaj, cowpy i collectGreetings miały tylko pojedyncze zadanie, więc to jest tylko pojedyncza linia. I mamy zarówno CPU i pamięć i czas trwania zadania, i były bardzo szybkie.

Jeśli używasz Seqera Platform, nawiasem mówiąc, dostajesz te same wykresy wbudowane w interfejs Platform bez konieczności robienia czegokolwiek. Więc zawsze masz te informacje pod ręką.

Dobra, więc możemy użyć tego raportu i na prawdziwym uruchomieniu, i poczuć, ile CPU i ile pamięci jest używane przez nasz pipeline i wrócić i umieścić te wartości z powrotem w naszym pliku konfiguracyjnym, tak żeby następnym razem może nie prosiliśmy o aż tyle. I możemy być trochę bardziej oszczędni.

Teraz możesz być naprawdę sprytny w konfigurowaniu plików konfiguracyjnych pipeline'u. I znowu, jeśli używasz Seqera Platform, szukaj małego przycisku, który wygląda jak żarówka. Bo jeśli go klikniesz, wygeneruje wysoce zoptymalizowany plik konfiguracyjny, który jest dostosowany specjalnie dla Twoich danych, Twojego uruchomienia i Twojego pipeline'u. Aby uruchomić go w najbardziej efektywny sposób.

Ale na razie powiem, że właściwie domyślna liczba CPU, którą Nextflow dawał, była w porządku, ale potrzebujemy tylko jednego gigabajta pamięci.

## 5.3. Ustaw alokację zasobów dla konkretnego procesu

Teraz w prawdziwym życiu dość niezwykłe jest, że wszystkie procesy w Twoim pipeline'ie będą potrzebować tych samych wymagań. Możesz mieć coś takiego jak MultiQC jako narzędzie raportowania, które potrzebuje bardzo mało zasobów i działa dość szybko.

A potem może masz coś, co indeksuje genom referencyjny lub wykonuje jakieś dopasowanie lub wykonuje jakąś inną pracę. Nieważne, co to jest, co wymaga dużo zasobów. I więc dla tych różnych przesłań zadań do harmonogramu chcesz dać różne ilości zasobów.

Pod tym zakresem process możemy zdefiniować konfigurację, która kieruje konkretne procesy na różne sposoby.

Tutaj używamy withName, możemy też używać etykiet, i mogą one używać wzorca do kierowania jednego lub wielu procesów. Tutaj po prostu mówimy, że dowolne procesy, które mają nazwę cowpy, ustawione na dwa gigabajty pamięci i dwa CPU, i ponieważ jest to bardziej specyficzny selektor niż na najwyższym poziomie process, jest to nadpisane w tych przypadkach, więc możesz zbudować ładny plik konfiguracyjny tutaj, który naprawdę dostosowuje wszystkie Twoje różne procesy w Twoim pipeline'ie, aby uczynić je naprawdę efektywnymi.

## 5.5. Dodaj limity zasobów

Teraz jako programista pipeline'u prawdopodobnie dość dobrze znam narzędzia i chcę, żeby wszystko działało tak szybko i tak dobrze, jak to możliwe. Więc może być tak, że umieszczam dość wysokie liczby dla niektórych z nich, ponieważ wiem, że będzie działać znacznie szybciej, jeśli dam cowpy 20 CPU zamiast dwóch.

To jest w porządku, dopóki nie przejdziesz do uruchamiania na swoim laptopie lub na testach ciągłej integracji GitHub Actions, lub jakimś innym systemie, który może nie mieć dostępnych 20 CPU.

Teraz kiedy próbujesz uruchomić pipeline, zawiedzie, ponieważ Nextflow powie, nie mogę przesłać tego zadania nigdzie. Nie mam dostępnych zasobów.

Teraz, aby uniknąć tego twardego błędu, możemy dodać trochę więcej konfiguracji, która jest specyficzna dla naszego systemu teraz, nazywanej limitami zasobów. I wygląda to tak. Jest pod zakresem process znowu.

I limity zasobów, możesz określić zasadniczo sufit tego, co masz dostępne. To mapa tutaj, i możesz, w tej mapie, możesz ustawić pamięć, CPU i czas.

Teraz to, co się dzieje, to gdy Nextflow przesyła zadanie z procesu, patrzy na to, co jest żądane i zasadniczo po prostu robi minimum między tym a tym. Więc jeśli zażądaliśmy 20 CPU, ale tylko cztery są dostępne, zażąda czterech. Pipeline nie zawiedzie i użyje możliwie najbliżej tego, co zostało zaprojektowane przez programistę pipeline'u.

## 6. Użyj profili, aby przełączać się między predefiniowanymi konfiguracjami

Dobra. Powiedziałem, że limity zasobów tutaj mogą być specyficzne dla systemu, i może mam plik Nextflow config w moim pipeline'ie, i wiem, że ludzie będą używać tego w szeregu różnych miejsc. Teraz, zamiast zmuszać wszystkich do tworzenia własnego pliku Nextflow config za każdym razem, to, co mogę zrobić, to pogrupować różne presety konfiguracji razem w profile konfiguracyjne.

Przewinę tutaj trochę w dół i tuż za params, ponieważ kolejność pliku konfiguracyjnego tutaj jest ważna, plik konfiguracyjny jest ładowany sekwencyjnie, więc umieszczę te profile po wszystkim innym, aby nadpisały wcześniej zdefiniowane params. I wkleję te profile z materiałów szkoleniowych.

Więc jest nowy zakres najwyższego poziomu o nazwie profiles. Możemy mieć dowolne nazwy tutaj. Więc mamy my_laptop i univ_hpc. I tutaj widzimy, że ustawiamy te same parametry konfiguracyjne, które mieliśmy wcześniej. Teraz tylko w profilu. Więc mamy lokalny executor do uruchamiania na moim laptopie i przesyłam do klastra SLURM na HPC.

Używam Docker'a lokalnie, condy na HPC, i system HPC ma znacznie wyższe limity zasobów.

Teraz mogę uruchomić pipeline z opcją CLI -profile, powiedzieć, którego profilu chcę użyć. Więc użyję my_laptop, i Nextflow zastosuje całą konfigurację w tym zakresie profilu. Więc mogę teraz spróbować. To jest to samo polecenie co wcześniej. Nextflow run hello config, i robię dash profile, pojedynczy myślnik, bo to podstawowa opcja Nextflow'a, dash profile my_laptop.

Teraz zamierza wsadowo zastosować całą tę opcję konfiguracyjną. O, i widzisz, mówiłem wcześniej, że to może się zdarzyć, że wymaganie procesu, prosił o cztery CPU, a mam tylko dwa na tej instancji Codespaces.

Więc to dobra okazja, żeby wypróbować limity zasobów procesu i powiedzieć, że mam tylko dwa CPU na moim laptopie, lub w tym Codespaces. Teraz jeśli uruchomimy to ponownie, powinno ograniczyć to wymaganie do dwóch i miejmy nadzieję, że pipeline się uruchomi. Świetnie.

## 6.2. Utwórz profil parametrów testowych

Zauważ, że te profile nie muszą zawierać tylko konfiguracji dotyczącej ich infrastruktury. Możesz mieć grupowania dowolnej konfiguracji tutaj, w tym parametrów.

Więc kolejną rzeczą, którą często zobaczysz w pipeline'ach ludzi, jest profil testowy, który zawiera parametry, które normalnie przesłałbyś na podstawie użytkownika. Ale tutaj mamy zasadniczo różne rozsądne wartości domyślne na wypadek, gdy chcę uruchomić przypadki testowe.

I to jest świetne, ponieważ nie muszę koniecznie iść i określać wszystkich tych rzeczy, które mogą być wymaganymi parametrami. W przeciwnym razie mogę po prostu powiedzieć dash profile test i po prostu uruchomi się od razu.

Teraz coś, co warto zauważyć, to to, że profile mogą być również łączone więcej niż jeden. Więc mogę zrobić profile my_laptop tutaj, a następnie również dodać test. Nie robię profile dwa razy. Po prostu robię listę oddzieloną przecinkami bez spacji. I zamierza zastosować te profile po kolei. Więc weźmie konfigurację z profilu my_laptop, a następnie zastosuje konfigurację test na wierzchu.

Naprawdę wygodne i widzisz, jak możesz ustawić wiele rozsądnych domyślnych grup tutaj, aby ułatwić uruchamianie Twojego pipeline'u.

## 6.3. Użyj nextflow config, aby zobaczyć rozwiązaną konfigurację

Miejmy nadzieję, przekonałem Cię, że rozwiązywanie konfiguracji Nextflow'a jest potężne, ale nie miałbym Ci za złe, gdybyś trochę mętniał w oczach w tym momencie po tym, jak powiedziałem około 20 różnych sposobów na dostarczenie konfiguracji i dałem wszystkie te różne warstwy jak skórkę cebuli.

Więc jeśli kiedykolwiek czujesz się niepewny co do tego, jaka jest ostateczna rozwiązana konfiguracja dla Nextflow'a, wiedz, że jest polecenie o nazwie "nextflow config", i możemy je uruchomić i powie nam, jaka jest rozwiązana konfiguracja w naszej obecnej lokalizacji.

Więc kiedy uruchamiam to tutaj, znajduje plik "nextflow.config" w bieżącym katalogu roboczym i przetwarza całą różną konfigurację i daje mi rozwiązane wyjście.

Zauważ, że plik konfiguracyjny Nextflow może również przyjąć opcję CLI profile. Więc jeśli powiem mu rozwiązać w profilach my_laptop i test, i widzisz, że również zastosował limity zasobów tutaj z opcji konfiguracji my_laptop i również ustawił params, które były w test.

Więc to fajny sposób po prostu na zbadanie, jak działa rozwiązywanie konfiguracji, jeśli jesteś w ogóle niepewny.

## Podsumowanie

Dobra, to wszystko. To jest konfiguracja Nextflow'a w pigułce. Możesz zrobić wiele rzeczy z konfiguracją. Jest naprawdę potężna. Ale to są większość typowych przypadków użycia, które będziesz robić, i te koncepcje stosują się do wszystkich różnych opcji.

Poklepie się po plecach, bo to koniec kursu szkoleniowego Hello Nextflow. Masz teraz, miejmy nadzieję, pewność zarówno w pisaniu własnego pipeline'u Nextflow'a od zera, konfigurowaniu go i uruchamianiu, i znasz wszystkie zawiłości i rzeczy, na które należy uważać.

Jest jeszcze jeden quiz, który możesz wypróbować na stronie szkoleniowej dotyczącej konfiguracji. Więc zejdź w dół i wypróbuj to i upewnij się, że zrozumiałeś wszystkie te części dotyczące konfiguracji.

I dołącz do nas w ostatnim wideo tylko na szybkie podsumowanie niektórych kolejnych kroków, które mogłyby być dobre do zrobienia po tym kursie szkoleniowym.

Dziękuję, że z nami wytrzymałeś. Dobra robota i do zobaczenia w następnym wideo.
