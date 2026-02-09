# Część 6: Hello Config - Transkrypcja wideo

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Tłumaczenie wspomagane przez AI - [dowiedz się więcej i zaproponuj poprawki](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

<div class="video-wrapper">
  <iframe width="560" height="315" src="https://www.youtube.com/embed/FcZTiE25TeA?si=tnXTi6mRkITY0zW_&amp;list=PLPZ8WHdZGxmWKozQuzr27jyMGqp9kElVK&amp;cc_load_policy=1" title="YouTube video player" frameborder="0" allow="accelerometer; autoplay; clipboard-write; encrypted-media; gyroscope; picture-in-picture; web-share" referrerpolicy="strict-origin-when-cross-origin" allowfullscreen></iframe>
</div>

!!!note "Ważne uwagi"

    Ta strona zawiera tylko transkrypcję. Pełne instrukcje krok po kroku znajdziesz w [materiałach szkoleniowych](../06_hello_config.md).

    Numery sekcji pokazane w transkrypcji mają charakter orientacyjny i mogą nie obejmować wszystkich numerów sekcji w materiałach.

## Powitanie

Cześć i witaj ponownie w części szóstej Hello Nextflow. Ta sekcja dotyczy plików konfiguracyjnych i jest ostatnią częścią tego kursu.

Nextflow szczególnie dobrze radzi sobie z dwoma rzeczami: odtwarzalnością i przenośnością. Pliki konfiguracyjne to miejsce, gdzie naprawdę widać tę drugą cechę. Możliwość skonfigurowania pipeline'u Nextflow'a do działania na różne sposoby i w różnych systemach, bez konieczności edytowania podstawowego kodu pipeline'u.

Ta supermoc pozwala na ponowne wykorzystanie pipeline'ów Nextflow'a przez inne osoby w różnych miejscach lub w różnych infrastrukturach, do których Ty sam możesz mieć dostęp.

Oznacza to, że możesz rozwijać kod pipeline'u na swoim laptopie, przesłać go do chmury, uruchomić na swoim HPC i to jest ten sam kod pipeline'u, który działa wszędzie.

W tej sekcji przejdziemy przez kilka tematów. Zaczniemy od tego, jak Nextflow obsługuje pliki konfiguracyjne, skąd je ładuje, jak je pisać i jak je strukturyzować, oraz tego rozdzielenia między samym pipeline'em a tym, co powinno znaleźć się w pliku konfiguracyjnym.

Następnie przejdziemy do kilku typowych przypadków użycia, takich jak zmiana miejsca przechowywania plików wyjściowych, a także jak sprawić, by pipeline działał w różnych infrastrukturach, zarówno przy użyciu różnych typów pakowania oprogramowania, jak i przesyłania zadań do różnych infrastruktur.

## Hierarchie plików konfiguracyjnych

Dobra, zaczynajmy. Jeśli chodzi o ładowanie plików konfiguracyjnych, Nextflow może pobierać je z wielu różnych miejsc, co jest dobrą rzeczą, ale może być też nieco ryzykowne, ponieważ czasami może być trochę trudno wiedzieć, skąd pobiera plik konfiguracyjny i w jakiej kolejności ładuje rzeczy.

Dlatego naprawdę polecam kliknięcie tego linku, który prowadzi do dokumentacji Nextflow'a. Na tej stronie konfiguracji wymienione są kluczowe miejsca, z których ładowana jest konfiguracja, i co ważne, kolejność, w jakiej te rzeczy są ładowane.

Możesz więc zobaczyć, że możesz umieścić plik konfiguracyjny w swoim katalogu domowym Nextflow'a, który zazwyczaj to ".nextflow" w Twoim katalogu domowym. I ten plik będzie zawsze ładowany przez każde uruchomienie Nextflow'a w Twoim systemie.

Następnym miejscem do sprawdzenia jest plik w katalogu głównym Twojego repozytorium lub katalogu pipeline'u o nazwie "nextflow.config".

Następnie kolejny plik o nazwie "nextflow.config", ale tym razem w katalogu, z którego uruchamiasz Nextflow'a: katalogu uruchomienia.

Na koniec możesz podać ścieżki do plików konfiguracyjnych w wierszu poleceń za pomocą argumentu "-c" i możesz to zrobić wiele razy. Są one stosowane w kolejności, w jakiej je określisz.

Możesz podać pliki konfiguracyjne we wszystkich tych lokalizacjach, jeśli chcesz, i będą one ładowane iteracyjnie, każdy nadpisując poprzedni tylko w tych zakresach konfiguracji, gdzie się pokrywają.

To naprawdę potężny system, ponieważ oznacza, że możesz ustawić rozsądne wartości domyślne, a następnie stopniowo być coraz bardziej szczegółowym, zawężając tę konfigurację.

## 0. Rozgrzewka: Uruchom hello-config.nf

Dobra, zamknijmy to i przejdźmy do naszego Codespaces i zacznijmy. Jak poprzednio, posprzątałem tutaj, usunąłem moje poprzednie katalogi wyników, mój Nextflow i katalogi robocze i tak dalej. Nie martw się, jeśli nadal masz te pliki. To tylko dlatego, że jestem bardzo przybliżony i w przeciwnym razie rzeczy szybko się zaśmiecają.

Będziemy pracować z hello-config.nf, ostatnim plikiem w naszym katalogu, i powinien on nawiązywać do miejsca, w którym skończyliśmy w poprzedniej sekcji.

Mamy więc nasze cztery różne procesy, które są dołączane z plików modułów. Mamy nasze parametry pipeline'u, nasz blok `workflow`, w którym wywołujemy różne procesy i łączymy kanały, publikujemy kanały wyjściowe, a następnie blok `output` na dole, gdzie definiujemy, gdzie te pliki powinny być przechowywane i jak powinny być kopiowane.

Mamy również już plik "nextflow.config" z ostatniego rozdziału, w którym włączamy Dockera, i będziemy rozbudowywać ten plik dzisiaj.

Jak poprzednio, zmieniliśmy ścieżkę wyjściową w tym głównym skrypcie na hello config, żeby nie kolidowała z poprzednimi wynikami, które wygenerowałeś.

Dobra, sprawdźmy szybko, czy wszystko nadal działa zgodnie z oczekiwaniami. Otwórzmy terminal i wykonajmy nextflow run hello-config.nf. Nextflow się ładuje. Powinien uruchomić nasze cztery różne procesy. Wygenerować trochę ładnej grafiki ASCII przy użyciu cowpy, a następnie zapisać nasze wyniki do naszych plików wyników w tym katalogu.

Mogę szybko zajrzeć tutaj, żeby upewnić się, że te pliki wyglądają zgodnie z oczekiwaniami, i rzeczywiście, oto nasz gigantyczny indyk. Świetnie.

## 1.1. Przenieś wartości domyślne do nextflow.config

Teraz pierwszą rzeczą, którą zrobimy, jest przeniesienie niektórych rzeczy z naszego skryptu do naszego pliku konfiguracyjnego.

A tym, na czym nam zależy, są głównie parametry na tym etapie. Chcemy przenieść wartości domyślne do pliku konfiguracyjnego, żeby było jaśniejsze, jakie są wartości domyślne i żeby łatwiej było je nadpisać.

Wezmę ten blok `params` ze skryptu i umieszczę go w pliku konfiguracyjnym. I musimy być tutaj trochę ostrożni, ponieważ w tej chwili składnia jest nieco inna między konfiguracją a skryptami. Plik konfiguracyjny nie może przyjmować deklaracji typów, ponieważ tak naprawdę nie definiujemy tych parametrów, tylko się do nich odwołujemy. Więc pozbędę się ich.

Ale poza tym jest bardzo podobnie. Mamy blok `params`, a następnie mamy nasze różne parametry wejściowe, parametr batch, parametr character.

Mogę teraz wrócić do mojego skryptu i nie muszę już definiować tych wartości domyślnych, ponieważ te wartości są teraz w moim pliku Nextflow config.

Jednak zostawiam nazwy parametrów i ich typy, aby Nextflow znał te informacje i mógł nadal wykonywać wszystkie sprawdzenia bezpieczeństwa typów i wszystko inne.

Dobra. Zapiszmy te pliki i szybko sprawdźmy, czy wszystko nadal działa tak samo jak wcześniej. Nie powinno być żadnych zmian. Zachowaliśmy te same wartości. Po prostu przenieśliśmy miejsce, w którym zostały zdefiniowane.

Świetnie.

## 1.2. Użyj pliku konfiguracyjnego specyficznego dla uruchomienia

Teraz do tej pory uruchamialiśmy Nextflow'a z tego samego katalogu, w którym mamy nasz skrypt pipeline'u. Więc nasz katalog uruchomienia i katalog pipeline'u to w zasadzie to samo.

Aby pokazać, jak możemy mieć różne pliki konfiguracyjne z różnymi katalogami uruchomienia, utworzymy teraz nowy podkatalog.

Więc powiem mkdir i nazwiemy go tux-run.

A następnie wykonam cd, zmienię katalog na tux-run. I zauważ, że jesteśmy teraz w naszym katalogu roboczym, który nie jest już w tym samym katalogu co skrypty pipeline'u.

Dobra, stwórzmy nowy plik "nextflow.config". Więc touch nextflow config i otwórzmy go w VS Code. Możesz również zobaczyć na pasku bocznym tutaj, że jesteśmy teraz w tym podkatalogu.

Teraz możemy wziąć ten sam blok `params`, który mieliśmy w nextflow.config najwyższego poziomu, skopiować go tutaj i teraz możemy zmienić te wartości.

Po pierwsze, dane są teraz inną ścieżką względną, ponieważ jesteśmy w podkatalogu, więc musimy to zaktualizować. A następnie zmienimy batch na experiment i zmienimy character z Turkey na tux.

Teraz kliknijmy zapisz i wypróbujmy to. Tak jak w przypadku danych, muszę teraz powiedzieć ../ aby dostać się do skryptu. Więc to Hello config. I naciskam enter.

Kod pipeline'u w ogóle się nie zmienił, ale teraz będziemy mieli dwa zestawy ładowania konfiguracji, a plik konfiguracyjny katalogu uruchomienia powinien nadpisać wartości domyślne, które zostały ustawione w nextflow.config pipeline'u, i powinniśmy otrzymać różne zestawy wyników.

Rzeczywiście, w naszym katalogu tutaj, w tux-run, możesz zobaczyć, że mamy katalog dot Nextflow i katalog work i to dlatego, że są one zawsze tworzone w Twoim katalogu uruchomienia. Więc są one inne niż katalogi work i results, które mieliśmy z wcześniejszych uruchomień.

Teraz, jeśli spojrzę w results, możemy zobaczyć nasz collected i tam jest nasza mała postać tux. Więc widać, że te parametry zostały prawidłowo zinterpretowane.

## 1.3. Użyj pliku parametrów

Dobra. Wcześniej, kiedy mówiłem o różnych plikach konfiguracyjnych, które mogą być ładowane, pominąłem jedno inne miejsce, z którego możemy pobrać konfigurację.

Możesz ją pobrać z wiersza poleceń, jak widzieliśmy z podwójnym myślnikiem i nazwami parametrów, ale możemy również dostarczyć plik YAML lub JSON, tylko z parametrami.

Plik konfiguracyjny może mieć wszystkie różne typy zakresów, ale te pliki to tylko parametry i jest to przyjazny dla użytkownika sposób dostarczania wielu parametrów naraz, a być może nieco bardziej odtwarzalny sposób, ponieważ zapisujesz je do pliku, więc łatwo je uzyskać na późniejszym etapie.

Więc wróćmy do naszego terminala i zanim zapomnimy, upewnijmy się, że wrócimy o katalog wyżej, więc nie jestem już w podkatalogu, i przyjrzę się plikowi YAML, który mamy tutaj o nazwie test-params.yaml.

Więc jeśli po prostu zrobię code test-params.yaml, możesz zobaczyć, że to jest po prostu zwykły plik YAML. Nic specjalnego. Z kluczami będącymi naszymi nazwami parametrów, z formatowaniem YAML, więc dwukropek tutaj, a następnie wartość.

Zauważ, że to nie jest kod Nextflow'a, więc nie możemy umieszczać tutaj rzeczy takich jak zmienne. To są tylko wartości statyczne.

Również ponieważ JSON faktycznie parsuje się jako YAML, możemy również mieć plik test-params.json, który wygląda bardzo podobnie. To po prostu inny format danych.

Więc mamy tutaj dwa różne pliki testowe i mamy nieco różne zmienne.

Dobra, więc jak przekazujemy je do Nextflow'a? To bardzo proste. Robimy Nextflow run hello config, jak poprzednio. I zamiast "-c" dla pliku konfiguracyjnego lub ładowania tych domyślnych nazw plików, robimy -params-file. Pojedynczy myślnik, ponieważ jest to podstawowa opcja Nextflow'a.

A następnie przekazujemy ścieżkę do tego pliku. Więc zrobię "-params-file test-params.yaml" i zobaczymy, czy są one prawidłowo załadowane.

Dobra. Uruchomił się. Przypomnijmy sobie, co było w tym pliku YAML. Więc batch był ustawiony na YAML, więc tak powinien być nazwany, i powinien mieć stegozaura. Więc przejdźmy w górę i spójrzmy w results. I mamy COLLECTED-yaml. Więc zobaczmy, czy mamy stegozaura. Fantastycznie, stegozaur w kapeluszu. To właśnie lubimy.

Więc to zadziałało naprawdę dobrze i jest dokładnie tak samo z plikiem JSON. Po prostu zamieniamy tutaj rozszerzenie pliku i Nextflow wie, jak to odczytać.

I w tym przypadku powinniśmy mieć batch o nazwie JSON i powinniśmy mieć żółwia. Spójrzmy. Wspaniale. Jedno z moich ulubionych narzędzi CLI.

## 2.1. Dostosuj katalog wyjściowy za pomocą -output-dir

Dobra, więc to było głównie myślenie o wejściach do pipeline'u i zmianie parametrów. A co z wyjściami?

Teraz, chociaż zmienialiśmy podkatalogi używając parametrów, mogłeś zauważyć, że wszystkie nasze pliki nadal trafiają do results.

Możemy zmienić ten katalog bazowy, do którego publikowane są wszystkie pliki, za pomocą flagi wiersza poleceń o nazwie -output-dir. Więc jeśli zrobię Nextflow run hello config, a następnie zrobię -output-dir i nazwiemy to "custom-outdir-cli". Nie mogę pisać. Po prostu żebyśmy pamiętali, skąd pochodzą te pliki.

To jest podstawowa opcja Nextflow'a i jest bardzo nowa. Została dodana dopiero niedawno i jest to jedna z rzeczy, które możemy zrobić z nowym parserem języka i wszystkim.

To trochę długo się pisze. Możesz również po prostu nazwać to "-o", jeśli chcesz. Więc jeśli po prostu wrócę. Mogę po prostu skrócić to do "-o", co jest nieco prostsze.

Dobra. Uruchamiamy to. Nie zmieniliśmy nic w naszym pipeline'ie ani nawet w naszej konfiguracji na tym etapie, i powinno to miejmy nadzieję zapisać wszystkie nasze wyniki do innego katalogu najwyższego poziomu. I możesz sobie wyobrazić, że możesz ustawić to na zasadniczo dowolną ścieżkę, którą chcesz.

Właśnie pojawił się na górze. Mamy custom-outdir-cli i wszystkie pliki są tam zorganizowane w dokładnie taki sam sposób, z tymi samymi podkatalogami i nazwami plików. Więc to naprawdę łatwy sposób, aby po prostu zmienić miejsce, w którym pipeline publikuje swoje wyniki, bez zastanawiania się zbytnio nad tym, jak te wyniki są zorganizowane.

## 2.1.2. Usuń zakodowane na stałe ścieżki z bloku output

Jeśli spojrzę do tego katalogu, możemy zobaczyć, że nadal mamy podkatalog o nazwie Hello Config, co wydaje się teraz trochę zbędne.

Więc po prostu załadujmy ponownie nasz skrypt i możemy teraz usunąć ten podkatalog z bloku `output` na dole. Ponieważ tak naprawdę już go nie potrzebujemy. Więc możemy to teraz zrobić, usunąć to stąd. A jeśli to jest tylko to, możesz albo całkowicie to usunąć, albo zostawić jako pusty ciąg znaków. Zostawię to jako pusty ciąg znaków na razie, ponieważ wrócimy i umieścimy tam różne rzeczy w przyszłości. Ale jeśli nie zależy Ci na podkatalogach, najczystsze jest po prostu całkowite usunięcie deklaracji path.

Dobra, zapiszmy. Szybko spróbujmy tego ponownie. Właściwie usunę mój katalog "custom-outdir-cli", żebyśmy nie byli zdezorientowani przez jakiekolwiek istniejące tam pliki. Ponieważ pamiętaj, kiedy publikujesz rzeczy, nie usuwa to plików, które tam już były. Po prostu dodaje nowe. Uruchommy to polecenie ponownie, custom-outdir-cli.

A teraz jeśli zrobisz "ls custom-outdir-cli", nie ma już tam katalogu o nazwie Hello Config.

## 2.2.1. Ustaw outputDir w pliku konfiguracyjnym

Dobra, flaga wiersza poleceń tutaj, "-o" lub "-output-dir" jest dobra. Ale co z ustawianiem wartości domyślnych dla tego w konfiguracji? Jak to robimy?

Otwieram plik "nextflow.config", zamykam wszystko inne i pozbywam się tego. Możemy dodać tutaj nową opcję konfiguracyjną, którą właśnie skopiowałem ze strony materiałów szkoleniowych, i nazywa się outputDir.

Nie jest pod żadnym zakresem. Nie jest pod params ani niczym innym. Jest na najwyższym poziomie i możemy ustawić to na ciąg znaków. Teraz prostą rzeczą do zrobienia jest po prostu zmienić to na cokolwiek innego niż results jako zakodowany na stałe ciąg znaków. Ale ponieważ to jest w pliku konfiguracyjnym Nextflow'a, możemy być tutaj trochę sprytni i również uwzględnić zmienne.

I możesz zobaczyć tutaj, że uwzględniliśmy zmienną params, params.batch, która jest częścią tego ciągu znaków. Oznacza to, że możemy ponownie wykorzystać zmienne, które pochodzą z innych miejsc. I w tym przypadku, jeśli zrobimy --batch, kiedy uruchamiamy Nextflow Pipeline, otrzymamy podkatalog w naszej niestandardowej ścieżce na podstawie tego, jaka była nazwa batch.

Dobra, więc wypróbujmy to i po prostu szybko zobaczmy, jak wyglądają wyniki. Więc jeśli zrobię Nextflow run hello config i --batch my_run. Przypomnijmy sobie, jak wyglądała konfiguracja. Więc to custom-outdir-config.

Tree custom-outdir-config. I możesz zobaczyć, że batch nazywał się my_run. A następnie mamy ten podkatalog o nazwie my_run. Więc ta dynamiczna ścieżka pliku zadziałała.

I nie tylko to, nie trafiło już do domyślnego katalogu results, i nie musiałem określać niczego w wierszu poleceń, aby zmienić katalog bazowy. Więc pomyślnie zresetowaliśmy wartość domyślną dla domyślnego outputDir.

## 2.2.2. Podkatalogi z nazwami batch i procesów

Dobra, pójdźmy z tym trochę dalej. To jest dynamiczna zmienna w pliku konfiguracyjnym. A co ze skryptem? Teraz do tej pory mieliśmy te ścieżki tutaj i one również mogą być dynamiczne. Więc zamiast po prostu kodować coś na stałe, możemy umieścić nawiasy klamrowe i umieścić coś dynamicznego.

Na przykład mamy nasze procesy o nazwie sayHello. Moglibyśmy zrobić sayHello.name, który jest atrybutem procesu, co jest trochę nudne, ponieważ to po prostu "sayHello" w tym przypadku. Ale jest zmienne.

Więc to daje Ci pomysł. Więc możemy to tutaj umieścić i powiedzieć convertToUpper.name, collectGreetings.name, collectGreetings.name ponownie i cowpy.

Teraz kiedy uruchomimy, katalog bazowy nadal będzie custom-outdir-config. I będzie w podkatalogu o nazwie params.batch, ale podkatalogi pod tym powinny być zorganizowane według nazwy procesu.

Wypróbujmy to i zobaczmy, czy to działa. Więc usunę poprzedni katalog, żebyśmy się nie pomylili, i po prostu użyję dokładnie tego samego polecenia Nextflow Run.

Powinno działać w ten sam sposób. Mógłbym używać dash resume na wszystkich tych, aby było trochę szybciej i używać wcześniej obliczonych wyników. Teraz, jeśli zrobię tree custom-outdir-config, możesz zobaczyć, że nie jest w results, jest w naszym katalogu bazowym z nazwą batch. I możesz zobaczyć, że wszystkie wyniki są teraz zorganizowane w podkatalogach nazwanych według procesu. Więc mamy dwa różne miejsca, w których definiujemy dynamiczne ścieżki wyjściowe tutaj.

Dobra. Ostatnia rzecz, dodajmy z powrotem te foldery pośrednie, które mieliśmy wcześniej, ponieważ były całkiem fajne. Intermediates.

I możemy również pomyśleć trochę o tym params.batch, może jako deweloper pipeline'u naprawdę lubiłem mieć to w podkatalogu, ale jeśli użytkownicy końcowi pipeline'u ustawiają "-o" lub -output-dir w CLI, całkowicie nadpisuje to całe wyrażenie i tracimy ten podkatalog.

Więc to, co możemy zrobić, to możemy wyjąć tę dynamiczną ścieżkę z konfiguracji outputDir, która byłaby nadpisana, i umieścić ją w ścieżce wyjściowej, która nie jest nadpisywana.

Więc możemy zrobić params.batch slash intermediates slash sayHello.name i zrobić to wszystko w ciągu znaków w podwójnych cudzysłowach, więc jest interpolowane przez Nextflow'a.

Mogę teraz skopiować, ups. Skopiować to do innych procesów. Pamiętaj, aby umieścić je wszystkie w cudzysłowach. I usunąć intermediates z tych konkretnych wyjść.

Dobra? Wygląda to teraz nieco bardziej skomplikowanie, ale możesz zobaczyć, że naprawdę zaczynamy budować ładnie zorganizowaną strukturę katalogów wyjściowych w naszym kodzie.

I co jest naprawdę fajne, to to, że ta dodatkowa złożoność w kodzie nie przechodzi do CLI. Więc możemy uruchomić nasze polecenie z -output-dir i jakimikolwiek zmiennymi batch, po prostu myśląc o tym, jak uruchomić pipeline i nie myśląc zbytnio o tym, co jest w kodzie. A nasze pliki wyjściowe będą skonstruowane naprawdę ładnie w bardzo dobrze zorganizowany sposób, co jest miłe dla osób korzystających z pipeline'u w zasadzie.

Świetnie. Kiedy to piszę, zdaję sobie sprawę, że popełniłem błąd. Zobaczmy, czy ktoś mnie złapał tutaj. Mamy collectGreetings.name, więc coś poszło nie tak. I tak, rzeczywiście, przypadkowo zapomniałem umieścić to w nawiasach klamrowych.

Więc pamiętaj, bądź ostrożny, kiedy piszesz swój kod i upewnij się, że mówisz Nextflow'owi, co jest zmienną, a co jest po prostu ciągiem znaków. Ponieważ zrobi dokładnie to, co mu powiesz. I nic więcej. Jak wszystkie dobre komputery. Dobra, to powinno to naprawić.

## 2.3. Ustaw tryb publikowania na poziomie workflow'u

Jest jedna część tego skryptu, której nadal nie lubię, a mianowicie fakt, że piszemy mode copy w kółko, a jeśli jest jedna rzecz, której nie lubimy, to powtarzanie się.

Więc możemy to trochę uporządkować, biorąc to i przenosząc do konfiguracji. I w rzeczywistości możemy ustawić to dla całego pipeline'u za jednym razem. Więc nie musimy mówić tego wiele razy.

Przechodzimy do naszego pliku konfiguracyjnego i mamy tutaj nowy zakres o nazwie workflow. I możemy albo zrobić nawiasy klamrowe, albo możemy użyć notacji kropkowej. Nie ma to żadnego znaczenia. Strona materiałów szkoleniowych używa notacji kropkowej. Mogę powiedzieć output i możemy mieszać i dopasowywać, więc mode equals copy. Świetnie.

I teraz możemy wrócić tutaj i usunąć te. Teraz moglibyśmy je zostawić na miejscu. Konfiguracja zasadniczo nadpisuje to, co jest tutaj napisane, ale ponieważ mamy to w konfiguracji na poziomie pipeline'u, a te dwa pliki są dostarczane razem, nie ma powodu, aby naprawdę robić to dwa razy.

Dobra. Po prostu sprawdźmy się, ponieważ najwyraźniej popełniamy błędy. Uruchommy to ponownie i po prostu sprawdźmy, czy prawidłowo używamy trybu copy do publikowania plików. Więc uruchomimy skrypt ponownie i tym razem umieściliśmy wyniki w katalogu o nazwie config-output-mode, zobaczmy, jak wyglądają tam pliki.

A następnie jeśli zrobię "ls -l", aby spojrzeć na batch, i możemy spojrzeć na cowpy, na przykład. I powinniśmy zobaczyć, tak, że to jest właściwy plik tutaj, który nie jest dowiązaniem symbolicznym, więc ten atrybut konfiguracji został prawidłowo zastosowany.

## 3. Wybierz technologię pakowania oprogramowania

Dobra. Do tej pory skupialiśmy się na wejściach i wyjściach, plikach, z którymi pracuje workflow. Ale co z infrastrukturą? Powiedziałem na początku, że Nextflow pozwala uruchomić ten sam pipeline w różnych konfiguracjach obliczeniowych. Więc jak to wygląda?

Aby to pokazać, przełączymy się z używania Dockera do uruchamiania cowpy i zamiast tego użyjemy Condy, aby zrobić to samo.

Mogę to zrobić bardzo prosto. Jeśli przejdę do code, "nextflow.config". Jeśli pamiętasz na górze, zdefiniowaliśmy docker.enabled wcześniej, w ostatnim rozdziale, abyśmy mogli użyć kontenera z cowpy.

Powiem Nextflow'owi, aby nie używał Dockera. Ustawię to na false. I powiem Conda enabled equals true. Więc powiem Nextflow'owi, proszę użyj Condy.

Teraz samo włączenie Condy nie wystarczy. Tak jak zrobiliśmy z Dockerem, musimy powiedzieć Nextflow'owi, skąd może pobrać potrzebne oprogramowanie.

Więc jeśli przejdziemy do modułów tutaj. I otworzymy skrypt cowpy. Możemy zobaczyć, że mamy deklarację container na górze. A kontener jest używany przez Dockera, ale także Singularity, Apptainer i wiele innych narzędzi programowych.

Ale nie może być używany dla Condy, więc mamy osobną deklarację o nazwie "conda" i moglibyśmy po prostu napisać "cowpy". I to pozostawi to rozwiązaniu pakietów conda, aby dowiedzieć się, jaki jest najlepszy sposób na rozwiązanie tego, zgodnie z Twoim lokalnym środowiskiem conda.

Lub dobrą praktyką jest zrobienie tego, co mówi strona materiałów szkoleniowych, czyli zdefiniowanie konkretnego kanału conda z jego notacją podwójnego dwukropka i zdecydowanie zdefiniowanie konkretnej wersji oprogramowania, aby każda osoba, która uruchamia pipeline, otrzymała tę samą wersję.

Zauważ, że kontenery są nieco lepsze pod tym względem, ponieważ kiedy instalujesz coś za pomocą Condy, nadal będzie rozwiązywać wszystkie zależności dla tego pakietu, a one mogą się zmieniać w czasie. Nazywa się to dryfem zależności.

Więc kontenery jednak blokują cały stos całych zależności oprogramowania aż do samego dołu, więc możesz być nieco bardziej pewny, że A, to zadziała, i B, będzie odtwarzalne.

Więc jeśli jesteś w stanie używać Dockera lub Singularity lub Apptainera, zdecydowanie bym to polecił.

Teraz co jest fajne w tym, to to, że plik modułu, który jest napisany przez dewelopera pipeline'u, ma teraz zarówno Container, jak i Conda, i więc mówimy osobie, która uruchamia ten pipeline, nie obchodzi nas, jakiego rozwiązania pakowania oprogramowania używasz. Będzie działać zarówno z Dockerem, jak i z Condą, i tutaj można pobrać oprogramowanie w obu przypadkach.

Możemy otworzyć terminal i wypróbujmy to. Więc Nextflow run hello config --batch conda. I za pierwszym razem, gdy to działa z condą, będzie trochę wolno, gdy dojdzie do tego konkretnego procesu, ponieważ musi uruchomić "conda install".

I tworzy specjalne środowisko conda tylko dla tego jednego procesu. Więc nie używa mojego globalnego środowiska conda, które mam w moim terminalu. Tworzy jedno tylko dla tego jednego procesu. To jest dobre, ponieważ unika rzeczy takich jak konflikty zależności między różnymi procesami w Twoim workflow'ie. Jeśli Twoje procesy mają narzędzia, które potrzebują różnych wersji Pythona lub rzeczy tego typu, to jest w porządku, ponieważ używają różnych środowisk conda.

Nextflow buforuje te środowiska conda lokalnie, możesz zobaczyć, że mówi Ci, gdzie jest ta ścieżka, jest w katalogu work tutaj. I więc następnym razem, gdy uruchomię ten skrypt z Condą, będzie znacznie szybciej, ponieważ znajdzie to istniejące środowisko conda i po prostu go ponownie użyje. Ale za pierwszym razem, gdy to robimy, musi pójść i pobrać to, rozwiązać to, pobrać wszystkie zależności i wszystko skonfigurować.

Dobra, świetnie, uruchomił się. Możemy po prostu przypomnieć sobie, do czego pipeline jest obecnie skonfigurowany. Jeśli spojrzymy w plik konfiguracyjny, to było "custom-outdir-config" w tej chwili dla mnie. Zobacz, czy przejdę do tego katalogu bazowego. I zrobiłem --batch conda. Tam jest nasz podkatalog conda. Więc to zadziałało i tam jest nasze wyjście cowpy.

Więc pobrało cowpy, zainstalowało je w moim lokalnym systemie używając condy i uruchomiło proces. I co jest świetne, to to, że jako użytkownik końcowy, nie musiałem w ogóle myśleć o żadnym zarządzaniu oprogramowaniem. Nextflow po prostu to dla mnie załatwił. Powiedziałem, muszę użyć condy w tym systemie. Deweloper pipeline'u powiedział, jakich pakietów potrzebuję. A Nextflow zrobił resztę. Bardzo potężne.

Zauważ, że możesz faktycznie używać mieszanki różnych technologii. Więc mogę włączyć Dockera dla konkretnych procesów i condę dla innych procesów, lub powiedzieć, że niektóre procesy powinny po prostu używać jakiegokolwiek lokalnego oprogramowania, które miałem zainstalowane. To jest dość nietypowe, ale jest możliwe, i w niektórych przypadkach, na przykład, jeśli używasz pewnego oprogramowania, które może być trudne do spakowania w Dockerze, masz wyjście awaryjne.

## 4. Wybierz platformę wykonawczą

Więc to jest pakowanie oprogramowania. Druga część przenośności do innych systemów to miejsce, gdzie faktycznie działają zadania. W tej chwili działam zasadniczo na moim laptopie lub w tym Codespaces, który jest pojedynczym komputerem. Nie ma nic wymyślnego. Nextflow jest trochę sprytny w paralelizowaniu zadań tak dobrze, jak może, ale wszystko jest w jednym systemie.

Teraz, jeśli działasz na HPC, prawdopodobnie masz jakiś harmonogram zadań, taki jak SLURM lub PBS lub coś takiego, i przesyłasz zadania do tego harmonogramu i on rozdziela wszystkie zadania do różnych węzłów obliczeniowych.

Innym sposobem działania jest chmura. Więc może używasz AWS Batch, lub Azure Cloud, lub Google. I wszystkie one działają w podobnym systemie, gdzie masz harmonogram i przesyłasz zadania i są one przesyłane do różnych miejsc do obliczenia.

Teraz w odległej przeszłości, kiedy zaczynałem zajmować się bioinformatyką, oprogramowanie wszystkich do prowadzenia analiz było bardzo związane z ich infrastrukturą obliczeniową, co sprawiało, że było prawie niemożliwe do replikacji.

Ale z tym rozdzieleniem konfiguracji w Nextflow'ie i ze zdolnością Nextflow'a do interakcji z bardzo wieloma różnymi backendami infrastruktury obliczeniowej, bardzo łatwo jest wziąć nasz pipeline bez modyfikowania kodu pipeline'u w ogóle i po prostu to zamienić.

## 4.1. Kierowanie na inny backend

Więc jeśli przejdziemy do naszego pliku "nextflow.config" i możemy teraz umieścić trochę konfiguracji na poziomie procesu. Więc jeśli umieszczę na górze zakres process i mogę ustawić executor, i tutaj jest ustawiony na local, co jest wartością domyślną.

Zauważ, że ponieważ to jest poziom procesu, możemy kierować rzeczy do różnych procesów. I więc możesz faktycznie skonfigurować executory tak, aby były specyficzne dla procesu i mieć hybrydowe wykonanie, gdzie niektóre zadania mogą działać lokalnie, gdziekolwiek zadanie Nextflow'a jest wykonywane. Niektóre są przesyłane do różnych HPC, a niektóre mogą być przesyłane do chmury. Możesz być tak sprytny, jak chcesz.

Teraz bardzo trudno jest to zademonstrować w środowisku szkoleniowym takim jak to, ponieważ nie mam HPC, do którego mógłbym przesłać. Ale to, co mogę zrobić, to jeśli wpiszę slurm, możemy trochę oszukać i możesz poczuć to.

I jest to najbardziej interesujące dla osób, które są przyzwyczajone do działania na SLURM i wiedzą, jak wyglądają nagłówki SLURM. Ale jeśli zrobię Nextflow run, hello config. To się nie powiedzie, ponieważ będzie próbowało przesłać zadania do klastra, który nie istnieje. Więc otrzymamy jakiś błąd o tym, że sbatch nie jest dostępny.

Tak, napisane. To jest narzędzie. To jest narzędzie CLI, którego używasz do przesyłania zadań do klastra slurm. Ale to, co możemy zrobić, to możemy przejść i spojrzeć w nasz katalog work tutaj przez command click, otworzyć ten katalog i spojrzeć na .command.run. I możesz zobaczyć na górze pliku .command.run, mamy nasze nagłówki sbatch, mówiące teoretycznemu klastrowi SLURM, jak obsługiwać to przesłanie zadania.

I więc możesz zobaczyć, że Nextflow jest sprytny, robi wszystkie właściwe rzeczy. Po prostu nie mieliśmy klastra, do którego moglibyśmy przesłać.

## 5. Kontroluj alokacje zasobów obliczeniowych

Co jeszcze jest różne między różnymi infrastrukturami obliczeniowymi? Inną rzeczą jest to, ile dostępnych zasobów masz, i w rzeczywistości w wielu środowiskach obliczeniowych jest to wymóg, że musisz określić, ile procesorów i ile pamięci potrzebuje zadanie.

Ponownie, Nextflow to dla nas abstrahuje, tak że nie jest już specyficzne dla pojedynczego typu środowiska obliczeniowego, i możemy wpisać w zakresie poziomu procesu tutaj. CPUs equals one, memory equals two gigabytes. Nasz pipeline nie jest bardzo wymagający, więc to powinno być w porządku.

Teraz po prostu zgadłem te liczby tutaj, ale skąd wiesz, jaka jest rozsądna ilość zasobów do użycia? To dość trudne zadanie, aby przejść i przeszukać wszystkie te różne procesy dużego pipeline'u wielu próbek i zrozumieć, jakie było wykorzystanie zasobów.

Więc dobrym podejściem do tego jest ustawienie tych wartości na wysokie liczby na początek, po prostu żeby Twój pipeline działał bez żadnych błędów, a następnie poprosić Nextflow'a o wygenerowanie raportu użycia dla Ciebie.

To jest super łatwe do zrobienia, więc wrócę do terminala. Och, muszę pamiętać, aby ustawić to z powrotem na local, żeby mój pipeline faktycznie działał. I powiem Nextflow run i użyję flagi wiersza poleceń -with-report.

I mogę zostawić to puste i da domyślną nazwę pliku, ale dam mu konkretną nazwę pliku, żeby to zostało zapisane w konkretnym miejscu.

Naciśnij Enter i pipeline działa dokładnie tak jak normalnie, ale kiedy się skończy, wygeneruje dla mnie ładny raport HTML.

Więc na pasku bocznym tutaj mam ten plik HTML. Gdybym uruchamiał to lokalnie, po prostu bym go otworzył. Ponieważ jestem w Codespaces, kliknę prawym przyciskiem myszy na to i kliknę download, co pobierze to na mój lokalny komputer. I mogę po prostu łatwo otworzyć to w przeglądarce internetowej.

Nextflow może wygenerować taki raport dla każdego pipeline'u i ma naprawdę fajne informacje. Więc dobrą praktyką jest zawsze zapisywanie tych rzeczy. Mówi nam, kiedy uruchomiliśmy, gdzie uruchomiliśmy, czy to się powiodło, czy nie, jakie parametry zostały użyte, jakie było polecenie CLI, rzeczy tego typu.

I są również te wykresy dotyczące wykorzystania zasobów. Więc mówi nam, jaki procent wywołań CPU został użyty dla każdego procesu jako wykres pudełkowy tutaj, ponieważ jest wiele zadań dla każdego procesu, więc możemy zobaczyć rozkład.

Możesz zobaczyć nasze procesy tutaj, cowpy i collectGreetings miały tylko jedno zadanie, więc to jest tylko pojedyncza linia. I mamy zarówno CPU, jak i pamięć, i czas trwania zadania, i były bardzo szybkie.

Jeśli używasz Seqera Platform, nawiasem mówiąc, otrzymujesz te same wykresy wbudowane w interfejs Platform bez konieczności robienia czegokolwiek. Więc zawsze masz te informacje pod ręką.

Dobra, więc możemy użyć tego raportu i na prawdziwym uruchomieniu i poczuć, ile procesorów i ile pamięci jest używane przez nasz pipeline i wrócić i umieścić te wartości z powrotem w naszym pliku konfiguracyjnym, żeby następnym razem może nie żądać aż tak dużo. I możemy być trochę bardziej oszczędni.

Teraz możesz być naprawdę sprytny w konfigurowaniu plików konfiguracyjnych pipeline'u. I ponownie, jeśli używasz Seqera Platform, poszukaj małego przycisku, który wygląda jak żarówka. Ponieważ jeśli na niego klikniesz, wygeneruje wysoce zoptymalizowany plik konfiguracyjny, który jest dostosowany specjalnie do Twoich danych, Twojego uruchomienia i Twojego pipeline'u. Aby uruchomić go w najbardziej efektywny sposób.

Ale na razie powiem, że faktycznie domyślna liczba procesorów, którą Nextflow dawał, była w porządku, ale potrzebujemy tylko jednego gigabajta pamięci.

## 5.3. Ustaw alokacje zasobów dla konkretnego procesu

Teraz w prawdziwym życiu dość nietypowe jest, że wszystkie procesy w Twoim pipeline'ie będą potrzebowały tych samych wymagań. Możesz mieć coś takiego jak MultiQC jako narzędzie raportujące, które potrzebuje bardzo mało pod względem zasobów i działa dość szybko.

A następnie może masz coś, co indeksuje genom referencyjny lub robi jakieś dopasowanie lub robi jakąś inną pracę. Nie ma znaczenia, co to jest, co wymaga dużo zasobów. I więc dla tych różnych przesłań zadań do harmonogramu chcesz dać różne ilości zasobów.

Pod tym zakresem process możemy zdefiniować konfigurację, która kieruje konkretne procesy na różne sposoby.

Tutaj używamy withName, możemy również używać etykiet, i mogą one używać wzorca do kierowania jednego lub wielu procesów. Tutaj po prostu mówimy, że wszystkie procesy, które mają nazwę cowpy, ustawiają na dwa gigabajty pamięci i dwa procesory, i ponieważ jest to bardziej specyficzny selektor niż proces najwyższego poziomu, jest to nadpisywane w tych przypadkach, więc możesz zbudować ładny plik konfiguracyjny tutaj, który naprawdę dostosowuje wszystkie Twoje różne procesy w Twoim pipeline'ie, aby uczynić je naprawdę wydajnymi.

## 5.5. Dodaj limity zasobów

Teraz jako deweloper pipeline'u prawdopodobnie znam narzędzia całkiem dobrze i chcę, aby wszystko działało tak szybko i tak dobrze, jak to możliwe. Więc może być tak, że umieszczam dość wysokie liczby dla niektórych z nich, ponieważ wiem, że będzie działać znacznie szybciej, jeśli dam cowpy 20 procesorów zamiast dwóch.

To jest w porządku, dopóki nie przejdziesz do uruchomienia na swoim laptopie lub na GitHub Actions Continuous Integration test, lub jakimś innym systemie, który może nie mieć dostępnych 20 procesorów.

Teraz kiedy próbujesz uruchomić pipeline, zawiedzie, ponieważ Nextflow powie, nie mogę przesłać tego zadania nigdzie. Nie mam dostępnych zasobów.

Teraz, aby uniknąć tego twardego awarii, możemy dodać trochę więcej konfiguracji, która jest teraz specyficzna dla naszego systemu, o nazwie limity zasobów. I to wygląda tak. Jest pod zakresem process ponownie.

I limity zasobów, możesz określić zasadniczo sufit tego, co masz dostępne. To jest mapa tutaj, i możesz w tej mapie ustawić pamięć, procesory i czas.

Teraz to, co się dzieje, to kiedy Nextflow przesyła zadanie z procesu, patrzy na to, co jest żądane i zasadniczo po prostu robi minimum między tym a tym. Więc jeśli zażądaliśmy 20 procesorów, ale tylko cztery są dostępne, zażąda czterech. Pipeline nie zawiedzie i używa tak blisko tego, co zostało zaprojektowane przez dewelopera pipeline'u, jak to możliwe.

## 6. Użyj profili, aby przełączać się między predefiniowanymi konfiguracjami

Dobra. Powiedziałem, że limity zasobów tutaj mogą być specyficzne dla systemu, i może mam plik Nextflow config w moim pipeline'ie i wiem, że ludzie będą używać tego w różnych miejscach. Teraz zamiast zmuszać wszystkich do tworzenia własnego pliku Nextflow config za każdym razem, to, co mogę zrobić, to mogę zgrupować różne presety konfiguracji razem w profile konfiguracyjne.

Przewinę trochę w dół tutaj, a następnie tuż za params, ponieważ kolejność pliku konfiguracyjnego tutaj jest ważna, plik konfiguracyjny jest ładowany sekwencyjnie, więc umieszczę te profile po wszystkim innym, żeby nadpisywały wcześniej zdefiniowane parametry. I wkleję te profile z materiałów szkoleniowych.

Więc jest nowy zakres najwyższego poziomu o nazwie profiles. Możemy mieć dowolne nazwy tutaj. Więc mamy my_laptop i univ_hpc. I tutaj możemy zobaczyć, że ustawiamy te same parametry konfiguracyjne, które mieliśmy wcześniej. Teraz tylko w profilu. Więc mamy lokalny executor do działania na moim laptopie i przesyłam do klastra SLURM na HPC.

Używam Dockera lokalnie, condy na HPC, a system HPC ma znacznie wyższe limity zasobów.

Teraz mogę uruchomić pipeline z opcją CLI -profile, powiedzieć, którego profilu chcę użyć. Więc użyję my_laptop, a Nextflow zastosuje całą konfigurację w tym zakresie profilu. Więc mogę to teraz spróbować. To jest to samo polecenie co wcześniej. Nextflow run hello config i robię dash profile, pojedynczy myślnik, ponieważ jest to podstawowa opcja Nextflow'a, dash profile my_laptop.

Teraz będzie stosować wsadowo całą tę opcję konfiguracyjną. Och, i możesz zobaczyć, powiedziałem wcześniej, że to może się zdarzyć, że wymaganie procesu, zażądało czterech procesorów, a mam tylko dwa w tej instancji Codespaces.

Więc to jest dobra okazja, aby po prostu wypróbować limity zasobów procesu i powiedzieć, że mam tylko dwa procesory na moim laptopie lub w tym Codespaces. Teraz jeśli uruchomimy to ponownie, powinno ograniczyć to wymaganie do dwóch i miejmy nadzieję, że pipeline będzie działał. Świetnie.

## 6.2. Utwórz profil parametrów testowych

Zauważ, że te profile nie muszą mieć tylko konfiguracji dotyczącej ich infrastruktury. Możesz mieć grupowania dowolnej konfiguracji tutaj, w tym parametrów.

Więc inną rzeczą, którą bardzo często zobaczysz w pipeline'ach ludzi, jest profil testowy, który zawiera parametry, które normalnie przesłałbyś na podstawie użytkownika. Ale tutaj mamy zasadniczo różne rozsądne wartości domyślne, kiedy chcę uruchomić przypadki testowe.

I to jest świetne, ponieważ niekoniecznie muszę iść i określać wszystkie te rzeczy, które mogą być wymaganymi parametrami. W przeciwnym razie mogę po prostu powiedzieć dash profile test i po prostu będzie działać od razu.

Teraz coś do zauważenia jest to, że profile mogą być również łączone więcej niż jeden. Więc mogę zrobić profile my_laptop tutaj, a następnie również dodać test. Nie robię profile dwa razy. Po prostu robię listę oddzieloną przecinkami tutaj bez spacji. I będzie stosować te profile w kolejności. Więc weźmie konfigurację z profilu my_laptop, a następnie zastosuje konfigurację test na wierzchu.

Naprawdę wygodne i możesz zobaczyć, jak możesz skonfigurować wiele rozsądnych domyślnych grup tutaj, aby ułatwić uruchomienie Twojego pipeline'u.

## 6.3. Użyj nextflow config, aby zobaczyć rozwiązaną konfigurację

Mam nadzieję, że przekonałem Cię, że rozwiązywanie konfiguracji Nextflow'a jest potężne, ale nie obwiniałbym Cię, gdybyś trochę krzyżował oczy w tym momencie po tym, jak powiedziałem około 20 różnych sposobów dostarczania konfiguracji i podania wszystkich tych różnych warstw jak skórka cebuli.

Więc jeśli kiedykolwiek czujesz się niepewny co do tego, jaka jest ostateczna rozwiązana konfiguracja dla Nextflow'a, wiedz, że jest polecenie o nazwie "nextflow config" i możemy to uruchomić i powie nam, jaka jest rozwiązana konfiguracja w naszej obecnej lokalizacji.

Więc kiedy uruchamiam to tutaj, znajduje plik "nextflow.config" w bieżącym katalogu roboczym i przetwarza całą różną konfigurację i daje mi rozwiązane wyjście.

Zauważ, że plik konfiguracyjny Nextflow'a może również przyjąć opcję CLI profilu. Więc jeśli powiem mu, aby rozwiązał w profilach my_laptop i test, i możesz zobaczyć, że również zastosował limity zasobów tutaj z opcji konfiguracyjnej my_laptop, a także ustawił parametry, które były w teście.

Więc to jest miły sposób, aby po prostu zbadać, jak działa rozwiązywanie konfiguracji, jeśli jesteś w ogóle niepewny.

## Podsumowanie

Dobra, to wszystko. To jest konfiguracja Nextflow'a w pigułce. Możesz zrobić wiele rzeczy z konfiguracją. Jest naprawdę potężna. Ale to są większość typowych przypadków użycia, które będziesz robić, i te koncepcje mają zastosowanie do wszystkich różnych opcji.

Poklepcie się po plecach, ponieważ to jest koniec kursu szkoleniowego Hello Nextflow. Mam nadzieję, że jesteś teraz pewny zarówno w pisaniu własnego pipeline'u Nextflow'a od podstaw, konfigurowaniu go i uruchamianiu go, i znasz wszystkie zawiłości i rzeczy, na które należy uważać.

Jest jeszcze jeden quiz, który możesz wypróbować na stronie szkoleniowej konfiguracji. Więc zejdź w dół i wypróbuj to i upewnij się, że zrozumiałeś wszystkie te części dotyczące konfiguracji.

I dołącz do nas w ostatnim wideo tylko na szybkie podsumowanie o niektórych następnych krokach, które mogą być dobre do zrobienia po tym kursie szkoleniowym.

Dziękuję za wytrwanie z nami. Dobra robota i do zobaczenia w następnym wideo.
