# Orientacja - Transkrypcja wideo

<div class="video-wrapper">
  <iframe width="560" height="315" src="https://www.youtube.com/embed/PIjOdFaYwWA?si=y8lAedhEHWaTV4zd&amp;list=PLPZ8WHdZGxmWKozQuzr27jyMGqp9kElVK&amp;cc_load_policy=1" title="YouTube video player" frameborder="0" allow="accelerometer; autoplay; clipboard-write; encrypted-media; gyroscope; picture-in-picture; web-share" referrerpolicy="strict-origin-when-cross-origin" allowfullscreen></iframe>
</div>

!!!note "Ważna uwaga"

    Ta strona zawiera tylko transkrypcję. Pełne instrukcje krok po kroku znajdziesz w [materiałach kursu](../00_orientation.md).

## Powitanie

Cześć i witaj w Hello Nextflow. Nazywam się Phil Ewels. Jestem kierownikiem ds. produktów oprogramowania open source w Seqera, firmie stojącej za Nextflow'em.

Ten kurs to praktyczne wprowadzenie do tworzenia workflow'ów przy pomocy Nextflow'a. Został zaprojektowany dla osób, które są zupełnie nowe w Nextflow'ie i chcą tworzyć własne pipeline'y.

Wszystkie przykłady dotyczą prostego przetwarzania tekstu, dzięki czemu możesz skupić się na koncepcjach Nextflow'a bez potrzeby posiadania specjalistycznej wiedzy dziedzinowej — wystarczy pewna znajomość linii poleceń.

Przejdziemy przez podstawy Nextflow'a: pisanie procesów, łączenie ich w wieloetapowe workflow'y, zarządzanie zależnościami programowymi za pomocą kontenerów oraz konfigurowanie pipeline'ów dla różnych środowisk obliczeniowych. Pod koniec zbudujesz działający pipeline od podstaw.

Ten kurs skupia się na _tworzeniu_ pipeline'ów. Jeśli chcesz tylko _uruchamiać_ istniejące pipeline'y bez zanurzania się zbytnio w kod, mamy krótszy kurs "Nextflow Run", który może lepiej Ci odpowiadać.

Gdy już opanujesz podstawy, mamy również kursy kontynuacyjne, które stosują te koncepcje do prawdziwych analiz naukowych. Nauczymy Cię, jak korzystać z pipeline'ów społeczności nf-core i najlepszych praktyk.

Jeśli utkniesz, przejdź do community.seqera.io. Tam znajduje się aktywne forum społeczności z sekcją dedykowaną wyłącznie pytaniom szkoleniowym. Możesz z niego korzystać w dowolnym momencie, jednak co kwartał organizujemy również tygodnie szkoleniowe z osobami specjalnie gotowymi do pomocy. Jeśli więc uczestniczysz w szkoleniu podczas jednego z nich, zdecydowanie nie wstydź się i poproś o pomoc.

Możesz również spróbować poprosić o pomoc Seqera AI. Świetnie radzi sobie z wyjaśnianiem kodu Nextflow'a i pomocą w debugowaniu.

Gdy będziesz gotowy do uruchamiania Nextflow'a na dużą skalę, Seqera Platform to najlepsze miejsce do tego. Działa na Twojej infrastrukturze bez żadnego uzależnienia od dostawcy, oferując wszystko od uruchamiania pipeline'ów przez monitorowanie w czasie rzeczywistym po interaktywne środowiska analityczne. Ale na razie skupmy się po prostu na podstawach.

No dobra, zaczynajmy.

## training.nextflow.io

Okej. Pierwszą rzeczą do zauważenia jest to, że wszystkie kursy szkoleniowe na training.nextflow.io są bardzo interaktywne. Chodzi o to, żebyś podążał za materiałem szkoleniowym i moimi instrukcjami, a my wspólnie przechodzimy przez materiał. Będziesz potrzebować dwóch rzeczy: swojego laptopa i tej strony internetowej. I to w zasadzie tyle.

To jest strona główna, tak jak wygląda dzisiaj, gdy nagrywam. Możesz zobaczyć przegląd różnych rzeczy, tło i różne kursy, które mamy, a lista ciągle się rozrasta.

Nextflow dla początkujących to miejsce, w którym jesteśmy. Są tu dwa kursy: Nextflow Run, który jest innym kursem, oraz Hello Nextflow, którym się zajmujemy.

Możesz też zobaczyć wszystkie różne kursy na pasku bocznym. Mogę przeskoczyć do Hello Nextflow i zobaczyć wszystkie różne rozdziały, które przejdziemy razem.

Jest kilka innych ważnych rzeczy do zauważenia. Po pierwsze, materiał szkoleniowy jest wersjonowany, więc możesz zobaczyć tutaj u góry. Mówi 3.0 latest, która w momencie nagrywania jest najnowszą stabilną wersją. To się zmieni z czasem. Publikujemy nowe kursy i aktualizujemy materiał na bieżąco. Więc jeśli jest 3.1 lub 3.2, nie martw się zbytnio. Jeśli jest 4.0, to prawdopodobnie jest nowe wideo i powinieneś je może znaleźć, ponieważ prawdopodobnie będą znaczące aktualizacje.

Kolejna rozwijana lista u góry to ta z językiem. To jest zupełnie nowe dla wersji 3.0. Wzięliśmy wcześniej przetłumaczony materiał, który został wykonany przez ludzi ręcznie, i przepuściliśmy go przez LLM i stworzyliśmy całą tę nową infrastrukturę do utrzymywania różnych tłumaczeń materiału szkoleniowego za pomocą tłumaczenia LLM.

Więc teraz mamy wszystkie te fantastyczne tłumaczenia tutaj. Więc jeśli chcesz słuchać po koreańsku, możesz załadować całą stronę po koreańsku i podążać tam. To samo dla wszystkich innych języków, hindi i niemieckiego i tak dalej. Ja będę podążał po angielsku. To jest główny język, w którym piszemy materiał.

Jest kilka innych przycisków, jeśli lubisz tryb jasny. Zamiast ciemnego trybu możesz śledzić stronę w trybie jasnym tutaj u góry.

A także wszystko, na co patrzymy, znajduje się w jednym repozytorium GitHub, które jest open source, nazywa się nextflow-io/training. I jeśli klikniesz ten przycisk w dowolnym momencie, przejdzie do repozytorium GitHub. Wrócimy do tego za chwilę.

## Konfiguracja GitHub Codespaces

Okej, więc teraz masz to otwarte w zakładce przeglądarki. Przejdźmy do Hello Nextflow i kliknijmy. Możesz zobaczyć na stronie wprowadzającej, że informuje nas o kilku wymaganiach, przeglądzie i planie lekcji tego, co w przybliżeniu omówimy, a następnie zagłębimy się w rozpoczęcie pracy.

Są różne sposoby, w które możesz wykonać ten interaktywny samouczek. Jeśli chcesz, możesz to zrobić lokalnie na swoim własnym komputerze z własną instalacją Nextflow'a. Jeśli klikniemy do Opcji środowiska, możesz zobaczyć więcej szczegółów, jak to zrobić używając lokalnych Devcontainerów lub możesz też po prostu zainstalować całe oprogramowanie lokalnie, z ręczną instalacją.

Pracujemy nad tym, aby to dobrze działało z Seqera Studios, więc to kolejna opcja. Ale najbardziej powszechną opcją teraz jest użycie GitHub Codespaces.

Codespaces konfiguruje środowisko piaskownicy na zdalnym serwerze prowadzonym przez GitHub. Jest to bezpłatne przez określoną ilość użytkowania, co zazwyczaj wystarcza na szkolenie. I skonfiguruje Ci instancję VS Code, IDE, gdzie możesz uzyskać dostęp do wszystkich plików z repozytorium, uruchamiać Nextflow'a i wszystko. I wstępnie skonfigurowaliśmy dla Ciebie Codespaces. Więc ma wszystko, czego potrzebujesz.

Piękno tego jest takie, że to tylko jedno kliknięcie, aby skonfigurować Codespace. Jest taki sam dla wszystkich i wiemy, że masz już zainstalowane wszystkie wymagania wstępne, więc jest miło i szybko.

Więc pierwszą rzeczą do zrobienia jest przejście do "Getting Started". Poszukaj tego przycisku, który mówi _Open in Codespaces_. Nacisnę cmd + klik, aby otworzyć to w nowej zakładce i przeniesie nas do GitHub.

Tak to wygląda. Możemy zobaczyć, że ustawiliśmy już wszystkie opcje dla Ciebie. Jeśli chcesz, możesz kliknąć opcje zmiany. Możesz tu zrobić kilka rzeczy. Możesz na przykład wybrać większą maszynę instancji, jeśli okaże się, że crashuje, bo kończy się pamięć lub coś takiego. Lub ustawić konkretne wersje materiału szkoleniowego. Ale zazwyczaj możesz po prostu użyć tego, co dla Ciebie skonfigurowaliśmy i możesz to zobaczyć. W tym przypadku używa wydania 3.0.

Więc kliknę create new Codespace. I to mnie przenosi.

Zauważ również, że mówi no Codespace to resume tam. Jeśli wcześniej utworzyłem Codespace, ponowne kliknięcie tego przycisku na materiale szkoleniowym przeniesie mnie na tę samą stronę i wyświetli wszystkie Codespaces, które już mam uruchomione. Wtedy możesz po prostu wskoczyć z powrotem do nich i kontynuować tam, gdzie skończyłeś. Więc nie ma znaczenia, jeśli zamknąłeś laptopa.

Automatycznie wyłączają się po kilku minutach braku aktywności, ale nie ma problemu. Możesz je po prostu zrestartować.

Gdy uruchomisz nowy Codespace, będzie się ładować na tej stronie przez dłuższy czas. Więc teraz jest dobry moment na krótką przerwę. Może zapomniałeś pójść do toalety lub chcesz filiżankę herbaty przed rozpoczęciem? Idź teraz, podczas czekania na to, bo będzie się kręcić przez chwilę.

Szybko, podczas gdy na to czekamy, przejdę również do github.com/codespaces i po prostu pokażę, że to jest strona przeglądu, gdzie możesz zobaczyć wszystkie różne Codespaces, które obecnie masz uruchomione.

Możesz zobaczyć, że mam tutaj jeden dla nextflow-io/training. Brak zmian, bo jeszcze nic w nim nie zrobiłem. Ilość zasobów, których używa, i możesz zobaczyć, że w tej chwili się konfiguruje. Mogę tutaj kliknąć tę małą rozwijalną listę i kliknąć usuń. Więc jeśli przypadkowo skonfigurujesz wiele Codespaces i nie używasz niektórych, możesz usunąć stare i posprzątać.

Na koniec, jeszcze jeden sposób, aby się tu dostać. Jeśli przejdziemy do repozytorium GitHub. I to działa dla każdego repozytorium GitHub. Kliknij code. Masz polecenia do klonowania repozytorium lokalnie. I jest zakładka o nazwie Codespaces. I znowu, możesz utworzyć nowy, i możesz zobaczyć te, które już działają.

Więc znowu, jeśli zapomnisz, jak utworzyłeś swój Codespace, zawsze możesz w ten sposób do niego wrócić.

## Interfejs VS Code

Okej, kompilator skończył i teraz zaczyna się ładowanie GitHub Codespaces. Nie zawsze trwa to tak długo, więc nie martw się. To tylko za pierwszym razem, gdy tworzysz Codespace. Jeśli wrócisz do tego, który już istnieje, jest znacznie szybciej.

Nie bądź zbyt niecierpliwy, jeśli to pierwszy raz, jeszcze się nie skończył, mimo że zaczyna nam dawać interfejs.

Ale podczas gdy czekamy na ostatnie rzeczy do skonfigurowania, po prostu pokażę Ci interfejs na wypadek, gdybyś był trochę nieznajomy z VS Code.

Po pierwsze, jest pasek boczny czatu dla rzeczy AI, który nie jest nam potrzebny. Więc zamknę to, pozbędę się tego i uwolnię trochę miejsca.

Po lewej stronie mamy eksplorator plików, który pokazuje nam wszystkie pliki w repozytorium Git, czyli obszarze roboczym, który utworzyliśmy. Zauważ, to nie są pliki lokalne. To wszystko jest na zdalnym serwerze, gdzie pracujemy. Możesz przeciągać i upuszczać lokalne pliki i rzeczy, ale w większości nie będziemy o tym myśleć dzisiaj. Pracujemy po prostu czysto zdalnie.

Są inne narzędzia w tym pasku bocznym, na przykład wyszukiwanie. Więc możesz przeszukać wszystkie pliki w repozytorium na raz. I gdybyśmy wykonywali pracę rozwojową nad repozytorium szkoleniowym, moglibyśmy wykonać integrację z kontrolą źródła z Git i debugowaniem i innymi rzeczami.

Inne rzeczy to, jest główne okno edycji kodu tutaj u góry, które właśnie załadowało podgląd readme, który jest dla materiału szkoleniowego. Więc w tym przypadku wyświetla markdown, ale normalnie będzie to edytor kodu.

A poniżej mamy terminal, w którym będziemy uruchamiać wszystkie nasze polecenia i wchodzić w interakcje bezpośrednio z Nextflow'em.

Wszystko w Codespace jest wstępnie zainstalowane, więc polecenie Nextflow'a jest już tam i tak dalej.

Okej. Gdy dojdziesz tak daleko, powinno być prawie gotowe. Możesz teraz zobaczyć, że pobrał serwer językowy Nextflow'a i skonfigurował dla nas kilka rozszerzeń w VS Code, w tym rozszerzenie Nextflow, które będzie użyteczne. Więc mogę to zamknąć i mogę zamknąć README.md.

A teraz możesz zobaczyć, że mam trochę więcej po lewej stronie. Jestem tu trochę przybliżony, ale jeśli pomniejszę, możesz zobaczyć, że jeden z przycisków mówi Nextflow z ikoną Nextflow'a i ma tam kilka fajnych rzeczy do eksploracji projektu i tak dalej, do których wrócimy później.

Okej. Na wypadek, gdybyś kiedykolwiek zgubił którykolwiek z tych paneli, te przyciski w prawym górnym rogu są naprawdę użyteczne i po prostu pokazują i ukrywają rzeczy. Więc to pokazuje i ukrywa Eksplorator, pokazuje i ukrywa terminal na dole. I tak dalej.

Będę ich używał dość często, ponieważ jestem bardzo przybliżony, więc próbuję pomóc Ci zobaczyć cały tekst na moim ekranie, więc przydatne jest móc przejść na pełny ekran z terminalem, a następnie ukryć go, gdy patrzymy na kod. Ale przez większość czasu możesz po prostu mieć wszystkie te rzeczy otwarte jednocześnie.

Okej, co jeszcze zobaczyć? Niewiele więcej. Zauważ, że Nextflow, jak mówię, jest zainstalowany. Więc mogę wpisać "nextflow -version" i powinno się pojawić mówiąc, którą wersję mamy zainstalowaną.

Jest również trochę innych rzeczy zainstalowanych tutaj. Na końcu każdego rozdziału mamy zestaw pytań quizowych, na przykład na stronie internetowej. I możesz je również wykonać w terminalu, jeśli chcesz, wpisując quiz.

Jest kilka innych skrótów klawiaturowych, których będę używał, na wypadek, gdybyś był ciekaw. Na przykład właśnie teraz nacisnąłem cmd+K na moim Macu i to wyczyściło terminal, aby pozbyć się całego poprzedniego wyjścia. Więc jest to miłe, aby zachować porządek. Jeśli zobaczysz, że to robię, tak to robię.

A także, jeśli jesteś nowy w terminalu, pamiętaj, że możesz używać tabulatora do automatycznego uzupełniania, którego będę używał dużo do automatycznego uzupełniania ścieżek.

Więc mogę zobaczyć po lewej stronie, że jest folder o nazwie Hello Nextflow, przez który będziemy pracować. Jeśli zrobię "ls", aby wyświetlić pliki, mogę zrobić "hel", nacisnąć tab, automatycznie uzupełnia. I tak to jest bardzo szybki sposób na uzupełnianie ścieżek.

## Otwieranie tylko folderu Hello Nextflow

Okej. To świetnie. Jest tutaj jednak dużo rzeczy w tym repozytorium.

Są wszystkie pliki do generowania strony internetowej i jest wiele różnych kursów tutaj, i możesz to zrobić z tego miejsca i po prostu kliknąć w folder "Hello Nextflow". Ale fajnie jest faktycznie po prostu skupić się wyłącznie na tym.

Możesz ustawić to jako swój obszar roboczy z grupą klikania tutaj i ustawiania katalogu projektu i tak dalej. Ale najłatwiejszym sposobem jest wpisać code, które jest poleceniem CLI do uruchamiania VS Code, a następnie "hello-nextflow".

To otworzy nową zakładkę przeglądarki i możesz zamknąć starą. I wygląda dokładnie tak samo. Ale teraz możesz zobaczyć, że jesteśmy w tym podkatalogu i wszystkie inne pliki są niewidoczne, i mamy czystszą konfigurację.

Możesz zobaczyć tutaj, że również bieżący katalog roboczy jest teraz w folderze Hello Nextflow. Więc miło i czysto. Nie musimy się martwić o bycie w złym miejscu. Okej.

## Nowa składnia Nextflow na 2026

Jest jedna specjalna rzecz, o której muszę wspomnieć w tym momencie. Teraz, na początku 2026 roku, zaczynamy wprowadzać różne funkcje do Nextflow'a, a jedną z dużych nowych jest nowy parser składni języka wewnątrz Nextflow'a.

Zasadniczo silnik, który odczytuje Twoje pliki Nextflow i rozumie to, dla środowiska uruchomieniowego. Są pewne zmiany w składni i naprawdę ważne jest, abyś używał Nextflow'a z poprawnie włączonym parserem składni.

Potrzebujesz do tego dwóch rzeczy. Potrzebujesz aktualnej wersji Nextflow'a i musisz upewnić się, że jest włączony.

Jeśli znowu zrobię "nextflow -version", zobaczysz, że Codespaces działa z 25.10.2, a 25.10 to minimalna wersja, aby móc używać tego.

Jeśli używasz 26.04, która dla mnie jeszcze nie wyszła, ale wkrótce wyjdzie, to będzie domyślnie uruchamiać nowy parser składni i nie musisz robić nic więcej.

Ale jeśli uruchamiasz 25.10, musisz włączyć ścisły parser składni, jak się go nazywa, lub parser składni v2.

Odbywa się to za pomocą zmiennej środowiskowej. Jest już ustawiona w Codespaces, więc nie musisz nic robić. Ale jeśli uruchamiasz lokalnie, musisz to ustawić, i mogę to zweryfikować robiąc "echo $NXF_SYNTAX_PARSER" i powinno być ustawione na v2.

Więc jeśli uruchamiasz lokalnie, po prostu zrób "export NXF_SYNTAX_PARSER=v2". Proste. Ale pamiętaj, aby to zrobić, bo w przeciwnym razie zobaczysz jakieś dziwne rozbieżności i błędy w trakcie.

Jeśli jesteś w ogóle niepewny co do którejkolwiek z tych rzeczy wokół wersji Nextflow'a i parsera składni, po pierwsze, pamiętaj, nie musisz się martwić, jeśli jesteś w Codespaces. Wszystko powinno być prawidłowo skonfigurowane. Ale po drugie, jeśli przejdziesz do materiału szkoleniowego Nextflow, jeśli zejdziesz w dół, porozmawiamy o wymaganiach wersji, jest link tutaj, który zabiera Cię na dół do strony pomocy dotyczącej eksplorowania wersji, i to przechodzi przez wszystko szczegółowo.

Warto to przeczytać, jeśli masz chwilę, bo pomaga wyjaśnić, jakie są niektóre z różnych terminów, które możesz usłyszeć, gdy zaczniesz używać Nextflow'a. Rzeczy takie jak DSL1, DSL2, parser składni jeden, parser składni dwa i tak dalej. Więc warto po prostu sprawdzić to, a to powtarza część tego, co właśnie powiedziałem.

Jest to również naprawdę użyteczne, jeśli wcześniej napisałeś kod Nextflow'a i wracasz na odświeżenie. Mówi Ci o niektórych rzeczach, które się zmieniły i linkuje Cię do części dokumentacji Nextflow, która mówi Ci, jak zaktualizować swój kod Nextflow'a.

## Pliki kursu

Okej. Ostatnią rzeczą do zaznajomienia się jest po prostu zobaczenie plików, które są w tym katalogu. Możesz albo spojrzeć w pasek boczny albo często w materiałach szkoleniowych używamy polecenia tree, -L, które jest liczbą poziomów do zajrzenia. Powiemy dwa, i jeśli zrobię to na pełnym ekranie, zobaczysz, że to zasadniczo odzwierciedla dokładnie to, co widzimy na pasku bocznym, ale wyklucza ukryte pliki, które zaczynają się od kropki.

Więc pliki \*.nf, to oznacza Nextflow. Więc to są pliki skryptowe Nextflow'a i jest tutaj plik startowy dla każdego z różnych rozdziałów materiału szkoleniowego, który otworzymy i zbadamy, a następnie edytujemy.

Będziemy zmieniać te pliki w trakcie, więc pod koniec każdego rozdziału pliki powinny wyglądać prawie tak samo jak na początku rozdziału dla następnego. Ale dajemy Ci te różne pliki, więc zawsze możesz zacząć od nowa i nie martwić się zbytnio zepsuciem składni.

Jeśli musisz porównać z czymś, co zdecydowanie powinno działać, możesz sprawdzić w folderze solutions, i to jest jak stan końcowy dla każdego z rozdziałów, więc możesz porównać to, co napisałeś, z tym, co tam jest.

Jest katalog data. Ma tylko plik greetings.csv, którego użyjemy jako przykładowych danych wejściowych w części kursu, i rzeczy takie jak plik konfiguracyjny i niektóre parametry, które opiszemy później w kursie.

## Podsumowanie

Okej, więc teraz mam nadzieję, że wszystko działa. Twój ekran wygląda tak samo jak mój i rozumiesz, jak dostać się do wszystkiego i czym są wszystkie różne pliki.

Jeśli przewiniesz w dół strony na getting started, małe pole wyboru powinno powiedzieć, że rozumiem, co robię. Moje środowisko jest uruchomione i działające, i Twój katalog roboczy jest prawidłowo ustawiony na folder "Hello Nextflow".

Jeśli zaznaczyłeś wszystkie te i wyglądają na zielone, możemy przejść do następnego wideo i następnego rozdziału, którym jest część pierwsza. Hello World. Do zobaczenia za chwilę.
