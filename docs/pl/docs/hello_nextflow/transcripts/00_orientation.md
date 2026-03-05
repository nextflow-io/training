# Orientacja - Transkrypcja wideo

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Tłumaczenie wspomagane przez AI - [dowiedz się więcej i zasugeruj ulepszenia](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

<div class="video-wrapper">
  <iframe width="560" height="315" src="https://www.youtube.com/embed/PIjOdFaYwWA?si=y8lAedhEHWaTV4zd&amp;list=PLPZ8WHdZGxmWKozQuzr27jyMGqp9kElVK&amp;cc_load_policy=1" title="YouTube video player" frameborder="0" allow="accelerometer; autoplay; clipboard-write; encrypted-media; gyroscope; picture-in-picture; web-share" referrerpolicy="strict-origin-when-cross-origin" allowfullscreen></iframe>
</div>

!!!note "Ważna uwaga"

    Ta strona zawiera tylko transkrypcję. Pełne instrukcje krok po kroku znajdziesz w [materiale szkoleniowym](../00_orientation.md).

## Powitanie

Cześć i witaj w Hello Nextflow. Nazywam się Phil Ewels. Jestem Product Managerem ds. Oprogramowania Open Source w Seqera, firmie stojącej za Nextflow'em.

Ten kurs to praktyczne wprowadzenie do budowania workflow'ów w Nextflow'ie. Jest przeznaczony dla osób, które są całkowicie nowe w Nextflow'ie i chcą tworzyć własne pipeline'y.

Wszystkie przykłady to proste przetwarzanie tekstu, dzięki czemu możesz skupić się na koncepcjach Nextflow'a bez potrzeby posiadania wiedzy specjalistycznej – wystarczy znajomość wiersza poleceń.

Przejdziemy przez podstawy Nextflow'a: pisanie procesów, łączenie ich w wieloetapowe workflow'y, zarządzanie zależnościami oprogramowania za pomocą kontenerów oraz konfigurowanie pipeline'ów dla różnych środowisk obliczeniowych. Pod koniec zbudujesz działający pipeline od podstaw.

Ten kurs koncentruje się na _tworzeniu_ pipeline'ów. Jeśli chcesz tylko _uruchamiać_ istniejące pipeline'y bez zagłębiania się zbytnio w kod, mamy krótszy kurs "Nextflow Run", który może Ci bardziej odpowiadać.

Gdy opanujesz tutaj podstawy, mamy również kursy kontynuacyjne, które stosują te koncepcje do rzeczywistej analizy naukowej. Nauczymy Cię, jak korzystać z pipeline'ów i najlepszych praktyk społeczności nf-core.

Jeśli utkniesz, odwiedź community.seqera.io. Jest tam aktywne forum społeczności z sekcją dedykowaną wyłącznie pytaniom szkoleniowym. Możesz z niego korzystać w dowolnym momencie, jednak organizujemy również kwartalne tygodnie szkoleniowe z osobami gotowymi specjalnie pomóc. Więc jeśli robisz szkolenie podczas jednego z nich, zdecydowanie nie krępuj się i proś o pomoc.

Możesz również spróbować poprosić o pomoc Seqera AI. Świetnie radzi sobie z wyjaśnianiem kodu Nextflow'a i pomaganiem w debugowaniu.

Gdy będziesz gotowy do uruchamiania Nextflow'a na dużą skalę, Seqera Platform to najlepsze miejsce do tego. Działa na Twojej infrastrukturze bez żadnego uzależnienia od dostawcy, ze wszystkim – od uruchamiania pipeline'ów po monitorowanie w czasie rzeczywistym, po interaktywne środowiska analityczne. Ale na razie skupmy się po prostu na podstawach.

No dobrze, zaczynajmy.

## training.nextflow.io

Okej. Pierwszą rzeczą do zauważenia jest to, że wszystkie kursy szkoleniowe na training.nextflow.io są bardzo interaktywne. Chodzi o to, żebyś podążał za materiałem szkoleniowym i moimi instrukcjami, a my przechodzimy przez materiał razem. Więc będziesz potrzebować dwóch rzeczy: Twojego laptopa i tej strony internetowej otwartej. I to w zasadzie wszystko.

To jest strona główna, jak wygląda dzisiaj, gdy to nagrywam. Widzisz, że jest przegląd różnych rzeczy, tło i różne kursy, które mamy, a lista cały czas rośnie.

Nextflow for newcomers to miejsce, w którym jesteśmy. Są tu dwa kursy: Nextflow Run, który jest innym kursem, i Hello Nextflow, który nas interesuje.

Możesz również zobaczyć wszystkie różne kursy na pasku bocznym. Mogę przejść do Hello Nextflow i zobaczyć wszystkie różne rozdziały, przez które będziemy pracować razem.

Jest kilka innych ważnych rzeczy do zauważenia. Po pierwsze, materiał szkoleniowy jest wersjonowany, więc widzisz tutaj na górze. Mówi 3.0 latest, co w momencie nagrywania jest najnowszą stabilną wersją. To będzie się zmieniać z czasem. Wypuszczamy nowe kursy i aktualizujemy materiał. Więc jeśli jest 3.1 lub 3.2, nie martw się zbytnio. Jeśli jest 4.0, to prawdopodobnie jest nowe wideo i powinieneś może je znaleźć, bo prawdopodobnie będą znaczące aktualizacje.

Kolejna lista rozwijana na górze to ta z językiem. To jest zupełnie nowe dla wersji 3.0. Wzięliśmy wcześniej przetłumaczony materiał, który był robiony przez ludzi, ręcznie, i przepuściliśmy go przez LLM i skonfigurowaliśmy całą tę nową infrastrukturę do utrzymywania różnych tłumaczeń materiału szkoleniowego przy użyciu tłumaczenia LLM.

Więc teraz mamy wszystkie te fantastyczne tłumaczenia tutaj. Więc jeśli chcesz słuchać po koreańsku, możesz załadować całą stronę po koreańsku i podążać tam. To samo dla wszystkich tych innych języków, hindi i niemieckiego i tak dalej. Będę podążał po angielsku. To jest główny język, w którym piszemy materiał.

Kilka innych przycisków – jeśli lubisz tryb jasny zamiast trybu ciemnego, możesz śledzić stronę w trybie jasnym na górze tutaj.

A potem wszystko, na co patrzymy, jest w jednym repozytorium GitHub, które jest open source, nazywa się nextflow-io/training. I jeśli klikniesz ten przycisk w dowolnym momencie, przejdzie do repozytorium GitHub. Wrócimy do tego za chwilę.

## Konfiguracja GitHub Codespaces

Okej, więc teraz masz to otwarte w karcie przeglądarki. Przejdźmy do Hello Nextflow i kliknijmy. Widzisz na stronie wprowadzającej, że mówi nam o kilku wymaganiach, przeglądzie i planie lekcji tego, co w przybliżeniu będziemy omawiać, a potem zanurzymy się w rozpoczęcie pracy.

Są różne sposoby, w jakie możesz zrobić ten interaktywny samouczek. Jeśli chcesz, możesz to zrobić lokalnie na swoim własnym komputerze z własną instalacją Nextflow'a. Jeśli klikniemy na Opcje środowiska, zobaczysz więcej szczegółów, jak to zrobić, używając lokalnych Devcontainerów lub możesz też po prostu zainstalować całe oprogramowanie lokalnie, z ręczną instalacją.

Pracujemy nad tym, żeby to działało ładnie z Seqera Studios, więc to kolejna opcja. Ale najbardziej powszechną teraz jest użycie GitHub Codespaces.

Codespaces konfiguruje środowisko piaskownicy na zdalnym serwerze prowadzonym przez GitHub. I jest darmowe dla pewnej ilości użycia, co zwykle wystarcza na szkolenie. I skonfiguruje Ci instancję VS Code, IDE, gdzie możesz uzyskać dostęp do wszystkich plików z repozytorium, uruchomić Nextflow'a i wszystko. I wstępnie skonfigurowaliśmy Codespaces dla Ciebie. Więc ma wszystko, czego potrzebujesz.

Piękno tego polega na tym, że to tylko jedno kliknięcie, aby skonfigurować Codespace. Jest takie samo dla wszystkich i wiemy, że masz już zainstalowane wszystkie wymagania wstępne, więc jest miło i szybko.

Więc pierwszą rzeczą do zrobienia jest przejście do "Pierwsze kroki". Poszukaj tego przycisku, który mówi _Otwórz w Codespaces_. Zrobię cmd + kliknięcie, aby otworzyć to w nowej karcie, i zabiera nas do GitHub.

Tak to wygląda. Widzimy, że ustawiliśmy wszystkie opcje tutaj dla Ciebie. Jeśli chcesz, możesz kliknąć zmień opcje. Kilka rzeczy, które możesz tutaj zrobić. Możesz dać większą maszynę instancji, na przykład, jeśli okaże się, że się zawiesza, bo kończy się pamięć lub coś takiego. Lub ustawić konkretne wersje materiału szkoleniowego. Ale zwykle możesz po prostu iść z tym, co skonfigurowaliśmy tutaj i widzisz to. W tym przypadku używa wydania 3.0.

Więc kliknę utwórz nowy Codespace. I to mnie zabiera.

Zauważ również, że mówi brak Codespace do wznowienia tam. Jeśli wcześniej utworzyłem Codespace, kliknięcie tego przycisku ponownie na materiale szkoleniowym zabierze mnie na tę samą stronę i wyświetli wszystkie Codespaces, które już mam uruchomione. Wtedy możesz po prostu wskoczyć z powrotem do nich i kontynuować tam, gdzie skończyłeś. Więc nie ma znaczenia, czy zamknąłeś laptopa.

Automatycznie się wyłączają po kilku minutach bezczynności, ale nie ma problemu. Możesz je po prostu zrestartować.

Gdy uruchomisz nowy Codespace, będzie siedział na tej stronie tak i będzie się ładował przez dłuższy czas. Więc teraz jest dobry czas na krótką przerwę. Może zapomniałeś pójść do toalety lub chcesz filiżankę herbaty, zanim zaczniemy? Idź teraz, podczas gdy czekasz na to, bo będzie się kręcić tam przez chwilę.

Tylko szybko, podczas gdy czekamy, aż się załaduje, również przejdę do github.com/codespaces i po prostu pokażę, że to jest strona przeglądu, gdzie możesz zobaczyć wszystkie różne Codespaces, które obecnie masz uruchomione.

Widzisz, że mam tutaj jeden dla nextflow-io/training. Brak zmian, bo jeszcze nic w nim nie zrobiłem. Ilość zasobów, których używa, i widzisz, że w tej chwili się konfiguruje. Mogę tutaj przejść, kliknąć tę małą listę rozwijaną i kliknąć usuń. Więc jeśli przypadkowo skonfigurujesz wiele Codespaces i nie używasz niektórych, możesz usunąć stare i posprzątać.

Na koniec, jeszcze jeden sposób, aby się tam dostać. Jeśli przejdziemy do repozytorium GitHub. I to działa dla każdego repozytorium GitHub. Kliknij code. Masz polecenia do klonowania repozytorium lokalnie. I jest zakładka o nazwie Codespaces. I znowu, możesz utworzyć nowy, i możesz zobaczyć te, które już działają.

Więc znowu, jeśli zapomnisz, jak utworzyłeś swój Codespace, zawsze możesz wrócić do niego w ten sposób.

## Interfejs VS Code

Okej, budowniczy skończył i teraz zaczyna ładować GitHub Codespaces. Nie zawsze trwa to tak długo, więc nie martw się. To tylko pierwszy raz, gdy tworzysz Codespace. Jeśli wskoczysz z powrotem do tego, który już istnieje, jest znacznie szybciej.

Nie bądź zbyt niecierpliwy, jeśli to pierwszy raz, jeszcze się nie skończył, mimo że zaczyna nam dawać interfejs.

Ale podczas gdy czekamy na ostatnie rzeczy do skonfigurowania, po prostu przeprowadzę Cię przez interfejs na wypadek, gdybyś był trochę nieznajomy z VS Code.

Po pierwsze, jest pasek boczny czatu dla rzeczy AI, którego nie potrzebujemy. Więc zamknę to, pozbędę się tego i zwolnię trochę miejsca.

Po lewej stronie mamy eksplorator plików, który pokazuje nam wszystkie pliki w repozytorium Git, które jest przestrzenią roboczą, którą utworzyliśmy. Zauważ, to nie są pliki lokalne. To wszystko jest na zdalnym serwerze, gdzie pracujemy. Możesz przeciągać i upuszczać pliki lokalne i rzeczy, ale w większości nie będziemy o tym myśleć dzisiaj. Po prostu pracujemy czysto zdalnie.

Są inne narzędzia na tym pasku bocznym, na przykład wyszukiwanie. Więc możesz przeszukać wszystkie pliki w repozytorium za jednym razem. I gdybyśmy robili pracę rozwojową nad repozytorium szkoleniowym, moglibyśmy zrobić integrację z kontrolą źródła z Git i debugowaniem i innymi rzeczami.

Inne rzeczy to, jest główne okno edycji kodu tutaj na górze, które właśnie załadowało podgląd readme, który jest dla materiału szkoleniowego. Więc w tym przypadku przegląda markdown, ale normalnie to będzie edytor kodu.

A potem poniżej tego mamy terminal, gdzie będziemy uruchamiać wszystkie nasze polecenia i wchodzić w interakcję bezpośrednio z Nextflow'em.

Wszystko w Codespace jest wstępnie zainstalowane, więc polecenie Nextflow jest już tam i tak dalej i tak dalej.

Okej. Gdy dojdziesz tak daleko, powinno być mniej więcej gotowe. Widzisz teraz, że pobrało serwer języka Nextflow i skonfigurowało dla nas kilka rozszerzeń w VS code, w tym rozszerzenie Nextflow, które będzie przydatne. Więc mogę to zamknąć i mogę zamknąć README.md.

I teraz widzisz, że mam więcej po lewej stronie. Jestem trochę przybliżony tutaj, ale jeśli oddalam, widzisz, że jeden z przycisków mówi Nextflow z ikoną Nextflow. I ma tam kilka fajnych rzeczy do eksplorowania projektu i rzeczy, do których wrócimy później.

Okej. Na wypadek, gdybyś kiedykolwiek stracił którykolwiek z tych paneli, te przyciski w prawym górnym rogu są naprawdę przydatne i po prostu pokazują i ukrywają rzeczy. Więc to pokazuje i ukrywa Eksplorator, pokazuje i ukrywa terminal na dole. I tak dalej.

Będę używał ich dość dużo, ponieważ jestem bardzo przybliżony, więc staram się pomóc Ci zobaczyć cały tekst na moim ekranie, i więc przydatne jest móc przejść na pełny ekran z terminalem, a potem ukryć go, gdy patrzymy na kod. Ale przez większość czasu możesz po prostu mieć wszystkie te rzeczy otwarte w tym samym czasie.

Okej, co jeszcze zobaczyć? Niewiele więcej. Zauważ, że Nextflow, jak mówię, jest zainstalowany. Więc mogę wpisać "nextflow -version" i powinno pojawić się, którą wersję mamy zainstalowaną.

Jest też inne rzeczy zainstalowane tutaj. Na końcu każdego rozdziału mamy zestaw pytań quizowych, na przykład, na stronie internetowej. I możesz również zrobić je w terminalu, jeśli chcesz, wpisując quiz.

Jest kilka innych skrótów klawiszowych, których będę używał, na wypadek, gdybyś był ciekawy. Na przykład, właśnie wtedy nacisnąłem cmd+K na moim Macu i to wyczyściło terminal, aby pozbyć się całego poprzedniego wyjścia. Więc to jest miłe, aby utrzymać rzeczy w czystości. Jeśli widzisz, że to robię, tak to robię.

A także, jeśli jesteś nowy w terminalu, pamiętaj, że możesz użyć tab do autouzupełniania, co będę robił dużo, aby autouzupełniać ścieżki.

Więc widzę po lewej stronie, że jest folder o nazwie Hello Nextflow, przez który będziemy pracować. Jeśli zrobię "ls", aby wyświetlić pliki, mogę zrobić "hel", nacisnąć tab, autouzupełnia. I więc to jest bardzo szybki sposób na uzupełnianie ścieżek.

## Otwieranie tylko folderu Hello Nextflow

Okej. To świetnie. Jest jednak dużo rzeczy w tym repozytorium.

Są wszystkie pliki do generowania strony internetowej i jest wiele różnych kursów tutaj, i możesz to zrobić z tego korzenia i po prostu kliknąć w folder "Hello Nextflow". Ale miło jest faktycznie skupić się wyłącznie na tym.

Możesz ustawić to jako swoją przestrzeń roboczą z kilkoma kliknięciami tutaj i ustawianiem katalogu projektu i rzeczy. Ale najłatwiejszym sposobem jest wpisanie code, które jest poleceniem CLI do uruchamiania VS Code, a potem "hello-nextflow".

To otworzy nową kartę przeglądarki i możesz zamknąć starą. I wygląda dokładnie tak samo. Ale teraz widzisz, że jesteśmy w tym podkatalogu i wszystkie inne pliki są niewidoczne, i mamy czystszą konfigurację.

Widzisz tutaj, że również bieżący katalog roboczy jest teraz w folderze Hello Nextflow. Więc miło i czysto. Nie musimy się martwić o bycie w złym miejscu. Okej.

## Nowa składnia Nextflow na 2026

Jest jedna specjalna rzecz, którą muszę wspomnieć w tym momencie. Teraz, na początku 2026 roku, zaczynamy wprowadzać różne funkcje do Nextflow'a, a jedną z dużych nowych jest nowy parser składni języka wewnątrz Nextflow'a.

Zasadniczo silnik, który czyta Twoje pliki Nextflow i rozumie to, dla środowiska uruchomieniowego. Są pewne zmiany w składni i naprawdę ważne jest, abyś używał Nextflow'a z włączonym poprawnym parserem składni.

Potrzebujesz dwóch rzeczy do tego. Potrzebujesz aktualnej wersji Nextflow'a i musisz upewnić się, że jest włączony.

Jeśli zrobię "nextflow -version" ponownie, zobaczysz, że Codespaces działa z 25.10.2 i 25.10 to minimalna wersja, aby móc używać tych rzeczy.

Jeśli używasz 26.04, który dla mnie jeszcze nie wyszedł, ale wkrótce wyjdzie. Wtedy to będzie uruchamiać nowy parser składni domyślnie i nie musisz robić nic więcej.

Ale jeśli uruchamiasz 25.10, musisz włączyć ścisły parser składni, jak to się nazywa, lub parser składni v2.

To jest zrobione za pomocą zmiennej środowiskowej. Jest już ustawiona w Codespaces, więc nie musisz nic robić. Ale jeśli uruchamiasz lokalnie, musisz to ustawić, i mogę to zweryfikować, robiąc "echo $NXF_SYNTAX_PARSER", i powinno być ustawione na v2.

Więc jeśli uruchamiasz lokalnie, po prostu zrób "export NXF_SYNTAX_PARSER=v2". Proste jak to. Ale pamiętaj, aby to zrobić, bo inaczej zobaczysz dziwne rozbieżności i błędy, gdy będziemy iść dalej.

Jeśli jesteś w ogóle niepewny co do którejkolwiek z tych rzeczy wokół wersji Nextflow'a i parsera składni, po pierwsze, pamiętaj, nie musisz się martwić, jeśli jesteś w Codespaces. Wszystko powinno być prawidłowo skonfigurowane. Ale po drugie, jeśli przejdziesz do materiału szkoleniowego Nextflow, jeśli zejdziesz w dół, porozmawiaj o wymaganiach wersji, jest link tutaj, który zabiera Cię w dół do strony pomocy wokół eksplorowania wersji, i to w pewnym sensie przechodzi przez to wszystko szczegółowo.

Warto to przeczytać, jeśli masz chwilę, bo pomaga wyjaśnić, jakie są niektóre z różnych terminów, które możesz usłyszeć, gdy zaczniesz używać Nextflow'a. Rzeczy takie jak DSL1, DSL2, parser składni jeden, parser składni dwa i tak dalej. Więc warto po prostu sprawdzić to i to powtarza część tego, co właśnie powiedziałem.

Jest również naprawdę przydatne, jeśli wcześniej pisałeś kod Nextflow i wracasz na odświeżenie. Mówi Ci o niektórych rzeczach, które się zmieniają i linkuje Cię do części dokumentacji Nextflow, która mówi Ci, jak zaktualizować Twój kod Nextflow.

## Pliki kursu

Okej. Ostatnia rzecz, aby zapoznać się, to po prostu zobaczyć pliki, które są w tym katalogu. Możesz albo spojrzeć na pasek boczny, albo często w materiale szkoleniowym używamy polecenia tree, -L, które jest liczbą poziomów do zajrzenia. Powiemy dwa, i jeśli zrobię to na pełnym ekranie, zobaczysz, że to zasadniczo odzwierciedla dokładnie to, co widzimy na pasku bocznym tam, ale wyklucza ukryte pliki, które zaczynają się od kropki.

Więc pliki \*.nf, oznacza Nextflow. Więc to są pliki skryptów Nextflow, i jest plik startowy tutaj dla każdego z różnych rozdziałów materiału szkoleniowego, który otworzymy i zbadamy, a potem edytujemy.

Będziemy zmieniać te pliki w miarę postępów, i więc pod koniec każdego rozdziału pliki powinny wyglądać mniej więcej tak samo jak początek rozdziału dla następnego. Ale dajemy Ci te różne pliki, więc zawsze możesz zacząć od nowa i nie martwić się zbytnio o psucie składni.

Jeśli musisz porównać do czegoś, co zdecydowanie powinno działać. Możesz sprawdzić w folderze solutions, i to jest jak stan końcowy dla każdego z rozdziałów, więc możesz porównać to, co napisałeś, z tym, co tam jest.

Jest katalog data. Ma tylko plik greetings.csv, którego użyjemy jako przykładowych danych wejściowych w części kursu, i rzeczy takie jak plik konfiguracyjny i niektóre parametry, które opiszemy później w kursie.

## Podsumowanie

Okej, więc teraz mam nadzieję, że wszystko działa. Twój ekran wygląda tak samo jak mój i rozumiesz, jak dostać się do wszystkiego i czym są wszystkie różne pliki.

Jeśli przewiniesz w dół do dołu strony na Pierwsze kroki, małe pole wyboru powinno powiedzieć, że rozumiem, co robię. Moje środowisko jest uruchomione i ustawiłeś, Twój katalog roboczy prawidłowo do folderu "Hello Nextflow".

Jeśli zaznaczyłeś wszystkie te i wyglądają na zielone. Możemy kontynuować do następnego wideo i następnego rozdziału, który jest częścią pierwszą. Hello World. Do zobaczenia za chwilę.
