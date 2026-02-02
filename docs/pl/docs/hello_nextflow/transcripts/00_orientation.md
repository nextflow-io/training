# Orientacja - Transkrypcja Wideo

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Tłumaczenie wspomagane przez AI - [dowiedz się więcej i zasugeruj ulepszenia](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

<div class="video-wrapper">
  <iframe width="560" height="315" src="https://www.youtube.com/embed/G3CV-FcV-rc?si=nyLvwhrSB2m1NPc5&amp;list=PLPZ8WHdZGxmXiHf8B26oB_fTfoKQdhlik" title="YouTube video player" frameborder="0" allow="accelerometer; autoplay; clipboard-write; encrypted-media; gyroscope; picture-in-picture; web-share" referrerpolicy="strict-origin-when-cross-origin" allowfullscreen></iframe>
</div>

!!!note "Ważna uwaga"

    Ta strona zawiera tylko transkrypcję. Aby uzyskać pełne instrukcje krok po kroku, wróć do [materiałów kursu](../00_orientation.md).

## Powitanie

Cześć, witaj w Hello Nextflow. Nazywam się Phil Ewels. Jestem Product Managerem ds. Open Source w Seqera i bardzo się cieszę, że mogę dzisiaj poprowadzić Cię przez ten pierwszy kurs szkoleniowy Nextflow.

Przejdziemy przez podstawy Nextflow, wyjaśniając jak pisać i uruchamiać pipelines oraz jak je konfigurować.

Zbudujesz własny prosty wieloetapowy pipeline. Omówimy terminologię taką jak operatory i fabryki kanałów, a pod koniec kursu będziesz gotowy, aby rozpocząć budowanie własnych pipelineów bioinformatycznych.

Jeśli masz jakiekolwiek pytania, skontaktuj się z nami na community.seqera.io. Mamy tam bardzo aktywną społeczność Nextflow, jest sekcja poświęcona szkoleniom, więc po prostu daj nam znać, gdzie utknąłeś, a ktoś będzie w stanie pomóc.

Dobrze. Zaczynajmy.

## Strona Szkoleniowa

Wszystkie materiały szkoleniowe dla kursów Nextflow znajdują się na training.nextflow.io. Możesz przejść do niej w Swojej przeglądarce internetowej. Uruchom ją teraz i możemy się rozejrzeć.

Będę korzystał z wersji 2.1.1. Wprowadzamy małe aktualizacje i poprawki co jakiś czas, więc nie martw się, jeśli będzie trochę inne, ale jeśli materiały zbyt bardzo się rozjechały, zawsze możesz użyć tego selektora wersji na górze, aby wybrać dokładną wersję materiałów, którą będę omawiał.

Jeśli wolisz jasny motyw, możesz zmienić motyw strony tutaj.

Zobacz tłumaczenia tutaj, chociaż w momencie nagrywania to tak naprawdę tylko angielski obejmuje te nowe materiały.

I zobacz także cały kod źródłowy strony szkoleniowej i wszystko, z czym będziemy pracować, na GitHub.

Strona główna tutaj wymienia wszystkie różne kursy materiałów szkoleniowych, które mamy. Jeśli przewinę w dół, zobaczymy Nextflow dla początkujących z kursem Hello Nextflow, który tutaj wykonamy. Możesz zobaczyć wszystkie inne kursy, które również mamy i które działają w podobny sposób.

## Konfiguracja Środowiska

Tak naprawdę zamierzam zacząć od tego pierwszego na górze, który jest wspólny dla wszystkich kursów szkoleniowych i dotyczy konkretnie konfigurowania naszego środowiska.

Klikam go, przenosi mnie do tej sekcji i widzimy instrukcje dotyczące lokalnego programowania. Jeśli chcesz użyć Swojego własnego laptopa z własną kopią VS Code i własnymi instalacjami oprogramowania, lub tego, czego oczekujemy od większości ludzi, czyli użycia czegoś, co nazywa się GitHub Codespaces.

Codespaces to usługa świadczona przez GitHub, w której uruchamiają serwer WWW w chmurze, z którym możesz się połączyć. Ten serwer ma zainstalowany VS Code, gdzie możesz go uruchomić w przeglądarce internetowej, lub jeśli wolisz, połączyć go z lokalną instalacją VS Code. Wszystkie obliczenia, wszystkie pliki, cała edycja odbywa się zdalnie, co oznacza, że wszystkie potrzebne oprogramowanie jest preinstalowane i jest takie samo dla wszystkich.

## Tworzenie GitHub Codespace

Aby utworzyć codespace ze wszystkim, czego potrzebujemy, poszukaj przycisków w materiałach dokumentacji, które mówią "Open in GitHub Codespaces". Kliknę to teraz, otwieram to w nowej karcie. I pokazano mi tę stronę internetową. Teraz widzisz, że jest to wstępnie skonfigurowane z nextflow-io training.

Mogę po prostu kliknąć create new codespace. Ale tak naprawdę zalecamy użycie nieco większej maszyny do szkolenia Nextflow z czterema procesorami zamiast dwóch. Możesz zmienić, której wersji materiałów używa. Więc to domyślnie ustawia się na 2.1.1, ponieważ to jest wersja dokumentacji, z której poszedłem za linkiem. Ale mogę też ustawić to na konkretną gałąź repozytorium, jeśli chcę.

Teraz kliknę create codespace. I zacznie konfigurować dla mnie środowisko.

## Tworzenie Codespace

Teraz, za pierwszym razem, kiedy to robisz, zajmie to dość dużo czasu, więc teraz jest dobry moment, aby pójść zrobić sobie herbatę. Usiądź wygodnie, porozmawiaj z osobą, przy której siedzisz.

Jeśli jesteś zainteresowany, możesz kliknąć building codespace tutaj na dole, aby zobaczyć logi konfiguracji. I możesz zobaczyć tutaj, że pobiera obraz Docker ze wszystkim, czego potrzebuję i konfiguruje środowisko.

Teraz musisz czekać w ten sposób tylko za pierwszym razem, gdy tworzysz codespace. Jeśli przejdziesz do github.com/codespaces tutaj, zobaczysz wszystkie różne Codespaces, które masz otwarte. Oto ten, który właśnie utworzyłem. Następnym razem, gdy to zrobisz, możesz przejść tutaj i możesz wybrać poprzedni codespace i po prostu od razu do niego wrócić. I to jest znacznie, znacznie szybszy proces rozgrzania tego istniejącego środowiska. To również zachowa wszystkie zmiany, które wprowadzałeś w VS Code i w plikach, więc nie stracisz postępu, jeśli wyjdziesz i wrócisz.

Możesz kliknąć trzy kropki tutaj, aby wykonać inne akcje. Na przykład, jeśli skonfigurowałeś go z dwoma procesorami, a teraz chcesz cztery, możesz zmienić typ maszyny. Lub jeśli chcesz zacząć od zera i od nowa, możesz usunąć codespace.

## Wprowadzenie do VS Code

Dobrze, Codespaces zakończył konfigurowanie mojego środowiska i teraz pokazano mi VS Code w przeglądarce internetowej.

Jeśli jesteś przyzwyczajony do VS Code, będzie to bardzo znajome. Jeśli nie używałeś go wcześniej, jest całkiem proste. Jest kilka różnych części strony, o których musisz wiedzieć.

Tutaj po lewej stronie mamy pasek boczny. Widzisz Eksplorator ustawiony ze wszystkimi różnymi plikami w repozytorium GitHub z repozytorium szkoleniowego.

Te przyciski na dole po lewej stronie mogą być różnymi narzędziami. W pasku bocznym mogę przeszukiwać wszystkie pliki we wszystkich projektach. Mogę pracować z Git, mogę pracować z GitHub, wszystkimi różnymi rzeczami tego typu.

Na górze tutaj jest główne menu. Eksplorator plików jest tym, który będziemy mieć najczęściej tutaj, i możesz kliknąć prawym przyciskiem myszy dowolny z tych plików i zrobić normalne rzeczy, których oczekujesz. Może być konieczne kliknięcie przez niektóre ostrzeżenia takie jak to, gdzie wytnij, kopiuj i możesz również pobrać na Swoją lokalną maszynę.

Gdy codespace się ładuje, daje nam podgląd pliku markdown w tym głównym obszarze tutaj. To jest dokładnie to samo, co renderuje się na github.com. Mogę to zamknąć, a jeśli kliknę dwukrotnie ten plik Readme, zobaczysz, że otwiera go jako kod w edytorze kodu i tak jak z każdym innym plikiem, możemy edytować ten kod bezpośrednio.

Na końcu tutaj na dole mamy okno terminala. Patrzyłem na logi, gdy się budował, więc to jest to, co aktualnie pokazuje. Mogę również nacisnąć ten przycisk plus, aby rozpocząć nową sesję terminala. To nie działa na mojej maszynie. Pamiętaj, to działa w chmurze, i jeśli wykonam tree trzy do głębokości dwóch, zobaczysz wszystkie te same pliki tutaj, które były po lewej stronie.

## Pokazywanie tylko plików "hello-nextflow"

To repozytorium GitHub zawiera wszystkie różne zestawy szkoleniowe, nie tylko ten, który robimy. Więc jeśli chcesz, możesz skupić się tylko na katalogu Hello Nextflow. Jednym ze sposobów, aby to trochę uporządkować, jest przejście do menu file, a następnie add folder to workspace.

Klikamy to, przechodzimy do training. Hello nextflow i klikamy add. Odświeży Twój ekran. A potem w Eksploratorze mamy teraz dwa różne przestrzenie robocze, tę, którą mieliśmy wcześniej dla training i jedną z tylko Hello Nextflow.

Jeśli chcesz, możesz kliknąć prawym przyciskiem myszy na training i kliknąć remove folder from workspace, aby całkowicie usunąć go z paska bocznego.

Teraz mamy tylko pliki dla tego konkretnego kursu szkoleniowego z boku. Mogę ukryć to ostrzeżenie i teraz mogę zrobić to samo w terminalu tutaj i wykonać CD dla zmiany katalogu. Hello, Nextflow. I znowu mamy te same pliki tutaj, które są w pasku bocznym.

## Hello Nextflow: pliki

Patrząc na te pliki dla kursu Hello Nextflow.

Mamy kilka plików .nf, które są dla Nextflow, i jest jeden z tych plików dla każdego z rozdziałów kursu szkoleniowego. Będziemy pracować przez te pliki i modyfikować je w ćwiczeniach.

Mamy również plik nextflow.config, który ma tylko podstawowe ustawienia konfiguracji do uruchamiania Nextflow w tym środowisku, którymi nie musisz się martwić w tym momencie. Plik greetings.csv, którego będziemy używać do przetwarzania danych, który zostanie wprowadzony w następnej części tego kursu, oraz plik test-params.json, który będzie używany w części szóstej i możesz zignorować go na razie.

Te pliki Nextflow to tylko początek każdego ćwiczenia. Jeśli chcesz zobaczyć, jak powinny wyglądać, gdy są skończone, możesz przejść do katalogu solutions i tam są odpowiedzi dla każdej części kursu szkoleniowego, więc możesz zobaczyć działającą wersję tego, do czego dążysz.

## Otwieranie terminala

Jeśli w którymkolwiek momencie zamkniesz terminal i nie możesz sobie przypomnieć, jak wrócić, nie martw się. Te przyciski na górze po prawej otwierają i zamykają różne panele w przestrzeni roboczej. Kliknij więc ten dla dolnego panelu i pojawi się ponownie. I upewnij się, że masz tutaj wybrany terminal. Możesz również kliknąć ten przycisk tutaj, strzałkę po prawej stronie terminala, aby zrobić pełny ekran.

Zobaczysz, jak robię to dość często, ponieważ mam VS Code powiększony, abyś mógł czytać tekst. W zależności od rozmiaru ekranu, możesz potrzebować lub nie potrzebować tego robić. To samo dotyczy minimalizowania panelu bocznego.

Dobrze. To wystarczy dla środowiska. Myślę, że jesteśmy gotowi, aby zacząć. Dołącz do mnie w następnym filmie dla rozdziału pierwszego.

[Następna transkrypcja wideo :octicons-arrow-right-24:](01_hello_world.md)
