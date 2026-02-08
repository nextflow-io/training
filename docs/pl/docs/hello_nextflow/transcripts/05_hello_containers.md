# Część 5: Hello Containers - Transkrypcja wideo

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Tłumaczenie wspomagane sztuczną inteligencją - [dowiedz się więcej i zaproponuj poprawki](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

<div class="video-wrapper">
  <iframe width="560" height="315" src="https://www.youtube.com/embed/Xqr--bKEN9U?si=QinuAnFwFj-Z8CrO&amp;list=PLPZ8WHdZGxmWKozQuzr27jyMGqp9kElVK&amp;cc_load_policy=1" title="YouTube video player" frameborder="0" allow="accelerometer; autoplay; clipboard-write; encrypted-media; gyroscope; picture-in-picture; web-share" referrerpolicy="strict-origin-when-cross-origin" allowfullscreen></iframe>
</div>

!!!note "Ważne uwagi"

    Ta strona zawiera tylko transkrypcję. Pełne instrukcje krok po kroku znajdziesz w [materiałach kursu](../05_hello_containers.md).

    Numery sekcji w transkrypcji mają charakter orientacyjny i mogą nie obejmować wszystkich numerów sekcji z materiałów.

## Powitanie i wprowadzenie

Cześć i witaj ponownie w Hello Nextflow. To jest część piąta, czyli Hello Containers. W tej części kursu omówimy szczegółowo sposób enkapsulacji wymagań programowych dla pipeline'u, aby osoby uruchamiające go nie musiały martwić się instalacją oprogramowania.

Jeśli pracujesz w bioinformatyce tak długo jak ja, być może pamiętasz to, co często nazywam złymi starymi czasami, kiedy chciałeś uruchomić czyjś pipeline albo odtworzyć jego pracę, a spędzałeś godziny lub dni próbując zainstalować wszystkie różne narzędzia programowe, których używali, w odpowiednich wersjach, próbując je skompilować na swoim komputerze - to był koszmar. Naprawdę trudny.

Jeśli pracowałeś na HPC, mogłeś używać modułów środowiskowych, gdzie administratorzy systemowi próbowali instalować oprogramowanie dla Ciebie, co było w porządku, ale wciąż niedoskonałe.

Ale teraz mamy lepsze sposoby, żeby to zrobić. Nextflow ma wbudowane wsparcie dla różnych technologii kontenerów programowych. Docker to najpopularniejsza. Tego będziemy używać dzisiaj. Dobrze działa w Codespaces, dobrze na Twoim lokalnym komputerze i dobrze w chmurze.

Są też Singularity czy Apptainer, które są bardzo popularne w systemach HPC i w zasadzie działają dokładnie tak samo. Albo Podman, Shifter - jest kilka innych, które są bardzo podobne.

Jest jeszcze jeden dodatkowy, który jest trochę podobny, ale niekoniecznie taki sam, który Nextflow obsługuje - to Conda. Nextflow może zarządzać środowiskami Conda dla Ciebie na poziomie każdego procesu, co jest dużo lepsze niż zarządzanie własnymi środowiskami Conda. I znowu, może być dostarczany wraz z pipeline'em.

Zaczniemy ten rozdział od krótkiego omówienia technologii kontenerowych, Dockera i tego, jak działają. Pierwszą połowę zrobimy ręcznie w Dockerze, żebyś zrozumiał, co dzieje się pod maską i jak to działa. Bo to jest naprawdę ważne, aby zrozumieć, co robi Nextflow i jak rozumieć, co robi Twój workflow'a, kiedy jest wykonywany.

Więc przeskoczmy do naszego Codespaces. Ponownie wszystko posprzątałem, ale jeśli przejdziesz do Hello Containers, powinieneś zobaczyć, że wszystkie nasze skrypty są tam takie same jak na końcu rozdziału o modułach. Mamy nasze różne moduły tutaj, które stworzyłem w katalogu modules.

Wciąż tam są. Muszą tam być, żeby to mogło działać. Workflow'a i wyjście są takie same, z wyjątkiem tego, że zmieniliśmy ścieżkę publikacji wyjść na Hello Containers, żeby Twoje pliki trafiały do tego katalogu.

Możemy uruchomić to teraz, żeby sprawdzić, czy działa, jeśli chcesz, albo możemy przejść do terminala.

## 1. Ręczne użycie kontenera

Będziemy używać Dockera do zarządzania naszymi kontenerami, mogę sprawdzić, czy jest zainstalowany w moim Codespaces, wykonując "docker -v", co pokazuje mi wersję, która jest zainstalowana i wszystko, i że działa prawidłowo.

Kontenery i Docker mają dwie koncepcje, które są naprawdę ważne. Jedna nazywa się obraz (image), a druga - kontener (container). Obraz to migawka, jeśli chcesz, całego systemu plików, którego będziesz używać, a kontener to działające środowisko. Tworzysz kontener przy użyciu obrazu.

Kiedy już jesteś w tym kontenerze, zwykle działa on jak cały system operacyjny. Jest odcięty od świata zewnętrznego. Jest oddzielony od wszystkiego innego, i to jest dobra rzecz. Tak właśnie uzyskujemy tak dobrą powtarzalność z Nextflow'em.

Ponieważ dla zadań uruchamianych wewnątrz kontenera nie mają wpływu żadne pliki konfiguracyjne w Twoim lokalnym systemie. Żadne inne zewnętrzne wpływy - działają w swoim własnym małym piaskownicy. Pliki są wtedy tworzone w bardzo, bardzo powtarzalny sposób, ponieważ używasz tych samych podstawowych bibliotek, wszystkich tych samych zależności, dokładnie tego samego oprogramowania dla każdej osoby uruchamiającej w każdym innym środowisku obliczeniowym. Co szczerze uważam za fantastyczne i niesamowite, że to działa. I nawet do dziś wciąż mnie to zachwycca, że jest to możliwe.

## 1.1. Pobranie obrazu kontenera

Więc zamierzamy wypróbować użycie niektórych obrazów Dockera, a Docker, gdy uruchamiasz go na swoim systemie, ma rejestr dockera na Twoim komputerze, lub w tym przypadku w Codespace, który śledzi wszystkie różne obrazy, które zostały pobrane i używane w przeszłości, oraz różne warstwy, z których są zbudowane.

Możemy zobaczyć, jakie obrazy mamy lokalnie w Dockerze, wykonując "docker image ls". W tym przypadku widzisz kilka obrazów Dockera tutaj, które wszystkie mają związek z konfiguracją tego Codespaces. Wszystkie związane z kontenerami deweloperskimi i tym podobne. Więc nie musisz się nimi zbytnio martwić, ale gdy będziemy dodawać więcej obrazów i pobierać je, w miarę postępu tego kursu, możesz sprawdzić tę listę i zobaczysz, że lokalny rejestr śledzi wszystkie te rzeczy, które pobraliśmy.

Ale zamierzamy pobrać nowy, wykonując "docker pull". To mówi Dockerowi, żeby pobrał nowy obraz z internetu.

Następnie wpisujemy URI dla tego kontenera. Może to być obraz dockera, który zbudowałeś lokalnie, a następnie wypchnąłeś do internetu. Może to być obraz, który ktoś inny zrobił. Jest wiele, wiele różnych sposobów tworzenia obrazów Dockera, ale prawdopodobnie jednym z najprostszych sposobów jest zlecenie tego na zewnątrz i pozwolenie komuś innemu to zrobić za Ciebie.

A tego, czego będziemy używać w tym tutorialu, to usługa od Seqera o nazwie Seqera Containers.

Seqera Containers jest całkowicie darmowa i używa oprogramowania open source, które opracowaliśmy, o nazwie Wave, które zostało zbudowane do zarządzania kontenerami w sposób uzupełniający Nextflow'a. Obsługuje wiele typowych przypadków użycia, z którymi mamy do czynienia przy Nextflow'ie.

Jest bardzo powszechne, że oprogramowanie, którego potrzebujemy, jest spakowane w Condzie, w kanałach Bioconda lub conda-forge, albo innych bardziej specjalistycznych kanałach. A Wave i Seqera Containers są naprawdę dobre w budowaniu obrazów z tego.

Mogę przejść do tego interfejsu webowego i będziemy majstrować przy pakiecie o nazwie "cowpy". Więc wpisuję nazwę pakietu, którego chcę. Przeszukuje, znalazł go w indeksie pakietów Pythona, więc mogę go użyć. Lub jeśli poczekam trochę dłużej, przeszukuje biocondę i conda-forge. Widzisz, mogę określić dowolny kanał conda tutaj. Jeśli chcesz znaleźć kanał Nvidia lub cokolwiek innego, to też powinno działać.

A potem mogę określić, czy chcę, żeby zbudował dla mnie obraz dockera czy obraz singularity, a także jaką architekturę CPU chcę użyć. Więc amd64 lub arm64.

A gdy wyniki biocondy są wylistowane, mogę teraz zobaczyć wszystkie różne wersje, które są również dostępne. Zamierzam to wstawić. I teraz mógłbym dalej szukać i zdobywać więcej pakietów z Condy, jeśli chcę, i komponować ten kontener tak, jak chcę, ale chcę tylko tego jednego. Więc kliknę Get Container.

Teraz ktoś inny już wcześniej prosił o ten sam kontener i jest zwracany z rejestru, więc po prostu dostajemy go natychmiast. Ale gdyby nikt nigdy nie poprosił o ten pakiet programowy lub tę kombinację pakietów programowych, Wave i Seqera Containers zbudowałyby go w locie dla nas.

Możemy skopiować ten URL i możemy również zobaczyć szczegóły budowy. To pokazuje nam, co usługa zrobiła na backendzie. Stworzyła plik środowiska condy. Plik dockera, a następnie to jest on, uruchamiający proces budowania dockera. Uruchomił też skanowanie, skanowanie bezpieczeństwa, więc możesz zobaczyć jakiekolwiek CVE. I mówi Ci, kiedy to zostało stworzone.

Wave i Seqera Containers mogą zrobić o wiele więcej niż to, ale to jest rodzaj prostego przypadku użycia, który jest najpopularniejszy. Powinienem powiedzieć, że te obrazy są hostowane przez co najmniej pięć lat. Możesz więc wbudować te URL-e w swoje pipeline'y i wiedzieć, że nie znikną w najbliższym czasie.

Więc mam mój URL dla mojego obrazu dockera dla cowpy.

Mogę teraz wykonać "docker pull" tego URL-a i pobierze wszystkie różne warstwy i pobierze ten obraz, aby był dostępny dla mnie lokalnie.

## 1.2. Użycie kontenera do uruchomienia cowpy jako jednorazowe polecenie

Dobra, teraz spróbujmy go faktycznie użyć. Więc teraz zamierzam użyć polecenia "docker run" zamiast "docker pull", i zamierzam użyć flagi "--rm", która po prostu mówi Dockerowi, żeby zamknął ten kontener, gdy skończy robić to, o co go prosiłem. A potem wpisuję identyfikator kontenera, który jest po prostu URI.

A potem na końcu określam polecenie, które chcę, żeby Docker uruchomił wewnątrz kontenera wygenerowanego z tego obrazu. Po prostu powiem cowpy, co jest nazwą narzędzia zainstalowanego z Conda Forge, które jest dostępne wewnątrz obrazu.

Nacisnę enter i proszę. Uruchomiliśmy cowpy na systemie. Mamy małą krowę dającą nam trochę informacji.

Zauważ, że cowpy nie jest zainstalowane w moim lokalnym systemie. Więc jeśli uruchomię to po prostu bez całego tego Dockera, mówi command not found. Więc to pobrało obraz. Stworzyło kontener przy użyciu Dockera, a następnie weszło do tego kontenera i uruchomiło to polecenie dla nas i zwróciło nam wyjście z powrotem do naszego terminala. Bardzo, bardzo fajne.

## 1.3. Interaktywne użycie kontenera do uruchomienia cowpy

Dobra, teraz pójdziemy o krok dalej i uruchomimy ten kontener interaktywnie i rozejrzymy się trochę, żebyśmy mogli zobaczyć, co się dzieje wewnątrz kontenera.

Więc jeśli wrócę i wezmę moje polecenie run i zamierzam pozbyć się cowpy na końcu, ponieważ tak naprawdę nie chcę uruchamiać cowpy. Chcę uruchomić terminal Bash.

A potem zamierzam cofnąć się tutaj i zamierzam zrobić "-it", co oznacza Interactive i Terminal lub TTY, i zamierzam nacisnąć enter.

I teraz widzisz, że prompt, część przed tym, co piszę, się zmienił. To był prompt Codespaces, gdzie pokazywał katalog, a teraz mówi base i roots i tmp. Więc jestem teraz wewnątrz kontenera, a jeśli zrobię "ls", zobaczysz, że pliki, które widzę w tym katalogu, są inne niż pliki, które mam w moim workspace.

I w rzeczywistości nie mogę zobaczyć żadnych plików z mojego lokalnego workspace'u codespaces ani mojego lokalnego dysku wewnątrz kontenera Dockera. Środowisko uruchomieniowe kontenera dockera jest całkowicie izolowane i nie może zapisywać ani odczytywać żadnych plików z hosta systemu plików na zewnątrz.

Mogę jednak zobaczyć oprogramowanie, które jest zainstalowane wewnątrz kontenera i je uruchomić. Więc mogę uruchomić cowpy i możemy zobaczyć trochę więcej o tym, jak używać cowpy. Tutaj mogę zrobić "cowpy 'Hello World'" i to mówi mu, żeby faktycznie umieścił mój cytat wewnątrz małego dymku mowy. Możesz również uruchamiać różne typy krów, więc nie musi to być krowa. Możesz zrobić "-c". A jestem w Szwecji, więc wybiorę łosia. Bardzo ładny. Dałem mu poroże.

I jest cała masa różnych, z którymi możesz się pobawić, które są opisane w dokumentacji szkoleniowej.

## 1.3.4. Montowanie danych do kontenera

Dobra. Byłoby miło, gdybyśmy mogli uruchomić cowpy na plikach w naszym systemie plików.

Oczywiście nie jest zbyt przydatne mieć tylko kontener i żaden dostęp do niczego. Może być bezpieczny i powtarzalny, ale nie jest zbyt użyteczny.

Więc jak to robimy? Zamierzam wyjść z tego kontenera Dockera, wpisując exit, i możesz zobaczyć, że prompt mówi nam, że jesteśmy teraz z powrotem w naszym zwykłym Codespaces.

Zamierzam uruchomić to samo polecenie ponownie. Ale tym razem zamierzam dodać kilka dodatkowych flag z powrotem tutaj. Ważna to "-v", co oznacza montowanie woluminu, który jest jak w zasadzie część przestrzeni dyskowej.

"-v" przyjmuje dwie części: jest taki ciąg znaków, a potem dwukropek i ciąg znaków. Pierwsza część to lokalny system plików, który powinien być zamontowany w kontenerze. A potem druga część to gdzie to powinno skończyć wewnątrz kontenera.

Teraz po prostu chcę załadować cały mój lokalny system plików tutaj. Więc "." to aktualny katalog roboczy. Więc po prostu zrobię "." a potem ":", a potem zamierzamy umieścić to w nowym katalogu wewnątrz kontenera o nazwie "my_project". To może naprawdę być nazwane jak chcesz.

A potem zamierzam uruchomić ponownie.

W katalogu roboczym, gdzie jestem rzucony, który jest /tmp, plików tam nie ma. Ale jeśli zrobię "ls my_project", proszę: wszystkie te same pliki, które mieliśmy lokalnie w naszym Codespaces, są teraz dostępne wewnątrz kontenera w tej ścieżce.

To jest dostęp do odczytu i zapisu, więc mogę tworzyć nowe pliki w tym katalogu i pojawią się w moim hostowym systemie plików. Ten konkretny katalog wtedy zachowuje się dokładnie tak, jakbym był poza kontenerem, więc mogę teraz czytać i zapisywać i robić rzeczy.

## 1.3.5. Użycie zamontowanych danych

Dobra, udowodnijmy, że możemy to zrobić. Robię "cat /my_project/data/greetings.csv". Jeśli pamiętasz, zawartość tego pliku wygląda tak. Mogę teraz przekierować to do cowpy, a krowa wydrukuje różne wyjścia tego pliku w swoim małym dymku mowy, co jest trochę zabawne.

Więc widzisz, możemy teraz używać oprogramowania w kontenerze do interakcji z plikami w naszym systemie hosta.

Dobra, wyskoczmy z powrotem i zajmiemy się resztą materiału szkoleniowego.

## 2. Użycie kontenerów w Nextflow'ie

Więc to jest naprawdę fajne, używanie kontenerów. Mam nadzieję, że to ma sens. I widzisz wartość tych kontenerów i dlaczego jest to przydatne do uruchamiania oprogramowania analitycznego.

Ale jak robimy cały ten sam proces wewnątrz Nextflow'a? Nie chcemy uruchamiać mnóstwa poleceń Dockera sami. Chcemy po prostu pozwolić Nextflow'owi obsłużyć to wszystko dla nas.

Więc przejdźmy przez to. Zamierzamy dodać nowy proces do naszego pipeline'u, żeby uruchomić cowpy. Dobra, więc stwórzmy nowy moduł dla naszego nowego procesu. Więc idź do modules, nazwijmy to cowPy.nf, a potem zamierzam skopiować kod z materiału szkoleniowego tutaj.

Ale widzisz, proces jest bardzo prosty. Wygląda podobnie do tych, które robiliśmy do tej pory, mamy blok wejściowy ze ścieżką, która jest naszym plikiem wejściowym, a także wartością tutaj, żeby to był znak, więc moglibyśmy użyć łosia ponownie, jeśli chcemy.

A potem wyjście, które jest pojedynczym plikiem tutaj, ścieżką, a potem skryptem. I robimy to samo, co robiliśmy interaktywnie wewnątrz kontenera: robimy "cat", żeby odczytać plik wejściowy. Przekierowujemy tę zawartość do cowpy. Wybieramy konkretny znak na podstawie tego wejścia, zapisujemy do pliku wyjściowego o nazwie cowpy input file, który jest następnie echowany na wyjście.

Świetnie. Dołączmy to. Więc include \{ cowpy \} from "./modules/cowpy.nf", czy nazwałem to cowpy? Tak.

A potem wywołajmy nasz nowy proces tutaj w dół w głównym bloku workflow'a. Więc uruchommy cowpy. I weźmiemy nasz nowy proces cowpy i zamierzamy powiedzieć collectGreetings.out.

A potem jeśli pamiętasz, były dwa wyjścia dla tego modułu. Jedno nazywało się outfile, a drugie report. Rozszerzenie VS Code automatycznie podpowiada nam to i chcemy .outfile.

Zawsze możesz wskoczyć do tego procesu tutaj. Albo najedź na niego i powinno Ci szybko pokazać, jakie były wyjścia. I możemy też command click w to i otworzy plik modułu, jeśli chcesz zobaczyć bardziej szczegółowo.

Więc proszę. To jest outfile tam, i to jest ścieżka. Więc to będzie teraz plik wejściowy dla naszego procesu cowpy. Fantastycznie.

Teraz jeśli pamiętasz, proces cowpy ma dwa wejścia. Mieliśmy również kanał wartości dla znaku. Więc możemy dodać tutaj "params.character". Mogłem to zakodować na stałe, jeśli chciałem, ale zróbmy to opcją CLI, żebyśmy mogli zrobić dash, dash character.

Dobra. Teraz muszę zdefiniować parametr wejściowy, który właśnie wywołaliśmy, i nadać mu wartość domyślną. Więc character, String. I lubię łosia, więc zamierzam ustawić go na moose domyślnie.

Dobra, spróbujmy to uruchomić. Więc jeśli zrobię Nextflow run hello containers, zobaczymy, co się stanie.

Mogłem użyć dash resume, gdybym miał stare katalogi robocze kręcące się dookoła. I znowu, te pierwsze procesy byłyby w pamięci cache i byłoby trochę szybciej, ale powinno być w zasadzie to samo.

Teraz możemy od razu zobaczyć, że zgłosił błąd, gdy doszedł do naszego nowego procesu, mówi nam tutaj, że był błąd wykonywania procesu cowpy i zakończył się ze statusem wyjścia 127. To jest polecenie, które próbował uruchomić. Wygląda dobrze, wygląda, jak się spodziewaliśmy. Bierze nazwę tego pliku wyjściowego, która wygląda w porządku, uruchamia to ze znakiem moose i próbuje zapisać.

Ale widzisz, błąd polecenia tutaj mówi, że polecenie cowpy nie zostało znalezione. I to ma sens, ponieważ tak naprawdę nie powiedzieliśmy jeszcze Nextflow'owi, żeby użył kontenera. Po prostu daliśmy mu polecenie cowpy. A jak powiedziałem wcześniej, cowpy nie jest zainstalowane w naszym lokalnym systemie. Więc kiedy próbował to uruchomić, zawiodło.

## 2.3.1. Określenie kontenera dla cowpy

Musimy powiedzieć Nextflow'owi, że jest dostępny kontener i może go użyć. Więc jak to robimy?

Jeśli wskoczmy do naszego modułu tutaj, zamierzamy dodać nową deklarację na górze o nazwie "container". I zamierzamy ustawić to na ciąg znaków.

Teraz, jeśli pamiętasz, w Seqera Containers mogę skopiować ten URL i po prostu wrzucam to w cudzysłowy tutaj.

Teraz wracam i próbuję uruchomić to ponownie.

Zobaczmy, czy zadziała tym razem.

Niestety zawodzi dokładnie w ten sam sposób, mimo że teraz zdefiniowaliśmy kontener dla procesu do uruchomienia. Więc żeby użyć naszego obrazu dockera, musimy powiedzieć Nextflow'owi, żeby włączył użycie Dockera, gdy uruchamiamy workflow'a.

I zamierzamy to zrobić, tworząc nowy plik konfiguracyjny. Więc zamierzam powiedzieć touch nextflow.config.

To jest specjalna nazwa pliku, gdzie jeśli jest w katalogu roboczym, podczas gdy uruchamiam pipeline'a, zostanie automatycznie załadowany. Więc jeśli wejdę do tego pliku Nextflow dot config, widzisz, że w rzeczywistości już istnieje, o czym zapomniałem. I mamy docker.enabled już tutaj, ale jest ustawiony na false, co jest domyślne.

Więc jeśli zmienię to na equals True zamiast tego, docker.enabled. I jest dokumentacja referencyjna dla wszystkich tych zakresów konfiguracji w dokumentacji Nextflow'a. A także możesz zobaczyć, że gdy najadę z rozszerzeniem VS Code, wciąga dokumentację specyficzną dla tego i mówi mi, co to znaczy i jak to ustawić.

Więc teraz ustawiliśmy to na true, i jeśli uruchomię Nextflow ponownie, Nextflow będzie teraz wiedział, żeby pobrać ten obraz dockera dla nas, jeśli jeszcze nie mamy go lokalnie, a następnie wykonać ten proces z tym środowiskiem kontenera.

I możemy zobaczyć, że uruchomił się pomyślnie i mamy małą ptaszkę obok cowpy. Fantastycznie. Jeśli przejdę w górę i spojrzę w katalogu wyników, pliku tam jeszcze nie ma. I to dlatego, że wciąż musimy opublikować ten plik wyjściowy tak samo jak wszystkie inne.

Więc idziemy do bloku publikowania wewnątrz workflow'a, mówimy mycowpy equals cowpy.out.

A potem tutaj w dole w bloku wyjściowym, mycowpy, klamrowe nawiasy path. Ups. Hello containers. Mode, copy.

Jeśli teraz uruchomię ponownie, powinien działać dokładnie w ten sam sposób. Mogłem użyć dash resume i zapominam za każdym razem. A potem idę w górę i teraz mamy nowy plik utworzony o nazwie cowpy-COLLECTED, i tam jest mój łoś mówiący BONJOUR, HELLO, HOLA Fantastycznie.

Oczywiście mogę również teraz przekazać "--character". Jakie są różne opcje? Myślę, że jest indyk? Więc mogę użyć character Turkey. Zamierza działać dokładnie w ten sam sposób. Przegapiłem kolejną okazję, żeby użyć dash resume, i teraz jeśli załadujemy nasz plik i teraz mamy indyka. Fantastycznie.

## 2.3.4. Sprawdzenie, jak Nextflow uruchomił konteneryzowane zadanie

Dobra. Ostatnia mała rzecz. Uruchommy to polecenie ponownie, tym razem resume, i rzućmy szybkie oko do katalogu roboczego, żeby zobaczyć, co to jest, co Nextflow robi pod maską, żeby sprawić, że to wszystko dla nas działa.

Tym razem jest super szybko, chodźmy do tego katalogu roboczego, cd work/. Teraz jeśli pamiętasz, mamy kilka plików dot tutaj, a ten, który nas interesuje w tym przypadku, to ten, o którym powiedziałem, że prawie nigdy nie musimy patrzeć, nazywany .command.run.

Jeśli zrobię code dot command run, otworzy to w edytorze. I mogę przeszukać ten plik i jeśli przewinę w dół, powinienem zobaczyć Docker run. I widzisz, Nextflow wykonuje polecenie docker run dla nas, gdy Docker jest włączony w konfiguracji. Ma całą masę różnych flag i rzeczy tutaj, ale widzisz flagę "-v", której używaliśmy sami, gdy uruchamialiśmy. I widzisz, że montuje lokalny katalog workspace do kontenera, żeby kontener mógł uzyskać dostęp do naszych plików wejściowych i zapisać wyjścia. A potem na końcu uruchamia również .command.sh, który jest wygenerowanym skryptem, który ma w sobie polecenie cowpy.

I możesz zobaczyć, że Nextflow bierze logikę workflow'a, która jest tym, na czym nam naprawdę zależy, która jest specyficzna dla naszej analizy, i robi wszystkie sprytne rzeczy za kulisami, żeby sprawić, że Docker działa w naszym systemie.

I robi to w naprawdę przenośny sposób, żeby użytkownik końcowy pipeline'u mógł zamienić technologię, której używa: Docker, Singularity, Apptainer, Conda. To naprawdę nie ma znaczenia dla logiki pipeline'u, ale Nextflow obsłuży wszystkie podstawowe potrzeby infrastruktury, żeby działało wszędzie.

I to jest naprawdę super moc Nextflow'a. To powtarzalność i przenośność. A z Nextflow'em możesz faktycznie udostępnić swój workflow'a i inne osoby mogą go uruchomić w swoich systemach i po prostu zadziała.

To jest naprawdę, naprawdę trudna rzecz do zrobienia, a teraz wiesz, jak to zrobić również ze swoimi workflow'ami.

Dobra, to tyle na ten rozdział. Jeśli przejdziesz na koniec kursu, znajdziesz quiz ponownie o niektórych kontenerach. Mam nadzieję, że to wszystko miało sens. To naprawdę fajny sposób pracy z analizą. A jeśli jesteś nowy w kontenerach, mam nadzieję, że przekonałem Cię, że to jest droga do podążania i nigdy nie spojrzysz wstecz.

Ale na tym zróbcie sobie małą przerwę może, i dołączycie do mnie za kilka minut, żeby przejść przez finałową część sześć Hello Nextflow, która jest cała o konfiguracji.

Dziękuję bardzo.
