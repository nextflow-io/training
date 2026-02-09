# Część 5: Hello Containers - Transkrypcja wideo

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Tłumaczenie wspomagane przez AI - [dowiedz się więcej i zaproponuj poprawki](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

<div class="video-wrapper">
  <iframe width="560" height="315" src="https://www.youtube.com/embed/Xqr--bKEN9U?si=QinuAnFwFj-Z8CrO&amp;list=PLPZ8WHdZGxmWKozQuzr27jyMGqp9kElVK&amp;cc_load_policy=1" title="YouTube video player" frameborder="0" allow="accelerometer; autoplay; clipboard-write; encrypted-media; gyroscope; picture-in-picture; web-share" referrerpolicy="strict-origin-when-cross-origin" allowfullscreen></iframe>
</div>

!!!note "Ważne uwagi"

    Ta strona zawiera tylko transkrypcję. Pełne instrukcje krok po kroku znajdziesz w [materiałach kursu](../05_hello_containers.md).

    Numery sekcji pokazane w transkrypcji mają charakter orientacyjny i mogą nie obejmować wszystkich numerów sekcji w materiałach.

## Powitanie i wprowadzenie

Cześć i witaj ponownie w Hello Nextflow. To część piąta zatytułowana Hello Containers. W tej części kursu porozmawiamy o tym, jak hermetyzować wymagania dotyczące oprogramowania dla pipeline'u, aby osoby uruchamiające pipeline nie musiały myśleć o instalowaniu oprogramowania.

Jeśli pracujesz w bioinformatyce tak długo jak ja, możesz pamiętać to, co często nazywam złymi starymi czasami, kiedy chciałeś uruchomić czyjś pipeline lub odtworzyć jego pracę, spędzałeś godziny lub dni próbując zainstalować wszystkie różne narzędzia programistyczne, których używali, w tych samych wersjach, próbując je skompilować na swojej maszynie, i to był koszmar. To było naprawdę trudne.

Jeśli pracowałeś na HPC, mogłeś używać modułów środowiskowych, gdzie administratorzy systemowi próbowali instalować dla Ciebie oprogramowanie, co było w porządku, ale wciąż niedoskonałe.

Ale teraz mamy lepsze sposoby, aby to zrobić. Nextflow ma wbudowane wsparcie dla różnych technologii kontenerów oprogramowania. Docker jest najpopularniejszy. To właśnie jego użyjemy dzisiaj. Dobrze działa w Codespaces. Dobrze działa na Twoim lokalnym komputerze i dobrze działa w chmurze.

Ale także Singularity lub Apptainer, które są bardzo popularne w systemach HPC i działają praktycznie w dokładnie taki sam sposób. Lub Podman, Shifter, jest wiele innych, które są bardzo podobne.

Jest jeszcze jeden dodatkowy, który jest trochę podobny, ale nie do końca, który Nextflow obsługuje - to Conda. Nextflow może zarządzać dla Ciebie środowiskami Conda dla każdego procesu osobno, co jest znacznie lepsze niż zarządzanie własnymi środowiskami Conda. I znowu, może być dostarczane z pipeline'em.

Zaczniemy ten rozdział od krótkiej rozmowy o technologiach kontenerów i Dockerze oraz o tym, jak działają. Pierwszą połowę zrobimy ręcznie w Dockerze, abyś zrozumiał, co dzieje się pod maską i jak to działa. Bo to naprawdę ważne, aby zrozumieć, co robi Nextflow i jak rozumieć, co robi Twój workflow podczas wykonywania.

Więc przeskoczmy do naszych Codespaces. Teraz wszystko posprzątałem, ale jeśli przejdziesz do Hello Containers, powinieneś zobaczyć, że wszystkie nasze skrypty i wszystko jest tam takie samo jak na końcu rozdziału o modułach. Mamy więc nasze różne moduły tutaj, które stworzyłem w katalogu modules.

Nadal tam są. Muszą tam być, aby mogło działać. Workflow i wyjście są takie same, z wyjątkiem tego, że zmieniliśmy ścieżkę publikowania wyjścia na Hello Containers, aby Twoje pliki trafiały do tego katalogu.

Możemy to teraz uruchomić, aby sprawdzić, czy działa, jeśli chcesz, lub możemy przejść do terminala.

## 1. Użyj kontenera 'ręcznie'

Będziemy używać Dockera do zarządzania naszymi kontenerami i mogę sprawdzić, czy jest zainstalowany w moich Codespaces, robiąc "docker -v", co pokazuje mi wersję, która jest zainstalowana i wszystko, i że działa poprawnie.

Teraz kontenery i Docker mają dwie koncepcje, które są naprawdę ważne. Jedna nazywa się obraz (image), a druga kontener (container). Obraz to migawka, jeśli chcesz, całego systemu plików, którego będziesz używać, a kontener to działające środowisko. Tworzysz kontener używając obrazu.

Gdy jesteś w tym kontenerze, zazwyczaj działa on jak cały system operacyjny. Jest odcięty od świata zewnętrznego. Jest oddzielony od wszystkiego innego, i to jest dobre. Tak właśnie uzyskujemy tak dobrą powtarzalność z Nextflow'em.

Ponieważ dla zadań uruchamianych wewnątrz kontenera nie są one skażone przez żadne pliki konfiguracyjne w Twoim lokalnym systemie. Żadne inne zewnętrzne wpływy, działają w swoim małym piaskownicy. Pliki są następnie tworzone w bardzo, bardzo powtarzalny sposób, ponieważ używasz tych samych podstawowych bibliotek, wszystkich tych samych zależności, dokładnie tego samego oprogramowania dla każdej osoby uruchamiającej w każdym innym środowisku obliczeniowym. Co szczerze mówiąc uważam za fantastyczne i niesamowite, że to działa. I nawet do dziś wciąż mnie to zadziwia, że to jest możliwe.

## 1.1. Pobierz obraz kontenera

Więc zamierzamy wypróbować używanie niektórych obrazów Dockera i Docker, gdy uruchamiasz go w swoim systemie, ma rejestr dockera na Twoim komputerze, lub w tym przypadku w Codespaces, który śledzi wszystkie różne obrazy, które zostały pobrane i używane w przeszłości, oraz różne warstwy, z których są zbudowane.

Możemy zobaczyć, jakie obrazy mamy lokalnie w Dockerze, robiąc "docker image ls". I w tym przypadku widzisz, że jest tutaj kilka obrazów Dockera, które wszystkie mają związek z konfiguracją tych Codespaces. Wszystkie związane z kontenerami deweloperskimi i tym podobnymi rzeczami. Więc nie musisz się nimi zbytnio przejmować, ale w miarę jak dodajemy więcej obrazów i je pobieramy, w miarę postępu tego kursu, możesz sprawdzić tę listę i zobaczysz, że lokalny rejestr śledzi wszystkie te rzeczy, które pobraliśmy.

Ale zamierzamy pobrać nowy, robiąc "docker pull". I to mówi Dockerowi, aby pobrał nowy obraz z sieci.

Następnie wpisujemy URI dla tego kontenera. Teraz może to być obraz dockera, który zbudowałeś lokalnie, a następnie wypchnąłeś do internetu. Może to być obraz, który ktoś inny zrobił. Jest wiele, wiele, wiele różnych sposobów tworzenia obrazów Dockera, ale prawdopodobnie jednym z najprostszych sposobów jest zlecenie tego na zewnątrz i pozwolenie komuś innemu to zrobić za Ciebie.

I to, czego będziemy używać w tym samouczku, to usługa od Seqera nazywana Seqera Containers.

Teraz Seqera Containers jest całkowicie darmowa i używa oprogramowania open source, które rozwijamy, nazywanego Wave, które zostało zbudowane do zarządzania kontenerami w sposób komplementarny do Nextflow'a. I obsługuje wiele typowych przypadków użycia, z którymi mamy do czynienia z Nextflow'em.

Bardzo często oprogramowanie, którego potrzebujemy, jest spakowane w Condzie, w kanałach Bioconda lub conda-forge lub innych bardziej specyficznych dla dziedziny kanałach. I Wave oraz Seqera Containers są naprawdę dobre w budowaniu obrazów z tego.

Więc mogę przejść do tego interfejsu webowego i będziemy bawić się pakietem o nazwie "cowpy". Więc wpisuję nazwę pakietu, którego chcę. Wyszukuje, znalazł go w indeksie pakietów Pythona, więc mogę go użyć. Lub jeśli poczekam trochę dłużej, przeszukuje biocondę i conda-forge. I widzisz, mogę określić dowolny kanał Condy tutaj. Więc jeśli chcesz znaleźć kanał Nvidia lub cokolwiek innego, to również powinno działać.

A następnie mogę określić, czy chcę, aby zbudował dla mnie obraz dockera czy obraz singularity, a także jaką architekturę CPU chcę użyć. Więc amd64 lub arm64.

I gdy wyniki z biocondy są wymienione, mogę teraz zobaczyć również wszystkie różne dostępne wersje. Zamierzam to dodać. I teraz mógłbym kontynuować wyszukiwanie i uzyskać więcej pakietów z Condy, jeśli chcę, i skomponować ten kontener, jak chcę, ale chcę tylko tego jednego. Więc kliknę Get Container.

Teraz ktoś inny już wcześniej poprosił o ten sam kontener i jest zwracany z rejestru, więc po prostu dostajemy go natychmiast. Ale gdyby nikt inny nigdy nie poprosił o ten pakiet oprogramowania lub tę kombinację pakietów oprogramowania, Wave i Seqera Containers zbudowałyby go w locie dla nas.

Możemy skopiować ten URL i możemy również zobaczyć szczegóły budowy. I to pokazuje nam, co usługa zrobiła na backendzie. Stworzyła plik środowiska condy. Plik dockera, a następnie to jest on, uruchamiający proces budowy dockera. Uruchomił również skan, skan bezpieczeństwa, więc możesz zobaczyć wszelkie CVE. I mówi Ci, kiedy to zostało utworzone.

Wave i Seqera Containers mogą zrobić znacznie więcej niż to, ale to jest rodzaj prostego przypadku użycia, który jest najczęstszy. I powinienem powiedzieć, że te obrazy są hostowane przez co najmniej pięć lat. Więc możesz wbudować te URL-e w swoje pipeline'y i wiedzieć, że nie znikną w najbliższym czasie.

Więc mam mój URL dla mojego obrazu dockera dla cowpy.

Mogę teraz zrobić "docker pull" tego URL-a, i pobierze wszystkie różne warstwy i pobierze ten obraz, aby był dostępny dla mnie lokalnie.

## 1.2. Użyj kontenera do uruchomienia cowpy jako jednorazowego polecenia

Okej, teraz spróbujmy faktycznie go użyć. Więc teraz zamierzam użyć polecenia "docker run" zamiast "docker pull", i zamierzam użyć flagi "--rm", która po prostu mówi Dockerowi, aby zamknął ten kontener, gdy skończy robić to, o co go poprosiłem. A następnie wpisuję identyfikator kontenera, który jest po prostu URI.

A następnie na końcu określam polecenie, które chcę, aby Docker uruchomił wewnątrz kontenera wygenerowanego z tego obrazu. Po prostu powiem cowpy, co jest nazwą narzędzia, które jest zainstalowane z Conda Forge, które jest dostępne wewnątrz obrazu.

Zamierzam nacisnąć enter i proszę. Uruchomiliśmy cowpy w systemie. Mamy małą krowę dającą nam pewne informacje.

Teraz zauważ, że cowpy nie jest zainstalowane w moim lokalnym systemie. Więc jeśli uruchomię to po prostu bez całego tego Dockera, mówi, polecenie nie znalezione. Więc to pobrało obraz. Stworzyło kontener używając Dockera, a następnie weszło do tego kontenera i uruchomiło to polecenie dla nas i dało nam wyjście z powrotem do naszego terminala. Bardzo, bardzo fajne.

## 1.3. Użyj kontenera do uruchomienia cowpy interaktywnie

Okej, pójdziemy teraz o krok dalej i uruchomimy ten kontener interaktywnie i pobawimy się trochę, abyśmy mogli zobaczyć, co dzieje się wewnątrz kontenera.

Więc jeśli wrócę i wezmę moje polecenie run i zamierzam pozbyć się cowpy na końcu, bo właściwie nie chcę uruchamiać cowpy. Chcę uruchomić terminal Bash.

A następnie zamierzam wrócić tutaj i zamierzam zrobić "-it", co oznacza Interactive i Terminal lub TTY, i zamierzam nacisnąć enter.

I teraz widzisz, że prompt, część przed tym, co wpisuję, się zmienił. To był prompt Codespaces, gdzie pokazywał katalog, a teraz mówi base i root i tmp. Więc jestem teraz wewnątrz kontenera, i jeśli zrobię "ls", zobaczysz, że pliki, które widzę w tym katalogu, są inne niż pliki, które mam w moim workspace.

I w rzeczywistości nie mogę zobaczyć żadnych plików z mojego lokalnego workspace'u Codespaces ani mojego lokalnego dysku wewnątrz kontenera Dockera. Środowisko uruchomieniowe kontenera dockera jest całkowicie odizolowane i nie może zapisywać ani odczytywać żadnych plików z systemu plików hosta na zewnątrz.

Mogę jednak zobaczyć oprogramowanie, które jest zainstalowane wewnątrz kontenera i je uruchomić. Więc mogę uruchomić cowpy i możemy zobaczyć trochę więcej o tym, jak używać cowpy. Tutaj mogę zrobić "cowpy 'Hello World'" i to mówi, mówi mu, aby faktycznie umieścił mój cytat w małym dymku mowy. I możesz również uruchamiać różne typy krów, więc nie musi to być krowa. Możesz zrobić "-c". I jestem w Szwecji, więc wybiorę łosia. Bardzo ładnie. Dałem mu poroże.

I jest cała masa różnych, z którymi możesz się pobawić, które możesz zobaczyć opisane w dokumentach szkoleniowych.

## 1.3.4. Zamontuj dane w kontenerze

Okej. Byłoby miło, gdybyśmy mogli uruchomić cowpy na plikach w naszym systemie plików.

Oczywiście nie jest zbyt przydatne mieć tylko kontener i brak dostępu do czegokolwiek w ogóle. Może to być bezpieczne i powtarzalne, ale nie jest zbyt przydatne.

Więc jak to robimy? Zamierzam wyjść z tego kontenera Dockera, wpisując exit, i widzisz, że prompt mówi nam, że jesteśmy teraz z powrotem w naszych zwykłych Codespaces.

I zamierzam uruchomić to samo polecenie ponownie. Ale tym razem zamierzam dodać kilka dodatkowych flag tutaj z powrotem. I ważna to "-v", co oznacza montowanie woluminu, który jest jak w zasadzie część, część przestrzeni dyskowej.

"-v" przyjmuje dwie części: jest jak ciąg znaków, a następnie dwukropek i ciąg znaków. I pierwsza część to lokalny system plików, który powinien być zamontowany w kontenerze. A następnie druga część to miejsce, gdzie to powinno się znaleźć wewnątrz kontenera.

Teraz chcę po prostu załadować cały mój lokalny system plików tutaj. Więc "." to bieżący katalog roboczy. Więc po prostu zrobię "." a następnie ":", a następnie zamierzamy umieścić to w nowym katalogu wewnątrz kontenera o nazwie "my_project". To naprawdę mogłoby być nazwane czymkolwiek.

A następnie zamierzam uruchomić ponownie.

W katalogu roboczym, gdzie jestem zrzucony, który jest /tmp, pliki tam nie są. Ale jeśli zrobię "ls my_project", mamy to: wszystkie te same pliki, które mieliśmy lokalnie w naszych Codespaces, są teraz dostępne wewnątrz kontenera w tej ścieżce.

To jest dostęp do odczytu i zapisu, więc mogę tworzyć nowe pliki w tym katalogu i pojawią się one w moim systemie plików hosta. Ten konkretny katalog zachowuje się wtedy dokładnie tak, jakbym był poza kontenerem, więc mogę teraz odczytywać i zapisywać i robić rzeczy.

## 1.3.5. Użyj zamontowanych danych

Okej, po prostu udowodnijmy, że możemy to zrobić. Robię "cat /my_project/data/greetings.csv". Jeśli pamiętasz, zawartość tego pliku wygląda tak. Mogę teraz przekazać to do cowpy i krowa wydrukuje różne wyjścia tego pliku w swoim małym dymku mowy, co jest trochę zabawne.

Więc widzisz, możemy teraz używać oprogramowania w kontenerze do interakcji z plikami w naszym systemie hosta.

Okej, wróćmy i przejdziemy do reszty materiału szkoleniowego.

## 2. Użyj kontenerów w Nextflow

Więc to naprawdę fajne używanie kontenerów. Mam nadzieję, że to ma sens. I widzisz wartość tych kontenerów i dlaczego to jest przydatne do uruchamiania oprogramowania analitycznego.

Ale jak robimy cały ten sam proces wewnątrz Nextflow'a? Nie chcemy uruchamiać mnóstwa poleceń Dockera sami. Chcemy po prostu pozwolić Nextflow'owi obsłużyć to wszystko za nas.

Więc przejdźmy przez to. Zamierzamy dodać nowy proces do naszego pipeline'u, aby uruchomić cowpy. Okej, więc stwórzmy nowy moduł dla naszego nowego procesu. Więc wchodzimy do modules, nazwijmy to cowPy.nf, a następnie zamierzam skopiować kod z materiału szkoleniowego tutaj.

Ale widzisz, proces jest bardzo prosty. Wygląda bardzo podobnie do tych, które robiliśmy do tej pory, mamy blok input ze ścieżką, która jest naszym plikiem wejściowym, a także wartością tutaj, więc to będzie znak, więc moglibyśmy użyć łosia ponownie, jeśli chcemy.

A następnie wyjście, które jest pojedynczym plikiem tutaj, ścieżką, a następnie skryptem. I robimy to samo, co robiliśmy interaktywnie wewnątrz kontenera: robimy "cat", aby odczytać, plik wejściowy. Przekazujemy tę zawartość do cowpy. Wybieramy konkretny znak na podstawie tego wejścia, zapisujemy do pliku wyjściowego o nazwie cowpy input file, który jest następnie echowany do wyjścia.

Świetnie. Dołączmy to. Więc include \{ cowpy \} from "./modules/cowpy.nf", czy nazwałem to cowpy? Tak.

A następnie wywołajmy nasz nowy proces tutaj w głównym bloku workflow'u. Więc uruchommy cowpy. I weźmiemy nasz nowy proces cowpy i zamierzamy powiedzieć collectGreetings.out.

A następnie, jeśli pamiętasz, były dwa wyjścia dla tego modułu. Jedno nazywało się outfile, a jedno report. Rozszerzenie VS Code automatycznie sugeruje nam te i chcemy .outfile.

Zawsze możesz wskoczyć do tego procesu tutaj. Albo najedź na niego i powinno pokazać Ci szybko, jakie były wyjścia. I możemy również kliknąć command w niego i otworzy plik modułu, jeśli chcesz zobaczyć bardziej szczegółowo.

Więc proszę. To jest outfile tam, i to jest ścieżka. Więc to będzie teraz plik wejściowy dla naszego procesu cowpy. Fantastycznie.

Teraz, jeśli pamiętasz, proces cowpy ma dwa wejścia. Mieliśmy również kanał wartości dla znaku. Więc możemy dodać tutaj "params.character". Mogłem to zakodować na stałe, jeśli chciałem, ale zróbmy to opcją CLI, abyśmy mogli zrobić dash, dash character.

Dobrze. Teraz muszę zdefiniować parametr wejściowy, który właśnie wywołaliśmy i nadać mu wartość domyślną. Więc character, String. I lubię łosia, więc zamierzam ustawić go na moose domyślnie.

Dobrze, spróbujmy to uruchomić. Więc jeśli zrobię Nextflow run hello containers, zobaczymy, co się stanie.

Mogłem użyć dash resume, gdybym miał stare katalogi robocze kręcące się, i znowu, te pierwsze procesy byłyby, buforowane i byłoby to trochę szybsze, ale powinno być w zasadzie to samo.

Teraz widzimy od razu, że rzucił błąd, gdy dotarł do naszego nowego procesu, mówi nam tutaj, że wystąpił błąd podczas wykonywania, procesu cowpy i zakończył się ze statusem wyjścia 127. To jest polecenie, które próbował uruchomić. Wygląda dobrze, wygląda, jak się spodziewaliśmy. Bierze tę nazwę pliku wyjściowego, która wygląda w porządku, uruchamia go ze znakiem moose i próbuje zapisać.

Ale widzisz, błąd polecenia tutaj mówi, że polecenie cowpy nie zostało znalezione. I to ma sens, ponieważ właściwie nie powiedzieliśmy jeszcze Nextflow'owi, aby użył kontenera. Po prostu daliśmy mu polecenie cowpy. I jak powiedziałem wcześniej, cowpy nie jest zainstalowane w naszym lokalnym systemie. Więc kiedy próbował to uruchomić, zawiodło.

## 2.3.1. Określ kontener dla cowpy

Musimy powiedzieć Nextflow'owi, że jest dostępny kontener i może go użyć. Więc jak to robimy?

Jeśli wskoczymy do naszego modułu tutaj, zamierzamy dodać nową deklarację na górze o nazwie "container". I zamierzamy następnie ustawić to na ciąg znaków.

Teraz, jeśli pamiętasz, w Seqera Containers mogę skopiować ten URL i po prostu upuszczam go w cudzysłowach tutaj.

Teraz wróć i spróbuj uruchomić to ponownie.

Zobaczmy, czy tym razem zadziała.

Niestety zawodzi w dokładnie taki sam sposób, mimo że teraz zdefiniowaliśmy kontener dla procesu do uruchomienia. Więc aby użyć naszego obrazu dockera, musimy powiedzieć Nextflow'owi, aby włączył użycie Dockera, gdy uruchamiamy workflow.

I zamierzamy to zrobić, tworząc nowy plik konfiguracyjny. Więc zamierzam powiedzieć touch nextflow.config.

To jest specjalna nazwa pliku, gdzie jeśli jest w katalogu roboczym, gdy uruchamiam pipeline, zostanie załadowany automatycznie. Więc jeśli wejdę do tego pliku Nextflow dot config, widzisz, że faktycznie już istnieje, o czym zapomniałem. I mamy docker.enabled tutaj już, ale jest ustawiony na false, co jest domyślne.

Więc jeśli zmienię to na equals True zamiast tego, docker.enabled. I są dokumenty referencyjne dla wszystkich tych zakresów konfiguracji w dokumentach Nextflow'a. A także widzisz, że gdy najadę z rozszerzeniem VS Code, pobiera dokumenty specyficzne dla tego i mówi mi, co to znaczy i jak to ustawić.

Więc teraz ustawiliśmy to na true, i jeśli uruchomię Nextflow ponownie, Nextflow będzie teraz wiedział, aby pobrać ten obraz dockera dla nas, jeśli jeszcze go nie mamy lokalnie, a następnie wykonać ten proces z tym środowiskiem kontenera.

I więc widzimy, że uruchomił się pomyślnie i mamy małą fajkę obok cowpy. Fantastycznie. Jeśli wejdę i spojrzę w katalogu results, plik jeszcze tam nie jest. I to dlatego, że wciąż musimy, opublikować ten plik wyjściowy tak samo jak, wszystkie inne.

Więc idziemy do bloku published w workflow'ie, mówimy mycowpy equals cowpy.out.

A następnie tutaj w bloku output, mycowpy, klamrowe nawiasy path. Ups. Hello containers. Mode, copy.

Jeśli uruchomię teraz ponownie, powinno uruchomić się w dokładnie taki sam sposób. Mogłem użyć dash resume i zapominam za każdym razem. A następnie wchodzę i teraz mamy nowy plik utworzony o nazwie cowpy-COLLECTED, i jest mój łoś mówiący BONJOUR, HELLO, HOLA Fantastycznie.

Teraz oczywiście mogłem również teraz przekazać "--character". Jakie są różne opcje? Myślę, że jest indyk? Więc mogę użyć character Turkey. Zamierza uruchomić się w dokładnie taki sam sposób. Przegapiłem kolejną okazję do użycia dash resume, i teraz jeśli załadujemy nasz plik i teraz mamy indyka. Fantastycznie.

## 2.3.4. Sprawdź, jak Nextflow uruchomił zadanie w kontenerze

Okej. Ostatnia mała rzecz. Po prostu szybko uruchommy to polecenie ponownie, resume tym razem, i rzućmy szybkie spojrzenie w katalogu roboczym, aby zobaczyć, co to jest, co Nextflow robi pod maską, aby wszystko to dla nas działało.

Tym razem jest super szybko, wejdźmy do tego katalogu roboczego, cd work/. Teraz, jeśli pamiętasz, mamy kilka plików kropkowych tutaj i ten, który nas interesuje w tym przypadku, to ten, o którym powiedziałem, że prawie nigdy nie musimy patrzeć, nazywany .command.run.

Jeśli zrobię code dot command run, otworzy się w edytorze. I mogę przeszukać ten plik i jeśli przewinę w dół, powinienem zobaczyć Docker run. I widzisz, Nextflow robi polecenie docker run dla nas, gdy Docker jest włączony w konfiguracji. Ma całą masę różnych, flag i rzeczy tutaj, ale widzisz flagę "-v", której użyliśmy sami, gdy uruchamialiśmy. I widzisz, że montuje lokalny, katalog workspace do kontenera, aby kontener mógł uzyskać dostęp do naszych plików wejściowych i zapisać wyjścia. A następnie na końcu uruchamia również .command.sh, który jest wygenerowanym skryptem, który ma polecenie cowpy w środku.

I więc widzisz, że Nextflow bierze logikę workflow'u, która jest rzeczą, na której nam naprawdę zależy, która jest specyficzna dla naszej analizy, i robi wszystkie sprytne rzeczy za kulisami, aby Docker działał w naszym systemie.

I robi to w naprawdę przenośny sposób, tak że użytkownik końcowy pipeline'u może zamienić technologię, której używa: Docker, Singularity, Apptainer, Conda. To naprawdę nie ma znaczenia dla logiki pipeline'u, ale Nextflow obsłuży wszystkie podstawowe potrzeby infrastruktury, tak że działa wszędzie.

I to jest naprawdę supermocy Nextflow'a. To powtarzalność i przenośność. I z Nextflow'em możesz faktycznie udostępnić swój workflow i inne osoby mogą go uruchomić w swoich systemach i po prostu zadziała.

To naprawdę, naprawdę trudna rzecz do zrobienia, a teraz wiesz, jak to zrobić również ze swoimi workflow'ami.

Okej, to wszystko na ten rozdział. Jeśli przejdziesz na koniec kursu, znajdziesz, quiz ponownie o niektórych kontenerach. Mam nadzieję, że to wszystko miało sens. To naprawdę fajny sposób pracy z analizą. I jeśli jesteś nowy w kontenerach, mam nadzieję, że przekonałem Cię, że to jest droga do pójścia, i nigdy nie będziesz oglądać się za siebie.

Ale z tym, zrób sobie małą przerwę może, i dołącz do mnie za kilka minut, aby przejść przez ostatnią część szóstą Hello Nextflow, która jest cała o konfiguracji.

Dziękuję bardzo.
