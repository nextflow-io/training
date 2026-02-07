# Część 5: Hello Containers - Transkrypcja

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Tłumaczenie wspomagane przez AI - [dowiedz się więcej i zasugeruj ulepszenia](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

<div class="video-wrapper">
  <iframe width="560" height="315" src="https://www.youtube.com/embed/5PyOWjKnNmg?si=QinuAnFwFj-Z8CrO&amp;list=PLPZ8WHdZGxmXiHf8B26oB_fTfoKQdhlik" title="YouTube video player" frameborder="0" allow="accelerometer; autoplay; clipboard-write; encrypted-media; gyroscope; picture-in-picture; web-share" referrerpolicy="strict-origin-when-cross-origin" allowfullscreen></iframe>
</div>

!!!note "Ważne informacje"

    Ta strona zawiera tylko transkrypcję. Pełne instrukcje krok po kroku znajdziesz w [materiale szkoleniowym](../05_hello_containers.md).

    Numery sekcji pokazane w transkrypcji mają charakter orientacyjny i mogą nie zawierać wszystkich numerów sekcji w materiałach.

## Witamy

Cześć, witamy w piątej części kursu szkoleniowego Hello Nextflow.

Ten rozdział nosi tytuł Hello Containers. Porozmawiamy o tym, jak Nextflow integruje się z narzędziami takimi jak Docker i Singularity, aby używać kontenerów oprogramowania do dostarczania narzędzi użytkownikom Twojego pipeline'u.

Oznacza to, że kiedy ludzie uruchamiają Twój pipeline, nie muszą sami instalować wszystkich różnych narzędzi. Nextflow zrobi to za nich.

Kontenery to niezwykle potężna technologia, kluczowa dla powtarzalności i łatwości użycia. Zaczniemy od krótkiego wprowadzenia do samych kontenerów, ręcznego uruchomienia kilku poleceń docker, a następnie użyjemy tych samych kontenerów w naszym pipeline'ie Nextflow'a.

Dobra. Zaczynajmy.

Tak jak poprzednio, zacznijmy od załadowania materiałów szkoleniowych. Wejdź na training.nextflow.io. Hello Nextflow, Rozdział 5, Hello Containers.

Przejdę do mojego środowiska Codespaces, a po lewej stronie widzimy hello containers dot nf.

Tak jak poprzednio, jest to ten sam skrypt, którym zakończyliśmy poprzedni rozdział 4, więc powinien wyglądać znajomo.

Mamy nasze parametry wiersza poleceń do określenia pliku wejściowego i nazwy partii. Dołączamy nasze trzy moduły i mamy nasz workflow, w którym uruchamiamy trzy procesy.

## 0. Rozgrzewka: Uruchom hello-containers.nf

Możesz uruchomić ten workflow ponownie i sprawdzić, czy produkuje oczekiwane wyniki. Na razie zamknę go i przejdę do terminala.

## 1. Użyj kontenera 'ręcznie'

Na początku tego rozdziału zrobimy podsumowanie technologii kontenerów. Jeśli jesteś bardzo zaznajomiony z docker, singularity lub innymi technologiami kontenerów, potraktuj to jako odświeżenie wiedzy lub pomiń tę część całkowicie.

Nextflow obsługuje wiele różnych typów technologii kontenerów. Obejmuje to Docker, Singularity, Podman, Shifter, Charliecloud i inne.

W tym szkoleniu skupimy się na Docker. Jest on preinstalowany w code spaces i jest jedną z najpopularniejszych technologii kontenerów, szczególnie jeśli rozwijasz na własnym komputerze.

Jeśli pracujesz w środowisku akademickim na współdzielonym HPC, możesz zauważyć, że dostępny jest Singularity, a nie Docker. To w porządku. Wszystkie koncepcje są dokładnie takie same. Kilka ręcznych poleceń jest różnych, ale jeśli rozumiesz Docker, zrozumiesz także singularity.

W rzeczywistości Singularity jest również zainstalowany w środowisku Code Spaces. Więc jeśli chcesz, możesz spróbować wykonać te same zadania używając Singularity zamiast Docker.

Dobra, czym jest technologia kontenerów? Ideą Docker jest to, że może pobrać obraz ze zdalnego źródła, ściągnąć go na lokalną maszynę, a następnie utworzyć kontener na podstawie tego obrazu.

Ten działający kontener jest trochę jak maszyna wirtualna działająca na Twoim komputerze. Jest odizolowany od Twojego środowiska i jest wstępnie zapakowany z systemem operacyjnym i zestawem dostępnego oprogramowania.

## 1.1. Pobierz obraz kontenera

Składnia, której potrzebujemy do pobrania istniejącego obrazu to "docker pull". Wpiszę to w moim terminalu, ale teraz potrzebujemy obrazu, z którym będziemy się bawić.

Możesz budować obrazy samodzielnie. Możesz je znaleźć w publicznych rejestrach takich jak Docker Hub czy quay.io. Ale naprawdę dobrym sposobem na szybkie uzyskanie obrazów jest użycie Seqera Containers.

To darmowa usługa społecznościowa, którą stworzyliśmy w 2024 roku i którą możesz używać bez logowania czy czegokolwiek innego.

Jeśli przejdziesz na seqera.io/containers lub klikniesz containers na górze, zobaczysz interfejs wyszukiwania i możesz wpisać nazwę dowolnego narzędzia dostępnego w Conda lub Python Package Index.

Domyślnie przeszukuje kanały Bioconda i Conda Forge, ale możesz poprzedzić dowolny kanał Conda, jeśli chcesz.

Dla zabawy użyjmy cowpy. Wpiszę cowpy. Daje mi wyniki z Python Package Index i Conda Forge. Kliknę to, aby dodać do mojego kontenera. Mógłbym dodać wiele pakietów, gdybym chciał. Wybieram Docker, wybieram linux/amd64 i klikam Get Container.

To buduje dla mnie obraz na żądanie, jeśli nie został jeszcze utworzony, i daje mi URL, który mogę skopiować.

Jeśli jesteś zainteresowany, możesz kliknąć view Build Details, co zabierze Cię na stronę pokazującą plik środowiska conda, który został użyty, oraz kompletny dziennik budowania wraz z wynikami skanowania bezpieczeństwa.

Jeśli wrócę do moich code spaces, mogę teraz wkleić tę nazwę kontenera i nacisnąć enter.

Docker teraz pobiera wszystkie różne warstwy w tym obrazie kontenera i teraz mówi nam, że ten obraz jest dostępny do użycia.

## Pobieranie obrazu Singularity

Jeśli używasz singularity, proces jest w zasadzie taki sam. Wybieramy nasze pakiety obrazów, wybieramy cowpy. Teraz wybieramy Singularity zamiast Docker i klikamy Get Container. To daje nam URL obrazu używający oras://. Lub jeśli wolisz, możesz użyć https:// zaznaczając to pole. Kopiuję ten URL. Teraz przechodzę do Code Spaces. W rzeczywistości mamy zainstalowany Apptainer w tej przestrzeni, który jest taki sam jak Singularity, ale są do siebie aliasowane. Więc zrobię apptainer pull, a następnie nazwę to cowpy sif, ale możesz nazwać to jak chcesz. Wklejam URL. I to pobierze dla mnie ten obraz.

Mógłbym zrobić ls -lh i zobaczyć cowpy.sif

Singularity różni się od Docker tym, że singularity przechowuje wszystkie obrazy w płaskich plikach, podczas gdy Docker ma rejestr, gdzie przechowuje wszystkie warstwy oddzielnie na Twojej maszynie hosta i ma działającego demona, aby śledzić to wszystko.

## 1.2. Użyj kontenera do uruchomienia cowpy jako jednorazowego polecenia

Dobra, wróćmy do Docker. Możemy teraz spróbować uruchomić ten obraz, który utworzyliśmy, robiąc docker run.

Zrobię dash dash rm, co po prostu wykonuje jednorazowe wykonanie obrazu. I wkleję URL obrazu. A na końcu kończysz to poleceniem, które chcesz uruchomić.

Obraz, który wygenerowaliśmy, miał zainstalowane cowpy, więc spróbujmy cowpy.

Proszę bardzo. Uruchomił nasze polecenie. Nie mam cowpy zainstalowanego lokalnie. Możesz zobaczyć, że jeśli spróbuję to uruchomić, nie istnieje. Jednak w tym poleceniu uruchomiłem to używając Docker i poprawnie wygenerowało to wyjście.

## 1.3. Użyj kontenera do uruchomienia cowpy interaktywnie

Możemy pójść dalej, jeśli chcemy, i uruchomić kontener interaktywnie i rozejrzeć się w środku. Ponownie robię "docker run dash dash rm". Teraz zrobię dash it, co mówi Docker, że chcemy interaktywny terminal. Ponownie robię URL obrazu, a tym razem, zamiast robić cowpy, zrobię bin bash, ponieważ polecenie, które chcemy uruchomić, to bash.

To zabiera nas do tego działającego kontenera i możesz zobaczyć, że prompt się teraz zmienił.

Jeśli zrobię LS slash, możesz zobaczyć, że katalogi tutaj są inne.

Jeśli otworzę drugi terminal tutaj po prawej stronie, który po prostu działa w GitHub Code Spaces i zrobię LS slash, widzisz, że mamy katalogi takie jak workspaces i temp, podczas gdy tutaj w Docker jest inaczej.

Więc to środowisko jest całkowicie oddzielne w Docker i odizolowane od mojego środowiska hosta. To dobra rzecz, ponieważ to izoluje wykonanie tego polecenia w obrazie Docker i utrzymuje to powtarzalne między różnymi ludźmi na różnych systemach hosta.

Jeśli chcesz użyć danych z Twojego systemu hosta w obrazie Docker, musisz jawnie zamontować to w kontenerze.

Za chwilę to zrobimy.

## 1.3.2. Uruchom żądane polecenia narzędzia

Najpierw jednak sprawdźmy, czy możemy uruchomić cowpy. Ponownie, polecenie jest teraz dostępne bezpośrednio w wierszu poleceń i możemy zacząć robić bardziej złożone rzeczy i przekazywać argumenty. Hello containers i zamiast krowy, zróbmy pingwina tux. Zobaczmy, co jeszcze mamy.

Zróbmy cheese. Wspaniale. Co powiesz na Dragon and Cow? Całkiem dobre.

## 1.3.3. Wyjdź z kontenera

Dobra. Nie mogę zrobić wiele więcej, ponieważ nie mam żadnych danych w tym kontenerze. Więc wyjdźmy z tego działającego obrazu i zobaczmy, czy możemy zamontować jakieś dane w kontenerze. Mogę to zrobić robiąc control D lub wpisując exit. Dobra, jestem teraz z powrotem w moim zwykłym GitHub code space.

## 1.3.4. Zamontuj dane w kontenerze

Aby zamontować jakieś dane w kontenerze Docker, muszę użyć dash V. Więc wezmę moje poprzednie polecenie docker, wrócę na początek i zrobię dash v. Zrobię "." dla bieżącego lokalnego katalogu roboczego, a następnie dwukropek, aby powiedzieć, gdzie to powinno być zamontowane w katalogu hosta i zrobię slash data. Więc to montuje ten konkretny katalog w kontenerze w slash data.

Teraz jeśli zrobię LS slash, możemy zobaczyć, że mamy nowy katalog o nazwie data, a jeśli zrobię LS data, możesz zobaczyć wszystkie pliki, które mamy na pasku bocznym tutaj. Fantastycznie.

## 1.3.5. Użyj zamontowanych danych

Teraz możemy zacząć używać niektórych plików, które są w systemie hosta w obrazie Docker. Więc mogę powiedzieć cat data greetings csv. Jeśli pamiętasz, to jest nasz plik CSV z naszymi różnymi powitaniami sprzed chwili, i mogę przekazać to do cowpy. Fantastycznie. Teraz gdzieś docieramy.

Dobra. Wystarczy uruchamiania Docker interaktywnie. Mam nadzieję, że masz teraz poczucie, czym w przybliżeniu jest Docker i jak go używać zarówno do uruchamiania polecenia jednorazowo, jak i do interaktywnego używania obrazu. Jeśli używasz singularity, polecenia są wszystkie bardzo podobne, z wyjątkiem tego, że robisz rzeczy takie jak apptainer exec lub apptainer run, lub singularity exec lub singularity run.

## 2. Użyj kontenerów w Nextflow

Następnie wrócimy do naszego workflow'a Nextflow'a i zobaczymy, jak używać tej technologii w pipeline'ie Nextflow'a.

Zamknijmy terminal i otwórzmy ponownie Hello Containers.

## 2.1. Napisz moduł cowpy

Aby trzymać się naszego przykładu cowpy, stwórzmy nowy proces w naszym workflow, który używa cowpy. Przejdźmy do modules, utwórzmy nowy plik i nazwijmy go cowpy nf. Teraz trochę oszukam i skopiuję kod z materiału szkoleniowego i nacisnę save. I spójrzmy.

To prosty proces. Mam nadzieję, że teraz rozumiesz, jak wyglądają elementy składowe procesu. Mamy ponownie nasz publishDir, idący do results. Mamy dwa wejścia, plik wejściowy i ciąg znaków o nazwie character. Mamy wyjście cowpy input file i mamy skrypt, który wygląda dokładnie tak samo jak to, co uruchamialiśmy ręcznie wewnątrz naszego obrazu docker przed chwilą: cat do wydrukowania pliku, przekazując to do cowpy, mówiąc, którego typu postaci cowpy chcemy użyć, i wyprowadzając to do pliku wyjściowego, który przekazujemy jako wyjście tutaj.

## 2.2. Dodaj cowpy do workflow'a

Dobra, wróćmy do naszego workflow'a, zaimportujmy ten nowy proces. Więc cowpy from modules cowpy nf. Stwórzmy nowy parametr, abyśmy mogli określić, którego znaku chcemy. Powiedzmy Turkey domyślnie. A następnie wywołajmy ten nowy proces na końcu workflow'a,

cowpy. I użyjmy wyjścia tutaj z Collect Greetings. Więc collect greetings out, out file tutaj. A następnie potrzebujemy drugiego argumentu, którym są nowe params, które właśnie stworzyliśmy. params dot character.

## 2.2.4. Uruchom workflow, aby sprawdzić, czy działa

Dobra, zobaczmy, czy nasz nowy proces działa. Nextflow run hello containers. To powinno uruchomić te pierwsze trzy procesy, a następnie spróbować uruchomić cowpy na końcu.

Mamy błąd. To, co tutaj mówi, cowpy miało błąd i miało status wyjścia 127 i rzeczywiście, polecenie sh cowpy polecenie nie znalezione.

Nie powiedzieliśmy Nextflow'owi, że mamy dostępny obraz Docker dla cowpy, więc próbował uruchomić to na naszym systemie hosta, a nie mamy cowpy zainstalowanego na naszym systemie hosta, więc wywołało to błąd.

## 2.3. Użyj kontenera, aby to uruchomić

Więc to, co musimy zrobić, to musimy powiedzieć Nextflow'owi, że mamy dostępny kontener. Przejdźmy do naszego procesu cowpy i dodajmy nową dyrektywę na górze procesu o nazwie container.

Następnie znajdujemy nasz obraz, kopiujemy URL i umieszczamy to w ciągu znaków.

To samo w sobie nie wystarcza, ponieważ pipeline X Flow może mieć kilka sposobów określania oprogramowania. Mógłbym również zrobić conda conda-forge cowpy, na przykład. I Nextflow musi wiedzieć, której z tych technologii chcesz użyć.

## 2.3.2. Włącz użycie Docker przez plik nextflow.config

Więc aby uruchomić z włączonym Docker, wyprzedzimy się trochę i użyjemy pliku Nextflow config, którym zajmiemy się bardziej szczegółowo w następnym rozdziale. Możesz zobaczyć w tym katalogu, że mamy plik o nazwie Nextflow Config, i tutaj już masz docker.enabled False.

Zmienimy to na True, aby włączyć Docker, a następnie możemy spróbować ponownie uruchomić workflow.

## 2.3.3. Uruchom workflow z włączonym Docker

Nextflow run hello containers nf i tym razem cowpy uruchomił się pomyślnie. Spójrzmy w Results. cowpy collected test i oto nasz Turkey. Wspaniale.

Więc w tle tam Nextflow wiedział, że ma dostępny kontener dla tego procesu.

Pobrał obraz i uruchomił dla nas polecenia.

## 2.3.4. Sprawdź, jak Nextflow uruchomił zadanie w kontenerze

Jeśli jesteś ciekawy, możemy faktycznie zobaczyć dokładnie, co to zrobiło, patrząc w katalog work. Jeśli zrobię code work, a następnie hash i następnie command run, który jeśli pamiętasz, to faktyczny plik, który jest wykonywany dla tego zadania, możemy wejść i możemy poszukać funkcji o nazwie NXF launch. I tutaj możesz zobaczyć dokładne polecenie docker, którego użył Nextflow, które wygląda podobnie do tego, co robiliśmy ręcznie w terminalu wcześniej. Docker run. Wiązanie tego katalogu hosta w kontenerze, a następnie określanie URL kontenera.

Więc nie ma tu magii. Po prostu Nextflow automatycznie robi ciężką pracę za Ciebie w sposób, który oznacza, że możesz łatwo określić kontenery w Swoim pipeline'ie, które są następnie łatwo dostępne dla każdego innego, kto uruchamia Twój workflow. I ci ludzie nie muszą już myśleć o zarządzaniu oprogramowaniem, aby uruchomić Twój pipeline analizy.

Bardzo, bardzo proste, bardzo wygodne, a także naprawdę powtarzalne. Dobre ze wszystkich stron.

Dobra, świetna robota. To koniec rozdziału 5. Dołącz do nas w następnym filmie w części 6, która jest ostatnią częścią tego szkolenia Hello Nextflow, gdzie porozmawiamy o konfiguracji Nextflow bardziej szczegółowo.

Do zobaczenia w następnym filmie.

[Następna transkrypcja wideo :octicons-arrow-right-24:](06_hello_config.md)
