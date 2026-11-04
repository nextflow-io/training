# Część 4: Hello Modules - transkrypcja wideo

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Tłumaczenie wspomagane przez sztuczną inteligencję - [dowiedz się więcej i zaproponuj ulepszenia](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

<div class="video-wrapper">
  <iframe width="560" height="315" src="https://www.youtube.com/embed/43Ot-f0iOME?si=0AWnXB7xqHAzJdJV&amp;list=PLPZ8WHdZGxmWKozQuzr27jyMGqp9kElVK&amp;cc_load_policy=1" title="YouTube video player" frameborder="0" allow="accelerometer; autoplay; clipboard-write; encrypted-media; gyroscope; picture-in-picture; web-share" referrerpolicy="strict-origin-when-cross-origin" allowfullscreen></iframe>
</div>

!!!note "Ważne uwagi"

    Ta strona zawiera tylko transkrypcję. Pełne instrukcje krok po kroku znajdziesz w [materiałach szkoleniowych](../04_hello_modules.md).

    Numery sekcji pokazane w transkrypcji służą jedynie celom orientacyjnym i mogą nie obejmować wszystkich numerów sekcji zawartych w materiałach.

## Powitanie

Cześć i witaj z powrotem w części czwartej Hello Nextflow. Ta sekcja dotyczy wyłącznie modułów i jest dość krótką częścią szkolenia. Właściwie nie napiszemy tu zbyt wiele kodu - chodzi bardziej o to, jak organizujemy kod w naszym pipeline'ie.

Do tej pory umieszczaliśmy wszystko w jednym pliku, co jest w porządku i tak właśnie budowaliśmy pipeline'y Nextflow'a w dawnych czasach.

Jednak w miarę wzrostu pipeline'u skrypt staje się coraz dłuższy, coraz trudniejszy w nawigacji i utrzymaniu, a także oznacza, że nie możemy naprawdę dzielić się żadnym z tego kodu.

Moduły Nextflow'a pozwalają nam wydzielić procesy z głównego skryptu, a następnie je zaimportować. Oznacza to, że kod jest łatwiejszy w nawigacji, a także możemy współdzielić kod modułów między różnymi pipeline'ami.

Ten mały diagram na głównej stronie dokumentacji ładnie ilustruje tę koncepcję. Zamiast jednego ogromnego skryptu będziemy dołączać te oddzielne pliki modułów z różnych skryptów modułowych, a wszystko zostanie włączone do workflow'a, ale nadal będzie działać dokładnie w ten sam sposób.

Przejdźmy więc do GitHub Codespaces i przyjrzyjmy się temu bliżej. Jak poprzednio, trochę posprzątałem moje środowisko pracy. Usunąłem stare katalogi Nextflow'a, katalog work i tak dalej. Ale nie ma problemu, jeśli nadal masz te pliki.

Zacznę pracę z pliku hello modules, który zasadniczo znajduje się w stanie, w jakim zostawiliśmy go na końcu poprzedniego rozdziału. Mamy tutaj nasze trzy procesy. Mamy kilka parametrów, blok workflow, w którym uruchamiamy te trzy procesy i łączymy je ze sobą za pomocą kanałów. Następnie publikujemy kanały wyjściowe i mamy blok output określający sposób publikacji tych plików.

## 1. Utwórz katalog do przechowywania modułów

Jak wspomniałem, tak naprawdę nie będziemy pisać ani edytować zbyt wiele kodu. Po prostu przeniesiemy kod, który już mamy. Pliki modułów Nextflow'a zazwyczaj zawierają pojedynczy proces, a zgodnie z konwencją normalnie przechowujemy je w katalogu o nazwie modules. Możesz nazwać go jednak jak chcesz. Ale ja utworzę katalog modules w moim repozytorium tutaj, a następnie utworzę jeden plik dla każdego procesu. Powiem więc new file, sayHello.nf.

## 2. Utwórz moduł dla sayHello()

Teraz wezmę mój proces i po prostu wybiorę ten kod, wytnę go z głównego pliku hello modules i wkleję tutaj.

Oczywiście samo to nic nie robi. Nasz główny skrypt nadal potrzebuje tego procesu, więc musimy go jakoś z powrotem włączyć. A robimy to za pomocą instrukcji include.

Wpisuję więc include i klamry, a następnie podaję nazwę procesu. I mówię from, a następnie podaję względną ścieżkę do pliku. Więc mówi, zaczyna się od ./, ponieważ jest względna w stosunku do miejsca, w którym zapisany jest ten skrypt. Czyli to modules sayHello.nf.

Zauważ, że rozszerzenie VS Code jest tutaj dość pomocne. Informuje nas, czy może znaleźć ten plik i czy może znaleźć proces, który wymieniam. Jeśli celowo wprowadzę tutaj literówkę, od razu pojawia się błąd i poinformuje mnie, że nie może znaleźć tego procesu, który próbuję zaimportować. Więc zwracaj uwagę na wszelkie błędy, które znajdziesz.

I to właściwie wszystko. Nadal mamy tutaj nasz proces. Nie są potrzebne żadne zmiany poniżej. Proces ma tę samą nazwę i jest wykonywany dokładnie w ten sam sposób. Po prostu rzeczywisty kod procesu znajduje się teraz w osobnym pliku.

Możemy ponownie uruchomić workflow Nextflow'a, będzie działać dokładnie w ten sam sposób. I to w zasadzie reszta tego rozdziału szkolenia - po prostu przenoszenie tych trzech procesów do ich własnych plików.

Zróbmy to teraz. Szybko utworzę nowy plik modułu dla drugiego procesu: convertToUpper.nf. Wytnę ten kod, wkleję go tam. A następnie zaimportuję ten. Świetnie.

A potem utworzę nowy plik dla collectGreetings.nf. Wytnę to.

Dużo wycinania i kopiowania oraz wklejania.

I teraz nasz główny skrypt workflow wygląda nagle znacznie, znacznie krócej, jest znacznie bardziej przystępny i o wiele łatwiejszy do przeczytania.

Widzisz, jak projekt zaczyna się teraz budować z naszymi różnymi plikami. Możemy zagłębić się w szczegóły tam, gdzie chcemy. Poruszać się w poszukiwaniu konkretnych kroków w pipeline'ie znacznie łatwiej i szybko uzyskać przegląd tego, co robi pipeline.

## Nawigacja po modułach za pomocą VS Code

Oczywiście wadą tego podejścia jest to, że jeśli masz duży pipeline, będziesz mieć wiele plików modułów i mogą być zorganizowane w wielu podkatalogach lub wszelkiego rodzaju rzeczach. Tutaj znowu mała wskazówka. Rozszerzenie VS Code jest całkiem dobre w nawigowaniu po Twojej bazie kodu i informowaniu Cię o kodzie tam znajdującym się.

Widzisz, że VS Code rozumie, czym jest ten proces i daje mi mały podgląd, kiedy najadę kursorem, więc mogę zobaczyć bez konieczności szukania kodu źródłowego, jakie są wejścia i wyjścia, co zazwyczaj jest najważniejszą rzeczą, gdy używam go w workflow.

A także jeśli przytrzymam command, jestem na Macu, i kliknę nazwę procesu, otworzy plik bezpośrednio od razu. Wciągnie go. Więc mogę tam po prostu skoczyć bez myślenia o tym, jakie są rzeczywiste ścieżki plików. I to działa wszędzie, mogę to zrobić również tam, gdzie procesy są wywoływane. To naprawdę szybkie.

## 4.4. Uruchom workflow

Dobra, po prostu sprawdźmy, czy pipeline nadal działa zgodnie z oczekiwaniami. Więc wywołam terminal. Zróbmy "nextflow run hello modules" i zobaczmy, czy wykona się bez problemów.

Mam nadzieję, że cały sens tego jest taki, że pipeline jest zasadniczo niezmieniony, więc tak naprawdę nie powinieneś zobaczyć żadnych zmian w porównaniu z tym, kiedy uruchamialiśmy go wcześniej. Wyjście tutaj wygląda dokładnie tak samo i widzisz nasz katalog results ze wszystkimi tymi samymi plikami, więc to świetnie. Brak zmian to dobrze.

## Uwaga o nf-core/modules

Zanim zakończymy, chcę krótko wspomnieć o sile współpracy, jeśli chodzi o moduły. Te pliki znajdują się w moim repozytorium, więc nie jest od razu oczywiste, jak moglibyśmy nad nimi współpracować. I jest wiele różnych sposobów, w które możesz to zrobić, ale prawdopodobnie największym i najbardziej znanym przykładem jest nf-core.

Jeśli przejdę na stronę nf-core, przejdę do resources, a następnie modules. Widzisz, że nf-core ma ogromną bibliotekę modułów, prawie niecałe 1700 modułów, kiedy to oglądam. Mogę więc wpisać nazwę dowolnego z moich ulubionych narzędzi, poszukać, czy ktoś inny już napisał dla niego moduł i zobaczyć ten wcześniej napisany proces modułu tutaj ze wszystkimi wejściami, wyjściami, kontenerami oprogramowania, wszystkimi tymi informacjami, a po boku widzisz, ile różnych pipeline'ów nf-core używa tego pojedynczego współdzielonego procesu.

To trochę ekstremalny przykład, ale widzisz, że to naprawdę ponowne wykorzystanie tego kodu. A jeśli kliknę, aby przejść do źródła GitHub dla tego, to dokładnie to samo, co robimy. To po prostu duży proces w pliku.

Po stronie nf-core robimy pewne sztuczki, aby móc udostępniać te pliki i wprowadzać je do różnych repozytoriów. A jeśli chcesz dowiedzieć się więcej na ten temat, sprawdź szkolenie, które mamy o używaniu i budowaniu z nf-core konkretnie. Ale chciałem tylko dać Ci wyobrażenie o tym, jak potężna może być ta koncepcja ponownego wykorzystania kodu.

## Podsumowanie

W porządku, to wszystko jeśli chodzi o moduły. Mówiłem Ci, że to krótka sekcja szkolenia. Sprawdź quiz, upewnij się, że rozumiesz to i upewnij się, że wszystko nadal działa prawidłowo. Zobaczę się w następnym filmie, który dotyczy kontenerów oprogramowania. Dziękuję bardzo.

I.
