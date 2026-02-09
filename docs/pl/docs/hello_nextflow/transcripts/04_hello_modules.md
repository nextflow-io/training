# Część 4: Hello Modules - Transkrypcja wideo

<div class="video-wrapper">
  <iframe width="560" height="315" src="https://www.youtube.com/embed/43Ot-f0iOME?si=0AWnXB7xqHAzJdJV&amp;list=PLPZ8WHdZGxmWKozQuzr27jyMGqp9kElVK&amp;cc_load_policy=1" title="YouTube video player" frameborder="0" allow="accelerometer; autoplay; clipboard-write; encrypted-media; gyroscope; picture-in-picture; web-share" referrerpolicy="strict-origin-when-cross-origin" allowfullscreen></iframe>
</div>

!!!note "Ważne uwagi"

    Ta strona zawiera tylko transkrypcję. Pełne instrukcje krok po kroku znajdziesz w [materiałach kursu](../04_hello_modules.md).

    Numery sekcji pokazane w transkrypcji służą wyłącznie celom orientacyjnym i mogą nie obejmować wszystkich numerów sekcji w materiałach.

## Powitanie

Cześć i witaj ponownie w czwartej części Hello Nextflow. Ta sekcja dotyczy modułów i jest dość krótką częścią kursu. Tak naprawdę nie będziemy pisać zbyt wiele kodu – chodzi bardziej o to, jak organizujemy kod w naszym pipeline'ie.

Do tej pory umieszczaliśmy wszystko w jednym pliku, co jest w porządku i tak właśnie budowaliśmy pipeline'y Nextflow'a w dawnych czasach.

Ale w miarę jak pipeline rośnie, skrypt staje się coraz dłuższy i trudniejszy w nawigacji, utrzymaniu, a także oznacza, że nie możemy tak naprawdę dzielić się żadnym kodem.

Moduły Nextflow'a pozwalają nam wyciągnąć procesy z głównego skryptu, a następnie je zaimportować. Oznacza to, że kod jest łatwiejszy w nawigacji, a także możemy dzielić się kodem modułów między różnymi pipeline'ami.

Ten mały diagram na głównej stronie dokumentacji ładnie pokazuje tę koncepcję. Zamiast jednego ogromnego skryptu, będziemy dołączać te oddzielne pliki modułów z różnych skryptów modułów, a wszystko zostanie wciągnięte do workflow'a, ale nadal będzie działać dokładnie w ten sam sposób.

Przejdźmy więc do GitHub Codespaces i rozejrzyjmy się trochę. Jak poprzednio, trochę posprzątałem moją przestrzeń roboczą. Usunąłem stare katalogi Nextflow'a i katalog work i tak dalej. Ale nie ma znaczenia, jeśli nadal masz te pliki.

Zacznę pracować w pliku hello modules, który jest w zasadzie tam, gdzie zostawiliśmy go na końcu poprzedniego rozdziału. Mamy tutaj nasze trzy procesy. Mamy kilka parametrów, blok workflow, w którym uruchamiamy te trzy procesy i łączymy je kanałami. Następnie publikujemy kanały wyjściowe i mamy blok output mówiący, jak publikować te pliki.

## 1. Utwórz katalog do przechowywania modułów

Teraz, jak mówię, tak naprawdę nie będziemy pisać ani edytować zbyt wiele kodu. Po prostu przeniesiemy kod, który już mamy. Pliki modułów Nextflow'a zazwyczaj zawierają pojedynczy proces, a zgodnie z konwencją normalnie przechowujemy je w katalogu o nazwie modules. Ale możesz nazwać go jak chcesz. Ja jednak zachowam katalog modules w moim repozytorium tutaj, a następnie utworzę jeden plik dla każdego procesu. Więc powiem new file, sayHello.nf.

## 2. Utwórz moduł dla sayHello()

Teraz wezmę mój proces i po prostu zaznaczę ten kod, wytnę go z głównego pliku hello modules i wkleję tutaj.

Oczywiście samo w sobie to nic nie robi. Nasz główny skrypt nadal potrzebuje tego procesu, więc musimy go jakoś z powrotem wciągnąć. I robimy to za pomocą instrukcji include.

Więc wpisuję include i klamrowe nawiasy, a następnie biorę nazwę procesu. I mówię from, a następnie podaję względną ścieżkę do pliku. Więc mówi, zaczyna się od ./, ponieważ jest względna od miejsca, gdzie zapisany jest ten skrypt. Więc to modules sayHello.nf.

Zauważ, że rozszerzenie VS Code jest tutaj dość pomocne. Mówi nam, czy może znaleźć ten plik i czy może znaleźć proces, który nazywam. Jeśli celowo wprowadzę tutaj literówkę, od razu dostanę błąd i powie mi, że nie może znaleźć tego procesu, który próbuję zaimportować. Więc po prostu zwracaj uwagę na wszelkie błędy, które znajdziesz.

I to właściwie wszystko. Nadal mamy tutaj nasz proces. Nie są potrzebne żadne zmiany tutaj na dole. Proces ma tę samą nazwę i jest wykonywany dokładnie w ten sam sposób. Po prostu faktyczny kod procesu jest teraz w osobnym pliku.

Możemy ponownie uruchomić workflow Nextflow'a, będzie działać dokładnie w ten sam sposób. I to w zasadzie reszta tego rozdziału kursu to po prostu przenoszenie tych trzech procesów do ich własnych plików.

Zróbmy to teraz. Szybko utworzę nowy plik modułu dla drugiego procesu: convertToUpper.nf. Wytnę ten kod, wkleję go tam. A następnie dołączę ten. Świetnie.

A potem utworzę nowy plik dla collectGreetings.nf. Wytnę to.

Dużo wycinania, kopiowania i wklejania.

I teraz nasz główny skrypt workflow'a nagle wygląda znacznie, znacznie krócej, jest bardziej przystępny i dużo łatwiejszy do odczytania.

I widzisz, jak projekt zaczyna się teraz budować z naszymi różnymi plikami. Możemy zagłębić się w szczegóły w miejscach, które nas interesują. Poruszać się, aby znaleźć konkretne kroki w pipeline'ie znacznie łatwiej i szybko uzyskać przegląd tego, co robi pipeline.

## Nawigacja po modułach za pomocą VS Code

Oczywiście wadą tego jest to, że jeśli masz duży pipeline, będziesz mieć wiele plików modułów i mogą być zorganizowane w wielu podkatalogach lub różnych rzeczach. Teraz, znowu, jedna mała wskazówka tutaj. Rozszerzenie VS Code jest całkiem dobre w nawigowaniu po Twojej bazie kodu i również informowaniu Cię o kodzie tam.

Widzisz, VS Code rozumie, czym jest ten proces i daje mi mały przegląd, gdy najadę kursorem, więc mogę zobaczyć bez konieczności szukania kodu źródłowego, jakie są wejścia i wyjścia, co jest zazwyczaj najważniejszą rzeczą, gdy używam go w workflow'ie.

A także, jeśli przytrzymam command, jestem na Macu, i kliknę nazwę procesu, otwiera plik bezpośrednio od razu. Wciąga go. Więc mogę po prostu tam przeskoczyć bez nawet myślenia o tym, jakie są faktyczne ścieżki plików. I to działa wszędzie, mogę to zrobić również tam, gdzie procesy są wywoływane. Więc to naprawdę przyspiesza.

## 4.4. Uruchom workflow

Dobra, sprawdźmy tylko, czy pipeline nadal działa tak, jak się spodziewamy. Więc wywołam terminal. Zróbmy "nextflow run hello modules" i zobaczmy, czy wykonuje się bez żadnych problemów.

Mam nadzieję, że cały sens tego jest taki, że pipeline jest w zasadzie niezmieniony, więc tak naprawdę nie powinieneś zobaczyć żadnych zmian w porównaniu do tego, kiedy uruchamialiśmy go wcześniej. Wyjście tutaj wygląda dokładnie tak samo i widzisz nasz katalog results ze wszystkimi tymi samymi plikami, więc to świetnie. Brak zmian to dobrze.

## Uwaga o nf-core/modules

Zanim zakończymy, chcę szybko poruszyć temat mocy współpracy, jeśli chodzi o moduły. Te pliki znajdują się w moim repozytorium, więc nie jest od razu oczywiste, jak moglibyśmy nad nimi współpracować. I jest wiele różnych sposobów, w jakie możesz to zrobić, ale prawdopodobnie największym i najbardziej znanym przykładem tego jest nf-core.

Jeśli przejdę na stronę nf-core, przejdę do resources, a następnie modules. Widzisz, że nf-core ma ogromną bibliotekę modułów, prawie niecałe 1700 modułów, kiedy to oglądam. I mogę wpisać nazwę dowolnego z moich ulubionych narzędzi, znaleźć, czy ktoś inny już napisał dla niego moduł i zobaczyć ten wcześniej napisany proces modułu tutaj ze wszystkimi wejściami, wyjściami, kontenerami oprogramowania, wszystkimi tymi informacjami, a po stronie możesz zobaczyć, ile różnych pipeline'ów nf-core używa tego pojedynczego współdzielonego procesu.

To jest trochę ekstremalny przykład, ale widzisz, że to naprawdę ponowne wykorzystanie tego kodu. A jeśli kliknę przez do źródła GitHub dla tego, to dokładnie to samo, co robimy. To po prostu duży proces w pliku.

Teraz po stronie nf-core robimy pewne sztuczki, aby móc dzielić się tymi plikami i wprowadzać je do różnych repozytoriów. A jeśli chcesz wiedzieć więcej na ten temat, sprawdź kurs, który mamy o używaniu i budowaniu z nf-core konkretnie. Ale chciałem tylko dać Ci wyobrażenie o tym, jak potężna może być ta koncepcja ponownego wykorzystania kodu.

## Podsumowanie

Dobra, to wszystko o modułach. Mówiłem Ci, że to krótka sekcja kursu. Sprawdź quiz, upewnij się, że rozumiesz i upewnij się, że wszystko nadal działa poprawnie. I zobaczę Cię z powrotem w następnym wideo, które dotyczy kontenerów oprogramowania. Dziękuję bardzo.
