# Część 4: Hello Modules - Transkrypcja

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Tłumaczenie wspomagane przez AI - [dowiedz się więcej i zasugeruj ulepszenia](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

<div class="video-wrapper">
  <iframe width="560" height="315" src="https://www.youtube.com/embed/Xxp_menS0E8?si=0AWnXB7xqHAzJdJV&amp;list=PLPZ8WHdZGxmXiHf8B26oB_fTfoKQdhlik" title="YouTube video player" frameborder="0" allow="accelerometer; autoplay; clipboard-write; encrypted-media; gyroscope; picture-in-picture; web-share" referrerpolicy="strict-origin-when-cross-origin" allowfullscreen></iframe>
</div>

!!!note "Ważne uwagi"

    Ta strona zawiera tylko transkrypcję. Aby uzyskać pełne instrukcje krok po kroku, wróć do [materiałów szkoleniowych](../04_hello_modules.md).

    Numery sekcji pokazane w transkrypcji są podane wyłącznie w celach poglądowych i mogą nie obejmować wszystkich numerów sekcji w materiałach.

## Powitanie

Cześć, witamy w czwartej części szkolenia Hello Nextflow.

Ten rozdział nazywa się Hello Modules i będziemy mówić o tym, jak modularyzować kod Nextflow. Zamierzamy wziąć nasz jeden skrypt i podzielić go na osobne pliki.

To sprawia, że kod jest łatwiejszy w nawigacji i utrzymaniu, gdy Twój workflow się powiększa, a także umożliwia udostępnianie modułów między pipeline'ami, dzięki czemu jeśli masz wiele pipeline'ów używających tego samego narzędzia, musisz napisać ten proces tylko raz.

Klasycznym przykładem tego jest repozytorium nf-core modules, które zawiera tysiące różnych gotowych do użycia narzędzi w modułach, które możesz zainstalować i użyć w Swoim workflow.

Nextflow może również pracować z sub workflows, które są jak moduły, ale mają wiele procesów. To wykracza poza zakres tego szkolenia, ale działa to w zasadzie w ten sam sposób.

Dobrze. Rzućmy okiem.

Jak zwykle, zacznij od przejścia do training.nextflow.io.

Przejdź do "Hello Nextflow" na pasku bocznym i robimy część czwartą: "Hello Modules".

Teraz przejdę do mojego środowiska GitHub Codespaces i przyjrzę się plikowi "hello-modules".

Tak jak wcześniej, zaczynamy od punktu końcowego poprzedniego rozdziału, więc ten skrypt powinien wyglądać znajomo. Mamy nasze trzy procesy, say hello, convert to upper i collect greetings, oraz prosty workflow, który uruchamia te trzy polecenia i emituje wiadomość na końcu. Mamy dwa parametry zwane greeting i batch, które określają nazwę, która jest używana dla zebranego pliku wyjściowego na końcu.

## 0. Rozgrzewka: Uruchom hello-modules.nf

Możemy zweryfikować, że ten workflow nadal działa zgodnie z oczekiwaniami, wykonując nextflow run hello-modules.

Świetnie. Uruchomiło trzy zadania z każdym z tych procesów, jedno zadanie collect i powiedziało nam, że w tej partii są trzy powitania. Jeśli wejdziemy do results, mamy nasze różne pliki wyjściowe tutaj, włącznie ze zebranym plikiem wyjściowym test batch.

## 1. Utwórz katalog do przechowywania modułów

Dobrze. Zróbmy trochę modularyzacji.

Generalnie dobrym pomysłem jest umieszczenie modułów w podfolderze w repozytorium Twojego pipeline'u, po prostu aby zachować porządek. Możesz nazwać to jak chcesz, ale zgodnie z konwencją zazwyczaj nazywamy to modules.

Więc chodźmy do terminala i zróbmy mkdir modules. Możesz zobaczyć, jak pojawia się w pasku bocznym VS Code tutaj.

## 2. Utwórz moduł dla sayHello()

Następnie utworzę nowy plik dla mojego pierwszego modułu. Możesz zrobić "touch" lub "code" lub możesz to zrobić na pasku bocznym, naprawdę to nie ma znaczenia. Więc zrobię code modules i nazwę go od nazwy procesu. Więc sayHello.nf. NF jest tradycyjnym rozszerzeniem pliku dla plików Nextflow.

Nacisnę zapisz tutaj i zobaczymy, jak pojawia się nasz nowy plik modułu.

## 2.2. Przenieś kod procesu sayHello do pliku modułu

Dobrze, teraz wezmę kod modułu z workflow. Wezmę również hash bang tutaj i najpierw go skopiuję, aby było jasne, że to plik Nextflow. A potem wezmę ten proces i wytnę. Więc usunę go z mojego głównego skryptu workflow i wkleję do tego nowego modułu.

To cała zawartość, którą będzie zawierał ten plik modułu. Tylko jeden proces, żadnego workflow, żadnej logiki, tylko sam proces.

Teraz mogę zamknąć ten plik.

## 2.3. Dodaj deklarację import przed blokiem workflow

Teraz mojemu workflow brakuje pierwszego procesu, więc musimy go przywrócić poprzez import. Składnia tego jest bardzo podobna do innych języków programowania, więc może być znajoma. Robimy include w klamrowych nawiasach, nazwa procesu, say hello, a następnie from i ścieżka pliku modules, say hello, nf. Fantastycznie.

Kilka sztuczek tutaj. Rozszerzenie VS Code jest w tym mądre. Rozpoznaje tę ścieżkę pliku i możesz najechać na nią kursorem i kliknąć follow link. Albo jestem na Mac, mogę nacisnąć option i kliknąć, a otworzy ten plik. Więc możemy szybko do niego przeskoczyć.

Ta nazwa procesu jest teraz używana przez workflow poniżej i możemy zrobić to samo tutaj. Pokazuje nam trochę informacji o tym procesie i znowu, mogę przytrzymać option, kliknąć na to i otworzy się w edytorze.

Więc to naprawdę szybki sposób, gdy masz wiele plików dla Swoich różnych procesów, aby szybko nawigować po bazie kodu w VS Code.

Dobrze. To w zasadzie wszystko w tym rozdziale. Teraz po prostu robimy to samo dla innych procesów.

## 3. Modularyzuj proces convertToUpper()

Więc stwórzmy nowy plik tutaj. Nazwijmy go Convert to upper nf. Znowu skopiujmy hash bang. A potem wytnijmy proces.

Skopiuj nazwę procesu tam, dodaj nową instrukcję include z nową nazwą procesu.

## 4. Modularyzuj proces collectGreetings()

A potem zróbmy to samo dla trzeciego procesu. Nowy plik, collect greetings,

zróbmy hash bang. Wytnijmy proces, wklejmy proces i zróbmy nową instrukcję include.

Teraz możesz zobaczyć tutaj, że mam podkreślenie błędu mówiące invalid include source. I to jest w rzeczywistości prawdziwy błąd, który popełniłem, ponieważ poruszałem się trochę za szybko. Jeśli spojrzysz uważnie, zobaczysz, że pominąłem T w convert to upper

Więc VS Code bardzo pomocnie powiedział mi, że popełniłem tam błąd. Jeśli naprawię tę nazwę pliku, błąd znika. To dobry przykład, dlaczego sprawdzanie błędów w VS Code jest tak przydatne do pisania kodu Nextflow. Inaczej bym tego nie zauważył i dowiedziałbym się o tym znacznie później, gdy próbowałbym uruchomić workflow.

Nasz główny skrypt pipeline wygląda teraz znacznie prościej. Nie ma w nim żadnych procesów, mamy tylko trzy instrukcje include i nasz workflow. Nie zmieniliśmy żadnej logiki workflow. Nie zmieniliśmy żadnego kodu procesu, więc miejmy nadzieję, że powinno działać dokładnie w ten sam sposób.

## 4.4. Uruchom workflow, aby zweryfikować, że robi to samo co wcześniej

Sprawdźmy. Otworzę terminal i uruchomię dokładnie to samo polecenie co wcześniej.

I rzeczywiście, uruchomiło nasze procesy, say hello, convert to upper, collect greetings i dało nam znowu trzy powitania.

Więc przenieśliśmy nasz kod, ale nie zmieniliśmy niczego w sposobie wykonywania workflow i jest całkowicie niezmieniony. Jedyna różnica polega na tym, że mamy teraz czystszy kod, łatwiejszy w utrzymaniu i łatwiejszy do udostępnienia innym.

I to wszystko. To był krótki rozdział. To prosta koncepcja, ale bardzo potężna i kluczowa dla sposobu, w jaki piszemy bardziej złożone workflow Nextflow. Więc ważne jest, abyś to rozumiał i wyrobił sobie nawyk korzystania z tego.

W następnym rozdziale zmienimy trochę tempo i przestaniemy myśleć tak dużo o składni pisania kodu Nextflow, a pomyślimy trochę o tym, jak używamy oprogramowania w samych procesach. Dołącz do nas w części piątej Hello Containers.

[Następna transkrypcja wideo :octicons-arrow-right-24:](05_hello_containers.md)
