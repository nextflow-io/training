# Część 4: Hello Modules - Transkrypcja

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Tłumaczenie wspomagane przez AI - [dowiedz się więcej i zasugeruj ulepszenia](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

<div class="video-wrapper">
  <iframe width="560" height="315" src="https://www.youtube.com/embed/Xxp_menS0E8?si=0AWnXB7xqHAzJdJV&amp;list=PLPZ8WHdZGxmXiHf8B26oB_fTfoKQdhlik" title="YouTube video player" frameborder="0" allow="accelerometer; autoplay; clipboard-write; encrypted-media; gyroscope; picture-in-picture; web-share" referrerpolicy="strict-origin-when-cross-origin" allowfullscreen></iframe>
</div>

!!!note "Ważne uwagi"

    Ta strona zawiera tylko transkrypcję. Aby uzyskać pełne instrukcje krok po kroku, wróć do [materiałów szkoleniowych](../04_hello_modules.md).

    Numery sekcji widoczne w transkrypcji są podane wyłącznie w celach poglądowych i mogą nie obejmować wszystkich numerów sekcji w materiałach.

## Powitanie

Cześć, witamy w czwartej części szkolenia Hello Nextflow.

Ten rozdział nazywa się Hello Modules i omówimy w nim sposób modularyzacji kodu Nextflow'a. Weźmiemy nasz workflow zapisany w jednym pliku i podzielimy go na osobne pliki.

Dzięki temu kod będzie łatwiejszy w nawigacji i utrzymaniu, gdy Twój workflow się rozrośnie. Umożliwia to również współdzielenie modułów między wieloma pipeline'ami, więc jeśli masz wiele pipeline'ów korzystających z tego samego narzędzia, musisz napisać ten proces tylko raz.

Klasycznym przykładem jest repozytorium nf-core modules, które zawiera tysiące różnych narzędzi w gotowych do użycia modułach, które możesz zainstalować i wykorzystać w swoim workflow.

Nextflow może również pracować z sub workflow'ami, które są jak moduły, ale zawierają wiele procesów. To wykracza poza zakres tego szkolenia, ale działają one w zasadzie w ten sam sposób.

Dobrze. Zobaczmy, jak to wygląda.

Jak zwykle zacznij od przejścia na training.nextflow.io.

Przejdź do "Hello Nextflow" na pasku bocznym i otwórz część czwartą: "Hello Modules".

Teraz przejdę do mojego środowiska GitHub Codespaces i przyjrzę się plikowi "hello-modules".

Tak jak wcześniej, zaczynamy od punktu końcowego poprzedniego rozdziału, więc ten skrypt powinien wyglądać znajomo. Mamy nasze trzy procesy: say hello, convert to upper i collect greetings, oraz prosty workflow, który uruchamia te trzy polecenia i emituje wiadomość na końcu. Mamy dwa parametry: greeting i batch, które określają nazwę zebranego pliku wyjściowego na końcu.

## 0. Rozgrzewka: Uruchom hello-modules.nf

Możemy zweryfikować, że ten workflow nadal działa zgodnie z oczekiwaniami, wykonując polecenie nextflow run hello-modules.

Świetnie. Uruchomiło trzy zadania dla każdego z tych procesów, jedno zadanie collect i powiedziało nam, że w tej partii są trzy powitania. Jeśli wejdziemy do katalogu results, mamy tam nasze różne pliki wyjściowe, włącznie ze zebranym plikiem test batch.

## 1. Utwórz katalog do przechowywania modułów

Dobrze. Zajmijmy się modularyzacją.

Generalnie dobrym pomysłem jest umieszczenie modułów w podkatalogu w repozytorium Twojego pipeline'u, aby zachować porządek. Możesz nazwać go jak chcesz, ale zgodnie z konwencją zazwyczaj nazywamy go modules.

Więc przejdźmy do terminala i wykonajmy mkdir modules. Możesz zobaczyć, jak pojawia się na pasku bocznym VS Code.

## 2. Utwórz moduł dla sayHello()

Następnie utworzę nowy plik dla mojego pierwszego modułu. Możesz użyć "touch" lub "code" albo zrobić to na pasku bocznym — naprawdę nie ma to znaczenia. Więc wykonam code modules i nazwę go od nazwy procesu. Tak więc sayHello.nf. NF jest tradycyjnym rozszerzeniem pliku dla plików Nextflow'a.

Nacisnę zapisz i zobaczymy, jak pojawia się nasz nowy plik modułu.

## 2.2. Przenieś kod procesu sayHello do pliku modułu

Dobrze, teraz wezmę kod procesu z workflow'a. Wezmę również hash bang i skopiuję go na początku, aby było jasne, że to plik Nextflow'a. A potem wezmę ten proces i wytnę go. Usunę go więc z mojego głównego skryptu workflow'a i wkleję do tego nowego modułu.

To cała zawartość tego pliku modułu. Tylko jeden proces, żadnego workflow'a, żadnej logiki — tylko sam proces.

Teraz mogę zamknąć ten plik.

## 2.3. Dodaj deklarację import przed blokiem workflow

Teraz mojemu workflow'owi brakuje tego pierwszego procesu, więc musimy go przywrócić poprzez import. Składnia jest bardzo podobna do innych języków programowania, więc może Ci się wydawać znajoma. Robimy include w klamrowych nawiasach, nazwa procesu say hello, a następnie from i ścieżka pliku modules, say hello, nf. Fantastycznie.

Kilka sztuczek. Rozszerzenie VS Code radzi sobie z tym inteligentnie. Rozpoznaje tę ścieżkę pliku i możesz najechać na nią kursorem i kliknąć follow link. Albo jestem na Mac, mogę nacisnąć option i kliknąć, a otworzy ten plik. Możemy więc szybko do niego przeskoczyć.

Ta nazwa procesu jest teraz używana przez workflow poniżej i możemy zrobić to samo tutaj. Pokazuje nam trochę informacji o tym procesie i znowu, mogę przytrzymać option, kliknąć na to i otworzy się w edytorze.

To naprawdę szybki sposób, gdy masz wiele plików dla różnych procesów, aby sprawnie nawigować po bazie kodu w VS Code.

Dobrze. To w zasadzie wszystko w tym rozdziale. Teraz po prostu robimy to samo dla pozostałych procesów.

## 3. Modularyzuj proces convertToUpper()

Więc stwórzmy tutaj nowy plik. Nazwijmy go Convert to upper nf. Znowu skopiujmy hash bang. A potem wytnijmy proces.

Skopiuj tam nazwę procesu i dodaj nową instrukcję include z nową nazwą procesu.

## 4. Modularyzuj proces collectGreetings()

A potem zróbmy to samo dla trzeciego procesu. Nowy plik, collect greetings,

dodajmy hash bang. Wytnijmy proces, wklejmy proces i dodajmy nową instrukcję include.

Teraz możesz zobaczyć tutaj, że mam podkreślenie błędu mówiące invalid include source. I to jest w rzeczywistości prawdziwy błąd, który popełniłem, ponieważ poruszałem się trochę za szybko. Jeśli spojrzysz uważnie, zobaczysz, że pominąłem T w convert to upper

Więc VS Code bardzo pomocnie powiedział mi, że popełniłem tam błąd. Jeśli naprawię tę nazwę pliku, błąd znika. To dobry przykład, dlaczego sprawdzanie błędów w VS Code jest tak przydatne do pisania kodu Nextflow'a. Inaczej bym tego nie zauważył i dowiedziałbym się o tym znacznie później, gdy próbowałbym uruchomić workflow'a.

Nasz główny skrypt pipeline'u wygląda teraz znacznie prościej. Nie ma w nim żadnych procesów, mamy tylko trzy instrukcje include i nasz workflow. Nie zmieniliśmy żadnej logiki workflow'a. Nie zmieniliśmy żadnego kodu procesu, więc miejmy nadzieję, że powinno działać dokładnie w ten sam sposób.

## 4.4. Uruchom workflow, aby zweryfikować, że robi to samo co wcześniej

Sprawdźmy. Otworzę terminal i uruchomię dokładnie to samo polecenie co wcześniej.

I rzeczywiście, uruchomiło nasze procesy: say hello, convert to upper, collect greetings i dało nam znowu trzy powitania.

Przenieśliśmy więc nasz kod, ale nie zmieniliśmy niczego w sposobie wykonywania workflow'a i pozostał on całkowicie niezmieniony. Jedyna różnica polega na tym, że mamy teraz czystszy kod, łatwiejszy w utrzymaniu i łatwiejszy do udostępniania innym.

I to wszystko. To był krótki rozdział. To prosta koncepcja, ale bardzo potężna i kluczowa dla sposobu, w jaki piszemy bardziej złożone workflow'e Nextflow'a. Więc ważne jest, abyś to zrozumiał i wyrobił sobie nawyk korzystania z tego.

W następnym rozdziale trochę zmienimy tempo i przestaniemy myśleć tak dużo o składni pisania kodu Nextflow'a, a pomyślimy o tym, jak używamy oprogramowania w samych procesach. Dołącz do nas w części piątej Hello Containers.

[Następna transkrypcja wideo :octicons-arrow-right-24:](05_hello_containers.md)
