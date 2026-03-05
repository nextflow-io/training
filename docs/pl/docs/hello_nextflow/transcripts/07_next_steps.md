# Kolejne kroki - transkrypcja wideo

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Tłumaczenie wspomagane przez AI - [dowiedz się więcej i zasugeruj ulepszenia](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

<div class="video-wrapper">
  <iframe width="560" height="315" src="https://www.youtube.com/embed/wnmUAfDL01E?si=Lp8hS8RdaMwbp5j5&amp;list=PLPZ8WHdZGxmWKozQuzr27jyMGqp9kElVK&amp;cc_load_policy=1" title="YouTube video player" frameborder="0" allow="accelerometer; autoplay; clipboard-write; encrypted-media; gyroscope; picture-in-picture; web-share" referrerpolicy="strict-origin-when-cross-origin" allowfullscreen></iframe>
</div>

!!!note "Ważne uwagi"

    Ta strona zawiera tylko transkrypcję. Pełne instrukcje krok po kroku znajdziesz w [materiale szkoleniowym](../next_steps.md).

## Powitanie

​

Gratulacje, udało Ci się!

Dotarłeś do końca i ukończyłeś szkolenie Hello Nextflow. Mamy nadzieję, że się podobało. Dziękujemy, że wytrwałeś z nami do samego końca. Naprawdę doceniamy czas i wysiłek, który włożyłeś w naukę Nextflow'a. Mamy nadzieję, że będzie przydatny w Twojej pracy.

## Inne szkolenia na training.nextflow.io

Nie zapomnij regularnie wracać na training.nextflow.io. Cały czas dodajemy nowe krótkie szkolenia, a także odświeżamy wiele materiałów, które już tam są. To szkolenie Hello Nextflow będzie aktualizowane na przestrzeni czasu.

Jest to szczególnie ważne, ponieważ aktualizujemy składnię w Nextflow'ie, a rok 2026 przyniesie sporo nowych funkcji, więc to szkolenie będzie wyglądać i działać nieco inaczej, gdy następnym razem przeprowadzimy je w 2027 roku.

W szczególności chcę zwrócić uwagę na stronę "Nextflow for Science". To krótkie szkolenia zaprojektowane jako kontynuacja tego szkolenia Hello Nextflow. Pokazują, jak używać Nextflow'a w różnych konkretnych przypadkach użycia, czy to genomika, RNAseq, czy wiele innych rzeczy. Staramy się cały czas dodawać więcej naukowych przypadków użycia.

Są też Side Quests. Kiedy tworzymy szkolenie takie jak Hello Nextflow, jest tak wiele rzeczy, które moglibyśmy omówić, i trudno jest utrzymać wszystko w zakresie. Więc jeśli jest konkretny temat, który uważamy za interesujący dla ludzi i zasługujący na głębsze omówienie, umieszczamy go w Side Quest.

Przejrzyj je i jeśli są różne rzeczy, które mogą być istotne dla Twojej pracy, jak nf-test czy różne operacje na metadanych i typowe wzorce skryptowe, sprawdź Side Quests i zobacz, czy warto dowiedzieć się więcej.

Jest też szkolenie o nf-core. Mamy nadzieję, że znasz już ten projekt, ale jeśli nie, sprawdź go. Jest tam prawie 150 różnych pipeline'ów do różnych typów analiz i różnych typów danych, więc całkiem możliwe, że będzie gotowy pipeline do użycia od razu dla rodzaju analizy danych, której potrzebujesz.

Co ważne, w nf-core są też komponenty, prawie 1700 różnych modułów, różnych procesów i wrapperów dla narzędzi. A dzięki narzędziom, które są w nf-core, możesz je mieszać i dopasowywać oraz budować własny pipeline jak z klocków Lego. Znacznie szybciej i bardziej powtarzalnie.

## Seqera Platform

W miarę jak zwiększasz swoje użycie Nextflow'a, sprawdź Seqera Platform - to najlepszy sposób na uruchamianie Nextflow'a. Możesz uruchamiać na własnej infrastrukturze, więc HPC lub AWS, Azure, Google Cloud, Oracle i więcej. Możesz też użyć naszego własnego Seqera Compute, jeśli w ogóle nie chcesz zarządzać infrastrukturą obliczeniową.

Seqera Platform naprawdę upraszcza konfigurację tych złożonych infrastruktur chmurowych dzięki funkcjom takim jak Batch Forge, który tworzy środowisko za Ciebie. Pomaga też naprawdę z obserwowalnością, logowaniem audytu i zgodnością.

Obiektywnie sprawia, że pipeline'y działają taniej i szybciej dzięki technologiom takim jak Fusion, które optymalizują dostęp do dysku i transfery danych. Jest też optymalizacja pipeline'ów, aby upewnić się, że konfiguracja Twoich pipeline'ów jest jak najlepiej dostrojona.

Są też zupełnie inne funkcje poza uruchamianiem pipeline'ów. Mamy Studios, gdzie możesz uruchamiać interaktywne analizy i tworzyć środowiska z dowolnego niestandardowego obrazu dockera, który stworzysz. I Data Explorer, który pomaga eksplorować różne systemy plików, gdziekolwiek się znajdują.

Jest darmowy tier Seqera Platform, więc możesz używać prawie wszystkich tych funkcji za darmo już teraz. A nawet damy Ci sto dolarów darmowego kredytu obliczeniowego z Seqera Compute, jeśli zarejestrujesz się swoim organizacyjnym adresem e-mail. Wreszcie, jest program akademicki, więc jeśli pracujesz na uniwersytecie, sprawdź stronę z cenami, znajdź tam formularz i daj nam znać, a zaktualizujemy Cię do Cloud Pro za darmo.

## Pomoc społeczności i wydarzenia

Dobra. Idąc dalej. Jeśli kiedykolwiek będziesz potrzebować wsparcia z Nextflow'em, sprawdź community.seqera.io. Jest naprawdę aktywne i mamy nadzieję zobaczyć Cię tam, aby omówić Twoje różne problemy i przypadki użycia, a może teraz będziesz mógł nawet pomóc innym ludziom.

Mamy też wiele wydarzeń. Mamy wydarzenia społeczności pochodzące z nf-core i Nextflow. Mamy online'owy i rozproszony hackathon nf-core w marcu - w zeszłym roku dołączyło ponad tysiąc osób z miejsc na całym świecie. Więc dołącz do nas, jeśli możesz.

Mamy też wydarzenia Nextflow Summit, jedno w Bostonie, a potem wydarzenie w Barcelonie i online. Fantastyczne prezentacje, gdzie możesz usłyszeć o ludziach używających Nextflow'a w naprawdę masowych, dzikich i ekscytujących różnych sposobach. Są też hackathony związane z nimi i szkolenia stacjonarne.

## Podcast i blog Nextflow

Jeśli chcesz być na bieżąco z tym, co dzieje się w ekosystemie Nextflow, sprawdź seqera.io/blog.

Jest tam sekcja dla Nextflow'a, gdzie możesz znaleźć posty blogowe społeczności od ludzi pracujących w społeczności, a także posty blogowe od Seqera o aktualizacjach Nextflow'a i innych narzędzi, które tworzymy.

Chciałbym też polecić mój ulubiony projekt, którym jest Nextflow Podcast. Sprawdź go na Spotify, Apple Music lub YouTube. Publikujemy nowe odcinki okresowo, gdzie rozmawiam z innymi ludźmi, albo pracującymi z Nextflow'em lub powiązanymi technologiami, albo ludźmi ze społeczności. Robimy naprawdę techniczne, głębokie zanurzenia w to, jak rzeczy działają i co ludzie robią. Więc jeśli jesteś zainteresowany, sprawdź je. Są naprawdę fajne.

## Podziękowania

Dobra, chciałbym złożyć podziękowania. Zespół szkoleniowy w Seqera jest odpowiedzialny za ten materiał. Siedzę przed kamerą, ale naprawdę całą ciężką pracę wykonali ci inni ludzie. Szczególne podziękowania dla Geraldine, która napisała i odświeżyła materiał szkoleniowy dla Hello Nextflow i innych. A także dla Jona, który naprawdę pomógł, szczególnie z aktualizacją składni dla nowej składni Nextflow'a, a także pisząc wiele szkoleń sam. Inni w zespole rozwoju naukowego, tacy jak Rike, Rob, Florian i wielu innych, mieli ogromny wkład w materiał, nad którym pracowaliśmy.

Chciałbym też podziękować ludziom ze społeczności. Nowe tłumaczenia, na przykład, które są bardzo świeże, były mocno wspierane przez ludzi z programu ambasadorów i innych miejsc. I naprawdę, otwartoźródłowy charakter materiału szkoleniowego oznacza, że mamy pull requesty i issues przychodzące dość często, co naprawdę nam pomaga.

## Ankieta

Teraz, gdy skończyłeś, jeśli jeszcze tego nie zrobiłeś, proszę, szybko wypełnij ankietę zwrotną. Jest na stronie training.nextflow.io tuż pod sekcją Hello Nextflow.

To tylko pięć pytań. Jest naprawdę, naprawdę szybka, ale pozwala nam śledzić mniej więcej, ile osób robi szkolenie, a także możesz nam powiedzieć, jak ulepszyć materiał szkoleniowy. Naprawdę sprawdzamy wszystkie odpowiedzi, więc naprawdę cenimy Twój wkład tam.

## Pożegnanie

Jeszcze raz, wielkie dzięki za dołączenie do nas na tym szkoleniu i na tej drodze. Zgłoś issue lub Pull Request na GitHubie, jeśli zauważyłeś cokolwiek w materiale szkoleniowym, co Twoim zdaniem można by poprawić. I naprawdę mam nadzieję zobaczyć Cię na innym szkoleniu Nextflow, lub może na hackathonie lub wydarzeniu. Jeszcze raz dziękuję.​
