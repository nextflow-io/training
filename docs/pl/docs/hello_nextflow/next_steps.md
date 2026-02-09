# Podsumowanie kursu

Gratulacje ukończenia kursu szkoleniowego Hello Nextflow! 🎉

<div class="video-wrapper">
  <iframe width="560" height="315" src="https://www.youtube.com/embed/wnmUAfDL01E?si=y8lAedhEHWaTV4zd&amp;list=PLPZ8WHdZGxmWKozQuzr27jyMGqp9kElVK&amp;cc_load_policy=1" title="YouTube video player" frameborder="0" allow="accelerometer; autoplay; clipboard-write; encrypted-media; gyroscope; picture-in-picture; web-share" referrerpolicy="strict-origin-when-cross-origin" allowfullscreen></iframe>
</div>

/// caption
:fontawesome-brands-youtube:{ .youtube } Zobacz [całą playlistę na kanale YouTube Nextflow](https://youtube.com/playlist?list=PLPZ8WHdZGxmWKozQuzr27jyMGqp9kElVK&si=eF7cLR62goy-lc6n).

:green_book: Możesz przeczytać [transkrypcję wideo](./transcripts/07_next_steps.md) równolegle z filmem.
///

## Twoja droga

Zacząłeś od bardzo prostego workflow'a, który uruchamiał zakodowane na stałe polecenie.
W ciągu sześciu części przekształciłeś ten podstawowy workflow w modularny, wieloetapowy pipeline, który wykorzystuje kluczowe funkcje Nextflow'a, w tym kanały, operatory, wbudowaną obsługę kontenerów oraz opcje konfiguracji.

### Co zbudowałeś

- Ostateczna forma workflow'a Hello przyjmuje jako wejście plik CSV zawierający tekstowe powitania.
- Cztery etapy są zaimplementowane jako procesy Nextflow (`sayHello`, `convertToUpper`, `collectGreetings` i `cowpy`) przechowywane w oddzielnych plikach modułów.
- Wyniki są publikowane do katalogu o nazwie `results/`.
- Końcowym wyjściem pipeline'u jest plik tekstowy zawierający grafikę ASCII postaci wypowiadającej powitania pisane wielkimi literami.

<figure class="excalidraw">
--8<-- "docs/en/docs/hello_nextflow/img/hello_pipeline_complete.svg"
</figure>

1. **`sayHello`:** Zapisuje każde powitanie do osobnego pliku wyjściowego (_np._ "Hello-output.txt")
2. **`convertToUpper`:** Konwertuje każde powitanie na wielkie litery (_np._ "HELLO")
3. **`collectGreetings`:** Zbiera wszystkie powitania pisane wielkimi literami do jednego pliku wsadowego
4. **`cowpy`:** Generuje grafikę ASCII przy użyciu narzędzia `cowpy`

Konfiguracja workflow'a umożliwia dostarczanie danych wejściowych i parametrów w elastyczny, powtarzalny sposób.

### Zdobyte umiejętności

Dzięki temu praktycznemu kursowi nauczyłeś się:

- Opisywać i wykorzystywać podstawowe komponenty Nextflow'a wystarczające do zbudowania prostego, wieloetapowego workflow'a
- Opisywać koncepcje kolejnego poziomu, takie jak operatory i fabryki kanałów
- Uruchamiać workflow Nextflow lokalnie
- Znajdować i interpretować wyniki oraz pliki dziennika generowane przez Nextflow'a
- Rozwiązywać podstawowe problemy

Posiadasz teraz fundamentalną wiedzę potrzebną do rozpoczęcia tworzenia własnych pipeline'ów w Nextflow.

## Kolejne kroki w rozwijaniu Twoich umiejętności

Oto nasze 3 najlepsze sugestie, co zrobić dalej:

- Zastosuj Nextflow'a do naukowego przypadku użycia analizy z kursem [Nextflow for Science](../nf4_science/index.md)
- Rozpocznij pracę z nf-core dzięki kursowi [Hello nf-core](../hello_nf-core/index.md)
- Poznaj bardziej zaawansowane funkcje Nextflow'a dzięki [Side Quests](../side_quests/index.md)

Na koniec polecamy zapoznać się z [**Seqera Platform**](https://seqera.io/) – platformą chmurową stworzoną przez twórców Nextflow'a, która jeszcze bardziej ułatwia uruchamianie i zarządzanie workflow'ami, a także zarządzanie danymi i interaktywne przeprowadzanie analiz w dowolnym środowisku.

## Ankieta zwrotna

Zanim przejdziesz dalej, poświęć chwilę na wypełnienie ankiety kursu! Twoja opinia pomaga nam ulepszać materiały szkoleniowe dla wszystkich.

[Wypełnij ankietę :material-arrow-right:](survey.md){ .md-button .md-button--primary }
