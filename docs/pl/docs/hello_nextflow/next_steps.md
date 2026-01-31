# Podsumowanie kursu

Gratulacje z okazji ukończenia kursu szkoleniowego Hello Nextflow! 🎉

<div class="video-wrapper">
  <iframe width="560" height="315" src="https://www.youtube.com/embed/xHOcx_4Ancg?si=Lp8hS8RdaMwbp5j5&amp;list=PLPZ8WHdZGxmXiHf8B26oB_fTfoKQdhlik" title="YouTube video player" frameborder="0" allow="accelerometer; autoplay; clipboard-write; encrypted-media; gyroscope; picture-in-picture; web-share" referrerpolicy="strict-origin-when-cross-origin" allowfullscreen></iframe>
</div>

/// caption
:fontawesome-brands-youtube:{ .youtube } Zobacz [całą playlistę na kanale YouTube Nextflow](https://www.youtube.com/playlist?list=PLPZ8WHdZGxmXiHf8B26oB_fTfoKQdhlik).

:green_book: Możesz przeczytać [transkrypcję wideo](./transcripts/07_next_steps.md) wraz z filmem.
///

## Twoja podróż

Zacząłeś od bardzo podstawowego workflow'u, który uruchamiał zakodowane na sztywno polecenie.
W ciągu sześciu części przekształciłeś ten podstawowy workflow w modularny, wieloetapowy pipeline, który wykorzystuje kluczowe funkcje Nextflow, w tym kanały, operatory, wbudowaną obsługę kontenerów i opcje konfiguracji.

### Co zbudowałeś

- Ostateczna forma workflow'u Hello przyjmuje jako wejście plik CSV zawierający tekstowe pozdrowienia.
- Cztery kroki są zaimplementowane jako procesy Nextflow (`sayHello`, `convertToUpper`, `collectGreetings` i `cowpy`) przechowywane w osobnych plikach modułów.
- Wyniki są publikowane do katalogu o nazwie `results/`.
- Końcowe wyjście pipeline'u to plik tekstowy zawierający grafikę ASCII postaci wypowiadającej pozdrowienia zapisane wielkimi literami.

<figure class="excalidraw">
--8<-- "docs/hello_nextflow/img/hello_pipeline_complete.svg"
</figure>

1. **`sayHello`:** Zapisuje każde pozdrowienie do własnego pliku wyjściowego (_np._ "Hello-output.txt")
2. **`convertToUpper`:** Konwertuje każde pozdrowienie na wielkie litery (_np._ "HELLO")
3. **`collectGreetings`:** Zbiera wszystkie pozdrowienia z wielkimi literami do jednego pliku batch
4. **`cowpy`:** Generuje grafikę ASCII za pomocą narzędzia `cowpy`

Konfiguracja workflow'u wspiera dostarczanie wejść i parametrów w elastyczny, powtarzalny sposób.

### Nabyte umiejętności

Przez ten praktyczny kurs nauczyłeś się:

- Opisywać i wykorzystywać podstawowe komponenty Nextflow wystarczające do budowania prostego, wieloetapowego workflow'u
- Opisywać koncepcje kolejnego kroku, takie jak operatory i channel factories
- Uruchamiać workflow Nextflow lokalnie
- Znajdować i interpretować wyjścia (wyniki) i pliki dziennika generowane przez Nextflow
- Rozwiązywać podstawowe problemy

Jesteś teraz wyposażony w fundamentalną wiedzę, aby zacząć tworzyć własne pipeline'y w Nextflow.

## Kolejne kroki do rozwijania umiejętności

Oto nasze 3 najlepsze sugestie, co robić dalej:

- Zastosuj Nextflow do naukowego przypadku analizy z [Nextflow dla nauki](../nf4_science/index.md)
- Rozpocznij pracę z nf-core z [Hello nf-core](../../hello_nf-core/index.md)
- Odkryj bardziej zaawansowane funkcje Nextflow z [Side Quests](../side_quests/index.md)

Na koniec polecamy zapoznać się z [**Seqera Platform**](https://seqera.io/), platformą chmurową opracowaną przez twórców Nextflow, która jeszcze bardziej ułatwia uruchamianie i zarządzanie workflow'ami, a także zarządzanie danymi i interaktywne uruchamianie analiz w dowolnym środowisku.

## Ankieta zwrotna

Zanim przejdziesz dalej, poświęć chwilę na wypełnienie ankiety kursu! Twoja opinia pomaga nam ulepszać materiały szkoleniowe dla wszystkich.

[Wypełnij ankietę :material-arrow-right:](survey.md){ .md-button .md-button--primary }
