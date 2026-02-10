# Podsumowanie kursu

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Tłumaczenie wspomagane przez AI - [dowiedz się więcej i zasugeruj ulepszenia](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

Gratulacje z okazji ukończenia kursu szkoleniowego Nextflow Run! 🎉

<!-- placeholder for video -->

## Twoja droga

Zacząłeś od bardzo prostego workflow'a i nauczyłeś się go uruchamiać, znajdować wyniki oraz zarządzać jego wykonaniem.
Następnie przepracowałeś coraz bardziej złożone wersje tego workflow'a i nauczyłeś się rozpoznawać podstawowe koncepcje i mechanizmy napędzające pipeline'y Nextflow, w tym kanały i operatory, modularyzację kodu oraz kontenery.
Na koniec nauczyłeś się dostosowywać konfigurację pipeline'a do Twoich preferencji i infrastruktury obliczeniowej.

### Czego się nauczyłeś

Jesteś teraz w stanie zarządzać wykonaniem pipeline'a Hello, opisać jego strukturę i zidentyfikować główne fragmenty zaangażowanego kodu.

- Ostateczna forma workflow'a Hello przyjmuje jako wejście plik CSV zawierający tekstowe powitania.
- Cztery kroki są zaimplementowane jako procesy Nextflow (`sayHello`, `convertToUpper`, `collectGreetings` i `cowpy`) przechowywane w oddzielnych plikach modułów.
- Wyniki są publikowane do katalogu o nazwie `results/`.
- Końcowym wyjściem pipeline'a jest plik tekstowy zawierający grafikę ASCII postaci wypowiadającej powitania pisane wielkimi literami.

<figure class="excalidraw">
--8<-- "docs/en/docs/hello_nextflow/img/hello_pipeline_complete.svg"
</figure>

1. **`sayHello`:** Zapisuje każde powitanie do własnego pliku wyjściowego (_np._ "Hello-output.txt")
2. **`convertToUpper`:** Konwertuje każde powitanie na wielkie litery (_np._ "HELLO")
3. **`collectGreetings`:** Zbiera wszystkie powitania pisane wielkimi literami do jednego pliku wsadowego
4. **`cowpy`:** Generuje grafikę ASCII przy użyciu narzędzia `cowpy`

Konfiguracja workflow'a umożliwia dostarczanie danych wejściowych i parametrów w elastyczny, powtarzalny sposób.

### Zdobyte umiejętności

Dzięki temu praktycznemu kursowi nauczyłeś się:

- Uruchamiać workflow Nextflow lokalnie
- Znajdować i interpretować wyniki oraz pliki dziennika generowane przez Nextflow
- Rozpoznawać podstawowe komponenty Nextflow tworzące prosty wieloetapowy workflow
- Opisywać zaawansowane koncepcje, takie jak operatory i fabryki kanałów
- Konfigurować pipeline'y dla różnych środowisk obliczeniowych

Jesteś teraz wyposażony w fundamentalną wiedzę potrzebną do rozpoczęcia integrowania istniejących pipeline'ów Nextflow z Twoją własną pracą.

## Kolejne kroki w rozwijaniu umiejętności

Oto nasze najlepsze sugestie, co zrobić dalej:

- Nie tylko uruchamiaj Nextflow, ale go pisz! Zostań programistą Nextflow dzięki [Hello Nextflow](../hello_nextflow/index.md)
- Zastosuj Nextflow do naukowego przypadku użycia z [Nextflow for Science](../nf4_science/index.md)
- Rozpocznij pracę z nf-core dzięki [Hello nf-core](../hello_nf-core/index.md)
- Naucz się technik rozwiązywania problemów dzięki [Debugging Side Quest](../side_quests/debugging.md)

Na koniec zalecamy zapoznanie się z [**Seqera Platform**](https://seqera.io/) — platformą chmurową opracowaną przez twórców Nextflow, która jeszcze bardziej ułatwia uruchamianie i zarządzanie workflow'ami, a także zarządzanie danymi i interaktywne przeprowadzanie analiz w dowolnym środowisku.

## Uzyskiwanie pomocy

Aby uzyskać zasoby pomocy i wsparcie społeczności, zobacz [stronę pomocy](../help.md).

## Ankieta zwrotna

Zanim przejdziesz dalej, poświęć chwilę na wypełnienie ankiety kursu! Twoja opinia pomaga nam ulepszać materiały szkoleniowe dla wszystkich.

[Wypełnij ankietę :material-arrow-right:](survey.md){ .md-button .md-button--primary }
