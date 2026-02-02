# Podsumowanie kursu

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Tłumaczenie wspomagane przez AI - [dowiedz się więcej i zasugeruj ulepszenia](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

Gratulacje z okazji ukończenia kursu szkoleniowego Nextflow Run!

<!-- placeholder for video -->

## Twoja droga

Zacząłeś od bardzo prostego workflow'u, ucząc się go uruchamiać, znajdować wyjścia i zarządzać wykonaniem.
Następnie pracowałeś przez coraz bardziej złożone warianty i poznałeś kluczowe koncepcje oraz mechanizmy napędzające pipeline'y Nextflow, w tym kanały i operatory, modularyzację kodu oraz kontenery.
Na koniec nauczyłeś się, jak dostosować konfigurację pipeline'u do swoich preferencji i infrastruktury obliczeniowej.

### Czego się nauczyłeś

Jesteś teraz w stanie zarządzać wykonaniem pipeline'u Hello, opisać, jak jest zbudowany, i zidentyfikować główne elementy kodu.

- Końcowa forma workflow'u Hello przyjmuje jako dane wejściowe plik CSV zawierający tekstowe powitania.
- Cztery kroki są zaimplementowane jako procesy Nextflow (`sayHello`, `convertToUpper`, `collectGreetings` i `cowpy`) przechowywane w oddzielnych plikach modułów.
- Wyniki są publikowane do katalogu o nazwie `results/`.
- Końcowym wyjściem pipeline'u jest plik tekstowy zawierający grafikę ASCII postaci mówiącej powitania wielkimi literami.

<figure class="excalidraw">
--8<-- "docs/en/docs/hello_nextflow/img/hello_pipeline_complete.svg"
</figure>

1. **`sayHello`:** Zapisuje każde powitanie do własnego pliku wyjściowego (_np._ "Hello-output.txt")
2. **`convertToUpper`:** Konwertuje każde powitanie na wielkie litery (_np._ "HELLO")
3. **`collectGreetings`:** Zbiera wszystkie powitania wielkimi literami do jednego pliku batch
4. **`cowpy`:** Generuje grafikę ASCII przy użyciu narzędzia `cowpy`

Konfiguracja workflow'u wspiera dostarczanie danych wejściowych i parametrów w elastyczny, odtwarzalny sposób.

### Zdobyte umiejętności

Dzięki temu praktycznemu kursowi nauczyłeś się:

- Uruchamiać workflow'y Nextflow lokalnie
- Znajdować i interpretować wyjścia (wyniki) i pliki dziennika generowane przez Nextflow
- Rozpoznawać podstawowe komponenty Nextflow, które tworzą prosty wieloetapowy workflow
- Opisywać koncepcje następnego poziomu, takie jak operatory i fabryki kanałów
- Konfigurować pipeline'y dla różnych środowisk obliczeniowych

Jesteś teraz wyposażony w podstawową wiedzę, aby zacząć integrować istniejące pipeline'y Nextflow ze Swoją własną pracą.

## Następne kroki do rozwijania umiejętności

Oto nasze najważniejsze sugestie, co robić dalej:

- Nie tylko uruchamiaj Nextflow, pisz go! Zostań programistą Nextflow z [Hello Nextflow](../hello_nextflow/index.md)
- Zastosuj Nextflow do naukowego przypadku użycia z [Nextflow for Science](../nf4_science/index.md)
- Rozpocznij pracę z nf-core z [Hello nf-core](../hello_nf-core/index.md)
- Poznaj techniki rozwiązywania problemów z [Debugging Side Quest](../side_quests/debugging.md)

Na koniec zalecamy zapoznanie się z [**Seqera Platform**](https://seqera.io/), platformą chmurową rozwijaną przez twórców Nextflow, która jeszcze bardziej ułatwia uruchamianie i zarządzanie workflow'ami, a także zarządzanie danymi i interaktywne uruchamianie analiz w dowolnym środowisku.

## Uzyskiwanie pomocy

Zasoby pomocy i wsparcie społeczności znajdziesz na [stronie Pomocy](../help.md).

## Ankieta zwrotna

Zanim przejdziesz dalej, poświęć chwilę na wypełnienie ankiety kursu! Twoja opinia pomaga nam ulepszać nasze materiały szkoleniowe dla wszystkich.

[Wypełnij ankietę :material-arrow-right:](survey.md){ .md-button .md-button--primary }
