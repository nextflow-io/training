# Podsumowanie kursu

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Tłumaczenie wspomagane przez AI - [dowiedz się więcej i zasugeruj ulepszenia](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

Gratulacje z okazji ukończenia kursu szkoleniowego Nextflow dla genomiki! 🎉

## Twoja droga

Zacząłeś od ręcznego uruchamiania narzędzi do wykrywania wariantów w terminalu, aby zrozumieć metodologię.
Następnie zbudowałeś pipeline Nextflow dla pojedynczej próbki, aby zautomatyzować ten proces, przeskalowałeś go do równoległej obsługi wielu próbek i dodałeś wspólne genotypowanie wielu próbek przy użyciu operatorów kanałów.

### Co zbudowałeś

- Pipeline do wykrywania wariantów, który przyjmuje pliki BAM jako wejście i produkuje wspólnie wywołane pliki VCF jako wyjście.
- Trzy procesy (`SAMTOOLS_INDEX`, `GATK_HAPLOTYPECALLER` i `GATK_JOINTGENOTYPING`) przechowywane w oddzielnych plikach modułów.
- Pipeline automatycznie paralelizuje przetwarzanie próbek wejściowych, wykorzystując paradygmat przepływu danych Nextflow.
- Wyniki są publikowane do katalogu o nazwie `results/`.

### Zdobyte umiejętności

Dzięki temu praktycznemu kursowi nauczyłeś się:

- Pisać liniowy workflow do zastosowania wykrywania wariantów dla pojedynczej próbki
- Odpowiednio obsługiwać pliki pomocnicze, takie jak pliki indeksów i zasoby genomu referencyjnego
- Wykorzystywać paradygmat przepływu danych Nextflow do paralelizacji wykrywania wariantów dla poszczególnych próbek
- Implementować wspólne wywoływanie wielu próbek przy użyciu odpowiednich operatorów kanałów

Jesteś teraz przygotowany, aby zacząć stosować Nextflow do workflow'ów analizy genomicznej we własnej pracy.

## Kolejne kroki w rozwijaniu Twoich umiejętności

Oto nasze najlepsze sugestie, co zrobić dalej:

- Zastosuj Nextflow do innych przypadków użycia analizy naukowej z [Nextflow for Science](../index.md)
- Rozpocznij pracę z nf-core dzięki [Hello nf-core](../../hello_nf-core/index.md)
- Poznaj bardziej zaawansowane funkcje Nextflow dzięki [Side Quests](../../side_quests/index.md)

Na koniec zalecamy zapoznanie się z [**Seqera Platform**](https://seqera.io/) — platformą opartą na chmurze, opracowaną przez twórców Nextflow, która jeszcze bardziej ułatwia uruchamianie i zarządzanie Twoimi workflow'ami, a także zarządzanie danymi i interaktywne przeprowadzanie analiz w dowolnym środowisku.

## Uzyskiwanie pomocy

Aby uzyskać zasoby pomocy i wsparcie społeczności, zobacz [stronę pomocy](../../help.md).

## Ankieta zwrotna

Zanim przejdziesz dalej, poświęć chwilę na wypełnienie ankiety dotyczącej kursu! Twoja opinia pomaga nam ulepszać nasze materiały szkoleniowe dla wszystkich.

[Wypełnij ankietę :material-arrow-right:](survey.md){ .md-button .md-button--primary }
