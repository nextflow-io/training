# Podsumowanie kursu

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Tłumaczenie wspomagane przez AI - [dowiedz się więcej i zasugeruj ulepszenia](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

Gratulacje ukończenia kursu szkoleniowego Nextflow dla genomiki! 🎉

## Twoja droga

Zacząłeś od ręcznego uruchamiania narzędzi do wykrywania wariantów w terminalu, aby zrozumieć metodologię.
Następnie zbudowałeś pipeline Nextflow'a dla pojedynczej próbki, aby zautomatyzować proces, przeskalowałeś go do obsługi wielu próbek równolegle i dodałeś wspólne genotypowanie wielu próbek przy użyciu operatorów kanałów.

### Co zbudowałeś

- Pipeline do wykrywania wariantów, który przyjmuje pliki BAM jako wejście i produkuje wspólnie wywołane pliki VCF jako wyjście.
- Trzy procesy (`SAMTOOLS_INDEX`, `GATK_HAPLOTYPECALLER` i `GATK_JOINTGENOTYPING`) przechowywane w oddzielnych plikach modułów.
- Pipeline skaluje się automatycznie do dowolnej liczby próbek wejściowych dzięki paradygmatowi przepływu danych Nextflow'a.
- Wyniki są publikowane do katalogu o nazwie `results/`.

### Zdobyte umiejętności

Dzięki temu praktycznemu kursowi nauczyłeś się:

- Pisać liniowy workflow do zastosowania wykrywania wariantów dla pojedynczej próbki
- Odpowiednio obsługiwać pliki pomocnicze, takie jak pliki indeksów i zasoby genomu referencyjnego
- Wykorzystywać paradygmat przepływu danych Nextflow'a do paralelizacji wykrywania wariantów dla poszczególnych próbek
- Implementować wspólne wywoływanie wielu próbek przy użyciu odpowiednich operatorów kanałów

Jesteś teraz przygotowany, aby zacząć stosować Nextflow'a do workflow'ów analizy genomicznej we własnej pracy.

## Kolejne kroki do rozwijania umiejętności

Oto nasze najlepsze sugestie, co zrobić dalej:

- Zastosuj Nextflow'a do innych przypadków użycia analizy naukowej z [Nextflow for Science](../index.md)
- Rozpocznij pracę z nf-core dzięki [Hello nf-core](../../hello_nf-core/index.md)
- Poznaj bardziej zaawansowane funkcje Nextflow'a dzięki [Side Quests](../../side_quests/index.md)

Na koniec zalecamy zapoznanie się z [**Seqera Platform**](https://seqera.io/), platformą chmurową opracowaną przez twórców Nextflow'a, która jeszcze bardziej ułatwia uruchamianie i zarządzanie workflow'ami, a także zarządzanie danymi i interaktywne przeprowadzanie analiz w dowolnym środowisku.

## Uzyskiwanie pomocy

Aby uzyskać zasoby pomocy i wsparcie społeczności, zobacz [stronę pomocy](../../help.md).

## Ankieta zwrotna

Zanim przejdziesz dalej, poświęć chwilę na wypełnienie ankiety kursu! Twoja opinia pomaga nam ulepszać materiały szkoleniowe dla wszystkich.

[Wypełnij ankietę :material-arrow-right:](survey.md){ .md-button .md-button--primary }
