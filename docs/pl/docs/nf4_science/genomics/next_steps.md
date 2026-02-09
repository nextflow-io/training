# Podsumowanie kursu

Gratulacje z okazji ukończenia kursu szkoleniowego Nextflow dla Genomiki! 🎉

## Twoja droga

Zaczęłaś/-eś od ręcznego uruchamiania narzędzi do wykrywania wariantów w terminalu, aby zrozumieć metodologię.
Następnie zbudowałaś/-eś pipeline Nextflow'a dla pojedynczej próbki, aby zautomatyzować proces, przeskalowałaś/-eś go do obsługi wielu próbek równolegle i dodałaś/-eś wielopróbkowe wspólne genotypowanie przy użyciu operatorów kanałów.

### Co zbudowałaś/-eś

- Pipeline do wykrywania wariantów, który przyjmuje pliki BAM jako wejście i produkuje wspólnie wywołane pliki VCF jako wyjście.
- Trzy procesy (`SAMTOOLS_INDEX`, `GATK_HAPLOTYPECALLER` i `GATK_JOINTGENOTYPING`) przechowywane w oddzielnych plikach modułów.
- Pipeline automatycznie skaluje się do dowolnej liczby próbek wejściowych dzięki paradygmatowi przepływu danych Nextflow'a.
- Wyniki są publikowane do katalogu o nazwie `results/`.

### Zdobyte umiejętności

Dzięki temu praktycznemu kursowi nauczyłaś/-eś się, jak:

- Napisać liniowy workflow, aby zastosować wykrywanie wariantów do pojedynczej próbki
- Odpowiednio obsługiwać pliki pomocnicze, takie jak pliki indeksu i zasoby genomu referencyjnego
- Wykorzystać paradygmat przepływu danych Nextflow'a do zrównoleglenia wykrywania wariantów per próbka
- Zaimplementować wielopróbkowe wspólne wywoływanie przy użyciu odpowiednich operatorów kanałów

Jesteś teraz wyposażona/-y, aby zacząć stosować Nextflow'a do workflow'ów analizy genomicznej we własnej pracy.

## Następne kroki w rozwijaniu umiejętności

Oto nasze najlepsze sugestie, co robić dalej:

- Zastosuj Nextflow'a do innych przypadków użycia w analizie naukowej z [Nextflow dla Nauki](../index.md)
- Rozpocznij pracę z nf-core dzięki [Hello nf-core](../../hello_nf-core/index.md)
- Odkryj bardziej zaawansowane funkcje Nextflow'a z [Side Quests](../../side_quests/index.md)

Na koniec zalecamy zapoznanie się z [**Seqera Platform**](https://seqera.io/) – platformą chmurową opracowaną przez twórców Nextflow'a, która jeszcze bardziej ułatwia uruchamianie i zarządzanie workflow'ami, a także zarządzanie danymi i przeprowadzanie analiz interaktywnie w dowolnym środowisku.

## Uzyskiwanie pomocy

W celu uzyskania zasobów pomocy i wsparcia społeczności zobacz [stronę Pomocy](../../help.md).

## Ankieta feedbackowa

Zanim przejdziesz dalej, poświęć chwilę na wypełnienie ankiety kursu! Twoja opinia pomaga nam ulepszać materiały szkoleniowe dla wszystkich.

[Wypełnij ankietę :material-arrow-right:](survey.md){ .md-button .md-button--primary }
