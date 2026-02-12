# Podsumowanie kursu

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Tłumaczenie wspomagane przez AI - [dowiedz się więcej i zasugeruj ulepszenia](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

Gratulujemy ukończenia kursu szkoleniowego Nextflow dla RNAseq!

## Twoja droga

Zacząłeś od ręcznego uruchamiania narzędzi do przetwarzania RNAseq w terminalu, aby zrozumieć metodologię.
Następnie zbudowałeś pipeline Nextflow dla pojedynczej próbki, aby zautomatyzować proces, przeskalowałeś go do obsługi wielu próbek równolegle i rozszerzyłeś, aby obsługiwał dane paired-end oraz agregował raporty QC dla wszystkich próbek.

### Co zbudowałeś

- Pipeline do przetwarzania RNAseq, który przyjmuje pliki FASTQ jako wejście i produkuje przycięte odczyty, dopasowania oraz zagregowane raporty QC jako wyjście.
- Procesy do przycinania (Trim Galore), dopasowania (HISAT2), kontroli jakości (FastQC) i agregacji raportów (MultiQC) przechowywane w oddzielnych plikach modułów.
- Pipeline automatycznie paralelizuje przetwarzanie próbek wejściowych, wykorzystując paradygmat przepływu danych Nextflow.
- Finalny pipeline obsługuje dane sekwencjonowania paired-end.

### Zdobyte umiejętności

Dzięki temu praktycznemu kursowi nauczyłeś się:

- Pisać liniowy workflow do zastosowania podstawowych metod przetwarzania i kontroli jakości RNAseq
- Odpowiednio obsługiwać pliki specyficzne dla dziedziny, takie jak FASTQ i zasoby genomu referencyjnego
- Obsługiwać dane sekwencjonowania single-end i paired-end
- Wykorzystywać paradygmat przepływu danych Nextflow do paralelizacji przetwarzania RNAseq dla poszczególnych próbek
- Agregować raporty QC z wielu kroków i próbek, używając odpowiednich operatorów kanałów

Jesteś teraz przygotowany, aby zacząć stosować Nextflow do workflow'ów analizy RNAseq we własnej pracy.

## Następne kroki w rozwijaniu umiejętności

Oto nasze najważniejsze sugestie dotyczące kolejnych kroków:

- Zastosuj Nextflow do innych przypadków użycia w analizie naukowej z [Nextflow for Science](../index.md)
- Rozpocznij pracę z nf-core dzięki [Hello nf-core](../../hello_nf-core/index.md)
- Poznaj bardziej zaawansowane funkcje Nextflow dzięki [Side Quests](../../side_quests/index.md)

Na koniec polecamy zapoznanie się z **[Seqera Platform](https://seqera.io/)**, platformą w chmurze opracowaną przez twórców Nextflow, która jeszcze bardziej ułatwia uruchamianie i zarządzanie workflow'ami, a także zarządzanie danymi i przeprowadzanie analiz interaktywnie w dowolnym środowisku.

## Uzyskiwanie pomocy

Aby znaleźć zasoby pomocy i wsparcie społeczności, zobacz [stronę Pomocy](../../help.md).

## Ankieta zwrotna

Zanim przejdziesz dalej, poświęć chwilę na wypełnienie ankiety kursu! Twoja opinia pomaga nam ulepszać materiały szkoleniowe dla wszystkich.

[Wypełnij ankietę :material-arrow-right:](survey.md){ .md-button .md-button--primary }
