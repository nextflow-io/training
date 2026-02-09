# Nextflow dla genomiki

**Praktyczny kurs stosowania Nextflow w rzeczywistym przypadku użycia z genomiki: wykrywanie wariantów za pomocą GATK.**

Ten kurs opiera się na szkoleniu dla początkujących [Hello Nextflow](../../hello_nextflow/) i pokazuje, jak używać Nextflow w specyficznym kontekście dziedziny genomiki.
Zaimplementujesz pipeline wykrywania wariantów za pomocą [GATK](https://gatk.broadinstitute.org/) (Genome Analysis Toolkit), szeroko stosowanego pakietu oprogramowania do analizy danych sekwencjonowania wysokoprzepustowego.

<!-- additional_information -->

## Przegląd kursu

Ten kurs ma charakter praktyczny, z ćwiczeniami ukierunkowanymi na cel, ustrukturyzowanymi tak, aby stopniowo wprowadzać informacje.

Zaczniesz od ręcznego uruchamiania narzędzi do wykrywania wariantów w terminalu, aby zrozumieć metodologię, a następnie stopniowo zbudujesz pipeline Nextflow, który automatyzuje i skaluje analizę.

### Plan lekcji

Podzieliliśmy to na trzy części, z których każda koncentruje się na konkretnych aspektach stosowania Nextflow w przypadku użycia z genomiki.

| Rozdział kursu                                                           | Podsumowanie                                                                                                  | Szacowany czas |
| ------------------------------------------------------------------------ | ------------------------------------------------------------------------------------------------------------- | -------------- |
| [Część 1: Przegląd metody](./01_method.md)                               | Zrozumienie metodologii wykrywania wariantów i ręczne uruchamianie narzędzi                                   | 30 min         |
| [Część 2: Wykrywanie wariantów dla pojedynczych próbek](./02_per_sample_variant_calling.md) | Budowanie pipeline'u, który indeksuje pliki BAM i wykrywa warianty, a następnie skalowanie do wielu próbek    | 60 min         |
| [Część 3: Wspólne wywoływanie na kohorcie](./03_joint_calling.md)        | Dodawanie wielopróbkowego wspólnego genotypowania przy użyciu operatorów kanałów do agregacji wyników próbek  | 45 min         |

Pod koniec tego kursu będziesz w stanie zastosować podstawowe koncepcje i narzędzia Nextflow do typowego przypadku użycia z genomiki.

Gotowy, aby rozpocząć kurs?

[Zacznij :material-arrow-right:](00_orientation.md){ .md-button .md-button--primary }

<!-- Clearfix for float -->
<div style="content: ''; clear: both; display: table;"></div>
