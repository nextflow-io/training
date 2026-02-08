# Część 2: Wspólne wywołanie wariantów dla kohorty

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Tłumaczenie wspomagane przez AI - [dowiedz się więcej i zasugeruj ulepszenia](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

W pierwszej części tego kursu zbudowałeś pipeline do wykrywania wariantów, który był całkowicie liniowy i przetwarzał dane każdej próbki niezależnie od innych.
Jednak w rzeczywistym przypadku użycia genomiki zazwyczaj będziesz musiał spojrzeć na wywołania wariantów wielu próbek razem.

W tej drugiej części pokażemy Ci, jak używać kanałów i operatorów kanałów do implementacji wspólnego wykrywania wariantów za pomocą GATK, budując na pipeline'ie z Części 1.

### Przegląd metody

Metoda GATK, której użyliśmy w pierwszej części tego kursu, generowała wyniki dla każdej próbki osobno.
Jest to w porządku, jeśli chcesz analizować warianty z każdej próbki indywidualnie, ale daje to ograniczone informacje.
Często bardziej interesujące jest porównanie wyników między wieloma próbkami, a GATK oferuje do tego alternatywną metodę zwaną wspólnym genotypowaniem, którą tutaj demonstrujemy.

Wspólne wywoływanie wariantów polega na wygenerowaniu specjalnego rodzaju wyjścia wariantów zwanego GVCF (Genomic VCF) dla każdej próbki, następnie połączeniu danych GVCF ze wszystkich próbek i wreszcie uruchomieniu statystycznej analizy 'wspólnego genotypowania'.

![Wspólna analiza](img/joint-calling.png)

To, co jest specjalnego w GVCF próbki, to to, że zawiera rekordy podsumowujące statystyki danych sekwencyjnych o wszystkich pozycjach w docelowym obszarze genomu, nie tylko pozycjach, gdzie program znalazł dowody zmienności.
Jest to kluczowe dla obliczenia wspólnego genotypowania ([dalsza lektura](https://gatk.broadinstitute.org/hc/en-us/articles/360035890431-The-logic-of-joint-calling-for-germline-short-variants)).

GVCF jest produkowany przez GATK HaplotypeCaller, to samo narzędzie, którego użyliśmy w Części 1, z dodatkowym parametrem (`-ERC GVCF`).
Łączenie GVCF odbywa się za pomocą GATK GenomicsDBImport, który łączy wywołania dla poszczególnych próbek w magazyn danych (analogiczny do bazy danych), a następnie właściwa analiza 'wspólnego genotypowania' jest wykonywana za pomocą GATK GenotypeGVCFs.

### Workflow

Podsumowując, w tej części kursu zamierzamy rozwinąć workflow, który wykonuje następujące czynności:

<figure class="excalidraw">
--8<-- "docs/en/docs/nf4_science/genomics/img/hello-gatk-2.svg"
</figure>

1. Wygenerować plik indeksu dla każdego pliku BAM wejściowego za pomocą Samtools
2. Uruchomić GATK HaplotypeCaller na każdym pliku BAM wejściowym, aby wygenerować GVCF wywołań wariantów genomowych dla każdej próbki
3. Zebrać wszystkie GVCF i połączyć je w magazyn danych GenomicsDB
4. Uruchomić wspólne genotypowanie na połączonym magazynie danych GVCF, aby wygenerować plik VCF na poziomie kohorty

Zastosujemy to do tego samego zestawu danych, co w Części 1.

---

## 0. Rozgrzewka: Uruchom Samtools i GATK bezpośrednio

Tak jak wcześniej, chcemy wypróbować polecenia ręcznie, zanim spróbujemy opakować je w workflow'ie.

!!! note

     Upewnij się, że jesteś we właściwym katalogu roboczym:
     `cd /workspaces/training/nf4-science/genomics`

### 0.1. Zaindeksuj plik BAM wejściowy za pomocą Samtools

Ten pierwszy krok jest taki sam jak w Części 1, więc powinien być bardzo znajomy, ale tym razem musimy to zrobić dla wszystkich trzech próbek.

!!! note

    Technicznie już wygenerowaliśmy pliki indeksu dla trzech próbek przez nasz pipeline, więc moglibyśmy je wydobyć z katalogu wyników. Jednak czystsze jest po prostu powtórzenie tego ręcznie, a zajmie to tylko minutę.

#### 0.1.1. Uruchom kontener Samtools interaktywnie

```bash
docker run -it -v ./data:/data community.wave.seqera.io/library/samtools:1.20--b5dfbd93de237464
```

<!--
??? success "Wyjście polecenia"

    ```console

    ```
-->

#### 0.1.2. Uruchom polecenie indeksowania dla trzech próbek

```bash
samtools index /data/bam/reads_mother.bam
samtools index /data/bam/reads_father.bam
samtools index /data/bam/reads_son.bam
```

Tak jak wcześniej, to powinno wytworzyć pliki indeksu w tym samym katalogu co odpowiadające pliki BAM.

??? abstract "Zawartość katalogu"

    ```console
    data/bam/
    ├── reads_father.bam
    ├── reads_father.bam.bai
    ├── reads_mother.bam
    ├── reads_mother.bam.bai
    ├── reads_son.bam
    └── reads_son.bam.bai
    ```

Teraz, gdy mamy pliki indeksu dla wszystkich trzech próbek, możemy przejść do generowania GVCF dla każdej z nich.

#### 0.1.3. Wyjdź z kontenera Samtools

```bash
exit
```

### 0.2. Wywołaj warianty za pomocą GATK HaplotypeCaller w trybie GVCF

Ten drugi krok jest bardzo podobny do tego, co zrobiliśmy w Części 1: Hello Genomics, ale teraz zamierzamy uruchomić GATK w 'trybie GVCF'.

#### 0.2.1. Uruchom kontener GATK interaktywnie

```bash
docker run -it -v ./data:/data community.wave.seqera.io/library/gatk4:4.5.0.0--730ee8817e436867
```

<!--
??? success "Wyjście polecenia"

    ```console

    ```
-->

#### 0.2.2. Uruchom polecenie wywoływania wariantów z opcją GVCF

Aby wytworzyć genomiczny VCF (GVCF), dodajemy opcję `-ERC GVCF` do podstawowego polecenia, co włącza tryb GVCF HaplotypeCaller'a.

Zmieniamy również rozszerzenie pliku wyjściowego z `.vcf` na `.g.vcf`.
Technicznie nie jest to wymagane, ale jest to zdecydowanie zalecana konwencja.

```bash
gatk HaplotypeCaller \
        -R /data/ref/ref.fasta \
        -I /data/bam/reads_mother.bam \
        -O reads_mother.g.vcf \
        -L /data/ref/intervals.bed \
        -ERC GVCF
```

<!--
??? success "Wyjście polecenia"

    ```console

    ```
-->

To tworzy plik wyjściowy GVCF `reads_mother.g.vcf` w bieżącym katalogu roboczym w kontenerze.

Jeśli użyjesz `cat`, aby wyświetlić zawartość, zobaczysz, że jest on znacznie dłuższy niż odpowiadający VCF, który wygenerowaliśmy w Części 1. Nie możesz nawet przewinąć do początku pliku, a większość linii wygląda zupełnie inaczej niż to, co widzieliśmy w VCF w Części 1.

```console title="Wyjście" linenums="1674"
20_10037292_10066351    14714   .       T       <NON_REF>       .       .       END=14718       GT:DP:GQ:MIN_DP:PL       0/0:37:99:37:0,99,1192
20_10037292_10066351    14719   .       T       <NON_REF>       .       .       END=14719       GT:DP:GQ:MIN_DP:PL       0/0:36:82:36:0,82,1087
20_10037292_10066351    14720   .       T       <NON_REF>       .       .       END=14737       GT:DP:GQ:MIN_DP:PL       0/0:42:99:37:0,100,1160
```

Reprezentują one regiony nie-wariantowe, gdzie wywołujący warianty nie znalazł dowodów zmienności, więc uchwycił pewne statystyki opisujące jego poziom pewności w braku zmienności. Umożliwia to rozróżnienie między dwoma bardzo różnymi przypadkami: (1) są dobre jakościowo dane pokazujące, że próbka jest homozygotyczna-referencyjna, i (2) nie ma wystarczająco dobrych danych dostępnych, aby dokonać określenia w jakikolwiek sposób.

W GVCF zazwyczaj jest wiele takich linii nie-wariantowych, z mniejszą liczbą rekordów wariantów rozproszonymi wśród nich. Spróbuj uruchomić `head -176` na GVCF, aby załadować tylko pierwsze 176 linii pliku i znaleźć rzeczywiste wywołanie wariantu.

```console title="Wyjście" linenums="174"
20_10037292_10066351    3479    .       T       <NON_REF>       .       .       END=3479        GT:DP:GQ:MIN_DP:PL       0/0:34:36:34:0,36,906
20_10037292_10066351    3480    .       C       CT,<NON_REF>    503.03  .       DP=23;ExcessHet=0.0000;MLEAC=2,0;MLEAF=1.00,0.00;RAW_MQandDP=82800,23    GT:AD:DP:GQ:PL:SB       1/1:0,18,0:18:54:517,54,0,517,54,517:0,0,7,11
20_10037292_10066351    3481    .       T       <NON_REF>       .       .       END=3481        GT:DP:GQ:MIN_DP:PL       0/0:21:51:21:0,51,765
```

Druga linia pokazuje pierwszy rekord wariantu w pliku, który odpowiada pierwszemu wariantowi w pliku VCF, na który patrzyliśmy w Części 1.

Tak jak oryginalny VCF, plik wyjściowy GVCF jest również dołączony do pliku indeksu, zwanego `reads_mother.g.vcf.idx`.

#### 0.2.3. Powtórz proces na dwóch pozostałych próbkach

Aby przetestować krok wspólnego genotypowania, potrzebujemy GVCF dla wszystkich trzech próbek, więc wygenerujmy je teraz ręcznie.

```bash
gatk HaplotypeCaller \
        -R /data/ref/ref.fasta \
        -I /data/bam/reads_father.bam \
        -O reads_father.g.vcf \
        -L /data/ref/intervals.bed \
        -ERC GVCF
```

<!--
??? success "Wyjście polecenia"

    ```console

    ```
-->

```bash
gatk HaplotypeCaller \
        -R /data/ref/ref.fasta \
        -I /data/bam/reads_son.bam \
        -O reads_son.g.vcf \
        -L /data/ref/intervals.bed \
        -ERC GVCF
```

<!--
??? success "Wyjście polecenia"

    ```console

    ```
-->

Po zakończeniu będziesz mieć trzy pliki kończące się na `.g.vcf` w Twoim bieżącym katalogu (jeden na próbkę) i ich odpowiednie pliki indeksu kończące się na `.g.vcf.idx`.

### 0.3. Uruchom wspólne genotypowanie

Teraz, gdy mamy wszystkie GVCF, możemy wreszcie wypróbować podejście wspólnego genotypowania do generowania wywołań wariantów dla kohorty próbek.
Jako przypomnienie, jest to dwuetapowa metoda, która polega na połączeniu danych ze wszystkich GVCF w magazyn danych, a następnie uruchomieniu właściwej analizy wspólnego genotypowania, aby wygenerować końcowy VCF wspólnie wywołanych wariantów.

#### 0.3.1. Połącz wszystkie GVCF dla poszczególnych próbek

Ten pierwszy krok używa innego narzędzia GATK, zwanego GenomicsDBImport, aby połączyć dane ze wszystkich GVCF w magazyn danych GenomicsDB.

```bash
gatk GenomicsDBImport \
    -V reads_mother.g.vcf \
    -V reads_father.g.vcf \
    -V reads_son.g.vcf \
    -L /data/ref/intervals.bed \
    --genomicsdb-workspace-path family_trio_gdb
```

<!--
??? success "Wyjście polecenia"

    ```console

    ```
-->

Wyjściem tego kroku jest faktycznie katalog zawierający zestaw dalej zagnieżdżonych katalogów przechowujących połączone dane wariantów w postaci wielu różnych plików.
Możesz w nim podziurać, ale szybko zobaczysz, że ten format magazynu danych nie jest przeznaczony do bezpośredniego czytania przez ludzi.

!!! note

    GATK zawiera narzędzia, które umożliwiają inspekcję i wydobycie danych wywołań wariantów z magazynu danych w razie potrzeby.

#### 0.3.2. Uruchom właściwą analizę wspólnego genotypowania

Ten drugi krok używa jeszcze innego narzędzia GATK, zwanego GenotypeGVCFs, aby przeliczyć statystyki wariantów i indywidualne genotypy w świetle danych dostępnych we wszystkich próbkach w kohorcie.

```bash
gatk GenotypeGVCFs \
    -R /data/ref/ref.fasta \
    -V gendb://family_trio_gdb \
    -O family_trio.vcf
```

<!--
??? success "Wyjście polecenia"

    ```console

    ```
-->

To tworzy plik wyjściowy VCF `family_trio.vcf` w bieżącym katalogu roboczym w kontenerze.
To kolejny dość mały plik, więc możesz użyć `cat` na tym pliku, aby wyświetlić jego zawartość i przewinąć w górę, aby znaleźć kilka pierwszych linii wariantów.

```console title="family_trio.vcf" linenums="40"
#CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO    FORMAT  reads_father    reads_mother    reads_son
20_10037292_10066351    3480    .       C       CT      1625.89 .       AC=5;AF=0.833;AN=6;BaseQRankSum=0.220;DP=85;ExcessHet=0.0000;FS=2.476;MLEAC=5;MLEAF=0.833;MQ=60.00;MQRankSum=0.00;QD=21.68;ReadPosRankSum=-1.147e+00;SOR=0.487    GT:AD:DP:GQ:PL  0/1:15,16:31:99:367,0,375       1/1:0,18:18:54:517,54,0 1/1:0,26:26:78:756,78,0
20_10037292_10066351    3520    .       AT      A       1678.89 .       AC=5;AF=0.833;AN=6;BaseQRankSum=1.03;DP=80;ExcessHet=0.0000;FS=2.290;MLEAC=5;MLEAF=0.833;MQ=60.00;MQRankSum=0.00;QD=22.39;ReadPosRankSum=0.701;SOR=0.730 GT:AD:DP:GQ:PL   0/1:18,13:31:99:296,0,424       1/1:0,18:18:54:623,54,0 1/1:0,26:26:78:774,78,0
20_10037292_10066351    3529    .       T       A       154.29  .       AC=1;AF=0.167;AN=6;BaseQRankSum=-5.440e-01;DP=104;ExcessHet=0.0000;FS=1.871;MLEAC=1;MLEAF=0.167;MQ=60.00;MQRankSum=0.00;QD=7.71;ReadPosRankSum=-1.158e+00;SOR=1.034       GT:AD:DP:GQ:PL  0/0:44,0:44:99:0,112,1347       0/1:12,8:20:99:163,0,328        0/0:39,0:39:99:0,105,1194
```

To wygląda bardziej jak oryginalny VCF, który wygenerowaliśmy w Części 1, z tym że tym razem mamy informacje na poziomie genotypu dla wszystkich trzech próbek.
Ostatnie trzy kolumny w pliku to bloki genotypu dla próbek, wymienione w porządku alfabetycznym.

Jeśli spojrzymy na genotypy wywołane dla naszego testowego tria rodzinnego dla pierwszego wariantu, widzimy, że ojciec jest heterozygotyczny-wariantowy (`0/1`), a matka i syn są obaj homozygotyczni-wariantowi (`1/1`).

To jest ostatecznie informacja, którą chcemy wydobyć z zestawu danych! Więc opakowujmy to wszystko w workflow Nextflow'a, abyśmy mogli robić to na dużą skalę.

#### 0.3.3. Wyjdź z kontenera GATK

```bash
exit
```

### Wnioski

Wiesz, jak uruchamiać poszczególne polecenia zaangażowane we wspólne wywoływanie wariantów w terminalu, aby sprawdzić, czy wyprodukują informacje, których chcesz.

### Co dalej?

Opakuj te polecenia w rzeczywisty pipeline.

---

## 1. Zmodyfikuj krok wywoływania wariantów dla poszczególnych próbek, aby wytwarzał GVCF

Dobra wiadomość jest taka, że nie musimy zaczynać od początku, ponieważ już napisaliśmy workflow, który wykonuje część tej pracy w Części 1.
Jednak ten pipeline produkuje pliki VCF, podczas gdy teraz chcemy plików GVCF, aby wykonać wspólne genotypowanie.
Więc musimy zacząć od włączenia trybu wywoływania wariantów GVCF i zaktualizowania rozszerzenia pliku wyjściowego.

!!! note

    Dla wygody będziemy pracować z nową kopią workflow'u GATK, jak stoi na końcu Części 1, ale pod inną nazwą: `genomics-2.nf`.

### 1.1. Powiedz HaplotypeCaller'owi, aby emitował GVCF i zaktualizuj rozszerzenie wyjściowe

Otwórzmy plik `genomics-2.nf` w edytorze kodu.
Powinien wyglądać bardzo znajomo, ale możesz go uruchomić, jeśli chcesz upewnić się, że działa zgodnie z oczekiwaniami.

Zamierzamy zacząć od wprowadzenia dwóch zmian:

- Dodać parametr `-ERC GVCF` do polecenia GATK HaplotypeCaller;
- Zaktualizować ścieżkę pliku wyjściowego, aby używała odpowiadającego rozszerzenia `.g.vcf`, zgodnie z konwencją GATK.

Upewnij się, że dodajesz ukośnik wsteczny (`\`) na końcu poprzedniej linii, gdy dodajesz `-ERC GVCF`.

=== "Po"

    ```groovy title="genomics-2.nf" linenums="56" hl_lines="4 6"
        """
        gatk HaplotypeCaller \
            -R ${ref_fasta} \
            -I ${input_bam} \
            -O ${input_bam}.g.vcf \
            -L ${interval_list} \
            -ERC GVCF
        """
    ```

=== "Przed"

    ```groovy title="genomics-2.nf" linenums="56" hl_lines="4"
        """
        gatk HaplotypeCaller \
            -R ${ref_fasta} \
            -I ${input_bam} \
            -O ${input_bam}.vcf \
            -L ${interval_list}
        """
    ```

I to wszystko, czego potrzeba, aby przełączyć HaplotypeCaller'a na generowanie GVCF zamiast VCF, prawda?

### 1.2. Uruchom pipeline, aby sprawdzić, czy możesz generować GVCF

Polecenie wykonania Nextflow'a jest takie samo jak wcześniej, z wyjątkiem samej nazwy pliku workflow'u.
Upewnij się, że odpowiednio zaktualizowałeś nazwę pliku.

```bash
nextflow run genomics-2.nf
```

??? success "Wyjście polecenia"

    ```console
    N E X T F L O W   ~  version 25.10.2

    ┃ Launching `genomics-2.nf` [crazy_venter] DSL2 - revision: a2d6f6f09f

    executor >  local (6)
    [f1/8d8486] SAMTOOLS_INDEX (1)       | 3 of 3 ✔
    [72/3249ca] GATK_HAPLOTYPECALLER (3) | 0 of 3
    ERROR ~ Error executing process > 'GATK_HAPLOTYPECALLER (2)'

    Caused by:
      Missing output file(s) `reads_son.bam.vcf` expected by process `GATK_HAPLOTYPECALLER (2)`

    Command executed:

      gatk HaplotypeCaller         -R ref.fasta         -I reads_son.bam         -O reads_son.bam.g.vcf         -L intervals.bed         -ERC GVCF
    ```

A wyjście jest... całe czerwone! O nie.

Polecenie, które zostało wykonane, jest poprawne, więc mieliśmy rację, że to wystarczyło do zmiany zachowania narzędzia GATK.
Ale spójrz na tę linię o brakującym pliku wyjściowym. Zauważyłeś coś?

Tak jest, zapomnieliśmy powiedzieć Nextflow'owi, aby oczekiwał nowej nazwy pliku. Ups.

### 1.3. Zaktualizuj również rozszerzenie pliku wyjściowego w bloku wyjść procesu

Ponieważ nie wystarczy po prostu zmienić rozszerzenie pliku w samym poleceniu narzędzia, musisz również powiedzieć Nextflow'owi, że oczekiwana nazwa pliku wyjściowego się zmieniła.

=== "Po"

    ```groovy title="genomics-2.nf" linenums="50" hl_lines="2 3"
        output:
        path "${input_bam}.g.vcf"     , emit: vcf
        path "${input_bam}.g.vcf.idx" , emit: idx
    ```

=== "Przed"

    ```groovy title="genomics-2.nf" linenums="50" hl_lines="2 3"
        output:
        path "${input_bam}.vcf"     , emit: vcf
        path "${input_bam}.vcf.idx" , emit: idx
    ```

### 1.4. Zaktualizuj cele publikacji dla nowych wyjść GVCF

Ponieważ teraz produkujemy GVCF zamiast VCF, powinniśmy zaktualizować sekcję `publish:` workflow'u, aby używać bardziej opisowych nazw.
Zorganizujemy również pliki GVCF w ich własnym podkatalogu dla jasności.

=== "Po"

    ```groovy title="genomics-2.nf" linenums="88" hl_lines="3 4"
        publish:
        indexed_bam = SAMTOOLS_INDEX.out
        gvcf = GATK_HAPLOTYPECALLER.out.vcf
        gvcf_idx = GATK_HAPLOTYPECALLER.out.idx
    ```

=== "Przed"

    ```groovy title="genomics-2.nf" linenums="88"
        publish:
        indexed_bam = SAMTOOLS_INDEX.out
        vcf = GATK_HAPLOTYPECALLER.out.vcf
        vcf_idx = GATK_HAPLOTYPECALLER.out.idx
    ```

### 1.5. Zaktualizuj blok wyjściowy dla nowej struktury katalogów

Musimy również zaktualizować blok `output`, aby umieścić pliki GVCF w podkatalogu `gvcf`.

=== "Po"

    ```groovy title="genomics-2.nf" linenums="94" hl_lines="3 5 6 8 9"
    output {
        indexed_bam {
            path 'indexed_bam'
        }
        gvcf {
            path 'gvcf'
        }
        gvcf_idx {
            path 'gvcf'
        }
    }
    ```

=== "Przed"

    ```groovy title="genomics-2.nf" linenums="94"
    output {
        indexed_bam {
            path '.'
        }
        vcf {
            path '.'
        }
        vcf_idx {
            path '.'
        }
    }
    ```

### 1.6. Uruchom pipeline ponownie

Uruchommy go tym razem z `-resume`.

```bash
nextflow run genomics-2.nf -resume
```

??? success "Wyjście polecenia"

    ```console
    N E X T F L O W   ~  version 25.10.2

    ┃ Launching `genomics-2.nf` [nostalgic_franklin] DSL2 - revision: f2c0a93c6a

    executor >  local (3)
    [cc/fbc705] SAMTOOLS_INDEX (3)       | 3 of 3, cached: 3 ✔
    [27/0d7eb9] GATK_HAPLOTYPECALLER (2) | 3 of 3 ✔
    ```

Tym razem działa.

Same wyjście Nextflow'a nie wygląda inaczej (w porównaniu do pomyślnego uruchomienia w normalnym trybie VCF), ale teraz możemy znaleźć pliki `.g.vcf` i ich odpowiednie pliki indeksu, dla wszystkich trzech próbek, zorganizowane w podkatalogach.

??? abstract "Zawartość katalogu (dowiązania symboliczne skrócone)"

    ```console
    results_genomics/
    ├── gvcf/
    │   ├── reads_father.bam.g.vcf -> */27/0d7eb9*/reads_father.bam.g.vcf
    │   ├── reads_father.bam.g.vcf.idx -> */27/0d7eb9*/reads_father.bam.g.vcf.idx
    │   ├── reads_mother.bam.g.vcf -> */e4/4ed55e*/reads_mother.bam.g.vcf
    │   ├── reads_mother.bam.g.vcf.idx -> */e4/4ed55e*/reads_mother.bam.g.vcf.idx
    │   ├── reads_son.bam.g.vcf -> */08/e95962*/reads_son.bam.g.vcf
    │   └── reads_son.bam.g.vcf.idx -> */08/e95962*/reads_son.bam.g.vcf.idx
    └── indexed_bam/
        ├── reads_father.bam -> */9a/c7a873*/reads_father.bam
        ├── reads_father.bam.bai -> */9a/c7a873*/reads_father.bam.bai
        ├── reads_mother.bam -> */f1/8d8486*/reads_mother.bam
        ├── reads_mother.bam.bai -> */f1/8d8486*/reads_mother.bam.bai
        ├── reads_son.bam -> */cc/fbc705*/reads_son.bam
        └── reads_son.bam.bai -> */cc/fbc705*/reads_son.bam.bai
    ```

Jeśli otworzysz jeden z plików GVCF i przewiniesz go, możesz sprawdzić, że GATK HaplotypeCaller wytworzyło pliki GVCF zgodnie z życzeniem.

### Wnioski

Okej, ten był minimalny pod względem nauki Nextflow'a...
Ale to była miła okazja, aby powtórzyć znaczenie bloku wyjściowego procesu!

### Co dalej?

Naucz się zbierać zawartość kanału i przekazywać je do następnego procesu jako pojedyncze wejście.

---

## 2. Zbierz i połącz dane GVCF ze wszystkich próbek

Teraz musimy połączyć dane ze wszystkich GVCF dla poszczególnych próbek w formę, która wspiera analizę wspólnego genotypowania, którą chcemy wykonać.

### 2.1. Zdefiniuj proces, który połączy GVCF

Jako przypomnienie tego, co zrobiliśmy wcześniej w sekcji rozgrzewkowej, łączenie GVCF jest zadaniem dla narzędzia GATK GenomicsDBImport, które wytworzy magazyn danych w tak zwanym formacie GenomicsDB.

Napiszmy nowy proces, aby zdefiniować, jak to będzie działać, na podstawie polecenia, którego użyliśmy wcześniej w sekcji rozgrzewkowej.

```groovy title="genomics-2.nf" linenums="66"
/*
 * Połącz GVCF w magazyn danych GenomicsDB
 */
process GATK_GENOMICSDB {

    container "community.wave.seqera.io/library/gatk4:4.5.0.0--730ee8817e436867"

    input:
    path all_gvcfs
    path all_idxs
    path interval_list
    val cohort_name

    output:
    path "${cohort_name}_gdb"

    script:
    """
    gatk GenomicsDBImport \
        -V ${all_gvcfs} \
        -L ${interval_list} \
        --genomicsdb-workspace-path ${cohort_name}_gdb
    """
}
```

Co myślisz, wygląda rozsądnie?

Podłączmy to i zobaczmy, co się stanie.

### 2.2. Dodaj parametr `cohort_name` z wartością domyślną

Musimy podać arbitralną nazwę dla kohorty.
Później w serii szkoleń nauczysz się używać metadanych próbek do tego rodzaju rzeczy, ale na razie po prostu deklarujemy parametr CLI używając `params` i dajemy mu wartość domyślną dla wygody.

```groovy title="genomics-2.nf" linenums="16"
    // Nazwa bazowa dla końcowego pliku wyjściowego
    cohort_name: String = "family_trio"
```

### 2.3. Zbierz wyjścia GATK_HAPLOTYPECALLER ze wszystkich próbek

Gdybyśmy mieli po prostu podłączyć kanał wyjściowy z procesu `GATK_HAPLOTYPECALLER` tak jak jest, Nextflow wywołałby proces na każdym GVCF próbki osobno.
Jednak chcemy zgrupować wszystkie trzy GVCF (i ich pliki indeksu) w taki sposób, aby Nextflow przekazał wszystkie razem do pojedynczego wywołania procesu.

Dobra wiadomość: możemy to zrobić za pomocą operatora kanału `collect()`. Dodajmy następujące linie do ciała `workflow`, zaraz po wywołaniu GATK_HAPLOTYPECALLER:

```groovy title="genomics-2.nf" linenums="118"
// Zbierz wyjścia wywoływania wariantów ze wszystkich próbek
all_gvcfs_ch = GATK_HAPLOTYPECALLER.out.vcf.collect()
all_idxs_ch = GATK_HAPLOTYPECALLER.out.idx.collect()
```

Czy to wydaje się trochę skomplikowane? Rozbijmy to i przetłumaczmy na prosty język.

1. Bierzemy kanał wyjściowy z procesu `GATK_HAPLOTYPECALLER`, odnoszony za pomocą właściwości `.out`.
2. Każdy 'element' wychodzący z kanału to para plików: GVCF i jego plik indeksu, w tej kolejności, ponieważ to jest kolejność, w jakiej są wymienione w bloku wyjściowym procesu. Wygodnie, ponieważ w ostatniej sesji nazwaliśmy wyjścia tego procesu (używając `emit:`), możemy wybrać GVCF z jednej strony, dodając `.vcf`, a pliki indeksu z drugiej strony, dodając `.idx` po właściwości `.out`. Gdybyśmy nie nazwali tych wyjść, musielibyśmy odnieść się do nich odpowiednio jako `.out[0]` i `.out[1]`.
3. Dołączamy operator kanału `collect()`, aby zgrupować wszystkie pliki GVCF razem w pojedynczy element w nowym kanale zwanym `all_gvcfs_ch`, i robimy to samo z plikami indeksu, aby utworzyć nowy kanał zwany `all_idxs_ch`.

!!! tip

    Jeśli masz trudności z wyobrażeniem sobie dokładnie, co się tutaj dzieje, pamiętaj, że możesz użyć operatora `view()`, aby sprawdzić zawartość kanałów przed i po zastosowaniu operatorów kanału.

Powstałe kanały `all_gvcfs_ch` i `all_idxs_ch` to te, które zamierzamy podłączyć do procesu `GATK_GENOMICSDB`, który właśnie napisaliśmy.

!!! note

    W przypadku, gdybyś się zastanawiał, zbieramy GVCF i ich pliki indeksu osobno, ponieważ polecenie GATK GenomicsDBImport chce widzieć tylko ścieżki plików GVCF. Na szczęście, ponieważ Nextflow będzie przygotowywał wszystkie pliki razem do wykonania, nie musimy martwić się o kolejność plików, jak to zrobiliśmy w przypadku BAM i ich indeksów w Części 1.

### 2.4. Dodaj wywołanie do bloku workflow, aby uruchomić GATK_GENOMICSDB

Mamy proces i mamy kanały wejściowe. Musimy tylko dodać wywołanie procesu.

```groovy title="genomics-2.nf" linenums="122"
    // Połącz GVCF w magazyn danych GenomicsDB
    GATK_GENOMICSDB(
        all_gvcfs_ch,
        all_idxs_ch,
        intervals_file,
        params.cohort_name
    )
```

Ok, wszystko jest podłączone.

### 2.5. Uruchom workflow

Zobaczmy, czy to działa.

```bash
nextflow run genomics-2.nf -resume
```

??? failure "Wyjście polecenia"

    ```console
    N E X T F L O W   ~  version 25.10.2

    ┃ Launching `genomics-2.nf` [disturbed_bell] DSL2 - revision: 57942246cc

    executor >  local (1)
    [f1/8d8486] SAMTOOLS_INDEX (1)       | 3 of 3, cached: 3 ✔
    [e4/4ed55e] GATK_HAPLOTYPECALLER (2) | 3 of 3, cached: 3 ✔
    [51/d350ea] GATK_GENOMICSDB          | 0 of 1
    ERROR ~ Error executing process > 'GATK_GENOMICSDB'

    Caused by:
      Process `GATK_GENOMICSDB` terminated with an error exit status (1)

    Command executed:

      gatk GenomicsDBImport         -V reads_son.bam.g.vcf reads_father.bam.g.vcf reads_mother.bam.g.vcf         -L intervals.bed         --genomicsdb-workspace-path family_trio_gdb
    ```

Uruchamia się dość szybko, ponieważ uruchamiamy z `-resume`, ale zawodzi!

Ach. Z jasnej strony, widzimy, że Nextflow pobrał proces `GATK_GENOMICSDB` i konkretnie wywołał go tylko raz.
To sugeruje, że podejście `collect()` zadziałało, do pewnego stopnia.
Ale, i to duże, wywołanie procesu nie powiodło się.

Kiedy wchodzimy w głąb wyjścia konsoli powyżej, możemy zobaczyć, że wykonane polecenie nie jest poprawne.

Czy możesz zauważyć błąd?
Spójrz na ten fragment: `-V reads_father.bam.g.vcf reads_son.bam.g.vcf reads_mother.bam.g.vcf`

Daliśmy `gatk GenomicsDBImport` wiele plików GVCF dla pojedynczego argumentu `-V`, ale narzędzie oczekuje oddzielnego argumentu `-V` dla każdego pliku GVCF.

Jako przypomnienie, to było polecenie, które uruchomiliśmy w kontenerze:

```bash
gatk GenomicsDBImport \
    -V reads_mother.g.vcf \
    -V reads_father.g.vcf \
    -V reads_son.g.vcf \
    -L /data/ref/intervals.bed \
    --genomicsdb-workspace-path family_trio_gdb
```

Więc to oznacza, że musimy jakoś przekształcić nasz pakiet plików GVCF w poprawnie sformatowany ciąg poleceń.

### 2.6. Skonstruuj linię poleceń z osobnym argumentem `-V` dla każdego wejściowego GVCF

To jest miejsce, gdzie Nextflow bazujący na Groovy okazuje się przydatny, ponieważ pozwoli nam użyć dość prostych manipulacji ciągami, aby skonstruować niezbędny ciąg poleceń.

Konkretnie, używając tej składni: `all_gvcfs.collect { gvcf -> "-V ${gvcf}" }.join(' ')`

Jeszcze raz rozbijmy to na składniki.

1. Najpierw bierzemy zawartość kanału wejściowego `all_gvcfs` i stosujemy na nim `.collect()` (tak jak wcześniej).
2. Pozwala nam to przekazać każdą indywidualną ścieżkę pliku GVCF w pakiecie do **domknięcia**, `{ gvcf -> "-V ${gvcf}" }`, gdzie `gvcf` odnosi się do tej ścieżki pliku GVCF.
   Domknięcie to mini-funkcja, której używamy, aby poprzedzić `-V ` do ścieżki pliku, w postaci `"-V ${gvcf}"`.
3. Następnie używamy `.join(' ')`, aby połączyć wszystkie trzy ciągi z pojedynczą spacją jako separatorem.

Z konkretnym przykładem wygląda to tak:

1. Mamy trzy pliki:

   `[A.ext, B.ext, C.ext]`

2. Domknięcie modyfikuje każdy, aby utworzyć ciągi:

   `"-V A.ext", "-V B.ext", "-V C.ext"`

3. Operacja `.join(' ')` generuje końcowy ciąg:

   `"-V A.ext -V B.ext -V C.ext"`

Kiedy mamy ten ciąg, możemy przypisać go do zmiennej lokalnej, `gvcfs_line`, zdefiniowanej słowem kluczowym `def`:

`def gvcfs_line = all_gvcfs.collect { gvcf -> "-V ${gvcf}" }.join(' ')`

Ok, mamy więc naszą rzecz do manipulacji ciągami. Gdzie ją umieścić?

Chcemy, aby to trafiło gdzieś do definicji procesu, ponieważ chcemy to zrobić _po_ tym, jak przekierowaliśmy ścieżki plików GVCF do procesu.
To dlatego, że Nextflow musi je zobaczyć jako ścieżki plików, aby poprawnie przygotować same pliki do wykonania.

Ale _gdzie_ w procesie możemy to dodać?

Ciekawy fakt: możesz dodać dowolny kod po `script:` i przed `"""` !

Świetnie, dodajmy więc tam naszą linię manipulacji ciągami, i zaktualizujmy polecenie `gatk GenomicsDBImport`, aby używało łączonego ciągu, który produkuje.

=== "Po"

    ```groovy title="genomics-2.nf" linenums="87"  hl_lines="2 5"
        script:
        def gvcfs_line = all_gvcfs.collect { gvcf -> "-V ${gvcf}" }.join(' ')
        """
        gatk GenomicsDBImport \
            ${gvcfs_line} \
            -L ${interval_list} \
            --genomicsdb-workspace-path ${cohort_name}_gdb
        """
    ```

=== "Przed"

    ```groovy title="genomics-2.nf" linenums="87"  hl_lines="4"
        script:
        """
        gatk GenomicsDBImport \
            -V ${all_gvcfs} \
            -L ${interval_list} \
            --genomicsdb-workspace-path ${cohort_name}_gdb
        """
    ```

To powinno być wszystko, czego potrzeba, aby prawidłowo dostarczyć wejścia do `gatk GenomicsDBImport`.

!!! tip

    Kiedy aktualizujesz polecenie `gatk GenomicsDBImport`, upewnij się, że usuniesz prefiks `-V `, gdy zamieniasz na zmienną `${gvcfs_line}`.

### 2.7. Uruchom workflow, aby sprawdzić, czy generuje wyjście GenomicsDB zgodnie z oczekiwaniami

W porządku, zobaczmy, czy to rozwiązało problem.

```bash
nextflow run genomics-2.nf -resume
```

??? success "Wyjście polecenia"

    ```console
    N E X T F L O W   ~  version 25.10.2

    ┃ Launching `genomics-2.nf` [peaceful_gates] DSL2 - revision: ca0bf847ed

    executor >  local (1)
    [cc/fbc705] SAMTOOLS_INDEX (3)       | 3 of 3, cached: 3 ✔
    [27/0d7eb9] GATK_HAPLOTYPECALLER (2) | 3 of 3, cached: 3 ✔
    [76/d13861] GATK_GENOMICSDB          | 1 of 1 ✔
    ```

Aha! Wygląda na to, że tym razem działa.

Pierwsze dwa kroki zostały pomyślnie pominięte, a trzeci krok zadziałał tym razem jak za dotknięciem różdżki.
Magazyn danych GenomicsDB jest tworzony w katalogu roboczym, ale nie jest publikowany do wyników, ponieważ to tylko format pośredni, którego użyjemy do wspólnego genotypowania.

Nawiasem mówiąc, nie musieliśmy robić nic specjalnego, aby obsłużyć wyjście będące katalogiem zamiast pojedynczym plikiem.

### Wnioski

Teraz wiesz, jak zbierać wyjścia z kanału i grupować je jako pojedyncze wejście do innego procesu.
Wiesz również, jak skonstruować linię poleceń, aby dostarczyć wejścia do danego narzędzia z odpowiednią składnią.

### Co dalej?

Naucz się, jak dodać drugie polecenie do tego samego procesu.

---

## 3. Uruchom krok wspólnego genotypowania jako część tego samego procesu

Teraz, gdy mamy połączone genomowe wywołania wariantów, możemy uruchomić narzędzie do wspólnego genotypowania, które wytworzy końcowe wyjście, na którym nam naprawdę zależy: VCF wywołań wariantów na poziomie kohorty.

Ze względów logistycznych decydujemy się włączyć wspólne genotypowanie do tego samego procesu.

### 3.1. Zmień nazwę procesu z GATK_GENOMICSDB na GATK_JOINTGENOTYPING

Ponieważ proces będzie uruchamiał więcej niż jedno narzędzie, zmieniamy jego nazwę, aby odnosiła się do całej operacji, a nie nazwy pojedynczego narzędzia.

=== "Po"

    ```groovy title="genomics-2.nf"
    /*
     * Połącz GVCF w magazyn danych GenomicsDB i uruchom wspólne genotypowanie, aby wytworzyć wywołania na poziomie kohorty
     */
    process GATK_JOINTGENOTYPING {
    ```

=== "Przed"

    ```groovy title="genomics-2.nf"
    /*
     * Połącz GVCF w magazyn danych GenomicsDB
     */
    process GATK_GENOMICSDB {
    ```

Pamiętaj, aby nazwy procesów były jak najbardziej opisowe, aby zmaksymalizować czytelność dla kolegów — i siebie w przyszłości!

### 3.2. Dodaj polecenie wspólnego genotypowania do procesu GATK_JOINTGENOTYPING

Po prostu dodaj drugie polecenie po pierwszym wewnątrz sekcji script.

=== "Po"

    ```groovy title="genomics-2.nf" linenums="89"  hl_lines="6-10"
        """
        gatk GenomicsDBImport \
            ${gvcfs_line} \
            -L ${interval_list} \
            --genomicsdb-workspace-path ${cohort_name}_gdb

        gatk GenotypeGVCFs \
            -R ${ref_fasta} \
            -V gendb://${cohort_name}_gdb \
            -L ${interval_list} \
            -O ${cohort_name}.joint.vcf
        """
    ```

=== "Przed"

    ```groovy title="genomics-2.nf" linenums="89"
        """
        gatk GenomicsDBImport \
            ${gvcfs_line} \
            -L ${interval_list} \
            --genomicsdb-workspace-path ${cohort_name}_gdb
        """
    ```

Dwa polecenia będą uruchamiane sekwencyjnie, w taki sam sposób, jak gdybyśmy mieli uruchomić je ręcznie w terminalu.

### 3.3. Dodaj pliki genomu referencyjnego do definicji wejść procesu GATK_JOINTGENOTYPING

Drugie polecenie wymaga plików genomu referencyjnego, więc musimy dodać je do wejść procesu.

=== "Po"

    ```groovy title="genomics-2.nf" linenums="78"  hl_lines="6-8"
        input:
        path all_gvcfs
        path all_idxs
        path interval_list
        val cohort_name
        path ref_fasta
        path ref_index
        path ref_dict
    ```

=== "Przed"

    ```groovy title="genomics-2.nf" linenums="78"
        input:
        path all_gvcfs
        path all_idxs
        path interval_list
        val cohort_name
    ```

Może wydawać się irytujące pisanie tego, ale pamiętaj, że piszesz to tylko raz, a następnie możesz uruchomić workflow milion razy. Warte tego?

### 3.4. Zaktualizuj definicję wyjścia procesu, aby emitowała VCF wywołań wariantów na poziomie kohorty

Naprawdę nie zależy nam na zapisaniu magazynu danych GenomicsDB, który jest tylko formatem pośrednim istniejącym z przyczyn logistycznych, więc możemy po prostu usunąć go z bloku wyjściowego, jeśli chcemy.

Wyjściem, na którym nam naprawdę zależy, jest VCF wytworzony przez polecenie wspólnego genotypowania.

=== "Po"

    ```groovy title="genomics-2.nf" linenums="87" hl_lines="2 3"
        output:
        path "${cohort_name}.joint.vcf"     , emit: vcf
        path "${cohort_name}.joint.vcf.idx" , emit: idx
    ```

=== "Przed"

    ```groovy title="genomics-2.nf" linenums="87"
        output:
        path "${cohort_name}_gdb"
    ```

Prawie skończyliśmy!

### 3.5. Zaktualizuj wywołanie procesu z GATK_GENOMICSDB na GATK_JOINTGENOTYPING

Nie zapomnijmy zmienić nazwy wywołania procesu w ciele workflow'u z GATK_GENOMICSDB na GATK_JOINTGENOTYPING. A skoro już przy tym jesteśmy, powinniśmy również dodać pliki genomu referencyjnego jako wejścia, ponieważ musimy je dostarczyć do narzędzia wspólnego genotypowania.

=== "Po"

    ```groovy title="genomics-2.nf" linenums="126"
    // Połącz GVCF w magazyn danych GenomicsDB i zastosuj wspólne genotypowanie
    GATK_JOINTGENOTYPING(
        all_gvcfs_ch,
        all_idxs_ch,
        intervals_file,
        params.cohort_name,
        ref_file,
        ref_index_file,
        ref_dict_file
    )
    ```

=== "Przed"

    ```groovy title="genomics-2.nf" linenums="126"
    // Połącz GVCF w magazyn danych GenomicsDB
    GATK_GENOMICSDB(
        all_gvcfs_ch,
        all_idxs_ch,
        intervals_file,
        params.cohort_name
    )
    ```

Teraz proces jest całkowicie podłączony.

### 3.6. Dodaj wspólny VCF do sekcji publikacji

Musimy opublikować wspólne wyjścia VCF z nowego procesu.
Dodaj te linie do sekcji `publish:` workflow'u:

```groovy title="genomics-2.nf" linenums="145"
    joint_vcf = GATK_JOINTGENOTYPING.out.vcf
    joint_vcf_idx = GATK_JOINTGENOTYPING.out.idx
```

### 3.7. Dodaj cele wspólnego VCF do bloku wyjściowego

Na koniec dodaj cele wyjściowe dla wspólnych plików VCF.
Umieścimy je w głównym katalogu wyników, ponieważ to jest końcowe wyjście.

```groovy title="genomics-2.nf" linenums="157"
    joint_vcf {
        path '.'
    }
    joint_vcf_idx {
        path '.'
    }
```

Teraz wszystko powinno być całkowicie podłączone.

### 3.8. Uruchom workflow

Wreszcie możemy uruchomić zmodyfikowany workflow...

```bash
nextflow run genomics-2.nf -resume
```

??? success "Wyjście polecenia"

    ```console
    N E X T F L O W   ~  version 25.10.2

    ┃ Launching `genomics-2.nf` [crazy_marconi] DSL2 - revision: 5da9afc841

    executor >  local (1)
    [9a/c7a873] SAMTOOLS_INDEX (2)       | 3 of 3, cached: 3 ✔
    [e4/4ed55e] GATK_HAPLOTYPECALLER (2) | 3 of 3, cached: 3 ✔
    [a6/7cc8ed] GATK_JOINTGENOTYPING     | 1 of 1 ✔
    ```

I to działa!

Znajdziesz końcowy plik wyjściowy, `family_trio.joint.vcf` (i jego indeks pliku), w katalogu wyników.

??? abstract "Zawartość katalogu (dowiązania symboliczne skrócone)"

    ```console
    results_genomics/
    ├── family_trio.joint.vcf -> */a6/7cc8ed*/family_trio.joint.vcf
    ├── family_trio.joint.vcf.idx -> */a6/7cc8ed*/family_trio.joint.vcf.idx
    ├── gvcf/
    │   ├── reads_father.bam.g.vcf -> */27/0d7eb9*/reads_father.bam.g.vcf
    │   ├── reads_father.bam.g.vcf.idx -> */27/0d7eb9*/reads_father.bam.g.vcf.idx
    │   ├── reads_mother.bam.g.vcf -> */e4/4ed55e*/reads_mother.bam.g.vcf
    │   ├── reads_mother.bam.g.vcf.idx -> */e4/4ed55e*/reads_mother.bam.g.vcf.idx
    │   ├── reads_son.bam.g.vcf -> */08/e95962*/reads_son.bam.g.vcf
    │   └── reads_son.bam.g.vcf.idx -> */08/e95962*/reads_son.bam.g.vcf.idx
    └── indexed_bam/
        ├── reads_father.bam -> */9a/c7a873*/reads_father.bam
        ├── reads_father.bam.bai -> */9a/c7a873*/reads_father.bam.bai
        ├── reads_mother.bam -> */f1/8d8486*/reads_mother.bam
        ├── reads_mother.bam.bai -> */f1/8d8486*/reads_mother.bam.bai
        ├── reads_son.bam -> */cc/fbc705*/reads_son.bam
        └── reads_son.bam.bai -> */cc/fbc705*/reads_son.bam.bai
    ```

Jeśli jesteś typem sceptycznym, możesz kliknąć na wspólny plik VCF, aby go otworzyć i sprawdzić, że workflow wygenerował te same wywołania wariantów, które otrzymałeś, uruchamiając narzędzia ręcznie na początku tej sekcji.

```console title="family_trio.joint.vcf" linenums="40"
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	reads_father	reads_mother	reads_son
20_10037292_10066351	3480	.	C	CT	1625.89	.	AC=5;AF=0.833;AN=6;BaseQRankSum=0.220;DP=85;ExcessHet=0.0000;FS=2.476;MLEAC=5;MLEAF=0.833;MQ=60.00;MQRankSum=0.00;QD=21.68;ReadPosRankSum=-1.147e+00;SOR=0.487	GT:AD:DP:GQ:PL	0/1:15,16:31:99:367,0,375	1/1:0,18:18:54:517,54,0	1/1:0,26:26:78:756,78,0
20_10037292_10066351	3520	.	AT	A	1678.89	.	AC=5;AF=0.833;AN=6;BaseQRankSum=1.03;DP=80;ExcessHet=0.0000;FS=2.290;MLEAC=5;MLEAF=0.833;MQ=60.00;MQRankSum=0.00;QD=22.39;ReadPosRankSum=0.701;SOR=0.730	GT:AD:DP:GQ:PL	0/1:18,13:31:99:296,0,424	1/1:0,18:18:54:623,54,0	1/1:0,26:26:78:774,78,0
20_10037292_10066351	3529	.	T	A	154.29	.	AC=1;AF=0.167;AN=6;BaseQRankSum=-5.440e-01;DP=104;ExcessHet=0.0000;FS=1.871;MLEAC=1;MLEAF=0.167;MQ=60.00;MQRankSum=0.00;QD=7.71;ReadPosRankSum=-1.158e+00;SOR=1.034	GT:AD:DP:GQ:PL	0/0:44,0:44:99:0,112,1347	0/1:12,8:20:99:163,0,328	0/0:39,0:39:99:0,105,1194
```

Masz teraz zautomatyzowany, w pełni odtwarzalny workflow wspólnego wywoływania wariantów!

!!! note

    Pamiętaj, że pliki danych, które Ci daliśmy, obejmują tylko maleńką część chromosomu 20.
    Rzeczywisty rozmiar zestawu wywołań wariantów byłby liczony w milionach wariantów.
    Dlatego używamy tylko małych podzbiorów danych do celów szkoleniowych!

### Wnioski

Wiesz, jak używać niektórych typowych operatorów, a także domknięć Groovy do kontrolowania przepływu danych w Twoim workflow'ie.

### Co dalej?

Świętuj Twój sukces i weź zasłużoną przerwę.

W następnej części tego kursu nauczysz się, jak modularyzować Twój workflow, wyodrębniając definicje procesów do modułów wielokrotnego użytku.
