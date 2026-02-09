# Część 1: Przegląd metod i ręczne testowanie

Wykrywanie wariantów to metoda analizy genomowej, która ma na celu identyfikację zmian w sekwencji genomu względem genomu referencyjnego.
Tutaj użyjemy narzędzi i metod zaprojektowanych do wykrywania krótkich wariantów zarodkowych, _tzn._ SNP i indeli, w danych z sekwencjonowania całego genomu.

![Pipeline GATK](img/gatk-pipeline.png)

Pełny pipeline wykrywania wariantów zazwyczaj obejmuje wiele kroków, w tym mapowanie do referencji (czasami określane jako dopasowanie genomu) oraz filtrowanie i priorytetyzację wariantów.
Dla uproszczenia w tym szkoleniu skupimy się tylko na części dotyczącej wykrywania wariantów.

### Metody

Pokażemy Ci dwa sposoby zastosowania wykrywania wariantów do próbek z sekwencjonowania całego genomu w celu identyfikacji zarodkowych SNP i indeli.
Najpierw zaczniemy od prostego **podejścia per-próbka**, które wykrywa warianty niezależnie dla każdej próbki.
Następnie pokażemy Ci bardziej zaawansowane **podejście z łącznym wykrywaniem (joint calling)**, które analizuje wiele próbek razem, dając dokładniejsze i bardziej informacyjne wyniki.

Zanim zaczniemy pisać jakikolwiek kod workflow'a dla któregokolwiek podejścia, przetestujemy polecenia ręcznie na danych testowych.

### Zbiór danych

Udostępniamy następujące dane i powiązane zasoby:

- **Genom referencyjny** składający się z małego regionu ludzkiego chromosomu 20 (z hg19/b37) oraz jego plików pomocniczych (indeks i słownik sekwencji).
- **Trzy próbki z sekwencjonowania całego genomu** odpowiadające trio rodzinnemu (matka, ojciec i syn), które zostały zawężone do małego fragmentu danych na chromosomie 20, aby zachować małe rozmiary plików.
  To dane z sekwencjonowania krótkich odczytów Illumina, które zostały już zmapowane do genomu referencyjnego, dostarczone w formacie [BAM](https://samtools.github.io/hts-specs/SAMv1.pdf) (Binary Alignment Map, skompresowana wersja SAM, Sequence Alignment Map).
- **Lista interwałów genomowych**, czyli koordynatów na genomie, gdzie nasze próbki mają dane odpowiednie do wykrywania wariantów, dostarczone w formacie BED.

### Oprogramowanie

Dwa główne narzędzia to [Samtools](https://www.htslib.org/), szeroko stosowany zestaw narzędzi do manipulacji plikami dopasowań sekwencji, oraz [GATK](https://gatk.broadinstitute.org/) (Genome Analysis Toolkit), zestaw narzędzi do wykrywania wariantów opracowany w Broad Institute.

Te narzędzia nie są zainstalowane w środowisku GitHub Codespaces, więc będziemy ich używać przez kontenery (zobacz [Hello Containers](../../hello_nextflow/05_hello_containers.md)).

!!! note "Uwaga"

     Upewnij się, że jesteś w katalogu `nf4-science/genomics`, aby ostatnia część ścieżki pokazanej po wpisaniu `pwd` to `genomics`.

---

## 1. Wykrywanie wariantów per-próbka

Wykrywanie wariantów per-próbka przetwarza każdą próbkę niezależnie: narzędzie wykrywające warianty analizuje dane sekwencjonowania dla jednej próbki na raz i identyfikuje pozycje, w których próbka różni się od referencji.

W tej sekcji testujemy dwa polecenia składające się na podejście wykrywania wariantów per-próbka: indeksowanie pliku BAM za pomocą Samtools oraz wykrywanie wariantów za pomocą GATK HaplotypeCaller.
To te polecenia, które opakujemy w workflow Nextflow w części 2 tego kursu.

1. Wygeneruj plik indeksu dla wejściowego pliku BAM za pomocą [Samtools](https://www.htslib.org/)
2. Uruchom GATK HaplotypeCaller na zindeksowanym pliku BAM, aby wygenerować wykrycia wariantów per-próbka w formacie VCF (Variant Call Format)

<figure class="excalidraw">
--8<-- "docs/en/docs/nf4_science/genomics/img/hello-gatk-1.svg"
</figure>

Zaczynamy od testowania dwóch poleceń na tylko jednej próbce.

### 1.1. Indeksowanie wejściowego pliku BAM za pomocą Samtools

Pliki indeksów są powszechną cechą formatów plików bioinformatycznych; zawierają informacje o strukturze pliku głównego, które pozwalają narzędziom takim jak GATK na dostęp do podzbioru danych bez konieczności odczytywania całego pliku.
Jest to ważne ze względu na to, jak duże mogą być te pliki.

Pliki BAM są często dostarczane bez indeksu, więc pierwszym krokiem w wielu workflow'ach analizy jest wygenerowanie go za pomocą `samtools index`.

Pobierzemy kontener Samtools, uruchomimy go interaktywnie i wykonamy polecenie `samtools index` na jednym z plików BAM.

#### 1.1.1. Pobierz kontener Samtools

Uruchom polecenie `docker pull`, aby pobrać obraz kontenera Samtools:

```bash
docker pull community.wave.seqera.io/library/samtools:1.20--b5dfbd93de237464
```

??? success "Wyjście polecenia"

    ```console
    1.20--b5dfbd93de237464: Pulling from library/samtools
    6360b3717211: Pull complete
    2ec3f7ad9b3c: Pull complete
    7716ca300600: Pull complete
    4f4fb700ef54: Pull complete
    8c61d418774c: Pull complete
    03dae77ff45c: Pull complete
    aab7f787139d: Pull complete
    4f4fb700ef54: Pull complete
    837d55536720: Pull complete
    897362c12ca7: Pull complete
    3893cbe24e91: Pull complete
    d1b61e94977b: Pull complete
    c72ff66fb90f: Pull complete
    0e0388f29b6d: Pull complete
    Digest: sha256:bbfc45b4f228975bde86cba95e303dd94ecf2fdacea5bfb2e2f34b0d7b141e41
    Status: Downloaded newer image for community.wave.seqera.io/library/samtools:1.20--b5dfbd93de237464
    community.wave.seqera.io/library/samtools:1.20--b5dfbd93de237464
    ```

Jeśli wcześniej nie pobierałeś/aś tego obrazu, ukończenie może zająć minutę.
Po zakończeniu będziesz mieć lokalną kopię obrazu kontenera.

#### 1.1.2. Uruchom kontener Samtools interaktywnie

Aby uruchomić kontener interaktywnie, użyj `docker run` z flagami `-it`.
Opcja `-v ./data:/data` montuje lokalny katalog `data` wewnątrz kontenera, aby narzędzia mogły uzyskać dostęp do plików wejściowych.

```bash
docker run -it -v ./data:/data community.wave.seqera.io/library/samtools:1.20--b5dfbd93de237464
```

Twój wiersz poleceń zmienia się na coś w stylu `(base) root@a1b2c3d4e5f6:/tmp#`, co wskazuje, że jesteś teraz wewnątrz kontenera.
Pliki danych są dostępne w katalogu `/data`.

#### 1.1.3. Uruchom polecenie indeksowania

[Dokumentacja Samtools](https://www.htslib.org/doc/samtools-index.html) podaje nam linię poleceń do uruchomienia w celu zindeksowania pliku BAM.

Musimy tylko podać plik wejściowy; narzędzie automatycznie wygeneruje nazwę dla wyjścia, dodając `.bai` do nazwy pliku wejściowego.

```bash
samtools index /data/bam/reads_mother.bam
```

??? abstract "Zawartość katalogu"

    ```console
    data/bam/
    ├── reads_father.bam
    ├── reads_mother.bam
    ├── reads_mother.bam.bai
    └── reads_son.bam
    ```

Teraz powinieneś/powinnaś zobaczyć plik o nazwie `reads_mother.bam.bai` w tym samym katalogu co oryginalny plik wejściowy BAM.

#### 1.1.4. Wyjdź z kontenera Samtools

Aby wyjść z kontenera, wpisz `exit`.

```bash
exit
```

Twój wiersz poleceń powinien teraz wrócić do tego, co było przed uruchomieniem kontenera.

### 1.2. Wykrywanie wariantów za pomocą GATK HaplotypeCaller

Pobierzemy kontener GATK, uruchomimy go interaktywnie i wykonamy polecenie `gatk HaplotypeCaller` na pliku BAM, który właśnie zindeksowaliśmy.

#### 1.2.1. Pobierz kontener GATK

Uruchom polecenie `docker pull`, aby pobrać obraz kontenera GATK:

```bash
docker pull community.wave.seqera.io/library/gatk4:4.5.0.0--730ee8817e436867
```

??? success "Wyjście polecenia"

    Niektóre warstwy pokazują `Already exists`, ponieważ są współdzielone z obrazem kontenera Samtools, który pobraliśmy wcześniej.

    ```console
    4.5.0.0--730ee8817e436867: Pulling from library/gatk4
    6360b3717211: Already exists
    2ec3f7ad9b3c: Already exists
    7716ca300600: Already exists
    4f4fb700ef54: Already exists
    8c61d418774c: Already exists
    03dae77ff45c: Already exists
    aab7f787139d: Already exists
    4f4fb700ef54: Already exists
    837d55536720: Already exists
    897362c12ca7: Already exists
    3893cbe24e91: Already exists
    d1b61e94977b: Already exists
    e5c558f54708: Pull complete
    087cce32d294: Pull complete
    Digest: sha256:e33413b9100f834fcc62fd5bc9edc1e881e820aafa606e09301eac2303d8724b
    Status: Downloaded newer image for community.wave.seqera.io/library/gatk4:4.5.0.0--730ee8817e436867
    community.wave.seqera.io/library/gatk4:4.5.0.0--730ee8817e436867
    ```

To powinno być szybsze niż pierwsze pobieranie, ponieważ oba obrazy kontenerów dzielą większość swoich warstw.

#### 1.2.2. Uruchom kontener GATK interaktywnie

Uruchom kontener GATK interaktywnie z zamontowanym katalogiem danych, tak jak zrobiliśmy to dla Samtools.

```bash
docker run -it -v ./data:/data community.wave.seqera.io/library/gatk4:4.5.0.0--730ee8817e436867
```

Twój wiersz poleceń zmienia się, wskazując, że jesteś teraz wewnątrz kontenera GATK.

#### 1.2.3. Uruchom polecenie wykrywania wariantów

[Dokumentacja GATK](https://gatk.broadinstitute.org/hc/en-us/articles/21905025322523-HaplotypeCaller) podaje nam linię poleceń do uruchomienia w celu wykonania wykrywania wariantów na pliku BAM.

Musimy podać plik wejściowy BAM (`-I`), a także genom referencyjny (`-R`), nazwę dla pliku wyjściowego (`-O`) oraz listę interwałów genomowych do analizy (`-L`).

Nie musimy jednak określać ścieżki do pliku indeksu; narzędzie automatycznie go wyszuka w tym samym katalogu, na podstawie ustalonej konwencji nazewnictwa i współlokalizacji.
To samo dotyczy plików pomocniczych genomu referencyjnego (pliki indeksu i słownika sekwencji, `*.fai` i `*.dict`).

```bash
gatk HaplotypeCaller \
        -R /data/ref/ref.fasta \
        -I /data/bam/reads_mother.bam \
        -O reads_mother.vcf \
        -L /data/ref/intervals.bed
```

??? success "Wyjście polecenia"

    Narzędzie generuje szczegółowy log wyjściowy. Podświetlone linie potwierdzają pomyślne zakończenie.

    ```console hl_lines="37 51 56 57"
    Using GATK jar /opt/conda/share/gatk4-4.5.0.0-0/gatk-package-4.5.0.0-local.jar
    Running:
        java -Dsamjdk.use_async_io_read_samtools=false -Dsamjdk.use_async_io_write_samtools=true -Dsamjdk.use_async_io_write_tribble=false -Dsamjdk.compression_level=2 -jar /opt/conda/share/gatk4-4.5.0.0-0/gatk-package-4.5.0.0-local.jar HaplotypeCaller -R /data/ref/ref.fasta -I /data/bam/reads_mother.bam -O reads_mother.vcf -L /data/ref/intervals.bed
    00:27:50.687 INFO  NativeLibraryLoader - Loading libgkl_compression.so from jar:file:/opt/conda/share/gatk4-4.5.0.0-0/gatk-package-4.5.0.0-local.jar!/com/intel/gkl/native/libgkl_compression.so
    00:27:50.854 INFO  HaplotypeCaller - ------------------------------------------------------------
    00:27:50.858 INFO  HaplotypeCaller - The Genome Analysis Toolkit (GATK) v4.5.0.0
    00:27:50.858 INFO  HaplotypeCaller - For support and documentation go to https://software.broadinstitute.org/gatk/
    00:27:50.858 INFO  HaplotypeCaller - Executing as root@a1fe8ff42d07 on Linux v6.10.14-linuxkit amd64
    00:27:50.858 INFO  HaplotypeCaller - Java runtime: OpenJDK 64-Bit Server VM v17.0.11-internal+0-adhoc..src
    00:27:50.859 INFO  HaplotypeCaller - Start Date/Time: February 8, 2026 at 12:27:50 AM GMT
    00:27:50.859 INFO  HaplotypeCaller - ------------------------------------------------------------
    00:27:50.859 INFO  HaplotypeCaller - ------------------------------------------------------------
    00:27:50.861 INFO  HaplotypeCaller - HTSJDK Version: 4.1.0
    00:27:50.861 INFO  HaplotypeCaller - Picard Version: 3.1.1
    00:27:50.861 INFO  HaplotypeCaller - Built for Spark Version: 3.5.0
    00:27:50.862 INFO  HaplotypeCaller - HTSJDK Defaults.COMPRESSION_LEVEL : 2
    00:27:50.862 INFO  HaplotypeCaller - HTSJDK Defaults.USE_ASYNC_IO_READ_FOR_SAMTOOLS : false
    00:27:50.862 INFO  HaplotypeCaller - HTSJDK Defaults.USE_ASYNC_IO_WRITE_FOR_SAMTOOLS : true
    00:27:50.863 INFO  HaplotypeCaller - HTSJDK Defaults.USE_ASYNC_IO_WRITE_FOR_TRIBBLE : false
    00:27:50.864 INFO  HaplotypeCaller - Deflater: IntelDeflater
    00:27:50.864 INFO  HaplotypeCaller - Inflater: IntelInflater
    00:27:50.864 INFO  HaplotypeCaller - GCS max retries/reopens: 20
    00:27:50.864 INFO  HaplotypeCaller - Requester pays: disabled
    00:27:50.865 INFO  HaplotypeCaller - Initializing engine
    00:27:50.991 INFO  FeatureManager - Using codec BEDCodec to read file file:///data/ref/intervals.bed
    00:27:51.016 INFO  IntervalArgumentCollection - Processing 6369 bp from intervals
    00:27:51.029 INFO  HaplotypeCaller - Done initializing engine
    00:27:51.040 INFO  NativeLibraryLoader - Loading libgkl_utils.so from jar:file:/opt/conda/share/gatk4-4.5.0.0-0/gatk-package-4.5.0.0-local.jar!/com/intel/gkl/native/libgkl_utils.so
    00:27:51.042 INFO  NativeLibraryLoader - Loading libgkl_smithwaterman.so from jar:file:/opt/conda/share/gatk4-4.5.0.0-0/gatk-package-4.5.0.0-local.jar!/com/intel/gkl/native/libgkl_smithwaterman.so
    00:27:51.042 INFO  SmithWatermanAligner - Using AVX accelerated SmithWaterman implementation
    00:27:51.046 INFO  HaplotypeCallerEngine - Disabling physical phasing, which is supported only for reference-model confidence output
    00:27:51.063 INFO  NativeLibraryLoader - Loading libgkl_pairhmm_omp.so from jar:file:/opt/conda/share/gatk4-4.5.0.0-0/gatk-package-4.5.0.0-local.jar!/com/intel/gkl/native/libgkl_pairhmm_omp.so
    00:27:51.085 INFO  IntelPairHmm - Flush-to-zero (FTZ) is enabled when running PairHMM
    00:27:51.086 INFO  IntelPairHmm - Available threads: 10
    00:27:51.086 INFO  IntelPairHmm - Requested threads: 4
    00:27:51.086 INFO  PairHMM - Using the OpenMP multi-threaded AVX-accelerated native PairHMM implementation
    00:27:51.128 INFO  ProgressMeter - Starting traversal
    00:27:51.136 INFO  ProgressMeter -        Current Locus  Elapsed Minutes     Regions Processed   Regions/Minute
    00:27:51.882 WARN  InbreedingCoeff - InbreedingCoeff will not be calculated at position 20_10037292_10066351:3480 and possibly subsequent; at least 10 samples must have called genotypes
    00:27:52.969 INFO  HaplotypeCaller - 7 read(s) filtered by: MappingQualityReadFilter
    0 read(s) filtered by: MappingQualityAvailableReadFilter
    0 read(s) filtered by: MappedReadFilter
    0 read(s) filtered by: NotSecondaryAlignmentReadFilter
    0 read(s) filtered by: NotDuplicateReadFilter
    0 read(s) filtered by: PassesVendorQualityCheckReadFilter
    0 read(s) filtered by: NonZeroReferenceLengthAlignmentReadFilter
    0 read(s) filtered by: GoodCigarReadFilter
    0 read(s) filtered by: WellformedReadFilter
    7 total reads filtered out of 1867 reads processed
    00:27:52.971 INFO  ProgressMeter - 20_10037292_10066351:13499              0.0                    35           1145.7
    00:27:52.971 INFO  ProgressMeter - Traversal complete. Processed 35 total regions in 0.0 minutes.
    00:27:52.976 INFO  VectorLoglessPairHMM - Time spent in setup for JNI call : 0.003346916
    00:27:52.976 INFO  PairHMM - Total compute time in PairHMM computeLogLikelihoods() : 0.045731709
    00:27:52.977 INFO  SmithWatermanAligner - Total compute time in native Smith-Waterman : 0.02 sec
    00:27:52.981 INFO  HaplotypeCaller - Shutting down engine
    [February 8, 2026 at 12:27:52 AM GMT] org.broadinstitute.hellbender.tools.walkers.haplotypecaller.HaplotypeCaller done. Elapsed time: 0.04 minutes.
    Runtime.totalMemory()=203423744
    ```

Plik wyjściowy `reads_mother.vcf` jest tworzony wewnątrz Twojego katalogu roboczego w kontenerze, więc nie zobaczysz go w eksploratorze plików VS Code, chyba że zmienisz ścieżkę pliku wyjściowego.
Jest to jednak mały plik testowy, więc możesz użyć `cat`, aby go otworzyć i zobaczyć zawartość.
Jeśli przewiniesz aż do początku pliku, znajdziesz nagłówek złożony z wielu linii metadanych, a następnie listę wykrytych wariantów, jeden w każdej linii.

??? abstract "Zawartość pliku"

    ```console title="reads_mother.vcf" linenums="26"
    #CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	reads_mother
    20_10037292_10066351	3480	.	C	CT	503.03	.	AC=2;AF=1.00;AN=2;DP=23;ExcessHet=0.0000;FS=0.000;MLEAC=2;MLEAF=1.00;MQ=60.00;QD=27.95;SOR=1.179	GT:AD:DP:GQ:PL	1/1:0,18:18:54:517,54,0
    20_10037292_10066351	3520	.	AT	A	609.03	.	AC=2;AF=1.00;AN=2;DP=18;ExcessHet=0.0000;FS=0.000;MLEAC=2;MLEAF=1.00;MQ=60.00;QD=33.83;SOR=0.693	GT:AD:DP:GQ:PL	1/1:0,18:18:54:623,54,0
    20_10037292_10066351	3529	.	T	A	155.64	.	AC=1;AF=0.500;AN=2;BaseQRankSum=-0.544;DP=21;ExcessHet=0.0000;FS=1.871;MLEAC=1;MLEAF=0.500;MQ=60.00;MQRankSum=0.000;QD=7.78;ReadPosRankSum=-1.158;SOR=1.034	GT:AD:DP:GQ:PL	0/1:12,8:20:99:163,0,328
    ```

Każda linia opisuje możliwy wariant zidentyfikowany w danych sekwencjonowania próbki. Aby uzyskać wskazówki dotyczące interpretacji formatu VCF, zobacz [ten pomocny artykuł](https://www.ebi.ac.uk/training/online/courses/human-genetic-variation-introduction/variant-identification-and-analysis/understanding-vcf-format/).

Plikowi wyjściowemu VCF towarzyszy plik indeksu o nazwie `reads_mother.vcf.idx`, który został automatycznie utworzony przez GATK.
Ma tę samą funkcję co plik indeksu BAM, aby umożliwić narzędziom wyszukiwanie i pobieranie podzbiorów danych bez ładowania całego pliku.

#### 1.2.4. Wyjdź z kontenera GATK

Aby wyjść z kontenera, wpisz `exit`.

```bash
exit
```

Twój wiersz poleceń powinien wrócić do normy.
To kończy test wykrywania wariantów per-próbka.

---

## 2. Łączne wykrywanie w kohorcie

Podejście wykrywania wariantów, którego właśnie użyliśmy, generuje wykrycia wariantów per-próbka.
Jest to w porządku do przeglądania wariantów z każdej próbki osobno, ale daje ograniczone informacje.
Często ciekawsze jest spojrzenie na to, jak wykrycia wariantów różnią się w wielu próbkach.
GATK oferuje alternatywną metodę zwaną łącznym wykrywaniem wariantów (joint variant calling) w tym celu.

Łączne wykrywanie wariantów obejmuje generowanie specjalnego rodzaju wyjścia wariantów zwanego GVCF (dla Genomic VCF) dla każdej próbki, następnie łączenie danych GVCF ze wszystkich próbek i uruchamianie analizy statystycznej „łącznego genotypowania".

![Analiza łączna](img/joint-calling.png)

To, co jest specjalne w GVCF próbki, to fakt, że zawiera rekordy podsumowujące statystyki danych sekwencjonowania dla wszystkich pozycji w docelowym obszarze genomu, nie tylko pozycji, w których program znalazł dowody zmienności.
Jest to kluczowe dla obliczeń łącznego genotypowania ([więcej informacji](https://gatk.broadinstitute.org/hc/en-us/articles/360035890431-The-logic-of-joint-calling-for-germline-short-variants)).

GVCF jest tworzony przez GATK HaplotypeCaller, to samo narzędzie, którego właśnie przetestowaliśmy, z dodatkowym parametrem (`-ERC GVCF`).
Łączenie plików GVCF odbywa się za pomocą GATK GenomicsDBImport, który łączy wykrycia per-próbka w magazyn danych (analogiczny do bazy danych).
Faktyczna analiza „łącznego genotypowania" jest następnie wykonywana za pomocą GATK GenotypeGVCFs.

Tutaj testujemy polecenia potrzebne do generowania GVCF i uruchomienia łącznego genotypowania.
To polecenia, które opakujemy w workflow Nextflow w części 3 tego kursu.

1. Wygeneruj plik indeksu dla każdego wejściowego pliku BAM za pomocą Samtools
2. Uruchom GATK HaplotypeCaller na każdym wejściowym pliku BAM, aby wygenerować GVCF wykrytych wariantów genomowych per-próbka
3. Zbierz wszystkie GVCF i połącz je w magazyn danych GenomicsDB
4. Uruchom łączne genotypowanie na połączonym magazynie danych GVCF, aby wytworzyć VCF na poziomie kohorty

<figure class="excalidraw">
--8<-- "docs/en/docs/nf4_science/genomics/img/hello-gatk-2.svg"
</figure>

Teraz musimy przetestować wszystkie te polecenia, zaczynając od zindeksowania wszystkich trzech plików BAM.

### 2.1. Indeksowanie plików BAM dla wszystkich trzech próbek

W pierwszej sekcji powyżej zindeksowaliśmy tylko jeden plik BAM.
Teraz musimy zindeksować wszystkie trzy próbki, aby GATK HaplotypeCaller mógł je przetworzyć.

#### 2.1.1. Uruchom kontener Samtools interaktywnie

Już pobraliśmy obraz kontenera Samtools, więc możemy go uruchomić bezpośrednio:

```bash
docker run -it -v ./data:/data community.wave.seqera.io/library/samtools:1.20--b5dfbd93de237464
```

Twój wiersz poleceń zmienia się, wskazując, że jesteś wewnątrz kontenera, z zamontowanym katalogiem danych jak wcześniej.

#### 2.1.2. Uruchom polecenie indeksowania na wszystkich trzech próbkach

Uruchom polecenie indeksowania na każdym z trzech plików BAM:

```bash
samtools index /data/bam/reads_mother.bam
samtools index /data/bam/reads_father.bam
samtools index /data/bam/reads_son.bam
```

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

To powinno wytworzyć pliki indeksów w tym samym katalogu co odpowiadające im pliki BAM.

#### 2.1.3. Wyjdź z kontenera Samtools

Aby wyjść z kontenera, wpisz `exit`.

```bash
exit
```

Twój wiersz poleceń powinien wrócić do normy.

### 2.2. Generowanie plików GVCF dla wszystkich trzech próbek

Aby uruchomić krok łącznego genotypowania, potrzebujemy plików GVCF dla wszystkich trzech próbek.

#### 2.2.1. Uruchom kontener GATK interaktywnie

Już pobraliśmy obraz kontenera GATK wcześniej, więc możemy go uruchomić bezpośrednio:

```bash
docker run -it -v ./data:/data community.wave.seqera.io/library/gatk4:4.5.0.0--730ee8817e436867
```

Twój wiersz poleceń zmienia się, wskazując, że jesteś wewnątrz kontenera GATK.

#### 2.2.2. Uruchom polecenie wykrywania wariantów z opcją GVCF

Aby wytworzyć genomowy VCF (GVCF), dodajemy opcję `-ERC GVCF` do podstawowego polecenia, która włącza tryb GVCF narzędzia HaplotypeCaller.

Zmieniamy też rozszerzenie pliku dla pliku wyjściowego z `.vcf` na `.g.vcf`.
Technicznie nie jest to wymóg, ale jest to silnie zalecana konwencja.

```bash
gatk HaplotypeCaller \
        -R /data/ref/ref.fasta \
        -I /data/bam/reads_mother.bam \
        -O reads_mother.g.vcf \
        -L /data/ref/intervals.bed \
        -ERC GVCF
```

??? success "Wyjście polecenia"

    ```console hl_lines="39 53 58 59"
    Using GATK jar /opt/conda/share/gatk4-4.5.0.0-0/gatk-package-4.5.0.0-local.jar
    Running:
        java -Dsamjdk.use_async_io_read_samtools=false -Dsamjdk.use_async_io_write_samtools=true -Dsamjdk.use_async_io_write_tribble=false -Dsamjdk.compression_level=2 -jar /opt/conda/share/gatk4-4.5.0.0-0/gatk-package-4.5.0.0-local.jar HaplotypeCaller -R /data/ref/ref.fasta -I /data/bam/reads_mother.bam -O reads_mother.g.vcf -L /data/ref/intervals.bed -ERC GVCF
    00:28:03.593 INFO  NativeLibraryLoader - Loading libgkl_compression.so from jar:file:/opt/conda/share/gatk4-4.5.0.0-0/gatk-package-4.5.0.0-local.jar!/com/intel/gkl/native/libgkl_compression.so
    00:28:03.765 INFO  HaplotypeCaller - ------------------------------------------------------------
    00:28:03.768 INFO  HaplotypeCaller - The Genome Analysis Toolkit (GATK) v4.5.0.0
    00:28:03.768 INFO  HaplotypeCaller - For support and documentation go to https://software.broadinstitute.org/gatk/
    00:28:03.768 INFO  HaplotypeCaller - Executing as root@8515e5a0598e on Linux v6.10.14-linuxkit amd64
    00:28:03.768 INFO  HaplotypeCaller - Java runtime: OpenJDK 64-Bit Server VM v17.0.11-internal+0-adhoc..src
    00:28:03.769 INFO  HaplotypeCaller - Start Date/Time: February 8, 2026 at 12:28:03 AM GMT
    00:28:03.769 INFO  HaplotypeCaller - ------------------------------------------------------------
    00:28:03.770 INFO  HaplotypeCaller - ------------------------------------------------------------
    00:28:03.772 INFO  HaplotypeCaller - HTSJDK Version: 4.1.0
    00:28:03.773 INFO  HaplotypeCaller - Picard Version: 3.1.1
    00:28:03.773 INFO  HaplotypeCaller - Built for Spark Version: 3.5.0
    00:28:03.773 INFO  HaplotypeCaller - HTSJDK Defaults.COMPRESSION_LEVEL : 2
    00:28:03.773 INFO  HaplotypeCaller - HTSJDK Defaults.USE_ASYNC_IO_READ_FOR_SAMTOOLS : false
    00:28:03.773 INFO  HaplotypeCaller - HTSJDK Defaults.USE_ASYNC_IO_WRITE_FOR_SAMTOOLS : true
    00:28:03.774 INFO  HaplotypeCaller - HTSJDK Defaults.USE_ASYNC_IO_WRITE_FOR_TRIBBLE : false
    00:28:03.774 INFO  HaplotypeCaller - Deflater: IntelDeflater
    00:28:03.774 INFO  HaplotypeCaller - Inflater: IntelInflater
    00:28:03.775 INFO  HaplotypeCaller - GCS max retries/reopens: 20
    00:28:03.775 INFO  HaplotypeCaller - Requester pays: disabled
    00:28:03.776 INFO  HaplotypeCaller - Initializing engine
    00:28:03.896 INFO  FeatureManager - Using codec BEDCodec to read file file:///data/ref/intervals.bed
    00:28:03.919 INFO  IntervalArgumentCollection - Processing 6369 bp from intervals
    00:28:03.934 INFO  HaplotypeCaller - Done initializing engine
    00:28:03.935 INFO  HaplotypeCallerEngine - Tool is in reference confidence mode and the annotation, the following changes will be made to any specified annotations: 'StrandBiasBySample' will be enabled. 'ChromosomeCounts', 'FisherStrand', 'StrandOddsRatio' and 'QualByDepth' annotations have been disabled
    00:28:03.943 INFO  NativeLibraryLoader - Loading libgkl_utils.so from jar:file:/opt/conda/share/gatk4-4.5.0.0-0/gatk-package-4.5.0.0-local.jar!/com/intel/gkl/native/libgkl_utils.so
    00:28:03.945 INFO  NativeLibraryLoader - Loading libgkl_smithwaterman.so from jar:file:/opt/conda/share/gatk4-4.5.0.0-0/gatk-package-4.5.0.0-local.jar!/com/intel/gkl/native/libgkl_smithwaterman.so
    00:28:03.946 INFO  SmithWatermanAligner - Using AVX accelerated SmithWaterman implementation
    00:28:03.955 INFO  HaplotypeCallerEngine - Standard Emitting and Calling confidence set to -0.0 for reference-model confidence output
    00:28:03.956 INFO  HaplotypeCallerEngine - All sites annotated with PLs forced to true for reference-model confidence output
    00:28:03.972 INFO  NativeLibraryLoader - Loading libgkl_pairhmm_omp.so from jar:file:/opt/conda/share/gatk4-4.5.0.0-0/gatk-package-4.5.0.0-local.jar!/com/intel/gkl/native/libgkl_pairhmm_omp.so
    00:28:03.993 INFO  IntelPairHmm - Flush-to-zero (FTZ) is enabled when running PairHMM
    00:28:03.994 INFO  IntelPairHmm - Available threads: 10
    00:28:03.994 INFO  IntelPairHmm - Requested threads: 4
    00:28:03.994 INFO  PairHMM - Using the OpenMP multi-threaded AVX-accelerated native PairHMM implementation
    00:28:04.044 INFO  ProgressMeter - Starting traversal
    00:28:04.070 INFO  ProgressMeter -        Current Locus  Elapsed Minutes     Regions Processed   Regions/Minute
    00:28:04.874 WARN  InbreedingCoeff - InbreedingCoeff will not be calculated at position 20_10037292_10066351:3480 and possibly subsequent; at least 10 samples must have called genotypes
    00:28:06.535 INFO  HaplotypeCaller - 7 read(s) filtered by: MappingQualityReadFilter
    0 read(s) filtered by: MappingQualityAvailableReadFilter
    0 read(s) filtered by: MappedReadFilter
    0 read(s) filtered by: NotSecondaryAlignmentReadFilter
    0 read(s) filtered by: NotDuplicateReadFilter
    0 read(s) filtered by: PassesVendorQualityCheckReadFilter
    0 read(s) filtered by: NonZeroReferenceLengthAlignmentReadFilter
    0 read(s) filtered by: GoodCigarReadFilter
    0 read(s) filtered by: WellformedReadFilter
    7 total reads filtered out of 1867 reads processed
    00:28:06.537 INFO  ProgressMeter - 20_10037292_10066351:13499              0.0                    35            851.6
    00:28:06.538 INFO  ProgressMeter - Traversal complete. Processed 35 total regions in 0.0 minutes.
    00:28:06.543 INFO  VectorLoglessPairHMM - Time spent in setup for JNI call : 0.003648749
    00:28:06.544 INFO  PairHMM - Total compute time in PairHMM computeLogLikelihoods() : 0.031498916
    00:28:06.544 INFO  SmithWatermanAligner - Total compute time in native Smith-Waterman : 0.02 sec
    00:28:06.547 INFO  HaplotypeCaller - Shutting down engine
    [February 8, 2026 at 12:28:06 AM GMT] org.broadinstitute.hellbender.tools.walkers.haplotypecaller.HaplotypeCaller done. Elapsed time: 0.05 minutes.
    Runtime.totalMemory()=281018368
    ```

To tworzy plik wyjściowy GVCF `reads_mother.g.vcf` w bieżącym katalogu roboczym w kontenerze.

Jeśli użyjesz `cat`, aby zobaczyć zawartość, zobaczysz, że jest on znacznie dłuższy niż równoważny VCF, który wygenerowaliśmy w sekcji 1. Nie możesz nawet przewinąć do początku pliku, a większość linii wygląda zupełnie inaczej niż to, co widzieliśmy w VCF.

??? abstract "Zawartość pliku"

    ```console title="reads_mother.g.vcf" linenums="1674"
    20_10037292_10066351    14714   .       T       <NON_REF>       .       .       END=14718       GT:DP:GQ:MIN_DP:PL       0/0:37:99:37:0,99,1192
    20_10037292_10066351    14719   .       T       <NON_REF>       .       .       END=14719       GT:DP:GQ:MIN_DP:PL       0/0:36:82:36:0,82,1087
    20_10037292_10066351    14720   .       T       <NON_REF>       .       .       END=14737       GT:DP:GQ:MIN_DP:PL       0/0:42:99:37:0,100,1160
    ```

Reprezentują regiony bez wariantów, w których narzędzie wykrywające warianty nie znalazło dowodów zmienności, więc przechwycono pewne statystyki opisujące jego poziom pewności w braku zmienności.
Umożliwia to rozróżnienie między dwoma bardzo różnymi przypadkami: (1) są dane dobrej jakości pokazujące, że próbka jest homozygotyczna względem referencji, i (2) nie ma wystarczającej ilości dobrych danych, aby w ogóle dokonać ustalenia.

W GVCF zazwyczaj jest wiele takich linii bez wariantów, z mniejszą liczbą rekordów wariantów rozsypanych wśród nich.
Spróbuj uruchomić `head -176` na pliku GVCF, aby załadować tylko pierwsze 176 linii pliku i znaleźć faktyczne wykrycie wariantu.

??? abstract "Zawartość pliku"

    ```console title="reads_mother.g.vcf" linenums="174"
    20_10037292_10066351    3479    .       T       <NON_REF>       .       .       END=3479        GT:DP:GQ:MIN_DP:PL       0/0:34:36:34:0,36,906
    20_10037292_10066351    3480    .       C       CT,<NON_REF>    503.03  .       DP=23;ExcessHet=0.0000;MLEAC=2,0;MLEAF=1.00,0.00;RAW_MQandDP=82800,23    GT:AD:DP:GQ:PL:SB       1/1:0,18,0:18:54:517,54,0,517,54,517:0,0,7,11
    20_10037292_10066351    3481    .       T       <NON_REF>       .       .       END=3481        GT:DP:GQ:MIN_DP:PL       0/0:21:51:21:0,51,765
    ```

Druga linia pokazuje pierwszy rekord wariantu w pliku, który odpowiada pierwszemu wariantowi w pliku VCF, na który patrzyliśmy wcześniej.

Tak jak oryginalny VCF, plikowi wyjściowemu GVCF również towarzyszy plik indeksu o nazwie `reads_mother.g.vcf.idx`.

#### 2.2.3. Powtórz proces dla pozostałych dwóch próbek

Wygeneruj pliki GVCF dla pozostałych dwóch próbek, uruchamiając poniższe polecenia, jedno po drugim.

```bash
gatk HaplotypeCaller \
        -R /data/ref/ref.fasta \
        -I /data/bam/reads_father.bam \
        -O reads_father.g.vcf \
        -L /data/ref/intervals.bed \
        -ERC GVCF
```

```bash
gatk HaplotypeCaller \
        -R /data/ref/ref.fasta \
        -I /data/bam/reads_son.bam \
        -O reads_son.g.vcf \
        -L /data/ref/intervals.bed \
        -ERC GVCF
```

Po zakończeniu powinieneś/powinnaś mieć trzy pliki kończące się na `.g.vcf` w swoim bieżącym katalogu (jeden na próbkę) oraz ich odpowiednie pliki indeksów kończące się na `.g.vcf.idx`.

Ale nie wychodź z kontenera!
Użyjemy tego samego kontenera w następnym kroku.

### 2.3. Uruchomienie łącznego genotypowania

Teraz, gdy mamy wszystkie pliki GVCF, możemy wypróbować podejście łącznego genotypowania do generowania wykrytych wariantów dla kohorty próbek.
Jest to metoda dwuetapowa, która składa się z połączenia danych ze wszystkich plików GVCF w magazyn danych, a następnie uruchomienia właściwej analizy łącznego genotypowania w celu wygenerowania końcowego VCF łącznie wykrytych wariantów.

#### 2.3.1. Połącz wszystkie pliki GVCF per-próbka

Ten pierwszy krok używa innego narzędzia GATK, zwanego GenomicsDBImport, do połączenia danych ze wszystkich plików GVCF w magazyn danych GenomicsDB.

```bash
gatk GenomicsDBImport \
    -V reads_mother.g.vcf \
    -V reads_father.g.vcf \
    -V reads_son.g.vcf \
    -L /data/ref/intervals.bed \
    --genomicsdb-workspace-path family_trio_gdb
```

??? success "Wyjście polecenia"

    ```console hl_lines="33 36 37 39 40"
    Using GATK jar /opt/conda/share/gatk4-4.5.0.0-0/gatk-package-4.5.0.0-local.jar
    Running:
        java -Dsamjdk.use_async_io_read_samtools=false -Dsamjdk.use_async_io_write_samtools=true -Dsamjdk.use_async_io_write_tribble=false -Dsamjdk.compression_level=2 -jar /opt/conda/share/gatk4-4.5.0.0-0/gatk-package-4.5.0.0-local.jar GenomicsDBImport -V reads_mother.g.vcf -V reads_father.g.vcf -V reads_son.g.vcf -L /data/ref/intervals.bed --genomicsdb-workspace-path family_trio_gdb
    00:28:20.772 INFO  NativeLibraryLoader - Loading libgkl_compression.so from jar:file:/opt/conda/share/gatk4-4.5.0.0-0/gatk-package-4.5.0.0-local.jar!/com/intel/gkl/native/libgkl_compression.so
    00:28:20.914 INFO  GenomicsDBImport - ------------------------------------------------------------
    00:28:20.917 INFO  GenomicsDBImport - The Genome Analysis Toolkit (GATK) v4.5.0.0
    00:28:20.917 INFO  GenomicsDBImport - For support and documentation go to https://software.broadinstitute.org/gatk/
    00:28:20.917 INFO  GenomicsDBImport - Executing as root@8515e5a0598e on Linux v6.10.14-linuxkit amd64
    00:28:20.917 INFO  GenomicsDBImport - Java runtime: OpenJDK 64-Bit Server VM v17.0.11-internal+0-adhoc..src
    00:28:20.918 INFO  GenomicsDBImport - Start Date/Time: February 8, 2026 at 12:28:20 AM GMT
    00:28:20.918 INFO  GenomicsDBImport - ------------------------------------------------------------
    00:28:20.918 INFO  GenomicsDBImport - ------------------------------------------------------------
    00:28:20.920 INFO  GenomicsDBImport - HTSJDK Version: 4.1.0
    00:28:20.921 INFO  GenomicsDBImport - Picard Version: 3.1.1
    00:28:20.921 INFO  GenomicsDBImport - Built for Spark Version: 3.5.0
    00:28:20.922 INFO  GenomicsDBImport - HTSJDK Defaults.COMPRESSION_LEVEL : 2
    00:28:20.923 INFO  GenomicsDBImport - HTSJDK Defaults.USE_ASYNC_IO_READ_FOR_SAMTOOLS : false
    00:28:20.923 INFO  GenomicsDBImport - HTSJDK Defaults.USE_ASYNC_IO_WRITE_FOR_SAMTOOLS : true
    00:28:20.923 INFO  GenomicsDBImport - HTSJDK Defaults.USE_ASYNC_IO_WRITE_FOR_TRIBBLE : false
    00:28:20.923 INFO  GenomicsDBImport - Deflater: IntelDeflater
    00:28:20.924 INFO  GenomicsDBImport - Inflater: IntelInflater
    00:28:20.924 INFO  GenomicsDBImport - GCS max retries/reopens: 20
    00:28:20.924 INFO  GenomicsDBImport - Requester pays: disabled
    00:28:20.925 INFO  GenomicsDBImport - Initializing engine
    00:28:21.144 INFO  FeatureManager - Using codec BEDCodec to read file file:///data/ref/intervals.bed
    00:28:21.152 INFO  IntervalArgumentCollection - Processing 6369 bp from intervals
    00:28:21.157 INFO  GenomicsDBImport - Done initializing engine
    00:28:21.287 INFO  GenomicsDBLibLoader - GenomicsDB native library version : 1.5.1-84e800e
    00:28:21.290 INFO  GenomicsDBImport - Vid Map JSON file will be written to /tmp/family_trio_gdb/vidmap.json
    00:28:21.290 INFO  GenomicsDBImport - Callset Map JSON file will be written to /tmp/family_trio_gdb/callset.json
    00:28:21.291 INFO  GenomicsDBImport - Complete VCF Header will be written to /tmp/family_trio_gdb/vcfheader.vcf
    00:28:21.291 INFO  GenomicsDBImport - Importing to workspace - /tmp/family_trio_gdb
    00:28:21.453 INFO  GenomicsDBImport - Importing batch 1 with 3 samples
    00:28:21.757 INFO  GenomicsDBImport - Importing batch 1 with 3 samples
    00:28:21.859 INFO  GenomicsDBImport - Importing batch 1 with 3 samples
    00:28:21.979 INFO  GenomicsDBImport - Done importing batch 1/1
    00:28:21.988 INFO  GenomicsDBImport - Import completed!
    00:28:21.988 INFO  GenomicsDBImport - Shutting down engine
    [February 8, 2026 at 12:28:21 AM GMT] org.broadinstitute.hellbender.tools.genomicsdb.GenomicsDBImport done. Elapsed time: 0.02 minutes.
    Runtime.totalMemory()=305135616
    ```

Wynikiem tego kroku jest efektywnie katalog zawierający zestaw dalszych zagnieżdżonych katalogów przechowujących połączone dane wariantów w postaci wielu różnych plików.
Możesz go przeglądać, ale szybko zobaczysz, że ten format magazynu danych nie jest przeznaczony do bezpośredniego odczytu przez ludzi.

!!! note "Uwaga"

    GATK zawiera narzędzia, które umożliwiają przeglądanie i wyodrębnianie danych wykrytych wariantów z magazynu danych w razie potrzeby.

#### 2.3.2. Uruchom właściwą analizę łącznego genotypowania

Ten drugi krok używa jeszcze innego narzędzia GATK, zwanego GenotypeGVCFs, do ponownego obliczenia statystyk wariantów i indywidualnych genotypów w świetle danych dostępnych we wszystkich próbkach w kohorcie.

```bash
gatk GenotypeGVCFs \
    -R /data/ref/ref.fasta \
    -V gendb://family_trio_gdb \
    -O family_trio.vcf
```

??? success "Wyjście polecenia"

    ```console hl_lines="30 35 37 38"
    Using GATK jar /opt/conda/share/gatk4-4.5.0.0-0/gatk-package-4.5.0.0-local.jar
    Running:
        java -Dsamjdk.use_async_io_read_samtools=false -Dsamjdk.use_async_io_write_samtools=true -Dsamjdk.use_async_io_write_tribble=false -Dsamjdk.compression_level=2 -jar /opt/conda/share/gatk4-4.5.0.0-0/gatk-package-4.5.0.0-local.jar GenotypeGVCFs -R /data/ref/ref.fasta -V gendb://family_trio_gdb -O family_trio.vcf
    00:28:24.625 INFO  NativeLibraryLoader - Loading libgkl_compression.so from jar:file:/opt/conda/share/gatk4-4.5.0.0-0/gatk-package-4.5.0.0-local.jar!/com/intel/gkl/native/libgkl_compression.so
    00:28:24.798 INFO  GenotypeGVCFs - ------------------------------------------------------------
    00:28:24.801 INFO  GenotypeGVCFs - The Genome Analysis Toolkit (GATK) v4.5.0.0
    00:28:24.801 INFO  GenotypeGVCFs - For support and documentation go to https://software.broadinstitute.org/gatk/
    00:28:24.801 INFO  GenotypeGVCFs - Executing as root@8515e5a0598e on Linux v6.10.14-linuxkit amd64
    00:28:24.801 INFO  GenotypeGVCFs - Java runtime: OpenJDK 64-Bit Server VM v17.0.11-internal+0-adhoc..src
    00:28:24.802 INFO  GenotypeGVCFs - Start Date/Time: February 8, 2026 at 12:28:24 AM GMT
    00:28:24.802 INFO  GenotypeGVCFs - ------------------------------------------------------------
    00:28:24.802 INFO  GenotypeGVCFs - ------------------------------------------------------------
    00:28:24.804 INFO  GenotypeGVCFs - HTSJDK Version: 4.1.0
    00:28:24.804 INFO  GenotypeGVCFs - Picard Version: 3.1.1
    00:28:24.804 INFO  GenotypeGVCFs - Built for Spark Version: 3.5.0
    00:28:24.805 INFO  GenotypeGVCFs - HTSJDK Defaults.COMPRESSION_LEVEL : 2
    00:28:24.805 INFO  GenotypeGVCFs - HTSJDK Defaults.USE_ASYNC_IO_READ_FOR_SAMTOOLS : false
    00:28:24.806 INFO  GenotypeGVCFs - HTSJDK Defaults.USE_ASYNC_IO_WRITE_FOR_SAMTOOLS : true
    00:28:24.806 INFO  GenotypeGVCFs - HTSJDK Defaults.USE_ASYNC_IO_WRITE_FOR_TRIBBLE : false
    00:28:24.806 INFO  GenotypeGVCFs - Deflater: IntelDeflater
    00:28:24.806 INFO  GenotypeGVCFs - Inflater: IntelInflater
    00:28:24.807 INFO  GenotypeGVCFs - GCS max retries/reopens: 20
    00:28:24.807 INFO  GenotypeGVCFs - Requester pays: disabled
    00:28:24.808 INFO  GenotypeGVCFs - Initializing engine
    00:28:25.023 INFO  GenomicsDBLibLoader - GenomicsDB native library version : 1.5.1-84e800e
    00:28:25.081 INFO  NativeGenomicsDB - pid=162 tid=163 No valid combination operation found for INFO field InbreedingCoeff  - the field will NOT be part of INFO fields in the generated VCF records
    00:28:25.082 INFO  NativeGenomicsDB - pid=162 tid=163 No valid combination operation found for INFO field MLEAC  - the field will NOT be part of INFO fields in the generated VCF records
    00:28:25.082 INFO  NativeGenomicsDB - pid=162 tid=163 No valid combination operation found for INFO field MLEAF  - the field will NOT be part of INFO fields in the generated VCF records
    00:28:25.109 INFO  GenotypeGVCFs - Done initializing engine
    00:28:25.184 INFO  ProgressMeter - Starting traversal
    00:28:25.187 INFO  ProgressMeter -        Current Locus  Elapsed Minutes    Variants Processed  Variants/Minute
    00:28:25.446 WARN  InbreedingCoeff - InbreedingCoeff will not be calculated at position 20_10037292_10066351:3480 and possibly subsequent; at least 10 samples must have called genotypes
    GENOMICSDB_TIMER,GenomicsDB iterator next() timer,Wall-clock time(s),0.15034835899999904,Cpu time(s),0.1355218420000006
    00:28:26.189 INFO  ProgressMeter - 20_10037292_10066351:13953              0.0                  3390         202994.0
    00:28:26.190 INFO  ProgressMeter - Traversal complete. Processed 3390 total variants in 0.0 minutes.
    00:28:26.194 INFO  GenotypeGVCFs - Shutting down engine
    [February 8, 2026 at 12:28:26 AM GMT] org.broadinstitute.hellbender.tools.walkers.GenotypeGVCFs done. Elapsed time: 0.03 minutes.
    Runtime.totalMemory()=296747008
    ```

To tworzy plik wyjściowy VCF `family_trio.vcf` w bieżącym katalogu roboczym w kontenerze.
Jest to kolejny dość mały plik, więc możesz użyć `cat`, aby zobaczyć jego zawartość, i przewinąć w górę, aby znaleźć pierwsze kilka linii wariantów.

??? abstract "Zawartość pliku"

    ```console title="family_trio.vcf" linenums="40"
    #CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO    FORMAT  reads_father    reads_mother    reads_son
    20_10037292_10066351    3480    .       C       CT      1625.89 .       AC=5;AF=0.833;AN=6;BaseQRankSum=0.220;DP=85;ExcessHet=0.0000;FS=2.476;MLEAC=5;MLEAF=0.833;MQ=60.00;MQRankSum=0.00;QD=21.68;ReadPosRankSum=-1.147e+00;SOR=0.487    GT:AD:DP:GQ:PL  0/1:15,16:31:99:367,0,375       1/1:0,18:18:54:517,54,0 1/1:0,26:26:78:756,78,0
    20_10037292_10066351    3520    .       AT      A       1678.89 .       AC=5;AF=0.833;AN=6;BaseQRankSum=1.03;DP=80;ExcessHet=0.0000;FS=2.290;MLEAC=5;MLEAF=0.833;MQ=60.00;MQRankSum=0.00;QD=22.39;ReadPosRankSum=0.701;SOR=0.730 GT:AD:DP:GQ:PL   0/1:18,13:31:99:296,0,424       1/1:0,18:18:54:623,54,0 1/1:0,26:26:78:774,78,0
    20_10037292_10066351    3529    .       T       A       154.29  .       AC=1;AF=0.167;AN=6;BaseQRankSum=-5.440e-01;DP=104;ExcessHet=0.0000;FS=1.871;MLEAC=1;MLEAF=0.167;MQ=60.00;MQRankSum=0.00;QD=7.71;ReadPosRankSum=-1.158e+00;SOR=1.034       GT:AD:DP:GQ:PL  0/0:44,0:44:99:0,112,1347       0/1:12,8:20:99:163,0,328        0/0:39,0:39:99:0,105,1194
    ```

Wygląda to podobnie do VCF, który wygenerowaliśmy wcześniej, z tą różnicą, że tym razem mamy informacje na poziomie genotypu dla wszystkich trzech próbek.
Ostatnie trzy kolumny w pliku to bloki genotypowe dla próbek, wymienione w porządku alfabetycznym.

Jeśli spojrzymy na genotypy wykryte dla naszego testowego trio rodzinnego dla pierwszego wariantu, widzimy, że ojciec jest heterozygotyczny-wariantowy (`0/1`), a matka i syn są obaj homozygotyczni-wariantowi (`1/1`).

To ostatecznie informacja, którą chcemy wydobyć ze zbioru danych!

#### 2.3.3. Wyjdź z kontenera GATK

Aby wyjść z kontenera, wpisz `exit`.

```bash
exit
```

Twój wiersz poleceń powinien wrócić do normy.
To kończy ręczne testowanie poleceń wykrywania wariantów.

---

### Podsumowanie

Wiesz, jak przetestować polecenia indeksowania Samtools i wykrywania wariantów GATK w ich odpowiednich kontenerach, w tym jak generować pliki GVCF i uruchamiać łączne genotypowanie na wielu próbkach.

### Co dalej?

Naucz się, jak opakować te same polecenia w workflow'y, które używają kontenerów do wykonania pracy.
