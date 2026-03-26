# 파트 1: 방법 개요

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } AI 지원 번역 - [자세히 알아보기 및 개선 제안](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

대량 RNAseq 데이터를 처리하고 분석하는 여러 가지 유효한 방법이 있습니다.
이 과정에서는 [Babraham Institute](https://www.bioinformatics.babraham.ac.uk/training/RNASeq_Course/Analysing%20RNA-Seq%20data%20Exercise.pdf)의 Simon Andrews 박사와 Laura Biggins 박사가 [여기](https://www.babraham.ac.uk/)에서 설명한 방법을 따릅니다.

우리의 목표는 다음 처리 단계를 구현하는 워크플로우를 개발하는 것입니다: 대량 RNAseq 샘플의 리드에 대한 초기 품질 관리 실행, 리드에서 어댑터 서열 트리밍, 참조 게놈에 리드 정렬, 그리고 포괄적인 품질 관리(QC) 보고서 생성.

<figure class="excalidraw">
--8<-- "docs/en/docs/nf4_science/rnaseq/img/rnaseq-method-03.svg"
</figure>

- **FASTQC:** FastQC를 사용하여 트리밍 전 리드 데이터에 대한 QC 수행
- **TRIM_GALORE:** Trim Galore를 사용하여 어댑터 서열 트리밍 및 트리밍 후 QC 수행 (Cutadapt 및 FastQC 번들)
- **HISAT2_ALIGN:** Hisat2를 사용하여 참조 게놈에 리드 정렬
- **MULTIQC:** MultiQC를 사용하여 포괄적인 QC 보고서 생성

### 방법

이러한 처리 단계를 두 단계로 나누어 적용하는 방법을 보여드리겠습니다.
먼저 하나의 샘플에 대해 QC, 트리밍 및 정렬 도구를 실행하는 **단일 샘플 처리**부터 시작합니다.
그런 다음 여러 샘플에 대해 동일한 도구를 실행하고 집계된 품질 관리 보고서를 생성하는 **다중 샘플 처리**로 확장합니다.

워크플로우 코드를 작성하기 전에 먼저 일부 테스트 데이터로 명령을 수동으로 시도해 봅니다.

### 데이터셋

다음 데이터 및 관련 리소스를 제공합니다:

- **RNAseq 데이터** (`reads/`): 6개 샘플의 FASTQ 파일, 파일 크기를 줄이기 위해 작은 영역으로 제한됨. 각 샘플은 paired-end 리드를 가지고 있지만(샘플당 2개 파일), 먼저 single-end 리드만 사용합니다.
- **참조 게놈** (`genome.fa`): 인간 염색체 20의 작은 영역 (hg19/b37에서 추출).
- **CSV 샘플시트** (`single-end.csv` 및 `paired-end.csv`): 예제 데이터 파일의 ID와 경로를 나열한 파일.

### 소프트웨어

관련된 네 가지 주요 도구는 품질 관리 메트릭 수집을 위한 [FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/), 어댑터 트리밍을 위한 [Trim Galore](https://www.bioinformatics.babraham.ac.uk/projects/trim_galore/) (Cutadapt 및 FastQC를 번들로 제공하여 트리밍 후 QC 수행), 참조 게놈에 대한 spliced 정렬을 위한 [HISAT2](http://daehwankimlab.github.io/hisat2/), 그리고 집계된 QC 보고서 생성을 위한 [MultiQC](https://multiqc.info/)입니다.

이러한 도구는 GitHub Codespaces 환경에 설치되어 있지 않으므로 Seqera Containers 서비스를 통해 검색한 컨테이너를 통해 사용할 것입니다 ([Hello Containers](../../hello_nextflow/05_hello_containers.md) 참조).

!!! tip "팁"

     `nf4-science/rnaseq` 디렉토리에 있는지 확인하십시오. `pwd`를 입력할 때 표시되는 경로의 마지막 부분은 `rnaseq`이어야 합니다.

---

## 1. 단일 샘플 처리

이 섹션에서는 단일 RNAseq 샘플을 처리하는 명령을 테스트합니다: 품질 관리, 어댑터 트리밍, 그리고 참조 게놈에 대한 정렬.
이것들이 이 과정의 파트 2에서 Nextflow 워크플로우로 적용할 명령입니다.

1. FastQC를 사용하여 FASTQ 파일에 대한 초기 QC 실행
2. Trim Galore를 사용하여 어댑터 서열 트리밍 및 트리밍 후 QC 실행
3. HISAT2를 사용하여 트리밍된 리드를 참조 게놈에 정렬

<figure class="excalidraw">
--8<-- "docs/en/docs/nf4_science/rnaseq/img/rnaseq-method-02.svg"
</figure>

먼저 하나의 샘플에 대해서만 이러한 명령을 테스트합니다.

### 1.1. QC 및 어댑터 트리밍

먼저 예제 데이터 파일 중 하나에서 QC 및 트리밍 명령을 실행합니다.

#### 1.1.1. 컨테이너 가져오기

`fastqc`와 `trim_galore`가 모두 설치된 컨테이너 이미지를 가져옵니다:

```bash
docker pull community.wave.seqera.io/library/trim-galore:0.6.10--1bf8ca4e1967cd18
```

??? success "명령 출력"

    ```console
    0.6.10--1bf8ca4e1967cd18: Pulling from library/trim-galore
    dafa2b0c44d2: Pull complete
    dec6b097362e: Pull complete
    f88da01cff0b: Pull complete
    4f4fb700ef54: Pull complete
    92dc97a3ef36: Pull complete
    403f74b0f85e: Pull complete
    10b8c00c10a5: Pull complete
    17dc7ea432cc: Pull complete
    bb36d6c3110d: Pull complete
    0ea1a16bbe82: Pull complete
    030a47592a0a: Pull complete
    32ec762be2d0: Pull complete
    d2cb90387285: Pull complete
    Digest: sha256:4f00e7b2a09f3c8d8a9ce955120e177152fb1e56f63a2a6e186088b1250d9907
    Status: Downloaded newer image for community.wave.seqera.io/library/trim-galore:0.6.10--1bf8ca4e1967cd18
    community.wave.seqera.io/library/trim-galore:0.6.10--1bf8ca4e1967cd18
    ```

이전에 이 이미지를 다운로드한 적이 없다면 완료하는 데 1분 정도 걸릴 수 있습니다.
완료되면 컨테이너 이미지의 로컬 복사본을 갖게 됩니다.

#### 1.1.2. 대화형으로 컨테이너 실행하기

컨테이너를 대화형으로 실행하려면 `-it` 플래그와 함께 `docker run`을 사용하십시오.
`-v ./data:/data` 옵션은 로컬 `data/` 디렉토리를 마운트하여 컨테이너 내부에서 입력 파일에 액세스할 수 있게 합니다.

```bash
docker run -it -v ./data:/data community.wave.seqera.io/library/trim-galore:0.6.10--1bf8ca4e1967cd18
```

??? success "명령 출력"

    ```console
    (base) root@b645838b3314:/tmp#
    ```

프롬프트가 `(base) root@b645838b3314:/tmp#`와 같이 변경되며, 이는 이제 컨테이너 내부에 있음을 나타냅니다.

`/data/reads` 아래에서 시퀀스 데이터 파일을 볼 수 있는지 확인하십시오:

```bash
ls /data/reads
```

??? abstract "디렉토리 내용"

    ```console
    ENCSR000COQ1_1.fastq.gz  ENCSR000COQ2_2.fastq.gz  ENCSR000COR2_1.fastq.gz  ENCSR000CPO1_2.fastq.gz
    ENCSR000COQ1_2.fastq.gz  ENCSR000COR1_1.fastq.gz  ENCSR000COR2_2.fastq.gz  ENCSR000CPO2_1.fastq.gz
    ENCSR000COQ2_1.fastq.gz  ENCSR000COR1_2.fastq.gz  ENCSR000CPO1_1.fastq.gz  ENCSR000CPO2_2.fastq.gz
    ```

이제 첫 번째 명령을 시도할 준비가 되었습니다.

#### 1.1.3. FastQC 명령 실행하기

위에서 참조한 방법은 단일 파일에 대해 QC를 실행하는 명령줄을 제공합니다.
입력 파일만 제공하면 됩니다. 도구는 원본 데이터와 동일한 디렉토리에 출력 파일을 자동으로 생성합니다.

하나의 데이터 파일에 대해 `fastqc` 명령을 실행하십시오:

```bash
fastqc /data/reads/ENCSR000COQ1_1.fastq.gz
```

??? success "명령 출력"

    ```console
    application/gzip
    Started analysis of ENCSR000COQ1_1.fastq.gz
    Approx 5% complete for ENCSR000COQ1_1.fastq.gz
    Approx 10% complete for ENCSR000COQ1_1.fastq.gz
    Approx 15% complete for ENCSR000COQ1_1.fastq.gz
    Approx 20% complete for ENCSR000COQ1_1.fastq.gz
    Approx 25% complete for ENCSR000COQ1_1.fastq.gz
    Approx 30% complete for ENCSR000COQ1_1.fastq.gz
    Approx 35% complete for ENCSR000COQ1_1.fastq.gz
    Approx 40% complete for ENCSR000COQ1_1.fastq.gz
    Approx 45% complete for ENCSR000COQ1_1.fastq.gz
    Approx 50% complete for ENCSR000COQ1_1.fastq.gz
    Approx 55% complete for ENCSR000COQ1_1.fastq.gz
    Approx 60% complete for ENCSR000COQ1_1.fastq.gz
    Approx 65% complete for ENCSR000COQ1_1.fastq.gz
    Approx 70% complete for ENCSR000COQ1_1.fastq.gz
    Approx 75% complete for ENCSR000COQ1_1.fastq.gz
    Approx 80% complete for ENCSR000COQ1_1.fastq.gz
    Approx 85% complete for ENCSR000COQ1_1.fastq.gz
    Approx 90% complete for ENCSR000COQ1_1.fastq.gz
    Approx 95% complete for ENCSR000COQ1_1.fastq.gz
    Analysis complete for ENCSR000COQ1_1.fastq.gz
    ```

이 명령은 매우 빠르게 실행됩니다.
원본 데이터와 동일한 디렉토리에서 출력 파일을 찾을 수 있습니다:

```bash
ls /data/reads/ENCSR000COQ1_1_fastqc*
```

??? abstract "디렉토리 내용"

    ```console
    /data/reads/ENCSR000COQ1_1_fastqc.html  /data/reads/ENCSR000COQ1_1_fastqc.zip
    ```

HTML 보고서와 QC 메트릭이 포함된 ZIP 아카이브가 표시됩니다.
이것으로 첫 번째 단계의 테스트가 완료되었습니다.

#### 1.1.4. Trim Galore로 어댑터 서열 트리밍하기

이제 Cutadapt 및 FastQC를 번들로 제공하는 `trim_galore`를 실행하여 어댑터 서열을 트리밍하고 트리밍 후 QC 메트릭을 수집합니다.
위에서 언급했듯이 소프트웨어는 동일한 컨테이너에 포함되어 있으므로 변경할 필요가 없습니다.

명령은 간단합니다. 트리밍이 완료된 후 QC 수집 단계를 자동으로 실행하도록 `--fastqc` 플래그를 추가하기만 하면 됩니다.

```bash
trim_galore --fastqc /data/reads/ENCSR000COQ1_1.fastq.gz
```

??? success "명령 출력"

    ```console hl_lines="54 55 56 58 59 60"
    Multicore support not enabled. Proceeding with single-core trimming.
    Path to Cutadapt set as: 'cutadapt' (default)
    Cutadapt seems to be working fine (tested command 'cutadapt --version')
    Cutadapt version: 4.9
    single-core operation.
    igzip command line interface 2.31.0
    igzip detected. Using igzip for decompressing

    No quality encoding type selected. Assuming that the data provided uses Sanger encoded Phred scores (default)



    AUTO-DETECTING ADAPTER TYPE
    ===========================
    Attempting to auto-detect adapter type from the first 1 million sequences of the first file (>> /data/reads/ENCSR000COQ1_1.fastq.gz <<)

    Found perfect matches for the following adapter sequences:
    Adapter type	Count	Sequence	Sequences analysed	Percentage
    Illumina	9	AGATCGGAAGAGC	27816	0.03
    smallRNA	0	TGGAATTCTCGG	27816	0.00
    Nextera	0	CTGTCTCTTATA	27816	0.00
    Using Illumina adapter for trimming (count: 9). Second best hit was smallRNA (count: 0)

    Writing report to 'ENCSR000COQ1_1.fastq.gz_trimming_report.txt'

    SUMMARISING RUN PARAMETERS
    ==========================
    Input filename: /data/reads/ENCSR000COQ1_1.fastq.gz
    Trimming mode: single-end
    Trim Galore version: 0.6.10
    Cutadapt version: 4.9
    Number of cores used for trimming: 1
    Quality Phred score cutoff: 20
    Quality encoding type selected: ASCII+33
    Adapter sequence: 'AGATCGGAAGAGC' (Illumina TruSeq, Sanger iPCR; auto-detected)
    Maximum trimming error rate: 0.1 (default)
    Minimum required adapter overlap (stringency): 1 bp
    Minimum required sequence length before a sequence gets removed: 20 bp
    Running FastQC on the data once trimming has completed
    Output file(s) will be GZIP compressed

    Cutadapt seems to be fairly up-to-date (version 4.9). Setting -j 1
    Writing final adapter and quality trimmed output to ENCSR000COQ1_1_trimmed.fq.gz


      >>> Now performing quality (cutoff '-q 20') and adapter trimming in a single pass for the adapter sequence: 'AGATCGGAAGAGC' from file /data/reads/ENCSR000COQ1_1.fastq.gz <<<
    This is cutadapt 4.9 with Python 3.12.7
    Command line parameters: -j 1 -e 0.1 -q 20 -O 1 -a AGATCGGAAGAGC /data/reads/ENCSR000COQ1_1.fastq.gz
    Processing single-end reads on 1 core ...
    Finished in 0.373 s (13.399 µs/read; 4.48 M reads/minute).

    === Summary ===

    Total reads processed:                  27,816
    Reads with adapters:                     9,173 (33.0%)
    Reads written (passing filters):        27,816 (100.0%)

    Total basepairs processed:     2,114,016 bp
    Quality-trimmed:                       0 bp (0.0%)
    Total written (filtered):      2,100,697 bp (99.4%)

    === Adapter 1 ===

    Sequence: AGATCGGAAGAGC; Type: regular 3'; Length: 13; Trimmed: 9173 times

    Minimum overlap: 1
    No. of allowed errors:
    1-9 bp: 0; 10-13 bp: 1

    Bases preceding removed adapters:
      A: 27.4%
      C: 37.4%
      G: 20.9%
      T: 14.3%
      none/other: 0.0%

    Overview of removed sequences
    length	count	expect	max.err	error counts
    1	6229	6954.0	0	6229
    2	2221	1738.5	0	2221
    3	581	434.6	0	581
    4	88	108.7	0	88
    5	33	27.2	0	33
    6	2	6.8	0	2
    7	1	1.7	0	1
    9	1	0.1	0	1
    10	2	0.0	1	2
    12	1	0.0	1	0 1
    14	4	0.0	1	3 1
    16	1	0.0	1	1
    19	1	0.0	1	1
    22	1	0.0	1	1
    29	4	0.0	1	0 4
    33	3	0.0	1	3

    RUN STATISTICS FOR INPUT FILE: /data/reads/ENCSR000COQ1_1.fastq.gz
    =============================================
    27816 sequences processed in total
    Sequences removed because they became shorter than the length cutoff of 20 bp:	0 (0.0%)


      >>> Now running FastQC on the data <<<

    application/gzip
    Started analysis of ENCSR000COQ1_1_trimmed.fq.gz
    Approx 5% complete for ENCSR000COQ1_1_trimmed.fq.gz
    Approx 10% complete for ENCSR000COQ1_1_trimmed.fq.gz
    Approx 15% complete for ENCSR000COQ1_1_trimmed.fq.gz
    Approx 20% complete for ENCSR000COQ1_1_trimmed.fq.gz
    Approx 25% complete for ENCSR000COQ1_1_trimmed.fq.gz
    Approx 30% complete for ENCSR000COQ1_1_trimmed.fq.gz
    Approx 35% complete for ENCSR000COQ1_1_trimmed.fq.gz
    Approx 40% complete for ENCSR000COQ1_1_trimmed.fq.gz
    Approx 45% complete for ENCSR000COQ1_1_trimmed.fq.gz
    Approx 50% complete for ENCSR000COQ1_1_trimmed.fq.gz
    Approx 55% complete for ENCSR000COQ1_1_trimmed.fq.gz
    Approx 60% complete for ENCSR000COQ1_1_trimmed.fq.gz
    Approx 65% complete for ENCSR000COQ1_1_trimmed.fq.gz
    Approx 70% complete for ENCSR000COQ1_1_trimmed.fq.gz
    Approx 75% complete for ENCSR000COQ1_1_trimmed.fq.gz
    Approx 80% complete for ENCSR000COQ1_1_trimmed.fq.gz
    Approx 85% complete for ENCSR000COQ1_1_trimmed.fq.gz
    Approx 90% complete for ENCSR000COQ1_1_trimmed.fq.gz
    Approx 95% complete for ENCSR000COQ1_1_trimmed.fq.gz
    Analysis complete for ENCSR000COQ1_1_trimmed.fq.gz
    ```

출력이 매우 상세하므로 위의 예제에서 가장 관련성 높은 줄을 강조 표시했습니다.
작업 디렉토리에서 출력 파일을 찾을 수 있습니다:

```bash
ls ENCSR000COQ1_1*
```

??? abstract "디렉토리 내용"

    ```console
    ENCSR000COQ1_1.fastq.gz_trimming_report.txt  ENCSR000COQ1_1_trimmed_fastqc.html
    ENCSR000COQ1_1_trimmed.fq.gz                 ENCSR000COQ1_1_trimmed_fastqc.zip
    ```

여기에는 트리밍된 리드, 트리밍 보고서 및 트리밍 후 QC 파일이 포함됩니다.

#### 1.1.5. 출력 파일 이동하기

컨테이너 내부에 남아 있는 모든 것은 향후 작업에 액세스할 수 없으므로 이러한 파일을 마운트된 파일시스템의 디렉토리로 이동해야 합니다.

```bash
mkdir /data/trimmed
mv ENCSR000COQ1_1* /data/trimmed
```

??? abstract "디렉토리 내용"

    ```console
    /data/trimmed
    ├── ENCSR000COQ1_1.fastq.gz_trimming_report.txt
    ├── ENCSR000COQ1_1_trimmed.fq.gz
    ├── ENCSR000COQ1_1_trimmed_fastqc.html
    └── ENCSR000COQ1_1_trimmed_fastqc.zip
    ```

이제 파일이 일반 파일시스템에서 액세스할 수 있습니다.

#### 1.1.6. 컨테이너 종료하기

컨테이너를 종료하려면 `exit`를 입력하십시오.

```bash
exit
```

프롬프트가 정상으로 돌아와야 합니다. 이것으로 처음 두 단계의 테스트가 완료되었습니다.

### 1.2. 참조 게놈에 리드 정렬하기

다음으로 트리밍된 RNAseq 리드를 참조 게놈에 정렬하는 정렬 명령을 실행합니다.

#### 1.2.1. 컨테이너 가져오기

`hisat2`와 `samtools`가 설치된 컨테이너 이미지를 가져옵니다:

```bash
docker pull community.wave.seqera.io/library/hisat2_samtools:5e49f68a37dc010e
```

??? success "명령 출력"

    ```console
    Unable to find image 'community.wave.seqera.io/library/hisat2_samtools:5e49f68a37dc010e' locally
    5e49f68a37dc010e: Pulling from library/hisat2_samtools
    dafa2b0c44d2: Already exists
    dec6b097362e: Already exists
    f88da01cff0b: Already exists
    4f4fb700ef54: Already exists
    92dc97a3ef36: Already exists
    403f74b0f85e: Already exists
    10b8c00c10a5: Already exists
    17dc7ea432cc: Already exists
    bb36d6c3110d: Already exists
    0ea1a16bbe82: Already exists
    030a47592a0a: Already exists
    e74ed5dd390b: Pull complete
    abfcf0185e51: Pull complete
    Digest: sha256:29d8e1a3172a2bdde7be813f7ebec22d331388194a7c0de872b4ccca4bed8f45
    Status: Downloaded newer image for community.wave.seqera.io/library/hisat2_samtools:5e49f68a37dc010e
    ```

일부 레이어가 이전에 가져온 Trim Galore 컨테이너 이미지와 공유되기 때문에 `Already exists`로 표시되는 것을 알 수 있습니다.
결과적으로 이 가져오기는 첫 번째 것보다 빠르게 진행됩니다.

#### 1.2.2. 대화형으로 컨테이너 실행하기

관련 컨테이너 URI로 교체하여 이전과 동일한 방법으로 컨테이너를 대화형으로 실행하십시오.

```bash
docker run -it -v ./data:/data community.wave.seqera.io/library/hisat2_samtools:5e49f68a37dc010e
```

프롬프트가 다시 변경되어 컨테이너 내부에 있음을 나타냅니다.

#### 1.2.3. 게놈 인덱스 파일 생성하기

HISAT2는 게놈 참조가 매우 특정한 형식으로 제공되어야 하며 우리가 제공하는 `genome.fa` FASTA 파일만으로는 사용할 수 없으므로 이 기회에 관련 리소스를 생성합니다.

```bash
hisat2-build /data/genome.fa genome_index
```

??? success "명령 출력"

    ```console hl_lines="1 2 218"
    Settings:
      Output files: "genome_index.*.ht2"
      Line rate: 6 (line is 64 bytes)
      Lines per side: 1 (side is 64 bytes)
      Offset rate: 4 (one in 16)
      FTable chars: 10
      Strings: unpacked
      Local offset rate: 3 (one in 8)
      Local fTable chars: 6
      Local sequence length: 57344
      Local sequence overlap between two consecutive indexes: 1024
      Endianness: little
      Actual local endianness: little
      Sanity checking: disabled
      Assertions: disabled
      Random seed: 0
      Sizeofs: void*:8, int:4, long:8, size_t:8
    Input files DNA, FASTA:
      /data/genome.fa
    Reading reference sizes
      Time reading reference sizes: 00:00:00
    Calculating joined length
    Writing header
    Reserving space for joined string
    Joining reference sequences
      Time to join reference sequences: 00:00:00
      Time to read SNPs and splice sites: 00:00:00
    Using parameters --bmax 6542727 --dcv 1024
      Doing ahead-of-time memory usage test
      Passed!  Constructing with these parameters: --bmax 6542727 --dcv 1024
    Constructing suffix-array element generator
    Building DifferenceCoverSample
      Building sPrime
      Building sPrimeOrder
      V-Sorting samples
      V-Sorting samples time: 00:00:01
      Allocating rank array
      Ranking v-sort output
      Ranking v-sort output time: 00:00:00
      Invoking Larsson-Sadakane on ranks
      Invoking Larsson-Sadakane on ranks time: 00:00:00
      Sanity-checking and returning
    Building samples
    Reserving space for 12 sample suffixes
    Generating random suffixes
    QSorting 12 sample offsets, eliminating duplicates
    QSorting sample offsets, eliminating duplicates time: 00:00:00
    Multikey QSorting 12 samples
      (Using difference cover)
      Multikey QSorting samples time: 00:00:00
    Calculating bucket sizes
    Splitting and merging
      Splitting and merging time: 00:00:00
    Split 1, merged 7; iterating...
    Splitting and merging
      Splitting and merging time: 00:00:00
    Avg bucket size: 4.98493e+06 (target: 6542726)
    Converting suffix-array elements to index image
    Allocating ftab, absorbFtab
    Entering GFM loop
    Getting block 1 of 7
      Reserving size (6542727) for bucket 1
      Calculating Z arrays for bucket 1
      Entering block accumulator loop for bucket 1:
      bucket 1: 10%
      bucket 1: 20%
      bucket 1: 30%
      bucket 1: 40%
      bucket 1: 50%
      bucket 1: 60%
      bucket 1: 70%
      bucket 1: 80%
      bucket 1: 90%
      bucket 1: 100%
      Sorting block of length 3540952 for bucket 1
      (Using difference cover)
      Sorting block time: 00:00:01
    Returning block of 3540953 for bucket 1
    Getting block 2 of 7
      Reserving size (6542727) for bucket 2
      Calculating Z arrays for bucket 2
      Entering block accumulator loop for bucket 2:
      bucket 2: 10%
      bucket 2: 20%
      bucket 2: 30%
      bucket 2: 40%
      bucket 2: 50%
      bucket 2: 60%
      bucket 2: 70%
      bucket 2: 80%
      bucket 2: 90%
      bucket 2: 100%
      Sorting block of length 6195795 for bucket 2
      (Using difference cover)
      Sorting block time: 00:00:01
    Returning block of 6195796 for bucket 2
    Getting block 3 of 7
      Reserving size (6542727) for bucket 3
      Calculating Z arrays for bucket 3
      Entering block accumulator loop for bucket 3:
      bucket 3: 10%
      bucket 3: 20%
      bucket 3: 30%
      bucket 3: 40%
      bucket 3: 50%
      bucket 3: 60%
      bucket 3: 70%
      bucket 3: 80%
      bucket 3: 90%
      bucket 3: 100%
      Sorting block of length 6199288 for bucket 3
      (Using difference cover)
      Sorting block time: 00:00:01
    Returning block of 6199289 for bucket 3
    Getting block 4 of 7
      Reserving size (6542727) for bucket 4
      Calculating Z arrays for bucket 4
      Entering block accumulator loop for bucket 4:
      bucket 4: 10%
      bucket 4: 20%
      bucket 4: 30%
      bucket 4: 40%
      bucket 4: 50%
      bucket 4: 60%
      bucket 4: 70%
      bucket 4: 80%
      bucket 4: 90%
      bucket 4: 100%
      Sorting block of length 6454986 for bucket 4
      (Using difference cover)
      Sorting block time: 00:00:00
    Returning block of 6454987 for bucket 4
    Getting block 5 of 7
      Reserving size (6542727) for bucket 5
      Calculating Z arrays for bucket 5
      Entering block accumulator loop for bucket 5:
      bucket 5: 10%
      bucket 5: 20%
      bucket 5: 30%
      bucket 5: 40%
      bucket 5: 50%
      bucket 5: 60%
      bucket 5: 70%
      bucket 5: 80%
      bucket 5: 90%
      bucket 5: 100%
      Sorting block of length 3493181 for bucket 5
      (Using difference cover)
      Sorting block time: 00:00:00
    Returning block of 3493182 for bucket 5
    Getting block 6 of 7
      Reserving size (6542727) for bucket 6
      Calculating Z arrays for bucket 6
      Entering block accumulator loop for bucket 6:
      bucket 6: 10%
      bucket 6: 20%
      bucket 6: 30%
      bucket 6: 40%
      bucket 6: 50%
      bucket 6: 60%
      bucket 6: 70%
      bucket 6: 80%
      bucket 6: 90%
      bucket 6: 100%
      Sorting block of length 5875908 for bucket 6
      (Using difference cover)
      Sorting block time: 00:00:00
    Returning block of 5875909 for bucket 6
    Getting block 7 of 7
      Reserving size (6542727) for bucket 7
      Calculating Z arrays for bucket 7
      Entering block accumulator loop for bucket 7:
      bucket 7: 10%
      bucket 7: 20%
      bucket 7: 30%
      bucket 7: 40%
      bucket 7: 50%
      bucket 7: 60%
      bucket 7: 70%
      bucket 7: 80%
      bucket 7: 90%
      bucket 7: 100%
      Sorting block of length 3134429 for bucket 7
      (Using difference cover)
      Sorting block time: 00:00:00
    Returning block of 3134430 for bucket 7
    Exited GFM loop
    fchr[A]: 0
    fchr[C]: 9094775
    fchr[G]: 17470759
    fchr[T]: 25839994
    fchr[$]: 34894545
    Exiting GFM::buildToDisk()
    Returning from initFromVector
    Wrote 15826295 bytes to primary GFM file: genome_index.1.ht2
    Wrote 8723644 bytes to secondary GFM file: genome_index.2.ht2
    Re-opening _in1 and _in2 as input streams
    Returning from GFM constructor
    Returning from initFromVector
    Wrote 15353415 bytes to primary GFM file: genome_index.5.ht2
    Wrote 8883598 bytes to secondary GFM file: genome_index.6.ht2
    Re-opening _in5 and _in5 as input streams
    Returning from HGFM constructor
    Headers:
        len: 34894545
        gbwtLen: 34894546
        nodes: 34894546
        sz: 8723637
        gbwtSz: 8723637
        lineRate: 6
        offRate: 4
        offMask: 0xfffffff0
        ftabChars: 10
        eftabLen: 0
        eftabSz: 0
        ftabLen: 1048577
        ftabSz: 4194308
        offsLen: 2180910
        offsSz: 8723640
        lineSz: 64
        sideSz: 64
        sideGbwtSz: 48
        sideGbwtLen: 192
        numSides: 181743
        numLines: 181743
        gbwtTotLen: 11631552
        gbwtTotSz: 11631552
        reverse: 0
        linearFM: Yes
    Total time for call to driver() for forward index: 00:00:12
    ```

출력이 매우 상세하므로 위의 예제에서 일부 관련 줄을 강조 표시했습니다.

이렇게 하면 작업 디렉토리에서 찾을 수 있는 여러 게놈 인덱스 파일이 생성됩니다.

```bash
ls genome_index.*
```

??? abstract "디렉토리 내용"

    ```console
    genome_index.1.ht2  genome_index.3.ht2  genome_index.5.ht2  genome_index.7.ht2
    genome_index.2.ht2  genome_index.4.ht2  genome_index.6.ht2  genome_index.8.ht2
    ```

나중에 이러한 파일이 필요하며, 이를 생성하는 것은 일반적으로 워크플로우의 일부로 수행하고 싶지 않은 작업이므로 필요에 따라 쉽게 전달할 수 있는 게놈 인덱스 파일이 포함된 gzip 압축된 tarball을 생성합니다.

```bash
tar -czvf /data/genome_index.tar.gz genome_index.*
```

??? success "명령 출력"

    ```console
    genome_index.1.ht2
    genome_index.2.ht2
    genome_index.3.ht2
    genome_index.4.ht2
    genome_index.5.ht2
    genome_index.6.ht2
    genome_index.7.ht2
    genome_index.8.ht2
    ```

몇 분 후에 게놈 인덱스 파일이 포함된 결과 `genome_index.tar.gz` tarball을 파일시스템의 `data/` 디렉토리로 이동할 것입니다.
이것은 이 과정의 파트 2에서 유용하게 사용됩니다.

#### 1.2.4. 정렬 명령 실행하기

이제 `hisat2`로 정렬 단계를 수행한 다음 출력을 `samtools`로 파이프하여 출력을 BAM 파일로 작성하는 정렬 명령을 실행할 수 있습니다.

리드 데이터 입력은 이전 단계에서 `trim_galore`로 생성한 `/data/trimmed/ENCSR000COQ1_1_trimmed.fq.gz` 파일입니다.

```bash
hisat2 -x genome_index -U /data/trimmed/ENCSR000COQ1_1_trimmed.fq.gz \
    --new-summary --summary-file ENCSR000COQ1_1_trimmed.hisat2.log | \
    samtools view -bS -o ENCSR000COQ1_1_trimmed.bam
```

??? success "명령 출력"

    ```console hl_lines="6"
    HISAT2 summary stats:
    	Total reads: 27816
    		Aligned 0 time: 1550 (5.57%)
    		Aligned 1 time: 25410 (91.35%)
    		Aligned >1 times: 856 (3.08%)
    	Overall alignment rate: 94.43%
    ```

매우 작은 테스트 파일이므로 거의 즉시 실행됩니다.
실제 규모에서는 훨씬 더 오래 걸릴 수 있습니다.

다시 한 번 작업 디렉토리에서 출력 파일을 찾을 수 있습니다:

```bash
ls ENCSR000COQ1_1*
```

??? abstract "디렉토리 내용"

    ```console title="Output"
    ENCSR000COQ1_1_trimmed.bam  ENCSR000COQ1_1_trimmed.hisat2.log
    ```

정렬은 BAM 파일과 정렬 통계가 포함된 로그 파일을 생성했습니다.

#### 1.2.5. 출력 파일 이동하기

이전과 마찬가지로 출력 파일을 마운트된 파일시스템의 디렉토리로 이동하여 컨테이너를 종료한 후에도 액세스할 수 있도록 합니다.

```bash
mkdir /data/aligned
mv ENCSR000COQ1_1* /data/aligned
```

이것으로 필요한 모든 것을 갖추었습니다.

#### 1.2.6. 컨테이너 종료하기

컨테이너를 종료하려면 `exit`를 입력하십시오.

```bash
exit
```

프롬프트가 정상으로 돌아와야 합니다.
이것으로 단일 샘플 처리 테스트 실행이 완료되었습니다.

!!! example "워크플로우로 작성하기!"

    바로 [파트 2](./02_single-sample.md)로 이동하여 이 분석을 Nextflow 워크플로우로 구현하기 시작하고 싶다면 자유롭게 진행하십시오.
    파트 3으로 넘어가기 전에 두 번째 테스트 라운드를 완료하기 위해 돌아오기만 하면 됩니다.

---

## 2. 다중 샘플 QC 집계

방금 테스트한 명령은 한 번에 하나의 샘플을 처리합니다.
실제로는 일반적으로 많은 샘플을 처리한 다음 모든 샘플에 걸쳐 QC 결과를 집계하여 전체 데이터셋의 품질을 평가해야 합니다.

[MultiQC](https://multiqc.info/)는 많은 일반적인 생물정보학 도구의 QC 보고서를 디렉토리에서 검색하고 이를 단일 포괄적인 HTML 보고서로 집계하는 도구입니다.
FastQC, Cutadapt (Trim Galore를 통해) 및 HISAT2 등의 출력을 인식할 수 있습니다.

여기서는 동일한 샘플별 도구를 통해 두 개의 추가 샘플을 처리한 다음 MultiQC를 사용하여 세 샘플 모두에 걸쳐 QC 보고서를 집계합니다.
이것들이 이 과정의 파트 3에서 Nextflow 워크플로우로 적용할 명령입니다.

1. Trim Galore를 사용하여 추가 샘플에 대한 QC 및 트리밍 실행
2. HISAT2를 사용하여 추가 샘플에 대한 정렬 실행
3. MultiQC를 사용하여 모든 QC 보고서를 포괄적인 보고서로 집계

<figure class="excalidraw">
--8<-- "docs/en/docs/nf4_science/rnaseq/img/rnaseq-method-03.svg"
</figure>

### 2.1. 추가 샘플 QC 및 트리밍

샘플별 QC 및 트리밍 명령은 섹션 1.1에서 실행한 것과 동일합니다.
이미 컨테이너 이미지를 가져왔으므로 직접 실행할 수 있습니다.

#### 2.1.1. 컨테이너 실행하기

섹션 1.1에서 이미 이 컨테이너 이미지를 가져왔으므로 직접 실행할 수 있습니다:

```bash
docker run -it -v ./data:/data community.wave.seqera.io/library/trim-galore:0.6.10--1bf8ca4e1967cd18
```

프롬프트가 변경되어 컨테이너 내부에 있음을 나타냅니다.

#### 2.1.2. 추가 샘플에 대한 QC 및 트리밍 실행

두 개의 추가 샘플에 대해 FastQC 및 Trim Galore를 차례로 실행하십시오.

```bash
fastqc /data/reads/ENCSR000COQ2_1.fastq.gz
fastqc /data/reads/ENCSR000COR1_1.fastq.gz
```

```bash
trim_galore --fastqc /data/reads/ENCSR000COQ2_1.fastq.gz
trim_galore --fastqc /data/reads/ENCSR000COR1_1.fastq.gz
```

완료되면 작업 디렉토리에 두 샘플 모두에 대한 Trim Galore 출력 파일이 있어야 합니다.

#### 2.1.3. 출력 파일 이동하기

Trim Galore 출력 파일을 섹션 1에서 사용한 동일한 디렉토리로 이동하십시오.

```bash
mv ENCSR000COQ2_1* ENCSR000COR1_1* /data/trimmed
```

??? abstract "디렉토리 내용"

    ```console
    /data/trimmed
    ├── ENCSR000COQ1_1.fastq.gz_trimming_report.txt
    ├── ENCSR000COQ1_1_trimmed.fq.gz
    ├── ENCSR000COQ1_1_trimmed_fastqc.html
    ├── ENCSR000COQ1_1_trimmed_fastqc.zip
    ├── ENCSR000COQ2_1.fastq.gz_trimming_report.txt
    ├── ENCSR000COQ2_1_trimmed.fq.gz
    ├── ENCSR000COQ2_1_trimmed_fastqc.html
    ├── ENCSR000COQ2_1_trimmed_fastqc.zip
    ├── ENCSR000COR1_1.fastq.gz_trimming_report.txt
    ├── ENCSR000COR1_1_trimmed.fq.gz
    ├── ENCSR000COR1_1_trimmed_fastqc.html
    └── ENCSR000COR1_1_trimmed_fastqc.zip
    ```

이제 파일이 일반 파일시스템에서 액세스할 수 있습니다.

#### 2.1.4. 컨테이너 종료하기

컨테이너를 종료하려면 `exit`를 입력하십시오.

```bash
exit
```

프롬프트가 정상으로 돌아와야 합니다.

### 2.2. 추가 샘플 정렬

정렬 명령은 섹션 1.2에서 실행한 것과 동일합니다.
원본 인덱스 파일이 더 이상 존재하지 않는 컨테이너 내부에서 생성되었으므로 이전에 저장한 tarball에서 게놈 인덱스를 추출해야 합니다.

#### 2.2.1. 컨테이너 실행하기

섹션 1.2에서 이미 이 컨테이너 이미지를 가져왔으므로 직접 실행할 수 있습니다:

```bash
docker run -it -v ./data:/data community.wave.seqera.io/library/hisat2_samtools:5e49f68a37dc010e
```

프롬프트가 변경되어 컨테이너 내부에 있음을 나타냅니다.

#### 2.2.2. 게놈 인덱스 추출하기

마운트된 파일시스템에 저장한 tarball에서 게놈 인덱스 파일을 추출하십시오:

```bash
tar -xzf /data/genome_index.tar.gz
```

이렇게 하면 작업 디렉토리에 `genome_index.*` 파일이 복원됩니다.

#### 2.2.3. 추가 샘플에 대한 정렬 실행

새로 트리밍된 두 샘플에 대해 HISAT2 정렬을 차례로 실행하십시오.

```bash
hisat2 -x genome_index -U /data/trimmed/ENCSR000COQ2_1_trimmed.fq.gz \
    --new-summary --summary-file ENCSR000COQ2_1_trimmed.hisat2.log | \
    samtools view -bS -o ENCSR000COQ2_1_trimmed.bam
```

??? success "명령 출력"

    ```console hl_lines="6"
    HISAT2 summary stats:
    	Total reads: 18736
    		Aligned 0 time: 1531 (8.17%)
    		Aligned 1 time: 16726 (89.27%)
    		Aligned >1 times: 479 (2.56%)
    	Overall alignment rate: 91.83%
    ```

```bash
hisat2 -x genome_index -U /data/trimmed/ENCSR000COR1_1_trimmed.fq.gz \
    --new-summary --summary-file ENCSR000COR1_1_trimmed.hisat2.log | \
    samtools view -bS -o ENCSR000COR1_1_trimmed.bam
```

??? success "명령 출력"

    ```console hl_lines="6"
    HISAT2 summary stats:
    	Total reads: 38056
    		Aligned 0 time: 2311 (6.07%)
    		Aligned 1 time: 33289 (87.47%)
    		Aligned >1 times: 2456 (6.45%)
    	Overall alignment rate: 93.93%
    ```

완료되면 작업 디렉토리에 두 샘플 모두에 대한 BAM 및 로그 파일이 있어야 합니다.

#### 2.2.4. 출력 파일 이동하기

정렬 출력 파일을 섹션 1에서 사용한 동일한 디렉토리로 이동하십시오.

```bash
mv ENCSR000COQ2_1* ENCSR000COR1_1* /data/aligned
```

??? abstract "디렉토리 내용"

    ```console
    /data/aligned
    ├── ENCSR000COQ1_1_trimmed.bam
    ├── ENCSR000COQ1_1_trimmed.hisat2.log
    ├── ENCSR000COQ2_1_trimmed.bam
    ├── ENCSR000COQ2_1_trimmed.hisat2.log
    ├── ENCSR000COR1_1_trimmed.bam
    └── ENCSR000COR1_1_trimmed.hisat2.log
    ```

이제 파일이 일반 파일시스템에서 액세스할 수 있습니다.

#### 2.2.5. 컨테이너 종료하기

컨테이너를 종료하려면 `exit`를 입력하십시오.

```bash
exit
```

프롬프트가 정상으로 돌아와야 합니다.

### 2.3. 포괄적인 QC 보고서 생성하기

이제 세 샘플에 대한 QC, 트리밍 및 정렬 출력이 있으므로 MultiQC를 사용하여 이를 단일 보고서로 집계할 수 있습니다.
MultiQC는 호환 가능한 QC 보고서를 찾기 위해 디렉토리를 검색하고 찾은 모든 것을 집계합니다.

#### 2.3.1. 컨테이너 가져오기

`multiqc`가 설치된 컨테이너 이미지를 가져옵니다:

```bash
docker pull community.wave.seqera.io/library/pip_multiqc:a3c26f6199d64b7c
```

??? success "명령 출력"

    ```console
    a3c26f6199d64b7c: Pulling from library/pip_multiqc
    dafa2b0c44d2: Already exists
    dec6b097362e: Already exists
    f88da01cff0b: Already exists
    4f4fb700ef54: Already exists
    92dc97a3ef36: Already exists
    403f74b0f85e: Already exists
    10b8c00c10a5: Already exists
    17dc7ea432cc: Already exists
    bb36d6c3110d: Already exists
    0ea1a16bbe82: Already exists
    030a47592a0a: Already exists
    2ed162b168e8: Pull complete
    ca06fe148f21: Pull complete
    Digest: sha256:af0e9de56896805aa2a065f7650362956f4213d99e95314f6fec472c6a3bf091
    Status: Downloaded newer image for community.wave.seqera.io/library/pip_multiqc:a3c26f6199d64b7c
    community.wave.seqera.io/library/pip_multiqc:a3c26f6199d64b7c
    ```

일부 레이어가 이전에 가져온 컨테이너 이미지와 공유되기 때문에 `Already exists`로 표시되는 것을 알 수 있습니다.
결과적으로 이 가져오기는 이전 것들보다 빠르게 진행됩니다.

#### 2.3.2. 대화형으로 컨테이너 실행하기

이전과 마찬가지로 데이터 디렉토리를 마운트하여 컨테이너를 대화형으로 실행하십시오.

```bash
docker run -it -v ./data:/data community.wave.seqera.io/library/pip_multiqc:a3c26f6199d64b7c
```

프롬프트가 변경되어 컨테이너 내부에 있음을 나타냅니다.

#### 2.3.3. MultiQC 명령 실행하기

세 샘플 모두에 대한 QC 관련 출력 파일을 저장한 디렉토리를 가리키며 `multiqc`를 실행하십시오.
`-n` 플래그는 출력 보고서의 이름을 설정합니다.

```bash
multiqc /data/reads /data/trimmed /data/aligned -n all_samples_QC
```

??? success "명령 출력"

    ```console hl_lines="8 9 10 11 12"

    /// MultiQC 🔍 v1.32

          file_search | Search path: /data/reads
          file_search | Search path: /data/trimmed
          file_search | Search path: /data/aligned
            searching | ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━ 100% 36/36
               hisat2 | Found 3 reports
             cutadapt | Found 3 reports
               fastqc | Found 3 reports
        write_results | Data        : all_samples_QC_data
        write_results | Report      : all_samples_QC.html
              multiqc | MultiQC complete
    ```

여기서 도구가 우리가 생성한 세 가지 QC 보고서를 모두 찾았음을 알 수 있습니다: `fastqc`로 수행한 초기 QC, (`trim_galore`를 통해 생성된) `cutadapt`의 트리밍 후 보고서, 그리고 `hisat2`가 생성한 정렬 후 QC입니다.

출력 파일은 작업 디렉토리에 있습니다:

```bash
ls all_samples_QC*
```

??? abstract "디렉토리 내용"

    ```console
    all_samples_QC.html

    all_samples_QC_data:
    cutadapt_filtered_reads_plot.txt                     multiqc.log
    cutadapt_trimmed_sequences_plot_3_Counts.txt         multiqc.parquet
    cutadapt_trimmed_sequences_plot_3_Obs_Exp.txt        multiqc_citations.txt
    fastqc-status-check-heatmap.txt                      multiqc_cutadapt.txt
    fastqc_adapter_content_plot.txt                      multiqc_data.json
    fastqc_overrepresented_sequences_plot.txt            multiqc_fastqc.txt
    fastqc_per_base_n_content_plot.txt                   multiqc_general_stats.txt
    fastqc_per_base_sequence_quality_plot.txt            multiqc_hisat2.txt
    fastqc_per_sequence_gc_content_plot_Counts.txt       multiqc_software_versions.txt
    fastqc_per_sequence_gc_content_plot_Percentages.txt  multiqc_sources.txt
    fastqc_per_sequence_quality_scores_plot.txt
    fastqc_sequence_counts_plot.txt
    fastqc_sequence_duplication_levels_plot.txt
    fastqc_top_overrepresented_sequences_table.txt
    hisat2_se_plot.txt
    llms-full.txt
    ```

주요 출력은 기본 메트릭이 포함된 데이터 디렉토리와 함께 제공되는 `all_samples_QC.html` 보고서입니다.

#### 2.3.4. 출력 파일 이동하기

보고서와 데이터 디렉토리를 마운트된 파일시스템으로 이동하십시오.

```bash
mkdir /data/multiqc
mv all_samples_QC* /data/multiqc
```

이제 파일이 일반 파일시스템에서 액세스할 수 있습니다.

#### 2.3.5. 컨테이너 종료하기

컨테이너를 종료하려면 `exit`를 입력하십시오.

```bash
exit
```

프롬프트가 정상으로 돌아와야 합니다.
이것으로 모든 RNAseq 처리 명령의 테스트가 완료되었습니다.

---

### 핵심 정리

각각의 컨테이너에서 FastQC, Trim Galore, HISAT2 및 MultiQC 명령을 실행하는 방법을 알게 되었습니다. 여기에는 여러 샘플을 처리하고 QC 보고서를 집계하는 방법이 포함됩니다.

### 다음 단계

잠시 휴식을 취한 다음 [파트 2](./02_single-sample.md)로 이동하여 동일한 명령을 컨테이너를 사용하여 작업을 실행하는 워크플로우로 적용하는 방법을 학습합니다.
