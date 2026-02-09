# 파트 1: 방법론 개요 및 수동 테스트

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice } AI 기반 번역 - [자세히 알아보고 개선 사항 제안하기](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

변이 호출(Variant calling)은 참조 게놈에 대비하여 게놈 서열의 변이를 식별하는 것을 목표로 하는 게놈 분석 방법입니다.
여기서는 전체 게놈 시퀀싱 데이터에서 짧은 배선 변이(short germline variants), 즉 SNP와 indel을 호출하기 위해 설계된 도구와 방법을 사용합니다.

![GATK 파이프라인](img/gatk-pipeline.png)

완전한 변이 호출 파이프라인은 일반적으로 참조 게놈에 대한 매핑(때로는 게놈 정렬이라고도 함)과 변이 필터링 및 우선순위 지정을 포함하여 많은 단계를 포함합니다.
이 과정에서는 간결성을 위해 변이 호출 부분에만 집중하겠습니다.

### 방법론

전체 게놈 시퀀싱 샘플에 변이 호출을 적용하여 배선 SNP와 indel을 식별하는 두 가지 방법을 보여드리겠습니다.
먼저 각 샘플에서 독립적으로 변이를 호출하는 간단한 **샘플별 접근법**부터 시작하겠습니다.
그런 다음 여러 샘플을 함께 분석하여 더 정확하고 유익한 결과를 생성하는 더 정교한 **공동 호출 접근법**을 보여드리겠습니다.

두 접근법 모두에 대한 워크플로우 코드를 작성하기 전에, 먼저 일부 테스트 데이터에서 명령을 수동으로 실행해 보겠습니다.

### 데이터셋

다음 데이터와 관련 리소스를 제공합니다:

- 인간 염색체 20(hg19/b37에서)의 작은 영역과 그 보조 파일(인덱스 및 서열 사전)로 구성된 **참조 게놈**.
- 가족 트리오(어머니, 아버지, 아들)에 해당하는 **세 개의 전체 게놈 시퀀싱 샘플**. 파일 크기를 작게 유지하기 위해 염색체 20의 작은 데이터 조각으로 부분집합화되었습니다.
  이것은 이미 참조 게놈에 매핑된 Illumina 단일 리드 시퀀싱 데이터로, [BAM](https://samtools.github.io/hts-specs/SAMv1.pdf) 형식(SAM(Sequence Alignment Map)의 압축 버전인 Binary Alignment Map)으로 제공됩니다.
- 샘플에 변이 호출에 적합한 데이터가 있는 게놈 상의 좌표인 **게놈 간격 목록**, BED 형식으로 제공됩니다.

### 소프트웨어

관련된 두 가지 주요 도구는 서열 정렬 파일을 조작하는 데 널리 사용되는 툴킷인 [Samtools](https://www.htslib.org/)와, Broad Institute에서 개발한 변이 발견을 위한 도구 세트인 [GATK](https://gatk.broadinstitute.org/)(Genome Analysis Toolkit)입니다.

이러한 도구는 GitHub Codespaces 환경에 설치되어 있지 않으므로 컨테이너를 통해 사용하겠습니다([Hello Containers](../../hello_nextflow/05_hello_containers.md) 참조).

!!! note "참고"

     `nf4-science/genomics` 디렉토리에 있는지 확인하세요. `pwd`를 입력할 때 표시되는 경로의 마지막 부분이 `genomics`이어야 합니다.

---

## 1. 샘플별 변이 호출

샘플별 변이 호출은 각 샘플을 독립적으로 처리합니다. 변이 호출기는 한 번에 한 샘플의 시퀀싱 데이터를 검사하고 샘플이 참조와 다른 위치를 식별합니다.

이 섹션에서는 샘플별 변이 호출 접근법을 구성하는 두 가지 명령을 테스트합니다: Samtools를 사용한 BAM 파일 인덱싱과 GATK HaplotypeCaller를 사용한 변이 호출.
이것들이 이 과정의 파트 2에서 Nextflow 워크플로우로 적용할 명령입니다.

1. [Samtools](https://www.htslib.org/)를 사용하여 BAM 입력 파일에 대한 인덱스 파일 생성
2. 인덱싱된 BAM 파일에서 GATK HaplotypeCaller를 실행하여 VCF(Variant Call Format)로 샘플별 변이 호출 생성

<figure class="excalidraw">
--8<-- "docs/en/docs/nf4_science/genomics/img/hello-gatk-1.svg"
</figure>

먼저 하나의 샘플에서만 두 명령을 테스트하는 것으로 시작합니다.

### 1.1. Samtools로 BAM 입력 파일 인덱싱

인덱스 파일은 생물정보학 파일 형식의 일반적인 기능입니다. 인덱스 파일은 메인 파일의 구조에 대한 정보를 포함하여 GATK와 같은 도구가 전체 파일을 읽지 않고도 데이터의 부분집합에 접근할 수 있게 합니다.
이것은 이러한 파일이 얼마나 커질 수 있는지를 고려할 때 중요합니다.

BAM 파일은 종종 인덱스 없이 제공되므로, 많은 분석 워크플로우의 첫 번째 단계는 `samtools index`를 사용하여 인덱스를 생성하는 것입니다.

Samtools 컨테이너를 가져와서 대화식으로 실행하고 BAM 파일 중 하나에서 `samtools index` 명령을 실행할 것입니다.

#### 1.1.1. Samtools 컨테이너 가져오기

`docker pull` 명령을 실행하여 Samtools 컨테이너 이미지를 다운로드하세요:

```bash
docker pull community.wave.seqera.io/library/samtools:1.20--b5dfbd93de237464
```

??? success "명령 출력"

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

이 이미지를 이전에 다운로드하지 않았다면 완료하는 데 1분 정도 걸릴 수 있습니다.
완료되면 컨테이너 이미지의 로컬 복사본이 생성됩니다.

#### 1.1.2. Samtools 컨테이너를 대화식으로 실행

컨테이너를 대화식으로 실행하려면 `-it` 플래그와 함께 `docker run`을 사용하세요.
`-v ./data:/data` 옵션은 로컬 `data` 디렉토리를 컨테이너에 마운트하여 도구가 입력 파일에 접근할 수 있도록 합니다.

```bash
docker run -it -v ./data:/data community.wave.seqera.io/library/samtools:1.20--b5dfbd93de237464
```

프롬프트가 `(base) root@a1b2c3d4e5f6:/tmp#`와 같이 변경되어 컨테이너 내부에 있음을 나타냅니다.
데이터 파일은 `/data` 아래에서 접근할 수 있습니다.

#### 1.1.3. 인덱싱 명령 실행

[Samtools 문서](https://www.htslib.org/doc/samtools-index.html)는 BAM 파일을 인덱싱하기 위해 실행할 명령줄을 제공합니다.

입력 파일만 제공하면 됩니다. 도구는 입력 파일명에 `.bai`를 추가하여 출력 파일명을 자동으로 생성합니다.

```bash
samtools index /data/bam/reads_mother.bam
```

??? abstract "디렉토리 내용"

    ```console
    data/bam/
    ├── reads_father.bam
    ├── reads_mother.bam
    ├── reads_mother.bam.bai
    └── reads_son.bam
    ```

이제 원본 BAM 입력 파일과 같은 디렉토리에 `reads_mother.bam.bai`라는 파일이 표시될 것입니다.

#### 1.1.4. Samtools 컨테이너 종료

컨테이너를 종료하려면 `exit`를 입력하세요.

```bash
exit
```

프롬프트가 컨테이너를 시작하기 전의 상태로 돌아가야 합니다.

### 1.2. GATK HaplotypeCaller로 변이 호출

GATK 컨테이너를 가져와서 대화식으로 실행하고 방금 인덱싱한 BAM 파일에서 `gatk HaplotypeCaller` 명령을 실행할 것입니다.

#### 1.2.1. GATK 컨테이너 가져오기

`docker pull` 명령을 실행하여 GATK 컨테이너 이미지를 다운로드하세요:

```bash
docker pull community.wave.seqera.io/library/gatk4:4.5.0.0--730ee8817e436867
```

??? success "명령 출력"

    일부 레이어는 앞서 가져온 Samtools 컨테이너 이미지와 공유되므로 `Already exists`로 표시됩니다.

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

두 컨테이너 이미지가 대부분의 레이어를 공유하기 때문에 첫 번째 pull보다 빨라야 합니다.

#### 1.2.2. GATK 컨테이너를 대화식으로 실행

Samtools에서 했던 것과 마찬가지로 data 디렉토리를 마운트하여 GATK 컨테이너를 대화식으로 실행하세요.

```bash
docker run -it -v ./data:/data community.wave.seqera.io/library/gatk4:4.5.0.0--730ee8817e436867
```

프롬프트가 변경되어 GATK 컨테이너 내부에 있음을 나타냅니다.

#### 1.2.3. 변이 호출 명령 실행

[GATK 문서](https://gatk.broadinstitute.org/hc/en-us/articles/21905025322523-HaplotypeCaller)는 BAM 파일에서 변이 호출을 수행하기 위해 실행할 명령줄을 제공합니다.

BAM 입력 파일(`-I`)과 참조 게놈(`-R`), 출력 파일 이름(`-O`), 분석할 게놈 간격 목록(`-L`)을 제공해야 합니다.

그러나 인덱스 파일의 경로를 지정할 필요는 없습니다. 도구는 확립된 명명 규칙과 공동 배치 규칙에 따라 같은 디렉토리에서 자동으로 찾습니다.
참조 게놈의 보조 파일(인덱스 및 서열 사전 파일, `*.fai` 및 `*.dict`)에도 동일하게 적용됩니다.

```bash
gatk HaplotypeCaller \
        -R /data/ref/ref.fasta \
        -I /data/bam/reads_mother.bam \
        -O reads_mother.vcf \
        -L /data/ref/intervals.bed
```

??? success "명령 출력"

    도구는 자세한 로깅 출력을 생성합니다. 강조 표시된 줄은 성공적으로 완료되었음을 확인합니다.

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

출력 파일 `reads_mother.vcf`는 컨테이너의 작업 디렉토리 내부에 생성되므로, 출력 파일 경로를 변경하지 않으면 VS Code 파일 탐색기에서 볼 수 없습니다.
그러나 작은 테스트 파일이므로 `cat` 명령으로 열어서 내용을 볼 수 있습니다.
파일의 시작 부분까지 스크롤하면 여러 줄의 메타데이터로 구성된 헤더와 그 뒤에 한 줄에 하나씩 나열된 변이 호출 목록을 찾을 수 있습니다.

??? abstract "파일 내용"

    ```console title="reads_mother.vcf" linenums="26"
    #CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	reads_mother
    20_10037292_10066351	3480	.	C	CT	503.03	.	AC=2;AF=1.00;AN=2;DP=23;ExcessHet=0.0000;FS=0.000;MLEAC=2;MLEAF=1.00;MQ=60.00;QD=27.95;SOR=1.179	GT:AD:DP:GQ:PL	1/1:0,18:18:54:517,54,0
    20_10037292_10066351	3520	.	AT	A	609.03	.	AC=2;AF=1.00;AN=2;DP=18;ExcessHet=0.0000;FS=0.000;MLEAC=2;MLEAF=1.00;MQ=60.00;QD=33.83;SOR=0.693	GT:AD:DP:GQ:PL	1/1:0,18:18:54:623,54,0
    20_10037292_10066351	3529	.	T	A	155.64	.	AC=1;AF=0.500;AN=2;BaseQRankSum=-0.544;DP=21;ExcessHet=0.0000;FS=1.871;MLEAC=1;MLEAF=0.500;MQ=60.00;MQRankSum=0.000;QD=7.78;ReadPosRankSum=-1.158;SOR=1.034	GT:AD:DP:GQ:PL	0/1:12,8:20:99:163,0,328
    ```

각 줄은 샘플의 시퀀싱 데이터에서 식별된 가능한 변이를 설명합니다. VCF 형식 해석에 대한 지침은 [이 유용한 문서](https://www.ebi.ac.uk/training/online/courses/human-genetic-variation-introduction/variant-identification-and-analysis/understanding-vcf-format/)를 참조하세요.

출력 VCF 파일에는 GATK에서 자동으로 생성한 `reads_mother.vcf.idx`라는 인덱스 파일이 함께 제공됩니다.
이것은 BAM 인덱스 파일과 동일한 기능을 가지고 있어, 도구가 전체 파일을 로드하지 않고도 데이터의 부분집합을 검색하고 가져올 수 있게 합니다.

#### 1.2.4. GATK 컨테이너 종료

컨테이너를 종료하려면 `exit`를 입력하세요.

```bash
exit
```

프롬프트가 정상으로 돌아가야 합니다.
이것으로 샘플별 변이 호출 테스트가 완료되었습니다.

---

## 2. 코호트에 대한 공동 호출

방금 사용한 변이 호출 접근법은 샘플별 변이 호출을 생성합니다.
이것은 각 샘플의 변이를 독립적으로 살펴보는 데는 좋지만, 제한적인 정보를 제공합니다.
여러 샘플에서 변이 호출이 어떻게 다른지 살펴보는 것이 더 흥미로운 경우가 많습니다.
GATK는 이 목적을 위해 공동 변이 호출이라는 대안 방법을 제공합니다.

공동 변이 호출은 각 샘플에 대해 GVCF(Genomic VCF)라는 특수한 종류의 변이 출력을 생성한 다음, 모든 샘플의 GVCF 데이터를 결합하고 '공동 유전형 분석(joint genotyping)' 통계 분석을 실행하는 것을 포함합니다.

![공동 분석](img/joint-calling.png)

샘플의 GVCF의 특별한 점은 프로그램이 변이의 증거를 찾은 위치뿐만 아니라 게놈의 대상 영역에 있는 모든 위치에 대한 서열 데이터 통계를 요약하는 레코드를 포함한다는 것입니다.
이것은 공동 유전형 분석 계산에 중요합니다([추가 정보](https://gatk.broadinstitute.org/hc/en-us/articles/360035890431-The-logic-of-joint-calling-for-germline-short-variants)).

GVCF는 추가 매개변수(`-ERC GVCF`)와 함께 방금 테스트한 것과 동일한 도구인 GATK HaplotypeCaller에 의해 생성됩니다.
GVCF를 결합하는 것은 GATK GenomicsDBImport로 수행되며, 이것은 샘플별 호출을 데이터 저장소(데이터베이스와 유사)로 결합합니다.
실제 '공동 유전형 분석' 분석은 GATK GenotypeGVCFs로 수행됩니다.

여기서는 GVCF를 생성하고 공동 유전형 분석을 실행하는 데 필요한 명령을 테스트합니다.
이것들이 이 과정의 파트 3에서 Nextflow 워크플로우로 적용할 명령입니다.

1. Samtools를 사용하여 각 BAM 입력 파일에 대한 인덱스 파일 생성
2. 각 BAM 입력 파일에서 GATK HaplotypeCaller를 실행하여 샘플별 게놈 변이 호출의 GVCF 생성
3. 모든 GVCF를 수집하여 GenomicsDB 데이터 저장소로 결합
4. 결합된 GVCF 데이터 저장소에서 공동 유전형 분석을 실행하여 코호트 수준 VCF 생성

<figure class="excalidraw">
--8<-- "docs/en/docs/nf4_science/genomics/img/hello-gatk-2.svg"
</figure>

이제 세 개의 BAM 파일을 모두 인덱싱하는 것부터 시작하여 이러한 모든 명령을 테스트해야 합니다.

### 2.1. 세 샘플 모두에 대한 BAM 파일 인덱싱

위의 첫 번째 섹션에서는 하나의 BAM 파일만 인덱싱했습니다.
이제 GATK HaplotypeCaller가 처리할 수 있도록 세 샘플 모두를 인덱싱해야 합니다.

#### 2.1.1. Samtools 컨테이너를 대화식으로 실행

이미 Samtools 컨테이너 이미지를 가져왔으므로 바로 실행할 수 있습니다:

```bash
docker run -it -v ./data:/data community.wave.seqera.io/library/samtools:1.20--b5dfbd93de237464
```

프롬프트가 변경되어 이전과 같이 data 디렉토리가 마운트된 컨테이너 내부에 있음을 나타냅니다.

#### 2.1.2. 세 샘플 모두에서 인덱싱 명령 실행

세 개의 BAM 파일 각각에서 인덱싱 명령을 실행하세요:

```bash
samtools index /data/bam/reads_mother.bam
samtools index /data/bam/reads_father.bam
samtools index /data/bam/reads_son.bam
```

??? abstract "디렉토리 내용"

    ```console
    data/bam/
    ├── reads_father.bam
    ├── reads_father.bam.bai
    ├── reads_mother.bam
    ├── reads_mother.bam.bai
    ├── reads_son.bam
    └── reads_son.bam.bai
    ```

이것은 해당 BAM 파일과 같은 디렉토리에 인덱스 파일을 생성해야 합니다.

#### 2.1.3. Samtools 컨테이너 종료

컨테이너를 종료하려면 `exit`를 입력하세요.

```bash
exit
```

프롬프트가 정상으로 돌아가야 합니다.

### 2.2. 세 샘플 모두에 대한 GVCF 생성

공동 유전형 분석 단계를 실행하려면 세 샘플 모두에 대한 GVCF가 필요합니다.

#### 2.2.1. GATK 컨테이너를 대화식으로 실행

이미 이전에 GATK 컨테이너 이미지를 가져왔으므로 바로 실행할 수 있습니다:

```bash
docker run -it -v ./data:/data community.wave.seqera.io/library/gatk4:4.5.0.0--730ee8817e436867
```

프롬프트가 변경되어 GATK 컨테이너 내부에 있음을 나타냅니다.

#### 2.2.2. GVCF 옵션과 함께 변이 호출 명령 실행

게놈 VCF(GVCF)를 생성하기 위해 기본 명령에 `-ERC GVCF` 옵션을 추가하여 HaplotypeCaller의 GVCF 모드를 켭니다.

또한 출력 파일의 파일 확장자를 `.vcf`에서 `.g.vcf`로 변경합니다.
이것은 기술적으로 요구 사항은 아니지만 강력히 권장되는 규칙입니다.

```bash
gatk HaplotypeCaller \
        -R /data/ref/ref.fasta \
        -I /data/bam/reads_mother.bam \
        -O reads_mother.g.vcf \
        -L /data/ref/intervals.bed \
        -ERC GVCF
```

??? success "명령 출력"

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

이것은 컨테이너의 현재 작업 디렉토리에 GVCF 출력 파일 `reads_mother.g.vcf`를 생성합니다.

내용을 보기 위해 `cat` 명령을 실행하면, 섹션 1에서 생성한 동등한 VCF보다 훨씬 더 긴 것을 볼 수 있습니다. 파일의 시작 부분까지 스크롤할 수도 없고, 대부분의 줄이 VCF에서 본 것과 상당히 다르게 보입니다.

??? abstract "파일 내용"

    ```console title="reads_mother.g.vcf" linenums="1674"
    20_10037292_10066351    14714   .       T       <NON_REF>       .       .       END=14718       GT:DP:GQ:MIN_DP:PL       0/0:37:99:37:0,99,1192
    20_10037292_10066351    14719   .       T       <NON_REF>       .       .       END=14719       GT:DP:GQ:MIN_DP:PL       0/0:36:82:36:0,82,1087
    20_10037292_10066351    14720   .       T       <NON_REF>       .       .       END=14737       GT:DP:GQ:MIN_DP:PL       0/0:42:99:37:0,100,1160
    ```

이것들은 변이 호출기가 변이의 증거를 찾지 못한 비변이 영역을 나타내므로, 변이가 없다는 것에 대한 신뢰 수준을 설명하는 일부 통계를 포착했습니다.
이를 통해 두 가지 매우 다른 경우를 구별할 수 있습니다: (1) 샘플이 동형접합 참조임을 보여주는 양질의 데이터가 있는 경우, (2) 어느 쪽으로든 결정하기에 충분한 양질의 데이터를 사용할 수 없는 경우.

GVCF에는 일반적으로 이러한 비변이 줄이 많이 있으며, 그 사이에 더 적은 수의 변이 레코드가 산재되어 있습니다.
GVCF에서 `head -176`을 실행하여 파일의 처음 176줄만 로드하여 실제 변이 호출을 찾아보세요.

??? abstract "파일 내용"

    ```console title="reads_mother.g.vcf" linenums="174"
    20_10037292_10066351    3479    .       T       <NON_REF>       .       .       END=3479        GT:DP:GQ:MIN_DP:PL       0/0:34:36:34:0,36,906
    20_10037292_10066351    3480    .       C       CT,<NON_REF>    503.03  .       DP=23;ExcessHet=0.0000;MLEAC=2,0;MLEAF=1.00,0.00;RAW_MQandDP=82800,23    GT:AD:DP:GQ:PL:SB       1/1:0,18,0:18:54:517,54,0,517,54,517:0,0,7,11
    20_10037292_10066351    3481    .       T       <NON_REF>       .       .       END=3481        GT:DP:GQ:MIN_DP:PL       0/0:21:51:21:0,51,765
    ```

두 번째 줄은 파일의 첫 번째 변이 레코드를 보여주며, 이것은 앞서 살펴본 VCF 파일의 첫 번째 변이에 해당합니다.

원본 VCF와 마찬가지로 출력 GVCF 파일에도 `reads_mother.g.vcf.idx`라는 인덱스 파일이 함께 제공됩니다.

#### 2.2.3. 나머지 두 샘플에서 프로세스 반복

아래 명령을 하나씩 실행하여 나머지 두 샘플에 대한 GVCF를 생성하세요.

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

완료되면 현재 디렉토리에 `.g.vcf`로 끝나는 세 개의 파일(샘플당 하나)과 `.g.vcf.idx`로 끝나는 각각의 인덱스 파일이 있어야 합니다.

하지만 컨테이너를 종료하지 마세요!
다음 단계에서 동일한 컨테이너를 사용할 것입니다.

### 2.3. 공동 유전형 분석 실행

이제 모든 GVCF가 있으므로 샘플 코호트에 대한 변이 호출을 생성하는 공동 유전형 분석 접근법을 시도해 볼 수 있습니다.
이것은 모든 GVCF의 데이터를 데이터 저장소로 결합한 다음 공동 유전형 분석을 실행하여 최종 공동 호출된 변이의 VCF를 생성하는 2단계 방법입니다.

#### 2.3.1. 모든 샘플별 GVCF 결합

이 첫 번째 단계는 GATK의 또 다른 도구인 GenomicsDBImport를 사용하여 모든 GVCF의 데이터를 GenomicsDB 데이터 저장소로 결합합니다.

```bash
gatk GenomicsDBImport \
    -V reads_mother.g.vcf \
    -V reads_father.g.vcf \
    -V reads_son.g.vcf \
    -L /data/ref/intervals.bed \
    --genomicsdb-workspace-path family_trio_gdb
```

??? success "명령 출력"

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

이 단계의 출력은 사실상 여러 개의 다른 파일 형태로 결합된 변이 데이터를 담고 있는 추가 중첩 디렉토리 세트를 포함하는 디렉토리입니다.
둘러볼 수는 있지만 이 데이터 저장소 형식은 사람이 직접 읽도록 만들어지지 않았다는 것을 빠르게 알 수 있습니다.

!!! note "참고"

    GATK에는 필요에 따라 데이터 저장소에서 변이 호출 데이터를 검사하고 추출할 수 있게 해주는 도구가 포함되어 있습니다.

#### 2.3.2. 실제 공동 유전형 분석 실행

이 두 번째 단계는 GATK의 또 다른 도구인 GenotypeGVCFs를 사용하여 코호트의 모든 샘플에서 사용 가능한 데이터를 고려하여 변이 통계와 개별 유전형을 재계산합니다.

```bash
gatk GenotypeGVCFs \
    -R /data/ref/ref.fasta \
    -V gendb://family_trio_gdb \
    -O family_trio.vcf
```

??? success "명령 출력"

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

이것은 컨테이너의 현재 작업 디렉토리에 VCF 출력 파일 `family_trio.vcf`를 생성합니다.
또 다른 적당히 작은 파일이므로 이 파일을 `cat`으로 내용을 보고 처음 몇 개의 변이 줄을 찾기 위해 위로 스크롤할 수 있습니다.

??? abstract "파일 내용"

    ```console title="family_trio.vcf" linenums="40"
    #CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO    FORMAT  reads_father    reads_mother    reads_son
    20_10037292_10066351    3480    .       C       CT      1625.89 .       AC=5;AF=0.833;AN=6;BaseQRankSum=0.220;DP=85;ExcessHet=0.0000;FS=2.476;MLEAC=5;MLEAF=0.833;MQ=60.00;MQRankSum=0.00;QD=21.68;ReadPosRankSum=-1.147e+00;SOR=0.487    GT:AD:DP:GQ:PL  0/1:15,16:31:99:367,0,375       1/1:0,18:18:54:517,54,0 1/1:0,26:26:78:756,78,0
    20_10037292_10066351    3520    .       AT      A       1678.89 .       AC=5;AF=0.833;AN=6;BaseQRankSum=1.03;DP=80;ExcessHet=0.0000;FS=2.290;MLEAC=5;MLEAF=0.833;MQ=60.00;MQRankSum=0.00;QD=22.39;ReadPosRankSum=0.701;SOR=0.730 GT:AD:DP:GQ:PL   0/1:18,13:31:99:296,0,424       1/1:0,18:18:54:623,54,0 1/1:0,26:26:78:774,78,0
    20_10037292_10066351    3529    .       T       A       154.29  .       AC=1;AF=0.167;AN=6;BaseQRankSum=-5.440e-01;DP=104;ExcessHet=0.0000;FS=1.871;MLEAC=1;MLEAF=0.167;MQ=60.00;MQRankSum=0.00;QD=7.71;ReadPosRankSum=-1.158e+00;SOR=1.034       GT:AD:DP:GQ:PL  0/0:44,0:44:99:0,112,1347       0/1:12,8:20:99:163,0,328        0/0:39,0:39:99:0,105,1194
    ```

이것은 앞서 생성한 VCF와 비슷해 보이지만, 이번에는 세 샘플 모두에 대한 유전형 수준 정보가 있습니다.
파일의 마지막 세 열은 샘플에 대한 유전형 블록이며, 알파벳 순서로 나열되어 있습니다.

첫 번째 변이에 대해 우리 테스트 가족 트리오에 대해 호출된 유전형을 보면, 아버지는 이형접합 변이(`0/1`)이고, 어머니와 아들은 모두 동형접합 변이(`1/1`)입니다.

이것이 궁극적으로 우리가 데이터셋에서 추출하려는 정보입니다!

#### 2.3.3. GATK 컨테이너 종료

컨테이너를 종료하려면 `exit`를 입력하세요.

```bash
exit
```

프롬프트가 정상으로 돌아가야 합니다.
이것으로 변이 호출 명령의 수동 테스트가 완료되었습니다.

---

### 핵심 정리

Samtools 인덱싱과 GATK 변이 호출 명령을 각각의 컨테이너에서 테스트하는 방법을 알게 되었습니다. 여기에는 GVCF를 생성하고 여러 샘플에 대해 공동 유전형 분석을 실행하는 방법도 포함됩니다.

### 다음 단계

컨테이너를 사용하여 작업을 실행하는 워크플로우에 동일한 명령을 적용하는 방법을 학습합니다.
