# 파트 1: 방법 개요 및 수동 테스트

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } AI 지원 번역 - [자세히 알아보기 및 개선 제안](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

대량 RNAseq 데이터를 처리하고 분석하는 여러 가지 유효한 방법이 있습니다.
이 과정에서는 [Babraham Institute](https://www.babraham.ac.uk/)의 Simon Andrews 박사와 Laura Biggins 박사가 [여기](https://www.bioinformatics.babraham.ac.uk/training/RNASeq_Course/Analysing%20RNA-Seq%20data%20Exercise.pdf)에서 설명한 방법을 따릅니다.

우리의 목표는 다음 처리 단계를 구현하는 워크플로우를 개발하는 것입니다: 대량 RNAseq 샘플의 리드에 대한 초기 품질 관리 실행, 리드에서 어댑터 서열 트리밍, 참조 게놈에 리드 정렬, 그리고 포괄적인 품질 관리(QC) 보고서 생성.

<figure class="excalidraw">
--8<-- "docs/nf4_science/rnaseq/img/preprocess.svg"
</figure>

- **FASTQC:** FastQC를 사용하여 트리밍 전 리드 데이터에 대한 QC 수행
- **TRIM_GALORE:** Trim Galore를 사용하여 어댑터 서열 트리밍 및 트리밍 후 QC 수행 (Cutadapt 및 FastQC 번들)
- **HISAT2_ALIGN:** Hisat2를 사용하여 참조 게놈에 리드 정렬
- **MULTIQC:** MultiQC를 사용하여 포괄적인 QC 보고서 생성

그러나 워크플로우 코드를 작성하기 전에 먼저 일부 테스트 데이터로 명령을 수동으로 시도해 보겠습니다.
필요한 도구가 GitHub Codespaces 환경에 설치되어 있지 않으므로 컨테이너를 통해 사용할 것입니다 ([Hello Containers](../../hello_nextflow/05_hello_containers.md) 참조).

!!! note "참고"

     `nf4-science/rnaseq` 디렉토리에 있는지 확인하십시오. `pwd`를 입력할 때 표시되는 경로의 마지막 부분은 `rnaseq`이어야 합니다.

---

## 1. 초기 QC 및 어댑터 트리밍

`fastqc`와 `trim_galore`가 모두 설치된 컨테이너 이미지를 가져와서 대화형으로 실행하고 예제 데이터 파일 중 하나에서 트리밍 및 QC 명령을 실행할 것입니다.

### 1.1. 컨테이너 가져오기

```bash
docker pull community.wave.seqera.io/library/trim-galore:0.6.10--1bf8ca4e1967cd18
```

시스템이 이미지를 다운로드하면 다음과 같은 콘솔 출력이 표시됩니다:

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

### 1.2. 대화형으로 컨테이너 실행하기

```bash
docker run -it -v ./data:/data community.wave.seqera.io/library/trim-galore:0.6.10--1bf8ca4e1967cd18
```

<!--
??? success "명령 출력"

    ```console

    ```
-->

프롬프트가 `(base) root@b645838b3314:/tmp#`와 같이 변경되며, 이는 이제 컨테이너 내부에 있음을 나타냅니다.

명령의 `-v ./data:/data` 부분은 컨테이너 내부에서 `data/` 디렉토리의 내용에 액세스할 수 있게 합니다.

```bash
ls /data/reads
```

??? success "명령 출력"

    ```console
    ENCSR000COQ1_1.fastq.gz  ENCSR000COQ2_2.fastq.gz  ENCSR000COR2_1.fastq.gz  ENCSR000CPO1_2.fastq.gz
    ENCSR000COQ1_2.fastq.gz  ENCSR000COR1_1.fastq.gz  ENCSR000COR2_2.fastq.gz  ENCSR000CPO2_1.fastq.gz
    ENCSR000COQ2_1.fastq.gz  ENCSR000COR1_2.fastq.gz  ENCSR000CPO1_1.fastq.gz  ENCSR000CPO2_2.fastq.gzO
    ```

### 1.3. 첫 번째 `fastqc` 명령 실행하기

`fastqc`를 실행하여 리드 데이터에 대한 품질 관리 메트릭을 수집해 보겠습니다.

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

<!-- switch to tree -->

```console title="출력"
/data/reads/ENCSR000COQ1_1_fastqc.html  /data/reads/ENCSR000COQ1_1_fastqc.zip
```

### 1.4. `trim_galore`로 어댑터 서열 트리밍하기

이제 Cutadapt 및 FastQC를 번들로 제공하는 `trim_galore`를 실행하여 어댑터 서열을 트리밍하고 트리밍 후 QC 메트릭을 수집해 보겠습니다.

```bash
trim_galore --fastqc /data/reads/ENCSR000COQ1_1.fastq.gz
```

`--fastqc` 플래그는 트리밍이 완료된 후 QC 수집 단계를 자동으로 실행하도록 합니다.

_출력이 매우 상세하므로 다음은 축약된 내용입니다._

??? success "명령 출력"

    ```console
    Multicore support not enabled. Proceeding with single-core trimming.
    Path to Cutadapt set as: 'cutadapt' (default)
    Cutadapt seems to be working fine (tested command 'cutadapt --version')
    Cutadapt version: 4.9
    single-core operation.
    igzip command line interface 2.31.0
    igzip detected. Using igzip for decompressing

    <...>

    Analysis complete for ENCSR000COQ1_1_trimmed.fq.gz
    ```

작업 디렉토리에서 출력 파일을 찾을 수 있습니다:

```bash
ls ENCSR000COQ1_1*
```

```console title="출력"
ENCSR000COQ1_1.fastq.gz_trimming_report.txt  ENCSR000COQ1_1_trimmed_fastqc.html
ENCSR000COQ1_1_trimmed.fq.gz                 ENCSR000COQ1_1_trimmed_fastqc.zip
```

### 1.5. 출력 파일을 컨테이너 외부의 파일시스템으로 이동하기

컨테이너 내부에 남아 있는 모든 것은 향후 작업에 액세스할 수 없으므로 새 디렉토리로 이동하겠습니다.

```bash
mkdir /data/trimmed
mv ENCSR000COQ1_1* /data/trimmed
```

### 1.6. 컨테이너 종료하기

```bash
exit
```

---

## 2. 참조 게놈에 리드 정렬하기

`hisat2`가 설치된 컨테이너 이미지를 가져와서 대화형으로 실행하고 RNAseq 데이터를 참조 게놈에 정렬하는 정렬 명령을 실행할 것입니다.

### 2.1. `hisat2` 컨테이너 가져오기

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

### 2.2. 대화형으로 `hisat2` 컨테이너 실행하기

```bash
docker run -it -v ./data:/data community.wave.seqera.io/library/hisat2_samtools:5e49f68a37dc010e
```

명령은 이전과 동일하며 관련 컨테이너 URI로 교체되었습니다.

### 2.3. Hisat2 게놈 인덱스 파일 생성하기

Hisat2는 게놈 참조가 매우 특정한 형식으로 제공되어야 하며 우리가 제공하는 `genome.fa` FASTA 파일만으로는 사용할 수 없으므로 이 기회에 관련 리소스를 생성하겠습니다.

```bash
hisat2-build /data/genome.fa genome_index
```

출력이 매우 상세하므로 다음은 축약된 내용입니다:

<!-- TODO: switch to full output -->

??? success "명령 출력"

    ```console
    Settings:
      Output files: "genome_index.*.ht2"
    <...>
    Total time for call to driver() for forward index: 00:00:16
    ```

이렇게 하면 작업 디렉토리에서 찾을 수 있는 여러 게놈 인덱스 파일이 생성됩니다.

```bash
ls genome_index.*
```

```console title="출력"
genome_index.1.ht2  genome_index.3.ht2  genome_index.5.ht2  genome_index.7.ht2
genome_index.2.ht2  genome_index.4.ht2  genome_index.6.ht2  genome_index.8.ht2
```

잠시 후에 사용할 것이지만 먼저 이러한 게놈 인덱스 파일로 gzip 압축된 tarball을 생성하겠습니다. 나중에 필요하며 이러한 파일을 생성하는 것은 일반적으로 워크플로우의 일부로 수행하고 싶지 않은 작업입니다.

```bash
tar -czvf /data/genome_index.tar.gz genome_index.*
```

이렇게 하면 게놈 인덱스 파일이 포함된 `genome_index.tar.gz` tarball이 파일시스템의 `data/` 디렉토리에 저장되며, 이 과정의 파트 2에서 유용하게 사용됩니다.

### 2.4. `hisat2` 명령 실행하기

이제 `hisat2`로 정렬 단계를 수행한 다음 출력을 `samtools`로 파이프하여 출력을 BAM 파일로 작성하는 정렬 명령을 실행할 수 있습니다.

리드 데이터 입력은 이전 단계에서 `trim_galore`로 생성한 `/data/trimmed/ENCSR000COQ1_1_trimmed.fq.gz` 파일입니다.

```bash
hisat2 -x genome_index -U /data/trimmed/ENCSR000COQ1_1_trimmed.fq.gz \
    --new-summary --summary-file ENCSR000COQ1_1_trimmed.hisat2.log | \
    samtools view -bS -o ENCSR000COQ1_1_trimmed.bam
```

??? success "명령 출력"

    ```console
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

```console title="출력"
ENCSR000COQ1_1_trimmed.bam  ENCSR000COQ1_1_trimmed.hisat2.log
```

### 2.5. 출력 파일을 컨테이너 외부의 파일시스템으로 이동하기

```bash
mkdir /data/aligned
mv ENCSR000COQ1_1* /data/aligned
```

### 2.6. 컨테이너 종료하기

```bash
exit
```

---

## 3. 포괄적인 QC 보고서 생성하기

`multiqc`가 설치된 컨테이너 이미지를 가져와서 대화형으로 실행하고 before/after FastQC 보고서 파일에서 보고서 생성 명령을 실행할 것입니다.

### 3.1. `multiqc` 컨테이너 가져오기

```bash
docker pull community.wave.seqera.io/library/pip_multiqc:a3c26f6199d64b7c
```

??? success "명령 출력"

    ```console
    ad8f247edb55897c: Pulling from library/pip_multiqc
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
    3f229294c69a: Pull complete
    5a5ad47fd84c: Pull complete
    Digest: sha256:0ebb1d9605395a7df49ad0eb366b21f46afd96a5090376b0d8941cf5294a895a
    Status: Downloaded newer image for community.wave.seqera.io/library/pip_multiqc:a3c26f6199d64b7c
    community.wave.seqera.io/library/pip_multiqc:a3c26f6199d64b7c
    ```

### 3.2. 대화형으로 `multiqc` 컨테이너 실행하기

```bash
docker run -it -v ./data:/data community.wave.seqera.io/library/pip_multiqc:a3c26f6199d64b7c
```

### 3.3. `multiqc` 명령 실행하기

```bash
multiqc /data/reads /data/trimmed /data/aligned -n ENCSR000COQ1_1_QC
```

??? success "명령 출력"

    ```console

    /// MultiQC 🔍 v1.27.1

          file_search | Search path: /data/reads
          file_search | Search path: /data/trimmed
          file_search | Search path: /data/aligned
            searching | ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━ 100% 20/20
                hisat2 | Found 1 reports
              cutadapt | Found 1 reports
                fastqc | Found 1 reports
        write_results | Data        : ENCSR000COQ1_1_QC_data
        write_results | Report      : ENCSR000COQ1_1_QC.html
              multiqc | MultiQC complete
    ```

MultiQC는 호환 가능한 QC 보고서를 찾기 위해 디렉토리를 검색할 수 있으며 찾은 모든 것을 집계합니다.

여기서 도구가 우리가 생성한 세 가지 QC 보고서를 모두 찾았음을 알 수 있습니다: `fastqc`로 수행한 초기 QC, (`trim_galore`를 통해 생성된) `cutadapt`의 트리밍 후 보고서, 그리고 `hisat2`가 생성한 정렬 후 QC입니다.

출력 파일은 다시 작업 디렉토리에 있습니다:

```bash
ls ENCSR000COQ1_1_QC*
```

```console title="출력"
ENCSR000COQ1_1_QC.html

ENCSR000COQ1_1_QC_data:
cutadapt_filtered_reads_plot.txt                     fastqc_top_overrepresented_sequences_table.txt
cutadapt_trimmed_sequences_plot_3_Counts.txt         hisat2_se_plot.txt
cutadapt_trimmed_sequences_plot_3_Obs_Exp.txt        multiqc.log
fastqc-status-check-heatmap.txt                      multiqc_citations.txt
fastqc_adapter_content_plot.txt                      multiqc_cutadapt.txt
fastqc_per_base_n_content_plot.txt                   multiqc_data.json
fastqc_per_base_sequence_quality_plot.txt            multiqc_fastqc.txt
fastqc_per_sequence_gc_content_plot_Counts.txt       multiqc_general_stats.txt
fastqc_per_sequence_gc_content_plot_Percentages.txt  multiqc_hisat2.txt
fastqc_per_sequence_quality_scores_plot.txt          multiqc_software_versions.txt
fastqc_sequence_counts_plot.txt                      multiqc_sources.txt
fastqc_sequence_duplication_levels_plot.txt
```

### 3.4. 출력 파일을 컨테이너 외부의 파일시스템으로 이동하기

```bash
mkdir /data/final_qc
mv ENCSR000COQ1_1_QC** /data/final_qc
```

### 3.5. 컨테이너 종료하기

```bash
exit
```

---

### 핵심 정리

관련 컨테이너에서 모든 개별 명령을 대화형으로 테스트했습니다.

### 다음 단계는?

동일한 명령을 컨테이너를 사용하여 작업을 실행하는 다단계 워크플로우로 적용하는 방법을 배웁니다.
