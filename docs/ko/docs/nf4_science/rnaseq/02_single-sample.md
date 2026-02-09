# 파트 2: 단일 샘플 구현

이 과정의 이번 파트에서는 파트 1에서 실행한 모든 명령을 자동화하는 가장 간단한 워크플로우를 작성하겠습니다. 한 번에 하나의 샘플만 처리하는 것을 목표로 합니다.

이 작업을 세 단계로 진행하겠습니다:

1. 초기 QC 단계를 실행하는 단일 단계 워크플로우 작성
2. 어댑터 트리밍 및 트리밍 후 QC 추가
3. 참조 게놈에 대한 정렬 추가

!!! warning "전제 조건"

    이 수업을 시작하기 전에 과정의 파트 1을 완료해야 합니다.
    특히, 섹션 2.1-3을 진행하면 이 수업의 정렬 단계에 필요한 게놈 인덱스 파일(`data/genome_index.tar.gz`)이 생성됩니다.

---

## 1. 초기 QC를 실행하는 단일 단계 워크플로우 작성

단일 엔드 RNAseq 리드가 포함된 FASTQ 파일에서 FastQC 도구를 실행하는 간단한 워크플로우를 작성하는 것부터 시작하겠습니다.

워크플로우의 주요 부분을 개요로 제공하는 워크플로우 파일 `rnaseq.nf`를 제공합니다.

```groovy title="rnaseq.nf" linenums="1"
#!/usr/bin/env nextflow

// Module INCLUDE statements

/*
 * Pipeline parameters
 */

// Primary input

workflow {

    // Create input channel

    // Call processes

}
```

이 워크플로우 코드는 올바르지만 기능적이지 않다는 점을 유념하세요. 실제 워크플로우를 작성하는 데 사용할 골격 역할만 합니다.

### 1.1. 모듈을 저장할 디렉토리 생성

각 프로세스에 대해 단독 모듈을 생성하여 관리 및 재사용을 쉽게 할 것이므로, 모듈을 저장할 디렉토리를 생성하겠습니다.

```bash
mkdir modules
```

### 1.2. QC 메트릭 수집 프로세스를 위한 모듈 생성

`FASTQC` 프로세스를 담을 `modules/fastqc.nf`라는 모듈 파일을 생성하겠습니다:

```bash
touch modules/fastqc.nf
```

코드 편집기에서 파일을 열고 다음 코드를 복사하세요:

```groovy title="modules/fastqc.nf" linenums="1"
#!/usr/bin/env nextflow

process FASTQC {

    container "community.wave.seqera.io/library/trim-galore:0.6.10--1bf8ca4e1967cd18"
    publishDir "results/fastqc", mode: 'symlink'

    input:
    path reads

    output:
    path "${reads.simpleName}_fastqc.zip", emit: zip
    path "${reads.simpleName}_fastqc.html", emit: html

    script:
    """
    fastqc $reads
    """
}
```

이 교육 시리즈의 파트 1과 파트 2에서 학습한 내용에서 모든 부분을 인식할 수 있을 것입니다. 유일하게 주목할 만한 변경 사항은 이번에는 `publishDir` 지시문에 `mode: symlink`를 사용하고 있으며, 매개변수를 사용하여 `publishDir`을 정의하고 있다는 점입니다.

!!! note "참고"

    여기서 사용하는 데이터 파일은 매우 작지만, 유전체학에서는 매우 클 수 있습니다. 교육 환경에서 시연 목적으로 불필요한 파일 복사를 피하기 위해 'symlink' 게시 모드를 사용하고 있습니다. `work` 디렉토리를 정리할 때 결과를 잃게 되므로 최종 워크플로우에서는 이렇게 하지 않아야 합니다.

### 1.3. 워크플로우 파일에 모듈 가져오기

`rnaseq.nf` 파일에 `include { FASTQC } from './modules/fastqc.nf'` 문을 추가하세요:

```groovy title="rnaseq.nf" linenums="3"
// Module INCLUDE statements
include { FASTQC } from './modules/fastqc.nf'
```

### 1.4. 입력 선언 추가

기본값을 가진 입력 매개변수를 선언하세요:

```groovy title="rnaseq.nf" linenums="10"
params {
    // Primary input
    reads: Path = "data/reads/ENCSR000COQ1_1.fastq.gz"
}
```

### 1.5. 워크플로우 블록에 입력 채널 생성

기본 `.fromPath()` 채널 팩토리를 사용하여 입력 채널을 생성하세요:

```groovy title="rnaseq.nf" linenums="13"
workflow {

    // Create input channel from a file path
    read_ch = channel.fromPath(params.reads)

    // Call processes

}
```

### 1.6. 입력 채널에서 `FASTQC` 프로세스 호출

```groovy title="rnaseq.nf" linenums="13"
workflow {

    // Create input channel from a file path
    read_ch = channel.fromPath(params.reads)

    // Initial quality control
    FASTQC(read_ch)

}
```

### 1.7. 워크플로우를 실행하여 작동하는지 테스트

`--reads` 매개변수를 사용하여 명령줄에서 입력을 지정할 수 있지만, 개발 중에는 설정한 테스트 기본값을 사용하면 편리합니다.

```bash
nextflow run rnaseq.nf
```

??? success "명령 출력"

    ```console
    N E X T F L O W   ~  version 24.10.0

    Launching `rnaseq.nf` [fabulous_snyder] DSL2 - revision: 3394c725ee

    executor >  local (1)
    [d6/d94c3a] FASTQC (1) [100%] 1 of 1 ✔
    ```

파트 1을 진행했고 이미 컨테이너를 가져왔다면 매우 빠르게 실행될 것입니다.
건너뛴 경우, Nextflow가 컨테이너를 자동으로 가져옵니다. 별도로 할 일은 없지만 최대 1분 정도 기다려야 할 수 있습니다.

`publishDir` 지시문에 의해 `FASTQC` 프로세스에서 지정한 대로 `results/fastqc` 아래에서 출력을 찾을 수 있습니다.

```bash
ls results/fastqc
```

```console title="Output"
ENCSR000COQ1_1_fastqc.html  ENCSR000COQ1_1_fastqc.zip
```

---

## 2. 어댑터 트리밍 및 트리밍 후 품질 관리 추가

트리밍 자체를 위한 Cutadapt와 트리밍 후 품질 관리를 위한 FastQC를 번들로 제공하는 Trim_Galore 래퍼를 사용하겠습니다.

### 2.1. 트리밍 및 QC 프로세스를 위한 모듈 생성

`TRIM_GALORE` 프로세스를 담을 `modules/trim_galore.nf`라는 모듈 파일을 생성하겠습니다:

```bash
touch modules/trim_galore.nf
```

코드 편집기에서 파일을 열고 다음 코드를 복사하세요:

```groovy title="modules/trim_galore.nf" linenums="1"
#!/usr/bin/env nextflow

process TRIM_GALORE {

    container "community.wave.seqera.io/library/trim-galore:0.6.10--1bf8ca4e1967cd18"
    publishDir "results/trimming", mode: 'symlink'

    input:
    path reads

    output:
    path "${reads.simpleName}_trimmed.fq.gz", emit: trimmed_reads
    path "${reads}_trimming_report.txt", emit: trimming_reports
    path "${reads.simpleName}_trimmed_fastqc.{zip,html}", emit: fastqc_reports

    script:
    """
    trim_galore --fastqc $reads
    """
}
```

### 2.2. 워크플로우 파일에 모듈 가져오기

`rnaseq.nf` 파일에 `include { TRIM_GALORE } from './modules/trim_galore.nf'` 문을 추가하세요:

```groovy title="rnaseq.nf" linenums="3"
// Module INCLUDE statements
include { FASTQC } from './modules/fastqc.nf'
include { TRIM_GALORE } from './modules/trim_galore.nf'
```

### 2.3. 입력 채널에서 프로세스 호출

```groovy title="rnaseq.nf" linenums="14"
workflow {

    // Create input channel from a file path
    read_ch = channel.fromPath(params.reads)

    // Initial quality control
    FASTQC(read_ch)

    // Adapter trimming and post-trimming QC
    TRIM_GALORE(read_ch)
}
```

### 2.4. 워크플로우를 실행하여 작동하는지 테스트

```bash
nextflow run rnaseq.nf
```

??? success "명령 출력"

    ```console
    N E X T F L O W   ~  version 24.10.0

    Launching `rnaseq.nf` [fabulous_snyder] DSL2 - revision: 3394c725ee

    executor >  local (1)
    [d6/d94c3a] FASTQC (1) [100%] 1 of 1 ✔
    [c2/e4a9bb] TRIM_GALORE (1)  [100%] 1 of 1 ✔
    ```

매우 작은 입력 파일에서 실행하고 있으므로 이것도 매우 빠르게 실행될 것입니다.

`publishDir` 지시문에 의해 `TRIM_GALORE` 프로세스에서 지정한 대로 `results/trimming` 아래에서 출력을 찾을 수 있습니다.

```bash
ls results/trimming
```

```console title="Output"
ENCSR000COQ1_1.fastq.gz_trimming_report.txt  ENCSR000COQ1_1_trimmed_fastqc.zip
ENCSR000COQ1_1_trimmed_fastqc.html           ENCSR000COQ1_1_trimmed.fq.gz
```

---

## 3. 참조 게놈에 리드 정렬

마지막으로 Hisat2를 사용하여 게놈 정렬 단계를 실행할 수 있습니다. 이는 FastQC 스타일의 품질 관리 메트릭도 출력합니다.

### 3.1. HiSat2 프로세스를 위한 모듈 생성

`HISAT2_ALIGN` 프로세스를 담을 `modules/hisat2_align.nf`라는 모듈 파일을 생성하겠습니다:

```bash
touch modules/hisat2_align.nf
```

코드 편집기에서 파일을 열고 다음 코드를 복사하세요:

```groovy title="modules/hisat2_align.nf" linenums="1"
#!/usr/bin/env nextflow

process HISAT2_ALIGN {

    container "community.wave.seqera.io/library/hisat2_samtools:5e49f68a37dc010e"
    publishDir "results/align", mode: 'symlink'

    input:
    path reads
    path index_zip

    output:
    path "${reads.simpleName}.bam", emit: bam
    path "${reads.simpleName}.hisat2.log", emit: log

    script:
    """
    tar -xzvf $index_zip
    hisat2 -x ${index_zip.simpleName} -U $reads \
        --new-summary --summary-file ${reads.simpleName}.hisat2.log | \
        samtools view -bS -o ${reads.simpleName}.bam
    """
}
```

### 3.2. 워크플로우 파일에 모듈 가져오기

`rnaseq.nf` 파일에 `include { HISAT2_ALIGN } from './modules/hisat2_align.nf'` 문을 추가하세요:

```groovy title="rnaseq.nf" linenums="3"
// Module INCLUDE statements
include { FASTQC } from './modules/fastqc.nf'
include { TRIM_GALORE } from './modules/trim_galore.nf'
include { HISAT2_ALIGN } from './modules/hisat2_align.nf'
```

### 3.3. 게놈 인덱스를 제공하기 위한 매개변수 선언 추가

기본값을 가진 입력 매개변수를 선언하세요:

```groovy title="rnaseq.nf" linenums="8"
params {
    // Primary input
    reads: Path = "data/reads/ENCSR000COQ1_1.fastq.gz"

    // Reference genome archive
    hisat2_index_zip: Path = "data/genome_index.tar.gz"
}
```

### 3.4. `TRIM_GALORE`에서 출력된 트리밍된 리드에서 `HISAT2_ALIGN` 프로세스 호출

트리밍된 리드는 이전 단계에서 출력된 `TRIM_GALORE.out.trimmed_reads` 채널에 있습니다.

또한 `file (params.hisat2_index_zip)`을 사용하여 Hisat2 도구에 압축된 게놈 인덱스 tarball을 제공합니다.

```groovy title="rnaseq.nf" linenums="16"
workflow {

    // Create input channel from a file path
    read_ch = channel.fromPath(params.reads)

    // Initial quality control
    FASTQC(read_ch)

    // Adapter trimming and post-trimming QC
    TRIM_GALORE(read_ch)

    // Alignment to a reference genome
    HISAT2_ALIGN(TRIM_GALORE.out.trimmed_reads, file (params.hisat2_index_zip))
}
```

### 3.5. 워크플로우를 실행하여 작동하는지 테스트

```bash
nextflow run rnaseq.nf
```

??? success "명령 출력"

    ```console
    N E X T F L O W   ~  version 24.10.0

    Launching `rnaseq.nf` [extravagant_khorana] DSL2 - revision: 701b41bd16

    executor >  local (3)
    [e4/d15ad4] FASTQC (1)       [100%] 1 of 1 ✔
    [c6/12b2be] TRIM_GALORE (1)  [100%] 1 of 1 ✔
    [c6/7a9f13] HISAT2_ALIGN (1) [100%] 1 of 1 ✔
    ```

`publishDir` 지시문에 의해 `HISAT2_ALIGN` 프로세스에서 지정한 대로 `results/align` 아래에서 출력을 찾을 수 있습니다.

```bash
ls results/align
```

```console title="Output"
ENCSR000COQ1_1_trimmed.bam  ENCSR000COQ1_1_trimmed.hisat2.log
```

이것으로 각 샘플에 적용해야 하는 기본 처리가 완료됩니다.

_한 번에 여러 샘플을 허용하도록 워크플로우를 수정한 후 파트 2에서 MultiQC 리포트 집계를 추가하겠습니다._

---

### 핵심 정리

단일 엔드 RNAseq 샘플을 개별적으로 처리하는 모든 핵심 단계를 적용하는 방법을 알게 되었습니다.

### 다음 단계

여러 샘플을 병렬로 처리하고, 모든 샘플의 모든 단계에 걸쳐 QC 리포트를 집계하며, 페어드 엔드 RNAseq 데이터에서 워크플로우를 실행할 수 있도록 워크플로우를 수정하는 방법을 학습하세요.
