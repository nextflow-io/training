# 파트 3: 다중 샘플 paired-end 구현

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } AI 지원 번역 - [자세히 알아보기 및 개선 제안](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

이 과정의 마지막 파트에서는 간단한 워크플로우를 한 단계 더 발전시켜 임의의 수의 샘플을 처리할 수 있는 강력한 배치 자동화 도구로 전환하겠습니다.
그리고 그 과정에서 더 최근 연구에서 일반적으로 사용되는 paired-end 데이터를 처리하도록 전환하겠습니다.

이를 세 단계로 진행하겠습니다:

1. 워크플로우가 여러 입력 샘플을 받고 병렬로 실행하도록 만들기
2. 포괄적인 QC 보고서 생성 추가
3. paired-end RNAseq 데이터로 전환

---

## 1. 워크플로우가 여러 입력 샘플을 받고 병렬로 실행하도록 만들기

입력을 관리하는 방식을 변경해야 합니다.

### 1.1. 기본 입력을 단일 파일 대신 파일 경로의 CSV로 변경

`data/` 디렉토리에 샘플 ID와 FASTQ 파일 경로가 포함된 CSV 파일을 제공합니다.
이 CSV 파일에는 헤더 줄이 포함되어 있습니다.
FASTQ 파일 경로는 절대 경로입니다.

```csv title="data/single-end.csv" linenums="1"
sample_id,fastq_path
ENCSR000COQ1,/workspaces/training/nf4-science/rnaseq/data/reads/ENCSR000COQ1_1.fastq.gz
ENCSR000COQ2,/workspaces/training/nf4-science/rnaseq/data/reads/ENCSR000COQ2_1.fastq.gz
ENCSR000COR1,/workspaces/training/nf4-science/rnaseq/data/reads/ENCSR000COR1_1.fastq.gz
ENCSR000COR2,/workspaces/training/nf4-science/rnaseq/data/reads/ENCSR000COR2_1.fastq.gz
ENCSR000CPO1,/workspaces/training/nf4-science/rnaseq/data/reads/ENCSR000CPO1_1.fastq.gz
ENCSR000CPO2,/workspaces/training/nf4-science/rnaseq/data/reads/ENCSR000CPO2_1.fastq.gz
```

기본 입력 매개변수의 이름을 `input_csv`로 바꾸고 기본값을 `single-end.csv` 파일의 경로로 변경하겠습니다.

```groovy title="rnaseq.nf" linenums="13"
params {
    // 기본 입력
    input_csv: Path = "data/single-end.csv"

    // 참조 게놈 아카이브
    hisat2_index_zip: Path = "data/genome_index.tar.gz"
}
```

### 1.2. CSV를 입력으로 처리하도록 입력 채널 팩토리 업데이트

파일 경로 자체가 아니라 파일의 내용을 채널에 로드하려고 하므로 `.splitCsv()` 연산자를 사용하여 CSV 형식을 파싱한 다음 `.map()` 연산자를 사용하여 원하는 특정 정보(FASTQ 파일 경로)를 가져옵니다.

```groovy title="rnaseq.nf" linenums="16"
    // CSV 파일의 내용에서 입력 채널 생성
    read_ch = channel.fromPath(params.input_csv)
        .splitCsv(header:true)
        .map { row -> file(row.fastq_path) }
```

### 1.3. 워크플로우를 실행하여 작동하는지 테스트

```bash
nextflow run rnaseq.nf
```

??? success "명령 출력"

    ```console
    N E X T F L O W   ~  version 24.10.0

    Launching `rnaseq.nf` [golden_curry] DSL2 - revision: 2a5ba5be1e

    executor >  local (18)
    [07/3ff9c5] FASTQC (6)       [100%] 6 of 6 ✔
    [cc/16859f] TRIM_GALORE (6)  [100%] 6 of 6 ✔
    [68/4c27b5] HISAT2_ALIGN (6) [100%] 6 of 6 ✔
    ```

이번에는 제공한 6개의 데이터 파일 각각에 대해 각 단계가 6번씩 실행되는 것을 볼 수 있습니다.

워크플로우가 여러 파일에서 실행되도록 하는 데 필요한 것은 이것이 전부입니다!
Nextflow가 모든 병렬 처리를 자동으로 처리합니다.

---

## 2. 전처리 QC 메트릭을 단일 MultiQC 보고서로 집계

이렇게 하면 많은 QC 보고서가 생성되며, 개별 보고서를 일일이 찾아보고 싶지 않습니다.
이것이 MultiQC 보고서 집계 단계를 추가하기에 완벽한 시점입니다!

### 2.1. QC 집계 프로세스를 위한 모듈 생성

`MULTIQC` 프로세스를 담을 `modules/multiqc.nf`라는 모듈 파일을 만들겠습니다:

```bash
touch modules/multiqc.nf
```

코드 편집기에서 파일을 열고 다음 코드를 복사합니다:

```groovy title="modules/multiqc.nf" linenums="1"
#!/usr/bin/env nextflow

process MULTIQC {

    container "community.wave.seqera.io/library/pip_multiqc:a3c26f6199d64b7c"
    publishDir "results/multiqc", mode: 'symlink'

    input:
    path '*'
    val output_name

    output:
    path "${output_name}.html", emit: report
    path "${output_name}_data", emit: data

    script:
    """
    multiqc . -n ${output_name}.html
    """
}
```

### 2.2. 워크플로우 파일에 모듈 import

`rnaseq.nf` 파일에 `include { MULTIQC } from './modules/multiqc.nf'` 문을 추가합니다:

```groovy title="rnaseq.nf" linenums="3"
// 모듈 INCLUDE 문
include { FASTQC } from './modules/fastqc.nf'
include { TRIM_GALORE } from './modules/trim_galore.nf'
include { HISAT2_ALIGN } from './modules/hisat2_align.nf'
include { MULTIQC } from './modules/multiqc.nf'
```

### 2.3. `report_id` 매개변수를 추가하고 적절한 기본값 지정

```groovy title="rnaseq.nf" linenums="9"
params {
    // 기본 입력
    input_csv: Path = "data/single-end.csv"

    // 참조 게놈 아카이브
    hisat2_index_zip: Path = "data/genome_index.tar.gz"

    // 보고서 ID
    report_id: String = "all_single-end"
}
```

### 2.4. 이전 단계의 출력에 대해 프로세스 호출

이전 단계의 모든 QC 관련 출력을 `MULTIQC` 프로세스에 제공해야 합니다.

이를 위해 여러 채널을 하나로 집계하는 `.mix()` 연산자를 사용하겠습니다.

각각 간단한 `.out` 채널을 가진 A, B, C, D라는 네 개의 프로세스가 있다면 구문은 다음과 같습니다: `A.out.mix( B.out, C.out, D.out )`. 보시다시피, 결합하려는 채널 중 첫 번째 채널에 적용하고(어느 것이든 상관없음) 그 뒤에 나오는 괄호 안에 쉼표로 구분하여 나머지를 모두 추가합니다.

우리 워크플로우의 경우 집계할 출력은 다음과 같습니다:

- `FASTQC.out.zip`
- `FASTQC.out.html`
- `TRIM_GALORE.out.trimming_reports`
- `TRIM_GALORE.out.fastqc_reports`
- `HISAT2_ALIGN.out.log`

따라서 구문 예제는 다음과 같이 됩니다:

```groovy title="MULTIQC 호출에서 .mix() 적용"
        FASTQC.out.zip.mix(
        FASTQC.out.html,
        TRIM_GALORE.out.trimming_reports,
        TRIM_GALORE.out.fastqc_reports,
        HISAT2_ALIGN.out.log
        )
```

이렇게 하면 샘플당 QC 보고서가 수집됩니다.
그러나 모든 샘플에 걸쳐 집계하려면 모든 샘플의 보고서를 단일 `MULTIQC` 호출로 가져오기 위해 `collect()` 연산자를 추가해야 합니다.
그리고 `report_id` 매개변수도 제공해야 합니다.

이렇게 하면 다음과 같이 됩니다:

```groovy title="완성된 MULTIQC 호출" linenums="33"
    // 종합 QC 보고서 생성
    MULTIQC(
        FASTQC.out.zip.mix(
        FASTQC.out.html,
        TRIM_GALORE.out.trimming_reports,
        TRIM_GALORE.out.fastqc_reports,
        HISAT2_ALIGN.out.log
        ).collect(),
        params.report_id
    )
```

전체 워크플로우 블록의 컨텍스트에서는 다음과 같이 보입니다:

```groovy title="rnaseq.nf" linenums="18"
workflow {
    // CSV 파일의 내용에서 입력 채널 생성
    read_ch = channel.fromPath(params.input_csv)
        .splitCsv(header:true)
        .map { row -> file(row.fastq_path) }

    /// 초기 품질 관리
    FASTQC(read_ch)

    // 어댑터 트리밍 및 트리밍 후 QC
    TRIM_GALORE(read_ch)

    // 참조 게놈에 정렬
    HISAT2_ALIGN(TRIM_GALORE.out.trimmed_reads, file (params.hisat2_index_zip))

    // 종합 QC 보고서 생성
    MULTIQC(
        FASTQC.out.zip.mix(
        FASTQC.out.html,
        TRIM_GALORE.out.trimming_reports,
        TRIM_GALORE.out.fastqc_reports,
        HISAT2_ALIGN.out.log
        ).collect(),
        params.report_id
    )
}
```

### 2.5. 워크플로우를 실행하여 작동하는지 테스트

```bash
nextflow run rnaseq.nf -resume
```

??? success "명령 출력"

    ```console
    N E X T F L O W   ~  version 24.10.0

    Launching `rnaseq.nf` [modest_pare] DSL2 - revision: fc724d3b49

    executor >  local (1)
    [07/3ff9c5] FASTQC (6)       [100%] 6 of 6, cached: 6 ✔
    [2c/8d8e1e] TRIM_GALORE (5)  [100%] 6 of 6, cached: 6 ✔
    [a4/7f9c44] HISAT2_ALIGN (6) [100%] 6 of 6, cached: 6 ✔
    [56/e1f102] MULTIQC          [100%] 1 of 1 ✔
    ```

이번에는 캐시된 프로세스 호출 후에 단일 MULTIQC 호출이 추가된 것을 볼 수 있습니다:

`TRIM_GALORE` 프로세스의 `publishDir` 지시문에 지정된 대로 `results/trimming` 아래에서 출력을 찾을 수 있습니다.

```bash
tree -L 2 results/multiqc
```

```console title="출력"
results/multiqc
├── all_single-end_data
│   ├── cutadapt_filtered_reads_plot.txt
│   ├── cutadapt_trimmed_sequences_plot_3_Counts.txt
│   ├── cutadapt_trimmed_sequences_plot_3_Obs_Exp.txt
│   ├── fastqc_adapter_content_plot.txt
│   ├── fastqc_overrepresented_sequences_plot.txt
│   ├── fastqc_per_base_n_content_plot.txt
│   ├── fastqc_per_base_sequence_quality_plot.txt
│   ├── fastqc_per_sequence_gc_content_plot_Counts.txt
│   ├── fastqc_per_sequence_gc_content_plot_Percentages.txt
│   ├── fastqc_per_sequence_quality_scores_plot.txt
│   ├── fastqc_sequence_counts_plot.txt
│   ├── fastqc_sequence_duplication_levels_plot.txt
│   ├── fastqc_sequence_length_distribution_plot.txt
│   ├── fastqc-status-check-heatmap.txt
│   ├── fastqc_top_overrepresented_sequences_table.txt
│   ├── hisat2_se_plot.txt
│   ├── multiqc_citations.txt
│   ├── multiqc_cutadapt.txt
│   ├── multiqc_data.json
│   ├── multiqc_fastqc.txt
│   ├── multiqc_general_stats.txt
│   ├── multiqc_hisat2.txt
│   ├── multiqc.log
│   ├── multiqc_software_versions.txt
│   └── multiqc_sources.txt
└── all_single-end.html
```

마지막 `all_single-end.html` 파일은 하나의 쉽게 탐색할 수 있는 HTML 파일로 편리하게 패키징된 전체 집계 보고서입니다.

---

## 3. paired-end RNAseq 데이터 처리 활성화

현재 우리의 워크플로우는 single-end RNAseq 데이터만 처리할 수 있습니다.
paired-end RNAseq 데이터를 보는 것이 점점 일반적이므로 이를 처리할 수 있기를 원합니다.

워크플로우를 데이터 타입에 완전히 독립적으로 만들려면 약간 더 고급 Nextflow 언어 기능을 사용해야 하므로 여기서는 그렇게 하지 않겠지만, 무엇을 조정해야 하는지 보여주기 위해 paired-end 처리 버전을 만들 수 있습니다.

### 3.1. `rnaseq_pe.nf`라는 워크플로우 복사본 만들기

```bash
cp rnaseq.nf rnaseq_pe.nf
```

### 3.2. 기본 `input_csv`를 paired-end 데이터를 가리키도록 수정

`data/` 디렉토리에 샘플 ID와 paired FASTQ 파일 경로가 포함된 두 번째 CSV 파일을 제공합니다

```csv title="data/paired-end.csv" linenums="1"
sample_id,fastq_1,fastq_2
ENCSR000COQ1,/workspaces/training/nf4-science/rnaseq/data/reads/ENCSR000COQ1_1.fastq.gz,/workspaces/training/nf4-science/rnaseq/data/reads/ENCSR000COQ1_2.fastq.gz
ENCSR000COQ2,/workspaces/training/nf4-science/rnaseq/data/reads/ENCSR000COQ2_1.fastq.gz,/workspaces/training/nf4-science/rnaseq/data/reads/ENCSR000COQ2_2.fastq.gz
ENCSR000COR1,/workspaces/training/nf4-science/rnaseq/data/reads/ENCSR000COR1_1.fastq.gz,/workspaces/training/nf4-science/rnaseq/data/reads/ENCSR000COR1_2.fastq.gz
ENCSR000COR2,/workspaces/training/nf4-science/rnaseq/data/reads/ENCSR000COR2_1.fastq.gz,/workspaces/training/nf4-science/rnaseq/data/reads/ENCSR000COR2_2.fastq.gz
ENCSR000CPO1,/workspaces/training/nf4-science/rnaseq/data/reads/ENCSR000CPO1_1.fastq.gz,/workspaces/training/nf4-science/rnaseq/data/reads/ENCSR000CPO1_2.fastq.gz
ENCSR000CPO2,/workspaces/training/nf4-science/rnaseq/data/reads/ENCSR000CPO2_1.fastq.gz,/workspaces/training/nf4-science/rnaseq/data/reads/ENCSR000CPO2_2.fastq.gz
```

`input_csv` 기본값을 `paired-end.csv` 파일의 경로로 변경하겠습니다.

```groovy title="rnaseq_pe.nf" linenums="15"
params {
    // 기본 입력
    input_csv: Path = "data/paired-end.csv"

    // 참조 게놈 아카이브
    hisat2_index_zip: Path = "data/genome_index.tar.gz"

    // 보고서 ID
    report_id: String = "all_single-end"
}
```

### 3.3. 채널 팩토리 업데이트

이제 두 FASTQ 파일 경로를 모두 가져오도록 `.map()` 연산자에 지시해야 합니다.

따라서 `row -> file(row.fastq_path)`는 `row -> [file(row.fastq_1), file(row.fastq_2)]`가 됩니다

```groovy title="rnaseq_pe.nf" linenums="19"
    // CSV 파일의 내용에서 입력 채널 생성
    read_ch = channel.fromPath(params.input_csv)
        .splitCsv(header:true)
        .map { row -> [file(row.fastq_1), file(row.fastq_2)] }
```

### 3.4. FASTQC 프로세스의 paired-end 버전 만들기

두 버전을 모두 사용할 수 있도록 모듈 복사본을 만들겠습니다.

```bash
cp modules/fastqc.nf modules/fastqc_pe.nf
```

코드 편집기에서 새 `fastqc_pe.nf` 모듈 파일을 열고 다음 코드 변경 사항을 적용합니다:

- `script` 블록(17번째 줄)에서 `fastqc $reads`를 `fastqc ${reads}`로 변경하여 `reads` 입력이 이제 단일 경로가 아니라 두 경로의 튜플이므로 압축이 해제되도록 합니다.
- 출력 파일을 개별적으로 처리하지 않도록 `${reads.simpleName}`을 와일드카드(`*`)로 바꿉니다.

```groovy title="modules/fastqc_pe.nf" linenums="8"
    input:
    path reads

    output:
    path "*_fastqc.zip", emit: zip
    path "*_fastqc.html", emit: html

    script:
    """
    fastqc ${reads}
    """
```

기술적으로 이것은 `FASTQC` 프로세스를 single-end 또는 paired-end RNAseq 데이터를 처리할 수 있도록 일반화합니다.

마지막으로, 모듈의 paired-end 버전을 사용하도록 모듈 import 문을 업데이트합니다.

```groovy title="rnaseq_pe.nf" linenums="4"
include { FASTQC } from './modules/fastqc_pe.nf'
```

### 3.5. TRIM_GALORE 프로세스의 paired-end 버전 만들기

두 버전을 모두 사용할 수 있도록 모듈 복사본을 만듭니다.

```bash
cp modules/trim_galore.nf modules/trim_galore_pe.nf
```

코드 편집기에서 새 `trim_galore_pe.nf` 모듈 파일을 열고 다음 코드 변경 사항을 적용합니다:

- 입력 선언을 `path reads`에서 `tuple path(read1), path(read2)`로 변경
- `script` 블록의 명령을 업데이트하여 `$reads`를 `--paired ${read1} ${read2}`로 바꿉니다
- 추가된 파일과 다른 명명 규칙을 반영하도록 출력 선언을 업데이트하고, 모든 것을 나열하지 않도록 와일드카드를 사용합니다.

```groovy title="modules/trim_galore_pe.nf" linenums="8"
    input:
    tuple path(read1), path(read2)

    output:
    tuple path("*_val_1.fq.gz"), path("*_val_2.fq.gz"), emit: trimmed_reads
    path "*_trimming_report.txt", emit: trimming_reports
    path "*_val_1_fastqc.{zip,html}", emit: fastqc_reports_1
    path "*_val_2_fastqc.{zip,html}", emit: fastqc_reports_2

    script:
    """
    trim_galore --fastqc --paired ${read1} ${read2}
    """
```

마지막으로, 모듈의 paired-end 버전을 사용하도록 모듈 import 문을 업데이트합니다.

```groovy title="rnaseq_pe.nf" linenums="5"
include { TRIM_GALORE } from './modules/trim_galore_pe.nf'
```

### 3.6. TRIM_GALORE에서 두 개의 보고서를 받을 것으로 예상하도록 MULTIQC 프로세스 호출 업데이트

`TRIM_GALORE` 프로세스가 이제 추가 출력 채널을 생성하므로 이를 MultiQC에 공급해야 합니다.

`TRIM_GALORE.out.fastqc_reports,`를 `TRIM_GALORE.out.fastqc_reports_1,`과 `TRIM_GALORE.out.fastqc_reports_2,`로 바꿉니다:

```groovy title="rnaseq_pe.nf" linenums="33"
    // 종합 QC 보고서 생성
    MULTIQC(
        FASTQC.out.zip.mix(
        FASTQC.out.html,
        TRIM_GALORE.out.trimming_reports,
        TRIM_GALORE.out.fastqc_reports_1,
        TRIM_GALORE.out.fastqc_reports_2,
        HISAT2_ALIGN.out.log
        ).collect(),
        params.report_id
    )
```

MultiQC에 대해 이야기하는 김에 `report_id` 매개변수 기본값도 `"all_single-end"`에서 `"all_paired-end"`로 업데이트하겠습니다.

```groovy title="rnaseq_pe.nf" linenums="9"
params {
    // 기본 입력
    input_csv: Path = "data/paired-end.csv"

    // 참조 게놈 아카이브
    hisat2_index_zip: Path = "data/genome_index.tar.gz"

    // 보고서 ID
    report_id: String = "all_paired-end"
}
```

### 3.7. HISAT2_ALIGN 프로세스의 paired-end 버전 만들기

두 버전을 모두 사용할 수 있도록 모듈 복사본을 만듭니다.

```bash
cp modules/hisat2_align.nf modules/hisat2_align_pe.nf
```

코드 편집기에서 새 `hisat2_align_pe.nf` 모듈 파일을 열고 다음 코드 변경 사항을 적용합니다:

- 입력 선언을 `path reads`에서 `tuple path(read1), path(read2)`로 변경
- `script` 블록의 명령을 업데이트하여 `-U $reads`를 `-1 ${read1} -2 ${read2}`로 바꿉니다
- `script` 블록의 명령과 출력 선언에서 `${reads.simpleName}`의 모든 인스턴스를 `${read1.simpleName}`으로 바꿉니다.

```groovy title="modules/hisat2_align_pe.nf" linenums="8"
    input:
    tuple path(read1), path(read2)
    path index_zip

    output:
    path "${read1.simpleName}.bam", emit: bam
    path "${read1.simpleName}.hisat2.log", emit: log

    script:
    """
    tar -xzvf $index_zip
    hisat2 -x ${index_zip.simpleName} -1 ${read1} -2 ${read2} \
        --new-summary --summary-file ${read1.simpleName}.hisat2.log | \
        samtools view -bS -o ${read1.simpleName}.bam
    """
```

마지막으로, 모듈의 paired-end 버전을 사용하도록 모듈 import 문을 업데이트합니다.

```groovy title="rnaseq_pe.nf" linenums="5"
include { HISAT2_ALIGN } from './modules/hisat2_align_pe.nf'
```

### 3.8. 워크플로우를 실행하여 작동하는지 테스트

캐시되지 않으므로 `-resume`을 사용하지 않으며, 이전보다 처리할 데이터가 두 배 더 많지만 1분 이내에 완료되어야 합니다.

```bash
nextflow run rnaseq_pe.nf
```

??? success "명령 출력"

    ```console
    N E X T F L O W   ~  version 24.10.0

    Launching `rnaseq_pe.nf` [reverent_kare] DSL2 - revision: 9c376cc219

    executor >  local (19)
    [c5/cbde15] FASTQC (5)       [100%] 6 of 6 ✔
    [e4/fa2784] TRIM_GALORE (5)  [100%] 6 of 6 ✔
    [3a/e23049] HISAT2_ALIGN (5) [100%] 6 of 6 ✔
    [e6/a3ccd9] MULTIQC          [100%] 1 of 1 ✔
    ```

이것으로 끝입니다! 이제 워크플로우의 약간 다른 두 버전이 있습니다. 하나는 single-end 읽기 데이터용이고 다른 하나는 paired-end 데이터용입니다.
다음 논리적 단계는 워크플로우가 즉시 두 데이터 타입을 모두 받을 수 있도록 만드는 것이며, 이는 이 과정의 범위를 벗어나지만 후속 과정에서 다룰 수 있습니다.

---

### 요점 정리

단일 샘플 워크플로우를 여러 샘플의 병렬 처리로 조정하고, 포괄적인 QC 보고서를 생성하며, 필요한 경우 paired-end 읽기 데이터를 사용하도록 워크플로우를 조정하는 방법을 알게 되었습니다.

### 다음 단계는?

축하합니다. Nextflow For RNAseq 단기 과정을 완료하셨습니다! 성공을 축하하고 충분한 휴식을 취하세요!

다음으로, 이 교육 과정에 대한 귀하의 경험에 대한 매우 짧은 설문조사를 완료해 주시기 바라며, 그런 다음 추가 교육 리소스 및 유용한 링크가 포함된 페이지로 안내해 드리겠습니다.
