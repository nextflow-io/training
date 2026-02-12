# 파트 3: 다중 샘플 paired-end 구현

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } AI 지원 번역 - [자세히 알아보기 및 개선 제안](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

이전에는 각 샘플의 데이터를 독립적으로 처리하는 샘플별 변이 호출 파이프라인을 구축했습니다.
이 과정의 이번 파트에서는 간단한 워크플로우를 한 단계 더 발전시켜 임의의 수의 샘플을 처리할 수 있는 강력한 배치 자동화 도구로 전환하겠습니다.
그리고 그 과정에서 더 최근 연구에서 일반적으로 사용되는 paired-end 데이터를 처리하도록 업데이트하겠습니다.

??? info "이 섹션을 시작하는 방법"

    이 과정의 이번 섹션은 [파트 1: 방법 개요](./01_method.md), [파트 2: 단일 샘플 구현](./02_single-sample.md)을 완료하고 모듈 파일이 채워진 작동하는 `rnaseq.nf` 파이프라인이 있다고 가정합니다.

    파트 2를 완료하지 않았거나 이번 파트를 새로 시작하고 싶다면 파트 2 해결책을 시작점으로 사용할 수 있습니다.
    `nf4-science/rnaseq/` 디렉토리 내에서 다음 명령을 실행하세요:

    ```bash
    cp solutions/part2/rnaseq-2.nf rnaseq.nf
    cp solutions/part2/modules/fastqc.nf modules/
    cp solutions/part2/modules/trim_galore.nf modules/
    cp solutions/part2/modules/hisat2_align.nf modules/
    cp solutions/part2/nextflow.config .
    ```

    이렇게 하면 완전한 단일 샘플 처리 워크플로우를 얻을 수 있습니다.
    성공적으로 실행되는지 테스트할 수 있습니다:

    ```bash
    nextflow run rnaseq.nf -profile test
    ```

## 과제

이 과정의 이번 파트에서는 워크플로우를 다음과 같이 확장하겠습니다:

1. CSV 샘플시트에서 샘플 정보 읽기
2. 모든 샘플에 대해 샘플별 QC, 트리밍, 정렬을 병렬로 실행
3. 모든 QC 보고서를 포괄적인 MultiQC 보고서로 집계

<figure class="excalidraw">
--8<-- "docs/en/docs/nf4_science/rnaseq/img/rnaseq-wf-03.svg"
</figure>

이는 [파트 1: 방법 개요](./01_method.md#2-multi-sample-qc-aggregation)의 두 번째 섹션에서 컨테이너에서 수동으로 실행했던 단계를 자동화합니다.

## 학습 계획

이를 세 단계로 나누었습니다:

1. **워크플로우가 여러 입력 샘플을 받도록 만들기.**
   단일 파일 경로에서 CSV 샘플시트로 전환하고, `splitCsv()`로 파싱하며, 모든 기존 프로세스를 여러 샘플에서 실행하는 것을 다룹니다.
2. **포괄적인 QC 보고서 생성 추가.**
   샘플 전체에 걸쳐 출력을 집계하는 `collect()` 연산자를 소개하고, 통합 보고서를 생성하는 MultiQC 프로세스를 추가합니다.
3. **paired-end RNAseq 데이터로 전환.**
   paired-end 입력을 위한 프로세스 조정(튜플 사용), paired-end 모듈 생성, 별도의 테스트 프로파일 설정을 다룹니다.

이는 [파트 1: 방법 개요](./01_method.md)에 설명된 방법(다중 샘플 사용 사례를 다루는 두 번째 섹션)을 구현하고 파트 2에서 생성된 워크플로우를 직접 기반으로 합니다.

!!! tip "팁"

     올바른 작업 디렉토리에 있는지 확인하세요:
     `cd /workspaces/training/nf4-science/rnaseq`

---

## 1. 워크플로우가 여러 입력 샘플을 받도록 만들기

여러 샘플에서 실행하려면 입력을 관리하는 방식을 변경해야 합니다. 단일 파일 경로를 제공하는 대신 CSV 파일에서 샘플 정보를 읽습니다.

`data/` 디렉토리에 샘플 ID와 FASTQ 파일 경로가 포함된 CSV 파일을 제공합니다.

```csv title="data/single-end.csv" linenums="1"
sample_id,fastq_path
ENCSR000COQ1,/workspaces/training/nf4-science/rnaseq/data/reads/ENCSR000COQ1_1.fastq.gz
ENCSR000COQ2,/workspaces/training/nf4-science/rnaseq/data/reads/ENCSR000COQ2_1.fastq.gz
ENCSR000COR1,/workspaces/training/nf4-science/rnaseq/data/reads/ENCSR000COR1_1.fastq.gz
ENCSR000COR2,/workspaces/training/nf4-science/rnaseq/data/reads/ENCSR000COR2_1.fastq.gz
ENCSR000CPO1,/workspaces/training/nf4-science/rnaseq/data/reads/ENCSR000CPO1_1.fastq.gz
ENCSR000CPO2,/workspaces/training/nf4-science/rnaseq/data/reads/ENCSR000CPO2_1.fastq.gz
```

이 CSV 파일에는 열 이름을 지정하는 헤더 줄이 포함되어 있습니다.

이것은 여전히 single-end 읽기 데이터입니다.

!!! warning "경고"

    CSV의 파일 경로는 환경과 일치해야 하는 절대 경로입니다.
    제공하는 교육 환경에서 실행하지 않는 경우 시스템에 맞게 경로를 업데이트해야 합니다.

### 1.1. 테스트 프로파일에서 기본 입력을 파일 경로의 CSV로 변경

먼저 `nextflow.config`의 테스트 프로파일을 업데이트하여 단일 FASTQ 경로 대신 CSV 파일 경로를 제공해야 합니다.

=== "후"

    ```groovy title="nextflow.config" linenums="1" hl_lines="5"
    docker.enabled = true

    profiles {
        test {
            params.input = "${projectDir}/data/single-end.csv"
            params.hisat2_index_zip = "${projectDir}/data/genome_index.tar.gz"
        }
    }
    ```

=== "전"

    ```groovy title="nextflow.config" linenums="1"
    docker.enabled = true

    profiles {
        test {
            params.input = "${projectDir}/data/reads/ENCSR000COQ1_1.fastq.gz"
            params.hisat2_index_zip = "${projectDir}/data/genome_index.tar.gz"
        }
    }
    ```

다음으로 이 CSV에서 읽도록 채널 생성을 업데이트해야 합니다.

### 1.2. CSV 입력을 파싱하도록 채널 팩토리 업데이트

파일 경로 자체가 아니라 파일의 내용을 채널에 로드해야 합니다.

[Hello Nextflow의 파트 2](../../hello_nextflow/02_hello_channels.md#42-use-the-splitcsv-operator-to-parse-the-file)에서 사용한 것과 동일한 패턴을 사용하여 이를 수행할 수 있습니다. 파일을 파싱하기 위해 [`splitCsv()`](https://nextflow.io/docs/latest/reference/operator.html#splitcsv) 연산자를 적용한 다음 `map` 연산을 사용하여 각 행에서 FASTQ 파일 경로를 추출합니다.

=== "후"

    ```groovy title="rnaseq.nf" linenums="19" hl_lines="4-7"
    workflow {

        main:
        // CSV 파일의 내용에서 입력 채널 생성
        read_ch = channel.fromPath(params.input)
            .splitCsv(header: true)
            .map { row -> file(row.fastq_path) }
    ```

=== "전"

    ```groovy title="rnaseq.nf" linenums="19"
    workflow {

        main:
        // 파일 경로에서 입력 채널 생성
        read_ch = channel.fromPath(params.input)
    ```

Hello Nextflow 과정에서 접한 것과 비교하여 새로운 점은 이 CSV에 헤더 줄이 있다는 것이므로 `splitCsv()` 호출에 `#!groovy header: true`를 추가합니다.
이를 통해 `map` 연산에서 이름으로 열을 참조할 수 있습니다. `#!groovy row.fastq_path`는 각 행의 `fastq_path` 열에서 파일 경로를 추출합니다.

입력 처리가 업데이트되었고 워크플로우를 테스트할 준비가 되었습니다.

### 1.3. 워크플로우 실행

이제 워크플로우는 CSV 파일에서 샘플 정보를 읽고 모든 샘플을 병렬로 처리합니다.

```bash
nextflow run rnaseq.nf -profile test
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

이번에는 CSV 파일의 각 샘플에 대해 각 단계가 6번씩 실행됩니다.

워크플로우가 여러 파일에서 실행되도록 하는 데 필요한 것은 이것이 전부입니다.
Nextflow가 모든 병렬 처리를 자동으로 처리합니다.

### 핵심 정리

단일 파일 입력에서 Nextflow가 병렬로 처리하는 CSV 기반 다중 샘플 입력으로 전환하는 방법을 알게 되었습니다.

### 다음 단계는?

모든 샘플의 메트릭을 결합하는 QC 보고서 집계 단계를 추가합니다.

---

## 2. 전처리 QC 메트릭을 단일 MultiQC 보고서로 집계

이렇게 하면 많은 QC 보고서가 생성되며, 개별 보고서를 일일이 찾아보고 싶지 않습니다.
이것이 MultiQC 보고서 집계 단계를 추가하기에 완벽한 시점입니다.

[파트 1](01_method.md)의 `multiqc` 명령을 기억하세요:

```bash
multiqc . -n <output_name>.html
```

이 명령은 현재 디렉토리에서 인식된 QC 출력 파일을 스캔하고 이를 단일 HTML 보고서로 집계합니다.
컨테이너 URI는 `community.wave.seqera.io/library/pip_multiqc:a3c26f6199d64b7c`였습니다.

추가 매개변수를 설정하고, 입력을 준비하고, 프로세스를 작성하고, 연결하고, 출력 처리를 업데이트해야 합니다.

### 2.1. 입력 설정

MultiQC 프로세스에는 보고서 이름 매개변수와 이전 모든 단계에서 수집된 QC 출력이 함께 번들로 필요합니다.

#### 2.1.1. `report_id` 매개변수 추가

출력 보고서의 이름을 지정하는 매개변수를 추가합니다.

=== "후"

    ```groovy title="rnaseq.nf" linenums="11" hl_lines="8-9"
    params {
        // 기본 입력
        input: Path

        // 참조 게놈 아카이브
        hisat2_index_zip: Path

        // 보고서 ID
        report_id: String
    }
    ```

=== "전"

    ```groovy title="rnaseq.nf" linenums="11"
    params {
        // 기본 입력
        input: Path

        // 참조 게놈 아카이브
        hisat2_index_zip: Path
    }
    ```

테스트 프로파일에 보고서 ID 기본값을 추가합니다:

=== "후"

    ```groovy title="nextflow.config" linenums="1" hl_lines="7"
    docker.enabled = true

    profiles {
        test {
            params.input = "${projectDir}/data/single-end.csv"
            params.hisat2_index_zip = "${projectDir}/data/genome_index.tar.gz"
            params.report_id = "all_single-end"
        }
    }
    ```

=== "전"

    ```groovy title="nextflow.config" linenums="1"
    docker.enabled = true

    profiles {
        test {
            params.input = "${projectDir}/data/single-end.csv"
            params.hisat2_index_zip = "${projectDir}/data/genome_index.tar.gz"
        }
    }
    ```

다음으로 MultiQC 프로세스의 입력을 준비해야 합니다.

#### 2.1.2. 이전 단계의 QC 출력 수집 및 결합

`MULTIQC` 프로세스에 이전 단계의 모든 QC 관련 출력을 함께 번들로 제공해야 합니다.

이를 위해 여러 채널을 하나로 집계하는 `.mix()` 연산자를 사용합니다.
`channel.empty()`에서 시작하여 결합하려는 모든 출력 채널을 혼합합니다.
이는 출력 채널 중 하나에 직접 `.mix()`를 연결하는 것보다 깔끔합니다. 모든 입력을 대칭적으로 처리하기 때문입니다.

우리 워크플로우에서 집계할 QC 관련 출력은 다음과 같습니다:

- `FASTQC.out.zip`
- `FASTQC.out.html`
- `TRIM_GALORE.out.trimming_reports`
- `TRIM_GALORE.out.fastqc_reports`
- `HISAT2_ALIGN.out.log`

이들을 단일 채널로 혼합한 다음 `.collect()`를 사용하여 모든 샘플의 보고서를 단일 리스트로 집계합니다.

`HISAT2_ALIGN` 호출 후 워크플로우 본문에 다음 줄을 추가합니다:

=== "후"

    ```groovy title="rnaseq.nf" linenums="38" hl_lines="4-13"
        // 참조 게놈에 정렬
        HISAT2_ALIGN(TRIM_GALORE.out.trimmed_reads, file(params.hisat2_index_zip))

        // 종합 QC 보고서 생성
        multiqc_files_ch = channel.empty().mix(
            FASTQC.out.zip,
            FASTQC.out.html,
            TRIM_GALORE.out.trimming_reports,
            TRIM_GALORE.out.fastqc_reports,
            HISAT2_ALIGN.out.log,
        )
        multiqc_files_list = multiqc_files_ch.collect()
    ```

=== "전"

    ```groovy title="rnaseq.nf" linenums="38"
        // 참조 게놈에 정렬
        HISAT2_ALIGN(TRIM_GALORE.out.trimmed_reads, file(params.hisat2_index_zip))
    ```

중간 변수를 사용하면 각 단계가 명확해집니다. `multiqc_files_ch`는 하나의 채널로 혼합된 모든 개별 QC 파일을 포함하고, `multiqc_files_list`는 MultiQC에 전달할 준비가 된 수집된 번들입니다.

### 2.2. QC 집계 프로세스 작성 및 워크플로우에서 호출

이전과 마찬가지로 프로세스 정의를 채우고, 모듈을 import하고, 프로세스 호출을 추가해야 합니다.

#### 2.2.1. QC 집계 프로세스의 모듈 채우기

`modules/multiqc.nf`를 열고 프로세스 정의의 개요를 검토합니다.

위에 제공된 정보를 사용하여 직접 프로세스 정의를 채운 다음 아래 "후" 탭의 해결책과 비교하여 작업을 확인하세요.

=== "전"

    ```groovy title="modules/multiqc.nf" linenums="1"
    #!/usr/bin/env nextflow

    /*
     * MultiQC로 QC 보고서 집계
     */
    process MULTIQC {

        container

        input:

        output:

        script:
        """

        """
    }
    ```

=== "후"

    ```groovy title="modules/multiqc.nf" linenums="1" hl_lines="8 11 12 15 16 20"
    #!/usr/bin/env nextflow

    /*
     * MultiQC로 QC 보고서 집계
     */
    process MULTIQC {

        container "community.wave.seqera.io/library/pip_multiqc:a3c26f6199d64b7c"

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

이 프로세스는 QC 파일의 입력 한정자로 `#!groovy path '*'`를 사용합니다.
`'*'` 와일드카드는 Nextflow에게 특정 이름을 요구하지 않고 수집된 모든 파일을 작업 디렉토리에 스테이징하도록 지시합니다.
`val output_name` 입력은 보고서 파일 이름을 제어하는 문자열입니다.

`multiqc .` 명령은 현재 디렉토리(모든 스테이징된 QC 파일이 있는 곳)를 스캔하고 보고서를 생성합니다.

이를 완료하면 프로세스를 사용할 준비가 됩니다.

#### 2.2.2. 모듈 include

`rnaseq.nf`에 import 문을 추가합니다:

=== "후"

    ```groovy title="rnaseq.nf" linenums="3" hl_lines="5"
    // 모듈 INCLUDE 문
    include { FASTQC } from './modules/fastqc.nf'
    include { TRIM_GALORE } from './modules/trim_galore.nf'
    include { HISAT2_ALIGN } from './modules/hisat2_align.nf'
    include { MULTIQC } from './modules/multiqc.nf'
    ```

=== "전"

    ```groovy title="rnaseq.nf" linenums="3"
    // 모듈 INCLUDE 문
    include { FASTQC } from './modules/fastqc.nf'
    include { TRIM_GALORE } from './modules/trim_galore.nf'
    include { HISAT2_ALIGN } from './modules/hisat2_align.nf'
    ```

이제 워크플로우에 프로세스 호출을 추가합니다.

#### 2.2.3. 프로세스 호출 추가

수집된 QC 파일과 보고서 ID를 `MULTIQC` 프로세스에 전달합니다:

=== "후"

    ```groovy title="rnaseq.nf" linenums="49" hl_lines="2"
        multiqc_files_list = multiqc_files_ch.collect()
        MULTIQC(multiqc_files_list, params.report_id)
    ```

=== "전"

    ```groovy title="rnaseq.nf" linenums="49"
        multiqc_files_list = multiqc_files_ch.collect()
    ```

MultiQC 프로세스가 이제 워크플로우에 연결되었습니다.

### 2.3. 출력 처리 업데이트

publish 선언에 MultiQC 출력을 추가하고 어디로 갈지 구성해야 합니다.

#### 2.3.1. MultiQC 출력에 대한 publish 대상 추가

`publish:` 섹션에 MultiQC 출력을 추가합니다:

=== "후"

    ```groovy title="rnaseq.nf" linenums="52" hl_lines="9-10"
        publish:
        fastqc_zip = FASTQC.out.zip
        fastqc_html = FASTQC.out.html
        trimmed_reads = TRIM_GALORE.out.trimmed_reads
        trimming_reports = TRIM_GALORE.out.trimming_reports
        trimming_fastqc = TRIM_GALORE.out.fastqc_reports
        bam = HISAT2_ALIGN.out.bam
        align_log = HISAT2_ALIGN.out.log
        multiqc_report = MULTIQC.out.report
        multiqc_data = MULTIQC.out.data
    }
    ```

=== "전"

    ```groovy title="rnaseq.nf" linenums="52"
        publish:
        fastqc_zip = FASTQC.out.zip
        fastqc_html = FASTQC.out.html
        trimmed_reads = TRIM_GALORE.out.trimmed_reads
        trimming_reports = TRIM_GALORE.out.trimming_reports
        trimming_fastqc = TRIM_GALORE.out.fastqc_reports
        bam = HISAT2_ALIGN.out.bam
        align_log = HISAT2_ALIGN.out.log
    }
    ```

다음으로 Nextflow에게 이러한 출력을 어디에 둘지 알려야 합니다.

#### 2.3.2. 새 출력 대상 구성

`output {}` 블록에 MultiQC 대상에 대한 항목을 추가하여 `multiqc/` 하위 디렉토리에 게시합니다:

=== "후"

    ```groovy title="rnaseq.nf" linenums="82" hl_lines="7-12"
        align_log {
            path 'align'
        }
        multiqc_report {
            path 'multiqc'
        }
        multiqc_data {
            path 'multiqc'
        }
    }
    ```

=== "전"

    ```groovy title="rnaseq.nf" linenums="82"
        align_log {
            path 'align'
        }
    }
    ```

출력 구성이 완료되었습니다.

### 2.4. 워크플로우 실행

이전 처리 단계가 캐시되고 새로운 MultiQC 단계만 실행되도록 `-resume`을 사용합니다.

```bash
nextflow run rnaseq.nf -profile test -resume
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

캐시된 프로세스 호출 후에 단일 MULTIQC 호출이 추가되었습니다.

결과 디렉토리에서 MultiQC 출력을 찾을 수 있습니다.

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

### 핵심 정리

여러 채널에서 출력을 수집하고, `.mix()`와 `.collect()`로 번들로 묶고, 집계 프로세스에 전달하는 방법을 알게 되었습니다.

### 다음 단계는?

paired-end RNAseq 데이터를 처리하도록 워크플로우를 조정합니다.

---

## 3. paired-end RNAseq 데이터 처리 활성화

현재 우리의 워크플로우는 single-end RNAseq 데이터만 처리할 수 있습니다.
paired-end RNAseq 데이터를 보는 것이 점점 일반적이므로 이를 처리할 수 있기를 원합니다.

워크플로우를 데이터 타입에 완전히 독립적으로 만들려면 약간 더 고급 Nextflow 언어 기능을 사용해야 하므로 여기서는 그렇게 하지 않겠지만, 무엇을 조정해야 하는지 보여주기 위해 paired-end 처리 버전을 만들 수 있습니다.

### 3.1. 워크플로우 복사 및 입력 업데이트

single-end 워크플로우 파일을 복사하고 paired-end 데이터용으로 업데이트하는 것으로 시작합니다.

#### 3.1.1. 워크플로우 파일 복사

paired-end 버전의 시작점으로 사용할 워크플로우 파일의 복사본을 만듭니다.

```bash
cp rnaseq.nf rnaseq_pe.nf
```

이제 새 파일에서 매개변수와 입력 처리를 업데이트합니다.

#### 3.1.2. paired-end 테스트 프로파일 추가

`data/` 디렉토리에 샘플 ID와 paired FASTQ 파일 경로가 포함된 두 번째 CSV 파일을 제공합니다.

```csv title="data/paired-end.csv" linenums="1"
sample_id,fastq_1,fastq_2
ENCSR000COQ1,/workspaces/training/nf4-science/rnaseq/data/reads/ENCSR000COQ1_1.fastq.gz,/workspaces/training/nf4-science/rnaseq/data/reads/ENCSR000COQ1_2.fastq.gz
ENCSR000COQ2,/workspaces/training/nf4-science/rnaseq/data/reads/ENCSR000COQ2_1.fastq.gz,/workspaces/training/nf4-science/rnaseq/data/reads/ENCSR000COQ2_2.fastq.gz
ENCSR000COR1,/workspaces/training/nf4-science/rnaseq/data/reads/ENCSR000COR1_1.fastq.gz,/workspaces/training/nf4-science/rnaseq/data/reads/ENCSR000COR1_2.fastq.gz
ENCSR000COR2,/workspaces/training/nf4-science/rnaseq/data/reads/ENCSR000COR2_1.fastq.gz,/workspaces/training/nf4-science/rnaseq/data/reads/ENCSR000COR2_2.fastq.gz
ENCSR000CPO1,/workspaces/training/nf4-science/rnaseq/data/reads/ENCSR000CPO1_1.fastq.gz,/workspaces/training/nf4-science/rnaseq/data/reads/ENCSR000CPO1_2.fastq.gz
ENCSR000CPO2,/workspaces/training/nf4-science/rnaseq/data/reads/ENCSR000CPO2_1.fastq.gz,/workspaces/training/nf4-science/rnaseq/data/reads/ENCSR000CPO2_2.fastq.gz
```

이 파일을 가리키고 paired-end 보고서 ID를 사용하는 `test_pe` 프로파일을 `nextflow.config`에 추가합니다.

=== "후"

    ```groovy title="nextflow.config" linenums="1" hl_lines="9-13"
    docker.enabled = true

    profiles {
        test {
            params.input = "${projectDir}/data/single-end.csv"
            params.hisat2_index_zip = "${projectDir}/data/genome_index.tar.gz"
            params.report_id = "all_single-end"
        }
        test_pe {
            params.input = "${projectDir}/data/paired-end.csv"
            params.hisat2_index_zip = "${projectDir}/data/genome_index.tar.gz"
            params.report_id = "all_paired-end"
        }
    }
    ```

=== "전"

    ```groovy title="nextflow.config" linenums="1"
    docker.enabled = true

    profiles {
        test {
            params.input = "${projectDir}/data/single-end.csv"
            params.hisat2_index_zip = "${projectDir}/data/genome_index.tar.gz"
            params.report_id = "all_single-end"
        }
    }
    ```

paired-end 데이터용 테스트 프로파일이 준비되었습니다.

#### 3.1.3. 채널 팩토리 업데이트

`.map()` 연산자가 두 FASTQ 파일 경로를 모두 가져와 리스트로 반환하도록 해야 합니다.

=== "후"

    ```groovy title="rnaseq_pe.nf" linenums="25" hl_lines="4"
        // CSV 파일의 내용에서 입력 채널 생성
        read_ch = channel.fromPath(params.input)
            .splitCsv(header: true)
            .map { row -> [file(row.fastq_1), file(row.fastq_2)] }
    ```

=== "전"

    ```groovy title="rnaseq_pe.nf" linenums="25"
        // CSV 파일의 내용에서 입력 채널 생성
        read_ch = channel.fromPath(params.input)
            .splitCsv(header: true)
            .map { row -> file(row.fastq_path) }
    ```

입력 처리가 paired-end 데이터용으로 구성되었습니다.

### 3.2. paired-end 데이터용 FASTQC 모듈 조정

paired-end 버전을 만들기 위해 모듈을 복사합니다:

```bash
cp modules/fastqc.nf modules/fastqc_pe.nf
```

FASTQC 프로세스 입력은 변경할 필요가 없습니다. Nextflow가 두 파일의 리스트를 받으면 둘 다 스테이징하고 `reads`가 두 파일 이름으로 확장됩니다.
필요한 유일한 변경 사항은 출력 블록에 있습니다. 이제 샘플당 두 개의 FastQC 보고서를 받으므로 `simpleName` 기반 패턴에서 와일드카드로 전환합니다.

=== "후"

    ```groovy title="modules/fastqc_pe.nf" linenums="10" hl_lines="2 3"
        output:
        path "*_fastqc.zip", emit: zip
        path "*_fastqc.html", emit: html
    ```

=== "전"

    ```groovy title="modules/fastqc_pe.nf" linenums="10"
        output:
        path "${reads.simpleName}_fastqc.zip", emit: zip
        path "${reads.simpleName}_fastqc.html", emit: html
    ```

이렇게 하면 프로세스가 single-end 또는 paired-end 데이터를 처리할 수 있도록 일반화됩니다.

`rnaseq_pe.nf`의 import를 paired-end 버전을 사용하도록 업데이트합니다:

=== "후"

    ```groovy title="rnaseq_pe.nf" linenums="4" hl_lines="1"
    include { FASTQC } from './modules/fastqc_pe.nf'
    ```

=== "전"

    ```groovy title="rnaseq_pe.nf" linenums="4"
    include { FASTQC } from './modules/fastqc.nf'
    ```

FASTQC 모듈과 import가 paired-end 데이터용으로 업데이트되었습니다.

### 3.3. paired-end 데이터용 TRIM_GALORE 모듈 조정

paired-end 버전을 만들기 위해 모듈을 복사합니다:

```bash
cp modules/trim_galore.nf modules/trim_galore_pe.nf
```

이 모듈은 더 실질적인 변경이 필요합니다:

- 입력이 단일 경로에서 두 경로의 튜플로 변경됩니다
- 명령에 `--paired` 플래그가 추가되고 두 읽기 파일을 모두 받습니다
- 출력이 Trim Galore의 paired-end 명명 규칙을 반영하도록 변경되어 각 읽기 파일에 대해 별도의 FastQC 보고서를 생성합니다

=== "후"

    ```groovy title="modules/trim_galore_pe.nf" linenums="8" hl_lines="2 5 7 8 12"
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

=== "전"

    ```groovy title="modules/trim_galore_pe.nf" linenums="8"
        input:
        path reads

        output:
        path "${reads.simpleName}_trimmed.fq.gz", emit: trimmed_reads
        path "${reads}_trimming_report.txt", emit: trimming_reports
        path "${reads.simpleName}_trimmed_fastqc.{zip,html}", emit: fastqc_reports

        script:
        """
        trim_galore --fastqc ${reads}
        """
    ```

`rnaseq_pe.nf`의 import를 업데이트합니다:

=== "후"

    ```groovy title="rnaseq_pe.nf" linenums="5" hl_lines="1"
    include { TRIM_GALORE } from './modules/trim_galore_pe.nf'
    ```

=== "전"

    ```groovy title="rnaseq_pe.nf" linenums="5"
    include { TRIM_GALORE } from './modules/trim_galore.nf'
    ```

TRIM_GALORE 모듈과 import가 paired-end 데이터용으로 업데이트되었습니다.

### 3.4. paired-end 데이터용 HISAT2_ALIGN 모듈 조정

paired-end 버전을 만들기 위해 모듈을 복사합니다:

```bash
cp modules/hisat2_align.nf modules/hisat2_align_pe.nf
```

이 모듈도 유사한 변경이 필요합니다:

- 입력이 단일 경로에서 두 경로의 튜플로 변경됩니다
- HISAT2 명령이 `-U`(unpaired)에서 `-1`과 `-2`(paired) 읽기 인수로 변경됩니다
- `reads.simpleName`의 모든 사용이 `read1.simpleName`으로 변경됩니다. 이제 쌍의 특정 멤버를 참조하기 때문입니다

=== "후"

    ```groovy title="modules/hisat2_align_pe.nf" linenums="8" hl_lines="2 6 7 12 13 14"
        input:
        tuple path(read1), path(read2)
        path index_zip

        output:
        path "${read1.simpleName}.bam", emit: bam
        path "${read1.simpleName}.hisat2.log", emit: log

        script:
        """
        tar -xzvf ${index_zip}
        hisat2 -x ${index_zip.simpleName} -1 ${read1} -2 ${read2} \
            --new-summary --summary-file ${read1.simpleName}.hisat2.log | \
            samtools view -bS -o ${read1.simpleName}.bam
        """
    ```

=== "전"

    ```groovy title="modules/hisat2_align_pe.nf" linenums="8"
        input:
        path reads
        path index_zip

        output:
        path "${reads.simpleName}.bam", emit: bam
        path "${reads.simpleName}.hisat2.log", emit: log

        script:
        """
        tar -xzvf ${index_zip}
        hisat2 -x ${index_zip.simpleName} -U ${reads} \
            --new-summary --summary-file ${reads.simpleName}.hisat2.log | \
            samtools view -bS -o ${reads.simpleName}.bam
        """
    ```

`rnaseq_pe.nf`의 import를 업데이트합니다:

=== "후"

    ```groovy title="rnaseq_pe.nf" linenums="6" hl_lines="1"
    include { HISAT2_ALIGN } from './modules/hisat2_align_pe.nf'
    ```

=== "전"

    ```groovy title="rnaseq_pe.nf" linenums="6"
    include { HISAT2_ALIGN } from './modules/hisat2_align.nf'
    ```

HISAT2_ALIGN 모듈과 import가 paired-end 데이터용으로 업데이트되었습니다.

### 3.5. paired-end 출력에 대한 MultiQC 집계 업데이트

paired-end `TRIM_GALORE` 프로세스는 이제 하나 대신 두 개의 별도 FastQC 보고서 채널(`fastqc_reports_1`과 `fastqc_reports_2`)을 생성합니다.
`rnaseq_pe.nf`의 `.mix()` 블록을 업데이트하여 둘 다 포함합니다:

=== "후"

    ```groovy title="rnaseq_pe.nf" linenums="40" hl_lines="5 6"
        multiqc_files_ch = channel.empty().mix(
            FASTQC.out.zip,
            FASTQC.out.html,
            TRIM_GALORE.out.trimming_reports,
            TRIM_GALORE.out.fastqc_reports_1,
            TRIM_GALORE.out.fastqc_reports_2,
            HISAT2_ALIGN.out.log,
        )
    ```

=== "전"

    ```groovy title="rnaseq_pe.nf" linenums="40"
        multiqc_files_ch = channel.empty().mix(
            FASTQC.out.zip,
            FASTQC.out.html,
            TRIM_GALORE.out.trimming_reports,
            TRIM_GALORE.out.fastqc_reports,
            HISAT2_ALIGN.out.log,
        )
    ```

MultiQC 집계가 이제 두 세트의 paired-end FastQC 보고서를 모두 포함합니다.

### 3.6. paired-end 출력에 대한 출력 처리 업데이트

`publish:` 섹션과 `output {}` 블록도 paired-end `TRIM_GALORE` 프로세스의 두 개의 별도 FastQC 보고서 채널을 반영해야 합니다.

`rnaseq_pe.nf`의 `publish:` 섹션을 업데이트합니다:

=== "후"

    ```groovy title="rnaseq_pe.nf" linenums="52" hl_lines="6-7"
        publish:
        fastqc_zip = FASTQC.out.zip
        fastqc_html = FASTQC.out.html
        trimmed_reads = TRIM_GALORE.out.trimmed_reads
        trimming_reports = TRIM_GALORE.out.trimming_reports
        trimming_fastqc_1 = TRIM_GALORE.out.fastqc_reports_1
        trimming_fastqc_2 = TRIM_GALORE.out.fastqc_reports_2
        bam = HISAT2_ALIGN.out.bam
        align_log = HISAT2_ALIGN.out.log
        multiqc_report = MULTIQC.out.report
        multiqc_data = MULTIQC.out.data
    }
    ```

=== "전"

    ```groovy title="rnaseq_pe.nf" linenums="52"
        publish:
        fastqc_zip = FASTQC.out.zip
        fastqc_html = FASTQC.out.html
        trimmed_reads = TRIM_GALORE.out.trimmed_reads
        trimming_reports = TRIM_GALORE.out.trimming_reports
        trimming_fastqc = TRIM_GALORE.out.fastqc_reports
        bam = HISAT2_ALIGN.out.bam
        align_log = HISAT2_ALIGN.out.log
        multiqc_report = MULTIQC.out.report
        multiqc_data = MULTIQC.out.data
    }
    ```

`output {}` 블록의 해당 항목을 업데이트합니다:

=== "후"

    ```groovy title="rnaseq_pe.nf" linenums="77" hl_lines="4-9"
        trimming_reports {
            path 'trimming'
        }
        trimming_fastqc_1 {
            path 'trimming'
        }
        trimming_fastqc_2 {
            path 'trimming'
        }
        bam {
            path 'align'
        }
    ```

=== "전"

    ```groovy title="rnaseq_pe.nf" linenums="77"
        trimming_reports {
            path 'trimming'
        }
        trimming_fastqc {
            path 'trimming'
        }
        bam {
            path 'align'
        }
    ```

paired-end 워크플로우가 이제 완전히 업데이트되어 실행할 준비가 되었습니다.

### 3.7. 워크플로우 실행

캐시되지 않으므로 `-resume`을 사용하지 않으며, 이전보다 처리할 데이터가 두 배 더 많지만 1분 이내에 완료되어야 합니다.

```bash
nextflow run rnaseq_pe.nf -profile test_pe
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

이제 워크플로우의 약간 다른 두 버전이 있습니다. 하나는 single-end 읽기 데이터용이고 다른 하나는 paired-end 데이터용입니다.
다음 논리적 단계는 워크플로우가 즉시 두 데이터 타입을 모두 받을 수 있도록 만드는 것이며, 이는 이 과정의 범위를 벗어나지만 후속 과정에서 다룰 수 있습니다.

---

### 핵심 정리

단일 샘플 워크플로우를 여러 샘플의 병렬 처리로 조정하고, 포괄적인 QC 보고서를 생성하며, paired-end 읽기 데이터를 사용하도록 워크플로우를 조정하는 방법을 알게 되었습니다.

### 다음 단계는?

스스로를 크게 칭찬하세요! Nextflow for RNAseq 과정을 완료하셨습니다.

최종 [과정 요약](./next_steps.md)으로 이동하여 학습한 내용을 검토하고 다음 단계를 알아보세요.
