# 파트 2: 단일 샘플 구현

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } AI 지원 번역 - [자세히 알아보기 및 개선 제안](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

이 과정의 이 부분에서는 파트 1에서 실행한 모든 명령을 적용하여 자동으로 실행하는 가장 간단한 워크플로우를 작성하며, 한 번에 하나의 샘플만 처리하는 것을 목표로 합니다.

!!! warning "필수 조건"

    이 레슨을 시작하기 전에 [파트 1: 방법 개요](./01_method.md)를 완료해야 합니다.
    특히, 섹션 1.2.3을 진행하면 이 레슨의 정렬 단계에 필요한 게놈 인덱스 파일(`data/genome_index.tar.gz`)이 생성됩니다.

## 과제

이 과정의 이 부분에서는 다음을 수행하는 워크플로우를 개발합니다:

1. 입력 리드에 대한 품질 관리 실행 (FastQC)
2. 어댑터 트리밍 및 트리밍 후 QC 실행 (Trim Galore)
3. 트리밍된 리드를 참조 게놈에 정렬 (HISAT2)

<figure class="excalidraw">
--8<-- "docs/en/docs/nf4_science/rnaseq/img/rnaseq-wf-02.svg"
</figure>

이는 [파트 1: 방법 개요](./01_method.md#1-single-sample-processing)의 첫 번째 섹션에서 컨테이너 내에서 수동으로 실행한 단계를 자동화합니다.

시작점으로 워크플로우의 주요 부분을 요약한 워크플로우 파일 `rnaseq.nf`와 `modules/` 디렉토리에 각 프로세스의 구조를 요약한 네 개의 모듈 파일(`fastqc.nf`, `trim_galore.nf`, `hisat2_align.nf`, `multiqc.nf`)을 제공합니다.

??? full-code "스캐폴드 파일"

    ```groovy title="rnaseq.nf"
    #!/usr/bin/env nextflow

    // 모듈 INCLUDE 문

    /*
     * Pipeline parameters
     */

    // Primary input

    workflow {

        main:
        // Create input channel

        // Call processes

        publish:
        // Declare outputs to publish
    }

    output {
        // Configure publish targets
    }
    ```

    ```groovy title="modules/fastqc.nf"
    #!/usr/bin/env nextflow

    /*
     * Run FastQC on input reads
     */
    process FASTQC {

        container

        input:

        output:

        script:
        """

        """
    }
    ```

    ```groovy title="modules/trim_galore.nf"
    #!/usr/bin/env nextflow

    /*
     * Trim adapters and run post-trimming QC
     */
    process TRIM_GALORE {

        container

        input:

        output:

        script:
        """

        """
    }
    ```

    ```groovy title="modules/hisat2_align.nf"
    #!/usr/bin/env nextflow

    /*
     * Align reads to a reference genome
     */
    process HISAT2_ALIGN {

        container

        input:

        output:

        script:
        """

        """
    }
    ```

    ```groovy title="modules/multiqc.nf"
    #!/usr/bin/env nextflow

    /*
     * Aggregate QC reports with MultiQC
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

이 파일들은 기능적이지 않으며, 코드의 중요한 부분을 채워 넣을 수 있는 스캐폴드 역할만 합니다.

## 레슨 계획

개발 과정을 더 교육적으로 만들기 위해 세 단계로 나누었습니다:

1. **초기 QC 단계를 실행하는 단일 단계 워크플로우 작성.**
   CLI 매개변수 설정, 입력 채널 생성, 프로세스 모듈 작성, 출력 게시 설정을 다룹니다.
2. **어댑터 트리밍 및 트리밍 후 QC 추가.**
   한 프로세스의 출력을 다른 프로세스의 입력에 연결하여 프로세스를 체인으로 연결하는 방법을 다룹니다.
3. **참조 게놈에 대한 정렬 추가.**
   추가 참조 입력 처리 및 압축 아카이브 작업을 다룹니다.

각 단계는 워크플로우 개발의 특정 측면에 초점을 맞춥니다.

!!! tip "팁"

     올바른 작업 디렉토리에 있는지 확인하세요:
     `cd /workspaces/training/nf4-science/rnaseq`

---

## 1. 초기 QC를 실행하는 단일 단계 워크플로우 작성

이 첫 번째 단계는 기본 사항에 초점을 맞춥니다: FASTQ 파일을 로드하고 품질 관리를 실행합니다.

[파트 1](01_method.md)의 `fastqc` 명령을 떠올려 보세요:

```bash
fastqc <reads>
```

이 명령은 FASTQ 파일을 입력으로 받아 `.zip` 아카이브와 `.html` 요약으로 품질 관리 보고서를 생성합니다.
컨테이너 URI는 `community.wave.seqera.io/library/trim-galore:0.6.10--1bf8ca4e1967cd18`였습니다.

이 정보를 가져와 세 단계로 Nextflow에 적용하겠습니다:

1. 입력 설정
2. QC 프로세스 작성 및 워크플로우에서 호출
3. 출력 처리 설정

### 1.1. 입력 설정

입력 매개변수를 선언하고, 편리한 기본값을 제공하는 테스트 프로필을 생성하고, 입력 채널을 생성해야 합니다.

#### 1.1.1. 입력 매개변수 선언 추가

`rnaseq.nf`의 `Pipeline parameters` 섹션 아래에 `Path` 타입의 `reads` 매개변수를 선언합니다.

=== "후"

    ```groovy title="rnaseq.nf" linenums="5" hl_lines="4-7"
    /*
     * Pipeline parameters
     */
    params {
        // Primary input
        input: Path
    }
    ```

=== "전"

    ```groovy title="rnaseq.nf" linenums="5"
    /*
     * Pipeline parameters
     */

    // Primary input
    ```

이것으로 CLI 매개변수가 설정되었지만, 개발 중에 워크플로우를 실행할 때마다 파일 경로를 입력하고 싶지 않습니다.
기본값을 제공하는 여러 옵션이 있으며, 여기서는 테스트 프로필을 사용합니다.

#### 1.1.2. `nextflow.config`에 기본값이 있는 테스트 프로필 생성

테스트 프로필은 명령줄에서 입력을 지정하지 않고 워크플로우를 시험해 볼 수 있는 편리한 기본값을 제공합니다.
이는 Nextflow 생태계에서 일반적인 관례입니다(자세한 내용은 [Hello Config](../../hello_nextflow/06_hello_config.md) 참조).

`nextflow.config`에 `profiles` 블록을 추가하고 `reads` 매개변수를 테스트 FASTQ 파일 중 하나로 설정하는 `test` 프로필을 추가합니다.

=== "후"

    ```groovy title="nextflow.config" linenums="1" hl_lines="3-7"
    docker.enabled = true

    profiles {
        test {
            params.input = "${projectDir}/data/reads/ENCSR000COQ1_1.fastq.gz"
        }
    }
    ```

=== "전"

    ```groovy title="nextflow.config" linenums="1"
    docker.enabled = true
    ```

여기서는 워크플로우 스크립트가 위치한 디렉토리를 가리키는 Nextflow 내장 변수인 `#!groovy ${projectDir}`를 사용합니다.
이를 통해 절대 경로를 하드코딩하지 않고도 데이터 파일 및 기타 리소스를 쉽게 참조할 수 있습니다.

이제 매개변수에 편리한 기본값이 있습니다. 다음으로 이로부터 채널을 생성해야 합니다.

#### 1.1.3. 입력 채널 설정

workflow 블록에서 `.fromPath` 채널 팩토리를 사용하여 매개변수 값으로부터 입력 채널을 생성합니다([Hello Channels](../../hello_nextflow/02_hello_channels.md)에서 사용됨).

=== "후"

    ```groovy title="rnaseq.nf" linenums="13" hl_lines="4-5"
    workflow {

        main:
        // 파일 경로에서 입력 채널 생성
        read_ch = channel.fromPath(params.input)

        // Call processes

        publish:
        // Declare outputs to publish
    }
    ```

=== "전"

    ```groovy title="rnaseq.nf" linenums="13"
    workflow {

        main:
        // Create input channel

        // Call processes

        publish:
        // Declare outputs to publish
    }
    ```

다음으로 이 입력에 대해 QC를 실행할 프로세스를 생성해야 합니다.

### 1.2. QC 프로세스 작성 및 워크플로우에서 호출

모듈 파일에서 프로세스 정의를 채우고, include 문을 사용하여 워크플로우로 가져오고, 입력에 대해 호출해야 합니다.

#### 1.2.1. QC 프로세스용 모듈 채우기

`modules/fastqc.nf`를 열고 프로세스 정의의 개요를 살펴보세요.
주요 구조 요소를 인식할 수 있어야 합니다. 그렇지 않다면 [Hello Nextflow](../../hello_nextflow/01_hello_world.md)를 읽어보는 것을 고려하세요.

위에 제공된 정보를 사용하여 직접 프로세스 정의를 채운 다음, 아래 "후" 탭의 솔루션과 비교하여 작업을 확인하세요.

=== "전"

    ```groovy title="modules/fastqc.nf" linenums="1"
    #!/usr/bin/env nextflow

    /*
     * Run FastQC on input reads
     */
    process FASTQC {

        container

        input:

        output:

        script:
        """

        """
    }
    ```

=== "후"

    ```groovy title="modules/fastqc.nf" linenums="1" hl_lines="8 11 14 15 19"
    #!/usr/bin/env nextflow

    /*
     * Run FastQC on input reads
     */
    process FASTQC {

        container "community.wave.seqera.io/library/trim-galore:0.6.10--1bf8ca4e1967cd18"

        input:
        path reads

        output:
        path "${reads.simpleName}_fastqc.zip", emit: zip
        path "${reads.simpleName}_fastqc.html", emit: html

        script:
        """
        fastqc ${reads}
        """
    }
    ```

`simpleName` 접근자는 파일명에서 모든 확장자를 제거하므로 `ENCSR000COQ1_1.fastq.gz`는 `ENCSR000COQ1_1`이 됩니다.
`emit:` 구문을 사용하여 각 출력 채널에 이름을 할당하며, 이는 출력을 publish 블록에 연결하는 데 유용합니다.

이것을 완료하면 프로세스가 완성됩니다.
워크플로우에서 사용하려면 모듈을 가져오고 프로세스 호출을 추가해야 합니다.

#### 1.2.2. 모듈 포함

`rnaseq.nf`에 `include` 문을 추가하여 프로세스를 워크플로우에서 사용할 수 있도록 합니다:

=== "후"

    ```groovy title="rnaseq.nf" linenums="3" hl_lines="2"
    // 모듈 INCLUDE 문
    include { FASTQC } from './modules/fastqc.nf'
    ```

=== "전"

    ```groovy title="rnaseq.nf" linenums="3"
    // 모듈 INCLUDE 문
    ```

이제 프로세스를 워크플로우 범위에서 사용할 수 있습니다.

#### 1.2.3. 입력에 대해 QC 프로세스 호출

workflow 블록에 `FASTQC` 호출을 추가하고 입력 채널을 인수로 전달합니다.

=== "후"

    ```groovy title="rnaseq.nf" linenums="14" hl_lines="7-8"
    workflow {

        main:
        // 파일 경로에서 입력 채널 생성
        read_ch = channel.fromPath(params.input)

        // 초기 품질 관리
        FASTQC(read_ch)

        publish:
        // Declare outputs to publish
    }
    ```

=== "전"

    ```groovy title="rnaseq.nf" linenums="14"
    workflow {

        main:
        // 파일 경로에서 입력 채널 생성
        read_ch = channel.fromPath(params.input)

        // Call processes

        publish:
        // Declare outputs to publish
    }
    ```

이제 워크플로우가 입력을 로드하고 QC 프로세스를 실행합니다.
다음으로 출력이 게시되는 방식을 설정해야 합니다.

### 1.3. 출력 처리 설정

게시할 프로세스 출력을 선언하고 출력 위치를 지정해야 합니다.

#### 1.3.1. `publish:` 섹션에서 출력 선언

workflow 블록 내부의 `publish:` 섹션은 게시할 프로세스 출력을 선언합니다.
`FASTQC`의 출력을 명명된 타겟에 할당합니다.

=== "후"

    ```groovy title="rnaseq.nf" linenums="23" hl_lines="2-3"
        publish:
        fastqc_zip = FASTQC.out.zip
        fastqc_html = FASTQC.out.html
    }
    ```

=== "전"

    ```groovy title="rnaseq.nf" linenums="23"
        publish:
        // Declare outputs to publish
    }
    ```

다음으로 Nextflow에 게시된 출력을 어디에 둘지 알려줘야 합니다.

#### 1.3.2. `output {}` 블록에서 출력 타겟 설정

`output {}` 블록은 workflow 외부에 위치하며 각 명명된 타겟이 게시되는 위치를 지정합니다.
두 타겟 모두 `fastqc/` 하위 디렉토리에 게시되도록 설정합니다.

=== "후"

    ```groovy title="rnaseq.nf" linenums="28" hl_lines="2-7"
    output {
        fastqc_zip {
            path 'fastqc'
        }
        fastqc_html {
            path 'fastqc'
        }
    }
    ```

=== "전"

    ```groovy title="rnaseq.nf" linenums="28"
    output {
        // Configure publish targets
    }
    ```

!!! note "참고"

    기본적으로 Nextflow는 출력 파일을 심볼릭 링크로 게시하여 불필요한 중복을 방지합니다.
    여기서 사용하는 데이터 파일은 매우 작지만, 유전체학에서는 매우 커질 수 있습니다.
    심볼릭 링크는 `work` 디렉토리를 정리할 때 깨지므로, 프로덕션 워크플로우에서는 기본 게시 모드를 `'copy'`로 재정의하는 것이 좋습니다.

### 1.4. 워크플로우 실행

이 시점에서 완전히 기능하는 1단계 QC 워크플로우가 있습니다.

테스트 프로필에 설정된 기본값을 사용하기 위해 `-profile test`로 실행하여 명령줄에 경로를 작성할 필요가 없습니다.

```bash
nextflow run rnaseq.nf -profile test
```

??? success "명령 출력"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `rnaseq.nf` [mad_lorenz] DSL2 - revision: 5846a164d2

    executor >  local (1)
    [7b/8ee79e] FASTQC (1) | 1 of 1 ✔
    ```

파트 1을 진행하고 이미 컨테이너를 가져왔다면 매우 빠르게 실행될 것입니다.
건너뛰었다면 Nextflow가 컨테이너를 가져올 것입니다. 이를 위해 아무것도 할 필요가 없지만 최대 1분 정도 기다려야 할 수 있습니다.

results 디렉토리에서 출력을 확인할 수 있습니다.

```bash
ls results/fastqc
```

```console title="출력"
ENCSR000COQ1_1_fastqc.html  ENCSR000COQ1_1_fastqc.zip
```

샘플의 QC 보고서가 이제 `fastqc/` 하위 디렉토리에 게시되었습니다.

### 핵심 정리

프로세스를 포함하는 모듈을 생성하고, 워크플로우로 가져오고, 입력 채널로 호출하고, 워크플로우 수준 output 블록을 사용하여 결과를 게시하는 방법을 배웠습니다.

### 다음 단계는?

워크플로우의 두 번째 단계로 트리밍 후 QC와 함께 어댑터 트리밍을 추가합니다.

---

## 2. 어댑터 트리밍 및 트리밍 후 QC 추가

이제 초기 QC가 준비되었으므로 내장된 트리밍 후 QC와 함께 어댑터 트리밍 단계를 추가할 수 있습니다.

[파트 1](01_method.md)의 `trim_galore` 명령을 떠올려 보세요:

```bash
trim_galore --fastqc <reads>
```

이 명령은 FASTQ 파일에서 어댑터를 트리밍하고 트리밍된 출력에 대해 FastQC를 실행합니다.
트리밍된 리드, 트리밍 보고서, 트리밍된 리드에 대한 FastQC 보고서를 생성합니다.
컨테이너 URI는 `community.wave.seqera.io/library/trim-galore:0.6.10--1bf8ca4e1967cd18`였습니다.

프로세스 정의를 작성하고, 가져오고, 워크플로우에서 호출하고, 출력 처리를 업데이트하기만 하면 됩니다.

### 2.1. 트리밍 프로세스 작성 및 워크플로우에서 호출

이전과 마찬가지로 프로세스 정의를 채우고, 모듈을 가져오고, 프로세스 호출을 추가해야 합니다.

#### 2.1.1. 트리밍 프로세스용 모듈 채우기

`modules/trim_galore.nf`를 열고 프로세스 정의의 개요를 살펴보세요.

위에 제공된 정보를 사용하여 직접 프로세스 정의를 채운 다음, 아래 "후" 탭의 솔루션과 비교하여 작업을 확인하세요.

=== "전"

    ```groovy title="modules/trim_galore.nf" linenums="1"
    #!/usr/bin/env nextflow

    /*
     * Trim adapters and run post-trimming QC
     */
    process TRIM_GALORE {

        container

        input:

        output:

        script:
        """

        """
    }
    ```

=== "후"

    ```groovy title="modules/trim_galore.nf" linenums="1" hl_lines="8 11 14 15 16 20"
    #!/usr/bin/env nextflow

    /*
     * Trim adapters and run post-trimming QC
     */
    process TRIM_GALORE {

        container "community.wave.seqera.io/library/trim-galore:0.6.10--1bf8ca4e1967cd18"

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
    }
    ```

이 프로세스는 세 개의 명명된 출력을 가집니다: 정렬 단계로 전달되는 트리밍된 리드, 트리밍 보고서, 트리밍 후 FastQC 보고서입니다.
`--fastqc` 플래그는 Trim Galore에 트리밍된 출력에 대해 자동으로 FastQC를 실행하도록 지시합니다.

#### 2.1.2. 모듈 포함

`rnaseq.nf`를 업데이트하여 새 모듈을 가져옵니다:

=== "후"

    ```groovy title="rnaseq.nf" linenums="3" hl_lines="3"
    // 모듈 INCLUDE 문
    include { FASTQC } from './modules/fastqc.nf'
    include { TRIM_GALORE } from './modules/trim_galore.nf'
    ```

=== "전"

    ```groovy title="rnaseq.nf" linenums="3"
    // 모듈 INCLUDE 문
    include { FASTQC } from './modules/fastqc.nf'
    ```

다음으로 워크플로우에 프로세스 호출을 추가하겠습니다.

#### 2.1.3. 입력에 대해 트리밍 프로세스 호출

workflow 블록에 프로세스 호출을 추가합니다:

=== "후"

    ```groovy title="rnaseq.nf" linenums="15" hl_lines="10-11"
    workflow {

        main:
        // 파일 경로에서 입력 채널 생성
        read_ch = channel.fromPath(params.input)

        // 초기 품질 관리
        FASTQC(read_ch)

        // 어댑터 트리밍 및 트리밍 후 QC
        TRIM_GALORE(read_ch)

        publish:
        fastqc_zip = FASTQC.out.zip
        fastqc_html = FASTQC.out.html
    }
    ```

=== "전"

    ```groovy title="rnaseq.nf" linenums="15"
    workflow {

        main:
        // 파일 경로에서 입력 채널 생성
        read_ch = channel.fromPath(params.input)

        // 초기 품질 관리
        FASTQC(read_ch)

        publish:
        fastqc_zip = FASTQC.out.zip
        fastqc_html = FASTQC.out.html
    }
    ```

이제 트리밍 프로세스가 워크플로우에 연결되었습니다.

### 2.2. 출력 처리 업데이트

트리밍 출력을 게시 선언에 추가하고 출력 위치를 설정해야 합니다.

#### 2.2.1. 트리밍 출력에 대한 게시 타겟 추가

`publish:` 섹션에 트리밍 출력을 추가합니다:

=== "후"

    ```groovy title="rnaseq.nf" linenums="27" hl_lines="4-6"
        publish:
        fastqc_zip = FASTQC.out.zip
        fastqc_html = FASTQC.out.html
        trimmed_reads = TRIM_GALORE.out.trimmed_reads
        trimming_reports = TRIM_GALORE.out.trimming_reports
        trimming_fastqc = TRIM_GALORE.out.fastqc_reports
    }
    ```

=== "전"

    ```groovy title="rnaseq.nf" linenums="27"
        publish:
        fastqc_zip = FASTQC.out.zip
        fastqc_html = FASTQC.out.html
    }
    ```

다음으로 Nextflow에 이러한 출력을 어디에 둘지 알려줘야 합니다.

#### 2.2.2. 새 출력 타겟 설정

`output {}` 블록에 트리밍 타겟에 대한 항목을 추가하여 `trimming/` 하위 디렉토리에 게시합니다:

=== "후"

    ```groovy title="rnaseq.nf" linenums="35" hl_lines="8-16"
    output {
        fastqc_zip {
            path 'fastqc'
        }
        fastqc_html {
            path 'fastqc'
        }
        trimmed_reads {
            path 'trimming'
        }
        trimming_reports {
            path 'trimming'
        }
        trimming_fastqc {
            path 'trimming'
        }
    }
    ```

=== "전"

    ```groovy title="rnaseq.nf" linenums="35"
    output {
        fastqc_zip {
            path 'fastqc'
        }
        fastqc_html {
            path 'fastqc'
        }
    }
    ```

출력 설정이 완료되었습니다.

### 2.3. 워크플로우 실행

이제 워크플로우에 초기 QC와 어댑터 트리밍이 모두 포함되었습니다.

```bash
nextflow run rnaseq.nf -profile test
```

??? success "명령 출력"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `rnaseq.nf` [gloomy_becquerel] DSL2 - revision: bb11055736

    executor >  local (2)
    [f6/c8ef2e] FASTQC (1)      | 1 of 1 ✔
    [58/c58d8a] TRIM_GALORE (1) | 1 of 1 ✔
    ```

매우 작은 입력 파일에서 실행하고 있으므로 이것도 매우 빠르게 실행될 것입니다.

results 디렉토리에서 트리밍 출력을 찾을 수 있습니다.

```bash
ls results/trimming
```

```console title="출력"
ENCSR000COQ1_1.fastq.gz_trimming_report.txt  ENCSR000COQ1_1_trimmed_fastqc.zip
ENCSR000COQ1_1_trimmed_fastqc.html           ENCSR000COQ1_1_trimmed.fq.gz
```

트리밍 출력과 트리밍 후 QC 보고서가 이제 `trimming/` 하위 디렉토리에 있습니다.

### 핵심 정리

동일한 입력에 대해 독립적으로 실행되어 여러 명명된 출력을 생성하는 두 번째 처리 단계를 추가하는 방법을 배웠습니다.

### 다음 단계는?

트리밍된 리드 출력에서 체인으로 연결되는 정렬 단계를 추가합니다.

---

## 3. 참조 게놈에 대한 정렬 추가

마지막으로 HISAT2를 사용하여 게놈 정렬 단계를 추가할 수 있습니다.

[파트 1](01_method.md)의 정렬 명령을 떠올려 보세요:

```bash
hisat2 -x <genome_index> -U <reads> \
    --new-summary --summary-file <reads>.hisat2.log | \
    samtools view -bS -o <reads>.bam
```

이 명령은 리드를 참조 게놈에 정렬하고 출력을 BAM 형식으로 변환합니다.
사전 구축된 게놈 인덱스 아카이브가 필요하며 BAM 파일과 정렬 요약 로그를 생성합니다.
컨테이너 URI는 `community.wave.seqera.io/library/hisat2_samtools:5e49f68a37dc010e`였습니다.

이 프로세스는 추가 입력(게놈 인덱스 아카이브)이 필요하므로 먼저 이를 설정한 다음 프로세스를 작성하고 연결해야 합니다.

### 3.1. 입력 설정

게놈 인덱스 아카이브에 대한 매개변수를 선언해야 합니다.

#### 3.1.1. 게놈 인덱스에 대한 매개변수 추가

`rnaseq.nf`에 게놈 인덱스 아카이브에 대한 매개변수 선언을 추가합니다:

=== "후"

    ```groovy title="rnaseq.nf" linenums="11" hl_lines="5-6"
    params {
        // Primary input
        input: Path

        // Reference genome archive
        hisat2_index_zip: Path
    }
    ```

=== "전"

    ```groovy title="rnaseq.nf" linenums="11"
    params {
        // Primary input
        input: Path
    }
    ```

#### 3.1.2. 테스트 프로필에 게놈 인덱스 기본값 추가

섹션 1.1.2에서 `reads`에 대해 했던 것처럼, `nextflow.config`의 테스트 프로필에 게놈 인덱스에 대한 기본값을 추가합니다:

=== "후"

    ```groovy title="nextflow.config" linenums="1" hl_lines="6"
    docker.enabled = true

    profiles {
        test {
            params.input = "${projectDir}/data/reads/ENCSR000COQ1_1.fastq.gz"
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
        }
    }
    ```

매개변수가 준비되었습니다. 이제 정렬 프로세스를 생성할 수 있습니다.

### 3.2. 정렬 프로세스 작성 및 워크플로우에서 호출

이전과 마찬가지로 프로세스 정의를 채우고, 모듈을 가져오고, 프로세스 호출을 추가해야 합니다.

#### 3.2.1. 정렬 프로세스용 모듈 채우기

`modules/hisat2_align.nf`를 열고 프로세스 정의의 개요를 살펴보세요.

위에 제공된 정보를 사용하여 직접 프로세스 정의를 채운 다음, 아래 "후" 탭의 솔루션과 비교하여 작업을 확인하세요.

=== "전"

    ```groovy title="modules/hisat2_align.nf" linenums="1"
    #!/usr/bin/env nextflow

    /*
     * Align reads to a reference genome
     */
    process HISAT2_ALIGN {

        container

        input:

        output:

        script:
        """

        """
    }
    ```

=== "후"

    ```groovy title="modules/hisat2_align.nf" linenums="1" hl_lines="8 11 12 15 16 20 21 22 23"
    #!/usr/bin/env nextflow

    /*
     * Align reads to a reference genome
     */
    process HISAT2_ALIGN {

        container "community.wave.seqera.io/library/hisat2_samtools:5e49f68a37dc010e"

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
    }
    ```

이 프로세스는 두 개의 입력을 받습니다: 리드와 게놈 인덱스 아카이브입니다.
script 블록은 먼저 아카이브에서 인덱스를 추출한 다음, HISAT2 정렬을 실행하고 `samtools view`로 파이프하여 출력을 BAM 형식으로 변환합니다.
`index_zip`의 `simpleName` 접근자는 아카이브의 기본 이름(`genome_index`)을 추출하여 인덱스 접두사로 사용합니다.

#### 3.2.2. 모듈 포함

`rnaseq.nf`를 업데이트하여 새 모듈을 가져옵니다:

=== "후"

    ```groovy title="rnaseq.nf" linenums="3" hl_lines="4"
    // 모듈 INCLUDE 문
    include { FASTQC } from './modules/fastqc.nf'
    include { TRIM_GALORE } from './modules/trim_galore.nf'
    include { HISAT2_ALIGN } from './modules/hisat2_align.nf'
    ```

=== "전"

    ```groovy title="rnaseq.nf" linenums="3"
    // 모듈 INCLUDE 문
    include { FASTQC } from './modules/fastqc.nf'
    include { TRIM_GALORE } from './modules/trim_galore.nf'
    ```

다음으로 워크플로우에 프로세스 호출을 추가하겠습니다.

#### 3.2.3. 정렬 프로세스 호출

트리밍된 리드는 이전 단계에서 출력된 `TRIM_GALORE.out.trimmed_reads` 채널에 있습니다.
`#!groovy file(params.hisat2_index_zip)`를 사용하여 게놈 인덱스 아카이브를 제공합니다.

=== "후"

    ```groovy title="rnaseq.nf" linenums="19" hl_lines="14-15"
    workflow {

        main:
        // 파일 경로에서 입력 채널 생성
        read_ch = channel.fromPath(params.input)

        // 초기 품질 관리
        FASTQC(read_ch)

        // 어댑터 트리밍 및 트리밍 후 QC
        TRIM_GALORE(read_ch)

        // 참조 게놈에 정렬
        HISAT2_ALIGN(TRIM_GALORE.out.trimmed_reads, file(params.hisat2_index_zip))
    ```

=== "전"

    ```groovy title="rnaseq.nf" linenums="19"
    workflow {

        main:
        // 파일 경로에서 입력 채널 생성
        read_ch = channel.fromPath(params.input)

        // 초기 품질 관리
        FASTQC(read_ch)

        // 어댑터 트리밍 및 트리밍 후 QC
        TRIM_GALORE(read_ch)
    ```

이제 정렬 프로세스가 워크플로우에 연결되었습니다.

### 3.3. 출력 처리 업데이트

정렬 출력을 게시 선언에 추가하고 출력 위치를 설정해야 합니다.

#### 3.3.1. 정렬 출력에 대한 게시 타겟 추가

`publish:` 섹션에 정렬 출력을 추가합니다:

=== "후"

    ```groovy title="rnaseq.nf" linenums="34" hl_lines="7-8"
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

=== "전"

    ```groovy title="rnaseq.nf" linenums="34"
        publish:
        fastqc_zip = FASTQC.out.zip
        fastqc_html = FASTQC.out.html
        trimmed_reads = TRIM_GALORE.out.trimmed_reads
        trimming_reports = TRIM_GALORE.out.trimming_reports
        trimming_fastqc = TRIM_GALORE.out.fastqc_reports
    }
    ```

다음으로 Nextflow에 이러한 출력을 어디에 둘지 알려줘야 합니다.

#### 3.3.2. 새 출력 타겟 설정

`output {}` 블록에 정렬 타겟에 대한 항목을 추가하여 `align/` 하위 디렉토리에 게시합니다:

=== "후"

    ```groovy title="rnaseq.nf" linenums="44" hl_lines="17-22"
    output {
        fastqc_zip {
            path 'fastqc'
        }
        fastqc_html {
            path 'fastqc'
        }
        trimmed_reads {
            path 'trimming'
        }
        trimming_reports {
            path 'trimming'
        }
        trimming_fastqc {
            path 'trimming'
        }
        bam {
            path 'align'
        }
        align_log {
            path 'align'
        }
    }
    ```

=== "전"

    ```groovy title="rnaseq.nf" linenums="44"
    output {
        fastqc_zip {
            path 'fastqc'
        }
        fastqc_html {
            path 'fastqc'
        }
        trimmed_reads {
            path 'trimming'
        }
        trimming_reports {
            path 'trimming'
        }
        trimming_fastqc {
            path 'trimming'
        }
    }
    ```

출력 설정이 완료되었습니다.

### 3.4. 워크플로우 실행

이제 워크플로우에 세 가지 처리 단계(QC, 트리밍, 정렬)가 모두 포함되었습니다.

```bash
nextflow run rnaseq.nf -profile test
```

??? success "명령 출력"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `rnaseq.nf` [elated_stonebraker] DSL2 - revision: e8e57d0cdd

    executor >  local (3)
    [e8/fa29d6] FASTQC (1)       | 1 of 1 ✔
    [ca/ffdde2] TRIM_GALORE (1)  | 1 of 1 ✔
    [b6/1c6ca3] HISAT2_ALIGN (1) | 1 of 1 ✔
    ```

results 디렉토리에서 정렬 출력을 찾을 수 있습니다.

```bash
ls results/align
```

```console title="출력"
ENCSR000COQ1_1_trimmed.bam  ENCSR000COQ1_1_trimmed.hisat2.log
```

이것으로 각 샘플에 적용해야 하는 기본 처리가 완료됩니다.

_워크플로우가 한 번에 여러 샘플을 허용하도록 만든 후 파트 3에서 MultiQC 보고서 집계를 추가할 것입니다._

---

### 핵심 정리

단일 말단 RNAseq 샘플을 개별적으로 처리하는 모든 핵심 단계를 적용하는 방법을 배웠습니다.

### 다음 단계는?

휴식을 취하세요! 많은 내용이었습니다.

기분이 상쾌해지면 [파트 3](./03_multi-sample.md)으로 이동하세요. 여기서는 여러 샘플을 병렬로 처리하도록 워크플로우를 수정하고, 모든 샘플의 모든 단계에서 QC 보고서를 집계하며, 페어드 엔드 RNAseq 데이터에서 워크플로우를 실행할 수 있도록 하는 방법을 배웁니다.
