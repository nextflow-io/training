# 파트 2: 샘플별 변이 호출

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } AI 지원 번역 - [자세히 알아보기 및 개선 제안](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

파트 1에서는 Samtools와 GATK 명령을 각각의 컨테이너에서 수동으로 테스트했습니다.
이제 동일한 명령을 Nextflow 워크플로우로 적용하겠습니다.

## 과제

이 과정의 이 파트에서는 다음을 수행하는 워크플로우를 개발하겠습니다:

1. [Samtools](https://www.htslib.org/)를 사용하여 각 BAM 입력 파일에 대한 인덱스 파일 생성
2. 각 BAM 입력 파일에 GATK HaplotypeCaller를 실행하여 VCF(Variant Call Format) 형식의 샘플별 변이 호출 생성

<figure class="excalidraw">
--8<-- "docs/en/docs/nf4_science/genomics/img/hello-gatk-1.svg"
</figure>

이는 파트 1에서 컨테이너에서 수동으로 실행했던 단계를 재현합니다.

시작점으로, 워크플로우의 주요 부분을 개요로 제시하는 워크플로우 파일 `genomics.nf`와 모듈의 구조를 개요로 제시하는 두 개의 모듈 파일 samtools_index.nf 및 gatk_haplotypecaller.nf를 제공합니다.
이 파일들은 기능적이지 않으며, 코드의 중요한 부분을 채워 넣을 수 있는 뼈대 역할만 합니다.

## 학습 계획

개발 과정을 더욱 교육적으로 만들기 위해 다음 네 단계로 나누었습니다:

1. **BAM 파일에 Samtools index를 실행하는 단일 단계 워크플로우 작성.**
   모듈 생성, 가져오기 및 워크플로우에서 호출하는 방법을 다룹니다.
2. **인덱싱된 BAM 파일에 GATK HaplotypeCaller를 실행하는 두 번째 프로세스 추가.**
   프로세스 출력을 입력으로 연결하고 보조 파일을 처리하는 방법을 다룹니다.
3. **샘플 배치에서 실행되도록 워크플로우 조정.**
   병렬 실행을 다루고 관련 파일을 함께 유지하기 위한 튜플을 소개합니다.
4. **입력 파일 배치를 포함하는 텍스트 파일을 받도록 워크플로우 수정.**
   대량으로 입력을 제공하는 일반적인 패턴을 보여줍니다.

각 단계는 워크플로우 개발의 특정 측면에 초점을 맞춥니다.

---

## 1. BAM 파일에 Samtools index를 실행하는 단일 단계 워크플로우 작성

첫 번째 단계는 기본 사항에 초점을 맞춥니다: BAM 파일을 로드하고 인덱스를 생성합니다.

[파트 1](01_method.md)의 `samtools index` 명령을 상기해 보세요:

```bash
samtools index '<input_bam>'
```

이 명령은 BAM 파일을 입력으로 받아 그 옆에 `.bai` 인덱스 파일을 생성합니다.
컨테이너 URI는 `community.wave.seqera.io/library/samtools:1.20--b5dfbd93de237464`였습니다.

이 정보를 가져와 세 단계로 Nextflow에 적용하겠습니다:

1. 입력 설정
2. 인덱싱 프로세스 작성 및 워크플로우에서 호출
3. 출력 처리 구성

### 1.1. 입력 설정

입력 매개변수를 선언하고, 편리한 기본값을 제공하는 테스트 프로파일을 생성하고, 입력 채널을 생성해야 합니다.

#### 1.1.1. 입력 매개변수 선언 추가

메인 워크플로우 파일 `genomics.nf`의 `Pipeline parameters` 섹션 아래에 `reads_bam`이라는 CLI 매개변수를 선언하세요.

=== "후"

    ```groovy title="genomics.nf" linenums="5" hl_lines="4-7"
    /*
     * Pipeline parameters
     */
    params {
        // Primary input
        reads_bam: Path
    }
    ```

=== "전"

    ```groovy title="genomics.nf" linenums="5"
    /*
     * Pipeline parameters
     */

    // Primary input
    ```

이것으로 CLI 매개변수가 설정되었지만, 개발 중에 워크플로우를 실행할 때마다 파일 경로를 입력하고 싶지는 않습니다.
기본값을 제공하는 여러 옵션이 있으며, 여기서는 테스트 프로파일을 사용합니다.

#### 1.1.2. `nextflow.config`에 기본값이 있는 테스트 프로파일 생성

테스트 프로파일은 명령줄에서 입력을 지정하지 않고 워크플로우를 시험해 볼 수 있는 편리한 기본값을 제공합니다.
이는 Nextflow 생태계에서 일반적인 관례입니다(자세한 내용은 [Hello Config](../../hello_nextflow/06_hello_config.md) 참조).

`nextflow.config`에 `profiles` 블록을 추가하고 `reads_bam` 매개변수를 테스트 BAM 파일 중 하나로 설정하는 `test` 프로파일을 추가하세요.

=== "후"

    ```groovy title="nextflow.config" linenums="1" hl_lines="3-7"
    docker.enabled = true

    profiles {
        test {
            params.reads_bam = "${projectDir}/data/bam/reads_mother.bam"
        }
    }
    ```

=== "전"

    ```groovy title="nextflow.config" linenums="1"
    docker.enabled = true
    ```

여기서는 워크플로우 스크립트가 위치한 디렉토리를 가리키는 내장 Nextflow 변수인 `${projectDir}`를 사용합니다.
이를 통해 절대 경로를 하드코딩하지 않고도 데이터 파일 및 기타 리소스를 쉽게 참조할 수 있습니다.

#### 1.1.3. 입력 채널 설정

워크플로우 블록에서 `.fromPath` 채널 팩토리를 사용하여 매개변수 값으로부터 입력 채널을 생성하세요([Hello Channels](../../hello_nextflow/02_hello_channels.md)에서 사용된 것처럼).

=== "후"

    ```groovy title="genomics.nf" linenums="13" hl_lines="4-5"
    workflow {

        main:
        // Create input channel (single file via CLI parameter)
        reads_ch = channel.fromPath(params.reads_bam)
    ```

=== "전"

    ```groovy title="genomics.nf" linenums="13"
    workflow {

        main:
        // Create input channel
    ```

이제 이 입력에 대해 인덱싱을 실행할 프로세스를 생성해야 합니다.

### 1.2. 인덱싱 프로세스 작성 및 워크플로우에서 호출

모듈 파일에 프로세스 정의를 작성하고, include 문을 사용하여 워크플로우로 가져오고, 입력에 대해 호출해야 합니다.

#### 1.2.1. 인덱싱 프로세스를 위한 모듈 작성

`modules/samtools_index.nf`를 열고 프로세스 정의의 개요를 살펴보세요.
주요 구조 요소를 인식할 수 있어야 합니다. 그렇지 않다면 [Hello Nextflow](../../hello_nextflow/01_hello_world.md)를 읽어보는 것을 고려하세요.

위에 제공된 정보를 사용하여 직접 프로세스 정의를 작성한 다음, 아래 "후" 탭의 솔루션과 비교하여 확인하세요.

=== "전"

    ```groovy title="modules/samtools_index.nf" linenums="1"
    #!/usr/bin/env nextflow

    /*
     * Generate BAM index file
     */
    process SAMTOOLS_INDEX {

        container

        input:

        output:

        script:
        """

        """
    }
    ```

=== "후"

    ```groovy title="modules/samtools_index.nf" linenums="1" hl_lines="8 11 14 18"
    #!/usr/bin/env nextflow

    /*
     * Generate BAM index file
     */
    process SAMTOOLS_INDEX {

        container 'community.wave.seqera.io/library/samtools:1.20--b5dfbd93de237464'

        input:
        path input_bam

        output:
        path "${input_bam}.bai"

        script:
        """
        samtools index '$input_bam'
        """
    }
    ```

이것으로 프로세스가 완성되었습니다.
워크플로우에서 사용하려면 모듈을 가져오고 프로세스 호출을 추가해야 합니다.

#### 1.2.2. 모듈 포함

`genomics.nf`에 `include` 문을 추가하여 워크플로우에서 프로세스를 사용할 수 있도록 하세요:

=== "후"

    ```groovy title="genomics.nf" linenums="3" hl_lines="2"
    // Module INCLUDE statements
    include { SAMTOOLS_INDEX } from './modules/samtools_index.nf'
    ```

=== "전"

    ```groovy title="genomics.nf" linenums="3"
    // Module INCLUDE statements
    ```

이제 프로세스를 워크플로우 범위에서 사용할 수 있습니다.

#### 1.2.3. 입력에 대해 인덱싱 프로세스 호출

이제 워크플로우 블록에 `SAMTOOLS_INDEX` 호출을 추가하고 입력 채널을 인수로 전달하세요.

=== "후"

    ```groovy title="genomics.nf" linenums="14" hl_lines="7-8"
    workflow {

        main:
        // Create input channel (single file via CLI parameter)
        reads_ch = channel.fromPath(params.reads_bam)

        // Create index file for input BAM file
        SAMTOOLS_INDEX(reads_ch)
    ```

=== "전"

    ```groovy title="genomics.nf" linenums="14"
    workflow {

        main:
        // Create input channel (single file via CLI parameter)
        reads_ch = channel.fromPath(params.reads_bam)

        // Call processes
    ```

이제 워크플로우가 입력을 로드하고 그에 대해 인덱싱 프로세스를 실행합니다.
다음으로 출력이 게시되는 방법을 구성해야 합니다.

### 1.3. 출력 처리 구성

게시할 프로세스 출력을 선언하고 출력 위치를 지정해야 합니다.

#### 1.3.1. `publish:` 섹션에서 출력 선언

워크플로우 블록 내부의 `publish:` 섹션은 게시할 프로세스 출력을 선언합니다.
`SAMTOOLS_INDEX`의 출력을 `bam_index`라는 이름의 대상에 할당하세요.

=== "후"

    ```groovy title="genomics.nf" linenums="22" hl_lines="2"
        publish:
        bam_index = SAMTOOLS_INDEX.out
    }
    ```

=== "전"

    ```groovy title="genomics.nf" linenums="22"
        publish:
        // Declare outputs to publish
    }
    ```

이제 Nextflow에 게시된 출력을 어디에 둘지 알려줘야 합니다.

#### 1.3.2. `output {}` 블록에서 출력 대상 구성

`output {}` 블록은 워크플로우 외부에 위치하며 각 이름이 지정된 대상이 게시되는 위치를 지정합니다.
`bam/` 하위 디렉토리에 게시하는 `bam_index` 대상을 추가하세요.

=== "후"

    ```groovy title="genomics.nf" linenums="26" hl_lines="2-4"
    output {
        bam_index {
            path 'bam'
        }
    }
    ```

=== "전"

    ```groovy title="genomics.nf" linenums="26"
    output {
        // Configure publish targets
    }
    ```

!!! note "참고"

    기본적으로 Nextflow는 출력 파일을 심볼릭 링크로 게시하여 불필요한 중복을 방지합니다.
    여기서 사용하는 데이터 파일은 매우 작지만, 유전체학에서는 매우 커질 수 있습니다.
    심볼릭 링크는 `work` 디렉토리를 정리할 때 깨지므로, 프로덕션 워크플로우의 경우 기본 게시 모드를 `'copy'`로 재정의할 수 있습니다.

### 1.4. 워크플로우 실행

이 시점에서 완전히 기능하는 1단계 인덱싱 워크플로우가 있습니다. 작동하는지 테스트해 봅시다!

테스트 프로파일에 설정된 기본값을 사용하고 명령줄에 경로를 작성하지 않도록 `-profile test`로 실행할 수 있습니다.

```bash
nextflow run genomics.nf -profile test
```

??? success "명령 출력"

    ```console
    N E X T F L O W   ~  version 25.10.2

    ┃ Launching `genomics.nf` [reverent_sinoussi] DSL2 - revision: 41d43ad7fe

    executor >  local (1)
    [2a/e69536] SAMTOOLS_INDEX (1) | 1 of 1 ✔
    ```

작업 디렉토리 또는 결과 디렉토리를 확인하여 인덱스 파일이 올바르게 생성되었는지 확인할 수 있습니다.

??? abstract "작업 디렉토리 내용"

    ```console
    work/2a/e695367b2f60df09cf826b07192dc3
    ├── reads_mother.bam -> /workspaces/training/nf4-science/genomics/data/bam/reads_mother.bam
    └── reads_mother.bam.bai
    ```

??? abstract "결과 디렉토리 내용"

    ```console
    results/
    └── bam/
        └── reads_mother.bam.bai -> ...
    ```

여기 있습니다!

### 핵심 정리

프로세스를 포함하는 모듈을 생성하고, 워크플로우로 가져오고, 입력 채널로 호출하고, 결과를 게시하는 방법을 알게 되었습니다.

### 다음 단계

인덱싱 프로세스의 출력을 가져와 변이 호출을 실행하는 두 번째 단계를 추가합니다.

---

## 2. 인덱싱된 BAM 파일에 GATK HaplotypeCaller를 실행하는 두 번째 프로세스 추가

이제 입력 파일에 대한 인덱스가 있으므로 변이 호출 단계 설정으로 넘어갈 수 있습니다.

[파트 1](01_method.md)의 `gatk HaplotypeCaller` 명령을 상기해 보세요:

```bash
gatk HaplotypeCaller \
        -R /data/ref/ref.fasta \
        -I /data/bam/reads_mother.bam \
        -O reads_mother.vcf \
        -L /data/ref/intervals.bed
```

이 명령은 BAM 파일(`-I`), 참조 게놈(`-R`), 인터벌 파일(`-L`)을 받아 VCF 파일(`-O`)과 그 인덱스를 생성합니다.
이 도구는 또한 BAM 인덱스, 참조 인덱스, 참조 딕셔너리가 각각의 파일과 함께 위치하기를 기대합니다.
컨테이너 URI는 `community.wave.seqera.io/library/gatk4:4.5.0.0--730ee8817e436867`였습니다.

이전과 동일한 세 단계를 따릅니다:

1. 입력 설정
2. 변이 호출 프로세스 작성 및 워크플로우에서 호출
3. 출력 처리 구성

### 2.1. 입력 설정

변이 호출 단계에는 여러 추가 입력 파일이 필요합니다.
이들에 대한 매개변수를 선언하고, 테스트 프로파일에 기본값을 추가하고, 이들을 로드할 변수를 생성해야 합니다.

#### 2.1.1. 보조 입력에 대한 매개변수 선언 추가

새 프로세스가 제공될 추가 파일을 기대하므로, `genomics.nf`의 `Pipeline parameters` 섹션 아래에 이들에 대한 매개변수 선언을 추가하세요:

=== "후"

    ```groovy title="genomics.nf" linenums="9" hl_lines="5-9"
    params {
        // Primary input
        reads_bam: Path

        // Accessory files
        reference: Path
        reference_index: Path
        reference_dict: Path
        intervals: Path
    }
    ```

=== "전"

    ```groovy title="genomics.nf" linenums="9"
    params {
        // Primary input
        reads_bam: Path
    }
    ```

이전과 마찬가지로 인라인이 아닌 테스트 프로파일을 통해 기본값을 제공합니다.

#### 2.1.2. 테스트 프로파일에 보조 파일 기본값 추가

섹션 1.1.2에서 `reads_bam`에 대해 했던 것처럼, `nextflow.config`의 테스트 프로파일에 보조 파일에 대한 기본값을 추가하세요:

=== "후"

    ```groovy title="nextflow.config" linenums="4" hl_lines="3-6"
    test {
        params.reads_bam = "${projectDir}/data/bam/reads_mother.bam"
        params.reference = "${projectDir}/data/ref/ref.fasta"
        params.reference_index = "${projectDir}/data/ref/ref.fasta.fai"
        params.reference_dict = "${projectDir}/data/ref/ref.dict"
        params.intervals = "${projectDir}/data/ref/intervals.bed"
    }
    ```

=== "전"

    ```groovy title="nextflow.config" linenums="4"
    test {
        params.reads_bam = "${projectDir}/data/bam/reads_mother.bam"
    }
    ```

이제 워크플로우에서 사용할 이러한 파일 경로를 로드하는 변수를 생성해야 합니다.

#### 2.1.3. 보조 파일에 대한 변수 생성

워크플로우 블록 내부에 보조 파일 경로에 대한 변수를 추가하세요:

=== "후"

    ```groovy title="genomics.nf" linenums="21" hl_lines="7-11"
    workflow {

        main:
        // Create input channel (single file via CLI parameter)
        reads_ch = channel.fromPath(params.reads_bam)

        // Load the file paths for the accessory files (reference and intervals)
        ref_file        = file(params.reference)
        ref_index_file  = file(params.reference_index)
        ref_dict_file   = file(params.reference_dict)
        intervals_file  = file(params.intervals)

        // Create index file for input BAM file
        SAMTOOLS_INDEX(reads_ch)
    ```

=== "전"

    ```groovy title="genomics.nf" linenums="21"
    workflow {

        main:
        // Create input channel (single file via CLI parameter)
        reads_ch = channel.fromPath(params.reads_bam)

        // Create index file for input BAM file
        SAMTOOLS_INDEX(reads_ch)
    ```

`file()` 구문은 Nextflow에 이러한 입력을 파일 경로로 명시적으로 처리하도록 지시합니다.
이에 대한 자세한 내용은 사이드 퀘스트 [파일 작업](../../side_quests/working_with_files.md)에서 확인할 수 있습니다.

### 2.2. 변이 호출 프로세스 작성 및 워크플로우에서 호출

모듈 파일에 프로세스 정의를 작성하고, include 문을 사용하여 워크플로우로 가져오고, 입력 리드와 인덱싱 단계의 출력 및 보조 파일에 대해 호출해야 합니다.

#### 2.2.1. 변이 호출 프로세스를 위한 모듈 작성

`modules/gatk_haplotypecaller.nf`를 열고 프로세스 정의의 개요를 살펴보세요.

위에 제공된 정보를 사용하여 직접 프로세스 정의를 작성한 다음, 아래 "후" 탭의 솔루션과 비교하여 확인하세요.

=== "전"

    ```groovy title="modules/gatk_haplotypecaller.nf" linenums="1"
    #!/usr/bin/env nextflow

    /*
     * Call variants with GATK HaplotypeCaller
     */
    process GATK_HAPLOTYPECALLER {

        container

        input:

        output:

        script:
        """

        """
    }
    ```

=== "후"

    ```groovy title="modules/gatk_haplotypecaller.nf" linenums="1" hl_lines="8 11-16 19-20 24-28"
    #!/usr/bin/env nextflow

    /*
     * Call variants with GATK HaplotypeCaller
     */
    process GATK_HAPLOTYPECALLER {

        container "community.wave.seqera.io/library/gatk4:4.5.0.0--730ee8817e436867"

        input:
        path input_bam
        path input_bam_index
        path ref_fasta
        path ref_index
        path ref_dict
        path interval_list

        output:
        path "${input_bam}.vcf"     , emit: vcf
        path "${input_bam}.vcf.idx" , emit: idx

        script:
        """
        gatk HaplotypeCaller \
            -R ${ref_fasta} \
            -I ${input_bam} \
            -O ${input_bam}.vcf \
            -L ${interval_list}
        """
    }
    ```

이 프로세스가 GATK 명령 자체가 요구하는 것보다 더 많은 입력을 가지고 있음을 알 수 있습니다.
GATK는 명명 규칙에 따라 BAM 인덱스 파일과 참조 게놈의 보조 파일을 찾을 위치를 알고 있지만, Nextflow는 도메인에 구애받지 않으며 이러한 규칙을 알지 못합니다.
Nextflow가 런타임에 작업 디렉토리에 이들을 스테이징하도록 명시적으로 나열해야 합니다. 그렇지 않으면 GATK가 누락된 파일에 대한 오류를 발생시킵니다.

마찬가지로, 후속 단계를 위해 Nextflow가 추적하도록 출력 VCF의 인덱스 파일(`"${input_bam}.vcf.idx"`)을 명시적으로 나열합니다.
각 출력 채널에 이름을 할당하기 위해 `emit:` 구문을 사용하며, 이는 출력을 게시 블록에 연결할 때 유용해집니다.

이것으로 프로세스가 완성되었습니다.
워크플로우에서 사용하려면 모듈을 가져오고 프로세스 호출을 추가해야 합니다.

#### 2.2.2. 새 모듈 가져오기

`genomics.nf`를 업데이트하여 새 모듈을 가져오세요:

=== "후"

    ```groovy title="genomics.nf" linenums="3" hl_lines="3"
    // Module INCLUDE statements
    include { SAMTOOLS_INDEX } from './modules/samtools_index.nf'
    include { GATK_HAPLOTYPECALLER } from './modules/gatk_haplotypecaller.nf'
    ```

=== "전"

    ```groovy title="genomics.nf" linenums="3"
    // Module INCLUDE statements
    include { SAMTOOLS_INDEX } from './modules/samtools_index.nf'
    ```

이제 프로세스를 워크플로우 범위에서 사용할 수 있습니다.

#### 2.2.3. 프로세스 호출 추가

워크플로우 본문의 `main:` 아래에 프로세스 호출을 추가하세요:

=== "후"

    ```groovy title="genomics.nf" linenums="33" hl_lines="4-12"
        // Create index file for input BAM file
        SAMTOOLS_INDEX(reads_ch)

        // Call variants from the indexed BAM file
        GATK_HAPLOTYPECALLER(
            reads_ch,
            SAMTOOLS_INDEX.out,
            ref_file,
            ref_index_file,
            ref_dict_file,
            intervals_file
        )
    ```

=== "전"

    ```groovy title="genomics.nf" linenums="33"
        // Create index file for input BAM file
        SAMTOOLS_INDEX(reads_ch)
    ```

Hello Nextflow 교육 시리즈에서 `*.out` 구문을 인식할 수 있어야 합니다. Nextflow에 `SAMTOOLS_INDEX`가 출력한 채널을 가져와 `GATK_HAPLOTYPECALLER` 프로세스 호출에 연결하도록 지시하고 있습니다.

!!! note "참고"

    입력이 프로세스 호출에서 프로세스의 입력 블록에 나열된 것과 정확히 동일한 순서로 제공된다는 점에 주목하세요.
    Nextflow에서 입력은 위치 기반이므로 동일한 순서를 _반드시_ 따라야 하며, 물론 동일한 수의 요소가 있어야 합니다.

### 2.3. 출력 처리 구성

게시 선언에 새 출력을 추가하고 출력 위치를 구성해야 합니다.

#### 2.3.1. 변이 호출 출력에 대한 게시 대상 추가

`publish:` 섹션에 VCF 및 인덱스 출력을 추가하세요:

=== "후"

    ```groovy title="genomics.nf" linenums="45" hl_lines="3-4"
        publish:
        bam_index = SAMTOOLS_INDEX.out
        vcf = GATK_HAPLOTYPECALLER.out.vcf
        vcf_idx = GATK_HAPLOTYPECALLER.out.idx
    }
    ```

=== "전"

    ```groovy title="genomics.nf" linenums="45"
        publish:
        bam_index = SAMTOOLS_INDEX.out
    }
    ```

이제 Nextflow에 새 출력을 어디에 둘지 알려줘야 합니다.

#### 2.3.2. 새 출력 대상 구성

`output {}` 블록에 `vcf` 및 `vcf_idx` 대상에 대한 항목을 추가하고, 둘 다 `vcf/` 하위 디렉토리에 게시하세요:

=== "후"

    ```groovy title="genomics.nf" linenums="51" hl_lines="5-10"
    output {
        bam_index {
            path 'bam'
        }
        vcf {
            path 'vcf'
        }
        vcf_idx {
            path 'vcf'
        }
    }
    ```

=== "전"

    ```groovy title="genomics.nf" linenums="49"
    output {
        bam_index {
            path 'bam'
        }
    }
    ```

VCF와 그 인덱스는 둘 다 `vcf/` 하위 디렉토리로 이동하는 별도의 대상으로 게시됩니다.

### 2.4. 워크플로우 실행

확장된 워크플로우를 실행하되, 이번에는 `-resume`을 추가하여 인덱싱 단계를 다시 실행하지 않도록 합니다.

```bash
nextflow run genomics.nf -profile test -resume
```

??? success "명령 출력"

    ```console
    N E X T F L O W   ~  version 25.10.2

    ┃ Launching `genomics.nf` [grave_volta] DSL2 - revision: 4790abc96a

    executor >  local (1)
    [2a/e69536] SAMTOOLS_INDEX (1)       | 1 of 1, cached: 1 ✔
    [53/e18e98] GATK_HAPLOTYPECALLER (1) | 1 of 1 ✔
    ```

이제 콘솔 출력을 보면 두 프로세스가 나열되어 있습니다.

예상대로 캐싱 덕분에 첫 번째 프로세스는 건너뛰었고, 두 번째 프로세스는 새로운 것이므로 실행되었습니다.

결과 디렉토리에서 출력 파일을 찾을 수 있습니다(작업 디렉토리에 대한 심볼릭 링크로).

??? abstract "디렉토리 내용"

    ```console
    results/
    ├── bam/
    │   └── reads_mother.bam.bai -> ...
    └── vcf/
        ├── reads_mother.bam.vcf -> ...
        └── reads_mother.bam.vcf.idx -> ...
    ```

VCF 파일을 열면 컨테이너에서 GATK 명령을 직접 실행하여 생성한 파일과 동일한 내용을 볼 수 있습니다.

??? abstract "파일 내용"

    ```console title="reads_mother.bam.vcf" linenums="26"
    #CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	reads_mother
    20_10037292_10066351	3480	.	C	CT	503.03	.	AC=2;AF=1.00;AN=2;DP=23;ExcessHet=0.0000;FS=0.000;MLEAC=2;MLEAF=1.00;MQ=60.00;QD=27.95;SOR=1.179	GT:AD:DP:GQ:PL	1/1:0,18:18:54:517,54,0
    20_10037292_10066351	3520	.	AT	A	609.03	.	AC=2;AF=1.00;AN=2;DP=18;ExcessHet=0.0000;FS=0.000;MLEAC=2;MLEAF=1.00;MQ=60.00;QD=33.83;SOR=0.693	GT:AD:DP:GQ:PL	1/1:0,18:18:54:623,54,0
    20_10037292_10066351	3529	.	T	A	155.64	.	AC=1;AF=0.500;AN=2;BaseQRankSum=-0.544;DP=21;ExcessHet=0.0000;FS=1.871;MLEAC=1;MLEAF=0.500;MQ=60.00;MQRankSum=0.000;QD=7.78;ReadPosRankSum=-1.158;SOR=1.034	GT:AD:DP:GQ:PL	0/1:12,8:20:99:163,0,328
    ```

이것이 연구의 각 샘플에 대해 생성하고자 하는 출력입니다.

### 핵심 정리

실제 분석 작업을 수행하고 보조 파일과 같은 유전체학 파일 형식의 특이성을 처리할 수 있는 2단계 모듈식 워크플로우를 만드는 방법을 알게 되었습니다.

### 다음 단계

워크플로우가 여러 샘플을 대량으로 처리하도록 만듭니다.

---

## 3. 샘플 배치에서 실행되도록 워크플로우 조정

단일 샘플에 대한 처리를 자동화할 수 있는 워크플로우를 갖는 것은 좋지만, 1000개의 샘플이 있다면 어떻게 될까요?
모든 샘플을 반복하는 bash 스크립트를 작성해야 할까요?

아니요, 다행히도 그렇지 않습니다! 코드를 약간만 수정하면 Nextflow가 그것도 처리해 줍니다.

### 3.1. 세 개의 샘플을 나열하도록 입력 업데이트

여러 샘플에서 실행하려면 테스트 프로파일을 업데이트하여 단일 파일 경로 대신 파일 경로 배열을 제공하세요.
이는 다중 샘플 실행을 테스트하는 빠른 방법입니다. 다음 단계에서는 입력 파일을 사용하는 보다 확장 가능한 접근 방식으로 전환하겠습니다.

먼저, 배열은 타입 선언을 사용할 수 없으므로 매개변수 선언에서 타입 주석을 주석 처리하세요:

=== "후"

    ```groovy title="genomics.nf" linenums="10" hl_lines="1-2"
        // Primary input (array of three samples)
        reads_bam //: Path
    ```

=== "전"

    ```groovy title="genomics.nf" linenums="10"
        // Primary input
        reads_bam: Path
    ```

그런 다음 테스트 프로파일을 업데이트하여 세 샘플을 모두 나열하세요:

=== "후"

    ```groovy title="nextflow.config" linenums="4" hl_lines="2-6"
    test {
        params.reads_bam = [
            "${projectDir}/data/bam/reads_mother.bam",
            "${projectDir}/data/bam/reads_father.bam",
            "${projectDir}/data/bam/reads_son.bam"
        ]
        params.reference = "${projectDir}/data/ref/ref.fasta"
        params.reference_index = "${projectDir}/data/ref/ref.fasta.fai"
        params.reference_dict = "${projectDir}/data/ref/ref.dict"
        params.intervals = "${projectDir}/data/ref/intervals.bed"
    }
    ```

=== "전"

    ```groovy title="nextflow.config" linenums="4"
    test {
        params.reads_bam = "${projectDir}/data/bam/reads_mother.bam"
        params.reference = "${projectDir}/data/ref/ref.fasta"
        params.reference_index = "${projectDir}/data/ref/ref.fasta.fai"
        params.reference_dict = "${projectDir}/data/ref/ref.dict"
        params.intervals = "${projectDir}/data/ref/intervals.bed"
    }
    ```

워크플로우 본문의 채널 팩토리(`.fromPath`)는 단일 파일 경로뿐만 아니라 여러 파일 경로도 받아들이므로 다른 변경은 필요하지 않습니다.

### 3.2. 워크플로우 실행

이제 세 개의 테스트 샘플 모두에서 실행되도록 배관이 설정되었으므로 워크플로우를 실행해 보세요.

```bash
nextflow run genomics.nf -profile test -resume
```

재미있는 점: 이것은 _작동할 수도_ 있고, _실패할 수도_ 있습니다. 예를 들어, 다음은 성공한 실행입니다:

??? success "명령 출력"

    ```console
    N E X T F L O W   ~  version 25.10.2

    ┃ Launching `genomics.nf` [peaceful_yalow] DSL2 - revision: a256d113ad

    executor >  local (6)
    [4f/7071b0] SAMTOOLS_INDEX (3)       | 3 of 3, cached: 1 ✔
    [7a/89bc43] GATK_HAPLOTYPECALLER (2) | 3 of 3, cached: 1 ✔
    ```

워크플로우 실행이 성공했다면 다음과 같은 오류가 발생할 때까지 다시 실행하세요:

??? failure "명령 출력"

    ```console
    N E X T F L O W   ~  version 25.10.2

    ┃ Launching `genomics.nf` [loving_pasteur] DSL2 - revision: d2a8e63076

    executor >  local (4)
    [01/eea165] SAMTOOLS_INDEX (2)       | 3 of 3, cached: 1 ✔
    [a5/fa9fd0] GATK_HAPLOTYPECALLER (3) | 1 of 3, cached: 1
    ERROR ~ Error executing process > 'GATK_HAPLOTYPECALLER (2)'

    Caused by:
      Process `GATK_HAPLOTYPECALLER (2)` terminated with an error exit status (2)

    Command executed:

      gatk HaplotypeCaller         -R ref.fasta         -I reads_father.bam         -O reads_father.bam.vcf         -L intervals.bed

    Command exit status:
      2

    Command error:
      ...
      A USER ERROR has occurred: Traversal by intervals was requested but some input files are not indexed.
      ...
    ```

GATK 명령 오류 출력을 보면 다음과 같은 줄이 있을 것입니다:

```console
A USER ERROR has occurred: Traversal by intervals was requested but some input files are not indexed.
```

워크플로우의 첫 번째 단계에서 BAM 파일을 명시적으로 인덱싱했는데 이상합니다. 배관에 문제가 있을 수 있을까요?

### 3.3. 문제 해결

작업 디렉토리를 검사하고 `view()` 연산자를 사용하여 무엇이 잘못되었는지 파악하겠습니다.

#### 3.3.1. 관련 호출에 대한 작업 디렉토리 확인

콘솔 출력에 나열된 실패한 `GATK_HAPLOTYPECALLER` 프로세스 호출의 작업 디렉토리 내부를 살펴보세요.

??? abstract "디렉토리 내용"

    ```console
    work/a5/fa9fd0994b6beede5fb9ea073596c2
    ├── intervals.bed -> /workspaces/training/nf4-science/genomics/data/ref/intervals.bed
    ├── reads_father.bam.bai -> /workspaces/training/nf4-science/genomics/work/01/eea16597bd6e810fb4cf89e60f8c2d/reads_father.bam.bai
    ├── reads_son.bam -> /workspaces/training/nf4-science/genomics/data/bam/reads_son.bam
    ├── reads_son.bam.vcf
    ├── reads_son.bam.vcf.idx
    ├── ref.dict -> /workspaces/training/nf4-science/genomics/data/ref/ref.dict
    ├── ref.fasta -> /workspaces/training/nf4-science/genomics/data/ref/ref.fasta
    └── ref.fasta.fai -> /workspaces/training/nf4-science/genomics/data/ref/ref.fasta.fai
    ```

이 디렉토리에 나열된 BAM 파일과 BAM 인덱스의 이름에 특히 주의하세요: `reads_son.bam`과 `reads_father.bam.bai`.

뭐라고요? Nextflow가 이 프로세스 호출의 작업 디렉토리에 인덱스 파일을 스테이징했지만 잘못된 것입니다. 어떻게 이런 일이 발생할 수 있을까요?

#### 3.3.2. [view() 연산자](https://www.nextflow.io/docs/latest/reference/operator.html#view)를 사용하여 채널 내용 검사

워크플로우 본문의 `GATK_HAPLOTYPECALLER` 프로세스 호출 전에 다음 두 줄을 추가하여 채널의 내용을 확인하세요:

=== "후"

    ```groovy title="genomics.nf" hl_lines="3-5"
        SAMTOOLS_INDEX(reads_ch)

        // temporary diagnostics
        reads_ch.view()
        SAMTOOLS_INDEX.out.view()

        // Call variants from the indexed BAM file
        GATK_HAPLOTYPECALLER(
    ```

=== "전"

    ```groovy title="genomics.nf"
        SAMTOOLS_INDEX(reads_ch)

        // Call variants from the indexed BAM file
        GATK_HAPLOTYPECALLER(
    ```

그런 다음 워크플로우 명령을 다시 실행하세요.

```bash
nextflow run genomics.nf -profile test
```

다시 한 번, 이것은 성공하거나 실패할 수 있습니다. 다음은 실패한 실행에 대한 두 `.view()` 호출의 출력입니다:

```console
/workspaces/training/nf4-science/genomics/data/bam/reads_mother.bam
/workspaces/training/nf4-science/genomics/data/bam/reads_father.bam
/workspaces/training/nf4-science/genomics/data/bam/reads_son.bam
/workspaces/training/nf4-science/genomics/work/9c/53492e3518447b75363e1cd951be4b/reads_father.bam.bai
/workspaces/training/nf4-science/genomics/work/cc/37894fffdf6cc84c3b0b47f9b536b7/reads_son.bam.bai
/workspaces/training/nf4-science/genomics/work/4d/dff681a3d137ba7d9866e3d9307bd0/reads_mother.bam.bai
```

처음 세 줄은 입력 채널에 해당하고 두 번째는 출력 채널에 해당합니다.
세 샘플의 BAM 파일과 인덱스 파일이 동일한 순서로 나열되지 않았음을 알 수 있습니다!

!!! note "참고"

    여러 요소를 포함하는 채널에서 Nextflow 프로세스를 호출하면 Nextflow는 가능한 한 많이 실행을 병렬화하려고 시도하며, 사용 가능해지는 순서대로 출력을 수집합니다.
    결과적으로 해당 출력은 원래 입력이 제공된 순서와 다른 순서로 수집될 수 있습니다.

현재 작성된 대로 워크플로우 스크립트는 인덱스 파일이 입력이 제공된 것과 동일한 mother/father/son 순서로 인덱싱 단계에서 나올 것으로 가정합니다.
그러나 이것이 항상 그런 것은 아니므로 때때로(항상은 아니지만) 두 번째 단계에서 잘못된 파일이 쌍을 이룹니다.

이를 수정하려면 BAM 파일과 그 인덱스 파일이 채널을 통해 함께 이동하도록 해야 합니다.

!!! tip "팁"

    워크플로우 코드의 `view()` 문은 아무것도 하지 않으므로 그대로 두어도 문제가 되지 않습니다.
    그러나 콘솔 출력을 어지럽히므로 문제 해결이 끝나면 제거하는 것이 좋습니다.

### 3.4. 인덱스 파일을 올바르게 처리하도록 워크플로우 업데이트

수정 방법은 각 BAM 파일을 그 인덱스와 함께 튜플로 묶은 다음, 다운스트림 프로세스와 워크플로우 배관을 일치하도록 업데이트하는 것입니다.

#### 3.4.1. SAMTOOLS_INDEX 모듈의 출력을 튜플로 변경

BAM 파일과 그 인덱스가 밀접하게 연관되도록 보장하는 가장 간단한 방법은 인덱스 작업에서 나오는 튜플로 함께 패키징하는 것입니다.

!!! note "참고"

    **튜플**은 함수에서 여러 값을 반환하는 데 일반적으로 사용되는 유한하고 순서가 지정된 요소 목록입니다. 튜플은 연관성과 순서를 유지하면서 프로세스 간에 여러 입력 또는 출력을 전달하는 데 특히 유용합니다.

`modules/samtools_index.nf`의 출력을 업데이트하여 BAM 파일을 포함하세요:

=== "후"

    ```groovy title="modules/samtools_index.nf" linenums="14" hl_lines="2"
        output:
        tuple path(input_bam), path("${input_bam}.bai")
    ```

=== "전"

    ```groovy title="modules/samtools_index.nf" linenums="14"
        output:
        path "${input_bam}.bai"
    ```

이렇게 하면 각 인덱스 파일이 원래 BAM 파일과 긴밀하게 결합되고, 인덱싱 단계의 전체 출력은 파일 쌍을 포함하는 단일 채널이 됩니다.

#### 3.4.2. GATK_HAPLOTYPECALLER 모듈의 입력을 튜플을 받도록 변경

첫 번째 프로세스의 출력 '형태'를 변경했으므로 두 번째 프로세스의 입력 정의를 일치하도록 업데이트해야 합니다.

`modules/gatk_haplotypecaller.nf`를 업데이트하세요:

=== "후"

    ```groovy title="modules/gatk_haplotypecaller.nf" linenums="11" hl_lines="2"
        input:
        tuple path(input_bam), path(input_bam_index)
    ```

=== "전"

    ```groovy title="modules/gatk_haplotypecaller.nf" linenums="11"
        input:
        path input_bam
        path input_bam_index
    ```

이제 프로세스 호출과 게시 대상에서 새 튜플 구조를 반영하도록 워크플로우를 업데이트해야 합니다.

#### 3.4.3. 워크플로우에서 GATK_HAPLOTYPECALLER 호출 업데이트

BAM 파일이 이제 `SAMTOOLS_INDEX`의 채널 출력에 번들로 포함되어 있으므로 더 이상 `GATK_HAPLOTYPECALLER` 프로세스에 원래 `reads_ch`를 제공할 필요가 없습니다.

`genomics.nf`의 호출을 업데이트하세요:

=== "후"

    ```groovy title="genomics.nf" linenums="42" hl_lines="2"
        GATK_HAPLOTYPECALLER(
            SAMTOOLS_INDEX.out,
    ```

=== "전"

    ```groovy title="genomics.nf" linenums="42"
        GATK_HAPLOTYPECALLER(
            reads_ch,
            SAMTOOLS_INDEX.out,
    ```

마지막으로 새 출력 구조를 반영하도록 게시 대상을 업데이트해야 합니다.

#### 3.4.4. 인덱싱된 BAM 출력에 대한 게시 대상 업데이트

SAMTOOLS_INDEX 출력이 이제 BAM 파일과 그 인덱스를 모두 포함하는 튜플이므로, 게시 대상의 이름을 `bam_index`에서 `indexed_bam`으로 변경하여 내용을 더 잘 반영하세요:

=== "후"

    ```groovy title="genomics.nf" linenums="46" hl_lines="2 8"
        publish:
        indexed_bam = SAMTOOLS_INDEX.out
        vcf = GATK_HAPLOTYPECALLER.out.vcf
        vcf_idx = GATK_HAPLOTYPECALLER.out.idx
    }

    output {
        indexed_bam {
            path 'bam'
        }
        vcf {
            path 'vcf'
        }
        vcf_idx {
            path 'vcf'
        }
    }
    ```

=== "전"

    ```groovy title="genomics.nf" linenums="46"
        publish:
        bam_index = SAMTOOLS_INDEX.out
        vcf = GATK_HAPLOTYPECALLER.out.vcf
        vcf_idx = GATK_HAPLOTYPECALLER.out.idx
    }

    output {
        bam_index {
            path 'bam'
        }
        vcf {
            path 'vcf'
        }
        vcf_idx {
            path 'vcf'
        }
    }
    ```

이러한 변경으로 BAM과 그 인덱스가 함께 이동하도록 보장되므로 쌍이 항상 올바릅니다.

### 3.5. 수정된 워크플로우 실행

앞으로 안정적으로 작동하는지 확인하기 위해 워크플로우를 다시 실행하세요.

```bash
nextflow run genomics.nf -profile test
```

이번에는(그리고 매번) 모든 것이 올바르게 실행되어야 합니다:

??? success "명령 출력"

    ```console
    N E X T F L O W   ~  version 25.10.2

    ┃ Launching `genomics.nf` [special_goldstine] DSL2 - revision: 4cbbf6ea3e

    executor >  local (6)
    [d6/10c2c4] SAMTOOLS_INDEX (1)       | 3 of 3 ✔
    [88/1783aa] GATK_HAPLOTYPECALLER (2) | 3 of 3 ✔
    ```

결과 디렉토리에는 이제 각 샘플에 대한 BAM 및 BAI 파일(튜플에서)과 VCF 출력이 모두 포함됩니다:

??? abstract "결과 디렉토리 내용"

    ```console
    results/
    ├── bam/
    │   ├── reads_father.bam -> ...
    │   ├── reads_father.bam.bai -> ...
    │   ├── reads_mother.bam -> ...
    │   ├── reads_mother.bam.bai -> ...
    │   ├── reads_son.bam -> ...
    │   └── reads_son.bam.bai -> ...
    └── vcf/
        ├── reads_father.bam.vcf -> ...
        ├── reads_father.bam.vcf.idx -> ...
        ├── reads_mother.bam.vcf -> ...
        ├── reads_mother.bam.vcf.idx -> ...
        ├── reads_son.bam.vcf -> ...
        └── reads_son.bam.vcf.idx -> ...
    ```

관련 파일을 튜플로 묶음으로써 올바른 파일이 항상 워크플로우를 통해 함께 이동하도록 보장했습니다.
이제 워크플로우는 임의의 수의 샘플을 안정적으로 처리하지만, 구성에서 개별적으로 나열하는 것은 확장성이 좋지 않습니다.
다음 단계에서는 파일에서 입력을 읽는 것으로 전환하겠습니다.

### 핵심 정리

워크플로우가 여러 샘플에서 (독립적으로) 실행되도록 만드는 방법을 알게 되었습니다.

### 다음 단계

샘플을 대량으로 처리하기 쉽게 만듭니다.

---

## 4. 입력 파일 배치를 포함하는 텍스트 파일을 받도록 워크플로우 수정

워크플로우에 여러 데이터 입력 파일을 제공하는 매우 일반적인 방법은 파일 경로를 포함하는 텍스트 파일로 제공하는 것입니다.
한 줄에 하나의 파일 경로만 나열하는 간단한 텍스트 파일일 수도 있고, 파일에 추가 메타데이터가 포함될 수도 있으며, 이 경우 종종 샘플시트라고 합니다.

여기서는 간단한 경우를 수행하는 방법을 보여드리겠습니다.

### 4.1. 입력 파일 경로를 나열하는 제공된 텍스트 파일 검사

`data/` 디렉토리에서 찾을 수 있는 `sample_bams.txt`라는 입력 파일 경로를 나열하는 텍스트 파일을 이미 만들었습니다.

```txt title="sample_bams.txt"
/workspaces/training/nf4-science/genomics/data/bam/reads_mother.bam
/workspaces/training/nf4-science/genomics/data/bam/reads_father.bam
/workspaces/training/nf4-science/genomics/data/bam/reads_son.bam
```

보시다시피 한 줄에 하나의 파일 경로를 나열했으며 절대 경로입니다.

!!! note "참고"

    여기서 사용하는 파일은 GitHub Codespaces의 로컬 파일 시스템에 있지만 클라우드 스토리지의 파일을 가리킬 수도 있습니다.
    제공된 Codespaces 환경을 사용하지 않는 경우 로컬 설정에 맞게 파일 경로를 조정해야 할 수 있습니다.

### 4.2. 매개변수 및 테스트 프로파일 업데이트

개별 샘플을 나열하는 대신 `sample_bams.txt` 파일을 가리키도록 `reads_bam` 매개변수를 전환하세요.

params 블록에서 타입 주석을 복원하세요(다시 단일 경로이므로):

=== "후"

    ```groovy title="genomics.nf" linenums="10" hl_lines="1-2"
        // Primary input (file of input files, one per line)
        reads_bam: Path
    ```

=== "전"

    ```groovy title="genomics.nf" linenums="10"
        // Primary input (array of three samples)
        reads_bam
    ```

그런 다음 테스트 프로파일을 업데이트하여 텍스트 파일을 가리키도록 하세요:

=== "후"

    ```groovy title="nextflow.config" linenums="4" hl_lines="2"
    test {
        params.reads_bam = "${projectDir}/data/sample_bams.txt"
        params.reference = "${projectDir}/data/ref/ref.fasta"
        params.reference_index = "${projectDir}/data/ref/ref.fasta.fai"
        params.reference_dict = "${projectDir}/data/ref/ref.dict"
        params.intervals = "${projectDir}/data/ref/intervals.bed"
    }
    ```

=== "전"

    ```groovy title="nextflow.config" linenums="4"
    test {
        params.reads_bam = [
            "${projectDir}/data/bam/reads_mother.bam",
            "${projectDir}/data/bam/reads_father.bam",
            "${projectDir}/data/bam/reads_son.bam"
        ]
        params.reference = "${projectDir}/data/ref/ref.fasta"
        params.reference_index = "${projectDir}/data/ref/ref.fasta.fai"
        params.reference_dict = "${projectDir}/data/ref/ref.dict"
        params.intervals = "${projectDir}/data/ref/intervals.bed"
    }
    ```

파일 목록이 더 이상 코드에 전혀 존재하지 않으며, 이는 올바른 방향으로 나아가는 큰 단계입니다.

### 4.3. 파일에서 줄을 읽도록 채널 팩토리 업데이트

현재 입력 채널 팩토리는 제공하는 모든 파일을 인덱싱 프로세스에 공급하려는 데이터 입력으로 처리합니다.
이제 입력 파일 경로를 나열하는 파일을 제공하고 있으므로 파일을 분석하고 포함된 파일 경로를 데이터 입력으로 처리하도록 동작을 변경해야 합니다.

[Hello Nextflow의 파트 2](../../hello_nextflow/02_hello_channels.md#42-use-the-splitcsv-operator-to-parse-the-file)에서 사용한 것과 동일한 패턴을 사용하여 이를 수행할 수 있습니다: [`splitCsv()`](https://nextflow.io/docs/latest/reference/operator.html#splitcsv) 연산자를 적용하여 파일을 분석한 다음 `map` 연산을 사용하여 각 줄의 첫 번째 필드를 선택합니다.

=== "후"

    ```groovy title="genomics.nf" linenums="24" hl_lines="1-4"
        // Create input channel from a CSV file listing input file paths
        reads_ch = Channel.fromPath(params.reads_bam)
                .splitCsv()
                .map { line -> file(line[0]) }
    ```

=== "전"

    ```groovy title="genomics.nf" linenums="24"
        // Create input channel (single file via CLI parameter)
        reads_ch = channel.fromPath(params.reads_bam)
    ```

기술적으로 입력 파일에 현재 파일 경로만 포함되어 있으므로 [`.splitText()`](https://www.nextflow.io/docs/latest/reference/operator.html#operator-splittext) 연산자를 사용하여 더 간단하게 수행할 수 있습니다.
그러나 더 다재다능한 `splitCsv` 연산자를 사용함으로써(`map`으로 보완), 파일 경로를 포함하는 파일에 메타데이터를 추가하기로 결정하는 경우 워크플로우를 미래에 대비할 수 있습니다.

!!! tip "팁"

    연산자가 여기서 무엇을 하는지 확신이 서지 않는다면, `.view()` 연산자를 사용하여 적용하기 전후의 채널 내용을 확인하는 것이 또 다른 좋은 기회입니다.

### 4.4. 워크플로우 실행

워크플로우를 한 번 더 실행하세요. 이전과 동일한 결과를 생성해야 하지 않을까요?

```bash
nextflow run genomics.nf -profile test -resume
```

??? success "명령 출력"

    ```console
    N E X T F L O W   ~  version 25.10.2

    ┃ Launching `genomics.nf` [sick_albattani] DSL2 - revision: 46d84642f6

    [18/23b4bb] SAMTOOLS_INDEX (1)       | 3 of 3, cached: 3 ✔
    [12/f727bb] GATK_HAPLOTYPECALLER (3) | 3 of 3, cached: 3 ✔
    ```

네! 실제로 Nextflow는 프로세스 호출이 정확히 동일하다는 것을 올바르게 감지하고, `-resume`으로 실행하고 있었으므로 모든 것을 다시 실행하지도 않습니다.

그리고 그게 전부입니다! 간단한 변이 호출 워크플로우에 원하는 모든 기본 기능이 있습니다.

### 핵심 정리

BAM 파일을 인덱싱하고 GATK를 사용하여 샘플별 변이 호출을 적용하는 다단계 모듈식 워크플로우를 만드는 방법을 알게 되었습니다.

더 일반적으로, 유전체학 파일 형식과 도구 요구 사항의 특이성을 고려하여 실제 작업을 수행하는 간단한 유전체학 파이프라인을 구축하기 위해 필수 Nextflow 구성 요소와 로직을 사용하는 방법을 배웠습니다.

### 다음 단계

성공을 축하하고 특별히 긴 휴식을 취하세요!

이 과정의 다음 파트에서는 이 간단한 샘플별 변이 호출 워크플로우를 데이터에 공동 변이 호출을 적용하도록 변환하는 방법을 배우게 됩니다.
