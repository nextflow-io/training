# 파트 3: 코드를 모듈로 이동하기

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } AI 지원 번역 - [자세히 알아보기 및 개선 제안](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

이 과정의 첫 번째 부분에서는 완전히 선형적이며 각 샘플의 데이터를 다른 샘플과 독립적으로 처리하는 변이 호출 파이프라인을 구축했습니다.

두 번째 부분에서는 채널과 채널 연산자를 사용하여 GATK로 공동 변이 호출을 구현하는 방법을 보여드렸으며, Part 1의 파이프라인을 기반으로 했습니다.

이 부분에서는 해당 워크플로우의 코드를 모듈로 변환하는 방법을 보여드리겠습니다. 이 교육 부분을 따라하려면 Part 1과 Part 2는 물론 모듈의 기본 사항을 다루는 [Hello Modules](../../../hello_nextflow/hello_modules.md)를 완료해야 합니다.

---

## 0. 준비 운동

워크플로우 개발을 시작할 때 모든 것을 하나의 코드 파일에 넣었습니다.
이제 코드를 **모듈화**할 시간입니다. 즉, 프로세스 정의를 모듈로 추출하는 것입니다.

Part 2와 동일한 워크플로우부터 시작하겠으며, `genomics-3.nf` 파일에 제공되어 있습니다.

!!! note "참고"

     올바른 작업 디렉토리에 있는지 확인하십시오:
     `cd /workspaces/training/nf4-science/genomics`

시작점을 확인하기 위해 워크플로우를 실행하십시오:

```bash
nextflow run genomics-3.nf -resume
```

```console title="출력"
 N E X T F L O W   ~  version 25.10.2

Launching `genomics-3.nf` [serene_borg] DSL2 - revision: 0cbebb67a1

executor >  local (7)
[6f/83ee72] SAMTOOLS_INDEX (3)       | 3 of 3 ✔
[53/b9d342] GATK_HAPLOTYPECALLER (1) | 3 of 3 ✔
[0c/fa6d15] GATK_JOINTGENOTYPING     | 1 of 1 ✔
```

이제 프로젝트 디렉토리 내에 `work` 디렉토리와 `results_genomics` 디렉토리가 생성됩니다.

### 요점 정리

워크플로우를 모듈화할 준비가 되었습니다.

### 다음 단계는?

Genomics 워크플로우의 프로세스를 모듈로 이동합니다.

---

## 1. 프로세스를 모듈로 이동하기

[Hello Modules](../../../hello_nextflow/hello_modules.md)에서 배운 것처럼, 프로세스 정의를 임의의 디렉토리에 있는 자체 파일로 복사하기만 하면 모듈을 만들 수 있으며, 해당 파일의 이름은 원하는 대로 지정할 수 있습니다.

나중에 명확해질 이유로 (특히 테스트를 할 때), 이 교육에서는 파일 이름을 `main.nf`로 지정하고 툴킷과 명령의 이름을 따서 명명된 디렉토리 구조에 배치하는 관례를 따르겠습니다.

### 1.1. `SAMTOOLS_INDEX` 프로세스를 위한 모듈 생성하기

`SAMTOOLS_INDEX` 프로세스의 경우, 'samtools'는 툴킷이고 'index'는 명령입니다. 따라서 `modules/samtools/index` 디렉토리 구조를 만들고 해당 디렉토리 내의 `main.nf` 파일에 `SAMTOOLS_INDEX` 프로세스 정의를 넣겠습니다.

```bash
mkdir -p modules/samtools/index
touch modules/samtools/index/main.nf
```

`main.nf` 파일을 열고 `SAMTOOLS_INDEX` 프로세스 정의를 복사합니다.

```groovy title="modules/samtools/index/main.nf" linenums="1"
/*
 * BAM 인덱스 파일 생성
 */
process SAMTOOLS_INDEX {

    container 'community.wave.seqera.io/library/samtools:1.20--b5dfbd93de237464'

    input:
    path input_bam

    output:
    tuple path(input_bam), path("${input_bam}.bai")

    script:
    """
    samtools index '$input_bam'
    """
}
```

그런 다음 `genomics-3.nf`에서 `SAMTOOLS_INDEX` 프로세스 정의를 제거하고, 다음 프로세스 정의 앞에 모듈에 대한 import 선언을 추가합니다:

=== "변경 후"

    ```groovy title="genomics-3.nf" linenums="1" hl_lines="1 2"
    // 모듈 포함
    include { SAMTOOLS_INDEX } from './modules/samtools/index/main.nf'

    /*
     * GATK HaplotypeCaller로 변이 호출
     */
    process GATK_HAPLOTYPECALLER {
    ```

=== "변경 전"

    ```groovy title="genomics-3.nf" linenums="1" hl_lines="1"
    /*
     * GATK HaplotypeCaller로 변이 호출
     */
    process GATK_HAPLOTYPECALLER {
    ```

이제 워크플로우를 다시 실행할 수 있으며, 이전과 동일한 방식으로 작동해야 합니다. `-resume` 플래그를 제공하면 새 작업을 실행할 필요조차 없어야 합니다:

```bash
nextflow run genomics-3.nf -resume
```

??? success "명령 출력"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `genomics-3.nf` [sleepy_snyder] DSL2 - revision: aa68d06c43

    [0f/71b55e] SAMTOOLS_INDEX (1)       | 3 of 3, cached: 3 ✔
    [f1/18971b] GATK_HAPLOTYPECALLER (3) | 3 of 3, cached: 3 ✔
    [0c/fa6d15] GATK_JOINTGENOTYPING     | 1 of 1, cached: 1 ✔
    ```

### 1.2. `GATK_HAPLOTYPECALLER` 및 `GATK_JOINTGENOTYPING` 프로세스를 위한 모듈 생성하기

나머지 프로세스에 대해서도 동일한 단계를 반복합니다.
각 프로세스에 대해:

1. 디렉토리 구조 생성 (`modules/gatk/haplotypecaller/` 및 `modules/gatk/jointgenotyping/`)
2. 프로세스 정의를 포함하는 `main.nf` 파일 생성
3. `genomics-3.nf`에서 프로세스 정의 제거
4. 모듈에 대한 import 선언 추가

완료되면 다음을 실행하여 modules 디렉토리 구조가 올바른지 확인하십시오:

```bash
tree modules/
```

??? abstract "디렉토리 내용"

    ```console
    modules/
    ├── gatk
    │   ├── haplotypecaller
    │   │   └── main.nf
    │   └── jointgenotyping
    │       └── main.nf
    └── samtools
        └── index
            └── main.nf

    5 directories, 3 files
    ```

또한 메인 워크플로우 파일의 매개변수 섹션 다음에 다음과 같은 내용이 있어야 합니다:

```
include { SAMTOOLS_INDEX } from './modules/samtools/index/main.nf'
include { GATK_HAPLOTYPECALLER } from './modules/gatk/haplotypecaller/main.nf'
include { GATK_JOINTGENOTYPING } from './modules/gatk/jointgenotyping/main.nf'

workflow {
```

### 요점 정리

Genomics 워크플로우를 예제로 워크플로우를 모듈화하는 연습을 했습니다.

### 다음 단계는?

모듈화된 워크플로우를 테스트합니다.

---

## 2. 모듈화된 워크플로우 테스트하기

모듈화된 워크플로우를 실행하여 모든 것이 여전히 작동하는지 확인합니다.

```bash
nextflow run genomics-3.nf -resume
```

```console title="출력"
 N E X T F L O W   ~  version 25.10.2

Launching `genomics-3.nf` [astonishing_venter] DSL2 - revision: ca27264c13

[6f/83ee72] SAMTOOLS_INDEX (3)       | 3 of 3, cached: 3 ✔
[53/b9d342] GATK_HAPLOTYPECALLER (3) | 3 of 3, cached: 3 ✔
[0c/fa6d15] GATK_JOINTGENOTYPING     | 1 of 1, cached: 1 ✔
```

파이프라인의 재개 가능성을 포함하여 모든 것이 여전히 작동합니다.
결과는 계속해서 `results_genomics` 디렉토리에 게시됩니다.

```console title="디렉토리 내용"
results_genomics/
├── family_trio.joint.vcf
├── family_trio.joint.vcf.idx
├── gvcf
│   ├── reads_father.bam.g.vcf
│   ├── reads_father.bam.g.vcf.idx
│   ├── reads_mother.bam.g.vcf
│   ├── reads_mother.bam.g.vcf.idx
│   ├── reads_son.bam.g.vcf
│   └── reads_son.bam.g.vcf.idx
└── indexed_bam
    ├── reads_father.bam
    ├── reads_father.bam.bai
    ├── reads_mother.bam
    ├── reads_mother.bam.bai
    ├── reads_son.bam
    └── reads_son.bam.bai
```

### 요점 정리

워크플로우를 모듈화했으며 이전과 동일한 방식으로 작동하는지 확인했습니다.

### 다음 단계는?

학습한 내용을 검토하고 테스트를 살펴봅니다.

---

## 3. 요약

워크플로우를 모듈화했으며 파이프라인 작동 방식에는 아무것도 변경되지 않았습니다.
이는 의도적인 것입니다: 기능에 영향을 주지 않고 코드를 재구성했습니다.

모듈은 프로세스 로직만 포함하여 깔끔하고 재사용 가능합니다.
메인 스크립트는 무엇이 어디에 게시되는지 제어하는 반면, 모듈은 계산 작업에 집중합니다.

코드를 더 쉽게 유지 관리할 수 있게 하는 기반을 마련했습니다.
예를 들어, 이제 nf-test 프레임워크를 사용하여 파이프라인에 테스트를 추가할 수 있습니다.
이것이 이 과정의 다음 부분에서 살펴볼 내용입니다.
