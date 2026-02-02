# 파트 1: 샘플별 변이 호출

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } AI 지원 번역 - [자세히 알아보기 및 개선 제안](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

이 과정의 첫 번째 파트에서는 개별 시퀀싱 샘플에 GATK 변이 호출을 적용하는 간단한 변이 호출 파이프라인을 구축하는 방법을 보여드립니다.

### 방법 개요

변이 호출은 참조 게놈에 대한 게놈 서열의 변이를 식별하는 것을 목표로 하는 게놈 분석 방법입니다.
여기서는 짧은 변이, 즉 SNP와 indel을 호출하기 위해 설계된 도구와 방법을 사용할 것입니다.

![GATK 파이프라인](img/gatk-pipeline.png)

완전한 변이 호출 파이프라인은 일반적으로 참조에 대한 매핑(게놈 정렬이라고도 함)과 변이 필터링 및 우선순위 지정을 포함한 많은 단계를 포함합니다.
간단하게 하기 위해 이 과정의 이 파트에서는 변이 호출 부분에만 집중하겠습니다.

### 데이터셋

다음과 같은 데이터와 관련 리소스를 제공합니다:

- 인간 염색체 20의 작은 영역(hg19/b37에서)과 그 보조 파일(인덱스 및 서열 사전)로 구성된 **참조 게놈**.
- 가족 3인조(어머니, 아버지, 아들)에 해당하는 **세 개의 전체 게놈 시퀀싱 샘플**로, 파일 크기를 작게 유지하기 위해 염색체 20의 작은 슬라이스 데이터로 부분 집합화되었습니다.
  이는 이미 참조 게놈에 매핑된 Illumina 단일 리드 시퀀싱 데이터이며, [BAM](https://samtools.github.io/hts-specs/SAMv1.pdf) 형식(Binary Alignment Map, SAM(Sequence Alignment Map)의 압축 버전)으로 제공됩니다.
- **게놈 간격 목록**, 즉 변이 호출에 적합한 데이터가 있는 샘플의 게놈 좌표로, BED 형식으로 제공됩니다.

### 워크플로우

이 과정의 이 파트에서는 다음을 수행하는 워크플로우를 개발할 것입니다:

1. [Samtools](https://www.htslib.org/)를 사용하여 각 BAM 입력 파일에 대한 인덱스 파일 생성
2. 각 BAM 입력 파일에 대해 GATK HaplotypeCaller를 실행하여 VCF(Variant Call Format)로 샘플별 변이 호출 생성

<figure class="excalidraw">
--8<-- "docs/en/docs/nf4_science/genomics/img/hello-gatk-1.svg"
</figure>

!!! note "참고"

    인덱스 파일은 생물정보학 파일 형식의 일반적인 기능입니다. 이들은 GATK와 같은 도구가 전체 파일을 읽지 않고도 데이터의 하위 집합에 액세스할 수 있도록 하는 메인 파일의 구조에 대한 정보를 포함합니다.
    이는 이러한 파일이 얼마나 커질 수 있는지 때문에 중요합니다.

---

## 0. 워밍업: Samtools 및 GATK 명령을 대화형으로 테스트합니다

먼저 워크플로우에 적용하기 전에 명령을 수동으로 시도해보고 싶습니다.
필요한 도구(Samtools 및 GATK)는 GitHub Codespaces 환경에 설치되어 있지 않으므로 컨테이너를 통해 사용할 것입니다([Hello Containers](../../hello_nextflow/05_hello_containers.md) 참조).

!!! note "참고"

     `pwd`를 입력할 때 표시되는 경로의 마지막 부분이 `genomics`가 되도록 `nf4-science/genomics` 디렉토리에 있는지 확인하십시오.

### 0.1. Samtools로 BAM 입력 파일 인덱싱

Samtools 컨테이너를 가져와서 대화형으로 실행하고 BAM 파일 중 하나에 `samtools index` 명령을 실행할 것입니다.

#### 0.1.1. Samtools 컨테이너 가져오기

```bash
docker pull community.wave.seqera.io/library/samtools:1.20--b5dfbd93de237464
```

<!--
??? success "명령 출력"

    ```console

    ```
-->

#### 0.1.2. Samtools 컨테이너를 대화형으로 실행

```bash
docker run -it -v ./data:/data community.wave.seqera.io/library/samtools:1.20--b5dfbd93de237464
```

<!--
??? success "명령 출력"

    ```console

    ```
-->

#### 0.1.3. 인덱싱 명령 실행

[Samtools 문서](https://www.htslib.org/doc/samtools-index.html)는 BAM 파일을 인덱싱하기 위해 실행할 명령줄을 제공합니다.

입력 파일만 제공하면 됩니다. 도구는 입력 파일 이름에 `.bai`를 추가하여 출력 이름을 자동으로 생성합니다.

```bash
samtools index /data/bam/reads_mother.bam
```

이것은 즉시 완료되어야 하며, 이제 원본 BAM 입력 파일과 같은 디렉토리에 `reads_mother.bam.bai`라는 파일이 표시되어야 합니다.

??? abstract "디렉토리 내용"

    ```console
    data/bam/
    ├── reads_father.bam
    ├── reads_mother.bam
    ├── reads_mother.bam.bai
    └── reads_son.bam
    ```

#### 0.1.4. Samtools 컨테이너 종료

```bash
exit
```

### 0.2. GATK HaplotypeCaller로 변이 호출

GATK 컨테이너를 가져와서 대화형으로 실행하고 방금 인덱싱한 BAM 파일에 `gatk HaplotypeCaller` 명령을 실행할 것입니다.

#### 0.2.1. GATK 컨테이너 가져오기

```bash
docker pull community.wave.seqera.io/library/gatk4:4.5.0.0--730ee8817e436867
```

<!--
??? success "명령 출력"

    ```console

    ```
-->

#### 0.2.2. GATK 컨테이너를 대화형으로 실행

```bash
docker run -it -v ./data:/data community.wave.seqera.io/library/gatk4:4.5.0.0--730ee8817e436867
```

<!--
??? success "명령 출력"

    ```console

    ```
-->

#### 0.2.3. 변이 호출 명령 실행

[GATK 문서](https://gatk.broadinstitute.org/hc/en-us/articles/21905025322523-HaplotypeCaller)는 BAM 파일에서 변이 호출을 수행하기 위해 실행할 명령줄을 제공합니다.

BAM 입력 파일(`-I`)과 참조 게놈(`-R`), 출력 파일 이름(`-O`) 및 분석할 게놈 간격 목록(`-L`)을 제공해야 합니다.

그러나 인덱스 파일의 경로를 지정할 필요는 없습니다. 도구는 확립된 명명 및 공동 배치 규칙에 따라 동일한 디렉토리에서 자동으로 찾습니다.
참조 게놈의 보조 파일(인덱스 및 서열 사전 파일, `*.fai` 및 `*.dict`)에도 동일하게 적용됩니다.

```bash
gatk HaplotypeCaller \
        -R /data/ref/ref.fasta \
        -I /data/bam/reads_mother.bam \
        -O reads_mother.vcf \
        -L /data/ref/intervals.bed
```

<!--
??? success "명령 출력"

    ```console

    ```
-->

출력 파일 `reads_mother.vcf`는 컨테이너 내의 작업 디렉토리에 생성되므로 출력 파일 경로를 변경하지 않는 한 VS Code 파일 탐색기에 표시되지 않습니다.
그러나 작은 테스트 파일이므로 `cat`으로 열어서 내용을 볼 수 있습니다.
파일의 시작 부분까지 스크롤하면 여러 줄의 메타데이터로 구성된 헤더가 있고, 그 다음에 한 줄에 하나씩 변이 호출 목록이 있습니다.

```console title="reads_mother.vcf" linenums="26"
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	reads_mother
20_10037292_10066351	3480	.	C	CT	503.03	.	AC=2;AF=1.00;AN=2;DP=23;ExcessHet=0.0000;FS=0.000;MLEAC=2;MLEAF=1.00;MQ=60.00;QD=27.95;SOR=1.179	GT:AD:DP:GQ:PL	1/1:0,18:18:54:517,54,0
20_10037292_10066351	3520	.	AT	A	609.03	.	AC=2;AF=1.00;AN=2;DP=18;ExcessHet=0.0000;FS=0.000;MLEAC=2;MLEAF=1.00;MQ=60.00;QD=33.83;SOR=0.693	GT:AD:DP:GQ:PL	1/1:0,18:18:54:623,54,0
20_10037292_10066351	3529	.	T	A	155.64	.	AC=1;AF=0.500;AN=2;BaseQRankSum=-0.544;DP=21;ExcessHet=0.0000;FS=1.871;MLEAC=1;MLEAF=0.500;MQ=60.00;MQRankSum=0.000;QD=7.78;ReadPosRankSum=-1.158;SOR=1.034	GT:AD:DP:GQ:PL	0/1:12,8:20:99:163,0,328
```

각 줄은 샘플의 시퀀싱 데이터에서 식별된 가능한 변이를 설명합니다. VCF 형식 해석에 대한 지침은 [이 유용한 문서](https://www.ebi.ac.uk/training/online/courses/human-genetic-variation-introduction/variant-identification-and-analysis/understanding-vcf-format/)를 참조하십시오.

출력 VCF 파일에는 GATK가 자동으로 생성한 `reads_mother.vcf.idx`라는 인덱스 파일이 함께 제공됩니다.
이것은 BAM 인덱스 파일과 같은 기능을 하며, 도구가 전체 파일을 로드하지 않고도 데이터의 하위 집합을 찾아서 검색할 수 있도록 합니다.

#### 0.2.4. GATK 컨테이너 종료

```bash
exit
```

### 핵심 사항

Samtools 인덱싱 및 GATK 변이 호출 명령을 각각의 컨테이너에서 테스트하는 방법을 알게 되었습니다.

### 다음은?

동일한 명령을 컨테이너를 사용하여 작업을 실행하는 2단계 워크플로우로 적용하는 방법을 배웁니다.

---

## 1. BAM 파일에 Samtools index를 실행하는 단일 단계 워크플로우 작성

워크플로우의 주요 부분을 개괄하는 워크플로우 파일 `genomics-1.nf`를 제공합니다.
기능적이지 않습니다. 그 목적은 실제 워크플로우를 작성하는 데 사용할 골격 역할을 하는 것입니다.

### 1.1. 인덱싱 프로세스 정의

인덱싱 작업을 설명하는 `SAMTOOLS_INDEX`라는 프로세스를 작성하는 것부터 시작하겠습니다.

```groovy title="genomics-1.nf" linenums="9"
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

이 교육 시리즈의 파트 1 & 파트 2에서 배운 모든 요소를 인식할 수 있을 것입니다.

이 프로세스는 `input_bam` 입력을 통해 파일 경로를 전달해야 하므로 다음에 설정하겠습니다.

### 1.2. 입력 매개변수 선언 추가

파일 상단의 `Pipeline parameters` 섹션 아래에 `reads_bam`이라는 CLI 매개변수를 선언하고 기본값을 지정합니다.
이렇게 하면 게으르게 파이프라인을 시작하는 명령을 입력할 때 입력을 지정하지 않을 수 있습니다(개발 목적).

```groovy title="genomics-1.nf" linenums="3"
/*
 * Pipeline parameters
 */
params {
    // 기본 입력
    reads_bam: Path = "${projectDir}/data/bam/reads_mother.bam"
}
```

이제 프로세스가 준비되었고 실행할 입력을 제공할 매개변수도 있으므로 이들을 함께 연결해 봅시다.

!!! note "참고"

    `${projectDir}`는 현재 Nextflow 워크플로우 스크립트(`genomics-1.nf`)가 위치한 디렉토리를 가리키는 내장 Nextflow 변수입니다.

    이를 통해 절대 경로를 하드코딩하지 않고도 워크플로우 리포지토리에 포함된 파일, 데이터 디렉토리 및 기타 리소스를 쉽게 참조할 수 있습니다.

### 1.3. SAMTOOLS_INDEX를 실행하기 위한 workflow 블록 추가

`workflow` 블록에서는 `SAMTOOLS_INDEX` 프로세스에 입력을 공급하기 위한 **채널**을 설정해야 합니다. 그런 다음 프로세스 자체를 호출하여 해당 채널의 내용에 대해 실행할 수 있습니다.

```groovy title="genomics-1.nf" linenums="24"
workflow {

    main:
    // 입력 채널 생성 (CLI 매개변수를 통한 단일 파일)
    reads_ch = channel.fromPath(params.reads_bam)

    // 입력 BAM 파일에 대한 인덱스 파일 생성
    SAMTOOLS_INDEX(reads_ch)

    publish:
    bam_index = SAMTOOLS_INDEX.out
}
```

workflow 블록에는 두 개의 섹션이 있습니다:

- `main:`에는 채널 작업 및 프로세스 호출이 포함됩니다
- `publish:`는 게시해야 하는 출력을 선언하고 명명된 타겟에 할당합니다

[Hello Channels](../../hello_nextflow/02_hello_channels.md)에서 사용한 것과 동일한 `.fromPath` channel factory를 사용하고 있음을 알 수 있습니다.
실제로 매우 유사한 작업을 하고 있습니다.
차이점은 Nextflow에게 내용을 읽어들이는 대신 파일 경로 자체를 입력 요소로 채널에 로드하도록 지시한다는 것입니다.

### 1.4. 결과가 게시되는 위치를 정의하기 위한 output 블록 추가

workflow 블록 이후에 워크플로우 출력을 게시할 위치를 지정하는 `output` 블록을 추가합니다.

```groovy title="genomics-1.nf" linenums="37"
output {
    bam_index {
        path '.'
    }
}
```

`publish:` 섹션의 각 명명된 타겟(예: `bam_index`)은 기본 출력 디렉토리에 상대적인 출력 경로를 구성할 수 있는 자체 블록을 갖습니다.

!!! note "참고"

    여기서 사용하는 데이터 파일은 매우 작지만 게놈학에서는 매우 클 수 있습니다.
    기본적으로 Nextflow는 게시 디렉토리의 출력 파일에 대한 심볼릭 링크를 생성하여 불필요한 파일 복사를 방지합니다.
    `mode` 옵션(예: `mode 'copy'`)을 사용하여 이 동작을 변경하여 실제 복사본을 대신 생성할 수 있습니다.
    `work` 디렉토리를 정리하면 심볼릭 링크가 끊어지므로 프로덕션 워크플로우의 경우 `mode 'copy'`를 사용할 수 있습니다.

### 1.5. 출력 디렉토리 구성

기본 출력 디렉토리는 `outputDir` 구성 옵션을 통해 설정됩니다. `nextflow.config`에 추가하십시오:

=== "변경 후"

    ```groovy title="nextflow.config" hl_lines="2"
    docker.enabled = true
    outputDir = 'results_genomics'
    ```

=== "변경 전"

    ```groovy title="nextflow.config"
    docker.enabled = true
    ```

### 1.6. 워크플로우를 실행하여 인덱싱 단계가 작동하는지 확인

워크플로우를 실행해 봅시다! 상기시키자면, 입력 매개변수를 선언할 때 입력의 기본값을 설정했기 때문에 명령줄에서 입력을 지정할 필요가 없습니다.

```bash
nextflow run genomics-1.nf
```

??? success "명령 출력"

    ```console
    N E X T F L O W   ~  version 25.10.2

    ┃ Launching `genomics-1.nf` [reverent_sinoussi] DSL2 - revision: 41d43ad7fe

    executor >  local (1)
    [2a/e69536] SAMTOOLS_INDEX (1) | 1 of 1 ✔
    ```

작업 디렉토리 또는 결과 디렉토리를 살펴보면 인덱스 파일이 올바르게 생성되었는지 확인할 수 있습니다.

??? abstract "작업 디렉토리 내용"

    ```console
    work/2a/e695367b2f60df09cf826b07192dc3
    ├── reads_mother.bam -> /workspaces/training/nf4-science/genomics/data/bam/reads_mother.bam
    └── reads_mother.bam.bai
    ```

??? abstract "결과 디렉토리 내용"

    ```console
    results_genomics/
    └── reads_mother.bam.bai -> */87/908bba*/reads_mother.bam.bai
    ```

여기 있습니다!

### 핵심 사항

게놈학 도구를 단일 단계 Nextflow 워크플로우로 적용하고 컨테이너를 사용하여 실행하는 방법을 알게 되었습니다.

### 다음은?

첫 번째 단계의 출력을 소비하는 두 번째 단계를 추가합니다.

---

## 2. 인덱싱된 BAM 파일에 GATK HaplotypeCaller를 실행하는 두 번째 프로세스 추가

이제 입력 파일에 대한 인덱스가 있으므로 워크플로우의 흥미로운 부분인 변이 호출 단계 설정으로 넘어갈 수 있습니다.

### 2.1. 변이 호출 프로세스 정의

변이 호출 작업을 설명하는 `GATK_HAPLOTYPECALLER`라는 프로세스를 작성하겠습니다.

```groovy title="genomics-1.nf" linenums="44"
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

여기서 각 출력 채널에 고유한 이름을 지정하기 위해 새로운 구문(`emit:`)을 도입했으며, 그 이유는 곧 명확해질 것입니다.

이 명령은 간단한 인덱싱 작업에 비해 GATK가 분석을 수행하기 위해 더 많은 정보가 필요하기 때문에 훨씬 더 많은 입력을 사용합니다.
그러나 입력 블록에 정의된 입력이 GATK 명령에 나열된 것보다 훨씬 더 많다는 것을 알 수 있습니다. 왜 그럴까요?

!!! note "참고"

    GATK는 해당 파일을 둘러싼 규칙을 인식하기 때문에 BAM 인덱스 파일과 참조 게놈의 보조 파일을 찾을 줄 압니다.
    그러나 Nextflow는 도메인에 구애받지 않도록 설계되었으며 생물정보학 파일 형식 요구 사항에 대해 아무것도 알지 못합니다.

런타임에 작업 디렉토리에 해당 파일을 스테이징해야 한다고 Nextflow에 명시적으로 알려야 합니다. 그렇지 않으면 수행하지 않으며 GATK는 인덱스 파일이 누락되었다는 오류를 (올바르게) 발생시킵니다.

마찬가지로 Nextflow가 후속 단계에서 필요한 경우 해당 파일을 추적하도록 출력 VCF의 인덱스 파일(`"${input_bam}.vcf.idx"` 파일)을 명시적으로 나열해야 합니다.

### 2.2. 보조 입력에 대한 정의 추가

새 프로세스는 추가 파일이 제공될 것으로 예상하므로 `Pipeline parameters` 섹션 아래에 일부 기본값과 함께 CLI 매개변수를 설정합니다(이전과 같은 이유).

```groovy title="genomics-1.nf" linenums="8"
    // 보조 파일
    reference: Path = "${projectDir}/data/ref/ref.fasta"
    reference_index: Path = "${projectDir}/data/ref/ref.fasta.fai"
    reference_dict: Path = "${projectDir}/data/ref/ref.dict"
    intervals: Path = "${projectDir}/data/ref/intervals.bed"
```

### 2.3. 보조 파일 경로를 보유할 변수 생성

메인 데이터 입력은 채널을 통해 동적으로 스트리밍되지만 보조 파일을 처리하는 두 가지 접근 방식이 있습니다. 권장되는 접근 방식은 명시적 채널을 생성하는 것으로, 데이터 흐름을 더 명확하고 일관되게 만듭니다. 또는 더 간단한 경우, 특히 여러 프로세스에서 동일한 파일을 참조해야 할 때 변수를 생성하기 위해 file() 함수를 사용할 수 있습니다. 다만 이것도 암묵적으로 채널을 생성한다는 점에 유의하십시오. <!-- TODO: Clarify: is this still necessary with typed inputs? -->

workflow 블록에 추가하십시오(`main:` 섹션 내부의 `reads_ch` 생성 후):

```groovy title="genomics-1.nf" linenums="79"
    // 보조 파일(참조 및 인터벌)에 대한 파일 경로 로드
    ref_file        = file(params.reference)
    ref_index_file  = file(params.reference_index)
    ref_dict_file   = file(params.reference_dict)
    intervals_file  = file(params.intervals)
```

이렇게 하면 보조 파일 경로를 필요한 프로세스에 입력으로 제공할 수 있습니다.

### 2.4. GATK_HAPLOTYPECALLER를 실행하기 위한 workflow 블록에 호출 추가

이제 두 번째 프로세스를 설정하고 모든 입력과 보조 파일이 준비되고 사용 가능하므로 워크플로우 본문에 `GATK_HAPLOTYPECALLER` 프로세스 호출을 추가할 수 있습니다.

```groovy title="genomics-1.nf" linenums="88"
    // 인덱싱된 BAM 파일에서 변이 호출
    GATK_HAPLOTYPECALLER(
        reads_ch,
        SAMTOOLS_INDEX.out,
        ref_file,
        ref_index_file,
        ref_dict_file,
        intervals_file
    )
```

이 교육 시리즈의 파트 1에서 `*.out` 구문을 인식해야 합니다. Nextflow에게 `SAMTOOLS_INDEX`가 출력한 채널을 가져와서 `GATK_HAPLOTYPECALLER` 프로세스 호출에 연결하도록 지시하고 있습니다.

!!! note "참고"

    입력이 프로세스 호출에서 제공되는 순서가 프로세스의 input 블록에 나열된 순서와 정확히 동일하다는 것을 알 수 있습니다.
    Nextflow에서 입력은 위치 지정적입니다. 즉, 동일한 순서를 _따라야_ 합니다. 그리고 물론 동일한 수의 요소가 있어야 합니다.

### 2.5. publish 섹션 및 output 블록 업데이트

VCF 출력을 포함하도록 `publish:` 섹션을 업데이트하고 `output` 블록에 해당 타겟을 추가해야 합니다.

```groovy title="genomics-1.nf" linenums="99"
    publish:
    bam_index = SAMTOOLS_INDEX.out
    vcf = GATK_HAPLOTYPECALLER.out.vcf
    vcf_idx = GATK_HAPLOTYPECALLER.out.idx
}

output {
    bam_index {
        path '.'
    }
    vcf {
        path '.'
    }
    vcf_idx {
        path '.'
    }
}
```

### 2.6. 워크플로우를 실행하여 변이 호출 단계가 작동하는지 확인

인덱싱 단계를 다시 실행할 필요가 없도록 `-resume`으로 확장된 워크플로우를 실행해 봅시다.

```bash
nextflow run genomics-1.nf -resume
```

??? success "명령 출력"

    ```console
    N E X T F L O W   ~  version 25.10.2

    ┃ Launching `genomics-1.nf` [grave_volta] DSL2 - revision: 4790abc96a

    executor >  local (1)
    [2a/e69536] SAMTOOLS_INDEX (1)       | 1 of 1, cached: 1 ✔
    [53/e18e98] GATK_HAPLOTYPECALLER (1) | 1 of 1 ✔
    ```

이제 콘솔 출력을 보면 두 프로세스가 나열됩니다.

예상대로 캐싱 덕분에 첫 번째 프로세스는 건너뛰었고, 두 번째 프로세스는 완전히 새로운 것이므로 실행되었습니다.

결과 디렉토리에서 출력 파일을 찾을 수 있습니다(작업 디렉토리에 대한 심볼릭 링크로).

??? abstract "디렉토리 내용"

    ```console
    results_genomics/
    ├── reads_mother.bam.bai -> */87/908bba*/reads_mother.bam.bai
    ├── reads_mother.bam.vcf -> */cf/36f756*/reads_mother.bam.vcf
    └── reads_mother.bam.vcf.idx -> */cf/36f756*/reads_mother.bam.vcf.idx
    ```

VCF 파일을 열면 컨테이너에서 GATK 명령을 직접 실행하여 생성한 파일과 동일한 내용을 볼 수 있습니다.

```console title="reads_mother.bam.vcf" linenums="26"
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	reads_mother
20_10037292_10066351	3480	.	C	CT	503.03	.	AC=2;AF=1.00;AN=2;DP=23;ExcessHet=0.0000;FS=0.000;MLEAC=2;MLEAF=1.00;MQ=60.00;QD=27.95;SOR=1.179	GT:AD:DP:GQ:PL	1/1:0,18:18:54:517,54,0
20_10037292_10066351	3520	.	AT	A	609.03	.	AC=2;AF=1.00;AN=2;DP=18;ExcessHet=0.0000;FS=0.000;MLEAC=2;MLEAF=1.00;MQ=60.00;QD=33.83;SOR=0.693	GT:AD:DP:GQ:PL	1/1:0,18:18:54:623,54,0
20_10037292_10066351	3529	.	T	A	155.64	.	AC=1;AF=0.500;AN=2;BaseQRankSum=-0.544;DP=21;ExcessHet=0.0000;FS=1.871;MLEAC=1;MLEAF=0.500;MQ=60.00;MQRankSum=0.000;QD=7.78;ReadPosRankSum=-1.158;SOR=1.034	GT:AD:DP:GQ:PL	0/1:12,8:20:99:163,0,328
```

이것이 연구의 각 샘플에 대해 생성하려는 출력입니다.

### 핵심 사항

보조 파일과 같은 게놈학 파일 형식의 특이성을 처리할 수 있는 실제 분석 작업을 수행하는 매우 기본적인 2단계 워크플로우를 만드는 방법을 알게 되었습니다.

### 다음은?

워크플로우가 여러 샘플을 일괄 처리하도록 만듭니다.

---

## 3. 샘플 배치에서 실행되도록 워크플로우 조정

단일 샘플에 대한 처리를 자동화할 수 있는 워크플로우를 갖는 것은 좋지만 1000개의 샘플이 있으면 어떻게 될까요?
모든 샘플을 반복하는 bash 스크립트를 작성해야 합니까?

아니요, 다행히도! 코드를 약간만 수정하면 Nextflow가 이를 처리해 줄 것입니다.

### 3.1. 입력 매개변수 선언을 세 개의 샘플을 나열하는 배열로 변환

`Pipeline parameters` 섹션 아래의 입력 BAM 파일 선언에서 해당 기본 파일 경로를 세 개의 테스트 샘플에 대한 파일 경로를 나열하는 배열로 바꿔 봅시다.

=== "변경 후"

    ```groovy title="genomics-1.nf" linenums="7"
    // 기본 입력 (세 샘플의 배열)
        reads_bam = [
            "${projectDir}/data/bam/reads_mother.bam",
            "${projectDir}/data/bam/reads_father.bam",
            "${projectDir}/data/bam/reads_son.bam"
        ]
    ```

=== "변경 전"

    ```groovy title="genomics-1.nf" linenums="7"
        // 기본 입력
        reads_bam: Path = "${projectDir}/data/bam/reads_mother.bam"
    ```

!!! note "참고"

    타입이 지정된 매개변수 선언(예: `reads_bam: Path`)을 사용하는 경우 배열 값을 할당할 수 없습니다.
    배열의 경우 타입 주석을 생략하십시오.

그리고 실제로 우리가 해야 할 일은 이것뿐입니다. 워크플로우 본문에서 사용하는 channel factory(`.fromPath`)는 입력 채널에 로드할 여러 파일 경로를 받는 것을 단일 파일 경로를 로드하는 것만큼 기쁘게 여기기 때문입니다.

!!! note "참고"

    일반적으로 샘플 목록을 워크플로우 파일에 하드코딩하고 싶지 않지만 여기서는 단순하게 유지하기 위해 그렇게 하고 있습니다.
    이 교육 시리즈의 후반부에서 입력을 처리하는 더 우아한 방법을 제시할 것입니다.

### 3.2. 워크플로우를 실행하여 세 샘플 모두에서 실행되는지 확인

이제 세 개의 테스트 샘플 모두에서 실행되도록 배관이 설정되었으니 워크플로우를 실행해 봅시다.

```bash
nextflow run genomics-1.nf -resume
```

재미있는 점: 이것이 _작동할 수도_ 있고 _실패할 수도_ 있습니다. 예를 들어 다음은 성공한 실행입니다:

??? success "명령 출력"

    ```console
    N E X T F L O W   ~  version 25.10.2

    ┃ Launching `genomics-1.nf` [peaceful_yalow] DSL2 - revision: a256d113ad

    executor >  local (6)
    [4f/7071b0] SAMTOOLS_INDEX (3)       | 3 of 3, cached: 1 ✔
    [7a/89bc43] GATK_HAPLOTYPECALLER (2) | 3 of 3, cached: 1 ✔
    ```

워크플로우 실행이 성공했다면 다음과 같은 오류가 발생할 때까지 다시 실행하십시오:

??? failure "명령 출력"

    ```console
    N E X T F L O W   ~  version 25.10.2

    ┃ Launching `genomics-1.nf` [loving_pasteur] DSL2 - revision: d2a8e63076

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

글쎄요, 워크플로우의 첫 번째 단계에서 BAM 파일을 명시적으로 인덱싱했다는 점을 고려하면 이상합니다. 배관에 문제가 있을까요?

#### 3.2.1. 관련 호출에 대한 작업 디렉토리 확인

콘솔 출력에 나열된 실패한 `GATK_HAPLOTYPECALLER` 프로세스 호출의 작업 디렉토리 내부를 살펴보겠습니다.

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

이 디렉토리에 나열된 BAM 파일과 BAM 인덱스의 이름에 특히 주의하십시오: `reads_son.bam` 및 `reads_father.bam.bai`.

뭐라고요? Nextflow가 이 프로세스 호출의 작업 디렉토리에 인덱스 파일을 스테이징했지만 잘못된 것입니다. 어떻게 이런 일이 일어날 수 있었을까요?

#### 3.2.2. [view() 연산자](https://www.nextflow.io/docs/latest/reference/operator.html#view)를 사용하여 채널 내용 검사

`GATK_HAPLOTYPER` 프로세스 호출 전에 workflow 본문에 이 두 줄을 추가하십시오:

```groovy title="genomics-1.nf" linenums="84"
    // 임시 진단
    reads_ch.view()
    SAMTOOLS_INDEX.out.view()
```

그런 다음 워크플로우 명령을 다시 실행하십시오.

```bash
nextflow run genomics-1.nf
```

다시 한 번, 이것이 성공하거나 실패할 수 있습니다. 다음은 성공적인 실행입니다:

??? success "명령 출력"

    ```console
    N E X T F L O W   ~  version 25.10.2

    ┃ Launching `genomics-1.nf` [fervent_pasteur] DSL2 - revision: a256d113ad

    /workspaces/training/nf4-science/genomics/data/bam/reads_mother.bam
    /workspaces/training/nf4-science/genomics/data/bam/reads_father.bam
    /workspaces/training/nf4-science/genomics/data/bam/reads_son.bam
    executor >  local (6)
    [4f/7071b0] SAMTOOLS_INDEX (3)       | 3 of 3 ✔
    /workspaces/training/nf4-science/genomics/work/b4/45a376f0e724be1dc626a6807f73d8/reads_mother.bam.bai
    /workspaces/training/nf4-science/genomics/work/4f/7071b082b45dd85b1c9b6b3b32cb69/reads_father.bam.bai
    /workspaces/training/nf4-science/genomics/work/3c/331645a9e20e67edae10da5ba17c7b/reads_son.bam.bai
    [a2/dbd8d5] GATK_HAPLOTYPECALLER (3) | 3 of 3 ✔
    ```

그리고 다음은 실패한 것입니다:

??? failure "명령 출력"

    ```console
    N E X T F L O W   ~  version 25.10.2

    ┃ Launching `genomics-1.nf` [angry_hamilton] DSL2 - revision: a256d113ad

    /workspaces/training/nf4-science/genomics/data/bam/reads_mother.bam
    /workspaces/training/nf4-science/genomics/data/bam/reads_father.bam
    /workspaces/training/nf4-science/genomics/data/bam/reads_son.bam
    executor >  local (6)
    [4f/7071b0] SAMTOOLS_INDEX (3)       | 3 of 3 ✔
    /workspaces/training/nf4-science/genomics/work/4f/7071b082b45dd85b1c9b6b3b32cb69/reads_father.bam.bai
    /workspaces/training/nf4-science/genomics/work/3c/331645a9e20e67edae10da5ba17c7b/reads_son.bam.bai
    /workspaces/training/nf4-science/genomics/work/b4/45a376f0e724be1dc626a6807f73d8/reads_mother.bam.bai
    [a3/cf3a89] GATK_HAPLOTYPECALLER (3) | 1 of 3
    ERROR ~ Error executing process > 'GATK_HAPLOTYPECALLER (3)'
    ...
    ```

다시 실패하려면 여러 번 실행해야 할 수 있습니다.
이 오류는 개별 프로세스 호출의 실행 시간에 일부 변동성에 의존하기 때문에 일관되게 재현되지 않습니다.

다음은 실패한 실행에 대해 추가한 두 개의 `.view()` 호출의 출력이 어떻게 보이는지입니다:

```console
/workspaces/training/nf4-science/genomics/data/bam/reads_mother.bam
/workspaces/training/nf4-science/genomics/data/bam/reads_father.bam
/workspaces/training/nf4-science/genomics/data/bam/reads_son.bam
/workspaces/training/nf4-science/genomics/work/9c/53492e3518447b75363e1cd951be4b/reads_father.bam.bai
/workspaces/training/nf4-science/genomics/work/cc/37894fffdf6cc84c3b0b47f9b536b7/reads_son.bam.bai
/workspaces/training/nf4-science/genomics/work/4d/dff681a3d137ba7d9866e3d9307bd0/reads_mother.bam.bai
```

처음 세 줄은 입력 채널에 해당하고 두 번째는 출력 채널에 해당합니다.
세 샘플의 BAM 파일과 인덱스 파일이 같은 순서로 나열되지 않았음을 알 수 있습니다!

!!! note "참고"

    여러 요소를 포함하는 채널에서 Nextflow 프로세스를 호출하면 Nextflow는 가능한 한 많이 병렬로 실행하려고 시도하고 사용 가능해지는 순서대로 출력을 수집합니다.
    결과적으로 해당 출력은 원래 입력이 제공된 것과 다른 순서로 수집될 수 있습니다.

현재 작성된 대로, 우리의 워크플로우 스크립트는 인덱스 파일이 입력이 제공된 것과 동일한 어머니/아버지/아들 순서로 인덱싱 단계에서 나올 것으로 가정합니다.
그러나 이것은 보장되지 않으며, 그래서 때때로(항상 그런 것은 아니지만) 잘못된 파일이 두 번째 단계에서 짝지어집니다.

이를 해결하려면 BAM 파일과 해당 인덱스 파일이 채널을 통해 함께 이동하도록 해야 합니다.

!!! tip "팁"

    워크플로우 코드의 `view()` 문은 아무것도 하지 않으므로 그대로 두는 것은 문제가 되지 않습니다.
    그러나 콘솔 출력이 복잡해지므로 문제 해결을 마치면 제거하는 것이 좋습니다.

### 3.3. SAMTOOLS_INDEX 프로세스의 출력을 입력 파일과 인덱스를 함께 유지하는 튜플로 변경

BAM 파일과 인덱스가 밀접하게 연관되도록 하는 가장 간단한 방법은 인덱스 작업에서 나오는 튜플로 함께 패키징하는 것입니다.

!!! note "참고"

    **튜플**은 함수에서 여러 값을 반환하는 데 일반적으로 사용되는 유한하고 순서가 지정된 요소 목록입니다. 튜플은 연관성과 순서를 유지하면서 프로세스 간에 여러 입력 또는 출력을 전달하는 데 특히 유용합니다.

먼저, `SAMTOOLS_INDEX` 프로세스의 출력을 출력 선언에 BAM 파일을 포함하도록 변경하겠습니다.

=== "변경 후"

    ```groovy title="genomics-1.nf" linenums="20"
        output:
        tuple path(input_bam), path("${input_bam}.bai")
    ```

=== "변경 전"

    ```groovy title="genomics-1.nf" linenums="20"
        output:
        path "${input_bam}.bai"
    ```

이렇게 하면 각 인덱스 파일이 원본 BAM 파일과 긴밀하게 결합되고 인덱싱 단계의 전체 출력은 파일 쌍을 포함하는 단일 채널이 됩니다.

### 3.4. GATK_HAPLOTYPECALLER 프로세스의 입력을 튜플로 변경

워크플로우의 첫 번째 프로세스 출력의 '형태'를 변경했으므로 일치하도록 두 번째 프로세스의 입력 정의를 업데이트해야 합니다.

특히 이전에 `GATK_HAPLOTYPECALLER` 프로세스의 input 블록에서 두 개의 별도 입력 경로를 선언했던 곳에서 이제 `SAMTOOLS_INDEX`가 방출한 튜플의 구조와 일치하는 단일 입력을 선언합니다.

=== "변경 후"

    ```groovy title="genomics-1.nf" linenums="51"
        input:
        tuple path(input_bam), path(input_bam_index)
    ```

=== "변경 전"

    ```groovy title="genomics-1.nf" linenums="51"
        input:
        path input_bam
        path input_bam_index
    ```

물론 이제 `GATK_HAPLOTYPECALLER`가 기대하는 입력의 형태를 변경했으므로 워크플로우 본문의 프로세스 호출을 그에 따라 업데이트해야 합니다.

### 3.5. workflow 블록에서 GATK_HAPLOTYPECALLER 호출 업데이트

BAM 파일이 이제 `SAMTOOLS_INDEX`의 채널 출력에 번들로 제공되므로 더 이상 원래 `reads_ch`를 `GATK_HAPLOTYPECALLER` 프로세스에 제공할 필요가 없습니다.

결과적으로 해당 줄을 삭제하기만 하면 됩니다.

=== "변경 후"

    ```groovy title="genomics-1.nf" linenums="88"
    GATK_HAPLOTYPECALLER(
        SAMTOOLS_INDEX.out,
    ```

=== "변경 전"

    ```groovy title="genomics-1.nf" linenums="88"
    GATK_HAPLOTYPECALLER(
        reads_ch,
        SAMTOOLS_INDEX.out,
    ```

인덱스 불일치 문제를 해결하기 위해 필요한 재배선은 이것뿐입니다.

### 3.6. 튜플에 대한 publish 섹션 및 output 블록 업데이트

`SAMTOOLS_INDEX.out`은 이제 BAM과 인덱스를 모두 포함하는 튜플이므로 두 파일이 함께 게시됩니다.
타겟 이름을 `bam_index`에서 `indexed_bam`으로 변경하여 이제 두 파일을 모두 포함한다는 것을 반영합니다.

=== "변경 후"

    ```groovy title="genomics-1.nf" hl_lines="2"
        publish:
        indexed_bam = SAMTOOLS_INDEX.out
    ```

=== "변경 전"

    ```groovy title="genomics-1.nf"
        publish:
        bam_index = SAMTOOLS_INDEX.out
    ```

또한 새 타겟 이름을 사용하도록 output 블록을 업데이트해야 합니다:

=== "변경 후"

    ```groovy title="genomics-1.nf" hl_lines="2"
    output {
        indexed_bam {
            path '.'
        }
    ```

=== "변경 전"

    ```groovy title="genomics-1.nf"
    output {
        bam_index {
            path '.'
        }
    ```

### 3.7. 워크플로우를 실행하여 매번 세 샘플 모두에서 올바르게 작동하는지 확인

물론 증명은 푸딩에 있으므로 앞으로 안정적으로 작동하는지 확인하기 위해 워크플로우를 몇 번 더 실행해 봅시다.

```bash
nextflow run genomics-1.nf
```

이번에는(그리고 매번) 모든 것이 올바르게 실행되어야 합니다:

??? success "명령 출력"

    ```console
    N E X T F L O W   ~  version 25.10.2

    ┃ Launching `genomics-1.nf` [special_goldstine] DSL2 - revision: 4cbbf6ea3e

    executor >  local (6)
    [d6/10c2c4] SAMTOOLS_INDEX (1)       | 3 of 3 ✔
    [88/1783aa] GATK_HAPLOTYPECALLER (2) | 3 of 3 ✔
    ```

결과 디렉토리에는 이제 VCF 출력과 함께 각 샘플에 대한 BAM 및 BAI 파일(튜플에서)이 모두 포함됩니다:

??? abstract "결과 디렉토리 내용"

    ```console
    results_genomics/
    ├── reads_father.bam -> */60/e2614c*/reads_father.bam
    ├── reads_father.bam.bai -> */60/e2614c*/reads_father.bam.bai
    ├── reads_father.bam.vcf -> */b8/91b3c8*/reads_father.bam.vcf
    ├── reads_father.bam.vcf.idx -> */b8/91b3c8*/reads_father.bam.vcf.idx
    ├── reads_mother.bam -> */3e/fededc*/reads_mother.bam
    ├── reads_mother.bam.bai -> */3e/fededc*/reads_mother.bam.bai
    ├── reads_mother.bam.vcf -> */32/5ca037*/reads_mother.bam.vcf
    ├── reads_mother.bam.vcf.idx -> */32/5ca037*/reads_mother.bam.vcf.idx
    ├── reads_son.bam -> */3c/36d1c2*/reads_son.bam
    ├── reads_son.bam.bai -> */3c/36d1c2*/reads_son.bam.bai
    ├── reads_son.bam.vcf -> */d7/a6b046*/reads_son.bam.vcf
    └── reads_son.bam.vcf.idx -> */d7/a6b046*/reads_son.bam.vcf.idx
    ```

원한다면 `.view()`를 다시 사용하여 `SAMTOOLS_INDEX` 출력 채널의 내용이 어떻게
