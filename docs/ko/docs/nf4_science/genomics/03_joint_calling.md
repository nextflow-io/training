# Part 3: 코호트에 대한 합동 호출

Part 2에서는 각 샘플의 데이터를 독립적으로 처리하는 샘플별 변이 호출 파이프라인을 구축했습니다.
이제 [Part 1](01_method.md)에서 다룬 합동 변이 호출을 구현하도록 파이프라인을 확장하겠습니다.

## 과제

이 과정의 이 부분에서는 다음을 수행하도록 워크플로우를 확장할 것입니다:

<figure class="excalidraw">
--8<-- "docs/en/docs/nf4_science/genomics/img/hello-gatk-2.svg"
</figure>

1. Samtools를 사용하여 각 BAM 입력 파일의 인덱스 파일 생성
2. 각 BAM 입력 파일에 대해 GATK HaplotypeCaller를 실행하여 샘플별 게놈 변이 호출의 GVCF 생성
3. 모든 GVCF를 수집하고 GenomicsDB 데이터 저장소로 결합
4. 결합된 GVCF 데이터 저장소에 대해 합동 유전자형 분석을 실행하여 코호트 수준 VCF 생성

이 부분은 Part 2에서 생성한 워크플로우를 직접 기반으로 합니다.

??? info "이 섹션부터 시작하는 방법"

    이 섹션은 [Part 2: 샘플별 변이 호출](./02_per_sample_variant_calling.md)을 완료하고 작동하는 `genomics.nf` 파이프라인을 가지고 있다고 가정합니다.

    Part 2를 완료하지 않았거나 이 부분을 새로 시작하고 싶다면, Part 2 해결책을 시작점으로 사용할 수 있습니다.
    `nf4-science/genomics/` 디렉토리 내부에서 다음 명령을 실행하세요:

    ```bash
    cp solutions/part2/genomics-2.nf genomics.nf
    cp solutions/part2/nextflow.config .
    cp solutions/part2/modules/* modules/
    ```

    이렇게 하면 완전한 샘플별 변이 호출 워크플로우를 얻을 수 있습니다.
    다음 명령을 실행하여 성공적으로 실행되는지 테스트할 수 있습니다:

    ```bash
    nextflow run genomics.nf -profile test
    ```

## 학습 계획

이를 두 단계로 나누었습니다:

1. **샘플별 변이 호출 단계를 GVCF를 생성하도록 수정합니다.**
   이는 프로세스 명령과 출력 업데이트를 다룹니다.
2. **샘플별 GVCF를 결합하고 유전자형을 분석하는 합동 유전자형 분석 단계를 추가합니다.**
   이는 `collect()` 연산자, 명령줄 구성을 위한 Groovy 클로저, 그리고 다중 명령 프로세스를 소개합니다.

!!! note

     올바른 작업 디렉토리에 있는지 확인하세요:
     `cd /workspaces/training/nf4-science/genomics`

---

## 1. 샘플별 변이 호출 단계를 GVCF를 생성하도록 수정

Part 2의 파이프라인은 VCF 파일을 생성하지만, 합동 호출에는 GVCF 파일이 필요합니다.
GVCF 변이 호출 모드를 활성화하고 출력 파일 확장자를 업데이트해야 합니다.

[Part 1](01_method.md)에서의 GVCF 변이 호출 명령을 상기해보세요:

```bash
gatk HaplotypeCaller \
        -R /data/ref/ref.fasta \
        -I /data/bam/reads_mother.bam \
        -O reads_mother.g.vcf \
        -L /data/ref/intervals.bed \
        -ERC GVCF
```

Part 2에서 래핑한 기본 HaplotypeCaller 명령과 비교하면, 차이점은 `-ERC GVCF` 매개변수와 `.g.vcf` 출력 확장자입니다.

### 1.1. HaplotypeCaller가 GVCF를 내보내도록 지시하고 출력 확장자 업데이트

`modules/gatk_haplotypecaller.nf` 모듈 파일을 열어 두 가지를 변경하세요:

- GATK HaplotypeCaller 명령에 `-ERC GVCF` 매개변수를 추가합니다;
- GATK 규칙에 따라 출력 파일 경로가 해당 `.g.vcf` 확장자를 사용하도록 업데이트합니다.

`-ERC GVCF`를 추가할 때 이전 줄 끝에 백슬래시(`\`)를 추가해야 합니다.

=== "후"

    ```groovy title="modules/gatk_haplotypecaller.nf" linenums="22" hl_lines="5 7"
        """
        gatk HaplotypeCaller \
            -R ${ref_fasta} \
            -I ${input_bam} \
            -O ${input_bam}.g.vcf \
            -L ${interval_list} \
            -ERC GVCF
        """
    ```

=== "전"

    ```groovy title="modules/gatk_haplotypecaller.nf" linenums="22" hl_lines="5"
        """
        gatk HaplotypeCaller \
            -R ${ref_fasta} \
            -I ${input_bam} \
            -O ${input_bam}.vcf \
            -L ${interval_list}
        """
    ```

새 파일 확장자와 일치하도록 output 블록도 업데이트해야 합니다.
명령 출력을 `.vcf`에서 `.g.vcf`로 변경했으므로, 프로세스 `output:` 블록도 동일한 변경 사항을 반영해야 합니다.

### 1.2. 프로세스 출력 블록에서 출력 파일 확장자 업데이트

=== "후"

    ```groovy title="modules/gatk_haplotypecaller.nf" linenums="17" hl_lines="2 3"
        output:
        path "${input_bam}.g.vcf"     , emit: vcf
        path "${input_bam}.g.vcf.idx" , emit: idx
    ```

=== "전"

    ```groovy title="modules/gatk_haplotypecaller.nf" linenums="17" hl_lines="2 3"
        output:
        path "${input_bam}.vcf"     , emit: vcf
        path "${input_bam}.vcf.idx" , emit: idx
    ```

새로운 GVCF 출력을 반영하도록 워크플로우의 publish 및 output 설정도 업데이트해야 합니다.

### 1.3. 새로운 GVCF 출력에 대한 publish 대상 업데이트

이제 VCF 대신 GVCF를 생성하므로, 더 설명적인 이름을 사용하도록 워크플로우의 `publish:` 섹션을 업데이트해야 합니다.
또한 명확성을 위해 GVCF 파일을 자체 하위 디렉토리로 구성할 것입니다.

=== "후"

    ```groovy title="genomics.nf" linenums="47" hl_lines="3 4"
        publish:
        indexed_bam = SAMTOOLS_INDEX.out
        gvcf = GATK_HAPLOTYPECALLER.out.vcf
        gvcf_idx = GATK_HAPLOTYPECALLER.out.idx
    ```

=== "전"

    ```groovy title="genomics.nf" linenums="47"
        publish:
        indexed_bam = SAMTOOLS_INDEX.out
        vcf = GATK_HAPLOTYPECALLER.out.vcf
        vcf_idx = GATK_HAPLOTYPECALLER.out.idx
    ```

이제 일치하도록 output 블록을 업데이트하세요.

### 1.4. 새로운 디렉토리 구조에 맞게 output 블록 업데이트

GVCF 파일을 `gvcf` 하위 디렉토리에 배치하도록 `output` 블록도 업데이트해야 합니다.

=== "후"

    ```groovy title="genomics.nf" linenums="53" hl_lines="3 5 6 8 9"
    output {
        indexed_bam {
            path 'indexed_bam'
        }
        gvcf {
            path 'gvcf'
        }
        gvcf_idx {
            path 'gvcf'
        }
    }
    ```

=== "전"

    ```groovy title="genomics.nf" linenums="53"
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

모듈, publish 대상, 그리고 output 블록이 모두 업데이트되었으므로, 변경 사항을 테스트할 수 있습니다.

### 1.5. 파이프라인 실행

변경 사항이 작동하는지 확인하기 위해 워크플로우를 실행하세요.

```bash
nextflow run genomics.nf
```

??? success "명령 출력"

    ```console
    N E X T F L O W   ~  version 25.10.2

    ┃ Launching `genomics.nf` [nostalgic_franklin] DSL2 - revision: f2c0a93c6a

    executor >  local (6)
    [cc/fbc705] SAMTOOLS_INDEX (3)       | 3 of 3 ✔
    [27/0d7eb9] GATK_HAPLOTYPECALLER (2) | 3 of 3 ✔
    ```

Nextflow 출력은 이전과 동일해 보이지만, `.g.vcf` 파일과 해당 인덱스 파일은 이제 하위 디렉토리로 구성되어 있습니다.

??? abstract "디렉토리 내용 (심볼릭 링크 축약됨)"

    ```console
    results/
    ├── gvcf/
    │   ├── reads_father.bam.g.vcf -> */27/0d7eb9*/reads_father.bam.g.vcf
    │   ├── reads_father.bam.g.vcf.idx -> */27/0d7eb9*/reads_father.bam.g.vcf.idx
    │   ├── reads_mother.bam.g.vcf -> */e4/4ed55e*/reads_mother.bam.g.vcf
    │   ├── reads_mother.bam.g.vcf.idx -> */e4/4ed55e*/reads_mother.bam.g.vcf.idx
    │   ├── reads_son.bam.g.vcf -> */08/e95962*/reads_son.bam.g.vcf
    │   └── reads_son.bam.g.vcf.idx -> */08/e95962*/reads_son.bam.g.vcf.idx
    └── indexed_bam/
        ├── reads_father.bam -> */9a/c7a873*/reads_father.bam
        ├── reads_father.bam.bai -> */9a/c7a873*/reads_father.bam.bai
        ├── reads_mother.bam -> */f1/8d8486*/reads_mother.bam
        ├── reads_mother.bam.bai -> */f1/8d8486*/reads_mother.bam.bai
        ├── reads_son.bam -> */cc/fbc705*/reads_son.bam
        └── reads_son.bam.bai -> */cc/fbc705*/reads_son.bam.bai
    ```

GVCF 파일 중 하나를 열어 스크롤하면, GATK HaplotypeCaller가 요청한 대로 GVCF 파일을 생성했음을 확인할 수 있습니다.

### 핵심 정리

도구 명령의 출력 파일 이름을 변경할 때는, 프로세스 `output:` 블록과 publish/output 설정이 일치하도록 업데이트되어야 합니다.

### 다음 단계

채널의 내용을 수집하고 단일 입력으로 다음 프로세스에 전달하는 방법을 학습합니다.

---

## 2. 합동 유전자형 분석 단계 추가

이제 샘플별 GVCF를 수집하고, GenomicsDB 데이터 저장소로 결합한 후, 합동 유전자형 분석을 실행하여 코호트 수준 VCF를 생성해야 합니다.
[Part 1](01_method.md)에서 다룬 것처럼, 이는 두 도구 작업입니다: GenomicsDBImport가 GVCF를 결합하고, 그런 다음 GenotypeGVCFs가 최종 변이 호출을 생성합니다.
`GATK_JOINTGENOTYPING`이라는 단일 프로세스에 두 도구를 모두 래핑하겠습니다.

[Part 1](01_method.md)의 두 명령을 상기해보세요:

```bash
gatk GenomicsDBImport \
    -V reads_mother.g.vcf \
    -V reads_father.g.vcf \
    -V reads_son.g.vcf \
    -L /data/ref/intervals.bed \
    --genomicsdb-workspace-path family_trio_gdb
```

```bash
gatk GenotypeGVCFs \
    -R /data/ref/ref.fasta \
    -V gendb://family_trio_gdb \
    -O family_trio.vcf
```

첫 번째 명령은 샘플별 GVCF와 intervals 파일을 받아 GenomicsDB 데이터 저장소를 생성합니다.
두 번째 명령은 해당 데이터 저장소와 참조 게놈을 받아 최종 코호트 수준 VCF를 생성합니다.
컨테이너 URI는 HaplotypeCaller와 동일합니다: `community.wave.seqera.io/library/gatk4:4.5.0.0--730ee8817e436867`.

### 2.1. 입력 설정

합동 유전자형 분석 프로세스에는 아직 없는 두 가지 종류의 입력이 필요합니다: 임의의 코호트 이름과 모든 샘플의 수집된 GVCF 출력이 함께 묶인 것입니다.

#### 2.1.1. `cohort_name` 매개변수 추가

코호트에 대한 임의의 이름을 제공해야 합니다.
나중에 교육 시리즈에서 이러한 종류의 작업에 샘플 메타데이터를 사용하는 방법을 배우겠지만, 지금은 `params`를 사용하여 CLI 매개변수를 선언하고 편의를 위해 기본값을 제공합니다.

=== "후"

    ```groovy title="genomics.nf" linenums="14" hl_lines="3-4"
        intervals: Path = "${projectDir}/data/ref/intervals.bed"

        // Base name for final output file
        cohort_name: String = "family_trio"
    }
    ```

=== "전"

    ```groovy title="genomics.nf" linenums="14"
        intervals: Path = "${projectDir}/data/ref/intervals.bed"
    }
    ```

#### 2.1.2. 샘플 간 HaplotypeCaller 출력 수집

`GATK_HAPLOTYPECALLER`의 출력 채널을 새 프로세스에 직접 연결하면, Nextflow는 각 샘플 GVCF에 대해 프로세스를 별도로 호출합니다.
세 개의 GVCF (및 해당 인덱스 파일)를 모두 묶어 Nextflow가 단일 프로세스 호출에 모두 함께 전달하도록 하려고 합니다.

`collect()` 채널 연산자를 사용하여 이를 수행할 수 있습니다.
GATK_HAPLOTYPECALLER 호출 바로 다음에 `workflow` 본문에 다음 줄을 추가하세요:

=== "후"

    ```groovy title="genomics.nf" hl_lines="4-6"
            intervals_file
        )

        // Collect variant calling outputs across samples
        all_gvcfs_ch = GATK_HAPLOTYPECALLER.out.vcf.collect()
        all_idxs_ch = GATK_HAPLOTYPECALLER.out.idx.collect()
    ```

=== "전"

    ```groovy title="genomics.nf"
            intervals_file
        )
    ```

이를 분석하면:

1. `.out` 속성을 사용하여 `GATK_HAPLOTYPECALLER`의 출력 채널을 가져옵니다.
2. 섹션 1에서 `emit:`을 사용하여 출력 이름을 지정했기 때문에, `.vcf`로 GVCF를, `.idx`로 인덱스 파일을 선택할 수 있습니다. 명명된 출력이 없다면 `.out[0]`와 `.out[1]`을 사용해야 합니다.
3. `collect()` 연산자는 모든 파일을 단일 요소로 묶으므로, `all_gvcfs_ch`에는 세 개의 GVCF가 모두 함께 포함되고, `all_idxs_ch`에는 세 개의 인덱스 파일이 모두 함께 포함됩니다.

Nextflow가 실행을 위해 모든 입력 파일을 함께 스테이징하므로 인덱스 파일이 GVCF와 함께 존재할 것이기 때문에, GVCF와 인덱스 파일을 (튜플로 함께 유지하는 대신) 별도로 수집할 수 있습니다.

!!! tip

    채널 연산자를 적용하기 전과 후에 `view()` 연산자를 사용하여 채널의 내용을 검사할 수 있습니다.

### 2.2. 합동 유전자형 분석 프로세스를 작성하고 워크플로우에서 호출

Part 2에서 사용한 것과 동일한 패턴을 따라, 모듈 파일에 프로세스 정의를 작성하고, 워크플로우로 가져온 후, 방금 준비한 입력으로 호출하겠습니다.

#### 2.2.1. 각 GVCF에 `-V` 인수를 제공하는 문자열 구성

프로세스 정의를 채우기 시작하기 전에, 해결해야 할 한 가지가 있습니다.
GenomicsDBImport 명령은 다음과 같이 각 GVCF 파일에 대해 별도의 `-V` 인수를 예상합니다:

```bash
gatk GenomicsDBImport \
    -V reads_mother.bam.g.vcf \
    -V reads_father.bam.g.vcf \
    -V reads_son.bam.g.vcf \
    ...
```

`-V ${all_gvcfs_ch}`로 작성하면, Nextflow는 단순히 파일 이름을 연결하여 명령의 해당 부분이 다음과 같이 보일 것입니다:

```groovy
-V reads_mother.bam.g.vcf reads_father.bam.g.vcf reads_son.bam.g.vcf
```

하지만 문자열이 다음과 같이 보여야 합니다:

```groovy
-V reads_mother.bam.g.vcf -V reads_father.bam.g.vcf -V reads_son.bam.g.vcf
```

중요한 것은, 수집된 채널에 있는 파일이 무엇이든 간에 이 문자열을 동적으로 구성해야 한다는 것입니다.
Nextflow는 (Groovy를 통해) 이를 수행하는 간결한 방법을 제공합니다:

```groovy
def gvcfs_line = all_gvcfs.collect { gvcf -> "-V ${gvcf}" }.join(' ')
```

이를 분석하면:

1. `all_gvcfs.collect { gvcf -> "-V ${gvcf}" }`는 각 파일 경로를 반복하고 `-V `를 앞에 추가하여 `["-V A.g.vcf", "-V B.g.vcf", "-V C.g.vcf"]`를 생성합니다.
2. `.join(' ')`는 공백으로 연결합니다: `"-V A.g.vcf -V B.g.vcf -V C.g.vcf"`.
3. 결과는 로컬 변수 `gvcfs_line` (`def`로 정의됨)에 할당되며, 명령 템플릿에 삽입할 수 있습니다.

이 줄은 프로세스의 `script:` 블록 내부, 명령 템플릿 전에 들어갑니다.
`script:`와 명령 템플릿의 여는 `"""` 사이에 임의의 Groovy 코드를 배치할 수 있습니다.

그러면 프로세스의 `script:` 블록에서 전체 문자열을 `gvcfs_line`으로 참조할 수 있습니다.

#### 2.2.2. 합동 유전자형 분석 프로세스의 모듈 작성

이제 전체 프로세스 작성에 착수할 수 있습니다.

`modules/gatk_jointgenotyping.nf`를 열고 프로세스 정의의 개요를 살펴보세요.

위에서 제공된 정보를 사용하여 프로세스 정의를 작성한 다음, 아래 "후" 탭의 해결책과 대조하여 작업을 확인하세요.

=== "전"

    ```groovy title="modules/gatk_jointgenotyping.nf" linenums="1"
    #!/usr/bin/env nextflow

    /*
     * GVCF를 GenomicsDB 데이터 저장소로 결합하고 합동 유전자형 분석을 실행하여 코호트 수준 호출 생성
     */
    process GATK_JOINTGENOTYPING {

        container

        input:

        output:

        script:
        """

        """
    }
    ```

=== "후"

    ```groovy title="modules/gatk_jointgenotyping.nf" linenums="1" hl_lines="8 11-17 20-21 24-25 29-33"
    #!/usr/bin/env nextflow

    /*
     * GVCF를 GenomicsDB 데이터 저장소로 결합하고 합동 유전자형 분석을 실행하여 코호트 수준 호출 생성
     */
    process GATK_JOINTGENOTYPING {

        container "community.wave.seqera.io/library/gatk4:4.5.0.0--730ee8817e436867"

        input:
        path all_gvcfs
        path all_idxs
        path interval_list
        val cohort_name
        path ref_fasta
        path ref_index
        path ref_dict

        output:
        path "${cohort_name}.joint.vcf"     , emit: vcf
        path "${cohort_name}.joint.vcf.idx" , emit: idx

        script:
        def gvcfs_line = all_gvcfs.collect { gvcf -> "-V ${gvcf}" }.join(' ')
        """
        gatk GenomicsDBImport \
            ${gvcfs_line} \
            -L ${interval_list} \
            --genomicsdb-workspace-path ${cohort_name}_gdb

        gatk GenotypeGVCFs \
            -R ${ref_fasta} \
            -V gendb://${cohort_name}_gdb \
            -L ${interval_list} \
            -O ${cohort_name}.joint.vcf
        """
    }
    ```

여기서 주목할 만한 몇 가지가 있습니다.

이전과 마찬가지로, 명령이 직접 참조하지 않더라도 여러 입력이 나열됩니다: `all_idxs`, `ref_index`, 그리고 `ref_dict`.
이들을 나열하면 Nextflow가 명령에 나타나는 파일과 함께 작업 디렉토리에 이러한 파일을 스테이징하며, GATK는 명명 규칙에 따라 이를 찾을 것으로 예상합니다.

`gvcfs_line` 변수는 GenomicsDBImport에 대한 `-V` 인수를 구성하기 위해 위에서 설명한 Groovy 클로저를 사용합니다.

이 프로세스는 터미널에서 수행하는 것처럼 두 개의 명령을 순차적으로 실행합니다.
GenomicsDBImport는 샘플별 GVCF를 데이터 저장소로 결합하고, 그런 다음 GenotypeGVCFs가 해당 데이터 저장소를 읽고 최종 코호트 수준 VCF를 생성합니다.
GenomicsDB 데이터 저장소 (`${cohort_name}_gdb`)는 프로세스 내에서만 사용되는 중간 산물입니다; output 블록에는 나타나지 않습니다.

이를 완료하면, 프로세스를 사용할 준비가 된 것입니다.
워크플로우에서 사용하려면, 모듈을 가져오고 프로세스 호출을 추가해야 합니다.

#### 2.2.3. 모듈 가져오기

기존 import 문 아래에 `genomics.nf`에 import 문을 추가하세요:

=== "후"

    ```groovy title="genomics.nf" linenums="21" hl_lines="3"
    include { SAMTOOLS_INDEX } from './modules/samtools_index.nf'
    include { GATK_HAPLOTYPECALLER } from './modules/gatk_haplotypecaller.nf'
    include { GATK_JOINTGENOTYPING } from './modules/gatk_jointgenotyping.nf'
    ```

=== "전"

    ```groovy title="genomics.nf" linenums="21"
    include { SAMTOOLS_INDEX } from './modules/samtools_index.nf'
    include { GATK_HAPLOTYPECALLER } from './modules/gatk_haplotypecaller.nf'
    ```

프로세스가 이제 워크플로우 범위에서 사용 가능합니다.

#### 2.2.4. 프로세스 호출 추가

`collect()` 줄 다음에 워크플로우 본문에 `GATK_JOINTGENOTYPING` 호출을 추가하세요:

=== "후"

    ```groovy title="genomics.nf" hl_lines="3-12"
        all_idxs_ch = GATK_HAPLOTYPECALLER.out.idx.collect()

        // Combine GVCFs into a GenomicsDB data store and apply joint genotyping
        GATK_JOINTGENOTYPING(
            all_gvcfs_ch,
            all_idxs_ch,
            intervals_file,
            params.cohort_name,
            ref_file,
            ref_index_file,
            ref_dict_file
        )
    ```

=== "전"

    ```groovy title="genomics.nf"
        all_idxs_ch = GATK_HAPLOTYPECALLER.out.idx.collect()
    ```

프로세스가 이제 완전히 연결되었습니다.
다음으로, 출력이 게시되는 방식을 설정합니다.

### 2.3. 출력 처리 설정

합동 VCF 출력을 게시해야 합니다.
합동 유전자형 분석 결과에 대한 publish 대상과 output 블록 항목을 추가하세요.

#### 2.3.1. 합동 VCF에 대한 publish 대상 추가

워크플로우의 `publish:` 섹션에 합동 VCF와 해당 인덱스를 추가하세요:

=== "후"

    ```groovy title="genomics.nf" hl_lines="5-6"
        publish:
        indexed_bam = SAMTOOLS_INDEX.out
        gvcf = GATK_HAPLOTYPECALLER.out.vcf
        gvcf_idx = GATK_HAPLOTYPECALLER.out.idx
        joint_vcf = GATK_JOINTGENOTYPING.out.vcf
        joint_vcf_idx = GATK_JOINTGENOTYPING.out.idx
    ```

=== "전"

    ```groovy title="genomics.nf"
        publish:
        indexed_bam = SAMTOOLS_INDEX.out
        gvcf = GATK_HAPLOTYPECALLER.out.vcf
        gvcf_idx = GATK_HAPLOTYPECALLER.out.idx
    ```

이제 일치하도록 output 블록을 업데이트하세요.

#### 2.3.2. 합동 VCF에 대한 output 블록 항목 추가

합동 VCF 파일에 대한 항목을 추가하세요.
최종 출력이므로 results 디렉토리의 루트에 배치할 것입니다.

=== "후"

    ```groovy title="genomics.nf" hl_lines="11-16"
    output {
        indexed_bam {
            path 'indexed_bam'
        }
        gvcf {
            path 'gvcf'
        }
        gvcf_idx {
            path 'gvcf'
        }
        joint_vcf {
            path '.'
        }
        joint_vcf_idx {
            path '.'
        }
    }
    ```

=== "전"

    ```groovy title="genomics.nf"
    output {
        indexed_bam {
            path 'indexed_bam'
        }
        gvcf {
            path 'gvcf'
        }
        gvcf_idx {
            path 'gvcf'
        }
    }
    ```

프로세스, publish 대상, 그리고 output 블록이 모두 준비되었으므로, 전체 워크플로우를 테스트할 수 있습니다.

### 2.4. 워크플로우 실행

모든 것이 작동하는지 확인하기 위해 워크플로우를 실행하세요.

```bash
nextflow run genomics.nf -resume
```

??? success "명령 출력"

    ```console
    N E X T F L O W   ~  version 25.10.2

    ┃ Launching `genomics.nf` [crazy_marconi] DSL2 - revision: 5da9afc841

    executor >  local (1)
    [9a/c7a873] SAMTOOLS_INDEX (2)       | 3 of 3, cached: 3 ✔
    [e4/4ed55e] GATK_HAPLOTYPECALLER (2) | 3 of 3, cached: 3 ✔
    [a6/7cc8ed] GATK_JOINTGENOTYPING     | 1 of 1 ✔
    ```

처음 두 단계는 이전 실행에서 캐시되고, 새로운 `GATK_JOINTGENOTYPING` 단계는 세 샘플 모두의 수집된 입력에 대해 한 번 실행됩니다.
최종 출력 파일인 `family_trio.joint.vcf` (및 해당 인덱스)는 results 디렉토리에 있습니다.

??? abstract "디렉토리 내용 (심볼릭 링크 축약됨)"

    ```console
    results/
    ├── family_trio.joint.vcf -> */a6/7cc8ed*/family_trio.joint.vcf
    ├── family_trio.joint.vcf.idx -> */a6/7cc8ed*/family_trio.joint.vcf.idx
    ├── gvcf/
    │   ├── reads_father.bam.g.vcf -> */27/0d7eb9*/reads_father.bam.g.vcf
    │   ├── reads_father.bam.g.vcf.idx -> */27/0d7eb9*/reads_father.bam.g.vcf.idx
    │   ├── reads_mother.bam.g.vcf -> */e4/4ed55e*/reads_mother.bam.g.vcf
    │   ├── reads_mother.bam.g.vcf.idx -> */e4/4ed55e*/reads_mother.bam.g.vcf.idx
    │   ├── reads_son.bam.g.vcf -> */08/e95962*/reads_son.bam.g.vcf
    │   └── reads_son.bam.g.vcf.idx -> */08/e95962*/reads_son.bam.g.vcf.idx
    └── indexed_bam/
        ├── reads_father.bam -> */9a/c7a873*/reads_father.bam
        ├── reads_father.bam.bai -> */9a/c7a873*/reads_father.bam.bai
        ├── reads_mother.bam -> */f1/8d8486*/reads_mother.bam
        ├── reads_mother.bam.bai -> */f1/8d8486*/reads_mother.bam.bai
        ├── reads_son.bam -> */cc/fbc705*/reads_son.bam
        └── reads_son.bam.bai -> */cc/fbc705*/reads_son.bam.bai
    ```

합동 VCF 파일을 열면, 워크플로우가 예상된 변이 호출을 생성했음을 확인할 수 있습니다.

```console title="family_trio.joint.vcf" linenums="40"
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	reads_father	reads_mother	reads_son
20_10037292_10066351	3480	.	C	CT	1625.89	.	AC=5;AF=0.833;AN=6;BaseQRankSum=0.220;DP=85;ExcessHet=0.0000;FS=2.476;MLEAC=5;MLEAF=0.833;MQ=60.00;MQRankSum=0.00;QD=21.68;ReadPosRankSum=-1.147e+00;SOR=0.487	GT:AD:DP:GQ:PL	0/1:15,16:31:99:367,0,375	1/1:0,18:18:54:517,54,0	1/1:0,26:26:78:756,78,0
20_10037292_10066351	3520	.	AT	A	1678.89	.	AC=5;AF=0.833;AN=6;BaseQRankSum=1.03;DP=80;ExcessHet=0.0000;FS=2.290;MLEAC=5;MLEAF=0.833;MQ=60.00;MQRankSum=0.00;QD=22.39;ReadPosRankSum=0.701;SOR=0.730	GT:AD:DP:GQ:PL	0/1:18,13:31:99:296,0,424	1/1:0,18:18:54:623,54,0	1/1:0,26:26:78:774,78,0
20_10037292_10066351	3529	.	T	A	154.29	.	AC=1;AF=0.167;AN=6;BaseQRankSum=-5.440e-01;DP=104;ExcessHet=0.0000;FS=1.871;MLEAC=1;MLEAF=0.167;MQ=60.00;MQRankSum=0.00;QD=7.71;ReadPosRankSum=-1.158e+00;SOR=1.034	GT:AD:DP:GQ:PL	0/0:44,0:44:99:0,112,1347	0/1:12,8:20:99:163,0,328	0/0:39,0:39:99:0,105,1194
```

이제 자동화되고 완전히 재현 가능한 합동 변이 호출 워크플로우를 갖게 되었습니다!

!!! note

    제공된 데이터 파일은 20번 염색체의 아주 작은 부분만 다룬다는 점을 기억하세요.
    실제 변이 호출 세트의 크기는 수백만 개의 변이로 계산됩니다.
    그래서 교육 목적으로 아주 작은 데이터 하위 집합만 사용합니다!

### 핵심 정리

채널에서 출력을 수집하고 다른 프로세스에 단일 입력으로 묶는 방법을 알게 되었습니다.
또한 Groovy 클로저를 사용하여 명령줄을 구성하는 방법과 단일 프로세스에서 여러 명령을 실행하는 방법도 알게 되었습니다.

### 다음 단계

스스로를 크게 칭찬하세요! Nextflow for Genomics 과정을 완료했습니다.

학습한 내용을 복습하고 다음 단계를 알아보려면 최종 [과정 요약](./next_steps.md)으로 이동하세요.
