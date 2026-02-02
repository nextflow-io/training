# 파트 2: 코호트에 대한 조인트 콜링

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } AI 지원 번역 - [자세히 알아보기 및 개선 제안](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

이 과정의 첫 번째 파트에서는 완전히 선형적이고 각 샘플의 데이터를 독립적으로 처리하는 변이 콜링 파이프라인을 구축했습니다.
그러나 실제 유전체학 사용 사례에서는 일반적으로 여러 샘플의 변이 콜을 함께 살펴봐야 합니다.

두 번째 파트에서는 채널과 채널 연산자를 사용하여 GATK로 조인트 변이 콜링을 구현하는 방법을 보여드리며, 파트 1의 파이프라인을 기반으로 합니다.

### 방법 개요

이 과정의 첫 번째 파트에서 사용한 GATK 변이 콜링 방법은 단순히 샘플별 변이 콜을 생성했습니다.
각 샘플의 변이를 개별적으로만 보려는 경우에는 괜찮지만, 제한적인 정보만 제공합니다.
여러 샘플에서 변이 콜이 어떻게 다른지 살펴보는 것이 더 흥미로운 경우가 많으며, 이를 위해 GATK는 조인트 변이 콜링이라는 대체 방법을 제공합니다.

조인트 변이 콜링은 각 샘플에 대해 GVCF(Genomic VCF)라는 특별한 종류의 변이 출력을 생성한 다음, 모든 샘플의 GVCF 데이터를 결합하고 마지막으로 '조인트 유전자형 분석(joint genotyping)' 통계 분석을 실행하는 과정을 포함합니다.

![조인트 분석](img/joint-calling.png)

샘플의 GVCF가 특별한 이유는 프로그램이 변이의 증거를 발견한 위치뿐만 아니라 유전체의 대상 영역에 있는 모든 위치에 대한 시퀀스 데이터 통계를 요약하는 레코드를 포함하기 때문입니다.
이는 조인트 유전자형 분석 계산에 매우 중요합니다([추가 자료](https://gatk.broadinstitute.org/hc/en-us/articles/360035890431-The-logic-of-joint-calling-for-germline-short-variants)).

GVCF는 파트 1에서 사용한 것과 동일한 도구인 GATK HaplotypeCaller에 추가 매개변수(`-ERC GVCF`)를 사용하여 생성됩니다.
GVCF를 결합하는 작업은 GATK GenomicsDBImport로 수행되며, 이는 샘플별 콜을 데이터 저장소(데이터베이스와 유사)로 결합한 다음, 실제 '조인트 유전자형 분석' 분석은 GATK GenotypeGVCFs로 수행됩니다.

### 워크플로

요약하자면, 이 과정의 이 파트에서는 다음을 수행하는 워크플로를 개발할 것입니다:

<figure class="excalidraw">
--8<-- "docs/en/docs/nf4_science/genomics/img/hello-gatk-2.svg"
</figure>

1. Samtools를 사용하여 각 BAM 입력 파일에 대한 인덱스 파일 생성
2. 각 BAM 입력 파일에 대해 GATK HaplotypeCaller를 실행하여 샘플별 유전체 변이 콜의 GVCF 생성
3. 모든 GVCF를 수집하고 GenomicsDB 데이터 저장소로 결합
4. 결합된 GVCF 데이터 저장소에 대해 조인트 유전자형 분석을 실행하여 코호트 수준의 VCF 생성

이를 파트 1과 동일한 데이터셋에 적용할 것입니다.

---

## 0. 준비 운동: Samtools 및 GATK를 직접 실행

이전과 마찬가지로, 워크플로로 적용하기 전에 명령을 수동으로 시도해보고자 합니다.

!!! note

     올바른 작업 디렉토리에 있는지 확인하세요:
     `cd /workspaces/training/nf4-science/genomics`

### 0.1. Samtools로 BAM 입력 파일 인덱싱

이 첫 번째 단계는 파트 1과 동일하므로 매우 익숙하게 느껴질 것입니다. 하지만 이번에는 세 개의 샘플 모두에 대해 수행해야 합니다.

!!! note

    파이프라인을 통해 이미 세 샘플에 대한 인덱스 파일을 생성했으므로, 결과 디렉토리에서 찾을 수 있습니다. 그러나 수동으로 다시 수행하는 것이 더 깔끔하며, 1분밖에 걸리지 않습니다.

#### 0.1.1. Samtools 컨테이너를 대화형으로 실행

```bash
docker run -it -v ./data:/data community.wave.seqera.io/library/samtools:1.20--b5dfbd93de237464
```

<!--
??? success "명령 출력"

    ```console

    ```
-->

#### 0.1.2. 세 샘플에 대해 인덱싱 명령 실행

```bash
samtools index /data/bam/reads_mother.bam
samtools index /data/bam/reads_father.bam
samtools index /data/bam/reads_son.bam
```

이전과 마찬가지로 해당 BAM 파일과 동일한 디렉토리에 인덱스 파일이 생성됩니다.

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

이제 세 샘플 모두에 대한 인덱스 파일이 있으므로 각 샘플에 대한 GVCF 생성을 진행할 수 있습니다.

#### 0.1.3. Samtools 컨테이너 종료

```bash
exit
```

### 0.2. GATK HaplotypeCaller로 GVCF 모드에서 변이 콜링

이 두 번째 단계는 파트 1: Hello Genomics에서 수행한 것과 매우 유사하지만, 이제 'GVCF 모드'에서 GATK를 실행할 것입니다.

#### 0.2.1. GATK 컨테이너를 대화형으로 실행

```bash
docker run -it -v ./data:/data community.wave.seqera.io/library/gatk4:4.5.0.0--730ee8817e436867
```

<!--
??? success "명령 출력"

    ```console

    ```
-->

#### 0.2.2. GVCF 옵션으로 변이 콜링 명령 실행

genomic VCF(GVCF)를 생성하기 위해 기본 명령에 `-ERC GVCF` 옵션을 추가하여 HaplotypeCaller의 GVCF 모드를 활성화합니다.

또한 출력 파일의 파일 확장자를 `.vcf`에서 `.g.vcf`로 변경합니다.
이는 기술적으로 요구 사항은 아니지만 강력히 권장되는 규칙입니다.

```bash
gatk HaplotypeCaller \
        -R /data/ref/ref.fasta \
        -I /data/bam/reads_mother.bam \
        -O reads_mother.g.vcf \
        -L /data/ref/intervals.bed \
        -ERC GVCF
```

<!--
??? success "명령 출력"

    ```console

    ```
-->

이렇게 하면 컨테이너의 현재 작업 디렉토리에 GVCF 출력 파일 `reads_mother.g.vcf`가 생성됩니다.

`cat`으로 내용을 보면 파트 1에서 생성한 동등한 VCF보다 훨씬 길다는 것을 알 수 있습니다. 파일의 시작 부분까지 스크롤할 수도 없고, 대부분의 줄이 파트 1의 VCF에서 본 것과 상당히 다릅니다.

```console title="출력" linenums="1674"
20_10037292_10066351    14714   .       T       <NON_REF>       .       .       END=14718       GT:DP:GQ:MIN_DP:PL       0/0:37:99:37:0,99,1192
20_10037292_10066351    14719   .       T       <NON_REF>       .       .       END=14719       GT:DP:GQ:MIN_DP:PL       0/0:36:82:36:0,82,1087
20_10037292_10066351    14720   .       T       <NON_REF>       .       .       END=14737       GT:DP:GQ:MIN_DP:PL       0/0:42:99:37:0,100,1160
```

이들은 변이 콜러가 변이의 증거를 발견하지 못한 비변이 영역을 나타내므로, 변이가 없다는 확신 수준을 설명하는 몇 가지 통계를 캡처했습니다. 이를 통해 두 가지 매우 다른 경우를 구별할 수 있습니다: (1) 샘플이 동형접합 참조임을 보여주는 양질의 데이터가 있는 경우, (2) 어느 쪽으로든 결정을 내리기에 충분한 좋은 데이터가 없는 경우.

GVCF에는 일반적으로 이러한 비변이 줄이 많이 있으며, 그 사이에 소수의 변이 레코드가 흩어져 있습니다. GVCF에서 `head -176`을 실행하여 파일의 처음 176줄만 로드하여 실제 변이 콜을 찾아보세요.

```console title="출력" linenums="174"
20_10037292_10066351    3479    .       T       <NON_REF>       .       .       END=3479        GT:DP:GQ:MIN_DP:PL       0/0:34:36:34:0,36,906
20_10037292_10066351    3480    .       C       CT,<NON_REF>    503.03  .       DP=23;ExcessHet=0.0000;MLEAC=2,0;MLEAF=1.00,0.00;RAW_MQandDP=82800,23    GT:AD:DP:GQ:PL:SB       1/1:0,18,0:18:54:517,54,0,517,54,517:0,0,7,11
20_10037292_10066351    3481    .       T       <NON_REF>       .       .       END=3481        GT:DP:GQ:MIN_DP:PL       0/0:21:51:21:0,51,765
```

두 번째 줄은 파일의 첫 번째 변이 레코드를 보여주며, 이는 파트 1에서 살펴본 VCF 파일의 첫 번째 변이에 해당합니다.

원래 VCF와 마찬가지로 출력 GVCF 파일에도 `reads_mother.g.vcf.idx`라는 인덱스 파일이 함께 제공됩니다.

#### 0.2.3. 다른 두 샘플에 대해 프로세스 반복

조인트 유전자형 분석 단계를 테스트하려면 세 샘플 모두에 대한 GVCF가 필요하므로 지금 수동으로 생성하겠습니다.

```bash
gatk HaplotypeCaller \
        -R /data/ref/ref.fasta \
        -I /data/bam/reads_father.bam \
        -O reads_father.g.vcf \
        -L /data/ref/intervals.bed \
        -ERC GVCF
```

<!--
??? success "명령 출력"

    ```console

    ```
-->

```bash
gatk HaplotypeCaller \
        -R /data/ref/ref.fasta \
        -I /data/bam/reads_son.bam \
        -O reads_son.g.vcf \
        -L /data/ref/intervals.bed \
        -ERC GVCF
```

<!--
??? success "명령 출력"

    ```console

    ```
-->

이 작업이 완료되면 현재 디렉토리에 `.g.vcf`로 끝나는 세 개의 파일(샘플당 하나)과 해당 인덱스 파일(`.g.vcf.idx`로 끝남)이 있어야 합니다.

### 0.3. 조인트 유전자형 분석 실행

이제 모든 GVCF가 있으므로 마침내 샘플 코호트에 대한 변이 콜을 생성하는 조인트 유전자형 분석 방법을 시도할 수 있습니다.
복습하자면, 모든 GVCF의 데이터를 데이터 저장소로 결합한 다음, 조인트 유전자형 분석 자체를 실행하여 조인트 콜 변이의 최종 VCF를 생성하는 2단계 방법입니다.

#### 0.3.1. 모든 샘플별 GVCF 결합

이 첫 번째 단계는 GATK의 또 다른 도구인 GenomicsDBImport를 사용하여 모든 GVCF의 데이터를 GenomicsDB 데이터 저장소로 결합합니다.

```bash
gatk GenomicsDBImport \
    -V reads_mother.g.vcf \
    -V reads_father.g.vcf \
    -V reads_son.g.vcf \
    -L /data/ref/intervals.bed \
    --genomicsdb-workspace-path family_trio_gdb
```

<!--
??? success "명령 출력"

    ```console

    ```
-->

이 단계의 출력은 실제로 여러 다른 파일 형태로 결합된 변이 데이터를 보유하는 추가 중첩 디렉토리 세트를 포함하는 디렉토리입니다.
둘러볼 수 있지만 이 데이터 저장소 형식은 사람이 직접 읽도록 의도되지 않았음을 빠르게 알 수 있습니다.

!!! note

    GATK에는 필요에 따라 데이터 저장소에서 변이 콜 데이터를 검사하고 추출할 수 있는 도구가 포함되어 있습니다.

#### 0.3.2. 조인트 유전자형 분석 자체 실행

이 두 번째 단계는 GATK의 또 다른 도구인 GenotypeGVCFs를 사용하여 코호트의 모든 샘플에서 사용 가능한 데이터를 고려하여 변이 통계 및 개별 유전자형을 다시 계산합니다.

```bash
gatk GenotypeGVCFs \
    -R /data/ref/ref.fasta \
    -V gendb://family_trio_gdb \
    -O family_trio.vcf
```

<!--
??? success "명령 출력"

    ```console

    ```
-->

이렇게 하면 컨테이너의 현재 작업 디렉토리에 VCF 출력 파일 `family_trio.vcf`가 생성됩니다.
또 다른 합리적으로 작은 파일이므로 이 파일을 `cat`하여 내용을 보고 위로 스크롤하여 처음 몇 개의 변이 줄을 찾을 수 있습니다.

```console title="family_trio.vcf" linenums="40"
#CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO    FORMAT  reads_father    reads_mother    reads_son
20_10037292_10066351    3480    .       C       CT      1625.89 .       AC=5;AF=0.833;AN=6;BaseQRankSum=0.220;DP=85;ExcessHet=0.0000;FS=2.476;MLEAC=5;MLEAF=0.833;MQ=60.00;MQRankSum=0.00;QD=21.68;ReadPosRankSum=-1.147e+00;SOR=0.487    GT:AD:DP:GQ:PL  0/1:15,16:31:99:367,0,375       1/1:0,18:18:54:517,54,0 1/1:0,26:26:78:756,78,0
20_10037292_10066351    3520    .       AT      A       1678.89 .       AC=5;AF=0.833;AN=6;BaseQRankSum=1.03;DP=80;ExcessHet=0.0000;FS=2.290;MLEAC=5;MLEAF=0.833;MQ=60.00;MQRankSum=0.00;QD=22.39;ReadPosRankSum=0.701;SOR=0.730 GT:AD:DP:GQ:PL   0/1:18,13:31:99:296,0,424       1/1:0,18:18:54:623,54,0 1/1:0,26:26:78:774,78,0
20_10037292_10066351    3529    .       T       A       154.29  .       AC=1;AF=0.167;AN=6;BaseQRankSum=-5.440e-01;DP=104;ExcessHet=0.0000;FS=1.871;MLEAC=1;MLEAF=0.167;MQ=60.00;MQRankSum=0.00;QD=7.71;ReadPosRankSum=-1.158e+00;SOR=1.034       GT:AD:DP:GQ:PL  0/0:44,0:44:99:0,112,1347       0/1:12,8:20:99:163,0,328        0/0:39,0:39:99:0,105,1194
```

이것은 파트 1에서 생성한 원래 VCF와 더 비슷해 보이지만, 이번에는 세 샘플 모두에 대한 유전자형 수준 정보가 있습니다.
파일의 마지막 세 열은 알파벳 순서로 나열된 샘플의 유전자형 블록입니다.

첫 번째 변이에 대해 테스트 가족 삼인조에 대해 호출된 유전자형을 살펴보면, 아버지는 이형접합 변이(`0/1`)이고 어머니와 아들은 모두 동형접합 변이(`1/1`)임을 알 수 있습니다.

이것이 궁극적으로 데이터셋에서 추출하려는 정보입니다! 이제 이 모든 것을 Nextflow 워크플로로 적용하여 대규모로 수행할 수 있습니다.

#### 0.3.3. GATK 컨테이너 종료

```bash
exit
```

### 핵심 내용

터미널에서 조인트 변이 콜링과 관련된 개별 명령을 실행하여 원하는 정보를 생성할 수 있는지 확인하는 방법을 알게 되었습니다.

### 다음 단계는?

이러한 명령을 실제 파이프라인으로 적용합니다.

---

## 1. 샘플별 변이 콜링 단계를 수정하여 GVCF 생성

좋은 소식은 파트 1에서 이미 이 작업의 일부를 수행하는 워크플로를 작성했기 때문에 처음부터 다시 시작할 필요가 없다는 것입니다.
그러나 해당 파이프라인은 VCF 파일을 생성하는 반면, 이제 조인트 유전자형 분석을 수행하기 위해 GVCF 파일을 원합니다.
따라서 GVCF 변이 콜링 모드를 활성화하고 출력 파일 확장자를 업데이트하는 것부터 시작해야 합니다.

!!! note

    편의를 위해 파트 1의 끝에 있는 GATK 워크플로의 새 복사본으로 작업하되 다른 이름인 `genomics-2.nf`를 사용합니다.

### 1.1. HaplotypeCaller에 GVCF 출력을 지시하고 출력 확장자 업데이트

코드 편집기에서 `genomics-2.nf` 파일을 열어봅시다.
매우 익숙하게 느껴질 것이지만, 예상대로 실행되는지 확인하고 싶다면 자유롭게 실행하세요.

두 가지 변경부터 시작하겠습니다:

- GATK HaplotypeCaller 명령에 `-ERC GVCF` 매개변수 추가
- GATK 규칙에 따라 해당 `.g.vcf` 확장자를 사용하도록 출력 파일 경로 업데이트

`-ERC GVCF`를 추가할 때 이전 줄 끝에 백슬래시(`\`)를 추가해야 합니다.

=== "수정 후"

    ```groovy title="genomics-2.nf" linenums="56" hl_lines="4 6"
        """
        gatk HaplotypeCaller \
            -R ${ref_fasta} \
            -I ${input_bam} \
            -O ${input_bam}.g.vcf \
            -L ${interval_list} \
            -ERC GVCF
        """
    ```

=== "수정 전"

    ```groovy title="genomics-2.nf" linenums="56" hl_lines="4"
        """
        gatk HaplotypeCaller \
            -R ${ref_fasta} \
            -I ${input_bam} \
            -O ${input_bam}.vcf \
            -L ${interval_list}
        """
    ```

그리고 이것이 HaplotypeCaller를 VCF 대신 GVCF를 생성하도록 전환하는 데 필요한 전부입니다, 맞죠?

### 1.2. 파이프라인을 실행하여 GVCF를 생성할 수 있는지 확인

Nextflow 실행 명령은 워크플로 파일 이름 자체를 제외하고 이전과 동일합니다.
적절하게 업데이트했는지 확인하세요.

```bash
nextflow run genomics-2.nf
```

??? success "명령 출력"

    ```console
    N E X T F L O W   ~  version 25.10.2

    ┃ Launching `genomics-2.nf` [crazy_venter] DSL2 - revision: a2d6f6f09f

    executor >  local (6)
    [f1/8d8486] SAMTOOLS_INDEX (1)       | 3 of 3 ✔
    [72/3249ca] GATK_HAPLOTYPECALLER (3) | 0 of 3
    ERROR ~ Error executing process > 'GATK_HAPLOTYPECALLER (2)'

    Caused by:
      Missing output file(s) `reads_son.bam.vcf` expected by process `GATK_HAPLOTYPECALLER (2)`

    Command executed:

      gatk HaplotypeCaller         -R ref.fasta         -I reads_son.bam         -O reads_son.bam.g.vcf         -L intervals.bed         -ERC GVCF
    ```

그리고 출력은... 모두 빨간색입니다! 오 이런.

실행된 명령은 올바르므로 GATK 도구의 동작을 변경하기에 충분하다는 것이 맞았습니다.
하지만 누락된 출력 파일에 대한 줄을 보세요. 뭔가 눈에 띄나요?

맞습니다, Nextflow에게 새 파일 이름을 예상하도록 알려주는 것을 잊었습니다. 죄송합니다.

### 1.3. process 출력 블록에서도 출력 파일 확장자 업데이트

도구 명령 자체의 파일 확장자만 변경하는 것으로는 충분하지 않습니다. 예상되는 출력 파일 이름이 변경되었음을 Nextflow에게도 알려야 합니다.

=== "수정 후"

    ```groovy title="genomics-2.nf" linenums="50" hl_lines="2 3"
        output:
        path "${input_bam}.g.vcf"     , emit: vcf
        path "${input_bam}.g.vcf.idx" , emit: idx
    ```

=== "수정 전"

    ```groovy title="genomics-2.nf" linenums="50" hl_lines="2 3"
        output:
        path "${input_bam}.vcf"     , emit: vcf
        path "${input_bam}.vcf.idx" , emit: idx
    ```

### 1.4. 새 GVCF 출력에 대한 게시 대상 업데이트

이제 VCF 대신 GVCF를 생성하므로 더 설명적인 이름을 사용하도록 워크플로의 `publish:` 섹션을 업데이트해야 합니다.
또한 명확성을 위해 GVCF 파일을 자체 하위 디렉토리로 구성합니다.

=== "수정 후"

    ```groovy title="genomics-2.nf" linenums="88" hl_lines="3 4"
        publish:
        indexed_bam = SAMTOOLS_INDEX.out
        gvcf = GATK_HAPLOTYPECALLER.out.vcf
        gvcf_idx = GATK_HAPLOTYPECALLER.out.idx
    ```

=== "수정 전"

    ```groovy title="genomics-2.nf" linenums="88"
        publish:
        indexed_bam = SAMTOOLS_INDEX.out
        vcf = GATK_HAPLOTYPECALLER.out.vcf
        vcf_idx = GATK_HAPLOTYPECALLER.out.idx
    ```

### 1.5. 새 디렉토리 구조에 맞게 output 블록 업데이트

또한 GVCF 파일을 `gvcf` 하위 디렉토리에 넣도록 `output` 블록을 업데이트해야 합니다.

=== "수정 후"

    ```groovy title="genomics-2.nf" linenums="94" hl_lines="3 5 6 8 9"
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

=== "수정 전"

    ```groovy title="genomics-2.nf" linenums="94"
    output {
        indexed_bam {
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

### 1.6. 파이프라인 다시 실행

이번에는 `-resume`과 함께 실행해봅시다.

```bash
nextflow run genomics-2.nf -resume
```

??? success "명령 출력"

    ```console
    N E X T F L O W   ~  version 25.10.2

    ┃ Launching `genomics-2.nf` [nostalgic_franklin] DSL2 - revision: f2c0a93c6a

    executor >  local (3)
    [cc/fbc705] SAMTOOLS_INDEX (3)       | 3 of 3, cached: 3 ✔
    [27/0d7eb9] GATK_HAPLOTYPECALLER (2) | 3 of 3 ✔
    ```

이번에는 작동합니다.

Nextflow 출력 자체는 다르게 보이지 않지만(일반 VCF 모드에서 성공적인 실행과 비교하여), 이제 세 샘플 모두에 대해 하위 디렉토리에 구성된 `.g.vcf` 파일과 해당 인덱스 파일을 찾을 수 있습니다.

??? abstract "디렉토리 내용 (심볼릭 링크 단축)"

    ```console
    results_genomics/
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

GVCF 파일 중 하나를 열고 스크롤하면 GATK HaplotypeCaller가 요청에 따라 GVCF 파일을 생성했음을 확인할 수 있습니다.

### 핵심 내용

좋아요, 이번 것은 Nextflow 학습 측면에서 최소한이었습니다...
하지만 process 출력 블록의 중요성을 반복할 좋은 기회였습니다!

### 다음 단계는?

채널의 내용을 수집하고 단일 입력으로 다음 process에 전달하는 방법을 배웁니다.

---

## 2. 모든 샘플에서 GVCF 데이터 수집 및 결합

이제 샘플별 모든 GVCF의 데이터를 우리가 수행하고자 하는 조인트 유전자형 분석을 지원하는 형태로 결합해야 합니다.

### 2.1. GVCF를 결합할 process 정의

준비 운동 섹션에서 이전에 수행한 작업을 상기시키자면, GVCF를 결합하는 것은 GATK 도구 GenomicsDBImport의 작업이며, 이는 소위 GenomicsDB 형식의 데이터 저장소를 생성합니다.

준비 운동 섹션에서 이전에 사용한 명령을 기반으로 이것이 어떻게 작동할지 정의하는 새 process를 작성해봅시다.

```groovy title="genomics-2.nf" linenums="66"
/*
 * GVCF를 GenomicsDB 데이터 저장소로 결합
 */
process GATK_GENOMICSDB {

    container "community.wave.seqera.io/library/gatk4:4.5.0.0--730ee8817e436867"

    input:
    path all_gvcfs
    path all_idxs
    path interval_list
    val cohort_name

    output:
    path "${cohort_name}_gdb"

    script:
    """
    gatk GenomicsDBImport \
        -V ${all_gvcfs} \
        -L ${interval_list} \
        --genomicsdb-workspace-path ${cohort_name}_gdb
    """
}
```

어떻게 생각하시나요, 합리적으로 보이나요?

연결해서 무슨 일이 일어나는지 봅시다.

### 2.2. 기본값과 함께 `cohort_name` 매개변수 추가

코호트에 대한 임의의 이름을 제공해야 합니다.
교육 시리즈 후반부에서 이러한 종류의 작업에 샘플 메타데이터를 사용하는 방법을 배우게 되지만, 지금은 `params`를 사용하여 CLI 매개변수를 선언하고 편의를 위해 기본값을 제공합니다.

```groovy title="genomics-2.nf" linenums="16"
    // 최종 출력 파일의 기본 이름
    cohort_name: String = "family_trio"
```

### 2.3. 샘플 전체에서 GATK_HAPLOTYPECALLER의 출력 수집

`GATK_HAPLOTYPECALLER` process의 출력 채널을 그대로 연결하면 Nextflow는 각 샘플 GVCF에 대해 process를 개별적으로 호출합니다.
그러나 세 개의 GVCF(및 인덱스 파일) 모두를 Nextflow가 하나의 process 호출에 함께 전달하는 방식으로 묶고 싶습니다.

좋은 소식: `collect()` 채널 연산자를 사용하여 이를 수행할 수 있습니다. GATK_HAPLOTYPECALLER 호출 직후 `workflow` 본문에 다음 줄을 추가해봅시다:

```groovy title="genomics-2.nf" linenums="118"
// 샘플 전체에서 변이 콜링 출력 수집
all_gvcfs_ch = GATK_HAPLOTYPECALLER.out.vcf.collect()
all_idxs_ch = GATK_HAPLOTYPECALLER.out.idx.collect()
```

좀 복잡해 보이나요? 이를 분해하여 일반 언어로 번역해봅시다.

1. `.out` 속성을 사용하여 참조되는 `GATK_HAPLOTYPECALLER` process의 출력 채널을 가져옵니다.
2. 채널에서 나오는 각 '요소'는 파일 쌍입니다: GVCF와 인덱스 파일이 process 출력 블록에 나열된 순서대로 있습니다. 편리하게도 마지막 세션에서 이 process의 출력 이름을 지정했기 때문에(`emit:` 사용), `.out` 속성 뒤에 `.vcf`를 추가하여 한편으로는 GVCF를, 다른 한편으로는 `.idx`를 추가하여 인덱스 파일을 선택할 수 있습니다. 이러한 출력 이름을 지정하지 않았다면 각각 `.out[0]`과 `.out[1]`로 참조해야 했을 것입니다.
3. `collect()` 채널 연산자를 추가하여 모든 GVCF 파일을 `all_gvcfs_ch`라는 새 채널의 단일 요소로 묶고, 인덱스 파일에 대해서도 동일하게 수행하여 `all_idxs_ch`라는 새 채널을 형성합니다.

!!! tip

    여기서 정확히 무슨 일이 일어나고 있는지 상상하기 어렵다면, `view()` 연산자를 사용하여 채널 연산자를 적용하기 전후에 채널 내용을 검사할 수 있다는 것을 기억하세요.

결과 `all_gvcfs_ch` 및 `all_idxs_ch` 채널은 방금 작성한 `GATK_GENOMICSDB` process에 연결할 것입니다.

!!! note

    궁금하실 경우를 대비하여, GATK GenomicsDBImport 명령이 GVCF 파일 경로만 보기를 원하기 때문에 GVCF와 인덱스 파일을 별도로 수집합니다. 다행히 Nextflow가 실행을 위해 모든 파일을 함께 스테이징하므로 파트 1에서 BAM과 인덱스에 대해 했던 것처럼 파일 순서를 걱정할 필요가 없습니다.

### 2.4. GATK_GENOMICSDB를 실행하기 위해 workflow 블록에 호출 추가

process가 있고 입력 채널이 있습니다. process 호출만 추가하면 됩니다.

```groovy title="genomics-2.nf" linenums="122"
    // GVCF를 GenomicsDB 데이터 저장소로 결합
    GATK_GENOMICSDB(
        all_gvcfs_ch,
        all_idxs_ch,
        intervals_file,
        params.cohort_name
    )
```

좋아요, 모든 것이 연결되었습니다.

### 2.5. 워크플로 실행

이것이 작동하는지 봅시다.

```bash
nextflow run genomics-2.nf -resume
```

??? failure "명령 출력"

    ```console
    N E X T F L O W   ~  version 25.10.2

    ┃ Launching `genomics-2.nf` [disturbed_bell] DSL2 - revision: 57942246cc

    executor >  local (1)
    [f1/8d8486] SAMTOOLS_INDEX (1)       | 3 of 3, cached: 3 ✔
    [e4/4ed55e] GATK_HAPLOTYPECALLER (2) | 3 of 3, cached: 3 ✔
    [51/d350ea] GATK_GENOMICSDB          | 0 of 1
    ERROR ~ Error executing process > 'GATK_GENOMICSDB'

    Caused by:
      Process `GATK_GENOMICSDB` terminated with an error exit status (1)

    Command executed:

      gatk GenomicsDBImport         -V reads_son.bam.g.vcf reads_father.bam.g.vcf reads_mother.bam.g.vcf         -L intervals.bed         --genomicsdb-workspace-path family_trio_gdb
    ```

`-resume`으로 실행하고 있으므로 상당히 빠르게 실행되지만 실패합니다!

아. 밝은 면에서, Nextflow가 `GATK_GENOMICSDB` process를 선택했고 특히 한 번만 호출했음을 볼 수 있습니다.
이는 `collect()` 접근 방식이 어느 정도 작동했음을 시사합니다.
하지만 큰 문제가 있습니다. process 호출이 실패했습니다.

위의 콘솔 출력을 자세히 살펴보면 실행된 명령이 올바르지 않음을 알 수 있습니다.

오류를 발견할 수 있나요?
이 부분을 보세요: `-V reads_father.bam.g.vcf reads_son.bam.g.vcf reads_mother.bam.g.vcf`

단일 `-V` 인수에 대해 `gatk GenomicsDBImport`에 여러 GVCF 파일을 제공했지만, 도구는 각 GVCF 파일에 대해 별도의 `-V` 인수를 예상합니다.

복습하자면, 컨테이너에서 실행한 명령은 다음과 같습니다:

```bash
gatk GenomicsDBImport \
    -V reads_mother.g.vcf \
    -V reads_father.g.vcf \
    -V reads_son.g.vcf \
    -L /data/ref/intervals.bed \
    --genomicsdb-workspace-path family_trio_gdb
```

즉, GVCF 파일 번들을 적절한 형식의 명령 문자열로 변환해야 합니다.

### 2.6. 각 입력 GVCF에 대해 별도의 `-V` 인수를 사용하여 명령줄 구성

여기서 Groovy를 기반으로 하는 Nextflow가 유용한데, 필요한 명령 문자열을 구성하기 위해 상당히 간단한 문자열 조작을 사용할 수 있기 때문입니다.

특히 이 구문을 사용합니다: `all_gvcfs.collect { gvcf -> "-V ${gvcf}" }.join(' ')`

다시 한 번, 구성 요소로 분해해봅시다.

1. 먼저 `all_gvcfs` 입력 채널의 내용을 가져와서 `.collect()`를 적용합니다(이전처럼).
2. 이를 통해 번들의 각 개별 GVCF 파일 경로를 **클로저** `{ gvcf -> "-V ${gvcf}" }`에 전달할 수 있습니다. 여기서 `gvcf`는 해당 GVCF 파일 경로를 나타냅니다.
   클로저는 파일 경로 앞에 `-V `를 추가하는 데 사용하는 미니 함수로, `"-V ${gvcf}"` 형태입니다.
3. 그런 다음 `.join(' ')`을 사용하여 세 문자열을 단일 공백을 구분자로 연결합니다.

구체적인 예를 들면 다음과 같습니다:

1. 세 개의 파일이 있습니다:

   `[A.ext, B.ext, C.ext]`

2. 클로저가 각각을 수정하여 문자열을 생성합니다:

   `"-V A.ext", "-V B.ext", "-V C.ext"`

3. `.join(' ')` 연산이 최종 문자열을 생성합니다:

   `"-V A.ext -V B.ext -V C.ext"`

이 문자열이 있으면 `def` 키워드로 정의된 로컬 변수 `gvcfs_line`에 할당할 수 있습니다:

`def gvcfs_line = all_gvcfs.collect { gvcf -> "-V ${gvcf}" }.join(' ')`

좋아요, 문자열 조작 작업이 있습니다. 어디에 넣을까요?

GVCF 파일 경로를 process로 채널링한 _후에_ 수행하고 싶기 때문에 process 정의 내부 어딘가에 이것을 넣고 싶습니다.
Nextflow가 파일 자체를 실행을 위해 올바르게 스테이징하려면 파일 경로로 보아야 하기 때문입니다.

하지만 process의 _어디에_ 이것을 추가할 수 있을까요?

재미있는 사실: `script:` 뒤와 `"""` 앞에 임의의 코드를 추가할 수 있습니다!

좋습니다, 거기에 문자열 조작 줄을 추가하고 생성하는 연결된 문자열을 사용하도록 `gatk GenomicsDBImport` 명령을 업데이트합시다.

=== "수정 후"

    ```groovy title="genomics-2.nf" linenums="87"  hl_lines="2 5"
        script:
        def gvcfs_line = all_gvcfs.collect { gvcf -> "-V ${gvcf}" }.join(' ')
        """
        gatk GenomicsDBImport \
            ${gvcfs_line} \
            -L ${interval_list} \
            --genomicsdb-workspace-path ${cohort_name}_gdb
        """
    ```

=== "수정 전"

    ```groovy title="genomics-2.nf" linenums="87"  hl_lines="4"
        script:
        """
        gatk GenomicsDBImport \
            -V ${all_gvcfs} \
            -L ${interval_list} \
            --genomicsdb-workspace-path ${cohort_name}_gdb
        """
    ```

이것이 `gatk GenomicsDBImport`에 입력을 올바르게 제공하는 데 필요한 전부일 것입니다.

!!! tip

    `gatk GenomicsDBImport` 명령을 업데이트할 때 `${gvcfs_line}` 변수를 교체할 때 `-V ` 접두사를 제거했는지 확인하세요.

### 2.7. 워크플로를 실행하여 예상대로 GenomicsDB 출력을 생성하는지 확인

좋습니다, 이것이 문제를 해결했는지 봅시다.

```bash
nextflow run genomics-2.nf -resume
```

??? success "명령 출력"

    ```console
    N E X T F L O W   ~  version 25.10.2

    ┃ Launching `genomics-2.nf` [peaceful_gates] DSL2 - revision: ca0bf847ed

    executor >  local (1)
    [cc/fbc705] SAMTOOLS_INDEX (3)       | 3 of 3, cached: 3 ✔
    [27/0d7eb9] GATK_HAPLOTYPECALLER (2) | 3 of 3, cached: 3 ✔
    [76/d13861] GATK_GENOMICSDB          | 1 of 1 ✔
    ```

아하! 이제 작동하는 것 같습니다.

처음 두 단계는 성공적으로 건너뛰었고 세 번째 단계는 이번에 완벽하게 작동했습니다.
GenomicsDB 데이터 저장소는 작업 디렉토리에 생성되지만 조인트 유전자형 분석에 사용할 중간 형식일 뿐이므로 결과에 게시되지 않습니다.

참고로, 출력이 단일 파일이 아닌 디렉토리인 것을 처리하기 위해 특별한 작업을 수행할 필요가 없었습니다.

### 핵심 내용

이제 채널에서 출력을 수집하고 다른 process에 단일 입력으로 묶는 방법을 알게 되었습니다.
또한 적절한 구문으로 주어진 도구에 입력을 제공하기 위한 명령줄을 구성하는 방법도 알게 되었습니다.

### 다음 단계는?

동일한 process의 일부로 두 번째 명령을 추가하는 방법을 배웁니다.

---

## 3. 동일한 process의 일부로 조인트 유전자형 분석 단계 실행

이제 결합된 유전체 변이 콜이 있으므로 조인트 유전자형 분석 도구를 실행할 수 있으며, 이는 실제로 우리가 관심 있는 최종 출력인 코호트 수준 변이 콜의 VCF를 생성합니다.

물류상의 이유로 동일한 process 내에 조인트 유전자형 분석을 포함하기로 결정했습니다.

### 3.1. process 이름을 GATK_GENOMICSDB에서 GATK_JOINTGENOTYPING으로 변경

process가 하나 이상의 도구를 실행할 것이므로 단일 도구 이름이 아닌 전체 작업을 참조하도록 이름을 변경합니다.

=== "수정 후"

    ```groovy title="genomics-2.nf"
    /*
     * GVCF를 GenomicsDB 데이터 저장소로 결합하고 조인트 유전자형 분석을 실행하여 코호트 수준 콜 생성
     */
    process GATK_JOINTGENOTYPING {
    ```

=== "수정 전"

    ```groovy title="genomics-2.nf"
    /*
     * GVCF를 GenomicsDB 데이터 저장소로 결합
     */
    process GATK_GENOMICSDB {
    ```

가독성을 극대화하기 위해 process 이름을 가능한 한 설명적으로 유지하는 것을 기억하세요. 동료와 미래의 자신을 위해서입니다!

### 3.2. GATK_JOINTGENOTYPING process에 조인트 유전자형 분석 명령 추가

script 섹션 내부의 첫 번째 명령 뒤에 두 번째 명령을 추가하기만 하면 됩니다.

=== "수정 후"

    ```groovy title="genomics-2.nf" linenums="89"  hl_lines="6-10"
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
    ```

=== "수정 전"

    ```groovy title="genomics-2.nf" linenums="89"
        """
        gatk GenomicsDBImport \
            ${gvcfs_line} \
            -L ${interval_list} \
            --genomicsdb-workspace-path ${cohort_name}_gdb
        """
    ```

두 명령은 터미널에서 수동으로 실행하는 것과 동일한 방식으로 순차적으로 실행됩니다.

### 3.3. GATK_JOINTGENOTYPING process 입력 정의에 참조 유전체 파일 추가

두 번째 명령은 참조 유전체 파일이 필요하므로 process 입력에 추가해야 합니다.

=== "수정 후"

    ```groovy title="genomics-2.nf" linenums="78"  hl_lines="6-8"
        input:
        path all_gvcfs
        path all_idxs
        path interval_list
        val cohort_name
        path ref_fasta
        path ref_index
        path ref_dict
    ```

=== "수정 전"

    ```groovy title="genomics-2.nf" linenums="78"
        input:
        path all_gvcfs
        path all_idxs
        path interval_list
        val cohort_name
    ```

입력하기 귀찮아 보일 수 있지만, 한 번만 입력하면 백만 번 워크플로를 실행할 수 있습니다. 그만한 가치가 있나요?

### 3.4. 코호트 수준 변이 콜의 VCF를 출력하도록 process 출력 정의 업데이트

물류상의 이유로만 존재하는 중간 형식인 GenomicsDB 데이터 저장소를 저장하는 데 실제로 관심이 없으므로 원하면 출력 블록에서 제거할 수 있습니다.

실제로 관심 있는 출력은 조인트 유전자형 분석 명령에 의해 생성된 VCF입니다.

=== "수정 후"

    ```groovy title="genomics-2.nf" linenums="87" hl_lines="2 3"
        output:
        path "${cohort_name}.joint.vcf"     , emit: vcf
        path "${cohort_name}.joint.vcf.idx" , emit: idx
    ```

=== "수정 전"

    ```groovy title="genomics-2.nf" linenums="87"
        output:
        path "${cohort_name}_gdb"
    ```

거의 다 됐습니다!

### 3.5. process 호출을 GATK_GENOMICSDB에서 GATK_JOINTGENOTYPING으로 업데이트

workflow 본문에서 process 호출 이름을 GATK_GENOMICSDB에서 GATK_JOINTGENOTYPING으로 변경하는 것을 잊지 맙시다. 그리고 조인트 유전자형 분석 도구에 제공해야 하므로 참조 유전체 파일도 입력으로 추가해야 합니다.

=== "수정 후"

    ```groovy title="genomics-2.nf" linenums="126"
    // GVCF를 GenomicsDB 데이터 저장소로 결합하고 조인트 유전자형 분석 적용
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

=== "수정 전"

    ```groovy title="genomics-2.nf" linenums="126"
    // GVCF를 GenomicsDB 데이터 저장소로 결합
    GATK_GENOMICSDB(
        all_gvcfs_ch,
        all_idxs_ch,
        intervals_file,
        params.cohort_name
    )
    ```

이제 process가 완전히 연결되었습니다.

### 3.6. 조인트 VCF를 publish 섹션에 추가

새 process의 조인트 VCF 출력을 게시해야 합니다.
워크플로의 `publish:` 섹션에 다음 줄을 추가하세요:

```groovy title="genomics-2.nf" linenums="145"
    joint_vcf = GATK_JOINTGENOTYPING.out.vcf
    joint_vcf_idx = GATK_JOINTGENOTYPING.out.idx
```

### 3.7. output 블록에 조인트 VCF 대상 추가

마지막
