# 파트 1: 데모 파이프라인 실행하기

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } AI 지원 번역 - [자세히 알아보기 및 개선 제안](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

Hello nf-core 교육 과정의 첫 번째 파트에서는 nf-core 파이프라인을 찾아서 사용해보고, 필요에 맞게 실행을 설정 및 맞춤화하며, 입력 검증이 일반적인 오류를 어떻게 방지하는지 학습합니다.

nf-core 프로젝트에서 코드 구조와 도구 작동을 시연하기 위한 파이프라인 모음의 일부로 유지 관리하는 nf-core/demo라는 파이프라인을 사용하겠습니다.

[시작하기](./00_orientation.md) 페이지의 안내에 따라 작업 디렉토리를 `hello-nf-core/`로 설정했는지 확인하십시오.

---

## 1. nf-core/demo 파이프라인 찾기 및 가져오기

[nf-co.re](https://nf-co.re) 프로젝트 웹사이트에서 nf-core/demo 파이프라인을 찾는 것부터 시작합니다. 이 웹사이트는 일반 문서 및 도움말 기사, 각 파이프라인에 대한 문서, 블로그 게시물, 이벤트 공지 등 모든 정보를 중앙에서 관리합니다.

### 1.1. 웹사이트에서 파이프라인 찾기

웹 브라우저에서 [https://nf-co.re/pipelines/](https://nf-co.re/pipelines/)로 이동하여 검색창에 `demo`를 입력하십시오.

![검색 결과](./img/search-results.png)

파이프라인 이름인 `demo`를 클릭하여 파이프라인 문서 페이지에 접근하십시오.

릴리스된 각 파이프라인에는 다음 문서 섹션이 포함된 전용 페이지가 있습니다:

- **Introduction:** 파이프라인의 소개 및 개요
- **Usage:** 파이프라인 실행 방법에 대한 설명
- **Parameters:** 설명과 함께 그룹화된 파이프라인 매개변수
- **Output:** 예상되는 출력 파일에 대한 설명 및 예제
- **Results:** 전체 테스트 데이터셋에서 생성된 출력 파일 예제
- **Releases & Statistics:** 파이프라인 버전 기록 및 통계

새로운 파이프라인 도입을 고려할 때는 실행을 시도하기 전에 파이프라인 문서를 주의 깊게 읽어 파이프라인이 무엇을 하는지, 어떻게 설정해야 하는지 이해해야 합니다.

지금 살펴보시고 다음을 찾을 수 있는지 확인해보십시오:

- 파이프라인이 실행할 도구들 (`Introduction` 탭 확인)
- 파이프라인이 받거나 요구하는 입력 및 매개변수 (`Parameters` 탭 확인)
- 파이프라인이 생성하는 출력 (`Output` 탭 확인)

#### 1.1.1. 파이프라인 개요

`Introduction` 탭은 시각적 표현(지하철 노선도라고 함)과 파이프라인의 일부로 실행되는 도구 목록을 포함한 파이프라인의 개요를 제공합니다.

![파이프라인 지하철 노선도](./img/nf-core-demo-subway-cropped.png)

1. Read QC (FASTQC)
2. Adapter and quality trimming (SEQTK_TRIM)
3. Present QC for raw reads (MULTIQC)

#### 1.1.2. 명령줄 예제

문서는 또한 입력 파일 예제(아래에서 자세히 설명)와 명령줄 예제를 제공합니다.

```bash
nextflow run nf-core/demo \
  -profile <docker/singularity/.../institute> \
  --input samplesheet.csv \
  --outdir <OUTDIR>
```

예제 명령이 워크플로우 파일을 지정하지 않고 파이프라인 저장소에 대한 참조인 `nf-core/demo`만 지정한다는 것을 알 수 있습니다.

이런 방식으로 호출되면 Nextflow는 코드가 특정 방식으로 구성되어 있다고 가정합니다.
이 구조를 살펴볼 수 있도록 코드를 가져오겠습니다.

### 1.2. 파이프라인 코드 가져오기

파이프라인이 우리의 목적에 적합해 보인다고 판단했다면 사용해보겠습니다.
다행히 Nextflow는 수동으로 다운로드할 필요 없이 올바르게 형식화된 저장소에서 파이프라인을 쉽게 가져올 수 있게 해줍니다.

#### 1.2.1. `nextflow pull` 사용하기

터미널로 돌아가서 다음을 실행합니다:

```bash
nextflow pull nf-core/demo
```

??? success "명령 출력"

    ```console
    Checking nf-core/demo ...
    downloaded from https://github.com/nf-core/demo.git - revision: 04060b4644 [master]
    ```

Nextflow는 파이프라인 코드를 `pull`하여 전체 저장소를 로컬 드라이브에 다운로드합니다.

명확히 하자면, nf-core 파이프라인뿐만 아니라 GitHub에서 적절하게 설정된 모든 Nextflow 파이프라인에 대해 이 작업을 수행할 수 있습니다.
그러나 nf-core는 가장 큰 오픈소스 Nextflow 파이프라인 모음입니다.

#### 1.2.2. `nextflow list` 사용하기

Nextflow에게 이런 방식으로 가져온 파이프라인 목록을 제공하도록 할 수 있습니다:

```bash
nextflow list
```

??? success "명령 출력"

    ```console
    nf-core/demo
    ```

다른 파이프라인을 몇 개 더 pull하여 여러 개가 있을 때 어떻게 나열되는지 확인해볼 수 있습니다.

#### 1.2.3. `$NXF_HOME/assets/`에서 파이프라인 찾기

파일들이 현재 작업 디렉토리에 없다는 것을 알 수 있습니다.
기본적으로 Nextflow는 `$NXF_HOME/assets`에 저장합니다.

```bash
tree -L 2 $NXF_HOME/assets/
```

```console title="Directory contents"
/workspaces/.nextflow/assets/
└── nf-core
    └── demo

2 directories, 0 files
```

!!! note "참고"

    교육 환경을 사용하지 않는 경우 시스템에서 전체 경로가 다를 수 있습니다.

Nextflow는 이러한 파이프라인을 직접 상호작용하는 코드가 아닌 라이브러리처럼 사용해야 한다는 원칙에 따라 다운로드된 소스 코드를 의도적으로 '방해가 되지 않는' 위치에 보관합니다.

#### 1.2.4. 소스 코드에 쉽게 접근하기 위한 심볼릭 링크 만들기

코드를 자세히 살펴보지는 않겠지만, 전체적인 구성이 어떻게 되어 있는지 간략히 확인해보겠습니다.

파이프라인 소스 코드를 더 쉽게 탐색할 수 있도록 assets 디렉토리에 대한 심볼릭 링크를 만드십시오:

```bash
ln -s $NXF_HOME/assets pipelines
```

이렇게 하면 `tree -L 2 pipelines`로 코드를 탐색하거나 파일을 직접 열 수 있는 바로가기가 만들어집니다.

#### 1.2.5. 코드 구성 개요

`tree`를 사용하거나 파일 탐색기를 사용하여 `nf-core/demo` 디렉토리를 찾아 열 수 있습니다.

```bash
tree -L 1 pipelines/nf-core/demo
```

??? abstract "디렉토리 내용"

    ```console
    pipelines/nf-core/demo
    ├── assets
    ├── CHANGELOG.md
    ├── CITATIONS.md
    ├── CODE_OF_CONDUCT.md
    ├── conf
    ├── docs
    ├── LICENSE
    ├── main.nf
    ├── modules
    ├── modules.json
    ├── nextflow.config
    ├── nextflow_schema.json
    ├── nf-test.config
    ├── README.md
    ├── ro-crate-metadata.json
    ├── subworkflows
    ├── tests
    ├── tower.yml
    └── workflows
    ```

보시다시피 많은 내용이 있지만, 대부분은 신경 쓰지 않아도 됩니다.

간략히 살펴보면, 최상위 레벨에서 요약 정보가 있는 README 파일과 라이선싱, 기여 가이드라인, 인용 및 행동 강령과 같은 프로젝트 정보를 요약하는 부속 파일을 찾을 수 있습니다.
상세한 파이프라인 문서는 `docs` 디렉토리에 있습니다.
이 모든 콘텐츠는 nf-core 웹사이트의 웹 페이지를 프로그래밍 방식으로 생성하는 데 사용되므로 항상 코드와 최신 상태를 유지합니다.

나머지 부분에서는 세 가지 기능적 그룹으로 구분할 수 있습니다:

1. 파이프라인 코드 구성요소 (`main.nf`, `workflows`, `subworkflows`, `modules`)
2. 파이프라인 설정
3. 파이프라인 매개변수 / 입력 및 검증

이 파트에서는 파이프라인 코드 구성요소를 자세히 다루지 않겠지만, nf-core 파이프라인의 최종 사용자로서 관련이 있을 설정 및 검증 요소를 살펴보겠습니다.

!!! tip "팁"

    nf-core 파이프라인의 소스 코드는 GitHub에서도 탐색할 수 있습니다. 예: [github.com/nf-core/demo](https://github.com/nf-core/demo).
    모든 nf-core 파이프라인은 동일한 디렉토리 구조를 따르므로, 구조를 한 번 파악하면 어떤 파이프라인에서도 같은 방식으로 설정 파일, 모듈, 워크플로우를 찾을 수 있습니다.

이제 파이프라인을 실행해보겠습니다!

### 핵심 정리

이제 nf-core 웹사이트를 통해 파이프라인을 찾고 소스 코드의 로컬 사본을 가져오는 방법을 알게 되었습니다.

### 다음 단계

최소한의 노력으로 nf-core 파이프라인을 사용해보는 방법을 배워보겠습니다.

---

## 2. 테스트 프로파일로 파이프라인 사용해보기

편리하게도 모든 nf-core 파이프라인에는 test 프로파일이 함께 제공됩니다.
이것은 [nf-core/test-datasets](https://github.com/nf-core/test-datasets) 저장소에서 호스팅되는 작은 테스트 데이터셋을 사용하여 파이프라인을 실행하기 위한 최소한의 설정 모음입니다.
작은 규모로 파이프라인을 빠르게 사용해볼 수 있는 좋은 방법입니다.

!!! note "참고"

    Nextflow의 configuration profile 시스템을 사용하면 다양한 컨테이너 엔진이나 실행 환경 간에 쉽게 전환할 수 있습니다.
    자세한 내용은 [Hello Nextflow Part 6: Configuration](../hello_nextflow/06_hello_config.md)을 참조하십시오.

### 2.1. 테스트 프로파일 살펴보기

파이프라인을 실행하기 전에 파이프라인의 test 프로파일이 무엇을 지정하는지 확인하는 것이 좋습니다.
`nf-core/demo`의 `test` 프로파일은 설정 파일 `conf/test.config`에 있습니다.
`nextflow pull`로 다운로드한 파이프라인 소스 내에서 로컬로 찾을 수 있습니다:

```bash
code $NXF_HOME/assets/nf-core/demo/conf/test.config
```

해당 파일의 내용은 다음과 같습니다:

```groovy title="conf/test.config" linenums="1" hl_lines="8 26"
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Nextflow config file for running minimal tests
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Defines input files and everything required to run a fast and simple pipeline test.

    Use as follows:
        nextflow run nf-core/demo -profile test,<docker/singularity> --outdir <OUTDIR>

----------------------------------------------------------------------------------------
*/

process {
    resourceLimits = [
        cpus: 2,
        memory: '4.GB',
        time: '1.h'
    ]
}

params {
    config_profile_name        = 'Test profile'
    config_profile_description = 'Minimal test dataset to check pipeline function'

    // 입력 데이터
    input  = 'https://raw.githubusercontent.com/nf-core/test-datasets/viralrecon/samplesheet/samplesheet_test_illumina_amplicon.csv'

}
```

상단의 주석 블록에 이 테스트 프로파일로 파이프라인을 실행하는 방법을 보여주는 사용 예제가 포함되어 있음을 바로 알 수 있습니다.

```groovy title="conf/test.config" linenums="7"
    Use as follows:
        nextflow run nf-core/demo -profile test,<docker/singularity> --outdir <OUTDIR>
```

제공해야 하는 것은 예제 명령에서 꺾쇠괄호 안에 표시된 것뿐입니다: `<docker/singularity>`와 `<OUTDIR>`.

다시 말하지만, `<docker/singularity>`는 컨테이너 시스템의 선택을 의미합니다. 모든 nf-core 파이프라인은 재현성을 보장하고 소프트웨어 설치 문제를 제거하기 위해 컨테이너(Docker, Singularity 등)와 함께 사용할 수 있도록 설계되었습니다.
따라서 파이프라인을 테스트하기 위해 Docker 또는 Singularity를 사용할지 지정해야 합니다.

`--outdir <OUTDIR>` 부분은 Nextflow가 파이프라인의 출력을 작성할 디렉토리를 의미합니다.
직접 만들 수 있는 이름을 제공해야 합니다.
이미 존재하지 않는 경우 Nextflow가 실행 시 생성합니다.

주석 블록 다음 섹션으로 이동하면, test 프로파일이 테스트를 위해 미리 설정된 내용을 보여줍니다: 가장 주목할 점은 `input` 매개변수가 이미 테스트 데이터셋을 가리키도록 설정되어 있으므로 자체 데이터를 제공할 필요가 없다는 것입니다.
미리 설정된 입력 링크를 따라가면 여러 실험 샘플에 대한 샘플 식별자와 파일 경로가 포함된 csv 파일임을 알 수 있습니다.

```csv title="samplesheet_test_illumina_amplicon.csv"
sample,fastq_1,fastq_2
SAMPLE1_PE,https://raw.githubusercontent.com/nf-core/test-datasets/viralrecon/illumina/amplicon/sample1_R1.fastq.gz,https://raw.githubusercontent.com/nf-core/test-datasets/viralrecon/illumina/amplicon/sample1_R2.fastq.gz
SAMPLE2_PE,https://raw.githubusercontent.com/nf-core/test-datasets/viralrecon/illumina/amplicon/sample2_R1.fastq.gz,https://raw.githubusercontent.com/nf-core/test-datasets/viralrecon/illumina/amplicon/sample2_R2.fastq.gz
SAMPLE3_SE,https://raw.githubusercontent.com/nf-core/test-datasets/viralrecon/illumina/amplicon/sample1_R1.fastq.gz,
SAMPLE3_SE,https://raw.githubusercontent.com/nf-core/test-datasets/viralrecon/illumina/amplicon/sample2_R1.fastq.gz,
```

이것은 샘플시트라고 하며, nf-core 파이프라인에 대한 가장 일반적인 입력 형식입니다.

!!! note "참고"

    데이터 형식과 유형에 익숙하지 않더라도 걱정하지 마십시오. 이후 내용에 중요하지 않습니다.

따라서 파이프라인을 사용해보는 데 필요한 모든 것이 있음을 확인했습니다.

### 2.2. 파이프라인 실행하기

컨테이너 시스템으로 Docker를, 출력 디렉토리로 `demo-results`를 사용하기로 결정했다면, 테스트 명령을 실행할 준비가 된 것입니다:

```bash
nextflow run nf-core/demo -profile docker,test --outdir demo-results
```

??? success "명령 출력"

    ```console
     N E X T F L O W   ~  version 25.10.4

    Launching `https://github.com/nf-core/demo` [magical_pauling] DSL2 - revision: 45904cb9d1 [master]


    ------------------------------------------------------
                                            ,--./,-.
            ___     __   __   __   ___     /,-._.--~'
      |\ | |__  __ /  ` /  \ |__) |__         }  {
      | \| |       \__, \__/ |  \ |___     \`-._,-`-,
                                            `._,._,'
      nf-core/demo 1.1.0
    ------------------------------------------------------
    Input/output options
      input                     : https://raw.githubusercontent.com/nf-core/test-datasets/viralrecon/samplesheet/samplesheet_test_illumina_amplicon.csv
      outdir                    : demo-results

    Institutional config options
      config_profile_name       : Test profile
      config_profile_description: Minimal test dataset to check pipeline function

    Generic options
      trace_report_suffix       : 2025-11-21_04-57-41

    Core Nextflow options
      revision                  : master
      runName                   : magical_pauling
      containerEngine           : docker
      launchDir                 : /workspaces/training/hello-nf-core
      workDir                   : /workspaces/training/hello-nf-core/work
      projectDir                : /workspaces/.nextflow/assets/nf-core/demo
      userName                  : root
      profile                   : docker,test
      configFiles               : /workspaces/.nextflow/assets/nf-core/demo/nextflow.config

    !! Only displaying parameters that differ from the pipeline defaults !!
    ------------------------------------------------------
    * The pipeline
        https://doi.org/10.5281/zenodo.12192442

    * The nf-core framework
        https://doi.org/10.1038/s41587-020-0439-x

    * Software dependencies
        https://github.com/nf-core/demo/blob/master/CITATIONS.md


    executor >  local (7)
    [ff/a6976b] NFCORE_DEMO:DEMO:FASTQC (SAMPLE3_SE)     | 3 of 3 ✔
    [39/731ab7] NFCORE_DEMO:DEMO:SEQTK_TRIM (SAMPLE3_SE) | 3 of 3 ✔
    [7c/78d96e] NFCORE_DEMO:DEMO:MULTIQC                 | 1 of 1 ✔
    -[nf-core/demo] Pipeline completed successfully-
    ```

출력이 일치한다면 축하합니다! 첫 번째 nf-core 파이프라인을 실행했습니다.

기본 Nextflow 파이프라인을 실행할 때보다 콘솔 출력이 훨씬 많다는 것을 알 수 있습니다.
파이프라인의 버전, 입력 및 출력, 그리고 몇 가지 설정 요소의 요약이 포함된 헤더가 있습니다.

!!! note "참고"

    출력에는 다른 타임스탬프, 실행 이름 및 파일 경로가 표시되지만 전체 구조와 프로세스 실행은 유사해야 합니다.

출력 상단 근처의 다음 줄을 확인하십시오:

```console
Launching `https://github.com/nf-core/demo` [magical_pauling] DSL2 - revision: 45904cb9d1 [master]
```

이 줄은 어떤 리비전의 파이프라인이 사용되었는지 알려줍니다.
버전을 지정하지 않았으므로 Nextflow는 `master`의 최신 커밋을 사용했습니다.
재현 가능한 실행을 위해서는 `-r` 플래그를 사용하여 특정 릴리스를 지정해야 합니다:

```bash
nextflow run nf-core/demo -r 1.1.0 -profile docker,test --outdir demo-results
```

이렇게 하면 새로운 커밋이나 릴리스에 관계없이 항상 동일한 파이프라인 코드가 사용됩니다.
이 교육에서는 간결함을 위해 `-r`을 생략하지만, 실제 운영 환경에서는 항상 지정해야 합니다.

실행 출력으로 이동하여 어떤 프로세스가 실행되었는지 알려주는 줄을 살펴보겠습니다:

```console
executor >  local (7)
[ff/a6976b] NFCORE_DEMO:DEMO:FASTQC (SAMPLE3_SE)     | 3 of 3 ✔
[39/731ab7] NFCORE_DEMO:DEMO:SEQTK_TRIM (SAMPLE3_SE) | 3 of 3 ✔
[7c/78d96e] NFCORE_DEMO:DEMO:MULTIQC                 | 1 of 1 ✔
-[nf-core/demo] Pipeline completed successfully-
```

이것은 nf-core 웹사이트의 파이프라인 문서 페이지에 표시된 세 가지 도구에 해당하는 세 개의 프로세스가 실행되었음을 알려줍니다: FASTQC, SEQTK_TRIM 및 MULTIQC.

여기에 표시된 `NFCORE_DEMO:DEMO:MULTIQC`와 같은 전체 프로세스 이름은 Hello Nextflow 입문 자료에서 본 것보다 깁니다.
여기에는 상위 워크플로우의 이름이 포함되어 있으며 파이프라인 코드의 모듈성을 반영합니다.
이 과정의 파트 2에서 이에 대해 자세히 설명하겠습니다.

### 2.3. 파이프라인 출력 살펴보기

마지막으로 파이프라인이 생성한 `demo-results` 디렉토리를 살펴보겠습니다.

```bash
tree -L 2 demo-results
```

??? abstract "디렉토리 내용"

    ```console
    demo-results
    ├── fastqc
    │   ├── SAMPLE1_PE
    │   ├── SAMPLE2_PE
    │   └── SAMPLE3_SE
    ├── fq
    │   ├── SAMPLE1_PE
    │   ├── SAMPLE2_PE
    │   └── SAMPLE3_SE
    ├── multiqc
    │   ├── multiqc_data
    │   ├── multiqc_plots
    │   └── multiqc_report.html
    └── pipeline_info
        ├── execution_report_2025-11-21_04-57-41.html
        ├── execution_timeline_2025-11-21_04-57-41.html
        ├── execution_trace_2025-11-21_04-57-41.txt
        ├── nf_core_demo_software_mqc_versions.yml
        ├── params_2025-11-21_04-57-46.json
        └── pipeline_dag_2025-11-21_04-57-41.html
    ```

많아 보일 수 있습니다.
`nf-core/demo` 파이프라인의 출력에 대해 자세히 알아보려면 [문서 페이지](https://nf-co.re/demo/1.1.0/docs/output/)를 확인하십시오.

현 단계에서 중요한 것은 결과가 모듈별로 구성되어 있고, 파이프라인 실행에 대한 다양한 타임스탬프가 있는 보고서가 포함된 `pipeline_info`라는 디렉토리가 추가로 있다는 것입니다.

예를 들어, `execution_timeline_*` 파일은 어떤 프로세스가 실행되었는지, 어떤 순서로, 얼마나 오래 실행되었는지 보여줍니다:

![실행 타임라인 보고서](./img/execution_timeline.png)

!!! note "참고"

    여기서 작업이 병렬로 실행되지 않은 이유는 Github Codespaces의 최소 사양 머신에서 실행하고 있기 때문입니다.
    병렬 실행을 보려면 codespace의 CPU 할당과 테스트 설정의 리소스 제한을 늘려보십시오.

이러한 보고서는 모든 nf-core 파이프라인에 대해 자동으로 생성됩니다.

### 핵심 정리

내장된 test 프로파일을 사용하여 nf-core 파이프라인을 실행하는 방법과 출력을 찾을 수 있는 위치를 알게 되었습니다.

### 다음 단계

파이프라인 실행을 맞춤화하기 위한 설정 방법을 배워보겠습니다.

---

## 3. 파이프라인 실행 설정하기

[Hello Config](../hello_nextflow/06_hello_config.md)에서 설명한 것처럼, 파이프라인 코드 자체를 변경하지 않고 파이프라인이 실행할 데이터와 실행 방식을 변경할 수 있어야 합니다.
이를 위해 Nextflow는 파이프라인 설정을 제어하는 여러 방법을 지원하는데, 처음에는 다소 복잡하게 느껴질 수 있습니다.

nf-core 프로젝트는 설정 요소를 구성하기 위한 규칙을 정의하며, 최상위 레벨에서 두 가지 종류의 설정을 구분합니다: **파이프라인 매개변수**와 엄밀한 의미의 **설정**.

- **파이프라인 매개변수** (`params` 시스템을 통해 설정)는 일반적으로 입력 파일, 도구 동작 플래그 및 분석 매개변수 등을 포함합니다.
- 엄밀한 의미의 **설정**은 파이프라인이 실행되는 방식의 세부 사항, 즉 executor, 컴퓨팅 리소스 할당 등을 의미합니다.

<figure class="excalidraw">
    --8<-- "docs/en/docs/hello_nf-core/img/params_vs_config.excalidraw.svg"
</figure>

파이프라인 매개변수부터 살펴본 다음, 엄밀한 의미의 설정을 다루겠습니다.

### 3.1. 파이프라인 매개변수

모든 nf-core 파이프라인에서 `--help` 플래그를 사용하여 명령줄에서 직접 파이프라인 매개변수의 전체 목록을 확인할 수 있습니다. `--help` 자체도 파이프라인 매개변수입니다.

#### 3.1.1. `--help`로 매개변수 목록 확인하기

데모 파이프라인에 대한 도움말 명령을 실행합니다:

```bash
nextflow run nf-core/demo --help
```

??? success "명령 출력"

    ```console
     N E X T F L O W   ~  version 25.10.4

    Launching `https://github.com/nf-core/demo` [run_name] DSL2 - revision: 45904cb9d1 [master]

    ----------------------------------------------------
                                            ,--./,-.
            ___     __   __   __   ___     /,-._.--~'
      |\ | |__  __ /  ` /  \ |__) |__         }  {
      | \| |       \__, \__/ |  \ |___     \`-._,-`-,
                                            `._,._,'
      nf-core/demo 1.1.0
    ----------------------------------------------------
    Typical pipeline command:

      nextflow run nf-core/demo -profile <docker/singularity/.../institute> --input samplesheet.csv --outdir <OUTDIR>

    Input/output options
      --input                       [string]           Path to a metadata file containing information about the samples in the experiment.
      --outdir                      [string]           The output directory where the results will be saved. You have to use absolute paths to storage on Cloud infrastructure.
      --email                       [string]           Email address for completion summary.
      --multiqc_title               [string]           MultiQC report title. Printed as page header, used for filename if not otherwise specified.

    Reference genome options
      --genome                      [string]           Name of iGenomes reference.
      --fasta                       [string]           Path to FASTA genome file.

    Process skipping options
      --skip_trim                   [boolean]          Skip trimming fastq files with seqtk

    Generic options
      --multiqc_methods_description [string]           Custom MultiQC yaml file containing HTML including a methods description.
      --help                        [boolean, string]  Display the help message.
      --help_full                   [boolean]          Display the full detailed help message.
      --show_hidden                 [boolean]          Display hidden parameters in the help message (only works when --help or --help_full are provided).
     !! Hiding 20 param(s), use the `--show_hidden` parameter to show them !!
    ----------------------------------------------------

    * The pipeline
        https://doi.org/10.5281/zenodo.12192442

    * The nf-core framework
        https://doi.org/10.1038/s41587-020-0439-x

    * Software dependencies
        https://github.com/nf-core/demo/blob/master/CITATIONS.md
    ```

보시다시피 출력은 매개변수를 카테고리별로 그룹화하여(입력/출력 옵션, 레퍼런스 게놈 옵션 등) 각각의 유형과 설명을 함께 표시합니다.

이 카테고리 분류는 스키마 파일에 의해 결정되며, 아래에서 자세히 다룹니다.
일반 Nextflow 파이프라인에서는 개발자가 직접 구현한 경우에만 `--help`가 작동합니다.

!!! tip "팁"

    `--help --show_hidden`을 사용하면 `--publish_dir_mode`나 `--monochrome_logs`와 같이 기본적으로 숨겨진 추가 매개변수를 확인할 수 있습니다.

#### 3.1.2. 매개변수 값 설정하기

[Hello Config](../hello_nextflow/06_hello_config.md)에서 다룬 것처럼, 명령줄에서 `--param_name`으로 매개변수 값을 설정하거나, YAML 파일에 매개변수 집합을 모아 `-params-file`로 전달할 수 있습니다.
두 방법 모두 nf-core 파이프라인에서 동일하게 작동합니다.

예를 들어, 트리밍 단계를 건너뛰려면:

```bash
nextflow run nf-core/demo -profile docker,test --outdir demo-results-notrim --skip_trim
```

??? success "명령 출력"

    ```console
    executor >  local (4)
    [3f/a82c91] NFCORE_DEMO:DEMO:FASTQC (SAMPLE3_SE) | 3 of 3 ✔
    [7d/c5e014] NFCORE_DEMO:DEMO:MULTIQC             | 1 of 1 ✔
    -[nf-core/demo] Pipeline completed successfully-
    ```

`SEQTK_TRIM` 프로세스가 출력에 더 이상 나타나지 않습니다.

!!! info "정보"

    `-c`로 전달하는 사용자 정의 설정 파일에서 파이프라인 매개변수를 설정하는 것이 기술적으로는 가능하지만, Nextflow의 설정 우선순위 규칙에 따라 파이프라인 자체의 `nextflow.config`에 이미 설정된 기본값을 재정의하지 못할 수 있습니다.
    명령줄에서 `--param_name`을 사용하거나 `-params-file`을 사용하는 것이 더 안정적입니다. 이 방법들은 항상 우선순위가 높습니다.

    **경험 법칙:** `--help` 출력에 나타나는 매개변수는 설정 파일이 아닌 명령줄이나 params 파일을 통해 설정하십시오.

#### 3.1.3. 매개변수 검증

흥미로운 사실: nf-core 프로젝트는 개발자가 모든 파이프라인 매개변수를 JSON 스키마 파일(`nextflow_schema.json`)에 공식적으로 정의하도록 요구하기 때문에 모든 nf-core 파이프라인에서 `--help` 명령이 작동합니다.
이 스키마는 각 매개변수의 유형, 설명, 기본값 및 그룹화를 기록합니다.

`--help` 출력을 지원하는 것 외에도, 스키마 파일은 실행 시 자동 검증을 가능하게 합니다.
즉, Nextflow는 전달한 모든 매개변수가 존재하고 적절한 값(적절한 유형, 허용된 값 범위 내 등)이 지정되었는지 확인할 수 있습니다.

이에 대해서는 [파트 5: 입력 검증](05_input_validation.md)에서 자세히 다루지만, 데모 파이프라인에 잘못된 매개변수 입력을 제공하여 이미 작동하는 것을 확인할 수 있습니다.

##### 3.1.3.1. 인식되지 않는 매개변수

존재하지 않는 매개변수를 전달해보십시오:

```bash
nextflow run nf-core/demo -profile docker,test --outdir demo-results --foobar "invalid"
```

콘솔 출력에 경고가 포함됩니다:

```console
WARN: The following invalid input values have been detected:

* --foobar: invalid
```

파이프라인은 계속 실행되지만, 경고를 통해 `--foobar`가 인식되지 않는 매개변수임을 즉시 알 수 있습니다.
이를 통해 출력이 잘못된 위치에 저장된 이유를 파악하느라 컴퓨팅 시간을 낭비하기 전에 `--outDir` 대신 `--outdir`과 같은 오타를 잡아낼 수 있습니다.

##### 3.1.3.2. 잘못된 매개변수 값

검증은 매개변수 **값**도 확인합니다.
`--skip_trim` 매개변수는 boolean 플래그이므로, 문자열 값을 전달하면 파이프라인이 즉시 실패합니다:

```bash
nextflow run nf-core/demo -profile docker,test --outdir demo-results --skip_trim yes
```

```console
ERROR ~ Validation of pipeline parameters failed!

The following invalid input values have been detected:

* --skip_trim (yes): Value is [string] but should be [boolean]
```

파이프라인은 어떤 프로세스도 실행하기 전에 중단되어, 실패하거나 잘못된 실행을 방지합니다.
Boolean 매개변수는 값 없이 플래그로 전달하거나(`--skip_trim`), params 파일에서 `true`/`false`로 설정해야 합니다.

#### 3.1.4. 입력 검증

동일한 검증 로직을 입력 파일의 유효성 확인에도 사용할 수 있습니다.
예를 들어, 파이프라인이 주요 데이터 입력으로 샘플시트를 기대하는 경우(많은 nf-core 파이프라인이 그렇습니다), 개발자는 입력 파일의 구조를 설명하는 입력 스키마(매개변수 스키마와는 별개)를 제공할 수 있습니다.

그러면 런타임에 Nextflow가 제공된 입력 파일이 유효한지 확인할 수 있습니다.

이에 대해서도 [파트 5: 입력 검증](05_input_validation.md)에서 자세히 다루지만, 데모 파이프라인에 잘못된 입력 샘플시트를 제공하여 이미 작동하는 것을 확인할 수 있습니다.

`nf-core/demo` 파이프라인은 `sample`, `fastq_1`, `fastq_2` 열이 있는 CSV 파일을 기대합니다.
이는 예상 구조, 열 유형 및 제약 조건을 지정하는 스키마 파일(`assets/schema_input.json`)에 정의되어 있습니다.

??? abstract "assets/schema_input.json"

    ```json title="assets/schema_input.json"
    {
        "$schema": "https://json-schema.org/draft/2020-12/schema",
        "$id": "https://raw.githubusercontent.com/nf-core/demo/master/assets/schema_input.json",
        "title": "nf-core/demo pipeline - params.input schema",
        "description": "Schema for the file provided with params.input",
        "type": "array",
        "items": {
            "type": "object",
            "properties": {
                "sample": {
                    "type": "string",
                    "pattern": "^\\S+$",
                    "errorMessage": "Sample name must be provided and cannot contain spaces",
                    "meta": ["id"]
                },
                "fastq_1": {
                    "type": "string",
                    "format": "file-path",
                    "exists": true,
                    "pattern": "^([\\S\\s]*\\/)?[^\\s\\/]+\\.f(ast)?q\\.gz$",
                    "errorMessage": "FastQ file for reads 1 must be provided, cannot contain spaces and must have extension '.fq.gz' or '.fastq.gz'"
                },
                "fastq_2": {
                    "type": "string",
                    "format": "file-path",
                    "exists": true,
                    "pattern": "^([\\S\\s]*\\/)?[^\\s\\/]+\\.f(ast)?q\\.gz$",
                    "errorMessage": "FastQ file for reads 2 cannot contain spaces and must have extension '.fq.gz' or '.fastq.gz'"
                }
            },
            "required": ["sample", "fastq_1"]
        }
    }
    ```

스키마는 `sample`과 `fastq_1`이 필수이고 `fastq_2`는 선택 사항임을 지정합니다(페어드 엔드 및 단일 엔드 데이터 모두 지원).
파일 경로는 존재 여부와 확장자 패턴이 검증됩니다.

##### 3.1.4.1. 잘못된 샘플시트 만들기

누락된 열과 존재하지 않는 파일 경로가 있는 샘플시트를 만드십시오:

```csv title="malformed_samplesheet.csv"
sample,fastq_2
SAMPLE1,/not/a/real/file.fastq.gz
```

이 샘플시트는 필수 `fastq_1` 열이 없고 `fastq_2`에 존재하지 않는 파일 경로가 있습니다.
두 문제 모두 다음 단계에서 검증 오류를 발생시킵니다.

##### 3.1.4.2. 잘못된 샘플시트로 데모 파이프라인 실행하기

`malformed_samplesheet.csv`를 입력으로 사용하여 데모 파이프라인을 실행합니다.

```bash
nextflow run nf-core/demo -profile docker,test --outdir demo-results --input malformed_samplesheet.csv
```

```console
ERROR ~ Validation of pipeline parameters failed!

The following invalid input values have been detected:

* --input (malformed_samplesheet.csv): Validation of file failed:
    -> Entry 1: Error for field 'fastq_2' (/not/a/real/file.fastq.gz): the file or directory
       '/not/a/real/file.fastq.gz' does not exist (FastQ file for reads 2 cannot contain spaces
       and must have extension '.fq.gz' or '.fastq.gz')
    -> Entry 1: Missing required field(s): fastq_1
```

보시다시피 파이프라인은 즉시 실패하고 **모든** 검증 오류를 한 번에 보고합니다.
nf-schema는 첫 번째 오류에서 멈추지 않고 모든 문제를 수집하여 함께 나열하므로, 문제를 하나씩 발견하는 대신 한 번에 모두 수정할 수 있습니다.

각 오류는 문제를 일으킨 정확한 항목과 필드를 식별하므로, 샘플시트를 수정한 후 Nextflow가 실제로 파일 경로에 접근할 때 나중에 실패하지 않을 것이라는 확신을 가지고 파이프라인을 다시 실행할 수 있습니다.

개발자를 위한 자세한 내용은 이 과정의 [파트 5](./05_input_validation.md)에서 다룹니다.

### 3.2. 설정

엄밀한 의미의 설정은 파이프라인이 **어떻게** 실행되는지를 제어합니다: 리소스 할당, 도구별 인자, 작업이 실행되는 위치, 사용할 소프트웨어 패키징 시스템 등입니다.

nf-core 파이프라인은 `nextflow.config`와 `conf/` 디렉토리에 기본 설정을 포함합니다.
무언가를 재정의하기 전에 기본값이 어디에 있는지 파악하는 것이 도움이 됩니다.

섹션 2.1에서 파이프라인 소스 코드가 `$NXF_HOME/assets`에 있다는 것을 이미 확인했습니다.
사용 가능한 설정 파일을 나열합니다:

```bash
ls $NXF_HOME/assets/nf-core/demo/conf/
```

```console
base.config  igenomes.config  igenomes_ignored.config  modules.config  test.config  test_full.config
```

<figure class="excalidraw">
--8<-- "docs/en/docs/hello_nf-core/img/nfcore_config_files.excalidraw.svg"
</figure>

가장 중요한 설정 파일은 다음과 같습니다:

- **`conf/base.config`**: 프로세스에 CPU, 메모리, 시간을 할당하는 리소스 레이블(`process_low`, `process_medium`, `process_high`)을 정의합니다. 프로세스가 예상보다 많은 리소스를 사용하는 경우, 이 파일에서 해당 기본값을 확인할 수 있습니다.
- **`conf/modules.config`**: 프로세스별 도구 인자(`ext.args`)와 출력 게시 설정(`publishDir`)을 지정합니다. 이 파일을 열어 각 도구가 기본적으로 받는 인자를 확인하십시오.
- **`conf/test.config`**: 섹션 2.1에서 사용한 test 프로파일로, `resourceLimits`를 통해 리소스를 제한하고 테스트 샘플시트를 설정합니다. `-profile test`로 활성화됩니다.
  전체 크기 테스트 데이터셋으로 실행하기 위한 `conf/test_full.config`도 있으며, 벤치마킹에 유용합니다.

중앙 `nextflow.config`는 위의 모든 파일을 로드하고 모든 항목에 대한 적절한 기본값을 설정합니다.

이러한 파일에 지정된 설정을 수정하려면 해당 파일을 직접 수정하지 마십시오.
대신 자체 설정 파일을 만들어 `-c`로 전달하십시오.
지정한 값이 다른 파일에 설정된 기본값을 재정의합니다.

실제로 이를 수행하는 몇 가지 연습을 진행해보겠습니다.

#### 3.2.1. 프로세스의 리소스 할당 변경하기

데모 파이프라인은 `base.config`에 정의된 레이블을 사용하여 리소스를 할당합니다.
예를 들어, `FASTQC`는 6개의 CPU와 36GB 메모리를 할당하는 `process_medium` 레이블을 사용합니다.

test 프로파일은 `resourceLimits`를 통해 리소스를 제한하지만, 특정 프로세스의 리소스를 재정의할 수도 있습니다.

`custom.config` 파일을 만드십시오:

```groovy title="custom.config" linenums="1"
process {
    withName: 'FASTQC' {
        cpus = 2
        memory = 4.GB
    }
}
```

사용자 정의 설정으로 파이프라인을 실행합니다:

```bash
nextflow run nf-core/demo -profile docker,test --outdir demo-results-custom -c custom.config
```

??? success "명령 출력"

    ```console
    executor >  local (7)
    [2a/f17b3e] NFCORE_DEMO:DEMO:FASTQC (SAMPLE3_SE)     | 3 of 3 ✔
    [9c/e4d028] NFCORE_DEMO:DEMO:SEQTK_TRIM (SAMPLE3_SE) | 3 of 3 ✔
    [5b/a93c71] NFCORE_DEMO:DEMO:MULTIQC                 | 1 of 1 ✔
    -[nf-core/demo] Pipeline completed successfully-
    ```

`-c` 플래그는 파이프라인의 내장 설정 위에 사용자 정의 설정을 추가합니다.

#### 3.2.2. `ext.args`로 도구 인자 값 설정하기

많은 명령줄 도구에는 필수가 아닌 인자가 있어, 매우 일반적으로 사용되지 않는 한 파이프라인 매개변수로 설정되지 않습니다.
이러한 도구 인자의 경우, nf-core 모듈은 설정 파일을 통해 기본 도구에 인자를 전달하기 위해 `ext.args`라는 Nextflow 규칙을 사용합니다.

예를 들어, `ext.args`를 사용하여 `SEQTK_TRIM` 모듈에 트리밍 인자를 추가해보겠습니다.

##### 3.2.2.1. 사용자 정의 설정 업데이트하기

`custom.config`를 업데이트합니다:

```groovy title="custom.config" linenums="1" hl_lines="6 7 8"
process {
    withName: 'FASTQC' {
        cpus = 2
        memory = 4.GB
    }
    withName: 'SEQTK_TRIM' {
        ext.args = '-b 5'
    }
}
```

이렇게 하면 `seqtk trimfq`가 품질 트리밍 외에 각 리드의 시작 부분에서 5개의 염기를 추가로 트리밍합니다.

##### 3.2.2.2. 파이프라인 실행하기

이 설정으로 파이프라인을 다시 실행하여 효과를 확인합니다:

```bash
nextflow run nf-core/demo -profile docker,test --outdir demo-results-extargs -c custom.config
```

??? success "명령 출력"

    ```console
    executor >  local (7)
    [1e/b7a392] NFCORE_DEMO:DEMO:FASTQC (SAMPLE3_SE)     | 3 of 3 ✔
    [ab/cd1234] NFCORE_DEMO:DEMO:SEQTK_TRIM (SAMPLE3_SE) | 3 of 3 ✔
    [4f/c8d105] NFCORE_DEMO:DEMO:MULTIQC                 | 1 of 1 ✔
    -[nf-core/demo] Pipeline completed successfully-
    ```

인자가 적용되었는지 확인하려면, 실행 출력에서 `SEQTK_TRIM` work 디렉토리 해시(예: `work/ab/cd1234...`)를 찾아 내부의 `.command.sh` 파일을 확인하십시오:

```bash
cat work/ab/cd1234/.command.sh
```

??? success "명령 출력"

    ```console
    #!/usr/bin/env bash
    ...
    seqtk trimfq -b 5 SAMPLE3_SE.fastq.gz | gzip -c > SAMPLE3_SE.trimmed.fastq.gz
    ```

`seqtk trimfq` 명령에 `-b 5`가 포함되어 있어 `ext.args` 재정의가 적용되었음을 확인할 수 있습니다.

##### 3.2.2.3. 기본값 재정의하기

일부 모듈에는 기본적으로 `ext.args`가 이미 설정되어 있습니다.
예를 들어, `FASTQC` 모듈은 기본적으로 `ext.args = '--quiet'`로 설정되어 있습니다(`conf/modules.config`에 정의됨).

```groovy title="conf/modules.config" linenums="21" hl_lines="2"
    withName: FASTQC {
        ext.args = '--quiet'
        publishDir = [
            path: { "${params.outdir}/fastqc/${meta.id}" },
            mode: params.publish_dir_mode,
            pattern: "*.{html,json}"
        ]
    }
```

사용자 정의 설정 파일을 통해 `ext.args` 값을 제공하면, 해당 값이 해당 프로세스에 설정된 기본값을 완전히 대체합니다.

예를 들어, 기본값이 `'--quiet'`이고 `ext.args = '--kmers 8'`로 설정하면 `--quiet` 플래그가 더 이상 적용되지 않습니다.
두 가지를 모두 유지하려면 `ext.args = '--quiet --kmers 8'`로 설정하십시오.

따라서 `ext.args`로 인자 값을 제공하려는 도구의 기본 설정을 직접 확인할 책임이 있습니다.

### 핵심 정리

nf-core 파이프라인에서 도움말을 얻는 방법, 매개변수를 설정하고 검증 방식을 이해하는 방법, 설정 파일을 통해 설정을 맞춤화하는 방법을 알게 되었습니다.

### 다음 단계

잠시 휴식을 취하십시오! 준비가 되었으면 파트 2로 이동하여 처음부터 nf-core 호환 파이프라인을 직접 만들어보겠습니다.
