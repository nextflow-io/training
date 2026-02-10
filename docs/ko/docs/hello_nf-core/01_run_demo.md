# 파트 1: 데모 파이프라인 실행하기

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } AI 지원 번역 - [자세히 알아보기 및 개선 제안](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

Hello nf-core 교육 과정의 첫 번째 파트에서는 nf-core 파이프라인을 찾아서 사용해보고, 코드 구성 방식을 이해하며, [Hello Nextflow](../hello_nextflow/index.md)에서 보여드린 일반 Nextflow 코드와 어떻게 다른지 알아보겠습니다.

nf-core 프로젝트에서 코드 구조와 도구 작동을 시연하기 위한 파이프라인 모음의 일부로 유지 관리하는 nf-core/demo라는 파이프라인을 사용하겠습니다.

[시작하기](./00_orientation.md) 페이지의 안내에 따라 작업 디렉토리를 `hello-nf-core/`로 설정했는지 확인하십시오.

---

## 1. nf-core/demo 파이프라인 찾기 및 가져오기

[nf-co.re](https://nf-co.re) 프로젝트 웹사이트에서 nf-core/demo 파이프라인을 찾는 것부터 시작하겠습니다. 이 웹사이트는 일반 문서 및 도움말 기사, 각 파이프라인에 대한 문서, 블로그 게시물, 이벤트 공지 등 모든 정보를 중앙에서 관리합니다.

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

새로운 파이프라인 도입을 고려할 때는 실행을 시도하기 전에 파이프라인 문서를 주의 깊게 읽어 파이프라인이 무엇을 하는지, 어떻게 구성해야 하는지 이해해야 합니다.

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

예제 명령이 workflow 파일을 지정하지 않고 파이프라인 저장소에 대한 참조인 `nf-core/demo`만 지정한다는 것을 알 수 있습니다.

이런 방식으로 호출되면 Nextflow는 코드가 특정 방식으로 구성되어 있다고 가정합니다.
이 구조를 살펴볼 수 있도록 코드를 가져오겠습니다.

### 1.2. 파이프라인 코드 가져오기

파이프라인이 우리의 목적에 적합해 보인다고 판단했다면 사용해보겠습니다.
다행히 Nextflow는 수동으로 다운로드할 필요 없이 올바르게 형식화된 저장소에서 파이프라인을 쉽게 가져올 수 있게 해줍니다.

터미널로 돌아가서 다음을 실행하겠습니다:

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

Nextflow에게 이런 방식으로 가져온 파이프라인 목록을 제공하도록 할 수 있습니다:

```bash
nextflow list
```

??? success "명령 출력"

    ```console
    nf-core/demo
    ```

파일들이 현재 작업 디렉토리에 없다는 것을 알 수 있습니다.
기본적으로 Nextflow는 `$NXF_HOME/assets`에 저장합니다.

```bash
tree -L 2 $NXF_HOME/assets/
```

```console title="디렉토리 내용"
/workspaces/.nextflow/assets/
└── nf-core
    └── demo

2 directories, 0 files
```

!!! note

    교육 환경을 사용하지 않는 경우 시스템에서 전체 경로가 다를 수 있습니다.

Nextflow는 이러한 파이프라인을 직접 상호작용하는 코드가 아닌 라이브러리처럼 사용해야 한다는 원칙에 따라 다운로드된 소스 코드를 의도적으로 '방해가 되지 않는' 위치에 보관합니다.

그러나 이 교육의 목적상 내부를 살펴보고 무엇이 있는지 확인하고 싶습니다.
따라서 쉽게 확인할 수 있도록 현재 작업 디렉토리에서 해당 위치로의 심볼릭 링크를 만들겠습니다.

```bash
ln -s $NXF_HOME/assets pipelines
```

이렇게 하면 방금 다운로드한 코드를 더 쉽게 탐색할 수 있는 바로가기가 만들어집니다.

```bash
tree -L 2 pipelines
```

```console title="디렉토리 내용"
pipelines
└── nf-core
    └── demo

2 directories, 0 files
```

이제 필요에 따라 소스 코드를 더 쉽게 살펴볼 수 있습니다.

하지만 먼저 첫 번째 nf-core 파이프라인을 실행해보겠습니다!

### 요약

이제 nf-core 웹사이트를 통해 파이프라인을 찾고 소스 코드의 로컬 사본을 가져오는 방법을 알게 되었습니다.

### 다음 단계

최소한의 노력으로 nf-core 파이프라인을 사용해보는 방법을 배워보겠습니다.

---

## 2. 테스트 프로파일로 파이프라인 사용해보기

편리하게도 모든 nf-core 파이프라인에는 test 프로파일이 함께 제공됩니다.
이것은 [nf-core/test-datasets](https://github.com/nf-core/test-datasets) 저장소에서 호스팅되는 작은 테스트 데이터셋을 사용하여 파이프라인을 실행하기 위한 최소한의 구성 설정 모음입니다.
작은 규모로 파이프라인을 빠르게 사용해볼 수 있는 좋은 방법입니다.

!!! note

    Nextflow의 configuration profile 시스템을 사용하면 다양한 컨테이너 엔진이나 실행 환경 간에 쉽게 전환할 수 있습니다.
    자세한 내용은 [Hello Nextflow Part 6: Configuration](../hello_nextflow/06_hello_config.md)을 참조하십시오.

### 2.1. 테스트 프로파일 살펴보기

파이프라인을 실행하기 전에 파이프라인의 test 프로파일이 무엇을 지정하는지 확인하는 것이 좋습니다.
`nf-core/demo`의 `test` 프로파일은 구성 파일 `conf/test.config`에 있으며 아래에 표시되어 있습니다.

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
        cpus: 4,
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

주석 블록 다음 섹션으로 이동하면, test 프로파일이 테스트를 위해 미리 구성된 내용을 보여줍니다: 가장 주목할 점은 `input` 매개변수가 이미 테스트 데이터셋을 가리키도록 설정되어 있으므로 자체 데이터를 제공할 필요가 없다는 것입니다.
미리 구성된 입력 링크를 따라가면 여러 실험 샘플에 대한 샘플 식별자와 파일 경로가 포함된 csv 파일임을 알 수 있습니다.

```csv title="samplesheet_test_illumina_amplicon.csv"
sample,fastq_1,fastq_2
SAMPLE1_PE,https://raw.githubusercontent.com/nf-core/test-datasets/viralrecon/illumina/amplicon/sample1_R1.fastq.gz,https://raw.githubusercontent.com/nf-core/test-datasets/viralrecon/illumina/amplicon/sample1_R2.fastq.gz
SAMPLE2_PE,https://raw.githubusercontent.com/nf-core/test-datasets/viralrecon/illumina/amplicon/sample2_R1.fastq.gz,https://raw.githubusercontent.com/nf-core/test-datasets/viralrecon/illumina/amplicon/sample2_R2.fastq.gz
SAMPLE3_SE,https://raw.githubusercontent.com/nf-core/test-datasets/viralrecon/illumina/amplicon/sample1_R1.fastq.gz,
SAMPLE3_SE,https://raw.githubusercontent.com/nf-core/test-datasets/viralrecon/illumina/amplicon/sample2_R1.fastq.gz,
```

이것은 샘플시트라고 하며, nf-core 파이프라인에 대한 가장 일반적인 입력 형식입니다.

!!! note

    데이터 형식과 유형에 익숙하지 않더라도 걱정하지 마십시오. 이후 내용에 중요하지 않습니다.

따라서 파이프라인을 사용해보는 데 필요한 모든 것이 있음을 확인했습니다.

### 2.2. 파이프라인 실행하기

컨테이너 시스템으로 Docker를, 출력 디렉토리로 `demo-results`를 사용하기로 결정했다면, 테스트 명령을 실행할 준비가 된 것입니다:

```bash
nextflow run nf-core/demo -profile docker,test --outdir demo-results
```

??? success "명령 출력"

    ```console
     N E X T F L O W   ~  version 25.04.3

    Launching `https://github.com/nf-core/demo` [magical_pauling] DSL2 - revision: db7f526ce1 [master]


    ------------------------------------------------------
                                            ,--./,-.
            ___     __   __   __   ___     /,-._.--~'
      |\ | |__  __ /  ` /  \ |__) |__         }  {
      | \| |       \__, \__/ |  \ |___     \`-._,-`-,
                                            `._,._,'
      nf-core/demo 1.0.2
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
파이프라인의 버전, 입력 및 출력, 그리고 몇 가지 구성 요소의 요약이 포함된 헤더가 있습니다.

!!! note

    출력에는 다른 타임스탬프, 실행 이름 및 파일 경로가 표시되지만 전체 구조와 프로세스 실행은 유사해야 합니다.

실행 출력으로 이동하여 어떤 프로세스가 실행되었는지 알려주는 줄을 살펴보겠습니다:

```console
    [ff/a6976b] NFCORE_DEMO:DEMO:FASTQC (SAMPLE3_SE)     | 3 of 3 ✔
    [39/731ab7] NFCORE_DEMO:DEMO:SEQTK_TRIM (SAMPLE3_SE) | 3 of 3 ✔
    [7c/78d96e] NFCORE_DEMO:DEMO:MULTIQC                 | 1 of 1 ✔
```

이것은 nf-core 웹사이트의 파이프라인 문서 페이지에 표시된 세 가지 도구에 해당하는 세 개의 프로세스가 실행되었음을 알려줍니다: FASTQC, SEQTK_TRIM 및 MULTIQC.

여기에 표시된 `NFCORE_DEMO:DEMO:MULTIQC`와 같은 전체 프로세스 이름은 Hello Nextflow 입문 자료에서 본 것보다 깁니다.
여기에는 상위 workflow의 이름이 포함되어 있으며 파이프라인 코드의 모듈성을 반영합니다.
조금 후에 이에 대해 자세히 설명하겠습니다.

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
`nf-core/demo` 파이프라인의 출력에 대해 자세히 알아보려면 [문서 페이지](https://nf-co.re/demo/1.0.2/docs/output/)를 확인하십시오.

현 단계에서 중요한 것은 결과가 모듈별로 구성되어 있고, 파이프라인 실행에 대한 다양한 타임스탬프가 있는 보고서가 포함된 `pipeline_info`라는 디렉토리가 추가로 있다는 것입니다.

예를 들어, `execution_timeline_*` 파일은 어떤 프로세스가 실행되었는지, 어떤 순서로, 얼마나 오래 실행되었는지 보여줍니다:

![실행 타임라인 보고서](./img/execution_timeline.png)

!!! note

    여기서 작업이 병렬로 실행되지 않은 이유는 Github Codespaces의 최소 사양 머신에서 실행하고 있기 때문입니다.
    병렬 실행을 보려면 codespace의 CPU 할당과 테스트 구성의 리소스 제한을 늘려보십시오.

이러한 보고서는 모든 nf-core 파이프라인에 대해 자동으로 생성됩니다.

### 요약

내장된 test 프로파일을 사용하여 nf-core 파이프라인을 실행하는 방법과 출력을 찾을 수 있는 위치를 알게 되었습니다.

### 다음 단계

파이프라인 코드가 어떻게 구성되어 있는지 배워보겠습니다.

---

## 3. 파이프라인 코드 구조 살펴보기

이제 사용자로서 파이프라인을 성공적으로 실행했으니, nf-core 파이프라인이 내부적으로 어떻게 구성되어 있는지 살펴보기 위해 관점을 전환해보겠습니다.

nf-core 프로젝트는 파이프라인 구조, 코드 구성 방법, 구성 및 문서화에 대한 강력한 가이드라인을 적용합니다.
이 모든 것이 어떻게 구성되어 있는지 이해하는 것이 이 과정의 Part 2에서 다룰 nf-core 호환 파이프라인을 개발하기 위한 첫 번째 단계입니다.

앞서 만든 `pipelines` 심볼릭 링크를 사용하여 `nf-core/demo` 저장소에서 파이프라인 코드가 어떻게 구성되어 있는지 살펴보겠습니다.

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

여기에는 많은 내용이 있으므로 단계별로 다루겠습니다.

먼저, 최상위 레벨에서 요약 정보가 있는 README 파일과 라이선싱, 기여 가이드라인, 인용 및 행동 강령과 같은 프로젝트 정보를 요약하는 보조 파일을 찾을 수 있습니다.
상세한 파이프라인 문서는 `docs` 디렉토리에 있습니다.
이 모든 콘텐츠는 nf-core 웹사이트의 웹 페이지를 프로그래밍 방식으로 생성하는 데 사용되므로 항상 코드와 최신 상태를 유지합니다.

이제 나머지 부분에 대해서는 세 단계로 나누어 살펴보겠습니다:

1. 파이프라인 코드 구성요소 (`main.nf`, `workflows`, `subworkflows`, `modules`)
2. 파이프라인 구성
3. 입력 및 검증

파이프라인 코드 구성요소부터 시작하겠습니다.
개별 파일 내의 코드를 자세히 살펴보기보다는 파일 계층 구조와 구조적 구성에 중점을 둘 것입니다.

### 3.1. 파이프라인 코드 구성요소

표준 nf-core 파이프라인 코드 구성은 [Hello Nextflow](../hello_nextflow/index.md) 과정의 Part 4인 [Hello Modules](../hello_nextflow/04_hello_modules.md)에서 소개된 것처럼 코드 재사용을 극대화하도록 설계된 모듈식 구조를 따릅니다. 다만 진정한 nf-core 방식으로 약간의 복잡성이 추가되어 구현됩니다.
특히, nf-core 파이프라인은 subworkflow, 즉 상위 workflow에서 가져오는 workflow 스크립트를 풍부하게 사용합니다.

조금 추상적으로 들릴 수 있으므로 `nf-core/demo` 파이프라인에서 실제로 어떻게 사용되는지 살펴보겠습니다.

!!! note

    이러한 모듈식 구성요소가 _어떻게_ 연결되는지에 대한 실제 코드는 다루지 않을 것입니다. subworkflow 사용과 관련된 복잡성이 혼란스러울 수 있고, 이를 이해하는 것이 교육의 현 단계에서는 필요하지 않기 때문입니다.
    지금은 전체적인 구성과 논리에 초점을 맞추겠습니다.

#### 3.1.1. 일반 개요

다음은 `nf-core/demo` 파이프라인에 대한 관련 코드 구성요소 간의 관계를 나타낸 것입니다:

<figure class="excalidraw">
    --8<-- "docs/en/docs/hello_nf-core/img/nf-core_demo_code_organization.svg"
</figure>

`main.nf`라는 _엔트리포인트_ 스크립트가 있으며, 이는 두 종류의 중첩된 workflow에 대한 래퍼 역할을 합니다: `workflows/` 아래에 있고 `demo.nf`라고 하는 실제 분석 로직을 포함하는 workflow와 `subworkflows/` 아래에 있는 관리용 workflow 집합입니다.
`demo.nf` workflow는 `modules/` 아래에 있는 **모듈**을 호출하며, 이들은 실제 분석 단계를 수행할 **프로세스**를 포함합니다.

!!! note

    Subworkflow는 관리 기능에만 국한되지 않으며 프로세스 모듈을 사용할 수 있습니다.

    여기에 표시된 `nf-core/demo` 파이프라인은 스펙트럼에서 단순한 편에 속하지만, 다른 nf-core 파이프라인(예: `nf-core/rnaseq`)은 실제 분석에 관여하는 subworkflow를 활용합니다.

이제 이러한 구성요소를 차례로 검토하겠습니다.

#### 3.1.2. 엔트리포인트 스크립트: `main.nf`

`main.nf` 스크립트는 `nextflow run nf-core/demo`를 실행할 때 Nextflow가 시작하는 엔트리포인트입니다.
즉, `nextflow run nf-core/demo`를 실행하여 파이프라인을 실행하면 Nextflow가 자동으로 `main.nf` 스크립트를 찾아 실행합니다.
이것은 nf-core 파이프라인뿐만 아니라 이러한 관례적인 명명 및 구조를 따르는 모든 Nextflow 파이프라인에 적용됩니다.

엔트리포인트 스크립트를 사용하면 실제 분석 스크립트가 실행되기 전후에 표준화된 '관리용' subworkflow를 쉽게 실행할 수 있습니다.
실제 분석 workflow와 그 모듈을 검토한 후에 이에 대해 살펴보겠습니다.

#### 3.1.3. 분석 스크립트: `workflows/demo.nf`

`workflows/demo.nf` workflow는 파이프라인의 핵심 로직이 저장된 곳입니다.
상위 workflow에서 호출되도록 설계되어 몇 가지 추가 기능이 필요하다는 점을 제외하면 일반 Nextflow workflow와 매우 유사하게 구성됩니다.
다음 과정 파트에서 Hello Nextflow의 간단한 Hello 파이프라인을 nf-core 호환 형태로 변환할 때 관련 차이점을 다루겠습니다.

`demo.nf` workflow는 다음에 검토할 `modules/` 아래에 있는 **모듈**을 호출합니다.

!!! note

    일부 nf-core 분석 workflow는 하위 레벨 subworkflow를 호출하여 추가 중첩 레벨을 표시합니다.
    이것은 주로 일반적으로 함께 사용되는 두 개 이상의 모듈을 쉽게 재사용 가능한 파이프라인 세그먼트로 적용하는 데 사용됩니다.
    nf-core 웹사이트에서 사용 가능한 [nf-core subworkflows](https://nf-co.re/subworkflows/)를 탐색하여 몇 가지 예제를 볼 수 있습니다.

    분석 스크립트가 subworkflow를 사용하는 경우 `subworkflows/` 디렉토리 아래에 저장됩니다.

#### 3.1.4. 모듈

모듈은 [Hello Nextflow 교육 과정의 Part 4](../hello_nextflow/04_hello_modules.md)에 설명된 대로 프로세스 코드가 있는 곳입니다.

nf-core 프로젝트에서 모듈은 출처와 내용을 모두 반영하는 다단계 중첩 구조를 사용하여 구성됩니다.
최상위 레벨에서 모듈은 `nf-core` 또는 `local`(nf-core 프로젝트의 일부가 아님)로 구분되며, 그런 다음 적용하는 도구의 이름을 따서 명명된 디렉토리에 배치됩니다.
도구가 툴킷(즉, 여러 도구를 포함하는 패키지)에 속하는 경우 툴킷 이름을 따서 명명된 중간 디렉토리 레벨이 있습니다.

`nf-core/demo` 파이프라인 모듈에 실제로 적용된 것을 볼 수 있습니다:

```bash
tree -L 3 pipelines/nf-core/demo/modules
```

??? abstract "디렉토리 내용"

    ```console
    pipelines/nf-core/demo/modules
    └── nf-core
        ├── fastqc
        │   ├── environment.yml
        │   ├── main.nf
        │   ├── meta.yml
        │   └── tests
        ├── multiqc
        │   ├── environment.yml
        │   ├── main.nf
        │   ├── meta.yml
        │   └── tests
        └── seqtk
            └── trim

    7 directories, 6 files
    ```

여기서 `fastqc` 및 `multiqc` 모듈은 `nf-core` 모듈 내 최상위 레벨에 있는 반면, `trim` 모듈은 속해 있는 툴킷인 `seqtk` 아래에 있습니다.
이 경우 `local` 모듈은 없습니다.

프로세스를 설명하는 모듈 코드 파일은 항상 `main.nf`라고 하며, 지금은 무시할 테스트 및 `.yml` 파일이 함께 제공됩니다.

엔트리포인트 workflow, 분석 workflow 및 모듈을 종합하면 파이프라인의 '흥미로운' 부분을 실행하기에 충분합니다.
그러나 관리용 subworkflow도 있다는 것을 알고 있으므로 이제 살펴보겠습니다.

#### 3.1.5. 관리용 subworkflow

모듈과 마찬가지로 subworkflow는 `local` 및 `nf-core` 디렉토리로 구분되며, 각 subworkflow에는 자체 `main.nf` 스크립트, 테스트 및 `.yml` 파일이 있는 고유한 중첩 디렉토리 구조가 있습니다.

```bash
tree -L 3 pipelines/nf-core/demo/subworkflows
```

??? abstract "디렉토리 내용"

    ```console
    pipelines/nf-core/demo/subworkflows
    ├── local
    │   └── utils_nfcore_demo_pipeline
    │       └── main.nf
    └── nf-core
        ├── utils_nextflow_pipeline
        │   ├── main.nf
        │   ├── meta.yml
        │   └── tests
        ├── utils_nfcore_pipeline
        │   ├── main.nf
        │   ├── meta.yml
        │   └── tests
        └── utils_nfschema_plugin
            ├── main.nf
            ├── meta.yml
            └── tests

    9 directories, 7 files
    ```

위에서 언급했듯이 `nf-core/demo` 파이프라인에는 분석 관련 subworkflow가 포함되어 있지 않으므로 여기에 표시된 모든 subworkflow는 이름의 `utils_` 접두사로 표시된 것처럼 소위 '관리용' 또는 '유틸리티' workflow입니다.
이러한 subworkflow는 다른 보조 기능 중에서도 콘솔 출력에 멋진 nf-core 헤더를 생성하는 것입니다.

!!! tip

    명명 패턴 외에도 이러한 subworkflow가 진정한 분석 관련 기능을 수행하지 않는다는 또 다른 표시는 프로세스를 전혀 호출하지 않는다는 것입니다.

이것으로 `nf-core/demo` 파이프라인을 구성하는 핵심 코드 구성요소의 정리가 완료되었습니다.
이제 개발에 뛰어들기 전에 조금 알아야 할 나머지 요소인 파이프라인 구성과 입력 검증을 살펴보겠습니다.

### 3.2. 파이프라인 구성

이전에 Nextflow가 입력 및 매개변수, 컴퓨팅 리소스 및 오케스트레이션의 기타 측면과 관련하여 파이프라인 실행을 구성하기 위한 많은 옵션을 제공한다는 것을 배웠습니다.
nf-core 프로젝트는 파이프라인 간 일관성과 유지 관리성을 제공하는 방식으로 Nextflow의 유연한 사용자 정의 옵션을 기반으로 구축하는 것을 목표로 하는 파이프라인 구성에 대한 고도로 표준화된 가이드라인을 적용합니다.

중앙 구성 파일 `nextflow.config`는 매개변수 및 기타 구성 옵션의 기본값을 설정하는 데 사용됩니다.
이러한 구성 옵션의 대부분은 기본적으로 적용되는 반면 다른 옵션(예: 소프트웨어 의존성 프로파일)은 선택적 프로파일로 포함됩니다.

`conf` 폴더에 저장된 여러 추가 구성 파일이 있으며, 기본적으로 또는 선택적으로 프로파일로 구성에 추가할 수 있습니다:

- `base.config`: 대부분의 고성능 컴퓨팅 환경에서 일반적으로 사용하기에 적합한 '백지 상태' 구성 파일입니다. 예를 들어 모듈에 적용하기 편리한 광범위한 리소스 사용 구간을 정의합니다.
- `modules.config`: 추가 모듈 지시문 및 인수입니다.
- `test.config`: 데모 파이프라인을 실행할 때 사용한 최소 테스트 데이터로 파이프라인을 실행하기 위한 프로파일입니다.
- `test_full.config`: 전체 크기 테스트 데이터셋으로 파이프라인을 실행하기 위한 프로파일입니다.

과정 후반부에서 이러한 파일 중 몇 가지를 다루겠습니다.

### 3.3. 입력 및 검증

앞서 `nf-core/demo` 파이프라인의 test 프로파일을 살펴보았을 때 언급했듯이, 파일 경로와 샘플 식별자가 포함된 샘플시트를 입력으로 받도록 설계되었습니다.
파일 경로는 `nf-core/test-datasets` 저장소에 있는 실제 데이터에 연결되어 있습니다.

샘플시트 예제도 `assets` 디렉토리에 제공되지만 이 파일의 경로는 실제가 아닙니다.

```csv title="assets/samplesheet.csv" linenums="1"
sample,fastq_1,fastq_2
SAMPLE_PAIRED_END,/path/to/fastq/files/AEG588A1_S1_L002_R1_001.fastq.gz,/path/to/fastq/files/AEG588A1_S1_L002_R2_001.fastq.gz
SAMPLE_SINGLE_END,/path/to/fastq/files/AEG588A4_S4_L003_R1_001.fastq.gz,

```

이 특정 샘플시트는 상당히 간단하지만, 일부 파이프라인은 기본 입력과 관련된 더 많은 메타데이터가 있는 더 복잡한 샘플시트에서 실행됩니다.

안타깝게도 이러한 파일은 눈으로 확인하기 어려울 수 있으므로 입력 데이터의 부적절한 형식 지정은 파이프라인 실패의 매우 일반적인 원인입니다.
관련 문제는 매개변수가 잘못 제공되는 경우입니다.

이러한 문제에 대한 해결책은 모든 입력 파일에 대해 예상되는 정보 유형이 올바르게 형식화되어 포함되어 있는지 확인하기 위해 자동 검증 확인을 실행하고, 매개변수가 예상 유형인지 확인하는 것입니다.
이것을 입력 검증이라고 하며, 파이프라인이 실패하여 입력에 문제가 있음을 알아낼 때까지 기다리기보다는 파이프라인을 실행하기 _전에_ 이상적으로 수행되어야 합니다.

구성과 마찬가지로 nf-core 프로젝트는 입력 검증에 대해 매우 의견이 강하며, Nextflow 파이프라인에 포괄적인 검증 기능을 제공하는 Nextflow 플러그인인 [nf-schema plugin](https://nextflow-io.github.io/nf-schema/latest/)의 사용을 권장합니다.

이 주제는 이 과정의 Part 5에서 더 자세히 다루겠습니다.
지금은 해당 목적으로 `nextflow_schema.json` 및 `assets/schema_input.json`이라는 두 개의 JSON 파일이 제공된다는 점만 알아두십시오.

`nextflow_schema.json`은 유형, 설명 및 도움말 텍스트를 포함한 파이프라인 매개변수에 대한 정보를 기계 판독 가능한 형식으로 저장하는 데 사용되는 파일입니다.
이것은 자동 매개변수 검증, 도움말 텍스트 생성, UI 인터페이스의 대화형 매개변수 양식 렌더링을 포함한 다양한 목적으로 사용됩니다.

`schema_input.json`은 입력 샘플시트 구조를 정의하는 데 사용되는 파일입니다.
각 열은 기계 판독 가능한 형식으로 유형, 패턴, 설명 및 도움말 텍스트를 가질 수 있습니다.
스키마는 자동 검증 및 유용한 오류 메시지 제공을 포함한 다양한 목적으로 사용됩니다.

### 요약

이제 nf-core 파이프라인의 주요 구성요소와 코드 구성 방법, 주요 구성 요소가 어디에 있는지, 입력 검증이 무엇을 위한 것인지 알게 되었습니다.

### 다음 단계

휴식을 취하십시오! 많은 내용이었습니다. 기분이 상쾌하고 준비가 되었을 때 다음 섹션으로 이동하여 배운 내용을 적용하여 nf-core 호환 파이프라인을 작성하십시오.

!!! tip

    다음 파트로 이동하기 전에 subworkflow를 사용하여 workflow를 구성하는 방법을 배우고 싶다면 [Workflows of Workflows](../side_quests/workflows_of_workflows.md) Side Quest를 확인하십시오.
