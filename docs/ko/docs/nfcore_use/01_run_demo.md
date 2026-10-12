# 파트 1: 데모 파이프라인 실행

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } AI 지원 번역 - [자세히 알아보기 및 개선 제안](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

Use nf-core 교육 과정의 첫 번째 파트에서는 nf-core 파이프라인을 찾고, 내장된 테스트 프로파일을 사용하여 실행하는 방법을 학습합니다.

여기서는 nf-core 프로젝트가 데모 및 교육 목적으로 관리하는 nf-core/demo 파이프라인을 사용합니다.

[시작하기](./00_orientation.md) 페이지의 안내에 따라 작업 디렉토리가 `nfcore-use/`로 설정되어 있는지 확인하세요.

---

## 1. nf-core/demo 파이프라인 찾기 및 가져오기

먼저 [nf-co.re](https://nf-co.re) 프로젝트 웹사이트에서 nf-core/demo 파이프라인을 찾습니다. 이 웹사이트는 일반 문서 및 도움말, 각 파이프라인 문서, 블로그 게시물, 이벤트 공지 등 모든 정보를 한곳에서 제공합니다.

### 1.1. 웹사이트에서 파이프라인 찾기

웹 브라우저에서 [https://nf-co.re/pipelines/](https://nf-co.re/pipelines/)로 이동하여 검색창에 `demo`를 입력합니다.

![검색 결과](./img/search-results.png)

파이프라인 이름 `demo`를 클릭하면 파이프라인 문서 페이지로 이동합니다.

릴리스된 각 파이프라인에는 다음과 같은 문서 섹션이 포함된 전용 페이지가 있습니다:

- **Introduction:** 파이프라인 소개 및 개요
- **Usage:** 파이프라인 실행 방법 설명
- **Parameters:** 설명이 포함된 파이프라인 매개변수 그룹
- **Output:** 예상 출력 파일에 대한 설명 및 예제
- **Results:** 전체 테스트 데이터셋에서 생성된 출력 파일 예제
- **Releases & Statistics:** 파이프라인 버전 이력 및 통계

새로운 파이프라인 도입을 고려할 때는 실행하기 전에 파이프라인 문서를 꼼꼼히 읽어 파이프라인이 무엇을 하는지, 어떻게 설정해야 하는지 이해해야 합니다.

지금 살펴보고 다음 내용을 확인해 보세요:

- 파이프라인이 실행할 도구 (`Introduction` 탭 확인)
- 파이프라인이 허용하거나 요구하는 입력 및 매개변수 (`Parameters` 탭 확인)
- 파이프라인이 생성하는 출력 (`Output` 탭 확인)

#### 1.1.1. 파이프라인 개요

`Introduction` 탭은 파이프라인의 개요를 제공하며, 시각적 표현(서브웨이 맵이라고 함)과 파이프라인의 일부로 실행되는 도구 목록이 포함되어 있습니다.

![파이프라인 서브웨이 맵](./img/nf-core-demo-subway-cropped.png)

1. 리드 QC ([FASTQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/))
2. 어댑터 및 품질 트리밍 ([SEQTK_TRIM](https://github.com/lh3/seqtk))
3. 원시 리드에 대한 QC 결과 표시 ([MULTIQC](http://multiqc.info/))
4. 소에서 재미있는 텍스트 메시지 생성 ([COWPY](https://github.com/jeffbuttars/cowpy))

#### 1.1.2. 예제 명령줄

문서에는 예제 입력 파일(아래에서 자세히 설명)과 예제 명령줄도 제공됩니다.

```bash
nextflow run nf-core/demo \
  -profile <docker/singularity/.../institute> \
  --input samplesheet.csv \
  --outdir <OUTDIR>
```

예제 명령에는 워크플로우 파일이 지정되지 않고, 파이프라인 저장소 참조인 `nf-core/demo`만 지정되어 있습니다.

이 방식으로 실행하면 Nextflow는 코드가 특정 방식으로 구성되어 있다고 가정합니다.
이 구조를 살펴보기 위해 코드를 가져옵니다.

### 1.2. 파이프라인 코드 가져오기

파이프라인이 목적에 적합하다고 판단되면 실행해 봅니다.
Nextflow는 올바른 형식의 저장소에서 파이프라인을 수동으로 다운로드하지 않고도 쉽게 가져올 수 있습니다.

#### 1.2.1. `nextflow pull` 사용

터미널로 돌아가서 다음 명령을 실행합니다:

```bash
nextflow pull nf-core/demo
```

??? success "명령 출력"

    ```console
    Checking nf-core/demo ...
     downloaded from https://github.com/nf-core/demo.git - revision: 32893afef8 [master]
    ```

Nextflow는 파이프라인 코드를 `pull`하여 전체 저장소를 로컬 드라이브에 다운로드합니다.

이 방법은 nf-core 파이프라인뿐만 아니라 GitHub에 적절히 설정된 모든 Nextflow 파이프라인에 사용할 수 있습니다.
다만 nf-core는 가장 큰 오픈소스 Nextflow 파이프라인 컬렉션입니다.

#### 1.2.2. `nextflow list` 사용

이 방식으로 가져온 파이프라인 목록을 확인할 수 있습니다:

```bash
nextflow list
```

??? success "명령 출력"

    ```console
    nf-core/demo
    ```

다른 파이프라인을 몇 개 더 pull하여 여러 개가 있을 때 어떻게 나열되는지 확인해 보세요.

#### 1.2.3. 파이프라인이 다운로드된 위치 확인

파일이 현재 작업 디렉토리에 없다는 것을 알 수 있습니다.
기본적으로 Nextflow는 pull한 파이프라인을 `$NXF_HOME/assets` 아래에 저장합니다.

특정 파이프라인의 위치를 확인하려면 Nextflow에 직접 물어보세요:

```bash
nextflow info nf-core/demo
```

??? success "명령 출력"

    ```console
     project name: nf-core/demo
     repository  : https://github.com/nf-core/demo
     local path  : /workspaces/.nextflow/assets/.repos/nf-core/demo
     main script : main.nf
     description : An nf-core demo pipeline
     revisions   :
       TEMPLATE
       bumper
       dev
       fix-nxfversion
       manually-merge-3_0_2
     > master (default)
       nf-core-template-merge-2.13.2.dev0
       nf-core-template-merge-2.14.0
       nf-core-template-merge-2.14.1
       nf-core-template-merge-3.0.0
       nf-core-template-merge-3.0.1
       nf-core-template-merge-3.0.2
       nf-core-template-merge-3.1.0
       nf-core-template-merge-3.1.2
       nf-core-template-merge-3.2.0
       nf-core-template-merge-3.2.1
       nf-core-template-merge-3.3.1
       nf-core-template-merge-3.3.2
       nf-core-template-merge-4.0.0
       nf-core-template-merge-4.0.3
       nf-core-template-merge-4.1.0
       nf-core-template-merge-4.1.0-2
       patch
       1.0.0 [t]
       1.0.1 [t]
       1.0.2 [t]
       1.1.0 [t]
     > 1.2.0 [t]
    ```

!!! info "정보"

    교육 환경을 사용하지 않는 경우 전체 경로가 다를 수 있습니다.

Nextflow는 다운로드된 소스 코드를 의도적으로 '눈에 띄지 않는 곳'에 보관합니다. 이는 이러한 파이프라인이 직접 상호작용하는 코드보다는 라이브러리처럼 사용되어야 한다는 원칙에 따른 것입니다.

내부적으로 Nextflow는 pull한 각 파이프라인을 `$NXF_HOME/assets/.repos/` 아래에 git 저장소로 저장하고, 각 리비전의 코드를 `clones/<commit>/` 하위 디렉토리에 체크아웃합니다.
`.repos`는 숨김 디렉토리이므로 `tree -L 2 $NXF_HOME/assets/`를 실행하면 비어 있는 것처럼 보입니다.

#### 1.2.4. 소스 코드에 쉽게 접근하기 위한 심볼릭 링크 생성

코드를 자세히 살펴보지는 않겠지만, 전체적인 구성이 어떻게 되어 있는지 간략히 확인해 봅니다.

파이프라인 소스 코드를 쉽게 탐색할 수 있도록 체크아웃된 파이프라인 복사본을 가리키는 심볼릭 링크를 생성합니다:

```bash
mkdir -p pipelines/nf-core
ln -s "$(echo $NXF_HOME/assets/.repos/nf-core/demo/clones/*/)" pipelines/nf-core/demo
```

이렇게 하면 `tree -L 2 pipelines/nf-core/demo`로 코드를 탐색하거나 파일을 직접 열 수 있는 바로가기가 생성됩니다.

#### 1.2.5. 코드 구성 개요

`tree`를 사용하거나 파일 탐색기(File Explorer)를 사용하여 `nf-core/demo` 디렉토리를 찾아 열 수 있습니다.

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

    7 directories, 12 files
    ```

보시다시피 많은 내용이 있지만, 대부분은 신경 쓰지 않아도 됩니다.

간략히 살펴보면, 최상위 레벨에는 요약 정보가 담긴 README 파일과 라이선스, 기여 가이드라인, 인용 및 행동 강령 등 프로젝트 정보를 요약한 부속 파일들이 있습니다.
자세한 파이프라인 문서는 `docs` 디렉토리에 있습니다.
이 모든 내용은 nf-core 웹사이트의 웹 페이지를 프로그래밍 방식으로 생성하는 데 사용되므로 항상 코드와 최신 상태를 유지합니다.

나머지 코드 파일은 세 가지 기능 그룹으로 구분할 수 있습니다:

1. 파이프라인 코드 구성 요소 (`main.nf`, `workflows`, `subworkflows`, `modules`)
2. 파이프라인 설정
3. 파이프라인 매개변수 / 입력 및 유효성 검사

이 파트에서는 파이프라인 코드 구성 요소를 다루지 않지만, nf-core 파이프라인의 최종 사용자에게 관련성이 높은 설정 및 유효성 검사 요소는 살펴볼 것입니다.

!!! tip "팁"

    nf-core 파이프라인의 소스 코드는 GitHub에서도 탐색할 수 있습니다. 예: [github.com/nf-core/demo](https://github.com/nf-core/demo).
    모든 nf-core 파이프라인은 동일한 디렉토리 구조를 따르므로, 구조를 한 번 파악하면 어떤 파이프라인에서도 같은 방식으로 설정 파일, 모듈, 워크플로우를 찾을 수 있습니다.

이제 파이프라인을 실행해 봅니다!

### 핵심 정리

nf-core 웹사이트를 통해 파이프라인을 찾고 소스 코드의 로컬 복사본을 가져오는 방법을 학습했습니다.

### 다음 단계

최소한의 노력으로 nf-core 파이프라인을 실행하는 방법을 학습합니다.

---

## 2. 테스트 프로파일로 파이프라인 실행해 보기

편리하게도 모든 nf-core 파이프라인에는 테스트 프로파일이 포함되어 있습니다.
이는 [nf-core/test-datasets](https://github.com/nf-core/test-datasets) 저장소에 호스팅된 소규모 테스트 데이터셋을 사용하여 파이프라인을 실행하기 위한 최소한의 설정 모음입니다.
소규모로 파이프라인을 빠르게 실행해 볼 수 있는 좋은 방법입니다.

!!! tip "팁"

    Nextflow의 설정 프로파일 시스템을 사용하면 다양한 컨테이너 엔진이나 실행 환경 간에 쉽게 전환할 수 있습니다.
    자세한 내용은 [Hello Nextflow 파트 6: 설정](../hello_nextflow/06_hello_config.md)을 참조하세요.

### 2.1. 테스트 프로파일 살펴보기

파이프라인을 실행하기 전에 테스트 프로파일이 무엇을 지정하는지 확인하는 것이 좋습니다.
`nf-core/demo`의 `test` 프로파일은 설정 파일 `conf/test.config`에 있습니다.
섹션 1.2.4에서 생성한 `pipelines` 심볼릭 링크를 통해 `nextflow pull`이 다운로드한 파이프라인 소스 내에서 로컬로 찾을 수 있습니다:

```bash
code pipelines/nf-core/demo/conf/test.config
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
        time: '1.h',
    ]
}

params {
    config_profile_name        = 'Test profile'
    config_profile_description = 'Minimal test dataset to check pipeline function'

    // 입력 데이터
    input                      = 'https://raw.githubusercontent.com/nf-core/test-datasets/viralrecon/samplesheet/samplesheet_test_illumina_amplicon.csv'
}
```

파일 상단의 주석 블록에 이 테스트 프로파일로 파이프라인을 실행하는 방법을 보여주는 사용 예제가 포함되어 있습니다.

```groovy title="conf/test.config" linenums="7"
    Use as follows:
        nextflow run nf-core/demo -profile test,<docker/singularity> --outdir <OUTDIR>
```

예제 명령에서 꺾쇠 괄호 사이에 표시된 내용만 지정하면 됩니다: `<docker/singularity>`와 `<OUTDIR>`.

`<docker/singularity>`는 컨테이너 시스템 선택을 의미합니다. 모든 nf-core 파이프라인은 재현성을 보장하고 소프트웨어 설치 문제를 없애기 위해 컨테이너(Docker, Singularity 등)와 함께 사용할 수 있도록 설계되어 있습니다.
따라서 파이프라인 테스트에 Docker와 Singularity 중 어느 것을 사용할지 지정해야 합니다.

`--outdir <OUTDIR>` 부분은 Nextflow가 파이프라인 출력을 저장할 디렉토리를 의미합니다.
원하는 이름을 지정하면 됩니다.
해당 디렉토리가 존재하지 않으면 Nextflow가 실행 시 자동으로 생성합니다.

주석 블록 다음 섹션을 보면 테스트를 위해 사전 설정된 내용을 확인할 수 있습니다. 특히 `input` 매개변수가 이미 테스트 데이터셋을 가리키도록 설정되어 있으므로 직접 데이터를 제공할 필요가 없습니다.
사전 설정된 입력 링크를 따라가면 여러 실험 샘플의 샘플 식별자와 파일 경로가 포함된 CSV 파일을 확인할 수 있습니다.

```csv title="samplesheet_test_illumina_amplicon.csv"
sample,fastq_1,fastq_2
SAMPLE1_PE,https://raw.githubusercontent.com/nf-core/test-datasets/viralrecon/illumina/amplicon/sample1_R1.fastq.gz,https://raw.githubusercontent.com/nf-core/test-datasets/viralrecon/illumina/amplicon/sample1_R2.fastq.gz
SAMPLE2_PE,https://raw.githubusercontent.com/nf-core/test-datasets/viralrecon/illumina/amplicon/sample2_R1.fastq.gz,https://raw.githubusercontent.com/nf-core/test-datasets/viralrecon/illumina/amplicon/sample2_R2.fastq.gz
SAMPLE3_SE,https://raw.githubusercontent.com/nf-core/test-datasets/viralrecon/illumina/amplicon/sample1_R1.fastq.gz,
SAMPLE3_SE,https://raw.githubusercontent.com/nf-core/test-datasets/viralrecon/illumina/amplicon/sample2_R1.fastq.gz,
```

이를 샘플시트(samplesheet)라고 하며, nf-core 파이프라인에서 가장 일반적인 입력 형식입니다.
데이터 형식과 유형에 익숙하지 않아도 괜찮습니다. 이후 내용에서는 중요하지 않습니다.

이제 파이프라인을 실행할 준비가 되었습니다.

### 2.2. 파이프라인 실행

위에서 언급했듯이 예제 테스트 명령을 거의 그대로 사용할 수 있습니다. 사용할 소프트웨어 패키징과 출력 디렉토리 이름만 지정하면 됩니다.
여기서는 컨테이너 시스템으로 Docker를, 출력 디렉토리 이름으로 `demo-results`를 사용합니다.

테스트 명령을 실행합니다:

```bash
nextflow run nf-core/demo -profile test,docker --outdir demo-results
```

??? success "명령 출력"

    ```console
     N E X T F L O W   ~  version 26.04.4

    Downloading plugin nf-schema@2.7.2
    Launching `https://github.com/nf-core/demo` [cranky_curry] revision: 32893afef8 [master]


    ------------------------------------------------------
                                            ,--./,-.
            ___     __   __   __   ___     /,-._.--~'
      |\ | |__  __ /  ` /  \ |__) |__         }  {
      | \| |       \__, \__/ |  \ |___     \`-._,-`-,
                                            `._,._,'
      nf-core/demo 1.2.0
    ------------------------------------------------------

    Input/output options
      input                     : https://raw.githubusercontent.com/nf-core/test-datasets/viralrecon/samplesheet/samplesheet_test_illumina_amplicon.csv
      outdir                    : demo-results

    Institutional config options
      config_profile_name       : Test profile
      config_profile_description: Minimal test dataset to check pipeline function

    Generic options
      trace_report_suffix       : 2026-07-03_21-31-35

    Core Nextflow options
      revision                  : master
      runName                   : cranky_curry
      containerEngine           : docker
      launchDir                 : /workspaces/training/nfcore-use
      workDir                   : /workspaces/training/nfcore-use/work
      projectDir                : /workspaces/.nextflow/assets/.repos/nf-core/demo/clones/32893afef8076a03a2767a020b3f0cab2e0b40b2
      userName                  : root
      profile                   : test,docker
      configFiles               : /workspaces/.nextflow/assets/.repos/nf-core/demo/clones/32893afef8076a03a2767a020b3f0cab2e0b40b2/nextflow.config

    !! Only displaying parameters that differ from the pipeline defaults !!
    ------------------------------------------------------

    * The pipeline
        https://doi.org/10.5281/zenodo.12192442

    * The nf-core framework
        https://doi.org/10.1038/s41587-020-0439-x

    * Software dependencies
        https://github.com/nf-core/demo/blob/master/CITATIONS.md

    executor >  local (8)
    [ca/5b0f3e] NFCORE_DEMO:DEMO:FASTQC (SAMPLE3_SE)     | 3 of 3 ✔
    [b7/cb6812] NFCORE_DEMO:DEMO:SEQTK_TRIM (SAMPLE3_SE) | 3 of 3 ✔
    [ff/6ebd98] NFCORE_DEMO:DEMO:COWPY                   | 1 of 1 ✔
    [09/bbd1b4] NFCORE_DEMO:DEMO:MULTIQC (demo)          | 1 of 1 ✔
    -[nf-core/demo] Pipeline completed successfully-
    ```

출력이 위와 같다면 축하합니다! 첫 번째 nf-core 파이프라인을 성공적으로 실행했습니다.

기본 Nextflow 파이프라인을 실행할 때보다 콘솔 출력이 훨씬 많다는 것을 알 수 있습니다.
파이프라인 버전, 입력 및 출력 요약, 일부 설정 요소가 포함된 헤더가 있습니다.

!!! info "정보"

    출력에는 다른 타임스탬프, 실행 이름, 파일 경로가 표시되지만 전체 구조와 프로세스 실행은 유사합니다.

출력 상단 근처의 다음 줄을 확인하세요:

```console
Launching `https://github.com/nf-core/demo` [cranky_curry] revision: 32893afef8 [master]
```

이 줄은 사용된 파이프라인의 리비전을 알려줍니다.
버전을 지정하지 않았으므로 Nextflow는 `master`의 최신 커밋을 사용했습니다.
재현 가능한 실행을 위해서는 `-r` 플래그를 사용하여 특정 릴리스를 고정해야 합니다:

```bash
nextflow run nf-core/demo -r 1.2.0 -profile test,docker --outdir demo-results
```

이렇게 하면 새로운 커밋이나 릴리스에 관계없이 항상 동일한 파이프라인 코드가 사용됩니다.
이 교육에서는 간단하게 `-r`을 생략하지만, 실제 운영 환경에서는 항상 지정해야 합니다.

실행 출력으로 넘어가서 실행된 프로세스를 알려주는 줄을 살펴봅니다:

```console
executor >  local (8)
[ca/5b0f3e] NFCORE_DEMO:DEMO:FASTQC (SAMPLE3_SE)     | 3 of 3 ✔
[b7/cb6812] NFCORE_DEMO:DEMO:SEQTK_TRIM (SAMPLE3_SE) | 3 of 3 ✔
[ff/6ebd98] NFCORE_DEMO:DEMO:COWPY                   | 1 of 1 ✔
[09/bbd1b4] NFCORE_DEMO:DEMO:MULTIQC (demo)          | 1 of 1 ✔
-[nf-core/demo] Pipeline completed successfully-
```

nf-core 웹사이트의 파이프라인 문서 페이지에 표시된 네 가지 도구인 `FASTQC`, `SEQTK_TRIM`, `MULTIQC`, `COWPY`에 해당하는 네 개의 프로세스가 실행되었습니다.

`NFCORE_DEMO:DEMO:MULTIQC`와 같이 여기에 표시된 전체 프로세스 이름은 Hello Nextflow 입문 자료에서 보았던 것보다 깁니다.
이는 상위 워크플로우의 이름을 포함하며 파이프라인 코드의 모듈성을 반영합니다.
nf-core 스타일의 파이프라인을 직접 개발하고 싶다면 [Build with nf-core](../nfcore_build/index.md) 과정을 참조하세요.

### 2.3. 파이프라인 출력 살펴보기

마지막으로 파이프라인이 생성한 `demo-results` 디렉토리를 살펴봅니다.

```bash
tree -L 2 demo-results
```

??? abstract "디렉토리 내용"

    ```console
    demo-results
    ├── cowpy
    │   └── cowpy.txt
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
    │   └── multiqc_report.html
    └── pipeline_info
        ├── execution_report_2026-07-03_21-31-35.html
        ├── execution_timeline_2026-07-03_21-31-35.html
        ├── execution_trace_2026-07-03_21-31-35.txt
        ├── nf_core_demo_software_mqc_versions.yml
        ├── params_2026-07-03_21-31-43.json
        └── pipeline_dag_2026-07-03_21-31-35.html

    12 directories, 8 files
    ```

꽤 많아 보일 수 있습니다.
`nf-core/demo` 파이프라인의 출력에 대해 자세히 알아보려면 [문서 페이지](https://nf-co.re/demo/1.2.0/docs/output/)를 확인하세요.

현재 단계에서 중요한 점은 결과가 모듈별로 구성되어 있으며, 파이프라인 실행에 관한 다양한 타임스탬프가 포함된 보고서가 담긴 `pipeline_info` 디렉토리가 추가로 있다는 것입니다.

예를 들어 `execution_timeline_*` 파일은 실행된 프로세스, 실행 순서, 소요 시간을 보여줍니다:

![실행 타임라인 보고서](./img/execution_timeline.png)

!!! info "정보"

    여기서는 Github Codespaces의 최소 사양 머신에서 실행하고 있기 때문에 작업이 병렬로 실행되지 않았습니다.
    병렬 실행을 확인하려면 코드스페이스의 CPU 할당과 테스트 설정의 리소스 제한을 늘려 보세요.

이러한 보고서는 모든 nf-core 파이프라인에서 자동으로 생성됩니다.

### 핵심 정리

내장된 테스트 프로파일을 사용하여 nf-core 파이프라인을 실행하는 방법과 출력 위치를 확인하는 방법을 학습했습니다.

### 다음 단계

[파트 2](./02_configure_execution.md)로 이동하여 파이프라인 실행 설정 방법을 학습합니다.

---

## 요약

이 파트에서 학습한 내용:

- nf-core 파이프라인을 찾고 가져와서 코드 구조 살펴보기
- 내장된 테스트 프로파일을 사용하여 파이프라인 실행하기
