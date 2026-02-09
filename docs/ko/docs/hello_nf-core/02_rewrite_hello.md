# 2부: nf-core용 Hello 재작성

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } AI 지원 번역 - [자세히 알아보기 및 개선 제안](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

Hello nf-core 교육 과정의 두 번째 부분에서는 [Hello Nextflow](../hello_nextflow/index.md) 초급자 과정에서 만든 파이프라인의 nf-core 호환 버전을 생성하는 방법을 보여드립니다.

교육의 첫 번째 섹션에서 nf-core 파이프라인이 많은 부속 파일과 함께 상당히 정교한 구조를 따른다는 것을 알아차리셨을 것입니다.
이 모든 것을 처음부터 만드는 것은 매우 지루할 것이므로, nf-core 커뮤니티는 프로세스를 시작하기 위해 템플릿에서 이를 수행하는 도구를 개발했습니다.

파이프라인 스캐폴드를 생성하기 위해 이 도구를 사용하는 방법을 보여드린 다음, 기존의 '일반' 파이프라인 코드를 nf-core 스캐폴드에 맞게 조정하겠습니다.

Hello 파이프라인에 익숙하지 않거나 복습이 필요하시면 [이 정보 페이지](../info/hello_pipeline.md)를 참조하세요.

---

## 1. 새 파이프라인 프로젝트 생성

먼저 새 파이프라인의 스캐폴드를 생성합니다.

!!! note "참고"

    터미널에서 `hello-nf-core` 디렉토리에 있는지 확인하세요.

### 1.1. 템플릿 기반 파이프라인 생성 도구 실행

`nf-core pipelines create` 명령으로 새 파이프라인을 생성하는 것부터 시작하겠습니다.
이는 파이프라인 이름, 설명, 작성자로 커스터마이징된 nf-core 기본 템플릿을 사용하여 새 파이프라인 스캐폴드를 생성합니다.

```bash
nf-core pipelines create
```

이 명령을 실행하면 파이프라인 생성을 위한 텍스트 사용자 인터페이스(TUI)가 열립니다:

<div style="text-align: center;">
    <iframe width="560" height="315" src="https://www.youtube.com/embed/VwjXNXONHlY?si=d0HkFSISnKn76TeI" title="YouTube video player" frameborder="0" allow="accelerometer; autoplay; clipboard-write; encrypted-media; gyroscope; picture-in-picture; web-share" allowfullscreen="" data-ruffle-polyfilled=""></iframe>
</div>

이 TUI는 파이프라인에 대한 기본 정보를 제공하도록 요청하며 파이프라인 스캐폴드에 포함하거나 제외할 기능을 선택할 수 있게 해줍니다.

- 시작 화면에서 **Let's go!**를 클릭합니다.
- `Choose pipeline type` 화면에서 **Custom**을 클릭합니다.
- 다음과 같이 파이프라인 세부 정보를 입력한 후(< YOUR NAME >을 자신의 이름으로 교체), **Next**를 클릭합니다.

```
[ ] GitHub organisation: core
[ ] Workflow name: hello
[ ] A short description of your pipeline: A basic nf-core style version of Hello Nextflow
[ ] Name of the main author(s): < YOUR NAME >
```

- Template features 화면에서 `Toggle all features`를 **off**로 설정한 다음, 다음 항목을 선택적으로 **활성화**합니다. 선택 사항을 확인하고 **Continue**를 클릭합니다.

```
[ ] Add testing profiles
[ ] Use nf-core components
[ ] Use nf-schema
[ ] Add configuration files
[ ] Add documentation
```

- `Final details` 화면에서 **Finish**를 클릭합니다. 파이프라인이 생성될 때까지 기다린 다음 **Continue**를 클릭합니다.
- Create GitHub repository 화면에서 **Finish without creating a repo**를 클릭합니다. 이렇게 하면 나중에 GitHub 저장소를 생성하기 위한 지침이 표시됩니다. 이를 무시하고 **Close**를 클릭합니다.

TUI가 닫히면 다음과 같은 콘솔 출력이 표시됩니다.

??? success "명령 출력"

    ```console
                                              ,--./,-.
              ___     __   __   __   ___     /,-._.--~\
        |\ | |__  __ /  ` /  \ |__) |__         }  {
        | \| |       \__, \__/ |  \ |___     \`-._,-`-,
                                              `._,._,'

        nf-core/tools version 3.4.1 - https://nf-co.re


    INFO     Launching interactive nf-core pipeline creation tool.
    ```

파이프라인 생성이 성공했다는 명시적인 확인은 콘솔 출력에 없지만, `core-hello`라는 새 디렉토리가 보일 것입니다.

새 디렉토리의 내용을 확인하여 템플릿을 사용함으로써 얼마나 많은 작업을 절약했는지 확인하세요.

```bash
tree core-hello
```

??? abstract "디렉토리 내용"

    ```console
    core-hello/
    ├── assets
    │   ├── samplesheet.csv
    │   └── schema_input.json
    ├── conf
    │   ├── base.config
    │   ├── modules.config
    │   ├── test.config
    │   └── test_full.config
    ├── docs
    │   ├── output.md
    │   ├── README.md
    │   └── usage.md
    ├── main.nf
    ├── modules.json
    ├── nextflow.config
    ├── nextflow_schema.json
    ├── README.md
    ├── subworkflows
    │   ├── local
    │   │   └── utils_nfcore_hello_pipeline
    │   │       └── main.nf
    │   └── nf-core
    │       ├── utils_nextflow_pipeline
    │       │   ├── main.nf
    │       │   ├── meta.yml
    │       │   └── tests
    │       │       ├── main.function.nf.test
    │       │       ├── main.function.nf.test.snap
    │       │       ├── main.workflow.nf.test
    │       │       └── nextflow.config
    │       ├── utils_nfcore_pipeline
    │       │   ├── main.nf
    │       │   ├── meta.yml
    │       │   └── tests
    │       │       ├── main.function.nf.test
    │       │       ├── main.function.nf.test.snap
    │       │       ├── main.workflow.nf.test
    │       │       ├── main.workflow.nf.test.snap
    │       │       └── nextflow.config
    │       └── utils_nfschema_plugin
    │           ├── main.nf
    │           ├── meta.yml
    │           └── tests
    │               ├── main.nf.test
    │               ├── nextflow.config
    │               └── nextflow_schema.json
    └── workflows
        └── hello.nf

    14 directories, 34 files
    ```

정말 많은 파일입니다!

`nf-core/demo` 파이프라인 구조를 탐색할 때 만났던 것과 동일한 많은 파일을 인식하실 것입니다.
하지만 아직 조금 헷갈리더라도 걱정하지 마세요. 이 교육 과정에서 중요한 부분들을 함께 살펴보겠습니다.

!!! note "참고"

    이 교육의 첫 번째 부분에서 살펴본 `nf-core/demo` 파이프라인과 비교했을 때 한 가지 중요한 차이점은 `modules` 디렉토리가 없다는 것입니다.
    이는 기본 nf-core 모듈을 포함하지 않기로 선택했기 때문입니다.

### 1.2. 스캐폴드가 작동하는지 테스트

믿기 어렵겠지만, 실제 작업을 수행할 모듈을 아직 추가하지 않았음에도 불구하고, 파이프라인 스캐폴드는 실제로 `nf-core/demo` 파이프라인을 실행했던 것과 동일한 방식으로 test 프로필을 사용하여 실행할 수 있습니다.

```bash
nextflow run ./core-hello -profile docker,test --outdir core-hello-results
```

??? success "명령 출력"

    ```console
    N E X T F L O W   ~  version 25.04.3

    Launching `./core-hello/main.nf` [scruffy_marconi] DSL2 - revision: b9e9b3b8de

    Downloading plugin nf-schema@2.5.1
    Input/output options
      input                     : https://raw.githubusercontent.com/nf-core/test-datasets/viralrecon/samplesheet/samplesheet_test_illumina_amplicon.csv
      outdir                    : core-hello-results

    Institutional config options
      config_profile_name       : Test profile
      config_profile_description: Minimal test dataset to check pipeline function

    Generic options
      trace_report_suffix       : 2025-11-21_04-47-18

    Core Nextflow options
      runName                   : scruffy_marconi
      containerEngine           : docker
      launchDir                 : /workspaces/training/hello-nf-core
      workDir                   : /workspaces/training/hello-nf-core/work
      projectDir                : /workspaces/training/hello-nf-core/core-hello
      userName                  : root
      profile                   : docker,test
      configFiles               : /workspaces/training/hello-nf-core/core-hello/nextflow.config

    !! Only displaying parameters that differ from the pipeline defaults !!
    ------------------------------------------------------
    -[core/hello] Pipeline completed successfully-
    ```

이는 모든 기본 배선이 제자리에 있음을 보여줍니다.
그렇다면 출력은 어디에 있을까요? 출력이 있기는 한가요?

실제로 표준 실행 보고서를 포함하는 `core-hello-results`라는 새 결과 디렉토리가 생성되었습니다:

```bash
tree core-hello-results
```

??? abstract "디렉토리 내용"

    ```console
    core-hello-results
    └── pipeline_info
        ├── execution_report_2025-11-21_04-47-18.html
        ├── execution_timeline_2025-11-21_04-47-18.html
        ├── execution_trace_2025-11-21_04-47-18.txt
        ├── hello_software_versions.yml
        ├── params_2025-11-21_04-47-18.json
        └── pipeline_dag_2025-11-21_04-47-18.html

    1 directory, 6 files
    ```

보고서를 살펴보면 무엇이 실행되었는지 확인할 수 있는데, 답은 아무것도 실행되지 않았다는 것입니다!

![비어있는 실행 타임라인 보고서](./img/execution_timeline_empty.png)

코드에 실제로 무엇이 있는지 살펴보겠습니다.

### 1.3. 플레이스홀더 workflow 검사

`main.nf` 파일 내부를 보면 `workflows/hello`에서 `HELLO`라는 workflow를 가져오는 것을 볼 수 있습니다.

이것은 1부에서 만난 `workflows/demo.nf` workflow와 동등하며, 이미 일부 nf-core 기능이 갖춰진 관심 workflow의 플레이스홀더 역할을 합니다.

```groovy title="core-hello/workflows/hello.nf" linenums="1" hl_lines="15 17 19 35"
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { paramsSummaryMap       } from 'plugin/nf-schema'
include { softwareVersionsToYAML } from '../subworkflows/nf-core/utils_nfcore_pipeline'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow HELLO {

    take:
    ch_samplesheet // channel: --input에서 읽은 샘플시트
    main:

    ch_versions = channel.empty()

    //
    // 소프트웨어 버전 수집 및 저장
    //
    softwareVersionsToYAML(ch_versions)
        .collectFile(
            storeDir: "${params.outdir}/pipeline_info",
            name:  'hello_software_'  + 'versions.yml',
            sort: true,
            newLine: true
        ).set { ch_collated_versions }


    emit:
    versions       = ch_versions                 // channel: [ path(versions.yml) ]

}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    끝
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
```

[Hello Nextflow](../hello_nextflow/index.md)에서 개발된 것과 같은 기본 Nextflow workflow와 비교하면, 여기에서 새로운 몇 가지 사항(위의 강조 표시된 줄)을 발견할 수 있습니다:

- workflow 블록에 이름이 있습니다
- workflow 입력은 `take:` 키워드를 사용하여 선언되고 채널 구성은 상위 workflow로 이동됩니다
- workflow 내용은 `main:` 블록 내부에 배치됩니다
- 출력은 `emit:` 키워드를 사용하여 선언됩니다

이것들은 workflow를 **구성 가능(composable)**하게 만드는 Nextflow의 선택적 기능으로, 다른 workflow 내에서 호출될 수 있음을 의미합니다.

!!! note "구성 가능한 workflow 심화"

    [Workflows of Workflows](../side_quests/workflows_of_workflows.md) Side Quest는 여러 workflow를 함께 구성하고 그들 사이의 복잡한 데이터 흐름을 관리하는 방법을 포함하여 workflow 구성을 훨씬 더 심도 있게 탐구합니다. 여기서 구성 가능성을 소개하는 이유는 nf-core 템플릿 아키텍처의 기본 요구 사항이기 때문입니다. nf-core 템플릿 아키텍처는 중첩된 workflow를 사용하여 파이프라인 초기화, 메인 분석 workflow, 완료 작업을 별도의 재사용 가능한 컴포넌트로 구성합니다.

우리는 관심 있는 workflow의 관련 로직을 해당 구조에 연결해야 합니다.
그를 위한 첫 번째 단계는 원래 workflow를 구성 가능하게 만드는 것입니다.

### 요약

이제 nf-core 도구를 사용하여 파이프라인 스캐폴드를 생성하는 방법을 알게 되었습니다.

### 다음 단계

nf-core와 호환되도록 만들기 위한 전제로 간단한 workflow를 구성 가능하게 만드는 방법을 배웁니다.

---

## 2. 원래 Hello Nextflow workflow를 구성 가능하게 만들기

이제 workflow를 nf-core 스캐폴드에 통합하는 작업을 시작할 때입니다.
상기하자면, 우리는 [Hello Nextflow](../hello_nextflow/index.md) 교육 과정에서 다룬 workflow로 작업하고 있습니다.

!!! tip "팁"

    해당 파이프라인에 익숙하지 않거나 복습이 필요하시면 [Hello 파이프라인](../info/hello_pipeline.md)을 참조하세요.

완성된 Hello Nextflow workflow의 깨끗하고 완전히 작동하는 복사본을 모듈 및 입력으로 사용할 것으로 예상되는 기본 CSV 파일과 함께 `original-hello` 디렉토리에 제공합니다.

```bash
tree original-hello/
```

??? abstract "디렉토리 내용"

    ```console
    original-hello/
    ├── hello.nf
    ├── modules
    │   ├── collectGreetings.nf
    │   ├── convertToUpper.nf
    │   ├── cowpy.nf
    │   └── sayHello.nf
    └── nextflow.config

    1 directory, 6 files
    ```

작동하는지 확인하기 위해 자유롭게 실행해 보세요:

```bash
nextflow run original-hello/hello.nf
```

??? success "명령 출력"

    ```console
    N E X T F L O W   ~  version 25.04.3

    Launching `original-hello/hello.nf` [goofy_babbage] DSL2 - revision: e9e72441e9

    executor >  local (8)
    [a4/081cec] sayHello (1)       | 3 of 3 ✔
    [e7/7e9058] convertToUpper (3) | 3 of 3 ✔
    [0c/17263b] collectGreetings   | 1 of 1 ✔
    [94/542280] cowpy              | 1 of 1 ✔
    ```

`hello.nf` workflow 파일을 열어 아래에 전체적으로 표시된 코드를 검사해 봅시다(모듈에 있는 프로세스는 제외):

```groovy title="original-hello/hello.nf" linenums="1"
#!/usr/bin/env nextflow

/*
* Pipeline parameters
*/
params.greeting = 'greetings.csv'
params.batch = 'test-batch'
params.character = 'turkey'

// 모듈 포함
include { sayHello } from './modules/sayHello.nf'
include { convertToUpper } from './modules/convertToUpper.nf'
include { collectGreetings } from './modules/collectGreetings.nf'
include { cowpy } from './modules/cowpy.nf'

workflow {

  // CSV 파일에서 입력용 채널 생성
  greeting_ch = channel.fromPath(params.greeting)
                      .splitCsv()
                      .map { line -> line[0] }

  // 인사말 출력
  sayHello(greeting_ch)

  // 인사말을 대문자로 변환
  convertToUpper(sayHello.out)

  // 모든 인사말을 하나의 파일에 수집
  collectGreetings(convertToUpper.out.collect(), params.batch)

  // cowpy로 인사말의 ASCII 아트 생성
  cowpy(collectGreetings.out.outfile, params.character)
}
```

보시다시피, 이 workflow는 단독으로 실행할 수 있는 간단한 이름 없는 workflow로 작성되었습니다.
nf-core 템플릿이 요구하는 대로 상위 workflow 내에서 실행 가능하게 만들려면 **구성 가능**하게 만들어야 합니다.

필요한 변경 사항을 하나씩 살펴보겠습니다.

### 2.1. workflow에 이름 지정

먼저 상위 workflow에서 참조할 수 있도록 workflow에 이름을 지정하겠습니다.

=== "변경 후"

    ```groovy title="original-hello/hello.nf" linenums="16"
    workflow HELLO {
    ```

=== "변경 전"

    ```groovy title="original-hello/hello.nf" linenums="16"
    workflow {
    ```

모듈 이름과 동일한 규칙이 workflow 이름에도 적용됩니다.

### 2.2. 채널 구성을 `take`로 교체

이제 채널 구성을 예상 입력을 선언하는 간단한 `take` 문으로 교체합니다.

=== "변경 후"

    ```groovy title="original-hello/hello.nf" linenums="18"
        take:
        // 인사말 채널
        greeting_ch
    ```

=== "변경 전"

    ```groovy title="original-hello/hello.nf" linenums="18"
        // CSV 파일에서 입력용 채널 생성
        greeting_ch = channel.fromPath(params.greeting)
                            .splitCsv()
                            .map { line -> line[0] }
    ```

이는 입력이 어떻게 제공되는지에 대한 세부 사항을 상위 workflow에 맡깁니다.

이 기회에 `params.greeting = 'greetings.csv'` 줄도 주석 처리할 수 있습니다.

=== "변경 후"

    ```groovy title="original-hello/hello.nf" linenums="3" hl_lines="4"
        /*
        * Pipeline parameters
        */
        //params.greeting = 'greetings.csv'
        params.batch = 'test-batch'
        params.character = 'turkey'
    ```

=== "변경 전"

    ```groovy title="original-hello/hello.nf" linenums="3" hl_lines="4"
        /*
        * Pipeline parameters
        */
        params.greeting = 'greetings.csv'
        params.batch = 'test-batch'
        params.character = 'turkey'
    ```

!!! note "참고"

    Nextflow 언어 서버 확장이 설치되어 있다면 구문 검사기가 코드에 빨간 물결선을 표시할 것입니다.
    이는 `take:` 문을 넣으면 `main:`도 있어야 하기 때문입니다.

    다음 단계에서 추가하겠습니다.

### 2.3. workflow 작업 앞에 `main` 문 추가

다음으로, workflow 본문에서 호출되는 나머지 작업 앞에 `main` 문을 추가합니다.

=== "변경 후"

    ```groovy title="original-hello/hello.nf" linenums="22" hl_lines="1"
        main:

        // 인사말 출력
        sayHello(greeting_ch)

        // 인사말을 대문자로 변환
        convertToUpper(sayHello.out)

        // 모든 인사말을 하나의 파일에 수집
        collectGreetings(convertToUpper.out.collect(), params.batch)

        // cowpy로 인사말의 ASCII 아트 생성
        cowpy(collectGreetings.out.outfile, params.character)
    ```

=== "변경 전"

    ```groovy title="original-hello/hello.nf" linenums="21"
        // 인사말 출력
        sayHello(greeting_ch)

        // 인사말을 대문자로 변환
        convertToUpper(sayHello.out)

        // 모든 인사말을 하나의 파일에 수집
        collectGreetings(convertToUpper.out.collect(), params.batch)

        // cowpy로 인사말의 ASCII 아트 생성
        cowpy(collectGreetings.out.outfile, params.character)
    ```

이것은 기본적으로 '이것이 이 workflow가 _하는_ 일'이라고 말하는 것입니다.

### 2.4. `emit` 문 추가

마지막으로 workflow의 최종 출력이 무엇인지 선언하는 `emit` 문을 추가합니다.

```groovy title="original-hello/hello.nf" linenums="35"
    emit:
    cowpy_hellos = cowpy.out
```

이것은 원래 workflow와 비교하여 코드에 새롭게 추가된 부분입니다.

### 2.5. 완료된 변경 사항 요약

설명대로 모든 변경을 수행했다면 workflow는 이제 다음과 같이 보일 것입니다:

```groovy title="original-hello/hello.nf" linenums="1" hl_lines="16 18-20 22 36-37"
#!/usr/bin/env nextflow

/*
* Pipeline parameters
*/
// params.greeting = 'greetings.csv'
params.batch = 'test-batch'
params.character = 'turkey'

// 모듈 포함
include { sayHello } from './modules/sayHello.nf'
include { convertToUpper } from './modules/convertToUpper.nf'
include { collectGreetings } from './modules/collectGreetings.nf'
include { cowpy } from './modules/cowpy.nf'

workflow HELLO {

    take:
    // 인사말 채널
    greeting_ch

    main:

    // 인사말 출력
    sayHello(greeting_ch)

    // 인사말을 대문자로 변환
    convertToUpper(sayHello.out)

    // 모든 인사말을 하나의 파일에 수집
    collectGreetings(convertToUpper.out.collect(), params.batch)

    // cowpy로 인사말의 ASCII 아트 생성
    cowpy(collectGreetings.out.outfile, params.character)

    emit:
    cowpy_hellos = cowpy.out
}
```

이는 입력 채널에 무엇을 공급할지를 제외하고 Nextflow가 필요로 하는 모든 것을 설명합니다.
그것은 **진입점(entrypoint)** workflow라고도 불리는 상위 workflow에서 정의될 것입니다.

### 2.6. 더미 진입점 workflow 만들기

구성 가능한 workflow를 복잡한 nf-core 스캐폴드에 통합하기 전에, 올바르게 작동하는지 확인해봅시다.
구성 가능한 workflow를 독립적으로 테스트하기 위해 간단한 더미 진입점 workflow를 만들 수 있습니다.

동일한 `original-hello` 디렉토리에 `main.nf`라는 빈 파일을 만듭니다.

```bash
touch original-hello/main.nf
```

다음 코드를 `main.nf` 파일에 복사합니다.

```groovy title="original-hello/main.nf" linenums="1"
#!/usr/bin/env nextflow

// hello.nf 파일에서 workflow 코드 가져오기
include { HELLO } from './hello.nf'

// 입력 매개변수 선언
params.greeting = 'greetings.csv'

workflow {
  // CSV 파일에서 입력용 채널 생성
  greeting_ch = channel.fromPath(params.greeting)
                      .splitCsv()
                      .map { line -> line[0] }

  // 인사말 채널에서 가져온 workflow 호출
  HELLO(greeting_ch)

  // workflow에서 방출된 출력 보기
  HELLO.out.view { output -> "Output: $output" }
}
```

여기서 두 가지 중요한 관찰 사항이 있습니다:

- 가져온 workflow를 호출하는 구문은 모듈을 호출하는 구문과 본질적으로 동일합니다.
- 입력을 workflow로 가져오는 것과 관련된 모든 것(입력 매개변수 및 채널 구성)은 이제 이 상위 workflow에 선언됩니다.

!!! note "참고"

    진입점 workflow 파일의 이름을 `main.nf`로 지정하는 것은 규칙이지 요구 사항이 아닙니다.

    이 규칙을 따르면 `nextflow run` 명령에서 workflow 파일 이름을 지정하지 않아도 됩니다.
    Nextflow는 실행 디렉토리에서 `main.nf`라는 파일을 자동으로 찾습니다.

    그러나 원한다면 진입점 workflow 파일의 이름을 다르게 지정할 수 있습니다.
    그 경우 `nextflow run` 명령에서 workflow 파일 이름을 반드시 지정해야 합니다.

### 2.7. workflow가 실행되는지 테스트

마침내 구성 가능한 workflow가 작동하는지 확인하는 데 필요한 모든 조각을 갖추게 되었습니다.
실행해 봅시다!

```bash
nextflow run ./original-hello
```

여기서 `main.nf` 명명 규칙을 사용하는 이점을 볼 수 있습니다.
진입점 workflow의 이름을 `something_else.nf`로 지정했다면 `nextflow run original-hello/something_else.nf`를 해야 했을 것입니다.

모든 변경 사항을 올바르게 수행했다면 완료될 때까지 실행되어야 합니다.

??? success "명령 출력"

    ```console
    N E X T F L O W   ~  version 25.04.3

    Launching `original-hello/main.nf` [friendly_wright] DSL2 - revision: 1ecd2d9c0a

    executor >  local (8)
    [24/c6c0d8] HELLO:sayHello (3)       | 3 of 3 ✔
    [dc/721042] HELLO:convertToUpper (3) | 3 of 3 ✔
    [48/5ab2df] HELLO:collectGreetings   | 1 of 1 ✔
    [e3/693b7e] HELLO:cowpy              | 1 of 1 ✔
    Output: /workspaces/training/hello-nf-core/work/e3/693b7e48dc119d0c54543e0634c2e7/cowpy-COLLECTED-test-batch-output.txt
    ```

이는 HELLO workflow를 구성 가능하게 업그레이드하는 데 성공했음을 의미합니다.

### 요약

이름을 지정하고 `take`, `main`, `emit` 문을 추가하여 workflow를 구성 가능하게 만드는 방법과 진입점 workflow에서 호출하는 방법을 알게 되었습니다.

### 다음 단계

nf-core 호환 가능하게 만들기 위한 전제로 기본 구성 가능한 workflow를 nf-core 스캐폴드에 접목하는 방법을 배웁니다.

---

## 3. 업데이트된 workflow 로직을 플레이스홀더 workflow에 맞추기

이제 구성 가능한 workflow가 올바르게 작동하는 것을 확인했으므로, 섹션 1에서 생성한 nf-core 파이프라인 스캐폴드로 돌아가겠습니다.
방금 개발한 구성 가능한 workflow를 nf-core 템플릿 구조에 통합하여 최종 결과가 다음과 같이 보이기를 원합니다.

<figure class="excalidraw">
--8<-- "docs/en/docs/hello_nf-core/img/core-hello.svg"
</figure>

그렇다면 어떻게 이를 실현할까요? `core-hello/workflows/hello.nf`(nf-core 스캐폴드)에 있는 `HELLO` workflow의 현재 내용을 살펴보겠습니다.

```groovy title="core-hello/workflows/hello.nf" linenums="1"
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { paramsSummaryMap       } from 'plugin/nf-schema'
include { softwareVersionsToYAML } from '../subworkflows/nf-core/utils_nfcore_pipeline'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow HELLO {

    take:
    ch_samplesheet // channel: --input에서 읽어온 samplesheet
    main:

    ch_versions = channel.empty()

    //
    // 소프트웨어 버전 수집 및 저장
    //
    softwareVersionsToYAML(ch_versions)
        .collectFile(
            storeDir: "${params.outdir}/pipeline_info",
            name:  'hello_software_'  + 'versions.yml',
            sort: true,
            newLine: true
        ).set { ch_collated_versions }


    emit:
    versions       = ch_versions                 // channel: [ path(versions.yml) ]

}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    끝
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
```

전반적으로 이 코드는 파이프라인에서 실행되는 소프트웨어 도구의 버전을 캡처하는 것과 관련된 일부 관리 작업을 제외하고는 거의 아무것도 하지 않습니다.

섹션 2에서 개발한 원래 workflow의 구성 가능한 버전에서 관련 코드를 추가해야 합니다.

다음 단계로 이 작업을 처리할 것입니다:

1. 모듈을 복사하고 모듈 import 설정
2. `take` 선언은 그대로 유지
3. workflow 로직을 `main` 블록에 추가
4. `emit` 블록 업데이트

!!! note "참고"

    이번 첫 번째 패스에서는 버전 캡처를 무시하고 나중에 이 교육의 뒷부분에서 이를 연결하는 방법을 살펴볼 것입니다.

### 3.1. 모듈 복사 및 모듈 import 설정

Hello Nextflow workflow의 네 가지 프로세스는 `original-hello/modules/`에 모듈로 저장됩니다.
이러한 모듈을 nf-core 프로젝트 구조(`core-hello/modules/local/` 아래)에 복사하고 nf-core workflow 파일에 import 문을 추가해야 합니다.

먼저 `original-hello/`에서 `core-hello/`로 모듈 파일을 복사하겠습니다:

```bash
mkdir -p core-hello/modules/local/
cp original-hello/modules/* core-hello/modules/local/.
```

이제 `core-hello/` 아래에 나열된 모듈 디렉토리를 볼 수 있어야 합니다.

```bash
tree core-hello/modules
```

??? abstract "디렉토리 내용"

    ```console
    core-hello/modules
    └── local
        ├── collectGreetings.nf
        ├── convertToUpper.nf
        ├── cowpy.nf
        └── sayHello.nf

    1 directory, 4 files
    ```

이제 모듈 import 문을 설정하겠습니다.

`original-hello/hello.nf` workflow의 import 문은 다음과 같았습니다:

```groovy title="original-hello/hello.nf" linenums="9"
// 모듈 포함
include { sayHello } from './modules/sayHello.nf'
include { convertToUpper } from './modules/convertToUpper.nf'
include { collectGreetings } from './modules/collectGreetings.nf'
include { cowpy } from './modules/cowpy.nf'
```

`core-hello/workflows/hello.nf` 파일을 열고 아래와 같이 이러한 import 문을 변환합니다.

=== "변경 후"

    ```groovy title="core-hello/workflows/hello.nf" linenums="1" hl_lines="8-11"
    /*
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    */
    include { paramsSummaryMap       } from 'plugin/nf-schema'
    include { softwareVersionsToYAML } from '../subworkflows/nf-core/utils_nfcore_pipeline'
    include { sayHello               } from '../modules/local/sayHello.nf'
    include { convertToUpper         } from '../modules/local/convertToUpper.nf'
    include { collectGreetings       } from '../modules/local/collectGreetings.nf'
    include { cowpy                  } from '../modules/local/cowpy.nf'
    ```

=== "변경 전"

    ```groovy title="core-hello/workflows/hello.nf" linenums="1"
    /*
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    */
    include { paramsSummaryMap       } from 'plugin/nf-schema'
    include { softwareVersionsToYAML } from '../subworkflows/nf-core/utils_nfcore_pipeline'
    ```

여기서 두 가지 더 흥미로운 관찰 사항이 있습니다:

- import 문의 형식을 nf-core 스타일 규칙을 따르도록 조정했습니다.
- 이제 다른 중첩 수준에 저장되어 있음을 반영하도록 모듈에 대한 상대 경로를 업데이트했습니다.

### 3.2. `take` 선언은 그대로 유지

nf-core 프로젝트에는 일반적으로 열 데이터를 포함하는 CSV 파일인 samplesheet 개념과 관련된 많은 사전 구축 기능이 있습니다.
본질적으로 우리의 `greetings.csv` 파일이 그러하므로 현재 `take` 선언을 그대로 유지하고 다음 단계에서 입력 채널의 이름만 업데이트하겠습니다.

```groovy title="core-hello/workflows/hello.nf" linenums="21"
    take:
    ch_samplesheet // channel: --input에서 읽어온 samplesheet
```

입력 처리는 이 workflow의 상위에서 수행됩니다(이 코드 파일에서가 아님).

### 3.3. workflow 로직을 `main` 블록에 추가

이제 모듈을 workflow에서 사용할 수 있으므로 workflow 로직을 `main` 블록에 연결할 수 있습니다.

참고로, 이것은 원래 workflow의 관련 코드이며, 구성 가능하게 만들 때 많이 변경되지 않았습니다(`main:` 줄만 추가했습니다):

```groovy title="original-hello/hello.nf" linenums="22"
    main:

    // 인사말 출력
    sayHello(greeting_ch)

    // 인사말을 대문자로 변환
    convertToUpper(sayHello.out)

    // 모든 인사말을 하나의 파일에 수집
    collectGreetings(convertToUpper.out.collect(), params.batch)

    // cowpy로 인사말의 ASCII 아트 생성
    cowpy(collectGreetings.out.outfile, params.character)
```

`main:` 뒤에 오는 코드를 workflow의 새 버전에 복사해야 합니다.

workflow를 실행하는 도구의 버전을 캡처하는 것과 관련된 일부 코드가 이미 있습니다. 지금은 그대로 두겠습니다(도구 버전은 나중에 처리하겠습니다).
맨 위에 `ch_versions = channel.empty()` 초기화를 유지한 다음 workflow 로직을 삽입하고 버전 수집 코드를 끝에 유지하겠습니다.
이러한 순서는 실제 파이프라인에서 프로세스가 workflow가 실행될 때 `ch_versions` 채널에 추가될 버전 정보를 방출하기 때문에 의미가 있습니다.

=== "변경 후"

    ```groovy title="core-hello/workflows/hello.nf" linenums="19" hl_lines="10-20"
    workflow HELLO {

        take:
        ch_samplesheet // channel: --input에서 읽어온 samplesheet

        main:

        ch_versions = Channel.empty()

        // 인사말 출력
        sayHello(greeting_ch)

        // 인사말을 대문자로 변환
        convertToUpper(sayHello.out)

        // 모든 인사말을 하나의 파일에 수집
        collectGreetings(convertToUpper.out.collect(), params.batch)

        // cowpy로 인사말의 ASCII 아트 생성
        cowpy(collectGreetings.out.outfile, params.character)

        //
        // 소프트웨어 버전 수집 및 저장
        //
        softwareVersionsToYAML(ch_versions)
            .collectFile(
                storeDir: "${params.outdir}/pipeline_info",
                name:  'hello_software_'  + 'versions.yml',
                sort: true,
                newLine: true
            ).set { ch_collated_versions }


        emit:
        versions       = ch_versions                 // channel: [ path(versions.yml) ]

    }
    ```

=== "변경 전"

    ```groovy title="core-hello/workflows/hello.nf" linenums="19"
    workflow HELLO {

        take:
    ch_samplesheet // channel: --input에서 읽어온 samplesheet
    main:

    ch_versions = channel.empty()

    //
    // 소프트웨어 버전 수집 및 저장
    //
    softwareVersionsToYAML(ch_versions)
            .collectFile(
                storeDir: "${params.outdir}/pipeline_info",
                name:  'hello_software_'  + 'versions.yml',
                sort: true,
                newLine: true
            ).set { ch_collated_versions }


        emit:
        versions       = ch_versions                 // channel: [ path(versions.yml) ]

    }
    ```

`main:` 앞에 빈 줄을 추가하여 코드를 더 읽기 쉽게 만들었습니다.

좋아 보이지만 `take:` 키워드 아래에 작성된 것과 일치하도록 아래와 같이 `sayHello()` 프로세스에 전달하는 채널의 이름을 `greeting_ch`에서 `ch_samplesheet`로 업데이트해야 합니다.

=== "변경 후"

    ```groovy title="core-hello/workflows/hello.nf" linenums="26"
        // 인사말 출력 (nf-core samplesheet 규칙을 사용하도록 업데이트됨)
        sayHello(ch_samplesheet)
    ```

=== "변경 전"

    ```groovy title="core-hello/workflows/hello.nf" linenums="26"
        // 인사말 출력
        sayHello(greeting_ch)
    ```

이제 workflow 로직이 올바르게 연결되었습니다.

### 3.4. `emit` 블록 업데이트

마지막으로 workflow의 최종 출력 선언을 포함하도록 `emit` 블록을 업데이트해야 합니다.

=== "변경 후"

    ```groovy title="core-hello/workflows/hello.nf" linenums="55" hl_lines="2"
        emit:
        cowpy_hellos   = cowpy.out
        versions       = ch_versions                 // channel: [ path(versions.yml) ]
    ```

=== "변경 전"

    ```groovy title="core-hello/workflows/hello.nf" linenums="55"
        emit:
        versions       = ch_versions                 // channel: [ path(versions.yml) ]
    ```

이것으로 HELLO workflow 자체에 대한 수정이 완료됩니다.
이 시점에서 우리는 구현하기로 설정한 전반적인 코드 구조를 달성했습니다.

### 요약

구성 가능한 workflow의 핵심 부분을 nf-core 플레이스홀더 workflow에 맞추는 방법을 알게 되었습니다.

### 다음 단계

nf-core 파이프라인 스캐폴드에서 입력이 처리되는 방식을 조정하는 방법을 배웁니다.

---

## 4. 입력 처리 조정

workflow 로직을 nf-core 스캐폴드에 성공적으로 통합했으므로 이제 한 가지 더 중요한 부분을 처리해야 합니다: 입력 데이터가 올바르게 처리되도록 하는 것입니다.
nf-core 템플릿은 복잡한 유전체학 데이터 세트를 위해 설계된 정교한 입력 처리와 함께 제공되므로 더 간단한 `greetings.csv` 파일과 작동하도록 조정해야 합니다.

### 4.1. 입력이 처리되는 위치 식별

첫 번째 단계는 입력 처리가 어디에서 수행되는지 파악하는 것입니다.

Hello Nextflow workflow를 구성 가능하게 다시 작성할 때 입력 매개변수 선언을 `main.nf` 진입점 workflow로 한 단계 위로 이동한 것을 기억하실 것입니다.
따라서 파이프라인 스캐폴드의 일부로 생성된 최상위 `main.nf` 진입점 workflow를 살펴보겠습니다:

```groovy title="core-hello/main.nf" linenums="1"
#!/usr/bin/env nextflow
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    core/hello
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Github : https://github.com/core/hello
----------------------------------------------------------------------------------------
*/

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT FUNCTIONS / MODULES / SUBWORKFLOWS / WORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { HELLO  } from './workflows/hello'
include { PIPELINE_INITIALISATION } from './subworkflows/local/utils_nfcore_hello_pipeline'
include { PIPELINE_COMPLETION     } from './subworkflows/local/utils_nfcore_hello_pipeline'
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    NAMED WORKFLOWS FOR PIPELINE
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// WORKFLOW: 입력 유형에 따라 메인 분석 파이프라인 실행
//
workflow CORE_HELLO {

    take:
    samplesheet // channel: --input에서 읽어온 samplesheet

    main:

    //
    // WORKFLOW: 파이프라인 실행
    //
    HELLO (
        samplesheet
    )
}
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    메인 WORKFLOW 실행
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow {

    main:
    //
    // SUBWORKFLOW: 초기화 작업 실행
    //
    PIPELINE_INITIALISATION (
        params.version,
        params.validate_params,
        params.monochrome_logs,
        args,
        params.outdir,
        params.input
    )

    //
    // WORKFLOW: 메인 워크플로우 실행
    //
    CORE_HELLO (
        PIPELINE_INITIALISATION.out.samplesheet
    )
    //
    // SUBWORKFLOW: 완료 작업 실행
    //
    PIPELINE_COMPLETION (
        params.outdir,
        params.monochrome_logs,
    )
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    끝
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
```

nf-core 프로젝트는 중첩된 subworkflow를 많이 사용하므로 이 부분은 처음 접근할 때 약간 혼란스러울 수 있습니다.

여기서 중요한 것은 두 개의 workflow가 정의되어 있다는 것입니다:

- `CORE_HELLO`는 `core-hello/workflows/hello.nf`에서 방금 조정을 완료한 HELLO workflow를 실행하기 위한 얇은 래퍼입니다.
- `CORE_HELLO`와 `PIPELINE_INITIALISATION`, `PIPELINE_COMPLETION`이라는 두 개의 다른 subworkflow를 호출하는 이름 없는 workflow입니다.

그들이 서로 어떻게 관련되어 있는지 보여주는 다이어그램이 있습니다:

<figure class="excalidraw">
--8<-- "docs/en/docs/hello_nf-core/img/hello-nested-workflows.svg"
</figure>

중요한 것은 이 수준에서 입력 채널을 구성하는 코드를 찾을 수 없고, `--input` 매개변수를 통해 제공되는 samplesheet에 대한 참조만 있다는 것입니다.

약간의 조사 결과 입력 처리는 적절하게도 `PIPELINE_INITIALISATION` subworkflow에서 수행되며, 이는 `core-hello/subworkflows/local/utils_nfcore_hello_pipeline/main.nf`에서 가져온다는 것이 드러납니다.

해당 파일을 열고 아래로 스크롤하면 다음 코드 덩어리에 도달합니다:

```groovy title="core-hello/subworkflows/local/utils_nfcore_hello_pipeline/main.nf" linenums="76"
    //
    // params.input을 통해 제공된 입력 파일에서 채널 생성
    //

    channel
        .fromList(samplesheetToList(params.input, "${projectDir}/assets/schema_input.json"))
        .map {
            meta, fastq_1, fastq_2 ->
                if (!fastq_2) {
                    return [ meta.id, meta + [ single_end:true ], [ fastq_1 ] ]
                } else {
                    return [ meta.id, meta + [ single_end:false ], [ fastq_1, fastq_2 ] ]
                }
        }
        .groupTuple()
        .map { samplesheet ->
            validateInputSamplesheet(samplesheet)
        }
        .map {
            meta, fastqs ->
                return [ meta, fastqs.flatten() ]
        }
        .set { ch_samplesheet }

    emit:
    samplesheet = ch_samplesheet
    versions    = ch_versions
```

이것은 samplesheet를 구문 분석하고 HELLO workflow에서 소비할 준비가 된 형태로 전달하는 channel factory입니다.

!!! note "참고"

    위의 구문은 이전에 사용한 것과 약간 다르지만 기본적으로 다음은:

    ```groovy
    channel.<...>.set { ch_samplesheet }
    ```

    다음과 동등합니다:

    ```groovy
    ch_samplesheet = channel.<...>
    ```

이 코드에는 nf-core 파이프라인 템플릿에 포함된 예제 samplesheet에 매우 특정적인 일부 구문 분석 및 검증 단계가 포함되어 있으며, 작성 당시 이는 매우 도메인별로 특정적이어서 우리의 간단한 파이프라인 프로젝트에 적합하지 않습니다.

### 4.2. 템플릿 입력 채널 코드 교체

좋은 소식은 우리 파이프라인의 요구 사항이 훨씬 더 간단하므로 원래 Hello Nextflow workflow에서 개발한 채널 구성 코드로 이 모든 것을 교체할 수 있다는 것입니다.

참고로 채널 구성은 다음과 같았습니다(solutions 디렉토리에서 확인):

```groovy title="solutions/composable-hello/main.nf" linenums="10" hl_lines="4"
    // CSV 파일에서 입력용 채널 생성
    greeting_ch = channel.fromPath(params.greeting)
        .splitCsv()
        .map { line -> line[0] }
```

따라서 약간의 변경을 통해 초기화 workflow에 연결하기만 하면 됩니다: 채널 이름을 `greeting_ch`에서 `ch_samplesheet`로, 매개변수 이름을 `params.greeting`에서 `params.input`으로 업데이트합니다(강조 표시된 줄 참조).

=== "변경 후"

    ```groovy title="core-hello/subworkflows/local/utils_nfcore_hello_pipeline/main.nf" linenums="76" hl_lines="5-7"
        //
        // params.input을 통해 제공된 입력 파일에서 채널 생성
        //

        ch_samplesheet = channel.fromPath(params.input)
            .splitCsv()
            .map { line -> line[0] }

        emit:
        samplesheet = ch_samplesheet
        versions    = ch_versions
    ```

=== "변경 전"

    ```groovy title="core-hello/subworkflows/local/utils_nfcore_hello_pipeline/main.nf" linenums="76" hl_lines="5-23"
        //
        // params.input을 통해 제공된 입력 파일에서 채널 생성
        //

        channel
            .fromList(samplesheetToList(params.input, "${projectDir}/assets/schema_input.json"))
            .map {
                meta, fastq_1, fastq_2 ->
                    if (!fastq_2) {
                        return [ meta.id, meta + [ single_end:true ], [ fastq_1 ] ]
                    } else {
                        return [ meta.id, meta + [ single_end:false ], [ fastq_1, fastq_2 ] ]
                    }
            }
            .groupTuple()
            .map { samplesheet ->
                validateInputSamplesheet(samplesheet)
            }
            .map {
                meta, fastqs ->
                    return [ meta, fastqs.flatten() ]
            }
            .set { ch_samplesheet }

        emit:
        samplesheet = ch_samplesheet
        versions    = ch_versions
    ```

이것으로 입력 처리가 작동하도록 만드는 데 필요한 변경 사항이 완료됩니다.

현재 형태에서는 스키마 검증을 위한 nf-core의 내장 기능을 활용할 수 없지만 나중에 추가할 수 있습니다.
지금은 테스트 데이터에서 성공적으로 실행할 수 있는 것을 얻기 위해 가능한 한 간단하게 유지하는 데 중점을 둡니다.

### 4.3. 테스트 프로필 업데이트

테스트 데이터와 매개변수에 대해 말하자면, 템플릿에 제공된 예제 samplesheet 대신 `greetings.csv` 미니 samplesheet를 사용하도록 이 파이프라인의 테스트 프로필을 업데이트하겠습니다.

`core-hello/conf` 아래에는 작은 데이터 샘플과 전체 크기 데이터를 테스트하기 위한 두 개의 템플릿 테스트 프로필인 `test.config`와 `test_full.config`가 있습니다.
우리 파이프라인의 목적을 고려할 때 전체 크기 테스트 프로필을 설정하는 데는 실제로 의미가 없으므로 `test_full.config`를 무시하거나 삭제하셔도 됩니다.
몇 가지 기본 매개변수로 `greetings.csv` 파일에서 실행되도록 `test.config`를 설정하는 데 집중하겠습니다.

#### 4.3.1. `greetings.csv` 파일 복사

먼저 파이프라인 프로젝트의 적절한 위치에 `greetings.csv` 파일을 복사해야 합니다.
일반적으로 작은 테스트 파일은 `assets` 디렉토리에 저장되므로 작업 디렉토리에서 파일을 복사하겠습니다.

```bash
cp greetings.csv core-hello/assets/.
```

이제 `greetings.csv` 파일을 테스트 입력으로 사용할 준비가 되었습니다.

#### 4.3.2. `test.config` 파일 업데이트

이제 다음과 같이 `test.config` 파일을 업데이트할 수 있습니다:

=== "변경 후"

    ```groovy title="core-hello/conf/test.config" linenums="21" hl_lines="6-10"
    params {
        config_profile_name        = 'Test profile'
        config_profile_description = 'Minimal test dataset to check pipeline function'

        // 입력 데이터
        input  = "${projectDir}/assets/greetings.csv"

        // 기타 매개변수
        batch     = 'test'
        character = 'tux'
    }
    ```

=== "변경 전"

    ```groovy title="core-hello/conf/test.config" linenums="21" hl_lines="6-8"
    params {
        config_profile_name        = 'Test profile'
        config_profile_description = 'Minimal test dataset to check pipeline function'

        // 입력 데이터
        // TODO nf-core: nf-core/test-datasets에 테스트 데이터 경로 지정
        // TODO nf-core: 명령줄 플래그가 필요하지 않도록 테스트에 필요한 매개변수 제공
        input  = params.pipelines_testdata_base_path + 'viralrecon/samplesheet/samplesheet_test_illumina_amplicon.csv'
    }
    ```

주요 요점:

- **`${projectDir}` 사용**: 이는 메인 workflow 스크립트가 위치한 디렉토리(파이프라인 루트)를 가리키는 Nextflow 암시적 변수입니다. 이를 사용하면 파이프라인이 어디서 실행되든 경로가 작동합니다.
- **절대 경로**: `${projectDir}`를 사용함으로써 절대 경로를 만들며, 이는 파이프라인과 함께 제공되는 테스트 데이터에 중요합니다.
- **테스트 데이터 위치**: nf-core 파이프라인은 일반적으로 작은 테스트 파일의 경우 파이프라인 저장소 내 `assets/` 디렉토리에, 더 큰 파일의 경우 외부 테스트 데이터셋에 테스트 데이터를 저장합니다.

또한 매우 기본적인 머신(Github Codespaces의 최소 VM과 같은)에서도 실행되도록 기본 리소스 제한을 강화하겠습니다:

=== "변경 후"

    ```groovy title="core-hello/config/test.config" linenums="13" hl_lines="3-4"
    process {
        resourceLimits = [
            cpus: 2,
            memory: '4.GB',
            time: '1.h'
        ]
    }
    ```

=== "변경 전"

    ```groovy title="core-hello/config/test.config" linenums="13" hl_lines="3-4"
    process {
        resourceLimits = [
            cpus: 4,
            memory: '15.GB',
            time: '1.h'
        ]
    }
    ```

이것으로 필요한 코드 수정이 완료되었습니다.

### 4.4. 테스트 프로필로 파이프라인 실행

많은 작업이었지만 이제 마침내 파이프라인을 실행해 볼 수 있습니다!
아직 검증을 설정하지 않았기 때문에 명령줄에 `--validate_params false`를 추가해야 합니다(이는 나중에 다룰 것입니다).

```bash
nextflow run core-hello --outdir core-hello-results -profile test,docker --validate_params false
```

모든 수정을 올바르게 수행했다면 완료될 때까지 실행되어야 합니다.

??? success "명령 출력"

    ```console
     N E X T F L O W   ~  version 25.04.3

    Launching `core-hello/main.nf` [condescending_allen] DSL2 - revision: b9e9b3b8de

    Input/output options
      input                     : /workspaces/training/hello-nf-core/core-hello/assets/greetings.csv
      outdir                    : core-hello-results

    Institutional config options
      config_profile_name       : Test profile
      config_profile_description: Minimal test dataset to check pipeline function

    Generic options
      validate_params           : false
      trace_report_suffix       : 2025-11-21_07-29-37

    Core Nextflow options
      runName                   : condescending_allen
      containerEngine           : docker
      launchDir                 : /workspaces/training/hello-nf-core
      workDir                   : /workspaces/training/hello-nf-core/work
      projectDir                : /workspaces/training/hello-nf-core/core-hello
      userName                  : root
      profile                   : test,docker
      configFiles               : /workspaces/training/hello-nf-core/core-hello/nextflow.config

    !! Only displaying parameters that differ from the pipeline defaults !!
    ------------------------------------------------------
    executor >  local (1)
    [ed/727b7e] CORE_HELLO:HELLO:sayHello (3)       [100%] 3 of 3 ✔
    [45/bb6096] CORE_HELLO:HELLO:convertToUpper (3) [100%] 3 of 3 ✔
    [81/7e2e34] CORE_HELLO:HELLO:collectGreetings   [100%] 1 of 1 ✔
    [96/9442a1] CORE_HELLO:HELLO:cowpy              [100%] 1 of 1 ✔
    -[core/hello] Pipeline completed successfully-
    ```

보시다시피, 초기화 subworkflow 덕분에 시작 부분에 일반적인 nf-core 요약이 생성되었으며, 각 모듈에 대한 줄은 이제 PIPELINE:WORKFLOW:module 전체 이름을 보여줍니다.

### 4.5. 파이프라인 출력 찾기

이제 질문은: 파이프라인의 출력은 어디에 있을까요?
그리고 답은 꽤 흥미롭습니다: 결과를 찾을 수 있는 두 개의 서로 다른 장소가 있습니다.

앞서 기억하실 수도 있듯이, 새로 생성된 workflow의 첫 번째 실행은 다양한 실행 보고서와 메타데이터를 포함하는 `core-hello-results/`라는 디렉토리를 생성했습니다.

```bash
tree core-hello-results
```

??? abstract "디렉토리 내용"

    ```console
    core-hello-results
    └── pipeline_info
        ├── execution_report_2025-11-21_04-47-18.html
        ├── execution_report_2025-11-21_07-29-37.html
        ├── execution_timeline_2025-11-21_04-47-18.html
        ├── execution_timeline_2025-11-21_07-29-37.html
        ├── execution_trace_2025-11-21_04-47-18.txt
        ├── execution_trace_2025-11-21_07-29-37.txt
        ├── hello_software_versions.yml
        ├── params_2025-11-21_04-47-13.json
        ├── params_2025-11-21_07-29-41.json
        └── pipeline_dag_2025-11-21_04-47-18.html
        └── pipeline_dag_2025-11-21_07-29-37.html

    1 directory, 12 files
    ```

첫 번째 실행에서 얻은 것 외에도 추가 실행 보고서 세트를 얻은 것을 볼 수 있습니다. 이때 workflow는 여전히 플레이스홀더일 뿐이었습니다.
이번에는 예상대로 실행된 모든 작업이 표시됩니다.

![Hello 파이프라인의 실행 타임라인 보고서](./img/execution_timeline_hello.png)

!!! note "참고"

    우리가 Github Codespaces에서 최소한의 머신에서 실행하고 있기 때문에 작업이 다시 병렬로 실행되지 않았습니다.
    이들이 병렬로 실행되는 것을 보려면 codespace의 CPU 할당량과 테스트 구성의 리소스 제한을 늘려보세요.

좋습니다만, 실제 파이프라인 결과는 거기에 없습니다!

다음과 같은 상황이 발생했습니다: 모듈 자체를 변경하지 않았으므로 모듈 수준 `publishDir` 지시문에 의해 처리되는 출력은 원래 파이프라인에 지정된 대로 여전히 `results` 디렉토리로 이동합니다.

```bash
tree results
```

??? abstract "디렉토리 내용"

    ```console
    results
    ├── Bonjour-output.txt
    ├── COLLECTED-test-batch-output.txt
    ├── COLLECTED-test-output.txt
    ├── cowpy-COLLECTED-test-batch-output.txt
    ├── cowpy-COLLECTED-test-output.txt
    ├── Hello-output.txt
    ├── Holà-output.txt
    ├── UPPER-Bonjour-output.txt
    ├── UPPER-Hello-output.txt
    └── UPPER-Holà-output.txt

    0 directories, 10 files
    ```

아, 거기 있군요. 원래 Hello 파이프라인의 이전 실행 결과와 섞여 있습니다.

데모 파이프라인의 출력처럼 깔끔하게 구성하려면 출력 게시 방법을 변경해야 합니다.
이 교육 과정의 뒷부분에서 그 방법을 보여드리겠습니다.

<!-- TODO: Hello Nextflow가 workflow 수준 출력을 사용하도록 업데이트되면 이를 업데이트하세요 -->

바로 이것입니다! 원래 파이프라인과 동일한 결과를 얻기 위해 많은 작업처럼 보일 수 있지만, 자동으로 생성되는 멋진 보고서를 모두 얻을 수 있으며, 이제 입력 검증과 나중 섹션에서 다룰 일부 메타데이터 처리 기능을 포함한 nf-core의 추가 기능을 활용할 수 있는 견고한 기반을 갖추게 되었습니다.

---

### 요약

nf-core 템플릿을 사용하여 일반 Nextflow 파이프라인을 nf-core 스타일 파이프라인으로 변환하는 방법을 알게 되었습니다.
그 과정에서 workflow를 구성 가능하게 만드는 방법과, 사용자 정의 nf-core 스타일 파이프라인을 개발할 때 가장 일반적으로 조정해야 하는 nf-core 템플릿 요소를 식별하는 방법을 배웠습니다.

### 다음 단계

잠시 휴식을 취하세요, 힘든 작업이었습니다! 준비가 되면 [3부: nf-core 모듈 사용하기](./03_use_module.md)로 넘어가서 nf-core/modules 저장소에서 커뮤니티가 관리하는 모듈을 활용하는 방법을 배우세요.
