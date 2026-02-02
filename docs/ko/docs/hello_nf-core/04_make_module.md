# 파트 4: nf-core 모듈 만들기

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } AI 지원 번역 - [자세히 알아보기 및 개선 제안](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

Hello nf-core 교육 과정의 네 번째 파트에서는 모듈을 이식 가능하고 유지보수 가능하게 만드는 주요 규칙을 적용하여 nf-core 모듈을 만드는 방법을 보여드립니다.

nf-core 프로젝트는 Part 2에서 워크플로에 사용했던 것과 유사하게 적절하게 구조화된 모듈 템플릿을 자동으로 생성하는 명령(`nf-core modules create`)을 제공합니다.
하지만 교육 목적으로, 먼저 수동으로 작업을 진행할 것입니다: `core-hello` 파이프라인의 로컬 `cowpy` 모듈을 단계별로 nf-core 스타일 모듈로 변환하겠습니다.
그 후에, 향후 더 효율적으로 작업하기 위해 템플릿 기반 모듈 생성을 사용하는 방법을 보여드리겠습니다.

??? info "이 섹션을 시작하는 방법"

    이 섹션은 [Part 3: nf-core 모듈 사용하기](./03_use_module.md)를 완료하고 `CAT_CAT` 모듈을 파이프라인에 통합했다고 가정합니다.

    Part 3을 완료하지 않았거나 이 파트를 새롭게 시작하려는 경우, `core-hello-part3` 솔루션을 시작점으로 사용할 수 있습니다.
    `hello-nf-core/` 디렉토리 내에서 다음 명령을 실행하십시오:

    ```bash
    cp -r solutions/core-hello-part3 core-hello
    cd core-hello
    ```

    이렇게 하면 `CAT_CAT` 모듈이 이미 통합된 파이프라인을 얻게 됩니다.
    다음 명령을 실행하여 성공적으로 실행되는지 테스트할 수 있습니다:

    ```bash
    nextflow run . --outdir core-hello-results -profile test,docker --validate_params false
    ```

---

## 1. `cowpy`를 nf-core 모듈로 변환하기

이 섹션에서는 nf-core 규칙을 `core-hello` 파이프라인의 로컬 `cowpy` 모듈에 적용하여 nf-core 커뮤니티 표준을 따르는 모듈로 변환합니다.

현재 `cowpy` 프로세스 모듈의 코드는 다음과 같습니다:

```groovy title="core-hello/modules/local/cowpy.nf" linenums="1"
#!/usr/bin/env nextflow

// cowpy로 ASCII 아트 생성 (https://github.com/jeffbuttars/cowpy)
process cowpy {

    publishDir 'results', mode: 'copy'

    container 'community.wave.seqera.io/library/cowpy:1.1.5--3db457ae1977a273'
    conda 'conda-forge::cowpy==1.1.5'

    input:
    path input_file
    val character

    output:
    path "cowpy-${input_file}"

    script:
    """
    cat $input_file | cowpy -c "$character" > cowpy-${input_file}
    """
}
```

다음 nf-core 규칙을 단계적으로 적용하겠습니다:

1. **프로세스 이름을 `COWPY`로 대문자화** - 규칙을 따르기 위함입니다.
2. **`COWPY`를 메타데이터 튜플을 사용하도록 업데이트** - 워크플로를 통해 샘플 메타데이터를 전달하기 위함입니다.
3. **`ext.args`로 도구 인수 구성 중앙화** - 인터페이스는 최소화하면서 모듈 다용성을 높이기 위함입니다.
4. **`ext.prefix`로 출력 이름 표준화** - 일관성을 촉진하기 위함입니다.
5. **발행 구성 중앙화** - 일관성을 촉진하기 위함입니다.

각 단계 후에는 모든 것이 예상대로 작동하는지 테스트하기 위해 파이프라인을 실행합니다.

!!! warning "작업 디렉토리"

    이 섹션의 모든 파일 편집과 명령 실행을 위해 `core-hello` 디렉토리(파이프라인 루트)에 있는지 확인하십시오.

    ```bash
    cd core-hello
    ```

### 1.1. 프로세스 이름 대문자화

이것은 순전히 스타일 규칙입니다(기술적 정당성은 없음) 하지만 nf-core 모듈의 표준이므로 준수하겠습니다.

세 가지 변경 사항을 적용해야 합니다:

1. 모듈에서 프로세스 이름 업데이트
2. 워크플로 헤더에서 모듈 import 문 업데이트
3. 워크플로 본문에서 프로세스 호출 및 emit 선언 업데이트

시작하겠습니다!

#### 1.1.1. 모듈에서 프로세스 이름 업데이트

`cowpy.nf` 모듈 파일(`core-hello/modules/local/` 아래)을 열고 프로세스 이름을 대문자로 수정하십시오:

=== "변경 후"

    ```groovy title="core-hello/modules/local/cowpy.nf" linenums="3" hl_lines="2"
    // cowpy로 ASCII 아트 생성 (https://github.com/jeffbuttars/cowpy)
    process COWPY {
    ```

=== "변경 전"

    ```groovy title="core-hello/modules/local/cowpy.nf" linenums="3" hl_lines="2"
    // cowpy로 ASCII 아트 생성 (https://github.com/jeffbuttars/cowpy)
    process cowpy {
    ```

이 경우 대문자화는 완전히 직관적입니다.

프로세스 이름이 여러 단어로 구성된 경우, 예를 들어 원래 camel case로 MyCowpyTool이라는 프로세스가 있었다면, nf-core 규칙은 밑줄을 사용하여 분리하여 MY_COWPY_TOOL이 됩니다.

#### 1.1.2. 모듈 import 문 업데이트

프로세스 이름은 대소문자를 구분하므로, 프로세스 이름을 변경했으므로 `hello.nf`의 워크플로 헤더에서 모듈 import 문을 업데이트해야 합니다:

=== "변경 후"

    ```groovy title="core-hello/workflows/hello.nf" linenums="1" hl_lines="10"
    /*
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    */
    include { paramsSummaryMap       } from 'plugin/nf-schema'
    include { softwareVersionsToYAML } from '../subworkflows/nf-core/utils_nfcore_pipeline'
    include { sayHello               } from '../modules/local/sayHello.nf'
    include { convertToUpper         } from '../modules/local/convertToUpper.nf'
    include { COWPY                  } from '../modules/local/cowpy/main.nf'
    include { CAT_CAT                } from '../modules/nf-core/cat/cat/main'
    ```

=== "변경 전"

    ```groovy title="core-hello/workflows/hello.nf" linenums="1" hl_lines="10"
    /*
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    */
    include { paramsSummaryMap       } from 'plugin/nf-schema'
    include { softwareVersionsToYAML } from '../subworkflows/nf-core/utils_nfcore_pipeline'
    include { sayHello               } from '../modules/local/sayHello.nf'
    include { convertToUpper         } from '../modules/local/convertToUpper.nf'
    include { cowpy                  } from '../modules/local/cowpy/main.nf'
    include { CAT_CAT                } from '../modules/nf-core/cat/cat/main'
    ```

import 문에서 별칭을 사용하여 프로세스 호출을 업데이트하지 않을 수도 있지만, 그렇게 하면 대문자화 규칙을 채택하는 목적이 다소 상실됩니다.

#### 1.1.3. 프로세스 호출 및 emit 선언 업데이트

이제 `hello.nf`의 workflow 블록에서 프로세스에 대한 두 참조를 업데이트하겠습니다:

=== "변경 후"

    ```groovy title="core-hello/workflows/hello.nf" linenums="43" hl_lines="2 17"
    // cowpy로 인사말의 ASCII 아트 생성
    COWPY(CAT_CAT.out.file_out)

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
    cowpy_hellos   = COWPY.out.cowpy_output
    versions       = ch_versions
    ```

=== "변경 전"

    ```groovy title="core-hello/workflows/hello.nf" linenums="43" hl_lines="2 17"
    // cowpy로 인사말의 ASCII 아트 생성
    cowpy(CAT_CAT.out.file_out)

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
    cowpy_hellos   = cowpy.out.cowpy_output
    versions       = ch_versions
    ```

**두 가지** 변경 사항을 모두 적용해야 합니다. 그렇지 않으면 실행할 때 오류가 발생합니다.

#### 1.1.4. 파이프라인을 실행하여 테스트

이러한 변경 후 모든 것이 올바르게 작동하는지 테스트하기 위해 워크플로를 실행해 봅시다.

```bash
nextflow run . --outdir core-hello-results -profile test,docker --validate_params false
```

??? success "명령 출력"

    ```console
     N E X T F L O W   ~  version 25.04.3

    Launching `./main.nf` [elegant_plateau] DSL2 - revision: b9e9b3b8de

    Input/output options
      input                     : /workspaces/training/hello-nf-core/core-hello/assets/greetings.csv
      outdir                    : core-hello-results

    Institutional config options
      config_profile_name       : Test profile
      config_profile_description: Minimal test dataset to check pipeline function

    Generic options
      validate_params           : false
      trace_report_suffix       : 2026-01-06_04-51-29

    Core Nextflow options
      runName                   : elegant_plateau
      containerEngine           : docker
      launchDir                 : /workspaces/training/hello-nf-core/core-hello
      workDir                   : /workspaces/training/hello-nf-core/core-hello/work
      projectDir                : /workspaces/training/hello-nf-core/core-hello
      userName                  : root
      profile                   : test,docker
      configFiles               : /workspaces/training/hello-nf-core/core-hello/nextflow.config

    !! Only displaying parameters that differ from the pipeline defaults !!
    ------------------------------------------------------
    executor >  local (8)
    [7b/66ceb5] CORE_HELLO:HELLO:sayHello (3)       | 3 of 3 ✔
    [8e/1bafb9] CORE_HELLO:HELLO:convertToUpper (3) | 3 of 3 ✔
    [bb/203575] CORE_HELLO:HELLO:CAT_CAT (test)     | 1 of 1 ✔
    [39/715489] CORE_HELLO:HELLO:COWPY              | 1 of 1 ✔
    -[core/hello] Pipeline completed successfully-
    ```

잘 작동합니다! 이제 더 실질적인 변경을 진행하겠습니다.

### 1.2. `COWPY`를 메타데이터 튜플을 사용하도록 업데이트

현재 버전의 `core-hello` 파이프라인에서는 아래 다이어그램의 상단에 표시된 것처럼 `CAT_CAT`의 출력 튜플에서 파일을 추출하여 `COWPY`에 전달하고 있습니다.

<figure class="excalidraw">
    --8<-- "docs/en/docs/hello_nf-core/img/cowpy-inputs.svg"
</figure>

다이어그램의 하단에 표시된 것처럼 `COWPY`가 메타데이터 튜플을 직접 받아들이도록 하여 워크플로를 통해 메타데이터가 계속 흐르도록 하는 것이 더 좋습니다.

이를 위해 다음 변경 사항을 적용해야 합니다:

1. 입력 및 출력 정의 업데이트
2. 워크플로에서 프로세스 호출 업데이트
3. 워크플로에서 emit 블록 업데이트

이 모든 작업을 완료한 후, 모든 것이 이전과 같이 작동하는지 확인하기 위해 파이프라인을 실행합니다.

#### 1.2.1. 입력 및 출력 정의 업데이트

`cowpy.nf` 모듈 파일로 돌아가서 아래와 같이 메타데이터 튜플을 받아들이도록 수정하십시오.

=== "변경 후"

    ```groovy title="core-hello/modules/local/cowpy.nf" linenums="11" hl_lines="2 6"
        input:
            tuple val(meta), path(input_file)
            val character

        output:
            tuple val(meta), path("cowpy-${input_file}"), emit: cowpy_output
    ```

=== "변경 전"

    ```groovy title="core-hello/modules/local/cowpy.nf" linenums="11" hl_lines="2 6"
        input:
            path input_file
            val character

        output:
            path "cowpy-${input_file}"
    ```

보시다시피, **주요 입력**과 **출력** 모두를 Part 3에서 소개한 `tuple val(meta), path(input_file)` 패턴을 따르는 튜플로 변경했습니다.
출력의 경우, 출력 채널에 설명적인 이름을 부여하기 위해 `emit: cowpy_output`을 추가하는 기회를 가졌습니다.

이제 프로세스가 예상하는 것을 변경했으므로, 프로세스 호출에서 제공하는 것을 업데이트해야 합니다.

#### 1.2.2. 워크플로에서 프로세스 호출 업데이트

좋은 소식은 이 변경으로 프로세스 호출이 단순화된다는 것입니다.
이제 `CAT_CAT`의 출력과 `COWPY`의 입력이 동일한 '형태', 즉 둘 다 `tuple val(meta), path(input_file)` 구조로 구성되므로, `CAT_CAT` 프로세스의 출력에서 파일을 명시적으로 추출하는 대신 직접 연결할 수 있습니다.

`hello.nf` 워크플로 파일(`core-hello/workflows/` 아래)을 열고 아래와 같이 `COWPY` 호출을 업데이트하십시오.

=== "변경 후"

    ```groovy title="core-hello/workflows/hello.nf" linenums="43" hl_lines="2"
        // cowpy로 인사말의 ASCII 아트 생성
        COWPY(CAT_CAT.out.file_out, params.character)
    ```

=== "변경 전"

    ```groovy title="core-hello/workflows/hello.nf" linenums="43" hl_lines="1-2 5"
        // cowpy가 아직 메타데이터를 사용하지 않으므로 튜플에서 파일 추출
        ch_for_cowpy = CAT_CAT.out.file_out.map{ meta, file -> file }

        // cowpy로 인사말의 ASCII 아트 생성
        COWPY(ch_for_cowpy, params.character)
    ```

이제 `CAT_CAT.out.file_out`에서 직접 `COWPY`를 호출합니다.

결과적으로 더 이상 `ch_for_cowpy` 채널을 구성할 필요가 없으므로 해당 줄(및 주석 줄)을 완전히 삭제할 수 있습니다.

#### 1.2.3. 워크플로에서 emit 블록 업데이트

`COWPY`가 이제 명명된 출력인 `cowpy_output`을 방출하므로, `hello.nf` 워크플로의 `emit:` 블록을 이를 사용하도록 업데이트할 수 있습니다.

=== "변경 후"

    ```groovy title="core-hello/workflows/hello.nf" linenums="60" hl_lines="2"
        emit:
        cowpy_hellos   = COWPY.out.cowpy_output
        versions       = ch_versions
    ```

=== "변경 전"

    ```groovy title="core-hello/workflows/hello.nf" linenums="60" hl_lines="2"
        emit:
        cowpy_hellos   = COWPY.out
        versions       = ch_versions
    ```

이것은 기술적으로 필요하지 않지만, 가능한 한 명명된 출력을 참조하는 것이 좋은 관행입니다.

#### 1.2.4. 파이프라인을 실행하여 테스트

이러한 변경 후 모든 것이 올바르게 작동하는지 테스트하기 위해 워크플로를 실행해 봅시다.

```bash
nextflow run . --outdir core-hello-results -profile test,docker --validate_params false
```

??? success "명령 출력"

    ```console
     N E X T F L O W   ~  version 25.04.3

    Launching `./main.nf` [modest_saha] DSL2 - revision: b9e9b3b8de

    Downloading plugin nf-schema@2.5.1
    Input/output options
      input                     : /workspaces/training/hello-nf-core/core-hello/assets/greetings.csv
      outdir                    : core-hello-results

    Institutional config options
      config_profile_name       : Test profile
      config_profile_description: Minimal test dataset to check pipeline function

    Generic options
      validate_params           : false
      trace_report_suffix       : 2025-12-27_06-16-55

    Core Nextflow options
      runName                   : modest_saha
      containerEngine           : docker
      launchDir                 : /workspaces/training/hello-nf-core/core-hello
      workDir                   : /workspaces/training/hello-nf-core/core-hello/work
      projectDir                : /workspaces/training/hello-nf-core/core-hello
      userName                  : root
      profile                   : test,docker
      configFiles               : /workspaces/training/hello-nf-core/core-hello/nextflow.config

    !! Only displaying parameters that differ from the pipeline defaults !!
    ------------------------------------------------------
    executor >  local (8)
    [a8/447993] CORE_HELLO:HELLO:sayHello (3)       [100%] 3 of 3 ✔
    [00/1fc59c] CORE_HELLO:HELLO:convertToUpper (3) [100%] 3 of 3 ✔
    [57/ac800d] CORE_HELLO:HELLO:CAT_CAT (test)     [100%] 1 of 1 ✔
    [b7/092f2b] CORE_HELLO:HELLO:COWPY              [100%] 1 of 1 ✔
    -[core/hello] Pipeline completed successfully-
    ```

파이프라인은 성공적으로 실행되어야 하며, 이제 메타데이터가 `CAT_CAT`에서 `COWPY`로 흐릅니다.

이것으로 `COWPY`가 메타데이터 튜플을 처리하도록 하기 위해 필요한 작업이 완료되었습니다.
이제 nf-core 모듈 패턴을 활용하기 위해 할 수 있는 다른 작업을 살펴보겠습니다.

### 1.3. `ext.args`로 도구 인수 구성 중앙화

현재 상태에서 `COWPY` 프로세스는 `character` 매개변수에 대한 값을 받을 것으로 예상합니다.
결과적으로, 도구에서 설정한 기본값으로 만족하더라도 프로세스를 호출할 때마다 값을 제공해야 합니다.
`COWPY`의 경우 이것이 큰 문제는 아니지만, 많은 선택적 매개변수가 있는 도구의 경우 상당히 번거로울 수 있습니다.

nf-core 프로젝트는 구성 파일을 통해 도구 인수를 더 편리하게 관리하기 위해 [`ext.args`](https://www.nextflow.io/docs/latest/reference/process.html#ext)라는 Nextflow 기능을 사용할 것을 권장합니다.

모든 도구 옵션에 대한 프로세스 입력을 선언하는 대신, 명령줄 구성에서 `ext.args`를 참조하도록 모듈을 작성합니다.
그런 다음 모든 모듈에 대한 구성 세부 정보를 통합하는 `modules.config` 파일에서 사용하려는 인수와 값을 보유하도록 `ext.args` 변수를 설정하기만 하면 됩니다.
Nextflow는 런타임에 해당 인수와 값을 도구 명령줄에 추가합니다.

이 접근 방식을 `COWPY` 모듈에 적용해 봅시다.
다음 변경 사항을 적용해야 합니다:

1. `COWPY` 모듈 업데이트
2. `modules.config` 파일에서 `ext.args` 구성
3. `hello.nf` 워크플로 업데이트

이 모든 작업을 완료한 후, 모든 것이 이전과 같이 작동하는지 확인하기 위해 파이프라인을 실행합니다.

#### 1.3.1. `COWPY` 모듈 업데이트

실행해 봅시다.
`cowpy.nf` 모듈 파일(`core-hello/modules/local/` 아래)을 열고 아래와 같이 `ext.args`를 참조하도록 수정하십시오.

=== "변경 후"

    ```groovy title="modules/local/cowpy.nf" linenums="1" hl_lines="18 20"
    #!/usr/bin/env nextflow

    // cowpy로 ASCII 아트 생성 (https://github.com/jeffbuttars/cowpy)
    process COWPY {

        publishDir 'results', mode: 'copy'

        container 'community.wave.seqera.io/library/cowpy:1.1.5--3db457ae1977a273'
        conda 'conda-forge::cowpy==1.1.5'

        input:
            tuple val(meta), path(input_file)

        output:
            tuple val(meta), path("cowpy-${input_file}"), emit: cowpy_output

        script:
        def args = task.ext.args ?: ''
        """
        cat $input_file | cowpy $args > cowpy-${input_file}
        """
    }
    ```

=== "변경 전"

    ```groovy title="core-hello/modules/local/cowpy.nf" linenums="1" hl_lines="13 20"
    #!/usr/bin/env nextflow

    // cowpy로 ASCII 아트 생성 (https://github.com/jeffbuttars/cowpy)
    process COWPY {

        publishDir 'results', mode: 'copy'

        container 'community.wave.seqera.io/library/cowpy:1.1.5--3db457ae1977a273'
        conda 'conda-forge::cowpy==1.1.5'

        input:
            tuple val(meta), path(input_file)
            val character

        output:
            tuple val(meta), path("cowpy-${input_file}"), emit: cowpy_output

        script:
        """
        cat $input_file | cowpy -c "$character" > cowpy-${input_file}
        """
    }
    ```

세 가지 변경 사항을 적용한 것을 볼 수 있습니다.

1. **`input:` 블록에서 `val character` 입력을 제거했습니다.**
   앞으로는 아래에 설명된 대로 `ext.args` 구성을 통해 해당 인수를 제공하겠습니다.

2. **`script:` 블록에서 `def args = task.ext.args ?: ''` 줄을 추가했습니다.**
   이 줄은 `?:` 연산자를 사용하여 `args` 변수의 값을 결정합니다: 비어 있지 않으면 `task.ext.args`의 내용, 비어 있으면 빈 문자열입니다.
   일반적으로 `ext.args`를 참조하지만, 이 코드는 모듈 수준 `ext.args` 구성을 가져오기 위해 `task.ext.args`를 참조해야 합니다.

3. **명령줄에서 `-c "$character"`를 `$args`로 대체했습니다.**
   여기에 Nextflow가 `modules.config` 파일의 `ext.args`에 설정된 모든 도구 인수를 주입합니다.

결과적으로 모듈 인터페이스가 더 간단해졌습니다: 필수 메타데이터와 파일 입력만 예상합니다.

!!! note

    `?:` 연산자는 옆으로 누운 Elvis Presley의 얼굴처럼 보이기 때문에 'Elvis 연산자'라고 불리며, `?` 문자가 그의 머리카락의 웨이브를 상징합니다.

#### 1.3.2. `modules.config` 파일에서 `ext.args` 구성

이제 모듈에서 `character` 선언을 제거했으므로, `modules.config` 구성 파일의 `ext.args`에 추가해야 합니다.

구체적으로, `process {}` 블록에 다음 코드 조각을 추가하겠습니다:

```groovy title="추가할 코드"
withName: 'COWPY' {
    ext.args = { "-c ${params.character}" }
}
```

`withName:` 구문은 이 구성을 `COWPY` 프로세스에만 할당하며, `ext.args = { "-c ${params.character}" }`는 단순히 `character` 매개변수의 값을 포함할 문자열을 구성합니다.
중괄호 사용에 주목하십시오. 이는 Nextflow에게 런타임에 매개변수 값을 평가하도록 지시합니다.

이해가 되시나요? 추가해 봅시다.

`conf/modules.config`를 열고 아래와 같이 `process {}` 블록 내에 구성 코드를 추가하십시오.

=== "변경 후"

    ```groovy title="core-hello/conf/modules.config" linenums="13" hl_lines="8-10"
    process {
        publishDir = [
            path: { "${params.outdir}/${task.process.tokenize(':')[-1].tokenize('_')[0].toLowerCase()}" },
            mode: params.publish_dir_mode,
            saveAs: { filename -> filename.equals('versions.yml') ? null : filename }
        ]

        withName: 'COWPY' {
            ext.args = { "-c ${params.character}" }
        }
    }
    ```

=== "변경 전"

    ```groovy title="core-hello/conf/modules.config" linenums="13"
    process {
        publishDir = [
            path: { "${params.outdir}/${task.process.tokenize(':')[-1].tokenize('_')[0].toLowerCase()}" },
            mode: params.publish_dir_mode,
            saveAs: { filename -> filename.equals('versions.yml') ? null : filename }
        ]
    }
    ```

파이프라인의 모든 모듈이 이 파일에서 `ext.args`를 지정하도록 하는 것을 상상할 수 있기를 바랍니다. 다음과 같은 이점이 있습니다:

- **모듈 인터페이스가 간단하게 유지됩니다** - 필수 메타데이터와 파일 입력만 받습니다
- **파이프라인은 여전히 `params.character`를 노출합니다** - 최종 사용자는 여전히 이전과 같이 구성할 수 있습니다
- **모듈은 이제 이식 가능합니다** - 특정 매개변수 이름을 기대하지 않고 다른 파이프라인에서 재사용할 수 있습니다
- 구성이 `modules.config`에 **중앙화되어** 워크플로 로직을 깔끔하게 유지합니다

모든 파이프라인이 모듈별 구성을 중앙화하는 장소로 `modules.config` 파일을 사용함으로써, 다른 파이프라인에서 모듈을 더 재사용 가능하게 만듭니다.

#### 1.3.3. `hello.nf` 워크플로 업데이트

`COWPY` 모듈이 더 이상 `character` 매개변수를 입력으로 요구하지 않으므로, 워크플로 호출을 그에 따라 업데이트해야 합니다.

`hello.nf` 워크플로 파일(`core-hello/workflows/` 아래)을 열고 아래와 같이 `COWPY` 호출을 업데이트하십시오.

=== "변경 후"

    ```groovy title="core-hello/workflows/hello.nf" linenums="39" hl_lines="2"
        // cowpy로 인사말의 ASCII 아트 생성
        COWPY(CAT_CAT.out.file_out)
    ```

=== "변경 전"

    ```groovy title="core-hello/workflows/hello.nf" linenums="39" hl_lines="2"
        // cowpy로 인사말의 ASCII 아트 생성
        COWPY(CAT_CAT.out.file_out, params.character)
    ```

워크플로 코드가 이제 더 깨끗해졌습니다: 프로세스에 `params.character`를 직접 전달할 필요가 없습니다.
모듈 인터페이스는 최소한으로 유지되어 더 이식 가능하게 만들면서, 파이프라인은 여전히 구성을 통해 명시적 옵션을 제공합니다.

#### 1.3.4. 파이프라인을 실행하여 테스트

워크플로가 여전히 예상대로 작동하는지 테스트해 봅시다. `ext.args` 구성이 작동하는지 확인하기 위해 다른 character를 지정합니다.

더 수수께끼 같은 옵션 중 하나인 `kosh`를 사용하여 이 명령을 실행하십시오:

```bash
nextflow run . --outdir core-hello-results -profile test,docker --validate_params false --character kosh
```

??? success "명령 출력"

    ```console
     N E X T F L O W   ~  version 25.04.3

    Launching `./main.nf` [exotic_planck] DSL2 - revision: b9e9b3b8de

    Input/output options
      input                     : /workspaces/training/hello-nf-core/core-hello/assets/greetings.csv
      outdir                    : core-hello-results

    Institutional config options
      config_profile_name       : Test profile
      config_profile_description: Minimal test dataset to check pipeline function

    Generic options
      validate_params           : false
      trace_report_suffix       : 2025-12-27_06-23-13

    Core Nextflow options
      runName                   : exotic_planck
      containerEngine           : docker
      launchDir                 : /workspaces/training/hello-nf-core/core-hello
      workDir                   : /workspaces/training/hello-nf-core/core-hello/work
      projectDir                : /workspaces/training/hello-nf-core/core-hello
      userName                  : root
      profile                   : test,docker
      configFiles               : /workspaces/training/hello-nf-core/core-hello/nextflow.config

    !! Only displaying parameters that differ from the pipeline defaults !!
    ------------------------------------------------------
    executor >  local (8)
    [13/9e3c0e] CORE_HELLO:HELLO:sayHello (3)       [100%] 3 of 3 ✔
    [e2/5b0ee5] CORE_HELLO:HELLO:convertToUpper (3) [100%] 3 of 3 ✔
    [b6/4fb569] CORE_HELLO:HELLO:CAT_CAT (test)     [100%] 1 of 1 ✔
    [38/eb29ea] CORE_HELLO:HELLO:COWPY              [100%] 1 of 1 ✔
    -[core/hello] Pipeline completed successfully-
    ```

이전처럼 성공적으로 실행되어야 합니다.

`ext.args` 구성이 작동했는지 출력을 확인하여 검증해 봅시다.
파일 브라우저에서 출력을 찾거나 작업 해시(위 예제의 `38/eb29ea` 부분)를 사용하여 출력 파일을 확인하십시오:

```bash
cat work/38/eb29ea*/cowpy-test.txt
```

??? success "명령 출력"

    ```console
    _________
    / HELLO   \
    | HOLà    |
    \ BONJOUR /
    ---------
        \
        \
          \
      ___       _____     ___
    /   \     /    /|   /   \
    |     |   /    / |  |     |
    |     |  /____/  |  |     |
    |     |  |    |  |  |     |
    |     |  | {} | /   |     |
    |     |  |____|/    |     |
    |     |    |==|     |     |
    |      \___________/      |
    |                         |
    |                         |
    ```

`kosh` character로 표시된 ASCII 아트를 확인할 수 있어야 하며, 이는 `ext.args` 구성이 작동했음을 확인합니다!

??? info "(선택 사항) 명령 파일 검사"

    구성이 정확히 어떻게 적용되었는지 확인하려면, `.command.sh` 파일을 검사할 수 있습니다:

    ```bash
    cat work/38/eb29ea*/.command.sh
    ```

    `-c kosh` 인수가 포함된 `cowpy` 명령을 볼 수 있습니다:

    ```console
    #!/usr/bin/env bash
    ...
    cat test.txt | cowpy -c kosh > cowpy-test.txt
    ```

    이는 `.command.sh` 파일이 `ext.args` 구성에 따라 올바르게 생성되었음을 보여줍니다.

여기서 우리가 달성한 것에 대해 생각해 보십시오.
이 접근 방식은 모듈 인터페이스를 필수 데이터(파일, 메타데이터 및 모든 필수 샘플별 매개변수)에 집중하도록 유지하면서, 도구의 동작을 제어하는 옵션은 구성을 통해 별도로 처리합니다.

이것은 `cowpy`와 같은 간단한 도구에는 불필요해 보일 수 있지만, 많은 선택적 인수가 있는 데이터 분석 도구의 경우 큰 차이를 만들 수 있습니다.

이 접근 방식의 이점을 요약하면:

- **깨끗한 인터페이스**: 모듈은 필수 데이터 입력(메타데이터 및 파일)에 집중합니다
- **유연성**: 사용자는 샘플별 값을 포함하여 구성을 통해 도구 인수를 지정할 수 있습니다
- **일관성**: 모든 nf-core 모듈이 이 패턴을 따릅니다
- **이식성**: 하드코딩된 도구 옵션 없이 모듈을 재사용할 수 있습니다
- **워크플로 변경 없음**: 도구 옵션을 추가하거나 변경해도 워크플로 코드를 업데이트할 필요가 없습니다

!!! note

    `ext.args` 시스템은 메타데이터에 따라 인수 값을 동적으로 전환하는 것을 포함하여 여기서 다루지 않은 강력한 추가 기능이 있습니다. 자세한 내용은 [nf-core 모듈 사양](https://nf-co.re/docs/guidelines/components/modules)을 참조하십시오.

### 1.4. `ext.prefix`로 출력 이름 표준화

이제 `COWPY` 프로세스가 메타맵에 액세스할 수 있게 되었으므로, 또 다른 유용한 nf-core 패턴을 활용할 수 있습니다: 메타데이터를 기반으로 출력 파일 이름을 지정하는 것입니다.

여기서는 `meta.id`(메타맵에 포함된 식별자)를 사용하여 모듈 전체에서 출력 파일 이름을 표준화하면서도 필요한 경우 개별적으로 모듈을 구성할 수 있도록 하는 `ext.prefix`라는 Nextflow 기능을 사용하겠습니다.

이것은 `ext.args`와 유사하지만, 진행하면서 자세히 설명할 몇 가지 차이점이 있습니다.

이 접근 방식을 `COWPY` 모듈에 적용해 봅시다.
다음 변경 사항을 적용해야 합니다:

1. `COWPY` 모듈 업데이트
2. `modules.config` 파일에서 `ext.prefix` 구성

(워크플로에는 변경이 필요하지 않습니다.)

이 작업을 완료한 후, 모든 것이 이전과 같이 작동하는지 확인하기 위해 파이프라인을 실행합니다.

#### 1.4.1. `COWPY` 모듈 업데이트

`cowpy.nf` 모듈 파일(`core-hello/modules/local/` 아래)을 열고 아래와 같이 `ext.prefix`를 참조하도록 수정하십시오.

=== "변경 후"

    ```groovy title="core-hello/modules/local/cowpy.nf" linenums="1" hl_lines="2 6 8"
        output:
            tuple val(meta), path("${prefix}.txt"), emit: cowpy_output

        script:
        def args = task.ext.args ?: ''
        prefix = task.ext.prefix ?: "${meta.id}"
        """
        cat $input_file | cowpy $args > ${prefix}.txt
        """
    }
    ```

=== "변경 전"

    ```groovy title="core-hello/modules/local/cowpy.nf" linenums="1" hl_lines="2 7"
        output:
            tuple val(meta), path("cowpy-${input_file}"), emit: cowpy_output

        script:
        def args = task.ext.args ?: ''
        """
        cat $input_file | cowpy $args > cowpy-${input_file}
        """
    }
    ```

세 가지 변경 사항을 적용한 것을 볼 수 있습니다.

1. **`script:` 블록에서 `prefix = task.ext.prefix ?: "${meta.id}"` 줄을 추가했습니다.**
   이 줄은 `?:` 연산자를 사용하여 `prefix` 변수의 값을 결정합니다: 비어 있지 않으면 `task.ext.prefix`의 내용, 비어 있으면 메타맵의 식별자(`meta.id`)입니다.
   일반적으로 `ext.prefix`를 참조하지만, 이 코드는 모듈 수준 `ext.prefix` 구성을 가져오기 위해 `task.ext.prefix`를 참조해야 합니다.

2. **명령줄에서 `cowpy-${input_file}`을 `${prefix}.txt`로 대체했습니다.**
   여기에 Nextflow가 위 줄에서 결정된 `prefix` 값을 주입합니다.

3. **`output:` 블록에서 `path("cowpy-${input_file}")`을 `path("${prefix}.txt")`로 대체했습니다.**
   이것은 단순히 명령줄에 작성된 내용에 따라 파일 경로가 무엇인지 반복합니다.

결과적으로 출력 파일 이름은 이제 적절한 파일 형식 확장자와 결합된 합리적인 기본값(메타맵의 식별자)을 사용하여 구성됩니다.

#### 1.4.2. `modules.config` 파일에서 `ext.prefix` 구성

이 경우 합리적인 기본값은 우리의 취향에 충분히 표현적이지 않습니다. 이전과 같이 도구 이름을 포함하는 사용자 정의 명명 패턴인 `cowpy-<id>.txt`를 사용하고 싶습니다.

`ext.args`로 `character` 매개변수를 처리했던 것처럼 `modules.config`에서 `ext.prefix`를 구성하여 이를 수행하겠습니다. 단, 이번에는 `withName: 'COWPY' {}` 블록이 이미 존재하므로 다음 줄만 추가하면 됩니다:

```groovy title="추가할 코드"
ext.prefix = { "cowpy-${meta.id}" }
```

이것은 원하는 문자열을 구성합니다.
다시 한 번 중괄호를 사용하여 Nextflow에게 런타임에 `meta.id` 값을 평가하도록 지시합니다.

추가해 봅시다.

`conf/modules.config`를 열고 아래와 같이 `process {}` 블록 내에 구성 코드를 추가하십시오.

=== "변경 후"

    ```groovy title="core-hello/conf/modules.config" linenums="21" hl_lines="3"
        withName: 'COWPY' {
            ext.args = { "-c ${params.character}" }
            ext.prefix = { "cowpy-${meta.id}" }
        }
    ```

=== "변경 전"

    ```groovy title="core-hello/conf/modules.config" linenums="21"
        withName: 'COWPY' {
            ext.args = { "-c ${params.character}" }
        }
    ```

궁금하실 수도 있는데, `ext.prefix` 클로저는 구성이 메타데이터를 사용할 수 있는 프로세스 실행 컨텍스트에서 평가되기 때문에 올바른 메타데이터 조각에 액세스할 수 있습니다.

#### 1.4.3. 파이프라인을 실행하여 테스트

워크플로가 여전히 예상대로 작동하는지 테스트해 봅시다.

```bash
nextflow run . --outdir core-hello-results -profile test,docker --validate_params false
```

??? success "명령 출력"

    ```console
     N E X T F L O W   ~  version 25.04.3

    Launching `./main.nf` [admiring_turing] DSL2 - revision: b9e9b3b8de

    Input/output options
      input                     : /workspaces/training/hello-nf-core/core-hello/assets/greetings.csv
      outdir                    : core-hello-results

    Institutional config options
      config_profile_name       : Test profile
      config_profile_description: Minimal test dataset to check pipeline function

    Generic options
      validate_params           : false
      trace_report_suffix       : 2025-12-27_06-29-02

    Core Nextflow options
      runName                   : admiring_turing
      containerEngine           : docker
      launchDir                 : /workspaces/training/hello-nf-core/core-hello
      workDir                   : /workspaces/training/hello-nf-core/core-hello/work
      projectDir                : /workspaces/training/hello-nf-core/core-hello
      userName                  : root
      profile                   : test,docker
      configFiles               : /workspaces/training/hello-nf-core/core-hello/nextflow.config

    !! Only displaying parameters that differ from the pipeline defaults !!
    ------------------------------------------------------
    executor >  local (8)
    [b2/e08524] CORE_HELLO:HELLO:sayHello (3)       [100%] 3 of 3 ✔
    [13/88939f] CORE_HELLO:HELLO:convertToUpper (3) [100%] 3 of 3 ✔
    [23/4554e1] CORE_HELLO:HELLO:CAT_CAT (test)     [100%] 1 of 1 ✔
    [a3/c6cbe9] CORE_HELLO:HELLO:COWPY              [100%] 1 of 1 ✔
    -[core/hello] Pipeline completed successfully-
    ```

결과 디렉토리의 출력을 살펴보십시오.
기본 배치 이름을 기반으로 이전과 동일한 이름을 가진 cowpy 출력 파일 `cowpy-test.txt`를 볼 수 있어야 합니다.

??? abstract "디렉토리 내용"

    ```console hl_lines="3"
    results
    ├── Bonjour-output.txt
    ├── cowpy-test.txt
    ├── Hello-output.txt
    ├── Holà-output.txt
    ├── UPPER-Bonjour-output.txt
    ├── UPPER-Hello-output.txt
    └── UPPER-Holà-output.txt
    ```

모듈이나 워크플로 코드를 변경하지 않고도 명명 패턴을 변경할 수 있다는 것을 확인하기 위해 `conf/modules.config`에서 `ext.prefix` 구성을 자유롭게 변경해 보십시오.

또는 명령줄에서 다른 `--batch` 매개변수를 지정하여 다시 실행하여 해당 부분이 여전히 즉시 사용자 정의 가능한지 확인할 수도 있습니다.

이것은 `ext.prefix`를 통해 모듈 인터페이스를 유연하게 유지하면서 선호하는 명명 규칙을 유지할 수 있는 방법을 보여줍니다.

이 접근 방식의 이점을 요약하면:

- **표준화된 명명**: 출력 파일은 일반적으로 메타데이터의 샘플 ID를 사용하여 명명됩니다
- **구성 가능**: 필요한 경우 사용자가 기본 명명을 재정의할 수 있습니다
- **일관성**: 모든 nf-core 모듈이 이 패턴을 따릅니다
- **예측 가능**: 출력 파일이 무엇으로 호출될지 쉽게 알 수 있습니다

꽤 좋죠?
그런데 nf-core 가이드라인에 맞게 모듈을 개선하기 위해 수행해야 할 중요한 변경 사항이 하나 더 있습니다.

### 1.5. 발행 구성 중앙화

이제 두 개의 서로 다른 디렉토리에 출력을 발행하고 있다는 것을 알아차렸을 것입니다:

- **`results`** — 로컬 모듈에 대해 처음부터 사용해 온 원래 출력 디렉토리로, 모듈별 `publishDir` 지시문을 사용하여 개별적으로 설정됨;
- **`core-hello-results`** — 명령줄에서 `--outdir`로 설정된 출력 디렉토리로, nf-core 로그와 `CAT_CAT`이 발행한 결과를 받고 있음.

이것은 지저분하고 최적이 아닙니다. 모든 것에 대해 하나의 위치를 갖는 것이 더 좋을 것입니다.
물론, 각 로컬 모듈로 이동하여 `publishDir` 지시문을 수동으로 업데이트하여 `core-hello-results` 디렉토리를 사용하도록 할 수 있지만, 다음에 출력 디렉토리를 변경하기로 결정하면 어떻게 될까요?

개별 모듈이 발행 결정을 내리도록 하는 것은 명백히 올바른 방법이 아닙니다. 특히 동일한 모듈이 다양한 요구 사항이나 선호도를 가진 사람들에 의해 많은 다른 파이프라인에서 사용될 수 있는 세상에서는 말입니다.
워크플로 구성 수준에서 출력이 발행되는 위치를 제어할 수 있어야 합니다.

"잠깐," 여러분이 말할 수 있습니다. "`CAT_CAT`은 `--outdir`로 출력을 보내고 있습니다. `publishDir` 지시문을 복사해야 할까요?"

네, 훌륭한 아이디어입니다.

단, `publishDir` 지시문이 없습니다. (계속해서 모듈 코드를 보십시오.)

nf-core 파이프라인은 개별 모듈에서 `publishDir`을 구성하는 대신 `conf/modules.config`에서 구성하여 워크플로 수준에서 제어를 중앙화하기 때문입니다.
구체적으로, nf-core 템플릿은 재정의 지시문이 제공되지 않는 한 모든 모듈에 적용되는 기본 `publishDir` 지시문(미리 정의된 디렉토리 구조 포함)을 선언합니다.

멋지게 들리지 않나요? 이 기본 지시문을 활용하기 위해 로컬 모듈에서 현재 `publishDir` 지시문을 제거하기만 하면 될까요?

`COWPY`에서 시도하여 어떤 일이 일어나는지 확인한 다음, 작동 방식을 이해하기 위해 기본 구성 코드를 살펴보겠습니다.

마지막으로, 원하는 경우 기본 동작을 재정의하는 방법을 보여드리겠습니다.

#### 1.5.1. `COWPY`에서 `publishDir` 지시문 제거

실행해 봅시다.
`cowpy.nf` 모듈 파일(`core-hello/modules/local/` 아래)을 열고 아래와 같이 `publishDir` 지시문을 제거하십시오.

=== "변경 후"

    ```groovy title="core-hello/modules/local/cowpy.nf (발췌)" linenums="1"
    #!/usr/bin/env nextflow

    // cowpy로 ASCII 아트 생성 (https://github.com/jeffbuttars/cowpy)
    process COWPY {

        container 'community.wave.seqera.io/library/cowpy:1.1.5--3db457ae1977a273'
        conda 'conda-forge::cowpy==1.1.5'
    ```

=== "변경 전"

    ```groovy title="core-hello/modules/local/cowpy.nf (발췌)" linenums="1" hl_lines="6"
    #!/usr/bin/env nextflow

    // cowpy로 ASCII 아트 생성 (https://github.com/jeffbuttars/cowpy)
    process COWPY {

        publishDir 'results', mode: 'copy'

        container 'community.wave.seqera.io/library/cowpy:1.1.5--3db457ae1977a273'
        conda 'conda-forge::cowpy==1.1.5'

    ```

이게 전부입니다!

#### 1.5.2. 파이프라인을 실행하여 테스트

이제 파이프라인을 실행하면 어떤 일이 일어나는지 살펴봅시다.

```bash
nextflow run . --outdir core-hello-results -profile test,docker --validate_params false
```

??? success "명령 출력"

    ```console
     N E X T F L O W   ~  version 25.04.3

    Launching `./main.nf` [silly_caravaggio] DSL2 - revision: b9e9b3b8de

    Input/output options
      input                     : /workspaces/training/hello-nf-core/core-hello/assets/greetings.csv
      outdir                    : core-hello-results

    Institutional config options
      config_profile_name       : Test profile
      config_profile_description: Minimal test dataset to check pipeline function

    Generic options
      validate_params           : false
      trace_report_suffix       : 2025-12-27_06-35-56

    Core Nextflow options
      runName                   : silly_caravaggio
      containerEngine           : docker
      launchDir                 : /workspaces/training/hello-nf-core/core-hello
      workDir                   : /worksp
