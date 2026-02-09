# 파트 3: nf-core 모듈 사용하기

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } AI 기반 번역 - [자세히 알아보고 개선 사항 제안하기](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

Hello nf-core 교육 과정의 세 번째 파트에서는 기존 nf-core 모듈을 찾아서 설치하고 파이프라인에서 사용하는 방법을 보여드립니다.

nf-core로 작업할 때 얻을 수 있는 큰 장점 중 하나는 [nf-core/modules](https://github.com/nf-core/modules) 저장소에서 미리 구축되고 테스트된 모듈을 활용할 수 있다는 것입니다.
모든 프로세스를 처음부터 작성하는 대신, 모범 사례를 따르는 커뮤니티 유지 관리 모듈을 설치하고 사용할 수 있습니다.

이것이 어떻게 작동하는지 보여드리기 위해, `core-hello` 파이프라인에서 커스텀 `collectGreetings` 모듈을 nf-core/modules의 `cat/cat` 모듈로 교체하겠습니다.

??? info "이 섹션을 시작하는 방법"

    이 과정 섹션은 [파트 2: nf-core를 위한 Hello 재작성](./02_rewrite_hello.md)을 완료하고 작동하는 `core-hello` 파이프라인이 있다고 가정합니다.

    파트 2를 완료하지 않았거나 이 파트를 새로 시작하고 싶다면, `core-hello-part2` 해결책을 시작점으로 사용할 수 있습니다.
    `hello-nf-core/` 디렉토리 내에서 다음 명령을 실행하세요:

    ```bash
    cp -r solutions/core-hello-part2 core-hello
    cd core-hello
    ```

    이렇게 하면 모듈을 추가할 준비가 된 완전히 작동하는 nf-core 파이프라인을 얻을 수 있습니다.
    다음 명령을 실행하여 성공적으로 실행되는지 테스트할 수 있습니다:

    ```bash
    nextflow run . --outdir core-hello-results -profile test,docker --validate_params false
    ```

---

## 1. 적합한 nf-core 모듈 찾기 및 설치하기

먼저, 기존 nf-core 모듈을 찾아서 파이프라인에 설치하는 방법을 학습하겠습니다.

여러 인사말 파일을 하나로 연결하기 위해 Unix `cat` 명령을 사용하는 `collectGreetings` 프로세스를 교체하는 것을 목표로 합니다.
파일 연결은 매우 일반적인 작업이므로, 이미 nf-core에 그 목적을 위해 설계된 모듈이 있을 것이라고 추론할 수 있습니다.

시작해 보겠습니다.

### 1.1. nf-core 웹사이트에서 사용 가능한 모듈 탐색하기

nf-core 프로젝트는 [https://nf-co.re/modules](https://nf-co.re/modules)에서 모듈의 중앙 집중식 카탈로그를 유지 관리합니다.

웹 브라우저에서 모듈 페이지로 이동하여 검색창을 사용해 'concatenate'를 검색하세요.

![모듈 검색 결과](./img/module-search-results.png)

보시다시피 상당히 많은 결과가 있으며, 그 중 많은 것들이 매우 특정한 유형의 파일을 연결하도록 설계된 모듈입니다.
그 중에서 범용 목적인 `cat_cat`이라는 모듈을 볼 수 있을 것입니다.

!!! note "모듈 명명 규칙"

    밑줄(`_`)은 모듈 이름에서 슬래시(`/`) 문자를 대신하여 사용됩니다.

    nf-core 모듈은 도구가 여러 명령을 제공할 때 `software/command` 명명 규칙을 따릅니다. 예를 들어 `samtools/view`(samtools 패키지, view 명령) 또는 `gatk/haplotypecaller`(GATK 패키지, HaplotypeCaller 명령)와 같습니다.
    하나의 주요 명령만 제공하는 도구의 경우, 모듈은 `fastqc` 또는 `multiqc`와 같이 단일 레벨을 사용합니다.

`cat_cat` 모듈 박스를 클릭하여 모듈 문서를 확인하세요.

모듈 페이지는 다음을 보여줍니다:

- 간단한 설명: "A module for concatenation of gzipped or uncompressed files"
- 설치 명령: `nf-core modules install cat/cat`
- 입력 및 출력 채널 구조
- 사용 가능한 매개변수

### 1.2. 명령줄에서 사용 가능한 모듈 나열하기

또는 nf-core 도구를 사용하여 명령줄에서 직접 모듈을 검색할 수도 있습니다.

```bash
nf-core modules list remote
```

이렇게 하면 nf-core/modules 저장소에서 사용 가능한 모든 모듈 목록이 표시되지만, 찾고 있는 모듈의 이름을 이미 알지 못한다면 조금 덜 편리합니다.
하지만 알고 있다면, 목록을 `grep`으로 파이프하여 특정 모듈을 찾을 수 있습니다:

```bash
nf-core modules list remote | grep 'cat/cat'
```

??? success "명령 출력"

    ```console
    │ cat/cat
    ```

`grep` 접근 방식은 이름에 검색어가 있는 결과만 추출한다는 점을 기억하세요. 이는 `cat_cat`에는 작동하지 않습니다.

### 1.3. 모듈에 대한 자세한 정보 얻기

명령줄에서 특정 모듈에 대한 자세한 정보를 보려면 `info` 명령을 사용하세요:

```bash
nf-core modules info cat/cat
```

이렇게 하면 입력, 출력 및 기본 사용 정보를 포함한 모듈에 대한 문서가 표시됩니다.

??? success "명령 출력"

    ```console

                                              ,--./,-.
              ___     __   __   __   ___     /,-._.--~\
        |\ | |__  __ /  ` /  \ |__) |__         }  {
        | \| |       \__, \__/ |  \ |___     \`-._,-`-,
                                              `._,._,'

        nf-core/tools version 3.4.1 - https://nf-co.re


    ╭─ Module: cat/cat  ─────────────────────────────────────────────────╮
    │ 🌐 Repository: https://github.com/nf-core/modules.git              │
    │ 🔧 Tools: cat                                                      │
    │ 📖 Description: A module for concatenation of gzipped or           │
    │ uncompressed files                                                 │
    ╰────────────────────────────────────────────────────────────────────╯
                      ╷                                          ╷
    📥 Inputs        │Description                               │Pattern
    ╺━━━━━━━━━━━━━━━━━┿━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┿━━━━━━━╸
    input[0]         │                                          │
    ╶─────────────────┼──────────────────────────────────────────┼───────╴
      meta  (map)     │Groovy Map containing sample information  │
                      │e.g. [ id:'test', single_end:false ]      │
    ╶─────────────────┼──────────────────────────────────────────┼───────╴
      files_in  (file)│List of compressed / uncompressed files   │      *
                      ╵                                          ╵
                          ╷                                 ╷
    📥 Outputs           │Description                      │     Pattern
    ╺━━━━━━━━━━━━━━━━━━━━━┿━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┿━━━━━━━━━━━━╸
    file_out             │                                 │
    ╶─────────────────────┼─────────────────────────────────┼────────────╴
      meta  (map)         │Groovy Map containing sample     │
                          │information                      │
    ╶─────────────────────┼─────────────────────────────────┼────────────╴
      ${prefix}  (file)   │Concatenated file. Will be       │ ${file_out}
                          │gzipped if file_out ends with    │
                          │".gz"                            │
    ╶─────────────────────┼─────────────────────────────────┼────────────╴
    versions             │                                 │
    ╶─────────────────────┼─────────────────────────────────┼────────────╴
      versions.yml  (file)│File containing software versions│versions.yml
                          ╵                                 ╵

    💻  Installation command: nf-core modules install cat/cat

    ```

이것은 웹사이트에서 찾을 수 있는 것과 정확히 동일한 정보입니다.

### 1.4. cat/cat 모듈 설치하기

이제 원하는 모듈을 찾았으므로, 파이프라인의 소스 코드에 추가해야 합니다.

좋은 소식은 nf-core 프로젝트에 이 부분을 쉽게 만드는 도구가 포함되어 있다는 것입니다.
특히, `nf-core modules install` 명령을 사용하면 코드를 검색하고 프로젝트에서 사용할 수 있도록 만드는 작업을 한 단계로 자동화할 수 있습니다.

파이프라인 디렉토리로 이동하여 설치 명령을 실행하세요:

```bash
cd core-hello
nf-core modules install cat/cat
```

도구가 먼저 저장소 유형을 지정하라는 메시지를 표시할 수 있습니다.
(그렇지 않다면, "마지막으로, 도구가 모듈 설치를 진행합니다."로 건너뛰세요.)

??? success "명령 출력"

    ```console

                                          ,--./,-.
          ___     __   __   __   ___     /,-._.--~\
    |\ | |__  __ /  ` /  \ |__) |__         }  {
    | \| |       \__, \__/ |  \ |___     \`-._,-`-,
                                          `._,._,'

    nf-core/tools version 3.4.1 - https://nf-co.re


    WARNING  'repository_type' not defined in .nf-core.yml
    ? Is this repository a pipeline or a modules repository? (Use arrow keys)
    » Pipeline
      Modules repository
    ```

그렇다면, Enter 키를 눌러 기본 응답(`Pipeline`)을 수락하고 계속하세요.

그런 다음 도구는 향후 이 프롬프트를 피하기 위해 프로젝트 구성을 수정할 것을 제안합니다.

??? success "명령 출력"

    ```console
        INFO     To avoid this prompt in the future, add the 'repository_type' key to your .nf-core.yml file.
        ? Would you like me to add this config now? [y/n] (y):
    ```

이 편리한 도구를 활용하는 것이 좋습니다!
Enter 키를 눌러 기본 응답(yes)을 수락하세요.

마지막으로, 도구가 모듈 설치를 진행합니다.

??? success "명령 출력"

    ```console
    INFO Config added to '.nf-core.yml'
    INFO Reinstalling modules found in 'modules.json' but missing from directory:
    INFO Installing 'cat/cat'
    INFO Use the following statement to include this module:

        include { CAT_CAT } from '../modules/nf-core/cat/cat/main'
    ```

명령은 자동으로:

- 모듈 파일을 `modules/nf-core/cat/cat/`에 다운로드합니다
- 설치된 모듈을 추적하기 위해 `modules.json`을 업데이트합니다
- 워크플로우에서 사용할 올바른 `include` 문을 제공합니다

!!! tip

    모듈 설치 명령을 실행하기 전에 현재 작업 디렉토리가 파이프라인 프로젝트의 루트인지 항상 확인하세요.

모듈이 올바르게 설치되었는지 확인해 보겠습니다:

```bash
tree -L 4 modules
```

??? abstract "디렉토리 내용"

    ```console
    modules
    ├── local
    │   ├── collectGreetings.nf
    │   ├── convertToUpper.nf
    │   ├── cowpy.nf
    │   └── sayHello.nf
    └── nf-core
        └── cat
            └── cat
                ├── environment.yml
                ├── main.nf
                ├── meta.yml
                └── tests

    5 directories, 7 files
    ```

nf-core 유틸리티에 로컬에 설치된 모듈을 나열하도록 요청하여 설치를 확인할 수도 있습니다:

```bash
nf-core modules list local
```

??? success "명령 출력"

    ```console
    INFO     Repository type: pipeline
    INFO     Modules installed in '.':

    ┏━━━━━━━━━━━━━┳━━━━━━━━━━━━━━━━━┳━━━━━━━━━━━━━┳━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┳━━━━━━━━━━━━┓
    ┃ Module Name ┃ Repository      ┃ Version SHA ┃ Message                                ┃ Date       ┃
    ┡━━━━━━━━━━━━━╇━━━━━━━━━━━━━━━━━╇━━━━━━━━━━━━━╇━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━╇━━━━━━━━━━━━┩
    │ cat/cat     │ nf-core/modules │ 41dfa3f     │ update meta.yml of all modules (#8747) │ 2025-07-07 │
    └─────────────┴─────────────────┴─────────────┴────────────────────────────────────────┴────────────┘
    ```

이것은 `cat/cat` 모듈이 이제 프로젝트의 소스 코드의 일부임을 확인합니다.

하지만 새 모듈을 실제로 사용하려면, 파이프라인으로 가져와야 합니다.

### 1.5. 모듈 가져오기 업데이트하기

`workflows/hello.nf` 워크플로우의 가져오기 섹션에서 `collectGreetings` 모듈에 대한 `include` 문을 `CAT_CAT`에 대한 것으로 교체하겠습니다.

참고로, 모듈 설치 도구가 사용할 정확한 문을 제공했습니다:

```groovy title="설치 명령으로 생성된 가져오기 문"
include { CAT_CAT } from '../modules/nf-core/cat/cat/main'`
```

nf-core 규칙은 모듈을 가져올 때 모듈 이름에 대문자를 사용하는 것입니다.

[core-hello/workflows/hello.nf](core-hello/workflows/hello.nf)를 열고 다음과 같이 교체하세요:

=== "후"

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
    include { CAT_CAT                } from '../modules/nf-core/cat/cat/main'
    include { cowpy                  } from '../modules/local/cowpy.nf'
    ```

=== "전"

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
    include { collectGreetings       } from '../modules/local/collectGreetings.nf'
    include { cowpy                  } from '../modules/local/cowpy.nf'
    ```

nf-core 모듈의 경로가 로컬 모듈과 어떻게 다른지 주목하세요:

- **nf-core 모듈**: `'../modules/nf-core/cat/cat/main'` (`main.nf` 참조)
- **로컬 모듈**: `'../modules/local/collectGreetings.nf'` (단일 파일 참조)

이제 모듈을 워크플로우에서 사용할 수 있으므로, `collectGreetings` 호출을 `CAT_CAT`을 사용하도록 교체하기만 하면 됩니다. 맞죠?

그렇게 빠르지 않습니다.

이 시점에서 코드 편집을 시작하고 싶을 수 있지만, 새 모듈이 무엇을 기대하고 무엇을 생성하는지 신중하게 검토하는 것이 좋습니다.

우리는 아직 다루지 않은 새로운 메커니즘인 메타데이터 맵을 포함하기 때문에 이것을 별도의 섹션으로 다루겠습니다.

!!! note

    선택적으로 `collectGreetings.nf` 파일을 삭제할 수 있습니다:

    ```bash
    rm modules/local/collectGreetings.nf
    ```

    하지만 로컬 모듈과 nf-core 모듈 간의 차이점을 이해하기 위한 참조로 유지하고 싶을 수 있습니다.

### 핵심 정리

nf-core 모듈을 찾아서 프로젝트에서 사용할 수 있도록 만드는 방법을 알게 되었습니다.

### 다음 단계

새 모듈이 요구하는 것을 평가하고 파이프라인에 통합하기 위해 필요한 중요한 변경 사항을 식별합니다.

---

## 2. 새 모듈의 요구 사항 평가하기

구체적으로, 모듈의 **인터페이스**, 즉 입력 및 출력 정의를 검토하고 교체하려는 모듈의 인터페이스와 비교해야 합니다.
이를 통해 새 모듈을 단순히 드롭인 교체로 처리할 수 있는지 또는 일부 배선을 조정해야 하는지 결정할 수 있습니다.

이상적으로는 모듈을 설치하기 _전에_ 이 작업을 수행해야 하지만, 늦었더라도 하는 것이 낫습니다.
(참고로, 더 이상 원하지 않는 모듈을 제거하는 `uninstall` 명령이 있습니다.)

!!! note

    CAT_CAT 프로세스에는 여기서 보여드리려는 것과 엄밀히 관련이 없는 다양한 압축 유형, 파일 확장자 등에 대한 상당히 영리한 처리가 포함되어 있으므로, 대부분을 무시하고 중요한 부분에만 집중하겠습니다.

### 2.1. 두 모듈의 인터페이스 비교하기

참고로, 우리의 `collectGreetings` 모듈 인터페이스는 다음과 같습니다:

```groovy title="modules/local/collectGreetings.nf (발췌)" linenums="1" hl_lines="6-7 10"
process collectGreetings {

    publishDir 'results', mode: 'copy'

    input:
        path input_files
        val batch_name

    output:
        path "COLLECTED-${batch_name}-output.txt" , emit: outfile
```

`collectGreetings` 모듈은 두 개의 입력을 받습니다:

- `input_files`는 처리할 하나 이상의 입력 파일을 포함합니다;
- `batch_name`은 출력 파일에 실행별 이름을 할당하는 데 사용하는 값으로, 메타데이터의 한 형태입니다.

완료되면, `collectGreetings`는 `outfile` 태그로 내보내진 단일 파일 경로를 출력합니다.

이에 비해, `cat/cat` 모듈의 인터페이스는 더 복잡합니다:

```groovy title="modules/nf-core/cat/cat/main.nf (발췌)" linenums="1" hl_lines="11 14"
process CAT_CAT {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/pigz:2.3.4' :
        'biocontainers/pigz:2.3.4' }"

    input:
    tuple val(meta), path(files_in)

    output:
    tuple val(meta), path("${prefix}"), emit: file_out
    path "versions.yml"               , emit: versions
```

CAT_CAT 모듈은 단일 입력을 받지만, 그 입력은 두 가지를 포함하는 튜플입니다:

- `meta`는 메타맵이라고 하는 메타데이터를 포함하는 구조입니다;
- `files_in`은 처리할 하나 이상의 입력 파일을 포함하며, `collectGreetings`의 `input_files`와 동등합니다.

완료되면, CAT_CAT은 출력을 두 부분으로 전달합니다:

- 메타맵과 연결된 출력 파일을 포함하는 또 다른 튜플로, `file_out` 태그로 내보내집니다;
- 사용된 소프트웨어 버전에 대한 정보를 캡처하는 `versions.yml` 파일로, `versions` 태그로 내보내집니다.

또한 기본적으로 출력 파일은 메타데이터의 일부인 식별자를 기반으로 이름이 지정됩니다(여기에는 코드가 표시되지 않음).

코드만 보면 추적하기 어려울 수 있으므로, 모든 것이 어떻게 함께 작동하는지 시각화하는 데 도움이 되는 다이어그램이 있습니다.

<figure class="excalidraw">
--8<-- "docs/en/docs/hello_nf-core/img/module_comparison.svg"
</figure>

두 모듈이 내용 측면에서 유사한 입력 요구 사항(입력 파일 세트와 일부 메타데이터)을 가지고 있지만, 해당 내용이 패키징되는 방식에 대한 기대는 매우 다르다는 것을 알 수 있습니다.
지금은 버전 파일을 무시하면, 주요 출력도 동등합니다(연결된 파일). 단, CAT_CAT은 출력 파일과 함께 메타맵도 내보냅니다.

곧 보시겠지만, 패키징 차이는 처리하기가 상당히 쉬울 것입니다.
하지만 메타맵 부분을 이해하려면 추가 컨텍스트를 소개해야 합니다.

### 2.2. 메타맵 이해하기

방금 CAT_CAT 모듈이 입력 튜플의 일부로 메타데이터 맵을 기대한다고 말씀드렸습니다.
이것이 무엇인지 자세히 살펴보는 데 몇 분을 할애하겠습니다.

**메타데이터 맵**은 종종 줄여서 **메타맵**이라고 하며, 데이터 단위에 대한 정보를 포함하는 Groovy 스타일 맵입니다.
Nextflow 파이프라인의 컨텍스트에서 데이터 단위는 원하는 모든 것이 될 수 있습니다: 개별 샘플, 샘플 배치 또는 전체 데이터 세트.

규칙에 따라, nf-core 메타맵은 `meta`라는 이름을 가지며 출력 이름 지정 및 데이터 단위 추적에 사용되는 필수 필드 `id`를 포함합니다.

예를 들어, 일반적인 메타데이터 맵은 다음과 같을 수 있습니다:

```groovy title="샘플 레벨 메타맵 예제"
[id: 'sample1', single_end: false, strandedness: 'forward']
```

또는 메타데이터가 배치 레벨에 첨부된 경우:

```groovy title="배치 레벨 메타맵 예제"
[id: 'batch1', date: '25.10.01']
```

이제 이것을 입력 파일이 메타맵과 함께 튜플로 패키징되기를 기대하고 출력 튜플의 일부로도 메타맵을 출력하는 `CAT_CAT` 프로세스의 컨텍스트에 넣어 보겠습니다.

```groovy title="modules/nf-core/cat/cat/main.nf (발췌)" linenums="1" hl_lines="2 5"
input:
tuple val(meta), path(files_in)

output:
tuple val(meta), path("${prefix}"), emit: file_out
```

결과적으로, 모든 데이터 단위는 관련 메타데이터가 첨부된 상태로 파이프라인을 통과합니다.
후속 프로세스도 해당 메타데이터에 쉽게 액세스할 수 있습니다.

`CAT_CAT`에서 출력된 파일이 메타데이터의 일부인 식별자를 기반으로 이름이 지정될 것이라고 말씀드렸던 것을 기억하시나요?
이것이 관련 코드입니다:

```groovy title="modules/nf-core/cat/cat/main.nf (발췌)" linenums="35"
prefix   = task.ext.prefix ?: "${meta.id}${getFileSuffix(file_list[0])}"
```

이것은 대략 다음과 같이 해석됩니다: 외부 작업 매개변수 시스템(`task.ext`)을 통해 `prefix`가 제공되면 그것을 사용하여 출력 파일 이름을 지정하고, 그렇지 않으면 메타맵의 `id` 필드에 해당하는 `${meta.id}`를 사용하여 생성합니다.

이 모듈로 들어오는 입력 채널이 다음과 같은 내용을 가지고 있다고 상상할 수 있습니다:

```groovy title="입력 채널 내용 예제"
ch_input = [[[id: 'batch1', date: '25.10.01'], ['file1A.txt', 'file1B.txt']],
            [[id: 'batch2', date: '25.10.26'], ['file2A.txt', 'file2B.txt']],
            [[id: 'batch3', date: '25.11.14'], ['file3A.txt', 'file3B.txt']]]
```

그러면 출력 채널 내용이 다음과 같이 나옵니다:

```groovy title="출력 채널 내용 예제"
ch_input = [[[id: 'batch1', date: '25.10.01'], 'batch1.txt'],
            [[id: 'batch2', date: '25.10.26'], 'batch2.txt'],
            [[id: 'batch3', date: '25.11.14'], 'batch3.txt']]
```

앞서 언급했듯이, `tuple val(meta), path(files_in)` 입력 설정은 모든 nf-core 모듈에서 사용되는 표준 패턴입니다.

이것이 얼마나 유용할 수 있는지 이해하기 시작하셨기를 바랍니다.
메타데이터를 기반으로 출력 이름을 지정할 수 있을 뿐만 아니라, 다른 매개변수 값을 적용하는 것과 같은 작업을 수행할 수 있으며, 특정 연산자와 결합하여 파이프라인을 통해 흐르는 데이터를 그룹화, 정렬 또는 필터링할 수도 있습니다.

!!! note "메타데이터에 대해 자세히 알아보기"

    샘플시트에서 메타데이터를 읽고 처리를 사용자 정의하는 데 사용하는 방법을 포함하여 Nextflow 워크플로우에서 메타데이터 작업에 대한 포괄적인 소개는 [워크플로우의 메타데이터](../side_quests/metadata) 사이드 퀘스트를 참조하세요.

### 2.3. 수행할 변경 사항 요약하기

검토한 내용을 바탕으로, `cat/cat` 모듈을 활용하기 위해 파이프라인에 수행해야 하는 주요 변경 사항은 다음과 같습니다:

- 배치 이름을 포함하는 메타맵 생성;
- 메타맵을 연결할 입력 파일 세트(`convertToUpper`에서 나오는)와 함께 튜플로 패키징;
- `collectGreetings()`에서 `CAT_CAT`으로 호출 전환;
- `CAT_CAT` 프로세스에서 생성된 튜플에서 출력 파일을 추출한 후 `cowpy`에 전달.

이것으로 충분할 것입니다! 이제 계획이 있으니 시작할 준비가 되었습니다.

### 핵심 정리

새 모듈의 입력 및 출력 인터페이스를 평가하여 요구 사항을 식별하는 방법을 알게 되었으며, 메타맵이 nf-core 파이프라인에서 파이프라인을 통해 흐르는 데이터와 메타데이터를 밀접하게 연결하는 데 어떻게 사용되는지 학습했습니다.

### 다음 단계

새 모듈을 워크플로우에 통합합니다.

---

## 3. CAT_CAT을 `hello.nf` 워크플로우에 통합하기

이제 메타맵에 대해 모든 것을 알게 되었으므로(또는 이 과정의 목적상 충분히), 위에서 설명한 변경 사항을 실제로 구현할 시간입니다.

명확성을 위해, 이것을 세분화하여 각 단계를 별도로 다루겠습니다.

!!! note

    아래 표시된 모든 변경 사항은 `core-hello/workflows/hello.nf` 워크플로우 파일의 `main` 블록에 있는 워크플로우 로직에 적용됩니다.

### 3.1. 메타데이터 맵 생성하기

먼저, nf-core 모듈이 메타맵에 최소한 `id` 필드를 요구한다는 점을 염두에 두고 `CAT_CAT`에 대한 메타데이터 맵을 생성해야 합니다.

다른 메타데이터가 필요하지 않으므로, 간단하게 유지하고 다음과 같은 것을 사용할 수 있습니다:

```groovy title="구문 예제"
def cat_meta = [id: 'test']
```

단, `id` 값을 하드코딩하고 싶지 않습니다. `params.batch` 매개변수의 값을 사용하고 싶습니다.
따라서 코드는 다음과 같이 됩니다:

```groovy title="구문 예제"
def cat_meta = [id: params.batch]
```

네, 기본 메타맵을 생성하는 것은 말 그대로 그렇게 간단합니다.

`convertToUpper` 호출 후에 이 줄을 추가하고, `collectGreetings` 호출을 제거하겠습니다:

=== "후"

    ```groovy title="core-hello/workflows/hello.nf" linenums="26" hl_lines="7-8"
        // 인사말을 내보냅니다
        sayHello(ch_samplesheet)

        // 인사말을 대문자로 변환합니다
        convertToUpper(sayHello.out)

        // 배치 이름을 ID로 사용하여 메타데이터 맵을 생성합니다
        def cat_meta = [ id: params.batch ]

        // cowpy로 인사말의 ASCII 아트를 생성합니다
        cowpy(collectGreetings.out.outfile, params.character)
    ```

=== "전"

    ```groovy title="core-hello/workflows/hello.nf" linenums="26" hl_lines="7-8"
        // 인사말을 내보냅니다
        sayHello(ch_samplesheet)

        // 인사말을 대문자로 변환합니다
        convertToUpper(sayHello.out)

        // 모든 인사말을 하나의 파일로 수집합니다
        collectGreetings(convertToUpper.out.collect(), params.batch)

        // cowpy로 인사말의 ASCII 아트를 생성합니다
        cowpy(collectGreetings.out.outfile, params.character)
    ```

이것은 `id`가 배치 이름으로 설정된 간단한 메타데이터 맵을 생성합니다(테스트 프로필을 사용할 때 `test`가 됩니다).

### 3.2. 메타데이터 튜플이 있는 채널 생성하기

다음으로, 파일 채널을 메타데이터와 파일을 포함하는 튜플 채널로 변환합니다:

=== "후"

    ```groovy title="core-hello/workflows/hello.nf" linenums="26" hl_lines="10-11"
        // 인사말을 내보냅니다
        sayHello(ch_samplesheet)

        // 인사말을 대문자로 변환합니다
        convertToUpper(sayHello.out)

        // 배치 이름을 ID로 사용하여 메타데이터 맵을 생성합니다
        def cat_meta = [ id: params.batch ]

        // 튜플 형식으로 메타데이터와 파일이 있는 채널을 생성합니다
        ch_for_cat = convertToUpper.out.collect().map { files -> tuple(cat_meta, files) }

        // cowpy로 인사말의 ASCII 아트를 생성합니다
        cowpy(collectGreetings.out.outfile, params.character)
    ```

=== "전"

    ```groovy title="core-hello/workflows/hello.nf" linenums="26"
        // 인사말을 내보냅니다
        sayHello(ch_samplesheet)

        // 인사말을 대문자로 변환합니다
        convertToUpper(sayHello.out)

        // 배치 이름을 ID로 사용하여 메타데이터 맵을 생성합니다
        def cat_meta = [ id: params.batch ]

        // cowpy로 인사말의 ASCII 아트를 생성합니다
        cowpy(collectGreetings.out.outfile, params.character)
    ```

추가한 줄은 두 가지를 달성합니다:

- `.collect()`는 `convertToUpper` 출력의 모든 파일을 단일 리스트로 수집합니다
- `.map { files -> tuple(cat_meta, files) }`는 `CAT_CAT`이 기대하는 형식인 `[metadata, files]`의 튜플을 생성합니다

이것이 `CAT_CAT`에 대한 입력 튜플을 설정하기 위해 해야 할 전부입니다.

### 3.3. CAT_CAT 모듈 호출하기

이제 새로 생성된 채널에서 `CAT_CAT`을 호출합니다:

=== "후"

    ```groovy title="core-hello/workflows/hello.nf" linenums="26" hl_lines="13-14"
        // 인사말을 내보냅니다
        sayHello(ch_samplesheet)

        // 인사말을 대문자로 변환합니다
        convertToUpper(sayHello.out)

        // 배치 이름을 ID로 사용하여 메타데이터 맵을 생성합니다
        def cat_meta = [ id: params.batch ]

        // 튜플 형식으로 메타데이터와 파일이 있는 채널을 생성합니다
        ch_for_cat = convertToUpper.out.collect().map { files -> tuple(cat_meta, files) }

        // nf-core cat/cat 모듈을 사용하여 파일을 연결합니다
        CAT_CAT(ch_for_cat)

        // cowpy로 인사말의 ASCII 아트를 생성합니다
        cowpy(collectGreetings.out.outfile, params.character)
    ```

=== "전"

    ```groovy title="core-hello/workflows/hello.nf" linenums="26"
        // 인사말을 내보냅니다
        sayHello(ch_samplesheet)

        // 인사말을 대문자로 변환합니다
        convertToUpper(sayHello.out)

        // 배치 이름을 ID로 사용하여 메타데이터 맵을 생성합니다
        def cat_meta = [ id: params.batch ]

        // 튜플 형식으로 메타데이터와 파일이 있는 채널을 생성합니다
        ch_for_cat = convertToUpper.out.collect().map { files -> tuple(cat_meta, files) }

        // cowpy로 인사말의 ASCII 아트를 생성합니다
        cowpy(collectGreetings.out.outfile, params.character)
    ```

이것으로 이 교체의 가장 까다로운 부분이 완료되었지만, 아직 완전히 끝나지 않았습니다. 연결된 출력을 `cowpy` 프로세스에 전달하는 방법을 업데이트해야 합니다.

### 3.4. `cowpy`를 위해 튜플에서 출력 파일 추출하기

이전에는 `collectGreetings` 프로세스가 `cowpy`에 직접 전달할 수 있는 파일만 생성했습니다.
하지만 `CAT_CAT` 프로세스는 출력 파일 외에 메타맵을 포함하는 튜플을 생성합니다.

`cowpy`가 아직 메타데이터 튜플을 받아들이지 않으므로(다음 과정 파트에서 수정하겠습니다), `CAT_CAT`에서 생성된 튜플에서 출력 파일을 추출한 후 `cowpy`에 전달해야 합니다:

=== "후"

    ```groovy title="core-hello/workflows/hello.nf" linenums="26" hl_lines="16-17 20"
        // 인사말을 내보냅니다
        sayHello(ch_samplesheet)

        // 인사말을 대문자로 변환합니다
        convertToUpper(sayHello.out)

        // 배치 이름을 ID로 사용하여 메타데이터 맵을 생성합니다
        def cat_meta = [ id: params.batch ]

        // 튜플 형식으로 메타데이터와 파일이 있는 채널을 생성합니다
        ch_for_cat = convertToUpper.out.collect().map { files -> tuple(cat_meta, files) }

        // 인사말을 연결합니다
        CAT_CAT(ch_for_cat)

        // cowpy가 아직 메타데이터를 사용하지 않으므로 튜플에서 파일을 추출합니다
        ch_for_cowpy = CAT_CAT.out.file_out.map{ meta, file -> file }

        // cowpy로 인사말의 ASCII 아트를 생성합니다
        cowpy(ch_for_cowpy, params.character)
    ```

=== "전"

    ```groovy title="core-hello/workflows/hello.nf" linenums="26" hl_lines="17"
        // 인사말을 내보냅니다
        sayHello(ch_samplesheet)

        // 인사말을 대문자로 변환합니다
        convertToUpper(sayHello.out)

        // 배치 이름을 ID로 사용하여 메타데이터 맵을 생성합니다
        def cat_meta = [ id: params.batch ]

        // 튜플 형식으로 메타데이터와 파일이 있는 채널을 생성합니다
        ch_for_cat = convertToUpper.out.collect().map { files -> tuple(cat_meta, files) }

        // 인사말을 연결합니다
        CAT_CAT(ch_for_cat)

        // cowpy로 인사말의 ASCII 아트를 생성합니다
        cowpy(collectGreetings.out.outfile, params.character)
    ```

`.map{ meta, file -> file }` 연산은 `CAT_CAT`에서 생성된 `[metadata, file]` 튜플에서 파일을 새 채널 `ch_for_cowpy`로 추출합니다.

그런 다음 마지막 줄에서 `collectGreetings.out.outfile` 대신 `ch_for_cowpy`를 `cowpy`에 전달하기만 하면 됩니다.

!!! note

    다음 과정 파트에서는 `cowpy`를 메타데이터 튜플과 직접 작동하도록 업데이트할 것이므로, 이 추출 단계는 더 이상 필요하지 않습니다.

### 3.5. 워크플로우 테스트하기

새로 통합된 `cat/cat` 모듈로 워크플로우가 작동하는지 테스트해 보겠습니다:

```bash
nextflow run . --outdir core-hello-results -profile test,docker --validate_params false
```

이것은 상당히 빠르게 실행되어야 합니다.

??? success "명령 출력"

    ```console
    N E X T F L O W ~ version 25.04.3

        Launching `./main.nf` [evil_pike] DSL2 - revision: b9e9b3b8de

        Input/output options
          input                     : /workspaces/training/hello-nf-core/core-hello/assets/greetings.csv
          outdir                    : core-hello-results

        Institutional config options
          config_profile_name       : Test profile
          config_profile_description: Minimal test dataset to check pipeline function

        Generic options
          validate_params           : false
          trace_report_suffix       : 2025-10-30_18-50-58

        Core Nextflow options
          runName                   : evil_pike
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
        [b3/f005fd] CORE_HELLO:HELLO:sayHello (3)       [100%] 3 of 3 ✔
        [08/f923d0] CORE_HELLO:HELLO:convertToUpper (3) [100%] 3 of 3 ✔
        [34/3729a9] CORE_HELLO:HELLO:CAT_CAT (test)     [100%] 1 of 1 ✔
        [24/df918a] CORE_HELLO:HELLO:cowpy              [100%] 1 of 1 ✔
        -[core/hello] Pipeline completed successfully-
    ```

`collectGreetings` 대신 `CAT_CAT`이 프로세스 실행 목록에 나타나는 것을 주목하세요.

그리고 그게 전부입니다! 이제 파이프라인의 해당 단계에 대해 커스텀 프로토타입 등급 코드 대신 강력한 커뮤니티 큐레이션 모듈을 사용하고 있습니다.

### 핵심 정리

이제 다음 방법을 알게 되었습니다:

- nf-core 모듈 찾기 및 설치하기
- nf-core 모듈의 요구 사항 평가하기
- nf-core 모듈과 함께 사용할 간단한 메타데이터 맵 생성하기
- nf-core 모듈을 워크플로우에 통합하기

### 다음 단계

로컬 모듈을 nf-core 규칙을 따르도록 조정하는 방법을 학습합니다.
또한 nf-core 도구를 사용하여 템플릿에서 새 nf-core 모듈을 생성하는 방법도 보여드리겠습니다.
