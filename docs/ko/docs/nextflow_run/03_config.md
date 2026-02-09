# 파트 3: 실행 구성

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } AI 기반 번역 - [자세히 알아보고 개선 사항 제안하기](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

이 섹션에서는 Nextflow 파이프라인의 구성을 관리하여 _워크플로우 코드 자체를 전혀 변경하지 않고도_ 동작을 사용자 정의하고, 다양한 환경에 적응시키며, 리소스 사용을 최적화하는 방법을 학습합니다.

이를 수행하는 방법은 여러 가지가 있으며, 이들을 조합하여 사용할 수 있고 [Configuration](https://nextflow.io/docs/latest/config.html) 문서에 설명된 우선순위에 따라 해석됩니다.

이 과정의 이 파트에서는 가장 간단하고 일반적인 구성 파일 메커니즘인 `nextflow.config` 파일을 소개합니다. 이 파일은 파트 2의 컨테이너 섹션에서 이미 접했습니다.

프로세스 지시문, 실행자, 프로필, 매개변수 파일과 같은 Nextflow 구성의 필수 구성 요소를 다룰 것입니다.
이러한 구성 옵션을 효과적으로 활용하는 방법을 학습함으로써 Nextflow 파이프라인의 유연성, 확장성 및 성능을 최대한 활용할 수 있습니다.

이러한 구성 요소를 연습하기 위해 이 교육 과정의 파트 2 마지막에 실행했던 워크플로우의 새 복사본을 `3-main.nf`로 이름을 변경하여 실행할 것입니다.

Hello 파이프라인에 익숙하지 않거나 복습이 필요한 경우 [이 정보 페이지](../info/hello_pipeline.md)를 참조하세요.

---

## 1. 워크플로우 입력 매개변수 관리

??? example "시나리오"

    파이프라인을 다운로드했고 동일한 입력 파일과 설정으로 반복적으로 실행하고 싶지만, 매번 모든 매개변수를 입력하고 싶지 않습니다.
    또는 명령줄 인수에 익숙하지 않은 동료를 위해 파이프라인을 설정하고 있을 수도 있습니다.

지금까지 작업해온 것의 확장인 구성 측면부터 시작하겠습니다: 입력 매개변수 관리입니다.

현재 워크플로우는 워크플로우 스크립트 자체의 `params` 블록에 선언된 여러 매개변수 값을 명령줄을 통해 받도록 설정되어 있습니다.
그 중 하나는 선언의 일부로 기본값이 설정되어 있습니다.

하지만 모든 매개변수에 대한 기본값을 설정하거나, 명령줄에서 매개변수를 지정하거나 원본 스크립트 파일을 수정하지 않고 기존 기본값을 재정의하고 싶을 수 있습니다.

이를 수행하는 방법은 여러 가지가 있습니다. 매우 일반적으로 사용되는 세 가지 기본 방법을 보여드리겠습니다.

### 1.1. `nextflow.config`에 값 설정

이것이 가장 간단한 접근 방식이지만, 메인 `nextflow.config` 파일은 매 실행마다 편집하고 싶은 것이 아니므로 유연성이 가장 낮을 수 있습니다.
하지만 워크플로우에서 매개변수를 _선언_하는 것(확실히 그곳에 속함)과 _기본값_을 제공하는 것(구성 파일에 더 적합함)의 관심사를 분리한다는 장점이 있습니다.

두 단계로 진행하겠습니다.

#### 1.1.1. 구성 파일에 `params` 블록 생성

`nextflow.config` 파일에서 다음 코드 변경을 수행하세요:

=== "후"

    ```groovy title="nextflow.config" linenums="1" hl_lines="3-10"
    docker.enabled = true

    /*
    * Pipeline parameters
    */
    params {
        input = 'data/greetings.csv'
        batch = 'batch'
        character = 'turkey'
    }
    ```

=== "전"

    ```groovy title="nextflow.config" linenums="1"
    docker.enabled = true
    ```

워크플로우의 `params` 블록을 구성 파일에 단순히 복사하지 않았다는 점에 유의하세요.
이미 기본값이 선언된 `batch` 매개변수의 경우 구문이 약간 다릅니다.
워크플로우 파일에서는 타입이 지정된 선언입니다.
구성에서는 값 할당입니다.

기술적으로 이것만으로도 워크플로우 파일에 여전히 지정된 기본값을 재정의하기에 충분합니다.
`batch`의 기본값을 수정하고 워크플로우를 실행하여 구성 파일에 설정된 값이 워크플로우 파일에 설정된 값을 재정의하는지 확인할 수 있습니다.

하지만 구성을 완전히 구성 파일로 이동하는 정신으로, 워크플로우 파일에서 해당 기본값을 완전히 제거하겠습니다.

#### 1.1.2. 워크플로우 파일에서 `batch`의 기본값 제거

`3-main.nf` 워크플로우 파일에서 다음 코드 변경을 수행하세요:

=== "후"

    ```groovy title="3-main.nf" linenums="9" hl_lines="6"
    /*
    * Pipeline parameters
    */
    params {
        input: Path
        batch: String
        character: String
    }
    ```

=== "전"

    ```groovy title="3-main.nf" linenums="9" hl_lines="6"
    /*
    * Pipeline parameters
    */
    params {
        input: Path
        batch: String = 'batch'
        character: String
    }
    ```

이제 워크플로우 파일 자체는 이러한 매개변수에 대한 기본값을 설정하지 않습니다.

#### 1.1.3. 파이프라인 실행

명령줄에서 매개변수를 지정하지 않고 올바르게 작동하는지 테스트해 봅시다.

```bash
nextflow run 3-main.nf
```

??? success "명령 출력"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `3-main.nf` [disturbed_einstein] DSL2 - revision: ede9037d02

    executor >  local (8)
    [f0/35723c] sayHello (2)       | 3 of 3 ✔
    [40/3efd1a] convertToUpper (3) | 3 of 3 ✔
    [17/e97d32] collectGreetings   | 1 of 1 ✔
    [98/c6b57b] cowpy              | 1 of 1 ✔
    ```

이것은 여전히 이전과 동일한 출력을 생성합니다.

최종 ASCII 아트 출력은 `results/3-main/` 디렉토리에 이전과 동일하게 `cowpy-COLLECTED-batch-output.txt`라는 이름으로 있습니다.

??? abstract "파일 내용"

    ```console title="results/3-main/cowpy-COLLECTED-batch-output.txt"
    _________
    / HOLà    \
    | HELLO   |
    \ BONJOUR /
    ---------
      \                                  ,+*^^*+___+++_
      \                           ,*^^^^              )
        \                       _+*                     ^**+_
        \                    +^       _ _++*+_+++_,         )
                  _+^^*+_    (     ,+*^ ^          \+_        )
                {       )  (    ,(    ,_+--+--,      ^)      ^\
                { (\@)    } f   ,(  ,+-^ __*_*_  ^^\_   ^\       )
              {:;-/    (_+*-+^^^^^+*+*<_ _++_)_    )    )      /
              ( /  (    (        ,___    ^*+_+* )   <    <      \
              U _/     )    *--<  ) ^\-----++__)   )    )       )
                (      )  _(^)^^))  )  )\^^^^^))^*+/    /       /
              (      /  (_))_^)) )  )  ))^^^^^))^^^)__/     +^^
            (     ,/    (^))^))  )  ) ))^^^^^^^))^^)       _)
              *+__+*       (_))^)  ) ) ))^^^^^^))^^^^^)____*^
              \             \_)^)_)) ))^^^^^^^^^^))^^^^)
              (_             ^\__^^^^^^^^^^^^))^^^^^^^)
                ^\___            ^\__^^^^^^))^^^^^^^^)\\
                      ^^^^^\uuu/^^\uuu/^^^^\^\^\^\^\^\^\^\
                        ___) >____) >___   ^\_\_\_\_\_\_\)
                        ^^^//\\_^^//\\_^       ^(\_\_\_\)
                          ^^^ ^^ ^^^ ^
    ```

기능적으로 이 이동은 아무것도 변경하지 않았지만, 개념적으로는 구성 파일에 기본값을 설정하는 것이 조금 더 깔끔합니다.

### 1.2. 실행별 구성 파일 사용

??? example "시나리오"

    메인 구성 파일을 수정하지 않고 다양한 설정을 실험하고 싶습니다.

실험을 위한 작업 디렉토리로 사용할 하위 디렉토리에 새 `nextflow.config` 파일을 생성하여 이를 수행할 수 있습니다.

#### 1.2.1. 빈 구성으로 작업 디렉토리 생성

먼저 새 디렉토리를 생성하고 그 안으로 이동합니다:

```bash
mkdir -p tux-run
cd tux-run
```

그런 다음 해당 디렉토리에 빈 구성 파일을 생성합니다:

```bash
touch nextflow.config
```

이렇게 하면 빈 파일이 생성됩니다.

#### 1.2.2. 실험 구성 설정

이제 새 파일을 열고 사용자 정의하려는 매개변수를 추가합니다:

```groovy title="tux-run/nextflow.config" linenums="1"
params {
    input = '../data/greetings.csv'
    batch = 'experiment'
    character = 'tux'
}
```

입력 파일의 경로는 디렉토리 구조를 반영해야 한다는 점에 유의하세요.

#### 1.2.3. 파이프라인 실행

이제 새 작업 디렉토리 내에서 파이프라인을 실행할 수 있습니다.
경로를 적절히 조정해야 합니다!

```bash
nextflow run ../3-main.nf
```

??? success "명령 출력"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `../3-main.nf` [trusting_escher] DSL2 - revision: 356df0818d

    executor >  local (8)
    [59/b66913] sayHello (2)       [100%] 3 of 3 ✔
    [ad/f06364] convertToUpper (3) [100%] 3 of 3 ✔
    [10/714895] collectGreetings   [100%] 1 of 1 ✔
    [88/3ece98] cowpy              [100%] 1 of 1 ✔
    ```

이렇게 하면 `tux-run/work/` 및 `tux-run/results/`를 포함하여 `tux-run/` 아래에 새로운 디렉토리 세트가 생성됩니다.

이 실행에서 Nextflow는 현재 디렉토리의 `nextflow.config`와 파이프라인 루트 디렉토리의 `nextflow.config`를 결합하여 기본 캐릭터(turkey)를 tux 캐릭터로 재정의합니다.

최종 출력 파일에는 인사말을 하는 tux 캐릭터가 포함되어야 합니다.

??? abstract "파일 내용"

    ```console title="tux-run/results/3-main/cowpy-COLLECTED-experiment-output.txt"
    _________
    / HELLO   \
    | BONJOUR |
    \ HOLà    /
    ---------
      \
        \
            .--.
          |o_o |
          |:_/ |
          //   \ \
        (|     | )
        /'\_   _/`\
        \___)=(___/

    ```

이제 '일반' 구성을 수정하지 않고 실험할 수 있는 공간이 생겼습니다.

!!! warning

    다음 섹션으로 이동하기 전에 이전 디렉토리로 돌아가야 합니다!

    ```bash
    cd ..
    ```

이제 매개변수 값을 설정하는 또 다른 유용한 방법을 살펴보겠습니다.

### 1.3. 매개변수 파일 사용

??? example "시나리오"

    공동 작업자와 정확한 실행 매개변수를 공유하거나 출판물을 위해 기록해야 합니다.

하위 디렉토리 접근 방식은 실험에 적합하지만 약간의 설정이 필요하고 경로를 적절히 조정해야 합니다.
특정 값 세트로 파이프라인을 실행하거나 다른 사람이 최소한의 노력으로 실행할 수 있도록 하려는 경우 더 간단한 접근 방식이 있습니다.

Nextflow를 사용하면 YAML 또는 JSON 형식의 [매개변수 파일](https://nextflow.io/docs/latest/config.html#parameter-file)을 통해 매개변수를 지정할 수 있으므로 대체 기본값 세트와 실행별 매개변수 값을 관리하고 배포하는 것이 매우 편리합니다.

#### 1.3.1. 예제 매개변수 파일 검토

이를 시연하기 위해 현재 디렉토리에 `test-params.yaml`이라는 예제 매개변수 파일을 제공합니다:

```yaml title="test-params.yaml" linenums="1"
input: "data/greetings.csv"
batch: "yaml"
character: "stegosaurus"
```

이 매개변수 파일에는 지정하려는 각 입력에 대한 키-값 쌍이 포함되어 있습니다.
구성 파일과 비교할 때 등호(`=`) 대신 콜론(`:`)을 사용한다는 점에 유의하세요.
구성 파일은 Groovy로 작성되는 반면, 매개변수 파일은 YAML로 작성됩니다.

!!! info

    예제로 매개변수 파일의 JSON 버전도 제공하지만 여기서는 실행하지 않을 것입니다.
    직접 시도해 보세요.

#### 1.3.2. 파이프라인 실행

이 매개변수 파일로 워크플로우를 실행하려면 기본 명령에 `-params-file <filename>`을 추가하기만 하면 됩니다.

```bash
nextflow run 3-main.nf -params-file test-params.yaml
```

??? success "명령 출력"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `3-main.nf` [disturbed_sammet] DSL2 - revision: ede9037d02

    executor >  local (8)
    [2b/9a7d1e] sayHello (2)       | 3 of 3 ✔
    [5c/8f3b2a] convertToUpper (3) | 3 of 3 ✔
    [a3/29d8fb] collectGreetings   | 1 of 1 ✔
    [b7/83ef12] cowpy              | 1 of 1 ✔
    ```

최종 출력 파일에는 인사말을 하는 stegosaurus 캐릭터가 포함되어야 합니다.

??? abstract "파일 내용"

    ```console title="results/3-main/cowpy-COLLECTED-yaml-output.txt"
    _________
    / HELLO   \
    | HOLà    |
    \ BONJOUR /
    ---------
    \                             .       .
    \                           / `.   .' "
      \                  .---.  <    > <    >  .---.
      \                 |    \  \ - ~ ~ - /  /    |
            _____          ..-~             ~-..-~
            |     |   \~~~\.'                    `./~~~/
          ---------   \__/                        \__/
          .'  O    \     /               /       \  "
        (_____,    `._.'               |         }  \/~~~/
          `----.          /       }     |        /    \__/
                `-.      |       /      |       /      `. ,~~|
                    ~-.__|      /_ - ~ ^|      /- _      `..-'
                        |     /        |     /     ~-.     `-. _  _  _
                        |_____|        |_____|         ~ - . _ _ _ _ _>
    ```

지정할 매개변수가 몇 개뿐일 때는 매개변수 파일을 사용하는 것이 과도해 보일 수 있지만, 일부 파이프라인은 수십 개의 매개변수를 예상합니다.
이러한 경우 매개변수 파일을 사용하면 대규모 명령줄을 입력하거나 워크플로우 스크립트를 수정하지 않고도 런타임에 매개변수 값을 제공할 수 있습니다.

또한 공동 작업자에게 매개변수 세트를 배포하거나 출판물의 보조 정보로 배포하는 것이 더 쉬워집니다.
이를 통해 다른 사람들이 작업을 더 재현 가능하게 만듭니다.

### 핵심 정리

워크플로우 입력 관리를 위한 주요 구성 옵션을 활용하는 방법을 알게 되었습니다.

### 다음 단계

워크플로우 출력이 게시되는 위치와 방법을 관리하는 방법을 학습합니다.

---

## 2. 워크플로우 출력 관리

??? example "시나리오"

    파이프라인이 하드코딩된 디렉토리에 출력을 게시하지만, 매번 워크플로우 코드를 편집하지 않고 프로젝트 또는 실험 이름별로 결과를 구성하고 싶습니다.

우리가 상속받은 워크플로우는 워크플로우 수준 출력 선언에 경로를 사용하는데, 이는 그다지 유연하지 않고 많은 반복이 포함됩니다.

이를 더 유연하게 구성할 수 있는 몇 가지 일반적인 방법을 살펴보겠습니다.

### 2.1. `outputDir` 디렉토리 이름 사용자 정의

지금까지 실행한 각 버전의 워크플로우는 출력 정의에 하드코딩된 다른 하위 디렉토리에 출력을 게시했습니다.

파트 1에서 `-output-dir` CLI 플래그를 사용하여 해당 하위 디렉토리의 위치를 변경했지만, 그것은 여전히 정적 문자열일 뿐입니다.
대신 구성 파일에서 이를 구성하여 더 복잡한 동적 경로를 정의할 수 있습니다.
이를 위해 완전히 새로운 매개변수를 만들 수도 있지만, 바로 거기에 있는 `batch` 매개변수를 사용하겠습니다.

#### 2.1.1. 구성 파일에 `outputDir` 값 설정

Nextflow가 출력 게시에 사용하는 경로는 `outputDir` 옵션으로 제어됩니다.
모든 출력의 경로를 변경하려면 `nextflow.config` 구성 파일에서 이 옵션의 값을 설정할 수 있습니다.

`nextflow.config` 파일에 다음 코드를 추가하세요:

=== "후"

    ```groovy title="nextflow.config" linenums="9" hl_lines="10-13"
    /*
    * Pipeline parameters
    */
    params {
        input = 'data/greetings.csv'
        batch = 'batch'
        character = 'turkey'
    }

    /*
    * Output settings
    */
    outputDir = "results_config/${params.batch}"
    ```

=== "전"

    ```groovy title="nextflow.config" linenums="9"
    /*
    * Pipeline parameters
    */
    params {
        input = 'data/greetings.csv'
        batch = 'batch'
        character = 'turkey'
    }
    ```

이렇게 하면 내장 기본 경로인 `results/`가 `results_config/`와 하위 디렉토리로서의 `batch` 매개변수 값으로 대체됩니다.

명령의 `-output-dir` 매개변수(`-o`로 줄임)를 사용하여 명령줄에서도 이 옵션을 설정할 수 있지만, 그러면 `batch` 매개변수 값을 사용할 수 없습니다.
CLI 플래그를 사용하면 구성에 설정된 경우 `outputDir`을 덮어씁니다.

#### 2.1.2. 하드코딩된 경로의 반복 부분 제거

출력 옵션에 여전히 하드코딩된 하위 디렉토리가 있으므로 이제 제거하겠습니다.

워크플로우 파일에서 다음 코드 변경을 수행하세요:

=== "후"

    ```groovy title="3-main.nf" linenums="42" hl_lines="3 7 11 15 19"
    output {
        first_output {
            path 'intermediates'
            mode 'copy'
        }
        uppercased {
            path 'intermediates'
            mode 'copy'
        }
        collected {
            path 'intermediates'
            mode 'copy'
        }
        batch_report {
            path ''
            mode 'copy'
        }
        cowpy_art {
            path ''
            mode 'copy'
        }
    }
    ```

=== "전"

    ```groovy title="3-main.nf" linenums="42" hl_lines="3 7 11 15 19"
    output {
        first_output {
            path '3-main/intermediates'
            mode 'copy'
        }
        uppercased {
            path '3-main/intermediates'
            mode 'copy'
        }
        collected {
            path '3-main/intermediates'
            mode 'copy'
        }
        batch_report {
            path '3-main'
            mode 'copy'
        }
        cowpy_art {
            path '3-main'
            mode 'copy'
        }
    }
    ```

`outputDir` 기본값을 수정하는 대신 각 경로에 `${params.batch}`를 추가할 수도 있었지만, 이것이 더 간결합니다.

#### 2.1.3. 파이프라인 실행

명령줄에서 배치 이름을 `outdir`로 설정하여 올바르게 작동하는지 테스트해 봅시다.

```bash
nextflow run 3-main.nf --batch outdir
```

??? success "명령 출력"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `3-main.nf` [amazing_church] DSL2 - revision: 6e18cd130e

    executor >  local (8)
    [9c/6a03ea] sayHello (2)       [100%] 3 of 3 ✔
    [11/9e58a6] convertToUpper (3) [100%] 3 of 3 ✔
    [c8/1977e5] collectGreetings   [100%] 1 of 1 ✔
    [38/f01eda] cowpy              [100%] 1 of 1 ✔
    ```

이것은 여전히 이전과 동일한 출력을 생성하지만, 이번에는 `results_config/outdir/` 아래에서 출력을 찾습니다.

??? abstract "디렉토리 내용"

    ```console
    results_config/outdir
    ├── cowpy-COLLECTED-outdir-output.txt
    ├── intermediates
    │   ├── Bonjour-output.txt
    │   ├── COLLECTED-outdir-output.txt
    │   ├── Hello-output.txt
    │   ├── Holà-output.txt
    │   ├── UPPER-Bonjour-output.txt
    │   ├── UPPER-Hello-output.txt
    │   └── UPPER-Holà-output.txt
    └── outdir-report.txt
    ```

이 접근 방식을 사용자 정의 경로 정의와 결합하여 원하는 디렉토리 계층 구조를 구성할 수 있습니다.

### 2.2. 프로세스별로 출력 구성

출력을 추가로 구성하는 인기 있는 방법 중 하나는 프로세스별로 수행하는 것입니다. _즉,_ 파이프라인에서 실행되는 각 프로세스에 대한 하위 디렉토리를 생성합니다.

#### 2.2.1. 출력 경로를 프로세스 이름 참조로 대체

출력 경로 선언에서 프로세스 이름을 `<process>.name`으로 참조하기만 하면 됩니다.

워크플로우 파일에서 다음 변경을 수행하세요:

=== "후"

    ```groovy title="3-main.nf" linenums="42" hl_lines="3 7 11 15 19"
    output {
        first_output {
            path { sayHello.name }
            mode 'copy'
        }
        uppercased {
            path { convertToUpper.name }
            mode 'copy'
        }
        collected {
            path { collectGreetings.name }
            mode 'copy'
        }
        batch_report {
            path { collectGreetings.name }
            mode 'copy'
        }
        cowpy_art {
            path { cowpy.name }
            mode 'copy'
        }
    }
    ```

=== "전"

    ```groovy title="3-main.nf" linenums="42" hl_lines="3 7 11 15 19"
    output {
        first_output {
            path 'intermediates'
            mode 'copy'
        }
        uppercased {
            path 'intermediates'
            mode 'copy'
        }
        collected {
            path 'intermediates'
            mode 'copy'
        }
        batch_report {
            path ''
            mode 'copy'
        }
        cowpy_art {
            path ''
            mode 'copy'
        }
    }
    ```

이렇게 하면 출력 경로 구성에서 나머지 하드코딩된 요소가 제거됩니다.

#### 2.2.2. 파이프라인 실행

명령줄에서 배치 이름을 `pnames`로 설정하여 올바르게 작동하는지 테스트해 봅시다.

```bash
nextflow run 3-main.nf --batch pnames
```

??? success "명령 출력"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `3-main.nf` [jovial_mcclintock] DSL2 - revision: ede9037d02

    executor >  local (8)
    [4a/c2e6b8] sayHello (2)       | 3 of 3 ✔
    [6f/d4a172] convertToUpper (3) | 3 of 3 ✔
    [e8/4f19d7] collectGreetings   | 1 of 1 ✔
    [f2/a85c36] cowpy              | 1 of 1 ✔
    ```

이것은 여전히 이전과 동일한 출력을 생성하지만, 이번에는 `results_config/pnames/` 아래에서 출력을 찾으며 프로세스별로 그룹화되어 있습니다.

??? abstract "디렉토리 내용"

    ```console
    results_config/pnames/
    ├── collectGreetings
    │   ├── COLLECTED-pnames-output.txt
    │   └── pnames-report.txt
    ├── convertToUpper
    │   ├── UPPER-Bonjour-output.txt
    │   ├── UPPER-Hello-output.txt
    │   └── UPPER-Holà-output.txt
    ├── cowpy
    │   └── cowpy-COLLECTED-pnames-output.txt
    └── sayHello
        ├── Bonjour-output.txt
        ├── Hello-output.txt
        └── Holà-output.txt
    ```

!!! note

    여기서는 `intermediates`와 최상위 수준의 최종 출력 간의 구분을 지웠습니다.
    이러한 접근 방식을 혼합하고 일치시킬 수 있으며, 예를 들어 첫 번째 출력의 경로를 `#!groovy "${params.batch}/intermediates/${sayHello.name}"`으로 설정하여 여러 변수를 포함할 수도 있습니다.

### 2.3. 워크플로우 수준에서 게시 모드 설정

마지막으로, 반복적인 코드의 양을 줄이는 정신으로 출력별 `mode` 선언을 구성의 단일 줄로 대체할 수 있습니다.

#### 2.3.1. 구성 파일에 `workflow.output.mode` 추가

`nextflow.config` 파일에 다음 코드를 추가하세요:

=== "후"

    ```groovy title="nextflow.config" linenums="2" hl_lines="5"
    /*
    * Output settings
    */
    outputDir = "results_config/${params.batch}"
    workflow.output.mode = 'copy'
    ```

=== "전"

    ```groovy title="nextflow.config" linenums="12"
    /*
    * Output settings
    */
    outputDir = "results_config/${params.batch}"
    ```

`outputDir` 옵션과 마찬가지로 구성 파일에서 `workflow.output.mode`에 값을 지정하면 워크플로우 파일에 설정된 것을 재정의하기에 충분하지만, 어쨌든 불필요한 코드를 제거하겠습니다.

#### 2.3.2. 워크플로우 파일에서 출력 모드 제거

워크플로우 파일에서 다음 변경을 수행하세요:

=== "후"

    ```groovy title="3-main.nf" linenums="42"
    output {
        first_output {
            path { sayHello.name }
        }
        uppercased {
            path { convertToUpper.name }
        }
        collected {
            path { collectGreetings.name }
        }
        batch_report {
            path { collectGreetings.name }
        }
        cowpy_art {
            path { cowpy.name }
        }
    }
    ```

=== "전"

    ```groovy title="3-main.nf" linenums="42" hl_lines="3 7 11 15 19"
    output {
        first_output {
            path { sayHello.name }
            mode 'copy'
        }
        uppercased {
            path { convertToUpper.name }
            mode 'copy'
        }
        collected {
            path { collectGreetings.name }
            mode 'copy'
        }
        batch_report {
            path { collectGreetings.name }
            mode 'copy'
        }
        cowpy_art {
            path { cowpy.name }
            mode 'copy'
        }
    }
    ```

더 간결하지 않나요?

#### 2.3.3. 파이프라인 실행

명령줄에서 배치 이름을 `outmode`로 설정하여 올바르게 작동하는지 테스트해 봅시다.

```bash
nextflow run 3-main.nf --batch outmode
```

??? success "명령 출력"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `3-main.nf` [rowdy_sagan] DSL2 - revision: ede9037d02

    executor >  local (8)
    [5b/d91e3c] sayHello (2)       | 3 of 3 ✔
    [8a/f6c241] convertToUpper (3) | 3 of 3 ✔
    [89/cd3a48] collectGreetings   | 1 of 1 ✔
    [9e/71fb52] cowpy              | 1 of 1 ✔
    ```

이것은 여전히 이전과 동일한 출력을 생성하지만, 이번에는 `results_config/outmode/` 아래에서 출력을 찾습니다.
모두 심볼릭 링크가 아닌 적절한 복사본입니다.

??? abstract "디렉토리 내용"

    ```console
    results_config/outmode/
    ├── collectGreetings
    │   ├── COLLECTED-outmode-output.txt
    │   └── outmode-report.txt
    ├── convertToUpper
    │   ├── UPPER-Bonjour-output.txt
    │   ├── UPPER-Hello-output.txt
    │   └── UPPER-Holà-output.txt
    ├── cowpy
    │   └── cowpy-COLLECTED-outmode-output.txt
    └── sayHello
        ├── Bonjour-output.txt
        ├── Hello-output.txt
        └── Holà-output.txt
    ```

출력별 모드 설정 방식을 여전히 사용하려는 주된 이유는 동일한 워크플로우 내에서 혼합하고 일치시키려는 경우입니다. _즉,_ 일부 출력은 복사하고 일부는 심볼릭 링크로 만들려는 경우입니다.

이러한 방식으로 사용자 정의할 수 있는 다른 많은 옵션이 있지만, 이것이 옵션의 범위와 선호도에 맞게 효과적으로 활용하는 방법에 대한 감각을 제공하기를 바랍니다.

### 핵심 정리

출력이 게시되는 디렉토리의 이름 지정 및 구조와 워크플로우 출력 게시 모드를 제어하는 방법을 알게 되었습니다.

### 다음 단계

소프트웨어 패키징 기술부터 시작하여 컴퓨팅 환경에 워크플로우 구성을 적응시키는 방법을 학습합니다.

---

## 3. 소프트웨어 패키징 기술 선택

지금까지 입력이 들어가는 방법과 출력이 나오는 위치를 제어하는 구성 요소를 살펴보았습니다. 이제 컴퓨팅 환경에 워크플로우 구성을 적응시키는 데 더 구체적으로 집중할 시간입니다.

그 경로의 첫 번째 단계는 각 단계에서 실행될 소프트웨어 패키지가 어디에서 올 것인지 지정하는 것입니다.
로컬 컴퓨팅 환경에 이미 설치되어 있나요?
컨테이너 시스템을 통해 이미지를 검색하고 실행해야 하나요?
아니면 Conda 패키지를 검색하고 로컬 Conda 환경을 구축해야 하나요?

이 교육 과정의 첫 번째 파트(파트 1-4)에서는 워크플로우에서 로컬에 설치된 소프트웨어만 사용했습니다.
그런 다음 파트 5에서 Docker 컨테이너와 `nextflow.config` 파일을 소개했으며, 이를 사용하여 Docker 컨테이너 사용을 활성화했습니다.

이제 `nextflow.config` 파일을 통해 대체 소프트웨어 패키징 옵션을 구성하는 방법을 살펴보겠습니다.

### 3.1. 구성 파일에서 Docker 비활성화 및 Conda 활성화

??? example "시나리오"

    보안상의 이유로 Docker가 허용되지 않는 HPC 클러스터로 파이프라인을 이동하고 있습니다.
    클러스터는 Singularity와 Conda를 지원하므로 그에 따라 구성을 전환해야 합니다.

앞서 언급했듯이 Nextflow는 HPC에서 더 널리 사용되는 Singularity를 포함한 여러 컨테이너 기술과 Conda와 같은 소프트웨어 패키지 관리자를 지원합니다.

Docker 대신 Conda를 사용하도록 구성 파일을 변경할 수 있습니다.
이를 위해 `docker.enabled`의 값을 `false`로 전환하고 Conda 사용을 활성화하는 지시문을 추가하겠습니다:

=== "후"

    ```groovy title="nextflow.config" linenums="1" hl_lines="1-2"
    docker.enabled = false
    conda.enabled = true
    ```

=== "전"

    ```groovy title="nextflow.config" linenums="1" hl_lines="1"
    docker.enabled = true
    ```

이렇게 하면 Nextflow가 Conda 패키지가 지정된 프로세스에 대한 Conda 환경을 생성하고 활용할 수 있습니다.
즉, 이제 `cowpy` 프로세스에 그 중 하나를 추가해야 합니다!

### 3.2. 프로세스 정의에 Conda 패키지 지정

`cowpy` 도구가 포함된 Conda 패키지의 URI를 이미 검색했습니다: `conda-forge::cowpy==1.1.5`

이제 `conda` 지시문을 사용하여 `cowpy` 프로세스 정의에 URI를 추가합니다:

=== "후"

    ```groovy title="modules/cowpy.nf" linenums="4" hl_lines="4"
    process cowpy {

        container 'community.wave.seqera.io/library/cowpy:1.1.5--3db457ae1977a273'
        conda 'conda-forge::cowpy==1.1.5'

        input:
    ```

=== "전"

    ```groovy title="modules/cowpy.nf" linenums="4"
    process cowpy {

        container 'community.wave.seqera.io/library/cowpy:1.1.5--3db457ae1977a273'

        input:
    ```

명확히 하자면, `docker` 지시문을 _대체_하는 것이 아니라 대체 옵션을 _추가_하는 것입니다.

!!! tip

    주어진 conda 패키지의 URI를 얻는 방법은 몇 가지가 있습니다.
    컨테이너를 만들 계획이 없더라도 복사하여 붙여넣을 수 있는 URI를 제공하는 [Seqera Containers](https://seqera.io/containers/) 검색 쿼리를 사용하는 것이 좋습니다.

### 3.3. 워크플로우를 실행하여 Conda를 사용할 수 있는지 확인

시도해 봅시다.

```bash
nextflow run 3-main.nf --batch conda
```

??? success "명령 출력"

    ```console title="출력"
    N E X T F L O W   ~  version 25.10.2

    Launching `3-main.nf` [trusting_lovelace] DSL2 - revision: 028a841db1

    executor >  local (8)
    [ee/4ca1f2] sayHello (3)       | 3 of 3 ✔
    [20/2596a7] convertToUpper (1) | 3 of 3 ✔
    [b3/e15de5] collectGreetings   | 1 of 1 ✔
    [c5/af5f88] cowpy              | 1 of 1 ✔
    ```

이것은 문제없이 작동하고 `results_config/conda` 아래에 이전과 동일한 출력을 생성해야 합니다.

백그라운드에서 Nextflow는 Conda 패키지를 검색하고 환경을 생성했으며, 이는 일반적으로 약간의 작업이 필요합니다. 따라서 우리가 직접 할 필요가 없다는 것이 좋습니다!

!!! info

    `cowpy` 패키지가 매우 작기 때문에 빠르게 실행되지만, 큰 패키지로 작업하는 경우 처음에는 평소보다 조금 더 오래 걸릴 수 있으며 콘솔 출력이 완료되기 전에 1분 정도 '멈춰' 있는 것을 볼 수 있습니다.
    이것은 정상이며 Nextflow가 새 패키지를 처음 사용할 때 수행하는 추가 작업 때문입니다.

우리 관점에서는 백엔드의 메커니즘이 약간 다르더라도 Docker로 실행하는 것과 정확히 동일하게 작동하는 것처럼 보입니다.

즉, 필요한 경우 Conda 환경으로 실행할 준비가 되었습니다.

??? info "Docker와 Conda 혼합 및 일치"

    이러한 지시문은 프로세스별로 할당되므로 '혼합 및 일치'가 가능합니다. _즉,_ 사용 중인 컴퓨팅 인프라가 둘 다 지원하는 경우 워크플로우의 일부 프로세스는 Docker로 실행하고 다른 프로세스는 Conda로 실행하도록 구성할 수 있습니다.
    이 경우 구성 파일에서 Docker와 Conda를 모두 활성화합니다.
    주어진 프로세스에 대해 둘 다 사용 가능한 경우 Nextflow는 컨테이너를 우선시합니다.

    그리고 앞서 언급했듯이 Nextflow는 여러 다른 소프트웨어 패키징 및 컨테이너 기술을 지원하므로 이 두 가지로만 제한되지 않습니다.

### 핵심 정리

각 프로세스가 사용해야 하는 소프트웨어 패키지를 구성하는 방법과 기술 간에 전환하는 방법을 알게 되었습니다.

### 다음 단계

Nextflow가 실제로 작업을 수행하는 데 사용하는 실행 플랫폼을 변경하는 방법을 학습합니다.

---

## 4. 실행 플랫폼 선택

??? example "시나리오"

    노트북에서 파이프라인을 개발하고 테스트해 왔지만 이제 수천 개의 샘플에서 실행해야 합니다.
    귀하의 기관에는 대신 사용하고 싶은 Slurm 스케줄러가 있는 HPC 클러스터가 있습니다.

지금까지 로컬 실행자로 파이프라인을 실행해 왔습니다.
이것은 Nextflow가 실행 중인 머신에서 각 작업을 실행합니다.
Nextflow가 시작되면 사용 가능한 CPU와 메모리를 확인합니다.
실행 준비가 된 작업의 리소스가 사용 가능한 리소스를 초과하면 Nextflow는 이전 작업 중 하나 이상이 완료되어 필요한 리소스가 확보될 때까지 마지막 작업을 실행에서 보류합니다.

로컬 실행자는 편리하고 효율적이지만 단일 머신으로 제한됩니다. 매우 큰 워크로드의 경우 사용 가능한 것보다 더 많은 리소스가 필요한 단일 작업이 있거나 단일 머신이 실행하기를 기다리는 데 너무 오래 걸리는 작업이 너무 많아 로컬 머신이 병목 현상이 될 수 있습니다.

Nextflow는 HPC 스케줄러(Slurm, LSF, SGE, PBS, Moab, OAR, Bridge, HTCondor 등)와 클라우드 실행 백엔드(AWS Batch, Google Cloud Batch, Azure Batch, Kubernetes 등)를 포함한 [다양한 실행 백엔드](https://nextflow.io/docs/latest/executor.html)를 지원합니다.

### 4.1. 다른 백엔드 대상 지정

실행자의 선택은 `executor`라는 프로세스 지시문으로 설정됩니다.
기본적으로 `local`로 설정되어 있으므로 다음 구성이 암시됩니다:

```groovy title="내장 구성"
process {
    executor = 'local'
}
```

다른 백엔드를 대상으로 실행자를 설정하려면 리소스 할당에 대해 위에서 설명한 것과 유사한 구문을 사용하여 원하는 실행자를 지정하기만 하면 됩니다(모든 옵션은 [Executors](https://nextflow.io/docs/latest/executor.html) 참조).

```groovy title="nextflow.config"
process {
    executor = 'slurm'
}
```

!!! warning

    교육 환경이 HPC에 연결되도록 설정되어 있지 않기 때문에 실제로 테스트할 수 없습니다.

### 4.2. 실행 매개변수에 대한 백엔드별 구문 처리

대부분의 고성능 컴퓨팅 플랫폼은 CPU 수 및 메모리와 같은 리소스 할당 요청 및 제한과 사용할 작업 큐 이름과 같은 특정 매개변수를 지정할 수 있도록(때로는 요구) 허용합니다.

불행히도 이러한 각 시스템은 작업을 정의하고 관련 스케줄러에 제출하는 방법을 정의하기 위해 서로 다른 기술, 구문 및 구성을 사용합니다.

??? abstract "예제"

    예를 들어, "my-science-work" 큐에서 실행될 8개의 CPU와 4GB의 RAM이 필요한 동일한 작업은 백엔드에 따라 다음과 같은 다른 방식으로 표현되어야 합니다.

    ```bash title="SLURM 구성 / sbatch를 사용하여 제출"
    #SBATCH -o /path/to/my/task/directory/my-task-1.log
    #SBATCH --no-requeue
    #SBATCH -c 8
    #SBATCH --mem 4096M
    #SBATCH -p my-science-work
    ```

    ```bash title="PBS 구성 / qsub를 사용하여 제출"
    #PBS -o /path/to/my/task/directory/my-task-1.log
    #PBS -j oe
    #PBS -q my-science-work
    #PBS -l nodes=1:ppn=5
    #PBS -l mem=4gb
    ```

    ```bash title="SGE 구성 / qsub를 사용하여 제출"
    #$ -o /path/to/my/task/directory/my-task-1.log
    #$ -j y
    #$ -terse
    #$ -notify
    #$ -q my-science-work
    #$ -l slots=5
    #$ -l h_rss=4096M,mem_free=4096M
    ```

다행히 Nextflow는 이 모든 것을 단순화합니다.
`cpus`, `memory`, `queue`와 같은 관련 속성을 한 번만 지정할 수 있도록 표준화된 구문을 제공합니다(사용 가능한 모든 옵션은 [Process directives](https://nextflow.io/docs/latest/reference/process.html#process-directives) 참조).
그런 다음 런타임에 Nextflow는 이러한 설정을 사용하여 실행자 설정에 따라 적절한 백엔드별 스크립트를 생성합니다.

다음 섹션에서 해당 표준화된 구문을 다룰 것입니다.

### 핵심 정리

이제 다양한 종류의 컴퓨팅 인프라를 사용하기 위해 실행자를 변경하는 방법을 알게 되었습니다.

### 다음 단계

Nextflow에서 리소스 할당 및 제한을 평가하고 표현하는 방법을 학습합니다.

---

## 5. 컴퓨팅 리소스 할당 제어

??? example "시나리오"

    파이프라인이 메모리 제한을 초과하여 작업이 종료되어 클러스터에서 계속 실패합니다.
    또는 사용하지 않는 리소스에 대해 요금이 청구되고 있으며 비용을 최적화하고 싶을 수 있습니다.

대부분의 고성능 컴퓨팅 플랫폼은 CPU 수 및 메모리와 같은 특정 리소스 할당 매개변수를 지정할 수 있도록(때로는 요구) 허용합니다.

기본적으로 Nextflow는 각 프로세스에 대해 단일 CPU와 2GB의 메모리를 사용합니다.
해당 프로세스 지시문은 `cpus` 및 `memory`라고 하므로 다음 구성이 암시됩니다:

```groovy title="내장 구성" linenums="1"
process {
    cpus = 1
    memory = 2.GB
}
```

구성 파일에서 추가 프로세스 지시문을 사용하여 모든 프로세스 또는 특정 명명된 프로세스에 대해 이러한 값을 수정할 수 있습니다.
Nextflow는 선택한 실행자에 대한 적절한 지침으로 변환합니다.

하지만 어떤 값을 사용해야 하는지 어떻게 알 수 있나요?

### 5.1. 워크플로우를 실행하여 리소스 사용률 보고서 생성

??? example "시나리오"

    프로세스에 필요한 메모리나 CPU의 양을 미리 알지 못하며 리소스를 낭비하거나 작업이 종료되는 것을 피하고 싶습니다.

프로세스가 필요로 할 가능성이 있는 CPU와 메모리의 양을 미리 알지 못하는 경우 리소스 프로파일링을 수행할 수 있습니다. 즉, 일부 기본 할당으로 워크플로우를 실행하고 각 프로세스가 사용한 양을 기록한 다음 거기에서 기본 할당을 조정하는 방법을 추정합니다.

편리하게도 Nextflow에는 이를 수행하기 위한 내장 도구가 포함되어 있으며 요청 시 기꺼이 보고서를 생성합니다.

이를 위해 명령줄에 `-with-report <filename>.html`을 추가합니다.

```bash
nextflow run 3-main.nf -with-report report-config-1.html
```

보고서는 html 파일이며 다운로드하여 브라우저에서 열 수 있습니다. 왼쪽의 파일 탐색기에서 마우스 오른쪽 버튼을 클릭하고 `Show preview`를 클릭하여 교육 환경에서 볼 수도 있습니다.

보고서를 살펴보고 리소스 조정 기회를 식별할 수 있는지 몇 분 정도 시간을 내세요.
할당된 것의 백분율로 사용률 결과를 보여주는 탭을 클릭해야 합니다.

사용 가능한 모든 기능에 대한 문서는 [Reports](https://nextflow.io/docs/latest/reports.html)를 참조하세요.

### 5.2. 모든 프로세스에 대한 리소스 할당 설정

프로파일링 결과 교육 워크플로우의 프로세스가 매우 가볍다는 것을 보여주므로 프로세스당 기본 메모리 할당을 1GB로 줄이겠습니다.

파이프라인 매개변수 섹션 앞에 `nextflow.config` 파일에 다음을 추가하세요:

```groovy title="nextflow.config" linenums="4"
/*
* Process settings
*/
process {
    memory = 1.GB
}
```

이렇게 하면 소비하는 컴퓨팅 양을 줄이는 데 도움이 됩니다.

### 5.3. 특정 프로세스에 대한 리소스 할당 설정

동시에 개별 프로세스에 대한 할당을 조정하는 방법을 시연하기 위해 `cowpy` 프로세스가 다른 프로세스보다 더 많은 리소스를 필요로 한다고 가정하겠습니다.

=== "후"

    ```groovy title="nextflow.config" linenums="4" hl_lines="6-9"
    /*
    * Process settings
    */
    process {
        memory = 1.GB
        withName: 'cowpy' {
            memory = 2.GB
            cpus = 2
        }
    }
    ```

=== "전"

    ```groovy title="nextflow.config" linenums="4"
    /*
    * Process settings
    */
    process {
        memory = 1.GB
    }
    ```

이 구성을 사용하면 모든 프로세스가 1GB의 메모리와 단일 CPU(암시된 기본값)를 요청하지만 `cowpy` 프로세스는 2GB와 2개의 CPU를 요청합니다.

!!! info

    CPU가 적은 머신이 있고 프로세스당 많은 수를 할당하면 프로세스 호출이 서로 뒤에 대기하는 것을 볼 수 있습니다.
    이는 Nextflow가 사용 가능한 것보다 더 많은 CPU를 요청하지 않도록 보장하기 때문입니다.

### 5.4. 업데이트된 구성으로 워크플로우 실행

구성 변경 전후의 성능을 비교할 수 있도록 프로파일링 보고서에 다른 파일 이름을 제공하여 시도해 봅시다.

```bash
nextflow run 3-main.nf -with-report report-config-2.html
```

이것은 매우 작은 워크로드이므로 실제 차이를 느끼지 못할 것이지만, 이것이 실제 워크플로우의 성능과 리소스 요구 사항을 분석하는 데 사용할 접근 방식입니다.

프로세스의 리소스 요구 사항이 다를 때 매우 유용합니다. 추측이 아닌 실제 데이터를 기반으로 각 프로세스에 대해 설정한 리소스 할당을 적절하게 조정할 수 있습니다.

!!! tip

    이것은 리소스 사용을 최적화하기 위해 할 수 있는 것의 작은 맛보기일 뿐입니다.
    Nextflow 자체에는 리소스 제한으로 인해 실패한 작업을 재시도하는 정말 멋진 [동적 재시도 로직](https://nextflow.io/docs/latest/process.html#dynamic-task-resources)이 내장되어 있습니다.
    또한 Seqera Platform은 리소스 할당을 자동으로 최적화하기 위한 AI 기반 도구도 제공합니다.

### 5.5. 리소스 제한 추가

사용 중인 컴퓨팅 실행자 및 컴퓨팅 인프라에 따라 할당할 수 있는(또는 해야 하는) 것에 대한 몇 가지 제약이 있을 수 있습니다.
예를 들어, 클러스터에서 특정 제한 내에 머물도록 요구할 수 있습니다.

`resourceLimits` 지시문을 사용하여 관련 제한을 설정할 수 있습니다. 프로세스 블록에 단독으로 있을 때 구문은 다음과 같습니다:

```groovy title="구문 예제"
process {
    resourceLimits = [
        memory: 750.GB,
        cpus: 200,
        time: 30.d
    ]
}
```

Nextflow는 지정한 실행자에 따라 이러한 값을 적절한 지침으로 변환합니다.

교육 환경에서 관련 인프라에 액세스할 수 없으므로 실행하지 않을 것입니다.
그러나 이러한 제한을 초과하는 리소스 할당으로 워크플로우를 실행한 다음 `.command.run` 스크립트 파일에서 `sbatch` 명령을 조회하면 실제로 실행자에게 전송되는 요청이 `resourceLimits`에 지정된 값으로 제한되는 것을 볼 수 있습니다.

??? info "기관 참조 구성"

    nf-core 프로젝트는 광범위한 HPC 및 클라우드 실행자를 다루는 전 세계 다양한 기관에서 공유한 [구성 파일 모음](https://nf-co.re/configs/)을 편집했습니다.

    이러한 공유 구성은 그곳에서 일하는 사람들이 기관의 구성을 즉시 활용할 수 있고, 자체 인프라에 대한 구성을 개발하려는 사람들을 위한 모델로서 가치가 있습니다.

### 핵심 정리

리소스 사용률을 평가하기 위한 프로파일링 보고서를 생성하는 방법과 모든 프로세스 및/또는 개별 프로세스에 대한 리소스 할당을 수정하는 방법, HPC에서 실행하기 위한 리소스 제한을 설정하는 방법을 알게 되었습니다.

### 다음 단계

사전 설정된 구성 프로필을 설정하고 런타임에 전환하는 방법을 학습합니다.

---

## 6. 프로필을 사용하여 사전 설정된 구성 간 전환

??? example "시나리오"

    개발을 위해 노트북에서 파이프라인을 실행하는 것과 프로덕션 실행을 위해 기관의 HPC에서 실행하는 것 사이를 정기적으로 전환합니다.
    환경을 전환할 때마다 구성 설정을 수동으로 변경하는 것이 지겹습니다.

작업 중인 프로젝트나 사용 중인 컴퓨팅 환경에 따라 파이프라인 구성을 사용자 정의할 수 있는 여러 가지 방법을 보여드렸습니다.

사용 중인 컴퓨팅 인프라에 따라 대체 설정 간에 전환하고 싶을 수 있습니다. 예를 들어, 노트북에서 로컬로 소규모 테스트를 개발하고 실행한 다음 HPC 또는 클라우드에서 전체 규모의 워크로드를 실행할 수 있습니다.

Nextflow를 사용하면 구성 파일 자체를 수정하는 대신 명령줄 인수를 사용하여 런타임에 선택할 수 있는 다양한 구성을 설명하는 [**프로필**](https://nextflow.io/docs/latest/config.html#profiles)을 원하는 수만큼 설정할 수 있습니다.

### 6.1. 로컬 개발과 HPC 실행 간 전환을 위한 프로필 생성

두 가지 대체 프로필을 설정하겠습니다. 하나는 Docker 컨테이너를 사용할 일반 컴퓨터에서 소규모 로드를 실행하기 위한 것이고, 다른 하나는 Conda 패키지를 사용할 Slurm 스케줄러가 있는 대학 HPC에서 실행하기 위한 것입니다.

#### 6.1.1. 프로필 설정

파이프라인 매개변수 섹션 뒤, 출력 설정 앞에 `nextflow.config` 파일에 다음을 추가하세요:

```groovy title="nextflow.config" linenums="24"
/*
* Profiles
*/
profiles {
    my_laptop {
        process.executor = 'local'
        docker.enabled = true
    }
    univ_hpc {
        process.executor = 'slurm'
        conda.enabled = true
        process.resourceLimits = [
            memory: 750.GB,
            cpus: 200,
            time: 30.d
        ]
    }
}
```

대학 HPC의 경우 리소스 제한도 지정하고 있습니다.

#### 6.1.2. 프로필로 워크플로우 실행

Nextflow 명령줄에서 프로필을 지정하려면 `-profile` 인수를 사용합니다.

`my_laptop` 구성으로 워크플로우를 실행해 봅시다.

```bash
nextflow run 3-main.nf -profile my_laptop
```

??? success "명령 출력"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `3-main.nf` [gigantic_brazil] DSL2 - revision: ede9037d02

    executor >  local (8)
    [58/da9437] sayHello (3)       | 3 of 3 ✔
    [35/9cbe77] convertToUpper (2) | 3 of 3 ✔
    [67/857d05] collectGreetings   | 1 of 1 ✔
    [37/7b51b5] cowpy              | 1 of 1 ✔
    ```

보시다시피 이를 통해 런타임에 구성 간을 매우 편리하게 전환할 수 있습니다.

!!! warning

    Slurm 스케줄러에 액세스할 수 없으므로 `univ_hpc` 프로필은 교육 환경에서 제대로 실행되지 않습니다.

향후 이러한 것과 항상 함께 발생하는 다른 구성 요소를 발견하면 해당 프로필에 간단히 추가할 수 있습니다.
함께 그룹화하려는 다른 구성 요소가 있는 경우 추가 프로필을 만들 수도 있습니다.

### 6.2. 테스트 매개변수 프로필 생성

??? example "시나리오"

    다른 사람들이 자신의 입력 데이터를 수집하지 않고도 파이프라인을 빠르게 시도할 수 있기를 원합니다.

프로필은 인프라 구성만을 위한 것이 아닙니다.
워크플로우 매개변수의 기본값을 설정하는 데도 사용할 수 있으므로 다른 사람들이 적절한 입력 값을 직접 수집하지 않고도 워크플로우를 더 쉽게 시도할 수 있습니다.
이것을 매개변수 파일을 사용하는 대안으로 간주할 수 있습니다.

#### 6.2.1. 프로필 설정

이 컨텍스트에서 기본값을 표현하는 구문은 `test`라는 이름의 프로필에 대해 다음과 같습니다:

```groovy title="구문 예제"
    test {
        params.<parameter1>
        params.<parameter2>
        ...
    }
```

워크플로우에 대한 테스트 프로필을 추가하면 `profiles` 블록은 다음과 같이 됩니다:

```groovy title="nextflow.config" linenums="24"
/*
* Profiles
*/
profiles {
    my_laptop {
        process.executor = 'local'
        docker.enabled = true
    }
    univ_hpc {
        process.executor = 'slurm'
        conda.enabled = true
        process.resourceLimits = [
            memory: 750.GB,
            cpus: 200,
            time: 30.d
        ]
    }
    test {
        params.input = 'data/greetings.csv'
        params.batch = 'test'
        params.character = 'dragonandcow'
    }
}
```

기술 구성 프로필과 마찬가지로 원하는 임의의 이름으로 매개변수를 지정하는 여러 다른 프로필을 설정할 수 있습니다.

#### 6.2.2. 테스트 프로필로 워크플로우를 로컬에서 실행

편리하게도 프로필은 상호 배타적이지 않으므로 다음 구문 `-profile <profile1>,<profile2>`(프로필 수에 관계없이)를 사용하여 명령줄에서 여러 프로필을 지정할 수 있습니다.

동일한 구성 요소에 대한 값을 설정하고 동일한 구성 파일에 설명된 프로필을 결합하면 Nextflow는 마지막으로 읽은 값(_즉,_ 파일에서 나중에 오는 것)을 사용하여 충돌을 해결합니다.
충돌하는 설정이 다른 구성 소스에 설정된 경우 기본 [우선순위](https://www.nextflow.io/docs/latest/config.html)가 적용됩니다.

이전 명령에 테스트 프로필을 추가해 봅시다:

```bash
nextflow run 3-main.nf -profile my_laptop,test
```

??? success "명령 출력"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `3-main.nf` [jovial_coulomb] DSL2 - revision: 46a6763141

    executor >  local (8)
    [9b/687cdc] sayHello (2)       | 3 of 3 ✔
    [ca/552187] convertToUpper (3) | 3 of 3 ✔
    [e8/83e306] collectGreetings   | 1 of 1 ✔
    [fd/e84fa9] cowpy              | 1 of 1 ✔
    ```

이렇게 하면 가능한 경우 Docker를 사용하고 `results_config/test` 아래에 출력을 생성하며, 이번에는 캐릭터가 코미디 듀오 `dragonandcow`입니다.

??? abstract "파일 내용"

    ```console title="results_config/test/"
     _________
    / HOLà    \
    | HELLO   |
    \ BONJOUR /
    ---------
                \                    ^    /^
                  \                  / \  // \
                  \   |\___/|      /   \//  .\
                    \  /O  O  \__  /    //  | \ \           *----*
                      /     /  \/_/    //   |  \  \          \   |
                      \@___\@`    \/_   //    |   \   \         \/\ \
                    0/0/|       \/_ //     |    \    \         \ \
                0/0/0/0/|        \///      |     \     \       | |
              0/0/0/0/0/_|_ /   (  //       |      \     _\     |  /
          0/0/0/0/0/0/`/,_ _ _/  ) ; -.    |    _ _\.-~       /   /
                      ,-}        _      *-.|.-~-.           .~    ~
      \     \__/        `/\      /                 ~-. _ .-~      /
      \____(oo)           *.   }            {                   /
      (    (--)          .----~-.\        \-`                 .~
      //__\\  \__ Ack!   ///.----..<        \             _ -~
      //    \\               ///-._ _ _ _ _ _ _{^ - - - - ~
    ```

즉, 워크플로우 코드와 함께 테스트 데이터 파일을 배포하는 한 누구나 명령줄이나 매개변수 파일을 통해 자신의 입력을 제공하지 않고도 워크플로우를 빠르게 시도할 수 있습니다.

!!! tip

    외부에 저장된 더 큰 파일에 대한 URL을 가리킬 수 있습니다.
    Nextflow는 열린 연결이 있는 한 자동으로 다운로드합니다.

    자세한 내용은 사이드 퀘스트 [Working with Files](../side_quests/working_with_files.md)를 참조하세요.

### 6.3. `nextflow config`를 사용하여 해결된 구성 확인

위에서 언급했듯이 때때로 동일한 매개변수가 결합하려는 프로필에서 다른 값으로 설정될 수 있습니다.
그리고 더 일반적으로 구성 요소를 저장할 수 있는 곳이 많으며 때때로 동일한 속성이 다른 곳에서 다른 값으로 설정될 수 있습니다.

Nextflow는 충돌을 해결하기 위해 설정된 [우선순위](https://nextflow.io/docs/latest/config.html#configuration-file)를 적용하지만 직접 결정하기는 까다로울 수 있습니다.
그리고 충돌하는 것이 없더라도 구성할 수 있는 모든 가능한 위치를 조회하는 것은 지루할 수 있습니다.

다행히 Nextflow에는 전체 프로세스를 자동화할 수 있는 `config`라는 편리한 유틸리티 도구가 포함되어 있습니다.

`config` 도구는 현재 작업 디렉토리의 모든 내용을 탐색하고 구성 파일을 수집하여 Nextflow가 워크플로우를 실행하는 데 사용할 완전히 해결된 구성을 생성합니다.
이를 통해 아무것도 시작하지 않고도 어떤 설정이 사용될지 알아낼 수 있습니다.

#### 6.3.1. 기본 구성 해결

기본적으로 적용될 구성을 해결하려면 이 명령을 실행하세요.

```bash
nextflow config
```

??? success "명령 출력"

    ```groovy
    params {
      input = 'data/greetings.csv'
      batch = 'batch'
      character = 'turkey'
    }

    docker {
      enabled = false
    }

    conda {
      enabled = true
    }

    process {
      memory = '1 GB'
      withName:cowpy {
          memory = '2 GB'
          cpus = 2
      }
    }

    outputDir = 'results_config/batch'

    workflow {
      output {
          mode = 'copy'
      }
    }
    ```

이것은 명령줄에서 추가로 지정하지 않으면 얻을 수 있는 기본 구성을 보여줍니다.

#### 6.3.2. 특정 설정이 활성화된 구성 해결

하나 이상의 프로필을 활성화하거나 매개변수 파일을 로드하는 등의 명령줄 매개변수를 제공하면 명령이 추가로 이를 고려합니다.

```bash
nextflow config -profile my_laptop,test
```

??? success "명령 출력"

    ```groovy
    params {
      input = 'data/greetings.csv'
      batch = 'test'
      character = 'dragonandcow'
    }

    docker {
      enabled = true
    }

    conda {
      enabled = true
    }

    process {
      memory = '1 GB'
      withName:cowpy {
          memory = '2 GB'
          cpus = 2
      }
      executor = 'local'
    }

    outputDir = 'results_config/test'

    workflow {
      output {
          mode = 'copy'
      }
    }
    ```

이것은 여러 구성 계층이 포함된 복잡한 프로젝트에 특히 유용합니다.

### 핵심 정리

프로필을 사용하여 최소한의 번거로움으로 런타임에 사전 설정된 구성을 선택하는 방법을 알게 되었습니다.
더 일반적으로 다양한 컴퓨팅 플랫폼에 맞게 워크플로우 실행을 구성하고 분석의 재현성을 향상시키는 방법을 알게 되었습니다.

### 다음 단계

GitHub와 같은 원격 저장소에서 직접 파이프라인을 실행하는 방법을 학습합니다.

---

## 7. 원격 저장소에서 파이프라인 실행

??? example "시나리오"

    코드를 직접 다운로드하고 관리하지 않고도 nf-core와 같은 잘 확립된 파이프라인을 실행하고 싶습니다.

지금까지 현재 디렉토리에 있는 워크플로우 스크립트를 실행해 왔습니다.
실제로는 GitHub와 같은 원격 저장소에 저장된 파이프라인을 실행하고 싶을 때가 많습니다.

Nextflow는 이를 간단하게 만듭니다. 수동으로 다운로드하지 않고도 Git 저장소 URL에서 직접 파이프라인을 실행할 수 있습니다.

### 7.1. GitHub에서 파이프라인 실행

원격 파이프라인을 실행하기 위한 기본 구문은 `nextflow run <repository>`이며, 여기서 `<repository>`는 `nextflow-io/hello`와 같은 GitHub 저장소 경로, 전체 URL 또는 GitLab, Bitbucket 또는 기타 Git 호스팅 서비스의 경로일 수 있습니다.

공식 Nextflow "hello" 데모 파이프라인을 실행해 보세요:

```bash
nextflow run nextflow-io/hello
```

??? success "명령 출력"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Pulling nextflow-io/hello ...
     downloaded from https://github.com/nextflow-io/hello.git
    Launching `https://github.com/nextflow-io/hello` [sleepy_swanson] DSL2 - revision: 2ce0b0e294 [master]

    executor >  local (4)
    [ba/08236d] sayHello (4) [100%] 4 of 4 ✔
    Ciao world!

    Hello world!

    Bonjour world!

    Hola world!
    ```

원격 파이프라인을 처음 실행하면 Nextflow가 다운로드하여 로컬에 캐시합니다.
후속 실행은 명시적으로 업데이트를 요청하지 않는 한 캐시된 버전을 사용합니다.

### 7.2. 재현성을 위한 버전 지정

기본적으로 Nextflow는 기본 브랜치의 최신 버전을 실행합니다.
`-r` 플래그를 사용하여 특정 버전(태그), 브랜치 또는 커밋을 지정할 수 있습니다:

```bash
nextflow run nextflow-io/hello -r v1.3
```

??? success "명령 출력"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `https://github.com/nextflow-io/hello` [sick_carson] DSL2 - revision: 2ce0b0e294 [v1.3]

    executor >  local (4)
    [61/e11f77] sayHello (4) [100%] 4 of 4 ✔
    Ciao world!

    Bonjour world!

    Hello world!

    Hola world!
    ```

정확한 버전을 지정하는 것은 재현성에 필수적입니다.

### 핵심 정리

GitHub 및 기타 원격 저장소에서 직접 파이프라인을 실행하는 방법과 재현성을 위해 버전을 지정하는 방법을 알게 되었습니다.

### 다음 단계

스스로에게 큰 박수를 보내세요!
Nextflow 파이프라인을 실행하고 관리하는 데 필요한 모든 것을 알게 되었습니다.

이것으로 이 과정을 마치지만, 계속 학습하고 싶다면 두 가지 주요 권장 사항이 있습니다:

- 자체 파이프라인 개발에 대해 더 깊이 파고들고 싶다면 채널과 연산자에 대해 훨씬 더 자세히 다루는 초보자를 위한 과정인 [Hello Nextflow](../hello_nextflow/index.md)를 살펴보세요.
- 코드를 더 깊이 파고들지 않고 Nextflow 파이프라인 실행에 대해 계속 학습하고 싶다면 매우 인기 있는 [nf-core](https://nf-co.re/) 프로젝트에서 파이프라인을 찾고 실행하기 위한 도구를 소개하는 [Hello nf-core](../hello_nf-core/index.md)의 첫 번째 파트를 살펴보세요.

즐거운 시간 되세요!

---

## 퀴즈

<quiz>
매개변수 값이 워크플로우 파일과 `nextflow.config` 모두에 설정된 경우 어느 것이 우선합니까?
- [ ] 워크플로우 파일 값
- [x] 구성 파일 값
- [ ] 처음 발견된 값
- [ ] 오류가 발생합니다

자세히 알아보기: [1.1. `nextflow.config`에 값 설정](#11-nextflowconfig에-값-설정)
</quiz>

<quiz>
워크플로우 파일과 구성 파일에서 매개변수 기본값을 설정하는 구문 차이는 무엇입니까?
- [ ] 동일한 구문을 사용합니다
- [x] 워크플로우는 타입이 지정된 선언(`#!groovy param: Type = value`)을 사용하고, 구성은 할당(`#!groovy param = value`)을 사용합니다
- [ ] 구성은 타입이 지정된 선언을 사용하고, 워크플로우는 할당을 사용합니다
- [ ] 구성 파일만 기본값을 설정할 수 있습니다

자세히 알아보기: [1.1. `nextflow.config`에 값 설정](#11-nextflowconfig에-값-설정)
</quiz>

<quiz>
워크플로우를 실행할 때 매개변수 파일을 지정하는 방법은 무엇입니까?
- [ ] `--params params.yaml`
- [ ] `-config params.yaml`
- [x] `-params-file params.yaml`
- [ ] `--input-params params.yaml`

자세히 알아보기: [1.3. 매개변수 파일 사용](#13-매개변수-파일-사용)
</quiz>

<quiz>
`outputDir` 구성 옵션은 무엇을 제어합니까?
- [ ] 작업 디렉토리의 위치
- [x] 워크플로우 출력이 게시되는 기본 경로
- [ ] 로그 파일의 디렉토리
- [ ] 모듈 파일의 위치

자세히 알아보기: [2.1. outputDir 디렉토리 이름 사용자 정의](#21-outputdir-디렉토리-이름-사용자-정의)
</quiz>

<quiz>
출력 경로 구성에서 프로세스 이름을 동적으로 참조하는 방법은 무엇입니까?
- [ ] `#!groovy ${processName}`
- [ ] `#!groovy path "<process>.name"`
- [x] `#!groovy path { <process>.name }`
- [ ] `@processName`

자세히 알아보기: [2.2. 프로세스별로 출력 구성](#22-프로세스별로-출력-구성)
</quiz>

<quiz>
Docker와 Conda가 모두 활성화되어 있고 프로세스에 두 지시문이 모두 있는 경우 어느 것이 우선합니까?
- [x] Docker (컨테이너)
- [ ] Conda
- [ ] 프로세스에서 먼저 정의된 것
- [ ] 오류가 발생합니다

자세히 알아보기: [3. 소프트웨어 패키징 기술 선택](#3-소프트웨어-패키징-기술-선택)
</quiz>

<quiz>
Nextflow의 기본 실행자는 무엇입니까?
- [x] `local`
- [ ] `slurm`
- [ ] `kubernetes`
- [ ] `aws`

자세히 알아보기: [4. 실행 플랫폼 선택](#4-실행-플랫폼-선택)
</quiz>

<quiz>
리소스 사용률 보고서를 생성하는 명령은 무엇입니까?
- [ ] `nextflow run workflow.nf -with-metrics`
- [ ] `nextflow run workflow.nf -with-stats`
- [x] `nextflow run workflow.nf -with-report report.html`
- [ ] `nextflow run workflow.nf -profile report`

자세히 알아보기: [5.1. 워크플로우를 실행하여 리소스 사용률 보고서 생성](#51-워크플로우를-실행하여-리소스-사용률-보고서-생성)
</quiz>

<quiz>
구성 파일에서 `cowpy`라는 특정 프로세스에 대한 리소스 요구 사항을 설정하는 방법은 무엇입니까?
- [ ] `#!groovy cowpy.memory = '2.GB'`
- [ ] `#!groovy process.cowpy.memory = '2.GB'`
- [x] `#!groovy process { withName: 'cowpy' { memory = '2.GB' } }`
- [ ] `#!groovy resources.cowpy.memory = '2.GB'`

자세히 알아보기: [5.3. 특정 프로세스에 대한 리소스 할당 설정](#53-특정-프로세스에-대한-리소스-할당-설정)
</quiz>

<quiz>
`resourceLimits` 지시문은 무엇을 합니까?
- [ ] 최소 리소스 요구 사항을 설정합니다
- [ ] 프로세스에 리소스를 할당합니다
- [x] 요청할 수 있는 최대 리소스를 제한합니다
- [ ] 실시간으로 리소스 사용을 모니터링합니다

자세히 알아보기: [5.5. 리소스 제한 추가](#55-리소스-제한-추가)
</quiz>

<quiz>
단일 명령에서 여러 프로필을 지정하는 방법은 무엇입니까?
- [ ] `-profile profile1 -profile profile2`
- [ ] `-profiles profile1,profile2`
- [x] `-profile profile1,profile2`
- [ ] `--profile profile1 --profile profile2`

자세히 알아보기: [6. 프로필을 사용하여 사전 설정된 구성 간 전환](#6-프로필을-사용하여-사전-설정된-구성-간-전환)
</quiz>

<quiz>
Nextflow가 사용할 완전히 해결된 구성을 보여주는 명령은 무엇입니까?
- [ ] `nextflow show-config`
- [ ] `nextflow settings`
- [x] `nextflow config`
- [ ] `nextflow resolve`

자세히 알아보기: [6.3. nextflow config를 사용하여 해결된 구성 확인](#63-nextflow-config를-사용하여-해결된-구성-확인)
</quiz>

<quiz>
프로필은 무엇에 사용할 수 있습니까? (해당하는 모든 항목 선택)
- [x] 인프라별 설정 정의(실행자, 컨테이너)
- [x] 다양한 환경에 대한 리소스 제한 설정
- [x] 쉬운 워크플로우 테스트를 위한 테스트 매개변수 제공
- [ ] 새 프로세스 정의

자세히 알아보기: [6. 프로필을 사용하여 사전 설정된 구성 간 전환](#6-프로필을-사용하여-사전-설정된-구성-간-전환)
</quiz>
