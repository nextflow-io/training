# 파트 6: Hello Config

<div class="video-wrapper">
  <iframe width="560" height="315" src="https://www.youtube.com/embed/FcZTiE25TeA?si=y8lAedhEHWaTV4zd&amp;list=PLPZ8WHdZGxmWKozQuzr27jyMGqp9kElVK&amp;cc_load_policy=1&amp;cc_lang_pref=ko" title="YouTube video player" frameborder="0" allow="accelerometer; autoplay; clipboard-write; encrypted-media; gyroscope; picture-in-picture; web-share" referrerpolicy="strict-origin-when-cross-origin" allowfullscreen></iframe>
</div>

/// caption
:fontawesome-brands-youtube:{ .youtube } Nextflow YouTube 채널에서 [전체 재생목록](https://youtube.com/playlist?list=PLPZ8WHdZGxmWKozQuzr27jyMGqp9kElVK&si=eF7cLR62goy-lc6n)을 확인하십시오.

:green_book: 비디오 스크립트는 [여기](./transcripts/06_hello_config.md)에서 확인하실 수 있습니다.
///

이 섹션에서는 _워크플로우 코드 자체를 한 줄도 변경하지 않고_ 동작을 사용자 정의하고, 다른 환경에 적응시키며, 리소스 사용을 최적화할 수 있도록 Nextflow 파이프라인의 구성을 설정하고 관리하는 방법을 다룹니다.

이를 수행하는 여러 방법이 있으며, 조합하여 사용할 수 있고 구성 문서에 설명된 [우선순위](https://nextflow.io/docs/latest/config.html)에 따라 해석됩니다.

이 파트에서는 파트 5: Hello Containers에서 이미 접한 가장 간단하고 일반적인 구성 파일 메커니즘인 [`nextflow.config`](https://nextflow.io/docs/latest/config.html) 파일을 보여드릴 것입니다.

Process 지시문, 실행자, 프로필 및 매개변수 파일과 같은 Nextflow 구성의 필수 구성 요소를 살펴볼 것입니다.
이러한 구성 옵션을 효과적으로 활용하는 방법을 배우면 파이프라인의 유연성, 확장성 및 성능을 향상시킬 수 있습니다.

??? info "이 섹션부터 시작하는 방법"

    이 섹션은 [Hello Nextflow](./index.md) 과정의 파트 1-5를 완료하고 완전히 작동하는 파이프라인이 있다고 가정합니다.

    이 시점부터 과정을 시작하는 경우 solutions에서 `modules` 디렉토리와 `nextflow.config` 파일을 복사해야 합니다:

    ```bash
    cp -r solutions/5-hello-containers/modules .
    cp solutions/5-hello-containers/nextflow.config .
    ```

    `nextflow.config` 파일에는 Docker 컨테이너 사용을 활성화하는 `docker.enabled = true` 줄이 포함되어 있습니다.

    Hello 파이프라인에 익숙하지 않거나 상기가 필요하면 [이 정보 페이지](../info/hello_pipeline.md)를 참조하십시오.

---

## 0. 준비 운동: `hello-config.nf` 실행

시작점으로 워크플로우 스크립트 `hello-config.nf`를 사용할 것입니다.
이 스크립트는 이 교육 과정의 파트 5를 완료하여 생성된 스크립트와 동일하지만, 출력 대상을 변경했습니다:

```groovy title="hello-config.nf" linenums="37" hl_lines="3 7 11 15"
output {
    first_output {
        path 'hello_config/intermediates'
        mode 'copy'
    }
    uppercased {
        path 'hello_config/intermediates'
        mode 'copy'
    }
    collected {
        path 'hello_config/intermediates'
        mode 'copy'
    }
    batch_report {
        path 'hello_config'
        mode 'copy'
    }
    cowpy_art {
        path 'hello_config'
        mode 'copy'
    }
}
```

변경을 시작하기 전에 모든 것이 제대로 작동하는지 확인하기 위해 스크립트를 한 번 실행하십시오:

```bash
nextflow run hello-config.nf
```

??? success "명령 출력"

    ```console
     N E X T F L O W   ~  version 25.10.2

    Launching `hello-config.nf` [nice_escher] DSL2 - revision: d5dfdc9872

    executor > local (7)
    [6a/bc46a6] sayHello (2) [100%] 3 of 3 ✔
    [33/67bc48] convertToUpper (3) [100%] 3 of 3 ✔
    [b5/de03ba] collectGreetings [100%] 1 of 1 ✔
    [98/c6b57b] cowpy              | 1 of 1 ✔
    ```

이전과 마찬가지로 `output` 블록에 지정된 디렉토리(`results/hello_config/`)에서 출력 파일을 찾을 수 있습니다.

??? abstract "디렉토리 내용"

    ```console
    results/hello_config/
    ├── cowpy-COLLECTED-batch-output.txt
    ├── intermediates
    │   ├── Bonjour-output.txt
    │   ├── COLLECTED-batch-output.txt
    │   ├── Hello-output.txt
    │   ├── Holà-output.txt
    │   ├── UPPER-Bonjour-output.txt
    │   ├── UPPER-Hello-output.txt
    │   └── UPPER-Holà-output.txt
    └── batch-report.txt
    ```

최종 ASCII 아트 출력은 `results/hello_config/` 디렉토리에 `cowpy-COLLECTED-batch-output.txt`라는 이름으로 있습니다.

??? abstract "파일 내용"

    ```console title="results/hello_config/cowpy-COLLECTED-batch-output.txt"
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

이것이 정상적으로 작동했다면 파이프라인을 구성하는 방법을 배울 준비가 되었습니다.

---

## 1. 워크플로우 입력 매개변수 관리

지금까지 작업해 온 것의 확장인 구성의 한 측면부터 시작할 것입니다: 입력 매개변수의 관리입니다.

현재 워크플로우는 명령줄을 통해 여러 매개변수 값을 허용하도록 설정되어 있으며, 기본값은 워크플로우 스크립트 자체의 `params` 블록에 설정되어 있습니다.
그러나 명령줄에서 매개변수를 지정하거나 원본 스크립트 파일을 수정하지 않고 해당 기본값을 재정의하고 싶을 수 있습니다.

이를 수행하는 여러 방법이 있습니다. 매우 일반적으로 사용되는 세 가지 기본 방법을 보여드리겠습니다.

### 1.1. 기본값을 `nextflow.config`로 이동

이것은 가장 간단한 접근 방식이지만, 실행할 때마다 편집하고 싶지 않은 기본 `nextflow.config` 파일이므로 아마도 가장 유연하지 않습니다.
그러나 워크플로우에서 매개변수를 _선언하는_ 것(확실히 거기에 속함)과 구성 파일에 더 적합한 *기본값*을 제공하는 것의 관심사를 분리하는 이점이 있습니다.

두 단계로 수행해 봅시다.

#### 1.1.1. 구성 파일에 `params` 블록 생성

`nextflow.config` 파일에서 다음 코드 변경을 수행하십시오:

=== "수정 후"

    ```groovy title="nextflow.config" linenums="1" hl_lines="3-10"
    docker.enabled = true

    /*
    * 파이프라인 매개변수
    */
    params {
        input = 'data/greetings.csv'
        batch = 'batch'
        character = 'turkey'
    }
    ```

=== "수정 전"

    ```groovy title="nextflow.config" linenums="1"
    docker.enabled = true
    ```

워크플로우에서 구성 파일로 `params` 블록을 단순히 복사한 것이 아닙니다.
구문이 약간 다릅니다.
워크플로우 파일에서는 타입이 지정된 선언입니다.
구성에서는 값 할당입니다.

기술적으로 이것만으로도 워크플로우 파일에 여전히 지정된 기본값을 재정의하기에 충분합니다.
예를 들어 캐릭터를 수정하고 워크플로우를 실행하여 구성 파일에 설정된 값이 워크플로우 파일에 설정된 값을 재정의하는지 확인할 수 있습니다.

그러나 구성을 완전히 구성 파일로 이동하는 정신으로, 워크플로우 파일에서 해당 값을 완전히 제거합시다.

#### 1.1.2. 워크플로우 파일의 `params` 블록에서 값 제거

`hello-config.nf` 워크플로우 파일에서 다음 코드 변경을 수행하십시오:

=== "수정 후"

    ```groovy title="hello-config.nf" linenums="9" hl_lines="5-7"
    /*
    * 파이프라인 매개변수
    */
    params {
        input: Path
        batch: String
        character: String
    }
    ```

=== "수정 전"

    ```groovy title="hello-config.nf" linenums="9" hl_lines="5-7"
    /*
    * 파이프라인 매개변수
    */
    params {
        input: Path = 'data/greetings.csv'
        batch: String = 'batch'
        character: String = 'turkey'
    }
    ```

이제 워크플로우 파일 자체는 이러한 매개변수에 대한 기본값을 설정하지 않습니다.

#### 1.1.3. 파이프라인 실행

올바르게 작동하는지 테스트해 봅시다.

```bash
nextflow run hello-config.nf
```

??? success "명령 출력"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `hello-config.nf` [disturbed_einstein] DSL2 - revision: ede9037d02

    executor >  local (8)
    [f0/35723c] sayHello (2)       | 3 of 3 ✔
    [40/3efd1a] convertToUpper (3) | 3 of 3 ✔
    [17/e97d32] collectGreetings   | 1 of 1 ✔
    [98/c6b57b] cowpy              | 1 of 1 ✔
    ```

이것은 이전과 동일한 출력을 생성합니다.

최종 ASCII 아트 출력은 이전과 마찬가지로 `results/hello_config/` 디렉토리에 `cowpy-COLLECTED-batch-output.txt`라는 이름으로 있습니다.

??? abstract "파일 내용"

    ```console title="results/hello_config/cowpy-COLLECTED-batch-output.txt"
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

기능적으로 이 이동은 아무것도 변경하지 않았지만, 개념적으로 구성 파일에 기본값을 설정하는 것이 약간 더 깔끔합니다.

### 1.2. 실행별 구성 파일 사용

좋습니다만, 때때로 메인 구성 파일을 건드리지 않고 다른 기본값으로 일부 임시 실험을 실행하고 싶을 수 있습니다.
실험을 위한 작업 디렉토리로 사용할 하위 디렉토리에 새 `nextflow.config` 파일을 생성하여 이를 수행할 수 있습니다.

#### 1.2.1. 빈 구성으로 작업 디렉토리 생성

새 디렉토리를 생성하고 그곳으로 이동하는 것으로 시작합시다:

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

이제 새 파일을 열고 사용자 정의하려는 매개변수를 추가하십시오:

```groovy title="tux-run/nextflow.config" linenums="1"
params {
    input = '../data/greetings.csv'
    batch = 'experiment'
    character = 'tux'
}
```

입력 파일 경로는 디렉토리 구조를 반영해야 합니다.

#### 1.2.3. 파이프라인 실행

이제 새 작업 디렉토리 내에서 파이프라인을 실행할 수 있습니다.
경로를 적절하게 조정해야 합니다!

```bash
nextflow run ../hello-config.nf
```

??? success "명령 출력"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `../hello-config.nf` [trusting_escher] DSL2 - revision: 356df0818d

    executor >  local (8)
    [59/b66913] sayHello (2)       [100%] 3 of 3 ✔
    [ad/f06364] convertToUpper (3) [100%] 3 of 3 ✔
    [10/714895] collectGreetings   [100%] 1 of 1 ✔
    [88/3ece98] cowpy              [100%] 1 of 1 ✔
    ```

이렇게 하면 `tux-run/work/` 및 `tux-run/results/`를 포함하여 `tux-run/` 아래에 새 디렉토리 세트가 생성됩니다.

이 실행에서 Nextflow는 현재 디렉토리의 `nextflow.config`를 파이프라인 루트 디렉토리의 `nextflow.config`와 결합하여 기본 캐릭터(칠면조)를 tux 캐릭터로 재정의합니다.

최종 출력 파일에는 인사말을 말하는 tux 캐릭터가 포함되어야 합니다.

??? abstract "파일 내용"

    ```console title="tux-run/results/hello_config/cowpy-COLLECTED-experiment-output.txt"
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

이제 '정상적인' 구성을 수정하지 않고 실험할 공간이 생겼습니다.

!!! warning "경고"

    다음 섹션으로 이동하기 전에 이전 디렉토리로 다시 변경해야 합니다!

    ```bash
    cd ..
    ```

이제 매개변수 값을 설정하는 또 다른 유용한 방법을 살펴봅시다.

### 1.3. 매개변수 파일 사용

하위 디렉토리 접근 방식은 실험에 적합하지만, 약간의 설정이 필요하고 그에 따라 경로를 조정해야 합니다.
특정 값 세트로 파이프라인을 실행하거나 다른 사람이 최소한의 노력으로 실행할 수 있도록 하려는 경우 더 간단한 접근 방식이 있습니다.

Nextflow를 사용하면 YAML 또는 JSON 형식의 [매개변수 파일](https://nextflow.io/docs/latest/config.html#params-file)을 통해 매개변수를 지정할 수 있으므로, 예를 들어 대체 기본값 세트와 실행별 매개변수 값을 관리하고 배포하는 것이 매우 편리합니다.

#### 1.3.1. 예제 매개변수 파일 검사

이를 시연하기 위해 현재 디렉토리에 `test-params.yaml`이라는 예제 매개변수 파일을 제공합니다:

```yaml title="test-params.yaml" linenums="1"
input: "data/greetings.csv"
batch: "yaml"
character: "stegosaurus"
```

이 매개변수 파일에는 지정하려는 각 입력에 대한 키-값 쌍이 포함되어 있습니다.
구성 파일과 비교하면 등호(`=`) 대신 콜론(`:`)을 사용합니다.
구성 파일은 Groovy로 작성되고, 매개변수 파일은 YAML로 작성됩니다.

!!! info "정보"

    예제로 매개변수 파일의 JSON 버전도 제공하지만 여기서는 실행하지 않을 것입니다.
    직접 시도해 보십시오.

#### 1.3.2. 파이프라인 실행

이 매개변수 파일로 워크플로우를 실행하려면 기본 명령에 `-params-file <filename>`을 추가하기만 하면 됩니다.

```bash
nextflow run hello-config.nf -params-file test-params.yaml
```

??? success "명령 출력"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `hello-config.nf` [disturbed_sammet] DSL2 - revision: ede9037d02

    executor >  local (8)
    [f0/35723c] sayHello (2)       | 3 of 3 ✔
    [40/3efd1a] convertToUpper (3) | 3 of 3 ✔
    [17/e97d32] collectGreetings   | 1 of 1 ✔
    [98/c6b57b] cowpy              | 1 of 1 ✔
    ```

최종 출력 파일에는 인사말을 말하는 스테고사우루스 캐릭터가 포함되어야 합니다.

??? abstract "파일 내용"

    ```console title="results/hello_config/cowpy-COLLECTED-yaml-output.txt"
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

몇 개의 매개변수만 지정해야 할 때 매개변수 파일을 사용하는 것은 과도해 보일 수 있지만, 일부 파이프라인은 수십 개의 매개변수를 예상합니다.
그런 경우 매개변수 파일을 사용하면 방대한 명령줄을 입력하거나 워크플로우 스크립트를 수정하지 않고도 런타임에 매개변수 값을 제공할 수 있습니다.

또한 예를 들어 공동 작업자에게 매개변수 세트를 배포하거나 출판물의 보충 정보로 제공하는 것이 더 쉬워집니다.
이렇게 하면 다른 사람들이 작업을 더 재현할 수 있습니다.

### 핵심 정리

워크플로우 입력 관리를 위한 주요 구성 옵션을 활용하는 방법을 알게 되었습니다.

### 다음 단계

워크플로우 출력이 게시되는 위치와 방법을 관리하는 방법을 배웁니다.

---

## 2. 워크플로우 출력 관리

지금까지 워크플로우 수준 출력 선언에 대한 모든 경로를 하드코딩해 왔으며, 여러 출력을 추가하기 시작할 때 언급했듯이 약간의 반복이 있을 수 있습니다.

이를 더 유연하게 구성할 수 있는 몇 가지 일반적인 방법을 살펴봅시다.

### 2.1. `-output-dir`로 출력 디렉토리 사용자 정의

'게시된' 출력이 구성되는 방식을 제어할 때 두 가지 별개의 우선순위가 있습니다:

- 최상위 출력 디렉토리
- 이 디렉토리 내에서 파일이 구성되는 방식

지금까지 기본 최상위 디렉토리인 `results`를 사용해 왔습니다.
이를 사용자 정의하는 것부터 시작합시다. `-output-dir` CLI 옵션을 사용합니다.

#### 2.1.1. `-output-dir`로 파이프라인 실행

`-output-dir` 옵션(단축형: `-o`)은 모든 워크플로우 출력에 대한 기본 출력 디렉토리(`results/`)를 재정의합니다.
이것은 출력이 게시되는 루트 경로를 제어하는 권장 방법입니다.

```bash
nextflow run hello-config.nf -output-dir custom-outdir-cli/
```

??? success "명령 출력"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `hello-config.nf` [prickly_kay] DSL2 - revision: 32ecc4fba2

    executor >  local (8)
    [9f/332636] sayHello (1)       [100%] 3 of 3 ✔
    [03/a55991] convertToUpper (3) [100%] 3 of 3 ✔
    [e5/ab7893] collectGreetings   [100%] 1 of 1 ✔
    [a8/97338e] cowpy              [100%] 1 of 1 ✔
    ```

이것은 `results/` 대신 `custom-outdir-cli/`에 출력을 게시합니다:

??? abstract "디렉토리 내용"

    ```console
    custom-outdir-cli/
    └── hello_config
        ├── batch-report.txt
        ├── cowpy-COLLECTED-batch-output.txt
        └── intermediates
            ├── Bonjour-output.txt
            ├── COLLECTED-batch-output.txt
            ├── Hello-output.txt
            ├── Holà-output.txt
            ├── UPPER-Bonjour-output.txt
            ├── UPPER-Hello-output.txt
            └── UPPER-Holà-output.txt
    ```

output 블록의 `path` 선언에서 `hello_config` 하위 디렉토리가 여전히 있음을 주목하십시오.
이를 정리합시다.

#### 2.1.2. output 블록에서 하드코딩된 경로 제거

`hello_config/` 접두사는 이전 장에서 하드코딩되었지만, 이제 출력 경로를 유연하게 구성하는 방법을 배우고 있으므로 이 하드코딩을 제거할 수 있습니다.
하위 디렉토리가 필요하지 않은 출력의 경우 `path` 지시문을 빈 문자열로 설정하거나 완전히 제거할 수 있습니다.

워크플로우 파일에서 다음 코드 변경을 수행하십시오:

=== "수정 후"

    ```groovy title="hello-config.nf" linenums="42" hl_lines="3 7 11 15 19"
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

=== "수정 전"

    ```groovy title="hello-config.nf" linenums="42" hl_lines="3 7 11 15 19"
    output {
        first_output {
            path 'hello_config/intermediates'
            mode 'copy'
        }
        uppercased {
            path 'hello_config/intermediates'
            mode 'copy'
        }
        collected {
            path 'hello_config/intermediates'
            mode 'copy'
        }
        batch_report {
            path 'hello_config'
            mode 'copy'
        }
        cowpy_art {
            path 'hello_config'
            mode 'copy'
        }
    }
    ```

파이프라인을 다시 실행하십시오:

```bash
nextflow run hello-config.nf -output-dir custom-outdir-cli-2/
```

이제 출력이 `hello_config` 하위 디렉토리 없이 `custom-outdir-cli-2/` 바로 아래에 게시됩니다:

??? abstract "디렉토리 내용"

    ```console
    custom-outdir-cli-2/
    ├── batch-report.txt
    ├── cowpy-COLLECTED-batch-output.txt
    └── intermediates
        ├── Bonjour-output.txt
        ├── COLLECTED-batch-output.txt
        ├── Hello-output.txt
        ├── Holà-output.txt
        ├── UPPER-Bonjour-output.txt
        ├── UPPER-Hello-output.txt
        └── UPPER-Holà-output.txt
    ```

!!! tip "팁"

    `-output-dir` 옵션은 출력이 _어디로_ 가는지 제어하는 데 사용되며, output 블록의 `path` 지시문은 _하위 디렉토리 구조_를 제어합니다.

### 2.2. 동적 출력 경로

CLI를 통해 출력 디렉토리를 변경하는 것 외에도, `outputDir`을 사용하여 구성 파일에서 사용자 정의 기본값을 설정할 수 있습니다.
이를 통해 디렉토리 경로를 동적으로 설정할 수 있습니다. 즉, 정적 문자열만 사용하는 것이 아닙니다.

#### 2.2.1. 구성 파일에서 `outputDir` 값 설정

`nextflow.config` 파일에 다음 코드를 추가하십시오:

=== "수정 후"

    ```groovy title="nextflow.config" linenums="9" hl_lines="10-13"
    /*
    * 파이프라인 매개변수
    */
    params {
        input = 'data/greetings.csv'
        batch = 'batch'
        character = 'turkey'
    }

    /*
    * 출력 설정
    */
    outputDir = "custom-outdir-config/${params.batch}"
    ```

=== "수정 전"

    ```groovy title="nextflow.config" linenums="9"
    /*
    * 파이프라인 매개변수
    */
    params {
        input = 'data/greetings.csv'
        batch = 'batch'
        character = 'turkey'
    }
    ```

이렇게 하면 출력 디렉토리가 `custom-outdir-config/` 더하기 하위 디렉토리로서 `batch` 매개변수 값으로 설정됩니다.
이제 `--batch` 매개변수를 설정하여 출력 위치를 변경할 수 있습니다:

```bash
nextflow run hello-config.nf --batch my_run
```

이것은 `custom-outdir-config/my_run/`에 출력을 게시합니다.

!!! note "참고"

    `-output-dir` CLI 옵션이 `outputDir` 구성 설정보다 우선합니다.
    설정되면 구성 옵션은 완전히 무시됩니다.

#### 2.2.2. batch 및 프로세스 이름을 사용한 하위 디렉토리

출력별로 하위 디렉토리 출력 `path` 선언을 동적으로 설정할 수도 있습니다.

예를 들어, 출력 경로 선언에서 `<process>.name`을 참조하여 프로세스별로 출력을 구성할 수 있습니다:

=== "수정 후"

    ```groovy title="hello-config.nf" linenums="42" hl_lines="3 7 11 15 19"
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

=== "수정 전"

    ```groovy title="hello-config.nf" linenums="42" hl_lines="3 7 11 15 19"
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

더 나아가 더 복잡한 하위 디렉토리 경로를 구성할 수 있습니다.

위 편집에서 `intermediates`와 최상위 수준의 최종 출력 간의 구분이 지워졌습니다.
이를 다시 추가하고 파일을 `params.batch` 하위 디렉토리에도 넣읍시다.

!!! tip "팁"

    output 블록 `path`에 `params.batch`를 포함하면, `outputDir` 구성 대신에 포함하므로 CLI에서 `-output-dir`로 재정의되지 않습니다.

먼저 `outputDir`에서 `${params.batch}`를 제거하도록 구성 파일을 업데이트하십시오(path 선언으로 이동하고 있으므로):

=== "수정 후"

    ```groovy title="nextflow.config" linenums="12" hl_lines="4"
    /*
    * 출력 설정
    */
    outputDir = "custom-outdir-config/"
    ```

=== "수정 전"

    ```groovy title="nextflow.config" linenums="12" hl_lines="4"
    /*
    * 출력 설정
    */
    outputDir = "custom-outdir-config/${params.batch}"
    ```

그런 다음 워크플로우 파일에서 다음 변경을 수행하십시오:

=== "수정 후"

    ```groovy title="hello-config.nf" linenums="42" hl_lines="3 7 11 15 19"
    output {
        first_output {
            path { "${params.batch}/intermediates/${sayHello.name}" }
            mode 'copy'
        }
        uppercased {
            path { "${params.batch}/intermediates/${convertToUpper.name}" }
            mode 'copy'
        }
        collected {
            path { "${params.batch}/intermediates/${collectGreetings.name}" }
            mode 'copy'
        }
        batch_report {
            path { "${params.batch}/${collectGreetings.name}" }
            mode 'copy'
        }
        cowpy_art {
            path { "${params.batch}/${cowpy.name}" }
            mode 'copy'
        }
    }
    ```

=== "수정 전"

    ```groovy title="hello-config.nf" linenums="42" hl_lines="3 7 11 15 19"
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

#### 2.2.3. 파이프라인 실행

실제로 어떻게 작동하는지 봅시다. 명령줄에서 `-output-dir`(또는 단축형 `-o`)를 `custom-outdir-config-2`로 설정하고 batch 이름을 `rep2`로 설정합니다:

```bash
nextflow run hello-config.nf -output-dir custom-outdir-config-2 --batch rep2
```

??? success "명령 출력"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `hello-config.nf` [mad_curry] DSL2 - revision: 668a98ccb9

    executor >  local (8)
    [9e/6095e0] sayHello (1)       [100%] 3 of 3 ✔
    [05/454d52] convertToUpper (3) [100%] 3 of 3 ✔
    [ed/e3ddfb] collectGreetings   [100%] 1 of 1 ✔
    [39/5e063a] cowpy              [100%] 1 of 1 ✔
    ```

이것은 지정된 기본 경로 _그리고_ batch 이름 하위 디렉토리 _그리고_ 프로세스별로 그룹화된 결과와 함께 `custom-outdir-config-2/rep2/`에 출력을 게시합니다:

??? abstract "디렉토리 내용"

    ```console
    custom-outdir-config-2
    └── rep2
        ├── collectGreetings
        │   └── rep2-report.txt
        ├── cowpy
        │   └── cowpy-COLLECTED-rep2-output.txt
        └── intermediates
            ├── collectGreetings
            │   └── COLLECTED-rep2-output.txt
            ├── convertToUpper
            │   ├── UPPER-Bonjour-output.txt
            │   ├── UPPER-Hello-output.txt
            │   └── UPPER-Holà-output.txt
            └── sayHello
                ├── Bonjour-output.txt
                ├── Hello-output.txt
                └── Holà-output.txt
    ```

### 2.3. 워크플로우 수준에서 게시 모드 설정

마지막으로, 반복적인 코드의 양을 줄이는 정신으로 출력별 `mode` 선언을 구성의 단일 줄로 대체할 수 있습니다.

#### 2.3.1. 구성 파일에 `workflow.output.mode` 추가

`nextflow.config` 파일에 다음 코드를 추가하십시오:

=== "수정 후"

    ```groovy title="nextflow.config" linenums="12" hl_lines="5"
    /*
    * 출력 설정
    */
    outputDir = "custom-outdir-config/"
    workflow.output.mode = 'copy'
    ```

=== "수정 전"

    ```groovy title="nextflow.config" linenums="12"
    /*
    * 출력 설정
    */
    outputDir = "custom-outdir-config/"
    ```

구성 파일에서 `workflow.output.mode`에 값을 지정하면 워크플로우 파일에 설정된 것을 재정의하기에 충분하지만, 불필요한 코드를 어쨌든 제거합시다.

#### 2.3.2. 워크플로우 파일에서 출력 모드 제거

워크플로우 파일에서 다음 변경을 수행하십시오:

=== "수정 후"

    ```groovy title="hello-config.nf" linenums="42"
    output {
        first_output {
            path { "${params.batch}/intermediates/${sayHello.name}" }
        }
        uppercased {
            path { "${params.batch}/intermediates/${convertToUpper.name}" }
        }
        collected {
            path { "${params.batch}/intermediates/${collectGreetings.name}" }
        }
        batch_report {
            path { "${params.batch}/${collectGreetings.name}" }
        }
        cowpy_art {
            path { "${params.batch}/${cowpy.name}" }
        }
    }
    ```

=== "수정 전"

    ```groovy title="hello-config.nf" linenums="42" hl_lines="4 8 12 16 20"
    output {
        first_output {
            path { "${params.batch}/intermediates/${sayHello.name}" }
            mode 'copy'
        }
        uppercased {
            path { "${params.batch}/intermediates/${convertToUpper.name}" }
            mode 'copy'
        }
        collected {
            path { "${params.batch}/intermediates/${collectGreetings.name}" }
            mode 'copy'
        }
        batch_report {
            path { "${params.batch}/${collectGreetings.name}" }
            mode 'copy'
        }
        cowpy_art {
            path { "${params.batch}/${cowpy.name}" }
            mode 'copy'
        }
    }
    ```

더 간결하지 않습니까?

#### 2.3.3. 파이프라인 실행

올바르게 작동하는지 테스트해 봅시다:

```bash
nextflow run hello-config.nf -output-dir config-output-mode
```

??? success "명령 출력"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `hello-config.nf` [small_stone] DSL2 - revision: 024d6361b5

    executor >  local (8)
    [e8/a0e93e] sayHello (1)       [100%] 3 of 3 ✔
    [14/176c9d] convertToUpper (3) [100%] 3 of 3 ✔
    [23/d667ca] collectGreetings   [100%] 1 of 1 ✔
    [e6/1dc80e] cowpy              [100%] 1 of 1 ✔
    ```

이것은 `config-output-mode/`에 출력을 게시하며, 여전히 심볼릭 링크가 아닌 적절한 복사본입니다.

??? abstract "디렉토리 내용"

    ```console
    config-output-mode
    └── batch
        ├── collectGreetings
        │   └── batch-report.txt
        ├── cowpy
        │   └── cowpy-COLLECTED-batch-output.txt
        └── intermediates
            ├── collectGreetings
            │   └── COLLECTED-batch-output.txt
            ├── convertToUpper
            │   ├── UPPER-Bonjour-output.txt
            │   ├── UPPER-Hello-output.txt
            │   └── UPPER-Holà-output.txt
            └── sayHello
                ├── Bonjour-output.txt
                ├── Hello-output.txt
                └── Holà-output.txt
    ```

출력별 모드 설정 방식을 여전히 사용하려는 주된 이유는 동일한 워크플로우 내에서 혼합하고 싶은 경우입니다. 즉, 일부 출력은 복사하고 일부는 심볼릭 링크로 설정합니다.

이 방식으로 사용자 정의할 수 있는 다른 옵션도 많이 있지만, 이를 통해 옵션의 범위와 선호도에 맞게 효과적으로 활용하는 방법에 대한 감각을 얻을 수 있기를 바랍니다.

### 핵심 정리

출력이 게시되는 디렉토리의 이름 지정 및 구조와 워크플로우 출력 게시 모드를 제어하는 방법을 알게 되었습니다.

### 다음 단계

소프트웨어 패키징 기술부터 시작하여 컴퓨팅 환경에 워크플로우 구성을 적응시키는 방법을 배웁니다.

---

## 3. 소프트웨어 패키징 기술 선택

지금까지 입력이 어떻게 들어가고 출력이 어디에서 나오는지를 제어하는 구성 요소를 살펴보았습니다. 이제 컴퓨팅 환경에 워크플로우 구성을 적응시키는 데 더 구체적으로 집중할 때입니다.

그 경로의 첫 번째 단계는 각 단계에서 실행될 소프트웨어 패키지가 어디에서 오는지 지정하는 것입니다.
로컬 컴퓨팅 환경에 이미 설치되어 있습니까?
이미지를 검색하고 컨테이너 시스템을 통해 실행해야 합니까?
아니면 Conda 패키지를 검색하고 로컬 Conda 환경을 빌드해야 합니까?

이 교육 과정의 첫 번째 부분(파트 1-4)에서는 워크플로우에서 로컬로 설치된 소프트웨어만 사용했습니다.
그런 다음 파트 5에서 Docker 컨테이너와 `nextflow.config` 파일을 도입하여 Docker 컨테이너 사용을 활성화했습니다.

이제 `nextflow.config` 파일을 통해 대체 소프트웨어 패키징 옵션을 구성하는 방법을 살펴봅시다.

### 3.1. 구성 파일에서 Docker 비활성화 및 Conda 활성화

보안상의 이유로 관리자가 Docker 사용을 허용하지 않는 HPC 클러스터에서 작업한다고 가정해 봅시다.
다행히도 Nextflow는 Singularity(HPC에서 더 널리 사용됨)를 포함한 여러 다른 컨테이너 기술과 Conda와 같은 소프트웨어 패키지 관리자를 지원합니다.

Docker 대신 [Conda](https://nextflow.io/docs/latest/conda.html)를 사용하도록 구성 파일을 변경할 수 있습니다.
그렇게 하려면 `docker.enabled` 값을 `false`로 전환하고 Conda 사용을 활성화하는 지시문을 추가합시다:

=== "수정 후"

    ```groovy title="nextflow.config" linenums="1" hl_lines="1-2"
    docker.enabled = false
    conda.enabled = true
    ```

=== "수정 전"

    ```groovy title="nextflow.config" linenums="1" hl_lines="1"
    docker.enabled = true
    ```

이렇게 하면 Nextflow가 Conda 패키지가 지정된 프로세스에 대해 Conda 환경을 생성하고 활용할 수 있습니다.
즉, 이제 `cowpy` 프로세스에 그 중 하나를 추가해야 합니다!

### 3.2. 프로세스 정의에서 Conda 패키지 지정

`cowpy` 도구가 포함된 Conda 패키지의 URI를 이미 검색했습니다: `conda-forge::cowpy==1.1.5`

이제 `conda` 지시문을 사용하여 `cowpy` 프로세스 정의에 URI를 추가합니다:

=== "수정 후"

    ```groovy title="modules/cowpy.nf" linenums="4" hl_lines="4"
    process cowpy {

        container 'community.wave.seqera.io/library/cowpy:1.1.5--3db457ae1977a273'
        conda 'conda-forge::cowpy==1.1.5'

        input:
    ```

=== "수정 전"

    ```groovy title="modules/cowpy.nf" linenums="4"
    process cowpy {

        container 'community.wave.seqera.io/library/cowpy:1.1.5--3db457ae1977a273'

        input:
    ```

명확히 하자면, `docker` 지시문을 *교체*하는 것이 아니라 대체 옵션을 *추가*하는 것입니다.

!!! tip "팁"

    주어진 conda 패키지에 대한 URI를 얻는 몇 가지 다른 방법이 있습니다.
    컨테이너를 생성할 계획이 없더라도 복사하여 붙여넣을 수 있는 URI를 제공하는 [Seqera Containers](https://seqera.io/containers/) 검색 쿼리를 사용하는 것이 좋습니다.

### 3.3. 워크플로우를 실행하여 Conda를 사용할 수 있는지 확인

시도해 봅시다.

```bash
nextflow run hello-config.nf --batch conda
```

??? success "명령 출력"

    ```console title="Output"
    N E X T F L O W   ~  version 25.10.2

    Launching `hello-config.nf` [friendly_lamport] DSL2 - revision: 024d6361b5

    executor >  local (8)
    [e8/91c116] sayHello (2)       [100%] 3 of 3 ✔
    [fe/6a70ce] convertToUpper (3) [100%] 3 of 3 ✔
    [99/7cc493] collectGreetings   [100%] 1 of 1 ✔
    [3c/09fb59] cowpy              [100%] 1 of 1 ✔
    ```

이것은 문제 없이 작동하고 `custom-outdir-config/conda` 아래에 이전과 동일한 출력을 생성해야 합니다.

백그라운드에서 Nextflow는 Conda 패키지를 검색하고 환경을 생성했으며, 이는 일반적으로 약간의 작업이 필요합니다. 그래서 우리가 직접 그 작업을 하지 않아도 되니 좋습니다!

!!! note "참고"

    `cowpy` 패키지가 상당히 작기 때문에 빠르게 실행되지만, 큰 패키지로 작업하는 경우 처음에는 평소보다 약간 더 오래 걸릴 수 있으며, 콘솔 출력이 완료되기 전에 1분 정도 '멈춰' 있는 것처럼 보일 수 있습니다.
    이것은 정상이며 새 패키지를 처음 사용할 때 Nextflow가 수행하는 추가 작업 때문입니다.

우리 관점에서 보면 백엔드의 메커니즘이 약간 다르더라도 Docker로 실행하는 것과 정확히 동일하게 작동하는 것 같습니다.

이것은 필요한 경우 Conda 환경으로 실행할 준비가 되었음을 의미합니다.

??? info "Docker와 Conda 혼합 사용"

    이러한 지시문은 프로세스별로 할당되므로 '혼합 사용'이 가능합니다. 즉, 사용 중인 컴퓨팅 인프라가 둘 다 지원하는 경우 예를 들어 워크플로우의 일부 프로세스는 Docker로 실행하고 다른 프로세스는 Conda로 실행하도록 구성할 수 있습니다.
    이 경우 구성 파일에서 Docker와 Conda를 모두 활성화합니다.
    주어진 프로세스에 대해 둘 다 사용 가능한 경우 Nextflow는 컨테이너를 우선시합니다.

    앞서 언급했듯이 Nextflow는 여러 다른 소프트웨어 패키징 및 컨테이너 기술을 지원하므로 이 두 가지에만 국한되지 않습니다.

### 핵심 정리

각 프로세스가 사용해야 하는 소프트웨어 패키지를 구성하는 방법과 기술 간 전환 방법을 알게 되었습니다.

### 다음 단계

Nextflow가 실제로 작업을 수행하는 데 사용하는 실행 플랫폼을 변경하는 방법을 배웁니다.

---

## 4. 실행 플랫폼 선택

지금까지 로컬 실행자로 파이프라인을 실행해 왔습니다.
이것은 Nextflow가 실행 중인 머신에서 각 작업을 실행합니다.
Nextflow가 시작되면 사용 가능한 CPU와 메모리를 확인합니다.
실행할 준비가 된 작업의 리소스가 사용 가능한 리소스를 초과하면 Nextflow는 하나 이상의 이전 작업이 완료되어 필요한 리소스가 확보될 때까지 마지막 작업을 실행에서 보류합니다.

로컬 실행자는 편리하고 효율적이지만 해당 단일 머신으로 제한됩니다. 매우 큰 워크로드의 경우 사용 가능한 것보다 더 많은 리소스가 필요한 단일 작업이 있거나 단일 머신에서 실행하기를 기다리는 데 너무 오래 걸릴 정도로 많은 작업이 있어 로컬 머신이 병목 현상을 일으키는 것을 발견할 수 있습니다.

Nextflow는 HPC 스케줄러(Slurm, LSF, SGE, PBS, Moab, OAR, Bridge, HTCondor 등)뿐만 아니라 클라우드 실행 백엔드(AWS Batch, Google Cloud Batch, Azure Batch, Kubernetes 등)를 포함한 [많은 다른 실행 백엔드](https://nextflow.io/docs/latest/executor.html)를 지원합니다.

### 4.1. 다른 백엔드 대상 지정

실행자 선택은 `executor`라는 프로세스 지시문으로 설정됩니다.
기본적으로 `local`로 설정되어 있으므로 다음 구성이 암시됩니다:

```groovy title="Built-in configuration"
process {
    executor = 'local'
}
```

다른 백엔드를 대상으로 하도록 실행자를 설정하려면 리소스 할당에 대해 위에서 설명한 것과 유사한 구문을 사용하여 원하는 실행자를 지정하기만 하면 됩니다(모든 옵션에 대해서는 [실행자 문서](https://nextflow.io/docs/latest/executor.html) 참조).

```groovy title="nextflow.config"
process {
    executor = 'slurm'
}
```

!!! warning "경고"

    HPC에 연결하도록 설정되지 않았기 때문에 교육 환경에서 이것을 실제로 테스트할 수 없습니다.

### 4.2. 실행 매개변수에 대한 백엔드별 구문 처리

대부분의 고성능 컴퓨팅 플랫폼은 리소스 할당 요청 및 제한(예: CPU 수 및 메모리) 및 사용할 작업 큐의 이름과 같은 특정 매개변수를 지정할 수 있도록 허용하고 때로는 요구합니다.

불행히도 이러한 각 시스템은 작업이 정의되고 관련 스케줄러에 제출되는 방법을 정의하는 데 다른 기술, 구문 및 구성을 사용합니다.

??? abstract "예제"

    예를 들어, "my-science-work" 큐에서 8개의 CPU와 4GB의 RAM이 필요한 동일한 작업은 백엔드에 따라 다른 방식으로 표현해야 합니다.

    ```bash title="Config for SLURM / submit using sbatch"
    #SBATCH -o /path/to/my/task/directory/my-task-1.log
    #SBATCH --no-requeue
    #SBATCH -c 8
    #SBATCH --mem 4096M
    #SBATCH -p my-science-work
    ```

    ```bash title="Config for PBS / submit using qsub"
    #PBS -o /path/to/my/task/directory/my-task-1.log
    #PBS -j oe
    #PBS -q my-science-work
    #PBS -l nodes=1:ppn=5
    #PBS -l mem=4gb
    ```

    ```bash title="Config for SGE / submit using qsub"
    #$ -o /path/to/my/task/directory/my-task-1.log
    #$ -j y
    #$ -terse
    #$ -notify
    #$ -q my-science-work
    #$ -l slots=5
    #$ -l h_rss=4096M,mem_free=4096M
    ```

다행히도 Nextflow는 이 모든 것을 단순화합니다.
[`cpus`](https://nextflow.io/docs/latest/reference/process.html#cpus), [`memory`](https://nextflow.io/docs/latest/reference/process.html#memory) 및 [`queue`](https://nextflow.io/docs/latest/reference/process.html#queue)(다른 속성은 [프로세스 지시문](https://nextflow.io/docs/latest/reference/process.html#process-directives) 참조)와 같은 관련 속성을 한 번만 지정할 수 있도록 표준화된 구문을 제공합니다.
그런 다음 런타임에 Nextflow는 해당 설정을 사용하여 실행자 설정에 따라 적절한 백엔드별 스크립트를 생성합니다.

다음 섹션에서 해당 표준화된 구문을 다룰 것입니다.

### 핵심 정리

이제 다른 종류의 컴퓨팅 인프라를 사용하도록 실행자를 변경하는 방법을 알게 되었습니다.

### 다음 단계

Nextflow에서 리소스 할당 및 제한을 평가하고 표현하는 방법을 배웁니다.

---

## 5. 컴퓨팅 리소스 할당 제어

대부분의 고성능 컴퓨팅 플랫폼은 CPU 수 및 메모리와 같은 특정 리소스 할당 매개변수를 지정할 수 있도록 허용하고 때로는 요구합니다.

기본적으로 Nextflow는 각 프로세스에 대해 단일 CPU와 2GB의 메모리를 사용합니다.
해당 프로세스 지시문은 `cpus` 및 `memory`라고 하므로 다음 구성이 암시됩니다:

```groovy title="Built-in configuration" linenums="1"
process {
    cpus = 1
    memory = 2.GB
}
```

구성 파일에서 추가 프로세스 지시문을 사용하여 모든 프로세스 또는 특정 명명된 프로세스에 대해 이러한 값을 수정할 수 있습니다.
Nextflow는 선택한 실행자에 적합한 명령으로 변환합니다.

그러나 어떤 값을 사용해야 하는지 어떻게 알 수 있습니까?

### 5.1. 리소스 활용 보고서를 생성하기 위해 워크플로우 실행

프로세스에 얼마나 많은 CPU와 메모리가 필요한지 미리 알지 못하는 경우 리소스 프로파일링을 수행할 수 있습니다. 즉, 일부 기본 할당으로 워크플로우를 실행하고, 각 프로세스가 사용한 양을 기록하고, 거기서 기본 할당을 조정하는 방법을 추정합니다.

편리하게도 Nextflow에는 이를 수행하기 위한 내장 도구가 포함되어 있으며 요청 시 보고서를 기꺼이 생성합니다.

그렇게 하려면 명령줄에 `-with-report <filename>.html`을 추가하십시오.

```bash
nextflow run hello-config.nf -with-report report-config-1.html
```

보고서는 html 파일이며 다운로드하여 브라우저에서 열 수 있습니다. 왼쪽의 파일 탐색기에서 마우스 오른쪽 버튼으로 클릭하고 `미리 보기 표시`를 클릭하여 교육 환경에서 볼 수도 있습니다.

보고서를 살펴보고 리소스를 조정할 수 있는 기회를 식별할 수 있는지 확인하는 데 몇 분을 투자하십시오.
할당된 것에 대한 백분율로 활용 결과를 보여주는 탭을 클릭해야 합니다.

사용 가능한 모든 기능을 설명하는 [보고서 문서](https://nextflow.io/docs/latest/reports.html)가 있습니다.

### 5.2. 모든 프로세스에 대한 리소스 할당 설정

프로파일링은 교육 워크플로우의 프로세스가 매우 가볍다는 것을 보여주므로 프로세스당 기본 메모리 할당을 1GB로 줄입시다.

파이프라인 매개변수 섹션 앞에 `nextflow.config` 파일에 다음을 추가하십시오:

=== "수정 후"

    ```groovy title="nextflow.config" linenums="1" hl_lines="4-9"
    docker.enabled = false
    conda.enabled = true

    /*
    * 프로세스 설정
    */
    process {
        memory = 1.GB
    }

    /*
    * 파이프라인 매개변수
    */
    params {
        input = 'data/greetings.csv'
        batch = 'batch'
        character = 'turkey'
    }
    ```

=== "수정 전"

    ```groovy title="nextflow.config" linenums="1"
    docker.enabled = false
    conda.enabled = true

    /*
    * 파이프라인 매개변수
    */
    params {
        input = 'data/greetings.csv'
        batch = 'batch'
        character = 'turkey'
    }
    ```

이렇게 하면 소비하는 컴퓨팅 양을 줄이는 데 도움이 됩니다.

### 5.3. 특정 프로세스에 대한 리소스 할당 설정

동시에 `cowpy` 프로세스가 다른 프로세스보다 더 많은 리소스를 필요로 한다고 가정하여 개별 프로세스에 대한 할당을 조정하는 방법을 시연하겠습니다.

=== "수정 후"

    ```groovy title="nextflow.config" linenums="4" hl_lines="6-9"
    /*
    * 프로세스 설정
    */
    process {
        memory = 1.GB
        withName: 'cowpy' {
            memory = 2.GB
            cpus = 2
        }
    }
    ```

=== "수정 전"

    ```groovy title="nextflow.config" linenums="4"
    /*
    * 프로세스 설정
    */
    process {
        memory = 1.GB
    }
    ```

이 구성으로 모든 프로세스는 1GB의 메모리와 단일 CPU(암시된 기본값)를 요청하지만, `cowpy` 프로세스는 2GB와 2개의 CPU를 요청합니다.

!!! tip "팁"

    CPU가 적은 머신이 있고 프로세스당 많은 수를 할당하면 프로세스 호출이 서로 뒤에 대기열에 있는 것을 볼 수 있습니다.
    이것은 Nextflow가 사용 가능한 것보다 더 많은 CPU를 요청하지 않도록 보장하기 때문입니다.

### 5.4. 업데이트된 구성으로 워크플로우 실행

시도해 봅시다. 구성 변경 전후의 성능을 비교할 수 있도록 프로파일링 보고서에 대해 다른 파일 이름을 제공합니다.

```bash
nextflow run hello-config.nf -with-report report-config-2.html
```

워크로드가 너무 작기 때문에 실제 차이를 느끼지 못할 수 있지만, 이것은 실제 워크플로우의 성능 및 리소스 요구 사항을 분석하는 데 사용할 접근 방식입니다.

프로세스에 다른 리소스 요구 사항이 있는 경우 매우 유용합니다. 추측이 아닌 실제 데이터를 기반으로 각 프로세스에 대해 설정하는 리소스 할당을 적절한 크기로 조정할 수 있습니다.

!!! tip "팁"

    이것은 리소스 사용을 최적화하기 위해 할 수 있는 것의 아주 작은 맛보기에 불과합니다.
    Nextflow 자체에는 리소스 제한으로 인해 실패한 작업을 재시도하는 정말 멋진 [동적 재시도 로직](https://nextflow.io/docs/latest/process.html#dynamic-task-resources)이 내장되어 있습니다.
    또한 Seqera Platform은 리소스 할당을 자동으로 최적화하기 위한 AI 기반 도구도 제공합니다.

### 5.5. 리소스 제한 추가

사용 중인 컴퓨팅 실행자와 컴퓨팅 인프라에 따라 할당할 수 있는(또는 해야 하는) 것에 대한 제약이 있을 수 있습니다.
예를 들어, 클러스터에서 특정 제한 내에 있어야 할 수 있습니다.

`resourceLimits` 지시문을 사용하여 관련 제한을 설정할 수 있습니다. process 블록에 단독으로 있을 때 구문은 다음과 같습니다:

```groovy title="Syntax example"
process {
    resourceLimits = [
        memory: 750.GB,
        cpus: 200,
        time: 30.d
    ]
}
```

Nextflow는 지정한 실행자에 따라 이러한 값을 적절한 명령으로 변환합니다.

교육 환경에서 관련 인프라에 액세스할 수 없으므로 이것을 실행하지 않을 것입니다.
그러나 이러한 제한을 초과하는 리소스 할당으로 워크플로우를 실행한 다음 `.command.run` 스크립트 파일에서 `sbatch` 명령을 조회하면 실행자에게 실제로 전송되는 요청이 `resourceLimits`에서 지정한 값으로 제한되어 있는 것을 볼 수 있습니다.

??? info "기관 참조 구성"

    nf-core 프로젝트는 전 세계 다양한 기관에서 공유하는 [구성 파일 모음](https://nf-co.re/configs/)을 컴파일하여 광범위한 HPC 및 클라우드 실행자를 다룹니다.

    이러한 공유 구성은 해당 기관에서 일하는 사람들이 자신의 기관 구성을 바로 활용할 수 있다는 점과 자체 인프라에 대한 구성을 개발하려는 사람들을 위한 모델로서 모두 가치가 있습니다.

### 핵심 정리

리소스 활용을 평가하기 위한 프로파일링 보고서를 생성하는 방법과 모든 프로세스 및/또는 개별 프로세스에 대한 리소스 할당을 수정하는 방법, HPC에서 실행하기 위한 리소스 제한을 설정하는 방법을 알게 되었습니다.

### 다음 단계

사전 설정 구성 프로필을 설정하고 런타임에 전환하는 방법을 배웁니다.

---

## 6. 프로필을 사용하여 사전 설정 구성 간 전환

작업 중인 프로젝트나 사용 중인 컴퓨팅 환경에 따라 파이프라인 구성을 사용자 정의할 수 있는 여러 방법을 보여드렸습니다.

사용 중인 컴퓨팅 인프라에 따라 대체 설정 간에 전환하고 싶을 수 있습니다. 예를 들어, 노트북에서 로컬로 개발하고 소규모 테스트를 실행한 다음 HPC 또는 클라우드에서 전체 규모 워크로드를 실행하고 싶을 수 있습니다.

Nextflow를 사용하면 다양한 구성을 설명하는 여러 [프로필](https://nextflow.io/docs/latest/config.html#config-profiles)을 설정할 수 있으며, 구성 파일 자체를 수정하지 않고 명령줄 인수를 사용하여 런타임에 선택할 수 있습니다.

### 6.1. 로컬 개발과 HPC에서 실행 간 전환을 위한 프로필 생성

두 가지 대체 프로필을 설정합시다. 하나는 Docker 컨테이너를 사용하는 일반 컴퓨터에서 소규모 로드를 실행하고, 하나는 Conda 패키지를 사용하는 Slurm 스케줄러가 있는 대학 HPC에서 실행하기 위한 것입니다.

#### 6.1.1. 프로필 설정

파이프라인 매개변수 섹션 뒤에 그리고 출력 설정 앞에 `nextflow.config` 파일에 다음을 추가하십시오:

=== "수정 후"

    ```groovy title="nextflow.config" linenums="15" hl_lines="10-27"
    /*
    * 파이프라인 매개변수
    */
    params {
        input = 'data/greetings.csv'
        batch = 'batch'
        character = 'turkey'
    }

    /*
    * 프로필
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

    /*
    * 출력 설정
    */
    outputDir = "custom-outdir-config/"
    workflow.output.mode = 'copy'
    ```

=== "수정 전"

    ```groovy title="nextflow.config" linenums="15"
    /*
    * 파이프라인 매개변수
    */
    params {
        input = 'data/greetings.csv'
        batch = 'batch'
        character = 'turkey'
    }

    /*
    * 출력 설정
    */
    outputDir = "custom-outdir-config/"
    workflow.output.mode = 'copy'
    ```

대학 HPC의 경우 리소스 제한도 지정하고 있습니다.

#### 6.1.2. 프로필로 워크플로우 실행

Nextflow 명령줄에서 프로필을 지정하려면 `-profile` 인수를 사용합니다.

`my_laptop` 구성으로 워크플로우를 실행해 봅시다.

```bash
nextflow run hello-config.nf -profile my_laptop
```

??? success "명령 출력"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `hello-config.nf` [hungry_sanger] DSL2 - revision: 024d6361b5

    executor >  local (8)
    [b0/fb2ec9] sayHello (3)       [100%] 3 of 3 ✔
    [4a/e039f0] convertToUpper (3) [100%] 3 of 3 ✔
    [6f/408fa9] collectGreetings   [100%] 1 of 1 ✔
    [f1/fd6520] cowpy              [100%] 1 of 1 ✔
    ```

보시다시피 런타임에 매우 편리하게 구성 간에 전환할 수 있습니다.

!!! warning "경고"

    Slurm 스케줄러에 액세스할 수 없으므로 교육 환경에서 `univ_hpc` 프로필은 제대로 실행되지 않습니다.

앞으로 이러한 프로필과 항상 함께 발생하는 다른 구성 요소를 찾으면 해당 프로필에 간단히 추가할 수 있습니다.
함께 그룹화하려는 다른 구성 요소가 있는 경우 추가 프로필을 만들 수도 있습니다.

### 6.2. 테스트 매개변수 프로필 생성

프로필은 인프라 구성만을 위한 것이 아닙니다.
다른 사람들이 적절한 입력 값을 직접 수집하지 않고도 워크플로우를 쉽게 시험해 볼 수 있도록 워크플로우 매개변수의 기본값을 설정하는 데도 사용할 수 있습니다.
이것을 매개변수 파일 사용의 대안으로 고려할 수 있습니다.

#### 6.2.1. 프로필 설정

이 컨텍스트에서 기본값을 표현하는 구문은 `test`라고 이름 지은 프로필에 대해 다음과 같습니다:

```groovy title="Syntax example"
    test {
        params.<parameter1>
        params.<parameter2>
        ...
    }
```

워크플로우에 대한 테스트 프로필을 추가하면 `profiles` 블록은 다음과 같이 됩니다:

```groovy title="nextflow.config" linenums="24" hl_lines="18-22"
/*
* 프로필
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

#### 6.2.2. 테스트 프로필로 로컬에서 워크플로우 실행

편리하게도 프로필은 상호 배타적이지 않으므로 다음 구문 `-profile <profile1>,<profile2>`를 사용하여 명령줄에서 여러 프로필을 지정할 수 있습니다(프로필 수에 관계없이).

동일한 구성 파일에 설명되어 있고 동일한 구성 요소에 대해 값을 설정하는 프로필을 결합하면 Nextflow는 마지막으로 읽은 값을 사용하여 충돌을 해결합니다(즉, 파일에서 나중에 나오는 것).
충돌하는 설정이 다른 구성 소스에 설정된 경우 기본 [우선순위](https://nextflow.io/docs/latest/config.html)가 적용됩니다.

이전 명령에 테스트 프로필을 추가해 봅시다:

```bash
nextflow run hello-config.nf -profile my_laptop,test
```

??? success "명령 출력"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `hello-config.nf` [modest_becquerel] DSL2 - revision: 024d6361b5

    executor >  local (8)
    [4c/fe2580] sayHello (1)       [100%] 3 of 3 ✔
    [fd/7d9017] convertToUpper (3) [100%] 3 of 3 ✔
    [13/1523bd] collectGreetings   [100%] 1 of 1 ✔
    [06/a1ee14] cowpy              [100%] 1 of 1 ✔
    ```

이것은 가능한 경우 Docker를 사용하고 `custom-outdir-config/test` 아래에 출력을 생성하며, 이번에 캐릭터는 코믹 듀오 `dragonandcow`입니다.

??? abstract "파일 내용"

    ```console title="custom-outdir-config/test/cowpy/cowpy-COLLECTED-test-output.txt"
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

이것은 워크플로우 코드와 함께 테스트 데이터 파일을 배포하는 한, 누구나 명령줄이나 매개변수 파일을 통해 자신의 입력을 제공하지 않고도 워크플로우를 빠르게 시험해 볼 수 있음을 의미합니다.

!!! tip "팁"

    외부에 저장된 더 큰 파일에 대한 URL을 가리킬 수 있습니다.
    열린 연결이 있는 한 Nextflow가 자동으로 다운로드합니다.

    자세한 내용은 사이드 퀘스트 [파일 작업](../side_quests/working_with_files.md)을 참조하십시오.

### 6.3. `nextflow config`를 사용하여 해결된 구성 확인

위에서 언급했듯이 때때로 동일한 매개변수가 결합하려는 프로필에서 다른 값으로 설정될 수 있습니다.
그리고 더 일반적으로, 구성 요소가 저장될 수 있는 수많은 장소가 있으며, 때때로 동일한 속성이 다른 장소에서 다른 값으로 설정될 수 있습니다.

Nextflow는 충돌을 해결하기 위해 설정된 [우선순위](https://nextflow.io/docs/latest/config.html)를 적용하지만 직접 결정하기 어려울 수 있습니다.
그리고 충돌이 없더라도 구성할 수 있는 모든 가능한 장소를 조회하는 것은 지루할 수 있습니다.

다행히도 Nextflow에는 전체 프로세스를 자동화할 수 있는 `config`라는 편리한 유틸리티 도구가 포함되어 있습니다.

`config` 도구는 현재 작업 디렉토리의 모든 내용을 탐색하고, 모든 구성 파일을 수집하고, Nextflow가 워크플로우를 실행하는 데 사용할 완전히 해결된 구성을 생성합니다.
이렇게 하면 아무것도 시작하지 않고도 어떤 설정이 사용될지 알 수 있습니다.

#### 6.3.1. 기본 구성 해결

기본적으로 적용될 구성을 해결하려면 이 명령을 실행하십시오.

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

    outputDir = 'custom-outdir-config/'

    workflow {
      output {
          mode = 'copy'
      }
    }
    ```

이것은 명령줄에서 추가 항목을 지정하지 않으면 얻는 기본 구성을 보여줍니다.

#### 6.3.2. 특정 설정이 활성화된 구성 해결

명령줄 매개변수(예: 하나 이상의 프로필 활성화 또는 매개변수 파일 로드)를 제공하면 명령이 이를 추가로 고려합니다.

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

    outputDir = 'custom-outdir-config/'

    workflow {
      output {
          mode = 'copy'
      }
    }
    ```

이것은 여러 계층의 구성을 포함하는 복잡한 프로젝트에 특히 유용합니다.

### 핵심 정리

프로필을 사용하여 최소한의 번거로움으로 런타임에 사전 설정 구성을 선택하는 방법을 알게 되었습니다.
더 일반적으로, 다른 컴퓨팅 플랫폼에 맞게 워크플로우 실행을 구성하고 분석의 재현성을 향상시키는 방법을 알게 되었습니다.

### 다음 단계

축하하고 자신에게 큰 칭찬을 해주십시오! 첫 번째 Nextflow 개발자 과정을 완료했습니다.

배운 내용을 검토하고 다음에 무엇이 올지 알아보려면 최종 [과정 요약](./next_steps.md)으로 이동하십시오.

---

## 퀴즈

<quiz>
Nextflow가 자동으로 로드하는 구성 파일의 이름은 무엇입니까?
- [ ] `config.nf`
- [ ] `pipeline.config`
- [x] `nextflow.config`
- [ ] `workflow.config`
</quiz>

<quiz>
동일한 매개변수가 구성 파일과 명령줄 모두에 설정된 경우 어떤 것이 우선합니까?
- [ ] 구성 파일 값
- [x] 명령줄 값
- [ ] 먼저 만난 값
- [ ] 둘 다 아님; 오류가 발생함

자세히 알아보기: [1.1. 기본값을 `nextflow.config`로 이동](#11-move-default-values-to-nextflowconfig)
</quiz>

<quiz>
동일한 구성에서 Docker와 Conda를 모두 활성화할 수 있습니까?
- [x] 예, Nextflow는 프로세스 지시문에 따라 둘 다 사용할 수 있음
- [ ] 아니요, 한 번에 하나만 활성화할 수 있음
- [ ] 예, 하지만 프로필에서만 가능
- [ ] 아니요, 상호 배타적임
</quiz>

<quiz>
Docker와 Conda가 모두 활성화되어 있고 프로세스에 두 지시문이 모두 있는 경우 어떤 것이 우선합니까?
- [x] Docker (컨테이너)
- [ ] Conda
- [ ] 먼저 정의된 것
- [ ] 오류가 발생함

자세히 알아보기: [3. 소프트웨어 패키징 기술 선택](#3-select-a-software-packaging-technology)
</quiz>

<quiz>
Nextflow 프로세스의 기본 메모리 할당은 얼마입니까?
- [ ] 1 GB
- [x] 2 GB
- [ ] 4 GB
- [ ] 제한 없음
</quiz>

<quiz>
구성 파일에서 특정 프로세스에 대한 리소스 요구 사항을 어떻게 설정합니까?
- [ ] `#!groovy processName.memory = '4 GB'`
- [ ] `#!groovy process.memory.processName = '4 GB'`
- [x] `#!groovy process { withName: 'processName' { memory = '4 GB' } }`
- [ ] `#!groovy resources.processName.memory = '4 GB'`

자세히 알아보기: [5.3. 특정 프로세스에 대한 리소스 할당 설정](#53-set-resource-allocations-for-a-specific-process)
</quiz>

<quiz>
리소스 활용 보고서를 생성하는 명령줄 옵션은 무엇입니까?
- [ ] `-with-metrics`
- [ ] `-with-stats`
- [x] `-with-report`
- [ ] `-with-profile`

자세히 알아보기: [5.1. 리소스 활용 보고서를 생성하기 위해 워크플로우 실행](#51-run-the-workflow-to-generate-a-resource-utilization-report)
</quiz>

<quiz>
`resourceLimits` 지시문은 무엇을 합니까?
- [ ] 최소 리소스 요구 사항 설정
- [ ] 프로세스에 리소스 할당
- [x] 요청할 수 있는 최대 리소스 제한
- [ ] 리소스 사용 모니터링

자세히 알아보기: [5.5. 리소스 제한 추가](#55-add-resource-limits)
</quiz>

<quiz>
Nextflow의 기본 실행자는 무엇입니까?
- [x] `local`
- [ ] `slurm`
- [ ] `kubernetes`
- [ ] `aws`

자세히 알아보기: [4. 실행 플랫폼 선택](#4-select-an-execution-platform)
</quiz>

<quiz>
Nextflow를 실행할 때 매개변수 파일을 어떻게 지정합니까?
- [ ] `--params params.json`
- [ ] `-config params.json`
- [x] `-params-file params.json`
- [ ] `--input params.json`

자세히 알아보기: [1.3. 매개변수 파일 사용](#13-use-a-parameter-file)
</quiz>

<quiz>
프로필은 무엇에 사용할 수 있습니까? (해당되는 것 모두 선택)
- [x] 인프라별 설정 정의
- [x] 다른 환경에 대한 리소스 제한 설정
- [x] 테스트 매개변수 제공
- [ ] 새 프로세스 정의

자세히 알아보기: [6. 프로필을 사용하여 사전 설정 구성 간 전환](#6-use-profiles-to-switch-between-preset-configurations)
</quiz>

<quiz>
단일 명령에서 여러 프로필을 어떻게 지정합니까?
- [ ] `-profile profile1 -profile profile2`
- [ ] `-profiles profile1,profile2`
- [x] `-profile profile1,profile2`
- [ ] `--profile profile1 --profile profile2`

자세히 알아보기: [6. 프로필을 사용하여 사전 설정 구성 간 전환](#6-use-profiles-to-switch-between-preset-configurations)
</quiz>
