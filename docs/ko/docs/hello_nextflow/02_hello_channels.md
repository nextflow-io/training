# 파트 2: Hello Channels

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } AI 지원 번역 - [자세히 알아보기 및 개선 제안](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

<div class="video-wrapper">
  <iframe width="560" height="315" src="https://www.youtube.com/embed/yDR66fzAMOg?si=xCItHLiOQWqoqBB9&amp;list=PLPZ8WHdZGxmWKozQuzr27jyMGqp9kElVK&amp;cc_load_policy=1&amp;cc_lang_pref=ko" title="YouTube video player" frameborder="0" allow="accelerometer; autoplay; clipboard-write; encrypted-media; gyroscope; picture-in-picture; web-share" referrerpolicy="strict-origin-when-cross-origin" allowfullscreen></iframe>
</div>

/// caption
:fontawesome-brands-youtube:{ .youtube } Nextflow YouTube 채널에서 [전체 재생목록](https://youtube.com/playlist?list=PLPZ8WHdZGxmWKozQuzr27jyMGqp9kElVK&si=eF7cLR62goy-lc6n)을 확인하십시오.

:green_book: 비디오 스크립트는 [여기](./transcripts/02_hello_channels.md)에서 확인하실 수 있습니다.
///

이 과정의 파트 1(Hello World)에서는 process 호출에서 직접 입력을 제공하여 process에 변수 입력을 전달하는 방법을 보여드렸습니다: `sayHello(params.input)`.
이것은 의도적으로 단순화된 접근 방식이었습니다.
실제로 이 접근 방식에는 주요한 제한 사항이 있습니다. 즉, 단일 값에 대해 process를 한 번만 실행하려는 매우 간단한 경우에만 작동합니다.
대부분의 실제 workflow 사용 사례에서는 여러 값(예: 여러 샘플의 실험 데이터)을 처리하려고 하므로 입력을 처리하는 더 정교한 방법이 필요합니다.

이것이 바로 Nextflow [**channel**](https://nextflow.io/docs/latest/channel.html)의 역할입니다.
Channel은 입력을 효율적으로 처리하고 다단계 workflow에서 한 단계에서 다른 단계로 입력을 전달하도록 설계된 큐이며, 내장된 병렬 처리와 많은 추가적인 이점을 제공합니다.

이 파트에서는 다양한 소스에서 여러 입력을 처리하기 위해 channel을 사용하는 방법을 배웁니다.
또한 필요에 따라 channel 내용을 변환하기 위해 [**연산자**](https://nextflow.io/docs/latest/reference/operator.html)를 사용하는 방법도 배웁니다.

??? info "이 섹션부터 시작하는 방법"

    이 섹션은 [Hello Nextflow](./index.md) 과정의 파트 1을 완료했다고 가정하지만, 해당 섹션에서 다룬 기본 사항에 익숙하다면 특별한 준비 없이 여기서 시작할 수 있습니다.

---

## 0. 준비 운동: `hello-channels.nf` 실행

시작점으로 workflow 스크립트 `hello-channels.nf`를 사용할 것입니다.
이 스크립트는 이 교육 과정의 파트 1을 완료하여 생성된 스크립트와 동일하지만, 출력 대상을 변경했습니다:

```groovy title="hello-channels.nf" linenums="37" hl_lines="3"
output {
    first_output {
        path 'hello_channels'
        mode 'copy'
    }
}
```

변경을 시작하기 전에 모든 것이 제대로 작동하는지 확인하기 위해 스크립트를 한 번 실행하십시오:

```bash
nextflow run hello-channels.nf --input 'Hello Channels!'
```

??? success "명령 출력"

    ```console
     N E X T F L O W   ~  version 25.10.2

    Launching `hello-channels.nf` [wise_jennings] DSL2 - revision: b24f4902d6

    executor >  local (1)
    [6f/824bc1] process > sayHello [100%] 1 of 1 ✔
    ```

이전과 마찬가지로, 위에 표시된 workflow 스크립트의 `output` 블록에 지정된 대로 `results/hello_channels` 디렉토리에서 `output.txt`라는 이름의 출력 파일을 찾을 수 있습니다.

??? abstract "디렉토리 내용"

    ```console title="results/hello_channels" hl_lines="2-3"
    results
    ├── hello_channels
    │   └── output.txt
    ├── hello_world
    │   └── output.txt
    └── output.txt -> /workspaces/training/hello-nextflow/work/8c/79499c11beea6e9d43605141f2817f/output.txt
    ```

??? abstract "파일 내용"

    ```console title="results/hello_channels/output.txt"
    Hello Channels!
    ```

이것이 정상적으로 작동했다면 channel에 대해 배울 준비가 되었습니다.

---

## 1. Channel을 통해 명시적으로 변수 입력 제공

암묵적 처리에 의존하는 대신, `sayHello()` process에 변수 입력을 전달하기 위해 **channel**을 생성할 것입니다. 암묵적 처리에는 특정 제한 사항이 있습니다.

### 1.1. 입력 channel 생성

Channel을 설정하는 데 사용할 수 있는 다양한 [**channel 팩토리**](https://nextflow.io/docs/latest/reference/channel.html)가 있습니다.
지금은 간단하게 유지하기 위해 단일 값을 포함하는 channel을 생성하는 [`channel.of`](https://nextflow.io/docs/latest/reference/channel.html#of)라는 가장 기본적인 channel 팩토리를 사용할 것입니다.
기능적으로는 이전에 설정한 방식과 유사하지만, Nextflow가 암묵적으로 channel을 생성하는 대신 명시적으로 수행합니다.

사용할 코드 라인입니다:

```console title="구문"
greeting_ch = channel.of('Hello Channels!')
```

이것은 `channel.of()` channel 팩토리를 사용하여 `greeting_ch`라는 channel을 생성하며, 간단한 queue channel을 설정하고 인사 값으로 사용할 문자열 `'Hello Channels!'`를 로드합니다.

<figure class="excalidraw">
--8<-- "docs/en/docs/hello_nextflow/img/hello-pipeline-channel.svg"
</figure>

!!! note

    가독성을 위해 CLI 매개변수 대신 임시로 하드코딩된 문자열로 전환합니다. Channel 수준에서 무슨 일이 일어나는지 다룬 후에 CLI 매개변수 사용으로 돌아갈 것입니다.

workflow 블록에서 channel 팩토리 코드를 추가하십시오:

=== "수정 후"

    ```groovy title="hello-channels.nf" linenums="27" hl_lines="4 5"
    workflow {

        main:
        // 입력용 channel 생성
        greeting_ch = channel.of('Hello Channels!')
        // 인사말을 내보냅니다
        sayHello(params.input)

        publish:
        first_output = sayHello.out
    }
    ```

=== "수정 전"

    ```groovy title="hello-channels.nf" linenums="27"
    workflow {

        main:
        // 인사말을 내보냅니다
        sayHello(params.input)

        publish:
        first_output = sayHello.out
    }
    ```

아직 process 호출의 입력을 전환하지 않았으므로 아직 작동하지 않습니다.

### 1.2. Process 호출에 channel을 입력으로 추가

이제 새로 생성한 channel을 `sayHello()` process 호출에 연결하여 이전에 직접 제공하던 CLI 매개변수를 대체해야 합니다.

workflow 블록에서 다음 코드 변경을 수행하십시오:

=== "수정 후"

    ```groovy title="hello-channels.nf" linenums="27" hl_lines="7"
    workflow {

        main:
        // 입력용 channel 생성
        greeting_ch = channel.of('Hello Channels!')
        // 인사말을 내보냅니다
        sayHello(greeting_ch)

        publish:
        first_output = sayHello.out
    }
    ```

=== "수정 전"

    ```groovy title="hello-channels.nf" linenums="27" hl_lines="7"
    workflow {

        main:
        // 입력용 channel 생성
        greeting_ch = channel.of('Hello Channels!')
        // 인사말을 내보냅니다
        sayHello(params.input)

        publish:
        first_output = sayHello.out
    }
    ```

이것은 Nextflow에게 `greeting_ch` channel의 내용에 대해 `sayHello` process를 실행하도록 지시합니다.

이제 workflow가 제대로 작동합니다. 이것은 `sayHello('Hello Channels!')`를 작성하는 것과 동일한 명시적 표현입니다.

### 1.3. Workflow 실행

실행해 봅시다!

```bash
nextflow run hello-channels.nf
```

??? success "명령 출력"

    ```console
     N E X T F L O W   ~  version 25.10.2

    Launching `hello-channels.nf` [fabulous_crick] DSL2 - revision: 23e20f76e8

    executor >  local (1)
    [c0/4f1872] process > sayHello (1) [100%] 1 of 1 ✔
    ```

두 편집을 모두 올바르게 수행했다면 성공적으로 실행될 것입니다.
results 디렉토리를 확인하여 결과가 이전과 동일한지 확인할 수 있습니다.

??? abstract "파일 내용"

    ```console title="results/hello_channels/output.txt"
    Hello Channels!
    ```

동일한 최종 결과를 얻으면서 workflow의 유연성을 높였습니다.
더 많은 코드를 작성하면서 실질적인 이점이 없는 것처럼 보일 수 있지만, 더 많은 입력을 처리하기 시작하면 그 가치가 명확해질 것입니다.

다음으로 넘어가기 전에 한 가지 더 살펴보겠습니다. 데이터 입력 관리를 위해 명시적 channel을 사용하는 것의 작지만 편리한 이점입니다.

### 1.4. `view()`를 사용하여 channel 내용 검사

Nextflow channel은 연산자를 사용하여 내용을 조작할 수 있도록 설계되어 있으며, 이에 대해서는 이 장의 후반부에서 자세히 다룰 것입니다.

지금은 channel의 내용을 검사하기 위해 [`view()`](https://www.nextflow.io/docs/latest/reference/operator.html#view)라는 매우 간단한 연산자를 사용하는 방법만 보여드리겠습니다.
`view()`를 Python의 `print()` 문이나 다른 언어의 동등한 것과 같은 디버깅 도구로 생각할 수 있습니다.

workflow 블록에 이 작은 라인을 추가하십시오:

=== "수정 후"

    ```groovy title="hello-channels.nf" linenums="27" hl_lines="7"
    workflow {

        main:
        // 입력용 channel 생성
        greeting_ch = channel.of('Hello Channels!')
                             .view()
        // 인사말을 내보냅니다
        sayHello(greeting_ch)

        publish:
        first_output = sayHello.out
    }
    ```

=== "수정 전"

    ```groovy title="hello-channels.nf" linenums="27"
    workflow {

        main:
        // 입력용 channel 생성
        greeting_ch = channel.of('Hello Channels!')
        // 인사말을 내보냅니다
        sayHello(greeting_ch)

        publish:
        first_output = sayHello.out
    }
    ```

정확한 공백 수는 4의 배수인 한 중요하지 않습니다. `.view()` 문의 시작을 channel 구성의 `.of()` 부분에 맞추는 것이 목표입니다.

이제 workflow를 다시 실행하십시오:

```bash
nextflow run hello-channels.nf
```

??? success "명령 출력"

    ```console hl_lines="7"
     N E X T F L O W   ~  version 25.10.2

    Launching `hello-channels.nf` [scruffy_shaw] DSL2 - revision: 2ede41e14a

    executor >  local (1)
    [ef/f7e40a] sayHello (1) [100%] 1 of 1 ✔
    Hello Channels!
    ```

보시다시피 channel 내용이 콘솔에 출력됩니다.
여기서는 요소가 하나뿐이지만, 다음 섹션에서 channel에 여러 값을 로드하기 시작하면 줄당 하나의 요소가 출력되도록 설정되어 있는 것을 볼 수 있습니다.

### 핵심 정리

기본 channel 팩토리를 사용하여 process에 입력을 제공하는 방법을 알게 되었습니다.

### 다음 단계

Channel을 사용하여 workflow가 여러 입력 값을 반복 처리하도록 하는 방법을 배웁니다.

---

## 2. Workflow를 수정하여 여러 입력 값에서 실행

Workflow는 일반적으로 대량으로 처리해야 하는 입력 배치에서 실행되므로, workflow를 업그레이드하여 여러 입력 값을 허용하도록 해야 합니다.

### 2.1. 입력 channel에 여러 인사말 로드

편리하게도, 사용해 온 `channel.of()` channel 팩토리는 둘 이상의 값을 기꺼이 허용하므로 전혀 수정할 필요가 없습니다.
Channel에 여러 값을 로드하기만 하면 됩니다.

`'Hello'`, `'Bonjour'`, `'Holà'`로 만들어 봅시다.

#### 2.1.1. 더 많은 인사말 추가

workflow 블록 전에 다음 코드 변경을 수행하십시오:

=== "수정 후"

    ```groovy title="hello-channels.nf" linenums="30" hl_lines="2"
    // 입력용 channel 생성
    greeting_ch = channel.of('Hello','Bonjour','Holà')
                         .view()
    ```

=== "수정 전"

    ```groovy title="hello-channels.nf" linenums="30" hl_lines="2"
    // 입력용 channel 생성
    greeting_ch = channel.of('Hello Channels')
                         .view()
    ```

문서에 따르면 이것이 작동해야 합니다. 정말 이렇게 간단할 수 있을까요?

#### 2.1.2. 명령을 실행하고 로그 출력 확인

시도해 봅시다.

```bash
nextflow run hello-channels.nf
```

??? success "명령 출력"

    ```console hl_lines="6"
     N E X T F L O W   ~  version 25.10.2

    Launching `hello-channels.nf` [amazing_crick] DSL2 - revision: 59a9a5888a

    executor >  local (3)
    [f4/c9962c] process > sayHello (1) [100%] 3 of 3 ✔
    Hello
    Bonjour
    Holà
    ```

확실히 잘 실행된 것 같습니다.
실행 모니터는 `sayHello` process에 대해 `3 of 3` 호출이 이루어졌음을 보여주며, 약속한 대로 `view()` 문에 의해 세 개의 인사말이 줄당 하나씩 나열되어 있습니다.

그러나 results 디렉토리에는 여전히 하나의 출력만 있습니다:

??? abstract "디렉토리 내용"

    ```console title="results/hello_channels" hl_lines="3"
    results
    ├── hello_channels
    │   └── output.txt
    ├── hello_world
    │   └── output.txt
    └── output.txt -> /workspaces/training/hello-nextflow/work/8c/79499c11beea6e9d43605141f2817f/output.txt
    ```

??? abstract "파일 내용"

    ```console title="results/hello_channels/output.txt"
    Holà
    ```

세 개의 인사말 중 하나가 있어야 하지만, 여기에 표시된 것과 다를 수 있습니다.
왜 그럴 수 있는지 생각해 보셨습니까?

<figure class="excalidraw">
--8<-- "docs/en/docs/hello_nextflow/img/hello-pipeline-channel-multi.svg"
</figure>

_다이어그램에서 channel은 녹색으로 표시되며, 요소의 순서는 파이프 안의 구슬처럼 표현됩니다. 처음 로드된 것이 오른쪽에, 두 번째가 중간에, 세 번째가 왼쪽에 있습니다._

실행 모니터를 다시 살펴보면, 하나의 하위 디렉토리 경로(`f4/c9962c`)만 제공했습니다.
그 안을 살펴봅시다.

??? abstract "디렉토리 내용"

    ```console hl_lines="9"
    work/f4/c9962ce91ef87480babcb86b2b9042/
    ├── .command.begin
    ├── .command.err
    ├── .command.log
    ├── .command.out
    ├── .command.run
    ├── .command.sh
    ├── .exitcode
    └── output.txt
    ```

??? abstract "파일 내용"

    ```console title="work/f4/c9962ce91ef87480babcb86b2b9042/output.txt"
    Hello
    ```

results 디렉토리에서 얻은 것과 같은 인사말도 아닙니다! 무슨 일이 일어나고 있는 걸까요?

이 시점에서 기본적으로 ANSI 로깅 시스템은 동일한 process에 대한 여러 호출의 로깅을 같은 줄에 기록한다는 것을 알려드려야 합니다.
따라서 sayHello() process에 대한 세 번의 호출 상태가 모두 같은 위치에 표시됩니다.

다행히도 전체 process 호출 목록을 볼 수 있도록 해당 동작을 비활성화할 수 있습니다.

#### 2.1.3. `-ansi-log false` 옵션으로 명령을 다시 실행

로깅을 확장하여 process 호출당 한 줄씩 표시하려면 명령에 `-ansi-log false`를 추가하십시오.

```bash
nextflow run hello-channels.nf -ansi-log false
```

??? success "명령 출력"

    ```console
     N E X T F L O W  ~  version 25.10.2
    Launching `hello-channels.nf` [desperate_monod] DSL2 - revision: 59a9a5888a
    Hello
    Bonjour
    Holà
    [23/871c7e] Submitted process > sayHello (2)
    [7f/21e2c2] Submitted process > sayHello (1)
    [f4/ea10a6] Submitted process > sayHello (3)
    ```

이번에는 세 개의 process 실행과 관련 작업 하위 디렉토리가 모두 출력에 나열됩니다.

훨씬 낫습니다. 적어도 간단한 workflow에서는요.
복잡한 workflow나 많은 수의 입력의 경우 전체 목록이 터미널에 출력되면 다소 압도적일 수 있습니다.
그래서 `-ansi-log false`가 기본 동작이 아닙니다.

!!! tip

    상태가 보고되는 방식은 두 로깅 모드 간에 약간 다릅니다.
    압축 모드에서 Nextflow는 호출이 성공적으로 완료되었는지 여부를 보고합니다.
    이 확장 모드에서는 제출되었다는 것만 보고합니다.

어쨌든 이제 각 process 호출의 하위 디렉토리가 있으므로 로그와 출력을 찾을 수 있습니다.

??? abstract "디렉토리 내용"

    ```console
    work/23/871c7ec3642a898ecd5e6090d21300/
    ├── .command.begin
    ├── .command.err
    ├── .command.log
    ├── .command.out
    ├── .command.run
    ├── .command.sh
    ├── .exitcode
    └── output.txt
    ```

    ```console
    work/7f/21e2c2f3cc8833ef3858b236e5575c/
    ├── .command.begin
    ├── .command.err
    ├── .command.log
    ├── .command.out
    ├── .command.run
    ├── .command.sh
    ├── .exitcode
    └── output.txt
    ```

    ```console
    work/f4/ea10a680d5687596d3eaa3fcf69272/
    ├── .command.begin
    ├── .command.err
    ├── .command.log
    ├── .command.out
    ├── .command.run
    ├── .command.sh
    ├── .exitcode
    └── output.txt
    ```

??? abstract "파일 내용"

    ```txt title="work/23/871c7ec3642a898ecd5e6090d21300/output.txt"
    Bonjour
    ```

    ```txt title="work/7f/21e2c2f3cc8833ef3858b236e5575c/output.txt"
    Hello
    ```

    ```txt title="work/f4/ea10a680d5687596d3eaa3fcf69272/output.txt"
    Holà
    ```

이것은 세 개의 process가 모두 성공적으로 실행되었음을 보여줍니다.

그렇긴 하지만, results 디렉토리에 여전히 하나의 출력 파일만 있다는 문제가 있습니다.

`sayHello` process의 출력 파일 이름을 하드코딩했기 때문에 세 번의 호출 모두 `output.txt`라는 파일을 생성했던 것을 기억하실 것입니다.

<figure class="excalidraw">
--8<-- "docs/en/docs/hello_nextflow/img/hello-channels-task-dirs.svg"
</figure>

출력 파일이 다른 process와 격리된 작업 하위 디렉토리에 있는 한 괜찮습니다.
그러나 동일한 results 디렉토리에 게시되면 먼저 복사된 것이 다음 것에 의해 덮어쓰여지고, 계속 반복됩니다.

### 2.2. 출력 파일 이름이 고유하도록 보장

모든 출력을 동일한 results 디렉토리에 계속 게시할 수 있지만 고유한 이름을 갖도록 해야 합니다.
구체적으로, 최종 파일 이름이 고유하도록 동적으로 파일 이름을 생성하도록 첫 번째 process를 수정해야 합니다.

그러면 파일 이름을 어떻게 고유하게 만들까요?
일반적인 방법은 입력(입력 channel에서 받은)의 고유한 메타데이터 조각을 출력 파일 이름의 일부로 사용하는 것입니다.
여기서는 편의상 인사말 자체가 짧은 문자열이므로 그것을 사용하여 기본 출력 파일명 앞에 추가하겠습니다.

<figure class="excalidraw">
--8<-- "docs/en/docs/hello_nextflow/img/hello-pipeline-channel-multi-unique.svg"
</figure>

#### 2.2.1. 동적 출력 파일 이름 구성

process 블록에서 다음 코드 변경을 수행하십시오:

=== "수정 후"

    ```groovy title="hello-channels.nf" linenums="6" hl_lines="7 11"
    process sayHello {

        input:
        val greeting

        output:
        path "${greeting}-output.txt"

        script:
        """
        echo '${greeting}' > '${greeting}-output.txt'
        """
    }
    ```

=== "수정 전"

    ```groovy title="hello-channels.nf" linenums="6" hl_lines="7 11"
    process sayHello {

        input:
        val greeting

        output:
        path 'output.txt'

        script:
        """
        echo '${greeting}' > output.txt
        """
    }
    ```

출력 정의와 `script:` 명령 블록 모두에서 `output.txt`를 교체해야 합니다.

!!! tip

    출력 정의에서 출력 파일명 표현식 주위에 큰따옴표를 사용해야 합니다(작은따옴표가 아님). 그렇지 않으면 실패합니다.

이렇게 하면 process가 호출될 때마다 고유한 출력 파일 이름이 생성되어 출력 디렉토리에서 동일한 process에 대한 다른 호출의 출력과 구별될 수 있습니다.

#### 2.2.2. Workflow 실행

실행해 봅시다. 기본 ANSI 로그 설정으로 다시 실행합니다.

```bash
nextflow run hello-channels.nf
```

??? success "명령 출력"

    ```console
     N E X T F L O W   ~  version 25.10.2

    Launching `hello-channels.nf` [sharp_minsky] DSL2 - revision: 16a291febe

    executor >  local (3)
    [e8/33ee64] sayHello (2) [100%] 3 of 3 ✔
    Hello
    Bonjour
    Holà
    ```

요약 보기로 돌아가면 출력이 다시 한 줄로 요약됩니다.
`results` 디렉토리를 확인하여 모든 출력 인사말이 있는지 확인하십시오.

??? abstract "디렉토리 내용"

    ```console
    results/hello_channels/
    ├── Bonjour-output.txt
    ├── Hello-output.txt
    ├── Holà-output.txt
    └── output.txt
    ```

그렇습니다! 각각 예상한 내용이 들어 있습니다.

??? abstract "파일 내용"

    ```console title="Bonjour-output.txt"
    Bonjour
    ```

    ```console title="Hello-output.txt"
    Hello
    ```

    ```console title="Holà-output.txt"
    Holà
    ```

성공입니다! 이제 출력 파일이 덮어쓰여지는 것에 대해 걱정하지 않고 원하는 만큼 많은 인사말을 추가할 수 있습니다.

!!! tip

    실제로 입력 데이터 자체를 기반으로 파일 이름을 지정하는 것은 거의 항상 비실용적입니다.
    동적 파일 이름을 생성하는 더 좋은 방법은 입력 파일과 함께 메타데이터를 process에 전달하는 것입니다.
    메타데이터는 일반적으로 'sample sheet' 또는 이와 동등한 것을 통해 제공됩니다.
    나중에 Nextflow 교육에서 이를 수행하는 방법을 배우게 됩니다([메타데이터 사이드 퀘스트](../side_quests/metadata.md) 참조).

### 핵심 정리

Channel을 통해 여러 입력 요소를 공급하는 방법을 알게 되었습니다.

### 다음 단계

연산자를 사용하여 channel의 내용을 변환하는 방법을 배웁니다.

---

## 3. 배열을 통해 여러 입력 제공

방금 channel 팩토리에 직접 하드코딩된 여러 입력 요소를 처리하는 방법을 보여드렸습니다.
다른 방식으로 여러 입력을 제공하려면 어떻게 해야 할까요?

예를 들어, 다음과 같이 요소 배열을 포함하는 입력 변수를 설정한다고 상상해 보십시오:

`greetings_array = ['Hello','Bonjour','Holà']`

이것을 출력 channel에 로드하고 작동할 것으로 기대할 수 있을까요?

<figure class="excalidraw">
--8<-- "docs/en/docs/hello_nextflow/img/hello-pipeline-multi-inputs-array.svg"
</figure>

알아봅시다.

### 3.1. Channel에 값 배열을 입력으로 제공

상식적으로 단일 값 대신 값 배열을 간단히 전달할 수 있어야 합니다.
시도해 봅시다. 입력 변수를 설정하고 channel 팩토리에 로드해야 합니다.

#### 3.1.1. 입력 변수 설정

방금 상상한 `greetings_array` 변수를 workflow 블록에 추가하여 현실로 만들어 봅시다:

=== "수정 후"

    ```groovy title="hello-channels.nf" linenums="27" hl_lines="4 5"
    workflow {

        main:
        // 입력 인사말 배열 선언
        greetings_array = ['Hello','Bonjour','Holà']
        // 입력용 channel 생성
        greeting_ch = channel.of('Hello','Bonjour','Holà')
                             .view()
        // 인사말을 내보냅니다
        sayHello(greeting_ch)

        publish:
        first_output = sayHello.out
    }
    ```

=== "수정 전"

    ```groovy title="hello-channels.nf" linenums="27"
    workflow {

        main:
        // 입력용 channel 생성
        greeting_ch = channel.of('Hello','Bonjour','Holà')
                             .view()
        // 인사말을 내보냅니다
        sayHello(greeting_ch)

        publish:
        first_output = sayHello.out
    }
    ```

아직 작동하지 않습니다. 배열에 대한 선언만 추가했습니다.

#### 3.1.2. Channel 팩토리에 인사말 배열을 입력으로 설정

이제 channel 팩토리에 현재 하드코딩된 값 `'Hello','Bonjour','Holà'`를 방금 생성한 `greetings_array`로 교체할 것입니다.

workflow 블록에서 다음 변경을 수행하십시오:

=== "수정 후"

    ```groovy title="hello-channels.nf" linenums="27" hl_lines="7"
    workflow {

        main:
        // 입력 인사말 배열 선언
        greetings_array = ['Hello','Bonjour','Holà']
        // 입력용 channel 생성
        greeting_ch = channel.of(greetings_array)
                             .view()
        // 인사말을 내보냅니다
        sayHello(greeting_ch)

        publish:
        first_output = sayHello.out
    }
    ```

=== "수정 전"

    ```groovy title="hello-channels.nf" linenums="27" hl_lines="7"
    workflow {

        main:
        // 입력 인사말 배열 선언
        greetings_array = ['Hello','Bonjour','Holà']
        // 입력용 channel 생성
        greeting_ch = channel.of('Hello','Bonjour','Holà')
                             .view()
        // 인사말을 내보냅니다
        sayHello(greeting_ch)

        publish:
        first_output = sayHello.out
    }
    ```

이제 작동해야 합니다.

#### 3.1.3. Workflow 실행

실행해 봅시다:

```bash
nextflow run hello-channels.nf
```

??? failure "명령 출력"

    ```console hl_lines="7 11 16"
    N E X T F L O W   ~  version 25.10.2

    Launching `hello-channels.nf` [friendly_koch] DSL2 - revision: 97256837a7

    executor >  local (1)
    [a8/1f6ead] sayHello (1) | 0 of 1
    [Hello, Bonjour, Holà]
    ERROR ~ Error executing process > 'sayHello (1)'

    Caused by:
      Missing output file(s) `[Hello, Bonjour, Holà]-output.txt` expected by process `sayHello (1)`


    Command executed:

      echo '[Hello, Bonjour, Holà]' > '[Hello, Bonjour, Holà]-output.txt'

    Command exit status:
      0

    Command output:
      (empty)

    Work dir:
      /workspaces/training/hello-nextflow/work/a8/1f6ead5f3fa30a3c508e2e7cf83ffb

    Tip: you can replicate the issue by changing to the process work dir and entering the command `bash .command.run`

    -- Check '.nextflow.log' file for details
    ```

오류가 발생했습니다!

`view()`의 출력과 오류 메시지를 살펴보십시오.

Nextflow가 배열의 세 문자열을 별도의 값으로 사용하는 대신 `[Hello, Bonjour, Holà]`를 단일 문자열 값으로 사용하여 단일 process 호출을 실행하려고 시도한 것 같습니다.

<figure class="excalidraw">
--8<-- "docs/en/docs/hello_nextflow/img/hello-channels-array-fail.svg"
</figure>

따라서 '패키징'이 문제를 일으키고 있습니다.
Nextflow가 배열을 풀고 개별 문자열을 channel에 로드하도록 하려면 어떻게 해야 할까요?

### 3.2. 연산자를 사용하여 channel 내용 변환

여기서 [**연산자**](https://nextflow.io/docs/latest/reference/operator.html)가 등장합니다.
내부 내용을 확인하는 `.view()` 연산자를 이미 사용했습니다.
이제 channel의 내용에 대해 작업할 수 있는 연산자를 살펴보겠습니다.

Nextflow 문서의 [연산자 목록](https://nextflow.io/docs/latest/reference/operator.html)을 훑어보면, 정확히 필요한 것을 수행하는 [`flatten()`](https://www.nextflow.io/docs/latest/reference/operator.html#flatten)을 찾을 수 있습니다. 배열의 내용을 풀고 개별 항목으로 내보냅니다.

#### 3.2.1. `flatten()` 연산자 추가

입력 channel에 `flatten()` 연산자를 적용하려면 channel 팩토리 선언에 추가합니다.

workflow 블록에서 다음 코드 변경을 수행하십시오:

=== "수정 후"

    ```groovy title="hello-channels.nf" linenums="27" hl_lines="9"
    workflow {

        main:
        // 입력 인사말 배열 선언
        greetings_array = ['Hello','Bonjour','Holà']
        // 입력용 channel 생성
        greeting_ch = channel.of(greetings_array)
                             .view()
                             .flatten()
        // 인사말을 내보냅니다
        sayHello(greeting_ch)

        publish:
        first_output = sayHello.out
    }
    ```

=== "수정 전"

    ```groovy title="hello-channels.nf" linenums="27"
    workflow {

        main:
        // 입력 인사말 배열 선언
        greetings_array = ['Hello','Bonjour','Holà']
        // 입력용 channel 생성
        greeting_ch = channel.of(greetings_array)
                             .view()
        // 인사말을 내보냅니다
        sayHello(greeting_ch)

        publish:
        first_output = sayHello.out
    }
    ```

여기서 가독성을 위해 다음 줄에 연산자를 추가했지만, 원한다면 다음과 같이 channel 팩토리와 같은 줄에 연산자를 추가할 수 있습니다:
`greeting_ch = channel.of(greetings_array).view().flatten()`

<figure class="excalidraw">
--8<-- "docs/en/docs/hello_nextflow/img/hello-channels-array-success.svg"
</figure>

#### 3.2.2. `view()` 문 개선

이것이 작동하는지 테스트하기 위해 바로 실행할 수 있지만, 그러는 김에 channel 내용을 검사하는 방식을 개선하겠습니다.

`flatten()` 연산자가 적용되기 전과 후의 내용이 어떻게 보이는지 대조할 수 있도록 두 번째 것을 추가하고, 출력에서 더 명확하게 레이블이 지정되도록 약간의 코드를 추가하겠습니다.

workflow 블록에서 다음 코드 변경을 수행하십시오:

=== "수정 후"

    ```groovy title="hello-channels.nf" linenums="27" hl_lines="8-10"
    workflow {

        main:
        // 입력 인사말 배열 선언
        greetings_array = ['Hello','Bonjour','Holà']
        // 입력용 channel 생성
        greeting_ch = channel.of(greetings_array)
                             .view { greeting -> "Before flatten: $greeting" }
                             .flatten()
                             .view { greeting -> "After flatten: $greeting" }
        // 인사말을 내보냅니다
        sayHello(greeting_ch)

        publish:
        first_output = sayHello.out
    }
    ```

=== "수정 전"

    ```groovy title="hello-channels.nf" linenums="27" hl_lines="8-9"
    workflow {

        main:
        // 입력 인사말 배열 선언
        greetings_array = ['Hello','Bonjour','Holà']
        // 입력용 channel 생성
        greeting_ch = channel.of(greetings_array)
                             .view()
                             .flatten()
        // 인사말을 내보냅니다
        sayHello(greeting_ch)

        publish:
        first_output = sayHello.out
    }
    ```

두 번째 `.view` 문을 추가했으며, 각각에 대해 빈 괄호(`()`)를 `{ greeting -> "Before flatten: $greeting" }`와 같은 코드가 포함된 중괄호로 교체했습니다.

이것을 *클로저*라고 합니다. 포함된 코드는 channel의 각 항목에 대해 실행됩니다.
해당 클로저의 범위 내에서만 사용되는 내부 값에 대한 임시 변수를 정의합니다. 여기서는 `greeting`이라고 합니다(임의의 이름일 수 있음).

이 예에서 `$greeting`은 channel에 로드된 각 개별 항목을 나타냅니다.
이렇게 하면 깔끔하게 레이블이 지정된 콘솔 출력이 생성됩니다.

!!! info

    일부 파이프라인에서는 연산자 클로저 내에서 `$it`라는 특수 변수가 사용되는 것을 볼 수 있습니다.
    이것은 `->`로 정의할 필요 없이 내부 변수에 대한 단축 액세스를 허용하는 _암묵적_ 변수입니다.

    코드 명확성을 위해 명시적으로 표현하는 것을 선호하므로, `$it` 구문은 권장되지 않으며 Nextflow 언어에서 점차 단계적으로 제거될 것입니다.

#### 3.2.3. Workflow 실행

마지막으로 workflow를 다시 실행해 볼 수 있습니다!

```bash
nextflow run hello-channels.nf
```

??? success "명령 출력"

    ```console hl_lines="7-10"
     N E X T F L O W   ~  version 25.10.2

    Launching `hello-channels.nf` [sleepy_gutenberg] DSL2 - revision: 1db4f760ee

    executor >  local (3)
    [b1/6a1e15] sayHello (2) [100%] 3 of 3 ✔
    Before flatten: [Hello, Bonjour, Holà]
    After flatten: Hello
    After flatten: Bonjour
    After flatten: Holà
    ```

이번에는 작동하며 `flatten()` 연산자를 실행하기 전과 후에 channel의 내용이 어떻게 보이는지에 대한 추가 통찰력을 제공합니다.

- 해당 시점에서 channel에는 원래 배열인 하나의 항목이 포함되어 있으므로 단일 `Before flatten:` 문이 표시됩니다.
- 그런 다음 이제 channel의 개별 항목인 각 인사말에 대해 세 개의 별도 `After flatten:` 문이 표시됩니다.

중요한 것은 이제 각 항목이 workflow에서 별도로 처리될 수 있다는 것입니다.

!!! tip

    작동에 암묵적 매핑 단계를 포함하는 다른 channel 팩토리인 [`channel.fromList`](https://nextflow.io/docs/latest/reference/channel.html#fromlist)를 사용하여 기술적으로 동일한 결과를 얻을 수 있습니다.
    여기서는 간단한 사용 사례에서 연산자 사용을 시연하기 위해 이를 사용하지 않기로 했습니다.

### 핵심 정리

`flatten()`과 같은 연산자를 사용하여 channel의 내용을 변환하는 방법과 연산자를 적용하기 전과 후에 channel 내용을 검사하기 위해 `view()` 연산자를 사용하는 방법을 알게 되었습니다.

### 다음 단계

Workflow가 파일을 입력 값의 소스로 사용하도록 하는 방법을 배웁니다.

---

## 4. CSV 파일에서 입력 값 읽기

현실적으로 값 배열에서 시작하는 경우는 거의 없습니다.
대부분의 경우 일종의 구조화된 형식으로 처리해야 하는 데이터가 포함된 하나 이상의 파일이 있을 것입니다.

실제 데이터 분석에서 처리하고자 하는 종류의 열 형식 데이터를 모방하여 여러 입력 인사말이 포함된 `greetings.csv`라는 CSV 파일을 `data/` 아래에 준비했습니다.
(숫자는 의미가 없으며 설명 목적으로만 있습니다.)

```csv title="data/greetings.csv" linenums="1"
Hello,English,123
Bonjour,French,456
Holà,Spanish,789
```

다음 작업은 이 파일에서 값을 읽도록 workflow를 조정하는 것입니다.

<figure class="excalidraw">
--8<-- "docs/en/docs/hello_nextflow/img/hello-pipeline-multi-inputs-csv.svg"
</figure>

이를 구현하는 방법을 살펴봅시다.

### 4.1. CSV 파일을 인사말 소스로 기대하도록 스크립트 수정

시작하려면 스크립트에 두 가지 주요 변경이 필요합니다:

- 입력 매개변수를 CSV 파일을 가리키도록 전환
- Channel 팩토리를 파일을 처리하도록 설계된 것으로 전환

#### 4.1.1. 입력 매개변수를 CSV 파일을 가리키도록 전환

파트 1에서 설정한 `params.input` 매개변수를 기억하십니까?
인사말이 포함된 CSV 파일을 가리키도록 업데이트하겠습니다.

매개변수 선언을 다음과 같이 편집하십시오:

=== "수정 후"

    ```groovy title="hello-channels.nf" linenums="20" hl_lines="5"
    /*
    * 파이프라인 매개변수
    */
    params {
        input: Path = 'data/greetings.csv'
    }
    ```

=== "수정 전"

    ```groovy title="hello-channels.nf" linenums="20" hl_lines="5"
    /*
     * 파이프라인 매개변수
     */
    input: String = 'Holà mundo!'
    ```

이것은 파일이 workflow 코드와 같은 위치에 있다고 가정합니다.
나중에 Nextflow 여정에서 다른 데이터 위치를 다루는 방법을 배우게 됩니다.

#### 4.1.2. 파일을 처리하도록 설계된 channel 팩토리로 전환

이제 입력으로 단순 문자열 대신 파일을 사용하고자 하므로 이전의 `channel.of()` channel 팩토리를 사용할 수 없습니다.
파일 경로 처리를 위한 일부 내장 기능이 있는 새로운 channel 팩토리 [`channel.fromPath()`](https://nextflow.io/docs/latest/reference/channel.html#frompath)를 사용하도록 전환해야 합니다.

workflow 블록에서 다음 코드 변경을 수행하십시오:

=== "수정 후"

    ```groovy title="hello-channels.nf" linenums="27" hl_lines="4-8"
    workflow {

        main:
        // CSV 파일에서 입력용 channel 생성
        greeting_ch = channel.fromPath(params.input)
                             .view { greeting -> "Before flatten: $greeting" }
                             // .flatten()
                             // .view { greeting -> "After flatten: $greeting" }
        // 인사말을 내보냅니다
        sayHello(greeting_ch)

        publish:
        first_output = sayHello.out
    }
    ```

=== "수정 전"

    ```groovy title="hello-channels.nf" linenums="27" hl_lines="4-8"
    workflow {

        main:
        // 입력 인사말 배열 선언
        greetings_array = ['Hello','Bonjour','Holà']
        // 입력용 channel 생성
        greeting_ch = channel.of(greetings_array)
                             .view { greeting -> "Before flatten: $greeting" }
                             .flatten()
                             .view { greeting -> "After flatten: $greeting" }
        // 인사말을 내보냅니다
        sayHello(greeting_ch)

        publish:
        first_output = sayHello.out
    }
    ```

channel 입력을 다시 `param.input`으로 전환하고, 더 이상 필요하지 않으므로 `greetings_array` 선언을 삭제했습니다.
또한 `flatten()`과 두 번째 `view()` 문을 주석 처리했습니다.

#### 4.1.3. Workflow 실행

새 channel 팩토리와 입력 파일로 workflow를 실행해 봅시다.

```bash
nextflow run hello-channels.nf
```

??? failure "명령 출력"

    ```console hl_lines="5 6 9 14"
     N E X T F L O W   ~  version 25.10.2

    Launching `hello-channels.nf` [peaceful_poisson] DSL2 - revision: a286c08ad5

    [-        ] sayHello [  0%] 0 of 1
    Before flatten: /workspaces/training/hello-nextflow/data/greetings.csv
    ERROR ~ Error executing process > 'sayHello (1)'

    Caused by:
      File `/workspaces/training/hello-nextflow/data/greetings.csv-output.txt` is outside the scope of the process work directory: /workspaces/training/hello-nextflow/work/30/e610cb4ea5ae8693f456ac3329c92f


    Command executed:

      echo '/workspaces/training/hello-nextflow/data/greetings.csv' > '/workspaces/training/hello-nextflow/data/greetings.csv-output.txt'

    Command exit status:
      -

    Command output:
      (empty)

    Work dir:
      /workspaces/training/hello-nextflow/work/30/e610cb4ea5ae8693f456ac3329c92f

    Tip: when you have fixed the problem you can continue the execution adding the option `-resume` to the run command line

    -- Check '.nextflow.log' file for details
    ```

작동하지 않습니다. 콘솔 출력의 시작과 오류 메시지를 살펴보십시오.
`Command executed:` 부분이 특히 도움이 됩니다.

약간 익숙해 보일 수 있습니다.
Nextflow가 파일 경로 자체를 문자열 값으로 사용하여 단일 process 호출을 실행하려고 시도한 것 같습니다.
파일 경로를 올바르게 확인했지만 실제로 원하는 내용을 구문 분석하지 않았습니다.

Nextflow가 파일을 열고 그 내용을 channel에 로드하도록 하려면 어떻게 해야 할까요?

또 다른 [연산자](https://nextflow.io/docs/latest/reference/operator.html)가 필요한 것 같습니다!

### 4.2. `splitCsv()` 연산자를 사용하여 파일 구문 분석

연산자 목록을 다시 살펴보면, CSV 형식 텍스트를 구문 분석하고 분할하도록 설계된 [`splitCsv()`](https://nextflow.io/docs/latest/reference/operator.html#splitcsv)를 찾을 수 있습니다.

#### 4.2.1. Channel에 `splitCsv()` 적용

연산자를 적용하려면 이전처럼 channel 팩토리 라인에 추가합니다.

workflow 블록에서 `flatten()`을 `splitcsv()`로 교체(주석 해제)하는 다음 코드 변경을 수행하십시오:

=== "수정 후"

    ```groovy title="hello-channels.nf" linenums="27" hl_lines="6-8"
    workflow {

        main:
        // CSV 파일에서 입력용 channel 생성
        greeting_ch = channel.fromPath(params.input)
                             .view { csv -> "Before splitCsv: $csv" }
                             .splitCsv()
                             .view { csv -> "After splitCsv: $csv" }
        // 인사말을 내보냅니다
        sayHello(greeting_ch)

        publish:
        first_output = sayHello.out
    }
    ```

=== "수정 전"

    ```groovy title="hello-channels.nf" linenums="27" hl_lines="6-8"
    workflow {

        main:
        // CSV 파일에서 입력용 channel 생성
        greeting_ch = channel.fromPath(params.input)
                             .view { greeting -> "Before flatten: $greeting" }
                             // .flatten()
                             // .view { greeting -> "After flatten: $greeting" }
        // 인사말을 내보냅니다
        sayHello(greeting_ch)

        publish:
        first_output = sayHello.out
    }
    ```

보시다시피 전후 `view()` 문도 업데이트했습니다.
기술적으로 동일한 변수 이름(`greeting`)을 사용할 수 있었지만, 다른 사람이 코드를 더 읽기 쉽게 하기 위해 더 적절한 것(`csv`)으로 업데이트했습니다.

#### 4.2.2. Workflow를 다시 실행

추가된 CSV 구문 분석 로직으로 workflow를 실행해 봅시다.

```bash
nextflow run hello-channels.nf
```

??? failure "명령 출력"

    ```console hl_lines="7-11 14 19"
     N E X T F L O W   ~  version 25.10.2

    Launching `hello-channels.nf` [insane_fermat] DSL2 - revision: 8e62fcbeb1

    executor >  local (3)
    [24/76da2f] sayHello (2) [  0%] 0 of 3 ✘
    Before splitCsv: /workspaces/training/hello-nextflow/data/greetings.csv
    After splitCsv: [Hello, English, 123]
    After splitCsv: [Bonjour, French, 456]
    After splitCsv: [Holà, Spanish, 789]
    ERROR ~ Error executing process > 'sayHello (2)'

    Caused by:
      Missing output file(s) `[Bonjour, French, 456]-output.txt` expected by process `sayHello (2)`


    Command executed:

      echo '[Bonjour, French, 456]' > '[Bonjour, French, 456]-output.txt'

    Command exit status:
      0

    Command output:
      (empty)

    Work dir:
      /workspaces/training/hello-nextflow/work/24/76da2fcc4876b61632749f99e26a50

    Tip: you can try to figure out what's wrong by changing to the process work dir and showing the script file named `.command.sh`

    -- Check '.nextflow.log' file for details
    ```

흥미롭게도 이것도 실패하지만 다른 오류가 발생합니다.
이번에 Nextflow는 파일 내용을 구문 분석했지만(야호!) 각 행을 배열로 로드했으며 각 배열이 channel의 요소입니다.

<figure class="excalidraw">
--8<-- "docs/en/docs/hello_nextflow/img/hello-channels-split-fail.svg"
</figure>

각 행에서 첫 번째 열만 가져오도록 지시해야 합니다.
그러면 이것을 어떻게 풀까요?

이전에 channel의 내용을 풀기 위해 `flatten()`을 사용했지만, flatten은 *모든 것*을 푸는 방식이므로 여기서는 작동하지 않습니다(직접 시도해 보고 싶다면 해 보십시오).

대신 Nextflow 파이프라인에서 정말 유용하고 자주 등장하는 `map()`이라는 다른 연산자를 사용할 것입니다.

### 4.3. `map()` 연산자를 사용하여 인사말 추출

[`map()`](https://nextflow.io/docs/latest/reference/operator.html#map) 연산자는 channel의 내용에 모든 종류의 매핑을 수행할 수 있게 해주는 매우 편리한 도구입니다.

이 경우 데이터 파일의 각 행에서 원하는 요소를 추출하는 데 사용할 것입니다.
구문은 다음과 같습니다:

```groovy title="구문"
.map { row -> row[0] }
```

이것은 'channel의 각 행에 대해 포함된 0번째(첫 번째) 항목을 가져오라'는 의미입니다.

이것을 CSV 구문 분석에 적용해 봅시다.

#### 4.3.1. Channel에 `map()` 적용

workflow 블록에서 다음 코드 변경을 수행하십시오:

=== "수정 후"

    ```groovy title="hello-channels.nf" linenums="27" hl_lines="9 10"
    workflow {

        main:
        // CSV 파일에서 입력용 channel 생성
        greeting_ch = channel.fromPath(params.input)
                             .view { csv -> "Before splitCsv: $csv" }
                             .splitCsv()
                             .view { csv -> "After splitCsv: $csv" }
                             .map { item -> item[0] }
                             .view { csv -> "After map: $csv" }
        // 인사말을 내보냅니다
        sayHello(greeting_ch)

        publish:
        first_output = sayHello.out
    }
    ```

=== "수정 전"

    ```groovy title="hello-channels.nf" linenums="27"
    workflow {

        main:
        // CSV 파일에서 입력용 channel 생성
        greeting_ch = channel.fromPath(params.input)
                             .view { csv -> "Before splitCsv: $csv" }
                             .splitCsv()
                             .view { csv -> "After splitCsv: $csv" }
        // 인사말을 내보냅니다
        sayHello(greeting_ch)

        publish:
        first_output = sayHello.out
    }
    ```

연산자가 예상대로 작동하는지 확인하기 위해 또 다른 `view()` 호출을 추가했습니다.

#### 4.3.2. Workflow 실행

한 번 더 실행해 봅시다:

```bash
nextflow run hello-channels.nf
```

??? success "명령 출력"

    ```console
     N E X T F L O W   ~  version 25.10.2

    Launching `hello-channels.nf` [focused_volhard] DSL2 - revision: de435e45be

    executor >  local (3)
    [54/6eebe3] sayHello (3) [100%] 3 of 3 ✔
    Before splitCsv: /workspaces/training/hello-nextflow/data/greetings.csv
    After splitCsv: [Hello, English, 123]
    After splitCsv: [Bonjour, French, 456]
    After splitCsv: [Holà, Spanish, 789]
    After map: Hello
    After map: Bonjour
    After map: Holà
    ```

이번에는 오류 없이 실행될 것입니다.

`view()` 문의 출력을 보면 다음을 볼 수 있습니다:

- 단일 `Before splitCsv:` 문: 해당 시점에서 channel에는 원래 파일 경로인 하나의 항목이 포함됩니다.
- 세 개의 별도 `After splitCsv:` 문: 각 인사말에 대해 하나씩이지만 각각 파일의 해당 줄에 해당하는 배열 내에 포함되어 있습니다.
- 세 개의 별도 `After map:` 문: 각 인사말에 대해 하나씩이며 이제 channel의 개별 요소입니다.

_출력에서 줄이 다른 순서로 나타날 수 있습니다._

<figure class="excalidraw">
--8<-- "docs/en/docs/hello_nextflow/img/hello-channels-split-and-map.svg"
</figure>

출력 파일을 확인하여 각 인사말이 올바르게 추출되어 workflow를 통해 처리되었는지 확인할 수도 있습니다.

이전과 동일한 결과를 얻었지만, 이제 코드를 수정하지 않고 입력 파일을 수정하여 처리하려는 인사말 channel에 더 많은 요소를 추가할 수 있는 훨씬 더 많은 유연성을 갖게 되었습니다.
나중 교육에서 복잡한 입력을 처리하는 더 정교한 접근 방식을 배우게 됩니다.

### 핵심 정리

`.fromPath()` channel 생성자와 `splitCsv()` 및 `map()` 연산자를 사용하여 입력 값 파일을 읽고 적절하게 처리하는 방법을 알게 되었습니다.

더 일반적으로, Nextflow가 process에 대한 입력을 관리하기 위해 **channel**을 사용하고 내용을 변환하기 위해 **연산자**를 사용하는 방법에 대한 기본적인 이해를 갖게 되었습니다.
또한 channel이 암묵적으로 병렬 실행을 처리하는 방법도 보았습니다.

<figure class="excalidraw">
--8<-- "docs/en/docs/hello_nextflow/img/hello-channels-parallel.svg"
</figure>

### 다음 단계

충분히 휴식을 취하십시오. 이 파트에서 열심히 작업했습니다!

준비가 되면 [**파트 3: Hello Workflow**](./03_hello_workflow.md)로 이동하여 더 많은 단계를 추가하고 적절한 workflow로 연결하는 방법을 배우십시오.

---

## 퀴즈

<quiz>
Nextflow에서 channel이란 무엇입니까?
- [ ] 파일 경로 지정
- [ ] Process 정의
- [x] Process 간 데이터 전달을 위한 큐와 같은 구조
- [ ] 구성 설정

자세히 알아보기: [1.1. 입력 channel 생성](#11-입력-channel-생성)
</quiz>

<quiz>
이 코드는 무엇을 출력합니까?

```groovy
channel.of('Hello', 'Bonjour', 'Hola')
    .view()
```

- [ ] `['Hello', 'Bonjour', 'Hola']` (단일 리스트)
- [x] 각 요소가 별도의 줄에: `Hello`, `Bonjour`, `Hola`
- [ ] 아무것도 (channel은 기본적으로 출력하지 않음)
- [ ] 오류 (잘못된 구문)

자세히 알아보기: [1.1. 입력 channel 생성](#11-입력-channel-생성)
</quiz>

<quiz>
Channel에 여러 값이 포함되어 있을 때 Nextflow는 process 실행을 어떻게 처리합니까?
- [ ] Process가 모든 값으로 한 번 실행됨
- [x] Process가 channel의 각 값에 대해 한 번씩 실행됨
- [ ] Process가 첫 번째 값으로만 실행됨
- [ ] Process가 마지막 값으로만 실행됨

자세히 알아보기: [2. Workflow를 수정하여 여러 입력 값에서 실행](#2-workflow를-수정하여-여러-입력-값에서-실행)
</quiz>

<quiz>
`flatten()` 연산자는 무엇을 합니까?
- [ ] 여러 channel을 하나로 결합
- [ ] Channel 요소 정렬
- [x] 배열을 개별 요소로 풀기
- [ ] 중복 요소 제거

자세히 알아보기: [3.2.1. `flatten()` 연산자 추가](#321-flatten-연산자-추가)
</quiz>

<quiz>
`view()` 연산자의 목적은 무엇입니까?
- [ ] Channel 내용 필터링
- [ ] Channel 요소 변환
- [x] Channel 내용 검사 및 디버깅
- [ ] Channel 내용을 파일에 저장

자세히 알아보기: [1.4. `view()`를 사용하여 channel 내용 검사](#14-view를-사용하여-channel-내용-검사)
</quiz>

<quiz>
`splitCsv()`는 무엇을 합니까?
- [ ] Channel 내용에서 CSV 파일 생성
- [ ] 문자열을 쉼표로 분할
- [x] CSV 파일을 각 행을 나타내는 배열로 구문 분석
- [ ] 여러 CSV 파일 병합

자세히 알아보기: [4.2. `splitCsv()` 연산자를 사용하여 파일 구문 분석](#42-splitcsv-연산자를-사용하여-파일-구문-분석)
</quiz>

<quiz>
`map()` 연산자의 목적은 무엇입니까?
- [ ] Channel에서 요소 필터링
- [ ] 여러 channel 결합
- [x] Channel의 각 요소 변환
- [ ] Channel의 요소 수 세기

자세히 알아보기: [4.3. `map()` 연산자를 사용하여 인사말 추출](#43-map-연산자를-사용하여-인사말-추출)
</quiz>

<quiz>
여러 입력을 처리할 때 동적 출력 파일 이름을 사용하는 것이 중요한 이유는 무엇입니까?
- [ ] 성능 향상
- [ ] 디스크 공간 절약
- [x] 출력 파일이 서로 덮어쓰이는 것을 방지
- [ ] 재개 기능 활성화

자세히 알아보기: [2.2. 출력 파일 이름이 고유하도록 보장](#22-출력-파일-이름이-고유하도록-보장)
</quiz>
