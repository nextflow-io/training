# 파트 3: Hello Workflow

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } AI 지원 번역 - [자세히 알아보기 및 개선 제안](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

<div class="video-wrapper">
  <iframe width="560" height="315" src="https://www.youtube.com/embed/_aO56V3iXGI?si=Irl9nAQniDyICp2b&amp;list=PLPZ8WHdZGxmWKozQuzr27jyMGqp9kElVK&amp;cc_load_policy=1&amp;cc_lang_pref=ko" title="YouTube video player" frameborder="0" allow="accelerometer; autoplay; clipboard-write; encrypted-media; gyroscope; picture-in-picture; web-share" referrerpolicy="strict-origin-when-cross-origin" allowfullscreen></iframe>
</div>

/// caption
:fontawesome-brands-youtube:{ .youtube } Nextflow YouTube 채널에서 [전체 재생목록](https://youtube.com/playlist?list=PLPZ8WHdZGxmWKozQuzr27jyMGqp9kElVK&si=eF7cLR62goy-lc6n)을 확인하십시오.

:green_book: 비디오 스크립트는 [여기](./transcripts/03_hello_workflow.md)에서 확인하실 수 있습니다.
///

대부분의 실제 워크플로우는 둘 이상의 단계를 포함합니다.
이 교육 모듈에서는 다단계 워크플로우에서 프로세스를 함께 연결하는 방법을 학습합니다.

다음과 같은 Nextflow 방식을 다룹니다:

1. 한 프로세스에서 다음 프로세스로 데이터 흐르게 하기
2. 여러 프로세스 호출의 출력을 단일 프로세스 호출로 수집하기
3. 프로세스에 추가 매개변수 전달하기
4. 프로세스에서 나오는 여러 출력 처리하기

파트 1과 2의 도메인에 구애받지 않는 Hello World 예제를 계속 구축할 것입니다.
이번에는 사람들이 실제 워크플로우를 구축하는 방식을 더 잘 반영하기 위해 다음과 같은 변경을 수행합니다:

1. 인사말을 대문자로 변환하는 두 번째 단계 추가
2. 모든 변환된 인사말을 수집하여 단일 파일에 기록하는 세 번째 단계 추가
3. 최종 출력 파일의 이름을 지정하는 매개변수를 추가하고 이를 수집 단계에 보조 입력으로 전달
4. 수집 단계에서 처리된 내용에 대한 간단한 통계도 보고하도록 설정

??? info "이 섹션부터 시작하는 방법"

    이 섹션은 [Hello Nextflow](./index.md) 과정의 파트 1-2를 완료했다고 가정하지만, 해당 섹션에서 다룬 기본 사항에 익숙하다면 특별한 준비 없이 여기서 시작할 수 있습니다.

---

## 0. 준비 운동: `hello-workflow.nf` 실행

시작점으로 워크플로우 스크립트 `hello-workflow.nf`를 사용할 것입니다.
이 스크립트는 이 교육 과정의 파트 2를 완료하여 생성된 스크립트와 동일하지만, `view()` 문을 제거하고 출력 대상을 변경했습니다:

```groovy title="hello-workflow.nf" linenums="37" hl_lines="3"
output {
    first_output {
        path 'hello_workflow'
        mode 'copy'
    }
}
```

이 다이어그램은 워크플로우의 현재 작동 방식을 요약합니다.
익숙해 보일 것입니다. 단, 이제 프로세스의 출력이 입력과 마찬가지로 채널로 패키징되는 것을 명시적으로 보여줍니다.
곧 이 출력 채널을 잘 활용할 것입니다.

<figure class="excalidraw">
--8<-- "docs/en/docs/hello_nextflow/img/hello-workflow-channels.svg"
</figure>

변경을 시작하기 전에 모든 것이 제대로 작동하는지 확인하기 위해 스크립트를 한 번 실행하십시오:

```bash
nextflow run hello-workflow.nf
```

??? success "명령 출력"

    ```console
     N E X T F L O W   ~  version 25.10.2

    Launching `hello-workflow.nf` [admiring_lamarr] DSL2 - revision: 4d4053520d

    executor >  local (3)
    [b1/5826b5] process > sayHello (2) [100%] 3 of 3 ✔
    ```

이전과 마찬가지로 `output` 블록에 지정된 위치에서 출력 파일을 찾을 수 있습니다.
이 장에서는 `results/hello_workflow/` 아래에 있습니다.

??? abstract "디렉토리 내용"

    ```console
    results/hello_workflow
    ├── Bonjour-output.txt
    ├── Hello-output.txt
    └── Holà-output.txt
    ```

이것이 정상적으로 작동했다면 다단계 워크플로우를 조립하는 방법을 학습할 준비가 되었습니다.

---

## 1. 워크플로우에 두 번째 단계 추가

각 인사말을 대문자로 변환하는 단계를 추가할 것입니다.

<figure class="excalidraw">
--8<-- "docs/en/docs/hello_nextflow/img/hello-multistep.svg"
</figure>

이를 위해 세 가지를 수행해야 합니다:

- 대문자 변환에 사용할 명령 정의
- 대문자 변환 명령을 감싸는 새 프로세스 작성
- 워크플로우 블록에서 새 프로세스를 호출하고 `sayHello()` 프로세스의 출력을 입력으로 받도록 설정

### 1.1. 대문자 변환 명령 정의 및 터미널에서 테스트

인사말을 대문자로 변환하기 위해 다음 구문과 함께 '텍스트 대체'를 위한 `tr`이라는 고전적인 UNIX 도구를 사용할 것입니다:

```bash title="구문"
tr '[a-z]' '[A-Z]'
```

이것은 악센트가 있는 문자를 고려하지 않는 매우 단순한 텍스트 대체 한 줄 명령이므로 예를 들어 'Holà'는 'HOLà'가 됩니다. 하지만 Nextflow 개념을 시연하기에는 충분히 작동하며 그것이 중요합니다.

테스트하려면 `echo 'Hello World'` 명령을 실행하고 그 출력을 `tr` 명령으로 파이프할 수 있습니다:

```bash
echo 'Hello World' | tr '[a-z]' '[A-Z]' > UPPER-output.txt
```

출력은 `Hello World` 문자열의 대문자 버전이 포함된 `UPPER-output.txt`라는 텍스트 파일입니다.

??? abstract "파일 내용"

    ```console title="UPPER-output.txt"
    HELLO WORLD
    ```

이것이 기본적으로 워크플로우로 수행하려는 것입니다.

### 1.2. 대문자 변환 단계를 Nextflow 프로세스로 작성

첫 번째 프로세스를 모델로 새 프로세스를 만들 수 있습니다. 동일한 모든 구성 요소를 사용하려고 하기 때문입니다.

다음 프로세스 정의를 첫 번째 바로 아래에 워크플로우 스크립트에 추가하십시오:

```groovy title="hello-workflow.nf" linenums="20"
/*
 * 텍스트 대체 도구를 사용하여 인사말을 대문자로 변환
 */
process convertToUpper {

    input:
    path input_file

    output:
    path "UPPER-${input_file}"

    script:
    """
    cat '${input_file}' | tr '[a-z]' '[A-Z]' > 'UPPER-${input_file}'
    """
}
```

이 프로세스에서는 첫 번째 프로세스의 출력에 대해 원래 수행했던 것과 유사하게 입력 파일 이름을 기반으로 두 번째 출력 파일 이름을 구성합니다.

### 1.3. 워크플로우 블록에 새 프로세스 호출 추가

이제 방금 정의한 프로세스를 실제로 호출하도록 Nextflow에 지시해야 합니다.

워크플로우 블록에서 다음 코드 변경을 수행하십시오:

=== "수정 후"

    ```groovy title="hello-workflow.nf" linenums="44" hl_lines="10-11"
    workflow {

        main:
        // CSV 파일에서 입력용 채널 생성
        greeting_ch = channel.fromPath(params.input)
                            .splitCsv()
                            .map { line -> line[0] }
        // 인사말을 내보냅니다
        sayHello(greeting_ch)
        // 인사말을 대문자로 변환
        convertToUpper()

        publish:
        first_output = sayHello.out
    }
    ```

=== "수정 전"

    ```groovy title="hello-workflow.nf" linenums="44"
    workflow {

        main:
        // CSV 파일에서 입력용 채널 생성
        greeting_ch = channel.fromPath(params.input)
                            .splitCsv()
                            .map { line -> line[0] }
        // 인사말을 내보냅니다
        sayHello(greeting_ch)

        publish:
        first_output = sayHello.out
    }
    ```

`convertToUpper()` 프로세스에 무엇을 입력해야 하는지 지정하지 않았기 때문에 아직 작동하지 않습니다.

### 1.4. 첫 번째 프로세스의 출력을 두 번째 프로세스로 전달

이제 `sayHello()` 프로세스의 출력이 `convertToUpper()` 프로세스로 흐르도록 해야 합니다.

편리하게도 준비 운동 섹션의 다이어그램에 표시된 것처럼 Nextflow는 자동으로 프로세스의 출력을 채널로 패키징합니다.
프로세스의 출력 채널을 `<process>.out`으로 참조할 수 있습니다.

따라서 `sayHello` 프로세스의 출력은 `sayHello.out`이라는 채널이며, 이를 `convertToUpper()` 호출에 바로 연결할 수 있습니다.

<figure class="excalidraw">
--8<-- "docs/en/docs/hello_nextflow/img/hello-multistep-connector.svg"
</figure>

워크플로우 블록에서 다음 코드 변경을 수행하십시오:

=== "수정 후"

    ```groovy title="hello-workflow.nf" linenums="53" hl_lines="2"
        // 인사말을 대문자로 변환
        convertToUpper(sayHello.out)
    ```

=== "수정 전"

    ```groovy title="hello-workflow.nf" linenums="53" hl_lines="2"
        // 인사말을 대문자로 변환
        convertToUpper()
    ```

이와 같이 간단한 경우(하나의 출력에서 하나의 입력으로)에는 두 프로세스를 연결하기 위해 이것만 하면 됩니다!

### 1.5. 워크플로우 출력 게시 설정

마지막으로, 두 번째 프로세스의 결과도 게시하도록 워크플로우 출력을 업데이트합시다.

#### 1.5.1. `workflow` 블록의 `publish:` 섹션 업데이트

`workflow` 블록에서 다음 코드 변경을 수행하십시오:

=== "수정 후"

    ```groovy title="hello-workflow.nf" linenums="56" hl_lines="3"
        publish:
        first_output = sayHello.out
        uppercased = convertToUpper.out
    }
    ```

=== "수정 전"

    ```groovy title="hello-workflow.nf" linenums="56"
        publish:
        first_output = sayHello.out
    }
    ```

로직은 이전과 동일합니다.

#### 1.5.2. `output` 블록 업데이트

`output` 블록에서 다음 코드 변경을 수행하십시오:

=== "수정 후"

    ```groovy title="hello-workflow.nf" linenums="61" hl_lines="6-9"
    output {
        first_output {
            path 'hello_workflow'
            mode 'copy'
        }
        uppercased {
            path 'hello_workflow'
            mode 'copy'
        }
    }
    ```

=== "수정 전"

    ```groovy title="hello-workflow.nf" linenums="61"
    output {
        first_output {
            path 'hello_workflow'
            mode 'copy'
        }
    }
    ```

다시 한 번, 로직은 이전과 동일합니다.

이것은 모든 개별 출력에 대해 매우 세분화된 수준에서 출력 설정을 제어할 수 있음을 보여줍니다.
프로세스 중 하나에 대해 경로나 게시 모드를 변경하여 어떤 일이 발생하는지 자유롭게 시도해 보십시오.

물론, 이것은 여기서 일부 정보를 반복하고 있다는 것을 의미하며, 모든 출력의 위치를 같은 방식으로 업데이트하려면 불편해질 수 있습니다.
과정의 후반부에서 구조화된 방식으로 여러 출력에 대한 이러한 설정을 구성하는 방법을 학습하게 됩니다.

### 1.6. `-resume`으로 워크플로우 실행

워크플로우의 첫 번째 단계를 이미 성공적으로 실행했으므로 `-resume` 플래그를 사용하여 테스트해 봅시다.

```bash
nextflow run hello-workflow.nf -resume
```

??? success "명령 출력"

    ```console
     N E X T F L O W   ~  version 25.10.2

    Launching `hello-workflow.nf` [high_cantor] DSL2 - revision: d746983511

    executor >  local (3)
    [ab/816321] process > sayHello (3)       [100%] 3 of 3, cached: 3 ✔
    [e0/ecf81b] process > convertToUpper (3) [100%] 3 of 3 ✔
    ```

방금 추가한 새 프로세스에 해당하는 추가 줄이 콘솔 출력에 나타납니다.

`output` 블록에 설정된 대로 `results/hello_workflow` 디렉토리에서 출력을 찾을 수 있습니다.

??? abstract "디렉토리 내용"

    ```console
    results/hello_workflow/
    ├── Bonjour-output.txt
    ├── Hello-output.txt
    ├── Holà-output.txt
    ├── UPPER-Bonjour-output.txt
    ├── UPPER-Hello-output.txt
    └── UPPER-Holà-output.txt
    ```

편리합니다! 하지만 두 번째 프로세스 호출 중 하나의 작업 디렉토리를 살펴볼 가치가 있습니다.

??? abstract "디렉토리 내용"

    ```console
    work/e0/ecf81b4cacc648b9b994218d5b29d7/
    ├── Holà-output.txt -> /workspaces/training/hello-nextflow/work/ab/81632178cd37e9e815959278808819/Holà-output.txt
    └── UPPER-Holà-output.txt
    ```

두 개의 `*-output` 파일이 있습니다: 첫 번째 프로세스의 출력과 두 번째의 출력입니다.

첫 번째 프로세스의 출력이 거기에 있는 이유는 Nextflow가 실행에 필요한 모든 것을 동일한 하위 디렉토리 내에 갖추기 위해 거기에 **스테이징**했기 때문입니다.

그러나 실제로는 첫 번째 프로세스 호출의 하위 디렉토리에 있는 원본 파일을 가리키는 심볼릭 링크입니다.
기본적으로, 여기서처럼 단일 머신에서 실행할 때 Nextflow는 입력 및 중간 파일을 스테이징하기 위해 복사본 대신 심볼릭 링크를 사용합니다.

이제, 계속 진행하기 전에, `sayHello`의 출력을 `convertToUpper`의 입력에 연결하기만 하면 두 프로세스가 직렬로 실행될 수 있었다는 점을 생각해 보십시오.
Nextflow가 개별 입력 및 출력 파일을 처리하고 두 명령 사이에 전달하는 어려운 작업을 대신해 주었습니다.

이것이 Nextflow 채널이 매우 강력한 이유 중 하나입니다: 워크플로우 단계를 함께 연결하는 데 관련된 반복 작업을 처리해 줍니다.

### 핵심 정리

한 단계의 출력을 다음 단계의 입력으로 제공하여 프로세스를 함께 연결하는 방법을 알게 되었습니다.

### 다음 단계

일괄 프로세스 호출에서 출력을 수집하고 단일 프로세스에 공급하는 방법을 학습합니다.

---

## 2. 모든 인사말을 수집하는 세 번째 단계 추가

여기서 여러 인사말에 대해 수행하는 것처럼 채널의 각 요소에 변환을 적용하기 위해 프로세스를 사용할 때, 때때로 해당 프로세스의 출력 채널에서 요소를 수집하여 일종의 분석 또는 요약을 수행하는 다른 프로세스에 공급하고 싶을 수 있습니다.

시연을 위해, `convertToUpper` 프로세스에서 생성된 모든 대문자 인사말을 수집하여 단일 파일에 기록하는 새 단계를 파이프라인에 추가할 것입니다.

<figure class="excalidraw">
--8<-- "docs/en/docs/hello_nextflow/img/hello-collect.svg"
</figure>

미리 스포일러를 하자면, 이것은 매우 유용한 연산자를 포함할 것입니다.

### 2.1. 수집 명령 정의 및 터미널에서 테스트

워크플로우에 추가하려는 수집 단계는 `cat` 명령을 사용하여 여러 대문자 인사말을 단일 파일로 연결합니다.

이전에 했던 것처럼 터미널에서 명령 자체를 실행하여 예상대로 작동하는지 확인해 봅시다.

터미널에서 다음을 실행하십시오:

```bash
echo 'Hello' | tr '[a-z]' '[A-Z]' > UPPER-Hello-output.txt
echo 'Bonjour' | tr '[a-z]' '[A-Z]' > UPPER-Bonjour-output.txt
echo 'Holà' | tr '[a-z]' '[A-Z]' > UPPER-Holà-output.txt
cat UPPER-Hello-output.txt UPPER-Bonjour-output.txt UPPER-Holà-output.txt > COLLECTED-output.txt
```

출력은 원래 인사말의 대문자 버전이 포함된 `COLLECTED-output.txt`라는 텍스트 파일입니다.

??? abstract "파일 내용"

    ```console title="COLLECTED-output.txt"
    HELLO
    BONJOUR
    HOLà
    ```

이것이 워크플로우로 달성하려는 결과입니다.

### 2.2. 수집 단계를 수행하는 새 프로세스 생성

새 프로세스를 생성하고 `collectGreetings()`라고 부릅시다.
지금까지 본 것을 기반으로 작성을 시작할 수 있습니다.

#### 2.2.1. 프로세스의 '명백한' 부분 작성

다음 프로세스 정의를 워크플로우 스크립트에 추가하십시오:

```groovy title="hello-workflow.nf" linenums="37"
/*
 * 대문자 인사말을 하나의 출력 파일에 수집
 */
process collectGreetings {

    input:
    ???

    output:
    path "COLLECTED-output.txt"

    script:
    """
    cat ??? > 'COLLECTED-output.txt'
    """
}
```

이것은 지금까지 배운 것을 기반으로 자신 있게 작성할 수 있는 부분입니다.
그러나 이것은 작동하지 않습니다!
입력 정의와 스크립트 명령의 첫 번째 절반을 생략했는데, 이것을 어떻게 작성해야 하는지 알아내야 합니다.

#### 2.2.2. `collectGreetings()`에 대한 입력 정의

`convertToUpper()` 프로세스에 대한 모든 호출에서 인사말을 수집해야 합니다.
워크플로우의 이전 단계에서 얻을 수 있는 것이 무엇인지 알고 있습니까?

`convertToUpper()`에서 출력되는 채널에는 대문자 인사말이 포함된 개별 파일의 경로가 포함됩니다.
이것은 하나의 입력 슬롯에 해당합니다; 단순함을 위해 `input_files`라고 부르겠습니다.

프로세스 블록에서 다음 코드 변경을 수행하십시오:

=== "수정 후"

    ```groovy title="hello-workflow.nf" linenums="42" hl_lines="2"
          input:
          path input_files
    ```

=== "수정 전"

    ```groovy title="hello-workflow.nf" linenums="42" hl_lines="2"
          input:
          ???
    ```

여러 파일이 포함될 것으로 예상하지만 `path` 접두사를 사용합니다.

#### 2.2.3. 연결 명령 구성

여기서 약간 까다로워질 수 있는데, 임의의 수의 입력 파일을 처리할 수 있어야 하기 때문입니다.
구체적으로, 명령을 미리 작성할 수 없으므로, 런타임에 프로세스로 흘러 들어오는 입력을 기반으로 명령을 구성하는 방법을 Nextflow에 알려줘야 합니다.

다시 말해, `[file1.txt, file2.txt, file3.txt]` 요소가 포함된 입력 채널이 있으면 Nextflow가 이를 `cat file1.txt file2.txt file3.txt`로 변환해야 합니다.

다행히도, 스크립트 명령에 `cat ${input_files}`라고 간단히 작성하면 Nextflow가 기꺼이 이를 수행합니다.

프로세스 블록에서 다음 코드 변경을 수행하십시오:

=== "수정 후"

    ```groovy title="hello-workflow.nf" linenums="54" hl_lines="3"
        script:
        """
        cat ${input_files} > 'COLLECTED-output.txt'
        """
    ```

=== "수정 전"

    ```groovy title="hello-workflow.nf" linenums="54"
        script:
        """
        cat ??? > 'COLLECTED-output.txt'
        """
    ```

이론적으로 이것은 임의의 수의 입력 파일을 처리해야 합니다.

!!! tip "팁"

    일부 명령줄 도구는 각 입력 파일에 대해 인수(예: `-input`)를 제공해야 합니다.
    그런 경우 명령을 구성하기 위해 약간의 추가 작업을 수행해야 합니다.
    [Nextflow for Genomics](../../nf4_science/genomics/) 교육 과정에서 이에 대한 예를 볼 수 있습니다.

### 2.3. 워크플로우에 수집 단계 추가

이제 대문자 변환 단계의 출력에서 수집 프로세스를 호출하기만 하면 됩니다.
그것도 `convertToUpper.out`이라는 채널입니다.

<figure class="excalidraw">
--8<-- "docs/en/docs/hello_nextflow/img/hello-collect-connector.svg"
</figure>

#### 2.3.1. 프로세스 호출 연결

워크플로우 블록에서 다음 코드 변경을 수행하십시오:

=== "수정 후"

    ```groovy title="hello-workflow.nf" linenums="75" hl_lines="4 5"
        // 인사말을 대문자로 변환
        convertToUpper(sayHello.out)

        // 모든 인사말을 하나의 파일에 수집
        collectGreetings(convertToUpper.out)
    }
    ```

=== "수정 전"

    ```groovy title="hello-workflow.nf" linenums="75"
        // 인사말을 대문자로 변환
        convertToUpper(sayHello.out)
    }
    ```

이것은 `convertToUpper()`의 출력을 `collectGreetings()`의 입력에 연결합니다.

#### 2.3.2. `-resume`으로 워크플로우 실행

시도해 봅시다.

```bash
nextflow run hello-workflow.nf -resume
```

??? success "명령 출력"

    ```console hl_lines="8"
    N E X T F L O W   ~  version 25.10.2

    Launching `hello-workflow.nf` [mad_gilbert] DSL2 - revision: 6acfd5e28d

    executor >  local (3)
    [79/33b2f0] sayHello (2)         | 3 of 3, cached: 3 ✔
    [99/79394f] convertToUpper (3)   | 3 of 3, cached: 3 ✔
    [47/50fe4a] collectGreetings (1) | 3 of 3 ✔
    ```

세 번째 단계를 포함하여 성공적으로 실행됩니다.

그러나 마지막 줄에서 `collectGreetings()`에 대한 호출 수를 보십시오.
하나만 예상했지만 세 개가 있습니다.

이제 최종 출력 파일의 내용을 살펴보십시오.

??? abstract "파일 내용"

    ```console title="results/COLLECTED-output.txt"
    Holà
    ```

이런. 수집 단계가 각 인사말에 대해 개별적으로 실행되었으며, 이것은 우리가 원했던 것이 아닙니다.

<figure class="excalidraw">
--8<-- "docs/en/docs/hello_nextflow/img/hello-collect-no-operator.svg"
</figure>

세 번째 단계가 `convertToUpper()`에서 출력된 채널의 모든 요소에서 실행되기를 원한다고 Nextflow에 명시적으로 알려야 합니다.

### 2.4. 연산자를 사용하여 인사말을 단일 입력으로 수집

네, 다시 한번 우리 문제에 대한 답은 연산자입니다.

구체적으로, 적절한 이름의 [`collect()`](https://nextflow.io/docs/latest/reference/operator.html#collect) 연산자를 사용할 것입니다.

#### 2.4.1. `collect()` 연산자 추가

이번에는 채널 팩토리의 컨텍스트에서 연산자를 추가하는 것이 아니라 출력 채널에 추가하기 때문에 약간 다르게 보일 것입니다.

`convertToUpper.out`을 가져와서 `collect()` 연산자를 추가하면 `convertToUpper.out.collect()`가 됩니다.
이를 `collectGreetings()` 프로세스 호출에 직접 연결할 수 있습니다.

워크플로우 블록에서 다음 코드 변경을 수행하십시오:

=== "수정 후"

    ```groovy title="hello-workflow.nf" linenums="73" hl_lines="2"
        // 모든 인사말을 하나의 파일에 수집
        collectGreetings(convertToUpper.out.collect())
    }
    ```

=== "수정 전"

    ```groovy title="hello-workflow.nf" linenums="73" hl_lines="2"
        // 모든 인사말을 하나의 파일에 수집
        collectGreetings(convertToUpper.out)
    }
    ```

#### 2.4.2. `view()` 문 추가

채널 내용의 전후 상태를 시각화하기 위해 몇 가지 `view()` 문도 포함합시다.

=== "수정 후"

    ```groovy title="hello-workflow.nf" linenums="73" hl_lines="4-6"
        // 모든 인사말을 하나의 파일에 수집
        collectGreetings(convertToUpper.out.collect())

        // 선택적 view 문
        convertToUpper.out.view { contents -> "Before collect: $contents" }
        convertToUpper.out.collect().view { contents -> "After collect: $contents" }
    }
    ```

=== "수정 전"

    ```groovy title="hello-workflow.nf" linenums="73"
        // 모든 인사말을 하나의 파일에 수집
        collectGreetings(convertToUpper.out.collect())
    }
    ```

`view()` 문은 원하는 곳에 배치할 수 있습니다; 가독성을 위해 호출 바로 뒤에 배치했습니다.

#### 2.4.3. `-resume`으로 워크플로우 다시 실행

시도해 봅시다:

```bash
nextflow run hello-workflow.nf -resume
```

??? success "명령 출력"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `hello-workflow.nf` [soggy_franklin] DSL2 - revision: bc8e1b2726

    [d6/cdf466] sayHello (1)       | 3 of 3, cached: 3 ✔
    [99/79394f] convertToUpper (2) | 3 of 3, cached: 3 ✔
    [1e/83586c] collectGreetings   | 1 of 1 ✔
    Before collect: /workspaces/training/hello-nextflow/work/b3/d52708edba8b864024589285cb3445/UPPER-Bonjour-output.txt
    Before collect: /workspaces/training/hello-nextflow/work/99/79394f549e3040dfc2440f69ede1fc/UPPER-Hello-output.txt
    Before collect: /workspaces/training/hello-nextflow/work/aa/56bfe7cf00239dc5badc1d04b60ac4/UPPER-Holà-output.txt
    After collect: [/workspaces/training/hello-nextflow/work/b3/d52708edba8b864024589285cb3445/UPPER-Bonjour-output.txt, /workspaces/training/hello-nextflow/work/99/79394f549e3040dfc2440f69ede1fc/UPPER-Hello-output.txt, /workspaces/training/hello-nextflow/work/aa/56bfe7cf00239dc5badc1d04b60ac4/UPPER-Holà-output.txt]
    ```

성공적으로 실행되지만, 로그 출력이 이것보다 약간 지저분하게 보일 수 있습니다(가독성을 위해 정리했습니다).

이번에는 세 번째 단계가 한 번만 호출되었습니다!
`view()` 문의 출력을 보면 다음을 볼 수 있습니다:

- 세 개의 `Before collect:` 문, 각 인사말에 대해 하나씩: 해당 시점에서 파일 경로는 채널의 개별 항목입니다.
- 단일 `After collect:` 문: 세 개의 파일 경로가 이제 단일 요소로 패키징되었습니다.

다음 다이어그램으로 요약할 수 있습니다:

<figure class="excalidraw">
--8<-- "docs/en/docs/hello_nextflow/img/hello-collect-WITH-operator.svg"
</figure>

마지막으로, 출력 파일의 내용을 살펴보고 모든 것이 올바르게 작동했는지 확인할 수 있습니다.

??? abstract "파일 내용"

    ```console title="results/COLLECTED-output.txt"
    BONJOUR
    HELLO
    HOLà
    ```

이번에는 최종 출력 파일에 세 개의 인사말이 모두 있습니다. 성공입니다!

!!! note "참고"

    `-resume` 없이 여러 번 실행하면 인사말의 순서가 실행마다 변경되는 것을 볼 수 있습니다.
    이것은 요소가 프로세스 호출을 통해 흐르는 순서가 일관되게 보장되지 않음을 보여줍니다.

#### 2.4.4. 가독성을 위해 `view()` 문 제거

다음 섹션으로 넘어가기 전에 콘솔 출력을 어지럽히지 않도록 `view()` 문을 삭제하는 것이 좋습니다.

=== "수정 후"

    ```groovy title="hello-workflow.nf" linenums="73"
        // 모든 인사말을 하나의 파일에 수집
        collectGreetings(convertToUpper.out.collect())
    ```

=== "수정 전"

    ```groovy title="hello-workflow.nf" linenums="73" hl_lines="4-6"
        // 모든 인사말을 하나의 파일에 수집
        collectGreetings(convertToUpper.out.collect())

        // 선택적 view 문
        convertToUpper.out.view { contents -> "Before collect: $contents" }
        convertToUpper.out.collect().view { contents -> "After collect: $contents" }
    ```

이것은 기본적으로 2.4.2 지점의 역순 작업입니다.

### 핵심 정리

일괄 프로세스 호출에서 출력을 수집하고 공동 분석 또는 요약 단계에 공급하는 방법을 알게 되었습니다.

요약하면, 지금까지 다음과 같이 구축했습니다:

<figure class="excalidraw">
--8<-- "docs/en/docs/hello_nextflow/img/hello-collect.svg"
</figure>

### 다음 단계

프로세스에 둘 이상의 입력을 전달하는 방법을 학습합니다.

---

## 3. 프로세스에 추가 매개변수 전달

최종 결과를 덮어쓰지 않고 후속 인사말 배치를 처리하기 위해 최종 출력 파일에 특정 이름을 지정할 수 있기를 원합니다.

이를 위해 워크플로우에 다음과 같은 개선을 수행할 것입니다:

- 수집기 프로세스가 출력 파일에 대해 사용자 정의 이름을 허용하도록 수정(`batch_name`)
- 워크플로우에 명령줄 매개변수를 추가하고(`--batch`) 수집기 프로세스에 전달

<figure class="excalidraw">
--8<-- "docs/en/docs/hello_nextflow/img/hello-collect-batch.svg"
</figure>

### 3.1. 수집기 프로세스 수정

추가 입력을 선언하고 출력 파일 이름에 통합해야 합니다.

#### 3.1.1. 추가 입력 선언

좋은 소식: 프로세스 정의에서 원하는 만큼 많은 입력 변수를 선언할 수 있습니다.
이것을 `batch_name`이라고 부르겠습니다.

프로세스 블록에서 다음 코드 변경을 수행하십시오:

=== "수정 후"

    ```groovy title="hello-workflow.nf" linenums="42" hl_lines="3"
        input:
        path input_files
        val batch_name
    ```

=== "수정 전"

    ```groovy title="hello-workflow.nf" linenums="42"
        input:
        path input_files
    ```

프로세스가 원하는 만큼 많은 입력을 예상하도록 설정할 수 있습니다.
현재 이것들은 모두 필수 입력으로 설정되어 있습니다; 워크플로우가 작동하려면 반드시 값을 제공해야 합니다.

Nextflow 여정의 후반부에서 필수 입력과 선택적 입력을 관리하는 방법을 학습하게 됩니다.

#### 3.1.2. 출력 파일 이름에 `batch_name` 변수 사용

이전에 동적 파일 이름을 구성한 것과 같은 방식으로 출력 파일 이름에 변수를 삽입할 수 있습니다.

프로세스 블록에서 다음 코드 변경을 수행하십시오:

=== "수정 후"

    ```groovy title="hello-workflow.nf" linenums="46" hl_lines="2 6"
        output:
        path "COLLECTED-${batch_name}-output.txt"

        script:
        """
        cat ${input_files} > 'COLLECTED-${batch_name}-output.txt'
        """
    ```

=== "수정 전"

    ```groovy title="hello-workflow.nf" linenums="46" hl_lines="2 6"
        output:
        path "COLLECTED-output.txt"

        script:
        """
        cat ${input_files} > 'COLLECTED-output.txt'
        """
    ```

이것은 프로세스가 `batch_name` 값을 사용하여 워크플로우의 최종 출력에 대한 특정 파일 이름을 생성하도록 설정합니다.

### 3.2. `batch` 명령줄 매개변수 추가

이제 `batch_name`에 대한 값을 제공하고 프로세스 호출에 공급하는 방법이 필요합니다.

#### 3.2.1. `params`를 사용하여 매개변수 설정

CLI 매개변수를 선언하기 위해 `params` 시스템을 사용하는 방법을 이미 알고 있습니다.
이를 사용하여 `batch` 매개변수를 선언합시다(게으르기 때문에 기본값과 함께).

파이프라인 매개변수 섹션에서 다음 코드 변경을 수행하십시오:

=== "수정 후"

    ```groovy title="hello-workflow.nf" linenums="55" hl_lines="6"
    /*
     * 파이프라인 매개변수
     */
    params {
        input: Path = 'data/greetings.csv'
        batch: String = 'batch'
    }
    ```

=== "수정 전"

    ```groovy title="hello-workflow.nf" linenums="55"
    /*
     * 파이프라인 매개변수
     */
    params {
        input: Path = 'data/greetings.csv'
    }
    ```

`--input`에 대해 시연한 것처럼 명령줄에서 `--batch`로 값을 지정하여 기본값을 재정의할 수 있습니다.

#### 3.2.2. `batch` 매개변수를 프로세스에 전달

매개변수 값을 프로세스에 제공하려면 프로세스 호출에 추가해야 합니다.

워크플로우 블록에서 다음 코드 변경을 수행하십시오:

=== "수정 후"

    ```groovy title="hello-workflow.nf" linenums="74" hl_lines="2"
        // 모든 인사말을 하나의 파일에 수집
        collectGreetings(convertToUpper.out.collect(), params.batch)
    ```

=== "수정 전"

    ```groovy title="hello-workflow.nf" linenums="74" hl_lines="2"
        // 모든 인사말을 하나의 파일에 수집
        collectGreetings(convertToUpper.out.collect())
    ```

프로세스에 여러 입력을 제공하려면 호출 괄호 안에 쉼표로 구분하여 나열하기만 하면 됩니다.

!!! warning "경고"

    프로세스에 대한 입력은 프로세스의 입력 정의 블록에 나열된 것과 정확히 동일한 순서로 제공해야 합니다.

### 3.3. 워크플로우 실행

명령줄에서 배치 이름으로 실행해 봅시다.

```bash
nextflow run hello-workflow.nf -resume --batch trio
```

??? success "명령 출력"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `hello-workflow.nf` [confident_rutherford] DSL2 - revision: bc58af409c

    executor >  local (1)
    [79/33b2f0] sayHello (2)       | 3 of 3, cached: 3 ✔
    [99/79394f] convertToUpper (2) | 3 of 3, cached: 3 ✔
    [b5/f19efe] collectGreetings   | 1 of 1 ✔
    ```

성공적으로 실행되고 원하는 출력을 생성합니다:

??? abstract "파일 내용"

    ```console title="results/COLLECTED-trio-output.txt"
    HELLO
    BONJOUR
    HOLà
    ```

이제 매개변수를 적절하게 지정하는 한 다른 입력 배치에 대한 후속 실행은 이전 결과를 덮어쓰지 않습니다.

### 핵심 정리

프로세스에 둘 이상의 입력을 전달하는 방법을 알게 되었습니다.

### 다음 단계

여러 출력을 내보내고 편리하게 처리하는 방법을 학습합니다.

---

## 4. 수집기 단계에 출력 추가

지금까지 하나의 출력만 생성하는 프로세스를 사용해 왔습니다.
출력을 다음 프로세스로 전달하는 컨텍스트(예: `convertToUpper(sayHello.out)`)와 `publish:` 섹션의 컨텍스트(예: `first_output = sayHello.out`) 모두에서 `<process>.out` 구문을 사용하여 각각의 출력에 매우 편리하게 액세스할 수 있었습니다.

프로세스가 둘 이상을 생성할 때는 어떻게 됩니까?
여러 출력을 어떻게 처리합니까?
특정 출력을 선택하여 사용할 수 있습니까?

모두 훌륭한 질문이며, 간단한 대답은 그렇습니다!

여러 출력은 별도의 채널로 패키징됩니다.
해당 출력 채널에 이름을 지정하여 나중에 개별적으로 쉽게 참조하거나 인덱스로 참조할 수 있습니다.

시연 목적으로, 주어진 입력 배치에 대해 수집되는 인사말 수를 세고 파일에 보고하고 싶다고 가정합시다.

### 4.1. 인사말 수를 세고 출력하도록 프로세스 수정

이것은 프로세스 정의에 두 가지 주요 변경이 필요합니다: 인사말을 세고 보고 파일을 작성하는 방법이 필요하고, 해당 보고 파일을 프로세스의 `output` 블록에 추가해야 합니다.

#### 4.1.1. 수집된 인사말 수 세기

편리하게도 Nextflow는 프로세스 정의의 `script:` 블록에 임의의 코드를 추가할 수 있으며, 이것은 이와 같은 작업을 수행하는 데 매우 유용합니다.

즉, Nextflow의 내장 `size()` 함수를 사용하여 `input_files` 배열의 파일 수를 가져오고 `echo` 명령으로 결과를 파일에 쓸 수 있습니다.

`collectGreetings` 프로세스 블록에서 다음 코드 변경을 수행하십시오:

=== "수정 후"

    ```groovy title="hello-workflow.nf" linenums="55" hl_lines="2 5"
        script:
        count_greetings = input_files.size()
        """
        cat ${input_files} > 'COLLECTED-${batch_name}-output.txt'
        echo 'There were ${count_greetings} greetings in this batch.' > '${batch_name}-report.txt'
        """
    ```

=== "수정 전"

    ```groovy title="hello-workflow.nf" linenums="55"
        script:
        """
        cat ${input_files} > 'COLLECTED-${batch_name}-output.txt'
        """
    ```

`count_greetings` 변수는 런타임에 계산됩니다.

#### 4.1.2. 보고 파일 내보내기 및 출력 이름 지정

원칙적으로 해야 할 일은 보고 파일을 `output:` 블록에 추가하는 것뿐입니다.

그러나 그러는 동안 출력 선언에 `emit:` 태그도 추가할 것입니다. 이렇게 하면 위치 인덱스를 사용하지 않고 이름으로 출력을 선택할 수 있습니다.

프로세스 블록에서 다음 코드 변경을 수행하십시오:

=== "수정 후"

    ```groovy title="hello-workflow.nf" linenums="46" hl_lines="2 3"
        output:
        path "COLLECTED-${batch_name}-output.txt", emit: outfile
        path "${batch_name}-report.txt", emit: report
    ```

=== "수정 전"

    ```groovy title="hello-workflow.nf" linenums="46"
        output:
        path "COLLECTED-${batch_name}-output.txt"
    ```

`emit:` 태그는 선택 사항이며 출력 중 하나에만 태그를 추가할 수도 있었습니다.
하지만 속담처럼, 둘 다 하면 안 될 이유가 있습니까?

!!! tip "팁"

    `emit:`를 사용하여 프로세스의 출력에 이름을 지정하지 않으면 각각의 (0 기반) 인덱스를 사용하여 개별적으로 액세스할 수 있습니다.
    예를 들어, `<process>.out[0]`을 사용하여 첫 번째 출력을 가져오고, `<process>.out[1]`을 사용하여 두 번째 출력을 가져오는 식입니다.

    출력에 이름을 지정하는 것을 선호하는 이유는 그렇지 않으면 특히 프로세스가 많은 출력을 생성할 때 실수로 잘못된 인덱스를 가져오기 너무 쉽기 때문입니다.

### 4.2. 워크플로우 출력 업데이트

이제 `collectGreetings` 프로세스에서 두 개의 출력이 나오므로 `collectGreetings.out` 출력에는 두 개의 채널이 포함됩니다:

- `collectGreetings.out.outfile`에는 최종 출력 파일이 포함됩니다
- `collectGreetings.out.report`에는 보고 파일이 포함됩니다

<figure class="excalidraw">
--8<-- "docs/en/docs/hello_nextflow/img/hello-collect-report.svg"
</figure>

이에 따라 워크플로우 출력을 업데이트해야 합니다.

#### 4.2.1. `publish:` 섹션 업데이트

`workflow` 블록에서 다음 코드 변경을 수행하십시오:

=== "수정 후"

    ```groovy title="hello-workflow.nf" linenums="80" hl_lines="4 5"
        publish:
        first_output = sayHello.out
        uppercased = convertToUpper.out
        collected = collectGreetings.out.outfile
        batch_report = collectGreetings.out.report
    ```

=== "수정 전"

    ```groovy title="hello-workflow.nf" linenums="80" hl_lines="4"
        publish:
        first_output = sayHello.out
        uppercased = convertToUpper.out
        collected = collectGreetings.out
    ```

보시다시피, 특정 프로세스 출력을 참조하는 것이 이제 간단해졌습니다.
파트 5(컨테이너)에서 파이프라인에 단계를 하나 더 추가할 때 `collectGreetings.out.outfile`을 쉽게 참조하고 새 프로세스에 전달할 수 있습니다(스포일러: 새 프로세스는 `cowpy`라고 합니다).

하지만 지금은 워크플로우 수준 출력 업데이트를 마무리합시다.

#### 4.2.2. `output` 블록 업데이트

`output` 블록에서 다음 코드 변경을 수행하십시오:

=== "수정 후"

    ```groovy title="hello-workflow.nf" linenums="86" hl_lines="14-17"
    output {
        first_output {
            path 'hello_workflow'
            mode 'copy'
        }
        uppercased {
            path 'hello_workflow'
            mode 'copy'
        }
        collected {
            path 'hello_workflow'
            mode 'copy'
        }
        batch_report {
            path 'hello_workflow'
            mode 'copy'
        }
    }
    ```

=== "수정 전"

    ```groovy title="hello-workflow.nf" linenums="80"
    output {
        first_output {
            path 'hello_workflow'
            mode 'copy'
        }
        uppercased {
            path 'hello_workflow'
            mode 'copy'
        }
        collected {
            path 'hello_workflow'
            mode 'copy'
        }
    }
    ```

`collected` 출력 정의는 이름이 변경되지 않았으므로 업데이트할 필요가 없습니다.
새 출력만 추가하면 됩니다.

### 4.3. 워크플로우 실행

현재 인사말 배치로 실행해 봅시다.

```bash
nextflow run hello-workflow.nf -resume --batch trio
```

??? success "명령 출력"

    ```console
     N E X T F L O W   ~  version 25.10.2

    Launching `hello-workflow.nf` [ecstatic_wilson] DSL2 - revision: c80285f8c8

    executor >  local (1)
    [c5/4c6ca9] sayHello (3)       [100%] 3 of 3, cached: 3 ✔
    [0e/6cbc59] convertToUpper (3) [100%] 3 of 3, cached: 3 ✔
    [02/61ead2] collectGreetings   [100%] 1 of 1 ✔
    ```

`results/hello_workflow/` 디렉토리를 보면 새 보고 파일인 `trio-report.txt`를 찾을 수 있습니다.
열어서 워크플로우가 처리된 인사말 수를 올바르게 보고했는지 확인하십시오.

??? abstract "파일 내용"

    ```txt title="trio-report.txt"
    There were 3 greetings in this batch.
    ```

<figure class="excalidraw">
--8<-- "docs/en/docs/hello_nextflow/img/hello-collect-4-way.svg"
</figure>

CSV에 더 많은 인사말을 추가하고 어떤 일이 발생하는지 테스트해 보십시오.

### 핵심 정리

프로세스가 여러 명명된 출력을 내보내도록 하는 방법과 워크플로우 수준에서 이를 적절하게 처리하는 방법을 알게 되었습니다.

더 일반적으로, 일반적인 방식으로 프로세스를 함께 연결하는 데 관련된 핵심 원칙을 이해하게 되었습니다.

### 다음 단계

충분히 긴 휴식을 취하십시오. 충분히 그럴 자격이 있습니다.

준비가 되면 [**파트 4: Hello Modules**](./04_hello_modules.md)로 이동하여 더 나은 유지보수성과 코드 효율성을 위해 코드를 모듈화하는 방법을 학습하십시오.

---

## 퀴즈

<quiz>
워크플로우 블록에서 프로세스의 출력에 어떻게 액세스합니까?
- [ ] `process.output`
- [ ] `output.processName`
- [x] `processName.out`
- [ ] `get(processName)`

자세히 알아보기: [1.4. 첫 번째 프로세스의 출력을 두 번째 프로세스로 전달](#14-첫-번째-프로세스의-출력을-두-번째-프로세스로-전달)
</quiz>

<quiz>
Nextflow에서 프로세스 실행 순서를 결정하는 것은 무엇입니까?
- [ ] 워크플로우 블록에서 프로세스가 작성된 순서
- [ ] 프로세스 이름의 알파벳 순서
- [x] 프로세스 간의 데이터 의존성
- [ ] 병렬 실행을 위한 무작위 순서

자세히 알아보기: [1.4. 첫 번째 프로세스의 출력을 두 번째 프로세스로 전달](#14-첫-번째-프로세스의-출력을-두-번째-프로세스로-전달)
</quiz>

<quiz>
다운스트림 프로세스를 위해 모든 출력을 단일 목록으로 수집하려면 `???`를 어떤 연산자로 교체해야 합니까?

```groovy hl_lines="4"
workflow {
    greetings_ch = Channel.of('Hello', 'Bonjour', 'Hola')
    SAYHELLO(greetings_ch)
    GATHER_ALL(SAYHELLO.out.???)
}
```

- [ ] `flatten()`
- [x] `collect()`
- [ ] `mix()`
- [ ] `join()`

자세히 알아보기: [2.4. 연산자를 사용하여 인사말을 단일 입력으로 수집](#24-연산자를-사용하여-인사말을-단일-입력으로-수집)
</quiz>

<quiz>
`collect()` 연산자는 언제 사용해야 합니까?
- [ ] 항목을 병렬로 처리하려고 할 때
- [ ] 채널 내용을 필터링해야 할 때
- [x] 다운스트림 프로세스가 업스트림 프로세스의 모든 항목을 필요로 할 때
- [ ] 여러 프로세스에 걸쳐 데이터를 분할하려고 할 때

자세히 알아보기: [2.4. 연산자를 사용하여 인사말을 단일 입력으로 수집](#24-연산자를-사용하여-인사말을-단일-입력으로-수집)
</quiz>

<quiz>
프로세스에서 명명된 출력에 어떻게 액세스합니까?
- [ ] `processName.outputName`
- [ ] `processName.get(outputName)`
- [x] `processName.out.outputName`
- [ ] `output.processName.outputName`

자세히 알아보기: [4.1.2. 보고 파일 내보내기 및 출력 이름 지정](#412-보고-파일-내보내기-및-출력-이름-지정)
</quiz>

<quiz>
프로세스에서 출력의 이름을 지정하는 올바른 구문은 무엇입니까?
- [ ] `name: outputName`
- [ ] `output: outputName`
- [x] `emit: outputName`
- [ ] `label: outputName`

자세히 알아보기: [4.1.2. 보고 파일 내보내기 및 출력 이름 지정](#412-보고-파일-내보내기-및-출력-이름-지정)
</quiz>

<quiz>
프로세스에 여러 입력을 제공할 때 무엇이 참이어야 합니까?
- [ ] 모든 입력이 동일한 유형이어야 함
- [ ] 입력은 알파벳 순서로 제공되어야 함
- [x] 입력 순서가 입력 블록에 정의된 순서와 일치해야 함
- [ ] 한 번에 두 개의 입력만 제공할 수 있음

자세히 알아보기: [3. 프로세스에 추가 매개변수 전달](#3-프로세스에-추가-매개변수-전달)
</quiz>
