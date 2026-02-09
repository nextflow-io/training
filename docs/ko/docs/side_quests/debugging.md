# 디버깅

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } AI 지원 번역 - [자세히 알아보기 및 개선 제안하기](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

Nextflow를 처음 시작할 때, 가장 큰 과제는 스크립트가 올바르게 실행되도록 하는 것입니다. Nextflow 문법에 익숙하지 않을 때는 스크립트 작성 중에 오류를 만들기 쉽습니다.

그러나 문법에 익숙해졌더라도, 문제는 완전히 사라지지 않습니다. 워크플로우를 구축할 때 논리적 오류는 여전히 문제가 될 수 있습니다.

이 사이드 퀘스트에서는 이러한 문제에 체계적으로 접근하는 방법을 배웁니다.
우리는 두 가지 주요 도구로 무장할 것입니다:

1. 일반적인 Nextflow 오류 유형 인식
2. 그것들을 식별하고 해결하는 체계적인 방법

이 두 도구는 디버깅 세션의 길이를 몇 시간에서 몇 분으로 줄이는 데 도움이 될 것입니다.

## 1. 문법 오류

가장 쉽게 식별할 수 있는 오류 유형은 문법 오류입니다. 이것은 Nextflow 스크립트에 문법적으로 올바르지 않은 무언가가 있을 때 발생합니다. 이 오류는 워크플로우 시작 시 즉시 나타나며, 우선 해결해야 할 문제입니다. 컴파일러는 이러한 오류를 찾는 데 도움이 되며, 워크플로우 실행 전에 표면화시킵니다.

### 1.1. 괄호 누락

가장 흔한 문법 오류 중 하나는 괄호나 중괄호를 빠뜨리는 것입니다. Nextflow는 DSL2에서 `{}` 중괄호를 다양한 곳에서 사용합니다: 프로세스 정의, 워크플로우 정의 및 로컬 스코프를 설정하는 모든 위치.

#### 파이프라인 실행

```bash
nextflow run bad_syntax.nf
```

??? failure "명령어 출력"

    ```console
    Error bad_syntax.nf:24:1: Unexpected input: '<EOF>'

    ERROR ~ Script compilation failed
    ```

#### 코드 확인

`<EOF>` 오류는 종종 누락된 괄호를 나타냅니다. 파일을 확인해 보겠습니다:

```groovy title="bad_syntax.nf" linenums="1"
#!/usr/bin/env nextflow

process PROCESS_FILES {
    input:
    tuple val(sample_name), path(file_name)

    output:
    path "${sample_name}_output.txt"

    script:
    """
    echo "Processing ${sample_name}"
    cat ${file_name} > ${sample_name}_output.txt
    """
// } <- 이 괄호가 누락됨!

workflow {
    input_ch = channel
        .fromPath("data/sample_data.csv")
        .splitCsv(header: true)
        .map { row -> [row.sample_id, file(row.file_path)] }

    PROCESS_FILES(input_ch)
}
```

누락된 것이 보이시나요? `PROCESS_FILES` 프로세스 정의의 닫는 중괄호(`}`)가 누락되었습니다. 이러한 상황에서는 일반적으로 괄호가 열린 위치로 돌아가서 닫는 괄호가 있는지 확인해야 합니다.

Nextflow가 괄호 일치를 확인하는 방법으로 인해, 오류 메시지는 종종 `Unexpected input: '<EOF>'`(예상치 못한 입력: 파일 끝)와 같은 일반적인 형태를 취합니다. 이는 파서가 파일 끝에 도달했지만 아직 닫히지 않은 구조가 있다는 것을 의미합니다.

#### 코드 수정

누락된 괄호를 추가합니다:

=== "수정 후"

    ```groovy title="bad_syntax.nf" hl_lines="14" linenums="1"
    #!/usr/bin/env nextflow

    process PROCESS_FILES {
        input:
        tuple val(sample_name), path(file_name)

        output:
        path "${sample_name}_output.txt"

        script:
        """
        echo "Processing ${sample_name}"
        cat ${file_name} > ${sample_name}_output.txt
        """
    } // 괄호 추가됨

    workflow {
        input_ch = channel
            .fromPath("data/sample_data.csv")
            .splitCsv(header: true)
            .map { row -> [row.sample_id, file(row.file_path)] }

        PROCESS_FILES(input_ch)
    }
    ```

=== "수정 전"

    ```groovy title="bad_syntax.nf" linenums="1"
    #!/usr/bin/env nextflow

    process PROCESS_FILES {
        input:
        tuple val(sample_name), path(file_name)

        output:
        path "${sample_name}_output.txt"

        script:
        """
        echo "Processing ${sample_name}"
        cat ${file_name} > ${sample_name}_output.txt
        """
    // } <- 이 괄호가 누락됨!

    workflow {
        input_ch = channel
            .fromPath("data/sample_data.csv")
            .splitCsv(header: true)
            .map { row -> [row.sample_id, file(row.file_path)] }

        PROCESS_FILES(input_ch)
    }
    ```

#### 파이프라인 실행

```bash
nextflow run bad_syntax.nf
```

??? success "명령어 출력"

    ```console
    N E X T F L O W  ~  version 25.10.2

    Launching `bad_syntax.nf` [angry_jennings] DSL2 - revision: 7e72f37c0d

    executor >  local (3)
    [8d/72cee4] process > PROCESS_FILES (1) [100%] 3 of 3 ✔
    ```

### 1.2. 쉘 변수를 제대로 이스케이프하지 않음

Nextflow에서 자주 발생하는 또 다른 오류는 스크립트 내에서 Bash 변수와 Nextflow 변수의 충돌입니다. Nextflow 프로세스는 스크립트 내의 `${variable}`을 Nextflow 변수로 해석합니다. Bash 변수를 사용하려면 달러 기호를 이스케이프해야 합니다.

#### 파이프라인 실행

```bash
nextflow run bad_var.nf
```

??? failure "명령어 출력"

    ```console
    N E X T F L O W  ~  version 25.10.2

    Launching `bad_var.nf` [boring_ampere] DSL2 - revision: 5c28df7302

    executor >  local (1)
    [79/50fcd8] process > COUNT_LINES (1)   [100%] 1 of 1, failed: 1

    Error executing process > 'COUNT_LINES (1)'

    Caused by:
      No such variable: variable1

    Command executed:

      for i in {1..5}; do
        echo "This is line $i"
      done > output.txt

      wc -l output.txt

    Command exit status:
      No such variable: variable1

    Command output:
      (empty)

    Work dir:
      /workspaces/training/side-quests/debugging/work/79/50fcd80a84c0c9b7eff21ad9d0ea42

    Tip: view the complete command output by changing to the process work dir and entering the command `cat .command.out`
    Tip: when you have fixed the problem you can continue the execution adding the option `-resume` to the run command line
     -- Check '.nextflow.log' file for details

    executor >  local (1)
    [79/50fcd8] process > COUNT_LINES (1)   [100%] 1 of 1, failed: 1
    ```

#### 코드 확인

오류는 `No such variable: variable1`인데, 이는 종종 Nextflow와 Bash 변수 문법이 충돌할 때 발생합니다. 코드를 확인해 보겠습니다:

```groovy title="bad_var.nf" linenums="8"
script:
"""
for i in {1..5}; do
  echo "This is line $i"
done > output.txt

wc -l output.txt
"""
```

문제는 `$i`가 Bash 변수로 사용되었지만, Nextflow는 이를 자체 변수로 해석하려고 한다는 것입니다. Nextflow는 `i`라는 변수를 찾으려고 하지만 존재하지 않아서 오류가 발생합니다.

#### 코드 수정

Bash 변수를 Nextflow 변수와 구분하기 위해 달러 기호를 이스케이프합니다:

=== "수정 후"

    ```groovy title="bad_var.nf" hl_lines="4" linenums="8"
    script:
    """
    for i in {1..5}; do
      echo "This is line \$i"
    done > output.txt

    wc -l output.txt
    """
    ```

=== "수정 전"

    ```groovy title="bad_var.nf" hl_lines="4" linenums="8"
    script:
    """
    for i in {1..5}; do
      echo "This is line $i"
    done > output.txt

    wc -l output.txt
    """
    ```

#### 파이프라인 실행

```bash
nextflow run bad_var.nf
```

??? success "명령어 출력"

    ```console
    N E X T F L O W  ~  version 25.10.2
    Launching `bad_var.nf` [sleepy_ampere] DSL2 - revision: 5c28df7302

    executor >  local (1)
    [dd/6a95fd] process > COUNT_LINES [100%] 1 of 1 ✔
    5
    ```

축하합니다! 이제 파이프라인이 성공적으로 실행되고 기대한 결과를 반환합니다.

!!! tip

    다음과 같은 상황에서 Bash 변수를 이스케이프해야 합니다:

    - 루프 변수 (예: `for i in {...}`)
    - 명령어 치환 (예: `\$(date)`)
    - 환경 변수 (예: `\${PATH}`)

    하지만 Nextflow 변수를 참조하려면 이스케이프하지 마세요:

    - 입력 값 (예: `${sample_id}`)
    - 파라미터 (예: `${params.input}`)
    - 정의한 변수 (예: `${my_variable}`)

### 1.3. 구문 강조 도구 사용하기

Nextflow 스크립트를 더 쉽게 작성하기 위해 구문 강조 기능이 있는 코드 에디터를 사용할 수 있습니다. VS Code를 사용하는 경우, [Nextflow 확장 프로그램](https://marketplace.visualstudio.com/items?itemName=nextflow.nextflow)을 설치할 수 있습니다.

구문 강조는 다음과 같은 일반적인 문제를 확인하는 데 도움이 됩니다:

- 괄호 불일치
- 오타 있는 키워드
- 변수 이름 오류
- 문법 오류

에디터의 구문 강조 기능을 사용하면 많은 문법 오류를 워크플로우를 실행하기 전에 잡아낼 수 있습니다.

### 중요 포인트

- 문법 오류는 워크플로우 시작 시 나타납니다
- 가장 흔한 문법 오류는 괄호 불일치입니다
- Bash 변수는 이스케이프하세요 (`\$variable`)
- IDE의 구문 강조 기능을 사용하여 문제를 조기에 식별하세요

문법 오류는 일반적으로 해결하기 가장 쉬운 오류이며, 코드를 자세히 살펴보면 대부분 찾을 수 있습니다. 하지만 워크플로우가 성공적으로 컴파일되더라도, 실행 중에 여전히 오류가 발생할 수 있습니다. 이제 이러한 런타임 오류를 살펴보겠습니다.

### 다음 단계

채널 구조 오류에 대해 알아봅시다.

---

## 2. 채널 구조 오류

문법 오류와 달리, 채널 구조 오류는 스크립트가 문법적으로는 올바르지만 채널 데이터의 형태나 사용 방식이 특정 컨텍스트에서 유효하지 않을 때 발생합니다. 이러한 오류는 프로세스를 실행하려고 할 때 나타납니다.

### 2.1. 튜플 구조 불일치

일반적인 오류는 프로세스의 기대하는 입력 형태와 제공된 채널의 형태가 일치하지 않는 경우입니다.

#### 파이프라인 실행

```bash
nextflow run bad_channel_shape.nf
```

??? failure "명령어 출력"

    ```console
    N E X T F L O W  ~  version 25.10.2
    Launching `bad_channel_shape.nf` [condescending_germain] DSL2 - revision: 9be4e32f67

    Process `PROCESS_FILES` declares 1 input channel but 2 were specified

    Error executing process > 'PROCESS_FILES'

    Caused by:
      Process `PROCESS_FILES` declares 1 input channel but 2 were specified

     -- Check script 'bad_channel_shape.nf' at line: 26 or see '.nextflow.log' file for more details
    ```

#### 코드 확인

오류 메시지는 입력 채널 수가 일치하지 않음을 나타냅니다. 스크립트를 살펴보겠습니다:

```groovy title="bad_channel_shape.nf" linenums="1"
#!/usr/bin/env nextflow

process PROCESS_FILES {
    input:
    val sample_name
    // ERROR: 이 프로세스는 하나의 입력만 선언하지만, 호출 시 두 개의 입력이 제공됨

    output:
    path "${sample_name}_output.txt"

    script:
    """
    echo "Processing ${sample_name}"
    echo "Done" > ${sample_name}_output.txt
    """
}

workflow {
    input_ch = channel
        .fromPath("data/sample_data.csv")
        .splitCsv(header: true)
        // ERROR: 이것은 [sample_id, file_path] 튜플을 생성하지만,
        // 프로세스는 단일 값만 예상함
        .map { row -> [row.sample_id, row.file_path] }

    PROCESS_FILES(input_ch)
}
```

이 코드에는 입력 불일치가 있습니다:

1. `PROCESS_FILES`는 단일 값(`val sample_name`)을 입력으로 기대합니다
2. `workflow`는 `[row.sample_id, row.file_path]` 튜플을 `input_ch`에 제공합니다
3. 프로세스를 호출할 때, 우리는 두 요소를 가진 튜플을 전달하고 있지만 프로세스는 하나의 값만 기대합니다

#### 코드 수정

문제를 해결하는 두 가지 접근 방식이 있습니다:

##### 옵션 1: 프로세스의 입력을 변경합니다.

=== "수정 후"

    ```groovy title="bad_channel_shape.nf" hl_lines="4 5" linenums="1"
    process PROCESS_FILES {
        input:
        // 두 값을 받도록 프로세스 입력을 수정
        val sample_name
        val file_path

        output:
        path "${sample_name}_output.txt"

        script:
        """
        echo "Processing ${sample_name}"
        echo "Done" > ${sample_name}_output.txt
        """
    }
    ```

=== "수정 전"

    ```groovy title="bad_channel_shape.nf" linenums="1"
    process PROCESS_FILES {
        input:
        val sample_name

        output:
        path "${sample_name}_output.txt"

        script:
        """
        echo "Processing ${sample_name}"
        echo "Done" > ${sample_name}_output.txt
        """
    }
    ```

##### 옵션 2: 채널의 출력을 변경합니다.

=== "수정 후"

    ```groovy title="bad_channel_shape.nf" hl_lines="5" linenums="21"
        .fromPath("data/sample_data.csv")
        .splitCsv(header: true)
        // 프로세스가 기대하는 대로 첫 번째 필드만 매핑
        .map { row -> row.sample_id }
    ```

=== "수정 전"

    ```groovy title="bad_channel_shape.nf" linenums="21"
        .fromPath("data/sample_data.csv")
        .splitCsv(header: true)
        // ERROR: 이것은 [sample_id, file_path] 튜플을 생성하지만,
        // 프로세스는 단일 값만 예상함
        .map { row -> [row.sample_id, row.file_path] }
    ```

두 번째 옵션을 선택해 보겠습니다.

#### 파이프라인 실행

```bash
nextflow run bad_channel_shape.nf
```

??? success "명령어 출력"

    ```console
    N E X T F L O W  ~  version 25.10.2
    Launching `bad_channel_shape.nf` [jolly_keller] DSL2 - revision: 2989b88b66

    executor >  local (3)
    [f8/bbc63d] process > PROCESS_FILES (sample3) [100%] 3 of 3 ✔
    ```

### 2.2. 튜플 구조 불일치 (복잡한 경우)

때로는 튜플 불일치 오류가 더 미묘할 수 있습니다. 특히 `.map`이나 다른 채널 연산자를 사용할 때 더욱 그렇습니다.

#### 파이프라인 실행

```bash
nextflow run bad_tuple.nf
```

??? failure "명령어 출력"

    ```console
    N E X T F L O W  ~  version 25.10.2
    Launching `bad_tuple.nf` [furious_feynman] DSL2 - revision: 6c4b5c77a9

    Process `PROCESS_FILES` expects [val sample_name, path file_name] as input

    Error executing process > 'PROCESS_FILES (1)'

    Caused by:
      Process `PROCESS_FILES` expects [val sample_name, path file_name] as input

     -- Check script 'bad_tuple.nf' at line: 27 or see '.nextflow.log' file for more details
    ```

#### 코드 확인

이 오류 메시지는 프로세스가 `[val sample_name, path file_name]` 튜플을 기대하지만 다른 것이 제공되었음을 나타냅니다:

```groovy title="bad_tuple.nf" linenums="1" hl_lines="28"
#!/usr/bin/env nextflow

process PROCESS_FILES {
    input:
    tuple val(sample_name), path(file_name)

    output:
    path "${sample_name}_output.txt"

    script:
    """
    echo "Processing ${sample_name}"
    cat ${file_name} > ${sample_name}_output.txt
    """
}

workflow {
    // 이 예에서 우리는 세 가지 샘플을 생성합니다
    input_ch = channel.of(
        [
            'sample1',
            file("${projectDir}/data/samples/sample1.txt")
        ],
        [
            'sample2',
            file("${projectDir}/data/samples/sample2.txt")
        ]
    )

    // 오류: 튜플이 평탄화되지 않음
    PROCESS_FILES(input_ch)
}
```

문제가 명확하지 않은데, 이것이 이런 종류의 오류가 어려운 이유입니다.

`input_ch`의 내용을 보면 어떻게 보일까요? 다음 수정을 확인해 봅시다:

=== "수정 후"

    ```groovy title="bad_tuple_viewed.nf" hl_lines="27-29" linenums="1"
    #!/usr/bin/env nextflow

    process PROCESS_FILES {
        input:
        tuple val(sample_name), path(file_name)

        output:
        path "${sample_name}_output.txt"

        script:
        """
        echo "Processing ${sample_name}"
        cat ${file_name} > ${sample_name}_output.txt
        """
    }

    workflow {
        // 이 예에서 우리는 세 가지 샘플을 생성합니다
        input_ch = channel.of(
            [
                'sample1',
                file("${projectDir}/data/samples/sample1.txt")
            ],
            [
                'sample2',
                file("${projectDir}/data/samples/sample2.txt")
            ]
        )
        input_ch.view() // 채널 내용을 확인합니다

        PROCESS_FILES(input_ch)
    }
    ```

=== "수정 전"

    ```groovy title="bad_tuple.nf" linenums="1"
    #!/usr/bin/env nextflow

    process PROCESS_FILES {
        input:
        tuple val(sample_name), path(file_name)

        output:
        path "${sample_name}_output.txt"

        script:
        """
        echo "Processing ${sample_name}"
        cat ${file_name} > ${sample_name}_output.txt
        """
    }

    workflow {
        // 이 예에서 우리는 세 가지 샘플을 생성합니다
        input_ch = channel.of(
            [
                'sample1',
                file("${projectDir}/data/samples/sample1.txt")
            ],
            [
                'sample2',
                file("${projectDir}/data/samples/sample2.txt")
            ]
        )

        // 오류: 튜플이 평탄화되지 않음
        PROCESS_FILES(input_ch)
    }
    ```

#### 파이프라인 실행

```bash
nextflow run bad_tuple_viewed.nf
```

??? failure "명령어 출력"

    ```console
    N E X T F L O W  ~  version 25.10.2
    Launching `bad_tuple_viewed.nf` [fervent_williams] DSL2 - revision: 1f7ffeab2e

    [[sample1, /workspaces/training/side-quests/debugging/data/samples/sample1.txt], [sample2, /workspaces/training/side-quests/debugging/data/samples/sample2.txt]]
    Process `PROCESS_FILES` expects [val sample_name, path file_name] as input

    Error executing process > 'PROCESS_FILES (1)'

    Caused by:
      Process `PROCESS_FILES` expects [val sample_name, path file_name] as input

     -- Check script 'bad_tuple_viewed.nf' at line: 31 or see '.nextflow.log' file for more details
    ```

`.view()` 메소드의 출력을 확인해 보세요. 출력은 중첩된 배열입니다! 이것은 `input_ch`가 각 항목이 튜플인 튜플임을 나타냅니다.

```
[[sample1, /workspaces/.../sample1.txt], [sample2, /workspaces/.../sample2.txt]]
```

그러나 프로세스는 중첩되지 않은 튜플 스트림을 기대합니다. 즉, 각 튜플이 샘플 이름과 파일 이름을 포함하는 스트림을 기대합니다:

```
[sample1, /workspaces/.../sample1.txt]
[sample2, /workspaces/.../sample2.txt]
```

#### 코드 수정

`channel.of`를 사용하여 배열 목록을 만들 때, 값 자체 대신 단일 배열로 래핑하면 Nextflow는 배열 자체를 단일 항목으로 취급할 수 있습니다. 우리는 `flatten()` 연산자를 사용하여 이 중첩 구조를 풀 수 있습니다:

=== "수정 후"

    ```groovy title="bad_tuple_viewed.nf" hl_lines="28" linenums="1"
    #!/usr/bin/env nextflow

    process PROCESS_FILES {
        input:
        tuple val(sample_name), path(file_name)

        output:
        path "${sample_name}_output.txt"

        script:
        """
        echo "Processing ${sample_name}"
        cat ${file_name} > ${sample_name}_output.txt
        """
    }

    workflow {
        input_ch = channel.of(
            [
                'sample1',
                file("${projectDir}/data/samples/sample1.txt")
            ],
            [
                'sample2',
                file("${projectDir}/data/samples/sample2.txt")
            ]
        )
        .flatten().collate(2)  // 리스트를 평탄화하고 2개씩 그룹화
        .view() // 수정된 채널 내용을 확인합니다

        PROCESS_FILES(input_ch)
    }
    ```

=== "수정 전"

    ```groovy title="bad_tuple_viewed.nf" linenums="1"
    #!/usr/bin/env nextflow

    process PROCESS_FILES {
        input:
        tuple val(sample_name), path(file_name)

        output:
        path "${sample_name}_output.txt"

        script:
        """
        echo "Processing ${sample_name}"
        cat ${file_name} > ${sample_name}_output.txt
        """
    }

    workflow {
        input_ch = channel.of(
            [
                'sample1',
                file("${projectDir}/data/samples/sample1.txt")
            ],
            [
                'sample2',
                file("${projectDir}/data/samples/sample2.txt")
            ]
        )
        input_ch.view() // 채널 내용을 확인합니다

        PROCESS_FILES(input_ch)
    }
    ```

또는, 채널을 더 명확하게 선언할 수도 있습니다:

```groovy
input_ch = channel.of(
    ['sample1', file("${projectDir}/data/samples/sample1.txt")],
    ['sample2', file("${projectDir}/data/samples/sample2.txt")]
)
```

#### 파이프라인 실행

```bash
nextflow run bad_tuple_viewed.nf
```

??? success "명령어 출력"

    ```console
    N E X T F L O W  ~  version 25.10.2
    Launching `bad_tuple_viewed.nf` [drunk_stone] DSL2 - revision: 69a6b0347e

    [sample1, /workspaces/training/side-quests/debugging/data/samples/sample1.txt]
    [sample2, /workspaces/training/side-quests/debugging/data/samples/sample2.txt]
    executor >  local (2)
    [6d/15b5df] process > PROCESS_FILES (sample2) [100%] 2 of 2 ✔
    ```

이제 `.view()` 출력은 예상되는 형태를 보여주고 워크플로우가 성공적으로 실행됩니다.

### 2.3. 채널 구조 오류 디버깅 기술

`.view()` 메소드는 채널 구조 오류를 디버깅하는 강력한 도구입니다. `.view()`는 채널의 내용을 보여줘서 예상과 다른 경우 확인할 수 있게 해줍니다.

실용적인 접근 방식으로, `.view()` 메소드를 사용하여 전체 데이터 흐름을 단계별로 검사할 수 있습니다. `bad_channel_shape_viewed.nf` 파일을 확인해 보겠습니다:

```groovy title="bad_channel_shape_viewed.nf" linenums="1"
#!/usr/bin/env nextflow

process PROCESS_FILES {
    input:
    val sample_name

    output:
    path "${sample_name}_output.txt"

    script:
    """
    echo "Processing ${sample_name}"
    echo "Done" > ${sample_name}_output.txt
    """
}

workflow {
    input_ch = channel.of(
        ['sample1', 'file1.txt'],
        ['sample2', 'file2.txt'],
        ['sample3', 'file3.txt']
        ) // [sample_name, file_name]
        .view { "Channel content: $it" }
        .map { tuple -> tuple[0] } // sample_name
        .view { "After mapping: $it" }

    PROCESS_FILES(input_ch)
}
```

#### 파이프라인 실행

```bash
nextflow run bad_channel_shape_viewed.nf
```

??? success "명령어 출력"

    ```console
    N E X T F L O W  ~  version 25.10.2

    Launching `bad_channel_shape_viewed.nf` [marvelous_koch] DSL2 - revision: 03e79cdbad

    executor >  local (3)
    [ff/d67cec] PROCESS_FILES (3) | 3 of 3 ✔
    Channel content: [sample1, file1.txt]
    Channel content: [sample2, file2.txt]
    Channel content: [sample3, file3.txt]
    After mapping: sample1
    After mapping: sample2
    After mapping: sample3
    ```

이것은 워크플로우가 복잡해지고 채널 구조가 더 불분명해짐에 따라 더 중요해질 것입니다.

#### 파이프라인 실행

```bash
nextflow run bad_channel_shape_viewed.nf
```

??? success "명령어 출력"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `bad_channel_shape_viewed.nf` [marvelous_koch] DSL2 - revision: 03e79cdbad

    executor >  local (3)
    [ff/d67cec] PROCESS_FILES (2) | 3 of 3 ✔
    Channel content: [sample1, file1.txt]
    Channel content: [sample2, file2.txt]
    Channel content: [sample3, file3.txt]
    After mapping: sample1
    After mapping: sample2
    After mapping: sample3
    ```

### 중요 포인트

많은 채널 구조 오류는 유효한 Nextflow 문법으로 생성될 수 있습니다. 데이터 흐름을 이해하고, 검사를 위한 `.view()` 연산자를 사용하여 채널 구조 오류를 디버깅할 수 있습니다. 그리고 예기치 않은 튜플 구조를 나타내는 대괄호와 같은 오류 메시지 패턴을 인식할 수 있습니다.

### 다음 단계

프로세스 정의에 의해 생성된 오류에 대해 알아봅시다.

---

## 3. 프로세스 구조 오류

프로세스와 관련하여 발생하는 대부분의 오류는 명령을 형성할 때 만든 실수나 기본 소프트웨어와 관련된 문제일 것입니다. 그렇다고 해도, 위의 채널 문제와 유사하게, 문법 오류는 아니지만 런타임에 오류를 일으키는 프로세스 정의에서 실수를 할 수 있습니다.

### 3.1. 출력 파일 누락

프로세스를 작성할 때 흔한 오류 중 하나는 프로세스가 기대하는 것과 생성되는 것 사이에 불일치를 초래하는 무언가를 하는 것입니다.

#### 파이프라인 실행

```bash
nextflow run missing_output.nf
```

??? failure "명령어 출력"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `missing_output.nf` [zen_stone] DSL2 - revision: 37ff61f926

    executor >  local (3)
    executor >  local (3)
    [fd/2642e9] process > PROCESS_FILES (2) [ 66%] 2 of 3, failed: 2
    ERROR ~ Error executing process > 'PROCESS_FILES (3)'

    Caused by:
      Missing output file(s) `sample3.txt` expected by process `PROCESS_FILES (3)`


    Command executed:

      echo "Processing sample3" > sample3_output.txt

    Command exit status:
      0

    Command output:
      (empty)

    Work dir:
      /workspaces/training/side-quests/debugging/work/02/9604d49fb8200a74d737c72a6c98ed

    Tip: when you have fixed the problem you can continue the execution adding the option `-resume` to the run command line

    -- Check '.nextflow.log' file for details
    ```

#### 코드 확인

오류 메시지는 프로세스가 `sample3.txt`라는 출력 파일을 기대했지만, 스크립트는 실제로 `sample3_output.txt`를 생성한다고 표시합니다. `missing_output.nf`의 프로세스 정의를 확인해 보겠습니다:

```groovy title="missing_output.nf" linenums="3" hl_lines="6 10"
process PROCESS_FILES {
    input:
    val sample_name

    output:
    path "${sample_name}.txt"  // 기대: sample3.txt

    script:
    """
    echo "Processing ${sample_name}" > ${sample_name}_output.txt  // 생성: sample3_output.txt
    """
}
```

출력 블록의 파일 이름과 스크립트에서 사용된 파일 이름 사이에 불일치가 있음을 알 수 있습니다. 이 불일치가 프로세스 실패의 원인입니다. 이런 종류의 오류가 발생하면, 프로세스 정의와 출력 블록 사이에 출력이 일치하는지 확인하세요.

문제가 여전히 명확하지 않다면, 작업 디렉토리 자체를 확인하여 실제로 생성된 출력 파일을 식별할 수 있습니다:

```bash
❯ ls -h work/02/9604d49fb8200a74d737c72a6c98ed
sample3_output.txt
```

이 예제에서는 우리의 `output:` 정의와 달리 출력 파일 이름에 `_output` 접미사가 포함되어 있음을 확인할 수 있습니다.

#### 코드 수정

출력 파일 이름을 일치시켜 불일치를 수정합니다:

=== "수정 후"

    ```groovy title="missing_output.nf" hl_lines="6 10" linenums="3"
    process PROCESS_FILES {
        input:
        val sample_name

        output:
        path "${sample_name}_output.txt"  // 수정: 스크립트 출력과 일치

        script:
        """
        echo "Processing ${sample_name}" > ${sample_name}_output.txt
        """
    }
    ```

=== "수정 전"

    ```groovy title="missing_output.nf" hl_lines="6 10" linenums="3"
    process PROCESS_FILES {
        input:
        val sample_name

        output:
        path "${sample_name}.txt"  // 기대: sample3.txt

        script:
        """
        echo "Processing ${sample_name}" > ${sample_name}_output.txt  // 생성: sample3_output.txt
        """
    }
    ```

#### 파이프라인 실행

```bash
nextflow run missing_output.nf
```

??? success "명령어 출력"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `missing_output.nf` [elated_hamilton] DSL2 - revision: 961938ee2b

    executor >  local (3)
    [16/1c437c] PROCESS_FILES (3) | 3 of 3 ✔
    ```

### 3.2. 소프트웨어 누락

다른 유형의 오류는 소프트웨어 프로비저닝 실수로 인해 발생합니다. `missing_software.nf`는 문법적으로 유효한 워크플로우이지만, 사용하는 `cowpy` 명령을 제공하는 외부 소프트웨어에 의존합니다.

#### 파이프라인 실행

```bash
nextflow run missing_software.nf
```

??? failure "명령어 출력"

    ```console hl_lines="12 18"
    ERROR ~ Error executing process > 'PROCESS_FILES (3)'

    Caused by:
      Process `PROCESS_FILES (3)` terminated with an error exit status (127)


    Command executed:

      cowpy sample3 > sample3_output.txt

    Command exit status:
      127

    Command output:
      (empty)

    Command error:
      .command.sh: line 2: cowpy: command not found

    Work dir:
      /workspaces/training/side-quests/debugging/work/82/42a5bfb60c9c6ee63ebdbc2d51aa6e

    Tip: you can try to figure out what's wrong by changing to the process work directory and showing the script file named `.command.sh`

    -- Check '.nextflow.log' file for details
    ```

프로세스가 우리가 지정한 명령에 접근할 수 없습니다. 때로는 이것이 워크플로우 `bin` 디렉토리에 스크립트가 있지만 실행 가능하지 않기 때문일 수 있습니다. 또 다른 경우에는 워크플로우가 실행되는 컨테이너나 환경에 소프트웨어가 설치되지 않았기 때문일 수 있습니다.

#### 코드 확인

그 `127` 종료 코드를 주목하세요 - 이것이 정확히 문제를 알려줍니다. `missing_software.nf`를 살펴보겠습니다:

```groovy title="missing_software.nf" linenums="3" hl_lines="3"
process PROCESS_FILES {

    container 'community.wave.seqera.io/library/cowpy:1.1.5--3db457ae1977a273'

    input:
    val sample_name

    output:
    path "${sample_name}_output.txt"

    script:
    """
    cowpy ${sample_name} > ${sample_name}_output.txt
    """
}
```

#### 코드 수정

여기서는 약간 부정직했습니다. 사실 코드에는 아무런 문제가 없습니다. 문제의 명령에 접근할 수 있도록 프로세스를 실행하는 데 필요한 구성을 지정해야 합니다. 이 경우 프로세스에는 컨테이너 정의가 있으므로, Docker를 활성화하여 워크플로우를 실행하기만 하면 됩니다.

#### 파이프라인 실행

`nextflow.config`에 Docker 프로필을 설정해 두었으므로 다음과 같이 워크플로우를 실행할 수 있습니다:

```bash
nextflow run missing_software.nf -profile docker
```

??? success "명령어 출력"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `missing_software.nf` [awesome_stonebraker] DSL2 - revision: 0296d12839

    executor >  local (3)
    [38/ab20d1] PROCESS_FILES (1) | 3 of 3 ✔
    ```

!!! note

    Nextflow가 컨테이너를 사용하는 방법에 대해 더 알아보려면 [Hello Nextflow](../hello_nextflow/05_hello_containers.md)를 참조하세요.

### 3.3. 잘못된 리소스 구성

프로덕션 사용에서는 프로세스에 리소스를 구성하게 될 것입니다. 예를 들어 `memory`는 프로세스에 사용 가능한 최대 메모리를 정의하며, 프로세스가 이 한도를 초과하면 일반적으로 스케줄러는 프로세스를 종료하고 종료 코드 `137`을 반환합니다. 우리는 `local` 실행자를 사용하고 있기 때문에 그것을 여기서 보여줄 수 없지만, `time`으로 비슷한 것을 보여줄 수 있습니다.

#### 파이프라인 실행

`bad_resources.nf`는 1밀리초라는 비현실적인 시간 제한이 있는 프로세스 구성을 가지고 있습니다:

```bash
nextflow run bad_resources.nf -profile docker
```

??? failure "명령어 출력"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `bad_resources.nf` [disturbed_elion] DSL2 - revision: 27d2066e86

    executor >  local (3)
    [c0/ded8e1] PROCESS_FILES (3) | 0 of 3 ✘
    ERROR ~ Error executing process > 'PROCESS_FILES (2)'

    Caused by:
      Process exceeded running time limit (1ms)

    Command executed:

      cowpy sample2 > sample2_output.txt

    Command exit status:
      -

    Command output:
      (empty)

    Work dir:
      /workspaces/training/side-quests/debugging/work/53/f0a4cc56d6b3dc2a6754ff326f1349

    Container:
      community.wave.seqera.io/library/cowpy:1.1.5--3db457ae1977a273

    Tip: you can replicate the issue by changing to the process work dir and entering the command `bash .command.run`

     -- Check '.nextflow.log' file for details
    ```

#### 코드 확인

`bad_resources.nf`를 살펴보겠습니다:

```groovy title="bad_resources.nf" linenums="3" hl_lines="3"
process PROCESS_FILES {

    time '1 ms'  // 오류: 비현실적인 시간 제한

    input:
    val sample_name

    output:
    path "${sample_name}_output.txt"

    script:
    """
    sleep 1  // 1초가 소요되지만, 시간 제한은 1ms
    cowpy ${sample_name} > ${sample_name}_output.txt
    """
}
```

프로세스가 1초 이상 걸릴 것이라는 것을 알고 있습니다(그것을 확실히 하기 위해 sleep을 추가했습니다), 하지만 프로세스는 1밀리초 후에 시간 초과되도록 설정되어 있습니다. 누군가 구성에서 약간 비현실적이었습니다!

#### 코드 수정

시간 제한을 현실적인 값으로 늘립니다:

=== "수정 후"

    ```groovy title="bad_resources.nf" hl_lines="3" linenums="3"
    process PROCESS_FILES {

        time '100 s'  // 수정: 현실적인 시간 제한

        input:
        val sample_name

        output:
        path "${sample_name}_output.txt"

        script:
        """
        sleep 1
        cowpy ${sample_name} > ${sample_name}_output.txt
        """
    }
    ```

=== "수정 전"

    ```groovy title="bad_resources.nf" hl_lines="3" linenums="3"
    process PROCESS_FILES {

        time '1 ms'  // 오류: 비현실적인 시간 제한

        input:
        val sample_name

        output:
        path "${sample_name}_output.txt"

        script:
        """
        sleep 1  // 1초가 소요되지만, 시간 제한은 1ms
        cowpy ${sample_name} > ${sample_name}_output.txt
        """
    }
    ```

#### 파이프라인 실행

```bash
nextflow run bad_resources.nf -profile docker
```

??? success "명령어 출력"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `bad_resources.nf` [friendly_mcclintock] DSL2 - revision: 381567d2c1

    executor >  local (3)
    [c2/9b4c41] PROCESS_FILES (3) | 3 of 3 ✔
    ```

오류 메시지를 잘 읽으면 이와 같은 실패에 대해 너무 오래 고민할 필요가 없을 것입니다. 하지만 실행 중인 명령의 리소스 요구 사항을 이해하여 리소스 지시어를 적절하게 구성할 수 있도록 하세요.

### 3.4. 프로세스 디버깅 기법

프로세스가 실패하거나 예상치 못한 동작을 보일 때, 무엇이 잘못되었는지 조사하기 위한 체계적인 기법이 필요합니다. 작업 디렉토리에는 프로세스 실행을 디버깅하는 데 필요한 모든 정보가 포함되어 있습니다.

#### 작업 디렉토리 검사 사용

프로세스를 위한 가장 강력한 디버깅 도구는 작업 디렉토리를 검사하는 것입니다. 프로세스가 실패하면, Nextflow는 해당 특정 프로세스 실행을 위한 작업 디렉토리를 생성하며, 이 디렉토리에는 무엇이 잘못되었는지 이해하는 데 필요한 모든 파일이 포함되어 있습니다.

#### 파이프라인 실행

작업 디렉토리 검사를 시연하기 위해 앞서 살펴본 `missing_output.nf` 예제를 사용해 보겠습니다(필요한 경우 출력 이름 불일치를 다시 생성하세요):

```bash
nextflow run missing_output.nf
```

??? failure "명령어 출력"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `missing_output.nf` [irreverent_payne] DSL2 - revision: 3d5117f7e2

    executor >  local (3)
    [5d/d544a4] PROCESS_FILES (2) | 0 of 3 ✘
    ERROR ~ Error executing process > 'PROCESS_FILES (1)'

    Caused by:
      Missing output file(s) `sample1.txt` expected by process `PROCESS_FILES (1)`

    Command executed:

      echo "Processing sample1" > sample1_output.txt

    Command exit status:
      0

    Command output:
      (empty)

    Work dir:
      /workspaces/training/side-quests/debugging/work/1e/2011154d0b0f001cd383d7364b5244

    Tip: you can replicate the issue by changing to the process work dir and entering the command `bash .command.run`
    Tip: when you have fixed the problem you can continue the execution adding the option `-resume` to the run command line
     -- Check '.nextflow.log' file for details

    executor >  local (1)
    [79/50fcd8] process > COUNT_LINES (1)   [100%] 1 of 1, failed: 1
    ```

#### 작업 디렉토리 확인

이 오류가 발생하면 작업 디렉토리에 모든 디버깅 정보가 포함되어 있습니다. 오류 메시지에서 작업 디렉토리 경로를 찾아 그 내용을 검사하세요:

```bash
# 오류 메시지에서 작업 디렉토리 찾기
ls work/02/9604d49fb8200a74d737c72a6c98ed/
```

그런 다음 핵심 파일을 검사할 수 있습니다:

##### 명령 스크립트 확인

`.command.sh` 파일은 정확히 어떤 명령이 실행되었는지 보여줍니다:

```bash
# 실행된 명령 보기
cat work/02/9604d49fb8200a74d737c72a6c98ed/.command.sh
```

이것은 다음을 보여줍니다:

- **변수 대체**: Nextflow 변수가 제대로 확장되었는지
- **파일 경로**: 입력 파일이 올바르게 위치했는지
- **명령 구조**: 스크립트 구문이 올바른지

찾아봐야 할 일반적인 문제:

- **따옴표 누락**: 공백이 포함된 변수는 적절한 따옴표가 필요함
- **잘못된 파일 경로**: 존재하지 않거나 잘못된 위치에 있는 입력 파일
- **잘못된 변수 이름**: 변수 참조의 오타
- **누락된 환경 설정**: 특정 환경에 의존하는 명령

##### 오류 출력 확인

`.command.err` 파일에는 실제 오류 메시지가 포함되어 있습니다:

```bash
# 오류 출력 보기
cat work/02/9604d49fb8200a74d737c72a6c98ed/.command.err
```

이 파일은 다음을 보여줍니다:

- **종료 코드**: 127(명령을 찾을 수 없음), 137(종료됨) 등
- **권한 오류**: 파일 액세스 문제
- **소프트웨어 오류**: 응용 프로그램별 오류 메시지
- **리소스 오류**: 메모리/시간 제한 초과

##### 표준 출력 확인

`.command.out` 파일은 명령이 생성한 내용을 보여줍니다:

```bash
# 표준 출력 보기
cat work/02/9604d49fb8200a74d737c72a6c98ed/.command.out
```

이것은 다음을 확인하는 데 도움이 됩니다:

- **예상 출력**: 명령이 올바른 결과를 생성했는지
- **부분 실행**: 명령이 시작되었지만 중간에 실패했는지
- **디버그 정보**: 스크립트의 진단 출력

##### 종료 코드 확인

`.exitcode` 파일에는 프로세스의 종료 코드가 포함되어 있습니다:

```bash
# 종료 코드 보기
cat work/*/*/.exitcode
```

일반적인 종료 코드와 그 의미:

- **종료 코드 127**: 명령을 찾을 수 없음 - 소프트웨어 설치 확인
- **종료 코드 137**: 프로세스 종료됨 - 메모리/시간 제한 확인

##### 파일 존재 확인

출력 파일이 누락되어 프로세스가 실패할 때, 실제로 어떤 파일이 생성되었는지 확인하세요:

```bash
# 작업 디렉토리의 모든 파일 나열
ls -la work/02/9604d49fb8200a74d737c72a6c98ed/
```

이것은 다음을 식별하는 데 도움이 됩니다:

- **파일 이름 불일치**: 예상과 다른 이름의 출력 파일
- **권한 문제**: 생성할 수 없는 파일
- **경로 문제**: 잘못된 디렉토리에 생성된 파일

앞서 예제에서, 이것은 예상된 `sample3.txt`가 없지만 `sample3_output.txt`가 있다는 것을 확인했습니다:

```bash
❯ ls -h work/02/9604d49fb8200a74d737c72a6c98ed
sample3_output.txt
```

### 중요 포인트

프로세스 디버깅은 무엇이 잘못되었는지 이해하기 위해 작업 디렉토리를 검사해야 합니다. 주요 파일에는 `.command.sh`(실행된 스크립트), `.command.err`(오류 메시지) 및 `.command.out`(표준 출력)이 포함됩니다. 127(명령을 찾을 수 없음) 및 137(프로세스 종료됨)과 같은 종료 코드는 실패 유형에 대한 즉각적인 진단 단서를 제공합니다.

### 다음 단계

Nextflow의 내장 디버깅 도구와 문제 해결을 위한 체계적인 접근 방식에 대해 알아봅시다.

---

## 4. 내장 디버깅 도구 및 고급 기법

Nextflow는 워크플로우 실행을 디버깅하고 분석하기 위한 몇 가지 강력한 내장 도구를 제공합니다. 이러한 도구는 무엇이 잘못되었는지, 어디서 잘못되었는지, 그리고 어떻게 효율적으로 수정할 수 있는지 이해하는 데 도움이 됩니다.

### 4.1. 실시간 프로세스 출력

때로는 실행 중인 프로세스 내부에서 무슨 일이 일어나고 있는지 확인해야 합니다. 실시간 프로세스 출력을 활성화하면 각 작업이 실행될 때 정확히 무엇을 하고 있는지 볼 수 있습니다.

#### 파이프라인 실행

앞서 예제에서 본 `bad_channel_shape_viewed.nf`는 `.view()`를 사용하여 채널 내용을 출력했지만, `debug` 지시어를 사용하여 프로세스 내부에서 변수를 출력할 수도 있으며, 이는 `bad_channel_shape_viewed_debug.nf`에서 보여줍니다. 워크플로우를 실행하세요:

```bash
nextflow run bad_channel_shape_viewed_debug.nf
```

??? success "명령어 출력"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `bad_channel_shape_viewed_debug.nf` [agitated_crick] DSL2 - revision: ea3676d9ec

    executor >  local (3)
    [c6/2dac51] process > PROCESS_FILES (3) [100%] 3 of 3 ✔
    Channel content: [sample1, file1.txt]
    Channel content: [sample2, file2.txt]
    Channel content: [sample3, file3.txt]
    After mapping: sample1
    After mapping: sample2
    After mapping: sample3
    Sample name inside process is sample2

    Sample name inside process is sample1

    Sample name inside process is sample3
    ```

#### 코드 확인

`debug` 지시어가 어떻게 작동하는지 `bad_channel_shape_viewed_debug.nf`를 살펴보겠습니다:

```groovy title="bad_channel_shape_viewed_debug.nf" linenums="3" hl_lines="2"
process PROCESS_FILES {
    debug true  // 실시간 출력 활성화

    input:
    val sample_name

    output:
    path "${sample_name}_output.txt"

    script:
    """
    echo "Sample name inside process is ${sample_name}"
    echo "Processing ${sample_name}" > ${sample_name}_output.txt
    """
}
```

`debug` 지시어는 프로세스의 환경을 이해하는 빠르고 편리한 방법이 될 수 있습니다.

### 4.2. 미리보기 모드

때로는 프로세스가 실행되기 전에 문제를 포착하고 싶을 수 있습니다. Nextflow는 이런 종류의 선제적 디버깅을 위한 플래그를 제공합니다: `-preview`.

#### 파이프라인 실행

미리보기 모드는 명령을 실행하지 않고 워크플로우 로직을 테스트할 수 있게 해줍니다. 이것은 실제 명령을 실행하지 않고도 워크플로우 구조를 빠르게 확인하고 프로세스가 올바르게 연결되어 있는지 확인하는 데 유용할 수 있습니다.

!!! note

    이전에 `bad_syntax.nf`를 수정했다면, 이 명령을 실행하기 전에 스크립트 블록 뒤의 닫는 중괄호를 제거하여 구문 오류를 다시 도입하세요.

다음 명령을 실행하세요:

```bash
nextflow run bad_syntax.nf -preview
```

??? failure "명령어 출력"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `bad_syntax.nf` [magical_mercator] DSL2 - revision: 550b9a8873

    Error bad_syntax.nf:24:1: Unexpected input: '<EOF>'

    ERROR ~ Script compilation failed

     -- Check '.nextflow.log' file for details
    ```

미리보기 모드는 특히 프로세스를 실행하지 않고 초기에 구문 오류를 포착하는 데 유용합니다. 이 모드는 실행 전에 워크플로우 구조와 프로세스 연결을 검증합니다.

### 4.3. 로직 테스트를 위한 스텁 실행

때로 명령이 너무 오래 걸리거나, 특별한 소프트웨어가 필요하거나, 복잡한 이유로 실패하기 때문에 오류를 디버깅하기 어려울 수 있습니다. 스텁 실행을 사용하면 실제 명령을 실행하지 않고 워크플로우 로직을 테스트할 수 있습니다.

#### 파이프라인 실행

Nextflow 프로세스를 개발할 때, `stub` 지시어를 사용하여 실제 명령을 실행하지 않고 올바른 형태의 출력을 생성하는 '더미' 명령을 정의할 수 있습니다. 이 접근 방식은 특히 실제 소프트웨어의 복잡성을 처리하기 전에 워크플로우 로직이 올바른지 확인하고자 할 때 가치가 있습니다.

예를 들어, 앞서 살펴본 `missing_software.nf`를 기억하시나요? `-profile docker`를 추가하기 전까지 워크플로우 실행을 방해했던 누락된 소프트웨어가 있던 것입니다? `missing_software_with_stub.nf`는 매우 유사한 워크플로우입니다. 같은 방식으로 실행하면 같은 오류가 발생합니다:

```bash
nextflow run missing_software_with_stub.nf
```

??? failure "명령어 출력"

    ```console hl_lines="12 18"
    ERROR ~ Error executing process > 'PROCESS_FILES (3)'

    Caused by:
      Process `PROCESS_FILES (3)` terminated with an error exit status (127)


    Command executed:

      cowpy sample3 > sample3_output.txt

    Command exit status:
      127

    Command output:
      (empty)

    Command error:
      .command.sh: line 2: cowpy: command not found

    Work dir:
      /workspaces/training/side-quests/debugging/work/82/42a5bfb60c9c6ee63ebdbc2d51aa6e

    Tip: you can try to figure out what's wrong by changing to the process work directory and showing the script file named `.command.sh`

    -- Check '.nextflow.log' file for details
    ```

그러나 이 워크플로우는 `docker` 프로필 없이도 `-stub-run`으로 실행하면 오류가 발생하지 않습니다:

```bash
nextflow run missing_software_with_stub.nf -stub-run
```

??? success "명령어 출력"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `missing_software_with_stub.nf` [astonishing_shockley] DSL2 - revision: f1f4f05d7d

    executor >  local (3)
    [b5/2517a3] PROCESS_FILES (3) | 3 of 3 ✔
    ```

#### 코드 확인

`missing_software_with_stub.nf`를 살펴보겠습니다:

```groovy title="missing_software.nf (with stub)" hl_lines="16-19" linenums="3"
process PROCESS_FILES {

    container 'community.wave.seqera.io/library/cowpy:1.1.5--3db457ae1977a273'

    input:
    val sample_name

    output:
    path "${sample_name}_output.txt"

    script:
    """
    cowpy ${sample_name} > ${sample_name}_output.txt
    """

    stub:
    """
    touch ${sample_name}_output.txt
    """
}
```

`missing_software.nf`에 비해, 이 프로세스는 Nextflow가 스텁 모드로 실행될 경우 `script:` 대신 사용될 명령을 지정하는 `stub:` 지시어가 있습니다.

여기서 사용하는 `touch` 명령은 어떤 소프트웨어나 적절한 입력에도 의존하지 않으며, 모든 상황에서 실행되어 프로세스 내부에 대해 걱정하지 않고 워크플로우 로직을 디버깅할 수 있습니다.

**스텁 실행이 디버깅에 도움을 주는 것**:

- 채널 구조 및 데이터 흐름
- 프로세스 연결 및 종속성
- 파라미터 전파
- 소프트웨어 종속성 없는 워크플로우 로직

### 4.4. 체계적인 디버깅 접근법

이제 트레이스 파일 및 작업 디렉토리부터 미리보기 모드, 스텁 실행 및 리소스 모니터링에 이르는 개별 디버깅 기법을 배웠으니, 그것들을 체계적인 방법론으로 결합해 보겠습니다. 구조화된 접근 방식을 갖추면 복잡한 오류에 압도당하지 않고 중요한 단서를 놓치지 않게 됩니다.

이 방법론은 우리가 다룬 모든 도구를 효율적인 워크플로우로 결합합니다:

**4단계 디버깅 방법**:

**1단계: 구문 오류 해결 (5분)**

1. VSCode나 IDE의 빨간색 밑줄 확인
2. `nextflow run workflow.nf -preview`를 실행하여 구문 문제 식별
3. 모든 구문 오류 수정 (누락된 괄호, 후행 쉼표 등)
4. 계속 진행하기 전에 워크플로우가 성공적으로 파싱되는지 확인

**2단계: 빠른 평가 (5분)**

1. 런타임 오류 메시지를 주의 깊게 읽기
2. 런타임, 로직 또는 리소스 오류인지 확인
3. 기본 워크플로우 로직을 테스트하기 위해 미리보기 모드 사용

**3단계: 상세 조사 (15-30분)**

1. 실패한 작업의 작업 디렉토리 찾기
2. 로그 파일 검사
3. 채널 검사를 위해 `.view()` 연산자 추가
4. 실행 없이 워크플로우 로직을 테스트하기 위해 `-stub-run` 사용

**4단계: 수정 및 검증 (15분)**

1. 최소한의 대상 수정 만들기
2. 재개로 테스트: `nextflow run workflow.nf -resume`
3. 완전한 워크플로우 실행 확인

!!! tip "효율적인 디버깅을 위한 Resume 사용"

    문제를 식별한 후에는 워크플로우의 성공적인 부분을 다시 실행하는 데 시간을 낭비하지 않고 수정 사항을 테스트할 효율적인 방법이 필요합니다. Nextflow의 `-resume` 기능은 디버깅에 매우 유용합니다.

    [Hello Nextflow](../hello_nextflow/)를 살펴봤다면 `-resume`를 접했을 것이며, 문제 프로세스 전에 있는 프로세스를 기다리지 않도록 디버깅할 때 이것을 잘 활용하는 것이 중요합니다.

    **Resume 디버깅 전략**:

    1. 실패할 때까지 워크플로우 실행
    2. 실패한 작업에 대한 작업 디렉토리 검사
    3. 특정 문제 수정
    4. 수정 사항만 테스트하기 위해 재개
    5. 워크플로우가 완료될 때까지 반복

#### 디버깅 구성 프로필

이 체계적인 접근 방식을 더욱 효율적으로 만들기 위해, 필요한 모든 도구를 자동으로 활성화하는 전용 디버깅 구성을 만들 수 있습니다:

```groovy title="nextflow.config (debug profile)" linenums="1"
profiles {
    debug {
        process {
            debug = true
            cleanup = false

            // 디버깅을 위한 보수적 리소스
            maxForks = 1
            memory = '2.GB'
            cpus = 1
        }
    }
}
```

그런 다음 이 프로필을 활성화하여 파이프라인을 실행할 수 있습니다:

```bash
nextflow run workflow.nf -profile debug
```

이 프로필은 실시간 출력을 활성화하고, 작업 디렉토리를 보존하며, 더 쉬운 디버깅을 위해 병렬화를 제한합니다.

### 4.5. 실용적인 디버깅 연습

이제 체계적인 디버깅 접근 방식을 실습해 볼 시간입니다. 워크플로우 `buggy_workflow.nf`는 실제 개발에서 마주칠 수 있는 유형의 문제를 나타내는 몇 가지 일반적인 오류를 포함하고 있습니다.

!!! exercise

    체계적인 디버깅 접근 방식을 사용하여 `buggy_workflow.nf`의 모든 오류를 식별하고 수정하세요. 이 워크플로우는 CSV 파일에서 샘플 데이터를 처리하려고 하지만, 일반적인 디버깅 시나리오를 나타내는 여러 의도적인 버그를 포함하고 있습니다.

    첫 번째 오류를 보기 위해 워크플로우를 실행하세요:

    ```bash
    nextflow run buggy_workflow.nf
    ```

    ??? failure "명령어 출력"

        ```console
        N E X T F L O W   ~  version 25.10.2

        Launching `buggy_workflow.nf` [wise_ramanujan] DSL2 - revision: d51a8e83fd

        ERROR ~ Range [11, 12) out of bounds for length 11

         -- Check '.nextflow.log' file for details
        ```

        이 암호 같은 오류는 `params{}` 블록 주변 11-12줄에서 파싱 문제를 나타냅니다. v2 파서는 구조적 문제를 일찍 잡아냅니다.

    배운 4단계 디버깅 방법을 적용하세요:

    **1단계: 구문 오류 해결**
    - VSCode나 IDE의 빨간색 밑줄 확인
    - `nextflow run workflow.nf -preview`를 실행하여 구문 문제 식별
    - 모든 구문 오류 수정 (누락된 괄호, 후행 쉼표 등)
    - 계속 진행하기 전에 워크플로우가 성공적으로 파싱되는지 확인

    **2단계: 빠른 평가**
    - 런타임 오류 메시지를 주의 깊게 읽기
    - 오류가 런타임, 로직 또는 리소스 관련인지 식별
    - 기본 워크플로우 로직을 테스트하기 위해 `-preview` 모드 사용

    **3단계: 상세 조사**
    - 실패한 작업에 대한 작업 디렉토리 검사
    - 채널을 검사하기 위해 `.view()` 연산자 추가
    - 작업 디렉토리의 로그 파일 확인
    - 실행 없이 워크플로우 로직을 테스트하기 위해 `-stub-run` 사용

    **4단계: 수정 및 검증**
    - 대상 수정 만들기
    - 수정 사항을 효율적으로 테스트하기 위해 `-resume` 사용
    - 완전한 워크플로우 실행 확인

    **사용할 수 있는 디버깅 도구**:
    ```bash
    # 구문 검사를 위한 미리보기 모드
    nextflow run buggy_workflow.nf -preview

    # 상세 출력을 위한 디버그 프로필
    nextflow run buggy_workflow.nf -profile debug

    # 로직 테스트를 위한 스텁 실행
    nextflow run buggy_workflow.nf -stub-run

    # 수정 후 재개
    nextflow run buggy_workflow.nf -resume
    ```

    ??? solution
        `buggy_workflow.nf`는 (계산 방법에 따라) 9-10개의 뚜렷한 오류를 포함하고 있으며, 모든 주요 디버깅 카테고리를 다룹니다. 각 오류와 수정 방법에 대한 체계적인 분석을 아래에서 확인할 수 있습니다.

        먼저 구문 오류부터 시작하겠습니다:

        **오류 1: 구문 오류 - 후행 쉼표**
        ```groovy linenums="21"
        output:
            path "${sample_id}_result.txt",  // 오류: 후행 쉼표
        ```
        **수정:** 후행 쉼표 제거
        ```groovy linenums="21"
        output:
            path "${sample_id}_result.txt"
        ```

        **오류 2: 구문 오류 - 닫는 괄호 누락**
        ```groovy linenums="24"
        script:
        """
        echo "Processing: ${sample}"
        cat ${input_file} > ${sample}_result.txt
        """
        // 오류: processFiles 프로세스에 대한 닫는 괄호 누락
        ```
        **수정:** 누락된 닫는 괄호 추가
        ```groovy linenums="29"
        """
        echo "Processing: ${sample_id}"
        cat ${input_file} > ${sample_id}_result.txt
        """
        }  // 누락된 닫는 괄호 추가
        ```

        **오류 3: 변수 이름 오류**
        ```groovy linenums="26"
        echo "Processing: ${sample}"     // 오류: sample_id여야 함
        cat ${input_file} > ${sample}_result.txt  // 오류: sample_id여야 함
        ```
        **수정:** 올바른 입력 변수 이름 사용
        ```groovy linenums="26"
        echo "Processing: ${sample_id}"
        cat ${input_file} > ${sample_id}_result.txt
        ```

        **오류 4: 정의되지 않은 변수 오류**
        ```groovy linenums="87"
        heavy_ch = heavyProcess(sample_ids)  // 오류: sample_ids 정의되지 않음
        ```
        **수정:** 올바른 채널 사용 및 샘플 ID 추출
        ```groovy linenums="87"
        heavy_ch = heavyProcess(input_ch)
        ```

        이 시점에서 워크플로우는 실행되지만, `processFiles`에서 `Path value cannot be null`과 같은 오류가 여전히 발생합니다. 이는 잘못된 채널 구조로 인한 것입니다.

        **오류 5: 채널 구조 오류 - 잘못된 맵 출력**
        ```groovy linenums="83"
        .map { row -> row.sample_id }  // 오류: processFiles는 튜플을 기대함
        ```
        **수정:** processFiles가 기대하는 튜플 구조 반환
        ```groovy linenums="83"
        .map { row -> [row.sample_id, file(row.fastq_path)] }
        ```

        하지만 이렇게 하면 위의 `heavyProcess()` 실행이 중단됩니다. 따라서 해당 프로세스에 샘플 ID만 전달하는 맵이 필요합니다:

        **오류 6: heavyProcess에 대한 잘못된 채널 구조**
        ```groovy linenums="87"
        heavy_ch = heavyProcess(input_ch)  // 오류: input_ch는 현재 방출당 2개의 요소를 가짐 - heavyProcess는 첫 번째 요소만 필요함
        ```
        **수정:** 올바른 채널 사용 및 샘플 ID 추출
        ```groovy linenums="87"
        heavy_ch = heavyProcess(input_ch.map{it[0]})
        ```

        이제 더 나아가지만 `No such variable: i`라는 오류가 발생합니다. 이는 Bash 변수를 이스케이프하지 않았기 때문입니다.

        **오류 7: Bash 변수 이스케이프 오류**
        ```groovy linenums="48"
        echo "Heavy computation $i for ${sample_id}"  // 오류: $i 이스케이프되지 않음
        ```
        **수정:** bash 변수 이스케이프
        ```groovy linenums="48"
        echo "Heavy computation \${i} for ${sample_id}"
        ```

        이제 `Process exceeded running time limit (1ms)`라는 오류가 발생하므로, 관련 프로세스의 실행 시간 제한을 수정해야 합니다:

        **오류 8: 리소스 구성 오류**
        ```groovy linenums="36"
        time '1 ms'  // 오류: 비현실적인 시간 제한
        ```
        **수정:** 현실적인 시간 제한으로 증가
        ```groovy linenums="36"
        time '100 s'
        ```

        다음으로 `Missing output file(s)` 오류를 해결해야 합니다:

        **오류 9: 출력 파일 이름 불일치**
        ```groovy linenums="49"
        done > ${sample_id}.txt  // 오류: 잘못된 파일 이름, 출력 선언과 일치해야 함
        ```
        **수정:** 출력 선언과 일치
        ```groovy linenums="49"
        done > ${sample_id}_heavy.txt
        ```

        처음 두 프로세스는 실행되었지만, 세 번째는 실행되지 않았습니다.

        **오류 10: 출력 파일 이름 불일치**
        ```groovy linenums="88"
        file_ch = channel.fromPath("*.txt") // 오류: pwd가 아닌 프로세스로부터 입력을 가져오려고 함
        handleFiles(file_ch)
        ```
        **수정:** 이전 프로세스의 출력 사용
        ```groovy linenums="88"
        handleFiles(heavyProcess.out)
        ```

        이것으로 전체 워크플로우가 실행되어야 합니다.

        **완전 수정된 워크플로우:**
        ```groovy linenums="1"
        #!/usr/bin/env nextflow

        /*
        * 디버깅 연습을 위한 버그가 있는 워크플로우
        * 이 워크플로우는 학습 목적으로 여러 의도적인 버그를 포함하고 있습니다
        */

        params{
            // 검증이 누락된 매개변수
            input: Path = 'data/sample_data.csv'
            output: String = 'results'
        }

        /*
        * 입력/출력 불일치가 있는 프로세스
        */
        process processFiles {
            publishDir "${params.output}/processed", mode: 'copy'

            input:
                tuple val(sample_id), path(input_file)

            output:
                path "${sample_id}_result.txt"

            script:
            """
            echo "Processing: ${sample_id}"
            cat ${input_file} > ${sample_id}_result.txt
            """
        }

        /*
        * 리소스 문제가 있는 프로세스
        */
        process heavyProcess {
            publishDir "${params.output}/heavy", mode: 'copy'

            time '100 s'

            input:
                val sample_id

            output:
                path "${sample_id}_heavy.txt"

            script:
            """
            # 무거운 계산 시뮬레이션
            for i in {1..1000000}; do
                echo "Heavy computation \$i for ${sample_id}"
            done > ${sample_id}_heavy.txt
            """
        }

        /*
        * 파일 처리 문제가 있는 프로세스
        */
        process handleFiles {
            publishDir "${params.output}/files", mode: 'copy'

            input:
                path input_file

            output:
                path "processed_${input_file}"

            script:
            """
            if [ -f "${input_file}" ]; then
                cp ${input_file} processed_${input_file}
            fi
            """
        }

        /*
        * 채널 문제가 있는 메인 워크플로우
        */
        workflow {

            // 잘못된 사용법이 있는 채널
            input_ch = channel
                .fromPath(params.input)
                .splitCsv(header: true)
                .map { row -> [row.sample_id, file(row.fastq_path)] }

            processed_ch = processFiles(input_ch)

            heavy_ch = heavyProcess(input_ch.map{it[0]})

            handleFiles(heavyProcess.out)
        }
        ```

    **다루어진 오류 카테고리:**

    - **구문 오류**: 누락된 괄호, 후행 쉼표, 정의되지 않은 변수
    - **채널 구조 오류**: 잘못된 데이터 형태, 정의되지 않은 채널
    - **프로세스 오류**: 출력 파일 불일치, 변수 이스케이프
    - **리소스 오류**: 비현실적인 시간 제한

    **핵심 디버깅 교훈:**

    1. **오류 메시지를 주의 깊게 읽기** - 종종 문제를 직접 가리킵니다
    2. **체계적인 접근 방식 사용** - 한 번에 하나의 오류를 수정하고 `-resume`으로 테스트
    3. **데이터 흐름 이해** - 채널 구조 오류는 종종 가장 미묘합니다
    4. **작업 디렉토리 확인** - 프로세스가 실패하면 로그가 정확히 무엇이 잘못되었는지 알려줍니다

---

## 요약

이 사이드 퀘스트에서는 Nextflow 워크플로우를 디버깅하기 위한 일련의 체계적인 기법을 배웠습니다.
이러한 기법을 자신의 작업에 적용하면 컴퓨터와 싸우는 데 시간을 덜 쓰고, 문제를 더 빨리 해결하며, 미래의 문제로부터 자신을 보호할 수 있습니다.

### 주요 패턴

**1. 구문 오류를 식별하고 수정하는 방법**:

- Nextflow 오류 메시지 해석 및 문제 위치 파악
- 일반적인 구문 오류: 누락된 괄호, 잘못된 키워드, 정의되지 않은 변수
- Nextflow(Groovy)와 Bash 변수 구분
- 초기 오류 감지를 위한 VS Code 확장 기능 사용

```groovy
// 누락된 괄호 - IDE에서 빨간색 밑줄 확인
process FOO {
    script:
    """
    echo "hello"
    """
// } <-- 누락됨!

// 잘못된 키워드
inputs:  // 'input:'이어야 함

// 정의되지 않은 변수 - Bash 변수는 백슬래시로 이스케이프
echo "${undefined_var}"      // Nextflow 변수 (정의되지 않은 경우 오류)
echo "\${bash_var}"          // Bash 변수 (이스케이프됨)
```

**2. 채널 구조 문제를 디버깅하는 방법**:

- 채널 카디널리티 및 소진 문제 이해
- 채널 내용 구조 불일치 디버깅
- 채널 검사를 위한 `.view()` 연산자 사용
- 출력에 대괄호와 같은 오류 패턴 인식

```groovy
// 채널 내용 검사
my_channel.view { "Content: $it" }

// 큐를 값 채널로 변환(소진 방지)
reference_ch = channel.value('ref.fa')
// 또는
reference_ch = channel.of('ref.fa').first()
```

**3. 프로세스 실행 문제를 해결하는 방법**:

- 누락된 출력 파일 오류 진단
- 종료 코드 이해(소프트웨어 누락 시 127, 메모리 문제 시 137)
- 작업 디렉토리 및 명령 파일 조사
- 적절한 리소스 구성

```bash
# 실제로 실행된 내용 확인
cat work/ab/cdef12/.command.sh

# 오류 출력 확인
cat work/ab/cdef12/.command.err

# 종료 코드 127 = 명령을 찾을 수 없음
# 종료 코드 137 = 종료됨 (메모리/시간 제한)
```

**4. Nextflow의 내장 디버깅 도구를 사용하는 방법**:

- 미리보기 모드 및 실시간 디버깅 활용
- 로직 테스트를 위한 스텁 실행 구현
- 효율적인 디버깅 사이클을 위한 재개 적용
- 4단계 체계적인 디버깅 방법론 따르기

!!! tip "빠른 디버깅 참조"

    **구문 오류?** → VSCode 경고 확인, `nextflow run workflow.nf -preview` 실행

    **채널 문제?** → `.view()`를 사용하여 내용 검사: `my_channel.view()`

    **프로세스 실패?** → 작업 디렉토리 파일 확인:

    - `.command.sh` - 실행된 스크립트
    - `.command.err` - 오류 메시지
    - `.exitcode` - 종료 상태(127 = 명령을 찾을 수 없음, 137 = 종료됨)

    **미스터리한 동작?** → 워크플로우 로직을 테스트하기 위해 `-stub-run`으로 실행

    **수정 완료?** → 시간을 절약하기 위해 `-resume` 사용: `nextflow run workflow.nf -resume`

---

### 추가 자료

- [Nextflow 문제 해결 가이드](https://www.nextflow.io/docs/latest/troubleshooting.html): 공식 문제 해결 문서
- [Nextflow 채널 이해](https://www.nextflow.io/docs/latest/channel.html): 채널 유형 및 동작에 대한 심층 설명
- [프로세스 지시어 참조](https://www.nextflow.io/docs/latest/process.html#directives): 사용 가능한 모든 프로세스 구성 옵션
- [nf-test](https://www.nf-test.com/): Nextflow 파이프라인을 위한 테스트 프레임워크
- [Nextflow Slack 커뮤니티](https://www.nextflow.io/slack-invite.html): 커뮤니티에서 도움 받기

프로덕션 워크플로우의 경우 다음을 고려하세요:

- [Seqera Platform](https://seqera.io/platform/) 설정으로 대규모 모니터링 및 디버깅
- 재현 가능한 소프트웨어 환경을 위한 [Wave 컨테이너](https://seqera.io/wave/) 사용

**기억하세요:** 효과적인 디버깅은 연습으로 향상되는 기술입니다. 여기에서 습득한 체계적인 방법론과 포괄적인 툴킷은 Nextflow 개발 여정 전반에 걸쳐 도움이 될 것입니다.

---

## 다음 단계

[사이드 퀘스트 메뉴](./index.md)로 돌아가거나 페이지 오른쪽 하단의 버튼을 클릭하여 목록의 다음 주제로 이동하세요.
