# 워크플로우 디버깅

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } AI 지원 번역 - [자세히 알아보기 및 개선 제안](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

디버깅은 수 시간의 좌절을 줄여주고 더 효과적인 Nextflow 개발자가 되도록 도와주는 핵심 기술입니다. 경력 전반에 걸쳐, 특히 시작 단계에서는 워크플로우를 구축하고 유지보수하면서 버그를 마주하게 될 것입니다. 체계적인 디버깅 접근법을 배우면 문제를 빠르게 식별하고 해결하는 데 도움이 됩니다.

### 학습 목표

이 사이드 퀘스트에서는 Nextflow 워크플로우를 위한 **체계적인 디버깅 기법**을 탐구합니다:

- **구문 오류 디버깅**: IDE 기능과 Nextflow 오류 메시지를 효과적으로 사용하기
- **채널 디버깅**: 데이터 흐름 문제와 채널 구조 문제 진단하기
- **프로세스 디버깅**: 실행 실패와 리소스 문제 조사하기
- **내장 디버깅 도구**: Nextflow의 미리보기 모드, stub 실행, 작업 디렉토리 활용하기
- **체계적 접근법**: 효율적인 디버깅을 위한 4단계 방법론

마지막에는 좌절스러운 오류 메시지를 해결책으로 가는 명확한 로드맵으로 변환하는 강력한 디버깅 방법론을 갖추게 될 것입니다.

### 사전 요구사항

이 사이드 퀘스트를 시작하기 전에 다음을 완료해야 합니다:

- [Hello Nextflow](../hello_nextflow/README.md) 튜토리얼 또는 동등한 초급 과정을 완료했어야 합니다.
- 기본 Nextflow 개념과 메커니즘(프로세스, 채널, 연산자) 사용에 익숙해야 합니다.

**선택 사항:** [IDE Features for Nextflow Development](./ide_features.md) 사이드 퀘스트를 먼저 완료하는 것을 권장합니다.
이는 여기서 많이 사용할 디버깅을 지원하는 IDE 기능(구문 강조, 오류 감지 등)에 대한 포괄적인 내용을 다룹니다.

---

## 0. 시작하기

#### 교육 코드스페이스 열기

아직 하지 않았다면, [환경 설정](../envsetup/index.md)에 설명된 대로 교육 환경을 여는지 확인하십시오.

[![Open in GitHub Codespaces](https://github.com/codespaces/badge.svg)](https://codespaces.new/nextflow-io/training?quickstart=1&ref=master)

#### 프로젝트 디렉토리로 이동

이 튜토리얼의 파일들이 위치한 디렉토리로 이동하겠습니다.

```bash
cd side-quests/debugging
```

VSCode가 이 디렉토리에 집중하도록 설정할 수 있습니다:

```bash
code .
```

#### 자료 검토

연습에 사용할 다양한 유형의 버그가 있는 예제 워크플로우 세트를 찾을 수 있습니다:

??? abstract "디렉토리 내용"

    ```console
    .
    ├── bad_bash_var.nf
    ├── bad_channel_shape.nf
    ├── bad_channel_shape_viewed_debug.nf
    ├── bad_channel_shape_viewed.nf
    ├── bad_number_inputs.nf
    ├── badpractice_syntax.nf
    ├── bad_resources.nf
    ├── bad_syntax.nf
    ├── buggy_workflow.nf
    ├── data
    │   ├── sample_001.fastq.gz
    │   ├── sample_002.fastq.gz
    │   ├── sample_003.fastq.gz
    │   ├── sample_004.fastq.gz
    │   ├── sample_005.fastq.gz
    │   └── sample_data.csv
    ├── exhausted.nf
    ├── invalid_process.nf
    ├── missing_output.nf
    ├── missing_software.nf
    ├── missing_software_with_stub.nf
    ├── nextflow.config
    └── no_such_var.nf
    ```

이 파일들은 실제 개발에서 마주치게 될 일반적인 디버깅 시나리오를 나타냅니다.

#### 과제 검토

여러분의 과제는 각 워크플로우를 실행하고, 오류를 식별하고, 수정하는 것입니다.

각 버그가 있는 워크플로우에 대해:

1. **워크플로우 실행**하고 오류 관찰하기
2. **오류 메시지 분석**: Nextflow가 무엇을 알려주고 있는가?
3. **문제 위치 찾기**: 제공된 단서를 사용하여 코드에서 찾기
4. **버그 수정**하고 해결책이 작동하는지 확인하기
5. 다음 섹션으로 넘어가기 전에 **파일 재설정**하기 (`git checkout <filename>` 사용)

연습은 간단한 구문 오류에서 더 미묘한 런타임 문제로 진행됩니다.
해결책은 인라인으로 논의되지만, 앞으로 읽기 전에 각각을 직접 해결해 보십시오.

#### 준비 체크리스트

뛰어들 준비가 되었다고 생각하십니까?

- [ ] 이 과정의 목표와 사전 요구사항을 이해합니다
- [ ] 코드스페이스가 실행 중입니다
- [ ] 작업 디렉토리를 적절하게 설정했습니다
- [ ] 과제를 이해합니다

모든 항목을 체크할 수 있다면, 시작할 준비가 되었습니다.

---

## 1. 구문 오류

구문 오류는 Nextflow 코드를 작성할 때 마주치게 될 가장 일반적인 유형의 오류입니다. 코드가 Nextflow DSL의 예상 구문 규칙을 준수하지 않을 때 발생합니다. 이러한 오류는 워크플로우가 전혀 실행되지 않도록 하므로, 빠르게 식별하고 수정하는 방법을 배우는 것이 중요합니다.

### 1.1. 중괄호 누락

가장 일반적인 구문 오류 중 하나이며, 때로는 디버그하기 더 복잡한 것 중 하나는 **괄호가 누락되거나 일치하지 않는 경우**입니다.

실용적인 예제로 시작하겠습니다.

#### 파이프라인 실행

```bash
nextflow run bad_syntax.nf
```

??? failure "명령 출력"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `bad_syntax.nf` [stupefied_bhabha] DSL2 - revision: ca6327fad2

    Error bad_syntax.nf:24:1: Unexpected input: '<EOF>'

    ERROR ~ Script compilation failed

     -- Check '.nextflow.log' file for details
    ```

**구문 오류 메시지의 주요 요소:**

- **파일과 위치**: 오류가 포함된 파일과 줄/열을 보여줍니다 (`bad_syntax.nf:24:1`)
- **오류 설명**: 분석기가 예상하지 않은 것을 찾았다고 설명합니다 (`Unexpected input: '<EOF>'`)
- **EOF 표시**: `<EOF>` (End Of File) 메시지는 분석기가 여전히 더 많은 내용을 기대하면서 파일 끝에 도달했음을 나타냅니다 - 닫히지 않은 중괄호의 전형적인 징후

#### 코드 확인

이제 `bad_syntax.nf`를 검토하여 오류의 원인을 이해하겠습니다:

```groovy title="bad_syntax.nf" hl_lines="14" linenums="1"
#!/usr/bin/env nextflow

process PROCESS_FILES {
    input:
    val sample_name

    output:
    path "${sample_name}_output.txt"

    script:
    """
    echo "Processing ${sample_name}" > ${sample_name}_output.txt
    """
// process에 닫는 중괄호가 누락됨

workflow {

    // 입력 채널 생성
    input_ch = channel.of('sample1', 'sample2', 'sample3')

    // 입력 채널로 process 호출
    PROCESS_FILES(input_ch)
}
```

이 예제의 목적을 위해 오류가 어디에 있는지 보여주는 주석을 남겨두었습니다. Nextflow VSCode 확장도 일치하지 않는 중괄호를 빨간색으로 표시하고 파일의 조기 종료를 강조하여 무엇이 잘못되었는지에 대한 힌트를 제공할 것입니다:

![잘못된 구문](img/bad_syntax.png)

**괄호 오류 디버깅 전략:**

1. VS Code의 괄호 일치 사용 (괄호 옆에 커서 배치)
2. Problems 패널에서 괄호 관련 메시지 확인
3. 각 여는 `{`에 대응하는 닫는 `}`가 있는지 확인

#### 코드 수정

주석을 누락된 닫는 중괄호로 교체합니다:

=== "수정 후"

    ```groovy title="bad_syntax.nf" hl_lines="14" linenums="1"
    #!/usr/bin/env nextflow

    process PROCESS_FILES {
        input:
        val sample_name

        output:
        path "${sample_name}_output.txt"

        script:
        """
        echo "Processing ${sample_name}" > ${sample_name}_output.txt
        """
    }  // 누락된 닫는 중괄호 추가

    workflow {

        // 입력 채널 생성
        input_ch = channel.of('sample1', 'sample2', 'sample3')

        // 입력 채널로 process 호출
        PROCESS_FILES(input_ch)
    }
    ```

=== "수정 전"

    ```groovy title="bad_syntax.nf" hl_lines="14" linenums="1"
    #!/usr/bin/env nextflow

    process PROCESS_FILES {
        input:
        val sample_name

        output:
        path "${sample_name}_output.txt"

        script:
        """
        echo "Processing ${sample_name}" > ${sample_name}_output.txt
        """
    // process에 닫는 중괄호가 누락됨

    workflow {

        // 입력 채널 생성
        input_ch = channel.of('sample1', 'sample2', 'sample3')

        // 입력 채널로 process 호출
        PROCESS_FILES(input_ch)
    }
    ```

#### 파이프라인 실행

이제 워크플로우를 다시 실행하여 작동하는지 확인합니다:

```bash
nextflow run bad_syntax.nf
```

??? success "명령 출력"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `bad_syntax.nf` [insane_faggin] DSL2 - revision: 961938ee2b

    executor >  local (3)
    [48/cd7f54] PROCESS_FILES (1) | 3 of 3 ✔
    ```

### 1.2. 잘못된 프로세스 키워드 또는 지시문 사용

또 다른 일반적인 구문 오류는 **잘못된 프로세스 정의**입니다. 이는 필수 블록을 정의하는 것을 잊거나 프로세스 정의에서 잘못된 지시문을 사용할 때 발생할 수 있습니다.

#### 파이프라인 실행

```bash
nextflow run invalid_process.nf
```

??? failure "명령 출력"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `invalid_process.nf` [nasty_jepsen] DSL2 - revision: da9758d614

    Error invalid_process.nf:3:1: Invalid process definition -- check for missing or out-of-order section labels
    │   3 | process PROCESS_FILES {
    │     | ^^^^^^^^^^^^^^^^^^^^^^^
    │   4 |     inputs:
    │   5 |     val sample_name
    │   6 |
    ╰   7 |     output:

    ERROR ~ Script compilation failed

     -- Check '.nextflow.log' file for details
    ```

#### 코드 확인

오류는 "Invalid process definition"을 나타내고 문제 주변의 컨텍스트를 보여줍니다. 3-7줄을 보면 4번 줄에 `inputs:`가 있는데, 이것이 문제입니다. `invalid_process.nf`를 살펴보겠습니다:

```groovy title="invalid_process.nf" hl_lines="4" linenums="1"
#!/usr/bin/env nextflow

process PROCESS_FILES {
    inputs:  // 오류: 'inputs'가 아니라 'input'이어야 함
    val sample_name

    output:
    path "${sample_name}_output.txt"

    script:
    """
    echo "Processing ${sample_name}" > ${sample_name}_output.txt
    """
}

workflow {

    // 입력 채널 생성
    input_ch = channel.of('sample1', 'sample2', 'sample3')

    // 입력 채널로 process 호출
    PROCESS_FILES(input_ch)
}
```

오류 컨텍스트의 4번 줄을 보면 문제를 발견할 수 있습니다: 올바른 `input` 지시문 대신 `inputs`를 사용하고 있습니다. Nextflow VSCode 확장도 이를 표시할 것입니다:

![잘못된 process 메시지](img/invalid_process_message.png)

#### 코드 수정

[문서](https://www.nextflow.io/docs/latest/process.html#)를 참조하여 잘못된 키워드를 올바른 것으로 교체합니다:

=== "수정 후"

    ```groovy title="invalid_process.nf" hl_lines="4" linenums="1"
    #!/usr/bin/env nextflow

    process PROCESS_FILES {
        input:  // 수정됨: 'inputs'를 'input'으로 변경
        val sample_name

        output:
        path "${sample_name}_output.txt"

        script:
        """
        echo "Processing ${sample_name}" > ${sample_name}_output.txt
        """
    }

    workflow {

        // 입력 채널 생성
        input_ch = channel.of('sample1', 'sample2', 'sample3')

        // 입력 채널로 process 호출
        PROCESS_FILES(input_ch)
    }
    ```

=== "수정 전"

    ```groovy title="invalid_process.nf" hl_lines="4" linenums="1"
    #!/usr/bin/env nextflow

    process PROCESS_FILES {
        inputs:  // 오류: 'inputs'가 아니라 'input'이어야 함
        val sample_name

        output:
        path "${sample_name}_output.txt"

        script:
        """
        echo "Processing ${sample_name}" > ${sample_name}_output.txt
        """
    }

    workflow {

        // 입력 채널 생성
        input_ch = channel.of('sample1', 'sample2', 'sample3')

        // 입력 채널로 process 호출
        PROCESS_FILES(input_ch)
    }
    ```

#### 파이프라인 실행

이제 워크플로우를 다시 실행하여 작동하는지 확인합니다:

```bash
nextflow run invalid_process.nf
```

??? success "명령 출력"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `invalid_process.nf` [silly_fermi] DSL2 - revision: 961938ee2b

    executor >  local (3)
    [b7/76cd9d] PROCESS_FILES (2) | 3 of 3 ✔
    ```

### 1.3. 잘못된 변수 이름 사용

스크립트 블록에서 사용하는 변수 이름은 입력에서 파생되거나 스크립트 전에 삽입된 groovy 코드에서 파생된 유효한 것이어야 합니다. 하지만 파이프라인 개발 초기에 복잡성을 다루다 보면 변수 이름에서 실수하기 쉽고, Nextflow가 빠르게 알려줄 것입니다.

#### 파이프라인 실행

```bash
nextflow run no_such_var.nf
```

??? failure "명령 출력"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `no_such_var.nf` [gloomy_meninsky] DSL2 - revision: 0c4d3bc28c

    Error no_such_var.nf:17:39: `undefined_var` is not defined
    │  17 |     echo "Using undefined variable: ${undefined_var}" >> ${output_pref
    ╰     |                                       ^^^^^^^^^^^^^

    ERROR ~ Script compilation failed

     -- Check '.nextflow.log' file for details
    ```

오류는 컴파일 시간에 포착되고 17번 줄의 정의되지 않은 변수를 직접 가리키며, 캐럿이 문제가 정확히 어디에 있는지 표시합니다.

#### 코드 확인

`no_such_var.nf`를 살펴보겠습니다:

```groovy title="no_such_var.nf" hl_lines="17" linenums="1"
#!/usr/bin/env nextflow

process PROCESS_FILES {
    input:
    val sample_name

    output:
    path "${sample_name}_processed.txt"

    script:
    // 스크립트 전에 Groovy 코드에서 변수 정의
    def output_prefix = "${sample_name}_processed"
    def timestamp = new Date().format("yyyy-MM-dd")

    """
    echo "Processing ${sample_name} on ${timestamp}" > ${output_prefix}.txt
    echo "Using undefined variable: ${undefined_var}" >> ${output_prefix}.txt  // 오류: undefined_var가 정의되지 않음
    """
}

workflow {
    input_ch = channel.of('sample1', 'sample2', 'sample3')
    PROCESS_FILES(input_ch)
}
```

오류 메시지는 변수가 스크립트 템플릿에서 인식되지 않음을 나타내며, 스크립트 블록에서 사용되지만 다른 곳에서 정의되지 않은 `${undefined_var}`를 볼 수 있어야 합니다.

#### 코드 수정

'No such variable' 오류가 발생하면 변수를 정의하거나(입력 변수 이름을 수정하거나 스크립트 전에 groovy 코드를 편집하여) 필요하지 않은 경우 스크립트 블록에서 제거하여 수정할 수 있습니다:

=== "수정 후"

    ```groovy title="no_such_var.nf" hl_lines="15-17" linenums="1"
    #!/usr/bin/env nextflow

    process PROCESS_FILES {
        input:
        val sample_name

        output:
        path "${sample_name}_output.txt"

        script:
        // 스크립트 전에 Groovy 코드에서 변수 정의
        def output_prefix = "${sample_name}_processed"
        def timestamp = new Date().format("yyyy-MM-dd")

        """
        echo "Processing ${sample_name} on ${timestamp}" > ${output_prefix}.txt
        """  // undefined_var가 있는 줄 제거됨
    }

    workflow {
        input_ch = channel.of('sample1', 'sample2', 'sample3')
        PROCESS_FILES(input_ch)
    }
    ```

=== "수정 전"

    ```groovy title="no_such_var.nf" hl_lines="17" linenums="1"
    #!/usr/bin/env nextflow

    process PROCESS_FILES {
        input:
        val sample_name

        output:
        path "${sample_name}_output.txt"

        script:
        // 스크립트 전에 Groovy 코드에서 변수 정의
        def output_prefix = "${sample_name}_processed"
        def timestamp = new Date().format("yyyy-MM-dd")

        """
        echo "Processing ${sample_name} on ${timestamp}" > ${output_prefix}.txt
        echo "Using undefined variable: ${undefined_var}" >> ${output_prefix}.txt  // 오류: undefined_var가 정의되지 않음
        """
    }

    workflow {
        input_ch = channel.of('sample1', 'sample2', 'sample3')
        PROCESS_FILES(input_ch)
    }
    ```

#### 파이프라인 실행

이제 워크플로우를 다시 실행하여 작동하는지 확인합니다:

```bash
nextflow run no_such_var.nf
```

??? success "명령 출력"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `no_such_var.nf` [suspicious_venter] DSL2 - revision: 6ba490f7c5

    executor >  local (3)
    [21/237300] PROCESS_FILES (2) | 3 of 3 ✔
    ```

### 1.4. Bash 변수의 잘못된 사용

Nextflow를 시작할 때 Nextflow(Groovy)와 Bash 변수의 차이를 이해하기 어려울 수 있습니다. 이는 스크립트 블록의 Bash 내용에서 변수를 사용하려고 할 때 나타나는 또 다른 형태의 잘못된 변수 오류를 생성할 수 있습니다.

#### 파이프라인 실행

```bash
nextflow run bad_bash_var.nf
```

??? failure "명령 출력"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `bad_bash_var.nf` [infallible_mandelbrot] DSL2 - revision: 0853c11080

    Error bad_bash_var.nf:13:42: `prefix` is not defined
    │  13 |     echo "Processing ${sample_name}" > ${prefix}.txt
    ╰     |                                          ^^^^^^

    ERROR ~ Script compilation failed

     -- Check '.nextflow.log' file for details
    ```

#### 코드 확인

오류는 `${prefix}`가 사용된 13번 줄을 가리킵니다. `bad_bash_var.nf`를 살펴보고 무엇이 문제를 일으키는지 확인하겠습니다:

```groovy title="bad_bash_var.nf" hl_lines="13" linenums="1"
#!/usr/bin/env nextflow

process PROCESS_FILES {
    input:
    val sample_name

    output:
    path "${sample_name}_output.txt"

    script:
    """
    prefix="${sample_name}_output"
    echo "Processing ${sample_name}" > ${prefix}.txt  # 오류: ${prefix}는 Bash가 아닌 Groovy 구문임
    """
}
```

이 예제에서 Bash에서 `prefix` 변수를 정의하고 있지만, Nextflow 프로세스에서 이를 참조하기 위해 사용한 `$` 구문(`${prefix}`)은 Bash가 아닌 Groovy 변수로 해석됩니다. 변수가 Groovy 컨텍스트에 존재하지 않으므로 'no such variable' 오류가 발생합니다.

#### 코드 수정

Bash 변수를 사용하려면 다음과 같이 달러 기호를 이스케이프해야 합니다:

=== "수정 후"

    ```groovy title="bad_bash_var.nf" hl_lines="13" linenums="1"
    #!/usr/bin/env nextflow

    process PROCESS_FILES {
        input:
        val sample_name

        output:
        path "${sample_name}_output.txt"

        script:
        """
        prefix="${sample_name}_output"
        echo "Processing ${sample_name}" > \${prefix}.txt  # 수정됨: 달러 기호를 이스케이프함
        """
    }

    workflow {
        input_ch = channel.of('sample1', 'sample2', 'sample3')
        PROCESS_FILES(input_ch)
    }
    ```

=== "수정 전"

    ```groovy title="bad_bash_var.nf" hl_lines="13" linenums="1"
    #!/usr/bin/env nextflow

    process PROCESS_FILES {
        input:
        val sample_name

        output:
        path "${sample_name}_output.txt"

        script:
        """
        prefix="${sample_name}_output"
        echo "Processing ${sample_name}" > ${prefix}.txt  # 오류: ${prefix}는 Bash가 아닌 Groovy 구문임
        """
    }
    ```

이는 Nextflow에게 이것을 Bash 변수로 해석하도록 알려줍니다.

#### 파이프라인 실행

이제 워크플로우를 다시 실행하여 작동하는지 확인합니다:

```bash
nextflow run bad_bash_var.nf
```

??? success "명령 출력"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `bad_bash_var.nf` [naughty_franklin] DSL2 - revision: 58c1c83709

    executor >  local (3)
    [4e/560285] PROCESS_FILES (2) | 3 of 3 ✔
    ```

!!! tip "Groovy vs Bash 변수"

    문자열 연결이나 접두사/접미사 작업과 같은 간단한 변수 조작의 경우, 스크립트 블록의 Bash 변수보다 스크립트 섹션의 Groovy 변수를 사용하는 것이 일반적으로 더 읽기 쉽습니다:

    ```groovy linenums="1"
    script:
    def output_prefix = "${sample_name}_processed"
    def output_file = "${output_prefix}.txt"
    """
    echo "Processing ${sample_name}" > ${output_file}
    """
    ```

    이 접근 방식은 달러 기호를 이스케이프할 필요를 피하고 코드를 더 읽기 쉽고 유지보수하기 쉽게 만듭니다.

### 1.5. Workflow 블록 외부의 문장

Nextflow VSCode 확장은 오류를 일으킬 코드 구조의 문제를 강조합니다. 일반적인 예는 `workflow {}` 블록 외부에서 채널을 정의하는 것입니다 - 이는 이제 구문 오류로 시행됩니다.

#### 파이프라인 실행

```bash
nextflow run badpractice_syntax.nf
```

??? failure "명령 출력"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `badpractice_syntax.nf` [intergalactic_colden] DSL2 - revision: 5e4b291bde

    Error badpractice_syntax.nf:3:1: Statements cannot be mixed with script declarations -- move statements into a process or workflow
    │   3 | input_ch = channel.of('sample1', 'sample2', 'sample3')
    ╰     | ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

    ERROR ~ Script compilation failed

     -- Check '.nextflow.log' file for details
    ```

오류 메시지는 문제를 명확하게 나타냅니다: 워크플로우나 프로세스 블록 외부에서 스크립트 선언과 문장(채널 정의와 같은)을 혼합할 수 없습니다.

#### 코드 확인

`badpractice_syntax.nf`를 살펴보고 오류의 원인을 확인하겠습니다:

```groovy title="badpractice_syntax.nf" hl_lines="3" linenums="1"
#!/usr/bin/env nextflow

input_ch = channel.of('sample1', 'sample2', 'sample3')  // 오류: workflow 외부에서 정의된 채널

process PROCESS_FILES {
    input:
    val sample_name

    output:
    path "${sample_name}_processed.txt"

    script:
    // 스크립트 전에 Groovy 코드에서 변수 정의
    def output_prefix = "${sample_name}_processed"
    def timestamp = new Date().format("yyyy-MM-dd")

    """
    echo "Processing ${sample_name} on ${timestamp}" > ${output_prefix}.txt
    """
}

workflow {
    PROCESS_FILES(input_ch)
}
```

VSCode 확장은 또한 워크플로우 블록 외부에서 정의된 `input_ch` 변수를 강조할 것입니다:

![치명적이지 않은 구문 오류](img/nonlethal.png)

#### 코드 수정

채널 정의를 워크플로우 블록 내부로 이동합니다:

=== "수정 후"

    ```groovy title="badpractice_syntax.nf" hl_lines="21" linenums="1"
    #!/usr/bin/env nextflow

    process PROCESS_FILES {
        input:
        val sample_name

        output:
        path "${sample_name}_processed.txt"

        script:
        // 스크립트 전에 Groovy 코드에서 변수 정의
        def output_prefix = "${sample_name}_processed"
        def timestamp = new Date().format("yyyy-MM-dd")

        """
        echo "Processing ${sample_name} on ${timestamp}" > ${output_prefix}.txt
        """
    }

    workflow {
        input_ch = channel.of('sample1', 'sample2', 'sample3')  // workflow 블록 내부로 이동됨
        PROCESS_FILES(input_ch)
    }
    ```

=== "수정 전"

    ```groovy title="badpractice_syntax.nf" hl_lines="3" linenums="1"
    #!/usr/bin/env nextflow

    input_ch = channel.of('sample1', 'sample2', 'sample3')  // 오류: workflow 외부에서 정의된 채널

    process PROCESS_FILES {
        input:
        val sample_name

        output:
        path "${sample_name}_processed.txt"

        script:
        // 스크립트 전에 Groovy 코드에서 변수 정의
        def output_prefix = "${sample_name}_processed"
        def timestamp = new Date().format("yyyy-MM-dd")

        """
        echo "Processing ${sample_name} on ${timestamp}" > ${output_prefix}.txt
        """
    }

    workflow {
        PROCESS_FILES(input_ch)
    }
    ```

#### 파이프라인 실행

워크플로우를 다시 실행하여 수정이 작동하는지 확인합니다:

```bash
nextflow run badpractice_syntax.nf
```

??? success "명령 출력"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `badpractice_syntax.nf` [naughty_ochoa] DSL2 - revision: 5e4b291bde

    executor >  local (3)
    [6a/84a608] PROCESS_FILES (2) | 3 of 3 ✔
    ```

입력 채널을 워크플로우 블록 내에 정의하고, 일반적으로 확장이 제공하는 다른 권장 사항을 따르십시오.

### 핵심 요약

Nextflow 오류 메시지와 IDE 시각적 표시를 사용하여 구문 오류를 체계적으로 식별하고 수정할 수 있습니다. 일반적인 구문 오류에는 중괄호 누락, 잘못된 프로세스 키워드, 정의되지 않은 변수, Bash와 Nextflow 변수의 부적절한 사용이 포함됩니다. VSCode 확장은 런타임 전에 이러한 많은 것들을 포착하는 데 도움이 됩니다. 이러한 구문 디버깅 기술을 도구 상자에 넣으면, 가장 일반적인 Nextflow 구문 오류를 빠르게 해결하고 더 복잡한 런타임 문제를 다루는 데 집중할 수 있습니다.

### 다음은?

구문이 올바른 경우에도 발생하는 더 복잡한 채널 구조 오류를 디버그하는 방법을 배웁니다.

---

## 2. 채널 구조 오류

채널 구조 오류는 코드가 구문적으로 올바르지만 데이터 형태가 프로세스가 예상하는 것과 일치하지 않기 때문에 구문 오류보다 더 미묘합니다. Nextflow는 파이프라인을 실행하려고 시도하지만, 입력 수가 예상과 일치하지 않음을 발견하고 실패할 수 있습니다. 이러한 오류는 일반적으로 런타임에만 나타나며 워크플로우를 통해 흐르는 데이터에 대한 이해가 필요합니다.

!!! tip "`.view()`로 채널 디버깅"

    이 섹션 전체에서 `.view()` 연산자를 사용하여 워크플로우의 어느 지점에서든 채널 내용을 검사할 수 있다는 것을 기억하십시오. 이것은 채널 구조 문제를 이해하기 위한 가장 강력한 디버깅 도구 중 하나입니다. 섹션 2.4에서 이 기법을 자세히 탐구할 것이지만, 예제를 작업하면서 자유롭게 사용하십시오.

    ```groovy
    my_channel.view()  // 채널을 통해 흐르는 내용을 표시
    ```

### 2.1. 잘못된 입력 채널 수

이 오류는 프로세스가 예상하는 것과 다른 수의 채널을 전달할 때 발생합니다.

#### 파이프라인 실행

```bash
nextflow run bad_number_inputs.nf
```

??? failure "명령 출력"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `bad_number_inputs.nf` [happy_swartz] DSL2 - revision: d83e58dcd3

    Error bad_number_inputs.nf:23:5: Incorrect number of call arguments, expected 1 but received 2
    │  23 |     PROCESS_FILES(samples_ch, files_ch)
    ╰     |     ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

    ERROR ~ Script compilation failed

     -- Check '.nextflow.log' file for details
    ```

#### 코드 확인

오류 메시지는 호출이 1개의 인수를 예상했지만 2개를 받았다고 명확히 설명하며 23번 줄을 가리킵니다. `bad_number_inputs.nf`를 살펴보겠습니다:

```groovy title="bad_number_inputs.nf" hl_lines="5 23" linenums="1"
#!/usr/bin/env nextflow

process PROCESS_FILES {
    input:
        val sample_name  // process는 1개의 입력만 예상함

    output:
        path "${sample_name}_output.txt"

    script:
    """
    echo "Processing ${sample_name}" > ${sample_name}_output.txt
    """
}

workflow {

    // 두 개의 별도 채널 생성
    samples_ch = channel.of('sample1', 'sample2', 'sample3')
    files_ch = channel.of('file1.txt', 'file2.txt', 'file3.txt')

    // 오류: 2개의 채널을 전달하지만 process는 1개만 예상함
    PROCESS_FILES(samples_ch, files_ch)
}
```

프로세스가 하나만 정의하는데 여러 입력 채널을 제공하는 일치하지 않는 `PROCESS_FILES` 호출을 볼 수 있어야 합니다. VSCode 확장도 프로세스 호출 아래에 빨간 줄을 표시하고 마우스를 올리면 진단 메시지를 제공할 것입니다:

![잘못된 인수 개수 메시지](img/incorrect_num_args.png)

#### 코드 수정

이 특정 예제의 경우, 프로세스는 단일 채널을 예상하고 두 번째 채널이 필요하지 않으므로 `samples_ch` 채널만 전달하여 수정할 수 있습니다:

=== "수정 후"

    ```groovy title="bad_number_inputs.nf" hl_lines="23" linenums="1"
    #!/usr/bin/env nextflow

    process PROCESS_FILES {
        input:
            val sample_name  // process는 1개의 입력만 예상함

        output:
            path "${sample_name}_output.txt"

        script:
        """
        echo "Processing ${sample_name}" > ${sample_name}_output.txt
        """
    }

    workflow {

        // 두 개의 별도 채널 생성
        samples_ch = channel.of('sample1', 'sample2', 'sample3')
        files_ch = channel.of('file1.txt', 'file2.txt', 'file3.txt')

        // 수정됨: process가 예상하는 채널만 전달
        PROCESS_FILES(samples_ch)
    }
    ```

=== "수정 전"

    ```groovy title="bad_number_inputs.nf" hl_lines="5 23" linenums="1"
    #!/usr/bin/env nextflow

    process PROCESS_FILES {
        input:
            val sample_name  // process는 1개의 입력만 예상함

        output:
            path "${sample_name}_output.txt"

        script:
        """
        echo "Processing ${sample_name}" > ${sample_name}_output.txt
        """
    }

    workflow {

        // 두 개의 별도 채널 생성
        samples_ch = channel.of('sample1', 'sample2', 'sample3')
        files_ch = channel.of('file1.txt', 'file2.txt', 'file3.txt')

        // 오류: 2개의 채널을 전달하지만 process는 1개만 예상함
        PROCESS_FILES(samples_ch, files_ch)
    }
    ```

#### 파이프라인 실행

```bash
nextflow run bad_number_inputs.nf
```

??? success "명령 출력"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `bad_number_inputs.nf` [big_euler] DSL2 - revision: e302bd87be

    executor >  local (3)
    [48/497f7b] PROCESS_FILES (3) | 3 of 3 ✔
    ```

이 예제보다 더 일반적으로, 프로세스에 추가 입력을 추가하고 워크플로우 호출을 그에 따라 업데이트하는 것을 잊어버릴 수 있으며, 이는 이러한 유형의 오류로 이어질 수 있습니다. 다행히도 오류 메시지가 불일치에 대해 매우 명확하므로 이것은 이해하고 수정하기 더 쉬운 오류 중 하나입니다.

### 2.2. 채널 소진 (프로세스가 예상보다 적게 실행됨)

일부 채널 구조 오류는 훨씬 더 미묘하고 전혀 오류를 생성하지 않습니다. 아마도 가장 일반적인 것은 queue 채널이 소진되어 항목이 부족할 수 있다는 것을 이해하는 데 새로운 Nextflow 사용자가 직면하는 문제를 반영하며, 워크플로우가 조기에 종료됩니다.

#### 파이프라인 실행

```bash
nextflow run exhausted.nf
```

??? success "명령 출력"

```console title="Exhausted channel output"
 N E X T F L O W   ~  version 25.10.2

Launching `exhausted.nf` [extravagant_gauss] DSL2 - revision: 08cff7ba2a

executor >  local (1)
[bd/f61fff] PROCESS_FILES (1) [100%] 1 of 1 ✔
```

이 워크플로우는 오류 없이 완료되지만, 단일 샘플만 처리합니다!

#### 코드 확인

`exhausted.nf`를 살펴보고 맞는지 확인하겠습니다:

```groovy title="exhausted.nf" hl_lines="23 24" linenums="1"
#!/usr/bin/env nextflow

process PROCESS_FILES {
    input:
    val reference
    val sample_name

    output:
    path "${output_prefix}.txt"

    script:
    // 스크립트 전에 Groovy 코드에서 변수 정의
    output_prefix = "${reference}_${sample_name}"
    def timestamp = new Date().format("yyyy-MM-dd")

    """
    echo "Processing ${sample_name} on ${timestamp}" > ${output_prefix}.txt
    """
}

workflow {

    reference_ch = channel.of('baseline_reference')
    input_ch = channel.of('sample1', 'sample2', 'sample3')

    PROCESS_FILES(reference_ch, input_ch)
}
```

프로세스가 세 번 대신 한 번만 실행되는 이유는 `reference_ch` 채널이 첫 번째 프로세스 실행 후 소진되는 queue 채널이기 때문입니다. 한 채널이 소진되면, 다른 채널에 여전히 항목이 있더라도 전체 프로세스가 중지됩니다.

이것은 여러 샘플에서 재사용해야 하는 단일 참조 파일이 있는 일반적인 패턴입니다. 해결책은 참조 채널을 무기한 재사용할 수 있는 value 채널로 변환하는 것입니다.

#### 코드 수정

영향을 받는 파일 수에 따라 이를 해결하는 몇 가지 방법이 있습니다.

**옵션 1**: 많이 재사용하는 단일 참조 파일이 있습니다. 반복해서 사용할 수 있는 value 채널 타입을 간단히 생성할 수 있습니다. 이를 수행하는 세 가지 방법이 있습니다:

**1a** `channel.value()` 사용:

```groovy title="exhausted.nf (fixed - Option 1a)" hl_lines="2" linenums="21"
workflow {
    reference_ch = channel.value('baseline_reference')  // value 채널은 재사용 가능
    input_ch = channel.of('sample1', 'sample2', 'sample3')

    PROCESS_FILES(reference_ch, input_ch)
}
```

**1b** `first()` [연산자](https://www.nextflow.io/docs/latest/reference/operator.html#first) 사용:

```groovy title="exhausted.nf (fixed - Option 1b)" hl_lines="2" linenums="21"
workflow {
    reference_ch = channel.of('baseline_reference').first()  // value 채널로 변환
    input_ch = channel.of('sample1', 'sample2', 'sample3')

    PROCESS_FILES(reference_ch, input_ch)
}
```

**1c.** `collect()` [연산자](https://www.nextflow.io/docs/latest/reference/operator.html#collect) 사용:

```groovy title="exhausted.nf (fixed - Option 1c)" hl_lines="2" linenums="21"
workflow {
    reference_ch = channel.of('baseline_reference').collect()  // value 채널로 변환
    input_ch = channel.of('sample1', 'sample2', 'sample3')

    PROCESS_FILES(reference_ch, input_ch)
}
```

**옵션 2**: 더 복잡한 시나리오에서, 아마도 샘플 채널의 모든 샘플에 대한 여러 참조 파일이 있는 경우, `combine` 연산자를 사용하여 두 채널을 튜플로 결합하는 새 채널을 생성할 수 있습니다:

```groovy title="exhausted.nf (fixed - Option 2)" hl_lines="4" linenums="21"
workflow {
    reference_ch = channel.of('baseline_reference','other_reference')
    input_ch = channel.of('sample1', 'sample2', 'sample3')
    combined_ch = reference_ch.combine(input_ch)  // 데카르트 곱 생성

    PROCESS_FILES(combined_ch)
}
```

`.combine()` 연산자는 두 채널의 데카르트 곱을 생성하므로 `reference_ch`의 각 항목이 `input_ch`의 각 항목과 쌍을 이룹니다. 이를 통해 프로세스가 참조를 계속 사용하면서 각 샘플에 대해 실행될 수 있습니다.

이를 위해서는 프로세스 입력을 조정해야 합니다. 우리 예제에서 프로세스 정의의 시작 부분은 다음과 같이 조정해야 합니다:

```groovy title="exhausted.nf (fixed - Option 2)" hl_lines="5" linenums="1"
#!/usr/bin/env nextflow

process PROCESS_FILES {
    input:
        tuple val(reference), val(sample_name)
```

이 접근 방식은 모든 상황에서 적합하지 않을 수 있습니다.

#### 파이프라인 실행

위의 수정 사항 중 하나를 시도하고 워크플로우를 다시 실행합니다:

```bash
nextflow run exhausted.nf
```

??? success "명령 출력"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `exhausted.nf` [maniac_leavitt] DSL2 - revision: f372a56a7d

    executor >  local (3)
    [80/0779e9] PROCESS_FILES (3) | 3 of 3 ✔
    ```

이제 하나 대신 세 개의 샘플이 모두 처리되는 것을 볼 수 있어야 합니다.

### 2.3. 잘못된 채널 내용 구조

워크플로우가 특정 수준의 복잡성에 도달하면 각 채널의 내부 구조를 추적하기가 조금 어려울 수 있으며, 사람들은 일반적으로 프로세스가 예상하는 것과 채널이 실제로 포함하는 것 사이에 불일치를 생성합니다. 이것은 채널 수가 잘못되었던 앞서 논의한 문제보다 더 미묘합니다. 이 경우 올바른 수의 입력 채널을 가질 수 있지만, 하나 이상의 채널의 내부 구조가 프로세스가 예상하는 것과 일치하지 않습니다.

#### 파이프라인 실행

```bash
nextflow run bad_channel_shape.nf
```

??? failure "명령 출력"

    ```console
    Launching `bad_channel_shape.nf` [hopeful_pare] DSL2 - revision: ffd66071a1

    executor >  local (3)
    executor >  local (3)
    [3f/c2dcb3] PROCESS_FILES (3) [  0%] 0 of 3 ✘
    ERROR ~ Error executing process > 'PROCESS_FILES (1)'

    Caused by:
      Missing output file(s) `[sample1, file1.txt]_output.txt` expected by process `PROCESS_FILES (1)`


    Command executed:

      echo "Processing [sample1, file1.txt]" > [sample1, file1.txt]_output.txt

    Command exit status:
      0

    Command output:
      (empty)

    Work dir:
      /workspaces/training/side-quests/debugging/work/d6/1fb69d1d93300bbc9d42f1875b981e

    Tip: when you have fixed the problem you can continue the execution adding the option `-resume` to the run command line

    -- Check '.nextflow.log' file for details
    ```

#### 코드 확인

오류 메시지의 대괄호가 여기서 단서를 제공합니다 - 프로세스가 튜플을 단일 값으로 처리하고 있으며, 이것은 우리가 원하는 것이 아닙니다. `bad_channel_shape.nf`를 살펴보겠습니다:

```groovy title="bad_channel_shape.nf" hl_lines="5 20-22" linenums="1"
#!/usr/bin/env nextflow

process PROCESS_FILES {
    input:
        val sample_name  // 단일 값을 예상하지만 튜플을 받음

    output:
        path "${sample_name}_output.txt"

    script:
    """
    echo "Processing ${sample_name}" > ${sample_name}_output.txt
    """
}

workflow {

    // 채널이 튜플을 내보내지만 process는 단일 값을 예상함
    input_ch = channel.of(
      ['sample1', 'file1.txt'],
      ['sample2', 'file2.txt'],
      ['sample3', 'file3.txt']
    )
    PROCESS_FILES(input_ch)
}
```

튜플로 구성된 채널을 생성하고 있음을 볼 수 있습니다: `['sample1', 'file1.txt']`, 하지만 프로세스는 단일 값 `val sample_name`을 예상합니다. 실행된 명령은 프로세스가 `[sample3, file3.txt]_output.txt`라는 파일을 생성하려고 시도하고 있음을 보여주며, 이는 의도한 출력이 아닙니다.

#### 코드 수정

이를 수정하려면 프로세스가 두 입력이 모두 필요한 경우 튜플을 수락하도록 프로세스를 조정할 수 있습니다:

=== "옵션 1: 프로세스에서 튜플 수락"

    === "수정 후"

        ```groovy title="bad_channel_shape.nf" hl_lines="5"  linenums="1"
        #!/usr/bin/env nextflow

        process PROCESS_FILES {
            input:
                tuple val(sample_name), val(file_name)  // 수정됨: 튜플 수락

            output:
                path "${sample_name}_output.txt"

            script:
            """
            echo "Processing ${sample_name}" > ${sample_name}_output.txt
            """
        }

        workflow {

            // 채널이 튜플을 내보내지만 process는 단일 값을 예상함
            input_ch = channel.of(
              ['sample1', 'file1.txt'],
              ['sample2', 'file2.txt'],
              ['sample3', 'file3.txt']
            )
            PROCESS_FILES(input_ch)
        }
        ```

    === "수정 전"

        ```groovy title="bad_channel_shape.nf" hl_lines="5" linenums="1"
        #!/usr/bin/env nextflow

        process PROCESS_FILES {
            input:
                val sample_name  // 단일 값을 예상하지만 튜플을 받음

            output:
                path "${sample_name}_output.txt"

            script:
            """
            echo "Processing ${sample_name}" > ${sample_name}_output.txt
            """
        }

        workflow {

            // 채널이 튜플을 내보내지만 process는 단일 값을 예상함
            input_ch = channel.of(
              ['sample1', 'file1.txt'],
              ['sample2', 'file2.txt'],
              ['sample3', 'file3.txt']
            )
            PROCESS_FILES(input_ch)
        }
        ```

=== "옵션 2: 첫 번째 요소 추출"

    === "수정 후"

        ```groovy title="bad_channel_shape.nf" hl_lines="9" linenums="16"
        workflow {

            // 채널이 튜플을 내보내지만 process는 단일 값을 예상함
            input_ch = channel.of(
              ['sample1', 'file1.txt'],
              ['sample2', 'file2.txt'],
              ['sample3', 'file3.txt']
            )
            PROCESS_FILES(input_ch.map { it[0] })  // 수정됨: 첫 번째 요소 추출
        }
        ```

    === "수정 전"

        ```groovy title="bad_channel_shape.nf" hl_lines="9" linenums="16"
        workflow {

            // 채널이 튜플을 내보내지만 process는 단일 값을 예상함
            input_ch = channel.of(
              ['sample1', 'file1.txt'],
              ['sample2', 'file2.txt'],
              ['sample3', 'file3.txt']
            )
            PROCESS_FILES(input_ch)
        }
        ```

#### 파이프라인 실행

해결책 중 하나를 선택하고 워크플로우를 다시 실행합니다:

```bash
nextflow run bad_channel_shape.nf
```

??? success "명령 출력"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `bad_channel_shape.nf` [clever_thompson] DSL2 - revision: 8cbcae3746

    executor >  local (3)
    [bb/80a958] PROCESS_FILES (2) | 3 of 3 ✔
    ```

### 2.4. 채널 디버깅 기법

#### 채널 검사를 위한 `.view()` 사용

채널을 위한 가장 강력한 디버깅 도구는 `.view()` 연산자입니다. `.view()`를 사용하면 디버깅에 도움이 되도록 모든 단계에서 채널의 형태를 이해할 수 있습니다.

#### 파이프라인 실행

`bad_channel_shape_viewed.nf`를 실행하여 이를 실제로 확인합니다:

```bash
nextflow run bad_channel_shape_viewed.nf
```

??? success "명령 출력"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `bad_channel_shape_viewed.nf` [maniac_poisson] DSL2 - revision: b4f24dc9da

    executor >  local (3)
    [c0/db76b3] PROCESS_FILES (3) [100%] 3 of 3 ✔
    Channel content: [sample1, file1.txt]
    Channel content: [sample2, file2.txt]
    Channel content: [sample3, file3.txt]
    After mapping: sample1
    After mapping: sample2
    After mapping: sample3
    ```

#### 코드 확인

`bad_channel_shape_viewed.nf`를 살펴보고 `.view()`가 어떻게 사용되는지 확인하겠습니다:

```groovy title="bad_channel_shape_viewed.nf" linenums="16" hl_lines="9 11"
workflow {

    // 채널이 튜플을 내보내지만 process는 단일 값을 예상함
    input_ch = channel.of(
      ['sample1', 'file1.txt'],
      ['sample2', 'file2.txt'],
      ['sample3', 'file3.txt']
    )
    .view { "Channel content: $it" }  // 디버그: 원래 채널 내용 표시
    .map { tuple -> tuple[0] }        // 변환: 첫 번째 요소 추출
    .view { "After mapping: $it" }    // 디버그: 변환된 채널 내용 표시

    PROCESS_FILES(input_ch)
}
```

#### 코드 수정

미래에 채널 내용을 이해하기 위해 과도하게 `.view()` 연산을 사용하는 것을 피하려면 도움이 되는 주석을 추가하는 것이 좋습니다:

```groovy title="bad_channel_shape_viewed.nf (with comments)" linenums="16" hl_lines="8 9"
workflow {

    // 채널이 튜플을 내보내지만 process는 단일 값을 예상함
    input_ch = channel.of(
            ['sample1', 'file1.txt'],
            ['sample2', 'file2.txt'],
            ['sample3', 'file3.txt'],
        ) // [sample_name, file_name]
        .map { tuple -> tuple[0] } // sample_name

    PROCESS_FILES(input_ch)
}
```

이것은 워크플로우가 복잡성이 증가하고
