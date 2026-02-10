# 워크플로우 디버깅

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } AI 지원 번역 - [자세히 알아보기 및 개선 제안](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

디버깅은 여러분의 시간을 절약하고 더 효과적인 Nextflow 개발자가 되도록 도와주는 중요한 기술입니다. 경력 전반에 걸쳐, 특히 시작 단계에서 워크플로우를 구축하고 유지 관리하는 동안 버그를 마주하게 될 것입니다. 체계적인 디버깅 접근법을 학습하면 문제를 빠르게 식별하고 해결하는 데 도움이 됩니다.

### 학습 목표

이 사이드 퀘스트에서는 Nextflow 워크플로우를 위한 **체계적인 디버깅 기법**을 학습합니다:

- **구문 오류 디버깅**: IDE 기능과 Nextflow 오류 메시지를 효과적으로 사용하기
- **채널 디버깅**: 데이터 흐름 문제와 채널 구조 문제 진단하기
- **프로세스 디버깅**: 실행 실패와 리소스 문제 조사하기
- **내장 디버깅 도구**: Nextflow의 미리보기 모드, 스텁 실행, 작업 디렉토리 활용하기
- **체계적 접근법**: 효율적인 디버깅을 위한 4단계 방법론

이 과정을 마치면 좌절스러운 오류 메시지를 명확한 해결 방법으로 전환하는 강력한 디버깅 방법론을 갖추게 됩니다.

### 사전 요구 사항

이 사이드 퀘스트를 시작하기 전에 다음을 충족해야 합니다:

- [Hello Nextflow](../hello_nextflow/README.md) 튜토리얼 또는 동등한 초급 과정을 완료했어야 합니다.
- 기본 Nextflow 개념과 메커니즘(프로세스, 채널, 연산자)을 편안하게 사용할 수 있어야 합니다.

**선택 사항:** [IDE Features for Nextflow Development](./ide_features.md) 사이드 퀘스트를 먼저 완료하는 것을 권장합니다.
해당 과정은 디버깅을 지원하는 IDE 기능(구문 강조, 오류 감지 등)을 포괄적으로 다루며, 여기서 이를 많이 사용할 것입니다.

---

## 0. 시작하기

#### 교육 코드스페이스 열기

아직 하지 않았다면 [환경 설정](../envsetup/index.md)에 설명된 대로 교육 환경을 열어야 합니다.

[![Open in GitHub Codespaces](https://github.com/codespaces/badge.svg)](https://codespaces.new/nextflow-io/training?quickstart=1&ref=master)

#### 프로젝트 디렉토리로 이동

이 튜토리얼의 파일이 있는 디렉토리로 이동하겠습니다.

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
2. **오류 메시지 분석**: Nextflow가 무엇을 알려주고 있나요?
3. **문제 위치 찾기**: 제공된 단서를 사용하여 코드에서 문제 찾기
4. **버그 수정**하고 솔루션이 작동하는지 확인하기
5. 다음 섹션으로 이동하기 전에 **파일 재설정**하기 (`git checkout <filename>` 사용)

연습은 간단한 구문 오류에서 더 미묘한 런타임 문제로 진행됩니다.
솔루션은 인라인으로 논의되지만, 앞으로 읽기 전에 각각을 직접 해결해 보세요.

#### 준비 체크리스트

시작할 준비가 되었다고 생각하시나요?

- [ ] 이 과정의 목표와 사전 요구 사항을 이해했습니다
- [ ] 코드스페이스가 실행 중입니다
- [ ] 작업 디렉토리를 적절하게 설정했습니다
- [ ] 과제를 이해했습니다

모든 항목을 체크할 수 있다면 시작할 준비가 된 것입니다.

---

## 1. 구문 오류

구문 오류는 Nextflow 코드를 작성할 때 마주치게 될 가장 일반적인 오류 유형입니다. 이는 코드가 Nextflow DSL의 예상 구문 규칙을 따르지 않을 때 발생합니다. 이러한 오류는 워크플로우가 전혀 실행되지 않도록 하므로, 빠르게 식별하고 수정하는 방법을 배우는 것이 중요합니다.

### 1.1. 중괄호 누락

가장 일반적인 구문 오류 중 하나이며, 때로는 디버깅하기 더 복잡한 오류 중 하나는 **중괄호 누락 또는 불일치**입니다.

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

- **파일과 위치**: 오류가 포함된 파일과 줄/열을 표시합니다 (`bad_syntax.nf:24:1`)
- **오류 설명**: 파서가 예상하지 못한 것을 발견했음을 설명합니다 (`Unexpected input: '<EOF>'`)
- **EOF 표시**: `<EOF>` (End Of File) 메시지는 파서가 더 많은 내용을 기대하면서 파일의 끝에 도달했음을 나타냅니다 - 닫히지 않은 중괄호의 전형적인 신호입니다

#### 코드 확인

이제 `bad_syntax.nf`를 검토하여 오류의 원인을 이해해 봅시다:

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
// 프로세스의 닫는 중괄호 누락

workflow {

    // 입력 채널 생성
    input_ch = channel.of('sample1', 'sample2', 'sample3')

    // 입력 채널로 프로세스 호출
    PROCESS_FILES(input_ch)
}
```

이 예제의 목적을 위해 오류가 어디에 있는지 보여주는 주석을 남겨두었습니다. Nextflow VSCode 확장 프로그램도 일치하지 않는 중괄호를 빨간색으로 표시하고 파일의 조기 종료를 강조하여 무엇이 잘못되었는지에 대한 힌트를 제공해야 합니다:

![Bad syntax](img/bad_syntax.png)

**중괄호 오류에 대한 디버깅 전략:**

1. VS Code의 중괄호 매칭 사용 (중괄호 옆에 커서 배치)
2. 중괄호 관련 메시지에 대한 Problems 패널 확인
3. 각 여는 `{`에 해당하는 닫는 `}`가 있는지 확인

#### 코드 수정

주석을 누락된 닫는 중괄호로 교체합니다:

=== "후"

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

        // 입력 채널로 프로세스 호출
        PROCESS_FILES(input_ch)
    }
    ```

=== "전"

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
    // 프로세스의 닫는 중괄호 누락

    workflow {

        // 입력 채널 생성
        input_ch = channel.of('sample1', 'sample2', 'sample3')

        // 입력 채널로 프로세스 호출
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

오류는 "Invalid process definition"을 나타내며 문제 주변의 컨텍스트를 보여줍니다. 3-7줄을 보면 4줄에 `inputs:`가 있는데, 이것이 문제입니다. `invalid_process.nf`를 검토해 봅시다:

```groovy title="invalid_process.nf" hl_lines="4" linenums="1"
#!/usr/bin/env nextflow

process PROCESS_FILES {
    inputs:  // 오류: 'inputs'가 아니라 'input'이어야 합니다
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

    // 입력 채널로 프로세스 호출
    PROCESS_FILES(input_ch)
}
```

오류 컨텍스트의 4줄을 보면 문제를 발견할 수 있습니다: 올바른 `input` 지시문 대신 `inputs`를 사용하고 있습니다. Nextflow VSCode 확장 프로그램도 이를 표시합니다:

![Invalid process message](img/invalid_process_message.png)

#### 코드 수정

[문서](https://www.nextflow.io/docs/latest/process.html#)를 참조하여 잘못된 키워드를 올바른 것으로 교체합니다:

=== "후"

    ```groovy title="invalid_process.nf" hl_lines="4" linenums="1"
    #!/usr/bin/env nextflow

    process PROCESS_FILES {
        input:  // 수정: 'inputs'를 'input'으로 변경
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

        // 입력 채널로 프로세스 호출
        PROCESS_FILES(input_ch)
    }
    ```

=== "전"

    ```groovy title="invalid_process.nf" hl_lines="4" linenums="1"
    #!/usr/bin/env nextflow

    process PROCESS_FILES {
        inputs:  // 오류: 'inputs'가 아니라 'input'이어야 합니다
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

        // 입력 채널로 프로세스 호출
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

스크립트 블록에서 사용하는 변수 이름은 유효해야 하며, 입력에서 파생되거나 스크립트 전에 삽입된 groovy 코드에서 파생되어야 합니다. 하지만 파이프라인 개발 초기에 복잡성을 다루다 보면 변수 이름 지정에서 실수하기 쉽고, Nextflow는 이를 빠르게 알려줄 것입니다.

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

오류는 컴파일 시점에 포착되며 17줄의 정의되지 않은 변수를 직접 가리키고, 캐럿이 정확히 문제가 있는 위치를 나타냅니다.

#### 코드 확인

`no_such_var.nf`를 검토해 봅시다:

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

오류 메시지는 변수가 스크립트 템플릿에서 인식되지 않음을 나타내며, 스크립트 블록에서 `${undefined_var}`가 사용되었지만 다른 곳에서 정의되지 않은 것을 볼 수 있습니다.

#### 코드 수정

'No such variable' 오류가 발생하면 변수를 정의하거나(입력 변수 이름을 수정하거나 스크립트 전에 groovy 코드를 편집하여) 필요하지 않은 경우 스크립트 블록에서 제거하여 수정할 수 있습니다:

=== "후"

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
        """  // undefined_var가 있는 줄 제거
    }

    workflow {
        input_ch = channel.of('sample1', 'sample2', 'sample3')
        PROCESS_FILES(input_ch)
    }
    ```

=== "전"

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

Nextflow를 시작할 때 Nextflow(Groovy) 변수와 Bash 변수의 차이를 이해하기 어려울 수 있습니다. 이는 스크립트 블록의 Bash 내용에서 변수를 사용하려고 할 때 나타나는 또 다른 형태의 잘못된 변수 오류를 생성할 수 있습니다.

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

오류는 `${prefix}`가 사용된 13줄을 가리킵니다. `bad_bash_var.nf`를 검토하여 무엇이 문제를 일으키는지 봅시다:

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
    echo "Processing ${sample_name}" > ${prefix}.txt  # 오류: ${prefix}는 Groovy 구문이지 Bash가 아닙니다
    """
}
```

이 예제에서는 Bash에서 `prefix` 변수를 정의하고 있지만, Nextflow 프로세스에서 이를 참조하는 데 사용한 `$` 구문(`${prefix}`)은 Bash가 아닌 Groovy 변수로 해석됩니다. 변수가 Groovy 컨텍스트에 존재하지 않으므로 'no such variable' 오류가 발생합니다.

#### 코드 수정

Bash 변수를 사용하려면 다음과 같이 달러 기호를 이스케이프해야 합니다:

=== "후"

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
        echo "Processing ${sample_name}" > \${prefix}.txt  # 수정: 달러 기호 이스케이프
        """
    }

    workflow {
        input_ch = channel.of('sample1', 'sample2', 'sample3')
        PROCESS_FILES(input_ch)
    }
    ```

=== "전"

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
        echo "Processing ${sample_name}" > ${prefix}.txt  # 오류: ${prefix}는 Groovy 구문이지 Bash가 아닙니다
        """
    }
    ```

이것은 Nextflow에게 이것을 Bash 변수로 해석하도록 알려줍니다.

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

!!! tip "Groovy 변수 vs Bash 변수"

    문자열 연결이나 접두사/접미사 작업과 같은 간단한 변수 조작의 경우, 스크립트 블록의 Bash 변수보다 스크립트 섹션의 Groovy 변수를 사용하는 것이 일반적으로 더 읽기 쉽습니다:

    ```groovy linenums="1"
    script:
    def output_prefix = "${sample_name}_processed"
    def output_file = "${output_prefix}.txt"
    """
    echo "Processing ${sample_name}" > ${output_file}
    """
    ```

    이 접근 방식은 달러 기호를 이스케이프할 필요가 없으며 코드를 더 읽기 쉽고 유지 관리하기 쉽게 만듭니다.

### 1.5. Workflow 블록 외부의 문장

Nextflow VSCode 확장 프로그램은 오류를 일으킬 코드 구조 문제를 강조합니다. 일반적인 예는 `workflow {}` 블록 외부에서 채널을 정의하는 것입니다 - 이것은 이제 구문 오류로 적용됩니다.

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

오류 메시지는 문제를 명확하게 나타냅니다: 문장(채널 정의와 같은)은 워크플로우나 프로세스 블록 외부의 스크립트 선언과 혼합될 수 없습니다.

#### 코드 확인

`badpractice_syntax.nf`를 검토하여 오류의 원인을 봅시다:

```groovy title="badpractice_syntax.nf" hl_lines="3" linenums="1"
#!/usr/bin/env nextflow

input_ch = channel.of('sample1', 'sample2', 'sample3')  // 오류: 워크플로우 외부에서 채널 정의

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

VSCode 확장 프로그램도 워크플로우 블록 외부에서 정의된 `input_ch` 변수를 강조합니다:

![Non-lethal syntax error](img/nonlethal.png)

#### 코드 수정

채널 정의를 워크플로우 블록 내부로 이동합니다:

=== "후"

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
        input_ch = channel.of('sample1', 'sample2', 'sample3')  // 워크플로우 블록 내부로 이동
        PROCESS_FILES(input_ch)
    }
    ```

=== "전"

    ```groovy title="badpractice_syntax.nf" hl_lines="3" linenums="1"
    #!/usr/bin/env nextflow

    input_ch = channel.of('sample1', 'sample2', 'sample3')  // 오류: 워크플로우 외부에서 채널 정의

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

입력 채널을 워크플로우 블록 내에 정의하고, 일반적으로 확장 프로그램이 제공하는 다른 권장 사항을 따르세요.

### 핵심 정리

Nextflow 오류 메시지와 IDE 시각적 표시기를 사용하여 구문 오류를 체계적으로 식별하고 수정할 수 있습니다. 일반적인 구문 오류에는 중괄호 누락, 잘못된 프로세스 키워드, 정의되지 않은 변수, Bash 변수와 Nextflow 변수의 부적절한 사용이 포함됩니다. VSCode 확장 프로그램은 런타임 전에 이러한 많은 오류를 포착하는 데 도움이 됩니다. 이러한 구문 디버깅 기술을 도구 상자에 넣으면 가장 일반적인 Nextflow 구문 오류를 빠르게 해결하고 더 복잡한 런타임 문제를 해결하는 데 집중할 수 있습니다.

### 다음 단계

구문이 올바른 경우에도 발생하는 더 복잡한 채널 구조 오류를 디버깅하는 방법을 학습합니다.

---

## 2. 채널 구조 오류

채널 구조 오류는 코드가 구문적으로 올바르지만 데이터 형태가 프로세스가 예상하는 것과 일치하지 않기 때문에 구문 오류보다 더 미묘합니다. Nextflow는 파이프라인을 실행하려고 시도하지만 입력 수가 예상과 일치하지 않는 것을 발견하고 실패할 수 있습니다. 이러한 오류는 일반적으로 런타임에만 나타나며 워크플로우를 통해 흐르는 데이터에 대한 이해가 필요합니다.

!!! tip "`.view()`로 채널 디버깅"

    이 섹션 전체에서 `.view()` 연산자를 사용하여 워크플로우의 어느 지점에서든 채널 내용을 검사할 수 있다는 것을 기억하세요. 이것은 채널 구조 문제를 이해하는 가장 강력한 디버깅 도구 중 하나입니다. 섹션 2.4에서 이 기법을 자세히 살펴보겠지만, 예제를 진행하면서 자유롭게 사용하세요.

    ```groovy
    my_channel.view()  // 채널을 통해 흐르는 것을 보여줍니다
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

오류 메시지는 호출이 1개의 인수를 예상했지만 2개를 받았다고 명확하게 설명하며 23줄을 가리킵니다. `bad_number_inputs.nf`를 검토해 봅시다:

```groovy title="bad_number_inputs.nf" hl_lines="5 23" linenums="1"
#!/usr/bin/env nextflow

process PROCESS_FILES {
    input:
        val sample_name  // 프로세스는 1개의 입력만 예상합니다

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

    // 오류: 2개의 채널을 전달하지만 프로세스는 1개만 예상합니다
    PROCESS_FILES(samples_ch, files_ch)
}
```

프로세스가 하나만 정의할 때 여러 입력 채널을 제공하는 일치하지 않는 `PROCESS_FILES` 호출을 볼 수 있어야 합니다. VSCode 확장 프로그램도 프로세스 호출에 빨간 밑줄을 그으며, 마우스를 올리면 진단 메시지를 제공합니다:

![Incorrect number of args message](img/incorrect_num_args.png)

#### 코드 수정

이 특정 예제의 경우 프로세스는 단일 채널을 예상하고 두 번째 채널이 필요하지 않으므로 `samples_ch` 채널만 전달하여 수정할 수 있습니다:

=== "후"

    ```groovy title="bad_number_inputs.nf" hl_lines="23" linenums="1"
    #!/usr/bin/env nextflow

    process PROCESS_FILES {
        input:
            val sample_name  // 프로세스는 1개의 입력만 예상합니다

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

        // 수정: 프로세스가 예상하는 채널만 전달
        PROCESS_FILES(samples_ch)
    }
    ```

=== "전"

    ```groovy title="bad_number_inputs.nf" hl_lines="5 23" linenums="1"
    #!/usr/bin/env nextflow

    process PROCESS_FILES {
        input:
            val sample_name  // 프로세스는 1개의 입력만 예상합니다

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

        // 오류: 2개의 채널을 전달하지만 프로세스는 1개만 예상합니다
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

이 예제보다 더 일반적으로, 프로세스에 추가 입력을 추가하고 그에 따라 워크플로우 호출을 업데이트하는 것을 잊어버려 이러한 유형의 오류가 발생할 수 있습니다. 다행히도 오류 메시지가 불일치에 대해 매우 명확하므로 이것은 이해하고 수정하기 쉬운 오류 중 하나입니다.

### 2.2. 채널 소진 (프로세스가 예상보다 적게 실행됨)

일부 채널 구조 오류는 훨씬 더 미묘하며 전혀 오류를 생성하지 않습니다. 아마도 이러한 오류 중 가장 일반적인 것은 큐 채널이 소진되어 항목이 부족해질 수 있다는 것을 이해하는 데 새로운 Nextflow 사용자가 직면하는 문제를 반영하며, 이는 워크플로우가 조기에 완료됨을 의미합니다.

#### 파이프라인 실행

```bash
nextflow run exhausted.nf
```

??? success "명령 출력"

```console title="소진된 채널 출력"
 N E X T F L O W   ~  version 25.10.2

Launching `exhausted.nf` [extravagant_gauss] DSL2 - revision: 08cff7ba2a

executor >  local (1)
[bd/f61fff] PROCESS_FILES (1) [100%] 1 of 1 ✔
```

이 워크플로우는 오류 없이 완료되지만 단일 샘플만 처리합니다!

#### 코드 확인

`exhausted.nf`를 검토하여 이것이 맞는지 봅시다:

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

프로세스가 세 번이 아닌 한 번만 실행되는 이유는 `reference_ch` 채널이 첫 번째 프로세스 실행 후 소진되는 큐 채널이기 때문입니다. 한 채널이 소진되면 다른 채널에 여전히 항목이 있더라도 전체 프로세스가 중지됩니다.

이것은 여러 샘플에서 재사용해야 하는 단일 참조 파일이 있는 일반적인 패턴입니다. 해결책은 참조 채널을 무한정 재사용할 수 있는 값 채널로 변환하는 것입니다.

#### 코드 수정

영향을 받는 파일 수에 따라 이를 해결하는 몇 가지 방법이 있습니다.

**옵션 1**: 많이 재사용하는 단일 참조 파일이 있습니다. 반복해서 사용할 수 있는 값 채널 유형을 간단히 생성할 수 있습니다. 이를 수행하는 세 가지 방법이 있습니다:

**1a** `channel.value()` 사용:

```groovy title="exhausted.nf (수정 - 옵션 1a)" hl_lines="2" linenums="21"
workflow {
    reference_ch = channel.value('baseline_reference')  // 값 채널은 재사용 가능
    input_ch = channel.of('sample1', 'sample2', 'sample3')

    PROCESS_FILES(reference_ch, input_ch)
}
```

**1b** `first()` [연산자](https://www.nextflow.io/docs/latest/reference/operator.html#first) 사용:

```groovy title="exhausted.nf (수정 - 옵션 1b)" hl_lines="2" linenums="21"
workflow {
    reference_ch = channel.of('baseline_reference').first()  // 값 채널로 변환
    input_ch = channel.of('sample1', 'sample2', 'sample3')

    PROCESS_FILES(reference_ch, input_ch)
}
```

**1c.** `collect()` [연산자](https://www.nextflow.io/docs/latest/reference/operator.html#collect) 사용:

```groovy title="exhausted.nf (수정 - 옵션 1c)" hl_lines="2" linenums="21"
workflow {
    reference_ch = channel.of('baseline_reference').collect()  // 값 채널로 변환
    input_ch = channel.of('sample1', 'sample2', 'sample3')

    PROCESS_FILES(reference_ch, input_ch)
}
```

**옵션 2**: 더 복잡한 시나리오에서, 아마도 샘플 채널의 모든 샘플에 대한 여러 참조 파일이 있는 경우, `combine` 연산자를 사용하여 두 채널을 튜플로 결합하는 새 채널을 생성할 수 있습니다:

```groovy title="exhausted.nf (수정 - 옵션 2)" hl_lines="4" linenums="21"
workflow {
    reference_ch = channel.of('baseline_reference','other_reference')
    input_ch = channel.of('sample1', 'sample2', 'sample3')
    combined_ch = reference_ch.combine(input_ch)  // 데카르트 곱 생성

    PROCESS_FILES(combined_ch)
}
```

`.combine()` 연산자는 두 채널의 데카르트 곱을 생성하므로 `reference_ch`의 각 항목이 `input_ch`의 각 항목과 쌍을 이룹니다. 이를 통해 프로세스가 참조를 사용하면서 각 샘플에 대해 실행될 수 있습니다.

이를 위해서는 프로세스 입력을 조정해야 합니다. 우리 예제에서 프로세스 정의의 시작 부분을 다음과 같이 조정해야 합니다:

```groovy title="exhausted.nf (수정 - 옵션 2)" hl_lines="5" linenums="1"
#!/usr/bin/env nextflow

process PROCESS_FILES {
    input:
        tuple val(reference), val(sample_name)
```

이 접근 방식은 모든 상황에 적합하지 않을 수 있습니다.

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

이제 하나가 아닌 세 개의 샘플이 모두 처리되는 것을 볼 수 있어야 합니다.

### 2.3. 잘못된 채널 내용 구조

워크플로우가 특정 수준의 복잡성에 도달하면 각 채널의 내부 구조를 추적하기가 조금 어려울 수 있으며, 사람들은 일반적으로 프로세스가 예상하는 것과 채널이 실제로 포함하는 것 사이에 불일치를 생성합니다. 이것은 채널 수가 잘못된 앞서 논의한 문제보다 더 미묘합니다. 이 경우 올바른 수의 입력 채널을 가질 수 있지만 하나 이상의 채널의 내부 구조가 프로세스가 예상하는 것과 일치하지 않습니다.

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

오류 메시지의 대괄호가 여기서 단서를 제공합니다 - 프로세스가 튜플을 단일 값으로 처리하고 있으며, 이것은 우리가 원하는 것이 아닙니다. `bad_channel_shape.nf`를 검토해 봅시다:

```groovy title="bad_channel_shape.nf" hl_lines="5 20-22" linenums="1"
#!/usr/bin/env nextflow

process PROCESS_FILES {
    input:
        val sample_name  // 단일 값을 예상하지만 튜플을 받습니다

    output:
        path "${sample_name}_output.txt"

    script:
    """
    echo "Processing ${sample_name}" > ${sample_name}_output.txt
    """
}

workflow {

    // 채널이 튜플을 내보내지만 프로세스는 단일 값을 예상합니다
    input_ch = channel.of(
      ['sample1', 'file1.txt'],
      ['sample2', 'file2.txt'],
      ['sample3', 'file3.txt']
    )
    PROCESS_FILES(input_ch)
}
```

튜플로 구성된 채널을 생성하고 있음을 볼 수 있습니다: `['sample1', 'file1.txt']`, 하지만 프로세스는 단일 값 `val sample_name`을 예상합니다. 실행된 명령은 프로세스가 `[sample3, file3.txt]_output.txt`라는 파일을 생성하려고 시도하고 있음을 보여주며, 이것은 의도된 출력이 아닙니다.

#### 코드 수정

이를 수정하려면 프로세스가 두 입력을 모두 필요로 하는 경우 튜플을 받도록 프로세스를 조정할 수 있습니다:

=== "옵션 1: 프로세스에서 튜플 받기"

    === "후"

        ```groovy title="bad_channel_shape.nf" hl_lines="5"  linenums="1"
        #!/usr/bin/env nextflow

        process PROCESS_FILES {
            input:
                tuple val(sample_name), val(file_name)  // 수정: 튜플 받기

            output:
                path "${sample_name}_output.txt"

            script:
            """
            echo "Processing ${sample_name}" > ${sample_name}_output.txt
            """
        }

        workflow {

            // 채널이 튜플을 내보내지만 프로세스는 단일 값을 예상합니다
            input_ch = channel.of(
              ['sample1', 'file1.txt'],
              ['sample2', 'file2.txt'],
              ['sample3', 'file3.txt']
            )
            PROCESS_FILES(input_ch)
        }
        ```

    === "전"

        ```groovy title="bad_channel_shape.nf" hl_lines="5" linenums="1"
        #!/usr/bin/env nextflow

        process PROCESS_FILES {
            input:
                val sample_name  // 단일 값을 예상하지만 튜플을 받습니다

            output:
                path "${sample_name}_output.txt"

            script:
            """
            echo "Processing ${sample_name}" > ${sample_name}_output.txt
            """
        }

        workflow {

            // 채널이 튜플을 내보내지만 프로세스는 단일 값을 예상합니다
            input_ch = channel.of(
              ['sample1', 'file1.txt'],
              ['sample2', 'file2.txt'],
              ['sample3', 'file3.txt']
            )
            PROCESS_FILES(input_ch)
        }
        ```

=== "옵션 2: 첫 번째 요소 추출"

    === "후"

        ```groovy title="bad_channel_shape.nf" hl_lines="9" linenums="16"
        workflow {

            // 채널이 튜플을 내보내지만 프로세스는 단일 값을 예상합니다
            input_ch = channel.of(
              ['sample1', 'file1.txt'],
              ['sample2', 'file2.txt'],
              ['sample3', 'file3.txt']
            )
            PROCESS_FILES(input_ch.map { it[0] })  // 수정: 첫 번째 요소 추출
        }
        ```

    === "전"

        ```groovy title="bad_channel_shape.nf" hl_lines="9" linenums="16"
        workflow {

            // 채널이 튜플을 내보내지만 프로세스는 단일 값을 예상합니다
            input_ch = channel.of(
              ['sample1', 'file1.txt'],
              ['sample2', 'file2.txt'],
              ['sample3', 'file3.txt']
            )
            PROCESS_FILES(input_ch)
        }
        ```

#### 파이프라인 실행

솔루션 중 하나를 선택하고 워크플로우를 다시 실행합니다:

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

`bad_channel_shape_viewed.nf`를 실행하여 이것이 작동하는 것을 확인합니다:

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

`bad_channel_shape_viewed.nf`를 검토하여 `.view()`가 어떻게 사용되는지 봅시다:

```groovy title="bad_channel_shape_viewed.nf" linenums="16" hl_lines="9 11"
workflow {

    // 채널이 튜플을 내보내지만 프로세스는 단일 값을 예상합니다
    input_ch = channel.of(
      ['sample1', 'file1.txt'],
      ['sample2', 'file2.txt'],
      ['sample3', 'file3.txt']
    )
    .view { "Channel content: $it" }  // 디버그: 원본 채널 내용 표시
    .map { tuple -> tuple[0] }        // 변환: 첫 번째 요소 추출
    .view { "After mapping: $it" }    // 디버그: 변환된 채널 내용 표시

    PROCESS_FILES(input_ch)
}
```

#### 코드 수정

향후 채널 내용을 이해하기 위해 `.view()` 작업을 과도하게 사용하지 않도록 하려면 도움이 되는 주석을 추가하는 것이 좋습니다:

```groovy title="bad_channel_shape_viewed.nf (주석 포함)" linenums="16" hl_lines="8 9"
workflow {

    // 채널이 튜플을 내보내지만 프로세스는 단일 값을 예상합니다
    input_ch = channel.of(
            ['sample1', 'file1.txt'],
            ['sample2', 'file2.txt'],
            ['sample3', 'file3.txt'],
        ) // [sample_name, file_name]
        .map { tuple -> tuple[0] } // sample_name

    PROCESS_FILES(input_ch)
}
```

이것은 워크플로우가 복잡성이 증가하고 채널 구조가 더 불투명해짐에 따라 더 중요해질 것입니다.

#### 파이프라인 실행

```bash
nextflow run bad_channel_shape_viewed.nf
```

??? success "명령 출력"

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

### 핵심 정리

많은 채널 구조 오류는 유효한 Nextflow 구문으로 생성될 수 있습니다. 데이터 흐름을 이해하고, 검사를 위해 `.view()` 연산자를 사용하고, 예상치 못한 튜플 구조를 나타내는 대괄호와 같은 오류 메시지 패턴을 인식하여 채널 구조 오류를 디버깅할 수 있습니다.

### 다음 단계

프로세스 정의로 인해 생성된 오류에 대해 학습합니다.

---

## 3. 프로세스 구조 오류

프로세스와 관련하여 마주치게 될 대부분의 오류는 명령을 형성하는 데 실수를 하거나 기본 소프트웨어와 관련된 문제와 관련이 있습니다. 그렇긴 하지만, 위의 채널 문제와 유사하게 구문 오류로 인정되지 않지만 런타임에 오류를 일으킬 프로세스 정의에서 실수를 할 수 있습니다.

### 3.1. 출력 파일 누락

프로세스를 작성할 때 일반적인 오류 중 하나는 프로세스가 예상하는 것과 생성되는 것 사이에 불일치를 생성하는 작업을 수행하는 것입니다.

#### 파이프라인 실행

```bash
nextflow run missing_output.nf
```

??? failure "명령 출력"

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

오류 메시지는 프로세스가 `sample3.txt`라는 출력 파일을 생성할 것으로 예상했지만 스크립트가 실제로 `sample3_output.txt`를 생성한다고 나타냅니다. `missing_output.nf`의 프로세스 정의를 검토해 봅시다:

```groovy title="missing_output.nf" linenums="3" hl_lines="6 10"
process PROCESS_FILES {
    input:
    val sample_name

    output:
    path "${sample_name}.txt"  // 예상: sample3.txt

    script:
    """
    echo "Processing ${sample_name}" > ${sample_name}_output.txt  // 생성: sample3_output.txt
    """
}
```

`output:` 블록의 출력 파일 이름과 스크립트에서 사용된 파일 이름 사이에 불일치가 있음을 볼 수 있어야 합니다. 이 불일치로 인해 프로세스가 실패합니다. 이러한 종류의 오류가 발생하면 돌아가서 프로세스 정의와 출력 블록 사이에 출력이 일치하는지 확인하세요.

문제가 여전히 명확하지 않으면 작업 디렉토리 자체를 확인하여 생성된 실제 출력 파일을 식별합니다:

```bash
❯ ls -h work/02/9604d49fb8200a74d737c72a6c98ed
sample3_output.txt
```

이 예제의 경우 `output:` 정의와 달리 `_output` 접미사가 출력 파일 이름에 통합되고 있음을 강조합니다.

#### 코드 수정

출력 파일 이름을 일관되게 만들어 불일치를 수정합니다:

=== "후"

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

=== "전"

    ```groovy title="missing_output.nf" hl_lines="6 10" linenums="3"
    process PROCESS_FILES {
        input:
        val sample_name

        output:
        path "${sample_name}.txt"  // 예상: sample3.txt

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

??? success "명령 출력"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `missing_output.nf` [elated_hamilton] DSL2 - revision: 961938ee2b

    executor >  local (3)
    [16/1c437c] PROCESS_FILES (3) | 3 of 3 ✔
    ```

### 3.2. 소프트웨어 누락

또 다른 오류 클래스는 소프트웨어 프로비저닝의 실수로 인해 발생합니다. `missing_software.nf`는 구문적으로 유효한 워크플로우이지만 사용하는 `cowpy` 명령을 제공하기 위해 일부 외부 소프트웨어에 의존합니다.

#### 파이프라인 실행

```bash
nextflow run missing_software.nf
```

??? failure "명령 출력"

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

프로세스가 지정하는 명령에 액세스할 수 없습니다. 때로는 스크립트가 워크플로우 `bin` 디렉토리에 있지만 실행 가능하게 만들어지지 않았기 때문입니다. 다른 경우에는 워크플로우가 실행되는 컨테이너나 환경에 소프트웨어가 설치되지 않았기 때문입니다.

#### 코드 확인

`127` 종료 코드를 주의하세요 - 정확히 문제를 알려줍니다. `missing_software.nf`를 검토해 봅시다:

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

여기서 우리는 조금 부정직했으며, 실제로 코드에는 아무런 문제가 없습니다. 문제의 명령에 액세스할 수 있는 방식으로 프로세스를 실행하는 데 필요한 구성을 지정하기만 하면 됩니다. 이 경우 프로세스에 컨테이너 정의가 있으므로 Docker를 활성화하여 워크플로우를 실행하기만 하면 됩니다.

#### 파이프라인 실행

`nextflow.config`에 Docker 프로파일을 설정해 두었으므로 다음과 같이 워크플로우를 실행할 수 있습니다:

```bash
nextflow run missing_software.nf -profile docker
```

??? success "명령 출력"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `missing_software.nf` [awesome_stonebraker] DSL2 - revision: 0296d12839

    executor >  local (3)
    [38/ab20d1] PROCESS_FILES (1) | 3 of 3 ✔
    ```

!!! note "참고"

    Nextflow가 컨테이너를 사용하는 방법에 대해 자세히 알아보려면 [Hello Nextflow](../hello_nextflow/05_hello_containers.md)를 참조하세요

### 3.3. 잘못된 리소스 구성

프로덕션 사용에서는 프로세스에 리소스를 구성하게 됩니다. 예를 들어 `memory`는 프로세스에 사용 가능한 최대 메모리 양을 정의하며, 프로세스가 이를 초과하면 스케줄러는 일반적으로 프로세스를 종료하고 종료 코드 `137`을 반환합니다. `local` 실행자를 사용하고 있기 때문에 여기서는 이를 시연할 수 없지만 `time`과 유사한 것을 보여줄 수 있습니다.

#### 파이프라인 실행

`bad_resources.nf`는 1밀리초라는 비현실적인 시간 제한이 있는 프로세스 구성을 가지고 있습니다:

```bash
nextflow run bad_resources.nf -profile docker
```

??? failure "명령 출력"

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

`bad_resources.nf`를 검토해 봅시다:

```groovy title="bad_resources.nf" linenums="3" hl_lines="3"
process PROCESS_FILES {

    time '1 ms'  // 오류: 비현실적인 시간 제한

    input:
    val sample_name

    output:
    path "${sample_name}_output.txt"

    script:
    """
    sleep 1  // 1초가 걸리지만 시간 제한은 1ms입니다
    cowpy ${sample_name} > ${sample_name}_output.txt
    """
}
```

프로세스가 1초 이상 걸릴 것을 알고 있지만(확실히 하기 위해 sleep을 추가했습니다), 프로세스는 1밀리초 후에 시간 초과되도록 설정되어 있습니다. 누군가 구성에 대해 조금 비현실적이었습니다!

#### 코드 수정

시간 제한을 현실적인 값으로 늘립니다:

=== "후"

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

=== "전"

    ```groovy title="bad_resources.nf" hl_lines="3" linenums="3"
    process PROCESS_FILES {

        time '1 ms'  // 오류: 비현실적인 시간 제한

        input:
        val sample_name

        output:
        path "${sample_name}_output.txt"

        script:
        """
        sleep 1  // 1초가 걸리지만 시간 제한은 1ms입니다
        cowpy ${sample_name} > ${sample_name}_output.txt
        """
    }
    ```

#### 파이프라인 실행

```bash
nextflow run bad_resources.nf -profile docker
```

??? success "명령 출력"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `bad_resources.nf` [friendly_mcclintock] DSL2 - revision: 381567d2c1

    executor >  local (3)
    [c2/9b4c41] PROCESS_FILES (3) | 3 of 3 ✔
    ```

오류 메시지를 주의 깊게 읽으면 이와 같은 실패가 오랫동안 당황스럽지 않을 것입니다. 하지만 리소스 지시문을 적절하게 구성할 수 있도록 실행 중인 명령의 리소스 요구 사항을 이해해야 합니다.

### 3.4. 프로세스 디버깅 기법

프로세스가 실패하거나 예상치 못하게 동작할 때 무엇이 잘못되었는지 조사하기 위한 체계적인 기법이 필요합니다. 작업 디렉토리에는 프로세스 실행을 디버깅하는 데 필요한 모든 정보가 포함되어 있습니다.

#### 작업 디렉토리 검사 사용

프로세스를 위한 가장 강력한 디버깅 도구는 작업 디렉토리를 검사하는 것입니다. 프로세스가 실패하면 Nextflow는 무슨 일이 일어났는지 이해하는 데 필요한 모든 파일을 포함하는 특정 프로세스 실행을 위한 작업 디렉토리를 생성합니다.

#### 파이프라인 실행

앞서의 `missing_output.nf` 예제를 사용하여 작업 디렉토리 검사를 시연해 봅시다(필요한 경우 출력 이름 불일치를 다시 생성하세요):

```bash
nextflow run missing_output.nf
```

??? failure "명령 출력"

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

     -- Check '.nextflow.log' file for details
    ```

#### 작업 디렉토리 확인

이 오류가 발생하면 작업 디렉토리에 모든 디버깅 정보가 포함됩니다. 오류 메시지에서 작업 디렉토리 경로를 찾아 내용을 검사합니다:

```bash
# 오류 메시지에서 작업 디렉토리 찾기
ls work/02/9604d49fb8200a74d737c72a6c98ed/
```

그런 다음 주요 파일을 검사할 수 있습니다:

##### 명령 스크립트 확인

`.command.sh` 파일은 정확히 어떤 명령이 실행되었는지 보여줍니다:

```bash
# 실행된 명령 보기
cat work/02/9604d49fb8200a74d737c72a6c98ed/.command.sh
```

이것은 다음을 보여줍니다:

- **변수 치환**: Nextflow 변수가 제대로 확장되었는지 여부
- **파일 경로**: 입력 파일이 올바르게 위치했는지 여부
- **명령 구조**: 스크립트 구문이 올바른지 여부

찾아야 할 일반적인 문제:

- **따옴표 누락**: 공백을 포함하는 변수는 적절한 따옴표가 필요합니다
- **잘못된 파일 경로**: 존재하지 않거나 잘못된 위치에 있는 입력 파일
- **잘못된 변수 이름**: 변수 참조의 오타
- **환경 설정 누락**: 특정 환경에 의존하는 명령

##### 오류 출력 확인

`.command.err` 파일에는 실제 오류 메시지가 포함되어 있습니다:

```bash
# 오류 출력 보기
cat work/02/9604d49fb8200a74d737c72a6c98ed/.command.err
```

이 파일은 다음을 보여줍니다:

- **종료 코드**: 127 (명령을 찾을 수 없음), 137 (종료됨) 등
- **권한 오류**: 파일 액세스 문제
- **소프트웨어 오류**: 애플리케이션별 오류 메시지
- **리소스 오류**: 메모리/시간 제한 초과

##### 표준 출력 확인

`.command.out` 파일은 명령이 생성한 것을 보여줍니다:

```bash
# 표준 출력 보기
cat work/02/9604d49fb8200a74d737c72a6c98ed/.command.out
```

이것은 다음을 확인하는 데 도움이 됩니다:

- **예상 출력**: 명령이 올바른 결과를 생성했는지 여부
- **부분 실행**: 명령이 시작되었지만 중간에 실패했는지 여부
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

출력 파일 누락으로 인해 프로세스가 실패하면 실제로 생성된 파일을 확인합니다:

```bash
# 작업 디렉토리의 모든 파일 나열
ls -la work/02/9604d49fb8200a74d737c72a6c98ed/
```

이것은 다음을 식별하는 데 도움이 됩니다:

- **파일 이름 불일치**: 예상과 다른 이름의 출력 파일
- **권한 문제**: 생성할 수 없는 파일
- **경로 문제**: 잘못된 디렉토리에 생성된 파일

앞서 예제에서 이것은 예상한 `sample3.txt`가 없지만 `sample3_output.txt`가 있음을 확인했습니다:

```bash
❯ ls -h work/02/9604d49fb8200a74d737c72a6c98ed
sample3_output.txt
```

### 핵심 정리

프로세스 디버깅은 무엇이 잘못되었는지 이해하기 위해 작업 디렉토리를 검사해야 합니다. 주요 파일에는 `.command.sh`(실행된 스크립트), `.command.err`(오류 메시지), `.command.out`(표준 출력)이 포함됩니다. 127(명령을 찾을 수 없음) 및 137(프로세스 종료됨)과 같은 종료 코드는 실패 유형에 대한 즉각적인 진단 단서를 제공합니다.

### 다음 단계

Nextflow의 내장 디버깅 도구와 문제 해결을 위한 체계적 접근법에 대해 학습합니다.

---

## 4. 내장 디버깅 도구 및 고급 기법

Nextflow는 워크플로우 실행을 디버깅하고 분석하기 위한 여러 강력한 내장 도구를 제공합니다. 이러한 도구는 무엇이 잘못되었는지, 어디서 잘못되었는지, 효율적으로 수정하는 방법을 이해하는 데 도움이 됩니다.

### 4.1. 실시간 프로세스 출력

때로는 실행 중인 프로세스 내부에서 무슨 일이 일어나고 있는지 확인해야 합니다. 실시간 프로세스 출력을 활성화할 수 있으며, 이는 각 작업이 실행될 때 정확히 무엇을 하고 있는지 보여줍니다.

#### 파이프라인 실행

앞서 예제의 `bad_channel_shape_viewed.nf`는 `.view()`를 사용하여 채널 내용을 출력했지만, `debug` 지시문을 사용하여 프로세스 자체 내에서 변수를 에코할 수도 있으며, 이를 `bad_channel_shape_viewed_debug.nf`에서 시연합니다. 워크플로우를 실행합니다:

```bash
nextflow run bad_channel_shape_viewed_debug.nf
```

??? success "명령 출력"

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

`bad_channel_shape_viewed_debug.nf`를 검토하여 `debug` 지시문이 어떻게 작동하는지 봅시다:

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

`debug` 지시문은 프로세스의 환경을 이해하는 빠르고 편리한 방법이 될 수 있습니다.

### 4.2. 미리보기 모드

때로는 프로세스가 실행되기 전에 문제를 포착하고 싶습니다. Nextflow는 이러한 종류의 사전 디버깅을 위한 플래그를 제공합니다: `-preview`.

#### 파이프라인 실행

미리보기 모드를 사용하면 명령을 실행하지 않고 워크플로우 로직을 테스트할 수 있습니다. 이것은 실제 명령을 실행하지 않고 워크플로우의 구조를 빠르게 확인하고 프로세스가 올바르게 연결되었는지 확인하는 데 매우 유용할 수 있습니다.

!!! note "참고"

    앞서 `bad_syntax.nf`를 수정했다면 이 명령을 실행하기 전에 스크립트 블록 뒤의 닫는 중괄호를 제거하여 구문 오류를 다시 도입하세요.

이 명령을 실행합니다:

```bash
nextflow run bad_syntax.nf -preview
```

??? failure "명령 출력"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `bad_syntax.nf` [magical_mercator] DSL2 - revision: 550b9a8873

    Error bad_syntax.nf:24:1: Unexpected input: '<EOF>'

    ERROR ~ Script compilation failed

     -- Check '.nextflow.log' file for details
    ```

미리보기 모드는 프로세스를 실행하지 않고 구문 오류를 조기에 포착하는 데 특히 유용합니다. 실행 전에 워크플로우 구조와 프로세스 연결을 검증합니다.

### 4.3. 로직 테스트를 위한 스텁 실행

때로는 명령이 너무 오래 걸리거나, 특수 소프트웨어가 필요하거나, 복잡한 이유로 실패하기 때문에 오류를 디버깅하기 어렵습니다. 스텁 실행을 사용하면 실제 명령을 실행하지 않고 워크플로우 로직을 테스트할 수 있습니다.

#### 파이프라인 실행

Nextflow 프로세스를 개발할 때 `stub` 지시문을 사용하여 실제 명령을 실행하지 않고 올바른 형태의 출력을 생성하는 '더미' 명령을 정의할 수 있습니다. 이 접근 방식은 실제 소프트웨어의 복잡성을 다루기 전에 워크플로우 로직이 올바른지 확인하려는 경우 특히 유용합니다.

예를 들어, 앞서의 `missing_software.nf`를 기억하시나요? `-profile docker`를 추가할 때까지 워크플로우가 실행되지 않도록 하는 소프트웨어가 누락된 것? `missing_software_with_stub.nf`는 매우 유사한 워크플로우입니다. 같은 방식으로 실행하면 같은 오류가 생성됩니다:

```bash
nextflow run missing_software_with_stub.nf
```

??? failure "명령 출력"

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

그러나 이 워크플로우는 `docker` 프로파일 없이도 `-stub-run`으로 실행하면 오류를 생성하지 않습니다:

```bash
nextflow run missing_software_with_stub.nf -stub-run
```

??? success "명령 출력"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `missing_software_with_stub.nf` [astonishing_shockley] DSL2 - revision: f1f4f05d7d

    executor >  local (3)
    [b5/2517a3] PROCESS_FILES (3) | 3 of 3 ✔
    ```

#### 코드 확인

`missing_software_with_stub.nf`를 검토해 봅시다:

```groovy title="missing_software.nf (스텁 포함)" hl_lines="16-19" linenums="3"
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

`missing_software.nf`에 비해 이 프로세스에는 Nextflow가 스텁 모드에서 실행되는 경우 `script:`에 지정된 명령 대신 사용할 명령을 지정하는 `stub:` 지시문이 있습니다.

여기서 사용하는 `touch` 명령은 소프트웨어나 적절한 입력에 의존하지 않으며 모든 상황에서 실행되므로 프로세스 내부에 대해 걱정하지 않고 워크플로우 로직을 디버깅할 수 있습니다.

**스텁 실행이 디버깅하는 데 도움이 되는 것:**

- 채널 구조 및 데이터 흐름
- 프로세스 연결 및 종속성
- 매개변수 전파
- 소프트웨어 종속성 없는 워크플로우 로직

### 4.4. 체계적 디버깅 접근법

이제 추적 파일과 작업 디렉토리에서 미리보기 모드, 스텁 실행, 리소스 모니터링에 이르기까지 개별 디버깅 기법을 배웠으므로 이를 체계적인 방법론으로 묶어봅시다. 구조화된 접근 방식을 갖추면 복잡한 오류에 압도되지 않고 중요한 단서를 놓치지 않을 수 있습니다.

이 방법론은 우리가 다룬 모든 도구를 효율적인 워크플로우로 결합합니다:

**4단계 디버깅 방법:**

**1단계: 구문 오류 해결 (5분)**

1. VSCode 또는 IDE에서 빨간 밑줄 확인
2. `nextflow run workflow.nf -preview`를 실행하여 구문 문제 식별
3. 모든 구문 오류 수정 (중괄호 누락, 후행 쉼표 등)
4. 진행하기 전에 워크플로우가 성공적으로 파싱되는지 확인

**2단계: 빠른 평가 (5분)**

1. 런타임 오류 메시지를 주의 깊게 읽기
2. 런타임, 로직 또는 리소스 오류인지 확인
3. 미리보기 모드를 사용하여 기본 워크플로우 로직 테스트

**3단계: 상세 조사 (15-30분)**

1. 실패한 작업의 작업 디렉토리 찾기
2. 로그 파일 검사
3. 채널을 검사하기 위해 `.view()` 연산자 추가
4. `-stub-run`을 사용하여 실행 없이 워크플로우 로직 테스트

**4단계: 수정 및 검증 (15분)**

1. 최소한의 대상 수정 수행
2. resume으로 테스트: `nextflow run workflow.nf -resume`
3. 완전한 워크플로우 실행 확인

!!! tip "효율적인 디버깅을 위한 Resume 사용"

    문제를 식별한 후에는 워크플로우의 성공적인 부분을 다시 실행하는 데 시간을 낭비하지 않고 수정 사항을 테스트하는 효율적인 방법이 필요합니다. Nextflow의 `-resume` 기능은 디버깅에 매우 유용합니다.

    [Hello Nextflow](../hello_nextflow/)를 진행했다면 `-resume`을 접했을 것이며, 문제 프로세스 전의 프로세스가 실행되는 동안 기다리는 시간을 절약하기 위해 디버깅할 때 이를 잘 활용하는 것이 중요합니다.

    **Resume 디버깅 전략:**

    1. 실패할 때까지 워크플로우 실행
    2. 실패한 작업의 작업 디렉토리 검사
    3. 특정 문제 수정
    4. Resume하여 수정 사항만 테스트
    5. 워크플로우가 완료될 때까지 반복

#### 디버깅 구성 프로파일

이 체계적인 접근 방식을 더욱 효율적으로 만들기 위해 필요한 모든 도구를 자동으로 활성화하는 전용 디버깅 구성을 생성할 수 있습니다:

```groovy title="nextflow.config (디버그 프로파일)" linenums="1"
profiles {
    debug {
        process {
            debug = true
            cleanup = false

            // 디버깅을 위한 보수적인 리소스
            maxForks = 1
            memory = '2.GB'
            cpus = 1
        }
    }
}
```

그런 다음 이 프로파일을 활성화하여 파이프라인을 실행할 수 있습니다:

```bash
nextflow run workflow.nf -profile debug
```

이 프로파일은 실시간 출력을 활성화하고, 작업 디렉토리를 보존하며, 더 쉬운 디버깅을 위해 병렬화를 제한합니다.

### 4.5. 실용적인 디버깅 연습

이제 체계적인 디버깅 접근법을 실제로 적용할 시간입니다. 워크플로우 `buggy_workflow.nf`에는 실제 개발에서 마주치게 될 유형의 문제를 나타내는 여러 일반적인 오류가 포함되어 있습니다.

!!! exercise "연습"

    체계적인 디버깅 접근법을 사용하여 `buggy_workflow.nf`의 모든 오류를 식별하고 수정하세요. 이 워크플로우는 CSV 파일에서 샘플 데이터를 처리하려고 시도하지만 일반적인 디버깅 시나리오를 나타내는 여러 의도적인 버그를 포함하고 있습니다.

    워크플로우를 실행하여 첫 번째 오류를 확인하는 것으로 시작하세요:

    ```bash
    nextflow run buggy_workflow.nf
    ```

    ??? failure "명령 출력"

        ```console
        N E X T F L O W   ~  version 25.10.2

        Launching `buggy_workflow.nf` [wise_ramanujan] DSL2 - revision: d51a8e83fd

        ERROR ~ Range [11, 12) out of bounds for length 11

         -- Check '.nextflow.log' file for details
        ```

        이 암호 같은 오류는 `params{}` 블록의 11-12줄 주변의 파싱 문제를 나타냅니다. v2 파서는 구조적 문제를 조기에 포착합니다.

    학습한 4단계 디버깅 방법을 적용하세요:

    **1단계: 구문 오류 해결**
    - VSCode 또는 IDE에서 빨간 밑줄 확인
    - `nextflow run workflow.nf -preview`를 실행하여 구문 문제 식별
    - 모든 구문 오류 수정 (중괄호 누락, 후행 쉼표 등)
    - 진행하기 전에 워크플로우가 성공적으로 파싱되는지 확인

    **2단계: 빠른 평가**
    - 런타임 오류 메시지를 주의 깊게 읽기
    - 오류가 런타임, 로직 또는 리소스 관련인지 식별
    - `-preview` 모드를 사용하여 기본 워크플로우 로직 테스트

    **3단계: 상세 조사**
    - 실패한 작업의 작업 디렉토리 검사
    - 채널을 검사하기 위해 `.view()` 연산자 추가
    - 작업 디렉토리의 로그 파일 확인
    - `-stub-run`을 사용하여 실행 없이 워크플로우 로직 테스트

    **4단계: 수정 및 검증**
    - 대상 수정 수행
    - `-resume`을 사용하여 수정 사항을 효율적으로 테스트
    - 완전한 워크플로우 실행 확인

    **사용 가능한 디버깅 도구:**
    ```bash
    # 구문 확인을 위한 미리보기 모드
    nextflow run buggy_workflow.nf -preview

    # 상세 출력을 위한 디버그 프로파일
    nextflow run buggy_workflow.nf -profile debug

    # 로직 테스트를 위한 스텁 실행
    nextflow run buggy_workflow.nf -stub-run

    # 수정 후 Resume
    nextflow run buggy_workflow.nf -resume
    ```

    ??? solution "해결책"
        `buggy_workflow.nf`에는 모든 주요 디버깅 범주를 다루는 9개 또는 10개의 개별 오류가 포함되어 있습니다(세는 방법에 따라). 다음은 각 오류와 수정 방법에 대한 체계적인 분석입니다

        구문 오류부터 시작하겠습니다:

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

        **오류 2: 구문 오류 - 닫는 중괄호 누락**
        ```groovy linenums="24"
        script:
        """
        echo "Processing: ${sample}"
        cat ${input_file} > ${sample}_result.txt
        """
        // 오류: processFiles 프로세스의 닫는 중괄호 누락
        ```
        **수정:** 누락된 닫는 중괄호 추가
        ```groovy linenums="29"
        """
        echo "Processing: ${sample_id}"
        cat ${input_file} > ${sample_id}_result.txt
        """
        }  // 누락된 닫는 중괄호 추가
        ```

        **오류 3: 변수 이름 오류**
        ```groovy linenums="26"
        echo "Processing: ${sample}"     // 오류: sample_id여야 합니다
        cat ${input_file} > ${sample}_result.txt  // 오류: sample_id여야 합니다
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

        이 시점에서 워크플로우가 실행되지만 여전히 오류가 발생합니다(예: `processFiles`에서 `Path value cannot be null`), 이는 잘못된 채널 구조로 인한 것입니다.

        **오류 5: 채널 구조 오류 - 잘못된 Map 출력**
        ```groovy linenums="83"
        .map { row -> row.sample_id }  // 오류: processFiles는 튜플을 예상합니다
        ```
        **수정:** processFiles가 예상하는 튜플 구조 반환
        ```groovy linenums="83"
        .map { row -> [row.sample_id, file(row.fastq_path)] }
        ```

        하지만 이것은 위의 `heavyProcess()` 실행을 위한 수정을 깨뜨리므로 해당 프로세스에 샘플 ID만 전달하기 위해 map을 사용해야 합니다:

        **오류 6: heavyProcess의 잘못된 채널 구조**
        ```groovy linenums="87"
        heavy_ch = heavyProcess(input_ch)  // 오류: input_ch는 이제 방출당 2개의 요소를 가집니다 - heavyProcess는 1개(첫 번째)만 필요합니다
        ```
        **수정:** 올바른 채널 사용 및 샘플 ID 추출
        ```groovy linenums="87"
        heavy_ch = heavyProcess(input_ch.map{it[0]})
        ```

        이제 조금 더 진행되지만 Bash 변수를 이스케이프하지 않았기 때문에 `No such variable: i`에 대한 오류를 받습니다.

        **오류 7: Bash 변수 이스케이프 오류**
        ```groovy linenums="48"
        echo "Heavy computation $i for ${sample_id}"  // 오류: $i가 이스케이프되지 않음
        ```
        **수정:** bash 변수 이스케이프
        ```groovy linenums="48"
        echo "Heavy computation \${i} for ${sample_id}"
        ```

        이제 `Process exceeded running time limit (1ms)`를 받으므로 관련 프로세스의 실행 시간 제한을 수정합니다:

        **오류 8: 리소스 구성 오류**
        ```groovy linenums="36"
        time '1 ms'  // 오류: 비현실적인 시간 제한
        ```
        **수정:** 현실적인 시간 제한으로 증가
        ```groovy linenums="36"
        time '100 s'
        ```

        다음으로 해결할 `Missing output file(s)` 오류가 있습니다:

        **오류 9: 출력 파일 이름 불일치**
        ```groovy linenums="49"
        done > ${sample_id}.txt  // 오류: 잘못된 파일 이름, 출력 선언과 일치해야 합니다
        ```
        **수정:** 출력 선언과 일치
        ```groovy linenums="49"
        done > ${sample_id}_heavy.txt
        ```

        처음 두 프로세스는 실행되었지만 세 번째는 실행되지 않았습니다.

        **오류 10: 출력 파일 이름 불일치**
        ```groovy linenums="88"
        file_ch = channel.fromPath("*.txt") // 오류: 프로세스가 아닌 pwd에서 입력을 가져오려고 시도
        handleFiles(file_ch)
        ```
        **수정:** 이전 프로세스의 출력 가져오기
        ```groovy linenums="88"
        handleFiles(heavyProcess.out)
        ```

        이것으로 전체 워크플로우가 실행되어야 합니다.

        **완전히 수정된 워크플로우:**
        ```groovy linenums="1"
        #!/usr/bin/env nextflow

        /*
        * 디버깅 연습을 위한 버그가 있는 워크플로우
        * 이 워크플로우에는 학습 목적을 위한 여러 의도적인 버그가 포함되어 있습니다
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

            // 잘못된 사용이 있는 채널
            input_ch = channel
                .fromPath(params.input)
                .splitCsv(header: true)
                .map { row -> [row.sample_id, file(row.fastq_path)] }

            processed_ch = processFiles(input_ch)

            heavy_ch = heavyProcess(input_ch.map{it[0]})

            handleFiles(heavyProcess.out)
        }
        ```

**다루는 오류 범주:**

- **구문 오류**: 중괄호 누락, 후행 쉼표, 정의되지 않은 변수
- **채널 구조 오류**: 잘못된 데이터 형태, 정의되지 않은 채널
- **프로세스 오류**: 출력 파일 불일치, 변수 이스케이프
- **리소스 오류**: 비현실적인 시간 제한

**주요 디버깅 교훈:**

1. **오류 메시지를 주의 깊게 읽기** - 종종 문제를 직접 가리킵니다
2. **체계적인 접근법 사용** - 한 번에 하나의 오류를 수정하고 `-resume`으로 테스트
3. **데이터 흐름 이해** - 채널 구조 오류가 종종 가장 미묘합니다
4. **작업 디렉토리 확인** - 프로세스가 실패하면 로그가 정확히 무엇이 잘못되었는지 알려줍니다

---

## 요약

이 사이드 퀘스트에서 Nextflow 워크플로우를 디버깅하기 위한 체계적인 기법 세트를 학습했습니다.
자신의 작업에 이러한 기법을 적용하면 컴퓨터와 싸우는 데 소비하는 시간을 줄이고, 문제를 더 빠르게 해결하며, 향후 문제로부터 자신을 보호할 수 있습니다.

### 주요 패턴

**1. 구문 오류를 식별하고 수정하는 방법**:

- Nextflow 오류 메시지를 해석하고 문제 위치 찾기
- 일반적인 구문 오류: 중괄호 누락, 잘못된 키워드, 정의되지 않은 변수
- Nextflow(Groovy) 변수와 Bash 변수 구별
- 조기 오류 감지를 위한 VS Code 확장 기능 사용

```groovy
// 중괄호 누락 - IDE에서 빨간 밑줄 찾기
process FOO {
    script:
    """
    echo "hello"
    """
// } <-- 누락!

// 잘못된 키워드
inputs:  // 'input:'이어야 합니다

// 정의되지 않은 변수 - Bash 변수의 경우 백슬래시로 이스케이프
echo "${undefined_var}"      // Nextflow 변수 (정의되지 않은 경우 오류)
echo "\${bash_var}"          // Bash 변수 (이스케이프됨)
```

**2. 채널 구조 문제를 디버깅하는 방법**:

- 채널 카디널리티 및 소진 문제 이해
- 채널 내용 구조 불일치 디버깅
- 채널 검사를 위한 `.view()` 연산자 사용
- 예상치 못한 튜플 구조를 나타내는 대괄호와 같은 오류 패턴 인식

````groovy
// 채널 내용 검사
my_channel.view { "Content: $it" }

// 큐를 값 채널로 변환 (소진 방지)
reference_ch = channel.value('ref.fa')
//```groovy
// 큐를 값 채널로 변환 (소진 방지)
reference_ch = channel.value('ref.fa')
// 또는
reference_ch = channel.of('ref.fa').first()
````

**3. 프로세스 실행 문제를 해결하는 방법**:

- 출력 파일 누락 오류 진단
- 종료 코드 이해 (소프트웨어 누락의 경우 127, 메모리 문제의 경우 137)
- 작업 디렉토리 및 명령 파일 조사
- 리소스를 적절하게 구성

```bash
# 실제로 실행된 것 확인
cat work/ab/cdef12/.command.sh

# 오류 출력 확인
cat work/ab/cdef12/.command.err

# 종료 코드 127 = 명령을 찾을 수 없음
# 종료 코드 137 = 종료됨 (메모리/시간 제한)
```

**4. Nextflow의 내장 디버깅 도구를 사용하는 방법**:

- 미리보기 모드 및 실시간 디버깅 활용
- 로직 테스트를 위한 스텁 실행 구현
- 효율적인 디버깅 사이클을 위한 resume 적용
- 4단계 체계적 디버깅 방법론 따르기

!!! tip "빠른 디버깅 참조"

    **구문 오류?** → VSCode 경고 확인, `nextflow run workflow.nf -preview` 실행

    **채널 문제?** → `.view()`를 사용하여 내용 검사: `my_channel.view()`

    **프로세스 실패?** → 작업 디렉토리 파일 확인:

    - `.command.sh` - 실행된 스크립트
    - `.command.err` - 오류 메시지
    - `.exitcode` - 종료 상태 (127 = 명령을 찾을 수 없음, 137 = 종료됨)

    **이상한 동작?** → `-stub-run`으로 워크플로우 로직 테스트

    **수정했나요?** → `-resume`을 사용하여 테스트 시간 절약: `nextflow run workflow.nf -resume`

---

### 추가 자료

- [Nextflow 문제 해결 가이드](https://www.nextflow.io/docs/latest/troubleshooting.html): 공식 문제 해결 문서
- [Nextflow 채널 이해](https://www.nextflow.io/docs/latest/channel.html): 채널 유형 및 동작에 대한 심층 분석
- [프로세스 지시문 참조](https://www.nextflow.io/docs/latest/process.html#directives): 사용 가능한 모든 프로세스 구성 옵션
- [nf-test](https://www.nf-test.com/): Nextflow 파이프라인을 위한 테스트 프레임워크
- [Nextflow Slack 커뮤니티](https://www.nextflow.io/slack-invite.html): 커뮤니티로부터 도움 받기

프로덕션 워크플로우의 경우 다음을 고려하세요:

- 대규모 모니터링 및 디버깅을 위한 [Seqera Platform](https://seqera.io/platform/) 설정
- 재현 가능한 소프트웨어 환경을 위한 [Wave 컨테이너](https://seqera.io/wave/) 사용

**기억하세요:** 효과적인 디버깅은 연습으로 향상되는 기술입니다. 여기서 습득한 체계적인 방법론과 포괄적인 도구 키트는 Nextflow 개발 여정 전반에 걸쳐 여러분에게 큰 도움이 될 것입니다.

---

## 다음 단계

[사이드 퀘스트 메뉴](./index.md)로 돌아가거나 페이지 오른쪽 하단의 버튼을 클릭하여 목록의 다음 주제로 이동하세요.
