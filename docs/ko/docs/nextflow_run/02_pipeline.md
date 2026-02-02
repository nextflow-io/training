# 파트 2: 실제 pipeline 실행

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } AI 지원 번역 - [자세히 알아보기 및 개선 제안](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

이 과정의 Part 1(기본 작업 실행)에서 코드 복잡성을 낮추기 위해 최소한의 기능만 있는 예제 workflow로 시작했습니다.
예를 들어 `1-hello.nf`는 명령줄 매개변수(`--input`)를 사용하여 한 번에 하나의 값을 제공했습니다.

그러나 대부분의 실제 pipeline은 대규모 데이터를 효율적으로 처리하고 때로는 복잡한 로직으로 연결된 여러 처리 단계를 적용하기 위해 더 정교한 기능을 사용합니다.

교육의 이 부분에서는 원래 Hello World pipeline의 확장된 버전을 시험해 보며 실제 pipeline의 핵심 기능을 시연합니다.

## 1. 파일에서 입력 데이터 처리

실제 pipeline에서는 일반적으로 하나 이상의 입력 파일에 포함된 여러 데이터 포인트(또는 데이터 시리즈)를 처리하려고 합니다.
그리고 가능하면 분석 대기 시간을 줄이기 위해 독립적인 데이터의 처리를 병렬로 실행하려고 합니다.

Nextflow가 이를 어떻게 수행하는지 보여주기 위해 실제 데이터 분석에서 처리할 수 있는 열 형식 데이터를 모방하여 여러 입력 인사말이 포함된 `greetings.csv`라는 CSV 파일을 준비했습니다.
숫자는 의미가 없으며 설명 목적으로만 있습니다.

```csv title="data/greetings.csv" linenums="1"
Hello,English,123
Bonjour,French,456
Holà,Spanish,789
```

원본 workflow의 개선된 버전인 `2a-inputs.nf`도 작성했습니다. 이 버전은 CSV 파일을 읽어 인사말을 추출하고 각각을 별도의 파일에 씁니다.

<figure class="excalidraw">
--8<-- "docs/en/docs/nextflow_run/img/hello-pipeline-multi-inputs.svg"
</figure>

먼저 workflow를 실행한 다음 관련 Nextflow 코드를 살펴보겠습니다.

### 1.1. Workflow 실행

터미널에서 다음 명령을 실행하세요.

```bash
nextflow run 2a-inputs.nf --input data/greetings.csv
```

??? success "명령 출력"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `2a-inputs.nf` [mighty_sammet] DSL2 - revision: 29fb5352b3

    executor >  local (3)
    [8e/0eb066] sayHello (2) [100%] 3 of 3 ✔
    ```

흥미롭게도 이것은 process에 대해 '3 of 3' 호출이 이루어졌음을 나타내며, 입력으로 제공한 CSV에 세 개의 데이터 행이 있었기 때문에 고무적입니다.
이것은 `sayHello()` process가 각 입력 행에 대해 한 번씩 세 번 호출되었음을 시사합니다.

### 1.2. `results` 디렉토리에서 게시된 출력 찾기

'results' 디렉토리를 살펴보고 workflow가 여전히 출력 사본을 거기에 쓰고 있는지 확인해 보겠습니다.

??? abstract "디렉토리 내용"

    ```console linenums="1" hl_lines="4-7"
    results
    ├── 1-hello
    |   └── output.txt
    └── 2a-inputs
        ├── Bonjour-output.txt
        ├── Hello-output.txt
        └── Holà-output.txt
    ```

네! 편리하게도 다른 이름을 가진 세 개의 출력 파일이 있는 `2a-inputs`라는 새 디렉토리가 보입니다.

각각을 열어 적절한 인사말 문자열이 포함되어 있는지 확인할 수 있습니다.

??? abstract "파일 내용"

    ```console title="results/2a-inputs/Hello-output.txt"
    Hello
    ```

    ```console title="results/2a-inputs/Bonjour-output.txt"
    Bonjour
    ```

    ```console title="results/2a-inputs/Holà-output.txt"
    Holà
    ```

이것은 입력 파일의 각 인사말이 적절하게 처리되었음을 확인합니다.

### 1.3. 원본 출력 및 로그 찾기

위의 콘솔 출력에서 하나의 작업 디렉토리만 참조된 것을 눈치채셨을 것입니다.
그것은 세 번의 `sayHello()` 호출이 모두 하나의 작업 디렉토리 내에서 실행되었다는 것을 의미하나요?

#### 1.3.1. 터미널에 제공된 작업 디렉토리 검토

`8e/0eb066` 작업 디렉토리 내부를 살펴보겠습니다.

??? abstract "디렉토리 내용"

    ```console title="8e/0eb066"
    work/8e/0eb066071cdb4123906b7b4ea8b047/
    └── Bonjour-output.txt
    ```

인사말 중 하나에 해당하는 출력만 찾을 수 있습니다(숨김 파일 표시를 활성화하면 부속 파일도 볼 수 있습니다).

무슨 일이 일어나고 있는 걸까요?

기본적으로 ANSI 로깅 시스템은 동일한 process에 대한 모든 호출의 상태 정보를 같은 줄에 씁니다.
결과적으로 콘솔 출력에서 세 개의 작업 디렉토리 경로 중 하나만(`8e/0eb066`) 표시했습니다.
거기에 나열되지 않은 두 개가 더 있습니다.

#### 1.3.2. 터미널에 더 많은 세부 정보 표시

다음과 같이 명령에 `-ansi-log false`를 추가하여 process 호출의 전체 목록을 보도록 로깅 동작을 수정할 수 있습니다:

```bash
nextflow run 2a-inputs.nf --input data/greetings.csv -ansi-log false
```

??? success "명령 출력"

    ```console linenums="1"
    N E X T F L O W  ~  version 25.10.2
    Launching `2a-inputs.nf` [pedantic_hamilton] DSL2 - revision: 6bbc42e49f
    [ab/1a8ece] Submitted process > sayHello (1)
    [0d/2cae24] Submitted process > sayHello (2)
    [b5/0df1d6] Submitted process > sayHello (3)
    ```

이번에는 출력에 세 가지 process 실행과 관련 작업 하위 디렉토리가 모두 나열됩니다.
ANSI 로깅을 비활성화하면 Nextflow가 터미널 출력에서 색상을 사용하지 못하게 됩니다.

두 로깅 모드 간에 상태 보고 방식이 약간 다릅니다.
압축 모드에서 Nextflow는 호출이 성공적으로 완료되었는지 여부를 보고합니다.
이 확장 모드에서는 제출되었다는 것만 보고합니다.

이것은 `sayHello()` process가 세 번 호출되고 각각에 대해 별도의 작업 디렉토리가 생성된다는 것을 확인합니다.

거기에 나열된 각 작업 디렉토리 내부를 살펴보면 각각이 인사말 중 하나에 해당하는지 확인할 수 있습니다.

??? abstract "디렉토리 내용"

    ```console title="ab/1a8ece"
    work/ab/1a8ece307e53f03fce689dde904b64/
    └── Hello-output.txt
    ```

    ```console title="0d/2cae24"
    work/0d/2cae2481a53593bc607077c80c9466/
    └── Bonjour-output.txt
    ```

    ```console title="b5/0df1d6"
    work/b5/0df1d642353269909c2ce23fc2a8fa/
    └── Holà-output.txt
    ```

이것은 각 process 호출이 다른 모든 것과 격리되어 실행된다는 것을 확인합니다.
이것은 process가 고유하지 않은 이름을 가진 중간 파일을 생성하는 경우 충돌을 피하는 것을 포함하여 많은 장점이 있습니다.

!!! tip

    복잡한 workflow나 많은 수의 입력의 경우 전체 목록을 터미널에 출력하면 약간 압도적일 수 있으므로 사람들은 일반적으로 일상적인 사용에서 `-ansi-log false`를 사용하지 않습니다.

### 1.4. Workflow 코드 검토

이 버전의 workflow는 입력의 CSV 파일을 읽고 입력을 별도로 처리하고 출력 이름을 고유하게 지정할 수 있습니다.

workflow 코드에서 무엇이 이를 가능하게 하는지 살펴보겠습니다.

??? full-code "전체 코드 파일"

    ```groovy title="2a-inputs.nf" linenums="1" hl_lines="31-33 35"
    #!/usr/bin/env nextflow

    /*
    * Use echo to print 'Hello World!' to a file
    */
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

    /*
    * Pipeline parameters
    */
    params {
        input: Path
    }

    workflow {

        main:
        // 입력용 채널 생성
        greeting_ch = channel.fromPath(params.input)
                            .splitCsv()
                            .map { line -> line[0] }
        // 인사말 출력
        sayHello(greeting_ch)

        publish:
        first_output = sayHello.out
    }

    output {
        first_output {
            path '2a-inputs'
            mode 'copy'
        }
    }
    ```

다시 한번, 코드 구문을 암기할 필요는 없지만 중요한 기능을 제공하는 workflow의 핵심 구성 요소를 인식하는 방법을 배우는 것이 좋습니다.

#### 1.4.1. CSV에서 입력 데이터 로드

이것이 가장 흥미로운 부분입니다: 명령줄에서 단일 값을 가져오는 것에서 CSV 파일을 가져와 구문 분석하고 포함된 개별 인사말을 처리하는 것으로 어떻게 전환했을까요?

Nextflow에서는 **channel**로 이를 수행합니다. channel은 입력을 효율적으로 처리하고 다단계 workflow에서 한 단계에서 다른 단계로 셔틀하도록 설계된 구조로, 내장된 병렬 처리와 많은 추가 이점을 제공합니다.

분석해 보겠습니다.

```groovy title="2a-inputs.nf" linenums="29" hl_lines="3-5"
    main:
    // 입력용 채널 생성
    greeting_ch = channel.fromPath(params.input)
                        .splitCsv()
                        .map { line -> line[0] }
    // 인사말 출력
    sayHello(greeting_ch)
```

이 코드는 CSV 파일을 읽고 구문 분석하고 각 행에서 첫 번째 열을 추출하는 `greeting_ch`라는 channel을 생성합니다.
결과는 `Hello`, `Bonjour`, `Holà`를 포함하는 channel입니다.

??? tip "이것이 어떻게 작동하나요?"

    이 줄이 일반 영어로 무엇을 의미하는지 설명하면:

    - `channel.fromPath`는 파일 경로에서 channel을 생성하는 **channel factory**입니다
    - `(params.input)`은 파일 경로가 명령줄의 `--input`에 의해 제공됨을 지정합니다

    즉, 그 줄은 Nextflow에게 `--input`으로 제공된 파일 경로를 가져와 그 내용을 입력 데이터로 처리할 준비를 하라고 말합니다.

    그런 다음 다음 두 줄은 파일의 실제 구문 분석과 적절한 데이터 구조로의 데이터 로드를 수행하는 **연산자**를 적용합니다:

    - `.splitCsv()`는 Nextflow에게 CSV 파일을 행과 열을 나타내는 배열로 구문 분석하라고 말합니다
    - `.map { line -> line[0] }`는 Nextflow에게 각 행에서 첫 번째 열의 요소만 가져오라고 말합니다

    실제로 다음 CSV 파일에서 시작하여:

    ```csv title="greetings.csv" linenums="1"
    Hello,English,123
    Bonjour,French,456
    Holà,Spanish,789
    ```

    이것을 다음과 같은 배열로 변환했습니다:

    ```txt title="배열 내용"
    [[Hello,English,123],[Bonjour,French,456],[Holà,Spanish,789]]
    ```

    그런 다음 세 행 각각에서 첫 번째 요소를 가져와 이제 `Hello`, `Bonjour`, `Holà`를 포함하는 Nextflow channel에 로드했습니다.

    channel과 연산자를 깊이 이해하고 직접 작성하는 방법을 포함하여 [Hello Nextflow Part 2: Hello Channels](../hello_nextflow/02_hello_channels.md#4-read-input-values-from-a-csv-file)를 참조하세요.

#### 1.4.2. 각 인사말에 대해 process 호출

다음으로 workflow의 `main:` 블록의 마지막 줄에서 로드된 `greeting_ch` channel을 `sayHello()` process에 입력으로 제공합니다.

```groovy title="2a-inputs.nf" linenums="29" hl_lines="7"
    main:
    // 입력용 채널 생성
    greeting_ch = channel.fromPath(params.input)
                        .splitCsv()
                        .map { line -> line[0] }
    // 인사말 출력
    sayHello(greeting_ch)
```

이것은 Nextflow에게 channel의 각 요소, 즉 각 인사말에 대해 process를 개별적으로 실행하라고 말합니다.
그리고 Nextflow는 똑똑하기 때문에 사용 가능한 컴퓨팅 인프라에 따라 가능한 경우 이러한 process 호출을 병렬로 실행합니다.

이것이 상대적으로 매우 적은 코드로 많은 데이터(많은 샘플, 데이터 포인트 또는 연구 단위)의 효율적이고 확장 가능한 처리를 달성할 수 있는 방법입니다.

#### 1.4.3. 출력 이름 지정 방법

마지막으로 출력 파일의 이름이 고유하게 지정되는 방법을 보기 위해 process 코드를 간략히 살펴볼 가치가 있습니다.

```groovy title="2a-inputs.nf" linenums="6" hl_lines="7 11"
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

`1-hello.nf`의 이 process 버전과 비교하여 출력 선언과 명령의 관련 부분이 출력 파일 이름에 인사말 값을 포함하도록 변경된 것을 볼 수 있습니다.

이것은 출력 파일 이름이 공통 results 디렉토리에 게시될 때 충돌하지 않도록 하는 한 가지 방법입니다.

그리고 이것이 process 선언 내에서 변경해야 했던 유일한 변경 사항입니다!

### 요약

channel과 연산자가 여러 입력을 효율적으로 처리하는 방법을 기본적인 수준에서 이해합니다.

### 다음 단계

다단계 workflow가 어떻게 구성되고 작동하는지 알아봅니다.

---

## 2. 다단계 workflow 실행

대부분의 실제 workflow에는 둘 이상의 단계가 포함됩니다.
channel에 대해 방금 배운 것을 바탕으로 Nextflow가 channel과 연산자를 사용하여 다단계 workflow에서 process를 함께 연결하는 방법을 살펴보겠습니다.

이를 위해 세 개의 별도 단계를 연결하고 다음을 시연하는 예제 workflow를 제공합니다:

1. 데이터가 한 process에서 다음 process로 흐르도록 만들기
2. 여러 process 호출의 출력을 단일 process 호출로 수집

구체적으로 각 입력 인사말을 가져와 대문자로 변환한 다음 모든 대문자 인사말을 단일 출력 파일로 수집하는 `2b-multistep.nf`라는 workflow의 확장된 버전을 만들었습니다.

<figure class="excalidraw">
--8<-- "docs/en/docs/nextflow_run/img/hello-pipeline-multi-steps.svg"
</figure>

이전과 마찬가지로 먼저 workflow를 실행한 다음 코드를 살펴보며 새로운 것이 무엇인지 확인합니다.

### 2.1. Workflow 실행

터미널에서 다음 명령을 실행하세요:

```bash
nextflow run 2b-multistep.nf --input data/greetings.csv
```

??? success "명령 출력"

    ```console linenums="1"
    N E X T F L O W   ~  version 25.10.2

    Launching `2b-multistep.nf` [soggy_franklin] DSL2 - revision: bc8e1b2726

    [d6/cdf466] sayHello (1)       | 3 of 3 ✔
    [99/79394f] convertToUpper (2) | 3 of 3 ✔
    [1e/83586c] collectGreetings   | 1 of 1 ✔
    ```

약속대로 workflow의 일부로 여러 단계가 실행된 것을 볼 수 있습니다. 처음 두 개(`sayHello`와 `convertToUpper`)는 아마도 각 개별 인사말에서 실행되었고, 세 번째(`collectGreetings`)는 세 개의 `convertToUpper` 호출 모두의 출력에서 한 번만 실행되었을 것입니다.

### 2.2. 출력 찾기

`results` 디렉토리를 살펴보며 실제로 그런 일이 일어났는지 확인해 보겠습니다.

??? abstract "디렉토리 내용"

    ```console linenums="1" hl_lines="8-16"
    results
    ├── 1-hello
    |   └── output.txt
    ├── 2a-inputs
    |   ├── Bonjour-output.txt
    |   ├── Hello-output.txt
    |   └── Holà-output.txt
    └── 2b-multistep
        ├── COLLECTED-batch-output.txt
        ├── batch-report.txt
        └── intermediates
            ├── Bonjour-output.txt
            ├── Hello-output.txt
            ├── Holà-output.txt
            ├── UPPER-Bonjour-output.txt
            ├── UPPER-Hello-output.txt
            └── UPPER-Holà-output.txt

    ```

보시다시피 `2b-multistep`이라는 새 디렉토리가 있으며 이전보다 훨씬 더 많은 파일이 포함되어 있습니다.
일부 파일은 `intermediates`라는 하위 디렉토리에 그룹화되어 있고 두 개의 파일은 최상위 수준에 있습니다.

그 두 개가 다단계 workflow의 최종 결과입니다.
파일 이름을 살펴보고 내용을 확인하여 예상대로인지 확인하세요.

??? abstract "파일 내용"

    ```txt title="results/2b-multistep/COLLECTED-batch-output.txt"
    HELLO
    BONJOUR
    HOLà
    ```

    ```txt title="results/2b-multistep/batch-report.txt"
    There were 3 greetings in this batch.
    ```

첫 번째는 약속대로 세 개의 인사말이 대문자로 변환되어 단일 파일로 수집된 것을 포함합니다.
두 번째는 실행에 대한 일부 정보를 요약하는 보고서 파일입니다.

### 2.3. 코드 검토

코드를 살펴보고 다단계 workflow의 핵심 패턴을 식별해 보겠습니다.

??? full-code "전체 코드 파일"

    ```groovy title="2b-multistep.nf" linenums="1" hl_lines="63 75-78 82-84"
    #!/usr/bin/env nextflow

    /*
    * Use echo to print 'Hello World!' to a file
    */
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

    /*
    * Use a text replacement tool to convert the greeting to uppercase
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

    /*
    * Collect uppercase greetings into a single output file
    */
    process collectGreetings {

        input:
        path input_files
        val batch_name

        output:
        path "COLLECTED-${batch_name}-output.txt", emit: outfile
        path "${batch_name}-report.txt", emit: report

        script:
        count_greetings = input_files.size()
        """
        cat ${input_files} > 'COLLECTED-${batch_name}-output.txt'
        echo 'There were ${count_greetings} greetings in this batch.' > '${batch_name}-report.txt'
        """
    }

    /*
    * Pipeline parameters
    */
    params {
        input: Path
        batch: String = 'batch'
    }

    workflow {

        main:
        // 입력용 채널 생성
        greeting_ch = channel.fromPath(params.input)
                            .splitCsv()
                            .map { line -> line[0] }
        // 인사말 출력
        sayHello(greeting_ch)
        // 인사말을 대문자로 변환
        convertToUpper(sayHello.out)
        // 모든 인사말을 하나의 파일에 수집
        collectGreetings(convertToUpper.out.collect(), params.batch)

        publish:
        first_output = sayHello.out
        uppercased = convertToUpper.out
        collected = collectGreetings.out.outfile
        batch_report = collectGreetings.out.report
    }

    output {
        first_output {
            path '2b-multistep/intermediates'
            mode 'copy'
        }
        uppercased {
            path '2b-multistep/intermediates'
            mode 'copy'
        }
        collected {
            path '2b-multistep'
            mode 'copy'
        }
        batch_report {
            path '2b-multistep'
            mode 'copy'
        }
    }
    ```

많은 일이 일어나고 있지만 이전 버전의 workflow와 비교했을 때 가장 명백한 차이점은 이제 여러 process 정의가 있고 해당하여 workflow 블록에 여러 process 호출이 있다는 것입니다.

자세히 살펴보고 가장 흥미로운 부분을 식별할 수 있는지 확인해 보겠습니다.

#### 2.3.1. Workflow 구조 시각화

Nextflow 확장이 있는 VSCode를 사용하는 경우 Nextflow 스크립트의 workflow 블록 바로 위에 표시되는 작은 `DAG preview` 링크를 클릭하여 process가 어떻게 연결되어 있는지에 대한 유용한 다이어그램을 얻을 수 있습니다.

<figure class="excalidraw">
--8<-- "docs/en/docs/nextflow_run/img/DAG-multistep.svg"
</figure>

이것은 process가 어떻게 연결되어 있고 무엇을 생성하는지에 대한 좋은 개요를 제공합니다.

원래 `sayHello` process 외에도 이제 `convertToUpper`와 `collectGreetings`도 있으며, 이는 콘솔 출력에서 본 process 이름과 일치합니다.
두 개의 새 process 정의는 `sayHello` process와 동일한 방식으로 구조화되어 있지만, `collectGreetings`는 `batch`라는 추가 입력 매개변수를 받고 두 개의 출력을 생성합니다.

각각의 코드에 대해 자세히 다루지 않겠지만, 궁금하시다면 [Hello Nextflow의 Part 2](../hello_nextflow/03_hello_workflow.md)에서 세부 사항을 찾아볼 수 있습니다.

지금은 process가 서로 어떻게 연결되어 있는지 살펴보겠습니다.

#### 2.3.2. Process 연결 방법

여기서 정말 흥미로운 것은 workflow의 `main:` 블록에서 process 호출이 어떻게 연결되어 있는지입니다.

```groovy title="2b-multistep.nf" linenums="68" hl_lines="9 11"
    main:
    // 입력용 채널 생성
    greeting_ch = channel.fromPath(params.input)
                        .splitCsv()
                        .map { line -> line[0] }
    // 인사말 출력
    sayHello(greeting_ch)
    // 인사말을 대문자로 변환
    convertToUpper(sayHello.out)
    // 모든 인사말을 하나의 파일에 수집
    collectGreetings(convertToUpper.out.collect(), params.batch)
```

첫 번째 process 호출인 `sayHello(greeting_ch)`는 변경되지 않은 것을 볼 수 있습니다.
그런 다음 `convertToUpper`에 대한 다음 process 호출은 `sayHello`의 출력을 `sayHello.out`으로 참조합니다.

패턴은 간단합니다: `processName.out`은 process의 출력 channel을 참조하며, 이를 다음 process에 직접 전달할 수 있습니다.
이것이 Nextflow에서 한 단계에서 다음 단계로 데이터를 셔틀하는 방법입니다.

#### 2.3.3. Process는 여러 입력을 받을 수 있습니다

세 번째 process 호출인 `collectGreetings`는 약간 다릅니다.

```groovy title="2b-multistep.nf" linenums="77"
    // 모든 인사말을 하나의 파일에 수집
    collectGreetings(convertToUpper.out.collect(), params.batch)
```

이 호출에 `convertToUpper.out.collect()`와 `params.batch`라는 두 개의 입력이 제공되는 것을 볼 수 있습니다.
지금은 `.collect()` 부분을 무시하고 이것을 `collectGreetings(input1, input2)`로 일반화할 수 있습니다.

이것은 process 모듈의 두 입력 선언과 일치합니다:

```groovy title="2b-multistep.nf" linenums="40"
process collectGreetings {

    input:
    path input_files
    val batch_name
```

Nextflow가 이것을 구문 분석할 때 호출의 첫 번째 입력을 `path input_files`에 할당하고 두 번째를 `val batch_name`에 할당합니다.

이제 process가 여러 입력을 받을 수 있다는 것과 workflow 블록에서 호출이 어떻게 보이는지 알게 되었습니다.

이제 첫 번째 입력인 `convertToUpper.out.collect()`를 자세히 살펴보겠습니다.

#### 2.3.4. `collectGreetings` 호출에서 `collect()`가 하는 일

`sayHello`의 출력을 `convertToUpper`에 전달하기 위해 `sayHello.out`으로 `sayHello`의 출력 channel을 간단히 참조했습니다. 그러나 다음 단계에서는 `convertToUpper.out.collect()`에 대한 참조를 보고 있습니다.

이 `collect()` 부분은 무엇이고 무엇을 하나요?

물론 연산자입니다. 앞서 만난 `splitCsv` 및 `map` 연산자와 마찬가지입니다.
이번에는 연산자가 `collect`이라고 하며 `convertToUpper`가 생성한 출력 channel에 적용됩니다.

`collect` 연산자는 동일한 process에 대한 여러 호출의 출력을 수집하고 단일 channel 요소로 패키징하는 데 사용됩니다.

이 workflow의 맥락에서 `convertToUpper.out` channel의 세 개의 대문자 인사말(세 개의 별도 channel 항목이며 일반적으로 다음 process에 의해 별도의 호출로 처리됨)을 가져와 단일 항목으로 패키징합니다.

더 실용적인 용어로: `convertToUpper()`의 출력을 `collectGreetings()`에 공급하기 전에 `collect()`를 적용하지 않으면 Nextflow는 단순히 각 인사말에 대해 `collectGreetings()`를 독립적으로 실행하여 목표를 달성하지 못합니다.

<figure class="excalidraw">
--8<-- "docs/en/docs/nextflow_run/img/without-collect-operator.svg"
</figure>

반대로 `collect()`를 사용하면 workflow의 두 번째 단계에서 생성된 모든 별도의 대문자 인사말을 가져와 pipeline의 세 번째 단계에서 단일 호출에 모두 함께 공급할 수 있습니다.

<figure class="excalidraw">
--8<-- "docs/en/docs/nextflow_run/img/with-collect-operator.svg"
</figure>

이것이 모든 인사말을 같은 파일로 다시 가져오는 방법입니다.

process 호출 간에 channel 내용에 변환을 적용할 수 있는 다른 많은 [연산자](https://www.nextflow.io/docs/latest/reference/operator.html#operator-page)가 있습니다.

이것은 pipeline 개발자에게 pipeline의 흐름 로직을 사용자 정의하는 데 많은 유연성을 제공합니다.
단점은 때때로 pipeline이 무엇을 하고 있는지 해독하기 어렵게 만들 수 있다는 것입니다.

#### 2.3.5. 입력 매개변수는 기본값을 가질 수 있습니다

`collectGreetings`가 두 번째 입력인 `params.batch`를 받는 것을 눈치채셨을 것입니다:

```groovy title="2b-multistep.nf" linenums="77"
    // 모든 인사말을 하나의 파일에 수집
    collectGreetings(convertToUpper.out.collect(), params.batch)
```

이것은 `--batch`라는 CLI 매개변수를 workflow에 전달합니다.
그러나 이전에 workflow를 시작할 때 `--batch` 매개변수를 지정하지 않았습니다.

무슨 일이 일어나고 있을까요?
`params` 블록을 살펴보세요:

```groovy title="2b-multistep.nf" linenums="61" hl_lines="3"
params {
    input: Path
    batch: String = 'batch'
}
```

workflow에 기본값이 구성되어 있으므로 제공할 필요가 없습니다.
하지만 명령줄에서 제공하면 지정한 값이 기본값 대신 사용됩니다.

시도해 보세요:

```bash
nextflow run 2b-multistep.nf --input data/greetings.csv --batch test
```

??? success "명령 출력"

    ```console linenums="1"
    N E X T F L O W   ~  version 25.10.2

    Launching `2b-multistep.nf` [soggy_franklin] DSL2 - revision: bc8e1b2726

    [a5/cdff26] sayHello (1)       | 3 of 3 ✔
    [c5/78794f] convertToUpper (2) | 3 of 3 ✔
    [d3/b4d86c] collectGreetings   | 1 of 1 ✔
    ```

사용자 정의 배치 이름으로 명명된 새 최종 출력을 볼 수 있습니다.

??? abstract "디렉토리 내용"

    ```console linenums="1" hl_lines="10 12"
    results
    ├── 1-hello
    |   └── output.txt
    ├── 2a-inputs
    |   ├── Bonjour-output.txt
    |   ├── Hello-output.txt
    |   └── Holà-output.txt
    └── 2b-multistep
        ├── COLLECTED-batch-output.txt
        ├── COLLECTED-test-output.txt
        ├── batch-report.txt
        ├── test-report.txt
        └── intermediates
            ├── Bonjour-output.txt
            ├── Hello-output.txt
            ├── Holà-output.txt
            ├── UPPER-Bonjour-output.txt
            ├── UPPER-Hello-output.txt
            └── UPPER-Holà-output.txt
    ```

이것은 Part 3에서 더 자세히 다룰 입력 구성의 한 측면이지만, 지금 중요한 것은 입력 매개변수에 기본값을 부여할 수 있다는 것입니다.

#### 2.3.6. Process는 여러 출력을 생성할 수 있습니다

`collectGreetings` process 정의에서 다음 출력 선언을 볼 수 있습니다:

```groovy title="2b-multistep.nf" linenums="46"
    output:
    path "COLLECTED-${batch_name}-output.txt", emit: outfile
    path "${batch_name}-report.txt", emit: report
```

이것은 `publish:` 블록에서 `emit:`으로 지정된 이름으로 참조됩니다:

```groovy title="2b-multistep.nf" linenums="80" hl_lines="4 5"
    publish:
    first_output = sayHello.out
    uppercased = convertToUpper.out
    collected = collectGreetings.out.outfile
    batch_report = collectGreetings.out.report
```

이렇게 하면 다양한 연산자와 함께 workflow의 다른 process에 특정 출력을 개별적으로 쉽게 전달할 수 있습니다.

#### 2.3.7. 게시된 출력을 구성할 수 있습니다

`output` 블록에서 workflow의 최종 출력만 쉽게 선택할 수 있도록 중간 결과를 그룹화하기 위해 사용자 정의 경로를 사용했습니다.

```groovy title="2b-multistep.nf" linenums="87" hl_lines="3 7 11 15"
output {
    first_output {
        path '2b-multistep/intermediates'
        mode 'copy'
    }
    uppercased {
        path '2b-multistep/intermediates'
        mode 'copy'
    }
    collected {
        path '2b-multistep'
        mode 'copy'
    }
    batch_report {
        path '2b-multistep'
        mode 'copy'
    }
}
```

게시된 출력을 구성하는 더 정교한 방법이 있습니다. 구성에 관한 부분에서 몇 가지를 다룰 것입니다.

!!! tip "Workflow 구축에 대해 더 알고 싶으신가요?"

    다단계 workflow 구축에 대한 자세한 내용은 [Hello Nextflow Part 3: Hello Workflow](../hello_nextflow/03_hello_workflow.md)를 참조하세요.

### 요약

channel과 연산자를 사용하여 다단계 workflow가 어떻게 구성되고 작동하는지 기본적인 수준에서 이해합니다.
또한 process가 여러 입력을 받고 여러 출력을 생성할 수 있으며 이것들이 구조화된 방식으로 게시될 수 있음을 보았습니다.

### 다음 단계

Nextflow pipeline이 코드 재사용과 유지 관리성을 촉진하기 위해 어떻게 모듈화될 수 있는지 배웁니다.

---

## 3. 모듈화된 pipeline 실행

지금까지 살펴본 모든 workflow는 모든 관련 코드를 포함하는 하나의 단일 workflow 파일로 구성되었습니다.

그러나 실제 pipeline은 일반적으로 *모듈화*의 이점을 누리며, 이는 코드가 다른 파일로 분할됨을 의미합니다.
이렇게 하면 개발 및 유지 관리가 더 효율적이고 지속 가능해질 수 있습니다.

여기서는 Nextflow에서 가장 일반적인 코드 모듈성 형태인 **모듈** 사용을 시연하겠습니다.

Nextflow에서 **모듈**은 독립 실행형 코드 파일에 자체적으로 캡슐화된 단일 process 정의입니다.
workflow에서 모듈을 사용하려면 workflow 코드 파일에 한 줄 import 문을 추가하기만 하면 됩니다. 그런 다음 일반적으로 하는 것과 같은 방식으로 process를 workflow에 통합할 수 있습니다.
이를 통해 코드의 여러 복사본을 생성하지 않고도 여러 workflow에서 process 정의를 재사용할 수 있습니다.

지금까지 모든 process가 단일 코드 파일에 포함된 workflow를 실행해 왔습니다.
이제 process가 개별 모듈에 저장될 때 어떻게 보이는지 살펴보겠습니다.

물론 `modules/` 디렉토리에 있는 모듈 세트와 함께 시연 목적으로 `2c-modules.nf`라는 적합한 workflow를 준비했습니다.

<figure class="excalidraw">
--8<-- "docs/en/docs/nextflow_run/img/modules.svg"
</figure>

??? abstract "디렉토리 내용"

    ```console
    modules/
    ├── collectGreetings.nf
    ├── convertToUpper.nf
    ├── cowpy.nf
    └── sayHello.nf
    ```

각각 process 중 하나의 이름을 딴 네 개의 Nextflow 파일이 있습니다.
지금은 `cowpy.nf` 파일을 무시해도 됩니다. 나중에 다룰 것입니다.

### 3.1. 코드 검토

이번에는 먼저 코드를 살펴보겠습니다.
`2c-modules.nf` workflow 파일을 열어 시작하세요.

??? full-code "전체 코드 파일"

    ```groovy title="2c-modules.nf" linenums="1"
    #!/usr/bin/env nextflow

    // 모듈 포함
    include { sayHello } from './modules/sayHello.nf'
    include { convertToUpper } from './modules/convertToUpper.nf'
    include { collectGreetings } from './modules/collectGreetings.nf'

    /*
    * Pipeline parameters
    */
    params {
        input: Path
        batch: String = 'batch'
    }

    workflow {

        main:
        // 입력용 채널 생성
        greeting_ch = channel.fromPath(params.input)
                            .splitCsv()
                            .map { line -> line[0] }
        // 인사말 출력
        sayHello(greeting_ch)
        // 인사말을 대문자로 변환
        convertToUpper(sayHello.out)
        // 모든 인사말을 하나의 파일에 수집
        collectGreetings(convertToUpper.out.collect(), params.batch)

        publish:
        first_output = sayHello.out
        uppercased = convertToUpper.out
        collected = collectGreetings.out.outfile
        batch_report = collectGreetings.out.report
    }

    output {
        first_output {
            path '2c-modules/intermediates'
            mode 'copy'
        }
        uppercased {
            path '2c-modules/intermediates'
            mode 'copy'
        }
        collected {
            path '2c-modules'
            mode 'copy'
        }
        batch_report {
            path '2c-modules'
            mode 'copy'
        }
    }
    ```

workflow 로직이 이전 버전의 workflow와 정확히 동일한 것을 볼 수 있습니다.
그러나 process 코드가 workflow 파일에서 사라지고 대신 `modules` 아래의 별도 파일을 가리키는 `include` 문이 있습니다.

```groovy title="hello-modules.nf" linenums="3"
// 모듈 포함
include { sayHello } from './modules/sayHello.nf'
include { convertToUpper } from './modules/convertToUpper.nf'
include { collectGreetings } from './modules/collectGreetings.nf'
```

해당 파일 중 하나를 열면 해당 process의 코드를 찾을 수 있습니다.

??? full-code "전체 코드 파일"

    ```groovy title="modules/sayHello.nf" linenums="1"
    #!/usr/bin/env nextflow

    /*
    * Use echo to print 'Hello World!' to a file
    */
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

보시다시피 process 코드는 변경되지 않았습니다. 기본 workflow 파일에 있는 대신 개별 모듈 파일에 복사되었을 뿐입니다.
다른 두 process에도 동일하게 적용됩니다.

이제 이 새 버전을 실행하면 어떻게 보이는지 살펴보겠습니다.

### 3.2. Workflow 실행

`-resume` 플래그와 함께 터미널에서 이 명령을 실행하세요:

```bash
nextflow run 2c-modules.nf --input data/greetings.csv -resume
```

??? success "명령 출력"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `2c-modules.nf` [soggy_franklin] DSL2 - revision: bc8e1b2726

    [j6/cdfa66] sayHello (1)       | 3 of 3, cached: ✔
    [95/79484f] convertToUpper (2) | 3 of 3, cached: ✔
    [5e/4358gc] collectGreetings   | 1 of 1, cached: ✔
    ```

코드가 분할되고 기본 workflow 파일 이름이 변경되었음에도 불구하고 process 실행이 모두 성공적으로 캐시되었음을 알 수 있습니다. 이는 Nextflow가 요청된 작업을 이미 수행했음을 인식했다는 것을 의미합니다.

Nextflow에게 중요한 것은 모든 코드가 함께 가져와지고 평가된 후 생성되는 작업 스크립트입니다.

!!! tip

    workflow의 한 섹션을 더 큰 pipeline에 가져올 수 있는 '하위 workflow'로 캡슐화하는 것도 가능하지만, 이 과정의 범위를 벗어납니다.

    [Workflows of Workflows](https://training.nextflow.io/latest/side_quests/workflows_of_workflows/)에 대한 Side Quest에서 조합 가능한 workflow 개발에 대해 더 자세히 알아볼 수 있습니다.

### 요약

process가 코드 재사용을 촉진하고 유지 관리성을 개선하기 위해 독립 실행형 모듈에 저장될 수 있음을 알게 되었습니다.

### 다음 단계

소프트웨어 의존성을 관리하기 위해 컨테이너를 사용하는 방법을 배웁니다.

---

## 4. 컨테이너화된 소프트웨어 사용

지금까지 예제로 사용한 workflow는 환경에서 사용 가능한 UNIX 도구를 사용하여 매우 기본적인 텍스트 처리 작업만 실행하면 되었습니다.

그러나 실제 pipeline은 일반적으로 대부분의 환경에 기본적으로 포함되지 않은 전문 도구와 패키지가 필요합니다.
일반적으로 이러한 도구를 설치하고 의존성을 관리하고 충돌을 해결해야 합니다.

이 모든 것은 매우 지루하고 짜증나는 일입니다.
이 문제를 해결하는 훨씬 더 좋은 방법은 **컨테이너**를 사용하는 것입니다.

**컨테이너**는 컨테이너 **이미지**에서 생성된 경량의 독립 실행형 실행 가능한 소프트웨어 단위로, 코드, 시스템 라이브러리 및 설정을 포함하여 애플리케이션을 실행하는 데 필요한 모든 것을 포함합니다.

!!! Tip

    [Docker](https://www.docker.com/get-started/) 기술을 사용하여 가르치지만, Nextflow는 [다른 여러 컨테이너 기술](https://www.nextflow.io/docs/latest/container.html#)도 지원합니다.

### 4.1. 컨테이너 직접 사용

먼저 컨테이너와 직접 상호 작용해 보겠습니다.
이렇게 하면 Nextflow에서 사용하기 전에 컨테이너가 무엇인지에 대한 이해를 확고히 하는 데 도움이 됩니다.

#### 4.1.1. 컨테이너 이미지 가져오기

컨테이너를 사용하려면 일반적으로 컨테이너 레지스트리에서 컨테이너 이미지를 다운로드하거나 "가져온" 다음 컨테이너 이미지를 실행하여 컨테이너 인스턴스를 생성합니다.

일반 구문은 다음과 같습니다:

```bash title="구문"
docker pull '<container>'
```

- `docker pull`은 컨테이너 시스템에 저장소에서 컨테이너 이미지를 가져오라는 지시입니다.
- `'<container>'`는 컨테이너 이미지의 URI 주소입니다.

예를 들어 임의의 텍스트 입력을 재미있는 방식으로 표시하기 위해 ASCII 아트를 생성하는 `cowsay`라는 도구의 Python 구현인 [cowpy](https://github.com/jeffbuttars/cowpy)가 포함된 컨테이너 이미지를 가져와 보겠습니다.

게시된 컨테이너를 찾을 수 있는 다양한 저장소가 있습니다.
[Seqera Containers](https://seqera.io/containers/) 서비스를 사용하여 `cowpy` Conda 패키지에서 이 Docker 컨테이너 이미지를 생성했습니다: `'community.wave.seqera.io/library/cowpy:1.1.5--3db457ae1977a273'`.

전체 가져오기 명령을 실행하세요:

```bash
docker pull 'community.wave.seqera.io/library/cowpy:1.1.5--3db457ae1977a273'
```

??? success "명령 출력"

    ```console
    Unable to find image 'community.wave.seqera.io/library/cowpy:1.1.5--3db457ae1977a273' locally
    131d6a1b707a8e65: Pulling from library/cowpy
    dafa2b0c44d2: Pull complete
    dec6b097362e: Pull complete
    f88da01cff0b: Pull complete
    4f4fb700ef54: Pull complete
    92dc97a3ef36: Pull complete
    403f74b0f85e: Pull complete
    10b8c00c10a5: Pull complete
    17dc7ea432cc: Pull complete
    bb36d6c3110d: Pull complete
    0ea1a16bbe82: Pull complete
    030a47592a0a: Pull complete
    622dd7f15040: Pull complete
    895fb5d0f4df: Pull complete
    Digest: sha256:fa50498b32534d83e0a89bb21fec0c47cc03933ac95c6b6587df82aaa9d68db3
    Status: Downloaded newer image for community.wave.seqera.io/library/cowpy:1.1.5--3db457ae1977a273
    community.wave.seqera.io/library/cowpy:1.1.5--3db457ae1977a273
    ```

이것은 시스템에 지정된 이미지를 다운로드하라고 지시합니다.
다운로드가 완료되면 컨테이너 이미지의 로컬 복사본이 생깁니다.

#### 4.1.2. 컨테이너 시작

컨테이너는 일회성 명령으로 실행할 수 있지만, 대화형으로 사용할 수도 있어 컨테이너 내부에서 셸 프롬프트를 제공하고 명령으로 놀 수 있습니다.

일반 구문은 다음과 같습니다:

```bash title="구문"
docker run --rm '<container>' [tool command]
```

- `docker run --rm '<container>'`는 컨테이너 시스템에 컨테이너 이미지에서 컨테이너 인스턴스를 시작하고 그 안에서 명령을 실행하라는 지시입니다.
- `--rm`은 명령이 완료된 후 컨테이너 인스턴스를 종료하라고 시스템에 지시합니다.

완전히 조립된 컨테이너 실행 명령은 다음과 같습니다:

```bash
docker run --rm -it 'community.wave.seqera.io/library/cowpy:1.1.5--3db457ae1977a273'
```

이 명령을 실행하면 프롬프트가 `(base) root@b645838b3314:/tmp#`과 같은 것으로 변경되어 이제 컨테이너 내부에 있음을 나타냅니다.

`ls`를 실행하여 디렉토리 내용을 나열하여 이를 확인할 수 있습니다:

```bash
ls /
```

??? success "명령 출력"

    ```console
    bin  boot  dev  etc  home  lib  lib64  media  mnt  opt  proc  root  run  sbin  srv  sys  tmp  usr  var
    ```

컨테이너 내부의 파일 시스템이 호스트 시스템의 파일 시스템과 다르다는 것을 볼 수 있습니다.

!!! Tip

    컨테이너를 실행할 때 기본적으로 호스트 시스템과 격리됩니다.
    이는 다음 구문을 사용하여 `docker run` 명령의 일부로 볼륨을 마운트하도록 명시적으로 허용하지 않는 한 컨테이너가 호스트 시스템의 파일에 액세스할 수 없음을 의미합니다:

    ```bash title="구문"
    -v <outside_path>:<inside_path>
    ```

    이것은 파일 시스템의 해당 부분에 액세스하는 데 사용할 수 있는 컨테이너 벽을 통해 터널을 효과적으로 설정합니다.

    이에 대한 자세한 내용은 [Hello Nextflow의 Part 5](../hello_nextflow/05_hello_containers.md)에서 다룹니다.

#### 4.1.3. `cowpy` 도구 실행

컨테이너 내부에서 `cowpy` 명령을 직접 실행할 수 있습니다.

```bash
cowpy "Hello Containers"
```

??? success "명령 출력"

    ```console
    ______________________________________________________
    < Hello Containers >
    ------------------------------------------------------
        \   ^__^
          \  (oo)\_______
            (__)\       )\/\
              ||----w |
              ||     ||
    ```

이것은 지정한 텍스트가 포함된 말풍선이 있는 기본 소 캐릭터(또는 'cowacter')의 ASCII 아트를 생성합니다.

이제 기본 사용법을 테스트했으니 몇 가지 매개변수를 제공해 볼 수 있습니다.
예를 들어 도구 문서에 따르면 `-c`로 캐릭터를 설정할 수 있습니다.

```bash
cowpy "Hello Containers" -c tux
```

??? success "명령 출력"

    ```console
    __________________
    < Hello Containers >
    ------------------
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

이번에는 `-c tux` 매개변수를 지정했기 때문에 ASCII 아트 출력에 Linux 펭귄 Tux가 표시됩니다.

컨테이너 내부에 있으므로 시스템 자체에 라이브러리를 설치할 필요 없이 입력 매개변수를 다양하게 변경하면서 cowpy 명령을 원하는 만큼 실행할 수 있습니다.

??? tip "다른 사용 가능한 캐릭터"

    다른 캐릭터를 선택하려면 '-c' 플래그를 사용하세요:

    `beavis`, `cheese`, `daemon`, `dragonandcow`, `ghostbusters`, `kitty`, `moose`, `milk`, `stegosaurus`, `turkey`, `turtle`, `tux`

마음껏 가지고 놀아 보세요.
완료되면 `exit` 명령을 사용하여 컨테이너를 종료하세요:

```bash
exit
```

일반 셸로 돌아옵니다.

### 4.2. Workflow에서 컨테이너 사용

pipeline을 실행할 때 각 단계에서 사용할 컨테이너를 Nextflow에 알려주고, 중요하게도 방금 수행한 모든 작업(컨테이너 가져오기, 시작, 명령 실행 및 완료 시 컨테이너 해제)을 처리하도록 하고 싶습니다.

좋은 소식: 바로 Nextflow가 우리를 위해 할 일입니다.
각 process에 대해 컨테이너를 지정하기만 하면 됩니다.

이것이 어떻게 작동하는지 보여주기 위해 세 번째 단계에서 생성된 수집된 인사말 파일에서 `cowpy`를 실행하는 workflow의 또 다른 버전을 만들었습니다.

<figure class="excalidraw">
--8<-- "docs/en/docs/nextflow_run/img/hello-pipeline-cowpy.svg"
</figure>

이것은 말풍선에 세 개의 인사말이 포함된 ASCII 아트가 포함된 파일을 출력해야 합니다.

#### 4.2.1. 코드 검토

workflow는 `cowpy`를 실행하는 추가 단계를 제외하고 이전 것과 매우 유사합니다.

??? full-code "전체 코드 파일"

    ```groovy title="2d-container.nf" linenums="1" hl_lines="7 15 32 39 59-62"
    #!/usr/bin/env nextflow

    // 모듈 포함
    include { sayHello } from './modules/sayHello.nf'
    include { convertToUpper } from './modules/convertToUpper.nf'
    include { collectGreetings } from './modules/collectGreetings.nf'
    include { cowpy } from './modules/cowpy.nf'

    /*
    * Pipeline parameters
    */
    params {
        input: Path
        batch: String = 'batch'
        character: String
    }

    workflow {

        main:
        // 입력용 채널 생성
        greeting_ch = channel.fromPath(params.input)
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

        publish:
        first_output = sayHello.out
        uppercased = convertToUpper.out
        collected = collectGreetings.out.outfile
        batch_report = collectGreetings.out.report
        cowpy_art = cowpy.out
    }

    output {
        first_output {
            path '2d-container/intermediates'
            mode 'copy'
        }
        uppercased {
            path '2d-container/intermediates'
            mode 'copy'
        }
        collected {
            path '2d-container/intermediates'
            mode 'copy'
        }
        batch_report {
            path '2d-container'
            mode 'copy'
        }
        cowpy_art {
            path '2d-container'
            mode 'copy'
        }
    }
    ```

이 workflow가 모듈 파일에서 `cowpy` process를 가져오고 `collectGreetings()` 호출의 출력과 `params.character`라는 입력 매개변수에서 호출하는 것을 볼 수 있습니다.

```groovy title="2d-container.nf" linenums="25"
// cowpy로 ASCII 아트 생성
cowpy(collectGreetings.out, params.character)
```

ASCII 아트를 생성하기 위해 cowpy 명령을 적용하는 `cowpy` process는 `cowpy.nf` 모듈에 정의되어 있습니다.

??? full-code "전체 코드 파일"

    ```groovy title="modules/cowpy.nf" linenums="1"
    #!/usr/bin/env nextflow

    // Generate ASCII art with cowpy (https://github.com/jeffbuttars/cowpy)
    process cowpy {

        container 'community.wave.seqera.io/library/cowpy:1.1.5--3db457ae1977a273'

        input:
        path input_file
        val character

        output:
        path "cowpy-${input_file}"

        script:
        """
        cat ${input_file} | cowpy -c "${character}" > cowpy-${input_file}
        """
    }
    ```

`cowpy` process에는 두 개의 입력이 필요합니다: 말풍선에 넣을 텍스트가 포함된 입력 파일의 경로(`input_file`)와 캐릭터 변수의 값입니다.

중요하게도 이전에 사용한 컨테이너 URI를 가리키는 `container 'community.wave.seqera.io/library/cowpy:1.1.5--3db457ae1977a273'` 줄도 포함합니다.

#### 4.2.2. 구성에서 Docker가 활성화되어 있는지 확인

`nextflow.config` 구성 파일을 소개하여 이 교육 과정의 Part 3를 약간 앞당기겠습니다. 이것은 Nextflow가 workflow 실행을 구성하는 주요 방법 중 하나입니다.
현재 디렉토리에 `nextflow.config`라는 이름의 파일이 있으면 Nextflow가 자동으로 로드하고 포함된 모든 구성을 적용합니다.

이를 위해 Docker를 활성화하는 한 줄의 코드가 포함된 `nextflow.config` 파일을 포함했습니다.

```groovy title="nextflow.config" linenums="1"
docker.enabled = true
```

이 구성은 호환되는 컨테이너를 지정하는 모든 process에 대해 Docker를 사용하도록 Nextflow에 지시합니다.

!!! tip

    기술적으로 `-with-docker <container>` 매개변수를 사용하여 실행별로 명령줄에서 Docker 실행을 활성화할 수 있습니다.
    그러나 이것은 전체 workflow에 대해 하나의 컨테이너만 지정할 수 있는 반면, 방금 보여드린 접근 방식은 process별로 다른 컨테이너를 지정할 수 있습니다.
    후자는 모듈성, 코드 유지 관리 및 재현성에 훨씬 더 좋습니다.

#### 4.2.3. Workflow 실행

요약하자면, 우리가 실행하려는 것은 다음과 같습니다:

<figure class="excalidraw">
--8<-- "docs/en/docs/hello_nextflow/img/hello_pipeline_complete.svg"
</figure>

작동할 것 같나요?

`-resume` 플래그와 함께 workflow를 실행하고 캐릭터를 turkey로 지정해 보겠습니다.

```bash
nextflow run 2d-container.nf --input data/greetings.csv --character turkey -resume
```

??? success "명령 출력"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `2d-container.nf` [elegant_brattain] DSL2 - revision: 028a841db1

    executor >  local (1)
    [95/fa0bac] sayHello (3)       | 3 of 3, cached: 3 ✔
    [92/32533f] convertToUpper (3) | 3 of 3, cached: 3 ✔
    [aa/e697a2] collectGreetings   | 1 of 1, cached: 1 ✔
    [7f/caf718] cowpy              | 1 of 1 ✔
    ```

처음 세 단계는 이전에 이미 실행했으므로 캐시되었지만 `cowpy` process는 새로운 것이므로 실제로 실행됩니다.

`results` 디렉토리에서 `cowpy` 단계의 출력을 찾을 수 있습니다.

??? abstract "파일 내용"

    ```console title="results/2d-container/cowpy-COLLECTED-batch-output.txt"
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

캐릭터가 수집된 대문자 인사말 파일에서 실행되었기 때문에 모든 인사말을 말하고 있는 것을 볼 수 있습니다.

더 중요한 것은 cowpy와 모든 의존성을 제대로 설치하지 않고도 pipeline의 일부로 이것을 실행할 수 있었다는 것입니다.
그리고 이제 pipeline을 협력자와 공유하고 위에서 언급한 Docker 또는 그 대안(예: Singularity/Apptainer) 외에는 아무것도 설치할 필요 없이 인프라에서 실행하도록 할 수 있습니다.

#### 4.2.4. Nextflow가 컨테이너화된 작업을 어떻게 시작했는지 검사

이 섹션의 마지막 코다로, `cowpy` process 호출 중 하나의 작업 하위 디렉토리를 살펴보고 Nextflow가 컨테이너와 함께 작동하는 방식에 대해 더 많은 통찰력을 얻어 보겠습니다.

`nextflow run` 명령의 출력을 확인하여 `cowpy` process의 작업 하위 디렉토리 경로를 찾으세요.
위에 표시된 실행에서 얻은 것을 보면 `cowpy` process의 콘솔 로그 줄은 `[7f/caf718]`로 시작합니다.
이것은 다음 축약된 디렉토리 경로에 해당합니다: `work/7f/caf718`.

그 디렉토리에서 Nextflow가 pipeline을 실행하는 과정에서 사용자를 대신하여 실행한 모든 명령이 포함된 `.command.run` 파일을 찾을 수 있습니다.

??? abstract "파일 내용"

    ```console title="work/7f/caf71890cce1667c094d880f4b6dcc/.command.run"
    #!/bin/bash
    ### ---
    ### name: 'cowpy'
    ### container: 'community.wave.seqera.io/library/cowpy:1.1.5--3db457ae1977a273'
    ### outputs:
    ### - 'cowpy-COLLECTED-batch-output.txt'
    ### ...
    set -e
    set -u
    NXF_DEBUG=${NXF_DEBUG:=0}; [[ $NXF_DEBUG > 1 ]] && set -x
    NXF_ENTRY=${1:-nxf_main}

    ...

    nxf_launch() {
        docker run -i --cpu-shares 1024 -e "NXF_TASK_WORKDIR" -v /workspaces/training/nextflow-run/work:/workspaces/training/nextflow-run/work -w "$NXF_TASK_WORKDIR" --name $NXF_BOXID community.wave.seqera.io/library/cowpy:1.1.5--3db457ae1977a273 /bin/bash -ue /workspaces/training/nextflow-run/work/7f/caf71890cce1667c094d880f4b6dcc/.command.sh
    }
    ```

이 파일에서 `nxf_launch`를 검색하면 다음과 같은 것을 볼 수 있습니다:

```console
nxf_launch() {
    docker run -i --cpu-shares 1024 -e "NXF_TASK_WORKDIR" -v /workspaces/training/nextflow-run/work:/workspaces/training/nextflow-run/work -w "$NXF_TASK_WORKDIR" --name $NXF_BOXID community.wave.seqera.io/library/pip_cowpy:131d6a1b707a8e65 /bin/bash -ue /workspaces/training/nextflow-run/work/7f/caf7189fca6c56ba627b75749edcb3/.command.sh
}
```

이 시작 명령은 Nextflow가 수동으로 실행할 때와 매우 유사한 `docker run` 명령을 사용하여 process 호출을 시작하고 있음을 보여줍니다.
또한 해당 작업 하위 디렉토리를 컨테이너에 마운트하고, 그에 따라 컨테이너 내부의 작업 디렉토리를 설정하고, `.command.sh` 파일의 템플릿화된 bash 스크립트를 실행합니다.

이것은 이전 섹션에서 수동으로 수행해야 했던 모든 힘든 작업이 이제 Nextflow에 의해 수행된다는 것을 확인합니다!

### 요약

컨테이너가 소프트웨어 도구 버전을 관리하고 재현성을 보장하는 데 어떤 역할을 하는지 이해합니다.

더 일반적으로 실제 Nextflow pipeline의 핵심 구성 요소가 무엇이고 어떻게 구성되어 있는지에 대한 기본적인 이해가 있습니다.
Nextflow가 여러 입력을 효율적으로 처리하고, 함께 연결된 여러 단계로 구성된 workflow를 실행하고, 모듈식 코드 구성 요소를 활용하고, 더 큰 재현성과 이식성을 위해 컨테이너를 활용하는 방법의 기본 사항을 알고 있습니다.

### 다음 단계

또 다른 휴식을 취하세요! Nextflow pipeline이 작동하는 방법에 대한 많은 정보였습니다.

교육의 마지막 섹션에서는 구성 주제에 대해 더 깊이 파고들 것입니다.
인프라에 맞게 pipeline 실행을 구성하고 입력 및 매개변수의 구성을 관리하는 방법을 배웁니다.

---

## 퀴즈

<quiz>
Nextflow가 각 process 호출에 대해 별도의 작업 디렉토리를 생성하는 이유는 무엇인가요?
- [ ] 실행 속도를 향상시키기 위해
- [ ] 메모리 사용량을 줄이기 위해
- [x] 실행을 격리하고 출력 간의 충돌을 피하기 위해
- [ ] 병렬 파일 압축을 활성화하기 위해

자세히 알아보기: [1.3. 원본 출력 및 로그 찾기](#13-원본-출력-및-로그-찾기)
</quiz>

<quiz>
workflow를 실행할 때 `-ansi-log false` 옵션은 무엇을 하나요?
- [ ] 모든 콘솔 출력을 비활성화합니다
- [x] 출력에서 색상을 제거합니다
- [x] 한 줄에 압축하는 대신 모든 작업 디렉토리 경로를 표시합니다
- [ ] 상세 디버깅 모드를 활성화합니다

자세히 알아보기: [1.3.2. 터미널에 더 많은 세부 정보 표시](#132-터미널에-더-많은-세부-정보-표시)

이 스타일을 선호하는 경우 다음 환경 변수 중 하나를 사용할 수도 있습니다:

```bash
export NXF_ANSI_LOG=0
# 또는
export NO_COLOR=1
```

</quiz>

<quiz>
`#!groovy channel.fromPath(params.input).splitCsv().map { line -> line[0] }` 코드에서 `#!groovy .map { line -> line[0] }`는 무엇을 하나요?
- [ ] 빈 줄을 필터링합니다
- [ ] 줄을 알파벳순으로 정렬합니다
- [x] 각 CSV 행에서 첫 번째 열을 추출합니다
- [ ] 줄 수를 셉니다

자세히 알아보기: [1.4.1. CSV에서 입력 데이터 로드](#141-csv에서-입력-데이터-로드)
</quiz>

<quiz>
출력 파일 이름에 입력 값을 포함하는 것이 중요한 이유는 무엇인가요(예: `#!groovy "${greeting}-output.txt"`)?
- [ ] 처리 속도를 향상시키기 위해
- [ ] resume 기능을 활성화하기 위해
- [x] 여러 입력을 처리할 때 출력 파일이 서로 덮어쓰지 않도록 방지하기 위해
- [ ] 파일을 더 쉽게 압축하기 위해

자세히 알아보기: [1.4.3. 출력 이름 지정 방법](#143-출력-이름-지정-방법)
</quiz>

<quiz>
모듈화된 workflow에서 `include` 문의 목적은 무엇인가요?
- [ ] process 코드를 workflow 파일에 복사합니다
- [x] 외부 모듈 파일에서 process 정의를 가져옵니다
- [ ] 구성 설정을 포함합니다
- [ ] 문서 주석을 추가합니다

자세히 알아보기: [3. 모듈화된 pipeline 실행](#3-모듈화된-pipeline-실행)
</quiz>

<quiz>
workflow를 모듈화하고 `-resume`으로 실행하면 어떻게 되나요?
- [ ] 모듈식 process에 대해 캐싱이 비활성화됩니다
- [ ] 모든 작업을 다시 실행해야 합니다
- [x] 생성된 작업 스크립트를 기반으로 캐싱이 정상적으로 작동합니다
- [ ] 기본 workflow 파일만 캐시됩니다

자세히 알아보기: [3.2. Workflow 실행](#32-workflow-실행)
</quiz>

<quiz>
process 정의의 `container` 지시문은 무엇을 지정하나요?
- [ ] process의 작업 디렉토리
- [ ] 최대 메모리 할당
- [x] process 실행에 사용할 컨테이너 이미지 URI
- [ ] 출력 파일 형식

자세히 알아보기: [4.2. Workflow에서 컨테이너 사용](#42-workflow에서-컨테이너-사용)
</quiz>

<quiz>
`.command.run` 파일에서 `nxf_launch` 함수에는 무엇이 포함되어 있나요?
- [ ] Nextflow 버전 정보
- [ ] workflow 매개변수
- [x] 볼륨 마운트 및 컨테이너 설정이 포함된 `docker run` 명령
- [ ] process 입력 선언

자세히 알아보기: [4.2.4. Nextflow가 컨테이너화된 작업을 어떻게 시작했는지 검사](#424-nextflow가-컨테이너화된-작업을-어떻게-시작했는지-검사)
</quiz>

<quiz>
Nextflow가 컨테이너화된 process를 실행할 때 자동으로 처리하는 것은 무엇인가요? (해당하는 모든 것을 선택하세요)
- [x] 필요한 경우 컨테이너 이미지 가져오기
- [x] 작업 디렉토리를 컨테이너에 마운트
- [x] 컨테이너 내에서 process 스크립트 실행
- [x] 실행 후 컨테이너 인스턴스 정리

자세히 알아보기: [4. 컨테이너화된 소프트웨어 사용](#4-컨테이너화된-소프트웨어-사용)
</quiz>
