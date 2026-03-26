# 파트 1: 플러그인 기초

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } AI 지원 번역 - [자세히 알아보기 및 개선 제안](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

이 섹션에서는 플러그인이 Nextflow를 어떻게 확장하는지 학습한 후, 세 가지 플러그인을 직접 사용해 봅니다.

---

## 1. 플러그인 작동 방식

플러그인은 여러 유형의 확장을 통해 Nextflow를 확장합니다:

| 확장 유형         | 기능                                             | 예시                         |
| ----------------- | ------------------------------------------------ | ---------------------------- |
| Functions         | 워크플로우에서 호출 가능한 사용자 정의 함수 추가 | `samplesheetToList()`        |
| Workflow monitors | 작업 완료 등의 이벤트에 응답                     | 사용자 정의 로깅, Slack 알림 |
| Executors         | 작업 실행 백엔드 추가                            | AWS Batch, Kubernetes        |
| Filesystems       | 스토리지 백엔드 추가                             | S3, Azure Blob               |

Functions와 workflow monitors("trace observers"라고도 함)는 플러그인 개발자들이 가장 많이 사용하는 유형입니다.
Executors와 filesystems는 일반적으로 플랫폼 벤더가 개발합니다.

다음 실습에서는 function 플러그인과 observer 플러그인을 사용해 보며 두 유형을 직접 확인합니다.

---

## 2. Function 플러그인 사용

Function 플러그인은 워크플로우에 import하여 호출할 수 있는 함수를 제공합니다.
nf-hello(간단한 예제)와 nf-schema(실제 환경에서 널리 사용되는 플러그인), 두 가지를 사용해 봅니다.
두 실습 모두 동일한 `hello.nf` 파이프라인을 수정하므로, 플러그인이 기존 워크플로우를 어떻게 개선하는지 확인할 수 있습니다.

### 2.1. nf-hello: 직접 작성한 코드 대체

[nf-hello](https://github.com/nextflow-io/nf-hello) 플러그인은 무작위 문자열을 생성하는 `randomString` 함수를 제공합니다.
파이프라인에는 이미 이 함수의 인라인 버전이 정의되어 있으며, 이를 플러그인의 함수로 교체합니다.

#### 2.1.1. 시작 상태 확인

파이프라인을 확인합니다:

```bash
cat hello.nf
```

```groovy title="Output"
#!/usr/bin/env nextflow

params.input = 'greetings.csv'

/**
 * 무작위 영숫자 문자열 생성
 */
def randomString(int length) {
    def chars = ('a'..'z') + ('A'..'Z') + ('0'..'9')
    def random = new Random()
    return (1..length).collect { chars[random.nextInt(chars.size())] }.join()
}

process SAY_HELLO {
    input:
        val greeting
    output:
        stdout
    script:
    """
    echo '$greeting'
    """
}

workflow {
    greeting_ch = channel.fromPath(params.input)
        .splitCsv(header: true)
        .map { row -> "${row.greeting}_${randomString(8)}" }
    SAY_HELLO(greeting_ch)
    SAY_HELLO.out.view { result -> "Output: ${result.trim()}" }
}
```

파이프라인은 `randomString` 함수를 인라인으로 정의하고, 이를 사용하여 각 인사말에 무작위 ID를 추가합니다.

실행합니다:

```bash
nextflow run hello.nf
```

```console title="Output"
Output: Hello_aBcDeFgH
Output: Bonjour_xYzWvUtS
Output: Holà_qRsPdMnK
Output: Ciao_jLhGfEcB
Output: Hallo_tNwOiAuR
```

출력 순서와 무작위 문자열은 다를 수 있으며, 스크립트를 다시 실행하면 다른 무작위 인사말이 생성됩니다.

#### 2.1.2. 플러그인 설정

인라인 함수를 플러그인의 함수로 교체합니다. `nextflow.config`에 다음을 추가합니다:

```groovy title="nextflow.config"
// 플러그인 개발 실습을 위한 설정
plugins {
    id 'nf-hello@0.5.0'
}
```

플러그인은 `nextflow.config`의 `plugins {}` 블록에 선언합니다.
Nextflow는 커뮤니티 및 공식 플러그인의 중앙 저장소인 [Nextflow Plugin Registry](https://registry.nextflow.io/)에서 플러그인을 자동으로 다운로드합니다.

#### 2.1.3. 플러그인 함수 사용

인라인 `randomString` 함수를 플러그인 버전으로 교체합니다:

=== "후"

    ```groovy title="hello.nf" hl_lines="3"
    #!/usr/bin/env nextflow

    include { randomString } from 'plugin/nf-hello'

    params.input = 'greetings.csv'

    process SAY_HELLO {
        input:
            val greeting
        output:
            stdout
        script:
        """
        echo '$greeting'
        """
    }

    workflow {
        greeting_ch = channel.fromPath(params.input)
            .splitCsv(header: true)
            .map { row -> "${row.greeting}_${randomString(8)}" }
        SAY_HELLO(greeting_ch)
        SAY_HELLO.out.view { result -> "Output: ${result.trim()}" }
    }
    ```

=== "전"

    ```groovy title="hello.nf" hl_lines="5-12"
    #!/usr/bin/env nextflow

    params.input = 'greetings.csv'

    /**
     * 무작위 영숫자 문자열 생성
     */
    def randomString(int length) {
        def chars = ('a'..'z') + ('A'..'Z') + ('0'..'9')
        def random = new Random()
        return (1..length).collect { chars[random.nextInt(chars.size())] }.join()
    }

    process SAY_HELLO {
        input:
            val greeting
        output:
            stdout
        script:
        """
        echo '$greeting'
        """
    }

    workflow {
        greeting_ch = channel.fromPath(params.input)
            .splitCsv(header: true)
            .map { row -> "${row.greeting}_${randomString(8)}" }
        SAY_HELLO(greeting_ch)
        SAY_HELLO.out.view { result -> "Output: ${result.trim()}" }
    }
    ```

`include` 구문은 검증되고 테스트되었으며, 더 많은 기여자들이 버그를 발견하고 수정하는 라이브러리에서 `randomString`을 가져옵니다.
각 파이프라인이 함수의 자체 복사본을 유지하는 대신, 플러그인을 사용하는 모든 파이프라인이 동일하게 검증된 구현을 사용합니다.
이를 통해 코드 중복과 그에 따른 유지보수 부담이 줄어듭니다.
`#!groovy include { function } from 'plugin/plugin-id'` 구문은 Nextflow 모듈에서 사용하는 `include`와 동일하며, `plugin/` 접두사가 붙습니다.
nf-hello 저장소의 GitHub에서 [`randomString`의 소스 코드](https://github.com/nextflow-io/nf-hello/blob/e67bddebfa589c7ae51f41bf780c92068dc09e93/plugins/nf-hello/src/main/nextflow/hello/HelloExtension.groovy#L110)를 확인할 수 있습니다.

#### 2.1.4. 실행

```bash
nextflow run hello.nf
```

```console title="Output"
Pipeline is starting! 🚀
Output: Hello_yqvtclcc
Output: Bonjour_vwwpyzcs
Output: Holà_wrghmgab
Output: Ciao_noniajuy
Output: Hallo_tvrtuxtp
Pipeline complete! 👋
```

(무작위 문자열은 다를 수 있습니다.)

출력에는 여전히 무작위 접미사가 있지만, 이제 `randomString`은 인라인 코드 대신 nf-hello 플러그인에서 제공됩니다.
"Pipeline is starting!"과 "Pipeline complete!" 메시지는 새로 추가된 것으로, 파트 5에서 살펴볼 플러그인의 observer 컴포넌트에서 출력됩니다.

Nextflow는 플러그인을 처음 사용할 때 자동으로 다운로드하므로, `nf-hello@0.5.0`을 선언한 모든 파이프라인은 프로젝트 간에 코드를 복사하지 않고도 동일하게 테스트된 `randomString` 함수를 사용할 수 있습니다.

이제 function 플러그인 사용의 세 단계를 확인했습니다: `nextflow.config`에 선언하고, `include`로 함수를 가져오고, 워크플로우에서 호출합니다.
다음 실습에서는 동일한 단계를 실제 환경의 플러그인에 적용합니다.

### 2.2. nf-schema: 검증된 CSV 파싱

[nf-schema](https://github.com/nextflow-io/nf-schema) 플러그인은 가장 널리 사용되는 Nextflow 플러그인 중 하나입니다.
예상 열과 유형을 정의하는 JSON 스키마를 사용하여 CSV/TSV 파일을 파싱하는 `samplesheetToList` 함수를 제공합니다.

파이프라인은 현재 `splitCsv`와 수동 `map`을 사용하여 `greetings.csv`를 읽고 있지만, nf-schema를 사용하면 검증된 스키마 기반 파싱으로 대체할 수 있습니다.
JSON 스키마 파일(`greetings_schema.json`)은 실습 디렉토리에 이미 제공되어 있습니다.

??? info "스키마란 무엇인가요?"

    스키마는 유효한 데이터의 형태를 공식적으로 정의한 것입니다.
    어떤 열이 필요한지, 각 값의 유형(문자열, 숫자 등)은 무엇인지, 어떤 필드가 필수인지 등을 정의합니다.

    일종의 계약으로 생각할 수 있습니다. 입력 데이터가 스키마와 일치하지 않으면, 파이프라인 후반부에서 혼란스러운 오류가 발생하기 전에 도구가 문제를 조기에 감지할 수 있습니다.

#### 2.2.1. 스키마 확인

```bash
cat greetings_schema.json
```

```json title="Output"
{
  "$schema": "https://json-schema.org/draft/2020-12/schema",
  "type": "array",
  "items": {
    "type": "object",
    "properties": {
      "greeting": {
        "type": "string",
        "description": "The greeting text"
      },
      "language": {
        "type": "string",
        "description": "The language of the greeting"
      }
    },
    "required": ["greeting"]
  }
}
```

스키마는 두 개의 열(`greeting`과 `language`)을 정의하고 `greeting`을 필수로 표시합니다.
`greeting` 열이 없는 CSV를 전달하면, nf-schema가 파이프라인 실행 전에 오류를 감지합니다.

#### 2.2.2. 설정에 nf-schema 추가

두 플러그인을 모두 포함하도록 `nextflow.config`를 업데이트합니다:

=== "후"

    ```groovy title="nextflow.config" hl_lines="3"
    plugins {
        id 'nf-hello@0.5.0'
        id 'nf-schema@2.6.1'
    }
    ```

=== "전"

    ```groovy title="nextflow.config"
    plugins {
        id 'nf-hello@0.5.0'
    }
    ```

#### 2.2.3. samplesheetToList를 사용하도록 hello.nf 업데이트

`splitCsv` 입력을 `samplesheetToList`로 교체합니다:

=== "후"

    ```groovy title="hello.nf" hl_lines="4 20 21 22"
    #!/usr/bin/env nextflow

    include { randomString } from 'plugin/nf-hello'
    include { samplesheetToList } from 'plugin/nf-schema'

    params.input = 'greetings.csv'

    process SAY_HELLO {
        input:
            val greeting
        output:
            stdout
        script:
        """
        echo '$greeting'
        """
    }

    workflow {
        def samplesheet_list = samplesheetToList(params.input, 'greetings_schema.json')
        greeting_ch = Channel.fromList(samplesheet_list)
            .map { row -> "${row[0]}_${randomString(8)}" }
        SAY_HELLO(greeting_ch)
        SAY_HELLO.out.view { result -> "Output: ${result.trim()}" }
    }
    ```

=== "전"

    ```groovy title="hello.nf" hl_lines="19 20 21"
    #!/usr/bin/env nextflow

    include { randomString } from 'plugin/nf-hello'

    params.input = 'greetings.csv'

    process SAY_HELLO {
        input:
            val greeting
        output:
            stdout
        script:
        """
        echo '$greeting'
        """
    }

    workflow {
        greeting_ch = channel.fromPath(params.input)
            .splitCsv(header: true)
            .map { row -> "${row.greeting}_${randomString(8)}" }
        SAY_HELLO(greeting_ch)
        SAY_HELLO.out.view { result -> "Output: ${result.trim()}" }
    }
    ```

직접 작성한 `splitCsv`와 `map` 파싱 코드가 `samplesheetToList`로 교체되었습니다. 이 함수는 검증되고 테스트된 함수로, 파이프라인 실행 전에 스키마에 대한 샘플시트 유효성 검사도 수행합니다.
직접 작성한 파싱 로직의 유지보수 부담을 줄이는 동시에, 입력이 예상 형식과 일치하지 않을 때 명확한 오류 메시지를 제공하여 파이프라인 사용자 경험을 개선합니다.
각 행은 열 순서대로 값의 목록이 되므로, `row[0]`은 인사말이고 `row[1]`은 언어입니다.

#### 2.2.4. 실행

```bash
nextflow run hello.nf
```

```console title="Output"
Pipeline is starting! 🚀
Output: Hello_diozjdwm
Output: Bonjour_speathmm
Output: Holà_dllxnzap
Output: Ciao_wzueddzc
Output: Hallo_hsxwrjbh
Pipeline complete! 👋
```

(무작위 문자열은 다를 수 있습니다.)

출력은 동일하지만, 이제 파이프라인 실행 전에 스키마가 CSV 구조를 검증합니다.
복잡한 샘플시트와 많은 열을 가진 실제 파이프라인에서는 이러한 유효성 검사가 수동 `splitCsv` + `map`으로는 놓칠 수 있는 오류를 방지합니다.

#### 2.2.5. 유효성 검사 동작 확인

스키마 유효성 검사가 무엇을 감지하는지 확인하기 위해 `greetings.csv`에 오류를 도입해 봅니다.

필수 열인 `greeting`의 이름을 `message`로 변경합니다:

```csv title="greetings.csv" hl_lines="1"
message,language
Hello,English
Bonjour,French
Holà,Spanish
Ciao,Italian
Hallo,German
```

파이프라인을 실행합니다:

```bash
nextflow run hello.nf
```

```console title="Output"
ERROR ~ Validation of samplesheet failed!

The following errors have been detected in greetings.csv:

-> Entry 1: Missing required field(s): greeting
-> Entry 2: Missing required field(s): greeting
-> Entry 3: Missing required field(s): greeting
-> Entry 4: Missing required field(s): greeting
-> Entry 5: Missing required field(s): greeting
```

스키마가 `greeting` 열을 요구하는데 찾을 수 없으므로 파이프라인 실행이 거부됩니다.

이제 필수 열을 복원하고 선택적 열인 `language`의 이름을 `lang`으로 변경합니다:

```csv title="greetings.csv" hl_lines="1"
greeting,lang
Hello,English
Bonjour,French
Holà,Spanish
Ciao,Italian
Hallo,German
```

```bash
nextflow run hello.nf
```

이번에는 파이프라인이 실행되지만 경고가 출력됩니다:

```console title="Output (partial)"
WARN: Found the following unidentified headers in greetings.csv:
	- lang
```

필수 열은 오류를 발생시키고, 선택적 열은 경고를 발생시킵니다.
이러한 조기 피드백은 수십 개의 열이 있는 실제 파이프라인에서 디버깅 시간을 절약해 줍니다.

#### 2.2.6. 유효성 검사 동작 설정

`lang`에 대한 경고는 유용하지만, 설정을 통해 심각도를 제어할 수 있습니다.
플러그인은 동작을 제어하는 자체 설정 스코프를 포함할 수 있습니다.
nf-schema 플러그인은 `validation` 설정 스코프를 포함하며, 여기서 설정을 수정하여 nf-schema의 동작을 변경할 수 있습니다.

인식되지 않는 헤더가 경고 대신 오류를 발생시키도록 `nextflow.config`에 `validation` 블록을 추가합니다:

=== "후"

    ```groovy title="nextflow.config" hl_lines="6-10"
    plugins {
        id 'nf-hello@0.5.0'
        id 'nf-schema@2.6.1'
    }

    validation {
        logging {
            unrecognisedHeaders = "error"
        }
    }
    ```

=== "전"

    ```groovy title="nextflow.config"
    plugins {
        id 'nf-hello@0.5.0'
        id 'nf-schema@2.6.1'
    }
    ```

`lang` 열이 그대로 있는 상태에서 파이프라인을 다시 실행합니다:

```bash
nextflow run hello.nf
```

```console title="Output (partial)"
Found the following unidentified headers in greetings.csv:
	- lang
 -- Check script 'hello.nf' at line: 20 or see '.nextflow.log' file for more details
```

이제 파이프라인은 경고 대신 오류로 실패합니다.
파이프라인 코드는 변경되지 않았으며, 설정만 변경되었습니다.

계속하기 전에 `greetings.csv`를 원래 상태로 복원하고 `validation` 블록을 제거합니다:

```csv title="greetings.csv"
greeting,language
Hello,English
Bonjour,French
Holà,Spanish
Ciao,Italian
Hallo,German
```

```groovy title="nextflow.config"
plugins {
    id 'nf-hello@0.5.0'
    id 'nf-schema@2.6.1'
}
```

nf-hello와 nf-schema는 모두 function 플러그인으로, `include`로 가져와 워크플로우 코드에서 호출하는 함수를 제공합니다.
다음 실습에서는 `include` 구문 없이도 작동하는 다른 유형의 플러그인을 살펴봅니다.

---

## 3. Observer 플러그인 사용: nf-co2footprint

모든 플러그인이 import할 함수를 제공하는 것은 아닙니다.
[nf-co2footprint](https://github.com/nextflow-io/nf-co2footprint) 플러그인은 **trace observer**를 사용하여 파이프라인의 리소스 사용량을 모니터링하고 탄소 발자국을 추정합니다.
파이프라인 코드를 변경할 필요 없이 설정에 추가하기만 하면 됩니다.

### 3.1. 설정에 nf-co2footprint 추가

`nextflow.config`를 업데이트합니다:

=== "후"

    ```groovy title="nextflow.config" hl_lines="4"
    plugins {
        id 'nf-hello@0.5.0'
        id 'nf-schema@2.6.1'
        id 'nf-co2footprint@1.2.0'
    }
    ```

=== "전"

    ```groovy title="nextflow.config"
    plugins {
        id 'nf-hello@0.5.0'
        id 'nf-schema@2.6.1'
    }
    ```

### 3.2. 파이프라인 실행

```bash
nextflow run hello.nf
```

플러그인은 실행 중에 여러 INFO 및 WARN 메시지를 출력합니다.
로컬 머신에서 실행되는 소규모 예제에서는 정상적인 동작입니다:

```console title="Output (partial)"
nf-co2footprint plugin  ~  version 1.2.0
WARN - [nf-co2footprint] Target zone null not found. Attempting to retrieve carbon intensity for fallback zone GLOBAL.
INFO - [nf-co2footprint] Using fallback carbon intensity from GLOBAL from CI table: 480.0 gCO₂eq/kWh.
WARN - [nf-co2footprint] Executor 'null' not mapped.
WARN - [nf-co2footprint] Fallback to: `machineType = null`, `pue = 1.0`. ...
...
WARN - [nf-co2footprint] No CPU model detected. Using default CPU power draw value (11.41 W).
WARN - [nf-co2footprint] 🔁 Requested memory is null for task 2. Using maximum consumed memory/`peak_rss` (0 GB) for CO₂e footprint computation.
```

zone, executor, CPU 모델, 메모리에 대한 경고는 플러그인이 로컬 교육 환경의 전체 하드웨어 정보를 감지할 수 없기 때문에 나타납니다.
프로덕션 환경(예: HPC 클러스터 또는 클라우드)에서는 이러한 값을 사용할 수 있어 추정치가 더 정확합니다.

마지막에 다음과 같은 줄을 확인합니다:

```console title="Output (partial)"
🌱 The workflow run used 126.76 uWh of electricity, resulting in the release of 60.84 ug of CO₂ equivalents into the atmosphere.
```

(수치는 다를 수 있습니다.)

### 3.3. 리포트 확인

플러그인은 작업 디렉토리에 출력 파일을 생성합니다:

```bash
ls co2footprint_*
```

```console title="Output"
co2footprint_report_<timestamp>.html
co2footprint_summary_<timestamp>.txt
co2footprint_trace_<timestamp>.txt
```

요약 파일을 확인합니다:

```bash
cat co2footprint_summary_*.txt
```

```console title="Output"
Total CO₂e footprint measures of this workflow run (including cached tasks):
  CO₂e emissions: 60.84 ug
  Energy consumption: 126.76 uWh
  CO₂e emissions (market): -

Which equals:
  - 3.48E-7 km travelled by car
  - It takes one tree 0.17s to sequester the equivalent amount of CO₂ from the atmosphere
  - 1.22E-7 % of a flight from Paris to London
```

(수치는 다를 수 있습니다.)

첫 번째 섹션은 원시 에너지 및 배출량 수치를 보여줍니다.
"Which equals" 섹션은 이 수치를 친숙한 단위로 변환하여 맥락을 제공합니다.
요약에는 플러그인의 설정 옵션 목록과 계산 방법의 기반이 되는 [Green Algorithms](https://doi.org/10.1002/advs.202100707) 연구 논문 인용도 포함되어 있습니다.

### 3.4. 플러그인 설정

3.2절의 "Target zone null" 경고는 플러그인에 위치가 설정되지 않았기 때문에 나타났습니다.
nf-co2footprint 플러그인은 지리적 위치를 설정할 수 있는 `co2footprint` 설정 스코프를 정의합니다.

`nextflow.config`에 `co2footprint` 블록을 추가합니다:

=== "후"

    ```groovy title="nextflow.config" hl_lines="7-9"
    plugins {
        id 'nf-hello@0.5.0'
        id 'nf-schema@2.6.1'
        id 'nf-co2footprint@1.2.0'
    }

    co2footprint {
        location = 'GB'
    }
    ```

=== "전"

    ```groovy title="nextflow.config"
    plugins {
        id 'nf-hello@0.5.0'
        id 'nf-schema@2.6.1'
        id 'nf-co2footprint@1.2.0'
    }
    ```

!!! tip "팁"

    원하는 경우 자신의 국가 코드를 사용하세요 (예: `'US'`, `'DE'`, `'FR'`).

파이프라인을 실행합니다:

```bash
nextflow run hello.nf
```

```console title="Output (partial)"
INFO - [nf-co2footprint] Using fallback carbon intensity from GB from CI table: 163.92 gCO₂eq/kWh.
```

zone 경고가 사라졌습니다.
플러그인이 이제 전역 기본값(480.0 gCO₂eq/kWh) 대신 GB 특정 탄소 집약도(163.92 gCO₂eq/kWh)를 사용합니다.

!!! note "참고"

    `WARN: Unrecognized config option 'co2footprint.location'` 메시지가 표시될 수도 있습니다.
    이는 외관상의 문제로 무시해도 됩니다. 플러그인은 여전히 값을 올바르게 읽습니다.

파트 6에서는 자신만의 플러그인을 위한 설정 스코프를 만들어 봅니다.

이 플러그인은 전적으로 observer 메커니즘을 통해 작동하며, 워크플로우 생명주기 이벤트에 연결하여 리소스 메트릭을 수집하고 파이프라인이 완료되면 리포트를 생성합니다.

이제 function 플러그인(`include`로 가져오기)과 observer 플러그인(설정만으로 활성화)을 모두 사용해 보았습니다.
이 두 가지가 가장 일반적인 확장 유형이지만, 섹션 1의 표에서 볼 수 있듯이 플러그인은 executor와 filesystem도 추가할 수 있습니다.

---

## 4. 플러그인 검색

[Nextflow Plugin Registry](https://registry.nextflow.io/)는 사용 가능한 플러그인을 찾을 수 있는 중앙 허브입니다.

![registry.nextflow.io의 nf-hello 플러그인 페이지](img/plugin-registry-nf-hello.png)

각 플러그인 페이지에는 설명, 사용 가능한 버전, 설치 방법, 문서 링크가 표시됩니다.

---

## 5. 플러그인 개발 준비

다음 섹션(파트 2-6)에서는 별도의 파이프라인 파일인 `greet.nf`를 사용합니다. 이 파일은 nf-schema에 의존하지만 nf-hello나 nf-co2footprint는 사용하지 않습니다.

nf-schema만 유지하도록 `nextflow.config`를 업데이트합니다:

```groovy title="nextflow.config"
// 플러그인 개발 실습을 위한 설정
plugins {
    id 'nf-schema@2.6.1'
}
```

co2footprint 출력 파일을 삭제합니다:

```bash
rm -f co2footprint_*
```

`hello.nf` 파일은 파트 1 작업 내용을 참조용으로 유지합니다. 이후에는 `greet.nf`로 작업합니다.

---

## 핵심 정리

세 가지 플러그인을 사용해 보았습니다:

- **nf-hello**: `randomString`을 제공하는 function 플러그인으로, `include`로 가져옵니다
- **nf-schema**: 스키마 검증 CSV 파싱을 위한 `samplesheetToList`를 제공하는 function 플러그인
- **nf-co2footprint**: `include` 없이 리소스 사용량을 자동으로 모니터링하는 observer 플러그인

주요 패턴:

- 플러그인은 `nextflow.config`에 `#!groovy plugins { id 'plugin-name@version' }`으로 선언합니다
- Function 플러그인은 `#!groovy include { function } from 'plugin/plugin-id'`가 필요합니다
- Observer 플러그인은 설정에 선언하면 자동으로 작동합니다
- 플러그인은 동작을 맞춤화하기 위한 설정 스코프를 정의할 수 있습니다 (예: `#!groovy validation {}`, `#!groovy co2footprint {}`)
- [Nextflow Plugin Registry](https://registry.nextflow.io/)에서 사용 가능한 플러그인 목록을 확인할 수 있습니다

---

## 다음 단계

다음 섹션에서는 직접 플러그인을 만드는 방법을 살펴봅니다.
플러그인 개발에 관심이 없다면 여기서 멈추거나 [요약](summary.md)으로 건너뛸 수 있습니다.

[파트 2로 계속 :material-arrow-right:](02_create_project.md){ .md-button .md-button--primary }
