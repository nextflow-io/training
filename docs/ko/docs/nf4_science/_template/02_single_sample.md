# 파트 2: 단일 샘플 처리

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } AI 지원 번역 - [자세히 알아보기 및 개선 제안](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

파트 1에서는 {TOOL_A}와 {TOOL_B} 명령을 각각의 컨테이너에서 수동으로 테스트했습니다.
이제 동일한 명령을 Nextflow 워크플로우로 적용하겠습니다.

## 과제

이 과정의 이 파트에서는 다음을 수행하는 워크플로우를 개발합니다:

1. {PROCESS_1_DESCRIPTION}
2. {PROCESS_2_DESCRIPTION}

<figure class="excalidraw">
--8<-- "docs/en/docs/nf4_science/{DOMAIN_DIR}/img/{DIAGRAM_1}.svg"
</figure>

이는 파트 1에서 컨테이너에서 이러한 명령을 수동으로 실행했던 단계를 재현합니다.

시작점으로, 워크플로우의 주요 부분을 개괄하는 워크플로우 파일 `{DOMAIN_DIR}.nf`와 모듈의 구조를 개괄하는 두 개의 모듈 파일 {TOOL_A_MODULE}.nf 및 {TOOL_B_MODULE}.nf를 제공합니다.
이 파일들은 기능적이지 않으며, 코드의 중요한 부분을 채워 넣을 수 있는 뼈대 역할만 합니다.

## 학습 계획

개발 과정을 더욱 교육적으로 만들기 위해 {N}단계로 나누었습니다:

1. **{TOOL_A_ACTION}을 실행하는 단일 단계 워크플로우 작성.**
   모듈 생성, 가져오기, 워크플로우에서 호출하는 방법을 다룹니다.
2. **{TOOL_B_ACTION}을 실행하는 두 번째 프로세스 추가.**
   프로세스 출력을 입력으로 연결하고 보조 파일을 처리하는 방법을 소개합니다.
3. **여러 샘플에서 실행되도록 워크플로우 조정.**
   병렬 실행을 다루고 관련 파일을 함께 유지하기 위한 튜플을 소개합니다.
4. **입력 파일 배치를 포함하는 샘플시트를 받아들이도록 워크플로우 수정.**
   대량으로 입력을 제공하는 일반적인 패턴을 보여줍니다.

각 단계는 워크플로우 개발의 특정 측면에 초점을 맞춥니다.

---

## 1. {TOOL_A_ACTION}을 실행하는 단일 단계 워크플로우 작성

이 첫 번째 단계는 기본 사항에 초점을 맞춥니다: {PRIMARY_INPUT_TYPE}을 로드하고 {TOOL_A_OUTPUT_DESCRIPTION}합니다.

[파트 1](01_method.md)의 `{TOOL_A_COMMAND_NAME}` 명령을 떠올려 보세요:

```bash
{TOOL_A_COMMAND_SUMMARY}
```

이 명령은 {INPUT_DESCRIPTION}을 받아 {OUTPUT_DESCRIPTION}을 생성합니다.
컨테이너 URI는 `{TOOL_A_CONTAINER_URI}`였습니다.

이 정보를 가져와서 세 단계로 Nextflow에 적용하겠습니다:

1. 입력 설정
2. 프로세스 작성 및 워크플로우에서 호출
3. 출력 처리 구성

### 1.1. 입력 설정

입력 매개변수를 선언하고, 편리한 기본값을 제공하는 테스트 프로파일을 생성하고, 입력 채널을 생성해야 합니다.

#### 1.1.1. 입력 매개변수 선언 추가

메인 워크플로우 파일 `{DOMAIN_DIR}.nf`의 `Pipeline parameters` 섹션 아래에 `{PRIMARY_PARAM_NAME}`이라는 CLI 매개변수를 선언하세요.

=== "후"

    ```groovy title="{DOMAIN_DIR}.nf" linenums="5" hl_lines="4-7"
    /*
     * Pipeline parameters
     */
    params {
        // Primary input
        {PRIMARY_PARAM_NAME}: Path
    }
    ```

=== "전"

    ```groovy title="{DOMAIN_DIR}.nf" linenums="5"
    /*
     * Pipeline parameters
     */

    // Primary input
    ```

이것은 CLI 매개변수를 설정하지만, 개발 중에 워크플로우를 실행할 때마다 파일 경로를 입력하고 싶지 않습니다.
기본값을 제공하는 여러 옵션이 있으며, 여기서는 테스트 프로파일을 사용합니다.

#### 1.1.2. `nextflow.config`에 기본값이 있는 테스트 프로파일 생성

테스트 프로파일은 명령줄에서 입력을 지정하지 않고 워크플로우를 시험해 볼 수 있는 편리한 기본값을 제공합니다.
이는 Nextflow 생태계에서 일반적인 관례입니다(자세한 내용은 [Hello Config](../../hello_nextflow/06_hello_config.md) 참조).

`nextflow.config`에 `profiles` 블록을 추가하고, `{PRIMARY_PARAM_NAME}` 매개변수를 테스트 {PRIMARY_INPUT_TYPE} 중 하나로 설정하는 `test` 프로파일을 추가하세요.

=== "후"

    ```groovy title="nextflow.config" linenums="1" hl_lines="3-7"
    docker.enabled = true

    profiles {
        test {
            params.{PRIMARY_PARAM_NAME} = "${projectDir}/data/{TEST_INPUT_PATH}"
        }
    }
    ```

=== "전"

    ```groovy title="nextflow.config" linenums="1"
    docker.enabled = true
    ```

여기서는 워크플로우 스크립트가 위치한 디렉토리를 가리키는 내장 Nextflow 변수인 `${projectDir}`을 사용합니다.
이를 통해 절대 경로를 하드코딩하지 않고도 데이터 파일 및 기타 리소스를 쉽게 참조할 수 있습니다.

#### 1.1.3. 입력 채널 설정

{INPUT_CHANNEL_INSTRUCTIONS}

### 1.2. {TOOL_A_NAME} 모듈 작성

{MODULE_INSTRUCTIONS_COVERING: container, input, output, script}

### 1.3. 워크플로우에서 모듈 가져오기 및 호출

{IMPORT_AND_CALL_INSTRUCTIONS}

### 1.4. 워크플로우 실행

이 시점에서 완전히 기능하는 단일 단계 워크플로우가 있습니다.

테스트 프로파일에 설정된 기본값을 사용하고 명령줄에 경로를 작성하지 않으려면 `-profile test`로 실행할 수 있습니다.

```bash
nextflow run {DOMAIN_DIR}.nf -profile test
```

??? success "명령 출력"

    ```console
    N E X T F L O W   ~  version 25.10.2

    ┃ Launching `{DOMAIN_DIR}.nf` [{RUN_NAME}] DSL2 - revision: {HASH}

    executor >  local (1)
    [{HASH}] {PROCESS_A_NAME} (1) | 1 of 1 ✔
    ```

{VERIFY_OUTPUT_INSTRUCTIONS}

### 핵심 정리

프로세스를 포함하는 모듈을 생성하고, 워크플로우로 가져오고, 입력 채널로 호출하고, 결과를 게시하는 방법을 알게 되었습니다.

### 다음 단계

추가 분석 단계를 연결하기 위해 두 번째 프로세스를 추가합니다.

---

## 2. {TOOL_B_ACTION}을 실행하는 두 번째 프로세스 추가

{DESCRIBE_WHAT_THIS_STEP_ADDS}

[파트 1](01_method.md)의 `{TOOL_B_COMMAND_NAME}` 명령을 떠올려 보세요:

```bash
{TOOL_B_COMMAND_SUMMARY}
```

### 2.1. {TOOL_B_NAME} 모듈 작성

{MODULE_INSTRUCTIONS}

### 2.2. 워크플로우에서 모듈 가져오기 및 호출

{IMPORT_AND_CALL_INSTRUCTIONS_CHAINING_OUTPUTS}

### 2.3. 워크플로우 실행

```bash
nextflow run {DOMAIN_DIR}.nf -profile test
```

??? success "명령 출력"

    ```console
    {EXPECTED_OUTPUT}
    ```

{VERIFY_OUTPUT_INSTRUCTIONS}

### 핵심 정리

프로세스 출력을 입력으로 연결하고 워크플로우에서 보조 파일을 처리하는 방법을 알게 되었습니다.

### 다음 단계

여러 샘플을 병렬로 처리하도록 워크플로우를 확장합니다.

---

## 3. 여러 샘플에서 실행되도록 워크플로우 조정

지금까지 워크플로우는 단일 샘플을 처리합니다.
여러 샘플을 처리하려면 입력 제공 방식을 수정하고 Nextflow의 데이터플로우 패러다임을 활용하여 실행을 병렬화해야 합니다.

### 3.1. {BATCH_INPUT_INSTRUCTIONS}

{INSTRUCTIONS_FOR_SWITCHING_TO_SAMPLESHEET_OR_BATCH_INPUT}

### 3.2. 여러 샘플에서 워크플로우 실행

```bash
nextflow run {DOMAIN_DIR}.nf -profile test
```

??? success "명령 출력"

    ```console
    {EXPECTED_OUTPUT_MULTIPLE_SAMPLES}
    ```

### 핵심 정리

Nextflow의 데이터플로우 패러다임을 활용하여 여러 입력 샘플에 걸쳐 샘플별 처리를 병렬화하는 방법을 알게 되었습니다.

### 다음 단계

[파트 3](03_multi_sample.md)에서는 모든 샘플의 결과를 결합하기 위해 다중 샘플 집계를 추가합니다.
