# 파트 1: 방법 개요 및 수동 테스트

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } AI 지원 번역 - [자세히 알아보기 및 개선 제안](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

{BRIEF_METHOD_DESCRIPTION}

![파이프라인 개요](img/{PIPELINE_DIAGRAM}.png)

{OPTIONAL_EXTENDED_METHOD_DESCRIPTION}

### 방법

{DESCRIBE_THE_TWO_APPROACHES: single-sample and multi-sample/aggregation}

워크플로우 코드를 작성하기 전에, 먼저 테스트 데이터에서 명령을 수동으로 실행해 보겠습니다.

### 데이터셋

다음 데이터 및 관련 리소스를 제공합니다:

- **{PRIMARY_INPUT_DESCRIPTION}**
- **{SAMPLE_DESCRIPTION}**
- **{ADDITIONAL_RESOURCES_DESCRIPTION}**

### 소프트웨어

주요 도구는 [{TOOL_A}]({TOOL_A_URL})와 [{TOOL_B}]({TOOL_B_URL})입니다.

이러한 도구는 GitHub Codespaces 환경에 설치되어 있지 않으므로, 컨테이너를 통해 사용하겠습니다([Hello Containers](../../hello_nextflow/05_hello_containers.md) 참조).

!!! note "참고"

    `nf4-science/{DOMAIN_DIR}` 디렉토리에 있는지 확인하세요. `pwd` 명령을 입력했을 때 경로의 마지막 부분이 `{DOMAIN_DIR}`이어야 합니다.

---

## 1. {SINGLE_SAMPLE_SECTION_TITLE}

{BRIEF_DESCRIPTION_OF_SINGLE_SAMPLE_APPROACH}

이 섹션에서는 단일 샘플 처리 방식을 구성하는 명령을 테스트합니다.
이 명령들은 이 과정의 파트 2에서 Nextflow 워크플로우로 작성할 것입니다.

1. {STEP_1_SUMMARY}
2. {STEP_2_SUMMARY}

<figure class="excalidraw">
--8<-- "docs/en/docs/nf4_science/{DOMAIN_DIR}/img/{DIAGRAM_1}.svg"
</figure>

먼저 하나의 샘플에서 명령을 테스트합니다.

### 1.1. {FIRST_TOOL_STEP_TITLE}

{BRIEF_DESCRIPTION_OF_STEP}

#### 1.1.1. 컨테이너 가져오기

`docker pull` 명령을 실행하여 컨테이너 이미지를 다운로드합니다:

```bash
docker pull {TOOL_A_CONTAINER_URI}
```

??? success "명령 출력"

    ```console
    {EXPECTED_PULL_OUTPUT}
    ```

#### 1.1.2. 컨테이너를 대화형으로 실행하기

컨테이너를 실행하고 `data` 디렉토리를 마운트하여 도구가 입력 파일에 접근할 수 있도록 합니다:

```bash
docker run -it -v ./data:/data {TOOL_A_CONTAINER_URI}
```

프롬프트가 변경되어 컨테이너 내부에 있음을 나타냅니다.

#### 1.1.3. 명령 실행하기

```bash
{TOOL_A_COMMAND}
```

??? success "명령 출력"

    ```console
    {EXPECTED_COMMAND_OUTPUT}
    ```

{DESCRIBE_EXPECTED_OUTPUT_FILES}

#### 1.1.4. 컨테이너 종료하기

컨테이너를 종료하려면 `exit`를 입력합니다.

```bash
exit
```

프롬프트가 정상으로 돌아와야 합니다.

### 1.2. {SECOND_TOOL_STEP_TITLE}

{BRIEF_DESCRIPTION_OF_STEP}

#### 1.2.1. 컨테이너 가져오기

```bash
docker pull {TOOL_B_CONTAINER_URI}
```

??? success "명령 출력"

    ```console
    {EXPECTED_PULL_OUTPUT}
    ```

#### 1.2.2. 컨테이너를 대화형으로 실행하기

```bash
docker run -it -v ./data:/data {TOOL_B_CONTAINER_URI}
```

#### 1.2.3. 명령 실행하기

```bash
{TOOL_B_COMMAND}
```

??? success "명령 출력"

    ```console
    {EXPECTED_COMMAND_OUTPUT}
    ```

{DESCRIBE_EXPECTED_OUTPUT_FILES}

#### 1.2.4. 컨테이너 종료하기

컨테이너를 종료하려면 `exit`를 입력합니다.

```bash
exit
```

프롬프트가 정상으로 돌아와야 합니다.
단일 샘플 처리 테스트가 완료되었습니다.

---

## 2. {MULTI_SAMPLE_SECTION_TITLE}

{BRIEF_DESCRIPTION_OF_MULTI_SAMPLE_APPROACH}

{EXPLAIN_WHY_MULTI_SAMPLE_IS_NEEDED}

이 섹션에서는 다중 샘플 처리에 필요한 추가 명령을 테스트합니다.
이 명령들은 이 과정의 파트 3에서 Nextflow 워크플로우로 작성할 것입니다.

1. {AGGREGATION_STEP_1_SUMMARY}
2. {AGGREGATION_STEP_2_SUMMARY}

<figure class="excalidraw">
--8<-- "docs/en/docs/nf4_science/{DOMAIN_DIR}/img/{DIAGRAM_2}.svg"
</figure>

### 2.1. {AGGREGATION_STEP_TITLE}

{STEP_INSTRUCTIONS_FOLLOWING_SAME_PATTERN_AS_ABOVE}

---

### 핵심 정리

각 컨테이너에서 {TOOL_A}와 {TOOL_B} 명령을 테스트하는 방법을 알게 되었습니다. 여기에는 {MULTI_SAMPLE_SUMMARY} 방법도 포함됩니다.

### 다음 단계

동일한 명령을 컨테이너를 사용하여 작업을 실행하는 워크플로우로 작성하는 방법을 학습합니다.
