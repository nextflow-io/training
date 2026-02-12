# 파트 3: 다중 샘플 집계

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } AI 지원 번역 - [자세히 알아보기 및 개선 제안](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

파트 2에서는 각 샘플을 독립적으로 처리하는 샘플별 처리 파이프라인을 구축했습니다.
이제 [파트 1](01_method.md)에서 다룬 다중 샘플 {AGGREGATION_METHOD}을 구현하도록 파이프라인을 확장하겠습니다.

## 과제

이 과정의 이번 파트에서는 다음을 수행하도록 워크플로우를 확장합니다:

<figure class="excalidraw">
--8<-- "docs/en/docs/nf4_science/{DOMAIN_DIR}/img/{DIAGRAM_2}.svg"
</figure>

1. {PER_SAMPLE_STEP_1}
2. {PER_SAMPLE_STEP_2}
3. {AGGREGATION_STEP_1}
4. {AGGREGATION_STEP_2}

이 파트는 파트 2에서 작성한 워크플로우를 직접 기반으로 합니다.

??? info "이 섹션을 시작하는 방법"

    이 과정의 이번 섹션은 [파트 2: 단일 샘플 처리](./02_single_sample.md)를 완료하고 작동하는 `{DOMAIN_DIR}.nf` 파이프라인이 있다고 가정합니다.

    파트 2를 완료하지 않았거나 이 파트를 새로 시작하려면, 파트 2 해결책을 시작점으로 사용할 수 있습니다.
    `nf4-science/{DOMAIN_DIR}/` 디렉토리 내에서 다음 명령을 실행하세요:

    ```bash
    cp solutions/part2/{DOMAIN_DIR}-2.nf {DOMAIN_DIR}.nf
    cp solutions/part2/nextflow.config .
    cp solutions/part2/modules/* modules/
    ```

    이렇게 하면 완전한 단일 샘플 처리 워크플로우를 얻을 수 있습니다.
    다음 명령을 실행하여 성공적으로 실행되는지 테스트할 수 있습니다:

    ```bash
    nextflow run {DOMAIN_DIR}.nf -profile test
    ```

## 학습 계획

이를 두 단계로 나누었습니다:

1. **{MODIFICATION_STEP_SUMMARY}.**
   프로세스 명령과 출력을 업데이트하는 내용을 다룹니다.
2. **{AGGREGATION_STEP_SUMMARY}.**
   `collect()` 연산자{AND_OTHER_CONCEPTS}를 소개합니다.

!!! note "참고"

     올바른 작업 디렉토리에 있는지 확인하세요:
     `cd /workspaces/training/nf4-science/{DOMAIN_DIR}`

---

## 1. {MODIFICATION_STEP_TITLE}

{DESCRIPTION_OF_MODIFICATION_TO_EXISTING_PROCESSES}

[파트 1](01_method.md)에서 수정한 명령을 다시 살펴보겠습니다:

```bash
{MODIFIED_COMMAND}
```

{EXPLAIN_DIFFERENCES_FROM_PART_2}

### 1.1. {MODIFICATION_SUBSTEP}

{INSTRUCTIONS_WITH_BEFORE_AFTER_TABS}

=== "후"

    ```groovy title="modules/{MODULE_FILE}.nf" linenums="{N}" hl_lines="{LINES}"
    {UPDATED_CODE}
    ```

=== "전"

    ```groovy title="modules/{MODULE_FILE}.nf" linenums="{N}" hl_lines="{LINES}"
    {ORIGINAL_CODE}
    ```

### 1.2. 워크플로우를 실행하여 수정 사항 확인

```bash
nextflow run {DOMAIN_DIR}.nf -profile test
```

??? success "명령 출력"

    ```console
    {EXPECTED_OUTPUT}
    ```

{VERIFY_OUTPUT_INSTRUCTIONS}

### 핵심 정리

워크플로우 동작을 조정하기 위해 프로세스 명령과 출력을 수정하는 방법을 학습했습니다.

### 다음 단계

다중 샘플 집계 단계를 추가합니다.

---

## 2. {AGGREGATION_STEP_TITLE}

{DESCRIPTION_OF_AGGREGATION}

### 2.1. 집계 모듈 작성

{MODULE_INSTRUCTIONS}

### 2.2. 샘플별 출력을 수집하여 집계 프로세스에 전달

{INSTRUCTIONS_USING_COLLECT_OPERATOR}

### 2.3. 완성된 워크플로우 실행

```bash
nextflow run {DOMAIN_DIR}.nf -profile test
```

??? success "명령 출력"

    ```console
    {EXPECTED_OUTPUT}
    ```

{VERIFY_OUTPUT_INSTRUCTIONS}

---

### 핵심 정리

샘플을 개별적으로 처리하고 모든 샘플에 걸쳐 결과를 집계하는 완전한 파이프라인을 구축했습니다.
`collect()`와 같은 채널 연산자를 사용하여 다중 샘플 분석을 위해 샘플별 출력을 집계하는 방법을 학습했습니다.

### 다음 단계

이 과정을 완료하신 것을 축하합니다! [과정 요약](next_steps.md)으로 이동하여 학습한 내용을 복습하고 다음 단계를 살펴보세요.
