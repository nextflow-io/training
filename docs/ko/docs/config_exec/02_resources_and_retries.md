# 파트 2: 컴퓨팅 리소스 및 실패 관리

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } AI 지원 번역 - [자세히 알아보기 및 개선 제안](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>


[파트 1](./01_packaging_and_execution.md)에서는 파이프라인의 작업이 실행되는 위치와 방식을 조정했습니다.
이번 파트에서는 각 작업에 할당되는 컴퓨팅 리소스의 양과, 최선의 추정으로 할당했음에도 작업이 실패했을 때의 처리 방법을 다룹니다.

---

## 1. 컴퓨팅 리소스 할당 제어

기본적으로 Nextflow는 `cpus` 지시문을 통해 각 프로세스에 CPU 1개를 할당하며, 별도로 설정하지 않는 한 메모리 제한을 적용하지 않습니다.

```groovy title="Built-in configuration"
process {
    cpus = 1
}
```

[Nextflow Run](../nextflow_run/index.md)에서 이미 학습했듯이, 이 파이프라인의 설정은 모든 프로세스에 대해 `memory`를 1 GB로 지정합니다.
그렇다면 실제로 사용할 값은 어떻게 결정해야 할까요?

### 1.1. 리소스 사용률 보고서 생성

[Nextflow Run](../nextflow_run/02_configure_pipeline.md)에서 `-with-report` 옵션으로 실행 보고서를 생성한 적이 있습니다.
프로세스가 실제로 필요로 하는 CPU와 메모리를 파악하는 방법도 동일합니다. 기본 할당값으로 워크플로우를 실행하고, 실제 사용량을 기록한 후 이를 기반으로 조정합니다.

```bash
nextflow run main.nf -with-report report-config-1.html
```

보고서는 브라우저에서 열 수 있는 HTML 파일입니다.
프로세스별 실행 시간과 리소스 사용률을 분석하여, 할당된 리소스 중 실제로 사용된 비율을 보여줍니다.
현재 기본값(CPU 1개, 메모리 1 GB)으로 `cowpy`를 실행했을 때의 결과는 다음과 같습니다.

| 지표             | 값     |
| ---------------- | ------ |
| CPU 사용률       | 116%   |
| 최대 메모리 사용 | 6.4 MB |
| 할당된 메모리    | 1 GB   |

`cowpy`는 1 GB 할당량의 1%도 사용하지 않습니다. `%cpu`가 100%를 초과하는 것은 컨테이너 내부에서 짧은 순간 동안 CPU 1개 이상의 처리 능력을 사용한다는 의미입니다.

사용 가능한 전체 기능 목록은 [Reports](https://nextflow.io/docs/latest/reports.html)를 참조하세요.

### 1.2. 특정 프로세스에 리소스 할당 설정

위 보고서에서 `cowpy`는 현재 할당량 내에서 여유롭게 실행되고 있습니다. 그러나 예를 들어 프로덕션 환경에서 더 큰 입력 데이터를 처리할 것으로 예상되는 경우, 더 넉넉한 여유를 주고 싶을 수 있습니다.
`withName`을 사용하면 단일 프로세스의 기본값을 재정의할 수 있습니다.

=== "후"

    ```groovy title="nextflow.config" linenums="6" hl_lines="5-6"
    process {
        memory = 1.GB
        withName: 'cowpy' {
            conda = 'conda-forge::cowpy==1.1.5'
            memory = 2.GB
            cpus = 2
        }
    }
    ```

=== "전"

    ```groovy title="nextflow.config" linenums="6"
    process {
        memory = 1.GB
        withName: 'cowpy' {
            conda = 'conda-forge::cowpy==1.1.5'
        }
    }
    ```

이 설정을 적용하면, `cowpy`를 제외한 모든 프로세스는 메모리 1 GB와 CPU 1개를 요청합니다. `cowpy`는 ([파트 1](./01_packaging_and_execution.md)의 `conda` 설정에 더하여) 메모리 2 GB와 CPU 2개를 요청합니다.

!!! info "정보"

    머신의 CPU 수가 적은데 프로세스당 높은 수의 CPU를 할당하면, Nextflow는 사용 가능한 CPU 수를 초과하여 요청하지 않으므로 작업 실행이 순차적으로 대기할 수 있습니다.

비교를 위해 다른 보고서 파일명을 지정하여 다시 실행합니다.

```bash
nextflow run main.nf -with-report report-config-2.html
```

??? success "명령 출력"

    ```console
    N E X T F L O W   ~  version 26.04.4

    Launching `main.nf` [voluminous_venter] revision: c3c85dec78

    executor >  local (8)
    [a1/0e96d4] sayHello (1)       | 3 of 3 ✔
    [a3/7173a3] convertToUpper (2) | 3 of 3 ✔
    [4f/a8ae3d] collectGreetings   | 1 of 1 ✔
    [91/3724f8] cowpy              | 1 of 1 ✔

    Outputs:

      /workspaces/training/config-exec/results

      first_output:
        - full_pipeline/intermediates/Bonjour-output.txt
        - full_pipeline/intermediates/Hola-output.txt
        - full_pipeline/intermediates/Hello-output.txt

      uppercased:
        - full_pipeline/intermediates/UPPER-Bonjour-output.txt
        - full_pipeline/intermediates/UPPER-Hello-output.txt
        - full_pipeline/intermediates/UPPER-Hola-output.txt

      collected: full_pipeline/intermediates/COLLECTED-batch-output.txt

      batch_report: full_pipeline/batch-report.txt

      cowpy_art: full_pipeline/cowpy-COLLECTED-batch-output.txt
    ```

`cowpy`에 대한 두 보고서를 비교하면 다음과 같습니다.

| 지표             | 전 (CPU 1개, 1 GB) | 후 (CPU 2개, 2 GB) |
| ---------------- | ------------------ | ------------------ |
| 최대 메모리 사용 | 6.4 MB             | 6.4 MB             |
| CPU 사용률       | 116%               | 118%               |

할당량을 두 배로 늘려도 실제 사용량은 전혀 변하지 않았습니다. 이는 원래의 1 GB / CPU 1개 할당이 이 간단한 워크로드에 이미 충분했음을 보여줍니다.
실제 파이프라인에서 의미 있는 데이터를 처리할 때는 프로세스 간에 수치가 크게 달라질 것입니다. 바로 이것이 할당량을 추측하는 대신 프로파일링을 먼저 수행해야 하는 이유입니다.

### 1.3. 리소스 한도 추가

컴퓨팅 인프라에 따라 요청 가능한 리소스에 클러스터 전체 상한선과 같은 엄격한 제약이 있을 수 있습니다.
`resourceLimits` 지시문을 사용하면 이러한 한도를 설정할 수 있습니다.

```groovy title="Syntax example"
process {
    resourceLimits = [
        memory: 750.GB,
        cpus: 200,
        time: 30.d
    ]
}
```

Nextflow는 이를 대상 executor가 요구하는 형식으로 변환합니다.
프로세스가 한도를 초과하는 리소스를 요청하면, 요청이 거부되는 대신 한도에 맞게 조정됩니다.

!!! warning "경고"

    이 기능은 HPC 인프라가 있어야 효과가 있으므로, 교육 환경에서는 실행할 수 없습니다.

??? info "기관별 참조 설정"

    nf-core 프로젝트는 전 세계 기관들이 공유하는 [설정 파일 모음](https://nf-co.re/configs/)을 관리하고 있으며, 다양한 HPC 및 클라우드 executor를 지원합니다.
    자신의 기관이 포함되어 있는지 여부와 관계없이 유용한 출발점이 됩니다.

### 핵심 정리

프로파일링 보고서를 생성하여 리소스 사용률을 평가하고, 특정 프로세스의 리소스 할당을 재정의하며, `resourceLimits`로 할당량에 상한선을 설정하는 방법을 학습했습니다.

### 다음 단계

리소스 할당 추정이 맞든 틀리든, 작업이 실패했을 때 파이프라인이 자동으로 복구되도록 만드는 방법을 학습합니다.

---

## 2. 재시도로 작업 실패 처리

프로파일링은 대부분의 경우 프로세스에 필요한 리소스를 알려주지만, 실제 워크로드는 다양합니다. 대부분의 입력에 충분한 할당량이 비정상적으로 큰 입력에는 부족할 수 있으며, 추정 자체가 틀릴 수도 있습니다.
단일 작업 실패로 전체 실행이 중단되는 대신, Nextflow는 실패한 작업을 자동으로 재시도할 수 있으며, 각 시도마다 더 많은 리소스를 제공하는 것도 가능합니다.

### 2.1. 실패한 작업 자동 재시도

이를 직접 확인하기 위해, `cowpy`의 메모리 할당을 실제 필요량보다 낮게 의도적으로 설정합니다. [1.1](#11-generate-a-resource-utilization-report)에서 최대 약 6.4 MB를 사용한다는 것을 확인했으므로, 6 MB는 충분하지 않을 것입니다.

=== "후"

    ```groovy title="nextflow.config" linenums="6" hl_lines="5-7"
    process {
        memory = 1.GB
        withName: 'cowpy' {
            conda = 'conda-forge::cowpy==1.1.5'
            memory = 6.MB
            errorStrategy = 'retry'
            maxRetries = 2
        }
    }
    ```

=== "전"

    ```groovy title="nextflow.config" linenums="6"
    process {
        memory = 1.GB
        withName: 'cowpy' {
            conda = 'conda-forge::cowpy==1.1.5'
            memory = 2.GB
            cpus = 2
        }
    }
    ```

`errorStrategy`는 작업이 실패했을 때 Nextflow가 취할 동작을 지정합니다. `'retry'`는 파이프라인 전체를 중단하는 대신 작업을 재제출합니다.
`maxRetries`는 Nextflow가 포기하기 전까지 허용되는 추가 시도 횟수를 제한합니다.

```bash
nextflow run main.nf
```

??? failure "명령 출력 (요약)"

    ```console
    N E X T F L O W   ~  version 26.04.4

    Launching `main.nf` [desperate_brazil] revision: c3c85dec78

    executor >  local (10)
    [67/fe1f49] sayHello (1)       | 3 of 3 ✔
    [8a/f13335] convertToUpper (1) | 3 of 3 ✔
    [39/1b24ed] collectGreetings   | 1 of 1 ✔
    [7a/d5eb6f] cowpy              | 0 of 1, retries: 2 ✘
    [9d/b79eb3] NOTE: Process `cowpy` terminated with an error exit status (137) -- Execution is retried (1)
    [6d/1d9d84] NOTE: Process `cowpy` terminated with an error exit status (137) -- Execution is retried (2)
    ERROR ~ Error executing process > 'cowpy'

    Caused by:
      Process `cowpy` terminated with an error exit status (137)

    Command executed:
      cat COLLECTED-batch-output.txt | cowpy -c "turkey" > cowpy-COLLECTED-batch-output.txt

    Command exit status:
      137

    Command output:
      (empty)

    Command error:
      /usr/local/bin/_activate_current_env.sh: line 35:    14 Killed                  micromamba activate "${ENV_NAME:-base}"

    Work dir:
      /workspaces/training/config-exec/work/7a/d5eb6feeac0eed18d95d3da7a7aeb4

    Tip: when you have fixed the problem you can continue the execution adding the option `-resume` to the run command line

    -- Check '.nextflow.log' file for details
    ```

종료 코드 137은 메모리 부족으로 인한 강제 종료를 나타내는 표준 신호입니다. 컨테이너에 `cowpy`를 실행하기에 충분한 메모리가 없었던 것입니다.
Nextflow는 작업을 두 번 재시도하여 총 세 번 시도했으며, 이는 `maxRetries = 2` 설정과 일치합니다.
시도 간에 메모리 할당이 변경되지 않았으므로 모든 시도가 동일한 한계에 부딪혔습니다. 재시도가 모두 소진되면 Nextflow는 실패를 전체적으로 보고하고 파이프라인을 중단하며, 0이 아닌 종료 상태로 종료됩니다.

근본적인 원인이 시도 간에 변경되지 않으면 재시도만으로는 아무것도 해결되지 않습니다.

### 2.2. 재시도마다 리소스 증가

프로세스 지시문 내에서 `task.attempt`는 현재 시도 번호를 나타내며, 1부터 시작합니다.
이를 closure에서 사용하면 재시도마다 리소스 할당을 단계적으로 늘릴 수 있습니다.

=== "후"

    ```groovy title="nextflow.config" linenums="6" hl_lines="5 7"
    process {
        memory = 1.GB
        withName: 'cowpy' {
            conda = 'conda-forge::cowpy==1.1.5'
            memory = { 6.MB * task.attempt }
            errorStrategy = 'retry'
            maxRetries = 3
        }
    }
    ```

=== "전"

    ```groovy title="nextflow.config" linenums="6"
    process {
        memory = 1.GB
        withName: 'cowpy' {
            conda = 'conda-forge::cowpy==1.1.5'
            memory = 6.MB
            errorStrategy = 'retry'
            maxRetries = 2
        }
    }
    ```

워크플로우를 다시 실행합니다.

```bash
nextflow run main.nf
```

??? success "명령 출력 (요약)"

    ```console
    N E X T F L O W   ~  version 26.04.4

    Launching `main.nf` [grave_joliot] revision: c3c85dec78

    executor >  local (9)
    [0f/211b8a] sayHello (2)       | 3 of 3 ✔
    [26/301ea2] convertToUpper (3) | 3 of 3 ✔
    [22/5a895b] collectGreetings   | 1 of 1 ✔
    [e1/beee86] cowpy              | 1 of 1, retries: 1 ✔
    [b6/7aed6a] NOTE: Process `cowpy` terminated with an error exit status (137) -- Execution is retried (1)

    Outputs:

      /workspaces/training/config-exec/results

      first_output:
        - full_pipeline/intermediates/Hola-output.txt
        - full_pipeline/intermediates/Hello-output.txt
        - full_pipeline/intermediates/Bonjour-output.txt

      uppercased:
        - full_pipeline/intermediates/UPPER-Hola-output.txt
        - full_pipeline/intermediates/UPPER-Hello-output.txt
        - full_pipeline/intermediates/UPPER-Bonjour-output.txt

      collected: full_pipeline/intermediates/COLLECTED-batch-output.txt

      batch_report: full_pipeline/batch-report.txt

      cowpy_art: full_pipeline/cowpy-COLLECTED-batch-output.txt
    ```

첫 번째 시도는 6 MB에서 여전히 실패하지만, 재시도는 12 MB(`6.MB * 2`)로 실행되어 성공하고, 파이프라인은 모든 출력을 게시하며 완료됩니다.

!!! warning "경고"

    파이프라인 전체가 성공했더라도 콘솔 출력에는 첫 번째 시도 실패를 보고하는 `NOTE:` 줄이 포함됩니다. Nextflow는 각 재시도를 개별적으로 기록하지만, 재시도된 실패는 전체 결과에 영향을 미치지 않습니다.
    실행이 실제로 성공했는지 확인하려면 `Outputs:` 요약 또는 명령의 종료 상태를 확인하세요.

더 고급 재시도 패턴(특정 오류 유형에 따른 확장 포함)은 Nextflow 문서의 [Dynamic computing resources](https://nextflow.io/docs/latest/process.html#dynamic-task-resources)를 참조하세요.

### 핵심 정리

파이프라인이 실패한 작업을 자동으로 재시도하도록 만드는 방법과, `task.attempt`를 사용하여 재시도마다 리소스 할당을 단계적으로 늘리는 방법을 학습했습니다.

### 다음 단계

[파트 3](./03_profiles.md)으로 이동하여, 이러한 설정을 전환 가능한 프로파일로 묶는 방법을 학습합니다.

---

## 요약

이번 파트에서 학습한 내용은 다음과 같습니다.

- 리소스 프로파일링 보고서를 생성하고 프로세스별 리소스 할당 설정
- `resourceLimits`로 리소스 요청에 상한선 적용
- `errorStrategy`와 `maxRetries`로 실패한 작업 자동 재시도
- `task.attempt`를 사용하여 재시도마다 리소스 할당 단계적 증가
