# 파트 1: 컴퓨팅 환경에 맞게 적용하기

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } AI 지원 번역 - [자세히 알아보기 및 개선 제안](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>


[Nextflow Run](../nextflow_run/index.md)에서는 파이프라인의 입력, 매개변수, 출력을 설정했습니다.
이 과정에서는 나머지 절반을 다룹니다. 워크플로우 코드를 변경하지 않고, 파이프라인이 실행되는 컴퓨팅 환경에 맞게 파이프라인 실행을 적용하는 방법입니다.

!!! example "시나리오"

    노트북에서 Docker를 사용하여 파이프라인을 개발하고 테스트했습니다.
    이제 이를 전달해야 합니다. 협력자는 Conda만 설정되어 있고, 소속 기관의 HPC 클러스터는 자체 스케줄러와 리소스 제한을 통해 작업을 처리합니다.
    이러한 상황에서 파이프라인 자체를 다시 작성할 필요는 없습니다.

동일한 파이프라인 코드가 이 모든 환경에서 실행될 수 있는 이유는, 해당 내용이 워크플로우에 고정되어 있지 않기 때문입니다.
소프트웨어 패키징, 실행 플랫폼, 리소스 할당은 모두 설정을 통해 제어되며, 코드 위에 계층적으로 적용됩니다. 이 과정에서는 코드가 아닌 설정을 변경하여 동일한 파이프라인을 새로운 환경에 적용하는 방법을 다룹니다.

---

## 1. 소프트웨어 패키징 기술 선택

[Nextflow Run](../nextflow_run/index.md)에서 `nextflow.config`에 Docker의 대안으로 이미 설정된 `conda` 프로파일을 확인했습니다.
여기서는 동일한 전환을 직접 구성하고, process를 Conda와 함께 실제로 사용하기 위해 필요한 것이 무엇인지 살펴봅니다.

### 1.1. Docker 비활성화 및 Conda 활성화

`docker.enabled`를 `false`로 전환하고 Conda를 활성화하는 지시문을 추가합니다.

=== "후"

    ```groovy title="nextflow.config" linenums="1" hl_lines="1-2"
    docker.enabled = false
    conda.enabled = true
    ```

=== "전"

    ```groovy title="nextflow.config" linenums="1"
    docker.enabled = true
    ```

이렇게 하면 Nextflow가 Conda 패키지가 지정된 모든 process에 대해 Conda 환경을 생성하고 사용할 수 있습니다.
`cowpy` process에는 아직 Conda 패키지가 없으므로, 설정만으로 추가해 보겠습니다.

### 1.2. 설정을 통해 Conda 패키지 추가

`conda` 지시문은 `modules/cowpy.nf`에서 `container`가 설정된 것처럼 process 정의 자체에 설정할 수 있지만, 반드시 그럴 필요는 없습니다. `withName`을 사용하면 `cowpy` process에만 범위를 지정하여 설정에서 지시문을 설정할 수 있습니다.

=== "후"

    ```groovy title="nextflow.config" linenums="6" hl_lines="3-5"
    process {
        memory = 1.GB
        withName: 'cowpy' {
            conda = 'conda-forge::cowpy==1.1.5'
        }
    }
    ```

=== "전"

    ```groovy title="nextflow.config" linenums="6"
    process {
        memory = 1.GB
    }
    ```

이는 파이프라인 코드에 이미 있는 `container` 지시문을 대체하는 것이 아니라, 해당 코드를 전혀 수정하지 않고 대안을 추가하는 것입니다.

!!! tip "팁"

    [Seqera Containers](https://seqera.io/containers/) 검색은 컨테이너를 빌드할 계획이 없더라도 특정 도구의 Conda 패키지 URI를 조회하는 편리한 방법입니다.

### 1.3. 워크플로를 실행하여 Conda 사용 가능 여부 확인

```bash
nextflow run main.nf --batch conda
```

??? success "명령 출력"

    ```console
    N E X T F L O W   ~  version 26.04.4

    Launching `main.nf` [romantic_pike] revision: c3c85dec78

    executor >  local (8)
    [6d/d48030] sayHello (2)       | 3 of 3 ✔
    [f5/7a9d76] convertToUpper (1) | 3 of 3 ✔
    [1c/79b693] collectGreetings   | 1 of 1 ✔
    Creating env using conda: conda-forge::cowpy==1.1.5 [cache /workspaces/training/config-exec/work/conda/env-898314d566668b6587ad714ae06b8520]
    [bb/64b67c] cowpy              | 1 of 1 ✔

    Outputs:

      /workspaces/training/config-exec/results

      first_output:
        - full_pipeline/intermediates/Hola-output.txt
        - full_pipeline/intermediates/Bonjour-output.txt
        - full_pipeline/intermediates/Hello-output.txt

      uppercased:
        - full_pipeline/intermediates/UPPER-Hello-output.txt
        - full_pipeline/intermediates/UPPER-Bonjour-output.txt
        - full_pipeline/intermediates/UPPER-Hola-output.txt

      collected: full_pipeline/intermediates/COLLECTED-conda-output.txt

      batch_report: full_pipeline/conda-report.txt

      cowpy_art: full_pipeline/cowpy-COLLECTED-conda-output.txt
    ```

내부 동작 방식은 다르지만 Docker로 실행했을 때와 동일한 출력이 생성됩니다. Nextflow는 컨테이너 이미지를 가져오는 대신 Conda 패키지를 가져와 환경을 빌드합니다.

!!! info "정보"

    새로운 Conda 환경을 빌드하는 것은 처음에 컨테이너를 가져오는 것보다 시간이 조금 더 걸릴 수 있지만, 여기서 사용하는 패키지는 작으므로 빠르게 완료됩니다.

이제 이 과정의 나머지 부분을 위해 Docker로 다시 전환합니다.

```groovy title="nextflow.config" linenums="1"
docker.enabled = true
```

??? tip "Docker와 Conda 혼합 사용"

    이 설정들은 process별로 범위가 지정되므로 혼합하여 사용할 수 있습니다. 각 도구에 사용 가능한 것에 따라 일부 process는 Docker를, 다른 process는 Conda를 사용할 수 있습니다.
    동일한 process에 대해 `container` 지시문(파이프라인 코드에서)과 `conda` 지시문(설정에서) 모두 설정되어 있고 두 패키징 시스템이 모두 활성화된 경우, Nextflow는 컨테이너를 우선시합니다.

### 핵심 정리

process가 사용해야 할 소프트웨어 패키징 기술을 설정하는 방법과 Docker와 Conda 간에 전환하는 방법을 학습했습니다.

### 다음 단계

Nextflow가 작업을 실제로 실행하는 데 사용하는 실행 플랫폼을 변경하는 방법을 학습합니다.

---

## 2. 실행 플랫폼 선택

지금까지 실행한 모든 파이프라인은 로컬 executor를 사용했습니다. 각 작업은 Nextflow 자체와 동일한 머신에서 실행됩니다.
Nextflow는 사용 가능한 CPU와 메모리를 확인하고, 충분한 리소스가 확보될 때까지 작업을 대기시킵니다.

로컬 executor는 편리하지만 단일 머신을 넘어서는 확장이 불가능합니다.
Nextflow는 HPC 스케줄러(Slurm, LSF, SGE, PBS 등)와 클라우드 플랫폼(AWS Batch, Google Cloud Batch, Azure Batch, Kubernetes 등)을 포함한 [다양한 실행 백엔드](https://nextflow.io/docs/latest/executor.html)를 지원합니다.

### 2.1. 다른 백엔드 지정

executor는 `executor`라는 process 지시문으로 설정됩니다.
기본값은 `local`이므로 다음 설정이 암묵적으로 적용됩니다.

```groovy title="Built-in configuration"
process {
    executor = 'local'
}
```

다른 백엔드를 지정하려면 지시문을 원하는 executor로 설정합니다.

```groovy title="nextflow.config"
process {
    executor = 'slurm'
}
```

!!! warning "경고"

    교육 환경은 HPC 클러스터에 연결되어 있지 않으므로, 여기서는 실행할 수 없습니다.

### 2.2. 백엔드별 문법의 추상화

대부분의 HPC 플랫폼은 작업 제출 시 CPU, 메모리, 큐 이름 등의 리소스 요청을 각자의 문법으로 지정해야 합니다.
`my-science-work`라는 큐에서 8개의 CPU와 4GB RAM을 요청하는 동일한 내용이 스케줄러에 따라 완전히 다르게 표현됩니다.

??? abstract "예제"

    ```bash title="Config for SLURM / submit using sbatch"
    #SBATCH -o /path/to/my/task/directory/my-task-1.log
    #SBATCH --no-requeue
    #SBATCH -c 8
    #SBATCH --mem 4096M
    #SBATCH -p my-science-work
    ```

    ```bash title="Config for PBS / submit using qsub"
    #PBS -o /path/to/my/task/directory/my-task-1.log
    #PBS -j oe
    #PBS -q my-science-work
    #PBS -l nodes=1:ppn=8
    #PBS -l mem=4gb
    ```

    ```bash title="Config for SGE / submit using qsub"
    #$ -o /path/to/my/task/directory/my-task-1.log
    #$ -j y
    #$ -terse
    #$ -notify
    #$ -q my-science-work
    #$ -l slots=8
    #$ -l h_rss=4096M,mem_free=4096M
    ```

Nextflow는 이 모든 것을 추상화합니다. `cpus`, `memory`, `queue`와 같은 표준화된 속성을 한 번만 지정하면([process 지시문](https://nextflow.io/docs/latest/reference/process.html#process-directives)에서 전체 목록 확인 가능), Nextflow가 런타임에 이를 적절한 백엔드별 스크립트로 변환합니다.

### 2.3. Nextflow가 실제로 실행하는 내용 확인

이 변환은 단순한 설정 파일의 편의 기능이 아닙니다. 로컬 executor를 사용하는 지금도 직접 확인할 수 있는 구체적인 결과물이 있습니다.
[Nextflow Run, 섹션 1.3](../nextflow_run/01_run_nextflow.md#13-explore-the-work-directory)에서 `work/` 아래의 작업 디렉토리를 살펴보고 Nextflow가 실행한 정확한 명령인 `.command.sh`를 확인했습니다.
동일한 디렉토리에는 아직 확인하지 않은 파일이 있습니다. 바로 `.command.run`입니다.

```bash
cat work/0a/0df4a1*/.command.run
```

??? success "명령 출력 (발췌)"

    ```console
    #!/bin/bash
    ### ---
    ### name: 'convertToUpper (3)'
    ### container: 'null'
    ### outputs:
    ### - 'UPPER-Bonjour-output.txt'
    ### ...
    set -e
    set -u
    ...
    nxf_launch() {
        /bin/bash -ue /workspaces/training/nextflow-run/work/0a/0df4a1028c2001758b1841cff92fc7/.command.sh
    }
    ...
    ```

`.command.run`은 Nextflow가 실행을 위해 전달하는 실제 스크립트입니다.
환경 설정, 입출력 staging, 결과를 Nextflow에 보고하는 등 실제 실행에 필요한 모든 것으로 `.command.sh`를 감쌉니다.
`local` executor에서는 Nextflow가 동일한 머신에서 이 스크립트를 단순히 실행합니다.

이것이 바로 다른 `executor`를 설정할 때 변경되는 부분입니다.
Slurm이나 PBS와 같은 HPC 스케줄러의 경우, Nextflow는 동일한 종류의 래퍼 스크립트를 생성하고, [2.2](#22-backend-specific-syntax-is-abstracted-away)에서 확인한 스케줄러별 헤더(`cpus`, `memory`, `queue` 설정에서 변환된)를 추가한 후, 결과를 해당 스케줄러의 제출 명령(예: Slurm의 경우 `sbatch`)에 전달합니다.
이후 Nextflow는 로컬 프로세스를 직접 감시하는 대신 스케줄러에 작업 상태를 폴링합니다.
클라우드 배치 백엔드는 제출 명령이 아닌 API 호출로 구동되므로 약간 다르게 동작하지만, 기본 개념은 동일합니다. 동일한 작업 스크립트가 실행되며, 변경되는 것은 실행 및 추적 방식뿐입니다.

### 핵심 정리

다양한 컴퓨팅 인프라를 대상으로 executor를 변경하는 방법, Nextflow가 백엔드별 제출 문법을 추상화하는 방식, 그리고 작업이 다른 백엔드에서 실행될 때 실제로 어떤 일이 일어나는지 학습했습니다.

### 다음 단계

[파트 2](./02_resources_and_retries.md)로 이동하여 컴퓨팅 리소스를 프로파일링하고 할당하는 방법, 그리고 재시도로 작업 실패를 처리하는 방법을 학습합니다.

---

## 요약

이 파트에서 학습한 내용:

- Docker와 Conda 간 소프트웨어 패키징 기술 전환
- process 정의에 `conda` 지시문 추가
- `executor` 지시문으로 실행 플랫폼 변경
- Nextflow가 작업에 대해 실제로 생성하고 실행하는 내용 확인, 그리고 executor에 따라 어떻게 변경되는지 이해
