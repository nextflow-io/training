# 파트 3: 프로파일을 사용하여 설정 전환하기

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } AI 지원 번역 - [자세히 알아보기 및 개선 제안](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>


[파트 1](./01_packaging_and_execution.md)과 [파트 2](./02_resources_and_retries.md)에 걸쳐 소프트웨어 패키징, 실행 플랫폼, 리소스 할당 등 몇 가지 설정 옵션을 추가했습니다.
실제 환경에서는 개발용 노트북과 운영용 HPC 클러스터처럼, 실행 환경에 따라 이러한 옵션 세트 전체를 전환해야 하는 경우가 많습니다.

Nextflow에서는 다양한 설정을 기술하는 [프로파일](https://nextflow.io/docs/latest/config.html#profiles)을 원하는 만큼 정의하고, 실행 시 단일 플래그로 하나 또는 여러 개를 선택할 수 있습니다.

이미 하나를 사용해 보셨습니다. [Nextflow Run](../nextflow_run/index.md)의 `test` 프로파일은 입력 매개변수를 소규모의 잘 정의된 세트로 재정의합니다.
이제 직접 인프라 프로파일을 만들고 이와 결합해 봅니다.

---

## 1. 다양한 환경을 위한 프로파일 생성

### 1.1. 프로파일 설정

`nextflow.config`에 두 개의 프로파일을 추가합니다. 하나는 Docker를 사용하는 일반 노트북용이고, 다른 하나는 Slurm 스케줄러와 Conda를 사용하는 대학 HPC 클러스터용입니다.

=== "후"

    ```groovy title="nextflow.config" linenums="35" hl_lines="10-19"
    profiles {
        test {
            params.input = 'data/greetings.csv'
            params.batch = 'test'
            params.character = 'tux'
        }
        conda {
            docker.enabled = false
            conda.enabled = true
        }
        my_laptop {
            process.executor = 'local'
            docker.enabled = true
        }
        univ_hpc {
            process.executor = 'slurm'
            conda.enabled = true
            process.resourceLimits = [
                memory: 750.GB,
                cpus: 200,
                time: 30.d
            ]
        }
    }
    ```

=== "전"

    ```groovy title="nextflow.config" linenums="35"
    profiles {
        test {
            params.input = 'data/greetings.csv'
            params.batch = 'test'
            params.character = 'tux'
        }
        conda {
            docker.enabled = false
            conda.enabled = true
        }
    }
    ```

`univ_hpc` 프로파일은 공유 HPC 인프라에서 일반적으로 요구되는 리소스 제한도 설정합니다.

### 1.2. 프로파일을 사용하여 워크플로우 실행

실행 시 `-profile`로 프로파일을 선택합니다.

```bash
nextflow run main.nf -profile my_laptop
```

??? success "명령 출력"

    ```console
    N E X T F L O W   ~  version 26.04.4

    Launching `main.nf` [nice_kimura] revision: c3c85dec78

    executor >  local (8)
    [df/b86b2f] sayHello (2)       | 3 of 3 ✔
    [b2/e2cb9c] convertToUpper (3) | 3 of 3 ✔
    [af/2ecff7] collectGreetings   | 1 of 1 ✔
    [c6/713c96] cowpy              | 1 of 1 ✔

    Outputs:

      /workspaces/training/config-exec/results

      first_output:
        - full_pipeline/intermediates/Bonjour-output.txt
        - full_pipeline/intermediates/Hello-output.txt
        - full_pipeline/intermediates/Hola-output.txt

      uppercased:
        - full_pipeline/intermediates/UPPER-Bonjour-output.txt
        - full_pipeline/intermediates/UPPER-Hello-output.txt
        - full_pipeline/intermediates/UPPER-Hola-output.txt

      collected: full_pipeline/intermediates/COLLECTED-batch-output.txt

      batch_report: full_pipeline/batch-report.txt

      cowpy_art: full_pipeline/cowpy-COLLECTED-batch-output.txt
    ```

!!! warning "경고"

    `univ_hpc` 프로파일은 교육 환경에서 실행할 수 없습니다. Slurm 스케줄러를 사용할 수 없기 때문입니다.

항상 함께 사용해야 하는 다른 설정이 있다면 해당 프로파일에 추가하세요.
필요한 다른 조합을 묶기 위해 추가 프로파일을 만들 수도 있습니다.

### 1.3. 여러 프로파일로 실행

프로파일은 상호 배타적이지 않습니다.
`-profile <profile1>,<profile2>` 형식으로 여러 프로파일을 동시에 활성화할 수 있습니다.
`my_laptop`과 Nextflow Run에서 이미 알고 있는 `test` 프로파일을 결합합니다.

```bash
nextflow run main.nf -profile my_laptop,test
```

??? success "명령 출력"

    ```console
    N E X T F L O W   ~  version 26.04.4

    Launching `main.nf` [romantic_volhard] revision: c3c85dec78

    executor >  local (8)
    [22/81de4f] sayHello (3)       | 3 of 3 ✔
    [93/e16b04] convertToUpper (1) | 3 of 3 ✔
    [ca/3e6a91] collectGreetings   | 1 of 1 ✔
    [16/0cc7e3] cowpy              | 1 of 1 ✔

    Outputs:

      /workspaces/training/config-exec/results

      first_output:
        - full_pipeline/intermediates/Hola-output.txt
        - full_pipeline/intermediates/Hello-output.txt
        - full_pipeline/intermediates/Bonjour-output.txt

      uppercased:
        - full_pipeline/intermediates/UPPER-Hello-output.txt
        - full_pipeline/intermediates/UPPER-Bonjour-output.txt
        - full_pipeline/intermediates/UPPER-Hola-output.txt

      collected: full_pipeline/intermediates/COLLECTED-test-output.txt

      batch_report: full_pipeline/test-report.txt

      cowpy_art: full_pipeline/cowpy-COLLECTED-test-output.txt
    ```

개별 파일 이름에 `test` 프로파일의 `batch = 'test'`가 올바르게 반영되어 있습니다(`COLLECTED-test-output.txt` 등).

동일한 옵션을 설정하는 프로파일을 결합하면, Nextflow는 파일에서 나중에 읽히는 값, 즉 파일에서 더 뒤에 위치한 값을 사용하여 충돌을 해결합니다.
충돌하는 설정이 완전히 다른 설정 소스에서 온 경우에는 표준 [우선순위](https://www.nextflow.io/docs/latest/config.html)가 적용됩니다.

### 핵심 정리

인프라별 설정을 묶는 프로파일을 정의하고, 실행 시 `-profile`로 선택하며, 단일 실행에서 여러 프로파일을 결합하고, 둘 이상의 프로파일이 동일한 옵션을 설정할 때 Nextflow가 충돌을 해결하는 방법을 학습했습니다.

### 다음 단계

실행 전에 최종적으로 결정된 설정을 검사하는 방법을 학습합니다.

---

## 2. 결정된 설정 검사

[Nextflow Run](../nextflow_run/02_configure_pipeline.md)에서 `nextflow config -profile test`를 사용하여 단일 프로파일이 어떻게 결정되는지 확인했습니다.
이 명령은 여러 프로파일을 결합할 때 특히 유용합니다. 방금 확인했듯이, 두 프로파일이 동일한 옵션을 설정하면 어떤 값이 실제로 적용되는지 직접 파악하기 어려울 수 있습니다.
`nextflow config` 명령은 파이프라인을 실행하지 않고도 이 모든 것을 해결해 줍니다.

### 2.1. 기본 설정 확인

```bash
nextflow config
```

??? success "명령 출력"

    ```groovy
    params {
       input = 'data/greetings.csv'
       batch = 'batch'
       character = 'turkey'
    }

    docker {
       enabled = true
    }

    process {
       memory = '1 GB'
       withName:cowpy {
          conda = 'conda-forge::cowpy==1.1.5'
          memory = '2 GB'
          cpus = 2
       }
    }
    ```

이것이 추가 플래그 없이 파이프라인을 실행했을 때 적용되는 설정입니다.

### 2.2. 프로파일이 활성화된 상태에서 설정 확인

실제 실행에 사용할 프로파일과 동일하게 추가합니다.

```bash
nextflow config -profile my_laptop,test
```

??? success "명령 출력"

    ```groovy
    params {
       input = 'data/greetings.csv'
       batch = 'test'
       character = 'tux'
    }

    docker {
       enabled = true
    }

    process {
       memory = '1 GB'
       withName:cowpy {
          conda = 'conda-forge::cowpy==1.1.5'
          memory = '2 GB'
          cpus = 2
       }
       executor = 'local'
    }
    ```

두 결과를 비교하면 변경된 내용을 확인할 수 있습니다. `params.batch`, `params.character`, `process.executor`가 모두 `my_laptop,test` 프로파일을 반영하고 있습니다.
설정 레이어가 많은 파이프라인에서는 결정된 설정을 직접 파악하는 것이 번거롭고 오류가 발생하기 쉬우므로, 이 기능이 특히 유용합니다.

### 핵심 정리

`nextflow config`를 사용하여 실행 전에 프로파일 조합에 대해 최종적으로 결정된 설정을 검사하는 방법을 학습했습니다.

### 다음 단계

Nextflow 파이프라인 설정의 핵심 내용을 모두 다루었습니다.
다음 단계는 [과정 요약](next_steps.md)을 참조하세요.

---

## 요약

이 파트에서 학습한 내용:

- 인프라별 설정을 묶는 프로파일 정의
- 단일 실행에서 여러 프로파일 결합 및 프로파일 간 충돌 해결 방식 이해
- `nextflow config`를 사용하여 최종적으로 결정된 설정 검사
