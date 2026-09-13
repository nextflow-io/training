# 파트 2: 파이프라인 설정

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } AI 지원 번역 - [자세히 알아보기 및 개선 제안](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

[파트 1](./01_run_nextflow.md)에서는 컨테이너를 사용하여 여러 입력을 병렬로 처리하는 완전한 다단계 파이프라인을 실행했습니다.
이제 `nextflow.config`를 사용하여 파이프라인 동작을 설정하는 방법을 살펴봅니다. 먼저 이미 제공된 설정 파일을 검토하고, 설정을 제공하는 다른 방법들을 살펴본 후, 출력 결과가 게시되는 방식과 위치를 제어하는 방법을 다룹니다.

---

## 1. 주요 설정 파일 검토

Nextflow는 작업 디렉토리에서 `nextflow.config`를 자동으로 인식하고, 모든 실행에 해당 설정을 적용합니다.

제공된 설정 파일은 소프트웨어 패키징, 프로세스 설정, 파이프라인 매개변수, 실행 프로파일의 네 가지 영역을 다룹니다.

??? full-code "nextflow.config"

    ```groovy title="nextflow.config" linenums="1"
    /*
     * 소프트웨어 패키징
     */
    docker.enabled = true

    /*
     * 프로세스 설정
     */
    process {
        cpus = 1
        memory = 1.GB
    }

    /*
     * 파이프라인 매개변수
     */
    params {
        input = 'data/greetings.csv'
        batch = 'batch'
        character = 'turkey'
    }

    /*
     * 프로파일
     */
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

각 항목을 살펴본 후, 프로파일을 실제로 사용하여 파이프라인을 실행합니다.

!!! note "참고"

    이 설정은 단일 머신에서의 로컬 실행을 다룹니다.
    Nextflow는 HPC 스케줄러(SLURM, PBS, LSF)와 클라우드 executor(AWS Batch, Google Cloud Batch, Azure Batch)도 지원하며, 모두 동일한 `nextflow.config` 메커니즘으로 설정합니다.
    이러한 옵션에 대한 전체 안내는 [Execution Config](../execution_config/01_packaging_and_execution.md) 과정의 [파트 1: 컴퓨팅 환경에 맞게 적용하기](../execution_config/index.md)를 참조하세요.

### 1.1. 소프트웨어 패키징

소프트웨어 패키징은 Nextflow가 프로세스에 필요한 실제 도구를 제공하는 방식으로, 컨테이너 이미지, Conda 환경 등을 사용할 수 있습니다.

```groovy title="nextflow.config" linenums="1"
/*
 * 소프트웨어 패키징
 */
docker.enabled = true
```

이 줄은 모든 프로세스에 Docker를 활성화합니다.
`container` 지시문을 선언한 프로세스는 지정된 이미지 내에서 실행됩니다.

### 1.2. 프로세스 설정

프로세스는 `sayHello`나 `cowpy`처럼 파이프라인의 단일 단계입니다.
Nextflow를 사용하면 각 프로세스의 실행 방식에 대한 여러 설정을 지정할 수 있습니다. CPU와 메모리 할당량, 사용할 컨테이너 또는 Conda 환경 등을 설정할 수 있습니다.

```groovy title="nextflow.config" linenums="6"
/*
 * 프로세스 설정
 */
process {
    cpus = 1
    memory = 1.GB
}
```

이 설정은 모든 프로세스를 단일 CPU와 1 GB 메모리로 제한합니다.

Nextflow에서는 개별 프로세스나 프로세스 그룹에 서로 다른 값을 설정할 수도 있습니다. 자세한 내용은 [Execution Config](../execution_config/02_resources_and_retries.md#12-set-resource-allocations-for-a-specific-process) 과정의 [파트 2: 컴퓨팅 리소스 및 오류 관리](../execution_config/index.md)에서 학습할 수 있습니다.

### 1.3. 파이프라인 매개변수

매개변수는 파이프라인의 명령줄 입력으로, 이미 명령줄에서 직접 설정해 본 `--input`, `--batch`, `--character` 플래그와 동일합니다.
여기서 기본값을 설정하면 매번 입력할 필요가 없습니다. 이 파트의 뒷부분에서 살펴보겠지만, 매개변수를 제공하는 다른 방법도 있습니다.

```groovy title="nextflow.config" linenums="14"
/*
 * 파이프라인 매개변수
 */
params {
    input = 'data/greetings.csv'
    batch = 'batch'
    character = 'turkey'
}
```

이 기본값은 명령줄에서 매개변수가 제공되지 않을 때 적용되므로, 플래그 없이 `nextflow run main.nf`를 실행해도 정상적으로 동작합니다.

### 1.4. 프로파일

프로파일을 사용하면 여러 설정을 하나의 이름으로 묶을 수 있습니다. 매번 값을 직접 변경하는 대신, 플래그 하나로 전체 설정을 전환할 수 있습니다.

```groovy title="nextflow.config" linenums="23"
/*
 * 프로파일
 */
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

`test` 프로파일은 세 가지 매개변수를 재정의하여 소규모의 명확한 입력 세트로 파이프라인을 실행합니다. 모든 nf-core 파이프라인에는 빠른 검증을 위한 이 프로파일이 포함되어 있으며, 직접 작성하는 파이프라인에도 따를 만한 좋은 관례입니다.

`conda` 프로파일은 소프트웨어 패키징을 Docker에서 Conda로 전환합니다.

프로파일은 명령줄에서 `-profile <name>`을 전달하여 활성화합니다.

`test` 프로파일을 사용해 봅니다.

```bash
nextflow run main.nf -profile test
```

??? success "명령 출력"

    ```console hl_lines="6 7 8 9"
    N E X T F L O W   ~  version 26.04.4

    Launching `main.nf` [reverent_heisenberg] revision: ce74f81996

    executor >  local (8)
    [3d/8a12c7] sayHello (3)       | 3 of 3 ✔
    [e7/e0934f] convertToUpper (2) | 3 of 3 ✔
    [1d/616569] collectGreetings   | 1 of 1 ✔
    [44/7d46cf] cowpy              | 1 of 1 ✔

    Outputs:

      /workspaces/training/nextflow-run/results

      first_output:
        - test/intermediates/Bonjour-output.txt
        - test/intermediates/Hello-output.txt
        - test/intermediates/Hola-output.txt

      uppercased:
        - test/intermediates/UPPER-Bonjour-output.txt
        - test/intermediates/UPPER-Hello-output.txt
        - test/intermediates/UPPER-Hola-output.txt

      collected: test/intermediates/COLLECTED-test-output.txt

      batch_report: test/test-report.txt

      cowpy_art: test/cowpy-COLLECTED-test-output.txt
    ```

파이프라인이 `batch = 'test'`와 `character = 'tux'`로 실행됩니다.
`results/test/`를 확인하면 배치 이름이 디렉토리 경로에 포함되어 있고, ASCII 아트에 터키 대신 tux 펭귄이 표시됩니다.

!!! note "참고"

    여러 프로파일을 동시에 활성화할 수 있으며, `nextflow config -profile <name>,<name>`을 사용하면 실행 전에 최종 적용된 결과를 확인할 수 있습니다.
    프로파일 조합 방법과 Nextflow가 프로파일 간 충돌을 해결하는 방식은 [Execution Config](../execution_config/03_profiles.md) 과정의 [파트 3: 프로파일을 사용하여 설정 전환하기](../execution_config/index.md)에서 자세히 다룹니다.

### 핵심 정리

`nextflow.config` 파일의 주요 구성 요소가 무엇을 하는지, 그리고 프로파일을 활성화하는 방법을 학습했습니다.

### 다음 단계

주요 `nextflow.config` 파일을 수정하지 않고 설정 값을 제공하는 다른 방법을 학습합니다. 개별 실행을 설정하거나 정확한 설정 세트를 다른 사람과 공유할 때 유용합니다.

---

## 2. 보조 파일을 통한 설정 제공

`nextflow.config`에서 기본값을 설정하는 방식은 거의 변경되지 않는 값에 적합합니다.
Nextflow는 두 가지 추가적인 방법도 제공합니다. 특정 환경에 맞게 실행을 조정하기 위한 실행별 설정 파일과, 정확한 입력 값 세트를 공동 작업자와 공유하기 위한 매개변수 파일입니다.

### 2.1. 실행별 설정 파일 사용

Docker가 없는 머신으로 파이프라인을 이전하면서 각 프로세스에 더 많은 리소스를 할당하고 싶다고 가정합니다.
필요한 재정의 설정만 포함하는 새 설정 파일을 생성합니다.

```groovy title="custom.config" linenums="1"
process {
    cpus = 2
    memory = 2.GB
}

docker.enabled = false
conda.enabled = true
```

`-c` 옵션으로 주요 파이프라인과 함께 전달합니다.

```bash
nextflow run main.nf -c custom.config
```

??? success "명령 출력"

    ```console hl_lines="6 7 8 9"
    N E X T F L O W   ~  version 26.04.4

    Launching `main.nf` [exotic_cray] revision: ce74f81996

    executor >  local (8)
    [77/e02315] sayHello (1)       | 3 of 3 ✔
    [a6/ccf44b] convertToUpper (3) | 3 of 3 ✔
    [56/fd1296] collectGreetings   | 1 of 1 ✔
    [2c/205a94] cowpy              | 1 of 1 ✔

    Outputs:

      /workspaces/training/nextflow-run/results

      first_output:
        - batch/intermediates/Bonjour-output.txt
        - batch/intermediates/Hello-output.txt
        - batch/intermediates/Hola-output.txt

      uppercased:
        - batch/intermediates/UPPER-Bonjour-output.txt
        - batch/intermediates/UPPER-Hello-output.txt
        - batch/intermediates/UPPER-Hola-output.txt

      collected: batch/intermediates/COLLECTED-batch-output.txt

      batch_report: batch/batch-report.txt

      cowpy_art: batch/cowpy-COLLECTED-batch-output.txt
    ```

Nextflow는 `custom.config`를 파이프라인의 `nextflow.config` 위에 병합하므로, 모든 프로세스가 기본값 대신 2개의 CPU와 2 GB 메모리를 사용하고 Docker 대신 Conda로 실행됩니다.
`cowpy`는 컨테이너와 함께 Conda 패키지가 선언된 유일한 프로세스이므로, Nextflow가 실제로 환경을 빌드하는 것을 확인할 수 있습니다.

```console
Creating env using conda: conda-forge::cowpy==1.1.5 [cache /path/to/work/conda/env-898314d566668b6587ad714ae06b8520]
```

파이프라인 매개변수는 건드리지 않고 리소스 할당과 패키징만 재정의하는 소규모 파일은 nf-core 파이프라인이 기관별 설정에서 기대하는 패턴과 정확히 일치합니다.
실제 사례는 [nf-core/configs](https://github.com/nf-core/configs) 저장소를 참조하세요.

이 방법을 사용하면 기존 설정을 변경하지 않고도 파이프라인을 새로운 환경에 유연하게 적용할 수 있습니다.

### 2.2. 매개변수 파일 사용

정확한 실행 매개변수 세트를 공동 작업자와 공유하거나 논문에 기록해야 하는 경우를 가정합니다.

Nextflow는 YAML 또는 JSON 형식의 [매개변수 파일](https://nextflow.io/docs/latest/config.html#parameter-file)을 지원합니다. 이는 정확하고 재현 가능한 값 세트를 배포하는 더 간단한 방법입니다.

`test-params.yaml`이라는 매개변수 파일이 작업 디렉토리에 이미 제공되어 있습니다.

```yaml title="test-params.yaml" linenums="1"
input: "data/greetings.csv"
batch: "yaml"
character: "stegosaurus"
```

이 파일은 Groovy가 아닌 일반 YAML이므로, `nextflow.config`에서 사용하는 등호(`=`) 대신 콜론(`:`)을 사용합니다.

!!! info "정보"

    JSON 버전인 `test-params.json`도 제공됩니다. 직접 사용해 보세요. 전달 방법은 동일합니다.

`-params-file`로 파일을 전달합니다.

```bash
nextflow run main.nf -params-file test-params.yaml
```

??? success "명령 출력"

    ```console
    N E X T F L O W   ~  version 26.04.4

    Launching `main.nf` [sharp_faraday] revision: ce74f81996

    executor >  local (8)
    [1c/9ff63e] sayHello (1)       | 3 of 3 ✔
    [3b/bb5691] convertToUpper (2) | 3 of 3 ✔
    [cd/2c1f6e] collectGreetings   | 1 of 1 ✔
    [89/c333bc] cowpy              | 1 of 1 ✔

    Outputs:

      /workspaces/training/nextflow-run/results

      first_output:
        - yaml/intermediates/Bonjour-output.txt
        - yaml/intermediates/Hello-output.txt
        - yaml/intermediates/Hola-output.txt

      uppercased:
        - yaml/intermediates/UPPER-Bonjour-output.txt
        - yaml/intermediates/UPPER-Hello-output.txt
        - yaml/intermediates/UPPER-Hola-output.txt

      collected: yaml/intermediates/COLLECTED-yaml-output.txt

      batch_report: yaml/yaml-report.txt

      cowpy_art: yaml/cowpy-COLLECTED-yaml-output.txt
    ```

??? abstract "파일 내용"

    ```console title="results/yaml/cowpy-COLLECTED-yaml-output.txt"
     _________
    / BONJOUR \
    | HOLA    |
    \ HELLO   /
     ---------
    \                             .       .
     \                           / `.   .' "
      \                  .---.  <    > <    >  .---.
       \                 |    \  \ - ~ ~ - /  /    |
             _____          ..-~             ~-..-~
            |     |   \~~~\.'                    `./~~~/
           ---------   \__/                        \__/
          .'  O    \     /               /       \  "
         (_____,    `._.'               |         }  \/~~~/
          `----.          /       }     |        /    \__/
                `-.      |       /      |       /      `. ,~~|
                    ~-.__|      /_ - ~ ^|      /- _      `..-'
                        |     /        |     /     ~-.     `-. _  _  _
                        |_____|        |_____|         ~ - . _ _ _ _ _>
    ```

매개변수 파일은 파이프라인에 매개변수가 많을 때 특히 유용합니다. 워크플로우 스크립트를 변경하거나 긴 명령줄을 입력하지 않고도 모든 매개변수를 한 번에 제공할 수 있으며, 결과물과 함께 배포하기도 쉽습니다.

### 핵심 정리

설정을 제공하는 두 가지 추가 방법을 학습했습니다. 새로운 환경에 맞게 실행을 조정하기 위한 실행별 설정 파일과, 정확하고 재현 가능한 입력 값을 공유하기 위한 매개변수 파일입니다.

### 다음 단계

파이프라인 출력 결과가 게시되는 방식과 위치를 제어하는 방법을 학습합니다.

---

## 3. 파이프라인 출력 관리

파이프라인 작성자가 코드에서 출력 구조를 결정하지만, 해당 코드를 수정하지 않고도 출력이 저장되는 위치와 방식을 제어할 수 있습니다.
Nextflow는 설정 수준에서 이를 제어하는 방법을 제공합니다. 기본 출력 디렉토리를 설정하고, 파일을 복사할지 심볼릭 링크로 연결할지 선택할 수 있습니다.

### 3.1. 출력 디렉토리 맞춤화

기본적으로 Nextflow는 `results/` 아래에 출력을 게시합니다.
`-output-dir`(또는 단축형 `-o`)을 사용하여 다른 위치를 지정합니다.

```bash
nextflow run main.nf -output-dir outputs
```

??? success "명령 출력"

    ```console hl_lines="6 7 8 9"
    N E X T F L O W   ~  version 26.04.4

    Launching `main.nf` [serene_kimura] revision: ce74f81996

    executor >  local (8)
    [31/df5c15] sayHello (3)       | 3 of 3 ✔
    [08/8bb2e5] convertToUpper (1) | 3 of 3 ✔
    [e5/5814da] collectGreetings   | 1 of 1 ✔
    [cf/8ab8c6] cowpy              | 1 of 1 ✔

    Outputs:

      /workspaces/training/nextflow-run/outputs

      first_output:
        - batch/intermediates/Bonjour-output.txt
        - batch/intermediates/Hello-output.txt
        - batch/intermediates/Hola-output.txt

      uppercased:
        - batch/intermediates/UPPER-Bonjour-output.txt
        - batch/intermediates/UPPER-Hello-output.txt
        - batch/intermediates/UPPER-Hola-output.txt

      collected: batch/intermediates/COLLECTED-batch-output.txt

      batch_report: batch/batch-report.txt

      cowpy_art: batch/cowpy-COLLECTED-batch-output.txt
    ```

??? abstract "디렉토리 내용"

    ```console
    outputs/batch
    ├── batch-report.txt
    ├── cowpy-COLLECTED-batch-output.txt
    └── intermediates
        ├── Bonjour-output.txt
        ├── COLLECTED-batch-output.txt
        ├── Hello-output.txt
        ├── Hola-output.txt
        ├── UPPER-Bonjour-output.txt
        ├── UPPER-Hello-output.txt
        └── UPPER-Hola-output.txt
    ```

출력이 기본값인 `results/batch/` 대신 `outputs/batch/` 아래에 저장됩니다.
`batch/`와 `intermediates/` 같은 하위 디렉토리 구조는 파이프라인 코드가 결정하며, `-output-dir`은 해당 구조가 시작되는 위치만 제어합니다.

`-output-dir`은 `outputDir` 설정 옵션의 명령줄 단축키이므로, 설정이 가능한 모든 위치에서 사용할 수 있습니다. `nextflow.config`에 직접 지정하거나, 프로파일 내부에 넣거나, 이 파트에서 앞서 사용한 것처럼 `-c` 오버레이 파일에 넣을 수 있습니다.
예를 들어, 다음 코드는 명령줄 대신 `nextflow.config`에 직접 동일한 설정을 적용하는 방법을 보여줍니다.

```groovy title="nextflow.config"
outputDir = 'outputs'
```

설정 옵션을 사용할 수 있는 전체 위치 목록은 Nextflow 참조 문서의 [Configuration file](https://nextflow.io/docs/latest/config.html)을 참조하세요.

### 3.2. 출력 게시 방식 선택

기본적으로 Nextflow는 실제 복사본이 아닌 `work/` 아래의 출력 위치를 가리키는 심볼릭 링크로 출력을 게시합니다.

```console
$ ls -l results/batch/intermediates/Hello-output.txt
lrwxr-xr-x  ...  Hello-output.txt -> /workspaces/training/nextflow-run/work/b7/b8c4d1.../Hello-output.txt
```

파이프라인 작성자는 워크플로우 코드에서 각 프로세스에 대해 게시 모드를 `'copy'` 또는 `'move'`로 설정할 수 있습니다.
일반적으로 파이프라인의 최종 출력에는 이 설정을 적용하고, 전체 파이프라인 실행 후 삭제 가능한 중간 파일에는 기본 `'symlink'` 동작을 유지합니다.

이 방식은 디스크의 데이터 중복을 방지하지만, `-resume` 기능을 사용하는 능력을 잃지 않으려면 `work/` 아래의 작업 디렉토리를 삭제할 수 없습니다.
모든 출력 파일을 실제로 복사하려면 파이프라인 설정에서 [`workflow.output.mode`](https://nextflow.io/docs/latest/reference/config.html#workflow)를 `'copy'`로 설정하세요. (`-output-dir`과 달리 이 설정에는 명령줄 플래그가 없으며, 설정 파일에서만 지정할 수 있습니다.)

`nextflow.config`에 설정해 봅니다.

=== "후"

    ```groovy title="nextflow.config" hl_lines="7"
    params {
        input = 'data/greetings.csv'
        batch = 'batch'
        character = 'turkey'
    }

    workflow.output.mode = 'copy'
    ```

=== "전"

    ```groovy title="nextflow.config"
    params {
        input = 'data/greetings.csv'
        batch = 'batch'
        character = 'turkey'
    }
    ```

그런 다음 출력의 차이를 확인할 수 있도록 배치 이름을 변경하여 파이프라인을 실행합니다.

```bash
nextflow run main.nf --batch withmode
```

??? success "명령 출력"

    ```console hl_lines="6 7 8 9"
    N E X T F L O W   ~  version 26.04.4

    Launching `main.nf` [angry_noether] revision: ce74f81996

    executor >  local (8)
    [41/2b478d] sayHello (3)       | 3 of 3 ✔
    [bf/dd2840] convertToUpper (2) | 3 of 3 ✔
    [ea/364e97] collectGreetings   | 1 of 1 ✔
    [10/76fe7b] cowpy              | 1 of 1 ✔

    Outputs:

      /workspaces/training/nextflow-run/results

      first_output:
        - withmode/intermediates/Bonjour-output.txt
        - withmode/intermediates/Hello-output.txt
        - withmode/intermediates/Hola-output.txt

      uppercased:
        - withmode/intermediates/UPPER-Bonjour-output.txt
        - withmode/intermediates/UPPER-Hello-output.txt
        - withmode/intermediates/UPPER-Hola-output.txt

      collected: withmode/intermediates/COLLECTED-withmode-output.txt

      batch_report: withmode/withmode-report.txt

      cowpy_art: withmode/cowpy-COLLECTED-withmode-output.txt
    ```

이전과 같이 출력 파일 중 하나를 확인합니다.

```console
$ ls -l results/withmode/intermediates/Hello-output.txt
-rw-r--r--  ...  Hello-output.txt
```

이제 `work/`가 정리되어도 계속 사용할 수 있는 실제 독립 파일입니다.

!!! warning "경고"

    `workflow.output.mode` 설정은 파이프라인 코드에서 아직 모드가 설정되지 않은 출력에 대한 기본값만 채웁니다.
    작성자가 하드코딩한 모드는 어떤 값을 설정하더라도 재정의할 수 없습니다.

### 핵심 정리

파이프라인 코드를 수정하지 않고도 기본 출력 디렉토리를 맞춤화하고, 복사된 출력과 심볼릭 링크 출력 중에서 선택하는 방법을 학습했습니다.

### 다음 단계

[파트 3](./03_manage_executions.md)으로 이동하여 과거 실행 기록을 검토하고, 실행 보고서를 생성하고, 오래된 work 디렉토리를 정리하는 방법을 학습합니다.

---

## 요약

이 파트에서 학습한 내용:

- `nextflow.config`와 프로파일을 사용하여 파이프라인 동작 설정하기
- 실행별 설정 파일 또는 매개변수 파일을 통해 설정 제공하기
- 출력 디렉토리 맞춤화 및 복사된 출력과 심볼릭 링크 출력 중 선택하기
