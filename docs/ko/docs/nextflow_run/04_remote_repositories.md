# 파트 4: 원격 파이프라인 실행

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } AI 지원 번역 - [자세히 알아보기 및 개선 제안](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

지금까지는 로컬에 저장된 워크플로우 스크립트를 실행했습니다.
실제 환경에서는 GitHub 등의 원격 저장소에 게시된 파이프라인을 직접 다운로드하지 않고 실행하는 경우가 많습니다.

Nextflow는 이를 간단하게 지원합니다. Git 저장소 URL에서 직접 파이프라인을 실행할 수 있습니다.

---

## 1. GitHub에서 파이프라인 실행

원격 파이프라인을 실행하는 기본 문법은 `nextflow run <repository>`입니다. 여기서 `<repository>`는 `nextflow-io/hello`와 같은 GitHub 저장소 경로, 전체 URL, 또는 GitLab, Bitbucket 등 다른 Git 호스팅 서비스의 경로가 될 수 있습니다.

### 1.1. 파이프라인 실행

Nextflow 공식 "hello" 데모 파이프라인을 실행합니다.
이 파이프라인은 이 과정에서 사용해 온 파이프라인과는 다른, 훨씬 단순한 파이프라인입니다. 이 교육 전반에서 사용된 "Hello" 파이프라인보다 이전에 만들어진 것으로, 몇 가지 하드코딩된 언어로 인사말을 출력하는 기능만 합니다. CSV 입력이나 ASCII 아트는 기대하지 마세요.

```bash
nextflow run nextflow-io/hello
```

??? success "명령 출력"

    ```console
    N E X T F L O W   ~  version 26.04.4

    Pulling nextflow-io/hello ...
     downloaded from https://github.com/nextflow-io/hello.git
    Launching `https://github.com/nextflow-io/hello` [sleepy_swanson] revision: 3c2cdc9823 [master]

    executor >  local (4)
    [ba/08236d] sayHello (4) | 4 of 4 ✔
    Ciao world!

    Hello world!

    Bonjour world!

    Hola world!
    ```

### 1.2. 파이프라인 캐시 위치 확인

원격 파이프라인을 처음 실행하면 Nextflow가 해당 파이프라인을 다운로드하여 로컬에 캐시합니다.
이후 실행 시에는 명시적으로 업데이트를 요청하지 않는 한 캐시된 버전을 재사용합니다.

기본적으로 Nextflow는 가져온 파이프라인을 `$NXF_HOME/assets` 아래에 저장합니다.
특정 파이프라인의 저장 위치와 사용 가능한 리비전을 확인하려면 Nextflow에 직접 조회합니다.

```bash
nextflow info nextflow-io/hello
```

??? success "명령 출력"

    ```console
     project name: nextflow-io/hello
     repository  : https://github.com/nextflow-io/hello
     local path  : /workspaces/.nextflow/assets/.repos/nextflow-io/hello
     main script : main.nf
     revisions   :
     > master (default)
       mybranch
       testing
       v1.1 [t]
       v1.2 [t]
       v1.3 [t]
    ```

    Nextflow는 이미 로컬에 체크아웃된 리비전에 `>`를 표시합니다. 나머지는 사용 가능하지만 아직 작업 복사본으로 가져오지 않은 상태입니다.

지금까지 가져온 모든 파이프라인 목록은 `nextflow list`로 확인할 수 있습니다.

```bash
nextflow list
```

??? success "명령 출력"

    ```console
    nextflow-io/hello
    ```

[Use nf-core](../nfcore_use/01_run_demo.md#12-retrieve-the-pipeline-code) 과정에서는 가져온 파이프라인의 소스 코드를 탐색하는 방법을 포함하여 이 캐싱 메커니즘을 더 자세히 다룹니다.

### 핵심 정리

GitHub 저장소에서 직접 파이프라인을 다운로드하지 않고 실행하는 방법과, 실행 후 로컬에서 파이프라인을 찾는 방법을 학습했습니다.

### 다음 단계

재현성을 위해 원격 파이프라인의 특정 버전을 고정하는 방법을 학습합니다.

---

## 2. 재현성을 위한 버전 지정

기본적으로 Nextflow는 기본 브랜치의 최신 리비전을 실행합니다.
`-r` 플래그를 사용하여 특정 버전(태그), 브랜치, 또는 커밋을 고정할 수 있습니다.

### 2.1. 특정 리비전 고정

```bash
nextflow run nextflow-io/hello -r v1.3
```

??? success "명령 출력"

    ```console
    N E X T F L O W   ~  version 26.04.4

    Pulling nextflow-io/hello:v1.3 ...
     downloaded from https://github.com/nextflow-io/hello.git
    Launching `https://github.com/nextflow-io/hello` [sick_carson] revision: 2ce0b0e294 [v1.3]

    executor >  local (4)
    [61/e11f77] sayHello (4) | 4 of 4 ✔
    Ciao world!

    Bonjour world!

    Hello world!

    Hola world!
    ```

Nextflow는 처음 요청할 때 해당 리비전을 가져오므로 `Pulling` 및 `downloaded from` 줄이 표시됩니다. 이후 동일한 리비전을 요청하면 바로 `Launching`으로 넘어갑니다.
정확한 리비전을 고정하는 것은 재현성을 위해 필수적입니다.
저장소에 변경 사항이 생기더라도 본인과 협업자가 동일한 파이프라인 코드를 실행한다는 것을 보장합니다.

### 2.2. 리비전은 해당 실행에만 적용됩니다

`-r`로 리비전을 고정하는 것은 해당 실행에만 영향을 미칩니다. 이후 `-r` 없이 `nextflow run`을 실행하면 이전에 고정한 리비전이 적용되지 않습니다.
`-r` 없이 파이프라인을 다시 실행합니다.

```bash
nextflow run nextflow-io/hello
```

??? success "명령 출력"

    ```console
    N E X T F L O W   ~  version 26.04.4

    Launching `https://github.com/nextflow-io/hello` [hungry_maxwell] revision: 3c2cdc9823 [master]

    executor >  local (4)
    [ba/08236d] sayHello (4) | 4 of 4 ✔
    Ciao world!

    Hello world!

    Bonjour world!

    Hola world!
    ```

이전 실행에서 `v1.3`을 명시적으로 고정했음에도 불구하고, 이번 실행은 기본 브랜치(`master`)로 돌아갑니다.
Nextflow는 사용한 각 리비전에 대해 별도의 로컬 작업 복사본을 유지합니다. 이것이 `nextflow info`에서 `>` 표시가 나타나는 이유입니다. 하지만 마지막으로 실행한 리비전을 기억하지는 않습니다.
파이프라인의 기본 브랜치 이름은 `nextflow info <pipeline>`을 실행하여 확인할 수 있습니다. `(default)`로 표시된 것이 기본 브랜치입니다.
재현성은 전적으로 사용자의 책임입니다. 이전 실행에서 고정한 리비전이 여전히 적용된다고 가정하지 말고, 필요할 때마다 항상 `-r`을 명시적으로 지정하세요.

### 핵심 정리

재현 가능한 실행을 위해 원격 파이프라인을 특정 버전, 브랜치, 또는 커밋으로 고정하는 방법과, 해당 고정이 그 실행에만 적용되며 이후 실행에는 적용되지 않는다는 점을 학습했습니다.

### 다음 단계

Nextflow 파이프라인 실행 및 관리의 기본 사항을 모두 학습했습니다.
이후 학습 방향은 [과정 요약](next_steps.md)을 참조하세요.

---

## 요약

이 파트에서는 다음 내용을 학습했습니다.

- GitHub 저장소에서 직접 파이프라인을 다운로드하지 않고 실행하기
- 재현 가능한 실행을 위해 원격 파이프라인을 특정 리비전으로 고정하기, 그리고 해당 고정이 그 실행에만 적용된다는 점 이해하기
