# 파트 3: 워크플로우 실행 관리

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } AI 지원 번역 - [자세히 알아보기 및 개선 제안](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

파이프라인을 반복적으로 실행하다 보면 실행 기록과 오래된 `work/` 디렉토리가 쌓이게 됩니다.
[파트 1](./01_run_nextflow.md#23-use-resume-to-skip-completed-work)에서 이미 완료된 작업을 건너뛰기 위해 `-resume`을 사용했습니다.
여기서는 실행 보고서를 생성하는 방법, [`nextflow log`](https://nextflow.io/docs/latest/reference/cli.html#log)로 과거 실행 기록을 확인하는 방법, 그리고 더 이상 필요하지 않은 오래된 work 디렉토리를 [`nextflow clean`](https://nextflow.io/docs/latest/reference/cli.html#clean)으로 삭제하는 방법을 학습합니다.

---

## 1. 파이프라인 보고서 생성

Nextflow는 실행에 대한 다양한 종류의 보고서를 생성할 수 있으며, 각각 고유한 `-with-*` 플래그로 추가합니다. 실행 보고서(`-with-report`), 실행 타임라인(`-with-timeline`), 작업 추적 파일(`-with-trace`), 워크플로우 다이어그램(`-with-dag`)이 있습니다.
여기서는 처음 두 가지를 생성합니다. 나머지는 Nextflow 참조 문서의 [Execution reports](https://nextflow.io/docs/latest/reports.html)를 확인하세요.

### 1.1. 실행 보고서 생성

`nextflow run` 명령에 `-with-report`를 추가하면 파이프라인 완료 후 HTML 보고서가 생성됩니다.

```bash
nextflow run main.nf --input data/greetings.csv --character turkey -with-report
```

??? success "명령 출력"

    ```console hl_lines="6 7 8 9"
    N E X T F L O W   ~  version 26.04.4

    Launching `main.nf` [intergalactic_dalembert] revision: ce74f81996

    executor >  local (8)
    [34/23f10f] sayHello (2)       | 3 of 3 ✔
    [af/cab69d] convertToUpper (3) | 3 of 3 ✔
    [9e/d73afb] collectGreetings   | 1 of 1 ✔
    [3c/392db0] cowpy              | 1 of 1 ✔

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

Nextflow는 현재 작업 디렉토리에 `report-<timestamp>.html`이라는 이름의 파일로 보고서를 저장합니다.
브라우저에서 열면 실행 요약, 각 작업의 상태 및 실행 시간이 포함된 표, 그리고 프로세스별로 분류된 리소스 사용량 차트를 확인할 수 있습니다.

**Tasks** 탭에는 파이프라인이 실행한 모든 작업이 프로세스 이름, 상태, 리소스 사용량과 함께 나열됩니다.

![실행 보고서 작업 표](img/execution_report_tasks.png)

이 보고서는 파이프라인이 예상보다 오래 걸리거나 작업이 실패했을 때 특히 유용합니다. 작업 표에서 시간이 어디에 소요되었는지, 어떤 작업이 성공하거나 실패했는지 정확히 확인할 수 있습니다.

### 1.2. 실행 타임라인 생성

실행 시 `-with-timeline`을 추가하면 각 작업이 언제 실행되었는지 간트 차트 형식으로 확인할 수 있습니다.

```bash
nextflow run main.nf --input data/greetings.csv --character turkey -with-timeline
```

??? success "명령 출력"

    ```console hl_lines="6 7 8 9"
    N E X T F L O W   ~  version 26.04.4

    Launching `main.nf` [jolly_noyce] revision: ce74f81996

    executor >  local (8)
    [ad/e92ef3] sayHello (3)       | 3 of 3 ✔
    [2a/df8a8d] convertToUpper (2) | 3 of 3 ✔
    [be/7fb72a] collectGreetings   | 1 of 1 ✔
    [63/dc9bd6] cowpy              | 1 of 1 ✔

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

Nextflow는 타임라인을 `timeline-<timestamp>.html`이라는 이름의 파일로 저장합니다.
브라우저에서 열면 각 작업이 언제 실행되었고 얼마나 걸렸는지를 막대 형태로 확인할 수 있습니다.

![실행 타임라인](img/execution_timeline.png)

타임라인을 통해 [파트 1](./01_run_nextflow.md#31-run-the-workflow)에서 살펴본 팬아웃-팬인 구조를 한눈에 확인할 수 있습니다. 세 개의 `sayHello` 작업이 병렬로 실행되고, 이어서 세 개의 `convertToUpper` 작업이 실행된 후, `collectGreetings`와 `cowpy`가 각각 이전 작업에 의존하기 때문에 순차적으로 실행됩니다.

### 핵심 정리

`-with-report`로 HTML 실행 보고서를 생성하고, `-with-timeline`으로 실행 타임라인을 생성하는 방법을 학습했습니다. 또한 Nextflow가 지원하는 다른 보고서 유형을 어디서 확인할 수 있는지도 알게 되었습니다.

### 다음 단계

과거 실행 기록을 확인하는 방법을 학습합니다.

---

## 2. 과거 실행 로그 확인

파이프라인을 개발하거나 운영 환경에서 실행할 때, 과거 실행에 대한 정보를 조회해야 하는 경우가 생깁니다.

### 2.1. 히스토리 파일

Nextflow 워크플로우를 실행할 때마다 현재 작업 디렉토리의 `.nextflow`라는 숨김 디렉토리 아래에 있는 `history`라는 로그 파일에 한 줄씩 기록됩니다.

??? abstract "파일 내용"

    ```txt title=".nextflow/history" linenums="1"
    2026-09-13 01:47:58	2.8s	determined_ramanujan	OK	ce74f81996b8ce999853fdbb25ba7969	03078950-5011-414f-992b-8fe1bd219ac2	nextflow run main.nf --input data/greetings.csv --character turkey
    2026-09-13 01:48:03	3s	prickly_cuvier	OK	ce74f81996b8ce999853fdbb25ba7969	5a9e8bf7-a7db-4c2c-b1d8-3bd4c7266528	nextflow run main.nf --input data/greetings.csv --character tux
    2026-09-13 01:48:09	2.3s	elegant_panini	OK	ce74f81996b8ce999853fdbb25ba7969	517ccc29-9db7-4d8d-8ea0-c8db6e9ac653	nextflow run main.nf --input data/greetings.csv --character stegosaurus
    2026-09-13 01:48:14	1.5s	deadly_lamport	OK	ce74f81996b8ce999853fdbb25ba7969	517ccc29-9db7-4d8d-8ea0-c8db6e9ac653	nextflow run main.nf --input data/greetings.csv --character stegosaurus -resume
    ```

각 줄에는 타임스탬프, 실행 시간, 실행 이름, 상태, 리비전 ID, 세션 ID, 그리고 해당 디렉토리에서 실행된 전체 명령이 포함됩니다.

마지막 두 줄을 살펴보면, 동일한 명령을 두 번 실행한 것입니다(하나는 일반 실행, 하나는 `-resume` 사용). 두 실행은 동일한 세션 ID를 공유합니다.
세션 ID는 완전히 새로운 실행을 시작할 때만 변경됩니다. `-resume`을 사용하면 세션 ID가 유지되며, Nextflow는 이를 통해 어떤 cache를 재사용할지 결정합니다.

### 2.2. `nextflow log`로 더 편리하게 확인하기

히스토리 파일을 직접 읽을 수도 있지만, `nextflow log`를 사용하면 헤더와 함께 동일한 정보를 보기 좋게 출력합니다.

```bash
nextflow log
```

??? success "명령 출력"

    ```console linenums="1"
    TIMESTAMP          	DURATION	RUN NAME            	STATUS	REVISION ID	SESSION ID                          	COMMAND
    2026-09-13 01:47:58	2.8s    	determined_ramanujan	OK    	ce74f81996 	03078950-5011-414f-992b-8fe1bd219ac2	nextflow run main.nf --input data/greetings.csv --character turkey
    2026-09-13 01:48:03	3s      	prickly_cuvier      	OK    	ce74f81996 	5a9e8bf7-a7db-4c2c-b1d8-3bd4c7266528	nextflow run main.nf --input data/greetings.csv --character tux
    2026-09-13 01:48:09	2.3s    	elegant_panini      	OK    	ce74f81996 	517ccc29-9db7-4d8d-8ea0-c8db6e9ac653	nextflow run main.nf --input data/greetings.csv --character stegosaurus
    2026-09-13 01:48:14	1.5s    	deadly_lamport      	OK    	ce74f81996 	517ccc29-9db7-4d8d-8ea0-c8db6e9ac653	nextflow run main.nf --input data/greetings.csv --character stegosaurus -resume
    ```

Nextflow는 `-resume`에 사용되는 caching 정보를 세션 ID를 키로 하여 `.nextflow/cache` 아래에 저장합니다.
따라서 과거 실행을 조사하거나 정리해야 할 때, 올바른 실행 이름이나 세션 ID를 여기서 먼저 확인하는 것이 첫 번째 단계입니다.

### 핵심 정리

Nextflow가 과거 실행 기록을 어디에 저장하는지, 그리고 `nextflow log`로 이를 확인하는 방법을 학습했습니다.

### 다음 단계

더 이상 필요하지 않은 오래된 work 디렉토리를 삭제하는 방법을 학습합니다.

---

## 3. 오래된 work 디렉토리 삭제

모든 실행은 결과물을 `results/`에 복사한 후에도 `work/` 아래에 작업 디렉토리를 남겨둡니다.
개발 중에 파이프라인을 여러 번 실행하다 보면 이 하위 디렉토리들이 쌓이게 됩니다. Nextflow는 더 이상 필요하지 않은 디렉토리를 삭제할 수 있도록 `nextflow clean`을 제공합니다.

### 3.1. 삭제 기준 결정

`nextflow clean`은 삭제할 항목을 선택하는 여러 방법을 지원합니다. 전체 목록은 [참조 문서](https://www.nextflow.io/docs/latest/reference/cli.html#clean)를 확인하세요.
여기서는 특정 실행 이름을 기준으로 그 이전의 모든 실행을 삭제합니다.

`nextflow log`를 사용하여 유지하고 싶은 가장 최근 실행을 확인합니다. [2.2의 예시](#22-use-nextflow-log-for-a-friendlier-view)에서는 `-resume` 실행 직전의 마지막 일반 실행인 `elegant_panini`가 해당됩니다.
실행 이름은 콘솔의 `Launching (...)` 줄이나 `nextflow log`의 `RUN NAME` 열에 표시되는 자동 생성된 두 단어 조합 문자열입니다.

### 3.2. 드라이 런 실행

실제로 삭제하기 전에 `-n`을 추가하여 어떤 항목이 삭제될지 먼저 확인합니다.

```bash
nextflow clean -before elegant_panini -n
```

??? success "명령 출력"

    ```console
    Would remove /workspaces/training/nextflow-run/work/e5/9ec3fe24deb5d1ffa478dff5e5e1d4
    Would remove /workspaces/training/nextflow-run/work/75/36b1ece86b213c331852d0a3dc659b
    Would remove /workspaces/training/nextflow-run/work/be/8157808f28699259b5ff37e4bf9f72
    Would remove /workspaces/training/nextflow-run/work/f3/d828b7ba6b315eb08e0f96168119f5
    Would remove /workspaces/training/nextflow-run/work/94/7a2d7ef67d81b6e66430f300816698
    Would remove /workspaces/training/nextflow-run/work/eb/e37c5c8638c45916de3d9794d6d22f
    Would remove /workspaces/training/nextflow-run/work/3f/f56a705a6370f6c615aff22d3f240a
    Would remove /workspaces/training/nextflow-run/work/25/4922fa35ba2c6087777ddce20e9f2f
    Would remove /workspaces/training/nextflow-run/work/a4/2de7bcc2f8fc974261f7280c1480f7
    Would remove /workspaces/training/nextflow-run/work/48/8c49bded446d4fcb78c8125318c5e3
    Would remove /workspaces/training/nextflow-run/work/8e/c9834a9bdd7c21de2991d959201a9d
    Would remove /workspaces/training/nextflow-run/work/d5/a9055dc8c817c6c1f436494685955c
    Would remove /workspaces/training/nextflow-run/work/27/9473df9e55e35c8869d9571b480926
    Would remove /workspaces/training/nextflow-run/work/d3/ebc066e1d844232b9f5ee0a7d87464
    Would remove /workspaces/training/nextflow-run/work/91/17ef1815ed06b8f177d0135d8ef0dc
    Would remove /workspaces/training/nextflow-run/work/12/0443f9fcff27a552d4bc73aaedf24a
    ```

총 16개의 작업 디렉토리입니다. `turkey` 실행의 8개 작업과 `tux` 실행의 8개 작업으로, 4개의 프로세스로 구성된 이 파이프라인을 두 번 완전히 실행했을 때 예상되는 수와 정확히 일치합니다.
`elegant_panini` 실행 자체와 `-resume` 실행에서 재사용된 캐시된 작업들은 그대로 유지됩니다.

출력에 표시되는 디렉토리 이름은 다를 수 있으며, 줄 수는 실행 횟수에 따라 달라집니다. 아무 줄도 표시되지 않는다면 실행 이름이 로그의 항목과 일치하지 않거나, 해당 실행 이전에 삭제할 항목이 없는 것입니다.

### 3.3. 실제 삭제 진행

드라이 런 결과가 올바르다면, 동일한 명령에서 `-n` 대신 `-f`를 사용하여 다시 실행합니다.

```bash
nextflow clean -before elegant_panini -f
```

??? success "명령 출력"

    ```console
    Removed /workspaces/training/nextflow-run/work/e5/9ec3fe24deb5d1ffa478dff5e5e1d4
    Removed /workspaces/training/nextflow-run/work/75/36b1ece86b213c331852d0a3dc659b
    Removed /workspaces/training/nextflow-run/work/be/8157808f28699259b5ff37e4bf9f72
    Removed /workspaces/training/nextflow-run/work/f3/d828b7ba6b315eb08e0f96168119f5
    Removed /workspaces/training/nextflow-run/work/94/7a2d7ef67d81b6e66430f300816698
    Removed /workspaces/training/nextflow-run/work/eb/e37c5c8638c45916de3d9794d6d22f
    Removed /workspaces/training/nextflow-run/work/3f/f56a705a6370f6c615aff22d3f240a
    Removed /workspaces/training/nextflow-run/work/25/4922fa35ba2c6087777ddce20e9f2f
    Removed /workspaces/training/nextflow-run/work/a4/2de7bcc2f8fc974261f7280c1480f7
    Removed /workspaces/training/nextflow-run/work/48/8c49bded446d4fcb78c8125318c5e3
    Removed /workspaces/training/nextflow-run/work/8e/c9834a9bdd7c21de2991d959201a9d
    Removed /workspaces/training/nextflow-run/work/d5/a9055dc8c817c6c1f436494685955c
    Removed /workspaces/training/nextflow-run/work/27/9473df9e55e35c8869d9571b480926
    Removed /workspaces/training/nextflow-run/work/d3/ebc066e1d844232b9f5ee0a7d87464
    Removed /workspaces/training/nextflow-run/work/91/17ef1815ed06b8f177d0135d8ef0dc
    Removed /workspaces/training/nextflow-run/work/12/0443f9fcff27a552d4bc73aaedf24a
    ```

`nextflow clean`은 작업 디렉토리의 내용을 삭제하지만, `e5/`와 같은 두 글자 상위 디렉토리는 그대로 남겨둡니다.

!!! warning "경고"

    과거 실행의 work 디렉토리를 삭제하면 Nextflow의 cache에서 해당 항목이 제거되고, 그곳에만 저장된 출력물도 삭제됩니다.
    이렇게 되면 해당 프로세스를 다시 실행하지 않고는 실행을 resume할 수 없게 됩니다. 따라서 resume이 필요하지 않다고 확신하는 실행만 정리하세요.
    이것이 바로 `work/` 디렉토리나 `symlink` publish 모드에 의존하지 않고, 중요한 결과물은 `mode 'copy'`를 사용하여 `results/`에 게시하는 것이 좋은 이유이기도 합니다.

### 핵심 정리

`nextflow clean`으로 오래된 work 디렉토리를 삭제하는 방법과, 삭제 시 해당 실행에서의 resume 기능을 포기하게 된다는 점을 학습했습니다.

### 다음 단계

[파트 4](./04_remote_repositories.md)에서 GitHub와 같은 원격 저장소에서 직접 파이프라인을 실행하는 방법을 학습합니다.

---

## 요약

이번 파트에서 학습한 내용:

- `-with-report`로 HTML 실행 보고서 생성하기, `-with-timeline`으로 실행 타임라인 생성하기
- `nextflow log`로 과거 실행 기록 확인하기
- `nextflow clean`으로 오래된 work 디렉토리 삭제하기, 그리고 이에 따른 resume 기능 포기에 대한 이해
