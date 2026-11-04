# 파트 2: 명령줄에서 파이프라인 실행

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } AI 지원 번역 - [자세히 알아보기 및 개선 제안](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

[파트 1](./01_run_with_seqera.md)에서는 Seqera 웹 인터페이스에서 nf-core/rnaseq를 실행했습니다.
이번에는 `tw` CLI를 사용하여 명령줄에서 동일한 작업을 수행하고, 워크스페이스에 새 파이프라인을 추가합니다.

---

## 1. 명령줄에서 파이프라인 실행

실행 보기에서 **Command line** 탭을 클릭합니다.
Platform이 사용자를 대신하여 구성하고 제출한 정확한 `nextflow run` 명령을 확인할 수 있습니다. 이는 Use nf-core 과정에서 수동으로 실행했던 것과 동일한 유형의 명령입니다.

Platform은 Nextflow를 대체하는 것이 아니라 이를 조율합니다.
웹 인터페이스를 통해 할 수 있는 모든 작업은 Platform API와 상호작용하는 명령줄 도구인 `tw` CLI를 사용하여 터미널에서도 수행할 수 있습니다.
이는 스크립트나 CI/CD 파이프라인에서 실행을 자동화할 때 유용합니다.

이전 과정에서 사용했던 동일한 Codespace에서 지금 바로 진행합니다.

### 1.1. tw CLI 설치

Codespace 터미널에서 다음 명령을 실행하여 `tw` 바이너리를 다운로드하고 설치합니다:

```bash
curl -fsSL https://github.com/seqeralabs/tower-cli/releases/latest/download/tw-linux-x86_64 -o tw
chmod +x tw
sudo mv tw /usr/local/bin/
```

설치를 확인합니다:

```bash
tw --version
```

??? success "명령 출력"

    ```console
    Tower CLI version 0.40.0 (build 26db579)
    ```

`tw` CLI가 설치되었으며 설정할 준비가 되었습니다.

### 1.2. 액세스 토큰 발급

`tw` CLI는 개인 액세스 토큰을 사용하여 Seqera에 인증합니다.

1. Seqera 웹 인터페이스에서 오른쪽 상단의 아바타를 클릭하고 **Your tokens**를 선택합니다.
2. **Add token**을 클릭하고 이름(예: `training`)을 입력한 후 **Add**를 클릭합니다.
3. 토큰 값을 복사합니다. 토큰은 한 번만 표시됩니다.
   바로 저장해 두지 않으면 새로 생성해야 합니다.

### 1.3. CLI 설정

편의를 위해 방금 생성한 액세스 토큰과 워크스페이스 식별자를 포함하는 설정 파일을 구성합니다.

편집기에서 이 디렉토리의 `.seqera_config` 파일을 열고 두 변수를 설정합니다:

- **`TOWER_ACCESS_TOKEN`**: 섹션 1.2에서 생성한 토큰
- **`TOWER_WORKSPACE_ID`**: 워크스페이스의 숫자 ID (섹션 1.4에서 실행하는 `tw workspaces list`의 `ID` 열)

값을 입력한 후 설정을 로드합니다:

```bash
source .seqera_config
```

연결을 확인합니다:

```bash
tw info
```

??? success "명령 출력"

    ```console
        Details
    -------------------------+-----------------------------
     Tower API endpoint      | https://api.cloud.seqera.io
     Tower API version       | 1.150.0
     Tower version           | 26.1.0-cycle54
     CLI version             | 0.30.0 (fde9dec)
     CLI minimum API version | 1.148.0
     Authenticated user      | <your-name>

    System health status
    ---------------------------------------+----
     Remote API server connection check    | OK
     Tower API version check               | OK
     Authentication API credential's token | OK
    ```

`tw` CLI가 인증되었으며 Seqera 계정에 연결되었습니다.
Codespace 세션을 시작할 때마다 `source .seqera_config`를 실행하여 설정을 다시 로드합니다.

!!! tip "팁"

    워크스페이스에 기본 컴퓨팅 환경이 설정되어 있지 않은 경우, 설정 파일에 `export TOWER_COMPUTE_ENV=<compute-env-name>`을 추가하여 기본값을 지정할 수 있습니다.
    모든 설정 값은 명령줄에서 플래그를 명시적으로 전달하여 재정의할 수 있습니다(예: `--compute-env other-env`).
    전체 옵션 및 환경 변수 목록은 [tw CLI 참조 문서](https://docs.seqera.io/platform/latest/cli/reference)를 확인하세요.

### 1.4. CLI에서 워크스페이스 탐색

접근 가능한 워크스페이스 목록을 확인합니다:

```bash
tw workspaces list
```

??? success "명령 출력"

    ```console
    Available workspaces:
    ID              | Name            | Full Name                | Visibility
    --------------- | --------------- | ------------------------ | ----------
    <workspace-id>  | my-workspace    | my-org/my-workspace      | PRIVATE
    ```

방금 실행한 nf-core/rnaseq 실행을 포함하여 워크스페이스의 실행 목록을 확인합니다:

```bash
tw runs list
```

??? success "명령 출력"

    ```console
    Pipeline runs for my-org/my-workspace workspace:
    ID        | Status   | Name          | Pipeline               | Run name
    --------- | -------- | ------------- | ---------------------- | --------
    <run-id>  | RUNNING  | ...           | nf-core/rnaseq         | happy_curie
    ```

웹 인터페이스에서 모니터링 중인 동일한 실행이 여기에도 표시됩니다.

!!! note "참고"

    `.seqera_config`에 `TOWER_WORKSPACE_ID`가 설정되어 있으므로 모든 `tw` 명령에서 `--workspace`를 생략할 수 있습니다.
    설정 파일 없이 사용하는 경우 다음과 같이 명시적으로 전달합니다:

    ```bash
    tw runs list --workspace <org>/<workspace>
    ```

웹 인터페이스에서 볼 수 있는 모든 것은 CLI에서도 접근할 수 있습니다.

### 1.5. CLI에서 nf-core/rnaseq 실행

[파트 1](./01_run_with_seqera.md)에서 워크스페이스에 추가한 파이프라인은 CLI에서 이름으로 사용할 수 있습니다.
`test` 프로파일로 실행합니다:

```bash
tw launch nf-core-rnaseq --profile test
```

??? success "명령 출력"

    ```console
    Launching pipeline nf-core-rnaseq
    Run name: focused_einstein
    https://cloud.seqera.io/orgs/my-org/workspaces/my-workspace/watch/<run-id>
    ```

브라우저에서 링크를 열고 **Runs** 패널에 실행이 표시되는지 확인합니다.

실행 중인 것을 확인했다면, CLI와 웹 인터페이스가 동일한 워크스페이스를 보는 두 가지 방식임을 확인한 것입니다.

!!! note "참고"

    파이프라인을 워크스페이스에 추가하지 않고도 `tw launch`에 GitHub URL을 직접 전달할 수 있습니다.
    그러나 실행 전에 파이프라인을 명시적으로 추가하는 것이 일반적으로 더 좋습니다. 이렇게 하면 향후 실행을 위한 파이프라인 설정이 저장되고, 이름으로 사용할 수 있으며, Launchpad에서 모든 워크스페이스 구성원이 볼 수 있습니다.

    `tw`를 사용하여 명령줄에서 직접 워크스페이스에 파이프라인을 추가할 수도 있습니다.
    다음 섹션에서는 nf-core/demo 파이프라인을 사용하여 이 방법을 설명합니다.

### 핵심 정리

`tw` CLI를 인증하고, 워크스페이스를 확인하며, 터미널에서 저장된 파이프라인을 실행하는 방법을 학습했습니다.

### 다음 단계

명령줄에서 워크스페이스에 새 파이프라인을 추가하고 실행합니다.

---

## 2. 새 파이프라인 추가 및 실행

GitHub에 있는 모든 Nextflow 파이프라인은 루트에 `main.nf` 진입점과 `nextflow.config`가 있는 한 `tw pipelines add`를 사용하여 워크스페이스에 추가할 수 있습니다.
nf-core/demo는 연습하기에 좋은 예입니다. Use nf-core 과정에서 이미 실행해 보았으므로 무엇을 하는지, 무엇을 기대할 수 있는지 알고 있습니다.

### 2.1. 워크스페이스에 nf-core/demo 추가

다음 명령을 실행하여 워크스페이스에 파이프라인을 등록합니다:

```bash
tw pipelines add \
    --name nf-core-demo \
    https://github.com/nf-core/demo
```

??? success "명령 출력"

    ```console
    New pipeline 'nf-core-demo' added at my-org/my-workspace workspace
    ```

파이프라인이 등록되었으며 Launchpad에 표시됩니다.

### 2.2. Launchpad에 표시되는지 확인

워크스페이스의 파이프라인 목록을 확인하여 추가되었는지 검증합니다:

```bash
tw pipelines list
```

??? success "명령 출력"

    ```console
    Pipelines at my-org/my-workspace:
    ID   | Name            | Repository
    ---- | --------------- | -----------------------------------------
    ...  | nf-core-demo    | https://github.com/nf-core/demo
    ...  | nf-core-rnaseq  | ...
    ```

브라우저에서 워크스페이스를 열고 **Launchpad**를 클릭하여 nf-core/demo가 nf-core/rnaseq와 함께 표시되는지 확인합니다.

!!! tip "팁"

    웹 인터페이스를 통해서도 파이프라인을 추가할 수 있습니다. 왼쪽 사이드바에서 **Launchpad**를 클릭한 후 **Add pipeline**을 클릭하고 양식을 작성합니다.

nf-core/demo 항목의 **Launch** 버튼을 클릭하여 실행 양식을 엽니다.
`input`과 `outdir` 매개변수가 빨간색으로 강조 표시된 것을 확인할 수 있습니다. 이는 기본값이 없는 필수 항목으로, `tw pipelines add`가 매개변수를 사전 설정하지 않고 파이프라인 소스만 등록하기 때문입니다.
다음 두 섹션에서는 이 값을 제공하는 방법을 설명합니다. 먼저 웹 양식을 통해, 그 다음 명령줄에서 진행합니다.

### 2.3. 웹 인터페이스에서 nf-core/demo 실행

실행 양식이 열린 상태에서 두 개의 필수 매개변수를 입력합니다.

`input`에는 nf-core/demo 테스트 프로파일의 테스트 샘플시트 URL을 입력합니다.
Use nf-core 과정에서 살펴본 파이프라인 저장소의 `conf/test.config`에서 확인할 수 있습니다:

```
https://raw.githubusercontent.com/nf-core/test-datasets/viralrecon/samplesheet/samplesheet_test_illumina_amplicon.csv
```

`outdir`에는 파이프라인이 결과를 저장할 클라우드 스토리지 경로를 입력합니다.
워크스페이스에 설정된 버킷을 사용하고, 실행을 체계적으로 관리하기 위해 하위 디렉토리를 지정합니다:

```
s3://my-bucket/demo-results
```

두 항목을 모두 입력한 후 파란색 **Launch** 버튼을 클릭합니다.

실행이 **Runs** 패널에 표시되며 테스트 데이터셋에서 몇 분 내에 완료됩니다.
실행을 클릭하여 작업 테이블과 실행 보고서를 확인합니다.

### 2.4. CLI에서 nf-core/demo 실행

`nextflow run`과 달리 `tw launch` 명령은 `--input`이나 `--outdir`과 같은 개별 매개변수 플래그를 허용하지 않습니다.
매개변수는 `--params-file`로 전달되는 YAML 또는 JSON 형식의 파일을 통해 제공해야 합니다.
이는 재현성을 높이기 위한 방식입니다. 저장된 매개변수 파일은 실행에 사용된 정확한 값을 기록하여 실행 설정을 쉽게 반복하거나 공유할 수 있게 합니다.

작업 디렉토리에 매개변수 파일을 생성합니다:

```bash
touch params.yaml
```

편집기에서 파일을 열고 출력 경로를 추가합니다:

```yaml title="params.yaml"
outdir: "s3://my-bucket/demo-results-cli"
```

이제 `test` 프로파일(`input` 샘플시트 제공)과 매개변수 파일(`outdir` 제공)을 사용하여 파이프라인을 실행합니다:

```bash
tw launch nf-core-demo --profile test --params-file params.yaml
```

??? success "명령 출력"

    ```console
    Launching pipeline nf-core-demo
    Run name: focused_feynman
    https://cloud.seqera.io/orgs/my-org/workspaces/my-workspace/watch/<run-id>
    ```

링크를 열어 **Runs** 패널에 실행이 표시되는지 확인합니다.

!!! tip "팁"

    초기 설정 단계에서 매개변수 파일을 포함하여 기본값을 설정하고, 웹 양식에서 설정한 내용과 일치하도록 추가 속성을 지정할 수 있습니다:

    ```bash
    tw pipelines add \
      --name nf-core-demo-gg \
      --description "Demo pipeline with defaults for testing" \
      --params-file params.yaml \
      --profile test \
      https://github.com/nf-core/demo
    ```

### 핵심 정리

GitHub에 호스팅된 Nextflow 파이프라인을 워크스페이스에 추가하고 실행하는 방법을 학습했습니다. 웹 인터페이스에서 매개변수를 수동으로 입력하는 방법과 `tw` CLI에서 프로파일과 매개변수 파일을 조합하는 방법 모두 확인했습니다.

---

## 요약

이번 파트에서 학습한 내용:

- `tw` CLI를 인증하고 터미널에서 저장된 파이프라인 실행
- CLI를 사용하여 GitHub에서 새 파이프라인을 추가하고 Launchpad에 표시되는지 확인
- Seqera 웹 인터페이스에서 필수 매개변수를 수동으로 입력하여 파이프라인 실행
- Nextflow 프로파일과 매개변수 파일을 사용하여 CLI에서 파이프라인 실행
