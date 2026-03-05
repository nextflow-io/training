# 시작하기

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } AI 지원 번역 - [자세히 알아보기 및 개선 제안](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

## 교육 환경 시작하기

GitHub Codespaces에서 제공하는 사전 구축된 환경을 사용하려면 아래의 "Open in GitHub Codespaces" 버튼을 클릭하십시오. 다른 옵션은 [환경 옵션](../envsetup/index.md)을 참조하십시오.

환경이 로드되는 동안 계속 읽을 수 있도록 새 브라우저 탭 또는 창에서 교육 환경을 여는 것을 권장합니다(장비에 따라 우클릭, ctrl-클릭 또는 cmd-클릭 사용).
과정을 진행하려면 이 지침을 계속 열어두어야 합니다.

[![Open in GitHub Codespaces](https://github.com/codespaces/badge.svg)](https://codespaces.new/nextflow-io/training?quickstart=1&ref=master)

### 환경 기본 사항

이 교육 환경에는 교육 과정을 진행하는 데 필요한 모든 소프트웨어, 코드 및 데이터가 포함되어 있으므로 직접 설치할 필요가 없습니다.

codespace는 VSCode 인터페이스로 설정되어 있으며, 파일 시스템 탐색기, 코드 편집기 및 터미널 셸이 포함되어 있습니다.
과정 중에 제공되는 모든 지침(예: '파일 열기', '코드 편집' 또는 '이 명령 실행')은 달리 명시되지 않는 한 VSCode 인터페이스의 이 세 부분을 참조합니다.

이 과정을 혼자 진행하는 경우 자세한 내용은 [환경 기본 사항](../envsetup/01_setup.md)을 숙지하십시오.

### 버전 요구 사항

이 교육은 **Nextflow 25.10.2** 또는 이후 버전에서 **v2 syntax parser가 비활성화된 상태**로 설계되었습니다.

#### 제공하는 교육 환경을 사용하는 경우:

더 진행하기 전에 다음 명령을 반드시 실행해야 합니다:

```bash
export NXF_SYNTAX_PARSER=v1
```

#### 로컬 또는 사용자 지정 환경을 사용하는 경우:

[여기](../info/nxf_versions.md)에 문서화된 올바른 설정을 사용하고 있는지 확인하십시오.

이 교육에는 추가로 **nf-core tools 3.4.1**이 필요합니다.
다른 버전의 nf-core 도구를 사용하는 경우 따라가기 어려울 수 있습니다.

`nf-core --version` 명령을 사용하여 환경에 설치된 버전을 확인할 수 있습니다.

## 작업 준비하기

codespace가 실행되면 교육에 들어가기 전에 두 가지 작업을 수행해야 합니다: 이 특정 과정의 작업 디렉토리를 설정하고 제공된 자료를 살펴보는 것입니다.

### 작업 디렉토리 설정하기

기본적으로 codespace는 모든 교육 과정의 루트에 작업 디렉토리가 설정된 상태로 열리지만, 이 과정에서는 `hello-nf-core/` 디렉토리에서 작업할 것입니다.

터미널에서 다음 명령을 실행하여 지금 디렉토리를 변경하십시오:

```bash
cd hello-nf-core/
```

!!! tip "팁"

    어떤 이유로든 이 디렉토리를 벗어나는 경우(예: codespace가 중지 모드로 전환됨) Github Codespaces 교육 환경 내에서 실행하는 경우 전체 경로를 사용하여 언제든지 돌아갈 수 있습니다:

    ```bash
    cd /workspaces/training/hello-nf-core
    ```

이제 이 디렉토리의 내용을 살펴보겠습니다.

### 제공된 자료 탐색하기

교육 작업 공간의 왼쪽에 있는 파일 탐색기를 사용하여 이 디렉토리의 내용을 탐색할 수 있습니다.
또는 `tree` 명령을 사용할 수 있습니다.

과정 전반에 걸쳐 `tree`의 출력을 사용하여 디렉토리 구조와 내용을 읽기 쉬운 형태로 표현하며, 때로는 명확성을 위해 약간 수정합니다.

여기서는 두 번째 레벨까지 목차를 생성합니다:

```bash
tree . -L 2
```

??? abstract "디렉토리 내용"

    ```console
    .
    ├── greetings.csv
    ├── original-hello
    │   ├── hello.nf
    │   ├── modules
    │   └── nextflow.config
    └── solutions
        ├── composable-hello
        ├── core-hello-part2
        ├── core-hello-part3
        ├── core-hello-part4
        ├── core-hello-part5
        └── core-hello-start
    ```

색상이 있는 상자를 클릭하여 섹션을 확장하고 내용을 확인하십시오.
예상되는 명령 출력을 간결한 방식으로 포함하기 위해 이와 같은 접을 수 있는 섹션을 사용합니다.

- **`greetings.csv` 파일**은 테스트 목적으로 사용하는 최소한의 열 데이터가 포함된 CSV입니다.

- **`original-hello` 디렉토리**에는 전체 Hello Nextflow 교육 시리즈를 진행하여 생성된 소스 코드의 사본이 포함되어 있습니다(Docker 활성화 상태).

- **`solutions` 디렉토리**에는 과정의 각 단계에서 생성되는 완성된 workflow 스크립트가 포함되어 있습니다.
  작업을 확인하고 문제를 해결하는 데 사용할 참조 자료로 제공됩니다.

## 준비 상태 체크리스트

시작할 준비가 되었다고 생각하십니까?

- [ ] 이 과정의 목표와 전제 조건을 이해했습니다
- [ ] 환경이 실행 중입니다
- [ ] syntax parser가 **v1**로 설정되어 있는지 확인했습니다
- [ ] 작업 디렉토리를 적절하게 설정했습니다

모든 항목을 체크할 수 있다면 준비가 완료된 것입니다.

**Part 1로 계속하려면 이 페이지의 오른쪽 하단 모서리에 있는 화살표를 클릭하십시오.**
