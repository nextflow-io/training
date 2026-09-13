# 시작하기

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } AI 지원 번역 - [자세히 알아보기 및 개선 제안](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

## 교육 환경 시작

GitHub Codespaces에서 제공하는 사전 구성된 환경을 사용하려면 아래의 "Open in GitHub Codespaces" 버튼을 클릭하세요. 다른 옵션은 [환경 설정](../envsetup/index.md)을 참조하세요.

환경이 로드되는 동안 계속 읽을 수 있도록, 교육 환경을 새 브라우저 탭이나 창에서 여는 것을 권장합니다(사용 중인 운영체제에 따라 마우스 오른쪽 클릭, ctrl-클릭 또는 cmd-클릭을 사용하세요).
과정을 진행하는 동안 이 안내 페이지를 함께 열어 두어야 합니다.

[![Open in GitHub Codespaces](https://github.com/codespaces/badge.svg)](https://codespaces.new/nextflow-io/training?quickstart=1&ref=master)

### 환경 기본 사항

이 교육 환경에는 교육 과정을 진행하는 데 필요한 모든 소프트웨어, 코드, 데이터가 포함되어 있으므로 별도로 설치할 필요가 없습니다.

Codespace는 VSCode 인터페이스로 구성되어 있으며, 파일 탐색기(File Explorer), 코드 편집기, 터미널 셸이 포함되어 있습니다.
과정 중 제공되는 모든 안내(예: '파일 열기', '코드 편집' 또는 '명령 실행')는 별도로 명시되지 않는 한 VSCode 인터페이스의 이 세 가지 구성 요소를 가리킵니다.

혼자 이 과정을 진행하는 경우, 자세한 내용은 [환경 기본 사항](../envsetup/01_setup.md)을 참고하세요.

### 버전 요구 사항

이 과정은 Nextflow 25.10.2 이상이 필요하며, v2 문법 분석기가 활성화되어 있어야 합니다(25.10 이상에서는 기본값입니다).
로컬 또는 커스텀 환경을 사용하는 경우, [여기](../info/nxf_versions.md)에 문서화된 올바른 설정을 사용하고 있는지 확인하세요.

## 작업 준비

Codespace가 실행되면, 본격적으로 시작하기 전에 두 가지를 먼저 진행합니다. 작업 디렉토리를 설정하고, 제공된 자료를 살펴보는 것입니다.

### 작업 디렉토리 설정

기본적으로 Codespace는 모든 교육 과정의 루트 디렉토리에서 열립니다.
이 과정에서는 `execution-config/` 디렉토리로 이동합니다:

```bash
cd execution-config/
```

그런 다음 VSCode가 이 디렉토리에 집중하도록 설정하여, 파일 탐색기 사이드바에 관련 파일만 표시되도록 합니다:

```bash
code .
```

!!! tip "팁"

    어떤 이유로든 이 디렉토리를 벗어난 경우(예: Codespace가 중지 상태가 된 경우), GitHub Codespaces 교육 환경 내에서 실행 중이라면 전체 경로를 사용하여 언제든지 돌아올 수 있습니다:

    ```bash
    cd /workspaces/training/execution-config
    ```

### 제공된 자료 살펴보기

왼쪽의 파일 탐색기를 사용하거나 `tree` 명령을 사용하여 과정 자료를 살펴볼 수 있습니다.
터미널에서 다음 명령을 실행하여 전체 구조를 확인합니다:

```bash
tree . -L 2
```

??? abstract "디렉토리 내용"

    ```console
    .
    ├── data
    │   └── greetings.csv
    ├── main.nf
    ├── modules
    │   ├── collectGreetings.nf
    │   ├── convertToUpper.nf
    │   ├── cowpy.nf
    │   └── sayHello.nf
    └── nextflow.config
    ```

**`main.nf`** 및 **`modules/`** 파일은 [Nextflow Run](../nextflow_run/index.md)에서 사용한 것과 동일한 다단계 파이프라인이며, **`nextflow.config`** 파일도 거기서 이미 확인한 것과 동일한 설정 파일입니다.
이 연습들을 진행하면서 두 파일 모두 확장하게 됩니다.

**`data/`** 디렉토리에는 파이프라인이 읽어 들이는 CSV 입력 파일이 포함되어 있습니다.

## 준비 확인 목록

시작할 준비가 되었나요?

- [ ] 이 과정의 목표와 사전 요구 사항을 이해했습니다
- [ ] 환경이 실행 중입니다
- [ ] 작업 디렉토리를 적절히 설정했습니다

모든 항목을 확인했다면 시작할 준비가 된 것입니다.

**[파트 1: 컴퓨팅 환경에 맞게 조정하기](./01_packaging_and_execution.md)로 계속하려면 이 페이지의 오른쪽 하단 모서리에 있는 화살표를 클릭하세요.**
