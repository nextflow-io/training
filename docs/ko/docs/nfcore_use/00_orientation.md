# 시작하기

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } AI 지원 번역 - [자세히 알아보기 및 개선 제안](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

## 교육 환경 시작

GitHub Codespaces에서 제공하는 사전 구성된 환경을 사용하려면 아래의 "Open in GitHub Codespaces" 버튼을 클릭하세요. 다른 옵션은 [환경 옵션](../envsetup/index.md)을 참조하세요.

환경이 로드되는 동안 계속 읽을 수 있도록, 교육 환경을 새 브라우저 탭이나 창에서 여는 것을 권장합니다(사용 중인 운영체제에 따라 마우스 오른쪽 클릭, ctrl-클릭 또는 cmd-클릭을 사용하세요). 과정을 진행하는 동안 이 안내 페이지를 함께 열어 두어야 합니다.

[![Open in GitHub Codespaces](https://github.com/codespaces/badge.svg)](https://codespaces.new/nextflow-io/training?quickstart=1&ref=master)

### 환경 기본 사항

이 교육 환경에는 교육 과정을 진행하는 데 필요한 모든 소프트웨어, 코드, 데이터가 포함되어 있으므로 별도로 설치할 필요가 없습니다.

Codespace는 VSCode 인터페이스로 구성되어 있으며, 파일 탐색기(File Explorer), 코드 편집기, 터미널 셸이 포함되어 있습니다. 과정 중에 제공되는 모든 안내(예: '파일 열기', '코드 편집' 또는 '명령 실행')는 별도로 명시되지 않는 한 VSCode 인터페이스의 이 세 가지 구성 요소를 가리킵니다.

혼자서 이 과정을 진행하는 경우, 자세한 내용은 [환경 기본 사항](../envsetup/01_setup.md)을 참조하세요.

### 버전 요구 사항

이 교육은 **v2 syntax parser**가 적용된 Nextflow 25.10.2 이상에서 작동하며, v2 parser는 Nextflow 26.04부터 기본값으로 설정됩니다. 교육 환경에서는 별도로 설정할 필요가 없습니다. v2 parser가 적용된 Nextflow 26.04.4가 실행됩니다. 로컬 또는 커스텀 환경을 사용하는 경우 [버전 참고 사항](../info/nxf_versions.md)을 확인하세요.

!!! warning "nf-core/demo는 Nextflow 25.10.4 이상이 필요합니다"

    파트 1에서 사용하는 `nf-core/demo` 파이프라인은 자체적으로 최소 Nextflow 버전(`>=25.10.4`)을 요구하며, 이는 교육의 일반 최소 버전인 25.10.2보다 더 엄격합니다.
    교육 환경은 이미 이 요구 사항을 충족합니다. 로컬 또는 커스텀 환경을 사용하는 경우 Nextflow 25.10.4 이상을 사용하고 있는지 확인하세요.

이 교육에는 추가로 **nf-core tools 4.0.2**가 필요합니다. 다른 버전의 nf-core 도구를 사용하는 경우 과정을 따라가는 데 어려움이 있을 수 있습니다.

`nf-core --version` 명령을 사용하여 환경에 설치된 버전을 확인할 수 있습니다.

!!! warning "v2 parser 호환성"

    많은 nf-core 파이프라인이 아직 v2 syntax parser를 지원하지 않습니다.
    이 과정에서 사용하는 파이프라인 외의 nf-core 파이프라인을 실행하다가 오류가 발생하면, `export NXF_SYNTAX_PARSER=v1`을 설정하여 v1 parser로 전환해야 할 수 있습니다.
    자세한 내용은 [버전 참고 사항](../info/nxf_versions.md)을 확인하세요.

## 작업 준비

Codespace가 실행되면, 교육을 시작하기 전에 두 가지를 먼저 해야 합니다. 이 과정의 작업 디렉토리를 설정하고, 제공된 자료를 살펴보는 것입니다.

### 작업 디렉토리 설정

기본적으로 Codespace는 모든 교육 과정의 루트 디렉토리에서 열리지만, 이 과정에서는 `nfcore-use/` 디렉토리에서 작업합니다.

터미널에서 다음 명령을 실행하여 디렉토리를 변경합니다.

```bash
cd nfcore-use/
```

!!! tip "팁"

    어떤 이유로든 이 디렉토리를 벗어난 경우(예: Codespace가 중지 상태가 된 경우), GitHub Codespaces 교육 환경 내에서 실행 중이라면 전체 경로를 사용하여 언제든지 돌아올 수 있습니다.

    ```bash
    cd /workspaces/training/nfcore-use
    ```

다음으로, 이 디렉토리의 내용을 살펴봅니다.

### 제공된 자료 살펴보기

교육 작업 공간 왼쪽의 파일 탐색기(File Explorer)를 사용하여 이 디렉토리의 내용을 살펴볼 수 있습니다. 또는 `tree` 명령을 사용할 수도 있습니다.

```bash
tree .
```

??? abstract "디렉토리 내용"

    ```console
    .
    ├── custom.config
    ├── laptop.config
    ├── malformed_samplesheet.csv
    └── my_params.yml
    ```

- **`laptop.config` 파일**은 프로덕션 규모의 파이프라인을 로컬에서 실행할 때 리소스 사용량을 제한하기 위해 섹션 4에서 사용할 설정 파일입니다. 그때까지는 무시해도 됩니다.
- **`my_params.yml`, `malformed_samplesheet.csv`, `custom.config` 파일**은 파트 2에서 파일로부터 매개변수 설정, 입력 유효성 검사, process 수준의 설정 재정의를 시연하는 데 사용됩니다. 이 파일들도 그때까지는 무시해도 됩니다.

## 준비 완료 체크리스트

시작할 준비가 되었나요?

- [ ] 이 과정의 목표와 사전 요구 사항을 이해했습니다
- [ ] 환경이 실행 중입니다
- [ ] nf-core tools 4.0.2를 사용하고 있습니다(`nf-core --version`으로 확인)
- [ ] 작업 디렉토리를 적절히 설정했습니다

모든 항목을 확인했다면 시작할 준비가 된 것입니다.

**파트 1로 계속하려면 이 페이지의 오른쪽 하단 모서리에 있는 화살표를 클릭하세요.**
