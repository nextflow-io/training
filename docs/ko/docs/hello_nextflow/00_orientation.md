# 시작하기

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } AI 지원 번역 - [자세히 알아보기 및 개선 제안](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

<div class="video-wrapper">
  <iframe width="560" height="315" src="https://www.youtube.com/embed/PIjOdFaYwWA?si=nyLvwhrSB2m1NPc5&amp;list=PLPZ8WHdZGxmWKozQuzr27jyMGqp9kElVK&amp;cc_load_policy=1&amp;cc_lang_pref=ko" title="YouTube video player" frameborder="0" allow="accelerometer; autoplay; clipboard-write; encrypted-media; gyroscope; picture-in-picture; web-share" referrerpolicy="strict-origin-when-cross-origin" allowfullscreen></iframe>
</div>

/// caption
:fontawesome-brands-youtube:{ .youtube } Nextflow YouTube 채널에서 [전체 재생 목록](https://youtube.com/playlist?list=PLPZ8WHdZGxmWKozQuzr27jyMGqp9kElVK&si=eF7cLR62goy-lc6n)을 확인하세요.

:green_book: 비디오 스크립트는 [여기](./transcripts/00_orientation.md)에서 확인할 수 있습니다.
///

!!! tip

    YouTube 비디오에는 몇 가지 강력한 기능이 있습니다!

    - :fontawesome-solid-closed-captioning: 고품질(수동 큐레이팅된) 자막. :material-subtitles: 아이콘으로 켤 수 있습니다
    - :material-bookmark: 페이지 제목에 해당하는 타임라인의 비디오 챕터

## 교육 환경 시작

GitHub Codespaces에서 제공하는 사전 구축된 환경을 사용하려면 아래 "Open in GitHub Codespaces" 버튼을 클릭하십시오. 다른 옵션은 [환경 옵션](../envsetup/index.md)을 참조하십시오.

환경이 로드되는 동안 계속 읽을 수 있도록 새 브라우저 탭이나 창에서 교육 환경을 여는 것이 좋습니다(장비에 따라 우클릭, ctrl-클릭 또는 cmd-클릭 사용).
과정을 진행하려면 이 지침을 병렬로 열어 두어야 합니다.

[![Open in GitHub Codespaces](https://github.com/codespaces/badge.svg)](https://codespaces.new/nextflow-io/training?quickstart=1&ref=master)

### 환경 기초

이 교육 환경에는 교육 과정을 진행하는 데 필요한 모든 소프트웨어, 코드 및 데이터가 포함되어 있으므로 직접 아무것도 설치할 필요가 없습니다.

codespace는 파일 시스템 탐색기, 코드 편집기 및 터미널 셸을 포함하는 VSCode 인터페이스로 설정되어 있습니다.
과정 중 제공되는 모든 지침(예: '파일 열기', '코드 편집' 또는 '이 명령 실행')은 별도로 지정하지 않는 한 VScode 인터페이스의 이 세 부분을 참조합니다.

직접 이 과정을 진행하는 경우 [환경 기초](../envsetup/01_setup.md)를 참조하여 자세한 내용을 확인하십시오.

### 버전 요구 사항

이 교육은 **v2 구문 분석기가 활성화된** Nextflow 25.10.2 이상용으로 설계되었습니다.
로컬 또는 사용자 정의 환경을 사용하는 경우 [여기](../info/nxf_versions.md)에 문서화된 대로 올바른 설정을 사용하고 있는지 확인하십시오.

## 작업 준비

codespace가 실행되면 교육에 들어가기 전에 두 가지를 해야 합니다: 이 특정 과정에 대한 작업 디렉토리 설정 및 제공된 자료 살펴보기.

### 작업 디렉토리 설정

기본적으로 codespace는 모든 교육 과정의 루트에서 작업 디렉토리가 설정된 상태로 열리지만, 이 과정에서는 `hello-nextflow/` 디렉토리에서 작업합니다.

터미널에서 다음 명령을 실행하여 지금 디렉토리를 변경하십시오:

```bash
cd hello-nextflow/
```

VSCode가 이 디렉토리에 초점을 맞추도록 설정하여 파일 탐색기 사이드바에 관련 파일만 표시되도록 할 수 있습니다:

```bash
code .
```

!!! tip

    어떤 이유로든 이 디렉토리에서 벗어난 경우(예: codespace가 중지 상태로 전환된 경우), Github Codespaces 교육 환경 내에서 실행한다고 가정하고 전체 경로를 사용하여 언제든지 돌아올 수 있습니다:

    ```bash
    cd /workspaces/training/hello-nextflow
    ```

이제 내용을 살펴보겠습니다.

### 제공된 자료 탐색

교육 작업 공간의 왼쪽에 있는 파일 탐색기를 사용하여 이 디렉토리의 내용을 탐색할 수 있습니다.
또는 `tree` 명령을 사용할 수 있습니다.

과정 전체에서 `tree`의 출력을 사용하여 디렉토리 구조와 내용을 읽기 쉬운 형태로 표현하며, 때로는 명확성을 위해 약간의 수정을 가합니다.

여기서는 두 번째 수준까지 목차를 생성합니다:

```bash
tree . -L 2
```

??? abstract "디렉토리 내용"

    ```console
    .
    ├── data
    │   └── greetings.csv
    ├── hello-channels.nf
    ├── hello-config.nf
    ├── hello-containers.nf
    ├── hello-modules.nf
    ├── hello-workflow.nf
    ├── hello-world.nf
    ├── nextflow.config
    ├── solutions
    │   ├── 1-hello-world
    │   ├── 2-hello-channels
    │   ├── 3-hello-workflow
    │   ├── 4-hello-modules
    │   ├── 5-hello-containers
    │   └── 6-hello-config
    ├── test-params.json
    └── test-params.yaml
    ```

색상이 있는 상자를 클릭하여 섹션을 확장하고 내용을 봅니다.
예상되는 명령 출력을 간결하게 포함하기 위해 이와 같은 접을 수 있는 섹션을 사용합니다.

- **`.nf` 파일**은 과정의 어느 부분에서 사용되는지에 따라 이름이 지정된 워크플로우 스크립트입니다.

- **`nextflow.config` 파일**은 최소한의 환경 속성을 설정하는 구성 파일입니다.
  지금은 무시해도 됩니다.

- **`data/` 아래의 `greetings.csv` 파일**에는 과정 대부분에서 사용할 입력 데이터가 포함되어 있습니다. Part 2(Channels)에서 처음 소개할 때 설명합니다.

- **`test-params.*` 파일**은 Part 6(Configuration)에서 사용할 구성 파일입니다. 지금은 무시해도 됩니다.

- **`solutions` 디렉토리**에는 과정의 각 단계에서 나온 완성된 워크플로우 스크립트가 포함되어 있습니다.
  작업을 확인하고 문제를 해결하기 위한 참조용입니다.

## 준비 완료 체크리스트

들어갈 준비가 되셨다고 생각하시나요?

- [ ] 이 과정의 목표와 전제 조건을 이해합니다
- [ ] 환경이 실행 중입니다
- [ ] 작업 디렉토리를 적절하게 설정했습니다

모든 항목을 체크할 수 있다면 준비 완료입니다.

**[Part 1: Hello World](./01_hello_world.md)로 계속하려면 이 페이지 오른쪽 하단의 화살표를 클릭하십시오.**
