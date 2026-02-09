# 시작하기

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } AI 지원 번역 - [자세히 알아보고 개선 사항 제안하기](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

## 교육 환경 시작하기

GitHub Codespaces에서 제공하는 사전 구축된 환경을 사용하려면 아래의 "Open in GitHub Codespaces" 버튼을 클릭하세요. 다른 옵션은 [환경 옵션](../../envsetup/index.md)을 참조하세요.

환경이 로드되는 동안 계속 읽을 수 있도록 새 브라우저 탭이나 창에서 교육 환경을 여는 것을 권장합니다(장비에 따라 마우스 오른쪽 클릭, ctrl-클릭 또는 cmd-클릭 사용).
과정을 진행하려면 이 지침을 병렬로 열어 두어야 합니다.

[![Open in GitHub Codespaces](https://github.com/codespaces/badge.svg)](https://codespaces.new/nextflow-io/training?quickstart=1&ref=master)

### 환경 기본 사항

이 교육 환경에는 교육 과정을 진행하는 데 필요한 모든 소프트웨어, 코드 및 데이터가 포함되어 있으므로 직접 설치할 필요가 없습니다.

Codespace는 VSCode 인터페이스로 설정되어 있으며, 파일 시스템 탐색기, 코드 편집기 및 터미널 셸이 포함되어 있습니다.
과정 중에 제공되는 모든 지침(예: '파일 열기', '코드 편집' 또는 '이 명령 실행')은 별도로 명시되지 않는 한 VSCode 인터페이스의 이 세 부분을 참조합니다.

이 과정을 혼자 진행하는 경우, 자세한 내용은 [환경 기본 사항](../../envsetup/01_setup.md)을 숙지하시기 바랍니다.

### 버전 요구 사항

이 교육은 Nextflow 25.10.2 이상 **및 v2 구문 분석기가 활성화된 상태**를 위해 설계되었습니다.
로컬 또는 사용자 정의 환경을 사용하는 경우, [여기](../../info/nxf_versions.md)에 문서화된 대로 올바른 설정을 사용하고 있는지 확인하세요.

## 작업 준비하기

Codespace가 실행되면 교육에 들어가기 전에 두 가지 작업을 수행해야 합니다. 이 특정 과정의 작업 디렉토리를 설정하고 제공된 자료를 살펴보는 것입니다.

### 작업 디렉토리 설정

기본적으로 Codespace는 모든 교육 과정의 루트에 작업 디렉토리가 설정된 상태로 열리지만, 이 과정에서는 `nf4-science/genomics/` 디렉토리에서 작업합니다.

터미널에서 다음 명령을 실행하여 지금 디렉토리를 변경하세요:

```bash
cd nf4-science/genomics/
```

파일 탐색기 사이드바에 관련 파일만 표시되도록 VSCode가 이 디렉토리에 집중하도록 설정할 수 있습니다:

```bash
code .
```

!!! tip "팁"

    어떤 이유로든 이 디렉토리에서 벗어난 경우(예: Codespace가 중지(sleep) 상태가 된 경우), Github Codespaces 교육 환경 내에서 실행한다고 가정하면 전체 경로를 사용하여 언제든지 돌아올 수 있습니다:

    ```bash
    cd /workspaces/training/nf4-science/genomics
    ```

이제 내용을 살펴보겠습니다.

### 제공된 자료 탐색

교육 작업 공간의 왼쪽에 있는 파일 탐색기를 사용하여 이 디렉토리의 내용을 탐색할 수 있습니다.
또는 `tree` 명령을 사용할 수도 있습니다.

과정 전반에 걸쳐 `tree`의 출력을 사용하여 디렉토리 구조와 내용을 읽기 쉬운 형태로 표현하며, 때로는 명확성을 위해 약간 수정합니다.

여기서는 두 번째 수준까지 목차를 생성합니다:

```bash
tree . -L 2
```

??? abstract "디렉토리 내용"

    ```console
    .
    ├── data
    │   ├── bam
    │   ├── ref
    │   ├── sample_bams.txt
    │   └── samplesheet.csv
    ├── genomics.nf
    ├── modules
    │   ├── gatk_haplotypecaller.nf
    │   └── samtools_index.nf
    ├── nextflow.config
    └── solutions
        ├── modules
        ├── nf-test.config
        ├── part2
        └── tests

    8 directories, 8 files
    ```

색상이 있는 상자를 클릭하여 섹션을 확장하고 내용을 확인하세요.
이와 같은 접을 수 있는 섹션을 사용하여 예상되는 명령 출력과 디렉토리 및 파일 내용을 간결하게 표시합니다.

- **`genomics.nf` 파일**은 과정을 진행하면서 구축할 워크플로우 스크립트입니다.

- **`modules` 디렉토리**에는 과정 중에 채울 스켈레톤 모듈 파일이 포함되어 있습니다.

- **`nextflow.config` 파일**은 최소한의 환경 속성을 설정하는 구성 파일입니다.
  지금은 무시해도 됩니다.

- **`data` 디렉토리**에는 입력 데이터와 관련 리소스가 포함되어 있으며, 과정 후반부에 설명됩니다.

- **`solutions` 디렉토리**에는 완성된 모듈 파일과 Part 3의 시작점으로 사용할 수 있는 Part 2 해결책이 포함되어 있습니다.
  이들은 작업을 확인하고 문제를 해결하는 데 참조용으로 사용하기 위한 것입니다.

## 준비 상태 체크리스트

준비가 되었다고 생각하시나요?

- [ ] 이 과정의 목표와 전제 조건을 이해했습니다
- [ ] 환경이 실행 중입니다
- [ ] 작업 디렉토리를 적절하게 설정했습니다

모든 항목을 체크할 수 있다면 시작할 준비가 된 것입니다.

**[Part 1: 메서드 개요 및 수동 테스트](./01_method.md)로 계속하려면 이 페이지의 오른쪽 하단 모서리에 있는 화살표를 클릭하세요.**
