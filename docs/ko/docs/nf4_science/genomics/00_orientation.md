# 시작하기

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } AI 지원 번역 - [자세히 알아보기 및 개선 제안](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

## 교육 환경 시작하기

GitHub Codespaces에서 제공하는 사전 구축된 환경을 사용하려면 아래의 "Open in GitHub Codespaces" 버튼을 클릭하십시오. 다른 옵션은 [환경 옵션](../../envsetup/index.md)을 참조하십시오.

환경이 로드되는 동안 계속 읽을 수 있도록 새 브라우저 탭이나 창에서 교육 환경을 여는 것을 권장합니다(사용하는 장비에 따라 마우스 오른쪽 버튼 클릭, ctrl-클릭 또는 cmd-클릭을 사용하십시오).
과정을 진행하려면 이 지침을 병렬로 열어 두어야 합니다.

[![Open in GitHub Codespaces](https://github.com/codespaces/badge.svg)](https://codespaces.new/nextflow-io/training?quickstart=1&ref=master)

### 환경 기본 사항

이 교육 환경에는 교육 과정을 진행하는 데 필요한 모든 소프트웨어, 코드 및 데이터가 포함되어 있으므로 직접 설치할 필요가 없습니다.

코드스페이스는 파일 시스템 탐색기, 코드 편집기, 터미널 셸이 포함된 VSCode 인터페이스로 설정되어 있습니다.
과정 중에 제공되는 모든 지침(예: '파일 열기', '코드 편집' 또는 '이 명령 실행')은 별도로 명시되지 않는 한 VScode 인터페이스의 이 세 부분을 참조합니다.

혼자서 이 과정을 진행하는 경우, 자세한 내용은 [환경 기본 사항](../../envsetup/01_setup.md)을 숙지하십시오.

### 버전 요구 사항

이 교육은 Nextflow 25.10.2 이상 **및 v2 구문 분석기가 활성화된** 버전을 위해 설계되었습니다.
로컬 또는 사용자 지정 환경을 사용하는 경우, [여기](../../info/nxf_versions.md)에 문서화된 대로 올바른 설정을 사용하고 있는지 확인하십시오.

## 작업 준비하기

코드스페이스가 실행되면, 교육에 들어가기 전에 두 가지를 해야 합니다: 이 특정 과정의 작업 디렉토리를 설정하고 제공된 자료를 살펴보는 것입니다.

### 작업 디렉토리 설정하기

기본적으로 코드스페이스는 모든 교육 과정의 루트에 작업 디렉토리가 설정된 상태로 열리지만, 이 과정에서는 `nf4-science/genomics/` 디렉토리에서 작업할 것입니다.

터미널에서 다음 명령을 실행하여 지금 디렉토리를 변경하십시오:

```bash
cd nf4-science/genomics/
```

VSCode가 이 디렉토리에 집중하도록 설정하여 파일 탐색기 사이드바에 관련 파일만 표시되도록 할 수 있습니다:

```bash
code .
```

!!! tip "팁"

    어떤 이유로든 이 디렉토리에서 이동한 경우(예: 코드스페이스가 중지(sleep) 상태가 된 경우), GitHub Codespaces 교육 환경 내에서 실행한다고 가정할 때 항상 전체 경로를 사용하여 돌아올 수 있습니다:

    ```bash
    cd /workspaces/training/nf4-science/genomics
    ```

이제 내용을 살펴보겠습니다.

### 제공된 자료 탐색하기

교육 작업 공간 왼쪽의 파일 탐색기를 사용하여 이 디렉토리의 내용을 탐색할 수 있습니다.
또는 `tree` 명령을 사용할 수도 있습니다.

과정 전반에 걸쳐 `tree`의 출력을 사용하여 디렉토리 구조와 내용을 읽기 쉬운 형태로 표현하며, 때로는 명확성을 위해 약간 수정합니다.

여기서는 두 번째 레벨까지 목차를 생성합니다:

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

색상이 있는 상자를 클릭하여 섹션을 확장하고 내용을 확인하십시오.
이와 같은 접을 수 있는 섹션을 사용하여 예상되는 명령 출력과 디렉토리 및 파일 내용을 간결한 방식으로 표시합니다.

- **`genomics.nf` 파일**은 과정을 진행하면서 구축할 워크플로우 스크립트입니다.

- **`modules` 디렉토리**에는 과정 중에 채워 넣을 스켈레톤 모듈 파일이 포함되어 있습니다.

- **`nextflow.config` 파일**은 최소한의 환경 속성을 설정하는 구성 파일입니다.
  지금은 무시하셔도 됩니다.

- **`data` 디렉토리**에는 과정 후반부에 설명될 입력 데이터와 관련 리소스가 포함되어 있습니다.

- **`solutions` 디렉토리**에는 완성된 모듈 파일과 Part 3의 시작점으로 사용할 수 있는 Part 2 해결책이 포함되어 있습니다.
  이는 작업을 확인하고 문제를 해결하기 위한 참조 자료로 사용됩니다.

## 준비 확인 목록

준비가 되셨다고 생각하십니까?

- [ ] 이 과정의 목표와 전제 조건을 이해했습니다
- [ ] 환경이 실행되고 있습니다
- [ ] 작업 디렉토리를 적절하게 설정했습니다

모든 항목을 확인할 수 있다면 시작할 준비가 된 것입니다.

**[Part 1: 방법 개요 및 수동 테스트](./01_method.md)로 계속 진행하려면 이 페이지 오른쪽 하단의 화살표를 클릭하십시오.**
