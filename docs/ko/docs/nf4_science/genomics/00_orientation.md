# 오리엔테이션

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } AI 지원 번역 - [자세히 알아보기 및 개선 제안](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

교육 환경에는 이 교육 과정을 진행하는 데 필요한 모든 소프트웨어, 코드 및 데이터가 포함되어 있으므로 직접 설치할 필요가 없습니다.
그러나 로그인하려면 (무료) 계정이 필요하며, 인터페이스에 익숙해지는 데 몇 분 정도 시간을 할애해야 합니다.

아직 완료하지 않으셨다면, 더 진행하기 전에 [이 링크](../../../envsetup/)를 따라가 주십시오.

## 제공되는 자료

이 교육 과정 전반에 걸쳐 `nf4-science/genomics/` 디렉토리에서 작업할 것이며, 교육 작업 공간을 열 때 이 디렉토리로 이동해야 합니다.
이 디렉토리에는 필요한 모든 코드 파일, 테스트 데이터 및 보조 파일이 포함되어 있습니다.

이 디렉토리의 내용을 자유롭게 탐색해 보십시오. 가장 쉬운 방법은 VSCode 인터페이스의 교육 작업 공간 왼쪽에 있는 파일 탐색기를 사용하는 것입니다.
또는 `tree` 명령을 사용할 수 있습니다.
과정 전반에 걸쳐 `tree`의 출력을 사용하여 디렉토리 구조와 내용을 읽기 쉬운 형태로 표현하며, 때로는 명확성을 위해 약간 수정합니다.

여기서는 두 번째 레벨까지 목차를 생성합니다:

```bash
tree . -L 2
```

`nf4-science/genomics` 내부에서 이 명령을 실행하면 다음과 같은 출력이 표시됩니다:

```console title="디렉토리 내용"

.
├── data
│   ├── bam
│   ├── ref
│   ├── sample_bams.txt
│   └── samplesheet.csv
├── genomics-1.nf
├── genomics-2.nf
├── genomics-3.nf
├── genomics-4.nf
├── nextflow.config
└── solutions
    ├── modules
    ├── nf-test.config
    └── tests

6 directories, 8 files

```

!!!note "참고"

    이것이 많아 보여도 걱정하지 마십시오. 과정의 각 단계에서 관련 부분을 다룰 것입니다.
    이것은 단지 개요를 제공하기 위한 것입니다.

**시작하기 위해 알아야 할 사항에 대한 요약입니다:**

- **`.nf` 파일**은 과정의 어느 부분에서 사용되는지에 따라 이름이 지정된 워크플로우 스크립트입니다.

- **`nextflow.config` 파일**은 최소한의 환경 속성을 설정하는 구성 파일입니다.
  지금은 무시하셔도 됩니다.

- **`data` 디렉토리**에는 입력 데이터와 관련 리소스가 포함되어 있으며, 과정 후반부에 설명됩니다.

- **`solutions` 디렉토리**에는 과정의 Part 3와 4에서 생성되는 모듈 파일과 테스트 구성이 포함되어 있습니다.
  이는 작업을 확인하고 문제를 해결하기 위한 참조 자료로 사용됩니다.

!!!tip "팁"

    어떤 이유로든 이 디렉토리에서 이동한 경우, 항상 다음 명령을 실행하여 돌아올 수 있습니다:

    ```bash
    cd /workspaces/training/nf4-science/genomics
    ```

이제 과정을 시작하려면 이 페이지 오른쪽 하단의 화살표를 클릭하십시오.
