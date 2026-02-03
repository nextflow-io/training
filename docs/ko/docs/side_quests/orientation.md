# 오리엔테이션

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } AI 지원 번역 - [자세히 알아보기 및 개선 제안](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

GitHub Codespaces 환경에는 이 교육 과정을 진행하는 데 필요한 모든 소프트웨어, 코드 및 데이터가 포함되어 있으므로, 별도로 설치할 필요가 없습니다.
그러나 로그인하려면 (무료) 계정이 필요하며, 인터페이스에 익숙해지는 데 몇 분 정도 시간을 할애하시기 바랍니다.

아직 수행하지 않으셨다면, 계속 진행하기 전에 [이 링크](../../envsetup/)를 따라가 주시기 바랍니다.

## 제공되는 자료

이 교육 과정을 진행하는 동안 `side-quests/` 디렉토리에서 작업하게 됩니다.
이 디렉토리에는 필요한 모든 코드 파일, 테스트 데이터 및 부속 파일이 포함되어 있습니다.

이 디렉토리의 내용을 자유롭게 탐색하시기 바랍니다. 가장 쉬운 방법은 GitHub Codespaces 작업 공간의 왼쪽에 있는 파일 탐색기를 사용하는 것입니다.
또는 `tree` 명령을 사용할 수도 있습니다.
과정 전반에 걸쳐 `tree`의 출력을 사용하여 디렉토리 구조와 내용을 읽기 쉬운 형태로 표현하며, 때로는 명확성을 위해 약간 수정하기도 합니다.

여기서는 두 번째 수준까지 목차를 생성합니다:

```bash
tree . -L 2
```

`side-quests` 내부에서 이를 실행하면 다음과 같은 출력이 표시됩니다:

```console title="디렉토리 내용"
.
├── metadata
├── nf-core
├── nf-test
├── solutions
├── splitting_and_grouping
└── workflows_of_workflows
```

**시작하기 위해 알아야 할 내용을 요약하면 다음과 같습니다:**

- **각 디렉토리는 개별 사이드 퀘스트에 해당합니다.**
  각 내용은 해당 사이드 퀘스트 페이지에 자세히 설명되어 있습니다.

- **`solutions` 디렉토리**에는 각 사이드 퀘스트의 다양한 단계를 진행하여 생성된 완성된 워크플로우 및/또는 모듈 스크립트가 포함되어 있습니다.
  이는 작업을 확인하고 문제를 해결하는 데 참조로 사용하도록 제공됩니다.

!!!tip

    어떤 이유로든 이 디렉토리를 벗어났다면, 다음 명령을 실행하여 언제든지 돌아올 수 있습니다:

    ```bash
    cd /workspaces/training/side-quests
    ```

이제 과정을 시작하려면 이 페이지의 오른쪽 하단에 있는 화살표를 클릭하시기 바랍니다.
