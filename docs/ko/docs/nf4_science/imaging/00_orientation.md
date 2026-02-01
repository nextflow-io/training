# 오리엔테이션

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } AI 지원 번역 - [자세히 알아보기 및 개선 제안](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

이 오리엔테이션은 "Open in GitHub Codespaces" 버튼을 클릭하여 이미 교육 환경을 열었다고 가정합니다.
아직 열지 않았다면, 지금 여십시오. 가능하면 이 지침을 참조할 수 있도록 두 번째 브라우저 창이나 탭에서 여는 것이 좋습니다.

[![Open in GitHub Codespaces](https://github.com/codespaces/badge.svg)](https://github.com/codespaces/new?hide_repo_select=true&ref=master&repo=290790519&skip_quickstart=true&machine=premiumLinux&devcontainer_path=.devcontainer%2Fdevcontainer.json)

!!!warning "머신 크기 요구사항"

    이 교육 과정을 위한 Codespace를 생성할 때 **8코어 머신**을 선택하십시오. 생물 이미징 워크플로우에는 추가 컴퓨팅 리소스가 필요합니다.

## GitHub Codespaces

GitHub Codespaces 환경에는 이 교육 과정을 진행하는 데 필요한 모든 소프트웨어, 코드 및 데이터가 포함되어 있으므로 별도로 설치할 필요가 없습니다.
다만, 로그인하려면 (무료) GitHub 계정이 필요하며, 인터페이스에 익숙하지 않으시다면 [GitHub Codespaces 오리엔테이션](../../envsetup/index.md) 단기 과정을 완료하여 몇 분간 익숙해지시기 바랍니다.

## Docker 이미지 사전 다운로드

Codespace를 열었으면, 이 교육 과정에 필요한 모든 Docker 이미지를 사전에 다운로드하겠습니다.
이렇게 하면 나중에 시간을 절약하고 워크플로우가 원활하게 실행됩니다.

새 터미널 탭을 열고 다음 명령을 실행하십시오:

```bash
nextflow run nf-core/molkart -profile docker,test -stub -resume --outdir results
```

이 명령은 백그라운드에서 필요한 모든 Docker 이미지를 다운로드합니다.
이 작업이 실행되는 동안 오리엔테이션의 나머지 부분을 계속 진행할 수 있습니다.

!!!tip "팁"

    `-stub` 플래그를 사용하면 실제 데이터를 처리하지 않고 파이프라인이 빠르게 실행되므로 이미지 다운로드에 완벽합니다. 터미널 탭에서 진행 상황을 모니터링할 수 있습니다.

## 작업 디렉토리

이 교육 과정 전체에서 `nf4-science/imaging/` 디렉토리에서 작업합니다.

터미널에서 다음 명령을 실행하여 지금 디렉토리를 변경하십시오:

```bash
cd nf4-science/imaging/
```

!!!tip "팁"

    어떤 이유로든 이 디렉토리에서 벗어났다면, GitHub Codespaces 교육 환경 내에서 실행 중이라고 가정할 때 항상 전체 경로를 사용하여 돌아올 수 있습니다:

    ```bash
    cd /workspaces/training/nf4-science/imaging
    ```

**이제 과정을 시작하려면 이 페이지의 오른쪽 하단 모서리에 있는 화살표를 클릭하십시오.**
