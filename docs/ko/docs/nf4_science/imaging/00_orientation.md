# 오리엔테이션

이 오리엔테이션은 "Open in GitHub Codespaces" 버튼을 클릭하여 교육 환경을 이미 열었다고 가정합니다.
아직 열지 않았다면, 지금 여세요. 가능하면 두 번째 브라우저 창이나 탭에서 열어서 이 지침을 참고할 수 있도록 하는 것이 좋습니다.

[![Open in GitHub Codespaces](https://github.com/codespaces/badge.svg)](https://github.com/codespaces/new?hide_repo_select=true&ref=master&repo=290790519&skip_quickstart=true&machine=premiumLinux&devcontainer_path=.devcontainer%2Fdevcontainer.json)

!!!warning "머신 크기 요구사항"

    이 교육 과정을 위해 Codespace를 생성할 때 반드시 **8-코어 머신**을 선택하세요. 생물 이미징 워크플로우는 추가 컴퓨팅 리소스가 필요합니다.

## GitHub Codespaces

GitHub Codespaces 환경에는 이 교육 과정을 진행하는 데 필요한 모든 소프트웨어, 코드 및 데이터가 포함되어 있으므로 직접 설치할 필요가 없습니다.
하지만 로그인하려면 (무료) GitHub 계정이 필요하며, 인터페이스에 익숙하지 않다면 [GitHub Codespaces 오리엔테이션](../../envsetup/index.md) 단기 과정을 완료하여 몇 분 동안 인터페이스에 익숙해지는 것이 좋습니다.

## Docker 이미지 사전 다운로드

Codespace를 열었다면, 이제 이 교육 과정에 필요한 모든 Docker 이미지를 사전에 다운로드하겠습니다.
이렇게 하면 나중에 시간을 절약하고 워크플로우가 원활하게 실행되도록 할 수 있습니다.

새 터미널 탭을 열고 다음 명령을 실행하세요:

```bash
nextflow run nf-core/molkart -profile docker,test -stub -resume --outdir results
```

이 명령은 필요한 모든 Docker 이미지를 백그라운드에서 다운로드합니다.
이 작업이 실행되는 동안 나머지 오리엔테이션을 계속 진행할 수 있습니다.

!!!tip

    `-stub` 플래그를 사용하면 실제 데이터를 처리하지 않고 파이프라인이 빠르게 실행되므로 이미지 다운로드에 적합합니다. 터미널 탭에서 진행 상황을 모니터링할 수 있습니다.

## 작업 디렉토리

이 교육 과정 전체에서 `nf4-science/imaging/` 디렉토리에서 작업할 것입니다.

터미널에서 다음 명령을 실행하여 지금 디렉토리를 변경하세요:

```bash
cd nf4-science/imaging/
```

!!!tip

    어떤 이유로든 이 디렉토리에서 벗어났다면, GitHub Codespaces 교육 환경에서 실행하는 경우 전체 경로를 사용하여 언제든지 돌아올 수 있습니다:

    ```bash
    cd /workspaces/training/nf4-science/imaging
    ```

**이제 과정을 시작하려면 이 페이지의 오른쪽 하단 모서리에 있는 화살표를 클릭하세요.**
