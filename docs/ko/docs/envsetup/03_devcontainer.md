# 로컬 Devcontainers

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } AI 지원 번역 - [자세히 알아보기 및 개선 제안](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

로컬 Docker 설치가 있거나 설치할 의향이 있다면, 이 자료로 로컬에서 작업하는 가장 쉬운 방법은 Visual Studio Code의 devcontainer 기능을 사용하는 것입니다. 이 방식은 수동 설치 없이 필요한 모든 도구와 의존성을 제공합니다.

## 요구 사항

로컬 devcontainer 설정을 사용하려면 다음이 필요합니다:

- [Visual Studio Code](https://code.visualstudio.com/)
- 로컬 Docker 설치, 예:
  - [Docker Desktop](https://docs.docker.com/get-docker/) (Windows/macOS용)
  - [Docker Engine](https://docs.docker.com/engine/install/) (Linux용)
  - [Colima](https://github.com/abiosoft/colima) (macOS 대안)
- [Docker Buildx](https://docs.docker.com/build/concepts/overview/#install-buildx) (Docker Desktop에 포함되어 있지만 다른 Docker 설정에서는 별도 설치가 필요할 수 있음)
- VS Code용 [Dev Containers 확장](https://marketplace.visualstudio.com/items?itemName=ms-vscode-remote.remote-containers)

devcontainer를 열기 전에 Docker 설치가 실행 중이어야 합니다.

Docker buildx가 사용 가능한지 확인하려면 다음을 실행하십시오:

```bash
docker buildx version
```

이 명령이 실패하면 진행하기 전에 buildx 확장을 설치해야 합니다.

## 설정 지침

VS Code devcontainers를 사용하여 로컬 환경을 설정하려면 다음 단계를 따르십시오:

### VS Code에 "Dev Containers" 확장 설치

- VS Code 열기
- Extensions로 이동 (Ctrl+Shift+X 또는 macOS에서 Cmd+Shift+X)
- "Dev Containers" 검색
- "Install" 클릭

![VS Code에서 Dev Containers 확장 프로그램 설치](img/install_extension.png)

### 저장소 복제:

```bash
git clone https://github.com/nextflow-io/training.git
cd training
```

### VS Code에서 저장소 열기:

- VS Code 실행
- 메뉴에서 **File -> Open Folder** 선택
- 방금 복제한 training 저장소 폴더로 이동하여 선택
- **Open** 클릭

### 컨테이너에서 다시 열기

VS Code에서 "Reopen in Container"라는 메시지가 표시되면 클릭합니다. 또는:

- F1 누르기 (또는 macOS에서 Ctrl+Shift+P / Cmd+Shift+P)
- "Dev Containers: Reopen in Container" 입력
- **중요**: 구성을 선택하라는 메시지가 표시되면 **local-dev** devcontainer 구성을 선택합니다

![컨테이너에서 다시 열기 프롬프트](img/reopen_prompt.png)

![로컬 구성 선택](img/select_local_config.png)

컨테이너가 빌드될 때까지 기다립니다. 처음에는 필요한 모든 구성 요소를 다운로드하고 설정해야 하므로 몇 분이 걸릴 수 있습니다.

컨테이너가 빌드되고 실행되면 다음을 포함한 필요한 모든 도구가 설치된 완전히 구성된 환경이 제공됩니다:

- Java
- Nextflow
- Docker
- Git
- 교육에 필요한 기타 모든 의존성

![devcontainer가 실행 중인 VS Code](img/running_container.png)

## Devcontainers 사용의 이점

devcontainer 방식을 사용하면 여러 가지 이점이 있습니다:

- **일관성**: 다른 머신에서도 일관된 개발 환경 보장
- **단순성**: 모든 의존성이 사전 설치되고 구성됨
- **격리**: 개발 환경이 로컬 시스템과 격리됨
- **재현성**: devcontainer를 사용하는 모든 사람이 동일한 설정을 얻음
- **수동 설치 불필요**: Java, Nextflow 및 기타 도구를 수동으로 설치할 필요 없음

## 환경 확인

devcontainer가 실행되면 다음을 실행하여 모든 것이 올바르게 설정되었는지 확인할 수 있습니다:

```bash
nextflow info
```

Nextflow 버전 및 런타임 정보가 표시되어 환경이 올바르게 구성되었음을 확인할 수 있습니다.

## 문제 해결

devcontainer 설정에 문제가 발생하면:

1. devcontainer를 열기 전에 Docker 설치(Docker Desktop, Colima, Docker Engine 등)가 실행 중인지 확인하십시오
2. 메시지가 표시되면 **local-dev** 구성을 선택했는지 확인하십시오
3. `docker buildx version`을 실행하여 Docker buildx가 설치되어 있고 작동하는지 확인하십시오
4. 컨테이너 빌드가 실패하면 "Dev Containers: Rebuild Container" 명령을 실행하여 다시 빌드해 보십시오
5. 지속적인 문제의 경우 [VS Code Dev Containers 문제 해결 가이드](https://code.visualstudio.com/docs/devcontainers/troubleshooting)를 참조하십시오
