# 수동 설치

자신의 로컬 환경에서 교육을 실행하는 데 필요한 모든 것을 수동으로 설치할 수 있습니다.

여기서는 표준 POSIX 호환 시스템(노트북과 같은 개인 컴퓨터 가정)에서 이를 수행하는 방법을 문서화했습니다.
특정 시스템에 따라 일부 세부 사항이 다를 수 있다는 점을 유념하세요.

!!! tip "팁"

    진행하기 전에 [Devcontainers 방식](03_devcontainer.md)을 고려해 보셨나요?
    수동 설치 없이 필요한 모든 도구와 의존성을 제공합니다.

## 일반 소프트웨어 요구 사항

Nextflow는 Java가 설치된 모든 POSIX 호환 시스템(Linux, macOS, Windows Subsystem for Linux 등)에서 사용할 수 있습니다.
우리의 교육 과정에는 몇 가지 추가 요구 사항이 있습니다.

총체적으로 다음 소프트웨어를 설치해야 합니다:

- Bash 또는 동등한 셸
- [Java 11 (또는 그 이상, 최대 21)](https://www.oracle.com/technetwork/java/javase/downloads/index.html)
- [Git](https://git-scm.com/)
- [Docker](https://docs.docker.com/get-docker/)
- [Conda](https://conda.io/) 4.5 (또는 그 이상)
- [VSCode](https://code.visualstudio.com)와 [Nextflow 확장](https://www.nextflow.io/docs/latest/developer-env.html#devenv-nextflow)

VSCode 애플리케이션은 기술적으로 선택 사항이지만, 과정을 진행하고 일반적인 Nextflow 개발 작업을 위해 사용하는 것을 강력히 권장합니다.

Nextflow 문서 매뉴얼은 [환경 설정](https://www.nextflow.io/docs/latest/developer-env.html)에서 이러한 의존성을 설치하는 방법에 대한 지침을 제공합니다.

## Nextflow 및 nf-core 도구

아래 링크된 문서에 자세히 설명된 대로 Nextflow 자체와 nf-core 도구를 설치해야 합니다:

- [Nextflow 설치](https://www.nextflow.io/docs/latest/install.html)
- [nf-core 도구](https://nf-co.re/docs/nf-core-tools/installation)

Nextflow는 자체 설치 옵션을, nf-core 도구는 PyPI 옵션을 사용하는 것을 권장합니다.

!!! warning "경고: 버전 호환성"

    <!-- 이 내용을 업데이트하면 홈 페이지에도 복사해야 합니다 -->
    **2026년 1월부터 모든 Nextflow 교육 과정은 별도로 명시되지 않는 한 Nextflow 버전 25.10.2 이상을 요구하며, strict v2 구문이 활성화되어야 합니다.**

    버전 요구 사항 및 strict v2 구문에 대한 자세한 내용은 [Nextflow 버전](../info/nxf_versions.md) 가이드를 참조하세요.

    이전 구문에 해당하는 교육 자료의 이전 버전은 이 웹페이지의 메뉴 바에 있는 버전 선택기를 통해 사용할 수 있습니다.

## 교육 자료

교육 자료를 다운로드하는 가장 쉬운 방법은 다음 명령을 사용하여 전체 저장소를 복제하는 것입니다:

```bash
git clone https://github.com/nextflow-io/training.git
```

각 과정에는 자체 디렉토리가 있습니다.
과정을 진행하려면 터미널 창을 열고(이상적으로는 VSCode 애플리케이션 내부에서) 관련 디렉토리로 `cd`하세요.

그런 다음 웹사이트에 제공된 과정 지침을 따를 수 있습니다.
