# 수동 설치

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } AI 지원 번역 - [자세히 알아보기 및 개선 제안](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

자체 로컬 환경에서 교육을 실행하는 데 필요한 모든 것을 수동으로 설치할 수 있습니다.

여기서는 표준 POSIX 호환 시스템(노트북과 같은 개인용 머신 가정)에서 이를 수행하는 방법을 문서화했습니다.
특정 시스템에 따라 일부 세부 사항이 다를 수 있다는 점에 유의하십시오.

!!! tip

    진행하기 전에 [Devcontainers 방식](03_devcontainer.md)을 고려해 보셨나요?
    수동 설치 없이 필요한 모든 도구와 의존성을 제공합니다.

## 일반 소프트웨어 요구 사항

Nextflow는 Java가 설치된 모든 POSIX 호환 시스템(Linux, macOS, Windows Subsystem for Linux 등)에서 사용할 수 있습니다.
저희 교육 과정에는 몇 가지 추가 요구 사항이 있습니다.

총합적으로 다음 소프트웨어가 설치되어 있어야 합니다:

- Bash 또는 동등한 셸
- [Java 11 (또는 최대 21까지의 이후 버전)](https://www.oracle.com/technetwork/java/javase/downloads/index.html)
- [Git](https://git-scm.com/)
- [Docker](https://docs.docker.com/get-docker/)
- [Conda](https://conda.io/) 4.5 (또는 이후 버전)
- [VSCode](https://code.visualstudio.com)와 [Nextflow 확장](https://www.nextflow.io/docs/latest/developer-env.html#devenv-nextflow)

VSCode 애플리케이션은 기술적으로 선택 사항이지만 과정을 진행하는 데뿐만 아니라 일반적인 Nextflow 개발 작업에도 사용하는 것을 강력히 권장합니다.

Nextflow 문서 매뉴얼에서는 [환경 설정](https://www.nextflow.io/docs/latest/developer-env.html)에서 이러한 의존성 설치에 대한 지침을 제공합니다.

## Nextflow 및 nf-core 도구

Nextflow 자체와 nf-core 도구를 설치해야 합니다. 아래 링크된 문서에 자세히 설명되어 있습니다:

- [Nextflow 설치](https://www.nextflow.io/docs/latest/install.html)
- [nf-core 도구](https://nf-co.re/docs/nf-core-tools/installation)

Nextflow의 경우 self-install 옵션을, nf-core 도구의 경우 PyPI 옵션을 사용하는 것을 권장합니다.

!!! warning "버전 호환성"

    <!-- Any update to this content needs to be copied to the home page -->
    **2026년 1월 기준, 별도의 언급이 없는 한 모든 Nextflow 교육 과정은 strict v2 구문이 활성화된 Nextflow 버전 25.10.2 이상을 필요로 합니다.**

    버전 요구 사항 및 strict v2 구문에 대한 자세한 내용은 [Nextflow 버전](../info/nxf_versions.md) 가이드를 참조하십시오.

    이전 구문에 해당하는 이전 버전의 교육 자료는 이 웹페이지 메뉴 바의 버전 선택기를 통해 사용할 수 있습니다.

## 교육 자료

교육 자료를 다운로드하는 가장 쉬운 방법은 다음 명령을 사용하여 전체 저장소를 복제하는 것입니다:

```bash
git clone https://github.com/nextflow-io/training.git
```

각 과정에는 자체 디렉토리가 있습니다.
과정을 진행하려면 터미널 창(이상적으로는 VSCode 애플리케이션 내부)을 열고 해당 디렉토리로 `cd`합니다.

그런 다음 웹사이트에 제공된 과정 지침을 따를 수 있습니다.
