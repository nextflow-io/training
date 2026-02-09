# GitHub Codespaces

GitHub Codespaces는 클라우드의 가상 머신을 기반으로 사전 구성된 교육 환경을 제공하는 웹 기반 플랫폼입니다.
이 플랫폼은 GitHub(Microsoft 소유)에서 운영하며, GitHub 계정이 있는 사람이라면 누구나 무료로(사용량 할당량 내에서) 접근할 수 있습니다.

!!! warning "경고"

    조직에 연결된 계정은 특정 추가 제한 사항이 적용될 수 있습니다.
    해당되는 경우, 독립적인 개인 계정을 사용하거나 로컬 설치를 대신 사용해야 할 수 있습니다.

## GitHub 계정 생성

[GitHub 홈페이지](https://github.com/)에서 무료 GitHub 계정을 생성할 수 있습니다.

## GitHub Codespace 실행하기

GitHub에 로그인한 후, 브라우저에서 다음 링크를 열어 Nextflow 교육 환경을 시작하세요: <https://codespaces.new/nextflow-io/training?quickstart=1&ref=master>

또는 아래에 표시된 버튼을 클릭할 수도 있습니다. 이 버튼은 각 교육 과정(일반적으로 오리엔테이션 페이지)에 반복적으로 표시됩니다.

[![Open in GitHub Codespaces](https://github.com/codespaces/badge.svg)](https://codespaces.new/nextflow-io/training?quickstart=1&ref=master)

새 GitHub Codespace를 생성할 수 있는 페이지가 표시됩니다:

![GitHub Codespace 생성](img/codespaces_create.png)

### 구성

일반적인 사용의 경우, 별도로 구성할 필요가 없습니다.
시작하려는 과정에서 별도로 명시하지 않는 한, 메인 버튼을 클릭하여 계속 진행하면 됩니다.

그러나 "Change options" 버튼을 클릭하여 환경을 사용자 정의할 수 있습니다.

??? info "구성 옵션"

    "Change options" 버튼을 클릭하면 다음을 사용자 정의할 수 있는 옵션이 제공됩니다:

    #### Branch

    교육 자료의 다른 버전을 선택할 수 있습니다.
    `master` 브랜치는 일반적으로 버그 수정과 최근 개발되어 승인되었지만 아직 웹사이트에 릴리스되지 않은 자료를 포함합니다.
    다른 브랜치는 완전히 작동하지 않을 수 있는 진행 중인 작업을 포함합니다.

    #### Machine type

    교육을 진행하는 데 사용할 가상 머신을 사용자 정의할 수 있습니다.

    더 많은 코어를 가진 머신을 사용하면 워크플로우 실행을 병렬화하는 Nextflow의 기능을 더 잘 활용할 수 있습니다.
    그러나 무료 할당량을 더 빨리 소비하므로, 수강하려는 과정의 지침에서 권장하지 않는 한 이 설정을 변경하지 않는 것이 좋습니다.

    할당량에 대한 자세한 내용은 아래의 'GitHub Codespaces 할당량'을 참조하세요.

### 시작 시간

새 GitHub Codespaces 환경을 처음 열 때는 시스템이 가상 머신을 설정해야 하므로 몇 분이 걸릴 수 있습니다. 대기 시간이 있어도 걱정하지 마세요.
그러나 5분 이상 걸리지는 않습니다.

## 교육 인터페이스 탐색

GitHub Codespaces가 로드되면 다음과 유사한 화면이 표시됩니다(계정 설정에 따라 라이트 모드로 열릴 수 있습니다):

![GitHub Codespaces 환영 화면](img/codespaces_welcome.png)

이것은 Nextflow 개발에 권장하는 인기 있는 코드 개발 애플리케이션인 VSCode IDE의 인터페이스입니다.

- **메인 편집기**는 Nextflow 코드와 기타 텍스트 파일이 열리는 곳입니다. 여기에서 코드를 편집합니다. Codespace를 열면 `README.md` 파일의 미리보기가 표시됩니다.
- 메인 편집기 아래의 **터미널**에서 명령을 실행할 수 있습니다. 과정 지침에 제공된 모든 명령줄을 여기에서 실행합니다.
- **사이드바**를 사용하면 환경을 사용자 정의하고 기본 작업(복사, 붙여넣기, 파일 열기, 검색, git 등)을 수행할 수 있습니다. 기본적으로 파일 탐색기로 열려 있으며, 저장소의 내용을 탐색할 수 있습니다. 탐색기에서 파일을 클릭하면 메인 편집기 창에서 열립니다.

원하는 대로 창 패널의 상대적 비율을 조정할 수 있습니다.

<!-- TODO (future) Link to development best practices side quest? -->

## GitHub Codespaces 사용에 대한 기타 참고 사항

### 세션 재개

환경을 생성한 후에는 쉽게 재개하거나 재시작하여 중단한 지점부터 계속할 수 있습니다.
환경은 30분 동안 활동이 없으면 중지(sleep) 상태가 되며, 최대 2주 동안 변경 사항을 저장합니다.

<https://github.com/codespaces/>에서 환경을 다시 열 수 있습니다.
이전 환경이 나열됩니다.
세션을 클릭하여 재개하세요.

![GitHub Codespace 세션 리스트](img/codespaces_list.png)

이전 GitHub Codespaces 환경의 URL을 저장한 경우, 브라우저에서 간단히 열 수 있습니다.
또는 처음 생성할 때 사용한 것과 동일한 버튼을 클릭하세요:

[![Open in GitHub Codespaces](https://github.com/codespaces/badge.svg)](https://codespaces.new/nextflow-io/training?quickstart=1&ref=master)

이전 세션이 표시되며, 기본 옵션은 재개하는 것입니다:

![GitHub Codespace 재개](img/codespaces_resume.png)

### 로컬 머신에 파일 저장

탐색기 패널에서 파일을 저장하려면 파일을 마우스 오른쪽 버튼으로 클릭하고 `Download`를 선택하세요.

### GitHub Codespaces 할당량 관리

GitHub Codespaces는 월 최대 15GB-월 스토리지와 월 120 코어-시간을 제공합니다.
이는 표준 작업 공간(2 코어, 8GB RAM, 32GB 스토리지)을 사용하는 기본 환경 런타임의 약 60시간에 해당합니다.

더 많은 리소스로 생성할 수 있지만(위의 설명 참조), 이는 무료 사용량을 더 빨리 소비하며 이 공간에 대한 접근 시간이 줄어듭니다.
예를 들어, 2 코어 기본값 대신 4 코어 머신을 선택하면 할당량이 절반의 시간에 소진됩니다.

선택적으로 더 많은 리소스에 대한 접근 권한을 구매할 수 있습니다.

자세한 내용은 GitHub 문서를 참조하세요:
[About billing for GitHub Codespaces](https://docs.github.com/en/billing/managing-billing-for-your-products/managing-billing-for-github-codespaces/about-billing-for-github-codespaces)
