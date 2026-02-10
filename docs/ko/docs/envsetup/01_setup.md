# GitHub Codespaces

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } AI 지원 번역 - [자세히 알아보기 및 개선 제안](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

GitHub Codespaces는 클라우드의 가상 머신에 의해 지원되는 사전 구성된 교육 환경을 제공할 수 있는 웹 기반 플랫폼입니다.
이 플랫폼은 Github(Microsoft 소유)에서 운영되며, Github 계정이 있는 누구나 무료(사용량 할당량 내에서)로 이용할 수 있습니다.

!!! warning

    조직에 연결된 계정은 특정 추가 제한이 적용될 수 있습니다.
    이 경우 독립적인 개인 계정을 사용하거나 대신 로컬 설치를 사용해야 할 수 있습니다.

## GitHub 계정 생성

[GitHub 홈페이지](https://github.com/)에서 무료 GitHub 계정을 만들 수 있습니다.

## GitHub Codespace 시작하기

GitHub에 로그인한 후 브라우저에서 이 링크를 열어 Nextflow 교육 환경을 엽니다: <https://codespaces.new/nextflow-io/training?quickstart=1&ref=master>

또는 아래 표시된 버튼을 클릭할 수 있습니다. 이 버튼은 각 교육 과정(일반적으로 Orientation 페이지)에서 반복됩니다.

[![Open in GitHub Codespaces](https://github.com/codespaces/badge.svg)](https://codespaces.new/nextflow-io/training?quickstart=1&ref=master)

새 GitHub Codespace를 생성할 수 있는 페이지가 표시됩니다:

![GitHub Codespace 생성](img/codespaces_create.png)

### 구성

일반적인 사용에서는 별도로 구성할 필요가 없습니다.
시작하려는 과정에서 별도로 지정하지 않는 한 기본 버튼을 클릭하여 계속 진행하면 됩니다.

그러나 "Change options" 버튼을 클릭하여 환경을 사용자 정의할 수 있습니다.

??? info "구성 옵션"

    "Change options" 버튼을 클릭하면 다음 사항을 사용자 정의할 수 있습니다:

    #### 브랜치

    다른 버전의 교육 자료를 선택할 수 있습니다.
    `master` 브랜치에는 일반적으로 버그 수정과 최근에 개발되어 승인되었지만 아직 웹사이트에 릴리스되지 않은 자료가 포함되어 있습니다.
    다른 브랜치에는 완전히 기능하지 않을 수 있는 진행 중인 작업이 포함되어 있습니다.

    #### 머신 유형

    교육을 진행하는 데 사용할 가상 머신을 사용자 정의할 수 있습니다.

    더 많은 코어가 있는 머신을 사용하면 Nextflow의 워크플로우 실행 병렬화 기능을 더 잘 활용할 수 있습니다.
    그러나 무료 할당량을 더 빨리 소비하므로 수강하려는 과정 지침에서 권장하지 않는 한 이 설정을 변경하는 것은 권장하지 않습니다.

    할당량에 대한 자세한 내용은 아래 'GitHub Codespaces 할당량'을 참조하십시오.

### 시작 시간

처음으로 새 GitHub Codespaces 환경을 여는 데 몇 분이 걸릴 수 있습니다. 시스템이 가상 머신을 설정해야 하므로 대기 시간이 있더라도 걱정하지 마십시오.
그러나 5분 이상 걸리지 않아야 합니다.

## 교육 인터페이스 탐색

GitHub Codespaces가 로드되면 다음과 유사한 화면이 표시됩니다(계정 기본 설정에 따라 라이트 모드로 열릴 수 있습니다):

![GitHub Codespaces 시작 화면](img/codespaces_welcome.png)

이것은 Nextflow 개발에 사용하기를 권장하는 인기 있는 코드 개발 애플리케이션인 VSCode IDE의 인터페이스입니다.

- **메인 에디터**는 Nextflow 코드 및 기타 텍스트 파일이 열리는 곳입니다. 여기에서 코드를 편집합니다. codespace를 열면 `README.md` 파일의 미리보기가 표시됩니다.
- 메인 에디터 아래의 **터미널**에서 명령을 실행할 수 있습니다. 과정 지침에 제공된 모든 명령줄은 여기에서 실행합니다.
- **사이드바**를 사용하면 환경을 사용자 정의하고 기본 작업(복사, 붙여넣기, 파일 열기, 검색, git 등)을 수행할 수 있습니다. 기본적으로 저장소의 내용을 탐색할 수 있는 파일 탐색기가 열려 있습니다. 탐색기에서 파일을 클릭하면 메인 에디터 창에서 열립니다.

원하는 대로 창 패널의 상대적 비율을 조정할 수 있습니다.

<!-- TODO (future) Link to development best practices side quest? -->

## GitHub Codespaces 사용에 대한 기타 참고 사항

### 세션 재개

환경을 생성한 후에는 쉽게 재개하거나 다시 시작하여 중단한 곳에서 계속할 수 있습니다.
환경은 30분간 비활동 후 타임아웃되며 최대 2주 동안 변경 사항을 저장합니다.

<https://github.com/codespaces/>에서 환경을 다시 열 수 있습니다.
이전 환경이 나열됩니다.
세션을 클릭하여 재개합니다.

![GitHub Codespace 세션 목록](img/codespaces_list.png)

이전 GitHub Codespaces 환경의 URL을 저장해 두었다면 브라우저에서 바로 열 수 있습니다.
또는 처음에 생성할 때 사용한 것과 동일한 버튼을 클릭하면 됩니다:

[![Open in GitHub Codespaces](https://github.com/codespaces/badge.svg)](https://codespaces.new/nextflow-io/training?quickstart=1&ref=master)

이전 세션이 표시되며 기본 옵션은 재개하는 것입니다:

![GitHub Codespace 재개](img/codespaces_resume.png)

### 로컬 머신에 파일 저장

탐색기 패널에서 파일을 저장하려면 파일을 마우스 오른쪽 버튼으로 클릭하고 `Download`를 선택합니다.

### GitHub Codespaces 할당량 관리

GitHub Codespaces는 월 최대 15GB-월 스토리지와 월 120 코어-시간을 제공합니다.
이는 표준 워크스페이스(2코어, 8GB RAM, 32GB 스토리지)를 사용하는 기본 환경 런타임의 약 60시간에 해당합니다.

더 많은 리소스로 생성할 수 있지만(위의 설명 참조), 무료 사용량을 더 빨리 소비하게 되어 이 공간에 대한 액세스 시간이 줄어듭니다.
예를 들어, 2코어 기본값 대신 4코어 머신을 선택하면 할당량이 절반의 시간 내에 소진됩니다.

선택적으로 더 많은 리소스에 대한 액세스를 구매할 수 있습니다.

자세한 내용은 GitHub 문서를 참조하십시오:
[GitHub Codespaces 청구에 대해](https://docs.github.com/en/billing/managing-billing-for-your-products/managing-billing-for-github-codespaces/about-billing-for-github-codespaces)
