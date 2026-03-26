---
title: Nextflow 버전
description: Nextflow 문법 버전의 변화를 이해하고 관리하기
hide:
  - toc
  - footer
---

## 현재 지원되는 Nextflow 문법 버전 및 요구 사항

교육 포털 버전 3.0부터, 과정 인덱스 페이지에 별도로 명시되지 않는 한 모든 교육 과정은 Nextflow 25.10.2 버전을 기준으로 합니다(더 이상 사용되지 않거나 아카이브된 자료는 버전 안내가 포함되지 않을 수 있습니다).

현재 과정에서는 워크플로우 수준의 타입이 지정된 입력과 워크플로우 수준의 출력 지시문을 사용하므로, V2 문법 분석기가 필요합니다.
[Github Codespaces](../envsetup/01_setup.md) 또는 [로컬 devcontainer](../envsetup/03_devcontainer.md)를 통해 제공되는 환경을 사용할 계획이라면, 과정 안내에 별도로 명시된 경우가 아닌 한 추가 작업은 필요하지 않습니다.
단, 직접 환경을 구성하여 교육을 진행할 계획이라면([수동 설치](../envsetup/02_local.md)), Nextflow 버전 25.10.2 이상을 사용하고 v2 문법 분석기를 활성화해야 합니다.

## 이전 버전의 교육 자료

교육 자료는 2025년 2월부터 버전 관리가 시작되었습니다.

**25.10.2 이전** Nextflow 버전과 호환되는 이전 교육 자료는, 모든 페이지 상단의 드롭다운 메뉴에서 교육 자료의 버전 번호를 선택하여 접근할 수 있습니다.
이전 버전의 교육 자료를 선택하면, 교육 환경 링크가 자동으로 해당 버전의 환경을 가리키도록 변경됩니다.

## Nextflow 문법 버전에 관한 추가 정보

Nextflow에는 혼동되기 쉬운 두 가지 버전 개념이 있습니다: **DSL 버전**과 **문법 분석기 버전**입니다.

**DSL1 vs DSL2**는 Nextflow 파이프라인을 작성하는 근본적으로 다른 방식을 의미합니다.
DSL1은 프로세스가 채널을 통해 암묵적으로 연결되던 기존 문법입니다.
Nextflow 20.07에서 도입된 DSL2는 모듈화 기능을 추가했습니다. 다른 파일에서 프로세스와 워크플로우를 가져오는 기능, 명시적인 `workflow` 블록, 그리고 이름이 지정된 프로세스 출력이 포함됩니다.
DSL1은 Nextflow 22.03에서 더 이상 사용되지 않게 되었고, 22.12에서 완전히 제거되었습니다.
현재 모든 Nextflow 코드는 DSL2를 사용합니다.

**문법 분석기 v1 vs v2**는 DSL2 코드를 처리하는 서로 다른 분석기를 의미합니다.
v1 분석기는 기존의 더 유연한 분석기입니다.
v2 분석기는 더 엄격하며, 정적 타이핑(타입이 지정된 입력 및 출력)과 워크플로우 수준의 출력 지시문 같은 새로운 언어 기능을 지원합니다.
또한 v2 분석기는 더 나은 오류 메시지를 제공하고, 런타임이 아닌 파싱 단계에서 더 많은 오류를 감지합니다.
v2 분석기는 Nextflow 26.04에서 기본값이 될 예정입니다.

요약하면, DSL2는 작성하는 언어이고, 문법 분석기 버전은 해당 언어를 얼마나 엄격하게 해석하고 어떤 고급 기능을 사용할 수 있는지를 결정합니다.

### Nextflow 버전 확인 및 설정

`nextflow --version` 명령어를 사용하여 시스템에 설치된 Nextflow 버전을 확인할 수 있습니다.

Nextflow 버전 업데이트 방법에 대한 자세한 내용은 [Nextflow 업데이트](https://www.nextflow.io/docs/latest/updating-nextflow.html) 참조 문서를 확인하세요.

### v2 문법 분석기 활성화

현재 세션에서 v2 문법 분석기를 **활성화**하려면, 터미널에서 다음 명령어를 실행하세요:

```bash
export NXF_SYNTAX_PARSER=v2
```

이 설정을 영구적으로 적용하려면(Nextflow 26.04에서 v2가 기본값이 될 때까지), 셸 프로파일(`~/.bashrc`, `~/.zshrc` 등)에 export 명령어를 추가하세요:

```bash
echo 'export NXF_SYNTAX_PARSER=v2' >> ~/.bashrc
source ~/.bashrc
```

`NXF_SYNTAX_PARSER=v2` 환경 변수는 임시 요구 사항입니다.
Nextflow 26.04부터는 v2 분석기가 기본값이 되므로 이 설정이 더 이상 필요하지 않습니다.

### v2 문법 분석기 비활성화

현재 세션에서 v2 문법 분석기를 **비활성화**하려면, 터미널에서 다음 명령어를 실행하세요:

```bash
export NXF_SYNTAX_PARSER=v1
```

<!-- 26.04 이후 버전에서도 비활성화가 가능한가요? -->

### 기존 코드 마이그레이션

기존 코드를 최신 Nextflow 버전에 맞게 마이그레이션하는 방법은 참조 문서의 [마이그레이션 노트](https://www.nextflow.io/docs/latest/migrations/index.html)를 참고하세요.

최신 릴리스로 마이그레이션할 때 특히 유용한 두 가지 문서는 다음과 같습니다:

- [워크플로우 출력으로 마이그레이션](https://www.nextflow.io/docs/latest/tutorials/workflow-outputs.html)
- [정적 타입으로 마이그레이션](https://www.nextflow.io/docs/latest/tutorials/static-types.html)

이 두 기능은 교육 자료 버전 3.0부터 초급 교육 과정의 일부로 다루어집니다.

마이그레이션하려는 Nextflow 코드의 세대에 따라, `nextflow lint -format` 명령어를 사용하는 Nextflow 린터로 대부분의 작업을 처리할 수 있습니다.
자세한 내용은 [`lint`](https://www.nextflow.io/docs/latest/reference/cli.html#lint) CLI 참조 문서를 확인하세요.

도움이 되셨으면 합니다.
도움이 필요하시면 Slack 또는 포럼에서 문의해 주세요.

---

<div markdown class="homepage_logos">

![Seqera](../assets/img/seqera_logo.png#only-light)

![Seqera](../assets/img/seqera_logo_dark.png#only-dark)

</div>
