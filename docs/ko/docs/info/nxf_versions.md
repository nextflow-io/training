---
title: Nextflow 버전
description: Nextflow 문법 버전의 발전을 이해하고 관리하기
hide:
  - toc
  - footer
---

## 현재 지원되는 Nextflow 문법 버전 및 요구 사항

교육 포털 버전 3.0부터 모든 교육 과정은 과정 인덱스 페이지에 별도로 명시되지 않는 한 Nextflow 25.10.2 버전 릴리스를 기반으로 합니다(버전 공지가 포함되지 않을 수 있는 더 이상 사용되지 않거나 아카이브된 자료 제외).

과정에서 이제 workflow 수준의 typed inputs과 workflow 수준 output directives를 사용하기 때문에 V2 syntax parser를 사용해야 합니다.
[Github Codespaces](../envsetup/01_setup.md) 또는 [로컬 devcontainer](../envsetup/03_devcontainer.md)를 통해 제공하는 환경을 사용할 계획이라면 과정 지침에 특별히 명시되지 않는 한 별도의 작업이 필요하지 않습니다.
그러나 자체 환경에서 교육을 진행할 계획이라면([수동 설치](../envsetup/02_local.md)), v2 syntax parser가 활성화된 Nextflow 버전 25.10.2 이상을 사용해야 합니다.

## 교육 자료의 이전 버전

교육 자료는 2025년 2월부터 버전 관리되고 있습니다.

모든 페이지 상단의 드롭다운 메뉴 항목을 통해 교육 자료의 번호가 매겨진 버전을 확인하고 **25.10.2 이전** Nextflow 버전에서 작동하는 이전 버전의 교육 자료에 접근할 수 있습니다.
이전 버전의 교육 자료를 선택하면 교육 환경 링크가 자동으로 해당 환경 버전을 지정합니다.

## Nextflow 문법 버전에 대한 기타 정보

Nextflow에는 때때로 혼동되는 두 가지 별개의 버전 관리 개념이 있습니다: **DSL 버전**과 **syntax parser 버전**.

**DSL1 대 DSL2**는 Nextflow pipeline을 작성하는 근본적으로 다른 방식을 나타냅니다.
DSL1은 process가 channel을 통해 암묵적으로 연결되는 원래 문법이었습니다.
Nextflow 20.07에 도입된 DSL2는 모듈성 기능을 추가했습니다: 다른 파일에서 process와 workflow를 가져오는 기능, 명시적인 `workflow` 블록, 명명된 process outputs.
DSL1은 Nextflow 22.03에서 더 이상 사용되지 않게 되었고 22.12에서 제거되었습니다.
모든 현대적인 Nextflow 코드는 DSL2를 사용합니다.

**Syntax parser v1 대 v2**는 둘 다 DSL2 코드와 함께 작동하는 서로 다른 parser를 나타냅니다.
v1 parser는 원래의 더 관대한 parser입니다.
v2 parser는 더 엄격하며 static typing(typed inputs 및 outputs)과 workflow 수준 output directives 같은 새로운 언어 기능을 활성화합니다.
v2 parser는 또한 더 나은 오류 메시지를 제공하고 런타임이 아닌 파싱 시점에 더 많은 오류를 잡아냅니다.
v2 parser는 Nextflow 26.04에서 기본값이 됩니다.

요약하면: DSL2는 작성하는 언어이고, syntax parser 버전은 해당 언어가 얼마나 엄격하게 해석되는지와 어떤 고급 기능을 사용할 수 있는지를 결정합니다.

### Nextflow 버전 확인 및 설정

`nextflow --version` 명령을 사용하여 시스템에 설치된 Nextflow 버전을 확인할 수 있습니다.

Nextflow 버전 업데이트 방법에 대한 자세한 내용은 [Nextflow 업데이트](https://www.nextflow.io/docs/latest/updating-nextflow.html)에 대한 참조 문서를 확인하십시오.

### v2 syntax parser 활성화

현재 세션에서 v2 syntax parser를 **활성화**하려면 터미널에서 다음 명령을 실행하십시오:

```bash
export NXF_SYNTAX_PARSER=v2
```

이를 영구적으로 만들려면(Nextflow 26.04에서 v2가 기본값이 될 때까지), export 명령을 shell 프로필(`~/.bashrc`, `~/.zshrc` 등)에 추가하십시오:

```bash
echo 'export NXF_SYNTAX_PARSER=v2' >> ~/.bashrc
source ~/.bashrc
```

`NXF_SYNTAX_PARSER=v2` 환경 변수는 임시 요구 사항입니다.
Nextflow 26.04부터 v2 parser가 기본값이 되며 이 설정은 더 이상 필요하지 않습니다.

### v2 syntax parser 비활성화

현재 세션에서 v2 syntax parser를 **비활성화**하려면 터미널에서 다음 명령을 실행하십시오:

```bash
export NXF_SYNTAX_PARSER=v1
```

<!-- Will it be possible to disable it in versions after 26.04? -->

### 기존 코드 마이그레이션

최신 Nextflow 버전에 맞게 기존 코드를 마이그레이션하는 방법에 대한 안내는 참조 문서의 [마이그레이션 노트](https://www.nextflow.io/docs/latest/migrations/index.html)를 참조하십시오.

다음 두 문서는 최신 릴리스로 마이그레이션하는 데 특히 유용합니다:

- [workflow outputs로 마이그레이션](https://www.nextflow.io/docs/latest/tutorials/workflow-outputs.html)
- [static types로 마이그레이션](https://www.nextflow.io/docs/latest/tutorials/static-types.html)

이 두 기능은 교육 자료 버전 3.0부터 시작하는 초급 교육의 일부로 다룹니다.

마이그레이션하려는 Nextflow 코드의 세대에 따라 `nextflow lint -format` 명령을 사용하여 Nextflow linter로 대부분의 작업을 완료할 수 있습니다.
자세한 내용은 [`lint`](https://www.nextflow.io/docs/latest/reference/cli.html#lint)에 대한 CLI 참조를 확인하십시오.

이것이 도움이 되기를 바랍니다.
도움이 필요하시면 Slack이나 포럼에서 연락해 주십시오.

---

<div markdown class="homepage_logos">

![Seqera](../assets/img/seqera_logo.png#only-light)

![Seqera](../assets/img/seqera_logo_dark.png#only-dark)

</div>
