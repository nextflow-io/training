---
title: Nextflow Training
description: Welcome to the Nextflow community training portal!
hide:
  - toc
  - footer
---

# Nextflow Training

Nextflow 커뮤니티 교육 포털에 오신 것을 환영합니다!

이 웹사이트에는 여러 가지 다양한 교육 과정이 준비되어 있습니다. 아래로 스크롤하여 여러분에게 맞는 과정을 찾아보세요!

아래 나열된 교육 과정들은 셀프 서비스 리소스로 활용할 수 있도록 설계되었으며, 언제든지 혼자서 학습을 진행할 수 있습니다 (자세한 내용은 Environment Setup을 참고하세요). 하지만 그룹 교육 행사에 참여하면 더욱 효과적으로 학습할 수 있습니다.

- nf-core 커뮤니티에서 정기적으로 무료 온라인 행사를 진행합니다. 자세한 내용은 [nf-core 행사 페이지](https://nf-co.re/events)를 확인하세요.
- Nextflow를 개발한 회사인 Seqera에서는 다양한 교육 행사를 운영하고 있습니다. [Seqera 행사 페이지](https://seqera.io/events/)를 방문하여 'Seqera Sessions' 및 'Nextflow Summit'을 확인해보세요.
- 저희 커뮤니티 팀은 제3자 기관이 주최하는 교육도 정기적으로 진행하고 있습니다. 해당 교육은 일반적으로 주최 측에서 공지 및 신청을 관리합니다.

준비가 되셨다면, 이 페이지나 선택한 과정의 인덱스 페이지에서 'Open in GitHub Codespaces' 버튼을 클릭하세요. 웹 기반의 교육 환경이 열립니다 (GitHub 무료 계정 필요).

[![Open in GitHub Codespaces](https://github.com/codespaces/badge.svg)](https://codespaces.new/nextflow-io/training?quickstart=1&ref=master)

## 교육 환경 설정

!!! exercise "환경 설정"

    !!! tip inline end ""

        :material-laptop:{.nextflow-primary} 처음 환경을 설정해보세요.

    모든 교육 과정에서 사용할 수 있도록 학습 환경을 설정하는 방법에 대한 안내입니다. GitHub Codespaces에 대한 기본 소개와 로컬 컴퓨터에서 직접 작업할 수 있는 대체 설치 방법도 포함되어 있습니다.

    [환경 설정 교육 시작하기 :material-arrow-right:](envsetup/index.md){ .md-button .md-button--primary }

## Nextflow 입문자를 위한 과정

이 과정들은 Nextflow를 완전히 처음 접하는 분들을 위한 기초 수업으로, 특정 분야의 전문 지식 없이도 학습할 수 있도록 구성되어 있습니다. 각 과정은 단계별로 실력을 쌓을 수 있도록 여러 개의 교육 모듈로 이루어져 있습니다.

!!! exercise "Hello Nextflow"

    !!! tip inline end ""

        :material-run-fast:{.nextflow-primary} Nextflow로 파이프라인 개발 배우기

        :fontawesome-brands-youtube:{.youtube} 동영상 자료 제공

    이 과정은 자신만의 파이프라인을 개발하고자 하는 입문자를 위한 교육입니다. Nextflow 언어의 핵심 구성 요소들을 충분히 상세하게 다루며, 단순하지만 완전한 기능을 갖춘 파이프라인을 개발할 수 있도록 도와줍니다. 또한 파이프라인 설계, 개발, 설정에 필요한 주요 개념들도 함께 다룹니다.

    [Hello Nextflow 교육 시작하기 :material-arrow-right:](hello_nextflow/index.md){ .md-button .md-button--primary }

!!! info ""

    **곧 출시 예정:** "Nextflow Run" — Nextflow 파이프라인 실행 배우기 (코드 개발 없이 실행만 다룸)

<!-- COMMENTED OUT UNTIL THIS IS READY
!!! exercise "Nextflow Run"

    !!! tip inline end ""

        :material-run-fast:{.nextflow-primary} Nextflow 파이프라인 실행 방법 배우기

    이 과정은 기존 파이프라인을 실행하는 방법을 배우고자 하는 입문자를 위한 교육입니다. 기존 파이프라인을 이해하고 실행할 수 있도록 Nextflow 언어의 필수 개념만을 간단히 다루며, 명령줄 환경에서 Nextflow 파이프라인을 설정하고 실행하는 방법을 설명합니다. 또한, 커뮤니티에서 큐레이션한 다양한 파이프라인을 제공하는 nf-core 프로젝트와, 대규모 파이프라인 실행을 관리할 수 있도록 Seqera(Nextflow 개발사)에서 운영하는 플랫폼 등 Nextflow 생태계의 중요한 구성 요소들도 함께 소개합니다.

    [Nextflow Run 교육 시작하기 :material-arrow-right:](nextflow_run/index.md){ .md-button .md-button--primary }
-->

## 과학 분야를 위한 Nextflow

이 과정들은 위의 'Hello Nextflow'에서 소개한 개념과 구성 요소들을 실제 과학 분야의 활용 사례에 적용하는 방법을 보여줍니다. 각 과정은 단계적으로 실력을 쌓아갈 수 있도록 여러 개의 교육 모듈로 구성되어 있습니다.

!!! exercise "유전체학을 위한 Nextflow"

    !!! tip inline end ""

        :material-dna:{.nextflow-primary} Nextflow로 유전체학 파이프라인 개발 배우기

    이 과정은 자신만의 유전체학 파이프라인을 개발하고자 하는 연구자를 위한 교육입니다. 변이 분석(variant calling) 사례를 통해 간단하면서도 기능적인 유전체학 파이프라인을 개발하는 방법을 보여줍니다.

    [유전체학 교육 시작하기 :material-arrow-right:](nf4_science/genomics/){ .md-button .md-button--primary }

!!! exercise "RNAseq를 위한 Nextflow"

    !!! tip inline end ""

        :material-flask:{.nextflow-primary} Nextflow로 RNAseq 데이터 처리 파이프라인 개발 배우기

    이 과정은 자신만의 RNAseq 파이프라인을 개발하고자 하는 연구자를 위한 교육입니다. 벌크 RNAseq 처리 사례를 통해 간단하면서도 기능적인 RNAseq 파이프라인을 개발하는 방법을 소개합니다.

    [RNAseq 교육 시작하기 :material-arrow-right:](nf4_science/rnaseq/){ .md-button .md-button--primary }

어떤 분야나 사용 사례가 추가되었으면 하는지 [커뮤니티 포럼의 교육 섹션](https://community.seqera.io/c/training/)에 글을 남겨 알려주세요.

## Nextflow 심화 교육

이 과정들은 Nextflow의 기능을 더 자세히, 또는 더 고급 수준에서 활용하는 방법을 보여줍니다. 각 과정은 관련 주제에 대한 실력을 더욱 깊이 있게 다듬을 수 있도록 하나 이상의 교육 모듈로 구성되어 있습니다

!!! exercise "사이드 퀘스트"

    !!! tip inline end ""

        :material-compass:{.nextflow-primary} 다양한 주제를 다루는 교육 모듈

    이 과정은 Nextflow 개발자가 기술의 폭을 넓히거나 깊이를 더하고자 할 때 적합합니다. 모듈은 순서대로 제공되지만, 학습자는 자유롭게 원하는 주제부터 선택하여 학습할 수 있습니다. ‘Hello Nextflow’ 과정을 넘어서는 내용이 필요한 경우, 각 모듈 개요에 그에 대한 의존성이 명시되어 있습니다.

    [사이드 퀘스트 교육 시작하기 :material-arrow-right:](side_quests/){ .md-button .md-button--primary }

!!! exercise "기초 교육"

    !!! quote inline end ""

        :octicons-mortar-board-16:{.nextflow-primary} Nextflow의 모든 기능을 아우르는 포괄적인 교육 자료

    이 기초 교육 자료는 Nextflow의 전반적인 내용을 다루며, 복잡한 워크플로우를 구축하려는 분들을 위한 참고 자료로 활용될 수 있습니다.

    [기초 교육 시작하기 :material-arrow-right:](basic_training/index.md){ .md-button .md-button--primary }

!!! exercise "고급 교육"

    !!! quote inline end ""

        :fontawesome-solid-hat-wizard:{.nextflow-primary} Nextflow 마스터를 위한 고급 교육 자료

    Nextflow 언어와 런타임의 고급 기능을 탐색하며, 이를 통해 데이터 집약적인 워크플로우를 효율적이고 확장 가능하게 작성하는 방법을 배울 수 있습니다.

    [고급 교육 시작하기 :material-arrow-right:](advanced/index.md){ .md-button .md-button--primary }

## 기타 / 실험적 과정

이 교육 과정들은 현재 적극적으로 유지되거나 진행되지 않고 있으며, 향후 다른 용도로 재구성되거나 삭제될 수 있습니다. 해당 자료들은 교육 환경에서는 제공되지 않으며, GitHub 저장소에서 직접 다운로드하여 로컬에서 사용할 수 있습니다.

- **nf-customize** — nf-core 파이프라인 설정 ([문서](other/nf_customize) / [코드](https://github.com/nextflow-io/training/tree/master/other/nf-customize))
- **troubleshoot** — 문제 해결 연습 ([문서](other/troubleshoot) / [코드](https://github.com/nextflow-io/training/tree/master/other/troubleshoot))
- **hands-on (rnaseq)** — 벌크 RNAseq 파이프라인 개발 (사용 중단됨) ([문서](other/hands_on) / [코드](https://github.com/nextflow-io/training/tree/master/other/hands-on))

## 자료 모음

유용한 참고 링크들을 빠르게 확인해보세요:

| 참고 자료                                                   | 커뮤니티                                                     |
| ----------------------------------------------------------- | ------------------------------------------------------------ |
| [Nextflow 문서](https://nextflow.io/docs/latest/index.html) | [Nextflow Slack](https://www.nextflow.io/slack-invite.html)  |
| [Nextflow 홈페이지](https://nextflow.io/)                   | [nf-core](https://nf-co.re/)                                 |
| [Seqera](https://seqera.io/)                                | [Seqera 커뮤니티 포럼](https://community.seqera.io)          |

어디서부터 시작할지 모르겠다면 [도움 받기](help.md) 페이지를 참고하세요.

## 크레딧 및 기여

[![크리에이티브 커먼즈 라이선스: CC BY-NC-SA 4.0](assets/img/cc_by-nc-nd.svg){ align=right }](https://creativecommons.org/licenses/by-nc-nd/4.0/)

이 교육 자료는 [Seqera](https://seqera.io)에서 개발 및 유지 관리하고 있으며, 커뮤니티의 이익을 위해 오픈소스 라이선스([CC BY-NC-ND](https://creativecommons.org/licenses/by-nc-nd/4.0/)) 하에 배포됩니다. 라이선스 조건에 따라 자유롭게 재사용하실 수 있습니다. 교육을 진행하는 강사분이라면, 사용 경험과 개선점을 저희에게 알려주시면 감사하겠습니다.

커뮤니티의 수정 및 개선 제안을 언제든지 환영합니다. 각 페이지 오른쪽 상단의 :material-file-edit-outline: 아이콘을 클릭하면 GitHub에서 Pull Request를 통해 제안할 수 있습니다.

<div markdown class="homepage_logos">

![Seqera](assets/img/seqera_logo.png#only-light)

![Seqera](assets/img/seqera_logo_dark.png#only-dark)

</div>
