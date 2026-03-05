---
title: 홈
description: Nextflow 커뮤니티 교육 포털에 오신 것을 환영합니다!
hide:
  - toc
  - footer
---

# Nextflow 교육

<div class="grid cards" markdown>

-   :material-book-open-variant:{ .lg .middle } __자율 학습 과정__

    ---

    **Nextflow 커뮤니티 교육 포털에 오신 것을 환영합니다!**

    아래 나열된 교육 과정은 자율 학습 자료로 사용할 수 있도록 설계되었습니다.
    Github Codespaces를 통해 제공되는 웹 기반 환경이나 자체 환경에서 언제든지 원하는 시간에 학습할 수 있습니다.

    [과정 둘러보기 :material-arrow-right:](#catalog-of-nextflow-training-courses){ .md-button .md-button--primary .mt-1 }

-   :material-information-outline:{ .lg .middle } __추가 정보__

    ---

    ??? warning "버전 호환성"

        <!-- Any update to this content needs to be copied to the local installation page -->
        **2026년 1월부터 별도로 명시되지 않는 한, 모든 Nextflow 교육 과정은 엄격한 구문이 활성화된 Nextflow 버전 25.10.2 이상을 필요로 합니다.**

        버전 요구 사항 및 엄격한 구문에 대한 자세한 내용은 [Nextflow 문서 마이그레이션 가이드](https://nextflow.io/docs/latest/strict-syntax.html)를 참조하세요.

        이전 구문에 해당하는 구버전 교육 자료는 이 웹페이지 메뉴 바의 버전 선택기를 통해 이용할 수 있습니다.

    ??? terminal "환경 옵션"

        교육에 필요한 모든 것이 사전 설치된 웹 기반 교육 환경을 Github Codespaces를 통해 제공합니다(무료 GitHub 계정 필요).

        [![Open in GitHub Codespaces](https://github.com/codespaces/badge.svg)](https://codespaces.new/nextflow-io/training?quickstart=1&ref=master)

        이 방법이 적합하지 않은 경우, 다른 [환경 옵션](./envsetup/index.md)을 참조하세요.

    ??? learning "교육 이벤트"

        구조화된 이벤트의 일부로 Nextflow 교육을 받고 싶으시다면, 다양한 기회가 있습니다. 다음 옵션을 확인하시기 바랍니다:

        - 커뮤니티 팀이 분기별로 주최하는 **[Training Weeks]()**
        - **[Seqera Events](https://seqera.io/events/)**에는 Seqera가 주최하는 대면 교육 이벤트가 포함됩니다('Seqera Sessions' 및 'Nextflow Summit' 검색)
        - **[Nextflow Ambassadors]()**는 지역 커뮤니티를 위한 이벤트를 주최합니다
        - **[nf-core events](https://nf-co.re/events)**에는 커뮤니티 해커톤이 포함됩니다

    ??? people "강사를 위한 정보"

        자체 교육을 진행하는 강사라면, 적절한 출처를 표시하는 한 교육 포털에서 직접 자료를 사용하실 수 있습니다. 자세한 내용은 아래 '크레딧 및 기여'를 참조하세요.

        또한, 교육 활동을 더 잘 지원할 수 있는 방법에 대한 의견을 듣고 싶습니다! [community@seqera.io](mailto:community@seqera.io) 또는 커뮤니티 포럼([도움말](help.md) 페이지 참조)으로 연락 주시기 바랍니다.

    ??? licensing "오픈소스 라이선스 및 기여 정책"

        [![Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International (CC BY-NC-SA 4.0)](assets/img/cc_by-nc-sa.svg){ align=right }](https://creativecommons.org/licenses/by-nc-sa/4.0/)

        이 교육 자료는 [Seqera](https://seqera.io)에서 개발 및 유지 관리하며, 커뮤니티의 이익을 위해 오픈소스 라이선스([CC BY-NC-SA](https://creativecommons.org/licenses/by-nc-sa/4.0/))로 공개됩니다. 이 자료를 라이선스 범위를 벗어나는 방식으로 사용하고자 하는 경우(상업적 사용 및 재배포에 대한 제한 사항 참고), [community@seqera.io](mailto:community@seqera.io)로 연락하여 요청 사항을 논의해 주시기 바랍니다.

        커뮤니티의 개선 사항, 수정 및 버그 리포트를 환영합니다. 모든 페이지의 오른쪽 상단에는 코드 저장소로 연결되는 :material-file-edit-outline: 아이콘이 있으며, 여기에서 문제를 보고하거나 풀 리퀘스트를 통해 교육 소스 자료에 대한 변경 사항을 제안할 수 있습니다. 자세한 내용은 저장소의 `README.md`를 참조하세요.

</div>

!!! note "AI 지원 번역"

    이 번역은 인공지능을 사용하여 생성되었으며 사람 번역자가 검토했습니다.
    피드백과 개선 제안을 환영합니다.
    자세한 내용은 [번역 가이드](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)를 참조하세요.

## Nextflow 교육 과정 목록

<div class="grid cards" markdown>

-   :material-walk:{ .lg .middle } __입문 트랙__

    ---

    ### :material-compass:{.nextflow-primary} Nextflow for Newcomers {.mt-1}

    Nextflow를 처음 접하는 분들을 위한 도메인 독립적 과정입니다. 각 과정은 학습자가 단계적으로 기술을 쌓을 수 있도록 설계된 일련의 교육 모듈로 구성됩니다.

    ??? courses "**Hello Nextflow:** 자체 파이프라인 개발 학습"

        이 과정은 간단하지만 완전히 기능하는 파이프라인을 개발할 수 있을 만큼 충분히 상세하게 Nextflow 언어의 핵심 구성 요소를 다루며, 파이프라인 설계, 개발 및 구성 관행의 주요 요소를 포함합니다.

        [Hello Nextflow 교육 시작하기 :material-arrow-right:](hello_nextflow/index.md){ .md-button .md-button--secondary }

    ??? courses "**Nextflow Run:** 기존 파이프라인 실행 학습"

        Hello Nextflow 개발자 과정을 기반으로 하되 코드에 대한 집중도를 낮춘, Nextflow 파이프라인 실행 및 구성에 대한 간결한 소개입니다. 실행, 출력, 기본 코드 구조 및 다양한 컴퓨팅 환경에 대한 구성을 다룹니다.

        [Nextflow Run 교육 시작하기 :material-arrow-right:](nextflow_run/index.md){ .md-button .md-button--secondary }

    ---

    ### :material-microscope:{.nextflow-primary} Nextflow for Science {.mt-1}

    'Hello Nextflow'에서 제시된 개념과 구성 요소를 특정 과학적 사용 사례에 적용하는 방법을 학습합니다.

    ??? courses "**Nextflow for Genomics** (변이 검출)"

        자체 유전체학 파이프라인을 개발하고자 하는 연구자를 위한 과정입니다. 변이 검출 사용 사례를 사용하여 간단하지만 기능적인 유전체학 파이프라인을 개발하는 방법을 보여줍니다.

        [Nextflow for Genomics 교육 시작하기 :material-arrow-right:](nf4_science/genomics/){ .md-button .md-button--secondary }

    ??? courses "**Nextflow for RNAseq** (bulk RNAseq)"

        자체 RNAseq 파이프라인을 개발하고자 하는 연구자를 위한 과정입니다. bulk RNAseq 처리 사용 사례를 사용하여 간단하지만 기능적인 RNAseq 파이프라인을 개발하는 방법을 보여줍니다.

        [Nextflow for RNAseq 교육 시작하기 :material-arrow-right:](nf4_science/rnaseq/){ .md-button .md-button--secondary }

    ??? courses "**Nextflow for Imaging** (공간 오믹스)"

        분석 파이프라인을 실행하고 사용자 정의하는 방법을 배우고자 하는 이미징 및 공간 오믹스 연구자를 위한 과정입니다. nf-core/molkart 파이프라인을 사용하여 생물학적으로 관련된 파이프라인을 제공하고, Nextflow 파이프라인 워크플로우를 실행, 구성 및 입력 관리하는 방법을 보여줍니다.

        [Nextflow for Imaging 교육 시작하기 :material-arrow-right:](nf4_science/imaging/){ .md-button .md-button--secondary }

-   :material-run:{ .lg .middle } __고급 트랙__

    ---

    ### :material-bridge:{.nextflow-primary} From Nextflow to nf-core {.mt-1}

    [nf-core](https://nf-co.re/) 커뮤니티 프로젝트의 코드와 모범 사례를 활용하는 방법을 학습합니다.

    이 과정은 Nextflow 기초부터 nf-core 모범 사례까지 학습할 수 있도록 도와줍니다.
    nf-core 커뮤니티가 파이프라인을 구축하는 방법과 이유를 이해하고, 기여하고 이러한 기술을 재사용하는 방법을 배웁니다.

    ??? courses "**Hello nf-core:** nf-core 시작하기"

        [nf-core](https://nf-co.re/) 준수 파이프라인을 실행하고 개발하고자 하는 개발자를 위한 과정입니다. nf-core 템플릿과 개발 모범 사례를 따르는 간단하지만 완전히 기능하는 파이프라인을 개발할 수 있을 만큼 충분히 상세하게 nf-core 파이프라인의 구조를 다루며, 기존 nf-core 모듈을 사용하는 방법도 포함합니다.

        [Hello nf-core 교육 시작하기 :material-arrow-right:](hello_nf-core/index.md){ .md-button .md-button--secondary }

    ---

    ### :material-rocket-launch:{.nextflow-primary} Advanced Nextflow Training {.mt-1}

    실제 사용 사례를 해결하기 위한 Nextflow 파이프라인 개발 및 배포를 위한 고급 개념과 메커니즘을 학습합니다.

    ??? courses "**Side Quests:** 독립적인 주제에 대한 심화 학습"

        특정 주제에 대한 범위를 넓히거나 기술을 심화하고자 하는 Nextflow 개발자를 위한 독립적인 단기 과정입니다. 선형적으로 제시되지만 원하는 순서대로 수강할 수 있습니다(각 단기 과정 개요의 의존성 참조).

        [Side Quests 둘러보기 :material-arrow-right:](side_quests/){ .md-button .md-button--secondary }

    ??? courses "**Training Collections:** Side Quests를 통한 권장 학습 경로"

        Training Collections는 특정 주제나 사용 사례를 중심으로 포괄적인 학습 경험을 제공하기 위해 여러 Side Quests를 결합합니다.

        [Training Collections 둘러보기 :material-arrow-right:](training_collections/){ .md-button .md-button--secondary }

</div>

!!! info "아카이브된 교육 자료를 찾고 계신가요?"

    구버전 교육 자료(Fundamentals Training, Advanced Training 및 기타 실험적 과정)는 Nextflow 3.0 엄격한 구문과 호환되지 않아 교육 포털에서 제거되었습니다.
    이러한 자료가 필요한 경우, 2026년 1월 이전의 [git 히스토리](https://github.com/nextflow-io/training)에서 이용할 수 있습니다.

---

<div markdown class="homepage_logos">

![Seqera](assets/img/seqera_logo.png#only-light)

![Seqera](assets/img/seqera_logo_dark.png#only-dark)

</div>
