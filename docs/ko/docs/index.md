---
title: 홈
description: Nextflow 커뮤니티 교육 포털에 오신 것을 환영합니다!
hide:
  - toc
  - footer
---

# Nextflow 교육

<div class="grid cards" markdown>

-   :material-book-open-variant:{ .lg .middle } __셀프 서비스 과정__

    ---

    **Nextflow 커뮤니티 교육 포털에 오신 것을 환영합니다!**

    아래에 나열된 교육 과정은 셀프 서비스 방식으로 활용할 수 있도록 설계되었습니다.
    GitHub Codespaces를 통해 제공되는 웹 기반 환경이나 직접 구성한 환경에서 언제든지 학습할 수 있습니다.

    [과정 살펴보기 :material-arrow-right:](#catalog-of-nextflow-training-courses){ .md-button .md-button--primary .mt-1 }

-   :material-information-outline:{ .lg .middle } __추가 정보__

    ---

    ??? warning "버전 호환성"

        <!-- 이 내용을 업데이트할 경우 로컬 설치 페이지에도 동일하게 반영해야 합니다 -->
        **2026년 1월부터 별도로 명시되지 않는 한, 모든 Nextflow 교육 과정은 strict syntax가 활성화된 Nextflow 버전 25.10.2 이상을 필요로 합니다.**

        버전 요구 사항 및 strict syntax에 대한 자세한 내용은 [Nextflow 문서 마이그레이션 가이드](https://nextflow.io/docs/latest/strict-syntax.html)를 참조하세요.

        이전 syntax에 해당하는 구버전 교육 자료는 이 웹페이지 메뉴 바의 버전 선택기를 통해 확인할 수 있습니다.

    ??? terminal "환경 옵션"

        교육에 필요한 모든 것이 사전 설치된 웹 기반 교육 환경을 제공합니다. GitHub Codespaces를 통해 이용할 수 있습니다(무료 GitHub 계정 필요).

        [![GitHub Codespaces에서 열기](https://github.com/codespaces/badge.svg)](https://codespaces.new/nextflow-io/training?quickstart=1&ref=master)

        이 방법이 적합하지 않은 경우, 다른 [환경 옵션](./envsetup/index.md)을 확인하세요.

    ??? learning "교육 이벤트"

        구조화된 이벤트 형태로 Nextflow 교육을 받고 싶다면 다양한 기회가 있습니다. 다음 옵션을 확인해 보세요:

        - 커뮤니티 팀이 분기별로 주관하는 **[Training Weeks]()**
        - Seqera가 주관하는 대면 교육 이벤트를 포함한 **[Seqera Events](https://seqera.io/events/)** ('Seqera Sessions' 및 'Nextflow Summit' 검색)
        - 지역 커뮤니티를 위한 이벤트를 주관하는 **[Nextflow Ambassadors]()**
        - 커뮤니티 해커톤을 포함한 **[nf-core events](https://nf-co.re/events)**

    ??? people "강사를 위한 정보"

        자체 교육을 진행하는 강사라면 출처를 적절히 표기하는 한 교육 포털의 자료를 직접 사용하실 수 있습니다. 자세한 내용은 아래 '크레딧 및 기여' 항목을 참조하세요.

        또한 교육 활동을 더 잘 지원할 수 있는 방법에 대한 의견을 환영합니다! [community@seqera.io](mailto:community@seqera.io) 또는 커뮤니티 포럼([도움말](help.md) 페이지 참조)을 통해 연락해 주세요.

    ??? licensing "오픈 소스 라이선스 및 기여 정책"

        [![Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International (CC BY-NC-SA 4.0)](assets/img/cc_by-nc-sa.svg){ align=right }](https://creativecommons.org/licenses/by-nc-sa/4.0/)

        이 교육 자료는 [Seqera](https://seqera.io)가 개발 및 유지 관리하며, 커뮤니티의 이익을 위해 오픈 소스 라이선스([CC BY-NC-SA](https://creativecommons.org/licenses/by-nc-sa/4.0/)) 하에 공개되었습니다. 라이선스 범위를 벗어나는 방식으로 이 자료를 사용하고자 하는 경우(상업적 이용 및 재배포 제한 사항 참고), [community@seqera.io](mailto:community@seqera.io)로 문의해 주세요.

        커뮤니티의 개선 사항, 수정 및 버그 리포트를 환영합니다. 각 페이지 오른쪽 상단의 :material-file-edit-outline: 아이콘을 클릭하면 코드 저장소로 이동하며, 이슈를 보고하거나 pull request를 통해 교육 자료 변경을 제안할 수 있습니다. 자세한 내용은 저장소의 `README.md`를 참조하세요.

</div>

!!! note "AI 지원 번역"

    이 번역은 인공지능을 사용하여 생성되었으며 사람 번역자가 검토했습니다.
    피드백과 개선 제안을 환영합니다.
    자세한 내용은 [번역 가이드](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)를 참조하세요.

## Nextflow 교육 과정 카탈로그

<div class="grid cards" markdown>

-   :material-walk:{ .lg .middle } __입문 트랙__

    ---

    ### :material-compass:{.nextflow-primary} Nextflow 입문자를 위한 과정 {.mt-1}

    Nextflow를 처음 접하는 분들을 위한 도메인 비특화 과정입니다. 각 과정은 학습자가 단계적으로 역량을 쌓을 수 있도록 설계된 일련의 교육 모듈로 구성됩니다.

    ??? courses "**Hello Nextflow:** 나만의 파이프라인 개발 학습"

        이 과정은 간단하지만 완전히 기능하는 파이프라인을 개발할 수 있을 만큼 Nextflow 언어의 핵심 구성 요소를 충분히 다루며, 파이프라인 설계, 개발 및 설정 방법의 핵심 요소도 포함합니다.

        [Hello Nextflow 교육 시작하기 :material-arrow-right:](hello_nextflow/index.md){ .md-button .md-button--secondary }

    ??? courses "**Nextflow Run:** 기존 파이프라인 실행 학습"

        Nextflow 파이프라인 실행 및 설정에 대한 간결한 입문 과정으로, Hello Nextflow 개발자 과정을 기반으로 하되 코드에 대한 비중이 적습니다. 실행, 출력, 기본 코드 구조, 다양한 컴퓨팅 환경을 위한 설정을 다룹니다.

        [Nextflow Run 교육 시작하기 :material-arrow-right:](nextflow_run/index.md){ .md-button .md-button--secondary }

    ---

    ### :material-microscope:{.nextflow-primary} 과학 연구를 위한 Nextflow {.mt-1}

    'Hello Nextflow'에서 소개된 개념과 구성 요소를 특정 과학적 활용 사례에 적용하는 방법을 학습합니다.

    ??? courses "**Nextflow for Genomics** (변이 분석)"

        자체 유전체 파이프라인 개발을 원하는 연구자를 위한 과정입니다. 변이 분석 사례를 활용하여 간단하지만 기능적인 유전체 파이프라인을 개발하는 방법을 학습합니다.

        [Nextflow for Genomics 교육 시작하기 :material-arrow-right:](nf4_science/genomics/){ .md-button .md-button--secondary }

    ??? courses "**Nextflow for RNAseq** (bulk RNAseq)"

        자체 RNAseq 파이프라인 개발을 원하는 연구자를 위한 과정입니다. bulk RNAseq 처리 사례를 활용하여 간단하지만 기능적인 RNAseq 파이프라인을 개발하는 방법을 학습합니다.

        [Nextflow for RNAseq 교육 시작하기 :material-arrow-right:](nf4_science/rnaseq/){ .md-button .md-button--secondary }

    ??? courses "**Nextflow for Imaging** (공간 오믹스)"

        이미징 및 공간 오믹스 분야 연구자 중 분석 파이프라인 실행 및 맞춤화를 원하는 분들을 위한 과정입니다. nf-core/molkart 파이프라인을 활용하여 생물학적으로 관련성 있는 파이프라인을 통해 Nextflow 파이프라인 워크플로우의 실행, 설정 및 입력 관리 방법을 학습합니다.

        [Nextflow for Imaging 교육 시작하기 :material-arrow-right:](nf4_science/imaging/){ .md-button .md-button--secondary }

-   :material-run:{ .lg .middle } __심화 트랙__

    ---

    ### :material-bridge:{.nextflow-primary} Nextflow에서 nf-core로 {.mt-1}

    [nf-core](https://nf-co.re/) 커뮤니티 프로젝트의 코드와 모범 사례를 활용하는 방법을 학습합니다.

    이 과정들은 Nextflow 기초부터 nf-core 모범 사례까지 단계적으로 학습할 수 있도록 구성되어 있습니다.
    nf-core 커뮤니티가 파이프라인을 구축하는 방법과 이유, 그리고 이러한 기법에 기여하고 재사용하는 방법을 학습합니다.

    ??? courses "**Hello nf-core:** nf-core 시작하기"

        [nf-core](https://nf-co.re/) 호환 파이프라인을 실행하고 개발하고자 하는 개발자를 위한 과정입니다. nf-core 템플릿과 개발 모범 사례를 따르는 간단하지만 완전히 기능하는 파이프라인을 개발하고, 기존 nf-core 모듈을 활용할 수 있을 만큼 nf-core 파이프라인 구조를 충분히 다룹니다.

        [Hello nf-core 교육 시작하기 :material-arrow-right:](hello_nf-core/index.md){ .md-button .md-button--secondary }

    ---

    ### :material-rocket-launch:{.nextflow-primary} Nextflow 심화 교육 {.mt-1}

    실제 활용 사례를 해결하기 위한 Nextflow 파이프라인 개발 및 배포의 고급 개념과 메커니즘을 학습합니다.

    ??? courses "**Side Quests:** 단독 주제 심층 학습"

        특정 주제에 대한 역량을 넓히거나 심화하고자 하는 Nextflow 개발자를 위한 단독 실행형 단기 과정입니다. 순서대로 제시되지만 어떤 순서로든 수강할 수 있습니다(각 단기 과정 개요의 의존성 참조).

        [Side Quests 살펴보기 :material-arrow-right:](side_quests/){ .md-button .md-button--secondary }

    ??? courses "**Training Collections:** Side Quests를 통한 추천 학습 경로"

        Training Collections는 특정 주제나 활용 사례를 중심으로 포괄적인 학습 경험을 제공하기 위해 여러 Side Quests를 결합합니다.

        [Training Collections 살펴보기 :material-arrow-right:](training_collections/){ .md-button .md-button--secondary }

</div>

!!! info "아카이브된 교육 자료를 찾고 계신가요?"

    구버전 교육 자료(Fundamentals Training, Advanced Training 및 기타 실험적 과정)는 Nextflow 3.0 strict syntax와 호환되지 않아 교육 포털에서 제거되었습니다.
    해당 자료가 필요한 경우 2026년 1월 이전의 [git history](https://github.com/nextflow-io/training)에서 확인할 수 있습니다.

---

<div markdown class="homepage_logos">

![Seqera](assets/img/seqera_logo.png#only-light)

![Seqera](assets/img/seqera_logo_dark.png#only-dark)

</div>
