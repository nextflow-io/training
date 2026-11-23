---
title: 홈
description: Nextflow 커뮤니티 교육 포털에 오신 것을 환영합니다!
hide:
  - toc
  - footer
---

# Nextflow 교육

<div class="grid cards" markdown>

-   :material-book-open-variant:{ .lg .middle } __자기 주도 학습 과정__

    ---

    **Nextflow 커뮤니티 교육 포털에 오신 것을 환영합니다!**

    아래 과정을 자신의 속도에 맞춰 진행하세요. 웹 기반 환경 또는 개인 환경에서 학습할 수 있습니다.
    각 과정은 실습 중심으로 구성되어 있으며, 목표 지향적인 연습 문제를 독립적으로 완료할 수 있습니다.

    [과정 살펴보기 :material-arrow-down:](#catalog-of-nextflow-training-courses){ .md-button .md-button--primary .mt-1 }

-   :material-account-group-outline:{ .lg .middle } __교육 이벤트__

    ---

    **자기 주도 학습 이상의 것을 찾고 계신가요?**

    구조화된 교육 이벤트, 자체 교육 운영을 위한 가이드, 오픈 소스 라이선스 및 기여 정책을 확인하세요.

    [교육 이벤트 보기 :material-arrow-right:](training_events.md){ .md-button .md-button--secondary .mt-1 }

</div>

!!! note "AI 지원 번역"

    이 번역은 인공지능을 사용하여 생성되었으며 사람 번역자가 검토했습니다.
    피드백과 개선 제안을 환영합니다.
    자세한 내용은 [번역 가이드](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)를 참조하세요.

## Nextflow 교육 과정 카탈로그

<div class="grid cards" markdown>

-   :material-account:{ .lg .middle } __사용자를 위한 과정__

    ---

    ### :material-play-circle:{.nextflow-primary} 파이프라인 실행 {.mt-1}

    코드 작성 없이 기존 파이프라인을 실행하는 방법을 학습합니다.

    ??? courses "**Nextflow Run:** Nextflow로 파이프라인 실행하기"

        코드 이해 없이 Nextflow 파이프라인을 실행하는 방법을 빠르게 학습합니다. 파이프라인 실행, 출력 결과 가져오기, 컨테이너 사용, 기본 수준의 실행 설정을 다룹니다.

        [교육 보기 :material-arrow-right:](nextflow_run/index.md){ .md-button .md-button--secondary }

    ??? courses "**Use nf-core:** 커뮤니티 큐레이션 파이프라인 찾기 및 실행하기"

        nf-core 커뮤니티 프로젝트에서 파이프라인을 찾고, 실행하고, 설정하는 방법을 빠르게 학습합니다. 최소한의 데모 파이프라인부터 시작하여 프로덕션 규모의 분석 파이프라인까지 단계적으로 확장합니다.

        [교육 보기 :material-arrow-right:](nfcore_use/index.md){ .md-button .md-button--secondary }

    ??? courses "**Scale with Seqera:** 대규모 파이프라인 실행 및 모니터링"

        웹 인터페이스와 명령줄 모두에서 Seqera Platform을 사용하여 Nextflow 파이프라인을 실행하고 모니터링하는 방법을 실습합니다.

        [교육 보기 :material-arrow-right:](seqera_scale/index.md){ .md-button .md-button--secondary }

    ---

    ### :material-tune:{.nextflow-primary} 실행 관리 {.mt-1}

    파이프라인 실행을 효과적으로 관리하는 방법을 학습합니다.

    ??? courses "**Configure Execution:** 리소스, 재시도, 실행 프로파일 설정하기"

        Nextflow 파이프라인 실행 설정을 실습합니다. 다양한 컴퓨팅 환경에 적응하고, 리소스 할당 및 재시도를 제어하며, 사전 설정된 설정 프로파일 간에 전환하는 방법을 다룹니다.

        [교육 보기 :material-arrow-right:](config_exec/index.md){ .md-button .md-button--secondary }

    !!! info compact "추가 주제 예정"

        성능 튜닝, HPC/클라우드 실행 등이 이 섹션에 추가될 예정입니다.
        다음에 다룰 주제에 대해 [간단한 관심 설문](https://seqera.typeform.com/to/JCs91e8v)에서 투표해 주세요.

-   :material-code-tags:{ .lg .middle } __개발자를 위한 과정__

    ---

    ### :material-wrench:{.nextflow-primary} 파이프라인 작성 {.mt-1}

    나만의 Nextflow 파이프라인을 개발하는 방법을 학습합니다.

    ??? courses "**Hello Nextflow:** 처음부터 나만의 파이프라인 개발하기"

        이 과정은 간단하지만 완전히 기능하는 파이프라인을 개발할 수 있을 만큼 Nextflow 언어의 핵심 구성 요소를 충분히 다룹니다. 파이프라인 설계, 개발 및 설정 방법의 핵심 요소도 포함됩니다.

        [교육 보기 :material-arrow-right:](hello_nextflow/index.md){ .md-button .md-button--secondary }

    ??? courses "**Build with nf-core:** nf-core 도구 및 규칙 활용하기"

        [nf-core](https://nf-co.re/) 호환 파이프라인을 개발하고자 하는 Nextflow 개발자를 위한 과정입니다.
        nf-core 파이프라인의 구조를 충분히 다루어 nf-core 템플릿과 개발 모범 사례를 활용하고 기존 nf-core 모듈을 사용하는 간단하지만 완전히 기능하는 파이프라인을 개발할 수 있도록 합니다.

        [교육 보기 :material-arrow-right:](nfcore_build/index.md){ .md-button .md-button--secondary }

    ??? catalog "**Side Quests:** 고급 Nextflow 주제 살펴보기"

        Nextflow 개발자가 특정 주제에 대한 역량을 넓히거나 심화하고자 할 때 활용할 수 있는 단독 실행형 단기 과정 모음입니다.
        순서대로 제시되어 있지만 어떤 순서로든 수강할 수 있습니다(각 단기 과정 개요의 의존성 참조).

        [Side Quests 살펴보기 :material-arrow-right:](side_quests/index.md){ .md-button .md-button--secondary }

    ---

    ### :material-microscope:{.nextflow-primary} 과학 연구를 위한 Nextflow {.mt-1}

    특정 과학 응용 분야를 위한 Nextflow 파이프라인 개발 방법을 학습합니다.

    ??? courses "**Genomics:** 변이 호출 파이프라인 개발하기"

        자체 유전체학 파이프라인을 개발하고자 하는 연구자를 위한 과정입니다. 변이 호출 사례를 활용하여 필수 Nextflow 개발 패턴을 학습합니다.

        [교육 보기 :material-arrow-right:](nf4_science/genomics/index.md){ .md-button .md-button--secondary }

    ??? courses "**RNAseq:** 벌크 RNAseq 처리 파이프라인 개발하기"

        자체 RNAseq 파이프라인을 개발하고자 하는 연구자를 위한 과정입니다. 벌크 RNAseq 처리 사례를 활용하여 필수 Nextflow 개발 패턴을 학습합니다.

        [교육 보기 :material-arrow-right:](nf4_science/rnaseq/index.md){ .md-button .md-button--secondary }

    ??? courses "**Bioimaging:** 이미징 파이프라인 실행 및 설정하기"

        바이오이미징 파이프라인을 실행하고 설정하는 방법을 학습하고자 하는 연구자를 위한 과정입니다. nf-core/molkart를 활용하여 필수 Nextflow 사용 패턴을 학습합니다.

        [교육 보기 :material-arrow-right:](nf4_science/imaging/index.md){ .md-button .md-button--secondary }

</div>

## 설정 및 도움말

<div class="grid cards mb-4" markdown>

-   :material-cog-outline:{ .lg .middle } __교육 환경__

    ---

    Nextflow 교육을 위한 환경 설정 옵션을 안내합니다.

    [교육 환경 보기 :material-arrow-right:](envsetup/index.md){ .md-button .md-button--secondary }

-   :material-tag-outline:{ .lg .middle } __Nextflow 버전__

    ---

    Nextflow 문법 버전의 변화를 이해하고 관리하는 방법을 안내합니다.

    [버전 요구사항 확인 :material-arrow-right:](info/nxf_versions.md){ .md-button .md-button--secondary }

-   :material-file-code-outline:{ .lg .middle } __Hello 파이프라인__

    ---

    Hello 파이프라인의 기능과 구조에 대한 요약을 확인합니다.

    [요약 읽기 :material-arrow-right:](info/hello_pipeline.md){ .md-button .md-button--secondary }

-   :material-lifebuoy:{ .lg .middle } __도움말__

    ---

    Nextflow 교육 중 문제가 발생했을 때 활용할 수 있는 유용한 자료를 안내합니다.

    [도움말 찾기 :material-arrow-right:](help.md){ .md-button .md-button--secondary }

</div>
