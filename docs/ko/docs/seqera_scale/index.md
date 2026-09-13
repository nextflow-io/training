---
title: Scale with Seqera
hide:
  - toc
page_type: index_page
index_type: course
additional_information:
  technical_requirements: true
  learning_objectives:
    - Sign up for Seqera Platform and explore the Community Showcase
    - Add a pipeline to a workspace and launch it from the web interface
    - Authenticate and launch pipelines from the command line with the `tw` CLI
    - Register a GitHub-hosted pipeline and launch it both ways
  audience_prerequisites:
    - "**Audience:** This course is designed for learners who want to run Nextflow pipelines at scale using Seqera Platform."
    - "**Skills:** Familiarity with running nf-core pipelines from the command line is assumed."
    - "**Courses:** Must have completed [Nextflow Run](../nextflow_run/index.md) and [Use nf-core](../nfcore_use/index.md), or otherwise be comfortable running local and `nf-core/rnaseq` pipelines."
---

---
title: Seqera로 확장하기
hide:
  - toc
page_type: index_page
index_type: course
additional_information:
  technical_requirements: true
  learning_objectives:
    - Seqera Platform에 가입하고 Community Showcase를 살펴봅니다
    - 워크스페이스에 파이프라인을 추가하고 웹 인터페이스에서 실행합니다
    - `tw` CLI로 인증하고 명령줄에서 파이프라인을 실행합니다
    - GitHub에 호스팅된 파이프라인을 등록하고 두 가지 방법으로 실행합니다
  audience_prerequisites:
    - "**대상:** 이 과정은 Seqera Platform을 사용하여 Nextflow 파이프라인을 대규모로 실행하려는 학습자를 위해 설계되었습니다."
    - "**기술:** 명령줄에서 nf-core 파이프라인을 실행한 경험이 있다고 가정합니다."
    - "**과정:** [Nextflow Run](../nextflow_run/index.md) 및 [Use nf-core](../nfcore_use/index.md)를 완료했거나, 로컬 및 `nf-core/rnaseq` 파이프라인 실행에 익숙해야 합니다."
---

# Seqera로 확장하기

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } AI 지원 번역 - [자세히 알아보기 및 개선 제안](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

**Scale with Seqera는 Seqera Platform을 사용하여 Nextflow 파이프라인을 실행하고 모니터링하는 실습 입문 과정입니다.**

실용적인 예제를 통해 Seqera Platform 접근 권한을 설정하고, 웹 인터페이스와 명령줄 모두에서 프로덕션 규모의 파이프라인을 실행하며, 워크스페이스에 새 파이프라인을 추가합니다.

이 과정을 마치면 Seqera Platform에서 자신의 파이프라인을 실행하고 모니터링할 수 있는 역량과 자신감을 갖추게 됩니다.

<!-- additional_information -->

## 과정 개요

이 과정은 실습 중심으로 구성되어 있으며, [Use nf-core](../nfcore_use/index.md)에서 이미 실행한 파이프라인을 기반으로 합니다.

먼저 Seqera Platform에 가입하고 웹 인터페이스에서 프로덕션 규모의 파이프라인인 `nf-core/rnaseq`를 실행합니다.
이후 `tw` 명령줄 도구로 전환하여 터미널에서 동일한 작업을 수행하고, 마지막으로 새 파이프라인인 `nf-core/demo`를 등록하여 두 가지 방법으로 실행합니다.

### 학습 계획

| 과정 챕터                                                                  | 요약                                                                                          | 예상 소요 시간 |
| -------------------------------------------------------------------------- | --------------------------------------------------------------------------------------------- | -------------- |
| [파트 1: 웹 인터페이스에서 파이프라인 실행](./01_run_with_seqera.md)       | Seqera Platform 접근 권한을 설정하고 웹 인터페이스에서 프로덕션 규모의 파이프라인을 실행합니다 | 20분           |
| [파트 2: 명령줄에서 파이프라인 실행](./02_launch_from_cli.md)              | `tw` CLI를 인증하고, 저장된 파이프라인을 실행하며, CLI에서 새 파이프라인을 등록합니다          | 25분           |

이 과정을 마치면 웹 인터페이스와 명령줄 중 어느 방법을 선호하든 Seqera Platform에서 Nextflow 파이프라인을 실행하고 모니터링하는 데 익숙해집니다.

과정을 시작할 준비가 되셨나요?

[학습 시작하기 :material-arrow-right:](00_orientation.md){ .md-button .md-button--primary }
