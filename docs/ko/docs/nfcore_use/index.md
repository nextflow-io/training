---
title: nf-core 사용하기
hide:
  - toc
page_type: index_page
index_type: course
additional_information:
  technical_requirements: true
  learning_objectives:
    - nf-core 커뮤니티 파이프라인을 찾고, 가져오고, 실행합니다
    - 매개변수와 설정 파일을 사용하여 파이프라인 실행을 설정합니다
    - nf-core 파이프라인이 매개변수와 입력 데이터를 검증하는 방법을 이해합니다
    - 프로덕션 규모의 파이프라인(nf-core/rnaseq)을 실행하고 기본 리소스 할당을 재정의합니다
  audience_prerequisites:
    - "**대상:** 이 과정은 로컬 Nextflow 파이프라인을 실행할 줄 알고 nf-core를 처음 접하는 학습자를 대상으로 하며, 기존 커뮤니티 파이프라인을 실행하고자 하는 분들을 위해 설계되었습니다."
    - "**기술:** 명령줄, 기본 스크립팅 개념 및 일반적인 파일 형식에 대한 기본적인 이해가 필요합니다."
    - "**과정:** [Nextflow Run](../nextflow_run/index.md)을 완료했거나 `nextflow run`으로 로컬 파이프라인을 실행하는 데 익숙해야 합니다."
    - "**도메인:** 실습에서 바이오인포매틱스 파이프라인을 사용하지만, 사전 과학 도메인 지식은 필요하지 않습니다."
---

# nf-core 사용하기

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } AI 지원 번역 - [자세히 알아보기 및 개선 제안](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

**nf-core 사용하기는 nf-core 커뮤니티 파이프라인을 찾고, 실행하고, 설정하는 방법을 실습 중심으로 학습하는 과정입니다.**

실용적인 예제와 유도형 실습을 통해 nf-core 파이프라인을 찾고 가져오는 방법, 내장된 테스트 프로파일을 사용하여 실행하는 방법, 매개변수와 설정 파일을 통해 실행을 맞춤화하는 방법을 학습합니다.

이 과정을 마치면 직접 분석에 nf-core 파이프라인을 활용할 수 있는 역량과 자신감을 갖추게 됩니다.

<!-- additional_information -->

## 과정 개요

이 과정은 실습 중심으로 구성되어 있으며, 목표 지향적인 실습을 통해 내용을 단계적으로 학습합니다.

교육 목적으로 nf-core 프로젝트가 관리하는 최소한의 파이프라인인 `nf-core/demo`부터 시작하여, 학습한 내용을 대규모 RNA 시퀀싱 분석에 널리 사용되는 프로덕션 파이프라인인 `nf-core/rnaseq`에 적용합니다.

이 과정은 파이프라인 실행에 초점을 맞춥니다.
nf-core 호환 파이프라인 개발 입문을 원하신다면 [Build with nf-core](../nfcore_build/index.md)를 참조하세요.

### 학습 계획

| 과정 챕터                                                                | 요약                                                                          | 예상 소요 시간 |
| ------------------------------------------------------------------------ | ----------------------------------------------------------------------------- | -------------- |
| [파트 1: 데모 파이프라인 실행](./01_run_demo.md)                         | nf-core 파이프라인을 찾고 가져온 후 테스트 프로파일을 사용하여 실행합니다    | 20분           |
| [파트 2: 파이프라인 실행 설정](./02_configure_execution.md)              | 매개변수를 설정하고, 검증을 이해하며, 리소스 할당과 도구 인자를 맞춤화합니다 | 20분           |
| [파트 3: 프로덕션 파이프라인 실행](./03_run_production_pipeline.md)      | nf-core/rnaseq를 가져와 실행하고 기본 리소스 할당을 재정의합니다             | 20분           |

이 과정을 마치면 nf-core 프로젝트가 제공하는 풍부한 커뮤니티 파이프라인을 충분히 활용할 수 있게 됩니다.

과정을 시작할 준비가 되셨나요?

[학습 시작하기 :material-arrow-right:](00_orientation.md){ .md-button .md-button--primary }
