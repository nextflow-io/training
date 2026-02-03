---
title: Hello nf-core
hide:
  - toc
page_type: index_page
index_type: course
additional_information:
  technical_requirements: true
  learning_objectives:
    - nf-core 파이프라인의 실행을 검색하고, 시작하고, 관리합니다
    - nf-core 파이프라인의 코드 구조와 프로젝트 구성을 설명합니다
    - 템플릿에서 기본 nf-core 호환 파이프라인을 생성합니다
    - 일반 Nextflow 워크플로우를 nf-core 표준에 맞게 업그레이드합니다
    - nf-core 호환 파이프라인에 nf-core 모듈을 추가합니다
    - 자신의 모듈을 nf-core에 기여합니다
    - nf-core 도구를 사용하여 입력과 매개변수를 검증합니다
  audience_prerequisites:
    - "**대상:** 이 과정은 기본 Nextflow에 이미 익숙하며 nf-core 리소스와 모범 사례 사용법을 배우고자 하는 학습자를 위해 설계되었습니다."
    - "**기술:** 명령줄, 기본 스크립팅 개념 및 일반적인 파일 형식에 대한 지식이 있다고 가정합니다."
    - "**과정:** [Hello Nextflow](../hello_nextflow/index.md) 과정 또는 이에 상응하는 과정을 완료해야 합니다."
    - "**영역:** 모든 연습은 영역에 구애받지 않으므로 사전 과학 지식이 필요하지 않습니다."
---

# Hello nf-core

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } AI 지원 번역 - [자세히 알아보기 및 개선 제안](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

**Hello nf-core는 nf-core 리소스와 모범 사례 사용에 대한 실습 입문 과정입니다.**

![nf-core logo](./img/nf-core-logo.png)

실용적인 예제와 안내된 연습을 통해 nf-core 호환 모듈과 파이프라인을 사용하고 개발하는 방법, 그리고 nf-core 도구를 효과적으로 활용하는 방법을 배우게 됩니다.

nf-core 모범 사례에 따라 파이프라인을 개발할 수 있는 기술과 자신감을 얻게 됩니다.

<!-- additional_information -->

## 과정 개요

이 과정은 실습 중심으로 설계되었으며, 목표 지향적 연습을 통해 정보를 단계적으로 학습합니다.

Nextflow를 사용하여 구축된 과학 파이프라인의 큐레이션된 세트를 개발하고 유지 관리하는 커뮤니티 활동인 [**nf-core**](https://nf-co.re/)와, 개방형 개발, 테스트 및 동료 검토를 촉진하는 관련 도구 및 가이드라인에 대해 소개합니다([Nat Biotechnol 38, 276–278 (2020)](https://www.nature.com/articles/s41587-020-0439-x), [Genome Biol 26, 228 (2025)](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-025-03673-9)).

nf-core 커뮤니티가 개발한 파이프라인은 모듈식이고, 확장 가능하며, 이식 가능하도록 설계되어 연구자들이 자신의 데이터와 컴퓨팅 리소스를 사용하여 쉽게 적응하고 실행할 수 있습니다.
프로젝트에서 시행하는 모범 사례 가이드라인은 파이프라인이 견고하고, 잘 문서화되며, 실제 데이터 세트에 대해 검증되도록 보장합니다.
이는 과학 분석의 신뢰성과 재현성을 높이는 데 도움이 되며, 궁극적으로 연구자들이 과학적 발견을 가속화할 수 있도록 합니다.

이 과정에서는 nf-core 파이프라인에 대한 모든 것을 다루지는 않습니다. nf-core는 커뮤니티가 수년에 걸쳐 개발한 많은 기능과 규칙을 포함하고 있기 때문입니다.
대신, 시작하는 데 도움이 되고 nf-core가 어떻게 작동하는지 이해하는 데 필요한 필수 개념에 초점을 맞출 것입니다.

### 수업 계획

이 과정은 nf-core 리소스 사용의 특정 측면에 각각 초점을 맞춘 다섯 부분으로 나뉩니다.

| 과정 챕터                                               | 요약                                                                                                                            | 예상 소요 시간 |
| ------------------------------------------------------- | ------------------------------------------------------------------------------------------------------------------------------- | -------------- |
| [파트 1: 데모 파이프라인 실행](./01_run_demo.md)        | 기존 nf-core 파이프라인을 실행하고 코드 구조를 검토하여 이러한 파이프라인이 기본 Nextflow 워크플로우와 무엇이 다른지 파악합니다 | 30분           |
| [파트 2: nf-core용 Hello 재작성](./02_rewrite_hello.md) | [Hello Nextflow](../hello_nextflow/index.md) 과정에서 생성한 간단한 워크플로우를 nf-core 템플릿 스캐폴드에 맞게 적응시킵니다    | 60분           |
| [파트 3: nf-core 모듈 사용](./03_use_module.md)         | 커뮤니티 모듈 라이브러리를 탐색하고 일반적인 생물정보학 도구를 적용하는 사전 구축되고 테스트된 모듈을 통합하는 방법을 배웁니다  | 30분           |
| [파트 4: nf-core 모듈 만들기](./04_make_module.md)      | nf-core가 정한 특정 구조, 명명 규칙 및 메타데이터 요구사항을 사용하여 자신만의 nf-core 스타일 모듈을 생성합니다                 | 30분           |
| [파트 5: 입력 검증 추가](./05_input_validation.md)      | nf-schema를 사용하여 명령줄 매개변수와 입력 데이터 파일 모두에 대한 입력 검증을 구현합니다                                      | 30분           |

이 과정을 마치면 nf-core 프로젝트가 제공하는 방대한 리소스를 활용할 수 있게 됩니다.

과정을 시작할 준비가 되셨습니까?

[학습 시작 :material-arrow-right:](00_orientation.md){ .md-button .md-button--primary }
