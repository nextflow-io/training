---
title: 실행 설정
hide:
  - toc
page_type: index_page
index_type: course
additional_information:
  technical_requirements: true
  learning_objectives:
    - Docker와 Conda 간 소프트웨어 패키징 기술 전환
    - 실행 플랫폼 선택 및 Nextflow가 작업 실행을 해당 플랫폼에 맞게 조정하는 방식 이해
    - 컴퓨팅 리소스 할당 제어 및 실패한 작업 자동 재시도
    - 사전 설정된 설정 간 전환을 위한 프로파일 정의 및 조합
  audience_prerequisites:
    - "**대상:** 이 과정은 로컬 Nextflow 파이프라인을 실행하는 방법을 이미 알고 있으며, 실행 설정을 더 깊이 이해하고자 하는 학습자를 위해 설계되었습니다."
    - "**기술:** 명령줄에 대한 기본적인 이해가 필요합니다."
    - "**선수 과정:** [Nextflow Run](../nextflow_run/index.md)을 완료했거나, `nextflow run`으로 로컬 파이프라인을 실행하는 데 익숙해야 합니다."
---

# 실행 설정

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } AI 지원 번역 - [자세히 알아보기 및 개선 제안](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>


**Execution Config는 Nextflow 파이프라인 실행을 다양한 컴퓨팅 환경에 맞게 조정하는 방법을 실습 중심으로 학습하는 과정입니다.**

목표 지향적인 실습을 통해 소프트웨어 패키징 기술 전환, 실행 플랫폼 선택, 컴퓨팅 리소스 할당 및 재시도 제어, 그리고 설정을 전환 가능한 프로파일로 묶는 방법을 학습합니다.

이 과정을 마치면 Nextflow 파이프라인 실행을 전문적으로 설정할 수 있는 역량과 자신감을 갖추게 됩니다.

<!-- additional_information -->

## 과정 개요

이 과정은 실습 중심으로 구성되어 있으며, [Nextflow Run](../nextflow_run/index.md)에서 다룬 기술을 기반으로 합니다.

해당 과정에서 사용한 동일한 다단계 파이프라인을 활용하여 다양한 컴퓨팅 환경에 맞게 설정을 단계적으로 조정하고, 런타임에 전환할 수 있는 프로파일로 모든 설정을 묶는 방법을 학습합니다.

### 학습 계획

| 과정 챕터                                                              | 요약                                                    | 예상 소요 시간 |
| ---------------------------------------------------------------------- | ------------------------------------------------------- | -------------- |
| [파트 1: 컴퓨팅 환경에 맞게 조정하기](./01_packaging_and_execution.md) | 소프트웨어 패키징 기술 전환 및 실행 플랫폼 선택         | 20분           |
| [파트 2: 컴퓨팅 리소스 및 실패 관리](./02_resources_and_retries.md)    | 리소스 할당 제어 및 실패한 작업 자동 재시도             | 15분           |
| [파트 3: 프로파일을 사용한 설정 전환](./03_profiles.md)                | 프로파일 정의 및 조합, 최종 완성된 설정 확인            | 15분           |

이 과정을 마치면 다양한 컴퓨팅 환경에 맞게 Nextflow 파이프라인을 설정하고, 최소한의 번거로움으로 환경 간 전환을 자유롭게 할 수 있습니다.

과정을 시작할 준비가 되셨나요?

[학습 시작하기 :material-arrow-right:](00_orientation.md){ .md-button .md-button--primary }
