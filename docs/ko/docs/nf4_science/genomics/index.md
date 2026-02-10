---
title: Nextflow for Genomics
hide:
  - toc
page_type: index_page
index_type: course
additional_information:
  technical_requirements: true
  learning_objectives:
    - 단일 샘플에 변이 호출을 적용하는 선형 워크플로우 작성하기
    - 인덱스 파일 및 참조 게놈 리소스와 같은 보조 파일을 적절하게 처리하기
    - Nextflow의 데이터플로우 패러다임을 활용하여 샘플별 변이 호출을 병렬화하기
    - 관련 채널 연산자를 사용하여 다중 샘플 공동 호출 구현하기
  audience_prerequisites:
    - "**대상:** 이 과정은 데이터 분석 파이프라인을 개발하거나 맞춤화하려는 유전체학 및 관련 분야의 연구자를 위해 설계되었습니다."
    - "**기술:** 명령줄, 기본 스크립팅 개념, 일반적인 유전체학 파일 형식에 대한 어느 정도의 친숙함을 가정합니다."
    - "**전제 조건:** [Hello Nextflow](../../hello_nextflow/)에서 다루는 기초적인 Nextflow 개념 및 도구에 대한 이해가 필요합니다."
---

# Nextflow for Genomics

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } AI 지원 번역 - [자세히 알아보기 및 개선 제안](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

**실제 유전체학 사용 사례에 Nextflow를 적용하는 실습 과정: GATK를 사용한 변이 호출.**

이 과정은 [Hello Nextflow](../../hello_nextflow/) 초급 교육을 기반으로 하며, 유전체학 도메인의 특정 맥락에서 Nextflow를 사용하는 방법을 보여줍니다.
고처리량 시퀀싱 데이터를 분석하는 데 널리 사용되는 소프트웨어 패키지인 [GATK](https://gatk.broadinstitute.org/) (Genome Analysis Toolkit)를 사용하여 변이 호출 파이프라인을 구현합니다.

<!-- additional_information -->

## 과정 개요

이 과정은 실습 중심이며, 정보를 단계적으로 학습할 수 있도록 목표 지향적인 연습으로 구성되어 있습니다.

먼저 터미널에서 변이 호출 도구를 수동으로 실행하여 방법론을 이해한 다음, 분석을 자동화하고 확장하는 Nextflow 파이프라인을 점진적으로 구축합니다.

### 학습 계획

이 과정은 유전체학 사용 사례에 Nextflow를 적용하는 특정 측면에 초점을 맞춘 세 부분으로 나뉩니다.

| 과정 챕터                                                                | 요약                                                                               | 예상 소요 시간 |
| ------------------------------------------------------------------------ | ---------------------------------------------------------------------------------- | -------------- |
| [Part 1: Method overview](./01_method.md)                                | 변이 호출 방법론 이해 및 도구를 수동으로 실행하기                                  | 30분           |
| [Part 2: Per-sample variant calling](./02_per_sample_variant_calling.md) | BAM 파일을 인덱싱하고 변이를 호출하는 파이프라인 구축, 그리고 여러 샘플로 확장하기 | 60분           |
| [Part 3: Joint calling on a cohort](./03_joint_calling.md)               | 채널 연산자를 사용하여 샘플별 출력을 집계하는 다중 샘플 공동 유전자형 분석 추가    | 45분           |

이 과정을 마치면 기초적인 Nextflow 개념과 도구를 일반적인 유전체학 사용 사례에 적용할 수 있습니다.

과정을 시작할 준비가 되셨나요?

[시작하기 :material-arrow-right:](00_orientation.md){ .md-button .md-button--primary }

<!-- Clearfix for float -->
<div style="content: ''; clear: both; display: table;"></div>
