---
title: 유전체학을 위한 Nextflow
hide:
  - toc
page_type: index_page
index_type: course
additional_information:
  technical_requirements: true
  learning_objectives:
    - 단일 샘플에 변이 호출을 적용하는 선형 workflow 작성하기
    - 인덱스 파일 및 참조 유전체 리소스와 같은 보조 파일을 적절하게 처리하기
    - Nextflow의 dataflow 패러다임을 활용하여 샘플별 변이 호출을 병렬화하기
    - 관련 channel 연산자를 사용하여 다중 샘플 변이 호출 구현하기
  audience_prerequisites:
    - "**대상:** 이 과정은 데이터 분석 파이프라인을 개발하거나 맞춤화하려는 유전체학 및 관련 분야의 연구자를 위해 설계되었습니다."
    - "**기술:** 명령줄, 기본 스크립팅 개념, 일반적인 유전체학 파일 형식에 대한 어느 정도의 친숙함을 가정합니다."
    - "**전제 조건:** [Hello Nextflow](../../hello_nextflow/)에서 다루는 기초적인 Nextflow 개념 및 도구"
---

# 유전체학을 위한 Nextflow

**실제 유전체학 사용 사례인 GATK를 사용한 변이 호출에 Nextflow를 적용하는 실습 과정입니다.**

이 과정은 [Hello Nextflow](../../hello_nextflow/) 초급 교육을 기반으로 하며, 유전체학 도메인의 특정 맥락에서 Nextflow를 사용하는 방법을 보여줍니다.
고처리량 시퀀싱 데이터를 분석하기 위해 널리 사용되는 소프트웨어 패키지인 [GATK](https://gatk.broadinstitute.org/) (Genome Analysis Toolkit)를 사용하여 변이 호출 파이프라인을 구현하게 됩니다.

<!-- additional_information -->

## 과정 개요

이 과정은 정보를 단계적으로 소개하도록 구성된 목표 지향적 연습을 통한 실습 과정입니다.

먼저 터미널에서 변이 호출 도구를 수동으로 실행하여 방법론을 이해한 다음, 분석을 자동화하고 확장하는 Nextflow 파이프라인을 점진적으로 구축하게 됩니다.

### 수업 계획

유전체학 사용 사례에 Nextflow를 적용하는 특정 측면에 각각 초점을 맞춘 세 부분으로 나누었습니다.

| 과정 챕터                                                      | 요약                                                                          | 예상 소요 시간 |
| -------------------------------------------------------------- | ----------------------------------------------------------------------------- | -------------- |
| [Part 1: 방법 개요](./01_method.md)                            | 변이 호출 방법론 이해 및 도구 수동 실행                                       | 30분           |
| [Part 2: 샘플별 변이 호출](./02_per_sample_variant_calling.md) | BAM 파일을 인덱싱하고 변이를 호출하는 파이프라인 구축, 여러 샘플로 확장       | 60분           |
| [Part 3: 코호트에 대한 결합 호출](./03_joint_calling.md)       | channel 연산자를 사용하여 샘플별 출력을 집계하는 다중 샘플 결합 유전자형 추가 | 45분           |

이 과정을 마치면 기초적인 Nextflow 개념과 도구를 일반적인 유전체학 사용 사례에 적용할 수 있게 됩니다.

과정을 시작할 준비가 되셨나요?

[시작하기 :material-arrow-right:](00_orientation.md){ .md-button .md-button--primary }

<!-- Clearfix for float -->
<div style="content: ''; clear: both; display: table;"></div>
