---
title: RNAseq
hide:
  - toc
page_type: index_page
index_type: course
additional_information:
  technical_requirements: true
  learning_objectives:
    - 기본 RNAseq 처리 및 QC 방법을 적용하는 선형 워크플로우 작성하기
    - FASTQ 및 참조 게놈 리소스와 같은 도메인별 파일을 적절하게 처리하기
    - 단일 말단 및 쌍 말단 시퀀싱 데이터 처리하기
    - Nextflow의 데이터 플로우 패러다임을 활용하여 샘플별 RNAseq 처리를 병렬화하기
    - 관련 채널 연산자를 사용하여 여러 단계와 샘플에 걸쳐 QC 리포트를 집계하기
  audience_prerequisites:
    - "**대상:** 이 과정은 데이터 분석 파이프라인을 개발하거나 맞춤화하는 데 관심이 있는 전사체학 및 관련 분야의 연구자를 위한 것입니다."
    - "**기술:** 명령줄, 기본 스크립팅 개념 및 일반적인 RNAseq 파일 형식에 대한 어느 정도의 친숙함을 가정합니다."
    - "**전제 조건:** [Hello Nextflow](../../hello_nextflow/)에서 다룬 기본적인 Nextflow 개념 및 도구."
---

# RNAseq을 위한 Nextflow

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } AI 지원 번역 - [자세히 알아보기 및 개선 제안](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

**실제 전사체학 사용 사례에 Nextflow를 적용하는 실습 과정: Trim Galore, HISAT2 및 FastQC를 사용한 bulk RNAseq 처리.**

이 과정은 [Hello Nextflow](../../hello_nextflow/) 초급 교육을 기반으로 하며 bulk RNAseq 분석의 특정 맥락에서 Nextflow를 사용하는 방법을 보여줍니다.
어댑터 서열을 트리밍하고, 리드를 게놈 참조에 정렬하며, 여러 단계에서 품질 관리(QC)를 수행하는 처리 파이프라인을 구현하게 됩니다.

<!-- additional_information -->

## 과정 개요

이 과정은 실습 중심이며, 정보를 단계적으로 학습할 수 있도록 구성된 목표 지향적 연습으로 이루어져 있습니다.

먼저 터미널에서 처리 도구를 수동으로 실행하여 방법론을 이해한 다음, 분석을 자동화하고 확장하는 Nextflow 파이프라인을 단계적으로 구축하게 됩니다.

### 학습 계획

이 과정은 RNAseq 사용 사례에 Nextflow를 적용하는 특정 측면에 초점을 맞춘 세 부분으로 나뉩니다.

| 과정 챕터                                              | 요약                                                                                        | 예상 소요 시간 |
| ------------------------------------------------------ | ------------------------------------------------------------------------------------------- | -------------- |
| [파트 1: 방법 개요](./01_method.md)                    | RNAseq 처리 방법론을 이해하고 도구를 수동으로 실행하기                                      | 30분           |
| [파트 2: 단일 샘플 구현](./02_single-sample.md)        | 단일 샘플을 트리밍, 정렬 및 QC하는 파이프라인을 구축한 다음 여러 샘플을 처리하도록 확장하기 | 60분           |
| [파트 3: 다중 샘플 쌍 말단 구현](./03_multi-sample.md) | 쌍 말단 데이터를 처리하고 샘플 전체에 걸쳐 QC 리포트를 집계하도록 파이프라인 확장하기       | 45분           |

이 과정을 마치면 일반적인 RNAseq 사용 사례에 기본적인 Nextflow 개념과 도구를 적용할 수 있게 됩니다.

과정을 시작할 준비가 되셨나요?

[시작하기 :material-arrow-right:](00_orientation.md){ .md-button .md-button--primary }

<!-- Clearfix for float -->
<div style="content: ''; clear: both; display: table;"></div>
