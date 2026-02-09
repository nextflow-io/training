# Nextflow for Genomics

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } AI 기반 번역 - [자세히 알아보기 및 개선 사항 제안](../../hello_nextflow/)</span>

**실제 유전체학 사례에 Nextflow를 적용하는 실습 과정: GATK를 사용한 변이 호출.**

이 과정은 [Hello Nextflow](../../hello_nextflow/) 초급 교육을 기반으로 하며, 유전체학 분야의 특정 맥락에서 Nextflow를 사용하는 방법을 다룹니다.
고처리량 시퀀싱 데이터 분석에 널리 사용되는 소프트웨어 패키지인 [GATK](https://gatk.broadinstitute.org/) (Genome Analysis Toolkit)를 사용하여 변이 호출 파이프라인을 구현합니다.

<!-- additional_information -->

## 과정 개요

이 과정은 실습 중심이며, 정보를 단계적으로 학습할 수 있도록 목표 지향적인 연습으로 구성되어 있습니다.

먼저 터미널에서 변이 호출 도구를 수동으로 실행하여 방법론을 이해한 다음, 분석을 자동화하고 확장하는 Nextflow 파이프라인을 점진적으로 구축합니다.

### 학습 계획

유전체학 사례에 Nextflow를 적용하는 특정 측면에 초점을 맞춘 세 부분으로 나누어져 있습니다.

| 과정 챕터                                                                | 요약                                                                      | 예상 소요 시간 |
| ------------------------------------------------------------------------ | ------------------------------------------------------------------------- | -------------- |
| [Part 1: Method overview](./01_method.md)                                | 변이 호출 방법론 이해 및 도구 수동 실행                                   | 30분           |
| [Part 2: Per-sample variant calling](./02_per_sample_variant_calling.md) | BAM 파일 인덱싱 및 변이 호출 파이프라인 구축, 여러 샘플로 확장            | 60분           |
| [Part 3: Joint calling on a cohort](./03_joint_calling.md)               | 채널 연산자를 사용하여 샘플별 출력을 집계하는 다중 샘플 공동 유전자형 추가 | 45분           |

이 과정을 마치면 일반적인 유전체학 사례에 기본적인 Nextflow 개념과 도구를 적용할 수 있습니다.

과정을 시작할 준비가 되셨나요?

[시작하기 :material-arrow-right:](00_orientation.md){ .md-button .md-button--primary }

<!-- Clearfix for float -->
<div style="content: ''; clear: both; display: table;"></div>
