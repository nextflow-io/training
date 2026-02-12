---
title: Nextflow for {DOMAIN}
hide:
  - toc
page_type: index_page
index_type: course
additional_information:
  technical_requirements: true
  learning_objectives:
    - Write a linear workflow to apply {METHOD} to a single sample
    - Handle accessory files such as {ACCESSORY_FILES} appropriately
    - Leverage Nextflow's dataflow paradigm to parallelize per-sample processing
    - Implement multi-sample aggregation using relevant channel operators
  audience_prerequisites:
    - "**Audience:** This course is designed for researchers in {DOMAIN} and related fields who want to develop or customize data analysis pipelines."
    - "**Skills:** Some familiarity with the command line, basic scripting concepts, and common {DOMAIN} file formats is assumed."
    - "**Prerequisites:** Foundational Nextflow concepts and tooling covered in [Hello Nextflow](../../hello_nextflow/)."
---

# Nextflow for {DOMAIN}

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } AI 지원 번역 - [자세히 알아보기 및 개선 제안](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

**{METHOD_SHORT_DESCRIPTION}를 다루는 실제 {DOMAIN} 사용 사례에 Nextflow를 적용하는 실습 과정입니다.**

이 과정은 [Hello Nextflow](../../hello_nextflow/) 초급 교육을 기반으로 하며, {DOMAIN} 분야의 특정 맥락에서 Nextflow를 사용하는 방법을 보여줍니다.
[{TOOL_A}]({TOOL_A_URL})와 [{TOOL_B}]({TOOL_B_URL})를 사용하여 {METHOD} 파이프라인을 구현합니다.

<!-- additional_information -->

## 과정 개요

이 과정은 실습 중심이며, 정보를 단계적으로 소개하도록 구성된 목표 지향적 연습으로 이루어져 있습니다.

먼저 터미널에서 분석 도구를 수동으로 실행하여 방법론을 이해한 다음, 분석을 자동화하고 확장하는 Nextflow 파이프라인을 점진적으로 구축합니다.

### 학습 계획

이 과정은 Nextflow를 {DOMAIN} 사용 사례에 적용하는 특정 측면에 초점을 맞춘 세 부분으로 나뉩니다.

| 과정 챕터                                                 | 요약                                                                                | 예상 소요 시간 |
| --------------------------------------------------------- | ----------------------------------------------------------------------------------- | -------------- |
| [Part 1: Method overview](./01_method.md)                 | {METHOD} 방법론을 이해하고 도구를 수동으로 실행하기                                 | 30분           |
| [Part 2: Single-sample processing](./02_single_sample.md) | {PART2_SUMMARY}하는 파이프라인을 구축한 다음 여러 샘플로 확장하기                   | 60분           |
| [Part 3: Multi-sample aggregation](./03_multi_sample.md)  | 채널 연산자를 사용하여 샘플별 출력을 집계하는 다중 샘플 {AGGREGATION_SUMMARY} 추가 | 45분           |

이 과정을 마치면 일반적인 {DOMAIN} 사용 사례에 기본적인 Nextflow 개념과 도구를 적용할 수 있습니다.

과정을 시작할 준비가 되셨나요?

[시작하기 :material-arrow-right:](00_orientation.md){ .md-button .md-button--primary }

<!-- Clearfix for float -->
<div style="content: ''; clear: both; display: table;"></div>
