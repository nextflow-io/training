---
title: Nextflow Run
hide:
  - toc
page_type: index_page
index_type: course
additional_information:
  technical_requirements: true
  learning_objectives:
    - Nextflow workflow 실행 및 관리
    - 출력(결과) 및 로그 파일 찾기 및 해석
    - 간단한 다단계 workflow에서 핵심 Nextflow 구성 요소 인식
    - HPC 및 클라우드를 포함한 일반적인 컴퓨팅 플랫폼에서 실행하도록 pipeline 구성
    - 코드 모듈화 및 소프트웨어 컨테이너를 포함하여 pipeline을 FAIR하게 만드는 재현성, 이식성 및 코드 재사용을 위한 모범 사례 요약
  audience_prerequisites:
    - "**대상:** 이 과정은 Nextflow를 처음 접하고 기존 pipeline을 실행하려는 학습자를 위해 설계되었습니다."
    - "**기술:** 명령줄, 기본 스크립팅 개념 및 일반적인 파일 형식에 대한 어느 정도의 친숙함이 필요합니다."
    - "**도메인:** 모든 연습은 도메인에 구애받지 않으므로 사전 과학 지식이 필요하지 않습니다."
---

# Nextflow Run

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } AI 지원 번역 - [자세히 알아보기 및 개선 제안](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

**Nextflow Run은 재현 가능하고 확장 가능한 데이터 분석 workflow 실행에 대한 실습 입문 과정입니다.**

실용적인 예제와 안내된 연습을 통해 pipeline 실행, 파일 및 소프트웨어 의존성 관리, 손쉬운 병렬 실행, 다양한 컴퓨팅 환경에서의 workflow 실행 방법을 포함한 Nextflow 사용의 기본 사항을 배웁니다.

Nextflow로 workflow를 실행하기 위한 기술과 자신감을 얻게 됩니다.

<!-- additional_information -->

## 과정 개요

### 학습 내용

이 과정은 실습 중심이며, 목표 지향적인 연습을 통해 정보를 단계적으로 학습합니다.

텍스트 입력을 처리하는 Nextflow pipeline의 여러 버전을 실행합니다.
단일 단계로 구성된 간단한 버전부터 시작하여, CSV 파일의 표 형식 텍스트 입력을 받아 몇 가지 변환 단계를 거친 후, 변환된 텍스트를 말하는 캐릭터의 ASCII 그림을 포함하는 단일 텍스트 파일을 출력하는 다단계 버전으로 진행합니다.

이 과정은 pipeline 실행에 중점을 둡니다(핵심 `nextflow run` 명령의 이름을 따서 명명됨).
Nextflow pipeline 개발에 대한 입문을 찾고 있다면 [Hello Nextflow](../hello_nextflow/index.md)를 참조하세요.

### 학습 계획

Nextflow로 작성된 pipeline을 실행하고 관리하는 특정 측면에 초점을 맞춘 세 부분으로 나누었습니다.

| 과정 챕터                                      | 요약                                                                      | 예상 소요 시간 |
| ---------------------------------------------- | ------------------------------------------------------------------------- | -------------- |
| [파트 1: 기본 작업 실행](./01_basics.md)       | 간단한 workflow의 실행 및 관리                                            | 30분           |
| [파트 2: 실제 pipeline 실행](./02_pipeline.md) | 복잡한 입력 처리, 다단계 workflow 실행, 컨테이너 사용 및 손쉬운 병렬 실행 | 60분           |
| [파트 3: 실행 구성](./03_config.md)            | pipeline 동작 사용자 정의 및 다양한 컴퓨팅 환경에서의 사용 최적화         | 60분           |

이 과정이 끝나면 과학 컴퓨팅 요구 사항에 맞는 재현 가능한 workflow를 실행하기 위한 다음 단계를 수행할 준비가 됩니다.

과정을 시작할 준비가 되셨나요?

[학습 시작 :material-arrow-right:](00_orientation.md){ .md-button .md-button--primary }
