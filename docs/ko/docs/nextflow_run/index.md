---
title: Nextflow Run
hide:
  - toc
page_type: index_page
index_type: course
additional_information:
  technical_requirements: true
  learning_objectives:
    - 명령줄에서 Nextflow 파이프라인 실행 및 관리
    - 채널과 연산자가 다중 입력, 다단계 워크플로우의 효율적인 처리를 가능하게 하는 방식 이해
    - 컨테이너를 사용하여 소프트웨어 의존성 관리 및 재현성 보장
    - 파이프라인 실행 및 출력 설정
    - 실행 보고서 생성, 이전 실행 기록 검사, 오래된 work 디렉토리 정리
    - GitHub 등 원격 저장소에서 직접 파이프라인 실행
  audience_prerequisites:
    - "**대상:** 이 과정은 Nextflow를 처음 접하고 기존 파이프라인을 실행하려는 학습자를 위해 설계되었습니다."
    - "**기술:** 명령줄, 기본 스크립팅 개념 및 일반적인 파일 형식에 대한 어느 정도의 친숙함이 필요합니다."
    - "**도메인:** 모든 연습은 도메인에 구애받지 않으므로 사전 과학 지식이 필요하지 않습니다."
---

# Nextflow Run

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } AI 지원 번역 - [자세히 알아보기 및 개선 제안](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

**Nextflow Run은 재현 가능하고 확장 가능한 데이터 분석 워크플로우 실행에 대한 실습 입문 과정입니다.**

목표 지향적인 연습을 통해 Nextflow 파이프라인 실행 및 관리의 기본 사항을 학습하고, 채널과 연산자가 다중 입력의 병렬 처리를 가능하게 하는 방식을 이해하며, 컨테이너를 사용하여 소프트웨어 의존성을 관리합니다.

Nextflow로 워크플로우를 실행하기 위한 기술과 자신감을 갖추게 됩니다.

<!-- additional_information -->

## 과정 개요

이 과정은 실습 중심이며, 목표 지향적인 연습을 통해 정보를 단계적으로 학습합니다.

텍스트 입력을 처리하는 Nextflow 파이프라인의 여러 버전을 실행합니다.
단일 단계로 구성된 간단한 버전부터 시작하여, CSV 파일의 입력을 받아 몇 가지 변환 단계를 거친 후, 컨테이너화된 도구가 생성한 ASCII 아트를 포함하는 단일 텍스트 파일을 출력하는 다단계 버전으로 진행합니다.

이 과정은 파이프라인 실행에 중점을 둡니다(핵심 `nextflow run` 명령의 이름을 따서 명명됨).
Nextflow 파이프라인 개발에 대한 입문을 찾고 있다면 [Hello Nextflow](../hello_nextflow/index.md)를 참조하세요.

!!! note "참고"

    이 과정의 이전 버전을 찾고 계신가요? 현재 페이지의 버전으로 대체되었지만, 교육 사이트의 [3.6.1 릴리스](https://training.nextflow.io/3.6.1/nextflow_run/)에서 계속 확인할 수 있습니다.

### 학습 계획

| 과정 챕터                                                   | 요약                                                             | 예상 소요 시간 |
| ----------------------------------------------------------- | ---------------------------------------------------------------- | -------------- |
| [파트 1: Nextflow 실행](./01_run_nextflow.md)               | Nextflow 파이프라인 실행 및 관리, 필수 워크플로우 메커니즘 이해  | 25분           |
| [파트 2: 파이프라인 설정](./02_configure_pipeline.md)       | `nextflow.config`를 사용하여 파이프라인 실행 및 출력 설정        | 20분           |
| [파트 3: 워크플로우 실행 관리](./03_manage_executions.md)   | 실행 보고서 생성, 이전 실행 기록 검사, 오래된 work 디렉토리 정리 | 10분           |
| [파트 4: 원격 파이프라인 실행](./04_remote_repositories.md) | GitHub에서 직접 파이프라인을 실행하고 특정 리비전으로 고정       | 10분           |

이 과정이 끝나면 과학 컴퓨팅 요구 사항에 맞는 재현 가능한 워크플로우를 실행하기 위한 다음 단계를 수행할 준비가 됩니다.

과정을 시작할 준비가 되셨나요?

[학습 시작 :material-arrow-right:](00_orientation.md){ .md-button .md-button--primary }
