---
title: Hello Nextflow
hide:
    - toc
page_type: index_page
index_type: course
additional_information:
    technical_requirements: true
    learning_objectives:
        - Nextflow 워크플로우 실행 시작 및 관리
        - Nextflow가 생성한 출력(결과) 및 로그 파일 찾기 및 해석
        - 기본적인 문제 해결
        - 핵심 Nextflow 구성 요소로 간단한 다단계 워크플로우 구축
        - 필수 유형의 채널 팩토리와 연산자를 구분하고 간단한 워크플로우에서 효과적으로 활용
        - HPC 및 클라우드를 포함한 일반적인 컴퓨팅 플랫폼에서 실행되도록 파이프라인 실행 구성
        - 코드 모듈화 및 소프트웨어 컨테이너를 포함하여 파이프라인을 FAIR하게 만드는 재현성, 이식성 및 코드 재사용을 위한 모범 사례 적용
    audience_prerequisites:
        - "**대상:** 이 과정은 Nextflow를 처음 접하고 자체 파이프라인을 개발하고자 하는 학습자를 위해 설계되었습니다."
        - "**기술:** 명령줄, 기본 스크립팅 개념 및 일반적인 파일 형식에 대한 약간의 친숙함이 있다고 가정합니다."
        - "**도메인:** 모든 연습은 도메인에 구애받지 않으므로 사전 과학적 지식이 필요하지 않습니다."
    videos_playlist: https://youtube.com/playlist?list=PLPZ8WHdZGxmWKozQuzr27jyMGqp9kElVK&si=eF7cLR62goy-lc6n
---

# Hello Nextflow

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } AI 지원 번역 - [자세히 알아보기 및 개선 제안](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

**Hello Nextflow는 재현 가능하고 확장 가능한 데이터 분석 워크플로우 구축에 대한 실습 입문 과정입니다.**

실제 예제와 안내 연습을 통해 Nextflow로 파이프라인을 개발하는 기초를 학습합니다. 프로세스 정의, 파이프라인 연결, 파일 및 소프트웨어 의존성 관리, 손쉬운 실행 병렬화, 다양한 컴퓨팅 환경에서 워크플로우 실행 방법 등을 다룹니다.

이 과정을 통해 Nextflow로 자체 워크플로우를 개발하고 실행하기 시작할 수 있는 기술과 자신감을 갖추게 됩니다.

<!-- additional_information -->

## 과정 개요

이 과정은 실습 중심으로 설계되었으며, 목표 지향적 연습을 통해 정보를 단계적으로 학습합니다.

텍스트 입력을 받아 몇 가지 변환 단계를 실행하고, 변환된 텍스트를 말하는 캐릭터의 ASCII 그림이 포함된 단일 텍스트 파일을 출력하는 간단한 Nextflow 파이프라인을 개발합니다.

### 강의 계획

개념과 코드로 압도하지 않기 위해 Nextflow로 파이프라인을 개발하는 특정 측면에 초점을 맞춘 여섯 부분으로 나누었습니다.

| 과정 챕터                                            | 요약                                                                    | 예상 소요 시간 |
| ---------------------------------------------------- | ----------------------------------------------------------------------- | -------------- |
| [파트 1: Hello World](./01_hello_world.md)           | Nextflow 워크플로우를 조립하고 실행하는 데 관련된 기본 구성 요소 및 원리 | 30분           |
| [파트 2: Hello Channels](./02_hello_channels.md)     | 입력 처리 및 손쉬운 실행 병렬화를 위한 채널과 연산자 사용                | 45분           |
| [파트 3: Hello Workflow](./03_hello_workflow.md)     | 채널을 사용하여 여러 단계를 연결하고 단계 간 데이터 전송 처리            | 60분           |
| [파트 4: Hello Modules](./04_hello_modules.md)       | 재사용성을 높이고 유지 보수 부담을 줄이기 위한 코드 모듈화 원칙 적용     | 20분           |
| [파트 5: Hello Containers](./05_hello_containers.md) | 소프트웨어 의존성 관리 및 재현성 향상을 위한 메커니즘으로 컨테이너 사용  | 60분           |
| [파트 6: Hello Config](./06_hello_config.md)         | 파이프라인 동작 사용자 정의 및 다양한 컴퓨팅 환경에서 사용 최적화        | 60분           |

이 과정을 마치면 과학적 컴퓨팅 요구 사항에 맞는 재현 가능한 워크플로우를 개발하는 다음 단계를 수행할 준비가 잘 되어 있을 것입니다.

과정을 수강할 준비가 되셨나요?

[시작하기 :material-arrow-right:](00_orientation.md){ .md-button .md-button--primary }

<!-- Clearfix for float -->
<div style="content: ''; clear: both; display: table;"></div>
