---
title: 유전체학을 위한 Nextflow
hide:
  - toc
---

# 유전체학을 위한 Nextflow

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } AI 지원 번역 - [자세히 알아보기 및 개선 제안](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

이 교육 과정은 데이터 분석 파이프라인을 개발하거나 맞춤화하는 데 관심이 있는 유전체학 및 관련 분야의 연구자를 위한 것입니다.
이 과정은 [Hello Nextflow](../../hello_nextflow/) 초급 교육을 기반으로 하며, 유전체학 도메인의 특정 맥락에서 Nextflow를 사용하는 방법을 보여줍니다.

구체적으로, 이 과정은 고처리량 시퀀싱 데이터를 분석하기 위해 널리 사용되는 소프트웨어 패키지인 [GATK](https://gatk.broadinstitute.org/) (Genome Analysis Toolkit)를 사용하여 간단한 변이 호출 파이프라인을 구현하는 방법을 보여줍니다.

시작해봅시다! 아래의 "Open in GitHub Codespaces" 버튼을 클릭하여 교육 환경을 실행하고(가급적 별도의 탭에서), 로딩되는 동안 계속 읽어주십시오.

[![Open in GitHub Codespaces](https://github.com/codespaces/badge.svg)](https://codespaces.new/nextflow-io/training?quickstart=1&ref=master)

## 학습 목표

이 과정을 통해 기초적인 Nextflow 개념과 도구를 일반적인 유전체학 사용 사례에 적용하는 방법을 배우게 됩니다.

이 워크숍을 마치면 다음을 수행할 수 있게 됩니다:

- 단일 샘플에 변이 호출을 적용하는 선형 workflow 작성하기
- 인덱스 파일 및 참조 유전체 리소스와 같은 보조 파일을 적절하게 처리하기
- Nextflow의 dataflow 패러다임을 활용하여 샘플별 변이 호출을 병렬화하기
- 관련 channel 연산자를 사용하여 다중 샘플 변이 호출 구현하기
- 유전체학 특유의 특성을 적절하게 처리하는 단계별 및 종단간 파이프라인 테스트 구현하기

<!-- TODO for future expansion: add metadata/samplesheet handling -->

## 전제 조건

이 과정은 다음에 대한 최소한의 친숙함을 가정합니다:

- 이 과학 분야에서 일반적으로 사용되는 도구 및 파일 형식
- 명령줄 경험
- [Hello Nextflow](../../hello_nextflow/) 초급 교육에서 다루는 기초적인 Nextflow 개념 및 도구

기술 요구 사항 및 환경 설정은 [환경 설정](../../envsetup/) 단기 과정을 참조하십시오.
