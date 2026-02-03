---
title: 이미징을 위한 Nextflow 실행
hide:
  - toc
---

# 이미징을 위한 Nextflow 실행

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } AI 지원 번역 - [자세히 알아보기 및 개선 제안](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

이 교육 과정은 데이터 분석 파이프라인을 실행하고 사용자 정의하는 데 관심이 있는 이미징 및 공간 생물학 연구자를 위한 것입니다.
이 과정은 Molecular Cartography 공간 전사체학 데이터를 처리하는 파이프라인인 [nf-core/molkart](https://nf-co.re/molkart)를 사용하여 워크플로우를 실행하고, 구성하고, 설정하는 것과 관련된 기본 Nextflow 개념을 가르칩니다.
여기서 배우는 기술은 모든 Nextflow 또는 nf-core 파이프라인에 적용할 수 있습니다.

시작해봅시다! 아래의 "Open in GitHub Codespaces" 버튼을 클릭하여 교육 환경을 실행하고(가급적 별도의 탭에서), 로드되는 동안 계속 읽어주십시오.

[![Open in GitHub Codespaces](https://github.com/codespaces/badge.svg)](https://codespaces.new/nextflow-io/training?quickstart=1&ref=master)

## 학습 목표

이 과정을 통해 이미징 분석 파이프라인을 실행하는 데 기본 Nextflow 개념과 도구를 적용하는 방법을 배우게 됩니다.

이 워크샵을 마치면 다음을 수행할 수 있습니다:

- 로컬에서 Nextflow 워크플로우를 실행하고 실행을 모니터링합니다
- Nextflow에서 생성된 출력(결과) 및 로그 파일을 찾고 해석합니다
- 테스트 데이터 및 사용자 정의 입력으로 nf-core 파이프라인을 실행합니다
- 프로파일 및 매개변수 파일을 사용하여 파이프라인 실행을 구성합니다
- 샘플시트 및 명령줄 매개변수를 사용하여 입력을 관리합니다

## 대상 및 전제 조건

이 과정은 다음에 대한 최소한의 지식을 전제로 합니다:

- 명령줄 사용 경험
- 이미징 파일 형식에 대한 기본적인 지식(TIFF 이미지, 표 형식 데이터)

기술적 요구 사항 및 환경 설정에 대해서는 [환경 설정](../../envsetup/) 단기 과정을 참조하십시오.
