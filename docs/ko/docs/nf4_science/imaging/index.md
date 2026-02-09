# Nextflow run for Imaging

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } AI 지원 번역 - [자세히 알아보기 및 개선 사항 제안](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

이 교육 과정은 데이터 분석 파이프라인을 실행하고 사용자 정의하는 데 관심이 있는 이미징 및 공간 생물학 연구자를 위한 과정입니다.
Molecular Cartography 공간 전사체학 데이터를 처리하는 파이프라인인 [nf-core/molkart](https://nf-co.re/molkart)를 사용하여 워크플로우를 실행하고, 구성하고, 설정하는 것과 관련된 기본적인 Nextflow 개념을 학습합니다.
여기서 배우는 기술은 모든 Nextflow 또는 nf-core 파이프라인에 적용할 수 있습니다.

시작해 봅시다! 아래의 "Open in GitHub Codespaces" 버튼을 클릭하여 교육 환경을 실행하세요(별도의 탭에서 여는 것을 권장합니다). 그런 다음 로딩되는 동안 계속 읽어주세요.

[![Open in GitHub Codespaces](https://github.com/codespaces/badge.svg)](https://codespaces.new/nextflow-io/training?quickstart=1&ref=master)

## 학습 목표

이 과정을 통해 이미징 분석 파이프라인을 실행하는 데 필요한 기본적인 Nextflow 개념과 도구를 적용하는 방법을 학습합니다.

이 워크숍을 마치면 다음을 수행할 수 있습니다:

- Nextflow 워크플로우를 로컬에서 실행하고 실행 과정을 모니터링하기
- Nextflow가 생성한 출력(결과) 및 로그 파일을 찾고 해석하기
- 테스트 데이터 및 사용자 정의 입력으로 nf-core 파이프라인 실행하기
- 프로파일 및 매개변수 파일을 사용하여 파이프라인 실행 구성하기
- 샘플시트 및 명령줄 매개변수를 사용하여 입력 관리하기

## 대상 및 사전 요구 사항

이 과정은 다음에 대한 최소한의 지식을 가정합니다:

- 명령줄 사용 경험
- 이미징 파일 형식에 대한 기본적인 이해(TIFF 이미지, 표 형식 데이터)

기술 요구 사항 및 환경 설정에 대해서는 [환경 설정](../../envsetup/) 단기 과정을 참조하세요.
