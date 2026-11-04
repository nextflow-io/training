# 과정 요약

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } AI 지원 번역 - [자세히 알아보기 및 개선 제안](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

nf-core 사용 교육 과정을 완료하신 것을 축하드립니다! 🎉

<!-- placeholder for video -->

## 학습 여정

`nf-core/demo` 파이프라인을 찾고 가져오는 것부터 시작하여, test 프로파일을 사용해 실행하고 출력 결과를 확인하는 방법을 학습했습니다.
이어서 파이프라인 매개변수와 설정 파일을 통해 실행을 설정하고, nf-core 파이프라인이 매개변수와 입력 데이터를 검증하는 방식을 살펴보았습니다.
마지막으로, 동일한 기술을 프로덕션 규모의 파이프라인인 `nf-core/rnaseq`에 적용하고, 사용 가능한 하드웨어에 맞게 기본 리소스 할당을 재정의하는 방법을 학습했습니다.

### 학습 내용

이제 nf-core 파이프라인을 찾고, 가져오고, 실행하고, 설정할 수 있습니다.

- nf-core 파이프라인은 `nextflow pull`로 가져오며, 표준 코드 구성을 따릅니다.
- 모든 nf-core 파이프라인에는 소규모 데이터셋에서 빠른 검증을 위한 `test` 프로파일이 포함되어 있습니다.
- 파이프라인 매개변수(`--param_name` 또는 `-params-file`로 설정)와 설정(`-c`로 설정)은 서로 다른 목적을 가집니다. 전자는 입력 및 분석 옵션을, 후자는 리소스 할당과 같은 실행 관련 설정을 담당합니다.
- nf-core 파이프라인은 매개변수와 입력 파일을 자동으로 검증하여, 작업이 시작되기 전에 오류를 감지합니다.
- 리소스 기본값은 `conf/base.config`에 정의된 레이블(`process_low`, `process_medium`, `process_high`)을 통해 할당되며, 사용자 정의 설정 파일로 재정의할 수 있습니다.

### 습득한 기술

이 실습 과정을 통해 다음을 수행하는 방법을 학습했습니다.

- nf-co.re 웹사이트에서 nf-core 파이프라인을 찾고 소스 코드 가져오기
- 내장된 test 프로파일을 사용하여 파이프라인을 실행하고 출력 결과 확인하기
- 도움말 확인, 매개변수 설정, 매개변수 및 입력 검증 이해하기
- 설정 파일을 통해 리소스 할당 및 도구 인자 맞춤화하기
- 프로덕션 규모의 파이프라인을 가져와 실행하고, 기본 리소스 레이블 재정의하기

이제 자신의 분석을 위해 nf-core 파이프라인을 실행하는 데 필요한 기초 지식을 갖추었습니다.

## 기술 향상을 위한 다음 단계

다음으로 진행할 수 있는 주요 제안 사항입니다.

- [Scale with Seqera](../seqera_scale/index.md)를 통해 대규모로 파이프라인을 실행하고 모니터링하기
- nf-core 파이프라인을 실행하는 것에 그치지 말고, 직접 개발해 보세요! [Build with nf-core](../nfcore_build/index.md)에서 nf-core 모범 사례를 학습하기
- Nextflow 자체가 처음이라면 [Nextflow Run](../nextflow_run/index.md)부터 시작하기
- [Nextflow for Science](../nf4_science/index.md)에서 과학적 분석 사례에 Nextflow 적용하기
- [Side Quests](../side_quests/index.md)에서 더 고급 Nextflow 기능 살펴보기

## 도움 받기

도움 자료 및 커뮤니티 지원은 [도움말 페이지](../help.md)를 참조하세요.

## 피드백 설문

다음 단계로 넘어가기 전에, 잠시 시간을 내어 과정 설문을 완료해 주세요! 여러분의 피드백은 모든 사람을 위한 교육 자료 개선에 도움이 됩니다.

[설문 참여하기 :material-arrow-right:](survey.md){ .md-button .md-button--primary }
