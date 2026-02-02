# 과정 요약

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } AI 지원 번역 - [자세히 알아보기 및 개선 제안](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

Nextflow Run 교육 과정을 완료하신 것을 축하합니다! 🎉

<!-- placeholder for video -->

## 여러분의 여정

매우 기본적인 workflow로 시작하여 실행하고, 출력을 찾고, 실행을 관리하는 방법을 배웠습니다.
그런 다음 해당 workflow의 점점 더 복잡한 버전을 통해 작업하고 channel과 연산자, 코드 모듈화, 컨테이너를 포함하여 Nextflow pipeline을 구동하는 필수 개념과 메커니즘을 인식하는 방법을 배웠습니다.
마지막으로 선호도와 컴퓨팅 인프라에 맞게 pipeline의 구성을 사용자 정의하는 방법을 배웠습니다.

### 배운 내용

이제 Hello pipeline의 실행을 관리하고, 구조화 방법을 설명하고, 관련된 주요 코드 조각을 식별할 수 있습니다.

- Hello workflow의 최종 형태는 텍스트 인사말이 포함된 CSV 파일을 입력으로 받습니다.
- 네 단계는 Nextflow process(`sayHello`, `convertToUpper`, `collectGreetings`, `cowpy`)로 구현되어 별도의 모듈 파일에 저장됩니다.
- 결과는 `results/`라는 디렉토리에 게시됩니다.
- pipeline의 최종 출력은 대문자 인사말을 말하는 캐릭터의 ASCII 아트가 포함된 일반 텍스트 파일입니다.

<figure class="excalidraw">
--8<-- "docs/en/docs/hello_nextflow/img/hello_pipeline_complete.svg"
</figure>

1. **`sayHello`:** 각 인사말을 자체 출력 파일에 씁니다(예: "Hello-output.txt")
2. **`convertToUpper`:** 각 인사말을 대문자로 변환합니다(예: "HELLO")
3. **`collectGreetings`:** 모든 대문자 인사말을 단일 배치 파일로 수집합니다
4. **`cowpy`:** `cowpy` 도구를 사용하여 ASCII 아트를 생성합니다

workflow 구성은 유연하고 재현 가능한 방식으로 입력 및 매개변수를 제공하는 것을 지원합니다.

### 습득한 기술

이 실습 과정을 통해 다음 방법을 배웠습니다:

- 로컬에서 Nextflow workflow 시작
- Nextflow가 생성한 출력(결과) 및 로그 파일 찾기 및 해석
- 간단한 다단계 workflow를 구성하는 핵심 Nextflow 구성 요소 인식
- 연산자 및 channel factory와 같은 다음 단계 개념 설명
- 다양한 컴퓨팅 환경에 맞게 pipeline 구성

이제 기존 Nextflow pipeline을 자신의 작업에 통합하기 위한 기초 지식을 갖추었습니다.

## 기술 향상을 위한 다음 단계

다음에 무엇을 할지에 대한 최고의 제안입니다:

- Nextflow를 실행만 하지 말고 작성하세요! [Hello Nextflow](../hello_nextflow/index.md)로 Nextflow 개발자가 되세요
- [Nextflow for Science](../nf4_science/index.md)로 과학적 분석 사용 사례에 Nextflow 적용
- [Hello nf-core](../hello_nf-core/index.md)로 nf-core 시작하기
- [디버깅 Side Quest](../side_quests/debugging.md)로 문제 해결 기술 배우기

마지막으로 [**Seqera Platform**](https://seqera.io/)을 살펴보시기를 권장합니다. Nextflow 제작자가 개발한 클라우드 기반 플랫폼으로 workflow를 시작하고 관리하고, 데이터를 관리하고, 모든 환경에서 대화형으로 분석을 실행하는 것을 더욱 쉽게 해줍니다.

## 도움 받기

도움 리소스 및 커뮤니티 지원은 [도움말 페이지](../help.md)를 참조하세요.

## 피드백 설문조사

다음으로 넘어가기 전에 과정 설문조사를 완료해 주세요! 여러분의 피드백은 모든 사람을 위한 교육 자료를 개선하는 데 도움이 됩니다.

[설문조사 참여 :material-arrow-right:](survey.md){ .md-button .md-button--primary }
