# 과정 요약

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } AI 지원 번역 - [자세히 알아보기 및 개선 제안](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

Hello Nextflow 교육 과정을 완료하신 것을 축하드립니다! 🎉

<!--
<div class="video-wrapper">
  <iframe width="560" height="315" src="https://www.youtube.com/embed/xHOcx_4Ancg?si=Lp8hS8RdaMwbp5j5&amp;list=PLPZ8WHdZGxmXiHf8B26oB_fTfoKQdhlik" title="YouTube video player" frameborder="0" allow="accelerometer; autoplay; clipboard-write; encrypted-media; gyroscope; picture-in-picture; web-share" referrerpolicy="strict-origin-when-cross-origin" allowfullscreen></iframe>
</div>

/// caption
:fontawesome-brands-youtube:{ .youtube } Nextflow YouTube 채널에서 [전체 재생 목록](https://www.youtube.com/playlist?list=PLPZ8WHdZGxmXiHf8B26oB_fTfoKQdhlik)을 확인하세요.

:green_book: 비디오와 함께 [비디오 스크립트](./transcripts/07_next_steps.md)를 읽을 수 있습니다.
///
-->

## 여러분의 여정

하드코딩된 명령을 실행하는 매우 기본적인 워크플로우로 시작했습니다.
여섯 부분에 걸쳐 그 기본 워크플로우를 채널, 연산자, 컨테이너에 대한 내장 지원 및 구성 옵션을 포함한 Nextflow의 핵심 기능을 활용하는 모듈식 다단계 파이프라인으로 변환했습니다.

### 구축한 것

- Hello 워크플로우의 최종 형태는 텍스트 인사말이 포함된 CSV 파일을 입력으로 받습니다.
- 네 단계는 별도의 모듈 파일에 저장된 Nextflow 프로세스(`sayHello`, `convertToUpper`, `collectGreetings`, `cowpy`)로 구현됩니다.
- 결과는 `results/`라는 디렉토리에 게시됩니다.
- 파이프라인의 최종 출력은 대문자로 변환된 인사말을 말하는 캐릭터의 ASCII 아트가 포함된 일반 텍스트 파일입니다.

<figure class="excalidraw">
--8<-- "docs/en/docs/hello_nextflow/img/hello_pipeline_complete.svg"
</figure>

1. **`sayHello`:** 각 인사말을 자체 출력 파일에 씁니다 (예: "Hello-output.txt")
2. **`convertToUpper`:** 각 인사말을 대문자로 변환합니다 (예: "HELLO")
3. **`collectGreetings`:** 모든 대문자 인사말을 단일 배치 파일로 수집합니다
4. **`cowpy`:** `cowpy` 도구를 사용하여 ASCII 아트를 생성합니다

워크플로우 구성은 유연하고 재현 가능한 방식으로 입력과 매개변수를 제공할 수 있도록 지원합니다.

### 습득한 기술

이 실습 과정을 통해 다음을 배웠습니다:

- 간단한 다단계 워크플로우를 구축하기에 충분한 핵심 Nextflow 구성 요소 설명 및 활용
- 연산자 및 채널 팩토리와 같은 다음 단계 개념 설명
- Nextflow 워크플로우를 로컬에서 실행
- Nextflow가 생성한 출력(결과) 및 로그 파일 찾기 및 해석
- 기본적인 문제 해결

이제 Nextflow에서 자체 파이프라인을 개발하기 시작할 수 있는 기초 지식을 갖추게 되었습니다.

## 기술 향상을 위한 다음 단계

다음에 할 일에 대한 상위 3가지 제안입니다:

- [과학을 위한 Nextflow](../nf4_science/index.md)로 과학적 분석 사용 사례에 Nextflow 적용
- [Hello nf-core](../../hello_nf-core/index.md)로 nf-core 시작하기
- [Side Quests](../side_quests/index.md)로 더 고급 Nextflow 기능 탐색

마지막으로, Nextflow 제작자가 개발한 클라우드 기반 플랫폼인 [**Seqera Platform**](https://seqera.io/)을 살펴보시기 바랍니다. 이 플랫폼은 워크플로우 시작 및 관리, 데이터 관리, 모든 환경에서 대화형 분석 실행을 더욱 쉽게 해줍니다.

## 피드백 설문조사

계속 진행하기 전에 과정 설문조사를 완료해 주세요! 여러분의 피드백은 모든 사람을 위한 교육 자료를 개선하는 데 도움이 됩니다.

[설문조사 참여하기 :material-arrow-right:](survey.md){ .md-button .md-button--primary }
