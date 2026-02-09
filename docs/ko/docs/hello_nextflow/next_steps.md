# 과정 요약

Hello Nextflow 교육 과정을 완료하신 것을 축하합니다! 🎉

<div class="video-wrapper">
  <iframe width="560" height="315" src="https://www.youtube.com/embed/wnmUAfDL01E?si=y8lAedhEHWaTV4zd&amp;list=PLPZ8WHdZGxmWKozQuzr27jyMGqp9kElVK&amp;cc_load_policy=1" title="YouTube video player" frameborder="0" allow="accelerometer; autoplay; clipboard-write; encrypted-media; gyroscope; picture-in-picture; web-share" referrerpolicy="strict-origin-when-cross-origin" allowfullscreen></iframe>
</div>

/// caption
:fontawesome-brands-youtube:{ .youtube } [Nextflow YouTube 채널의 전체 재생목록](https://youtube.com/playlist?list=PLPZ8WHdZGxmWKozQuzr27jyMGqp9kElVK&si=eF7cLR62goy-lc6n)을 확인하세요.

:green_book: 영상과 함께 [영상 스크립트](./transcripts/07_next_steps.md)를 읽을 수 있습니다.
///

## 여러분의 학습 여정

하드코딩된 명령을 실행하는 매우 기본적인 워크플로우로 시작했습니다.
6개 파트를 거치면서, 그 기본 워크플로우를 채널, 연산자, 컨테이너에 대한 내장 지원, 구성 옵션을 포함한 Nextflow의 주요 기능을 활용하는 모듈식 다단계 파이프라인으로 변환했습니다.

### 구축한 내용

- Hello 워크플로우의 최종 형태는 텍스트 인사말이 포함된 CSV 파일을 입력으로 받습니다.
- 4개의 단계는 별도의 모듈 파일에 저장된 Nextflow 프로세스(`sayHello`, `convertToUpper`, `collectGreetings`, `cowpy`)로 구현됩니다.
- 결과는 `results/`라는 디렉토리에 게시됩니다.
- 파이프라인의 최종 출력은 대문자로 변환된 인사말을 말하는 캐릭터의 ASCII 아트가 포함된 일반 텍스트 파일입니다.

<figure class="excalidraw">
--8<-- "docs/en/docs/hello_nextflow/img/hello_pipeline_complete.svg"
</figure>

1. **`sayHello`:** 각 인사말을 개별 출력 파일에 작성합니다 (_예:_ "Hello-output.txt")
2. **`convertToUpper`:** 각 인사말을 대문자로 변환합니다 (_예:_ "HELLO")
3. **`collectGreetings`:** 모든 대문자 인사말을 단일 배치 파일로 수집합니다
4. **`cowpy`:** `cowpy` 도구를 사용하여 ASCII 아트를 생성합니다

워크플로우 구성은 유연하고 재현 가능한 방식으로 입력과 매개변수를 제공하는 것을 지원합니다.

### 습득한 기술

이 실습 과정을 통해 다음을 학습했습니다:

- 간단한 다단계 워크플로우를 구축하기에 충분한 핵심 Nextflow 구성 요소를 설명하고 활용하기
- 연산자 및 채널 팩토리와 같은 다음 단계 개념 설명하기
- Nextflow 워크플로우를 로컬에서 실행하기
- Nextflow가 생성한 출력(결과) 및 로그 파일을 찾고 해석하기
- 기본적인 문제 해결하기

이제 Nextflow에서 자신만의 파이프라인을 개발하기 시작할 수 있는 기초 지식을 갖추게 되었습니다.

## 기술을 향상시키기 위한 다음 단계

다음에 할 일에 대한 상위 3가지 제안은 다음과 같습니다:

- [Nextflow for Science](../nf4_science/index.md)로 과학 분석 사용 사례에 Nextflow 적용하기
- [Hello nf-core](../hello_nf-core/index.md)로 nf-core 시작하기
- [Side Quests](../side_quests/index.md)로 더 고급 Nextflow 기능 탐색하기

마지막으로, Nextflow 제작자가 개발한 클라우드 기반 플랫폼인 [**Seqera Platform**](https://seqera.io/)을 살펴보시기를 권장합니다. 이 플랫폼은 워크플로우를 더욱 쉽게 실행하고 관리할 수 있게 하며, 데이터를 관리하고 모든 환경에서 대화형으로 분석을 실행할 수 있게 합니다.

## 피드백 설문조사

다음 단계로 넘어가기 전에, 잠시 시간을 내어 과정 설문조사를 완료해 주세요! 여러분의 피드백은 모두를 위한 교육 자료를 개선하는 데 도움이 됩니다.

[설문조사 참여하기 :material-arrow-right:](survey.md){ .md-button .md-button--primary }
