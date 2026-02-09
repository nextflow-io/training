# 과정 요약

Nextflow Run 교육 과정을 완료하신 것을 축하합니다! 🎉

<!-- placeholder for video -->

## 여러분의 학습 여정

여러분은 매우 기본적인 워크플로우로 시작하여 이를 실행하고, 출력을 찾고, 실행을 관리하는 방법을 학습했습니다.
그런 다음 점점 더 복잡한 버전의 워크플로우를 다루면서 채널과 연산자, 코드 모듈화, 컨테이너를 포함하여 Nextflow 파이프라인을 구동하는 핵심 개념과 메커니즘을 인식하는 방법을 학습했습니다.
마지막으로 여러분의 선호도와 컴퓨팅 인프라에 맞게 파이프라인의 구성을 사용자 정의하는 방법을 학습했습니다.

### 학습한 내용

이제 여러분은 Hello 파이프라인의 실행을 관리하고, 구조를 설명하며, 관련된 주요 코드 부분을 식별할 수 있습니다.

- Hello 워크플로우의 최종 형태는 텍스트 인사말이 포함된 CSV 파일을 입력으로 받습니다.
- 네 단계는 별도의 모듈 파일에 저장된 Nextflow 프로세스(`sayHello`, `convertToUpper`, `collectGreetings`, `cowpy`)로 구현됩니다.
- 결과는 `results/`라는 디렉토리에 게시됩니다.
- 파이프라인의 최종 출력은 대문자로 변환된 인사말을 말하는 캐릭터의 ASCII 아트가 포함된 일반 텍스트 파일입니다.

<figure class="excalidraw">
--8<-- "docs/en/docs/hello_nextflow/img/hello_pipeline_complete.svg"
</figure>

1. **`sayHello`:** 각 인사말을 자체 출력 파일에 작성합니다(_예:_ "Hello-output.txt")
2. **`convertToUpper`:** 각 인사말을 대문자로 변환합니다(_예:_ "HELLO")
3. **`collectGreetings`:** 모든 대문자 인사말을 단일 배치 파일로 수집합니다
4. **`cowpy`:** `cowpy` 도구를 사용하여 ASCII 아트를 생성합니다

워크플로우 구성은 유연하고 재현 가능한 방식으로 입력과 매개변수를 제공하는 것을 지원합니다.

### 습득한 기술

이 실습 과정을 통해 다음 방법을 학습했습니다:

- 로컬에서 Nextflow 워크플로우 실행하기
- Nextflow가 생성한 출력(결과) 및 로그 파일 찾기 및 해석하기
- 간단한 다단계 워크플로우를 구성하는 핵심 Nextflow 구성 요소 인식하기
- 연산자 및 채널 팩토리와 같은 다음 단계 개념 설명하기
- 다양한 컴퓨팅 환경에 맞게 파이프라인 구성하기

이제 여러분은 기존 Nextflow 파이프라인을 자신의 작업에 통합하기 시작할 수 있는 기초 지식을 갖추었습니다.

## 기술을 향상시키기 위한 다음 단계

다음에 할 일에 대한 주요 제안 사항은 다음과 같습니다:

- Nextflow를 실행만 하지 말고 작성해 보세요! [Hello Nextflow](../hello_nextflow/index.md)로 Nextflow 개발자가 되어 보세요
- [Nextflow for Science](../nf4_science/index.md)로 과학적 분석 사례에 Nextflow를 적용해 보세요
- [Hello nf-core](../hello_nf-core/index.md)로 nf-core를 시작해 보세요
- [Debugging Side Quest](../side_quests/debugging.md)로 문제 해결 기법을 학습해 보세요

마지막으로, Nextflow 제작자가 개발한 클라우드 기반 플랫폼인 [**Seqera Platform**](https://seqera.io/)을 살펴보시기를 권장합니다. 이 플랫폼은 워크플로우를 더욱 쉽게 실행하고 관리할 수 있게 해주며, 데이터를 관리하고 모든 환경에서 대화형으로 분석을 실행할 수 있게 해줍니다.

## 도움 받기

도움 리소스 및 커뮤니티 지원에 대해서는 [도움말 페이지](../help.md)를 참조하세요.

## 피드백 설문조사

다음 단계로 넘어가기 전에 잠시 시간을 내어 과정 설문조사를 완료해 주세요! 여러분의 피드백은 모든 사람을 위한 교육 자료를 개선하는 데 도움이 됩니다.

[설문조사 참여하기 :material-arrow-right:](survey.md){ .md-button .md-button--primary }
