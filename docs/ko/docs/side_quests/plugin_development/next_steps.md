# 다음 단계

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } AI 지원 번역 - [자세히 알아보기 및 개선 제안](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

**플러그인 개발** 교육 과정을 완료하신 것을 축하드립니다.

---

## 1. 플러그인 여정을 이어가는 3가지 방법

다음 단계로 나아가기 위한 세 가지 권장 사항을 소개합니다.

### 1.1. 플러그인 레지스트리 살펴보기

[Nextflow Plugin Registry](https://registry.nextflow.io/)에서 사용 가능한 플러그인을 확인합니다.
입력 유효성 검사, 클라우드 플랫폼 통합, 알림, 출처 추적 등 다양한 플러그인을 찾아볼 수 있습니다.

### 1.2. 기존 플러그인 소스 코드 분석하기

인기 있는 플러그인의 소스 코드를 살펴보며 구조를 파악합니다.
[nf-hello](https://github.com/nextflow-io/nf-hello) 플러그인은 의도적으로 단순하게 설계되어 있으며 문서화가 잘 되어 있습니다.
[공식 플러그인 저장소](https://github.com/nextflow-io/plugins)에는 executor 및 파일 시스템 구현을 포함한 더 복잡한 예제들이 있습니다.

### 1.3. 유용한 플러그인 직접 만들기

Nextflow에서 필요했지만 찾지 못했던 기능을 생각해 보세요.
몇 가지 아이디어를 제안합니다:

- 조직에 맞는 커스텀 로깅 또는 리포팅
- 내부 도구 또는 API와의 통합
- 도메인 특화 유틸리티 함수
- 선호하는 플랫폼으로의 알림 전송

작게 시작하여 동작하도록 만든 후, 반복적으로 개선합니다.

---

## 2. 커뮤니티에서 도움 받기

- [Nextflow Slack](https://www.nextflow.io/slack-invite.html): 질문하고 작업 결과를 공유합니다.
- [커뮤니티 포럼](https://community.seqera.io/): 아이디어를 논의하고 조언을 구합니다.
- [GitHub discussions](https://github.com/nextflow-io/nextflow/discussions): 기술적인 질문 및 기능 요청을 남깁니다.

유용한 플러그인을 만들었다면, 플러그인 레지스트리를 통해 커뮤니티와 공유하는 것을 고려해 보세요.

---

## 3. Nextflow 교육 계속하기

아직 수강하지 않으셨다면, 다음 교육 과정들을 확인해 보세요:

- **[Hello Nextflow](../../hello_nextflow/index.md)**: Nextflow의 기초 개념
- **[Hello nf-core](../../hello_nf-core/index.md)**: nf-core 파이프라인 및 모범 사례
- **[Side Quests](../index.md)**: 특정 주제에 대한 심층 학습
