# 교육 과정 요약

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } AI 지원 번역 - [자세히 알아보기 및 개선 제안](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

Hello nf-core 교육 과정을 완료하신 것을 축하드립니다! 🎉

<!-- placeholder for video -->

## 학습 여정

데모 pipeline을 가져와서 실행하는 방법을 배우는 것으로 시작하여, 간단한 Nextflow workflow를 nf-core pipeline으로 변환하는 작업을 수행하셨습니다.
template을 사용하여 pipeline scaffold를 생성하는 방법을 배우고 기존 pipeline을 해당 scaffold에 접목시켰습니다.
그런 다음 local module 중 하나를 nf-core module로 교체하고, 다른 local module을 nf-core 표준에 맞게 변환하며, 입력 검증을 추가하여 pipeline을 단계적으로 개선하셨습니다.

### 구축한 결과물

최종 `core-hello` pipeline은 다음을 갖추고 있습니다:

- **표준화된 구조**: nf-core template을 사용하여 workflow, subworkflow, module 및 구성을 위한 체계적인 디렉토리 구조
- **커뮤니티 module**: nf-core 저장소의 module(`cat/cat`)과 함께 사용자 정의 module
- **포괄적인 검증**: pipeline 실행 전에 매개변수와 입력 데이터를 모두 확인
- **전문적인 구성**: 다양한 실행 환경을 위한 프로파일
- **완전한 문서화**: nf-core 규칙을 따르는 메타데이터

### 습득한 핵심 기술

이 실습 교육 과정을 통해 다음을 배우셨습니다:

1. 기존 pipeline을 탐색하여 nf-core pipeline 구조를 **탐색하고 이해하기**
2. nf-core template 내에서 조합 가능하도록 workflow를 **재구성하기**
3. 커뮤니티 저장소에서 사전 구축된 module을 **찾아서 통합하기**
4. 명명, 구조 및 메타데이터에 대한 nf-core 표준을 따라 **사용자 정의 module 생성하기**
5. nf-schema를 사용하여 명확한 피드백으로 오류를 조기에 발견하는 **검증 구현하기**

이제 커뮤니티 모범 사례를 따르는 프로덕션 준비가 된 nf-core pipeline을 구축하기 위한 기초 지식을 갖추셨습니다.

## 기술 향상을 위한 다음 단계

다음 단계로 추천하는 상위 3가지 사항입니다:

- [Nextflow for Science](../nf4_science/index.md)를 통해 과학적 분석 사용 사례에 Nextflow 적용하기
- [Side Quests](../side_quests/index.md)를 통해 더 고급 Nextflow 기능 탐색하기
- [nf-core 커뮤니티 참여](https://nf-co.re/join)를 통해 참여하기

마지막으로, Nextflow 제작자들이 개발한 클라우드 기반 플랫폼인 [**Seqera Platform**](https://seqera.io/)을 살펴보시길 권장합니다. 이 플랫폼을 사용하면 workflow를 더 쉽게 시작하고 관리할 수 있으며, 데이터를 관리하고 모든 환경에서 대화형으로 분석을 실행할 수 있습니다.

## 피드백 설문조사

다음 단계로 넘어가기 전에, 교육 과정 설문조사를 작성하는 데 1분만 투자해 주십시오! 여러분의 피드백은 모두를 위한 교육 자료를 개선하는 데 도움이 됩니다.

[설문조사 참여하기 :material-arrow-right:](survey.md){ .md-button .md-button--primary }
