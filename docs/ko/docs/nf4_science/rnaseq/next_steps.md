# 과정 요약

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } AI 지원 번역 - [자세히 알아보기 및 개선 제안](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

Nextflow for RNAseq 교육 과정을 완료하신 것을 축하드립니다!

## 여러분의 여정

터미널에서 RNAseq 처리 도구를 수동으로 실행하여 방법론을 이해하는 것으로 시작했습니다.
그런 다음 프로세스를 자동화하기 위해 단일 샘플 Nextflow 파이프라인을 구축하고, 여러 샘플을 병렬로 처리하도록 확장했으며, paired-end 데이터를 처리하고 샘플 전체에 걸쳐 QC 보고서를 집계하도록 확장했습니다.

### 구축한 것

- FASTQ 파일을 입력으로 받아 트리밍된 리드, 정렬 및 집계된 QC 보고서를 출력으로 생성하는 RNAseq 처리 파이프라인
- 별도의 모듈 파일에 저장된 트리밍(Trim Galore), 정렬(HISAT2), 품질 관리(FastQC) 및 보고서 집계(MultiQC)를 위한 프로세스
- Nextflow의 데이터플로우 패러다임을 사용하여 입력 샘플 처리를 자동으로 병렬화하는 파이프라인
- paired-end 시퀀싱 데이터를 처리하는 최종 파이프라인

### 습득한 기술

이 실습 과정을 통해 다음 방법을 배웠습니다:

- 기본 RNAseq 처리 및 QC 방법을 적용하기 위한 선형 워크플로우 작성
- FASTQ 및 참조 게놈 리소스와 같은 도메인별 파일을 적절하게 처리
- single-end 및 paired-end 시퀀싱 데이터 처리
- Nextflow의 데이터플로우 패러다임을 활용하여 샘플별 RNAseq 처리 병렬화
- 관련 채널 연산자를 사용하여 여러 단계 및 샘플에 걸쳐 QC 보고서 집계

이제 여러분은 자신의 작업에서 RNAseq 분석 워크플로우에 Nextflow를 적용할 준비가 되었습니다.

## 기술 향상을 위한 다음 단계

다음에 수행할 작업에 대한 주요 권장 사항은 다음과 같습니다:

- [Nextflow for Science](../index.md)를 통해 다른 과학적 분석 사용 사례에 Nextflow 적용
- [Hello nf-core](../../hello_nf-core/index.md)로 nf-core 시작하기
- [Side Quests](../../side_quests/index.md)로 더 고급 Nextflow 기능 탐색

마지막으로, Nextflow 창시자들이 개발한 클라우드 기반 플랫폼인 [**Seqera Platform**](https://seqera.io/)을 살펴보시기를 권장합니다. 이 플랫폼은 워크플로우를 더욱 쉽게 시작하고 관리할 수 있으며, 모든 환경에서 데이터를 관리하고 대화형으로 분석을 실행할 수 있습니다.

## 도움 받기

도움 리소스 및 커뮤니티 지원에 대해서는 [도움말 페이지](../../help.md)를 참조하십시오.

## 피드백 설문조사

다음 단계로 넘어가기 전에 잠시 시간을 내어 과정 설문조사를 완료해 주십시오! 여러분의 피드백은 모두를 위한 교육 자료를 개선하는 데 도움이 됩니다.

[설문조사 참여하기 :material-arrow-right:](survey.md){ .md-button .md-button--primary }
