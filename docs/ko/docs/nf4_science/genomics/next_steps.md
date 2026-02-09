# 과정 요약

Nextflow for Genomics 교육 과정을 완료하신 것을 축하드립니다! 🎉

## 여러분의 여정

터미널에서 variant calling 도구를 수동으로 실행하여 방법론을 이해하는 것으로 시작했습니다.
그런 다음 프로세스를 자동화하기 위해 단일 샘플 Nextflow 파이프라인을 구축하고, 이를 확장하여 여러 샘플을 병렬로 처리하도록 했으며, 채널 연산자를 사용하여 다중 샘플 joint genotyping을 추가했습니다.

### 구축한 것

- BAM 파일을 입력으로 받아 joint-called VCF를 출력으로 생성하는 variant calling 파이프라인.
- 별도의 모듈 파일에 저장된 세 개의 프로세스(`SAMTOOLS_INDEX`, `GATK_HAPLOTYPECALLER`, `GATK_JOINTGENOTYPING`).
- Nextflow의 데이터 흐름 패러다임을 사용하여 입력 샘플 수에 관계없이 자동으로 확장되는 파이프라인.
- `results/` 디렉토리에 게시되는 결과.

### 습득한 기술

이 실습 과정을 통해 다음 방법을 학습했습니다.

- 단일 샘플에 variant calling을 적용하는 선형 워크플로우 작성
- 인덱스 파일 및 참조 유전체 리소스와 같은 보조 파일을 적절하게 처리
- Nextflow의 데이터 흐름 패러다임을 활용하여 샘플별 variant calling 병렬화
- 관련 채널 연산자를 사용하여 다중 샘플 joint calling 구현

이제 여러분은 자신의 작업에서 Nextflow를 유전체학 분석 워크플로우에 적용하기 시작할 준비가 되었습니다.

## 기술 향상을 위한 다음 단계

다음에 수행할 작업에 대한 주요 권장 사항입니다.

- [Nextflow for Science](../index.md)로 다른 과학적 분석 사용 사례에 Nextflow 적용하기
- [Hello nf-core](../../hello_nf-core/index.md)로 nf-core 시작하기
- [Side Quests](../../side_quests/index.md)로 더욱 고급 Nextflow 기능 탐색하기

마지막으로, Nextflow의 제작자가 개발한 클라우드 기반 플랫폼인 **[Seqera Platform](https://seqera.io/)**을 살펴보실 것을 권장합니다. 이 플랫폼은 워크플로우를 시작하고 관리하며, 데이터를 관리하고 모든 환경에서 대화형으로 분석을 실행하는 것을 훨씬 쉽게 만들어줍니다.

## 도움 받기

도움 리소스 및 커뮤니티 지원은 [도움말 페이지](../../help.md)를 참조하십시오.

## 피드백 설문조사

다음 단계로 넘어가기 전에, 잠시 시간을 내어 과정 설문조사를 완료해 주십시오! 여러분의 피드백은 모든 사람을 위한 교육 자료를 개선하는 데 도움이 됩니다.

[설문조사 참여하기 :material-arrow-right:](survey.md){ .md-button .md-button--primary }
