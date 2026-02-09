# Nextflow for RNAseq

이 교육 과정은 데이터 분석 파이프라인을 개발하거나 커스터마이징하는 데 관심이 있는 전사체학 및 관련 분야의 연구자를 대상으로 합니다.
이 과정은 [Hello Nextflow](../../hello_nextflow/) 초급 교육을 기반으로 하며, 대량 RNAseq 분석의 특정 맥락에서 Nextflow를 사용하는 방법을 보여줍니다.

구체적으로, 이 과정은 어댑터 서열을 트리밍하고, 리드를 게놈 참조에 정렬하며, 여러 단계에서 품질 관리(QC)를 수행하는 간단한 대량 RNAseq 처리 파이프라인을 구현하는 방법을 보여줍니다.

시작해 봅시다! 아래의 "Open in GitHub Codespaces" 버튼을 클릭하여 교육 환경을 실행하세요(별도의 탭에서 여는 것을 권장합니다). 그런 다음 로딩되는 동안 계속 읽어주세요.

[![Open in GitHub Codespaces](https://github.com/codespaces/badge.svg)](https://codespaces.new/nextflow-io/training?quickstart=1&ref=master)

## 학습 목표

이 과정을 통해 일반적인 RNAseq 사용 사례에 기본적인 Nextflow 개념과 도구를 적용하는 방법을 학습하게 됩니다.

이 워크숍을 마치면 다음을 수행할 수 있습니다:

- 기본 RNAseq 처리 및 QC 방법을 적용하는 선형 워크플로우 작성하기
- FASTQ 및 참조 게놈 리소스와 같은 도메인별 파일을 적절하게 처리하기
- 단일 말단 및 쌍 말단 시퀀싱 데이터 처리하기
- Nextflow의 데이터플로우 패러다임을 활용하여 샘플별 RNAseq 처리를 병렬화하기
- 관련 채널 연산자를 사용하여 여러 단계와 샘플에 걸쳐 QC 보고서 집계하기

<!-- TODO
- Configure pipeline execution and manage and optimize resource allocations
- Implement per-step and end-to-end pipeline tests that handle RNAseq-specific idiosyncrasies appropriately
-->
<!-- TODO for future expansion: add metadata/samplesheet handling -->

## 전제 조건

이 과정은 다음에 대한 최소한의 지식을 가정합니다:

- 이 과학 분야에서 일반적으로 사용되는 도구 및 파일 형식
- 커맨드 라인 사용 경험
- [Hello Nextflow](../../hello_nextflow/) 초급 교육에서 다루는 기본적인 Nextflow 개념 및 도구

기술 요구 사항 및 환경 설정에 대해서는 [환경 설정](../../envsetup/) 단기 과정을 참조하세요.
