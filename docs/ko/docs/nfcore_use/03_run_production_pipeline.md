# 파트 3: 프로덕션 파이프라인 실행

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } AI 지원 번역 - [자세히 알아보기 및 개선 제안](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

[파트 2](./02_configure_execution.md)에서는 nf-core/demo의 매개변수 설정 및 실행 설정 맞춤화 방법을 학습했습니다.
이제 학습한 내용을 실제 프로덕션 파이프라인인 nf-core/rnaseq에 적용합니다.

---

## 1. nf-core/rnaseq 가져오기 및 실행

지금까지는 교육용으로 설계된 최소한의 파이프라인인 `nf-core/demo`를 사용했습니다.
이제 실제 프로덕션 파이프라인을 가져와 테스트 프로파일로 실행합니다.

`nf-core/rnaseq` 파이프라인은 벌크 RNA 시퀀싱 분석의 핵심 단계인 품질 관리, 어댑터 트리밍, 리드 정렬, 유전자 수준 정량화를 수행합니다.
현재까지 가장 널리 사용되는 nf-core 파이프라인입니다.

### 1.1. 파이프라인 가져오기

다음 명령을 실행하여 파이프라인을 다운로드합니다.

```bash
nextflow pull nf-core/rnaseq
```

??? success "명령 출력"

    ```console
    Checking nf-core/rnaseq ...
     downloaded from https://github.com/nf-core/rnaseq.git - revision: 1f03b53ef7 [master]
    ```

파이프라인이 로컬에 캐시되어 실행할 준비가 완료되었습니다.

### 1.2. 테스트 프로파일 실행

테스트 프로파일과 Docker로 실행합니다.

```bash
nextflow run nf-core/rnaseq -profile test,docker --outdir rnaseq-results
```

??? failure "명령 출력"

    ```console
     N E X T F L O W   ~  version 26.04.4

    Launching `https://github.com/nf-core/rnaseq` [suspicious_dijkstra] revision: 1f03b53ef7 [master]

    ------------------------------------------------------
                                            ,--./,-.
            ___     __   __   __   ___     /,-._.--~'
      |\ | |__  __ /  ` /  \ |__) |__         }  {
      | \| |       \__, \__/ |  \ |___     \`-._,-`-,
                                            `._,._,'
      nf-core/rnaseq 3.26.0
    ------------------------------------------------------
    ...
    [-        ] NFCORE_RNASEQ:RNASEQ:FASTQ_QC_TRIM_FILTER_SETSTRANDEDNESS:FQ_LINT  -
    Plus 47 more processes waiting for tasks…

    Execution cancelled -- Finishing pending tasks before exit
    -[nf-core/rnaseq] Pipeline completed with errors-
    ERROR ~ Error executing process > 'NFCORE_RNASEQ:RNASEQ:FASTQ_QC_TRIM_FILTER_SETSTRANDEDNESS:FQ_LINT (WT_REP2)'

    Caused by:
      Process requirement exceeds available memory -- req: 12 GB; avail: 7.7 GB


    Command executed:

      fq lint \
          --disable-validator P001 \
          SRR6357072_1.fastq.gz SRR6357072_2.fastq.gz > WT_REP2.fq_lint.txt

    Command exit status:
      -

    Command output:
      (empty)

    Work dir:
      /workspaces/training/nfcore-use/work/xx/xxxxxxxxxxxxxxxxxxxxxx

    Container:
      quay.io/biocontainers/fq:0.12.0--h9ee0642_0

    Tip: you can replicate the issue by changing to the process work dir and entering the command `bash .command.run`

     -- Check '.nextflow.log' file for details
    ERROR ~ Pipeline failed. Please refer to troubleshooting docs for common issues: https://nf-co.re/docs/running/troubleshooting
    ```

오류의 핵심 내용은 다음과 같습니다.

```console
Process requirement exceeds available memory -- req: 12 GB; avail: 7.7 GB
```

기본 Codespaces 머신은 8 GB의 RAM을 제공하며, 이는 Docker Desktop의 일반적인 기본값이기도 합니다.
파이프라인의 `FQ_LINT` 프로세스가 12 GB를 요청하고 있어 머신이 제공할 수 있는 용량을 초과합니다.

이 12 GB는 `conf/base.config`에 정의된 `process_low` 리소스 레이블에서 비롯됩니다.

```groovy title="conf/base.config"
withLabel:process_low {
    cpus   = { 2     * task.attempt }
    memory = { 12.GB * task.attempt }
    time   = { 4.h   * task.attempt }
}
```

더 큰 머신 유형을 사용하는 방법도 있지만, 테스트 목적으로는 사용 가능한 하드웨어에서 실행할 수 있어야 합니다.
더 나은 방법은 사용자 정의 설정 파일로 기본 리소스 값을 재정의하는 것입니다.

### 1.3. 사용자 정의 설정으로 재실행

레이블 기반 리소스 기본값을 재정의하는 사용자 정의 설정 파일을 제공합니다.

??? full-code "laptop.config"

    ```groovy title="laptop.config"
    process {
        withLabel: 'process_low' {
            cpus   = 2
            memory = 6.GB
        }
        withLabel: 'process_medium' {
            cpus   = 4
            memory = 6.GB
        }
        withLabel: 'process_high' {
            cpus   = 6
            memory = 6.GB
        }
        withLabel: 'process_high_memory' {
            memory = 6.GB
        }
    }
    ```

[파트 2](./02_configure_execution.md)에서는 `withName:`을 사용하여 단일 프로세스를 이름으로 지정하는 방법을 다뤘습니다.
여기서는 `withLabel:`을 사용하여 동일한 레이블을 공유하는 모든 프로세스를 한 번에 지정합니다.

이 파일은 이미 작업 디렉토리에 있습니다.
`-c` 옵션으로 전달하여 재정의를 적용합니다.

```bash
nextflow run nf-core/rnaseq -profile test,docker -c laptop.config --outdir rnaseq-results
```

??? success "명령 출력 (파이프라인 시작)"

    ```console
     N E X T F L O W   ~  version 26.04.4

    Launching `https://github.com/nf-core/rnaseq` [romantic_faraday] revision: 1f03b53ef7 [master]

    ------------------------------------------------------
                                            ,--./,-.
            ___     __   __   __   ___     /,-._.--~'
      |\ | |__  __ /  ` /  \ |__) |__         }  {
      | \| |       \__, \__/ |  \ |___     \`-._,-`-,
                                            `._,._,'
      nf-core/rnaseq 3.26.0
    ------------------------------------------------------
    ...
    executor >  local (7)
    [xx/xxxxxx] NFCORE_RNASEQ:RNASEQ:FQ_LINT (RAP1_IAA_30M_REP1)   | 3 of 5
    [xx/xxxxxx] NFCORE_RNASEQ:RNASEQ:FASTQC (RAP1_IAA_30M_REP1)    | 2 of 5
    ...
    ```

파이프라인이 실행 중이며, 작업이 하나씩 완료되는 것을 확인할 수 있습니다.
이 최소한의 테스트 데이터셋에서는 총 200개 이상의 작업을 실행하며 15~20분 내에 완료됩니다.

실제 RNA-seq 실험은 일반적으로 수십 개의 샘플을 포함하며 몇 시간 또는 며칠이 걸립니다.
Nextflow는 HPC 스케줄러(SLURM, PBS, LSF)와 클라우드 플랫폼(AWS, Google Cloud, Azure)을 지원하며, 여러 노드에 작업을 분산하여 실제 소요 시간을 크게 줄일 수 있습니다.
그러나 이러한 환경을 설정하는 데는 상당한 복잡성이 따릅니다.

Seqera 플랫폼(Nextflow 개발팀이 만든 서비스)은 HPC 또는 클라우드 인프라에서 Nextflow 파이프라인을 실행하기 위한 웹 기반 인터페이스를 제공합니다. 자체 인프라 또는 관리형 인프라를 사용할 수 있으며, 대규모 파이프라인 실행 과정을 간소화하는 컴퓨팅 및 데이터 관리 기능을 갖추고 있습니다.

!!! tip "팁"

    학술 연구자는 [Seqera 학술 프로그램](https://seqera.io/academic-program/)을 통해 Seqera Platform을 무료로 이용할 수 있습니다.

### 핵심 정리

`nf-core/rnaseq`를 가져오고, nf-core 리소스 레이블의 작동 방식을 확인했으며, 사용자 정의 설정 파일로 이를 재정의하는 방법을 학습했습니다.
더 중요한 것은, 로컬 실행이 실제 규모의 분석에서는 출발점일 뿐 최종 목적지가 아님을 확인했다는 점입니다.

### 다음 단계

nf-core 파이프라인 실행의 기본 사항을 모두 다뤘습니다.
앞으로 나아갈 방향은 [다음 단계](next_steps.md)를 참조하세요.

---

## 요약

이 파트에서는 다음을 학습했습니다.

- 프로덕션 규모의 파이프라인(nf-core/rnaseq)을 가져와 실행하고, 기본 리소스 레이블을 재정의하는 방법
