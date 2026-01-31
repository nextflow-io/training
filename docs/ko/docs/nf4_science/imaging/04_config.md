# 파트 4: 구성

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } AI 지원 번역 - [자세히 알아보기 및 개선 제안](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

파트 1-3에서는 Nextflow를 실행하고, nf-core 파이프라인을 실행하며, 매개변수 파일과 샘플시트로 입력을 관리하는 방법을 배웠습니다.
이제 **구성 파일**과 **프로파일**을 사용하여 다양한 컴퓨팅 환경에 맞게 파이프라인을 구성하는 방법을 살펴보겠습니다.

## 학습 목표

이 파트를 마치면 다음을 수행할 수 있습니다:

- 여러 소스에서 Nextflow가 구성을 해석하는 방법 이해하기
- 컨테이너 및 테스트를 위한 nf-core 내장 프로파일 사용하기
- 다양한 컴퓨팅 환경을 위한 사용자 정의 프로파일 생성하기
- 프로세스 레이블을 사용하여 리소스 요청 사용자 정의하기
- 제한된 환경에서 리소스 제한 관리하기
- `nextflow config`로 해석된 구성 검사하기

---

## 1. Nextflow 구성 이해하기

### 1.1. 구성 파일이란 무엇인가요?

Nextflow는 구성 파일을 사용하여 **워크플로우 로직**(무엇을 할 것인가)과 **실행 설정**(어떻게 그리고 어디서 할 것인가)을 분리합니다.

구성 파일이 제어하는 것:

- 컨테이너 엔진(Docker, Singularity, Conda)
- 컴퓨팅 리소스(CPU, 메모리, 시간)
- 실행 플랫폼(로컬, HPC, 클라우드)
- 파이프라인 매개변수

### 1.2. 구성 우선순위

Nextflow는 여러 소스에서 구성을 로드하며, 나중에 로드된 소스가 이전 소스를 재정의합니다:

1. **파이프라인 구성**: 파이프라인 저장소의 `nextflow.config`
2. **디렉토리 구성**: 현재 작업 디렉토리의 `nextflow.config`
3. **사용자 구성**: `~/.nextflow/config`
4. **명령줄**: 직접 전달된 매개변수 및 옵션

이러한 계층적 접근 방식을 통해 파이프라인의 기본값을 유지하고, 사용자별 설정으로 재정의하며, 명령줄에서 빠른 조정을 할 수 있습니다.

### 1.3. 현재 구성

지금까지 사용한 구성을 살펴보겠습니다:

```groovy title="nextflow.config"
docker.enabled = true
process {
    resourceLimits = [
        cpus: 2,
        memory: '7.GB',
    ]
}

```

파트 2에서 추가한 `docker.enabled = true` 줄을 주석 처리하거나 다시 변경하고, 대신 molkart에서 프로파일을 사용하여 동일한 결과를 얻을 수 있는 방법을 알아보겠습니다.

---

## 2. 프로파일 사용하기

### 2.1. 프로파일이란 무엇인가요?

프로파일은 `nextflow run` 명령을 통해 `-profile` 플래그로 활성화할 수 있는 명명된 구성 세트입니다.
프로파일을 사용하면 구성 파일을 편집하지 않고도 다양한 컴퓨팅 시나리오 간에 쉽게 전환할 수 있습니다.

모든 nf-core 파이프라인에는 우리가 사용할 수 있는 여러 기본 프로파일이 함께 제공됩니다.

### 2.2. 내장 프로파일 검사하기

파이프라인 코드베이스와 관련된 `molkart/nextflow.config` 파일에서 이들을 검사해 보겠습니다:

```bash
code molkart/nextflow.config
```

`profiles` 블록을 검색하십시오:

```groovy title="molkart/nextflow.config (발췌)"
profiles {
    docker {
        docker.enabled          = true
        singularity.enabled     = false
        conda.enabled           = false
    }
    singularity {
        singularity.enabled     = true
        singularity.autoMounts  = true
        docker.enabled          = false
    }
    conda {
        conda.enabled           = true
        docker.enabled          = false
        conda.channels          = ['conda-forge', 'bioconda']
    }
}
```

일반적인 컨테이너 프로파일:

- `docker`: Docker 컨테이너 사용(로컬 개발에 가장 일반적)
- `singularity`: Singularity/Apptainer 사용(HPC에서 일반적)
- `conda`: Conda 환경 사용
- `apptainer`: Apptainer 컨테이너 사용

### 2.3. nextflow.config 대신 프로파일로 다시 실행하기

이제 로컬 `nextflow.config` 파일에서 docker 구성을 비활성화하고 프로파일을 이해했으므로, `-profile` 플래그를 사용하여 파이프라인을 다시 실행해 보겠습니다.

이전 파트 3에서는 사용자 정의 매개변수가 포함된 `params.yaml` 파일을 만들었습니다.
이제 이를 내장 Docker 프로파일과 결합할 수 있습니다:

```bash
nextflow run ./molkart \
  -profile docker \
  -params-file params.yaml \
  -resume
```

각 플래그가 수행하는 작업을 살펴보겠습니다:

- `-profile docker`: molkart의 `nextflow.config`에서 Docker 프로파일을 활성화하여 `docker.enabled = true`로 설정합니다
- `-params-file params.yaml`: YAML 파일에서 모든 파이프라인 매개변수를 로드합니다
- `-resume`: 이전 실행에서 캐시된 결과를 재사용합니다

`-resume`을 사용하기 때문에 Nextflow는 마지막 실행 이후 변경 사항이 있는지 확인합니다.
매개변수, 입력 및 코드가 동일하면 모든 작업이 캐시에서 검색되고 파이프라인이 거의 즉시 완료됩니다.

```console title="출력 (발췌)"
executor >  local (12)
...
[1a/2b3c4d] NFCORE_MOLKART:MOLKART:MINDAGAP_MINDAGAP (mem_only)   [100%] 2 of 2, cached: 2 ✔
[5e/6f7g8h] NFCORE_MOLKART:MOLKART:CLAHE (mem_only)               [100%] 2 of 2, cached: 2 ✔
...
-[nf-core/molkart] Pipeline completed successfully-
```

모든 프로세스가 `cached: 2` 또는 `cached: 1`을 표시하는 것을 확인하십시오 - 아무것도 재실행되지 않았습니다!

### 2.4. 테스트 프로파일

테스트 프로파일은 파이프라인이 작동하는지 확인할 수 있도록 기본 입력 매개변수와 데이터 파일을 지정하는 빠른 방법을 제공합니다.
nf-core 파이프라인은 항상 최소 두 개의 테스트 프로파일을 포함합니다:

- `test`: 빠른 테스트를 위한 작은 데이터셋과 빠른 매개변수
- `test_full`: 더 큰 데이터로 더 포괄적인 테스트

`includeConfig` 지시문을 사용하여 포함되는 molkart의 `test` 프로파일을 자세히 살펴보겠습니다:

```groovy title="molkart/nextflow.config (발췌)"
profiles {
  ...
    test      { includeConfig 'conf/test.config'      }
}
```

이는 `-profile test`로 파이프라인을 실행할 때마다 Nextflow가 `conf/test.config`에서 구성을 로드한다는 의미입니다.

```groovy title="molkart/conf/test.config (발췌)"
params {
    config_profile_name        = 'Test profile'
    config_profile_description = 'Minimal test dataset to check pipeline function'

    input = 'https://raw.githubusercontent.com/nf-core/test-datasets/molkart/test_data/samplesheets/samplesheet_membrane.csv'
    mindagap_tilesize = 90
    mindagap_boxsize = 7
    mindagap_loopnum = 100
    clahe_pyramid_tile = 368
    segmentation_method = "mesmer,cellpose,stardist"
}

process {
    resourceLimits = [
        cpus: 4,
        memory: '15.GB',
        time: '1.h'
    ]
}
```

이 프로파일에는 이전에 `params.yaml` 파일에서 사용했던 것과 동일한 매개변수가 포함되어 있습니다.

쉼표로 구분하여 여러 프로파일을 활성화할 수 있습니다.
이를 사용하여 params 파일 없이 파이프라인을 테스트해 보겠습니다:

```bash
nextflow run ./molkart -profile docker,test --outdir results -resume
```

이것은 다음을 결합합니다:

- `docker`: Docker 컨테이너 활성화
- `test`: 테스트 데이터셋 및 매개변수 사용

프로파일은 왼쪽에서 오른쪽으로 적용되므로, 동일한 값을 설정하는 경우 나중 프로파일이 이전 프로파일을 재정의합니다.

### 요점

nf-core 파이프라인에는 컨테이너, 테스트 및 특수 환경을 위한 내장 프로파일이 함께 제공됩니다.
여러 프로파일을 결합하여 필요한 구성을 만들 수 있습니다.

### 다음 단계

다양한 컴퓨팅 환경을 위한 자체 사용자 정의 프로파일을 만드는 방법을 배웁니다.

---

## 3. 사용자 정의 프로파일 만들기

### 3.1. 로컬 개발과 HPC 실행 간 전환을 위한 프로파일 만들기

두 가지 시나리오에 대한 사용자 정의 프로파일을 만들어 보겠습니다:

1. Docker를 사용한 로컬 개발
2. Slurm 스케줄러와 Singularity를 사용하는 대학 HPC

`nextflow.config`에 다음을 추가하십시오:

```groovy title="nextflow.config"
profiles {
    local_dev {
        docker.enabled = true
        process.executor = 'local'
    }

    hpc_cluster {
        singularity.enabled = true
        process.executor = 'slurm'
        process.queue = 'standard_queue'
        singularity.cacheDir = '/shared/containers'
    }
}
```

이제 환경 간에 쉽게 전환할 수 있습니다:

```bash
# 로컬 개발용
nextflow run ./molkart -profile local_dev --input data/samplesheet.csv --outdir results

# HPC용 (사용 가능한 경우)
nextflow run ./molkart -profile hpc_cluster --input data/samplesheet.csv --outdir results
```

!!! Note "참고"

    Slurm 스케줄러에 접근할 수 없으므로 이 교육 환경에서는 HPC 프로파일을 테스트할 수 없습니다.
    하지만 이것은 실제 사용을 위해 구성하는 방법을 보여줍니다.

### 3.2. `nextflow config`를 사용하여 구성 검사하기

`nextflow config` 명령은 파이프라인을 실행하지 않고 완전히 해석된 구성을 보여줍니다.

기본 구성 보기:

```bash
nextflow config ./molkart
```

특정 프로파일로 구성 보기:

```bash
nextflow config -profile local_dev ./molkart
```

이것은 다음에 매우 유용합니다:

- 구성 문제 디버깅
- 실제로 사용될 값 이해
- 여러 프로파일이 상호 작용하는 방식 확인

### 요점

사용자 정의 프로파일을 사용하면 단일 명령줄 플래그로 다양한 컴퓨팅 환경 간에 전환할 수 있습니다.
`nextflow config`를 사용하여 실행 전에 해석된 구성을 검사하십시오.

### 다음 단계

nf-core의 프로세스 레이블 시스템을 사용하여 개별 프로세스에 대한 리소스 요청을 사용자 정의하는 방법을 배웁니다.

---

## 4. 리소스 요청 사용자 정의하기

### 4.1. nf-core 파이프라인의 프로세스 레이블 이해하기

간단히 하기 위해, nf-core 파이프라인은 [**프로세스 레이블**](https://www.nextflow.io/docs/latest/reference/process.html#process-label)을 사용하여 모든 파이프라인에서 리소스 할당을 표준화합니다.
각 프로세스에는 낮음, 중간 또는 높은 컴퓨팅 리소스 요구 사항을 설명하는 `process_low`, `process_medium` 또는 `process_high`와 같은 레이블이 태그됩니다.
이러한 레이블은 파이프라인의 `conf/` 디렉토리에 있는 구성 파일 중 하나에서 특정 리소스 요청으로 변환됩니다.

```groovy title="molkart/conf/base.config (발췌)"
process {
    cpus   = { 1      * task.attempt }
    memory = { 6.GB   * task.attempt }
    time   = { 4.h    * task.attempt }

    errorStrategy = { task.exitStatus in ((130..145) + 104) ? 'retry' : 'finish' }
    maxRetries    = 1
    maxErrors     = '-1'

    withLabel:process_single {
        cpus   = { 1                   }
        memory = { 6.GB * task.attempt }
        time   = { 4.h  * task.attempt }
    }
    withLabel:process_low {
        cpus   = { 2     * task.attempt }
        memory = { 12.GB * task.attempt }
        time   = { 4.h   * task.attempt }
    }
    withLabel:process_medium {
        cpus   = { 6     * task.attempt }
        memory = { 36.GB * task.attempt }
        time   = { 8.h   * task.attempt }
    }
    withLabel:process_high {
        cpus   = { 12    * task.attempt }
        memory = { 72.GB * task.attempt }
        time   = { 16.h  * task.attempt }
    }
}
```

`task.attempt` 승수를 주목하십시오 - 이를 통해 파이프라인이 `process.maxRetries > 1`로 설정된 경우 후속 작업 재시도가 더 많은 리소스를 요청할 수 있습니다.

### 4.2. 특정 프로세스에 대한 리소스 재정의

세밀한 제어를 위해 이름으로 개별 프로세스를 대상으로 지정하십시오:

```groovy title="nextflow.config"
process {
    withName: 'NFCORE_MOLKART:MOLKART:CELLPOSE' {
        cpus   = 16
        memory = 32.GB
    }
}
```

위의 재정의로 이 파이프라인을 실행하려고 하면 `CELLPOSE` 프로세스는 해당 레이블로 정의된 기본값 대신 16개의 CPU와 32GB 메모리를 요청합니다.
사용 가능한 RAM이 충분하지 않기 때문에 현재 환경에서는 파이프라인이 실패합니다.
다음 섹션에서 이러한 유형의 실패를 방지하는 방법을 배웁니다.

!!! Tip "팁"

    프로세스 이름을 찾으려면 파이프라인 실행 출력을 살펴보거나 `.nextflow.log`를 확인하십시오.
    프로세스 이름은 `WORKFLOW:SUBWORKFLOW:PROCESS` 패턴을 따릅니다.

### 요점

nf-core 파이프라인은 프로세스 레이블을 사용하여 리소스 할당을 표준화합니다.
레이블별로(여러 프로세스에 영향) 또는 이름별로(특정 프로세스 하나에 영향) 리소스를 재정의할 수 있습니다.

### 다음 단계

GitHub Codespaces와 같이 제한된 환경에서 리소스 제한을 관리하는 방법을 배웁니다.

---

## 5. 제한된 환경에서 리소스 관리하기

### 5.1. 리소스 제한 문제

섹션 4.2에 표시된 것처럼 16개의 CPU와 32GB 메모리를 요청하는 프로세스로 molkart를 실행하려고 하면 사용 가능한 리소스가 충분하지 않아 현재 환경에서 실패합니다.
더 큰 노드가 있는 클러스터 환경에서는 이러한 요청이 스케줄러에 제출됩니다.

GitHub Codespaces와 같이 제한된 환경에서 제한이 없으면 Nextflow는 사용 가능한 리소스를 초과하는 프로세스의 실행을 거부합니다.

### 5.2. 리소스 제한 설정하기

`resourceLimits` 지시문은 지정된 값으로 리소스 요청을 제한합니다:

```groovy title="nextflow.config"
process {
    resourceLimits = [ cpus: 2, memory: 7.GB ]
}
```

이것은 Nextflow에게 다음과 같이 알려줍니다: "프로세스가 2개 이상의 CPU 또는 7GB 이상의 메모리를 요청하면 대신 이러한 제한으로 제한하십시오."

### 5.3. 사용자 정의 프로파일에 리소스 제한 추가하기

적절한 제한을 포함하도록 사용자 정의 프로파일을 업데이트하십시오:

```groovy title="nextflow.config"
profiles {
    local_dev {
        docker.enabled = true
        process.executor = 'local'
        process.resourceLimits = [
            cpus: 2,
            memory: 7.GB
        ]
    }

    hpc_cluster {
        singularity.enabled = true
        process.executor = 'slurm'
        process.queue = 'batch'
        process.resourceLimits = [
            cpus: 32,
            memory: 128.GB,
            time: 24.h
        ]
    }
}
```

!!! Warning "경고"

    리소스 제한을 너무 낮게 설정하면 프로세스가 실패하거나 느리게 실행될 수 있습니다.
    파이프라인은 메모리 집약도가 낮은 알고리즘을 사용하거나 더 작은 청크로 데이터를 처리해야 할 수 있습니다.

### 요점

프로세스 리소스 요청을 제한하여 리소스 제약 환경에서 파이프라인을 실행하려면 `resourceLimits`를 사용하십시오.
서로 다른 프로파일은 해당 환경에 적합한 서로 다른 제한을 가질 수 있습니다.

### 다음 단계

바이오이미징을 위한 핵심 Nextflow 교육을 완료했습니다!

---

## 결론

이제 다양한 컴퓨팅 환경에 맞게 Nextflow 파이프라인을 구성하는 방법을 이해했습니다.

배운 주요 기술:

- **구성 우선순위**: 여러 소스에서 Nextflow가 설정을 해석하는 방법
- **nf-core 프로파일**: 컨테이너, 테스트 및 유틸리티를 위한 내장 프로파일 사용
- **사용자 정의 프로파일**: 다양한 환경을 위한 자체 프로파일 생성
- **프로세스 레이블**: 레이블별 리소스 요청 이해 및 재정의
- **리소스 제한**: `resourceLimits`로 제약된 환경 관리
- **구성 검사**: 설정을 디버그하고 확인하기 위해 `nextflow config` 사용

이러한 구성 기술은 모든 Nextflow 파이프라인에 적용 가능하며 로컬 머신, HPC 클러스터 및 클라우드 플랫폼에서 워크플로우를 효율적으로 실행하는 데 도움이 됩니다.

### 다음 단계

바이오이미징을 위한 Nextflow 과정을 완료한 것을 축하합니다!

다음 단계:

- 피드백을 제공하기 위해 과정 설문조사를 작성하십시오
- 워크플로우 개발에 대해 더 자세히 알아보려면 [Hello Nextflow](../hello_nextflow/index.md)를 확인하십시오
- nf-core 도구에 대해 더 깊이 이해하려면 [Hello nf-core](../hello_nf-core/index.md)를 살펴보십시오
- [교육 컬렉션](../training_collections/index.md)에서 다른 과정을 찾아보십시오
