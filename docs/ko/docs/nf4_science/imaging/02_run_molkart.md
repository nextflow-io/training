# 파트 2: nf-core/molkart 실행하기

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } AI 지원 번역 - [자세히 알아보기 및 개선 제안](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

Part 1에서는 간단한 Hello World 워크플로를 실행하여 Nextflow 실행의 기본 사항을 이해했습니다.
이제 실제 바이오이미징 파이프라인인 **nf-core/molkart**를 실행하겠습니다.

이 파이프라인은 Resolve Bioscience의 Molecular Cartography 공간 전사체학 데이터를 처리합니다.
하지만 여기서 배우는 Nextflow 패턴은 모든 nf-core 파이프라인이나 프로덕션 워크플로에 적용됩니다.

## 1. nf-core 파이프라인 이해하기

파이프라인을 실행하기 전에 nf-core가 무엇이며 워크플로 실행에 왜 중요한지 이해하겠습니다.

### 1.1. nf-core란 무엇인가?

[nf-core](https://nf-co.re/)는 고품질 Nextflow 파이프라인의 커뮤니티 기반 컬렉션입니다.
모든 nf-core 파이프라인은 동일한 구조와 규칙을 따르므로, 하나를 실행하는 방법을 배우면 모든 파이프라인을 실행할 수 있습니다.

nf-core 파이프라인의 주요 특징:

- **표준화된 구조**: 모든 파이프라인이 일관된 매개변수 이름과 사용 패턴을 가짐
- **내장된 테스트 데이터**: 모든 파이프라인에 빠른 검증을 위한 테스트 프로파일 포함
- **포괄적인 문서**: 상세한 사용 지침 및 매개변수 설명
- **품질 관리**: MultiQC를 사용한 자동 QC 보고서
- **컨테이너 지원**: 재현성을 위해 미리 빌드된 컨테이너 제공

!!! tip "nf-core에 대해 더 알고 싶으신가요?"

    nf-core 파이프라인 개발에 대한 심층적인 소개를 보려면 [Hello nf-core](../../hello_nf-core/index.md) 교육 과정을 확인하세요.
    처음부터 nf-core 파이프라인을 생성하고 사용자 정의하는 방법을 다룹니다.

### 1.2. molkart 파이프라인

![nf-core/molkart 파이프라인](img/molkart.png)

[nf-core/molkart](https://nf-co.re/molkart) 파이프라인은 공간 전사체학 이미징 데이터를 여러 단계로 처리합니다:

1. **이미지 전처리**: 그리드 패턴 채우기 및 선택적 대비 향상
2. **세포 분할**: 다중 알고리즘 옵션 (Cellpose, Mesmer, ilastik, Stardist)
3. **스팟 할당**: 전사체 스팟을 분할된 세포에 할당
4. **품질 관리**: 포괄적인 QC 보고서 생성

주요 출력물:

- 세포별 전사체 카운트 테이블
- 분할 마스크
- MultiQC 품질 관리 보고서

---

## 2. 테스트 데이터로 molkart 실행하기

시작하기 전에 molkart 저장소를 로컬에 클론하여 코드를 검사할 수 있도록 하겠습니다:

```bash
cd /workspaces/training/nf4-science/imaging
git clone --branch 1.2.0 --depth 1 https://github.com/nf-core/molkart
```

이렇게 하면 전체 파이프라인 소스 코드가 포함된 `molkart/` 디렉토리가 생성됩니다.

!!! note "왜 로컬에 클론하나요?"

    일반적으로 `nextflow run nf-core/molkart -r 1.2.0`을 사용하여 GitHub에서 직접 nf-core 파이프라인을 실행합니다.
    Nextflow는 요청된 파이프라인 버전을 자동으로 `$HOME/.nextflow/assets/nf-core/molkart`에 다운로드하고 거기서 실행합니다.
    하지만 이 교육에서는 코드를 더 쉽게 검사할 수 있도록 파이프라인을 다른 로컬 디렉토리에 클론하고 있습니다.

### 2.1. 컨테이너 요구사항 이해하기

전체 파이프라인을 실행하기 전에 컨테이너가 nf-core 파이프라인에 왜 필수적인지 알아보겠습니다.

molkart 테스트 구성의 테스트 데이터셋과 매개변수를 사용하여 파이프라인을 실행해 보겠습니다:

```bash
nextflow run ./molkart \
  --input 'data/samplesheet.csv' \
  --mindagap_tilesize 90 \
  --mindagap_boxsize 7 \
  --mindagap_loopnum 100 \
  --clahe_pyramid_tile 368 \
  --segmentation_method "mesmer,cellpose,stardist" \
  --outdir results
```

이 매개변수들을 분석해 보겠습니다:

- `--input`: 샘플 메타데이터가 포함된 샘플시트 경로
- `--mindagap_tilesize`, `--mindagap_boxsize`, `--mindagap_loopnum`: 그리드 패턴 채우기를 위한 매개변수
- `--clahe_pyramid_tile`: 대비 향상을 위한 커널 크기
- `--segmentation_method`: 세포 분할에 사용할 알고리즘
- `--outdir`: 결과를 저장할 위치

!!! Warning "이 명령은 실패할 것입니다 - 의도적입니다!"

    컨테이너가 필요한 이유를 보여주기 위해 의도적으로 컨테이너 없이 실행하고 있습니다.

잠시 후 다음과 같은 오류가 표시됩니다:

??? failure "명령 출력"

    ```console
    ERROR ~ Error executing process > 'NFCORE_MOLKART:MOLKART:MINDAGAP_DUPLICATEFINDER (mem_only)'

    Caused by:
      Process `NFCORE_MOLKART:MOLKART:MINDAGAP_DUPLICATEFINDER (mem_only)` terminated with an error exit status (127)

    Command executed:

      duplicate_finder.py \
          spots.txt \
          90

    Command exit status:
      127

    Command error:
      .command.sh: line 3: duplicate_finder.py: command not found
    ```

**무슨 일이 일어나고 있나요?**

`command not found` 오류 (종료 상태 127)는 Nextflow가 `duplicate_finder.py`를 실행하려고 했지만 시스템에서 찾을 수 없다는 의미입니다.
이는 다음과 같은 이유 때문입니다:

1. 파이프라인이 특수한 바이오인포매틱스 소프트웨어가 설치되어 있을 것으로 예상합니다
2. 이러한 도구들 (예: `duplicate_finder.py`, `apply_clahe.dask.py` 등)은 표준 Linux 배포판에 포함되어 있지 않습니다
3. 컨테이너 없이 Nextflow는 로컬 머신에서 직접 명령을 실행하려고 합니다

**이 도구들은 어디서 와야 하나요?**

프로세스 모듈 중 하나를 검사하여 소프트웨어 요구사항을 어떻게 선언하는지 살펴보겠습니다.

CLAHE 전처리 모듈을 열어보세요:

```bash
code molkart/modules/local/clahe/main.nf
```

5번째 줄을 보면 다음과 같이 표시됩니다:

```groovy
container 'ghcr.io/schapirolabor/molkart-local:v0.0.4'
```

이 줄은 Nextflow에게 다음과 같이 알려줍니다: "이 프로세스를 실행하려면 필요한 모든 소프트웨어가 포함된 Docker 이미지 `ghcr.io/schapirolabor/molkart-local:v0.0.4`를 사용하세요."

각 프로세스는 필요한 도구를 제공하는 컨테이너 이미지를 선언합니다.
하지만 Nextflow는 사용자가 지시할 때만 이러한 컨테이너를 사용합니다!

**해결책: 구성에서 Docker 활성화**

### 2.2. Docker 구성 및 파이프라인 실행

Docker를 활성화하려면 `nextflow.config` 파일에서 `docker.enabled`를 `false`에서 `true`로 변경해야 합니다.

구성 파일을 여세요:

```bash
code nextflow.config
```

`docker.enabled = false`를 `docker.enabled = true`로 변경하세요:

```groovy
docker.enabled = true
process {
    resourceLimits = [
        cpus: 2,
        memory: '7.GB',
    ]
}
```

이제 동일한 명령으로 파이프라인을 다시 실행하세요:

```bash
nextflow run ./molkart \
  --input 'data/samplesheet.csv' \
  --mindagap_tilesize 90 \
  --mindagap_boxsize 7 \
  --mindagap_loopnum 100 \
  --clahe_pyramid_tile 368 \
  --segmentation_method "cellpose,mesmer,stardist" \
  --outdir results
```

이번에는 Nextflow가:

1. 구성에서 `docker.enabled = true` 설정을 읽습니다
2. 필요한 Docker 이미지를 가져옵니다 (처음만)
3. 각 프로세스를 지정된 컨테이너 내에서 실행합니다
4. 모든 도구가 컨테이너 내에서 사용 가능하므로 성공적으로 실행됩니다

!!! Tip "컨테이너가 중요한 이유"

    대부분의 nf-core 파이프라인은 컨테이너화 (Docker, Singularity, Podman 등)를 **필요로** 합니다. 왜냐하면:

    - 표준 환경에서 사용할 수 없는 특수한 바이오인포매틱스 소프트웨어를 사용합니다
    - 컨테이너는 재현성을 보장합니다 - 정확히 동일한 소프트웨어 버전이 모든 곳에서 실행됩니다
    - 수십 개의 도구와 그 의존성을 수동으로 설치할 필요가 없습니다

    Nextflow의 컨테이너에 대한 자세한 내용은 Hello Nextflow 교육의 [Hello Containers](../../hello_nextflow/05_hello_containers.md)를 참조하세요.

### 2.3. 실행 모니터링

파이프라인이 실행되는 동안 다음과 유사한 출력이 표시됩니다:

??? success "명령 출력"

    ```console
    Nextflow 25.04.8 is available - Please consider updating your version to it

    N E X T F L O W   ~  version 25.04.3

    Launching `https://github.com/nf-core/molkart` [soggy_kalam] DSL2 - revision: 5e54b29cb3 [dev]


    ------------------------------------------------------
                                            ,--./,-.
            ___     __   __   __   ___     /,-._.--~'
      |\ | |__  __ /  ` /  \ |__) |__         }  {
      | \| |       \__, \__/ |  \ |___     \`-._,-`-,
                                            `._,._,'
      nf-core/molkart 1.2.0dev
    ------------------------------------------------------
    Segmentation methods and options
      segmentation_method       : mesmer,cellpose,stardist

    Image preprocessing
      mindagap_boxsize          : 7
      mindagap_loopnum          : 100
      clahe_kernel              : 25
      mindagap_tilesize         : 90
      clahe_pyramid_tile        : 368

    Input/output options
      input                     : https://raw.githubusercontent.com/nf-core/test-datasets/molkart/test_data/samplesheets/samplesheet_membrane.csv
      outdir                    : results

    Institutional config options
      config_profile_name       : Test profile
      config_profile_description: Minimal test dataset to check pipeline function

    Generic options
      trace_report_suffix       : 2025-10-18_22-22-21

    Core Nextflow options
      revision                  : dev
      runName                   : soggy_kalam
      containerEngine           : docker
      launchDir                 : /workspaces/training/nf4-science/imaging
      workDir                   : /workspaces/training/nf4-science/imaging/work
      projectDir                : /workspaces/.nextflow/assets/nf-core/molkart
      userName                  : root
      profile                   : docker,test
      configFiles               :

    !! Only displaying parameters that differ from the pipeline defaults !!
    ------------------------------------------------------
    * The pipeline
        https://doi.org/10.5281/zenodo.10650748

    * The nf-core framework
        https://doi.org/10.1038/s41587-020-0439-x

    * Software dependencies
        https://github.com/nf-core/molkart/blob/master/CITATIONS.md

    executor >  local (22)
    [c1/da5009] NFCORE_MOLKART:MOLKART:MINDAGAP_MINDAGAP (mem_only)        [100%] 2 of 2 ✔
    [73/8f5e8a] NFCORE_MOLKART:MOLKART:CLAHE (mem_only)                    [100%] 2 of 2 ✔
    [ec/8f84d5] NFCORE_MOLKART:MOLKART:CREATE_STACK (mem_only)             [100%] 1 of 1 ✔
    [a2/99349b] NFCORE_MOLKART:MOLKART:MINDAGAP_DUPLICATEFINDER (mem_only) [100%] 1 of 1 ✔
    [95/c9b4b1] NFCORE_MOLKART:MOLKART:DEEPCELL_MESMER (mem_only)          [100%] 1 of 1 ✔
    [d4/1ebd1e] NFCORE_MOLKART:MOLKART:STARDIST (mem_only)                 [100%] 1 of 1 ✔
    [3e/3c0736] NFCORE_MOLKART:MOLKART:CELLPOSE (mem_only)                 [100%] 1 of 1 ✔
    [a0/415c6a] NFCORE_MOLKART:MOLKART:MASKFILTER (mem_only)               [100%] 3 of 3 ✔
    [14/a830c9] NFCORE_MOLKART:MOLKART:SPOT2CELL (mem_only)                [100%] 3 of 3 ✔
    [b5/391836] NFCORE_MOLKART:MOLKART:CREATE_ANNDATA (mem_only)           [100%] 3 of 3 ✔
    [77/aed558] NFCORE_MOLKART:MOLKART:MOLKARTQC (mem_only)                [100%] 3 of 3 ✔
    [e6/b81475] NFCORE_MOLKART:MOLKART:MULTIQC                             [100%] 1 of 1 ✔
    -[nf-core/molkart] Pipeline completed successfully-
    Completed at: 19-Oct-2025 22:23:01
    Duration    : 2m 52s
    CPU hours   : 0.1
    Succeeded   : 22
    ```

파이프라인이 따르는 nf-core 규칙 때문에 이 출력이 Hello World 예제보다 더 상세하다는 것을 알 수 있습니다:

- 파이프라인이 버전과 로고를 표시합니다
- 구성 매개변수가 표시됩니다
- 여러 프로세스가 병렬로 실행됩니다 (여러 프로세스 라인으로 표시)
- 프로세스 이름에는 전체 모듈 경로가 포함됩니다 (예: `NFCORE_MOLKART:MOLKART:MINDAGAP_MINDAGAP`)

### 2.4. 프로세스 실행 이해하기

executor 라인 `executor > local (22)`는 다음을 알려줍니다:

- **executor**: 사용 중인 컴퓨팅 환경 (`local` = 사용자의 머신)
- **(22)**: 실행된 총 작업 수

각 프로세스 라인은 다음을 보여줍니다:

- **해시** (`[1a/2b3c4d]`): 작업 디렉토리 식별자 (이전과 같음)
- **프로세스 이름**: 전체 모듈 경로 및 프로세스 이름
- **입력 식별자**: 괄호 안의 샘플 이름
- **진행률**: 완료 백분율 및 개수 (예: `1 of 1 ✔`)

### 핵심 사항

테스트 데이터로 nf-core 파이프라인을 실행하고 실행 출력을 해석하는 방법을 알게 되었습니다.

### 다음 단계

결과를 찾는 위치와 해석하는 방법을 배웁니다.

---

## 3. 출력물 찾기 및 검사하기

파이프라인이 성공적으로 완료되면 완료 메시지와 실행 요약이 표시됩니다.

### 3.1. 결과 디렉토리 찾기

기본적으로 nf-core 파이프라인은 `outdir` 매개변수로 지정된 디렉토리에 출력을 작성하며, 우리는 이를 `results/`로 설정했습니다.

내용을 나열하세요:

```bash
tree results/
```

여러 하위 디렉토리가 표시됩니다:

```console title="results/"
results/
├── anndata/
├── clahe/
├── mindagap/
├── molkartqc/
├── multiqc/
├── pipeline_info/
├── segmentation/
├── spot2cell/
└── stack/
```

각 하위 디렉토리에는 파이프라인의 특정 단계에서 나온 출력물이 포함되어 있습니다:

- **mindagap/**: MindaGap 전처리 단계의 그리드 채워진 이미지
- **clahe/**: CLAHE 전처리의 대비 향상된 이미지
- **stack/**: 분할을 위해 생성된 다중 채널 이미지 스택
- **segmentation/**: 다양한 알고리즘의 분할 결과 (cellpose/, mesmer/, stardist/, filtered_masks/)
- **spot2cell/**: 세포별 전사체 카운트 테이블
- **anndata/**: 세포별 전사체 매트릭스 및 공간 좌표를 포함하는 AnnData 객체
- **molkartqc/**: 스팟 할당을 위한 품질 관리 메트릭
- **multiqc/**: 포괄적인 품질 관리 보고서
- **pipeline_info/**: 실행 보고서 및 로그

### 3.2. MultiQC 보고서 검사하기

MultiQC 보고서는 모든 파이프라인 단계의 품질 메트릭을 집계하는 포괄적인 HTML 파일입니다.

파일 브라우저에서 보고서를 연 다음 "Show Preview" 버튼을 클릭하여 VS Code에서 직접 렌더링된 내용을 확인하세요.

보고서에는 다음이 포함됩니다:

- 모든 샘플에 대한 일반 통계
- 전처리 메트릭
- 분할 품질 메트릭
- 감지된 세포 및 스팟 수

!!! Tip

    MultiQC 보고서는 일반적으로 모든 nf-core 파이프라인에 포함됩니다.
    항상 파이프라인 실행 및 데이터 품질에 대한 높은 수준의 개요를 제공합니다.

### 3.3. 세포별 전사체 테이블 검사하기

가장 중요한 과학적 출력물은 세포별 전사체 카운트 테이블입니다.
이는 각 세포에서 각 전사체가 몇 개 감지되었는지 알려줍니다.

spot2cell 디렉토리로 이동하세요:

```bash
ls results/spot2cell/
```

다음과 같은 파일들을 찾을 수 있습니다:

- `cellxgene_mem_only_cellpose.csv`: Cellpose 분할을 사용한 세포별 전사체 테이블
- `cellxgene_mem_only_mesmer.csv`: Mesmer 분할을 사용한 세포별 전사체 테이블
- `cellxgene_mem_only_stardist.csv`: Stardist 분할을 사용한 세포별 전사체 테이블

이 테스트 데이터셋에서는 1개의 샘플만 실행했지만, 실제 실험에서는 각 샘플에 대해 이러한 테이블이 있을 것입니다.
Nextflow가 여러 분할 방법을 병렬로 처리하여 결과를 쉽게 비교할 수 있도록 하는 방법을 주목하세요.

### 3.4. 실행 보고서 보기

Nextflow는 여러 실행 보고서를 자동으로 생성합니다.

pipeline_info 디렉토리를 확인하세요:

```bash
ls results/pipeline_info/
```

주요 파일:

- **execution_report.html**: 타임라인 및 리소스 사용량 시각화
- **execution_timeline.html**: 프로세스 실행의 Gantt 차트
- **execution_trace.txt**: 상세한 작업 실행 메트릭
- **pipeline_dag.html**: 워크플로 구조를 보여주는 방향성 비순환 그래프

실행 보고서를 열어 리소스 사용량을 확인하세요:

```bash
code results/pipeline_info/execution_report.html
```

다음을 보여줍니다:

- 각 프로세스가 소요된 시간
- CPU 및 메모리 사용량
- 캐시된 작업과 실행된 작업

!!! Tip

    이러한 보고서는 리소스 할당을 최적화하고 성능 문제를 해결하는 데 매우 유용합니다.

### 핵심 사항

파이프라인 출력물을 찾고, 품질 관리 보고서를 검사하고, 실행 메트릭에 액세스하는 방법을 알게 되었습니다.

### 다음 단계

작업 디렉토리와 Nextflow가 중간 파일을 관리하는 방법에 대해 배웁니다.

---

## 4. 작업 디렉토리 탐색하기

Hello World 예제와 마찬가지로 모든 실제 작업은 `work/` 디렉토리에서 발생합니다.

### 4.1. 작업 디렉토리 구조 이해하기

작업 디렉토리에는 실행된 각 작업에 대한 하위 디렉토리가 포함됩니다.
22개의 작업이 있는 이 파이프라인의 경우 22개의 작업 하위 디렉토리가 있습니다.

작업 디렉토리를 나열하세요:

```bash
ls -d work/*/*/ | head -5
```

처음 5개의 작업 디렉토리가 표시됩니다.

### 4.2. 작업 디렉토리 검사하기

콘솔 출력의 분할 프로세스 해시 중 하나를 선택하고 (예: `[3m/4n5o6p]`) 내부를 살펴보세요:

```bash
ls -la work/3m/4n5o6p*/
```

다음을 볼 수 있습니다:

- **.command.\* 파일**: Nextflow 실행 스크립트 및 로그 (이전과 같음)
- **스테이징된 입력 파일**: 실제 입력 파일에 대한 심볼릭 링크
- **출력 파일**: 분할 마스크, 중간 결과 등

Hello World와의 주요 차이점:

- 실제 파이프라인은 대용량 입력 파일 (이미지, 참조 데이터)을 스테이징합니다
- 출력 파일이 상당히 클 수 있습니다 (분할 마스크, 처리된 이미지)
- 작업당 여러 입력 및 출력 파일

!!! Tip

    프로세스가 실패하면 해당 작업 디렉토리로 이동하여 `.command.err`에서 오류 메시지를 검사하고, 문제를 디버그하기 위해 `.command.sh`를 수동으로 다시 실행할 수도 있습니다.

### 4.3. 작업 디렉토리 정리

여러 파이프라인 실행을 거치면서 작업 디렉토리가 상당히 커질 수 있습니다.
Part 1에서 배운 것처럼 `nextflow clean`을 사용하여 이전 실행의 작업 디렉토리를 제거할 수 있습니다.

하지만 대용량 중간 파일이 있는 nf-core 파이프라인의 경우 정기적으로 정리하는 것이 특히 중요합니다.

### 핵심 사항

nf-core 파이프라인이 작업 디렉토리를 구성하는 방법과 디버깅을 위해 개별 작업을 검사하는 방법을 이해했습니다.

### 다음 단계

Nextflow 캐시와 실패한 파이프라인 실행을 재개하는 방법에 대해 배웁니다.

---

## 5. 파이프라인 실행 재개하기

Nextflow의 가장 강력한 기능 중 하나는 실패 지점부터 파이프라인을 재개할 수 있는 능력입니다.

### 5.1. 캐시 메커니즘

`-resume`으로 파이프라인을 실행하면 Nextflow는:

1. 각 작업에 대해 캐시를 확인합니다
2. 입력, 코드 및 매개변수가 동일하면 캐시된 결과를 재사용합니다
3. 변경되거나 실패한 작업만 다시 실행합니다

이는 실행 후반에 실패가 발생할 수 있는 장시간 실행 파이프라인에 필수적입니다.

### 5.2. molkart로 재개 시도하기

동일한 명령을 다시 실행하되 `-resume`을 추가하세요:

```bash
nextflow run ./molkart \
  --input 'data/samplesheet.csv' \
  --mindagap_tilesize 90 \
  --mindagap_boxsize 7 \
  --mindagap_loopnum 100 \
  --clahe_pyramid_tile 368 \
  --segmentation_method "cellpose" \
  --outdir results \
  -resume
```

다음과 같은 출력이 표시됩니다: <!-- TODO: full output -->

```console
executor >  local (0)
[1a/2b3c4d] NFCORE_MOLKART:MOLKART:MINDAGAP_MINDAGAP (mem_only)        [100%] 2 of 2, cached: 2 ✔
[5e/6f7g8h] NFCORE_MOLKART:MOLKART:CLAHE (mem_only)                    [100%] 2 of 2, cached: 2 ✔
[7f/8g9h0i] NFCORE_MOLKART:MOLKART:CREATE_STACK (mem_only)             [100%] 1 of 1, cached: 1 ✔
[9h/0i1j2k] NFCORE_MOLKART:MOLKART:MINDAGAP_DUPLICATEFINDER (mem_only) [100%] 1 of 1, cached: 1 ✔
[2k/3l4m5n] NFCORE_MOLKART:MOLKART:CELLPOSE (mem_only)                 [100%] 1 of 1, cached: 1 ✔
...
```

각 프로세스에 대해 `cached: 2` 또는 `cached: 1`이 표시되는 것을 주목하세요 - 아무것도 다시 실행되지 않았습니다!

### 5.3. 재개가 유용한 경우

재개는 다음과 같은 경우에 특히 유용합니다:

- 리소스 제한으로 인해 파이프라인이 실패한 경우 (메모리 부족, 시간 제한 초과)
- 업스트림 단계를 다시 실행하지 않고 다운스트림 프로세스를 수정해야 하는 경우
- 데이터 다운로드 중 네트워크 연결이 끊긴 경우
- 계산을 다시 하지 않고 추가 출력물을 추가하려는 경우

!!! Warning

    입력 데이터, 파이프라인 코드 또는 매개변수를 변경하지 않은 경우에만 재개가 작동합니다.
    이 중 하나를 변경하면 Nextflow는 영향을 받는 작업을 올바르게 다시 실행합니다.

### 핵심 사항

성공한 작업을 반복하지 않고 파이프라인을 효율적으로 다시 실행하기 위해 `-resume`을 사용하는 방법을 알게 되었습니다.

### 다음 단계

이제 테스트 데이터로 nf-core/molkart를 실행할 수 있으므로 자신의 데이터셋에 맞게 구성하는 방법을 배울 준비가 되었습니다.
