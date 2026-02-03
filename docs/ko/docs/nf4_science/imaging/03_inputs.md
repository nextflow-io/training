# 파트 3: 입력 구성하기

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } AI 지원 번역 - [자세히 알아보기 및 개선 제안](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

Part 2에서는 명령줄에서 여러 매개변수를 사용하여 molkart를 실행했습니다.
이제 입력을 관리하는 두 가지 더 나은 방법인 **매개변수 파일**과 **샘플시트**를 배우겠습니다.

## 1. 매개변수 파일 사용하기

### 1.1. 긴 명령줄의 문제점

Part 2의 명령을 다시 살펴보겠습니다:

```bash
nextflow run ./molkart \
  --input 'data/samplesheet.csv' \
  --mindagap_tilesize 90 \
  --mindagap_boxsize 7 \
  --mindagap_loopnum 100 \
  --clahe_pyramid_tile 368 \
  --segmentation_method "cellpose" \
  --outdir results
```

이 방법은 작동하지만, 재현하거나 공유하거나 수정하기 어렵습니다.
다음 달에 동일한 분석을 다시 실행해야 한다면 어떻게 하시겠습니까?
동료가 정확히 같은 설정을 사용하고 싶어한다면 어떻게 하시겠습니까?

### 1.2. 해결책: 매개변수 파일 사용

`params.yaml`이라는 파일을 생성하십시오:

```yaml title="params.yaml"
input: "data/samplesheet.csv"
outdir: "results"
mindagap_tilesize: 90
mindagap_boxsize: 7
mindagap_loopnum: 100
clahe_pyramid_tile: 368
segmentation_method: "cellpose"
```

이제 명령이 다음과 같이 간단해집니다:

```bash
nextflow run ./molkart -params-file params.yaml -resume
```

이것으로 끝입니다! 매개변수 파일은 정확한 구성을 문서화하고 재실행이나 공유를 쉽게 만듭니다.

### 1.3. 매개변수 재정의하기

명령줄에서 특정 매개변수를 여전히 재정의할 수 있습니다:

```bash
nextflow run ./molkart -params-file params.yaml --segmentation_method "stardist" --outdir stardist_results -resume
```

위 명령은 `segmentation_method`를 `stardist`로 변경하고 `--outdir` 이름을 `params.yaml` 파일의 매개변수 대신 `stardist_results`로 변경합니다.
또한 `-resume` 플래그가 이전 실행의 전처리 결과를 재사용하여 시간을 절약했음을 확인할 수 있습니다.
이 패턴을 사용하여 파이프라인의 다양한 변형을 빠르게 테스트할 수 있습니다.

### 요점

매개변수 파일은 분석을 재현 가능하고 공유하기 쉽게 만듭니다.
실제 분석 작업에는 매개변수 파일을 사용하십시오.

### 다음 단계

샘플시트가 여러 샘플에 대한 정보를 어떻게 구성하는지 알아보겠습니다.

---

## 2. 샘플시트 패턴

### 2.1. 샘플시트란 무엇입니까?

샘플시트는 입력 샘플을 설명하는 CSV 파일입니다.
각 행은 샘플이며, 열은 해당 샘플의 파일과 메타데이터를 지정합니다.

지금까지 사용해온 샘플시트를 살펴보겠습니다:

```bash
cat data/samplesheet.csv
```

```csv title="samplesheet.csv"
sample,nuclear_image,spot_table,membrane_image
mem_only,https://raw.githubusercontent.com/nf-core/test-datasets/molkart/test_data/input_data/nuclear.tiff,https://raw.githubusercontent.com/nf-core/test-datasets/molkart/test_data/input_data/spots.txt,https://raw.githubusercontent.com/nf-core/test-datasets/molkart/test_data/input_data/membrane.tiff
```

열은 다음과 같습니다:

- `sample`: 고유한 샘플 식별자
- `nuclear_image`: 핵 염색 이미지 (TIFF)
- `spot_table`: 전사체 스팟 (TXT)
- `membrane_image`: 막 염색 이미지 (TIFF, 선택사항)

### 2.2. 파일 경로

샘플시트는 여러 경로 유형을 허용합니다:

- **URL**: Nextflow가 자동으로 다운로드합니다 (위에 표시된 대로)
- **로컬 경로**: `data/nuclear.tiff` 또는 `/absolute/path/to/nuclear.tiff`
- **클라우드 스토리지**: `s3://bucket/nuclear.tiff`, `gs://bucket/nuclear.tiff`, `az://container/nuclear.tiff`

동일한 샘플시트에서 경로 유형을 혼합할 수 있습니다.

### 2.3. 자신만의 샘플시트 만들기

먼저 테스트 데이터 파일을 로컬로 다운로드하겠습니다:

```bash
cd /workspaces/training/nf4-science/imaging/data
curl -O https://raw.githubusercontent.com/nf-core/test-datasets/molkart/test_data/input_data/nuclear.tiff
curl -O https://raw.githubusercontent.com/nf-core/test-datasets/molkart/test_data/input_data/spots.txt
curl -O https://raw.githubusercontent.com/nf-core/test-datasets/molkart/test_data/input_data/membrane.tiff
cd ..
```

이제 이러한 로컬 파일을 참조하도록 샘플시트를 수정하겠습니다:

```csv title="samplesheet.csv"
sample,nuclear_image,spot_table,membrane_image
mem_only,data/nuclear.tiff,data/spots.txt,data/membrane.tiff
```

!!! warning "경고"

    샘플시트의 경로는 샘플시트가 위치한 곳이 아니라 Nextflow를 **실행**하는 위치에 상대적입니다.

마지막으로 로컬 파일 경로가 있는 샘플시트로 nf-core/molkart를 한 번 더 실행하겠습니다:

`nextflow run ./molkart -params-file params.yaml -resume`

보시다시피 Nextflow는 파일이 Github에서 다운로드되었을 때와 유사하게 이 실행을 수행합니다. 이것은 Nextflow의 훌륭한 기능 중 하나로, 데이터가 어디에 위치하든 상관없이 적절하게 스테이징합니다.

### 요점

샘플시트는 파일 경로와 함께 메타데이터를 명시적으로 정의할 수 있는 방식으로 다중 샘플 데이터셋을 구성합니다.
대부분의 nf-core 파이프라인은 이 패턴을 사용합니다.

### 다음 단계

입력에 대해 다루었으니, 이제 다양한 컴퓨팅 환경에 맞게 Nextflow 파이프라인을 구성하는 방법을 살펴보겠습니다.
