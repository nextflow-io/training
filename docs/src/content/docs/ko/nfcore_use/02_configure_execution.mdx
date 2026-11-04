# 파트 2: 파이프라인 실행 설정

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } AI 지원 번역 - [자세히 알아보기 및 개선 제안](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

[파트 1](./01_run_demo.md)에서는 nf-core/demo 파이프라인을 찾아 test 프로파일로 실행했습니다.
이제 파이프라인 실행 설정 방법을 살펴봅니다. 매개변수 설정, 유효성 검사 이해, 리소스 할당 및 도구 인자 맞춤화를 다룹니다.

[Hello Config](../hello_nextflow/06_hello_config.md)에서 설명한 것처럼, 파이프라인 코드 자체를 변경하지 않고도 파이프라인이 실행할 데이터와 실행 방식을 변경할 수 있어야 합니다.
이를 위해 Nextflow는 파이프라인 설정을 제어하는 여러 방법을 지원하는데, 처음에는 다소 복잡하게 느껴질 수 있습니다.

nf-core 프로젝트는 설정 요소를 구성하는 규칙을 정의하며, 최상위 수준에서 두 가지 종류의 설정을 구분합니다. 바로 **파이프라인 매개변수**와 엄밀한 의미의 **설정(configuration)**입니다.

- **파이프라인 매개변수** (`params` 시스템을 통해 설정)는 일반적으로 입력 파일, 도구 동작 플래그, 분석 매개변수 등을 포함합니다.
- 엄밀한 의미의 **설정**은 파이프라인이 실행되는 방식, 즉 executor, 컴퓨팅 리소스 할당 등의 운영 측면을 의미합니다.

<figure class="excalidraw">
    --8<-- "docs/en/docs/nfcore_use/img/params_vs_config.excalidraw.svg"
</figure>

먼저 파이프라인 매개변수를 살펴본 후, 엄밀한 의미의 설정을 확인합니다.

---

## 1. 파이프라인 매개변수

모든 nf-core 파이프라인에서는 `--help` 플래그를 사용하여 명령줄에서 직접 파이프라인 매개변수의 전체 목록을 확인할 수 있습니다. `--help` 자체도 파이프라인 매개변수입니다.

### 1.1. `--help`로 매개변수 목록 확인하기

demo 파이프라인의 도움말 명령을 실행합니다.

```bash
nextflow run nf-core/demo --help
```

??? success "명령 출력"

    ```console
     N E X T F L O W   ~  version 26.04.4

    Launching `https://github.com/nf-core/demo` [adoring_meucci] revision: 32893afef8 [master]


    ------------------------------------------------------
                                            ,--./,-.
            ___     __   __   __   ___     /,-._.--~'
      |\ | |__  __ /  ` /  \ |__) |__         }  {
      | \| |       \__, \__/ |  \ |___     \`-._,-`-,
                                            `._,._,'
      nf-core/demo 1.2.0
    ------------------------------------------------------
    Typical pipeline command:

      nextflow run nf-core/demo -profile <docker/singularity/.../institute> --input samplesheet.csv --outdir <OUTDIR>


    Input/output options
      --input                       [string] Path to a metadata file containing information about the samples in the experiment.
      --outdir                      [string] The output directory where the results will be saved. You have to use absolute paths to storage on Cloud infrastructure.
      --email                       [string] Email address for completion summary.
      --multiqc_title               [string] MultiQC report title. Printed as page header, used for filename if not otherwise specified.

    Reference genome options
      --genome                      [string] Name of iGenomes reference.
      --fasta                       [string] Path to FASTA genome file.

    Process skipping options
      --skip_trim                   [boolean] Skip trimming fastq files with seqtk

    Generic options
      --multiqc_methods_description [string]          Custom MultiQC yaml file containing HTML including a methods description.
      --help                        [boolean, string] Display the help message.
      --help_full                   [boolean]         Display the full detailed help message.
      --show_hidden                 [boolean]         Display hidden parameters in the help message (only works when --help or --help_full are provided).
     !! Hiding 19 param(s), use the `--showHidden` parameter to show them !!
    ------------------------------------------------------

    * The pipeline
        https://doi.org/10.5281/zenodo.12192442

    * The nf-core framework
        https://doi.org/10.1038/s41587-020-0439-x

    * Software dependencies
        https://github.com/nf-core/demo/blob/master/CITATIONS.md
    ```

출력 결과에서 매개변수가 카테고리별로 그룹화되어 있으며(입력/출력 옵션, 레퍼런스 게놈 옵션 등), 각 매개변수의 타입과 설명이 함께 표시됩니다.

이 카테고리 분류는 스키마 파일에 의해 결정되며, 자세한 내용은 아래에서 다룹니다.
일반 Nextflow 파이프라인에서는 개발자가 직접 구현한 경우에만 `--help`가 동작합니다.

!!! tip "팁"

    `--help --show_hidden`을 사용하면 기본적으로 숨겨진 추가 매개변수(예: `--publish_dir_mode`, `--monochrome_logs`)를 확인할 수 있습니다.

### 1.2. 매개변수 값 설정하기

[Hello Config](../hello_nextflow/06_hello_config.md)에서 다룬 것처럼, 명령줄에서 `--param_name`으로 매개변수 값을 설정하거나, YAML 파일에 매개변수 집합을 모아 `-params-file`로 전달할 수 있습니다.
두 방법 모두 nf-core 파이프라인에서 동일하게 동작합니다.

예를 들어, 트리밍 단계를 건너뛰려면 boolean 매개변수 `skip_trim`을 `true`로 설정해야 합니다.
작업 디렉토리에 해당 값이 이미 설정된 `my_params.yml` 파일이 제공되어 있습니다.

```yaml title="my_params.yml"
skip_trim: true
```

`-params-file`로 전달합니다.

```bash
nextflow run nf-core/demo -profile test,docker --outdir demo-results-notrim -params-file my_params.yml
```

??? success "명령 출력"

    ```console
     N E X T F L O W   ~  version 26.04.4

    Launching `https://github.com/nf-core/demo` [focused_heisenberg] revision: 32893afef8 [master]


    ------------------------------------------------------
                                            ,--./,-.
            ___     __   __   __   ___     /,-._.--~'
      |\ | |__  __ /  ` /  \ |__) |__         }  {
      | \| |       \__, \__/ |  \ |___     \`-._,-`-,
                                            `._,._,'
      nf-core/demo 1.2.0
    ------------------------------------------------------

    Input/output options
      input                     : https://raw.githubusercontent.com/nf-core/test-datasets/viralrecon/samplesheet/samplesheet_test_illumina_amplicon.csv
      outdir                    : demo-results-notrim

    Process skipping options
      skip_trim                 : true

    Institutional config options
      config_profile_name       : Test profile
      config_profile_description: Minimal test dataset to check pipeline function

    Generic options
      trace_report_suffix       : 2026-07-03_22-08-47

    Core Nextflow options
      revision                  : master
      runName                   : focused_heisenberg
      containerEngine           : docker
      launchDir                 : /workspaces/training/nfcore-use
      workDir                   : /workspaces/training/nfcore-use/work
      projectDir                : /workspaces/.nextflow/assets/.repos/nf-core/demo/clones/32893afef8076a03a2767a020b3f0cab2e0b40b2
      userName                  : root
      profile                   : test,docker
      configFiles               : /workspaces/.nextflow/assets/.repos/nf-core/demo/clones/32893afef8076a03a2767a020b3f0cab2e0b40b2/nextflow.config

    !! Only displaying parameters that differ from the pipeline defaults !!
    ------------------------------------------------------

    * The pipeline
        https://doi.org/10.5281/zenodo.12192442

    * The nf-core framework
        https://doi.org/10.1038/s41587-020-0439-x

    * Software dependencies
        https://github.com/nf-core/demo/blob/master/CITATIONS.md

    executor >  local (5)
    [7a/f3599e] NFCORE_DEMO:DEMO:FASTQC (SAMPLE3_SE) | 3 of 3 ✔
    [b0/2f0bdc] NFCORE_DEMO:DEMO:COWPY               | 1 of 1 ✔
    [c3/3c2278] NFCORE_DEMO:DEMO:MULTIQC (demo)      | 1 of 1 ✔
    -[nf-core/demo] Pipeline completed successfully-
    ```

출력 결과에서 `SEQTK_TRIM` 프로세스가 더 이상 나타나지 않습니다.

!!! warning "매개변수 입력에 관한 중요한 제한 사항"

    **명령줄에서 boolean 매개변수 설정하기**

    Nextflow 버전 26.04부터 명령줄에서 전달되는 모든 값은 string으로 처리됩니다.
    `skip_trim`과 같은 boolean 매개변수를 단독 플래그(`--skip_trim`)나 `--skip_trim true`로 전달하면 **string** `"true"`로 평가되어 스키마 유효성 검사에 실패합니다.

    ```console
    * --skip_trim (true): Value is [string] but should be [boolean]
    ```

    boolean 매개변수를 실제 `true`/`false` 값으로 설정하려면 위에서 보여준 것처럼 `-params-file`을 사용하거나 설정 파일에서 설정하세요.
    string, integer, 파일 경로 매개변수는 영향을 받지 않으며 명령줄에서 직접 설정할 수 있습니다.
    이 과정에서는 boolean 매개변수에 대해 이 방식을 일관되게 사용합니다.

    **사용자 정의 설정 파일 사용하기**

    `-c`로 전달하는 사용자 정의 설정 파일에서 파이프라인 매개변수를 설정하는 것이 기술적으로는 가능하지만, Nextflow의 설정 우선순위 규칙에 따라 파이프라인 자체의 `nextflow.config`에 이미 설정된 기본값을 재정의하지 못할 수 있습니다.
    명령줄에서 `--param_name`을 사용하거나 `-params-file`을 사용하는 것이 더 안정적입니다. 이 방법들은 항상 우선순위가 높습니다.

    경험상 원칙: `--help` 출력에 나타나는 매개변수는 설정 파일이 아닌 명령줄이나 params 파일을 통해 설정하세요.

### 1.3. 매개변수 유효성 검사

흥미로운 사실: nf-core 프로젝트는 개발자가 모든 파이프라인 매개변수를 JSON 스키마 파일(`nextflow_schema.json`)에 공식적으로 정의하도록 요구하기 때문에, `--help` 명령이 모든 nf-core 파이프라인에서 동작합니다.
이 스키마는 각 매개변수의 타입, 설명, 기본값, 그룹화 정보를 기록합니다.

`--help` 출력을 지원하는 것 외에도, 스키마 파일은 실행 시 자동 유효성 검사를 가능하게 합니다.
즉, Nextflow는 전달된 모든 매개변수가 존재하는지, 적절한 값(적절한 타입, 허용된 값 범위 내)이 지정되었는지 확인할 수 있습니다.

이에 대한 자세한 내용은 [입력 유효성 검사 섹션](../nfcore_build/04_input_validation.md)에서 다루지만, demo 파이프라인에 잘못된 매개변수 입력을 제공하여 이미 동작을 확인할 수 있습니다.

#### 1.3.1. 인식되지 않는 매개변수

존재하지 않는 매개변수를 전달해 봅니다.

```bash
nextflow run nf-core/demo -profile test,docker --outdir demo-results --foobar "invalid"
```

콘솔 출력에 경고가 포함됩니다.

```console
WARN: The following invalid input values have been detected:

* --foobar: invalid
```

파이프라인은 계속 실행되지만, `--foobar`가 인식되지 않는 매개변수임을 즉시 경고합니다.
이는 `--outDir`을 `--outdir` 대신 사용하는 것과 같이 실행에 영향을 주지 않는 오타를 발견하는 데 도움을 주어, 시간과 컴퓨팅 자원을 낭비하지 않도록 합니다.

#### 1.3.2. 잘못된 매개변수 값

유효성 검사는 매개변수 **값**도 확인합니다.
`--skip_trim` 매개변수는 boolean 플래그이므로, string 값을 전달하면 파이프라인이 즉시 실패합니다.

```bash
nextflow run nf-core/demo -profile test,docker --outdir demo-results --skip_trim yes
```

```console
ERROR ~ Validation of pipeline parameters failed!

 -- Check '.nextflow.log' file for details
The following invalid input values have been detected:

* --skip_trim (yes): Value is [string] but should be [boolean]

 -- Check script '/workspaces/.nextflow/assets/.repos/nf-core/demo/clones/32893afef8076a03a2767a020b3f0cab2e0b40b2/subworkflows/nf-core/utils_nfschema_plugin/main.nf' at line: 68 or see '.nextflow.log' file for more details
```

어떤 프로세스도 실행되기 전에 파이프라인이 중단되어, 실패하거나 잘못된 실행을 방지합니다.
[1.2](#12-set-parameter-values)에서 언급한 것처럼, boolean 매개변수는 명령줄 값이 string으로 처리되므로 명령줄에서 전달하는 대신 params 파일에서 실제 `true`/`false` 값으로 설정해야 합니다.

### 1.4. 입력 유효성 검사

동일한 유효성 검사 로직을 입력 파일의 유효성 확인에도 사용할 수 있습니다.
예를 들어, 파이프라인이 주요 데이터 입력으로 샘플시트를 기대하는 경우(많은 nf-core 파이프라인이 이에 해당), 개발자는 입력 파일의 구조를 설명하는 입력 스키마(매개변수 스키마와는 별개)를 제공할 수 있습니다.

그러면 런타임에 Nextflow가 제공된 입력 파일이 유효한지 확인합니다.

이에 대한 자세한 내용도 [입력 유효성 검사 섹션](../nfcore_build/04_input_validation.md)에서 다루지만, demo 파이프라인에 잘못된 입력 샘플시트를 제공하여 이미 동작을 확인할 수 있습니다.

`nf-core/demo` 파이프라인은 `sample`, `fastq_1`, `fastq_2` 열을 가진 CSV 파일을 기대합니다.
이는 예상 구조, 열 타입, 제약 조건을 지정하는 스키마 파일(`assets/schema_input.json`)에 정의되어 있습니다.

??? abstract "입력에 대한 스키마 파일"

    ```json title="assets/schema_input.json"
    {
        "$schema": "https://json-schema.org/draft/2020-12/schema",
        "$id": "https://raw.githubusercontent.com/nf-core/demo/master/assets/schema_input.json",
        "title": "nf-core/demo pipeline - params.input schema",
        "description": "Schema for the file provided with params.input",
        "type": "array",
        "items": {
            "type": "object",
            "properties": {
                "sample": {
                    "type": "string",
                    "pattern": "^\\S+$",
                    "errorMessage": "Sample name must be provided and cannot contain spaces",
                    "meta": ["id"]
                },
                "fastq_1": {
                    "type": "string",
                    "format": "file-path",
                    "exists": true,
                    "pattern": "^([\\S\\s]*\\/)?[^\\s\\/]+\\.f(ast)?q\\.gz$",
                    "errorMessage": "FastQ file for reads 1 must be provided, cannot contain spaces and must have extension '.fq.gz' or '.fastq.gz'"
                },
                "fastq_2": {
                    "type": "string",
                    "format": "file-path",
                    "exists": true,
                    "pattern": "^([\\S\\s]*\\/)?[^\\s\\/]+\\.f(ast)?q\\.gz$",
                    "errorMessage": "FastQ file for reads 2 cannot contain spaces and must have extension '.fq.gz' or '.fastq.gz'"
                }
            },
            "required": ["sample", "fastq_1"]
        }
    }
    ```

스키마는 `sample`과 `fastq_1`이 필수이며, `fastq_2`는 선택 사항(페어드 엔드 및 단일 엔드 데이터 모두 지원)임을 지정합니다.
파일 경로는 존재 여부와 확장자 패턴에 대해 유효성이 검사됩니다.

이를 확인하기 위해 작업 디렉토리에 `malformed_samplesheet.csv`라는 잘못된 형식의 샘플시트가 제공되어 있습니다.

```csv title="malformed_samplesheet.csv"
sample,fastq_2
SAMPLE1,/not/a/real/file.fastq.gz
```

이 샘플시트는 필수 `fastq_1` 열이 없고, `fastq_2`에 존재하지 않는 파일 경로가 있습니다.

`malformed_samplesheet.csv`를 입력으로 사용하여 demo 파이프라인을 실행합니다.

```bash
nextflow run nf-core/demo -profile test,docker --outdir demo-results --input malformed_samplesheet.csv
```

```console
ERROR ~ Validation of pipeline parameters failed!

 -- Check '.nextflow.log' file for details
The following invalid input values have been detected:

* --input (malformed_samplesheet.csv): Validation of file failed:
    -> Entry 1: Error for field 'fastq_2' (/not/a/real/file.fastq.gz): the file or directory '/not/a/real/file.fastq.gz' does not exist (FastQ file for reads 2 cannot contain spaces and must have extension '.fq.gz' or '.fastq.gz')
    -> Entry 1: Missing required field(s): fastq_1

 -- Check script '/workspaces/.nextflow/assets/.repos/nf-core/demo/clones/32893afef8076a03a2767a020b3f0cab2e0b40b2/subworkflows/nf-core/utils_nfschema_plugin/main.nf' at line: 68 or see '.nextflow.log' file for more details
```

파이프라인이 즉시 실패하고 **모든** 유효성 검사 오류를 한 번에 보고합니다.
nf-schema는 첫 번째 오류에서 멈추지 않고 모든 문제를 수집하여 함께 나열하므로, 문제를 하나씩 발견하는 대신 한 번에 모두 수정할 수 있습니다.

각 오류는 문제를 일으킨 정확한 항목과 필드를 식별하므로, 샘플시트를 수정한 후 Nextflow가 실제로 파일 경로에 접근할 때 나중에 실패하지 않을 것이라는 확신을 가지고 파이프라인을 다시 실행할 수 있습니다.

개발자를 위한 자세한 내용은 [Build with nf-core 파트 4](../nfcore_build/04_input_validation.md)에서 다룹니다.

### 핵심 정리

`--help`로 파이프라인의 전체 매개변수 목록을 확인하고, 명령줄이나 params 파일을 통해 설정하는 방법을 학습했습니다. 또한 Nextflow가 파이프라인의 스키마에 대해 매개변수 값과 입력 파일 모두의 유효성을 검사하는 방법도 확인했습니다.

### 다음 단계

다른 종류의 설정인 파이프라인 실행 방식, 즉 리소스 할당과 도구 인자에 대해 학습합니다.

---

## 2. 설정(Configuration)

엄밀한 의미의 설정은 파이프라인이 **어떻게** 실행되는지를 제어합니다. 리소스 할당, 도구별 인자, 작업이 실행되는 위치, 사용할 소프트웨어 패키징 시스템 등이 이에 해당합니다.

nf-core 파이프라인은 `nextflow.config`와 `conf/` 디렉토리에 기본 설정을 포함합니다.
무언가를 재정의하기 전에 기본값이 어디에 있는지 파악하는 것이 도움이 됩니다.

### 2.1. 설정 파일 살펴보기

[파트 1](./01_run_demo.md)에서 파이프라인 소스 코드가 `$NXF_HOME/assets` 아래에 있다는 것을 확인했습니다.
[파트 1](./01_run_demo.md)에서 생성한 `pipelines` 심볼릭 링크를 사용하여 설정 파일 목록을 확인합니다.

```bash
ls pipelines/nf-core/demo/conf/
```

```console
base.config
containers_conda_lock_files_amd64.config
containers_conda_lock_files_arm64.config
containers_docker_amd64.config
containers_docker_arm64.config
containers_singularity_https_amd64.config
containers_singularity_https_arm64.config
containers_singularity_oras_amd64.config
containers_singularity_oras_arm64.config
igenomes.config
igenomes_ignored.config
modules.config
test.config
test_full.config
```

<figure class="excalidraw">
--8<-- "docs/en/docs/nfcore_use/img/nfcore_config_files.excalidraw.svg"
</figure>

가장 중요한 설정 파일은 다음과 같습니다.

- **`conf/base.config`**: 프로세스에 CPU, 메모리, 시간을 할당하는 리소스 레이블(`process_low`, `process_medium`, `process_high`)을 정의합니다. 프로세스가 예상보다 많은 리소스를 사용하는 경우, 이 파일에서 해당 기본값을 확인할 수 있습니다.
- **`conf/modules.config`**: 프로세스별 도구 인자(`ext.args`)와 출력 게시 설정(`publishDir`)을 지정합니다. 이 파일을 열어 각 도구가 기본적으로 받는 인자를 확인하세요.
- **`conf/test.config`**: [파트 1](./01_run_demo.md)에서 사용한 test 프로파일로, `resourceLimits`를 통해 리소스를 제한하고 테스트 샘플시트를 설정합니다. `-profile test`로 활성화됩니다.
  전체 크기 테스트 데이터셋으로 실행하기 위한 `conf/test_full.config`도 있으며, 벤치마킹에 유용합니다.

중앙의 `nextflow.config`는 위의 모든 파일을 로드하고 모든 항목에 대한 적절한 기본값을 설정합니다.

이 파일들에 지정된 설정을 수정하려면 해당 파일을 직접 수정하지 마세요.
대신 자체 설정 파일을 만들어 `-c`로 전달하세요.
지정한 값이 다른 파일에 설정된 기본값을 재정의합니다.

실제로 적용해 봅니다.

### 2.2. 프로세스 리소스 및 도구 인자 맞춤화하기

nf-core 모듈은 두 가지 일반적인 설정 재정의 유형을 지원합니다. **리소스 할당**(CPU, 메모리, 시간)과 `ext.args`를 통한 **도구 인자**입니다.

많은 명령줄 도구에는 파이프라인 매개변수로 노출하기에는 사용 빈도가 낮은 인자들이 있습니다.
`ext.args` 규칙을 사용하면 이러한 인자를 설정 파일을 통해 기본 도구에 전달할 수 있습니다.

작업 디렉토리에 제공된 `custom.config` 파일은 두 가지 재정의를 모두 보여줍니다.

```groovy title="custom.config" linenums="1"
process {
    withName: 'FASTQC' {
        cpus = 2
        memory = 4.GB
    }
    withName: 'SEQTK_TRIM' {
        ext.args = '-b 5'
    }
}
```

첫 번째 블록은 `FASTQC` 리소스 할당을 재정의합니다.
기본적으로 `FASTQC`는 `base.config`의 `process_medium` 레이블을 사용하여 6개의 CPU와 36GB 메모리를 할당합니다. 여기서는 2개의 CPU와 4GB로 제한합니다.

두 번째 블록은 `ext.args`를 통해 `SEQTK_TRIM`에 추가 인자를 전달합니다.
`-b 5` 플래그는 `seqtk trimfq`에게 품질 트리밍 외에도 각 리드의 시작 부분에서 5개의 염기를 트리밍하도록 지시합니다.

이 설정으로 파이프라인을 실행합니다.

```bash
nextflow run nf-core/demo -profile test,docker --outdir demo-results-custom -c custom.config
```

??? success "명령 출력"

    ```console
    executor >  local (8)
    [95/b32876] NFCORE_DEMO:DEMO:FASTQC (SAMPLE1_PE)     | 3 of 3 ✔
    [17/428668] NFCORE_DEMO:DEMO:SEQTK_TRIM (SAMPLE1_PE) | 3 of 3 ✔
    [cf/85991a] NFCORE_DEMO:DEMO:COWPY                   | 1 of 1 ✔
    [3c/94a7a0] NFCORE_DEMO:DEMO:MULTIQC (demo)          | 1 of 1 ✔
    -[nf-core/demo] Pipeline completed successfully-
    ```

`-c` 플래그는 파이프라인의 내장 설정 위에 사용자 설정을 추가합니다.

`ext.args` 재정의가 적용되었는지 확인하려면, 실행 출력에서 `SEQTK_TRIM` work 디렉토리 해시(예: `work/17/428668...`)를 찾아 내부의 `.command.sh` 파일을 확인합니다.

```bash
cat work/17/428668/.command.sh
```

??? success "명령 출력"

    ```console
    #!/usr/bin/env bash -e -u -o pipefail
    printf "%s\n" sample1_R1.fastq.gz sample1_R2.fastq.gz | while read f;
    do
        seqtk \
            trimfq \
            -b 5 \
            $f \
            | gzip --no-name > SAMPLE1_PE_$(basename $f)
    done
    ...
    ```

`seqtk trimfq` 명령에서 `-b 5`를 확인할 수 있습니다.

`ext.args`에 대해 알아야 할 중요한 사항이 있습니다. 모듈에 이미 기본값이 설정되어 있는 경우, 사용자의 값이 기존 값에 추가되는 것이 아니라 **완전히 대체**됩니다.
예를 들어, `FASTQC`는 `conf/modules.config`에 기본적으로 `ext.args = '--quiet'`가 설정되어 있습니다.

```groovy title="conf/modules.config" linenums="21" hl_lines="2"
    withName: FASTQC {
        ext.args   = '--quiet'
        publishDir = [
            path: { "${params.outdir}/fastqc/${meta.id}" },
            mode: params.publish_dir_mode,
            pattern: "*.{html,json}",
        ]
    }
```

`FASTQC`에 `ext.args = '--kmers 8'`을 설정하면 `--quiet` 플래그가 더 이상 적용되지 않습니다.
두 플래그를 모두 유지하려면 `ext.args = '--quiet --kmers 8'`로 설정하세요.

`ext.args`를 재정의하기 전에 항상 모듈의 기본 설정을 확인해야 합니다.

### 핵심 정리

nf-core 파이프라인 설정 기본값이 어디에 있는지 파악하고, 사용자 정의 설정 파일로 리소스 할당과 도구 인자를 재정의하는 방법을 학습했습니다.

### 다음 단계

[파트 3](./03_run_production_pipeline.md)으로 이동하여 지금까지 학습한 내용을 실제 프로덕션 파이프라인에 적용합니다.

---

## 요약

이 파트에서 학습한 내용:

- 도움말 확인, 매개변수 설정, 매개변수 및 입력 유효성 검사 이해
- 설정 파일을 통한 리소스 할당 및 도구 인자 맞춤화
