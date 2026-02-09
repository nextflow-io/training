# 파일 입력 처리

과학 분석 워크플로우는 종종 대량의 파일을 처리하는 작업을 포함합니다.
Nextflow는 파일을 효율적으로 처리할 수 있는 강력한 도구를 제공하여 최소한의 코드로 데이터를 구성하고 처리할 수 있도록 도와줍니다.

### 학습 목표

이 사이드 퀘스트에서는 기본 파일 작업부터 파일 컬렉션을 다루는 고급 기법까지 Nextflow가 파일을 처리하는 방법을 학습합니다.
과학 분석 파이프라인에서 일반적으로 요구되는 파일명에서 메타데이터를 추출하는 방법을 배우게 됩니다.

이 사이드 퀘스트를 마치면 다음을 수행할 수 있습니다:

- Nextflow의 `file()` 메서드를 사용하여 파일 경로 문자열에서 Path 객체 생성하기
- 이름, 확장자, 상위 디렉토리와 같은 파일 속성에 접근하기
- URI를 사용하여 로컬 및 원격 파일을 투명하게 처리하기
- `channel.fromPath()` 및 `channel.fromFilePairs()`를 사용하여 파일 처리 자동화하기
- 문자열 조작을 사용하여 파일명에서 메타데이터 추출 및 구조화하기
- 패턴 매칭 및 glob 표현식을 사용하여 관련 파일 그룹화하기
- 적절한 입력 처리를 통해 파일 작업을 Nextflow 프로세스에 통합하기
- 메타데이터 기반 디렉토리 구조를 사용하여 프로세스 출력 구성하기

이러한 기술은 다양한 종류의 파일 입력을 뛰어난 유연성으로 처리할 수 있는 워크플로우를 구축하는 데 도움이 됩니다.

### 사전 요구 사항

이 사이드 퀘스트를 시작하기 전에 다음을 준비해야 합니다:

- [Hello Nextflow](../../hello_nextflow/) 튜토리얼 또는 이에 상응하는 초급 과정을 완료했어야 합니다.
- 기본 Nextflow 개념 및 메커니즘(프로세스, 채널, 연산자)을 편안하게 사용할 수 있어야 합니다.

---

## 0. 시작하기

#### 교육 코드스페이스 열기

아직 열지 않았다면 [환경 설정](../envsetup/index.md)에 설명된 대로 교육 환경을 열어야 합니다.

[![Open in GitHub Codespaces](https://github.com/codespaces/badge.svg)](https://codespaces.new/nextflow-io/training?quickstart=1&ref=master)

#### 프로젝트 디렉토리로 이동

이 튜토리얼의 파일이 있는 디렉토리로 이동하겠습니다.

```bash
cd side-quests/working_with_files
```

VSCode가 이 디렉토리에 집중하도록 설정할 수 있습니다:

```bash
code .
```

#### 자료 검토

`main.nf`라는 간단한 워크플로우 파일, 두 개의 모듈 파일이 포함된 `modules` 디렉토리, 그리고 일부 예제 데이터 파일이 포함된 `data` 디렉토리를 찾을 수 있습니다.

??? abstract "디렉토리 내용"

    ```console
    .
    ├── data
    │   ├── patientA_rep1_normal_R1_001.fastq.gz
    │   ├── patientA_rep1_normal_R2_001.fastq.gz
    │   ├── patientA_rep1_tumor_R1_001.fastq.gz
    │   ├── patientA_rep1_tumor_R2_001.fastq.gz
    │   ├── patientA_rep2_normal_R1_001.fastq.gz
    │   ├── patientA_rep2_normal_R2_001.fastq.gz
    │   ├── patientA_rep2_tumor_R1_001.fastq.gz
    │   ├── patientA_rep2_tumor_R2_001.fastq.gz
    │   ├── patientB_rep1_normal_R1_001.fastq.gz
    │   ├── patientB_rep1_normal_R2_001.fastq.gz
    │   ├── patientB_rep1_tumor_R1_001.fastq.gz
    │   ├── patientB_rep1_tumor_R2_001.fastq.gz
    │   ├── patientC_rep1_normal_R1_001.fastq.gz
    │   ├── patientC_rep1_normal_R2_001.fastq.gz
    │   ├── patientC_rep1_tumor_R1_001.fastq.gz
    │   └── patientC_rep1_tumor_R2_001.fastq.gz
    ├── main.nf
    └── modules
        ├── analyze_reads.nf
        └── count_lines.nf
    ```

이 디렉토리에는 세 명의 환자(A, B, C)로부터 얻은 페어드 엔드 시퀀싱 데이터가 포함되어 있습니다.

각 환자에 대해 `tumor`(일반적으로 종양 생검에서 유래) 또는 `normal`(건강한 조직이나 혈액에서 채취) 유형의 샘플이 있습니다.
암 분석에 익숙하지 않다면, 이것이 대조 분석을 수행하기 위해 페어드 종양/정상 샘플을 사용하는 실험 모델에 해당한다는 것만 알아두세요.

특히 환자 A의 경우 두 세트의 기술적 반복(replicate)이 있습니다.

시퀀싱 데이터 파일은 '정방향 리드'와 '역방향 리드'로 알려진 것에 대한 일반적인 `_R1_` 및 `_R2_` 규칙으로 명명됩니다.

_이 실험 설계에 익숙하지 않더라도 걱정하지 마세요. 이 튜토리얼을 이해하는 데 중요하지 않습니다._

#### 과제 검토

여러분의 과제는 다음을 수행하는 Nextflow 워크플로우를 작성하는 것입니다:

1. Nextflow의 파일 처리 메서드를 사용하여 입력 파일 **로드**
2. 파일명 구조에서 메타데이터(환자 ID, 반복, 샘플 유형) **추출**
3. `channel.fromFilePairs()`를 사용하여 페어드 파일(R1/R2) **그룹화**
4. 제공된 분석 모듈로 파일 **처리**
5. 추출된 메타데이터를 기반으로 디렉토리 구조로 출력 **구성**

#### 준비 체크리스트

시작할 준비가 되었다고 생각하시나요?

- [ ] 이 과정의 목표와 사전 요구 사항을 이해했습니다
- [ ] 코드스페이스가 실행 중입니다
- [ ] 작업 디렉토리를 적절하게 설정했습니다
- [ ] 과제를 이해했습니다

모든 항목을 체크할 수 있다면 시작할 준비가 된 것입니다.

---

## 1. 기본 파일 작업

### 1.1. `.class`로 객체 유형 식별하기

워크플로우 파일 `main.nf`를 살펴보세요:

```groovy title="main.nf" linenums="1"
#!/usr/bin/env nextflow

workflow {

    // Create a Path object from a string path
    myFile = 'data/patientA_rep1_normal_R1_001.fastq.gz'

    println "${myFile} is of class ${myFile.class}"
}
```

이것은 워크플로우에서 단일 파일 경로를 참조한 다음 해당 클래스와 함께 콘솔에 출력하는 미니 워크플로우(프로세스 없음)입니다.

??? info "`.class`란 무엇인가요?"

    Nextflow에서 `.class`는 우리가 작업하고 있는 객체의 유형을 알려줍니다. 이것은 "이것이 어떤 종류의 것인가?"라고 물어보는 것과 같으며, 문자열인지, 숫자인지, 파일인지, 아니면 다른 것인지 알아내는 데 도움이 됩니다.
    이것은 다음 섹션에서 일반 문자열과 Path 객체의 차이를 설명하는 데 도움이 될 것입니다.

워크플로우를 실행해 봅시다:

```bash
nextflow run main.nf
```

??? success "명령 출력"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `main.nf` [romantic_chandrasekhar] DSL2 - revision: 5a4a89bc3a

    data/patientA_rep1_normal_R1_001.fastq.gz is of class java.lang.String
    ```

보시다시피 Nextflow는 우리가 작성한 대로 정확히 문자열 경로를 출력했습니다.

이것은 단지 텍스트 출력일 뿐입니다. Nextflow는 아직 이것으로 특별한 작업을 하지 않았습니다.
또한 Nextflow의 관점에서 이것은 단지 문자열(`java.lang.String` 클래스)일 뿐이라는 것을 확인했습니다.
이것은 아직 파일에 해당한다고 Nextflow에 알리지 않았기 때문에 당연합니다.

### 1.2. file()로 Path 객체 생성하기

경로 문자열에서 [Path 객체](https://www.nextflow.io/docs/latest/reference/stdlib-types.html#path)를 생성하여 Nextflow에 파일을 처리하는 방법을 알려줄 수 있습니다.

워크플로우에서 `file()` 메서드를 사용하여 문자열 경로 `data/patientA_rep1_normal_R1_001.fastq.gz`를 Path 객체로 변환할 수 있으며, 이는 파일 속성 및 작업에 대한 접근을 제공합니다.

다음과 같이 `main.nf`를 편집하여 문자열을 `file()`로 감싸세요:

=== "후"

    ```groovy title="main.nf" linenums="5" hl_lines="2"
        // Create a Path object from a string path
        myFile = file('data/patientA_rep1_normal_R1_001.fastq.gz')

        println "${myFile} is of class ${myFile.class}"
    ```

=== "전"

    ```groovy title="main.nf" linenums="5" hl_lines="2"
        // Create a Path object from a string path
        myFile = 'data/patientA_rep1_normal_R1_001.fastq.gz'

        println "${myFile} is of class ${myFile.class}"
    ```

이제 워크플로우를 다시 실행하세요:

```bash
nextflow run main.nf
```

??? success "명령 출력"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `main.nf` [kickass_coulomb] DSL2 - revision: 5af44b1b59

    /workspaces/training/side-quests/working_with_files/data/patientA_rep1_normal_R1_001.fastq.gz is of class class sun.nio.fs.UnixPath
    ```

이번에는 입력으로 제공한 상대 경로 대신 전체 절대 경로가 표시됩니다.

Nextflow는 우리의 문자열을 Path 객체로 변환하고 시스템의 실제 파일 위치로 해석했습니다.
파일 경로는 이제 `/workspaces/training/side-quests/working_with_files/data/patientA_rep1_normal_R1_001.fastq.gz`와 같이 절대 경로가 됩니다.

또한 Path 객체 클래스가 `sun.nio.fs.UnixPath`임을 주목하세요. 이것은 Nextflow가 로컬 파일을 나타내는 방식입니다.
나중에 보겠지만, 원격 파일은 다른 클래스 이름(예: HTTP 파일의 경우 `nextflow.file.http.XPath`)을 가지지만, 모두 정확히 같은 방식으로 작동하며 워크플로우에서 동일하게 사용할 수 있습니다.

!!! tip

    **주요 차이점:**

    - **경로 문자열**: Nextflow가 문자로 취급하는 단순한 텍스트
    - **Path 객체**: Nextflow가 작업할 수 있는 스마트 파일 참조

    이렇게 생각해 보세요: 경로 문자열은 종이에 주소를 쓰는 것과 같고, Path 객체는 GPS 장치에 주소를 로드하여 그곳으로 이동하는 방법을 알고 여정에 대한 세부 정보를 알려줄 수 있는 것과 같습니다.

### 1.3. 파일 속성 접근하기

왜 이것이 유용할까요? 이제 Nextflow가 `myFile`이 단순한 문자열이 아니라 Path 객체라는 것을 이해하므로 Path 객체의 다양한 속성에 접근할 수 있습니다.

워크플로우를 업데이트하여 내장 파일 속성을 출력해 봅시다:

=== "후"

    ```groovy title="main.nf" linenums="5" hl_lines="4-9"
        // Create a Path object from a string path
        myFile = file('data/patientA_rep1_normal_R1_001.fastq.gz')

        // Print file attributes
        println "File object class: ${myFile.class}"
        println "File name: ${myFile.name}"
        println "Simple name: ${myFile.simpleName}"
        println "Extension: ${myFile.extension}"
        println "Parent directory: ${myFile.parent}"
    ```

=== "전"

    ```groovy title="main.nf" linenums="5" hl_lines="4"
        // Create a Path object from a string path
        myFile = file('data/patientA_rep1_normal_R1_001.fastq.gz')

        println "${myFile} is of class ${myFile.class}"
    ```

워크플로우를 실행하세요:

```bash
nextflow run main.nf
```

??? success "명령 출력"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `main.nf` [ecstatic_ampere] DSL2 - revision: f3fa3dcb48

    File object class: sun.nio.fs.UnixPath
    File name: patientA_rep1_normal_R1_001.fastq.gz
    Simple name: patientA_rep1_normal_R1_001
    Extension: gz
    Parent directory: /workspaces/training/side-quests/working_with_files/data
    ```

위의 콘솔에 다양한 파일 속성이 출력된 것을 볼 수 있습니다.

### 1.4. 프로세스에 파일 전달하기

문자열과 Path 객체의 차이는 프로세스가 있는 실제 워크플로우를 구축하기 시작할 때 중요해집니다.
지금까지 Nextflow가 입력 파일을 파일로 취급하고 있음을 확인했지만, 프로세스에서 해당 파일에 대해 실제로 무언가를 실행할 수 있는지 확인해 봅시다.

#### 1.4.1. 프로세스 가져오기 및 코드 검토

파일 입력을 받아 몇 줄이 있는지 세는 `COUNT_LINES`라는 미리 작성된 프로세스 모듈을 제공합니다.

워크플로우에서 프로세스를 사용하려면 워크플로우 블록 앞에 include 문을 추가하기만 하면 됩니다:

=== "후"

    ```groovy title="main.nf" linenums="1" hl_lines="3"
    #!/usr/bin/env nextflow

    include { COUNT_LINES } from './modules/count_lines.nf'

    workflow {
    ```

=== "전"

    ```groovy title="main.nf" linenums="1"
    #!/usr/bin/env nextflow

    workflow {
    ```

모듈 파일을 열어 코드를 검토할 수 있습니다:

```groovy title="modules/count_lines.nf" linenums="1"
#!/usr/bin/env nextflow

process COUNT_LINES {
    debug true

    input:
    path input_file

    script:
    """
    set -o pipefail
    echo "Processing file: $input_file"
    gzip -dc $input_file | wc -l
    """
}
```

보시다시피, 파일의 압축을 풀고 포함된 줄 수를 세는 매우 간단한 작은 스크립트입니다.

??? info "`debug true`는 무엇을 하나요?"

    프로세스 정의의 `debug true` 지시문은 Nextflow가 스크립트의 출력(예: 줄 수 "40")을 실행 로그에 직접 출력하도록 합니다.
    이것이 없으면 프로세스 실행 상태만 볼 수 있고 스크립트의 실제 출력은 볼 수 없습니다.

    Nextflow 프로세스 디버깅에 대한 자세한 내용은 [Debugging Nextflow Workflows](debugging.md) 사이드 퀘스트를 참조하세요.

#### 1.4.2. `COUNT_LINES` 호출 추가

이제 프로세스를 워크플로우에서 사용할 수 있으므로 입력 파일에서 실행하기 위해 `COUNT_LINES` 프로세스에 대한 호출을 추가할 수 있습니다.

워크플로우를 다음과 같이 편집하세요:

=== "후"

    ```groovy title="main.nf" linenums="7" hl_lines="11-12"
        // Create a Path object from a string path
        myFile = file('data/patientA_rep1_normal_R1_001.fastq.gz')

        // Print file attributes
        println "File object class: ${myFile.class}"
        println "File name: ${myFile.name}"
        println "Simple name: ${myFile.simpleName}"
        println "Extension: ${myFile.extension}"
        println "Parent directory: ${myFile.parent}"

        // Count the lines in the file
        COUNT_LINES(myFile)
    ```

=== "전"

    ```groovy title="main.nf" linenums="7" hl_lines="4-9"
        // Create a Path object from a string path
        myFile = file('data/patientA_rep1_normal_R1_001.fastq.gz')

        // Print file attributes
        println "File object class: ${myFile.class}"
        println "File name: ${myFile.name}"
        println "Simple name: ${myFile.simpleName}"
        println "Extension: ${myFile.extension}"
        println "Parent directory: ${myFile.parent}"
    ```

이제 워크플로우를 실행하세요:

```bash
nextflow run main.nf
```

??? success "명령 출력"

    ```console
     N E X T F L O W   ~  version 25.10.2

    Launching `main.nf` [cheeky_hypatia] DSL2 - revision: 281d13c414

    File object class: class sun.nio.fs.UnixPath
    File name: patientA_rep1_normal_R1_001.fastq.gz
    Simple name: patientA_rep1_normal_R1_001
    Extension: gz
    Parent directory: /workspaces/training/side-quests/working_with_files/data
    executor >  local (1)
    [e9/341c05] COUNT_LINES [100%] 1 of 1 ✔
    Processing file: /workspaces/training/side-quests/working_with_files/data/patientA_rep1_normal_R1_001.fastq.gz
    40
    ```

이것은 프로세스 내에서 파일을 적절하게 작업할 수 있음을 보여줍니다.

구체적으로 Nextflow는 다음 작업을 성공적으로 수행했습니다:

- 작업 디렉토리에 파일 스테이징
- .gz 파일 압축 해제
- 줄 수 세기(이 경우 40줄)
- 오류 없이 완료

이 원활한 작업의 핵심은 입력이 파일이며 그렇게 취급되어야 한다고 Nextflow에 명시적으로 알리는 것입니다.

### 1.5. 기본 파일 입력 오류 문제 해결

이것은 종종 Nextflow 초보자를 혼란스럽게 하므로 잘못했을 때 어떤 일이 발생하는지 몇 분 동안 살펴보겠습니다.

파일 처리를 잘못할 수 있는 두 가지 주요 위치가 있습니다: 워크플로우 수준과 프로세스 수준입니다.

#### 1.5.1. 워크플로우 수준 오류

워크플로우 블록에서 입력을 지정할 때 파일을 문자열로 취급하도록 되돌리면 어떤 일이 발생하는지 살펴봅시다.

워크플로우를 다음과 같이 편집하고 경로별 출력 문을 주석 처리하세요:

=== "후"

    ```groovy title="main.nf" linenums="7" hl_lines="2 6-11"
        // Create a Path object from a string path
        myFile = 'data/patientA_rep1_normal_R1_001.fastq.gz'

        // Print file attributes
        println "File object class: ${myFile.class}"
        /*
        println "File name: ${myFile.name}"
        println "Simple name: ${myFile.simpleName}"
        println "Extension: ${myFile.extension}"
        println "Parent directory: ${myFile.parent}"
        */

        // Count the lines in the file
        COUNT_LINES(myFile)
    ```

=== "전"

    ```groovy title="main.nf" linenums="7" hl_lines="4-9"
        // Create a Path object from a string path
        myFile = file('data/patientA_rep1_normal_R1_001.fastq.gz')

        // Print file attributes
        println "File object class: ${myFile.class}"
        println "File name: ${myFile.name}"
        println "Simple name: ${myFile.simpleName}"
        println "Extension: ${myFile.extension}"
        println "Parent directory: ${myFile.parent}"

        // Count the lines in the file
        COUNT_LINES(myFile)
    ```

이제 워크플로우를 실행하세요:

```bash
nextflow run main.nf
```

??? failure "명령 출력"

    ```console
     N E X T F L O W   ~  version 25.10.2

    Launching `main.nf` [friendly_goodall] DSL2 - revision: ae50609b20

    [-        ] COUNT_LINES -
    ERROR ~ Error executing process > 'COUNT_LINES'

    Caused by:
      Not a valid path value: 'data/patientA_rep1_normal_R1_001.fastq.gz'



    Tip: view the complete command output by changing to the process work dir and entering the command `cat .command.out`

    -- Check '.nextflow.log' file for details
    ```

중요한 부분은 다음과 같습니다:

```console
Not a valid path value: 'data/patientA_rep1_normal_R1_001.fastq.gz'
```

`path` 입력을 지정하면 Nextflow는 단순한 문자열이 아닌 실제 파일 참조를 전달하고 있는지 검증합니다.
이 오류는 `'data/patientA_rep1_normal_R1_001.fastq.gz'`가 Path 객체가 아니라 문자열이기 때문에 유효한 경로 값이 아니라고 알려줍니다.

Nextflow는 즉시 문제를 감지하고 프로세스를 시작하기 전에 중지했습니다.

#### 1.5.2. 프로세스 수준 오류

입력을 파일로 취급하도록 Nextflow에 지정하는 것을 잊을 수 있는 다른 위치는 프로세스 정의입니다.

!!! warning "1.5.1의 워크플로우 오류 유지"

    이 테스트가 올바르게 작동하려면 워크플로우를 손상된 상태로 유지하세요(`file()` 대신 일반 문자열 사용).
    프로세스의 `val`과 결합하면 아래에 표시된 오류가 발생합니다.

모듈을 다음과 같이 편집하세요:

=== "후"

    ```groovy title="modules/count_lines.nf" linenums="3" hl_lines="5"
    process COUNT_LINES {
        debug true

        input:
        val input_file
    ```

=== "전"

    ```groovy title="modules/count_lines.nf" linenums="3" hl_lines="5"
    process COUNT_LINES {
        debug true

        input:
        path input_file
    ```

이제 워크플로우를 다시 실행하세요:

```bash
nextflow run main.nf
```

??? failure "명령 출력"

    ```console
     N E X T F L O W   ~  version 25.10.2

    Launching `main.nf` [soggy_golick] DSL2 - revision: ae50609b20

    executor >  local (1)
    [b3/b3023c] COUNT_LINES [  0%] 0 of 1 ✘
    ERROR ~ Error executing process > 'COUNT_LINES'

    Caused by:
      Process `COUNT_LINES` terminated with an error exit status (1)


    Command executed:

      set -o pipefail
      echo "Processing file: data/patientA_rep1_normal_R1_001.fastq.gz"
      gzip -dc data/patientA_rep1_normal_R1_001.fastq.gz | wc -l

    Command exit status:
      1

    Command output:
      Processing file: data/patientA_rep1_normal_R1_001.fastq.gz
      0

    Command error:
      Processing file: data/patientA_rep1_normal_R1_001.fastq.gz
      gzip: data/patientA_rep1_normal_R1_001.fastq.gz: No such file or directory
      0

    Work dir:
      /workspaces/training/side-quests/working_with_files/work/b3/b3023cb2ccb986851301d8e369e79f

    Tip: when you have fixed the problem you can continue the execution adding the option `-resume` to the run command line

    -- Check '.nextflow.log' file for details
    ```

프로세스가 위에서 언급한 대로 디버깅 정보를 출력하도록 설정되어 있기 때문에 오류에 대한 많은 세부 정보가 표시됩니다.

가장 관련성이 높은 섹션은 다음과 같습니다:

```console
Command executed:

  set -o pipefail
  echo "Processing file: data/patientA_rep1_normal_R1_001.fastq.gz"
  gzip -dc data/patientA_rep1_normal_R1_001.fastq.gz | wc -l
```

```console
Command error:
  Processing file: data/patientA_rep1_normal_R1_001.fastq.gz
  gzip: data/patientA_rep1_normal_R1_001.fastq.gz: No such file or directory
  0
```

시스템이 파일을 찾을 수 없다고 나와 있지만, 경로를 확인하면 해당 위치에 해당 이름의 파일이 있습니다.

이것을 실행했을 때 Nextflow는 문자열 값을 스크립트로 전달했지만 실제 파일을 작업 디렉토리에 _스테이징_하지 않았습니다.
따라서 프로세스는 상대 문자열 `data/patientA_rep1_normal_R1_001.fastq.gz`를 사용하려고 했지만 해당 파일은 프로세스 작업 디렉토리 내에 존재하지 않습니다.

이 두 예제를 종합하면 입력을 파일로 처리해야 한다고 Nextflow에 알리는 것이 얼마나 중요한지 알 수 있습니다.

!!! note

    다음 섹션으로 계속하기 전에 두 가지 의도적인 오류를 모두 수정하세요.

### 핵심 정리

- 경로 문자열 vs Path 객체: 문자열은 단순한 텍스트이고, Path 객체는 스마트 파일 참조입니다
- `file()` 메서드는 문자열 경로를 Nextflow가 작업할 수 있는 Path 객체로 변환합니다
- `name`, `simpleName`, `extension`, `parent`와 같은 파일 속성에 [파일 속성을 사용하여](https://www.nextflow.io/docs/latest/working-with-files.html#getting-file-attributes) 접근할 수 있습니다
- 문자열 대신 Path 객체를 사용하면 Nextflow가 워크플로우에서 파일을 적절하게 관리할 수 있습니다
- 프로세스 입력 결과: 적절한 파일 처리를 위해서는 문자열이 아닌 Path 객체가 필요하며, 이를 통해 파일이 올바르게 스테이징되고 프로세스에서 사용할 수 있습니다

---

## 2. 원격 파일 사용하기

Nextflow의 주요 기능 중 하나는 로컬 파일(같은 머신에 있는)에서 인터넷을 통해 접근 가능한 원격 파일로 원활하게 전환할 수 있는 능력입니다.

올바르게 수행하면 다른 위치에서 오는 파일을 수용하기 위해 워크플로우의 로직을 변경할 필요가 없습니다.
원격 파일을 사용하기 위해 해야 할 일은 워크플로우에 제공할 때 파일 경로에 적절한 접두사를 지정하는 것뿐입니다.

예를 들어, `/path/to/data`에는 접두사가 없어 '일반' 로컬 파일 경로임을 나타내는 반면, `s3://path/to/data`에는 `s3://` 접두사가 포함되어 Amazon의 S3 객체 스토리지에 위치함을 나타냅니다.

다양한 프로토콜이 지원됩니다:

- HTTP(S)/FTP (http://, https://, ftp://)
- Amazon S3 (s3://)
- Azure Blob Storage (az://)
- Google Cloud Storage (gs://)

이들 중 하나를 사용하려면 문자열에 관련 접두사를 지정하기만 하면 되며, 이는 기술적으로 파일 경로 대신 URI(Uniform Resource Identifier)라고 불립니다.
Nextflow는 인증 및 파일을 올바른 위치로 스테이징하고, 다운로드 또는 업로드 및 예상되는 모든 기타 파일 작업을 처리합니다.

이 시스템의 주요 강점은 파이프라인 로직을 변경하지 않고도 환경 간에 전환할 수 있다는 것입니다.
예를 들어, URI를 변경하기만 하면 작은 로컬 테스트 세트로 개발한 후 원격 스토리지에 있는 전체 규모 테스트 세트로 전환할 수 있습니다.

### 2.1. 인터넷의 파일 사용하기

Github에 저장된 동일한 데이터의 복사본을 가리키는 HTTPS 경로로 워크플로우에 제공하는 로컬 경로를 전환하여 이를 테스트해 봅시다.

!!! warning

    이것은 활성 인터넷 연결이 있는 경우에만 작동합니다.

`main.nf`를 다시 열고 다음과 같이 입력 경로를 변경하세요:

=== "후"

    ```groovy title="main.nf" linenums="2" hl_lines="2"
        // Using a remote file from the internet
        myFile = file('https://raw.github.com/nextflow-io/training/master/side-quests/working_with_files/data/patientA_rep1_normal_R1_001.fastq.gz')

        // Print file attributes
        println "File object class: ${myFile.class}"
        println "File name: ${myFile.name}"
        println "Simple name: ${myFile.simpleName}"
        println "Extension: ${myFile.extension}"
        println "Parent directory: ${myFile.parent}"
    ```

=== "전"

    ```groovy title="main.nf" linenums="2" hl_lines="2"
        // Create a Path object from a string path
        myFile = file('data/patientA_rep1_normal_R1_001.fastq.gz')

        // Print file attributes
        println "File object class: ${myFile.class}"
        println "File name: ${myFile.name}"
        println "Simple name: ${myFile.simpleName}"
        println "Extension: ${myFile.extension}"
        println "Parent directory: ${myFile.parent}"
    ```

워크플로우를 실행해 봅시다:

```bash
nextflow run main.nf
```

??? success "명령 출력"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `main.nf` [insane_swartz] DSL2 - revision: fff18abe6d

    File object class: class nextflow.file.http.XPath
    File name: patientA_rep1_normal_R1_001.fastq.gz
    Simple name: patientA_rep1_normal_R1_001
    Extension: gz
    Parent directory: /nextflow-io/training/master/side-quests/working_with_files/data
    executor >  local (1)
    [8a/2ab7ca] COUNT_LINES [100%] 1 of 1 ✔
    Processing file: patientA_rep1_normal_R1_001.fastq.gz
    40
    ```

작동합니다! 콘솔 출력에서 거의 변경되지 않은 것을 볼 수 있습니다.

콘솔 출력의 한 가지 차이점은 경로 객체 클래스가 이제 `nextflow.file.http.XPath`인 반면, 로컬 경로의 경우 클래스가 `sun.nio.fs.UnixPath`였다는 것입니다.
이러한 클래스를 기억할 필요는 없습니다. Nextflow가 다른 위치를 적절하게 식별하고 처리한다는 것을 보여주기 위해 언급할 뿐입니다.

내부적으로 Nextflow는 작업 디렉토리 내에 있는 스테이징 디렉토리로 파일을 다운로드했습니다.
그런 다음 스테이징된 파일을 로컬 파일로 취급하고 관련 프로세스 디렉토리에 심볼릭 링크할 수 있습니다.

프로세스의 해시 값에 있는 작업 디렉토리의 내용을 살펴보면 이것이 발생했음을 확인할 수 있습니다.

??? abstract "작업 디렉토리 내용"

    프로세스 해시가 `8a/2ab7ca`인 경우 작업 디렉토리를 탐색할 수 있습니다:

    ```console
    $ ls -la work/8a/2ab7ca*/
    total 16
    drwxr-xr-x  6 user  staff   192 Jan 28 10:00 .
    drwxr-xr-x  3 user  staff    96 Jan 28 10:00 ..
    -rw-r--r--  1 user  staff     0 Jan 28 10:00 .command.begin
    -rw-r--r--  1 user  staff   127 Jan 28 10:00 .command.sh
    lrwxr-xr-x  1 user  staff    89 Jan 28 10:00 patientA_rep1_normal_R1_001.fastq.gz -> /path/to/work/stage/.../patientA_rep1_normal_R1_001.fastq.gz
    ```

    심볼릭 링크는 Nextflow가 자동으로 다운로드한 원격 파일의 스테이징된 복사본을 가리킵니다.

더 큰 파일의 경우 다운로드 단계가 로컬 파일에서 실행하는 것보다 추가 시간이 걸립니다.
그러나 Nextflow는 불필요한 다운로드를 피하기 위해 이미 스테이징된 복사본이 있는지 확인합니다.
따라서 동일한 파일에서 다시 실행하고 스테이징된 파일을 삭제하지 않은 경우 Nextflow는 스테이징된 복사본을 사용합니다.

이것은 Nextflow를 사용하여 로컬 및 원격 데이터 간에 전환하는 것이 얼마나 쉬운지 보여주며, 이는 Nextflow의 주요 기능입니다.

!!! note

    이 원칙의 한 가지 중요한 예외는 HTTPS가 여러 파일을 나열할 수 없기 때문에 HTTPS에서는 glob 패턴이나 디렉토리 경로를 사용할 수 없으며 정확한 파일 URL을 지정해야 한다는 것입니다.
    그러나 blob 스토리지(`s3://`, `az://`, `gs://`)와 같은 다른 스토리지 프로토콜은 glob 및 디렉토리 경로를 모두 사용할 수 있습니다.

    클라우드 스토리지에서 glob 패턴을 사용하는 방법은 다음과 같습니다:

    ```groovy title="Cloud storage examples (not runnable in this environment)"
    // S3 with glob patterns - would match multiple files
    ch_s3_files = channel.fromPath('s3://my-bucket/data/*.fastq.gz')

    // Azure Blob Storage with glob patterns
    ch_azure_files = channel.fromPath('az://container/data/patient*_R{1,2}.fastq.gz')

    // Google Cloud Storage with glob patterns
    ch_gcs_files = channel.fromPath('gs://bucket/data/sample_*.fastq.gz')
    ```

    다음 섹션에서 glob을 실제로 사용하는 방법을 보여드리겠습니다.

### 2.2. 로컬 파일로 다시 전환

이 사이드 퀘스트의 나머지 부분에서는 로컬 예제 파일을 사용할 것이므로 워크플로우 입력을 원래 파일로 다시 전환하겠습니다:

=== "후"

    ```groovy title="main.nf" linenums="2" hl_lines="2"
        // Create a Path object from a string path
        myFile = file('data/patientA_rep1_normal_R1_001.fastq.gz')

        // Print file attributes
        println "File object class: ${myFile.class}"
        println "File name: ${myFile.name}"
        println "Simple name: ${myFile.simpleName}"
        println "Extension: ${myFile.extension}"
        println "Parent directory: ${myFile.parent}"
    ```

=== "전"

    ```groovy title="main.nf" linenums="2" hl_lines="2"
        // Create a Path object from a string path
        myFile = file('https://raw.github.com/nextflow-io/training/master/side-quests/working_with_files/data/patientA_rep1_normal_R1_001.fastq.gz')

        // Print file attributes
        println "File object class: ${myFile.class}"
        println "File name: ${myFile.name}"
        println "Simple name: ${myFile.simpleName}"
        println "Extension: ${myFile.extension}"
        println "Parent directory: ${myFile.parent}"
    ```

### 핵심 정리

- 원격 데이터는 URI(HTTP, FTP, S3, Azure, Google Cloud)를 사용하여 접근됩니다
- 이러한 경로가 프로세스에 전달되는 한 Nextflow는 자동으로 데이터를 올바른 위치로 다운로드하고 스테이징합니다
- 원격 파일을 다운로드하거나 업로드하는 로직을 작성하지 마세요!
- 로컬 및 원격 파일은 다른 객체 유형을 생성하지만 동일하게 작동합니다
- **중요**: HTTP/HTTPS는 단일 파일에서만 작동합니다(glob 패턴 없음)
- 클라우드 스토리지(S3, Azure, GCS)는 단일 파일과 glob 패턴을 모두 지원합니다
- 프로토콜이 필요한 작업을 지원하는 한 코드 로직을 변경하지 않고도 로컬 및 원격 데이터 소스 간에 원활하게 전환할 수 있습니다

---

## 3. `fromPath()` 채널 팩토리 사용하기

지금까지 한 번에 하나의 파일로 작업했지만 Nextflow에서는 일반적으로 처리할 여러 입력 파일이 있는 입력 채널을 생성하려고 합니다.

이를 수행하는 순진한 방법은 다음과 같이 `file()` 메서드를 [`channel.of()`](https://www.nextflow.io/docs/latest/reference/channel.html#of)와 결합하는 것입니다:

```groovy title="Syntax example"
ch_files = channel.of([file('data/patientA_rep1_normal_R1_001.fastq.gz')],
                      [file('data/patientA_rep1_normal_R1_001.fastq.gz')])
```

이것은 작동하지만 번거롭습니다.

!!! tip "`file()` vs `channel.fromPath()` 사용 시기"

    - 단일 Path 객체를 직접 조작해야 할 때(파일이 존재하는지 확인, 속성 읽기 또는 단일 프로세스 호출에 전달) `file()`을 사용하세요
    - 여러 파일을 보유할 수 있는 채널이 필요할 때, 특히 glob 패턴을 사용하거나 파일이 여러 프로세스를 통해 흐를 때 `channel.fromPath()`를 사용하세요

여기서 [`channel.fromPath()`](https://www.nextflow.io/docs/latest/reference/channel.html#frompath)가 등장합니다: 하나 이상의 정적 파일 문자열과 glob 패턴에서 채널을 생성하는 데 필요한 모든 기능을 번들로 제공하는 편리한 채널 팩토리입니다.

### 3.1. 채널 팩토리 추가

워크플로우를 업데이트하여 `channel.fromPath`를 사용해 봅시다.

=== "후"

    ```groovy title="main.nf" linenums="7" hl_lines="1-3"
        // Load files with channel.fromPath
        ch_files = channel.fromPath('data/patientA_rep1_normal_R1_001.fastq.gz')
        ch_files.view { myFile -> "Found file: $myFile" }

        // Print file attributes
        /* 지금은 주석 처리하고 나중에 다시 돌아오겠습니다!
        println "File object class: ${myFile.class}"
        println "File name: ${myFile.name}"
        println "Simple name: ${myFile.simpleName}"
        println "Extension: ${myFile.extension}"
        println "Parent directory: ${myFile.parent}"
        */

        // Count the lines in the file
        // COUNT_LINES(myFile)
    ```

=== "전"

    ```groovy title="main.nf" linenums="7" hl_lines="1-2"
        // Create a Path object from a string path
        myFile = file('data/patientA_rep1_normal_R1_001.fastq.gz')

        // Print file attributes
        println "File object class: ${myFile.class}"
        println "File name: ${myFile.name}"
        println "Simple name: ${myFile.simpleName}"
        println "Extension: ${myFile.extension}"
        println "Parent directory: ${myFile.parent}"

        // Count the lines in the file
        COUNT_LINES(myFile)
    ```

지금은 속성을 출력하는 코드를 주석 처리하고 대신 파일명만 출력하는 `.view` 문을 추가했습니다.

워크플로우를 실행하세요:

```bash
nextflow run main.nf
```

??? success "명령 출력"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `main.nf` [grave_meucci] DSL2 - revision: b09964a583

    Found file: /workspaces/training/side-quests/working_with_files/data/patientA_rep1_normal_R1_001.fastq.gz
    ```

보시다시피 파일 경로가 채널에 `Path` 유형 객체로 로드되고 있습니다.
이것은 `file()`이 했을 것과 유사하지만 이제 원하는 경우 더 많은 파일을 로드할 수 있는 채널이 있습니다.

`channel.fromPath()`를 사용하는 것은 파일 목록으로 채워진 새 채널을 생성하는 편리한 방법입니다.

### 3.2. 채널의 파일 속성 보기

채널 팩토리를 처음 사용할 때 코드를 단순화하고 파일 이름만 출력했습니다.

전체 파일 속성을 출력하는 것으로 돌아가 봅시다:

=== "후"

    ```groovy title="main.nf" linenums="7" hl_lines="3-9 12"
        // Load files with channel.fromPath
        ch_files = channel.fromPath('data/patientA_rep1_normal_R1_001.fastq.gz')
        ch_files.view { myFile ->
            println "File object class: ${myFile.class}"
            println "File name: ${myFile.name}"
            println "Simple name: ${myFile.simpleName}"
            println "Extension: ${myFile.extension}"
            println "Parent directory: ${myFile.parent}"
        }

        // Count the lines in the file
        COUNT_LINES(ch_files)
    ```

=== "전"

    ```groovy title="main.nf" linenums="7" hl_lines="3"
        // Load files with channel.fromPath
        ch_files = channel.fromPath('data/patientA_rep1_normal_R1_001.fastq.gz')
        ch_files.view { myFile -> "Found file: $myFile" }

        // Count the lines in the file
        // COUNT_LINES(ch_files)
    ```

또한 채널 기반 접근 방식으로 파일 처리가 여전히 올바르게 작동하는지 확인하기 위해 `COUNT_LINES` 프로세스 호출을 다시 활성화하고 있습니다.

워크플로우를 실행하세요:

```bash
nextflow run main.nf
```

??? success "명령 출력"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `main.nf` [furious_swanson] DSL2 - revision: c35c34950d

    executor >  local (1)
    [9d/6701a6] COUNT_LINES (1) [100%] 1 of 1 ✔
    File object class: sun.nio.fs.UnixPath
    File name: patientA_rep1_normal_R1_001.fastq.gz
    Simple name: patientA_rep1_normal_R1_001
    Extension: gz
    Parent directory: /workspaces/training/side-quests/working_with_files/data
    Processing file: patientA_rep1_normal_R1_001.fastq.gz
    40
    ```

이전과 동일한 결과이지만 이제 파일이 채널에 있으므로 더 추가할 수 있습니다.

### 3.3. glob을 사용하여 여러 파일 매칭

채널에 더 많은 파일을 로드할 수 있는 여러 가지 방법이 있습니다.
여기서는 와일드카드 문자를 기반으로 파일 및 디렉토리 이름을 매칭하고 검색하는 편리한 방법인 glob 패턴을 사용하는 방법을 보여드리겠습니다.
이러한 패턴을 매칭하는 프로세스를 "globbing" 또는 "파일명 확장"이라고 합니다.

!!! note

    앞서 언급했듯이 Nextflow는 HTTPS 파일 경로를 제외한 대부분의 경우 입력 및 출력 파일을 관리하기 위해 globbing을 지원합니다. HTTPS는 여러 파일을 나열할 수 없기 때문입니다.

특정 환자 `patientA`와 관련된 파일 쌍의 두 파일을 모두 검색하고 싶다고 가정해 봅시다:

```console
patientA_rep1_normal_R1_001.fastq.gz
patientA_rep1_normal_R2_001.fastq.gz
```

파일명 간의 유일한 차이점이 `R` 뒤의 반복 번호이므로 와일드카드 문자 `*`를 사용하여 다음과 같이 숫자를 대신할 수 있습니다:

```console
patientA_rep1_normal_R*_001.fastq.gz
```

이것이 우리가 필요한 glob 패턴입니다.

이제 채널 팩토리의 파일 경로를 다음과 같이 해당 glob 패턴을 사용하도록 업데이트하기만 하면 됩니다:

=== "후"

    ```groovy title="main.nf" linenums="7"
      // Load files with channel.fromPath
        ch_files = channel.fromPath('data/patientA_rep1_normal_R*_001.fastq.gz')
    ```

=== "전"

    ```groovy title="main.nf" linenums="7"
      // Load files with channel.fromPath
        ch_files = channel.fromPath('data/patientA_rep1_normal_R1_001.fastq.gz')
    ```

Nextflow는 이것이 glob 패턴임을 자동으로 인식하고 적절하게 처리합니다.

워크플로우를 실행하여 테스트하세요:

```bash
nextflow run main.nf
```

??? success "명령 출력"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `main.nf` [boring_sammet] DSL2 - revision: d2aa789c9a

    executor >  local (2)
    [3c/a65de5] COUNT_LINES (2) [100%] 2 of 2 ✔
    File object class: class sun.nio.fs.UnixPath
    File name: patientA_rep1_normal_R1_001.fastq.gz
    Simple name: patientA_rep1_normal_R1_001
    Extension: gz
    Parent directory: /workspaces/training/side-quests/working_with_files/data
    File object class: class sun.nio.fs.UnixPath
    File name: patientA_rep1_normal_R2_001.fastq.gz
    Simple name: patientA_rep1_normal_R2_001
    Extension: gz
    Parent directory: /workspaces/training/side-quests/working_with_files/data
    Processing file: patientA_rep1_normal_R1_001.fastq.gz
    40

    Processing file: patientA_rep1_normal_R2_001.fastq.gz
    40
    ```

보시다시피 이제 채널에 두 개의 Path 객체가 있으며, 이는 Nextflow가 파일명 확장을 올바르게 수행했고 예상대로 두 파일을 모두 로드하고 처리했음을 보여줍니다.

이 방법을 사용하면 glob 패턴을 변경하기만 하면 원하는 만큼 많거나 적은 파일을 검색할 수 있습니다. 예를 들어 파일명의 모든 가변 부분을 `*`로 대체하여 더 관대하게 만들면(_예:_ `data/patient*_rep*_*_R*_001.fastq.gz`) `data` 디렉토리의 모든 예제 파일을 가져올 수 있습니다.

### 핵심 정리

- `channel.fromPath()`는 패턴과 일치하는 파일로 채널을 생성합니다
- 각 파일은 채널의 별도 요소로 방출됩니다
- glob 패턴을 사용하여 여러 파일을 매칭할 수 있습니다
- 파일은 전체 속성을 가진 Path 객체로 자동 변환됩니다
- `.view()` 메서드를 사용하면 채널 내용을 검사할 수 있습니다

---

## 4. 파일명에서 기본 메타데이터 추출하기

대부분의 과학 분야에서 데이터를 포함하는 파일 이름에 메타데이터가 인코딩되어 있는 것이 매우 일반적입니다.
예를 들어, 생물정보학에서 시퀀싱 데이터를 포함하는 파일은 종종 샘플, 조건, 반복 및 리드 번호에 대한 정보를 인코딩하는 방식으로 명명됩니다.

파일명이 일관된 규칙에 따라 구성된 경우 표준화된 방식으로 해당 메타데이터를 추출하여 분석 과정에서 사용할 수 있습니다.
물론 이것은 큰 '만약'이며, 파일명 구조에 의존할 때마다 매우 신중해야 합니다. 하지만 현실은 이 접근 방식이 매우 널리 사용되므로 Nextflow에서 어떻게 수행되는지 살펴보겠습니다.

예제 데이터의 경우 파일명에 일관되게 구조화된 메타데이터가 포함되어 있다는 것을 알고 있습니다.
예를 들어, 파일명 `patientA_rep1_normal_R2_001`은 다음을 인코딩합니다:

- 환자 ID: `patientA`
- 반복 ID: `rep1`
- 샘플 유형: `normal`(`tumor`와 반대)
- 리드 세트: `R1`(`R2`와 반대)

세 단계로 이 정보를 검색하도록 워크플로우를 수정할 것입니다:

1. 메타데이터를 포함하는 파일의 `simpleName` 검색
2. `tokenize()`라는 메서드를 사용하여 메타데이터 분리
3. 맵을 사용하여 메타데이터 구성

!!! warning

    환자 이름이나 기타 식별 특성과 같은 민감한 정보를 파일명에 인코딩해서는 안 됩니다. 이는 환자 개인 정보 보호 또는 기타 관련 보안 제한을 손상시킬 수 있습니다.

### 4.1. `simpleName` 검색

`simpleName`은 경로와 확장자가 제거된 파일명에 해당하는 파일 속성입니다.

워크플로우를 다음과 같이 편집하세요:

=== "후"

    ```groovy title="main.nf" linenums="7" hl_lines="3-6"
        // Load files with channel.fromPath
        ch_files = channel.fromPath('data/patientA_rep1_normal_R*_001.fastq.gz')
        ch_files.map { myFile ->
            [ myFile.simpleName, myFile ]
        }
        .view()
    ```

=== "전"

    ```groovy title="main.nf" linenums="7" hl_lines="3-9"
        // Load files with channel.fromPath
        ch_files = channel.fromPath('data/patientA_rep1_normal_R*_001.fastq.gz')
        ch_files.view { myFile ->
            println "File object class: ${myFile.class}"
            println "File name: ${myFile.name}"
            println "Simple name: ${myFile.simpleName}"
            println "Extension: ${myFile.extension}"
            println "Parent directory: ${myFile.parent}"
        }
    ```

이것은 `simpleName`을 검색하고 `map()` 작업을 사용하여 전체 파일 객체와 연결합니다.

워크플로우를 실행하여 작동하는지 테스트하세요:

```bash
nextflow run main.nf
```

??? success "명령 출력"

    ```console hl_lines="7-8"
    N E X T F L O W   ~  version 25.10.2

    Launching `main.nf` [suspicious_mahavira] DSL2 - revision: ae8edc4e48

    executor >  local (2)
    [e9/55774b] COUNT_LINES (2) [100%] 2 of 2 ✔
    [patientA_rep1_normal_R2_001, /workspaces/training/side-quests/working_with_files/data/patientA_rep1_normal_R2_001.fastq.gz]
    [patientA_rep1_normal_R1_001, /workspaces/training/side-quests/working_with_files/data/patientA_rep1_normal_R1_001.fastq.gz]
    Processing file: patientA_rep1_normal_R1_001.fastq.gz
    40

    Processing file: patientA_rep1_normal_R2_001.fastq.gz
    40
    ```

채널의 각 요소는 이제 `simpleName`과 원래 파일 객체를 포함하는 튜플입니다.

### 4.2. `simplename`에서 메타데이터 추출

이 시점에서 원하는 메타데이터는 `simplename`에 포함되어 있지만 개별 항목에 직접 접근할 수 없습니다.
따라서 `simplename`을 구성 요소로 분할해야 합니다.
다행히도 이러한 구성 요소는 원래 파일명에서 단순히 밑줄로 구분되어 있으므로 이 작업에 완벽한 `tokenize()`라는 일반적인 Nextflow 메서드를 적용할 수 있습니다.

워크플로우를 다음과 같이 편집하세요:

=== "후"

    ```groovy title="main.nf" linenums="7" hl_lines="4"
        // Load files with channel.fromPath
        ch_files = channel.fromPath('data/patientA_rep1_normal_R*_001.fastq.gz')
        ch_files.map { myFile ->
            [ myFile.simpleName.tokenize('_'), myFile ]
        }
    ```

=== "전"

    ```groovy title="main.nf" linenums="7" hl_lines="4"
        // Load files with channel.fromPath
        ch_files = channel.fromPath('data/patientA_rep1_normal_R*_001.fastq.gz')
        ch_files.map { myFile ->
            [ myFile.simpleName, myFile ]
        }
    ```

`tokenize()` 메서드는 밑줄을 찾을 때마다 `simpleName` 문자열을 분할하고 부분 문자열을 포함하는 리스트를 반환합니다.

워크플로우를 실행하세요:

```bash
nextflow run main.nf
```

??? success "명령 출력"

    ```console hl_lines="7-8"
    N E X T F L O W   ~  version 25.10.2

    Launching `main.nf` [gigantic_gauss] DSL2 - revision: a39baabb57

    executor >  local (2)
    [e7/da2f4b] COUNT_LINES (2) [100%] 2 of 2 ✔
    [[patientA, rep1, normal, R2, 001], /workspaces/training/side-quests/working_with_files/data/patientA_rep1_normal_R2_001.fastq.gz]
    [[patientA, rep1, normal, R1, 001], /workspaces/training/side-quests/working_with_files/data/patientA_rep1_normal_R1_001.fastq.gz]
    Processing file: patientA_rep1_normal_R2_001.fastq.gz
    40

    Processing file: patientA_rep1_normal_R1_001.fastq.gz
    40
    ```

이제 채널의 각 요소에 대한 튜플에는 메타데이터 리스트(_예:_ `[patientA, rep1, normal, R1, 001]`)와 원래 파일 객체가 포함됩니다.

훌륭합니다!
환자 정보를 단일 문자열에서 문자열 리스트로 분해했습니다.
이제 환자 정보의 각 부분을 별도로 처리할 수 있습니다.

### 4.3. 맵을 사용하여 메타데이터 구성

현재 메타데이터는 단순한 평면 리스트입니다.
사용하기는 쉽지만 읽기는 어렵습니다.

```console
[patientA, rep1, normal, R1, 001]
```

인덱스 3의 항목은 무엇인가요? 메타데이터 구조에 대한 원래 설명을 다시 참조하지 않고 알 수 있나요?

이것은 키-값 저장소를 사용할 좋은 기회입니다. 여기서 모든 항목에는 키 세트와 관련 값이 있으므로 각 키를 쉽게 참조하여 해당 값을 얻을 수 있습니다.

예제에서 이것은 다음 구성에서:

```groovy
data = [patientA, 1, normal, R1]

println data[3]
```

다음 구성으로 이동하는 것을 의미합니다:

```groovy
data = [id: patientA, replicate: 1, type: normal, readNum: 1]

println data.readNum
```

Nextflow에서 이것을 [맵](https://nextflow.io/docs/latest/script.html#maps)이라고 합니다.

이제 평면 리스트를 맵으로 변환해 봅시다.
워크플로우를 다음과 같이 편집하세요:

=== "후"

    ```groovy title="main.nf" linenums="7" hl_lines="4-13"
        // Load files with channel.fromPath
        ch_files = channel.fromPath('data/patientA_rep1_normal_R*_001.fastq.gz')
        ch_files.map { myFile ->
            def (patient, replicate, type, readNum) = myFile.simpleName.tokenize('_')
            [
              [
                id: patient,
                replicate: replicate.replace('rep', ''),
                type: type,
                readNum: readNum.replace('R', ''),
              ],
              myFile
            ]
        }
    ```

=== "전"

    ```groovy title="main.nf" linenums="7" hl_lines="4"
        // Load files with channel.fromPath
        ch_files = channel.fromPath('data/patientA_rep1_normal_R*_001.fastq.gz')
        ch_files.map { myFile ->
            [ myFile.simpleName.tokenize('_'), myFile ]
        }
    ```

여기서 주요 변경 사항은 다음과 같습니다:

- **구조 분해 할당**: `def (patient, replicate, type, readNum) = ...`는 토큰화된 값을 한 줄에 명명된 변수로 추출합니다
- **맵 리터럴 구문**: `[id: patient, replicate: ...]`는 각 키(예: `id`)가 값(예: `patient`)과 연결된 맵을 생성합니다
- **중첩 구조**: 외부 리스트 `[..., myFile]`는 메타데이터 맵을 원래 파일 객체와 쌍으로 만듭니다

또한 `replace()`라는 문자열 대체 메서드를 사용하여 불필요한 일부 문자를 제거하여 메타데이터 문자열 몇 개를 단순화했습니다(_예:_ 반복 ID에서 숫자만 유지하기 위해 `replicate.replace('rep', '')`).

워크플로우를 다시 실행해 봅시다:

```bash
nextflow run main.nf
```

??? success "명령 출력"

    ```console hl_lines="7-8"
    N E X T F L O W   ~  version 25.10.2

    Launching `main.nf` [infallible_swartz] DSL2 - revision: 7f4e68c0cb

    executor >  local (2)
    [1b/e7fb27] COUNT_LINES (1) [100%] 2 of 2 ✔
    [[id:patientA, replicate:1, type:normal, readNum:2], /workspaces/training/side-quests/working_with_files/data/patientA_rep1_normal_R2_001.fastq.gz]
    [[id:patientA, replicate:1, type:normal, readNum:1], /workspaces/training/side-quests/working_with_files/data/patientA_rep1_normal_R1_001.fastq.gz]
    Processing file: patientA_rep1_normal_R2_001.fastq.gz
    40

    Processing file: patientA_rep1_normal_R1_001.fastq.gz
    40
    ```

이제 메타데이터가 깔끔하게 레이블이 지정되어(_예:_ `[id:patientA, replicate:1, type:normal, readNum:2]`) 무엇이 무엇인지 훨씬 쉽게 알 수 있습니다.

또한 워크플로우에서 메타데이터 요소를 실제로 사용하는 것이 훨씬 쉬워지고 코드를 더 읽기 쉽고 유지 관리하기 쉽게 만들 것입니다.

### 핵심 정리

- Nextflow에서 전체 프로그래밍 언어의 기능으로 파일명을 처리할 수 있습니다
- 파일명을 문자열로 취급하여 관련 정보를 추출할 수 있습니다
- `tokenize()` 및 `replace()`와 같은 메서드를 사용하면 파일명의 문자열을 조작할 수 있습니다
- `.map()` 작업은 구조를 유지하면서 채널 요소를 변환합니다
- 구조화된 메타데이터(맵)는 위치 리스트보다 코드를 더 읽기 쉽고 유지 관리하기 쉽게 만듭니다

다음으로 페어드 데이터 파일을 처리하는 방법을 살펴보겠습니다.

---

## 5. 페어드 데이터 파일 처리

많은 실험 설계는 명시적으로 페어드 방식으로 처리하는 것이 유리한 페어드 데이터 파일을 생성합니다.
예를 들어, 생물정보학에서 시퀀싱 데이터는 종종 페어드 리드 형태로 생성됩니다. 즉, 동일한 DNA 조각에서 유래한 시퀀스 문자열(종종 반대쪽 끝에서 읽히기 때문에 '정방향' 및 '역방향'이라고 함)을 의미합니다.

이것이 예제 데이터의 경우이며, 여기서 R1과 R2는 두 세트의 리드를 나타냅니다.

```console
data/patientA_rep1_normal_R1_001.fastq.gz
data/patientA_rep1_normal_R2_001.fastq.gz
```

Nextflow는 이와 같은 페어드 파일 작업을 위한 특수 채널 팩토리인 `channel.fromFilePairs()`를 제공하며, 이는 공유 명명 패턴을 기반으로 파일을 자동으로 그룹화합니다. 이를 통해 더 적은 노력으로 페어드 파일을 더 긴밀하게 연결할 수 있습니다.

이를 활용하도록 워크플로우를 수정할 것입니다.
두 단계가 필요합니다:

1. 채널 팩토리를 `channel.fromFilePairs()`로 전환
2. 메타데이터 추출 및 매핑

### 5.1. 채널 팩토리를 `channel.fromFilePairs()`로 전환

`channel.fromFilePairs`를 사용하려면 Nextflow가 쌍의 두 멤버를 식별하는 데 사용해야 하는 패턴을 지정해야 합니다.

예제 데이터로 돌아가서 명명 패턴을 다음과 같이 공식화할 수 있습니다:

```console
data/patientA_rep1_normal_R{1,2}_001.fastq.gz
```

이것은 이전에 사용한 glob 패턴과 유사하지만 쌍의 두 멤버를 식별하는 부분 문자열(R 바로 뒤에 오는 `1` 또는 `2`)을 구체적으로 열거합니다.

워크플로우 `main.nf`를 그에 따라 업데이트해 봅시다:

=== "후"

    ```groovy title="main.nf" linenums="7" hl_lines="1-2"
        // Load files with channel.fromFilePairs
        ch_files = channel.fromFilePairs('data/patientA_rep1_normal_R{1,2}_001.fastq.gz')
        /* 지금은 매핑을 주석 처리하고 나중에 다시 돌아오겠습니다!
        ch_files.map { myFile ->
            def (sample, replicate, type, readNum) = myFile.simpleName.tokenize('_')
            [
                [
                    id: sample,
                    replicate: replicate.replace('rep', ''),
                    type: type,
                    readNum: readNum,
                ],
                myFile
            ]
        }
        */
        .view()
    ```

=== "전"

    ```groovy title="main.nf" linenums="7" hl_lines="1-2"
        // Load files with channel.fromFilePairs
        ch_files = channel.fromPath('data/patientA_rep1_normal_R*_001.fastq.gz')
        ch_files.map { myFile ->
            def (sample, replicate, type, readNum) = myFile.simpleName.tokenize('_')
            [
                [
                    id: sample,
                    replicate: replicate.replace('rep', ''),
                    type: type,
                    readNum: readNum,
                ],
                myFile
            ]
        }
        .view()
    ```

채널 팩토리를 전환하고 파일 매칭 패턴을 조정했으며, 그 과정에서 맵 작업을 주석 처리했습니다.
나중에 몇 가지 수정과 함께 다시 추가하겠습니다.

워크플로우를 실행하여 테스트하세요:

```bash
nextflow run main.nf
```

??? failure "명령 출력"

    ```console hl_lines="7-8"
     N E X T F L O W   ~  version 25.10.2

    Launching `main.nf` [angry_koch] DSL2 - revision: 44fdf66105

    [-        ] COUNT_LINES -
    [-        ] COUNT_LINES -
    [patientA_rep1_normal_R, [/workspaces/training/side-quests/working_with_files/data/patientA_rep1_normal_R1_001.fastq.gz, /workspaces/training/side-quests/working_with_files/data/patientA_rep1_normal_R2_001.fastq.gz]]
    ERROR ~ Error executing process > 'COUNT_LINES (1)'

    Caused by:
      Not a valid path value: 'patientA_rep1_normal_R'



    Tip: when you have fixed the problem you can continue the execution adding the option `-resume` to the run command line

    -- Check '.nextflow.log' file for details
    ```

어머, 이번에는 실행이 실패했습니다!

오류 메시지의 관련 부분은 다음과 같습니다:

```console
Not a valid path value: 'patientA_rep1_normal_R'
```

채널 팩토리를 변경했기 때문입니다.
지금까지 원래 입력 채널에는 파일 경로만 포함되어 있었습니다.
우리가 수행한 모든 메타데이터 조작은 실제로 채널 내용에 영향을 미치지 않았습니다.

이제 `.fromFilePairs` 채널 팩토리를 사용하고 있으므로 결과 채널의 내용이 다릅니다.
두 항목을 포함하는 튜플로 구성된 하나의 채널 요소만 볼 수 있습니다: 두 파일이 공유하는 `simpleName`의 일부(식별자 역할을 함)와 두 파일 객체를 포함하는 튜플, `id, [ file1, file2 ]` 형식입니다.

Nextflow가 공유 접두사를 검사하고 환자 식별자로 사용하여 환자 이름을 추출하는 어려운 작업을 수행했기 때문에 훌륭합니다.

그러나 현재 워크플로우가 중단됩니다.
프로세스를 변경하지 않고 동일한 방식으로 `COUNT_LINES`를 계속 실행하려면 파일 경로를 추출하기 위해 매핑 작업을 적용해야 합니다.
하지만 우리는 그렇게 하지 않을 것입니다. 왜냐하면 우리의 궁극적인 목표는 파일 쌍을 적절하게 처리하는 다른 프로세스인 `ANALYZE_READS`를 사용하는 것이기 때문입니다.

따라서 `COUNT_LINES` 호출을 주석 처리(또는 삭제)하고 계속 진행하겠습니다.

=== "후"

    ```groovy title="main.nf" linenums="26" hl_lines="2"
        // Count the lines in the file
        // COUNT_LINES(ch_files)
    ```

=== "전"

    ```groovy title="main.nf" linenums="26" hl_lines="2"
        // Count the lines in the file
        COUNT_LINES(ch_files)
    ```

`COUNT_LINES` include 문도 주석 처리하거나 삭제할 수 있지만 기능적 효과는 없습니다.

이제 워크플로우를 다시 실행해 봅시다:

```bash
nextflow run main.nf
```

??? success "명령 출력"

    ```console hl_lines="5"
     N E X T F L O W   ~  version 25.10.2

    Launching `main.nf` [fabulous_davinci] DSL2 - revision: 22b53268dc

    [patientA_rep1_normal_R, [/workspaces/training/side-quests/working_with_files/data/patientA_rep1_normal_R1_001.fastq.gz, /workspaces/training/side-quests/working_with_files/data/patientA_rep1_normal_R2_001.fastq.gz]]
    ```

좋습니다, 이번에는 워크플로우가 성공했습니다!

그러나 여전히 `id` 필드에서 나머지 메타데이터를 가져와야 합니다.

### 5.2. 파일 쌍에서 메타데이터 추출 및 구성

이전의 `map` 작업은 데이터 구조와 일치하지 않기 때문에 작동하지 않지만 작동하도록 수정할 수 있습니다.

`fromFilePairs()`가 식별자로 사용한 문자열에서 실제 환자 식별자에 이미 접근할 수 있으므로 이전처럼 Path 객체에서 `simpleName`을 가져오지 않고도 메타데이터를 추출하는 데 사용할 수 있습니다.

워크플로우에서 맵 작업의 주석을 해제하고 다음과 같이 편집하세요:

=== "후"

    ```groovy title="main.nf" linenums="7" hl_lines="3-4 9 11 13"
        // Load files with channel.fromFilePairs
        ch_files = channel.fromFilePairs('data/patientA_rep1_normal_R{1,2}_001.fastq.gz')
        ch_files.map { id, files ->
            def (sample, replicate, type) = id.tokenize('_')
            [
                [
                    id: sample,
                    replicate: replicate.replace('rep', ''),
                    type: type
                ],
                files
            ]
        }
        .view()
    ```

=== "전"

    ```groovy title="main.nf" linenums="7" hl_lines="3-5 11 13"
        // Load files with channel.fromFilePairs
        ch_files = channel.fromFilePairs('data/patientA_rep1_normal_R{1,2}_001.fastq.gz')
        /* 지금은 매핑을 주석 처리하고 나중에 다시 돌아오겠습니다!
        ch_files.map { myFile ->
            def (sample, replicate, type, readNum) = myFile.simpleName.tokenize('_')
            [
                [
                    id: sample,
                    replicate: replicate.replace('rep', ''),
                    type: type,
                    readNum: readNum,
                ],
                myFile
            ]
        }
        */
        .view()
    ```

이번에는 맵이 단순히 `myFile` 대신 `id, files`에서 시작하고 `tokenize()`가 `myFile.simpleName` 대신 `id`에 적용됩니다.

또한 `tokenize()` 줄에서 `readNum`을 삭제했습니다. 구체적으로 명명하지 않은 부분 문자열(왼쪽에서 시작)은 자동으로 삭제됩니다.
페어드 파일이 이제 긴밀하게 연결되어 있으므로 메타데이터 맵에 더 이상 `readNum`이 필요하지 않기 때문에 이렇게 할 수 있습니다.

워크플로우를 실행해 봅시다:

```bash
nextflow run main.nf
```

??? success "명령 출력"

    ```console

    N E X T F L O W   ~  version 25.10.2

    Launching `main.nf` [prickly_stonebraker] DSL2 - revision: f62ab10a3f

    [[id:patientA, replicate:1, type:normal], [/workspaces/training/side-quests/working_with_files/data/patientA_rep1_normal_R1_001.fastq.gz, /workspaces/training/side-quests/working_with_files/data/patientA_rep1_normal_R2_001.fastq.gz]]
    ```

그리고 여기 있습니다: 출력 튜플의 첫 번째 위치에 메타데이터 맵(`[id:patientA, replicate:1, type:normal]`)이 있고 그 뒤에 의도한 대로 페어드 파일의 튜플이 있습니다.

물론 이것은 특정 파일 쌍만 선택하고 처리합니다.
여러 쌍을 처리하는 것을 실험하고 싶다면 입력 패턴에 와일드카드를 추가하고 어떤 일이 발생하는지 확인할 수 있습니다.
예를 들어 `data/patientA_rep1_*_R{1,2}_001.fastq.gz`를 사용해 보세요.

### 핵심 정리

- [`channel.fromFilePairs()`는 관련 파일을 자동으로 찾아 쌍으로 만듭니다](https://www.nextflow.io/docs/latest/reference/channel.html#fromfilepairs)
- 이것은 파이프라인에서 페어드 엔드 리드 처리를 단순화합니다
- 페어드 파일은 `[id, [file1, file2]]` 튜플로 그룹화할 수 있습니다
- 메타데이터 추출은 개별 파일이 아닌 페어드 파일 ID에서 수행할 수 있습니다

---

## 6. 프로세스에서 파일 작업 사용하기

이제 이 모든 것을 간단한 프로세스에 통합하여 Nextflow 프로세스 내에서 파일 작업을 사용하는 방법을 강화해 봅시다.

메타데이터 튜플과 입력 파일 쌍을 받아 분석하는 `ANALYZE_READS`라는 미리 작성된 프로세스 모듈을 제공합니다.
이것이 시퀀스 정렬, 변이 호출 또는 이 데이터 유형에 적합한 다른 단계를 수행한다고 상상할 수 있습니다.

시작해 봅시다.

### 6.1. 프로세스 가져오기 및 코드 검토

워크플로우에서 이 프로세스를 사용하려면 워크플로우 블록 앞에 모듈 include 문을 추가하기만 하면 됩니다.

워크플로우를 다음과 같이 편집하세요:

=== "후"

    ```groovy title="main.nf" linenums="1" hl_lines="3"
    #!/usr/bin/env nextflow

    include { ANALYZE_READS } from './modules/analyze_reads.nf'

    workflow {
    ```

=== "전"

    ```groovy title="main.nf" linenums="1"
    #!/usr/bin/env nextflow

    workflow {
    ```

모듈 파일을 열어 코드를 검토할 수 있습니다:

```groovy title="modules/analyze_reads.nf - process example" linenums="1"
#!/usr/bin/env nextflow

process ANALYZE_READS {
    tag { meta.id }

    publishDir { "results/${meta.id}" }, mode: 'copy'

    input:
    tuple val(meta), path(files)

    output:
    tuple val(meta.id), path("${meta.id}_stats.txt")

    script:
    """
    echo "Sample metadata: ${meta.id}" > ${meta.id}_stats.txt
    echo "Replicate: ${meta.replicate}" >> ${meta.id}_stats.txt
    echo "Type: ${meta.type}" >> ${meta.id}_stats.txt
    echo "Read 1: ${files[0]}" >> ${meta.id}_stats.txt
    echo "Read 2: ${files[1]}" >> ${meta.id}_stats.txt
    echo "Read 1 size: \$(gunzip -dc ${files[0]} | wc -l | awk '{print \$1/4}') reads" >> ${meta.id}_stats.txt
    echo "Read 2 size: \$(gunzip -dc ${files[1]} | wc -l | awk '{print \$1/4}') reads" >> ${meta.id}_stats.txt
    """
}
```

!!! note

    `tag` 및 `publishDir` 지시문은 문자열 보간(`"${...}"`) 대신 클로저 구문(`{ ... }`)을 사용합니다.
    이는 이러한 지시문이 런타임까지 사용할 수 없는 입력 변수(`meta`)를 참조하기 때문입니다.
    클로저 구문은 프로세스가 실제로 실행될 때까지 평가를 연기합니다.

!!! note

    관례에 따라 메타데이터 맵을 `meta`라고 부릅니다.
    메타 맵에 대한 심층 분석은 [Metadata and meta maps](./metadata.md) 사이드 퀘스트를 참조하세요.

### 6.2. 워크플로우에서 프로세스 호출

이제 프로세스를 워크플로우에서 사용할 수 있으므로 실행하기 위해 `ANALYZE_READS` 프로세스에 대한 호출을 추가할 수 있습니다.

예제 데이터에서 실행하려면 두 가지 작업을 수행해야 합니다:

1. 재매핑된 채널에 이름 지정
2. 프로세스 호출 추가

#### 6.2.1. 재매핑된 입력 채널 이름 지정

이전에 입력 채널에 직접 매핑 조작을 적용했습니다.
재매핑된 내용을 `ANALYZE_READS` 프로세스에 전달하기 위해(그리고 명확하고 읽기 쉬운 방식으로 수행하기 위해) `ch_samples`라는 새 채널을 생성하려고 합니다.

[`set`](https://www.nextflow.io/docs/latest/reference/operator.html#set) 연산자를 사용하여 이를 수행할 수 있습니다.

메인 워크플로우에서 `.view()` 연산자를 `.set { ch_samples }`로 바꾸고 이름으로 채널을 참조할 수 있는지 테스트하는 줄을 추가하세요.

=== "후"

    ```groovy title="main.nf" linenums="7" hl_lines="14 16-17"
        // Load files with channel.fromFilePairs
        ch_files = channel.fromFilePairs('data/patientA_rep1_normal_R{1,2}_001.fastq.gz')
        ch_files.map { id,  files ->
           def (sample, replicate, type, readNum) = id.tokenize('_')
           [
               [
                   id: sample,
                   replicate: replicate.replace('rep', ''),
                   type: type
               ],
               files
           ]
        }
            .set { ch_samples }

        // Temporary: peek into ch_samples
        ch_samples.view()
    ```

=== "전"

    ```groovy title="main.nf" linenums="7" hl_lines="14"
        // Load files with channel.fromFilePairs
        ch_files = channel.fromFilePairs('data/patientA_rep1_normal_R{1,2}_001.fastq.gz')
        ch_files.map { id,  files ->
           def (sample, replicate, type, readNum) = id.tokenize('_')
           [
               [
                   id: sample,
                   replicate: replicate.replace('rep', ''),
                   type: type
               ],
               files
           ]
        }
        .view()
    }
    ```

이것을 실행해 봅시다:

```bash
nextflow run main.nf
```

??? success "명령 출력"

    ```console
     N E X T F L O W   ~  version 25.10.2

    Launching `main.nf` [goofy_kirch] DSL2 - revision: 3313283e42

    [[id:patientA, replicate:1, type:normal], [/workspaces/training/side-quests/working_with_files/data/patientA_rep1_normal_R1_001.fastq.gz, /workspaces/training/side-quests/working_with_files/data/patientA_rep1_normal_R2_001.fastq.gz]]
    ```

이제 이름으로 채널을 참조할 수 있음을 확인합니다.

#### 6.2.2. 데이터에서 프로세스 호출

이제 `ch_samples` 채널에서 `ANALYZE_READS` 프로세스를 실제로 호출해 봅시다.

메인 워크플로우에서 다음 코드 변경을 수행하세요:

=== "후"

    ```groovy title="main.nf" linenums="23"
        // Run the analysis
        ANALYZE_READS(ch_samples)
    ```

=== "전"

    ```groovy title="main.nf" linenums="23"
        // Temporary: peek into ch_samples
        ch_samples.view()
    ```

이것을 실행해 봅시다:

```bash
nextflow run main.nf
```

??? success "명령 출력"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `./main.nf` [shrivelled_cori] DSL2 - revision: b546a31769

    executor >  local (1)
    [b5/110360] process > ANALYZE_READS (patientA) [100%] 1 of 1 ✔
    ```

이 프로세스는 출력을 `results` 디렉토리에 게시하도록 설정되어 있으므로 거기를 살펴보세요.

??? abstract "디렉토리 및 파일 내용"

    ```console
    results
    └── patientA
        └── patientA_stats.txt
    ```

    ```txt title="patientA_stats.txt"
    Sample metadata: patientA
    Replicate: 1
    Type: normal
    Read 1: patientA_rep1_normal_R1_001.fastq.gz
    Read 2: patientA_rep1_normal_R2_001.fastq.gz
    Read 1 size: 10 reads
    Read 2 size: 10 reads
    ```

프로세스는 입력을 받아 설계된 대로 환자 메타데이터를 포함하는 새 파일을 생성했습니다.
훌륭합니다!

### 6.3. 더 많은 환자 포함

물론 이것은 단일 환자에 대한 단일 파일 쌍만 처리하는 것으로, Nextflow로 얻고자 하는 높은 처리량과는 거리가 멉니다.
아마도 한 번에 훨씬 더 많은 데이터를 처리하고 싶을 것입니다.

`channel.fromPath()`는 입력으로 _glob_을 허용하므로 패턴과 일치하는 파일을 원하는 만큼 받을 수 있습니다.
따라서 모든 환자를 포함하려면 더 많은 환자를 포함하도록 입력 문자열을 수정하기만 하면 됩니다. 앞서 간단히 언급한 대로입니다.

가능한 한 욕심을 부리고 싶다고 가정해 봅시다.
워크플로우를 다음과 같이 편집하세요:

=== "후"

    ```groovy title="main.nf" linenums="7" hl_lines="2"
        // Load files with channel.fromFilePairs
        ch_files = channel.fromFilePairs('data/*_R{1,2}_001.fastq.gz')
    ```

=== "전"

    ```groovy title="main.nf" linenums="7" hl_lines="2"
        // Load files with channel.fromFilePairs
        ch_files = channel.fromFilePairs('data/patientA_rep1_normal_R{1,2}_001.fastq.gz')
    ```

파이프라인을 다시 실행하세요:

```bash
nextflow run main.nf
```

??? success "명령 출력"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `./main.nf` [big_stonebraker] DSL2 - revision: f7f9b8a76c

    executor >  local (8)
    [d5/441891] process > ANALYZE_READS (patientC) [100%] 8 of 8 ✔
    ```

결과 디렉토리에는 이제 사용 가능한 모든 데이터에 대한 결과가 포함되어야 합니다.

??? abstract "디렉토리 내용"

    ```console
    results
    ├── patientA
    │   └── patientA_stats.txt
    ├── patientB
    │   └── patientB_stats.txt
    └── patientC
        └── patientC_stats.txt
    ```

성공! 한 번에 모든 환자를 분석했습니다! 맞나요?

아닐 수도 있습니다.
더 자세히 살펴보면 문제가 있습니다: patientA에 대한 두 개의 반복이 있지만 출력 파일은 하나뿐입니다!
매번 출력 파일을 덮어쓰고 있습니다.

### 6.4. 게시된 파일을 고유하게 만들기

환자 메타데이터에 접근할 수 있으므로 디렉토리 구조나 파일명 자체에 차별화 메타데이터를 포함하여 게시된 파일을 고유하게 만들 수 있습니다.

워크플로우를 다음과 같이 변경하세요:

=== "후"

    ```groovy title="modules/analyze_reads.nf" linenums="6"
        publishDir { "results/${meta.type}/${meta.id}/${meta.replicate}" }, mode: 'copy'
    ```

=== "전"

    ```groovy title="modules/analyze_reads.nf" linenums="6"
        publishDir { "results/${meta.id}" }, mode: 'copy'
    ```

여기서는 샘플 유형과 반복을 고려하기 위해 추가 디렉토리 수준을 사용하는 옵션을 보여주지만 파일명 수준에서도 실험해 볼 수 있습니다.

이제 파이프라인을 한 번 더 실행하되, 깨끗한 작업 공간을 제공하기 위해 먼저 결과 디렉토리를 제거하세요:

```bash
rm -r results
nextflow run main.nf
```

??? success "명령 출력"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `./main.nf` [insane_swartz] DSL2 - revision: fff18abe6d

    executor >  local (8)
    [e3/449081] process > ANALYZE_READS (patientC) [100%] 8 of 8 ✔
    ```

이제 결과 디렉토리를 확인하세요:

??? abstract "디렉토리 내용"

    ```console
    results/
    ├── normal
    │   ├── patientA
    │   │   ├── 1
    │   │   │   └── patientA_stats.txt
    │   │   └── 2
    │   │       └── patientA_stats.txt
    │   ├── patientB
    │   │   └── 1
    │   │       └── patientB_stats.txt
    │   └── patientC
    │       └── 1
    │           └── patientC_stats.txt
    └── tumor
        ├── patientA
        │   ├── 1
        │   │   └── patientA_stats.txt
        │   └── 2
        │       └── patientA_stats.txt
        ├── patientB
        │   └── 1
        │       └── patientB_stats.txt
        └── patientC
            └── 1
                └── patientC_stats.txt
    ```

그리고 여기 있습니다, 모든 메타데이터가 깔끔하게 구성되어 있습니다. 성공입니다!

이와 같이 맵에 메타데이터를 로드하면 훨씬 더 많은 작업을 수행할 수 있습니다:

1. 환자 속성을 기반으로 구성된 출력 디렉토리 생성
2. 환자 속성을 기반으로 프로세스에서 결정 내리기
3. 메타데이터 값을 기반으로 데이터 분할, 결합 및 재결합

메타데이터를 명시적으로 유지하고 데이터에 첨부하는 이 패턴(파일명에 인코딩하는 대신)은 강력하고 유지 관리 가능한 분석 워크플로우를 구축할 수 있게 하는 Nextflow의 핵심 모범 사례입니다.
이에 대한 자세한 내용은 [Metadata and meta maps](./metadata.md) 사이드 퀘스트에서 확인할 수 있습니다.

### 핵심 정리

- `publishDir` 지시문은 메타데이터 값을 기반으로 출력을 구성할 수 있습니다
- 튜플의 메타데이터는 결과의 구조화된 구성을 가능하게 합니다
- 이 접근 방식은 명확한 데이터 출처를 가진 유지 관리 가능한 워크플로우를 생성합니다
- 프로세스는 메타데이터와 파일의 튜플을 입력으로 받을 수 있습니다
- `tag` 지시문은 실행 로그에서 프로세스 식별을 제공합니다
- 워크플로우 구조는 채널 생성과 프로세스 실행을 분리합니다

---

## 요약

이 사이드 퀘스트에서는 기본 작업부터 파일 컬렉션을 처리하는 고급 기법까지 Nextflow에서 파일을 다루는 방법을 배웠습니다.

자신의 작업에 이러한 기법을 적용하면 특히 복잡한 명명 규칙을 가진 대량의 파일을 다룰 때 더 효율적이고 유지 관리 가능한 워크플로우를 구축할 수 있습니다.

### 주요 패턴

1.  **기본 파일 작업:** `file()`로 Path 객체를 생성하고 이름, 확장자, 상위 디렉토리와 같은 파일 속성에 접근하여 문자열과 Path 객체의 차이를 학습했습니다.

    - `file()`로 Path 객체 생성

    ```groovy
    myFile = file('path/to/file.txt')
    ```

    - 파일 속성 가져오기

    ```groovy
    println myFile.name       // file.txt
    println myFile.baseName   // file
    println myFile.extension  // txt
    println myFile.parent     // path/to
    ```

2.  **원격 파일 사용**: URI를 사용하여 로컬 및 원격 파일 간에 투명하게 전환하는 방법을 배웠으며, 워크플로우 로직을 변경하지 않고도 다양한 소스의 파일을 처리하는 Nextflow의 능력을 보여주었습니다.

    - 로컬 파일

    ```groovy
    myFile = file('path/to/file.txt')
    ```

    - FTP

    ```groovy
    myFile = file('ftp://path/to/file.txt')
    ```

    - HTTPS

    ```groovy
    myFile = file('https://path/to/file.txt')
    ```

    - Amazon S3

    ```groovy
    myFile = file('s3://path/to/file.txt')
    ```

    - Azure Blob Storage

    ```groovy
    myFile = file('az://path/to/file.txt')
    ```

    - Google Cloud Storage

    ```groovy
    myFile = file('gs://path/to/file.txt')
    ```

3.  **`fromPath()` 채널 팩토리를 사용하여 파일 로드:** `channel.fromPath()`로 파일 패턴에서 채널을 생성하고 객체 유형을 포함한 파일 속성을 확인했습니다.

    - 파일 패턴에서 채널 생성

    ```groovy
     ch_files = channel.fromPath('data/*.fastq.gz')
    ```

    - 파일 속성 가져오기

    ```groovy
     ch_files.view { myFile ->
        println "File object class: ${myFile.class}"
        println "File name: ${myFile.name}"
        println "Simple name: ${myFile.simpleName}"
        println "Extension: ${myFile.extension}"
        println "Parent directory: ${myFile.parent}"
    }
    ```

4.  **파일명에서 환자 메타데이터 추출:** `tokenize()` 및 `replace()`를 사용하여 파일명에서 메타데이터를 추출하고 구조화하여 구성된 맵으로 변환했습니다.

    ```groovy
    def name = file.name.tokenize('_')
    def patientId = name[0]
    def replicate = name[1].replace('rep', '')
    def type = name[2]
    def readNum = name[3].replace('R', '')
    ```

5.  **channel.fromFilePairs로 단순화:** `channel.fromFilePairs()`를 사용하여 관련 파일을 자동으로 쌍으로 만들고 페어드 파일 ID에서 메타데이터를 추출했습니다.

    ```groovy
    ch_pairs = channel.fromFilePairs('data/*_R{1,2}_001.fastq.gz')
    ```

6.  **프로세스에서 파일 작업 사용:** 적절한 입력 처리를 통해 파일 작업을 Nextflow 프로세스에 통합하고 `publishDir`를 사용하여 메타데이터를 기반으로 출력을 구성했습니다.

    - 프로세스 입력과 메타 맵 연결

    ```groovy
    ch_files = channel.fromFilePairs('data/patientA_rep1_normal_R{1,2}_001.fastq.gz')
    ch_files.map { id,  files ->
        def (sample, replicate, type, readNum) = id.tokenize('_')
        [
            [
                id: sample,
                replicate: replicate.replace('rep', ''),
                type: type
            ],
             files
        ]
    }
        .set { ch_samples }

    ANALYZE_READS(ch_samples)
    ```

    - 메타데이터를 기반으로 출력 구성

    ```groovy
    publishDir { "results/${meta.type}/${meta.id}/${meta.replicate}" }, mode: 'copy'
    ```

### 추가 자료

- [Nextflow 문서: 파일 작업](https://www.nextflow.io/docs/latest/working-with-files.html)
- [channel.fromPath](https://www.nextflow.io/docs/latest/reference/channel.html#frompath)
- [channel.fromFilePairs](https://www.nextflow.io/docs/latest/reference/channel.html#fromfilepairs)

---

## 다음 단계

[사이드 퀘스트 메뉴](./index.md)로 돌아가거나 페이지 오른쪽 하단의 버튼을 클릭하여 목록의 다음 주제로 이동하세요.
