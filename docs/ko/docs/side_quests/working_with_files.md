# 파일 입력 처리

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } AI 지원 번역 - [자세히 알아보기 및 개선 제안](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

과학 분석 워크플로우는 종종 대량의 파일 처리를 포함합니다.
Nextflow는 파일을 효율적으로 처리할 수 있는 강력한 도구를 제공하여, 최소한의 코드로 데이터를 구성하고 처리할 수 있도록 돕습니다.

### 학습 목표

이 사이드 퀘스트에서는 기본 파일 작업부터 파일 컬렉션 작업을 위한 고급 기술까지 Nextflow가 파일을 처리하는 방법을 살펴봅니다.
과학 분석 파이프라인에서 일반적으로 요구되는 파일 이름에서 메타데이터를 추출하는 방법을 배우게 됩니다.

이 사이드 퀘스트를 마치면 다음을 수행할 수 있습니다:

- Nextflow의 `file()` 메서드를 사용하여 파일 경로 문자열에서 Path 객체 생성하기
- 이름, 확장자, 상위 디렉토리와 같은 파일 속성에 접근하기
- URI를 사용하여 로컬 파일과 원격 파일을 투명하게 처리하기
- `channel.fromPath()` 및 `channel.fromFilePairs()`를 사용하여 채널로 파일 처리 자동화하기
- 문자열 조작을 사용하여 파일 이름에서 메타데이터 추출 및 구조화하기
- 패턴 매칭과 glob 표현식을 사용하여 관련 파일 그룹화하기
- 적절한 입력 처리로 파일 작업을 Nextflow process에 통합하기
- 메타데이터 기반 디렉토리 구조를 사용하여 process 출력 구성하기

이러한 기술은 다양한 종류의 파일 입력을 뛰어난 유연성으로 처리할 수 있는 워크플로우를 구축하는 데 도움이 됩니다.

### 사전 요구 사항

이 사이드 퀘스트를 시작하기 전에 다음을 완료해야 합니다:

- [Hello Nextflow](../../hello_nextflow/) 튜토리얼 또는 이에 상응하는 초급 과정을 완료했어야 합니다.
- 기본 Nextflow 개념과 메커니즘(process, channel, operator)을 편안하게 사용할 수 있어야 합니다.

---

## 0. 시작하기

#### 교육 코드스페이스 열기

아직 수행하지 않았다면, [환경 설정](../envsetup/index.md)에 설명된 대로 교육 환경을 열어야 합니다.

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

각 환자에 대해 `tumor`(일반적으로 종양 생검에서 유래) 또는 `normal`(건강한 조직이나 혈액에서 채취) 타입의 샘플이 있습니다.
암 분석에 익숙하지 않다면, 이것이 대조 분석을 수행하기 위해 페어드 tumor/normal 샘플을 사용하는 실험 모델에 해당한다는 것만 알아두시면 됩니다.

특히 환자 A의 경우, 두 세트의 기술적 복제(반복)가 있습니다.

시퀀싱 데이터 파일은 '정방향 리드'와 '역방향 리드'로 알려진 것에 대한 일반적인 `_R1_` 및 `_R2_` 규칙으로 명명되어 있습니다.

_이 실험 설계에 익숙하지 않더라도 걱정하지 마십시오. 이 튜토리얼을 이해하는 데 중요하지 않습니다._

#### 과제 검토

여러분의 과제는 다음을 수행하는 Nextflow 워크플로우를 작성하는 것입니다:

1. Nextflow의 파일 처리 메서드를 사용하여 입력 파일 **로드**
2. 파일 이름 구조에서 메타데이터(환자 ID, 복제, 샘플 타입) **추출**
3. `channel.fromFilePairs()`를 사용하여 페어드 파일(R1/R2) **그룹화**
4. 제공된 분석 모듈로 파일 **처리**
5. 추출된 메타데이터를 기반으로 디렉토리 구조로 출력 **구성**

#### 준비 체크리스트

준비가 되었다고 생각하십니까?

- [ ] 이 과정의 목표와 사전 요구 사항을 이해합니다
- [ ] 코드스페이스가 실행 중입니다
- [ ] 작업 디렉토리를 적절하게 설정했습니다
- [ ] 과제를 이해합니다

모든 항목을 체크할 수 있다면 시작할 준비가 된 것입니다.

---

## 1. 기본 파일 작업

### 1.1. `.class`로 객체의 타입 식별

워크플로우 파일 `main.nf`를 살펴보십시오:

```groovy title="main.nf" linenums="1"
#!/usr/bin/env nextflow

workflow {

    // 문자열 경로에서 Path 객체 생성
    myFile = 'data/patientA_rep1_normal_R1_001.fastq.gz'

    println "${myFile} is of class ${myFile.class}"
}
```

이것은 워크플로우에서 단일 파일 경로를 참조한 다음, 해당 클래스와 함께 콘솔에 출력하는 미니 워크플로우(process 없이)입니다.

??? info "`.class`란 무엇입니까?"

    Nextflow에서 `.class`는 우리가 다루고 있는 객체의 타입을 알려줍니다. 이것은 문자열인지, 숫자인지, 파일인지, 아니면 다른 것인지 알아보기 위해 "이것이 무엇인가요?"라고 묻는 것과 같습니다.
    이것은 다음 섹션에서 일반 문자열과 Path 객체의 차이를 설명하는 데 도움이 됩니다.

워크플로우를 실행해 보겠습니다:

```bash
nextflow run main.nf
```

??? success "명령 출력"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `main.nf` [romantic_chandrasekhar] DSL2 - revision: 5a4a89bc3a

    data/patientA_rep1_normal_R1_001.fastq.gz is of class java.lang.String
    ```

보시다시피, Nextflow는 우리가 작성한 대로 문자열 경로를 정확히 출력했습니다.

이것은 단지 텍스트 출력입니다. Nextflow는 아직 이것으로 특별한 작업을 수행하지 않았습니다.
또한 Nextflow에 관한 한 이것은 문자열(`java.lang.String` 클래스)일 뿐이라는 것을 확인했습니다.
이것은 이해가 됩니다. 우리는 아직 Nextflow에게 이것이 파일에 해당한다고 말하지 않았기 때문입니다.

### 1.2. file()로 Path 객체 생성

경로 문자열에서 [Path 객체](https://www.nextflow.io/docs/latest/reference/stdlib-types.html#path)를 생성하여 Nextflow에게 파일을 처리하는 방법을 알려줄 수 있습니다.

워크플로우에서 `file()` 메서드를 사용하여 문자열 경로 `data/patientA_rep1_normal_R1_001.fastq.gz`를 Path 객체로 변환할 수 있으며, 이는 파일 속성 및 작업에 대한 접근을 제공합니다.

`main.nf`를 다음과 같이 편집하여 문자열을 `file()`로 감싸십시오:

=== "수정 후"

    ```groovy title="main.nf" linenums="5" hl_lines="2"
        // 문자열 경로에서 Path 객체 생성
        myFile = file('data/patientA_rep1_normal_R1_001.fastq.gz')

        println "${myFile} is of class ${myFile.class}"
    ```

=== "수정 전"

    ```groovy title="main.nf" linenums="5" hl_lines="2"
        // 문자열 경로에서 Path 객체 생성
        myFile = 'data/patientA_rep1_normal_R1_001.fastq.gz'

        println "${myFile} is of class ${myFile.class}"
    ```

이제 워크플로우를 다시 실행하십시오:

```bash
nextflow run main.nf
```

??? success "명령 출력"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `main.nf` [kickass_coulomb] DSL2 - revision: 5af44b1b59

    /workspaces/training/side-quests/working_with_files/data/patientA_rep1_normal_R1_001.fastq.gz is of class class sun.nio.fs.UnixPath
    ```

이번에는 입력으로 제공한 상대 경로 대신 전체 절대 경로를 볼 수 있습니다.

Nextflow는 우리의 문자열을 Path 객체로 변환하고 시스템의 실제 파일 위치로 확인했습니다.
파일 경로는 이제 `/workspaces/training/side-quests/working_with_files/data/patientA_rep1_normal_R1_001.fastq.gz`와 같이 절대 경로가 됩니다.

또한 Path 객체 클래스가 `sun.nio.fs.UnixPath`임을 주목하십시오. 이것은 Nextflow가 로컬 파일을 표현하는 방식입니다.
나중에 보겠지만, 원격 파일은 다른 클래스 이름(예: HTTP 파일의 경우 `nextflow.file.http.XPath`)을 가지지만, 모두 정확히 동일한 방식으로 작동하며 워크플로우에서 동일하게 사용할 수 있습니다.

!!! tip

    **주요 차이점:**

    - **경로 문자열**: Nextflow가 문자로 처리하는 텍스트일 뿐입니다
    - **Path 객체**: Nextflow가 작업할 수 있는 스마트 파일 참조입니다

    이렇게 생각하십시오: 경로 문자열은 종이에 주소를 쓰는 것과 같고, Path 객체는 해당 위치로 이동하는 방법을 알고 여정에 대한 세부 정보를 알려줄 수 있는 GPS 장치에 주소를 로드한 것과 같습니다.

### 1.3. 파일 속성 접근

이것이 왜 유용할까요? 이제 Nextflow가 `myFile`이 단순한 문자열이 아니라 Path 객체라는 것을 이해하므로, Path 객체의 다양한 속성에 접근할 수 있습니다.

워크플로우를 업데이트하여 내장 파일 속성을 출력해 봅시다:

=== "수정 후"

    ```groovy title="main.nf" linenums="5" hl_lines="4-9"
        // 문자열 경로에서 Path 객체 생성
        myFile = file('data/patientA_rep1_normal_R1_001.fastq.gz')

        // 파일 속성 출력
        println "File object class: ${myFile.class}"
        println "File name: ${myFile.name}"
        println "Simple name: ${myFile.simpleName}"
        println "Extension: ${myFile.extension}"
        println "Parent directory: ${myFile.parent}"
    ```

=== "수정 전"

    ```groovy title="main.nf" linenums="5" hl_lines="4"
        // 문자열 경로에서 Path 객체 생성
        myFile = file('data/patientA_rep1_normal_R1_001.fastq.gz')

        println "${myFile} is of class ${myFile.class}"
    ```

워크플로우를 실행하십시오:

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

위 콘솔에 출력된 다양한 파일 속성을 볼 수 있습니다.

### 1.4. process에 파일 제공

문자열과 Path 객체의 차이는 process를 사용하여 실제 워크플로우를 구축하기 시작할 때 중요해집니다.
지금까지 Nextflow가 입력 파일을 파일로 처리하고 있다는 것을 확인했지만, process에서 실제로 해당 파일에 대해 무언가를 실행할 수 있는지 살펴봅시다.

#### 1.4.1. process 가져오기 및 코드 검토

파일 입력을 받아 포함된 줄 수를 세는 `COUNT_LINES`라는 미리 작성된 process 모듈을 제공합니다.

워크플로우에서 process를 사용하려면 workflow 블록 앞에 include 문을 추가하기만 하면 됩니다:

=== "수정 후"

    ```groovy title="main.nf" linenums="1" hl_lines="3"
    #!/usr/bin/env nextflow

    include { COUNT_LINES } from './modules/count_lines.nf'

    workflow {
    ```

=== "수정 전"

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

보시다시피, 파일의 압축을 풀고 포함된 줄 수를 세는 상당히 간단한 작은 스크립트입니다.

??? info "`debug true`는 무엇을 합니까?"

    process 정의의 `debug true` 지시문은 Nextflow가 스크립트의 출력(줄 수 "40"과 같은)을 실행 로그에 직접 출력하도록 합니다.
    이것이 없으면 process 실행 상태만 볼 수 있고 스크립트의 실제 출력은 볼 수 없습니다.

    Nextflow process 디버깅에 대한 자세한 내용은 [Debugging Nextflow Workflows](debugging.md) 사이드 퀘스트를 참조하십시오.

#### 1.4.2. `COUNT_LINES` 호출 추가

이제 process를 워크플로우에서 사용할 수 있으므로, 입력 파일에서 실행하기 위해 `COUNT_LINES` process에 대한 호출을 추가할 수 있습니다.

워크플로우를 다음과 같이 편집하십시오:

=== "수정 후"

    ```groovy title="main.nf" linenums="7" hl_lines="11-12"
        // 문자열 경로에서 Path 객체 생성
        myFile = file('data/patientA_rep1_normal_R1_001.fastq.gz')

        // 파일 속성 출력
        println "File object class: ${myFile.class}"
        println "File name: ${myFile.name}"
        println "Simple name: ${myFile.simpleName}"
        println "Extension: ${myFile.extension}"
        println "Parent directory: ${myFile.parent}"

        // 파일의 줄 수 계산
        COUNT_LINES(myFile)
    ```

=== "수정 전"

    ```groovy title="main.nf" linenums="7" hl_lines="4-9"
        // 문자열 경로에서 Path 객체 생성
        myFile = file('data/patientA_rep1_normal_R1_001.fastq.gz')

        // 파일 속성 출력
        println "File object class: ${myFile.class}"
        println "File name: ${myFile.name}"
        println "Simple name: ${myFile.simpleName}"
        println "Extension: ${myFile.extension}"
        println "Parent directory: ${myFile.parent}"
    ```

이제 워크플로우를 실행하십시오:

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

이것은 process 내부에서 파일을 적절하게 조작할 수 있음을 보여줍니다.

구체적으로, Nextflow는 다음 작업을 성공적으로 수행했습니다:

- 작업 디렉토리에 파일 스테이징
- .gz 파일 압축 해제
- 줄 수 세기(이 경우 40줄)
- 오류 없이 완료

이 원활한 작업의 핵심은 Nextflow에게 입력이 파일이며 그렇게 처리되어야 함을 명시적으로 알려주는 것입니다.

### 1.5. 기본 파일 입력 오류 문제 해결

이것은 종종 Nextflow 초보자를 당황하게 하므로, 잘못 수행했을 때 어떤 일이 발생하는지 몇 분 동안 살펴보겠습니다.

파일 처리를 잘못할 수 있는 주요 위치가 두 곳 있습니다: 워크플로우 수준과 process 수준입니다.

#### 1.5.1. 워크플로우 수준 오류

workflow 블록에서 입력을 지정할 때 파일을 문자열로 처리하도록 되돌리면 어떤 일이 발생하는지 살펴봅시다.

워크플로우를 다음과 같이 편집하고, 경로 관련 print 문은 주석 처리하십시오:

=== "수정 후"

    ```groovy title="main.nf" linenums="7" hl_lines="2 6-11"
        // 문자열 경로에서 Path 객체 생성
        myFile = 'data/patientA_rep1_normal_R1_001.fastq.gz'

        // 파일 속성 출력
        println "File object class: ${myFile.class}"
        /*
        println "File name: ${myFile.name}"
        println "Simple name: ${myFile.simpleName}"
        println "Extension: ${myFile.extension}"
        println "Parent directory: ${myFile.parent}"
        */

        // 파일의 줄 수 계산
        COUNT_LINES(myFile)
    ```

=== "수정 전"

    ```groovy title="main.nf" linenums="7" hl_lines="4-9"
        // 문자열 경로에서 Path 객체 생성
        myFile = file('data/patientA_rep1_normal_R1_001.fastq.gz')

        // 파일 속성 출력
        println "File object class: ${myFile.class}"
        println "File name: ${myFile.name}"
        println "Simple name: ${myFile.simpleName}"
        println "Extension: ${myFile.extension}"
        println "Parent directory: ${myFile.parent}"

        // 파일의 줄 수 계산
        COUNT_LINES(myFile)
    ```

이제 워크플로우를 실행하십시오:

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

`path` 입력을 지정하면 Nextflow는 단순한 문자열이 아닌 실제 파일 참조를 전달하는지 확인합니다.
이 오류는 `'data/patientA_rep1_normal_R1_001.fastq.gz'`가 Path 객체가 아니라 문자열이기 때문에 유효한 경로 값이 아니라고 알려줍니다.

Nextflow는 즉시 문제를 감지하고 process를 시작하기 전에 중지했습니다.

#### 1.5.2. Process 수준 오류

입력을 파일로 처리하도록 Nextflow에게 지정하는 것을 잊을 수 있는 다른 위치는 process 정의입니다.

!!! warning "1.5.1의 워크플로우 오류 유지"

    이 테스트가 올바르게 작동하려면 워크플로우를 손상된 상태로 유지하십시오(`file()` 대신 일반 문자열 사용).
    process에서 `val`과 결합하면 아래에 표시된 오류가 발생합니다.

모듈을 다음과 같이 편집하십시오:

=== "수정 후"

    ```groovy title="modules/count_lines.nf" linenums="3" hl_lines="5"
    process COUNT_LINES {
        debug true

        input:
        val input_file
    ```

=== "수정 전"

    ```groovy title="modules/count_lines.nf" linenums="3" hl_lines="5"
    process COUNT_LINES {
        debug true

        input:
        path input_file
    ```

이제 워크플로우를 다시 실행하십시오:

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

이것은 위에서 언급한 대로 process가 디버깅 정보를 출력하도록 설정되어 있기 때문에 오류에 대한 많은 세부 정보를 보여줍니다.

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

시스템이 파일을 찾을 수 없다고 말합니다. 그러나 경로를 찾아보면 해당 위치에 그 이름의 파일이 있습니다.

이것을 실행했을 때, Nextflow는 문자열 값을 스크립트로 전달했지만 작업 디렉토리에 실제 파일을 *스테이징*하지 않았습니다.
따라서 process는 상대 문자열 `data/patientA_rep1_normal_R1_001.fastq.gz`를 사용하려고 했지만, 해당 파일은 process 작업 디렉토리 내에 존재하지 않습니다.

이 두 예제를 종합하면, 입력을 파일로 처리해야 한다고 Nextflow에게 알려주는 것이 얼마나 중요한지 알 수 있습니다.

!!! note

    다음 섹션으로 계속하기 전에 두 의도적 오류를 모두 수정하십시오.

### 요점 정리

- 경로 문자열 vs Path 객체: 문자열은 단순한 텍스트이고, Path 객체는 스마트 파일 참조입니다
- `file()` 메서드는 문자열 경로를 Nextflow가 작업할 수 있는 Path 객체로 변환합니다
- `name`, `simpleName`, `extension`, `parent`와 같은 파일 속성에 [파일 속성을 사용하여](https://www.nextflow.io/docs/latest/working-with-files.html#getting-file-attributes) 접근할 수 있습니다
- 문자열 대신 Path 객체를 사용하면 Nextflow가 워크플로우에서 파일을 적절하게 관리할 수 있습니다
- Process 입력 결과: 적절한 파일 처리는 파일이 올바르게 스테이징되고 process에서 사용할 수 있도록 문자열이 아닌 Path 객체가 필요합니다

---

## 2. 원격 파일 사용

Nextflow의 주요 기능 중 하나는 로컬 파일(동일한 머신에 있는)에서 인터넷을 통해 접근할 수 있는 원격 파일로 원활하게 전환할 수 있다는 것입니다.

올바르게 수행하면 다른 위치에서 오는 파일을 수용하기 위해 워크플로우의 로직을 변경할 필요가 없어야 합니다.
원격 파일을 사용하기 위해 해야 할 일은 워크플로우에 파일을 제공할 때 파일 경로에 적절한 접두사를 지정하는 것뿐입니다.

예를 들어, `/path/to/data`에는 접두사가 없어 '일반' 로컬 파일 경로임을 나타내는 반면, `s3://path/to/data`에는 Amazon의 S3 객체 스토리지에 위치함을 나타내는 `s3://` 접두사가 포함되어 있습니다.

다양한 프로토콜이 지원됩니다:

- HTTP(S)/FTP (http://, https://, ftp://)
- Amazon S3 (s3://)
- Azure Blob Storage (az://)
- Google Cloud Storage (gs://)

이들 중 하나를 사용하려면 문자열에 관련 접두사를 지정하기만 하면 되며, 이는 기술적으로 파일 경로 대신 URI(Uniform Resource Identifier)라고 불립니다.
Nextflow는 인증 및 파일을 적절한 위치로 스테이징, 다운로드 또는 업로드 및 예상되는 모든 다른 파일 작업을 처리합니다.

이 시스템의 주요 강점은 환경 간 전환 시 파이프라인 로직을 변경하지 않고도 가능하다는 것입니다.
예를 들어, URI만 변경하여 작고 로컬인 테스트 세트로 개발한 후 원격 스토리지에 있는 전체 규모 테스트 세트로 전환할 수 있습니다.

### 2.1. 인터넷의 파일 사용

워크플로우에 제공하는 로컬 경로를 Github에 저장된 동일한 데이터의 사본을 가리키는 HTTPS 경로로 전환하여 이것을 테스트해 봅시다.

!!! warning

    이것은 활성 인터넷 연결이 있는 경우에만 작동합니다.

`main.nf`를 다시 열고 다음과 같이 입력 경로를 변경하십시오:

=== "수정 후"

    ```groovy title="main.nf" linenums="2" hl_lines="2"
        // 인터넷에서 원격 파일 사용
        myFile = file('https://raw.github.com/nextflow-io/training/master/side-quests/working_with_files/data/patientA_rep1_normal_R1_001.fastq.gz')

        // 파일 속성 출력
        println "File object class: ${myFile.class}"
        println "File name: ${myFile.name}"
        println "Simple name: ${myFile.simpleName}"
        println "Extension: ${myFile.extension}"
        println "Parent directory: ${myFile.parent}"
    ```

=== "수정 전"

    ```groovy title="main.nf" linenums="2" hl_lines="2"
        // 문자열 경로에서 Path 객체 생성
        myFile = file('data/patientA_rep1_normal_R1_001.fastq.gz')

        // 파일 속성 출력
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

작동합니다! 콘솔 출력에서 변경된 것이 거의 없음을 볼 수 있습니다.

콘솔 출력의 한 가지 차이점은 경로 객체 클래스가 이제 `nextflow.file.http.XPath`인 반면 로컬 경로의 경우 클래스는 `sun.nio.fs.UnixPath`였다는 것입니다.
이러한 클래스를 기억할 필요는 없습니다. Nextflow가 다른 위치를 적절하게 식별하고 처리한다는 것을 보여주기 위해 언급하는 것입니다.

백그라운드에서 Nextflow는 work 디렉토리 내에 있는 스테이징 디렉토리로 파일을 다운로드했습니다.
그런 다음 스테이징된 파일은 로컬 파일로 처리되어 관련 process 디렉토리에 심볼릭 링크될 수 있습니다.

process의 해시 값에 있는 작업 디렉토리의 내용을 살펴보면 이것이 발생했음을 확인할 수 있습니다.

??? abstract "작업 디렉토리 내용"

    process 해시가 `8a/2ab7ca`인 경우 작업 디렉토리를 탐색할 수 있습니다:

    ```console
    $ ls -la work/8a/2ab7ca*/
    total 16
    drwxr-xr-x  6 user  staff   192 Jan 28 10:00 .
    drwxr-xr-x  3 user  staff    96 Jan 28 10:00 ..
    -rw-r--r--  1 user  staff     0 Jan 28 10:00 .command.begin
    -rw-r--r--  1 user  staff   127 Jan 28 10:00 .command.sh
    lrwxr-xr-x  1 user  staff    89 Jan 28 10:00 patientA_rep1_normal_R1_001.fastq.gz -> /path/to/work/stage/.../patientA_rep1_normal_R1_001.fastq.gz
    ```

    심볼릭 링크는 Nextflow가 자동으로 다운로드한 원격 파일의 스테이징된 사본을 가리킵니다.

더 큰 파일의 경우 다운로드 단계가 로컬 파일에서 실행하는 것에 비해 추가 시간이 걸립니다.
그러나 Nextflow는 불필요한 다운로드를 피하기 위해 이미 스테이징된 사본이 있는지 확인합니다.
따라서 동일한 파일에서 다시 실행하고 스테이징된 파일을 삭제하지 않았다면 Nextflow는 스테이징된 사본을 사용합니다.

이것은 Nextflow를 사용하여 로컬 및 원격 데이터 간을 전환하는 것이 얼마나 쉬운지 보여주며, 이것이 Nextflow의 주요 기능입니다.

!!! note

    이 원칙에 대한 한 가지 중요한 예외는 HTTPS가 여러 파일을 나열할 수 없기 때문에 HTTPS에서 glob 패턴이나 디렉토리 경로를 사용할 수 없다는 것입니다. 따라서 정확한 파일 URL을 지정해야 합니다.
    그러나 blob 스토리지(`s3://`, `az://`, `gs://`)와 같은 다른 스토리지 프로토콜은 glob과 디렉토리 경로를 모두 사용할 수 있습니다.

    클라우드 스토리지에서 glob 패턴을 사용하는 방법은 다음과 같습니다:

    ```groovy title="클라우드 스토리지 예제 (이 환경에서 실행 불가)"
    // glob 패턴을 사용한 S3 - 여러 파일과 일치
    ch_s3_files = channel.fromPath('s3://my-bucket/data/*.fastq.gz')

    // glob 패턴을 사용한 Azure Blob Storage
    ch_azure_files = channel.fromPath('az://container/data/patient*_R{1,2}.fastq.gz')

    // glob 패턴을 사용한 Google Cloud Storage
    ch_gcs_files = channel.fromPath('gs://bucket/data/sample_*.fastq.gz')
    ```

    다음 섹션에서 glob를 실제로 사용하는 방법을 보여드리겠습니다.

### 2.2. 로컬 파일로 다시 전환

이 사이드 퀘스트의 나머지 부분에서는 로컬 예제 파일을 사용할 것이므로 워크플로우 입력을 원래 파일로 다시 전환하겠습니다:

=== "수정 후"

    ```groovy title="main.nf" linenums="2" hl_lines="2"
        // 문자열 경로에서 Path 객체 생성
        myFile = file('data/patientA_rep1_normal_R1_001.fastq.gz')

        // 파일 속성 출력
        println "File object class: ${myFile.class}"
        println "File name: ${myFile.name}"
        println "Simple name: ${myFile.simpleName}"
        println "Extension: ${myFile.extension}"
        println "Parent directory: ${myFile.parent}"
    ```

=== "수정 전"

    ```groovy title="main.nf" linenums="2" hl_lines="2"
        // 문자열 경로에서 Path 객체 생성
        myFile = file('https://raw.github.com/nextflow-io/training/master/side-quests/working_with_files/data/patientA_rep1_normal_R1_001.fastq.gz')

        // 파일 속성 출력
        println "File object class: ${myFile.class}"
        println "File name: ${myFile.name}"
        println "Simple name: ${myFile.simpleName}"
        println "Extension: ${myFile.extension}"
        println "Parent directory: ${myFile.parent}"
    ```

### 요점 정리

- 원격 데이터는 URI를 사용하여 접근합니다(HTTP, FTP, S3, Azure, Google Cloud)
- 이러한 경로를 process에 제공하는 한 Nextflow는 자동으로 데이터를 다운로드하고 적절한 위치로 스테이징합니다
- 원격 파일을 다운로드하거나 업로드하는 로직을 작성하지 마십시오!
- 로컬 및 원격 파일은 다른 객체 타입을 생성하지만 동일하게 작동합니다
- **중요**: HTTP/HTTPS는 단일 파일에서만 작동합니다(glob 패턴 없음)
- 클라우드 스토리지(S3, Azure, GCS)는 단일 파일과 glob 패턴을 모두 지원합니다
- 프로토콜이 필요한 작업을 지원하는 한 코드 로직을 변경하지 않고 로컬 및 원격 데이터 소스 간을 원활하게 전환할 수 있습니다

---

## 3. `fromPath()` 채널 팩토리 사용

지금까지 우리는 한 번에 하나의 파일로 작업했지만, Nextflow에서는 일반적으로 처리할 여러 입력 파일이 있는 입력 채널을 생성하려고 합니다.

이를 수행하는 단순한 방법은 `file()` 메서드를 다음과 같이 [`channel.of()`](https://www.nextflow.io/docs/latest/reference/channel.html#of)와 결합하는 것입니다:

```groovy title="구문 예제"
ch_files = channel.of([file('data/patientA_rep1_normal_R1_001.fastq.gz')],
                      [file('data/patientA_rep1_normal_R1_001.fastq.gz')])
```

이것은 작동하지만 투박합니다.

!!! tip "`file()` vs `channel.fromPath()`를 사용할 때"

    - 단일 Path 객체를 직접 조작해야 할 때(파일이 존재하는지 확인, 속성 읽기 또는 단일 process 호출로 전달) `file()`을 사용하십시오
    - 여러 파일을 보유할 수 있는 채널이 필요할 때, 특히 glob 패턴을 사용하거나 파일이 여러 process를 통해 흐를 때 `channel.fromPath()`를 사용하십시오

여기서 [`channel.fromPath()`](https://www.nextflow.io/docs/latest/reference/channel.html#frompath)가 등장합니다: 하나 이상의 정적 파일 문자열뿐만 아니라 glob 패턴에서 채널을 생성하는 데 필요한 모든 기능을 번들로 제공하는 편리한 채널 팩토리입니다.

### 3.1. 채널 팩토리 추가

`channel.fromPath`를 사용하도록 워크플로우를 업데이트해 봅시다.

=== "수정 후"

    ```groovy title="main.nf" linenums="7" hl_lines="1-3"
        // channel.fromPath로 파일 로드
        ch_files = channel.fromPath('data/patientA_rep1_normal_R1_001.fastq.gz')
        ch_files.view { myFile -> "Found file: $myFile" }

        // 파일 속성 출력
        /* Comment these out for now, we'll come back to them!
        println "File object class: ${myFile.class}"
        println "File name: ${myFile.name}"
        println "Simple name: ${myFile.simpleName}"
        println "Extension: ${myFile.extension}"
        println "Parent directory: ${myFile.parent}"
        */

        // 파일의 줄 수 계산
        // COUNT_LINES(myFile)
    ```

=== "수정 전"

    ```groovy title="main.nf" linenums="7" hl_lines="1-2"
        // 문자열 경로에서 Path 객체 생성
        myFile = file('data/patientA_rep1_normal_R1_001.fastq.gz')

        // 파일 속성 출력
        println "File object class: ${myFile.class}"
        println "File name: ${myFile.name}"
        println "Simple name: ${myFile.simpleName}"
        println "Extension: ${myFile.extension}"
        println "Parent directory: ${myFile.parent}"

        // 파일의 줄 수 계산
        COUNT_LINES(myFile)
    ```

또한 지금은 속성을 출력하는 코드를 주석 처리했으며 대신 파일 이름만 출력하는 `.view` 문을 추가했습니다.

워크플로우를 실행하십시오:

```bash
nextflow run main.nf
```

??? success "명령 출력"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `main.nf` [grave_meucci] DSL2 - revision: b09964a583

    Found file: /workspaces/training/side-quests/working_with_files/data/patientA_rep1_normal_R1_001.fastq.gz
    ```

보시다시피, 파일 경로가 채널에 `Path` 타입 객체로 로드되고 있습니다.
이것은 `file()`이 했을 것과 유사하지만, 이제 원한다면 더 많은 파일을 로드할 수 있는 채널이 있습니다.

`channel.fromPath()`를 사용하는 것은 파일 목록으로 채워진 새 채널을 생성하는 편리한 방법입니다.

### 3.2. 채널의 파일 속성 보기

채널 팩토리를 사용한 첫 번째 시도에서는 코드를 단순화하고 파일 이름만 출력했습니다.

전체 파일 속성을 출력하는 것으로 돌아가 봅시다:

=== "수정 후"

    ```groovy title="main.nf" linenums="7" hl_lines="3-9 12"
        // channel.fromPath로 파일 로드
        ch_files = channel.fromPath('data/patientA_rep1_normal_R1_001.fastq.gz')
        ch_files.view { myFile ->
            println "File object class: ${myFile.class}"
            println "File name: ${myFile.name}"
            println "Simple name: ${myFile.simpleName}"
            println "Extension: ${myFile.extension}"
            println "Parent directory: ${myFile.parent}"
        }

        // 파일의 줄 수 계산
        COUNT_LINES(ch_files)
    ```

=== "수정 전"

    ```groovy title="main.nf" linenums="7" hl_lines="3"
        // channel.fromPath로 파일 로드
        ch_files = channel.fromPath('data/patientA_rep1_normal_R1_001.fastq.gz')
        ch_files.view { myFile -> "Found file: $myFile" }

        // 파일의 줄 수 계산
        // COUNT_LINES(ch_files)
    ```

또한 채널 기반 접근 방식으로 파일 처리가 여전히 올바르게 작동하는지 확인하기 위해 `COUNT_LINES` process 호출을 다시 활성화하고 있습니다.

워크플로우를 실행하십시오:

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

이제 파일이 채널에 있으므로 이전과 동일한 결과를 얻었으며, 더 많이 추가할 수 있습니다.

### 3.3. glob을 사용하여 여러 파일 매칭

채널에 더 많은 파일을 로드할 수 있는 여러 방법이 있습니다.
여기서는 와일드카드 문자를 기반으로 파일 및 디렉토리 이름을 매칭하고 검색하는 편리한 방법인 glob 패턴을 사용하는 방법을 보여드리겠습니다.
이러한 패턴을 매칭하는 프로세스를 "globbing" 또는 "filename expansion"이라고 합니다.

!!! note

    이전에 언급했듯이, HTTPS가 여러 파일을 나열할 수 없기 때문에 HTTPS 파일 경로를 제외한 대부분의 경우 Nextflow는 입력 및 출력 파일을 관리하기 위한 globbing을 지원합니다.

특정 환자 `patientA`와 관련된 파일 쌍의 두 파일을 모두 검색하고 싶다고 가정해 봅시다:

```console
patientA_rep1_normal_R1_001.fastq.gz
patientA_rep1_normal_R2_001.fastq.gz
```

파일 이름 간의 유일한 차이점은 복제 번호, 즉 `R` 뒤의 숫자이므로 다음과 같이 와일드카드 문자 `*`를 사용하여 숫자를 대신할 수 있습니다:

```console
patientA_rep1_normal_R*_001.fastq.gz
```

이것이 우리에게 필요한 glob 패턴입니다.

이제 채널 팩토리의 파일 경로를 다음과 같이 해당 glob 패턴을 사용하도록 업데이트하기만 하면 됩니다:

=== "수정 후"

    ```groovy title="main.nf" linenums="7"
      // channel.fromPath로 파일 로드
        ch_files = channel.fromPath('data/patientA_rep1_normal_R*_001.fastq.gz')
    ```

=== "수정 전"

    ```groovy title="main.nf" linenums="7"
      // channel.fromPath로 파일 로드
        ch_files = channel.fromPath('data/patientA_rep1_normal_R1_001.fastq.gz')
    ```

Nextflow는 이것이 glob 패턴임을 자동으로 인식하고 적절하게 처리합니다.

워크플로우를 실행하여 테스트하십시오:

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

보시다시피, 이제 채널에 두 개의 Path 객체가 있으며, 이는 Nextflow가 파일 이름 확장을 올바르게 수행했고 예상대로 두 파일을 모두 로드하고 처리했음을 보여줍니다.

이 방법을 사용하면 glob 패턴을 변경하기만 하면 원하는 만큼 많은 파일을 검색할 수 있습니다. 예를 들어 파일 이름의 모든 가변 부분을 `*`로 바꾸어 더 관대하게 만들면(_예:_ `data/patient*_rep*_*_R*_001.fastq.gz`) `data` 디렉토리의 모든 예제 파일을 가져올 수 있습니다.

### 요점 정리

- `channel.fromPath()`는 패턴과 매칭하는 파일로 채널을 생성합니다
- 각 파일은 채널의 별도 요소로 방출됩니다
- glob 패턴을 사용하여 여러 파일을 매칭할 수 있습니다
- 파일은 전체 속성을 가진 Path 객체로 자동 변환됩니다
- `.view()` 메서드를 사용하면 채널 내용을 검사할 수 있습니다

---

## 4. 파일 이름에서 기본 메타데이터 추출

대부분의 과학 분야에서 데이터를 포함하는 파일 이름
