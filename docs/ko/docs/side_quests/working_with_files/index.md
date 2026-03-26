# 파일 입력 처리

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } AI 지원 번역 - [자세히 알아보기 및 개선 제안](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

과학적 분석 워크플로우는 대개 많은 수의 파일을 처리하는 작업을 포함합니다.
Nextflow는 파일을 효율적으로 처리할 수 있는 강력한 도구를 제공하여, 최소한의 코드로 데이터를 구성하고 처리할 수 있도록 도와줍니다.

### 학습 목표

이 사이드 퀘스트에서는 기본적인 파일 작업부터 파일 컬렉션을 다루는 고급 기법까지, Nextflow가 파일을 처리하는 방법을 살펴봅니다.
파일 이름에서 메타데이터를 추출하는 방법도 학습합니다. 이는 과학적 분석 파이프라인에서 자주 요구되는 작업입니다.

이 사이드 퀘스트를 마치면 다음을 수행할 수 있습니다:

- Nextflow의 `file()` 메서드를 사용하여 파일 경로 문자열에서 Path 객체 생성
- 파일 이름, 확장자, 상위 디렉토리 등 파일 속성 접근
- URI를 사용하여 로컬 파일과 원격 파일을 동일하게 처리
- `channel.fromPath()`와 `channel.fromFilePairs()`를 사용하여 파일 처리 자동화
- 문자열 조작을 통해 파일 이름에서 메타데이터 추출 및 구조화
- 패턴 매칭과 glob 표현식을 사용하여 관련 파일 그룹화
- 적절한 입력 처리를 통해 Nextflow 프로세스에 파일 작업 통합
- 메타데이터 기반 디렉토리 구조를 사용하여 프로세스 출력 구성

이러한 기술을 익히면 다양한 종류의 파일 입력을 유연하게 처리하는 워크플로우를 구축할 수 있습니다.

### 사전 요구 사항

이 사이드 퀘스트를 시작하기 전에 다음을 완료해야 합니다:

- [Hello Nextflow](../../hello_nextflow/) 튜토리얼 또는 동급의 입문 과정 완료
- 기본적인 Nextflow 개념과 메커니즘(프로세스, 채널, 연산자)에 익숙할 것

<!-- I removed the suggestion to do the metamaps SQ first because that works more naturally after -->

---

## 0. 시작하기

#### 교육 코드스페이스 열기

아직 열지 않았다면, [환경 설정](../envsetup/index.md)에 설명된 대로 교육 환경을 열어 주세요.

[![Open in GitHub Codespaces](https://github.com/codespaces/badge.svg)](https://codespaces.new/nextflow-io/training?quickstart=1&ref=master)

#### 프로젝트 디렉토리로 이동

이 튜토리얼의 파일이 있는 디렉토리로 이동합니다.

```bash
cd side-quests/working_with_files
```

VSCode에서 이 디렉토리를 포커스로 설정할 수 있습니다:

```bash
code .
```

#### 자료 검토

`main.nf`라는 간단한 워크플로우 파일, 두 개의 모듈 파일이 있는 `modules` 디렉토리, 그리고 예제 데이터 파일이 있는 `data` 디렉토리를 확인할 수 있습니다.

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
암 분석에 익숙하지 않더라도, 이것이 대조 분석을 수행하기 위해 종양/정상 샘플 쌍을 사용하는 실험 모델에 해당한다는 것만 알면 됩니다.

특히 환자 A의 경우, 두 세트의 기술적 반복(리플리케이트)이 있습니다.

시퀀싱 데이터 파일은 'forward reads'와 'reverse reads'로 알려진 것에 대해 일반적인 `_R1_` 및 `_R2_` 규칙으로 명명되어 있습니다.

_이 실험 설계에 익숙하지 않더라도 이 튜토리얼을 이해하는 데 중요하지 않으니 걱정하지 마세요._

#### 과제 검토

여러분의 과제는 다음을 수행하는 Nextflow 워크플로우를 작성하는 것입니다:

1. Nextflow의 파일 처리 메서드를 사용하여 입력 파일 **로드**
2. 파일 이름 구조에서 메타데이터(환자 ID, 리플리케이트, 샘플 유형) **추출**
3. `channel.fromFilePairs()`를 사용하여 페어드 파일(R1/R2) **그룹화**
4. 제공된 분석 모듈로 파일 **처리**
5. 추출된 메타데이터를 기반으로 디렉토리 구조에 출력 **구성**

#### 준비 체크리스트

시작할 준비가 되었나요?

- [ ] 이 과정의 목표와 사전 요구 사항을 이해했습니다
- [ ] 코드스페이스가 실행 중입니다
- [ ] 작업 디렉토리를 적절히 설정했습니다
- [ ] 과제를 이해했습니다

모든 항목을 확인했다면 시작할 준비가 된 것입니다.

---

## 1. 기본 파일 작업

### 1.1. `.class`로 객체 유형 확인

워크플로우 파일 `main.nf`를 살펴봅니다:

```groovy title="main.nf" linenums="1"
#!/usr/bin/env nextflow

workflow {

    // 문자열 경로에서 Path 객체 생성
    myFile = 'data/patientA_rep1_normal_R1_001.fastq.gz'

    println "${myFile} is of class ${myFile.class}"
}
```

이것은 워크플로우에서 단일 파일 경로를 참조한 다음, 해당 클래스와 함께 콘솔에 출력하는 간단한 미니 워크플로우(프로세스 없음)입니다.

??? info "`.class`란 무엇인가요?"

    Nextflow에서 `.class`는 우리가 다루고 있는 객체의 유형을 알려줍니다. "이것은 어떤 종류의 것인가?"라고 묻는 것과 같습니다. 문자열인지, 숫자인지, 파일인지, 아니면 다른 무언가인지 확인할 수 있습니다.
    이를 통해 다음 섹션에서 일반 문자열과 Path 객체의 차이를 설명할 것입니다.

워크플로우를 실행합니다:

```bash
nextflow run main.nf
```

??? success "명령 출력"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `main.nf` [romantic_chandrasekhar] DSL2 - revision: 5a4a89bc3a

    data/patientA_rep1_normal_R1_001.fastq.gz is of class java.lang.String
    ```

보시다시피, Nextflow는 우리가 작성한 그대로 문자열 경로를 출력했습니다.

이것은 단순한 텍스트 출력이며, Nextflow는 아직 이것으로 특별한 작업을 수행하지 않았습니다.
또한 Nextflow 입장에서 이것은 단순한 문자열(`java.lang.String` 클래스)임을 확인했습니다.
아직 Nextflow에게 이것이 파일에 해당한다고 알려주지 않았으므로 당연한 결과입니다.

### 1.2. file()로 Path 객체 생성

[Path 객체](https://www.nextflow.io/docs/latest/reference/stdlib-types.html#path)를 경로 문자열에서 생성하여 Nextflow가 파일을 처리하는 방법을 알려줄 수 있습니다.

워크플로우에서 `file()` 메서드를 사용하여 문자열 경로 `data/patientA_rep1_normal_R1_001.fastq.gz`를 Path 객체로 변환할 수 있습니다. 이를 통해 파일 속성과 작업에 접근할 수 있습니다.

다음과 같이 `main.nf`를 편집하여 문자열을 `file()`로 감쌉니다:

=== "후"

    ```groovy title="main.nf" linenums="5" hl_lines="2"
        // 문자열 경로에서 Path 객체 생성
        myFile = file('data/patientA_rep1_normal_R1_001.fastq.gz')

        println "${myFile} is of class ${myFile.class}"
    ```

=== "전"

    ```groovy title="main.nf" linenums="5" hl_lines="2"
        // 문자열 경로에서 Path 객체 생성
        myFile = 'data/patientA_rep1_normal_R1_001.fastq.gz'

        println "${myFile} is of class ${myFile.class}"
    ```

워크플로우를 다시 실행합니다:

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

Nextflow가 문자열을 Path 객체로 변환하고 시스템의 실제 파일 위치로 확인했습니다.
파일 경로는 이제 `/workspaces/training/side-quests/working_with_files/data/patientA_rep1_normal_R1_001.fastq.gz`와 같이 절대 경로가 됩니다.

Path 객체의 클래스가 `sun.nio.fs.UnixPath`임을 확인할 수 있습니다. 이것은 Nextflow가 로컬 파일을 표현하는 방식입니다.
나중에 살펴보겠지만, 원격 파일은 다른 클래스 이름을 가집니다(예: HTTP 파일의 경우 `nextflow.file.http.XPath`). 하지만 모두 동일한 방식으로 작동하며 워크플로우에서 동일하게 사용할 수 있습니다.

!!! tip "팁"

    **핵심 차이점:**

    - **Path 문자열**: Nextflow가 문자로 처리하는 단순한 텍스트
    - **Path 객체**: Nextflow가 작업할 수 있는 스마트한 파일 참조

    이렇게 생각해 보세요: 경로 문자열은 종이에 주소를 적는 것과 같고, Path 객체는 GPS 장치에 주소를 입력하여 경로를 안내받고 여정에 대한 세부 정보를 알 수 있는 것과 같습니다.

### 1.3. 파일 속성 접근

이것이 왜 유용할까요? Nextflow가 `myFile`이 단순한 문자열이 아닌 Path 객체임을 이해하면, Path 객체의 다양한 속성에 접근할 수 있습니다.

내장 파일 속성을 출력하도록 워크플로우를 업데이트합니다:

=== "후"

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

=== "전"

    ```groovy title="main.nf" linenums="5" hl_lines="4"
        // 문자열 경로에서 Path 객체 생성
        myFile = file('data/patientA_rep1_normal_R1_001.fastq.gz')

        println "${myFile} is of class ${myFile.class}"
    ```

워크플로우를 실행합니다:

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

위에서 다양한 파일 속성이 콘솔에 출력된 것을 확인할 수 있습니다.

### 1.4. 프로세스에 파일 전달

문자열과 Path 객체의 차이는 실제 워크플로우를 프로세스와 함께 구축할 때 중요해집니다.
지금까지 Nextflow가 입력 파일을 파일로 처리하고 있음을 확인했지만, 실제로 프로세스에서 해당 파일에 대해 무언가를 실행할 수 있는지 확인해 봅니다.

#### 1.4.1. 프로세스 가져오기 및 코드 검토

파일 입력을 받아 파일의 줄 수를 세는 `COUNT_LINES`라는 미리 작성된 프로세스 모듈을 제공합니다.

워크플로우에서 이 프로세스를 사용하려면 workflow 블록 앞에 include 문을 추가하면 됩니다:

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

보시다시피, 파일의 압축을 풀고 줄 수를 세는 간단한 스크립트입니다.

??? info "`debug true`는 무엇을 하나요?"

    프로세스 정의의 `debug true` 지시문은 Nextflow가 스크립트의 출력(예: 줄 수 "40")을 실행 로그에 직접 출력하도록 합니다.
    이 지시문이 없으면 프로세스 실행 상태만 볼 수 있고 스크립트의 실제 출력은 볼 수 없습니다.

    Nextflow 프로세스 디버깅에 대한 자세한 내용은 [Nextflow 워크플로우 디버깅](debugging.md) 사이드 퀘스트를 참조하세요.

#### 1.4.2. `COUNT_LINES` 실행 추가

프로세스를 워크플로우에서 사용할 수 있게 되었으므로, 입력 파일에 대해 `COUNT_LINES` 프로세스를 실행하는 호출을 추가할 수 있습니다.

워크플로우를 다음과 같이 편집합니다:

=== "후"

    ```groovy title="main.nf" linenums="7" hl_lines="11-12"
        // 문자열 경로에서 Path 객체 생성
        myFile = file('data/patientA_rep1_normal_R1_001.fastq.gz')

        // 파일 속성 출력
        println "File object class: ${myFile.class}"
        println "File name: ${myFile.name}"
        println "Simple name: ${myFile.simpleName}"
        println "Extension: ${myFile.extension}"
        println "Parent directory: ${myFile.parent}"

        // 파일의 줄 수 세기
        COUNT_LINES(myFile)
    ```

=== "전"

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

워크플로우를 실행합니다:

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

프로세스 내에서 파일을 적절히 처리할 수 있음을 보여줍니다.

구체적으로, Nextflow는 다음 작업을 성공적으로 수행했습니다:

- 파일을 work directory에 스테이징
- .gz 파일 압축 해제
- 줄 수 계산(이 경우 40줄)
- 오류 없이 완료

이 원활한 작업의 핵심은 입력이 파일임을 Nextflow에 명시적으로 알려주는 것입니다.

### 1.5. 기본 파일 입력 오류 해결

이것은 Nextflow 초보자들이 자주 실수하는 부분이므로, 잘못 처리했을 때 어떤 일이 발생하는지 살펴봅니다.

파일 처리를 잘못할 수 있는 두 가지 주요 위치가 있습니다: 워크플로우 수준과 프로세스 수준입니다.

#### 1.5.1. 워크플로우 수준 오류

워크플로우 블록에서 입력을 지정할 때 파일을 문자열로 처리하도록 되돌리면 어떻게 되는지 확인합니다.

경로별 출력 문을 주석 처리하면서 워크플로우를 다음과 같이 편집합니다:

=== "후"

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

        // 파일의 줄 수 세기
        COUNT_LINES(myFile)
    ```

=== "전"

    ```groovy title="main.nf" linenums="7" hl_lines="4-9"
        // 문자열 경로에서 Path 객체 생성
        myFile = file('data/patientA_rep1_normal_R1_001.fastq.gz')

        // 파일 속성 출력
        println "File object class: ${myFile.class}"
        println "File name: ${myFile.name}"
        println "Simple name: ${myFile.simpleName}"
        println "Extension: ${myFile.extension}"
        println "Parent directory: ${myFile.parent}"

        // 파일의 줄 수 세기
        COUNT_LINES(myFile)
    ```

워크플로우를 실행합니다:

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

중요한 부분은 다음입니다:

```console
Not a valid path value: 'data/patientA_rep1_normal_R1_001.fastq.gz'
```

`path` 입력을 지정하면 Nextflow는 실제 파일 참조를 전달하고 있는지 검증합니다. 단순한 문자열은 허용되지 않습니다.
이 오류는 `'data/patientA_rep1_normal_R1_001.fastq.gz'`가 Path 객체가 아닌 문자열이기 때문에 유효한 경로 값이 아님을 알려줍니다.

Nextflow는 즉시 문제를 감지하고 프로세스를 시작하기 전에 중단했습니다.

#### 1.5.2. 프로세스 수준 오류

입력을 파일로 처리하도록 지정하는 것을 잊을 수 있는 또 다른 위치는 프로세스 정의입니다.

!!! warning "1.5.1의 워크플로우 오류 유지"

    이 테스트가 올바르게 작동하려면 워크플로우를 오류 상태(`file()` 대신 일반 문자열 사용)로 유지하세요.
    프로세스에서 `val`과 결합하면 아래에 표시된 오류가 발생합니다.

모듈을 다음과 같이 편집합니다:

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

워크플로우를 다시 실행합니다:

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

위에서 언급한 것처럼 프로세스가 디버깅 정보를 출력하도록 설정되어 있어 오류에 대한 많은 세부 정보가 표시됩니다.

가장 관련성 높은 섹션은 다음과 같습니다:

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

시스템이 파일을 찾을 수 없다고 하지만, 해당 경로를 확인하면 그 위치에 해당 이름의 파일이 있습니다.

이 경우 Nextflow는 문자열 값을 스크립트에 전달했지만, 실제 파일을 work directory에 *스테이징*하지 않았습니다.
따라서 프로세스는 상대 문자열 `data/patientA_rep1_normal_R1_001.fastq.gz`를 사용하려 했지만, 해당 파일이 프로세스 work directory 내에 존재하지 않습니다.

이 두 가지 예시를 통해 입력을 파일로 처리해야 한다고 Nextflow에 알려주는 것이 얼마나 중요한지 알 수 있습니다.

!!! note "참고"

    다음 섹션으로 진행하기 전에 두 가지 의도적인 오류를 모두 수정하세요.

### 핵심 정리

- Path 문자열 vs Path 객체: 문자열은 단순한 텍스트이고, Path 객체는 스마트한 파일 참조입니다
- `file()` 메서드는 문자열 경로를 Nextflow가 작업할 수 있는 Path 객체로 변환합니다
- [파일 속성 사용](https://www.nextflow.io/docs/latest/working-with-files.html#getting-file-attributes)을 통해 `name`, `simpleName`, `extension`, `parent`와 같은 파일 속성에 접근할 수 있습니다
- 문자열 대신 Path 객체를 사용하면 Nextflow가 워크플로우에서 파일을 올바르게 관리할 수 있습니다
- 프로세스 입력 결과: 올바른 파일 처리를 위해서는 문자열이 아닌 Path 객체가 필요하며, 이를 통해 파일이 올바르게 스테이징되고 프로세스에서 접근 가능합니다.

---

## 2. 원격 파일 사용

Nextflow의 주요 기능 중 하나는 로컬 파일(동일한 머신)과 인터넷을 통해 접근 가능한 원격 파일 간에 원활하게 전환할 수 있다는 것입니다.

올바르게 사용한다면, 다른 위치에서 오는 파일을 처리하기 위해 워크플로우 로직을 변경할 필요가 없습니다.
원격 파일을 사용하려면 워크플로우에 파일 경로를 제공할 때 적절한 접두사를 지정하기만 하면 됩니다.

예를 들어, `/path/to/data`는 접두사가 없어 '일반' 로컬 파일 경로임을 나타내고, `s3://path/to/data`는 `s3://` 접두사를 포함하여 Amazon의 S3 객체 스토리지에 위치함을 나타냅니다.

다양한 프로토콜이 지원됩니다:

- HTTP(S)/FTP (http://, https://, ftp://)
- Amazon S3 (s3://)
- Azure Blob Storage (az://)
- Google Cloud Storage (gs://)

이 중 하나를 사용하려면 문자열에 관련 접두사를 지정하면 됩니다. 이것은 기술적으로 파일 경로 대신 URI(Uniform Resource Identifier)라고 합니다.
Nextflow는 인증을 처리하고 파일을 올바른 위치에 스테이징하며, 다운로드, 업로드 및 기타 모든 파일 작업을 처리합니다.

이 시스템의 핵심 강점은 파이프라인 로직을 변경하지 않고 환경 간에 전환할 수 있다는 것입니다.
예를 들어, URI만 변경하여 소규모 로컬 테스트 세트로 개발한 후 원격 스토리지에 있는 전체 규모 테스트 세트로 전환할 수 있습니다.

### 2.1. 인터넷의 파일 사용

로컬 경로를 Github에 저장된 동일한 데이터의 복사본을 가리키는 HTTPS 경로로 전환하여 테스트합니다.

!!! warning "경고"

    인터넷 연결이 활성화된 경우에만 작동합니다.

`main.nf`를 다시 열고 입력 경로를 다음과 같이 변경합니다:

=== "후"

    ```groovy title="main.nf" linenums="2" hl_lines="2"
        // 인터넷의 원격 파일 사용
        myFile = file('https://raw.github.com/nextflow-io/training/master/side-quests/working_with_files/data/patientA_rep1_normal_R1_001.fastq.gz')

        // 파일 속성 출력
        println "File object class: ${myFile.class}"
        println "File name: ${myFile.name}"
        println "Simple name: ${myFile.simpleName}"
        println "Extension: ${myFile.extension}"
        println "Parent directory: ${myFile.parent}"
    ```

=== "전"

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

워크플로우를 실행합니다:

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

작동합니다! 거의 변경된 것이 없음을 확인할 수 있습니다.

콘솔 출력에서 한 가지 차이점은 경로 객체 클래스가 이제 `nextflow.file.http.XPath`인 반면, 로컬 경로의 경우 클래스는 `sun.nio.fs.UnixPath`였습니다.
이 클래스들을 기억할 필요는 없습니다. Nextflow가 다른 위치를 적절히 식별하고 처리한다는 것을 보여주기 위해 언급한 것입니다.

내부적으로 Nextflow는 work directory 내의 스테이징 디렉토리에 파일을 다운로드했습니다.
그런 다음 스테이징된 파일은 로컬 파일로 처리되어 관련 프로세스 디렉토리에 심볼릭 링크로 연결됩니다.

프로세스의 해시 값에 위치한 work directory의 내용을 확인하여 이를 검증할 수 있습니다.

??? abstract "작업 디렉토리 내용"

    프로세스 해시가 `8a/2ab7ca`인 경우 work directory를 탐색할 수 있습니다:

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

대용량 파일의 경우 로컬 파일 실행에 비해 다운로드 단계에 추가 시간이 소요됩니다.
그러나 Nextflow는 불필요한 다운로드를 방지하기 위해 이미 스테이징된 복사본이 있는지 확인합니다.
따라서 동일한 파일로 다시 실행하고 스테이징된 파일을 삭제하지 않았다면, Nextflow는 스테이징된 복사본을 사용합니다.

이것은 Nextflow의 핵심 기능인 로컬 데이터와 원격 데이터 간에 얼마나 쉽게 전환할 수 있는지를 보여줍니다.

!!! note "참고"

    이 원칙의 한 가지 중요한 예외는 HTTPS는 여러 파일을 나열할 수 없기 때문에 glob 패턴이나 디렉토리 경로를 사용할 수 없다는 것입니다. 정확한 파일 URL을 지정해야 합니다.
    그러나 blob 스토리지(`s3://`, `az://`, `gs://`)와 같은 다른 스토리지 프로토콜은 glob과 디렉토리 경로를 모두 사용할 수 있습니다.

    클라우드 스토리지에서 glob 패턴을 사용하는 방법은 다음과 같습니다:

    ```groovy title="Cloud storage examples (not runnable in this environment)"
    // glob 패턴을 사용한 S3 - 여러 파일과 매칭됩니다
    ch_s3_files = channel.fromPath('s3://my-bucket/data/*.fastq.gz')

    // glob 패턴을 사용한 Azure Blob Storage
    ch_azure_files = channel.fromPath('az://container/data/patient*_R{1,2}.fastq.gz')

    // glob 패턴을 사용한 Google Cloud Storage
    ch_gcs_files = channel.fromPath('gs://bucket/data/sample_*.fastq.gz')
    ```

    다음 섹션에서 실제로 glob을 사용하는 방법을 보여드리겠습니다.

### 2.2. 로컬 파일로 다시 전환

이 사이드 퀘스트의 나머지 부분에서는 로컬 예제 파일을 계속 사용할 것이므로, 워크플로우 입력을 원래 파일로 다시 전환합니다:

=== "후"

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

=== "전"

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

### 핵심 정리

- 원격 데이터는 URI(HTTP, FTP, S3, Azure, Google Cloud)를 사용하여 접근합니다
- Nextflow는 이러한 경로가 프로세스에 전달되는 한 데이터를 자동으로 다운로드하고 올바른 위치에 스테이징합니다
- 원격 파일을 다운로드하거나 업로드하는 로직을 직접 작성하지 마세요!
- 로컬 파일과 원격 파일은 다른 객체 유형을 생성하지만 동일하게 작동합니다
- **중요**: HTTP/HTTPS는 단일 파일에서만 작동합니다(glob 패턴 불가)
- 클라우드 스토리지(S3, Azure, GCS)는 단일 파일과 glob 패턴을 모두 지원합니다
- 코드 로직을 변경하지 않고 로컬 데이터 소스와 원격 데이터 소스 간에 원활하게 전환할 수 있습니다(프로토콜이 필요한 작업을 지원하는 한)

---

## 3. `fromPath()` 채널 팩토리 사용

지금까지는 한 번에 하나의 파일로 작업했지만, Nextflow에서는 일반적으로 처리할 여러 입력 파일이 있는 입력 채널을 생성하려고 합니다.

단순한 방법으로는 `file()` 메서드와 [`channel.of()`](https://www.nextflow.io/docs/latest/reference/channel.html#of)를 다음과 같이 결합할 수 있습니다:

```groovy title="Syntax example"
ch_files = channel.of([file('data/patientA_rep1_normal_R1_001.fastq.gz')],
                      [file('data/patientA_rep1_normal_R1_001.fastq.gz')])
```

작동하기는 하지만 번거롭습니다.

!!! tip "`file()` vs `channel.fromPath()` 사용 시기"

    - 단일 Path 객체가 필요한 경우(파일 존재 여부 확인, 속성 읽기, 단일 프로세스 실행에 전달)에는 `file()`을 사용하세요
    - 여러 파일을 담을 수 있는 채널이 필요한 경우, 특히 glob 패턴을 사용하거나 파일이 여러 프로세스를 통해 흐를 때는 `channel.fromPath()`를 사용하세요

여기서 [`channel.fromPath()`](https://www.nextflow.io/docs/latest/reference/channel.html#frompath)가 등장합니다. 하나 이상의 정적 파일 문자열과 glob 패턴에서 채널을 생성하는 데 필요한 모든 기능을 번들로 제공하는 편리한 채널 팩토리입니다.

### 3.1. 채널 팩토리 추가

`channel.fromPath`를 사용하도록 워크플로우를 업데이트합니다.

=== "후"

    ```groovy title="main.nf" linenums="7" hl_lines="1-3"
        // channel.fromPath로 파일 로드
        ch_files = channel.fromPath('data/patientA_rep1_normal_R1_001.fastq.gz')
        ch_files.view { myFile -> "Found file: $myFile" }

        // 파일 속성 출력
        /* 지금은 주석 처리합니다. 나중에 다시 살펴볼 것입니다!
        println "File object class: ${myFile.class}"
        println "File name: ${myFile.name}"
        println "Simple name: ${myFile.simpleName}"
        println "Extension: ${myFile.extension}"
        println "Parent directory: ${myFile.parent}"
        */

        // 파일의 줄 수 세기
        // COUNT_LINES(myFile)
    ```

=== "전"

    ```groovy title="main.nf" linenums="7" hl_lines="1-2"
        // 문자열 경로에서 Path 객체 생성
        myFile = file('data/patientA_rep1_normal_R1_001.fastq.gz')

        // 파일 속성 출력
        println "File object class: ${myFile.class}"
        println "File name: ${myFile.name}"
        println "Simple name: ${myFile.simpleName}"
        println "Extension: ${myFile.extension}"
        println "Parent directory: ${myFile.parent}"

        // 파일의 줄 수 세기
        COUNT_LINES(myFile)
    ```

속성을 출력하는 코드는 지금은 주석 처리하고, 대신 파일 이름만 출력하는 `.view` 문을 추가했습니다.

워크플로우를 실행합니다:

```bash
nextflow run main.nf
```

??? success "명령 출력"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `main.nf` [grave_meucci] DSL2 - revision: b09964a583

    Found file: /workspaces/training/side-quests/working_with_files/data/patientA_rep1_normal_R1_001.fastq.gz
    ```

보시다시피, 파일 경로가 채널에서 `Path` 유형 객체로 로드되고 있습니다.
이것은 `file()`이 했을 것과 유사하지만, 이제 원하는 경우 더 많은 파일을 로드할 수 있는 채널이 생겼습니다.

`channel.fromPath()`는 파일 목록으로 채워진 새 채널을 생성하는 편리한 방법입니다.

### 3.2. 채널의 파일 속성 보기

채널 팩토리를 처음 사용할 때는 코드를 단순화하여 파일 이름만 출력했습니다.

전체 파일 속성을 다시 출력하도록 돌아갑니다:

=== "후"

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

        // 파일의 줄 수 세기
        COUNT_LINES(ch_files)
    ```

=== "전"

    ```groovy title="main.nf" linenums="7" hl_lines="3"
        // channel.fromPath로 파일 로드
        ch_files = channel.fromPath('data/patientA_rep1_normal_R1_001.fastq.gz')
        ch_files.view { myFile -> "Found file: $myFile" }

        // 파일의 줄 수 세기
        // COUNT_LINES(ch_files)
    ```

채널 기반 접근 방식에서도 파일 처리가 올바르게 작동하는지 확인하기 위해 `COUNT_LINES` 프로세스 실행도 다시 활성화했습니다.

워크플로우를 실행합니다:

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

이전과 동일한 결과이지만, 이제 파일이 채널에 있으므로 더 많은 파일을 추가할 수 있습니다.

### 3.3. glob을 사용하여 여러 파일 매칭

채널에 더 많은 파일을 로드하는 방법은 여러 가지가 있습니다.
여기서는 와일드카드 문자를 기반으로 파일 및 디렉토리 이름을 매칭하고 검색하는 편리한 방법인 glob 패턴을 사용하는 방법을 보여드립니다.
이러한 패턴을 매칭하는 과정을 "globbing" 또는 "파일 이름 확장"이라고 합니다.

!!! note "참고"

    앞서 언급한 것처럼, Nextflow는 HTTPS가 여러 파일을 나열할 수 없기 때문에 HTTPS 파일 경로를 제외한 대부분의 경우에 입력 및 출력 파일을 관리하기 위한 globbing을 지원합니다.

특정 환자 `patientA`와 관련된 파일 쌍의 두 파일을 모두 검색하려고 한다고 가정합니다:

```console
patientA_rep1_normal_R1_001.fastq.gz
patientA_rep1_normal_R2_001.fastq.gz
```

파일 이름의 유일한 차이점은 리플리케이트 번호, 즉 `R` 뒤의 숫자이므로, 와일드카드 문자 `*`를 사용하여 다음과 같이 숫자를 대체할 수 있습니다:

```console
patientA_rep1_normal_R*_001.fastq.gz
```

이것이 필요한 glob 패턴입니다.

이제 채널 팩토리의 파일 경로를 해당 glob 패턴을 사용하도록 업데이트하면 됩니다:

=== "후"

    ```groovy title="main.nf" linenums="7"
      // channel.fromPath로 파일 로드
        ch_files = channel.fromPath('data/patientA_rep1_normal_R*_001.fastq.gz')
    ```

=== "전"

    ```groovy title="main.nf" linenums="7"
      // channel.fromPath로 파일 로드
        ch_files = channel.fromPath('data/patientA_rep1_normal_R1_001.fastq.gz')
    ```

Nextflow는 이것이 glob 패턴임을 자동으로 인식하고 적절히 처리합니다.

워크플로우를 실행하여 테스트합니다:

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

보시다시피, 이제 채널에 두 개의 Path 객체가 있습니다. Nextflow가 파일 이름 확장을 올바르게 수행하고 두 파일을 모두 로드하고 처리했음을 보여줍니다.

이 방법을 사용하면 glob 패턴을 변경하는 것만으로 원하는 만큼 많거나 적은 파일을 검색할 수 있습니다. 패턴을 더 넓게 만들면, 예를 들어 파일 이름의 모든 가변 부분을 `*`로 대체하면(_예:_ `data/patient*_rep*_*_R*_001.fastq.gz`) `data` 디렉토리의 모든 예제 파일을 가져올 수 있습니다.

### 핵심 정리

- `channel.fromPath()`는 패턴과 일치하는 파일로 채널을 생성합니다
- 각 파일은 채널의 별도 요소로 내보내집니다
- glob 패턴을 사용하여 여러 파일을 매칭할 수 있습니다
- 파일은 자동으로 전체 속성을 가진 Path 객체로 변환됩니다
- `.view()` 메서드를 사용하여 채널 내용을 검사할 수 있습니다

---

## 4. 파일 이름에서 기본 메타데이터 추출

대부분의 과학 분야에서 데이터를 포함하는 파일의 이름에 메타데이터가 인코딩되는 것은 매우 일반적입니다.
예를 들어, 생물정보학에서 시퀀싱 데이터를 포함하는 파일은 종종 샘플, 조건, 리플리케이트 및 리드 번호에 대한 정보를 인코딩하는 방식으로 명명됩니다.

파일 이름이 일관된 규칙에 따라 구성된 경우, 표준화된 방식으로 해당 메타데이터를 추출하여 분석 과정에서 사용할 수 있습니다.
물론 이것은 큰 '만약'이며, 파일 이름 구조에 의존할 때는 매우 주의해야 합니다. 하지만 이 접근 방식이 매우 널리 사용되는 것이 현실이므로, Nextflow에서 어떻게 수행되는지 살펴봅니다.

예제 데이터의 경우, 파일 이름에 일관되게 구조화된 메타데이터가 포함되어 있음을 알고 있습니다.
예를 들어, 파일 이름 `patientA_rep1_normal_R2_001`은 다음을 인코딩합니다:

- 환자 ID: `patientA`
- 리플리케이트 ID: `rep1`
- 샘플 유형: `normal` (`tumor`와 반대)
- 리드 세트: `R1` (`R2`와 반대)

세 단계로 이 정보를 검색하도록 워크플로우를 수정합니다:

1. 메타데이터를 포함하는 파일의 `simpleName` 검색
2. `tokenize()` 메서드를 사용하여 메타데이터 분리
3. map을 사용하여 메타데이터 구성

!!! warning "경고"

    환자 이름이나 기타 식별 특성과 같은 민감한 정보를 파일 이름에 인코딩해서는 안 됩니다. 이는 환자 개인 정보나 기타 관련 보안 제한을 침해할 수 있습니다.

### 4.1. `simpleName` 검색

`simpleName`은 경로와 확장자가 제거된 파일 이름에 해당하는 파일 속성입니다.

워크플로우를 다음과 같이 편집합니다:

=== "후"

    ```groovy title="main.nf" linenums="7" hl_lines="3-6"
        // channel.fromPath로 파일 로드
        ch_files = channel.fromPath('data/patientA_rep1_normal_R*_001.fastq.gz')
        ch_files.map { myFile ->
            [ myFile.simpleName, myFile ]
        }
        .view()
    ```

=== "전"

    ```groovy title="main.nf" linenums="7" hl_lines="3-9"
        // channel.fromPath로 파일 로드
        ch_files = channel.fromPath('data/patientA_rep1_normal_R*_001.fastq.gz')
        ch_files.view { myFile ->
            println "File object class: ${myFile.class}"
            println "File name: ${myFile.name}"
            println "Simple name: ${myFile.simpleName}"
            println "Extension: ${myFile.extension}"
            println "Parent directory: ${myFile.parent}"
        }
    ```

이것은 `map()` 작업을 사용하여 `simpleName`을 검색하고 전체 파일 객체와 연결합니다.

워크플로우를 실행하여 작동하는지 테스트합니다:

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

### 4.2. `simpleName`에서 메타데이터 추출

이 시점에서 원하는 메타데이터는 `simpleName`에 내포되어 있지만, 개별 항목에 직접 접근할 수 없습니다.
따라서 `simpleName`을 구성 요소로 분리해야 합니다.
다행히 원래 파일 이름에서 해당 구성 요소들은 단순히 밑줄로 구분되어 있으므로, 이 작업에 완벽한 `tokenize()`라는 일반적인 Nextflow 메서드를 적용할 수 있습니다.

워크플로우를 다음과 같이 편집합니다:

=== "후"

    ```groovy title="main.nf" linenums="7" hl_lines="4"
        // channel.fromPath로 파일 로드
        ch_files = channel.fromPath('data/patientA_rep1_normal_R*_001.fastq.gz')
        ch_files.map { myFile ->
            [ myFile.simpleName.tokenize('_'), myFile ]
        }
    ```

=== "전"

    ```groovy title="main.nf" linenums="7" hl_lines="4"
        // channel.fromPath로 파일 로드
        ch_files = channel.fromPath('data/patientA_rep1_normal_R*_001.fastq.gz')
        ch_files.map { myFile ->
            [ myFile.simpleName, myFile ]
        }
    ```

`tokenize()` 메서드는 밑줄을 찾을 때마다 `simpleName` 문자열을 분리하고, 부분 문자열을 포함하는 목록을 반환합니다.

워크플로우를 실행합니다:

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

이제 채널의 각 요소에 대한 튜플에는 메타데이터 목록(_예:_ `[patientA, rep1, normal, R1, 001]`)과 원래 파일 객체가 포함됩니다.

훌륭합니다!
환자 정보를 단일 문자열에서 문자열 목록으로 분해했습니다.
이제 환자 정보의 각 부분을 별도로 처리할 수 있습니다.

### 4.3. map을 사용하여 메타데이터 구성

현재 메타데이터는 단순한 평면 목록입니다.
사용하기는 쉽지만 읽기가 어렵습니다.

```console
[patientA, rep1, normal, R1, 001]
```

인덱스 3의 항목은 무엇인가요? 메타데이터 구조에 대한 원래 설명을 참조하지 않고 알 수 있나요?

이것은 키-값 저장소를 사용하기에 좋은 기회입니다. 모든 항목에 키 세트와 관련 값이 있어 각 키를 쉽게 참조하여 해당 값을 얻을 수 있습니다.

예시에서는 다음 구성에서:

```groovy
data = [patientA, 1, normal, R1]

println data[3]
```

다음으로 전환하는 것을 의미합니다:

```groovy
data = [id: patientA, replicate: 1, type: normal, readNum: 1]

println data.readNum
```

Nextflow에서 이것을 [map](https://nextflow.io/docs/latest/script.html#maps)이라고 합니다.

이제 평면 목록을 map으로 변환합니다.
워크플로우를 다음과 같이 편집합니다:

=== "후"

    ```groovy title="main.nf" linenums="7" hl_lines="4-13"
        // channel.fromPath로 파일 로드
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
        // channel.fromPath로 파일 로드
        ch_files = channel.fromPath('data/patientA_rep1_normal_R*_001.fastq.gz')
        ch_files.map { myFile ->
            [ myFile.simpleName.tokenize('_'), myFile ]
        }
    ```

여기서 주요 변경 사항은 다음과 같습니다:

- **구조 분해 할당**: `def (patient, replicate, type, readNum) = ...`는 토큰화된 값을 한 줄에서 명명된 변수로 추출합니다
- **map 리터럴 문법**: `[id: patient, replicate: ...]`는 각 키(예: `id`)가 값(예: `patient`)과 연결된 map을 생성합니다
- **내포된 구조**: 외부 목록 `[..., myFile]`은 메타데이터 map과 원래 파일 객체를 쌍으로 묶습니다

또한 불필요한 문자를 제거하기 위해 `replace()`라는 문자열 대체 메서드를 사용하여 일부 메타데이터 문자열을 단순화했습니다(_예:_ `replicate.replace('rep', '')`로 리플리케이트 ID에서 숫자만 유지).

워크플로우를 다시 실행합니다:

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

이제 메타데이터가 명확하게 레이블링되어(_예:_ `[id:patientA, replicate:1, type:normal, readNum:2]`) 무엇이 무엇인지 훨씬 쉽게 알 수 있습니다.

워크플로우에서 메타데이터 요소를 실제로 활용하기도 훨씬 쉬워지고, 코드를 더 읽기 쉽고 유지 관리하기 쉽게 만들 것입니다.

### 핵심 정리

- Nextflow에서 파일 이름을 완전한 프로그래밍 언어의 기능으로 처리할 수 있습니다
- 파일 이름을 문자열로 처리하여 관련 정보를 추출할 수 있습니다
- `tokenize()` 및 `replace()`와 같은 메서드를 사용하여 파일 이름의 문자열을 조작할 수 있습니다
- `.map()` 작업은 구조를 유지하면서 채널 요소를 변환합니다
- 구조화된 메타데이터(map)는 위치 목록보다 코드를 더 읽기 쉽고 유지 관리하기 쉽게 만듭니다

다음으로, 페어드 데이터 파일을 처리하는 방법을 살펴봅니다.

---

## 5. 페어드 데이터 파일 처리

많은 실험 설계에서 명시적으로 쌍을 이루어 처리하면 이점이 있는 페어드 데이터 파일이 생성됩니다.
예를 들어, 생물정보학에서 시퀀싱 데이터는 종종 페어드 리드 형태로 생성됩니다. 이는 동일한 DNA 단편에서 유래한 시퀀스 문자열을 의미합니다(반대쪽 끝에서 읽기 때문에 종종 'forward'와 'reverse'라고 불립니다).

이것이 R1과 R2가 두 세트의 리드를 나타내는 예제 데이터의 경우입니다.

```console
data/patientA_rep1_normal_R1_001.fastq.gz
data/patientA_rep1_normal_R2_001.fastq.gz
```

Nextflow는 이와 같은 페어드 파일을 처리하기 위한 특수 채널 팩토리인 `channel.fromFilePairs()`를 제공합니다. 이것은 공유 명명 패턴을 기반으로 파일을 자동으로 그룹화합니다. 이를 통해 더 적은 노력으로 페어드 파일을 더 긴밀하게 연결할 수 있습니다.

이를 활용하도록 워크플로우를 수정합니다.
두 단계가 필요합니다:

1. 채널 팩토리를 `channel.fromFilePairs()`로 전환
2. 메타데이터 추출 및 매핑

### 5.1. 채널 팩토리를 `channel.fromFilePairs()`로 전환

`channel.fromFilePairs`를 사용하려면 Nextflow가 쌍의 두 구성원을 식별하는 데 사용할 패턴을 지정해야 합니다.

예제 데이터로 돌아가서, 명명 패턴을 다음과 같이 공식화할 수 있습니다:

```console
data/patientA_rep1_normal_R{1,2}_001.fastq.gz
```

이것은 앞서 사용한 glob 패턴과 유사하지만, 쌍의 두 구성원을 식별하는 부분 문자열(R 바로 뒤에 오는 `1` 또는 `2`)을 구체적으로 열거합니다.

워크플로우 `main.nf`를 다음과 같이 업데이트합니다:

=== "후"

    ```groovy title="main.nf" linenums="7" hl_lines="1-2"
        // channel.fromFilePairs로 파일 로드
        ch_files = channel.fromFilePairs('data/patientA_rep1_normal_R{1,2}_001.fastq.gz')
        /* 지금은 매핑을 주석 처리합니다. 나중에 다시 살펴볼 것입니다!
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
        // channel.fromFilePairs로 파일 로드
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

채널 팩토리를 전환하고 파일 매칭 패턴을 조정했으며, map 작업을 주석 처리했습니다.
나중에 몇 가지 수정과 함께 다시 추가할 것입니다.

워크플로우를 실행하여 테스트합니다:

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

이런, 이번에는 실행이 실패했습니다!

오류 메시지의 관련 부분은 다음과 같습니다:

```console
Not a valid path value: 'patientA_rep1_normal_R'
```

채널 팩토리를 변경했기 때문입니다.
지금까지 원래 입력 채널에는 파일 경로만 포함되어 있었습니다.
수행한 모든 메타데이터 조작은 실제로 채널 내용에 영향을 미치지 않았습니다.

이제 `.fromFilePairs` 채널 팩토리를 사용하면 결과 채널의 내용이 다릅니다.
두 파일이 공유하는 `simpleName`의 일부(식별자 역할)와 두 파일 객체를 포함하는 튜플로 구성된 하나의 채널 요소만 볼 수 있습니다. 형식은 `id, [ file1, file2 ]`입니다.

Nextflow가 공유 접두사를 검사하여 환자 이름을 추출하고 환자 식별자로 사용하는 어려운 작업을 수행했기 때문에 훌륭합니다.

그러나 현재 워크플로우가 중단됩니다.
프로세스를 변경하지 않고 `COUNT_LINES`를 동일한 방식으로 계속 실행하려면 파일 경로를 추출하기 위한 매핑 작업을 적용해야 합니다.
하지만 우리의 궁극적인 목표는 파일 쌍을 적절히 처리하는 다른 프로세스인 `ANALYZE_READS`를 사용하는 것이므로 그렇게 하지 않을 것입니다.

따라서 `COUNT_LINES` 실행을 주석 처리(또는 삭제)하고 계속 진행합니다.

=== "후"

    ```groovy title="main.nf" linenums="26" hl_lines="2"
        // 파일의 줄 수 세기
        // COUNT_LINES(ch_files)
    ```

=== "전"

    ```groovy title="main.nf" linenums="26" hl_lines="2"
        // 파일의 줄 수 세기
        COUNT_LINES(ch_files)
    ```

`COUNT_LINES` include 문도 주석 처리하거나 삭제할 수 있지만, 기능적인 영향은 없습니다.

워크플로우를 다시 실행합니다:

```bash
nextflow run main.nf
```

??? success "명령 출력"

    ```console hl_lines="5"
     N E X T F L O W   ~  version 25.10.2

    Launching `main.nf` [fabulous_davinci] DSL2 - revision: 22b53268dc

    [patientA_rep1_normal_R, [/workspaces/training/side-quests/working_with_files/data/patientA_rep1_normal_R1_001.fastq.gz, /workspaces/training/side-quests/working_with_files/data/patientA_rep1_normal_R2_001.fastq.gz]]
    ```

이번에는 워크플로우가 성공했습니다!

그러나 `id` 필드에서 나머지 메타데이터를 추출해야 합니다.

### 5.2. 파일 쌍에서 메타데이터 추출 및 구성

이전의 `map` 작업은 데이터 구조와 일치하지 않아 작동하지 않지만, 수정하여 작동하게 할 수 있습니다.

`fromFilePairs()`가 식별자로 사용한 문자열에서 실제 환자 식별자에 이미 접근할 수 있으므로, 이전처럼 Path 객체에서 `simpleName`을 가져오지 않고 이를 사용하여 메타데이터를 추출할 수 있습니다.

워크플로우에서 map 작업의 주석을 해제하고 다음과 같이 편집합니다:

=== "후"

    ```groovy title="main.nf" linenums="7" hl_lines="3-4 9 11 13"
        // channel.fromFilePairs로 파일 로드
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
        // channel.fromFilePairs로 파일 로드
        ch_files = channel.fromFilePairs('data/patientA_rep1_normal_R{1,2}_001.fastq.gz')
        /* 지금은 매핑을 주석 처리합니다. 나중에 다시 살펴볼 것입니다!
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

이번에는 map이 `myFile` 대신 `id, files`에서 시작하고, `tokenize()`가 `myFile.simpleName` 대신 `id`에 적용됩니다.

또한 `tokenize()` 줄에서 `readNum`을 제거했습니다. 구체적으로 명명하지 않은 부분 문자열(왼쪽부터 시작)은 자동으로 삭제됩니다.
페어드 파일이 이제 긴밀하게 연결되어 있으므로 메타데이터 map에서 `readNum`이 더 이상 필요하지 않기 때문에 이렇게 할 수 있습니다.

워크플로우를 실행합니다:

```bash
nextflow run main.nf
```

??? success "명령 출력"

    ```console

    N E X T F L O W   ~  version 25.10.2

    Launching `main.nf` [prickly_stonebraker] DSL2 - revision: f62ab10a3f

    [[id:patientA, replicate:1, type:normal], [/workspaces/training/side-quests/working_with_files/data/patientA_rep1_normal_R1_001.fastq.gz, /workspaces/training/side-quests/working_with_files/data/patientA_rep1_normal_R2_001.fastq.gz]]
    ```

결과가 나왔습니다: 출력 튜플의 첫 번째 위치에 메타데이터 map(`[id:patientA, replicate:1, type:normal]`)이 있고, 그 뒤에 의도한 대로 페어드 파일의 튜플이 있습니다.

물론, 이것은 특정 파일 쌍만 선택하고 처리합니다.
여러 쌍을 처리하는 실험을 하려면 입력 패턴에 와일드카드를 추가하고 어떻게 되는지 확인해 보세요.
예를 들어, `data/patientA_rep1_*_R{1,2}_001.fastq.gz`를 사용해 보세요.

### 핵심 정리

- [`channel.fromFilePairs()`는 관련 파일을 자동으로 찾아 쌍으로 묶습니다](https://www.nextflow.io/docs/latest/reference/channel.html#fromfilepairs)
- 이를 통해 파이프라인에서 페어드 엔드 리드 처리가 단순화됩니다
- 페어드 파일은 `[id, [file1, file2]]` 튜플로 그룹화할 수 있습니다
- 메타데이터 추출은 개별 파일이 아닌 페어드 파일 ID에서 수행할 수 있습니다

---

## 6. 프로세스에서 파일 작업 사용

이제 이 모든 것을 간단한 프로세스에 통합하여 Nextflow 프로세스 내에서 파일 작업을 사용하는 방법을 강화합니다.

메타데이터 튜플과 한 쌍의 입력 파일을 받아 분석하는 `ANALYZE_READS`라는 미리 작성된 프로세스 모듈을 제공합니다.
이것이 시퀀스 정렬, 변이 호출 또는 이 데이터 유형에 적합한 다른 단계를 수행한다고 상상할 수 있습니다.

시작합니다.

### 6.1. 프로세스 가져오기 및 코드 검토

워크플로우에서 이 프로세스를 사용하려면 workflow 블록 앞에 모듈 include 문을 추가하면 됩니다.

워크플로우를 다음과 같이 편집합니다:

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

!!! note "참고"

    `tag`와 `publishDir` 지시문은 문자열 보간(`"${...}"`) 대신 closure 문법(`{ ... }`)을 사용합니다.
    이는 이러한 지시문이 런타임까지 사용할 수 없는 입력 변수(`meta`)를 참조하기 때문입니다.
    closure 문법은 프로세스가 실제로 실행될 때까지 평가를 지연시킵니다.

!!! note "참고"

    관례에 따라 메타데이터 map을 `meta`라고 부릅니다.
    meta map에 대한 자세한 내용은 [메타데이터와 meta map](../metadata/) 사이드 퀘스트를 참조하세요.

### 6.2. 워크플로우에서 프로세스 실행

프로세스를 워크플로우에서 사용할 수 있게 되었으므로, `ANALYZE_READS` 프로세스를 실행하는 호출을 추가할 수 있습니다.

예제 데이터에서 실행하려면 두 가지 작업이 필요합니다:

1. 재매핑된 채널에 이름 지정
2. 프로세스 실행 추가

#### 6.2.1. 재매핑된 입력 채널 이름 지정

이전에는 매핑 조작을 입력 채널에 직접 적용했습니다.
재매핑된 내용을 `ANALYZE_READS` 프로세스에 전달하기 위해(명확하고 읽기 쉬운 방식으로) `ch_samples`라는 새 채널을 생성합니다.

[`set`](https://www.nextflow.io/docs/latest/reference/operator.html#set) 연산자를 사용하여 이를 수행할 수 있습니다.

메인 워크플로우에서 `.view()` 연산자를 `.set { ch_samples }`로 교체하고, 이름으로 채널을 참조할 수 있는지 테스트하는 줄을 추가합니다.

=== "후"

    ```groovy title="main.nf" linenums="7" hl_lines="14 16-17"
        // channel.fromFilePairs로 파일 로드
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

        // 임시: ch_samples 내용 확인
        ch_samples.view()
    ```

=== "전"

    ```groovy title="main.nf" linenums="7" hl_lines="14"
        // channel.fromFilePairs로 파일 로드
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

실행합니다:

```bash
nextflow run main.nf
```

??? success "명령 출력"

    ```console
     N E X T F L O W   ~  version 25.10.2

    Launching `main.nf` [goofy_kirch] DSL2 - revision: 3313283e42

    [[id:patientA, replicate:1, type:normal], [/workspaces/training/side-quests/working_with_files/data/patientA_rep1_normal_R1_001.fastq.gz, /workspaces/training/side-quests/working_with_files/data/patientA_rep1_normal_R2_001.fastq.gz]]
    ```

이제 이름으로 채널을 참조할 수 있음을 확인했습니다.

#### 6.2.2. 데이터에 프로세스 실행

이제 `ch_samples` 채널에 `ANALYZE_READS` 프로세스를 실제로 실행합니다.

메인 워크플로우에서 다음과 같이 코드를 변경합니다:

=== "후"

    ```groovy title="main.nf" linenums="23"
        // 분석 실행
        ANALYZE_READS(ch_samples)
    ```

=== "전"

    ```groovy title="main.nf" linenums="23"
        // 임시: ch_samples 내용 확인
        ch_samples.view()
    ```

실행합니다:

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

이 프로세스는 출력을 `results` 디렉토리에 게시하도록 설정되어 있으므로 확인해 봅니다.

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

프로세스가 입력을 받아 설계된 대로 환자 메타데이터를 포함하는 새 파일을 생성했습니다.
훌륭합니다!

### 6.3. 더 많은 환자 포함

물론, 이것은 단일 환자의 단일 파일 쌍만 처리하는 것으로, Nextflow로 기대하는 높은 처리량은 아닙니다.
한 번에 훨씬 더 많은 데이터를 처리하고 싶을 것입니다.

`channel.fromPath()`는 *glob*을 입력으로 허용한다는 것을 기억하세요. 즉, 패턴과 일치하는 임의의 수의 파일을 허용할 수 있습니다.
따라서 모든 환자를 포함하려면 앞서 언급한 것처럼 더 많은 환자를 포함하도록 입력 문자열을 수정하면 됩니다.

가능한 한 많은 파일을 포함하고 싶다고 가정합니다.
워크플로우를 다음과 같이 편집합니다:

=== "후"

    ```groovy title="main.nf" linenums="7" hl_lines="2"
        // channel.fromFilePairs로 파일 로드
        ch_files = channel.fromFilePairs('data/*_R{1,2}_001.fastq.gz')
    ```

=== "전"

    ```groovy title="main.nf" linenums="7" hl_lines="2"
        // channel.fromFilePairs로 파일 로드
        ch_files = channel.fromFilePairs('data/patientA_rep1_normal_R{1,2}_001.fastq.gz')
    ```

파이프라인을 다시 실행합니다:

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

이제 results 디렉토리에 사용 가능한 모든 데이터에 대한 결과가 포함되어야 합니다.

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

성공! 모든 환자를 한 번에 분석했습니다! 맞죠?

아마도 아닐 것입니다.
더 자세히 살펴보면 문제가 있습니다: patientA에 두 개의 리플리케이트가 있지만 출력 파일은 하나뿐입니다!
매번 출력 파일을 덮어쓰고 있습니다.

### 6.4. 게시된 파일을 고유하게 만들기

환자 메타데이터에 접근할 수 있으므로, 디렉토리 구조나 파일 이름 자체에 차별화 메타데이터를 포함하여 게시된 파일을 고유하게 만들 수 있습니다.

워크플로우를 다음과 같이 변경합니다:

=== "후"

    ```groovy title="modules/analyze_reads.nf" linenums="6"
        publishDir { "results/${meta.type}/${meta.id}/${meta.replicate}" }, mode: 'copy'
    ```

=== "전"

    ```groovy title="modules/analyze_reads.nf" linenums="6"
        publishDir { "results/${meta.id}" }, mode: 'copy'
    ```

여기서는 샘플 유형과 리플리케이트를 처리하기 위해 추가 디렉토리 수준을 사용하는 옵션을 보여주지만, 파일 이름 수준에서도 실험해 볼 수 있습니다.

이제 파이프라인을 한 번 더 실행하되, 깨끗한 작업 공간을 위해 먼저 results 디렉토리를 삭제하세요:

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

이제 results 디렉토리를 확인합니다:

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

모든 메타데이터가 깔끔하게 구성되었습니다. 성공입니다!

이와 같이 메타데이터를 map에 로드하면 훨씬 더 많은 작업을 수행할 수 있습니다:

1. 환자 속성을 기반으로 구성된 출력 디렉토리 생성
2. 환자 속성을 기반으로 프로세스에서 결정 내리기
3. 메타데이터 값을 기반으로 데이터 분리, 결합 및 재결합

메타데이터를 명시적으로 유지하고 데이터에 첨부하는 이 패턴(파일 이름에 인코딩하는 대신)은 견고하고 유지 관리 가능한 분석 워크플로우를 구축할 수 있게 하는 Nextflow의 핵심 모범 사례입니다.
이에 대한 자세한 내용은 [메타데이터와 meta map](../metadata/) 사이드 퀘스트에서 확인할 수 있습니다.

### 핵심 정리

- `publishDir` 지시문은 메타데이터 값을 기반으로 출력을 구성할 수 있습니다
- 튜플의 메타데이터는 결과의 구조화된 구성을 가능하게 합니다
- 이 접근 방식은 명확한 데이터 출처를 가진 유지 관리 가능한 워크플로우를 생성합니다
- 프로세스는 메타데이터와 파일의 튜플을 입력으로 받을 수 있습니다
- `tag` 지시문은 실행 로그에서 프로세스 식별을 제공합니다
- 워크플로우 구조는 채널 생성과 프로세스 실행을 분리합니다

---

## 요약

이 사이드 퀘스트에서는 기본 작업부터 파일 컬렉션을 처리하는 고급 기법까지 Nextflow에서 파일을 다루는 방법을 학습했습니다.

이러한 기법을 실제 작업에 적용하면 복잡한 명명 규칙을 가진 많은 수의 파일을 처리할 때 특히 더 효율적이고 유지 관리 가능한 워크플로우를 구축할 수 있습니다.

### 핵심 패턴

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

2.  **원격 파일 사용**: URI를 사용하여 로컬 파일과 원격 파일 간에 투명하게 전환하는 방법을 학습하여, 워크플로우 로직을 변경하지 않고 다양한 소스의 파일을 처리하는 Nextflow의 기능을 보여주었습니다.

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

3.  **`fromPath()` 채널 팩토리를 사용한 파일 로드:** `channel.fromPath()`로 파일 패턴에서 채널을 생성하고 객체 유형을 포함한 파일 속성을 확인했습니다.

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

4.  **파일 이름에서 환자 메타데이터 추출:** `tokenize()`와 `replace()`를 사용하여 파일 이름에서 메타데이터를 추출하고 구조화하여 구성된 map으로 변환했습니다.

    ```groovy
    def name = file.name.tokenize('_')
    def patientId = name[0]
    def replicate = name[1].replace('rep', '')
    def type = name[2]
    def readNum = name[3].replace('R', '')
    ```

5.  **channel.fromFilePairs로 단순화:** `channel.fromFilePairs()`를 사용하여 관련 파일을 자동으로 쌍으로 묶고 페어드 파일 ID에서 메타데이터를 추출했습니다.

    ```groovy
    ch_pairs = channel.fromFilePairs('data/*_R{1,2}_001.fastq.gz')
    ```

6.  **프로세스에서 파일 작업 사용:** 적절한 입력 처리를 통해 파일 작업을 Nextflow 프로세스에 통합하고, `publishDir`을 사용하여 메타데이터를 기반으로 출력을 구성했습니다.

    - 프로세스 입력에 meta map 연결

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

[사이드 퀘스트 메뉴](../)로 돌아가거나 페이지 오른쪽 하단의 버튼을 클릭하여 목록의 다음 주제로 이동하세요.
