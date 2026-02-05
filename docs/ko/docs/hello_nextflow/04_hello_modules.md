# 파트 4: Hello Modules

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } AI 지원 번역 - [자세히 알아보기 및 개선 제안](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

<!--
<div class="video-wrapper">
  <iframe width="560" height="315" src="https://www.youtube.com/embed/Xxp_menS0E8?si=0AWnXB7xqHAzJdJV&amp;list=PLPZ8WHdZGxmXiHf8B26oB_fTfoKQdhlik" title="YouTube video player" frameborder="0" allow="accelerometer; autoplay; clipboard-write; encrypted-media; gyroscope; picture-in-picture; web-share" referrerpolicy="strict-origin-when-cross-origin" allowfullscreen></iframe>
</div>

/// caption
:fontawesome-brands-youtube:{ .youtube } Nextflow YouTube 채널에서 [전체 재생목록](https://www.youtube.com/playlist?list=PLPZ8WHdZGxmXiHf8B26oB_fTfoKQdhlik)을 확인하십시오.

:green_book: 비디오 스크립트는 [여기](./transcripts/04_hello_modules.md)에서 확인하실 수 있습니다.
///
-->

이 섹션에서는 파이프라인의 개발 및 유지보수를 더 효율적이고 지속 가능하게 만들기 위해 workflow 코드를 구성하는 방법을 다룹니다.
구체적으로, [**모듈**](https://nextflow.io/docs/latest/module.html)을 사용하는 방법을 시연할 것입니다.

Nextflow에서 **모듈**은 독립적인 코드 파일이며, 종종 단일 process 정의를 캡슐화합니다.
Workflow에서 모듈을 사용하려면 workflow 코드 파일에 한 줄의 `include` 문만 추가하면 됩니다. 그런 다음 일반적으로 하는 것과 같은 방식으로 process를 workflow에 통합할 수 있습니다.
이렇게 하면 코드의 여러 복사본을 생성하지 않고도 여러 workflow에서 process 정의를 재사용할 수 있습니다.

Workflow 개발을 시작했을 때, 우리는 모든 것을 하나의 단일 코드 파일에 작성했습니다.
이제 process를 개별 모듈로 이동할 것입니다.

<figure class="excalidraw">
--8<-- "docs/en/docs/hello_nextflow/img/modules.svg"
</figure>

이렇게 하면 코드가 더 공유 가능하고, 유연하며, 유지보수하기 쉬워집니다.

??? info "이 섹션부터 시작하는 방법"

    이 섹션은 [Hello Nextflow](./index.md) 과정의 파트 1-3을 완료했다고 가정하지만, 해당 섹션에서 다룬 기본 사항에 익숙하다면 특별한 준비 없이 여기서 시작할 수 있습니다.

---

## 0. 준비 운동: `hello-modules.nf` 실행

시작점으로 workflow 스크립트 `hello-modules.nf`를 사용할 것입니다.
이 스크립트는 이 교육 과정의 파트 3을 완료하여 생성된 스크립트와 동일하지만, 출력 대상을 변경했습니다:

```groovy title="hello-modules.nf" linenums="37" hl_lines="3 7 11 15"
output {
    first_output {
        path 'hello_modules'
        mode 'copy'
    }
    uppercased {
        path 'hello_modules'
        mode 'copy'
    }
    collected {
        path 'hello_modules'
        mode 'copy'
    }
    batch_report {
        path 'hello_modules'
        mode 'copy'
    }
}
```

변경을 시작하기 전에 모든 것이 제대로 작동하는지 확인하기 위해 스크립트를 한 번 실행하십시오:

```bash
nextflow run hello-modules.nf
```

??? success "명령 출력"

    ```console
     N E X T F L O W   ~  version 25.10.2

    Launching `hello-modules.nf` [hopeful_avogadro] DSL2 - revision: b09af1237d

    executor >  local (7)
    [0f/8795c9] sayHello (3)       [100%] 3 of 3 ✔
    [6a/eb2510] convertToUpper (3) [100%] 3 of 3 ✔
    [af/479117] collectGreetings   [100%] 1 of 1 ✔
    ```

이전과 마찬가지로 `output` 블록에 지정된 디렉토리(여기서는 `results/hello_modules/`)에서 출력 파일을 찾을 수 있습니다.

??? abstract "디렉토리 내용"

    ```console
    results/hello_modules/
    ├── Bonjour-output.txt
    ├── COLLECTED-batch-output.txt
    ├── Hello-output.txt
    ├── Holà-output.txt
    ├── batch-report.txt
    ├── UPPER-Bonjour-output.txt
    ├── UPPER-Hello-output.txt
    └── UPPER-Holà-output.txt
    ```

이것이 정상적으로 작동했다면 workflow 코드를 모듈화하는 방법을 배울 준비가 되었습니다.

---

## 1. 모듈을 저장할 디렉토리 생성

모듈을 특정 디렉토리에 저장하는 것이 모범 사례입니다.
해당 디렉토리의 이름은 원하는 대로 지을 수 있지만, 관례상 `modules/`라고 부릅니다.

```bash
mkdir modules
```

---

## 2. `sayHello()`용 모듈 생성

가장 간단한 형태에서, 기존 process를 모듈로 전환하는 것은 복사-붙여넣기 작업에 불과합니다.
모듈에 대한 파일 스텁을 생성하고, 관련 코드를 복사한 다음 메인 workflow 파일에서 삭제할 것입니다.

그런 다음 `include` 문만 추가하면 Nextflow가 런타임에 관련 코드를 가져올 것입니다.

### 2.1. 새 모듈에 대한 파일 스텁 생성

`sayHello.nf`라는 모듈에 대한 빈 파일을 생성합시다.

```bash
touch modules/sayHello.nf
```

이것은 process 코드를 넣을 장소를 제공합니다.

### 2.2. `sayHello` process 코드를 모듈 파일로 이동

전체 process 정의를 workflow 파일에서 모듈 파일로 복사하고, `#!/usr/bin/env nextflow` shebang도 함께 복사해야 합니다.

```groovy title="modules/sayHello.nf" linenums="1"
#!/usr/bin/env nextflow

/*
 * echo를 사용하여 'Hello World!'를 파일에 출력
 */
process sayHello {

    input:
    val greeting

    output:
    path "${greeting}-output.txt"

    script:
    """
    echo '${greeting}' > '${greeting}-output.txt'
    """
}
```

완료되면 workflow 파일에서 process 정의를 삭제하되, shebang은 그대로 두십시오.

### 2.3. Workflow 블록 전에 include 선언 추가

모듈에서 process를 가져오는 구문은 매우 간단합니다:

```groovy title="구문: Include 선언"
include { <PROCESS_NAME> } from '<path_to_module>'
```

`params` 블록 위에 이것을 삽입하고 적절하게 채웁시다.

=== "수정 후"

    ```groovy title="hello-modules.nf" linenums="44" hl_lines="1 2"
    // 모듈 포함
    include { sayHello } from './modules/sayHello.nf'

    /*
    * 파이프라인 매개변수
    */
    params {
        greeting: Path = 'data/greetings.csv'
        batch: String = 'batch'
    }
    ```

=== "수정 전"

    ```groovy title="hello-modules.nf" linenums="44"
    /*
    * 파이프라인 매개변수
    */
    params {
        greeting: Path = 'data/greetings.csv'
        batch: String = 'batch'
    }
    ```

process 이름인 `sayHello`와 모듈 코드가 포함된 파일 경로인 `./modules/sayHello.nf`를 채웠습니다.

### 2.4. Workflow 실행

이전과 본질적으로 동일한 코드와 입력으로 workflow를 실행하고 있으므로, `-resume` 플래그로 실행하고 어떻게 되는지 살펴봅시다.

```bash
nextflow run hello-modules.nf -resume
```

??? success "명령 출력"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `hello-modules.nf` [romantic_poisson] DSL2 - revision: 96edfa9ad3

    [f6/cc0107] sayHello (1)       | 3 of 3, cached: 3 ✔
    [3c/4058ba] convertToUpper (2) | 3 of 3, cached: 3 ✔
    [1a/bc5901] collectGreetings   | 1 of 1, cached: 1 ✔
    ```

모든 것이 캐시되어 있으므로 매우 빠르게 실행되어야 합니다.
게시된 출력을 자유롭게 확인하십시오.

Nextflow는 코드가 여러 파일로 분할되어 있더라도 여전히 동일한 작업이 수행되어야 함을 인식했습니다.

### 핵심 정리

Process를 로컬 모듈로 추출하는 방법과 이렇게 해도 workflow의 재개 가능성이 깨지지 않는다는 것을 알게 되었습니다.

### 다음 단계

더 많은 모듈을 만드는 연습을 합니다.
하나를 만들면 백만 개 더 만들 수 있습니다...
하지만 지금은 두 개만 더 만들겠습니다.

---

## 3. `convertToUpper()` process 모듈화

### 3.1. 새 모듈에 대한 파일 스텁 생성

`convertToUpper.nf`라는 모듈에 대한 빈 파일을 생성합니다.

```bash
touch modules/convertToUpper.nf
```

### 3.2. `convertToUpper` process 코드를 모듈 파일로 이동

전체 process 정의를 workflow 파일에서 모듈 파일로 복사하고, `#!/usr/bin/env nextflow` shebang도 함께 복사해야 합니다.

```groovy title="modules/convertToUpper.nf" linenums="1"
#!/usr/bin/env nextflow

/*
 * 텍스트 대체 도구를 사용하여 인사말을 대문자로 변환
 */
process convertToUpper {

    input:
    path input_file

    output:
    path "UPPER-${input_file}"

    script:
    """
    cat '${input_file}' | tr '[a-z]' '[A-Z]' > 'UPPER-${input_file}'
    """
}
```

완료되면 workflow 파일에서 process 정의를 삭제하되, shebang은 그대로 두십시오.

### 3.3. `params` 블록 전에 include 선언 추가

`params` 블록 위에 include 선언을 삽입하고 적절하게 채웁니다.

=== "수정 후"

    ```groovy title="hello-modules.nf" linenums="23" hl_lines="3"
    // 모듈 포함
    include { sayHello } from './modules/sayHello.nf'
    include { convertToUpper } from './modules/convertToUpper.nf'

    /*
    * 파이프라인 매개변수
    */
    params {
        greeting: Path = 'data/greetings.csv'
        batch: String = 'batch'
    }
    ```

=== "수정 전"

    ```groovy title="hello-modules.nf" linenums="23"
    // 모듈 포함
    include { sayHello } from './modules/sayHello.nf'

    /*
    * 파이프라인 매개변수
    */
    params {
        greeting: Path = 'data/greetings.csv'
        batch: String = 'batch'
    }
    ```

이것이 매우 익숙해지기 시작해야 합니다.

### 3.4. Workflow 다시 실행

`-resume` 플래그로 실행합니다.

```bash
nextflow run hello-modules.nf -resume
```

??? success "명령 출력"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `hello-modules.nf` [nauseous_heisenberg] DSL2 - revision: a04a9f2da0

    [c9/763d42] sayHello (3)       | 3 of 3, cached: 3 ✔
    [60/bc6831] convertToUpper (3) | 3 of 3, cached: 3 ✔
    [1a/bc5901] collectGreetings   | 1 of 1, cached: 1 ✔
    ```

이것은 여전히 이전과 동일한 출력을 생성해야 합니다.

두 개 완료, 하나 더 남았습니다!

---

## 4. `collectGreetings()` process 모듈화

### 4.1. 새 모듈에 대한 파일 스텁 생성

`collectGreetings.nf`라는 모듈에 대한 빈 파일을 생성합니다.

```bash
touch modules/collectGreetings.nf
```

### 4.2. `collectGreetings` process 코드를 모듈 파일로 이동

전체 process 정의를 workflow 파일에서 모듈 파일로 복사하고, `#!/usr/bin/env nextflow` shebang도 함께 복사해야 합니다.

```groovy title="modules/collectGreetings.nf" linenums="1"
#!/usr/bin/env nextflow

/*
 * 대문자 인사말을 하나의 출력 파일에 수집
 */
process collectGreetings {

    input:
    path input_files
    val batch_name

    output:
    path "COLLECTED-${batch_name}-output.txt", emit: outfile
    path "${batch_name}-report.txt", emit: report

    script:
    count_greetings = input_files.size()
    """
    cat ${input_files} > 'COLLECTED-${batch_name}-output.txt'
    echo 'There were ${count_greetings} greetings in this batch.' > '${batch_name}-report.txt'
    """
}
```

완료되면 workflow 파일에서 process 정의를 삭제하되, shebang은 그대로 두십시오.

### 4.3. `params` 블록 전에 include 선언 추가

`params` 블록 위에 include 선언을 삽입하고 적절하게 채웁니다.

=== "수정 후"

    ```groovy title="hello-modules.nf" linenums="3" hl_lines="4"
    // 모듈 포함
    include { sayHello } from './modules/sayHello.nf'
    include { convertToUpper } from './modules/convertToUpper.nf'
    include { collectGreetings } from './modules/collectGreetings.nf'

    /*
    * 파이프라인 매개변수
    */
    params {
        greeting: Path = 'data/greetings.csv'
        batch: String = 'batch'
    }
    ```

=== "수정 전"

    ```groovy title="hello-modules.nf" linenums="3"
    // 모듈 포함
    include { sayHello } from './modules/sayHello.nf'
    include { convertToUpper } from './modules/convertToUpper.nf'

    /*
    * 파이프라인 매개변수
    */
    params {
        greeting: Path = 'data/greetings.csv'
        batch: String = 'batch'
    }
    ```

마지막입니다!

### 4.4. Workflow 실행

`-resume` 플래그로 실행합니다.

```bash
nextflow run hello-modules.nf -resume
```

??? success "명령 출력"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `hello-modules.nf` [friendly_coulomb] DSL2 - revision: 7aa2b9bc0f

    [f6/cc0107] sayHello (1)       | 3 of 3, cached: 3 ✔
    [3c/4058ba] convertToUpper (2) | 3 of 3, cached: 3 ✔
    [1a/bc5901] collectGreetings   | 1 of 1, cached: 1 ✔
    ```

이것은 여전히 이전과 동일한 출력을 생성해야 합니다.

### 핵심 정리

Workflow에서 여러 process를 모듈화하는 방법을 알게 되었습니다.

축하합니다, 이 모든 작업을 수행했지만 파이프라인 작동 방식에는 아무것도 변경되지 않았습니다!

농담은 제쳐두고, 이제 코드가 더 모듈화되었으며, 해당 process 중 하나를 호출하는 다른 파이프라인을 작성하기로 결정하면 관련 모듈을 사용하기 위해 짧은 `include` 문 하나만 입력하면 됩니다.
이것은 코드를 복사-붙여넣기하는 것보다 낫습니다. 나중에 모듈을 개선하기로 결정하면 모든 파이프라인이 개선 사항을 상속받기 때문입니다.

### 다음 단계

원한다면 잠시 휴식을 취하십시오.

준비가 되면 [**파트 5: Hello Containers**](./05_hello_containers.md)로 이동하여 소프트웨어 의존성을 더 편리하고 재현 가능하게 관리하기 위해 컨테이너를 사용하는 방법을 배우십시오.

---

## 퀴즈

<quiz>
Nextflow에서 모듈이란 무엇입니까?
- [ ] 구성 파일
- [x] Process 정의를 포함할 수 있는 독립 파일
- [ ] Workflow 정의
- [ ] Channel 연산자

자세히 알아보기: [2. `sayHello()`용 모듈 생성](#2-sayhello용-모듈-생성)
</quiz>

<quiz>
모듈 파일을 저장하는 데 일반적으로 사용되는 규칙은 무엇입니까?
- [ ] Workflow와 같은 디렉토리에
- [ ] `bin/` 디렉토리에
- [x] `modules/` 디렉토리에
- [ ] `lib/` 디렉토리에

자세히 알아보기: [1. 모듈을 저장할 디렉토리 생성](#1-모듈을-저장할-디렉토리-생성)
</quiz>

<quiz>
모듈을 사용하는 올바른 구문은 무엇입니까?

- [ ] `#!groovy import { SAYHELLO } from './modules/sayhello.nf'`
- [ ] `#!groovy require { SAYHELLO } from './modules/sayhello.nf'`
- [x] `#!groovy include { SAYHELLO } from './modules/sayhello.nf'`
- [ ] `#!groovy load { SAYHELLO } from './modules/sayhello.nf'`

자세히 알아보기: [2.3. include 선언 추가](#23-workflow-블록-전에-include-선언-추가)
</quiz>

<quiz>
모듈을 사용할 때 `-resume` 기능은 어떻게 됩니까?
- [ ] 더 이상 작동하지 않음
- [ ] 추가 구성이 필요함
- [x] 이전과 동일하게 작동함
- [ ] 로컬 모듈에서만 작동함
</quiz>

<quiz>
모듈 사용의 이점은 무엇입니까? (해당되는 것 모두 선택)
- [x] Workflow 간 코드 재사용 가능
- [x] 더 쉬운 유지보수
- [x] Workflow 코드의 더 나은 구성
- [ ] 더 빠른 실행 속도
</quiz>
