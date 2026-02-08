# 파트 5: Hello Containers

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } AI 지원 번역 - [자세히 알아보기 및 개선 제안](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

<div class="video-wrapper">
  <iframe width="560" height="315" src="https://www.youtube.com/embed/Xqr--bKEN9U?si=QinuAnFwFj-Z8CrO&amp;list=PLPZ8WHdZGxmWKozQuzr27jyMGqp9kElVK&amp;cc_load_policy=1&amp;cc_lang_pref=ko" title="YouTube video player" frameborder="0" allow="accelerometer; autoplay; clipboard-write; encrypted-media; gyroscope; picture-in-picture; web-share" referrerpolicy="strict-origin-when-cross-origin" allowfullscreen></iframe>
</div>

/// caption
:fontawesome-brands-youtube:{ .youtube } Nextflow YouTube 채널에서 [전체 재생목록](https://youtube.com/playlist?list=PLPZ8WHdZGxmWKozQuzr27jyMGqp9kElVK&si=eF7cLR62goy-lc6n)을 확인하십시오.

:green_book: 비디오 스크립트는 [여기](./transcripts/05_hello_containers.md)에서 확인하실 수 있습니다.
///

이 교육 과정의 파트 1-4에서는 Nextflow의 기본 구성 요소를 사용하여 텍스트를 처리하고, 여러 입력이 있는 경우 실행을 병렬화하며, 추가 처리를 위해 결과를 수집할 수 있는 간단한 워크플로우를 조립하는 방법을 배웠습니다.

그러나 환경에서 사용할 수 있는 기본 UNIX 도구로 제한되었습니다.
실제 작업에는 기본적으로 포함되지 않은 다양한 도구와 패키지가 필요한 경우가 많습니다.
일반적으로 이러한 도구를 설치하고, 의존성을 관리하며, 충돌을 해결해야 합니다.

이 모든 것이 매우 지루하고 귀찮으므로, 이 문제를 훨씬 더 편리하게 해결하기 위해 **컨테이너**를 사용하는 방법을 보여드리겠습니다.

**컨테이너**는 코드, 시스템 라이브러리 및 설정을 포함하여 애플리케이션을 실행하는 데 필요한 모든 것이 포함된 컨테이너 **이미지**에서 생성된 경량의 단독 실행형 소프트웨어 단위입니다.
예상하시다시피, 이것은 파이프라인을 더 재현 가능하게 만드는 데 매우 도움이 될 것입니다.

여기서는 [Docker](https://www.docker.com/get-started/)를 사용하여 이를 다루지만, Nextflow가 [여러 다른 컨테이너 기술](https://www.nextflow.io/docs/latest/container.html#)도 지원한다는 점을 유의하십시오.

??? info "이 섹션부터 시작하는 방법"

    이 섹션은 [Hello Nextflow](./index.md) 과정의 파트 1-4를 완료하고 완전히 작동하는 파이프라인이 있다고 가정합니다.

    이 시점부터 과정을 시작하는 경우 solutions에서 `modules` 디렉토리를 복사해야 합니다:

    ```bash
    cp -r solutions/4-hello-modules/modules .
    ```

---

## 0. 준비 운동: `hello-containers.nf` 실행

시작점으로 워크플로우 스크립트 `hello-containers.nf`를 사용할 것입니다.
이 스크립트는 이 교육 과정의 파트 4를 완료하여 생성된 스크립트와 동일하지만, 출력 대상을 변경했습니다:

```groovy title="hello-containers.nf" linenums="37" hl_lines="3 7 11 15"
output {
    first_output {
        path 'hello_containers'
        mode 'copy'
    }
    uppercased {
        path 'hello_containers'
        mode 'copy'
    }
    collected {
        path 'hello_containers'
        mode 'copy'
    }
    batch_report {
        path 'hello_containers'
        mode 'copy'
    }
}
```

변경을 시작하기 전에 모든 것이 제대로 작동하는지 확인하기 위해 스크립트를 한 번 실행하십시오:

```bash
nextflow run hello-containers.nf
```

??? success "명령 출력"

    ```console
     N E X T F L O W   ~  version 25.10.2

    Launching `hello-containers.nf` [nice_escher] DSL2 - revision: d5dfdc9872

    executor > local (7)
    [5a/ec1fa1] sayHello (2) [100%] 3 of 3 ✔
    [30/32b5b8] convertToUpper (3) [100%] 3 of 3 ✔
    [d3/be01bc] collectGreetings [100%] 1 of 1 ✔

    ```

이전과 마찬가지로 `output` 블록에 지정된 디렉토리(`results/hello_containers/`)에서 출력 파일을 찾을 수 있습니다.

??? abstract "디렉토리 내용"

    ```console
    results/hello_containers/
    ├── Bonjour-output.txt
    ├── COLLECTED-batch-output.txt
    ├── Hello-output.txt
    ├── Holà-output.txt
    ├── batch-report.txt
    ├── UPPER-Bonjour-output.txt
    ├── UPPER-Hello-output.txt
    └── UPPER-Holà-output.txt
    ```

이것이 정상적으로 작동했다면 컨테이너를 사용하는 방법을 배울 준비가 되었습니다.

---

## 1. 컨테이너를 '수동으로' 사용

우리가 하고자 하는 것은 실행에 컨테이너를 사용하는 단계를 워크플로우에 추가하는 것입니다.

그러나 먼저 Nextflow에서 컨테이너를 사용하기 전에 컨테이너가 무엇인지에 대한 이해를 확고히 하기 위해 몇 가지 기본 개념과 작업을 살펴볼 것입니다.

### 1.1. 컨테이너 이미지 풀(Pull)

컨테이너를 사용하려면 일반적으로 컨테이너 레지스트리에서 컨테이너 이미지를 다운로드하거나 *풀(pull)*한 다음, 컨테이너 이미지를 실행하여 컨테이너 인스턴스를 생성합니다.

일반적인 구문은 다음과 같습니다:

```bash title="구문"
docker pull '<container>'
```

`docker pull` 부분은 저장소에서 컨테이너 이미지를 풀(pull)하라는 컨테이너 시스템에 대한 명령입니다.

`'<container>'` 부분은 컨테이너 이미지의 URI 주소입니다.

예를 들어, 임의의 텍스트 입력을 재미있는 방식으로 표시하는 ASCII 아트를 생성하는 `cowsay`라는 도구의 Python 구현인 [cowpy](https://github.com/jeffbuttars/cowpy)가 포함된 컨테이너 이미지를 풀해 봅시다.

```txt title="예제"
 ________________________
< Are we having fun yet? >
 ------------------------
    \                                  ___-------___
     \                             _-~~             ~~-_
      \                         _-~                    /~-_
             /^\__/^\         /~  \                   /    \
           /|  O|| O|        /      \_______________/        \
          | |___||__|      /       /                \          \
          |          \    /      /                    \          \
          |   (_______) /______/                        \_________ \
          |         / /         \                      /            \
           \         \^\\         \                  /               \     /
             \         ||           \______________/      _-_       //\__//
               \       ||------_-~~-_ ------------- \ --/~   ~\    || __/
                 ~-----||====/~     |==================|       |/~~~~~
                  (_(__/  ./     /                    \_\      \.
                         (_(___/                         \_____)_)
```

게시된 컨테이너를 찾을 수 있는 다양한 저장소가 있습니다.
우리는 [Seqera Containers](https://seqera.io/containers/) 서비스를 사용하여 `cowpy` Conda 패키지에서 이 Docker 컨테이너 이미지를 생성했습니다: `'community.wave.seqera.io/library/cowpy:1.1.5--3db457ae1977a273'`.

전체 풀 명령을 실행하십시오:

```bash
docker pull 'community.wave.seqera.io/library/cowpy:1.1.5--3db457ae1977a273'
```

??? success "명령 출력"

    ```console
    1.1.5--3db457ae1977a273: Pulling from library/cowpy
    dafa2b0c44d2: Pull complete
    dec6b097362e: Pull complete
    f88da01cff0b: Pull complete
    4f4fb700ef54: Pull complete
    92dc97a3ef36: Pull complete
    403f74b0f85e: Pull complete
    10b8c00c10a5: Pull complete
    17dc7ea432cc: Pull complete
    bb36d6c3110d: Pull complete
    0ea1a16bbe82: Pull complete
    030a47592a0a: Pull complete
    c23bdb422167: Pull complete
    e1686ff32a11: Pull complete
    Digest: sha256:1ebc0043e8cafa61203bf42d29fd05bd14e7b4298e5e8cf986504c15f5aa4160
    Status: Downloaded newer image for community.wave.seqera.io/library/cowpy:1.1.5--3db457ae1977a273
    community.wave.seqera.io/library/cowpy:1.1.5--3db457ae1977a273
    ```

이전에 이미지를 다운로드한 적이 없다면 완료하는 데 1분 정도 걸릴 수 있습니다.
완료되면 컨테이너 이미지의 로컬 복사본이 생깁니다.

### 1.2. 컨테이너를 사용하여 `cowpy`를 일회성 명령으로 실행

사람들이 컨테이너를 사용하는 매우 일반적인 방법 중 하나는 직접 실행하는 것, 즉 비대화식으로 실행하는 것입니다.
이것은 일회성 명령을 실행하는 데 적합합니다.

일반적인 구문은 다음과 같습니다:

```bash title="구문"
docker run --rm '<container>' [tool command]
```

`docker run --rm '<container>'` 부분은 컨테이너 이미지에서 컨테이너 인스턴스를 생성하고 그 안에서 명령을 실행하라는 컨테이너 시스템에 대한 명령입니다.
`--rm` 플래그는 명령이 완료된 후 컨테이너 인스턴스를 종료하도록 시스템에 지시합니다.

`[tool command]` 구문은 사용 중인 도구와 컨테이너 설정 방식에 따라 다릅니다.
`cowpy`로 시작해 봅시다.

완전히 조립된 컨테이너 실행 명령은 다음과 같습니다. 계속 진행하여 실행하십시오.

```bash
docker run --rm 'community.wave.seqera.io/library/cowpy:1.1.5--3db457ae1977a273' cowpy
```

??? success "명령 출력"

    ```console
    ______________________________________________________
    < Cowacter, eyes:default, tongue:False, thoughts:False >
    ------------------------------------------------------
        \   ^__^
          \  (oo)\_______
            (__)\       )\/\
              ||----w |
              ||     ||
    ```

시스템이 컨테이너를 생성하고, 매개변수와 함께 `cowpy` 명령을 실행하고, 출력을 콘솔로 보내고, 마지막으로 컨테이너 인스턴스를 종료했습니다.

### 1.3. 컨테이너를 대화식으로 사용하여 `cowpy` 실행

컨테이너를 대화식으로 실행할 수도 있으며, 이렇게 하면 컨테이너 내부에서 셸 프롬프트를 받고 명령으로 작업할 수 있습니다.

#### 1.3.1. 컨테이너 생성

대화식으로 실행하려면 `docker run` 명령에 `-it`를 추가하기만 하면 됩니다.
선택적으로, 명령에 예를 들어 `/bin/bash`를 추가하여 컨테이너 내부에서 사용할 셸을 지정할 수 있습니다.

```bash
docker run --rm -it 'community.wave.seqera.io/library/cowpy:1.1.5--3db457ae1977a273' /bin/bash
```

프롬프트가 `(base) root@b645838b3314:/tmp#`과 같은 것으로 변경되는 것을 확인하십시오. 이는 이제 컨테이너 내부에 있음을 나타냅니다.

파일 시스템의 루트에서 디렉토리 내용을 나열하는 `ls /`를 실행하여 이를 확인할 수 있습니다:

```bash
ls /
```

??? abstract "명령 출력"

    ```console
    bin  boot  dev  etc  home  lib  lib64  media  mnt  opt  proc  root  run  sbin  srv  sys  tmp  usr  var
    ```

`tree` 유틸리티가 이 컨테이너에서 사용할 수 없기 때문에 여기서 `tree` 대신 `ls`를 사용합니다.
컨테이너 내부의 파일 시스템이 호스트 시스템의 파일 시스템과 다르다는 것을 알 수 있습니다.

방금 수행한 작업의 한 가지 제한 사항은 컨테이너가 기본적으로 호스트 시스템과 완전히 격리되어 있다는 것입니다.
이것은 명시적으로 허용하지 않는 한 컨테이너가 호스트 시스템의 파일에 액세스할 수 없음을 의미합니다.

잠시 후에 그 방법을 보여드리겠습니다.

#### 1.3.2. 원하는 도구 명령 실행

이제 컨테이너 내부에 있으므로 `cowpy` 명령을 직접 실행하고 몇 가지 매개변수를 제공할 수 있습니다.
예를 들어, 도구 문서에 따르면 `-c`로 캐릭터('cowacter')를 변경할 수 있습니다.

```bash
cowpy "Hello Containers" -c tux
```

??? success "명령 출력"

    ```console
    __________________
    < Hello Containers >
    ------------------
      \
        \
            .--.
          |o_o |
          |:_/ |
          //   \ \
        (|     | )
        /'\_   _/`\
        \___)=(___/
    ```

이제 출력에 기본 소 대신 리눅스 펭귄 Tux가 표시됩니다. `-c tux` 매개변수를 지정했기 때문입니다.

컨테이너 내부에 있으므로 Docker 명령을 다룰 필요 없이 입력 매개변수를 변경하면서 `cowpy` 명령을 원하는 만큼 실행할 수 있습니다.

!!! Tip "팁"

    '-c' 플래그를 사용하여 다른 캐릭터를 선택할 수 있습니다:
    `beavis`, `cheese`, `daemon`, `dragonandcow`, `ghostbusters`, `kitty`, `moose`, `milk`, `stegosaurus`, `turkey`, `turtle`, `tux`

좋습니다. 더 좋은 것은 `greetings.csv`를 입력으로 제공할 수 있다면 좋겠습니다.
그러나 파일 시스템에 액세스할 수 없으므로 그렇게 할 수 없습니다.

수정해 봅시다.

#### 1.3.3. 컨테이너 종료

컨테이너를 종료하려면 프롬프트에서 `exit`를 입력하거나 ++ctrl+d++ 키보드 단축키를 사용할 수 있습니다.

```bash
exit
```

이제 프롬프트가 컨테이너를 시작하기 전 상태로 돌아가야 합니다.

#### 1.3.4. 컨테이너에 데이터 마운트

앞서 언급했듯이 컨테이너는 기본적으로 호스트 시스템과 격리되어 있습니다.

컨테이너가 호스트 파일 시스템에 액세스할 수 있도록 하려면 다음 구문을 사용하여 호스트 시스템에서 컨테이너로 **볼륨**을 **마운트**할 수 있습니다:

```bash title="구문"
-v <outside_path>:<inside_path>
```

우리의 경우 `<outside_path>`는 현재 작업 디렉토리이므로 점(`.`)만 사용할 수 있고, `<inside_path>`는 우리가 만든 별칭일 뿐입니다. `/my_project`라고 부르겠습니다(내부 경로는 절대 경로여야 함).

볼륨을 마운트하려면 경로를 교체하고 다음과 같이 docker run 명령에 볼륨 마운트 인수를 추가합니다:

```bash
docker run --rm -it -v .:/my_project 'community.wave.seqera.io/library/cowpy:1.1.5--3db457ae1977a273' /bin/bash
```

이렇게 하면 현재 작업 디렉토리가 컨테이너 내부의 `/my_project` 아래에서 액세스할 수 있는 볼륨으로 마운트됩니다.

`/my_project`의 내용을 나열하여 작동하는지 확인할 수 있습니다:

```bash
ls /my_project
```

??? success "명령 출력"

    ```console
    data               hello-config.nf      hello-modules.nf   hello-world.nf  nextflow.config  solutions         work
    hello-channels.nf  hello-containers.nf  hello-workflow.nf  modules         results          test-params.json
    ```

이제 컨테이너 내부에서 `data/` 아래의 `greetings.csv` 파일을 포함하여 작업 디렉토리의 내용을 볼 수 있습니다.

이것은 효과적으로 파일 시스템의 해당 부분에 액세스하는 데 사용할 수 있는 컨테이너 벽을 통한 터널을 설정했습니다.

#### 1.3.5. 마운트된 데이터 사용

이제 작업 디렉토리를 컨테이너에 마운트했으므로 `cowpy` 명령을 사용하여 `greetings.csv` 파일의 내용을 표시할 수 있습니다.

이를 위해 `cat /my_project/data/greetings.csv | `를 사용하여 CSV 파일의 내용을 `cowpy` 명령으로 파이프합니다.

```bash
cat /my_project/data/greetings.csv | cowpy -c turkey
```

??? success "명령 출력"

    ```console title="data/greetings.csv"
     ____________________
    / Hello,English,123  \
    | Bonjour,French,456 |
    \ Holà,Spanish,789   /
    --------------------
      \                                  ,+*^^*+___+++_
      \                           ,*^^^^              )
        \                       _+*                     ^**+_
        \                    +^       _ _++*+_+++_,         )
                  _+^^*+_    (     ,+*^ ^          \+_        )
                {       )  (    ,(    ,_+--+--,      ^)      ^\
                { (\@)    } f   ,(  ,+-^ __*_*_  ^^\_   ^\       )
              {:;-/    (_+*-+^^^^^+*+*<_ _++_)_    )    )      /
              ( /  (    (        ,___    ^*+_+* )   <    <      \
              U _/     )    *--<  ) ^\-----++__)   )    )       )
                (      )  _(^)^^))  )  )\^^^^^))^*+/    /       /
              (      /  (_))_^)) )  )  ))^^^^^))^^^)__/     +^^
            (     ,/    (^))^))  )  ) ))^^^^^^^))^^)       _)
              *+__+*       (_))^)  ) ) ))^^^^^^))^^^^^)____*^
              \             \_)^)_)) ))^^^^^^^^^^))^^^^)
              (_             ^\__^^^^^^^^^^^^))^^^^^^^)
                ^\___            ^\__^^^^^^))^^^^^^^^)\\
                      ^^^^^\uuu/^^\uuu/^^^^\^\^\^\^\^\^\^\
                        ___) >____) >___   ^\_\_\_\_\_\_\)
                        ^^^//\\_^^//\\_^       ^(\_\_\_\)
                          ^^^ ^^ ^^^ ^
    ```

이것은 예제 인사말을 읊는 칠면조의 ASCII 아트를 생성합니다!
단, 여기서 칠면조는 인사말만이 아닌 전체 행을 반복하고 있습니다.
우리는 Nextflow 워크플로우가 더 잘할 것이라는 것을 이미 알고 있습니다!

이 명령으로 자유롭게 작업해 보십시오.
완료되면 이전처럼 컨테이너를 종료하십시오:

```bash
exit
```

일반 셸로 돌아갑니다.

### 핵심 정리

컨테이너를 풀하고 일회성 또는 대화식으로 실행하는 방법을 알게 되었습니다. 또한 컨테이너 내부에서 데이터에 액세스할 수 있도록 하는 방법을 알게 되었으며, 이를 통해 시스템에 소프트웨어를 설치하지 않고도 관심 있는 도구를 실제 데이터로 시험해 볼 수 있습니다.

### 다음 단계

Nextflow 프로세스 실행에 컨테이너를 사용하는 방법을 배웁니다.

---

## 2. Nextflow에서 컨테이너 사용

Nextflow는 컴퓨팅 환경에 설치되지 않은 도구를 실행할 수 있도록 컨테이너 내부에서 프로세스를 실행하는 기능을 기본적으로 지원합니다.
즉, 프로세스를 실행하는 데 원하는 컨테이너 이미지를 사용할 수 있으며, Nextflow가 이미지 풀, 데이터 마운트 및 그 안에서 프로세스 실행을 처리합니다.

이를 시연하기 위해 개발해 온 파이프라인에 `collectGreetings` 단계 이후에 `cowpy` 단계를 추가할 것입니다.

<figure class="excalidraw">
--8<-- "docs/en/docs/nextflow_run/img/hello-pipeline-cowpy.svg"
</figure>

### 2.1. `cowpy` 모듈 작성

먼저 `cowpy` 프로세스 모듈을 생성합시다.

#### 2.1.1. 새 모듈에 대한 파일 스텁 생성

`cowpy.nf`라는 모듈에 대한 빈 파일을 생성합니다.

```bash
touch modules/cowpy.nf
```

이것은 프로세스 코드를 넣을 장소를 제공합니다.

#### 2.1.2. 모듈 파일에 `cowpy` 프로세스 코드 복사

이전에 작성한 다른 프로세스를 모델로 `cowpy` 프로세스를 만들 수 있습니다.

```groovy title="modules/cowpy.nf" linenums="1"
#!/usr/bin/env nextflow

// cowpy로 ASCII 아트 생성
process cowpy {

    input:
    path input_file
    val character

    output:
    path "cowpy-${input_file}"

    script:
    """
    cat ${input_file} | cowpy -c "${character}" > cowpy-${input_file}
    """

}
```

프로세스는 인사말이 포함된 `input_file`과 `character` 값을 예상합니다.

출력은 `cowpy` 도구에서 생성된 ASCII 아트가 포함된 새 텍스트 파일입니다.

### 2.2. 워크플로우에 cowpy 추가

이제 모듈을 가져오고 프로세스를 호출해야 합니다.

#### 2.2.1. `hello-containers.nf`에 `cowpy` 프로세스 가져오기

워크플로우 블록 위에 import 선언을 삽입하고 적절하게 채웁니다.

=== "수정 후"

    ```groovy title="hello-containers.nf" linenums="3" hl_lines="5"
    // 모듈 포함
    include { sayHello } from './modules/sayHello.nf'
    include { convertToUpper } from './modules/convertToUpper.nf'
    include { collectGreetings } from './modules/collectGreetings.nf'
    include { cowpy } from './modules/cowpy.nf'
    ```

=== "수정 전"

    ```groovy title="hello-containers.nf" linenums="3"
    // 모듈 포함
    include { sayHello } from './modules/sayHello.nf'
    include { convertToUpper } from './modules/convertToUpper.nf'
    include { collectGreetings } from './modules/collectGreetings.nf'
    ```

이제 `cowpy` 모듈을 워크플로우에서 사용할 수 있습니다.

#### 2.2.2. 워크플로우에 `cowpy` 프로세스 호출 추가

기억하시다시피 두 개의 출력을 생성하는 `collectGreetings()` 프로세스의 출력에 `cowpy()` 프로세스를 연결합시다:

- `collectGreetings.out.outfile`에는 출력 파일이 포함됩니다 <--_우리가 원하는 것_
- `collectGreetings.out.report`에는 배치당 인사말 수가 포함된 보고 파일이 포함됩니다

워크플로우 블록에서 다음 코드 변경을 수행하십시오:

=== "수정 후"

    ```groovy title="hello-containers.nf" linenums="19" hl_lines="12-13"
        main:
        // CSV 파일에서 입력용 채널 생성
        greeting_ch = channel.fromPath(params.input)
                            .splitCsv()
                            .map { line -> line[0] }
        // 인사말을 내보냅니다
        sayHello(greeting_ch)
        // 인사말을 대문자로 변환
        convertToUpper(sayHello.out)
        // 모든 인사말을 하나의 파일에 수집
        collectGreetings(convertToUpper.out.collect(), params.batch)
        // cowpy로 인사말의 ASCII 아트 생성
        cowpy(collectGreetings.out.outfile, params.character)
    ```

=== "수정 전"

    ```groovy title="hello-containers.nf" linenums="19"
        main:
        // CSV 파일에서 입력용 채널 생성
        greeting_ch = channel.fromPath(params.input)
                            .splitCsv()
                            .map { line -> line[0] }
        // 인사말을 내보냅니다
        sayHello(greeting_ch)
        // 인사말을 대문자로 변환
        convertToUpper(sayHello.out)
        // 모든 인사말을 하나의 파일에 수집
        collectGreetings(convertToUpper.out.collect(), params.batch)
    ```

인사말을 말할 캐릭터를 지정하기 위해 새 CLI 매개변수 `params.character`를 선언했습니다.

#### 2.2.3. `params` 블록에 `character` 매개변수 추가

이것은 기술적으로 선택 사항이지만 권장되는 관행이며, 그러는 동안 캐릭터의 기본값을 설정할 기회입니다.

=== "수정 후"

    ```groovy title="hello-containers.nf" linenums="9" hl_lines="7"
    /*
    * 파이프라인 매개변수
    */
    params {
        input: Path = 'data/greetings.csv'
        batch: String = 'batch'
        character: String = 'turkey'
    }
    ```

=== "수정 전"

    ```groovy title="hello-containers.nf" linenums="9"
    /*
    * 파이프라인 매개변수
    */
    params {
        input: Path = 'data/greetings.csv'
        batch: String = 'batch'
    }
    ```

이제 게으르게 명령줄에서 캐릭터 매개변수 입력을 건너뛸 수 있습니다.

#### 2.2.4. 워크플로우 출력 업데이트

`cowpy` 프로세스의 출력을 게시하도록 워크플로우 출력을 업데이트해야 합니다.

##### 2.2.4.1. `publish:` 섹션 업데이트

`workflow 블록`에서 다음 코드 변경을 수행하십시오:

=== "수정 후"

    ```groovy title="hello-containers.nf" linenums="34" hl_lines="6"
        publish:
        first_output = sayHello.out
        uppercased = convertToUpper.out
        collected = collectGreetings.out.outfile
        batch_report = collectGreetings.out.report
        cowpy_art = cowpy.out
    ```

=== "수정 전"

    ```groovy title="hello-containers.nf" linenums="34"
        publish:
        first_output = sayHello.out
        uppercased = convertToUpper.out
        collected = collectGreetings.out.outfile
        batch_report = collectGreetings.out.report
    ```

`cowpy` 프로세스는 하나의 출력만 생성하므로 `.out`을 추가하여 일반적인 방식으로 참조할 수 있습니다.

하지만 지금은 워크플로우 수준 출력 업데이트를 마무리합시다.

##### 2.2.4.2. `output` 블록 업데이트

`output` 블록에 최종 `cowpy_art` 출력을 추가해야 합니다. 그러는 동안 파이프라인이 완료되어 실제로 중요한 출력이 무엇인지 알게 되었으므로 게시 대상도 편집합시다.

`output` 블록에서 다음 코드 변경을 수행하십시오:

=== "수정 후"

    ```groovy title="hello-containers.nf" linenums="42" hl_lines="3 7 11 15 18-21"
    output {
        first_output {
            path 'hello_containers/intermediates'
            mode 'copy'
        }
        uppercased {
            path 'hello_containers/intermediates'
            mode 'copy'
        }
        collected {
            path 'hello_containers/intermediates'
            mode 'copy'
        }
        batch_report {
            path 'hello_containers'
            mode 'copy'
        }
        cowpy_art {
            path 'hello_containers'
            mode 'copy'
        }
    }
    ```

=== "수정 전"

    ```groovy title="hello-containers.nf" linenums="42" hl_lines="3 7 11 15"
    output {
        first_output {
            path 'hello_containers'
            mode 'copy'
        }
        uppercased {
            path 'hello_containers'
            mode 'copy'
        }
        collected {
            path 'hello_containers'
            mode 'copy'
        }
        batch_report {
            path 'hello_containers'
            mode 'copy'
        }
    }
    ```

이제 게시된 출력이 약간 더 정리됩니다.

#### 2.2.5. 워크플로우 실행

요약하자면, 우리가 목표로 하는 것은 다음과 같습니다:

<figure class="excalidraw">
--8<-- "docs/en/docs/hello_nextflow/img/hello_pipeline_complete.svg"
</figure>

작동할 것 같습니까?

이전에 게시된 출력을 삭제하여 깨끗한 상태로 만들고 `-resume` 플래그로 워크플로우를 실행합시다.

```bash
rm -r hello_containers/
nextflow run hello-containers.nf -resume
```

??? failure "명령 출력 (명확성을 위해 편집됨)"

    ```console hl_lines="10 13 20-21 26-27"
     N E X T F L O W   ~  version 25.10.2

    Launching `hello-containers.nf` [lonely_woese] DSL2 - revision: abf1dccf7f

    executor >  local (1)
    [c9/f5c686] sayHello (3)       [100%] 3 of 3, cached: 3 ✔
    [ef/3135a8] convertToUpper (3) [100%] 3 of 3, cached: 3 ✔
    [7f/f435e3] collectGreetings   [100%] 1 of 1, cached: 1 ✔
    [9b/02e776] cowpy              [  0%] 0 of 1 ✘
    ERROR ~ Error executing process > 'cowpy'

    Caused by:
      Process `cowpy` terminated with an error exit status (127)


    Command executed:

      cat COLLECTED-batch-output.txt | cowpy -c "turkey" > cowpy-COLLECTED-batch-output.txt

    Command exit status:
      127

    Command output:
      (empty)

    Command error:
      .command.sh: line 2: cowpy: command not found

    Work dir:
      /workspaces/training/hello-nextflow/work/9b/02e7761db848f82db3c3e59ff3a9b6

    Tip: when you have fixed the problem you can continue the execution adding the option `-resume` to the run command line

    -- Check '.nextflow.log' file for details
    ERROR ~ Cannot access first() element from an empty List

    -- Check '.nextflow.log' file for details
    ```

오류가 발생했습니다!
`error exit status (127)`에서 제공된 오류 코드는 요청한 실행 파일을 찾을 수 없음을 의미합니다.

`cowpy` 도구를 호출하고 있지만 아직 컨테이너를 지정하지 않았기 때문에(이런) 당연히 그렇습니다.

### 2.3. 컨테이너를 사용하여 `cowpy` 프로세스 실행

컨테이너를 지정하고 Nextflow에게 `cowpy()` 프로세스에 사용하도록 지시해야 합니다.

#### 2.3.1. `cowpy`에 대한 컨테이너 지정

이 튜토리얼의 첫 번째 섹션에서 직접 사용했던 동일한 이미지를 사용할 수 있습니다.

다음과 같이 프로세스 정의에 `container` 지시문을 추가하도록 `cowpy.nf` 모듈을 편집하십시오:

=== "수정 후"

    ```groovy title="modules/cowpy.nf" linenums="4" hl_lines="3"
    process cowpy {

        container 'community.wave.seqera.io/library/cowpy:1.1.5--3db457ae1977a273'

        input:
        path input_file
        val character

        output:
        path "cowpy-${input_file}"

        script:
        """
        cat ${input_file} | cowpy -c "${character}" > cowpy-${input_file}
        """
    }
    ```

=== "수정 전"

    ```groovy title="modules/cowpy.nf" linenums="4"
    process cowpy {

        input:
        path input_file
        val character

        output:
        path "cowpy-${input_file}"

        script:
        """
        cat ${input_file} | cowpy -c "${character}" > cowpy-${input_file}
        """
    }
    ```

이것은 _Docker 사용이 활성화된 경우_, Nextflow가 여기에 지정된 컨테이너 이미지를 사용하여 프로세스를 실행해야 함을 알려줍니다.

#### 2.3.2. `nextflow.config` 파일을 통해 Docker 사용 활성화

*'Docker 사용이 활성화된 경우'*라고 말했다는 것을 주목하십시오. 기본적으로 활성화되어 있지 않으므로 Nextflow에게 Docker를 사용할 수 있도록 허용해야 합니다.
이를 위해 구성을 다루는 이 과정의 다음이자 마지막 파트(파트 6)의 주제를 약간 앞당기겠습니다.

Nextflow가 워크플로우 실행을 구성하기 위해 제공하는 주요 방법 중 하나는 `nextflow.config` 파일을 사용하는 것입니다.
현재 디렉토리에 이러한 파일이 있으면 Nextflow가 자동으로 로드하고 포함된 모든 구성을 적용합니다.

Docker를 명시적으로 비활성화하는 한 줄의 코드가 포함된 `nextflow.config` 파일을 제공했습니다: `docker.enabled = false`.

이제 Docker를 활성화하기 위해 `true`로 전환합시다:

=== "수정 후"

    ```console title="nextflow.config" linenums="1" hl_lines="1"
    docker.enabled = true
    ```

=== "수정 전"

    ```console title="nextflow.config" linenums="1" hl_lines="1"
    docker.enabled = false
    ```

!!! tip "팁"

    `-with-docker <container>` 매개변수를 사용하여 명령줄에서 실행별로 Docker 실행을 활성화할 수 있습니다.
    그러나 이것은 전체 워크플로우에 대해 하나의 컨테이너만 지정할 수 있는 반면, 방금 보여드린 접근 방식은 프로세스당 다른 컨테이너를 지정할 수 있습니다.
    이것은 모듈성, 코드 유지보수 및 재현성에 더 좋습니다.

#### 2.3.3. Docker가 활성화된 상태로 워크플로우 실행

`-resume` 플래그로 워크플로우를 실행하십시오:

```bash
nextflow run hello-containers.nf -resume
```

??? success "명령 출력"

    ```console
     N E X T F L O W   ~  version 25.10.2

    Launching `hello-containers.nf` [drunk_perlman] DSL2 - revision: abf1dccf7f

    executor >  local (1)
    [c9/f5c686] sayHello (3)       [100%] 3 of 3, cached: 3 ✔
    [ef/3135a8] convertToUpper (3) [100%] 3 of 3, cached: 3 ✔
    [7f/f435e3] collectGreetings   [100%] 1 of 1, cached: 1 ✔
    [98/656c6c] cowpy              [100%] 1 of 1 ✔
    ```

이번에는 정말로 작동합니다!
평소와 같이 해당 results 디렉토리에서 워크플로우 출력을 찾을 수 있지만, 이번에는 보고서와 최종 출력만 최상위 수준에 있고 모든 중간 파일은 하위 디렉토리로 이동되어 약간 더 깔끔하게 정리되어 있습니다.

??? abstract "디렉토리 내용"

    ```console
    results/hello_containers/
    ├── cowpy-COLLECTED-batch-output.txt
    ├── intermediates
    │   ├── Bonjour-output.txt
    │   ├── COLLECTED-batch-output.txt
    │   ├── Hello-output.txt
    │   ├── Holà-output.txt
    │   ├── UPPER-Bonjour-output.txt
    │   ├── UPPER-Hello-output.txt
    │   └── UPPER-Holà-output.txt
    └── batch-report.txt
    ```

최종 ASCII 아트 출력은 `results/hello_containers/` 디렉토리에 `cowpy-COLLECTED-batch-output.txt`라는 이름으로 있습니다.

??? abstract "파일 내용"

    ```console title="results/hello_containers/cowpy-COLLECTED-batch-output.txt"
    _________
    / HOLà    \
    | HELLO   |
    \ BONJOUR /
    ---------
      \                                  ,+*^^*+___+++_
      \                           ,*^^^^              )
        \                       _+*                     ^**+_
        \                    +^       _ _++*+_+++_,         )
                  _+^^*+_    (     ,+*^ ^          \+_        )
                {       )  (    ,(    ,_+--+--,      ^)      ^\
                { (\@)    } f   ,(  ,+-^ __*_*_  ^^\_   ^\       )
              {:;-/    (_+*-+^^^^^+*+*<_ _++_)_    )    )      /
              ( /  (    (        ,___    ^*+_+* )   <    <      \
              U _/     )    *--<  ) ^\-----++__)   )    )       )
                (      )  _(^)^^))  )  )\^^^^^))^*+/    /       /
              (      /  (_))_^)) )  )  ))^^^^^))^^^)__/     +^^
            (     ,/    (^))^))  )  ) ))^^^^^^^))^^)       _)
              *+__+*       (_))^)  ) ) ))^^^^^^))^^^^^)____*^
              \             \_)^)_)) ))^^^^^^^^^^))^^^^)
              (_             ^\__^^^^^^^^^^^^))^^^^^^^)
                ^\___            ^\__^^^^^^))^^^^^^^^)\\
                      ^^^^^\uuu/^^\uuu/^^^^\^\^\^\^\^\^\^\
                        ___) >____) >___   ^\_\_\_\_\_\_\)
                        ^^^//\\_^^//\\_^       ^(\_\_\_\)
                          ^^^ ^^ ^^^ ^
    ```

우리의 아름다운 칠면조가 원하는 대로 인사말을 말하고 있습니다.

#### 2.3.4. Nextflow가 컨테이너화된 작업을 시작한 방법 검사

이 섹션의 마지막 코다로, Nextflow가 내부적으로 컨테이너와 어떻게 작동하는지에 대해 조금 더 통찰력을 얻기 위해 `cowpy` 프로세스 호출 중 하나의 작업 하위 디렉토리를 살펴봅시다.

`nextflow run` 명령의 출력을 확인하여 `cowpy` 프로세스의 작업 하위 디렉토리 경로를 찾으십시오.
위에 표시된 실행에서 얻은 것을 보면, `cowpy` 프로세스의 콘솔 로그 라인은 `[98/656c6c]`로 시작합니다.
이것은 다음 잘린 디렉토리 경로에 해당합니다: `work/98/656c6c`.

해당 디렉토리에서 파이프라인을 실행하는 과정에서 Nextflow가 대신 실행한 모든 명령이 포함된 `.command.run` 파일을 찾을 수 있습니다.

??? abstract "파일 내용"

    ```console title="work/98/656c6c90cce1667c094d880f4b6dcc/.command.run"
    #!/bin/bash
    ### ---
    ### name: 'cowpy'
    ### container: 'community.wave.seqera.io/library/cowpy:1.1.5--3db457ae1977a273'
    ### outputs:
    ### - 'cowpy-COLLECTED-batch-output.txt'
    ### ...
    set -e
    set -u
    NXF_DEBUG=${NXF_DEBUG:=0}; [[ $NXF_DEBUG > 1 ]] && set -x
    NXF_ENTRY=${1:-nxf_main}


    nxf_sleep() {
      sleep $1 2>/dev/null || sleep 1;
    }

    nxf_date() {
        local ts=$(date +%s%3N);
        if [[ ${#ts} == 10 ]]; then echo ${ts}000
        elif [[ $ts == *%3N ]]; then echo ${ts/\%3N/000}
        elif [[ $ts == *3N ]]; then echo ${ts/3N/000}
        elif [[ ${#ts} == 13 ]]; then echo $ts
        else echo "Unexpected timestamp value: $ts"; exit 1
        fi
    }

    nxf_env() {
        echo '============= task environment ============='
        env | sort | sed "s/\(.*\)AWS\(.*\)=\(.\{6\}\).*/\1AWS\2=\3xxxxxxxxxxxxx/"
        echo '============= task output =================='
    }

    nxf_kill() {
        declare -a children
        while read P PP;do
            children[$PP]+=" $P"
        done < <(ps -e -o pid= -o ppid=)

        kill_all() {
            [[ $1 != $$ ]] && kill $1 2>/dev/null || true
            for i in ${children[$1]:=}; do kill_all $i; done
        }

        kill_all $1
    }

    nxf_mktemp() {
        local base=${1:-/tmp}
        mkdir -p "$base"
        if [[ $(uname) = Darwin ]]; then mktemp -d $base/nxf.XXXXXXXXXX
        else TMPDIR="$base" mktemp -d -t nxf.XXXXXXXXXX
        fi
    }

    nxf_fs_copy() {
      local source=$1
      local target=$2
      local basedir=$(dirname $1)
      mkdir -p $target/$basedir
      cp -fRL $source $target/$basedir
    }

    nxf_fs_move() {
      local source=$1
      local target=$2
      local basedir=$(dirname $1)
      mkdir -p $target/$basedir
      mv -f $source $target/$basedir
    }

    nxf_fs_rsync() {
      rsync -rRl $1 $2
    }

    nxf_fs_rclone() {
      rclone copyto $1 $2/$1
    }

    nxf_fs_fcp() {
      fcp $1 $2/$1
    }

    on_exit() {
        local last_err=$?
        local exit_status=${nxf_main_ret:=0}
        [[ ${exit_status} -eq 0 && ${nxf_unstage_ret:=0} -ne 0 ]] && exit_status=${nxf_unstage_ret:=0}
        [[ ${exit_status} -eq 0 && ${last_err} -ne 0 ]] && exit_status=${last_err}
        printf -- $exit_status > /workspaces/training/hello-nextflow/work/98/656c6c90cce1667c094d880f4b6dcc/.exitcode
        set +u
        docker rm $NXF_BOXID &>/dev/null || true
        exit $exit_status
    }

    on_term() {
        set +e
        docker stop $NXF_BOXID
    }

    nxf_launch() {
        docker run -i --cpu-shares 1024 -e "NXF_TASK_WORKDIR" -v /workspaces/training/hello-nextflow/work:/workspaces/training/hello-nextflow/work -w "$NXF_TASK_WORKDIR" --name $NXF_BOXID community.wave.seqera.io/library/cowpy:1.1.5--3db457ae1977a273 /bin/bash -ue /workspaces/training/hello-nextflow/work/98/656c6c90cce1667c094d880f4b6dcc/.command.sh
    }

    nxf_stage() {
        true
        # stage input files
        rm -f COLLECTED-batch-output.txt
        ln -s /workspaces/training/hello-nextflow/work/7f/f435e3f2cf95979b5f3d7647ae6696/COLLECTED-batch-output.txt COLLECTED-batch-output.txt
    }

    nxf_unstage_outputs() {
        true
    }

    nxf_unstage_controls() {
        true
    }

    nxf_unstage() {
        if [[ ${nxf_main_ret:=0} == 0 ]]; then
            (set -e -o pipefail; (nxf_unstage_outputs | tee -a .command.out) 3>&1 1>&2 2>&3 | tee -a .command.err)
            nxf_unstage_ret=$?
        fi
        nxf_unstage_controls
    }

    nxf_main() {
        trap on_exit EXIT
        trap on_term TERM INT USR2
        trap '' USR1

        [[ "${NXF_CHDIR:-}" ]] && cd "$NXF_CHDIR"
        export NXF_BOXID="nxf-$(dd bs=18 count=1 if=/dev/urandom 2>/dev/null | base64 | tr +/ 0A | tr -d '\r\n')"
        NXF_SCRATCH=''
        [[ $NXF_DEBUG > 0 ]] && nxf_env
        touch /workspaces/training/hello-nextflow/work/98/656c6c90cce1667c094d880f4b6dcc/.command.begin
        set +u
        set -u
        [[ $NXF_SCRATCH ]] && cd $NXF_SCRATCH
        export NXF_TASK_WORKDIR="$PWD"
        nxf_stage

        set +e
        (set -o pipefail; (nxf_launch | tee .command.out) 3>&1 1>&2 2>&3 | tee .command.err) &
        pid=$!
        wait $pid || nxf_main_ret=$?
        nxf_unstage
    }

    $NXF_ENTRY

    ```

이 파일에서 `nxf_launch`를 검색하면 다음과 같은 것을 볼 수 있습니다:

```console
nxf_launch() {
    docker run -i --cpu-shares 1024 -e "NXF_TASK_WORKDIR" -v /workspaces/training/hello-nextflow/work:/workspaces/training/hello-nextflow/work -w "$NXF_TASK_WORKDIR" --name $NXF_BOXID community.wave.seqera.io/library/cowpy:1.1.5--3db457ae1977a273 /bin/bash -ue /workspaces/training/hello-nextflow/work/98/656c6c90cce1667c094d880f4b6dcc/.command.sh
}
```

보시다시피, Nextflow는 `docker run` 명령을 사용하여 프로세스 호출을 시작합니다.
또한 해당 작업 하위 디렉토리를 컨테이너에 마운트하고, 컨테이너 내부의 작업 디렉토리를 적절하게 설정하고, `.command.sh` 파일에서 템플릿된 bash 스크립트를 실행합니다.

첫 번째 섹션에서 수동으로 해야 했던 모든 어려운 작업? Nextflow가 백그라운드에서 우리를 위해 해줍니다!

```txt
 _______________________
< Hurray for robots...! >
 -----------------------
                                   ,-----.
                                   |     |
                                ,--|     |-.
                         __,----|  |     | |
                       ,;::     |  `_____' |
                       `._______|    i^i   |
                                `----| |---'| .
                           ,-------._| |== ||//
                           |       |_|P`.  /'/
                           `-------' 'Y Y/'/'
                                     .==\ /_\
   ^__^                             /   /'|  `i
   (oo)\_______                   /'   /  |   |
   (__)\       )\/\             /'    /   |   `i
       ||----w |           ___,;`----'.___L_,-'`\__
       ||     ||          i_____;----\.____i""\____\
```

### 핵심 정리

프로세스를 실행하기 위해 Nextflow에서 컨테이너를 사용하는 방법을 알게 되었습니다.

### 다음 단계

휴식을 취하십시오!

준비가 되면 [**파트 6: Hello Config**](./06_hello_config.md)로 이동하여 인프라에 맞게 파이프라인 실행을 구성하고 입력 및 매개변수 구성을 관리하는 방법을 배우십시오.

마지막 파트이며, 이것을 완료하면 이 과정을 마치게 됩니다!

---

## 퀴즈

<quiz>
컨테이너란 무엇입니까?
- [ ] 가상 머신의 일종
- [ ] 파일 압축 형식
- [x] 애플리케이션 실행에 필요한 모든 것이 포함된 경량의 단독 실행형 단위
- [ ] 네트워크 프로토콜
</quiz>

<quiz>
컨테이너 이미지와 컨테이너 인스턴스의 차이점은 무엇입니까?
- [ ] 동일한 것임
- [x] 이미지는 템플릿이고 인스턴스는 해당 이미지에서 생성된 실행 중인 컨테이너임
- [ ] 인스턴스는 템플릿이고 이미지는 실행 중인 컨테이너임
- [ ] 이미지는 Docker용이고 인스턴스는 Singularity용임
</quiz>

<quiz>
`docker run` 명령에서 `-v` 플래그는 무엇을 합니까?
- [ ] 자세한 출력 활성화
- [ ] 컨테이너 유효성 검사
- [x] 호스트 시스템에서 컨테이너로 볼륨 마운트
- [ ] 컨테이너 버전 지정

자세히 알아보기: [1.3.4. 컨테이너에 데이터 마운트](#134-컨테이너에-데이터-마운트)
</quiz>

<quiz>
컨테이너를 사용할 때 볼륨을 마운트해야 하는 이유는 무엇입니까?
- [ ] 컨테이너 성능 향상
- [ ] 디스크 공간 절약
- [x] 컨테이너가 기본적으로 호스트 파일 시스템과 격리되어 있기 때문
- [ ] 네트워킹 활성화

자세히 알아보기: [1.3.4. 컨테이너에 데이터 마운트](#134-컨테이너에-데이터-마운트)
</quiz>

<quiz>
Nextflow 프로세스에 컨테이너를 어떻게 지정합니까?
- [ ] `docker 'container-uri'`
- [ ] `image 'container-uri'`
- [x] `container 'container-uri'`
- [ ] `use 'container-uri'`

자세히 알아보기: [2.3.1. cowpy에 대한 컨테이너 지정](#231-cowpy에-대한-컨테이너-지정)
</quiz>

<quiz>
어떤 `nextflow.config` 설정이 워크플로우에 대해 Docker를 활성화합니까?
- [ ] `#!groovy process.docker = true`
- [x] `#!groovy docker.enabled = true`
- [ ] `#!groovy container.engine = 'docker'`
- [ ] `#!groovy docker.activate = true`

자세히 알아보기: [2.3.2. `nextflow.config` 파일을 통해 Docker 사용 활성화](#232-nextflowconfig-파일을-통해-docker-사용-활성화)
</quiz>

<quiz>
컨테이너에서 프로세스를 실행할 때 Nextflow가 자동으로 처리하는 것은 무엇입니까? (해당되는 것 모두 선택)
- [x] 필요한 경우 컨테이너 이미지 풀
- [x] 작업 디렉토리 마운트
- [x] 컨테이너 내부에서 프로세스 스크립트 실행
- [x] 실행 후 컨테이너 인스턴스 정리

자세히 알아보기: [2.3.4. Nextflow가 컨테이너화된 작업을 시작한 방법 검사](#234-nextflow가-컨테이너화된-작업을-시작한-방법-검사)
</quiz>
