# 파트 1: 컨테이너 더 알아보기

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } AI 지원 번역 - [자세히 알아보기 및 개선 제안](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

[TODO]

---

## 1. 컨테이너 이미지를 찾거나 만드는 방법

일부 소프트웨어 개발자는 Docker Hub와 같은 컨테이너 레지스트리에서 사용할 수 있는 소프트웨어용 컨테이너 이미지를 제공하지만, 많은 경우 제공하지 않습니다.
이 선택적 섹션에서는 Nextflow 파이프라인에서 사용하려는 도구의 컨테이너 이미지를 얻는 두 가지 방법을 보여드리겠습니다: Seqera Containers를 사용하는 방법과 컨테이너 이미지를 직접 빌드하는 방법입니다.

이 섹션의 마지막 연습에서 사용할 `quote` pip 패키지의 컨테이너 이미지를 얻거나 빌드하게 됩니다.

### 1.1. Seqera Containers에서 컨테이너 이미지 얻기

Seqera Containers는 pip 및 conda(bioconda 포함)로 설치 가능한 도구의 컨테이너 이미지를 빌드하는 무료 서비스입니다.
[Seqera Containers](https://www.seqera.io/containers/)로 이동하여 `quote` pip 패키지를 검색하세요.

![Seqera Containers](img/seqera-containers-1.png)

"+Add"를 클릭한 다음 "Get Container"를 클릭하여 `quote` pip 패키지의 컨테이너 이미지를 요청하세요.

![Seqera Containers](img/seqera-containers-2.png)

이 버전의 패키지에 대해 커뮤니티 컨테이너가 처음 빌드되는 경우, 완료하는 데 몇 분이 걸릴 수 있습니다.
생성된 컨테이너 이미지의 URI(예: `community.wave.seqera.io/library/pip_quote:ae07804021465ee9`)를 복사하려면 클릭하세요.

이제 컨테이너 이미지를 사용하여 `quote` 명령을 실행하고 Grace Hopper의 무작위 명언을 얻을 수 있습니다.

```bash
docker run --rm community.wave.seqera.io/library/pip_quote:ae07804021465ee9 quote "Grace Hopper"
```

출력:

```console title="Output"
Humans are allergic to change. They love to say, 'We've always done it
this way.' I try to fight that. That's why I have a clock on my wall
that runs counter-clockwise.
```

### 1.2. 컨테이너 이미지 직접 빌드하기

Seqera Containers 웹사이트의 빌드 세부 정보를 사용하여 `quote` pip 패키지의 컨테이너 이미지를 직접 빌드해 보겠습니다.
Seqera Containers 웹사이트로 돌아가서 "Build Details" 버튼을 클릭하세요.

먼저 살펴볼 항목은 `Dockerfile`입니다. 이는 컨테이너 이미지를 빌드하는 데 필요한 모든 명령이 포함된 스크립트 파일의 한 유형입니다.
각 부분이 무엇을 하는지 이해하는 데 도움이 되도록 아래 Dockerfile에 설명 주석을 추가했습니다.

```Dockerfile title="Dockerfile"
# micromamba 베이스 도커 이미지에서 시작
FROM mambaorg/micromamba:1.5.10-noble
# conda.yml 파일을 컨테이너로 복사
COPY --chown=$MAMBA_USER:$MAMBA_USER conda.yml /tmp/conda.yml
# Nextflow가 사용할 다양한 유틸리티와 conda.yml 파일의 패키지를 설치
RUN micromamba install -y -n base -f /tmp/conda.yml \
    && micromamba install -y -n base conda-forge::procps-ng \
    && micromamba env export --name base --explicit > environment.lock \
    && echo ">> CONDA_LOCK_START" \
    && cat environment.lock \
    && echo "<< CONDA_LOCK_END" \
    && micromamba clean -a -y
# 컨테이너를 root 사용자로 실행
USER root
# PATH 환경 변수를 micromamba 설치 디렉토리를 포함하도록 설정
ENV PATH="$MAMBA_ROOT_PREFIX/bin:$PATH"
```

두 번째로 살펴볼 항목은 `conda.yml` 파일입니다. 이 파일에는 컨테이너 이미지에 설치해야 하는 패키지 목록이 포함되어 있습니다.

```conda.yml title="conda.yml"
channels:
- conda-forge
- bioconda
dependencies:
- pip
- pip:
  - quote==3.0.0 #
```

이 파일들의 내용을 `containers/build` 디렉토리에 있는 스텁에 복사한 다음, 다음 명령을 실행하여 컨테이너 이미지를 직접 빌드하세요.

!!! Note "참고"

    `-t quote:latest` 플래그를 사용하여 컨테이너 이미지에 `quote`라는 이름과 `latest` 태그를 지정합니다.
    이 시스템에서 컨테이너 이미지를 실행할 때 이 태그를 사용하여 참조할 수 있습니다.

```bash
docker build -t quote:latest containers/build
```

빌드가 완료되면 방금 빌드한 컨테이너 이미지를 실행할 수 있습니다.

```bash
docker run --rm quote:latest quote "Margaret Oakley Dayhoff"
```

### 핵심 정리

Nextflow 파이프라인에서 사용하려는 도구의 컨테이너 이미지를 얻는 두 가지 방법을 학습했습니다: Seqera Containers를 사용하는 방법과 컨테이너 이미지를 직접 빌드하는 방법입니다.

### 다음 단계

이 교육 시리즈의 [다음 챕터](./04_hello_genomics.md)로 계속 진행하는 데 필요한 모든 것을 갖추었습니다.
또한 `quote` 컨테이너를 사용하여 컴퓨터/생물학 선구자들의 명언을 가져오고 `cowsay` 컨테이너를 사용하여 출력하는 선택적 연습을 계속할 수도 있습니다.

---

## 2. 소가 유명한 과학자들을 인용하게 만들기

이 섹션에는 지금까지 학습한 내용을 연습할 수 있는 심화 연습이 포함되어 있습니다.
이러한 연습을 수행하는 것은 교육의 후반부를 이해하는 데 _필수는 아니지만_, 소가 유명한 과학자들을 인용하게 만드는 방법을 알아내면서 학습한 내용을 강화하는 재미있는 방법을 제공합니다.

```console title="cowsay-output-Grace-Hopper.txt"
  _________________________________________________
 /                                                 \
| Humans are allergic to change. They love to       |
| say, 'We've always done it this way.' I try to fi |
| ght that. That's why I have a clock on my wall th |
| at runs counter-clockwise.                        |
| -Grace Hopper                                     |
 \                                                 /
  =================================================
                                                 \
                                                  \
                                                    ^__^
                                                    (oo)\_______
                                                    (__)\       )\/\
                                                        ||----w |
                                                        ||     ||
```

### 2.1. getQuote 프로세스를 사용하도록 `hello-containers.nf` 스크립트 수정하기

`containers/data/pioneers.csv` 파일에 컴퓨터 및 생물학 선구자 목록이 있습니다.
이 연습을 완료하려면 높은 수준에서 다음을 수행해야 합니다:

- 기본 `params.input_file`을 `pioneers.csv` 파일을 가리키도록 수정합니다.
- `quote` 컨테이너를 사용하여 각 입력에 대한 명언을 가져오는 `getQuote` 프로세스를 생성합니다.
- `getQuote` 프로세스의 출력을 `cowsay` 프로세스에 연결하여 명언을 표시합니다.

`quote` 컨테이너 이미지의 경우, 이전 심화 연습에서 직접 빌드한 것을 사용하거나 Seqera Containers에서 얻은 것을 사용할 수 있습니다.

!!! Hint "힌트"

    getQuote 프로세스의 `script` 블록에 적합한 선택은 다음과 같습니다:
        ```groovy
        script:
            def safe_author = author.tokenize(' ').join('-')
            """
            quote "$author" > quote-${safe_author}.txt
            echo "-${author}" >> quote-${safe_author}.txt
            """
        ```

이 연습의 해결책은 `containers/solutions/hello-containers-4.1.nf`에서 찾을 수 있습니다.

### 2.2. `quote` 및 `sayHello` 모드에서 실행할 수 있도록 Nextflow 파이프라인 수정하기

파이프라인에 분기 로직을 추가하여 `quote`와 `sayHello` 모두를 위한 입력을 받을 수 있도록 하세요.
다음은 Nextflow 워크플로우에서 `if` 문을 사용하는 방법의 예입니다:

```groovy title="hello-containers.nf"
workflow {
    if (params.quote) {
        ...
    }
    else {
        ...
    }
    cowSay(text_ch)
}
```

!!! Hint "힌트"

    `new_ch = processName.out`을 사용하여 프로세스의 출력 채널에 이름을 할당할 수 있습니다.

이 연습의 해결책은 `containers/solutions/hello-containers-4.2.nf`에서 찾을 수 있습니다.

### 핵심 정리

Nextflow에서 컨테이너를 사용하여 프로세스를 실행하는 방법과 파이프라인에 분기 로직을 구축하는 방법을 알게 되었습니다!

### 다음 단계

축하하고, 스트레칭 휴식을 취하고 물을 마시세요!

준비가 되면 이 교육 시리즈의 파트 3으로 이동하여 지금까지 학습한 내용을 보다 현실적인 데이터 분석 사용 사례에 적용하는 방법을 학습하세요.
