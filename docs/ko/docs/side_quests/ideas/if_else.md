# 파트 2: If - Else

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } AI 지원 번역 - [자세히 알아보기 및 개선 제안](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

[TODO]

---

## 1. 소를 유명한 과학자들의 인용구를 말하게 만들기

이 섹션에는 지금까지 배운 내용을 연습할 수 있는 심화 연습 문제가 포함되어 있습니다.
이러한 연습 문제를 수행하는 것은 교육의 이후 부분을 이해하는 데 _필수는 아니지만_, 소가 유명한 과학자들의 말을 인용하도록 만드는 방법을 알아내면서 학습한 내용을 강화하는 재미있는 방법을 제공합니다.

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

### 1.1. `hello-containers.nf` 스크립트를 수정하여 getQuote 프로세스 사용하기

`containers/data/pioneers.csv` 파일에 컴퓨터 및 생물학 분야 선구자들의 목록이 있습니다.
높은 수준에서 이 연습 문제를 완료하려면 다음이 필요합니다:

- 기본 `params.input_file`을 `pioneers.csv` 파일을 가리키도록 수정합니다.
- `quote` 컨테이너를 사용하여 각 입력에 대한 인용구를 가져오는 `getQuote` 프로세스를 생성합니다.
- `getQuote` 프로세스의 출력을 `cowsay` 프로세스에 연결하여 인용구를 표시합니다.

`quote` 컨테이너 이미지의 경우, 이전 심화 연습 문제에서 직접 빌드한 것을 사용하거나 Seqera Containers에서 가져온 것을 사용할 수 있습니다.

!!! Hint

    getQuote 프로세스의 `script` 블록에 적합한 선택은 다음과 같습니다:
        ```groovy
        script:
            def safe_author = author.tokenize(' ').join('-')
            """
            quote "$author" > quote-${safe_author}.txt
            echo "-${author}" >> quote-${safe_author}.txt
            """
        ```

이 연습 문제의 해결책은 `containers/solutions/hello-containers-4.1.nf`에서 찾을 수 있습니다.

### 1.2. Nextflow 파이프라인을 수정하여 `quote` 및 `sayHello` 모드에서 실행되도록 하기

파이프라인에 분기 로직을 추가하여 `quote`와 `sayHello` 모두를 위한 입력을 받을 수 있도록 합니다.
다음은 Nextflow 워크플로우에서 `if` 문을 사용하는 방법의 예시입니다:

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

!!! Hint

    `new_ch = processName.out`을 사용하여 프로세스의 출력 채널에 이름을 할당할 수 있습니다.

이 연습 문제의 해결책은 `containers/solutions/hello-containers-4.2.nf`에서 찾을 수 있습니다.

### 요점

Nextflow에서 컨테이너를 사용하여 프로세스를 실행하는 방법과 파이프라인에 분기 로직을 구축하는 방법을 알게 되었습니다!

### 다음 단계

축하하고, 스트레칭 휴식을 취하고 물을 마시십시오!

준비가 되면, 지금까지 배운 내용을 보다 현실적인 데이터 분석 사용 사례에 적용하는 방법을 배우기 위해 이 교육 시리즈의 Part 3로 이동하십시오.
