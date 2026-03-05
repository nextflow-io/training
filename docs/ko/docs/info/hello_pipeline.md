---
title: Hello pipeline
description: Hello pipeline이 무엇을 하고 어떻게 구성되어 있는지에 대한 요약.
hide:
  - toc
  - footer
---

# Hello pipeline

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } AI 지원 번역 - [자세히 알아보기 및 개선 제안](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

대부분의 교육 과정에서는 Nextflow 개념과 메커니즘을 시연하기 위해 간단한 도메인 비의존적 pipeline을 사용합니다.
Hello Nextflow 과정에서는 모든 설계 및 구현 결정을 설명하면서 이 pipeline을 단계별로 개발하는 방법을 보여줍니다.
다른 교육에서는 이 pipeline 또는 그 일부를 시작점으로 사용합니다.

이 페이지는 Hello Nextflow 과정 완료 시점의 pipeline 상태를 요약합니다.

### 요약 설명

Hello workflow는 인사말이 포함된 CSV 파일을 받아 별도의 파일에 작성하고, 각각을 대문자로 변환하고, 다시 모아서 인사말을 말하는 재미있는 캐릭터의 ASCII 그림이 포함된 단일 텍스트 파일을 출력합니다.

### Workflow 단계 (process)

네 단계는 별도의 모듈 파일에 저장된 Nextflow process(`sayHello`, `convertToUpper`, `collectGreetings`, `cowpy`)로 구현됩니다.

1. **`sayHello`:** 각 인사말을 자체 출력 파일에 작성합니다 (예: "Hello-output.txt")
2. **`convertToUpper`:** 각 인사말을 대문자로 변환합니다 (예: "HELLO")
3. **`collectGreetings`:** 모든 대문자 인사말을 단일 배치 파일로 수집합니다
4. **`cowpy`:** `cowpy` 도구를 사용하여 ASCII 아트를 생성합니다

### 다이어그램

<figure class="excalidraw">
--8<-- "docs/en/docs/hello_nextflow/img/hello_pipeline_complete.svg"
</figure>

### 결과

결과는 `results/`라는 디렉토리에 게시되며, pipeline의 최종 출력(기본 매개변수로 실행 시)은 대문자 인사말을 말하는 칠면조의 ASCII 아트가 포함된 일반 텍스트 파일입니다.

```txt title="results/cowpy-COLLECTED-test-batch-output.txt"
  _________
/ BONJOUR \
| HELLO   |
\ HOLà    /
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

pipeline이 포함된 과정에 따라 세부 사항에서 약간의 차이가 있을 수 있습니다.

---

<div markdown class="homepage_logos">

![Seqera](../assets/img/seqera_logo.png#only-light)

![Seqera](../assets/img/seqera_logo_dark.png#only-dark)

</div>
