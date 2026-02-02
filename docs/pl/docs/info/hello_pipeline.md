---
title: Pipeline Hello
description: Podsumowanie działania pipeline'u Hello i jego struktury.
hide:
  - toc
  - footer
---

# Pipeline Hello

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Tłumaczenie wspomagane przez AI - [dowiedz się więcej i zasugeruj ulepszenia](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

Większość naszych kursów szkoleniowych używa prostego, domenowo-agnostycznego pipeline'u do demonstracji koncepcji i mechanizmów Nextflow.
Kurs Hello Nextflow pokazuje, jak rozwijać ten pipeline krok po kroku, wyjaśniając każdą decyzję projektową i implementacyjną.
Inne szkolenia używają tego pipeline'u lub jego części jako punktu wyjścia.

Ta strona podsumowuje stan pipeline'u po ukończeniu kursu Hello Nextflow.

### Krótki opis

Workflow Hello przyjmuje plik CSV zawierający teksty powitalne. Zapisuje je do oddzielnych plików i konwertuje na wielkie litery. Następnie zbiera je z powrotem razem i generuje pojedynczy plik tekstowy z obrazkiem ASCII zabawnej postaci wypowiadającej te teksty.

### Kroki workflow (procesy)

Cztery kroki są zaimplementowane jako procesy Nextflow (`sayHello`, `convertToUpper`, `collectGreetings` i `cowpy`) przechowywane w oddzielnych plikach modułów.

1. **`sayHello`:** Zapisuje każdy tekst powitalny do własnego pliku wyjściowego (np. "Hello-output.txt")
2. **`convertToUpper`:** Konwertuje każdy wpis na wielkie litery (np. "HELLO")
3. **`collectGreetings`:** Zbiera wszystkie pozdrowienia wielkimi literami do pojedynczego pliku wsadowego
4. **`cowpy`:** Generuje grafikę ASCII za pomocą narzędzia `cowpy`

### Diagram

<figure class="excalidraw">
--8<-- "docs/en/docs/hello_nextflow/img/hello_pipeline_complete.svg"
</figure>

### Wyniki

Wyniki są publikowane do katalogu o nazwie `results/`, a końcowe wyjście pipeline'u (przy uruchomieniu z domyślnymi parametrami) to plik tekstowy zawierający grafikę ASCII indyka wypowiadającego pozdrowienia wielkimi literami.

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

Możesz napotkać pewne różnice w szczegółach w zależności od kursu, w którym pipeline jest prezentowany.

---

<div markdown class="homepage_logos">

![Seqera](../assets/img/seqera_logo.png#only-light)

![Seqera](../assets/img/seqera_logo_dark.png#only-dark)

</div>
