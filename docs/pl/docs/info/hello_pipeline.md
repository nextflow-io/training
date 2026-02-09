---
title: Pipeline Hello
description: Podsumowanie działania i struktury pipeline'u Hello.
hide:
  - toc
  - footer
---

# Pipeline Hello

Większość naszych kursów szkoleniowych wykorzystuje prosty, niezależny od dziedziny pipeline do demonstracji koncepcji i mechanizmów Nextflow'a.
Kurs Hello Nextflow pokazuje, jak rozwijać ten pipeline krok po kroku, wyjaśniając każdą decyzję projektową i implementacyjną.
Inne szkolenia wykorzystują ten pipeline lub jego fragmenty jako punkt wyjścia.

Ta strona podsumowuje stan pipeline'u po ukończeniu kursu Hello Nextflow.

### Opis ogólny

Workflow Hello przyjmuje plik CSV zawierający powitania, zapisuje je do osobnych plików, konwertuje każde na wielkie litery, zbiera je z powrotem i wyświetla pojedynczy plik tekstowy zawierający obrazek ASCII zabawnej postaci wypowiadającej powitania.

### Kroki workflow'a (procesy)

Cztery kroki są zaimplementowane jako procesy Nextflow'a (`sayHello`, `convertToUpper`, `collectGreetings` i `cowpy`) przechowywane w osobnych plikach modułów.

1. **`sayHello`:** Zapisuje każde powitanie do własnego pliku wyjściowego (np. „Hello-output.txt")
2. **`convertToUpper`:** Konwertuje każde powitanie na wielkie litery (np. „HELLO")
3. **`collectGreetings`:** Zbiera wszystkie powitania pisane wielkimi literami do pojedynczego pliku wsadowego
4. **`cowpy`:** Generuje grafikę ASCII przy użyciu narzędzia `cowpy`

### Diagram

<figure class="excalidraw">
--8<-- "docs/en/docs/hello_nextflow/img/hello_pipeline_complete.svg"
</figure>

### Wyniki

Wyniki są publikowane do katalogu o nazwie `results/`, a końcowym rezultatem pipeline'u (przy uruchomieniu z domyślnymi parametrami) jest zwykły plik tekstowy zawierający grafikę ASCII indyka wypowiadającego powitania pisane wielkimi literami.

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

W zależności od kursu, w którym pipeline jest prezentowany, możesz napotkać pewne różnice w szczegółach.

---

<div markdown class="homepage_logos">

![Seqera](../assets/img/seqera_logo.png#only-light)

![Seqera](../assets/img/seqera_logo_dark.png#only-dark)

</div>
