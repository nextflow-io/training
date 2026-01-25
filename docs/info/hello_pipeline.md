---
title: The Hello pipeline
description: Recap of what the Hello pipeline does and how it is structured.
  - toc
  - footer
---

# The Hello pipeline

Most of our training courses use a simple domain-agnostic pipeline to demonstrate Nextflow concepts and mechanisms.
The Hello Nextflow course shows how to develop this pipeline in a step-by step way that explains every design and implementation decision.
Other trainings use this pipeline, or parts of it, as a starting point.

This page summarizes the state of the pipeline as it stands at the completion of the Hello Nextflow course.

### 1. Summary description

The Hello workflow takes a CSV file containing greetings, writes them to separate files, converts each to uppercase, collects them back together again and outputs a single text file containing an ASCII picture of a fun character saying the greetings.

### 2. Workflow steps (processes)

The four steps are implemented as Nextflow processes (`sayHello`, `convertToUpper`, `collectGreetings`, and `cowpy`) stored in separate module files.

1. **`sayHello`:** Writes each greeting to its own output file (e.g., "Hello-output.txt")
2. **`convertToUpper`:** Converts each greeting to uppercase (e.g., "HELLO")
3. **`collectGreetings`:** Collects all uppercase greetings into a single batch file
4. **`cowpy`:** Generates ASCII art using the `cowpy` tool

### 3. Diagram

<figure class="excalidraw">
--8<-- "docs/hello_nextflow/img/hello_pipeline_complete.svg"
</figure>

### 4. Results

The results are published to a directory called `results/`, and the final output of the pipeline (when run with default parameters) is a plain text file containing ASCII art of a turkey saying the uppercased greetings.

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

You may encounter some variations in the specifics depending on the course the pipeline is featured in.

---

<div markdown class="homepage_logos">

![Seqera](assets/img/seqera_logo.png#only-light)

![Seqera](assets/img/seqera_logo_dark.png#only-dark)

</div>
