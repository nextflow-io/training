---
title: El pipeline Hello
description: Resum del que fa el pipeline Hello i com està estructurat.
hide:
  - toc
  - footer
---

# El pipeline Hello

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Traducció assistida per IA - [més informació i suggeriments](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

La majoria dels nostres cursos de formació utilitzen un pipeline senzill independent del domini per demostrar conceptes i mecanismes de Nextflow.
El curs Hello Nextflow mostra com desenvolupar aquest pipeline pas a pas, explicant cada decisió de disseny i implementació.
Altres formacions utilitzen aquest pipeline, o parts d'ell, com a punt de partida.

Aquesta pàgina resumeix l'estat del pipeline tal com queda en completar el curs Hello Nextflow.

### Descripció resumida

El workflow Hello pren un fitxer CSV que conté salutacions, les escriu en fitxers separats, converteix cadascuna a majúscules, les recull de nou juntes i genera un únic fitxer de text que conté una imatge ASCII d'un personatge divertit dient les salutacions.

### Passos del workflow (processos)

Els quatre passos s'implementen com a processos de Nextflow (`sayHello`, `convertToUpper`, `collectGreetings` i `cowpy`) emmagatzemats en fitxers de mòdul separats.

1. **`sayHello`:** Escriu cada salutació al seu propi fitxer de sortida (p. ex., "Hello-output.txt")
2. **`convertToUpper`:** Converteix cada salutació a majúscules (p. ex., "HELLO")
3. **`collectGreetings`:** Recull totes les salutacions en majúscules en un únic fitxer per lots
4. **`cowpy`:** Genera art ASCII utilitzant l'eina `cowpy`

### Diagrama

<figure class="excalidraw">
--8<-- "docs/en/docs/hello_nextflow/img/hello_pipeline_complete.svg"
</figure>

### Resultats

Els resultats es publiquen en un directori anomenat `results/`, i la sortida final del pipeline (quan s'executa amb paràmetres per defecte) és un fitxer de text pla que conté art ASCII d'un gall dindi dient les salutacions en majúscules.

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

Podeu trobar algunes variacions en els detalls específics depenent del curs en què aparegui el pipeline.

---

<div markdown class="homepage_logos">

![Seqera](../assets/img/seqera_logo.png#only-light)

![Seqera](../assets/img/seqera_logo_dark.png#only-dark)

</div>
