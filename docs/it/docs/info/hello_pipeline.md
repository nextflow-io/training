---
title: La pipeline Hello
description: Riepilogo di cosa fa la pipeline Hello e come è strutturata.
hide:
  - toc
  - footer
---

# La pipeline Hello

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Traduzione assistita da IA - [scopri di più e suggerisci miglioramenti](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

La maggior parte dei nostri corsi di formazione utilizza una semplice pipeline agnostica rispetto al dominio per dimostrare i concetti e i meccanismi di Nextflow.
Il corso Hello Nextflow mostra come sviluppare questa pipeline passo dopo passo, spiegando ogni decisione di progettazione e implementazione.
Altri corsi di formazione utilizzano questa pipeline, o parti di essa, come punto di partenza.

Questa pagina riassume lo stato della pipeline al completamento del corso Hello Nextflow.

### Descrizione di riepilogo

Il workflow Hello prende un file CSV contenente saluti, li scrive in file separati, converte ciascuno in maiuscolo, li raccoglie nuovamente insieme e produce un singolo file di testo contenente un'immagine ASCII di un personaggio divertente che pronuncia i saluti.

### Passaggi del workflow (process)

I quattro passaggi sono implementati come process Nextflow (`sayHello`, `convertToUpper`, `collectGreetings` e `cowpy`) memorizzati in file modulo separati.

1. **`sayHello`:** Scrive ogni saluto nel proprio file di output (es. "Hello-output.txt")
2. **`convertToUpper`:** Converte ogni saluto in maiuscolo (es. "HELLO")
3. **`collectGreetings`:** Raccoglie tutti i saluti in maiuscolo in un singolo file batch
4. **`cowpy`:** Genera arte ASCII utilizzando lo strumento `cowpy`

### Diagramma

<figure class="excalidraw">
--8<-- "docs/en/docs/hello_nextflow/img/hello_pipeline_complete.svg"
</figure>

### Risultati

I risultati vengono pubblicati in una directory chiamata `results/`, e l'output finale della pipeline (quando eseguita con i parametri predefiniti) è un file di testo semplice contenente arte ASCII di un tacchino che pronuncia i saluti in maiuscolo.

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

Potrebbe incontrare alcune variazioni nei dettagli a seconda del corso in cui la pipeline è presente.

---

<div markdown class="homepage_logos">

![Seqera](../assets/img/seqera_logo.png#only-light)

![Seqera](../assets/img/seqera_logo_dark.png#only-dark)

</div>
