---
title: O pipeline Hello
description: Recapitulação do que o pipeline Hello faz e como ele é estruturado.
hide:
  - toc
  - footer
---

# O pipeline Hello

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Tradução assistida por IA - [saiba mais e sugira melhorias](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

A maioria dos nossos cursos de treinamento usa um pipeline simples e agnóstico de domínio para demonstrar conceitos e mecanismos do Nextflow.
O curso Hello Nextflow mostra como desenvolver este pipeline passo a passo, explicando cada decisão de design e implementação.
Outros treinamentos usam este pipeline, ou partes dele, como ponto de partida.

Esta página resume o estado do pipeline como está ao final do curso Hello Nextflow.

### Descrição resumida

O fluxo de trabalho Hello recebe um arquivo CSV contendo saudações, escreve-as em arquivos separados, converte cada uma para maiúsculas, coleta-as de volta juntas e gera um único arquivo de texto contendo uma imagem ASCII de um personagem divertido dizendo as saudações.

### Etapas do fluxo de trabalho (processos)

As quatro etapas são implementadas como processos Nextflow (`sayHello`, `convertToUpper`, `collectGreetings` e `cowpy`) armazenados em arquivos de módulo separados.

1. **`sayHello`:** Escreve cada saudação em seu próprio arquivo de saída (ex: "Hello-output.txt")
2. **`convertToUpper`:** Converte cada saudação para maiúsculas (ex: "HELLO")
3. **`collectGreetings`:** Coleta todas as saudações em maiúsculas em um único arquivo de lote
4. **`cowpy`:** Gera arte ASCII usando a ferramenta `cowpy`

### Diagrama

<figure class="excalidraw">
--8<-- "docs/en/docs/hello_nextflow/img/hello_pipeline_complete.svg"
</figure>

### Resultados

Os resultados são publicados em um diretório chamado `results/`, e a saída final do pipeline (quando executado com parâmetros padrão) é um arquivo de texto contendo arte ASCII de um peru dizendo as saudações em maiúsculas.

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

Você pode encontrar algumas variações nos detalhes dependendo do curso em que o pipeline é apresentado.

---

<div markdown class="homepage_logos">

![Seqera](../assets/img/seqera_logo.png#only-light)

![Seqera](../assets/img/seqera_logo_dark.png#only-dark)

</div>
