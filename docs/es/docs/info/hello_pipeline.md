---
title: El pipeline Hello
description: Resumen de lo que hace el pipeline Hello y cómo está estructurado.
hide:
  - toc
  - footer
---

# El pipeline Hello

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Traducción asistida por IA - [más información y sugerencias](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

La mayoría de nuestros cursos de entrenamiento utilizan un pipeline simple e independiente del dominio para demostrar conceptos y mecanismos de Nextflow.
El curso Hello Nextflow muestra cómo desarrollar este pipeline paso a paso, explicando cada decisión de diseño e implementación.
Otros entrenamientos utilizan este pipeline, o partes de él, como punto de partida.

Esta página resume el estado del pipeline tal como está al completar el curso Hello Nextflow.

### Descripción resumida

El workflow Hello toma un archivo CSV que contiene saludos, los escribe en archivos separados, convierte cada uno a mayúsculas, los recopila nuevamente y genera un único archivo de texto que contiene una imagen ASCII de un personaje divertido diciendo los saludos.

### Pasos del workflow (procesos)

Los cuatro pasos están implementados como processes de Nextflow (`sayHello`, `convertToUpper`, `collectGreetings` y `cowpy`) almacenados en archivos de módulo separados.

1. **`sayHello`:** Escribe cada saludo en su propio archivo de salida (por ejemplo, "Hello-output.txt")
2. **`convertToUpper`:** Convierte cada saludo a mayúsculas (por ejemplo, "HELLO")
3. **`collectGreetings`:** Recopila todos los saludos en mayúsculas en un único archivo de lote
4. **`cowpy`:** Genera arte ASCII usando la herramienta `cowpy`

### Diagrama

<figure class="excalidraw">
--8<-- "docs/en/docs/hello_nextflow/img/hello_pipeline_complete.svg"
</figure>

### Resultados

Los resultados se publican en un directorio llamado `results/`, y la salida final del pipeline (cuando se ejecuta con parámetros predeterminados) es un archivo de texto plano que contiene arte ASCII de un pavo diciendo los saludos en mayúsculas.

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

Puede encontrar algunas variaciones en los detalles específicos dependiendo del curso en el que se presenta el pipeline.

---

<div markdown class="homepage_logos">

![Seqera](../assets/img/seqera_logo.png#only-light)

![Seqera](../assets/img/seqera_logo_dark.png#only-dark)

</div>
