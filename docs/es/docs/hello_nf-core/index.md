---
title: Hello nf-core
hide:
  - toc
page_type: index_page
index_type: course
additional_information:
  technical_requirements: true
  learning_objectives:
    - Recuperar, ejecutar y gestionar la ejecución de pipelines de nf-core
    - Describir la estructura del código y la organización del proyecto de pipelines de nf-core
    - Crear un pipeline básico compatible con nf-core desde una plantilla
    - Actualizar un workflow de Nextflow básico para cumplir con los estándares de nf-core
    - Añadir módulos de nf-core a un pipeline compatible con nf-core
    - Contribuir con sus propios módulos a nf-core
    - Validar entradas y parámetros utilizando las herramientas de nf-core
  audience_prerequisites:
    - "**Audiencia:** Este curso está diseñado para estudiantes que ya están familiarizados con Nextflow básico y desean aprender a usar recursos y mejores prácticas de nf-core."
    - "**Habilidades:** Se asume familiaridad con la línea de comandos, conceptos básicos de scripting y formatos de archivo comunes."
    - "**Cursos:** Debe haber completado el curso [Hello Nextflow](../hello_nextflow/index.md) o equivalente."
    - "**Dominio:** Los ejercicios son todos agnósticos al dominio, por lo que no se requiere conocimiento científico previo."
---

# Hello nf-core

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Traducción asistida por IA - [más información y sugerencias](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

**Hello nf-core es una introducción práctica al uso de recursos y mejores prácticas de nf-core.**

![nf-core logo](./img/nf-core-logo.png)

Trabajando a través de ejemplos prácticos y ejercicios guiados, aprenderá a usar y desarrollar módulos y pipelines compatibles con nf-core, y a utilizar las herramientas de nf-core de manera efectiva.

Obtendrá las habilidades y la confianza para comenzar a desarrollar pipelines de acuerdo con las mejores prácticas de nf-core.

<!-- additional_information -->

## Descripción general del curso

Este curso está diseñado para ser práctico, con ejercicios orientados a objetivos estructurados para introducir información gradualmente.

Se le presentará [**nf-core**](https://nf-co.re/), un esfuerzo comunitario para desarrollar y mantener un conjunto curado de pipelines científicos construidos usando Nextflow, así como herramientas y directrices relevantes que promueven el desarrollo abierto, las pruebas y la revisión por pares ([Nat Biotechnol 38, 276–278 (2020)](https://www.nature.com/articles/s41587-020-0439-x), [Genome Biol 26, 228 (2025)](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-025-03673-9)).

Los pipelines desarrollados por la comunidad nf-core están diseñados para ser modulares, escalables y portables, permitiendo a los investigadores adaptarlos y ejecutarlos fácilmente usando sus propios datos y recursos de cómputo.
Las directrices de mejores prácticas aplicadas por el proyecto aseguran además que los pipelines sean robustos, estén bien documentados y validados contra conjuntos de datos del mundo real.
Esto ayuda a aumentar la confiabilidad y reproducibilidad de los análisis científicos y, en última instancia, permite a los investigadores acelerar sus descubrimientos científicos.

No cubriremos todo lo que hay que saber sobre pipelines de nf-core en este curso, porque nf-core abarca muchas características y convenciones desarrolladas por la comunidad a lo largo de años.
En su lugar, nos centraremos en los conceptos esenciales que le ayudarán a comenzar y entender cómo funciona nf-core.

### Plan de lecciones

Hemos dividido esto en cinco partes que se enfocarán cada una en aspectos específicos del uso de recursos de nf-core.

| Capítulo del curso                                                | Resumen                                                                                                                                                                           | Duración estimada |
| ----------------------------------------------------------------- | --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- | ----------------- |
| [Parte 1: Ejecutar un pipeline de demostración](./01_run_demo.md) | Ejecutar un pipeline existente de nf-core y examinar su estructura de código para tener una idea de lo que hace diferentes a estos pipelines de los workflows básicos de Nextflow | 30 mins           |
| [Parte 2: Reescribir Hello para nf-core](./02_rewrite_hello.md)   | Adaptar un workflow existente a la estructura de plantilla de nf-core, comenzando desde el workflow simple producido en el curso [Hello Nextflow](../hello_nextflow/index.md)     | 60 mins           |
| [Parte 3: Usar un módulo de nf-core](./03_use_module.md)          | Explorar la biblioteca de módulos de la comunidad y aprender a integrar módulos preconstruidos y probados que envuelven herramientas bioinformáticas comunes                      | 30 mins           |
| [Parte 4: Crear un módulo de nf-core](./04_make_module.md)        | Crear su propio módulo al estilo nf-core usando la estructura específica, convenciones de nomenclatura y requisitos de metadatos establecidos por nf-core                         | 30 mins           |
| [Parte 5: Añadir validación de entrada](./05_input_validation.md) | Implementar validación de entrada tanto para parámetros de línea de comandos como para archivos de datos de entrada usando nf-schema                                              | 30 mins           |

Al final de este curso, podrá aprovechar la enorme riqueza de recursos ofrecidos por el proyecto nf-core.

¿Listo para tomar el curso?

[Comenzar a aprender :material-arrow-right:](00_orientation.md){ .md-button .md-button--primary }
