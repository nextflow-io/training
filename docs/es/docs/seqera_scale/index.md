---
title: Escala con Seqera
hide:
  - toc
page_type: index_page
index_type: course
additional_information:
  technical_requirements: true
  learning_objectives:
    - Registrarse en Seqera Platform y explorar el Community Showcase
    - Agregar un pipeline a un workspace y ejecutarlo desde la interfaz web
    - Autenticarse y ejecutar pipelines desde la línea de comandos con el CLI `tw`
    - Registrar un pipeline alojado en GitHub y ejecutarlo de ambas formas
  audience_prerequisites:
    - "**Audiencia:** Este curso está diseñado para personas que desean ejecutar pipelines de Nextflow a escala usando Seqera Platform."
    - "**Habilidades:** Se asume familiaridad con la ejecución de pipelines de nf-core desde la línea de comandos."
    - "**Cursos:** Es necesario haber completado [Nextflow Run](../nextflow_run/index.md) y [Use nf-core](../nfcore_use/index.md), o tener experiencia ejecutando pipelines locales y el pipeline `nf-core/rnaseq`."
---

# Escala con Seqera

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Traducción asistida por IA - [más información y sugerencias](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

**Scale with Seqera es una introducción práctica para ejecutar y monitorear pipelines de Nextflow con Seqera Platform.**

A través de ejemplos prácticos, configurará el acceso a Seqera Platform, ejecutará un pipeline a escala de producción tanto desde la interfaz web como desde la línea de comandos, y agregará un nuevo pipeline a su workspace.

Al finalizar, contará con las habilidades y la confianza necesarias para ejecutar y monitorear sus propios pipelines en Seqera Platform.

<!-- additional_information -->

## Descripción general del curso

Este curso es práctico y se basa en los pipelines que ya ejecutó en [Use nf-core](../nfcore_use/index.md).

Comenzará registrándose en Seqera Platform y ejecutando `nf-core/rnaseq`, un pipeline a escala de producción, desde la interfaz web.
Luego cambiará a la herramienta de línea de comandos `tw` para hacer lo mismo desde una terminal y, finalmente, registrará un nuevo pipeline, `nf-core/demo`, y lo ejecutará de ambas formas.

### Plan de lecciones

| Capítulo del curso                                                                          | Resumen                                                                                                          | Duración estimada |
| ------------------------------------------------------------------------------------------- | ---------------------------------------------------------------------------------------------------------------- | ----------------- |
| [Parte 1: Ejecutar pipelines desde la interfaz web](./01_run_with_seqera.md)                | Configurar el acceso a Seqera Platform y ejecutar un pipeline a escala de producción desde la interfaz web       | 20 min            |
| [Parte 2: Ejecutar pipelines desde la línea de comandos](./02_launch_from_cli.md)           | Autenticar el CLI `tw`, ejecutar un pipeline guardado y registrar un nuevo pipeline desde la línea de comandos   | 25 min            |

Al finalizar este curso, se sentirá cómodo/a ejecutando y monitoreando pipelines de Nextflow en Seqera Platform, ya sea desde la interfaz web o desde la línea de comandos.

¿Listo/a para comenzar el curso?

[Comenzar a aprender :material-arrow-right:](00_orientation.md){ .md-button .md-button--primary }
