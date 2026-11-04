---
title: Configuración de Ejecución
hide:
  - toc
page_type: index_page
index_type: course
additional_information:
  technical_requirements: true
  learning_objectives:
    - Cambiar la tecnología de empaquetado de software entre Docker y Conda
    - Seleccionar una plataforma de ejecución y comprender cómo Nextflow adapta la ejecución de tareas a ella
    - Controlar la asignación de recursos de cómputo y reintentar automáticamente las tareas que fallen
    - Definir y combinar perfiles para alternar entre configuraciones predefinidas
  audience_prerequisites:
    - "**Audiencia:** Este curso está diseñado para personas que ya saben cómo ejecutar pipelines de Nextflow localmente y desean configurar la ejecución con mayor profundidad."
    - "**Habilidades:** Se asume cierta familiaridad con la línea de comandos."
    - "**Cursos:** Debe haber completado [Nextflow Run](../nextflow_run/index.md) o estar familiarizado/a con la ejecución de un pipeline local con `nextflow run`."
---

# Configuración de Ejecución

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Traducción asistida por IA - [más información y sugerencias](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

**Configure Execution es una introducción práctica a la adaptación de la ejecución de pipelines de Nextflow a diferentes entornos de cómputo.**

A través de ejercicios orientados a objetivos, aprenderá a cambiar la tecnología de empaquetado de software, seleccionar una plataforma de ejecución, controlar la asignación de recursos de cómputo y los reintentos, y agrupar la configuración en perfiles intercambiables.

Adquirirá las habilidades y la confianza para configurar la ejecución de pipelines de Nextflow como un profesional.

<!-- additional_information -->

## Descripción general del curso

Este curso es práctico y se basa en las habilidades cubiertas en [Nextflow Run](../nextflow_run/index.md).

Tomará el mismo pipeline de múltiples pasos de ese curso y adaptará progresivamente su configuración a diferentes entornos de cómputo, para luego agrupar todo en perfiles entre los que podrá alternar en tiempo de ejecución.

### Plan de lecciones

| Capítulo del curso                                                                | Resumen                                                                                    | Duración estimada |
| --------------------------------------------------------------------------------- | ------------------------------------------------------------------------------------------ | ----------------- |
| [Parte 1: Adaptación al entorno de cómputo](./01_packaging_and_execution.md)      | Cambiar la tecnología de empaquetado de software y seleccionar una plataforma de ejecución | 20 min            |
| [Parte 2: Gestión de recursos de cómputo y fallos](./02_resources_and_retries.md) | Controlar la asignación de recursos y reintentar automáticamente las tareas que fallen     | 15 min            |
| [Parte 3: Uso de perfiles para alternar configuraciones](./03_profiles.md)        | Definir y combinar perfiles, e inspeccionar la configuración completamente resuelta        | 15 min            |

Al finalizar este curso, estará en condiciones de configurar pipelines de Nextflow para una variedad de entornos de cómputo y alternar entre ellos con el mínimo esfuerzo.

¿Listo/a para comenzar el curso?

[Comenzar a aprender :material-arrow-right:](00_orientation.md){ .md-button .md-button--primary }
