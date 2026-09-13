---
title: Nextflow Run
hide:
  - toc
page_type: index_page
index_type: course
additional_information:
  technical_requirements: true
  learning_objectives:
    - Lanzar y gestionar pipelines de Nextflow desde la línea de comandos
    - Comprender cómo los canales y operadores permiten workflows eficientes con múltiples entradas y múltiples pasos
    - Usar contenedores para gestionar dependencias de software y garantizar la reproducibilidad
    - Configurar la ejecución de pipelines y sus salidas
    - Generar reportes de ejecución, inspeccionar el historial de ejecuciones anteriores y limpiar directorios de trabajo antiguos
    - Ejecutar pipelines directamente desde repositorios remotos como GitHub
  audience_prerequisites:
    - "**Audiencia:** Este curso está diseñado para quienes son completamente nuevos en Nextflow y desean ejecutar pipelines existentes."
    - "**Habilidades:** Se asume cierta familiaridad con la línea de comandos, conceptos básicos de scripting y formatos de archivo comunes."
    - "**Dominio:** Los ejercicios son independientes del dominio, por lo que no se requiere conocimiento científico previo."
---

# Nextflow Run

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Traducción asistida por IA - [más información y sugerencias](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

**Nextflow Run es una introducción práctica a la ejecución de workflows de análisis de datos reproducibles y escalables.**

A través de una serie de ejercicios orientados a objetivos, aprenderá los fundamentos para lanzar y gestionar pipelines de Nextflow, comprenderá cómo los canales y operadores permiten el procesamiento paralelo de múltiples entradas, y usará contenedores para gestionar dependencias de software.

Se llevará las habilidades y la confianza para comenzar a ejecutar workflows con Nextflow.

<!-- additional_information -->

## Descripción general del curso

Este curso es práctico, con ejercicios orientados a objetivos estructurados para introducir información gradualmente.

Ejecutará varias versiones de un pipeline de Nextflow que procesa entradas de texto, comenzando con una versión simple de un solo paso y progresando a una versión de múltiples pasos que toma un archivo CSV de entradas, ejecuta algunos pasos de transformación y produce un único archivo de texto con arte ASCII generado por una herramienta en contenedor.

Este curso se enfoca en ejecutar pipelines (nombrado así por el comando principal `nextflow run`).
Si busca una introducción al desarrollo de pipelines de Nextflow, consulte [Hello Nextflow](../hello_nextflow/index.md).

!!! note "Nota"

    ¿Busca la versión anterior de este curso? Ha sido reemplazada por la versión en esta página, pero aún puede consultarse en la [versión 3.6.1](https://training.nextflow.io/3.6.1/nextflow_run/) del sitio de capacitación.

### Plan de lecciones

| Capítulo del curso                                                       | Resumen                                                                                                                      | Duración estimada |
| ------------------------------------------------------------------------ | ---------------------------------------------------------------------------------------------------------------------------- | ----------------- |
| [Parte 1: Ejecutar Nextflow](./01_run_nextflow.md)                       | Lanzar y gestionar pipelines de Nextflow, y comprender los mecanismos esenciales del workflow                                | 25 mins           |
| [Parte 2: Configurar el pipeline](./02_configure_pipeline.md)            | Configurar la ejecución del pipeline y sus salidas usando `nextflow.config`                                                  | 20 mins           |
| [Parte 3: Gestionar ejecuciones del workflow](./03_manage_executions.md) | Generar reportes de ejecución, inspeccionar el historial de ejecuciones anteriores y limpiar directorios de trabajo antiguos | 10 mins           |
| [Parte 4: Ejecutar pipelines remotos](./04_remote_repositories.md)       | Ejecutar un pipeline directamente desde GitHub y fijarlo a una revisión específica                                           | 10 mins           |

Al final de este curso, estará bien preparado para abordar los próximos pasos en su camino para ejecutar workflows reproducibles para sus necesidades de computación científica.

¿Listo para tomar el curso?

[Comenzar a aprender :material-arrow-right:](00_orientation.md){ .md-button .md-button--primary }
