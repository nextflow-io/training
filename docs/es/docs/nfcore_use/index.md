---
title: Usar nf-core
hide:
  - toc
page_type: index_page
index_type: course
additional_information:
  technical_requirements: true
  learning_objectives:
    - Encontrar, obtener y ejecutar pipelines comunitarios de nf-core
    - Configurar la ejecución de pipelines usando parámetros y archivos de configuración
    - Comprender cómo los pipelines de nf-core validan parámetros y datos de entrada
    - Ejecutar un pipeline a escala de producción (nf-core/rnaseq) y modificar sus asignaciones de recursos predeterminadas
  audience_prerequisites:
    - "**Audiencia:** Este curso está diseñado para personas que ya saben cómo ejecutar pipelines locales de Nextflow, son nuevas en nf-core y desean ejecutar pipelines comunitarios existentes."
    - "**Habilidades:** Se asume cierta familiaridad con la línea de comandos, conceptos básicos de scripting y formatos de archivo comunes."
    - "**Cursos:** Debe haber completado [Nextflow Run](../nextflow_run/index.md) o sentirse cómodo/a ejecutando un pipeline local con `nextflow run`."
    - "**Dominio:** Los ejercicios utilizan pipelines de bioinformática, pero no se requiere conocimiento científico previo del dominio."
---

# Usar nf-core

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Traducción asistida por IA - [más información y sugerencias](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

**Usar nf-core es una introducción práctica para encontrar, ejecutar y configurar pipelines comunitarios de nf-core.**

A través de ejemplos prácticos y ejercicios guiados, aprenderá a encontrar y obtener pipelines de nf-core, ejecutarlos usando sus perfiles de prueba integrados, y personalizar su ejecución mediante parámetros y archivos de configuración.

Al finalizar, contará con las habilidades y la confianza necesarias para comenzar a ejecutar pipelines de nf-core en sus propios análisis.

<!-- additional_information -->

## Descripción general del curso

Este curso es práctico, con ejercicios orientados a objetivos estructurados para introducir la información de forma gradual.

Comenzará con `nf-core/demo`, un pipeline mínimo mantenido por el proyecto nf-core con fines de capacitación, y luego aplicará lo aprendido a `nf-core/rnaseq`, un pipeline de producción ampliamente utilizado para el análisis de secuenciación de RNA en bulk.

Este curso se enfoca en la ejecución de pipelines.
Si busca una introducción al desarrollo de pipelines compatibles con nf-core, consulte [Build with nf-core](../nfcore_build/index.md).

### Plan de lecciones

| Capítulo del curso                                                             | Resumen                                                                                                                   | Duración estimada |
| ------------------------------------------------------------------------------ | ------------------------------------------------------------------------------------------------------------------------- | ----------------- |
| [Parte 1: Ejecutar un pipeline de demostración](./01_run_demo.md)              | Encontrar y obtener un pipeline de nf-core y ejecutarlo usando su perfil de prueba                                        | 20 min            |
| [Parte 2: Configurar la ejecución del pipeline](./02_configure_execution.md)   | Establecer parámetros, comprender la validación y personalizar la asignación de recursos y los argumentos de herramientas | 20 min            |
| [Parte 3: Ejecutar un pipeline de producción](./03_run_production_pipeline.md) | Obtener y ejecutar nf-core/rnaseq, y modificar sus asignaciones de recursos predeterminadas                               | 20 min            |

Al finalizar este curso, podrá aprovechar la gran cantidad de pipelines comunitarios que ofrece el proyecto nf-core.

¿Listo/a para comenzar el curso?

[Comenzar a aprender :material-arrow-right:](00_orientation.md){ .md-button .md-button--primary }
