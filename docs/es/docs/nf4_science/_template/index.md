---
title: Nextflow for {DOMAIN}
hide:
  - toc
page_type: index_page
index_type: course
additional_information:
  technical_requirements: true
  learning_objectives:
    - Write a linear workflow to apply {METHOD} to a single sample
    - Handle accessory files such as {ACCESSORY_FILES} appropriately
    - Leverage Nextflow's dataflow paradigm to parallelize per-sample processing
    - Implement multi-sample aggregation using relevant channel operators
  audience_prerequisites:
    - "**Audience:** This course is designed for researchers in {DOMAIN} and related fields who want to develop or customize data analysis pipelines."
    - "**Skills:** Some familiarity with the command line, basic scripting concepts, and common {DOMAIN} file formats is assumed."
    - "**Prerequisites:** Foundational Nextflow concepts and tooling covered in [Hello Nextflow](../../hello_nextflow/)."
---

# Nextflow para {DOMAIN}

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Traducción asistida por IA - [más información y sugerencias](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

**Un curso práctico que aplica Nextflow a un caso de uso real de {DOMAIN}: {METHOD_SHORT_DESCRIPTION}.**

Este curso se basa en la capacitación para principiantes [Hello Nextflow](../../hello_nextflow/) y demuestra cómo usar Nextflow en el contexto específico del dominio de {DOMAIN}.
Implementarás un pipeline de {METHOD} con [{TOOL_A}]({TOOL_A_URL}) y [{TOOL_B}]({TOOL_B_URL}).

<!-- additional_information -->

## Descripción general del curso

Este curso es práctico, con ejercicios orientados a objetivos estructurados para introducir información gradualmente.

Comenzarás ejecutando las herramientas de análisis manualmente en la terminal para comprender la metodología, luego construirás progresivamente un pipeline de Nextflow que automatiza y escala el análisis.

### Plan de lecciones

Hemos dividido esto en tres partes que se enfocan en aspectos específicos de la aplicación de Nextflow a un caso de uso de {DOMAIN}.

| Capítulo del curso                                                | Resumen                                                                                                                  | Duración estimada |
| ----------------------------------------------------------------- | ------------------------------------------------------------------------------------------------------------------------ | ----------------- |
| [Parte 1: Descripción general del método](./01_method.md)         | Comprender la metodología de {METHOD} y ejecutar las herramientas manualmente                                           | 30 mins           |
| [Parte 2: Procesamiento de muestra única](./02_single_sample.md)  | Construir un pipeline que {PART2_SUMMARY}, luego escalar a múltiples muestras                                           | 60 mins           |
| [Parte 3: Agregación de múltiples muestras](./03_multi_sample.md) | Agregar {AGGREGATION_SUMMARY} de múltiples muestras usando operadores de canal para agregar salidas por muestra         | 45 mins           |

Al final de este curso, podrás aplicar conceptos y herramientas fundamentales de Nextflow a un caso de uso típico de {DOMAIN}.

¿Listo para tomar el curso?

[Comenzar :material-arrow-right:](00_orientation.md){ .md-button .md-button--primary }

<!-- Clearfix for float -->
<div style="content: ''; clear: both; display: table;"></div>
