---
title: Nextflow para Genómica
hide:
  - toc
page_type: index_page
index_type: course
additional_information:
  technical_requirements: true
  learning_objectives:
    - Escribir un workflow lineal para aplicar llamado de variantes a una sola muestra
    - Manejar archivos accesorios como archivos de índice y recursos del genoma de referencia de manera apropiada
    - Aprovechar el paradigma de flujo de datos de Nextflow para paralelizar el llamado de variantes por muestra
    - Implementar llamado conjunto de múltiples muestras usando operadores de canal relevantes
  audience_prerequisites:
    - "**Audiencia:** Este curso está diseñado para investigadores en genómica y campos relacionados que desean desarrollar o personalizar pipelines de análisis de datos."
    - "**Habilidades:** Se asume cierta familiaridad con la línea de comandos, conceptos básicos de scripting y formatos comunes de archivos genómicos."
    - "**Requisitos previos:** Conceptos fundamentales de Nextflow y herramientas cubiertas en [Hello Nextflow](../../hello_nextflow/)."
---

# Nextflow para Genómica

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Traducción asistida por IA - [más información y sugerencias](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

**Un curso práctico que aplica Nextflow a un caso de uso real en genómica: llamado de variantes con GATK.**

Este curso se basa en la capacitación para principiantes [Hello Nextflow](../../hello_nextflow/) y demuestra cómo usar Nextflow en el contexto específico del dominio de la genómica.
Implementarás un pipeline de llamado de variantes con [GATK](https://gatk.broadinstitute.org/) (Genome Analysis Toolkit), un paquete de software ampliamente utilizado para analizar datos de secuenciación de alto rendimiento.

<!-- additional_information -->

## Descripción general del curso

Este curso es práctico, con ejercicios orientados a objetivos estructurados para introducir información gradualmente.

Comenzarás ejecutando las herramientas de llamado de variantes manualmente en la terminal para comprender la metodología, luego construirás progresivamente un pipeline de Nextflow que automatiza y escala el análisis.

### Plan de lecciones

Hemos dividido esto en tres partes que se enfocan en aspectos específicos de aplicar Nextflow a un caso de uso en genómica.

| Capítulo del curso                                                              | Resumen                                                                                                       | Duración estimada |
| ------------------------------------------------------------------------------- | ------------------------------------------------------------------------------------------------------------- | ----------------- |
| [Parte 1: Descripción del método](./01_method.md)                               | Comprender la metodología de llamado de variantes y ejecutar las herramientas manualmente                     | 30 mins           |
| [Parte 2: Llamado de variantes por muestra](./02_per_sample_variant_calling.md) | Construir un pipeline que indexa archivos BAM y llama variantes, luego escalar a múltiples muestras           | 60 mins           |
| [Parte 3: Llamado conjunto en una cohorte](./03_joint_calling.md)               | Agregar genotipado conjunto de múltiples muestras usando operadores de canal para agregar salidas por muestra | 45 mins           |

Al final de este curso, podrás aplicar conceptos fundamentales de Nextflow y herramientas a un caso de uso típico en genómica.

¿Listo/a para tomar el curso?

[Comenzar :material-arrow-right:](00_orientation.md){ .md-button .md-button--primary }

<!-- Clearfix for float -->
<div style="content: ''; clear: both; display: table;"></div>
