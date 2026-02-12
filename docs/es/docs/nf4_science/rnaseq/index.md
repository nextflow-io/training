---
title: RNAseq
hide:
  - toc
page_type: index_page
index_type: course
additional_information:
  technical_requirements: true
  learning_objectives:
    - Escribir un workflow lineal para aplicar métodos básicos de procesamiento de RNAseq y QC
    - Manejar apropiadamente archivos específicos del dominio como FASTQ y recursos de genoma de referencia
    - Manejar datos de secuenciación de extremo simple y extremo pareado
    - Aprovechar el paradigma de flujo de datos de Nextflow para paralelizar el procesamiento de RNAseq por muestra
    - Agregar reportes de QC a través de múltiples pasos y muestras usando operadores de canal relevantes
  audience_prerequisites:
    - "**Audiencia:** Este curso está diseñado para investigadores en transcriptómica y campos relacionados que deseen desarrollar o personalizar pipelines de análisis de datos."
    - "**Habilidades:** Se asume cierta familiaridad con la línea de comandos, conceptos básicos de scripting y formatos de archivo comunes de RNAseq."
    - "**Requisitos previos:** Conceptos fundamentales de Nextflow y herramientas cubiertas en [Hello Nextflow](../../hello_nextflow/)."
---

# Nextflow para RNAseq

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Traducción asistida por IA - [más información y sugerencias](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

**Un curso práctico que aplica Nextflow a un caso de uso real de transcriptómica: procesamiento de RNAseq masivo con Trim Galore, HISAT2 y FastQC.**

Este curso se basa en el entrenamiento para principiantes [Hello Nextflow](../../hello_nextflow/) y demuestra cómo usar Nextflow en el contexto específico del análisis de RNAseq masivo.
Implementará un pipeline de procesamiento que recorta secuencias adaptadoras, alinea las lecturas a un genoma de referencia y realiza control de calidad (QC) en varias etapas.

<!-- additional_information -->

## Descripción general del curso

Este curso es práctico, con ejercicios orientados a objetivos estructurados para introducir información gradualmente.

Comenzará ejecutando las herramientas de procesamiento manualmente en la terminal para comprender la metodología, luego construirá progresivamente un pipeline de Nextflow que automatiza y escala el análisis.

### Plan de lecciones

Hemos dividido esto en tres partes que se enfocan en aspectos específicos de aplicar Nextflow a un caso de uso de RNAseq.

| Capítulo del curso                                                               | Resumen                                                                                                                  | Duración estimada |
| -------------------------------------------------------------------------------- | ------------------------------------------------------------------------------------------------------------------------ | ----------------- |
| [Parte 1: Descripción del método](./01_method.md)                                | Comprender la metodología de procesamiento de RNAseq y ejecutar las herramientas manualmente                             | 30 mins           |
| [Parte 2: Implementación de muestra única](./02_single-sample.md)                | Construir un pipeline que recorta, alinea y realiza QC de una sola muestra, luego escala para manejar múltiples muestras | 60 mins           |
| [Parte 3: Implementación multi-muestra de extremo pareado](./03_multi-sample.md) | Extender el pipeline para manejar datos de extremo pareado y agregar reportes de QC a través de muestras                 | 45 mins           |

Al finalizar este curso, podrá aplicar conceptos fundamentales de Nextflow y herramientas a un caso de uso típico de RNAseq.

¿Listo para tomar el curso?

[Comenzar :material-arrow-right:](00_orientation.md){ .md-button .md-button--primary }

<!-- Clearfix for float -->
<div style="content: ''; clear: both; display: table;"></div>
