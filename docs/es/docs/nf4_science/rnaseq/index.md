# Nextflow para RNAseq

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Traducción asistida por IA - [más información y sugerencias](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

Este curso de entrenamiento está dirigido a investigadores en transcriptómica y campos relacionados que estén interesados en desarrollar o personalizar pipelines de análisis de datos.
Se basa en el entrenamiento para principiantes [Hello Nextflow](../../hello_nextflow/) y demuestra cómo usar Nextflow en el contexto específico del análisis de RNAseq masivo.

Específicamente, este curso demuestra cómo implementar un pipeline simple de procesamiento de RNAseq masivo para recortar secuencias adaptadoras, alinear las lecturas a un genoma de referencia y realizar control de calidad (QC) en varias etapas.

¡Comencemos! Haga clic en el botón "Open in GitHub Codespaces" a continuación para iniciar el entorno de entrenamiento (preferiblemente en una pestaña separada), luego continúe leyendo mientras se carga.

[![Open in GitHub Codespaces](https://github.com/codespaces/badge.svg)](https://codespaces.new/nextflow-io/training?quickstart=1&ref=master)

## Objetivos de aprendizaje

Al trabajar en este curso, aprenderá a aplicar conceptos fundamentales de Nextflow y herramientas a un caso de uso típico de RNAseq.

Al finalizar este taller, usted será capaz de:

- Escribir un workflow lineal para aplicar métodos básicos de procesamiento de RNAseq y QC
- Manejar apropiadamente archivos específicos del dominio como FASTQ y recursos de genoma de referencia
- Manejar datos de secuenciación de extremo simple y extremo pareado
- Aprovechar el paradigma de flujo de datos de Nextflow para paralelizar el procesamiento de RNAseq por muestra
- Agregar reportes de QC a través de múltiples pasos y muestras usando operadores de channel relevantes

<!-- TODO
- Configure pipeline execution and manage and optimize resource allocations
- Implement per-step and end-to-end pipeline tests that handle RNAseq-specific idiosyncrasies appropriately
-->
<!-- TODO for future expansion: add metadata/samplesheet handling -->

## Requisitos previos

El curso asume cierta familiaridad mínima con lo siguiente:

- Herramientas y formatos de archivo comúnmente usados en este dominio científico
- Experiencia con la línea de comandos
- Conceptos fundamentales de Nextflow y herramientas cubiertas en el entrenamiento para principiantes [Hello Nextflow](../../hello_nextflow/).

Para requisitos técnicos y configuración del entorno, consulte el mini-curso [Configuración del Entorno](../../envsetup/).
