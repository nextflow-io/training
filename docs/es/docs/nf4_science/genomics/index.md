# Nextflow para Genómica

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Traducción asistida por IA - [más información y sugerencias](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

Este curso de entrenamiento está dirigido a investigadores en genómica y campos relacionados que estén interesados en desarrollar o personalizar pipelines de análisis de datos.
Se basa en el entrenamiento para principiantes [Hello Nextflow](../../hello_nextflow/) y demuestra cómo usar Nextflow en el contexto específico del dominio de la genómica.

Específicamente, este curso demuestra cómo implementar un pipeline simple de llamado de variantes con [GATK](https://gatk.broadinstitute.org/) (Genome Analysis Toolkit), un paquete de software ampliamente utilizado para analizar datos de secuenciación de alto rendimiento.

¡Comencemos! Haga clic en el botón "Open in GitHub Codespaces" a continuación para iniciar el entorno de entrenamiento (preferiblemente en una pestaña separada), luego continúe leyendo mientras se carga.

[![Open in GitHub Codespaces](https://github.com/codespaces/badge.svg)](https://codespaces.new/nextflow-io/training?quickstart=1&ref=master)

## Objetivos de aprendizaje

Al completar este curso, aprenderá cómo aplicar conceptos y herramientas fundamentales de Nextflow a un caso de uso típico de genómica.

Al final de este taller será capaz de:

- Escribir un workflow lineal para aplicar llamado de variantes a una única muestra
- Manejar archivos accesorios como archivos de índice y recursos del genoma de referencia de manera apropiada
- Aprovechar el paradigma de flujo de datos de Nextflow para paralelizar el llamado de variantes por muestra
- Implementar llamado de variantes de múltiples muestras utilizando operadores de canal relevantes
- Implementar pruebas por paso y de extremo a extremo del pipeline que manejen idiosincrasias específicas de genómica de manera apropiada

<!-- TODO for future expansion: add metadata/samplesheet handling -->

## Requisitos previos

El curso asume cierta familiaridad mínima con lo siguiente:

- Herramientas y formatos de archivo comúnmente utilizados en este dominio científico
- Experiencia con la línea de comandos
- Conceptos fundamentales de Nextflow y herramientas cubiertas en el entrenamiento para principiantes [Hello Nextflow](../../hello_nextflow/)

Para requisitos técnicos y configuración del entorno, consulte el mini-curso [Configuración del Entorno](../../envsetup/).
