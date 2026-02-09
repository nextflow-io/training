---
title: Nextflow run for Imaging
hide:
  - toc
---

# Nextflow run for Imaging

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Traducción asistida por IA - [más información y sugerencias](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

Este curso de capacitación está dirigido a investigadores en imagenología y biología espacial que estén interesados en ejecutar y personalizar pipelines de análisis de datos.
Enseña conceptos fundamentales de Nextflow relacionados con la ejecución, organización y configuración de workflows utilizando [nf-core/molkart](https://nf-co.re/molkart), un pipeline para procesar datos de transcriptómica espacial de Molecular Cartography.
Las habilidades que aprenderás aquí son transferibles a cualquier pipeline de Nextflow o nf-core.

¡Comencemos! Haz clic en el botón "Open in GitHub Codespaces" a continuación para iniciar el entorno de capacitación (preferiblemente en una pestaña separada), luego continúa leyendo mientras se carga.

[![Open in GitHub Codespaces](https://github.com/codespaces/badge.svg)](https://codespaces.new/nextflow-io/training?quickstart=1&ref=master)

## Objetivos de aprendizaje

Al trabajar en este curso, aprenderás a aplicar conceptos y herramientas fundamentales de Nextflow para ejecutar pipelines de análisis de imagenología.

Al finalizar este taller serás capaz de:

- Lanzar un workflow de Nextflow localmente y monitorear su ejecución
- Encontrar e interpretar las salidas (resultados) y archivos de registro generados por Nextflow
- Ejecutar un pipeline de nf-core con datos de prueba y entradas personalizadas
- Configurar la ejecución del pipeline usando perfiles y archivos de parámetros
- Gestionar entradas usando samplesheets y parámetros de línea de comandos

## Audiencia y requisitos previos

Este curso asume cierta familiaridad mínima con lo siguiente:

- Experiencia con la línea de comandos
- Familiaridad básica con formatos de archivos de imagenología (imágenes TIFF, datos tabulares)

Para los requisitos técnicos y la configuración del entorno, consulta el mini-curso de [Configuración del entorno](../../envsetup/).
