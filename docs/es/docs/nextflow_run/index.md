---
title: Nextflow Run
hide:
  - toc
page_type: index_page
index_type: course
additional_information:
  technical_requirements: true
  learning_objectives:
    - Lanzar y gestionar la ejecución de workflows de Nextflow
    - Encontrar e interpretar salidas (resultados) y archivos de registro
    - Reconocer los componentes principales de Nextflow en un workflow simple de múltiples pasos
    - Configurar la ejecución de pipelines para ejecutarse en plataformas de computación comunes, incluyendo HPC y nube
    - Resumir las mejores prácticas para reproducibilidad, portabilidad y reutilización de código que hacen los pipelines FAIR, incluyendo modularidad de código y contenedores de software
  audience_prerequisites:
    - "**Audiencia:** Este curso está diseñado para quienes son completamente nuevos en Nextflow y desean ejecutar pipelines existentes."
    - "**Habilidades:** Se asume cierta familiaridad con la línea de comandos, conceptos básicos de scripting y formatos de archivo comunes."
    - "**Dominio:** Los ejercicios son independientes del dominio, por lo que no se requiere conocimiento científico previo."
---

# Nextflow Run

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Traducción asistida por IA - [más información y sugerencias](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

**Nextflow Run es una introducción práctica a la ejecución de workflows de análisis de datos reproducibles y escalables.**

A través de ejemplos prácticos y ejercicios guiados, aprenderá los fundamentos del uso de Nextflow, incluyendo cómo ejecutar pipelines, gestionar archivos y dependencias de software, paralelizar la ejecución sin esfuerzo y ejecutar workflows en diferentes entornos de computación.

Se llevará las habilidades y la confianza para comenzar a ejecutar workflows con Nextflow.

<!-- additional_information -->

## Descripción general del curso

### Qué hará

Este curso es práctico, con ejercicios orientados a objetivos estructurados para introducir información gradualmente.

Ejecutará varias versiones de un pipeline de Nextflow que procesa entradas de texto.
Comenzará con una versión simple que consiste en un solo paso, y eventualmente progresará a una versión de múltiples pasos que toma un archivo CSV de entradas de texto tabulares, ejecuta algunos pasos de transformación y produce un único archivo de texto que contiene una imagen ASCII de un personaje diciendo el texto transformado.

Este curso se enfoca en ejecutar pipelines (nombrado así por el comando principal `nextflow run`).
Si busca una introducción al desarrollo de pipelines de Nextflow, consulte [Hello Nextflow](../hello_nextflow/index.md).

### Plan de lecciones

Hemos dividido esto en tres partes que se enfocarán en aspectos específicos de la ejecución y gestión de pipelines escritos en Nextflow.

| Capítulo del curso                                      | Resumen                                                                                                                       | Duración estimada |
| ------------------------------------------------------- | ----------------------------------------------------------------------------------------------------------------------------- | ----------------- |
| [Parte 1: Ejecutar operaciones básicas](./01_basics.md) | Lanzar y gestionar la ejecución de un workflow simple                                                                         | 30 mins           |
| [Parte 2: Ejecutar pipelines reales](./02_pipeline.md)  | Procesar entradas complejas, ejecutar workflows de múltiples pasos, usar contenedores y paralelizar la ejecución sin esfuerzo | 60 mins           |
| [Parte 3: Configuración de ejecución](./03_config.md)   | Personalizar el comportamiento del pipeline y optimizar el uso en diferentes entornos computacionales                         | 60 mins           |

Al final de este curso, estará bien preparado para abordar los próximos pasos en su camino para ejecutar workflows reproducibles para sus necesidades de computación científica.

¿Listo para tomar el curso?

[Comenzar a aprender :material-arrow-right:](00_orientation.md){ .md-button .md-button--primary }
