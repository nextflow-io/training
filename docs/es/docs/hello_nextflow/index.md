---
title: Hello Nextflow
hide:
  - toc
page_type: index_page
index_type: course
additional_information:
  technical_requirements: true
  learning_objectives:
    - Iniciar y gestionar la ejecución de flujos de trabajo de Nextflow
    - Encontrar e interpretar salidas (resultados) y archivos de registro generados por Nextflow
    - Solucionar problemas básicos
    - Construir un flujo de trabajo simple de múltiples pasos a partir de componentes principales de Nextflow
    - Distinguir entre tipos esenciales de channel factories y operadores y utilizarlos efectivamente en un flujo de trabajo simple
    - Configurar la ejecución de pipelines para ejecutar en plataformas de cómputo comunes incluyendo HPC y nube
    - Aplicar mejores prácticas de reproducibilidad, portabilidad y reutilización de código que hacen los pipelines FAIR, incluyendo modularidad del código y contenedores de software
  audience_prerequisites:
    - "**Audiencia:** Este curso está diseñado para estudiantes que son completamente nuevos en Nextflow y quieren desarrollar sus propios pipelines."
    - "**Habilidades:** Se asume cierta familiaridad con la línea de comandos, conceptos básicos de scripting y formatos de archivo comunes."
    - "**Dominio:** Los ejercicios son todos independientes del dominio, por lo que no se requiere conocimiento científico previo."
  videos_playlist: https://youtube.com/playlist?list=PLPZ8WHdZGxmWKozQuzr27jyMGqp9kElVK&si=eF7cLR62goy-lc6n
---

# Hello Nextflow

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Traducción asistida por IA - [más información y sugerencias](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

**Hello Nextflow es una introducción práctica a la construcción de flujos de trabajo de análisis de datos reproducibles y escalables.**

Trabajando a través de ejemplos prácticos y ejercicios guiados, aprenderá los fundamentos del desarrollo de pipelines con Nextflow, incluyendo cómo definir procesos, conectarlos en pipelines, gestionar archivos y dependencias de software, paralelizar la ejecución sin esfuerzo y ejecutar flujos de trabajo en diferentes entornos de cómputo.

Se llevará las habilidades y la confianza para comenzar a desarrollar y ejecutar sus propios flujos de trabajo con Nextflow.

<!-- additional_information -->

## Descripción general del curso

Este curso está diseñado para ser práctico, con ejercicios orientados a objetivos estructurados para introducir información gradualmente.

Desarrollará un pipeline simple de Nextflow que toma algunas entradas de texto, ejecuta algunos pasos de transformación, y produce un único archivo de texto que contiene una imagen ASCII de un personaje diciendo el texto transformado.

### Plan de lecciones

Para evitar abrumarle con conceptos y código, hemos dividido esto en seis partes que se enfocarán cada una en aspectos específicos del desarrollo de pipelines con Nextflow.

| Capítulo del curso                                    | Resumen                                                                                                          | Duración estimada |
| ----------------------------------------------------- | ---------------------------------------------------------------------------------------------------------------- | ----------------- |
| [Parte 1: Hello World](./01_hello_world.md)           | Componentes básicos y principios involucrados en ensamblar y ejecutar un flujo de trabajo de Nextflow            | 30 mins           |
| [Parte 2: Hello Channels](./02_hello_channels.md)     | Usar canales y operadores para procesar entradas y paralelizar la ejecución sin esfuerzo                         | 45 mins           |
| [Parte 3: Hello Workflow](./03_hello_workflow.md)     | Usar canales para encadenar múltiples pasos juntos y manejar la transferencia de datos entre pasos               | 60 mins           |
| [Parte 4: Hello Modules](./04_hello_modules.md)       | Aplicar principios de modularidad de código para aumentar la reutilización y disminuir la carga de mantenimiento | 20 mins           |
| [Parte 5: Hello Containers](./05_hello_containers.md) | Usar contenedores como mecanismo para gestionar dependencias de software y aumentar la reproducibilidad          | 60 mins           |
| [Parte 6: Hello Config](./06_hello_config.md)         | Personalizar el comportamiento del pipeline y optimizar el uso en diferentes entornos computacionales            | 60 mins           |

Al final de este curso, estará bien preparado para abordar los próximos pasos en su viaje para desarrollar flujos de trabajo reproducibles para sus necesidades de computación científica.

¿Listo para tomar el curso?

[Comenzar :material-arrow-right:](00_orientation.md){ .md-button .md-button--primary }

<!-- Clearfix for float -->
<div style="content: ''; clear: both; display: table;"></div>
