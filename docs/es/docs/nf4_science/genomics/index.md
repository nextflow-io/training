# Nextflow para Genómica

**Un curso práctico que aplica Nextflow a un caso de uso real de genómica: llamado de variantes con GATK.**

Este curso se basa en el entrenamiento para principiantes [Hello Nextflow](../../hello_nextflow/) y demuestra cómo usar Nextflow en el contexto específico del dominio de la genómica.
Implementará un pipeline de llamado de variantes con [GATK](https://gatk.broadinstitute.org/) (Genome Analysis Toolkit), un paquete de software ampliamente utilizado para analizar datos de secuenciación de alto rendimiento.

<!-- additional_information -->

## Descripción del curso

Este curso es práctico, con ejercicios orientados a objetivos estructurados para introducir información gradualmente.

Comenzará ejecutando las herramientas de llamado de variantes manualmente en la terminal para comprender la metodología, luego construirá progresivamente un pipeline de Nextflow que automatiza y escala el análisis.

### Plan de lecciones

Hemos dividido esto en tres partes que se enfocan en aspectos específicos de aplicar Nextflow a un caso de uso de genómica.

| Capítulo del curso                                                              | Resumen                                                                                                       | Duración estimada |
| ------------------------------------------------------------------------------- | ------------------------------------------------------------------------------------------------------------- | ----------------- |
| [Parte 1: Descripción del método](./01_method.md)                               | Comprender la metodología de llamado de variantes y ejecutar las herramientas manualmente                     | 30 min            |
| [Parte 2: Llamado de variantes por muestra](./02_per_sample_variant_calling.md) | Construir un pipeline que indexa archivos BAM y llama variantes, luego escalar a múltiples muestras           | 60 min            |
| [Parte 3: Llamado conjunto en una cohorte](./03_joint_calling.md)               | Agregar genotipado conjunto de múltiples muestras usando operadores de canal para agregar salidas por muestra | 45 min            |

Al finalizar este curso, podrá aplicar conceptos y herramientas fundamentales de Nextflow a un caso de uso típico de genómica.

¿Listo para tomar el curso?

[Comenzar :material-arrow-right:](00_orientation.md){ .md-button .md-button--primary }

<!-- Clearfix for float -->
<div style="content: ''; clear: both; display: table;"></div>
