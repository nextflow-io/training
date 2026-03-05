# Resumen del curso

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Traducción asistida por IA - [más información y sugerencias](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

¡Felicitaciones por completar el curso de entrenamiento de Nextflow para RNAseq!

## Tu recorrido

Comenzaste ejecutando herramientas de procesamiento de RNAseq manualmente en la terminal para comprender la metodología.
Luego construiste un pipeline de Nextflow para una sola muestra para automatizar el proceso, lo escalaste para manejar múltiples muestras en paralelo y lo extendiste para manejar datos paired-end y agregar reportes de QC entre muestras.

### Lo que construiste

- Un pipeline de procesamiento de RNAseq que toma archivos FASTQ como entrada y produce lecturas recortadas, alineamientos y reportes de QC agregados como salida.
- Procesos para recorte (Trim Galore), alineamiento (HISAT2), control de calidad (FastQC) y agregación de reportes (MultiQC) almacenados en archivos de módulos separados.
- El pipeline paraleliza automáticamente el procesamiento de muestras de entrada usando el paradigma de flujo de datos de Nextflow.
- El pipeline final maneja datos de secuenciación paired-end.

### Habilidades adquiridas

A través de este curso práctico, has aprendido cómo:

- Escribir un workflow lineal para aplicar métodos básicos de procesamiento y QC de RNAseq
- Manejar archivos específicos del dominio como FASTQ y recursos de genoma de referencia apropiadamente
- Manejar datos de secuenciación single-end y paired-end
- Aprovechar el paradigma de flujo de datos de Nextflow para paralelizar el procesamiento de RNAseq por muestra
- Agregar reportes de QC a través de múltiples pasos y muestras usando operadores de canal relevantes

Ahora estás equipado para comenzar a aplicar Nextflow a workflows de análisis de RNAseq en tu propio trabajo.

## Próximos pasos para desarrollar tus habilidades

Aquí están nuestras principales sugerencias sobre qué hacer a continuación:

- Aplica Nextflow a otros casos de uso de análisis científico con [Nextflow para Ciencia](../index.md)
- Comienza con nf-core con [Hello nf-core](../../hello_nf-core/index.md)
- Explora características más avanzadas de Nextflow con las [Misiones Secundarias](../../side_quests/index.md)

Finalmente, te recomendamos que eches un vistazo a [**Seqera Platform**](https://seqera.io/), una plataforma basada en la nube desarrollada por los creadores de Nextflow que facilita aún más el lanzamiento y gestión de tus workflows, así como gestionar tus datos y ejecutar análisis de forma interactiva en cualquier entorno.

## Obtener ayuda

Para recursos de ayuda y soporte de la comunidad, consulta la [página de Ayuda](../../help.md).

## Encuesta de retroalimentación

Antes de continuar, ¡por favor toma un minuto para completar la encuesta del curso! Tu retroalimentación nos ayuda a mejorar nuestros materiales de capacitación para todos.

[Realizar la encuesta :material-arrow-right:](survey.md){ .md-button .md-button--primary }
