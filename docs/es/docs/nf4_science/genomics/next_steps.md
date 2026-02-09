# Resumen del curso

¡Felicidades por completar el curso de capacitación de Nextflow para Genómica! 🎉

## Tu recorrido

Comenzaste ejecutando herramientas de llamado de variantes manualmente en la terminal para comprender la metodología.
Luego construiste un pipeline de Nextflow para una sola muestra para automatizar el proceso, lo escalaste para manejar múltiples muestras en paralelo y agregaste genotipado conjunto de múltiples muestras usando operadores de canal.

### Lo que construiste

- Un pipeline de llamado de variantes que toma archivos BAM como entrada y produce VCFs con llamado conjunto como salida.
- Tres procesos (`SAMTOOLS_INDEX`, `GATK_HAPLOTYPECALLER` y `GATK_JOINTGENOTYPING`) almacenados en archivos de módulo separados.
- El pipeline escala automáticamente a cualquier número de muestras de entrada usando el paradigma de flujo de datos de Nextflow.
- Los resultados se publican en un directorio llamado `results/`.

### Habilidades adquiridas

A través de este curso práctico, has aprendido cómo:

- Escribir un workflow lineal para aplicar llamado de variantes a una sola muestra
- Manejar archivos accesorios como archivos de índice y recursos del genoma de referencia apropiadamente
- Aprovechar el paradigma de flujo de datos de Nextflow para paralelizar el llamado de variantes por muestra
- Implementar llamado conjunto de múltiples muestras usando operadores de canal relevantes
  Ahora estás equipado/a para comenzar a aplicar Nextflow a workflows de análisis genómico en tu propio trabajo.

## Próximos pasos para desarrollar tus habilidades

Aquí están nuestras principales sugerencias sobre qué hacer a continuación:

- Aplica Nextflow a otros casos de uso de análisis científico con [Nextflow for Science](../index.md)
- Comienza con nf-core con [Hello nf-core](../../hello_nf-core/index.md)
- Explora características más avanzadas de Nextflow con las [Side Quests](../../side_quests/index.md)

Finalmente, te recomendamos que eches un vistazo a [**Seqera Platform**](https://seqera.io/), una plataforma basada en la nube desarrollada por los creadores de Nextflow que hace aún más fácil lanzar y gestionar tus workflows, así como administrar tus datos y ejecutar análisis de forma interactiva en cualquier entorno.

## Obtener ayuda

Para recursos de ayuda y soporte de la comunidad, consulta la [página de Ayuda](../../help.md).

## Encuesta de retroalimentación

Antes de continuar, ¡por favor toma un minuto para completar la encuesta del curso! Tu retroalimentación nos ayuda a mejorar nuestros materiales de capacitación para todos.

[Realizar la encuesta :material-arrow-right:](survey.md){ .md-button .md-button--primary }
