---
title: Inicio
description: ¡Bienvenido/a al portal de capacitación de la comunidad Nextflow!
hide:
  - toc
  - footer
---

# Capacitación en Nextflow

<div class="grid cards" markdown>

-   :material-book-open-variant:{ .lg .middle } __Cursos de autoservicio__

    ---

    **¡Bienvenido/a al portal de capacitación de la comunidad Nextflow!**

    Avance a su propio ritmo a través de los cursos que se presentan a continuación, en nuestro entorno basado en la web o en el suyo propio.
    Cada curso es práctico, con ejercicios orientados a objetivos que puede completar de forma independiente.

    [Explorar los cursos :material-arrow-down:](#catalog-of-nextflow-training-courses){ .md-button .md-button--primary .mt-1 }

-   :material-account-group-outline:{ .lg .middle } __Eventos de capacitación__

    ---

    **¿Busca algo más allá del autoservicio?**

    Encuentre eventos de capacitación estructurados, orientación para organizar sus propios talleres y nuestra licencia de código abierto y política de contribución.

    [Ver eventos de capacitación :material-arrow-right:](training_events.md){ .md-button .md-button--secondary .mt-1 }

</div>

!!! note "Traducción asistida por IA"

    Esta traducción fue creada utilizando inteligencia artificial y revisada por traductores humanos.
    Agradecemos sus comentarios y sugerencias de mejora.
    Consulte nuestra [guía de traducción](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md) para más información.

## Catalog of Nextflow training courses

<div class="grid cards" markdown>

-   :material-account:{ .lg .middle } __Para usuarios__

    ---

    ### :material-play-circle:{.nextflow-primary} Ejecutar pipelines {.mt-1}

    Aprenda a ejecutar pipelines existentes sin escribir ningún código.

    ??? courses "**Nextflow Run:** Ejecutar pipelines con Nextflow"

        Una introducción rápida a la ejecución de pipelines de Nextflow que no requiere comprensión del código. Cubre el lanzamiento de pipelines, la recuperación de salidas, el uso de contenedores y la configuración de la ejecución a un nivel básico.

        [Ver la capacitación :material-arrow-right:](nextflow_run/index.md){ .md-button .md-button--secondary }

    ??? courses "**Use nf-core:** Encontrar y ejecutar pipelines curados por la comunidad"

        Una introducción rápida para encontrar, ejecutar y configurar pipelines del proyecto comunitario nf-core, comenzando con un pipeline de demostración mínimo y escalando hasta un pipeline de análisis a escala de producción.

        [Ver la capacitación :material-arrow-right:](nfcore_use/index.md){ .md-button .md-button--secondary }

    ??? courses "**Scale with Seqera:** Lanzar y monitorear pipelines a escala"

        Una introducción práctica al lanzamiento y monitoreo de pipelines de Nextflow con Seqera Platform, tanto desde la interfaz web como desde la línea de comandos.

        [Ver la capacitación :material-arrow-right:](seqera_scale/index.md){ .md-button .md-button--secondary }

    ---

    ### :material-tune:{.nextflow-primary} Gestionar la ejecución {.mt-1}

    Aprenda a gestionar la ejecución de pipelines de manera efectiva.

    ??? courses "**Execution Config:** Configurar pipelines como un profesional"

        Una introducción práctica a la configuración de la ejecución de pipelines de Nextflow: adaptación a diferentes entornos de cómputo, control de asignaciones de recursos y reintentos, y cambio entre perfiles de configuración predefinidos.

        [Ver la capacitación :material-arrow-right:](execution_config/index.md){ .md-button .md-button--secondary }

    !!! info compact "Más temas próximamente"

        La optimización del rendimiento, la ejecución en HPC/nube y más están planificados para esta sección.
        Vote sobre qué temas cubrir a continuación en nuestra [breve encuesta de interés](https://seqera.typeform.com/to/JCs91e8v).

-   :material-code-tags:{ .lg .middle } __Para desarrolladores__

    ---

    ### :material-wrench:{.nextflow-primary} Escribir pipelines {.mt-1}

    Aprenda a desarrollar sus propios pipelines de Nextflow.

    ??? courses "**Hello Nextflow:** Desarrollar sus propios pipelines desde cero"

        Este curso cubre los componentes principales del lenguaje Nextflow con suficiente detalle para permitir el desarrollo de pipelines simples pero completamente funcionales, además de elementos clave de diseño, desarrollo y prácticas de configuración de pipelines.

        [Ver la capacitación :material-arrow-right:](hello_nextflow/index.md){ .md-button .md-button--secondary }

    ??? courses "**Build with nf-core:** Usar las herramientas y reglas de nf-core"

        Para desarrolladores de Nextflow que deseen aprender a desarrollar pipelines compatibles con [nf-core](https://nf-co.re/).
        El curso cubre la estructura de los pipelines de nf-core con suficiente detalle para permitir el desarrollo de pipelines simples pero completamente funcionales que aprovechen la plantilla de nf-core y las mejores prácticas de desarrollo, así como el uso de módulos de nf-core existentes.

        [Ver la capacitación :material-arrow-right:](nfcore_build/index.md){ .md-button .md-button--secondary }

    ??? catalog "**Side Quests:** Profundizar en temas avanzados de Nextflow"

        Una colección de mini-cursos independientes destinados a desarrolladores de Nextflow que deseen ampliar su alcance y/o profundizar sus habilidades en temas particulares.
        Se presentan de forma lineal, pero pueden tomarse en cualquier orden (consulte las dependencias en la descripción general de cada mini-curso).

        [Explorar los Side Quests :material-arrow-right:](side_quests/index.md){ .md-button .md-button--secondary }

    ---

    ### :material-microscope:{.nextflow-primary} Nextflow para la Ciencia {.mt-1}

    Aprenda a desarrollar pipelines de Nextflow para aplicaciones científicas específicas.

    ??? courses "**Genomics:** Desarrollar un pipeline de llamado de variantes"

        Un curso para investigadores que deseen aprender a desarrollar sus propios pipelines de genómica, utilizando un caso de uso de llamado de variantes para demostrar los patrones esenciales de desarrollo en Nextflow.

        [Ver la capacitación :material-arrow-right:](nf4_science/genomics/index.md){ .md-button .md-button--secondary }

    ??? courses "**RNAseq:** Desarrollar un pipeline de procesamiento de RNAseq masivo"

        Un curso para investigadores que deseen aprender a desarrollar sus propios pipelines de RNAseq, utilizando un caso de uso de procesamiento de RNAseq masivo para demostrar los patrones esenciales de desarrollo en Nextflow.

        [Ver la capacitación :material-arrow-right:](nf4_science/rnaseq/index.md){ .md-button .md-button--secondary }

    ??? courses "**Bioimaging:** Ejecutar y configurar pipelines de imágenes"

        Un curso para investigadores que deseen aprender a ejecutar y configurar pipelines de bioimagen, utilizando nf-core/molkart para demostrar los patrones esenciales de uso de Nextflow.

        [Ver la capacitación :material-arrow-right:](nf4_science/imaging/index.md){ .md-button .md-button--secondary }

</div>

## Configuración y ayuda

<div class="grid cards mb-4" markdown>

-   :material-cog-outline:{ .lg .middle } __Entorno de capacitación__

    ---

    Opciones para configurar su entorno para las capacitaciones de Nextflow.

    [Ver los entornos de capacitación :material-arrow-right:](envsetup/index.md){ .md-button .md-button--secondary }

-   :material-tag-outline:{ .lg .middle } __Versiones de Nextflow__

    ---

    Comprensión y gestión de la evolución de las versiones de sintaxis de Nextflow.

    [Verificar los requisitos de versión :material-arrow-right:](info/nxf_versions.md){ .md-button .md-button--secondary }

-   :material-file-code-outline:{ .lg .middle } __El pipeline Hello__

    ---

    Resumen de lo que hace el pipeline Hello y cómo está estructurado.

    [Leer el resumen :material-arrow-right:](info/hello_pipeline.md){ .md-button .md-button--secondary }

-   :material-lifebuoy:{ .lg .middle } __Obtener ayuda__

    ---

    Recursos útiles cuando tenga algún problema con la capacitación de Nextflow.

    [Encontrar ayuda :material-arrow-right:](help.md){ .md-button .md-button--secondary }

</div>
