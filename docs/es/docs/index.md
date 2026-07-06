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

    Los cursos de capacitación que se enumeran a continuación están diseñados para usarse como un recurso de autoservicio.
    Puede trabajar en ellos por su cuenta en cualquier momento, ya sea en el entorno basado en web que proporcionamos a través de Github Codespaces o en su propio entorno.

    [Explorar los cursos :material-arrow-right:](#catalog-of-nextflow-training-courses){ .md-button .md-button--primary .mt-1 }

-   :material-information-outline:{ .lg .middle } __Información adicional__

    ---

    ??? warning "Compatibilidad de versiones"

        <!-- Cualquier actualización de este contenido debe copiarse en la página de instalación local -->
        **A partir de enero de 2026, todos nuestros cursos de capacitación de Nextflow requieren la versión 25.10.2 de Nextflow o posterior, con la sintaxis estricta activada, a menos que se indique lo contrario.**

        Para obtener más información sobre los requisitos de versión y la sintaxis estricta, consulte la [guía de migración de la documentación de Nextflow](https://nextflow.io/docs/latest/strict-syntax.html).

        Las versiones anteriores del material de capacitación correspondientes a la sintaxis previa están disponibles a través del selector de versiones en la barra de menú de esta página web.

    ??? terminal "Opciones de entorno"

        Proporcionamos un entorno de capacitación basado en web donde todo lo que necesita para tomar la capacitación está preinstalado, disponible a través de Github Codespaces (requiere una cuenta gratuita de GitHub).

        [![Open in GitHub Codespaces](https://github.com/codespaces/badge.svg)](https://codespaces.new/nextflow-io/training?quickstart=1&ref=master)

        Si esto no se adapta a sus necesidades, consulte las otras [opciones de entorno](./envsetup/index.md).

    ??? learning "Eventos de capacitación"

        Si prefiere tomar la capacitación de Nextflow como parte de un evento estructurado, hay muchas oportunidades para hacerlo. Recomendamos consultar las siguientes opciones:

        - **[Training Weeks]()** organizadas trimestralmente por el equipo de la Comunidad
        - **[Seqera Events](https://seqera.io/events/)** incluyen eventos de capacitación presenciales organizados por Seqera (busque 'Seqera Sessions' y 'Nextflow Summit')
        - **[Nextflow Ambassadors]()** organizan eventos para su comunidad local
        - **[nf-core events](https://nf-co.re/events)** incluyen hackatones comunitarios

    ??? people "Información para instructores"

        Si es un instructor que organiza sus propias capacitaciones, puede usar nuestros materiales directamente desde el portal de capacitación siempre que otorgue el crédito correspondiente. Consulte 'Créditos y contribuciones' a continuación para más detalles.

        Además, ¡nos encantaría saber cómo podríamos apoyar mejor sus esfuerzos de capacitación! Contáctenos en [community@seqera.io](mailto:community@seqera.io) o en el foro de la comunidad (consulte la página de [Ayuda](help.md)).

    ??? licensing "Licencia de código abierto y política de contribución"

        [![Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International (CC BY-NC-SA 4.0)](assets/img/cc_by-nc-sa.svg){ align=right }](https://creativecommons.org/licenses/by-nc-sa/4.0/)

        Este material de capacitación es desarrollado y mantenido por [Seqera](https://seqera.io) y publicado bajo una licencia de código abierto ([CC BY-NC-SA](https://creativecommons.org/licenses/by-nc-sa/4.0/)) para el beneficio de la comunidad. Si desea utilizar este material de una manera que esté fuera del alcance de la licencia (tenga en cuenta las limitaciones sobre uso comercial y redistribución), contáctenos en [community@seqera.io](mailto:community@seqera.io) para analizar su solicitud.

        Agradecemos mejoras, correcciones e informes de errores de la comunidad. Cada página tiene un ícono :material-file-edit-outline: en la parte superior derecha que enlaza al repositorio de código, donde puede reportar problemas o proponer cambios al material fuente de capacitación mediante un pull request. Consulte el `README.md` en el repositorio para más detalles.

</div>

!!! note "Traducción asistida por IA"

    Esta traducción fue creada utilizando inteligencia artificial y revisada por traductores humanos.
    Agradecemos sus comentarios y sugerencias de mejora.
    Consulte nuestra [guía de traducción](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md) para más información.

## Catálogo de cursos de capacitación de Nextflow

<div class="grid cards" markdown>

-   :material-walk:{ .lg .middle } __Nivel introductorio__

    ---

    ### :material-compass:{.nextflow-primary} Nextflow para principiantes {.mt-1}

    Cursos independientes del dominio destinados a quienes son completamente nuevos en Nextflow. Cada curso consiste en una serie de módulos de capacitación diseñados para ayudar a los participantes a desarrollar sus habilidades de forma progresiva.

    ??? courses "**Hello Nextflow:** Aprenda a desarrollar sus propios pipelines"

        Este curso cubre los componentes principales del lenguaje Nextflow con suficiente detalle para permitir el desarrollo de pipelines simples pero completamente funcionales, además de elementos clave de diseño, desarrollo y prácticas de configuración de pipelines.

        [Comenzar la capacitación Hello Nextflow :material-arrow-right:](hello_nextflow/index.md){ .md-button .md-button--secondary }

    ??? courses "**Nextflow Run:** Aprenda a ejecutar pipelines existentes"

        Una introducción concisa para ejecutar y configurar pipelines de Nextflow, basada en el curso para desarrolladores Hello Nextflow pero con menos énfasis en el código. Cubre la ejecución, las salidas, la estructura básica del código y la configuración para diferentes entornos de cómputo.

        [Comenzar la capacitación Nextflow Run :material-arrow-right:](nextflow_run/index.md){ .md-button .md-button--secondary }

    ---

    ### :material-microscope:{.nextflow-primary} Nextflow para la ciencia {.mt-1}

    Aprenda a aplicar los conceptos y componentes presentados en 'Hello Nextflow' a casos de uso científicos específicos.

    ??? courses "**Nextflow para Genómica** (llamado de variantes)"

        Para investigadores que deseen aprender a desarrollar sus propios pipelines de genómica. El curso utiliza un caso de uso de llamado de variantes para demostrar cómo desarrollar un pipeline de genómica simple pero funcional.

        [Comenzar la capacitación Nextflow para Genómica :material-arrow-right:](nf4_science/genomics/index.md){ .md-button .md-button--secondary }

    ??? courses "**Nextflow para RNAseq** (RNAseq masivo)"

        Para investigadores que deseen aprender a desarrollar sus propios pipelines de RNAseq. El curso utiliza un caso de uso de procesamiento de RNAseq masivo para demostrar cómo desarrollar un pipeline de RNAseq simple pero funcional.

        [Comenzar la capacitación Nextflow para RNAseq :material-arrow-right:](nf4_science/rnaseq/index.md){ .md-button .md-button--secondary }

    ??? courses "**Nextflow para Imágenes** (ómica espacial)"

        Para investigadores en imágenes y ómica espacial que deseen aprender a ejecutar y personalizar pipelines de análisis. El curso utiliza el pipeline nf-core/molkart para proporcionar un pipeline biológicamente relevante que demuestra cómo ejecutar, configurar y gestionar entradas para workflows de pipelines de Nextflow.

        [Comenzar la capacitación Nextflow para Imágenes :material-arrow-right:](nf4_science/imaging/index.md){ .md-button .md-button--secondary }

-   :material-run:{ .lg .middle } __Nivel avanzado__

    ---

    ### :material-bridge:{.nextflow-primary} De Nextflow a nf-core {.mt-1}

    Aprenda a utilizar código y buenas prácticas del proyecto comunitario [nf-core](https://nf-co.re/).

    Estos cursos le ayudan a pasar de los fundamentos de Nextflow a las buenas prácticas de nf-core.
    Comprenda cómo y por qué la comunidad nf-core construye pipelines, y cómo puede contribuir y reutilizar estas técnicas.

    ??? courses "**Hello nf-core:** Primeros pasos con nf-core"

        Para desarrolladores que deseen aprender a ejecutar y desarrollar pipelines compatibles con [nf-core](https://nf-co.re/). El curso cubre la estructura de los pipelines de nf-core con suficiente detalle para permitir el desarrollo de pipelines simples pero completamente funcionales que sigan la plantilla y las buenas prácticas de desarrollo de nf-core, así como el uso de módulos nf-core existentes.

        [Comenzar la capacitación Hello nf-core :material-arrow-right:](hello_nf-core/index.md){ .md-button .md-button--secondary }

    ---

    ### :material-rocket-launch:{.nextflow-primary} Capacitación avanzada en Nextflow {.mt-1}

    Aprenda conceptos y mecanismos avanzados para desarrollar e implementar pipelines de Nextflow que aborden casos de uso del mundo real.

    ??? courses "**Side Quests:** Exploraciones en profundidad de temas independientes"

        Mini-cursos independientes destinados a desarrolladores de Nextflow que deseen ampliar su alcance y/o profundizar sus habilidades en temas particulares. Se presentan de forma lineal pero pueden tomarse en cualquier orden (consulte las dependencias en la descripción general de cada mini-curso).

        [Explorar los Side Quests :material-arrow-right:](side_quests/index.md){ .md-button .md-button--secondary }

    ??? courses "**Training Collections:** Rutas de aprendizaje recomendadas a través de los Side Quests"

        Las Training Collections combinan múltiples Side Quests para proporcionar una experiencia de aprendizaje integral en torno a un tema o caso de uso particular.

        [Explorar las Training Collections :material-arrow-right:](training_collections/index.md){ .md-button .md-button--secondary }

</div>

!!! info "¿Busca materiales de capacitación archivados?"

    Los materiales de capacitación anteriores (Fundamentals Training, Advanced Training y otros cursos experimentales) han sido eliminados del portal de capacitación ya que son incompatibles con la sintaxis estricta de Nextflow 3.0.
    Si necesita acceder a estos materiales, están disponibles en el [historial de git](https://github.com/nextflow-io/training) previo a enero de 2026.

---

<div markdown class="homepage_logos">

![Seqera](assets/img/seqera_logo.png#only-light)

![Seqera](assets/img/seqera_logo_dark.png#only-dark)

</div>
