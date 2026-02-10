---
title: Inicio
description: ¡Te damos la bienvenida al portal de capacitación de la comunidad Nextflow!
hide:
  - toc
  - footer
---

# Capacitación Nextflow

<div class="grid cards" markdown>

-   :material-book-open-variant:{ .lg .middle } __Cursos de autoaprendizaje__

    ---

    **¡Te damos la bienvenida al portal de capacitación de la comunidad Nextflow!**

    Los cursos de capacitación que se enumeran a continuación están diseñados para ser utilizados como un recurso de autoaprendizaje.
    Puedes trabajar en ellos por tu cuenta en cualquier momento, ya sea en el entorno basado en web que proporcionamos a través de Github Codespaces o en tu propio entorno.

    [Explorar los cursos :material-arrow-right:](#catalog-of-nextflow-training-courses){ .md-button .md-button--primary .mt-1 }

-   :material-information-outline:{ .lg .middle } __Información adicional__

    ---

    ??? warning "Compatibilidad de versiones"

        <!-- Any update to this content needs to be copied to the local installation page -->
        **A partir de enero de 2026, todos nuestros cursos de capacitación Nextflow requieren Nextflow versión 25.10.2 o posterior, con sintaxis estricta activada, a menos que se indique lo contrario.**

        Para obtener más información sobre los requisitos de versión y la sintaxis estricta, consulta la [guía de migración de la documentación de Nextflow](https://nextflow.io/docs/latest/strict-syntax.html).

        Las versiones anteriores del material de capacitación correspondientes a la sintaxis previa están disponibles a través del selector de versiones en la barra de menú de esta página web.

    ??? terminal "Opciones de entorno"

        Proporcionamos un entorno de capacitación basado en web donde todo lo que necesitas para tomar la capacitación está preinstalado, disponible a través de Github Codespaces (requiere una cuenta gratuita de GitHub).

        [![Abrir en GitHub Codespaces](https://github.com/codespaces/badge.svg)](https://codespaces.new/nextflow-io/training?quickstart=1&ref=master)

        Si esto no se ajusta a tus necesidades, consulta las otras [Opciones de entorno](./envsetup/index.md).

    ??? learning "Eventos de capacitación"

        Si prefieres tomar la capacitación de Nextflow como parte de un evento estructurado, hay muchas oportunidades para hacerlo. Te recomendamos consultar las siguientes opciones:

        - **[Training Weeks]()** organizadas trimestralmente por el equipo de la Comunidad
        - **[Seqera Events](https://seqera.io/events/)** incluyen eventos de capacitación presenciales organizados por Seqera (busca 'Seqera Sessions' y 'Nextflow Summit')
        - **[Nextflow Ambassadors]()** organizan eventos para su comunidad local
        - **[nf-core events](https://nf-co.re/events)** incluyen hackathons de la comunidad

    ??? people "Información para instructores"

        Si eres un instructor que organiza sus propias capacitaciones, puedes utilizar nuestros materiales directamente desde el portal de capacitación siempre que atribuyas el crédito apropiado. Consulta 'Créditos y contribuciones' a continuación para más detalles.

        Además, ¡nos encantaría saber de ti sobre cómo podríamos apoyar mejor tus esfuerzos de capacitación! Contáctanos en [community@seqera.io](mailto:community@seqera.io) o en el foro de la comunidad (consulta la página de [Ayuda](help.md)).

    ??? licensing "Licencia de código abierto y política de contribución"

        [![Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International (CC BY-NC-SA 4.0)](assets/img/cc_by-nc-sa.svg){ align=right }](https://creativecommons.org/licenses/by-nc-sa/4.0/)

        Este material de capacitación es desarrollado y mantenido por [Seqera](https://seqera.io) y publicado bajo una licencia de código abierto ([CC BY-NC-SA](https://creativecommons.org/licenses/by-nc-sa/4.0/)) para el beneficio de la comunidad. Si deseas utilizar este material de una manera que esté fuera del alcance de la licencia (ten en cuenta las limitaciones sobre el uso comercial y la redistribución), contáctanos en [community@seqera.io](mailto:community@seqera.io) para discutir tu solicitud.

        Damos la bienvenida a mejoras, correcciones y reportes de errores de la comunidad. Cada página tiene un ícono :material-file-edit-outline: en la parte superior derecha de la página que enlaza al repositorio de código, donde puedes reportar problemas o proponer cambios al material fuente de capacitación a través de un pull request. Consulta el `README.md` en el repositorio para más detalles.

</div>

!!! note "Traducción asistida por IA"

    Esta traducción fue creada utilizando inteligencia artificial y revisada por traductores humanos.
    Agradecemos sus comentarios y sugerencias de mejora.
    Consulte nuestra [guía de traducción](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md) para más información.

## Catálogo de cursos de capacitación Nextflow

<div class="grid cards" markdown>

-   :material-walk:{ .lg .middle } __Ruta introductoria__

    ---

    ### :material-compass:{.nextflow-primary} Nextflow para Principiantes {.mt-1}

    Cursos independientes del dominio destinados a aquellos que son completamente nuevos en Nextflow. Cada curso consiste en una serie de módulos de capacitación diseñados para ayudar a los estudiantes a desarrollar sus habilidades progresivamente.

    ??? courses "**Hello Nextflow:** Aprende a desarrollar tus propios pipelines"

        Este curso cubre los componentes principales del lenguaje Nextflow con suficiente detalle para permitir el desarrollo de pipelines simples pero completamente funcionales, además de elementos clave de diseño, desarrollo y prácticas de configuración de pipelines.

        [Comenzar la capacitación Hello Nextflow :material-arrow-right:](hello_nextflow/index.md){ .md-button .md-button--secondary }

    ??? courses "**Nextflow Run:** Aprende a ejecutar pipelines existentes"

        Una introducción concisa a la ejecución y configuración de pipelines Nextflow, basada en el curso para desarrolladores Hello Nextflow pero con menos enfoque en el código. Cubre ejecución, salidas, estructura básica del código y configuración para diferentes entornos de cómputo.

        [Comenzar la capacitación Nextflow Run :material-arrow-right:](nextflow_run/index.md){ .md-button .md-button--secondary }

    ---

    ### :material-microscope:{.nextflow-primary} Nextflow para Ciencia {.mt-1}

    Aprende a aplicar los conceptos y componentes presentados en 'Hello Nextflow' a casos de uso científicos específicos.

    ??? courses "**Nextflow para Genómica** (llamado de variantes)"

        Para investigadores que desean aprender a desarrollar sus propios pipelines de genómica. El curso utiliza un caso de uso de llamado de variantes para demostrar cómo desarrollar un pipeline de genómica simple pero funcional.

        [Comenzar la capacitación Nextflow para Genómica :material-arrow-right:](nf4_science/genomics/){ .md-button .md-button--secondary }

    ??? courses "**Nextflow para RNAseq** (RNAseq bulk)"

        Para investigadores que desean aprender a desarrollar sus propios pipelines de RNAseq. El curso utiliza un caso de uso de procesamiento de RNAseq bulk para demostrar cómo desarrollar un pipeline de RNAseq simple pero funcional.

        [Comenzar la capacitación Nextflow para RNAseq :material-arrow-right:](nf4_science/rnaseq/){ .md-button .md-button--secondary }

    ??? courses "**Nextflow para Imágenes** (ómicas espaciales)"

        Para investigadores en imágenes y ómicas espaciales que desean aprender a ejecutar y personalizar pipelines de análisis. El curso utiliza el pipeline nf-core/molkart para proporcionar un pipeline biológicamente relevante que demuestra cómo ejecutar, configurar y gestionar entradas para workflows de pipelines Nextflow.

        [Comenzar la capacitación Nextflow para Imágenes :material-arrow-right:](nf4_science/imaging/){ .md-button .md-button--secondary }

-   :material-run:{ .lg .middle } __Ruta avanzada__

    ---

    ### :material-bridge:{.nextflow-primary} De Nextflow a nf-core {.mt-1}

    Aprende a utilizar código y mejores prácticas del proyecto comunitario [nf-core](https://nf-co.re/).

    Estos cursos te ayudan a pasar de los fundamentos de Nextflow a las mejores prácticas de nf-core.
    Comprende cómo y por qué la comunidad nf-core construye pipelines, y cómo puedes contribuir y reutilizar estas técnicas.

    ??? courses "**Hello nf-core:** Comienza con nf-core"

        Para desarrolladores que desean aprender a ejecutar y desarrollar pipelines compatibles con [nf-core](https://nf-co.re/). El curso cubre la estructura de los pipelines nf-core con suficiente detalle para permitir el desarrollo de pipelines simples pero completamente funcionales que siguen la plantilla nf-core y las mejores prácticas de desarrollo, así como el uso de módulos nf-core existentes.

        [Comenzar la capacitación Hello nf-core :material-arrow-right:](hello_nf-core/index.md){ .md-button .md-button--secondary }

    ---

    ### :material-rocket-launch:{.nextflow-primary} Capacitación Avanzada en Nextflow {.mt-1}

    Aprende conceptos y mecanismos avanzados para desarrollar e implementar pipelines Nextflow para abordar casos de uso del mundo real.

    ??? courses "**Side Quests:** Inmersiones profundas en temas independientes"

        Mini-cursos independientes destinados a desarrolladores de Nextflow que desean ampliar su rango y/o profundizar sus habilidades en temas particulares. Se presentan linealmente pero pueden tomarse en cualquier orden (consulta las dependencias en cada descripción general del mini-curso).

        [Explorar los Side Quests :material-arrow-right:](side_quests/){ .md-button .md-button--secondary }

    ??? courses "**Training Collections:** Rutas de aprendizaje recomendadas a través de los Side Quests"

        Las Training Collections combinan múltiples Side Quests para proporcionar una experiencia de aprendizaje integral en torno a un tema o caso de uso particular.

        [Explorar las Training Collections :material-arrow-right:](training_collections/){ .md-button .md-button--secondary }

</div>

!!! info "¿Buscas materiales de capacitación archivados?"

    Los materiales de capacitación más antiguos (Fundamentals Training, Advanced Training y otros cursos experimentales) han sido eliminados del portal de capacitación ya que son incompatibles con la sintaxis estricta de Nextflow 3.0.
    Si necesitas acceso a estos materiales, están disponibles en el [historial de git](https://github.com/nextflow-io/training) antes de enero de 2026.

---

<div markdown class="homepage_logos">

![Seqera](assets/img/seqera_logo.png#only-light)

![Seqera](assets/img/seqera_logo_dark.png#only-dark)

</div>
