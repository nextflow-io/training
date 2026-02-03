---
title: Inicio
description: ¡Bienvenido al portal de entrenamiento de la comunidad Nextflow!
hide:
  - toc
  - footer
---

# Entrenamiento de Nextflow

<div class="grid cards" markdown>

-   :material-book-open-variant:{ .lg .middle } __Cursos de autoservicio__

    ---

    **¡Bienvenido al portal de entrenamiento de la comunidad Nextflow!**

    Los cursos de entrenamiento listados a continuación están diseñados para ser utilizados como recurso de autoservicio.
    Puede trabajar con ellos por su cuenta en cualquier momento, ya sea en el entorno web que proporcionamos a través de Github Codespaces o en su propio entorno.

    [Explorar los cursos :material-arrow-right:](#catalogo-de-cursos-de-entrenamiento-de-nextflow){ .md-button .md-button--primary .mt-1 }

-   :material-information-outline:{ .lg .middle } __Información adicional__

    ---

    ??? warning "Compatibilidad de versiones"

        <!-- Cualquier actualización de este contenido debe copiarse a la página de instalación local -->
        **A partir de enero de 2026, todos nuestros cursos de entrenamiento de Nextflow requieren la versión 25.10.2 o posterior de Nextflow, con la sintaxis estricta activada, a menos que se indique lo contrario.**

        Para más información sobre los requisitos de versión y la sintaxis estricta, consulte la [guía de migración de la documentación de Nextflow](https://nextflow.io/docs/latest/strict-syntax.html).

        Las versiones anteriores del material de entrenamiento correspondientes a la sintaxis anterior están disponibles a través del selector de versiones en la barra de menú de esta página web.

    ??? terminal "Opciones de entorno"

        Proporcionamos un entorno de entrenamiento basado en web donde todo lo que necesita para realizar el entrenamiento está preinstalado, disponible a través de Github Codespaces (requiere una cuenta gratuita de GitHub).

        [![Open in GitHub Codespaces](https://github.com/codespaces/badge.svg)](https://codespaces.new/nextflow-io/training?quickstart=1&ref=master)

        Si esto no se adapta a sus necesidades, consulte las otras [Opciones de entorno](./envsetup/index.md).

    ??? learning "Eventos de entrenamiento"

        Si prefiere tomar el entrenamiento de Nextflow como parte de un evento estructurado, hay muchas oportunidades para hacerlo. Le recomendamos consultar las siguientes opciones:

        - **[Semanas de Entrenamiento]()** organizadas trimestralmente por el equipo de la Comunidad
        - **[Eventos de Seqera](https://seqera.io/events/)** incluyen eventos de entrenamiento presencial organizados por Seqera (busque 'Seqera Sessions' y 'Nextflow Summit')
        - **[Embajadores de Nextflow]()** organizan eventos para su comunidad local
        - **[Eventos de nf-core](https://nf-co.re/events)** incluyen hackathons de la comunidad

    ??? people "Información para instructores"

        Si usted es un instructor que dirige sus propios entrenamientos, puede utilizar nuestros materiales directamente desde el portal de entrenamiento siempre que atribuya el crédito correspondiente. Consulte 'Créditos y contribuciones' a continuación para más detalles.

        Además, nos encantaría saber de usted cómo podríamos apoyar mejor sus esfuerzos de entrenamiento. Por favor contáctenos en [community@seqera.io](mailto:community@seqera.io) o en el foro de la comunidad (consulte la página de [Ayuda](help.md)).

    ??? licensing "Licencia de código abierto y política de contribuciones"

        [![Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International (CC BY-NC-SA 4.0)](assets/img/cc_by-nc-sa.svg){ align=right }](https://creativecommons.org/licenses/by-nc-sa/4.0/)

        Este material de entrenamiento está desarrollado y mantenido por [Seqera](https://seqera.io) y publicado bajo una licencia de código abierto ([CC BY-NC-SA](https://creativecommons.org/licenses/by-nc-sa/4.0/)) para beneficio de la comunidad. Si desea utilizar este material de una manera que quede fuera del alcance de la licencia (tenga en cuenta las limitaciones sobre uso comercial y redistribución), por favor contáctenos en [community@seqera.io](mailto:community@seqera.io) para discutir su solicitud.

        Damos la bienvenida a mejoras, correcciones y reportes de errores de la comunidad. Cada página tiene un ícono :material-file-edit-outline: en la parte superior derecha de la página que enlaza al repositorio de código, donde puede reportar problemas o proponer cambios al material de entrenamiento fuente a través de un pull request. Consulte el `README.md` en el repositorio para más detalles.

</div>

!!! note "Traducción asistida por IA"

    Esta traducción fue creada utilizando inteligencia artificial y revisada por traductores humanos.
    Agradecemos tus comentarios y sugerencias de mejora.
    Consulta nuestra [guía de traducción](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md) para más información.

## Catálogo de cursos de entrenamiento de Nextflow

<div class="grid cards" markdown>

-   :material-walk:{ .lg .middle } __Pista introductoria__

    ---

    ### :material-compass:{.nextflow-primary} Nextflow para principiantes {.mt-1}

    Cursos independientes del dominio destinados a quienes son completamente nuevos en Nextflow. Cada curso consiste en una serie de módulos de entrenamiento diseñados para ayudar a los estudiantes a desarrollar sus habilidades progresivamente.

    ??? courses "**Hello Nextflow:** Aprenda a desarrollar sus propios pipelines"

        Este curso cubre los componentes principales del lenguaje Nextflow con suficiente detalle para permitir el desarrollo de pipelines simples pero completamente funcionales, además de elementos clave de diseño de pipelines, prácticas de desarrollo y configuración.

        [Comenzar el entrenamiento Hello Nextflow :material-arrow-right:](hello_nextflow/index.md){ .md-button .md-button--secondary }

    ??? courses "**Nextflow Run:** Aprenda a ejecutar pipelines existentes"

        Una introducción concisa a la ejecución y configuración de pipelines de Nextflow, basada en el curso de desarrolladores Hello Nextflow pero con menos enfoque en el código. Cubre la ejecución, salidas, estructura básica del código y configuración para diferentes entornos de cómputo.

        [Comenzar el entrenamiento Nextflow Run :material-arrow-right:](nextflow_run/index.md){ .md-button .md-button--secondary }

    ---

    ### :material-microscope:{.nextflow-primary} Nextflow para Ciencia {.mt-1}

    Aprenda a aplicar los conceptos y componentes presentados en 'Hello Nextflow' a casos de uso científicos específicos.

    ??? courses "**Nextflow para Genómica** (llamado de variantes)"

        Para investigadores que desean aprender a desarrollar sus propios pipelines de genómica. El curso utiliza un caso de uso de llamado de variantes para demostrar cómo desarrollar un pipeline de genómica simple pero funcional.

        [Comenzar el entrenamiento Nextflow para Genómica :material-arrow-right:](nf4_science/genomics/){ .md-button .md-button--secondary }

    ??? courses "**Nextflow para RNAseq** (RNAseq bulk)"

        Para investigadores que desean aprender a desarrollar sus propios pipelines de RNAseq. El curso utiliza un caso de uso de procesamiento de RNAseq bulk para demostrar cómo desarrollar un pipeline de RNAseq simple pero funcional.

        [Comenzar el entrenamiento Nextflow para RNAseq :material-arrow-right:](nf4_science/rnaseq/){ .md-button .md-button--secondary }

    ??? courses "**Nextflow para Imaging** (ómicas espaciales)"

        Para investigadores en imaging y ómicas espaciales que desean aprender a ejecutar y personalizar pipelines de análisis. El curso utiliza el pipeline nf-core/molkart para proporcionar un pipeline biológicamente relevante y demostrar cómo ejecutar, configurar y gestionar entradas para flujos de trabajo de Nextflow.

        [Comenzar el entrenamiento Nextflow para Imaging :material-arrow-right:](nf4_science/imaging/){ .md-button .md-button--secondary }

-   :material-run:{ .lg .middle } __Pista avanzada__

    ---

    ### :material-bridge:{.nextflow-primary} De Nextflow a nf-core {.mt-1}

    Aprenda a utilizar código y mejores prácticas del proyecto de comunidad [nf-core](https://nf-co.re/).

    Estos cursos le ayudan a pasar de los fundamentos de Nextflow a las mejores prácticas de nf-core.
    Entienda cómo y por qué la comunidad nf-core construye pipelines, y cómo puede contribuir y reutilizar estas técnicas.

    ??? courses "**Hello nf-core:** Comience con nf-core"

        Para desarrolladores que desean aprender a ejecutar y desarrollar pipelines compatibles con [nf-core](https://nf-co.re/). El curso cubre la estructura de los pipelines nf-core con suficiente detalle para permitir el desarrollo de pipelines simples pero completamente funcionales que siguen la plantilla de nf-core y las mejores prácticas de desarrollo, así como el uso de módulos nf-core existentes.

        [Comenzar el entrenamiento Hello nf-core :material-arrow-right:](hello_nf-core/index.md){ .md-button .md-button--secondary }

    ---

    ### :material-rocket-launch:{.nextflow-primary} Entrenamiento Avanzado de Nextflow {.mt-1}

    Aprenda conceptos y mecanismos avanzados para desarrollar e implementar pipelines de Nextflow para abordar casos de uso del mundo real.

    ??? courses "**Side Quests:** Profundizaciones en temas independientes"

        Mini-cursos independientes destinados a desarrolladores de Nextflow que desean ampliar su alcance y/o profundizar sus habilidades en temas particulares. Se presentan de forma lineal pero pueden tomarse en cualquier orden (consulte las dependencias en la descripción general de cada mini-curso).

        [Explorar los Side Quests :material-arrow-right:](side_quests/){ .md-button .md-button--secondary }

    ??? courses "**Colecciones de Entrenamiento:** Rutas de aprendizaje recomendadas a través de los Side Quests"

        Las Colecciones de Entrenamiento combinan múltiples Side Quests para proporcionar una experiencia de aprendizaje integral en torno a un tema o caso de uso particular.

        [Explorar las Colecciones de Entrenamiento :material-arrow-right:](training_collections/){ .md-button .md-button--secondary }

</div>

!!! info "¿Busca materiales de entrenamiento archivados?"

    Los materiales de entrenamiento más antiguos (Entrenamiento Fundamental, Entrenamiento Avanzado y otros cursos experimentales) han sido eliminados del portal de entrenamiento ya que son incompatibles con la sintaxis estricta de Nextflow 3.0.
    Si necesita acceso a estos materiales, están disponibles en el [historial de git](https://github.com/nextflow-io/training) antes de enero de 2026.

---

<div markdown class="homepage_logos">

![Seqera](assets/img/seqera_logo.png#only-light)

![Seqera](assets/img/seqera_logo_dark.png#only-dark)

</div>
