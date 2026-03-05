# Instalación manual

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Traducción asistida por IA - [más información y sugerencias](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

Es posible instalar todo lo que necesita para ejecutar el entrenamiento en su propio entorno local manualmente.

Aquí hemos documentado cómo hacerlo en sistemas estándar compatibles con POSIX (asumiendo una máquina personal como una laptop).
Tenga en cuenta que algunos detalles pueden ser diferentes dependiendo de su sistema específico.

!!! tip "Consejo"

    Antes de continuar, ¿ha considerado usar el [enfoque de Devcontainers](03_devcontainer.md)?
    Proporciona todas las herramientas y dependencias necesarias sin requerir instalación manual.

## Requisitos generales de software

Nextflow se puede usar en cualquier sistema compatible con POSIX (Linux, macOS, Windows Subsystem for Linux, etc.) con Java instalado.
Nuestros cursos de entrenamiento tienen algunos requisitos adicionales.

En total, necesitará tener instalado el siguiente software:

- Bash o shell equivalente
- [Java 11 (o posterior, hasta 21)](https://www.oracle.com/technetwork/java/javase/downloads/index.html)
- [Git](https://git-scm.com/)
- [Docker](https://docs.docker.com/get-docker/)
- [Conda](https://conda.io/) 4.5 (o posterior)
- [VSCode](https://code.visualstudio.com) con la [extensión de Nextflow](https://www.nextflow.io/docs/latest/developer-env.html#devenv-nextflow)

La aplicación VSCode es técnicamente opcional pero recomendamos encarecidamente que la use para trabajar en los cursos así como para su trabajo de desarrollo de Nextflow en general.

El manual de documentación de Nextflow proporciona instrucciones para instalar estas dependencias en [Configuración del entorno](https://www.nextflow.io/docs/latest/developer-env.html).

## Nextflow y herramientas de nf-core

Necesitará instalar Nextflow en sí, además de las herramientas de nf-core, como se detalla en los artículos enlazados a continuación:

- [Instalación de Nextflow](https://www.nextflow.io/docs/latest/install.html)
- [Herramientas de nf-core](https://nf-co.re/docs/nf-core-tools/installation)

Recomendamos usar la opción de autoinstalación para Nextflow y la opción de PyPI para las herramientas de nf-core.

!!! warning "Compatibilidad de versiones"

    <!-- Cualquier actualización de este contenido debe copiarse a la página principal -->
    **A partir de enero de 2026, todos nuestros cursos de entrenamiento de Nextflow requieren la versión 25.10.2 o posterior de Nextflow, con la sintaxis estricta v2 activada, a menos que se indique lo contrario.**

    Para más información sobre los requisitos de versión y la sintaxis estricta v2, consulte la guía de [versiones de Nextflow](../info/nxf_versions.md).

    Las versiones anteriores del material de entrenamiento correspondientes a la sintaxis anterior están disponibles a través del selector de versiones en la barra de menú de esta página web.

## Materiales de entrenamiento

La forma más fácil de descargar los materiales de entrenamiento es clonar todo el repositorio usando este comando:

```bash
git clone https://github.com/nextflow-io/training.git
```

Cada curso tiene su propio directorio.
Para trabajar en un curso, abra una ventana de terminal (idealmente, desde dentro de la aplicación VSCode) y use `cd` para entrar en el directorio relevante.

Luego puede seguir las instrucciones del curso proporcionadas en el sitio web.
