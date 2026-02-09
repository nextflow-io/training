# Instalación manual

Es posible instalar todo lo que necesitas para ejecutar la capacitación en tu propio entorno local de forma manual.

Aquí hemos documentado cómo hacerlo en sistemas estándar compatibles con POSIX (asumiendo una máquina personal como una laptop).
Ten en cuenta que algunos detalles pueden ser diferentes dependiendo de tu sistema específico.

!!! tip "Consejo"

    Antes de continuar, ¿has considerado usar el [enfoque de Devcontainers](03_devcontainer.md)?
    Proporciona todas las herramientas y dependencias necesarias sin requerir instalación manual.

## Requisitos generales de software

Nextflow puede usarse en cualquier sistema compatible con POSIX (Linux, macOS, Windows Subsystem for Linux, etc.) con Java instalado.
Nuestros cursos de capacitación tienen algunos requisitos adicionales.

En total, necesitarás tener el siguiente software instalado:

- Bash o shell equivalente
- [Java 11 (o posterior, hasta 21)](https://www.oracle.com/technetwork/java/javase/downloads/index.html)
- [Git](https://git-scm.com/)
- [Docker](https://docs.docker.com/get-docker/)
- [Conda](https://conda.io/) 4.5 (o posterior)
- [VSCode](https://code.visualstudio.com) con la [extensión de Nextflow](https://www.nextflow.io/docs/latest/developer-env.html#devenv-nextflow)

La aplicación VSCode es técnicamente opcional, pero recomendamos encarecidamente que la uses para trabajar en los cursos, así como para tu trabajo de desarrollo con Nextflow en general.

El manual de documentación de Nextflow proporciona instrucciones para instalar estas dependencias en [Configuración del entorno](https://www.nextflow.io/docs/latest/developer-env.html).

## Nextflow y herramientas nf-core

Necesitarás instalar Nextflow mismo, además de las herramientas nf-core, como se detalla en los artículos vinculados a continuación:

- [Instalación de Nextflow](https://www.nextflow.io/docs/latest/install.html)
- [Herramientas nf-core](https://nf-co.re/docs/nf-core-tools/installation)

Recomendamos usar la opción de autoinstalación para Nextflow y la opción PyPI para las herramientas nf-core.

!!! warning "Advertencia sobre compatibilidad de versiones"

    <!-- Any update to this content needs to be copied to the home page -->
    **A partir de enero de 2026, todos nuestros cursos de capacitación de Nextflow requieren Nextflow versión 25.10.2 o posterior, con sintaxis estricta v2 activada, a menos que se indique lo contrario.**

    Para más información sobre los requisitos de versión y la sintaxis estricta v2, consulte la guía de [versiones de Nextflow](../info/nxf_versions.md).

    Las versiones anteriores del material de capacitación correspondientes a sintaxis previas están disponibles a través del selector de versión en la barra de menú de esta página web.

## Materiales de capacitación

La forma más fácil de descargar los materiales de capacitación es clonar todo el repositorio usando este comando:

```bash
git clone https://github.com/nextflow-io/training.git
```

Cada curso tiene su propio directorio.
Para trabajar en un curso, abre una ventana de terminal (idealmente, desde dentro de la aplicación VSCode) y usa `cd` para entrar al directorio relevante.

Luego puedes seguir las instrucciones del curso proporcionadas en el sitio web.
