# Opciones de entorno

Nuestro objetivo es proporcionar un entorno consistente y completamente probado que permita a los estudiantes enfocarse en aprender Nextflow sin tener que dedicar tiempo y esfuerzo a la gestión de software.
Con ese fin, hemos desarrollado un entorno en contenedor que contiene todo el software necesario, archivos de código y datos de ejemplo para trabajar en todos nuestros cursos.

Este entorno en contenedor puede ejecutarse directamente en Github Codespaces o localmente en VS Code con la extensión Devcontainers.

<div class="grid cards" markdown>

-   :material-cloud-outline:{ .lg .middle } __Github Codespaces__

    ---

    GitHub Codespaces es un servicio basado en la web que nos permite proporcionar un entorno preconfigurado para capacitación, con todas las herramientas y datos incluidos, respaldado por máquinas virtuales en la nube. Es accesible de forma gratuita para cualquier persona con una cuenta de Github.

    [Usar Github Codespaces:material-arrow-right:](01_setup.md){ .md-button .md-button--primary .mt-1 }

-   :material-laptop:{ .lg .middle } __Devcontainers locales__

    ---

    VS Code con Devcontainers proporciona un entorno de desarrollo en contenedor que se ejecuta localmente con todas las herramientas de capacitación preconfiguradas. Ofrece el mismo entorno preconfigurado que Codespaces pero ejecutándose completamente en tu hardware local.

    [Usar Devcontainers localmente :material-arrow-right:](03_devcontainer.md){ .md-button .md-button--primary .mt-1 }

</div>

## Instrucciones para instalación manual

Si ninguna de las opciones anteriores se ajusta a tus necesidades, puedes replicar este entorno en tu propio sistema local instalando las dependencias de software manualmente y clonando el repositorio de capacitación.

[Instalación manual :material-arrow-right:](02_local.md){ .md-button .md-button--primary .mt-1 }

---

!!! info "Discontinuación de Gitpod"

    Nextflow Training solía usar [Gitpod](https://gitpod.io) hasta febrero de 2025.
    Sin embargo, los creadores de Gitpod decidieron retirar la funcionalidad gratuita en favor del sistema [Gitpod Flex](https://www.gitpod.io/blog/introducing-gitpod-flex).
    Por esa razón, cambiamos a usar GitHub Codespaces, que también ofrece un entorno de desarrollo con un solo clic sin configuración previa.

    Dependiendo de cuándo te registraste en Gitpod y cuándo exactamente retiren el servicio, es posible que aún puedas iniciar la capacitación en su antiguo IDE en la nube, aunque no podemos garantizar un acceso confiable en el futuro:
    [Abrir en Gitpod](https://gitpod.io/#https://github.com/nextflow-io/training).
