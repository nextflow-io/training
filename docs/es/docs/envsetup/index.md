---
title: Opciones de entorno
description: Opciones para configurar su entorno para los entrenamientos de Nextflow
hide:
  - toc
  - footer
---

# Opciones de entorno

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Traducción asistida por IA - [más información y sugerencias](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

Nuestro objetivo es proporcionar un entorno consistente y exhaustivamente probado que permita a los estudiantes concentrarse en aprender Nextflow sin tener que dedicar tiempo y esfuerzo a gestionar software.
Para ese fin, hemos desarrollado un entorno en contenedor que contiene todo el software necesario, archivos de código y datos de ejemplo para trabajar con todos nuestros cursos.

Este entorno en contenedor puede ejecutarse directamente en Github Codespaces o localmente en VS Code con la extensión Devcontainers.

<div class="grid cards" markdown>

- :material-cloud-outline:{ .lg .middle } **Github Codespaces**

  ***

  GitHub Codespaces es un servicio basado en web que nos permite proporcionar un entorno preconstruido para el entrenamiento, con todas las herramientas y datos incluidos, respaldado por máquinas virtuales en la nube. Es accesible de forma gratuita para cualquier persona con una cuenta de Github.

  [Usar Github Codespaces:material-arrow-right:](01_setup.md){ .md-button .md-button--primary .mt-1 }

- :material-laptop:{ .lg .middle } **Devcontainers locales**

  ***

  VS Code con Devcontainers proporciona un entorno de desarrollo en contenedor ejecutado localmente con todas las herramientas de entrenamiento preconfiguradas. Ofrece el mismo entorno preconstruido que Codespaces pero ejecutándose completamente en su hardware local.

  [Usar Devcontainers localmente :material-arrow-right:](03_devcontainer.md){ .md-button .md-button--primary .mt-1 }

</div>

## Instrucciones para instalación manual

Si ninguna de las opciones anteriores se adapta a sus necesidades, puede replicar este entorno en su propio sistema local instalando las dependencias de software manualmente y clonando el repositorio de entrenamiento.

[Instalación manual :material-arrow-right:](02_local.md){ .md-button .md-button--primary .mt-1 }

---

!!! info "Obsolescencia de Gitpod"

    El entrenamiento de Nextflow solía usar [Gitpod](https://gitpod.io) hasta febrero de 2025.
    Sin embargo, los creadores de Gitpod decidieron retirar la funcionalidad gratuita en favor del sistema [Gitpod Flex](https://www.gitpod.io/blog/introducing-gitpod-flex).
    Por esa razón, cambiamos a usar GitHub Codespaces, que también ofrece un entorno de desarrollo con un solo clic sin configuración previa.

    Dependiendo de cuándo se registró en Gitpod y cuándo exactamente retiren el servicio, es posible que aún pueda iniciar el entrenamiento en su antiguo IDE en la nube, aunque no podemos garantizar un acceso confiable en el futuro:
    [Abrir en Gitpod](https://gitpod.io/#https://github.com/nextflow-io/training).
