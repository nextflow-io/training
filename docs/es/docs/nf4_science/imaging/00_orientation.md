# Orientación

Esta orientación asume que ya has abierto el entorno de capacitación haciendo clic en el botón "Open in GitHub Codespaces".
Si no lo has hecho, por favor hazlo ahora, idealmente en una segunda ventana o pestaña del navegador para que puedas consultar estas instrucciones.

[![Open in GitHub Codespaces](https://github.com/codespaces/badge.svg)](https://github.com/codespaces/new?hide_repo_select=true&ref=master&repo=290790519&skip_quickstart=true&machine=premiumLinux&devcontainer_path=.devcontainer%2Fdevcontainer.json)

!!!warning "Requisito de tamaño de máquina"

    Asegúrate de seleccionar una **máquina de 8 núcleos** al crear tu Codespace para este curso de capacitación. Los workflows de bioimagen requieren recursos computacionales adicionales.

## GitHub Codespaces

El entorno de GitHub Codespaces contiene todo el software, código y datos necesarios para trabajar en este curso de capacitación, por lo que no necesitas instalar nada por tu cuenta.
Sin embargo, sí necesitas una cuenta de GitHub (gratuita) para iniciar sesión, y si no estás familiarizado con la interfaz, deberías tomarte unos minutos para familiarizarte con ella completando el mini-curso de [Orientación de GitHub Codespaces](../../envsetup/index.md).

## Pre-descarga de imágenes Docker

Una vez que hayas abierto tu Codespace, vamos a pre-descargar todas las imágenes Docker que necesitaremos para este curso de capacitación.
Esto ahorrará tiempo más adelante y garantizará una ejecución fluida de los workflows.

Abre una nueva pestaña de terminal y ejecuta el siguiente comando:

```bash
nextflow run nf-core/molkart -profile docker,test -stub -resume --outdir results
```

Este comando descargará todas las imágenes Docker necesarias en segundo plano.
Puedes continuar con el resto de la orientación mientras esto se ejecuta.

!!!tip

    La bandera `-stub` permite que el pipeline se ejecute rápidamente sin procesar datos reales, lo cual es perfecto para descargar imágenes. Puedes monitorear el progreso en la pestaña del terminal.

## Directorio de trabajo

A lo largo de este curso de capacitación, trabajaremos en el directorio `nf4-science/imaging/`.

Cambia de directorio ahora ejecutando este comando en el terminal:

```bash
cd nf4-science/imaging/
```

!!!tip

    Si por alguna razón sales de este directorio, siempre puedes usar la ruta completa para regresar a él, asumiendo que estás ejecutando esto dentro del entorno de capacitación de GitHub Codespaces:

    ```bash
    cd /workspaces/training/nf4-science/imaging
    ```

**Ahora, para comenzar el curso, haz clic en la flecha en la esquina inferior derecha de esta página.**
