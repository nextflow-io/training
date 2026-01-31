# Orientación

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Traducción asistida por IA - [más información y sugerencias](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

Esta orientación asume que ya ha abierto el entorno de entrenamiento haciendo clic en el botón "Open in GitHub Codespaces".
Si no lo ha hecho, hágalo ahora, idealmente en una segunda ventana o pestaña del navegador para poder consultar estas instrucciones.

[![Open in GitHub Codespaces](https://github.com/codespaces/badge.svg)](https://github.com/codespaces/new?hide_repo_select=true&ref=master&repo=290790519&skip_quickstart=true&machine=premiumLinux&devcontainer_path=.devcontainer%2Fdevcontainer.json)

!!!warning "Requisito de tamaño de máquina"

    Asegúrese de seleccionar una **máquina de 8 núcleos** al crear su Codespace para este curso de entrenamiento. Los flujos de trabajo de bioimagen requieren recursos computacionales adicionales.

## GitHub Codespaces

El entorno de GitHub Codespaces contiene todo el software, código y datos necesarios para trabajar en este curso de entrenamiento, por lo que no necesita instalar nada usted mismo.
Sin embargo, sí necesita una cuenta (gratuita) de GitHub para iniciar sesión, y si no está familiarizado con la interfaz, debería tomarse unos minutos para familiarizarse con ella completando el mini-curso [GitHub Codespaces Orientation](../../envsetup/index.md).

## Pre-descargar imágenes de Docker

Una vez que haya abierto su Codespace, pre-descarguemos todas las imágenes de Docker que necesitaremos para este curso de entrenamiento.
Esto ahorrará tiempo más adelante y asegurará una ejecución fluida de los flujos de trabajo.

Abra una nueva pestaña de terminal y ejecute el siguiente comando:

```bash
nextflow run nf-core/molkart -profile docker,test -stub -resume --outdir results
```

Este comando descargará todas las imágenes de Docker necesarias en segundo plano.
Puede continuar con el resto de la orientación mientras esto se ejecuta.

!!!tip "Consejo"

    La bandera `-stub` permite que el pipeline se ejecute rápidamente sin procesar datos reales, lo cual es perfecto para descargar imágenes. Puede monitorear el progreso en la pestaña de terminal.

## Directorio de trabajo

A lo largo de este curso de entrenamiento, trabajaremos en el directorio `nf4-science/imaging/`.

Cambie de directorio ahora ejecutando este comando en la terminal:

```bash
cd nf4-science/imaging/
```

!!!tip "Consejo"

    Si por alguna razón se mueve fuera de este directorio, siempre puede usar la ruta completa para regresar a él, asumiendo que está ejecutando esto dentro del entorno de entrenamiento de GitHub Codespaces:

    ```bash
    cd /workspaces/training/nf4-science/imaging
    ```

**Ahora, para comenzar el curso, haga clic en la flecha en la esquina inferior derecha de esta página.**
