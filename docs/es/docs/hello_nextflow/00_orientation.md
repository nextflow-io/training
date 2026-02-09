# Primeros pasos

<div class="video-wrapper">
  <iframe width="560" height="315" src="https://www.youtube.com/embed/PIjOdFaYwWA?si=gZxlXgkVxuLEzOsC" title="YouTube video player" frameborder="0" allow="accelerometer; autoplay; clipboard-write; encrypted-media; gyroscope; picture-in-picture; web-share" referrerpolicy="strict-origin-when-cross-origin" allowfullscreen></iframe>
</div>

/// caption
:fontawesome-brands-youtube:{ .youtube } Consulte [la lista de reproducción completa](https://youtube.com/playlist?list=PLPZ8WHdZGxmWKozQuzr27jyMGqp9kElVK&si=eF7cLR62goy-lc6n) en el canal de YouTube de Nextflow.

:green_book: La transcripción del video está disponible [aquí](./transcripts/00_orientation.md).
///

!!! tip

    ¡Los videos de YouTube tienen súper poderes!

    - :fontawesome-solid-closed-captioning: Subtítulos de alta calidad (curados manualmente). Actívelos con el ícono :material-subtitles:
    - :material-bookmark: Capítulos de video en la línea de tiempo que corresponden a los encabezados de la página.

## Iniciar un entorno de capacitación

Para usar el entorno preconfigurado que proporcionamos en GitHub Codespaces, haga clic en el botón "Open in GitHub Codespaces" a continuación. Para otras opciones, consulte [Opciones de entorno](../envsetup/index.md).

Recomendamos abrir el entorno de capacitación en una nueva pestaña o ventana del navegador (use clic derecho, ctrl-clic o cmd-clic según su equipo) para que pueda continuar leyendo mientras se carga el entorno.
Necesitará mantener estas instrucciones abiertas en paralelo para trabajar en el curso.

[![Open in GitHub Codespaces](https://github.com/codespaces/badge.svg)](https://codespaces.new/nextflow-io/training?quickstart=1&ref=master)

### Conceptos básicos del entorno

Este entorno de capacitación contiene todo el software, código y datos necesarios para trabajar en el curso de capacitación, por lo que no necesita instalar nada usted mismo.

El codespace está configurado con una interfaz VSCode, que incluye un explorador de archivos, un editor de código y un terminal shell.
Todas las instrucciones dadas durante el curso (por ejemplo, 'abra el archivo', 'edite el código' o 'ejecute este comando') se refieren a esas tres partes de la interfaz VSCode a menos que se especifique lo contrario.

Si está trabajando en este curso por su cuenta, familiarícese con los [conceptos básicos del entorno](../envsetup/01_setup.md) para obtener más detalles.

### Requisitos de versión

Esta capacitación está diseñada para Nextflow 25.10.2 o posterior **con el analizador de sintaxis v2 HABILITADO**.
Si está utilizando un entorno local o personalizado, asegúrese de estar usando la configuración correcta como se documenta [aquí](../info/nxf_versions.md).

## Prepárese para trabajar

Una vez que su codespace esté en ejecución, hay dos cosas que debe hacer antes de sumergirse en la capacitación: establecer su directorio de trabajo para este curso específico y revisar los materiales proporcionados.

### Establecer el directorio de trabajo

Por defecto, el codespace se abre con el directorio de trabajo establecido en la raíz de todos los cursos de capacitación, pero para este curso, trabajaremos en el directorio `hello-nextflow/`.

Cambie de directorio ahora ejecutando este comando en el terminal:

```bash
cd hello-nextflow/
```

Puede configurar VSCode para enfocarse en este directorio, de modo que solo los archivos relevantes se muestren en la barra lateral del explorador de archivos:

```bash
code .
```

!!! tip

    Si por alguna razón sale de este directorio (por ejemplo, su codespace se suspende), siempre puede usar la ruta completa para regresar a él, asumiendo que está ejecutando esto dentro del entorno de capacitación de Github Codespaces:

    ```bash
    cd /workspaces/training/hello-nextflow
    ```

Ahora echemos un vistazo al contenido.

### Explorar los materiales proporcionados

Puede explorar el contenido de este directorio usando el explorador de archivos en el lado izquierdo del espacio de trabajo de capacitación.
Alternativamente, puede usar el comando `tree`.

A lo largo del curso, usamos la salida de `tree` para representar la estructura y el contenido del directorio en una forma legible, a veces con modificaciones menores para mayor claridad.

Aquí generamos una tabla de contenidos hasta el segundo nivel:

```bash
tree . -L 2
```

??? abstract "Contenido del directorio"

    ```console
    .
    ├── data
    │   └── greetings.csv
    ├── hello-channels.nf
    ├── hello-config.nf
    ├── hello-containers.nf
    ├── hello-modules.nf
    ├── hello-workflow.nf
    ├── hello-world.nf
    ├── nextflow.config
    ├── solutions
    │   ├── 1-hello-world
    │   ├── 2-hello-channels
    │   ├── 3-hello-workflow
    │   ├── 4-hello-modules
    │   ├── 5-hello-containers
    │   └── 6-hello-config
    ├── test-params.json
    └── test-params.yaml
    ```

Haga clic en el cuadro de color para expandir la sección y ver su contenido.
Usamos secciones plegables como esta para incluir la salida esperada del comando de manera concisa.

- **Los archivos `.nf`** son scripts de workflow que se nombran según la parte del curso en la que se usan.

- **El archivo `nextflow.config`** es un archivo de configuración que establece propiedades mínimas del entorno.
  Puede ignorarlo por ahora.

- **El archivo `greetings.csv`** bajo `data/` contiene datos de entrada que usaremos en la mayor parte del curso. Se describe en la Parte 2 (Canales), cuando lo introducimos por primera vez.

- **Los archivos `test-params.*`** son archivos de configuración que usaremos en la Parte 6 (Configuración). Puede ignorarlos por ahora.

- **El directorio `solutions`** contiene los scripts de workflow completados que resultan de cada paso del curso.
  Están destinados a ser utilizados como referencia para verificar su trabajo y solucionar cualquier problema.

## Lista de verificación de preparación

¿Cree que está listo para comenzar?

- [ ] Entiendo el objetivo de este curso y sus requisitos previos
- [ ] Mi entorno está en funcionamiento
- [ ] He establecido mi directorio de trabajo apropiadamente

Si puede marcar todas las casillas, está listo para comenzar.

**Para continuar a [Parte 1: Hello World](./01_hello_world.md), haga clic en la flecha en la esquina inferior derecha de esta página.**
