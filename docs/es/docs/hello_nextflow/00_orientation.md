# Comenzando

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Traducción asistida por IA - [más información y sugerencias](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

<!--
<div class="video-wrapper">
  <iframe width="560" height="315" src="https://www.youtube.com/embed/G3CV-FcV-rc?si=nyLvwhrSB2m1NPc5&amp;list=PLPZ8WHdZGxmXiHf8B26oB_fTfoKQdhlik" title="YouTube video player" frameborder="0" allow="accelerometer; autoplay; clipboard-write; encrypted-media; gyroscope; picture-in-picture; web-share" referrerpolicy="strict-origin-when-cross-origin" allowfullscreen></iframe>
</div>

/// caption
:fontawesome-brands-youtube:{ .youtube } Vea [la lista de reproducción completa](https://www.youtube.com/playlist?list=PLPZ8WHdZGxmXiHf8B26oB_fTfoKQdhlik) en el canal de YouTube de Nextflow.

:green_book: La transcripción del video está disponible [aquí](./transcripts/00_orientation.md).
///

!!! tip "Consejo"

    ¡Los videos de YouTube tienen algunas funciones especiales!

    - :fontawesome-solid-closed-captioning: Subtítulos de alta calidad (curados manualmente). Actívelos con el ícono :material-subtitles:
    - :material-bookmark: Capítulos de video en la línea de tiempo que corresponden a los encabezados de la página.

-->

## Iniciar un entorno de entrenamiento

Para usar el entorno preconstruido que proporcionamos en GitHub Codespaces, haga clic en el botón "Open in GitHub Codespaces" a continuación. Para otras opciones, consulte [Opciones de entorno](../envsetup/index.md).

Recomendamos abrir el entorno de entrenamiento en una nueva pestaña o ventana del navegador (use clic derecho, ctrl-clic o cmd-clic dependiendo de su equipo) para que pueda seguir leyendo mientras el entorno carga.
Necesitará mantener estas instrucciones abiertas en paralelo para trabajar en el curso.

[![Open in GitHub Codespaces](https://github.com/codespaces/badge.svg)](https://codespaces.new/nextflow-io/training?quickstart=1&ref=master)

### Conceptos básicos del entorno

Este entorno de entrenamiento contiene todo el software, código y datos necesarios para trabajar en el curso de entrenamiento, por lo que no necesita instalar nada usted mismo.

El codespace está configurado con una interfaz de VSCode, que incluye un explorador de sistema de archivos, un editor de código y una terminal.
Todas las instrucciones dadas durante el curso (por ejemplo, 'abra el archivo', 'edite el código' o 'ejecute este comando') se refieren a esas tres partes de la interfaz de VSCode a menos que se especifique lo contrario.

Si está trabajando en este curso por su cuenta, familiarícese con los [conceptos básicos del entorno](../envsetup/01_setup.md) para más detalles.

### Requisitos de versión

Este entrenamiento está diseñado para Nextflow 25.10.2 o posterior **con el analizador de sintaxis v2 HABILITADO**.
Si está usando un entorno local o personalizado, asegúrese de estar usando la configuración correcta como se documenta [aquí](../info/nxf_versions.md).

## Prepárese para trabajar

Una vez que su codespace esté ejecutándose, hay dos cosas que debe hacer antes de sumergirse en el entrenamiento: establecer su directorio de trabajo para este curso específico, y echar un vistazo a los materiales proporcionados.

### Establecer el directorio de trabajo

Por defecto, el codespace se abre con el directorio de trabajo establecido en la raíz de todos los cursos de entrenamiento, pero para este curso, trabajaremos en el directorio `hello-nextflow/`.

Cambie de directorio ahora ejecutando este comando en el terminal:

```bash
cd hello-nextflow/
```

Puede configurar VSCode para enfocarse en este directorio, de modo que solo los archivos relevantes se muestren en la barra lateral del explorador de archivos:

```bash
code .
```

!!! tip "Consejo"

    Si por cualquier razón sale de este directorio (por ejemplo, su codespace se suspende), siempre puede usar la ruta completa para volver a él, asumiendo que está ejecutando esto dentro del entorno de entrenamiento de Github Codespaces:

    ```bash
    cd /workspaces/training/hello-nextflow
    ```

Ahora echemos un vistazo a los contenidos.

### Explorar los materiales proporcionados

Puede explorar los contenidos de este directorio usando el explorador de archivos en el lado izquierdo del espacio de trabajo de entrenamiento.
Alternativamente, puede usar el comando `tree`.

A lo largo del curso, usamos la salida de `tree` para representar la estructura y contenidos del directorio de forma legible, a veces con modificaciones menores para mayor claridad.

Aquí generamos una tabla de contenidos hasta el segundo nivel:

```bash
tree . -L 2
```

??? abstract "Contenidos del directorio"

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

Haga clic en el cuadro coloreado para expandir la sección y ver sus contenidos.
Usamos secciones colapsables como esta para incluir la salida esperada de comandos de forma concisa.

- **Los archivos `.nf`** son scripts de flujo de trabajo que se nombran según la parte del curso en la que se usan.

- **El archivo `nextflow.config`** es un archivo de configuración que establece propiedades mínimas del entorno.
  Puede ignorarlo por ahora.

- **El archivo `greetings.csv`** bajo `data/` contiene datos de entrada que usaremos en la mayor parte del curso. Se describe en la Parte 2 (Channels), cuando lo introducimos por primera vez.

- **Los archivos `test-params.*`** son archivos de configuración que usaremos en la Parte 6 (Configuration). Puede ignorarlos por ahora.

- **El directorio `solutions`** contiene los scripts de flujo de trabajo completados que resultan de cada paso del curso.
  Están destinados a ser usados como referencia para verificar su trabajo y solucionar cualquier problema.

## Lista de verificación de preparación

¿Cree que está listo para sumergirse?

- [ ] Entiendo el objetivo de este curso y sus prerrequisitos
- [ ] Mi entorno está funcionando
- [ ] He establecido mi directorio de trabajo apropiadamente

Si puede marcar todas las casillas, está listo para continuar.

**Para continuar a [Parte 1: Hello World](./01_hello_world.md), haga clic en la flecha en la esquina inferior derecha de esta página.**
