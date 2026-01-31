# Primeros pasos

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Traducción asistida por IA - [más información y sugerencias](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

## Iniciar un entorno de entrenamiento

Para usar el entorno preconfigurado que proporcionamos en GitHub Codespaces, haga clic en el botón "Open in GitHub Codespaces" a continuación. Para otras opciones, consulte [Opciones de entorno](../envsetup/index.md).

Recomendamos abrir el entorno de entrenamiento en una nueva pestaña o ventana del navegador (use clic derecho, ctrl-clic o cmd-clic dependiendo de su equipo) para que pueda seguir leyendo mientras se carga el entorno.
Necesitará mantener estas instrucciones abiertas en paralelo para trabajar a través del curso.

[![Open in GitHub Codespaces](https://github.com/codespaces/badge.svg)](https://codespaces.new/nextflow-io/training?quickstart=1&ref=master)

### Conceptos básicos del entorno

Este entorno de entrenamiento contiene todo el software, código y datos necesarios para trabajar a través del curso de entrenamiento, por lo que no necesita instalar nada usted mismo.

El codespace está configurado con una interfaz VSCode, que incluye un explorador de archivos, un editor de código y una terminal.
Todas las instrucciones dadas durante el curso (por ejemplo, 'abrir el archivo', 'editar el código' o 'ejecutar este comando') se refieren a esas tres partes de la interfaz VSCode a menos que se especifique lo contrario.

Si está trabajando en este curso por su cuenta, familiarícese con los [conceptos básicos del entorno](../envsetup/01_setup.md) para más detalles.

### Requisitos de versión

Este entrenamiento está diseñado para Nextflow 25.10.2 o posterior **con el analizador de sintaxis v2 HABILITADO**.
Si está usando un entorno local o personalizado, asegúrese de usar la configuración correcta como se documenta [aquí](../info/nxf_versions.md).

## Prepararse para trabajar

Una vez que su codespace esté ejecutándose, hay dos cosas que necesita hacer antes de sumergirse en el entrenamiento: establecer su directorio de trabajo para este curso específico y echar un vistazo a los materiales proporcionados.

### Establecer el directorio de trabajo

Por defecto, el codespace se abre con el directorio de trabajo establecido en la raíz de todos los cursos de entrenamiento, pero para este curso, trabajaremos en el directorio `nextflow-run/`.

Cambie de directorio ahora ejecutando este comando en la terminal:

```bash
cd nextflow-run/
```

Puede configurar VSCode para enfocarse en este directorio, de modo que solo los archivos relevantes se muestren en la barra lateral del explorador de archivos:

```bash
code .
```

!!! tip "Consejo"

    Si por cualquier razón sale de este directorio (por ejemplo, si su codespace entra en suspensión), siempre puede usar la ruta completa para volver a él, asumiendo que está ejecutando esto dentro del entorno de entrenamiento de GitHub Codespaces:

    ```bash
    cd /workspaces/training/nextflow-run
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
    ├── 1-hello.nf
    ├── 2a-inputs.nf
    ├── 2b-multistep.nf
    ├── 2c-modules.nf
    ├── 2d-container.nf
    ├── 3-main.nf
    ├── data
    │   └── greetings.csv
    ├── modules
    │   ├── collectGreetings.nf
    │   ├── convertToUpper.nf
    │   ├── cowpy.nf
    │   └── sayHello.nf
    ├── nextflow.config
    ├── solutions
    │   ├── 3-main.nf
    │   ├── modules
    │   └── nextflow.config
    ├── test-params.json
    └── test-params.yaml
    ```

Haga clic en el cuadro coloreado para expandir la sección y ver su contenido.
Usamos secciones colapsables como esta para mostrar la salida esperada de comandos así como contenidos de directorios y archivos de manera concisa.

- **Los archivos `.nf`** son scripts de workflow que están numerados según la parte del curso en la que se usan.

- **El archivo `nextflow.config`** es un archivo de configuración que establece propiedades mínimas del entorno.
  Puede ignorarlo por ahora.

- **El archivo `greetings.csv`** bajo `data/` contiene datos de entrada que usaremos en la mayor parte del curso. Se describe en la Parte 2 (Ejecutar pipelines), cuando lo introducimos por primera vez.

- **Los archivos `test-params.*`** son archivos de configuración que usaremos en la Parte 3 (Configuración). Puede ignorarlos por ahora.

- **El directorio `solutions`** contiene el estado final del workflow y sus archivos accesorios (config y modules) que resultan de completar el curso.
  Están destinados a ser usados como referencia para verificar su trabajo y solucionar cualquier problema.

## Lista de verificación de preparación

¿Cree que está listo para comenzar?

- [ ] Entiendo el objetivo de este curso y sus requisitos previos
- [ ] Mi entorno está funcionando
- [ ] He establecido mi directorio de trabajo apropiadamente

Si puede marcar todas las casillas, está listo para comenzar.

**Para continuar a [Parte 1: Ejecutar operaciones básicas](./01_basics.md), haga clic en la flecha en la esquina inferior derecha de esta página.**
