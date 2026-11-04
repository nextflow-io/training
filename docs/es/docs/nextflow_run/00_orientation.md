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

Este curso requiere Nextflow 25.10.2 o posterior, con el analizador de sintaxis v2 habilitado (el valor predeterminado en 25.10+).
Si está usando un entorno local o personalizado, asegúrese de usar la configuración correcta como se documenta [aquí](../info/nxf_versions.md).

## Prepararse para trabajar

Una vez que su codespace esté ejecutándose, hay dos cosas que necesita hacer antes de sumergirse en el entrenamiento: establecer su directorio de trabajo y echar un vistazo a los materiales proporcionados.

### Establecer el directorio de trabajo

Por defecto, el codespace se abre en la raíz de todos los cursos de entrenamiento.
Para este curso, cambie al directorio `nextflow-run/`:

```bash
cd nextflow-run/
```

Luego configure VSCode para enfocarse en este directorio, de modo que solo los archivos relevantes se muestren en la barra lateral del explorador de archivos:

```bash
code .
```

!!! tip "Consejo"

    Si por cualquier razón sale de este directorio (por ejemplo, si su codespace entra en suspensión), siempre puede usar la ruta completa para volver a él, asumiendo que está ejecutando esto dentro del entorno de entrenamiento de GitHub Codespaces:

    ```bash
    cd /workspaces/training/nextflow-run
    ```

### Explorar los materiales proporcionados

Puede explorar los materiales del curso usando el explorador de archivos a la izquierda, o con el comando `tree`.
Ejecute lo siguiente desde la terminal para ver la estructura completa:

```bash
tree . -L 2
```

??? abstract "Contenido del directorio"

    ```console
    .
    ├── 1-hello.nf
    ├── 2-inputs.nf
    ├── data
    │   ├── greetings-extended.csv
    │   └── greetings.csv
    ├── main.nf
    ├── modules
    │   ├── collectGreetings.nf
    │   ├── convertToUpper.nf
    │   ├── cowpy.nf
    │   └── sayHello.nf
    └── nextflow.config
    ```

Los **archivos `.nf`** son scripts de workflow de complejidad creciente, usados en ese orden a lo largo del curso.

El **directorio `data/`** contiene los archivos CSV de entrada que usaremos a partir de la sección 2.

El **directorio `modules/`** contiene las definiciones de procesos usadas por `main.nf`.

El **archivo `nextflow.config`** es un archivo de configuración que establece propiedades mínimas del entorno. Puede ignorarlo por ahora; lo revisaremos en la sección 4.

## Lista de verificación de preparación

¿Cree que está listo para comenzar?

- [ ] Entiendo el objetivo de este curso y sus requisitos previos
- [ ] Mi entorno está funcionando
- [ ] He establecido mi directorio de trabajo apropiadamente

Si puede marcar todas las casillas, está listo para comenzar.

**Para continuar a [Parte 1: Ejecutar Nextflow](./01_run_nextflow.md), haga clic en la flecha en la esquina inferior derecha de esta página.**
