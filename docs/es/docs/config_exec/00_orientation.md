# Primeros pasos

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Traducción asistida por IA - [más información y sugerencias](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>


## Iniciar un entorno de capacitación

Para usar el entorno preconfigurado que proporcionamos en GitHub Codespaces, haga clic en el botón "Open in GitHub Codespaces" que aparece a continuación. Para otras opciones, consulte [Opciones de entorno](../envsetup/index.md).

Recomendamos abrir el entorno de capacitación en una nueva pestaña o ventana del navegador (use clic derecho, ctrl+clic o cmd+clic según su equipo) para que pueda seguir leyendo mientras el entorno se carga.
Deberá mantener estas instrucciones abiertas en paralelo para trabajar a lo largo del curso.

[![Open in GitHub Codespaces](https://github.com/codespaces/badge.svg)](https://codespaces.new/nextflow-io/training?quickstart=1&ref=master)

### Conceptos básicos del entorno

Este entorno de capacitación contiene todo el software, el código y los datos necesarios para trabajar a lo largo del curso, por lo que no necesita instalar nada por su cuenta.

El codespace está configurado con una interfaz de VSCode, que incluye un explorador de archivos, un editor de código y una terminal.
Todas las instrucciones dadas durante el curso (por ejemplo, "abra el archivo", "edite el código" o "ejecute este comando") hacen referencia a esas tres partes de la interfaz de VSCode, a menos que se indique lo contrario.

Si está trabajando en este curso por su cuenta, familiarícese con los [conceptos básicos del entorno](../envsetup/01_setup.md) para obtener más detalles.

### Requisitos de versión

Este curso requiere Nextflow 25.10.2 o posterior, con el analizador de sintaxis v2 habilitado (el predeterminado en 25.10+).
Si está usando un entorno local o personalizado, asegúrese de estar usando la configuración correcta según se documenta [aquí](../info/nxf_versions.md).

## Prepárese para trabajar

Una vez que su codespace esté en ejecución, hay dos cosas que hacer antes de comenzar: establecer su directorio de trabajo y revisar los materiales proporcionados.

### Establecer el directorio de trabajo

De forma predeterminada, el codespace se abre en la raíz de todos los cursos de capacitación.
Para este curso, cambie al directorio `config-exec/`:

```bash
cd config-exec/
```

Luego configure VSCode para que se enfoque en este directorio, de modo que solo los archivos relevantes aparezcan en la barra lateral del explorador de archivos:

```bash
code .
```

!!! tip "Consejo"

    Si por alguna razón sale de este directorio (por ejemplo, si su codespace entra en suspensión), siempre puede usar la ruta completa para regresar, asumiendo que está ejecutando esto dentro del entorno de capacitación de GitHub Codespaces:

    ```bash
    cd /workspaces/training/config-exec
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
    ├── data
    │   └── greetings.csv
    ├── main.nf
    ├── modules
    │   ├── collectGreetings.nf
    │   ├── convertToUpper.nf
    │   ├── cowpy.nf
    │   └── sayHello.nf
    └── nextflow.config
    ```

Los archivos **`main.nf`** y **`modules/`** son el mismo pipeline de múltiples pasos de [Nextflow Run](../nextflow_run/index.md), y el archivo **`nextflow.config`** es la misma configuración que ya vio allí.
Extenderá ambos a lo largo de estos ejercicios.

El directorio **`data/`** contiene el archivo de entrada CSV que el pipeline lee.

## Lista de verificación de preparación

¿Cree que está listo para comenzar?

- [ ] Entiendo el objetivo de este curso y sus requisitos previos
- [ ] Mi entorno está en funcionamiento
- [ ] He establecido mi directorio de trabajo correctamente

Si puede marcar todas las casillas, está listo para continuar.

**Para continuar a [Parte 1: Adaptarse a su entorno de cómputo](./01_packaging_and_execution.md), haga clic en la flecha en la esquina inferior derecha de esta página.**
