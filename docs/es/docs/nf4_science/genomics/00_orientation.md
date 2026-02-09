# Primeros pasos

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Traducción asistida por IA - [más información y sugerencias](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

## Iniciar un entorno de capacitación

Para usar el entorno preconfigurado que proporcionamos en GitHub Codespaces, haz clic en el botón "Open in GitHub Codespaces" a continuación. Para otras opciones, consulta [Opciones de entorno](../../envsetup/index.md).

Recomendamos abrir el entorno de capacitación en una nueva pestaña o ventana del navegador (usa clic derecho, ctrl-clic o cmd-clic según tu equipo) para que puedas seguir leyendo mientras se carga el entorno.
Necesitarás mantener estas instrucciones abiertas en paralelo para trabajar en el curso.

[![Open in GitHub Codespaces](https://github.com/codespaces/badge.svg)](https://codespaces.new/nextflow-io/training?quickstart=1&ref=master)

### Conceptos básicos del entorno

Este entorno de capacitación contiene todo el software, código y datos necesarios para trabajar en el curso de capacitación, por lo que no necesitas instalar nada por tu cuenta.

El codespace está configurado con una interfaz VSCode, que incluye un explorador de archivos, un editor de código y una terminal shell.
Todas las instrucciones dadas durante el curso (por ejemplo, 'abre el archivo', 'edita el código' o 'ejecuta este comando') se refieren a esas tres partes de la interfaz VSCode a menos que se especifique lo contrario.

Si estás trabajando en este curso por tu cuenta, familiarízate con los [conceptos básicos del entorno](../../envsetup/01_setup.md) para obtener más detalles.

### Requisitos de versión

Esta capacitación está diseñada para Nextflow 25.10.2 o posterior **con el analizador de sintaxis v2 HABILITADO**.
Si estás usando un entorno local o personalizado, asegúrate de estar usando la configuración correcta como se documenta [aquí](../../info/nxf_versions.md).

## Prepárate para trabajar

Una vez que tu codespace esté en ejecución, hay dos cosas que debes hacer antes de sumergirte en la capacitación: establecer tu directorio de trabajo para este curso específico y revisar los materiales proporcionados.

### Establecer el directorio de trabajo

Por defecto, el codespace se abre con el directorio de trabajo establecido en la raíz de todos los cursos de capacitación, pero para este curso, trabajaremos en el directorio `nf4-science/genomics/`.

Cambia de directorio ahora ejecutando este comando en la terminal:

```bash
cd nf4-science/genomics/
```

Puedes configurar VSCode para enfocarse en este directorio, de modo que solo los archivos relevantes se muestren en la barra lateral del explorador de archivos:

```bash
code .
```

!!! tip "Consejo"

    Si por alguna razón sales de este directorio (por ejemplo, tu codespace se suspende), siempre puedes usar la ruta completa para regresar a él, asumiendo que estás ejecutando esto dentro del entorno de capacitación de Github Codespaces:

    ```bash
    cd /workspaces/training/nf4-science/genomics
    ```

Ahora echemos un vistazo a los contenidos.

### Explorar los materiales proporcionados

Puedes explorar los contenidos de este directorio usando el explorador de archivos en el lado izquierdo del espacio de trabajo de capacitación.
Alternativamente, puedes usar el comando `tree`.

A lo largo del curso, usamos la salida de `tree` para representar la estructura y contenidos del directorio de forma legible, a veces con modificaciones menores para mayor claridad.

Aquí generamos una tabla de contenidos hasta el segundo nivel:

```bash
tree . -L 2
```

??? abstract "Contenidos del directorio"

    ```console
    .
    ├── data
    │   ├── bam
    │   ├── ref
    │   ├── sample_bams.txt
    │   └── samplesheet.csv
    ├── genomics.nf
    ├── modules
    │   ├── gatk_haplotypecaller.nf
    │   └── samtools_index.nf
    ├── nextflow.config
    └── solutions
        ├── modules
        ├── nf-test.config
        ├── part2
        └── tests

    8 directories, 8 files
    ```

Haz clic en el cuadro de color para expandir la sección y ver sus contenidos.
Usamos secciones plegables como esta para mostrar la salida esperada de comandos, así como los contenidos de directorios y archivos de manera concisa.

- **El archivo `genomics.nf`** es un script de workflow que construirás a lo largo del curso.

- **El directorio `modules`** contiene archivos de módulos esqueleto que completarás durante el curso.

- **El archivo `nextflow.config`** es un archivo de configuración que establece propiedades mínimas del entorno.
  Puedes ignorarlo por ahora.

- **El directorio `data`** contiene datos de entrada y recursos relacionados, descritos más adelante en el curso.

- **El directorio `solutions`** contiene archivos de módulos completados y una solución de la Parte 2 que puede servir como punto de partida para la Parte 3.
  Están destinados a ser usados como referencia para verificar tu trabajo y solucionar cualquier problema.

## Lista de verificación de preparación

¿Crees que estás listo para comenzar?

- [ ] Entiendo el objetivo de este curso y sus requisitos previos
- [ ] Mi entorno está funcionando
- [ ] He establecido mi directorio de trabajo apropiadamente

Si puedes marcar todas las casillas, estás listo para comenzar.

**Para continuar a la [Parte 1: Descripción general del método y pruebas manuales](./01_method.md), haz clic en la flecha en la esquina inferior derecha de esta página.**
