# Primeros pasos

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Traducción asistida por IA - [más información y sugerencias](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

## Iniciar un entorno de capacitación

Para usar el entorno preconfigurado que proporcionamos en GitHub Codespaces, haga clic en el botón "Open in GitHub Codespaces" a continuación. Para otras opciones, consulte [Opciones de entorno](../../envsetup/index.md).

Recomendamos abrir el entorno de capacitación en una nueva pestaña o ventana del navegador (use clic derecho, ctrl-clic o cmd-clic según su equipo) para que pueda leer mientras se carga el entorno.
Necesitará mantener estas instrucciones abiertas en paralelo para trabajar en el curso.

[![Open in GitHub Codespaces](https://github.com/codespaces/badge.svg)](https://codespaces.new/nextflow-io/training?quickstart=1&ref=master)

### Conceptos básicos del entorno

Este entorno de capacitación contiene todo el software, código y datos necesarios para trabajar en el curso de capacitación, por lo que no necesita instalar nada por su cuenta.

El codespace está configurado con una interfaz VSCode, que incluye un explorador de archivos, un editor de código y un terminal shell.
Todas las instrucciones dadas durante el curso (por ejemplo, 'abrir el archivo', 'editar el código' o 'ejecutar este comando') se refieren a esas tres partes de la interfaz de VScode a menos que se especifique lo contrario.

Si está trabajando en este curso por su cuenta, por favor familiarícese con los [conceptos básicos del entorno](../../envsetup/01_setup.md) para más detalles.

### Requisitos de versión

Esta capacitación está diseñada para Nextflow 25.10.2 o posterior **con el analizador de sintaxis v2 HABILITADO**.
Si está usando un entorno local o personalizado, por favor asegúrese de estar usando la configuración correcta como está documentado [aquí](../../info/nxf_versions.md).

## Prepárese para trabajar

Una vez que su codespace esté en ejecución, hay dos cosas que necesita hacer antes de sumergirse en la capacitación: establecer su directorio de trabajo para este curso específico y revisar los materiales proporcionados.

### Establecer el directorio de trabajo

Por defecto, el codespace se abre con el directorio de trabajo establecido en la raíz de todos los cursos de capacitación, pero para este curso, trabajaremos en el directorio `nf4-science/genomics/`.

Cambie de directorio ahora ejecutando este comando en el terminal:

```bash
cd nf4-science/genomics/
```

Puede configurar VSCode para enfocarse en este directorio, de modo que solo los archivos relevantes se muestren en la barra lateral del explorador de archivos:

```bash
code .
```

!!! tip "Consejo"

    Si por cualquier razón sale de este directorio (por ejemplo, su codespace se suspende), siempre puede usar la ruta completa para regresar a él, asumiendo que está ejecutando esto dentro del entorno de capacitación de Github Codespaces:

    ```bash
    cd /workspaces/training/nf4-science/genomics
    ```

Ahora echemos un vistazo a los contenidos.

### Explorar los materiales proporcionados

Puede explorar los contenidos de este directorio usando el explorador de archivos en el lado izquierdo del espacio de trabajo de capacitación.
Alternativamente, puede usar el comando `tree`.

A lo largo del curso, usamos la salida de `tree` para representar la estructura y contenido del directorio de una forma legible, a veces con modificaciones menores para mayor claridad.

Aquí generamos una tabla de contenidos hasta el segundo nivel:

```bash
tree . -L 2
```

??? abstract "Contenido del directorio"

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

Haga clic en el cuadro de color para expandir la sección y ver su contenido.
Usamos secciones plegables como esta para mostrar la salida esperada de comandos, así como el contenido de directorios y archivos de forma concisa.

- **El archivo `genomics.nf`** es un script de workflow que construirá a lo largo del curso.

- **El directorio `modules`** contiene archivos de módulo esqueleto que completará durante el curso.

- **El archivo `nextflow.config`** es un archivo de configuración que establece propiedades mínimas del entorno.
  Puede ignorarlo por ahora.

- **El directorio `data`** contiene datos de entrada y recursos relacionados, descritos más adelante en el curso.

- **El directorio `solutions`** contiene archivos de módulo completados y una solución de la Parte 2 que puede servir como punto de partida para la Parte 3.
  Están destinados a ser usados como referencia para verificar su trabajo y solucionar cualquier problema.

## Lista de verificación de preparación

¿Cree que está listo para comenzar?

- [ ] Entiendo el objetivo de este curso y sus requisitos previos
- [ ] Mi entorno está en funcionamiento
- [ ] He establecido mi directorio de trabajo apropiadamente

Si puede marcar todas las casillas, está listo para comenzar.

**Para continuar a la [Parte 1: Descripción general del método y pruebas manuales](./01_method.md), haga clic en la flecha en la esquina inferior derecha de esta página.**
