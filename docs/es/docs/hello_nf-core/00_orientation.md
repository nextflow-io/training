# Primeros pasos

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Traducción asistida por IA - [más información y sugerencias](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

## Iniciar un entorno de entrenamiento

Para usar el entorno preconfigurado que proporcionamos en GitHub Codespaces, haga clic en el botón "Open in GitHub Codespaces" a continuación. Para otras opciones, consulte [Opciones de entorno](../envsetup/index.md).

Recomendamos abrir el entorno de entrenamiento en una nueva pestaña o ventana del navegador (use clic derecho, ctrl-clic o cmd-clic según su equipo) para que pueda continuar leyendo mientras el entorno se carga.
Necesitará mantener estas instrucciones abiertas en paralelo para trabajar a través del curso.

[![Open in GitHub Codespaces](https://github.com/codespaces/badge.svg)](https://codespaces.new/nextflow-io/training?quickstart=1&ref=master)

### Conceptos básicos del entorno

Este entorno de entrenamiento contiene todo el software, código y datos necesarios para trabajar a través del curso de entrenamiento, por lo que no necesita instalar nada por su cuenta.

El codespace está configurado con una interfaz VSCode, que incluye un explorador de sistema de archivos, un editor de código y un terminal shell.
Todas las instrucciones dadas durante el curso (por ejemplo, 'abrir el archivo', 'editar el código' o 'ejecutar este comando') se refieren a estas tres partes de la interfaz VSCode a menos que se especifique lo contrario.

Si está trabajando en este curso por su cuenta, por favor familiarícese con los [conceptos básicos del entorno](../envsetup/01_setup.md) para más detalles.

### Requisitos de versión

Este entrenamiento funciona con **Nextflow 25.10.2** o posterior **con el analizador de sintaxis v2**, que es el predeterminado a partir de Nextflow 26.04.
En nuestro entorno de entrenamiento no necesita hacer nada: ejecuta Nextflow 26.04.4 con el analizador v2. Si está usando un entorno local o personalizado, consulte las [notas de versión](../info/nxf_versions.md).

El entrenamiento además requiere **nf-core tools 4.0.2**.
Si usa una versión diferente de las herramientas nf-core, puede tener dificultades para seguir el curso.

Puede verificar qué versión está instalada en su entorno usando el comando `nf-core --version`.

!!! warning "Compatibilidad con el analizador v2"

    Muchos pipelines de nf-core aún no son compatibles con el analizador de sintaxis v2.
    Si ejecuta un pipeline de nf-core distinto de los usados en este curso y encuentra errores, puede que necesite cambiar al analizador v1 configurando `export NXF_SYNTAX_PARSER=v1`.
    Consulte las [notas de versión](../info/nxf_versions.md) para más detalles.

## Prepararse para trabajar

Una vez que su codespace esté funcionando, hay dos cosas que necesita hacer antes de sumergirse en la capacitación: establecer su directorio de trabajo para este curso específico y echar un vistazo a los materiales proporcionados.

### Establecer el directorio de trabajo

Por defecto, el codespace se abre con el directorio de trabajo establecido en la raíz de todos los cursos de capacitación, pero para este curso, trabajaremos en el directorio `hello-nf-core/`.

Cambie de directorio ahora ejecutando este comando en el terminal:

```bash
cd hello-nf-core/
```

!!! tip "Consejo"

    Si por alguna razón sale de este directorio (por ejemplo, su codespace se suspende), siempre puede usar la ruta completa para volver a él, asumiendo que está ejecutando esto dentro del entorno de capacitación de Github Codespaces:

    ```bash
    cd /workspaces/training/hello-nf-core
    ```

A continuación, explore el contenido de este directorio.

### Explorar los materiales proporcionados

Puede explorar el contenido de este directorio usando el explorador de archivos en el lado izquierdo del espacio de trabajo de capacitación.
Alternativamente, puede usar el comando `tree`.

A lo largo del curso, usamos la salida de `tree` para representar la estructura y contenidos del directorio en una forma legible, a veces con modificaciones menores para mayor claridad.

Aquí generamos una tabla de contenidos hasta el segundo nivel:

```bash
tree . -L 2
```

??? abstract "Contenido del directorio"

    ```console
    .
    ├── custom.config
    ├── greetings.csv
    ├── malformed_samplesheet.csv
    ├── my_params.yml
    ├── original-hello
    │   ├── hello.nf
    │   ├── modules
    │   └── nextflow.config
    └── solutions
        ├── composable-hello
        ├── core-hello-part2
        ├── core-hello-part3
        ├── core-hello-part4
        ├── core-hello-part5
        └── core-hello-start
    ```

Haga clic en el cuadro de color para expandir la sección y ver su contenido.
Usamos secciones colapsables como esta para incluir la salida esperada de comandos de manera concisa.

- **El archivo `greetings.csv`** es un CSV que contiene algunos datos columnares mínimos que usamos con fines de prueba.

- **El archivo `custom.config`** es un ejemplo de archivo de configuración de Nextflow usado en la Parte 1 para demostrar la sobreescritura de recursos de procesos y `ext.args`.

- **El archivo `malformed_samplesheet.csv`** es una hoja de muestras intencionalmente incorrecta usada en la Parte 1 para demostrar la validación de entradas.

- **El archivo `my_params.yml`** es un ejemplo de archivo de parámetros usado en la Parte 1 para demostrar cómo pasar parámetros booleanos a un pipeline.

- **El directorio `original-hello`** contiene una copia del código fuente producido al trabajar a través de la serie completa de capacitación Hello Nextflow (con Docker habilitado).

- **El directorio `solutions`** contiene los scripts de workflow completados que resultan de cada paso del curso.
  Están destinados a ser usados como referencia para verificar su trabajo y solucionar cualquier problema.

## Lista de verificación de preparación

¿Cree que está listo para sumergirse?

- [ ] Entiendo el objetivo de este curso y sus requisitos previos
- [ ] Mi entorno está funcionando
- [ ] Estoy usando nf-core tools 4.0.2 (verifique con `nf-core --version`)
- [ ] He establecido mi directorio de trabajo apropiadamente

Si puede marcar todas las casillas, está listo para comenzar.

**Para continuar a la Parte 1, haga clic en la flecha en la esquina inferior derecha de esta página.**
