# Primeros pasos

## Iniciar un entorno de capacitación

Para usar el entorno preconfigurado que proporcionamos en GitHub Codespaces, haz clic en el botón "Open in GitHub Codespaces" a continuación. Para otras opciones, consulta [Opciones de entorno](../envsetup/index.md).

Te recomendamos abrir el entorno de capacitación en una nueva pestaña o ventana del navegador (usa clic derecho, ctrl-clic o cmd-clic según tu equipo) para que puedas seguir leyendo mientras se carga el entorno.
Necesitarás mantener estas instrucciones abiertas en paralelo para trabajar en el curso.

[![Open in GitHub Codespaces](https://github.com/codespaces/badge.svg)](https://codespaces.new/nextflow-io/training?quickstart=1&ref=master)

### Conceptos básicos del entorno

Este entorno de capacitación contiene todo el software, código y datos necesarios para trabajar en el curso de capacitación, por lo que no necesitas instalar nada por tu cuenta.

El codespace está configurado con una interfaz VSCode, que incluye un explorador de archivos, un editor de código y una terminal shell.
Todas las instrucciones dadas durante el curso (por ejemplo, 'abre el archivo', 'edita el código' o 'ejecuta este comando') se refieren a esas tres partes de la interfaz VSCode a menos que se especifique lo contrario.

Si estás trabajando en este curso por tu cuenta, por favor familiarízate con los [conceptos básicos del entorno](../envsetup/01_setup.md) para más detalles.

### Requisitos de versión

Esta capacitación está diseñada para **Nextflow 25.10.2** o posterior **con el analizador de sintaxis v2 DESHABILITADO**.

#### Si estás usando nuestro entorno de capacitación:

DEBES ejecutar el siguiente comando antes de continuar:

```bash
export NXF_SYNTAX_PARSER=v1
```

#### Si estás usando un entorno local o personalizado:

Por favor asegúrate de estar usando la configuración correcta como se documenta [aquí](../info/nxf_versions.md).

La capacitación además requiere **nf-core tools 3.4.1**.
Si usas una versión diferente de las herramientas nf-core, podrías tener dificultades para seguir el curso.

Puedes verificar qué versión está instalada en tu entorno usando el comando `nf-core --version`.

## Prepárate para trabajar

Una vez que tu codespace esté ejecutándose, hay dos cosas que necesitas hacer antes de sumergirte en la capacitación: establecer tu directorio de trabajo para este curso específico y revisar los materiales proporcionados.

### Establecer el directorio de trabajo

Por defecto, el codespace se abre con el directorio de trabajo establecido en la raíz de todos los cursos de capacitación, pero para este curso, trabajaremos en el directorio `hello-nf-core/`.

Cambia de directorio ahora ejecutando este comando en la terminal:

```bash
cd hello-nf-core/
```

!!! tip "Consejo"

    Si por alguna razón sales de este directorio (por ejemplo, tu codespace se suspende), siempre puedes usar la ruta completa para regresar a él, asumiendo que estás ejecutando esto dentro del entorno de capacitación de Github Codespaces:

    ```bash
    cd /workspaces/training/hello-nf-core
    ```

Ahora echemos un vistazo al contenido de este directorio.

### Explorar los materiales proporcionados

Puedes explorar el contenido de este directorio usando el explorador de archivos en el lado izquierdo del espacio de trabajo de capacitación.
Alternativamente, puedes usar el comando `tree`.

A lo largo del curso, usamos la salida de `tree` para representar la estructura y contenido del directorio de forma legible, a veces con modificaciones menores para mayor claridad.

Aquí generamos una tabla de contenidos hasta el segundo nivel:

```bash
tree . -L 2
```

??? abstract "Contenido del directorio"

    ```console
    .
    ├── greetings.csv
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

Haz clic en el cuadro de color para expandir la sección y ver su contenido.
Usamos secciones plegables como esta para incluir la salida esperada de comandos de manera concisa.

- **El archivo `greetings.csv`** es un CSV que contiene algunos datos columnares mínimos que usamos para propósitos de prueba.

- **El directorio `original-hello`** contiene una copia del código fuente producido al trabajar en la serie completa de capacitación Hello Nextflow (con Docker habilitado).

- **El directorio `solutions`** contiene los scripts de workflow completados que resultan de cada paso del curso.
  Están destinados a ser usados como referencia para verificar tu trabajo y solucionar cualquier problema.

## Lista de verificación de preparación

¿Crees que estás listo para comenzar?

- [ ] Entiendo el objetivo de este curso y sus requisitos previos
- [ ] Mi entorno está funcionando
- [ ] Me he asegurado de que el analizador de sintaxis esté configurado en **v1**
- [ ] He establecido mi directorio de trabajo apropiadamente

Si puedes marcar todas las casillas, estás listo para comenzar.

**Para continuar a la Parte 1, haz clic en la flecha en la esquina inferior derecha de esta página.**
