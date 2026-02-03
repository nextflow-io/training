# Orientación

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Traducción asistida por IA - [más información y sugerencias](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

El entorno de entrenamiento contiene todo el software, código y datos necesarios para trabajar en este curso de entrenamiento, por lo que no necesita instalar nada por su cuenta.
Sin embargo, sí necesita una cuenta (gratuita) para iniciar sesión, y debería tomarse unos minutos para familiarizarse con la interfaz.

Si aún no lo ha hecho, por favor siga [este enlace](../../../envsetup/) antes de continuar.

## Materiales proporcionados

A lo largo de este curso de entrenamiento, trabajaremos en el directorio `nf4-science/genomics/`, al cual debe moverse cuando abra el espacio de trabajo de entrenamiento.
Este directorio contiene todos los archivos de código, datos de prueba y archivos accesorios que necesitará.

Siéntase libre de explorar el contenido de este directorio; la forma más fácil de hacerlo es usar el explorador de archivos en el lado izquierdo del espacio de trabajo de entrenamiento en la interfaz de VSCode.
Alternativamente, puede usar el comando `tree`.
A lo largo del curso, usamos la salida de `tree` para representar la estructura y contenido del directorio de una forma legible, a veces con modificaciones menores para mayor claridad.

Aquí generamos una tabla de contenidos hasta el segundo nivel:

```bash
tree . -L 2
```

Si ejecuta esto dentro de `nf4-science/genomics`, debería ver la siguiente salida:

```console title="Contenido del directorio"

.
├── data
│   ├── bam
│   ├── ref
│   ├── sample_bams.txt
│   └── samplesheet.csv
├── genomics-1.nf
├── genomics-2.nf
├── genomics-3.nf
├── genomics-4.nf
├── nextflow.config
└── solutions
    ├── modules
    ├── nf-test.config
    └── tests

6 directories, 8 files

```

!!!note "Nota"

    No se preocupe si esto parece mucho; revisaremos las partes relevantes en cada paso del curso.
    Esto está destinado solo a darle una visión general.

**Aquí hay un resumen de lo que debe saber para comenzar:**

- **Los archivos `.nf`** son scripts de workflow que se nombran según la parte del curso en la que se utilizan.

- **El archivo `nextflow.config`** es un archivo de configuración que establece propiedades mínimas del entorno.
  Puede ignorarlo por ahora.

- **El directorio `data`** contiene datos de entrada y recursos relacionados, descritos más adelante en el curso.

- **El directorio `solutions`** contiene archivos de módulo y configuraciones de prueba que resultan de las Partes 3 y 4 del curso.
  Están destinados a ser usados como referencia para verificar su trabajo y solucionar cualquier problema.

!!!tip "Consejo"

    Si por cualquier razón sale de este directorio, siempre puede ejecutar este comando para regresar a él:

    ```bash
    cd /workspaces/training/nf4-science/genomics
    ```

Ahora, para comenzar el curso, haga clic en la flecha en la esquina inferior derecha de esta página.
