# Orientación

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Traducción asistida por IA - [más información y sugerencias](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

El entorno de entrenamiento contiene todo el software, código y datos necesarios para trabajar en este curso de entrenamiento, por lo que no necesita instalar nada usted mismo.
Sin embargo, sí necesita una cuenta (gratuita) para iniciar sesión, y debe tomarse unos minutos para familiarizarse con la interfaz.

Si aún no lo ha hecho, por favor complete el mini-curso de [Configuración del Entorno](../../envsetup/) antes de continuar.

## Materiales proporcionados

A lo largo de este curso de entrenamiento, trabajaremos en el directorio `nf4-science/rnaseq/`, al cual debe moverse cuando abra el espacio de trabajo de entrenamiento.
Este directorio contiene todos los archivos de código, datos de prueba y archivos accesorios que necesitará.

Siéntase libre de explorar el contenido de este directorio; la forma más fácil de hacerlo es usar el explorador de archivos en el lado izquierdo del espacio de trabajo de entrenamiento en la interfaz de VSCode.
Alternativamente, puede usar el comando `tree`.
A lo largo del curso, usamos la salida de `tree` para representar la estructura y contenido del directorio en una forma legible, a veces con modificaciones menores para mayor claridad.

Aquí generamos una tabla de contenidos hasta el segundo nivel:

```bash
tree . -L 3
```

??? success "Contenido del directorio"

    ```console
    rnaseq
    ├── data
    │   ├── genome.fa
    │   ├── paired-end.csv
    │   ├── reads
    │   │   ├── ENCSR000COQ1_1.fastq.gz
    │   │   ├── ENCSR000COQ1_2.fastq.gz
    │   │   ├── ENCSR000COQ2_1.fastq.gz
    │   │   ├── ENCSR000COQ2_2.fastq.gz
    │   │   ├── ENCSR000COR1_1.fastq.gz
    │   │   ├── ENCSR000COR1_2.fastq.gz
    │   │   ├── ENCSR000COR2_1.fastq.gz
    │   │   ├── ENCSR000COR2_2.fastq.gz
    │   │   ├── ENCSR000CPO1_1.fastq.gz
    │   │   ├── ENCSR000CPO1_2.fastq.gz
    │   │   ├── ENCSR000CPO2_1.fastq.gz
    │   │   └── ENCSR000CPO2_2.fastq.gz
    │   └── single-end.csv
    ├── nextflow.config
    ├── rnaseq.nf
    └── solutions
        ├── modules
        │   ├── fastqc.nf
        │   ├── fastqc_pe.nf
        │   ├── hisat2_align.nf
        │   ├── hisat2_align_pe.nf
        │   ├── multiqc.nf
        │   ├── trim_galore.nf
        │   └── trim_galore_pe.nf
        ├── rnaseq-2.1.nf
        ├── rnaseq-2.2.nf
        ├── rnaseq-2.3.nf
        ├── rnaseq-3.1.nf
        ├── rnaseq-3.2.nf
        └── rnaseq_pe-3.3.nf
    ```

!!!note

    No se preocupe si esto parece mucho; revisaremos las partes relevantes en cada paso del curso.
    Esto es solo para darle una visión general.

**Aquí hay un resumen de lo que debe saber para comenzar:**

- **El archivo `rnaseq.nf`** es el esquema del script de flujo de trabajo que desarrollaremos.

- **El archivo `nextflow.config`** es un archivo de configuración que establece propiedades mínimas del entorno. Puede ignorarlo por ahora.

- **El directorio `data`** contiene datos de entrada y recursos relacionados:

  - _Un genoma de referencia_ llamado `genome.fa` que consiste en una pequeña región del cromosoma 20 humano (de hg19/b37).
  - _Datos de RNAseq_ que han sido reducidos a una pequeña región para mantener los tamaños de archivo pequeños, en el directorio `reads/`.
  - _Archivos CSV_ que listan los IDs y rutas de los archivos de datos de ejemplo, para procesamiento en lotes.

- **El directorio `solutions`** contiene los scripts de flujo de trabajo y módulos completados que resultan de cada paso del curso.
  Están destinados a ser usados como referencia para verificar su trabajo y solucionar cualquier problema.
  El número en el nombre del archivo corresponde al paso de la parte relevante del curso.

!!!tip

    Si por cualquier razón sale de este directorio, siempre puede ejecutar este comando para regresar:

    ```bash
    cd /workspaces/training/nf4-science/rnaseq
    ```

Ahora, para comenzar el curso, haga clic en la flecha en la esquina inferior derecha de esta página.
