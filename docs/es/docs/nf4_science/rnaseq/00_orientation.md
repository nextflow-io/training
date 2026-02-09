# Orientación

El entorno de capacitación contiene todo el software, código y datos necesarios para trabajar en este curso de capacitación, por lo que no necesitas instalar nada por tu cuenta.
Sin embargo, sí necesitas una cuenta (gratuita) para iniciar sesión, y deberías tomarte unos minutos para familiarizarte con la interfaz.

Si aún no lo has hecho, por favor completa el mini-curso de [Configuración del Entorno](../../envsetup/) antes de continuar.

## Materiales proporcionados

A lo largo de este curso de capacitación, trabajaremos en el directorio `nf4-science/rnaseq/`, al cual debes moverte cuando abras el espacio de trabajo de capacitación.
Este directorio contiene todos los archivos de código, datos de prueba y archivos accesorios que necesitarás.

Siéntete libre de explorar el contenido de este directorio; la forma más fácil de hacerlo es usar el explorador de archivos en el lado izquierdo del espacio de trabajo de capacitación en la interfaz de VSCode.
Alternativamente, puedes usar el comando `tree`.
A lo largo del curso, usamos la salida de `tree` para representar la estructura y contenido del directorio de forma legible, a veces con modificaciones menores para mayor claridad.

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

    No te preocupes si esto parece mucho; revisaremos las partes relevantes en cada paso del curso.
    Esto es solo para darte una visión general.

**Aquí hay un resumen de lo que debes saber para comenzar:**

- **El archivo `rnaseq.nf`** es el esquema del script de workflow que trabajaremos para desarrollar.

- **El archivo `nextflow.config`** es un archivo de configuración que establece propiedades mínimas del entorno. Puedes ignorarlo por ahora.

- **El directorio `data`** contiene datos de entrada y recursos relacionados:

  - _Un genoma de referencia_ llamado `genome.fa` que consiste en una pequeña región del cromosoma 20 humano (de hg19/b37).
  - _Datos de RNAseq_ que han sido reducidos a una pequeña región para mantener los tamaños de archivo pequeños, en el directorio `reads/`.
  - _Archivos CSV_ que listan los IDs y rutas de los archivos de datos de ejemplo, para procesamiento por lotes.

- **El directorio `solutions`** contiene los scripts de workflow completos y módulos que resultan de cada paso del curso.
  Están destinados a ser usados como referencia para verificar tu trabajo y solucionar cualquier problema.
  El número en el nombre del archivo corresponde al paso de la parte relevante del curso.

!!!tip

    Si por alguna razón te sales de este directorio, siempre puedes ejecutar este comando para regresar a él:

    ```bash
    cd /workspaces/training/nf4-science/rnaseq
    ```

Ahora, para comenzar el curso, haz clic en la flecha en la esquina inferior derecha de esta página.
