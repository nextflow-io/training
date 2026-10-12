# Parte 3: Ejecutar un pipeline de producción

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Traducción asistida por IA - [más información y sugerencias](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

En la [Parte 2](./02_configure_execution.md), aprendió a configurar parámetros y personalizar la configuración para nf-core/demo.
Ahora aplicamos lo aprendido a un pipeline de producción real, nf-core/rnaseq.

---

## 1. Descargar y ejecutar nf-core/rnaseq

Hasta ahora hemos usado `nf-core/demo`, que es un pipeline mínimo diseñado para capacitación.
Ahora descargamos un pipeline de producción real y lo ejecutamos con su perfil de prueba.

El pipeline `nf-core/rnaseq` realiza los pasos principales del análisis de secuenciación de RNA en bulk: control de calidad, recorte de adaptadores, alineamiento de lecturas y cuantificación a nivel de genes.
Es probablemente el pipeline de nf-core más utilizado hasta la fecha.

### 1.1. Descargar el pipeline

Ejecute el siguiente comando para descargarlo.

```bash
nextflow pull nf-core/rnaseq
```

??? success "Salida del comando"

    ```console
    Checking nf-core/rnaseq ...
     downloaded from https://github.com/nf-core/rnaseq.git - revision: 1f03b53ef7 [master]
    ```

El pipeline ahora está almacenado en caché localmente y listo para ejecutarse.

### 1.2. Ejecutar el perfil de prueba

Ejecútelo con el perfil de prueba y Docker:

```bash
nextflow run nf-core/rnaseq -profile test,docker --outdir rnaseq-results
```

??? failure "Salida del comando"

    ```console
     N E X T F L O W   ~  version 26.04.4

    Launching `https://github.com/nf-core/rnaseq` [suspicious_dijkstra] revision: 1f03b53ef7 [master]

    ------------------------------------------------------
                                            ,--./,-.
            ___     __   __   __   ___     /,-._.--~'
      |\ | |__  __ /  ` /  \ |__) |__         }  {
      | \| |       \__, \__/ |  \ |___     \`-._,-`-,
                                            `._,._,'
      nf-core/rnaseq 3.26.0
    ------------------------------------------------------
    ...
    [-        ] NFCORE_RNASEQ:RNASEQ:FASTQ_QC_TRIM_FILTER_SETSTRANDEDNESS:FQ_LINT  -
    Plus 47 more processes waiting for tasks…

    Execution cancelled -- Finishing pending tasks before exit
    -[nf-core/rnaseq] Pipeline completed with errors-
    ERROR ~ Error executing process > 'NFCORE_RNASEQ:RNASEQ:FASTQ_QC_TRIM_FILTER_SETSTRANDEDNESS:FQ_LINT (WT_REP2)'

    Caused by:
      Process requirement exceeds available memory -- req: 12 GB; avail: 7.7 GB


    Command executed:

      fq lint \
          --disable-validator P001 \
          SRR6357072_1.fastq.gz SRR6357072_2.fastq.gz > WT_REP2.fq_lint.txt

    Command exit status:
      -

    Command output:
      (empty)

    Work dir:
      /workspaces/training/nfcore-use/work/xx/xxxxxxxxxxxxxxxxxxxxxx

    Container:
      quay.io/biocontainers/fq:0.12.0--h9ee0642_0

    Tip: you can replicate the issue by changing to the process work dir and entering the command `bash .command.run`

     -- Check '.nextflow.log' file for details
    ERROR ~ Pipeline failed. Please refer to troubleshooting docs for common issues: https://nf-co.re/docs/running/troubleshooting
    ```

La línea clave en ese error es:

```console
Process requirement exceeds available memory -- req: 12 GB; avail: 7.7 GB
```

La máquina de Codespaces predeterminada tiene 8 GB de RAM, que también es el valor predeterminado típico para Docker Desktop.
El pipeline está solicitando 12 GB para el proceso `FQ_LINT` — más de lo que la máquina puede proporcionar.

Esos 12 GB provienen de la etiqueta de recursos `process_low` definida en `conf/base.config`:

```groovy title="conf/base.config"
withLabel:process_low {
    cpus   = { 2     * task.attempt }
    memory = { 12.GB * task.attempt }
    time   = { 4.h   * task.attempt }
}
```

Una opción sería usar un tipo de máquina más grande, pero para propósitos de prueba queremos poder ejecutar en cualquier hardware disponible.
El mejor enfoque es sobrescribir los valores predeterminados de recursos en un archivo de configuración personalizado.

### 1.3. Volver a ejecutar con una configuración personalizada

Le proporcionamos un archivo de configuración personalizado que sobrescribe los valores predeterminados de recursos basados en etiquetas.

??? full-code "laptop.config"

    ```groovy title="laptop.config"
    process {
        withLabel: 'process_low' {
            cpus   = 2
            memory = 6.GB
        }
        withLabel: 'process_medium' {
            cpus   = 4
            memory = 6.GB
        }
        withLabel: 'process_high' {
            cpus   = 6
            memory = 6.GB
        }
        withLabel: 'process_high_memory' {
            memory = 6.GB
        }
    }
    ```

La [Parte 2](./02_configure_execution.md) introdujo `withName:` para apuntar a un único proceso por nombre.
Aquí usamos `withLabel:` para apuntar a todos los procesos que comparten una etiqueta a la vez.

Este archivo ya está presente en su directorio de trabajo.
Páselo con `-c` para aplicar las sobrescrituras:

```bash
nextflow run nf-core/rnaseq -profile test,docker -c laptop.config --outdir rnaseq-results
```

??? success "Salida del comando (pipeline iniciando)"

    ```console
     N E X T F L O W   ~  version 26.04.4

    Launching `https://github.com/nf-core/rnaseq` [romantic_faraday] revision: 1f03b53ef7 [master]

    ------------------------------------------------------
                                            ,--./,-.
            ___     __   __   __   ___     /,-._.--~'
      |\ | |__  __ /  ` /  \ |__) |__         }  {
      | \| |       \__, \__/ |  \ |___     \`-._,-`-,
                                            `._,._,'
      nf-core/rnaseq 3.26.0
    ------------------------------------------------------
    ...
    executor >  local (7)
    [xx/xxxxxx] NFCORE_RNASEQ:RNASEQ:FQ_LINT (RAP1_IAA_30M_REP1)   | 3 of 5
    [xx/xxxxxx] NFCORE_RNASEQ:RNASEQ:FASTQC (RAP1_IAA_30M_REP1)    | 2 of 5
    ...
    ```

El pipeline ahora está en ejecución y puede observar cómo las tareas se completan una por una.
Con este conjunto de datos de prueba mínimo, se completará en 15–20 minutos, ejecutando más de 200 tareas en total.

Los experimentos reales de RNA-seq típicamente involucran decenas de muestras y se ejecutan durante horas o días.
Nextflow es compatible con schedulers de HPC (SLURM, PBS, LSF) y plataformas en la nube (AWS, Google Cloud, Azure), que pueden reducir drásticamente el tiempo de ejecución al distribuir el trabajo entre muchos nodos.
Sin embargo, configurar esos entornos añade una complejidad considerable.

La plataforma Seqera (desarrollada por los creadores de Nextflow) proporciona una interfaz web para lanzar pipelines de Nextflow en infraestructura HPC o en la nube (ya sea la propia o una gestionada para usted), con capacidades de gestión de cómputo y datos que simplifican el proceso de ejecutar pipelines a escala.

!!! tip "Consejo"

    Los investigadores académicos pueden acceder a Seqera Platform de forma gratuita a través del [programa académico de Seqera](https://seqera.io/academic-program/).

### Conclusión

Ha descargado `nf-core/rnaseq`, visto cómo funcionan las etiquetas de recursos de nf-core y aprendido a sobrescribirlas con un archivo de configuración personalizado.
Más importante aún, ha comprendido por qué la ejecución local es un punto de partida y no un destino para análisis a escala real.

### ¿Qué sigue?

Ha cubierto los fundamentos de la ejecución de pipelines de nf-core.
Consulte los [Próximos pasos](next_steps.md) para saber hacia dónde continuar.

---

## Resumen

En esta parte aprendió a:

- Descargar y ejecutar un pipeline a escala de producción (nf-core/rnaseq) y sobrescribir sus etiquetas de recursos predeterminadas
