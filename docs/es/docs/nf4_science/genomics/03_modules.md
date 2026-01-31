# Parte 3: Mover código a módulos

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Traducción asistida por IA - [más información y sugerencias](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

En la primera parte de este curso, construiste un pipeline de llamado de variantes que era completamente lineal y procesaba los datos de cada muestra independientemente de las demás.

En la segunda parte, te mostramos cómo usar canales y operadores de canal para implementar el llamado conjunto de variantes con GATK, construyendo sobre el pipeline de la Parte 1.

En esta parte, te mostraremos cómo convertir el código de ese workflow en módulos. Para seguir esta parte del entrenamiento, debes haber completado la Parte 1 y la Parte 2, así como [Hello Modules](../../../hello_nextflow/hello_modules.md), que cubre los conceptos básicos de los módulos.

---

## 0. Calentamiento

Cuando comenzamos a desarrollar nuestro workflow, pusimos todo en un único archivo de código.
Ahora es momento de abordar la **modularización** de nuestro código, _es decir_, extraer las definiciones de procesos en módulos.

Vamos a comenzar con el mismo workflow que en la Parte 2, que hemos proporcionado en el archivo `genomics-3.nf`.

!!! note "Nota"

     Asegúrate de estar en el directorio de trabajo correcto:
     `cd /workspaces/training/nf4-science/genomics`

Ejecuta el workflow para verificar el punto de partida:

```bash
nextflow run genomics-3.nf -resume
```

```console title="Salida"
 N E X T F L O W   ~  version 25.10.2

Launching `genomics-3.nf` [serene_borg] DSL2 - revision: 0cbebb67a1

executor >  local (7)
[6f/83ee72] SAMTOOLS_INDEX (3)       | 3 of 3 ✔
[53/b9d342] GATK_HAPLOTYPECALLER (1) | 3 of 3 ✔
[0c/fa6d15] GATK_JOINTGENOTYPING     | 1 of 1 ✔
```

Ahora habrá un directorio `work` y un directorio `results_genomics` dentro de tu directorio del proyecto.

### Conclusión

Estás listo para comenzar a modularizar tu workflow.

### ¿Qué sigue?

Mover los procesos del workflow de Genómica a módulos.

---

## 1. Mover procesos a módulos

Como aprendiste en [Hello Modules](../../../hello_nextflow/hello_modules.md), puedes crear un módulo simplemente copiando la definición del proceso en su propio archivo, en cualquier directorio, y puedes nombrar ese archivo como quieras.

Por razones que quedarán claras más adelante (en particular cuando lleguemos a las pruebas), en este entrenamiento seguiremos la convención de nombrar el archivo `main.nf`, y colocarlo en una estructura de directorios nombrada según el conjunto de herramientas y el comando.

### 1.1. Crear un módulo para el proceso `SAMTOOLS_INDEX`

En el caso del proceso `SAMTOOLS_INDEX`, 'samtools' es el conjunto de herramientas e 'index' es el comando. Entonces, crearemos una estructura de directorios `modules/samtools/index` y pondremos la definición del proceso `SAMTOOLS_INDEX` en el archivo `main.nf` dentro de ese directorio.

```bash
mkdir -p modules/samtools/index
touch modules/samtools/index/main.nf
```

Abre el archivo `main.nf` y copia la definición del proceso `SAMTOOLS_INDEX` en él.

```groovy title="modules/samtools/index/main.nf" linenums="1"
/*
 * Generar archivo de índice BAM
 */
process SAMTOOLS_INDEX {

    container 'community.wave.seqera.io/library/samtools:1.20--b5dfbd93de237464'

    input:
    path input_bam

    output:
    tuple path(input_bam), path("${input_bam}.bai")

    script:
    """
    samtools index '$input_bam'
    """
}
```

Luego, elimina la definición del proceso `SAMTOOLS_INDEX` de `genomics-3.nf`, y agrega una declaración de importación para el módulo antes de la siguiente definición de proceso, así:

=== "Después"

    ```groovy title="genomics-3.nf" linenums="1" hl_lines="1 2"
    // Incluir módulos
    include { SAMTOOLS_INDEX } from './modules/samtools/index/main.nf'

    /*
     * Llamar variantes con GATK HaplotypeCaller
     */
    process GATK_HAPLOTYPECALLER {
    ```

=== "Antes"

    ```groovy title="genomics-3.nf" linenums="1" hl_lines="1"
    /*
     * Llamar variantes con GATK HaplotypeCaller
     */
    process GATK_HAPLOTYPECALLER {
    ```

Ahora puedes ejecutar el workflow nuevamente, y debería funcionar de la misma manera que antes. Si proporcionas la bandera `-resume`, ni siquiera deberían ejecutarse nuevas tareas:

```bash
nextflow run genomics-3.nf -resume
```

??? success "Salida del comando"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `genomics-3.nf` [sleepy_snyder] DSL2 - revision: aa68d06c43

    [0f/71b55e] SAMTOOLS_INDEX (1)       | 3 of 3, cached: 3 ✔
    [f1/18971b] GATK_HAPLOTYPECALLER (3) | 3 of 3, cached: 3 ✔
    [0c/fa6d15] GATK_JOINTGENOTYPING     | 1 of 1, cached: 1 ✔
    ```

### 1.2. Crear módulos para los procesos `GATK_HAPLOTYPECALLER` y `GATK_JOINTGENOTYPING`

Repite los mismos pasos para los procesos restantes.
Para cada proceso:

1. Crea la estructura de directorios (`modules/gatk/haplotypecaller/` y `modules/gatk/jointgenotyping/`)
2. Crea un archivo `main.nf` que contenga la definición del proceso
3. Elimina la definición del proceso de `genomics-3.nf`
4. Agrega una declaración de importación para el módulo

Una vez que hayas terminado, verifica que la estructura de directorios de tus módulos sea correcta ejecutando:

```bash
tree modules/
```

??? abstract "Contenidos del directorio"

    ```console
    modules/
    ├── gatk
    │   ├── haplotypecaller
    │   │   └── main.nf
    │   └── jointgenotyping
    │       └── main.nf
    └── samtools
        └── index
            └── main.nf

    5 directories, 3 files
    ```

También deberías tener algo como esto en el archivo principal del workflow, después de la sección de parámetros:

```
include { SAMTOOLS_INDEX } from './modules/samtools/index/main.nf'
include { GATK_HAPLOTYPECALLER } from './modules/gatk/haplotypecaller/main.nf'
include { GATK_JOINTGENOTYPING } from './modules/gatk/jointgenotyping/main.nf'

workflow {
```

### Conclusión

Has practicado la modularización de un workflow, con el workflow de genómica como ejemplo.

### ¿Qué sigue?

Probar el workflow modularizado.

---

## 2. Probar el workflow modularizado

Ejecuta el workflow modularizado para verificar que todo siga funcionando.

```bash
nextflow run genomics-3.nf -resume
```

```console title="Salida"
 N E X T F L O W   ~  version 25.10.2

Launching `genomics-3.nf` [astonishing_venter] DSL2 - revision: ca27264c13

[6f/83ee72] SAMTOOLS_INDEX (3)       | 3 of 3, cached: 3 ✔
[53/b9d342] GATK_HAPLOTYPECALLER (3) | 3 of 3, cached: 3 ✔
[0c/fa6d15] GATK_JOINTGENOTYPING     | 1 of 1, cached: 1 ✔
```

Todo sigue funcionando, incluyendo la capacidad de reanudación del pipeline.
Los resultados continúan siendo publicados en el directorio `results_genomics`.

```console title="Contenidos del directorio"
results_genomics/
├── family_trio.joint.vcf
├── family_trio.joint.vcf.idx
├── gvcf
│   ├── reads_father.bam.g.vcf
│   ├── reads_father.bam.g.vcf.idx
│   ├── reads_mother.bam.g.vcf
│   ├── reads_mother.bam.g.vcf.idx
│   ├── reads_son.bam.g.vcf
│   └── reads_son.bam.g.vcf.idx
└── indexed_bam
    ├── reads_father.bam
    ├── reads_father.bam.bai
    ├── reads_mother.bam
    ├── reads_mother.bam.bai
    ├── reads_son.bam
    └── reads_son.bam.bai
```

### Conclusión

Has modularizado un workflow y verificado que sigue funcionando de la misma manera que antes.

### ¿Qué sigue?

Revisar lo que has aprendido y mirar hacia adelante a las pruebas.

---

## 3. Resumen

Has modularizado el workflow, y nada ha cambiado en cómo funciona el pipeline.
Esto es intencional: has reestructurado el código sin impactar su función.

Los módulos contienen solo la lógica del proceso, haciéndolos limpios y reutilizables.
El script principal controla qué se publica y dónde, mientras que los módulos permanecen enfocados en su tarea computacional.

Has sentado las bases para cosas que harán tu código más fácil de mantener.
Por ejemplo, ahora puedes agregar pruebas a tu pipeline usando el framework nf-test.
Esto es lo que veremos en la siguiente parte de este curso.
