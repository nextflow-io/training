# División y Agrupación

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Traducción asistida por IA - [más información y sugerencias](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

Nextflow proporciona herramientas poderosas para trabajar con datos de manera flexible. Una capacidad clave es dividir datos en diferentes flujos y luego agrupar elementos relacionados. Esto es especialmente valioso en workflows de bioinformática donde necesitas procesar diferentes tipos de muestras por separado antes de combinar los resultados para el análisis.

Piénsalo como clasificar el correo: separas las cartas por destino, procesas cada pila de manera diferente y luego recombinas los elementos que van a la misma persona. Nextflow usa operadores especiales para lograr esto con datos científicos. Este enfoque también se conoce comúnmente como el patrón **scatter/gather** en computación distribuida y workflows de bioinformática.

El sistema de canales de Nextflow es el corazón de esta flexibilidad. Los canales conectan diferentes partes de tu workflow, permitiendo que los datos fluyan a través de tu análisis. Puedes crear múltiples canales desde una sola fuente de datos, procesar cada canal de manera diferente y luego combinar los canales cuando sea necesario. Este enfoque te permite diseñar workflows que reflejen naturalmente los caminos ramificados y convergentes de los análisis complejos de bioinformática.

### Objetivos de aprendizaje

En esta misión secundaria, aprenderás a dividir y agrupar datos usando los operadores de canales de Nextflow.
Comenzaremos con un archivo CSV que contiene información de muestras y archivos de datos asociados, luego manipularemos y reorganizaremos estos datos.

Al final de esta misión secundaria, podrás separar y combinar flujos de datos de manera efectiva, usando las siguientes técnicas:

- Leer datos de archivos usando `splitCsv`
- Filtrar y transformar datos con `filter` y `map`
- Combinar datos relacionados usando `join` y `groupTuple`
- Crear combinaciones de datos con `combine` para procesamiento paralelo
- Optimizar la estructura de datos usando `subMap` y estrategias de deduplicación
- Construir funciones reutilizables con closures con nombre para manipular estructuras de canales

Estas habilidades te ayudarán a construir workflows que puedan manejar múltiples archivos de entrada y diferentes tipos de datos de manera eficiente, manteniendo una estructura de código limpia y fácil de mantener.

### Requisitos previos

Antes de comenzar esta misión secundaria, debes:

- Haber completado el tutorial [Hello Nextflow](../hello_nextflow/README.md) o un curso equivalente para principiantes.
- Estar familiarizado con los conceptos y mecanismos básicos de Nextflow (procesos, canales, operadores, trabajo con archivos, metadatos)

**Opcional:** Recomendamos completar primero la misión secundaria [Metadata in workflows](../metadata/).
Esta cubre los fundamentos de lectura de archivos CSV con `splitCsv` y la creación de mapas de metadatos, que usaremos ampliamente aquí.

---

## 0. Primeros pasos

#### Abrir el codespace de capacitación

Si aún no lo has hecho, asegúrate de abrir el entorno de capacitación como se describe en la [Configuración del entorno](../envsetup/index.md).

[![Open in GitHub Codespaces](https://github.com/codespaces/badge.svg)](https://codespaces.new/nextflow-io/training?quickstart=1&ref=master)

#### Moverse al directorio del proyecto

Vamos a movernos al directorio donde se encuentran los archivos de este tutorial.

```bash
cd side-quests/splitting_and_grouping
```

Puedes configurar VSCode para que se enfoque en este directorio:

```bash
code .
```

#### Revisar los materiales

Encontrarás un archivo de workflow principal y un directorio `data` que contiene una hoja de muestras llamada `samplesheet.csv`.

```console title="Directory contents"
.
├── data
│   └── samplesheet.csv
└── main.nf
```

La hoja de muestras contiene información sobre muestras de diferentes pacientes, incluyendo el ID del paciente, el número de repetición de la muestra, el tipo (normal o tumor) y las rutas a archivos de datos hipotéticos (que en realidad no existen, pero fingiremos que sí).

```console title="samplesheet.csv"
id,repeat,type,bam
patientA,1,normal,patientA_rep1_normal.bam
patientA,1,tumor,patientA_rep1_tumor.bam
patientA,2,normal,patientA_rep2_normal.bam
patientA,2,tumor,patientA_rep2_tumor.bam
patientB,1,normal,patientB_rep1_normal.bam
patientB,1,tumor,patientB_rep1_tumor.bam
patientC,1,normal,patientC_rep1_normal.bam
patientC,1,tumor,patientC_rep1_tumor.bam
```

Esta hoja de muestras lista ocho muestras de tres pacientes (A, B, C).

Para cada paciente, tenemos muestras de tipo `tumor` (que típicamente provienen de biopsias tumorales) o `normal` (tomadas de tejido sano o sangre).
Si no estás familiarizado con el análisis de cáncer, solo debes saber que esto corresponde a un modelo experimental que usa pares de muestras tumor/normal para realizar análisis contrastivos.

Para el paciente A específicamente, tenemos dos conjuntos de réplicas técnicas (repeticiones).

!!! note "Nota"

    No te preocupes si no estás familiarizado con este diseño experimental, no es fundamental para entender este tutorial.

#### Revisar la tarea

Tu desafío es escribir un workflow de Nextflow que:

1. **Lea** datos de muestras desde un archivo CSV y los estructure con mapas de metadatos
2. **Separe** las muestras en diferentes canales según el tipo (normal vs tumor)
3. **Una** los pares tumor/normal coincidentes por ID de paciente y número de réplica
4. **Distribuya** las muestras a través de intervalos genómicos para procesamiento paralelo
5. **Agrupe** las muestras relacionadas para el análisis posterior

Esto representa un patrón común en bioinformática donde necesitas dividir datos para procesamiento independiente y luego recombinar elementos relacionados para análisis comparativo.

#### Lista de verificación de preparación

¿Crees que estás listo para comenzar?

- [ ] Entiendo el objetivo de este curso y sus requisitos previos
- [ ] Mi codespace está en funcionamiento
- [ ] He configurado mi directorio de trabajo correctamente
- [ ] Entiendo la tarea

Si puedes marcar todas las casillas, estás listo para continuar.

---

## 1. Leer datos de muestras

### 1.1. Leer datos de muestras con `splitCsv` y crear mapas de metadatos

Comencemos leyendo los datos de muestras con `splitCsv` y organizándolos en el patrón de mapa de metadatos. En el archivo `main.nf`, verás que ya hemos comenzado el workflow.

```groovy title="main.nf" linenums="1" hl_lines="2"
workflow {
    ch_samplesheet = channel.fromPath("./data/samplesheet.csv")
}
```

!!! note "Nota"

    A lo largo de este tutorial, usaremos el prefijo `ch_` para todas las variables de canal para indicar claramente que son canales de Nextflow.

Si completaste la misión secundaria [Metadata in workflows](../metadata/), reconocerás este patrón. Usaremos `splitCsv` para leer el CSV e inmediatamente estructurar los datos con un mapa de metadatos para separar los metadatos de las rutas de archivos.

!!! info "Info"

    Encontraremos dos conceptos diferentes llamados `map` en esta capacitación:

    - **Estructura de datos**: El mapa de Groovy (equivalente a diccionarios/hashes en otros lenguajes) que almacena pares clave-valor
    - **Operador de canal**: El operador `.map()` que transforma elementos en un canal

    Aclararemos cuál de los dos queremos decir según el contexto, pero esta distinción es importante de entender cuando se trabaja con Nextflow.

Aplica estos cambios a `main.nf`:

=== "Después"

    ```groovy title="main.nf" linenums="2" hl_lines="2-6"
        ch_samples = channel.fromPath("./data/samplesheet.csv")
            .splitCsv(header: true)
            .map{ row ->
              [[id:row.id, repeat:row.repeat, type:row.type], row.bam]
            }
            .view()
    ```

=== "Antes"

    ```groovy title="main.nf" linenums="2" hl_lines="1"
        ch_samplesheet = channel.fromPath("./data/samplesheet.csv")
    ```

Esto combina la operación `splitCsv` (lectura del CSV con encabezados) y la operación `map` (estructuración de datos como tuplas `[meta, archivo]`) en un solo paso. Aplica ese cambio y ejecuta el pipeline:

```bash
nextflow run main.nf
```

??? success "Salida del comando"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `main.nf` [deadly_mercator] DSL2 - revision: bd6b0224e9

    [[id:patientA, repeat:1, type:normal], patientA_rep1_normal.bam]
    [[id:patientA, repeat:1, type:tumor], patientA_rep1_tumor.bam]
    [[id:patientA, repeat:2, type:normal], patientA_rep2_normal.bam]
    [[id:patientA, repeat:2, type:tumor], patientA_rep2_tumor.bam]
    [[id:patientB, repeat:1, type:normal], patientB_rep1_normal.bam]
    [[id:patientB, repeat:1, type:tumor], patientB_rep1_tumor.bam]
    [[id:patientC, repeat:1, type:normal], patientC_rep1_normal.bam]
    [[id:patientC, repeat:1, type:tumor], patientC_rep1_tumor.bam]
    ```

Ahora tenemos un canal donde cada elemento es una tupla `[meta, archivo]`: los metadatos separados de las rutas de archivos. Esta estructura nos permite dividir y agrupar nuestra carga de trabajo según los campos de metadatos.

---

## 2. Filtrar y transformar datos

### 2.1. Filtrar datos con `filter`

Podemos usar el [operador `filter`](https://www.nextflow.io/docs/latest/operator.html#filter) para filtrar los datos según una condición. Digamos que solo queremos procesar muestras normales. Podemos hacer esto filtrando los datos según el campo `type`. Insertemos esto antes del operador `view`.

=== "Después"

    ```groovy title="main.nf" linenums="2" hl_lines="6"
        ch_samples = channel.fromPath("./data/samplesheet.csv")
            .splitCsv(header: true)
            .map{ row ->
              [[id:row.id, repeat:row.repeat, type:row.type], row.bam]
            }
            .filter { meta, file -> meta.type == 'normal' }
            .view()
    ```

=== "Antes"

    ```groovy title="main.nf" linenums="2"
        ch_samples = channel.fromPath("./data/samplesheet.csv")
            .splitCsv(header: true)
            .map{ row ->
              [[id:row.id, repeat:row.repeat, type:row.type], row.bam]
            }
            .view()
    ```

Ejecuta el workflow nuevamente para ver el resultado filtrado:

```bash
nextflow run main.nf
```

??? success "Salida del comando"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `main.nf` [admiring_brown] DSL2 - revision: 194d61704d

    [[id:patientA, repeat:1, type:normal], patientA_rep1_normal.bam]
    [[id:patientA, repeat:2, type:normal], patientA_rep2_normal.bam]
    [[id:patientB, repeat:1, type:normal], patientB_rep1_normal.bam]
    [[id:patientC, repeat:1, type:normal], patientC_rep1_normal.bam]
    ```

Hemos filtrado exitosamente los datos para incluir solo muestras normales. Repasemos cómo funciona esto.

El operador `filter` toma un closure que se aplica a cada elemento del canal. Si el closure devuelve `true`, el elemento se incluye; si devuelve `false`, el elemento se excluye.

En nuestro caso, queremos conservar solo las muestras donde `meta.type == 'normal'`. El closure usa la tupla `meta,file` para referirse a cada muestra, accede al tipo de muestra con `meta.type` y verifica si es igual a `'normal'`.

Esto se logra con el único closure que introdujimos anteriormente:

```groovy title="main.nf" linenums="7"
    .filter { meta, file -> meta.type == 'normal' }
```

### 2.2. Crear canales filtrados separados

Actualmente estamos aplicando el filtro al canal creado directamente desde el CSV, pero queremos filtrar de más de una manera, así que reescribamos la lógica para crear un canal filtrado separado para las muestras normales:

=== "Después"

    ```groovy title="main.nf" linenums="2" hl_lines="6 8"
        ch_samples = channel.fromPath("./data/samplesheet.csv")
            .splitCsv(header: true)
            .map{ row ->
                [[id:row.id, repeat:row.repeat, type:row.type], row.bam]
            }
        ch_normal_samples = ch_samples
            .filter { meta, file -> meta.type == 'normal' }
        ch_normal_samples
            .view()
    ```

=== "Antes"

    ```groovy title="main.nf" linenums="2"
        ch_samples = channel.fromPath("./data/samplesheet.csv")
            .splitCsv(header: true)
            .map{ row ->
              [[id:row.id, repeat:row.repeat, type:row.type], row.bam]
            }
            .filter { meta, file -> meta.type == 'normal' }
            .view()
    ```

Ejecuta el pipeline para ver los resultados:

```bash
nextflow run main.nf
```

??? success "Salida del comando"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `main.nf` [trusting_poisson] DSL2 - revision: 639186ee74

    [[id:patientA, repeat:1, type:normal], patientA_rep1_normal.bam]
    [[id:patientA, repeat:2, type:normal], patientA_rep2_normal.bam]
    [[id:patientB, repeat:1, type:normal], patientB_rep1_normal.bam]
    [[id:patientC, repeat:1, type:normal], patientC_rep1_normal.bam]
    ```

Hemos filtrado exitosamente los datos y creado un canal separado para las muestras normales.

Creemos también un canal filtrado para las muestras tumorales:

=== "Después"

    ```groovy title="main.nf" linenums="7" hl_lines="3-8"
        ch_normal_samples = ch_samples
            .filter { meta, file -> meta.type == 'normal' }
        ch_tumor_samples = ch_samples
            .filter { meta, file -> meta.type == 'tumor' }
        ch_normal_samples
            .view{'Normal sample: ' + it}
        ch_tumor_samples
            .view{'Tumor sample: ' + it}
    ```

=== "Antes"

    ```groovy title="main.nf" linenums="7" hl_lines="3 4"
        ch_normal_samples = ch_samples
            .filter { meta, file -> meta.type == 'normal' }
        ch_normal_samples
            .view()
    ```

```bash
nextflow run main.nf
```

??? success "Salida del comando"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `main.nf` [maniac_boltzmann] DSL2 - revision: 3636b6576b

    Tumor sample: [[id:patientA, repeat:1, type:tumor], patientA_rep1_tumor.bam]
    Tumor sample: [[id:patientA, repeat:2, type:tumor], patientA_rep2_tumor.bam]
    Normal sample: [[id:patientA, repeat:1, type:normal], patientA_rep1_normal.bam]
    Normal sample: [[id:patientA, repeat:2, type:normal], patientA_rep2_normal.bam]
    Normal sample: [[id:patientB, repeat:1, type:normal], patientB_rep1_normal.bam]
    Normal sample: [[id:patientC, repeat:1, type:normal], patientC_rep1_normal.bam]
    Tumor sample: [[id:patientB, repeat:1, type:tumor], patientB_rep1_tumor.bam]
    Tumor sample: [[id:patientC, repeat:1, type:tumor], patientC_rep1_tumor.bam]
    ```

Hemos separado las muestras normales y tumorales en dos canales diferentes, y usamos un closure proporcionado a `view()` para etiquetarlas de manera diferente en la salida: `ch_tumor_samples.view{'Tumor sample: ' + it}`.

### Conclusión

En esta sección, has aprendido:

- **Filtrar datos**: Cómo filtrar datos con `filter`
- **Dividir datos**: Cómo dividir datos en diferentes canales según una condición
- **Visualizar datos**: Cómo usar `view` para imprimir los datos y etiquetar la salida de diferentes canales

Ahora hemos separado las muestras normales y tumorales en dos canales diferentes. A continuación, uniremos las muestras normales y tumorales por el campo `id`.

---

## 3. Unir canales por identificadores

En la sección anterior, separamos las muestras normales y tumorales en dos canales diferentes. Estas podrían procesarse de forma independiente usando procesos o workflows específicos según su tipo. Pero ¿qué sucede cuando queremos comparar las muestras normales y tumorales del mismo paciente? En este punto, necesitamos unirlas nuevamente asegurándonos de hacer coincidir las muestras según su campo `id`.

Nextflow incluye muchos métodos para combinar canales, pero en este caso el operador más apropiado es [`join`](https://www.nextflow.io/docs/latest/operator.html#join). Si estás familiarizado con SQL, actúa como la operación `JOIN`, donde especificamos la clave para unir y el tipo de unión a realizar.

### 3.1. Usar `map` y `join` para combinar por ID de paciente

Si revisamos la documentación de [`join`](https://www.nextflow.io/docs/latest/operator.html#join), podemos ver que por defecto une dos canales basándose en el primer elemento de cada tupla.

#### 3.1.1. Verificar la estructura de datos

Si ya no tienes disponible la salida de la consola, ejecutemos el pipeline para verificar nuestra estructura de datos y ver cómo necesitamos modificarla para unir por el campo `id`.

```bash
nextflow run main.nf
```

??? success "Salida del comando"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `main.nf` [maniac_boltzmann] DSL2 - revision: 3636b6576b

    Tumor sample: [[id:patientA, repeat:1, type:tumor], patientA_rep1_tumor.bam]
    Tumor sample: [[id:patientA, repeat:2, type:tumor], patientA_rep2_tumor.bam]
    Normal sample: [[id:patientA, repeat:1, type:normal], patientA_rep1_normal.bam]
    Normal sample: [[id:patientA, repeat:2, type:normal], patientA_rep2_normal.bam]
    Normal sample: [[id:patientB, repeat:1, type:normal], patientB_rep1_normal.bam]
    Normal sample: [[id:patientC, repeat:1, type:normal], patientC_rep1_normal.bam]
    Tumor sample: [[id:patientB, repeat:1, type:tumor], patientB_rep1_tumor.bam]
    Tumor sample: [[id:patientC, repeat:1, type:tumor], patientC_rep1_tumor.bam]
    ```

Podemos ver que el campo `id` es el primer elemento en cada mapa de metadatos. Para que `join` funcione, debemos aislar el campo `id` en cada tupla. Después de eso, podemos simplemente usar el operador `join` para combinar los dos canales.

#### 3.1.2. Aislar el campo `id`

Para aislar el campo `id`, podemos usar el [operador `map`](https://www.nextflow.io/docs/latest/operator.html#map) para crear una nueva tupla con el campo `id` como primer elemento.

=== "Después"

    ```groovy title="main.nf" linenums="7" hl_lines="3 6"
        ch_normal_samples = ch_samples
            .filter { meta, file -> meta.type == 'normal' }
            .map { meta, file -> [meta.id, meta, file] }
        ch_tumor_samples = ch_samples
            .filter { meta, file -> meta.type == 'tumor' }
            .map { meta, file -> [meta.id, meta, file] }
        ch_normal_samples
            .view{'Normal sample: ' + it}
        ch_tumor_samples
            .view{'Tumor sample: ' + it}
    ```

=== "Antes"

    ```groovy title="main.nf" linenums="7"
        ch_normal_samples = ch_samples
            .filter { meta, file -> meta.type == 'normal' }
        ch_tumor_samples = ch_samples
            .filter { meta, file -> meta.type == 'tumor' }
        ch_normal_samples
            .view{'Normal sample: ' + it}
        ch_tumor_samples
            .view{'Tumor sample: ' + it}
    ```

```bash
nextflow run main.nf
```

??? success "Salida del comando"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `main.nf` [mad_lagrange] DSL2 - revision: 9940b3f23d

    Tumor sample: [patientA, [id:patientA, repeat:1, type:tumor], patientA_rep1_tumor.bam]
    Tumor sample: [patientA, [id:patientA, repeat:2, type:tumor], patientA_rep2_tumor.bam]
    Normal sample: [patientA, [id:patientA, repeat:1, type:normal], patientA_rep1_normal.bam]
    Normal sample: [patientA, [id:patientA, repeat:2, type:normal], patientA_rep2_normal.bam]
    Tumor sample: [patientB, [id:patientB, repeat:1, type:tumor], patientB_rep1_tumor.bam]
    Tumor sample: [patientC, [id:patientC, repeat:1, type:tumor], patientC_rep1_tumor.bam]
    Normal sample: [patientB, [id:patientB, repeat:1, type:normal], patientB_rep1_normal.bam]
    Normal sample: [patientC, [id:patientC, repeat:1, type:normal], patientC_rep1_normal.bam]
    ```

Puede ser sutil, pero deberías poder ver que el primer elemento en cada tupla es el campo `id`.

#### 3.1.3. Combinar los dos canales

Ahora podemos usar el operador `join` para combinar los dos canales basándonos en el campo `id`.

Una vez más, usaremos `view` para imprimir las salidas unidas.

=== "Después"

    ```groovy title="main.nf" linenums="7" hl_lines="7-9"
        ch_normal_samples = ch_samples
            .filter { meta, file -> meta.type == 'normal' }
            .map { meta, file -> [meta.id, meta, file] }
        ch_tumor_samples = ch_samples
            .filter { meta, file -> meta.type == 'tumor' }
            .map { meta, file -> [meta.id, meta, file] }
        ch_joined_samples = ch_normal_samples
            .join(ch_tumor_samples)
        ch_joined_samples.view()
    ```

=== "Antes"

    ```groovy title="main.nf" linenums="7" hl_lines="7-10"
        ch_normal_samples = ch_samples
            .filter { meta, file -> meta.type == 'normal' }
            .map { meta, file -> [meta.id, meta, file] }
        ch_tumor_samples = ch_samples
            .filter { meta, file -> meta.type == 'tumor' }
            .map { meta, file -> [meta.id, meta, file] }
        ch_normal_samples
            .view{'Normal sample: ' + it}
        ch_tumor_samples
            .view{'Tumor sample: ' + it}
    ```

```bash
nextflow run main.nf
```

??? success "Salida del comando"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `main.nf` [soggy_wiles] DSL2 - revision: 3bc1979889

    [patientA, [id:patientA, repeat:1, type:normal], patientA_rep1_normal.bam, [id:patientA, repeat:1, type:tumor], patientA_rep1_tumor.bam]
    [patientA, [id:patientA, repeat:2, type:normal], patientA_rep2_normal.bam, [id:patientA, repeat:2, type:tumor], patientA_rep2_tumor.bam]
    [patientB, [id:patientB, repeat:1, type:normal], patientB_rep1_normal.bam, [id:patientB, repeat:1, type:tumor], patientB_rep1_tumor.bam]
    [patientC, [id:patientC, repeat:1, type:normal], patientC_rep1_normal.bam, [id:patientC, repeat:1, type:tumor], patientC_rep1_tumor.bam]
    ```

Es un poco difícil de leer porque es muy ancho, pero deberías poder ver que las muestras se han unido por el campo `id`. Cada tupla ahora tiene el formato:

- `id`: El ID de la muestra
- `normal_meta_map`: Los metadatos de la muestra normal incluyendo tipo, réplica y ruta al archivo BAM
- `normal_sample_file`: El archivo de la muestra normal
- `tumor_meta_map`: Los metadatos de la muestra tumoral incluyendo tipo, réplica y ruta al archivo BAM
- `tumor_sample`: La muestra tumoral incluyendo tipo, réplica y ruta al archivo BAM

!!! warning "Advertencia"

    El operador `join` descartará cualquier tupla sin coincidencia. En este ejemplo, nos aseguramos de que todas las muestras tuvieran coincidencia para tumor y normal, pero si esto no es así, debes usar el parámetro `remainder: true` para conservar las tuplas sin coincidencia. Consulta la [documentación](https://www.nextflow.io/docs/latest/operator.html#join) para más detalles.

Ahora sabes cómo usar `map` para aislar un campo en una tupla, y cómo usar `join` para combinar tuplas basándose en el primer campo.
Con este conocimiento, podemos combinar exitosamente canales basándonos en un campo compartido.

A continuación, consideraremos la situación en la que quieres unir por múltiples campos.

### 3.2. Unir por múltiples campos

Tenemos 2 réplicas para la muestra A, pero solo 1 para las muestras B y C. En este caso pudimos unirlas efectivamente usando el campo `id`, pero ¿qué pasaría si estuvieran desincronizadas? ¡Podríamos mezclar las muestras normales y tumorales de diferentes réplicas!

Para evitar esto, podemos unir por múltiples campos. En realidad hay múltiples formas de lograr esto, pero nos centraremos en crear una nueva clave de unión que incluya tanto el `id` de la muestra como el número de `replicate`.

Comencemos creando una nueva clave de unión. Podemos hacerlo de la misma manera que antes, usando el [operador `map`](https://www.nextflow.io/docs/latest/operator.html#map) para crear una nueva tupla con los campos `id` y `repeat` como primer elemento.

=== "Después"

    ```groovy title="main.nf" linenums="7" hl_lines="3 6"
        ch_normal_samples = ch_samples
            .filter { meta, file -> meta.type == 'normal' }
            .map { meta, file -> [[meta.id, meta.repeat], meta, file] }
        ch_tumor_samples = ch_samples
            .filter { meta, file -> meta.type == 'tumor' }
            .map { meta, file -> [[meta.id, meta.repeat], meta, file] }
    ```

=== "Antes"

    ```groovy title="main.nf" linenums="7" hl_lines="3 6"
        ch_normal_samples = ch_samples
            .filter { meta, file -> meta.type == 'normal' }
            .map { meta, file -> [meta.id, meta, file] }
        ch_tumor_samples = ch_samples
            .filter { meta, file -> meta.type == 'tumor' }
            .map { meta, file -> [meta.id, meta, file] }
    ```

Ahora deberíamos ver que la unión ocurre pero usando tanto los campos `id` como `repeat`. Ejecuta el workflow:

```bash
nextflow run main.nf
```

??? success "Salida del comando"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `main.nf` [prickly_wing] DSL2 - revision: 3bebf22dee

    [[patientA, 1], [id:patientA, repeat:1, type:normal], patientA_rep1_normal.bam, [id:patientA, repeat:1, type:tumor], patientA_rep1_tumor.bam]
    [[patientA, 2], [id:patientA, repeat:2, type:normal], patientA_rep2_normal.bam, [id:patientA, repeat:2, type:tumor], patientA_rep2_tumor.bam]
    [[patientB, 1], [id:patientB, repeat:1, type:normal], patientB_rep1_normal.bam, [id:patientB, repeat:1, type:tumor], patientB_rep1_tumor.bam]
    [[patientC, 1], [id:patientC, repeat:1, type:normal], patientC_rep1_normal.bam, [id:patientC, repeat:1, type:tumor], patientC_rep1_tumor.bam]
    ```

Observa cómo tenemos una tupla de dos elementos (campos `id` y `repeat`) como primer elemento de cada resultado unido. Esto demuestra cómo se pueden usar elementos complejos como clave de unión, lo que permite una coincidencia bastante intrincada entre muestras de las mismas condiciones.

Si quieres explorar más formas de unir por diferentes claves, consulta la [documentación del operador join](https://www.nextflow.io/docs/latest/operator.html#join) para opciones y ejemplos adicionales.

### 3.3. Usar `subMap` para crear una nueva clave de unión

El enfoque anterior pierde los nombres de campo de nuestra clave de unión: los campos `id` y `repeat` se convierten en solo una lista de valores. Para conservar los nombres de campo para acceso posterior, podemos usar el [método `subMap`](<https://docs.groovy-lang.org/latest/html/groovy-jdk/java/util/Map.html#subMap(java.util.Collection)>).

El método `subMap` extrae solo los pares clave-valor especificados de un mapa. Aquí extraeremos solo los campos `id` y `repeat` para crear nuestra clave de unión.

=== "Después"

    ```groovy title="main.nf" linenums="7" hl_lines="3 6"
        ch_normal_samples = ch_samples
            .filter { meta, file -> meta.type == 'normal' }
            .map { meta, file -> [meta.subMap(['id', 'repeat']), meta, file] }
        ch_tumor_samples = ch_samples
            .filter { meta, file -> meta.type == 'tumor' }
            .map { meta, file -> [meta.subMap(['id', 'repeat']), meta, file] }
    ```

=== "Antes"

    ```groovy title="main.nf" linenums="7" hl_lines="3 6"
        ch_normal_samples = ch_samples
            .filter { meta, file -> meta.type == 'normal' }
            .map { meta, file -> [[meta.id, meta.repeat], meta, file] }
        ch_tumor_samples = ch_samples
            .filter { meta, file -> meta.type == 'tumor' }
            .map { meta, file -> [[meta.id, meta.repeat], meta, file] }
    ```

```bash
nextflow run main.nf
```

??? success "Salida del comando"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `main.nf` [reverent_wing] DSL2 - revision: 847016c3b7

    [[id:patientA, repeat:1], [id:patientA, repeat:1, type:normal], patientA_rep1_normal.bam, [id:patientA, repeat:1, type:tumor], patientA_rep1_tumor.bam]
    [[id:patientA, repeat:2], [id:patientA, repeat:2, type:normal], patientA_rep2_normal.bam, [id:patientA, repeat:2, type:tumor], patientA_rep2_tumor.bam]
    [[id:patientB, repeat:1], [id:patientB, repeat:1, type:normal], patientB_rep1_normal.bam, [id:patientB, repeat:1, type:tumor], patientB_rep1_tumor.bam]
    [[id:patientC, repeat:1], [id:patientC, repeat:1, type:normal], patientC_rep1_normal.bam, [id:patientC, repeat:1, type:tumor], patientC_rep1_tumor.bam]
    ```

Ahora tenemos una nueva clave de unión que no solo incluye los campos `id` y `repeat`, sino que también conserva los nombres de campo para que podamos acceder a ellos más tarde por nombre, por ejemplo `meta.id` y `meta.repeat`.

### 3.4. Usar un closure con nombre en map

Para evitar la duplicación y reducir errores, podemos usar un closure con nombre. Un closure con nombre nos permite crear una función reutilizable que podemos llamar en múltiples lugares.

Para hacerlo, primero definimos el closure como una nueva variable:

=== "Después"

    ```groovy title="main.nf" linenums="2" hl_lines="7"
        ch_samples = channel.fromPath("./data/samplesheet.csv")
            .splitCsv(header: true)
            .map{ row ->
                [[id:row.id, repeat:row.repeat, type:row.type], row.bam]
            }

        getSampleIdAndReplicate = { meta, bam -> [ meta.subMap(['id', 'repeat']), meta, file(bam) ] }

        ch_normal_samples = ch_samples
            .filter { meta, file -> meta.type == 'normal' }
    ```

=== "Antes"

    ```groovy title="main.nf" linenums="2"
        ch_samples = channel.fromPath("./data/samplesheet.csv")
            .splitCsv(header: true)
            .map{ row ->
                [[id:row.id, repeat:row.repeat, type:row.type], row.bam]
            }
        ch_normal_samples = ch_samples
            .filter { meta, file -> meta.type == 'normal' }
    ```

Hemos definido la transformación de map como una variable con nombre que podemos reutilizar.

Ten en cuenta que también convertimos la ruta del archivo a un objeto Path usando `file()` para que cualquier proceso que reciba este canal pueda manejar el archivo correctamente (para más información consulta [Working with files](../working_with_files/)).

Implementemos el closure en nuestro workflow:

=== "Después"

    ```groovy title="main.nf" linenums="10" hl_lines="3 6"
        ch_normal_samples = ch_samples
            .filter { meta, file -> meta.type == 'normal' }
             .map ( getSampleIdAndReplicate )
        ch_tumor_samples = ch_samples
            .filter { meta, file -> meta.type == 'tumor' }
             .map ( getSampleIdAndReplicate )

    ```

=== "Antes"

    ```groovy title="main.nf" linenums="10" hl_lines="3 6"
        ch_normal_samples = ch_samples
            .filter { meta, file -> meta.type == 'normal' }
            .map { meta, file -> [meta.subMap(['id', 'repeat']), meta, file] }
        ch_tumor_samples = ch_samples
            .filter { meta, file -> meta.type == 'tumor' }
            .map { meta, file -> [meta.subMap(['id', 'repeat']), meta, file] }
    ```

!!! note "Nota"

    El operador `map` ha cambiado de usar `{ }` a usar `( )` para pasar el closure como argumento. Esto se debe a que el operador `map` espera un closure como argumento y `{ }` se usa para definir un closure anónimo. Al llamar a un closure con nombre, usa la sintaxis `( )`.

Ejecuta el workflow una vez más para verificar que todo sigue funcionando:

```bash
nextflow run main.nf
```

??? success "Salida del comando"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `main.nf` [angry_meninsky] DSL2 - revision: 2edc226b1d

    [[id:patientA, repeat:1], [id:patientA, repeat:1, type:normal], patientA_rep1_normal.bam, [id:patientA, repeat:1, type:tumor], patientA_rep1_tumor.bam]
    [[id:patientA, repeat:2], [id:patientA, repeat:2, type:normal], patientA_rep2_normal.bam, [id:patientA, repeat:2, type:tumor], patientA_rep2_tumor.bam]
    [[id:patientB, repeat:1], [id:patientB, repeat:1, type:normal], patientB_rep1_normal.bam, [id:patientB, repeat:1, type:tumor], patientB_rep1_tumor.bam]
    [[id:patientC, repeat:1], [id:patientC, repeat:1, type:normal], patientC_rep1_normal.bam, [id:patientC, repeat:1, type:tumor], patientC_rep1_tumor.bam]
    ```

Usar un closure con nombre nos permite reutilizar la misma transformación en múltiples lugares, reduciendo el riesgo de errores y haciendo el código más legible y fácil de mantener.

### 3.5. Reducir la duplicación de datos

Tenemos muchos datos duplicados en nuestro workflow. Cada elemento en las muestras unidas repite los campos `id` y `repeat`. Dado que esta información ya está disponible en la clave de agrupación, podemos evitar esta redundancia. Como recordatorio, nuestra estructura de datos actual se ve así:

```groovy
[
  [
    "id": "sampleC",
    "repeat": "1",
  ],
  [
    "id": "sampleC",
    "repeat": "1",
    "type": "normal",
  ],
  "sampleC_rep1_normal.bam"
  [
    "id": "sampleC",
    "repeat": "1",
    "type": "tumor",
  ],
  "sampleC_rep1_tumor.bam"
]
```

Dado que los campos `id` y `repeat` están disponibles en la clave de agrupación, eliminémoslos del resto de cada elemento del canal para evitar la duplicación. Podemos hacer esto usando el método `subMap` para crear un nuevo mapa con solo el campo `type`. Este enfoque nos permite mantener toda la información necesaria mientras eliminamos la redundancia en nuestra estructura de datos.

=== "Después"

    ```groovy title="main.nf" linenums="8" hl_lines="1"
        getSampleIdAndReplicate = { meta, bam -> [ meta.subMap(['id', 'repeat']), meta.subMap(['type']), file(bam) ] }
    ```

=== "Antes"

    ```groovy title="main.nf" linenums="8" hl_lines="1"
        getSampleIdAndReplicate = { meta, bam -> [ meta.subMap(['id', 'repeat']), meta, file(bam) ] }
    ```

Ahora el closure devuelve una tupla donde el primer elemento contiene los campos `id` y `repeat`, y el segundo elemento contiene solo el campo `type`. Esto elimina la redundancia almacenando la información de `id` y `repeat` una sola vez en la clave de agrupación, mientras se mantiene toda la información necesaria.

Ejecuta el workflow para ver cómo se ve esto:

```bash
nextflow run main.nf
```

??? success "Salida del comando"

    ```console
    [[id:patientA, repeat:1], [type:normal], /workspaces/training/side-quests/splitting_and_grouping/patientA_rep1_normal.bam, [type:tumor], /workspaces/training/side-quests/splitting_and_grouping/patientA_rep1_tumor.bam]
    [[id:patientA, repeat:2], [type:normal], /workspaces/training/side-quests/splitting_and_grouping/patientA_rep2_normal.bam, [type:tumor], /workspaces/training/side-quests/splitting_and_grouping/patientA_rep2_tumor.bam]
    [[id:patientB, repeat:1], [type:normal], /workspaces/training/side-quests/splitting_and_grouping/patientB_rep1_normal.bam, [type:tumor], /workspaces/training/side-quests/splitting_and_grouping/patientB_rep1_tumor.bam]
    [[id:patientC, repeat:1], [type:normal], /workspaces/training/side-quests/splitting_and_grouping/patientC_rep1_normal.bam, [type:tumor], /workspaces/training/side-quests/splitting_and_grouping/patientC_rep1_tumor.bam]
    ```

Podemos ver que solo indicamos los campos `id` y `repeat` una vez en la clave de agrupación y tenemos el campo `type` en los datos de la muestra. No hemos perdido ninguna información pero logramos hacer el contenido de nuestro canal más conciso.

### 3.6. Eliminar información redundante

Eliminamos la información duplicada anteriormente, pero todavía tenemos otra información redundante en nuestros canales.

Al principio, separamos las muestras normales y tumorales usando `filter`, luego las unimos basándonos en las claves `id` y `repeat`. El operador `join` preserva el orden en que se fusionan las tuplas, por lo que en nuestro caso, con las muestras normales en el lado izquierdo y las tumorales en el derecho, el canal resultante mantiene esta estructura: `id, <elementos normales>, <elementos tumorales>`.

Dado que conocemos la posición de cada elemento en nuestro canal, podemos simplificar aún más la estructura eliminando los metadatos `[type:normal]` y `[type:tumor]`.

=== "Después"

    ```groovy title="main.nf" linenums="8" hl_lines="1"
        getSampleIdAndReplicate = { meta, file -> [ meta.subMap(['id', 'repeat']), file ] }
    ```

=== "Antes"

    ```groovy title="main.nf" linenums="8" hl_lines="1"
        getSampleIdAndReplicate = { meta, file -> [ meta.subMap(['id', 'repeat']), meta.subMap(['type']), file ] }
    ```

Ejecuta nuevamente para ver el resultado:

```bash
nextflow run main.nf
```

??? success "Salida del comando"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `main.nf` [confident_leavitt] DSL2 - revision: a2303895bd

    [[id:patientA, repeat:1], patientA_rep1_normal.bam, patientA_rep1_tumor.bam]
    [[id:patientA, repeat:2], patientA_rep2_normal.bam, patientA_rep2_tumor.bam]
    [[id:patientB, repeat:1], patientB_rep1_normal.bam, patientB_rep1_tumor.bam]
    [[id:patientC, repeat:1], patientC_rep1_normal.bam, patientC_rep1_tumor.bam]
    ```

### Conclusión

En esta sección, has aprendido:

- **Manipular tuplas**: Cómo usar `map` para aislar un campo en una tupla
- **Unir tuplas**: Cómo usar `join` para combinar tuplas basándose en el primer campo
- **Crear claves de unión**: Cómo usar `subMap` para crear una nueva clave de unión
- **Closures con nombre**: Cómo usar un closure con nombre en map
- **Unión por múltiples campos**: Cómo unir por múltiples campos para una coincidencia más precisa
- **Optimización de la estructura de datos**: Cómo simplificar la estructura del canal eliminando información redundante

Ahora tienes un workflow que puede dividir una hoja de muestras, filtrar las muestras normales y tumorales, unirlas por ID de muestra y número de réplica, y luego imprimir los resultados.

Este es un patrón común en workflows de bioinformática donde necesitas hacer coincidir muestras u otros tipos de datos después de procesarlos de forma independiente, por lo que es una habilidad muy útil. A continuación, veremos cómo repetir una muestra múltiples veces.

## 4. Distribuir muestras en intervalos

Un patrón clave en los workflows de bioinformática es distribuir el análisis a través de regiones genómicas. Por ejemplo, la llamada de variantes puede paralelizarse dividiendo el genoma en intervalos (como cromosomas o regiones más pequeñas). Esta estrategia de paralelización mejora significativamente la eficiencia del pipeline al distribuir la carga computacional entre múltiples núcleos o nodos, reduciendo el tiempo total de ejecución.

En la siguiente sección, demostraremos cómo distribuir nuestros datos de muestras a través de múltiples intervalos genómicos. Emparejaremos cada muestra con cada intervalo, permitiendo el procesamiento paralelo de diferentes regiones genómicas. Esto multiplicará el tamaño de nuestro conjunto de datos por el número de intervalos, creando múltiples unidades de análisis independientes que pueden reunirse más adelante.

### 4.1. Distribuir muestras en intervalos usando `combine`

Comencemos creando un canal de intervalos. Para simplificar, usaremos solo 3 intervalos que definiremos manualmente. En un workflow real, podrías leerlos desde un archivo de entrada o incluso crear un canal con muchos archivos de intervalos.

=== "Después"

    ```groovy title="main.nf" linenums="17" hl_lines="2"
            .join(ch_tumor_samples)
        ch_intervals = channel.of('chr1', 'chr2', 'chr3')
    ```

=== "Antes"

    ```groovy title="main.nf" linenums="17" hl_lines="2"
            .join(ch_tumor_samples)
        ch_joined_samples.view()
    ```

Recuerda que queremos repetir cada muestra para cada intervalo. Esto a veces se denomina el producto cartesiano de las muestras e intervalos. Podemos lograrlo usando el [operador `combine`](https://www.nextflow.io/docs/latest/operator.html#combine). Este tomará cada elemento del canal 1 y lo repetirá para cada elemento del canal 2. Agreguemos un operador combine a nuestro workflow:

=== "Después"

    ```groovy title="main.nf" linenums="18" hl_lines="3-5"
        ch_intervals = channel.of('chr1', 'chr2', 'chr3')

        ch_combined_samples = ch_joined_samples
            .combine(ch_intervals)
            .view()
    ```

=== "Antes"

    ```groovy title="main.nf" linenums="18"
        ch_intervals = channel.of('chr1', 'chr2', 'chr3')
    ```

Ahora ejecutémoslo y veamos qué sucede:

```bash
nextflow run main.nf
```

??? success "Salida del comando"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `main.nf` [mighty_tesla] DSL2 - revision: ae013ab70b

    [[id:patientA, repeat:1], patientA_rep1_normal.bam, patientA_rep1_tumor.bam, chr1]
    [[id:patientA, repeat:1], patientA_rep1_normal.bam, patientA_rep1_tumor.bam, chr2]
    [[id:patientA, repeat:1], patientA_rep1_normal.bam, patientA_rep1_tumor.bam, chr3]
    [[id:patientA, repeat:2], patientA_rep2_normal.bam, patientA_rep2_tumor.bam, chr1]
    [[id:patientA, repeat:2], patientA_rep2_normal.bam, patientA_rep2_tumor.bam, chr2]
    [[id:patientA, repeat:2], patientA_rep2_normal.bam, patientA_rep2_tumor.bam, chr3]
    [[id:patientB, repeat:1], patientB_rep1_normal.bam, patientB_rep1_tumor.bam, chr1]
    [[id:patientB, repeat:1], patientB_rep1_normal.bam, patientB_rep1_tumor.bam, chr2]
    [[id:patientB, repeat:1], patientB_rep1_normal.bam, patientB_rep1_tumor.bam, chr3]
    [[id:patientC, repeat:1], patientC_rep1_normal.bam, patientC_rep1_tumor.bam, chr1]
    [[id:patientC, repeat:1], patientC_rep1_normal.bam, patientC_rep1_tumor.bam, chr2]
    [[id:patientC, repeat:1], patientC_rep1_normal.bam, patientC_rep1_tumor.bam, chr3]
    ```

¡Éxito! Hemos repetido cada muestra para cada intervalo en nuestra lista de 3 intervalos. Hemos triplicado efectivamente el número de elementos en nuestro canal.

Es un poco difícil de leer, así que en la siguiente sección lo organizaremos mejor.

### 4.2. Organizar el canal

Podemos usar el operador `map` para ordenar y refactorizar nuestros datos de muestras para que sean más fáciles de entender. Movamos la cadena de intervalos al mapa de agrupación en el primer elemento.

=== "Después"

    ```groovy title="main.nf" linenums="20" hl_lines="3-9"
        ch_combined_samples = ch_joined_samples
            .combine(ch_intervals)
            .map { grouping_key, normal, tumor, interval ->
                [
                    grouping_key + [interval: interval],
                    normal,
                    tumor
                ]
            }
            .view()
    ```

=== "Antes"

    ```groovy title="main.nf" linenums="20"
        ch_combined_samples = ch_joined_samples
            .combine(ch_intervals)
            .view()
    ```

Analicemos paso a paso lo que hace esta operación map.

Primero, usamos parámetros con nombre para hacer el código más legible. Al usar los nombres `grouping_key`, `normal`, `tumor` e `interval`, podemos referirnos a los elementos de la tupla por nombre en lugar de por índice:

```groovy
        .map { grouping_key, normal, tumor, interval ->
```

A continuación, combinamos el `grouping_key` con el campo `interval`. El `grouping_key` es un mapa que contiene los campos `id` y `repeat`. Creamos un nuevo mapa con el `interval` y los fusionamos usando la adición de mapas de Groovy (`+`):

```groovy
                grouping_key + [interval: interval],
```

Finalmente, devolvemos esto como una tupla con tres elementos: el mapa de metadatos combinado, el archivo de la muestra normal y el archivo de la muestra tumoral:

```groovy
            [
                grouping_key + [interval: interval],
                normal,
                tumor
            ]
```

Ejecutémoslo nuevamente y verifiquemos el contenido del canal:

```bash
nextflow run main.nf
```

??? success "Salida del comando"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `main.nf` [sad_hawking] DSL2 - revision: 1f6f6250cd

    [[id:patientA, repeat:1, interval:chr1], patientA_rep1_normal.bam, patientA_rep1_tumor.bam]
    [[id:patientA, repeat:1, interval:chr2], patientA_rep1_normal.bam, patientA_rep1_tumor.bam]
    [[id:patientA, repeat:1, interval:chr3], patientA_rep1_normal.bam, patientA_rep1_tumor.bam]
    [[id:patientA, repeat:2, interval:chr1], patientA_rep2_normal.bam, patientA_rep2_tumor.bam]
    [[id:patientA, repeat:2, interval:chr2], patientA_rep2_normal.bam, patientA_rep2_tumor.bam]
    [[id:patientA, repeat:2, interval:chr3], patientA_rep2_normal.bam, patientA_rep2_tumor.bam]
    [[id:patientB, repeat:1, interval:chr1], patientB_rep1_normal.bam, patientB_rep1_tumor.bam]
    [[id:patientB, repeat:1, interval:chr2], patientB_rep1_normal.bam, patientB_rep1_tumor.bam]
    [[id:patientB, repeat:1, interval:chr3], patientB_rep1_normal.bam, patientB_rep1_tumor.bam]
    [[id:patientC, repeat:1, interval:chr1], patientC_rep1_normal.bam, patientC_rep1_tumor.bam]
    [[id:patientC, repeat:1, interval:chr2], patientC_rep1_normal.bam, patientC_rep1_tumor.bam]
    [[id:patientC, repeat:1, interval:chr3], patientC_rep1_normal.bam, patientC_rep1_tumor.bam]
    ```

Usar `map` para transformar tus datos en la estructura correcta puede ser complicado, pero es fundamental para una manipulación efectiva de los datos.

Ahora tenemos cada muestra repetida a través de todos los intervalos genómicos, creando múltiples unidades de análisis independientes que pueden procesarse en paralelo. Pero ¿qué pasa si queremos reunir muestras relacionadas? En la siguiente sección, aprenderemos cómo agrupar muestras que comparten atributos comunes.

### Conclusión

En esta sección, has aprendido:

- **Distribuir muestras en intervalos**: Cómo usar `combine` para repetir muestras en intervalos
- **Crear productos cartesianos**: Cómo generar todas las combinaciones de muestras e intervalos
- **Organizar la estructura del canal**: Cómo usar `map` para reestructurar datos para mejor legibilidad
- **Preparación para procesamiento paralelo**: Cómo configurar datos para análisis distribuido

## 5. Agregar muestras usando `groupTuple`

En las secciones anteriores, aprendimos cómo dividir datos de un archivo de entrada y filtrar por campos específicos (en nuestro caso muestras normales y tumorales). Pero esto solo cubre un tipo de unión. ¿Qué pasa si queremos agrupar muestras por un atributo específico? Por ejemplo, en lugar de unir pares normal-tumor coincidentes, podríamos querer procesar todas las muestras de "sampleA" juntas independientemente de su tipo. Este patrón es común en workflows de bioinformática donde puedes querer procesar muestras relacionadas por separado por razones de eficiencia antes de comparar o combinar los resultados al final.

Nextflow incluye métodos integrados para hacer esto, el principal que veremos es `groupTuple`.

Comencemos agrupando todas nuestras muestras que tienen los mismos campos `id` e `interval`; esto sería típico de un análisis donde queremos agrupar réplicas técnicas pero mantener separadas las muestras significativamente diferentes.

Para hacer esto, debemos separar nuestras variables de agrupación para poder usarlas de forma aislada.

El primer paso es similar a lo que hicimos en la sección anterior. Debemos aislar nuestra variable de agrupación como el primer elemento de la tupla. Recuerda que nuestro primer elemento es actualmente un mapa de los campos `id`, `repeat` e `interval`:

```groovy title="main.nf" linenums="1"
{
  "id": "sampleA",
  "repeat": "1",
  "interval": "chr1"
}
```

Podemos reutilizar el método `subMap` de antes para aislar nuestros campos `id` e `interval` del mapa. Como antes, usaremos el operador `map` para aplicar el método `subMap` al primer elemento de la tupla para cada muestra.

=== "Después"

    ```groovy title="main.nf" linenums="20" hl_lines="11-19"
        ch_combined_samples = ch_joined_samples
            .combine(ch_intervals)
            .map { grouping_key, normal, tumor, interval ->
                [
                    grouping_key + [interval: interval],
                    normal,
                    tumor
                ]
            }

        ch_grouped_samples = ch_combined_samples
            .map { grouping_key, normal, tumor ->
                [
                    grouping_key.subMap('id', 'interval'),
                    normal,
                    tumor
                ]
              }
              .view()
    ```

=== "Antes"

    ```groovy title="main.nf" linenums="20" hl_lines="10"
        ch_combined_samples = ch_joined_samples
            .combine(ch_intervals)
            .map { grouping_key, normal, tumor, interval ->
                [
                    grouping_key + [interval: interval],
                    normal,
                    tumor
                ]
            }
            .view()
    ```

Ejecutémoslo nuevamente y verifiquemos el contenido del canal:

```bash
nextflow run main.nf
```

??? success "Salida del comando"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `main.nf` [hopeful_brenner] DSL2 - revision: 7f4f7fea76

    [[id:patientA, interval:chr1], patientA_rep1_normal.bam, patientA_rep1_tumor.bam]
    [[id:patientA, interval:chr2], patientA_rep1_normal.bam, patientA_rep1_tumor.bam]
    [[id:patientA, interval:chr3], patientA_rep1_normal.bam, patientA_rep1_tumor.bam]
    [[id:patientA, interval:chr1], patientA_rep2_normal.bam, patientA_rep2_tumor.bam]
    [[id:patientA, interval:chr2], patientA_rep2_normal.bam, patientA_rep2_tumor.bam]
    [[id:patientA, interval:chr3], patientA_rep2_normal.bam, patientA_rep2_tumor.bam]
    [[id:patientB, interval:chr1], patientB_rep1_normal.bam, patientB_rep1_tumor.bam]
    [[id:patientB, interval:chr2], patientB_rep1_normal.bam, patientB_rep1_tumor.bam]
    [[id:patientB, interval:chr3], patientB_rep1_normal.bam, patientB_rep1_tumor.bam]
    [[id:patientC, interval:chr1], patientC_rep1_normal.bam, patientC_rep1_tumor.bam]
    [[id:patientC, interval:chr2], patientC_rep1_normal.bam, patientC_rep1_tumor.bam]
    [[id:patientC, interval:chr3], patientC_rep1_normal.bam, patientC_rep1_tumor.bam]
    ```

Podemos ver que hemos aislado exitosamente los campos `id` e `interval`, pero aún no hemos agrupado las muestras.

!!! note "Nota"

    Estamos descartando el campo `replicate` aquí. Esto se debe a que no lo necesitamos para el procesamiento posterior. Después de completar este tutorial, ¡intenta incluirlo sin afectar la agrupación posterior!

Ahora agrupemos las muestras por este nuevo elemento de agrupación, usando el [operador `groupTuple`](https://www.nextflow.io/docs/latest/operator.html#grouptuple).

=== "Después"

    ```groovy title="main.nf" linenums="30" hl_lines="9"
        ch_grouped_samples = ch_combined_samples
            .map { grouping_key, normal, tumor ->
                [
                    grouping_key.subMap('id', 'interval'),
                    normal,
                    tumor
                ]
              }
              .groupTuple()
              .view()
    ```

=== "Antes"

    ```groovy title="main.nf" linenums="30"
        ch_grouped_samples = ch_combined_samples
            .map { grouping_key, normal, tumor ->
                [
                    grouping_key.subMap('id', 'interval'),
                    normal,
                    tumor
                ]
              }
              .view()
    ```

¡Eso es todo! Solo agregamos una línea de código. Veamos qué sucede cuando lo ejecutamos:

```bash
nextflow run main.nf
```

??? success "Salida del comando"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `main.nf` [friendly_jang] DSL2 - revision: a1bee1c55d

    [[id:patientA, interval:chr1], [patientA_rep1_normal.bam, patientA_rep2_normal.bam], [patientA_rep1_tumor.bam, patientA_rep2_tumor.bam]]
    [[id:patientA, interval:chr2], [patientA_rep1_normal.bam, patientA_rep2_normal.bam], [patientA_rep1_tumor.bam, patientA_rep2_tumor.bam]]
    [[id:patientA, interval:chr3], [patientA_rep1_normal.bam, patientA_rep2_normal.bam], [patientA_rep1_tumor.bam, patientA_rep2_tumor.bam]]
    [[id:patientB, interval:chr1], [patientB_rep1_normal.bam], [patientB_rep1_tumor.bam]]
    [[id:patientB, interval:chr2], [patientB_rep1_normal.bam], [patientB_rep1_tumor.bam]]
    [[id:patientB, interval:chr3], [patientB_rep1_normal.bam], [patientB_rep1_tumor.bam]]
    [[id:patientC, interval:chr1], [patientC_rep1_normal.bam], [patientC_rep1_tumor.bam]]
    [[id:patientC, interval:chr2], [patientC_rep1_normal.bam], [patientC_rep1_tumor.bam]]
    [[id:patientC, interval:chr3], [patientC_rep1_normal.bam], [patientC_rep1_tumor.bam]]
    ```

Observa que nuestros datos han cambiado de estructura y dentro de cada elemento del canal los archivos ahora están contenidos en tuplas como `[patientA_rep1_normal.bam, patientA_rep2_normal.bam]`. Esto se debe a que cuando usamos `groupTuple`, Nextflow combina los archivos individuales para cada muestra de un grupo. Esto es importante recordarlo al intentar manejar los datos posteriormente.

!!! note "Nota"

    [`transpose`](https://www.nextflow.io/docs/latest/reference/operator.html#transpose) es lo opuesto de groupTuple. Desempaqueta los elementos en un canal y los aplana. ¡Intenta agregar `transpose` y deshacer la agrupación que realizamos anteriormente!

### Conclusión

En esta sección, has aprendido:

- **Agrupar muestras relacionadas**: Cómo usar `groupTuple` para agregar muestras por atributos comunes
- **Aislar claves de agrupación**: Cómo usar `subMap` para extraer campos específicos para la agrupación
- **Manejar estructuras de datos agrupadas**: Cómo trabajar con la estructura anidada creada por `groupTuple`
- **Manejo de réplicas técnicas**: Cómo agrupar muestras que comparten las mismas condiciones experimentales

---

## Resumen

En esta misión secundaria, has aprendido cómo dividir y agrupar datos usando canales.

Al modificar los datos a medida que fluyen a través del pipeline, puedes construir un pipeline escalable sin usar bucles o sentencias while, lo que ofrece varias ventajas sobre los enfoques más tradicionales:

- Podemos escalar a tantas o tan pocas entradas como queramos sin código adicional
- Nos enfocamos en manejar el flujo de datos a través del pipeline, en lugar de la iteración
- Podemos ser tan complejos o simples como sea necesario
- El pipeline se vuelve más declarativo, enfocándose en qué debe suceder en lugar de cómo debe suceder
- Nextflow optimizará la ejecución por nosotros ejecutando operaciones independientes en paralelo

Dominar estas operaciones de canal te permitirá construir pipelines flexibles y escalables que manejen relaciones de datos complejas sin recurrir a bucles o programación iterativa, permitiendo que Nextflow optimice la ejecución y paralelice las operaciones independientes automáticamente.

### Patrones clave

1.  **Crear datos de entrada estructurados:** Comenzando desde un archivo CSV con mapas de metadatos (basándose en patrones de [Metadata in workflows](../metadata/))

    ```groovy
    ch_samples = channel.fromPath("./data/samplesheet.csv")
        .splitCsv(header: true)
        .map{ row ->
          [[id:row.id, repeat:row.repeat, type:row.type], row.bam]
        }
    ```

2.  **Dividir datos en canales separados:** Usamos `filter` para dividir datos en flujos independientes basados en el campo `type`

    ```groovy
    channel.filter { it.type == 'tumor' }
    ```

3.  **Unir muestras coincidentes:** Usamos `join` para recombinar muestras relacionadas basándonos en los campos `id` y `repeat`

    - Unir dos canales por clave (primer elemento de la tupla)

    ```groovy
    tumor_ch.join(normal_ch)
    ```

    - Extraer clave de unión y unir por este valor

    ```groovy
    tumor_ch.map { meta, file -> [meta.id, meta, file] }
        .join(
          normal_ch.map { meta, file -> [meta.id, meta, file] }
        )
    ```

    - Unir por múltiples campos usando subMap

    ```groovy
    tumor_ch.map { meta, file -> [meta.subMap(['id', 'repeat']), meta, file] }
        .join(
          normal_ch.map { meta, file -> [meta.subMap(['id', 'repeat']), meta, file] }
        )
    ```

4.  **Distribuir a través de intervalos:** Usamos `combine` para crear productos cartesianos de muestras con intervalos genómicos para procesamiento paralelo.

    ```groovy
    samples_ch.combine(intervals_ch)
    ```

5.  **Agregar por claves de agrupación:** Usamos `groupTuple` para agrupar por el primer elemento en cada tupla, recopilando así las muestras que comparten los campos `id` e `interval` y fusionando las réplicas técnicas.

    ```groovy
    channel.groupTuple()
    ```

6.  **Optimizar la estructura de datos:** Usamos `subMap` para extraer campos específicos y creamos un closure con nombre para hacer las transformaciones reutilizables.

    - Extraer campos específicos de un mapa

    ```groovy
    meta.subMap(['id', 'repeat'])
    ```

    - Usar closure con nombre para transformaciones reutilizables

    ```groovy
    getSampleIdAndReplicate = { meta, file -> [meta.subMap(['id', 'repeat']), file] }
    channel.map(getSampleIdAndReplicate)
    ```

### Recursos adicionales

- [filter](https://www.nextflow.io/docs/latest/operator.html#filter)
- [map](https://www.nextflow.io/docs/latest/operator.html#map)
- [join](https://www.nextflow.io/docs/latest/operator.html#join)
- [groupTuple](https://www.nextflow.io/docs/latest/operator.html#grouptuple)
- [combine](https://www.nextflow.io/docs/latest/operator.html#combine)

---

## ¿Qué sigue?

Regresa al [menú de misiones secundarias](../) o haz clic en el botón en la parte inferior derecha de la página para continuar con el siguiente tema de la lista.
