# Procesamiento de archivos de entrada

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Traducción asistida por IA - [más información y sugerencias](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

Los workflows de análisis científico frecuentemente implican el procesamiento de grandes cantidades de archivos.
Nextflow proporciona herramientas poderosas para manejar archivos de manera eficiente, ayudándole a organizar y procesar sus datos con un mínimo de código.

### Objetivos de aprendizaje

En esta misión secundaria, exploraremos cómo Nextflow maneja los archivos, desde operaciones básicas hasta técnicas más avanzadas para trabajar con colecciones de archivos.
Aprenderá a extraer metadatos de los nombres de archivo, lo cual es un requisito común en los pipelines de análisis científico.

Al finalizar esta misión secundaria, podrá:

- Crear objetos Path a partir de cadenas de rutas de archivo usando el método `file()` de Nextflow
- Acceder a atributos de archivos como nombre, extensión y directorio padre
- Manejar archivos locales y remotos de forma transparente usando URIs
- Usar canales para automatizar el manejo de archivos con `channel.fromPath()` y `channel.fromFilePairs()`
- Extraer y estructurar metadatos de nombres de archivo usando manipulación de cadenas
- Agrupar archivos relacionados usando coincidencia de patrones y expresiones glob
- Integrar operaciones de archivos en procesos de Nextflow con el manejo adecuado de entradas
- Organizar las salidas de los procesos usando estructuras de directorios basadas en metadatos

Estas habilidades le ayudarán a construir workflows que puedan manejar diferentes tipos de archivos de entrada con gran flexibilidad.

### Requisitos previos

Antes de comenzar esta misión secundaria, debe:

- Haber completado el tutorial [Hello Nextflow](../../hello_nextflow/) o un curso equivalente para principiantes.
- Estar familiarizado con los conceptos y mecanismos básicos de Nextflow (procesos, canales, operadores)

<!-- I removed the suggestion to do the metamaps SQ first because that works more naturally after -->

---

## 0. Primeros pasos

#### Abrir el codespace de capacitación

Si aún no lo ha hecho, asegúrese de abrir el entorno de capacitación como se describe en la [Configuración del entorno](../envsetup/index.md).

[![Open in GitHub Codespaces](https://github.com/codespaces/badge.svg)](https://codespaces.new/nextflow-io/training?quickstart=1&ref=master)

#### Moverse al directorio del proyecto

Vamos a movernos al directorio donde se encuentran los archivos de este tutorial.

```bash
cd side-quests/working_with_files
```

Puede configurar VSCode para que se enfoque en este directorio:

```bash
code .
```

#### Revisar los materiales

Encontrará un archivo de workflow simple llamado `main.nf`, un directorio `modules` que contiene dos archivos de módulo, y un directorio `data` que contiene algunos archivos de datos de ejemplo.

??? abstract "Contenido del directorio"

    ```console
    .
    ├── data
    │   ├── patientA_rep1_normal_R1_001.fastq.gz
    │   ├── patientA_rep1_normal_R2_001.fastq.gz
    │   ├── patientA_rep1_tumor_R1_001.fastq.gz
    │   ├── patientA_rep1_tumor_R2_001.fastq.gz
    │   ├── patientA_rep2_normal_R1_001.fastq.gz
    │   ├── patientA_rep2_normal_R2_001.fastq.gz
    │   ├── patientA_rep2_tumor_R1_001.fastq.gz
    │   ├── patientA_rep2_tumor_R2_001.fastq.gz
    │   ├── patientB_rep1_normal_R1_001.fastq.gz
    │   ├── patientB_rep1_normal_R2_001.fastq.gz
    │   ├── patientB_rep1_tumor_R1_001.fastq.gz
    │   ├── patientB_rep1_tumor_R2_001.fastq.gz
    │   ├── patientC_rep1_normal_R1_001.fastq.gz
    │   ├── patientC_rep1_normal_R2_001.fastq.gz
    │   ├── patientC_rep1_tumor_R1_001.fastq.gz
    │   └── patientC_rep1_tumor_R2_001.fastq.gz
    ├── main.nf
    └── modules
        ├── analyze_reads.nf
        └── count_lines.nf
    ```

Este directorio contiene datos de secuenciación de extremos pareados de tres pacientes (A, B, C).

Para cada paciente, tenemos muestras de tipo `tumor` (que típicamente provienen de biopsias tumorales) o `normal` (tomadas de tejido sano o sangre).
Si no está familiarizado con el análisis de cáncer, sepa que esto corresponde a un modelo experimental que utiliza muestras pareadas tumor/normal para realizar análisis contrastivos.

Para el paciente A específicamente, tenemos dos conjuntos de réplicas técnicas (repeticiones).

Los archivos de datos de secuenciación están nombrados con la convención típica `_R1_` y `_R2_` para lo que se conoce como 'lecturas directas' y 'lecturas inversas'.

_No se preocupe si no está familiarizado con este diseño experimental, no es fundamental para entender este tutorial._

#### Revisar la tarea

Su desafío es escribir un workflow de Nextflow que:

1. **Cargue** los archivos de entrada usando los métodos de manejo de archivos de Nextflow
2. **Extraiga** metadatos (ID de paciente, réplica, tipo de muestra) de la estructura del nombre de archivo
3. **Agrupe** los archivos pareados (R1/R2) usando `channel.fromFilePairs()`
4. **Procese** los archivos con un módulo de análisis proporcionado
5. **Organice** las salidas en una estructura de directorios basada en los metadatos extraídos

#### Lista de verificación de preparación

¿Cree que está listo para comenzar?

- [ ] Entiendo el objetivo de este curso y sus requisitos previos
- [ ] Mi codespace está en funcionamiento
- [ ] He configurado mi directorio de trabajo correctamente
- [ ] Entiendo la tarea

Si puede marcar todas las casillas, está listo para continuar.

---

## 1. Operaciones básicas con archivos

### 1.1. Identificar el tipo de un objeto con `.class`

Eche un vistazo al archivo de workflow `main.nf`:

```groovy title="main.nf" linenums="1"
#!/usr/bin/env nextflow

workflow {

    // Crear un objeto Path a partir de una ruta de archivo en cadena
    myFile = 'data/patientA_rep1_normal_R1_001.fastq.gz'

    println "${myFile} is of class ${myFile.class}"
}
```

Este es un mini-workflow (sin ningún proceso) que hace referencia a una única ruta de archivo en su workflow, luego la imprime en la consola junto con su clase.

??? info "¿Qué es `.class`?"

    En Nextflow, `.class` nos indica con qué tipo de objeto estamos trabajando. Es como preguntar "¿qué tipo de cosa es esto?" para saber si es una cadena, un número, un archivo u otra cosa.
    Esto nos ayudará a ilustrar la diferencia entre una cadena simple y un objeto Path en las siguientes secciones.

Ejecutemos el workflow:

```bash
nextflow run main.nf
```

??? success "Salida del comando"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `main.nf` [romantic_chandrasekhar] DSL2 - revision: 5a4a89bc3a

    data/patientA_rep1_normal_R1_001.fastq.gz is of class java.lang.String
    ```

Como puede ver, Nextflow imprimió la ruta en cadena exactamente como la escribimos.

Esta es solo una salida de texto; Nextflow aún no ha hecho nada especial con ella.
También hemos confirmado que, en lo que respecta a Nextflow, esto es solo una cadena (de clase `java.lang.String`).
Tiene sentido, ya que aún no le hemos indicado a Nextflow que corresponde a un archivo.

### 1.2. Crear un objeto Path con file()

Podemos indicarle a Nextflow cómo manejar los archivos creando [objetos Path](https://www.nextflow.io/docs/latest/reference/stdlib-types.html#path) a partir de cadenas de rutas.

En nuestro workflow, podemos convertir la cadena de ruta `data/patientA_rep1_normal_R1_001.fastq.gz` en un objeto Path usando el método `file()`, que proporciona acceso a las propiedades y operaciones del archivo.

Edite el `main.nf` para envolver la cadena con `file()` de la siguiente manera:

=== "Después"

    ```groovy title="main.nf" linenums="5" hl_lines="2"
        // Crear un objeto Path a partir de una ruta de archivo en cadena
        myFile = file('data/patientA_rep1_normal_R1_001.fastq.gz')

        println "${myFile} is of class ${myFile.class}"
    ```

=== "Antes"

    ```groovy title="main.nf" linenums="5" hl_lines="2"
        // Crear un objeto Path a partir de una ruta de archivo en cadena
        myFile = 'data/patientA_rep1_normal_R1_001.fastq.gz'

        println "${myFile} is of class ${myFile.class}"
    ```

Ahora ejecute el workflow nuevamente:

```bash
nextflow run main.nf
```

??? success "Salida del comando"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `main.nf` [kickass_coulomb] DSL2 - revision: 5af44b1b59

    /workspaces/training/side-quests/working_with_files/data/patientA_rep1_normal_R1_001.fastq.gz is of class class sun.nio.fs.UnixPath
    ```

Esta vez, verá la ruta absoluta completa en lugar de la ruta relativa que proporcionamos como entrada.

Nextflow ha convertido nuestra cadena en un objeto Path y la ha resuelto a la ubicación real del archivo en el sistema.
La ruta del archivo ahora será absoluta, como en `/workspaces/training/side-quests/working_with_files/data/patientA_rep1_normal_R1_001.fastq.gz`.

Observe también que la clase del objeto Path es `sun.nio.fs.UnixPath`: esta es la forma en que Nextflow representa los archivos locales.
Como veremos más adelante, los archivos remotos tendrán nombres de clase diferentes (como `nextflow.file.http.XPath` para archivos HTTP), pero todos funcionan exactamente de la misma manera y pueden usarse de forma idéntica en sus workflows.

!!! tip "Consejo"

    **La diferencia clave:**

    - **Cadena de ruta**: Solo texto que Nextflow trata como caracteres
    - **Objeto Path**: Una referencia de archivo inteligente con la que Nextflow puede trabajar

    Piénselo así: una cadena de ruta es como escribir una dirección en papel, mientras que un objeto Path es como tener la dirección cargada en un dispositivo GPS que sabe cómo navegar hasta allí y puede darle detalles sobre el trayecto.

### 1.3. Acceder a los atributos del archivo

¿Por qué es útil esto? Ahora que Nextflow entiende que `myFile` es un objeto Path y no solo una cadena, podemos acceder a los diversos atributos del objeto Path.

Actualicemos nuestro workflow para imprimir los atributos de archivo integrados:

=== "Después"

    ```groovy title="main.nf" linenums="5" hl_lines="4-9"
        // Crear un objeto Path a partir de una ruta de archivo en cadena
        myFile = file('data/patientA_rep1_normal_R1_001.fastq.gz')

        // Imprimir atributos del archivo
        println "File object class: ${myFile.class}"
        println "File name: ${myFile.name}"
        println "Simple name: ${myFile.simpleName}"
        println "Extension: ${myFile.extension}"
        println "Parent directory: ${myFile.parent}"
    ```

=== "Antes"

    ```groovy title="main.nf" linenums="5" hl_lines="4"
        // Crear un objeto Path a partir de una ruta de archivo en cadena
        myFile = file('data/patientA_rep1_normal_R1_001.fastq.gz')

        println "${myFile} is of class ${myFile.class}"
    ```

Ejecute el workflow:

```bash
nextflow run main.nf
```

??? success "Salida del comando"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `main.nf` [ecstatic_ampere] DSL2 - revision: f3fa3dcb48

    File object class: sun.nio.fs.UnixPath
    File name: patientA_rep1_normal_R1_001.fastq.gz
    Simple name: patientA_rep1_normal_R1_001
    Extension: gz
    Parent directory: /workspaces/training/side-quests/working_with_files/data
    ```

Puede ver los diversos atributos del archivo impresos en la consola.

### 1.4. Pasar el archivo a un proceso

La diferencia entre cadenas y objetos Path se vuelve crítica cuando comienza a construir workflows reales con procesos.
Hasta ahora hemos verificado que Nextflow está tratando nuestro archivo de entrada como un archivo, pero veamos si realmente podemos ejecutar algo sobre ese archivo en un proceso.

#### 1.4.1. Importar el proceso y examinar el código

Le proporcionamos un módulo de proceso preescrito llamado `COUNT_LINES` que toma un archivo de entrada y cuenta cuántas líneas contiene.

Para usar el proceso en el workflow, solo necesita agregar una declaración include antes del bloque workflow:

=== "Después"

    ```groovy title="main.nf" linenums="1" hl_lines="3"
    #!/usr/bin/env nextflow

    include { COUNT_LINES } from './modules/count_lines.nf'

    workflow {
    ```

=== "Antes"

    ```groovy title="main.nf" linenums="1"
    #!/usr/bin/env nextflow

    workflow {
    ```

Puede abrir el archivo del módulo para examinar su código:

```groovy title="modules/count_lines.nf" linenums="1"
#!/usr/bin/env nextflow

process COUNT_LINES {
    debug true

    input:
    path input_file

    script:
    """
    set -o pipefail
    echo "Processing file: $input_file"
    gzip -dc $input_file | wc -l
    """
}
```

Como puede ver, es un script bastante sencillo que descomprime el archivo y cuenta cuántas líneas contiene.

??? info "¿Qué hace `debug true`?"

    La directiva `debug true` en la definición del proceso hace que Nextflow imprima la salida de su script (como el conteo de líneas "40") directamente en el registro de ejecución.
    Sin esto, solo vería el estado de ejecución del proceso pero no la salida real de su script.

    Para más información sobre la depuración de procesos de Nextflow, consulte la misión secundaria [Debugging Nextflow Workflows](debugging.md).

#### 1.4.2. Agregar una llamada a `COUNT_LINES`

Ahora que el proceso está disponible para el workflow, podemos agregar una llamada al proceso `COUNT_LINES` para ejecutarlo sobre el archivo de entrada.

Realice los siguientes cambios en el workflow:

=== "Después"

    ```groovy title="main.nf" linenums="7" hl_lines="11-12"
        // Crear un objeto Path a partir de una ruta de archivo en cadena
        myFile = file('data/patientA_rep1_normal_R1_001.fastq.gz')

        // Imprimir atributos del archivo
        println "File object class: ${myFile.class}"
        println "File name: ${myFile.name}"
        println "Simple name: ${myFile.simpleName}"
        println "Extension: ${myFile.extension}"
        println "Parent directory: ${myFile.parent}"

        // Contar las líneas del archivo
        COUNT_LINES(myFile)
    ```

=== "Antes"

    ```groovy title="main.nf" linenums="7" hl_lines="4-9"
        // Crear un objeto Path a partir de una ruta de archivo en cadena
        myFile = file('data/patientA_rep1_normal_R1_001.fastq.gz')

        // Imprimir atributos del archivo
        println "File object class: ${myFile.class}"
        println "File name: ${myFile.name}"
        println "Simple name: ${myFile.simpleName}"
        println "Extension: ${myFile.extension}"
        println "Parent directory: ${myFile.parent}"
    ```

Ahora ejecute el workflow:

```bash
nextflow run main.nf
```

??? success "Salida del comando"

    ```console
     N E X T F L O W   ~  version 25.10.2

    Launching `main.nf` [cheeky_hypatia] DSL2 - revision: 281d13c414

    File object class: class sun.nio.fs.UnixPath
    File name: patientA_rep1_normal_R1_001.fastq.gz
    Simple name: patientA_rep1_normal_R1_001
    Extension: gz
    Parent directory: /workspaces/training/side-quests/working_with_files/data
    executor >  local (1)
    [e9/341c05] COUNT_LINES [100%] 1 of 1 ✔
    Processing file: /workspaces/training/side-quests/working_with_files/data/patientA_rep1_normal_R1_001.fastq.gz
    40
    ```

Esto muestra que podemos operar sobre el archivo correctamente dentro de un proceso.

Específicamente, Nextflow realizó las siguientes operaciones con éxito:

- Preparó el archivo en el directorio de trabajo (staging)
- Descomprimió el archivo .gz
- Contó las líneas (40 líneas en este caso)
- Completó sin errores

La clave de esta operación fluida es que le estamos indicando explícitamente a Nextflow que nuestra entrada es un archivo y debe tratarse como tal.

### 1.5. Solucionar errores básicos de entrada de archivos

Esto suele confundir a los recién llegados a Nextflow, así que tomemos unos minutos para ver qué sucede cuando se hace incorrectamente.

Hay dos lugares principales donde puede cometer errores en el manejo de archivos: a nivel del workflow y a nivel del proceso.

#### 1.5.1. Error a nivel del workflow

Veamos qué sucede si volvemos a tratar el archivo como una cadena cuando especificamos la entrada en el bloque workflow.

Realice los siguientes cambios en el workflow, asegurándose de comentar las declaraciones de impresión específicas de la ruta:

=== "Después"

    ```groovy title="main.nf" linenums="7" hl_lines="2 6-11"
        // Crear un objeto Path a partir de una ruta de archivo en cadena
        myFile = 'data/patientA_rep1_normal_R1_001.fastq.gz'

        // Imprimir atributos del archivo
        println "File object class: ${myFile.class}"
        /*
        println "File name: ${myFile.name}"
        println "Simple name: ${myFile.simpleName}"
        println "Extension: ${myFile.extension}"
        println "Parent directory: ${myFile.parent}"
        */

        // Contar las líneas del archivo
        COUNT_LINES(myFile)
    ```

=== "Antes"

    ```groovy title="main.nf" linenums="7" hl_lines="4-9"
        // Crear un objeto Path a partir de una ruta de archivo en cadena
        myFile = file('data/patientA_rep1_normal_R1_001.fastq.gz')

        // Imprimir atributos del archivo
        println "File object class: ${myFile.class}"
        println "File name: ${myFile.name}"
        println "Simple name: ${myFile.simpleName}"
        println "Extension: ${myFile.extension}"
        println "Parent directory: ${myFile.parent}"

        // Contar las líneas del archivo
        COUNT_LINES(myFile)
    ```

Ahora ejecute el workflow:

```bash
nextflow run main.nf
```

??? failure "Salida del comando"

    ```console
     N E X T F L O W   ~  version 25.10.2

    Launching `main.nf` [friendly_goodall] DSL2 - revision: ae50609b20

    [-        ] COUNT_LINES -
    ERROR ~ Error executing process > 'COUNT_LINES'

    Caused by:
      Not a valid path value: 'data/patientA_rep1_normal_R1_001.fastq.gz'



    Tip: view the complete command output by changing to the process work dir and entering the command `cat .command.out`

    -- Check '.nextflow.log' file for details
    ```

Esta es la parte importante:

```console
Not a valid path value: 'data/patientA_rep1_normal_R1_001.fastq.gz'
```

Cuando especifica una entrada `path`, Nextflow valida que está pasando referencias de archivos reales, no solo cadenas.
Este error le indica que `'data/patientA_rep1_normal_R1_001.fastq.gz'` no es un valor de ruta válido porque es una cadena, no un objeto Path.

Nextflow detectó el problema de inmediato y se detuvo antes de iniciar el proceso.

#### 1.5.2. Error a nivel del proceso

El otro lugar donde podríamos olvidar especificar que queremos que Nextflow trate la entrada como un archivo es en la definición del proceso.

!!! warning "Advertencia"

    **Mantenga el error del workflow del paso 1.5.1**

    Para que esta prueba funcione correctamente, mantenga el workflow en su estado incorrecto (usando una cadena simple en lugar de `file()`).
    Cuando se combina con `val` en el proceso, esto produce el error que se muestra a continuación.

Realice el siguiente cambio en el módulo:

=== "Después"

    ```groovy title="modules/count_lines.nf" linenums="3" hl_lines="5"
    process COUNT_LINES {
        debug true

        input:
        val input_file
    ```

=== "Antes"

    ```groovy title="modules/count_lines.nf" linenums="3" hl_lines="5"
    process COUNT_LINES {
        debug true

        input:
        path input_file
    ```

Ahora ejecute el workflow nuevamente:

```bash
nextflow run main.nf
```

??? failure "Salida del comando"

    ```console
     N E X T F L O W   ~  version 25.10.2

    Launching `main.nf` [soggy_golick] DSL2 - revision: ae50609b20

    executor >  local (1)
    [b3/b3023c] COUNT_LINES [  0%] 0 of 1 ✘
    ERROR ~ Error executing process > 'COUNT_LINES'

    Caused by:
      Process `COUNT_LINES` terminated with an error exit status (1)


    Command executed:

      set -o pipefail
      echo "Processing file: data/patientA_rep1_normal_R1_001.fastq.gz"
      gzip -dc data/patientA_rep1_normal_R1_001.fastq.gz | wc -l

    Command exit status:
      1

    Command output:
      Processing file: data/patientA_rep1_normal_R1_001.fastq.gz
      0

    Command error:
      Processing file: data/patientA_rep1_normal_R1_001.fastq.gz
      gzip: data/patientA_rep1_normal_R1_001.fastq.gz: No such file or directory
      0

    Work dir:
      /workspaces/training/side-quests/working_with_files/work/b3/b3023cb2ccb986851301d8e369e79f

    Tip: when you have fixed the problem you can continue the execution adding the option `-resume` to the run command line

    -- Check '.nextflow.log' file for details
    ```

Esto muestra muchos detalles sobre el error porque el proceso está configurado para mostrar información de depuración, como se mencionó anteriormente.

Estas son las secciones más relevantes:

```console
Command executed:

  set -o pipefail
  echo "Processing file: data/patientA_rep1_normal_R1_001.fastq.gz"
  gzip -dc data/patientA_rep1_normal_R1_001.fastq.gz | wc -l
```

```console
Command error:
  Processing file: data/patientA_rep1_normal_R1_001.fastq.gz
  gzip: data/patientA_rep1_normal_R1_001.fastq.gz: No such file or directory
  0
```

Esto indica que el sistema no pudo encontrar el archivo; sin embargo, si busca la ruta, hay un archivo con ese nombre en esa ubicación.

Cuando ejecutamos esto, Nextflow pasó el valor de la cadena al script, pero no preparó (_staged_) el archivo real en el directorio de trabajo.
Por lo tanto, el proceso intentó usar la cadena relativa, `data/patientA_rep1_normal_R1_001.fastq.gz`, pero ese archivo no existe dentro del directorio de trabajo del proceso.

En conjunto, estos dos ejemplos muestran lo importante que es indicarle a Nextflow si una entrada debe manejarse como un archivo.

!!! note "Nota"

    Asegúrese de volver y corregir ambos errores intencionales antes de continuar con la siguiente sección.

### Conclusión

- Cadenas de ruta vs objetos Path: Las cadenas son solo texto, los objetos Path son referencias de archivo inteligentes
- El método `file()` convierte una cadena de ruta en un objeto Path con el que Nextflow puede trabajar
- Puede acceder a propiedades del archivo como `name`, `simpleName`, `extension` y `parent` [usando atributos de archivo](https://www.nextflow.io/docs/latest/working-with-files.html#getting-file-attributes)
- Usar objetos Path en lugar de cadenas permite a Nextflow gestionar correctamente los archivos en su workflow
- Resultados del manejo de entradas en procesos: El manejo adecuado de archivos requiere objetos Path, no cadenas, para garantizar que los archivos sean correctamente preparados y accesibles para su uso por los procesos.

---

## 2. Uso de archivos remotos

Una de las características clave de Nextflow es la capacidad de cambiar sin problemas entre archivos locales (en la misma máquina) y archivos remotos accesibles a través de internet.

Si lo hace correctamente, nunca debería necesitar cambiar la lógica de su workflow para acomodar archivos provenientes de diferentes ubicaciones.
Todo lo que necesita hacer para usar un archivo remoto es especificar el prefijo apropiado en la ruta del archivo cuando lo proporciona al workflow.

Por ejemplo, `/path/to/data` no tiene prefijo, lo que indica que es una ruta de archivo local 'normal', mientras que `s3://path/to/data` incluye el prefijo `s3://`, lo que indica que está ubicado en el almacenamiento de objetos S3 de Amazon.

Se admiten muchos protocolos diferentes:

- HTTP(S)/FTP (http://, https://, ftp://)
- Amazon S3 (s3://)
- Azure Blob Storage (az://)
- Google Cloud Storage (gs://)

Para usar cualquiera de estos, simplemente especifique el prefijo relevante en la cadena, que técnicamente se denomina Identificador Uniforme de Recursos (URI) en lugar de ruta de archivo.
Nextflow manejará la autenticación y la preparación de los archivos en el lugar correcto, descargando o subiendo y realizando todas las demás operaciones de archivo que esperaría.

La principal fortaleza de este sistema es que nos permite cambiar entre entornos sin modificar ninguna lógica del pipeline.
Por ejemplo, puede desarrollar con un conjunto de pruebas pequeño y local antes de cambiar a un conjunto de pruebas a escala completa ubicado en almacenamiento remoto simplemente cambiando el URI.

### 2.1. Usar un archivo de internet

Probemos esto cambiando la ruta local que estamos proporcionando a nuestro workflow por una ruta HTTPS que apunta a una copia de los mismos datos almacenados en Github.

!!! warning "Advertencia"

    Esto solo funcionará si tiene una conexión a internet activa.

Abra `main.nf` nuevamente y cambie la ruta de entrada de la siguiente manera:

=== "Después"

    ```groovy title="main.nf" linenums="2" hl_lines="2"
        // Usar un archivo remoto de internet
        myFile = file('https://raw.github.com/nextflow-io/training/master/side-quests/working_with_files/data/patientA_rep1_normal_R1_001.fastq.gz')

        // Imprimir atributos del archivo
        println "File object class: ${myFile.class}"
        println "File name: ${myFile.name}"
        println "Simple name: ${myFile.simpleName}"
        println "Extension: ${myFile.extension}"
        println "Parent directory: ${myFile.parent}"
    ```

=== "Antes"

    ```groovy title="main.nf" linenums="2" hl_lines="2"
        // Crear un objeto Path a partir de una ruta de archivo en cadena
        myFile = file('data/patientA_rep1_normal_R1_001.fastq.gz')

        // Imprimir atributos del archivo
        println "File object class: ${myFile.class}"
        println "File name: ${myFile.name}"
        println "Simple name: ${myFile.simpleName}"
        println "Extension: ${myFile.extension}"
        println "Parent directory: ${myFile.parent}"
    ```

Ejecutemos el workflow:

```bash
nextflow run main.nf
```

??? success "Salida del comando"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `main.nf` [insane_swartz] DSL2 - revision: fff18abe6d

    File object class: class nextflow.file.http.XPath
    File name: patientA_rep1_normal_R1_001.fastq.gz
    Simple name: patientA_rep1_normal_R1_001
    Extension: gz
    Parent directory: /nextflow-io/training/master/side-quests/working_with_files/data
    executor >  local (1)
    [8a/2ab7ca] COUNT_LINES [100%] 1 of 1 ✔
    Processing file: patientA_rep1_normal_R1_001.fastq.gz
    40
    ```

¡Funciona! Puede ver que muy poco ha cambiado.

La única diferencia en la salida de la consola es que la clase del objeto Path ahora es `nextflow.file.http.XPath`, mientras que para la ruta local la clase era `sun.nio.fs.UnixPath`.
No necesita recordar estas clases; solo las mencionamos para demostrar que Nextflow identifica y maneja las diferentes ubicaciones de manera apropiada.

En segundo plano, Nextflow descargó el archivo a un directorio de preparación ubicado dentro del directorio de trabajo.
Ese archivo preparado puede tratarse como un archivo local y vincularse simbólicamente al directorio del proceso correspondiente.

Puede verificar que eso ocurrió aquí mirando el contenido del directorio de trabajo ubicado en el valor hash del proceso.

??? abstract "Contenido del directorio work"

    Si el hash del proceso fue `8a/2ab7ca`, podría explorar el directorio de trabajo:

    ```console
    $ ls -la work/8a/2ab7ca*/
    total 16
    drwxr-xr-x  6 user  staff   192 Jan 28 10:00 .
    drwxr-xr-x  3 user  staff    96 Jan 28 10:00 ..
    -rw-r--r--  1 user  staff     0 Jan 28 10:00 .command.begin
    -rw-r--r--  1 user  staff   127 Jan 28 10:00 .command.sh
    lrwxr-xr-x  1 user  staff    89 Jan 28 10:00 patientA_rep1_normal_R1_001.fastq.gz -> /path/to/work/stage/.../patientA_rep1_normal_R1_001.fastq.gz
    ```

    El enlace simbólico apunta a una copia preparada del archivo remoto que Nextflow descargó automáticamente.

Tenga en cuenta que para archivos más grandes, el paso de descarga tomará tiempo adicional en comparación con la ejecución sobre archivos locales.
Sin embargo, Nextflow verifica si ya tiene una copia preparada para evitar descargas innecesarias.
Por lo tanto, si ejecuta nuevamente sobre el mismo archivo y no ha eliminado el archivo preparado, Nextflow usará la copia preparada.

Esto muestra lo fácil que es cambiar entre datos locales y remotos usando Nextflow, lo cual es una característica clave de Nextflow.

!!! note "Nota"

    La única excepción importante a este principio es que no puede usar patrones glob o rutas de directorio con HTTPS porque HTTPS no puede listar múltiples archivos, por lo que debe especificar URLs de archivo exactas.
    Sin embargo, otros protocolos de almacenamiento como el almacenamiento en la nube (`s3://`, `az://`, `gs://`) pueden usar tanto globs como rutas de directorio.

    Así es como podría usar patrones glob con almacenamiento en la nube:

    ```groovy title="Cloud storage examples (not runnable in this environment)"
    // S3 con patrones glob - coincidiría con múltiples archivos
    ch_s3_files = channel.fromPath('s3://my-bucket/data/*.fastq.gz')

    // Azure Blob Storage con patrones glob
    ch_azure_files = channel.fromPath('az://container/data/patient*_R{1,2}.fastq.gz')

    // Google Cloud Storage con patrones glob
    ch_gcs_files = channel.fromPath('gs://bucket/data/sample_*.fastq.gz')
    ```

    Le mostraremos cómo trabajar con globs en la práctica en la siguiente sección.

### 2.2. Volver al archivo local

Vamos a volver a usar nuestros archivos de ejemplo locales para el resto de esta misión secundaria, así que cambiemos la entrada del workflow de vuelta al archivo original:

=== "Después"

    ```groovy title="main.nf" linenums="2" hl_lines="2"
        // Crear un objeto Path a partir de una ruta de archivo en cadena
        myFile = file('data/patientA_rep1_normal_R1_001.fastq.gz')

        // Imprimir atributos del archivo
        println "File object class: ${myFile.class}"
        println "File name: ${myFile.name}"
        println "Simple name: ${myFile.simpleName}"
        println "Extension: ${myFile.extension}"
        println "Parent directory: ${myFile.parent}"
    ```

=== "Antes"

    ```groovy title="main.nf" linenums="2" hl_lines="2"
        // Crear un objeto Path a partir de una ruta de archivo en cadena
        myFile = file('https://raw.github.com/nextflow-io/training/master/side-quests/working_with_files/data/patientA_rep1_normal_R1_001.fastq.gz')

        // Imprimir atributos del archivo
        println "File object class: ${myFile.class}"
        println "File name: ${myFile.name}"
        println "Simple name: ${myFile.simpleName}"
        println "Extension: ${myFile.extension}"
        println "Parent directory: ${myFile.parent}"
    ```

### Conclusión

- Los datos remotos se acceden usando un URI (HTTP, FTP, S3, Azure, Google Cloud)
- Nextflow descargará y preparará automáticamente los datos en el lugar correcto, siempre que estas rutas se estén pasando a los procesos
- ¡No escriba lógica para descargar o subir archivos remotos!
- Los archivos locales y remotos producen diferentes tipos de objetos pero funcionan de manera idéntica
- **Importante**: HTTP/HTTPS solo funciona con archivos individuales (sin patrones glob)
- El almacenamiento en la nube (S3, Azure, GCS) admite tanto archivos individuales como patrones glob
- Puede cambiar sin problemas entre fuentes de datos locales y remotas sin cambiar la lógica del código (siempre que el protocolo admita las operaciones requeridas)

---

## 3. Uso de la fábrica de canales `fromPath()`

Hasta ahora hemos estado trabajando con un solo archivo a la vez, pero en Nextflow, típicamente querremos crear un canal de entrada con múltiples archivos de entrada para procesar.

Una forma ingenua de hacerlo sería combinar el método `file()` con [`channel.of()`](https://www.nextflow.io/docs/latest/reference/channel.html#of) así:

```groovy title="Syntax example"
ch_files = channel.of([file('data/patientA_rep1_normal_R1_001.fastq.gz')],
                      [file('data/patientA_rep1_normal_R1_001.fastq.gz')])
```

Eso funciona, pero es engorroso.

!!! tip "Consejo: Cuándo usar `file()` vs `channel.fromPath()`"

    - Use `file()` cuando necesite un único objeto Path para manipulación directa (verificar si un archivo existe, leer sus atributos, o pasarlo a una única invocación de proceso)
    - Use `channel.fromPath()` cuando necesite un canal que pueda contener múltiples archivos, especialmente con patrones glob, o cuando los archivos fluirán a través de múltiples procesos

Aquí es donde entra [`channel.fromPath()`](https://www.nextflow.io/docs/latest/reference/channel.html#frompath): una fábrica de canales conveniente que agrupa toda la funcionalidad que necesitamos para generar un canal a partir de una o más cadenas de archivos estáticas, así como patrones glob.

### 3.1. Agregar la fábrica de canales

Actualicemos nuestro workflow para usar `channel.fromPath`.

=== "Después"

    ```groovy title="main.nf" linenums="7" hl_lines="1-3"
        // Cargar archivos con channel.fromPath
        ch_files = channel.fromPath('data/patientA_rep1_normal_R1_001.fastq.gz')
        ch_files.view { myFile -> "Found file: $myFile" }

        // Imprimir atributos del archivo
        /* ¡Comentamos esto por ahora, volveremos a ello!
        println "File object class: ${myFile.class}"
        println "File name: ${myFile.name}"
        println "Simple name: ${myFile.simpleName}"
        println "Extension: ${myFile.extension}"
        println "Parent directory: ${myFile.parent}"
        */

        // Contar las líneas del archivo
        // COUNT_LINES(myFile)
    ```

=== "Antes"

    ```groovy title="main.nf" linenums="7" hl_lines="1-2"
        // Crear un objeto Path a partir de una ruta de archivo en cadena
        myFile = file('data/patientA_rep1_normal_R1_001.fastq.gz')

        // Imprimir atributos del archivo
        println "File object class: ${myFile.class}"
        println "File name: ${myFile.name}"
        println "Simple name: ${myFile.simpleName}"
        println "Extension: ${myFile.extension}"
        println "Parent directory: ${myFile.parent}"

        // Contar las líneas del archivo
        COUNT_LINES(myFile)
    ```

También hemos comentado el código que imprime los atributos por ahora, y hemos agregado una declaración `.view` para imprimir solo el nombre del archivo.

Ejecute el workflow:

```bash
nextflow run main.nf
```

??? success "Salida del comando"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `main.nf` [grave_meucci] DSL2 - revision: b09964a583

    Found file: /workspaces/training/side-quests/working_with_files/data/patientA_rep1_normal_R1_001.fastq.gz
    ```

Como puede ver, la ruta del archivo se está cargando como un objeto de tipo `Path` en el canal.
Esto es similar a lo que habría hecho `file()`, excepto que ahora tenemos un canal en el que podemos cargar más archivos si queremos.

Usar `channel.fromPath()` es una forma conveniente de crear un nuevo canal poblado con una lista de archivos.

### 3.2. Ver atributos de los archivos en el canal

En nuestra primera aproximación al uso de la fábrica de canales, simplificamos el código y solo imprimimos el nombre del archivo.

Volvamos a imprimir los atributos completos del archivo:

=== "Después"

    ```groovy title="main.nf" linenums="7" hl_lines="3-9 12"
        // Cargar archivos con channel.fromPath
        ch_files = channel.fromPath('data/patientA_rep1_normal_R1_001.fastq.gz')
        ch_files.view { myFile ->
            println "File object class: ${myFile.class}"
            println "File name: ${myFile.name}"
            println "Simple name: ${myFile.simpleName}"
            println "Extension: ${myFile.extension}"
            println "Parent directory: ${myFile.parent}"
        }

        // Contar las líneas del archivo
        COUNT_LINES(ch_files)
    ```

=== "Antes"

    ```groovy title="main.nf" linenums="7" hl_lines="3"
        // Cargar archivos con channel.fromPath
        ch_files = channel.fromPath('data/patientA_rep1_normal_R1_001.fastq.gz')
        ch_files.view { myFile -> "Found file: $myFile" }

        // Contar las líneas del archivo
        // COUNT_LINES(ch_files)
    ```

También estamos reactivando la llamada al proceso `COUNT_LINES` para verificar que el procesamiento de archivos sigue funcionando correctamente con nuestro enfoque basado en canales.

Ejecute el workflow:

```bash
nextflow run main.nf
```

??? success "Salida del comando"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `main.nf` [furious_swanson] DSL2 - revision: c35c34950d

    executor >  local (1)
    [9d/6701a6] COUNT_LINES (1) [100%] 1 of 1 ✔
    File object class: sun.nio.fs.UnixPath
    File name: patientA_rep1_normal_R1_001.fastq.gz
    Simple name: patientA_rep1_normal_R1_001
    Extension: gz
    Parent directory: /workspaces/training/side-quests/working_with_files/data
    Processing file: patientA_rep1_normal_R1_001.fastq.gz
    40
    ```

Y ahí está, los mismos resultados que antes pero ahora tenemos el archivo en un canal, por lo que podemos agregar más.

### 3.3. Usar un glob para coincidir con múltiples archivos

Hay varias formas de cargar más archivos en el canal.
Aquí le mostraremos cómo usar patrones glob, que son una forma conveniente de coincidir y recuperar nombres de archivos y directorios basados en caracteres comodín.
El proceso de coincidir estos patrones se llama "globbing" o "expansión de nombres de archivo".

!!! note "Nota"

    Como se mencionó anteriormente, Nextflow admite globbing para gestionar archivos de entrada y salida en la mayoría de los casos, excepto con rutas HTTPS porque HTTPS no puede listar múltiples archivos.

Digamos que queremos recuperar ambos archivos en un par de archivos asociados con un paciente dado, `patientA`:

```console
patientA_rep1_normal_R1_001.fastq.gz
patientA_rep1_normal_R2_001.fastq.gz
```

Dado que la única diferencia entre los nombres de archivo es el número de réplica, _es decir_, el número después de `R`, podemos usar el carácter comodín `*` para sustituir el número de la siguiente manera:

```console
patientA_rep1_normal_R*_001.fastq.gz
```

Ese es el patrón glob que necesitamos.

Ahora todo lo que necesitamos hacer es actualizar la ruta del archivo en la fábrica de canales para usar ese patrón glob de la siguiente manera:

=== "Después"

    ```groovy title="main.nf" linenums="7"
      // Cargar archivos con channel.fromPath
        ch_files = channel.fromPath('data/patientA_rep1_normal_R*_001.fastq.gz')
    ```

=== "Antes"

    ```groovy title="main.nf" linenums="7"
      // Cargar archivos con channel.fromPath
        ch_files = channel.fromPath('data/patientA_rep1_normal_R1_001.fastq.gz')
    ```

Nextflow reconocerá automáticamente que este es un patrón glob y lo manejará apropiadamente.

Ejecute el workflow para probarlo:

```bash
nextflow run main.nf
```

??? success "Salida del comando"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `main.nf` [boring_sammet] DSL2 - revision: d2aa789c9a

    executor >  local (2)
    [3c/a65de5] COUNT_LINES (2) [100%] 2 of 2 ✔
    File object class: class sun.nio.fs.UnixPath
    File name: patientA_rep1_normal_R1_001.fastq.gz
    Simple name: patientA_rep1_normal_R1_001
    Extension: gz
    Parent directory: /workspaces/training/side-quests/working_with_files/data
    File object class: class sun.nio.fs.UnixPath
    File name: patientA_rep1_normal_R2_001.fastq.gz
    Simple name: patientA_rep1_normal_R2_001
    Extension: gz
    Parent directory: /workspaces/training/side-quests/working_with_files/data
    Processing file: patientA_rep1_normal_R1_001.fastq.gz
    40

    Processing file: patientA_rep1_normal_R2_001.fastq.gz
    40
    ```

Como puede ver, ahora tenemos dos objetos Path en nuestro canal, lo que muestra que Nextflow ha realizado la expansión de nombres de archivo correctamente, y ha cargado y procesado ambos archivos como se esperaba.

Usando este método, podemos recuperar tantos o tan pocos archivos como queramos simplemente cambiando el patrón glob. Si lo hacemos más amplio, por ejemplo reemplazando todas las partes variables de los nombres de archivo por `*` (_p. ej._ `data/patient*_rep*_*_R*_001.fastq.gz`) podríamos obtener todos los archivos de ejemplo en el directorio `data`.

### Conclusión

- `channel.fromPath()` crea un canal con archivos que coinciden con un patrón
- Cada archivo se emite como un elemento separado en el canal
- Podemos usar un patrón glob para coincidir con múltiples archivos
- Los archivos se convierten automáticamente en objetos Path con todos sus atributos
- El método `.view()` permite inspeccionar el contenido del canal

---

## 4. Extracción de metadatos básicos de los nombres de archivo

En la mayoría de los dominios científicos, es muy común tener metadatos codificados en los nombres de los archivos que contienen los datos.
Por ejemplo, en bioinformática, los archivos que contienen datos de secuenciación a menudo se nombran de una manera que codifica información sobre la muestra, la condición, la réplica y el número de lectura.

Si los nombres de archivo se construyen según una convención consistente, puede extraer esos metadatos de manera estandarizada y usarlos en el transcurso de su análisis.
Eso es un gran 'si', por supuesto, y debe ser muy cauteloso cuando dependa de la estructura del nombre de archivo; pero la realidad es que este enfoque es muy ampliamente utilizado, así que veamos cómo se hace en Nextflow.

En el caso de nuestros datos de ejemplo, sabemos que los nombres de archivo incluyen metadatos estructurados de manera consistente.
Por ejemplo, el nombre de archivo `patientA_rep1_normal_R2_001` codifica lo siguiente:

- ID de paciente: `patientA`
- ID de réplica: `rep1`
- tipo de muestra: `normal` (en contraposición a `tumor`)
- conjunto de lecturas: `R1` (en contraposición a `R2`)

Vamos a modificar nuestro workflow para recuperar esta información en tres pasos:

1. Recuperar el `simpleName` del archivo, que incluye los metadatos
2. Separar los metadatos usando un método llamado `tokenize()`
3. Usar un map para organizar los metadatos

!!! warning "Advertencia"

    Nunca debe codificar información sensible en los nombres de archivo, como nombres de pacientes u otras características identificativas, ya que eso puede comprometer la privacidad del paciente u otras restricciones de seguridad relevantes.

### 4.1. Recuperar el `simpleName`

El `simpleName` es un atributo de archivo que corresponde al nombre del archivo sin su ruta ni extensión.

Realice los siguientes cambios en el workflow:

=== "Después"

    ```groovy title="main.nf" linenums="7" hl_lines="3-6"
        // Cargar archivos con channel.fromPath
        ch_files = channel.fromPath('data/patientA_rep1_normal_R*_001.fastq.gz')
        ch_files.map { myFile ->
            [ myFile.simpleName, myFile ]
        }
        .view()
    ```

=== "Antes"

    ```groovy title="main.nf" linenums="7" hl_lines="3-9"
        // Cargar archivos con channel.fromPath
        ch_files = channel.fromPath('data/patientA_rep1_normal_R*_001.fastq.gz')
        ch_files.view { myFile ->
            println "File object class: ${myFile.class}"
            println "File name: ${myFile.name}"
            println "Simple name: ${myFile.simpleName}"
            println "Extension: ${myFile.extension}"
            println "Parent directory: ${myFile.parent}"
        }
    ```

Esto recupera el `simpleName` y lo asocia con el objeto de archivo completo usando una operación `map()`.

Ejecute el workflow para verificar que funciona:

```bash
nextflow run main.nf
```

??? success "Salida del comando"

    ```console hl_lines="7-8"
    N E X T F L O W   ~  version 25.10.2

    Launching `main.nf` [suspicious_mahavira] DSL2 - revision: ae8edc4e48

    executor >  local (2)
    [e9/55774b] COUNT_LINES (2) [100%] 2 of 2 ✔
    [patientA_rep1_normal_R2_001, /workspaces/training/side-quests/working_with_files/data/patientA_rep1_normal_R2_001.fastq.gz]
    [patientA_rep1_normal_R1_001, /workspaces/training/side-quests/working_with_files/data/patientA_rep1_normal_R1_001.fastq.gz]
    Processing file: patientA_rep1_normal_R1_001.fastq.gz
    40

    Processing file: patientA_rep1_normal_R2_001.fastq.gz
    40
    ```

Cada elemento en el canal es ahora una tupla que contiene el `simpleName` y el objeto de archivo original.

### 4.2. Extraer los metadatos del `simpleName`

En este punto, los metadatos que queremos están integrados en el `simpleName`, pero no podemos acceder a los elementos individuales directamente.
Por lo tanto, necesitamos dividir el `simpleName` en sus componentes.
Afortunadamente, esos componentes están simplemente separados por guiones bajos en el nombre de archivo original, por lo que podemos aplicar un método común de Nextflow llamado `tokenize()` que es perfecto para esta tarea.

Realice los siguientes cambios en el workflow:

=== "Después"

    ```groovy title="main.nf" linenums="7" hl_lines="4"
        // Cargar archivos con channel.fromPath
        ch_files = channel.fromPath('data/patientA_rep1_normal_R*_001.fastq.gz')
        ch_files.map { myFile ->
            [ myFile.simpleName.tokenize('_'), myFile ]
        }
    ```

=== "Antes"

    ```groovy title="main.nf" linenums="7" hl_lines="4"
        // Cargar archivos con channel.fromPath
        ch_files = channel.fromPath('data/patientA_rep1_normal_R*_001.fastq.gz')
        ch_files.map { myFile ->
            [ myFile.simpleName, myFile ]
        }
    ```

El método `tokenize()` dividirá la cadena `simpleName` donde encuentre guiones bajos, y devolverá una lista que contiene las subcadenas.

Ejecute el workflow:

```bash
nextflow run main.nf
```

??? success "Salida del comando"

    ```console hl_lines="7-8"
    N E X T F L O W   ~  version 25.10.2

    Launching `main.nf` [gigantic_gauss] DSL2 - revision: a39baabb57

    executor >  local (2)
    [e7/da2f4b] COUNT_LINES (2) [100%] 2 of 2 ✔
    [[patientA, rep1, normal, R2, 001], /workspaces/training/side-quests/working_with_files/data/patientA_rep1_normal_R2_001.fastq.gz]
    [[patientA, rep1, normal, R1, 001], /workspaces/training/side-quests/working_with_files/data/patientA_rep1_normal_R1_001.fastq.gz]
    Processing file: patientA_rep1_normal_R2_001.fastq.gz
    40

    Processing file: patientA_rep1_normal_R1_001.fastq.gz
    40
    ```

Ahora la tupla para cada elemento en nuestro canal contiene la lista de metadatos (_p. ej._ `[patientA, rep1, normal, R1, 001]`) y el objeto de archivo original.

¡Excelente!
Hemos desglosado la información de nuestro paciente de una sola cadena en una lista de cadenas.
Ahora podemos manejar cada parte de la información del paciente por separado.

### 4.3. Usar un map para organizar los metadatos

Nuestros metadatos son solo una lista plana en este momento.
Es suficientemente fácil de usar pero difícil de leer.

```console
[patientA, rep1, normal, R1, 001]
```

¿Cuál es el elemento en el índice 3? ¿Puede saberlo sin consultar la explicación original de la estructura de metadatos?

Esta es una gran oportunidad para usar un almacén de clave-valor, donde cada elemento tiene un conjunto de claves y sus valores asociados, para que pueda referirse fácilmente a cada clave para obtener el valor correspondiente.

En nuestro ejemplo, eso significa pasar de esta organización:

```groovy
data = [patientA, 1, normal, R1]

println data[3]
```

A esta:

```groovy
data = [id: patientA, replicate: 1, type: normal, readNum: 1]

println data.readNum
```

En Nextflow, eso se llama un [map](https://nextflow.io/docs/latest/script.html#maps).

Convirtamos nuestra lista plana en un map ahora.
Realice los siguientes cambios en el workflow:

=== "Después"

    ```groovy title="main.nf" linenums="7" hl_lines="4-13"
        // Cargar archivos con channel.fromPath
        ch_files = channel.fromPath('data/patientA_rep1_normal_R*_001.fastq.gz')
        ch_files.map { myFile ->
            def (patient, replicate, type, readNum) = myFile.simpleName.tokenize('_')
            [
              [
                id: patient,
                replicate: replicate.replace('rep', ''),
                type: type,
                readNum: readNum.replace('R', ''),
              ],
              myFile
            ]
        }
    ```

=== "Antes"

    ```groovy title="main.nf" linenums="7" hl_lines="4"
        // Cargar archivos con channel.fromPath
        ch_files = channel.fromPath('data/patientA_rep1_normal_R*_001.fastq.gz')
        ch_files.map { myFile ->
            [ myFile.simpleName.tokenize('_'), myFile ]
        }
    ```

Los cambios clave aquí son:

- **Asignación por desestructuración**: `def (patient, replicate, type, readNum) = ...` extrae los valores tokenizados en variables con nombre en una sola línea
- **Sintaxis literal de map**: `[id: patient, replicate: ...]` crea un map donde cada clave (como `id`) está asociada con un valor (como `patient`)
- **Estructura anidada**: La lista exterior `[..., myFile]` empareja el map de metadatos con el objeto de archivo original

También simplificamos un par de las cadenas de metadatos usando un método de reemplazo de cadenas llamado `replace()` para eliminar algunos caracteres innecesarios (_p. ej._ `replicate.replace('rep', '')` para conservar solo el número de los IDs de réplica).

Ejecutemos el workflow nuevamente:

```bash
nextflow run main.nf
```

??? success "Salida del comando"

    ```console hl_lines="7-8"
    N E X T F L O W   ~  version 25.10.2

    Launching `main.nf` [infallible_swartz] DSL2 - revision: 7f4e68c0cb

    executor >  local (2)
    [1b/e7fb27] COUNT_LINES (1) [100%] 2 of 2 ✔
    [[id:patientA, replicate:1, type:normal, readNum:2], /workspaces/training/side-quests/working_with_files/data/patientA_rep1_normal_R2_001.fastq.gz]
    [[id:patientA, replicate:1, type:normal, readNum:1], /workspaces/training/side-quests/working_with_files/data/patientA_rep1_normal_R1_001.fastq.gz]
    Processing file: patientA_rep1_normal_R2_001.fastq.gz
    40

    Processing file: patientA_rep1_normal_R1_001.fastq.gz
    40
    ```

Ahora los metadatos están claramente etiquetados (_p. ej._ `[id:patientA, replicate:1, type:normal, readNum:2]`) por lo que es mucho más fácil saber qué es qué.

También será mucho más fácil hacer uso real de los elementos de metadatos en el workflow, y hará que nuestro código sea más fácil de leer y mantener.

### Conclusión

- Podemos manejar nombres de archivo en Nextflow con el poder de un lenguaje de programación completo
- Podemos tratar los nombres de archivo como cadenas para extraer información relevante
- El uso de métodos como `tokenize()` y `replace()` nos permite manipular cadenas en el nombre de archivo
- La operación `.map()` transforma los elementos del canal preservando la estructura
- Los metadatos estructurados (maps) hacen que el código sea más legible y mantenible que las listas posicionales

A continuación, veremos cómo manejar archivos de datos pareados.

---

## 5. Manejo de archivos de datos pareados

Muchos diseños experimentales producen archivos de datos pareados que se benefician de ser manejados de manera explícitamente pareada.
Por ejemplo, en bioinformática, los datos de secuenciación a menudo se generan en forma de lecturas pareadas, es decir, cadenas de secuencias que se originan del mismo fragmento de ADN (a menudo llamadas 'directas' e 'inversas' porque se leen desde extremos opuestos).

Ese es el caso de nuestros datos de ejemplo, donde R1 y R2 se refieren a los dos conjuntos de lecturas.

```console
data/patientA_rep1_normal_R1_001.fastq.gz
data/patientA_rep1_normal_R2_001.fastq.gz
```

Nextflow proporciona una fábrica de canales especializada para trabajar con archivos pareados como este llamada `channel.fromFilePairs()`, que agrupa automáticamente los archivos basándose en un patrón de nomenclatura compartido. Eso le permite asociar los archivos pareados de manera más estrecha con menos esfuerzo.

Vamos a modificar nuestro workflow para aprovechar esto.
Tomará dos pasos:

1. Cambiar la fábrica de canales a `channel.fromFilePairs()`
2. Extraer y mapear los metadatos

### 5.1. Cambiar la fábrica de canales a `channel.fromFilePairs()`

Para usar `channel.fromFilePairs`, necesitamos especificar el patrón que Nextflow debe usar para identificar los dos miembros de un par.

Volviendo a nuestros datos de ejemplo, podemos formalizar el patrón de nomenclatura de la siguiente manera:

```console
data/patientA_rep1_normal_R{1,2}_001.fastq.gz
```

Esto es similar al patrón glob que usamos anteriormente, excepto que este enumera específicamente las subcadenas (ya sea `1` o `2` justo después de la R) que identifican a los dos miembros del par.

Actualicemos el workflow `main.nf` en consecuencia:

=== "Después"

    ```groovy title="main.nf" linenums="7" hl_lines="1-2"
        // Cargar archivos con channel.fromFilePairs
        ch_files = channel.fromFilePairs('data/patientA_rep1_normal_R{1,2}_001.fastq.gz')
        /* ¡Comentamos el mapeo por ahora, volveremos a él!
        ch_files.map { myFile ->
            def (sample, replicate, type, readNum) = myFile.simpleName.tokenize('_')
            [
                [
                    id: sample,
                    replicate: replicate.replace('rep', ''),
                    type: type,
                    readNum: readNum,
                ],
                myFile
            ]
        }
        */
        .view()
    ```

=== "Antes"

    ```groovy title="main.nf" linenums="7" hl_lines="1-2"
        // Cargar archivos con channel.fromFilePairs
        ch_files = channel.fromPath('data/patientA_rep1_normal_R*_001.fastq.gz')
        ch_files.map { myFile ->
            def (sample, replicate, type, readNum) = myFile.simpleName.tokenize('_')
            [
                [
                    id: sample,
                    replicate: replicate.replace('rep', ''),
                    type: type,
                    readNum: readNum,
                ],
                myFile
            ]
        }
        .view()
    ```

Hemos cambiado la fábrica de canales y adaptado el patrón de coincidencia de archivos, y de paso, hemos comentado la operación map.
La agregaremos de nuevo más adelante, con algunas modificaciones.

Ejecute el workflow para probarlo:

```bash
nextflow run main.nf
```

??? failure "Salida del comando"

    ```console hl_lines="7-8"
     N E X T F L O W   ~  version 25.10.2

    Launching `main.nf` [angry_koch] DSL2 - revision: 44fdf66105

    [-        ] COUNT_LINES -
    [-        ] COUNT_LINES -
    [patientA_rep1_normal_R, [/workspaces/training/side-quests/working_with_files/data/patientA_rep1_normal_R1_001.fastq.gz, /workspaces/training/side-quests/working_with_files/data/patientA_rep1_normal_R2_001.fastq.gz]]
    ERROR ~ Error executing process > 'COUNT_LINES (1)'

    Caused by:
      Not a valid path value: 'patientA_rep1_normal_R'



    Tip: when you have fixed the problem you can continue the execution adding the option `-resume` to the run command line

    -- Check '.nextflow.log' file for details
    ```

¡Vaya, esta vez la ejecución falló!

La parte relevante del mensaje de error está aquí:

```console
Not a valid path value: 'patientA_rep1_normal_R'
```

Eso se debe a que hemos cambiado la fábrica de canales.
Hasta ahora, el canal de entrada original solo contenía las rutas de archivo.
Toda la manipulación de metadatos que hemos estado haciendo en realidad no afectó el contenido del canal.

Ahora que estamos usando la fábrica de canales `.fromFilePairs`, el contenido del canal resultante es diferente.
Vemos solo un elemento del canal, compuesto por una tupla que contiene dos elementos: la parte del `simpleName` compartida por los dos archivos, que sirve como identificador, y una tupla que contiene los dos objetos de archivo, en el formato `id, [ file1, file2 ]`.

Eso es excelente, porque Nextflow ha hecho el trabajo difícil de extraer el nombre del paciente examinando el prefijo compartido y usándolo como identificador del paciente.

Sin embargo, sí rompe nuestro workflow actual.
Si quisiéramos seguir ejecutando `COUNT_LINES` de la misma manera sin cambiar el proceso, tendríamos que aplicar una operación de mapeo para extraer las rutas de archivo.
Pero no vamos a hacer eso, porque nuestro objetivo final es usar un proceso diferente, `ANALYZE_READS`, que maneja los pares de archivos apropiadamente.

Así que simplemente comentemos (o eliminemos) la llamada a `COUNT_LINES` y continuemos.

=== "Después"

    ```groovy title="main.nf" linenums="26" hl_lines="2"
        // Contar las líneas del archivo
        // COUNT_LINES(ch_files)
    ```

=== "Antes"

    ```groovy title="main.nf" linenums="26" hl_lines="2"
        // Contar las líneas del archivo
        COUNT_LINES(ch_files)
    ```

También puede comentar o eliminar la declaración include de `COUNT_LINES`, pero eso no tendrá ningún efecto funcional.

Ahora ejecutemos el workflow nuevamente:

```bash
nextflow run main.nf
```

??? success "Salida del comando"

    ```console hl_lines="5"
     N E X T F L O W   ~  version 25.10.2

    Launching `main.nf` [fabulous_davinci] DSL2 - revision: 22b53268dc

    [patientA_rep1_normal_R, [/workspaces/training/side-quests/working_with_files/data/patientA_rep1_normal_R1_001.fastq.gz, /workspaces/training/side-quests/working_with_files/data/patientA_rep1_normal_R2_001.fastq.gz]]
    ```

¡Esta vez el workflow se ejecuta con éxito!

Sin embargo, aún necesitamos extraer el resto de los metadatos del campo `id`.

### 5.2. Extraer y organizar metadatos de los pares de archivos

Nuestra operación `map` anterior no funcionará porque no coincide con la estructura de datos, pero podemos modificarla para que funcione.

Ya tenemos acceso al identificador real del paciente en la cadena que `fromFilePairs()` usó como identificador, por lo que podemos usarla para extraer los metadatos sin obtener el `simpleName` del objeto Path como hicimos antes.

Descomente la operación map en el workflow y realice los siguientes cambios:

=== "Después"

    ```groovy title="main.nf" linenums="7" hl_lines="3-4 9 11 13"
        // Cargar archivos con channel.fromFilePairs
        ch_files = channel.fromFilePairs('data/patientA_rep1_normal_R{1,2}_001.fastq.gz')
        ch_files.map { id, files ->
            def (sample, replicate, type) = id.tokenize('_')
            [
                [
                    id: sample,
                    replicate: replicate.replace('rep', ''),
                    type: type
                ],
                files
            ]
        }
        .view()
    ```

=== "Antes"

    ```groovy title="main.nf" linenums="7" hl_lines="3-5 11 13"
        // Cargar archivos con channel.fromFilePairs
        ch_files = channel.fromFilePairs('data/patientA_rep1_normal_R{1,2}_001.fastq.gz')
        /* ¡Comentamos el mapeo por ahora, volveremos a él!
        ch_files.map { myFile ->
            def (sample, replicate, type, readNum) = myFile.simpleName.tokenize('_')
            [
                [
                    id: sample,
                    replicate: replicate.replace('rep', ''),
                    type: type,
                    readNum: readNum,
                ],
                myFile
            ]
        }
        */
        .view()
    ```

Esta vez el map comienza desde `id, files` en lugar de solo `myFile`, y `tokenize()` se aplica a `id` en lugar de a `myFile.simpleName`.

Observe también que hemos eliminado `readNum` de la línea `tokenize()`; cualquier subcadena que no nombremos específicamente (comenzando desde la izquierda) se descartará silenciosamente.
Podemos hacer esto porque los archivos pareados ahora están estrechamente asociados, por lo que ya no necesitamos `readNum` en el map de metadatos.

Ejecutemos el workflow:

```bash
nextflow run main.nf
```

??? success "Salida del comando"

    ```console

    N E X T F L O W   ~  version 25.10.2

    Launching `main.nf` [prickly_stonebraker] DSL2 - revision: f62ab10a3f

    [[id:patientA, replicate:1, type:normal], [/workspaces/training/side-quests/working_with_files/data/patientA_rep1_normal_R1_001.fastq.gz, /workspaces/training/side-quests/working_with_files/data/patientA_rep1_normal_R2_001.fastq.gz]]
    ```

Y ahí está: tenemos el map de metadatos (`[id:patientA, replicate:1, type:normal]`) en la primera posición de la tupla de salida, seguido de la tupla de archivos pareados, como se pretendía.

Por supuesto, esto solo recogerá y procesará ese par específico de archivos.
Si desea experimentar con el procesamiento de múltiples pares, puede intentar agregar comodines al patrón de entrada y ver qué sucede.
Por ejemplo, pruebe usando `data/patientA_rep1_*_R{1,2}_001.fastq.gz`

### Conclusión

- [`channel.fromFilePairs()` encuentra y empareja automáticamente archivos relacionados](https://www.nextflow.io/docs/latest/reference/channel.html#fromfilepairs)
- Esto simplifica el manejo de lecturas de extremos pareados en su pipeline
- Los archivos pareados pueden agruparse como tuplas `[id, [file1, file2]]`
- La extracción de metadatos puede hacerse desde el ID del par de archivos en lugar de archivos individuales

---

## 6. Uso de operaciones de archivos en procesos

Ahora pongamos todo esto junto en un proceso simple para reforzar cómo usar operaciones de archivos dentro de un proceso de Nextflow.

Le proporcionamos un módulo de proceso preescrito llamado `ANALYZE_READS` que toma una tupla de metadatos y un par de archivos de entrada y los analiza.
Podríamos imaginar que esto está realizando alineamiento de secuencias, llamada de variantes o cualquier otro paso que tenga sentido para este tipo de datos.

Comencemos.

### 6.1. Importar el proceso y examinar el código

Para usar este proceso en el workflow, solo necesitamos agregar una declaración include del módulo antes del bloque workflow.

Realice el siguiente cambio en el workflow:

=== "Después"

    ```groovy title="main.nf" linenums="1" hl_lines="3"
    #!/usr/bin/env nextflow

    include { ANALYZE_READS } from './modules/analyze_reads.nf'

    workflow {
    ```

=== "Antes"

    ```groovy title="main.nf" linenums="1"
    #!/usr/bin/env nextflow

    workflow {
    ```

Puede abrir el archivo del módulo para examinar su código:

```groovy title="modules/analyze_reads.nf - process example" linenums="1"
#!/usr/bin/env nextflow

process ANALYZE_READS {
    tag { meta.id }

    publishDir { "results/${meta.id}" }, mode: 'copy'

    input:
    tuple val(meta), path(files)

    output:
    tuple val(meta.id), path("${meta.id}_stats.txt")

    script:
    """
    echo "Sample metadata: ${meta.id}" > ${meta.id}_stats.txt
    echo "Replicate: ${meta.replicate}" >> ${meta.id}_stats.txt
    echo "Type: ${meta.type}" >> ${meta.id}_stats.txt
    echo "Read 1: ${files[0]}" >> ${meta.id}_stats.txt
    echo "Read 2: ${files[1]}" >> ${meta.id}_stats.txt
    echo "Read 1 size: \$(gunzip -dc ${files[0]} | wc -l | awk '{print \$1/4}') reads" >> ${meta.id}_stats.txt
    echo "Read 2 size: \$(gunzip -dc ${files[1]} | wc -l | awk '{print \$1/4}') reads" >> ${meta.id}_stats.txt
    """
}
```

!!! note "Nota"

    Las directivas `tag` y `publishDir` usan sintaxis de closure (`{ ... }`) en lugar de interpolación de cadenas (`"${...}"`).
    Esto se debe a que estas directivas hacen referencia a variables de entrada (`meta`) que no están disponibles hasta el tiempo de ejecución.
    La sintaxis de closure difiere la evaluación hasta que el proceso realmente se ejecuta.

!!! note "Nota"

    Llamamos a nuestro map de metadatos `meta` por convención.
    Para una exploración más profunda de los meta maps, consulte la misión secundaria [Metadata and meta maps](../metadata/).

### 6.2. Llamar al proceso en el workflow

Ahora que el proceso está disponible para el workflow, podemos agregar una llamada al proceso `ANALYZE_READS` para ejecutarlo.

Para ejecutarlo sobre nuestros datos de ejemplo, necesitaremos hacer dos cosas:

1. Dar un nombre al canal remapeado
2. Agregar una llamada al proceso

#### 6.2.1. Nombrar el canal de entrada remapeado

Anteriormente aplicamos las manipulaciones de mapeo directamente al canal de entrada.
Para alimentar el contenido remapeado al proceso `ANALYZE_READS` (y hacerlo de una manera clara y fácil de leer) queremos crear un nuevo canal llamado `ch_samples`.

Podemos hacer eso usando el operador [`set`](https://www.nextflow.io/docs/latest/reference/operator.html#set).

En el workflow principal, reemplace el operador `.view()` con `.set { ch_samples }`, y agregue una línea para verificar que podemos referirnos al canal por nombre.

=== "Después"

    ```groovy title="main.nf" linenums="7" hl_lines="14 16-17"
        // Cargar archivos con channel.fromFilePairs
        ch_files = channel.fromFilePairs('data/patientA_rep1_normal_R{1,2}_001.fastq.gz')
        ch_files.map { id,  files ->
           def (sample, replicate, type, readNum) = id.tokenize('_')
           [
               [
                   id: sample,
                   replicate: replicate.replace('rep', ''),
                   type: type
               ],
               files
           ]
        }
            .set { ch_samples }

        // Temporal: inspeccionar ch_samples
        ch_samples.view()
    ```

=== "Antes"

    ```groovy title="main.nf" linenums="7" hl_lines="14"
        // Cargar archivos con channel.fromFilePairs
        ch_files = channel.fromFilePairs('data/patientA_rep1_normal_R{1,2}_001.fastq.gz')
        ch_files.map { id,  files ->
           def (sample, replicate, type, readNum) = id.tokenize('_')
           [
               [
                   id: sample,
                   replicate: replicate.replace('rep', ''),
                   type: type
               ],
               files
           ]
        }
        .view()
    }
    ```

Ejecutemos esto:

```bash
nextflow run main.nf
```

??? success "Salida del comando"

    ```console
     N E X T F L O W   ~  version 25.10.2

    Launching `main.nf` [goofy_kirch] DSL2 - revision: 3313283e42

    [[id:patientA, replicate:1, type:normal], [/workspaces/training/side-quests/working_with_files/data/patientA_rep1_normal_R1_001.fastq.gz, /workspaces/training/side-quests/working_with_files/data/patientA_rep1_normal_R2_001.fastq.gz]]
    ```

Esto confirma que ahora podemos referirnos al canal por nombre.

#### 6.2.2. Llamar al proceso sobre los datos

Ahora llamemos realmente al proceso `ANALYZE_READS` sobre el canal `ch_samples`.

En el workflow principal, realice los siguientes cambios en el código:

=== "Después"

    ```groovy title="main.nf" linenums="23"
        // Ejecutar el análisis
        ANALYZE_READS(ch_samples)
    ```

=== "Antes"

    ```groovy title="main.nf" linenums="23"
        // Temporal: inspeccionar ch_samples
        ch_samples.view()
    ```

Ejecutemos esto:

```bash
nextflow run main.nf
```

??? success "Salida del comando"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `./main.nf` [shrivelled_cori] DSL2 - revision: b546a31769

    executor >  local (1)
    [b5/110360] process > ANALYZE_READS (patientA) [100%] 1 of 1 ✔
    ```

Este proceso está configurado para publicar sus salidas en un directorio `results`, así que eche un vistazo allí.

??? abstract "Contenido de directorio y archivo"

    ```console
    results
    └── patientA
        └── patientA_stats.txt
    ```

    ```txt title="patientA_stats.txt"
    Sample metadata: patientA
    Replicate: 1
    Type: normal
    Read 1: patientA_rep1_normal_R1_001.fastq.gz
    Read 2: patientA_rep1_normal_R2_001.fastq.gz
    Read 1 size: 10 reads
    Read 2 size: 10 reads
    ```

El proceso tomó nuestras entradas y creó un nuevo archivo que contiene los metadatos del paciente, tal como fue diseñado.
¡Espléndido!

### 6.3. Incluir muchos más pacientes

Por supuesto, esto solo está procesando un par de archivos para un solo paciente, lo cual no es exactamente el tipo de alto rendimiento que espera obtener con Nextflow.
Probablemente querrá procesar muchos más datos a la vez.

Recuerde que `channel.fromPath()` acepta un _glob_ como entrada, lo que significa que puede aceptar cualquier número de archivos que coincidan con el patrón.
Por lo tanto, si queremos incluir a todos los pacientes, simplemente podemos modificar la cadena de entrada para incluir más pacientes, como se mencionó de pasada anteriormente.

Supongamos que queremos ser lo más inclusivos posible.
Realice los siguientes cambios en el workflow:

=== "Después"

    ```groovy title="main.nf" linenums="7" hl_lines="2"
        // Cargar archivos con channel.fromFilePairs
        ch_files = channel.fromFilePairs('data/*_R{1,2}_001.fastq.gz')
    ```

=== "Antes"

    ```groovy title="main.nf" linenums="7" hl_lines="2"
        // Cargar archivos con channel.fromFilePairs
        ch_files = channel.fromFilePairs('data/patientA_rep1_normal_R{1,2}_001.fastq.gz')
    ```

Ejecute el pipeline nuevamente:

```bash
nextflow run main.nf
```

??? success "Salida del comando"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `./main.nf` [big_stonebraker] DSL2 - revision: f7f9b8a76c

    executor >  local (8)
    [d5/441891] process > ANALYZE_READS (patientC) [100%] 8 of 8 ✔
    ```

El directorio de resultados ahora debería contener resultados para todos los datos disponibles.

??? abstract "Contenido del directorio"

    ```console
    results
    ├── patientA
    │   └── patientA_stats.txt
    ├── patientB
    │   └── patientB_stats.txt
    └── patientC
        └── patientC_stats.txt
    ```

¡Éxito! ¡Hemos analizado a todos los pacientes de una sola vez! ¿Verdad?

Quizás no.
Si mira más de cerca, tenemos un problema: tenemos dos réplicas para patientA, ¡pero solo un archivo de salida!
Estamos sobreescribiendo el archivo de salida cada vez.

### 6.4. Hacer únicos los archivos publicados

Dado que tenemos acceso a los metadatos del paciente, podemos usarlos para hacer únicos los archivos publicados incluyendo metadatos diferenciadores, ya sea en la estructura de directorios o en los propios nombres de archivo.

Realice el siguiente cambio en el workflow:

=== "Después"

    ```groovy title="modules/analyze_reads.nf" linenums="6"
        publishDir { "results/${meta.type}/${meta.id}/${meta.replicate}" }, mode: 'copy'
    ```

=== "Antes"

    ```groovy title="modules/analyze_reads.nf" linenums="6"
        publishDir { "results/${meta.id}" }, mode: 'copy'
    ```

Aquí mostramos la opción de usar niveles de directorio adicionales para dar cuenta de los tipos de muestra y las réplicas, pero también podría experimentar haciéndolo a nivel del nombre de archivo.

Ahora ejecute el pipeline una vez más, pero asegúrese de eliminar el directorio de resultados primero para tener un espacio de trabajo limpio:

```bash
rm -r results
nextflow run main.nf
```

??? success "Salida del comando"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `./main.nf` [insane_swartz] DSL2 - revision: fff18abe6d

    executor >  local (8)
    [e3/449081] process > ANALYZE_READS (patientC) [100%] 8 of 8 ✔
    ```

Revise el directorio de resultados ahora:

??? abstract "Contenido del directorio"

    ```console
    results/
    ├── normal
    │   ├── patientA
    │   │   ├── 1
    │   │   │   └── patientA_stats.txt
    │   │   └── 2
    │   │       └── patientA_stats.txt
    │   ├── patientB
    │   │   └── 1
    │   │       └── patientB_stats.txt
    │   └── patientC
    │       └── 1
    │           └── patientC_stats.txt
    └── tumor
        ├── patientA
        │   ├── 1
        │   │   └── patientA_stats.txt
        │   └── 2
        │       └── patientA_stats.txt
        ├── patientB
        │   └── 1
        │       └── patientB_stats.txt
        └── patientC
            └── 1
                └── patientC_stats.txt
    ```

Y ahí está, todos nuestros metadatos, perfectamente organizados. ¡Eso es un éxito!

Hay mucho más que puede hacer una vez que tiene sus metadatos cargados en un map como este:

1. Crear directorios de salida organizados basados en los atributos del paciente
2. Tomar decisiones en los procesos basadas en las propiedades del paciente
3. Dividir, unir y recombinar datos basándose en valores de metadatos

Este patrón de mantener los metadatos explícitos y adjuntos a los datos (en lugar de codificados en los nombres de archivo) es una práctica recomendada fundamental en Nextflow que permite construir workflows de análisis robustos y mantenibles.
Puede aprender más sobre esto en la misión secundaria [Metadata and meta maps](../metadata/).

### Conclusión

- La directiva `publishDir` puede organizar las salidas basándose en valores de metadatos
- Los metadatos en tuplas permiten la organización estructurada de los resultados
- Este enfoque crea workflows mantenibles con una clara procedencia de los datos
- Los procesos pueden tomar tuplas de metadatos y archivos como entrada
- La directiva `tag` proporciona identificación del proceso en los registros de ejecución
- La estructura del workflow separa la creación de canales de la ejecución de procesos

---

## Resumen

En esta misión secundaria, ha aprendido cómo trabajar con archivos en Nextflow, desde operaciones básicas hasta técnicas más avanzadas para manejar colecciones de archivos.

Aplicar estas técnicas en su propio trabajo le permitirá construir workflows más eficientes y mantenibles, especialmente cuando trabaje con grandes cantidades de archivos con convenciones de nomenclatura complejas.

### Patrones clave

1.  **Operaciones básicas con archivos:** Creamos objetos Path con `file()` y accedimos a atributos del archivo como nombre, extensión y directorio padre, aprendiendo la diferencia entre cadenas y objetos Path.

    - Crear un objeto Path con `file()`

    ```groovy
    myFile = file('path/to/file.txt')
    ```

    - Obtener atributos del archivo

    ```groovy
    println myFile.name       // file.txt
    println myFile.baseName   // file
    println myFile.extension  // txt
    println myFile.parent     // path/to
    ```

2.  **Uso de archivos remotos**: Aprendimos cómo cambiar de manera transparente entre archivos locales y remotos usando URIs, demostrando la capacidad de Nextflow para manejar archivos de diversas fuentes sin cambiar la lógica del workflow.

    - Archivo local

    ```groovy
    myFile = file('path/to/file.txt')
    ```

    - FTP

    ```groovy
    myFile = file('ftp://path/to/file.txt')
    ```

    - HTTPS

    ```groovy
    myFile = file('https://path/to/file.txt')
    ```

    - Amazon S3

    ```groovy
    myFile = file('s3://path/to/file.txt')
    ```

    - Azure Blob Storage

    ```groovy
    myFile = file('az://path/to/file.txt')
    ```

    - Google Cloud Storage

    ```groovy
    myFile = file('gs://path/to/file.txt')
    ```

3.  **Carga de archivos usando la fábrica de canales `fromPath()`:** Creamos canales a partir de patrones de archivos con `channel.fromPath()` y visualizamos sus atributos de archivo, incluyendo los tipos de objetos.

    - Crear un canal a partir de un patrón de archivos

    ```groovy
     ch_files = channel.fromPath('data/*.fastq.gz')
    ```

    - Obtener atributos del archivo

    ```groovy
     ch_files.view { myFile ->
        println "File object class: ${myFile.class}"
        println "File name: ${myFile.name}"
        println "Simple name: ${myFile.simpleName}"
        println "Extension: ${myFile.extension}"
        println "Parent directory: ${myFile.parent}"
    }
    ```

4.  **Extracción de metadatos del paciente de los nombres de archivo:** Usamos `tokenize()` y `replace()` para extraer y estructurar metadatos de los nombres de archivo, convirtiéndolos en maps organizados.

    ```groovy
    def name = file.name.tokenize('_')
    def patientId = name[0]
    def replicate = name[1].replace('rep', '')
    def type = name[2]
    def readNum = name[3].replace('R', '')
    ```

5.  **Simplificación con channel.fromFilePairs:** Usamos `channel.fromFilePairs()` para emparejar automáticamente archivos relacionados y extraer metadatos de los IDs de los pares de archivos.

    ```groovy
    ch_pairs = channel.fromFilePairs('data/*_R{1,2}_001.fastq.gz')
    ```

6.  **Uso de operaciones de archivos en procesos:** Integramos operaciones de archivos en procesos de Nextflow con el manejo adecuado de entradas, usando `publishDir` para organizar las salidas basándose en metadatos.

    - Asociar un meta map con las entradas del proceso

    ```groovy
    ch_files = channel.fromFilePairs('data/patientA_rep1_normal_R{1,2}_001.fastq.gz')
    ch_files.map { id,  files ->
        def (sample, replicate, type, readNum) = id.tokenize('_')
        [
            [
                id: sample,
                replicate: replicate.replace('rep', ''),
                type: type
            ],
             files
        ]
    }
        .set { ch_samples }

    ANALYZE_READS(ch_samples)
    ```

    - Organizar las salidas basándose en metadatos

    ```groovy
    publishDir { "results/${meta.type}/${meta.id}/${meta.replicate}" }, mode: 'copy'
    ```

### Recursos adicionales

- [Documentación de Nextflow: Trabajando con archivos](https://www.nextflow.io/docs/latest/working-with-files.html)
- [channel.fromPath](https://www.nextflow.io/docs/latest/reference/channel.html#frompath)
- [channel.fromFilePairs](https://www.nextflow.io/docs/latest/reference/channel.html#fromfilepairs)

---

## ¿Qué sigue?

Regrese al [menú de misiones secundarias](../) o haga clic en el botón en la parte inferior derecha de la página para continuar con el siguiente tema de la lista.
