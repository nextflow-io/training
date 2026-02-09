# Procesamiento de archivos de entrada

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Traducción asistida por IA - [más información y sugerencias](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

Los flujos de trabajo de análisis científico a menudo implican el procesamiento de grandes cantidades de archivos.
Nextflow proporciona herramientas poderosas para manejar archivos de manera eficiente, ayudándole a organizar y procesar sus datos con un código mínimo.

### Objetivos de aprendizaje

En esta misión secundaria, exploraremos cómo Nextflow maneja archivos, desde operaciones básicas hasta técnicas más avanzadas para trabajar con colecciones de archivos.
Aprenderá cómo extraer metadata de los nombres de archivos, que es un requisito común en los pipelines de análisis científico.

Al final de esta misión secundaria, podrá:

- Crear objetos Path a partir de cadenas de rutas de archivo usando el método `file()` de Nextflow
- Acceder a atributos de archivos como nombre, extensión y directorio padre
- Manejar tanto archivos locales como remotos de manera transparente usando URIs
- Usar channels para automatizar el manejo de archivos con `channel.fromPath()` y `channel.fromFilePairs()`
- Extraer y estructurar metadata de los nombres de archivos usando manipulación de cadenas
- Agrupar archivos relacionados usando coincidencia de patrones y expresiones glob
- Integrar operaciones de archivos en procesos de Nextflow con el manejo apropiado de entradas
- Organizar salidas de procesos usando estructuras de directorios basadas en metadata

Estas habilidades le ayudarán a construir flujos de trabajo que puedan manejar diferentes tipos de entradas de archivos con gran flexibilidad.

### Requisitos previos

Antes de emprender esta misión secundaria, debería:

- Haber completado el tutorial [Hello Nextflow](../../hello_nextflow/) o un curso equivalente para principiantes.
- Sentirse cómodo usando conceptos y mecanismos básicos de Nextflow (procesos, channels, operadores)

<!-- I removed the suggestion to do the metamaps SQ first because that works more naturally after -->

---

## 0. Primeros pasos

#### Abrir el codespace de entrenamiento

Si aún no lo ha hecho, asegúrese de abrir el entorno de entrenamiento como se describe en [Configuración del Entorno](../envsetup/index.md).

[![Open in GitHub Codespaces](https://github.com/codespaces/badge.svg)](https://codespaces.new/nextflow-io/training?quickstart=1&ref=master)

#### Moverse al directorio del proyecto

Vamos a movernos al directorio donde se encuentran los archivos para este tutorial.

```bash
cd side-quests/working_with_files
```

Puede configurar VSCode para enfocarse en este directorio:

```bash
code .
```

#### Revisar los materiales

Encontrará un archivo de workflow simple llamado `main.nf`, un directorio `modules` que contiene dos archivos de módulos, y un directorio `data` que contiene algunos archivos de datos de ejemplo.

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

Este directorio contiene datos de secuenciación de extremos emparejados de tres pacientes (A, B, C).

Para cada paciente, tenemos muestras que son de tipo `tumor` (típicamente originadas de biopsias de tumores) o `normal` (tomadas de tejido sano o sangre).
Si no está familiarizado con el análisis de cáncer, solo sepa que esto corresponde a un modelo experimental que usa muestras emparejadas tumor/normal para realizar análisis contrastivos.

Para el paciente A específicamente, tenemos dos conjuntos de réplicas técnicas (repeticiones).

Los archivos de datos de secuenciación se nombran con una convención típica `_R1_` y `_R2_` para lo que se conoce como 'lecturas forward' y 'lecturas reverse'.

_No se preocupe si no está familiarizado con este diseño experimental, no es crítico para entender este tutorial._

#### Revisar la asignación

Su desafío es escribir un workflow de Nextflow que:

1. **Cargue** archivos de entrada usando los métodos de manejo de archivos de Nextflow
2. **Extraiga** metadata (ID del paciente, réplica, tipo de muestra) de la estructura del nombre de archivo
3. **Agrupe** archivos emparejados (R1/R2) juntos usando `channel.fromFilePairs()`
4. **Procese** los archivos con un módulo de análisis proporcionado
5. **Organice** las salidas en una estructura de directorios basada en la metadata extraída

#### Lista de verificación de preparación

¿Cree que está listo para comenzar?

- [ ] Entiendo el objetivo de este curso y sus requisitos previos
- [ ] Mi codespace está funcionando
- [ ] He establecido mi directorio de trabajo apropiadamente
- [ ] Entiendo la asignación

Si puede marcar todas las casillas, está listo para comenzar.

---

## 1. Operaciones básicas de archivos

### 1.1. Identificar el tipo de un objeto con `.class`

Eche un vistazo al archivo de workflow `main.nf`:

```groovy title="main.nf" linenums="1"
#!/usr/bin/env nextflow

workflow {

    // Crear un objeto Path a partir de una ruta de texto
    myFile = 'data/patientA_rep1_normal_R1_001.fastq.gz'

    println "${myFile} is of class ${myFile.class}"
}
```

Este es un mini-workflow (sin procesos) que hace referencia a una única ruta de archivo en su workflow, luego la imprime en la consola, junto con su clase.

??? info "¿Qué es `.class`?"

    En Nextflow, `.class` nos dice qué tipo de objeto estamos manejando. Es como preguntar "¿qué tipo de cosa es esto?" para averiguar si es una cadena, un número, un archivo o algo más.
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

Como puede ver, Nextflow imprimió la ruta de cadena exactamente como la escribimos.

Esto es solo salida de texto; Nextflow aún no ha hecho nada especial con ella.
También hemos confirmado que, hasta donde Nextflow tiene conocimiento, esto es solo una cadena (de clase `java.lang.String`).
Eso tiene sentido, ya que aún no le hemos dicho a Nextflow que corresponde a un archivo.

### 1.2. Crear un objeto Path con file()

Podemos indicarle a Nextflow cómo manejar archivos creando [objetos Path](https://www.nextflow.io/docs/latest/reference/stdlib-types.html#path) a partir de cadenas de rutas.

En nuestro workflow, podemos convertir la cadena de ruta `data/patientA_rep1_normal_R1_001.fastq.gz` a un objeto Path usando el método `file()`, que proporciona acceso a propiedades y operaciones de archivos.

Edite el `main.nf` para envolver la cadena con `file()` de la siguiente manera:

=== "Después"

    ```groovy title="main.nf" linenums="5" hl_lines="2"
        // Crear un objeto Path a partir de una ruta de texto
        myFile = file('data/patientA_rep1_normal_R1_001.fastq.gz')

        println "${myFile} is of class ${myFile.class}"
    ```

=== "Antes"

    ```groovy title="main.nf" linenums="5" hl_lines="2"
        // Crear un objeto Path a partir de una ruta de texto
        myFile = 'data/patientA_rep1_normal_R1_001.fastq.gz'

        println "${myFile} is of class ${myFile.class}"
    ```

Ahora ejecute el workflow de nuevo:

```bash
nextflow run main.nf
```

??? success "Salida del comando"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `main.nf` [kickass_coulomb] DSL2 - revision: 5af44b1b59

    /workspaces/training/side-quests/working_with_files/data/patientA_rep1_normal_R1_001.fastq.gz is of class class sun.nio.fs.UnixPath
    ```

Esta vez, ve la ruta absoluta completa en lugar de la ruta relativa que proporcionamos como entrada.

Nextflow ha convertido nuestra cadena en un objeto Path y la ha resuelto a la ubicación real del archivo en el sistema.
La ruta del archivo ahora será absoluta, como en `/workspaces/training/side-quests/working_with_files/data/patientA_rep1_normal_R1_001.fastq.gz`.

Note también que la clase del objeto Path es `sun.nio.fs.UnixPath`: esta es la forma de Nextflow de representar archivos locales.
Como veremos más adelante, los archivos remotos tendrán nombres de clase diferentes (como `nextflow.file.http.XPath` para archivos HTTP), pero todos funcionan exactamente de la misma manera y pueden usarse de manera idéntica en sus flujos de trabajo.

!!! tip

    **La diferencia clave:**

    - **Cadena de ruta**: Solo texto que Nextflow trata como caracteres
    - **Objeto Path**: Una referencia de archivo inteligente con la que Nextflow puede trabajar

    Piénselo así: una cadena de ruta es como escribir una dirección en papel, mientras que un objeto Path es como tener la dirección cargada en un dispositivo GPS que sabe cómo navegar hasta allí y puede decirle detalles sobre el viaje.

### 1.3. Acceder a atributos de archivo

¿Por qué es esto útil? Bueno, ahora que Nextflow entiende que `myFile` es un objeto Path y no solo una cadena, podemos acceder a los diversos atributos del objeto Path.

Actualicemos nuestro workflow para imprimir los atributos de archivo incorporados:

=== "Después"

    ```groovy title="main.nf" linenums="5" hl_lines="4-9"
        // Crear un objeto Path a partir de una ruta de texto
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
        // Crear un objeto Path a partir de una ruta de texto
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

Ve los diversos atributos de archivo impresos en la consola arriba.

### 1.4. Alimentar el archivo a un proceso

La diferencia entre cadenas y objetos Path se vuelve crítica cuando comienza a construir workflows reales con procesos.
Hasta ahora hemos verificado que Nextflow ahora está tratando nuestro archivo de entrada como un archivo, pero veamos si podemos realmente ejecutar algo en ese archivo en un proceso.

#### 1.4.1. Importar el proceso y examinar el código

Le proporcionamos un módulo de proceso pre-escrito llamado `COUNT_LINES` que toma una entrada de archivo y cuenta cuántas líneas contiene.

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

Puede abrir el archivo de módulo para examinar su código:

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

Como puede ver, es un pequeño script bastante directo que descomprime el archivo y cuenta cuántas líneas contiene.

??? info "¿Qué hace `debug true`?"

    La directiva `debug true` en la definición del proceso hace que Nextflow imprima la salida de su script (como el conteo de líneas "40") directamente en el registro de ejecución.
    Sin esto, solo vería el estado de ejecución del proceso pero no la salida real de su script.

    Para más información sobre depuración de procesos de Nextflow, vea la misión secundaria [Debugging Nextflow Workflows](debugging.md).

#### 1.4.2. Agregar una llamada a `COUNT_LINES`

Ahora que el proceso está disponible para el workflow, podemos agregar una llamada al proceso `COUNT_LINES` para ejecutarlo en el archivo de entrada.

Haga las siguientes ediciones al workflow:

=== "Después"

    ```groovy title="main.nf" linenums="7" hl_lines="11-12"
        // Crear un objeto Path a partir de una ruta de texto
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
        // Crear un objeto Path a partir de una ruta de texto
        myFile = file('data/patientA_rep1_normal_R1_001.fastq.gz')

        // Imprimir atributos del archivo
        println "File object class: ${myFile.class}"
        println "File name: ${myFile.name}"
        println "Simple name: ${myFile.simpleName}"
        println "Extension: ${myFile.extension}"
        println "Parent directory: ${myFile.parent}"
    ```

Y ahora ejecute el workflow:

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

Esto muestra que somos capaces de operar en el archivo apropiadamente dentro de un proceso.

Específicamente, Nextflow llevó a cabo las siguientes operaciones exitosamente:

- Organizó el archivo en el directorio de trabajo
- Descomprimió el archivo .gz
- Contó las líneas (40 líneas en este caso)
- Completó sin error

La clave de esta operación fluida es que estamos diciéndole explícitamente a Nextflow que nuestra entrada es un archivo y debe ser tratada como tal.

### 1.5. Solucionar errores básicos de entrada de archivo

Esto a menudo confunde a los recién llegados a Nextflow, así que tomemos unos minutos para ver qué sucede cuando lo hace mal.

Hay dos lugares principales donde puede manejar incorrectamente los archivos: al nivel del workflow y al nivel del proceso.

#### 1.5.1. Error a nivel de workflow

Veamos qué sucede si volvemos a tratar el archivo como una cadena cuando especificamos la entrada en el bloque workflow.

Haga las siguientes ediciones al workflow, asegurándose de comentar las declaraciones print específicas de path:

=== "Después"

    ```groovy title="main.nf" linenums="7" hl_lines="2 6-11"
        // Crear un objeto Path a partir de una ruta de texto
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
        // Crear un objeto Path a partir de una ruta de texto
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

Y ahora ejecute el workflow:

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

Cuando especifica una entrada `path`, Nextflow valida que esté pasando referencias de archivo reales, no solo cadenas.
Este error le está diciendo que `'data/patientA_rep1_normal_R1_001.fastq.gz'` no es un valor de ruta válido porque es una cadena, no un objeto Path.

Nextflow detectó inmediatamente el problema y se detuvo antes de siquiera iniciar el proceso.

#### 1.5.2. Error a nivel de proceso

El otro lugar donde podríamos olvidar especificar que queremos que Nextflow trate la entrada como un archivo es en la definición del proceso.

!!! warning "Mantenga el error del workflow de 1.5.1"

    Para que esta prueba funcione correctamente, mantenga el workflow en su estado roto (usando una cadena simple en lugar de `file()`).
    Cuando se combina con `val` en el proceso, esto produce el error que se muestra a continuación.

Haga la siguiente edición al módulo:

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

Y ahora ejecute el workflow de nuevo:

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

Esto muestra muchos detalles sobre el error porque el proceso está configurado para generar información de depuración, como se señaló anteriormente.

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

Esto dice que el sistema no pudo encontrar el archivo; sin embargo, si busca la ruta, hay un archivo con ese nombre en esa ubicación.

Cuando ejecutamos esto, Nextflow pasó el valor de cadena al script, pero no _organizó_ el archivo real en el directorio de trabajo.
Entonces el proceso intentó usar la cadena relativa, `data/patientA_rep1_normal_R1_001.fastq.gz`, pero ese archivo no existe dentro del directorio de trabajo del proceso.

En conjunto, estos dos ejemplos le muestran cuán importante es decirle a Nextflow si una entrada debe ser manejada como un archivo.

!!! note

    Asegúrese de volver atrás y corregir ambos errores intencionales antes de continuar a la siguiente sección.

### Conclusión

- Cadenas de ruta vs objetos Path: Las cadenas son solo texto, los objetos Path son referencias de archivo inteligentes
- El método `file()` convierte una cadena de ruta en un objeto Path con el que Nextflow puede trabajar
- Puede acceder a propiedades de archivo como `name`, `simpleName`, `extension` y `parent` [usando atributos de archivo](https://www.nextflow.io/docs/latest/working-with-files.html#getting-file-attributes)
- Usar objetos Path en lugar de cadenas permite a Nextflow gestionar apropiadamente los archivos en su workflow
- Manejo de Entradas de Proceso: El manejo apropiado de archivos requiere objetos Path, no cadenas, para asegurar que los archivos se organicen correctamente y sean accesibles para su uso por los procesos.

---

## 2. Usar archivos remotos

Una de las características clave de Nextflow es la capacidad de cambiar sin problemas entre archivos locales (en la misma máquina) y archivos remotos accesibles a través de internet.

Si lo está haciendo correctamente, nunca debería necesitar cambiar la lógica de su workflow para acomodar archivos provenientes de diferentes ubicaciones.
Todo lo que necesita hacer para usar un archivo remoto es especificar el prefijo apropiado en la ruta del archivo cuando lo suministra al workflow.

Por ejemplo, `/path/to/data` no tiene prefijo, indicando que es una ruta de archivo local 'normal', mientras que `s3://path/to/data` incluye el prefijo `s3://`, indicando que está ubicado en el almacenamiento de objetos S3 de Amazon.

Muchos protocolos diferentes son soportados:

- HTTP(S)/FTP (http://, https://, ftp://)
- Amazon S3 (s3://)
- Azure Blob Storage (az://)
- Google Cloud Storage (gs://)

Para usar cualquiera de estos, simplemente especifique el prefijo relevante en la cadena, que entonces técnicamente se llama Identificador Uniforme de Recursos (URI) en lugar de ruta de archivo.
Nextflow manejará la autenticación y la organización de los archivos en el lugar correcto, descargando o subiendo y todas las demás operaciones de archivos que esperaría.

La fortaleza clave de este sistema es que nos permite cambiar entre entornos sin cambiar ninguna lógica del pipeline.
Por ejemplo, puede desarrollar con un conjunto de prueba pequeño y local antes de cambiar a un conjunto de prueba a gran escala ubicado en almacenamiento remoto simplemente cambiando el URI.

### 2.1. Usar un archivo de internet

Probemos esto cambiando la ruta local que estamos proporcionando a nuestro workflow con una ruta HTTPS que apunta a una copia de los mismos datos que está almacenada en Github.

!!! warning

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
        // Crear un objeto Path a partir de una ruta de texto
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

La única diferencia en la salida de la consola es que la clase del objeto path ahora es `nextflow.file.http.XPath`, mientras que para la ruta local la clase era `sun.nio.fs.UnixPath`.
No necesita recordar estas clases; solo mencionamos esto para demostrar que Nextflow identifica y maneja las diferentes ubicaciones apropiadamente.

Entre bastidores, Nextflow descargó el archivo a un directorio de organización ubicado dentro del directorio de trabajo.
Ese archivo organizado puede entonces ser tratado como un archivo local y vinculado simbólicamente en el directorio del proceso relevante.

Puede verificar que eso sucedió aquí mirando el contenido del directorio de trabajo ubicado en el valor hash del proceso.

??? abstract "Contenido del directorio de trabajo"

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

    El enlace simbólico apunta a una copia organizada del archivo remoto que Nextflow descargó automáticamente.

Note que para archivos más grandes, el paso de descarga tomará tiempo adicional en comparación con ejecutar en archivos locales.
Sin embargo, Nextflow verifica si ya tiene una copia organizada para evitar descargas innecesarias.
Entonces, si ejecuta nuevamente en el mismo archivo y no ha eliminado el archivo organizado, Nextflow usará la copia organizada.

Esto muestra cuán fácil es cambiar entre datos locales y remotos usando Nextflow, que es una característica clave de Nextflow.

!!! note

    La única excepción importante a este principio es que no puede usar patrones glob o rutas de directorio con HTTPS porque HTTPS no puede listar múltiples archivos, por lo que debe especificar URLs de archivos exactas.
    Sin embargo, otros protocolos de almacenamiento como blob storage (`s3://`, `az://`, `gs://`) pueden usar tanto globs como rutas de directorio.

    Aquí está cómo podría usar patrones glob con almacenamiento en la nube:

    ```groovy title="Ejemplos de almacenamiento en la nube (no ejecutables en este entorno)"
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
        // Crear un objeto Path a partir de una ruta de texto
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
        // Usar un archivo remoto de internet
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
- Nextflow descargará y organizará automáticamente los datos en el lugar correcto, siempre que estas rutas se alimenten a los procesos
- ¡No escriba lógica para descargar o subir archivos remotos!
- Los archivos locales y remotos producen diferentes tipos de objetos pero funcionan de manera idéntica
- **Importante**: HTTP/HTTPS solo funcionan con archivos individuales (sin patrones glob)
- El almacenamiento en la nube (S3, Azure, GCS) soporta tanto archivos individuales como patrones glob
- Puede cambiar sin problemas entre fuentes de datos locales y remotas sin cambiar la lógica del código (siempre que el protocolo soporte sus operaciones requeridas)

---

## 3. Usar el channel factory `fromPath()`

Hasta ahora hemos estado trabajando con un solo archivo a la vez, pero en Nextflow, típicamente vamos a querer crear un canal de entrada con múltiples archivos de entrada para procesar.

Una forma ingenua de hacer eso sería combinar el método `file()` con [`channel.of()`](https://www.nextflow.io/docs/latest/reference/channel.html#of) así:

```groovy title="Ejemplo de sintaxis"
ch_files = channel.of([file('data/patientA_rep1_normal_R1_001.fastq.gz')],
                      [file('data/patientA_rep1_normal_R1_001.fastq.gz')])
```

Eso funciona, pero es torpe.

!!! tip "Cuándo usar `file()` vs `channel.fromPath()`"

    - Use `file()` cuando necesite un único objeto Path para manipulación directa (verificar si existe un archivo, leer sus atributos, o pasar a una única invocación de proceso)
    - Use `channel.fromPath()` cuando necesite un canal que pueda contener múltiples archivos, especialmente con patrones glob, o cuando los archivos fluirán a través de múltiples procesos

Aquí es donde entra [`channel.fromPath()`](https://www.nextflow.io/docs/latest/reference/channel.html#frompath): un conveniente channel factory que agrupa toda la funcionalidad que necesitamos para generar un canal a partir de una o más cadenas de archivos estáticas así como patrones glob.

### 3.1. Agregar el channel factory

Actualicemos nuestro workflow para usar `channel.fromPath`.

=== "Después"

    ```groovy title="main.nf" linenums="7" hl_lines="1-3"
        // Cargar archivos con channel.fromPath
        ch_files = channel.fromPath('data/patientA_rep1_normal_R1_001.fastq.gz')
        ch_files.view { myFile -> "Found file: $myFile" }

        // Imprimir atributos del archivo
        /* Comentar esto por ahora, volveremos a ello más adelante!
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
        // Crear un objeto Path a partir de una ruta de texto
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

También hemos comentado el código que imprime los atributos por ahora, y agregamos una declaración `.view` para imprimir solo el nombre del archivo en su lugar.

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
Esto es similar a lo que `file()` habría hecho, excepto que ahora tenemos un canal en el que podemos cargar más archivos si queremos.

Usar `channel.fromPath()` es una forma conveniente de crear un nuevo canal poblado por una lista de archivos.

### 3.2. Ver atributos de archivos en el canal

En nuestro primer intento de usar el channel factory, simplificamos el código y solo imprimimos el nombre del archivo.

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

También estamos re-habilitando la llamada al proceso `COUNT_LINES` para verificar que el procesamiento de archivos todavía funciona correctamente con nuestro enfoque basado en canales.

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

Y ahí está, mismos resultados que antes pero ahora tenemos el archivo en un canal, así que podemos agregar más.

### 3.3. Usar un glob para coincidir con múltiples archivos

Hay varias formas en que podríamos cargar más archivos en el canal.
Aquí vamos a mostrarle cómo usar patrones glob, que son una forma conveniente de coincidir y recuperar nombres de archivos y directorios basados en caracteres comodín.
El proceso de coincidir estos patrones se llama "globbing" o "expansión de nombres de archivo".

!!! note

    Como se señaló anteriormente, Nextflow soporta globbing para gestionar archivos de entrada y salida en la mayoría de los casos, excepto con rutas de archivo HTTPS porque HTTPS no puede listar múltiples archivos.

Digamos que queremos recuperar ambos archivos en un par de archivos asociados con un paciente dado, `patientA`:

```console
patientA_rep1_normal_R1_001.fastq.gz
patientA_rep1_normal_R2_001.fastq.gz
```

Ya que la única diferencia entre los nombres de archivos es el número de réplica, _es decir_, el número después de `R`, podemos usar el carácter comodín `*` para representar el número de la siguiente manera:

```console
patientA_rep1_normal_R*_001.fastq.gz
```

Ese es el patrón glob que necesitamos.

Ahora todo lo que necesitamos hacer es actualizar la ruta del archivo en el channel factory para usar ese patrón glob de la siguiente manera:

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

Como puede ver, ahora tenemos dos objetos Path en nuestro canal, lo que muestra que Nextflow ha hecho la expansión de nombres de archivo correctamente, y ha cargado y procesado ambos archivos como se esperaba.

Usando este método, podemos recuperar tantos o tan pocos archivos como queramos simplemente cambiando el patrón glob. Si lo hiciéramos más generoso, por ejemplo reemplazando todas las partes variables de los nombres de archivo por `*` (_e.g._ `data/patient*_rep*_*_R*_001.fastq.gz`) podríamos tomar todos los archivos de ejemplo en el directorio `data`.

### Conclusión

- `channel.fromPath()` crea un canal con archivos que coinciden con un patrón
- Cada archivo se emite como un elemento separado en el canal
- Podemos usar un patrón glob para coincidir con múltiples archivos
- Los archivos se convierten automáticamente en objetos Path con atributos completos
- El método `.view()` permite inspección del contenido del canal

---

## 4. Extraer metadata básica de nombres de archivos

En la mayoría de dominios científicos, es muy común tener metadata codificada en los nombres de los archivos que contienen los datos.
Por ejemplo, en bioinformática, los archivos que contienen datos de secuenciación a menudo se nombran de una manera que codifica información sobre la muestra, condición, réplica y número de lectura.

Si los nombres de archivos se construyen de acuerdo con una convención consistente, puede extraer esa metadata de manera estandarizada y usarla en el curso de su análisis.
Ese es un gran 'si', por supuesto, y debe ser muy cauteloso cada vez que confíe en la estructura del nombre de archivo; pero la realidad es que este enfoque es muy ampliamente utilizado, así que veamos cómo se hace en Nextflow.

En el caso de nuestros datos de ejemplo, sabemos que los nombres de archivos incluyen metadata estructurada consistentemente.
Por ejemplo, el nombre de archivo `patientA_rep1_normal_R2_001` codifica lo siguiente:

- ID del paciente: `patientA`
- ID de réplica: `rep1`
- tipo de muestra: `normal` (en oposición a `tumor`)
- conjunto de lectura: `R1` (en oposición a `R2`)

Vamos a modificar nuestro workflow para recuperar esta información en tres pasos:

1. Recuperar el `simpleName` del archivo, que incluye la metadata
2. Separar la metadata usando un método llamado `tokenize()`
3. Usar un map para organizar la metadata

!!! warning

    Nunca debe codificar información sensible en nombres de archivos, como nombres de pacientes u otras características identificatorias, ya que eso puede comprometer la privacidad del paciente u otras restricciones de seguridad relevantes.

### 4.1. Recuperar el `simpleName`

El `simpleName` es un atributo de archivo que corresponde al nombre de archivo sin su ruta y extensión.

Haga las siguientes ediciones al workflow:

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

Ejecute el workflow para probar que funciona:

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

Cada elemento en el canal ahora es una tupla que contiene el `simpleName` y el objeto de archivo original.

### 4.2. Extraer la metadata del `simplename`

En este punto, la metadata que queremos está incrustada en el `simplename`, pero no podemos acceder a elementos individuales directamente.
Entonces necesitamos dividir el `simplename` en sus componentes.
Afortunadamente, esos componentes están simplemente separados por guiones bajos en el nombre de archivo original, así que podemos aplicar un método común de Nextflow llamado `tokenize()` que es perfecto para esta tarea.

Haga las siguientes ediciones al workflow:

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

Ahora la tupla para cada elemento en nuestro canal contiene la lista de metadata (_e.g._ `[patientA, rep1, normal, R1, 001]`) y el objeto de archivo original.

¡Eso es genial!
Hemos desglosado nuestra información de paciente de una sola cadena en una lista de cadenas.
Ahora podemos manejar cada parte de la información del paciente por separado.

### 4.3. Usar un map para organizar la metadata

Nuestra metadata es solo una lista plana en este momento.
Es suficientemente fácil de usar pero difícil de leer.

```console
[patientA, rep1, normal, R1, 001]
```

¿Qué es el elemento en el índice 3? ¿Puede decirlo sin referirse a la explicación original de la estructura de metadata?

Esta es una gran oportunidad para usar un almacén clave-valor, donde cada elemento tiene un conjunto de claves y sus valores asociados, para que pueda referirse fácilmente a cada clave para obtener el valor correspondiente.

En nuestro ejemplo, eso significa ir de esta organización:

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
Haga las siguientes ediciones al workflow:

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

- **Asignación por desestructuración**: `def (patient, replicate, type, readNum) = ...` extrae los valores tokenizados en variables nombradas en una línea
- **Sintaxis literal de map**: `[id: patient, replicate: ...]` crea un map donde cada clave (como `id`) está asociada con un valor (como `patient`)
- **Estructura anidada**: La lista externa `[..., myFile]` empareja el map de metadata con el objeto de archivo original

También simplificamos un par de cadenas de metadata usando un método de reemplazo de cadenas llamado `replace()` para eliminar algunos caracteres que son innecesarios (_e.g._ `replicate.replace('rep', '')` para mantener solo el número de los IDs de réplica).

Ejecutemos el workflow de nuevo:

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

Ahora la metadata está claramente etiquetada (_e.g._ `[id:patientA, replicate:1, type:normal, readNum:2]`) así que es mucho más fácil decir qué es qué.

También será mucho más fácil hacer uso de elementos de metadata en el workflow, y hará nuestro código más fácil de leer y más mantenible.

### Conclusión

- Podemos manejar nombres de archivos en Nextflow con el poder de un lenguaje de programación completo
- Podemos tratar los nombres de archivos como cadenas para extraer información relevante
- El uso de métodos como `tokenize()` y `replace()` nos permite manipular cadenas en el nombre de archivo
- La operación `.map()` transforma elementos del canal mientras preserva la estructura
- La metadata estructurada (maps) hace el código más legible y mantenible que las listas posicionales

A continuación, veremos cómo manejar archivos de datos emparejados.

---

## 5. Manejar archivos de datos emparejados

Muchos diseños experimentales producen archivos de datos emparejados que se benefician de ser manejados de una manera explícitamente emparejada.
Por ejemplo, en bioinformática, los datos de secuenciación a menudo se generan en forma de lecturas emparejadas, que significa cadenas de secuencia que se originan del mismo fragmento de ADN (a menudo llamadas 'forward' y 'reverse' porque se leen desde extremos opuestos).

Ese es el caso de nuestros datos de ejemplo, donde R1 y R2 se refieren a los dos conjuntos de lecturas.

```console
data/patientA_rep1_normal_R1_001.fastq.gz
data/patientA_rep1_normal_R2_001.fastq.gz
```

Nextflow proporciona un channel factory especializado para trabajar con archivos emparejados como este llamado `channel.fromFilePairs()`, que automáticamente agrupa archivos basándose en un patrón de nomenclatura compartido. Eso le permite asociar los archivos emparejados más estrechamente con menos esfuerzo.

Vamos a modificar nuestro workflow para aprovechar esto.
Va a tomar dos pasos:

1. Cambiar el channel factory a `channel.fromFilePairs()`
2. Extraer y mapear la metadata

### 5.1. Cambiar el channel factory a `channel.fromFilePairs()`

Para usar `channel.fromFilePairs`, necesitamos especificar el patrón que Nextflow debe usar para identificar los dos miembros en un par.

Volviendo a nuestros datos de ejemplo, podemos formalizar el patrón de nomenclatura de la siguiente manera:

```console
data/patientA_rep1_normal_R{1,2}_001.fastq.gz
```

Esto es similar al patrón glob que usamos anteriormente, excepto que específicamente enumera las subcadenas (ya sea `1` o `2` viniendo justo después de la R) que identifican los dos miembros del par.

Actualicemos el workflow `main.nf` en consecuencia:

=== "Después"

    ```groovy title="main.nf" linenums="7" hl_lines="1-2"
        // Cargar archivos con channel.fromFilePairs
        ch_files = channel.fromFilePairs('data/patientA_rep1_normal_R{1,2}_001.fastq.gz')
        /* Comentar el mapeo por ahora, volveremos a él!
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
        // Cargar archivos con channel.fromPath
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

Hemos cambiado el channel factory y adaptado el patrón de coincidencia de archivos, y mientras estábamos en ello, comentamos la operación map.
Agregaremos eso de vuelta más tarde, con algunas modificaciones.

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

Uh-oh, ¡esta vez la ejecución falló!

La parte relevante del mensaje de error está aquí:

```console
Not a valid path value: 'patientA_rep1_normal_R'
```

Eso es porque hemos cambiado el channel factory.
Hasta ahora, el canal de entrada original solo contenía las rutas de archivo.
Toda la manipulación de metadata que hemos estado haciendo no afectó realmente el contenido del canal.

Ahora que estamos usando el channel factory `.fromFilePairs`, el contenido del canal resultante es diferente.
Vemos solo un elemento de canal, compuesto de una tupla que contiene dos elementos: la parte del `simpleName` compartida por los dos archivos, que sirve como identificador, y una tupla que contiene los dos objetos de archivo, en el formato `id, [ file1, file2 ]`.

Eso es genial, porque Nextflow ha hecho el trabajo difícil de extraer el nombre del paciente examinando el prefijo compartido y usándolo como identificador del paciente.

Sin embargo, sí rompe nuestro workflow actual.
Si quisiéramos todavía ejecutar `COUNT_LINES` de la misma manera sin cambiar el proceso, tendríamos que aplicar una operación de mapeo para extraer las rutas de archivo.
Pero no vamos a hacer eso, porque nuestro objetivo final es usar un proceso diferente, `ANALYZE_READS`, que maneja pares de archivos apropiadamente.

Así que simplemente comentemos (o eliminemos) la llamada a `COUNT_LINES` y sigamos adelante.

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

También puede comentar o eliminar la declaración include de `COUNT_LINES`, pero eso no tendrá efecto funcional.

Ahora ejecutemos el workflow de nuevo:

```bash
nextflow run main.nf
```

??? success "Salida del comando"

    ```console hl_lines="5"
     N E X T F L O W   ~  version 25.10.2

    Launching `main.nf` [fabulous_davinci] DSL2 - revision: 22b53268dc

    [patientA_rep1_normal_R, [/workspaces/training/side-quests/working_with_files/data/patientA_rep1_normal_R1_001.fastq.gz, /workspaces/training/side-quests/working_with_files/data/patientA_rep1_normal_R2_001.fastq.gz]]
    ```

¡Yay, esta vez el workflow tiene éxito!

Sin embargo, todavía necesitamos obtener el resto de la metadata del campo `id`.

### 5.2. Extraer y organizar metadata de pares de archivos

Nuestra operación `map` de antes no funcionará porque no coincide con la estructura de datos, pero podemos modificarla para que funcione.

Ya tenemos acceso al identificador del paciente real en la cadena que `fromFilePairs()` usó como identificador, así que podemos usar eso para extraer la metadata sin obtener el `simpleName` del objeto Path como lo hicimos antes.

Descomente la operación map en el workflow y haga las siguientes ediciones:

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
        /* Comentar el mapeo por ahora, volveremos a él!
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

Note también que hemos eliminado `readNum` de la línea `tokenize()`; cualquier subcadena que no nombremos específicamente (comenzando desde la izquierda) será eliminada silenciosamente.
Podemos hacer esto porque los archivos emparejados ahora están estrechamente asociados, así que ya no necesitamos `readNum` en el map de metadata.

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

Y ahí está: tenemos el map de metadata (`[id:patientA, replicate:1, type:normal]`) en la primera posición de la tupla de salida, seguido de la tupla de archivos emparejados, como se pretendía.

Por supuesto, esto solo recogerá y procesará ese par específico de archivos.
Si desea experimentar con el procesamiento de múltiples pares, puede intentar agregar comodines al patrón de entrada y ver qué sucede.
Por ejemplo, intente usar `data/patientA_rep1_*_R{1,2}_001.fastq.gz`

### Conclusión

- [`channel.fromFilePairs()` encuentra y empareja automáticamente archivos relacionados](https://www.nextflow.io/docs/latest/reference/channel.html#fromfilepairs)
- Esto simplifica el manejo de lecturas de extremos emparejados en su pipeline
- Los archivos emparejados pueden agruparse como tuplas `[id, [file1, file2]]`
- La extracción de metadata puede hacerse desde el ID de archivo emparejado en lugar de desde archivos individuales

---

## 6. Usar operaciones de archivos en procesos

Ahora pongamos todo esto junto en un proceso simple para reforzar cómo usar operaciones de archivos dentro de un proceso de Nextflow.

Le proporcionamos un módulo de proceso pre-escrito llamado `ANALYZE_READS` que toma una tupla de metadata y un par de archivos de entrada y los analiza.
Podríamos imaginar que esto está haciendo alineamiento de secuencias, o llamado de variantes o cualquier otro paso que tenga sentido para este tipo de datos.

Comencemos.

### 6.1. Importar el proceso y examinar el código

Para usar este proceso en el workflow, solo necesitamos agregar una declaración include de módulo antes del bloque workflow.

Haga la siguiente edición al workflow:

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

Puede abrir el archivo de módulo para examinar su código:

```groovy title="modules/analyze_reads.nf - ejemplo de proceso" linenums="1"
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

!!! note

    Las directivas `tag` y `publishDir` usan sintaxis de closure (`{ ... }`) en lugar de interpolación de cadenas (`"${...}"`).
    Esto es porque estas directivas referencian variables de entrada (`meta`) que no están disponibles hasta el tiempo de ejecución.
    La sintaxis de closure difiere la evaluación hasta que el proceso realmente se ejecuta.

!!! note

    Estamos llamando a nuestro map de metadata `meta` por convención.
    Para una inmersión más profunda en meta maps, vea la misión secundaria [Metadata and meta maps](./metadata.md).

### 6.2. Llamar al proceso en el workflow

Ahora que el proceso está disponible para el workflow, podemos agregar una llamada al proceso `ANALYZE_READS` para ejecutarlo.

Para ejecutarlo en nuestros datos de ejemplo, necesitaremos hacer dos cosas:

1. Dar un nombre al canal remapeado
2. Agregar una llamada al proceso

#### 6.2.1. Nombrar el canal de entrada remapeado

Anteriormente aplicamos las manipulaciones de mapeo directamente al canal de entrada.
Para alimentar el contenido remapeado al proceso `ANALYZE_READS` (y hacerlo de una manera que sea clara y fácil de leer) queremos crear un nuevo canal llamado `ch_samples`.

Podemos hacer eso usando el operador [`set`](https://www.nextflow.io/docs/latest/reference/operator.html#set).

En el workflow principal, reemplace el operador `.view()` con `.set { ch_samples }`, y agregue una línea probando que podemos referirnos al canal por nombre.

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

        // Temporal: mirar dentro de ch_samples
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

#### 6.2.2. Llamar al proceso en los datos

Ahora llamemos realmente al proceso `ANALYZE_READS` en el canal `ch_samples`.

En el workflow principal, haga los siguientes cambios de código:

=== "Después"

    ```groovy title="main.nf" linenums="23"
        // Ejecutar el análisis
        ANALYZE_READS(ch_samples)
    ```

=== "Antes"

    ```groovy title="main.nf" linenums="23"
        // Temporal: mirar dentro de ch_samples
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

El proceso tomó nuestras entradas y creó un nuevo archivo conteniendo la metadata del paciente, como fue diseñado.
¡Espléndido!

### 6.3. Incluir muchos más pacientes

Por supuesto, esto solo está procesando un único par de archivos para un único paciente, que no es exactamente el tipo de alto rendimiento que espera obtener con Nextflow.
Probablemente querrá procesar muchos más datos a la vez.

Recuerde que `channel.fromPath()` acepta un _glob_ como entrada, lo que significa que puede aceptar cualquier número de archivos que coincidan con el patrón.
Por lo tanto, si queremos incluir todos los pacientes, simplemente podemos modificar la cadena de entrada para incluir más pacientes, como se señaló de pasada anteriormente.

Pretendamos que queremos ser tan codiciosos como sea posible.
Haga las siguientes ediciones al workflow:

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

Ejecute el pipeline de nuevo:

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

¡Éxito! ¡Hemos analizado todos los pacientes de una vez! ¿Verdad?

Quizás no.
Si mira más de cerca, tenemos un problema: ¡tenemos dos réplicas para patientA, pero solo un archivo de salida!
Estamos sobrescribiendo el archivo de salida cada vez.

### 6.4. Hacer los archivos publicados únicos

Dado que tenemos acceso a la metadata del paciente, podemos usarla para hacer los archivos publicados únicos incluyendo metadata diferenciadora, ya sea en la estructura de directorios o en los nombres de archivo mismos.

Haga el siguiente cambio al workflow:

=== "Después"

    ```groovy title="modules/analyze_reads.nf" linenums="6"
        publishDir { "results/${meta.type}/${meta.id}/${meta.replicate}" }, mode: 'copy'
    ```

=== "Antes"

    ```groovy title="modules/analyze_reads.nf" linenums="6"
        publishDir { "results/${meta.id}" }, mode: 'copy'
    ```

Aquí mostramos la opción de usar niveles de directorio adicionales para dar cuenta de tipos de muestra y réplicas, pero podría experimentar haciéndolo al nivel del nombre de archivo también.

Ahora ejecute el pipeline una vez más, pero asegúrese de eliminar el directorio de resultados primero para darse un espacio de trabajo limpio:

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

Verifique el directorio de resultados ahora:

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

Y ahí está, toda nuestra metadata, ordenadamente organizada. ¡Eso es éxito!

Hay mucho más que puede hacer una vez que tiene su metadata cargada en un map como este:

1. Crear directorios de salida organizados basados en atributos del paciente
2. Tomar decisiones en procesos basadas en propiedades del paciente
3. Dividir, unir y recombinar datos basados en valores de metadata

Este patrón de mantener la metadata explícita y adjunta a los datos (en lugar de codificada en nombres de archivo) es una práctica fundamental en Nextflow que permite construir flujos de trabajo de análisis robustos y mantenibles.
Puede aprender más sobre esto en la misión secundaria [Metadata and meta maps](./metadata.md).

### Conclusión

- La directiva `publishDir` puede organizar salidas basadas en valores de metadata
- La metadata en tuplas permite la organización estructurada de resultados
- Este enfoque crea flujos de trabajo mantenibles con proveniencia de datos clara
- Los procesos pueden tomar tuplas de metadata y archivos como entrada
- La directiva `tag` proporciona identificación de procesos en los registros de ejecución
- La estructura del workflow separa la creación de canales de la ejecución de procesos

---

## Resumen

En esta misión secundaria, ha aprendido cómo trabajar con archivos en Nextflow, desde operaciones básicas hasta técnicas más avanzadas para manejar colecciones de archivos.

Aplicar estas técnicas en su propio trabajo le permitirá construir flujos de trabajo más eficientes y mantenibles, especialmente cuando trabaje con grandes números de archivos con convenciones de nomenclatura complejas.

### Patrones clave

1.  **Operaciones Básicas de Archivos:** Creamos objetos Path con `file()` y accedimos a atributos de archivo como nombre, extensión y directorio padre, aprendiendo la diferencia entre cadenas y objetos Path.

    - Crear un objeto Path con `file()`

    ```groovy
    myFile = file('path/to/file.txt')
    ```

    - Obtener atributos de archivo

    ```groovy
    println myFile.name       // file.txt
    println myFile.baseName   // file
    println myFile.extension  // txt
    println myFile.parent     // path/to
    ```

2.  **Usar Archivos Remotos**: Aprendimos cómo cambiar transparentemente entre archivos locales y remotos usando URIs, demostrando la capacidad de Nextflow para manejar archivos de varias fuentes sin cambiar la lógica del workflow.

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

3.  **Cargar archivos usando el channel factory `fromPath()`:** Creamos canales a partir de patrones de archivo con `channel.fromPath()` y vimos sus atributos de archivo, incluyendo tipos de objetos.

    - Crear un canal a partir de un patrón de archivo

    ```groovy
     ch_files = channel.fromPath('data/*.fastq.gz')
    ```

    - Obtener atributos de archivo

    ```groovy
     ch_files.view { myFile ->
        println "File object class: ${myFile.class}"
        println "File name: ${myFile.name}"
        println "Simple name: ${myFile.simpleName}"
        println "Extension: ${myFile.extension}"
        println "Parent directory: ${myFile.parent}"
    }
    ```

4.  **Extraer Metadata de Pacientes de Nombres de Archivos:** Usamos `tokenize()` y `replace()` para extraer y estructurar metadata de nombres de archivos, convirtiéndolos en maps organizados.

    ```groovy
    def name = file.name.tokenize('_')
    def patientId = name[0]
    def replicate = name[1].replace('rep', '')
    def type = name[2]
    def readNum = name[3].replace('R', '')
    ```

5.  **Simplificar con channel.fromFilePairs:** Usamos `channel.fromFilePairs()` para emparejar automáticamente archivos relacionados y extraer metadata de IDs de archivo emparejados.

    ```groovy
    ch_pairs = channel.fromFilePairs('data/*_R{1,2}_001.fastq.gz')
    ```

6.  **Usar Operaciones de Archivos en Procesos:** Integramos operaciones de archivos en procesos de Nextflow con el manejo apropiado de entradas, usando `publishDir` para organizar salidas basadas en metadata.

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

    - Organizar salidas basadas en metadata

    ```groovy
    publishDir { "results/${meta.type}/${meta.id}/${meta.replicate}" }, mode: 'copy'
    ```

### Recursos adicionales

- [Documentación de Nextflow: Working with Files](https://www.nextflow.io/docs/latest/working-with-files.html)
- [channel.fromPath](https://www.nextflow.io/docs/latest/reference/channel.html#frompath)
- [channel.fromFilePairs](https://www.nextflow.io/docs/latest/reference/channel.html#fromfilepairs)

---

## ¿Qué sigue?

Regrese al [menú de Misiones Secundarias](./index.md) o haga clic en el botón en la esquina inferior derecha de la página para pasar al siguiente tema en la lista.
