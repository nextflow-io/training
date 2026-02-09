# Parte 4: Crear un módulo nf-core

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Traducción asistida por IA - [más información y sugerencias](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

En esta cuarta parte del curso de capacitación Hello nf-core, te mostramos cómo crear un módulo nf-core aplicando las convenciones clave que hacen que los módulos sean portables y mantenibles.

El proyecto nf-core proporciona un comando (`nf-core modules create`) que genera plantillas de módulos estructuradas correctamente de forma automática, similar a lo que usamos para el workflow en la Parte 2.
Sin embargo, con fines didácticos, vamos a comenzar haciéndolo manualmente: transformando el módulo local `cowpy` en tu pipeline `core-hello` en un módulo estilo nf-core paso a paso.
Después de eso, te mostraremos cómo usar la creación de módulos basada en plantillas para trabajar de manera más eficiente en el futuro.

??? info "Cómo comenzar desde esta sección"

    Esta sección asume que has completado la [Parte 3: Usar un módulo nf-core](./03_use_module.md) y has integrado el módulo `CAT_CAT` en tu pipeline.

    Si no completaste la Parte 3 o quieres comenzar de nuevo para esta parte, puedes usar la solución `core-hello-part3` como punto de partida.
    Ejecuta estos comandos desde dentro del directorio `hello-nf-core/`:

    ```bash
    cp -r solutions/core-hello-part3 core-hello
    cd core-hello
    ```

    Esto te proporciona un pipeline con el módulo `CAT_CAT` ya integrado.
    Puedes verificar que se ejecuta correctamente ejecutando el siguiente comando:

    ```bash
    nextflow run . --outdir core-hello-results -profile test,docker --validate_params false
    ```

---

## 1. Transformar `cowpy` en un módulo nf-core

En esta sección, aplicaremos las convenciones nf-core al módulo local `cowpy` en tu pipeline `core-hello`, transformándolo en un módulo que sigue los estándares de la comunidad nf-core.

Este es el código actual para el módulo de proceso `cowpy`:

```groovy title="core-hello/modules/local/cowpy.nf" linenums="1"
#!/usr/bin/env nextflow

// Generate ASCII art with cowpy (https://github.com/jeffbuttars/cowpy)
process cowpy {

    publishDir 'results', mode: 'copy'

    container 'community.wave.seqera.io/library/cowpy:1.1.5--3db457ae1977a273'
    conda 'conda-forge::cowpy==1.1.5'

    input:
    path input_file
    val character

    output:
    path "cowpy-${input_file}"

    script:
    """
    cat $input_file | cowpy -c "$character" > cowpy-${input_file}
    """
}
```

Aplicaremos las siguientes convenciones nf-core de forma incremental:

1. **Poner en mayúsculas el nombre del proceso a `COWPY`** para seguir la convención.
2. **Actualizar `COWPY` para usar tuplas de metadata** para propagar metadata de muestra a través del workflow.
3. **Centralizar la configuración de argumentos de herramientas con `ext.args`** para aumentar la versatilidad del módulo mientras se mantiene la interfaz mínima.
4. **Estandarizar el nombramiento de salidas con `ext.prefix`** para promover la consistencia.
5. **Centralizar la configuración de publicación** para promover la consistencia.

Después de cada paso, ejecutaremos el pipeline para verificar que todo funciona como se espera.

!!! warning "Directorio de trabajo"

    Asegúrate de estar en el directorio `core-hello` (la raíz de tu pipeline) para todas las ediciones de archivos y ejecuciones de comandos en esta sección.

    ```bash
    cd core-hello
    ```

### 1.1. Poner en mayúsculas el nombre del proceso

Esta es puramente una convención estilística (no hay justificación técnica) pero como es la norma para los módulos nf-core, cumplamos con ella.

Necesitamos hacer tres conjuntos de cambios:

1. Actualizar el nombre del proceso en el módulo
2. Actualizar la declaración de importación del módulo en el encabezado del workflow
3. Actualizar la llamada al proceso y la declaración emit en el cuerpo del workflow

¡Comencemos!

#### 1.1.1. Actualizar el nombre del proceso en el módulo

Abre el archivo del módulo `cowpy.nf` (bajo `core-hello/modules/local/`) y modifica el nombre del proceso a mayúsculas:

=== "Después"

    ```groovy title="core-hello/modules/local/cowpy.nf" linenums="3" hl_lines="2"
    // Generate ASCII art with cowpy (https://github.com/jeffbuttars/cowpy)
    process COWPY {
    ```

=== "Antes"

    ```groovy title="core-hello/modules/local/cowpy.nf" linenums="3" hl_lines="2"
    // Generate ASCII art with cowpy (https://github.com/jeffbuttars/cowpy)
    process cowpy {
    ```

En este caso, poner en mayúsculas es completamente directo.

Si el nombre del proceso estuviera compuesto de varias palabras, por ejemplo si tuviéramos un proceso llamado MyCowpyTool originalmente en camel case, la convención nf-core sería usar guiones bajos para separarlas, resultando en MY_COWPY_TOOL.

#### 1.1.2. Actualizar la declaración de importación del módulo

Los nombres de procesos distinguen entre mayúsculas y minúsculas, así que ahora que hemos cambiado el nombre del proceso, necesitamos actualizar la declaración de importación del módulo en consecuencia en el encabezado del workflow de `hello.nf`:

=== "Después"

    ```groovy title="core-hello/workflows/hello.nf" linenums="1" hl_lines="10"
    /*
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    */
    include { paramsSummaryMap       } from 'plugin/nf-schema'
    include { softwareVersionsToYAML } from '../subworkflows/nf-core/utils_nfcore_pipeline'
    include { sayHello               } from '../modules/local/sayHello.nf'
    include { convertToUpper         } from '../modules/local/convertToUpper.nf'
    include { COWPY                  } from '../modules/local/cowpy/main.nf'
    include { CAT_CAT                } from '../modules/nf-core/cat/cat/main'
    ```

=== "Antes"

    ```groovy title="core-hello/workflows/hello.nf" linenums="1" hl_lines="10"
    /*
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    */
    include { paramsSummaryMap       } from 'plugin/nf-schema'
    include { softwareVersionsToYAML } from '../subworkflows/nf-core/utils_nfcore_pipeline'
    include { sayHello               } from '../modules/local/sayHello.nf'
    include { convertToUpper         } from '../modules/local/convertToUpper.nf'
    include { cowpy                  } from '../modules/local/cowpy/main.nf'
    include { CAT_CAT                } from '../modules/nf-core/cat/cat/main'
    ```

Podríamos usar un alias en la declaración de importación para evitar tener que actualizar las llamadas al proceso, pero eso de alguna manera anularía el propósito de adoptar la convención de mayúsculas.

#### 1.1.3. Actualizar la llamada al proceso y la declaración emit

Así que ahora actualicemos las dos referencias al proceso en el bloque workflow de `hello.nf`:

=== "Después"

    ```groovy title="core-hello/workflows/hello.nf" linenums="43" hl_lines="2 17"
    // generate ASCII art of the greetings with cowpy
    COWPY(CAT_CAT.out.file_out)

    //
    // Collate and save software versions
    //
    softwareVersionsToYAML(ch_versions)
        .collectFile(
            storeDir: "${params.outdir}/pipeline_info",
            name:  'hello_software_'  + 'versions.yml',
            sort: true,
            newLine: true
        ).set { ch_collated_versions }


    emit:
    cowpy_hellos   = COWPY.out.cowpy_output
    versions       = ch_versions
    ```

=== "Antes"

    ```groovy title="core-hello/workflows/hello.nf" linenums="43" hl_lines="2 17"
    // generate ASCII art of the greetings with cowpy
    cowpy(CAT_CAT.out.file_out)

    //
    // Collate and save software versions
    //
    softwareVersionsToYAML(ch_versions)
        .collectFile(
            storeDir: "${params.outdir}/pipeline_info",
            name:  'hello_software_'  + 'versions.yml',
            sort: true,
            newLine: true
        ).set { ch_collated_versions }


    emit:
    cowpy_hellos   = cowpy.out.cowpy_output
    versions       = ch_versions
    ```

Asegúrate de hacer **ambos** cambios, de lo contrario obtendrás un error cuando ejecutes esto.

#### 1.1.4. Ejecutar el pipeline para probarlo

Ejecutemos el workflow para verificar que todo funciona correctamente después de estos cambios.

```bash
nextflow run . --outdir core-hello-results -profile test,docker --validate_params false
```

??? success "Salida del comando"

    ```console
     N E X T F L O W   ~  version 25.04.3

    Launching `./main.nf` [elegant_plateau] DSL2 - revision: b9e9b3b8de

    Input/output options
      input                     : /workspaces/training/hello-nf-core/core-hello/assets/greetings.csv
      outdir                    : core-hello-results

    Institutional config options
      config_profile_name       : Test profile
      config_profile_description: Minimal test dataset to check pipeline function

    Generic options
      validate_params           : false
      trace_report_suffix       : 2026-01-06_04-51-29

    Core Nextflow options
      runName                   : elegant_plateau
      containerEngine           : docker
      launchDir                 : /workspaces/training/hello-nf-core/core-hello
      workDir                   : /workspaces/training/hello-nf-core/core-hello/work
      projectDir                : /workspaces/training/hello-nf-core/core-hello
      userName                  : root
      profile                   : test,docker
      configFiles               : /workspaces/training/hello-nf-core/core-hello/nextflow.config

    !! Only displaying parameters that differ from the pipeline defaults !!
    ------------------------------------------------------
    executor >  local (8)
    [7b/66ceb5] CORE_HELLO:HELLO:sayHello (3)       | 3 of 3 ✔
    [8e/1bafb9] CORE_HELLO:HELLO:convertToUpper (3) | 3 of 3 ✔
    [bb/203575] CORE_HELLO:HELLO:CAT_CAT (test)     | 1 of 1 ✔
    [39/715489] CORE_HELLO:HELLO:COWPY              | 1 of 1 ✔
    -[core/hello] Pipeline completed successfully-
    ```

¡Muy bien, esto funciona! Ahora pasemos a hacer cambios más sustanciales.

### 1.2. Actualizar `COWPY` para usar tuplas de metadata

En la versión actual del pipeline `core-hello`, estamos extrayendo el archivo de la tupla de salida de `CAT_CAT` para pasarlo a `COWPY`, como se muestra en la mitad superior del diagrama a continuación.

<figure class="excalidraw">
    --8<-- "docs/en/docs/hello_nf-core/img/cowpy-inputs.svg"
</figure>

Sería mejor que `COWPY` aceptara tuplas de metadata directamente, permitiendo que la metadata fluya a través del workflow, como se muestra en la mitad inferior del diagrama.

Para lograrlo, necesitaremos hacer los siguientes cambios:

1. Actualizar las definiciones de entrada y salida
2. Actualizar la llamada al proceso en el workflow
3. Actualizar el bloque emit en el workflow

Una vez que hayamos hecho todo eso, ejecutaremos el pipeline para verificar que todo sigue funcionando como antes.

#### 1.2.1. Actualizar las definiciones de entrada y salida

Regresa al archivo del módulo `cowpy.nf` y modifícalo para aceptar tuplas de metadata como se muestra a continuación.

=== "Después"

    ```groovy title="core-hello/modules/local/cowpy.nf" linenums="11" hl_lines="2 6"
        input:
            tuple val(meta), path(input_file)
            val character

        output:
            tuple val(meta), path("cowpy-${input_file}"), emit: cowpy_output
    ```

=== "Antes"

    ```groovy title="core-hello/modules/local/cowpy.nf" linenums="11" hl_lines="2 6"
        input:
            path input_file
            val character

        output:
            path "cowpy-${input_file}"
    ```

Como puedes ver, cambiamos tanto la **entrada principal** como la **salida** a una tupla que sigue el patrón `tuple val(meta), path(input_file)` introducido en la Parte 3 de esta capacitación.
Para la salida, también aprovechamos esta oportunidad para agregar `emit: cowpy_output` con el fin de darle un nombre descriptivo al canal de salida.

Ahora que hemos cambiado lo que el proceso espera, necesitamos actualizar lo que le proporcionamos en la llamada al proceso.

#### 1.2.2. Actualizar la llamada al proceso en el workflow

La buena noticia es que este cambio simplificará la llamada al proceso.
Ahora que la salida de `CAT_CAT` y la entrada de `COWPY` tienen la misma 'forma', es decir, ambas consisten en una estructura `tuple val(meta), path(input_file)`, podemos simplemente conectarlas directamente en lugar de tener que extraer el archivo explícitamente de la salida del proceso `CAT_CAT`.

Abre el archivo del workflow `hello.nf` (bajo `core-hello/workflows/`) y actualiza la llamada a `COWPY` como se muestra a continuación.

=== "Después"

    ```groovy title="core-hello/workflows/hello.nf" linenums="43" hl_lines="2"
        // generate ASCII art of the greetings with cowpy
        COWPY(CAT_CAT.out.file_out, params.character)
    ```

=== "Antes"

    ```groovy title="core-hello/workflows/hello.nf" linenums="43" hl_lines="1-2 5"
        // extract the file from the tuple since cowpy doesn't use metadata yet
        ch_for_cowpy = CAT_CAT.out.file_out.map{ meta, file -> file }

        // generate ASCII art of the greetings with cowpy
        COWPY(ch_for_cowpy, params.character)
    ```

Ahora llamamos a `COWPY` directamente sobre `CAT_CAT.out.file_out`.

Como resultado, ya no necesitamos construir el canal `ch_for_cowpy`, por lo que esa línea (y su línea de comentario) puede eliminarse por completo.

#### 1.2.3. Actualizar el bloque emit en el workflow

Dado que `COWPY` ahora emite una salida nombrada, `cowpy_output`, podemos actualizar el bloque `emit:` del workflow `hello.nf` para usarla.

=== "Después"

    ```groovy title="core-hello/workflows/hello.nf" linenums="60" hl_lines="2"
        emit:
        cowpy_hellos   = COWPY.out.cowpy_output
        versions       = ch_versions
    ```

=== "Antes"

    ```groovy title="core-hello/workflows/hello.nf" linenums="60" hl_lines="2"
        emit:
        cowpy_hellos   = COWPY.out
        versions       = ch_versions
    ```

Técnicamente esto no es requerido, pero es una buena práctica referirse a salidas nombradas siempre que sea posible.

#### 1.2.4. Ejecutar el pipeline para probarlo

Ejecutemos el workflow para verificar que todo funciona correctamente después de estos cambios.

```bash
nextflow run . --outdir core-hello-results -profile test,docker --validate_params false
```

??? success "Salida del comando"

    ```console
     N E X T F L O W   ~  version 25.04.3

    Launching `./main.nf` [modest_saha] DSL2 - revision: b9e9b3b8de

    Downloading plugin nf-schema@2.5.1
    Input/output options
      input                     : /workspaces/training/hello-nf-core/core-hello/assets/greetings.csv
      outdir                    : core-hello-results

    Institutional config options
      config_profile_name       : Test profile
      config_profile_description: Minimal test dataset to check pipeline function

    Generic options
      validate_params           : false
      trace_report_suffix       : 2025-12-27_06-16-55

    Core Nextflow options
      runName                   : modest_saha
      containerEngine           : docker
      launchDir                 : /workspaces/training/hello-nf-core/core-hello
      workDir                   : /workspaces/training/hello-nf-core/core-hello/work
      projectDir                : /workspaces/training/hello-nf-core/core-hello
      userName                  : root
      profile                   : test,docker
      configFiles               : /workspaces/training/hello-nf-core/core-hello/nextflow.config

    !! Only displaying parameters that differ from the pipeline defaults !!
    ------------------------------------------------------
    executor >  local (8)
    [a8/447993] CORE_HELLO:HELLO:sayHello (3)       [100%] 3 of 3 ✔
    [00/1fc59c] CORE_HELLO:HELLO:convertToUpper (3) [100%] 3 of 3 ✔
    [57/ac800d] CORE_HELLO:HELLO:CAT_CAT (test)     [100%] 1 of 1 ✔
    [b7/092f2b] CORE_HELLO:HELLO:COWPY              [100%] 1 of 1 ✔
    -[core/hello] Pipeline completed successfully-
    ```

El pipeline debería ejecutarse exitosamente, con la metadata ahora fluyendo desde `CAT_CAT` a través de `COWPY`.

Eso completa lo que necesitábamos hacer para que `COWPY` maneje tuplas de metadata.
Ahora, veamos qué más podemos hacer para aprovechar los patrones de módulos nf-core.

### 1.3. Centralizar la configuración de argumentos de herramientas con `ext.args`

En su estado actual, el proceso `COWPY` espera recibir un valor para el parámetro `character`.
Como resultado, tenemos que proporcionar un valor cada vez que llamamos al proceso, incluso si estaríamos contentos con los valores predeterminados establecidos por la herramienta.
Para `COWPY` esto admitidamente no es un gran problema, pero para herramientas con muchos parámetros opcionales, puede volverse bastante engorroso.

El proyecto nf-core recomienda usar una característica de Nextflow llamada [`ext.args`](https://www.nextflow.io/docs/latest/reference/process.html#ext) para gestionar los argumentos de herramientas de manera más conveniente a través de archivos de configuración.

En lugar de declarar entradas de proceso para cada opción de herramienta, escribes el módulo para referenciar `ext.args` en la construcción de su línea de comandos.
Luego es solo cuestión de configurar la variable `ext.args` para contener los argumentos y valores que quieres usar en el archivo `modules.config`, que consolida los detalles de configuración para todos los módulos.
Nextflow agregará esos argumentos con sus valores en la línea de comandos de la herramienta en tiempo de ejecución.

Apliquemos este enfoque al módulo `COWPY`.
Vamos a necesitar hacer los siguientes cambios:

1. Actualizar el módulo `COWPY`
2. Configurar `ext.args` en el archivo `modules.config`
3. Actualizar el workflow `hello.nf`

Una vez que hayamos hecho todo eso, ejecutaremos el pipeline para verificar que todo sigue funcionando como antes.

#### 1.3.1. Actualizar el módulo `COWPY`

Hagámoslo.
Abre el archivo del módulo `cowpy.nf` (bajo `core-hello/modules/local/`) y modifícalo para referenciar `ext.args` como se muestra a continuación.

=== "Después"

    ```groovy title="modules/local/cowpy.nf" linenums="1" hl_lines="18 20"
    #!/usr/bin/env nextflow

    // Generate ASCII art with cowpy (https://github.com/jeffbuttars/cowpy)
    process COWPY {

        publishDir 'results', mode: 'copy'

        container 'community.wave.seqera.io/library/cowpy:1.1.5--3db457ae1977a273'
        conda 'conda-forge::cowpy==1.1.5'

        input:
            tuple val(meta), path(input_file)

        output:
            tuple val(meta), path("cowpy-${input_file}"), emit: cowpy_output

        script:
        def args = task.ext.args ?: ''
        """
        cat $input_file | cowpy $args > cowpy-${input_file}
        """
    }
    ```

=== "Antes"

    ```groovy title="core-hello/modules/local/cowpy.nf" linenums="1" hl_lines="13 20"
    #!/usr/bin/env nextflow

    // Generate ASCII art with cowpy (https://github.com/jeffbuttars/cowpy)
    process COWPY {

        publishDir 'results', mode: 'copy'

        container 'community.wave.seqera.io/library/cowpy:1.1.5--3db457ae1977a273'
        conda 'conda-forge::cowpy==1.1.5'

        input:
            tuple val(meta), path(input_file)
            val character

        output:
            tuple val(meta), path("cowpy-${input_file}"), emit: cowpy_output

        script:
        """
        cat $input_file | cowpy -c "$character" > cowpy-${input_file}
        """
    }
    ```

Puedes ver que hicimos tres cambios.

1. **En el bloque `input:`, eliminamos la entrada `val character`.**
   En adelante, proporcionaremos ese argumento a través de la configuración `ext.args` como se describe más adelante.

2. **En el bloque `script:`, agregamos la línea `def args = task.ext.args ?: ''`.**
   Esa línea usa el operador `?:` para determinar el valor de la variable `args`: el contenido de `task.ext.args` si no está vacío, o una cadena vacía si lo está.
   Ten en cuenta que aunque generalmente nos referimos a `ext.args`, este código debe referenciar `task.ext.args` para extraer la configuración `ext.args` a nivel de módulo.

3. **En la línea de comandos, reemplazamos `-c "$character"` con `$args`.**
   Aquí es donde Nextflow inyectará cualquier argumento de herramienta establecido en `ext.args` en el archivo `modules.config`.

Como resultado, la interfaz del módulo ahora es más simple: solo espera las entradas esenciales de metadata y archivo.

!!! note

    El operador `?:` a menudo se llama 'operador Elvis' porque parece una cara de Elvis Presley de lado, con el carácter `?` simbolizando la onda en su cabello.

#### 1.3.2. Configurar `ext.args` en el archivo `modules.config`

Ahora que hemos sacado la declaración `character` del módulo, tenemos que agregarla a `ext.args` en el archivo de configuración `modules.config`.

Específicamente, vamos a agregar este pequeño fragmento de código al bloque `process {}`:

```groovy title="Código a agregar"
withName: 'COWPY' {
    ext.args = { "-c ${params.character}" }
}
```

La sintaxis `withName:` asigna esta configuración solo al proceso `COWPY`, y `ext.args = { "-c ${params.character}" }` simplemente compone una cadena que incluirá el valor del parámetro `character`.
Ten en cuenta el uso de llaves, que le dicen a Nextflow que evalúe el valor del parámetro en tiempo de ejecución.

¿Tiene sentido? Agreguémoslo.

Abre `conf/modules.config` y agrega el código de configuración dentro del bloque `process {}` como se muestra a continuación.

=== "Después"

    ```groovy title="core-hello/conf/modules.config" linenums="13" hl_lines="8-10"
    process {
        publishDir = [
            path: { "${params.outdir}/${task.process.tokenize(':')[-1].tokenize('_')[0].toLowerCase()}" },
            mode: params.publish_dir_mode,
            saveAs: { filename -> filename.equals('versions.yml') ? null : filename }
        ]

        withName: 'COWPY' {
            ext.args = { "-c ${params.character}" }
        }
    }
    ```

=== "Antes"

    ```groovy title="core-hello/conf/modules.config" linenums="13"
    process {
        publishDir = [
            path: { "${params.outdir}/${task.process.tokenize(':')[-1].tokenize('_')[0].toLowerCase()}" },
            mode: params.publish_dir_mode,
            saveAs: { filename -> filename.equals('versions.yml') ? null : filename }
        ]
    }
    ```

Esperamos que puedas imaginar tener todos los módulos en un pipeline con sus `ext.args` especificados en este archivo, con los siguientes beneficios:

- La **interfaz del módulo se mantiene simple** - Solo acepta las entradas esenciales de metadata y archivos
- El **pipeline aún expone `params.character`** - Los usuarios finales aún pueden configurarlo como antes
- El **módulo ahora es portable** - Puede reutilizarse en otros pipelines sin esperar un nombre de parámetro específico
- La configuración está **centralizada** en `modules.config`, manteniendo limpia la lógica del workflow

Al usar el archivo `modules.config` como el lugar donde todos los pipelines centralizan la configuración por módulo, hacemos que nuestros módulos sean más reutilizables en diferentes pipelines.

#### 1.3.3. Actualizar el workflow `hello.nf`

Dado que el módulo `COWPY` ya no requiere el parámetro `character` como entrada, necesitamos actualizar la llamada del workflow en consecuencia.

Abre el archivo del workflow `hello.nf` (bajo `core-hello/workflows/`) y actualiza la llamada a `COWPY` como se muestra a continuación.

=== "Después"

    ```groovy title="core-hello/workflows/hello.nf" linenums="39" hl_lines="2"
        // generate ASCII art of the greetings with cowpy
        COWPY(CAT_CAT.out.file_out)
    ```

=== "Antes"

    ```groovy title="core-hello/workflows/hello.nf" linenums="39" hl_lines="2"
        // generate ASCII art of the greetings with cowpy
        COWPY(CAT_CAT.out.file_out, params.character)
    ```

El código del workflow ahora es más limpio: no necesitamos pasar `params.character` directamente al proceso.
La interfaz del módulo se mantiene mínima, haciéndola más portable, mientras que el pipeline aún proporciona la opción explícita a través de la configuración.

#### 1.3.4. Ejecutar el pipeline para probarlo

Probemos que el workflow aún funciona como se espera, especificando un carácter diferente para verificar que la configuración `ext.args` está funcionando.

Ejecuta este comando usando `kosh`, una de las opciones más... enigmáticas:

```bash
nextflow run . --outdir core-hello-results -profile test,docker --validate_params false --character kosh
```

??? success "Salida del comando"

    ```console
     N E X T F L O W   ~  version 25.04.3

    Launching `./main.nf` [exotic_planck] DSL2 - revision: b9e9b3b8de

    Input/output options
      input                     : /workspaces/training/hello-nf-core/core-hello/assets/greetings.csv
      outdir                    : core-hello-results

    Institutional config options
      config_profile_name       : Test profile
      config_profile_description: Minimal test dataset to check pipeline function

    Generic options
      validate_params           : false
      trace_report_suffix       : 2025-12-27_06-23-13

    Core Nextflow options
      runName                   : exotic_planck
      containerEngine           : docker
      launchDir                 : /workspaces/training/hello-nf-core/core-hello
      workDir                   : /workspaces/training/hello-nf-core/core-hello/work
      projectDir                : /workspaces/training/hello-nf-core/core-hello
      userName                  : root
      profile                   : test,docker
      configFiles               : /workspaces/training/hello-nf-core/core-hello/nextflow.config

    !! Only displaying parameters that differ from the pipeline defaults !!
    ------------------------------------------------------
    executor >  local (8)
    [13/9e3c0e] CORE_HELLO:HELLO:sayHello (3)       [100%] 3 of 3 ✔
    [e2/5b0ee5] CORE_HELLO:HELLO:convertToUpper (3) [100%] 3 of 3 ✔
    [b6/4fb569] CORE_HELLO:HELLO:CAT_CAT (test)     [100%] 1 of 1 ✔
    [38/eb29ea] CORE_HELLO:HELLO:COWPY              [100%] 1 of 1 ✔
    -[core/hello] Pipeline completed successfully-
    ```

Esto debería ejecutarse exitosamente como antes.

Verifiquemos que la configuración `ext.args` funcionó revisando la salida.
Encuentra la salida en el explorador de archivos o usa el hash de la tarea (la parte `38/eb29ea` en el ejemplo anterior) para ver el archivo de salida:

```bash
cat work/38/eb29ea*/cowpy-test.txt
```

??? success "Salida del comando"

    ```console
    _________
    / HELLO   \
    | HOLà    |
    \ BONJOUR /
    ---------
        \
        \
          \
      ___       _____     ___
    /   \     /    /|   /   \
    |     |   /    / |  |     |
    |     |  /____/  |  |     |
    |     |  |    |  |  |     |
    |     |  | {} | /   |     |
    |     |  |____|/    |     |
    |     |    |==|     |     |
    |      \___________/      |
    |                         |
    |                         |
    ```

¡Deberías ver el arte ASCII mostrado con el carácter `kosh`, confirmando que la configuración `ext.args` funcionó!

??? info "(Opcional) Inspeccionar el archivo de comando"

    Si quieres ver exactamente cómo se aplicó la configuración, puedes inspeccionar el archivo `.command.sh`:

    ```bash
    cat work/38/eb29ea*/.command.sh
    ```

    Verás el comando `cowpy` con el argumento `-c kosh`:

    ```console
    #!/usr/bin/env bash
    ...
    cat test.txt | cowpy -c kosh > cowpy-test.txt
    ```

    Esto muestra que el archivo `.command.sh` se generó correctamente basándose en la configuración `ext.args`.

Tómate un momento para pensar en lo que logramos aquí.
Este enfoque mantiene la interfaz del módulo enfocada en datos esenciales (archivos, metadata y cualquier parámetro obligatorio por muestra), mientras que las opciones que controlan el comportamiento de la herramienta se manejan por separado a través de la configuración.

Esto puede parecer innecesario para una herramienta simple como `cowpy`, pero puede marcar una gran diferencia para herramientas de análisis de datos que tienen muchos argumentos opcionales.

Para resumir los beneficios de este enfoque:

- **Interfaz limpia**: El módulo se enfoca en entradas de datos esenciales (metadata y archivos)
- **Flexibilidad**: Los usuarios pueden especificar argumentos de herramientas a través de configuración, incluyendo valores específicos por muestra
- **Consistencia**: Todos los módulos nf-core siguen este patrón
- **Portabilidad**: Los módulos pueden reutilizarse sin opciones de herramientas codificadas
- **Sin cambios en el workflow**: Agregar o cambiar opciones de herramientas no requiere actualizar el código del workflow

!!! note

    El sistema `ext.args` tiene capacidades adicionales poderosas no cubiertas aquí, incluyendo cambiar valores de argumentos dinámicamente basándose en metadata. Consulta las [especificaciones de módulos nf-core](https://nf-co.re/docs/guidelines/components/modules) para más detalles.

### 1.4. Estandarizar el nombramiento de salidas con `ext.prefix`

Ahora que le hemos dado al proceso `COWPY` acceso al metamap, podemos comenzar a aprovechar otro patrón útil de nf-core: nombrar archivos de salida basándose en metadata.

Aquí vamos a usar una característica de Nextflow llamada `ext.prefix` que nos permitirá estandarizar el nombramiento de archivos de salida en todos los módulos usando `meta.id` (el identificador incluido en el metamap), mientras aún podemos configurar módulos individualmente si se desea.

Esto será similar a lo que hicimos con `ext.args`, con algunas diferencias que detallaremos a medida que avancemos.

Apliquemos este enfoque al módulo `COWPY`.
Vamos a necesitar hacer los siguientes cambios:

1. Actualizar el módulo `COWPY`
2. Configurar `ext.prefix` en el archivo `modules.config`

(No se necesitan cambios en el workflow.)

Una vez que hayamos hecho eso, ejecutaremos el pipeline para verificar que todo sigue funcionando como antes.

#### 1.4.1. Actualizar el módulo `COWPY`

Abre el archivo del módulo `cowpy.nf` (bajo `core-hello/modules/local/`) y modifícalo para referenciar `ext.prefix` como se muestra a continuación.

=== "Después"

    ```groovy title="core-hello/modules/local/cowpy.nf" linenums="1" hl_lines="2 6 8"
        output:
            tuple val(meta), path("${prefix}.txt"), emit: cowpy_output

        script:
        def args = task.ext.args ?: ''
        prefix = task.ext.prefix ?: "${meta.id}"
        """
        cat $input_file | cowpy $args > ${prefix}.txt
        """
    }
    ```

=== "Antes"

    ```groovy title="core-hello/modules/local/cowpy.nf" linenums="1" hl_lines="2 7"
        output:
            tuple val(meta), path("cowpy-${input_file}"), emit: cowpy_output

        script:
        def args = task.ext.args ?: ''
        """
        cat $input_file | cowpy $args > cowpy-${input_file}
        """
    }
    ```

Puedes ver que hicimos tres cambios.

1. **En el bloque `script:`, agregamos la línea `prefix = task.ext.prefix ?: "${meta.id}"`.**
   Esa línea usa el operador `?:` para determinar el valor de la variable `prefix`: el contenido de `task.ext.prefix` si no está vacío, o el identificador del metamap (`meta.id`) si lo está.
   Ten en cuenta que aunque generalmente nos referimos a `ext.prefix`, este código debe referenciar `task.ext.prefix` para extraer la configuración `ext.prefix` a nivel de módulo.

2. **En la línea de comandos, reemplazamos `cowpy-${input_file}` con `${prefix}.txt`.**
   Aquí es donde Nextflow inyectará el valor de `prefix` determinado por la línea anterior.

3. **En el bloque `output:`, reemplazamos `path("cowpy-${input_file}")` con `path("${prefix}.txt")`.**
   Esto simplemente reitera cuál será la ruta del archivo según lo escrito en la línea de comandos.

Como resultado, el nombre del archivo de salida ahora se construye usando un valor predeterminado sensato (el identificador del metamap) combinado con la extensión de formato de archivo apropiada.

#### 1.4.2. Configurar `ext.prefix` en el archivo `modules.config`

En este caso, el valor predeterminado sensato no es suficientemente expresivo para nuestro gusto; queremos usar un patrón de nombramiento personalizado que incluya el nombre de la herramienta, `cowpy-<id>.txt`, como teníamos antes.

Lo haremos configurando `ext.prefix` en `modules.config`, tal como hicimos para el parámetro `character` con `ext.args`, excepto que esta vez el bloque `withName: 'COWPY' {}` ya existe, y solo necesitamos agregar la siguiente línea:

```groovy title="Código a agregar"
ext.prefix = { "cowpy-${meta.id}" }
```

Esto compondrá la cadena que queremos.
Ten en cuenta que una vez más usamos llaves, esta vez para decirle a Nextflow que evalúe el valor de `meta.id` en tiempo de ejecución.

Agreguémoslo.

Abre `conf/modules.config` y agrega el código de configuración dentro del bloque `process {}` como se muestra a continuación.

=== "Después"

    ```groovy title="core-hello/conf/modules.config" linenums="21" hl_lines="3"
        withName: 'COWPY' {
            ext.args = { "-c ${params.character}" }
            ext.prefix = { "cowpy-${meta.id}" }
        }
    ```

=== "Antes"

    ```groovy title="core-hello/conf/modules.config" linenums="21"
        withName: 'COWPY' {
            ext.args = { "-c ${params.character}" }
        }
    ```

En caso de que te lo estés preguntando, el closure `ext.prefix` tiene acceso a la pieza correcta de metadata porque la configuración se evalúa en el contexto de la ejecución del proceso, donde la metadata está disponible.

#### 1.4.3. Ejecutar el pipeline para probarlo

Probemos que el workflow aún funciona como se espera.

```bash
nextflow run . --outdir core-hello-results -profile test,docker --validate_params false
```

??? success "Salida del comando"

    ```console
     N E X T F L O W   ~  version 25.04.3

    Launching `./main.nf` [admiring_turing] DSL2 - revision: b9e9b3b8de

    Input/output options
      input                     : /workspaces/training/hello-nf-core/core-hello/assets/greetings.csv
      outdir                    : core-hello-results

    Institutional config options
      config_profile_name       : Test profile
      config_profile_description: Minimal test dataset to check pipeline function

    Generic options
      validate_params           : false
      trace_report_suffix       : 2025-12-27_06-29-02

    Core Nextflow options
      runName                   : admiring_turing
      containerEngine           : docker
      launchDir                 : /workspaces/training/hello-nf-core/core-hello
      workDir                   : /workspaces/training/hello-nf-core/core-hello/work
      projectDir                : /workspaces/training/hello-nf-core/core-hello
      userName                  : root
      profile                   : test,docker
      configFiles               : /workspaces/training/hello-nf-core/core-hello/nextflow.config

    !! Only displaying parameters that differ from the pipeline defaults !!
    ------------------------------------------------------
    executor >  local (8)
    [b2/e08524] CORE_HELLO:HELLO:sayHello (3)       [100%] 3 of 3 ✔
    [13/88939f] CORE_HELLO:HELLO:convertToUpper (3) [100%] 3 of 3 ✔
    [23/4554e1] CORE_HELLO:HELLO:CAT_CAT (test)     [100%] 1 of 1 ✔
    [a3/c6cbe9] CORE_HELLO:HELLO:COWPY              [100%] 1 of 1 ✔
    -[core/hello] Pipeline completed successfully-
    ```

Echa un vistazo a la salida en el directorio de resultados.
Deberías ver el archivo de salida de cowpy con el mismo nombramiento que antes: `cowpy-test.txt`, basado en el nombre de lote predeterminado.

??? abstract "Contenido del directorio"

    ```console hl_lines="3"
    results
    ├── Bonjour-output.txt
    ├── cowpy-test.txt
    ├── Hello-output.txt
    ├── Holà-output.txt
    ├── UPPER-Bonjour-output.txt
    ├── UPPER-Hello-output.txt
    └── UPPER-Holà-output.txt
    ```

Siéntete libre de cambiar la configuración `ext.prefix` en `conf/modules.config` para satisfacerte de que puedes cambiar el patrón de nombramiento sin tener que hacer ningún cambio en el código del módulo o del workflow.

Alternativamente, también puedes intentar ejecutar esto nuevamente con un parámetro `--batch` diferente especificado en la línea de comandos para satisfacerte de que esa parte aún es personalizable sobre la marcha.

Esto demuestra cómo `ext.prefix` te permite mantener tu convención de nombramiento preferida mientras mantienes flexible la interfaz del módulo.

Para resumir los beneficios de este enfoque:

- **Nombramiento estandarizado**: Los archivos de salida típicamente se nombran usando IDs de muestra de la metadata
- **Configurable**: Los usuarios pueden anular el nombramiento predeterminado si es necesario
- **Consistente**: Todos los módulos nf-core siguen este patrón
- **Predecible**: Fácil saber cómo se llamarán los archivos de salida

¿Bastante bien, verdad?
Bueno, hay un cambio más importante que necesitamos hacer para mejorar nuestro módulo para que se ajuste a las directrices nf-core.

### 1.5. Centralizar la configuración de publicación

Puede que hayas notado que hemos estado publicando salidas en dos directorios diferentes:

- **`results`** — El directorio de salida original que hemos estado usando desde el principio para nuestros módulos locales, establecido individualmente usando directivas `publishDir` por módulo;
- **`core-hello-results`** — El directorio de salida establecido con `--outdir` en la línea de comandos, que ha estado recibiendo los logs de nf-core y los resultados publicados por `CAT_CAT`.

Esto es desordenado y subóptimo; sería mejor tener una ubicación para todo.
Por supuesto, podríamos ir a cada uno de nuestros módulos locales y actualizar la directiva `publishDir` manualmente para usar el directorio `core-hello-results`, pero ¿qué pasa la próxima vez que decidamos cambiar el directorio de salida?

Tener módulos individuales tomando decisiones de publicación claramente no es el camino a seguir, especialmente en un mundo donde el mismo módulo podría usarse en muchos pipelines diferentes, por personas que tienen diferentes necesidades o preferencias.
Queremos poder controlar dónde se publican las salidas a nivel de la configuración del workflow.

"Oye," podrías decir, "`CAT_CAT` está enviando sus salidas al `--outdir`. ¿Tal vez deberíamos copiar su directiva `publishDir`?"

Sí, esa es una gran idea.

Excepto que no tiene una directiva `publishDir`. (Adelante, mira el código del módulo.)

Eso es porque los pipelines nf-core centralizan el control a nivel del workflow configurando `publishDir` en `conf/modules.config` en lugar de en módulos individuales.
Específicamente, la plantilla nf-core declara una directiva `publishDir` predeterminada (con una estructura de directorios predefinida) que se aplica a todos los módulos a menos que se proporcione una directiva de anulación.

¿No suena increíble? ¿Podría ser que para aprovechar esta directiva predeterminada, todo lo que necesitamos hacer es eliminar la directiva `publishDir` actual de nuestros módulos locales?

Probemos eso en `COWPY` para ver qué sucede, luego veremos el código para la configuración predeterminada para entender cómo funciona.

Finalmente, demostraremos cómo anular el comportamiento predeterminado si se desea.

#### 1.5.1. Eliminar la directiva `publishDir` de `COWPY`

Hagamos esto.
Abre el archivo del módulo `cowpy.nf` (bajo `core-hello/modules/local/`) y elimina la directiva `publishDir` como se muestra a continuación.

=== "Después"

    ```groovy title="core-hello/modules/local/cowpy.nf (extracto)" linenums="1"
    #!/usr/bin/env nextflow

    // Generate ASCII art with cowpy (https://github.com/jeffbuttars/cowpy)
    process COWPY {

        container 'community.wave.seqera.io/library/cowpy:1.1.5--3db457ae1977a273'
        conda 'conda-forge::cowpy==1.1.5'
    ```

=== "Antes"

    ```groovy title="core-hello/modules/local/cowpy.nf (extracto)" linenums="1" hl_lines="6"
    #!/usr/bin/env nextflow

    // Generate ASCII art with cowpy (https://github.com/jeffbuttars/cowpy)
    process COWPY {

        publishDir 'results', mode: 'copy'

        container 'community.wave.seqera.io/library/cowpy:1.1.5--3db457ae1977a273'
        conda 'conda-forge::cowpy==1.1.5'

    ```

¡Eso es todo!

#### 1.5.2. Ejecutar el pipeline para probarlo

Echemos un vistazo a qué sucede si ejecutamos el pipeline ahora.

```bash
nextflow run . --outdir core-hello-results -profile test,docker --validate_params false
```

??? success "Salida del comando"

    ```console
     N E X T F L O W   ~  version 25.04.3

    Launching `./main.nf` [silly_caravaggio] DSL2 - revision: b9e9b3b8de

    Input/output options
      input                     : /workspaces/training/hello-nf-core/core-hello/assets/greetings.csv
      outdir                    : core-hello-results

    Institutional config options
      config_profile_name       : Test profile
      config_profile_description: Minimal test dataset to check pipeline function

    Generic options
      validate_params           : false
      trace_report_suffix       : 2025-12-27_06-35-56

    Core Nextflow options
      runName                   : silly_caravaggio
      containerEngine           : docker
      launchDir                 : /workspaces/training/hello-nf-core/core-hello
      workDir                   : /workspaces/training/hello-nf-core/core-hello/work
      projectDir                : /workspaces/training/hello-nf-core/core-hello
      userName                  : root
      profile                   : test,docker
      configFiles               : /workspaces/training/hello-nf-core/core-hello/nextflow.config

    !! Only displaying parameters that differ from the pipeline defaults !!
    ------------------------------------------------------
    executor >  local (8)
    [db/39978e] CORE_HELLO:HELLO:sayHello (3)       [100%] 3 of 3 ✔
    [b5/bf6a8d] CORE_HELLO:HELLO:convertToUpper (3) [100%] 3 of 3 ✔
    [b7/c61842] CORE_HELLO:HELLO:CAT_CAT (test)     [100%] 1 of 1 ✔
    [46/5839d6] CORE_HELLO:HELLO:COWPY              [100%] 1 of 1 ✔
    -[core/hello] Pipeline completed successfully-
    ```

Echa un vistazo a tu directorio de trabajo actual.
Ahora el `core-hello-results` también contiene las salidas del módulo `COWPY`.

??? abstract "Contenido del directorio"

    ```console hl_lines="4-5"
    core-hello-results/
    ├── cat
    │   └── test.txt
    ├── cowpy
    │   └── cowpy-test.txt
    └── pipeline_info
        ├── execution_report_2025-12-27_06-16-55.html
        ├── execution_report_2025-12-27_06-23-13.html
        ├── execution_report_2025-12-27_06-29-02.html
        ├── execution_report_2025-12-27_06-35-56.html
        ├── execution_timeline_2025-12-27_06-16-55.html
        ├── execution_timeline_2025-12-27_06-23-13.html
        ├── execution_timeline_2025-12-27_06-29-02.html
        ├── execution_timeline_2025-12-27_06-35-56.html
        ├── execution_trace_2025-12-27_06-16-55.txt
        ├── execution_trace_2025-12-27_06-23-13.txt
        ├── execution_trace_2025-12-27_06-29-02.txt
        ├── execution_trace_2025-12-27_06-35-56.txt
        ├── hello_software_versions.yml
        ├── params_2025-12-27_06-17-00.json
        ├── params_2025-12-27_06-23-17.json
        ├── params_2025-12-27_06-29-07.json
        ├── params_2025-12-27_06-36-01.json
        ├── pipeline_dag_2025-12-27_06-16-55.html
        ├── pipeline_dag_2025-12-27_06-23-13.html
        ├── pipeline_dag_2025-12-27_06-29-02.html
        └── pipeline_dag_2025-12-27_06-35-56.html
    ```

Puedes ver que Nextflow creó esta jerarquía de directorios basándose en los nombres del workflow y del módulo.

El código responsable vive en el archivo `conf/modules.config`.
Esta es la configuración `publishDir` predeterminada que es parte de la plantilla nf-core y se aplica a todos los procesos:

```groovy
process {
    publishDir = [
        path: { "${params.outdir}/${task.process.tokenize(':')[-1].tokenize('_')[0].toLowerCase()}" },
        mode: params.publish_dir_mode,
        saveAs: { filename -> filename.equals('versions.yml') ? null : filename }
    ]
}
```

Esto puede parecer complicado, así que veamos cada uno de los tres componentes:

- **`path:`** Determina el directorio de salida basándose en el nombre del proceso.
  El nombre completo de un proceso contenido en `task.process` incluye la jerarquía de importaciones de workflow y módulo (como `CORE_HELLO:HELLO:CAT_CAT`).
  Las operaciones `tokenize` eliminan esa jerarquía para obtener solo el nombre del proceso, luego toman la primera parte antes de cualquier guión bajo (si aplica), y lo convierten a minúsculas.
  Esto es lo que determina que los resultados de `CAT_CAT` se publiquen en `${params.outdir}/cat/`.
- **`mode:`** Controla cómo se publican los archivos (copia, enlace simbólico, etc.).
  Esto es configurable a través del parámetro `params.publish_dir_mode`.
- **`saveAs:`** Filtra qué archivos publicar.
  Este ejemplo excluye archivos `versions.yml` devolviendo `null` para ellos, evitando que se publiquen.

Esto proporciona una lógica consistente para organizar salidas.

La salida se ve aún mejor cuando todos los módulos en un pipeline adoptan esta convención, así que siéntete libre de ir a eliminar las directivas `publishDir` de los otros módulos en tu pipeline.
Este valor predeterminado se aplicará incluso a módulos que no modificamos explícitamente para seguir las directrices nf-core.

Dicho esto, puedes decidir que quieres organizar tus entradas de manera diferente, y la buena noticia es que es fácil hacerlo.

#### 1.5.3. Anular el valor predeterminado

Para anular la directiva `publishDir` predeterminada, simplemente puedes agregar tus propias directivas al archivo `conf/modules.config`.

Por ejemplo, podrías anular el valor predeterminado para un solo proceso usando el selector `withName:`, como en este ejemplo donde agregamos una directiva `publishDir` personalizada para el proceso 'COWPY'.

```groovy title="core-hello/conf/modules.config" linenums="13" hl_lines="8-10"
process {
    publishDir = [
        path: { "${params.outdir}/${task.process.tokenize(':')[-1].tokenize('_')[0].toLowerCase()}" },
        mode: params.publish_dir_mode,
        saveAs: { filename -> filename.equals('versions.yml') ? null : filename }
    ]

    withName: 'COWPY' {
        ext.args = { "-c ${params.character}" }
        publishDir = [
            path: 'my_custom_results'
        ]
    }
}
```

En realidad no vamos a hacer ese cambio, pero siéntete libre de jugar con esto y ver qué lógica puedes implementar.

El punto es que este sistema te da lo mejor de ambos mundos: consistencia por defecto y la flexibilidad para personalizar la configuración bajo demanda.

Para resumir, obtienes:

- **Fuente única de verdad**: Toda la configuración de publicación vive en `modules.config`
- **Valor predeterminado útil**: Los procesos funcionan de inmediato sin configuración por módulo
- **Personalización fácil**: Anula el comportamiento de publicación en la configuración, no en el código del módulo
- **Módulos portables**: Los módulos no codifican ubicaciones de salida

Esto completa el conjunto de características de módulos nf-core que absolutamente debes aprender a usar, pero hay otras que puedes leer en las [especificaciones de módulos nf-core](https://nf-co.re/docs/guidelines/components/modules).

### Conclusión

Ahora sabes cómo adaptar módulos locales para seguir las convenciones nf-core:

- Diseña tus módulos para aceptar y propagar tuplas de metadata;
- Usa `ext.args` para mantener las interfaces de módulos mínimas y portables;
- Usa `ext.prefix` para nombramiento de archivos de salida configurable y estandarizado;
- Adopta la directiva `publishDir` centralizada predeterminada para una estructura de directorio de resultados consistente.

### ¿Qué sigue?

Aprende cómo usar las herramientas integradas de nf-core basadas en plantillas para crear módulos de manera fácil.

---

## 2. Crear un módulo con las herramientas nf-core

Ahora que has aprendido los patrones de módulos nf-core aplicándolos manualmente, veamos cómo crearías módulos en la práctica.

### 2.1. Generar un esqueleto de módulo desde una plantilla

Similar a lo que existe para crear pipelines, el proyecto nf-core proporciona herramientas para generar módulos estructurados correctamente basados en una plantilla, con todos estos patrones incorporados desde el inicio.

#### 2.1.1. Ejecutar el comando de creación de módulo

El comando `nf-core modules create` genera una plantilla de módulo que ya sigue todas las convenciones que has aprendido.

Creemos una nueva versión del módulo `COWPY` con una plantilla mínima ejecutando este comando:

```bash
nf-core modules create --empty-template COWPY
```

La bandera `--empty-template` crea una plantilla inicial limpia sin código extra, facilitando ver la estructura esencial.

El comando se ejecuta de forma interactiva, guiándote a través de la configuración.
Busca automáticamente información de herramientas de repositorios de paquetes como Bioconda y bio.tools para pre-poblar metadata.

Se te pedirán varias opciones de configuración:

- **Información del autor**: Tu nombre de usuario de GitHub para atribución
- **Etiqueta de recursos**: Un conjunto predefinido de requisitos computacionales.
  El proyecto nf-core proporciona etiquetas estándar como `process_single` para herramientas ligeras y `process_high` para las exigentes.
  Estas etiquetas ayudan a gestionar la asignación de recursos en diferentes entornos de ejecución.
- **Requisito de metadata**: Si el módulo necesita información específica de muestra a través de un mapa `meta` (generalmente sí para módulos de procesamiento de datos).

La herramienta maneja la complejidad de encontrar información de paquetes y configurar la estructura, permitiéndote enfocarte en implementar la lógica específica de la herramienta.

#### 2.1.2. Examinar el esqueleto del módulo

La herramienta crea una estructura de módulo completa en `modules/local/` (o `modules/nf-core/` si estás en el repositorio nf-core/modules):

??? abstract "Contenido del directorio"

    ```console
    modules/local/cowpy
    ├── environment.yml
    ├── main.nf
    ├── meta.yml
    └── tests
        └── main.nf.test
    ```

Cada archivo tiene un propósito específico:

- **`main.nf`**: Definición del proceso con todos los patrones nf-core incorporados
- **`meta.yml`**: Documentación del módulo describiendo entradas, salidas y la herramienta
- **`environment.yml`**: Especificación del entorno Conda para dependencias
- **`tests/main.nf.test`**: Casos de prueba nf-test para validar que el módulo funciona

!!! tip "Aprende más sobre pruebas"

    El archivo de prueba generado usa nf-test, un framework de pruebas para pipelines y módulos Nextflow. Para aprender cómo escribir y ejecutar estas pruebas, consulta la [misión secundaria nf-test](../side_quests/nf-test.md).

El `main.nf` generado incluye todos los patrones que acabas de aprender, más algunas características adicionales:

```groovy title="modules/local/cowpy/main.nf" hl_lines="11 21 22"
process COWPY {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/YOUR-TOOL-HERE':
        'biocontainers/YOUR-TOOL-HERE' }"

    input:
    tuple val(meta), path(input)        // Patrón 1: Tuplas de metadata ✓

    output:
    tuple val(meta), path("*"), emit: output
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''                // Patrón 2: ext.args ✓
    def prefix = task.ext.prefix ?: "${meta.id}"  // Patrón 3: ext.prefix ✓

    """
    // Add your tool command here

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        COWPY: \$(cowpy --version)
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    echo $args
    touch ${prefix}.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        COWPY: \$(cowpy --version)
    END_VERSIONS
    """
}
```

¡Nota cómo todos los patrones que aplicaste manualmente arriba ya están ahí!

La plantilla también incluye varias convenciones adicionales de nf-core.
Algunas de estas funcionan de inmediato, mientras que otras son marcadores de posición que necesitaremos completar, como se describe a continuación.

**Características que funcionan tal cual:**

- **`tag "$meta.id"`**: Agrega el ID de muestra a los nombres de procesos en los logs para un seguimiento más fácil
- **`label 'process_single'`**: Etiqueta de recursos para configurar requisitos de CPU/memoria
- **Bloque `when:`**: Permite ejecución condicional a través de la configuración `task.ext.when`

Estas características ya son funcionales y hacen que los módulos sean más mantenibles.

**Marcadores de posición que personalizaremos a continuación:**

- **Bloques `input:` y `output:`**: Declaraciones genéricas que actualizaremos para que coincidan con nuestra herramienta
- **Bloque `script:`**: Contiene un comentario donde agregaremos el comando `cowpy`
- **Bloque `stub:`**: Plantilla que actualizaremos para producir las salidas correctas
- **Container y environment**: Marcadores de posición que completaremos con información de paquetes

Las siguientes secciones recorren la finalización de estas personalizaciones.

### 2.2. Configurar el contenedor y el entorno conda

Las directrices nf-core requieren que especifiquemos tanto un contenedor como un entorno Conda como parte del módulo.

#### 2.2.1. Contenedor

Para el contenedor, puedes usar [Seqera Containers](https://seqera.io/containers/) para construir automáticamente un contenedor desde cualquier paquete Conda, incluyendo paquetes conda-forge.
En este caso estamos usando el mismo contenedor precompilado que antes.

El código predeterminado ofrece alternar entre Docker y Singularity, pero vamos a simplificar esa línea y solo especificar el contenedor Docker que obtuvimos de Seqera Containers arriba.

=== "Después"

```groovy title="modules/local/cowpy/main.nf" linenums="3" hl_lines="6"
process COWPY {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "community.wave.seqera.io/library/cowpy:1.1.5--3db457ae1977a273"
```

=== "Antes"

```groovy title="modules/local/cowpy/main.nf" linenums="3" hl_lines="6"
process COWPY {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/YOUR-TOOL-HERE':
        'biocontainers/YOUR-TOOL-HERE' }"
```

#### 2.2.2. Entorno Conda

Para el entorno Conda, el código del módulo especifica `conda "${moduleDir}/environment.yml"` lo que significa que debe configurarse en el archivo `environment.yml`.

La herramienta de creación de módulos nos advirtió que no pudo encontrar el paquete `cowpy` en Bioconda (el canal principal para herramientas bioinformáticas).
Sin embargo, `cowpy` está disponible en conda-forge, así que puedes completar el `environment.yml` así:

=== "Después"

    ```yaml title="modules/local/cowpy/environment.yml"  linenums="1" hl_lines="1 3 5"
    name: COWPY
    channels:
      - conda-forge
    dependencies:
      - cowpy=1.1.5
    ```

=== "Antes"

    ```yaml title="modules/local/cowpy/environment.yml" linenums="1"
    ---
    # yaml-language-server: $schema=https://raw.githubusercontent.com/nf-core/modules/master/modules/environment-schema.json
    channels:
      - conda-forge
      - bioconda
    dependencies:
      # TODO nf-core: List required Conda package(s).
      #               Software MUST be pinned to channel (i.e. "bioconda"), version (i.e. "1.10").
      #               For Conda, the build (i.e. "h9402c20_2") must be EXCLUDED to support installation on different operating systems.
      - "YOUR-TOOL-HERE"
    ```

Para envío a nf-core, tendríamos que seguir los valores predeterminados más de cerca, pero para nuestro propio uso podemos simplificar el código de esta manera.

!!! tip "Paquetes Bioconda vs conda-forge"

    - **Paquetes Bioconda**: Automáticamente obtienen BioContainers construidos, proporcionando contenedores listos para usar
    - **Paquetes conda-forge**: Pueden usar Seqera Containers para construir contenedores bajo demanda desde la receta Conda

    La mayoría de las herramientas bioinformáticas están en Bioconda, pero para herramientas conda-forge, Seqera Containers proporciona una solución fácil para la contenedorización.

### 2.3. Conectar la lógica de `COWPY`

Ahora actualicemos los elementos de código que son específicos de lo que hace el proceso `COWPY`: las entradas y salidas, y el bloque script.

#### 2.3.1. Entradas y salidas

La plantilla generada incluye declaraciones genéricas de entrada y salida que necesitarás personalizar para tu herramienta específica.
Mirando hacia atrás a nuestro módulo `COWPY` manual de la sección 1, podemos usarlo como guía.

Actualiza los bloques de entrada y salida:

=== "Después"

    ```groovy title="modules/local/cowpy/main.nf" linenums="8" hl_lines="2 5"
    input:
    tuple val(meta), path(input_file)

    output:
    tuple val(meta), path("${prefix}.txt"), emit: cowpy_output
    path "versions.yml"           , emit: versions
    ```

=== "Antes"

    ```groovy title="modules/local/cowpy/main.nf" linenums="8" hl_lines="2 5"
    input:
    tuple val(meta), path(input)

    output:
    tuple val(meta), path("*"), emit: output
    path "versions.yml"           , emit: versions
    ```

Esto especifica:

- El nombre del parámetro del archivo de entrada (`input_file` en lugar de `input` genérico)
- El nombre del archivo de salida usando el patrón de prefijo configurable (`${prefix}.txt` en lugar del comodín `*`)
- Un nombre de emit descriptivo (`cowpy_output` en lugar de `output` genérico)

Si estás usando el servidor de lenguaje Nextflow para validar sintaxis, la parte `${prefix}` se marcará como un error en esta etapa porque aún no la hemos agregado al bloque script.
Hagámoslo ahora.

#### 2.3.2. El bloque script

La plantilla proporciona un comentario marcador de posición en el bloque script donde debes agregar el comando real de la herramienta.

Basándonos en el módulo que escribimos manualmente antes, deberíamos hacer las siguientes ediciones:

=== "Después"

    ```groovy title="modules/local/cowpy/main.nf" linenums="15" hl_lines="3 6"
    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"

    """
    cat $input_file | cowpy $args > ${prefix}.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        COWPY: \$(cowpy --version)
    END_VERSIONS
    """
    ```

=== "Antes"

    ```groovy title="modules/local/cowpy/main.nf" linenums="15" hl_lines="6"
    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    // Add your tool command here

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        COWPY: \$(cowpy --version)
    END_VERSIONS
    """
    ```

Cambios clave:

- Cambiar `def prefix` a solo `prefix` (sin `def`) para hacerlo accesible en el bloque output
- Reemplazar el comentario con el comando `cowpy` real que usa tanto `$args` como `${prefix}.txt`

Ten en cuenta que si no hubiéramos hecho ya el trabajo de agregar la configuración `ext.args` y `ext.prefix` para el proceso `COWPY` al archivo `modules.config`, necesitaríamos hacerlo ahora.

#### 2.3.3. Implementar el bloque stub

En el contexto de Nextflow, un bloque [stub](https://www.nextflow.io/docs/latest/process.html#stub) te permite definir un script ligero y ficticio usado para prototipado rápido y pruebas de la lógica de un pipeline sin ejecutar el comando real.

<!-- TODO (futuro) Esto está muy resumido pero realmente debería explicarse o al menos enlazar a una explicación sobre stubs (el doc de referencia tampoco es terriblemente útil). Ahora mismo esto probablemente será mayormente sin sentido para cualquiera que no sepa ya sobre stubs. -->

No te preocupes demasiado si esto parece misterioso; incluimos esto por completitud pero también puedes simplemente eliminar la sección stub si no quieres lidiar con ella, ya que es completamente opcional.

=== "Después"

    ```groovy title="modules/local/cowpy/main.nf" linenums="27" hl_lines="3 6"
    stub:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"

    """
    touch ${prefix}.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        COWPY: \$(cowpy --version)
    END_VERSIONS
    """
    ```

=== "Antes"

    ```groovy title="modules/local/cowpy/main.nf" linenums="27" hl_lines="3 6"
    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    echo $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        COWPY: \$(cowpy --version)
    END_VERSIONS
    """
    ```

Cambios clave:

- Cambiar `def prefix` a solo `prefix` para coincidir con el bloque script
- Eliminar la línea `echo $args` (que era solo código marcador de posición de plantilla)
- El stub crea un archivo `${prefix}.txt` vacío que coincide con lo que produce el bloque script

Esto te permite probar la lógica del workflow y el manejo de archivos sin esperar a que la herramienta real se ejecute.

Una vez que hayas completado la configuración del entorno (sección 2.2), entradas/salidas (sección 2.3.1), bloque script (sección 2.3.2) y bloque stub (sección 2.3.3), ¡el módulo está listo para probar!

### 2.4. Intercambiar el nuevo módulo `COWPY` y ejecutar el pipeline

Todo lo que necesitamos hacer para probar esta nueva versión del módulo `COWPY` es cambiar la declaración de importación en el archivo del workflow `hello.nf` para que apunte al nuevo archivo.

=== "Después"

    ```groovy title="workflows/hello.nf" linenums="1" hl_lines="10"
    /*
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    */
    include { paramsSummaryMap       } from 'plugin/nf-schema'
    include { softwareVersionsToYAML } from '../subworkflows/nf-core/utils_nfcore_pipeline'
    include { sayHello               } from '../modules/local/sayHello.nf'
    include { convertToUpper         } from '../modules/local/convertToUpper.nf'
    include { COWPY                  } from '../modules/local/cowpy/main.nf'
    include { CAT_CAT                } from '../modules/nf-core/cat/cat/main'
    ```

=== "Antes"

    ```groovy title="modules/local/cowpy/main.nf" linenums="1" hl_lines="10"
    /*
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    */
    include { paramsSummaryMap       } from 'plugin/nf-schema'
    include { softwareVersionsToYAML } from '../subworkflows/nf-core/utils_nfcore_pipeline'
    include { sayHello               } from '../modules/local/sayHello.nf'
    include { convertToUpper         } from '../modules/local/convertToUpper.nf'
    include { COWPY                  } from '../modules/local/cowpy.nf'
    include { CAT_CAT                } from '../modules/nf-core/cat/cat/main'
    ```

Ejecutemos el pipeline para probarlo.

```bash
nextflow run . --outdir core-hello-results -profile test,docker --validate_params false
```

??? success "Salida del comando"

    ```console hl_lines="33"
      N E X T F L O W   ~  version 25.04.3

    Launching `./main.nf` [prickly_neumann] DSL2 - revision: b9e9b3b8de

    Input/output options
      input                     : /workspaces/training/hello-nf-core/core-hello/assets/greetings.csv
      outdir                    : core-hello-results

    Institutional config options
      config_profile_name       : Test profile
      config_profile_description: Minimal test dataset to check pipeline function

    Generic options
      validate_params           : false
      trace_report_suffix       : 2025-12-27_08-23-51

    Core Nextflow options
      runName                   : prickly_neumann
      containerEngine           : docker
      launchDir                 : /workspaces/training/hello-nf-core/core-hello
      workDir                   : /workspaces/training/hello-nf-core/core-hello/work
      projectDir                : /workspaces/training/hello-nf-core/core-hello
      userName                  : root
      profile                   : test,docker
      configFiles               : /workspaces/training/hello-nf-core/core-hello/nextflow.config

    !! Only displaying parameters that differ from the pipeline defaults !!
    ------------------------------------------------------
    executor >  local (8)
    [e9/008ede] CORE_HELLO:HELLO:sayHello (3)       [100%] 3 of 3 ✔
    [f0/d70cfe] CORE_HELLO:HELLO:convertToUpper (3) [100%] 3 of 3 ✔
    [be/0ecc58] CORE_HELLO:HELLO:CAT_CAT (test)     [100%] 1 of 1 ✔
    [11/8e082f] CORE_HELLO:HELLO:COWPY (test)       [100%] 1 of 1 ✔
    -[core/hello] Pipeline completed successfully-
    ```

Esto produce los mismos resultados que antes.

### Conclusión

Ahora sabes cómo usar las herramientas integradas de nf-core para crear módulos eficientemente usando plantillas en lugar de escribir todo desde cero.

### ¿Qué sigue?

Aprende cuáles son los beneficios de contribuir módulos a nf-core y cuáles son los principales pasos y requisitos involucrados.

---

## 3. Contribuir módulos de vuelta a nf-core

El repositorio [nf-core/modules](https://github.com/nf-core/modules) da la bienvenida a contribuciones de módulos bien probados y estandarizados.

### 3.1. ¿Por qué contribuir?

Contribuir tus módulos a nf-core:

- Hace que tus herramientas estén disponibles para toda la comunidad nf-core a través del catálogo de módulos en [nf-co.re/modules](https://nf-co.re/modules)
- Asegura mantenimiento y mejoras continuas de la comunidad
- Proporciona aseguramiento de calidad a través de revisión de código y pruebas automatizadas
- Da visibilidad y reconocimiento a tu trabajo

### 3.2. Lista de verificación del contribuidor

Para contribuir un módulo a nf-core, necesitarás pasar por los siguientes pasos:

1. Verificar si ya existe en [nf-co.re/modules](https://nf-co.re/modules)
2. Hacer fork del repositorio [nf-core/modules](https://github.com/nf-core/modules)
3. Usar `nf-core modules create` para generar la plantilla
4. Completar la lógica del módulo y las pruebas
5. Probar con `nf-core modules test tool/subtool`
6. Hacer lint con `nf-core modules lint tool/subtool`
7. Enviar un pull request

Para instrucciones detalladas, consulta el [tutorial de componentes nf-core](https://nf-co.re/docs/tutorials/nf-core_components/components).

### 3.3. Recursos

- **Tutorial de componentes**: [Guía completa para crear y contribuir módulos](https://nf-co.re/docs/tutorials/nf-core_components/components)
- **Especificaciones de módulos**: [Requisitos técnicos y directrices](https://nf-co.re/docs/guidelines/components/modules)
- **Soporte de la comunidad**: [Slack de nf-core](https://nf-co.re/join) - Únete al canal `#modules`

### Conclusión

¡Ahora sabes cómo crear módulos nf-core! Aprendiste los cuatro patrones clave que hacen que los módulos sean portables y mantenibles:

- **Tuplas de metadata** propagan metadata a través del workflow
- **`ext.args`** simplifica las interfaces de módulos manejando argumentos opcionales a través de configuración
- **`ext.prefix`** estandariza el nombramiento de archivos de salida
- **Publicación centralizada** a través de `publishDir` configurado en `modules.config` en lugar de codificado en módulos

Al transformar `COWPY` paso a paso, desarrollaste una comprensión profunda de estos patrones, haciéndote capaz de trabajar con, depurar y crear módulos nf-core.
En la práctica, usarás `nf-core modules create` para generar módulos estructurados correctamente con estos patrones incorporados desde el inicio.

Finalmente, aprendiste cómo contribuir módulos a la comunidad nf-core, haciendo herramientas disponibles para investigadores en todo el mundo mientras te beneficias del mantenimiento continuo de la comunidad.

### ¿Qué sigue?

Cuando estés listo, continúa a la [Parte 5: Validación de entrada](./05_input_validation.md) para aprender cómo agregar validación de entrada basada en esquemas a tu pipeline.
