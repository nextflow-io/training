# Parte 4: Crear un módulo nf-core

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Traducción asistida por IA - [más información y sugerencias](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

En esta cuarta parte del curso de entrenamiento Hello nf-core, le mostramos cómo crear un módulo nf-core aplicando las convenciones clave que hacen que los módulos sean portables y mantenibles.

El proyecto nf-core proporciona un comando (`nf-core modules create`) que genera plantillas de módulos estructuradas correctamente de forma automática, similar a lo que usamos para el flujo de trabajo en la Parte 2.
Sin embargo, con fines didácticos, vamos a comenzar haciéndolo manualmente: transformando el módulo local `cowpy` en su pipeline `core-hello` en un módulo de estilo nf-core paso a paso.
Después de eso, le mostraremos cómo usar la creación de módulos basada en plantillas para trabajar de manera más eficiente en el futuro.

??? info "Cómo comenzar desde esta sección"

    Esta sección asume que ha completado la [Parte 3: Usar un módulo nf-core](./03_use_module.md) e ha integrado el módulo `CAT_CAT` en su pipeline.

    Si no completó la Parte 3 o desea comenzar de nuevo para esta parte, puede usar la solución `core-hello-part3` como punto de partida.
    Ejecute estos comandos desde dentro del directorio `hello-nf-core/`:

    ```bash
    cp -r solutions/core-hello-part3 core-hello
    cd core-hello
    ```

    Esto le proporciona un pipeline con el módulo `CAT_CAT` ya integrado.
    Puede probar que se ejecuta correctamente ejecutando el siguiente comando:

    ```bash
    nextflow run . --outdir core-hello-results -profile test,docker --validate_params false
    ```

---

## 1. Transformar `cowpy` en un módulo nf-core

En esta sección, aplicaremos las convenciones de nf-core al módulo local `cowpy` en su pipeline `core-hello`, transformándolo en un módulo que sigue los estándares de la comunidad nf-core.

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

Aplicaremos las siguientes convenciones de nf-core de forma incremental:

1. **Cambiar el nombre del proceso a mayúsculas `COWPY`** para seguir la convención.
2. **Actualizar `COWPY` para usar tuplas de metadatos** para propagar los metadatos de muestra a través del flujo de trabajo.
3. **Centralizar la configuración de argumentos de la herramienta con `ext.args`** para aumentar la versatilidad del módulo mientras se mantiene la interfaz mínima.
4. **Estandarizar el nombre de salida con `ext.prefix`** para promover la consistencia.
5. **Centralizar la configuración de publicación** para promover la consistencia.

Después de cada paso, ejecutaremos el pipeline para probar que todo funciona como se espera.

!!! warning "Directorio de trabajo"

    Asegúrese de estar en el directorio `core-hello` (la raíz de su pipeline) para todas las ediciones de archivos y ejecuciones de comandos en esta sección.

    ```bash
    cd core-hello
    ```

### 1.1. Cambiar el nombre del proceso a mayúsculas

Esta es puramente una convención estilística (no hay justificación técnica), pero dado que es la norma para los módulos nf-core, vamos a cumplirla.

Necesitamos hacer tres conjuntos de cambios:

1. Actualizar el nombre del proceso en el módulo
2. Actualizar la declaración de importación del módulo en el encabezado del flujo de trabajo
3. Actualizar la llamada al proceso y la declaración emit en el cuerpo del flujo de trabajo

¡Comencemos!

#### 1.1.1. Actualizar el nombre del proceso en el módulo

Abra el archivo del módulo `cowpy.nf` (en `core-hello/modules/local/`) y modifique el nombre del proceso a mayúsculas:

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

En este caso, el cambio a mayúsculas es completamente directo.

Si el nombre del proceso estuviera compuesto por varias palabras, por ejemplo si tuviéramos un proceso llamado MyCowpyTool originalmente en camelCase, la convención de nf-core sería usar guiones bajos para separarlas, resultando en MY_COWPY_TOOL.

#### 1.1.2. Actualizar la declaración de importación del módulo

Los nombres de procesos distinguen entre mayúsculas y minúsculas, así que ahora que hemos cambiado el nombre del proceso, necesitamos actualizar la declaración de importación del módulo en consecuencia en el encabezado del flujo de trabajo de `hello.nf`:

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
    // generar arte ASCII de los saludos con cowpy
    COWPY(CAT_CAT.out.file_out)

    //
    // Recopilar y guardar versiones de software
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
    // generar arte ASCII de los saludos con cowpy
    cowpy(CAT_CAT.out.file_out)

    //
    // Recopilar y guardar versiones de software
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

Asegúrese de hacer **ambos** cambios, de lo contrario obtendrá un error cuando ejecute esto.

#### 1.1.4. Ejecutar el pipeline para probarlo

Ejecutemos el flujo de trabajo para probar que todo está funcionando correctamente después de estos cambios.

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

### 1.2. Actualizar `COWPY` para usar tuplas de metadatos

En la versión actual del pipeline `core-hello`, estamos extrayendo el archivo de la tupla de salida de `CAT_CAT` para pasarlo a `COWPY`, como se muestra en la mitad superior del diagrama a continuación.

<figure class="excalidraw">
    --8<-- "docs/en/docs/hello_nf-core/img/cowpy-inputs.svg"
</figure>

Sería mejor que `COWPY` aceptara tuplas de metadatos directamente, permitiendo que los metadatos fluyan a través del flujo de trabajo, como se muestra en la mitad inferior del diagrama.

Para lograrlo, necesitaremos hacer los siguientes cambios:

1. Actualizar las definiciones de entrada y salida
2. Actualizar la llamada al proceso en el flujo de trabajo
3. Actualizar el bloque emit en el flujo de trabajo

Una vez que hayamos hecho todo eso, ejecutaremos el pipeline para probar que todo sigue funcionando como antes.

#### 1.2.1. Actualizar las definiciones de entrada y salida

Regrese al archivo del módulo `cowpy.nf` y modifíquelo para aceptar tuplas de metadatos como se muestra a continuación.

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

Como puede ver, cambiamos tanto la **entrada principal** como la **salida** a una tupla que sigue el patrón `tuple val(meta), path(input_file)` introducido en la Parte 3 de este entrenamiento.
Para la salida, también aprovechamos esta oportunidad para agregar `emit: cowpy_output` con el fin de dar al canal de salida un nombre descriptivo.

Ahora que hemos cambiado lo que el proceso espera, necesitamos actualizar lo que le proporcionamos en la llamada al proceso.

#### 1.2.2. Actualizar la llamada al proceso en el flujo de trabajo

La buena noticia es que este cambio simplificará la llamada al proceso.
Ahora que la salida de `CAT_CAT` y la entrada de `COWPY` tienen la misma 'forma', es decir, ambas consisten en una estructura `tuple val(meta), path(input_file)`, simplemente podemos conectarlas directamente en lugar de tener que extraer el archivo explícitamente de la salida del proceso `CAT_CAT`.

Abra el archivo de flujo de trabajo `hello.nf` (en `core-hello/workflows/`) y actualice la llamada a `COWPY` como se muestra a continuación.

=== "Después"

    ```groovy title="core-hello/workflows/hello.nf" linenums="43" hl_lines="2"
        // generar arte ASCII de los saludos con cowpy
        COWPY(CAT_CAT.out.file_out, params.character)
    ```

=== "Antes"

    ```groovy title="core-hello/workflows/hello.nf" linenums="43" hl_lines="1-2 5"
        // extract the file from the tuple since cowpy doesn't use metadata yet
        ch_for_cowpy = CAT_CAT.out.file_out.map{ meta, file -> file }

        // generar arte ASCII de los saludos con cowpy
        COWPY(ch_for_cowpy, params.character)
    ```

Ahora llamamos a `COWPY` directamente en `CAT_CAT.out.file_out`.

Como resultado, ya no necesitamos construir el canal `ch_for_cowpy`, por lo que esa línea (y su línea de comentario) se pueden eliminar por completo.

#### 1.2.3. Actualizar el bloque emit en el flujo de trabajo

Dado que `COWPY` ahora emite una salida nombrada, `cowpy_output`, podemos actualizar el bloque `emit:` del flujo de trabajo `hello.nf` para usar eso.

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

Ejecutemos el flujo de trabajo para probar que todo está funcionando correctamente después de estos cambios.

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

El pipeline debería ejecutarse exitosamente, con los metadatos ahora fluyendo desde `CAT_CAT` a través de `COWPY`.

Eso completa lo que necesitábamos hacer para que `COWPY` maneje tuplas de metadatos.
Ahora, veamos qué más podemos hacer para aprovechar los patrones de módulos de nf-core.

### 1.3. Centralizar la configuración de argumentos de herramientas con `ext.args`

En su estado actual, el proceso `COWPY` espera recibir un valor para el parámetro `character`.
Como resultado, tenemos que proporcionar un valor cada vez que llamamos al proceso, incluso si estaríamos contentos con los valores predeterminados establecidos por la herramienta.
Para `COWPY` esto admitidamente no es un gran problema, pero para herramientas con muchos parámetros opcionales, puede volverse bastante engorroso.

El proyecto nf-core recomienda usar una característica de Nextflow llamada [`ext.args`](https://www.nextflow.io/docs/latest/reference/process.html#ext) para gestionar los argumentos de herramientas de manera más conveniente a través de archivos de configuración.

En lugar de declarar entradas de proceso para cada opción de herramienta, usted escribe el módulo para referenciar `ext.args` en la construcción de su línea de comandos.
Luego es solo cuestión de configurar la variable `ext.args` para contener los argumentos y valores que desea usar en el archivo `modules.config`, que consolida los detalles de configuración para todos los módulos.
Nextflow agregará esos argumentos con sus valores a la línea de comandos de la herramienta en tiempo de ejecución.

Apliquemos este enfoque al módulo `COWPY`.
Vamos a necesitar hacer los siguientes cambios:

1. Actualizar el módulo `COWPY`
2. Configurar `ext.args` en el archivo `modules.config`
3. Actualizar el flujo de trabajo `hello.nf`

Una vez que hayamos hecho todo eso, ejecutaremos el pipeline para probar que todo sigue funcionando como antes.

#### 1.3.1. Actualizar el módulo `COWPY`

Hagámoslo.
Abra el archivo del módulo `cowpy.nf` (en `core-hello/modules/local/`) y modifíquelo para referenciar `ext.args` como se muestra a continuación.

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

Puede ver que hicimos tres cambios.

1. **En el bloque `input:`, eliminamos la entrada `val character`.**
   De ahora en adelante, proporcionaremos ese argumento a través de la configuración `ext.args` como se describe más adelante.

2. **En el bloque `script:`, agregamos la línea `def args = task.ext.args ?: ''`.**
   Esa línea usa el operador `?:` para determinar el valor de la variable `args`: el contenido de `task.ext.args` si no está vacío, o una cadena vacía si lo está.
   Tenga en cuenta que aunque generalmente nos referimos a `ext.args`, este código debe referenciar `task.ext.args` para extraer la configuración de `ext.args` a nivel de módulo.

3. **En la línea de comandos, reemplazamos `-c "$character"` con `$args`.**
   Aquí es donde Nextflow inyectará cualquier argumento de herramienta configurado en `ext.args` en el archivo `modules.config`.

Como resultado, la interfaz del módulo ahora es más simple: solo espera las entradas esenciales de metadatos y archivos.

!!! note

    El operador `?:` a menudo se llama 'operador Elvis' porque parece una cara de Elvis Presley de lado, con el carácter `?` simbolizando la onda en su cabello.

#### 1.3.2. Configurar `ext.args` en el archivo `modules.config`

Ahora que hemos sacado la declaración de `character` del módulo, tenemos que agregarla a `ext.args` en el archivo de configuración `modules.config`.

Específicamente, vamos a agregar este pequeño fragmento de código al bloque `process {}`:

```groovy title="Código a agregar"
withName: 'COWPY' {
    ext.args = { "-c ${params.character}" }
}
```

La sintaxis `withName:` asigna esta configuración solo al proceso `COWPY`, y `ext.args = { "-c ${params.character}" }` simplemente compone una cadena que incluirá el valor del parámetro `character`.
Tenga en cuenta el uso de llaves, que le dicen a Nextflow que evalúe el valor del parámetro en tiempo de ejecución.

¿Tiene sentido? Agreguémoslo.

Abra `conf/modules.config` y agregue el código de configuración dentro del bloque `process {}` como se muestra a continuación.

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

Esperamos que pueda imaginar tener todos los módulos en un pipeline con sus `ext.args` especificados en este archivo, con los siguientes beneficios:

- La **interfaz del módulo se mantiene simple** - Solo acepta las entradas esenciales de metadatos y archivos
- El **pipeline todavía expone `params.character`** - Los usuarios finales aún pueden configurarlo como antes
- El **módulo ahora es portable** - Puede reutilizarse en otros pipelines sin esperar un nombre de parámetro específico
- La configuración está **centralizada** en `modules.config`, manteniendo limpia la lógica del flujo de trabajo

Al usar el archivo `modules.config` como el lugar donde todos los pipelines centralizan la configuración por módulo, hacemos que nuestros módulos sean más reutilizables en diferentes pipelines.

#### 1.3.3. Actualizar el flujo de trabajo `hello.nf`

Dado que el módulo `COWPY` ya no requiere el parámetro `character` como entrada, necesitamos actualizar la llamada del flujo de trabajo en consecuencia.

Abra el archivo de flujo de trabajo `hello.nf` (en `core-hello/workflows/`) y actualice la llamada a `COWPY` como se muestra a continuación.

=== "Después"

    ```groovy title="core-hello/workflows/hello.nf" linenums="39" hl_lines="2"
        // generar arte ASCII de los saludos con cowpy
        COWPY(CAT_CAT.out.file_out)
    ```

=== "Antes"

    ```groovy title="core-hello/workflows/hello.nf" linenums="39" hl_lines="2"
        // generar arte ASCII de los saludos con cowpy
        COWPY(CAT_CAT.out.file_out, params.character)
    ```

El código del flujo de trabajo ahora es más limpio: no necesitamos pasar `params.character` directamente al proceso.
La interfaz del módulo se mantiene mínima, haciéndola más portable, mientras que el pipeline todavía proporciona la opción explícita a través de la configuración.

#### 1.3.4. Ejecutar el pipeline para probarlo

Probemos que el flujo de trabajo todavía funciona como se espera, especificando un personaje diferente para verificar que la configuración de `ext.args` está funcionando.

Ejecute este comando usando `kosh`, una de las opciones más... enigmáticas:

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

Esto debería ejecutarse exitosamente como anteriormente.

Verifiquemos que la configuración de `ext.args` funcionó verificando la salida.
Encuentre la salida en el navegador de archivos o use el hash de la tarea (la parte `38/eb29ea` en el ejemplo anterior) para ver el archivo de salida:

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

¡Debería ver el arte ASCII mostrado con el personaje `kosh`, confirmando que la configuración de `ext.args` funcionó!

??? info "(Opcional) Inspeccionar el archivo de comando"

    Si desea ver exactamente cómo se aplicó la configuración, puede inspeccionar el archivo `.command.sh`:

    ```bash
    cat work/38/eb29ea*/.command.sh
    ```

    Verá el comando `cowpy` con el argumento `-c kosh`:

    ```console
    #!/usr/bin/env bash
    ...
    cat test.txt | cowpy -c kosh > cowpy-test.txt
    ```

    Esto muestra que el archivo `.command.sh` se generó correctamente basándose en la configuración de `ext.args`.

Tómese un momento para pensar en lo que logramos aquí.
Este enfoque mantiene la interfaz del módulo enfocada en datos esenciales (archivos, metadatos y cualquier parámetro obligatorio por muestra), mientras que las opciones que controlan el comportamiento de la herramienta se manejan por separado a través de la configuración.

Esto puede parecer innecesario para una herramienta simple como `cowpy`, pero puede marcar una gran diferencia para herramientas de análisis de datos que tienen muchos argumentos opcionales.

Para resumir los beneficios de este enfoque:

- **Interfaz limpia**: El módulo se centra en entradas de datos esenciales (metadatos y archivos)
- **Flexibilidad**: Los usuarios pueden especificar argumentos de herramientas a través de la configuración, incluyendo valores específicos por muestra
- **Consistencia**: Todos los módulos nf-core siguen este patrón
- **Portabilidad**: Los módulos pueden reutilizarse sin opciones de herramientas codificadas
- **Sin cambios en el flujo de trabajo**: Agregar o cambiar opciones de herramientas no requiere actualizar el código del flujo de trabajo

!!! note

    El sistema `ext.args` tiene capacidades adicionales poderosas no cubiertas aquí, incluyendo el cambio de valores de argumentos dinámicamente basado en metadatos. Consulte las [especificaciones de módulos nf-core](https://nf-co.re/docs/guidelines/components/modules) para más detalles.

### 1.4. Estandarizar el nombre de salida con `ext.prefix`

Ahora que le hemos dado al proceso `COWPY` acceso al metamap, podemos comenzar a aprovechar otro patrón útil de nf-core: nombrar archivos de salida basándose en metadatos.

Aquí vamos a usar una característica de Nextflow llamada `ext.prefix` que nos permitirá estandarizar el nombre de archivos de salida en todos los módulos usando `meta.id` (el identificador incluido en el metamap), mientras seguimos pudiendo configurar módulos individualmente si se desea.

Esto será similar a lo que hicimos con `ext.args`, con algunas diferencias que detallaremos a medida que avancemos.

Apliquemos este enfoque al módulo `COWPY`.
Vamos a necesitar hacer los siguientes cambios:

1. Actualizar el módulo `COWPY`
2. Configurar `ext.prefix` en el archivo `modules.config`

(No se necesitan cambios en el flujo de trabajo.)

Una vez que hayamos hecho eso, ejecutaremos el pipeline para probar que todo sigue funcionando como antes.

#### 1.4.1. Actualizar el módulo `COWPY`

Abra el archivo del módulo `cowpy.nf` (en `core-hello/modules/local/`) y modifíquelo para referenciar `ext.prefix` como se muestra a continuación.

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

Puede ver que hicimos tres cambios.

1. **En el bloque `script:`, agregamos la línea `prefix = task.ext.prefix ?: "${meta.id}"`.**
   Esa línea usa el operador `?:` para determinar el valor de la variable `prefix`: el contenido de `task.ext.prefix` si no está vacío, o el identificador del metamap (`meta.id`) si lo está.
   Tenga en cuenta que aunque generalmente nos referimos a `ext.prefix`, este código debe referenciar `task.ext.prefix` para extraer la configuración de `ext.prefix` a nivel de módulo.

2. **En la línea de comandos, reemplazamos `cowpy-${input_file}` con `${prefix}.txt`.**
   Aquí es donde Nextflow inyectará el valor de `prefix` determinado por la línea anterior.

3. **En el bloque `output:`, reemplazamos `path("cowpy-${input_file}")` con `path("${prefix}.txt")`.**
   Esto simplemente reitera cuál será la ruta del archivo según lo escrito en la línea de comandos.

Como resultado, el nombre del archivo de salida ahora se construye usando un valor predeterminado sensato (el identificador del metamap) combinado con la extensión de formato de archivo apropiada.

#### 1.4.2. Configurar `ext.prefix` en el archivo `modules.config`

En este caso, el valor predeterminado sensato no es suficientemente expresivo para nuestro gusto; queremos usar un patrón de nomenclatura personalizado que incluya el nombre de la herramienta, `cowpy-<id>.txt`, como teníamos antes.

Lo haremos configurando `ext.prefix` en `modules.config`, tal como lo hicimos para el parámetro `character` con `ext.args`, excepto que esta vez el bloque `withName: 'COWPY' {}` ya existe, y solo necesitamos agregar la siguiente línea:

```groovy title="Código a agregar"
ext.prefix = { "cowpy-${meta.id}" }
```

Esto compondrá la cadena que queremos.
Tenga en cuenta que una vez más usamos llaves, esta vez para decirle a Nextflow que evalúe el valor de `meta.id` en tiempo de ejecución.

Agreguémoslo.

Abra `conf/modules.config` y agregue el código de configuración dentro del bloque `process {}` como se muestra a continuación.

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

En caso de que se esté preguntando, el closure `ext.prefix` tiene acceso a la pieza correcta de metadatos porque la configuración se evalúa en el contexto de la ejecución del proceso, donde los metadatos están disponibles.

#### 1.4.3. Ejecutar el pipeline para probarlo

Probemos que el flujo de trabajo todavía funciona como se espera.

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

Eche un vistazo a la salida en el directorio de resultados.
Debería ver el archivo de salida de cowpy con el mismo nombre que antes: `cowpy-test.txt`, basado en el nombre de lote predeterminado.

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

Siéntase libre de cambiar la configuración de `ext.prefix` en `conf/modules.config` para satisfacerse de que puede cambiar el patrón de nomenclatura sin tener que hacer ningún cambio en el código del módulo o del flujo de trabajo.

Alternativamente, también puede intentar ejecutar esto nuevamente con un parámetro `--batch` diferente especificado en la línea de comandos para satisfacerse de que esa parte todavía es personalizable sobre la marcha.

Esto demuestra cómo `ext.prefix` le permite mantener su convención de nomenclatura preferida mientras mantiene flexible la interfaz del módulo.

Para resumir los beneficios de este enfoque:

- **Nomenclatura estandarizada**: Los archivos de salida se nombran típicamente usando IDs de muestra de metadatos
- **Configurable**: Los usuarios pueden sobrescribir la nomenclatura predeterminada si es necesario
- **Consistente**: Todos los módulos nf-core siguen este patrón
- **Predecible**: Es fácil saber cómo se llamarán los archivos de salida

¿Bastante bien, verdad?
Bueno, hay un cambio más importante que necesitamos hacer para mejorar nuestro módulo para que se ajuste a las directrices de nf-core.

### 1.5. Centralizar la configuración de publicación

Puede haber notado que hemos estado publicando salidas en dos directorios diferentes:

- **`results`** — El directorio de salida original que hemos estado usando desde el principio para nuestros módulos locales, establecido individualmente usando directivas `publishDir` por módulo;
- **`core-hello-results`** — El directorio de salida establecido con `--outdir` en la línea de comandos, que ha estado recibiendo los registros de nf-core y los resultados publicados por `CAT_CAT`.

Esto es desordenado y subóptimo; sería mejor tener una ubicación para todo.
Por supuesto, podríamos ir a cada uno de nuestros módulos locales y actualizar la directiva `publishDir` manualmente para usar el directorio `core-hello-results`, pero ¿qué pasa la próxima vez que decidamos cambiar el directorio de salida?

Tener módulos individuales tomando decisiones de publicación claramente no es el camino a seguir, especialmente en un mundo donde el mismo módulo podría usarse en muchos pipelines diferentes, por personas que tienen diferentes necesidades o preferencias.
Queremos poder controlar dónde se publican las salidas al nivel de la configuración del flujo de trabajo.

"Oye", podría decir, "`CAT_CAT` está enviando sus salidas a `--outdir`. ¿Quizás deberíamos copiar su directiva `publishDir`?"

Sí, esa es una gran idea.

Excepto que no tiene una directiva `publishDir`. (Adelante, mire el código del módulo.)

Eso es porque los pipelines nf-core centralizan el control al nivel del flujo de trabajo configurando `publishDir` en `conf/modules.config` en lugar de en módulos individuales.
Específicamente, la plantilla nf-core declara una directiva `publishDir` predeterminada (con una estructura de directorio predefinida) que se aplica a todos los módulos a menos que se proporcione una directiva de sobrescritura.

¿No suena increíble? ¿Podría ser que para aprovechar esta directiva predeterminada, todo lo que necesitamos hacer es eliminar la directiva `publishDir` actual de nuestros módulos locales?

Probemos eso en `COWPY` para ver qué sucede, luego veremos el código para la configuración predeterminada para entender cómo funciona.

Finalmente, demostraremos cómo sobrescribir el comportamiento predeterminado si se desea.

#### 1.5.1. Eliminar la directiva `publishDir` de `COWPY`

Hagamos esto.
Abra el archivo del módulo `cowpy.nf` (en `core-hello/modules/local/`) y elimine la directiva `publishDir` como se muestra a continuación.

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

Echemos un vistazo a lo que sucede si ejecutamos el pipeline ahora.

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

Eche un vistazo a su directorio de trabajo actual.
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

Puede ver que Nextflow creó esta jerarquía de directorios basada en los nombres del flujo de trabajo y del módulo.

El código responsable vive en el archivo `conf/modules.config`.
Esta es la configuración `publishDir` predeterminada que forma parte de la plantilla nf-core y se aplica a todos los procesos:

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

- **`path:`** Determina el directorio de salida basado en el nombre del proceso.
  El nombre completo de un proceso contenido en `task.process` incluye la jerarquía de importaciones de flujo de trabajo y módulo (como `CORE_HELLO:HELLO:CAT_CAT`).
  Las operaciones `tokenize` eliminan esa jerarquía para obtener solo el nombre del proceso, luego toman la primera parte antes de cualquier guion bajo (si corresponde), y la convierten a minúsculas.
  Esto es lo que determina que los resultados de `CAT_CAT` se publiquen en `${params.outdir}/cat/`.
- **`mode:`** Controla cómo se publican los archivos (copia, enlace simbólico, etc.).
  Esto es configurable a través del parámetro `params.publish_dir_mode`.
- **`saveAs:`** Filtra qué archivos publicar.
  Este ejemplo excluye archivos `versions.yml` devolviendo `null` para ellos, evitando que se publiquen.

Esto proporciona una lógica consistente para organizar salidas.

La salida se ve aún mejor cuando todos los módulos en un pipeline adoptan esta convención, así que siéntase libre de ir a eliminar las directivas `publishDir` de los otros módulos en su pipeline.
Este valor predeterminado se aplicará incluso a módulos que no modificamos explícitamente para seguir las directrices de nf-core.

Dicho esto, puede decidir que desea organizar sus entradas de manera diferente, y la buena noticia es que es fácil hacerlo.

#### 1.5.3. Sobrescribir el valor predeterminado

Para sobrescribir la directiva `publishDir` predeterminada, simplemente puede agregar sus propias directivas al archivo `conf/modules.config`.

Por ejemplo, podría sobrescribir el valor predeterminado para un solo proceso usando el selector `withName:`, como en este ejemplo donde agregamos una directiva `publishDir` personalizada para el proceso 'COWPY'.

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

En realidad no vamos a hacer ese cambio, pero siéntase libre de jugar con esto y ver qué lógica puede implementar.

El punto es que este sistema permite le da lo mejor de ambos mundos: consistencia por defecto y la flexibilidad para personalizar la configuración bajo demanda.

Para resumir, obtiene:

- **Fuente única de verdad**: Toda la configuración de publicación vive en `modules.config`
- **Valor predeterminado útil**: Los procesos funcionan de inmediato sin configuración por módulo
- **Personalización fácil**: Sobrescriba el comportamiento de publicación en configuración, no en código de módulo
- **Módulos portables**: Los módulos no codifican ubicaciones de salida

Esto completa el conjunto de características de módulos nf-core que absolutamente debe aprender a usar, pero hay otras que puede leer en las [especificaciones de módulos nf-core](https://nf-co.re/docs/guidelines/components/modules).

### Conclusión

Ahora sabe cómo adaptar módulos locales para seguir las convenciones de nf-core:

- Diseñe sus módulos para aceptar y propagar tuplas de metadatos;
- Use `ext.args` para mantener las interfaces de módulos mínimas y portables;
- Use `ext.prefix` para nombres de archivos de salida configurables y estandarizados;
- Adopte la directiva `publishDir` centralizada predeterminada para una estructura de directorio de resultados consistente.

### ¿Qué sigue?

Aprenda a usar las herramientas integradas basadas en plantillas de nf-core para crear módulos de la manera fácil.

---

## 2. Crear un módulo con las herramientas nf-core

Ahora que ha aprendido los patrones de módulos nf-core aplicándolos manualmente, veamos cómo crearía módulos en la práctica.

### 2.1. Generar un esqueleto de módulo desde una plantilla

Similar a lo que existe para crear pipelines, el proyecto nf-core proporciona herramientas para generar módulos estructurados correctamente basados en una plantilla, con todos estos patrones incorporados desde el principio.

#### 2.1.1. Ejecutar el comando de creación de módulo

El comando `nf-core modules create` genera una plantilla de módulo que ya sigue todas las convenciones que ha aprendido.

Creemos una nueva versión del módulo `COWPY` con una plantilla mínima ejecutando este comando:

```bash
nf-core modules create --empty-template COWPY
```

La bandera `--empty-template` crea una plantilla inicial limpia sin código extra, facilitando ver la estructura esencial.

El comando se ejecuta de forma interactiva, guiándolo a través de la configuración.
Busca automáticamente información de herramientas de repositorios de paquetes como Bioconda y bio.tools para prepoblar metadatos.

Se le solicitarán varias opciones de configuración:

- **Información del autor**: Su nombre de usuario de GitHub para atribución
- **Etiqueta de recurso**: Un conjunto predefinido de requisitos computacionales.
  El proyecto nf-core proporciona etiquetas estándar como `process_single` para herramientas ligeras y `process_high` para las exigentes.
  Estas etiquetas ayudan a gestionar la asignación de recursos en diferentes entornos de ejecución.
- **Requisito de metadatos**: Si el módulo necesita información específica de muestra a través de un mapa `meta` (generalmente sí para módulos de procesamiento de datos).

La herramienta maneja la complejidad de encontrar información de paquetes y configurar la estructura, permitiéndole concentrarse en implementar la lógica específica de la herramienta.

#### 2.1.2. Examinar el esqueleto del módulo

La herramienta crea una estructura de módulo completa en `modules/local/` (o `modules/nf-core/` si está en el repositorio nf-core/modules):

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

- **`main.nf`**: Definición de proceso con todos los patrones nf-core incorporados
- **`meta.yml`**: Documentación del módulo describiendo entradas, salidas y la herramienta
- **`environment.yml`**: Especificación de entorno Conda para dependencias
- **`tests/main.nf.test`**: Casos de prueba nf-test para validar que el módulo funciona

!!! tip "Aprenda más sobre pruebas"

    El archivo de prueba generado usa nf-test, un framework de pruebas para pipelines y módulos de Nextflow. Para aprender cómo escribir y ejecutar estas pruebas, consulte la [misión secundaria nf-test](../side_quests/nf-test.md).

El `main.nf` generado incluye todos los patrones que acaba de aprender, más algunas características adicionales:

```groovy title="modules/local/cowpy/main.nf" hl_lines="11 21 22"
process COWPY {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/YOUR-TOOL-HERE':
        'biocontainers/
```
