# Depuración de Flujos de Trabajo

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Traducción asistida por IA - [más información y sugerencias](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

La depuración es una habilidad crítica que puede ahorrarle horas de frustración y ayudarle a convertirse en un desarrollador de Nextflow más eficaz. A lo largo de su carrera, especialmente cuando está comenzando, encontrará errores mientras construye y mantiene sus flujos de trabajo. Aprender enfoques sistemáticos de depuración le ayudará a identificar y resolver problemas rápidamente.

### Objetivos de aprendizaje

En esta misión secundaria, exploraremos **técnicas sistemáticas de depuración** para flujos de trabajo de Nextflow:

- **Depuración de errores de sintaxis**: Uso efectivo de características del IDE y mensajes de error de Nextflow
- **Depuración de canales**: Diagnóstico de problemas de flujo de datos y problemas de estructura de canales
- **Depuración de procesos**: Investigación de fallas de ejecución y problemas de recursos
- **Herramientas de depuración integradas**: Aprovechamiento del modo preview, ejecución stub y directorios de trabajo de Nextflow
- **Enfoques sistemáticos**: Una metodología de cuatro fases para depuración eficiente

Al final, tendrá una metodología robusta de depuración que transforma mensajes de error frustrantes en hojas de ruta claras hacia soluciones.

### Requisitos previos

Antes de emprender esta misión secundaria, debería:

- Haber completado el tutorial [Hello Nextflow](../hello_nextflow/README.md) o un curso equivalente para principiantes.
- Sentirse cómodo usando conceptos y mecanismos básicos de Nextflow (procesos, canales, operadores)

**Opcional:** Recomendamos completar primero la misión secundaria [Características del IDE para Desarrollo con Nextflow](./ide_features.md).
Esa cubre características completas del IDE que apoyan la depuración (resaltado de sintaxis, detección de errores, etc.), que usaremos intensivamente aquí.

---

## 0. Comenzar

#### Abrir el codespace de entrenamiento

Si aún no lo ha hecho, asegúrese de abrir el entorno de entrenamiento como se describe en [Configuración del Entorno](../envsetup/index.md).

[![Abrir en GitHub Codespaces](https://github.com/codespaces/badge.svg)](https://codespaces.new/nextflow-io/training?quickstart=1&ref=master)

#### Moverse al directorio del proyecto

Vamos a movernos al directorio donde se encuentran los archivos para este tutorial.

```bash
cd side-quests/debugging
```

Puede configurar VSCode para enfocarse en este directorio:

```bash
code .
```

#### Revisar los materiales

Encontrará un conjunto de flujos de trabajo de ejemplo con varios tipos de errores que usaremos para practicar:

??? abstract "Contenido del directorio"

    ```console
    .
    ├── bad_bash_var.nf
    ├── bad_channel_shape.nf
    ├── bad_channel_shape_viewed_debug.nf
    ├── bad_channel_shape_viewed.nf
    ├── bad_number_inputs.nf
    ├── badpractice_syntax.nf
    ├── bad_resources.nf
    ├── bad_syntax.nf
    ├── buggy_workflow.nf
    ├── data
    │   ├── sample_001.fastq.gz
    │   ├── sample_002.fastq.gz
    │   ├── sample_003.fastq.gz
    │   ├── sample_004.fastq.gz
    │   ├── sample_005.fastq.gz
    │   └── sample_data.csv
    ├── exhausted.nf
    ├── invalid_process.nf
    ├── missing_output.nf
    ├── missing_software.nf
    ├── missing_software_with_stub.nf
    ├── nextflow.config
    └── no_such_var.nf
    ```

Estos archivos representan escenarios comunes de depuración que encontrará en el desarrollo del mundo real.

#### Revisar la asignación

Su desafío es ejecutar cada flujo de trabajo, identificar el(los) error(es) y corregirlos.

Para cada flujo de trabajo con errores:

1. **Ejecutar el flujo de trabajo** y observar el error
2. **Analizar el mensaje de error**: ¿qué le está diciendo Nextflow?
3. **Localizar el problema** en el código usando las pistas proporcionadas
4. **Corregir el error** y verificar que su solución funciona
5. **Restablecer el archivo** antes de pasar a la siguiente sección (use `git checkout <filename>`)

Los ejercicios progresan desde errores de sintaxis simples hasta problemas de tiempo de ejecución más sutiles.
Las soluciones se discuten en línea, pero intente resolver cada uno usted mismo antes de leer más adelante.

#### Lista de verificación de preparación

¿Cree que está listo para comenzar?

- [ ] Entiendo el objetivo de este curso y sus requisitos previos
- [ ] Mi codespace está funcionando
- [ ] He establecido mi directorio de trabajo apropiadamente
- [ ] Entiendo la asignación

Si puede marcar todas las casillas, está listo para comenzar.

---

## 1. Errores de Sintaxis

Los errores de sintaxis son el tipo más común de error que encontrará al escribir código Nextflow. Ocurren cuando el código no se ajusta a las reglas de sintaxis esperadas del DSL de Nextflow. Estos errores evitan que su flujo de trabajo se ejecute en absoluto, por lo que es importante aprender a identificarlos y corregirlos rápidamente.

### 1.1. Llaves faltantes

Uno de los errores de sintaxis más comunes, y a veces uno de los más complejos de depurar, son **corchetes faltantes o desemparejados**.

Comencemos con un ejemplo práctico.

#### Ejecutar el pipeline

```bash
nextflow run bad_syntax.nf
```

??? failure "Salida del comando"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `bad_syntax.nf` [stupefied_bhabha] DSL2 - revision: ca6327fad2

    Error bad_syntax.nf:24:1: Unexpected input: '<EOF>'

    ERROR ~ Script compilation failed

     -- Check '.nextflow.log' file for details
    ```

**Elementos clave de los mensajes de error de sintaxis:**

- **Archivo y ubicación**: Muestra qué archivo y línea/columna contienen el error (`bad_syntax.nf:24:1`)
- **Descripción del error**: Explica lo que el analizador encontró que no esperaba (`Unexpected input: '<EOF>'`)
- **Indicador EOF**: El mensaje `<EOF>` (End Of File - Fin de Archivo) indica que el analizador llegó al final del archivo mientras todavía esperaba más contenido - una señal clásica de llaves sin cerrar

#### Verificar el código

Ahora, examinemos `bad_syntax.nf` para entender qué está causando el error:

```groovy title="bad_syntax.nf" hl_lines="14" linenums="1"
#!/usr/bin/env nextflow

process PROCESS_FILES {
    input:
    val sample_name

    output:
    path "${sample_name}_output.txt"

    script:
    """
    echo "Processing ${sample_name}" > ${sample_name}_output.txt
    """
// Falta la llave de cierre para el proceso

workflow {

    // Crear canal de entrada
    input_ch = channel.of('sample1', 'sample2', 'sample3')

    // Llamar al proceso con el canal de entrada
    PROCESS_FILES(input_ch)
}
```

Para el propósito de este ejemplo, hemos dejado un comentario para mostrarle dónde está el error. La extensión de Nextflow para VSCode también debería estar dándole algunas pistas sobre lo que podría estar mal, poniendo la llave desemparejada en rojo y resaltando el final prematuro del archivo:

![Bad syntax](img/bad_syntax.png)

**Estrategia de depuración para errores de corchetes:**

1. Use el emparejamiento de corchetes de VS Code (coloque el cursor junto a un corchete)
2. Revise el panel de Problemas para mensajes relacionados con corchetes
3. Asegúrese de que cada `{` de apertura tenga su correspondiente `}` de cierre

#### Corregir el código

Reemplace el comentario con la llave de cierre faltante:

=== "Después"

    ```groovy title="bad_syntax.nf" hl_lines="14" linenums="1"
    #!/usr/bin/env nextflow

    process PROCESS_FILES {
        input:
        val sample_name

        output:
        path "${sample_name}_output.txt"

        script:
        """
        echo "Processing ${sample_name}" > ${sample_name}_output.txt
        """
    }  // Agregar la llave de cierre faltante

    workflow {

        // Crear canal de entrada
        input_ch = channel.of('sample1', 'sample2', 'sample3')

        // Llamar al proceso con el canal de entrada
        PROCESS_FILES(input_ch)
    }
    ```

=== "Antes"

    ```groovy title="bad_syntax.nf" hl_lines="14" linenums="1"
    #!/usr/bin/env nextflow

    process PROCESS_FILES {
        input:
        val sample_name

        output:
        path "${sample_name}_output.txt"

        script:
        """
        echo "Processing ${sample_name}" > ${sample_name}_output.txt
        """
    // Falta la llave de cierre para el proceso

    workflow {

        // Crear canal de entrada
        input_ch = channel.of('sample1', 'sample2', 'sample3')

        // Llamar al proceso con el canal de entrada
        PROCESS_FILES(input_ch)
    }
    ```

#### Ejecutar el pipeline

Ahora ejecute el flujo de trabajo nuevamente para confirmar que funciona:

```bash
nextflow run bad_syntax.nf
```

??? success "Salida del comando"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `bad_syntax.nf` [insane_faggin] DSL2 - revision: 961938ee2b

    executor >  local (3)
    [48/cd7f54] PROCESS_FILES (1) | 3 of 3 ✔
    ```

### 1.2. Uso de palabras clave o directivas de proceso incorrectas

Otro error de sintaxis común es una **definición de proceso inválida**. Esto puede suceder si olvida definir bloques requeridos o usa directivas incorrectas en la definición del proceso.

#### Ejecutar el pipeline

```bash
nextflow run invalid_process.nf
```

??? failure "Salida del comando"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `invalid_process.nf` [nasty_jepsen] DSL2 - revision: da9758d614

    Error invalid_process.nf:3:1: Invalid process definition -- check for missing or out-of-order section labels
    │   3 | process PROCESS_FILES {
    │     | ^^^^^^^^^^^^^^^^^^^^^^^
    │   4 |     inputs:
    │   5 |     val sample_name
    │   6 |
    ╰   7 |     output:

    ERROR ~ Script compilation failed

     -- Check '.nextflow.log' file for details
    ```

#### Verificar el código

El error indica una "Definición de proceso inválida" y muestra el contexto alrededor del problema. Mirando las líneas 3-7, podemos ver `inputs:` en la línea 4, que es el problema. Examinemos `invalid_process.nf`:

```groovy title="invalid_process.nf" hl_lines="4" linenums="1"
#!/usr/bin/env nextflow

process PROCESS_FILES {
    inputs:  // ERROR: Debería ser 'input' no 'inputs'
    val sample_name

    output:
    path "${sample_name}_output.txt"

    script:
    """
    echo "Processing ${sample_name}" > ${sample_name}_output.txt
    """
}

workflow {

    // Crear canal de entrada
    input_ch = channel.of('sample1', 'sample2', 'sample3')

    // Llamar al proceso con el canal de entrada
    PROCESS_FILES(input_ch)
}
```

Mirando la línea 4 en el contexto del error, podemos identificar el problema: estamos usando `inputs` en lugar de la directiva correcta `input`. La extensión de Nextflow para VSCode también marcará esto:

![Invalid process message](img/invalid_process_message.png)

#### Corregir el código

Reemplace la palabra clave incorrecta con la correcta consultando [la documentación](https://www.nextflow.io/docs/latest/process.html#):

=== "Después"

    ```groovy title="invalid_process.nf" hl_lines="4" linenums="1"
    #!/usr/bin/env nextflow

    process PROCESS_FILES {
        input:  // Corregido: Cambiado 'inputs' a 'input'
        val sample_name

        output:
        path "${sample_name}_output.txt"

        script:
        """
        echo "Processing ${sample_name}" > ${sample_name}_output.txt
        """
    }

    workflow {

        // Crear canal de entrada
        input_ch = channel.of('sample1', 'sample2', 'sample3')

        // Llamar al proceso con el canal de entrada
        PROCESS_FILES(input_ch)
    }
    ```

=== "Antes"

    ```groovy title="invalid_process.nf" hl_lines="4" linenums="1"
    #!/usr/bin/env nextflow

    process PROCESS_FILES {
        inputs:  // ERROR: Debería ser 'input' no 'inputs'
        val sample_name

        output:
        path "${sample_name}_output.txt"

        script:
        """
        echo "Processing ${sample_name}" > ${sample_name}_output.txt
        """
    }

    workflow {

        // Crear canal de entrada
        input_ch = channel.of('sample1', 'sample2', 'sample3')

        // Llamar al proceso con el canal de entrada
        PROCESS_FILES(input_ch)
    }
    ```

#### Ejecutar el pipeline

Ahora ejecute el flujo de trabajo nuevamente para confirmar que funciona:

```bash
nextflow run invalid_process.nf
```

??? success "Salida del comando"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `invalid_process.nf` [silly_fermi] DSL2 - revision: 961938ee2b

    executor >  local (3)
    [b7/76cd9d] PROCESS_FILES (2) | 3 of 3 ✔
    ```

### 1.3. Uso de nombres de variable incorrectos

Los nombres de variable que usa en sus bloques de script deben ser válidos, derivados ya sea de entradas o de código groovy insertado antes del script. Pero cuando está manejando complejidad al inicio del desarrollo del pipeline, es fácil cometer errores en el nombramiento de variables, y Nextflow se lo hará saber rápidamente.

#### Ejecutar el pipeline

```bash
nextflow run no_such_var.nf
```

??? failure "Salida del comando"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `no_such_var.nf` [gloomy_meninsky] DSL2 - revision: 0c4d3bc28c

    Error no_such_var.nf:17:39: `undefined_var` is not defined
    │  17 |     echo "Using undefined variable: ${undefined_var}" >> ${output_pref
    ╰     |                                       ^^^^^^^^^^^^^

    ERROR ~ Script compilation failed

     -- Check '.nextflow.log' file for details
    ```

El error se detecta en tiempo de compilación y apunta directamente a la variable no definida en la línea 17, con un acento circunflejo indicando exactamente dónde está el problema.

#### Verificar el código

Examinemos `no_such_var.nf`:

```groovy title="no_such_var.nf" hl_lines="17" linenums="1"
#!/usr/bin/env nextflow

process PROCESS_FILES {
    input:
    val sample_name

    output:
    path "${sample_name}_processed.txt"

    script:
    // Definir variables en código Groovy antes del script
    def output_prefix = "${sample_name}_processed"
    def timestamp = new Date().format("yyyy-MM-dd")

    """
    echo "Processing ${sample_name} on ${timestamp}" > ${output_prefix}.txt
    echo "Using undefined variable: ${undefined_var}" >> ${output_prefix}.txt  // ERROR: undefined_var no definida
    """
}

workflow {
    input_ch = channel.of('sample1', 'sample2', 'sample3')
    PROCESS_FILES(input_ch)
}
```

El mensaje de error indica que la variable no se reconoce en la plantilla del script, y ahí está: puede ver `${undefined_var}` usado en el bloque de script, pero no definido en otro lugar.

#### Corregir el código

Si obtiene un error de 'No existe tal variable', puede corregirlo definiendo la variable (corrigiendo nombres de variables de entrada o editando código groovy antes del script), o eliminándola del bloque de script si no es necesaria:

=== "Después"

    ```groovy title="no_such_var.nf" hl_lines="15-17" linenums="1"
    #!/usr/bin/env nextflow

    process PROCESS_FILES {
        input:
        val sample_name

        output:
        path "${sample_name}_output.txt"

        script:
        // Definir variables en código Groovy antes del script
        def output_prefix = "${sample_name}_processed"
        def timestamp = new Date().format("yyyy-MM-dd")

        """
        echo "Processing ${sample_name} on ${timestamp}" > ${output_prefix}.txt
        """  // Eliminada la línea con undefined_var
    }

    workflow {
        input_ch = channel.of('sample1', 'sample2', 'sample3')
        PROCESS_FILES(input_ch)
    }
    ```

=== "Antes"

    ```groovy title="no_such_var.nf" hl_lines="17" linenums="1"
    #!/usr/bin/env nextflow

    process PROCESS_FILES {
        input:
        val sample_name

        output:
        path "${sample_name}_output.txt"

        script:
        // Definir variables en código Groovy antes del script
        def output_prefix = "${sample_name}_processed"
        def timestamp = new Date().format("yyyy-MM-dd")

        """
        echo "Processing ${sample_name} on ${timestamp}" > ${output_prefix}.txt
        echo "Using undefined variable: ${undefined_var}" >> ${output_prefix}.txt  // ERROR: undefined_var no definida
        """
    }

    workflow {
        input_ch = channel.of('sample1', 'sample2', 'sample3')
        PROCESS_FILES(input_ch)
    }
    ```

#### Ejecutar el pipeline

Ahora ejecute el flujo de trabajo nuevamente para confirmar que funciona:

```bash
nextflow run no_such_var.nf
```

??? success "Salida del comando"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `no_such_var.nf` [suspicious_venter] DSL2 - revision: 6ba490f7c5

    executor >  local (3)
    [21/237300] PROCESS_FILES (2) | 3 of 3 ✔
    ```

### 1.4. Mal uso de variables de Bash

Comenzando en Nextflow, puede ser difícil entender la diferencia entre variables de Nextflow (Groovy) y Bash. Esto puede generar otra forma del error de variable incorrecta que aparece al intentar usar variables en el contenido Bash del bloque de script.

#### Ejecutar el pipeline

```bash
nextflow run bad_bash_var.nf
```

??? failure "Salida del comando"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `bad_bash_var.nf` [infallible_mandelbrot] DSL2 - revision: 0853c11080

    Error bad_bash_var.nf:13:42: `prefix` is not defined
    │  13 |     echo "Processing ${sample_name}" > ${prefix}.txt
    ╰     |                                          ^^^^^^

    ERROR ~ Script compilation failed

     -- Check '.nextflow.log' file for details
    ```

#### Verificar el código

El error apunta a la línea 13 donde se usa `${prefix}`. Examinemos `bad_bash_var.nf` para ver qué está causando el problema:

```groovy title="bad_bash_var.nf" hl_lines="13" linenums="1"
#!/usr/bin/env nextflow

process PROCESS_FILES {
    input:
    val sample_name

    output:
    path "${sample_name}_output.txt"

    script:
    """
    prefix="${sample_name}_output"
    echo "Processing ${sample_name}" > ${prefix}.txt  # ERROR: ${prefix} es sintaxis Groovy, no Bash
    """
}
```

En este ejemplo, estamos definiendo la variable `prefix` en Bash, pero en un proceso de Nextflow la sintaxis `$` que usamos para referirnos a ella (`${prefix}`) se interpreta como una variable Groovy, no Bash. La variable no existe en el contexto Groovy, por lo que obtenemos un error de 'no existe tal variable'.

#### Corregir el código

Si quiere usar una variable de Bash, debe escapar el signo de dólar así:

=== "Después"

    ```groovy title="bad_bash_var.nf" hl_lines="13" linenums="1"
    #!/usr/bin/env nextflow

    process PROCESS_FILES {
        input:
        val sample_name

        output:
        path "${sample_name}_output.txt"

        script:
        """
        prefix="${sample_name}_output"
        echo "Processing ${sample_name}" > \${prefix}.txt  # Corregido: Escapado el signo de dólar
        """
    }

    workflow {
        input_ch = channel.of('sample1', 'sample2', 'sample3')
        PROCESS_FILES(input_ch)
    }
    ```

=== "Antes"

    ```groovy title="bad_bash_var.nf" hl_lines="13" linenums="1"
    #!/usr/bin/env nextflow

    process PROCESS_FILES {
        input:
        val sample_name

        output:
        path "${sample_name}_output.txt"

        script:
        """
        prefix="${sample_name}_output"
        echo "Processing ${sample_name}" > ${prefix}.txt  # ERROR: ${prefix} es sintaxis Groovy, no Bash
        """
    }
    ```

Esto le dice a Nextflow que interprete esto como una variable de Bash.

#### Ejecutar el pipeline

Ahora ejecute el flujo de trabajo nuevamente para confirmar que funciona:

```bash
nextflow run bad_bash_var.nf
```

??? success "Salida del comando"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `bad_bash_var.nf` [naughty_franklin] DSL2 - revision: 58c1c83709

    executor >  local (3)
    [4e/560285] PROCESS_FILES (2) | 3 of 3 ✔
    ```

!!! tip "Variables Groovy vs Bash"

    Para manipulaciones de variables simples como concatenación de strings u operaciones de prefijo/sufijo, generalmente es más legible usar variables Groovy en la sección de script en lugar de variables Bash en el bloque de script:

    ```groovy linenums="1"
    script:
    def output_prefix = "${sample_name}_processed"
    def output_file = "${output_prefix}.txt"
    """
    echo "Processing ${sample_name}" > ${output_file}
    """
    ```

    Este enfoque evita la necesidad de escapar signos de dólar y hace que el código sea más fácil de leer y mantener.

### 1.5. Declaraciones Fuera del Bloque Workflow

La extensión de Nextflow para VSCode resalta problemas con la estructura del código que causarán errores. Un ejemplo común es definir canales fuera del bloque `workflow {}` - esto ahora se aplica como un error de sintaxis.

#### Ejecutar el pipeline

```bash
nextflow run badpractice_syntax.nf
```

??? failure "Salida del comando"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `badpractice_syntax.nf` [intergalactic_colden] DSL2 - revision: 5e4b291bde

    Error badpractice_syntax.nf:3:1: Statements cannot be mixed with script declarations -- move statements into a process or workflow
    │   3 | input_ch = channel.of('sample1', 'sample2', 'sample3')
    ╰     | ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

    ERROR ~ Script compilation failed

     -- Check '.nextflow.log' file for details
    ```

El mensaje de error indica claramente el problema: las declaraciones (como definiciones de canales) no pueden mezclarse con declaraciones de script fuera de un bloque workflow o process.

#### Verificar el código

Examinemos `badpractice_syntax.nf` para ver qué está causando el error:

```groovy title="badpractice_syntax.nf" hl_lines="3" linenums="1"
#!/usr/bin/env nextflow

input_ch = channel.of('sample1', 'sample2', 'sample3')  // ERROR: Canal definido fuera del workflow

process PROCESS_FILES {
    input:
    val sample_name

    output:
    path "${sample_name}_processed.txt"

    script:
    // Definir variables en código Groovy antes del script
    def output_prefix = "${sample_name}_processed"
    def timestamp = new Date().format("yyyy-MM-dd")

    """
    echo "Processing ${sample_name} on ${timestamp}" > ${output_prefix}.txt
    """
}

workflow {
    PROCESS_FILES(input_ch)
}
```

La extensión de VSCode también resaltará la variable `input_ch` como definida fuera del bloque workflow:

![Non-lethal syntax error](img/nonlethal.png)

#### Corregir el código

Mueva la definición del canal dentro del bloque workflow:

=== "Después"

    ```groovy title="badpractice_syntax.nf" hl_lines="21" linenums="1"
    #!/usr/bin/env nextflow

    process PROCESS_FILES {
        input:
        val sample_name

        output:
        path "${sample_name}_processed.txt"

        script:
        // Definir variables en código Groovy antes del script
        def output_prefix = "${sample_name}_processed"
        def timestamp = new Date().format("yyyy-MM-dd")

        """
        echo "Processing ${sample_name} on ${timestamp}" > ${output_prefix}.txt
        """
    }

    workflow {
        input_ch = channel.of('sample1', 'sample2', 'sample3')  // Movido dentro del bloque workflow
        PROCESS_FILES(input_ch)
    }
    ```

=== "Antes"

    ```groovy title="badpractice_syntax.nf" hl_lines="3" linenums="1"
    #!/usr/bin/env nextflow

    input_ch = channel.of('sample1', 'sample2', 'sample3')  // ERROR: Canal definido fuera del workflow

    process PROCESS_FILES {
        input:
        val sample_name

        output:
        path "${sample_name}_processed.txt"

        script:
        // Definir variables en código Groovy antes del script
        def output_prefix = "${sample_name}_processed"
        def timestamp = new Date().format("yyyy-MM-dd")

        """
        echo "Processing ${sample_name} on ${timestamp}" > ${output_prefix}.txt
        """
    }

    workflow {
        PROCESS_FILES(input_ch)
    }
    ```

#### Ejecutar el pipeline

Ejecute el flujo de trabajo nuevamente para confirmar que la corrección funciona:

```bash
nextflow run badpractice_syntax.nf
```

??? success "Salida del comando"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `badpractice_syntax.nf` [naughty_ochoa] DSL2 - revision: 5e4b291bde

    executor >  local (3)
    [6a/84a608] PROCESS_FILES (2) | 3 of 3 ✔
    ```

Mantenga sus canales de entrada definidos dentro del bloque workflow, y en general siga cualquier otra recomendación que haga la extensión.

### Conclusión

Puede identificar y corregir errores de sintaxis sistemáticamente usando mensajes de error de Nextflow e indicadores visuales del IDE. Los errores de sintaxis comunes incluyen llaves faltantes, palabras clave de proceso incorrectas, variables no definidas y uso inadecuado de variables de Bash vs. Nextflow. La extensión de VSCode ayuda a detectar muchos de estos antes del tiempo de ejecución. Con estas habilidades de depuración de sintaxis en su caja de herramientas, podrá resolver rápidamente los errores de sintaxis más comunes de Nextflow y pasar a abordar problemas de tiempo de ejecución más complejos.

### ¿Qué sigue?

Aprenda a depurar errores de estructura de canal más complejos que ocurren incluso cuando la sintaxis es correcta.

---

## 2. Errores de Estructura de Canal

Los errores de estructura de canal son más sutiles que los errores de sintaxis porque el código es sintácticamente correcto, pero las formas de los datos no coinciden con lo que los procesos esperan. Nextflow intentará ejecutar el pipeline, pero podría encontrar que el número de entradas no coincide con lo que espera y fallar. Estos errores típicamente solo aparecen en tiempo de ejecución y requieren una comprensión de los datos que fluyen a través de su flujo de trabajo.

!!! tip "Depuración de Canales con `.view()`"

    A lo largo de esta sección, recuerde que puede usar el operador `.view()` para inspeccionar el contenido del canal en cualquier punto de su flujo de trabajo. Esta es una de las herramientas de depuración más poderosas para entender problemas de estructura de canal. Exploraremos esta técnica en detalle en la sección 2.4, pero siéntase libre de usarla mientras trabaja en los ejemplos.

    ```groovy
    my_channel.view()  // Muestra lo que está fluyendo a través del canal
    ```

### 2.1. Número Incorrecto de Canales de Entrada

Este error ocurre cuando pasa un número diferente de canales del que un proceso espera.

#### Ejecutar el pipeline

```bash
nextflow run bad_number_inputs.nf
```

??? failure "Salida del comando"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `bad_number_inputs.nf` [happy_swartz] DSL2 - revision: d83e58dcd3

    Error bad_number_inputs.nf:23:5: Incorrect number of call arguments, expected 1 but received 2
    │  23 |     PROCESS_FILES(samples_ch, files_ch)
    ╰     |     ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

    ERROR ~ Script compilation failed

     -- Check '.nextflow.log' file for details
    ```

#### Verificar el código

El mensaje de error indica claramente que la llamada esperaba 1 argumento pero recibió 2, y apunta a la línea 23. Examinemos `bad_number_inputs.nf`:

```groovy title="bad_number_inputs.nf" hl_lines="5 23" linenums="1"
#!/usr/bin/env nextflow

process PROCESS_FILES {
    input:
        val sample_name  // El proceso espera solo 1 entrada

    output:
        path "${sample_name}_output.txt"

    script:
    """
    echo "Processing ${sample_name}" > ${sample_name}_output.txt
    """
}

workflow {

    // Crear dos canales separados
    samples_ch = channel.of('sample1', 'sample2', 'sample3')
    files_ch = channel.of('file1.txt', 'file2.txt', 'file3.txt')

    // ERROR: Pasando 2 canales pero el proceso espera solo 1
    PROCESS_FILES(samples_ch, files_ch)
}
```

Debería ver la llamada desemparejada `PROCESS_FILES`, suministrando múltiples canales de entrada cuando el proceso solo define uno. La extensión de VSCode también subrayará la llamada al proceso en rojo y suministrará un mensaje de diagnóstico cuando pase el mouse:

![Incorrect number of args message](img/incorrect_num_args.png)

#### Corregir el código

Para este ejemplo específico, el proceso espera un solo canal y no requiere el segundo canal, por lo que podemos corregirlo pasando solo el canal `samples_ch`:

=== "Después"

    ```groovy title="bad_number_inputs.nf" hl_lines="23" linenums="1"
    #!/usr/bin/env nextflow

    process PROCESS_FILES {
        input:
            val sample_name  // El proceso espera solo 1 entrada

        output:
            path "${sample_name}_output.txt"

        script:
        """
        echo "Processing ${sample_name}" > ${sample_name}_output.txt
        """
    }

    workflow {

        // Crear dos canales separados
        samples_ch = channel.of('sample1', 'sample2', 'sample3')
        files_ch = channel.of('file1.txt', 'file2.txt', 'file3.txt')

        // Corregido: Pasar solo el canal que el proceso espera
        PROCESS_FILES(samples_ch)
    }
    ```

=== "Antes"

    ```groovy title="bad_number_inputs.nf" hl_lines="5 23" linenums="1"
    #!/usr/bin/env nextflow

    process PROCESS_FILES {
        input:
            val sample_name  // El proceso espera solo 1 entrada

        output:
            path "${sample_name}_output.txt"

        script:
        """
        echo "Processing ${sample_name}" > ${sample_name}_output.txt
        """
    }

    workflow {

        // Crear dos canales separados
        samples_ch = channel.of('sample1', 'sample2', 'sample3')
        files_ch = channel.of('file1.txt', 'file2.txt', 'file3.txt')

        // ERROR: Pasando 2 canales pero el proceso espera solo 1
        PROCESS_FILES(samples_ch, files_ch)
    }
    ```

#### Ejecutar el pipeline

```bash
nextflow run bad_number_inputs.nf
```

??? success "Salida del comando"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `bad_number_inputs.nf` [big_euler] DSL2 - revision: e302bd87be

    executor >  local (3)
    [48/497f7b] PROCESS_FILES (3) | 3 of 3 ✔
    ```

Más comúnmente que este ejemplo, podría agregar entradas adicionales a un proceso y olvidar actualizar la llamada del workflow en consecuencia, lo que puede conducir a este tipo de error. Afortunadamente, este es uno de los errores más fáciles de entender y corregir, ya que el mensaje de error es bastante claro sobre el desajuste.

### 2.2. Agotamiento de Canal (El Proceso se Ejecuta Menos Veces de lo Esperado)

Algunos errores de estructura de canal son mucho más sutiles y no producen errores en absoluto. Probablemente el más común de estos refleja un desafío que los nuevos usuarios de Nextflow enfrentan al entender que los canales de cola pueden agotarse y quedarse sin elementos, lo que significa que el flujo de trabajo termina prematuramente.

#### Ejecutar el pipeline

```bash
nextflow run exhausted.nf
```

??? success "Salida del comando"

```console title="Salida de canal agotado"
 N E X T F L O W   ~  version 25.10.2

Launching `exhausted.nf` [extravagant_gauss] DSL2 - revision: 08cff7ba2a

executor >  local (1)
[bd/f61fff] PROCESS_FILES (1) [100%] 1 of 1 ✔
```

¡Este flujo de trabajo se completa sin errores, pero solo procesa una sola muestra!

#### Verificar el código

Examinemos `exhausted.nf` para ver si eso es correcto:

```groovy title="exhausted.nf" hl_lines="23 24" linenums="1"
#!/usr/bin/env nextflow

process PROCESS_FILES {
    input:
    val reference
    val sample_name

    output:
    path "${output_prefix}.txt"

    script:
    // Definir variables en código Groovy antes del script
    output_prefix = "${reference}_${sample_name}"
    def timestamp = new Date().format("yyyy-MM-dd")

    """
    echo "Processing ${sample_name} on ${timestamp}" > ${output_prefix}.txt
    """
}

workflow {

    reference_ch = channel.of('baseline_reference')
    input_ch = channel.of('sample1', 'sample2', 'sample3')

    PROCESS_FILES(reference_ch, input_ch)
}
```

El proceso solo se ejecuta una vez en lugar de tres veces porque el canal `reference_ch` es un canal de cola que se agota después de la primera ejecución del proceso. Cuando un canal se agota, todo el proceso se detiene, incluso si otros canales todavía tienen elementos.

Este es un patrón común donde tiene un archivo de referencia único que necesita ser reutilizado en múltiples muestras. La solución es convertir el canal de referencia en un canal de valor que puede ser reutilizado indefinidamente.

#### Corregir el código

Hay un par de formas de abordar esto dependiendo de cuántos archivos están afectados.

**Opción 1**: Tiene un solo archivo de referencia que está reutilizando mucho. Puede simplemente crear un tipo de canal de valor, que puede usarse una y otra vez. Hay tres formas de hacer esto:

**1a** Usar `channel.value()`:

```groovy title="exhausted.nf (corregido - Opción 1a)" hl_lines="2" linenums="21"
workflow {
    reference_ch = channel.value('baseline_reference')  // El canal de valor puede ser reutilizado
    input_ch = channel.of('sample1', 'sample2', 'sample3')

    PROCESS_FILES(reference_ch, input_ch)
}
```

**1b** Usar el [operador](https://www.nextflow.io/docs/latest/reference/operator.html#first) `first()`:

```groovy title="exhausted.nf (corregido - Opción 1b)" hl_lines="2" linenums="21"
workflow {
    reference_ch = channel.of('baseline_reference').first()  // Convertir a canal de valor
    input_ch = channel.of('sample1', 'sample2', 'sample3')

    PROCESS_FILES(reference_ch, input_ch)
}
```

**1c.** Usar el [operador](https://www.nextflow.io/docs/latest/reference/operator.html#collect) `collect()`:

```groovy title="exhausted.nf (corregido - Opción 1c)" hl_lines="2" linenums="21"
workflow {
    reference_ch = channel.of('baseline_reference').collect()  // Convertir a canal de valor
    input_ch = channel.of('sample1', 'sample2', 'sample3')

    PROCESS_FILES(reference_ch, input_ch)
}
```

**Opción 2**: En escenarios más complejos, quizás donde tiene múltiples archivos de referencia para todas las muestras en el canal de muestras, puede usar el operador `combine` para crear un nuevo canal que combine los dos canales en tuplas:

```groovy title="exhausted.nf (corregido - Opción 2)" hl_lines="4" linenums="21"
workflow {
    reference_ch = channel.of('baseline_reference','other_reference')
    input_ch = channel.of('sample1', 'sample2', 'sample3')
    combined_ch = reference_ch.combine(input_ch)  // Crea producto cartesiano

    PROCESS_FILES(combined_ch)
}
```

El operador `.combine()` genera un producto cartesiano de los dos canales, por lo que cada elemento en `reference_ch` se emparejará con cada elemento en `input_ch`. Esto permite que el proceso se ejecute para cada muestra mientras sigue usando la referencia.

Esto requiere que la entrada del proceso sea ajustada. En nuestro ejemplo, el inicio de la definición del proceso necesitaría ser ajustado de la siguiente manera:

```groovy title="exhausted.nf (corregido - Opción 2)" hl_lines="5" linenums="1"
#!/usr/bin/env nextflow

process PROCESS_FILES {
    input:
        tuple val(reference), val(sample_name)
```

Este enfoque puede no ser adecuado en todas las situaciones.

#### Ejecutar el pipeline

Pruebe una de las correcciones anteriores y ejecute el flujo de trabajo nuevamente:

```bash
nextflow run exhausted.nf
```

??? success "Salida del comando"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `exhausted.nf` [maniac_leavitt] DSL2 - revision: f372a56a7d

    executor >  local (3)
    [80/0779e9] PROCESS_FILES (3) | 3 of 3 ✔
    ```

Ahora debería ver las tres muestras siendo procesadas en lugar de solo una.

### 2.3. Estructura de Contenido de Canal Incorrecta

Cuando los flujos de trabajo alcanzan un cierto nivel de complejidad, puede ser un poco difícil hacer un seguimiento de las estructuras internas de cada canal, y las personas comúnmente generan desajustes entre lo que el proceso espera y lo que el canal realmente contiene. Esto es más sutil que el problema que discutimos anteriormente, donde el número de canales era incorrecto. En este caso, puede tener el número correcto de canales de entrada, pero la estructura interna de uno o más de esos canales no coincide con lo que el proceso espera.

#### Ejecutar el pipeline

```bash
nextflow run bad_channel_shape.nf
```

??? failure "Salida del comando"

    ```console
    Launching `bad_channel_shape.nf` [hopeful_pare] DSL2 - revision: ffd66071a1

    executor >  local (3)
    executor >  local (3)
    [3f/c2dcb3] PROCESS_FILES (3) [  0%] 0 of 3 ✘
    ERROR ~ Error executing process > 'PROCESS_FILES (1)'

    Caused by:
      Missing output file(s) `[sample1, file1.txt]_output.txt` expected by process `PROCESS_FILES (1)`


    Command executed:

      echo "Processing [sample1, file1.txt]" > [sample1, file1.txt]_output.txt

    Command exit status:
      0

    Command output:
      (empty)

    Work dir:
      /workspaces/training/side-quests/debugging/work/d6/1fb69d1d93300bbc9d42f1875b981e

    Tip: when you have fixed the problem you can continue the execution adding the option `-resume` to the run command line

    -- Check '.nextflow.log' file for details
    ```

#### Verificar el código

Los corchetes en el mensaje de error proporcionan la pista aquí - el proceso está tratando la tupla como un solo valor, lo cual no es lo que queremos. Examinemos `bad_channel_shape.nf`:

```groovy title="bad_channel_shape.nf" hl_lines="5 20-22" linenums="1"
#!/usr/bin/env nextflow

process PROCESS_FILES {
    input:
        val sample_name  // Espera un solo valor, obtiene una tupla

    output:
        path "${sample_name}_output.txt"

    script:
    """
    echo "Processing ${sample_name}" > ${sample_name}_output.txt
    """
}

workflow {

    // El canal emite tuplas, pero el proceso espera valores únicos
    input_ch = channel.of(
      ['sample1', 'file1.txt'],
      ['sample2', 'file2.txt'],
      ['sample3', 'file3.txt']
    )
    PROCESS_FILES(input_ch)
}
```

Puede ver que estamos generando un canal compuesto de tuplas: `['sample1', 'file1.txt']`, pero el proceso espera un solo valor, `val sample_name`. El comando ejecutado muestra que el proceso está intentando crear un archivo llamado `[sample3, file3.txt]_output.txt`, que no es la salida prevista.

#### Corregir el código

Para corregir esto, si el proceso requiere ambas entradas podríamos ajustar el proceso para aceptar una tupla:

=== "Opción 1: Aceptar tupla en el proceso"

    === "Después"

        ```groovy title="bad_channel_shape.nf" hl_lines="5"  linenums="1"
        #!/usr/bin/env nextflow

        process PROCESS_FILES {
            input:
                tuple val(sample_name), val(file_name)  // Corregido: Aceptar tupla

            output:
                path "${sample_name}_output.txt"

            script:
            """
            echo "Processing ${sample_name}" > ${sample_name}_output.txt
            """
        }

        workflow {

            // El canal emite tuplas, pero el proceso espera valores únicos
            input_ch = channel.of(
              ['sample1', 'file1.txt'],
              ['sample2', 'file2.txt'],
              ['sample3', 'file3.txt']
            )
            PROCESS_FILES(input_ch)
        }
        ```

    === "Antes"

        ```groovy title="bad_channel_shape.nf" hl_lines="5" linenums="1"
        #!/usr/bin/env nextflow

        process PROCESS_FILES {
            input:
                val sample_name  // Espera un solo valor, obtiene una tupla

            output:
                path "${sample_name}_output.txt"

            script:
            """
            echo "Processing ${sample_name}" > ${sample_name}_output.txt
            """
        }

        workflow {

            // El canal emite tuplas, pero el proceso espera valores únicos
            input_ch = channel.of(
              ['sample1', 'file1.txt'],
              ['sample2', 'file2.txt'],
              ['sample3', 'file3.txt']
            )
            PROCESS_FILES(input_ch)
        }
        ```

=== "Opción 2: Extraer primer elemento"

    === "Después"

        ```groovy title="bad_channel_shape.nf" hl_lines="9" linenums="16"
        workflow {

            // El canal emite tuplas, pero el proceso espera valores únicos
            input_ch = channel.of(
              ['sample1', 'file1.txt'],
              ['sample2', 'file2.txt'],
              ['sample3', 'file3.txt']
            )
            PROCESS_FILES(input_ch.map { it[0] })  // Corregido: Extraer primer elemento
        }
        ```

    === "Antes"

        ```groovy title="bad_channel_shape.nf" hl_lines="9" linenums="16"
        workflow {

            // El canal emite tuplas, pero el proceso espera valores únicos
            input_ch = channel.of(
              ['sample1', 'file1.txt'],
              ['sample2', 'file2.txt'],
              ['sample3', 'file3.txt']
            )
            PROCESS_FILES(input_ch)
        }
        ```

#### Ejecutar el pipeline

Elija una de las soluciones y vuelva a ejecutar el flujo de trabajo:

```bash
nextflow run bad_channel_shape.nf
```

??? success "Salida del comando"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `bad_channel_shape.nf` [clever_thompson] DSL2 - revision: 8cbcae3746

    executor >  local (3)
    [bb/80a958] PROCESS_FILES (2) | 3 of 3 ✔
    ```

### 2.4. Técnicas de Depuración de Canales

#### Uso de `.view()` para Inspección de Canales

La herramienta de depuración más poderosa para canales es el operador `.view()`. Con `.view()`, puede entender la forma de sus canales en todas las etapas para ayudar con la depuración.

#### Ejecutar el pipeline

Ejecute `bad_channel_shape_viewed.nf` para ver esto en acción:

```bash
nextflow run bad_channel_shape_viewed.nf
```

??? success "Salida del comando"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `bad_channel_shape_viewed.nf` [maniac_poisson] DSL2 - revision: b4f24dc9da

    executor >  local (3)
    [c0/db76b3] PROCESS_FILES (3) [100%] 3 of 3 ✔
    Channel content: [sample1, file1.txt]
    Channel content: [sample2, file2.txt]
    Channel content: [sample3, file3.txt]
    After mapping: sample1
    After mapping: sample2
    After mapping: sample3
    ```

#### Verificar el código

Examinemos `bad_channel_shape_viewed.nf` para ver cómo se usa `.view()`:

```groovy title="bad_channel_shape_viewed.nf" linenums="16" hl_lines="9 11"
workflow {

    // El canal emite tuplas, pero el proceso espera valores únicos
    input_ch = channel.of(
      ['sample1', 'file1.txt'],
      ['sample2', 'file2.txt'],
      ['sample3', 'file3.txt']
    )
    .view { "Channel content: $it" }  // Debug: Mostrar contenido original del canal
    .map { tuple -> tuple[0] }        // Transform: Extraer primer elemento
    .view { "After mapping: $it" }    // Debug: Mostrar contenido transformado del canal

    PROCESS_FILES(input_ch)
}
```

#### Corregir el código

Para evitar usar operaciones `.view()` excesivamente en el futuro para entender el contenido del canal, es aconsejable agregar algunos comentarios para ayudar:

```groovy title="bad_channel_shape_viewed.nf (con comentarios)" linenums="16" hl_lines="8 9"
workflow {

    // El canal emite tuplas, pero el proceso espera valores únicos
    input_ch = channel.of(
            ['sample1', 'file1.txt'],
            ['sample2', 'file2.txt'],
            ['sample3', 'file3.txt'],
        ) // [sample_name, file_name]
        .map { tuple -> tuple[0] } // sample_name

    PROCESS_FILES(input_ch)
}
```

Esto se volverá más importante a medida que sus flujos de trabajo crezcan en complejidad y la estructura del canal se vuelva más opaca.

#### Ejecutar el pipeline

```bash
nextflow run bad_channel_shape_viewed.nf
```

??? success "Salida del comando"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `bad_channel_shape_viewed.nf` [marvelous_koch] DSL2 - revision: 03e79cdbad

    executor >  local (3)
    [ff/d67cec] PROCESS_FILES (2) | 3 of 3 ✔
    Channel content: [sample1, file1.txt]
    Channel content: [sample2, file2.txt]
    Channel content: [sample3, file3.txt]
    After mapping: sample1
    After mapping: sample2
    After mapping: sample3
    ```

### Conclusión

Muchos errores de estructura de canal pueden crearse con sintaxis de Nextflow válida. Puede depurar errores de estructura de canal entendiendo el flujo de datos, usando operadores `.view()` para inspección y reconociendo patrones de mensajes de error como corchetes que indican estructuras de tupla inesperadas.

### ¿Qué sigue?

Aprenda sobre errores creados por definiciones de procesos.

---

## 3. Errores de Estructura de Proceso

La mayoría de los errores que encuentre relacionados con procesos estarán relacionados con errores que ha cometido al formar el comando, o con problemas relacionados con el software subyacente. Dicho esto, de manera similar a los problemas de canal anteriores, puede cometer errores en la definición del proceso que no califican como errores de sintaxis, pero que causarán errores en tiempo de ejecución.

### 3.1. Archivos de Salida Faltantes

Un error común al escribir procesos es hacer algo que genera un desajuste entre lo que el proceso espera y lo que se genera.

#### Ejecutar el pipeline

```bash
nextflow run missing_output.nf
```

??? failure "Salida del comando"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `missing_output.nf` [zen_stone] DSL2 - revision: 37ff61f926

    executor >  local (3)
    executor >  local (3)
    [fd/2642e9] process > PROCESS_FILES (2) [ 66%] 2 of 3, failed: 2
    ERROR ~ Error executing process > 'PROCESS_FILES (3)'

    Caused by:
      Missing output file(s) `sample3.txt` expected by process `PROCESS_FILES (3)`


    Command executed:

      echo "Processing sample3" > sample3_output.txt

    Command exit status:
      0

    Command output:
      (empty)

    Work dir:
      /workspaces/training/side-quests/debugging/work/02/9604d49fb8200a74d737c72a6c98ed

    Tip: when you have fixed the problem you can continue the execution adding the option `-resume` to the run command line

    -- Check '.nextflow.log' file for details
    ```

#### Verificar el código

El mensaje de error indica que el proceso esperaba producir un archivo de salida llamado `sample3.txt`, pero el script realmente crea `sample3_output.txt`. Examinemos la definición del proceso en `missing_output.nf`:

```groovy title="missing_output.nf" linenums="3" hl_lines="6 10"
process PROCESS_FILES {
    input:
    val sample_name

    output:
    path "${sample_name}.txt"  // Espera: sample3.txt

    script:
    """
    echo "Processing ${sample_name}" > ${sample_name}_output.txt  // Crea: sample3_output.txt
    """
}
```

Debería ver que hay un desajuste entre el nombre del archivo de salida en el bloque `output:`, y el utilizado en el script. Este desajuste hace que el proceso falle. Si encuentra este tipo de error, revise que las salidas coincidan entre su definición de proceso y su bloque de salida.

Si el problema aún no está claro, verifique el directorio de trabajo mismo para identificar los archivos de salida reales creados:

```bash
❯ ls -h work/02/9604d49fb8200a74d737c72a6c98ed
sample3_output.txt
```

Para este ejemplo, esto nos resaltaría que se está incorporando un sufijo `_output` en el nombre del archivo de salida, contrario a nuestra definición `output:`.

#### Corregir el código

Corrija el desajuste haciendo que el nombre del archivo de salida sea consistente:

=== "Después"

    ```groovy title="missing_output.nf" hl_lines="6 10" linenums="3"
    process PROCESS_FILES {
        input:
        val sample_name

        output:
        path "${sample_name}_output.txt"  // Corregido: Coincidir con la salida del script

        script:
        """
        echo "Processing ${sample_name}" > ${sample_name}_output.txt
        """
    }
    ```

=== "Antes"

    ```groovy title="missing_output.nf" hl_lines="6 10" linenums="3"
    process PROCESS_FILES {
        input:
        val sample_name

        output:
        path "${sample_name}.txt"  // Espera: sample3.txt

        script:
        """
        echo "Processing ${sample_name}" > ${sample_name}_output.txt  // Crea: sample3_output.txt
        """
    }
    ```

#### Ejecutar el pipeline

```bash
nextflow run missing_output.nf
```

??? success "Salida del comando"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `missing_output.nf` [elated_hamilton] DSL2 - revision: 961938ee2b

    executor >  local (3)
    [16/1c437c] PROCESS_FILES (3) | 3 of 3 ✔
    ```

### 3.2. Software faltante

Otra clase de errores ocurre debido a errores en el aprovisionamiento de software. `missing_software.nf` es un flujo de trabajo sintácticamente válido, pero depende de algún software externo para proporcionar el comando `cowpy` que utiliza.

#### Ejecutar el pipeline

```bash
nextflow run missing_software.nf
```

??? failure "Salida del comando"

    ```console hl_lines="12 18"
    ERROR ~ Error executing process > 'PROCESS_FILES (3)'

    Caused by:
      Process `PROCESS_FILES (3)` terminated with an error exit status (127)


    Command executed:

      cowpy sample3 > sample3_output.txt

    Command exit status:
      127

    Command output:
      (empty)

    Command error:
      .command.sh: line 2: cowpy: command not found

    Work dir:
      /workspaces/training/side-quests/debugging/work/82/42a5bfb60c9c6ee63ebdbc2d51aa6e

    Tip: you can try to figure out what's wrong by changing to the process work directory and showing the script file named `.command.sh`

    -- Check '.nextflow.log' file for details
    ```

El proceso no tiene acceso al comando que estamos especificando. A veces esto es porque un script está presente en el directorio `bin` del flujo de trabajo, pero no se ha hecho ejecutable. Otras veces es porque el software no está instalado en el contenedor o entorno donde se está ejecutando el flujo de trabajo.

#### Verificar el código

Esté atento a ese código de salida `127` - le dice exactamente el problema. Examinemos `missing_software.nf`:

```groovy title="missing_software.nf" linenums="3" hl_lines="3"
process PROCESS_FILES {

    container 'community.wave.seqera.io/library/cowpy:1.1.5--3db457ae1977a273'

    input:
    val sample_name

    output:
    path "${sample_name}_output.txt"

    script:
    """
    cowpy ${sample_name} > ${sample_name}_output.txt
    """
}
```

#### Corregir el código

Hemos sido un poco engañosos aquí, y en realidad no hay nada malo con el código. Solo necesitamos especificar la configuración necesaria para ejecutar el proceso de tal manera que tenga acceso al comando en cuestión. En este caso, el proceso tiene una definición de contenedor, por lo que todo lo que necesitamos hacer es ejecutar el flujo de trabajo con Docker habilitado.

#### Ejecutar el pipeline

Hemos configurado un perfil de Docker para usted en `nextflow.config`, por lo que puede ejecutar el flujo de trabajo con:

```bash
nextflow run missing_software.nf -profile docker
```

??? success "Salida del comando"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `missing_software.nf` [awesome_stonebraker] DSL2 - revision: 0296d12839

    executor >  local (3)
    [38/ab20d1] PROCESS_FILES (1) | 3 of 3 ✔
    ```

!!! note

    Para aprender más sobre cómo Nextflow usa contenedores, vea [Hello Nextflow](../hello_nextflow/05_hello_containers.md)

### 3.3. Mala configuración de recursos

En uso de producción, estará configurando recursos en sus procesos. Por ejemplo, `memory` define la cantidad máxima de memoria disponible para su proceso, y si el proceso excede eso, su planificador típicamente matará el proceso y devolverá un código de salida de `137`. No podemos demostrar eso aquí porque estamos usando el executor `local`, pero podemos mostrar algo similar con `time`.

#### Ejecutar el pipeline

`bad_resources.nf` tiene una configuración de proceso con un límite de tiempo poco realista de 1 milisegundo:

```bash
nextflow run bad_resources.nf -profile docker
```

??? failure "Salida del comando"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `bad_resources.nf` [disturbed_elion] DSL2 - revision: 27d2066e86

    executor >  local (3)
    [c0/ded8e1] PROCESS_FILES (3) | 0 of 3 ✘
    ERROR ~ Error executing process > 'PROCESS_FILES (2)'

    Caused by:
      Process exceeded running time limit (1ms)

    Command executed:

      cowpy sample2 > sample2_output.txt

    Command exit status:
      -

    Command output:
      (empty)

    Work dir:
      /workspaces/training/side-quests/debugging/work/53/f0a4cc56d6b3dc2a6754ff326f1349

    Container:
      community.wave.seqera.io/library/cowpy:1.1.5--3db457ae1977a273

    Tip: you can replicate the issue by changing to the process work dir and entering the command `bash .command.run`

     -- Check '.nextflow.log' file for details
    ```

#### Verificar el código

Examinemos `bad_resources.nf`:

```groovy title="bad_resources.nf" linenums="3" hl_lines="3"
process PROCESS_FILES {

    time '1 ms'  // ERROR: Límite de tiempo poco realista

    input:
    val sample_name

    output:
    path "${sample_name}_output.txt"

    script:
    """
    sleep 1  // Tarda 1 segundo, pero el límite de tiempo es 1ms
    cowpy ${sample_name} > ${sample_name}_output.txt
    """
}
```

Sabemos que el proceso tardará más de un segundo (hemos añadido un sleep ahí para asegurarnos), pero el proceso está configurado para expirar después de 1 milisegundo. ¡Alguien ha sido un poco poco realista con su configuración!

#### Corregir el código

Aumente el límite de tiempo a un valor realista:

=== "Después"

    ```groovy title="bad_resources.nf" hl_lines="3" linenums="3"
    process PROCESS_FILES {

        time '100 s'  // Corregido: Límite de tiempo realista

        input:
        val sample_name

        output:
        path "${sample_name}_output.txt"

        script:
        """
        sleep 1
        cowpy ${sample_name} > ${sample_name}_output.txt
        """
    }
    ```

=== "Antes"

    ```groovy title="bad_resources.nf" hl_lines="3" linenums="3"
    process PROCESS_FILES {

        time '1 ms'  // ERROR: Límite de tiempo poco realista

        input:
        val sample_name

        output:
        path "${sample_name}_output.txt"

        script:
        """
        sleep 1  // Tarda 1 segundo, pero el límite de tiempo es 1ms
        cowpy ${sample_name} > ${sample_name}_output.txt
        """
    }
    ```

#### Ejecutar el pipeline

```bash
nextflow run bad_resources.nf -profile docker
```

??? success "Salida del comando"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `bad_resources.nf` [friendly_mcclintock] DSL2 - revision: 381567d2c1

    executor >  local (3)
    [c2/9b4c41] PROCESS_FILES (3) | 3 of 3 ✔
    ```

Si se asegura de leer sus mensajes de error, fallos como este no deberían confundirle por mucho tiempo. Pero asegúrese de entender los requisitos de recursos de los comandos que está ejecutando para poder configurar sus directivas de recursos apropiadamente.

### 3.4. Técnicas de Depuración de Procesos

Cuando los procesos fallan o se comportan inesperadamente, necesita técnicas sistemáticas para investigar qué salió mal. El directorio de trabajo contiene toda la información que necesita para depurar la ejecución del proceso.

#### Uso de la Inspección del Directorio de Trabajo

La herramienta de depuración más poderosa para procesos es examinar el directorio de trabajo. Cuando un proceso falla, Nextflow crea un directorio de trabajo para esa ejecución específica del proceso que contiene todos los archivos necesarios para entender qué pasó.

#### Ejecutar el pipeline

Usemos el ejemplo `missing_output.nf` de antes para demostrar la inspección del directorio de trabajo (regenere un desajuste de nombres de salida si es necesario):

```bash
nextflow run missing_output.nf
```

??? failure "Salida del comando"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `missing_output.nf` [irreverent_payne] DSL2 - revision: 3d5117f7e2

    executor >  local (3)
    [5d/d544a4] PROCESS_FILES (2) | 0 of 3 ✘
    ERROR ~ Error executing process > 'PROCESS_FILES (1)'

    Caused by:
      Missing output file(s) `sample1.txt` expected by process `PROCESS_FILES (1)`

    Command executed:

      echo "Processing sample1" > sample1_output.txt

    Command exit status:
      0

    Command output:
      (empty)

    Work dir:
      /workspaces/training/side-quests/debugging/work/1e/2011154d0b0f001cd383d7364b5244

    Tip: you can replicate the issue by changing to the process work dir and entering the command `bash .command.run`

     -- Check '.nextflow.log' file for details
    ```

#### Verificar el directorio de trabajo

Cuando obtiene este error, el directorio de trabajo contiene toda la información de depuración. Encuentre la ruta del directorio de trabajo del mensaje de error y examine su contenido:

```bash
# Encuentre el directorio de trabajo del mensaje de error
ls work/02/9604d49fb8200a74d737c72a6c98ed/
```

Luego puede examinar los archivos clave:

##### Verificar el Script de Comando

El archivo `.command.sh` muestra exactamente qué comando se ejecutó:

```bash
# Ver el comando ejecutado
cat work/02/9604d49fb8200a74d737c72a6c98ed/.command.sh
```

Esto revela:

- **Sustitución de variables**: Si las variables de Nextflow fueron expandidas correctamente
- **Rutas de archivos**: Si los archivos de entrada fueron localizados correctamente
- **Estructura del comando**: Si la sintaxis del script es correcta

Problemas comunes a buscar:

- **Comillas faltantes**: Variables que contienen espacios necesitan comillas adecuadas
- **Rutas de archivos incorrectas**: Archivos de entrada que no existen o están en ubicaciones incorrectas
- **Nombres de variables incorrectos**: Errores tipográficos en referencias de variables
- **Configuración de entorno faltante**: Comandos que dependen de entornos específicos

##### Verificar la Salida de Error

El archivo `.command.err` contiene los mensajes de error reales:

```bash
# Ver salida de error
cat work/02/9604d49fb8200a74d737c72a6c98ed/.command.err
```

Este archivo mostrará:

- **Códigos de salida**: 127 (comando no encontrado), 137 (terminado), etc.
- **Errores de permisos**: Problemas de acceso a archivos
- **Errores de software**: Mensajes de error específicos de la aplicación
- **Errores de recursos**: Límite de memoria/tiempo excedido

##### Verificar la Salida Estándar

El archivo `.command.out` muestra lo que produjo su comando:

```bash
# Ver salida estándar
cat work/02/9604d49fb8200a74d737c72a6c98ed/.command.out
```

Esto ayuda a verificar:

- **Salida esperada**: Si el comando produjo los resultados correctos
- **Ejecución parcial**: Si el comando comenzó pero falló a mitad de camino
- **Información de depuración**: Cualquier salida de diagnóstico de su script

##### Verificar el Código de Salida

El archivo `.exitcode` contiene el código de salida del proceso:

```bash
# Ver código de salida
cat work/*/*/.exitcode
```

Códigos de salida comunes y sus significados:

- **Código de salida 127**: Comando no encontrado - verifique la instalación del software
- **Código de salida 137**: Proceso terminado - verifique los límites de memoria/tiempo

##### Verificar la Existencia de Archivos

Cuando los procesos fallan debido a archivos de salida faltantes, verifique qué archivos fueron realmente creados:

```bash
# Listar todos los archivos en el directorio de trabajo
ls -la work/02/9604d49fb8200a74d737c72a6c98ed/
```

Esto ayuda a identificar:

- **Desajustes de nombres de archivo**: Archivos de salida con nombres diferentes a los esperados
- **Problemas de permisos**: Archivos que no pudieron ser creados
- **Problemas de ruta**: Archivos creados en directorios incorrectos

En nuestro ejemplo anterior, esto nos confirmó que mientras nuestro esperado `sample3.txt` no estaba presente, `sample3_output.txt` sí lo estaba:

```bash
❯ ls -h work/02/9604d49fb8200a74d737c72a6c98ed
sample3_output.txt
```

### Conclusión

La depuración de procesos requiere examinar los directorios de trabajo para entender qué salió mal. Los archivos clave incluyen `.command.sh` (el script ejecutado), `.command.err` (mensajes de error) y `.command.out` (salida estándar). Los códigos de salida como 127 (comando no encontrado) y 137 (proceso terminado) proporcionan pistas diagnósticas inmediatas sobre el tipo de fallo.

### ¿Qué sigue?

Aprenda sobre las herramientas de depuración integradas de Nextflow y los enfoques sistemáticos para la resolución de problemas.

---

## 4. Herramientas de Depuración Integradas y Técnicas Avanzadas

Nextflow proporciona varias herramientas integradas poderosas para depurar y analizar la ejecución de flujos de trabajo. Estas herramientas le ayudan a entender qué salió mal, dónde salió mal y cómo solucionarlo eficientemente.

### 4.1. Salida de Proceso en Tiempo Real

A veces necesita ver qué está pasando dentro de los procesos en ejecución. Puede habilitar la salida de proceso en tiempo real, que le muestra exactamente qué está haciendo cada tarea mientras se ejecuta.

#### Ejecutar el pipeline

`bad_channel_shape_viewed.nf` de nuestros ejemplos anteriores imprimió el contenido del canal usando `.view()`, pero también podemos usar la directiva `debug` para mostrar variables desde dentro del proceso mismo, lo cual demostramos en `bad_channel_shape_viewed_debug.nf`. Ejecute el flujo de trabajo:

```bash
nextflow run bad_channel_shape_viewed_debug.nf
```

??? success "Salida del comando"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `bad_channel_shape_viewed_debug.nf` [agitated_crick] DSL2 - revision: ea3676d9ec

    executor >  local (3)
    [c6/2dac51] process > PROCESS_FILES (3) [100%] 3 of 3 ✔
    Channel content: [sample1, file1.txt]
    Channel content: [sample2, file2.txt]
    Channel content: [sample3, file3.txt]
    After mapping: sample1
    After mapping: sample2
    After mapping: sample3
    Sample name inside process is sample2

    Sample name inside process is sample1

    Sample name inside process is sample3
    ```

#### Verificar el código

Examinemos `bad_channel_shape_viewed_debug.nf` para ver cómo funciona la directiva `debug`:

```groovy title="bad_channel_shape_viewed_debug.nf" linenums="3" hl_lines="2"
process PROCESS_FILES {
    debug true  // Enable real-time output

    input:
    val sample_name

    output:
    path "${sample_name}_output.txt"

    script:
    """
    echo "Sample name inside process is ${sample_name}"
    echo "Processing ${sample_name}" > ${sample_name}_output.txt
    """
}
```

La directiva `debug` puede ser una forma rápida y conveniente de entender el entorno de un proceso.

### 4.2. Modo Preview

A veces quiere detectar problemas antes de que se ejecute cualquier proceso. Nextflow proporciona un flag para este tipo de depuración proactiva: `-preview`.

#### Ejecutar el pipeline

El modo preview le permite probar la lógica del flujo de trabajo sin ejecutar comandos. Esto puede ser muy útil para verificar rápidamente la estructura de su flujo de trabajo y asegurar que los procesos estén conectados correctamente sin ejecutar ningún comando real.

!!! note

    Si corrigió `bad_syntax.nf` anteriormente, reintroduzca el error de sintaxis eliminando la llave de cierre después del bloque script antes de ejecutar este comando.

Ejecute este comando:

```bash
nextflow run bad_syntax.nf -preview
```

??? failure "Salida del comando"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `bad_syntax.nf` [magical_mercator] DSL2 - revision: 550b9a8873

    Error bad_syntax.nf:24:1: Unexpected input: '<EOF>'

    ERROR ~ Script compilation failed

     -- Check '.nextflow.log' file for details
    ```

El modo preview es particularmente útil para detectar errores de sintaxis tempranamente sin ejecutar ningún proceso. Valida la estructura del flujo de trabajo y las conexiones de procesos antes de la ejecución.

### 4.3. Ejecución Stub para Pruebas de Lógica

A veces los errores son difíciles de depurar porque los comandos tardan demasiado, requieren software especial o fallan por razones complejas. La ejecución stub le permite probar la lógica del flujo de trabajo sin ejecutar los comandos reales.

#### Ejecutar el pipeline

Cuando está desarrollando un proceso de Nextflow, puede usar la directiva `stub` para definir comandos 'ficticios' que generan salidas de la forma correcta sin ejecutar el comando real. Este enfoque es particularmente valioso cuando quiere verificar que la lógica de su flujo de trabajo es correcta antes de lidiar con las complejidades del software real.

Por ejemplo, ¿recuerda nuestro `missing_software.nf` de antes? ¿El que tenía software faltante que impedía que el flujo de trabajo se ejecutara hasta que añadimos `-profile docker`? `missing_software_with_stub.nf` es un flujo de trabajo muy similar. Si lo ejecutamos de la misma manera, generaremos el mismo error:

```bash
nextflow run missing_software_with_stub.nf
```

??? failure "Salida del comando"

    ```console hl_lines="12 18"
    ERROR ~ Error executing process > 'PROCESS_FILES (3)'

    Caused by:
      Process `PROCESS_FILES (3)` terminated with an error exit status (127)


    Command executed:

      cowpy sample3 > sample3_output.txt

    Command exit status:
      127

    Command output:
      (empty)

    Command error:
      .command.sh: line 2: cowpy: command not found

    Work dir:
      /workspaces/training/side-quests/debugging/work/82/42a5bfb60c9c6ee63ebdbc2d51aa6e

    Tip: you can try to figure out what's wrong by changing to the process work directory and showing the script file named `.command.sh`

    -- Check '.nextflow.log' file for details
    ```

Sin embargo, este flujo de trabajo no producirá errores si lo ejecutamos con `-stub-run`, incluso sin el perfil `docker`:

```bash
nextflow run missing_software_with_stub.nf -stub-run
```

??? success "Salida del comando"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `missing_software_with_stub.nf` [astonishing_shockley] DSL2 - revision: f1f4f05d7d

    executor >  local (3)
    [b5/2517a3] PROCESS_FILES (3) | 3 of 3 ✔
    ```

#### Verificar el código

Examinemos `missing_software_with_stub.nf`:

```groovy title="missing_software.nf (with stub)" hl_lines="16-19" linenums="3"
process PROCESS_FILES {

    container 'community.wave.seqera.io/library/cowpy:1.1.5--3db457ae1977a273'

    input:
    val sample_name

    output:
    path "${sample_name}_output.txt"

    script:
    """
    cowpy ${sample_name} > ${sample_name}_output.txt
    """

    stub:
    """
    touch ${sample_name}_output.txt
    """
}
```

Respecto a `missing_software.nf`, este proceso tiene una directiva `stub:` que especifica un comando a usar en lugar del especificado en `script:`, en el caso de que Nextflow se ejecute en modo stub.

El comando `touch` que estamos usando aquí no depende de ningún software o entradas apropiadas, y se ejecutará en todas las situaciones, permitiéndonos depurar la lógica del flujo de trabajo sin preocuparnos por los internos del proceso.

**La ejecución stub ayuda a depurar:**

- Estructura de canales y flujo de datos
- Conexiones y dependencias de procesos
- Propagación de parámetros
- Lógica del flujo de trabajo sin dependencias de software

### 4.4. Enfoque Sistemático de Depuración

Ahora que ha aprendido técnicas de depuración individuales - desde archivos de traza y directorios de trabajo hasta modo preview, ejecución stub y monitoreo de recursos - unámoslas en una metodología sistemática. Tener un enfoque estructurado evita que se sienta abrumado por errores complejos y asegura que no pierda pistas importantes.

Esta metodología combina todas las herramientas que hemos cubierto en un flujo de trabajo eficiente:

**Método de Depuración en Cuatro Fases:**

**Fase 1: Resolución de Errores de Sintaxis (5 minutos)**

1. Verifique subrayados rojos en VSCode o su IDE
2. Ejecute `nextflow run workflow.nf -preview` para identificar problemas de sintaxis
3. Corrija todos los errores de sintaxis (llaves faltantes, comas finales, etc.)
4. Asegúrese de que el flujo de trabajo se analice exitosamente antes de continuar

**Fase 2: Evaluación Rápida (5 minutos)**

1. Lea los mensajes de error de tiempo de ejecución cuidadosamente
2. Verifique si es un error de tiempo de ejecución, lógica o recursos
3. Use el modo preview para probar la lógica básica del flujo de trabajo

**Fase 3: Investigación Detallada (15-30 minutos)**

1. Encuentre el directorio de trabajo de la tarea fallida
2. Examine los archivos de registro
3. Añada operadores `.view()` para inspeccionar canales
4. Use `-stub-run` para probar la lógica del flujo de trabajo sin ejecución

**Fase 4: Corregir y Validar (15 minutos)**

1. Haga correcciones mínimas y dirigidas
2. Pruebe con resume: `nextflow run workflow.nf -resume`
3. Verifique la ejecución completa del flujo de trabajo

!!! tip "Uso de Resume para Depuración Eficiente"

    Una vez que ha identificado un problema, necesita una forma eficiente de probar sus correcciones sin perder tiempo re-ejecutando partes exitosas de su flujo de trabajo. La funcionalidad `-resume` de Nextflow es invaluable para la depuración.

    Habrá encontrado `-resume` si ha trabajado con [Hello Nextflow](../hello_nextflow/), y es importante que lo use bien al depurar para ahorrarse esperar mientras los procesos antes de su proceso problemático se ejecutan.

    **Estrategia de depuración con resume:**

    1. Ejecute el flujo de trabajo hasta el fallo
    2. Examine el directorio de trabajo de la tarea fallida
    3. Corrija el problema específico
    4. Reanude para probar solo la corrección
    5. Repita hasta que el flujo de trabajo se complete

#### Perfil de Configuración de Depuración

Para hacer este enfoque sistemático aún más eficiente, puede crear una configuración de depuración dedicada que habilite automáticamente todas las herramientas que necesita:

```groovy title="nextflow.config (debug profile)" linenums="1"
profiles {
    debug {
        process {
            debug = true
            cleanup = false

            // Conservative resources for debugging
            maxForks = 1
            memory = '2.GB'
            cpus = 1
        }
    }
}
```

Luego puede ejecutar el pipeline con este perfil habilitado:

```bash
nextflow run workflow.nf -profile debug
```

Este perfil habilita la salida en tiempo real, preserva los directorios de trabajo y limita la paralelización para una depuración más fácil.

### 4.5. Ejercicio Práctico de Depuración

Ahora es momento de poner en práctica el enfoque sistemático de depuración. El flujo de trabajo `buggy_workflow.nf` contiene varios errores comunes que representan los tipos de problemas que encontrará en el desarrollo del mundo real.

!!! exercise

    Use el enfoque sistemático de depuración para identificar y corregir todos los errores en `buggy_workflow.nf`. Este flujo de trabajo intenta procesar datos de muestra de un archivo CSV pero contiene múltiples errores intencionales que representan escenarios comunes de depuración.

    Comience ejecutando el flujo de trabajo para ver el primer error:

    ```bash
    nextflow run buggy_workflow.nf
    ```

    ??? failure "Salida del comando"

        ```console
        N E X T F L O W   ~  version 25.10.2

        Launching `buggy_workflow.nf` [wise_ramanujan] DSL2 - revision: d51a8e83fd

        ERROR ~ Range [11, 12) out of bounds for length 11

         -- Check '.nextflow.log' file for details
        ```

        Este error críptico indica un problema de análisis alrededor de las líneas 11-12 en el bloque `params{}`. El analizador v2 detecta problemas estructurales tempranamente.

    Aplique el método de depuración en cuatro fases que ha aprendido:

    **Fase 1: Resolución de Errores de Sintaxis**
    - Verifique subrayados rojos en VSCode o su IDE
    - Ejecute `nextflow run workflow.nf -preview` para identificar problemas de sintaxis
    - Corrija todos los errores de sintaxis (llaves faltantes, comas finales, etc.)
    - Asegúrese de que el flujo de trabajo se analice exitosamente antes de continuar

    **Fase 2: Evaluación Rápida**
    - Lea los mensajes de error de tiempo de ejecución cuidadosamente
    - Identifique si los errores son de tiempo de ejecución, lógica o recursos
    - Use el modo `-preview` para probar la lógica básica del flujo de trabajo

    **Fase 3: Investigación Detallada**
    - Examine los directorios de trabajo de las tareas fallidas
    - Añada operadores `.view()` para inspeccionar canales
    - Verifique los archivos de registro en los directorios de trabajo
    - Use `-stub-run` para probar la lógica del flujo de trabajo sin ejecución

    **Fase 4: Corregir y Validar**
    - Haga correcciones dirigidas
    - Use `-resume` para probar las correcciones eficientemente
    - Verifique la ejecución completa del flujo de trabajo

    **Herramientas de Depuración a Su Disposición:**
    ```bash
    # Modo preview para verificación de sintaxis
    nextflow run buggy_workflow.nf -preview

    # Perfil debug para salida detallada
    nextflow run buggy_workflow.nf -profile debug

    # Ejecución stub para pruebas de lógica
    nextflow run buggy_workflow.nf -stub-run

    # Resume después de correcciones
    nextflow run buggy_workflow.nf -resume
    ```

    ??? solution
        El `buggy_workflow.nf` contiene 9 o 10 errores distintos (dependiendo de cómo cuente) que cubren todas las categorías principales de depuración. Aquí hay un desglose sistemático de cada error y cómo corregirlo

        Comencemos con esos errores de sintaxis:

        **Error 1: Error de Sintaxis - Coma Final**
        ```groovy linenums="21"
        output:
            path "${sample_id}_result.txt",  // ERROR: Coma final
        ```
        **Corrección:** Elimine la coma final
        ```groovy linenums="21"
        output:
            path "${sample_id}_result.txt"
        ```

        **Error 2: Error de Sintaxis - Llave de Cierre Faltante**
        ```groovy linenums="24"
        script:
        """
        echo "Processing: ${sample}"
        cat ${input_file} > ${sample}_result.txt
        """
        // ERROR: Falta la llave de cierre para el proceso processFiles
        ```
        **Corrección:** Añada la llave de cierre faltante
        ```groovy linenums="29"
        """
        echo "Processing: ${sample_id}"
        cat ${input_file} > ${sample_id}_result.txt
        """
        }  // Añadir llave de cierre faltante
        ```

        **Error 3: Error de Nombre de Variable**
        ```groovy linenums="26"
        echo "Processing: ${sample}"     // ERROR: debería ser sample_id
        cat ${input_file} > ${sample}_result.txt  // ERROR: debería ser sample_id
        ```
        **Corrección:** Use el nombre correcto de la variable de entrada
        ```groovy linenums="26"
        echo "Processing: ${sample_id}"
        cat ${input_file} > ${sample_id}_result.txt
        ```

        **Error 4: Error de Variable No Definida**
        ```groovy linenums="87"
        heavy_ch = heavyProcess(sample_ids)  // ERROR: sample_ids no definido
        ```
        **Corrección:** Use el canal correcto y extraiga los IDs de muestra
        ```groovy linenums="87"
        heavy_ch = heavyProcess(input_ch)
        ```

        En este punto el flujo de trabajo se ejecutará, pero seguiremos obteniendo errores (ej. `Path value cannot be null` en `processFiles`), causados por una estructura de canal incorrecta.

        **Error 5: Error de Estructura de Canal - Salida Map Incorrecta**
        ```groovy linenums="83"
        .map { row -> row.sample_id }  // ERROR: processFiles espera tupla
        ```
        **Corrección:** Devuelva la estructura de tupla que processFiles espera
        ```groovy linenums="83"
        .map { row -> [row.sample_id, file(row.fastq_path)] }
        ```

        Pero esto romperá nuestra solución para ejecutar `heavyProcess()` arriba, así que necesitaremos usar un map para pasar solo los IDs de muestra a ese proceso:

        **Error 6: Estructura de canal incorrecta para heavyProcess**
        ```groovy linenums="87"
        heavy_ch = heavyProcess(input_ch)  // ERROR: input_ch ahora tiene 2 elementos por emisión - heavyProcess solo necesita 1 (el primero)
        ```
        **Corrección:** Use el canal correcto y extraiga los IDs de muestra
        ```groovy linenums="87"
        heavy_ch = heavyProcess(input_ch.map{it[0]})
        ```

        Ahora avanzamos un poco más pero recibimos un error sobre `No such variable: i`, porque no escapamos una variable de Bash.

        **Error 7: Error de Escapado de Variable Bash**
        ```groovy linenums="48"
        echo "Heavy computation $i for ${sample_id}"  // ERROR: $i no escapado
        ```
        **Corrección:** Escape la variable bash
        ```groovy linenums="48"
        echo "Heavy computation \${i} for ${sample_id}"
        ```

        Ahora obtenemos `Process exceeded running time limit (1ms)`, así que corregimos el límite de tiempo de ejecución para el proceso relevante:

        **Error 8: Error de Configuración de Recursos**
        ```groovy linenums="36"
        time '1 ms'  // ERROR: Límite de tiempo poco realista
        ```
        **Corrección:** Aumente a un límite de tiempo realista
        ```groovy linenums="36"
        time '100 s'
        ```

        A continuación tenemos un error `Missing output file(s)` para resolver:

        **Error 9: Desajuste de Nombre de Archivo de Salida**
        ```groovy linenums="49"
        done > ${sample_id}.txt  // ERROR: Nombre de archivo incorrecto, debería coincidir con la declaración de salida
        ```
        **Corrección:** Coincida con la declaración de salida
        ```groovy linenums="49"
        done > ${sample_id}_heavy.txt
        ```

        Los primeros dos procesos se ejecutaron, pero no el tercero.

        **Error 10: Desajuste de Nombre de Archivo de Salida**
        ```groovy linenums="88"
        file_ch = channel.fromPath("*.txt") // Error: intentando tomar entrada del pwd en vez de un proceso
        handleFiles(file_ch)
        ```
        **Corrección:** Tome la salida del proceso anterior
        ```groovy linenums="88"
        handleFiles(heavyProcess.out)
        ```

        Con eso, todo el flujo de trabajo debería ejecutarse.

        **Flujo de Trabajo Corregido Completo:**
        ```groovy linenums="1"
        #!/usr/bin/env nextflow

        /*
        * Buggy workflow for debugging exercises
        * This workflow contains several intentional bugs for learning purposes
        */

        params{
            // Parameters with missing validation
            input: Path = 'data/sample_data.csv'
            output: String = 'results'
        }

        /*
        * Process with input/output mismatch
        */
        process processFiles {
            publishDir "${params.output}/processed", mode: 'copy'

            input:
                tuple val(sample_id), path(input_file)

            output:
                path "${sample_id}_result.txt"

            script:
            """
            echo "Processing: ${sample_id}"
            cat ${input_file} > ${sample_id}_result.txt
            """
        }

        /*
        * Process with resource issues
        */
        process heavyProcess {
            publishDir "${params.output}/heavy", mode: 'copy'

            time '100 s'

            input:
                val sample_id

            output:
                path "${sample_id}_heavy.txt"

            script:
            """
            # Simulate heavy computation
            for i in {1..1000000}; do
                echo "Heavy computation \$i for ${sample_id}"
            done > ${sample_id}_heavy.txt
            """
        }

        /*
        * Process with file handling issues
        */
        process handleFiles {
            publishDir "${params.output}/files", mode: 'copy'

            input:
                path input_file

            output:
                path "processed_${input_file}"

            script:
            """
            if [ -f "${input_file}" ]; then
                cp ${input_file} processed_${input_file}
            fi
            """
        }

        /*
        * Main workflow with channel issues
        */
        workflow {

            // Channel with incorrect usage
            input_ch = channel
                .fromPath(params.input)
                .splitCsv(header: true)
                .map { row -> [row.sample_id, file(row.fastq_path)] }

            processed_ch = processFiles(input_ch)

            heavy_ch = heavyProcess(input_ch.map{it[0]})

            handleFiles(heavyProcess.out)
        }
        ```

**Categorías de Error Cubiertas:**

- **Errores de sintaxis**: Llaves faltantes, comas finales, variables no definidas
- **Errores de estructura de canal**: Formas de datos incorrectas, canales no definidos
- **Errores de proceso**: Desajustes de archivos de salida, escapado de variables
- **Errores de recursos**: Límites de tiempo poco realistas

**Lecciones Clave de Depuración:**

1. **Lea los mensajes de error cuidadosamente** - a menudo apuntan directamente al problema
2. **Use enfoques sistemáticos** - corrija un error a la vez y pruebe con `-resume`
3. **Entienda el flujo de datos** - los errores de estructura de canal son a menudo los más sutiles
4. **Verifique los directorios de trabajo** - cuando los procesos fallan, los registros le dicen exactamente qué salió mal

---

## Resumen

En esta misión secundaria, ha aprendido un conjunto de técnicas sistemáticas para depurar flujos de trabajo de Nextflow.
Aplicar estas técnicas en su propio trabajo le permitirá pasar menos tiempo luchando con su computadora, resolver problemas más rápido y protegerse de problemas futuros.

### Patrones clave

**1. Cómo identificar y corregir errores de sintaxis**:

- Interpretación de mensajes de error de Nextflow y localización de problemas
- Errores de sintaxis comunes: llaves faltantes, palabras clave incorrectas, variables no definidas
- Distinción entre variables de Nextflow (Groovy) y Bash
- Uso de las características de la extensión VS Code para detección temprana de errores

```groovy
// Llave faltante - busque subrayados rojos en el IDE
process FOO {
    script:
    """
    echo "hello"
    """
// } <-- ¡faltante!

// Palabra clave incorrecta
inputs:  // Debería ser 'input:'

// Variable no definida - escape con barra invertida para variables Bash
echo "${undefined_var}"      // Variable Nextflow (error si no está definida)
echo "\${bash_var}"          // Variable Bash (escapada)
```

**2. Cómo depurar problemas de estructura de canal**:

- Comprensión de la cardinalidad de canales y problemas de agotamiento
- Depuración de desajustes de estructura de contenido de canal
- Uso de operadores `.view()` para inspección de canales
- Reconocimiento de patrones de error como corchetes en la salida

```groovy
// Inspeccionar contenido del canal
my_channel.view { "Content: $it" }

// Convertir canal de cola a canal de valor (previene agotamiento)
reference_ch = channel.value('ref.fa')
// o
reference_ch = channel.of('ref.fa').first()
```

**3. Cómo solucionar problemas de ejecución de procesos**:

- Diagnóstico de errores de archivos de salida faltantes
- Comprensión de códigos de salida (127 para software faltante, 137 para problemas de memoria)
- Investigación de directorios de trabajo y archivos de comando
- Configuración apropiada de recursos

```bash
# Verificar qué se ejecutó realmente
cat work/ab/cdef12/.command.sh

# Verificar salida de error
cat work/ab/cdef12/.command.err

# Código de salida 127 = comando no encontrado
# Código de salida 137 = terminado (límite de memoria/tiempo)
```

**4. Cómo usar las herramientas de depuración integradas de Nextflow**:

- Aprovechamiento del modo preview y depuración en tiempo real
- Implementación de ejecución stub para pruebas de lógica
- Aplicación de resume para ciclos de depuración eficientes
- Seguimiento de una metodología sistemática de depuración en cuatro fases

!!! tip "Referencia Rápida de Depuración"

    **¿Errores de sintaxis?** → Verifique advertencias de VSCode, ejecute `nextflow run workflow.nf -preview`

    **¿Problemas de canal?** → Use `.view()` para inspeccionar contenido: `my_channel.view()`

    **¿Fallos de proceso?** → Verifique archivos del directorio de trabajo:

    - `.command.sh` - el script ejecutado
    - `.command.err` - mensajes de error
    - `.exitcode` - estado de salida (127 = comando no encontrado, 137 = terminado)

    **¿Comportamiento misterioso?** → Ejecute con `-stub-run` para probar la lógica del flujo de trabajo

    **¿Hizo correcciones?** → Use `-resume` para ahorrar tiempo probando: `nextflow run workflow.nf -resume`

---

### Recursos adicionales

- [Guía de resolución de problemas de Nextflow](https://www.nextflow.io/docs/latest/troubleshooting.html): Documentación oficial de resolución de problemas
- [Comprender los canales de Nextflow](https://www.nextflow.io/docs/latest/channel.html): Inmersión profunda en tipos y comportamiento de canales
- [Referencia de directivas de proceso](https://www.nextflow.io/docs/latest/process.html#directives): Todas las opciones de configuración de proceso disponibles
- [nf-test](https://www.nf-test.com/): Framework de pruebas para pipelines de Nextflow
- [Comunidad Slack de Nextflow](https://www.nextflow.io/slack-invite.html): Obtenga ayuda de la comunidad

Para flujos de trabajo de producción, considere:

- Configurar [Seqera Platform](https://seqera.io/platform/) para monitoreo y depuración a escala
- Usar [Wave containers](https://seqera.io/wave/) para entornos de software reproducibles

**Recuerde:** La depuración efectiva es una habilidad que mejora con la práctica. La metodología sistemática y el conjunto de herramientas integral que ha adquirido aquí le servirán bien a lo largo de su viaje de desarrollo con Nextflow.

---

## ¿Qué sigue?

Regrese al [menú de Misiones Secundarias](./index.md) o haga clic en el botón en la parte inferior derecha de la página para avanzar al siguiente tema de la lista.
