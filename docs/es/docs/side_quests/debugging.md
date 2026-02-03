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

En uso de producción, estará configurando recursos en sus procesos. Por ejemplo, `memory` define la cantidad máxima de memoria disponible para su proceso, y si el proceso excede eso, su planificador típicamente matará el proceso y devolverá un código de salida de `137`. No podemos demostrar eso aquí porque estamos usando el executor `local`, pero podemos mostrar
