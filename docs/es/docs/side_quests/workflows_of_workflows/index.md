# Workflows de Workflows

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Traducción asistida por IA - [más información y sugerencias](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

Cuando desarrollas un pipeline, con frecuencia te encuentras creando secuencias similares de procesos para diferentes tipos de datos o pasos de análisis.
Podrías terminar copiando y pegando estas secuencias de procesos, lo que genera código duplicado difícil de mantener; o podrías crear un workflow masivo que sea difícil de entender y modificar.

Una de las características más poderosas de Nextflow es su capacidad para componer pipelines complejos a partir de módulos de workflow más pequeños y reutilizables.
Este enfoque modular hace que los pipelines sean más fáciles de desarrollar, probar y mantener.

### Objetivos de aprendizaje

En esta misión secundaria, exploraremos cómo desarrollar módulos de workflow que puedan probarse y usarse de forma independiente, componer esos módulos en un pipeline más grande y gestionar el flujo de datos entre módulos.

Al finalizar esta misión secundaria, podrás:

- Dividir pipelines complejos en unidades lógicas y reutilizables
- Probar cada módulo de workflow de forma independiente
- Combinar workflows para crear nuevos pipelines
- Compartir módulos de workflow comunes entre diferentes pipelines
- Hacer que tu código sea más fácil de mantener y entender

Estas habilidades te ayudarán a construir pipelines complejos manteniendo una estructura de código limpia y fácil de mantener.

### Requisitos previos

Antes de comenzar esta misión secundaria deberías:

- Haber completado el tutorial [Hello Nextflow](../../hello_nextflow/index.md) o un curso equivalente para principiantes.
- Sentirte cómodo/a usando los conceptos y mecanismos básicos de Nextflow (procesos, canales, operadores, módulos)

---

## 0. Primeros pasos

#### Abre el codespace de capacitación

Si aún no lo has hecho, asegúrate de abrir el entorno de capacitación como se describe en la [Configuración del entorno](../../envsetup/index.md).

[![Open in GitHub Codespaces](https://github.com/codespaces/badge.svg)](https://codespaces.new/nextflow-io/training?quickstart=1&ref=master)

#### Muévete al directorio del proyecto

Muévete al directorio donde se encuentran los archivos de este tutorial.

```bash
cd side-quests/workflows_of_workflows
```

Puedes configurar VSCode para que se enfoque en este directorio:

```bash
code .
```

El editor se abre con el directorio del proyecto en foco.

#### Revisa los materiales

Encontrarás un directorio `modules` con definiciones de procesos, un directorio `workflows` con dos scripts de workflow preescritos, y un archivo `main.nf` que irás actualizando progresivamente:

```console title="Directory contents"
├── main.nf
├── workflows/
│   ├── greeting.nf              # Standalone greeting workflow (to be made composable)
│   └── transform.nf             # Standalone transform workflow (to be made composable)
└── modules/
    ├── say_hello.nf             # Creates a greeting (from Hello Nextflow)
    ├── say_hello_upper.nf       # Converts to uppercase (from Hello Nextflow)
    ├── timestamp_greeting.nf    # Adds timestamps to greetings
    ├── validate_name.nf         # Validates input names
    └── reverse_text.nf          # Reverses text content
```

El directorio `modules/` contiene las definiciones individuales de procesos, y el directorio `workflows/` contiene los dos scripts de workflow preescritos con los que trabajarás en esta misión secundaria.

#### Revisa la tarea

Tu desafío es ensamblar estos módulos en dos workflows separados que luego compondremos en un workflow principal:

- Un `GREETING_WORKFLOW` que valida nombres, crea saludos y agrega marcas de tiempo
- Un `TRANSFORM_WORKFLOW` que convierte texto a mayúsculas y lo invierte

El pipeline terminado encadena ambos workflows en un único flujo de datos:

1. **Validar**: verificar que cada nombre esté bien formado
2. **Saludar**: generar un saludo para cada nombre válido
3. **Marca de tiempo**: registrar cuándo se creó cada saludo
4. **Mayúsculas**: convertir el saludo con marca de tiempo a mayúsculas
5. **Invertir**: invertir el texto en mayúsculas

Los primeros tres pasos pertenecen a `GREETING_WORKFLOW`, y los últimos dos pertenecen a `TRANSFORM_WORKFLOW`.
Construirlos como workflows separados y componibles te permite desarrollar y probar cada etapa de forma independiente antes de conectarlos entre sí.

#### Lista de verificación de preparación

¿Crees que estás listo/a para comenzar?

- [ ] Entiendo el objetivo de este curso y sus requisitos previos
- [ ] Mi codespace está en funcionamiento
- [ ] He configurado mi directorio de trabajo correctamente
- [ ] Entiendo la tarea

Si puedes marcar todas las casillas, estás listo/a para continuar.

---

## 1. Agregar el greeting workflow al pipeline

El greeting workflow valida nombres y genera saludos con marcas de tiempo.

### 1.1. Revisar y ejecutar el greeting workflow

Abre `workflows/greeting.nf` y examina el código:

```groovy title="workflows/greeting.nf" linenums="1"
#!/usr/bin/env nextflow

include { VALIDATE_NAME } from '../modules/validate_name'
include { SAY_HELLO } from '../modules/say_hello'
include { TIMESTAMP_GREETING } from '../modules/timestamp_greeting'

workflow {
    main:
    names_ch = channel.of('Alice', 'Bob', 'Charlie')

    // Encadena los procesos: validar -> crear saludo -> agregar marca de tiempo
    validated_ch = VALIDATE_NAME(names_ch)
    greetings_ch = SAY_HELLO(validated_ch)
    timestamped_ch = TIMESTAMP_GREETING(greetings_ch)

    publish:
    greetings = greetings_ch
    timestamped = timestamped_ch
}

output {
    greetings {
    }
    timestamped {
    }
}
```

Este es un workflow completo y autónomo con la misma estructura que viste en el tutorial 'Hello Nextflow'.
Tiene los nombres de entrada codificados directamente, encadena tres procesos y publica dos salidas.

Ejecútalo para verificar que todo funciona:

```bash
nextflow run workflows/greeting.nf
```

??? success "Salida del comando"

    ```console
    N E X T F L O W   ~  version 26.04.4
    Launching `workflows/greeting.nf` [sharp_kirch] revision: 04c0a33f79
    executor >  local (9)
    [c5/9df3c8] VAL…TE_NAME (validating Alice) | 3 of 3 ✔
    [e0/243874] SAY_HELLO (greeting Alice)     | 3 of 3 ✔
    [c5/e40c8d] TIM…ing timestamp to greeting) | 3 of 3 ✔

    Outputs:

      /workspaces/training/side-quests/workflows_of_workflows/results

      greetings:
        - Charlie-output.txt
        - Bob-output.txt
        - Alice-output.txt

      timestamped:
        - timestamped_Bob-output.txt
        - timestamped_Charlie-output.txt
        - timestamped_Alice-output.txt
    ```

Para hacerlo componible con otros workflows, hay algunas cosas que necesitan cambiar.

### 1.2. Hacer el workflow componible

Para hacer un workflow componible, tres cosas deben cambiar:
el workflow recibe un nombre, las entradas se mueven a un bloque `take:` y las salidas se mueven a un bloque `emit:`
(reemplazando los bloques independientes `publish:`/`output {}`, que pertenecen al entry workflow).

Las siguientes secciones cubren estos cambios uno por uno.

#### 1.2.1. Nombrar el workflow

Dale un nombre al workflow para que pueda importarse desde un workflow padre.

=== "Después"

    ```groovy title="workflows/greeting.nf" linenums="5" hl_lines="1"
    workflow GREETING_WORKFLOW {
    ```

=== "Antes"

    ```groovy title="workflows/greeting.nf" linenums="5" hl_lines="1"
    workflow {
    ```

Con un nombre, el workflow puede importarse en otros scripts.

#### 1.2.2. Declarar entradas con `take:`

Reemplaza la declaración de canal codificada directamente con un bloque `take:` que declare qué entradas espera el workflow.
El bloque `take:` va antes de `main:`, y se elimina la línea `names_ch = channel.of(...)`.

=== "Después"

    ```groovy title="workflows/greeting.nf" linenums="5" hl_lines="2 3 5"
    workflow GREETING_WORKFLOW {
        take:
        names_ch // Canal de entrada con nombres

        main:
        // Encadena los procesos: validar -> crear saludo -> agregar marca de tiempo
        validated_ch = VALIDATE_NAME(names_ch)
        greetings_ch = SAY_HELLO(validated_ch)
        timestamped_ch = TIMESTAMP_GREETING(greetings_ch)
    ```

=== "Antes"

    ```groovy title="workflows/greeting.nf" linenums="5"
    workflow GREETING_WORKFLOW {
        main:
        names_ch = channel.of('Alice', 'Bob', 'Charlie')

        // Encadena los procesos: validar -> crear saludo -> agregar marca de tiempo
        validated_ch = VALIDATE_NAME(names_ch)
        greetings_ch = SAY_HELLO(validated_ch)
        timestamped_ch = TIMESTAMP_GREETING(greetings_ch)
    ```

El bloque `take:` declara el canal solo por nombre. El workflow padre define qué contiene.

#### 1.2.3. Declarar salidas con `emit:`

Reemplaza la sección `publish:` y elimina el bloque `output {}`, sustituyéndolos por un bloque `emit:` que nombra las salidas.

=== "Después"

    ```groovy title="workflows/greeting.nf" linenums="14" hl_lines="2 3 4"

        emit:
        greetings = greetings_ch // Saludos originales
        timestamped = timestamped_ch // Saludos con marca de tiempo
    }
    ```

=== "Antes"

    ```groovy title="workflows/greeting.nf" linenums="14"

        publish:
        greetings = greetings_ch
        timestamped = timestamped_ch
    }

    output {
        greetings {
        }
        timestamped {
        }
    }
    ```

El bloque `emit:` expone salidas con nombre a las que los workflows padre pueden acceder mediante `GREETING_WORKFLOW.out.greetings` y `GREETING_WORKFLOW.out.timestamped`.

#### 1.2.4. Verificar el resultado y probarlo

Después de los tres cambios, el archivo completo debería verse así:

```groovy title="workflows/greeting.nf" linenums="1" hl_lines="7 8 9 11 17 18 19"
#!/usr/bin/env nextflow

include { VALIDATE_NAME } from '../modules/validate_name'
include { SAY_HELLO } from '../modules/say_hello'
include { TIMESTAMP_GREETING } from '../modules/timestamp_greeting'

workflow GREETING_WORKFLOW {
    take:
    names_ch // Canal de entrada con nombres

    main:
    // Encadena los procesos: validar -> crear saludo -> agregar marca de tiempo
    validated_ch = VALIDATE_NAME(names_ch)
    greetings_ch = SAY_HELLO(validated_ch)
    timestamped_ch = TIMESTAMP_GREETING(greetings_ch)

    emit:
    greetings = greetings_ch // Saludos originales
    timestamped = timestamped_ch // Saludos con marca de tiempo
}
```

Ahora intenta ejecutarlo directamente:

```bash
nextflow run workflows/greeting.nf
```

??? failure "Salida del comando"

    ```console
    N E X T F L O W   ~  version 26.04.4
    Launching `workflows/greeting.nf` [intergalactic_leavitt] revision: 89ac88c6c2
    No entry workflow specified
    ```

Esto introduce un concepto clave: el **entry workflow**.
Nextflow usa un bloque `workflow {}` sin nombre como punto de entrada cuando ejecutas un script directamente.
`GREETING_WORKFLOW` tiene nombre, por lo que Nextflow no sabe cómo ejecutarlo por sí solo.

Esto es intencional. Los workflows componibles están diseñados para ser llamados desde un entry workflow, no para ejecutarse directamente.
La solución es un entry workflow en `main.nf` que importe y llame a `GREETING_WORKFLOW`.

### 1.3. Actualizar y probar el workflow principal

Ahora actualiza el workflow principal para llamar al greeting workflow.

#### 1.3.1. Incluir el greeting workflow y llamarlo

Agrega la sentencia `include`, actualiza el cuerpo del workflow para llamar a `GREETING_WORKFLOW` y reemplaza el marcador de posición `channel.empty()` en `publish:`:

=== "Después"

    ```groovy title="main.nf" linenums="1" hl_lines="3 9 10 13"
    #!/usr/bin/env nextflow

    include { GREETING_WORKFLOW } from './workflows/greeting'

    workflow {
        main:
        names = channel.of('Alice', 'Bob', 'Charlie')

        // Ejecuta el greeting workflow
        GREETING_WORKFLOW(names)

        publish:
        greetings = GREETING_WORKFLOW.out.greetings
    }
    ```

=== "Antes"

    ```groovy title="main.nf" linenums="1"
    #!/usr/bin/env nextflow

    workflow {
        main:
        names = channel.of('Alice', 'Bob', 'Charlie')

        publish:
        greetings = channel.empty()
    }
    ```

El entry workflow permanece sin nombre para que Nextflow lo use como punto de entrada del pipeline.

#### 1.3.2. Actualizar el bloque output

Agrega una directiva `path` para enrutar los saludos publicados hacia un subdirectorio `greetings/`:

=== "Después"

    ```groovy title="main.nf" linenums="14" hl_lines="3"
    output {
        greetings {
            path 'greetings'
        }
    }
    ```

=== "Antes"

    ```groovy title="main.nf" linenums="14" hl_lines="2 3"
    output {
        greetings {
        }
    }
    ```

#### 1.3.3. Ejecutar el workflow

Ejecuta el workflow para verificar que funciona:

```bash
nextflow run main.nf
```

??? success "Salida del comando"

    ```console
    N E X T F L O W   ~  version 26.04.4
    Launching `main.nf` [agitated_jones] revision: 6b0e98dd05
    executor >  local (9)
    [00/5ac4ce] GRE…DATE_NAME (validating Bob) | 3 of 3 ✔
    [7d/dde45c] GRE…SAY_HELLO (greeting Alice) | 3 of 3 ✔
    [fd/a044cc] GRE…ing timestamp to greeting) | 3 of 3 ✔

    Outputs:

      /workspaces/training/side-quests/workflows_of_workflows/results

      greetings:
        - greetings/Bob-output.txt
        - greetings/Charlie-output.txt
        - greetings/Alice-output.txt
    ```

??? abstract "Nuevo contenido agregado en `results/`"

    ```console
    results/
    └── greetings
        ├── Alice-output.txt
        ├── Bob-output.txt
        └── Charlie-output.txt
    ```

    Los seis archivos de la ejecución autónoma en la sección 1.1 (`Alice-output.txt`, `timestamped_Alice-output.txt`, etc.) siguen estando en `results/` junto a este nuevo subdirectorio `greetings/`.
    Esto es esperado: nada en esta lección los elimina, y la sección 2.1 los lee directamente como entrada.

??? abstract "Contenido del archivo"

    ```console title="results/greetings/Alice-output.txt"
    Hello, Alice!
    ```

Los archivos de saludo se publican en `results/greetings/`.
El workflow principal llama a `GREETING_WORKFLOW` y conecta su salida directamente a la sección `publish:`.

`GREETING_WORKFLOW.out.timestamped` no está conectado a `publish:` aquí.
A partir de la sección 2, ese canal se convierte en la entrada de `TRANSFORM_WORKFLOW` en lugar de ser una salida publicada por sí misma.
Contrasta esto con la ejecución autónoma en la sección 1.1, donde `timestamped` no tenía nada más a qué alimentar, por lo que se publicó directamente.

### Conclusión

En esta sección, has aprendido varios conceptos importantes:

- **Workflows con nombre**: Crear un workflow con nombre (`GREETING_WORKFLOW`) que puede importarse y reutilizarse
- **Interfaces de workflow**: Definir entradas claras con `take:` y salidas con `emit:` para crear un workflow componible
- **Entry points**: Entender que Nextflow necesita un entry workflow sin nombre para ejecutar un script
- **Composición de workflows**: Importar y usar un workflow con nombre dentro de otro workflow
- **Namespaces de workflow**: Acceder a las salidas del workflow usando la notación `.out` (`GREETING_WORKFLOW.out.greetings`)

Ahora tienes un greeting workflow funcional que:

- Recibe un canal de nombres como entrada
- Valida cada nombre
- Crea un saludo para cada nombre válido
- Agrega marcas de tiempo a los saludos
- Emite tanto los saludos originales como los saludos con marca de tiempo, aunque solo los saludos originales se publican en `results/` en esta etapa

Este enfoque modular te permite probar el greeting workflow de forma independiente o usarlo como componente en pipelines más grandes.

---

## 2. Agregar el transform workflow al pipeline

El transform workflow aplica transformaciones de texto a los saludos con marca de tiempo.

### 2.1. Revisar y ejecutar el workflow

Abre `workflows/transform.nf` y examina el código:

```groovy title="workflows/transform.nf" linenums="1"
#!/usr/bin/env nextflow

include { SAY_HELLO_UPPER } from '../modules/say_hello_upper'
include { REVERSE_TEXT } from '../modules/reverse_text'

workflow {
    main:
    input_ch = channel.fromPath('results/timestamped_*.txt')

    // Aplica las transformaciones en secuencia
    upper_ch = SAY_HELLO_UPPER(input_ch)
    reversed_ch = REVERSE_TEXT(upper_ch)

    publish:
    upper = upper_ch
    reversed = reversed_ch
}

output {
    upper {
    }
    reversed {
    }
}
```

Este workflow independiente lee archivos de saludos con marca de tiempo del directorio `results/` producido por `greeting.nf`, los convierte a mayúsculas y luego invierte el texto.

Ejecútalo para verificar que funciona con los resultados del greeting de la sección 1.1:

```bash
nextflow run workflows/transform.nf
```

??? success "Salida del comando"

    ```console
    N E X T F L O W   ~  version 26.04.4
    Launching `workflows/transform.nf` [evil_meninsky] revision: 4e35790d32
    executor >  local (6)
    [7d/4ed7de] SAY…estamped_Alice-output.txt) | 3 of 3 ✔
    [70/ddaafa] REV…imestamped_Bob-output.txt) | 3 of 3 ✔

    Outputs:

      /workspaces/training/side-quests/workflows_of_workflows/results

      upper:
        - UPPER-timestamped_Bob-output.txt
        - UPPER-timestamped_Charlie-output.txt
        - UPPER-timestamped_Alice-output.txt

      reversed:
        - REVERSED-UPPER-timestamped_Charlie-output.txt
        - REVERSED-UPPER-timestamped_Alice-output.txt
        - REVERSED-UPPER-timestamped_Bob-output.txt
    ```

Para hacerlo componible con `GREETING_WORKFLOW`, se aplican los mismos tres cambios de la sección 1.2.

### 2.2. Hacerlo componible

Aplica los mismos tres cambios que en la sección 1.2: nombra el workflow, reemplaza la entrada codificada directamente con `take:`, y reemplaza `publish:`/`output {}` con `emit:`.

El archivo terminado debería verse así:

```groovy title="workflows/transform.nf" linenums="1" hl_lines="6 7 8 10 15 16 17"
#!/usr/bin/env nextflow

include { SAY_HELLO_UPPER } from '../modules/say_hello_upper'
include { REVERSE_TEXT } from '../modules/reverse_text'

workflow TRANSFORM_WORKFLOW {
    take:
    input_ch // Canal de entrada con saludos

    main:
    // Aplica las transformaciones en secuencia
    upper_ch = SAY_HELLO_UPPER(input_ch)
    reversed_ch = REVERSE_TEXT(upper_ch)

    emit:
    upper = upper_ch // Saludos en mayúsculas
    reversed = reversed_ch // Saludos en mayúsculas invertidos
}
```

El transform workflow ahora es componible y está listo para importarse en el workflow principal.

### 2.3. Actualizar y probar el workflow principal

Ahora actualiza el workflow principal para llamar al transform workflow.

#### 2.3.1. Incluir el transform workflow y llamarlo

Agrega la sentencia include, una llamada a `TRANSFORM_WORKFLOW` encadenada sobre los saludos con marca de tiempo, y las dos nuevas entradas en `publish:`:

=== "Después"

    ```groovy title="main.nf" linenums="1" hl_lines="4 13 14 18 19"
    #!/usr/bin/env nextflow

    include { GREETING_WORKFLOW } from './workflows/greeting'
    include { TRANSFORM_WORKFLOW } from './workflows/transform'

    workflow {
        main:
        names = channel.of('Alice', 'Bob', 'Charlie')

        // Ejecuta el greeting workflow
        GREETING_WORKFLOW(names)

        // Ejecuta el transform workflow
        TRANSFORM_WORKFLOW(GREETING_WORKFLOW.out.timestamped)

        publish:
        greetings = GREETING_WORKFLOW.out.greetings
        upper = TRANSFORM_WORKFLOW.out.upper
        reversed = TRANSFORM_WORKFLOW.out.reversed
    }
    ```

=== "Antes"

    ```groovy title="main.nf" linenums="1"
    #!/usr/bin/env nextflow

    include { GREETING_WORKFLOW } from './workflows/greeting'

    workflow {
        main:
        names = channel.of('Alice', 'Bob', 'Charlie')

        // Ejecuta el greeting workflow
        GREETING_WORKFLOW(names)

        publish:
        greetings = GREETING_WORKFLOW.out.greetings
    }
    ```

Esto ejecutará el transform workflow sobre los saludos con marca de tiempo.

#### 2.3.2. Actualizar el bloque output

Agrega las entradas `upper` y `reversed` al bloque `output {}`, cada una con una directiva `path` para su subdirectorio:

=== "Después"

    ```groovy title="main.nf" linenums="20" hl_lines="5 6 7 8 9 10"
    output {
        greetings {
            path 'greetings'
        }
        upper {
            path 'upper'
        }
        reversed {
            path 'reversed'
        }
    }
    ```

=== "Antes"

    ```groovy title="main.nf" linenums="20" hl_lines="2 3 4 5"
    output {
        greetings {
            path 'greetings'
        }
    }
    ```

Esto publicará las salidas finales en los directorios correspondientes.

#### 2.3.3. Ejecutar el pipeline completo

Ejecuta el pipeline para verificar que todo funciona:

```bash
nextflow run main.nf
```

??? success "Salida del comando"

    ```console
    N E X T F L O W   ~  version 26.04.4
    Launching `main.nf` [special_magritte] revision: 2a0e54fecd
    executor >  local (15)
    [e7/e051a5] GRE…TE_NAME (validating Alice) | 3 of 3 ✔
    [ca/cf1b18] GRE…SAY_HELLO (greeting Alice) | 3 of 3 ✔
    [e5/aa6be2] GRE…ing timestamp to greeting) | 3 of 3 ✔
    [d2/7afc88] TRA…imestamped_Bob-output.txt) | 3 of 3 ✔
    [5b/70f190] TRA…estamped_Alice-output.txt) | 3 of 3 ✔

    Outputs:

      /workspaces/training/side-quests/workflows_of_workflows/results

      greetings:
        - greetings/Bob-output.txt
        - greetings/Charlie-output.txt
        - greetings/Alice-output.txt

      upper:
        - upper/UPPER-timestamped_Charlie-output.txt
        - upper/UPPER-timestamped_Bob-output.txt
        - upper/UPPER-timestamped_Alice-output.txt

      reversed:
        - reversed/REVERSED-UPPER-timestamped_Charlie-output.txt
        - reversed/REVERSED-UPPER-timestamped_Bob-output.txt
        - reversed/REVERSED-UPPER-timestamped_Alice-output.txt
    ```

??? abstract "Nuevo contenido agregado en `results/`"

    ```console
    results/
    ├── greetings
    │   ├── Alice-output.txt
    │   ├── Bob-output.txt
    │   └── Charlie-output.txt
    ├── reversed
    │   ├── REVERSED-UPPER-timestamped_Alice-output.txt
    │   ├── REVERSED-UPPER-timestamped_Bob-output.txt
    │   └── REVERSED-UPPER-timestamped_Charlie-output.txt
    └── upper
        ├── UPPER-timestamped_Alice-output.txt
        ├── UPPER-timestamped_Bob-output.txt
        └── UPPER-timestamped_Charlie-output.txt
    ```

    La raíz de `results/` también sigue conteniendo los archivos sobrantes de las ejecuciones autónomas en las secciones 1.1 y 2.1 (`Alice-output.txt`, `timestamped_Alice-output.txt`, `UPPER-timestamped_Alice-output.txt`, etc.).
    Esto es esperado: esta lección nunca indica que debes eliminarlos.

??? abstract "Contenido del archivo"

    ```console title="results/reversed/REVERSED-UPPER-timestamped_Alice-output.txt"
    !ECILA ,OLLEH ]54:84:81 31-90-6202[
    ```

El pipeline funciona de extremo a extremo: el saludo ha sido convertido a mayúsculas e invertido.

### Conclusión

Ahora deberías tener un pipeline completo que:

- Procesa nombres a través del greeting workflow
- Alimenta los saludos con marca de tiempo al transform workflow
- Produce versiones en mayúsculas e invertidas de los saludos

---

## Resumen

En esta misión secundaria, hemos explorado el poderoso concepto de composición de workflows en Nextflow, que nos permite construir pipelines complejos a partir de componentes más pequeños y reutilizables.

Este enfoque modular ofrece varias ventajas sobre los pipelines monolíticos:

- Cada workflow puede desarrollarse, probarse y depurarse de forma independiente
- Los workflows pueden reutilizarse en diferentes pipelines
- La estructura general del pipeline se vuelve más legible y fácil de mantener
- Los cambios en un workflow no necesariamente afectan a otros si las interfaces permanecen consistentes
- Los entry points pueden configurarse para ejecutar diferentes partes de tu pipeline según sea necesario

Llamar a un workflow es similar a llamar a un proceso, pero no es lo mismo.
No puedes ejecutar un workflow N veces pasándole un canal de tamaño N: pasas el canal una vez y el workflow itera sobre él internamente.

Aplicar estas técnicas en tu propio trabajo te permitirá construir pipelines de Nextflow más sofisticados que puedan manejar tareas complejas de procesamiento de datos mientras se mantienen fáciles de mantener y escalar.

### Patrones clave

1.  **Estructura del workflow**: Definimos entradas y salidas claras para cada workflow usando la sintaxis `take:` y `emit:`, creando interfaces bien definidas entre componentes, y envolvimos la lógica del workflow dentro del bloque `main:`.

    ```groovy
    workflow EXAMPLE_WORKFLOW {
        take:
        // Los canales de entrada se declaran aquí
        input_ch

        main:
        // La lógica del workflow va aquí
        // Aquí es donde se llaman los procesos y se manipulan los canales
        result_ch = SOME_PROCESS(input_ch)

        emit:
        // Los canales de salida se declaran aquí
        output_ch = result_ch
    }
    ```

2.  **Importaciones de workflow**: Construimos dos módulos de workflow independientes y los importamos a un pipeline principal con sentencias `include`.

    - Incluir un único workflow

    ```groovy
    include { WORKFLOW_NAME } from './path/to/workflow'
    ```

    - Incluir múltiples workflows

    ```groovy
    include { WORKFLOW_A; WORKFLOW_B } from './path/to/workflows'
    ```

    - Incluir con alias para evitar conflictos de nombres

    ```groovy
    include { WORKFLOW_A as WORKFLOW_A_ALIAS } from './path/to/workflow'
    ```

3.  **Entry points**: Nextflow requiere un entry workflow sin nombre para saber dónde comenzar la ejecución.
    Este entry workflow llama a tus workflows con nombre.

    - Workflow sin nombre (entry point)

    ```groovy
    workflow {
        // Este es el entry point cuando se ejecuta el script
        NAMED_WORKFLOW(input_ch)
    }
    ```

    - Workflow con nombre (llamado desde el entry workflow)

    ```groovy
    workflow NAMED_WORKFLOW {
        // Debe ser llamado desde el entry workflow
    }
    ```

4.  **Gestión del flujo de datos**: Aprendimos cómo acceder a las salidas del workflow usando la notación de namespace (`WORKFLOW_NAME.out.channel_name`) y pasarlas a otros workflows.

    ```groovy
    WORKFLOW_A(input_ch)
    WORKFLOW_B(WORKFLOW_A.out.some_channel)
    ```

### Recursos adicionales

- [Documentación de Workflows de Nextflow](https://www.nextflow.io/docs/latest/workflow.html)
- [Referencia de operadores de canales](https://www.nextflow.io/docs/latest/operator.html)
- [Documentación de Error Strategy](https://www.nextflow.io/docs/latest/process.html#errorstrategy)

---

## ¿Qué sigue?

Regresa al [menú de misiones secundarias](../index.md) o haz clic en el botón en la parte inferior derecha de la página para continuar con el siguiente tema de la lista.
