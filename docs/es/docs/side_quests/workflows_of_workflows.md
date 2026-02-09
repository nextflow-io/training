# Workflows de Workflows

Cuando estás desarrollando un pipeline, a menudo te encuentras creando secuencias similares de procesos para diferentes tipos de datos o pasos de análisis. Podrías terminar copiando y pegando estas secuencias de procesos, lo que lleva a código duplicado que es difícil de mantener; o podrías crear un workflow masivo que es difícil de entender y modificar.

Una de las características más poderosas de Nextflow es su capacidad para componer pipelines complejos a partir de módulos de workflow más pequeños y reutilizables. Este enfoque modular hace que los pipelines sean más fáciles de desarrollar, probar y mantener.

### Objetivos de aprendizaje

En esta misión secundaria, exploraremos cómo desarrollar módulos de workflow que pueden ser probados y usados por separado, componer esos módulos en un pipeline más grande, y gestionar el flujo de datos entre módulos.

Al final de esta misión secundaria, serás capaz de:

- Dividir pipelines complejos en unidades lógicas y reutilizables
- Probar cada módulo de workflow de forma independiente
- Mezclar y combinar workflows para crear nuevos pipelines
- Compartir módulos de workflow comunes entre diferentes pipelines
- Hacer tu código más mantenible y fácil de entender

Estas habilidades te ayudarán a construir pipelines complejos mientras mantienes una estructura de código limpia y mantenible.

### Requisitos previos

Antes de emprender esta misión secundaria deberías:

- Haber completado el tutorial [Hello Nextflow](../hello_nextflow/README.md) o un curso equivalente para principiantes.
- Sentirte cómodo usando conceptos y mecanismos básicos de Nextflow (procesos, canales, operadores, módulos)

---

## 0. Primeros pasos

#### Abrir el codespace de capacitación

Si aún no lo has hecho, asegúrate de abrir el entorno de capacitación como se describe en [Configuración del entorno](../envsetup/index.md).

[![Open in GitHub Codespaces](https://github.com/codespaces/badge.svg)](https://codespaces.new/nextflow-io/training?quickstart=1&ref=master)

#### Moverse al directorio del proyecto

Vamos a movernos al directorio donde se encuentran los archivos para este tutorial.

```bash
cd side-quests/workflows_of_workflows
```

Puedes configurar VSCode para enfocarse en este directorio:

```bash
code .
```

#### Revisar los materiales

Encontrarás un directorio `modules` que contiene varias definiciones de procesos que se basan en lo que aprendiste en 'Hello Nextflow':

```console title="Directory contents"
modules/
├── say_hello.nf             # Creates a greeting (from Hello Nextflow)
├── say_hello_upper.nf       # Converts to uppercase (from Hello Nextflow)
├── timestamp_greeting.nf    # Adds timestamps to greetings
├── validate_name.nf         # Validates input names
└── reverse_text.nf          # Reverses text content
```

#### Revisar la asignación

Tu desafío es ensamblar estos módulos en dos workflows separados que luego compondremos en un workflow principal:

- Un `GREETING_WORKFLOW` que valida nombres, crea saludos y añade marcas de tiempo
- Un `TRANSFORM_WORKFLOW` que convierte texto a mayúsculas y lo invierte

#### Lista de verificación de preparación

¿Crees que estás listo para comenzar?

- [ ] Entiendo el objetivo de este curso y sus requisitos previos
- [ ] Mi codespace está funcionando
- [ ] He configurado mi directorio de trabajo apropiadamente
- [ ] Entiendo la asignación

Si puedes marcar todas las casillas, estás listo para comenzar.

---

## 1. Crear el Greeting Workflow

Comencemos creando un workflow que valida nombres y genera saludos con marcas de tiempo.

### 1.1. Crear la estructura del workflow

```bash title="Create workflow directory and file"
mkdir -p workflows
touch workflows/greeting.nf
```

### 1.2. Añadir el código del primer (sub)workflow

Añade este código a `workflows/greeting.nf`:

```groovy title="workflows/greeting.nf" linenums="1"
include { VALIDATE_NAME } from '../modules/validate_name'
include { SAY_HELLO } from '../modules/say_hello'
include { TIMESTAMP_GREETING } from '../modules/timestamp_greeting'

workflow {

    names_ch = channel.of('Alice', 'Bob', 'Charlie')

    // Chain processes: validate -> create greeting -> add timestamp
    validated_ch = VALIDATE_NAME(names_ch)
    greetings_ch = SAY_HELLO(validated_ch)
    timestamped_ch = TIMESTAMP_GREETING(greetings_ch)
}
```

Este es un workflow completo, con una estructura similar a los que viste en el tutorial 'Hello Nextflow', que podemos probar de forma independiente. Probémoslo ahora:

```bash
nextflow run workflows/greeting.nf
```

??? success "Salida del comando"

    ```console
    N E X T F L O W  ~  version 24.10.0
    Launching `workflows/greeting.nf` [peaceful_montalcini] DSL2 - revision: 90f61b7093
    executor >  local (9)
    [51/4f980f] process > VALIDATE_NAME (validating Bob)                    [100%] 3 of 3 ✔
    [2b/dd8dc2] process > SAY_HELLO (greeting Bob)                          [100%] 3 of 3 ✔
    [8e/882565] process > TIMESTAMP_GREETING (adding timestamp to greeting) [100%] 3 of 3 ✔
    ```

Esto funciona como se esperaba, pero para hacerlo componible hay algunas cosas que necesitamos cambiar.

### 1.3. Hacer el workflow componible

Los workflows componibles tienen algunas diferencias con respecto a los que viste en el tutorial 'Hello Nextflow':

- El bloque workflow necesita tener un nombre
- Las entradas se declaran usando la palabra clave `take:`
- El contenido del workflow se coloca dentro del bloque `main:`
- Las salidas se declaran usando la palabra clave `emit:`

Actualicemos el greeting workflow para que coincida con esta estructura. Cambia el código a lo siguiente:

```groovy title="workflows/greeting.nf" linenums="1" hl_lines="6 7 9 15 16 17"
include { VALIDATE_NAME } from '../modules/validate_name'
include { SAY_HELLO } from '../modules/say_hello'
include { TIMESTAMP_GREETING } from '../modules/timestamp_greeting'

workflow GREETING_WORKFLOW {
    take:
        names_ch        // Input channel with names

    main:
        // Chain processes: validate -> create greeting -> add timestamp
        validated_ch = VALIDATE_NAME(names_ch)
        greetings_ch = SAY_HELLO(validated_ch)
        timestamped_ch = TIMESTAMP_GREETING(greetings_ch)

    emit:
        greetings = greetings_ch      // Original greetings
        timestamped = timestamped_ch  // Timestamped greetings
}
```

Puedes ver que el workflow ahora tiene nombre y tiene un bloque `take:` y `emit:`, y estas son las conexiones que usaremos para componer un workflow de nivel superior.
El contenido del workflow también se coloca dentro del bloque `main:`. Nota también que hemos eliminado la declaración del canal de entrada `names_ch`, ya que ahora se pasa como argumento al workflow.

Probemos el workflow nuevamente para ver si funciona como se esperaba:

```bash
nextflow run workflows/greeting.nf
```

??? failure "Salida del comando"

    ```console
    N E X T F L O W  ~  version 24.10.0
    Launching `workflows/greeting.nf` [high_brahmagupta] DSL2 - revision: 8f5857af25
    No entry workflow specified
    ```

Esto te informa sobre otro concepto nuevo, un 'entry workflow'. El entry workflow es el workflow que se llama cuando ejecutas un script de Nextflow. Por defecto, Nextflow usará un workflow sin nombre como entry workflow, cuando esté presente, y eso es lo que has estado haciendo hasta ahora, con bloques workflow que comienzan así:

```groovy title="hello.nf" linenums="1"
workflow {
```

Pero nuestro greeting workflow no tiene un workflow sin nombre, sino que tenemos un workflow con nombre:

```groovy title="workflows/greeting.nf" linenums="1"
workflow GREETING_WORKFLOW {
```

Por eso Nextflow lanzó un error y no hizo lo que queríamos.

No añadimos la sintaxis `take:`/`emit:` para poder llamar al workflow directamente - lo hicimos para poder componerlo con otros workflows. La solución es crear un script principal con un entry workflow sin nombre que importe y llame a nuestro workflow con nombre.

### 1.4. Crear y probar el workflow principal

Ahora crearemos un workflow principal que importa y usa el workflow `greeting`.

Crea `main.nf`:

```groovy title="main.nf" linenums="1"
include { GREETING_WORKFLOW } from './workflows/greeting'

workflow {
    names = channel.of('Alice', 'Bob', 'Charlie')
    GREETING_WORKFLOW(names)

    GREETING_WORKFLOW.out.greetings.view { "Original: $it" }
    GREETING_WORKFLOW.out.timestamped.view { "Timestamped: $it" }
}

```

Nota que nuestra entrada de workflow en este archivo no tiene nombre, y eso es porque vamos a usarla como entry workflow.

Ejecuta esto y observa la salida:

```bash
nextflow run main.nf
```

??? success "Salida del comando"

    ```console
    N E X T F L O W  ~  version 24.10.0
    Launching `main.nf` [goofy_mayer] DSL2 - revision: 543f8742fe
    executor >  local (9)
    [05/3cc752] process > GREETING_WORKFLOW:VALIDATE_NAME (validating Char... [100%] 3 of 3 ✔
    [b1/b56ecf] process > GREETING_WORKFLOW:SAY_HELLO (greeting Charlie)      [100%] 3 of 3 ✔
    [ea/342168] process > GREETING_WORKFLOW:TIMESTAMP_GREETING (adding tim... [100%] 3 of 3 ✔
    Original: /workspaces/training/side_quests/workflows_of_workflows/work/bb/c8aff3df0ebc15a4d7d35f736db44c/Alice-output.txt
    Original: /workspaces/training/side_quests/workflows_of_workflows/work/fb/fa877776e8a5d90b537b1bcd3b6f5b/Bob-output.txt
    Original: /workspaces/training/side_quests/workflows_of_workflows/work/b1/b56ecf938fda8bcbec211847c8f0be/Charlie-output.txt
    Timestamped: /workspaces/training/side_quests/workflows_of_workflows/work/06/877bc909f140bbf8223343450cea36/timestamped_Alice-output.txt
    Timestamped: /workspaces/training/side_quests/workflows_of_workflows/work/aa/bd31b71cdb745b7c155ca7f8837b8a/timestamped_Bob-output.txt
    Timestamped: /workspaces/training/side_quests/workflows_of_workflows/work/ea/342168d4ba04cc899a89c56cbfd9b0/timestamped_Charlie-output.txt
    ```

¡Funciona! Hemos envuelto el greeting workflow con nombre en un workflow principal con un bloque `workflow` de entrada sin nombre. El workflow principal está usando el workflow `GREETING_WORKFLOW` casi (no exactamente) como un proceso, y pasando el canal `names` como argumento.

### Conclusión

En esta sección, has aprendido varios conceptos importantes:

- **Workflows con nombre**: Crear un workflow con nombre (`GREETING_WORKFLOW`) que puede ser importado y reutilizado
- **Interfaces de workflow**: Definir entradas claras con `take:` y salidas con `emit:` para crear un workflow componible
- **Puntos de entrada**: Entender que Nextflow necesita un entry workflow sin nombre para ejecutar un script
- **Composición de workflows**: Importar y usar un workflow con nombre dentro de otro workflow
- **Espacios de nombres de workflow**: Acceder a las salidas del workflow usando el espacio de nombres `.out` (`GREETING_WORKFLOW.out.greetings`)

Ahora tienes un greeting workflow funcional que:

- Toma un canal de nombres como entrada
- Valida cada nombre
- Crea un saludo para cada nombre válido
- Añade marcas de tiempo a los saludos
- Expone tanto los saludos originales como los que tienen marca de tiempo como salidas

Este enfoque modular te permite probar el greeting workflow de forma independiente o usarlo como componente en pipelines más grandes.

---

## 2. Añadir el Transform Workflow

Ahora creemos un workflow que aplica transformaciones de texto a los saludos.

### 2.1. Crear el archivo del workflow

```bash
touch workflows/transform.nf
```

### 2.2. Añadir el código del workflow

Añade este código a `workflows/transform.nf`:

```groovy title="workflows/transform.nf" linenums="1"
include { SAY_HELLO_UPPER } from '../modules/say_hello_upper'
include { REVERSE_TEXT } from '../modules/reverse_text'

workflow TRANSFORM_WORKFLOW {
    take:
        input_ch         // Input channel with messages

    main:
        // Apply transformations in sequence
        upper_ch = SAY_HELLO_UPPER(input_ch)
        reversed_ch = REVERSE_TEXT(upper_ch)

    emit:
        upper = upper_ch        // Uppercase greetings
        reversed = reversed_ch  // Reversed uppercase greetings
}
```

No repetiremos la explicación de la sintaxis componible aquí, pero nota que el workflow con nombre se declara nuevamente con un bloque `take:` y `emit:`, y el contenido del workflow se coloca dentro del bloque `main:`.

### 2.3. Actualizar el workflow principal

Actualiza `main.nf` para usar ambos workflows:

```groovy title="main.nf" linenums="1"
include { GREETING_WORKFLOW } from './workflows/greeting'
include { TRANSFORM_WORKFLOW } from './workflows/transform'

workflow {
    names = channel.of('Alice', 'Bob', 'Charlie')

    // Run the greeting workflow
    GREETING_WORKFLOW(names)

    // Run the transform workflow
    TRANSFORM_WORKFLOW(GREETING_WORKFLOW.out.timestamped)

    // View results
    TRANSFORM_WORKFLOW.out.upper.view { "Uppercase: $it" }
    TRANSFORM_WORKFLOW.out.reversed.view { "Reversed: $it" }
}
```

Ejecuta el pipeline completo:

```bash
nextflow run main.nf
```

??? success "Salida del comando"

    ```console
    N E X T F L O W  ~  version 24.10.0
    Launching `main.nf` [sick_kimura] DSL2 - revision: 8dc45fc6a8
    executor >  local (13)
    executor >  local (15)
    [83/1b51f4] process > GREETING_WORKFLOW:VALIDATE_NAME (validating Alice)  [100%] 3 of 3 ✔
    [68/556150] process > GREETING_WORKFLOW:SAY_HELLO (greeting Alice)        [100%] 3 of 3 ✔
    [de/511abd] process > GREETING_WORKFLOW:TIMESTAMP_GREETING (adding tim... [100%] 3 of 3 ✔
    [cd/e6a7e0] process > TRANSFORM_WORKFLOW:SAY_HELLO_UPPER (converting t... [100%] 3 of 3 ✔
    [f0/74ba4a] process > TRANSFORM_WORKFLOW:REVERSE_TEXT (reversing UPPER... [100%] 3 of 3 ✔
    Uppercase: /workspaces/training/side_quests/workflows_of_workflows/work/a0/d4f5df4d6344604498fa47a6084a11/UPPER-timestamped_Bob-output.txt
    Uppercase: /workspaces/training/side_quests/workflows_of_workflows/work/69/b5e37f6c79c2fd38adb75d0eca8f87/UPPER-timestamped_Charlie-output.txt
    Uppercase: /workspaces/training/side_quests/workflows_of_workflows/work/cd/e6a7e0b17e7d5a2f71bb8123cd53a7/UPPER-timestamped_Alice-output.txt
    Reversed: /workspaces/training/side_quests/workflows_of_workflows/work/7a/7a222f7957b35d1d121338566a24ac/REVERSED-UPPER-timestamped_Bob-output.txt
    Reversed: /workspaces/training/side_quests/workflows_of_workflows/work/46/8d19af62e33a5a6417c773496e0f90/REVERSED-UPPER-timestamped_Charlie-output.txt
    Reversed: /workspaces/training/side_quests/workflows_of_workflows/work/f0/74ba4a10d9ef5c82f829d1c154d0f6/REVERSED-UPPER-timestamped_Alice-output.txt
    ```

Si echas un vistazo a uno de esos archivos invertidos, verás que es la versión en mayúsculas del saludo invertido:

```bash
cat /workspaces/training/side_quests/workflows_of_workflows/work/f0/74ba4a10d9ef5c82f829d1c154d0f6/REVERSED-UPPER-timestamped_Alice-output.txt
```

```console title="Reversed file content"
!ECILA ,OLLEH ]04:50:71 60-30-5202[
```

### Conclusión

Ahora deberías tener un pipeline completo que:

- Procesa nombres a través del greeting workflow
- Alimenta los saludos con marca de tiempo al transform workflow
- Produce versiones tanto en mayúsculas como invertidas de los saludos

---

## Resumen

En esta misión secundaria, hemos explorado el poderoso concepto de composición de workflows en Nextflow, que nos permite construir pipelines complejos a partir de componentes más pequeños y reutilizables.

Este enfoque modular ofrece varias ventajas sobre los pipelines monolíticos:

- Cada workflow puede ser desarrollado, probado y depurado de forma independiente
- Los workflows pueden ser reutilizados en diferentes pipelines
- La estructura general del pipeline se vuelve más legible y mantenible
- Los cambios en un workflow no necesariamente afectan a otros si las interfaces permanecen consistentes
- Los puntos de entrada pueden configurarse para ejecutar diferentes partes de tu pipeline según sea necesario

_Es importante notar sin embargo que aunque llamar workflows es un poco como llamar procesos, no es realmente lo mismo. No puedes, por ejemplo, ejecutar un workflow N veces llamándolo con un canal de tamaño N - necesitarías pasar un canal de tamaño N al workflow e iterar internamente._

Aplicar estas técnicas en tu propio trabajo te permitirá construir pipelines de Nextflow más sofisticados que pueden manejar tareas bioinformáticas complejas mientras permanecen mantenibles y escalables.

### Patrones clave

1.  **Estructura de workflow**: Definimos entradas y salidas claras para cada workflow usando la sintaxis `take:` y `emit:`, creando interfaces bien definidas entre componentes, y envolvimos la lógica del workflow dentro del bloque `main:`.

    ```groovy
    workflow EXAMPLE_WORKFLOW {
        take:
            // Input channels are declared here
            input_ch

        main:
            // Workflow logic goes here
            // This is where processes are called and channels are manipulated
            result_ch = SOME_PROCESS(input_ch)

        emit:
            // Output channels are declared here
            output_ch = result_ch
    }
    ```

2.  **Importaciones de workflow:** Construimos dos módulos de workflow independientes y los importamos en un pipeline principal con declaraciones include.

    - Incluir un solo workflow

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

3.  **Puntos de entrada**: Nextflow requiere un entry workflow sin nombre para saber dónde comenzar la ejecución. Este entry workflow llama a tus workflows con nombre.

    - Workflow sin nombre (punto de entrada)

    ```groovy
    workflow {
        // This is the entry point when the script is run
        NAMED_WORKFLOW(input_ch)
    }
    ```

    - Workflow con nombre (llamado desde el entry workflow)

    ```groovy
    workflow NAMED_WORKFLOW {
        // Must be called from the entry workflow
    }
    ```

4.  **Gestión del flujo de datos:** Aprendimos cómo acceder a las salidas del workflow usando la notación de espacio de nombres (`WORKFLOW_NAME.out.channel_name`) y pasarlas a otros workflows.

    ```nextflow
    WORKFLOW_A(input_ch)
    WORKFLOW_B(WORKFLOW_A.out.some_channel)
    ```

### Recursos adicionales

- [Documentación de Workflow de Nextflow](https://www.nextflow.io/docs/latest/workflow.html)
- [Referencia de Operadores de Canal](https://www.nextflow.io/docs/latest/operator.html)
- [Documentación de Error Strategy](https://www.nextflow.io/docs/latest/process.html#errorstrategy)

---

## ¿Qué sigue?

Regresa al [menú de Misiones Secundarias](./index.md) o haz clic en el botón en la parte inferior derecha de la página para continuar con el siguiente tema de la lista.
