# Flujos de trabajo de flujos de trabajo

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Traducción asistida por IA - [más información y sugerencias](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

Cuando está desarrollando un pipeline, a menudo se encuentra creando secuencias similares de procesos para diferentes tipos de datos o pasos de análisis. Podría terminar copiando y pegando estas secuencias de procesos, lo que lleva a código duplicado que es difícil de mantener; o podría crear un flujo de trabajo masivo que es difícil de entender y modificar.

Una de las características más poderosas de Nextflow es su capacidad para componer pipelines complejos a partir de módulos de flujos de trabajo más pequeños y reutilizables. Este enfoque modular hace que los pipelines sean más fáciles de desarrollar, probar y mantener.

### Objetivos de aprendizaje

En esta misión secundaria, exploraremos cómo desarrollar módulos de flujos de trabajo que se pueden probar y usar por separado, componer esos módulos en un pipeline más grande y gestionar el flujo de datos entre módulos.

Al final de esta misión secundaria, usted podrá:

- Dividir pipelines complejos en unidades lógicas y reutilizables
- Probar cada módulo de flujo de trabajo de forma independiente
- Mezclar y combinar flujos de trabajo para crear nuevos pipelines
- Compartir módulos de flujos de trabajo comunes entre diferentes pipelines
- Hacer su código más mantenible y más fácil de entender

Estas habilidades le ayudarán a construir pipelines complejos mientras mantiene una estructura de código limpia y mantenible.

### Requisitos previos

Antes de asumir esta misión secundaria debería:

- Haber completado el tutorial [Hello Nextflow](../hello_nextflow/README.md) o un curso equivalente para principiantes.
- Sentirse cómodo usando conceptos y mecanismos básicos de Nextflow (procesos, canales, operadores, módulos)

---

## 0. Comenzar

#### Abrir el codespace de entrenamiento

Si aún no lo ha hecho, asegúrese de abrir el entorno de entrenamiento como se describe en [Configuración del Entorno](../envsetup/index.md).

[![Open in GitHub Codespaces](https://github.com/codespaces/badge.svg)](https://codespaces.new/nextflow-io/training?quickstart=1&ref=master)

#### Moverse al directorio del proyecto

Vamos a movernos al directorio donde se encuentran los archivos para este tutorial.

```bash
cd side-quests/workflows_of_workflows
```

Puede configurar VSCode para enfocarse en este directorio:

```bash
code .
```

#### Revisar los materiales

Encontrará un directorio `modules` que contiene varias definiciones de procesos que se basan en lo que aprendió en 'Hello Nextflow':

```console title="Contenido del directorio"
modules/
├── say_hello.nf             # Crea un saludo (de Hello Nextflow)
├── say_hello_upper.nf       # Convierte a mayúsculas (de Hello Nextflow)
├── timestamp_greeting.nf    # Agrega marcas de tiempo a los saludos
├── validate_name.nf         # Valida nombres de entrada
└── reverse_text.nf          # Invierte el contenido del texto
```

#### Revisar la asignación

Su desafío es ensamblar estos módulos en dos flujos de trabajo separados que luego compondremos en un flujo de trabajo principal:

- Un `GREETING_WORKFLOW` que valida nombres, crea saludos y agrega marcas de tiempo
- Un `TRANSFORM_WORKFLOW` que convierte texto a mayúsculas y lo invierte

#### Lista de verificación de preparación

¿Cree que está listo para comenzar?

- [ ] Entiendo el objetivo de este curso y sus requisitos previos
- [ ] Mi codespace está funcionando
- [ ] He configurado mi directorio de trabajo apropiadamente
- [ ] Entiendo la asignación

Si puede marcar todas las casillas, está listo para comenzar.

---

## 1. Crear el flujo de trabajo de saludo

Comencemos creando un flujo de trabajo que valide nombres y genere saludos con marcas de tiempo.

### 1.1. Crear la estructura del flujo de trabajo

```bash title="Crear directorio de flujo de trabajo y archivo"
mkdir -p workflows
touch workflows/greeting.nf
```

### 1.2. Agregar el código del primer (sub)flujo de trabajo

Agregue este código a `workflows/greeting.nf`:

```groovy title="workflows/greeting.nf" linenums="1"
include { VALIDATE_NAME } from '../modules/validate_name'
include { SAY_HELLO } from '../modules/say_hello'
include { TIMESTAMP_GREETING } from '../modules/timestamp_greeting'

workflow {

    names_ch = channel.of('Alice', 'Bob', 'Charlie')

    // Encadenar procesos: validar -> crear saludo -> agregar marca de tiempo
    validated_ch = VALIDATE_NAME(names_ch)
    greetings_ch = SAY_HELLO(validated_ch)
    timestamped_ch = TIMESTAMP_GREETING(greetings_ch)
}
```

Este es un flujo de trabajo completo, con una estructura similar a los que vio en el tutorial 'Hello Nextflow', que podemos probar de forma independiente. Probémoslo ahora:

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

### 1.3. Hacer el flujo de trabajo componible

Los flujos de trabajo componibles tienen algunas diferencias con respecto a los que vio en el tutorial 'Hello Nextflow':

- El bloque workflow necesita tener un nombre
- Las entradas se declaran usando la palabra clave `take:`
- El contenido del flujo de trabajo se coloca dentro del bloque `main:`
- Las salidas se declaran usando la palabra clave `emit:`

Actualicemos el flujo de trabajo de saludo para que coincida con esta estructura. Cambie el código a lo siguiente:

```groovy title="workflows/greeting.nf" linenums="1" hl_lines="6 7 9 15 16 17"
include { VALIDATE_NAME } from '../modules/validate_name'
include { SAY_HELLO } from '../modules/say_hello'
include { TIMESTAMP_GREETING } from '../modules/timestamp_greeting'

workflow GREETING_WORKFLOW {
    take:
        names_ch        // Canal de entrada con nombres

    main:
        // Encadenar procesos: validar -> crear saludo -> agregar marca de tiempo
        validated_ch = VALIDATE_NAME(names_ch)
        greetings_ch = SAY_HELLO(validated_ch)
        timestamped_ch = TIMESTAMP_GREETING(greetings_ch)

    emit:
        greetings = greetings_ch      // Saludos originales
        timestamped = timestamped_ch  // Saludos con marca de tiempo
}
```

Puede ver que el flujo de trabajo ahora tiene un nombre y tiene un bloque `take:` y `emit:`, y estas son las conexiones que usaremos para componer un flujo de trabajo de nivel superior.
El contenido del flujo de trabajo también se coloca dentro del bloque `main:`. Note también que hemos eliminado la declaración del canal de entrada `names_ch`, ya que ahora se pasa como argumento al flujo de trabajo.

Probemos el flujo de trabajo nuevamente para ver si funciona como se esperaba:

```bash
nextflow run workflows/greeting.nf
```

??? failure "Salida del comando"

    ```console
    N E X T F L O W  ~  version 24.10.0
    Launching `workflows/greeting.nf` [high_brahmagupta] DSL2 - revision: 8f5857af25
    No entry workflow specified
    ```

Esto le informa sobre otro concepto nuevo, un 'flujo de trabajo de entrada'. El flujo de trabajo de entrada es el flujo de trabajo que se llama cuando ejecuta un script de Nextflow. Por defecto, Nextflow usará un flujo de trabajo sin nombre como flujo de trabajo de entrada, cuando esté presente, y eso es lo que ha estado haciendo hasta ahora, con bloques de flujo de trabajo que comienzan así:

```groovy title="hello.nf" linenums="1"
workflow {
```

Pero nuestro flujo de trabajo de saludo no tiene un flujo de trabajo sin nombre, sino que tenemos un flujo de trabajo con nombre:

```groovy title="workflows/greeting.nf" linenums="1"
workflow GREETING_WORKFLOW {
```

Por eso Nextflow arrojó un error y no hizo lo que queríamos.

No agregamos la sintaxis `take:`/`emit:` para poder llamar al flujo de trabajo directamente - lo hicimos para poder componerlo con otros flujos de trabajo. La solución es crear un script principal con un flujo de trabajo de entrada sin nombre que importe y llame a nuestro flujo de trabajo con nombre.

### 1.4. Crear y probar el flujo de trabajo principal

Ahora crearemos un flujo de trabajo principal que importa y usa el flujo de trabajo de saludo.

Crear `main.nf`:

```groovy title="main.nf" linenums="1"
include { GREETING_WORKFLOW } from './workflows/greeting'

workflow {
    names = channel.of('Alice', 'Bob', 'Charlie')
    GREETING_WORKFLOW(names)

    GREETING_WORKFLOW.out.greetings.view { "Original: $it" }
    GREETING_WORKFLOW.out.timestamped.view { "Timestamped: $it" }
}

```

Note que nuestra entrada de flujo de trabajo en este archivo no tiene nombre, y eso es porque vamos a usarla como un flujo de trabajo de entrada.

Ejecute esto y vea la salida:

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

¡Funciona! Hemos envuelto el flujo de trabajo de saludo con nombre en un flujo de trabajo principal con un bloque de entrada `workflow` sin nombre. El flujo de trabajo principal está usando el flujo de trabajo `GREETING_WORKFLOW` casi (no exactamente) como un proceso, y pasando el canal `names` como argumento.

### Conclusión

En esta sección, ha aprendido varios conceptos importantes:

- **Flujos de trabajo con nombre**: Crear un flujo de trabajo con nombre (`GREETING_WORKFLOW`) que puede importarse y reutilizarse
- **Interfaces de flujo de trabajo**: Definir entradas claras con `take:` y salidas con `emit:` para crear un flujo de trabajo componible
- **Puntos de entrada**: Entender que Nextflow necesita un flujo de trabajo de entrada sin nombre para ejecutar un script
- **Composición de flujos de trabajo**: Importar y usar un flujo de trabajo con nombre dentro de otro flujo de trabajo
- **Espacios de nombres de flujo de trabajo**: Acceder a las salidas del flujo de trabajo usando el espacio de nombres `.out` (`GREETING_WORKFLOW.out.greetings`)

Ahora tiene un flujo de trabajo de saludo funcional que:

- Toma un canal de nombres como entrada
- Valida cada nombre
- Crea un saludo para cada nombre válido
- Agrega marcas de tiempo a los saludos
- Expone tanto los saludos originales como los que tienen marca de tiempo como salidas

Este enfoque modular le permite probar el flujo de trabajo de saludo de forma independiente o usarlo como componente en pipelines más grandes.

---

## 2. Agregar el flujo de trabajo de transformación

Ahora creemos un flujo de trabajo que aplica transformaciones de texto a los saludos.

### 2.1. Crear el archivo de flujo de trabajo

```bash
touch workflows/transform.nf
```

### 2.2. Agregar el código del flujo de trabajo

Agregue este código a `workflows/transform.nf`:

```groovy title="workflows/transform.nf" linenums="1"
include { SAY_HELLO_UPPER } from '../modules/say_hello_upper'
include { REVERSE_TEXT } from '../modules/reverse_text'

workflow TRANSFORM_WORKFLOW {
    take:
        input_ch         // Canal de entrada con mensajes

    main:
        // Aplicar transformaciones en secuencia
        upper_ch = SAY_HELLO_UPPER(input_ch)
        reversed_ch = REVERSE_TEXT(upper_ch)

    emit:
        upper = upper_ch        // Saludos en mayúsculas
        reversed = reversed_ch  // Saludos en mayúsculas invertidos
}
```

No repetiremos la explicación de la sintaxis componible aquí, pero note que el flujo de trabajo con nombre se declara nuevamente con un bloque `take:` y `emit:`, y el contenido del flujo de trabajo se coloca dentro del bloque `main:`.

### 2.3. Actualizar el flujo de trabajo principal

Actualice `main.nf` para usar ambos flujos de trabajo:

```groovy title="main.nf" linenums="1"
include { GREETING_WORKFLOW } from './workflows/greeting'
include { TRANSFORM_WORKFLOW } from './workflows/transform'

workflow {
    names = channel.of('Alice', 'Bob', 'Charlie')

    // Ejecutar el flujo de trabajo de saludo
    GREETING_WORKFLOW(names)

    // Ejecutar el flujo de trabajo de transformación
    TRANSFORM_WORKFLOW(GREETING_WORKFLOW.out.timestamped)

    // Ver resultados
    TRANSFORM_WORKFLOW.out.upper.view { "Uppercase: $it" }
    TRANSFORM_WORKFLOW.out.reversed.view { "Reversed: $it" }
}
```

Ejecute el pipeline completo:

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

Si echa un vistazo a uno de esos archivos invertidos, verá que es la versión en mayúsculas del saludo invertido:

```bash
cat /workspaces/training/side_quests/workflows_of_workflows/work/f0/74ba4a10d9ef5c82f829d1c154d0f6/REVERSED-UPPER-timestamped_Alice-output.txt
```

```console title="Contenido del archivo invertido"
!ECILA ,OLLEH ]04:50:71 60-30-5202[
```

### Conclusión

Ahora debería tener un pipeline completo que:

- Procesa nombres a través del flujo de trabajo de saludo
- Alimenta los saludos con marca de tiempo al flujo de trabajo de transformación
- Produce versiones tanto en mayúsculas como invertidas de los saludos

---

## Resumen

En esta misión secundaria, hemos explorado el poderoso concepto de composición de flujos de trabajo en Nextflow, que nos permite construir pipelines complejos a partir de componentes más pequeños y reutilizables.

Este enfoque modular ofrece varias ventajas sobre los pipelines monolíticos:

- Cada flujo de trabajo puede desarrollarse, probarse y depurarse de forma independiente
- Los flujos de trabajo pueden reutilizarse en diferentes pipelines
- La estructura general del pipeline se vuelve más legible y mantenible
- Los cambios en un flujo de trabajo no necesariamente afectan a otros si las interfaces permanecen consistentes
- Los puntos de entrada pueden configurarse para ejecutar diferentes partes de su pipeline según sea necesario

_Es importante notar sin embargo que aunque llamar flujos de trabajo es un poco como llamar procesos, en realidad no es lo mismo. No puede, por ejemplo, ejecutar un flujo de trabajo N veces llamándolo con un canal de tamaño N - necesitaría pasar un canal de tamaño N al flujo de trabajo e iterar internamente._

Aplicar estas técnicas en su propio trabajo le permitirá construir pipelines de Nextflow más sofisticados que puedan manejar tareas bioinformáticas complejas mientras permanecen mantenibles y escalables.

### Patrones clave

1.  **Estructura de flujo de trabajo**: Definimos entradas y salidas claras para cada flujo de trabajo usando la sintaxis `take:` y `emit:`, creando interfaces bien definidas entre componentes, y envolvimos la lógica del flujo de trabajo dentro del bloque `main:`.

    ```groovy
    workflow EXAMPLE_WORKFLOW {
        take:
            // Los canales de entrada se declaran aquí
            input_ch

        main:
            // La lógica del flujo de trabajo va aquí
            // Aquí es donde se llaman los procesos y se manipulan los canales
            result_ch = SOME_PROCESS(input_ch)

        emit:
            // Los canales de salida se declaran aquí
            output_ch = result_ch
    }
    ```

2.  **Importaciones de flujos de trabajo:** Construimos dos módulos de flujos de trabajo independientes y los importamos en un pipeline principal con declaraciones include.

    - Incluir un único flujo de trabajo

    ```groovy
    include { WORKFLOW_NAME } from './path/to/workflow'
    ```

    - Incluir múltiples flujos de trabajo

    ```groovy
    include { WORKFLOW_A; WORKFLOW_B } from './path/to/workflows'
    ```

    - Incluir con alias para evitar conflictos de nombres

    ```groovy
    include { WORKFLOW_A as WORKFLOW_A_ALIAS } from './path/to/workflow'
    ```

3.  **Puntos de entrada**: Nextflow requiere un flujo de trabajo de entrada sin nombre para saber dónde comenzar la ejecución. Este flujo de trabajo de entrada llama a sus flujos de trabajo con nombre.

    - Flujo de trabajo sin nombre (punto de entrada)

    ```groovy
    workflow {
        // Este es el punto de entrada cuando se ejecuta el script
        NAMED_WORKFLOW(input_ch)
    }
    ```

    - Flujo de trabajo con nombre (llamado desde el flujo de trabajo de entrada)

    ```groovy
    workflow NAMED_WORKFLOW {
        // Debe ser llamado desde el flujo de trabajo de entrada
    }
    ```

4.  **Gestión del flujo de datos:** Aprendimos cómo acceder a las salidas del flujo de trabajo usando la notación de espacio de nombres (`WORKFLOW_NAME.out.channel_name`) y pasarlas a otros flujos de trabajo.

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

Regrese al [menú de Misiones Secundarias](./index.md) o haga clic en el botón en la parte inferior derecha de la página para pasar al siguiente tema de la lista.
