# Pruebas con nf-test

Poder probar sistemáticamente que cada parte de tu flujo de trabajo está haciendo lo que se supone que debe hacer es fundamental para la reproducibilidad y el mantenimiento a largo plazo, y puede ser de gran ayuda durante el proceso de desarrollo.

Tomemos un minuto para hablar sobre por qué las pruebas son tan importantes. Si estás desarrollando un flujo de trabajo, una de las primeras cosas que harás es tomar algunos datos de prueba que sabes que son válidos y deberían producir un resultado. Agregas el primer proceso al pipeline y lo conectas a tus entradas para que funcione. Luego, para verificar que todo está funcionando, lo ejecutas con los datos de prueba. Asumiendo que funciona, pasas al siguiente proceso y ejecutas los datos de prueba nuevamente. Repites este proceso hasta que tienes un pipeline con el que estás satisfecho.

Luego, tal vez agregues un parámetro simple de verdadero o falso como `--skip_process`. Ahora debes ejecutar el pipeline dos veces, una con cada parámetro para asegurarte de que funciona como se espera. Pero espera, ¿cómo verificamos si `--skip_process` realmente omite el proceso? ¡Tenemos que revisar las salidas o verificar los archivos de registro! Esto es tedioso y propenso a errores.

A medida que desarrollas tu pipeline, rápidamente se volverá tan complejo que probar manualmente cada iteración será lento y propenso a errores. Además, si encuentras un error, será muy difícil identificar exactamente de dónde en tu pipeline proviene el error. Aquí es donde entran las pruebas.

Las pruebas te permiten verificar sistemáticamente que cada parte de tu pipeline está funcionando como se espera. Los beneficios para un desarrollador de pruebas bien escritas son enormes:

- **Confianza**: Debido a que las pruebas cubren todo el pipeline, puedes estar seguro de que cambiar algo no afecta nada más
- **Confiabilidad**: Cuando múltiples desarrolladores trabajan en el pipeline, saben que los otros desarrolladores no han roto el pipeline ni ningún componente.
- **Transparencia**: Las pruebas muestran dónde está fallando un pipeline y facilitan el rastreo del problema. También funcionan como una forma de documentación, mostrando cómo ejecutar un proceso o workflow.
- **Velocidad**: Debido a que las pruebas están automatizadas, se pueden ejecutar muy rápidamente y repetidamente. Puedes iterar rápidamente con menos temor de introducir nuevos errores.

Hay muchos tipos diferentes de pruebas que podemos escribir:

1. **Pruebas a nivel de módulo**: Para procesos individuales
2. **Pruebas a nivel de workflow**: Para un solo workflow
3. **Pruebas a nivel de pipeline**: Para el pipeline en su conjunto
4. **Pruebas de rendimiento**: Para la velocidad y eficiencia del pipeline
5. **Pruebas de estrés**: Evaluando el rendimiento del pipeline bajo condiciones extremas para determinar sus límites

Probar procesos individuales es análogo a las pruebas unitarias en otros lenguajes. Probar el workflow o el pipeline completo es análogo a lo que se llama pruebas de integración en otros lenguajes, donde probamos las interacciones de los componentes.

[**nf-test**](https://www.nf-test.com/) es una herramienta que te permite escribir pruebas a nivel de módulo, workflow y pipeline. En resumen, te permite verificar sistemáticamente que cada parte individual del pipeline está funcionando como se espera, _de forma aislada_.

### Objetivos de aprendizaje

En esta misión secundaria, aprenderás a usar nf-test para escribir una prueba a nivel de workflow para el pipeline, así como pruebas a nivel de módulo para los tres procesos que invoca.

Al final de esta misión secundaria, podrás usar las siguientes técnicas de manera efectiva:

- Inicializar nf-test en tu proyecto
- Generar pruebas a nivel de módulo y workflow
- Agregar tipos comunes de aserciones
- Comprender cuándo usar snapshots vs. aserciones de contenido
- Ejecutar pruebas para un proyecto completo

Estas habilidades te ayudarán a implementar una estrategia de pruebas integral en tus proyectos de pipeline, asegurando que sean más robustos y mantenibles.

### Requisitos previos

Antes de emprender esta misión secundaria, deberías:

- Haber completado el tutorial [Hello Nextflow](../hello_nextflow/README.md) o un curso equivalente para principiantes.
- Sentirte cómodo usando conceptos y mecanismos básicos de Nextflow (procesos, canales, operadores, trabajar con archivos, metadatos)

---

## 0. Primeros pasos

#### Abrir el codespace de capacitación

Si aún no lo has hecho, asegúrate de abrir el entorno de capacitación como se describe en [Configuración del entorno](../envsetup/index.md).

[![Open in GitHub Codespaces](https://github.com/codespaces/badge.svg)](https://codespaces.new/nextflow-io/training?quickstart=1&ref=master)

#### Moverse al directorio del proyecto

Vamos a movernos al directorio donde se encuentran los archivos para este tutorial.

```bash
cd side-quests/nf-test
```

Puedes configurar VSCode para enfocarse en este directorio:

```bash
code .
```

#### Revisar los materiales

Encontrarás un archivo de workflow principal y un archivo CSV llamado `greetings.csv` que contiene la entrada al pipeline.

```console title="Directory contents"
.
├── greetings.csv
└── main.nf
```

Para una descripción detallada de los archivos, consulta el [calentamiento de Hello Nextflow](../hello_nextflow/00_orientation.md).

El workflow que probaremos es un subconjunto del workflow Hello construido en [Hello Workflow](../hello_nextflow/03_hello_workflow.md).

??? example "¿Qué hace el workflow Hello Nextflow?"

    Si no has realizado la capacitación [Hello Nextflow](../hello_nextflow/index.md), aquí hay una descripción rápida de lo que hace este workflow simple.

    El workflow toma un archivo CSV que contiene saludos, ejecuta cuatro pasos de transformación consecutivos sobre ellos y genera un único archivo de texto que contiene una imagen ASCII de un personaje divertido diciendo los saludos.

    Los cuatro pasos se implementan como procesos de Nextflow (`sayHello`, `convertToUpper`, `collectGreetings` y `cowpy`) almacenados en archivos de módulo separados.

    1. **`sayHello`:** Escribe cada saludo en su propio archivo de salida (por ejemplo, "Hello-output.txt")
    2. **`convertToUpper`:** Convierte cada saludo a mayúsculas (por ejemplo, "HELLO")
    3. **`collectGreetings`:** Recopila todos los saludos en mayúsculas en un único archivo por lotes
    4. **`cowpy`:** Genera arte ASCII usando la herramienta `cowpy`

    Los resultados se publican en un directorio llamado `results/`, y la salida final del pipeline (cuando se ejecuta con parámetros predeterminados) es un archivo de texto plano que contiene arte ASCII de un personaje diciendo los saludos en mayúsculas.

    En esta misión secundaria, usamos una forma intermedia del workflow Hello que solo contiene los dos primeros procesos. <!-- TODO: change this to use the full finished workflow as suggested in https://github.com/nextflow-io/training/issues/735 -->

El subconjunto con el que trabajaremos está compuesto por dos procesos: `sayHello` y `convertToUpper`.
Puedes ver el código completo del workflow a continuación.

??? example "Código del workflow"

    ```groovy title="main.nf"
    /*
    * Parámetros del pipeline
    */
    params.input_file = "greetings.csv"

    /*
    * Usa echo para imprimir 'Hello World!' a la salida estándar
    */
    process sayHello {

        publishDir 'results', mode: 'copy'

        input:
            val greeting

        output:
            path "${greeting}-output.txt"

        script:
        """
        echo '$greeting' > '$greeting-output.txt'
        """
    }

    /*
    * Usa una utilidad de reemplazo de texto para convertir el saludo a mayúsculas
    */
    process convertToUpper {

        publishDir 'results', mode: 'copy'

        input:
            path input_file

        output:
            path "UPPER-${input_file}"

        script:
        """
        cat '$input_file' | tr '[a-z]' '[A-Z]' > UPPER-${input_file}
        """
    }

    workflow {

        // create a channel for inputs from a CSV file
        greeting_ch = channel.fromPath(params.input_file).splitCsv().flatten()

        // emit a greeting
        sayHello(greeting_ch)

        // convert the greeting to uppercase
        convertToUpper(sayHello.out)
    }
    ```

#### Ejecutar el workflow

Ejecutemos el workflow para asegurarnos de que está funcionando como se espera.

```bash
nextflow run main.nf
```

```console title="Result of running the workflow"
 N E X T F L O W   ~  version 24.10.2

Launching `main.nf` [soggy_linnaeus] DSL2 - revision: bbf79d5c31

executor >  local (6)
[f7/c3be66] sayHello (3)       | 3 of 3 ✔
[cd/e15303] convertToUpper (3) | 3 of 3 ✔
```

¡FELICIDADES! ¡Acabas de ejecutar una prueba!

"Espera, ¿qué? ¡Solo ejecuté el workflow y funcionó! ¿Cómo es eso una prueba?"

¡Buena pregunta!

Analicemos lo que acaba de suceder.

Ejecutaste el workflow con los parámetros predeterminados, confirmaste que funcionó y estás satisfecho con los resultados. Esta es la esencia de las pruebas. Si trabajaste en el curso de capacitación Hello Nextflow, habrás notado que siempre comenzamos cada sección ejecutando el workflow que estábamos usando como punto de partida, para confirmar que todo está configurado correctamente.

Probar software esencialmente hace este proceso por nosotros.

#### Revisar la asignación

Tu desafío es agregar pruebas estandarizadas a este workflow usando nf-test, para facilitar la verificación de que cada parte continúa funcionando como se espera en caso de que se realicen más cambios.

<!-- TODO: give a bit more details, similar to how it's done in the Metadata side quest -->

#### Lista de verificación de preparación

¿Crees que estás listo para comenzar?

- [ ] Entiendo el objetivo de este curso y sus requisitos previos
- [ ] Mi codespace está funcionando
- [ ] He configurado mi directorio de trabajo apropiadamente
- [ ] He ejecutado el workflow exitosamente
- [ ] Entiendo la asignación

Si puedes marcar todas las casillas, estás listo para comenzar.

---

## 1. Inicializar `nf-test`

El paquete `nf-test` proporciona un comando de inicialización que configura algunas cosas para que podamos comenzar a desarrollar pruebas para nuestro proyecto.

```bash
nf-test init
```

Esto debería producir la siguiente salida:

```bash
🚀 nf-test 0.9.3
https://www.nf-test.com
(c) 2021 - 2024 Lukas Forer and Sebastian Schoenherr

Project configured. Configuration is stored in nf-test.config
```

También crea un directorio `tests` que contiene un esqueleto de archivo de configuración.

### 1.1. Generar un nf-test

`nf-test` viene con un conjunto de herramientas para construir archivos nf-test, ahorrándonos la mayor parte del trabajo. Estas vienen bajo el subcomando `generate`. Generemos una prueba para el pipeline:

```bash
nf-test generate pipeline main.nf
```

```console title="Output"
> nf-test generate pipeline main.nf

Load source file '/workspaces/training/side-quests/nf-test/main.nf'
Wrote pipeline test file '/workspaces/training/side-quests/nf-test/tests/main.nf.test

SUCCESS: Generated 1 test files.
```

Esto creará un archivo `main.nf.test` dentro del directorio `tests`. Este es nuestro archivo de prueba a nivel de pipeline. Si ejecutas `tree tests/` deberías ver algo como esto:

```console title="Test directory contents"
tests/
├── main.nf.test
└── nextflow.config
```

El archivo `main.nf.test` es nuestro archivo de prueba a nivel de pipeline. Ábrelo y echemos un vistazo al contenido.

```groovy title="tests/main.nf.test"
nextflow_pipeline {

    name "Test Workflow main.nf"
    script "main.nf"

    test("Should run without failures") {

        when {
            params {
                // define parameters here. Example:
                // outdir = "tests/results"
            }
        }

        then {
            assert workflow.success
        }

    }

}
```

Tomaremos un segundo para entender la estructura del archivo de prueba.

El bloque `nextflow_pipeline` es el punto de entrada para todas las pruebas a nivel de pipeline. Contiene lo siguiente:

- `name`: El nombre de la prueba.
- `script`: La ruta al script del pipeline.

El bloque `test` es la prueba real. Contiene lo siguiente:

- `when`: Las condiciones bajo las cuales se debe ejecutar la prueba. Esto incluye los parámetros que se usarán para ejecutar el pipeline.
- `then`: Las aserciones que se deben hacer. Esto incluye los resultados esperados del pipeline.

En lenguaje sencillo, la lógica de la prueba se lee de la siguiente manera:
"**Cuando** estos _parámetros_ se proporcionan a este _pipeline_, **entonces** esperamos ver estos resultados."

Esta no es una prueba funcional, demostraremos cómo convertirla en una en la siguiente sección.

### Una nota sobre los nombres de las pruebas

En el ejemplo anterior, usamos el nombre predeterminado "Should run without failures" que es apropiado para una prueba básica que solo verifica si el pipeline se ejecuta exitosamente. Sin embargo, a medida que agregamos casos de prueba más específicos, debemos usar nombres más descriptivos que indiquen lo que realmente estamos probando. Por ejemplo:

- "Should convert input to uppercase" - cuando se prueba funcionalidad específica
- "Should handle empty input gracefully" - cuando se prueban casos extremos
- "Should respect max memory parameter" - cuando se prueban restricciones de recursos
- "Should create expected output files" - cuando se prueba la generación de archivos

Los buenos nombres de prueba deben:

1. Comenzar con "Should" para dejar claro cuál es el comportamiento esperado
2. Describir la funcionalidad o escenario específico que se está probando
3. Ser lo suficientemente claros para que si la prueba falla, sepas qué funcionalidad está rota

A medida que agreguemos más aserciones y casos de prueba específicos más adelante, usaremos estos nombres más descriptivos para dejar claro qué está verificando cada prueba.

### 1.2. Ejecutar la prueba

Ejecutemos la prueba para ver qué sucede.

```bash
nf-test test tests/main.nf.test
```

```console title="nf-test pipeline fail"
> nf-test test tests/main.nf.test

🚀 nf-test 0.9.3
https://www.nf-test.com
(c) 2021 - 2024 Lukas Forer and Sebastian Schoenherr


Test Workflow main.nf

  Test [693ba951] 'Should run without failures' FAILED (4.652s)

  Assertion failed:

  assert workflow.success
         |        |
         workflow false

  Nextflow stdout:

  ERROR ~ No such file or directory: /workspaces/training/side-quests/nf-test/.nf-test/tests/693ba951a20fec36a5a9292ed1cc8a9f/greetings.csv

   -- Check '/workspaces/training/side-quests/nf-test/.nf-test/tests/693ba951a20fec36a5a9292ed1cc8a9f/meta/nextflow.log' file for details
  Nextflow stderr:

FAILURE: Executed 1 tests in 4.679s (1 failed)
```

¡La prueba falla! ¿Qué sucedió?

1. nf-test intentó ejecutar el pipeline tal como está, usando la configuración en el bloque `when`:

```groovy title="tests/main.nf.test"
when {
    params {
        // define parameters here. Example:
        // outdir = "tests/results"
    }
}
```

2. nf-test verificó el estado del pipeline y lo comparó con el bloque `when`:

```groovy title="tests/main.nf.test"
then {
    assert workflow.success
}
```

Nota cómo nf-test ha reportado que el pipeline falló y proporcionó el mensaje de error de Nextflow:

```console title="Error"
ERROR ~ No such file or directory: /workspaces/training/side-quests/nf-test/.nf-test/tests/693ba951a20fec36a5a9292ed1cc8a9f/greetings.csv
```

Entonces, ¿cuál fue el problema? Recuerda que el pipeline tiene un archivo greetings.csv en el directorio del proyecto. Cuando nf-test ejecuta el pipeline, buscará este archivo, pero no puede encontrarlo. El archivo está ahí, ¿qué está pasando? Bueno, si miramos la ruta podemos ver que la prueba está ocurriendo en la ruta `./nf-test/tests/longHashString/`. Al igual que Nextflow, nf-test crea un nuevo directorio para cada prueba para mantener todo aislado. El archivo de datos no está ubicado allí, por lo que debemos corregir la ruta al archivo en la prueba original.

Puede que te estés preguntando cómo vamos a apuntar a la raíz del pipeline en la prueba. Dado que esta es una situación común, nf-test tiene un rango de variables globales que podemos usar para facilitarnos la vida. Puedes encontrar la lista completa [aquí](https://www.nf-test.com/docs/testcases/global_variables/) pero mientras tanto usaremos la variable `projectDir`, que significa la raíz del proyecto del pipeline.

_Antes:_

```groovy title="tests/main.nf.test" linenums="1" hl_lines="3 4"
when {
    params {
        // define parameters here. Example:
        // outdir = "tests/results"
    }
}
```

_Después:_

```groovy title="tests/main.nf.test" linenums="1" hl_lines="3"
when {
    params {
        input_file = "${projectDir}/greetings.csv"
    }
}
```

Ejecutemos la prueba nuevamente para ver si funciona.

```bash title="nf-test pipeline pass"
nf-test test tests/main.nf.test
```

```console title="Pipeline passes"
> nf-test test tests/main.nf.test

🚀 nf-test 0.9.3
https://www.nf-test.com
(c) 2021 - 2024 Lukas Forer and Sebastian Schoenherr


Test Workflow main.nf

  Test [1d4aaf12] 'Should run without failures' PASSED (1.619s)


SUCCESS: Executed 1 tests in 1.626s
```

¡Éxito! El pipeline se ejecuta exitosamente y la prueba pasa. ¡Ejecútalo tantas veces como quieras y siempre obtendrás el mismo resultado!

Por defecto, la salida de Nextflow está oculta, pero para convencerte de que nf-test definitivamente está ejecutando el workflow, puedes usar la bandera `--verbose`:

```bash
nf-test test tests/main.nf.test --verbose
```

```console title="Pipeline runs all processes"
> nf-test test tests/main.nf.test

🚀 nf-test 0.9.3
https://www.nf-test.com
(c) 2021 - 2024 Lukas Forer and Sebastian Schoenherr


Test Workflow main.nf

  Test [693ba951] 'Should run without failures'
    > Nextflow 24.10.4 is available - Please consider updating your version to it
    > N E X T F L O W  ~  version 24.10.0
    > Launching `/workspaces/training/side-quests/nf-test/main.nf` [zen_ampere] DSL2 - revision: bbf79d5c31
    > [2b/61e453] Submitted process > sayHello (2)
    > [31/4e1606] Submitted process > sayHello (1)
    > [bb/5209ee] Submitted process > sayHello (3)
    > [83/83db6f] Submitted process > convertToUpper (2)
    > [9b/3428b1] Submitted process > convertToUpper (1)
    > [ca/0ba51b] Submitted process > convertToUpper (3)
    PASSED (5.206s)


SUCCESS: Executed 1 tests in 5.239s
```

### 1.3. Agregar aserciones

Una verificación simple es asegurar que nuestro pipeline está ejecutando todos los procesos que esperamos y no está omitiendo ninguno silenciosamente. Recuerda que nuestro pipeline ejecuta 6 procesos, uno llamado `sayHello` y uno llamado `convertToUpper` para cada uno de los 3 saludos.

Agreguemos una aserción a nuestra prueba para verificar que el pipeline ejecuta el número esperado de procesos. También actualizaremos el nombre de nuestra prueba para reflejar mejor lo que estamos probando.

**Antes:**

```groovy title="tests/main.nf.test" linenums="1" hl_lines="1"
    test("Should run without failures") {

        when {
            params {
                input_file = "${projectDir}/greetings.csv"
            }
        }

        then {
            assert workflow.success
        }

    }
```

**Después:**

```groovy title="tests/main.nf.test" linenums="1" hl_lines="1 11"
    test("Should run successfully with correct number of processes") {

        when {
            params {
                input_file = "${projectDir}/greetings.csv"
            }
        }

        then {
            assert workflow.success
            assert workflow.trace.tasks().size() == 6
        }

    }
```

El nombre de la prueba ahora refleja mejor lo que realmente estamos verificando: no solo que el pipeline se ejecuta sin fallar, sino que ejecuta el número esperado de procesos.

Ejecutemos la prueba nuevamente para ver si funciona.

```bash title="nf-test pipeline pass"
nf-test test tests/main.nf.test
```

```console title="Pipeline passes with assertions"
🚀 nf-test 0.9.3
https://www.nf-test.com
(c) 2021 - 2024 Lukas Forer and Sebastian Schoenherr


Test Workflow main.nf

  Test [1d4aaf12] 'Should run successfully with correct number of processes' PASSED (1.567s)


SUCCESS: Executed 1 tests in 1.588s
```

¡Éxito! El pipeline se ejecuta exitosamente y la prueba pasa. Ahora hemos comenzado a probar los detalles del pipeline, así como el estado general.

### 1.4. Probar la salida

Agreguemos una aserción a nuestra prueba para verificar que se creó el archivo de salida. Lo agregaremos como una prueba separada, con un nombre informativo, para que los resultados sean más fáciles de interpretar.

**Antes:**

```groovy title="tests/main.nf.test" linenums="1" hl_lines="14"
    test("Should run successfully with correct number of processes") {

        when {
            params {
                input_file = "${projectDir}/greetings.csv"
            }
        }

        then {
            assert workflow.success
            assert workflow.trace.tasks().size() == 6
        }

    }
```

**Después:**

```groovy title="tests/main.nf.test" linenums="1" hl_lines="14-33"
    test("Should run successfully with correct number of processes") {

        when {
            params {
                input_file = "${projectDir}/greetings.csv"
            }
        }

        then {
            assert workflow.success
            assert workflow.trace.tasks().size() == 6
        }

    }

    test("Should produce correct output files") {

        when {
            params {
                input_file = "${projectDir}/greetings.csv"
            }
        }

        then {
            assert file("$launchDir/results/Bonjour-output.txt").exists()
            assert file("$launchDir/results/Hello-output.txt").exists()
            assert file("$launchDir/results/Holà-output.txt").exists()
            assert file("$launchDir/results/UPPER-Bonjour-output.txt").exists()
            assert file("$launchDir/results/UPPER-Hello-output.txt").exists()
            assert file("$launchDir/results/UPPER-Holà-output.txt").exists()
        }

    }
```

Ejecuta la prueba nuevamente para ver si funciona.

```bash title="nf-test pipeline pass"
nf-test test tests/main.nf.test
```

```console title="Pipeline passes with file assertions"
> nf-test test tests/main.nf.test

🚀 nf-test 0.9.3
https://www.nf-test.com
(c) 2021 - 2024 Lukas Forer and Sebastian Schoenherr


Test Workflow main.nf

  Test [f0e08a68] 'Should run successfully with correct number of processes' PASSED (8.144s)
  Test [d7e32a32] 'Should produce correct output files' PASSED (6.994s)


SUCCESS: Executed 2 tests in 15.165s
```

¡Éxito! Las pruebas pasan porque el pipeline se completó exitosamente, se ejecutó el número correcto de procesos y se crearon los archivos de salida. Esto también debería mostrarte lo útil que es proporcionar esos nombres informativos para tus pruebas.

Esto es solo la superficie, podemos seguir escribiendo aserciones para verificar los detalles del pipeline, pero por ahora pasemos a probar los componentes internos del pipeline.

### Conclusión

Sabes cómo escribir un nf-test para un pipeline.

### ¿Qué sigue?

Aprende cómo probar un proceso de Nextflow.

---

## 2. Probar un proceso de Nextflow

No tenemos que escribir pruebas para cada parte del pipeline, pero cuantas más pruebas tengamos, más completos podemos ser sobre el pipeline y más confiados podemos estar de que está funcionando como se espera. En esta sección vamos a probar ambos procesos en el pipeline como unidades individuales.

### 2.1. Probar el proceso `sayHello`

Comencemos con el proceso `sayHello`.

Usemos el comando `nf-test generate` nuevamente para generar pruebas para el proceso.

```bash
nf-test generate process main.nf
```

```console title="Output"
> nf-test generate process main.nf

Load source file '/workspaces/training/side-quests/nf-test/main.nf'
Wrote process test file '/workspaces/training/side-quests/nf-test/tests/main.sayhello.nf.test
Wrote process test file '/workspaces/training/side-quests/nf-test/tests/main.converttoupper.nf.test

SUCCESS: Generated 2 test files.
```

Enfoquémonos por ahora en el proceso `sayhello` en el archivo `main.sayhello.nf.test`.

Abramos el archivo y echemos un vistazo al contenido.

```groovy title="tests/main.sayhello.nf.test"
nextflow_process {

    name "Test Process sayHello"
    script "main.nf"
    process "sayHello"

    test("Should run without failures") {

        when {
            params {
                // define parameters here. Example:
                // outdir = "tests/results"
            }
            process {
                """
                // define inputs of the process here. Example:
                // input[0] = file("test-file.txt")
                """
            }
        }

        then {
            assert process.success
            assert snapshot(process.out).match()
        }

    }

}
```

Como antes, comenzamos con los detalles de la prueba, seguidos de los bloques `when` y `then`. Sin embargo, también tenemos un bloque `process` adicional que nos permite definir las entradas al proceso.

Ejecutemos la prueba para ver si funciona.

```bash title="nf-test pipeline pass"
nf-test test tests/main.sayhello.nf.test
```

```console title="Process test fails"
> nf-test test tests/main.sayhello.nf.test

🚀 nf-test 0.9.3
https://www.nf-test.com
(c) 2021 - 2024 Lukas Forer and Sebastian Schoenherr


Test Process sayHello

  Test [1eaad118] 'Should run without failures' FAILED (4.876s)

  Assertion failed:

  assert process.success
         |       |
         |       false
         sayHello

  Nextflow stdout:

  Process `sayHello` declares 1 input but was called with 0 arguments
  Nextflow stderr:

FAILURE: Executed 1 tests in 4.884s (1 failed)
```

La prueba falla porque el proceso `sayHello` declara 1 entrada pero fue llamado con 0 argumentos. Arreglemos eso agregando una entrada al proceso. Recuerda de [Hello Workflow](../hello_nextflow/03_hello_workflow.md) (y la sección de calentamiento anterior) que nuestro proceso `sayHello` toma una única entrada de valor, que necesitaremos proporcionar. También deberíamos corregir el nombre de la prueba para reflejar mejor lo que estamos probando.

**Antes:**

```groovy title="tests/main.sayhello.nf.test" linenums="1" hl_lines="1 10 11"
    test("Should run without failures") {

        when {
            params {
                // define parameters here. Example:
                // outdir = "tests/results"
            }
            process {
                """
                // define inputs of the process here. Example:
                // input[0] = file("test-file.txt")
                """
            }
        }

        then {
            assert process.success
            assert snapshot(process.out).match()
        }

    }
```

**Después:**

```groovy title="tests/main.sayhello.nf.test" linenums="1" hl_lines="1 10"
    test("Should run without failures and produce correct output") {

        when {
            params {
                // define parameters here. Example:
                // outdir = "tests/results"
            }
            process {
                """
                input[0] = "hello"
                """
            }
        }

        then {
            assert process.success
            assert snapshot(process.out).match()
        }

    }
```

Ejecutemos la prueba nuevamente para ver si funciona.

```console title="nf-test pipeline pass"
> nf-test test tests/main.sayhello.nf.test

🚀 nf-test 0.9.3
https://www.nf-test.com
(c) 2021 - 2024 Lukas Forer and Sebastian Schoenherr


Test Process sayHello

  Test [f91a1bcd] 'Should run without failures and produce correct output' PASSED (1.604s)
  Snapshots:
    1 created [Should run without failures and produce correct output]


Snapshot Summary:
  1 created

SUCCESS: Executed 1 tests in 1.611s
```

¡Éxito! La prueba pasa porque el proceso `sayHello` se ejecutó exitosamente y se creó la salida.

### 2.2. Revisar el snapshot creado por la prueba

Si miramos el archivo `tests/main.sayhello.nf.test`, podemos ver que usa un método `snapshot()` en el bloque de aserción:

```groovy title="tests/main.sayhello.nf.test"
assert snapshot(process.out).match()
```

Esto le está diciendo a nf-test que cree un snapshot de la salida del proceso `sayHello`. Echemos un vistazo al contenido del archivo de snapshot.

```console title="Snapshot file contents"
code tests/main.sayhello.nf.test.snap
```

No lo imprimiremos aquí, pero deberías ver un archivo JSON que contiene detalles del proceso y las salidas del proceso. En particular, podemos ver una línea que se ve así:

```json title="Snapshot file contents"
"0": [
    "hello-output.txt:md5,b1946ac92492d2347c6235b4d2611184"
]
```

Esto representa las salidas creadas por el proceso `sayHello`, que estamos probando explícitamente. Si volvemos a ejecutar la prueba, el programa verificará que la nueva salida coincida con la salida que se registró originalmente. Esta es una forma rápida y simple de probar que las salidas del proceso no cambian, por lo que nf-test la proporciona como predeterminada.

!!!warning "Advertencia"

    ¡Eso significa que debemos estar seguros de que la salida que registramos en la ejecución original es correcta!

Si, en el curso del desarrollo futuro, algo en el código cambia que hace que la salida sea diferente, la prueba fallará y tendremos que determinar si el cambio es esperado o no.

- Si resulta que algo en el código se rompió, tendremos que arreglarlo, con la expectativa de que el código corregido pasará la prueba.
- Si es un cambio esperado (por ejemplo, la herramienta ha sido mejorada y los resultados son mejores) entonces necesitaremos actualizar el snapshot para aceptar la nueva salida como la referencia a coincidir. nf-test tiene un parámetro `--update-snapshot` para este propósito.

Podemos ejecutar la prueba nuevamente y ver que la prueba debería pasar:

```console title="nf-test process pass with snapshot"
> nf-test test tests/main.sayhello.nf.test

🚀 nf-test 0.9.3
https://www.nf-test.com
(c) 2021 - 2024 Lukas Forer and Sebastian Schoenherr


Test Process sayHello

  Test [f91a1bcd] 'Should run without failures and produce correct output' PASSED (1.675s)


SUCCESS: Executed 1 tests in 1.685s
```

¡Éxito! La prueba pasa porque el proceso `sayHello` se ejecutó exitosamente y la salida coincidió con el snapshot.

### 2.3. Alternativa a los Snapshots: Aserciones de Contenido Directo

Si bien los snapshots son excelentes para detectar cualquier cambio en la salida, a veces quieres verificar contenido específico sin ser tan estricto sobre que todo el archivo coincida. Por ejemplo:

- Cuando partes de la salida pueden cambiar (marcas de tiempo, IDs aleatorios, etc.) pero cierto contenido clave debe estar presente
- Cuando quieres verificar patrones o valores específicos en la salida
- Cuando quieres hacer la prueba más explícita sobre qué constituye el éxito

Aquí está cómo podríamos modificar nuestra prueba para verificar contenido específico:

**Antes:**

```groovy title="tests/main.sayhello.nf.test" linenums="1" hl_lines="1 5 6 17"
    test("Should run without failures and produce correct output") {

        when {
            params {
                // define parameters here. Example:
                // outdir = "tests/results"
            }
            process {
                """
                input[0] = "hello"
                """
            }
        }

        then {
            assert process.success
            assert snapshot(process.out).match()
        }

    }
```

**Después:**

```groovy title="tests/main.sayhello.nf.test" linenums="1" hl_lines="1 5 16 17"
    test("Should run without failures and contain expected greeting") {

        when {
            params {
                // define parameters here
            }
            process {
                """
                input[0] = "hello"
                """
            }
        }

        then {
            assert process.success
            assert path(process.out[0][0]).readLines().contains('hello')
            assert !path(process.out[0][0]).readLines().contains('HELLO')
        }

    }
```

Nota que nf-test ve las salidas del proceso como una lista de listas, por lo que `process.out[0][0]` está obteniendo la primera parte del primer elemento del canal (o 'emisión') de este proceso.

Este enfoque:

- Deja claro exactamente qué esperamos en la salida
- Es más resistente a cambios irrelevantes en la salida
- Proporciona mejores mensajes de error cuando las pruebas fallan
- Permite validaciones más complejas (patrones regex, comparaciones numéricas, etc.)

Ejecutemos la prueba para ver si funciona.

```bash title="nf-test pipeline pass"
nf-test test tests/main.sayhello.nf.test
```

```console title="Process test fails"
> nf-test test tests/main.sayhello.nf.test

🚀 nf-test 0.9.3
https://www.nf-test.com
(c) 2021 - 2024 Lukas Forer and Sebastian Schoenherr


Test Process sayHello

  Test [58df4e4b] 'Should run without failures and contain expected greeting' PASSED (7.196s)


SUCCESS: Executed 1 tests in 7.208s
```

### 2.4. Probar el proceso `convertToUpper`

Abramos el archivo `tests/main.converttoupper.nf.test` y echemos un vistazo al contenido:

```groovy title="tests/main.converttoupper.nf.test"
nextflow_process {

    name "Test Process convertToUpper"
    script "main.nf"
    process "convertToUpper"

    test("Should run without failures") {

        when {
            params {
                // define parameters here. Example:
                // outdir = "tests/results"
            }
            process {
                """
                // define inputs of the process here. Example:
                // input[0] = file("test-file.txt")
                """
            }
        }

        then {
            assert process.success
            assert snapshot(process.out).match()
        }

    }

}
```

Esta es una prueba similar al proceso `sayHello`, pero está probando el proceso `convertToUpper`. Sabemos que esta fallará porque al igual que con `sayHello`, el proceso `convertToUpper` toma una única entrada de ruta, pero no hemos especificado una.

Ahora necesitamos proporcionar un único archivo de entrada al proceso convertToUpper, que incluya algo de texto que queremos convertir a mayúsculas. Hay muchas formas en que podríamos hacer esto:

- Podríamos crear un archivo dedicado para probar
- Podríamos reutilizar el archivo existente data/greetings.csv
- Podríamos crearlo sobre la marcha dentro de la prueba

Por ahora, reutilicemos el archivo existente data/greetings.csv usando el ejemplo que usamos con la prueba a nivel de pipeline. Como antes, podemos nombrar la prueba para reflejar mejor lo que estamos probando, pero esta vez dejemos que haga un 'snapshot' del contenido en lugar de verificar cadenas específicas (como hicimos en el otro proceso).

**Antes:**

```groovy title="tests/main.converttoupper.nf.test" linenums="1" hl_lines="1 10 11"
    test("Should run without failures") {

        when {
            params {
                // define parameters here. Example:
                // outdir = "tests/results"
            }
            process {
                """
                // define inputs of the process here. Example:
                // input[0] = file("test-file.txt")
                """
            }
        }

        then {
            assert process.success
            assert snapshot(process.out).match()
        }

    }
```

**Después:**

```groovy title="tests/main.converttoupper.nf.test" linenums="1" hl_lines="1 10"
    test("Should run without failures and produce correct output") {

        when {
            params {
                // define parameters here. Example:
                // outdir = "tests/results"
            }
            process {
                """
                input[0] = "${projectDir}/greetings.csv"
                """
            }
        }

        then {
            assert process.success
            assert snapshot(process.out).match()
        }

    }
```

¡Y ejecuta la prueba!

```bash title="nf-test pipeline pass"
nf-test test tests/main.converttoupper.nf.test
```

```console title="nf-test process convertToUpper pass"
> nf-test test tests/main.converttoupper.nf.test

🚀 nf-test 0.9.3
https://www.nf-test.com
(c) 2021 - 2024 Lukas Forer and Sebastian Schoenherr


Test Process convertToUpper

  Test [c59b6044] 'Should run without failures and produce correct output' PASSED (1.755s)
  Snapshots:
    1 created [Should run without failures and produce correct output]


Snapshot Summary:
  1 created

SUCCESS: Executed 1 tests in 1.764s
```

Nota, hemos creado un archivo de snapshot para el proceso `convertToUpper` en `tests/main.converttoupper.nf.test.snap`. Si ejecutamos la prueba nuevamente, deberíamos ver que nf-test pasa nuevamente.

```bash title="nf-test process convertToUpper pass"
nf-test test tests/main.converttoupper.nf.test
```

```console title="nf-test process convertToUpper pass"
> nf-test test tests/main.converttoupper.nf.test

🚀 nf-test 0.9.3
https://www.nf-test.com
(c) 2021 - 2024 Lukas Forer and Sebastian Schoenherr


Test Process convertToUpper

  Test [c59b6044] 'Should run without failures and produce correct output' PASSED (1.798s)


SUCCESS: Executed 1 tests in 1.811s
```

### Conclusión

Sabes cómo escribir pruebas para un proceso de Nextflow y ejecutarlas.

### ¿Qué sigue?

¡Aprende cómo ejecutar pruebas para todo a la vez!

## 3. Ejecutar pruebas para todo el repositorio

Ejecutar nf-test en cada componente está bien, pero es laborioso y propenso a errores. ¿No podemos simplemente probar todo a la vez?

¡Sí podemos!

Ejecutemos nf-test en todo el repositorio.

### 3.1. Ejecutar nf-test en todo el repositorio

Podemos ejecutar nf-test en todo el repositorio ejecutando el comando `nf-test test`.

```bash
nf-test test .
```

Nota, solo estamos usando el `.` para ejecutar todo desde nuestro directorio actual. ¡Esto incluirá cada prueba!

```console title="nf-test repo pass"
> nf-test test .

🚀 nf-test 0.9.3
https://www.nf-test.com
(c) 2021 - 2024 Lukas Forer and Sebastian Schoenherr


Test Process convertToUpper

  Test [3d26d9af] 'Should run without failures and produce correct output' PASSED (4.155s)

Test Workflow main.nf

  Test [f183df37] 'Should run successfully with correct number of processes' PASSED (3.33s)
  Test [d7e32a32] 'Should produce correct output files' PASSED (3.102s)

Test Process sayHello

  Test [58df4e4b] 'Should run without failures and contain expected greeting' PASSED (2.614s)


SUCCESS: Executed 4 tests in 13.481s
```

¡Mira eso! Ejecutamos 4 pruebas, 1 para cada proceso y 2 para todo el pipeline con un solo comando. ¡Imagina lo poderoso que es esto en una base de código grande!

---

## Resumen

En esta misión secundaria, has aprendido a aprovechar las características de nf-test para crear y ejecutar pruebas para procesos individuales, así como pruebas de extremo a extremo para todo el pipeline.
Ahora conoces los dos enfoques principales para la validación de salidas, snapshots y aserciones de contenido directo, y cuándo usar cada uno.
También sabes cómo ejecutar pruebas una por una o para un proyecto completo.

Aplicar estas técnicas en tu propio trabajo te permitirá asegurar que:

- Tu código funciona como se espera
- Los cambios no rompen la funcionalidad existente
- Otros desarrolladores pueden contribuir con confianza
- Los problemas pueden identificarse y corregirse rápidamente
- El contenido de salida coincide con las expectativas

### Patrones clave

<!-- TODO: Can we add snippets of code below to illustrate? -->

1. Pruebas a nivel de pipeline:
   - Pruebas básicas de éxito
   - Verificación del conteo de procesos
   - Verificaciones de existencia de archivos de salida
2. Pruebas a nivel de proceso
3. Dos enfoques para la validación de salidas:
   - Usar snapshots para verificación completa de salida
   - Usar aserciones de contenido directo para verificaciones de contenido específico
4. Ejecutar todas las pruebas en un repositorio con un solo comando

### Recursos adicionales

Consulta la [documentación de nf-test](https://www.nf-test.com/) para características de prueba más avanzadas y mejores prácticas. Puede que quieras:

- Agregar aserciones más completas a tus pruebas
- Escribir pruebas para casos extremos y condiciones de error
- Configurar integración continua para ejecutar pruebas automáticamente
- Aprender sobre otros tipos de pruebas como pruebas de workflow y módulo
- Explorar técnicas de validación de contenido más avanzadas

**Recuerda:** Las pruebas son documentación viva de cómo debe comportarse tu código. Cuantas más pruebas escribas, y más específicas sean tus aserciones, más confiado puedes estar en la confiabilidad de tu pipeline.

---

## ¿Qué sigue?

Regresa al [menú de Misiones Secundarias](./index.md) o haz clic en el botón en la parte inferior derecha de la página para pasar al siguiente tema de la lista.
