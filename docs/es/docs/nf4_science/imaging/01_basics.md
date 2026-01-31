# Parte 1: Ejecutar operaciones básicas

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Traducción asistida por IA - [más información y sugerencias](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

En esta primera parte del curso de entrenamiento de Nextflow para Bioimagen, utilizaremos un ejemplo muy básico de Hello World independiente del dominio para demostrar las operaciones esenciales y señalar los componentes de código de Nextflow correspondientes.

## 1. Ejecutar el workflow

Le proporcionamos un script de workflow llamado `hello-world.nf` que recibe una entrada mediante un argumento de línea de comandos llamado `--greeting` y produce un archivo de texto que contiene ese saludo.
Todavía no vamos a ver el código; primero veamos cómo se ve ejecutarlo.

### 1.1. Lanzar el workflow y monitorear la ejecución

En la terminal, ejecute el siguiente comando:

```bash
nextflow run hello-world.nf --greeting 'Hello World!'
```

La salida de su consola debería verse algo así:

```console title="Output" linenums="1"
 N E X T F L O W   ~  version 25.04.3

Launching `hello-world.nf` [goofy_torvalds] DSL2 - revision: c33d41f479

executor >  local (1)
[a3/7be2fa] sayHello | 1 of 1 ✔
```

¡Felicidades, acaba de ejecutar su primer workflow de Nextflow!

La salida más importante aquí es la última línea (línea 6):

```console title="Output" linenums="6"
[a3/7be2fa] sayHello | 1 of 1 ✔
```

Esto nos dice que el proceso `sayHello` fue ejecutado exitosamente una vez (`1 of 1 ✔`).

Eso es genial, pero puede estar preguntándose: ¿dónde está la salida?

### 1.2. Encontrar el archivo de salida en el directorio `results`

Este workflow está configurado para publicar su salida en un directorio llamado `results`.
Si observa su directorio actual, verá que cuando ejecutó el workflow, Nextflow creó un nuevo directorio llamado `results`, que contiene un archivo llamado `output.txt`.

```console title="results/" linenums="1"
results
└── output.txt
```

Abra el archivo; el contenido debería coincidir con el saludo que especificó en la línea de comandos.

<details>
  <summary>Contenido del archivo</summary>

```console title="results/output.txt" linenums="1"
Hello World!
```

</details>

¡Eso es genial, nuestro workflow hizo lo que se suponía que debía hacer!

Sin embargo, tenga en cuenta que el resultado 'publicado' es una copia (o en algunos casos un enlace simbólico) de la salida real producida por Nextflow cuando ejecutó el workflow.

Así que ahora, vamos a mirar bajo el capó para ver dónde Nextflow realmente ejecutó el trabajo.

!!! warning "Advertencia"

    No todos los workflows estarán configurados para publicar las salidas en un directorio de resultados, y/o el nombre del directorio puede ser diferente.
    Un poco más adelante en esta sección, le mostraremos cómo averiguar dónde se especifica este comportamiento.

### 1.3. Encontrar la salida original y los logs en el directorio `work/`

Cuando ejecuta un workflow, Nextflow crea un 'directorio de tarea' distinto para cada invocación de cada proceso en el workflow (=cada paso en el pipeline).
Para cada uno, preparará las entradas necesarias, ejecutará la(s) instrucción(es) relevante(s) y escribirá las salidas y archivos de log dentro de ese único directorio, que se nombra automáticamente usando un hash para hacerlo único.

Todos estos directorios de tarea vivirán bajo un directorio llamado `work` dentro de su directorio actual (donde está ejecutando el comando).

Eso puede sonar confuso, así que veamos cómo se ve eso en la práctica.

Volviendo a la salida de la consola para el workflow que ejecutamos anteriormente, teníamos esta línea:

```console title="Extracto de la salida del comando" linenums="6"
[a3/7be2fa] sayHello | 1 of 1 ✔
```

¿Ve cómo la línea comienza con `[a3/7be2fa]`?
Esa es una forma truncada de la ruta del directorio de tarea para esa llamada de proceso, y le dice dónde encontrar la salida de la llamada al proceso `sayHello` dentro de la ruta del directorio `work/`.

Puede encontrar la ruta completa escribiendo el siguiente comando (reemplazando `a3/7be2fa` con lo que ve en su propia terminal) y presionando la tecla tab para autocompletar la ruta o agregando un asterisco:

```bash
tree work/a3/7be2fa*
```

Esto debería producir la ruta completa del directorio: `work/a3/7be2fa7be2fad5e71e5f49998f795677fd68`

Echemos un vistazo a lo que hay ahí.

!!! Tip "Consejo"

    Si navega por el contenido del subdirectorio de tarea en el explorador de archivos de VSCode, verá todos los archivos de inmediato.
    Sin embargo, los archivos de log están configurados para ser invisibles en la terminal, así que si desea usar `ls` o `tree` para verlos, deberá establecer la opción relevante para mostrar archivos invisibles.

    ```bash
    tree -a work
    ```

Los nombres exactos de los subdirectorios serán diferentes en su sistema.

<details>
  <summary>Contenido del directorio</summary>

```console title="work/"
work
└── a3
    └── 7be2fad5e71e5f49998f795677fd68
        ├── .command.begin
        ├── .command.err
        ├── .command.log
        ├── .command.out
        ├── .command.run
        ├── .command.sh
        ├── .exitcode
        └── output.txt
```

</details>

Debería reconocer inmediatamente el archivo `output.txt`, que de hecho es la salida original del proceso `sayHello` que fue publicada en el directorio `results`.
Si lo abre, encontrará el saludo `Hello World!` nuevamente.

<details>
  <summary>Contenido del archivo output.txt</summary>

```console title="work/a3/7be2fa7be2fad5e71e5f49998f795677fd68/output.txt" linenums="1"
Hello World!
```

</details>

Entonces, ¿qué hay de todos esos otros archivos?

Estos son los archivos auxiliares y de log que Nextflow escribió como parte de la ejecución de la tarea:

- **`.command.begin`**: Archivo centinela creado tan pronto como se lanza la tarea.
- **`.command.err`**: Mensajes de error (`stderr`) emitidos por la llamada al proceso
- **`.command.log`**: Salida de log completa emitida por la llamada al proceso
- **`.command.out`**: Salida regular (`stdout`) de la llamada al proceso
- **`.command.run`**: Script completo ejecutado por Nextflow para ejecutar la llamada al proceso
- **`.command.sh`**: El comando que fue realmente ejecutado por la llamada al proceso
- **`.exitcode`**: El código de salida resultante del comando

El archivo `.command.sh` es especialmente útil porque muestra el comando principal que Nextflow ejecutó sin incluir toda la contabilidad y configuración de tarea/entorno.

<details>
  <summary>Contenido del archivo</summary>

```console title="work/a3/7be2fa7be2fad5e71e5f49998f795677fd68/.command.sh" linenums="1"
#!/bin/bash -ue
echo 'Hello World!' > output.txt

```

</details>

!!! Tip "Consejo"

    Cuando algo sale mal y necesita solucionar problemas sobre qué sucedió, puede ser útil mirar el script `command.sh` para verificar exactamente qué comando compuso Nextflow basándose en las instrucciones del workflow, interpolación de variables y demás.

### 1.4. Ejercicio opcional: volver a ejecutar con diferentes saludos

Intente volver a ejecutar el workflow varias veces con diferentes valores para el argumento `--greeting`, luego observe tanto el contenido del directorio `results/` como los directorios de tarea.

Observe cómo se preservan las salidas y logs de los directorios de tarea aislados, mientras que el contenido del directorio `results` es sobrescrito por la salida de ejecuciones posteriores.

### Conclusión

Usted sabe cómo ejecutar un script simple de Nextflow, monitorear su ejecución y encontrar sus salidas.

### ¿Qué sigue?

Aprenda a leer un script básico de Nextflow e identificar cómo sus componentes se relacionan con su funcionalidad.

---

## 2. Examinar el script inicial del workflow Hello World

Lo que hicimos allí fue básicamente tratar el script del workflow como una caja negra.
Ahora que hemos visto qué hace, abramos la caja y miremos dentro.

_El objetivo aquí no es memorizar la sintaxis del código de Nextflow, sino formar algo de intuición básica sobre cuáles son los componentes principales y cómo están organizados._

### 2.1. Examinar la estructura general del código

Abramos el script `hello-world.nf` en el panel del editor.

<details>
  <summary>Código</summary>

```groovy title="hello-world.nf" linenums="1"
#!/usr/bin/env nextflow

/*
 * Usar echo para imprimir un saludo a un archivo
 */
process sayHello {

    publishDir 'results', mode: 'copy'

    input:
    val greeting

    output:
    path 'output.txt'

    script:
    """
    echo '$greeting' > output.txt
    """
}

workflow {

    // emitir un saludo
    sayHello(params.greeting)
}
```

</details>

Un script de Nextflow involucra dos tipos principales de componentes centrales: uno o más **processes**, y el **workflow** en sí.
Cada **process** describe qué operación(es) debe realizar el paso correspondiente en el pipeline, mientras que el **workflow** describe la lógica de flujo de datos que conecta los diversos pasos.

Veamos más de cerca el bloque **process** primero, luego veremos el bloque **workflow**.

### 2.2. La definición del `process`

El primer bloque de código describe un **process**.
La definición del proceso comienza con la palabra clave `process`, seguida del nombre del proceso y finalmente el cuerpo del proceso delimitado por llaves.
El cuerpo del proceso debe contener un bloque script que especifica el comando a ejecutar, que puede ser cualquier cosa que pueda ejecutar en una terminal de línea de comandos.

Aquí tenemos un **process** llamado `sayHello` que toma una variable de **input** llamada `greeting` y escribe su **output** en un archivo llamado `output.txt`.

<details>
  <summary>Código</summary>

```groovy title="hello-world.nf" linenums="3"
/*
 * Usar echo para imprimir un saludo a un archivo
 */
process sayHello {

    publishDir 'results', mode: 'copy'

    input:
    val greeting

    output:
    path 'output.txt'

    script:
    """
    echo '$greeting' > output.txt
    """
}
```

</details>

Esta es una definición de proceso muy mínima que solo contiene una definición de `input`, una definición de `output` y el `script` a ejecutar.

La definición de `input` incluye el calificador `val`, que le dice a Nextflow que espere un valor de algún tipo (puede ser una cadena, un número, lo que sea).

La definición de `output` incluye el calificador `path`, que le dice a Nextflow que esto debe manejarse como una ruta (incluye tanto rutas de directorio como archivos).

!!! Tip "Consejo"

    La definición de salida no _determina_ qué salida se creará.
    Simplemente _declara_ dónde encontrar el(los) archivo(s) de salida esperado(s), para que Nextflow pueda buscarlo una vez que la ejecución esté completa.

    Esto es necesario para verificar que el comando se ejecutó exitosamente y para pasar la salida a los procesos posteriores si es necesario.
    La salida producida que no coincida con lo declarado en el bloque de salida no se pasará a los procesos posteriores.

En un pipeline del mundo real, un proceso generalmente contiene información adicional como directivas de proceso, que presentaremos en un momento.

### 2.3. La definición del `workflow`

El segundo bloque de código describe el **workflow** en sí.
La definición del workflow comienza con la palabra clave `workflow`, seguida de un nombre opcional, luego el cuerpo del workflow delimitado por llaves.

Aquí tenemos un **workflow** que consiste en una llamada al proceso `sayHello`, que toma una entrada, `params.greeting`, que contiene el valor que dimos al parámetro `--greeting`.

```groovy title="hello-world.nf" linenums="22"
workflow {

    // emitir un saludo
    sayHello(params.greeting)
}
```

Esta es una definición de **workflow** muy mínima.
En un pipeline del mundo real, el workflow típicamente contiene múltiples llamadas a **processes** conectados por **channels**, y puede haber valores predeterminados configurados para las entradas de variables.

Veremos esto en acción cuando ejecutemos nf-core/molkart en la Parte 2 del curso.

### 2.4. El sistema `params` de parámetros de línea de comandos

El `params.greeting` que proporcionamos a la llamada del proceso `sayHello()` es un fragmento ingenioso de código de Nextflow y vale la pena dedicarle un minuto extra.

Como se mencionó anteriormente, así es como pasamos el valor del parámetro de línea de comandos `--greeting` a la llamada del proceso `sayHello()`.
De hecho, simplemente declarar `params.someParameterName` nos permitirá dar al workflow un parámetro llamado `--someParameterName` desde la línea de comandos.

!!! Tip "Consejo"

    Estos parámetros de workflow declarados usando el sistema `params` siempre toman dos guiones (`--`).
    Esto los distingue de los parámetros a nivel de Nextflow, que solo toman un guion (`-`).

### Conclusión

Ahora sabe cómo está estructurado un workflow simple de Nextflow, y cómo los componentes básicos se relacionan con su funcionalidad.

### ¿Qué sigue?

Aprenda a gestionar sus ejecuciones de workflow convenientemente.

---

## 3. Gestionar ejecuciones de workflow

Saber cómo lanzar workflows y recuperar salidas es genial, pero rápidamente encontrará que hay algunos otros aspectos de la gestión de workflows que harán su vida más fácil.

Aquí le mostramos cómo aprovechar la función `resume` para cuando necesite relanzar el mismo workflow, cómo inspeccionar los logs de ejecución con `nextflow log`, y cómo eliminar directorios de trabajo más antiguos con `nextflow clean`.

### 3.1. Relanzar un workflow con `-resume`

A veces, va a querer volver a ejecutar un pipeline que ya lanzó previamente sin rehacer ningún trabajo que ya se haya completado exitosamente.

Nextflow tiene una opción llamada `-resume` que le permite hacer esto.
Específicamente, en este modo, cualquier proceso que ya se haya ejecutado con exactamente el mismo código, configuraciones y entradas será omitido.
Esto significa que Nextflow solo ejecutará los procesos que haya agregado o modificado desde la última ejecución, o a los que esté proporcionando nuevas configuraciones o entradas.

Hay dos ventajas clave al hacer esto:

- Si está en medio del desarrollo de un pipeline, puede iterar más rápidamente ya que solo tiene que ejecutar el(los) proceso(s) en los que está trabajando activamente para probar sus cambios.
- Si está ejecutando un pipeline en producción y algo sale mal, en muchos casos puede corregir el problema y relanzar el pipeline, y se reanudará ejecutándose desde el punto de falla, lo que puede ahorrarle mucho tiempo y cómputo.

Para usarlo, simplemente agregue `-resume` a su comando y ejecútelo:

```bash
nextflow run hello-world.nf --greeting 'Hello World!' -resume
```

??? success "Salida del comando"

    ```console
    N E X T F L O W   ~  version 25.04.3

    Launching `hello-world.nf` [tiny_noyce] DSL2 - revision: c33d41f479

    [a3/7be2fa] process > sayHello [100%] 1 of 1, cached: 1 ✔
    ```

Busque el bit `cached:` que se ha agregado en la línea de estado del proceso (línea 5), lo que significa que Nextflow ha reconocido que ya ha hecho este trabajo y simplemente reutilizó el resultado de la ejecución exitosa anterior.

También puede ver que el hash del subdirectorio de trabajo es el mismo que en la ejecución anterior.
Nextflow literalmente le está señalando la ejecución anterior y diciendo "Ya hice eso allí".

!!! Tip "Consejo"

    Cuando vuelve a ejecutar un pipeline con `resume`, Nextflow no sobrescribe ningún archivo escrito en un directorio `publishDir` por ninguna llamada de proceso que se haya ejecutado previamente de manera exitosa.

### 3.2. Inspeccionar el log de ejecuciones pasadas

Cada vez que lanza un workflow de nextflow, se escribe una línea en un archivo de log llamado `history`, bajo un directorio oculto llamado `.nextflow` en el directorio de trabajo actual.

Una forma más conveniente de acceder a esta información es usar el comando `nextflow log`.

```bash
nextflow log
```

Esto mostrará el contenido del archivo de log en la terminal, mostrándole la marca de tiempo, nombre de ejecución, estado y línea de comandos completa para cada ejecución de Nextflow que se haya lanzado desde dentro del directorio de trabajo actual.

### 3.3. Eliminar directorios de trabajo más antiguos

Durante el proceso de desarrollo, típicamente ejecutará sus pipelines en borrador un gran número de veces, lo que puede llevar a una acumulación de muchos archivos en muchos subdirectorios.
Dado que los subdirectorios se nombran aleatoriamente, es difícil distinguir por sus nombres cuáles son ejecuciones más antiguas vs. más recientes.

Nextflow incluye un conveniente subcomando `clean` que puede eliminar automáticamente los subdirectorios de trabajo de ejecuciones pasadas que ya no le interesan, con varias [opciones](https://www.nextflow.io/docs/latest/reference/cli.html#clean) para controlar qué se eliminará.

Puede usar el log de Nextflow para buscar una ejecución basándose en su marca de tiempo y/o línea de comandos, luego usar `nextflow clean -before <run_name> -f` para eliminar directorios de trabajo de ejecuciones anteriores.

!!! Warning "Advertencia"

    Eliminar subdirectorios de trabajo de ejecuciones pasadas los elimina del caché de Nextflow y elimina cualquier salida que se haya almacenado en esos directorios.
    Esto significa que rompe la capacidad de Nextflow de reanudar la ejecución sin volver a ejecutar los procesos correspondientes.

    ¡Usted es responsable de guardar cualquier salida que le importe o en la que planee confiar! Si está usando la directiva `publishDir` para ese propósito, asegúrese de usar el modo `copy`, no el modo `symlink`.

### Conclusión

Usted sabe cómo relanzar un pipeline sin repetir pasos que ya se ejecutaron de manera idéntica, inspeccionar el log de ejecución, y usar el comando `nextflow clean` para limpiar directorios de trabajo antiguos.

### ¿Qué sigue?

Ahora que entiende las operaciones básicas de Nextflow, está listo para ejecutar un pipeline real de bioimagen con nf-core/molkart.
