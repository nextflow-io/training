# Parte 4: Hello Modules

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Traducción asistida por IA - [más información y sugerencias](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

<!--
<div class="video-wrapper">
  <iframe width="560" height="315" src="https://www.youtube.com/embed/Xxp_menS0E8?si=0AWnXB7xqHAzJdJV&amp;list=PLPZ8WHdZGxmXiHf8B26oB_fTfoKQdhlik" title="YouTube video player" frameborder="0" allow="accelerometer; autoplay; clipboard-write; encrypted-media; gyroscope; picture-in-picture; web-share" referrerpolicy="strict-origin-when-cross-origin" allowfullscreen></iframe>
</div>

/// caption
:fontawesome-brands-youtube:{ .youtube } Vea [la lista de reproducción completa](https://www.youtube.com/playlist?list=PLPZ8WHdZGxmXiHf8B26oB_fTfoKQdhlik) en el canal de YouTube de Nextflow.

:green_book: La transcripción del video está disponible [aquí](./transcripts/04_hello_modules.md).
///
-->

Esta sección cubre cómo organizar el código de su flujo de trabajo para hacer el desarrollo y mantenimiento de su pipeline más eficiente y sostenible.
Específicamente, vamos a demostrar cómo usar **módulos**.

En Nextflow, un **módulo** es una única definición de proceso que está encapsulada por sí misma en un archivo de código independiente.
Para usar un módulo en un flujo de trabajo, solo agrega una declaración de importación de una línea a su archivo de código de workflow; luego puede integrar el proceso en el flujo de trabajo de la misma manera que normalmente lo haría.
Eso hace posible reutilizar definiciones de proceso en múltiples flujos de trabajo sin producir múltiples copias del código.

Cuando comenzamos a desarrollar nuestro flujo de trabajo, escribimos todo en un único archivo de código.
Ahora vamos a mover los procesos a módulos individuales.

<figure class="excalidraw">
--8<-- "docs/en/docs/hello_nextflow/img/modules.svg"
</figure>

Esto hará nuestro código más compartible, flexible y mantenible.

??? info "Cómo comenzar desde esta sección"

    Esta sección del curso asume que ha completado las Partes 1-3 del curso [Hello Nextflow](./index.md), pero si se siente cómodo con los conceptos básicos cubiertos en esas secciones, puede comenzar desde aquí sin hacer nada especial.

---

## 0. Calentamiento: Ejecutar `hello-modules.nf`

Vamos a usar el script de workflow `hello-modules.nf` como punto de partida.
Es equivalente al script producido al trabajar en la Parte 3 de este curso de entrenamiento, excepto que hemos cambiado los destinos de salida:

```groovy title="hello-modules.nf" linenums="37" hl_lines="3 7 11 15"
output {
    first_output {
        path 'hello_modules'
        mode 'copy'
    }
    uppercased {
        path 'hello_modules'
        mode 'copy'
    }
    collected {
        path 'hello_modules'
        mode 'copy'
    }
    batch_report {
        path 'hello_modules'
        mode 'copy'
    }
}
```

Solo para asegurarse de que todo funciona, ejecute el script una vez antes de hacer cualquier cambio:

```bash
nextflow run hello-modules.nf
```

??? success "Salida del comando"

    ```console
     N E X T F L O W   ~  version 25.10.2

    Launching `hello-modules.nf` [hopeful_avogadro] DSL2 - revision: b09af1237d

    executor >  local (7)
    [0f/8795c9] sayHello (3)       [100%] 3 of 3 ✔
    [6a/eb2510] convertToUpper (3) [100%] 3 of 3 ✔
    [af/479117] collectGreetings   [100%] 1 of 1 ✔
    ```

Como anteriormente, encontrará los archivos de salida en el directorio especificado en el bloque `output` (aquí, `results/hello_modules/`).

??? abstract "Contenido del directorio"

    ```console
    results/hello_modules/
    ├── Bonjour-output.txt
    ├── COLLECTED-batch-output.txt
    ├── Hello-output.txt
    ├── Holà-output.txt
    ├── batch-report.txt
    ├── UPPER-Bonjour-output.txt
    ├── UPPER-Hello-output.txt
    └── UPPER-Holà-output.txt
    ```

Si eso funcionó para usted, está listo para aprender cómo modularizar el código de su flujo de trabajo.

---

## 1. Crear un directorio para almacenar módulos

Es una buena práctica almacenar sus módulos en un directorio específico.
Puede llamar a ese directorio como quiera, pero la convención es llamarlo `modules/`.

```bash
mkdir modules
```

!!! tip "Consejo"

    Aquí le estamos mostrando cómo usar **módulos locales**, es decir, módulos almacenados localmente en el mismo repositorio que el resto del código del flujo de trabajo, en contraste con módulos remotos, que se almacenan en otros repositorios (remotos).
    Para más información sobre **módulos remotos**, vea la [documentación](https://www.nextflow.io/docs/latest/module.html).

---

## 2. Crear un módulo para `sayHello()`

En su forma más simple, convertir un proceso existente en un módulo es poco más que una operación de copiar y pegar.
Vamos a crear un archivo stub para el módulo, copiar el código relevante y luego eliminarlo del archivo de flujo de trabajo principal.

Luego todo lo que necesitaremos hacer es agregar una declaración de importación para que Nextflow sepa que debe traer el código relevante en tiempo de ejecución.

### 2.1. Crear un archivo stub para el nuevo módulo

Creemos un archivo vacío para el módulo llamado `sayHello.nf`.

```bash
touch modules/sayHello.nf
```

Esto nos da un lugar para poner el código del proceso.

### 2.2. Mover el código del proceso `sayHello` al archivo del módulo

Copie toda la definición del proceso desde el archivo de workflow al archivo del módulo, asegurándose de copiar también el shebang `#!/usr/bin/env nextflow`.

```groovy title="modules/sayHello.nf" linenums="1"
#!/usr/bin/env nextflow

/*
 * Usar echo para imprimir 'Hello World!' a un archivo
 */
process sayHello {

    input:
    val greeting

    output:
    path "${greeting}-output.txt"

    script:
    """
    echo '${greeting}' > '${greeting}-output.txt'
    """
}
```

Una vez hecho eso, elimine la definición del proceso del archivo de workflow, pero asegúrese de dejar el shebang en su lugar.

### 2.3. Agregar una declaración de importación antes del bloque workflow

La sintaxis para importar un módulo local es bastante sencilla:

```groovy title="Sintaxis: Declaración de importación"
include { <MODULE_NAME> } from '<path_to_module>'
```

Insertemos eso arriba del bloque `params` y completémoslo apropiadamente.

=== "Después"

    ```groovy title="hello-modules.nf" linenums="44" hl_lines="1 2"
    // Incluir módulos
    include { sayHello } from './modules/sayHello.nf'

    /*
    * Parámetros del pipeline
    */
    params {
        greeting: Path = 'data/greetings.csv'
        batch: String = 'batch'
    }
    ```

=== "Antes"

    ```groovy title="hello-modules.nf" linenums="44"
    /*
    * Parámetros del pipeline
    */
    params {
        greeting: Path = 'data/greetings.csv'
        batch: String = 'batch'
    }
    ```

Verá que hemos completado el nombre del módulo, `sayHello`, y la ruta al archivo que contiene el código del módulo, `./modules/sayHello.nf`.

### 2.4. Ejecutar el flujo de trabajo

Estamos ejecutando el flujo de trabajo con esencialmente el mismo código y entradas que antes, así que ejecutemos con la bandera `-resume` y veamos qué sucede.

```bash
nextflow run hello-modules.nf -resume
```

??? success "Salida del comando"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `hello-modules.nf` [romantic_poisson] DSL2 - revision: 96edfa9ad3

    [f6/cc0107] sayHello (1)       | 3 of 3, cached: 3 ✔
    [3c/4058ba] convertToUpper (2) | 3 of 3, cached: 3 ✔
    [1a/bc5901] collectGreetings   | 1 of 1, cached: 1 ✔
    ```

Esto debería ejecutarse muy rápidamente porque todo está en caché.
Siéntase libre de verificar las salidas publicadas.

Nextflow reconoció que sigue siendo el mismo trabajo por hacer, incluso si el código está dividido en múltiples archivos.

### Conclusión

Sabe cómo extraer un proceso en un módulo local y sabe que hacer esto no rompe la capacidad de reanudar del flujo de trabajo.

### ¿Qué sigue?

Practicar haciendo más módulos.
Una vez que ha hecho uno, puede hacer un millón más...
Pero por ahora hagamos solo dos más.

---

## 3. Modularizar el proceso `convertToUpper()`

### 3.1. Crear un archivo stub para el nuevo módulo

Cree un archivo vacío para el módulo llamado `convertToUpper.nf`.

```bash
touch modules/convertToUpper.nf
```

### 3.2. Mover el código del proceso `convertToUpper` al archivo del módulo

Copie toda la definición del proceso desde el archivo de workflow al archivo del módulo, asegurándose de copiar también el shebang `#!/usr/bin/env nextflow`.

```groovy title="modules/convertToUpper.nf" linenums="1"
#!/usr/bin/env nextflow

/*
 * Usar una herramienta de reemplazo de texto para convertir el saludo a mayúsculas
 */
process convertToUpper {

    input:
    path input_file

    output:
    path "UPPER-${input_file}"

    script:
    """
    cat '${input_file}' | tr '[a-z]' '[A-Z]' > 'UPPER-${input_file}'
    """
}
```

Una vez hecho eso, elimine la definición del proceso del archivo de workflow, pero asegúrese de dejar el shebang en su lugar.

### 3.3. Agregar una declaración de importación antes del bloque `params`

Inserte la declaración de importación arriba del bloque `params` y complétela apropiadamente.

=== "Después"

    ```groovy title="hello-modules.nf" linenums="23" hl_lines="3"
    // Incluir módulos
    include { sayHello } from './modules/sayHello.nf'
    include { convertToUpper } from './modules/convertToUpper.nf'

    /*
    * Parámetros del pipeline
    */
    params {
        greeting: Path = 'data/greetings.csv'
        batch: String = 'batch'
    }
    ```

=== "Antes"

    ```groovy title="hello-modules.nf" linenums="23"
    // Incluir módulos
    include { sayHello } from './modules/sayHello.nf'

    /*
    * Parámetros del pipeline
    */
    params {
        greeting: Path = 'data/greetings.csv'
        batch: String = 'batch'
    }
    ```

Esto debería empezar a verse muy familiar.

### 3.4. Ejecutar el flujo de trabajo nuevamente

Ejecute esto con la bandera `-resume`.

```bash
nextflow run hello-modules.nf -resume
```

??? success "Salida del comando"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `hello-modules.nf` [nauseous_heisenberg] DSL2 - revision: a04a9f2da0

    [c9/763d42] sayHello (3)       | 3 of 3, cached: 3 ✔
    [60/bc6831] convertToUpper (3) | 3 of 3, cached: 3 ✔
    [1a/bc5901] collectGreetings   | 1 of 1, cached: 1 ✔
    ```

Esto debería seguir produciendo la misma salida que anteriormente.

¡Dos listos, uno más por hacer!

---

## 4. Modularizar el proceso `collectGreetings()`

### 4.1. Crear un archivo stub para el nuevo módulo

Cree un archivo vacío para el módulo llamado `collectGreetings.nf`.

```bash
touch modules/collectGreetings.nf
```

### 4.2. Mover el código del proceso `collectGreetings` al archivo del módulo

Copie toda la definición del proceso desde el archivo de workflow al archivo del módulo, asegurándose de copiar también el shebang `#!/usr/bin/env nextflow`.

```groovy title="modules/collectGreetings.nf" linenums="1"
#!/usr/bin/env nextflow

/*
 * Recopilar saludos en mayúsculas en un único archivo de salida
 */
process collectGreetings {

    input:
    path input_files
    val batch_name

    output:
    path "COLLECTED-${batch_name}-output.txt", emit: outfile
    path "${batch_name}-report.txt", emit: report

    script:
    count_greetings = input_files.size()
    """
    cat ${input_files} > 'COLLECTED-${batch_name}-output.txt'
    echo 'There were ${count_greetings} greetings in this batch.' > '${batch_name}-report.txt'
    """
}
```

Una vez hecho eso, elimine la definición del proceso del archivo de workflow, pero asegúrese de dejar el shebang en su lugar.

### 4.3. Agregar una declaración de importación antes del bloque `params`

Inserte la declaración de importación arriba del bloque `params` y complétela apropiadamente.

=== "Después"

    ```groovy title="hello-modules.nf" linenums="3" hl_lines="4"
    // Incluir módulos
    include { sayHello } from './modules/sayHello.nf'
    include { convertToUpper } from './modules/convertToUpper.nf'
    include { collectGreetings } from './modules/collectGreetings.nf'

    /*
    * Parámetros del pipeline
    */
    params {
        greeting: Path = 'data/greetings.csv'
        batch: String = 'batch'
    }
    ```

=== "Antes"

    ```groovy title="hello-modules.nf" linenums="3"
    // Incluir módulos
    include { sayHello } from './modules/sayHello.nf'
    include { convertToUpper } from './modules/convertToUpper.nf'

    /*
    * Parámetros del pipeline
    */
    params {
        greeting: Path = 'data/greetings.csv'
        batch: String = 'batch'
    }
    ```

¡El último!

### 4.4. Ejecutar el flujo de trabajo

Ejecute esto con la bandera `-resume`.

```bash
nextflow run hello-modules.nf -resume
```

??? success "Salida del comando"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `hello-modules.nf` [friendly_coulomb] DSL2 - revision: 7aa2b9bc0f

    [f6/cc0107] sayHello (1)       | 3 of 3, cached: 3 ✔
    [3c/4058ba] convertToUpper (2) | 3 of 3, cached: 3 ✔
    [1a/bc5901] collectGreetings   | 1 of 1, cached: 1 ✔
    ```

Esto debería seguir produciendo la misma salida que anteriormente.

### Conclusión

Sabe cómo modularizar múltiples procesos en un flujo de trabajo.

¡Felicitaciones, ha hecho todo este trabajo y absolutamente nada ha cambiado en cómo funciona el pipeline!

Bromas aparte, ahora su código es más modular, y si decide escribir otro pipeline que llame a uno de esos procesos, solo necesita escribir una corta declaración de importación para usar el módulo relevante.
Esto es mejor que copiar y pegar el código, porque si más tarde decide mejorar el módulo, todos sus pipelines heredarán las mejoras.

### ¿Qué sigue?

Tome un pequeño descanso si lo desea.

Cuando esté listo, continúe con [**Parte 5: Hola Containers**](./05_hello_containers.md) para aprender cómo usar contenedores para gestionar dependencias de software de manera más conveniente y reproducible.

---

## Cuestionario

<quiz>
¿Qué es un módulo en Nextflow?
- [ ] Un archivo de configuración
- [x] Un archivo independiente que contiene una única definición de proceso
- [ ] Una definición de workflow
- [ ] Un operador de channel

Aprenda más: [2. Crear un módulo para `sayHello()`](#2-crear-un-modulo-para-sayhello)
</quiz>

<quiz>
¿Cuál es la convención de nombres recomendada para archivos de módulos?
- [ ] `module_processName.nf`
- [ ] `processName_module.nf`
- [x] `processName.nf`
- [ ] `mod_processName.nf`
</quiz>

<quiz>
¿Dónde deberían almacenarse los archivos de módulos?
- [ ] En el mismo directorio que el workflow
- [ ] En un directorio `bin/`
- [x] En un directorio `modules/`
- [ ] En un directorio `lib/`

Aprenda más: [1. Crear un directorio para almacenar módulos](#1-crear-un-directorio-para-almacenar-modulos)
</quiz>

<quiz>
¿Cuál es la sintaxis correcta para importar un módulo?

- [ ] `#!groovy import { SAYHELLO } from './modules/sayhello.nf'`
- [ ] `#!groovy require { SAYHELLO } from './modules/sayhello.nf'`
- [x] `#!groovy include { SAYHELLO } from './modules/sayhello.nf'`
- [ ] `#!groovy load { SAYHELLO } from './modules/sayhello.nf'`

Aprenda más: [2.3. Agregar una declaración de importación](#23-agregar-una-declaracion-de-importacion-antes-del-bloque-workflow)
</quiz>

<quiz>
¿Qué sucede con la funcionalidad `-resume` cuando se usan módulos?
- [ ] Ya no funciona
- [ ] Requiere configuración adicional
- [x] Funciona igual que antes
- [ ] Solo funciona para módulos locales
</quiz>

<quiz>
¿Cuáles son los beneficios de usar módulos? (Seleccione todos los que apliquen)
- [x] Reutilización de código entre flujos de trabajo
- [x] Mantenimiento más fácil
- [x] Mejor organización del código del flujo de trabajo
- [ ] Velocidad de ejecución más rápida
</quiz>
