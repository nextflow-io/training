# Parte 3: Configuración de ejecución

Esta sección explorará cómo gestionar la configuración de un pipeline de Nextflow para personalizar su comportamiento, adaptarlo a diferentes entornos y optimizar el uso de recursos _sin alterar una sola línea del código del workflow_.

Hay múltiples formas de hacer esto, que pueden usarse en combinación y se interpretan según el orden de precedencia descrito en la documentación de [Configuration](https://nextflow.io/docs/latest/config.html).

En esta parte del curso, vamos a mostrarte el mecanismo de archivo de configuración más simple y común, el archivo `nextflow.config`, que ya encontraste en la sección sobre contenedores en la Parte 2.

Revisaremos componentes esenciales de la configuración de Nextflow como directivas de proceso, ejecutores, perfiles y archivos de parámetros.
Al aprender a utilizar estas opciones de configuración de manera efectiva, podrás aprovechar al máximo la flexibilidad, escalabilidad y rendimiento de los pipelines de Nextflow.

Para ejercitar estos elementos de configuración, vamos a ejecutar una copia nueva del workflow que ejecutamos por última vez al final de la Parte 2 de este curso de capacitación, renombrado `3-main.nf`.

Si no estás familiarizado con el pipeline Hello o necesitas un recordatorio, consulta [esta página de información](../info/hello_pipeline.md).

---

## 1. Gestionar parámetros de entrada del workflow

??? example "Escenario"

    Has descargado un pipeline y quieres ejecutarlo repetidamente con los mismos archivos de entrada y configuraciones, pero no quieres escribir todos los parámetros cada vez.
    O tal vez estás configurando el pipeline para un colega que no se siente cómodo con argumentos de línea de comandos.

Vamos a comenzar con un aspecto de la configuración que es simplemente una extensión de lo que hemos estado trabajando hasta ahora: la gestión de parámetros de entrada.

Actualmente, nuestro workflow está configurado para aceptar varios valores de parámetros a través de la línea de comandos, declarados en un bloque `params` en el propio script del workflow.
Uno tiene un valor predeterminado establecido como parte de su declaración.

Sin embargo, es posible que desees establecer valores predeterminados para todos ellos, o anular el predeterminado existente sin tener que especificar parámetros en la línea de comandos o modificar el archivo de script original.

Hay múltiples formas de hacer eso; vamos a mostrarte tres formas básicas que se usan muy comúnmente.

### 1.1. Configurar valores en `nextflow.config`

Este es el enfoque más simple, aunque posiblemente el menos flexible ya que el archivo principal `nextflow.config` no es algo que quieras estar editando para cada ejecución.
Pero tiene la ventaja de separar las preocupaciones de _declarar_ los parámetros en el workflow (que definitivamente pertenece allí) versus proporcionar _valores predeterminados_, que están más en casa en un archivo de configuración.

Hagamos esto en dos pasos.

#### 1.1.1. Crear un bloque `params` en el archivo de configuración

Realiza los siguientes cambios de código en el archivo `nextflow.config`:

=== "Después"

    ```groovy title="nextflow.config" linenums="1" hl_lines="3-10"
    docker.enabled = true

    /*
    * Pipeline parameters
    */
    params {
        input = 'data/greetings.csv'
        batch = 'batch'
        character = 'turkey'
    }
    ```

=== "Antes"

    ```groovy title="nextflow.config" linenums="1"
    docker.enabled = true
    ```

Nota que no simplemente copiamos el bloque `params` del workflow al archivo de configuración.
Para el parámetro `batch` que ya tenía un valor predeterminado declarado, la sintaxis es un poco diferente.
En el archivo del workflow, esa es una declaración tipada.
En la configuración, esas son asignaciones de valores.

Técnicamente, esto es suficiente para anular los valores predeterminados aún especificados en el archivo del workflow.
Podrías modificar el valor predeterminado para `batch` y ejecutar el workflow para verificar que el valor establecido en el archivo de configuración anula el establecido en el archivo del workflow.

Pero en el espíritu de mover la configuración completamente al archivo de configuración, eliminemos ese valor predeterminado del archivo del workflow por completo.

#### 1.1.2. Eliminar el valor predeterminado para `batch` en el archivo del workflow

Realiza el siguiente cambio de código en el archivo del workflow `3-main.nf`:

=== "Después"

    ```groovy title="3-main.nf" linenums="9" hl_lines="6"
    /*
    * Pipeline parameters
    */
    params {
        input: Path
        batch: String
        character: String
    }
    ```

=== "Antes"

    ```groovy title="3-main.nf" linenums="9" hl_lines="6"
    /*
    * Pipeline parameters
    */
    params {
        input: Path
        batch: String = 'batch'
        character: String
    }
    ```

Ahora el archivo del workflow en sí no establece ningún valor predeterminado para estos parámetros.

#### 1.1.3. Ejecutar el pipeline

Probemos que funciona correctamente sin especificar ningún parámetro en la línea de comandos.

```bash
nextflow run 3-main.nf
```

??? success "Salida del comando"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `3-main.nf` [disturbed_einstein] DSL2 - revision: ede9037d02

    executor >  local (8)
    [f0/35723c] sayHello (2)       | 3 of 3 ✔
    [40/3efd1a] convertToUpper (3) | 3 of 3 ✔
    [17/e97d32] collectGreetings   | 1 of 1 ✔
    [98/c6b57b] cowpy              | 1 of 1 ✔
    ```

Esto aún produce la misma salida que anteriormente.

La salida final de arte ASCII está en el directorio `results/3-main/`, bajo el nombre `cowpy-COLLECTED-batch-output.txt`, igual que antes.

??? abstract "Contenido del archivo"

    ```console title="results/3-main/cowpy-COLLECTED-batch-output.txt"
    _________
    / HOLà    \
    | HELLO   |
    \ BONJOUR /
    ---------
      \                                  ,+*^^*+___+++_
      \                           ,*^^^^              )
        \                       _+*                     ^**+_
        \                    +^       _ _++*+_+++_,         )
                  _+^^*+_    (     ,+*^ ^          \+_        )
                {       )  (    ,(    ,_+--+--,      ^)      ^\
                { (\@)    } f   ,(  ,+-^ __*_*_  ^^\_   ^\       )
              {:;-/    (_+*-+^^^^^+*+*<_ _++_)_    )    )      /
              ( /  (    (        ,___    ^*+_+* )   <    <      \
              U _/     )    *--<  ) ^\-----++__)   )    )       )
                (      )  _(^)^^))  )  )\^^^^^))^*+/    /       /
              (      /  (_))_^)) )  )  ))^^^^^))^^^)__/     +^^
            (     ,/    (^))^))  )  ) ))^^^^^^^))^^)       _)
              *+__+*       (_))^)  ) ) ))^^^^^^))^^^^^)____*^
              \             \_)^)_)) ))^^^^^^^^^^))^^^^)
              (_             ^\__^^^^^^^^^^^^))^^^^^^^)
                ^\___            ^\__^^^^^^))^^^^^^^^)\\
                      ^^^^^\uuu/^^\uuu/^^^^\^\^\^\^\^\^\^\
                        ___) >____) >___   ^\_\_\_\_\_\_\)
                        ^^^//\\_^^//\\_^       ^(\_\_\_\)
                          ^^^ ^^ ^^^ ^
    ```

Funcionalmente, este movimiento no ha cambiado nada, pero conceptualmente es un poco más limpio tener los valores predeterminados establecidos en el archivo de configuración.

### 1.2. Usar un archivo de configuración específico para la ejecución

??? example "Escenario"

    Quieres experimentar con diferentes configuraciones sin modificar tu archivo de configuración principal.

Puedes hacer eso creando un nuevo archivo `nextflow.config` en un subdirectorio que usarás como directorio de trabajo para tus experimentos.

#### 1.2.1. Crear el directorio de trabajo con una configuración en blanco

Comencemos creando un nuevo directorio y moviéndonos a él:

```bash
mkdir -p tux-run
cd tux-run
```

Luego, crea un archivo de configuración en blanco en ese directorio:

```bash
touch nextflow.config
```

Esto produce un archivo vacío.

#### 1.2.2. Configurar la configuración experimental

Ahora abre el nuevo archivo y agrega los parámetros que deseas personalizar:

```groovy title="tux-run/nextflow.config" linenums="1"
params {
    input = '../data/greetings.csv'
    batch = 'experiment'
    character = 'tux'
}
```

Nota que la ruta al archivo de entrada debe reflejar la estructura de directorios.

#### 1.2.3. Ejecutar el pipeline

Ahora podemos ejecutar nuestro pipeline desde dentro de nuestro nuevo directorio de trabajo.
¡Asegúrate de adaptar la ruta en consecuencia!

```bash
nextflow run ../3-main.nf
```

??? success "Salida del comando"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `../3-main.nf` [trusting_escher] DSL2 - revision: 356df0818d

    executor >  local (8)
    [59/b66913] sayHello (2)       [100%] 3 of 3 ✔
    [ad/f06364] convertToUpper (3) [100%] 3 of 3 ✔
    [10/714895] collectGreetings   [100%] 1 of 1 ✔
    [88/3ece98] cowpy              [100%] 1 of 1 ✔
    ```

Esto creará un nuevo conjunto de directorios bajo `tux-run/` incluyendo `tux-run/work/` y `tux-run/results/`.

En esta ejecución, Nextflow combina el `nextflow.config` en nuestro directorio actual con el `nextflow.config` en el directorio raíz del pipeline, y por lo tanto anula el carácter predeterminado (turkey) con el carácter tux.

El archivo de salida final debe contener el carácter tux diciendo los saludos.

??? abstract "Contenido del archivo"

    ```console title="tux-run/results/3-main/cowpy-COLLECTED-experiment-output.txt"
    _________
    / HELLO   \
    | BONJOUR |
    \ HOLà    /
    ---------
      \
        \
            .--.
          |o_o |
          |:_/ |
          //   \ \
        (|     | )
        /'\_   _/`\
        \___)=(___/

    ```

Eso es todo; ahora tienes un espacio para experimentar sin modificar tu configuración 'normal'.

!!! warning

    ¡Asegúrate de volver al directorio anterior antes de pasar a la siguiente sección!

    ```bash
    cd ..
    ```

Ahora veamos otra forma útil de establecer valores de parámetros.

### 1.3. Usar un archivo de parámetros

??? example "Escenario"

    Necesitas compartir parámetros de ejecución exactos con un colaborador, o registrarlos para una publicación.

El enfoque de subdirectorio funciona muy bien para experimentar, pero implica un poco de configuración y requiere que adaptes las rutas en consecuencia.
Hay un enfoque más simple para cuando quieres ejecutar tu pipeline con un conjunto específico de valores, o permitir que alguien más lo haga con un esfuerzo mínimo.

Nextflow nos permite especificar parámetros a través de un [archivo de parámetros](https://nextflow.io/docs/latest/config.html#parameter-file) en formato YAML o JSON, lo que hace muy conveniente gestionar y distribuir conjuntos alternativos de valores predeterminados, por ejemplo, así como valores de parámetros específicos de ejecución.

#### 1.3.1. Examinar el archivo de parámetros de ejemplo

Para demostrar esto, proporcionamos un archivo de parámetros de ejemplo en el directorio actual, llamado `test-params.yaml`:

```yaml title="test-params.yaml" linenums="1"
input: "data/greetings.csv"
batch: "yaml"
character: "stegosaurus"
```

Este archivo de parámetros contiene un par clave-valor para cada una de las entradas que queremos especificar.
Nota el uso de dos puntos (`:`) en lugar de signos de igual (`=`) si comparas la sintaxis con el archivo de configuración.
El archivo de configuración está escrito en Groovy, mientras que el archivo de parámetros está escrito en YAML.

!!! info

    También proporcionamos una versión JSON del archivo de parámetros como ejemplo, pero no vamos a ejecutarla aquí.
    Siéntete libre de probar esa por tu cuenta.

#### 1.3.2. Ejecutar el pipeline

Para ejecutar el workflow con este archivo de parámetros, simplemente agrega `-params-file <nombre_archivo>` al comando base.

```bash
nextflow run 3-main.nf -params-file test-params.yaml
```

??? success "Salida del comando"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `3-main.nf` [disturbed_sammet] DSL2 - revision: ede9037d02

    executor >  local (8)
    [2b/9a7d1e] sayHello (2)       | 3 of 3 ✔
    [5c/8f3b2a] convertToUpper (3) | 3 of 3 ✔
    [a3/29d8fb] collectGreetings   | 1 of 1 ✔
    [b7/83ef12] cowpy              | 1 of 1 ✔
    ```

El archivo de salida final debe contener el carácter stegosaurus diciendo los saludos.

??? abstract "Contenido del archivo"

    ```console title="results/3-main/cowpy-COLLECTED-yaml-output.txt"
    _________
    / HELLO   \
    | HOLà    |
    \ BONJOUR /
    ---------
    \                             .       .
    \                           / `.   .' "
      \                  .---.  <    > <    >  .---.
      \                 |    \  \ - ~ ~ - /  /    |
            _____          ..-~             ~-..-~
            |     |   \~~~\.'                    `./~~~/
          ---------   \__/                        \__/
          .'  O    \     /               /       \  "
        (_____,    `._.'               |         }  \/~~~/
          `----.          /       }     |        /    \__/
                `-.      |       /      |       /      `. ,~~|
                    ~-.__|      /_ - ~ ^|      /- _      `..-'
                        |     /        |     /     ~-.     `-. _  _  _
                        |_____|        |_____|         ~ - . _ _ _ _ _>
    ```

Usar un archivo de parámetros puede parecer excesivo cuando solo tienes unos pocos parámetros para especificar, pero algunos pipelines esperan docenas de parámetros.
En esos casos, usar un archivo de parámetros nos permitirá proporcionar valores de parámetros en tiempo de ejecución sin tener que escribir líneas de comandos masivas y sin modificar el script del workflow.

También facilita distribuir conjuntos de parámetros a colaboradores, o como información de apoyo para una publicación, por ejemplo.
Esto hace que tu trabajo sea más reproducible por otros.

### Conclusión

Sabes cómo aprovechar las opciones clave de configuración para gestionar las entradas del workflow.

### ¿Qué sigue?

Aprende cómo gestionar dónde y cómo se publican las salidas de tu workflow.

---

## 2. Gestionar salidas del workflow

??? example "Escenario"

    Tu pipeline publica salidas en un directorio codificado de forma fija, pero quieres organizar los resultados por proyecto o nombre de experimento sin editar el código del workflow cada vez.

El workflow que heredamos usa rutas para declaraciones de salida a nivel de workflow, lo cual no es terriblemente flexible e implica mucha repetición.

Veamos algunas formas comunes en las que podrías configurar esto para que sea más flexible.

### 2.1. Personalizar el nombre del directorio `outputDir`

Cada versión del workflow que hemos ejecutado hasta ahora ha publicado sus salidas en un subdirectorio diferente codificado en las definiciones de salida.

Cambiamos dónde estaba ese subdirectorio en la Parte 1 usando la bandera CLI `-output-dir`, pero eso sigue siendo solo una cadena estática.
En su lugar, configuremos esto en un archivo de configuración, donde podemos definir rutas dinámicas más complejas.
Podríamos crear un parámetro completamente nuevo para esto, pero usemos el parámetro `batch` ya que está justo ahí.

#### 2.1.1. Establecer un valor para `outputDir` en el archivo de configuración

La ruta que Nextflow usa para publicar salidas está controlada por la opción `outputDir`.
Para cambiar la ruta para todas las salidas, puedes establecer un valor para esta opción en el archivo de configuración `nextflow.config`.

Agrega el siguiente código al archivo `nextflow.config`:

=== "Después"

    ```groovy title="nextflow.config" linenums="9" hl_lines="10-13"
    /*
    * Pipeline parameters
    */
    params {
        input = 'data/greetings.csv'
        batch = 'batch'
        character = 'turkey'
    }

    /*
    * Output settings
    */
    outputDir = "results_config/${params.batch}"
    ```

=== "Antes"

    ```groovy title="nextflow.config" linenums="9"
    /*
    * Pipeline parameters
    */
    params {
        input = 'data/greetings.csv'
        batch = 'batch'
        character = 'turkey'
    }
    ```

Esto reemplazará la ruta predeterminada incorporada, `results/`, con `results_config/` más el valor del parámetro `batch` como subdirectorio.

Recuerda que también puedes establecer esta opción desde la línea de comandos usando el parámetro `-output-dir` en tu comando (`-o` para abreviar), pero entonces no podrías usar el valor del parámetro `batch`.
Usar la bandera CLI sobrescribirá `outputDir` en la configuración si está establecido.

#### 2.1.2. Eliminar la parte repetida de la ruta codificada

Todavía tenemos un subdirectorio codificado en las opciones de salida, así que eliminémoslo ahora.

Realiza los siguientes cambios de código en el archivo del workflow:

=== "Después"

    ```groovy title="3-main.nf" linenums="42" hl_lines="3 7 11 15 19"
    output {
        first_output {
            path 'intermediates'
            mode 'copy'
        }
        uppercased {
            path 'intermediates'
            mode 'copy'
        }
        collected {
            path 'intermediates'
            mode 'copy'
        }
        batch_report {
            path ''
            mode 'copy'
        }
        cowpy_art {
            path ''
            mode 'copy'
        }
    }
    ```

=== "Antes"

    ```groovy title="3-main.nf" linenums="42" hl_lines="3 7 11 15 19"
    output {
        first_output {
            path '3-main/intermediates'
            mode 'copy'
        }
        uppercased {
            path '3-main/intermediates'
            mode 'copy'
        }
        collected {
            path '3-main/intermediates'
            mode 'copy'
        }
        batch_report {
            path '3-main'
            mode 'copy'
        }
        cowpy_art {
            path '3-main'
            mode 'copy'
        }
    }
    ```

También podríamos haber simplemente agregado `${params.batch}` a cada ruta en lugar de modificar el `outputDir` predeterminado, pero esto es más conciso.

#### 2.1.3. Ejecutar el pipeline

Probemos que funciona correctamente, estableciendo el nombre del lote en `outdir` desde la línea de comandos.

```bash
nextflow run 3-main.nf --batch outdir
```

??? success "Salida del comando"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `3-main.nf` [amazing_church] DSL2 - revision: 6e18cd130e

    executor >  local (8)
    [9c/6a03ea] sayHello (2)       [100%] 3 of 3 ✔
    [11/9e58a6] convertToUpper (3) [100%] 3 of 3 ✔
    [c8/1977e5] collectGreetings   [100%] 1 of 1 ✔
    [38/f01eda] cowpy              [100%] 1 of 1 ✔
    ```

Esto aún produce la misma salida que anteriormente, excepto que esta vez encontramos nuestras salidas bajo `results_config/outdir/`.

??? abstract "Contenido del directorio"

    ```console
    results_config/outdir
    ├── cowpy-COLLECTED-outdir-output.txt
    ├── intermediates
    │   ├── Bonjour-output.txt
    │   ├── COLLECTED-outdir-output.txt
    │   ├── Hello-output.txt
    │   ├── Holà-output.txt
    │   ├── UPPER-Bonjour-output.txt
    │   ├── UPPER-Hello-output.txt
    │   └── UPPER-Holà-output.txt
    └── outdir-report.txt
    ```

Puedes combinar este enfoque con definiciones de ruta personalizadas para construir cualquier jerarquía de directorios que desees.

### 2.2. Organizar salidas por proceso

Una forma popular de organizar las salidas aún más es hacerlo por proceso, _es decir_, crear subdirectorios para cada proceso ejecutado en el pipeline.

#### 2.2.1. Reemplazar las rutas de salida por una referencia a nombres de proceso

Todo lo que necesitas hacer es hacer referencia al nombre del proceso como `<proceso>.name` en la declaración de ruta de salida.

Realiza los siguientes cambios en el archivo del workflow:

=== "Después"

    ```groovy title="3-main.nf" linenums="42" hl_lines="3 7 11 15 19"
    output {
        first_output {
            path { sayHello.name }
            mode 'copy'
        }
        uppercased {
            path { convertToUpper.name }
            mode 'copy'
        }
        collected {
            path { collectGreetings.name }
            mode 'copy'
        }
        batch_report {
            path { collectGreetings.name }
            mode 'copy'
        }
        cowpy_art {
            path { cowpy.name }
            mode 'copy'
        }
    }
    ```

=== "Antes"

    ```groovy title="3-main.nf" linenums="42" hl_lines="3 7 11 15 19"
    output {
        first_output {
            path 'intermediates'
            mode 'copy'
        }
        uppercased {
            path 'intermediates'
            mode 'copy'
        }
        collected {
            path 'intermediates'
            mode 'copy'
        }
        batch_report {
            path ''
            mode 'copy'
        }
        cowpy_art {
            path ''
            mode 'copy'
        }
    }
    ```

Esto elimina los elementos codificados restantes de la configuración de ruta de salida.

#### 2.2.2. Ejecutar el pipeline

Probemos que funciona correctamente, estableciendo el nombre del lote en `pnames` desde la línea de comandos.

```bash
nextflow run 3-main.nf --batch pnames
```

??? success "Salida del comando"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `3-main.nf` [jovial_mcclintock] DSL2 - revision: ede9037d02

    executor >  local (8)
    [4a/c2e6b8] sayHello (2)       | 3 of 3 ✔
    [6f/d4a172] convertToUpper (3) | 3 of 3 ✔
    [e8/4f19d7] collectGreetings   | 1 of 1 ✔
    [f2/a85c36] cowpy              | 1 of 1 ✔
    ```

Esto aún produce la misma salida que anteriormente, excepto que esta vez encontramos nuestras salidas bajo `results_config/pnames/`, y están agrupadas por proceso.

??? abstract "Contenido del directorio"

    ```console
    results_config/pnames/
    ├── collectGreetings
    │   ├── COLLECTED-pnames-output.txt
    │   └── pnames-report.txt
    ├── convertToUpper
    │   ├── UPPER-Bonjour-output.txt
    │   ├── UPPER-Hello-output.txt
    │   └── UPPER-Holà-output.txt
    ├── cowpy
    │   └── cowpy-COLLECTED-pnames-output.txt
    └── sayHello
        ├── Bonjour-output.txt
        ├── Hello-output.txt
        └── Holà-output.txt
    ```

!!! note

    Nota que aquí hemos borrado la distinción entre `intermediates` versus salidas finales en el nivel superior.
    Puedes mezclar y combinar estos enfoques e incluso incluir múltiples variables, por ejemplo estableciendo la ruta de la primera salida como `#!groovy "${params.batch}/intermediates/${sayHello.name}"`

### 2.3. Establecer el modo de publicación a nivel de workflow

Finalmente, en el espíritu de reducir la cantidad de código repetitivo, podemos reemplazar las declaraciones de `mode` por salida con una sola línea en la configuración.

#### 2.3.1. Agregar `workflow.output.mode` al archivo de configuración

Agrega el siguiente código al archivo `nextflow.config`:

=== "Después"

    ```groovy title="nextflow.config" linenums="2" hl_lines="5"
    /*
    * Output settings
    */
    outputDir = "results_config/${params.batch}"
    workflow.output.mode = 'copy'
    ```

=== "Antes"

    ```groovy title="nextflow.config" linenums="12"
    /*
    * Output settings
    */
    outputDir = "results_config/${params.batch}"
    ```

Al igual que la opción `outputDir`, dar a `workflow.output.mode` un valor en el archivo de configuración sería suficiente para anular lo que está establecido en el archivo del workflow, pero eliminemos el código innecesario de todos modos.

#### 2.3.2. Eliminar el modo de salida del archivo del workflow

Realiza los siguientes cambios en el archivo del workflow:

=== "Después"

    ```groovy title="3-main.nf" linenums="42"
    output {
        first_output {
            path { sayHello.name }
        }
        uppercased {
            path { convertToUpper.name }
        }
        collected {
            path { collectGreetings.name }
        }
        batch_report {
            path { collectGreetings.name }
        }
        cowpy_art {
            path { cowpy.name }
        }
    }
    ```

=== "Antes"

    ```groovy title="3-main.nf" linenums="42" hl_lines="3 7 11 15 19"
    output {
        first_output {
            path { sayHello.name }
            mode 'copy'
        }
        uppercased {
            path { convertToUpper.name }
            mode 'copy'
        }
        collected {
            path { collectGreetings.name }
            mode 'copy'
        }
        batch_report {
            path { collectGreetings.name }
            mode 'copy'
        }
        cowpy_art {
            path { cowpy.name }
            mode 'copy'
        }
    }
    ```

Eso es más conciso, ¿no es así?

#### 2.3.3. Ejecutar el pipeline

Probemos que funciona correctamente, estableciendo el nombre del lote en `outmode` desde la línea de comandos.

```bash
nextflow run 3-main.nf --batch outmode
```

??? success "Salida del comando"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `3-main.nf` [rowdy_sagan] DSL2 - revision: ede9037d02

    executor >  local (8)
    [5b/d91e3c] sayHello (2)       | 3 of 3 ✔
    [8a/f6c241] convertToUpper (3) | 3 of 3 ✔
    [89/cd3a48] collectGreetings   | 1 of 1 ✔
    [9e/71fb52] cowpy              | 1 of 1 ✔
    ```

Esto aún produce la misma salida que anteriormente, excepto que esta vez encontramos nuestras salidas bajo `results_config/outmode/`.
Todas siguen siendo copias adecuadas, no enlaces simbólicos.

??? abstract "Contenido del directorio"

    ```console
    results_config/outmode/
    ├── collectGreetings
    │   ├── COLLECTED-outmode-output.txt
    │   └── outmode-report.txt
    ├── convertToUpper
    │   ├── UPPER-Bonjour-output.txt
    │   ├── UPPER-Hello-output.txt
    │   └── UPPER-Holà-output.txt
    ├── cowpy
    │   └── cowpy-COLLECTED-outmode-output.txt
    └── sayHello
        ├── Bonjour-output.txt
        ├── Hello-output.txt
        └── Holà-output.txt
    ```

La razón principal por la que aún podrías querer usar la forma por salida de establecer el modo es si quieres mezclar y combinar dentro del mismo workflow, _es decir_, tener algunas salidas copiadas y algunas como enlaces simbólicos.

Hay muchas otras opciones que puedes personalizar de esta manera, pero con suerte esto te da una idea del rango de opciones y cómo utilizarlas efectivamente para adaptarse a tus preferencias.

### Conclusión

Sabes cómo controlar el nombre y la estructura de los directorios donde se publican tus salidas, así como el modo de publicación de salida del workflow.

### ¿Qué sigue?

Aprende cómo adaptar la configuración de tu workflow a tu entorno de cómputo, comenzando con la tecnología de empaquetado de software.

---

## 3. Seleccionar una tecnología de empaquetado de software

Hasta ahora hemos estado viendo elementos de configuración que controlan cómo entran las entradas y dónde salen las salidas. Ahora es momento de enfocarnos más específicamente en adaptar la configuración de tu workflow a tu entorno de cómputo.

El primer paso en ese camino es especificar de dónde vendrán los paquetes de software que se ejecutarán en cada paso.
¿Ya están instalados en el entorno de cómputo local?
¿Necesitamos recuperar imágenes y ejecutarlas a través de un sistema de contenedores?
¿O necesitamos recuperar paquetes de Conda y construir un entorno Conda local?

En la primera parte de este curso de capacitación (Partes 1-4) simplemente usamos software instalado localmente en nuestro workflow.
Luego en la Parte 5, introdujimos contenedores Docker y el archivo `nextflow.config`, que usamos para habilitar el uso de contenedores Docker.

Ahora veamos cómo podemos configurar una opción alternativa de empaquetado de software a través del archivo `nextflow.config`.

### 3.1. Deshabilitar Docker y habilitar Conda en el archivo de configuración

??? example "Escenario"

    Estás moviendo tu pipeline a un clúster HPC donde Docker no está permitido por razones de seguridad.
    El clúster soporta Singularity y Conda, así que necesitas cambiar tu configuración en consecuencia.

Como se señaló anteriormente, Nextflow soporta múltiples tecnologías de contenedores incluyendo Singularity (que se usa más ampliamente en HPC), así como gestores de paquetes de software como Conda.

Podemos cambiar nuestro archivo de configuración para usar Conda en lugar de Docker.
Para hacerlo, cambiemos el valor de `docker.enabled` a `false`, y agreguemos una directiva habilitando el uso de Conda:

=== "Después"

    ```groovy title="nextflow.config" linenums="1" hl_lines="1-2"
    docker.enabled = false
    conda.enabled = true
    ```

=== "Antes"

    ```groovy title="nextflow.config" linenums="1" hl_lines="1"
    docker.enabled = true
    ```

Esto permitirá a Nextflow crear y utilizar entornos Conda para procesos que tienen paquetes Conda especificados.
¡Lo que significa que ahora necesitamos agregar uno de esos a nuestro proceso `cowpy`!

### 3.2. Especificar un paquete Conda en la definición del proceso

Ya hemos recuperado el URI para un paquete Conda que contiene la herramienta `cowpy`: `conda-forge::cowpy==1.1.5`

Ahora agregamos el URI a la definición del proceso `cowpy` usando la directiva `conda`:

=== "Después"

    ```groovy title="modules/cowpy.nf" linenums="4" hl_lines="4"
    process cowpy {

        container 'community.wave.seqera.io/library/cowpy:1.1.5--3db457ae1977a273'
        conda 'conda-forge::cowpy==1.1.5'

        input:
    ```

=== "Antes"

    ```groovy title="modules/cowpy.nf" linenums="4"
    process cowpy {

        container 'community.wave.seqera.io/library/cowpy:1.1.5--3db457ae1977a273'

        input:
    ```

Para ser claros, no estamos _reemplazando_ la directiva `docker`, estamos _agregando_ una opción alternativa.

!!! tip

    Hay algunas formas diferentes de obtener el URI para un paquete conda dado.
    Recomendamos usar la consulta de búsqueda de [Seqera Containers](https://seqera.io/containers/), que te dará un URI que puedes copiar y pegar, incluso si no planeas crear un contenedor a partir de él.

### 3.3. Ejecutar el workflow para verificar que puede usar Conda

Probémoslo.

```bash
nextflow run 3-main.nf --batch conda
```

??? success "Salida del comando"

    ```console title="Output"
    N E X T F L O W   ~  version 25.10.2

    Launching `3-main.nf` [trusting_lovelace] DSL2 - revision: 028a841db1

    executor >  local (8)
    [ee/4ca1f2] sayHello (3)       | 3 of 3 ✔
    [20/2596a7] convertToUpper (1) | 3 of 3 ✔
    [b3/e15de5] collectGreetings   | 1 of 1 ✔
    [c5/af5f88] cowpy              | 1 of 1 ✔
    ```

Esto debería funcionar sin problemas y producir las mismas salidas que anteriormente bajo `results_config/conda`.

Detrás de escena, Nextflow ha recuperado los paquetes Conda y creado el entorno, lo que normalmente requiere un poco de trabajo; ¡así que es agradable que no tengamos que hacer nada de eso nosotros mismos!

!!! info

    Esto se ejecuta rápidamente porque el paquete `cowpy` es bastante pequeño, pero si estás trabajando con paquetes grandes, puede tomar un poco más de tiempo de lo habitual la primera vez, y podrías ver la salida de la consola quedarse 'atascada' durante un minuto más o menos antes de completarse.
    Esto es normal y se debe al trabajo extra que Nextflow hace la primera vez que usas un nuevo paquete.

Desde nuestro punto de vista, parece que funciona exactamente igual que ejecutar con Docker, aunque en el backend la mecánica es un poco diferente.

Esto significa que estamos listos para ejecutar con entornos Conda si es necesario.

??? info "Mezclar y combinar Docker y Conda"

    Dado que estas directivas se asignan por proceso, es posible 'mezclar y combinar', _es decir_, configurar algunos de los procesos en tu workflow para ejecutarse con Docker y otros con Conda, por ejemplo, si la infraestructura de cómputo que estás usando soporta ambos.
    En ese caso, habilitarías tanto Docker como Conda en tu archivo de configuración.
    Si ambos están disponibles para un proceso dado, Nextflow priorizará los contenedores.

    Y como se señaló anteriormente, Nextflow soporta múltiples otras tecnologías de empaquetado de software y contenedores, por lo que no estás limitado a solo esas dos.

### Conclusión

Sabes cómo configurar qué paquete de software debe usar cada proceso, y cómo cambiar entre tecnologías.

### ¿Qué sigue?

Aprende cómo cambiar la plataforma de ejecución utilizada por Nextflow para hacer realmente el trabajo.

---

## 4. Seleccionar una plataforma de ejecución

??? example "Escenario"

    Has estado desarrollando y probando tu pipeline en tu laptop, pero ahora necesitas ejecutarlo en miles de muestras.
    Tu institución tiene un clúster HPC con un planificador Slurm que te gustaría usar en su lugar.

Hasta ahora, hemos estado ejecutando nuestro pipeline con el ejecutor local.
Esto ejecuta cada tarea en la máquina en la que se está ejecutando Nextflow.
Cuando Nextflow comienza, examina los CPUs y la memoria disponibles.
Si los recursos de las tareas listas para ejecutarse exceden los recursos disponibles, Nextflow retendrá las últimas tareas de la ejecución hasta que una o más de las tareas anteriores hayan terminado, liberando los recursos necesarios.

El ejecutor local es conveniente y eficiente, pero está limitado a esa única máquina. Para cargas de trabajo muy grandes, puedes descubrir que tu máquina local es un cuello de botella, ya sea porque tienes una sola tarea que requiere más recursos de los que tienes disponibles, o porque tienes tantas tareas que esperar a que una sola máquina las ejecute tomaría demasiado tiempo.

Nextflow soporta [muchos backends de ejecución diferentes](https://nextflow.io/docs/latest/executor.html), incluyendo planificadores HPC (Slurm, LSF, SGE, PBS, Moab, OAR, Bridge, HTCondor y otros) así como backends de ejecución en la nube (AWS Batch, Google Cloud Batch, Azure Batch, Kubernetes y más).

### 4.1. Apuntar a un backend diferente

La elección del ejecutor se establece mediante una directiva de proceso llamada `executor`.
Por defecto está establecido en `local`, por lo que la siguiente configuración está implícita:

```groovy title="Built-in configuration"
process {
    executor = 'local'
}
```

Para establecer el ejecutor para apuntar a un backend diferente, simplemente especificarías el ejecutor que deseas usando una sintaxis similar a la descrita anteriormente para asignaciones de recursos (consulta [Executors](https://nextflow.io/docs/latest/executor.html) para todas las opciones).

```groovy title="nextflow.config"
process {
    executor = 'slurm'
}
```

!!! warning

    En realidad no podemos probar esto en el entorno de capacitación porque no está configurado para conectarse a un HPC.

### 4.2. Lidiar con sintaxis específica del backend para parámetros de ejecución

La mayoría de las plataformas de computación de alto rendimiento permiten (y a veces requieren) que especifiques ciertos parámetros como solicitudes y limitaciones de asignación de recursos (por ejemplo, número de CPUs y memoria) y nombre de la cola de trabajos a usar.

Desafortunadamente, cada uno de estos sistemas usa diferentes tecnologías, sintaxis y configuraciones para definir cómo debe definirse y enviarse un trabajo al planificador relevante.

??? abstract "Ejemplos"

    Por ejemplo, el mismo trabajo que requiere 8 CPUs y 4GB de RAM para ejecutarse en la cola "my-science-work" necesita expresarse de las siguientes formas diferentes dependiendo del backend.

    ```bash title="Config for SLURM / submit using sbatch"
    #SBATCH -o /path/to/my/task/directory/my-task-1.log
    #SBATCH --no-requeue
    #SBATCH -c 8
    #SBATCH --mem 4096M
    #SBATCH -p my-science-work
    ```

    ```bash title="Config for PBS / submit using qsub"
    #PBS -o /path/to/my/task/directory/my-task-1.log
    #PBS -j oe
    #PBS -q my-science-work
    #PBS -l nodes=1:ppn=5
    #PBS -l mem=4gb
    ```

    ```bash title="Config for SGE / submit using qsub"
    #$ -o /path/to/my/task/directory/my-task-1.log
    #$ -j y
    #$ -terse
    #$ -notify
    #$ -q my-science-work
    #$ -l slots=5
    #$ -l h_rss=4096M,mem_free=4096M
    ```

Afortunadamente, Nextflow simplifica todo esto.
Proporciona una sintaxis estandarizada para que puedas especificar las propiedades relevantes como `cpus`, `memory` y `queue` solo una vez (consulta [Process directives](https://nextflow.io/docs/latest/reference/process.html#process-directives) para todas las opciones disponibles).
Luego, en tiempo de ejecución, Nextflow usará esas configuraciones para generar los scripts específicos del backend apropiados basados en la configuración del ejecutor.

Cubriremos esa sintaxis estandarizada en la siguiente sección.

### Conclusión

Ahora sabes cómo cambiar el ejecutor para usar diferentes tipos de infraestructura de cómputo.

### ¿Qué sigue?

Aprende cómo evaluar y expresar asignaciones y limitaciones de recursos en Nextflow.

---

## 5. Controlar asignaciones de recursos de cómputo

??? example "Escenario"

    Tu pipeline sigue fallando en el clúster porque las tareas están siendo eliminadas por exceder los límites de memoria.
    O tal vez te están cobrando por recursos que no estás usando y quieres optimizar costos.

La mayoría de las plataformas de computación de alto rendimiento permiten (y a veces requieren) que especifiques ciertos parámetros de asignación de recursos como número de CPUs y memoria.

Por defecto, Nextflow usará un solo CPU y 2GB de memoria para cada proceso.
Las directivas de proceso correspondientes se llaman `cpus` y `memory`, por lo que la siguiente configuración está implícita:

```groovy title="Built-in configuration" linenums="1"
process {
    cpus = 1
    memory = 2.GB
}
```

Puedes modificar estos valores, ya sea para todos los procesos o para procesos específicos nombrados, usando directivas de proceso adicionales en tu archivo de configuración.
Nextflow las traducirá en las instrucciones apropiadas para el ejecutor elegido.

Pero, ¿cómo sabes qué valores usar?

### 5.1. Ejecutar el workflow para generar un informe de utilización de recursos

??? example "Escenario"

    No sabes cuánta memoria o CPU necesitan tus procesos y quieres evitar desperdiciar recursos o que los trabajos sean eliminados.

Si no sabes de antemano cuánto CPU y memoria es probable que necesiten tus procesos, puedes hacer un perfilado de recursos, lo que significa que ejecutas el workflow con algunas asignaciones predeterminadas, registras cuánto usó cada proceso y, a partir de ahí, estimas cómo ajustar las asignaciones base.

Convenientemente, Nextflow incluye herramientas incorporadas para hacer esto, y generará felizmente un informe para ti bajo solicitud.

Para hacerlo, agrega `-with-report <nombre_archivo>.html` a tu línea de comandos.

```bash
nextflow run 3-main.nf -with-report report-config-1.html
```

El informe es un archivo html, que puedes descargar y abrir en tu navegador. También puedes hacer clic derecho en él en el explorador de archivos a la izquierda y hacer clic en `Show preview` para verlo en el entorno de capacitación.

Tómate unos minutos para revisar el informe y ver si puedes identificar algunas oportunidades para ajustar recursos.
Asegúrate de hacer clic en las pestañas que muestran los resultados de utilización como un porcentaje de lo que fue asignado.

Consulta [Reports](https://nextflow.io/docs/latest/reports.html) para documentación sobre todas las características disponibles.

### 5.2. Establecer asignaciones de recursos para todos los procesos

El perfilado muestra que los procesos en nuestro workflow de capacitación son muy ligeros, así que reduzcamos la asignación de memoria predeterminada a 1GB por proceso.

Agrega lo siguiente a tu archivo `nextflow.config`, antes de la sección de parámetros del pipeline:

```groovy title="nextflow.config" linenums="4"
/*
* Process settings
*/
process {
    memory = 1.GB
}
```

Eso ayudará a reducir la cantidad de cómputo que consumimos.

### 5.3. Establecer asignaciones de recursos para un proceso específico

Al mismo tiempo, vamos a pretender que el proceso `cowpy` requiere más recursos que los demás, solo para que podamos demostrar cómo ajustar asignaciones para un proceso individual.

=== "Después"

    ```groovy title="nextflow.config" linenums="4" hl_lines="6-9"
    /*
    * Process settings
    */
    process {
        memory = 1.GB
        withName: 'cowpy' {
            memory = 2.GB
            cpus = 2
        }
    }
    ```

=== "Antes"

    ```groovy title="nextflow.config" linenums="4"
    /*
    * Process settings
    */
    process {
        memory = 1.GB
    }
    ```

Con esta configuración, todos los procesos solicitarán 1GB de memoria y un solo CPU (el predeterminado implícito), excepto el proceso `cowpy`, que solicitará 2GB y 2 CPUs.

!!! info

    Si tienes una máquina con pocos CPUs y asignas un número alto por proceso, podrías ver llamadas de proceso poniéndose en cola una detrás de otra.
    Esto es porque Nextflow asegura que no solicitemos más CPUs de los que están disponibles.

### 5.4. Ejecutar el workflow con la configuración actualizada

Probemos eso, proporcionando un nombre de archivo diferente para el informe de perfilado para que podamos comparar el rendimiento antes y después de los cambios de configuración.

```bash
nextflow run 3-main.nf -with-report report-config-2.html
```

Probablemente no notarás ninguna diferencia real ya que esta es una carga de trabajo tan pequeña, pero este es el enfoque que usarías para analizar el rendimiento y los requisitos de recursos de un workflow del mundo real.

Es muy útil cuando tus procesos tienen diferentes requisitos de recursos. Te permite dimensionar correctamente las asignaciones de recursos que configuras para cada proceso basándote en datos reales, no en conjeturas.

!!! tip

    Esto es solo una pequeña muestra de lo que puedes hacer para optimizar tu uso de recursos.
    Nextflow en sí tiene una [lógica de reintento dinámico](https://nextflow.io/docs/latest/process.html#dynamic-task-resources) realmente ingeniosa incorporada para reintentar trabajos que fallan debido a limitaciones de recursos.
    Además, la Seqera Platform ofrece herramientas impulsadas por IA para optimizar tus asignaciones de recursos automáticamente también.

### 5.5. Agregar límites de recursos

Dependiendo de qué ejecutor de cómputo e infraestructura de cómputo estés usando, puede haber algunas restricciones sobre lo que puedes (o debes) asignar.
Por ejemplo, tu clúster puede requerir que te mantengas dentro de ciertos límites.

Puedes usar la directiva `resourceLimits` para establecer las limitaciones relevantes. La sintaxis se ve así cuando está sola en un bloque de proceso:

```groovy title="Syntax example"
process {
    resourceLimits = [
        memory: 750.GB,
        cpus: 200,
        time: 30.d
    ]
}
```

Nextflow traducirá estos valores en las instrucciones apropiadas dependiendo del ejecutor que especificaste.

No vamos a ejecutar esto, ya que no tenemos acceso a infraestructura relevante en el entorno de capacitación.
Sin embargo, si intentaras ejecutar el workflow con asignaciones de recursos que excedan estos límites, luego buscar el comando `sbatch` en el archivo de script `.command.run`, verías que las solicitudes que realmente se envían al ejecutor están limitadas a los valores especificados por `resourceLimits`.

??? info "Configuraciones de referencia institucionales"

    El proyecto nf-core ha compilado una [colección de archivos de configuración](https://nf-co.re/configs/) compartidos por varias instituciones alrededor del mundo, cubriendo una amplia gama de ejecutores HPC y en la nube.

    Esas configuraciones compartidas son valiosas tanto para personas que trabajan allí y por lo tanto pueden simplemente utilizar la configuración de su institución de forma inmediata, como modelo para personas que buscan desarrollar una configuración para su propia infraestructura.

### Conclusión

Sabes cómo generar un informe de perfilado para evaluar la utilización de recursos y cómo modificar asignaciones de recursos para todos los procesos y/o para procesos individuales, así como establecer limitaciones de recursos para ejecutar en HPC.

### ¿Qué sigue?

Aprende cómo configurar perfiles de configuración preestablecidos y cambiar entre ellos en tiempo de ejecución.

---

## 6. Usar perfiles para cambiar entre configuraciones preestablecidas

??? example "Escenario"

    Cambias regularmente entre ejecutar pipelines en tu laptop para desarrollo y en el HPC de tu institución para ejecuciones de producción.
    Estás cansado de cambiar manualmente las configuraciones cada vez que cambias de entorno.

Te hemos mostrado varias formas en las que puedes personalizar la configuración de tu pipeline dependiendo del proyecto en el que estés trabajando o del entorno de cómputo que estés usando.

Es posible que desees cambiar entre configuraciones alternativas dependiendo de qué infraestructura de cómputo estés usando. Por ejemplo, podrías querer desarrollar y ejecutar pruebas a pequeña escala localmente en tu laptop, luego ejecutar cargas de trabajo a escala completa en HPC o en la nube.

Nextflow te permite configurar cualquier número de [**perfiles**](https://nextflow.io/docs/latest/config.html#profiles) que describen diferentes configuraciones, que luego puedes seleccionar en tiempo de ejecución usando un argumento de línea de comandos, en lugar de tener que modificar el archivo de configuración en sí.

### 6.1. Crear perfiles para cambiar entre desarrollo local y ejecución en HPC

Configuremos dos perfiles alternativos; uno para ejecutar cargas a pequeña escala en una computadora regular, donde usaremos contenedores Docker, y uno para ejecutar en un HPC universitario con un planificador Slurm, donde usaremos paquetes Conda.

#### 6.1.1. Configurar los perfiles

Agrega lo siguiente a tu archivo `nextflow.config`, después de la sección de parámetros del pipeline pero antes de la configuración de salida:

```groovy title="nextflow.config" linenums="24"
/*
* Profiles
*/
profiles {
    my_laptop {
        process.executor = 'local'
        docker.enabled = true
    }
    univ_hpc {
        process.executor = 'slurm'
        conda.enabled = true
        process.resourceLimits = [
            memory: 750.GB,
            cpus: 200,
            time: 30.d
        ]
    }
}
```

Ves que para el HPC universitario, también estamos especificando limitaciones de recursos.

#### 6.1.2. Ejecutar el workflow con un perfil

Para especificar un perfil en nuestra línea de comandos de Nextflow, usamos el argumento `-profile`.

Intentemos ejecutar el workflow con la configuración `my_laptop`.

```bash
nextflow run 3-main.nf -profile my_laptop
```

??? success "Salida del comando"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `3-main.nf` [gigantic_brazil] DSL2 - revision: ede9037d02

    executor >  local (8)
    [58/da9437] sayHello (3)       | 3 of 3 ✔
    [35/9cbe77] convertToUpper (2) | 3 of 3 ✔
    [67/857d05] collectGreetings   | 1 of 1 ✔
    [37/7b51b5] cowpy              | 1 of 1 ✔
    ```

Como puedes ver, esto nos permite alternar entre configuraciones muy convenientemente en tiempo de ejecución.

!!! warning

    El perfil `univ_hpc` no se ejecutará correctamente en el entorno de capacitación ya que no tenemos acceso a un planificador Slurm.

Si en el futuro encontramos otros elementos de configuración que siempre co-ocurren con estos, podemos simplemente agregarlos al(los) perfil(es) correspondiente(s).
También podemos crear perfiles adicionales si hay otros elementos de configuración que queremos agrupar.

### 6.2. Crear un perfil de parámetros de prueba

??? example "Escenario"

    Quieres que otros puedan probar tu pipeline rápidamente sin reunir sus propios datos de entrada.

Los perfiles no son solo para configuración de infraestructura.
También podemos usarlos para establecer valores predeterminados para parámetros del workflow, para facilitar que otros prueben el workflow sin tener que reunir valores de entrada apropiados ellos mismos.
Puedes considerar esto una alternativa a usar un archivo de parámetros.

#### 6.2.1. Configurar el perfil

La sintaxis para expresar valores predeterminados en este contexto se ve así, para un perfil que nombramos `test`:

```groovy title="Syntax example"
    test {
        params.<parameter1>
        params.<parameter2>
        ...
    }
```

Si agregamos un perfil de prueba para nuestro workflow, el bloque `profiles` se convierte en:

```groovy title="nextflow.config" linenums="24"
/*
* Profiles
*/
profiles {
    my_laptop {
        process.executor = 'local'
        docker.enabled = true
    }
    univ_hpc {
        process.executor = 'slurm'
        conda.enabled = true
        process.resourceLimits = [
            memory: 750.GB,
            cpus: 200,
            time: 30.d
        ]
    }
    test {
        params.input = 'data/greetings.csv'
        params.batch = 'test'
        params.character = 'dragonandcow'
    }
}
```

Al igual que para los perfiles de configuración técnica, puedes configurar múltiples perfiles diferentes especificando parámetros bajo cualquier nombre arbitrario que desees.

#### 6.2.2. Ejecutar el workflow localmente con el perfil de prueba

Convenientemente, los perfiles no son mutuamente excluyentes, por lo que podemos especificar múltiples perfiles en nuestra línea de comandos usando la siguiente sintaxis `-profile <perfil1>,<perfil2>` (para cualquier número de perfiles).

Si combinas perfiles que establecen valores para los mismos elementos de configuración y están descritos en el mismo archivo de configuración, Nextflow resolverá el conflicto usando el valor que leyó al final (_es decir_, lo que viene después en el archivo).
Si las configuraciones en conflicto se establecen en diferentes fuentes de configuración, se aplica el [orden de precedencia](https://www.nextflow.io/docs/latest/config.html) predeterminado.

Intentemos agregar el perfil de prueba a nuestro comando anterior:

```bash
nextflow run 3-main.nf -profile my_laptop,test
```

??? success "Salida del comando"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `3-main.nf` [jovial_coulomb] DSL2 - revision: 46a6763141

    executor >  local (8)
    [9b/687cdc] sayHello (2)       | 3 of 3 ✔
    [ca/552187] convertToUpper (3) | 3 of 3 ✔
    [e8/83e306] collectGreetings   | 1 of 1 ✔
    [fd/e84fa9] cowpy              | 1 of 1 ✔
    ```

Esto usará Docker donde sea posible y producirá salidas bajo `results_config/test`, y esta vez el carácter es el dúo cómico `dragonandcow`.

??? abstract "Contenido del archivo"

    ```console title="results_config/test/"
     _________
    / HOLà    \
    | HELLO   |
    \ BONJOUR /
    ---------
                \                    ^    /^
                  \                  / \  // \
                  \   |\___/|      /   \//  .\
                    \  /O  O  \__  /    //  | \ \           *----*
                      /     /  \/_/    //   |  \  \          \   |
                      \@___\@`    \/_   //    |   \   \         \/\ \
                    0/0/|       \/_ //     |    \    \         \ \
                0/0/0/0/|        \///      |     \     \       | |
              0/0/0/0/0/_|_ /   (  //       |      \     _\     |  /
          0/0/0/0/0/0/`/,_ _ _/  ) ; -.    |    _ _\.-~       /   /
                      ,-}        _      *-.|.-~-.           .~    ~
      \     \__/        `/\      /                 ~-. _ .-~      /
      \____(oo)           *.   }            {                   /
      (    (--)          .----~-.\        \-`                 .~
      //__\\  \__ Ack!   ///.----..<        \             _ -~
      //    \\               ///-._ _ _ _ _ _ _{^ - - - - ~
    ```

Esto significa que mientras distribuyamos cualquier archivo de datos de prueba con el código del workflow, cualquiera puede probar rápidamente el workflow sin tener que proporcionar sus propias entradas a través de la línea de comandos o un archivo de parámetros.

!!! tip

    Podemos apuntar a URLs para archivos más grandes que están almacenados externamente.
    Nextflow los descargará automáticamente siempre que haya una conexión abierta.

    Para más detalles, consulta la Side Quest [Working with Files](../side_quests/working_with_files.md)

### 6.3. Usar `nextflow config` para ver la configuración resuelta

Como se señaló anteriormente, a veces el mismo parámetro puede establecerse en diferentes valores en perfiles que deseas combinar.
Y más generalmente, hay numerosos lugares donde pueden almacenarse elementos de configuración, y a veces las mismas propiedades pueden establecerse en diferentes valores en diferentes lugares.

Nextflow aplica un [orden de precedencia](https://nextflow.io/docs/latest/config.html#configuration-file) establecido para resolver cualquier conflicto, pero eso puede ser difícil de determinar por ti mismo.
E incluso si nada está en conflicto, puede ser tedioso buscar todos los lugares posibles donde las cosas podrían estar configuradas.

Afortunadamente, Nextflow incluye una herramienta de utilidad conveniente llamada `config` que puede automatizar todo ese proceso para ti.

La herramienta `config` explorará todos los contenidos en tu directorio de trabajo actual, recogerá cualquier archivo de configuración y producirá la configuración completamente resuelta que Nextflow usaría para ejecutar el workflow.
Esto te permite averiguar qué configuraciones se usarán sin tener que lanzar nada.

#### 6.3.1. Resolver la configuración predeterminada

Ejecuta este comando para resolver la configuración que se aplicaría por defecto.

```bash
nextflow config
```

??? success "Salida del comando"

    ```groovy
    params {
      input = 'data/greetings.csv'
      batch = 'batch'
      character = 'turkey'
    }

    docker {
      enabled = false
    }

    conda {
      enabled = true
    }

    process {
      memory = '1 GB'
      withName:cowpy {
          memory = '2 GB'
          cpus = 2
      }
    }

    outputDir = 'results_config/batch'

    workflow {
      output {
          mode = 'copy'
      }
    }
    ```

Esto te muestra la configuración base que obtienes si no especificas nada extra en la línea de comandos.

#### 6.3.2. Resolver la configuración con configuraciones específicas activadas

Si proporcionas parámetros de línea de comandos, por ejemplo, habilitando uno o más perfiles o cargando un archivo de parámetros, el comando adicionalmente los tomará en cuenta.

```bash
nextflow config -profile my_laptop,test
```

??? success "Salida del comando"

    ```groovy
    params {
      input = 'data/greetings.csv'
      batch = 'test'
      character = 'dragonandcow'
    }

    docker {
      enabled = true
    }

    conda {
      enabled = true
    }

    process {
      memory = '1 GB'
      withName:cowpy {
          memory = '2 GB'
          cpus = 2
      }
      executor = 'local'
    }

    outputDir = 'results_config/test'

    workflow {
      output {
          mode = 'copy'
      }
    }
    ```

Esto se vuelve especialmente útil para proyectos complejos que involucran múltiples capas de configuración.

### Conclusión

Sabes cómo usar perfiles para seleccionar una configuración preestablecida en tiempo de ejecución con un mínimo de complicaciones.
Más generalmente, sabes cómo configurar tus ejecuciones de workflow para adaptarse a diferentes plataformas de cómputo y mejorar la reproducibilidad de tus análisis.

### ¿Qué sigue?

Aprende cómo ejecutar pipelines directamente desde repositorios remotos como GitHub.

---

## 7. Ejecutar pipelines desde repositorios remotos

??? example "Escenario"

    Quieres ejecutar un pipeline bien establecido como los de nf-core sin tener que descargar y gestionar el código tú mismo.

Hasta ahora hemos estado ejecutando scripts de workflow ubicados en el directorio actual.
En la práctica, a menudo querrás ejecutar pipelines almacenados en repositorios remotos, como GitHub.

Nextflow hace esto sencillo: puedes ejecutar cualquier pipeline directamente desde una URL de repositorio Git sin descargarlo manualmente primero.

### 7.1. Ejecutar un pipeline desde GitHub

La sintaxis básica para ejecutar un pipeline remoto es `nextflow run <repositorio>`, donde `<repositorio>` puede ser una ruta de repositorio de GitHub como `nextflow-io/hello`, una URL completa, o una ruta a GitLab, Bitbucket u otros servicios de alojamiento Git.

Intenta ejecutar el pipeline de demostración oficial "hello" de Nextflow:

```bash
nextflow run nextflow-io/hello
```

??? success "Salida del comando"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Pulling nextflow-io/hello ...
     downloaded from https://github.com/nextflow-io/hello.git
    Launching `https://github.com/nextflow-io/hello` [sleepy_swanson] DSL2 - revision: 2ce0b0e294 [master]

    executor >  local (4)
    [ba/08236d] sayHello (4) [100%] 4 of 4 ✔
    Ciao world!

    Hello world!

    Bonjour world!

    Hola world!
    ```

La primera vez que ejecutas un pipeline remoto, Nextflow lo descarga y lo almacena en caché localmente.
Las ejecuciones subsiguientes usan la versión en caché a menos que solicites explícitamente una actualización.

### 7.2. Especificar una versión para reproducibilidad

Por defecto, Nextflow ejecuta la última versión de la rama predeterminada.
Puedes especificar una versión particular (etiqueta), rama o commit usando la bandera `-r`:

```bash
nextflow run nextflow-io/hello -r v1.3
```

??? success "Salida del comando"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `https://github.com/nextflow-io/hello` [sick_carson] DSL2 - revision: 2ce0b0e294 [v1.3]

    executor >  local (4)
    [61/e11f77] sayHello (4) [100%] 4 of 4 ✔
    Ciao world!

    Bonjour world!

    Hello world!

    Hola world!
    ```

Especificar versiones exactas es esencial para la reproducibilidad.

### Conclusión

Sabes cómo ejecutar pipelines directamente desde GitHub y otros repositorios remotos, y cómo especificar versiones para reproducibilidad.

### ¿Qué sigue?

¡Date una gran palmadita en la espalda!
Sabes todo lo que necesitas saber para comenzar a ejecutar y gestionar pipelines de Nextflow.

Eso concluye este curso, pero si estás ansioso por seguir aprendiendo, tenemos dos recomendaciones principales:

- Si quieres profundizar en el desarrollo de tus propios pipelines, echa un vistazo a [Hello Nextflow](../hello_nextflow/index.md), un curso para principiantes que cubre la misma progresión general que este pero entra en mucho más detalle sobre canales y operadores.
- Si te gustaría continuar aprendiendo cómo ejecutar pipelines de Nextflow sin profundizar más en el código, echa un vistazo a la primera parte de [Hello nf-core](../hello_nf-core/index.md), que introduce las herramientas para encontrar y ejecutar pipelines del enormemente popular proyecto [nf-core](https://nf-co.re/).

¡Diviértete!

---

## Quiz

<quiz>
Cuando los valores de parámetros se establecen tanto en el archivo del workflow como en `nextflow.config`, ¿cuál tiene precedencia?
- [ ] El valor del archivo del workflow
- [x] El valor del archivo de configuración
- [ ] El primer valor encontrado
- [ ] Causa un error

Aprende más: [1.1. Configurar valores en `nextflow.config`](#11-set-up-values-in-nextflowconfig)
</quiz>

<quiz>
¿Cuál es la diferencia de sintaxis entre establecer un parámetro predeterminado en un archivo de workflow vs. un archivo de configuración?
- [ ] Usan la misma sintaxis
- [x] El workflow usa declaración tipada (`#!groovy param: Type = value`), la configuración usa asignación (`#!groovy param = value`)
- [ ] La configuración usa declaración tipada, el workflow usa asignación
- [ ] Solo los archivos de configuración pueden establecer valores predeterminados

Aprende más: [1.1. Configurar valores en `nextflow.config`](#11-set-up-values-in-nextflowconfig)
</quiz>

<quiz>
¿Cómo especificas un archivo de parámetros al ejecutar un workflow?
- [ ] `--params params.yaml`
- [ ] `-config params.yaml`
- [x] `-params-file params.yaml`
- [ ] `--input-params params.yaml`

Aprende más: [1.3. Usar un archivo de parámetros](#13-use-a-parameter-file)
</quiz>

<quiz>
¿Qué controla la opción de configuración `outputDir`?
- [ ] La ubicación del directorio de trabajo
- [x] La ruta base donde se publican las salidas del workflow
- [ ] El directorio para archivos de registro
- [ ] La ubicación de archivos de módulos

Aprende más: [2.1. Personalizar el nombre del directorio outputDir](#21-customize-the-outputdir-directory-name)
</quiz>

<quiz>
¿Cómo haces referencia a un nombre de proceso dinámicamente en la configuración de ruta de salida?
- [ ] `#!groovy ${processName}`
- [ ] `#!groovy path "<proceso>.name"`
- [x] `#!groovy path { <proceso>.name }`
- [ ] `@processName`

Aprende más: [2.2. Organizar salidas por proceso](#22-organize-outputs-by-process)
</quiz>

<quiz>
Si tanto Docker como Conda están habilitados y un proceso tiene ambas directivas, ¿cuál se prioriza?
- [x] Docker (contenedores)
- [ ] Conda
- [ ] La primera definida en el proceso
- [ ] Causa un error

Aprende más: [3. Seleccionar una tecnología de empaquetado de software](#3-select-a-software-packaging-technology)
</quiz>

<quiz>
¿Cuál es el ejecutor predeterminado en Nextflow?
- [x] `local`
- [ ] `slurm`
- [ ] `kubernetes`
- [ ] `aws`

Aprende más: [4. Seleccionar una plataforma de ejecución](#4-select-an-execution-platform)
</quiz>

<quiz>
¿Qué comando genera un informe de utilización de recursos?
- [ ] `nextflow run workflow.nf -with-metrics`
- [ ] `nextflow run workflow.nf -with-stats`
- [x] `nextflow run workflow.nf -with-report report.html`
- [ ] `nextflow run workflow.nf -profile report`

Aprende más: [5.1. Ejecutar el workflow para generar un informe de utilización de recursos](#51-run-the-workflow-to-generate-a-resource-utilization-report)
</quiz>

<quiz>
¿Cómo estableces requisitos de recursos para un proceso específico llamado `cowpy` en el archivo de configuración?
- [ ] `#!groovy cowpy.memory = '2.GB'`
- [ ] `#!groovy process.cowpy.memory = '2.GB'`
- [x] `#!groovy process { withName: 'cowpy' { memory = '2.GB' } }`
- [ ] `#!groovy resources.cowpy.memory = '2.GB'`

Aprende más: [5.3. Establecer asignaciones de recursos para un proceso específico](#53-set-resource-allocations-for-a-specific-process)
</quiz>

<quiz>
¿Qué hace la directiva `resourceLimits`?
- [ ] Establece requisitos mínimos de recursos
- [ ] Asigna recursos a procesos
- [x] Limita los recursos máximos que pueden solicitarse
- [ ] Monitorea el uso de recursos en tiempo real

Aprende más: [5.5. Agregar límites de recursos](#55-add-resource-limits)
</quiz>

<quiz>
¿Cómo especificas múltiples perfiles en un solo comando?
- [ ] `-profile perfil1 -profile perfil2`
- [ ] `-profiles perfil1,perfil2`
- [x] `-profile perfil1,perfil2`
- [ ] `--profile perfil1 --profile perfil2`

Aprende más: [6. Usar perfiles para cambiar entre configuraciones preestablecidas](#6-use-profiles-to-switch-between-preset-configurations)
</quiz>

<quiz>
¿Qué comando muestra la configuración completamente resuelta que Nextflow usaría?
- [ ] `nextflow show-config`
- [ ] `nextflow settings`
- [x] `nextflow config`
- [ ] `nextflow resolve`

Aprende más: [6.3. Usar `nextflow config` para ver la configuración resuelta](#63-use-nextflow-config-to-see-the-resolved-configuration)
</quiz>

<quiz>
¿Para qué se pueden usar los perfiles? (Selecciona todas las que apliquen)
- [x] Definir configuraciones específicas de infraestructura (ejecutores, contenedores)
- [x] Establecer límites de recursos para diferentes entornos
- [x] Proporcionar parámetros de prueba para facilitar las pruebas del workflow
- [ ] Definir nuevos procesos

Aprende más: [6. Usar perfiles para cambiar entre configuraciones preestablecidas](#6-use-profiles-to-switch-between-preset-configurations)
</quiz>
