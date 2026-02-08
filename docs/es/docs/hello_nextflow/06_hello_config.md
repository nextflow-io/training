# Parte 6: Hello Config

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Traducción asistida por IA - [más información y sugerencias](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

<div class="video-wrapper">
  <iframe width="560" height="315" src="https://www.youtube.com/embed/FcZTiE25TeA?si=tnXTi6mRkITY0zW_&amp;list=PLPZ8WHdZGxmWKozQuzr27jyMGqp9kElVK&amp;cc_load_policy=1&amp;cc_lang_pref=es" title="YouTube video player" frameborder="0" allow="accelerometer; autoplay; clipboard-write; encrypted-media; gyroscope; picture-in-picture; web-share" referrerpolicy="strict-origin-when-cross-origin" allowfullscreen></iframe>
</div>

/// caption
:fontawesome-brands-youtube:{ .youtube } Vea [la lista de reproducción completa](https://youtube.com/playlist?list=PLPZ8WHdZGxmWKozQuzr27jyMGqp9kElVK&si=eF7cLR62goy-lc6n) en el canal de YouTube de Nextflow.

:green_book: La transcripción del video está disponible [aquí](./transcripts/06_hello_config.md).
///

Esta sección explorará cómo configurar y gestionar la configuración de su pipeline de Nextflow para que pueda personalizar su comportamiento, adaptarlo a diferentes entornos y optimizar el uso de recursos _sin alterar una sola línea del código del flujo de trabajo en sí_.

Hay múltiples formas de hacer esto, que pueden usarse en combinación y se interpretan según el [orden de precedencia](https://www.nextflow.io/docs/latest/config.html) descrito en la documentación de configuración.

En esta parte del curso, vamos a mostrarle el mecanismo de archivo de configuración más simple y común, el archivo [`nextflow.config`](https://www.nextflow.io/docs/latest/config.html), que ya encontró en la Parte 5: Hello Containers.

Repasaremos los componentes esenciales de la configuración de Nextflow como directivas de proceso, executors, perfiles y archivos de parámetros.
Al aprender a utilizar estas opciones de configuración efectivamente, puede mejorar la flexibilidad, escalabilidad y rendimiento de sus pipelines.

??? info "Cómo comenzar desde esta sección"

    Esta sección del curso asume que ha completado las Partes 1-5 del curso [Hello Nextflow](./index.md) y tiene un pipeline completo funcionando.

    Si está comenzando el curso desde este punto, necesitará copiar el directorio `modules` y el archivo `nextflow.config` desde las soluciones:

    ```bash
    cp -r solutions/5-hello-containers/modules .
    cp solutions/5-hello-containers/nextflow.config .
    ```

    El archivo `nextflow.config` contiene la línea `docker.enabled = true` que habilita el uso de contenedores Docker.

    Si no está familiarizado con el pipeline Hello o podría usar un recordatorio, vea [esta página de información](../info/hello_pipeline.md).

---

## 0. Calentamiento: Ejecutar `hello-config.nf`

Vamos a usar el script de workflow `hello-config.nf` como punto de partida.
Es equivalente al script producido al trabajar en la Parte 5 de este curso de capacitación, excepto que hemos cambiado los destinos de salida:

```groovy title="hello-config.nf" linenums="37" hl_lines="3 7 11 15"
output {
    first_output {
        path 'hello_config/intermediates'
        mode 'copy'
    }
    uppercased {
        path 'hello_config/intermediates'
        mode 'copy'
    }
    collected {
        path 'hello_config/intermediates'
        mode 'copy'
    }
    batch_report {
        path 'hello_config'
        mode 'copy'
    }
    cowpy_art {
        path 'hello_config'
        mode 'copy'
    }
}
```

Solo para asegurarse de que todo funciona, ejecute el script una vez antes de hacer cualquier cambio:

```bash
nextflow run hello-config.nf
```

??? success "Salida del comando"

    ```console
     N E X T F L O W   ~  version 25.10.2

    Launching `hello-config.nf` [nice_escher] DSL2 - revision: d5dfdc9872

    executor > local (7)
    [6a/bc46a6] sayHello (2) [100%] 3 of 3 ✔
    [33/67bc48] convertToUpper (3) [100%] 3 of 3 ✔
    [b5/de03ba] collectGreetings [100%] 1 of 1 ✔
    [98/c6b57b] cowpy              | 1 of 1 ✔
    ```

Como anteriormente, encontrará los archivos de salida en el directorio especificado en el bloque `output` (`results/hello_config/`).

??? abstract "Contenido del directorio"

    ```console
    results/hello_config/
    ├── cowpy-COLLECTED-batch-output.txt
    ├── intermediates
    │   ├── Bonjour-output.txt
    │   ├── COLLECTED-batch-output.txt
    │   ├── Hello-output.txt
    │   ├── Holà-output.txt
    │   ├── UPPER-Bonjour-output.txt
    │   ├── UPPER-Hello-output.txt
    │   └── UPPER-Holà-output.txt
    └── batch-report.txt
    ```

La salida final de arte ASCII está en el directorio `results/hello_config/`, bajo el nombre `cowpy-COLLECTED-batch-output.txt`.

??? abstract "Contenido del archivo"

    ```console title="results/hello_config/cowpy-COLLECTED-batch-output.txt"
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

Si eso funcionó para usted, está listo para aprender cómo configurar sus pipelines.

---

## 1. Gestionar parámetros de entrada del flujo de trabajo

Vamos a comenzar con un aspecto de la configuración que es simplemente una extensión de lo que hemos estado trabajando hasta ahora: la gestión de parámetros de entrada.

Actualmente, nuestro flujo de trabajo está configurado para aceptar varios valores de parámetros a través de la línea de comandos, con valores predeterminados establecidos en un bloque `params` en el script del flujo de trabajo mismo.
Sin embargo, podría querer sobrescribir esos valores predeterminados sin tener que especificar parámetros en la línea de comandos ni modificar el archivo de script original.

Hay múltiples formas de hacer eso; vamos a mostrarle tres formas básicas que son muy comúnmente usadas.

### 1.1. Mover valores predeterminados a `nextflow.config`

Este es el enfoque más simple, aunque posiblemente sea el menos flexible ya que el archivo `nextflow.config` principal no es algo que quiera estar editando para cada ejecución.
Pero tiene la ventaja de separar las preocupaciones de _declarar_ los parámetros en el flujo de trabajo (que definitivamente pertenece allí) versus suministrar _valores predeterminados_, que están más en casa en un archivo de configuración.

Hagamos esto en dos pasos.

#### 1.1.1. Crear un bloque `params` en el archivo de configuración

Haga los siguientes cambios de código en el archivo `nextflow.config`:

=== "Después"

    ```groovy title="nextflow.config" linenums="1" hl_lines="3-10"
    docker.enabled = true

    /*
    * Parámetros del pipeline
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

Note que no simplemente copiamos el bloque `params` del flujo de trabajo al archivo de configuración.
La sintaxis es un poco diferente.
En el archivo de flujo de trabajo, esas son declaraciones tipadas.
En la configuración, esas son asignaciones de valores.

Técnicamente, esto es suficiente para sobrescribir los valores predeterminados aún especificados en el archivo de flujo de trabajo.
Podría modificar el personaje, por ejemplo, y ejecutar el flujo de trabajo para comprobar que el valor establecido en el archivo de configuración sobrescribe el establecido en el archivo de flujo de trabajo.

Pero en el espíritu de mover la configuración completamente al archivo de configuración, eliminemos esos valores del archivo de flujo de trabajo por completo.

#### 1.1.2. Eliminar los valores del bloque `params` en el archivo de flujo de trabajo

Haga los siguientes cambios de código en el archivo de flujo de trabajo `hello-config.nf`:

=== "Después"

    ```groovy title="hello-config.nf" linenums="9" hl_lines="5-7"
    /*
    * Parámetros del pipeline
    */
    params {
        input: Path
        batch: String
        character: String
    }
    ```

=== "Antes"

    ```groovy title="hello-config.nf" linenums="9" hl_lines="5-7"
    /*
    * Parámetros del pipeline
    */
    params {
        input: Path = 'data/greetings.csv'
        batch: String = 'batch'
        character: String = 'turkey'
    }
    ```

Ahora el archivo de flujo de trabajo mismo no establece ningún valor predeterminado para estos parámetros.

#### 1.1.3. Ejecutar el pipeline

Probemos que funciona correctamente.

```bash
nextflow run hello-config.nf
```

??? success "Salida del comando"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `hello-config.nf` [disturbed_einstein] DSL2 - revision: ede9037d02

    executor >  local (8)
    [f0/35723c] sayHello (2)       | 3 of 3 ✔
    [40/3efd1a] convertToUpper (3) | 3 of 3 ✔
    [17/e97d32] collectGreetings   | 1 of 1 ✔
    [98/c6b57b] cowpy              | 1 of 1 ✔
    ```

Esto todavía produce la misma salida que anteriormente.

La salida final de arte ASCII está en el directorio `results/hello_config/`, bajo el nombre `cowpy-COLLECTED-batch-output.txt`, igual que antes.

??? abstract "Contenido del archivo"

    ```console title="results/hello_config/cowpy-COLLECTED-batch-output.txt"
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

### 1.2. Usar un archivo de configuración específico de ejecución

Eso es genial, pero a veces podría querer ejecutar algunos experimentos temporales con diferentes valores predeterminados sin tocar el archivo de configuración principal.
Puede hacer eso creando un nuevo archivo `nextflow.config` en un subdirectorio que usará como directorio de trabajo para sus experimentos.

#### 1.2.1. Crear el directorio de trabajo con una configuración en blanco

Comencemos creando un nuevo directorio y entrando en él:

```bash
mkdir -p tux-run
cd tux-run
```

Luego, cree un archivo de configuración en blanco en ese directorio:

```bash
touch nextflow.config
```

Esto produce un archivo vacío.

#### 1.2.2. Configurar la configuración experimental

Ahora abra el nuevo archivo y agregue los parámetros que quiere personalizar:

```groovy title="tux-run/nextflow.config" linenums="1"
params {
    input = '../data/greetings.csv'
    batch = 'experiment'
    character = 'tux'
}
```

Note que la ruta al archivo de entrada debe reflejar la estructura de directorios.

#### 1.2.3. Ejecutar el pipeline

Ahora podemos ejecutar nuestro pipeline desde dentro de nuestro nuevo directorio de trabajo.
¡Asegúrese de adaptar la ruta en consecuencia!

```bash
nextflow run ../hello-config.nf
```

??? success "Salida del comando"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `../hello-config.nf` [trusting_escher] DSL2 - revision: 356df0818d

    executor >  local (8)
    [59/b66913] sayHello (2)       [100%] 3 of 3 ✔
    [ad/f06364] convertToUpper (3) [100%] 3 of 3 ✔
    [10/714895] collectGreetings   [100%] 1 of 1 ✔
    [88/3ece98] cowpy              [100%] 1 of 1 ✔
    ```

Esto creará un nuevo conjunto de directorios bajo `tux-run/` incluyendo `tux-run/work/` y `tux-run/results/`.

En esta ejecución, Nextflow combina el `nextflow.config` en nuestro directorio actual con el `nextflow.config` en el directorio raíz del pipeline, y así sobrescribe el personaje predeterminado (turkey) con el personaje tux.

El archivo de salida final debería contener el personaje tux diciendo los saludos.

??? abstract "Contenido del archivo"

    ```console title="tux-run/results/hello_config/cowpy-COLLECTED-experiment-output.txt"
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

Eso es todo; ahora tiene un espacio para experimentar sin modificar su configuración 'normal'.

!!! warning "Advertencia"

    ¡Asegúrese de volver al directorio anterior antes de pasar a la siguiente sección!

    ```bash
    cd ..
    ```

Ahora veamos otra forma útil de establecer valores de parámetros.

### 1.3. Usar un archivo de parámetros

El enfoque del subdirectorio funciona muy bien para experimentar, pero involucra un poco de configuración y requiere que adapte las rutas en consecuencia.
Hay un enfoque más simple para cuando quiere ejecutar su pipeline con un conjunto específico de valores, o permitir que alguien más lo haga con mínimo esfuerzo.

Nextflow nos permite especificar parámetros a través de un [archivo de parámetros](https://nextflow.io/docs/latest/config.html#params-file) en formato YAML o JSON, lo que hace muy conveniente gestionar y distribuir conjuntos alternativos de valores predeterminados, por ejemplo, así como valores de parámetros específicos de ejecución.

#### 1.3.1. Examinar el archivo de parámetros de ejemplo

Para demostrar esto, proporcionamos un archivo de parámetros de ejemplo en el directorio actual, llamado `test-params.yaml`:

```yaml title="test-params.yaml" linenums="1"
input: "data/greetings.csv"
batch: "yaml"
character: "stegosaurus"
```

Este archivo de parámetros contiene un par clave-valor para cada una de las entradas que queremos especificar.
Note el uso de dos puntos (`:`) en lugar de signos de igual (`=`) si compara la sintaxis con el archivo de configuración.
El archivo de configuración está escrito en Groovy, mientras que el archivo de parámetros está escrito en YAML.

!!! info "Información"

    También proporcionamos una versión JSON del archivo de parámetros como ejemplo pero no vamos a ejecutarla aquí.
    Siéntase libre de probar esa por su cuenta.

#### 1.3.2. Ejecutar el pipeline

Para ejecutar el flujo de trabajo con este archivo de parámetros, simplemente agregue `-params-file <filename>` al comando base.

```bash
nextflow run hello-config.nf -params-file test-params.yaml
```

??? success "Salida del comando"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `hello-config.nf` [disturbed_sammet] DSL2 - revision: ede9037d02

    executor >  local (8)
    [f0/35723c] sayHello (2)       | 3 of 3 ✔
    [40/3efd1a] convertToUpper (3) | 3 of 3 ✔
    [17/e97d32] collectGreetings   | 1 of 1 ✔
    [98/c6b57b] cowpy              | 1 of 1 ✔
    ```

El archivo de salida final debería contener el personaje stegosaurus diciendo los saludos.

??? abstract "Contenido del archivo"

    ```console title="results/hello_config/cowpy-COLLECTED-yaml-output.txt"
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

Usar un archivo de parámetros puede parecer excesivo cuando solo tiene unos pocos parámetros para especificar, pero algunos pipelines esperan docenas de parámetros.
En esos casos, usar un archivo de parámetros nos permitirá proporcionar valores de parámetros en tiempo de ejecución sin tener que escribir líneas de comando masivas y sin modificar el script de flujo de trabajo.

También hace más fácil distribuir conjuntos de parámetros a colaboradores, o como información de soporte para una publicación, por ejemplo.
Esto hace que su trabajo sea más reproducible por otros.

### Conclusión

Sabe cómo aprovechar las opciones de configuración clave para gestionar entradas de flujo de trabajo.

### ¿Qué sigue?

Aprender cómo gestionar dónde y cómo se publican las salidas de su flujo de trabajo.

---

## 2. Gestionar salidas del flujo de trabajo

Hasta ahora hemos estado codificando todas las rutas para las declaraciones de salida a nivel de flujo de trabajo, y como notamos cuando comenzamos a agregar múltiples salidas, puede haber un poco de repetición involucrada.

Veamos algunas formas comunes en que podría configurar esto para ser más flexible.

### 2.1. Personalizar el directorio de salida con `-output-dir`

Cuando estamos controlando cómo se organizan nuestras salidas 'publicadas', tenemos dos prioridades distintas:

- El directorio de salida de nivel superior
- Cómo se organizan los archivos dentro de este directorio

Hemos estado usando el directorio de nivel superior predeterminado hasta ahora: `results`.
Comencemos personalizando eso, usando la opción CLI `-output-dir`.

#### 2.1.1. Ejecutar el pipeline con `-output-dir`

La opción `-output-dir` (forma corta: `-o`) sobrescribe el directorio de salida predeterminado (`results/`) para todas las salidas del flujo de trabajo.
Esta es la forma recomendada de controlar la ruta raíz donde se publican las salidas.

```bash
nextflow run hello-config.nf -output-dir custom-outdir-cli/
```

??? success "Salida del comando"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `hello-config.nf` [prickly_kay] DSL2 - revision: 32ecc4fba2

    executor >  local (8)
    [9f/332636] sayHello (1)       [100%] 3 of 3 ✔
    [03/a55991] convertToUpper (3) [100%] 3 of 3 ✔
    [e5/ab7893] collectGreetings   [100%] 1 of 1 ✔
    [a8/97338e] cowpy              [100%] 1 of 1 ✔
    ```

Esto publica las salidas en `custom-outdir-cli/` en lugar de `results/`:

??? abstract "Contenido del directorio"

    ```console
    custom-outdir-cli/
    └── hello_config
        ├── batch-report.txt
        ├── cowpy-COLLECTED-batch-output.txt
        └── intermediates
            ├── Bonjour-output.txt
            ├── COLLECTED-batch-output.txt
            ├── Hello-output.txt
            ├── Holà-output.txt
            ├── UPPER-Bonjour-output.txt
            ├── UPPER-Hello-output.txt
            └── UPPER-Holà-output.txt
    ```

Note que todavía tenemos el subdirectorio `hello_config` de las declaraciones `path` en el bloque output.
Limpiemos eso.

#### 2.1.2. Eliminar rutas codificadas del bloque output

El prefijo `hello_config/` fue codificado en capítulos anteriores, pero dado que ahora estamos aprendiendo a configurar rutas de salida de manera flexible, podemos eliminar esta codificación.
Para salidas que no necesitan un subdirectorio podemos establecer la directiva `path` a una cadena vacía, o eliminarla por completo.

Haga los siguientes cambios de código en el archivo de flujo de trabajo:

=== "Después"

    ```groovy title="hello-config.nf" linenums="42" hl_lines="3 7 11 15 19"
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

    ```groovy title="hello-config.nf" linenums="42" hl_lines="3 7 11 15 19"
    output {
        first_output {
            path 'hello_config/intermediates'
            mode 'copy'
        }
        uppercased {
            path 'hello_config/intermediates'
            mode 'copy'
        }
        collected {
            path 'hello_config/intermediates'
            mode 'copy'
        }
        batch_report {
            path 'hello_config'
            mode 'copy'
        }
        cowpy_art {
            path 'hello_config'
            mode 'copy'
        }
    }
    ```

Ejecute el pipeline nuevamente:

```bash
nextflow run hello-config.nf -output-dir custom-outdir-cli-2/
```

Ahora las salidas se publican directamente bajo `custom-outdir-cli-2/`, sin el subdirectorio `hello_config`:

??? abstract "Contenido del directorio"

    ```console
    custom-outdir-cli-2/
    ├── batch-report.txt
    ├── cowpy-COLLECTED-batch-output.txt
    └── intermediates
        ├── Bonjour-output.txt
        ├── COLLECTED-batch-output.txt
        ├── Hello-output.txt
        ├── Holà-output.txt
        ├── UPPER-Bonjour-output.txt
        ├── UPPER-Hello-output.txt
        └── UPPER-Holà-output.txt
    ```

!!! tip "Consejo"

    La opción `-output-dir` se usa para controlar _dónde_ van las salidas, mientras que la directiva `path` en el bloque output controla la _estructura de subdirectorios_.

### 2.2. Rutas de salida dinámicas

Además de cambiar el directorio de salida a través del CLI, también podemos establecer un valor predeterminado personalizado en el archivo de configuración usando `outputDir`.
Esto nos permite establecer la ruta del directorio dinámicamente - no solo usando cadenas estáticas.

#### 2.2.1. Establecer `outputDir` en el archivo de configuración

Agregue el siguiente código al archivo `nextflow.config`:

=== "Después"

    ```groovy title="nextflow.config" linenums="9" hl_lines="10-13"
    /*
    * Parámetros del pipeline
    */
    params {
        input = 'data/greetings.csv'
        batch = 'batch'
        character = 'turkey'
    }

    /*
    * Configuración de salida
    */
    outputDir = "custom-outdir-config/${params.batch}"
    ```

=== "Antes"

    ```groovy title="nextflow.config" linenums="9"
    /*
    * Parámetros del pipeline
    */
    params {
        input = 'data/greetings.csv'
        batch = 'batch'
        character = 'turkey'
    }
    ```

Esto establece el directorio de salida a `custom-outdir-config/` más el valor del parámetro `batch` como subdirectorio.
Ahora puede cambiar la ubicación de salida estableciendo el parámetro `--batch`:

```bash
nextflow run hello-config.nf --batch my_run
```

Esto publica las salidas en `custom-outdir-config/my_run/`.

!!! note "Nota"

    La opción CLI `-output-dir` tiene precedencia sobre la configuración `outputDir`.
    Si se establece, la opción de configuración será ignorada por completo.

#### 2.2.2. Subdirectorios con nombres de batch y proceso

También podemos establecer declaraciones de `path` de salida de subdirectorio dinámicamente, por salida individual.

Por ejemplo, podemos organizar nuestras salidas por proceso referenciando `<process>.name` en la declaración de ruta de salida:

=== "Después"

    ```groovy title="hello-config.nf" linenums="42" hl_lines="3 7 11 15 19"
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

    ```groovy title="hello-config.nf" linenums="42" hl_lines="3 7 11 15 19"
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

Podemos ir más allá y componer rutas de subdirectorio más complejas.

En la edición anterior borramos la distinción entre `intermediates` versus salidas finales estando en el nivel superior.
Recuperemos eso, y también pongamos los archivos en un subdirectorio `params.batch`.

!!! tip "Consejo"

    Incluir `params.batch` en el `path` del bloque output, en lugar del `outputDir` config, significa que no será sobrescrito con `-output-dir` en el CLI.

Primero, actualice el archivo de configuración para eliminar `${params.batch}` de `outputDir` (ya que lo estamos moviendo a las declaraciones de ruta):

=== "Después"

    ```groovy title="nextflow.config" linenums="12" hl_lines="4"
    /*
    * Configuración de salida
    */
    outputDir = "custom-outdir-config/"
    ```

=== "Antes"

    ```groovy title="nextflow.config" linenums="12" hl_lines="4"
    /*
    * Configuración de salida
    */
    outputDir = "custom-outdir-config/${params.batch}"
    ```

Luego, haga los siguientes cambios en el archivo de flujo de trabajo:

=== "Después"

    ```groovy title="hello-config.nf" linenums="42" hl_lines="3 7 11 15 19"
    output {
        first_output {
            path { "${params.batch}/intermediates/${sayHello.name}" }
            mode 'copy'
        }
        uppercased {
            path { "${params.batch}/intermediates/${convertToUpper.name}" }
            mode 'copy'
        }
        collected {
            path { "${params.batch}/intermediates/${collectGreetings.name}" }
            mode 'copy'
        }
        batch_report {
            path { "${params.batch}/${collectGreetings.name}" }
            mode 'copy'
        }
        cowpy_art {
            path { "${params.batch}/${cowpy.name}" }
            mode 'copy'
        }
    }
    ```

=== "Antes"

    ```groovy title="hello-config.nf" linenums="42" hl_lines="3 7 11 15 19"
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

#### 2.2.3. Ejecutar el pipeline

Veamos cómo funciona esto en la práctica, estableciendo tanto `-output-dir` (o `-o` para abreviar) a `custom-outdir-config-2` como el nombre de batch a `rep2` desde la línea de comandos:

```bash
nextflow run hello-config.nf -output-dir custom-outdir-config-2 --batch rep2
```

??? success "Salida del comando"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `hello-config.nf` [mad_curry] DSL2 - revision: 668a98ccb9

    executor >  local (8)
    [9e/6095e0] sayHello (1)       [100%] 3 of 3 ✔
    [05/454d52] convertToUpper (3) [100%] 3 of 3 ✔
    [ed/e3ddfb] collectGreetings   [100%] 1 of 1 ✔
    [39/5e063a] cowpy              [100%] 1 of 1 ✔
    ```

Esto publica las salidas en `custom-outdir-config-2/rep2/`, con la ruta base especificada _y_ el subdirectorio del nombre de batch _y_ resultados agrupados por proceso:

??? abstract "Contenido del directorio"

    ```console
    custom-outdir-config-2
    └── rep2
        ├── collectGreetings
        │   └── rep2-report.txt
        ├── cowpy
        │   └── cowpy-COLLECTED-rep2-output.txt
        └── intermediates
            ├── collectGreetings
            │   └── COLLECTED-rep2-output.txt
            ├── convertToUpper
            │   ├── UPPER-Bonjour-output.txt
            │   ├── UPPER-Hello-output.txt
            │   └── UPPER-Holà-output.txt
            └── sayHello
                ├── Bonjour-output.txt
                ├── Hello-output.txt
                └── Holà-output.txt
    ```

### 2.3. Establecer el modo de publicación a nivel de flujo de trabajo

Finalmente, en el espíritu de reducir la cantidad de código repetitivo, podemos reemplazar las declaraciones `mode` por salida con una única línea en la configuración.

#### 2.3.1. Agregar `workflow.output.mode` al archivo de configuración

Agregue el siguiente código al archivo `nextflow.config`:

=== "Después"

    ```groovy title="nextflow.config" linenums="12" hl_lines="5"
    /*
    * Configuración de salida
    */
    outputDir = "custom-outdir-config/"
    workflow.output.mode = 'copy'
    ```

=== "Antes"

    ```groovy title="nextflow.config" linenums="12"
    /*
    * Configuración de salida
    */
    outputDir = "custom-outdir-config/"
    ```

Establecer `workflow.output.mode` en el archivo de configuración es suficiente para sobrescribir lo que está establecido en el archivo de flujo de trabajo, pero eliminemos el código innecesario de todos modos.

#### 2.3.2. Eliminar el modo de salida del archivo de flujo de trabajo

Haga los siguientes cambios en el archivo de flujo de trabajo:

=== "Después"

    ```groovy title="hello-config.nf" linenums="42"
    output {
        first_output {
            path { "${params.batch}/intermediates/${sayHello.name}" }
        }
        uppercased {
            path { "${params.batch}/intermediates/${convertToUpper.name}" }
        }
        collected {
            path { "${params.batch}/intermediates/${collectGreetings.name}" }
        }
        batch_report {
            path { "${params.batch}/${collectGreetings.name}" }
        }
        cowpy_art {
            path { "${params.batch}/${cowpy.name}" }
        }
    }
    ```

=== "Antes"

    ```groovy title="hello-config.nf" linenums="42" hl_lines="4 8 12 16 20"
    output {
        first_output {
            path { "${params.batch}/intermediates/${sayHello.name}" }
            mode 'copy'
        }
        uppercased {
            path { "${params.batch}/intermediates/${convertToUpper.name}" }
            mode 'copy'
        }
        collected {
            path { "${params.batch}/intermediates/${collectGreetings.name}" }
            mode 'copy'
        }
        batch_report {
            path { "${params.batch}/${collectGreetings.name}" }
            mode 'copy'
        }
        cowpy_art {
            path { "${params.batch}/${cowpy.name}" }
            mode 'copy'
        }
    }
    ```

Eso es más conciso, ¿no?

#### 2.3.3. Ejecutar el pipeline

Probemos que funciona correctamente:

```bash
nextflow run hello-config.nf -output-dir config-output-mode
```

??? success "Salida del comando"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `hello-config.nf` [small_stone] DSL2 - revision: 024d6361b5

    executor >  local (8)
    [e8/a0e93e] sayHello (1)       [100%] 3 of 3 ✔
    [14/176c9d] convertToUpper (3) [100%] 3 of 3 ✔
    [23/d667ca] collectGreetings   [100%] 1 of 1 ✔
    [e6/1dc80e] cowpy              [100%] 1 of 1 ✔
    ```

Esto publica las salidas en `config-output-mode/`, y todavía son todas copias apropiadas, no enlaces simbólicos.

??? abstract "Contenido del directorio"

    ```console
    config-output-mode
    └── batch
        ├── collectGreetings
        │   └── batch-report.txt
        ├── cowpy
        │   └── cowpy-COLLECTED-batch-output.txt
        └── intermediates
            ├── collectGreetings
            │   └── COLLECTED-batch-output.txt
            ├── convertToUpper
            │   ├── UPPER-Bonjour-output.txt
            │   ├── UPPER-Hello-output.txt
            │   └── UPPER-Holà-output.txt
            └── sayHello
                ├── Bonjour-output.txt
                ├── Hello-output.txt
                └── Holà-output.txt
    ```

La razón principal por la que aún podría querer usar la forma por salida de establecer el modo es si quiere mezclar y combinar dentro del mismo flujo de trabajo, _es decir_, tener algunas salidas copiadas y algunas enlazadas simbólicamente.

Hay muchas otras opciones que puede personalizar de esta manera, pero esperamos que esto le dé una idea del rango de opciones y cómo utilizarlas efectivamente para adaptarse a sus preferencias.

### Conclusión

Sabe cómo controlar el nombre y la estructura de los directorios donde se publican sus salidas, así como el modo de publicación de salida del flujo de trabajo.

### ¿Qué sigue?

Aprender cómo adaptar la configuración de su flujo de trabajo a su entorno de cómputo, comenzando con la tecnología de empaquetado de software.

---

## 3. Seleccionar una tecnología de empaquetado de software

Hasta ahora hemos estado viendo elementos de configuración que controlan cómo entran las entradas y de dónde salen las salidas. Ahora es tiempo de enfocarnos más específicamente en adaptar la configuración de su flujo de trabajo a su entorno de cómputo.

El primer paso en ese camino es especificar de dónde van a venir los paquetes de software que se ejecutarán en cada paso.
¿Ya están instalados en el entorno de cómputo local?
¿Necesitamos recuperar imágenes y ejecutarlas a través de un sistema de contenedores?
¿O necesitamos recuperar paquetes Conda y construir un entorno Conda local?

En la primera parte de este curso de capacitación (Partes 1-4) solo usamos software instalado localmente en nuestro flujo de trabajo.
Luego en la Parte 5, introdujimos contenedores Docker y el archivo `nextflow.config`, que usamos para habilitar el uso de contenedores Docker.

Ahora veamos cómo podemos configurar una opción alternativa de empaquetado de software a través del archivo `nextflow.config`.

### 3.1. Deshabilitar Docker y habilitar Conda en el archivo de configuración

Imaginemos que estamos trabajando en un clúster HPC y el administrador no permite el uso de Docker por razones de seguridad.
Afortunadamente para nosotros, Nextflow soporta múltiples otras tecnologías de contenedores incluyendo Singularity (que es más ampliamente usado en HPC), y gestores de paquetes de software como Conda.

Podemos cambiar nuestro archivo de configuración para usar [Conda](https://nextflow.io/docs/latest/conda.html) en lugar de Docker.
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
Lo que significa que ahora necesitamos agregar uno de esos a nuestro proceso `cowpy`!

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

!!! tip "Consejo"

    Hay algunas formas diferentes de obtener el URI para un paquete conda dado.
    Recomendamos usar la consulta de búsqueda de [Seqera Containers](https://seqera.io/containers/), que le dará un URI que puede copiar y pegar, incluso si no planea crear un contenedor a partir de él.

### 3.3. Ejecutar el flujo de trabajo para verificar que puede usar Conda

Probémoslo.

```bash
nextflow run hello-config.nf --batch conda
```

??? success "Salida del comando"

    ```console title="Output"
    N E X T F L O W   ~  version 25.10.2

    Launching `hello-config.nf` [friendly_lamport] DSL2 - revision: 024d6361b5

    executor >  local (8)
    [e8/91c116] sayHello (2)       [100%] 3 of 3 ✔
    [fe/6a70ce] convertToUpper (3) [100%] 3 of 3 ✔
    [99/7cc493] collectGreetings   [100%] 1 of 1 ✔
    [3c/09fb59] cowpy              [100%] 1 of 1 ✔
    ```

Esto debería funcionar sin problemas y producir las mismas salidas que anteriormente bajo `custom-outdir-config/conda`.

Detrás de escenas, Nextflow ha recuperado los paquetes Conda y creado el entorno, lo cual normalmente toma un poco de trabajo; ¡así que es bueno que no tengamos que hacer nada de eso nosotros mismos!

!!! note "Nota"

    Esto se ejecuta rápidamente porque el paquete `cowpy` es bastante pequeño, pero si está trabajando con paquetes grandes, puede tomar un poco más de lo usual la primera vez, y podría ver la salida de la consola quedarse 'atascada' por un minuto más o menos antes de completarse.
    Esto es normal y se debe al trabajo extra que Nextflow hace la primera vez que usa un nuevo paquete.

Desde nuestro punto de vista, parece que funciona exactamente igual que ejecutar con Docker, aunque en el backend la mecánica es un poco diferente.

Esto significa que estamos listos para ejecutar con entornos Conda si es necesario.

??? info "Mezclar y combinar Docker y Conda"

    Ya que estas directivas se asignan por proceso, es posible 'mezclar y combinar', _es decir_, configurar algunos de los procesos en su flujo de trabajo para ejecutar con Docker y otros con Conda, por ejemplo, si la infraestructura de cómputo que está usando soporta ambos.
    En ese caso, habilitaría tanto Docker como Conda en su archivo de configuración.
    Si ambos están disponibles para un proceso dado, Nextflow priorizará contenedores.

    Y como se señaló anteriormente, Nextflow soporta múltiples otras tecnologías de empaquetado de software y contenedores, así que no está limitado a solo esas dos.

### Conclusión

Sabe cómo configurar qué paquete de software debería usar cada proceso, y cómo cambiar entre tecnologías.

### ¿Qué sigue?

Aprender cómo cambiar la plataforma de ejecución usada por Nextflow para realmente hacer el trabajo.

---

## 4. Seleccionar una plataforma de ejecución

Hasta ahora, hemos estado ejecutando nuestro pipeline con el executor local.
Esto ejecuta cada tarea en la máquina donde Nextflow está corriendo.
Cuando Nextflow comienza, mira los CPUs y memoria disponibles.
Si los recursos de las tareas listas para ejecutar exceden los recursos disponibles, Nextflow retendrá las últimas tareas de la ejecución hasta que una o más de las tareas anteriores hayan terminado, liberando los recursos necesarios.

El executor local es conveniente y eficiente, pero está limitado a esa única máquina. Para cargas de trabajo muy grandes, puede descubrir que su máquina local es un cuello de botella, ya sea porque tiene una única tarea que requiere más recursos de los que tiene disponibles, o porque tiene tantas tareas que esperar a que una sola máquina las ejecute tomaría demasiado tiempo.

Nextflow soporta [muchos executors diferentes](https://www.nextflow.io/docs/latest/executor.html), incluyendo programadores HPC (Slurm, LSF, SGE, PBS, Moab, OAR, Bridge, HTCondor y otros) así como backends de ejecución en la nube (AWS Batch, Google Cloud Batch, Azure Batch, Kubernetes y más).

### 4.1. Apuntar a un backend diferente

La elección del executor se establece por una directiva de proceso llamada `executor`.
Por defecto está establecido a `local`, así que la siguiente configuración está implícita:

```groovy title="Configuración incorporada"
process {
    executor = 'local'
}
```

Para establecer el executor para apuntar a un backend diferente, simplemente especificaría el executor que quiere usando sintaxis similar a la descrita arriba para asignaciones de recursos (vea la [documentación de executors](https://www.nextflow.io/docs/latest/executor.html) para todas las opciones).

```groovy title="nextflow.config"
process {
    executor = 'slurm'
}
```

!!! warning "Advertencia"

    En realidad no podemos probar esto en el entorno de capacitación porque no está configurado para conectarse a un HPC.

### 4.2. Lidiar con sintaxis específica del backend para parámetros de ejecución

La mayoría de las plataformas de computación de alto rendimiento permiten (y a veces requieren) que especifique ciertos parámetros como solicitudes y limitaciones de asignación de recursos (por ej. número de CPUs y memoria) y nombre de la cola de trabajos a usar.

Desafortunadamente, cada uno de estos sistemas usa diferentes tecnologías, sintaxis y configuraciones para definir cómo un trabajo debería ser definido y enviado al programador relevante.

??? abstract "Ejemplos"

    Por ejemplo, el mismo trabajo que requiere 8 CPUs y 4GB de RAM para ser ejecutado en la cola "my-science-work" necesita ser expresado de diferentes maneras dependiendo del backend.

    ```bash title="Config para SLURM / enviar usando sbatch"
    #SBATCH -o /path/to/my/task/directory/my-task-1.log
    #SBATCH --no-requeue
    #SBATCH -c 8
    #SBATCH --mem 4096M
    #SBATCH -p my-science-work
    ```

    ```bash title="Config para PBS / enviar usando qsub"
    #PBS -o /path/to/my/task/directory/my-task-1.log
    #PBS -j oe
    #PBS -q my-science-work
    #PBS -l nodes=1:ppn=5
    #PBS -l mem=4gb
    ```

    ```bash title="Config para SGE / enviar usando qsub"
    #$ -o /path/to/my/task/directory/my-task-1.log
    #$ -j y
    #$ -terse
    #$ -notify
    #$ -q my-science-work
    #$ -l slots=5
    #$ -l h_rss=4096M,mem_free=4096M
    ```

Afortunadamente, Nextflow simplifica todo esto.
Proporciona una sintaxis estandarizada para que pueda especificar las propiedades relevantes como [`cpus`](https://nextflow.io/docs/latest/reference/process.html#cpus), [`memory`](https://nextflow.io/docs/latest/reference/process.html#memory) y [`queue`](https://nextflow.io/docs/latest/reference/process.html#queue) (vea las [directivas de proceso](https://nextflow.io/docs/latest/reference/process.html#process-directives) para otras propiedades) solo una vez.
Luego, en tiempo de ejecución, Nextflow usará esas configuraciones para generar los scripts específicos del backend apropiados basados en la configuración del executor.

Cubriremos esa sintaxis estandarizada en la siguiente sección.

### Conclusión

Ahora sabe cómo cambiar el executor para usar diferentes tipos de infraestructura de cómputo.

### ¿Qué sigue?

Aprender cómo evaluar y expresar asignaciones y limitaciones de recursos en Nextflow.

---

## 5. Controlar asignaciones de recursos de cómputo

La mayoría de las plataformas de computación de alto rendimiento permiten (y a veces requieren) que especifique ciertos parámetros de asignación de recursos como número de CPUs y memoria.

Por defecto, Nextflow usará un único CPU y 2GB de memoria para cada proceso.
Las directivas de proceso correspondientes se llaman `cpus` y `memory`, así que la siguiente configuración está implícita:

```groovy title="Configuración incorporada" linenums="1"
process {
    cpus = 1
    memory = 2.GB
}
```

Puede modificar estos valores, ya sea para todos los procesos o para procesos nombrados específicos, usando directivas de proceso adicionales en su archivo de configuración.
Nextflow las traducirá a las instrucciones apropiadas para el executor elegido.

¿Pero cómo sabe qué valores usar?

### 5.1. Ejecutar el flujo de trabajo para generar un reporte de utilización de recursos

Si no sabe de antemano cuánto CPU y memoria es probable que necesiten sus procesos, puede hacer algo de perfilado de recursos, lo que significa que ejecuta el flujo de trabajo con algunas asignaciones predeterminadas, registra cuánto usó cada proceso, y de ahí, estima cómo ajustar las asignaciones base.

Convenientemente, Nextflow incluye herramientas incorporadas para hacer esto, y felizmente generará un reporte para usted cuando lo solicite.

Para hacerlo, agregue `-with-report <filename>.html` a su línea de comandos.

```bash
nextflow run hello-config.nf -with-report report-config-1.html
```

El reporte es un archivo html, que puede descargar y abrir en su navegador. También puede hacer clic derecho en él en el explorador de archivos a la izquierda y hacer clic en `Show preview` para verlo en el entorno de capacitación.

Tómese unos minutos para revisar el reporte y ver si puede identificar algunas oportunidades para ajustar recursos.
Asegúrese de hacer clic en las pestañas que muestran los resultados de utilización como porcentaje de lo que fue asignado.

Vea [Reports](https://www.nextflow.io/docs/latest/reports.html) para documentación sobre todas las características disponibles.

### 5.2. Establecer asignaciones de recursos para todos los procesos

El perfilado muestra que los procesos en nuestro flujo de trabajo de capacitación son muy ligeros, así que reduzcamos la asignación de memoria predeterminada a 1GB por proceso.

Agregue lo siguiente a su archivo `nextflow.config`, antes de la sección de parámetros del pipeline:

=== "Después"

    ```groovy title="nextflow.config" linenums="1" hl_lines="4-9"
    docker.enabled = false
    conda.enabled = true

    /*
    * Configuración de procesos
    */
    process {
        memory = 1.GB
    }

    /*
    * Parámetros del pipeline
    */
    params {
        input = 'data/greetings.csv'
        batch = 'batch'
        character = 'turkey'
    }
    ```

=== "Antes"

    ```groovy title="nextflow.config" linenums="1"
    docker.enabled = false
    conda.enabled = true

    /*
    * Parámetros del pipeline
    */
    params {
        input = 'data/greetings.csv'
        batch = 'batch'
        character = 'turkey'
    }
    ```

Eso ayudará a reducir la cantidad de cómputo que consumimos.

### 5.3. Establecer asignaciones de recursos para un proceso específico

Al mismo tiempo, vamos a pretender que el proceso `cowpy` requiere más recursos que los otros, solo para que podamos demostrar cómo ajustar asignaciones para un proceso individual.

=== "Después"

    ```groovy title="nextflow.config" linenums="4" hl_lines="6-9"
    /*
    * Configuración de procesos
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
    * Configuración de procesos
    */
    process {
        memory = 1.GB
    }
    ```

Con esta configuración, todos los procesos solicitarán 1GB de memoria y un único CPU (el predeterminado implícito), excepto el proceso `cowpy`, que solicitará 2GB y 2 CPUs.

!!! tip "Consejo"

    Si tiene una máquina con pocos CPUs y asigna un número alto por proceso, podría ver llamadas de proceso siendo encoladas detrás de otras.
    Esto es porque Nextflow asegura que no solicitemos más CPUs de los que están disponibles.

### 5.4. Ejecutar el flujo de trabajo con la configuración actualizada

Probemos eso, suministrando un nombre de archivo diferente para el reporte de perfilado para que podamos comparar el rendimiento antes y después de los cambios de configuración.

```bash
nextflow run hello-config.nf -with-report report-config-2.html
```

Probablemente no notará ninguna diferencia real ya que esta es una carga de trabajo tan pequeña, pero este es el enfoque que usaría para analizar el rendimiento y los requisitos de recursos de un flujo de trabajo del mundo real.

Es muy útil cuando sus procesos tienen diferentes requisitos de recursos. Le permite dimensionar correctamente las asignaciones de recursos que configura para cada proceso basándose en datos reales, no en suposiciones.

!!! tip "Consejo"

    Esto es solo un pequeño adelanto de lo que puede hacer para optimizar su uso de recursos.
    Nextflow mismo tiene una [lógica de reintento dinámica](https://www.nextflow.io/docs/latest/process.html#dynamic-task-resources) realmente elegante incorporada para reintentar trabajos que fallan debido a limitaciones de recursos.
    Adicionalmente, la Plataforma Seqera ofrece herramientas impulsadas por IA para optimizar sus asignaciones de recursos automáticamente también.

### 5.5. Agregar límites de recursos

Dependiendo de qué executor de cómputo e infraestructura de cómputo esté usando, puede haber algunas restricciones sobre lo que puede (o debe) asignar.
Por ejemplo, su clúster puede requerir que permanezca dentro de ciertos límites.

Puede usar la directiva `resourceLimits` para establecer las limitaciones relevantes. La sintaxis se ve así cuando está por sí sola en un bloque process:

```groovy title="Ejemplo de sintaxis"
process {
    resourceLimits = [
        memory: 750.GB,
        cpus: 200,
        time: 30.d
    ]
}
```

Nextflow traducirá estos valores a las instrucciones apropiadas dependiendo del executor que especificó.

No vamos a ejecutar esto, ya que no tenemos acceso a infraestructura relevante en el entorno de capacitación.
Sin embargo, si intentara ejecutar el flujo de trabajo con asignaciones de recursos que excedan estos límites, luego buscara el comando `sbatch` en el archivo de script `.command.run`, vería que las solicitudes que realmente se envían al executor están limitadas a los valores especificados por `resourceLimits`.

??? info "Configuraciones de referencia institucionales"

    El proyecto nf-core ha compilado una [colección de archivos de configuración](https://nf-co.re/configs/) compartidos por varias instituciones alrededor del mundo, cubriendo una amplia gama de executors HPC y de nube.

    Esas configuraciones compartidas son valiosas tanto para las personas que trabajan allí y por lo tanto pueden simplemente utilizar la configuración de su institución directamente, como un modelo para personas que buscan desarrollar una configuración para su propia infraestructura.

### Conclusión

Sabe cómo generar un reporte de perfilado para evaluar la utilización de recursos y cómo modificar las asignaciones de recursos para todos los procesos y/o para procesos individuales, así como establecer limitaciones de recursos para ejecutar en HPC.

### ¿Qué sigue?

Aprender cómo configurar perfiles de configuración preestablecidos y cambiar entre ellos en tiempo de ejecución.

---

## 6. Usar perfiles para cambiar entre configuraciones preestablecidas

Le hemos mostrado varias formas en que puede personalizar la configuración de su pipeline dependiendo del proyecto en el que está trabajando o el entorno de cómputo que está usando.

Puede querer cambiar entre configuraciones alternativas dependiendo de qué infraestructura de cómputo está usando. Por ejemplo, podría querer desarrollar y ejecutar pruebas de pequeña escala localmente en su laptop, luego ejecutar cargas de trabajo a escala completa en HPC o nube.

Nextflow le permite configurar cualquier número de [perfiles](https://nextflow.io/docs/latest/config.html#config-profiles) que describen diferentes configuraciones, que luego puede seleccionar en tiempo de ejecución usando un argumento de línea de comandos, en lugar de tener que modificar el archivo de configuración mismo.

### 6.1. Crear perfiles para cambiar entre desarrollo local y ejecución en HPC

Configuremos dos perfiles alternativos; uno para ejecutar cargas pequeñas en una computadora regular, donde usaremos contenedores Docker, y uno para ejecutar en un HPC universitario con un programador Slurm, donde usaremos paquetes Conda.

#### 6.1.1. Configurar los perfiles

Agregue lo siguiente a su archivo `nextflow.config`, después de la sección de parámetros del pipeline pero antes de la configuración de salida:

=== "Después"

    ```groovy title="nextflow.config" linenums="15" hl_lines="10-27"
    /*
    * Parámetros del pipeline
    */
    params {
        input = 'data/greetings.csv'
        batch = 'batch'
        character = 'turkey'
    }

    /*
    * Perfiles
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

    /*
    * Configuración de salida
    */
    outputDir = "custom-outdir-config/"
    workflow.output.mode = 'copy'
    ```

=== "Antes"

    ```groovy title="nextflow.config" linenums="15"
    /*
    * Parámetros del pipeline
    */
    params {
        input = 'data/greetings.csv'
        batch = 'batch'
        character = 'turkey'
    }

    /*
    * Configuración de salida
    */
    outputDir = "custom-outdir-config/"
    workflow.output.mode = 'copy'
    ```

Verá que para el HPC universitario, también estamos especificando limitaciones de recursos.

#### 6.1.2. Ejecutar el flujo de trabajo con un perfil

Para especificar un perfil en nuestra línea de comandos de Nextflow, usamos el argumento `-profile`.

Intentemos ejecutar el flujo de trabajo con la configuración `my_laptop`.

```bash
nextflow run hello-config.nf -profile my_laptop
```

??? success "Salida del comando"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `hello-config.nf` [hungry_sanger] DSL2 - revision: 024d6361b5

    executor >  local (8)
    [b0/fb2ec9] sayHello (3)       [100%] 3 of 3 ✔
    [4a/e039f0] convertToUpper (3) [100%] 3 of 3 ✔
    [6f/408fa9] collectGreetings   [100%] 1 of 1 ✔
    [f1/fd6520] cowpy              [100%] 1 of 1 ✔
    ```

Como puede ver, esto nos permite alternar entre configuraciones muy convenientemente en tiempo de ejecución.

!!! warning "Advertencia"

    El perfil `univ_hpc` no se ejecutará correctamente en el entorno de capacitación ya que no tenemos acceso a un programador Slurm.

Si en el futuro encontramos otros elementos de configuración que siempre co-ocurren con estos, simplemente podemos agregarlos al(los) perfil(es) correspondiente(s).
También podemos crear perfiles adicionales si hay otros elementos de configuración que queramos agrupar.

### 6.2. Crear un perfil de parámetros de prueba

Los perfiles no son solo para configuración de infraestructura.
También podemos usarlos para establecer valores predeterminados para parámetros de flujo de trabajo, para hacer más fácil para otros probar el flujo de trabajo sin tener que reunir valores de entrada apropiados ellos mismos.
Puede considerar esto una alternativa a usar un archivo de parámetros.

#### 6.2.1. Configurar el perfil

La sintaxis para expresar valores predeterminados en este contexto se ve así, para un perfil que nombramos `test`:

```groovy title="Ejemplo de sintaxis"
    test {
        params.<parameter1>
        params.<parameter2>
        ...
    }
```

Si agregamos un perfil de prueba para nuestro flujo de trabajo, el bloque `profiles` se convierte en:

```groovy title="nextflow.config" linenums="24" hl_lines="18-22"
/*
* Perfiles
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

Al igual que para los perfiles de configuración técnica, puede configurar múltiples perfiles diferentes especificando parámetros bajo cualquier nombre arbitrario que desee.

#### 6.2.2. Ejecutar el flujo de trabajo localmente con el perfil de prueba

Convenientemente, los perfiles no son mutuamente excluyentes, así que podemos especificar múltiples perfiles en nuestra línea de comandos usando la siguiente sintaxis `-profile <profile1>,<profile2>` (para cualquier número de perfiles).

Si combina perfiles que establecen valores para los mismos elementos de configuración y están descritos en el mismo archivo de configuración, Nextflow resolverá el conflicto usando cualquier valor que haya leído último (_es decir_, lo que viene después en el archivo).
Si las configuraciones en conflicto están establecidas en diferentes fuentes de configuración, se aplica el [orden de precedencia](https://www.nextflow.io/docs/latest/config.html) predeterminado.

Intentemos agregar el perfil de prueba a nuestro comando anterior:

```bash
nextflow run hello-config.nf -profile my_laptop,test
```

??? success "Salida del comando"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `hello-config.nf` [modest_becquerel] DSL2 - revision: 024d6361b5

    executor >  local (8)
    [4c/fe2580] sayHello (1)       [100%] 3 of 3 ✔
    [fd/7d9017] convertToUpper (3) [100%] 3 of 3 ✔
    [13/1523bd] collectGreetings   [100%] 1 of 1 ✔
    [06/a1ee14] cowpy              [100%] 1 of 1 ✔
    ```

Esto usará Docker donde sea posible y producirá salidas bajo `custom-outdir-config/test`, y esta vez el personaje es el dúo cómico `dragonandcow`.

??? abstract "Contenido del archivo"

    ```console title="custom-outdir-config/test/cowpy/cowpy-COLLECTED-test-output.txt"
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

Esto significa que mientras distribuyamos cualquier archivo de datos de prueba con el código del flujo de trabajo, cualquiera puede probar rápidamente el flujo de trabajo sin tener que suministrar sus propias entradas a través de la línea de comandos o un archivo de parámetros.

!!! tip "Consejo"

    Podemos apuntar a URLs para archivos más grandes que están almacenados externamente.
    Nextflow los descargará automáticamente mientras haya una conexión abierta.

    Para más detalles, vea la Misión Secundaria [Trabajando con Archivos](../side_quests/working_with_files.md)

### 6.3. Usar `nextflow config` para ver la configuración resuelta

Como se señaló arriba, a veces el mismo parámetro puede establecerse a diferentes valores en perfiles que quiere combinar.
Y más generalmente, hay numerosos lugares donde los elementos de configuración pueden almacenarse, y a veces las mismas propiedades pueden establecerse a diferentes valores en diferentes lugares.

Nextflow aplica un [orden de precedencia](https://www.nextflow.io/docs/latest/config.html) establecido para resolver cualquier conflicto, pero eso puede ser complicado de determinar usted mismo.
E incluso si nada está en conflicto, puede ser tedioso buscar todos los lugares posibles donde las cosas podrían estar configuradas.

Afortunadamente, Nextflow incluye una herramienta de utilidad conveniente llamada `config` que puede automatizar todo ese proceso por usted.

La herramienta `config` explorará todos los contenidos en su directorio de trabajo actual, aspirará cualquier archivo de configuración, y producirá la configuración completamente resuelta que Nextflow usaría para ejecutar el flujo de trabajo.
Esto le permite averiguar qué configuraciones se usarán sin tener que lanzar nada.

#### 6.3.1. Resolver la configuración predeterminada

Ejecute este comando para resolver la configuración que se aplicaría por defecto.

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

    outputDir = 'custom-outdir-config/'

    workflow {
      output {
          mode = 'copy'
      }
    }
    ```

Esto le muestra la configuración base que obtiene si no especifica nada extra en la línea de comandos.

#### 6.3.2. Resolver la configuración con configuraciones específicas activadas

Si proporciona parámetros de línea de comandos, ej. habilitando uno o más perfiles o cargando un archivo de parámetros, el comando adicionalmente tomará esos en cuenta.

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

    outputDir = 'custom-outdir-config/'

    workflow {
      output {
          mode = 'copy'
      }
    }
    ```

Esto se vuelve especialmente útil para proyectos complejos que involucran múltiples capas de configuración.

### Conclusión

Sabe cómo usar perfiles para seleccionar una configuración preestablecida en tiempo de ejecución con mínimo esfuerzo.
Más generalmente, sabe cómo configurar las ejecuciones de su flujo de trabajo para adaptarse a diferentes plataformas de cómputo y mejorar la reproducibilidad de sus análisis.

### ¿Qué sigue?

¡Celebre y dese una gran palmada en la espalda! Ha completado su primer curso de desarrollador de Nextflow.

Diríjase al [resumen final del curso](./next_steps.md) para revisar lo que aprendió y descubrir qué viene después.

---

## Cuestionario

<quiz>
¿Cuál es el nombre del archivo de configuración que Nextflow carga automáticamente?
- [ ] `config.nf`
- [ ] `pipeline.config`
- [x] `nextflow.config`
- [ ] `workflow.config`
</quiz>

<quiz>
¿Qué tiene precedencia cuando el mismo parámetro está establecido tanto en el archivo de configuración como en la línea de comandos?
- [ ] El valor del archivo de configuración
- [x] El valor de la línea de comandos
- [ ] El primer valor encontrado
- [ ] Ninguno; causa un error

Aprenda más: [1.1. Mover valores predeterminados a `nextflow.config`](#11-mover-valores-predeterminados-a-nextflowconfig)
</quiz>

<quiz>
¿Puede tener tanto Docker como Conda habilitados en la misma configuración?
- [x] Sí, Nextflow puede usar ambos dependiendo de las directivas del proceso
- [ ] No, solo uno puede estar habilitado a la vez
- [ ] Sí, pero solo en perfiles
- [ ] No, son mutuamente excluyentes
</quiz>

<quiz>
Si tanto Docker como Conda están habilitados y un proceso tiene ambas directivas, ¿cuál se prioriza?
- [x] Docker (contenedores)
- [ ] Conda
- [ ] El primero definido
- [ ] Causa un error

Aprenda más: [3. Seleccionar una tecnología de empaquetado de software](#3-seleccionar-una-tecnologia-de-empaquetado-de-software)
</quiz>

<quiz>
¿Cuál es la asignación de memoria predeterminada para los procesos de Nextflow?
- [ ] 1 GB
- [x] 2 GB
- [ ] 4 GB
- [ ] Sin límite
</quiz>

<quiz>
¿Cómo establece los requisitos de recursos para un proceso específico en el archivo de configuración?
- [ ] `#!groovy processName.memory = '4 GB'`
- [ ] `#!groovy process.memory.processName = '4 GB'`
- [x] `#!groovy process { withName: 'processName' { memory = '4 GB' } }`
- [ ] `#!groovy resources.processName.memory = '4 GB'`

Aprenda más: [5.3. Establecer asignaciones de recursos para un proceso específico](#53-establecer-asignaciones-de-recursos-para-un-proceso-especifico)
</quiz>

<quiz>
¿Qué opción de línea de comandos genera un reporte de utilización de recursos?
- [ ] `-with-metrics`
- [ ] `-with-stats`
- [x] `-with-report`
- [ ] `-with-profile`

Aprenda más: [5.1. Ejecutar el flujo de trabajo para generar un reporte de utilización de recursos](#51-ejecutar-el-flujo-de-trabajo-para-generar-un-reporte-de-utilizacion-de-recursos)
</quiz>

<quiz>
¿Qué hace la directiva `resourceLimits`?
- [ ] Establece requisitos mínimos de recursos
- [ ] Asigna recursos a procesos
- [x] Limita los recursos máximos que pueden solicitarse
- [ ] Monitorea el uso de recursos

Aprenda más: [5.5. Agregar límites de recursos](#55-agregar-limites-de-recursos)
</quiz>

<quiz>
¿Cuál es el executor predeterminado en Nextflow?
- [x] `local`
- [ ] `slurm`
- [ ] `kubernetes`
- [ ] `aws`

Aprenda más: [4. Seleccionar una plataforma de ejecución](#4-seleccionar-una-plataforma-de-ejecucion)
</quiz>

<quiz>
¿Cómo especifica un archivo de parámetros al ejecutar Nextflow?
- [ ] `--params params.json`
- [ ] `-config params.json`
- [x] `-params-file params.json`
- [ ] `--input params.json`

Aprenda más: [1.3. Usar un archivo de parámetros](#13-usar-un-archivo-de-parametros)
</quiz>

<quiz>
¿Para qué pueden usarse los perfiles? (Seleccione todos los que apliquen)
- [x] Definir configuraciones específicas de infraestructura
- [x] Establecer límites de recursos para diferentes entornos
- [x] Proporcionar parámetros de prueba
- [ ] Definir nuevos procesos

Aprenda más: [6. Usar perfiles para cambiar entre configuraciones preestablecidas](#6-usar-perfiles-para-cambiar-entre-configuraciones-preestablecidas)
</quiz>

<quiz>
¿Cómo especifica múltiples perfiles en un único comando?
- [ ] `-profile profile1 -profile profile2`
- [ ] `-profiles profile1,profile2`
- [x] `-profile profile1,profile2`
- [ ] `--profile profile1 --profile profile2`

Aprenda más: [6. Usar perfiles para cambiar entre configuraciones preestablecidas](#6-usar-perfiles-para-cambiar-entre-configuraciones-preestablecidas)
</quiz>
