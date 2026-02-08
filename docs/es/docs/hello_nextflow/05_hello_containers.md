# Parte 5: Hello Containers

<div class="video-wrapper">
  <iframe width="560" height="315" src="https://www.youtube.com/embed/Xqr--bKEN9U?si=QinuAnFwFj-Z8CrO&amp;list=PLPZ8WHdZGxmWKozQuzr27jyMGqp9kElVK&amp;cc_load_policy=1&amp;cc_lang_pref=es" title="YouTube video player" frameborder="0" allow="accelerometer; autoplay; clipboard-write; encrypted-media; gyroscope; picture-in-picture; web-share" referrerpolicy="strict-origin-when-cross-origin" allowfullscreen></iframe>
</div>

/// caption
:fontawesome-brands-youtube:{ .youtube } Vea [la lista de reproducción completa](https://youtube.com/playlist?list=PLPZ8WHdZGxmWKozQuzr27jyMGqp9kElVK&si=eF7cLR62goy-lc6n) en el canal de YouTube de Nextflow.

:green_book: La transcripción del video está disponible [aquí](./transcripts/05_hello_containers.md).
///

En las Partes 1-4 de este curso de capacitación, aprendió cómo usar los bloques de construcción básicos de Nextflow para ensamblar un workflow simple capaz de procesar algo de texto, paralelizar la ejecución si había múltiples entradas, y recopilar los resultados para procesamiento adicional.

Sin embargo, estaba limitado a herramientas UNIX básicas disponibles en su entorno.
Las tareas del mundo real a menudo requieren varias herramientas y paquetes que no están incluidos por defecto.
Típicamente, necesitaría instalar estas herramientas, gestionar sus dependencias y resolver cualquier conflicto.

Todo eso es muy tedioso y molesto, así que vamos a mostrarle cómo usar **contenedores** para resolver este problema de manera mucho más conveniente.

Un **contenedor** es una unidad de software ligera, independiente y ejecutable creada a partir de una **imagen** de contenedor que incluye todo lo necesario para ejecutar una aplicación incluyendo código, bibliotecas del sistema y configuraciones.
Como puede imaginar, eso va a ser muy útil para hacer sus pipelines más reproducibles.

Note que enseñaremos esto usando [Docker](https://www.docker.com/get-started/), pero tenga en cuenta que Nextflow soporta [varias otras tecnologías de contenedores](https://nextflow.io/docs/latest/container.html) también.

??? info "Cómo comenzar desde esta sección"

    Esta sección del curso asume que ha completado las Partes 1-4 del curso [Hello Nextflow](./index.md) y tiene un pipeline completo funcionando.

    Si está comenzando el curso desde este punto, necesitará copiar el directorio `modules` desde las soluciones:

    ```bash
    cp -r solutions/4-hello-modules/modules .
    ```

---

## 0. Calentamiento: Ejecutar `hello-containers.nf`

Vamos a usar el script de workflow `hello-containers.nf` como punto de partida.
Es equivalente al script producido al trabajar en la Parte 4 de este curso de capacitación, excepto que hemos cambiado los destinos de salida:

```groovy title="hello-containers.nf" linenums="37" hl_lines="3 7 11 15"
output {
    first_output {
        path 'hello_containers'
        mode 'copy'
    }
    uppercased {
        path 'hello_containers'
        mode 'copy'
    }
    collected {
        path 'hello_containers'
        mode 'copy'
    }
    batch_report {
        path 'hello_containers'
        mode 'copy'
    }
}
```

Solo para asegurarse de que todo funciona, ejecute el script una vez antes de hacer cualquier cambio:

```bash
nextflow run hello-containers.nf
```

??? success "Salida del comando"

    ```console
     N E X T F L O W   ~  version 25.10.2

    Launching `hello-containers.nf` [nice_escher] DSL2 - revision: d5dfdc9872

    executor > local (7)
    [5a/ec1fa1] sayHello (2) [100%] 3 of 3 ✔
    [30/32b5b8] convertToUpper (3) [100%] 3 of 3 ✔
    [d3/be01bc] collectGreetings [100%] 1 of 1 ✔

    ```

Como anteriormente, encontrará los archivos de salida en el directorio especificado en el bloque `output` (`results/hello_containers/`).

??? abstract "Contenido del directorio"

    ```console
    results/hello_containers/
    ├── Bonjour-output.txt
    ├── COLLECTED-batch-output.txt
    ├── Hello-output.txt
    ├── Holà-output.txt
    ├── batch-report.txt
    ├── UPPER-Bonjour-output.txt
    ├── UPPER-Hello-output.txt
    └── UPPER-Holà-output.txt
    ```

Si eso funcionó para usted, está listo para aprender cómo usar contenedores.

---

## 1. Usar un contenedor 'manualmente'

Lo que queremos hacer es agregar un paso a nuestro workflow que usará un contenedor para la ejecución.

Sin embargo, primero vamos a repasar algunos conceptos y operaciones básicas para solidificar su comprensión de qué son los contenedores antes de comenzar a usarlos en Nextflow.

### 1.1. Descargar la imagen del contenedor

Para usar un contenedor, usualmente descarga o _pull_ una imagen de contenedor de un registro de contenedores, y luego ejecuta la imagen del contenedor para crear una instancia de contenedor.

La sintaxis general es la siguiente:

```bash title="Syntax"
docker pull '<container>'
```

La parte `docker pull` es la instrucción al sistema de contenedores para descargar una imagen de contenedor de un repositorio.

La parte `'<container>'` es la dirección URI de la imagen del contenedor.

Como ejemplo, descarguemos una imagen de contenedor que contiene [cowpy](https://github.com/jeffbuttars/cowpy), una implementación en Python de una herramienta llamada `cowsay` que genera arte ASCII para mostrar entradas de texto arbitrarias de una manera divertida.

```txt title="Example"
 ________________________
< Are we having fun yet? >
 ------------------------
    \                                  ___-------___
     \                             _-~~             ~~-_
      \                         _-~                    /~-_
             /^\__/^\         /~  \                   /    \
           /|  O|| O|        /      \_______________/        \
          | |___||__|      /       /                \          \
          |          \    /      /                    \          \
          |   (_______) /______/                        \_________ \
          |         / /         \                      /            \
           \         \^\\         \                  /               \     /
             \         ||           \______________/      _-_       //\__//
               \       ||------_-~~-_ ------------- \ --/~   ~\    || __/
                 ~-----||====/~     |==================|       |/~~~~~
                  (_(__/  ./     /                    \_\      \.
                         (_(___/                         \_____)_)
```

Hay varios repositorios donde puede encontrar contenedores publicados.
Usamos el servicio [Seqera Containers](https://seqera.io/containers/) para generar esta imagen de contenedor Docker desde el paquete Conda `cowpy`: `'community.wave.seqera.io/library/cowpy:1.1.5--3db457ae1977a273'`.

Ejecute el comando de descarga completo:

```bash
docker pull 'community.wave.seqera.io/library/cowpy:1.1.5--3db457ae1977a273'
```

??? success "Salida del comando"

    ```console
    1.1.5--3db457ae1977a273: Pulling from library/cowpy
    dafa2b0c44d2: Pull complete
    dec6b097362e: Pull complete
    f88da01cff0b: Pull complete
    4f4fb700ef54: Pull complete
    92dc97a3ef36: Pull complete
    403f74b0f85e: Pull complete
    10b8c00c10a5: Pull complete
    17dc7ea432cc: Pull complete
    bb36d6c3110d: Pull complete
    0ea1a16bbe82: Pull complete
    030a47592a0a: Pull complete
    c23bdb422167: Pull complete
    e1686ff32a11: Pull complete
    Digest: sha256:1ebc0043e8cafa61203bf42d29fd05bd14e7b4298e5e8cf986504c15f5aa4160
    Status: Downloaded newer image for community.wave.seqera.io/library/cowpy:1.1.5--3db457ae1977a273
    community.wave.seqera.io/library/cowpy:1.1.5--3db457ae1977a273
    ```

Si nunca ha descargado la imagen antes, esto puede tardar un minuto en completarse.
Una vez que esté hecho, tiene una copia local de la imagen del contenedor.

### 1.2. Usar el contenedor para ejecutar `cowpy` como un comando único

Una forma muy común en que las personas usan contenedores es ejecutarlos directamente, _es decir_, no interactivamente.
Esto es genial para ejecutar comandos únicos.

La sintaxis general es la siguiente:

```bash title="Syntax"
docker run --rm '<container>' [tool command]
```

La parte `docker run --rm '<container>'` es la instrucción al sistema de contenedores para iniciar una instancia de contenedor desde una imagen de contenedor y ejecutar un comando en ella.
La bandera `--rm` le dice al sistema que apague la instancia del contenedor después de que el comando se haya completado.

La sintaxis `[tool command]` depende de la herramienta que esté usando y cómo esté configurado el contenedor.
Comencemos simplemente con `cowpy`.

Completamente ensamblado, el comando de ejecución del contenedor se ve así; adelante y ejecútelo.

```bash
docker run --rm 'community.wave.seqera.io/library/cowpy:1.1.5--3db457ae1977a273' cowpy
```

??? success "Salida del comando"

    ```console
    ______________________________________________________
    < Cowacter, eyes:default, tongue:False, thoughts:False >
    ------------------------------------------------------
        \   ^__^
          \  (oo)\_______
            (__)\       )\/\
              ||----w |
              ||     ||
    ```

El sistema inició el contenedor, ejecutó el comando `cowpy` con sus parámetros, envió la salida a la consola y finalmente, apagó la instancia del contenedor.

### 1.3. Usar el contenedor para ejecutar `cowpy` interactivamente

También puede ejecutar un contenedor interactivamente, lo que le da un prompt de shell dentro del contenedor y le permite jugar con el comando.

#### 1.3.1. Iniciar el contenedor

Para ejecutar interactivamente, solo agregamos `-it` al comando `docker run`.
Opcionalmente, podemos especificar el shell que queremos usar dentro del contenedor agregando _ej._ `/bin/bash` al comando.

```bash
docker run --rm -it 'community.wave.seqera.io/library/cowpy:1.1.5--3db457ae1977a273' /bin/bash
```

Note que su prompt cambia a algo como `(base) root@b645838b3314:/tmp#`, lo que indica que ahora está dentro del contenedor.

Puede verificar esto ejecutando `ls /` para listar el contenido del directorio desde la raíz del sistema de archivos:

```bash
ls /
```

??? abstract "Salida del comando"

    ```console
    bin  boot  dev  etc  home  lib  lib64  media  mnt  opt  proc  root  run  sbin  srv  sys  tmp  usr  var
    ```

Usamos `ls` aquí en lugar de `tree` porque la utilidad `tree` no está disponible en este contenedor.
Puede ver que el sistema de archivos dentro del contenedor es diferente del sistema de archivos en su sistema host.

Una limitación de lo que acabamos de hacer es que el contenedor está completamente aislado del sistema host por defecto.
Esto significa que el contenedor no puede acceder a ningún archivo en el sistema host a menos que explícitamente se lo permita.

Le mostraremos cómo hacer eso en un minuto.

#### 1.3.2. Ejecutar el/los comando(s) de la herramienta deseada

Ahora que está dentro del contenedor, puede ejecutar el comando `cowpy` directamente y darle algunos parámetros.
Por ejemplo, la documentación de la herramienta dice que podemos cambiar el personaje ('cowacter') con `-c`.

```bash
cowpy "Hello Containers" -c tux
```

??? success "Salida del comando"

    ```console
    __________________
    < Hello Containers >
    ------------------
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

Ahora la salida muestra el pingüino de Linux, Tux, en lugar de la vaca predeterminada, porque especificamos el parámetro `-c tux`.

Como está dentro del contenedor, puede ejecutar el comando `cowpy` tantas veces como quiera, variando los parámetros de entrada, sin tener que molestarse con los comandos de Docker.

!!! Tip "Consejo"

    Use la bandera '-c' para elegir un personaje diferente, incluyendo:
    `beavis`, `cheese`, `daemon`, `dragonandcow`, `ghostbusters`, `kitty`, `moose`, `milk`, `stegosaurus`, `turkey`, `turtle`, `tux`

Esto es genial. Lo que sería aún más genial es si pudiéramos alimentar nuestro `greetings.csv` como entrada en esto.
Pero como no tenemos acceso al sistema de archivos, no podemos.

Arreglemos eso.

#### 1.3.3. Salir del contenedor

Para salir del contenedor, puede escribir `exit` en el prompt o usar el atajo de teclado ++ctrl+d++.

```bash
exit
```

Su prompt ahora debería estar de vuelta a lo que era antes de que iniciara el contenedor.

#### 1.3.4. Montar datos en el contenedor

Como se señaló anteriormente, el contenedor está aislado del sistema host por defecto.

Para permitir que el contenedor acceda al sistema de archivos del host, puede **montar** un **volumen** desde el sistema host en el contenedor usando la siguiente sintaxis:

```bash title="Syntax"
-v <outside_path>:<inside_path>
```

En nuestro caso `<outside_path>` será el directorio de trabajo actual, así que podemos simplemente usar un punto (`.`), y `<inside_path>` es solo un alias que inventamos; llamémoslo `/my_project` (la ruta interna debe ser absoluta).

Para montar un volumen, reemplazamos las rutas y agregamos el argumento de montaje de volumen al comando docker run de la siguiente manera:

```bash
docker run --rm -it -v .:/my_project 'community.wave.seqera.io/library/cowpy:1.1.5--3db457ae1977a273' /bin/bash
```

Esto monta el directorio de trabajo actual como un volumen que será accesible bajo `/my_project` dentro del contenedor.

Puede verificar que funciona listando el contenido de `/my_project`:

```bash
ls /my_project
```

??? success "Salida del comando"

    ```console
    data               hello-config.nf      hello-modules.nf   hello-world.nf  nextflow.config  solutions         work
    hello-channels.nf  hello-containers.nf  hello-workflow.nf  modules         results          test-params.json
    ```

Ahora puede ver el contenido del directorio de trabajo desde dentro del contenedor, incluyendo el archivo `greetings.csv` bajo `data/`.

Esto efectivamente estableció un túnel a través de la pared del contenedor que puede usar para acceder a esa parte de su sistema de archivos.

#### 1.3.5. Usar los datos montados

Ahora que hemos montado el directorio de trabajo en el contenedor, podemos usar el comando `cowpy` para mostrar el contenido del archivo `greetings.csv`.

Para hacer esto, usaremos `cat /my_project/data/greetings.csv | ` para canalizar el contenido del archivo CSV al comando `cowpy`.

```bash
cat /my_project/data/greetings.csv | cowpy -c turkey
```

??? success "Salida del comando"

    ```console title="data/greetings.csv"
     ____________________
    / Hello,English,123  \
    | Bonjour,French,456 |
    \ Holà,Spanish,789   /
    --------------------
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

¡Esto produce el arte ASCII deseado de un pavo recitando nuestros saludos de ejemplo!
Excepto que aquí el pavo está repitiendo las filas completas en lugar de solo los saludos.
¡Ya sabemos que nuestro workflow de Nextflow hará un mejor trabajo!

Siéntase libre de jugar con este comando.
Cuando haya terminado, salga del contenedor como anteriormente:

```bash
exit
```

Se encontrará de vuelta en su shell normal.

### Conclusión

Sabe cómo descargar un contenedor y ejecutarlo ya sea como un comando único o interactivamente. También sabe cómo hacer que sus datos sean accesibles desde dentro de su contenedor, lo que le permite probar cualquier herramienta que le interese en datos reales sin tener que instalar ningún software en su sistema.

### ¿Qué sigue?

Aprender cómo usar contenedores para la ejecución de procesos de Nextflow.

---

## 2. Usar contenedores en Nextflow

Nextflow tiene soporte integrado para ejecutar procesos dentro de contenedores para permitirle ejecutar herramientas que no tiene instaladas en su entorno de cómputo.
Esto significa que puede usar cualquier imagen de contenedor que desee para ejecutar sus procesos, y Nextflow se encargará de descargar la imagen, montar los datos y ejecutar el proceso dentro de ella.

Para demostrar esto, vamos a agregar un paso `cowpy` al pipeline que hemos estado desarrollando, después del paso `collectGreetings`.

<figure class="excalidraw">
--8<-- "docs/en/docs/hello_nextflow/img/hello-pipeline-cowpy.svg"
</figure>

### 2.1. Escribir un módulo `cowpy`

Primero, creemos el módulo del proceso `cowpy`.

#### 2.1.1. Crear un archivo stub para el nuevo módulo

Cree un archivo vacío para el módulo llamado `cowpy.nf`.

```bash
touch modules/cowpy.nf
```

Esto nos da un lugar para poner el código del proceso.

#### 2.1.2. Copiar el código del proceso `cowpy` en el archivo del módulo

Podemos modelar nuestro proceso `cowpy` en los otros procesos que hemos escrito anteriormente.

```groovy title="modules/cowpy.nf" linenums="1"
#!/usr/bin/env nextflow

// Generate ASCII art with cowpy
process cowpy {

    input:
    path input_file
    val character

    output:
    path "cowpy-${input_file}"

    script:
    """
    cat ${input_file} | cowpy -c "${character}" > cowpy-${input_file}
    """

}
```

El proceso espera un `input_file` que contiene los saludos así como un valor `character`.

La salida será un nuevo archivo de texto que contiene el arte ASCII generado por la herramienta `cowpy`.

### 2.2. Agregar cowpy al workflow

Ahora necesitamos importar el módulo y llamar al proceso.

#### 2.2.1. Importar el proceso `cowpy` en `hello-containers.nf`

Inserte la declaración de importación arriba del bloque workflow y complétela apropiadamente.

=== "Después"

    ```groovy title="hello-containers.nf" linenums="3" hl_lines="5"
    // Include modules
    include { sayHello } from './modules/sayHello.nf'
    include { convertToUpper } from './modules/convertToUpper.nf'
    include { collectGreetings } from './modules/collectGreetings.nf'
    include { cowpy } from './modules/cowpy.nf'
    ```

=== "Antes"

    ```groovy title="hello-containers.nf" linenums="3"
    // Include modules
    include { sayHello } from './modules/sayHello.nf'
    include { convertToUpper } from './modules/convertToUpper.nf'
    include { collectGreetings } from './modules/collectGreetings.nf'
    ```

Ahora el módulo `cowpy` está disponible para usar en el workflow.

#### 2.2.2. Agregar una llamada al proceso `cowpy` en el workflow

Conectemos el proceso `cowpy()` a la salida del proceso `collectGreetings()`, que como recordará produce dos salidas:

- `collectGreetings.out.outfile` contiene el archivo de salida <--_lo que queremos_
- `collectGreetings.out.report` contiene el archivo de reporte con el conteo de saludos por lote

En el bloque workflow, haga el siguiente cambio de código:

=== "Después"

    ```groovy title="hello-containers.nf" linenums="19" hl_lines="12-13"
        main:
        // create a channel for inputs from a CSV file
        greeting_ch = channel.fromPath(params.input)
                            .splitCsv()
                            .map { line -> line[0] }
        // emit a greeting
        sayHello(greeting_ch)
        // convert the greeting to uppercase
        convertToUpper(sayHello.out)
        // collect all the greetings into one file
        collectGreetings(convertToUpper.out.collect(), params.batch)
        // generate ASCII art of the greetings with cowpy
        cowpy(collectGreetings.out.outfile, params.character)
    ```

=== "Antes"

    ```groovy title="hello-containers.nf" linenums="19"
        main:
        // create a channel for inputs from a CSV file
        greeting_ch = channel.fromPath(params.input)
                            .splitCsv()
                            .map { line -> line[0] }
        // emit a greeting
        sayHello(greeting_ch)
        // convert the greeting to uppercase
        convertToUpper(sayHello.out)
        // collect all the greetings into one file
        collectGreetings(convertToUpper.out.collect(), params.batch)
    ```

Note que declaramos un nuevo parámetro CLI, `params.character`, para especificar qué personaje queremos que diga los saludos.

#### 2.2.3. Agregar el parámetro `character` al bloque `params`

Esto es técnicamente opcional pero es la práctica recomendada y es una oportunidad para establecer un valor predeterminado para el personaje mientras estamos en ello.

=== "Después"

    ```groovy title="hello-containers.nf" linenums="9" hl_lines="7"
    /*
    * Parámetros del pipeline
    */
    params {
        input: Path = 'data/greetings.csv'
        batch: String = 'batch'
        character: String = 'turkey'
    }
    ```

=== "Antes"

    ```groovy title="hello-containers.nf" linenums="9"
    /*
    * Parámetros del pipeline
    */
    params {
        input: Path = 'data/greetings.csv'
        batch: String = 'batch'
    }
    ```

Ahora podemos ser perezosos y omitir escribir el parámetro character en nuestras líneas de comando.

#### 2.2.4. Actualizar las salidas del workflow

Necesitamos actualizar las salidas del workflow para publicar la salida del proceso `cowpy`.

##### 2.2.4.1. Actualizar la sección `publish:`

En el bloque `workflow`, haga el siguiente cambio de código:

=== "Después"

    ```groovy title="hello-containers.nf" linenums="34" hl_lines="6"
        publish:
        first_output = sayHello.out
        uppercased = convertToUpper.out
        collected = collectGreetings.out.outfile
        batch_report = collectGreetings.out.report
        cowpy_art = cowpy.out
    ```

=== "Antes"

    ```groovy title="hello-containers.nf" linenums="34"
        publish:
        first_output = sayHello.out
        uppercased = convertToUpper.out
        collected = collectGreetings.out.outfile
        batch_report = collectGreetings.out.report
    ```

El proceso `cowpy` solo produce una salida así que podemos referirnos a ella de la manera usual agregando `.out`.

Pero por ahora, terminemos de actualizar las salidas a nivel de workflow.

##### 2.2.4.2. Actualizar el bloque `output`

Necesitamos agregar la salida final `cowpy_art` al bloque `output`. Mientras estamos en ello, también editemos los destinos de publicación ya que ahora nuestro pipeline está completo y sabemos qué salidas realmente nos importan.

En el bloque `output`, haga los siguientes cambios de código:

=== "Después"

    ```groovy title="hello-containers.nf" linenums="42" hl_lines="3 7 11 15 18-21"
    output {
        first_output {
            path 'hello_containers/intermediates'
            mode 'copy'
        }
        uppercased {
            path 'hello_containers/intermediates'
            mode 'copy'
        }
        collected {
            path 'hello_containers/intermediates'
            mode 'copy'
        }
        batch_report {
            path 'hello_containers'
            mode 'copy'
        }
        cowpy_art {
            path 'hello_containers'
            mode 'copy'
        }
    }
    ```

=== "Antes"

    ```groovy title="hello-containers.nf" linenums="42" hl_lines="3 7 11 15"
    output {
        first_output {
            path 'hello_containers'
            mode 'copy'
        }
        uppercased {
            path 'hello_containers'
            mode 'copy'
        }
        collected {
            path 'hello_containers'
            mode 'copy'
        }
        batch_report {
            path 'hello_containers'
            mode 'copy'
        }
    }
    ```

Ahora las salidas publicadas estarán un poco más organizadas.

#### 2.2.5. Ejecutar el workflow

Solo para recapitular, esto es lo que estamos buscando:

<figure class="excalidraw">
--8<-- "docs/en/docs/hello_nextflow/img/hello_pipeline_complete.svg"
</figure>

¿Cree que va a funcionar?

Eliminemos las salidas publicadas anteriores para tener una pizarra limpia, y ejecutemos el workflow con la bandera `-resume`.

```bash
rm -r hello_containers/
nextflow run hello-containers.nf -resume
```

??? failure "Salida del comando (editada para claridad)"

    ```console hl_lines="10 13 20-21 26-27"
     N E X T F L O W   ~  version 25.10.2

    Launching `hello-containers.nf` [lonely_woese] DSL2 - revision: abf1dccf7f

    executor >  local (1)
    [c9/f5c686] sayHello (3)       [100%] 3 of 3, cached: 3 ✔
    [ef/3135a8] convertToUpper (3) [100%] 3 of 3, cached: 3 ✔
    [7f/f435e3] collectGreetings   [100%] 1 of 1, cached: 1 ✔
    [9b/02e776] cowpy              [  0%] 0 of 1 ✘
    ERROR ~ Error executing process > 'cowpy'

    Caused by:
      Process `cowpy` terminated with an error exit status (127)


    Command executed:

      cat COLLECTED-batch-output.txt | cowpy -c "turkey" > cowpy-COLLECTED-batch-output.txt

    Command exit status:
      127

    Command output:
      (empty)

    Command error:
      .command.sh: line 2: cowpy: command not found

    Work dir:
      /workspaces/training/hello-nextflow/work/9b/02e7761db848f82db3c3e59ff3a9b6

    Tip: when you have fixed the problem you can continue the execution adding the option `-resume` to the run command line

    -- Check '.nextflow.log' file for details
    ERROR ~ Cannot access first() element from an empty List

    -- Check '.nextflow.log' file for details
    ```

¡Oh no, hay un error!
El código de error dado por `error exit status (127)` significa que el ejecutable que pedimos no fue encontrado.

Eso tiene sentido, ya que estamos llamando a la herramienta `cowpy` pero en realidad no hemos especificado un contenedor todavía (ups).

### 2.3. Usar un contenedor para ejecutar el proceso `cowpy`

Necesitamos especificar un contenedor y decirle a Nextflow que lo use para el proceso `cowpy()`.

#### 2.3.1. Especificar un contenedor para `cowpy`

Podemos usar la misma imagen que estábamos usando directamente en la primera sección de este tutorial.

Edite el módulo `cowpy.nf` para agregar la directiva `container` a la definición del proceso de la siguiente manera:

=== "Después"

    ```groovy title="modules/cowpy.nf" linenums="4" hl_lines="3"
    process cowpy {

        container 'community.wave.seqera.io/library/cowpy:1.1.5--3db457ae1977a273'

        input:
        path input_file
        val character

        output:
        path "cowpy-${input_file}"

        script:
        """
        cat ${input_file} | cowpy -c "${character}" > cowpy-${input_file}
        """
    }
    ```

=== "Antes"

    ```groovy title="modules/cowpy.nf" linenums="4"
    process cowpy {

        input:
        path input_file
        val character

        output:
        path "cowpy-${input_file}"

        script:
        """
        cat ${input_file} | cowpy -c "${character}" > cowpy-${input_file}
        """
    }
    ```

Esto le dice a Nextflow que _si el uso de Docker está habilitado_, debería usar la imagen de contenedor especificada aquí para ejecutar el proceso.

#### 2.3.2. Habilitar el uso de Docker a través del archivo `nextflow.config`

Note que dijimos _'si el uso de Docker está habilitado'_. Por defecto, no lo está, así que necesitamos decirle a Nextflow que está permitido usar Docker.
Para ese fin, vamos a anticipar ligeramente el tema de la siguiente y última parte de este curso (Parte 6), que cubre la configuración.

Una de las principales formas que ofrece Nextflow para configurar la ejecución del workflow es usar un archivo `nextflow.config`.
Cuando tal archivo está presente en el directorio actual, Nextflow lo cargará automáticamente y aplicará cualquier configuración que contenga.

Proporcionamos un archivo `nextflow.config` con una única línea de código que explícitamente deshabilita Docker: `docker.enabled = false`.

Ahora, cambiemos eso a `true` para habilitar Docker:

=== "Después"

    ```console title="nextflow.config" linenums="1" hl_lines="1"
    docker.enabled = true
    ```

=== "Antes"

    ```console title="nextflow.config" linenums="1" hl_lines="1"
    docker.enabled = false
    ```

!!! tip "Consejo"

    Es posible habilitar la ejecución de Docker desde la línea de comandos, por ejecución, usando el parámetro `-with-docker <container>`.
    Sin embargo, eso solo nos permite especificar un contenedor para todo el workflow, mientras que el enfoque que acabamos de mostrarle nos permite especificar un contenedor diferente por proceso.
    Esto es mejor para la modularidad, el mantenimiento del código y la reproducibilidad.

#### 2.3.3. Ejecutar el workflow con Docker habilitado

Ejecute el workflow con la bandera `-resume`:

```bash
nextflow run hello-containers.nf -resume
```

??? success "Salida del comando"

    ```console
     N E X T F L O W   ~  version 25.10.2

    Launching `hello-containers.nf` [drunk_perlman] DSL2 - revision: abf1dccf7f

    executor >  local (1)
    [c9/f5c686] sayHello (3)       [100%] 3 of 3, cached: 3 ✔
    [ef/3135a8] convertToUpper (3) [100%] 3 of 3, cached: 3 ✔
    [7f/f435e3] collectGreetings   [100%] 1 of 1, cached: 1 ✔
    [98/656c6c] cowpy              [100%] 1 of 1 ✔
    ```

¡Esta vez sí funciona!
Como siempre puede encontrar las salidas del workflow en el directorio de resultados correspondiente, aunque esta vez están un poco más organizadas, con solo el reporte y la salida final en el nivel superior, y todos los archivos intermedios apartados en un subdirectorio.

??? abstract "Contenido del directorio"

    ```console
    results/hello_containers/
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

La salida final de arte ASCII está en el directorio `results/hello_containers/`, bajo el nombre `cowpy-COLLECTED-batch-output.txt`.

??? abstract "Contenido del archivo"

    ```console title="results/hello_containers/cowpy-COLLECTED-batch-output.txt"
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

Y ahí está, nuestro hermoso pavo diciendo los saludos como se deseaba.

#### 2.3.4. Inspeccionar cómo Nextflow lanzó la tarea containerizada

Como coda final a esta sección, echemos un vistazo al subdirectorio de trabajo para una de las llamadas del proceso `cowpy` para obtener un poco más de información sobre cómo Nextflow trabaja con contenedores bajo el capó.

Verifique la salida de su comando `nextflow run` para encontrar la ruta al subdirectorio de trabajo para el proceso `cowpy`.
Mirando lo que obtuvimos para la ejecución mostrada arriba, la línea del log de consola para el proceso `cowpy` comienza con `[98/656c6c]`.
Eso corresponde a la siguiente ruta de directorio truncada: `work/98/656c6c`.

En ese directorio, encontrará el archivo `.command.run` que contiene todos los comandos que Nextflow ejecutó en su nombre durante el curso de la ejecución del pipeline.

??? abstract "Contenido del archivo"

    ```console title="work/98/656c6c90cce1667c094d880f4b6dcc/.command.run"
    #!/bin/bash
    ### ---
    ### name: 'cowpy'
    ### container: 'community.wave.seqera.io/library/cowpy:1.1.5--3db457ae1977a273'
    ### outputs:
    ### - 'cowpy-COLLECTED-batch-output.txt'
    ### ...
    set -e
    set -u
    NXF_DEBUG=${NXF_DEBUG:=0}; [[ $NXF_DEBUG > 1 ]] && set -x
    NXF_ENTRY=${1:-nxf_main}


    nxf_sleep() {
      sleep $1 2>/dev/null || sleep 1;
    }

    nxf_date() {
        local ts=$(date +%s%3N);
        if [[ ${#ts} == 10 ]]; then echo ${ts}000
        elif [[ $ts == *%3N ]]; then echo ${ts/\%3N/000}
        elif [[ $ts == *3N ]]; then echo ${ts/3N/000}
        elif [[ ${#ts} == 13 ]]; then echo $ts
        else echo "Unexpected timestamp value: $ts"; exit 1
        fi
    }

    nxf_env() {
        echo '============= task environment ============='
        env | sort | sed "s/\(.*\)AWS\(.*\)=\(.\{6\}\).*/\1AWS\2=\3xxxxxxxxxxxxx/"
        echo '============= task output =================='
    }

    nxf_kill() {
        declare -a children
        while read P PP;do
            children[$PP]+=" $P"
        done < <(ps -e -o pid= -o ppid=)

        kill_all() {
            [[ $1 != $$ ]] && kill $1 2>/dev/null || true
            for i in ${children[$1]:=}; do kill_all $i; done
        }

        kill_all $1
    }

    nxf_mktemp() {
        local base=${1:-/tmp}
        mkdir -p "$base"
        if [[ $(uname) = Darwin ]]; then mktemp -d $base/nxf.XXXXXXXXXX
        else TMPDIR="$base" mktemp -d -t nxf.XXXXXXXXXX
        fi
    }

    nxf_fs_copy() {
      local source=$1
      local target=$2
      local basedir=$(dirname $1)
      mkdir -p $target/$basedir
      cp -fRL $source $target/$basedir
    }

    nxf_fs_move() {
      local source=$1
      local target=$2
      local basedir=$(dirname $1)
      mkdir -p $target/$basedir
      mv -f $source $target/$basedir
    }

    nxf_fs_rsync() {
      rsync -rRl $1 $2
    }

    nxf_fs_rclone() {
      rclone copyto $1 $2/$1
    }

    nxf_fs_fcp() {
      fcp $1 $2/$1
    }

    on_exit() {
        local last_err=$?
        local exit_status=${nxf_main_ret:=0}
        [[ ${exit_status} -eq 0 && ${nxf_unstage_ret:=0} -ne 0 ]] && exit_status=${nxf_unstage_ret:=0}
        [[ ${exit_status} -eq 0 && ${last_err} -ne 0 ]] && exit_status=${last_err}
        printf -- $exit_status > /workspaces/training/hello-nextflow/work/98/656c6c90cce1667c094d880f4b6dcc/.exitcode
        set +u
        docker rm $NXF_BOXID &>/dev/null || true
        exit $exit_status
    }

    on_term() {
        set +e
        docker stop $NXF_BOXID
    }

    nxf_launch() {
        docker run -i --cpu-shares 1024 -e "NXF_TASK_WORKDIR" -v /workspaces/training/hello-nextflow/work:/workspaces/training/hello-nextflow/work -w "$NXF_TASK_WORKDIR" --name $NXF_BOXID community.wave.seqera.io/library/cowpy:1.1.5--3db457ae1977a273 /bin/bash -ue /workspaces/training/hello-nextflow/work/98/656c6c90cce1667c094d880f4b6dcc/.command.sh
    }

    nxf_stage() {
        true
        # stage input files
        rm -f COLLECTED-batch-output.txt
        ln -s /workspaces/training/hello-nextflow/work/7f/f435e3f2cf95979b5f3d7647ae6696/COLLECTED-batch-output.txt COLLECTED-batch-output.txt
    }

    nxf_unstage_outputs() {
        true
    }

    nxf_unstage_controls() {
        true
    }

    nxf_unstage() {
        if [[ ${nxf_main_ret:=0} == 0 ]]; then
            (set -e -o pipefail; (nxf_unstage_outputs | tee -a .command.out) 3>&1 1>&2 2>&3 | tee -a .command.err)
            nxf_unstage_ret=$?
        fi
        nxf_unstage_controls
    }

    nxf_main() {
        trap on_exit EXIT
        trap on_term TERM INT USR2
        trap '' USR1

        [[ "${NXF_CHDIR:-}" ]] && cd "$NXF_CHDIR"
        export NXF_BOXID="nxf-$(dd bs=18 count=1 if=/dev/urandom 2>/dev/null | base64 | tr +/ 0A | tr -d '\r\n')"
        NXF_SCRATCH=''
        [[ $NXF_DEBUG > 0 ]] && nxf_env
        touch /workspaces/training/hello-nextflow/work/98/656c6c90cce1667c094d880f4b6dcc/.command.begin
        set +u
        set -u
        [[ $NXF_SCRATCH ]] && cd $NXF_SCRATCH
        export NXF_TASK_WORKDIR="$PWD"
        nxf_stage

        set +e
        (set -o pipefail; (nxf_launch | tee .command.out) 3>&1 1>&2 2>&3 | tee .command.err) &
        pid=$!
        wait $pid || nxf_main_ret=$?
        nxf_unstage
    }

    $NXF_ENTRY

    ```

Si busca `nxf_launch` en este archivo, debería ver algo como esto:

```console
nxf_launch() {
    docker run -i --cpu-shares 1024 -e "NXF_TASK_WORKDIR" -v /workspaces/training/hello-nextflow/work:/workspaces/training/hello-nextflow/work -w "$NXF_TASK_WORKDIR" --name $NXF_BOXID community.wave.seqera.io/library/cowpy:1.1.5--3db457ae1977a273 /bin/bash -ue /workspaces/training/hello-nextflow/work/98/656c6c90cce1667c094d880f4b6dcc/.command.sh
}
```

Como puede ver, Nextflow está usando el comando `docker run` para lanzar la llamada del proceso.
También monta el subdirectorio de trabajo correspondiente en el contenedor, establece el directorio de trabajo dentro del contenedor en consecuencia, y ejecuta nuestro script bash plantillado en el archivo `.command.sh`.

¡Todo el trabajo duro que tuvimos que hacer manualmente en la primera sección? ¡Nextflow lo hace por nosotros detrás de escenas!

```txt
 _______________________
< Hurray for robots...! >
 -----------------------
                                   ,-----.
                                   |     |
                                ,--|     |-.
                         __,----|  |     | |
                       ,;::     |  `_____' |
                       `._______|    i^i   |
                                `----| |---'| .
                           ,-------._| |== ||//
                           |       |_|P`.  /'/
                           `-------' 'Y Y/'/'
                                     .==\ /_\
   ^__^                             /   /'|  `i
   (oo)\_______                   /'   /  |   |
   (__)\       )\/\             /'    /   |   `i
       ||----w |           ___,;`----'.___L_,-'`\__
       ||     ||          i_____;----\.____i""\____\
```

### Conclusión

Sabe cómo usar contenedores en Nextflow para ejecutar procesos.

### ¿Qué sigue?

¡Tome un descanso!

Cuando esté listo, continúe con [**Parte 6: Hello Config**](./06_hello_config.md) para aprender cómo configurar la ejecución de su pipeline para adaptarse a su infraestructura así como gestionar la configuración de entradas y parámetros.

¡Es la última parte, y luego habrá terminado con este curso!

---

## Cuestionario

<quiz>
¿Qué es un contenedor?
- [ ] Un tipo de máquina virtual
- [ ] Un formato de compresión de archivos
- [x] Una unidad ejecutable ligera e independiente que incluye todo lo necesario para ejecutar una aplicación
- [ ] Un protocolo de red
</quiz>

<quiz>
¿Cuál es la diferencia entre una imagen de contenedor y una instancia de contenedor?
- [ ] Son lo mismo
- [x] Una imagen es una plantilla; una instancia es un contenedor en ejecución creado a partir de esa imagen
- [ ] Una instancia es una plantilla; una imagen es un contenedor en ejecución
- [ ] Las imágenes son para Docker; las instancias son para Singularity
</quiz>

<quiz>
¿Qué hace la bandera `-v` en un comando `docker run`?
- [ ] Habilita salida verbosa
- [ ] Valida el contenedor
- [x] Monta un volumen del sistema host en el contenedor
- [ ] Especifica la versión del contenedor

Aprenda más: [1.3.4. Montar datos en el contenedor](#134-mount-data-into-the-container)
</quiz>

<quiz>
¿Por qué necesita montar volúmenes cuando usa contenedores?
- [ ] Para mejorar el rendimiento del contenedor
- [ ] Para ahorrar espacio en disco
- [x] Porque los contenedores están aislados del sistema de archivos del host por defecto
- [ ] Para habilitar la red

Aprenda más: [1.3.4. Montar datos en el contenedor](#134-mount-data-into-the-container)
</quiz>

<quiz>
¿Cómo especifica un contenedor para un proceso de Nextflow?
- [ ] `docker 'container-uri'`
- [ ] `image 'container-uri'`
- [x] `container 'container-uri'`
- [ ] `use 'container-uri'`

Aprenda más: [2.3.1. Especificar un contenedor para cowpy](#231-specify-a-container-for-cowpy)
</quiz>

<quiz>
¿Qué configuración de `nextflow.config` habilita Docker para su workflow?
- [ ] `#!groovy process.docker = true`
- [x] `#!groovy docker.enabled = true`
- [ ] `#!groovy container.engine = 'docker'`
- [ ] `#!groovy docker.activate = true`

Aprenda más: [2.3.2. Habilitar el uso de Docker a través del archivo `nextflow.config`](#232-enable-use-of-docker-via-the-nextflowconfig-file)
</quiz>

<quiz>
¿Qué maneja automáticamente Nextflow cuando ejecuta un proceso en un contenedor? (Seleccione todos los que apliquen)
- [x] Descargar la imagen del contenedor si es necesario
- [x] Montar el directorio de trabajo
- [x] Ejecutar el script del proceso dentro del contenedor
- [x] Limpiar la instancia del contenedor después de la ejecución

Aprenda más: [2.3.4. Inspeccionar cómo Nextflow lanzó la tarea containerizada](#234-inspect-how-nextflow-launched-the-containerized-task)
</quiz>
