# Parte 3: Hello Workflow

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Traducción asistida por IA - [más información y sugerencias](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

<div class="video-wrapper">
  <iframe width="560" height="315" src="https://www.youtube.com/embed/_aO56V3iXGI?si=y8lAedhEHWaTV4zd&amp;list=PLPZ8WHdZGxmWKozQuzr27jyMGqp9kElVK&amp;cc_load_policy=1" title="YouTube video player" frameborder="0" allow="accelerometer; autoplay; clipboard-write; encrypted-media; gyroscope; picture-in-picture; web-share" referrerpolicy="strict-origin-when-cross-origin" allowfullscreen></iframe>
</div>

/// caption
:fontawesome-brands-youtube:{ .youtube } Consulte [la lista de reproducción completa](https://youtube.com/playlist?list=PLPZ8WHdZGxmWKozQuzr27jyMGqp9kElVK&si=eF7cLR62goy-lc6n) en el canal de YouTube de Nextflow.

:green_book: La transcripción del video está disponible [aquí](./transcripts/03_hello_workflow.md).
///

La mayoría de los workflows del mundo real involucran más de un paso.
En este módulo de capacitación, aprenderá a conectar procesos en un workflow de múltiples pasos.

Esto le enseñará la forma Nextflow de lograr lo siguiente:

1. Hacer que los datos fluyan de un proceso al siguiente
2. Recopilar salidas de múltiples llamadas a procesos en una sola llamada a proceso
3. Pasar parámetros adicionales a un proceso
4. Manejar múltiples salidas provenientes de un proceso

Para demostrarlo, continuaremos construyendo sobre el ejemplo Hello World independiente del dominio de las Partes 1 y 2.
Esta vez, vamos a realizar los siguientes cambios a nuestro workflow para reflejar mejor cómo las personas construyen workflows reales:

1. Agregar un segundo paso que convierta el saludo a mayúsculas.
2. Agregar un tercer paso que recopile todos los saludos transformados y los escriba en un solo archivo.
3. Agregar un parámetro para nombrar el archivo de salida final y pasarlo como entrada secundaria al paso de recopilación.
4. Hacer que el paso de recopilación también reporte una estadística simple sobre lo que se procesó.

??? info "Cómo comenzar desde esta sección"

    Esta sección del curso asume que ha completado las Partes 1-2 del curso [Hello Nextflow](./index.md), pero si se siente cómodo con los conceptos básicos cubiertos en esas secciones, puede comenzar desde aquí sin hacer nada especial.

---

## 0. Calentamiento: Ejecute `hello-workflow.nf`

Vamos a usar el script de workflow `hello-workflow.nf` como punto de partida.
Es equivalente al script producido al trabajar en la Parte 2 de este curso de capacitación, excepto que hemos eliminado las declaraciones `view()` y cambiado el destino de salida:

```groovy title="hello-workflow.nf" linenums="37" hl_lines="3"
output {
    first_output {
        path 'hello_workflow'
        mode 'copy'
    }
}
```

Este diagrama resume la operación actual del workflow.
Debería verse familiar, excepto que ahora estamos mostrando explícitamente que las salidas del proceso están empaquetadas en un canal, al igual que las entradas.
Vamos a darle un buen uso a ese canal de salida en un minuto.

<figure class="excalidraw">
--8<-- "docs/en/docs/hello_nextflow/img/hello-workflow-channels.svg"
</figure>

Solo para asegurarnos de que todo funciona, ejecute el script una vez antes de realizar cualquier cambio:

```bash
nextflow run hello-workflow.nf
```

??? success "Salida del comando"

    ```console
     N E X T F L O W   ~  version 25.10.2

    Launching `hello-workflow.nf` [admiring_lamarr] DSL2 - revision: 4d4053520d

    executor >  local (3)
    [b1/5826b5] process > sayHello (2) [100%] 3 of 3 ✔
    ```

Como anteriormente, encontrará los archivos de salida en la ubicación especificada en el bloque `output`.
Para este capítulo, está bajo `results/hello_workflow/`.

??? abstract "Contenido del directorio"

    ```console
    results/hello_workflow
    ├── Bonjour-output.txt
    ├── Hello-output.txt
    └── Holà-output.txt
    ```

Si eso funcionó para usted, está listo para aprender cómo ensamblar un workflow de múltiples pasos.

---

## 1. Agregue un segundo paso al workflow

Vamos a agregar un paso para convertir cada saludo a mayúsculas.

<figure class="excalidraw">
--8<-- "docs/en/docs/hello_nextflow/img/hello-multistep.svg"
</figure>

Para ello, necesitamos hacer tres cosas:

- Definir el comando que vamos a usar para hacer la conversión a mayúsculas.
- Escribir un nuevo proceso que envuelva el comando de conversión a mayúsculas.
- Llamar al nuevo proceso en el bloque workflow y configurarlo para que tome la salida del proceso `sayHello()` como entrada.

### 1.1. Defina el comando de conversión a mayúsculas y pruébelo en la terminal

Para hacer la conversión de los saludos a mayúsculas, vamos a usar una herramienta UNIX clásica llamada `tr` para 'reemplazo de texto', con la siguiente sintaxis:

```bash title="Sintaxis"
tr '[a-z]' '[A-Z]'
```

Este es un reemplazo de texto muy ingenuo que no tiene en cuenta las letras acentuadas, por lo que, por ejemplo, 'Holà' se convertirá en 'HOLà', pero hará un trabajo lo suficientemente bueno para demostrar los conceptos de Nextflow y eso es lo que importa.

Para probarlo, podemos ejecutar el comando `echo 'Hello World'` y canalizar su salida al comando `tr`:

```bash
echo 'Hello World' | tr '[a-z]' '[A-Z]' > UPPER-output.txt
```

La salida es un archivo de texto llamado `UPPER-output.txt` que contiene la versión en mayúsculas de la cadena `Hello World`.

??? abstract "Contenido del archivo"

    ```console title="UPPER-output.txt"
    HELLO WORLD
    ```

Eso es básicamente lo que vamos a intentar hacer con nuestro workflow.

### 1.2. Escriba el paso de conversión a mayúsculas como un proceso de Nextflow

Podemos modelar nuestro nuevo proceso en el primero, ya que queremos usar todos los mismos componentes.

Agregue la siguiente definición de proceso al script de workflow, justo debajo del primero:

```groovy title="hello-workflow.nf" linenums="20"
/*
 * Usa una herramienta de reemplazo de texto para convertir el saludo a mayúsculas
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

En este, componemos el segundo nombre de archivo de salida basado en el nombre del archivo de entrada, de manera similar a lo que hicimos originalmente para la salida del primer proceso.

### 1.3. Agregue una llamada al nuevo proceso en el bloque workflow

Ahora necesitamos decirle a Nextflow que realmente llame al proceso que acabamos de definir.

En el bloque workflow, realice el siguiente cambio de código:

=== "Después"

    ```groovy title="hello-workflow.nf" linenums="44" hl_lines="10-11"
    workflow {

        main:
        // crea un canal para entradas desde un archivo CSV
        greeting_ch = channel.fromPath(params.input)
                            .splitCsv()
                            .map { line -> line[0] }
        // emite un saludo
        sayHello(greeting_ch)
        // convierte el saludo a mayúsculas
        convertToUpper()

        publish:
        first_output = sayHello.out
    }
    ```

=== "Antes"

    ```groovy title="hello-workflow.nf" linenums="44"
    workflow {

        main:
        // crea un canal para entradas desde un archivo CSV
        greeting_ch = channel.fromPath(params.input)
                            .splitCsv()
                            .map { line -> line[0] }
        // emite un saludo
        sayHello(greeting_ch)

        publish:
        first_output = sayHello.out
    }
    ```

Esto aún no es funcional porque no hemos especificado qué debe ser la entrada del proceso `convertToUpper()`.

### 1.4. Pase la salida del primer proceso al segundo proceso

Ahora necesitamos hacer que la salida del proceso `sayHello()` fluya hacia el proceso `convertToUpper()`.

Convenientemente, Nextflow empaqueta automáticamente la salida de un proceso en un canal, como se muestra en el diagrama en la sección de calentamiento.
Podemos referirnos al canal de salida de un proceso como `<proceso>.out`.

Entonces, la salida del proceso `sayHello` es un canal llamado `sayHello.out`, que podemos conectar directamente a la llamada a `convertToUpper()`.

<figure class="excalidraw">
--8<-- "docs/en/docs/hello_nextflow/img/hello-multistep-connector.svg"
</figure>

En el bloque workflow, realice el siguiente cambio de código:

=== "Después"

    ```groovy title="hello-workflow.nf" linenums="53" hl_lines="2"
        // convierte el saludo a mayúsculas
        convertToUpper(sayHello.out)
    ```

=== "Antes"

    ```groovy title="hello-workflow.nf" linenums="53" hl_lines="2"
        // convierte el saludo a mayúsculas
        convertToUpper()
    ```

Para un caso simple como este (una salida a una entrada), ¡eso es todo lo que necesitamos hacer para conectar dos procesos!

### 1.5. Configure la publicación de salida del workflow

Finalmente, actualicemos las salidas del workflow para publicar también los resultados del segundo proceso.

#### 1.5.1. Actualice la sección `publish:` del bloque `workflow`

En el bloque `workflow`, realice el siguiente cambio de código:

=== "Después"

    ```groovy title="hello-workflow.nf" linenums="56" hl_lines="3"
        publish:
        first_output = sayHello.out
        uppercased = convertToUpper.out
    }
    ```

=== "Antes"

    ```groovy title="hello-workflow.nf" linenums="56"
        publish:
        first_output = sayHello.out
    }
    ```

La lógica es la misma que anteriormente.

#### 1.5.2. Actualice el bloque `output`

En el bloque `output`, realice el siguiente cambio de código:

=== "Después"

    ```groovy title="hello-workflow.nf" linenums="61" hl_lines="6-9"
    output {
        first_output {
            path 'hello_workflow'
            mode 'copy'
        }
        uppercased {
            path 'hello_workflow'
            mode 'copy'
        }
    }
    ```

=== "Antes"

    ```groovy title="hello-workflow.nf" linenums="61"
    output {
        first_output {
            path 'hello_workflow'
            mode 'copy'
        }
    }
    ```

Una vez más, la lógica es la misma que antes.

Esto le muestra que puede controlar la configuración de salida a un nivel muy granular, para cada salida individual.
Siéntase libre de intentar cambiar las rutas o el modo de publicación para uno de los procesos para ver qué sucede.

Por supuesto, eso significa que estamos repitiendo cierta información aquí, lo que podría volverse inconveniente si quisiéramos actualizar la ubicación para todas las salidas de la misma manera.
Más adelante en el curso, aprenderá cómo configurar estos ajustes para múltiples salidas de manera estructurada.

### 1.6. Ejecute el workflow con `-resume`

Probemos esto usando la bandera `-resume`, ya que ya hemos ejecutado el primer paso del workflow con éxito.

```bash
nextflow run hello-workflow.nf -resume
```

??? success "Salida del comando"

    ```console
     N E X T F L O W   ~  version 25.10.2

    Launching `hello-workflow.nf` [high_cantor] DSL2 - revision: d746983511

    executor >  local (3)
    [ab/816321] process > sayHello (3)       [100%] 3 of 3, cached: 3 ✔
    [e0/ecf81b] process > convertToUpper (3) [100%] 3 of 3 ✔
    ```

Ahora hay una línea adicional en la salida de la consola que corresponde al nuevo proceso que acabamos de agregar.

Encontrará las salidas en el directorio `results/hello_workflow` como se estableció en el bloque `output`.

??? abstract "Contenido del directorio"

    ```console
    results/hello_workflow/
    ├── Bonjour-output.txt
    ├── Hello-output.txt
    ├── Holà-output.txt
    ├── UPPER-Bonjour-output.txt
    ├── UPPER-Hello-output.txt
    └── UPPER-Holà-output.txt
    ```

¡Eso es conveniente! Pero aún vale la pena echar un vistazo dentro del directorio de trabajo de una de las llamadas al segundo proceso.

??? abstract "Contenido del directorio"

    ```console
    work/e0/ecf81b4cacc648b9b994218d5b29d7/
    ├── Holà-output.txt -> /workspaces/training/hello-nextflow/work/ab/81632178cd37e9e815959278808819/Holà-output.txt
    └── UPPER-Holà-output.txt
    ```

Note que hay dos archivos `*-output`: la salida del primer proceso así como la salida del segundo.

La salida del primer proceso está ahí porque Nextflow la **preparó** allí para tener todo lo necesario para la ejecución dentro del mismo subdirectorio.

Sin embargo, en realidad es un enlace simbólico que apunta al archivo original en el subdirectorio de la primera llamada al proceso.
Por defecto, cuando se ejecuta en una sola máquina como lo estamos haciendo aquí, Nextflow usa enlaces simbólicos en lugar de copias para preparar archivos de entrada e intermedios.

Ahora, antes de continuar, piense en cómo todo lo que hicimos fue conectar la salida de `sayHello` a la entrada de `convertToUpper` y los dos procesos pudieron ejecutarse en serie.
Nextflow hizo el trabajo pesado de manejar archivos de entrada y salida individuales y pasarlos entre los dos comandos por nosotros.

Esta es una de las razones por las que los canales de Nextflow son tan poderosos: se encargan del trabajo tedioso involucrado en conectar los pasos del workflow.

### Conclusión

Sabe cómo encadenar procesos conectando la salida de un paso como entrada al siguiente paso.

### ¿Qué sigue?

Aprenda cómo recopilar salidas de llamadas a procesos por lotes y alimentarlas en un solo proceso.

---

## 2. Agregue un tercer paso para recopilar todos los saludos

Cuando usamos un proceso para aplicar una transformación a cada uno de los elementos en un canal, como lo estamos haciendo aquí con los múltiples saludos, a veces queremos recopilar elementos del canal de salida de ese proceso y alimentarlos en otro proceso que realiza algún tipo de análisis o suma.

Para demostrarlo, agregaremos un nuevo paso a nuestro pipeline que recopila todos los saludos en mayúsculas producidos por el proceso `convertToUpper` y los escribe en un solo archivo.

<figure class="excalidraw">
--8<-- "docs/en/docs/hello_nextflow/img/hello-collect.svg"
</figure>

Sin arruinar la sorpresa, esto va a involucrar un operador muy útil.

### 2.1. Defina el comando de recopilación y pruébelo en la terminal

El paso de recopilación que queremos agregar a nuestro workflow usará el comando `cat` para concatenar múltiples saludos en mayúsculas en un solo archivo.

Ejecutemos el comando por sí solo en la terminal para verificar que funciona como se espera, tal como lo hemos hecho anteriormente.

Ejecute lo siguiente en su terminal:

```bash
echo 'Hello' | tr '[a-z]' '[A-Z]' > UPPER-Hello-output.txt
echo 'Bonjour' | tr '[a-z]' '[A-Z]' > UPPER-Bonjour-output.txt
echo 'Holà' | tr '[a-z]' '[A-Z]' > UPPER-Holà-output.txt
cat UPPER-Hello-output.txt UPPER-Bonjour-output.txt UPPER-Holà-output.txt > COLLECTED-output.txt
```

La salida es un archivo de texto llamado `COLLECTED-output.txt` que contiene las versiones en mayúsculas de los saludos originales.

??? abstract "Contenido del archivo"

    ```console title="COLLECTED-output.txt"
    HELLO
    BONJOUR
    HOLà
    ```

Ese es el resultado que queremos lograr con nuestro workflow.

### 2.2. Cree un nuevo proceso para hacer el paso de recopilación

Creemos un nuevo proceso y llamémoslo `collectGreetings()`.
Podemos comenzar a escribirlo basándonos en lo que hemos visto antes.

#### 2.2.1. Escriba las partes 'obvias' del proceso

Agregue la siguiente definición de proceso al script de workflow:

```groovy title="hello-workflow.nf" linenums="37"
/*
 * Recopila saludos en mayúsculas en un solo archivo de salida
 */
process collectGreetings {

    input:
    ???

    output:
    path "COLLECTED-output.txt"

    script:
    """
    cat ??? > 'COLLECTED-output.txt'
    """
}
```

Esto es lo que podemos escribir con confianza basándonos en lo que ha aprendido hasta ahora.
¡Pero esto no es funcional!
Omite la(s) definición(es) de entrada y la primera mitad del comando del script porque necesitamos averiguar cómo escribir eso.

#### 2.2.2. Defina las entradas a `collectGreetings()`

Necesitamos recopilar los saludos de todas las llamadas al proceso `convertToUpper()`.
¿Qué sabemos que podemos obtener del paso anterior en el workflow?

El canal de salida de `convertToUpper()` contendrá las rutas a los archivos individuales que contienen los saludos en mayúsculas.
Eso equivale a un espacio de entrada; llamémoslo `input_files` por simplicidad.

En el bloque de proceso, realice el siguiente cambio de código:

=== "Después"

    ```groovy title="hello-workflow.nf" linenums="42" hl_lines="2"
          input:
          path input_files
    ```

=== "Antes"

    ```groovy title="hello-workflow.nf" linenums="42" hl_lines="2"
          input:
          ???
    ```

Note que usamos el prefijo `path` aunque esperamos que esto contenga múltiples archivos.

#### 2.2.3. Componga el comando de concatenación

Aquí es donde las cosas podrían complicarse un poco, porque necesitamos poder manejar un número arbitrario de archivos de entrada.
Específicamente, no podemos escribir el comando por adelantado, por lo que necesitamos decirle a Nextflow cómo componerlo en tiempo de ejecución basándose en qué entradas fluyen hacia el proceso.

En otras palabras, si tenemos un canal de entrada que contiene el elemento `[file1.txt, file2.txt, file3.txt]`, necesitamos que Nextflow lo convierta en `cat file1.txt file2.txt file3.txt`.

Afortunadamente, Nextflow está muy dispuesto a hacer eso por nosotros si simplemente escribimos `cat ${input_files}` en el comando del script.

En el bloque de proceso, realice el siguiente cambio de código:

=== "Después"

    ```groovy title="hello-workflow.nf" linenums="54" hl_lines="3"
        script:
        """
        cat ${input_files} > 'COLLECTED-output.txt'
        """
    ```

=== "Antes"

    ```groovy title="hello-workflow.nf" linenums="54"
        script:
        """
        cat ??? > 'COLLECTED-output.txt'
        """
    ```

En teoría, esto debería manejar cualquier número arbitrario de archivos de entrada.

!!! tip

    Algunas herramientas de línea de comandos requieren proporcionar un argumento (como `-input`) para cada archivo de entrada.
    En ese caso, tendríamos que hacer un poco de trabajo extra para componer el comando.
    Puede ver un ejemplo de esto en el curso de capacitación [Nextflow for Genomics](../../nf4_science/genomics/).

### 2.3. Agregue el paso de recopilación al workflow

Ahora solo deberíamos necesitar llamar al proceso de recopilación en la salida del paso de conversión a mayúsculas.
Ese también es un canal, llamado `convertToUpper.out`.

<figure class="excalidraw">
--8<-- "docs/en/docs/hello_nextflow/img/hello-collect-connector.svg"
</figure>

#### 2.3.1. Conecte las llamadas a procesos

En el bloque workflow, realice el siguiente cambio de código:

=== "Después"

    ```groovy title="hello-workflow.nf" linenums="75" hl_lines="4 5"
        // convierte el saludo a mayúsculas
        convertToUpper(sayHello.out)

        // recopila todos los saludos en un archivo
        collectGreetings(convertToUpper.out)
    }
    ```

=== "Antes"

    ```groovy title="hello-workflow.nf" linenums="75"
        // convierte el saludo a mayúsculas
        convertToUpper(sayHello.out)
    }
    ```

Esto conecta la salida de `convertToUpper()` a la entrada de `collectGreetings()`.

#### 2.3.2. Ejecute el workflow con `-resume`

Intentémoslo.

```bash
nextflow run hello-workflow.nf -resume
```

??? success "Salida del comando"

    ```console hl_lines="8"
    N E X T F L O W   ~  version 25.10.2

    Launching `hello-workflow.nf` [mad_gilbert] DSL2 - revision: 6acfd5e28d

    executor >  local (3)
    [79/33b2f0] sayHello (2)         | 3 of 3, cached: 3 ✔
    [99/79394f] convertToUpper (3)   | 3 of 3, cached: 3 ✔
    [47/50fe4a] collectGreetings (1) | 3 of 3 ✔
    ```

Se ejecuta con éxito, incluyendo el tercer paso.

Sin embargo, observe el número de llamadas para `collectGreetings()` en la última línea.
Solo esperábamos una, pero hay tres.

Ahora eche un vistazo al contenido del archivo de salida final.

??? abstract "Contenido del archivo"

    ```console title="results/COLLECTED-output.txt"
    Holà
    ```

Oh no. El paso de recopilación se ejecutó individualmente en cada saludo, lo cual NO es lo que queríamos.

<figure class="excalidraw">
--8<-- "docs/en/docs/hello_nextflow/img/hello-collect-no-operator.svg"
</figure>

Necesitamos hacer algo para decirle a Nextflow explícitamente que queremos que ese tercer paso se ejecute en todos los elementos del canal de salida de `convertToUpper()`.

### 2.4. Use un operador para recopilar los saludos en una sola entrada

Sí, una vez más la respuesta a nuestro problema es un operador.

Específicamente, vamos a usar el operador [`collect()`](https://nextflow.io/docs/latest/reference/operator.html#collect), con un nombre muy apropiado.

#### 2.4.1. Agregue el operador `collect()`

Esta vez se verá un poco diferente porque no estamos agregando el operador en el contexto de una fábrica de canales; lo estamos agregando a un canal de salida.

Tomamos el `convertToUpper.out` y agregamos el operador `collect()`, lo que nos da `convertToUpper.out.collect()`.
Podemos conectar eso directamente en la llamada al proceso `collectGreetings()`.

En el bloque workflow, realice el siguiente cambio de código:

=== "Después"

    ```groovy title="hello-workflow.nf" linenums="73" hl_lines="2"
        // recopila todos los saludos en un archivo
        collectGreetings(convertToUpper.out.collect())
    }
    ```

=== "Antes"

    ```groovy title="hello-workflow.nf" linenums="73" hl_lines="2"
        // recopila todos los saludos en un archivo
        collectGreetings(convertToUpper.out)
    }
    ```

#### 2.4.2. Agregue algunas declaraciones `view()`

También incluyamos un par de declaraciones `view()` para visualizar los estados antes y después del contenido del canal.

=== "Después"

    ```groovy title="hello-workflow.nf" linenums="73" hl_lines="4-6"
        // recopila todos los saludos en un archivo
        collectGreetings(convertToUpper.out.collect())

        // declaraciones view opcionales
        convertToUpper.out.view { contents -> "Antes de collect: $contents" }
        convertToUpper.out.collect().view { contents -> "Después de collect: $contents" }
    }
    ```

=== "Antes"

    ```groovy title="hello-workflow.nf" linenums="73"
        // recopila todos los saludos en un archivo
        collectGreetings(convertToUpper.out.collect())
    }
    ```

Las declaraciones `view()` pueden ir donde quiera; las pusimos justo después de la llamada para mayor legibilidad.

#### 2.4.3. Ejecute el workflow nuevamente con `-resume`

Intentémoslo:

```bash
nextflow run hello-workflow.nf -resume
```

??? success "Salida del comando"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `hello-workflow.nf` [soggy_franklin] DSL2 - revision: bc8e1b2726

    [d6/cdf466] sayHello (1)       | 3 of 3, cached: 3 ✔
    [99/79394f] convertToUpper (2) | 3 of 3, cached: 3 ✔
    [1e/83586c] collectGreetings   | 1 of 1 ✔
    Antes de collect: /workspaces/training/hello-nextflow/work/b3/d52708edba8b864024589285cb3445/UPPER-Bonjour-output.txt
    Antes de collect: /workspaces/training/hello-nextflow/work/99/79394f549e3040dfc2440f69ede1fc/UPPER-Hello-output.txt
    Antes de collect: /workspaces/training/hello-nextflow/work/aa/56bfe7cf00239dc5badc1d04b60ac4/UPPER-Holà-output.txt
    Después de collect: [/workspaces/training/hello-nextflow/work/b3/d52708edba8b864024589285cb3445/UPPER-Bonjour-output.txt, /workspaces/training/hello-nextflow/work/99/79394f549e3040dfc2440f69ede1fc/UPPER-Hello-output.txt, /workspaces/training/hello-nextflow/work/aa/56bfe7cf00239dc5badc1d04b60ac4/UPPER-Holà-output.txt]
    ```

Se ejecuta con éxito, aunque la salida del registro puede verse un poco más desordenada que esto (la limpiamos para mayor legibilidad).

¡Esta vez el tercer paso solo se llamó una vez!
Mirando la salida de las declaraciones `view()`, vemos lo siguiente:

- Tres declaraciones `Antes de collect:`, una para cada saludo: en ese punto las rutas de archivo son elementos individuales en el canal.
- Una sola declaración `Después de collect:`: las tres rutas de archivo ahora están empaquetadas en un solo elemento.

Podemos resumir eso con el siguiente diagrama:

<figure class="excalidraw">
--8<-- "docs/en/docs/hello_nextflow/img/hello-collect-WITH-operator.svg"
</figure>

Finalmente, puede echar un vistazo al contenido del archivo de salida para satisfacerse de que todo funcionó correctamente.

??? abstract "Contenido del archivo"

    ```console title="results/COLLECTED-output.txt"
    BONJOUR
    HELLO
    HOLà
    ```

Esta vez tenemos los tres saludos en el archivo de salida final. ¡Éxito!

!!! note

    Si ejecuta esto varias veces sin `-resume`, verá que el orden de los saludos cambia de una ejecución a otra.
    Esto le muestra que no se garantiza que el orden en que los elementos fluyen a través de las llamadas a procesos sea consistente.

#### 2.4.4. Elimine las declaraciones `view()` para mayor legibilidad

Antes de pasar a la siguiente sección, le recomendamos que elimine las declaraciones `view()` para evitar saturar la salida de la consola.

=== "Después"

    ```groovy title="hello-workflow.nf" linenums="73"
        // recopila todos los saludos en un archivo
        collectGreetings(convertToUpper.out.collect())
    ```

=== "Antes"

    ```groovy title="hello-workflow.nf" linenums="73" hl_lines="4-6"
        // recopila todos los saludos en un archivo
        collectGreetings(convertToUpper.out.collect())

        // declaraciones view opcionales
        convertToUpper.out.view { contents -> "Antes de collect: $contents" }
        convertToUpper.out.collect().view { contents -> "Después de collect: $contents" }
    ```

Esta es básicamente la operación inversa del punto 2.4.2.

### Conclusión

Sabe cómo recopilar salidas de un lote de llamadas a procesos y alimentarlas en un paso de análisis o suma conjunta.

Para recapitular, esto es lo que ha construido hasta ahora:

<figure class="excalidraw">
--8<-- "docs/en/docs/hello_nextflow/img/hello-collect.svg"
</figure>

### ¿Qué sigue?

Aprenda cómo pasar más de una entrada a un proceso.

---

## 3. Pase parámetros adicionales a un proceso

Queremos poder nombrar el archivo de salida final algo específico para procesar lotes posteriores de saludos sin sobrescribir los resultados finales.

Para ello, vamos a hacer los siguientes refinamientos al workflow:

- Modificar el proceso recopilador para aceptar un nombre definido por el usuario para el archivo de salida (`batch_name`)
- Agregar un parámetro de línea de comandos al workflow (`--batch`) y pasarlo al proceso recopilador

<figure class="excalidraw">
--8<-- "docs/en/docs/hello_nextflow/img/hello-collect-batch.svg"
</figure>

### 3.1. Modifique el proceso recopilador

Vamos a necesitar declarar la entrada adicional e integrarla en el nombre del archivo de salida.

#### 3.1.1. Declare la entrada adicional

Buenas noticias: podemos declarar tantas variables de entrada como queramos en la definición del proceso.
Llamemos a esta `batch_name`.

En el bloque de proceso, realice el siguiente cambio de código:

=== "Después"

    ```groovy title="hello-workflow.nf" linenums="42" hl_lines="3"
        input:
        path input_files
        val batch_name
    ```

=== "Antes"

    ```groovy title="hello-workflow.nf" linenums="42"
        input:
        path input_files
    ```

Puede configurar sus procesos para esperar tantas entradas como desee.
En este momento, todas están configuradas para ser entradas requeridas; _debe_ proporcionar un valor para que el workflow funcione.

Aprenderá cómo administrar entradas requeridas vs. opcionales más adelante en su viaje con Nextflow.

#### 3.1.2. Use la variable `batch_name` en el nombre del archivo de salida

Podemos insertar la variable en el nombre del archivo de salida de la misma manera que hemos compuesto nombres de archivo dinámicos antes.

En el bloque de proceso, realice el siguiente cambio de código:

=== "Después"

    ```groovy title="hello-workflow.nf" linenums="46" hl_lines="2 6"
        output:
        path "COLLECTED-${batch_name}-output.txt"

        script:
        """
        cat ${input_files} > 'COLLECTED-${batch_name}-output.txt'
        """
    ```

=== "Antes"

    ```groovy title="hello-workflow.nf" linenums="46" hl_lines="2 6"
        output:
        path "COLLECTED-output.txt"

        script:
        """
        cat ${input_files} > 'COLLECTED-output.txt'
        """
    ```

Esto configura el proceso para usar el valor de `batch_name` para generar un nombre de archivo específico para la salida final del workflow.

### 3.2. Agregue un parámetro de línea de comandos `batch`

Ahora necesitamos una forma de proporcionar el valor para `batch_name` y alimentarlo a la llamada del proceso.

#### 3.2.1. Use `params` para configurar el parámetro

Ya sabe cómo usar el sistema `params` para declarar parámetros CLI.
Usemos eso para declarar un parámetro `batch` (con un valor predeterminado porque somos perezosos).

En la sección de parámetros del pipeline, realice los siguientes cambios de código:

=== "Después"

    ```groovy title="hello-workflow.nf" linenums="55" hl_lines="6"
    /*
     * Parámetros del pipeline
     */
    params {
        input: Path = 'data/greetings.csv'
        batch: String = 'batch'
    }
    ```

=== "Antes"

    ```groovy title="hello-workflow.nf" linenums="55"
    /*
     * Parámetros del pipeline
     */
    params {
        input: Path = 'data/greetings.csv'
    }
    ```

Al igual que demostramos para `--input`, puede anular ese valor predeterminado especificando un valor con `--batch` en la línea de comandos.

#### 3.2.2. Pase el parámetro `batch` al proceso

Para proporcionar el valor del parámetro al proceso, necesitamos agregarlo en la llamada al proceso.

En el bloque workflow, realice el siguiente cambio de código:

=== "Después"

    ```groovy title="hello-workflow.nf" linenums="74" hl_lines="2"
        // recopila todos los saludos en un archivo
        collectGreetings(convertToUpper.out.collect(), params.batch)
    ```

=== "Antes"

    ```groovy title="hello-workflow.nf" linenums="74" hl_lines="2"
        // recopila todos los saludos en un archivo
        collectGreetings(convertToUpper.out.collect())
    ```

Ve que para proporcionar múltiples entradas a un proceso, simplemente las enumera en los paréntesis de la llamada, separadas por comas.

!!! warning

    DEBE proporcionar las entradas al proceso en el MISMO ORDEN EXACTO en que están listadas en el bloque de definición de entrada del proceso.

### 3.3. Ejecute el workflow

Intentemos ejecutar esto con un nombre de lote en la línea de comandos.

```bash
nextflow run hello-workflow.nf -resume --batch trio
```

??? success "Salida del comando"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `hello-workflow.nf` [confident_rutherford] DSL2 - revision: bc58af409c

    executor >  local (1)
    [79/33b2f0] sayHello (2)       | 3 of 3, cached: 3 ✔
    [99/79394f] convertToUpper (2) | 3 of 3, cached: 3 ✔
    [b5/f19efe] collectGreetings   | 1 of 1 ✔
    ```

Se ejecuta con éxito y produce la salida deseada:

??? abstract "Contenido del archivo"

    ```console title="results/COLLECTED-trio-output.txt"
    HELLO
    BONJOUR
    HOLà
    ```

Ahora, siempre que especifiquemos el parámetro apropiadamente, las ejecuciones posteriores en otros lotes de entradas no sobrescribirán los resultados anteriores.

### Conclusión

Sabe cómo pasar más de una entrada a un proceso.

### ¿Qué sigue?

Aprenda cómo emitir múltiples salidas y manejarlas convenientemente.

---

## 4. Agregue una salida al paso recopilador

Hasta ahora hemos estado usando procesos que solo producían una salida cada uno.
Pudimos acceder a sus respectivas salidas muy convenientemente usando la sintaxis `<proceso>.out`, que usamos tanto en el contexto de pasar una salida al siguiente proceso (por ejemplo, `convertToUpper(sayHello.out)`) como en el contexto de la sección `publish:` (por ejemplo, `first_output = sayHello.out`).

¿Qué sucede cuando un proceso produce más de una?
¿Cómo manejamos las múltiples salidas?
¿Podemos seleccionar y usar una salida específica?

Todas son excelentes preguntas, y la respuesta corta es ¡sí podemos!

Las múltiples salidas se empaquetarán en canales separados.
Podemos elegir darles nombres a esos canales de salida, lo que facilita referirse a ellos individualmente más adelante, o podemos referirnos a ellos por índice.

Para fines de demostración, digamos que queremos contar el número de saludos que se están recopilando para un lote dado de entradas e informarlo en un archivo.

### 4.1. Modifique el proceso para contar y emitir el número de saludos

Esto requerirá dos cambios clave en la definición del proceso: necesitamos una forma de contar los saludos y escribir un archivo de informe, luego necesitamos agregar ese archivo de informe al bloque `output` del proceso.

#### 4.1.1. Cuente el número de saludos recopilados

Convenientemente, Nextflow nos permite agregar código arbitrario en el bloque `script:` de la definición del proceso, lo que resulta muy útil para hacer cosas como esta.

Eso significa que podemos usar la función incorporada `size()` de Nextflow para obtener el número de archivos en el array `input_files`, y escribir el resultado en un archivo con un comando `echo`.

En el bloque de proceso `collectGreetings`, realice los siguientes cambios de código:

=== "Después"

    ```groovy title="hello-workflow.nf" linenums="55" hl_lines="2 5"
        script:
        count_greetings = input_files.size()
        """
        cat ${input_files} > 'COLLECTED-${batch_name}-output.txt'
        echo 'Hubo ${count_greetings} saludos en este lote.' > '${batch_name}-report.txt'
        """
    ```

=== "Antes"

    ```groovy title="hello-workflow.nf" linenums="55"
        script:
        """
        cat ${input_files} > 'COLLECTED-${batch_name}-output.txt'
        """
    ```

La variable `count_greetings` se calculará en tiempo de ejecución.

#### 4.1.2. Emita el archivo de informe y nombre las salidas

En principio, todo lo que necesitamos hacer es agregar el archivo de informe al bloque `output:`.

Sin embargo, mientras lo hacemos, también vamos a agregar algunas etiquetas `emit:` a nuestras declaraciones de salida. Estas nos permitirán seleccionar las salidas por nombre en lugar de tener que usar índices posicionales.

En el bloque de proceso, realice el siguiente cambio de código:

=== "Después"

    ```groovy title="hello-workflow.nf" linenums="46" hl_lines="2 3"
        output:
        path "COLLECTED-${batch_name}-output.txt", emit: outfile
        path "${batch_name}-report.txt", emit: report
    ```

=== "Antes"

    ```groovy title="hello-workflow.nf" linenums="46"
        output:
        path "COLLECTED-${batch_name}-output.txt"
    ```

Las etiquetas `emit:` son opcionales, y podríamos haber agregado una etiqueta a solo una de las salidas.
Pero como dice el refrán, ¿por qué no ambas?

!!! tip

    Si no nombra las salidas de un proceso usando `emit:`, aún puede acceder a ellas individualmente usando su respectivo índice (basado en cero).
    Por ejemplo, usaría `<proceso>.out[0]` para obtener la primera salida, `<proceso>.out[1]` para obtener la segunda salida, y así sucesivamente.

    Preferimos nombrar las salidas porque de lo contrario, es demasiado fácil tomar el índice incorrecto por error, especialmente cuando el proceso produce muchas salidas.

### 4.2. Actualice las salidas del workflow

Ahora que tenemos dos salidas provenientes del proceso `collectGreetings`, la salida `collectGreetings.out` contiene dos canales:

- `collectGreetings.out.outfile` contiene el archivo de salida final
- `collectGreetings.out.report` contiene el archivo de informe

<figure class="excalidraw">
--8<-- "docs/en/docs/hello_nextflow/img/hello-collect-report.svg"
</figure>

Necesitamos actualizar las salidas del workflow en consecuencia.

#### 4.2.1. Actualice la sección `publish:`

En el bloque `workflow`, realice el siguiente cambio de código:

=== "Después"

    ```groovy title="hello-workflow.nf" linenums="80" hl_lines="4 5"
        publish:
        first_output = sayHello.out
        uppercased = convertToUpper.out
        collected = collectGreetings.out.outfile
        batch_report = collectGreetings.out.report
    ```

=== "Antes"

    ```groovy title="hello-workflow.nf" linenums="80" hl_lines="4"
        publish:
        first_output = sayHello.out
        uppercased = convertToUpper.out
        collected = collectGreetings.out
    ```

Como puede ver, referirse a salidas de procesos específicas ahora es trivial.
Cuando vayamos a agregar un paso más a nuestro pipeline en la Parte 5 (Contenedores), podremos referirnos fácilmente a `collectGreetings.out.outfile` y pasarlo al nuevo proceso (spoiler: el nuevo proceso se llama `cowpy`).

Pero por ahora, terminemos de actualizar las salidas a nivel de workflow.

#### 4.2.2. Actualice el bloque `output`

En el bloque `output`, realice el siguiente cambio de código:

=== "Después"

    ```groovy title="hello-workflow.nf" linenums="86" hl_lines="14-17"
    output {
        first_output {
            path 'hello_workflow'
            mode 'copy'
        }
        uppercased {
            path 'hello_workflow'
            mode 'copy'
        }
        collected {
            path 'hello_workflow'
            mode 'copy'
        }
        batch_report {
            path 'hello_workflow'
            mode 'copy'
        }
    }
    ```

=== "Antes"

    ```groovy title="hello-workflow.nf" linenums="80"
    output {
        first_output {
            path 'hello_workflow'
            mode 'copy'
        }
        uppercased {
            path 'hello_workflow'
            mode 'copy'
        }
        collected {
            path 'hello_workflow'
            mode 'copy'
        }
    }
    ```

No necesitamos actualizar la definición de salida `collected` ya que ese nombre no ha cambiado.
Solo necesitamos agregar la nueva salida.

### 4.3. Ejecute el workflow

Intentemos ejecutar esto con el lote actual de saludos.

```bash
nextflow run hello-workflow.nf -resume --batch trio
```

??? success "Salida del comando"

    ```console
     N E X T F L O W   ~  version 25.10.2

    Launching `hello-workflow.nf` [ecstatic_wilson] DSL2 - revision: c80285f8c8

    executor >  local (1)
    [c5/4c6ca9] sayHello (3)       [100%] 3 of 3, cached: 3 ✔
    [0e/6cbc59] convertToUpper (3) [100%] 3 of 3, cached: 3 ✔
    [02/61ead2] collectGreetings   [100%] 1 of 1 ✔
    ```

Si mira en el directorio `results/hello_workflow/`, encontrará el nuevo archivo de informe, `trio-report.txt`.
Ábralo para verificar que el workflow informó correctamente el recuento de saludos que se procesaron.

??? abstract "Contenido del archivo"

    ```txt title="trio-report.txt"
    Hubo 3 saludos en este lote.
    ```

<figure class="excalidraw">
--8<-- "docs/en/docs/hello_nextflow/img/hello-collect-4-way.svg"
</figure>

Siéntase libre de agregar más saludos al CSV y probar qué sucede.

### Conclusión

Sabe cómo hacer que un proceso emita múltiples salidas nombradas y cómo manejarlas apropiadamente a nivel de workflow.

Más generalmente, comprende los principios clave involucrados en conectar procesos de maneras comunes.

### ¿Qué sigue?

Tómese un descanso extra largo, se lo ha ganado.

Cuando esté listo, continúe con [**Parte 4: Hello Modules**](./04_hello_modules.md) para aprender cómo modularizar su código para mejor mantenibilidad y eficiencia del código.

---

## Cuestionario

<quiz>
¿Cómo accede a la salida de un proceso en el bloque workflow?
- [ ] `process.output`
- [ ] `output.processName`
- [x] `processName.out`
- [ ] `get(processName)`

Más información: [1.4. Pase la salida del primer proceso al segundo proceso](#14-pase-la-salida-del-primer-proceso-al-segundo-proceso)
</quiz>

<quiz>
¿Qué determina el orden de ejecución de procesos en Nextflow?
- [ ] El orden en que los procesos están escritos en el bloque workflow
- [ ] Orden alfabético por nombre de proceso
- [x] Dependencias de datos entre procesos
- [ ] Orden aleatorio para ejecución paralela

Más información: [1.4. Pase la salida del primer proceso al segundo proceso](#14-pase-la-salida-del-primer-proceso-al-segundo-proceso)
</quiz>

<quiz>
¿Qué operador debe reemplazar `???` para reunir todas las salidas en una sola lista para el proceso descendente?

```groovy hl_lines="4"
workflow {
    greetings_ch = Channel.of('Hello', 'Bonjour', 'Hola')
    SAYHELLO(greetings_ch)
    GATHER_ALL(SAYHELLO.out.???)
}
```

- [ ] `flatten()`
- [x] `collect()`
- [ ] `mix()`
- [ ] `join()`

Más información: [2.4. Use un operador para recopilar los saludos en una sola entrada](#24-use-un-operador-para-recopilar-los-saludos-en-una-sola-entrada)
</quiz>

<quiz>
¿Cuándo debe usar el operador `collect()`?
- [ ] Cuando quiere procesar elementos en paralelo
- [ ] Cuando necesita filtrar el contenido del canal
- [x] Cuando un proceso descendente necesita todos los elementos de un proceso ascendente
- [ ] Cuando quiere dividir datos entre múltiples procesos

Más información: [2.4. Use un operador para recopilar los saludos en una sola entrada](#24-use-un-operador-para-recopilar-los-saludos-en-una-sola-entrada)
</quiz>

<quiz>
¿Cómo accede a una salida nombrada de un proceso?
- [ ] `processName.outputName`
- [ ] `processName.get(outputName)`
- [x] `processName.out.outputName`
- [ ] `output.processName.outputName`

Más información: [4.1.2. Emita el archivo de informe y nombre las salidas](#412-emita-el-archivo-de-informe-y-nombre-las-salidas)
</quiz>

<quiz>
¿Cuál es la sintaxis correcta para nombrar una salida en un proceso?
- [ ] `name: outputName`
- [ ] `output: outputName`
- [x] `emit: outputName`
- [ ] `label: outputName`

Más información: [4.1.2. Emita el archivo de informe y nombre las salidas](#412-emita-el-archivo-de-informe-y-nombre-las-salidas)
</quiz>

<quiz>
Al proporcionar múltiples entradas a un proceso, ¿qué debe ser cierto?
- [ ] Todas las entradas deben ser del mismo tipo
- [ ] Las entradas deben proporcionarse en orden alfabético
- [x] El orden de las entradas debe coincidir con el orden definido en el bloque de entrada
- [ ] Solo se pueden proporcionar dos entradas a la vez

Más información: [3. Pase más de una entrada a un proceso](#3-pase-parametros-adicionales-a-un-proceso)
</quiz>
