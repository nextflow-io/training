---
descripción: Primeros pasos con Nextflow
---

# Introducción

## Conceptos básicos

Nextflow es un orquestador de flujos de trabajo y un lenguaje de dominio específico (DSL) que facilita la escritura de flujos computacionales con uso intensivo de datos.

Está diseñado en torno a la idea de que la plataforma Linux es la _lingua franca_ de la ciencia de datos. Linux proporciona muchas herramientas de linea de comandos simples pero poderosas que, cuando se encadenan, facilitan la manipulación de datos complejos.

Nextflow amplía este enfoque, agregando la capacidad de definir interacciones de programas complejas y un entorno computacional paralelo de alto nivel, basado en el modelo de programación de flujo de datos. Las características principales de Nextflow son:

- Portabilidad y reproducibilidad del flujo de trabajo
- Escalabilidad de paralelización y despliegue
- Integración de herramientas, sistemas y estándares industriales existentes

### Procesos y Canales

En la práctica, un flujo de Nextflow se realiza uniendo diferentes procesos. Cada `proceso` se puede escribir en cualquier lenguaje de programación que pueda ejecutar la plataforma Linux (Bash, Perl, Ruby, Python, etc.).

Los procesos se ejecutan de forma independiente y están aislados entre sí, es decir, no comparten un estado común (de escritura). La única forma en que pueden comunicarse es a través de colas FIFO (primero en entrar, primero en salir) asincrónicas, llamadas "canales".

Cualquier `proceso` puede definir uno o más `canales` como `entrada` y `salida`. La interacción entre estos procesos y, en última instancia, el propio flujo de ejecución del flujo, se define implícitamente mediante estas declaraciones de "entrada" y "salida".

<figure class="excalidraw">
--8<-- "docs/basic_training/img/channel-process.excalidraw.svg"
</figure>

### Abstracción de ejecución

Mientras que un 'proceso' define _qué_ comando o 'script' debe ejecutarse, el ejecutor determina _cómo_ se ejecuta ese 'script' en la plataforma de destino.

Si no se especifica lo contrario, los procesos se ejecutan en la computadora local. El ejecutor local es muy útil para desarrollar flujos y con fines de prueba; sin embargo, para los flujos computacionales del mundo real, a menudo se requiere una plataforma de alto rendimiento (HPC) o la nube.

En otras palabras, Nextflow proporciona una abstracción entre la lógica funcional del flujo y el sistema de ejecución (o runtime) subyacente. Por lo tanto, es posible escribir un flujo que se ejecute sin problemas en su computadora, un clúster o la nube, sin ser modificada. Simplemente defina la plataforma de ejecución de destino en el archivo de configuración.

<figure markdown>

![Execution abstraction](img/execution_abstraction.png)

</figure>

### Lenguaje de escritura

Nextflow implementa un DSL declarativo que simplifica la escritura de flujos de trabajo de análisis de datos complejos como una extensión de un lenguaje de programación de propósito general.

Este enfoque hace que Nextflow sea flexible: brinda los beneficios de un DSL conciso para el manejo de casos de uso frecuente con facilidad **y** da la flexibilidad y el poder de un lenguaje de programación de propósito general para manejar casos más complejos en el mismo entorno. Esto sería difícil de implementar utilizando un enfoque puramente declarativo.

En términos prácticos, las secuencias de comandos de Nextflow son una extensión del [lenguaje de programación Groovy] (https://groovy-lang.org/) que, a su vez, es un superconjunto del lenguaje de programación Java. Groovy puede considerarse como "Python para Java", ya que simplifica la escritura de código y es más accesible.

## Tu primer script

Aquí ejecutará su primer script de Nextflow (`hello.nf`), que revisaremos línea por línea.

En este ejemplo sencillo, el script toma una cadena de entrada (un parámetro llamado `params.saludo`) y la divide en partes de seis caracteres en el primer proceso. El segundo proceso convierte los caracteres a mayúsculas. El resultado finalmente se muestra en la pantalla.

### Nextflow code

<!-- NOTE: (Phil, Jan 2023)
We can dynamically include external files using mkdocs, as follows:

```groovy title="nf-training/hello.nf" linenums="1"
--8<-- "nf-training/hello.nf"
```

This inserts a code snippet identical to the one below, and we don't have to worry about keeping the two in sync.

HOWEVER - currently the line annotations cannot be added for external files. So for now, we still need to copy the scripts.

TODO: Maybe either:
    - Rewrite docs to not use loads of annotations
    - Wait for future versions to allow annotations with external files
-->

!!! info

    Click the :material-plus-circle: icons in the code for explanations.

```groovy title="nf-training/hello.nf" linenums="1"
#!/usr/bin/env nextflow
// (1)!

params.saludo = '¡Hola mundo!' // (2)!
saludo_ch = Channel.of(params.saludo) // (3)!

process SEPARALETRAS { // (4)!
    input: // (5)!
    val x // (6)!

    output: // (7)!
    path 'segmento_*' // (8)!

    // (9)!
    """
    printf '$x' | split -b 6 - segmento_
    """
} // (10)!

process MAYUSCULAS { // (11)!
    input: // (12)!
    path y // (13)!

    output: // (14)!
    stdout // (15)!

    // (16)!
    """
    cat $y | tr '[a-z]' '[A-Z]'
    """
} // (17)!

workflow { // (18)!
    letras_ch = SEPARALETRAS(saludo_ch) // (19)!
    resultados = MAYUSCULAS(letras_ch.flatten()) // (20)!
    resultados.view{ it } // (21)!
} // (22)!
```

1. El código comienza con un "shebang", que declara a Nextflow como intérprete.
2. Declara un parámetro 'saludo' que se inicializa con el valor 'Hello world!'.
3. Inicializa un `canal` etiquetado como `saludo_ch`, que contiene el valor de `params.saludo`. Los canales son el tipo de entrada para los procesos en Nextflow.
4. Comienza el primer bloque de proceso, definido como `SEPARALETRAS`.
5. Declaración de entrada para el proceso `SEPARALETRAS`. Las entradas pueden ser valores (`val`), archivos o rutas (`path`), u otros calificadores ([ver aquí](https://www.nextflow.io/docs/latest/process.html#inputs)).
6. Le dice al `proceso` que espere un valor de entrada (`val`), que asignamos a la variable 'x'.
7. Declaración de salida para el proceso `SEPARALETRAS`.
8. Le dice al proceso que espere un archivo de salida (`ruta`), con un nombre de archivo que comience con 'segmento\_\*', como salida del script. El proceso envía la salida como un canal.
9. Tres comillas dobles inician y finalizan el bloque de código para ejecutar este `proceso`.
   Dentro está el código para ejecutar: imprimir el valor de `entrada` x (llamado usando el prefijo del símbolo de dólar [$]), dividir la cadena en fragmentos con una longitud de 6 caracteres ("¡Hola " y "mundo!"), y guardar cada uno a un archivo (segmento_aa y segmento_ab).
10. Fin del primer bloque de proceso.
11. Comienza el segundo bloque de proceso, definido como `MAYUSCULAS`.
12. Declaración de entrada para el `proceso` `MAYUSCULAS`.
13. Le dice al `proceso` que espere un archivo(s) de `entrada` (`ruta`; es decir, segmento_aa y segmento_ab), que asignamos a la variable 'y'.
14. Declaración de salida para el proceso `MAYUSCULAS`.
15. Le dice al proceso que espere la salida como salida estándar (stdout) y envía esta salida como un canal.
16. Tres comillas dobles inician y finalizan el bloque de código para ejecutar este `proceso`.
    Dentro del bloque hay una secuencia de comandos para leer archivos (cat) usando la variable de entrada '$y', luego realizar la conversión a mayúsculas y generar la salida estándar.
17. Fin del segundo bloque `proceso`.
18. Inicio del alcance del flujo de trabajo donde se puede llamar a cada proceso.
19. Ejecute el `proceso` `SEPARALETRAS` en `saludo_ch` (también conocido como canal de saludo), y almacene la salida en el canal `letras_ch`.
20. Ejecute el `proceso` `MAYUSCULAS` en el canal de letras `letras_ch`, que se aplana usando el operador `.flatten()`. Esto transforma el canal de entrada de tal manera que cada elemento es un elemento separado. Almacenamos la salida en el canal `resultados_ch`.
21. El resultado final (en el canal `resultados_ch`) se imprime en la pantalla usando el operador `view` (junto al nombre del canal).
22. Fin del flujo de trabajo.

El operador `.flatten()` se usa aquí para dividir los dos archivos en dos elementos separados para pasar por el siguiente proceso (de lo contrario, serían tratados como un solo elemento).

### En la práctica

Ahora copie el ejemplo anterior en su editor de texto favorito y guárdelo en un archivo llamado `hello.nf`.

!!! advertencia

    Para el tutorial de Gitpod, asegúrese de estar en la carpeta llamada `nf-training`

Ejecute el script ingresando el siguiente comando en su terminal:

```bash
nextflow run hola.nf
```

El resultado será similar al texto que se muestra a continuación:

```linenums="1"
N E X T F L O W  ~  version 22.04.5
Launching `hola.nf` [gigantic_poitras] DSL2 - revision: 197a0e289a
executor >  local (3)
[c8/c36893] process > SEPARALETRAS (1)   [100%] 1 of 1 ✔
[1a/3c54ed] process > MAYUSCULAS (2) [100%] 2 of 2 ✔
MUNDO!
¡HOLA
```

La salida estándar muestra (línea por línea):

1. La versión de Nextflow que se ejecutó.
2. Los nombres de la secuencia de comandos y la versión.
3. El ejecutor utilizado (en el caso anterior: local).
4. El primer `proceso` se ejecuta una vez. La línea comienza con un valor hexadecimal único (consulte el CONSEJO a continuación) y finaliza con el porcentaje y la información de finalización del trabajo.
5. El segundo proceso se ejecuta dos veces (una vez para segmento_aa y otra para segmento_ab).
6. Se imprime la cadena de resultado de stdout.

!!! info

    Los números hexadecimales, como `c8/c36893`, identifican una ejecución única del proceso. Estos números son también el prefijo de los directorios donde se ejecuta cada proceso. Puede inspeccionar los archivos producidos cambiando al directorio `$PWD/work` y usando estos números para encontrar la ruta de ejecución específica del proceso.

!!! tip

    El segundo proceso se ejecuta dos veces, ejecutándose en dos directorios de trabajo diferentes para cada archivo de entrada. La salida de registro [ANSI](https://en.wikipedia.org/wiki/ANSI_escape_code) de Nextflow se actualiza dinámicamente a medida que se ejecuta el flujo; en el ejemplo anterior, el directorio de trabajo `[1a/3c54ed]` es el segundo de los dos directorios que se procesaron (sobrescribiendo el registro con el primero). Para imprimir todas las rutas relevantes a la pantalla, deshabilite la salida de registro ANSI usando el indicador `-ansi-log` (por ejemplo, `nextflow run hello.nf -ansi-log false`).

Vale la pena señalar que el proceso `MAYUSCULAS` se ejecuta en paralelo, por lo que no hay garantía de que la instancia que procesa el primer segmento (el fragmento '¡Hola ') se ejecute antes que la que procesa el segundo (el fragmento 'mundo!').

Por lo tanto, podría ser que su resultado final se imprima en un orden diferente:

```
MUNDO!
¡HOLA
```

## Modificar y reanudar

Nextflow realiza un seguimiento de todos los procesos ejecutados en su flujo. Si modifica algunas partes de su secuencia de comandos, solo se volverán a ejecutar los procesos modificados. Se omitirá la ejecución de los procesos que no se modifican y, en su lugar, se utilizará el resultado almacenado en caché.

Esto permite probar o modificar parte del script sin tener que volver a ejecutarla desde cero.

Modifique el proceso `MAYUSCULAS` en el ejemplo anterior, reemplazando el script del proceso con la cadena `rev $y`, para que el proceso se vea así:

```groovy
process MAYUSCULAS {
    input:
    path y

    output:
    stdout

    """
    rev $y
    """
}
```

Luego guarde el archivo con el mismo nombre y ejecútelo agregando la opción `-resume` a la línea de comando:


```console
$ nextflow run hello.nf -resume

N E X T F L O W  ~  version 22.04.5
Launching `hola.nf` [amazing_becquerel] DSL2 - revision: 525206806b
executor >  local (2)
[c8/c36893] process > SEPARALETRAS (1)   [100%] 1 of 1, cached: 1 ✔
[77/cf83b6] process > MAYUSCULAS (1) [100%] 2 of 2 ✔
!odnum
 aloH¡
```

Verá que se omite la ejecución del proceso `SEPARALETRAS` (el ID del proceso es el mismo que en la primera salida); sus resultados se recuperan de la caché. El segundo proceso se ejecuta como se esperaba, imprimiendo las cadenas invertidas.

!!! info

    Los resultados se almacenan en caché de forma predeterminada en el directorio `$PWD/work`. Dependiendo de su trabajo, esta carpeta puede ocupar mucho espacio en el disco. Si está seguro de que no necesitará reanudar la ejecución de su flujo, limpie esta carpeta periódicamente.

## Parámetros de la pipeline

Los parámetros se declaran simplemente anteponiendo el prefijo `params` al nombre de una variable, separados por un carácter de punto. Su valor se puede especificar en la línea de comando anteponiendo el nombre del parámetro con un carácter de doble guión, es decir, `--paramName`.

Ahora, intentemos ejecutar el ejemplo anterior especificando un parámetro de cadena de entrada diferente, como se muestra a continuación:

```bash
nextflow run hello.nf --saludo 'Bonjour le monde!'
```

La cadena especificada en la línea de comando anulará el valor predeterminado del parámetro. La salida se verá así:

```
N E X T F L O W  ~  version 22.04.5
Launching `hello.nf` [fervent_galileo] DSL2 - revision: 525206806b
executor >  local (4)
[e9/139d7d] process > SEPARALETRAS (1)   [100%] 1 of 1 ✔
[bb/fc8548] process > MAYUSCULAS (1) [100%] 3 of 3 ✔
m el r
!edno
uojnoB
```

### En formato similar a DAG

Para comprender mejor cómo Nextflow maneja los datos en este caso, a continuación se muestra una figura similar al DAG para visualizar todas las "entradas", "salidas", "canales" y "procesos":

<figure markdown>

![Hello world diagram](img/helloworlddiagram.png)

</figure>
