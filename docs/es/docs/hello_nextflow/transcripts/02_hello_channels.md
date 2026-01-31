# Parte 2: Hola Channels - Transcripción

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Traducción asistida por IA - [más información y sugerencias](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

<div class="video-wrapper">
  <iframe width="560" height="315" src="https://www.youtube.com/embed/lJ41WMMm44M?si=xCItHLiOQWqoqBB9&amp;list=PLPZ8WHdZGxmXiHf8B26oB_fTfoKQdhlik" title="YouTube video player" frameborder="0" allow="accelerometer; autoplay; clipboard-write; encrypted-media; gyroscope; picture-in-picture; web-share" referrerpolicy="strict-origin-when-cross-origin" allowfullscreen></iframe>
</div>

!!!note "Notas importantes"

    Esta página muestra únicamente la transcripción. Para instrucciones completas paso a paso, regrese al [material del curso](../02_hello_channels.md).

    Los números de sección mostrados en la transcripción se proporcionan solo con fines indicativos y pueden no incluir todos los números de sección de los materiales.

## Bienvenida

Hola, bienvenidos a la parte dos de Hello Nextflow.

Este capítulo se llama Hello Channels. Vamos a hablar sobre esta parte fundamental de Nextflow.

Los channels son las cosas que conectan los diferentes pasos en su pipeline, la forma en que sus datos y lógica fluyen a través de su flujo de trabajo.

Está bien, comencemos.

Comencemos yendo a training.nextflow.io

Hello Nextflow en la barra lateral y haciendo clic en la parte dos. Hello Channels.

Todo el material está escrito aquí abajo para que pueda seguir a su propio ritmo y capturar cualquier cosa que haya podido perder.

Una vez que tenga el sitio web abierto, puede cargar Codespaces y continuaremos desde donde estábamos al final del último capítulo.

## 0. Calentamiento: Ejecutar hello-channels.nf

Para este capítulo, vamos a editar un archivo diferente. Este se llama Hello Channels, así que puede encontrarlo en la barra lateral, haga doble clic para abrirlo.

Ahora, si acaba de venir del capítulo uno, este archivo le resultará muy familiar. El punto de partida aquí es básicamente donde terminamos el capítulo uno, con nuestro proceso llamado sayHello, nuestra entrada, salida, nuestro publishDir y nuestro params.greeting, y nuestro flujo de trabajo simple.

Estamos comenzando con un archivo nuevo, por lo que es un terreno nivelado para todos, pero puede continuar con su archivo anterior si lo prefiere.

Nota, también he eliminado todos los archivos .nextflow\* y los directorios work aquí, solo para que sea un punto de partida limpio. No importa si hace eso o no, depende de usted.

Está bien. Comencemos verificando que este pipeline todavía funciona como esperamos. Voy a mostrar el terminal aquí.

Hacer "nextflow run hello-channels.nf" y presionar enter.

Va a ejecutar ese pequeño flujo de trabajo, ejecuta nuestro paso sayHello, genera un directorio work con ese hash, y aquí está nuestra carpeta de resultados y ahí está nuestro archivo de salida, tal como esperábamos de nuestro params.greeting predeterminado.

Eso es genial. Exactamente lo mismo que el capítulo uno, funcionando como esperamos.

## 1. Proporcionar entradas variables mediante un channel explícitamente

En el capítulo uno, ya estaba usando channels, simplemente no se daba cuenta. Cuando especificamos una cadena aquí, Nextflow automáticamente creó un channel alrededor de esa cadena para nosotros, solo porque sabía que estábamos llamando a un proceso, así que necesitábamos un canal de entrada.

Lo primero que vamos a hacer es hacerlo explícito escribiendo realmente el channel en sí.

## 1.1. Crear un canal de entrada

Entonces voy a ir al workflow aquí en la parte inferior del script, y voy a decir greeting_ch. Esta es una convención que a menudo usamos en código Nextflow de tener un guión bajo ch al final de un nombre de variable cuando es un channel, solo para que sea fácil identificar que es un channel, pero no tiene que hacer eso. Igual a channel of Hello Channels.

Lo que acabamos de usar es algo llamado "Channel Factory" en el lenguaje de Nextflow. Esto es esto aquí, estamos estableciendo esta variable a un nuevo channel, y este channel factory aquí está creando un channel para nosotros de una manera particular.

Hay un puñado de diferentes channel factories que Nextflow tiene, para crear channels desde diferentes tipos de entradas. Dot of es el más simple, y simplemente toma cualquier cadena que le demos.

Note que cuando paso el cursor sobre estas palabras en VS Code, la extensión de Nextflow me está dando una ventana emergente explicando lo que hace esta sintaxis, y también hay un texto de leer más en la parte inferior de esa ventana emergente.

Si hago clic en eso, abrirá los documentos de Nextflow en una nueva pestaña y me llevará directamente a la documentación para esta cosa específica. En este caso para channel.of.

## 1.2. Agregar el channel como entrada a la llamada del proceso

Note que la extensión también nos está dando una advertencia, diciendo que hemos creado un nuevo channel aquí, pero no está siendo usado por nada.

Entonces, arreglemos eso. Voy a tomar el nuevo nombre del channel y voy a reemplazar este params.greeting con nuestro nuevo channel.

Note que ya no estamos usando la bandera de línea de comandos --greeting ahora, params.greeting no se está usando, estamos volviendo a codificar esta cadena de forma fija. Está bien. Solo estoy tratando de mantener las cosas simples. Volveremos más tarde y usaremos los params nuevamente.

## 1.3. Ejecutar el comando workflow nuevamente

Está bien, verifiquemos que esto funciona. Abro el terminal y noto de nuevo. Nextflow run hello channels. Verifico output.txt, y ahí está.

Gran ejemplo un poco aburrido, haciendo exactamente lo mismo que hicimos antes, pero ahora al menos la lógica es un poco más clara. Estamos siendo explícitos sobre escribir un nuevo channel.

Efectivamente acabamos de escribir más código para hacer lo mismo. Pero esto comenzará a tener más sentido a medida que nos volvamos un poco más complicados con la forma en que creamos nuestros channels.

## 2. Modificar el workflow para ejecutar en múltiples valores de entrada

Está bien, hagamos esto un poco más interesante. Es muy raro que quiera ejecutar un pipeline de Nextflow en una sola entrada, así que démosle varias entradas.

## 2.1. Cargar múltiples saludos en el canal de entrada

Desde los documentos aquí. Voy a copiar estas diferentes cadenas, tres de ellas. Hello, Bonjour, Olà. Oh, obtengo Esperanza. Copilot está sugiriendo un par más. Entonces tabulemos e ingresemos esos.

Los documentos de Nextflow aquí nos dicen que podemos dar múltiples valores a este operador, por lo que debería funcionar, pero probémoslo y veamos qué sucede.

## 2.1.2. Ejecutar el comando y observar la salida del registro

Bueno. Sí y no. Veamos. Dice que cinco de cinco tareas se han ejecutado aquí, pero solo nos muestra un hash, lo cual es un poco extraño. Está bien. Todo está como se esperaba aquí. Por defecto. Nextflow usa un tipo especial de salida a un terminal llamado códigos de control ANSI, lo que significa que sobrescribe ciertas líneas para dar una vista comprimida agradable de todos los diferentes procesos que se están ejecutando.

Esto tiene mucho más sentido cuando tiene flujos de trabajo más grandes y está ejecutando cientos o miles de muestras diferentes. Puede generar tanta salida en el terminal que es imposible de ver, mientras que esta vista de actualización le da un progreso en tiempo real.

## 2.1.3. Ejecutar el comando nuevamente con la opción -ansi-log false

Si lo desea, puede ejecutarlo de nuevo, y esta vez voy a usar un argumento central adicional de Nextflow con un solo guión diciendo, "-ansi-log false". Esto usa la versión anterior de la salida de registro de Nextflow. Y aquí puede ver todos los procesos individuales que se han lanzado.

Depende de usted si hace esto o no. La salida de Nextflow es exactamente la misma en ambos casos.

## 2.2. Asegurar que los nombres de los archivos de salida sean únicos

Está bien, echemos un vistazo a los archivos de salida, luego iremos a results. Pero solo hay un único archivo de salida. ¿Qué pasó? Vimos que el proceso se había ejecutado muchas veces. Podemos ir al directorio work y ver todos los diferentes hashes, todas las tareas se ejecutaron correctamente. Pero si recuerda en nuestro proceso aquí, estamos guardando todo en un archivo output.txt y luego publicando eso en este directorio.

Entonces, el mismo archivo fue creado cinco veces, y luego fue sobrescrito cinco veces. Y solo tenemos la tarea que sucedió ejecutarse al último.

## 2.2.1. Construir un nombre de archivo de salida dinámico

La forma en que arreglamos esto es usando un nombre de archivo de salida dinámico. Aquí ya tenemos una variable llamada greeting dentro del proceso, por lo que podemos usar eso en el nombre del archivo de salida. Copio eso y hago $greeting-output.txt.

Voy a rodear esto entre comillas, solo para que bash no se confunda con ningún espacio que pueda colarse aquí. Y luego voy a tomar el mismo nombre de archivo y actualizar la salida aquí.

Es realmente importante que la salida coincida con esto, porque de lo contrario, este archivo no se encontrará y Nextflow se bloqueará.

Voy a hacer una edición más realmente importante, que es que voy a cambiar estas comillas simples por comillas dobles. Note que el color del código cambió cuando hice eso. Esta variable solo se expande si usamos comillas dobles. Si uso comillas simples aquí, se usa como un valor literal, y obtendría un solo archivo llamado $greeting-output, que no es lo que quiero.

## 2.2.2. Ejecutar el workflow

Entonces pongamos las comillas dobles de nuevo y probémoslo.

Solo voy a limpiar mi directorio antes de comenzar, para que sea fácil ver los nuevos archivos. Voy a eliminar cualquier cosa llamada .nextflow, work y results.

Y voy a ejecutar ese comando Nextflow de nuevo y veamos qué archivos se crean. Entonces ejecuta los cinco procesos allí. Si estaba mirando muy de cerca, podría haber visto que la línea se actualizaba mientras se estaba ejecutando.

Y ahora podemos ir al directorio results, y efectivamente, tenemos cinco salidas diferentes, y todas tienen el prefijo del saludo diferente.

Si abro cada uno de estos, veremos que cada uno contiene el saludo correspondiente. Fantástico. Eso es lo que queremos.

## 3. Usar un operador para transformar el contenido de un channel

Está bien, entonces ahora sabemos qué son los channels y sabemos qué son los channel factories. ¿Qué hay de los operadores? Este es otro término para parte del lenguaje Nextflow, que es una serie de funciones que nos permiten operar en channels para hacer ciertas cosas con ellos. Nextflow viene con un conjunto de operadores, que nos permiten manipular channels de diversas maneras.

## 3.1. Proporcionar un array de valores como entrada al channel

Trabajemos a través de esto con un ejemplo. Digamos que queremos tomar estas cadenas de entrada, pero en lugar de simplemente ponerlas directamente en un channel factory, queremos definirlas como un array.

## 3.1.1. Configurar la variable de entrada

Entonces voy a tomar estas y hacer eso como una nueva línea arriba y decir, greetings, array.

Ahí vamos. Voy a tomar esa variable de array y ponerla en el channel.of, y presiono guardar.

## 3.1.3. Ejecutar el workflow

Ahora, veamos qué pasa. Vuelvo a mi terminal. Voy a limpiar todos esos archivos temporales de nuevo. Y ejecutemos el workflow.

No está bien. Está bien. Se rompió. Está bien. Esperaba que se rompiera esta vez. Depurar qué sale mal cuando un workflow de Nextflow falla es una parte clave de ser un desarrollador de Nextflow. Esto sucederá mucho y es importante entender qué dice el mensaje de error y cómo tratarlo.

Los mensajes de error de Nextflow son en realidad bastante estructurados. Nos dice qué proceso salió mal. Nos da un mensaje de error por una razón. Dice cuál fue el comando que intentó ejecutar dentro de esa tarea particular, cuál fue el estado de salida, cuál fue la salida en donde estaba ese directorio work de la tarea.

Note que puedo hacer clic con opción en esto en VS Code y lo abre en una barra lateral para que pueda ir directamente allí y ver todos estos archivos ocultos, de los que hablamos en el capítulo anterior, incluido el archivo .command.sh. Este puede ver que es el mismo que los comandos que se ejecutaron aquí.

Al mirar este archivo, podemos tener una idea de qué podría haber salido mal aquí en lugar de ejecutar una sola tarea para cada elemento en el array como lo hizo la última vez, simplemente proporcionó todo el array de una vez como una cadena. Entonces necesitamos desempaquetar ese array en valores individuales antes de pasarlo al channel. Volvamos y veamos si podemos hacer eso usando un operador.

## 3.2. Usar un operador para transformar el contenido del channel

En este caso, no vamos a cambiar el array antes de pasarlo al channel. Vamos a ajustar el channel para que se comporte de la manera que esperamos. Vamos a hacer eso usando el operador flatten puede hacer punto comenzar a escribir y podemos ver que la extensión de VS Code comienza a sugerir todos los diferentes operadores que tenemos disponibles.

## 3.2.1. Agregar el operador flatten()

Y voy a seleccionar flatten. Note que el espacio en blanco no importa en este contexto para Nextflow. Entonces puede poner estos operadores en una nueva línea si lo desea. Entonces puedo bajar esto aquí e indentarlo para que quede debajo de ".of" y verá que las personas a menudo encadenan muchos operadores como este en un channel e indentan de esta manera para que sea más fácil de leer.

También puede ver, como antes, puedo pasar el cursor sobre esto y leer lo que está haciendo el operador flatten, y también seguir un enlace a la documentación si quiero.

Entonces, este operador está tomando este channel, que tiene un solo array dentro de él, y separando los valores del array.

## 3.2.2. Agregar view() para inspeccionar el contenido del channel

Podemos echar un vistazo a los channels usando el operador especial view, y voy a agregar un par de ellos aquí. Esto es un poco como usar declaraciones de impresión en otros lenguajes. Entonces voy a hacer dot view y luego voy a usar estas llaves onduladas.

Esto se llama un closure. Básicamente esto da código adicional al operador view, que ejecutará en cada elemento dentro del channel. En este caso, voy a decir greeting before flatten. Greeting.

Estoy definiendo una variable aquí, que está solo dentro del alcance de este closure. Entonces esta variable solo se usa aquí y podría llamarla como quisiera. No importa realmente. Solo estoy usando greeting para que sea fácil de leer.

En algunos pipelines de Nextflow, puede ver que las personas usan una variable implícita especial llamada "$it". Así. Esta es una variable especial dentro del código Nextflow, que es una forma abreviada para que no tenga que hacer la pequeña definición de una variable. Sin embargo, con el tiempo estamos pensando que esto no es muy claro para las personas que son nuevas en Nextflow, y desalentamos el uso de "$it" ahora.

Entonces voy a mantenerme con el comportamiento anterior de greeting y usarlo así porque eso es más explícito y es más claro sobre lo que está sucediendo.

Entonces voy a copiar esta línea y hacer exactamente lo mismo de nuevo después de los argumentos flatten. El operador view es un poco especial porque hace algo en los elementos, pero también simplemente continúa pasándolos al siguiente operador para que podamos encadenarlo en medio de una cadena de operaciones como esta, e imprimirá el estado allí y seguirá adelante. Entonces, con suerte, esto nos mostrará cómo se ve el channel antes y después del operador flatten.

## 3.2.3. Ejecutar el workflow

Probémoslo. Limpio. Limpio todo en el espacio de trabajo. Ejecuto el pipeline de nuevo.

Está bien, entonces podemos ver que ejecutó nuestros cinco procesos. De nuevo, no se bloqueó con un error, así que eso es definitivamente bueno. Y ahora tenemos el before flatten y efectivamente tenemos nuestro array y tenemos after flatten, impreso cinco veces una vez para cada elemento del array. Eso es exactamente lo que esperábamos. Entonces esas son realmente buenas noticias. Y eso encaja exactamente con lo que esperaríamos del código.

Ya no necesitamos estas declaraciones de depuración, así que puedo comentarlas o eliminarlas. Voy a eliminarlas solo para mantener mi código agradable y limpio. Está bien, genial. Este ejemplo ahora está funcionando bien y podemos comenzar a ver cómo los channels pueden hacer una lógica un poco más complicada.

## 4. Usar un operador para analizar valores de entrada desde un archivo CSV

Ahora vamos a intentar hacer esto usando un archivo con una serie de entradas en su lugar. Esta es una forma muy común de escribir pipelines de Nextflow usando una hoja de muestras o un CSV de metadatos.

## 4.1. Modificar el script para esperar un archivo CSV como la fuente de saludos

Si voy a la barra lateral, puede ver greetings.csv en el repositorio de ejemplo, y este es un archivo CSV muy, muy simple que solo contiene tres líneas con tres saludos diferentes. Veamos si podemos usar este archivo CSV dentro de nuestro workflow.

Ahora voy a volver a usar params como lo hicimos en el capítulo uno, para que podamos tener una entrada de línea de comandos.

Voy a eliminar este array greetings.

## 4.1.1. Cambiar el parámetro de entrada para apuntar al archivo CSV

Voy a establecer params greeting al nombre del archivo, que es greetings.csv, y voy a usar esta variable especial para generar el channel. Voy a poner eso ahí, y los errores desaparecen. Recuerde que esto está estableciendo esta variable por defecto ahora. Entonces, si ejecuto el pipeline sin ningún argumento, usará greetings.csv, pero podría hacer --greeting para sobrescribir esta variable si quisiera.

## 4.1.2. Cambiar a un channel factory diseñado para manejar un archivo

Está bien, estamos pasando un archivo ahora en lugar de una cadena o un array de cadenas, por lo que probablemente necesitamos un channel factory diferente.

Nos vamos a deshacer de "of" que hemos estado usando hasta ahora, y en su lugar usar .fromPath. Esto hace exactamente lo que suena. Crea un channel con rutas en lugar de valores, usando un nombre de archivo de cadena o glob. También voy a eliminar el operador flatten ya que ya no necesitamos esto, ahora que estamos pasando un archivo.

## 4.1.3. Ejecutar el workflow

Voy a presionar guardar, abrir el terminal, ejecutar el workflow, y luego ver qué sucede.

Está bien. Se bloqueó de nuevo. No se preocupe. También esperaba este. Echemos un vistazo al mensaje de error y veamos si podemos averiguar qué está saliendo mal. Aquí podemos ver el comando ejecutado, y un poco como antes donde teníamos todo el array impreso. Ahora tenemos la ruta del archivo siendo repetida en el comando, en lugar de recorrer el contenido del archivo.

## 4.2. Usar el operador splitCsv() para analizar el archivo

Entonces, para usar el contenido del archivo en su lugar, necesitamos otro operador. El operador que vamos a usar para este es llamado splitCsv. Tiene sentido, porque es un archivo CSV que estamos cargando.

## 4.2.1. Aplicar splitCsv() al channel

Ok, entonces splitCsv. Cerrar paréntesis. No necesitamos ningún argumento aquí. Y de nuevo, voy a usar algunos operadores view para dar una idea de lo que está pasando aquí.

.view csv after splitCsv. Before split Cv.s

## 4.2.2. Ejecutar el workflow nuevamente

Está bien, intentemos ejecutar esto y ver qué sucede.

Está bien, tenemos un poco más de salida esta vez, pero todavía falló. Podemos mirar las declaraciones view, y aquí puede ver before split CSV, y tenemos una ruta de archivo como vimos en el mensaje de error anterior. After split CSV, ahora tenemos tres valores correspondientes a las tres líneas en el archivo CSV.

Sin embargo, puede ver que cada uno de estos valores está rodeado por corchetes. Entonces cada uno de esos era un array en sí mismo, y eso nos ha dado el mismo área que teníamos antes donde está intentando repetir un array en lugar de solo una cadena única.

Si pensamos en un archivo CSV, esto tiene sentido. Típicamente, un archivo CSV tendrá filas y columnas, por lo que split CSV hace un array bidimensional. La primera dimensión del array es cada fila, y luego hay una segunda dimensión, que es cada columna para cada fila.

Entonces aquí solo tenemos un solo valor en cada línea, por lo que tenemos una sola columna, por lo que tenemos un array de un elemento para cada línea del archivo.

Está bien. Solo necesitamos otro operador para colapsar ese array para cada línea del archivo CSV analizado. Limpiemos esto. Deshagámonos de un terminal y veamos qué podemos hacer.

## 4.3. Usar el operador map() para extraer los saludos

Ahora podríamos usar el operador flatten de nuevo, que usamos antes. Hemos visto cómo puede colapsar un array en una serie de valores, lo que funcionaría muy bien aquí. Pero voy a usar la oportunidad para demostrar otro operador, que es muy común dentro de los flujos de trabajo llamado el operador map.

## 4.3.1. Aplicar map() al channel

Voy a hacer dot map y voy a hacer item item[0].

Si escribe mucho código en otros lenguajes, puede estar familiarizado con el operador map ya. Toma un iterable, como un array o un channel, y hace alguna operación en cada valor de eso.

Aquí estamos diciendo que deberíamos definir una variable llamada item dentro del alcance de este closure, y luego queremos devolver, solo el primer valor en ese array. Entonces item índice cero.

Esto está efectivamente aplanando el array. Puede ver cómo podríamos extender esto para ser más complejo, sin embargo: si nuestro archivo CSV tuviera seis columnas, pero solo estamos interesados en la cuarta columna, podríamos acceder a un índice específico aquí. O hacer cualquier otro tipo de operación en el valor antes de pasarlo al procesamiento posterior.

Entonces el operador map es extremadamente flexible y muy poderoso para modificar channels en vuelo. Pongamos otra declaración view solo para que podamos ver lo que está haciendo en nuestra ejecución. Puedo adjudicar esa línea y moverla hacia abajo. Y after map.

## 4.3.2. Ejecutar el workflow una vez más

Abramos el terminal e intentemos ejecutar el workflow.

Está bien, no hay errores esta vez. Esa es una buena señal. Ahora podemos revisar todas estas diferentes salidas de las declaraciones view. Before split CSV, teníamos una sola ruta. After split CSV, teníamos los arrays de un solo valor, y luego after map, tenemos solo los valores sin ninguna sintaxis de array. Vayamos al directorio results, y aquí están nuestros archivos de salida comportándose exactamente como queríamos.

Hay un pequeño bonus aquí. Puede ver realmente que los operadores view están ligeramente mezclados en el orden en que han hecho la salida. Esto se debe a que Nextflow está haciendo la paralelización de estas diferentes tareas. Entonces, después de que dividió el CSV, hay tres elementos en este channel, y está manejando el procesamiento de esos tres elementos en paralelo automáticamente. Eso significa que el orden de las salidas es estocástico y puede variar. En este caso, simplemente sucedió que algunos de los operadores view regresaron después de que el paso subsiguiente se había completado, y por eso vino en este orden.

Si ejecuto el mismo workflow de nuevo. Entonces, efectivamente, ha venido en un orden diferente y esta vez tenemos los split CSV y los maps en el orden que esperaríamos.

Entonces solo tenga en cuenta, no puede confiar en el orden de salidas de una tarea de proceso porque Nextflow está manejando esta paralelización para usted automáticamente. Nextflow hace eso por usted con su lógica de flujo de datos, y ese es el verdadero poder de Nextflow.

Está bien, este es probablemente uno de los capítulos más importantes de todo el entrenamiento. Una vez que comprenda los channels, los channel factories y los operadores, comenzará a entrar en la fortaleza de Nextflow y lo que lo hace único como lenguaje de programación. Esta funcionalidad permite a Nextflow paralelizar todos sus flujos de trabajo para usted y generar lógica de flujo de trabajo extremadamente compleja con una sintaxis muy limpia y un modelo de flujo de datos de empuje. Puede ser un concepto un poco extraño al principio, pero una vez que se acostumbre a escribir código como este, rápidamente se sentirá natural y antes de que se dé cuenta, estará escribiendo flujos de trabajo fantásticos.

Tome un descanso, una taza de té, camine y pasemos al capítulo tres, donde comenzamos a extender estos conceptos a flujos de trabajo más complejos. Nos vemos en el próximo video.

[Transcripción del siguiente video :octicons-arrow-right-24:](03_hello_workflow.md)
