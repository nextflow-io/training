# Parte 3: Hola Workflow - Transcripción

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Traducción asistida por IA - [más información y sugerencias](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

<div class="video-wrapper">
  <iframe width="560" height="315" src="https://www.youtube.com/embed/zJP7cUYPEbA?si=Irl9nAQniDyICp2b&amp;list=PLPZ8WHdZGxmXiHf8B26oB_fTfoKQdhlik" title="YouTube video player" frameborder="0" allow="accelerometer; autoplay; clipboard-write; encrypted-media; gyroscope; picture-in-picture; web-share" referrerpolicy="strict-origin-when-cross-origin" allowfullscreen></iframe>
</div>

!!!note "Notas importantes"

    Esta página muestra únicamente la transcripción. Para instrucciones completas paso a paso, regrese al [material del curso](../03_hello_workflow.md).

    Los números de sección mostrados en la transcripción se proporcionan solo con fines indicativos y pueden no incluir todos los números de sección de los materiales.

## Bienvenida

Hola, bienvenidos a la parte tres del curso de entrenamiento "Hola Nextflow".

Este capítulo se llama "Hola Workflow".

En el capítulo dos, construimos un flujo de trabajo simple de un proceso, pero en realidad, los pipelines son útiles porque pueden encadenar múltiples pasos de análisis juntos.

En este capítulo, vamos a tomar ese ejemplo inicial y ampliarlo para que sea un poco más realista.

Vamos a agregar algunos pasos adicionales y vamos a ver cómo usamos los canales para conectar esos pasos.

Vamos a ver múltiples tareas, que pueden colapsar en un solo proceso y vamos a ver procesos que pueden tener múltiples entradas y múltiples salidas.

Bien, comencemos.

Así que empecemos. Igual que antes. Vamos a training.nextflow.io. Hola Nextflow, capítulo tres. Hola Workflow. Y abramos nuestro espacio de trabajo. He limpiado todos mis archivos de trabajo de mis capítulos anteriores y voy a abrir Hola Workflow.

Ahora este es el mismo archivo en el que hemos estado trabajando hasta ahora, así que debería verse familiar. Tenemos nuestro proceso say hello. Tenemos nuestro params.greeting con su archivo greetings CSV, y tenemos nuestro workflow en la parte inferior, que carga ese archivo CSV, crea el canal y lo pasa a nuestro proceso.

## 0. Calentamiento: Ejecutar hello-workflow.nf

Si lo desea, podemos probar esto y verificar dos veces que esté funcionando como esperamos. Abra una terminal para nextflow run hello workflow nf y presione enter.

Bien, genial. Nuestros tres procesos se ejecutan. Tenemos nuestro directorio de resultados con nuestras tres salidas. Bonjour. Hello. Holà. Así que cerremos esos archivos, cerremos la terminal, volvamos al script.

## 1. Agregar un segundo paso al flujo de trabajo

Bien. Para nuestro ejemplo, nos mantenemos básicos e intentamos mantenernos agnósticos del dominio. Así que nuestro segundo proceso solo va a manipular estas cadenas, estas palabras, de una manera simple. Vamos a usar el comando Unix translate para tomar estos archivos y ponerlos todos en mayúsculas. Lo hacemos con el comando "tr".

## 1.1. Definir el comando de conversión a mayúsculas y probarlo en la terminal

Podemos probar esto solo en la terminal bash, y ver si funciona. Así que haces echo, Hello World, y luego pasas eso con el carácter pipe a tr, y le damos un patrón de reconocimiento, a a z y a qué debería traducirse. A a Z en mayúsculas.

Esto es muy simple porque está literalmente haciendo los caracteres de A a Z. Así que no funcionará con nada que tenga acentos o algo así. Pero para los propósitos del ejemplo, debería entender la idea.

Voy a presionar enter y se imprime en una terminal, HELLO WORLD en mayúsculas. Y como antes, podríamos redirigir esto a un archivo si quisiéramos. Outfile.

Bien. Limpiemos esto.

## 1.1. Escribir el paso de conversión a mayúsculas como un proceso de Nextflow

Volvamos a nuestro script y escribamos un nuevo proceso para manejar este comando bash. Voy a copiar el proceso anterior, pegarlo debajo, y llamarlo convert to upper. Para mayúsculas. Voy a usar el mismo publishDir results, pero voy a hacer algunos cambios aquí. En lugar de tomar un val, voy a tomar un path input file, y voy a tener un prefijo aquí upper, para que nuestros archivos de salida no sobrescriban la salida. Y voy a usar el nombre de variable de la entrada. Y luego voy a cambiar un script aquí abajo, y en su lugar voy a usar cat en el archivo de entrada y al igual que hicimos en Bash TR, a-z, upper input file .txt. Bien, hagamos clic en guardar.

## 1.2. Agregar una llamada al nuevo proceso en el bloque workflow

Ahora si me desplazo hacia abajo, necesitamos realmente llamar a este proceso. Solo agregar el proceso al script no es suficiente. Tenemos que decirle a Nextflow que necesitamos ejecutar este proceso y dónde hacerlo.

Así que voy a hacer aquí, convert to upper y

bien, estamos obteniendo un error aquí que dice que espera un argumento. Por supuesto, necesitamos pasar algo a este proceso para que realmente tenga algo que hacer.

## 1.3. Pasar la salida del primer proceso al segundo proceso

Lo que vamos a hacer es vamos a tomar la salida de este proceso. Así que tomo el nombre, say hello, y cuando hago dot out.

Para un ejemplo simple como este, donde tenemos un proceso que tiene solo una salida y estamos pasando eso a un nuevo proceso, por lo que tiene una entrada, eso debería ser todo lo que necesitamos. Así que voy a hacer clic en guardar, abrir la terminal, y tratemos de ejecutar esto nuevamente.

## 1.4. Ejecutar el flujo de trabajo nuevamente

Ahora, no he limpiado mi directorio de trabajo de la última vez que ejecuté este flujo de trabajo. Voy a ejecutarlo nuevamente y voy a usar esto como una oportunidad para mostrar cómo funciona el almacenamiento en caché parcial. Así que si hago single dash resume. Con suerte debería reutilizar las salidas de ese primer proceso, que eran exactamente las mismas que la última vez que ejecuté. Pero ahora tenemos un nuevo proceso aquí que no se ha ejecutado antes, que se ejecuta desde cero. Y efectivamente, puede ver que el primer proceso usó las salidas de caché, y la segunda salida ejecutó tres de tres. También puede ver que tenemos ambos procesos aquí ahora, nuestro primer proceso, say hello, se ejecutó tres veces, y nuestro segundo proceso convert to upper se ejecutó tres veces.

Si ejecuto esto nuevamente, como recordatorio, con -ansi-log false, deberíamos ver que se ejecutaron seis tareas de proceso diferentes, tres para cada una de ellas. Así que esto está haciendo exactamente lo que esperábamos. El primer proceso se está ejecutando tres veces, pasando esas salidas a un segundo proceso, que luego se está ejecutando tres veces.

Así que echemos un vistazo dentro del directorio de trabajo y veamos cómo Nextflow está manejando estas entradas de archivos. Si tomo este directorio hash aquí del segundo proceso, podemos usar un comando tree nuevamente con -a solo para mirar estos archivos. Puede ver aquí que tenemos nuestro archivo de entrada, que es el archivo Bonjour-output.txt, y eso es en realidad un enlace simbólico. Eso es lo que nos muestra esta flecha, y está apuntando al archivo en el directorio de trabajo anterior.

Esto tiene sentido. Nextflow maneja la ejecución de cada tarea en su propio directorio encapsulado, por lo que está completamente autocontenido. Sin embargo, necesita proporcionar los archivos de pasos anteriores como entrada. En lugar de salir fuera del directorio de trabajo para obtener esos archivos, Nextflow los prepara en el directorio de trabajo.

Si tenemos un sistema de archivos compartido como aquí, hace eso usando un enlace simbólico para que no use ningún espacio de archivo adicional. Si usamos almacenamiento en la nube con buckets en diferentes ubicaciones, buscaría esos archivos y realmente los copiaría en el directorio de trabajo.

Echemos un vistazo al archivo command sh. Si hago code work, command sh, puede ver, efectivamente, está accediendo a ese archivo desde el directorio local. Así que todo está muy autocontenido y limpio.

También podemos verificar el directorio de resultados y asegurarnos de que estos archivos se hayan producido correctamente. Y efectivamente, en resultados, podemos ver todos los archivos de salida del primer proceso y todos los archivos de salida del segundo. Y todos están en mayúsculas como esperábamos.

Aquí es donde el poder de Nextflow comienza a brillar. Con un código muy mínimo y Nextflow manejó la ejecución en paralelo de estas tareas con una encapsulación limpia dentro de directorios de trabajo separados y preparación de archivos de entrada y salida y publicación de archivos, todo automáticamente para nosotros, justo desde el principio. Así que puede ver cómo, a medida que escalamos esta complejidad de nuestros flujos de trabajo de análisis, esta funcionalidad es realmente, realmente valiosa.

## 2. Agregar un tercer paso para recolectar todos los saludos

Bien. Estos pasos fueron uno a uno. Tuvimos una salida del primer proceso yendo a una entrada para el segundo proceso. A continuación, vamos a hablar sobre cómo recolectar estas diferentes salidas en una sola tarea de proceso, lo cual es nuevamente algo muy común de hacer. Así que rápidamente abramos la terminal y hagamos una ejecución de prueba de esto.

## 2.1. Definir el comando de recolección y probarlo en la terminal

Voy a hacer trampa y copiar el código bash de ejemplo del material de entrenamiento y simplemente presionar enter.

Lo que podemos ver aquí es que ejecutamos este comando echo tres veces en tres archivos de salida diferentes, que puedo ver aquí. Y luego usamos el comando cat para imprimir la salida de cada uno de estos tres archivos diferentes, y redirigir eso a un único archivo recolectado.

Y si hago "cat COLLECTED-output", puede ver que tiene el contenido de esos tres archivos diferentes, ahora en un solo archivo.

## 2.2. Crear un nuevo proceso para hacer el paso de recolección

Así que veamos si podemos replicar lo mismo dentro de nuestro pipeline de Nextflow.

Desplacémonos hacia arriba y creemos un tercer proceso. Voy a copiar este anterior, y esta vez lo voy a llamar Collect Greetings.

En la terminal bash, lo llamamos collected output txt. Así que voy a decir lo mismo path output aquí. Y voy a hacer la redirección aquí, para que se guarde de la misma manera.

Bien. Necesitamos cambiar lo que sucede al principio de ese comando, y necesitamos pensar en cuál es el archivo de entrada aquí. De hecho, este proceso va a tomar múltiples archivos de entrada. Voy a mantener path y voy a cambiar esto a una nueva variable llamada input files, en plural.

Luego voy a nuevamente, hacer cat con ellos como hicimos en nuestro script bash. Y voy a usar la variable aquí.

Ahora, podrías pensar que esto no funcionaría. Hemos visto fallas anteriormente donde se pasó un arreglo de cadenas o un arreglo de rutas a un proceso y eso causó un error. Pero de hecho, aquí Nextflow va a manejar esto automáticamente para nosotros de la manera correcta. Va a tomar varios archivos de entrada diferentes, y solo va a imprimir las diferentes rutas de archivos aquí.

Por supuesto, ayuda que el comando cat pueda tomar una serie de nombres de archivos como este. Si estuviera usando un comando diferente que requiriera un argumento antes de cada ruta de archivo o algo así, tendríamos que tener un poco más de código aquí y lógica para poder manejar la iteración de estas rutas de archivo. Pero en este caso, debería funcionar simplemente.

## 2.3. Agregar el paso de recolección al flujo de trabajo

Bien, bajemos al workflow y agreguemos nuestro nuevo proceso. Collect greetings. Y nuevamente, tomemos la salida de convert to upper out. Guardemos esto.

Démosle una oportunidad. nextflow run hello workflow.

Bien, el flujo de trabajo se ejecutó, pero algo es un poco extraño aquí. Tenemos tres ejecuciones del primer paso, lo cual esperamos. Tres tareas para el segundo, pero también tenemos tres tareas al final cuando esperábamos tener solo una sola tarea aquí fusionando todas las salidas.

Si vamos a nuestro directorio de resultados. También vemos que la salida recolectada solo tiene un valor único en lugar de los tres. Esto es porque ese archivo de salida fue sobrescrito tres veces con tres valores diferentes.

Esto tiene sentido porque pasamos una salida a una entrada aquí de la misma manera que hicimos en el paso anterior.

## 2.4. Usar un operador para recolectar los saludos en una sola entrada

Así que necesitamos un operador aquí para tomar este canal con tres elementos y colapsarlos a un solo elemento, para que ese proceso final solo se ejecute una vez.

Para hacer eso, vamos a usar el operador collect. Puedo hacer esto directamente dentro del workflow. Puedo hacer .out y encadenar un operador aquí al final .collect.

Presionar guardar. Y luego para los propósitos de este entrenamiento, también voy a hacer algunos operadores view como hicimos antes, para que podamos echar un vistazo a este canal antes y después de usar el operador collect, para que podamos entender qué está sucediendo.

Voy a tomar este canal, deshacerme del collect y dot view greetings, y luego voy a duplicar esta línea, agregar el operador collect. Y cambiar eso a after.

Esto es separado de donde estamos llamando esto, pero eso está bien porque estamos usando las mismas llamadas de operador en el mismo canal de salida.

Bien, guardemos y probémoslo en la terminal. Voy a ejecutar nextflow run. Hello, workflow. Volver a ejecutar nuestro script.

Bien. Esto se ve mejor. Como antes, podemos ver que los dos primeros procesos se ejecutan tres veces y ahora nuestro proceso final solo se ejecutó una vez.

Si miramos lo que fue impreso por el operador view, aquí abajo, dijimos before collect, que es esta salida aquí, y eso se imprimió tres veces. Y puede ver que hay una sola ruta para cada uno de esos. Y luego after collect, puede ver que tenemos este arreglo de tres rutas. Así que eso es como esperamos.

Bien, revisemos el archivo de resultados y veamos si es lo que esperamos esta vez. Efectivamente, ahora hay tres líneas en el archivo - eso concatenó exitosamente estas tres salidas en un solo archivo de salida. Fantástico.

Bien, voy a limpiar y pasemos al siguiente paso. Y voy a eliminar estas declaraciones view solo para mantener las cosas limpias.

## 3. Pasar más de una entrada a un proceso para nombrar el archivo de salida final de forma única

Bien. Hasta ahora, todos nuestros procesos solo han tomado una sola entrada. Ahora vamos a hacer un ejercicio donde agregamos más de una entrada a un proceso para ver cómo funciona esto. Para hacer esto, vamos a usar este ejemplo collect greetings.

Cada vez que ejecuté el flujo de trabajo, sobrescribió ese archivo en el directorio de resultados, lo cual puede no ser lo que queremos.

## 3.1. Modificar el proceso recolector para aceptar un nombre definido por el usuario para el archivo de salida

Así que para este ejemplo, vamos a pasar un parámetro adicional para que podamos personalizar el nombre del archivo de salida.

Agregar una segunda entrada a un proceso es muy simple. Solo agrego una segunda línea en el bloque input. Esta vez va a ser un valor, en lugar de una ruta, porque queremos pasar una cadena y voy a llamarlo batch underscore name.

Ahora puedo usar esta variable en el bloque script, y voy a decir collected dash dollar batch name.

Estoy usando llaves aquí alrededor del nombre de la variable. Eso es solo para mantenerlo separado del resto de una cadena, y probablemente no sea necesario en este caso, pero creo que hace que sea más fácil de leer.

Bien. Finalmente, recuerde actualizar la ruta de salida porque ahora el nombre del archivo ha cambiado, así que voy a hacer lo mismo y poner el batch name en la salida de path como se esperaba.

## 3.2. Agregar un parámetro batch de línea de comandos

Ahora necesitamos pasar un nombre de lote desde algún lugar, y voy a crear un segundo parámetro para hacer esto para que podamos hacerlo en la línea de comandos cuando ejecutemos el flujo de trabajo.

Así que voy a hacer params batch name, y por defecto, llamemos a esto test batch. Ahora puedo usar esta variable de parámetros especiales abajo, donde llamamos al proceso.

Y efectivamente VS Code nos está diciendo que no hay suficientes argumentos para este proceso ahora, y que espera una segunda entrada.

Simplemente hago coma y paso nuestra nueva variable y el error desaparece.

Note que el orden de las entradas aquí es realmente importante. La primera entrada del proceso era la ruta, y la segunda entrada es el nombre. Si cambio el orden aquí, también debo cambiar el orden cuando llamo al proceso. De lo contrario. Siguiente, pasaremos el canal equivocado a la entrada equivocada.

## 3.3. Ejecutar el flujo de trabajo

Bien, probémoslo y veamos si funciona. Hagamos "nextflow run hello- workflow. Bien, se ejecutó como antes. Echemos un vistazo al directorio de resultados.

Efectivamente, nuestro nombre de archivo aquí ahora se llama "collected test batch output txt". Fantástico.

Y ahora veamos si podemos sobrescribir eso ejecutando nuevamente. Esta vez voy a hacer --batch_name para que coincida con ese nombre de variable de parámetro especial aquí. Y voy a llamarlo demo output.

Ejecute el flujo de trabajo nuevamente y veremos si algo sucede.

Bien, ahora tenemos un collected demo output .txt. Y debido a que este nombre de archivo es diferente a ese, no lo sobrescribió. Ambos están ahora presentes en el directorio de resultados.

## 4. Agregar una salida al paso recolector

Bien, así que ahí mostramos dar múltiples entradas a un proceso, pero ¿qué hay de múltiples salidas? Para este ejemplo, vamos a calcular el número de saludos que se procesan y producir eso como una salida secundaria para este paso collect greeting.

## 4.1. Modificar el proceso para contar y producir el número de saludos

Vamos a hacer un poco de truco aquí. Los procesos de Nextflow tienen este bloque script con una cadena multilínea, y eso se pasa como salida bash al dot comando dot sh. Pero en realidad podemos escribir cualquier código personalizado arriba de eso, y eso se ejecutará como parte de una tarea pero no se incluirá dentro del script bash.

Una de las funciones incorporadas en la sintaxis de Nextflow se llama size. Así que voy a tomar la entrada de ruta, y voy a decir count underscore greetings, solo para definir un nombre de variable. Voy a tomar los archivos de entrada y voy a llamar "size" en él.

Esta función contará el tamaño de este canal de entrada y lo asignará a una variable.

Ahora podemos devolver esa variable como parte del bloque output. Así que decimos, val, porque es un valor, no un archivo. Y count greetings.

Ahora esto es suficiente por sí mismo, y ahora podríamos acceder a estas diferentes salidas de este proceso. Sin embargo, tendríamos que acceder a ellas de manera posicional. Así que usando una clave de índice como cero y uno.

Para que sea un poco más fácil obtener las salidas, podemos nombrarlas y hacemos eso usando una declaración emit.

Así que hacemos coma emit out file o como quiera llamar a esto. Y hago aquí emit count. Esto es básicamente solo un decorador, que solo nos ayuda a escribir un código un poco más limpio para que podamos fácilmente referenciar las salidas específicas más adelante en el bloque workflow.

## 4.2. Reportar la salida al final del flujo de trabajo

Bien. Si me desplazo hacia abajo al bloque workflow, ahora puedo tomar las salidas de collect greetings, hacer collect greetings, dot out, y podemos ver nuestras dos salidas nombradas se sugieren aquí por la extensión de VS Code. Muy útil.

Así que voy a hacer dot count para obtener el valor de conteo que acabamos de crear, y voy a hacer view, para que se imprima en la línea de comandos. Para que podamos verlo cuando ejecutemos el flujo de trabajo.

Escribamos algo en el closure aquí solo para hacerlo un poco más agradable. num greetings, there were greetings greetings.

Y en realidad no nos importa la otra salida porque no la estamos usando como entrada para ningún otro proceso. Pero puede ver cómo podríamos fácilmente pasar esto como entrada a otro proceso si quisiéramos, posteriormente.

## 4.3. Ejecutar el flujo de trabajo

Vamos a hacer clic en guardar. Echemos un vistazo a la terminal y probémoslo.

Bien, fantástico. Aquí vamos. Hay tres saludos. Eso es exactamente correcto.

Bien, gran trabajo. Ese es el final de este capítulo. Hemos terminado por llegar hasta aquí. Ahora está comenzando a construir un flujo de trabajo bastante realista, donde somos capaces de manejar entradas y salidas y lógica dentro de nuestro flujo de trabajo.

A medida que estos archivos de flujo de trabajo se hacen más largos, comienzan a volverse un poco difíciles de manejar. Así que en el próximo capítulo, veremos cómo podemos modularizar el código de Nextflow en archivos separados para que sea más fácil encontrar y mantener el código dentro del flujo de trabajo.

Únase a nosotros en el próximo video para el capítulo cuatro. Hola Modules.

[Siguiente transcripción de video :octicons-arrow-right-24:](04_hello_modules.md)
