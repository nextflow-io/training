# Parte 1: Hola Mundo - Transcripción

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Traducción asistida por IA - [más información y sugerencias](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

<div class="video-wrapper">
  <iframe width="560" height="315" src="https://www.youtube.com/embed/8X2hHI-9vms?si=F0t9LFYLjAWoyRXj&amp;list=PLPZ8WHdZGxmXiHf8B26oB_fTfoKQdhlik" title="YouTube video player" frameborder="0" allow="accelerometer; autoplay; clipboard-write; encrypted-media; gyroscope; picture-in-picture; web-share" referrerpolicy="strict-origin-when-cross-origin" allowfullscreen></iframe>
</div>

!!!note "Notas importantes"

    Esta página muestra solo la transcripción. Para instrucciones completas paso a paso, regrese al [material del curso](../01_hello_world.md).

    Los números de sección mostrados en la transcripción se proporcionan solo con fines indicativos y pueden no incluir todos los números de sección en los materiales.

## Bienvenida

Hola, bienvenidos al Capítulo Uno de Hola Nextflow.

En esta primera parte de un curso de seis partes, vamos a ver lo más básico de Nextflow. Comenzaremos ejecutando algunos comandos en una terminal, y luego tomaremos esos comandos Bash y veremos cómo construirlos en un script de Nextflow.

Intentaremos ejecutar ese primer pipeline de Nextflow, veremos qué hace Nextflow, dónde se ejecuta, qué archivos crea, y cuál es el propósito de esos archivos.

Muy bien, comencemos.

## training.nextflow.io

Primero que nada, vayan a training.nextflow.io. Al igual que antes, todo el material está escrito aquí, y lo trabajaré paso a paso. Mostraré mi pantalla mientras realizo los pasos del entrenamiento, pero todo lo que estoy diciendo está en el material de entrenamiento para que puedan seguirlo a su propio ritmo, y pueden encontrarlo todo escrito allí.

Este video también tiene subtítulos habilitados, así que siéntanse libres de activarlos y seguir exactamente lo que digo mientras lo digo.

Bien, vayamos a Hola Nextflow. Ese es el curso que vamos a hacer hoy, y ya hicimos la orientación en el primer video, así que iremos directamente a la parte uno. Hola Mundo.

Bien, voy a dejar este material de entrenamiento ahora y saltar a mi entorno de Code Spaces. Esto es lo que configuramos en el primer video. Espero que tengan algo que se vea muy similar a esto en su propio sistema. Estoy usando VS Code y estoy viendo el material de entrenamiento y he cambiado de directorios al directorio hello Nextflow.

## 0. Calentamiento: Ejecutar Hola Mundo directamente

Bien. Comencemos con un par de conceptos básicos, que espero sean familiares para todos. Voy a comenzar simplemente escribiendo un comando muy básico en la terminal. Aquí abajo voy a decir 'echo Hola Mundo!' presiono enter y, sin sorpresas, la terminal hace lo que le pido y devuelve esa cadena. Hola mundo.

Bien, luego voy a presionar arriba para obtener ese comando y editarlo un poco más. Esta vez vamos a redirigir esa salida a un archivo. Voy a escribirlo en output.txt y presionar enter, nada en la terminal esta vez porque la salida no vino a la terminal. Fue a ese archivo.

Luego puedo leer ese archivo haciendo 'cat output.txt' presiono tab allí para expandir automáticamente el nombre del archivo y ahí está. El archivo está ahí.

También puedo ver ese archivo en la barra lateral en el explorador de archivos en VS Code. Puedo hacer doble clic en él y abrirlo aquí. Si quieren abrirlo en VS Code sin hacer clic en nada, también pueden hacer "code" y luego "output.txt" y hace lo mismo.

Genial. Ese es el primer paso. Muy simple.

## 1. Examinar el script inicial del workflow Hola Mundo

Bien. Ahora vamos a hacer exactamente lo mismo, pero en Nextflow, en lugar de directamente en la terminal.

Vamos a usar el primer script de ejemplo para comenzar, este archivo se llama Hello World. Puedo hacer "ls" para verlo en una terminal, y estoy en Mac, así que puedo hacer command clic para abrir ese archivo, o podría haber simplemente hecho doble clic en la barra lateral aquí.

Hay algunas cosas que podemos ver en este archivo. Justo arriba, hay una declaración hash que dice que este es un archivo Nextflow y así es como podría ejecutarse. Hay algunos comentarios aquí, solo comentarios regulares de código en gris claro, que no afectan la ejecución, y solo nos ayudan a leer el script.

Y luego hay dos estructuras principales. Hay un process aquí y un workflow.

Los processes en Nextflow son los pasos del pipeline. Son las partes que realmente hacen la lógica y realizan el procesamiento.

El workflow entonces en la parte inferior une estos procesos juntos y gobierna la lógica del workflow, cómo todo se conecta entre sí.

Vamos a comenzar viendo un process. Volveremos al workflow en un momento.

## 1.2 La definición del process

Entonces cada process comienza con una palabra clave process. Tiene un nombre y luego tiene algunos corchetes y todo dentro de esos corchetes es ese único process.

Un process debe tener una sección script, y contenido aquí hay un fragmento de bash en una cadena multilínea, que es la parte del código que realmente se ejecuta en el entorno de cómputo.

También tenemos una declaración output aquí, que le dice a Nextflow, qué archivos se espera que sean creados por el script. Note que el output aquí tiene una palabra clave path, que le dice a Nextflow que esto es un archivo, no un valor, o una cadena.

Dentro del bloque script, esto es solo una declaración bash regular, y es exactamente lo mismo que escribimos en la terminal. Estamos haciendo echo de hola mundo a un archivo llamado output.txt. Este output.txt es luego recogido por la definición de output. La definición de output en realidad no está haciendo nada. Solo está diciéndole a Nextflow qué esperar, y si este archivo no se creara, Nextflow lanzaría un error.

Note que este ejemplo no es muy bueno porque hemos codificado el nombre del archivo aquí, output.txt y output.txt. Si cualquiera de estos se cambiara, eso causaría un error en nuestro workflow.

Hay una mejor manera de hacer esto con variables, que cubriremos en un minuto.

## 1.3 La definición del workflow

Bien. Bajando al workflow, podemos ver que tenemos un comentario y luego ejecutamos el process llamado sayHello. Esta es la misma palabra clave que está aquí arriba. Esto es tan simple como puede ser un workflow. Solo estamos llamando a un único process sin entrada variable, así que no lo estamos conectando a nada más. En la parte posterior de este curso, hablaremos sobre cómo hacer esto más poderoso usando entradas variables y conectando cosas con channels.

## 2. Ejecutar el workflow

Bien, esto es todo lo que necesitamos. Veamos si podemos ejecutarlo y ver qué sucede. Voy a simplemente limpiar la terminal y luego voy a hacer "nextflow run", y voy a llamar al nombre del archivo, que es hello-world.nf. Eso es todo lo que necesitamos para ejecutar un pipeline de Nextflow. Este pipeline no toma ninguna entrada, así que no necesitamos ningún otro argumento.

Presionemos enter y veamos qué pasa.

Bien. Esperemos que deberían tener alguna salida, que se vea así. Tenemos algunos fragmentos de información que nos dicen que Nextflow se ejecutó y qué versión estaba usando. Nos dice qué script se lanzó y nos da un nombre generado aleatoriamente para esta ejecución particular del workflow. En este caso, el mío se llamó "gloomy_crick".

La parte más importante de esto, sin embargo, es que nos dice qué pasos se ejecutaron en el pipeline. Pueden ver que nuestro process llamado sayHello se ejecutó, y se ejecutó una vez y estuvo cien por ciento completo.

Esta parte aquí es el hash para esa tarea particular del workflow. Cada process se ejecuta una o más veces, y cada una de esas ejecuciones se llama una tarea.

## 2.2. Encontrar la salida y los logs en el directorio work

Cada tarea obtiene su propio directorio aislado donde se ejecuta, por lo que está separado del resto de la ejecución del workflow. Este hash corresponde a la estructura de archivos dentro del directorio work. Si hago "tree work", podemos ver a0, y luego una versión más larga del hash corto, y luego nuestro archivo output.txt. También pueden verlo en una barra lateral.

Pueden ver en la barra lateral que hay algunos archivos adicionales aquí. La razón por la que estos no aparecieron en una terminal es porque son archivos ocultos, comienzan con un punto. Y de hecho, si hago "tree -a" para todos, y "work", podemos verlos aquí.

Estos archivos con punto están presentes en cada directorio work que Nextflow crea, y cada uno tiene una tarea ligeramente diferente. Primero .command.begin simplemente incluye algunas instrucciones para Nextflow que configuran la tarea antes de que se ejecute. .command.run son las instrucciones reales ejecutadas por Nextflow mismo. Luego .command.sh es probablemente el más interesante. Este es el script que fue resuelto desde nuestro bloque script del process.

Si lo abro, pueden ver que tenemos nuestro "echo Hello World" al archivo output.txt. Esto es exactamente lo mismo que nuestro process en este caso, pero si tenemos alguna variable dentro de nuestro código Nextflow, cada tarea tendrá un .command.sh diferente, y pueden ver cómo se resolvieron esas variables.

Los otros archivos tienen que ver con cómo se ejecutó la tarea. Entonces .command.err, .log y .out son el error estándar, salida estándar y los dos combinados. Y .exitcode le dice a Nextflow cómo se ejecutó esta tarea con qué código de salida, si fue exitosa o no.

Finalmente, tenemos nuestro archivo output.txt y efectivamente, "Hello World" esto es lo que estábamos esperando y esto es lo que se creó.

Bien, genial. Esa fue su primera ejecución de Nextflow. Felicitaciones. Realmente es así de simple.

A continuación, vamos a ver cómo hacer esto un poco más conveniente para que no tengamos que editar el código cada vez que queramos hacer un cambio en cómo se ejecuta el pipeline.

## 3. Gestionar ejecuciones del workflow

Esta estructura de directorios es genial para mantener todas las tareas separadas y todo organizado, pero por supuesto, no es muy conveniente para encontrar sus archivos de salida. No quieren estar excavando a través de montones de directorios anidados tratando de encontrar los resultados de su pipeline.

## 3.1. Publicar salidas

La buena noticia es que no se supone que lo hagan. Los directorios work son realmente solo para que Nextflow los use. Así que lo que vamos a hacer es usar una función para Nextflow llamada "publishDir".

Volvemos a nuestro workflow, vamos al process. Podemos agregar una nueva declaración aquí llamada directiva. Esto es lo que Nextflow llama estas cosas en la parte superior de los procesos que aumentan cómo funciona la funcionalidad, y la que vamos a usar se llama publishDir.

Pueden ver que he comenzado a escribir aquí y la extensión de Nextflow para VS Code me ha sugerido la directiva, así que solo puedo presionar enter.

Bien. Voy a seguir esto con un directorio llamado "results" y vamos a decirle que copie los archivos de salida allí. Así que voy a decir mode copy. Genial. Voy a guardar y ejecutemos el workflow nuevamente.

nextflow run hello-world.nf

Se ejecuta exactamente igual. Aunque note que tenemos un hash ligeramente diferente esta vez. Nextflow, usará un hash diferente cada vez que ejecute el workflow. Y tenemos un conjunto diferente de directorios work como resultado. Áreas, uno llamado EB en su lugar, pero pueden ver que todos los archivos son los mismos. Sin embargo, lo que es nuevo esta vez es que también tenemos un directorio llamado "results".

Dentro de "results" aquí tenemos nuestro archivo de salida. Eso es lo que le dijimos a Nextflow que hiciera. Dijimos, guarda los archivos de resultados en un directorio llamado "results" y cópialos allí. Y así esto ahora es mucho más fácil de encontrar. Está justo allí junto a donde lanzamos un workflow y todos los diferentes archivos se pueden organizar allí como deseemos, independientemente de dónde o cómo Nextflow ejecutó la ejecución real.

Note que publishDir puede manejar enlaces simbólicos, lo cual es bueno si está trabajando en un sistema de archivos compartido y desea ahorrar espacio. Y también no tiene que definir todos los archivos que son creados por un process como una salida.

Nextflow solo copiará las cosas que están definidas en este bloque output. Así que si tiene archivos intermedios creados por el paso, que no son necesarios más adelante de este process, simplemente no los define en output y no aparecerán en el publishDir. Entonces esta es una forma de mantener sus archivos de salida de un pipeline limpios y eliminar fácilmente archivos intermedios una vez que el lugar de trabajo haya terminado.

Una nota rápida aquí. Hay una nueva sintaxis de Nextflow que viene llamada workflow output definitions, que eventualmente reemplazará publishDir. Esto nos da una forma de definir todas las salidas de un workflow a nivel de pipeline en el bloque workflow. Esto se describe en los documentos de Nextflow si quieren probarlo. Pero por ahora, publishDir estará presente por un tiempo, así que todavía lo tenemos en un entrenamiento para 2025.

## 3.2. Relanzar un workflow con -resume

Bien. Mencioné que el directorio work aquí ahora tiene dos conjuntos de resultados con un hash diferente de cada vez que ejecutamos el workflow. Eso es bueno. Sin embargo, a veces no queremos volver a calcular pasos cada vez si no lo necesitamos.

Tal vez estén construyendo iterativamente su workflow y estén agregando pasos y quieren que los primeros pasos simplemente reutilicen las versiones en caché. O tal vez algo salió mal en su sistema de cómputo a mitad de camino a través de su workflow y quieren que continúe desde donde se quedó, pero omita los pasos que ya había completado.

Nextflow tiene funcionalidad incorporada para esto llamada resume. Probémoslo. Entonces primero, solo voy a echar un vistazo al directorio work para que podamos recordar qué había allí.

Y luego voy a hacer "nextflow run hello-world.nf" y voy a agregar un único comando aquí, "-resume".

Note, guión simple, eso es realmente importante. Voy a ejecutarlo y la salida va a verse básicamente exactamente igual, con un par de pequeñas diferencias.

Note aquí dice "cached" en gris. Eso significa que Nextflow no ejecutó la tarea. Esta vez encontró algo que coincidía con lo que eran los requisitos y reutilizó esas salidas directamente en lugar de volver a ejecutar el paso.

Y efectivamente, si miran el hash aquí, pueden ver que esto corresponde al hash existente que teníamos de una ejecución anterior.

## 3.3. Eliminar directorios work antiguos

Bien. Pero si están desarrollando iterativamente, van a acumular muchos de estos archivos de workflow. Eso puede ser un problema si pueden tener poco espacio.

Nextflow puede ayudarnos a limpiar estos directorios work con un par de comandos de ayuda. Si hago "nextflow log". Eso me dará una lista de todas las diferentes ejecuciones de workflow que he hecho en este directorio, y tienen los nombres de ejecución aquí. Pueden ver el gloomy quick que fue el primero que ejecutamos, y luego estos dos nuevos.

Ahora podemos tomar ese nombre y usarlos con el comando "nextflow clean". Puedo especificar un único nombre de ejecución. O incluso mejor, puedo decirle a Nextflow que elimine todo antes de un único nombre de workflow con "-before", y voy a poner "stupefied_shaw". Esa fue mi ejecución más reciente, "-n".

El comando "-n" le dijo a Nextflow que lo hiciera como una ejecución en seco sin realmente eliminar nada de verdad, y nos dice cuáles de los directorios hash habrían sido eliminados. Efectivamente, es solo ese del primer ejecución. Ambas ejecuciones segundas usan el mismo directorio hash.

Voy a ejecutarlo de nuevo, pero ahora en lugar de "-n" para ejecución en seco, voy a hacer "-f" para forzar y ha eliminado ese directorio hash. Ahora si hago "tree work", podemos ver, solo tenemos este archivo de salida restante.

Genial. Así que hemos logrado limpiar un montón de espacio en disco allí.

Un par de cosas a tener en cuenta al eliminar directorios work, si enlazan simbólicamente cosas a su directorio de resultados, esas fuentes de enlace simbólico ahora serán eliminadas y sus resultados se habrán ido para siempre. Por eso usar el modo copy es una cosa más segura de hacer, y generalmente lo que recomendamos.

En segundo lugar, la funcionalidad resume de Nextflow depende de estos directorios work. Entonces, si los eliminan y ejecutan Nextflow nuevamente, la funcionalidad resume ya no funcionará. Así que depende de ustedes hacer un seguimiento de qué cosas pueden necesitar o pueden no necesitar, y solo eliminar cosas cuando estén seguros de que es seguro hacerlo.

La otra cosa que podemos hacer es simplemente eliminar todo el directorio work si hemos terminado nuestra ejecución del workflow y estamos seguros de que ya no lo necesitamos.

Entonces puedo hacer "rm -r work". Sé que no había nada importante allí. Tengo mis resultados que me importan en el directorio results donde los copiamos. Y así fue seguro eliminar el directorio work. Depende de ustedes cuál de estos enfoques usan.

## 4. Usar una entrada variable pasada en la línea de comandos

Bien, ¿qué sigue? Mencioné que habíamos codificado algunos de los valores en nuestro script de workflow aquí, el archivo output.txt, y que podría haber una mejor manera de hacer eso.

Hagamos un comienzo en esto. Lo que vamos a hacer son tres cosas. Vamos a agregar una nueva entrada al process. Vamos a decirle al script del process cómo usar esa entrada, y luego vamos a conectarlo en el workflow para que podamos usarlo dinámicamente con una bandera de línea de comandos al ejecutar Nextflow.

Entonces primero que nada. Agreguemos un bloque input aquí. Igual que output. Esta es una nueva sección para el process, y voy a decir, "val greeting".

Note aquí, estoy diciendo "val", lo que dice que esto es una variable, no un path.

Luego puedo bajar al script y luego puedo sacar este texto codificado aquí y hacer $greeting. Esto funciona como cualquier otro lenguaje de programación. Estamos definiendo una variable aquí y la estamos referenciando dentro de este bloque script. Cuando Nextflow ejecute este process, la variable será interpolada. Y cuando vayamos y miremos ese archivo .command.sh, veremos la cadena real codificada aquí en su lugar.

## 4.1.3. Configurar un parámetro CLI y proporcionarlo como entrada a la llamada del process

Bien, pero ¿dónde proporcionamos la variable? A continuación bajamos a la sección workflow, y pueden ver que la extensión aquí está diciendo, ahora esperamos una entrada, y me ha dado una advertencia.

Ahora, lo más simple que podríamos hacer es simplemente codificarlo. Podría escribir "Hello World" y proporcionar esa entrada de cadena al process. Pero nuevamente, eso realmente no resolvería ningún problema. Todavía tendríamos que volver y editar el código del pipeline cada vez que quisiéramos cambiar algo, lo cual no es bueno.

La buena noticia es que Nextflow tiene un sistema incorporado para manejar argumentos de línea de comandos llamados parámetros. Entonces en su lugar, puedo usar una de estas variables especiales llamada params y puedo llamarla como quiera, pero voy a decir greeting para que coincida con la lógica del workflow.

Guardar y veamos qué podemos hacer con esto.

Entonces si vuelvo a la terminal. Entonces hacemos "nextflow run hello-world.nf". Justo como antes, pero la diferencia clave es que hacemos --greeting

Note, hay dos guiones aquí porque este es un parámetro. Cuando reanudamos el workflow antes, eso fue un guión simple. Eso es porque resume es una opción central de Nextflow, y este es un parámetro que es específico para nuestro pipeline.

No mezclen los dos. Es fácil hacer eso. Si hicieran --resume en lugar de solo un guión, entonces eso sería "params.resume", lo que no haría nada. Del mismo modo, si hicieran un guión simple aquí, Nextflow no lo reconocería como un argumento clave.

Entonces es --greeting, que corresponde a parameters greeting.

Ahora puedo seguir eso con cualquier texto que quiera. Entonces estoy en Suecia en este momento, así que voy a decir, "Hej världen".

Entonces ejecutémoslo, veamos qué sucede, momento de la verdad.

Bien, entonces pueden ver que el process se ejecutó nuevamente, justo como antes, sayHello con una única ejecución.

Esto habrá sobrescrito el archivo que estaba en el directorio publishDir "results". Y así tengan cuidado cuando estén volviendo a ejecutar los archivos porque las cosas en el directorio publicado serán sobrescritas.

Ahora puedo hacer "code results/output.txt", y efectivamente, nuestra salida ha sido actualizada y ahora dice "Hej världen".

## 4.2. Usar valores predeterminados para parámetros de línea de comandos

Bien, eso es genial. Pero el problema ahora es que nuestro workflow depende de que siempre definamos este parámetro, y es bueno tener valores predeterminados sensatos para que las cosas se ejecuten de una manera sensata para su workflow a menos que sobrescriban los valores predeterminados.

Entonces la forma en que hacemos eso es estableciendo un valor predeterminado para el parámetro en nuestro script de workflow.

Entonces si vuelvo a mi archivo hello-world.nf, puedo ir al script justo encima del workflow, escribir "params.greeting" y definirlo como cualquier otra variable. Entonces pongamos una cadena aquí y digamos "¡Hola mundo!"

Ahora este parámetro tiene un valor predeterminado definido, que se usará aquí, o todavía podemos sobrescribirlo en la línea de comandos con --greeting, justo como lo hicimos antes.

Entonces verifiquemos que funcione. "nextflow run hello-world.nf"

Sin argumentos de línea de comandos esta vez, y verifiquemos si hizo lo correcto.

"code results/output.txt". Y ahí está. Obtuvimos nuestro valor predeterminado.

Bien, intentemos de nuevo, solo para verificar que no les esté diciendo ninguna mentira. Ejecutémoslo de nuevo, pero hagamos --greeting, y usemos el ejemplo del material de entrenamiento, digamos "¡Konnichiwa!"

Vuelve a ejecutar, el workflow, y efectivamente, nuestro archivo de salida en la parte superior acaba de actualizarse con el nuevo valor que proporcionamos en la línea de comandos.

Genial. Este es un aspecto realmente central para escribir cualquier workflow de Nextflow. Definir valores predeterminados sensatos en su código de pipeline, pero hacerlo muy fácil de configurar para el usuario final teniendo argumentos de línea de comandos en la terminal.

Note que el usuario final puede sobrescribir la configuración en múltiples lugares diferentes. Pueden tener un archivo de configuración en su directorio principal, que se aplica a cada ejecución de Nextflow que hagan. Pueden tener un archivo de configuración en un directorio de lanzamiento. Pueden tener un archivo de configuración en un directorio de pipeline. Todas estas diferentes ubicaciones de configuración se cargan en un orden específico, que se describe en los documentos de Nextflow.

Bien, ese es el final de la sección uno. Hemos tenido nuestro primer script de workflow en Nextflow con un process y un workflow. Hemos visto entradas, salidas, scripts y publicación, y cómo conectar parámetros y un canal de entrada a nuestro process.

Felicitaciones, su primer paso hacia escribir código de Nextflow está completo.

Tomen un pequeño descanso y los veré de vuelta en unos minutos para el capítulo dos.

[Siguiente transcripción de video :octicons-arrow-right-24:](02_hello_channels.md)
