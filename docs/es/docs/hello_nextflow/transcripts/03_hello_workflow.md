# Parte 3: Hello Workflow - Transcripción del Video

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Traducción asistida por IA - [más información y sugerencias](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

<div class="video-wrapper">
  <iframe width="560" height="315" src="https://www.youtube.com/embed/_aO56V3iXGI?si=Irl9nAQniDyICp2b&amp;list=PLPZ8WHdZGxmWKozQuzr27jyMGqp9kElVK&amp;cc_load_policy=1" title="YouTube video player" frameborder="0" allow="accelerometer; autoplay; clipboard-write; encrypted-media; gyroscope; picture-in-picture; web-share" referrerpolicy="strict-origin-when-cross-origin" allowfullscreen></iframe>
</div>

!!!note "Notas importantes"

    Esta página muestra únicamente la transcripción. Para instrucciones completas paso a paso, regrese al [material del curso](../03_hello_workflow.md).

    Los números de sección mostrados en la transcripción se proporcionan solo con fines indicativos y pueden no incluir todos los números de sección en los materiales.

## Bienvenida y recapitulación

Hola, y bienvenido/a de nuevo a la parte tres de Hello Nextflow. Esta parte se llama Hello Workflow, y es en esta parte del curso donde realmente empezamos a justificar el nombre pipeline o workflow.

Vamos a tomar nuestro script de pipeline simple hasta ahora con su único proceso, y vamos a empezar a agregar procesos adicionales y ver cómo Nextflow maneja esta orquestación y el flujo de datos a través del pipeline.

Regresemos a nuestros code spaces. Verán que he eliminado todos mis directorios .nextflow\* y los directorios work y todo para intentar mantenerlo limpio. No se preocupen si todavía tienen esos archivos de partes anteriores del curso.

Vamos a trabajar desde un archivo llamado hello-workflow.nf. Como antes, esto básicamente representa el script que hemos construido hasta este punto, y nos da un punto de partida limpio. Y nuevamente, abajo en la salida podemos ver que la ruta ahora es hello_workflow. Así que los archivos publicados deberían ir a un subdirectorio diferente en su carpeta results.

Para recapitular dónde estamos hasta ahora, tenemos un solo proceso aquí, con una entrada greeting, una salida greeting file. Y luego el script Bash simple, que solo hace un comando echo a un archivo.

Tenemos una sola entrada de workflow, el bloque params aquí, donde decimos que espera una ruta, y el valor predeterminado es data/greetings.csv, que es este archivo aquí arriba.

Luego en el workflow mismo, tenemos un bloque main. Estamos creando un canal. Estamos analizando el CSV en filas y luego tomando el primer elemento de cada array, y estamos pasando ese canal a ese proceso, que luego está generando tres tareas, y estamos publicando desde el workflow, las salidas de ese proceso.

Y luego finalmente, en el bloque output, le estamos diciendo a Nextflow que publique estos archivos desde este canal al directorio llamado hello_workflow. Y que copie esos archivos en lugar de crear enlaces simbólicos.

## 1. Agregar un segundo paso al workflow

Bien, en esta parte vamos a agregar un segundo proceso a nuestro workflow. Tomaremos las salidas del proceso sayHello, y las procesaremos en un segundo paso, que va a convertir todas las letras dentro de esos archivos convertToUppercase.

Este es solo un ejemplo tonto, es solo un procesamiento de cadenas simple nuevamente, pero les muestra cómo podemos tomar la lógica, dentro del workflow.

Vamos a hacer un comando bash llamado "tr" para esto, que es la abreviatura de translate. Es un comando Unix que ha existido desde siempre. Si no están familiarizados con él, no los culpo. No creo que lo haya usado nunca antes de la capacitación, pero pueden probarlo muy rápidamente en la terminal. Si hago "echo 'hello world'" y luego pipe a 'tr' y luego entre comillas dicen rango de caracteres, así que A a Z, minúsculas, y luego quieren hacer A a Z mayúsculas. Y solo dice, traduce estas letras a estas letras.

Y cuando presiono enter, pueden ver que ahora ha puesto todo en mayúsculas. Muy bueno si les gusta gritarle a la gente.

Así que ese es un estilo muy simple de comando bash que vamos a usar en nuestro segundo proceso.

## 1.2. Escribir el paso de conversión a mayúsculas como un proceso de Nextflow

Entonces, si regreso a mi script, voy a hacer un poco de trampa y simplemente copiar el código de, de los documentos de la capacitación. Pero pueden ver exactamente lo que está pasando.

Tenemos un nuevo proceso aquí. Este lo hemos llamado convertToUpper, pero podríamos llamarlo como queramos.

Tenemos una sola entrada path, como lo hicimos antes. No es un canal de valor, es un canal de ruta. Y luego una sola salida.

En el bloque script hacemos "cat" en el archivo de entrada. Y podemos poner esto entre llaves si queremos. y que toma esa variable. Y ejecutamos ese mismo comando bash en el pipe y escribimos los resultados a un archivo con este nombre de archivo, y eso es capturado por la ruta de salida.

Ahora necesitamos hacer algo con este nuevo proceso. Así que vamos a ir al workflow donde construimos la lógica diferente de un workflow, y después de ese primer proceso, vamos a ejecutar nuestro segundo proceso. Así que convertToUpper es el nombre del proceso aquí.

Toma una entrada así que no podemos simplemente llamarlo por sí mismo. Queremos procesar la salida del primer proceso. Así que al igual que hicimos con esto, sayHello out donde estamos publicando esos resultados. Queremos usar esos mismos resultados aquí como la entrada, así que podemos copiarlos y ponerlos ahí.

Queremos el proceso sayHello ".out", y Nextflow sabe que esto significa un registro de salida simple único aquí, que es este archivo. Así que eso luego será pasado como una entrada a un segundo proceso.

## 1.5. Configurar la publicación de salida del workflow

Bien. Y finalmente, para que realmente guardemos los resultados de este segundo proceso, también necesitamos publicarlos desde el workflow, y luego definirlos en el bloque output, la misma sintaxis que antes. Así que podemos copiar esto y decir second outputs, o como quieran llamarlo.

Tomar el nombre del proceso que nos interesa, convertToUpper out, y luego aquí abajo en el bloque output. Agregar esto y podríamos hacer los mismos atributos aquí. Así que también queremos estos archivos en el subdirectorio Hello Workflow, y también queremos copiarlos.

Genial. Intentemos ejecutarlo. Así que si abro la terminal y hago "nextflow run hello-workflow.nf", y veremos qué hace. A ver si se ve diferente a las partes anteriores.

Así que lanza Nextflow. En los documentos, dice hacer esto con "-resume", pero eliminé todo mi directorio work, así que no habría hecho ninguna diferencia aquí. Pero si lo hicieron, entonces eso también funcionará.

Y se ve casi exactamente igual. Pero pueden ver ahora que hay una segunda línea de salida aquí, donde pueden ver el nombre del segundo proceso que acabamos de agregar. Y efectivamente, pueden ver que se ejecutó tres veces exitosamente.

Brillante. Si tuviera mis directorios work anteriores y hubiera hecho esto con "-resume", estos habrían sido, almacenados en caché solo el primer paso en el pipeline. Porque esas salidas eran exactamente las mismas, así que Nextflow habría sabido reutilizarlas nuevamente.

Y así pueden ver cómo pueden usar -resume para construir iterativamente su workflow, paso a paso, si lo necesitan.

Bien, echemos un vistazo al directorio results aquí arriba y veamos si funcionó. Podemos ver que tenemos algunos archivos más aquí arriba. Tenemos nuestros archivos originales como antes del primer proceso. Y efectivamente, tenemos nuestros archivos upper y las letras están todas en mayúsculas, así que funcionó. Es realmente agradable de ver.

También es interesante solo revisar dentro de estos directorios work. Como antes el, hash aquí corresponde a los, directorios work. Así que si miro en "ls work", y luego expando eso, veremos los diferentes archivos aquí.

Vemos el archivo de salida del primer proceso, que ha sido traído aquí como la entrada. Y podemos ver el nuevo archivo de salida que fue generado.

Ahora si hago esto con "-la" para listar y mostrar todos los archivos, veremos algunas cosas más. Primero, verán que este archivo es en realidad un enlace simbólico al primer proceso. Esto es básicamente siempre un enlace simbólico si puede serlo, para ahorrar espacio de archivo. No estamos publicando los archivos aquí y solo hace referencia a ese archivo de una primera tarea a una segunda tarea para que todo esté encapsulado dentro de ese directorio de trabajo, y seguro y aislado de todo lo demás.

Y eso necesita estar ahí porque si miramos el archivo .command.sh, así que si hago "cat work/b8/56\*", pueden ver que las, partes de archivo aquí son relativas, así que está haciendo cat de ese archivo de entrada, que ha sido enlazado simbólicamente al mismo directorio de trabajo.

Así que así es como se verá cada directorio work. Cuando lo miren en Nextflow, tendrán todos los archivos de entrada ahí preparados en ese directorio de trabajo. Y luego también tendrán cualquier archivo de salida que fue creado. Así que eso es genial. Eso se ve como esperábamos.

## 2.1. Definir el comando de recolección y probarlo en la terminal

Bien, regresemos a nuestro workflow. ¿Cuál es el siguiente paso que queremos hacer?

Tenemos dos procesos ahora y están tomando este archivo CSV, analizándolo y dividiéndolo. Y luego tenemos tres tareas para cada uno de estos procesos y Nextflow maneja la paralelización de todo eso, así que todo se ejecuta lado a lado donde sea posible.

Esa forma de dividir el trabajo para ejecutar cosas en paralelo es muy común. Y lo inverso de eso es luego reunir todo de nuevo. Así que eso es lo que vamos a hacer con nuestro proceso final en el workflow es que tendremos un tercero aquí, que toma estas tres salidas diferentes y las combina todas en un solo archivo.

Podemos hacer esto bastante simplemente en una terminal, solo para tener una idea de cómo se verá esto.

Si voy a la carpeta results. Así que, "cd results/hello_workflow/", y tenemos todos los archivos UPPER aquí. Puedo simplemente usar "cat", que usamos para imprimir el contenido de ese archivo, y pueden dar múltiples archivos a "cat" y leerá uno tras otro.

Así que puedo decir "UPPER-\*", que me da la misma lista de tres nombres de archivo con expansión Bash. Y puedo decir combined.txt. Creo que en los documentos, lista los nombres exactos de archivo, pero está haciendo lo mismo.

Ahora, si uso "cat combined.txt", podemos ver que tenemos el contenido de archivo de los tres archivos.

Así que eso es básicamente todo lo que este proceso va a hacer es que vamos a intentar darle todos los diferentes archivos de salida de un proceso anterior en una sola tarea de proceso, y luego vamos a hacer "cat" de ellos juntos y guardar el archivo de salida.

## 2.2. Crear un nuevo proceso para hacer el paso de recolección

Bien, así que agreguemos nuestro nuevo proceso. Voy a pegar esto de los materiales de capacitación, y pueden ver que nos ha dejado un poco de ejercicio para el lector aquí con estos signos de interrogación. Pero pueden ver el esquema general del proceso es básicamente lo que acabamos de hacer en la terminal, donde estamos haciendo "cat" de un montón de archivos de entrada y escribiéndolo a un archivo de salida aquí llamado collected, y luego la salida espera esa ruta única nuevamente.

Así que necesitamos algún tipo de entrada aquí y van a ser un conjunto de rutas. Así que nuevamente, definimos un canal de entrada path y llamémoslo input_files. Ahora, esto anteriormente nos ha dado una sola ruta aquí, pero una ruta también puede tener múltiples archivos aquí, aunque todavía es una sola declaración.

Voy a copiar eso aquí abajo porque queremos hacer "cat" de estos archivos. Y podrían pensar que tenemos algunos problemas aquí con imprimir un array o cosas así, pero Nextflow es generalmente bastante sensato cuando se trata de esto. Y si se le da un canal con múltiples archivos en él como este, los, pondrá todos juntos con separadores de espacio. Así que esto nos dará la sintaxis correcta.

Eso es genial. Así que ahora conectemos nuestro nuevo proceso. Voy al workflow. Voy a decir combine the outputs, el nuevo nombre del proceso, y justo igual que antes. Voy a tomar este proceso anterior, convertToUpper y hacer ".out".

Genial. Probémoslo y veamos si funciona en la terminal. Si solo regreso un par de directorios y luego vuelvo a ejecutar el comando Nextflow, y veremos qué pasa.

Así que el workflow se ha lanzado y ahora pueden ver que tenemos tres nombres de proceso diferentes, lo cual es genial. Los primeros dos se ven iguales que antes, y el tercero nuevo se ejecuta, lo cual es bueno.

Sin embargo, hay algo un poco extraño aquí. Queríamos combinar esos archivos de salida en un solo archivo, y sin embargo este proceso podemos ver que se ha ejecutado tres veces, no una.

Efectivamente, si vamos a uno de estos directorios work. Y hacemos "cat work/" "collected", entonces veremos. Solo hay una sola palabra aquí, no tres.

Y entonces lo que ha pasado es que Nextflow ha continuado esa paralelización justo como lo hizo en los pasos anteriores. Y este proceso nos dio un canal con tres elementos, y esos tres elementos de canal fueron pasados a nuestro proceso descendente, que generó tres tareas de proceso.

Básicamente intentó recolectar tres veces separadas y cada vez solo tenía un solo archivo, así que solo hizo cat de un solo archivo a una salida, y de hecho, podemos ver eso en el archivo .command.sh también.

Si hago .command.sh, podemos ver que solo tiene un solo nombre de archivo aquí y solo un solo archivo fue preparado en ese directorio de trabajo.

## 2.3. Agregar el paso de recolección al workflow

Así que de alguna manera necesitamos decirle a Nextflow que reúna todas esas salidas de un proceso anterior y se las dé a este proceso descendente como un solo elemento de canal, en lugar de tres.

Hacemos eso con un operador de canal llamado _collect_.

Este es un operador súper útil, que verán en pipelines de Nextflow todo el tiempo. Este es un canal aquí, este canal de salida, justo igual que el que creamos arriba. Y así podemos agregar operadores de canal a él justo como lo hicimos antes. Podemos simplemente hacer punto, y luego en este caso, collect, paréntesis.

Y eso es todo lo que necesitamos. Eso va a manipular este canal antes de que sea pasado a este proceso.

Si quieren ver qué le está pasando, también podemos verlo aquí. Así que aquí, esto no está relacionado con ejecutar este proceso en absoluto, así que podría ponerlo en cualquier punto después de ejecutar ese proceso. Pero tomamos el mismo, canal de salida, y lo estamos mirando con .view, y luego lo estamos mirando nuevamente con .collect.view.

Y cuando ejecutemos esto, nos mostrará las dos estructuras diferentes de ese canal, antes y después de collect. Así que probemos eso ahora. Bien, acabo de alejar un poco porque algunas de las salidas son bastante largas, pero si ejecuto el pipeline, veremos si funciona.

Espero que un tercer proceso se ejecute solo una vez, porque está recolectando las salidas y efectivamente, pueden ver collectGreetings como uno de uno. Así que eso ejecutó solo una tarea.

Y luego si miramos las declaraciones view, tenemos tres declaraciones view para los tres elementos de antes, con una ruta de archivo en cada una.

Y luego después de esa declaración collect, eso solo se activó una vez porque hay un solo elemento en ese canal. Y ahora tenemos esta, lista de tres rutas de archivo diferentes.

Eso es exactamente lo que esperábamos. Y pueden ver con suerte, esto es básicamente lo inverso de ese operador "map" que hicimos para ir de los arrays CSV a elementos de canal separados. Ahora estamos tomando elementos de canal separados y poniéndolos de vuelta en un solo array.

Genial, podemos limpiar estas declaraciones view. Ya no las necesitamos. Podemos pasar al siguiente paso.

Antes de ir más lejos, y antes de que olvide, voy a agregar una nueva declaración publish aquí. Third output. Pueden llamar a esto algo más semántico y descriptivo en su workflow. Y luego voy a agregar eso al bloque output nuevamente y decir path 'hello_workflow' mode 'copy'. Solo para que el archivo de salida generado por este proceso se guarde en nuestra carpeta results aquí arriba.

Solo para verificar rápidamente que funciona. Debería ser un poco más limpio ahora porque no tenemos esas declaraciones view. Y, veremos si obtenemos nuestro nuevo archivo de salida aquí arriba. Una de, una tarea se ejecutó, obtuvimos un nuevo archivo llamado collected, y ahora tenemos las tres palabras. Fantástico. ¿Qué sigue?

## 3. Pasar parámetros adicionales a un proceso

Bien. A continuación vamos a ver cómo manejar múltiples entradas en un solo proceso. Hasta ahora pueden ver que todos nuestros procesos solo están tomando una cosa como entrada. Todos tienen una sola línea bajo su entrada.

Vamos a demostrar esto permitiendo que Nextflow especifique un identificador de lote diferente para que tal vez ejecuten este, workflow múltiples veces y puedan darle un ID de lote diferente cada vez.

Simplemente voy a agregar una segunda línea en la entrada aquí para collectGreetings. Y voy a llamarlo "val", porque esto es una cadena. Ahora es un valor, no una ruta, y voy a llamarlo "batch_name".

Luego voy a editar el script aquí abajo para usar esta variable, y voy a intentar ponerlo en el mismo lugar que el material de capacitación. Así que lo pongo en medio de esta ruta de archivo COLLECTED-$\{batch_name\}-output.

Todavía no terminamos. Recuerden que tenemos que decirle a Nextflow cuáles van a ser los nombres de archivo de salida. Así que también tenemos que hacer lo mismo aquí arriba: COLLECTED-$\{batch_name\}-output.txt".

Fantástico. Nextflow ahora está obteniendo una segunda entrada de variable y la está interpolando en el script y la salida.

Una última cosa, ahora tenemos que encontrar dónde se está llamando esto, y tenemos que pasar la segunda entrada al proceso. Esto es como cualquier otra entrada a una función en cualquier otro lenguaje.

Justo como hicimos antes en la capacitación, voy a usar el "params" especial aquí, y vamos a llamarlo "params.batch" para que podamos tener una opción CLI --batch. Y ahora pueden ver que nuestro proceso aquí tiene dos entradas separadas solo separadas por comas, que están siendo pasadas.

Es realmente importante obtener el orden correcto, así que el orden de argumentos aquí para channel y luego el param debe coincidir. El, channel y el batch name ahí. Esto es solo coincidencia posicional.

Bien. Puedo ejecutar este pipeline ahora directamente con --batch, pero primero hagamos lo correcto y definámoslo en la entrada aquí en Params. Así que voy a agregarlo a batch y luego vamos a decir que es una cadena y démosle un valor predeterminado. Así que simplemente llamémoslo batch. ¿Bien? Ahora intentemos ejecutar el workflow.

--batch Trio. Creo que dice en el material de capacitación, pero podríamos usar cualquier cadena que queramos ahí. Y con suerte veremos que ese archivo de salida de resultados aparezca aquí.

Y efectivamente, COLLECTED-trio-output - eso ha funcionado correctamente. Ha renombrado nuestro archivo. Y pueden imaginar ahora que esto es útil porque si ejecuto eso nuevamente con un nombre de lote diferente, como replicate_two, entonces nos va a dar un nombre de lote diferente aquí arriba.

Y y no va a sobrescribir los archivos de salida en este caso. Así que eso es bueno.

## 4. Agregar una salida al paso de recolección

Bien, así que ahora tenemos múltiples entradas a nuestro proceso aquí. Pero ¿qué pasa si queremos crear múltiples salidas? Nuestro ejemplo aquí entonces es que vamos a crear un reporte para este proceso, solo diciendo cuántos archivos fueron recolectados.

Y haremos eso con un comando echo aquí. Así que podemos decir echo. There were, voy a copiar esto del material de capacitación, para que no tengan que verme escribirlo.

There were $\{count_greetings\} greetings in this batch, y guardar eso en un nuevo archivo ahora llamado $\{batch_name\}, así que la misma variable, podemos reutilizarla tantas veces como queramos, report.txt.

## 4.1.1. Contar el número de saludos recolectados

Necesitamos calcular eso de alguna manera. Podríamos hacer esa lógica en el script Bash si quisiéramos, usando lógica Bash. Sin embargo, también podemos simplemente hacer scripting directamente dentro del código Nextflow, siempre y cuando esté dentro del bloque script en el proceso y arriba de la sección entre comillas.

Cualquier cosa aquí no será incluida en el script renderizado final, y solo será ejecutada por Nextflow cuando renderice una tarea.

Así que aquí solo estamos haciendo algo de lógica. Estamos creando una nueva variable llamada count_greetings. Tomamos el canal de archivos de entrada aquí, y estamos llamando .size() en él.

Bien, esa función me va a dar un número aquí en esta variable, y ahora nuestra advertencia se ha ido porque esta variable está siendo definida.

Bien, así que estamos creando ese segundo archivo en el directorio work, pero necesitamos decirle a Nextflow que lo espere como una salida publicada de este proceso. Así que hacemos eso con exactamente la misma sintaxis que hicimos para el primer archivo.

Decimos path porque es, nuevamente, podríamos estar publicando una variable aquí si quisiéramos con "val", pero vamos a decir "path". Y luego el, nombre de archivo esperado. Noten que no está resaltado aquí. Eso es porque usé comillas simples. Tengo que usar comillas dobles.

## 4.1.2. Emitir el archivo de reporte y nombrar salidas

Bien, eso es genial. Y ahora podríamos empezar a acceder a estas salidas aquí abajo justo como hice aquí. Pero ahora es un array de diferentes objetos, así que podría hacer collectGreetings.out[0] para obtener el primero, o uno para obtener el segundo, que es nuestro nuevo reporte.

Pero realmente no me gusta hacer eso mucho porque es bastante fácil equivocarse con el conteo de índices. Y se sientan ahí contando líneas mucho y agregan una nueva salida y de repente todo se rompe. Así que

es mucho más agradable referenciar todo por nombre en su lugar. Y podemos hacer eso con una clave especial aquí llamada "emit".

Así que podemos llamar a esto como queramos. Llamémoslo emit outfile, y emit reports. Si definen estos y pueden hacerlo en uno o muchos, depende de ustedes. Ahora puedo ir aquí abajo y en su lugar puedo ir punto out punto reports y simplemente llamarlo por nombre, lo cual es mucho más fácil de entender su código cuando lo leen, y es más seguro ante cambios en el código.

He, agregado el .out.report aquí, pero en realidad necesito tener dos salidas diferentes siendo publicadas. Así que voy a renombrar como algo más interesante como collected y report y ¿es así como lo llamé? Lo llamé out file, perdón. Así que ese nombre emit aquí outfile y report. porque estamos publicando dos canales de salida diferentes y así necesitamos referenciar ambos en el bloque publish.

Luego también necesitamos definir estos en el bloque output. Así que renombré eso collected, y nuevamente, para reports, un poco verboso aquí, pero es realmente útil cuando entran a leer un nuevo workflow, ver todas las diferentes salidas aquí, todos los diferentes canales listados lado a lado, y hay formas de hacer esto menos verboso, lo cual tocaremos más adelante.

Bien, probémoslo y ejecutemos nuestro workflow y veamos qué pasa.

Con suerte ahora debería ejecutarse básicamente igual que antes. Y vamos a obtener un nuevo archivo de salida aquí arriba llamado replicate_two, report. Y ahí está. Se abrió y dice que hay tres saludos en el lote, que es lo que esperábamos, así que es perfecto.

Si voy al directorio work aquí solo para probarles que fue ejecutado en el, código Nextflow en lugar del script bash, puedo ir a cat work/ command.sh, y verán aquí que solo está haciendo echo de esta cadena directamente. There were three greetings in this batch, y así esa variable fue interpolada por Nextflow. Fue calculada en el bloque script antes de que escribiera el archivo .command.sh. Así que el cálculo de variable resultante está básicamente codificado en esto antes de que sea ejecutado en su entorno de cómputo en este caso.

Y así pueden ver esa separación entre el script. Bloque aquí y cualquier cosa arriba de él. Espero que eso tenga sentido.

## Conclusión y cuestionario

Bien, ese es el final de esta parte de Hello Nextflow. Así que como antes, vayan y revisen el cuestionario. Háganlo en la página web o en el CLI, repasen algunas de las preguntas y solo verifiquen que han entendido algo del material que hemos cubierto. Vean si hay algo ahí que resalte algo que podrían no haber entendido. No demasiadas preguntas. Agradable y fácil de hacer. O pueden hacerlo en la página web aquí abajo también.

Y tomen un pequeño descanso, un pequeño paseo y regresen y únanse a nosotros en la parte cuatro de Hello, Nextflow, donde hablaremos sobre módulos. Muchas gracias.
