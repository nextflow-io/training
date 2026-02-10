# Parte 2: Hello Channels - Transcripción del Video

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Traducción asistida por IA - [más información y sugerencias](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

<div class="video-wrapper">
  <iframe width="560" height="315" src="https://www.youtube.com/embed/yDR66fzAMOg?si=xCItHLiOQWqoqBB9&amp;list=PLPZ8WHdZGxmWKozQuzr27jyMGqp9kElVK&amp;cc_load_policy=1" title="YouTube video player" frameborder="0" allow="accelerometer; autoplay; clipboard-write; encrypted-media; gyroscope; picture-in-picture; web-share" referrerpolicy="strict-origin-when-cross-origin" allowfullscreen></iframe>
</div>

!!!note "Notas importantes"

    Esta página muestra únicamente la transcripción. Para instrucciones completas paso a paso, regrese al [material del curso](../02_hello_channels.md).

    Los números de sección mostrados en la transcripción se proporcionan solo con fines indicativos y pueden no incluir todos los números de sección en los materiales.

## Bienvenida

Hola y bienvenido/a de nuevo a la Parte 2 de Hello Nextflow. Este capítulo se llama Hello Channels.

Los canales son como el pegamento en tu pipeline de Nextflow. Son las piezas que mantienen unidos todos los diferentes procesos, que Nextflow usa para pasar toda la información y orquestar tu workflow.

Hay otra parte de los canales que son los operadores. Estos son básicamente funciones que podemos usar en los canales para modificar el contenido. Vamos a sumergirnos en VS Code y ver dónde estamos.

Estoy muy ampliado en este VS Code, así que para mantener las cosas limpias y ordenadas, he eliminado todos los archivos _.nextflow\*_ y el directorio _work/_ y el _results/_ y todo del Capítulo Uno. Y simplemente estoy empezando de cero aquí. Pero no te preocupes demasiado por eso. Si no quieres, puedes dejar esos archivos ahí. No causarán ningún problema.

Vamos a comenzar trabajando en _hello-channels.nf_ para este capítulo, y si abro esto, debería verse muy similar al archivo en el que estábamos trabajando anteriormente. Puede ser que diferentes partes estén en diferentes partes del script, pero todo debería ser básicamente lo mismo.

Una cosa que es diferente es que la ruta en el bloque output aquí ahora es _hello_channels_ para esta parte, lo que significa que los archivos de resultado se almacenarán en un subdirectorio diferente en tus resultados si todavía lo tienes ahí. Así que debería ser un lugar agradable y limpio para comenzar sin confundirse con las salidas.

Bien, entonces recordemos rápidamente qué hace este script cuando ejecutamos este workflow. Hacemos _"nextflow run hello-channels.nf"_. Podemos hacer _"--input myinput"_, y cuando ejecutamos esto, va a usar este parámetro, params.input, que se pasó como la variable para el proceso sayHello aquí arriba, que va a greeting y se guarda en output.txt. Y podemos ver eso en el archivo de resultados. Genial.

## 1. Proporcionar entradas variables a través de un canal explícitamente

Eso es bueno. Pero es bastante simplista. Tenemos una variable en este parámetro, que va a un proceso que se ejecuta una vez, y realmente no escala. Y no podemos darle muchos archivos diferentes para crear aquí. No podemos darle muchos saludos diferentes. Solo tenemos uno.

En realidad, Nextflow se trata de escalar tu análisis. Así que probablemente quieras que haga más de una cosa. Y hacemos eso con _canales_.

Los canales son un concepto un poco único para muchas personas que están aprendiendo Nextflow. Viene de estos conceptos de programación funcional, y puede tomar un poco de tiempo entenderlo, pero una vez que haces clic, realmente desbloquean el poder de Nextflow y es clave para cómo escribes tus workflows.

## 1.1. Crear un canal de entrada

Comencemos tomando este script y haciendo que use un _canal_ en lugar de solo un _parámetro_.

Vamos al workflow, que es donde está toda nuestra lógica de workflow sobre cómo unir las cosas. Y voy a entrar aquí y voy a crear un nuevo canal.

Crear un nuevo canal.

Y lo voy a llamar "_greeting_ch"_. Esta es la convención de hacer "_\_ch"_ así, solo para que puedas recordar que esta variable es un canal. Pero puedes llamarlo como quieras.

Y luego voy a decir igual, y voy a hacer _"channel.of"._

Channel es como el espacio de nombres para todo lo relacionado con canales. "c" minúscula si has estado usando Nextflow antes. Y el _".of"_ es algo llamado una fábrica de canales, que es básicamente una forma de crear un canal.

Hay muchas fábricas de canales diferentes. Si solo hago "." aquí, puedes ver que VS Code está sugiriendo muchas de ellas, pero _".of"_ es la más simple y solo toma una entrada aquí.

Así que puedo hacer unos paréntesis y voy a decir _"Hello Channels!"_.

Genial. Tengo un canal. Fantástico. Puedo guardar, podría ejecutarlo de nuevo, pero nada interesante va a pasar. VS Code me ha dado una línea de advertencia naranja aquí y me ha dicho que esto está configurado: has creado esto, pero nunca lo has usado realmente para nada. Este canal no está siendo consumido.

Bien, entonces ¿cómo lo usamos? Muy simple. Voy a tomar esto, copiarlo, y voy a eliminar _params.input_ y voy a poner _"greeting_ch"_ aquí en su lugar. Así que vamos a pasar este canal como la entrada a sayHello.

Ten en cuenta que he codificado esta cadena por ahora. Esto es un poco un paso atrás después de nuestro bonito parámetro que usamos al final del último capítulo, pero simplemente mantiene las cosas simples para empezar para que puedas ver la lógica.

Bien, voy a ir a mi terminal y voy a ejecutar el workflow de nuevo. Sin ningún _"--input"_ esta vez, y va a ejecutarse y va a usar ese canal que hemos creado y con suerte deberíamos tener un archivo aquí arriba en _results/hello_channels/_ y ahora dice "Hello Channels!". Fantástico. Así que eso es lo que esperábamos de nuestro canal aquí. Genial.

## 1.4. Usar view() para inspeccionar el contenido del canal

Una cosa más para agregar aquí, solo una introducción rápida a otra función que podemos usar en canales llamada "_.view"_.

Esto es análogo al comando _print_ en Python u otros lenguajes que podrías estar acostumbrado, y simplemente vuelca el contenido de este canal a la terminal cuando lo ejecutamos.

Así que hago "_.view"_, y luego si vuelvo a ejecutar el workflow, debería imprimir en la terminal cuál es el contenido de ese canal, en el momento en que lo creamos.

Efectivamente, puedes ver que está impreso en la terminal aquí. _"Hello Channels!"_.

Ten en cuenta que puedes dividir estas cosas en líneas si quieres, y de hecho, el formateador automático de Nextflow intentará hacer eso por ti. El espacio en blanco no es realmente importante aquí, así que puedes encadenar estas cosas una tras otra.

## 2. Modificar el workflow para ejecutarse con múltiples valores de entrada

Bien, entonces nuestro canal tiene una cosa que es buena, pero es básicamente lo mismo que antes. Así que hagámoslo un poco más complicado. Agreguemos algunas cosas más a nuestro canal.

La fábrica de canales "_.of()"_ puede tomar múltiples elementos, así que escribamos algunos más. Haremos _Hello, Bonjour, Hej_. Y luego podemos ejecutar este workflow de nuevo y veremos qué pasa.

Debería ejecutarse de nuevo. Y ahora hemos impreso. _"Hello", "Bonjour"_ y _"Hej"_ a la terminal con nuestra declaración view. Fantástico.

## 2.1.2. Ejecutar el comando y observar la salida del registro

Podrías pensar que hemos terminado en este punto. Pero en realidad hay un pequeño problema aquí, que nos va a hacer tropezar. Si miramos nuestro archivo de salida aquí. Puedes ver que tiene _"Hello"_ adentro, pero no tiene ninguna de las otras salidas. De hecho, es solo este.

Si ejecutamos este workflow múltiples veces, incluso podríamos ver que a veces tiene _"Bonjour"_, a veces tiene _"Hej"_. Es un poco aleatorio.

Si miramos la terminal, podemos ver que se ejecutó tres veces y podemos ver las diferentes salidas de view. Pero si voy al directorio work, puedo hacer _"cat work"_. Poner este hash y expandir eso y _output.txt_. Puedes ver que este archivo en el directorio work es diferente al directorio results, y este es _"Hej"._ Así que hay algo que no está funcionando bien aquí.

Y la clave es que tenemos tres tareas que se ejecutaron. La salida de Nextflow intenta resumir eso a medida que avanza el procesamiento, para que no se apodere completamente de toda tu terminal, y ese registro ANSI usa códigos de escape ANSI, básicamente ha sobrescrito las otras tareas. Así que solo te muestra la última que resultó ser actualizada.

## 2.1.3. Ejecutar el comando nuevamente con la opción -ansi-log false

Hay algunas cosas que podemos hacer para entender esto un poco mejor. Podemos mirar en el directorio work mismo y puedes ver todos los diferentes directorios work allí, pero eso es un poco confuso porque se mezclará con diferentes ejecuciones de Nextflow.

O podemos decirle a Nextflow que no use los códigos de escape ANSI.

Así que si ejecuto el comando de nuevo, pero esta vez digo _"-ansi-log false"_ para apagarlo, también podría usar las variables de entorno _$NO_COLOR_ o _"$NXF_ANSI_LOG=false"_. Entonces usa el estilo más antiguo de registro de Nextflow sin ninguno de estos códigos de escape. Simplemente imprime directamente a una terminal sin actualizaciones inteligentes.

Y ahora podemos ver los tres procesos que se ejecutaron. Y cada uno de ellos su propio hash de tarea. Y si entramos en estos directorios work, veremos los tres saludos diferentes que especificamos.

Así que eso tiene un poco más de sentido ahora. Con suerte entiendes que Nextflow estaba haciendo esto, solo estaba siendo un poco inteligente con lo que te mostraba en la terminal con esos directorios work.

Sin embargo, esto está arreglado para un problema con los directorios work, pero no ha arreglado un problema con el archivo de salida. Todavía solo tenemos un archivo de salida que dice _"Hello"_.

## 2.2. Asegurar que los nombres de archivo de salida sean únicos

Ahora para entender esto, necesitamos volver a nuestro script de workflow. Estamos generando nuestro canal aquí, lo estamos pasando a nuestro proceso, y si miramos el proceso, estamos escribiendo el saludo en un archivo llamado _"output.txt"_ y pasando ese archivo de salida de vuelta al bloque output aquí abajo, publicándolo.

Sin embargo, cada tres veces que este proceso ejecuta estas tres tareas diferentes. Todas generan un archivo llamado _"output.txt"_, todos esos archivos de salida se publican en el directorio results, y todos se sobrescriben entre sí. Así que cualquier archivo de resultado que obtengas allí es solo el último que se generó, pero aplastó todos los demás. Eso no es realmente lo que queremos.

## 2.2.1. Construir un nombre de archivo de salida dinámico

Hay diferentes formas de manejar esto, pero la más simple por ahora es simplemente crear diferentes nombres de archivo únicos. Así que cada vez que la tarea se ejecute con un saludo diferente, generará un archivo de salida diferente, que ya no chocará cuando se publique. Y entonces obtendremos tres archivos de salida únicos.

Hacemos esto exactamente de la misma manera. Podemos usar esta variable en cualquier lugar dentro del bloque script y podemos usarla múltiples veces.

Así que puedo pegarlo aquí, _"$\{greeting\}\_output.txt"_, y luego también necesito pegarlo aquí arriba porque ya no estamos creando un archivo llamado _output.txt_. Así que si no actualizo esto, Nextflow fallará con un error diciendo que esperaba un archivo, que nunca se generó.

Así que necesito hacer lo mismo allí y necesito usar comillas dobles, no comillas simples, para que esta variable sea entendida.

Bien, probémoslo y veamos si funcionó. Vamos a ejecutar el workflow de nuevo. Con suerte nos mostrará las tres tareas diferentes dentro de los tres directorios work diferentes. Y efectivamente, puedes ver aquí arriba en la carpeta results aquí a la izquierda. Ahora tenemos tres archivos diferentes con tres nombres de archivo diferentes y cada uno con el contenido diferente que esperamos. Así que los archivos ya no se están aplastando entre sí, y todo está ahí como esperamos.

Esta es una configuración un poco trivial por la que hemos pasado aquí, pero subraya algunos de los conceptos clave que necesitas entender sobre cómo funciona la publicación de archivos, y algunas de las cosas en las que podrías caer como trampas. Así que con suerte puedes evitar eso en tus propios workflows.

También vale la pena señalar que lo que hemos hecho aquí es un poco poco práctico en situaciones de la vida real. Hemos tomado algunos datos de entrada y estamos usando esos datos, pero también estamos nombrando el archivo después de esos datos, lo cual generalmente no puedes hacer.

Así que en pipelines de Nextflow más maduros y reales, a menudo pasarás un objeto meta con todos los metadatos asociados con una muestra dada. Luego puedes crear nombres de archivo dinámicos basados en eso, lo cual es mucho más práctico.

Si estás interesado en cómo hacer esto con las mejores prácticas, hay una misión secundaria en _training.nextflow.io_, que trata específicamente sobre metadatos y mapas meta, así que puedes profundizar allí para más detalles.

## 3. Proporcionar múltiples entradas a través de un array

Bien. A continuación vamos a explorar un poco sobre cómo están estructurados los canales y cómo difieren de otros tipos de estructuras de datos en el lenguaje de codificación. Y voy a pensar un poco sobre cómo podría potencialmente usar un array, que podría ser un concepto familiar si vienes de otros lenguajes.

¿Puedo usar un array en un canal? Intentémoslo. Voy a crear un array, y he copiado esto de los documentos, _"greetings_array"_ y _"Hello", "Bonjour"_ y _"Holà"_. Y luego voy a poner eso aquí en lugar de mis cadenas codificadas. Así que voy a decir "channel.of" _"greetings_array"_, pasando este array a un canal. Intentémoslo.

Abrir la terminal, y ejecutar el pipeline.

Bien. Puedes ver que la declaración view aquí imprimió nuestro array como se esperaba, pero luego todo este texto rojo, o no será rojo si todavía tienes _"-ansi-log"_ apagado, pero todo este texto rojo nos está diciendo que algo salió mal.

Ya no tenemos una marca verde aquí. Tenemos una cruz roja, y si solo hago esto un poco más ancho para que sea más fácil de leer, Nextflow nos está diciendo qué salió mal.

Así que desglosemos esto sección por sección. Dice que el error fue causado por, y luego la razón del error, que son archivos de salida faltantes. Así que básicamente ese bloque output dijo que este archivo debería ser creado y no lo fue. A continuación dice que este es el comando que se ejecutó. Así que esto es básicamente el contenido de ese archivo _.command.sh_. Así es como se veía después de que todas esas variables se pusieran.

Y puedes ver aquí que nuestro comando echo en realidad solo se ha ejecutado una vez y ha usado el array completo, pero en una representación de cadena, que no es realmente lo que queríamos.

Y luego el comando salió así, y ese fue el directorio work donde podemos ir y ver los archivos para entender un poco más.

Bien. Entonces lo que pasó fue que Nextflow simplemente pasó este array completo como un solo elemento de canal al proceso, lo que significó que el proceso solo se ejecutó una vez. Tuvo una tarea y no usó los datos en una estructura que esperábamos.

## 3.2. Usar un operador para transformar el contenido del canal

Así que necesitamos hacer algo a este canal primero, antes de que pueda ser usado. Y esto está preparando el escenario para usar operadores, que son funciones especiales que podemos usar en canales para manipular el contenido del canal.

En este caso, vamos a usar algo llamado _flatten_. Que pasamos al final del canal aquí. Así que creamos el canal y luego ejecutamos _flatten_. Y de nuevo, si pasamos el cursor sobre él, nos muestra la documentación para este comando directamente en VS Code, lo cual es muy útil. También puedes encontrar todos estos documentos en el sitio web de Nextflow, la documentación.

Podría simplemente ejecutar este código ahora y ver si funciona, pero también es una buena oportunidad para introducir cómo hacer código dinámico dentro de operadores y dentro del código de Nextflow, que se llaman closures.

Así que voy a agregar de nuevo un comando view aquí antes de que ejecutemos _flatten_. Y aquí este tiene estas llaves onduladas, que es el closure dinámico. Y solo hay algo de código arbitrario dentro de aquí que se ejecutará, dentro del contexto de un operador view.

Aquí, esto está diciendo toma el saludo, que es la entrada del operador view, y eso está aquí. Podría llamar a esto como quisiera, podría llamar a esto _"foo"_ y solo necesito referirme a él como _"foo"_ más tarde. Y luego digo con esto, devuelve esto.

Y luego establezco devolviendo una cadena que dice antes del flatten para una variable. muy simple.

Ahora voy a agregar otro de estos exactamente igual, pero voy a decir después de _flatten_.

Así que lo que esto hace, porque esto se ejecuta en secuencia, vas a ver cómo se ve el canal antes de que ejecutemos _flatten_, y luego de nuevo después de que ejecutemos _flatten_.

Y luego este canal greeting todavía está creado, así que todavía va a ser pasado al proceso. Y con suerte ahora el workflow se ejecutará. Probémoslo.

Genial. Así que lo primero es que el pipeline no falló esta vez. Tuvimos tres procesos que se ejecutaron correctamente y tenemos una pequeña marca de verificación. Y luego podemos ver que nuestras declaraciones view funcionaron.

Tenemos antes de _flatten_, que es ese array que vimos antes del fallo, y luego tenemos tres veces el después de _flatten_ fue llamado donde tenemos _"Hello", "Bonjour",_ y todos esos otros tres elementos separados en el array, que ahora son como esperábamos, tres elementos separados en el canal.

Y puedes ver que el operador _view_ se ejecutó tres veces. Y eso es porque este canal después de _flatten_ ahora tiene tres elementos. Y entonces el operador se llama tres veces.

Muy rápidamente, solo mencionaría que cuando estaba creando fábricas de canales antes, hice _"."_, y luego vimos que había muchas formas diferentes de crear canales, y una de ellas se llama "_fromList"_. Y eso está específicamente diseñado para hacer esta misma operación. Así que podríamos haber hecho simplemente from list greetings away, y eso funcionará. Es una sintaxis ligeramente más limpia y agradable. Pero para los propósitos de esta demostración, queríamos hacerlo un poco más paso a paso para que pudieras ver cómo se está manipulando el canal y cómo diferentes operadores pueden cambiar lo que está en el contenido de un canal.

## 4. Leer valores de entrada desde un archivo CSV

Bien, ¿cómo podemos hacer esto un poco más realista? Probablemente no vas a querer estar creando mucho código en tu pipeline de Nextflow con arrays codificados. Probablemente vas a querer tomar los datos de afuera cuando lances, y esos datos casi seguramente van a estar en archivos.

Así que lo siguiente que vamos a hacer es que vamos a replicar esto, pero en lugar de tomar los datos de un solo parámetro CLI o de una cadena o array codificado, vamos a tomarlo de un archivo.

Así que deshagámonos de nuestro greetings away. Y ahora vamos a cambiar esta fábrica de canales de nuevo. Acabo de decir que había un montón para elegir y hay uno llamado _".fromPath"_. Y voy a decirle que, en este caso, tome _params.input_, que está volviendo a nuestra entrada que estábamos usando antes.

Ahora ese parámetro realmente no está listo para ser usado todavía. Todavía estamos diciendo que es una cadena y está codificado aquí con un valor predeterminado, pero podríamos sobrescribir esa cadena. Ahora queremos que esto sea un archivo en su lugar. Así que el tipo es diferente. Ya no es un _String_. Es un _Path_.

Y luego podemos establecer el valor predeterminado si queremos, de nuevo, a un Path. Y si miro en explorar a la izquierda, puedes ver en este repositorio, en este directorio de trabajo, tengo un directorio llamado data. Tengo un archivo allí llamado _"greetings.csv"._

Así que puedo simplemente establecer el valor predeterminado aquí a _"data/greetings.csv"_. Ahora, cuando ejecute este pipeline de nuevo sin ninguna opción de línea de comandos, usará este valor predeterminado. Sabe que es una ruta, así que sabe que debe manejar eso como una ruta y no una cadena.

Y luego va a pasar eso a una fábrica de canales desde este _params.input_ y crear nuestro canal, que luego va a ser usado en este proceso llamado _sayHello_. Probémoslo.

Bien. Falló. No te preocupes. Esto era esperado. Y si estás siguiendo el material de capacitación, verás que también era esperado allí. Veamos qué está pasando aquí.

Ha intentado ejecutar el pipeline. Ha intentado ejecutar el proceso, y tiene un error bastante similar al que vimos antes.

Aquí dice: intentamos ejecutar _echo_, pero en lugar de hacer eco del contenido de este archivo CSV, solo hizo eco de la ruta. Y puedes ver que es la ruta absoluta completa aquí a este archivo CSV.

Y luego, efectivamente, porque intentó escribir eso en esta ruta realmente complicada, realmente no sabía qué hacer. Y estaba fuera del alcance del directorio work del proceso.

Mencioné al principio que Nextflow encapsula cada tarea ejecutada dentro de un directorio work especial. Y si intentas escribir en datos, que están fuera de ese directorio work, Nextflow te detendrá como precaución de seguridad. Y eso es lo que ha pasado aquí. Intentamos escribir en una ruta absoluta y Nextflow falló y nos impidió.

## 4.2. Usar el operador splitCsv() para analizar el archivo

Bien, echemos un vistazo a este canal y veamos cómo se ve. Podemos hacer _".view",_ y he copiado esto del sitio web. Así que _.view_, y tenemos un closure dinámico aquí y decimos un nombre de variable "_csv"_ como la entrada. Así que ese es el contenido del canal, y decimos antes de splitCsv, y así es como se ve.

Si lo ejecuto de nuevo, todavía fallará, pero nos mostrará qué hay dentro de este canal. No es particularmente emocionante. Es esa variable _path_. Así que puedes ver que es solo una cadena aquí porque se está imprimiendo en una terminal, pero es un objeto _path_, que contiene la información y metadatos sobre este archivo.

No queremos pasar los metadatos del archivo a la entrada. Queremos pasar el contenido de ese archivo. Si miramos el archivo _greetings.csv_, puedes ver aquí que tiene estas diferentes variables aquí. _Hello, Bonjour, Holà_ de nuevo. Y estas son las cosas que realmente queremos estar pasando a nuestro proceso, no solo el archivo en sí como un solo objeto.

Así que necesitamos analizar este archivo CSV. Necesitamos desempaquetarlo, llegar al contenido del archivo CSV, y luego pasar el contenido dentro del canal al proceso.

Como probablemente puedas decir por el mensaje de registro, queremos usar el _splitCsv_, que es otro operador, otro operador de canal. Así que si hago "_dot" "s"_, y luego puedes ver que está auto sugerido. Ups, _splitCsv_ y algunos paréntesis.

Y luego después de _splitCsv_, voy a poner otra declaración _view_ solo para que podamos ver cómo se ve después. Ejecutemos el pipeline y veamos qué tenemos.

Bien. Todavía falló, pero de una manera nueva y emocionante, que es progreso.

Esta vez de nuevo, tenemos algún problema con nuestro script, que ha sido renderizado. Ahora. Ya no tenemos la ruta final, pero tenemos un array de variables, que se parece mucho al error que tuvimos antes cuando estábamos pasando un array como una entrada fija.

Con nuestro registro del operador view, podemos ver que antes de _splitCsv_ era la ruta. Y efectivamente, después de _splitCsv_, tenemos tres salidas diferentes y cada una de esas salidas se parece mucho a cada una de las filas del archivo _greetings.csv_, lo cual tiene sentido.

Así que lo que ha pasado aquí es que Nextflow ha analizado este archivo CSV dándonos tres objetos, un array para cada línea del archivo CSV. Así que luego tres veces hemos pasado un array de variables al canal en lugar de un solo valor de cadena.

Bien, entonces la última vez que tuvimos este problema, usamos _flatten_. Intentemos muy rápidamente. Probar flatten y ver qué pasa.

Puedo llamar a estas variables como sea. Así que voy a llamarlo _myarray_ porque ya no es realmente un CSV. Intentemos ejecutarlo de nuevo y veamos qué pasa con _flatten_.

Así que esta vez vamos a ejecutar, analizamos el CSV en tres objetos de array, y luego lo aplanamos. Y esta vez, pasó. Y el pipeline de Nextflow se ejecutó. Sin embargo, puedes ver que _flatten_ realmente se pone serio y aplana todo. Y entonces obtenemos tres entradas de array independientes para cada fila. Y entonces ejecutó el proceso tres veces cada fila de un CSV. Y ahora tenemos un montón de archivos de resultados, y 123, 456, y todo tipo de cosas, no solo esa primera columna del CSV, que es lo que realmente queríamos.

## 4.3. Usar el operador map() para extraer los saludos

Entonces, ¿cómo llegamos solo a la primera columna? Si flatten es demasiado simplista aquí, necesitamos un operador más complejo donde realmente podamos personalizar y decirle lo que queremos del CSV.

Para hacer eso, vamos a usar _map_. Básicamente _map_ solo dice, ejecuta algún código, alguna función sobre cada elemento que me den y haz algún tipo de transformación en él. Y porque es tan flexible, lo verás aparecer en el código de Nextflow todo el tiempo.

Por sí mismo, no hace nada. Así que no queremos paréntesis regulares, queremos un closure aquí y necesitamos decirle qué hacer. Así que voy a decir _"row"_, porque eso está recibiendo filas del CSV, así que es un nombre de variable lógico. Es la entrada. Y quiero devolver solo el primer elemento de ese array.

Los arrays en Nextflow están basados en cero, así que vamos a decir solo el primer elemento, que es la fila cero. Si quisiéramos la segunda columna, podría ser uno o la tercera columna ser dos, y así sucesivamente. Podemos devolver lo que queramos aquí, pero voy a devolver solo el primer valor.

Y ahora, podemos ejecutar el pipeline de nuevo y ver si hace lo que esperamos.

Efectivamente, después de _splitCsv_ tenemos nuestros arrays, y luego después del _map,_ tenemos nuestras cadenas limpias y agradables, solo _"Hello", "Bonjour"_ y _"Holà"_. Y el pipeline ahora está haciendo lo que queremos. Fantástico.

Así que podemos deshacernos de todos estos comandos view ahora. Ya no los necesitamos.

## Resumen

Terminamos nuestra especie de depuración y este es el código con el que terminamos. Tomando nuestro parámetro CLI llamado _input_, que está clasificado como un _Path_. Nextflow encuentra la ruta, la carga y entiende el archivo CSV. Devuelve todas las diferentes filas. Y luego mapeamos solo el primer elemento de esa fila en el canal que nos da el contenido del canal, que se pasa al proceso.

Y el proceso se ejecuta sobre cada elemento en el canal, que son tres. Y ejecuta el proceso tres veces, dándole tres tareas. Y esos resultados luego se publican desde el workflow, recogidos por la salida del proceso. Publicados desde un workflow y guardados en el bloque output en un subdirectorio llamado _"hello_channels"_.

Bastante genial. Ahora estamos llegando a algo que se parece más a un pipeline de Nextflow de la vida real que podrías ejecutar para algún análisis real.

## Conclusión

Bien. Con suerte ahora estás teniendo una idea de qué son los canales y operadores de Nextflow y cómo los operadores trabajan en los canales y cómo puedes crearlos.

Los canales, como dije al principio de este video, son el pegamento de Nextflow. Y puedes ver aquí que podemos tomar diferentes entradas y manipularlas y tomar esos datos y luego pasarlos a la lógica de workflow descendente.

Y este bloque workflow aquí es realmente donde construyes toda esa paralelización y toda la lógica inteligente, y explicas a Nextflow cómo construir tu DAG de workflow, y cómo orquestar tu pipeline.

Los canales no son el concepto más fácil de entender. Así que toma un descanso, piensa un poco sobre esto, tal vez lee el material de nuevo, y realmente asegúrate de que tienes estos conceptos claros porque esto es clave para tu comprensión de Nextflow y cuanto mejor entiendas los canales y los diferentes operadores de canal y las diferentes fábricas de canales. Más diversión tendrás escribiendo Nextflow y más poderosos serán tus pipelines.

Esto no es lo mismo que la programación regular en Python u otros lenguajes. No estamos usando declaraciones _if_ aquí, esto es programación de flujo funcional usando canales y operadores. Así que es un poco diferente, pero también es súper poderoso.

Ese es el final de este capítulo. Ve y toma un descanso rápido y te veré en el próximo video para la parte tres donde vamos a repasar Hello Workflow, y hablar un poco más sobre los workflows.

Al igual que el capítulo anterior, hay algunas preguntas de cuestionario en la parte inferior de la página web aquí, así que puedes hacer una revisión rápida de estas y asegurarte de que entiendes todas las diferentes partes del material que acabamos de hacer. Y aparte de eso, te veré en el próximo video. Muchas gracias.

Bien.

​
