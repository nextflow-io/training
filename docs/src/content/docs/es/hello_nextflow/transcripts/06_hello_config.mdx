# Parte 6: Hello Config - Transcripción del video

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Traducción asistida por IA - [más información y sugerencias](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

<div class="video-wrapper">
  <iframe width="560" height="315" src="https://www.youtube.com/embed/FcZTiE25TeA?si=tnXTi6mRkITY0zW_&amp;list=PLPZ8WHdZGxmWKozQuzr27jyMGqp9kElVK&amp;cc_load_policy=1" title="YouTube video player" frameborder="0" allow="accelerometer; autoplay; clipboard-write; encrypted-media; gyroscope; picture-in-picture; web-share" referrerpolicy="strict-origin-when-cross-origin" allowfullscreen></iframe>
</div>

!!!note "Notas importantes"

    Esta página muestra solo la transcripción. Para instrucciones completas paso a paso, regrese al [material del curso](../06_hello_config.md).

    Los números de sección mostrados en la transcripción se proporcionan solo con fines indicativos y pueden no incluir todos los números de sección en los materiales.

## Bienvenida

Hola, y bienvenido/a de nuevo a la Parte Seis de Hello Nextflow. Esta sección trata sobre configuraciones, y es la última parte de este curso.

Nextflow es particularmente bueno en dos cosas: reproducibilidad y portabilidad. Las configuraciones es donde vemos brillar realmente la segunda de estas. La capacidad de configurar un pipeline de Nextflow para ejecutarse de diferentes maneras y funcionar en diferentes sistemas, sin tener que editar el código subyacente del pipeline.

Este superpoder permite que los pipelines de Nextflow sean reutilizados por otras personas en diferentes lugares, o en diferentes infraestructuras a las que tú mismo podrías tener acceso.

Significa que puedes desarrollar código de pipeline en tu computadora portátil, enviarlo a la nube, ejecutarlo en tu HPC, y es el mismo código de pipeline y se ejecuta en todas partes.

En esta sección, vamos a repasar algunos temas. Comenzaremos con cómo Nextflow maneja los archivos de configuración, de dónde los carga, y cómo los escribes y cómo los estructuras, y esa separación entre el pipeline en sí y lo que debería ir en un archivo de configuración.

Luego pasaremos a algunos casos de uso comunes como cambiar dónde se almacenan los archivos de salida, y también cómo hacer que el pipeline funcione en diferentes infraestructuras, tanto usando diferentes tipos de empaquetado de software como enviando trabajos a diferentes infraestructuras.

## Jerarquías de archivos de configuración

Bien, comencemos. Cuando se trata de cargar archivos de configuración, Nextflow puede obtenerlos de muchos lugares diferentes, lo cual es algo bueno y también puede ser algo un poco arriesgado porque a veces puede ser un poco difícil saber de dónde está obteniendo un archivo de configuración y en qué orden carga las cosas.

Así que realmente recomiendo que hagas clic en este enlace aquí, que nos lleva a la documentación de Nextflow. Y en esta página de configuración, enumera los lugares clave de donde se carga la configuración, y lo que es importante, el orden en que estas cosas se cargan.

Entonces puedes ver, puedes poner un archivo de configuración en tu directorio home de Nextflow, que típicamente es ".nextflow" en tu directorio home. Y ese archivo siempre será cargado por cada ejecución de Nextflow en tu sistema.

El siguiente lugar a buscar es un archivo en la raíz de tu repositorio o directorio de pipeline llamado "nextflow.config".

Luego después de eso, otro archivo llamado "nextflow.config", pero esta vez en el directorio desde donde estás lanzando Nextflow: el directorio de lanzamiento.

Finalmente, puedes proporcionar rutas de archivos de configuración en la línea de comandos con un argumento "-c", y puedes hacer eso múltiples veces. Y se aplican en el orden en que los especificas.

Puedes proporcionar archivos de configuración en todos estos lugares si quieres, y se cargarán iterativamente, cada uno sobrescribiendo el anterior solo en los ámbitos de configuración donde entren en conflicto.

Este es un sistema realmente poderoso porque significa que puedes establecer valores predeterminados sensatos y luego hacerte gradualmente más y más específico a medida que te concentras en esa configuración.

## 0. Calentamiento: Ejecutar hello-config.nf

Bien, cerremos esto y saltemos a nuestros Codespaces y comencemos. Como antes, he limpiado aquí, he eliminado mis directorios de resultados anteriores, mis directorios de Nextflow y work, y demás. No te preocupes si todavía tienes esos archivos por ahí. Es solo porque estoy muy acercado y las cosas se ensucian muy rápidamente de otra manera.

Vamos a estar trabajando con hello-config.nf, el último archivo en nuestro directorio, y esto debería seguir desde donde lo dejamos en la sección anterior.

Entonces tenemos nuestros cuatro procesos diferentes, que están siendo incluidos desde archivos de módulos. Tenemos nuestros parámetros del pipeline, nuestro bloque workflow donde estamos llamando a los diferentes procesos y uniendo los canales, publicando los canales de salida, y luego el bloque output en la parte inferior donde definimos dónde esos archivos deberían ser almacenados y cómo deberían ser copiados.

También ya tenemos un archivo "nextflow.config" del último capítulo, donde habilitamos Docker, y vamos a estar construyendo sobre este archivo hoy.

Como antes, hemos cambiado la ruta de salida en este script principal a hello config, solo para que no entre en conflicto con resultados anteriores que has generado.

Bien, vamos a verificar rápidamente que todo sigue funcionando como esperamos. Abro un terminal y hago nextflow run hello-config.nf. Nextflow se carga. Debería ejecutar nuestros cuatro procesos diferentes. Generar algo de bonito arte ASCII usando cowpy y luego guardar nuestros resultados en nuestros archivos de resultados en ese directorio.

Puedo echar un vistazo rápido aquí solo para asegurarme de que estos archivos se vean como esperamos, y efectivamente, ahí está nuestro pavo gigante. Genial.

## 1.1. Mover valores predeterminados a nextflow.config

Ahora lo primero que vamos a hacer es comenzar a mover algunas cosas de nuestro script a nuestro archivo de configuración.

Y lo que nos importa son principalmente los parámetros en esta etapa. Queremos llevar los valores predeterminados al archivo de configuración, para que esté más claro cuáles son los valores predeterminados y sea más fácil para las personas sobrescribirlos.

Voy a tomar este bloque params aquí del script y ponerlo en el archivo de configuración. Y necesitamos tener un poco de cuidado aquí, porque ahora mismo la sintaxis es ligeramente diferente entre config y scripts. El archivo de configuración no puede tomar declaraciones de tipo porque realmente no estamos definiendo estos params, solo los estamos referenciando. Así que voy a deshacerme de esos.

Pero por lo demás es muy parecido. Tenemos un bloque params y luego tenemos nuestros diferentes parámetros de entrada, parámetro batch, parámetro character.

Ahora puedo volver a mi script y ya no necesito definir estos valores predeterminados porque estos valores ahora están en mi archivo Nextflow config.

Sin embargo, sí dejo los nombres de los parámetros y sus tipos, para que Nextflow conozca esa información y aún pueda hacer toda la seguridad de tipos y todo eso.

Bien. Guardamos esos archivos y verificamos rápidamente que todo sigue funcionando igual que antes. No debería haber ningún cambio aquí. Hemos mantenido los valores iguales. Solo hemos movido dónde han sido definidos.

Genial.

## 1.2. Usar un archivo de configuración específico de ejecución

Ahora, hasta ahora hemos estado lanzando Nextflow desde el mismo directorio donde tenemos nuestro script del pipeline. Así que nuestro directorio de lanzamiento y nuestro directorio del pipeline son básicamente lo mismo.

Para mostrar cómo podemos tener diferentes archivos de configuración con diferentes directorios de lanzamiento, vamos a crear un nuevo subdirectorio ahora.

Así que voy a decir mkdir, y vamos a llamarlo tux-run.

Y luego voy a hacer cd, cambiar directorio a tux-run. Y nota que ahora estamos en nuestro directorio de trabajo ya no está en el mismo directorio que los scripts del pipeline.

Bien, vamos a crear un nuevo archivo "nextflow.config". Así que touch Nextflow config, y vamos a abrirlo en VS Code. También puedes ver en la barra lateral aquí que ahora estamos en este subdirectorio.

Ahora podemos tomar el mismo bloque params que teníamos en nuestro nextflow.config de nivel superior, copiar esto y ahora podemos cambiar estos valores.

Primero, los data ahora son una ruta relativa diferente porque estamos en un subdirectorio, así que necesitamos actualizar eso. Y luego vamos a cambiar batch a experiment, y vamos a cambiar el character de Turkey a tux.

Ahora hago clic en guardar ahí, y vamos a probarlo. Al igual que con data, ahora necesito decir ../ para llegar al script. Así que es Hello config. Y presiono enter.

El código del pipeline no ha cambiado en absoluto, pero ahora vamos a tener dos conjuntos de configuración cargándose, y el archivo de configuración del directorio de lanzamiento debería sobrescribir los valores predeterminados, que fueron establecidos en el nextflow.config del pipeline, y deberíamos obtener diferentes conjuntos de resultados.

Efectivamente, dentro de nuestro directorio aquí, dentro de tux-run, puedes ver que tenemos un directorio punto Nextflow y un directorio work y eso es porque estos se crean siempre en tu directorio de lanzamiento. Así que estos son diferentes a los work y results que teníamos de ejecuciones anteriores.

Ahora, si miro en results, podemos ver nuestro collected y ahí está nuestro pequeño personaje tux. Así que puedes ver que esos parámetros fueron interpretados correctamente.

## 1.3. Usar un archivo de parámetros

Bien. Antes cuando estaba hablando sobre los diferentes archivos de configuración que podían cargarse, omití otro lugar del que podemos obtener configuración.

Puedes obtenerla desde una línea de comandos como hemos visto con guion guion nombres de parámetros, pero también podemos proporcionar un archivo YAML o un archivo JSON, solo de params.

El archivo de configuración puede tener todos los diferentes tipos de ámbitos, pero estos archivos son solo parámetros, y es una forma agradable y amigable para el usuario de proporcionar muchos parámetros a la vez, y quizás una forma un poco más reproducible porque los escribes en un archivo, por lo que es fácil obtenerlos en una etapa posterior.

Así que volvamos a nuestro terminal y justo antes de que olvidemos, asegurémonos de volver a subir un directorio, así que ya no estoy en el subdirectorio, y voy a mirar el archivo YAML que tenemos aquí llamado test-params.yaml.

Así que si solo hago code test-params.yaml, puedes ver que esto es solo un archivo YAML regular. Nada especial al respecto. Con las claves siendo nuestros nombres de parámetros, con el formato YAML así que dos puntos aquí, y luego un valor.

Nota que esto no es código Nextflow, así que no podemos poner cosas como variables aquí. Estos son solo valores estáticos.

También porque JSON en realidad se analiza como YAML, también podemos tener un archivo test-params.json, que se ve muy similar. Es solo un formato de datos diferente.

Así que tenemos dos archivos de prueba diferentes aquí y tenemos variables ligeramente diferentes.

Bien, entonces ¿cómo damos estos a Nextflow? Es muy simple. Hacemos Nextflow run hello config, como antes. Y en lugar de "-c" para archivo de configuración, o cargar esos nombres de archivo predeterminados, hacemos -params-file. Un solo guion porque es una opción central de Nextflow.

Y luego pasamos la ruta para ese archivo. Así que voy a hacer "-params-file test-params.yaml", y veremos si esos se cargan correctamente.

Bien. Se ejecutó. Vamos a recordarnos qué había en este archivo YAML. Entonces el batch estaba establecido en YAML, así que así es como debería llamarse, y debería tener un stegosaurus. Así que vamos arriba y miremos en results. Y tenemos COLLECTED-yaml. Así que veamos si tenemos un Stegosaurus. Fantástico, un Stegosaurus usando un sombrero. Eso es lo que nos gusta.

Así que eso funcionó realmente bien, y es exactamente lo mismo con el archivo JSON. Solo cambiamos la extensión del archivo aquí y Nextflow sabe cómo leerlo.

Y en este caso, deberíamos tener un batch llamado JSON y deberíamos tener una tortuga. Vamos a echar un vistazo. Maravilloso. Una de mis herramientas CLI favoritas.

## 2.1. Personalizar el directorio de salida con -output-dir

Bien, así que eso ha sido principalmente pensar en las entradas al pipeline y cambiar parámetros. ¿Qué hay de las salidas?

Ahora, aunque hemos estado cambiando los subdirectorios usando params, podrías haber notado que todos nuestros archivos siguen yendo a results.

Podemos cambiar ese directorio base al que se publican todos los archivos con un indicador de línea de comandos llamado -output-dir. Así que si hago Nextflow run hello config, y luego hago -output-dir, y vamos a llamarlo "custom-outdir-cli". No puedo escribir. Solo para que recordemos de dónde vinieron estos archivos.

Esta es una opción central de Nextflow y es muy nueva. Esto fue añadido recientemente, y esta es una de las cosas que podemos hacer con el nuevo analizador de lenguaje y todo.

Es un poco difícil de escribir. También puedes simplemente llamarlo "-o" si quieres. Así que si solo voy atrás. Puedo simplemente acortar eso a "-o", que es un poco más simple.

Bien. Ejecutamos eso. No hemos cambiado nada en nuestro pipeline o incluso en nuestro config en este punto, y debería guardar con suerte todos nuestros resultados en un directorio de nivel superior diferente. Y puedes imaginar que puedes establecer esto a básicamente cualquier ruta que quieras.

Acaba de llegar arriba. Tenemos un custom-outdir-cli, y todos los archivos están organizados ahí de exactamente la misma manera, con sus mismos subdirectorios y nombres de archivo. Así que esta es una manera realmente fácil de simplemente cambiar dónde el pipeline publica sus resultados, sin pensar demasiado en cómo esos resultados están organizados.

## 2.1.2. Eliminar rutas codificadas del bloque output

Si miro en este directorio, podemos ver que todavía tenemos un subdirectorio llamado Hello Config, que se siente un poco redundante ahora.

Así que simplemente carguemos nuestro script de nuevo y ahora podemos eliminar ese subdirectorio del bloque output en la parte inferior. Porque realmente ya no lo necesitamos. Así que podemos hacer eso ahora, eliminar eso de aquí. Y luego si es solo esto, puedes eliminar eso completamente o dejarlo como una cadena vacía. Voy a dejarlo como una cadena vacía por ahora, porque vamos a volver y poner algunas cosas diferentes en su lugar en el futuro. Pero si no te importan los subdirectorios, es más limpio simplemente eliminar completamente la declaración de path ahí.

Bien, vamos a guardar. Solo probemos rápidamente de nuevo. De hecho voy a eliminar mi directorio "custom-outdir-cli" para que no nos confundan los archivos existentes ahí. Porque recuerda, cuando publicas cosas, no elimina los archivos que ya estaban ahí. Solo añade nuevos. Vamos a ejecutar ese comando de nuevo, custom-outdir-cli.

Y ahora si haces "ls custom-outdir-cli", ya no hay más directorio ahí llamado Hello Config.

## 2.2.1. Establecer outputDir en el archivo de configuración

Bien, el indicador de línea de comandos aquí, "-o" o "-output-dir" es bueno. Pero ¿qué hay de establecer valores predeterminados para esto en config? ¿Cómo hacemos eso?

Abro el archivo "nextflow.config", cierro todo lo demás y me deshago de eso. Podemos añadir una nueva opción de configuración aquí, que he copiado del sitio web del material de capacitación, y se llama outputDir.

No está bajo ningún ámbito. No está bajo params ni nada. Es de nivel superior, y podemos establecer esto en una cadena. Ahora una cosa simple de hacer es simplemente cambiarlo a cualquier cosa que no sea results como una cadena codificada. Pero porque esto está en un archivo Nextflow config, podemos ser un poco inteligentes aquí y también incluir variables.

Y puedes ver aquí que hemos incluido una variable params, params.batch, que es parte de esta cadena. Esto significa que podemos reutilizar variables que están viniendo de otros lugares. Y en este caso, si hacemos --batch, cuando ejecutamos Nextflow Pipeline, vamos a obtener un subdirectorio en nuestra ruta personalizada basado en cuál era el nombre del batch.

Bien, así que probemos esto y solo echemos un vistazo rápido para ver cómo se ven los resultados. Así que si hago Nextflow run hello config y --batch my_run. Recordemos cómo se veía el config. Así que es custom-outdir-config.

Tree custom-outdir-config. Y puedes ver que el batch se llamaba my_run. Y luego tenemos ese subdirectorio llamado my_run. Así que esa ruta de archivo dinámica funcionó.

Y no solo eso, ya no fue a un directorio results predeterminado, y no tuve que especificar nada en la línea de comandos para cambiar el directorio base. Así que hemos restablecido exitosamente el valor predeterminado para el outputDir predeterminado.

## 2.2.2. Subdirectorios con nombres de batch y proceso

Bien, llevemos eso un poco más allá. Esa es una variable dinámica dentro del archivo de configuración. ¿Qué hay del script en sí? Ahora, hasta ahora hemos tenido estas rutas aquí y estas también pueden ser dinámicas. Así que en lugar de solo codificar algo, podemos poner algunas llaves onduladas y poner algo dinámico.

Así que por ejemplo, tenemos nuestros procesos llamados sayHello. Podríamos hacer sayHello.name, que es un atributo del proceso, que es un poco aburrido porque es solo "sayHello" en este caso. Pero es variable.

Así que esto te da una idea. Así que podemos poner esto aquí y decir convertToUpper.name, collectGreetings.name, collectGreetings.name de nuevo, y cowpy.

Ahora cuando ejecutemos, el directorio base todavía va a ser custom-outdir-config. Y va a estar en un subdirectorio llamado params.batch, pero los subdirectorios bajo eso deberían estar organizados por nombre de proceso.

Vamos a probar eso y ver si funciona. Así que voy a eliminar el directorio anterior para que no nos confundamos, y solo usar exactamente el mismo comando Nextflow Run.

Debería ejecutarse de la misma manera. Podría estar usando dash resume en todos estos para hacerlo un poco más rápido y usar los resultados previamente calculados. Ahora, si hago tree custom-outdir-config, puedes ver que no está en results, está en nuestro directorio base con el nombre del batch. Y puedes ver que todos los resultados ahora están organizados dentro de subdirectorios nombrados después del proceso. Así que tenemos dos lugares diferentes donde estamos definiendo rutas de salida dinámicas aquí.

Bien. Última cosa, vamos a añadir de vuelta esas carpetas intermedias, que teníamos antes porque eran un poco agradables. Intermediates.

Y también podemos pensar un poco sobre este params.batch, tal vez como desarrollador del pipeline realmente me gustaba tener eso en el subdirectorio, pero si los usuarios finales del pipeline están estableciendo "-o" o -output-dir en el CLI, está sobrescribiendo completamente toda esta declaración, y perdemos ese subdirectorio.

Así que lo que podemos hacer es que podemos sacar esa ruta dinámica del outputDir config, que sería machacado, y ponerlo en la ruta de salida, que no es machacada.

Así que podemos hacer params.batch barra intermediates barra sayHello.name, y hacer todo esto en una cadena entre comillas dobles, para que sea interpolado por Nextflow.

Ahora puedo copiar, ups. Copiar estos hacia abajo a los otros procesos. Recuerda ponerlos todos entre comillas. Y eliminar intermediates de estas salidas particulares.

¿Bien? Se ve ligeramente más complejo ahora, pero puedes ver que realmente estamos empezando a construir una bonita estructura de directorio de salida organizada en nuestro código.

Y lo que es realmente agradable es que esta complejidad extra en el código no pasa a través del CLI. Así que podemos ejecutar nuestro comando con -output-dir y cualquier variable batch, solo pensando en cómo ejecutar el pipeline y sin pensar realmente demasiado en lo que está en el código. Y nuestros archivos de salida van a ser construidos realmente bien de una manera muy bien organizada, lo cual es agradable para las personas que usan el pipeline básicamente.

Genial. Mientras escribo esto, me doy cuenta de que cometí un error. Veamos si alguien me atrapó aquí. Tenemos collectGreetings.name, así que algo ha ido un poco mal. Y sí, efectivamente, accidentalmente olvidé poner estos entre llaves onduladas.

Así que recuerda, ten cuidado cuando estés escribiendo tu código y asegúrate de decirle a Nextflow qué es una variable y qué es solo una cadena. Porque hará exactamente lo que le digas que haga. Y nada más. Como todos los buenos computadores. Bien, eso debería arreglarlo.

## 2.3. Establecer el modo de publicación a nivel de workflow

Hay un poco de este script, que todavía no me encanta, que es el hecho de que estamos escribiendo mode copy una y otra vez, y si hay una cosa que no nos gusta, es repetirnos.

Así que podemos limpiar esto un poco tomando esto y moviéndolo al config. Y de hecho, podemos establecerlo para todo el pipeline de una sola vez. Así que no tenemos que decirlo múltiples veces.

Vamos a nuestro archivo de configuración y tenemos un nuevo ámbito aquí llamado workflow. Y podemos hacer llaves onduladas o podemos hacer notación de punto. No hace ninguna diferencia. El sitio web del material de capacitación usa notación de punto. Puedo decir output y podemos mezclar y combinar, así que mode equals copy. Genial.

Y ahora podemos volver aquí y eliminar estos. Ahora podríamos dejarlos en su lugar. El config está básicamente sobrescribiendo lo que está escrito aquí, pero como lo tenemos en el config a nivel del pipeline, y estos dos archivos se envían juntos, realmente no hay razón para hacerlo dos veces.

Bien. Solo una verificación de cordura, porque aparentemente sí cometemos errores. Vamos a ejecutar eso de nuevo y solo verificar que estamos usando correctamente el modo copy para publicar archivos. Así que vamos a ejecutar el script de nuevo y esta vez hemos puesto los resultados en un directorio llamado config-output-mode, veamos cómo se ven los archivos ahí.

Y luego si hago "ls -l" para mirar batch, y podemos mirar cowpy, por ejemplo. Y deberíamos ver, sí, que este es un archivo propio aquí, que no es un enlace simbólico, así que ese atributo de configuración ha sido aplicado correctamente.

## 3. Seleccionar una tecnología de empaquetado de software

Bien. Hasta ahora nos hemos estado enfocando en las entradas y las salidas, los archivos con los que el workflow está ejecutándose. Pero ¿qué hay de la infraestructura? Dije al principio que Nextflow te permite ejecutar el mismo pipeline en diferentes configuraciones computacionales. Entonces, ¿cómo se ve eso?

Para mostrar esto, vamos a cambiar de usar Docker para ejecutar cowpy, y en su lugar usaremos Conda para hacer lo mismo.

Puedo hacer esto muy simplemente. Si voy a code, "nextflow.config". Si recuerdas en la parte superior, definimos docker.enabled anteriormente, y en el último capítulo para que pudiéramos usar el contenedor con cowpy dentro.

Voy a decirle a Nextflow que no use Docker. Establecer eso en false. Y voy a decir Conda enabled equals true. Así que decirle a Nextflow, por favor usa Conda.

Ahora simplemente habilitar Conda no es suficiente por sí solo. Justo como lo hicimos con Docker, tenemos que decirle a Nextflow dónde puede obtener el software que necesita.

Así que si saltamos a los modules aquí. Y abrimos el script cowpy. Podemos ver que tenemos una declaración container en la parte superior. Y el container es usado por Docker, pero también Singularity, Apptainer, y muchas de las otras herramientas de software.

Pero no puede ser usado para Conda, así que tenemos una declaración separada llamada "conda", y podríamos simplemente escribir "cowpy". Y eso lo dejará a la resolución de paquetes de conda para descubrir la mejor manera de resolver eso, según tu entorno conda local.

O es buena práctica hacer lo que el sitio web del material de capacitación dice hacer, que es definir un canal conda específico con su notación de doble dos puntos, y definitivamente definir una versión específica del software para que cada persona que ejecute el pipeline obtenga la misma versión.

Nota que los contenedores son un poco superiores en este aspecto, porque cuando instalas algo con Conda, todavía va a resolver todas las dependencias para ese paquete, y pueden cambiar con el tiempo. Llamado deriva de dependencias.

Así que los contenedores, sin embargo, bloquean toda la pila de todas las dependencias de software hasta el fondo, así que puedes estar un poco más confiado de que A, va a funcionar, y B, será reproducible.

Así que si puedes usar Docker o Singularity o Apptainer, definitivamente lo recomendaría.

Ahora lo que es agradable de esto es que el archivo de módulo, que está escrito por el desarrollador del pipeline, ahora tiene tanto Container como Conda, y así estamos diciéndole a la persona que está ejecutando este pipeline, no nos importa qué solución de empaquetado de software uses. Funcionará tanto con Docker como con Conda, y aquí es donde obtener el software en ambos casos.

Podemos abrir el terminal y vamos a darle una oportunidad. Así que Nextflow run hello config --batch conda. Y la primera vez que esto se ejecute con conda, va a ser un poco lento cuando llegue a ese proceso particular, porque tiene que ejecutar "conda install".

Y está creando un entorno conda especial solo para este un proceso. Así que no está usando mi entorno conda global, que tengo en mi terminal. Está creando uno solo para ese un proceso. Esto es bueno porque evita cosas como conflictos de dependencias entre diferentes procesos en tu workflow. Si tus procesos tienen herramientas que necesitan diferentes versiones de Python o cosas así, está bien porque están usando diferentes entornos conda.

Nextflow almacena en caché estos entornos conda localmente, puedes ver que te dice dónde está esa ruta, está en el directorio work aquí. Y así la próxima vez que ejecute este script con Conda, será mucho más rápido porque encontrará ese entorno conda existente y simplemente lo reutilizará. Pero la primera vez que lo hacemos, tiene que ir y buscarlo, resolverlo, descargar todas las dependencias, y configurar todo.

Bien, genial, se ejecutó. Podemos simplemente recordarnos qué está configurado actualmente el pipeline para usar. Si miramos en el archivo de configuración, era "custom-outdir-config" ahora mismo para mí. Veamos si voy hasta ese directorio base. E hice --batch conda. Ahí está nuestro subdirectorio conda. Así que funcionó y ahí está nuestra salida de cowpy.

Así que buscó cowpy, lo instaló en mi sistema local usando conda, y ejecutó el proceso. Y lo que es genial es que, como ese usuario final, no tuve que pensar en absoluto sobre ninguna de la gestión de software ahí. Nextflow simplemente lo resolvió por mí. Dije, necesito usar conda en este sistema. El desarrollador del pipeline dijo qué paquetes necesitaba. Y Nextflow hizo el resto. Muy poderoso.

Nota que en realidad puedes usar una mezcla de diferentes tecnologías. Así que puedo habilitar Docker para procesos específicos, y conda para otros procesos, o decir que algunos procesos deberían simplemente usar cualquier software local que tuviera instalado. Esto es bastante inusual, pero es posible, y en algunos casos, por ejemplo, si estás usando cierto software que podría ser difícil de empaquetar en Docker, tienes una escapatoria.

## 4. Seleccionar una plataforma de ejecución

Así que eso es empaquetado de software. La otra parte de la portabilidad a otros sistemas es dónde se ejecutan realmente los trabajos. En este momento, estoy ejecutando básicamente en mi computadora portátil o en estos Codespaces, que es una sola computadora. No hay nada elegante. Nextflow está siendo un poco inteligente sobre paralelizar los trabajos lo mejor que puede, pero todo está en un sistema.

Ahora, si estás ejecutando en un HPC, probablemente tienes algún tipo de programador de trabajos como SLURM o PBS o algo, y enviarás trabajos a ese programador y distribuirá todos los trabajos a diferentes nodos de cómputo.

Otra forma de ejecutar es en la nube. Así que tal vez estés usando AWS Batch, o Azure Cloud, o Google. Y todos estos funcionan en un sistema similar donde tienes un programador y envías trabajos y se envían a diferentes lugares para ser computados.

Ahora en el pasado lejano cuando comencé a hacer bioinformática, el software de todos para ejecutar análisis estaba muy atado a su infraestructura computacional, lo que lo hacía casi imposible de replicar.

Pero con esta separación de configuración en Nextflow, y con la capacidad de Nextflow de interactuar con muchos backends de infraestructura de cómputo diferentes, es muy simple tomar nuestro pipeline sin modificar el código del pipeline en absoluto y simplemente cambiar eso.

## 4.1. Apuntando a un backend diferente

Así que si vamos a nuestro archivo "nextflow.config", y ahora podemos poner algo de configuración a nivel de proceso. Así que si pongo en la parte superior ámbito process y puedo establecer el executor, y aquí está establecido en local, que es el predeterminado.

Nota que porque esto es a nivel de proceso, podemos apuntar cosas a diferentes procesos. Y así puedes en realidad configurar executors para ser específicos de proceso y tener una ejecución híbrida, donde algunos trabajos podrían ejecutarse localmente, donde sea que el trabajo de Nextflow esté siendo ejecutado. Algunos se envían a diferentes HPC y algunos podrían enviarse a la nube. Puedes ser tan inteligente como quieras.

Ahora, es muy difícil demostrar esto en un entorno de capacitación como este porque no tengo un HPC al cual enviar. Pero lo que puedo hacer es si escribo slurm, podemos hacer un poco de trampa y puedes tener una idea de esto.

Y esto es realmente más interesante para personas que están acostumbradas a ejecutar en SLURM y conocen cómo se ven los encabezados SLURM. Pero si hago Nextflow run, hello config. Va a fallar porque va a intentar enviar trabajos a un clúster que no existe. Así que obtendremos algún tipo de error sobre sbatch no estar disponible.

Sí, escrito. Esa es la herramienta. Esa es la herramienta CLI que usas para enviar trabajos a un clúster slurm. Pero lo que podemos hacer es que podemos ir y mirar en nuestro directorio work aquí con command clic, abrir ese directorio y mirar el .command.run. Y puedes ver en la parte superior del archivo .command.run, tenemos nuestros encabezados sbatch, diciéndole a un clúster SLURM teórico cómo manejar este envío de trabajo.

Y así puedes ver que Nextflow está siendo inteligente, está haciendo todas las cosas correctas. Solo que no teníamos un clúster al cual enviar.

## 5. Controlar asignaciones de recursos computacionales

¿Qué más es diferente entre diferentes infraestructuras computacionales? Otra cosa es cuántos recursos disponibles tienes, y de hecho, en muchos entornos de cómputo, es un requisito que tengas que especificar cuántos CPUs y cuánta memoria necesita un trabajo.

De nuevo, Nextflow abstrae esto para nosotros, para que ya no sea específico de un solo tipo de entorno de cómputo, y podemos escribir en el ámbito a nivel de proceso aquí. CPUs equals one, memory equals two gigabytes. Nuestro pipeline no es muy exigente, así que eso debería estar bien.

Ahora, solo he adivinado estos números aquí, pero ¿cómo sabes qué es una cantidad sensata de recursos para usar? Es un trabajo bastante difícil ir y escarbar a través de todos estos diferentes procesos de un gran pipeline de muchas muestras y entender cuál fue la utilización de recursos.

Así que un buen enfoque para esto es establecer estos valores en números altos para empezar, solo para que tu pipeline se ejecute sin ningún error, y luego pedirle a Nextflow que genere un reporte de uso para ti.

Esto es súper fácil de hacer, así que voy a volver a un terminal. Oh, necesito recordar establecer eso de vuelta a local para que mi pipeline en realidad se ejecute. Y voy a decir Nextflow run, y voy a usar un indicador de línea de comandos -with-report.

Y puedo dejar eso en blanco y dará un nombre de archivo predeterminado, pero voy a darle un nombre de archivo específico para que eso se guarde en un lugar específico.

Presiono Enter, y el pipeline se ejecuta exactamente como normal, pero cuando termine, va a generar un bonito reporte HTML para mí.

Así que en la barra lateral aquí, tengo este archivo HTML. Si estuviera ejecutando esto localmente, simplemente lo abriría. Estoy, porque estoy en Codespaces, voy a hacer clic derecho en eso y hacer clic en descargar, lo que va a descargarlo a mi computadora local. Y puedo simplemente abrirlo fácilmente en el navegador web.

Nextflow puede generar un reporte como este para cualquier pipeline y tiene información realmente agradable. Así que es una buena práctica siempre guardar estas cosas. Nos dice cuándo ejecutamos, dónde ejecutamos, si fue exitoso o no, qué parámetros se usaron, cuál fue el comando CLI, cosas así.

Y también hay estos gráficos sobre uso de recursos. Así que nos dice qué porcentaje de llamadas de CPU se usaron para cada proceso como un diagrama de caja aquí, porque hay muchas tareas para cada proceso, así que podemos ver la distribución.

Puedes ver nuestros procesos aquí, cowpy y collectGreetings solo tenían una sola tarea, así que es solo una sola línea. Y tenemos tanto CPU como memoria y duración del trabajo, y fueron muy rápidos.

Si estás usando Seqera Platform, por cierto, obtienes los mismos gráficos incorporados en la interfaz de Platform sin tener que hacer nada. Así que siempre tienes esta información al alcance de tu mano.

Bien, así que podemos usar este reporte y en una ejecución real, y tener una idea de cuántos CPUs y cuánta memoria está siendo usada por nuestro pipeline y volver y poner esos valores de vuelta en nuestro archivo de configuración, para que la próxima vez tal vez no solicitemos tanto. Y podamos ser un poco más eficientes.

Ahora puedes ponerte realmente inteligente sobre configurar archivos de configuración del pipeline. Y de nuevo, si estás usando Seqera Platform, busca un pequeño botón que se vea como una bombilla. Porque si haces clic en eso, generará un archivo de configuración altamente optimizado, que está adaptado específicamente para tus datos, tu ejecución y tu pipeline. Para ejecutarlo de la manera más eficiente posible.

Pero por ahora, voy a decir que en realidad el número predeterminado de CPUs que Nextflow estaba dando estaba bien y solo necesito un gigabyte de memoria.

## 5.3. Establecer asignaciones de recursos para un proceso específico

Ahora, en la vida real, es bastante inusual que todos los procesos en tu pipeline vayan a necesitar los mismos requisitos. Podrías tener algo como MultiQC como una herramienta de reportes, que necesita muy poco en términos de recursos y se ejecuta bastante rápido.

Y luego tal vez tienes algo que está indexando un genoma de referencia o haciendo algún alineamiento o haciendo algún otro trabajo. No importa qué es, que toma muchos recursos. Y así para estos diferentes envíos de trabajos a un programador, quieres dar diferentes cantidades de recursos.

Bajo este ámbito process, podemos definir un config, que apunta a procesos específicos de diferentes maneras.

Aquí estamos usando withName, también podemos usar labels, y estos pueden usar un patrón para apuntar a uno o múltiples procesos. Aquí solo estamos diciendo cualquier proceso que tenga un nombre cowpy establecido a dos gigabytes de memoria y dos CPUs, y porque este es un selector más específico que el process de nivel superior, esto es sobrescrito en estos casos, así que puedes construir un bonito archivo de configuración aquí, que realmente adapta todos tus diferentes procesos en tu pipeline para hacerlos realmente eficientes.

## 5.5. Añadir límites de recursos

Ahora como desarrollador de pipeline, probablemente conozco las herramientas bastante bien, y quiero que todo se ejecute lo más rápido y lo mejor posible. Así que podría ser que ponga números bastante altos para algunos de estos porque sé que se ejecutará mucho más rápido si le doy a cowpy 20 CPUs en lugar de dos.

Eso está bien hasta que vas a ejecutar en tu computadora portátil o en GitHub Actions Continuous Integration test, o algún otro sistema, que tal vez no tenga 20 CPUs disponibles.

Ahora cuando intentas ejecutar el pipeline, fallará porque Nextflow dirá, no puedo enviar este trabajo a ninguna parte. No tengo los recursos disponibles.

Ahora para evitar ese fallo duro, podemos añadir un poco más de configuración, que es específica de nuestro sistema ahora, llamada límites de recursos. Y eso se ve así. Está bajo el ámbito process de nuevo.

Y límites de recursos, puedes especificar básicamente el techo de lo que tienes disponible. Es un mapa aquí, y puedes, dentro de este mapa, puedes establecer la memoria, los CPUs, y el tiempo.

Ahora lo que sucede es cuando Nextflow envía una tarea de un proceso, mira lo que se solicita y básicamente solo hace un mínimo entre eso y eso. Así que si solicitamos 20 CPUs, pero solo cuatro están disponibles, solicitará cuatro. El pipeline no falla y usa lo más cercano a lo que fue diseñado por el desarrollador del pipeline como sea posible.

## 6. Usar perfiles para cambiar entre configuraciones preestablecidas

Bien. Dije que los límites de recursos aquí podrían ser específicos del sistema, y tal vez tengo un archivo Nextflow config en mi pipeline, y sé que la gente va a estar usando esto en una variedad de lugares diferentes. Ahora, en lugar de forzar a todos a crear su propio archivo Nextflow config cada vez, lo que puedo hacer es que puedo agrupar diferentes preajustes de configuración juntos en perfiles de configuración.

Voy a desplazarme un poco hacia abajo aquí y luego justo después de params, porque el orden del archivo de configuración aquí es importante, el archivo de configuración se carga secuencialmente, así que voy a poner estos perfiles después de todo lo demás para que sobrescriba los params previamente definidos. Y voy a pegar estos perfiles del material de capacitación.

Así que hay un nuevo ámbito de nivel superior llamado profiles. Podemos tener nombres arbitrarios aquí. Así que tenemos my_laptop y univ_hpc. Y aquí podemos ver que estamos estableciendo los otros mismos parámetros de configuración que teníamos antes. Ahora dentro de solo un perfil. Así que tenemos un executor local para ejecutar en my_laptop y estoy enviando a un clúster SLURM en el HPC.

Estoy usando Docker localmente, conda en el HPC, y el sistema HPC tiene límites de recursos mucho más altos.

Ahora puedo ejecutar el pipeline con la opción CLI -profile, decir qué perfil quiero usar. Así que voy a usar my_laptop, y Nextflow aplicará toda la configuración dentro de ese ámbito de perfil. Así que puedo intentar eso ahora. Es el mismo comando de antes. Nextflow run hello config, y hago guion profile, guion simple porque es la opción central de Nextflow, guion profile my_laptop.

Ahora va a aplicar toda esa opción de configuración por lotes. Oh, y puedes ver, dije antes que esto podría suceder que el requisito del proceso, solicitó cuatro CPUs y solo tengo dos en esta instancia de Codespaces.

Así que esta es una buena oportunidad solo para probar los límites de recursos del proceso, y decir que solo tengo dos CPUs en my_laptop, o en estos Codespaces. Ahora si lo ejecutamos de nuevo, debería limitar ese requisito a dos y con suerte el pipeline se ejecutará. Genial.

## 6.2. Crear un perfil de parámetros de prueba

Nota que estos perfiles no tienen que tener solo configuración sobre su infraestructura. Puedes tener agrupaciones de cualquier configuración aquí, incluyendo parámetros.

Así que otra cosa que verás muy a menudo en los pipelines de la gente es un perfil test, que incluye parámetros, que normalmente enviarías por usuario. Pero aquí tenemos, básicamente diferentes valores predeterminados sensatos para cuando quiero ejecutar casos de prueba.

Y esto es genial porque no tengo que necesariamente ir y especificar todas estas cosas, que podrían ser parámetros requeridos. De otra manera puedo simplemente decir guion profile test y simplemente se ejecutará sin problemas.

Ahora algo a notar es que los perfiles también pueden ser combinados más de uno. Así que puedo hacer profile my_laptop aquí, y luego también añadir test. No hago profile dos veces. Solo hago una lista separada por comas aquí sin espacios. Y va a aplicar estos perfiles en orden. Así que tomará la configuración del perfil my_laptop, y luego aplicará la configuración de test encima.

Realmente conveniente y puedes ver cómo puedes configurar muchos grupos de valores predeterminados sensatos aquí para hacer que sea fácil ejecutar tu pipeline.

## 6.3. Usar nextflow config para ver la configuración resuelta

Con suerte, te he convencido de que la resolución de configuración de Nextflow es poderosa, pero no te culparía si te estás poniendo un poco bizco en este punto después de que he dicho unas 20 formas diferentes de proporcionar configuración y dar todas estas diferentes capas como una piel de cebolla.

Así que si alguna vez te sientes inseguro sobre cuál es la configuración final resuelta para Nextflow, debes saber que hay un comando llamado "nextflow config", y podemos ejecutar eso y nos dirá cuál es la configuración resuelta en nuestra ubicación actual.

Así que cuando lo ejecuto aquí, encuentra el archivo "nextflow.config" en el directorio de trabajo actual, y procesa toda la configuración diferente, y me da la salida resuelta.

Nota que el archivo Nextflow config también puede tomar la opción CLI de perfil. Así que si le digo que resuelva en los perfiles my_laptop y test, y puedes ver que también aplicó los límites de recursos aquí de la opción de configuración my_laptop y también estableció los params, que estaban en el test.

Así que esta es una buena manera solo para explorar cómo está funcionando la resolución de configuración, si tienes alguna duda.

## Conclusión

Bien, eso es todo. Eso es Nextflow config en pocas palabras. Puedes hacer muchas cosas con config. Es realmente poderoso. Pero estos son la mayoría de los casos de uso comunes que te encontrarás haciendo, y estos conceptos se aplican a todas las diferentes opciones.

Date una palmadita en la espalda porque este es el final del curso de capacitación Hello Nextflow. Con suerte ahora estás confiado tanto escribiendo tu propio pipeline de Nextflow desde cero, configurándolo y ejecutándolo, y conoces todos los detalles y las cosas a tener en cuenta.

Hay un cuestionario más que puedes probar en la página de capacitación de configuración. Así que baja y prueba eso y asegúrate de haber entendido todas estas partes sobre la configuración.

Y únete a nosotros en el último video solo para una conclusión rápida sobre algunos de los próximos pasos que podrían ser buenos después de este curso de capacitación.

Gracias por quedarte con nosotros. Bien hecho y te veré en el próximo video.
