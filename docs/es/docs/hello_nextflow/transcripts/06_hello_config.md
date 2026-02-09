# Parte 6: Hello Config - Transcripción del Video

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Traducción asistida por IA - [más información y sugerencias](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

<div class="video-wrapper">
  <iframe width="560" height="315" src="https://www.youtube.com/embed/FcZTiE25TeA?si=tnXTi6mRkITY0zW_&amp;list=PLPZ8WHdZGxmWKozQuzr27jyMGqp9kElVK&amp;cc_load_policy=1" title="YouTube video player" frameborder="0" allow="accelerometer; autoplay; clipboard-write; encrypted-media; gyroscope; picture-in-picture; web-share" referrerpolicy="strict-origin-when-cross-origin" allowfullscreen></iframe>
</div>

!!!note "Notas importantes"

    Esta página muestra únicamente la transcripción. Para instrucciones completas paso a paso, regrese al [material del curso](../06_hello_config.md).

    Los números de sección mostrados en la transcripción se proporcionan solo con fines indicativos y pueden no incluir todos los números de sección en los materiales.

## Bienvenida

Hola, y bienvenido/a de nuevo a la Parte Seis de Hello Nextflow. Esta sección trata sobre configuraciones, y es la última parte de este curso.

Nextflow es particularmente bueno en dos cosas: reproducibilidad y portabilidad. Las configuraciones es donde vemos brillar realmente la segunda de estas. La capacidad de configurar un pipeline de Nextflow para ejecutarse de diferentes maneras y funcionar en diferentes sistemas, sin tener que editar el código subyacente del pipeline.

Este superpoder permite que los pipelines de Nextflow sean reutilizados por otras personas en diferentes lugares, o a través de diferentes infraestructuras a las que usted mismo podría tener acceso.

Significa que puede desarrollar código de pipeline en su computadora portátil, enviarlo a la nube, ejecutarlo en su HPC, y es el mismo código de pipeline y se ejecuta en todas partes.

En esta sección, vamos a repasar algunos temas. Comenzaremos con cómo Nextflow maneja los archivos de configuración, de dónde los carga, y cómo los escribe y cómo los estructura, y esa separación entre el pipeline mismo y lo que debería ir en un archivo de configuración.

Luego pasaremos a algunos casos de uso comunes como cambiar dónde se almacenan los archivos de salida, y también cómo hacer que el pipeline funcione en diferentes infraestructuras, tanto usando diferentes tipos de empaquetado de software como enviando trabajos a diferentes infraestructuras.

## Jerarquías de archivos de configuración

Bien, comencemos. Cuando se trata de cargar archivos de configuración, Nextflow puede obtenerlos de muchos lugares diferentes, lo cual es algo bueno y también puede ser algo un poco arriesgado porque a veces puede ser un poco difícil saber de dónde está obteniendo un archivo de configuración y en qué orden carga las cosas.

Así que realmente recomiendo que haga clic en este enlace aquí, que nos lleva a la documentación de Nextflow. Y en esta página de configuración, enumera los lugares clave desde donde se carga la configuración, y lo que es importante, el orden en que se cargan estas cosas.

Entonces puede ver, puede poner un archivo de configuración en su directorio home de Nextflow, que típicamente es ".nextflow" en su directorio home. Y ese archivo siempre será cargado por cada ejecución de Nextflow en su sistema.

El siguiente lugar a buscar es un archivo en la raíz de su repositorio o directorio de pipeline llamado "nextflow.config".

Luego después de eso, otro archivo llamado "nextflow.config", pero esta vez en el directorio desde el cual está lanzando Nextflow: el directorio de lanzamiento.

Finalmente, puede proporcionar rutas de archivos de configuración en la línea de comandos con un argumento "-c", y puede hacer eso múltiples veces. Y se aplican en el orden en que los especifica.

Puede proporcionar archivos de configuración en todas estas ubicaciones si lo desea, y se cargarán iterativamente, cada uno sobrescribiendo el anterior solo en los ámbitos de configuración donde entren en conflicto.

Este es un sistema realmente poderoso porque significa que puede establecer valores predeterminados sensatos y luego volverse gradualmente más y más específico a medida que se enfoca en esa configuración.

## 0. Calentamiento: Ejecutar hello-config.nf

Bien, cerremos esto y saltemos a nuestros Codespaces y comencemos. Como antes, he limpiado aquí, he eliminado mis directorios de resultados anteriores, mi Nextflow y mis directorios de trabajo y demás. No se preocupe si todavía tiene esos archivos por ahí. Es solo porque estoy muy ampliado y entonces las cosas se ensucian muy rápidamente de lo contrario.

Vamos a trabajar con hello-config.nf, el último archivo en nuestro directorio, y esto debería continuar desde donde lo dejamos en la sección anterior.

Entonces tenemos nuestros cuatro procesos diferentes, que están siendo incluidos desde archivos de módulo. Tenemos nuestros parámetros de pipeline, nuestro bloque de workflow donde estamos llamando a los diferentes procesos y uniendo los canales, publicando los canales de salida, y luego el bloque de output en la parte inferior donde definimos dónde deben almacenarse esos archivos y cómo deben copiarse.

También ya tenemos un archivo "nextflow.config" del último capítulo, donde habilitamos Docker, y vamos a construir sobre este archivo hoy.

Como antes, hemos cambiado la ruta de salida en este script principal a hello config, solo para que no entre en conflicto con resultados anteriores que haya generado.

Bien, solo verifiquemos rápidamente que todo sigue funcionando como esperamos. Abro una terminal y hago nextflow run hello-config.nf. Nextflow se carga. Debería ejecutar nuestros cuatro procesos diferentes. Generar algo de arte ascii agradable usando cowpy y luego guardar nuestros resultados en nuestros archivos de resultados en ese directorio.

Puedo echar un vistazo rápido aquí solo para asegurarme de que estos archivos se vean como esperamos, y efectivamente, ahí está nuestro pavo gigante. Genial.

## 1.1. Mover valores predeterminados a nextflow.config

Ahora lo primero que vamos a hacer es comenzar a mover algunas cosas de nuestro script a nuestro archivo de configuración.

Y lo que nos importa son principalmente los parámetros en esta etapa. Queremos llevar los valores predeterminados al archivo de configuración, para que sea más claro cuáles son los valores predeterminados y sea más fácil para las personas sobrescribirlos.

Voy a tomar este bloque de params aquí del script y ponerlo en el archivo de configuración. Y necesitamos tener un poco de cuidado aquí, porque ahora mismo la sintaxis es ligeramente diferente entre configuración y scripts. El archivo de configuración no puede tomar declaraciones de tipo porque realmente no estamos definiendo estos params, solo los estamos referenciando. Así que voy a deshacerme de esos.

Pero por lo demás es muy similar. Tenemos un bloque de params y luego tenemos nuestros diferentes parámetros de entrada, parámetro batch, parámetro character.

Ahora puedo volver a mi script y ya no necesito definir estos valores predeterminados porque estos valores ahora están en mi archivo Nextflow config.

Sin embargo, dejo los nombres de los parámetros y sus tipos, para que Nextflow conozca esa información y aún pueda hacer toda la seguridad de tipos y todo.

Bien. Guardamos esos archivos y verificamos rápidamente que todo sigue funcionando igual que antes. No debería haber cambios aquí. Hemos mantenido los valores iguales. Solo hemos movido dónde se han definido.

Genial.

## 1.2. Usar un archivo de configuración específico de ejecución

Ahora, hasta ahora hemos estado lanzando Nextflow desde el mismo directorio donde tenemos nuestro script de pipeline. Entonces nuestro directorio de lanzamiento y nuestro directorio de pipeline son básicamente lo mismo.

Para mostrar cómo podemos tener diferentes archivos de configuración con diferentes directorios de lanzamiento, vamos a crear un nuevo subdirectorio ahora.

Así que voy a decir mkdir, y lo vamos a llamar tux-run.

Y luego voy a hacer cd, cambiar de directorio a tux-run. Y note que ahora estamos en nuestro directorio de trabajo ya no está en el mismo directorio que los scripts de pipeline.

Bien, creemos un nuevo archivo "nextflow.config". Entonces touch Nextflow config, y simplemente ábralo en VS Code. También puede ver en la barra lateral aquí que ahora estamos en este subdirectorio.

Ahora podemos tomar el mismo bloque de params que teníamos en nuestro nextflow.config de nivel superior, copiar esto y ahora podemos cambiar estos valores.

Primero, los datos ahora son una ruta relativa diferente porque estamos en un subdirectorio, así que necesitamos actualizar eso. Y luego vamos a cambiar batch a experiment, y vamos a cambiar el character de Turkey a tux.

Ahora hago clic en guardar ahí, y probémoslo. Al igual que con data, ahora necesito decir ../ para llegar al script. Entonces es Hello config. Y presiono enter.

El código del pipeline no ha cambiado en absoluto, pero ahora vamos a tener dos conjuntos de configuración cargándose, y el archivo de configuración del directorio de lanzamiento debería sobrescribir los valores predeterminados, que se establecieron en el nextflow.config del pipeline, y deberíamos obtener diferentes conjuntos de resultados.

Efectivamente, dentro de nuestro directorio aquí, dentro de tux-run, puede ver que tenemos un directorio punto Nextflow y un directorio work y eso es porque estos se crean siempre en su directorio de lanzamiento. Entonces estos son diferentes a los directorios work y results que teníamos de ejecuciones anteriores.

Ahora, si miro en results, podemos ver nuestro collected y ahí está nuestro pequeño personaje tux. Entonces puede ver que esos parámetros fueron interpretados correctamente.

## 1.3. Usar un archivo de parámetros

Bien. Antes cuando estaba hablando sobre los diferentes archivos de configuración que podían cargarse, omití otro lugar del que podemos obtener configuración.

Puede obtenerla de una línea de comandos como hemos visto con guión guión nombres de parámetros, pero también podemos proporcionar un archivo YAML o JSON, solo de params.

El archivo de configuración puede tener todos los diferentes tipos de ámbitos, pero estos archivos son solo parámetros, y es una forma agradable y amigable para el usuario de proporcionar muchos parámetros a la vez, y quizás una forma un poco más reproducible porque los escribe en un archivo, por lo que es fácil obtenerlos en una etapa posterior.

Así que volvamos a nuestra terminal y justo antes de que olvidemos, asegurémonos de retroceder un directorio, así que ya no estoy en el subdirectorio, y voy a mirar el archivo YAML que tenemos aquí llamado test-params.yaml.

Entonces si solo hago code test-params.yaml, puede ver que esto es solo un archivo YAML regular. Nada especial al respecto. Con las claves siendo nuestros nombres de parámetros, con el formato YAML así que dos puntos aquí, y luego un valor.

Note que esto no es código Nextflow, así que no podemos poner cosas como variables aquí. Estos son solo valores estáticos.

También porque JSON en realidad se analiza como YAML, también podemos tener un archivo test-params.json, que se ve muy similar. Es solo un formato de datos diferente.

Entonces tenemos dos archivos de prueba diferentes aquí y tenemos variables ligeramente diferentes.

Bien, entonces ¿cómo le damos estos a Nextflow? Es muy simple. Hacemos Nextflow run hello config, como antes. Y en lugar de "-c" para archivo de configuración, o cargar esos nombres de archivo predeterminados, hacemos -params-file. Guión simple porque es una opción central de Nextflow.

Y luego pasamos la ruta para ese archivo. Así que voy a hacer "-params-file test-params.yaml", y veremos si esos se cargan correctamente.

Bien. Se ejecutó. Recordémonos qué había en este archivo YAML. Entonces el batch se estableció en YAML, así que así es como debería llamarse, y debería tener un stegosaurus. Así que subamos y miremos en results. Y tenemos COLLECTED-yaml. Así que veamos si tenemos un Stegosaurus. Fantástico, un Stegosaurus usando un sombrero. Eso es lo que nos gusta.

Entonces eso funcionó muy bien, y es exactamente lo mismo con el archivo JSON. Solo cambiamos la extensión del archivo aquí y Nextflow sabe cómo leer eso.

Y en este caso, deberíamos tener un batch llamado JSON y deberíamos tener una tortuga. Echemos un vistazo. Maravilloso. Una de mis herramientas CLI favoritas.

## 2.1. Personalizar el directorio de salida con -output-dir

Bien, entonces eso ha sido principalmente pensar en las entradas al pipeline y cambiar parámetros. ¿Qué hay de las salidas?

Ahora, aunque hemos estado cambiando los subdirectorios usando params, puede haber notado que todos nuestros archivos todavía van a results.

Podemos cambiar ese directorio base al que se publican todos los archivos con una bandera de línea de comandos llamada -output-dir. Entonces si hago Nextflow run hello config, y luego hago -output-dir, y lo vamos a llamar "custom-outdir-cli". No puedo escribir. Solo para que recordemos de dónde vinieron estos archivos.

Esta es una opción central de Nextflow y es muy nueva. Esto solo se agregó recientemente, y esta es una de las cosas que podemos hacer con el nuevo analizador de lenguaje y todo.

Es un poco largo de escribir. También puede simplemente llamarlo "-o" si lo desea. Entonces si solo retrocedo. Puedo simplemente acortar eso a "-o", que es un poco más simple.

Bien. Ejecutamos eso. No hemos cambiado nada en nuestro pipeline o incluso en nuestra configuración en este punto, y con suerte debería guardar todos nuestros resultados en un directorio de nivel superior diferente. Y puede imaginar que puede establecer esto en básicamente cualquier ruta que desee.

Acaba de llegar arriba. Tenemos un custom-outdir-cli, y todos los archivos están organizados allí exactamente de la misma manera, con sus mismos subdirectorios y nombres de archivo. Entonces esta es una forma realmente fácil de simplemente cambiar dónde el pipeline publica sus resultados, sin pensar demasiado en cómo se organizan esos resultados.

## 2.1.2. Eliminar rutas codificadas del bloque output

Si miro en este directorio, podemos ver que todavía tenemos un subdirectorio llamado Hello Config, que se siente un poco redundante ahora.

Así que simplemente carguemos nuestro script nuevamente y ahora podemos eliminar ese subdirectorio del bloque output en la parte inferior. Porque realmente ya no lo necesitamos. Así que podemos hacer eso ahora, eliminar eso de aquí. Y luego si es solo esto, puede eliminarlo completamente o dejarlo como una cadena vacía. Voy a dejarlo como una cadena vacía por ahora, porque vamos a volver y poner algunas cosas diferentes en su lugar en el futuro. Pero si no le importan los subdirectorios, es más limpio simplemente eliminar completamente la declaración de path allí.

Bien, guardemos. Solo probemos rápidamente de nuevo. En realidad voy a eliminar mi directorio "custom-outdir-cli" para que no nos confundamos con ningún archivo existente allí. Porque recuerde, cuando publica cosas, no elimina los archivos que ya estaban allí. Solo agrega nuevos. Ejecutemos ese comando nuevamente, custom-outdir-cli.

Y ahora si hace "ls custom-outdir-cli", ya no hay más directorio allí llamado Hello Config.

## 2.2.1. Establecer outputDir en el archivo de configuración

Bien, la bandera de línea de comandos aquí, "-o" o "-output-dir" es buena. ¿Pero qué hay de establecer valores predeterminados para esto en la configuración? ¿Cómo hacemos eso?

Abro el archivo "nextflow.config", cierro todo lo demás y me deshago de eso. Podemos agregar una nueva opción de configuración aquí, que acabo de copiar del sitio web del material de capacitación, y se llama outputDir.

No está bajo ningún ámbito. No está bajo params ni nada. Es de nivel superior, y podemos establecer esto en una cadena. Ahora algo simple de hacer es simplemente cambiarlo a cualquier cosa que no sea results como una cadena codificada. Pero como esto está en un archivo de configuración de Nextflow, podemos ser un poco inteligentes aquí y también incluir variables.

Y puede ver aquí que hemos incluido una variable params, params.batch, que es parte de esta cadena. Esto significa que podemos reutilizar variables que vienen de otros lugares. Y en este caso, si hacemos --batch, cuando ejecutamos Nextflow Pipeline, vamos a obtener un subdirectorio en nuestra ruta personalizada basado en cuál era el nombre del batch.

Bien, entonces probemos esto y solo echemos un vistazo rápido para ver cómo se ven los resultados. Entonces si hago Nextflow run hello config y --batch my_run. Recordémonos cómo se veía la configuración. Entonces es custom-outdir-config.

Tree custom-outdir-config. Y puede ver que el batch se llamaba my_run. Y luego tenemos ese subdirectorio llamado my_run. Entonces esa ruta de archivo dinámica funcionó.

Y no solo eso, ya no fue a un directorio results predeterminado, y no tuve que especificar nada en la línea de comandos para cambiar el directorio base. Entonces hemos restablecido exitosamente el valor predeterminado para el outputDir predeterminado.

## 2.2.2. Subdirectorios con nombres de batch y proceso

Bien, llevemos eso un poco más lejos. Esa es una variable dinámica dentro del archivo de configuración. ¿Qué hay del script mismo? Ahora, hasta ahora hemos tenido estas rutas aquí y estas también pueden ser dinámicas. Entonces en lugar de solo codificar algo, podemos poner algunas llaves onduladas y poner algo dinámico.

Entonces, por ejemplo, tenemos nuestros procesos llamados sayHello. Podríamos hacer sayHello.name, que es un atributo del proceso, que es un poco aburrido porque es solo "sayHello" en este caso. Pero es variable.

Entonces esto le da una idea. Entonces podemos poner esto aquí y decir convertToUpper.name, collectGreetings.name, collectGreetings.name nuevamente, y cowpy.

Ahora cuando ejecutemos, el directorio base todavía va a ser custom-outdir-config. Y va a estar en un subdirectorio llamado params.batch, pero los subdirectorios bajo eso deberían estar organizados por nombre de proceso.

Probemos eso y veamos si funciona. Entonces voy a eliminar el directorio anterior para que no nos confundamos, y solo usar exactamente el mismo comando Nextflow Run.

Debería ejecutarse de la misma manera. Podría estar usando dash resume en todos estos para hacerlo un poco más rápido y usar los resultados calculados previamente. Ahora, si hago tree custom-outdir-config, puede ver que no está en results, está en nuestro directorio base con el nombre del batch. Y puede ver que todos los resultados ahora están organizados dentro de subdirectorios nombrados según el proceso. Entonces tenemos dos lugares diferentes donde estamos definiendo rutas de salida dinámicas aquí.

Bien. Última cosa, agreguemos de nuevo esas carpetas intermedias, que teníamos antes porque eran algo agradables. Intermediates.

Y también podemos pensar un poco sobre este params.batch, tal vez como desarrollador de pipeline realmente me gustaba tener eso en el subdirectorio, pero si los usuarios finales del pipeline están estableciendo "-o" o -output-dir en el CLI, está sobrescribiendo completamente esta declaración completa, y perdemos ese subdirectorio.

Entonces lo que podemos hacer es tomar esa ruta dinámica fuera del outputDir config, que sería eliminada, y ponerla en la ruta de output, que no es eliminada.

Entonces podemos hacer params.batch barra intermediates barra sayHello.name, y hacer todo esto en una cadena entre comillas dobles, para que sea interpolada por Nextflow.

Ahora puedo copiar, ups. Copiar estos hacia abajo a los otros procesos. Recuerde ponerlos todos entre comillas. Y eliminar intermediates de estas salidas particulares.

¿Bien? Se está viendo un poco más complejo ahora, pero puede ver que realmente estamos comenzando a construir una estructura de directorio de salida bien organizada en nuestro código.

Y lo que es realmente agradable es que esta complejidad adicional en el código no pasa al CLI. Entonces podemos ejecutar nuestro comando con -output-dir y cualquier variable batch, solo pensando en cómo ejecutar el pipeline y sin pensar realmente demasiado en lo que hay en el código. Y nuestros archivos de salida van a ser construidos muy bien de una manera muy bien organizada, lo cual es agradable para las personas que usan el pipeline básicamente.

Genial. Mientras escribo esto, me doy cuenta de que cometí un error. Veamos si alguien me atrapó aquí. Tenemos collectGreetings.name, entonces algo salió un poco mal. Y sí, efectivamente, accidentalmente olvidé poner estos entre llaves onduladas.

Entonces recuerde, tenga cuidado cuando esté escribiendo su código y asegúrese de decirle a Nextflow qué es una variable y qué es solo una cadena. Porque hará exactamente lo que le diga que haga. Y nada más. Como todas las buenas computadoras. Bien, eso debería arreglarlo.

## 2.3. Establecer el modo de publicación a nivel de workflow

Hay una parte de este script que todavía no me encanta, que es el hecho de que estamos escribiendo mode copy una y otra vez, y si hay algo que no nos gusta, es repetirnos.

Entonces podemos limpiar esto un poco tomando esto y moviéndolo a la configuración. Y de hecho, podemos establecerlo para todo el pipeline de una sola vez. Entonces no tenemos que decirlo múltiples veces.

Vamos a nuestro archivo de configuración y tenemos un nuevo ámbito aquí llamado workflow. Y podemos hacer llaves onduladas o podemos hacer notación de punto. No hace ninguna diferencia. El sitio web del material de capacitación usa notación de punto. Puedo decir output y podemos mezclar y combinar, entonces mode equals copy. Genial.

Y ahora podemos volver aquí y eliminar estos. Ahora podríamos dejarlos en su lugar. La configuración básicamente está sobrescribiendo lo que está escrito aquí, pero como lo tenemos en la configuración a nivel de pipeline, y estos dos archivos se envían juntos, realmente no hay razón para hacerlo dos veces.

Bien. Solo una verificación de cordura, porque aparentemente cometemos errores. Ejecutemos eso de nuevo y solo verifiquemos que estamos usando correctamente el modo copy para publicar archivos. Entonces vamos a ejecutar el script nuevamente y esta vez hemos puesto los resultados en un directorio llamado config-output-mode, veamos cómo se ven los archivos allí.

Y luego si hago "ls -l" para mirar batch, y podemos mirar cowpy, por ejemplo. Y deberíamos ver, sí, que este es un archivo apropiado aquí, que no es un enlace simbólico, entonces ese atributo de configuración se ha aplicado correctamente.

## 3. Seleccionar una tecnología de empaquetado de software

Bien. Hasta ahora nos hemos estado enfocando en las entradas y las salidas, los archivos con los que el workflow está ejecutándose. ¿Pero qué hay de la infraestructura? Dije al principio que Nextflow le permite ejecutar el mismo pipeline en diferentes configuraciones de computación. Entonces, ¿cómo se ve eso?

Para mostrar esto, vamos a cambiar de usar Docker para ejecutar cowpy, y en su lugar usaremos Conda para hacer lo mismo.

Puedo hacer esto muy simplemente. Si voy a code, "nextflow.config". Si recuerda en la parte superior, definimos docker.enabled anteriormente, y en el último capítulo para que pudiéramos usar el contenedor con cowpy.

Voy a decirle a Nextflow que no use Docker. Establecer eso en false. Y voy a decir Conda enabled equals true. Entonces decirle a Nextflow, por favor use Conda.

Ahora solo habilitar Conda no es suficiente por sí mismo. Así como hicimos con Docker, tenemos que decirle a Nextflow dónde puede obtener el software que necesita.

Entonces si saltamos a los módulos aquí. Y abrimos el script cowpy. Podemos ver que tenemos una declaración de container en la parte superior. Y el contenedor es usado por Docker, pero también Singularity, Apptainer, y muchas de las otras herramientas de software.

Pero no puede ser usado para Conda, así que tenemos una declaración separada llamada "conda", y podríamos simplemente escribir "cowpy". Y eso dejará que la resolución de paquetes de conda descubra la mejor manera de resolver eso según su entorno conda local.

O es una buena práctica hacer lo que dice el sitio web del material de capacitación, que es definir un canal conda específico con su notación de doble dos puntos, y definitivamente definir una versión específica del software para que cada persona que ejecute el pipeline obtenga la misma versión.

Note que los contenedores son un poco superiores en este aspecto, porque cuando instala algo con Conda, todavía va a descubrir todas las dependencias para ese paquete, y pueden cambiar con el tiempo. Llamado deriva de dependencias.

Entonces los contenedores, sin embargo, bloquean toda la pila de dependencias de software hasta el final, por lo que puede estar un poco más seguro de que A, va a funcionar, y B, será reproducible.

Entonces si puede usar Docker o Singularity o Apptainer, definitivamente lo recomendaría.

Ahora lo que es agradable de esto es que el archivo de módulo, que está escrito por el desarrollador del pipeline, ahora tiene tanto Container como Conda, y entonces le estamos diciendo a la persona que está ejecutando este pipeline, no nos importa qué solución de empaquetado de software use. Funcionará tanto con Docker como con Conda, y aquí es donde obtener el software en ambos casos.

Podemos abrir la terminal y probemos esto. Entonces Nextflow run hello config --batch conda. Y la primera vez que esto se ejecuta con conda, va a ser un poco lento cuando llegue a ese proceso particular, porque tiene que ejecutar "conda install".

Y está creando un entorno conda especial solo para este proceso. Entonces no está usando mi entorno conda global, que tengo en mi terminal. Está creando uno solo para ese proceso. Esto es bueno porque evita cosas como conflictos de dependencias entre diferentes procesos en su workflow. Si sus procesos tienen herramientas que necesitan diferentes versiones de Python o cosas así, está bien porque están usando diferentes entornos conda.

Nextflow almacena en caché estos entornos conda localmente, puede ver que le dice dónde está esa ruta, está en el directorio work aquí. Y entonces la próxima vez que ejecute este script con Conda, será mucho más rápido porque encontrará ese entorno conda existente y simplemente lo reutilizará. Pero la primera vez que lo hacemos, tiene que ir y buscarlo, resolverlo, descargar todas las dependencias y configurar todo.

Bien, genial, se ejecutó. Podemos simplemente recordarnos para qué está configurado actualmente el pipeline. Si miramos en el archivo de configuración, era "custom-outdir-config" ahora mismo para mí. Veamos si subo a ese directorio base. E hice --batch conda. Ahí está nuestro subdirectorio conda. Entonces funcionó y ahí está nuestra salida de cowpy.

Entonces buscó cowpy, lo instaló en mi sistema local usando conda, y ejecutó el proceso. Y lo que es genial es que, como ese usuario final, no tuve que pensar en absoluto sobre ninguna gestión de software allí. Nextflow simplemente lo resolvió para mí. Dije, necesito usar conda en este sistema. El desarrollador del pipeline dijo qué paquetes necesitaba. Y Nextflow hizo el resto. Muy poderoso.

Note que en realidad puede usar una mezcla de diferentes tecnologías. Entonces puedo habilitar Docker para procesos específicos, y conda para otros procesos, o decir que algunos procesos deberían simplemente usar cualquier software local que tuviera instalado. Esto es bastante inusual, pero es posible, y en algunos casos, por ejemplo, si está usando cierto software que podría ser difícil de empaquetar en Docker, tiene una escapatoria.

## 4. Seleccionar una plataforma de ejecución

Entonces eso es empaquetado de software. La otra parte de la portabilidad a otros sistemas es dónde se ejecutan realmente los trabajos. En este momento, estoy ejecutando básicamente en mi computadora portátil o en estos Codespaces, que es una sola computadora. No hay nada sofisticado. Nextflow está siendo un poco inteligente sobre paralelizar los trabajos lo mejor que puede, pero todo está en un sistema.

Ahora, si está ejecutando en un HPC, probablemente tenga algún tipo de programador de trabajos como SLURM o PBS o algo, y enviará trabajos a ese programador y distribuirá todos los trabajos a diferentes nodos de cómputo.

Otra forma de ejecutar es en la nube. Entonces tal vez esté usando AWS Batch, o Azure Cloud, o Google. Y todos estos funcionan en un sistema similar donde tiene un programador y envía trabajos y se envían a diferentes lugares para ser computados.

Ahora en el pasado lejano cuando comencé a hacer bioinformática, el software de todos para ejecutar análisis estaba muy vinculado a su infraestructura computacional, lo que lo hacía casi imposible de replicar.

Pero con esta separación de configuración en Nextflow, y con la capacidad de Nextflow de interactuar con muchos backends de infraestructura de cómputo diferentes, es muy simple tomar nuestro pipeline sin modificar el código del pipeline en absoluto y simplemente cambiar eso.

## 4.1. Apuntar a un backend diferente

Entonces si vamos a nuestro archivo "nextflow.config", y ahora podemos poner algo de configuración a nivel de proceso. Entonces si pongo en la parte superior el ámbito process y puedo establecer el executor, y aquí está establecido en local, que es el predeterminado.

Note que como esto es a nivel de proceso, podemos apuntar cosas a diferentes procesos. Y entonces en realidad puede configurar executors para ser específicos del proceso y tener una ejecución híbrida, donde algunos trabajos podrían ejecutarse localmente, donde sea que se esté ejecutando el trabajo de Nextflow. Algunos se envían a diferentes HPC y algunos podrían enviarse a la nube. Puede ser tan inteligente como desee.

Ahora, es muy difícil demostrar esto en un entorno de capacitación como este porque no tengo un HPC al que enviar. Pero lo que puedo hacer es si escribo slurm, podemos hacer trampa un poco y puede tener una idea de esto.

Y esto es realmente más interesante para las personas que están acostumbradas a ejecutar en SLURM y saben cómo se ven los encabezados de SLURM. Pero si hago Nextflow run, hello config. Va a fallar porque va a intentar enviar trabajos a un clúster que no existe. Entonces obtendremos algún tipo de error sobre sbatch no estar disponible.

Sí, escrito. Esa es la herramienta. Esa es la herramienta CLI que usa para enviar trabajos a un clúster slurm. Pero lo que podemos hacer es ir y mirar en nuestro directorio work aquí haciendo clic con comando, abrir ese directorio y mirar el .command.run. Y puede ver en la parte superior del archivo .command.run, tenemos nuestros encabezados sbatch, diciéndole a un clúster SLURM teórico cómo manejar este envío de trabajo.

Y entonces puede ver que Nextflow está siendo inteligente, está haciendo todas las cosas correctas. Es solo que no teníamos un clúster al que enviar.

## 5. Controlar asignaciones de recursos de cómputo

¿Qué más es diferente entre diferentes infraestructuras de computación? Otra cosa es cuántos recursos disponibles tiene, y de hecho, en muchos entornos de cómputo, es un requisito que tenga que especificar cuántos CPUs y cuánta memoria necesita un trabajo.

Nuevamente, Nextflow abstrae esto para nosotros, de modo que ya no es específico de un solo tipo de entorno de cómputo, y podemos escribir en el ámbito de nivel de proceso aquí. CPUs equals uno, memory equals dos gigabytes. Nuestro pipeline no es muy exigente, así que eso debería estar bien.

Ahora, acabo de adivinar estos números aquí, pero ¿cómo sabe qué es una cantidad sensata de recursos para usar? Es un trabajo bastante difícil ir y revisar todos estos diferentes procesos de un gran pipeline de muchas muestras y entender cuál fue la utilización de recursos.

Entonces un buen enfoque para esto es establecer estos valores en números altos para comenzar, solo para que su pipeline se ejecute sin errores, y luego pedirle a Nextflow que genere un reporte de uso para usted.

Esto es súper fácil de hacer, así que voy a volver a una terminal. Oh, necesito recordar establecer eso de nuevo a local para que mi pipeline realmente se ejecute. Y voy a decir Nextflow run, y voy a usar una bandera de línea de comandos -with-report.

Y puedo dejar eso en blanco y dará un nombre de archivo predeterminado, pero voy a darle un nombre de archivo específico para que se guarde en un lugar específico.

Presiono Enter, y el pipeline se ejecuta exactamente como normal, pero cuando termine, va a generar un bonito reporte HTML para mí.

Entonces en la barra lateral aquí, tengo este archivo HTML. Si estuviera ejecutando esto localmente, simplemente lo abriría. Como estoy en Codespaces, voy a hacer clic derecho en eso y hacer clic en descargar, lo que va a descargarlo a mi computadora local. Y puedo simplemente abrirlo fácilmente en el navegador web.

Nextflow puede generar un reporte como este para cualquier pipeline y tiene información realmente agradable. Entonces es una buena práctica guardar siempre estas cosas. Nos dice cuándo ejecutamos, dónde ejecutamos, si fue exitoso o no, qué parámetros se usaron, cuál fue el comando CLI, cosas así.

Y también hay estos gráficos sobre el uso de recursos. Entonces nos dice qué porcentaje de llamadas de CPU se usaron para cada proceso como un diagrama de caja aquí, porque hay muchas tareas para cada proceso, por lo que podemos ver la distribución.

Puede ver nuestros procesos aquí, cowpy y collectGreetings solo tenían una sola tarea, por lo que es solo una sola línea. Y tenemos tanto CPU como memoria y duración del trabajo, y fueron muy rápidos.

Si está usando Seqera Platform, por cierto, obtiene los mismos gráficos integrados en la interfaz de Platform sin tener que hacer nada. Entonces siempre obtiene esta información al alcance de su mano.

Bien, entonces podemos usar este reporte y en una ejecución real, y tener una idea de cuántos CPUs y cuánta memoria está siendo usada por nuestro pipeline y volver y poner esos valores de nuevo en nuestro archivo de configuración, para que la próxima vez tal vez no solicitemos tanto. Y podemos ser un poco más eficientes.

Ahora puede volverse realmente inteligente sobre configurar archivos de configuración de pipeline. Y nuevamente, si está usando Seqera Platform, busque un pequeño botón que parezca una bombilla. Porque si hace clic en eso, generará un archivo de configuración altamente optimizado, que está adaptado específicamente para sus datos, su ejecución y su pipeline. Para ejecutarlo de la manera más eficiente posible.

Pero por ahora, voy a decir que en realidad el número predeterminado de CPUs que Nextflow estaba dando estaba bien y pero solo necesito un gigabyte de memoria.

## 5.3. Establecer asignaciones de recursos para un proceso específico

Ahora, en la vida real, es bastante inusual que todos los procesos en su pipeline vayan a necesitar los mismos requisitos. Podría tener algo como MultiQC como una herramienta de reporte, que necesita muy poco en términos de recursos y se ejecuta bastante rápido.

Y luego tal vez tenga algo que está indexando un genoma de referencia o haciendo algún alineamiento o haciendo algún otro trabajo. No importa qué sea, que toma muchos recursos. Y entonces para estos diferentes envíos de trabajos a un programador, desea dar diferentes cantidades de recursos.

Bajo este ámbito de proceso, podemos definir una configuración, que apunta a procesos específicos de diferentes maneras.

Aquí estamos usando withName, también podemos usar labels, y estos pueden usar un patrón para apuntar a uno o múltiples procesos. Aquí solo estamos diciendo que cualquier proceso que tenga un nombre cowpy se establezca en dos gigabytes de memoria y dos CPUs, y como este es un selector más específico que el proceso de nivel superior, esto se sobrescribe en estos casos, por lo que puede construir un bonito archivo de configuración aquí, que realmente adapta todos sus diferentes procesos en su pipeline para hacerlos realmente eficientes.

## 5.5. Agregar límites de recursos

Ahora como desarrollador de pipeline, probablemente conozco las herramientas bastante bien, y quiero que todo se ejecute lo más rápido y lo mejor posible. Entonces podría ser que ponga números bastante altos para algunos de estos porque sé que se ejecutará mucho más rápido si le doy a cowpy 20 CPUs en lugar de dos.

Eso está bien hasta que vaya a ejecutar en su computadora portátil o en pruebas de Integración Continua de GitHub Actions, o algún otro sistema, que tal vez no tenga 20 CPUs disponibles.

Ahora cuando intente ejecutar el pipeline, fallará porque Nextflow dirá, no puedo enviar este trabajo a ningún lado. No tengo los recursos disponibles.

Ahora para evitar ese fallo duro, podemos agregar un poco más de configuración, que es específica de nuestro sistema ahora, llamada límites de recursos. Y eso se ve así. Está bajo el ámbito de proceso nuevamente.

Y límites de recursos, puede especificar básicamente el techo de lo que tiene disponible. Es un mapa aquí, y puede, dentro de este mapa, puede establecer la memoria, los CPUs y el tiempo.

Ahora lo que sucede es cuando Nextflow envía una tarea desde un proceso, mira lo que se solicita y básicamente solo hace un mínimo entre eso y eso. Entonces si solicitamos 20 CPUs, pero solo cuatro están disponibles, solicitará cuatro. El pipeline no falla y usa lo más cercano a lo que fue diseñado por el desarrollador del pipeline como sea posible.

## 6. Usar perfiles para cambiar entre configuraciones preestablecidas

Bien. Dije que los límites de recursos aquí podrían ser específicos del sistema, y tal vez tengo un archivo Nextflow config en mi pipeline, y sé que las personas van a estar usando esto en una variedad de lugares diferentes. Ahora, en lugar de forzar a todos a crear su propio archivo Nextflow config cada vez, lo que puedo hacer es agrupar diferentes preajustes de configuración juntos en perfiles de configuración.

Voy a desplazarme un poco hacia abajo aquí y luego justo después de params, porque el orden del archivo de configuración aquí es importante, el archivo de configuración se carga secuencialmente, así que voy a poner estos perfiles después de todo lo demás para que sobrescriba los params definidos previamente. Y voy a pegar estos perfiles del material de capacitación.

Entonces hay un nuevo ámbito de nivel superior llamado profiles. Podemos tener nombres arbitrarios aquí. Entonces tenemos my_laptop y univ_hpc. Y aquí podemos ver que estamos estableciendo los otros mismos parámetros de configuración que teníamos antes. Ahora dentro de solo un perfil. Entonces tenemos un executor local para ejecutar en my_laptop y estoy enviando a un clúster SLURM en el HPC.

Estoy usando Docker localmente, conda en el HPC, y el sistema HPC tiene límites de recursos mucho más altos.

Ahora puedo ejecutar el pipeline con la opción CLI -profile, decir qué perfil quiero usar. Entonces voy a usar my_laptop, y Nextflow aplicará toda la configuración dentro de ese ámbito de perfil. Entonces puedo probar eso ahora. Es el mismo comando que antes. Nextflow run hello config, y hago guión profile, guión simple porque es la opción central de Nextflow, guión profile my_laptop.

Ahora va a aplicar en lote toda esa opción de configuración. Oh, y puede ver, dije antes que esto podría suceder que el requisito del proceso, pidió cuatro CPUs y solo tengo dos en esta instancia de Codespaces.

Entonces esta es una buena oportunidad solo para probar los límites de recursos del proceso, y decir que solo tengo dos CPUs en my_laptop, o en estos Codespaces. Ahora si lo ejecutamos de nuevo, debería limitar ese requisito a dos y con suerte el pipeline se ejecutará. Genial.

## 6.2. Crear un perfil de parámetros de prueba

Note que estos perfiles no tienen que tener solo configuración sobre su infraestructura. Puede tener agrupaciones de cualquier configuración aquí, incluyendo parámetros.

Entonces otra cosa que verá muy a menudo en los pipelines de las personas es un perfil de prueba, que incluye parámetros, que normalmente enviaría por usuario. Pero aquí tenemos, básicamente diferentes valores predeterminados sensatos para cuando quiero ejecutar casos de prueba.

Y esto es genial porque no tengo que necesariamente ir y especificar todas estas cosas, que podrían ser parámetros requeridos. De lo contrario, puedo simplemente decir guión profile test y simplemente se ejecutará sin problemas.

Ahora algo a notar es que los perfiles también pueden combinarse más de uno. Entonces puedo hacer profile my_laptop aquí, y luego también agregar test. No hago profile dos veces. Solo hago una lista separada por comas aquí sin espacios. Y va a aplicar estos perfiles en orden. Entonces tomará la configuración del perfil my_laptop, y luego aplicará la configuración de test encima.

Realmente conveniente y puede ver cómo puede configurar muchos grupos predeterminados sensatos aquí para facilitar la ejecución de su pipeline.

## 6.3. Usar nextflow config para ver la configuración resuelta

Con suerte, los he convencido de que la resolución de configuración de Nextflow es poderosa, pero no los culparía si están un poco bizcos en este punto después de que he dicho unas 20 formas diferentes de proporcionar configuración y dar todas estas diferentes capas como una piel de cebolla.

Entonces si alguna vez se siente inseguro sobre cuál es la configuración resuelta final para Nextflow, sepa que hay un comando llamado "nextflow config", y podemos ejecutar eso y nos dirá cuál es la configuración resuelta en nuestra ubicación actual.

Entonces cuando lo ejecuto aquí, encuentra el archivo "nextflow.config" en el directorio de trabajo actual, y procesa toda la configuración diferente, y me da la salida resuelta.

Note que el archivo de configuración de Nextflow también puede tomar la opción CLI de perfil. Entonces si le digo que resuelva en los perfiles my_laptop y test, y puede ver que también ha aplicado los límites de recursos aquí de la opción de configuración my_laptop y también estableció los params, que estaban en el test.

Entonces esta es una buena manera de simplemente explorar cómo está funcionando la resolución de configuración, si tiene alguna duda.

## Conclusión

Bien, eso es todo. Eso es la configuración de Nextflow en pocas palabras. Puede hacer muchas cosas con la configuración. Es realmente poderoso. Pero estos son la mayoría de los casos de uso comunes que se encontrará haciendo, y estos conceptos se aplican a todas las diferentes opciones.

Dese una palmadita en la espalda porque este es el final del curso de capacitación Hello Nextflow. Con suerte ahora está seguro tanto de escribir su propio pipeline de Nextflow desde cero, configurarlo y ejecutarlo, y conoce todos los detalles y las cosas a tener en cuenta.

Hay un cuestionario más que puede probar en la página de capacitación de configuración. Entonces baje y pruebe eso y asegúrese de haber entendido todas estas partes sobre la configuración.

Y únase a nosotros en el último video solo para una conclusión rápida sobre algunos de los próximos pasos que podrían ser buenos para hacer después de este curso de capacitación.

Gracias por quedarse con nosotros. Bien hecho y lo veré en el próximo video.
