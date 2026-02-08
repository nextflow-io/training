# Parte 5: Hello Containers - Transcripción del video

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Traducción asistida por IA - [más información y sugerencias](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

<div class="video-wrapper">
  <iframe width="560" height="315" src="https://www.youtube.com/embed/Xqr--bKEN9U?si=QinuAnFwFj-Z8CrO&amp;list=PLPZ8WHdZGxmWKozQuzr27jyMGqp9kElVK&amp;cc_load_policy=1" title="YouTube video player" frameborder="0" allow="accelerometer; autoplay; clipboard-write; encrypted-media; gyroscope; picture-in-picture; web-share" referrerpolicy="strict-origin-when-cross-origin" allowfullscreen></iframe>
</div>

!!!note "Notas importantes"

    Esta página muestra únicamente la transcripción. Para instrucciones completas paso a paso, regrese al [material del curso](../05_hello_containers.md).

    Los números de sección mostrados en la transcripción se proporcionan únicamente con fines indicativos y pueden no incluir todos los números de sección en los materiales.

## Bienvenida y antecedentes

Hola, y bienvenidos de nuevo a Hello Nextflow. Esta es la parte cinco llamada Hello Containers. Y en esta parte del curso, vamos a hablar todo sobre cómo encapsular los requisitos de software para un pipeline de manera que las personas que ejecuten el pipeline no tengan que pensar en instalar el software.

Si has estado trabajando en bioinformática tanto tiempo como yo, quizás recuerdes lo que a menudo llamo los viejos malos tiempos, donde cuando querías ejecutar el pipeline de otra persona o replicar su trabajo, pasabas horas o días tratando de instalar todas las diferentes herramientas de software que usaron, en las mismas versiones, tratando de compilarlas en tu máquina, y era una pesadilla. Era realmente difícil.

Si estabas ejecutando en un HPC, quizás usaste módulos de entorno donde los administradores de sistemas trataban de instalar software para ti, lo cual estaba bien, pero todavía era imperfecto.

Pero ahora tenemos mejores formas de hacer esto. Nextflow tiene soporte integrado para diferentes tecnologías de contenedores de software. Docker es la más común. Esa es la que vamos a usar hoy. Funciona bien en Codespaces. Funciona bien en tu computadora local y funciona bien en la nube.

Pero también Singularity o Apptainer, que son muy comunes en sistemas HPC y efectivamente funcionan exactamente de la misma manera. O Poman, Shifter, hay un montón de otros que son todos muy similares.

El único extra, que es algo similar pero no exactamente, que Nextflow soporta es Conda. Y Nextflow puede gestionar entornos Conda para ti por proceso, lo cual es mucho mejor que hacer tus propios entornos Conda. Y nuevamente, puede enviarse con un pipeline.

Vamos a comenzar este capítulo hablando un poco sobre las tecnologías de contenedores y Docker y cómo funcionan. Y vamos a hacer la primera mitad solo manualmente en Docker para que entiendas qué está pasando bajo el capó y cómo funciona esto. Porque eso es realmente importante para entender qué está haciendo Nextflow y cómo entender qué está haciendo tu workflow cuando se está ejecutando.

Así que. Saltemos a nuestros Codespaces. Ahora he limpiado todo de nuevo, pero si saltamos a Hello Containers, deberías ver que todos nuestros scripts y todo están ahí igual que al final del capítulo de módulos. Así que tenemos nuestros diferentes módulos aquí, que creé en el directorio de módulos.

Todavía están ahí. Necesitan estar ahí para que pueda ejecutarse. y el workflow y la salida son todos iguales excepto que hemos cambiado la ruta de publicación de salida a Hello Containers, para que tus archivos terminen en ese directorio.

Podemos ejecutar esto ahora para verificar que funciona si quieres, o podemos continuar con la terminal.

## 1. Usar un contenedor 'manualmente'

Vamos a usar Docker para gestionar nuestros contenedores, y puedo verificar que está instalado en mi Codespaces haciendo "docker -v", lo cual me muestra la versión que está instalada y todo, y que está funcionando correctamente.

Ahora los contenedores y Docker tienen dos conceptos que son realmente importantes. Uno se llama imagen, y uno se llama contenedor. La imagen es la instantánea, si quieres, de todo el sistema de archivos que estarás usando, y el contenedor es el entorno en ejecución. Así que creas un contenedor usando una imagen.

Una vez que estás en ese contenedor, típicamente funciona como un sistema operativo completo. Está aislado del mundo exterior. Está separado de todo lo demás, y eso es algo bueno. Así es como obtenemos tan buena reproducibilidad con Nextflow.

Porque para tareas ejecutadas dentro de un contenedor, no están contaminadas por ningún archivo de configuración en tu sistema local. Ninguna otra influencia externa, se ejecutan en su propia pequeña caja de arena. Los archivos se producen entonces de una manera muy, muy reproducible porque estás usando las mismas bibliotecas subyacentes, todas las mismas dependencias, exactamente el mismo software para cada persona ejecutando en cada entorno de computación diferente. Lo cual francamente creo que es fantástico y asombroso que funcione. E incluso, incluso hasta el día de hoy todavía me sorprende que esto sea posible.

## 1.1. Descargar la imagen del contenedor

Así que vamos a probar usando algunas imágenes Docker y Docker, cuando lo ejecutas en tu sistema, tiene un registro docker en tu computadora, o en este caso, en el codespace, que lleva un registro de todas las diferentes imágenes que han sido descargadas y usadas en el pasado, y las diferentes capas de las cuales están construidas.

Podemos ver qué imágenes tenemos localmente con Docker haciendo "docker image ls". Y en este caso puedes ver que hay un montón de imágenes Docker aquí, que todas tienen que ver con configurar este Codespaces. Todas tienen que ver con contenedores de desarrollo y cosas así. Así que no necesitas preocuparte mucho por ellas, pero a medida que agreguemos más imágenes y las descarguemos, a medida que avanza este curso, puedes verificar esa lista y verás que el registro local está llevando un registro de todas estas cosas que hemos descargado.

Pero vamos a agarrar una nueva haciendo "docker pull". Y eso le dice a Docker que busque una nueva imagen desde la web.

Luego ponemos el URI para ese contenedor. Ahora esto podría ser una imagen docker que has construido localmente y luego enviado a internet. podría ser una imagen que otra persona ha hecho. Hay muchas, muchas, muchas formas diferentes de hacer imágenes Docker, pero posiblemente una de las formas más simples es externalizar eso, y hacer que alguien más lo haga por ti.

Y lo que vamos a usar en este, tutorial es un servicio de Seqera llamado Seqera Containers.

Ahora, Seqera Containers es totalmente gratis, y usa una pieza de software de código abierto que desarrollamos llamado Wave, que fue construido para gestionar contenedores de manera complementaria a Nextflow. Y maneja muchos de los casos de uso comunes con los que nos encontramos lidiando con Nextflow.

Es muy común que el software que necesitamos esté empaquetado en Conda, en canales Bioconda, o conda-forge o canales más específicos de dominio. Y Wave y Seqera Containers es realmente bueno construyendo imágenes a partir de eso.

Así que puedo ir a esta interfaz web y vamos a jugar con el paquete llamado "cowpy". Así que escribo el nombre del paquete que quiero. Busca, lo ha encontrado en el índice de paquetes Python, así que puedo usar eso. O si espero un poco más, está buscando bioconda y conda-forge. Y puedes ver, puedo especificar cualquier canal conda aquí. Así que si quieres encontrar un canal Nvidia o cualquier otra cosa, eso también debería funcionar.

Y luego puedo especificar si quiero que construya una imagen docker para mí o una imagen singularity y también qué arquitectura de CPU quiero usar. Así que amd64 o arm64.

Y una vez que los resultados de bioconda están listados, ahora puedo ver todas las diferentes versiones que están disponibles también. Voy a agregar eso. Y ahora podría seguir buscando y obtener más paquetes de Conda si quiero y componer este contenedor como quiera, pero solo quiero ese. Así que voy a hacer clic en Get Container.

Ahora, alguien más ya ha solicitado el mismo contenedor antes y se devuelve desde un registro, así que lo obtenemos inmediatamente. Pero si nadie más hubiera pedido este paquete de software o esta combinación de paquetes de software, Wave y Seqera Containers lo construiría sobre la marcha para nosotros.

Podemos copiar esta URL y también podemos ver los detalles de construcción. Y esto nos muestra qué hizo el servicio en el backend. Creó un archivo de entorno conda. Un archivo docker, y luego esto es, ejecutando el proceso de construcción docker. También ejecutó un escaneo, un escaneo de seguridad, así que puedes ver cualquier CVE. Y te dice cuándo se creó esto.

Wave y Seqera Containers pueden hacer mucho más que esto, pero este es un caso de uso simple, que es el más común. Y debería decir que estas imágenes se alojan durante al menos cinco años. Así que puedes construir estas URLs en tus pipelines y saber que no van a desaparecer pronto.

Así que tengo mi URL para mi imagen docker para cowpy.

Ahora puedo hacer "docker pull" de esa URL, y buscará todas las diferentes capas y descargará esta imagen para que esté disponible localmente para mí.

## 1.2. Usar el contenedor para ejecutar cowpy como un comando único

Bien, ahora intentemos realmente usarlo. Así que ahora voy a usar un comando "docker run" en lugar de "docker pull", y voy a usar la bandera "--rm", que solo le dice a Docker que apague este contenedor una vez que haya terminado de hacer lo que le he pedido. Y luego pongo el identificador del contenedor, que es solo un URI.

Y luego al final, especifico el comando que quiero que Docker ejecute dentro del contenedor generado desde esta imagen. Solo voy a decir cowpy, que es el nombre de la herramienta que está instalada desde Conda Forge, que está disponible dentro de la imagen.

Voy a presionar enter y ahí lo tienes. Hemos ejecutado cowpy en un sistema. Tenemos una pequeña vaca dándonos algo de información.

Ahora nota que cowpy no está instalado en mi sistema local. Así que si lo ejecuto solo sin todas las cosas de Docker, dice, comando no encontrado. Así que esto ha descargado una imagen. Ha creado un contenedor usando Docker, y luego ha ido a ese contenedor y ha ejecutado este comando para nosotros y nos ha devuelto la salida a nuestra terminal. Muy, muy genial.

## 1.3. Usar el contenedor para ejecutar cowpy interactivamente

Bien, vamos a ir un paso más allá ahora y ejecutar este contenedor interactivamente y echar un vistazo, para que podamos ver qué está pasando dentro del contenedor.

Así que si vuelvo atrás y tomo mi comando run y voy a deshacerme de cowpy al final ahí, porque realmente no quiero ejecutar cowpy. Quiero ejecutar una terminal Bash.

Y luego voy a volver aquí y voy a hacer "-it", que significa Interactivo y Terminal o TTY, y voy a presionar enter.

Y ahora puedes ver que el prompt, la parte antes de que escriba, ha cambiado. Este era el prompt de Codespaces donde decía el directorio, y ahora dice base y root y tmp. Así que ahora estoy dentro del contenedor, y si hago "ls", verás que los archivos que veo en este directorio son diferentes a los archivos que tengo en mi espacio de trabajo.

Y de hecho, no puedo ver ninguno de los archivos de mi espacio de trabajo local de codespaces o mi unidad local dentro del contenedor Docker. El tiempo de ejecución del contenedor docker, está completamente aislado y no puede escribir o leer ningún archivo de un sistema de archivos host externo.

Sin embargo, puedo ver el software que está instalado dentro del contenedor y ejecutarlo. Así que puedo ejecutar cowpy y podemos ver un poco más sobre cómo usar cowpy. Aquí puedo hacer "cowpy 'Hello World'" y eso le dice, le dice que realmente ponga mi cita dentro de una pequeña burbuja de diálogo. Y también puedes ejecutar diferentes tipos de vacas, así que no tiene que ser una vaca. Puedes hacer un "-c". Y estoy en Suecia, así que voy a elegir un alce. Muy bonito. Le di algunos cuernos.

Y hay un montón de diferentes que puedes explorar, que puedes ver descritos en los documentos de capacitación.

## 1.3.4. Montar datos en el contenedor

Bien. Sería bueno si pudiéramos ejecutar cowpy en los archivos de nuestro sistema de archivos.

Por supuesto, no es muy útil solo tener el contenedor y sin acceso a nada en absoluto. Podría ser seguro y reproducible, pero no es muy útil.

Entonces, ¿cómo hacemos eso? Voy a salir de este contenedor Docker escribiendo exit, y puedes ver que el prompt nos dice que ahora estamos de vuelta en nuestro Codespaces regular de nuevo.

Y voy a ejecutar el mismo comando de nuevo. Pero esta vez voy a agregar algunas banderas adicionales aquí. Y la importante es "-v", que significa montar un volumen, que es básicamente como una parte, parte de un espacio de disco.

El "-v" toma dos partes: hay como una cadena y luego dos puntos y una cadena. Y la primera parte es el sistema de archivos local, que debería montarse en el contenedor. Y luego la segunda parte es dónde eso debería terminar dentro del contenedor.

Ahora solo quiero cargar todo mi sistema de archivos local aquí. Así que "." es el directorio de trabajo actual. Así que solo voy a hacer "." y luego ":", y luego vamos a poner eso en un nuevo directorio dentro del contenedor llamado "my_project". Esto realmente podría llamarse cualquier cosa.

Y luego voy a ejecutar de nuevo.

En el directorio de trabajo donde me deja, que es /tmp, los archivos no están ahí. Pero si hago "ls my_project", ahí lo tenemos: todos los mismos archivos que teníamos localmente en nuestros Codespaces ahora están disponibles dentro del contenedor en esa ruta.

Este es acceso de lectura y escritura, así que puedo crear nuevos archivos en este directorio y aparecerán en mi sistema de archivos host. Este directorio particular, entonces se comporta exactamente como si estuviera fuera del contenedor, así que ahora puedo leer y escribir y hacer cosas.

## 1.3.5. Usar los datos montados

Bien, vamos a probar que podemos hacer esto. Hago "cat /my_project/data/greetings.csv". Si recuerdas, el contenido de este archivo se ve así. Ahora puedo canalizar eso a cowpy y la vaca imprimirá las diferentes salidas de ese archivo en su pequeña burbuja de diálogo, lo cual es algo divertido.

Así que puedes ver, ahora podemos usar el software en el contenedor para interactuar con los archivos en nuestro sistema host.

Bien, salgamos de vuelta y continuaremos con el resto del material de capacitación.

## 2. Usar contenedores en Nextflow

Así que eso es realmente genial usar contenedores. Espero que tenga sentido. Y puedes ver el valor de estos contenedores y por qué eso es útil para ejecutar software de análisis.

Pero, ¿cómo hacemos todo este mismo proceso dentro de Nextflow? No queremos estar ejecutando montones de comandos Docker nosotros mismos. Queremos simplemente dejar que Nextflow maneje todo eso por nosotros.

Así que trabajemos en esto. Vamos a agregar un nuevo proceso a nuestro pipeline, para ejecutar cowpy. Bien, así que creemos un nuevo módulo para nuestro nuevo proceso. Así que vamos a módulos, llamémoslo cowPy.nf, y luego voy a copiar el código del material de capacitación aquí.

Pero puedes ver que el proceso es muy simple. Se parece mucho a los que hemos hecho hasta ahora, tenemos un bloque de entrada con un path, que es nuestro archivo de entrada, y también un valor aquí para que esto sea un carácter, así que podríamos usar un alce de nuevo si queremos.

Y luego una salida, que es un solo archivo aquí, un path y luego un script. Y estamos haciendo lo mismo que hicimos interactivamente dentro del contenedor: estamos haciendo "cat" para leer el archivo de entrada. Estamos canalizando ese contenido a cowpy. Estamos eligiendo un carácter específico basado en esa entrada, estamos escribiendo a un archivo de salida llamado cowpy input file, que luego se hace eco a la salida.

Genial. Incluyamos eso. Así que include \{ cowpy \} from "./modules/cowpy.nf", ¿lo llamé cowpy? Sí.

Y luego llamemos a nuestro nuevo proceso aquí abajo en el bloque principal del workflow. Así que ejecutemos cowpy. Y tomaremos nuestro nuevo proceso cowpy y vamos a decir collectGreetings.out.

Y luego si recuerdas, había dos salidas para este módulo. Una llamada outfile y una llamada report. La extensión VS Code está auto-sugiriendo estas para nosotros y queremos el .outfile.

Siempre puedes saltar a este proceso aquí. O pasa el cursor sobre él y debería mostrarte rápidamente cuáles eran las salidas. Y también podemos hacer command clic en él y abrirá el archivo del módulo si quieres ver con más detalle.

Así que aquí vamos. Ese es el outfile ahí, y ese es el path. Así que ahora ese será el archivo de entrada para nuestro proceso cowpy. Fantástico.

Ahora si recuerdas, un proceso cowpy tiene dos entradas. También teníamos el canal de valor para el carácter. Así que podemos agregar "params.character" aquí. Podría haber codificado esto si hubiera querido, pero hagámoslo una opción CLI para que podamos hacer dash, dash character.

Bien. Ahora necesito definir el parámetro de entrada que acabamos de llamar y darle un valor predeterminado. Así que character, String. Y me gusta el alce, así que voy a establecerlo en moose por defecto.

Bien, intentemos ejecutarlo. Así que si hago Nextflow run hello containers, veremos qué pasa.

Podría haber usado dash resume si tuviera los viejos directorios de trabajo dando vueltas. Y nuevamente, estos primeros procesos habrían sido, almacenados en caché y habría sido un poco más rápido, pero debería ser básicamente lo mismo.

Ahora podemos ver de inmediato que ha lanzado un error cuando llegó a nuestro nuevo proceso, nos está diciendo aquí que hubo un error ejecutando el proceso cowpy y salió con un estado de salida 127. Este es el comando que intentó ejecutar. Se ve bien, se ve como esperábamos. Está tomando ese nombre de archivo de salida, que se ve correcto, lo está ejecutando con un carácter moose y tratando de guardarlo.

Pero puedes ver que el error del comando aquí está diciendo que el comando cowpy no se encontró. Y eso tiene sentido porque en realidad no le hemos dicho a Nextflow que use un contenedor todavía. Solo le hemos dado el comando cowpy. Y como dije antes, cowpy no está instalado en nuestro sistema local. Así que cuando intentó ejecutarlo, falló.

## 2.3.1. Especificar un contenedor para cowpy

Necesitamos decirle a Nextflow que hay un contenedor disponible y que puede usarlo. Entonces, ¿cómo hacemos eso?

Si entramos en nuestro módulo aquí, vamos a agregar una nueva declaración en la parte superior llamada "container". Y luego vamos a establecer eso en una cadena.

Ahora, si recuerdas, allá en Seqera Containers, puedo copiar esa URL y simplemente la coloco entre comillas aquí.

Ahora regresa e intenta ejecutarlo de nuevo.

Déjame ver si funciona esta vez.

Desafortunadamente, falla exactamente de la misma manera, a pesar de que ahora hemos definido un contenedor para que se ejecute el proceso. Así que para usar nuestra imagen docker, necesitamos decirle a Nextflow que habilite el uso de Docker cuando ejecutemos el workflow.

Y vamos a hacer eso creando un nuevo archivo de configuración. Así que voy a decir touch nextflow.config.

Este es un nombre de archivo especial donde si está en el directorio de trabajo mientras lanzo el pipeline, se cargará automáticamente. Así que si voy a este archivo Nextflow dot config, puedes ver que en realidad ya existe, lo cual había olvidado. Y tenemos docker.enabled aquí ya, pero está establecido en false, que es el valor predeterminado.

Así que si cambio eso a equals True en su lugar, docker.enabled. Y hay documentos de referencia para todos estos ámbitos de configuración en los documentos de Nextflow. Y también puedes ver que cuando paso el cursor con una extensión VS Code, trae los documentos específicos para esto y me dice qué significa y cómo configurarlo.

Así que ahora lo hemos establecido en true, y si ejecuto Nextflow de nuevo, Nextflow ahora sabrá buscar esa imagen docker para nosotros si aún no la tenemos localmente, y luego ejecutar ese proceso con ese entorno de contenedor.

Y así podemos ver que se ha ejecutado exitosamente y tenemos un pequeño tic junto a cowpy. Fantástico. Si subo y miro en el directorio de resultados, el archivo no está ahí todavía. Y eso es porque todavía necesitamos, publicar este archivo de salida igual que todos los demás.

Así que vamos al bloque publish dentro del workflow, digamos mycowpy equals cowpy.out.

Y luego aquí abajo en el bloque output, mycowpy, llaves onduladas path. Ups. Hello containers. Mode, copy.

Si ejecuto de nuevo ahora, debería ejecutarse exactamente de la misma manera. Podría haber usado dash resume y olvido cada vez. Y luego subo y ahora tenemos un nuevo archivo creado llamado cowpy-COLLECTED, y ahí está mi alce diciendo BONJOUR, HELLO, HOLà Fantástico.

Ahora, por supuesto, también podría pasar ahora "--character". ¿Cuáles son las diferentes opciones? ¿Creo que hay un pavo? Así que puedo usar character Turkey. Va a ejecutarse exactamente de la misma manera. Perdí otra oportunidad de usar dash resume, y ahora si cargamos nuestro archivo y ahora tenemos un pavo. Fantástico.

## 2.3.4. Inspeccionar cómo Nextflow lanzó la tarea en contenedor

Bien. Última pequeña cosa. Ejecutemos rápidamente este comando de nuevo, resume esta vez, y echemos un vistazo rápido en el directorio de trabajo para ver qué es lo que Nextflow está haciendo bajo el capó para hacer que todo esto funcione para nosotros.

Esta vez es súper rápido, vamos a ese directorio de trabajo, cd work/. Ahora si recuerdas tenemos un montón de archivos punto aquí y el que nos interesa en este caso es el que dije que casi nunca necesitamos mirar, llamado .command.run.

Si hago code dot command run, lo va a abrir en el editor. Y puedo buscar en este archivo y si me desplazo hacia abajo debería ver Docker run. Y puedes ver que Nextflow está haciendo el comando docker run para nosotros, cuando Docker está habilitado en una configuración. Tiene un montón de diferentes, banderas y cosas aquí, pero puedes ver la bandera "-v" que usamos nosotros mismos cuando estábamos ejecutando. Y puedes ver que está montando el directorio del espacio de trabajo local en el contenedor, para que el contenedor pueda acceder a nuestros archivos de entrada y guardar las salidas. Y luego al final, también está ejecutando .command.sh, que es el script generado, que tiene el comando cowpy dentro.

Y así puedes ver que Nextflow está tomando la lógica del workflow, que es lo que realmente nos importa, que es específica para nuestro análisis, y está haciendo todo lo inteligente detrás de escena para hacer que Docker funcione en nuestro sistema.

Y lo está haciendo de una manera realmente portátil para que un usuario final del pipeline pueda cambiar la tecnología que está usando: Docker, Singularity, Apptainer, Conda. Eso realmente no importa para la lógica del pipeline, pero Nextflow manejará todas las necesidades de infraestructura subyacentes, para que se ejecute en cualquier lugar.

Y ese es realmente el superpoder de Nextflow. Es reproducibilidad y portabilidad. Y con Nextflow puedes realmente compartir tu workflow y otras personas pueden ejecutarlo en sus sistemas y simplemente funcionará.

Esa es una cosa realmente, realmente difícil de hacer, y ahora tú también sabes cómo hacerlo con tus workflows.

Bien, eso es todo para este capítulo. si bajas al final del curso, encontrarás un, un cuestionario de nuevo sobre algunos contenedores. Espero que todo haya tenido sentido. Es una forma realmente genial de trabajar con análisis. Y si eres nuevo en contenedores, espero haberte convencido de que es el camino a seguir, y nunca mirarás atrás.

Pero con eso, toma un pequeño descanso quizás, y te unes a mí en un par de minutos para revisar la parte seis final del Hello Nextflow, que es todo sobre configuración.

Muchas gracias.
