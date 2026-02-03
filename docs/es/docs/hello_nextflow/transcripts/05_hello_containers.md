# Parte 5: Hola Contenedores - Transcripción

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Traducción asistida por IA - [más información y sugerencias](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

<div class="video-wrapper">
  <iframe width="560" height="315" src="https://www.youtube.com/embed/5PyOWjKnNmg?si=QinuAnFwFj-Z8CrO&amp;list=PLPZ8WHdZGxmXiHf8B26oB_fTfoKQdhlik" title="YouTube video player" frameborder="0" allow="accelerometer; autoplay; clipboard-write; encrypted-media; gyroscope; picture-in-picture; web-share" referrerpolicy="strict-origin-when-cross-origin" allowfullscreen></iframe>
</div>

!!!note "Notas importantes"

    Esta página muestra únicamente la transcripción. Para instrucciones completas paso a paso, regrese al [material del curso](../05_hello_containers.md).

    Los números de sección mostrados en la transcripción se proporcionan únicamente con fines indicativos y pueden no incluir todos los números de sección en los materiales.

## Bienvenida

Hola, bienvenido a la Parte Cinco del curso de entrenamiento Hello Nextflow.

Este capítulo se llama Hola Contenedores. Vamos a hablar sobre cómo Nextflow se integra con herramientas como Docker y Singularity para usar contenedores de software para proporcionar software a los usuarios de su pipeline.

Esto significa que cuando las personas ejecutan su pipeline, no tienen que ir e instalar todas las diferentes herramientas ellos mismos. Nextflow lo hará por ellos.

Los contenedores son una tecnología extremadamente poderosa y crucial en reproducibilidad y facilidad de uso. Vamos a comenzar haciendo una breve introducción a los contenedores mismos, ejecutando algunos comandos de docker manualmente, y luego tomaremos esos mismos contenedores y los pondremos en nuestro pipeline de Nextflow.

Muy bien. Empecemos.

Entonces, como antes, comencemos cargando el material de entrenamiento. Vaya a training.nextflow.io. Hello Nextflow, Capítulo Cinco, Hola Contenedores.

Voy a entrar en mi entorno de Codespaces y a la izquierda aquí vemos hello containers punto nf.

Como antes, este es el mismo script con el que terminamos el capítulo cuatro anterior, así que debería verse familiar.

Tenemos nuestros parámetros de línea de comandos para especificar el archivo de entrada y el nombre del lote. Estamos incluyendo nuestros tres módulos, y tenemos nuestro workflow donde ejecutamos los tres procesos.

## 0. Calentamiento: Ejecute hello-containers.nf

Siéntase libre de ejecutar este workflow nuevamente y verificar que está produciendo las salidas que espera. Por ahora, en realidad voy a cerrarlo y sumergirme en la terminal.

## 1. Use un contenedor 'manualmente'

Para comenzar este capítulo, vamos a hacer un poco de recapitulación sobre la tecnología de contenedores. Si está muy acostumbrado a docker o singularity u otras tecnologías de contenedores, entonces trate esto como un repaso, o siéntase libre de omitirlo completamente.

Nextflow soporta muchos tipos diferentes de tecnologías de contenedores. Eso incluye Docker, Singularity, Podman, Shifter, Charliecloud, y más.

En este entrenamiento, vamos a enfocarnos en Docker. Ese viene preinstalado en los code spaces y es una de las tecnologías de contenedores más populares, especialmente si está desarrollando en su propia computadora o en su propia laptop.

Si está trabajando en un entorno académico en un HPC compartido, puede que encuentre que Singularity está disponible y no Docker. Está bien. Todos los conceptos son exactamente los mismos. Algunos de los comandos manuales son diferentes, pero si entiende Docker, también entenderá singularity.

De hecho, Singularity también está instalado en el entorno de Code Spaces. Entonces si gusta, puede intentar hacer las mismas tareas usando Singularity en lugar de Docker.

Bien, entonces ¿qué es la tecnología de contenedores? La idea detrás de Docker es que puede obtener una imagen de una fuente remota. Descargarla a su máquina local y luego crear un contenedor basado en esa imagen.

Este contenedor en ejecución es un poco como una máquina virtual ejecutándose en su computadora. Está aislado de su entorno, y viene preempaquetado con un sistema operativo y un conjunto de software disponible.

## 1.1. Descargue la imagen del contenedor

La sintaxis que necesitamos para obtener una imagen preexistente es "docker pull". Entonces voy a escribir eso en mi terminal, pero ahora necesitamos una imagen con la cual jugar.

Puede construir imágenes usted mismo. Puede encontrarlas en registros públicos como Docker Hub o quay.io. Pero una forma realmente buena de obtener imágenes rápidamente es usando Seqera Containers.

Este es un servicio comunitario gratuito que hemos construido en 2024, que puede usar sin inicio de sesión ni nada.

Si va a seqera.io/containers o hace clic en containers en la parte superior aquí, se le presenta una interfaz de búsqueda y puede escribir el nombre de cualquier herramienta disponible en Conda o en el Python package Index.

Por defecto, busca en los canales de Bioconda y Conda Forge, pero puede prefijar cualquier canal de Conda. Estoy aquí si lo desea.

Por diversión, usemos cowpy. Voy a escribir cowpy. Me da resultados de Python Package Index y Conda Forge. Voy a hacer clic en eso para agregarlo a mi contenedor. Podría agregar múltiples paquetes aquí si quisiera. Selecciono Docker, selecciono linux/amd64, y hago clic en Get Container.

Esto construye la imagen para mí bajo demanda si no ha sido creada ya, y me da una URL que puedo copiar.

Si está interesado, puede hacer clic en view Build Details, y eso lo lleva a una página que muestra el archivo de entorno conda que fue usado y el registro completo de construcción para la construcción, junto con los resultados del escaneo de seguridad.

Si regreso a mis code spaces, ahora puedo pegar este nombre de contenedor y presionar enter.

Docker ahora descarga todas las diferentes capas dentro de esta imagen de contenedor, y ahora nos dice que esta imagen está disponible para usar.

## Descargando una imagen de Singularity

Si está usando singularity, el proceso es básicamente el mismo. Seleccionamos nuestros paquetes de imagen, seleccionamos cowpy. Ahora elegimos Singularity en lugar de Docker y hacemos clic en Get Container. Eso nos da una URL de imagen usando oras://. O si prefiere, puede usar https:// marcando esa casilla. Copie esa URL. Ahora vaya a Code Spaces. En realidad tenemos Apptainer instalado en este espacio, que es lo mismo que Singularity, pero están aliados entre sí. Entonces voy a hacer apptainer pull y luego voy a llamarlo cowpy sif, pero puede llamarlo como quiera. Pegue la URL. Y eso va a descargar esa imagen para mí.

Podría hacer ls -lh y ver cowpy.sif

Singularity es diferente a Docker, en que singularity almacena todas las imágenes en archivos planos, mientras que Docker tiene un registro donde mantiene todas las capas por separado en su máquina host, y tiene un demonio ejecutándose para hacer seguimiento de todo eso.

## 1.2. Use el contenedor para ejecutar cowpy como un comando único

Bien, volvamos a Docker. Ahora podemos intentar ejecutar esta imagen que creamos haciendo docker run.

Voy a hacer dash dash rm, que simplemente hace una ejecución única de la imagen. Y voy a pegar la URL de la imagen. Y luego finalmente, termina esto con un comando que desea ejecutar.

La imagen que generamos tenía cowpy instalado, así que intentemos cowpy.

Ahí está. Ejecutó nuestro comando. No tengo cowpy instalado localmente. Puede ver que si intento ejecutarlo, no existe. Sin embargo, en este comando, lo ejecuté usando Docker y correctamente generó esta salida.

## 1.3. Use el contenedor para ejecutar cowpy interactivamente

Podemos ir más allá de esto si queremos y activar un contenedor interactivamente y mirar alrededor adentro. Nuevamente, hago "docker run dash dash rm". Ahora voy a hacer dash it, que le dice a Docker que queremos una terminal interactiva. Hago la URL de la imagen nuevamente, y esta vez, en lugar de hacer cowpy, voy a hacer bin bash porque el comando que queremos ejecutar es bash.

Esto nos lleva dentro de este contenedor en ejecución y puede ver que el prompt ha cambiado ahora.

Si hago LS slash puede ver que los directorios aquí son diferentes.

Si abro una segunda terminal aquí a la derecha, que simplemente está ejecutándose en GitHub Code Spaces y hago LS slash, ve que tenemos directorios como workspaces y temp, mientras que aquí en Docker es diferente.

Entonces este entorno está completamente separado dentro de Docker y aislado de mi entorno host. Eso es algo bueno, porque eso aísla la ejecución de este comando en la imagen de Docker y lo mantiene reproducible entre diferentes personas en diferentes sistemas host.

Si desea usar datos de su sistema host dentro de la imagen de Docker, tiene que montar eso explícitamente en el contenedor.

Vamos a hacer eso en un segundo.

## 1.3.2. Ejecute el(los) comando(s) de la herramienta deseada

Primero, sin embargo, veamos si podemos ejecutar cowpy. Ahí nuevamente, el comando está disponible ahora directamente en la línea de comandos, y podemos comenzar a hacer cosas más complejas y pasar argumentos. Hello containers y en lugar de la vaca, hagamos el pingüino tux. Veamos qué más tenemos.

Hagamos cheese. Maravilloso. ¿Qué tal Dragon y Cow? Bastante bien.

## 1.3.3. Salga del contenedor

Bien. No puedo hacer mucho más porque no tengo ningún dato en este contenedor. Así que salgamos de esta imagen en ejecución y veamos si podemos montar algunos datos en el contenedor. Puedo hacer eso haciendo control D o escribiendo exit. Bien, ahora estoy de vuelta en mi code space regular de GitHub.

## 1.3.4. Monte datos en el contenedor

Para montar algunos datos en el contenedor de Docker, necesito usar dash V. Entonces voy a tomar mi comando docker anterior, ir al inicio hacer dash v. Voy a hacer "." para el directorio de trabajo local actual, y luego dos puntos para decir dónde debería montarse eso en el directorio host y hacer slash data. Entonces eso está montando este directorio particular en el contenedor en slash data.

Ahora si hago LS slash podemos ver que tenemos un nuevo directorio llamado data, y si hago LS data, puede ver todos los archivos que tenemos en la barra lateral aquí. Fantástico.

## 1.3.5. Use los datos montados

Ahora podemos comenzar a usar algunos de los archivos que están en el sistema host dentro de la imagen de Docker. Entonces puedo decir cat data greetings csv. Si recuerda, este es nuestro archivo CSV con nuestros diferentes saludos de antes, y puedo canalizarlo a cowpy. Fantástico. Ahora estamos llegando a algún lado.

Bien. Eso es suficiente para ejecutar Docker interactivamente. Esperamos que ahora tenga una idea de aproximadamente qué es Docker y cómo usarlo tanto para ejecutar un comando de manera única, como también para usar una imagen interactivamente. Si está usando singularity. Los comandos son todos muy similares excepto que hace cosas como apptainer exec o apptainer run, o singularity exec o singularity run.

## 2. Use contenedores en Nextflow

A continuación vamos a volver a nuestro workflow de Nextflow y ver cómo usar esta tecnología dentro del pipeline de Nextflow.

Cerremos la terminal y abramos Hello Containers nuevamente.

## 2.1. Escriba un módulo cowpy

Para continuar con nuestro ejemplo de cowpy, creemos un nuevo proceso en nuestro workflow, que use cowpy. Vayamos a modules, creemos un nuevo archivo y llamémoslo cowpy nf. Ahora voy a hacer trampa un poco y copiar el código del material de entrenamiento y presionar guardar. Y echemos un vistazo.

Entonces este es un proceso simple. Esperamos que ahora entienda cómo lucen los bloques de construcción de un proceso. Tenemos nuestro publishDir nuevamente, yendo a results. Tenemos dos entradas, un archivo de entrada y un string llamado character. Tenemos una salida cowpy input file, y tenemos un script que se ve exactamente igual a lo que ejecutamos manualmente dentro de nuestra imagen de docker hace un segundo: cat para imprimir un archivo, canalizándolo a cowpy, diciendo qué tipo de carácter de cowpy queremos usar, y enviando eso al archivo de salida, que pasamos como la salida aquí.

## 2.2. Agregue cowpy al workflow

Bien, volvamos a nuestro workflow, importemos este nuevo proceso. Entonces cowpy from modules cowpy nf. Creemos un nuevo parámetro para que podamos especificar qué carácter queríamos. Digamos Turkey por defecto. Y luego llamemos a este nuevo proceso al final del workflow,

cowpy. Y usemos la salida aquí de Collect Greetings. Entonces collect greetings out, out file aquí. Y luego necesitamos un segundo argumento, que es el nuevo params que acabamos de hacer. params dot character.

## 2.2.4. Ejecute el workflow para verificar que funciona

Bien, veamos si nuestro nuevo proceso funciona. Nextflow run hello containers. Esto debería ejecutar esos primeros tres procesos y luego intentar ejecutar cowpy al final.

Tuvimos un error. Lo que está diciendo aquí, cowpy tuvo un error y tuvo un estado de salida 127 y efectivamente, comando sh cowpy comando no encontrado.

No le dijimos a Nextflow que tenemos una imagen de Docker disponible para cowpy, así que intentó ejecutarlo en nuestro sistema host y no tenemos cowpy instalado en nuestro sistema host, así que desencadenó un error.

## 2.3. Use un contenedor para ejecutarlo

Entonces lo que necesitamos hacer es que necesitamos decirle a Nextflow que tenemos un contenedor disponible. Vayamos a nuestro proceso cowpy y vamos a agregar una nueva directiva en la parte superior del proceso llamada container.

Luego encontramos nuestra imagen, copiamos la URL, y la ponemos en un string.

Esto no es suficiente por sí mismo porque un pipeline de Nextflow puede tener varias formas de especificar software. También podría hacer conda conda-forge cowpy, por ejemplo. Y Nextflow necesita saber cuál de estas tecnologías desea usar.

## 2.3.2. Habilite el uso de Docker a través del archivo nextflow.config

Entonces para ejecutar con Docker habilitado, vamos a adelantarnos un poco y usar el archivo Nextflow config, que es algo que vamos a cubrir con más detalle en el próximo capítulo. Puede ver en este directorio que tenemos un archivo llamado Nextflow Config, y aquí ya tiene docker.enabled False.

Vamos a cambiar eso a True para habilitar Docker, y luego podemos intentar ejecutar el workflow nuevamente.

## 2.3.3. Ejecute el workflow con Docker habilitado

Nextflow run hello containers nf y esta vez cowpy se ejecutó exitosamente. Miremos en Results. cowpy collected test y ahí está nuestro Turkey. Maravilloso.

Entonces en segundo plano ahí, Nextflow sabía que tenía un contenedor disponible para ese proceso.

Obtuvo la imagen y ejecutó los comandos por nosotros.

## 2.3.4. Inspeccione cómo Nextflow lanzó la tarea contenerizada

Si tiene curiosidad, en realidad podemos ver exactamente qué hizo mirando en el directorio de trabajo. Si hago code work, y luego el hash y luego command run, que si recuerda es el archivo real que se ejecuta para esa tarea, podemos entrar y podemos buscar una función llamada NXF launch. Y aquí puede ver el comando docker exacto que Nextflow usó, que se parece mucho a lo que estábamos haciendo manualmente en la terminal antes. Docker run. Vinculando este directorio host en el contenedor, y luego especificando la URL del contenedor.

Entonces no hay magia aquí. Es solo que Nextflow está haciendo automáticamente el trabajo pesado por usted de una manera que significa que puede especificar fácilmente contenedores en su pipeline, que luego están fácilmente disponibles para cualquier otra persona que ejecute su workflow. Y esas personas ya no tienen que pensar en gestionar software para ejecutar su pipeline de análisis.

Muy, muy simple, muy conveniente, y también realmente reproducible. Bueno en general.

Bien, buen trabajo. Ese es el final del Capítulo Cinco. Únase a nosotros en el próximo video para la parte seis, que es la parte final de este entrenamiento Hello Nextflow, donde hablaremos sobre la configuración de Nextflow con más detalle.

Nos vemos en el próximo video.

[Transcripción del próximo video :octicons-arrow-right-24:](06_hello_config.md)
