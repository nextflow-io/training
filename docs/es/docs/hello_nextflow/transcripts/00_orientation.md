# Orientación - Transcripción del Video

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Traducción asistida por IA - [más información y sugerencias](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

<div class="video-wrapper">
  <iframe width="560" height="315" src="https://www.youtube.com/embed/G3CV-FcV-rc?si=nyLvwhrSB2m1NPc5&amp;list=PLPZ8WHdZGxmXiHf8B26oB_fTfoKQdhlik" title="YouTube video player" frameborder="0" allow="accelerometer; autoplay; clipboard-write; encrypted-media; gyroscope; picture-in-picture; web-share" referrerpolicy="strict-origin-when-cross-origin" allowfullscreen></iframe>
</div>

!!!note "Nota importante"

    Esta página muestra únicamente la transcripción. Para instrucciones completas paso a paso, regrese al [material del curso](../00_orientation.md).

## Bienvenida

Hola, bienvenido a Hello Nextflow. Mi nombre es Phil Ewels. Soy Product Manager para Open Source en Seqera, y estoy encantado de estar aquí hoy para guiarlo a través de este primer curso de entrenamiento de Nextflow.

Vamos a repasar los conceptos básicos de Nextflow, explicando cómo escribir y ejecutar pipelines y configurarlos.

Y usted va a construir su propio pipeline simple de múltiples pasos. Cubriremos terminología como operadores y channel factories, y al final del curso, estará listo para comenzar a construir sus propios pipelines de bioinformática.

Si tiene alguna pregunta, por favor comuníquese en community.seqera.io. Tenemos una comunidad de Nextflow muy activa, hay una sección dedicada al entrenamiento, así que simplemente déjenos saber dónde está atascado y alguien podrá ayudar.

Bien. Empecemos.

## Sitio Web de Entrenamiento

Todo el material de entrenamiento para los cursos de Nextflow se encuentra en training.nextflow.io. Puede acceder en su navegador web. Así que ábralo ahora y podemos echar un vistazo.

Estaré ejecutando esto con la versión 2.1.1. Publicamos pequeñas actualizaciones y correcciones de vez en cuando, así que no se preocupe si es un poco diferente, pero si el material ha cambiado demasiado, siempre puede usar este selector de versión en la parte superior para elegir la versión exacta de los materiales de los que voy a hablar.

Si prefiere el modo claro, puede cambiar el tema del sitio web aquí.

Vea las traducciones aquí, aunque al momento de esta grabación, realmente solo está el inglés, que cubre este nuevo material.

Y también vea todo el código fuente del sitio web de entrenamiento y todo con lo que trabajaremos en GitHub.

La página de inicio aquí enumera todos los diferentes cursos de material de entrenamiento que tenemos. Así que desplazo hacia abajo, veremos Nextflow para principiantes con el curso Hello Nextflow que haremos aquí. Puede ver todos los demás cursos que también tenemos, que funcionan de manera similar.

## Configuración del Entorno

De hecho, voy a comenzar usando este primero en la parte superior, que es común para todos los cursos de entrenamiento, y específicamente se trata de configurar nuestro entorno.

Hago clic, me lleva a esta sección, y podemos ver instrucciones para desarrollar localmente. Si desea usar su propia computadora portátil con su propia copia de VS Code y sus propias instalaciones de software, o lo que esperamos que haga la mayoría de la gente, que es usar algo llamado GitHub Codespaces.

Codespaces es un servicio proporcionado por GitHub donde ejecutan un servidor web en la nube, al que puede conectarse. Ese servidor tiene VS Code instalado, donde puede ejecutarlo en su navegador web, o si lo prefiere, conectarlo a su instalación local de VS Code. Toda la computación, todos los archivos, toda la edición ocurre de forma remota, lo que significa que todo el software que necesita viene preinstalado y es el mismo para todos.

## Crear un GitHub Codespace

Para crear el codespace con todo lo que necesitamos, busque los botones en el material de documentación, que dicen "Open in GitHub Codespaces". Voy a hacer clic en eso ahora, abrirlo en una nueva pestaña. Y se me presenta esta página web. Ahora puede ver que esto está preconfigurado para establecerse con nextflow-io training.

Simplemente puedo hacer clic en crear nuevo codespace. Pero en realidad recomendamos que usemos una máquina un poco más grande para el entrenamiento de Nextflow con cuatro CPUs en lugar de dos. Puede cambiar qué versión del material usa. Así que esto está predeterminado a 2.1.1 porque esa es la versión de la documentación desde la que seguí el enlace. Pero también podría establecerlo en una rama específica del repositorio si lo deseo.

Ahora voy a hacer clic en crear codespace. Y va a comenzar a configurar el entorno para mí.

## Creación de Codespace

Ahora, la primera vez que haga esto, va a tomar bastante tiempo, así que ahora es un buen momento para ir y tomar una taza de té. Póngase cómodo, charle con la persona sentada a su lado.

Si está interesado, puede hacer clic en building codespace aquí abajo para ver los registros de la configuración. Y puede ver aquí que está descargando una imagen de Docker con todo lo que necesito y configurando el entorno.

Ahora, solo tiene que esperar así la primera vez que cree un codespace. Si va a github.com/codespaces aquí, verá todos los diferentes Codespaces que tiene abiertos. Aquí está el que acabo de crear. La próxima vez que haga esto, puede ir aquí y puede seleccionar el codespace anterior y simplemente volver directamente a él. Y es un proceso mucho, mucho más rápido para calentar ese entorno existente. Eso también mantendrá todos los cambios que haya realizado en VS Code y en los archivos, por lo que no perderá su progreso si se va y regresa.

Puede hacer clic en los tres puntos aquí para hacer otras acciones. Por ejemplo, si lo configuró con dos CPUs y ahora quiere cuatro, puede cambiar el tipo de máquina. O si quiere comenzar desde cero y fresco, puede eliminar el codespace.

## Introducción a VS Code

Bien, Codespaces ha terminado de configurar mi entorno y ahora se presenta con VS Code en el navegador web.

Si está acostumbrado a VS Code. Esto se sentirá muy familiar si no lo ha usado antes, es bastante simple. Hay algunas partes diferentes de la página de las que debe estar al tanto.

Aquí a la izquierda, tenemos la barra lateral. Puede ver el Explorador configurado con todos los diferentes archivos en el repositorio de GitHub del repositorio de entrenamiento.

En estos botones en la parte inferior izquierda, pueden haber diferentes herramientas. En la barra lateral. Puedo buscar todos los archivos en todo el proyecto. Puedo trabajar con Git, puedo trabajar con GitHub, todo tipo de cosas así.

En la parte superior aquí está el menú principal. El explorador de archivos es el que tendremos más aquí, y puede hacer clic derecho en cualquiera de estos archivos y hacer las cosas normales que esperaría. Es posible que deba hacer clic en algunas advertencias como esta donde puede cortar, copiar y también puede descargar a su máquina local.

Cuando se carga el codespace, nos da una vista previa del archivo markdown en esta área principal aquí. Este es el mismo que el que se renderiza en github.com. Puedo cerrar eso y si hago doble clic en ese archivo Readme, verá que lo abre como código en el editor de código y al igual que con cualquier otro archivo, podemos editar este código directamente.

Finalmente, aquí abajo, tenemos la ventana de terminal. Estaba viendo los registros mientras se construía, así que eso es lo que está mostrando actualmente. También puedo presionar este botón más para iniciar una nueva sesión de terminal. Esto no se está ejecutando en mi máquina. Recuerde, esto se está ejecutando en la nube, y si hago tree tres a profundidad de dos, verá todos los mismos archivos aquí, que estaban a la izquierda.

## Mostrar solo archivos "hello-nextflow"

Este repositorio de GitHub contiene todos los diferentes conjuntos de entrenamiento, no solo el que estamos haciendo. Así que si lo desea, puede enfocarse solo en la carpeta Hello Nextflow. Una forma de limpiar esto un poco es ir al menú archivo y luego agregar carpeta al espacio de trabajo.

Hacemos clic en eso, vamos a training. Hello nextflow, y hacemos clic en agregar. Actualizará su pantalla. Y luego en el Explorador, ahora tenemos dos espacios de trabajo diferentes, el que teníamos antes para training y uno con solo Hello Nextflow.

Si lo desea, puede hacer clic derecho en training y hacer clic en eliminar carpeta del espacio de trabajo para deshacerse de ella de la barra lateral por completo.

Ahora tenemos solo archivos para este curso de entrenamiento en particular en el lateral. Puedo ocultar esa advertencia y ahora puedo hacer lo mismo en el terminal aquí y hacer CD para cambiar directorio. Hello, Nextflow. Y nuevamente, tenemos los mismos archivos aquí, que están en la barra lateral.

## Hello Nextflow: archivos

Mirando estos archivos para el curso Hello Nextflow.

Tenemos un montón de archivos .nf, que son para Nextflow, y hay uno de estos archivos para cada uno de los capítulos del curso de entrenamiento. Trabajaremos a través de estos archivos y los modificaremos en los ejercicios.

También tenemos un archivo nextflow.config, que solo tiene configuraciones básicas de configuración para ejecutar Nextflow en este entorno, de las que realmente no necesita preocuparse en este punto. Un archivo greetings.csv, que usaremos para procesar datos, que se introducirá en la siguiente parte de este curso, y un archivo test-params.json, que se usará en la parte seis y puede ignorar por ahora.

Estos archivos de Nextflow son solo el comienzo de cada ejercicio. Si desea ver cómo deberían verse cuando estén terminados, puede ir a un directorio solutions y están las respuestas para cada parte del curso de entrenamiento, por lo que puede ver una versión funcional de lo que está apuntando.

## Abrir un terminal

Si en algún momento cierra el terminal y no puede recordar cómo volver, no se preocupe por eso. Estos botones en la parte superior derecha abren y cierran diferentes paneles en el espacio de trabajo. Así que haga clic en este para el panel inferior y reaparecerá. Y solo asegúrese de que tenga terminal seleccionado aquí. También puede hacer clic en este botón aquí, la flecha en el lado derecho de un terminal para ponerlo en pantalla completa.

Me verá haciendo eso bastante porque tengo VS Code ampliado para que pueda leer el texto. Dependiendo del tamaño de su pantalla, es posible que necesite o no hacer esto. Lo mismo ocurre con minimizar el panel lateral.

Bien. Eso es suficiente para el entorno. Creo que estamos listos para comenzar. Acompáñeme de vuelta en el siguiente video para el capítulo uno.

[Siguiente transcripción del video :octicons-arrow-right-24:](01_hello_world.md)
