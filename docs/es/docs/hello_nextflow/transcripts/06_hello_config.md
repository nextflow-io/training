# Parte 6: Hola Config - Transcripción

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Traducción asistida por IA - [más información y sugerencias](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

<div class="video-wrapper">
  <iframe width="560" height="315" src="https://www.youtube.com/embed/IuDO2HeKvXk?si=tnXTi6mRkITY0zW_&amp;list=PLPZ8WHdZGxmXiHf8B26oB_fTfoKQdhlik" title="YouTube video player" frameborder="0" allow="accelerometer; autoplay; clipboard-write; encrypted-media; gyroscope; picture-in-picture; web-share" referrerpolicy="strict-origin-when-cross-origin" allowfullscreen></iframe>
</div>

!!!note "Notas importantes"

    Esta página muestra solo la transcripción. Para instrucciones completas paso a paso, regrese al [material del curso](../06_hello_config.md).

    Los números de sección mostrados en la transcripción se proporcionan solo con fines indicativos y pueden no incluir todos los números de sección en los materiales.

## Bienvenida

Hola, bienvenido a la parte seis del curso de entrenamiento Hello Nextflow.

Este capítulo se llama Hola Config, y es la parte final de nuestro curso de entrenamiento.

En este capítulo, vamos a hablar sobre la configuración de Nextflow. La configuración de Nextflow es realmente poderosa. Nos permite ejecutar el mismo pipeline en múltiples infraestructuras de cómputo diferentes con diferente aprovisionamiento de software y diferentes opciones en el pipeline mismo.

Esto significa que puede tomar pipelines de Nextflow construidos por otras personas y ejecutarlos en su sistema, aunque pueden haber sido construidos para una infraestructura completamente diferente. Esta capacidad de configurar Nextflow hace que los workflows sean verdaderamente portables y compartibles.

En este capítulo, usaremos el workflow que hemos construido en partes anteriores, pero no vamos a editar el código del workflow en absoluto. Solo vamos a mirar nuestro archivo de configuración de Nextflow y ver cómo cambiar la configuración altera la forma en que Nextflow se ejecuta.

Bien, comencemos.

Como antes, comencemos yendo a training.nextflow.io. Vaya a la izquierda en Hello Nextflow y capítulo seis, Hola config. Ahora voy a entrar en mi entorno de GitHub Codespaces y verificar el script que usaremos.

## 0. Calentamiento: Verificar que Docker esté habilitado y ejecutar el workflow Hello Config

Este se llama Hello Config, y está comenzando desde donde estábamos antes. Entonces se ve exactamente igual con nuestros tres parámetros: greetings para el archivo CSV, batch para el nombre de lote de salida y character para el nombre de cowpy. Tenemos nuestras cuatro importaciones de los diferentes procesos, y luego tenemos un workflow donde los encadenamos.

De hecho, voy a cerrar este archivo ahora porque no vamos a tocar el archivo de Nextflow en absoluto en este capítulo. Vamos a trabajar puramente dentro del archivo de configuración. Si miro el archivo nextflow.config que brevemente vimos en el capítulo cinco anterior, podemos ver que tenemos una sola declaración aquí: docker.enabled = true, que le está diciendo a Nextflow que use Docker cuando ejecute este workflow.

Estoy usando nextflow.config en la raíz del pipeline aquí, que se carga automáticamente cuando ejecuto Nextflow. Pero recuerde, Nextflow puede cargar archivos de configuración desde múltiples lugares.

Si verifico con la documentación de Nextflow, voy a Configuración, puede ver una lista de estos lugares y una prioridad en la que se cargan.

Bien. Verifiquemos que nuestro workflow se está ejecutando como esperamos. Abro una terminal. Hago nextflow run hello-config y presiono enter. Deberíamos tener esos cuatro procesos ejecutándose, terminando con un comando cowpy. Efectivamente, esto funcionó correctamente. Tenía Docker habilitado, descargó Docker y ejecutó cowpy para mí, tal como lo hizo al final del capítulo cinco.

## 1. Determinar qué tecnología de empaquetado de software usar

Bien. Digamos que estoy ejecutando en un HPC y no tengo Docker instalado. Lo mejor que podría hacer en este escenario sería usar Singularity o Apptainer. Si fuera a hacer eso, iría al módulo cowpy y cambiaría este contenedor para usar la imagen de singularity como mostré en el capítulo anterior, con un oras://, que también puede obtener de Seqera Containers.

Luego iría a nextflow.config, establecería docker.enabled en false y haría singularity.enabled = true. O, si uso Apptainer, apptainer.enabled = true y eso funcionaría.

Nextflow también admite otras tecnologías además de contenedores, algo con lo que podría estar familiarizado es conda. Aquí podemos hacer conda.enabled = true y establecer Docker en false. conda no usa la misma directiva de contenedor. En su lugar, podemos agregar una nueva aquí llamada conda. Luego especificamos el paquete conda que queremos usar. Es una buena práctica ser lo más específico posible para tratar de hacer el pipeline lo más reproducible posible. Así que voy a especificar el canal conda, conda-forge, y luego cowpy, y la versión exacta, que fue 1.1.5.

También podría simplemente escribir cowpy si quisiera, pero eso podría resolver a una versión diferente de cowpy en diferentes ejecuciones del pipeline.

Lo bueno de esto es que no he tocado la directiva docker en absoluto. Esta imagen de Docker todavía está ahí. Solo estoy proporcionando dos alternativas ahora, y estas pueden activarse o desactivarse usando solo un archivo de configuración.

## 1.3. Ejecutar el workflow para verificar que pueda usar Conda

Conda ahora está habilitado, así que probémoslo.

Genial. Se está ejecutando y puede ver que hay un mensaje de Nextflow aquí diciendo que Nextflow está creando un entorno conda para mí, y está usando esta ubicación de caché.

En segundo plano, Nextflow está ejecutando comandos "conda create" para mí para crear un nuevo entorno conda aislado con solo los paquetes que quiero, y luego instalando y obteniendo esos paquetes conda para que pueda ejecutar el proceso.

Puede ver que tomó un poco de tiempo allí porque estaba creando el entorno e instalando el software por primera vez. Sin embargo, ha almacenado en caché este entorno, así que si ejecuto el mismo comando de Nextflow nuevamente, debería ser mucho más rápido porque reutilizará el mismo entorno conda.

Una de las cosas geniales de esto es que estas directivas pueden especificarse a nivel de proceso, no solo para todo el workflow. Entonces, si lo desea, puede mezclar y combinar qué tecnología se usa para diferentes procesos.

## 2. Asignar recursos de cómputo con directivas de proceso

El archivo de configuración de Nextflow puede hacer mucho más que solo empaquetado de software. También podemos decirle a Nextflow cómo ejecutar realmente los pasos en el pipeline. Un ejemplo es decirle al sistema host qué recursos deben estar disponibles para cada tarea en ejecución.

Por defecto, Nextflow no proporciona mucho. Proporciona un solo CPU y solo dos gigabytes de memoria para cada proceso.

Esto es probablemente algo que querríamos cambiar, para que los procesos que toman mucho tiempo para ejecutarse puedan tener más recursos y ejecutarse más rápidamente, pero puede ser difícil saber qué asignar a un proceso. Nextflow tiene algunos buenos trucos bajo la manga para ayudarlo con esto.

## 2.1. Ejecutar el workflow para generar un reporte de utilización de recursos

Ejecutemos el workflow nuevamente. Esta vez, voy a agregar un argumento adicional, que es -with-report. Es una opción central de Nextflow, así que es un solo guion. Y luego cualquier nombre de archivo que me guste. En este caso, voy a llamarlo report-config-1.html.

Voy a ejecutar el workflow nuevamente. Va a ejecutarse exactamente como antes, pero me va a dar un reporte de ayuda adicional, que puede ver que ahora ha aparecido aquí en la barra lateral.

Voy a hacer clic derecho en este archivo, hacer clic en descargar, que lo descarga de GitHub Codespaces a mi sistema local, para que pueda verlo fácilmente en el navegador web aquí arriba.

Este reporte puede generarse para cualquier ejecución de Nextflow, y tiene mucha información. Comienza en la parte superior con algunos metadatos sobre qué comando se usó, cuándo se ejecutó el workflow, cuánto tiempo tomó, pero a medida que se desplaza hacia abajo, obtenemos información más detallada sobre los recursos que fueron utilizados por cada paso en el pipeline.

Debido a que cada proceso se ejecuta múltiples veces para diferentes tareas, tenemos un diagrama de caja que muestra la variación de los recursos que usamos para cada proceso.

Si me desplazo un poco más hacia abajo, veo información similar sobre la memoria utilizada y la duración del trabajo. También lectura y escritura de disco.

Puede imaginar que para un pipeline grande con tareas de larga ejecución, esto puede ser muy informativo sobre cómo ajustar la configuración de los recursos que está solicitando para que no solicite de más, pero también para que pueda proporcionar suficiente para que se ejecute rápidamente.

Si sigo desplazándome hacia abajo en el reporte, también vemos una tabla de tareas, que nos muestra información detallada sobre cada tarea individual que se ejecutó en el workflow. Esto incluye información como el script resuelto, que se ejecutó.

Bien, volvamos a nuestro archivo de configuración. Vimos que realmente no necesitábamos mucho para nuestro workflow, así que digámosle a Nextflow que solo necesitamos un gigabyte de memoria para cada proceso en el workflow.

Ahora cuando lo definimos así a nivel de proceso, esto se aplica a cada proceso individual en el pipeline.

## 2.3. Establecer asignaciones de recursos para un proceso individual

Por el bien del argumento, supongamos que cowpy está realmente haciendo mucho trabajo pesado y necesita más recursos que las otras tareas. Podemos definir un bloque adicional de configuración aquí, que se aplica solo a ese proceso usando withName: COWPY.

Esto se llama un selector de configuración, y podemos definir diferentes patrones aquí para coincidir con diferentes procesos. Por ejemplo, podría hacer COW\*. Luego sigo eso con algunas llaves y démosle dos gigabytes de memoria en lugar de uno y digamos dos CPUs.

Ahora Nextflow proporcionará a cada proceso en el workflow un gigabyte, aparte de esta solicitud, que es más específica. Así que la sobrescribe. Y solo para cualquier proceso que se llame COWPY, obtendrá dos gigabytes de memoria y dos CPUs.

Tenga en cuenta que Nextflow es inteligente sobre la utilización de recursos. Entonces, si comienza a poner estos números en valores más altos, verá que Nextflow comienza a poner en cola los envíos de trabajos uno tras otro, en lugar de ejecutarlos todos en paralelo, para que no solicite en exceso los recursos que están disponibles.

## 2.4. Ejecutar el workflow con la configuración modificada

Intentemos ejecutar un workflow nuevamente y guardemos un nuevo reporte esta vez.

Bien, podemos descargar este archivo y echarle un vistazo.

Sí, como era de esperar, se ve básicamente exactamente igual porque este es un workflow ficticio, que no está haciendo nada real. Pero puede imaginar cómo este enfoque iterativo de definir límites y hacer workflows de la vida real con este tipo de reportes le permite hacer un enfoque basado en evidencia para establecer la configuración apropiada y realmente aprovechar al máximo los recursos computacionales que tiene disponibles.

Puede comenzar a ser realmente inteligente sobre esto. Nextflow tiene una capacidad incorporada para reintentar fallas, y puede aprovecharla en su archivo de configuración usando un closure como este y estableciendo dinámicamente los recursos que se ponen a disposición. Entonces aquí le he dicho a Nextflow que multiplique esos dos gigabytes por el intento de reintento. Entonces, el segundo reintento obtendrá cuatro gigabytes, el tercer reintento obtendrá seis gigabytes y así sucesivamente. Esto está un poco más allá del alcance de este curso de entrenamiento, pero si está interesado, consulte la documentación de Nextflow, que tiene una buena sección sobre lógica de reintento dinámico.

## 2.5. Agregar límites de recursos

Ahora, una cosa que podría notar sobre esto es que este tipo de cosas puede hacer que sea bastante fácil ir accidentalmente más allá de los recursos disponibles en su sistema. Si solicita más recursos de los que están disponibles, Nextflow generará un error sobre su configuración y detendrá la ejecución. Para evitar eso, puede usar algo llamado límites de recursos.

Bajo el ámbito process, en nuestro workflow, podemos definir límites de recursos como este, que toma un array, y podemos especificar la memoria máxima, CPUs y tiempo que están disponibles en este sistema.

Establecer valores altos aquí no aumenta la cantidad de recursos que se solicitan. Todavía vamos a usar un gigabyte en nuestras solicitudes, pero significa que si alguna de estas solicitudes llega a 750, alcanzarán ese límite y no se solicitará más que eso, lo que significa que Nextflow continuará ejecutándose y no fallará debido a recursos no disponibles.

Entonces, esta es una buena salvaguarda para usar, especialmente si está usando lógica dinámica con su asignación de recursos.

La otra situación donde esto es realmente útil es si está usando pipelines que son públicos y no están controlados por usted. Pueden venir con valores predeterminados de configuración, y Nextflow automáticamente tomará el enfoque correcto de limitar cualquier solicitud de recursos para ejecutarse en su sistema.

Bien, genial. Hemos hablado sobre software. Hemos hablado sobre la asignación de recursos, y hemos descrito diferentes ámbitos de configuración, tanto para todos los procesos como para procesos específicos.

## 3. Usar un archivo de parámetros para almacenar parámetros del workflow

Bien, a continuación vamos a dirigir nuestra atención a los parámetros. Podemos definir parámetros en el archivo de configuración tal como lo hicimos antes en el script de Nextflow. Entonces params.greeting = 'hello' o usar el ámbito params y establecer foo = 'bar'.

Y eso es genial para establecer valores predeterminados para su workflow. Sin embargo, cuando está ejecutando pipelines, puede ser bueno especificar parámetros en un archivo JSON o YAML.

Usar un archivo como este es mucho mejor que especificar opciones de línea de comandos con --. Como cuando ejecuta un workflow, podría tener que especificar muchos parámetros y puede ser tedioso escribirlos todos en una sola CLI y propenso a errores. Además, es poco probable que recuerde todos los parámetros que usó, así que si lo codifica en un archivo, es más fácil lanzar el workflow nuevamente, usando los mismos parámetros en el futuro.

Tenemos un archivo de ejemplo aquí llamado test-params, y puede ver que esto especifica los tres parámetros que tenemos en nuestro workflow con tres valores diferentes. Personalmente, encuentro que YAML es más fácil de escribir que JSON. Entonces, solo para demostrar que funciona, voy a crear un nuevo archivo llamado test.yaml y copiar estos, deshacerme de las comillas. Y guardar.

Estos archivos JSON y YAML pueden ser más fáciles de escribir ya que son sintaxis más familiar. Pero tenga en cuenta que estos son solo para parámetros y solo toman sintaxis de clave-valor como esta.

## 3.1. Ejecutar el workflow usando un archivo de parámetros

Probémoslo. Hago el mismo comando que antes. Me deshago del reporte y voy a hacer -params-file test-params.yaml.

No, esta es una opción central de Nextflow, así que es un solo guion.

Bien. Ejecutó el workflow y usó los parámetros en ese archivo YAML en lugar de que yo los especifique todos en la línea de comandos. Puede parecer exagerado solo para este ejemplo simple, pero puede imaginar si tiene 10 o 20 parámetros diferentes, puede ser una molestia escribirlos manualmente, y esto es simplemente mucho más fácil de editar en un editor de código y conservar para fines de reproducibilidad.

## 3. Determinar qué executor(es) deben usarse para hacer el trabajo

Bien. Hemos hablado sobre empaquetado de software con Docker y conda. Hemos hablado sobre requisitos de recursos de proceso con CPUs y memoria. Y hablamos un poco sobre cómo especificar parámetros al ejecutar workflows.

Las partes finales de la configuración realmente son la ejecución, la infraestructura de cómputo subyacente misma, y esta es la verdadera joya de la corona de Nextflow: que podemos ejecutar estos mismos workflows a través de múltiples infraestructuras de cómputo diferentes.

De hecho, voy a cambiar al material de entrenamiento escrito por un segundo. Bajo esta parte del entrenamiento, podemos ver algunos ejemplos diferentes de cómo diferentes executors, en este caso, planificadores HPC, definen los requisitos de recursos necesarios para enviar un trabajo.

Entonces, para Slurm, tiene estos encabezados SBATCH, que definen --mem y el número de CPU. Si está usando PBS, tiene encabezados diferentes, y si usa Grid Engine, tiene encabezados diferentes nuevamente.

Puede imaginar que es aún más diferente si quiere ejecutar en la nube, ya sea AWS Batch, Google Cloud, Azure o más.

Cada una de estas infraestructuras de cómputo subyacentes se llama un executor y Nextflow sabe cómo hablar con todos estos diferentes executors para enviar trabajos con la sintaxis correcta.

La buena noticia es que no tiene que saber sobre esto. Todo lo que tiene que hacer es decirle a Nextflow qué executor usar.

## 3.1. Apuntar a un backend diferente

Volvemos a nuestro archivo de configuración y en process hacemos executor, y voy a escribir 'local'.

Local es en realidad el predeterminado, si no especifica ningún otro executor, local es lo que se usará, y eso simplemente significa su sistema host, donde sea que haya lanzado Nextflow.

En su lugar, podría especificar slurm. Y eso enviaría trabajos de Slurm, o podría decir awsbatch, y eso enviaría trabajos a AWS Batch.

Necesita alguna configuración adicional en algunos casos, por ejemplo, ejecutar en la nube necesitará ciertas credenciales, pero realmente este es el núcleo de ello, y puede ser tan simple como una o dos líneas de configuración para ejecutar su workflow en un entorno de cómputo completamente diferente.

Aunque estamos ejecutando en un sistema simple dentro de Codespaces, todavía puedo jugar con esto un poco y pretender que estamos ejecutando en Slurm. Si luego lanzo el workflow nuevamente, nextflow run hello-config. Fallará porque no podrá enviar trabajos a Slurm. Pero aún podemos ir a los directorios de trabajo y ver qué hizo Nextflow. Entonces, si vamos a este directorio de trabajo y miramos .command.run. Puede ver en la parte superior de este archivo, ahora tenemos estas líneas de encabezado sbatch, que intentaron especificar los recursos necesarios para el trabajo de Slurm.

## 4. Usar perfiles para seleccionar configuraciones preestablecidas

Bien, ya casi llegamos. La parte final de este capítulo es hablar sobre perfiles de configuración. Si está ejecutando su pipeline en varios sistemas diferentes, podría ser molesto tener todos estos diferentes archivos de configuración de Nextflow, que necesita especificar cada vez.

En su lugar, puede codificar agrupaciones de configuración dentro de su archivo nextflow.config, y activar y desactivar esos grupos usando una bandera de perfil. Veamos cómo se ve eso.

## 4.1. Crear perfiles para cambiar entre desarrollo local y ejecución en HPC

Vamos a crear dos perfiles en nuestro ejemplo aquí, uno para mi laptop y uno para un sistema HPC más pesado. Voy a hacer trampa un poco y simplemente copiar el código del material de entrenamiento y ponerlo aquí.

Tenemos un nuevo ámbito llamado profiles, y luego tenemos un nombre para cada perfil, que puede ser cualquier cosa. Y dentro de eso tenemos configuración, que se ve exactamente igual que la configuración de nivel superior que ya escribimos. Entonces, nuevamente, tenemos ámbito process, ámbito docker.

En el perfil llamado my_laptop, estoy diciendo que se ejecute usando el executor local, entonces en mi sistema host y que use Docker.

En el perfil university_hpc aquí estoy diciendo que use Slurm para enviar trabajos, que use conda en lugar de Docker, y estoy especificando diferentes límites de recursos, que pueden coincidir con el tamaño del sistema de nodos en el HPC que estoy usando.

Por defecto, ninguna de esta configuración se usará cuando ejecute Nextflow, tengo que especificar que quiero usar uno de estos perfiles.

## 4.2. Ejecutar el workflow con un perfil

Hagamos nextflow run hello-config. Y voy a hacer -profile, un solo guion porque es una opción central de Nextflow. Y luego el nombre que le di, que es my_laptop. Nextflow ahora debería usar el bloque de configuración que se especificó dentro de ese perfil de configuración, y aplicarlo cuando ejecute Nextflow. Si quisiera usar el otro bloque de configuración, solo tengo que cambiar ese nombre de perfil. Mucho más fácil de recordar. Mucho más fácil de usar.

## 4.3. Crear un perfil de prueba

Tenga en cuenta, los perfiles pueden tener cualquier tipo de configuración, así que no tiene que estar relacionado con su entorno de ejecución. Por ejemplo, creemos un nuevo perfil aquí, que tenga un conjunto de parámetros. Podemos cambiar esto a tux y cambiar a My Profile, y ahora cuando hagamos profile test, va a especificar estos parámetros, que sobrescribirán los parámetros que se especifican en el nivel superior del workflow.

Cuando ejecute Nextflow, puede encadenar múltiples perfiles y se aplicarán en secuencia.

## 4.4. Ejecutar el workflow localmente con el perfil de prueba

Entonces puedo tomar el comando anterior y hacer coma test. Eso aplicará la configuración my_laptop primero, y luego aplicará la configuración test. Si hay alguna superposición, entonces el perfil a la derecha sobrescribirá cualquier configuración en perfiles anteriores. Si presiono enter, veamos qué pasa.

Bien, tenemos un nuevo archivo de resultados aquí. Puede ver el My Profile, que especifiqué como una de las opciones. Y también podemos ver cowpy, my_profile, y efectivamente, ahí está tux. Entonces eso ha funcionado.

## Conclusión

¡Bien! Increíble. Eso es todo. Lo ha logrado hasta el final del curso. Obtiene un poco de confeti de celebración. Bien hecho por terminar este capítulo.

[Siguiente transcripción de video :octicons-arrow-right-24:](07_next_steps.md)
