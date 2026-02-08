# Parte 1: Hello World - Transcripción del Video

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Traducción asistida por IA - [más información y sugerencias](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

<div class="video-wrapper">
  <iframe width="560" height="315" src="https://www.youtube.com/embed/tOukLxWCHiA?si=F0t9LFYLjAWoyRXj&amp;list=PLPZ8WHdZGxmWKozQuzr27jyMGqp9kElVK&amp;cc_load_policy=1" title="YouTube video player" frameborder="0" allow="accelerometer; autoplay; clipboard-write; encrypted-media; gyroscope; picture-in-picture; web-share" referrerpolicy="strict-origin-when-cross-origin" allowfullscreen></iframe>
</div>

!!!note "Notas importantes"

    Esta página muestra únicamente la transcripción. Para instrucciones completas paso a paso, regrese al [material del curso](../01_hello_world.md).

    Los números de sección mostrados en la transcripción se proporcionan solo con fines indicativos y pueden no incluir todos los números de sección de los materiales.

## Bienvenida

Hola, y bienvenido/a de nuevo.

Ahora estás en la Parte Uno del curso "Hello Nextflow" llamada "Hello World". En este capítulo, vamos a empezar a construir comprensión de los conceptos más básicos de Nextflow.

Así que con suerte ya estás configurado en Codespaces o en algún lugar equivalente con VS Code ejecutándose, y tienes tu carpeta Hello Nextflow en el espacio de trabajo en el Explorador con todos estos diferentes archivos aquí.

Vamos a empezar haciendo algunas cosas muy básicas en la terminal usando Bash, y luego veremos si podemos hacer las mismas cosas dentro de Nextflow para que te familiarices con la sintaxis.

## 0. Calentamiento

Así que empecemos muy simple. Comencemos solo con "echo", para imprimir algo en la terminal. "Hello World". Presiono enter y eso va a la terminal. Hello World. Con suerte eso no es una sorpresa para nadie viendo este curso.

Bien, hagamos algo con esto. En lugar de solo imprimirlo en la terminal, escribámoslo en un archivo. Voy a presionar la flecha arriba de mi teclado, que recorre el historial de Bash, así que me da mi último comando, y voy a agregar al final de él, un pequeño símbolo de mayor que, que redirige la salida de este comando a un archivo, y voy a llamarlo output.txt.

Enter de nuevo, para ejecutar ese comando, nada en la terminal esta vez, pero podemos ver en el lado izquierdo, el nuevo archivo ha aparecido aquí, llamado output.txt.

Podemos ver eso en la terminal con algo como cat. Entonces cat output.txt y efectivamente dice "Hello World". También podemos hacer doble clic en él y se abre en el editor de código en VS Code.

## 1.1. Examinar el código

Muy bien. Te dije que era simple. ¿Qué sigue? Intentemos tomar este proceso y hacerlo de nuevo, pero esta vez, hagámoslo dentro de Nextflow.

Como dije, todos los diferentes capítulos en este curso comienzan con un script y este se llama, Hello World. Así que voy a encontrar Hello World. Lo previsualiza si hago clic una vez, voy a hacer doble clic para abrirlo en el editor aquí. Y voy a deshacerme rápidamente de la terminal.

Ahora este es un script muy simple, tan simple como puede ser. Solo tiene 22 líneas de largo, y hace básicamente lo mismo. De hecho. Algo de esto debería verse familiar. Es lo que acabamos de escribir. Podemos ver nuestro comando bash redirigiendo a un archivo allí.

Bien. ¿Qué más? También, en este archivo, podemos empezar a ver algunos de los conceptos centrales de Nextflow. Tenemos un **process** en rojo aquí y un **workflow**. Estas son palabras clave especiales y terminología especial en Nextflow.

## 1.1.1. La definición del proceso

Diferentes **process** dentro de un **workflow** envuelven diferentes unidades lógicas de tu workflow. Cada **process** hace una cosa.

Cuando lo ejecutamos, genera una tarea o múltiples tareas, que son pasos reales de ejecución de un pipeline. Todos los **process** son entonces orquestados dentro de un bloque **workflow**, que vemos en la parte inferior, y en este caso, solo ejecuta ese único **process**.

El nombre del **process** sigue esta palabra clave aquí, y esto puede ser básicamente cualquier cosa. Y luego los contenidos del **process** están dentro de estas llaves.

Solo hay realmente un requisito para un **process**, que es que incluye algún tipo de bloque script o exec. Este está en las triples comillas aquí, y este es el script bash que se escribe en el directorio de trabajo cuando ejecutamos el pipeline y es la cosa que realmente se ejecuta en tu computadora o servidor.

Esto es bash típicamente, pero también puedes poner un shebang diferente aquí en la parte superior, y podría ser un script Python o un script R. No importa. Lo que esté en este script será ejecutado.

Hay otra cosa que hemos agregado a este **process** aquí, que es la declaración output. Esto le dice a Nextflow que este **process** está esperando un archivo de salida llamado output.txt. Dice que es un path, así que debe ser manejado como un archivo, no digamos, si esto fuera val, diría que es como una variable o valor.

Nota que esto no está creando este archivo. En realidad no lo está generando. Eso lo hace el script aquí abajo. Solo le está diciendo a Nextflow que espere un archivo de salida con este nombre de archivo.

## 1.1.2. La definición del workflow

Bien. Y luego en la parte inferior tenemos un **workflow** aquí, y de nuevo, tenemos una declaración. Este se llama Main. Esta es la parte equivalente del **workflow** a un bloque script, si quieres. Es la parte del **workflow** que hace algo. Y en este caso, estamos diciendo, llamar al **process** llamado sayHello.

Normalmente, por supuesto, tu pipeline se verá mucho más complejo que esto. Probablemente tendrás más de un **process**, y usarás canales para orquestar el flujo de datos entre ellos. Vamos a llegar a eso en las próximas partes de este curso, pero por ahora, esto es suficiente. Este es un pipeline válido, que debería funcionar.

Incluso puedo hacer clic en preview DAG aquí en VS Code. El DAG es una representación de una estructura de flujo de datos en el pipeline, y podemos verlo renderizado al lado como un diagrama mermaid. En este caso es muy simple. Hay una caja, que es el **workflow** y un **process**, que se llama sayHello, pero eso podría verse más interesante a medida que avanzamos.

## 1.2. Ejecutar el workflow

Bien, intentemos ejecutar este **workflow** y veamos qué sucede.

Voy a abrir la terminal de nuevo en la parte inferior, limpiar la salida, y voy a escribir Nextflow Run. Y luego solo voy a escribir el nombre del script, que es hello-world.nf. Y voy a presionar enter.

Bien, tiene algunas cosas estándar en la parte superior, que nos dice que Nextflow se ejecutó y qué versión estaba ejecutándose y cuál era el nombre del script y todo eso.

Y realmente lo importante que estamos buscando aquí es _aquí_, que es un resumen de las diferentes tareas que se ejecutaron.

Si el tuyo se ve así con una pequeña marca verde, entonces bien hecho. Acabas de ejecutar tu primer pipeline. Fantástico.

Nos dice aquí el nombre del **process**, que se ejecutó, que se llamaba Say Hello, y nos dijo que se ejecutó una vez y que fue exitoso. Esto se actualiza a medida que avanzas, así que cuando estés ejecutando un pipeline más grande, verás el progreso representado aquí. Pero como esto es tan pequeño, se ejecuta básicamente de inmediato.

## 1.2.2. Encontrar la salida y los registros en el directorio work

Ahora cuando ejecutas un pipeline de Nextflow, cada uno de esos **process** se une, y cada **process**, como dije antes, puede generar tareas una o múltiples. Así que en este caso, teníamos una sola tarea de este **process**. Solo se ejecutó una vez y eso se hizo bajo este _hash_ de tarea.

Nextflow no trabaja con los archivos en tu directorio de trabajo directamente, crea una carpeta especial llamada work. Y si hago "ls", veremos que ha aparecido aquí: _work_, y dentro de aquí hay subdirectorios para cada tarea que se ejecuta. Y eso coincide con este hash. Así que puedes ver si voy a "ls work/c4", y luego está truncado, pero comienza 203, y ese es el directorio de trabajo, que fue creado por este **process** cuando ejecutamos el pipeline. Y puedes verlo en el lado también.

Cuando listo esos archivos, puedes ver que el archivo output.txt fue generado. Puedes verlo aquí también. Y hay un montón de archivos ocultos, que no se muestran con mi "ls" regular.

Si hago clic en output.txt, efectivamente, tenemos nuestra salida. Fantástico. Así que el pipeline funcionó.

Puede parecer mucho código repetitivo para ejecutar lo que era esencialmente un script bash de una línea, pero tendrá más sentido a medida que nuestros **process** se vuelvan más complicados. Y este directorio work con Nextflow y estos archivos, que se crean, es realmente la columna vertebral de lo que hace a Nextflow tan poderoso.

Cada tarea, cada elemento de un pipeline está aislado de todas las demás tareas. Es reproducible. No entran en conflicto entre sí, y todo puede ejecutarse en paralelo. En realidad es una forma muy agradable cuando te acostumbras a ella, debido a este aislamiento, que puedes entrar y ver exactamente qué sucedió para una sola tarea y depurar.

Echemos un vistazo rápido a estos otros archivos en el directorio work. De arriba a abajo, tenemos un archivo llamado _.command.begin_. Este está vacío. Es solo lo que se llama un archivo centinela, creado por Nextflow diciendo, bien, estoy comenzando la tarea. Nada interesante allí.

Luego está _.command.error_, _.command.log_ y _.command.out_. Todas estas son salidas del comando bash o este script que se ejecutó. Esto es error estándar. Esto es salida estándar, y este es la combinación de los dos como salieron. Así que obtienes el orden lógico.

Bien, esos también estaban todos vacíos para esto, así que no muy interesantes, pero las cosas se vuelven más interesantes cuando llegas a _.command.run_.

Este es típicamente un script muy largo. Y esto es lo que Nextflow realmente ejecuta. Si entras aquí, comenzarás a ver toda la lógica interna de Nextflow y ver qué está haciendo y cómo está ejecutando tu **process**. Esto dependerá de dónde estés ejecutando, si estamos ejecutando localmente o enviándolo como un trabajo a SLURM, en cuyo caso tendremos encabezados SLURM en la parte superior. Todas estas diferentes configuraciones.

Generalmente, realmente no necesitas mirar este archivo. Es autogenerado por Nextflow y no hay nada realmente particular de tu pipeline que esté en él. Pero eso es realmente el núcleo de lo que se está ejecutando.

El siguiente es mucho más interesante. _.command.sh_ es el script generado, que vino de tu **process**, y aquí puedes ver que Nextflow agregó el encabezado Bash, y luego ejecutó nuestro comando, que estaba en nuestro bloque script.

Y eso es todo lo que hace el archivo _.command.run_, simplemente ejecuta este archivo _.command.sh_.

Este es uno muy útil, que es el que generalmente terminas mirando más cuando estás tratando de depurar algo y verificar que la lógica de tu pipeline de Nextflow esté haciendo lo que esperas que haga.

Finalmente, tenemos un archivo llamado _.exitcode_, y este solo captura el código de salida de una tarea, que en este caso fue exitoso. Así que el código de salida fue cero.

Si algo sale mal, te quedas sin memoria o algo más y falla, entonces esto es muy útil para entender qué salió mal.

## 1.3. Ejecutar el workflow nuevamente

Una cosa más para entender sobre los directorios work es que si sigo ejecutando este pipeline repetidamente, entonces si hago _"nextflow run hello-world.nf"_, va a hacer exactamente lo mismo, pero esta vez tendrá un nuevo id de tarea. Puedes ver que este hash aquí es diferente, y ahora si miro en work, hay dos directorios hash. Y estos están, de nuevo, separados entre sí.

Así que cada vez que ejecutas un **workflow** de Nextflow, a menos que uses el resume, que usa la caché, tocaremos más tarde, va a reejecutar esos **process** en nuevos directorios work, que están separados unos de otros. No tendrás colisiones de nombres de archivos, no tendrás ningún problema como ese. Todo está aislado y limpio.

Y si entramos en este directorio, puedes ver todos los mismos archivos y el mismo _output.txt_, que ha sido recreado desde cero.

## 2. Publicar salidas

Bien, eso es genial para Nextflow por sí mismo, mientras está ejecutando tu pipeline para que todas las cosas estén separadas entre sí y limpias y puedan ser administradas.

Pero no es súper útil si eres una persona tratando de explorar tus resultados. Realmente no quieres estar buscando a través de miles y miles de diferentes directorios work tratando de encontrar tus archivos de resultados. Y realmente no se supone que lo hagas. Los directorios work no están destinados a ser el estado final de dónde se crean tus archivos.

Hacemos esto publicando nuestros archivos.

## 2.1.1. Declarar la salida del proceso sayHello

Entonces, si vuelvo a nuestro script, vamos a trabajar en nuestro bloque **workflow** aquí. Vamos a decirle qué archivos esperar, qué archivos nos importan, y luego vamos a crear un nuevo bloque debajo llamado el bloque output.

Esta es la nueva sintaxis, que vino con el analizador sintáctico y será por defecto en la versión 26.04 de Nextflow. Entonces, si has usado Nextflow un poco antes, esta es una de las cosas que es nueva.

Así que tenemos el bloque main, y a continuación voy a decir publish y voy a decirle a Nextflow qué esperar de la publicación. Vamos a llamarlo _first_output_, y vamos a llamarlo, _sayHello.out_.

Accidentalmente cometí un error tipográfico ahí, pero esta es una buena oportunidad para también señalar algunas de las características de la extensión Nextflow VS Code. Puedes ver que de inmediato me dio una pequeña línea ondulada roja debajo de esto diciendo que algo está mal. Y si paso el cursor sobre ella, me va a decir que esta variable no está definida. No sé qué es.

Es bastante obvio en este caso, cometí un error tipográfico. Quise escribir, sayHello, y luego la línea ondulada desaparece.

Ahora es púrpura. El analizador sintáctico de Nextflow sabe que esto es un **process** y cuando paso el cursor sobre él, me da una representación reducida de cómo se ve este **process**. Así que puedo ver muy rápidamente de un vistazo que no toma ninguna entrada y nos da esta salida. Entonces trabajar en VS Code con esta extensión te da mucha información contextual mientras escribes código.

Nota que podemos referirnos a la salida de este **process** con la sintaxis _.out_. Y por el momento podemos llamar a esto como queramos, es solo un nombre de variable arbitrario.

## 2.1.2. Agregar un bloque output: al script

Donde se vuelve importante es cuando hacemos nuestro nuevo bloque aquí, y esto está debajo del bloque **workflow** ahora, ya no estamos dentro de **workflow**. Llaves de nuevo. Y aquí es donde simplemente le decimos a Nextflow dónde poner todos los archivos, que son creados por el **workflow**.

Ahora voy a tomar este nombre de variable, que creé aquí, y voy a ponerlo allí y poner algunas llaves para esto. Y voy a decirle a Nextflow que use un path. Ups. Path, entre comillas. Y voy a usar punto. Eso simplemente le dice a Nextflow que ponga el archivo en la raíz del directorio results. Así que no en subdirectorios ni nada.

Intentemos ejecutar nuestro **workflow** de nuevo. Si hago _"nextflow run hello-world.nf"_, entonces con suerte debería verse básicamente exactamente igual. Nada ha cambiado realmente con Nextflow aquí. Está ejecutando las mismas cosas. Solo las está haciendo en directorios work de nuevo.

Pero ahora si hago _"ls results/"_, verás que hay un nuevo directorio aquí que ha sido creado llamado results, que es el directorio base por defecto para la publicación del **workflow**. Y dentro hay un archivo llamado _output.txt_.

Si hago _"ls -l results"_, verás que esto en realidad es un enlace simbólico al directorio work. Así que este no es un archivo real, está enlazado al directorio work y ha recopilado todos los archivos allí para nosotros.

## 2.2. Establecer una ubicación personalizada

"Results" es el nombre por defecto para esta ruta. Si ejecuto el **workflow** de nuevo, y esta vez hago _dash_ guion simple, esto es, porque es una opción central de Nextflow. _" -output-dir **my** results"._ También podría simplemente hacer _"-o"_ para abreviar. Entonces va a establecer un directorio base diferente para dónde se almacenan los archivos y una vez más, aquí arriba en _myresults/_, ahora tenemos un _output.txt_.

Eso es genial, pero probablemente no queremos todos los archivos solo en la raíz. Queremos algo de organización, así que también podemos crear un subdirectorio aquí llamado como queramos. Digamos _"path 'hello_world'"_, y simplemente ejecuto esto de nuevo. _"nextflow run hello-world.nf"_. Debería ir al directorio results en un subdirectorio y efectivamente, ahora bajo results aquí en la parte superior tenemos _hello_world/_ y tenemos _output.txt_.

Cosa importante a notar, el viejo archivo _output.txt_ todavía está allí. El directorio results no se borra cuando haces esto. Solo los nuevos archivos se copian allí. Sobrescribirán archivos que ya están allí si tienen el mismo nombre de archivo, pero no limpiarán los viejos. Así que necesitas tener un poco de cuidado sobre cuándo reejecutas pipelines. Si no quieres que estén encima de los archivos que ya están allí. Asegúrate de usar un directorio vacío y limpio.

## 2.3. Establecer el modo de publicación para copiar

Bien, mencioné que estos archivos son enlaces simbólicos, así que si hago _"ls -l results/hello_world/"_, puedes ver que está enlazando simbólicamente al directorio work. Eso generalmente es algo bueno si estás trabajando en algo como HPC, y estos son archivos realmente enormes y no quieres duplicarlos, porque significa que los archivos solo se almacenan una vez en el sistema de archivos.

Sin embargo, significa que si eliminas el directorio work: si hago _"rm -r work"_ y limpio todos esos archivos intermedios que se crearon. Ahora, si intento leer este archivo _"results/hello_world/"_. Va a estar apuntando como un enlace simbólico a un archivo que ya no existe y los datos se han ido para siempre y son irrecuperables, lo cual tal vez no sea genial.

Así que generalmente decimos que es una buena práctica copiar los archivos en lugar de enlazarlos simbólicamente si puedes, porque es más seguro. Solo ten en cuenta que usará el doble de espacio en disco a menos que elimines esos directorios work.

Para hacer eso con el bloque output, voy a ir al primer output aquí. Establecí la ruta antes y ahora voy a establecer el mode y puedes ver mientras escribo, la extensión VS code está, sugiriendo cosas que sabe que es una directiva de salida aquí. Y voy a decir copy. Presiono guardar.

Reejecutemos el **workflow**. Va a crear los archivos de nuevo, nuevo directorio work.

Ahora, si voy a _"ls -l results/hello_world/"_ puedes ver que este es un archivo real y ya no es un enlace simbólico, y Nextflow lo copió. Bueno saberlo. Así que path y mode son cosas que te encontrarás escribiendo bastante.

Ahora, por supuesto, esto es muy simple. Vamos a hacer esto más complejo y poderoso a medida que avanzamos, y verás cómo hacer estas cosas dinámicas y no demasiado verbosas.

## 2.4. Nota sobre las directivas publishDir a nivel de proceso

Ahora, dije cuando comenzamos esto, que esta es una forma de sintaxis bastante nueva. Solo está disponible en las últimas versiones de Nextflow mientras grabo esto, y se llama Workflow Outputs.

Si usas esto, es genial. Desbloquea muchas otras características geniales dentro de Nextflow, como, Nextflow Lineage para ayudar a rastrear el origen de estos archivos a medida que se crean, y pronto será el predeterminado en 26.04. Y en una fecha posterior en el futuro, esta será la única forma de escribir tus **workflows**.

Sin embargo, como estamos en esta fase de transición ahora mismo, bien podrías ver pipelines en la práctica, que usan algo llamado publishDir, que es la forma antigua de hacerlo, y esto se define no a nivel de **workflow** y output, sino que se define a nivel de **process**.

Y esta declaración dice básicamente lo mismo. Dice, publica los archivos de resultados en un directorio llamado results, y usa un mode copy. Así que puedes ver que la sintaxis es muy similar. Pero cuando estés escribiendo nuevos pipelines ahora, intenta no usar esta directiva publishDir, incluso si la ves, en resultados de IA o en documentación u otros pipelines, porque esa es la forma antigua de hacerlo.

En 2026 todos deberíamos estar usando workflow outputs.

Todo esto está documentado, si estás haciendo esto y has usado Nextflow antes, puedes ir a los documentos de Nextflow aquí, nextflow.io/docs/. Y si desplazo hacia abajo a tutoriales, hay un tutorial llamado _Migrating to Workflow Outputs_.

Es realmente bueno. Recorre toda la sintaxis, cómo es equivalente a la sintaxis antigua, por qué la cambiamos, y, tiene una línea de tiempo y todo. Y recorre todos los diferentes escenarios con montones y montones de ejemplos. Así que puedes convertir fácilmente el código Nextflow existente a la nueva sintaxis.

## 3.1. Cambiar el proceso sayHello para esperar una entrada variable

Bien, así que tenemos nuestro script simple, que está ejecutando un **process**, creando un archivo, diciéndole a Nextflow que es una salida, y luego le estamos diciendo a Nextflow dónde guardar ese archivo. Ese es un buen comienzo.

Pero sería más interesante si no estuviera todo codificado. Así que a continuación, pensemos en cómo decirle a Nextflow que este **process** puede tomar una entrada variable, que es algo que podemos controlar en tiempo de ejecución cuando lanzamos un **workflow**.

Necesitamos hacer algunas cosas diferentes para que esto suceda.

Primero, necesitamos decirle a este **process** que puede aceptar una variable de entrada y escribimos _input_ aquí como un nuevo bloque de declaración. Y vamos a llamar a esto _"val greeting"_.

La parte val es el equivalente de un path aquí abajo. Le dice a Nextflow que esto es una variable, como un string en este caso. Y si pasas el cursor sobre ella de nuevo, te dice de la extensión qué significa esto.

A continuación vamos a decirle a Nextflow qué hacer con esto. No es suficiente solo decir que hay una variable. Tienes que decir en el script cómo usar esa variable. Y así que voy a deshacerme de este string codificado aquí, y voy a poner una variable.

Voy a hacerlo rápidamente sin llaves solo para mostrarte que esto está permitido, y esta es la forma antigua de hacerlo. Pero ahora con la nueva sintaxis, realmente recomendamos ponerlo dentro de llaves así, y deja muy claro que esto está siendo interpolado por Nextflow aquí.

Genial. Así que _"input greeting"_ va a _$\{greeting\}._ Lo último es que necesitamos decirle a Nextflow a nivel de **workflow** que este **process** ahora toma una entrada. Y para hacer eso, básicamente vamos a darle una variable.

## 3.2. Configurar un parámetro de línea de comandos para capturar entrada del usuario

Podríamos codificarlo de nuevo, como Hello World, y eso funcionaría bien, pero obviamente no nos da realmente ninguna ventaja. Queríamos poder configurar esto en tiempo de ejecución, así que queremos poder hacerlo en la CLI, cuando lanzas Nextflow.

Y la forma en que hacemos eso es un concepto especial de Nextflow llamado _params_. Vamos a llamar a eso _params.input_.

Lo que esto hace es exponer esta variable input en la CLI y ahí es donde usamos un doble guion cuando lanzamos Nextflow.

Puedo llamar a esto como quiera, puedo llamarlo _hello, greeting_. No importa. Lo que haga allí se expondrá como una opción CLI cuando lancemos un pipeline. Y este es un verdadero truco de magia de Nextflow porque significa que puedes construir tu script de **workflow** muy rápidamente con estos parámetros, y esencialmente estás construyendo una CLI personalizada para tu pipeline, haciendo que sea realmente fácil personalizar diferentes opciones sobre la marcha cuando lanzas.

Entonces. Probémoslo. Volvamos a nuestra terminal. Tenemos nuestro comando _"nextflow run"_ aquí. Y ahora voy a hacer _"--input"_, que coincide con _"params.input"_ que vimos antes. Creo que en los documentos está en francés. A Geraldine le gusta hablar francés. Voy a hacerlo en sueco porque vivo en Suecia. así que voy a decir, "_Hej Världen_" y presionar enter.

Puedo usar comillas simples o dobles, solo afecta cómo Bash lo interpreta.

Ejecuta el pipeline de Nextflow exactamente de la misma manera. Puedes ver el directorio de trabajo y todo es lo mismo. Pero ahora si voy arriba a _"results/hello_world/output"_. Podemos ver nuestro bonito sueco aquí en su lugar.

Así que hemos pasado dinámicamente una entrada desde una CLI a un parámetro. Hemos pasado eso como una entrada al **process** y el **process** lo interpretó y lo puso en un bloque script, que luego ha cambiado dinámicamente la salida de ese resultado del script. Bastante genial.

Lógica bastante compleja con muy poca sintaxis aquí. Y con suerte puedes ver cómo esto ahora comienza a escalar. Y así es como realmente construimos la lógica y la capacidad de personalización de nuestros pipelines en el script de Nextflow.

## 3.4. Usar valores predeterminados para parámetros de línea de comandos

Bien, eso es genial. El problema ahora es que cada vez que ejecuto este pipeline, necesito hacer dash, input para que se ejecute.

Si intento ejecutar sin este parámetro, ahora Nextflow va a lanzar un error diciendo que necesitaba este parámetro y no fue establecido. y así no supo qué hacer.

Esta es una cosa genial nueva, por cierto. En el pasado, Nextflow simplemente se habría ejecutado con un string vacío, y habrías tenido todo tipo de errores extraños, que habrían sido difíciles de entender. Pero en el nuevo analizador sintáctico de Nextflow, es un poco más cuidadoso y te lo dice de inmediato.

Así que no siempre queremos especificar cada opción. Es una buena práctica especificar valores predeterminados sensatos. Entonces, ¿cómo hacemos eso en nuestro script?

Notarás que cuando escribimos esto, simplemente pusimos _params.input_ directamente donde lo estamos usando. Así que la solución obvia es definir un valor predeterminado, y hacemos eso en la parte superior del script aquí en un bloque params especial en el **workflow**. Esto está en el script del **workflow** aquí.

De nuevo, algo de sintaxis nueva aquí, así que presta atención. Esto es realmente genial. Tenemos el nombre del parámetro, que se esperará aquí.

Y luego después de este carácter de dos puntos, estamos definiendo un tipo de la variable. No tienes que hacer esto, puedes simplemente dejarlo en blanco, pero es realmente agradable. Le dice a Nextflow que estamos esperando un string y tratarlo como tal.

Si queremos un número en su lugar, por ejemplo, podríamos escribir float, y eso diría que queremos un número de punto flotante. Y si tratamos de ejecutar con eso, entonces lanzará un error. Si le damos un string, que no es un float. Y también lo pasará como tal. Como si hacemos string, entonces sabe que es un string. E incluso si tiene ceros iniciales y es todo numérico, aún lo pasará como un string real.

Así que esa seguridad de tipo es una característica muy nueva de Nextflow, pero realmente poderosa para hacer tu código más seguro de escribir y ejecutar.

Luego después de eso tenemos un símbolo de igual y luego el valor predeterminado aquí. Nextflow fue escrito en Barcelona originalmente, así que parece apropiado que tengamos algo de español aquí, _"Holà mundo!"_ como predeterminado.

Bien, voy a guardar ese script, volver, ejecutar el script de nuevo sin _--input_. Y esta vez debería ejecutarse y creará nuestro nuevo archivo en _results_. Y en este archivo ahora dice _"Holà mundo!"_.

Este es solo un valor predeterminado, así que no significa que no podamos hacer lo mismo que antes. Si vuelvo y encuentro mi viejo script aquí, _"Hej Världen"_, porque hago _--input_ en la línea de comandos, eso sobrescribirá ese valor predeterminado y lo usará de nuevo en el archivo output.txt.

Así que esto en el script es solo el valor predeterminado que estoy estableciendo.

A medida que construyamos nuestro **workflow** para que sea más complejo e incluya más parámetros, este bloque params en la parte superior del script comenzará a recopilarlos todos en un solo lugar.

Y terminas con esta simetría bastante agradable en tu script, donde efectivamente tienes todas tus entradas del **workflow** aquí y tus salidas del **workflow** en la parte inferior. Y es muy claro cuál es la interfaz de tu **workflow** al mundo exterior. Así que puedes tomar un nuevo pipeline muy rápidamente con la nueva sintaxis y entender cómo usarlo.

Una última cosa genial. No tenemos que establecer un valor predeterminado con esto. Si hacemos params input pero no establecemos un valor predeterminado, entonces le dice a Nextflow que este parámetro es requerido, y de nuevo, el pipeline fallará en ejecutarse sin él, pero te dará un mensaje de error más útil en lugar de algo sobre que es nulo.

Así que dice que estamos esperando que su entrada es requerida, pero no fue especificada en la línea de comandos. Muy bien.

Bien, así que con suerte ahora está claro cómo configurar tu pipeline de Nextflow con entradas variables y parámetros, cómo establecer el valor predeterminado, establecer los tipos, podría ser un booleano verdadero falso flag o un entero o diferentes tipos aquí. Cómo pasarlos a tu **workflow**, dónde va, y luego interpola en tu **process**. Y también sabes cómo personalizar esos en la línea de comandos cuando lanzas Nextflow. Esto está comenzando a verse más interesante que nuestro simple comando bash.

## 4. Administrar ejecuciones del workflow

Bien. ¿Qué sigue? Para la parte final de este capítulo, vamos a hablar un poco sobre cómo administrar todas las diferentes ejecuciones del **workflow**. Si miras en mi barra lateral aquí y el Explorador debajo de work, verás que he ejecutado un montón de pipelines diferentes y estos directorios work se están volviendo bastante largos, hay muchos de ellos.

Y la otra cosa es, como dije antes, cada vez que reejecuto este pipeline, está creando un nuevo conjunto de directorios work, y está reejecutando todos los **process** desde cero, lo cual es algo bueno. Ese es el comportamiento previsto. Es reproducible y está regenerando todo fresco. Pero obviamente, si estás ejecutando **process** que tardan mucho en ejecutarse, es molesto tener que siempre empezar tu pipeline desde el principio si se bloqueó a mitad de camino, o si cambias algo al final del pipeline.

## 4.1. Relanzar un workflow con -resume

Afortunadamente, Nextflow es realmente bueno sabiendo qué se ha ejecutado previamente y qué está disponible, y para reutilizar esos viejos resultados es muy simple. Simplemente agregamos un nuevo flag al final del comando _"-resume"_.

Ahora, nota que hay dos guiones en input porque ese es el parámetro. Solo hay un guion en resume porque esa es una opción central de Nextflow.

Esto confunde a la gente todo el tiempo, incluso si has estado usando Nextflow durante mucho tiempo. Así que siempre recuerda uno o dos guiones. Depende si es una opción central de Nextflow.

Bien, entonces ahora hago _-resume_ y ejecuto exactamente el mismo **workflow** de nuevo. Y esta vez debería verse prácticamente exactamente igual con una diferencia clave.

En la salida aquí, puedes ver que los resultados fueron cacheados. Y de hecho, este hash de tarea aquí es exactamente el mismo que la ejecución anterior, y simplemente reutilizó ese directorio work en su totalidad. Las entradas y las salidas y el script estaban todos sin modificar. Y así que simplemente toma ese archivo de allí y si hay pasos siguientes en el **process**, los pasaría al siguiente paso en el pipeline.

Así que todavía está ejecutando todo el pipeline de principio a fin, pero está usando resultados cacheados para cada una de esas tareas, donde puede.

Ahora, cuando haces _-resume_, simplemente reanuda la última ejecución del pipeline en tu directorio de trabajo, lo que sea que haya sido. Pero en realidad puedes reanudar desde cualquier ejecución anterior que hayas hecho allí. Y hemos hecho bastantes ahora.

## 4.2. Inspeccionar el registro de ejecuciones pasadas

Para ver todas ellas, podemos hacer _"nextflow log"_ en lugar de _"nextflow run"_, y eso nos dará una salida agradable mostrando todas estas diferentes.. Necesito hacer mi pantalla un poco más pequeña para que podamos verlo, todas estas diferentes ejecuciones cuando las hicimos, el id de sesión, el comando y todo.

Y podemos mirar aquí y podemos tomar el nombre de ejecución de cualquiera de estas y luego reanudar una de esas específicas. Así que puedo volver y puedo reanudar esa llamada _hungry_ekeblad_. Y simplemente pongo eso después del _resume_.

Si tienes curiosidad, por cierto, todos estos adjetivos y nombres de científicos están en el código fuente de Nextflow. Es una muy buena manera de obtener tu primer pull request a Nextflow yendo y encontrándolo y agregando tu científico favorito.

Y de todos modos, entonces hice eso y volvió y miró los resultados cacheados de esta ejecución del **workflow**, se dio cuenta de que todavía podía reutilizarlos, y lo hizo. Así que obtuve los resultados cacheados de nuevo.

## 4.3. Eliminar directorios work más antiguos

Eso es genial. ¿Qué pasa si quiero limpiar estos directorios work? Hay montones de ellos aquí. Hay montones de archivos. Tal vez sé de hecho que quiero reanudar desde las últimas dos ejecuciones del pipeline, pero no me importan todas las anteriores a esas.

Entonces puedo elegir uno aquí y puedo usar otro comando de Nextflow, que es _"nextflow clean"_, y puedo hacer _"nextflow clean"_, voy a hacer _"-before"_, y el nombre de ejecución particular, que en este caso fue _reverent_pike_ y voy a hacer _"-n"_, que le dice a Nextflow que solo haga una ejecución en seco. Así que solo me dice qué eliminará. Sin hacer nada en realidad, así que eliminaría estos directorios work.

Eso parece sensato. Así que voy a hacer el mismo comando de nuevo, pero en lugar de _"-n"_ haré _"-f"_ para hacer la limpieza realmente. Y esta vez realmente ha eliminado todos estos directorios. Y si entro y miro los directorios work, ahora se ve mucho más ligero. Fantástico.

Así que así es como limpiar todos tus directorios work locales de una manera bastante segura sin destruir completamente la caché. Así que todavía puedes reanudar si quieres.

Si alguna vez olvidas qué son estos flags para cada comando de Nextflow puedes hacer _"nextflow help"_, y luego el nombre del comando. Así que si hago _"nextflow help clean"_, puedes ver todas las diferentes opciones: _-after, -before, -but_, todas las diferentes formas de configurar este comportamiento de limpieza. Bastante genial.

## Conclusión

Bien, ese es el final de la parte uno de Hello Nextflow. Es un comienzo bastante intenso del curso, pero con suerte ahora tienes una comprensión bastante buena de cómo se ve un script de Nextflow; con diferentes partes clave, los **process**, los **workflows**, las salidas, y los parámetros. Sabes cómo configurarlos con sobrescrituras básicas desde la línea de comandos, cómo hacer un bloque de entrada dinámico con un script dinámico y sabes cómo administrar todas tus ejecuciones de carga de trabajo: ver qué ya has ejecutado, reanudar, limpiar. Hay muchas cosas. Has avanzado mucho. Así que si quieres tomar un descanso y dar un paseo rápido y tomar una taza de té, ahora es probablemente un buen momento. Te lo has ganado.

De aquí en adelante, estamos básicamente construyendo sobre esta base. ¿Cómo podemos hacer esto más complejo, más poderoso? ¿Cómo podemos hacerlo más flexible? Hacer las cosas que queremos hacer en nuestro análisis a escala.

## Cuestionario

Ahora si te desplazas hacia abajo a la parte uno, hello world, en la página web verás un pequeño cuestionario y esto es algo nuevo que hemos hecho para esta versión de la capacitación de Nextflow. Y puedes pasar y ponerte a prueba para verificar que hayas entendido todo el material que hemos hecho en este capítulo.

Esto no se nos envía ni nada, solo se almacena en tu navegador. Así que no sabemos cuáles son tus respuestas, pero es solo una pequeña autoevaluación para asegurarte de que no te hayas perdido nada o malinterpretado algo. Y puedes intentarlo tantas veces como quieras.

Si eres como yo, tal vez quieras quedarte en la terminal en tu instancia de VS Code, en cuyo caso puedes escribir el comando _quiz_ y luego simplemente decirle en qué capítulo estás. Así que hacemos _"Hello World"_, y luego puedes hacer exactamente las mismas preguntas del cuestionario, que están en el navegador web, pero solo en tu terminal.

Genial. Bien. Espero que disfrutes eso. Que te diviertas un poco y, te veremos en el próximo capítulo en solo un minuto para hablar todo sobre los canales de Nextflow.
