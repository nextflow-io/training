# Orientación - Transcripción del Video

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Traducción asistida por IA - [más información y sugerencias](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

<div class="video-wrapper">
  <iframe width="560" height="315" src="https://www.youtube.com/embed/PIjOdFaYwWA?si=y8lAedhEHWaTV4zd&amp;list=PLPZ8WHdZGxmWKozQuzr27jyMGqp9kElVK&amp;cc_load_policy=1" title="YouTube video player" frameborder="0" allow="accelerometer; autoplay; clipboard-write; encrypted-media; gyroscope; picture-in-picture; web-share" referrerpolicy="strict-origin-when-cross-origin" allowfullscreen></iframe>
</div>

!!!note "Nota importante"

    Esta página muestra únicamente la transcripción. Para instrucciones completas paso a paso, regrese al [material del curso](../00_orientation.md).

## Bienvenida

Hola, y bienvenido/a a Hello Nextflow. Mi nombre es Phil Ewels. Soy Gerente de Producto para Software de Código Abierto en Seqera, la empresa detrás de Nextflow.

Este curso es una introducción práctica a la construcción de workflows con Nextflow. Está diseñado para personas que son completamente nuevas en Nextflow y quieren desarrollar sus propios pipelines.

Todos los ejemplos son procesamiento simple de texto, para que pueda enfocarse en los conceptos de Nextflow sin necesitar experiencia en el dominio, solo cierta familiaridad con la línea de comandos.

Vamos a repasar los fundamentos de Nextflow: escribir procesos, conectarlos en workflows de múltiples pasos, gestionar dependencias de software con contenedores, y configurar pipelines para diferentes entornos de cómputo. Al final, habrá construido un pipeline funcional desde cero.

Este curso se enfoca en _desarrollar_ pipelines. Si solo quiere _ejecutar_ pipelines existentes sin profundizar demasiado en el código, tenemos un curso más corto "Nextflow Run" que podría ser más adecuado para usted.

Una vez que domine los fundamentos aquí, también tenemos cursos de seguimiento que aplican estos conceptos a análisis científicos reales. Le enseñaremos cómo usar los pipelines y mejores prácticas de la comunidad nf-core.

Si se queda atascado, diríjase a community.seqera.io. Hay un foro comunitario activo allí con una sección dedicada solo para preguntas de capacitación. Puede usarlo en cualquier momento, sin embargo, también realizamos semanas de capacitación trimestrales con personas disponibles específicamente para ayudar. Así que si está haciendo la capacitación durante una de esas, definitivamente no sea tímido y pida ayuda.

También puede intentar pedir ayuda a Seqera AI. Es excelente para explicar código de Nextflow y ayudarle con la depuración.

Cuando esté listo para ejecutar Nextflow a escala, Seqera Platform es el mejor lugar para hacerlo. Se ejecuta en su infraestructura sin ningún bloqueo de proveedor, con todo desde lanzamiento de pipelines hasta monitoreo en tiempo real, hasta entornos de análisis interactivos. Pero por ahora, enfoquémonos solo en los fundamentos.

Bien, comencemos.

## training.nextflow.io

Bien. Lo primero a notar es que todos los cursos de capacitación en training.nextflow.io son muy interactivos. La idea es que siga el material de capacitación y mis instrucciones, y repasemos el material de capacitación juntos. Así que necesitará dos cosas: necesitará su computadora y necesitará este sitio web abierto. Y eso es prácticamente todo.

Así que esta es la página de inicio como se ve hoy cuando grabo esto. Puede ver que hay una descripción general de las diferentes cosas, el contexto, y los diferentes cursos que tenemos, cuya lista está creciendo todo el tiempo.

Nextflow para principiantes es donde estamos. Hay dos cursos dentro de aquí, Nextflow Run, que es un curso diferente, y el Hello Nextflow, que es el que nos interesa.

Y también puede ver todos los diferentes cursos en la barra lateral. Puedo saltar a Hello Nextflow, y podemos ver todos los diferentes capítulos que vamos a trabajar juntos.

Hay un par de otras cosas importantes a notar aquí. Primero, el material de capacitación está versionado, así que puede ver aquí arriba. Dice 3.0 latest, que al momento de grabar es la última versión estable. Esto cambiará con el tiempo. Publicamos nuevos cursos y actualizamos el material con el tiempo. Así que si es 3.1 o 3.2, no se preocupe demasiado. Si es 4.0, entonces probablemente hay un nuevo video, y debería tal vez ir y encontrarlo porque probablemente habrá actualizaciones significativas.

Otro menú desplegable en la parte superior es este, de idioma. Ahora esto es completamente nuevo para la versión 3.0. Hemos tomado el material previamente traducido, que fue hecho por humanos, a mano, y lo hemos pasado a un LLM y configurado toda esta nueva infraestructura para mantener diferentes traducciones del material de capacitación usando traducción LLM.

Así que ahora tenemos todas estas fantásticas traducciones aquí. Así que si quiere escuchar en coreano, puede cargar todo el sitio web en coreano. Y seguir allí. Lo mismo para todos estos otros idiomas, hindi y alemán y demás. Voy a seguir en inglés. Ese es como el idioma principal en el que escribimos el material.

Un par de otros botones si le gusta tener el modo claro. En lugar de ese modo, puede seguir el sitio web en modo claro aquí arriba.

Y luego también todo lo que miramos está en un único repositorio de GitHub, que es de código abierto, llamado nextflow-io/training. Y si hace clic en este botón en cualquier momento, irá al repositorio de GitHub. Volveremos a eso en un minuto.

## Configuración de GitHub Codespaces

Bien, ahora que tiene esto abierto en la pestaña del navegador. Vamos a Hello Nextflow y hagamos clic. Puede ver en la página de introducción, nos dice algunos de los requisitos, la descripción general, y el plan de lecciones de aproximadamente lo que vamos a cubrir, y luego vamos a sumergirnos en los primeros pasos.

Hay diferentes formas en que puede hacer este tutorial interactivo. Si está dispuesto, es bienvenido a hacer esto localmente en su propia computadora con su propia instalación de Nextflow. Si hacemos clic en Opciones de Entorno, puede ver que hay más detalles sobre cómo hacer esto ya sea usando Devcontainers locales o también puede simplemente instalar todo el software localmente, con instalación manual.

Estamos trabajando en hacer que esto funcione bien con Seqera Studios, así que esa es otra opción. Pero la más común ahora mismo es usar GitHub Codespaces.

Codespaces, configura un entorno sandbox en un servidor remoto ejecutado por GitHub. Y es gratuito para una cierta cantidad de uso, que generalmente está bien para capacitación. Y le configurará una instancia de VS Code, un IDE donde puede acceder a todos los archivos del repositorio, ejecutar Nextflow y todo. Y hemos preconfigurado Codespaces para usted. Así que tiene todo lo que necesita.

La belleza de esto es que es solo un clic para configurar un Codespace. Es lo mismo para todos, y sabemos que ya tiene todos los requisitos previos instalados, así que es rápido y agradable.

Así que lo primero que debe hacer es ir a "Primeros Pasos". Busque este botón, que dice, _Abrir en Codespaces_. Voy a hacer comando \+ clic para abrirlo en una nueva pestaña, y nos lleva a GitHub.

Así es como se ve. Podemos ver, hemos configurado todas las opciones aquí para usted. Si quiere, puede hacer clic en cambiar opciones. Algunas cosas que puede hacer aquí. Puede dar una máquina de instancia más grande, por ejemplo, si encuentra que se bloquea porque se queda sin memoria o algo así. O establecer versiones específicas del material de capacitación. Pero generalmente puede simplemente ir con lo que hemos configurado aquí y puede verlo. En este caso está, usando la versión 3.0.

Así que voy a hacer clic en crear nuevo Codespace. Y eso me lleva adentro.

Note también, dice no hay Codespace para reanudar allí. Si he creado previamente un Codespace, hacer clic en ese botón nuevamente en el material de capacitación me llevará a la misma página y listará todos los Codespaces que ya tengo ejecutándose. Entonces puede simplemente volver directamente a ellos y continuar donde lo dejó. Así que no importa si cerró su computadora portátil.

Se apagan automáticamente después de unos minutos de inactividad, pero no hay problema. Simplemente puede reiniciarlos.

Una vez que inicie un nuevo Codespace, va a quedarse en esta página así y va a cargar por bastante tiempo. Así que ahora es un buen momento para tomar un breve descanso. ¿Tal vez olvidó ir al baño o quiere una taza de té antes de comenzar? Vaya ahora mientras espera esto, porque va a girar allí por un tiempo.

Rápidamente mientras esperamos que cargue, también voy a ir a github.com/codespaces y solo mostrar esta es la página de descripción general donde puede ver todos los diferentes Codespaces que tiene ejecutándose actualmente.

Puede ver que tengo uno aquí para nextflow-io/training. Sin cambios, porque no he hecho nada en él todavía. La cantidad de recursos que está usando, y puede ver en este momento se está configurando. Puedo ir aquí, hacer clic en este pequeño menú desplegable y hacer clic en eliminar. Así que si accidentalmente configuró múltiples Codespaces y no está usando algunos, puede eliminar los antiguos y limpiar.

Finalmente, una forma más de entrar a esto. Si vamos al repositorio de GitHub. Y esto funciona para cualquier repositorio de GitHub. Haga clic en code. Puede tener comandos para clonar el repositorio localmente. Y hay una pestaña llamada Codespaces. Y nuevamente, puede crear uno nuevo, y puede ver los que ya están ejecutándose.

Así que nuevamente, si olvida cómo creó su Codespace, siempre puede volver a él de esta manera.

## La interfaz de VS Code

Bien, los constructores terminaron y ahora está comenzando a cargar los GitHub Codespaces. No siempre toma tanto tiempo, así que no se preocupe. Es solo la primera vez que crea el Codespace. Si vuelve a uno que ya existe, es mucho más rápido.

No sea demasiado impaciente si esta es la primera vez, aún no ha terminado, aunque está comenzando a darnos una interfaz.

Pero mientras esperamos que se configuren las cosas finales, solo le mostraré la interfaz en caso de que no esté familiarizado con VS Code.

Primero, está la barra lateral de chat para cosas de IA, que no necesitamos. Así que voy a cerrar eso, deshacerme de eso y liberar algo de espacio.

A la izquierda, tenemos el explorador de archivos que nos muestra todos los archivos en el repositorio de Git, que es el espacio de trabajo que hemos creado. Note, estos no son archivos locales. Todo esto está en el servidor remoto donde estamos trabajando. Puede arrastrar y soltar archivos locales y cosas, pero en su mayor parte, no vamos a pensar en eso hoy. Solo estamos trabajando puramente de forma remota.

Hay otras herramientas en esta barra lateral, por ejemplo, búsqueda. Así que puede buscar todos los archivos en un repositorio de una vez. Y si estuviéramos haciendo trabajo de desarrollo en el repositorio de capacitación, podríamos hacer integración con control de fuente con Git y depuración y otras cosas.

Otras cosas son, hay una ventana principal de edición de código aquí arriba, que acaba de cargar una vista previa del readme, que es para el material de capacitación. Así que en este caso está viendo markdown, pero normalmente esto será un editor de código.

Y luego debajo de eso tenemos la terminal, que es donde vamos a estar ejecutando todos nuestros comandos e interactuando directamente con Nextflow.

Todo en el Codespace está preinstalado, así que el comando Nextflow ya está allí y demás.

Bien. Cuando llegue hasta aquí, debería estar casi listo. Puede ver ahora que ha descargado el servidor de lenguaje Nextflow y ha configurado algunas extensiones para nosotros en VS code, incluyendo la extensión de Nextflow, que va a ser útil. Así que puedo cerrar eso y puedo cerrar el README.md.

Y ahora puede ver que tengo algo más en el lado izquierdo. Estoy un poco ampliado aquí, pero si reduzco el zoom puede ver que uno de los botones dice Nextflow con el ícono de Nextflow. y eso tiene algunas cosas agradables aquí para explorar el proyecto y cosas, a las que volveremos más tarde.

Bien. en caso de que pierda alguno de estos paneles, estos botones en la parte superior derecha son realmente útiles y estos solo muestran y ocultan cosas. Así que eso muestra y oculta el Explorador, muestra y oculta la terminal en la parte inferior. Y demás.

Voy a estar usando estos bastante porque estoy muy ampliado, así que trato de ayudarle a ver todo el texto en mi pantalla, y por lo tanto es útil poder ir a pantalla completa con la terminal y luego ocultarla cuando estemos mirando código. Pero la mayor parte del tiempo puede simplemente tener todas estas cosas abiertas al mismo tiempo.

Bien, ¿qué más ver? No mucho más. Note que Nextflow, como digo, está instalado. Así que puedo escribir "nextflow -version" y debería aparecer diciendo qué versión tenemos instalada.

Hay algunas otras cosas instaladas aquí también. Al final de cada capítulo, tenemos un conjunto de preguntas de cuestionario, por ejemplo, en el sitio web. Y también puede hacer esas en la terminal si quiere escribiendo quiz.

Hay algunos otros atajos de teclado que voy a estar usando, solo en caso de que tenga curiosidad. Por ejemplo, justo entonces presioné cmd\+K en mi Mac y eso limpió la terminal, para deshacerse de toda la salida anterior. Así que eso es agradable para mantener las cosas limpias. Si me ve haciendo eso, así es como lo estoy haciendo.

Y también si es nuevo en la terminal, recuerde que puede usar tab para autocompletar, lo cual estaré haciendo mucho para autocompletar rutas.

Así que puedo ver en el lado izquierdo aquí hay una carpeta llamada Hello Nextflow, que es lo que vamos a estar trabajando. Si hago "ls" para listar archivos, puedo hacer "hel", presionar tab, autocompleta. Y así esta es una forma muy rápida de completar rutas.

## Abriendo solo la carpeta Hello Nextflow

Bien. Esto es genial. Sin embargo, hay muchas cosas en este repositorio.

Están todos los archivos para generar el sitio web, y hay múltiples cursos diferentes aquí, y puede hacerlo desde esta ruta y simplemente hacer clic en la carpeta "Hello Nextflow". Pero es agradable enfocarse puramente en esto.

Puede establecer esto como su espacio de trabajo con un montón de clics por aquí y configurando un directorio de proyecto y cosas. Pero la forma más fácil es escribir code, que es el comando CLI para lanzar VS Code, y luego "hello-nextflow".

Eso abrirá una nueva pestaña del navegador y puede cerrar la antigua. Y se ve exactamente igual. Pero ahora puede ver que estamos en este subdirectorio y todos los otros archivos son invisibles, y tenemos una configuración más limpia.

Puede ver aquí que también el directorio de trabajo actual ahora está dentro de la carpeta Hello Nextflow. Así que agradable y limpio. No necesitamos preocuparnos por estar en el lugar equivocado. Bien.

## Nueva Sintaxis de Nextflow para 2026

hay una cosa especial que necesito mencionar en este punto. Ahora mismo, a principios de 2026, estamos comenzando a traer diferentes características a Nextflow, y una de las grandes nuevas es un nuevo analizador de sintaxis de lenguaje dentro de Nextflow.

Básicamente el motor que lee sus archivos de Nextflow y entiende eso, para tiempo de ejecución. Hay algunos cambios en la sintaxis, y es realmente importante que use Nextflow con el analizador de sintaxis correcto habilitado.

Dos cosas que necesita para esto. Necesita una versión actualizada de Nextflow y necesita asegurarse de que esté habilitado.

Si hago "nextflow -version" nuevamente, verá que los Codespaces están ejecutándose con 25.10.2 y 25.10 es la versión mínima para poder estar usando estas cosas.

Si está usando 26.04, que para mí aún no ha salido, pero lo hará pronto. Entonces esto estará ejecutando el nuevo analizador de sintaxis por defecto, y no tiene que hacer nada más.

Pero si está ejecutando 25.10, necesita habilitar el analizador de sintaxis estricto, como se llama, o analizador de sintaxis v2.

Esto se hace con una variable de entorno. Ya está configurada en los Codespaces, así que no necesita hacer nada. Pero si está ejecutando localmente, necesita configurar esto, y puedo verificar esto haciendo "echo $NXF_SYNTAX_PARSER", y debería estar configurado en v2.

Así que si está ejecutando localmente, solo haga "export NXF_SYNTAX_PARSER=v2". Así de simple. Pero recuerde hacer eso. porque de lo contrario va a ver algunas discrepancias raras y, errores a medida que avanzamos.

Si tiene alguna duda sobre cualquiera de estas cosas alrededor de la versión de Nextflow y el analizador de sintaxis, primero, recuerde, no necesita preocuparse si está en Codespaces. Todo debería estar configurado correctamente. Pero segundo, si va al material de capacitación de Nextflow, si baja, habla sobre requisitos de versión, hay un enlace aquí que lo lleva a la página de ayuda alrededor de explorar versiones, y esto básicamente repasa todo en detalle.

Vale la pena leer esto si tiene un momento. porque ayuda a aclarar cuáles son algunos de los diferentes términos, que podría escuchar cuando comience a usar Nextflow. Cosas como DSL1, DSL2, analizador de sintaxis uno, analizador de sintaxis dos, y demás. Así que vale la pena solo echar un vistazo a eso y eso repite algo de lo que acabo de decir.

También es realmente útil si ha escrito previamente código de Nextflow y está regresando para un repaso. Le dice algunas de las cosas que cambian y lo enlaza a partes de la documentación de Nextflow, que le dice cómo actualizar su código de Nextflow.

## Archivos del curso

Bien. Lo último para familiarizarnos es solo ver los archivos, que están en este directorio. Puede mirar en la barra lateral o a menudo en el material de capacitación, usamos el comando tree, -L, que es el número de niveles para mirar. Diremos dos, y si hago esto en pantalla completa, verá que esto básicamente refleja exactamente lo que vemos en la barra lateral allí, pero excluye archivos ocultos, que comienzan con un punto.

Así que los archivos \*.nf, significa Nextflow. Así que estos son los archivos de script de Nextflow, y hay un archivo de inicio aquí para cada uno de los diferentes capítulos del material de capacitación, que abriremos y exploraremos y luego editaremos.

Cambiaremos estos archivos a medida que avancemos, y así al final de cada capítulo, los archivos deberían verse prácticamente igual que el inicio del capítulo para el siguiente. Pero le damos estos diferentes archivos para que siempre pueda comenzar de nuevo y no preocuparse demasiado por estropear la sintaxis.

Si necesita comparar con algo que definitivamente debería funcionar. Puede verificar en la carpeta de soluciones, y este es como un estado final para cada uno de los capítulos, así que puede comparar lo que ha escrito contra lo que está allí.

Hay un directorio de datos. Este tiene solo un archivo greetings.csv, que usaremos como ejemplo, datos de entrada en parte del curso, y cosas como un archivo de configuración y algunos parámetros, que describiremos más adelante en el curso.

## Conclusión

Bien, así que ahora esperamos que todo esté funcionando. Su pantalla se ve igual que la mía y entiende cómo acceder a todo y cuáles son todos los diferentes archivos.

Si se desplaza hacia abajo hasta la parte inferior de la página en primeros pasos, pequeña casilla de verificación debería decir que entiendo lo que estoy haciendo. Mi entorno está funcionando y ha configurado, su directorio de trabajo correctamente a la carpeta "Hello Nextflow".

Si ha marcado todos esos y se ven verdes. Podemos continuar al siguiente video y el siguiente capítulo, que es la parte uno. Hello World. Nos vemos en un momento.
