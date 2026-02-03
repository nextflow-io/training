# Parte 4: Hola Módulos - Transcripción

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Traducción asistida por IA - [más información y sugerencias](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

<div class="video-wrapper">
  <iframe width="560" height="315" src="https://www.youtube.com/embed/Xxp_menS0E8?si=0AWnXB7xqHAzJdJV&amp;list=PLPZ8WHdZGxmXiHf8B26oB_fTfoKQdhlik" title="YouTube video player" frameborder="0" allow="accelerometer; autoplay; clipboard-write; encrypted-media; gyroscope; picture-in-picture; web-share" referrerpolicy="strict-origin-when-cross-origin" allowfullscreen></iframe>
</div>

!!!note "Notas importantes"

    Esta página muestra solo la transcripción. Para instrucciones completas paso a paso, regrese al [material del curso](../04_hello_modules.md).

    Los números de sección mostrados en la transcripción se proporcionan solo con fines indicativos y pueden no incluir todos los números de sección en los materiales.

## Bienvenida

Hola, bienvenido a la Parte Cuatro del curso de entrenamiento Hello Nextflow.

Este capítulo se llama Hola Módulos, y hablaremos sobre cómo modularizar código Nextflow. Lo que vamos a hacer es tomar nuestro script de workflow y dividirlo en archivos separados.

Esto hace que el código sea más fácil de navegar y mantener a medida que su workflow se hace más grande, y también hace posible compartir módulos entre pipelines para que si tiene múltiples pipelines usando la misma herramienta, solo necesite escribir ese process una vez.

Un ejemplo clásico de esto es el repositorio de módulos nf-core, que tiene miles de herramientas diferentes en módulos listos para usar, los cuales puede instalar y utilizar en su workflow.

Nextflow también puede trabajar con sub workflows, que son como módulos, pero tienen múltiples processes. Eso está fuera del alcance de este entrenamiento, pero funciona básicamente de la misma manera.

Muy bien. Echemos un vistazo.

Como siempre, comience yendo a training.nextflow.io.

Vaya a "Hello Nextflow" en la barra lateral, y estamos en la parte cuatro: "Hello Modules".

Ahora voy a saltar a mi entorno de GitHub Code Spaces y echar un vistazo al archivo "hello-modules".

Como antes, estamos comenzando en el punto final del capítulo anterior, así que este script debería resultarle familiar. Tenemos nuestros tres processes, say hello, convert to upper y collect greetings, y en un simple workflow, que ejecuta estos tres comandos y emite un mensaje al final. Tenemos dos parámetros llamados greeting y batch, que especifica el nombre, que se utiliza para el archivo de salida recopilado al final.

## 0. Calentamiento: Ejecutar hello-modules.nf

Podemos verificar que este workflow todavía funciona como esperamos haciendo nextflow run hello, modules.

Genial. Ejecutó tres tareas con cada uno de estos processes, una tarea recopilada, y nos dijo que hay tres saludos en este lote. Si entramos en results, tenemos nuestros diferentes archivos de salida aquí, incluida la salida recopilada test batch.

## 1. Crear un directorio para almacenar módulos

Bien. Hagamos algo de modularización.

Generalmente es una buena idea poner los módulos en una subcarpeta en su repositorio de pipeline, solo para mantener las cosas ordenadas. Puede llamarla como quiera, pero por convención usualmente la llamamos modules.

Así que vamos adelante, vaya a una terminal y haga make the modules. Puede verla aparecer en la barra lateral y VS Code aquí.

## 2. Crear un módulo para sayHello()

Entonces voy a crear un nuevo archivo para mi primer módulo. Puede hacer "touch" o "code" o puede hacerlo en la barra lateral, realmente no importa. Así que voy a hacer code modules y voy a nombrarlo según el process. Así que sayHello.nf. NF es una extensión de archivo tradicional para archivos Nextflow.

Voy a presionar guardar aquí y podemos ver que aparece nuestro nuevo archivo de módulo.

## 2.2. Mover el código del process sayHello al archivo del módulo

Bien, a continuación voy a tomar el código del módulo del workflow. También voy a tomar el hash bang aquí y copiarlo primero para que sea claramente un archivo Nextflow. Y luego voy a tomar este process y voy a cortar. Así que voy a eliminarlo de mi script principal de workflow y voy a pegarlo en este nuevo módulo.

Ese es todo el contenido que este archivo de módulo va a contener. Solo un process único, no workflow, no lógica, solo un process solo.

Ahora puedo cerrar este archivo.

## 2.3. Agregar una declaración de importación antes del bloque workflow

Ahora mi workflow está faltando ese primer process, así que necesitamos traerlo de vuelta importándolo. La sintaxis para esto es muy similar a otros lenguajes de programación, así que puede sentirse familiar. Hacemos include llaves, el nombre del process, say hello, y luego from la ruta del archivo modules, say hello, nf. Fantástico.

Un par de trucos aquí. La extensión de VS Code es inteligente con esto. Reconoce esta ruta de archivo y puede pasar el cursor sobre ella y hacer follow link. O estoy en Mac, puedo hacer option click y abre este archivo. Así que podemos saltar rápidamente a él.

Este nombre de process ahora está siendo usado por el workflow aquí abajo, y podemos hacer lo mismo aquí. Nos muestra un poco de información sobre ese process, y nuevamente, puedo mantener presionado option, hacer clic en él, y lo abrirá en el editor.

Así que es una forma realmente rápida cuando tiene muchos archivos para sus diferentes processes de navegar rápidamente por su base de código en VS Code.

Bien. Eso es básicamente todo para este capítulo. Ahora solo hacemos lo mismo nuevamente para los otros processes.

## 3. Modularizar el process convertToUpper()

Así que vamos a crear un nuevo archivo aquí. Llamémoslo Convert to upper nf. Nuevamente, copie el hash bang. Y luego corte el process.

Copie el nombre del process allí, incluya una nueva declaración include con el nuevo nombre de process.

## 4. Modularizar el process collectGreetings()

Y luego haga lo mismo para el tercer process. Nuevo archivo, connect. Greetings,

haga el hash bang. Corte el process, pegue el process, y haga una nueva declaración include.

Ahora puede ver aquí tengo un subrayado de error aquí diciendo invalid include source. Y este es en realidad un error genuino que cometí porque me estaba moviendo un poco demasiado rápido. Si mira de cerca, puede ver que me perdí la T en convert to upper

Así que VS Code muy útilmente me ha dicho que cometí un error allí. Si corrijo ese nombre de archivo, el error desaparece. Es un buen ejemplo de por qué la verificación de errores dentro de VS Code es tan útil para escribir código Nextflow. De lo contrario no lo habría detectado y solo me habría dado cuenta mucho más tarde cuando intenté ejecutar el workflow.

Nuestro script principal de pipeline ahora se ve mucho más simple. No tiene ningún process, solo tenemos tres declaraciones include y nuestro workflow. No hemos cambiado ninguna de la lógica del workflow. No hemos cambiado ningún código del process, así que con suerte debería funcionar exactamente de la misma manera.

## 4.4. Ejecutar el workflow para verificar que hace lo mismo que antes

Vamos a verificar. Voy a abrir una terminal y voy a ejecutar exactamente el mismo comando que antes.

Efectivamente, ha ejecutado nuestros processes, say hello, convert to upper collect greetings, y nos dio tres saludos nuevamente.

Así que hemos movido nuestro código, pero no hemos cambiado nada sobre cómo se ejecuta el workflow y está completamente sin cambios. La única diferencia es que ahora tenemos código más limpio, más fácil de mantener, y más fácil de compartir con otros.

Y eso es todo. Fue un capítulo corto. Es un concepto simple, pero es muy poderoso y clave para cómo escribimos workflows Nextflow más complejos. Así que es importante que lo entienda y adquiera el hábito de usarlo.

En el siguiente capítulo, vamos a tener un poco de cambio de ritmo y dejar de pensar tanto en la sintaxis de escribir código Nextflow, y pensar un poco sobre cómo usamos software en los processes mismos. Únase a nosotros en la parte cinco para Hola Contenedores.

[Siguiente transcripción de video :octicons-arrow-right-24:](05_hello_containers.md)
