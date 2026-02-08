# Parte 4: Hello Modules - Transcripción del video

<div class="video-wrapper">
  <iframe width="560" height="315" src="https://www.youtube.com/embed/43Ot-f0iOME?si=0AWnXB7xqHAzJdJV&amp;list=PLPZ8WHdZGxmWKozQuzr27jyMGqp9kElVK&amp;cc_load_policy=1" title="YouTube video player" frameborder="0" allow="accelerometer; autoplay; clipboard-write; encrypted-media; gyroscope; picture-in-picture; web-share" referrerpolicy="strict-origin-when-cross-origin" allowfullscreen></iframe>
</div>

!!!note "Notas importantes"

    Esta página muestra únicamente la transcripción. Para instrucciones paso a paso completas, regrese al [material del curso](../04_hello_modules.md).

    Los números de sección mostrados en la transcripción se proporcionan solo con fines indicativos y pueden no incluir todos los números de sección en los materiales.

## Bienvenida

Hola, y bienvenidos de nuevo a la parte cuatro de Hello Nextflow. Esta sección trata completamente sobre módulos, y es una sección bastante corta del curso. En realidad no vamos a escribir mucho código, es más sobre cómo organizamos el código en nuestro pipeline.

Hasta ahora, hemos estado colocando todo en un único archivo, lo cual está bien, y así es como solíamos construir pipelines de Nextflow en los viejos tiempos.

Pero a medida que ese pipeline crece, el script se vuelve cada vez más largo y más difícil de navegar, de mantener, y también significa que realmente no podemos compartir nada del código.

Los módulos de Nextflow nos permiten extraer los procesos de ese script principal y luego importarlos. Eso significa que el código es más fácil de navegar y también significa que podemos compartir ese código de módulo entre diferentes pipelines.

Este pequeño diagrama en la página principal de la documentación muestra el concepto claramente. En lugar de un script enorme, vamos a incluir estos archivos de módulo separados, de diferentes scripts de módulo, y todo se va a integrar en el workflow, pero aún va a ejecutarse exactamente de la misma manera.

Así que saltemos a GitHub Codespaces y echemos un vistazo. Como antes, he limpiado un poco mi espacio de trabajo aquí. Eliminé los directorios antiguos de Nextflow y el directorio work y demás. Pero no importa si todavía tienes esos archivos.

Voy a comenzar a trabajar en el archivo hello modules, que es básicamente donde lo dejamos al final del capítulo anterior. Tenemos nuestros tres procesos aquí. Tenemos un par de params, el bloque workflow, donde estamos ejecutando esos tres procesos y conectándolos con canales. Luego publicamos los canales de salida y tenemos el bloque output diciendo cómo publicar esos archivos.

## 1. Crear un directorio para almacenar módulos

Ahora, como dije, realmente no vamos a escribir o editar mucho código. Solo vamos a reorganizar el código que ya tenemos. Los archivos de módulo de Nextflow típicamente tienen un único proceso en ellos, y por convención normalmente los mantenemos en un directorio llamado modules. Pero puedes llamarlo como quieras. Pero voy a mantener un directorio modules en mi repositorio aquí, y luego voy a crear un archivo para cada proceso. Así que voy a decir archivo nuevo, sayHello.nf.

## 2. Crear un módulo para sayHello()

Ahora voy a tomar mi proceso y simplemente voy a seleccionar este código, cortarlo del archivo principal hello modules y pegarlo aquí.

Obviamente eso no hace nada por sí solo. Nuestro script principal todavía necesita ese proceso, así que necesitamos traerlo de vuelta de alguna manera. Y hacemos eso con la declaración include.

Entonces escribo include y unas llaves, y luego tomo el nombre del proceso. Y digo from, y luego le doy una ruta de archivo relativa. Así que dice, comienza con ./ porque es relativa desde donde este script está guardado. Entonces es modules sayHello.nf.

Noten que la extensión de VS Code es bastante útil aquí. Me dice, si puede encontrar este archivo y si puede encontrar un proceso, que estoy nombrando. Si deliberadamente pongo un error tipográfico aquí, me da un error de inmediato y me dirá que no puede encontrar este proceso que estoy tratando de importar. Así que mantén un ojo en cualquier error que encuentres.

Y eso es realmente todo. Todavía tenemos nuestro proceso aquí. No se necesitan cambios aquí abajo. El proceso tiene el mismo nombre y se ejecuta exactamente de la misma manera. Solo que el código real del proceso ahora está en un archivo separado.

Podemos ejecutar el workflow de Nextflow nuevamente, va a funcionar exactamente de la misma manera. Y eso es básicamente el resto de este capítulo del curso, solo mover estos tres procesos a sus propios archivos.

Así que hagámoslo ahora. Voy a crear rápidamente un nuevo archivo de módulo para el segundo proceso: convertToUpper.nf. Voy a cortar ese código, pegarlo ahí. Y luego voy a incluir ese. Vamos, genial.

Y luego voy a crear un nuevo archivo para collectGreetings.nf. Cortar eso.

Mucho cortar, cortar y copiar y pegar.

Y ahora nuestro script de workflow principal de repente se ve mucho, mucho más corto, mucho más accesible y mucho más fácil de leer.

Y pueden ver cómo el proyecto ahora comienza a construirse con nuestros diferentes archivos. Podemos profundizar en los detalles en los lugares que queremos. Navegar para encontrar pasos específicos en el pipeline mucho más fácilmente, y obtener una visión general de lo que está haciendo el pipeline rápidamente.

## Navegar módulos con VS Code

Ahora, por supuesto, la desventaja de hacer esto es que si tienes un pipeline grande, tendrás muchos archivos de módulo y podrían estar organizados en múltiples subdirectorios o todo tipo de cosas. Ahora, de nuevo, un pequeño consejo aquí. La extensión de VS Code es bastante buena para navegar tu base de código por ti y también informarte sobre el código ahí.

Pueden ver que VS Code entiende qué es este proceso y me da una pequeña descripción general cuando paso el cursor, así puedo ver sin tener que ir y encontrar el código fuente, cuáles son las entradas y las salidas, que típicamente es lo más importante cuando lo estoy usando en un workflow.

Y también si mantengo presionado command, estoy en una Mac, y hago clic en el nombre del proceso, abre el archivo directamente de inmediato. Lo trae. Así que puedo saltar directamente allí sin siquiera pensar en cuáles son las rutas de archivo reales. Y eso funciona en cualquier lugar, también puedo hacer eso, donde sea que se estén llamando los procesos. Así que eso lo hace realmente rápido.

## 4.4. Ejecutar el workflow

Bien, solo verifiquemos que el pipeline todavía se ejecuta como esperamos. Así que abro la terminal. Hagamos "nextflow run hello modules", y veamos si se ejecuta sin ningún problema.

Esperemos que todo el punto de esto sea que el pipeline está básicamente sin cambios, así que realmente no deberías ver ningún cambio desde cuando lo ejecutamos antes. La salida aquí se ve exactamente igual, y puedes ver nuestro directorio results con todos los mismos archivos, así que eso es genial. Ningún cambio es bueno.

## Una nota sobre nf-core/modules

Justo antes de terminar, quiero tocar rápidamente el poder de la colaboración cuando se trata de módulos. Estos archivos están en mi repositorio, por lo que no es obvio de inmediato cómo podríamos colaborar en ellos. Y hay muchas formas diferentes en que puedes hacer esto, pero probablemente el ejemplo más grande y mejor conocido de esto es nf-core.

Si voy al sitio web de nf-core, voy a recursos, y módulos. Pueden ver que nf-core tiene una biblioteca enorme de módulos, casi 1700 módulos cuando veo esto. Y así puedo escribir el nombre de cualquiera de mis herramientas favoritas, ir y encontrar si alguien más ya ha escrito un módulo para ella, y ver este proceso de módulo pre-escrito aquí con todas las entradas, las salidas, los contenedores de software, toda esta información, y pueden ver en el lado aquí, cuántos pipelines diferentes de nf-core están todos usando este único proceso compartido.

Este es un ejemplo un poco extremo, pero pueden ver que esto es realmente reutilizar este código. Y si hago clic para ir al código fuente de GitHub para esto, es exactamente lo mismo que lo que estamos haciendo. Es solo un gran proceso en un archivo.

Ahora en el lado de nf-core, hacemos algunos trucos para poder compartir esos archivos y traerlos a diferentes repositorios. Y si quieres saber más sobre eso, ve y revisa el curso que tenemos sobre usar y construir con nf-core específicamente. Pero quería solo darte una idea de qué tan poderoso puede ser este concepto de reutilización de código.

## Conclusión

Bien, eso es todo para módulos. Les dije que era una sección corta del curso. Revisen el cuestionario, asegúrense de entenderlo y asegúrense de que todo todavía esté funcionando correctamente. Y los veré de vuelta en el siguiente video, que trata completamente sobre contenedores de software. Muchas gracias.

I.
