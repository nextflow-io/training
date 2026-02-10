# Resumen del curso

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Traducción asistida por IA - [más información y sugerencias](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

¡Felicidades por completar el curso de capacitación Hello Nextflow! 🎉

<div class="video-wrapper">
  <iframe width="560" height="315" src="https://www.youtube.com/embed/wnmUAfDL01E?si=y8lAedhEHWaTV4zd&amp;list=PLPZ8WHdZGxmWKozQuzr27jyMGqp9kElVK&amp;cc_load_policy=1" title="YouTube video player" frameborder="0" allow="accelerometer; autoplay; clipboard-write; encrypted-media; gyroscope; picture-in-picture; web-share" referrerpolicy="strict-origin-when-cross-origin" allowfullscreen></iframe>
</div>

/// caption
:fontawesome-brands-youtube:{ .youtube } Vea la [lista de reproducción completa en el canal de YouTube de Nextflow](https://youtube.com/playlist?list=PLPZ8WHdZGxmWKozQuzr27jyMGqp9kElVK&si=eF7cLR62goy-lc6n).

:green_book: Puede leer la [transcripción del video](./transcripts/07_next_steps.md) junto con el video.
///

## Su recorrido

Comenzó con un workflow muy básico que ejecutaba un comando codificado de forma fija.
A lo largo de seis partes, transformó ese workflow básico en un pipeline modular de múltiples pasos que ejercita características clave de Nextflow, incluyendo canales, operadores, soporte integrado para contenedores y opciones de configuración.

### Lo que construyó

- La forma final del workflow Hello toma como entrada un archivo CSV que contiene saludos de texto.
- Los cuatro pasos están implementados como procesos de Nextflow (`sayHello`, `convertToUpper`, `collectGreetings` y `cowpy`) almacenados en archivos de módulo separados.
- Los resultados se publican en un directorio llamado `results/`.
- La salida final del pipeline es un archivo de texto plano que contiene arte ASCII de un personaje diciendo los saludos en mayúsculas.

<figure class="excalidraw">
--8<-- "docs/en/docs/hello_nextflow/img/hello_pipeline_complete.svg"
</figure>

1. **`sayHello`:** Escribe cada saludo en su propio archivo de salida (_p. ej._ "Hello-output.txt")
2. **`convertToUpper`:** Convierte cada saludo a mayúsculas (_p. ej._ "HELLO")
3. **`collectGreetings`:** Recopila todos los saludos en mayúsculas en un único archivo por lotes
4. **`cowpy`:** Genera arte ASCII usando la herramienta `cowpy`

La configuración del workflow permite proporcionar entradas y parámetros de manera flexible y reproducible.

### Habilidades adquiridas

A través de este curso práctico, ha aprendido a:

- Describir y utilizar componentes centrales de Nextflow suficientes para construir un workflow simple de múltiples pasos
- Describir conceptos del siguiente nivel como operadores y factorías de canales
- Lanzar un workflow de Nextflow localmente
- Encontrar e interpretar salidas (resultados) y archivos de registro generados por Nextflow
- Solucionar problemas básicos

Ahora está equipado con el conocimiento fundamental para comenzar a desarrollar sus propios pipelines en Nextflow.

## Próximos pasos para desarrollar sus habilidades

Aquí están nuestras 3 principales sugerencias sobre qué hacer a continuación:

- Aplique Nextflow a un caso de uso de análisis científico con [Nextflow for Science](../nf4_science/index.md)
- Comience con nf-core con [Hello nf-core](../hello_nf-core/index.md)
- Explore características más avanzadas de Nextflow con las [Side Quests](../side_quests/index.md)

Finalmente, le recomendamos que eche un vistazo a [**Seqera Platform**](https://seqera.io/), una plataforma basada en la nube desarrollada por los creadores de Nextflow que hace aún más fácil lanzar y gestionar sus workflows, así como administrar sus datos y ejecutar análisis de forma interactiva en cualquier entorno.

## Encuesta de retroalimentación

Antes de continuar, ¡por favor tome un minuto para completar la encuesta del curso! Su retroalimentación nos ayuda a mejorar nuestros materiales de capacitación para todos.

[Realizar la encuesta :material-arrow-right:](survey.md){ .md-button .md-button--primary }
