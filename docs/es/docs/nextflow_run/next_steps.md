# Resumen del curso

¡Felicidades por completar el curso de capacitación Nextflow Run! 🎉

<!-- placeholder for video -->

## Tu recorrido

Comenzaste con un workflow muy básico y aprendiste a ejecutarlo, encontrar las salidas y gestionar su ejecución.
Luego, trabajaste con versiones cada vez más complejas de ese workflow y aprendiste a reconocer los conceptos y mecanismos esenciales que impulsan los pipelines de Nextflow, incluyendo canales y operadores, modularización de código y contenedores.
Finalmente, aprendiste cómo personalizar la configuración de un pipeline para adaptarlo a tus preferencias y tu infraestructura computacional.

### Lo que aprendiste

Ahora eres capaz de gestionar la ejecución del pipeline Hello, describir cómo está estructurado e identificar las principales piezas de código involucradas.

- La forma final del workflow Hello toma como entrada un archivo CSV que contiene saludos de texto.
- Los cuatro pasos están implementados como procesos de Nextflow (`sayHello`, `convertToUpper`, `collectGreetings` y `cowpy`) almacenados en archivos de módulos separados.
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

A través de este curso práctico, has aprendido cómo:

- Ejecutar un workflow de Nextflow localmente
- Encontrar e interpretar salidas (resultados) y archivos de registro generados por Nextflow
- Reconocer los componentes centrales de Nextflow que constituyen un workflow simple de múltiples pasos
- Describir conceptos de siguiente nivel como operadores y factorías de canales
- Configurar pipelines para diferentes entornos computacionales

Ahora estás equipado/a con el conocimiento fundamental para comenzar a integrar pipelines de Nextflow existentes en tu propio trabajo.

## Próximos pasos para desarrollar tus habilidades

Aquí están nuestras principales sugerencias sobre qué hacer a continuación:

- ¡No solo ejecutes Nextflow, escríbelo! Conviértete en desarrollador/a de Nextflow con [Hello Nextflow](../hello_nextflow/index.md)
- Aplica Nextflow a un caso de uso de análisis científico con [Nextflow for Science](../nf4_science/index.md)
- Comienza con nf-core con [Hello nf-core](../hello_nf-core/index.md)
- Aprende técnicas de resolución de problemas con la [Misión Secundaria de Depuración](../side_quests/debugging.md)

Finalmente, te recomendamos que eches un vistazo a [**Seqera Platform**](https://seqera.io/), una plataforma basada en la nube desarrollada por los creadores de Nextflow que facilita aún más el lanzamiento y la gestión de tus workflows, así como la gestión de tus datos y la ejecución de análisis de forma interactiva en cualquier entorno.

## Obtener ayuda

Para recursos de ayuda y soporte de la comunidad, consulta la [página de Ayuda](../help.md).

## Encuesta de retroalimentación

Antes de continuar, ¡por favor tómate un minuto para completar la encuesta del curso! Tus comentarios nos ayudan a mejorar nuestros materiales de capacitación para todos.

[Realizar la encuesta :material-arrow-right:](survey.md){ .md-button .md-button--primary }
