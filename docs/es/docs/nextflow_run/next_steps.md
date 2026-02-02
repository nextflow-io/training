# Resumen del curso

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Traducción asistida por IA - [más información y sugerencias](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

¡Felicidades por completar el curso de entrenamiento Nextflow Run!

<!-- placeholder for video -->

## Su recorrido

Comenzó con un workflow muy básico, y aprendió a ejecutarlo, encontrar las salidas y gestionar su ejecución.
Luego, trabajó a través de versiones cada vez más complejas de ese workflow y aprendió a reconocer los conceptos y mecanismos esenciales que impulsan los pipelines de Nextflow, incluyendo channels y operadores, modularización de código y contenedores.
Finalmente, aprendió cómo personalizar la configuración de un pipeline para adaptarlo a sus preferencias y su infraestructura computacional.

### Lo que aprendió

Ahora puede gestionar la ejecución del pipeline Hello, describir cómo está estructurado e identificar las piezas principales de código involucradas.

- La forma final del workflow Hello toma como entrada un archivo CSV que contiene saludos de texto.
- Los cuatro pasos están implementados como processes de Nextflow (`sayHello`, `convertToUpper`, `collectGreetings` y `cowpy`) almacenados en archivos de módulos separados.
- Los resultados se publican en un directorio llamado `results/`.
- La salida final del pipeline es un archivo de texto plano que contiene arte ASCII de un personaje diciendo los saludos en mayúsculas.

<figure class="excalidraw">
--8<-- "docs/en/docs/hello_nextflow/img/hello_pipeline_complete.svg"
</figure>

1. **`sayHello`:** Escribe cada saludo en su propio archivo de salida (_ej._ "Hello-output.txt")
2. **`convertToUpper`:** Convierte cada saludo a mayúsculas (_ej._ "HELLO")
3. **`collectGreetings`:** Recolecta todos los saludos en mayúsculas en un único archivo de lote
4. **`cowpy`:** Genera arte ASCII usando la herramienta `cowpy`

La configuración del workflow soporta proporcionar entradas y parámetros de manera flexible y reproducible.

### Habilidades adquiridas

A través de este curso práctico, ha aprendido cómo:

- Lanzar un workflow de Nextflow localmente
- Encontrar e interpretar salidas (resultados) y archivos de registro generados por Nextflow
- Reconocer los componentes principales de Nextflow que constituyen un workflow simple de múltiples pasos
- Describir conceptos del siguiente nivel como operadores y fábricas de channels
- Configurar pipelines para diferentes entornos de cómputo

Ahora está equipado con el conocimiento fundamental para comenzar a integrar pipelines de Nextflow existentes en su propio trabajo.

## Próximos pasos para desarrollar sus habilidades

Aquí están nuestras principales sugerencias de qué hacer a continuación:

- ¡No solo ejecute Nextflow, escríbalo! Conviértase en un desarrollador de Nextflow con [Hello Nextflow](../hello_nextflow/index.md)
- Aplique Nextflow a un caso de uso de análisis científico con [Nextflow for Science](../nf4_science/index.md)
- Comience con nf-core con [Hello nf-core](../hello_nf-core/index.md)
- Aprenda técnicas de solución de problemas con el [Debugging Side Quest](../side_quests/debugging.md)

Finalmente, le recomendamos que eche un vistazo a [**Seqera Platform**](https://seqera.io/), una plataforma basada en la nube desarrollada por los creadores de Nextflow que hace aún más fácil lanzar y gestionar sus workflows, así como gestionar sus datos y ejecutar análisis interactivamente en cualquier entorno.

## Obtener ayuda

Para recursos de ayuda y soporte de la comunidad, vea la [página de Ayuda](../help.md).

## Encuesta de retroalimentación

Antes de continuar, ¡tómese un minuto para completar la encuesta del curso! Su retroalimentación nos ayuda a mejorar nuestros materiales de entrenamiento para todos.

[Tomar la encuesta :material-arrow-right:](survey.md){ .md-button .md-button--primary }
