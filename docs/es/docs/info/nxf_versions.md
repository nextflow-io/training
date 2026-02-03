---
title: Versiones de Nextflow
description: Comprensión y gestión de la evolución de las versiones de sintaxis de Nextflow
hide:
  - toc
  - footer
---

## Versión de sintaxis de Nextflow actualmente soportada y requisitos

A partir de la versión 3.0 del portal de entrenamiento, todos nuestros cursos de entrenamiento se basan en la versión 25.10.2 de Nextflow, a menos que se especifique lo contrario en la página de índice del curso (excepto materiales obsoletos o archivados que pueden no incluir un aviso de versión).

Debido a que los cursos ahora utilizan entradas tipadas a nivel de workflow así como directivas de salida a nivel de workflow, requieren el uso del parser de sintaxis V2.
Si planea utilizar el entorno que proporcionamos a través de [Github Codespaces](../envsetup/01_setup.md) o [devcontainers locales](../envsetup/03_devcontainer.md), no necesita hacer nada a menos que se indique específicamente en las instrucciones del curso.
Sin embargo, si planea trabajar en los entrenamientos en su propio entorno ([Instalación manual](../envsetup/02_local.md)), deberá asegurarse de usar Nextflow versión 25.10.2 o posterior con el parser de sintaxis v2 habilitado.

## Versiones anteriores de los materiales de entrenamiento

Nuestros materiales de entrenamiento han sido versionados desde febrero de 2025.

Puede acceder a versiones anteriores de los materiales de entrenamiento que funcionan con versiones de Nextflow **anteriores a 25.10.2** a través del elemento del menú desplegable en la parte superior de cada página que muestra la versión numerada de los materiales de entrenamiento.
Cuando seleccione una versión anterior de los materiales de entrenamiento, los enlaces al entorno de entrenamiento especificarán automáticamente la versión correspondiente del entorno.

## Otra información sobre las versiones de sintaxis de Nextflow

Nextflow tiene dos conceptos de versionado distintos que a veces se confunden: **versiones DSL** y **versiones del parser de sintaxis**.

**DSL1 vs DSL2** se refiere a formas fundamentalmente diferentes de escribir pipelines de Nextflow.
DSL1 era la sintaxis original donde los processes se conectaban implícitamente a través de channels.
DSL2, introducido en Nextflow 20.07, agregó características de modularidad: la capacidad de importar processes y workflows desde otros archivos, bloques `workflow` explícitos y salidas de process nombradas.
DSL1 fue declarado obsoleto en Nextflow 22.03 y eliminado en 22.12.
Todo el código moderno de Nextflow utiliza DSL2.

**Parser de sintaxis v1 vs v2** se refiere a diferentes parsers que ambos funcionan con código DSL2.
El parser v1 es el original, más permisivo.
El parser v2 es más estricto y habilita nuevas características del lenguaje como tipado estático (entradas y salidas tipadas) y directivas de salida a nivel de workflow.
El parser v2 también proporciona mejores mensajes de error y detecta más errores en tiempo de análisis en lugar de en tiempo de ejecución.
El parser v2 se convertirá en el predeterminado en Nextflow 26.04.

En resumen: DSL2 es el lenguaje que escribe; la versión del parser de sintaxis determina cuán estrictamente se interpreta ese lenguaje y qué características avanzadas están disponibles.

### Verificar y establecer la versión de Nextflow

Puede verificar qué versión de Nextflow está instalada en su sistema usando el comando `nextflow --version`.

Para más información sobre cómo actualizar su versión de Nextflow, consulte la documentación de referencia sobre [Updating Nextflow](https://www.nextflow.io/docs/latest/updating-nextflow.html).

### Habilitar el parser de sintaxis v2

Para **habilitar** el parser de sintaxis v2 para su sesión actual, ejecute el siguiente comando en su terminal:

```bash
export NXF_SYNTAX_PARSER=v2
```

Para hacer esto permanente (hasta que v2 se convierta en el predeterminado en Nextflow 26.04), agregue el comando export a su perfil de shell (`~/.bashrc`, `~/.zshrc`, etc.):

```bash
echo 'export NXF_SYNTAX_PARSER=v2' >> ~/.bashrc
source ~/.bashrc
```

Tenga en cuenta que la variable de entorno `NXF_SYNTAX_PARSER=v2` es un requisito temporal.
A partir de Nextflow 26.04, el parser v2 se convertirá en el predeterminado y esta configuración ya no será necesaria.

### Deshabilitar el parser de sintaxis v2

Para **deshabilitar** el parser de sintaxis v2 para su sesión actual, ejecute el siguiente comando en su terminal:

```bash
export NXF_SYNTAX_PARSER=v1
```

<!-- Will it be possible to disable it in versions after 26.04? -->

### Migrar código existente

Para obtener orientación sobre la migración de código existente para cumplir con versiones más recientes de Nextflow, consulte las [Migration Notes](https://www.nextflow.io/docs/latest/migrations/index.html) en la documentación de referencia.

Estos dos artículos son particularmente útiles para migrar a la versión más reciente:

- [Migrating to workflow outputs](https://www.nextflow.io/docs/latest/tutorials/workflow-outputs.html)
- [Migrating to static types](https://www.nextflow.io/docs/latest/tutorials/static-types.html)

Ambas características se cubren como parte del entrenamiento para principiantes a partir de la versión 3.0 de los materiales de entrenamiento.

Dependiendo de la generación de código Nextflow que pretenda migrar, es posible que pueda realizar la mayor parte del trabajo con el linter de Nextflow usando el comando `nextflow lint -format`.
Consulte la referencia de CLI para [`lint`](https://www.nextflow.io/docs/latest/reference/cli.html#lint) para más detalles.

Esperamos que esto sea útil.
Si necesita ayuda, comuníquese en Slack o en el foro.

---

<div markdown class="homepage_logos">

![Seqera](../assets/img/seqera_logo.png#only-light)

![Seqera](../assets/img/seqera_logo_dark.png#only-dark)

</div>
