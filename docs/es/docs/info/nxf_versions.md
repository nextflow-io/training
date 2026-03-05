---
title: Versiones de Nextflow
description: Comprender y gestionar la evolución de las versiones de sintaxis de Nextflow
hide:
  - toc
  - footer
---

## Versión de sintaxis de Nextflow actualmente compatible y requisitos

A partir de la versión 3.0 del portal de capacitación, todos nuestros cursos de capacitación se basan en la versión 25.10.2 de Nextflow, a menos que se especifique lo contrario en la página de índice del curso (excepto los materiales obsoletos o archivados que pueden no incluir un aviso de versión).

Debido a que los cursos ahora utilizan entradas tipadas a nivel de workflow, así como directivas de salida a nivel de workflow, requieren el uso del analizador de sintaxis V2.
Si planeas usar el entorno que proporcionamos a través de [Github Codespaces](../envsetup/01_setup.md) o [devcontainers locales](../envsetup/03_devcontainer.md), no necesitas hacer nada a menos que se indique específicamente en las instrucciones del curso.
Sin embargo, si planeas trabajar en las capacitaciones en tu propio entorno ([Instalación manual](../envsetup/02_local.md)), deberás asegurarte de usar Nextflow versión 25.10.2 o posterior con el analizador de sintaxis v2 habilitado.

## Versiones anteriores de los materiales de capacitación

Nuestros materiales de capacitación han sido versionados desde febrero de 2025.

Puedes acceder a versiones anteriores de los materiales de capacitación que funcionan con versiones de Nextflow **anteriores a 25.10.2** a través del menú desplegable en la parte superior de cada página que muestra la versión numerada de los materiales de capacitación.
Cuando seleccionas una versión anterior de los materiales de capacitación, los enlaces al entorno de capacitación especificarán automáticamente la versión correspondiente del entorno.

## Otra información sobre las versiones de sintaxis de Nextflow

Nextflow tiene dos conceptos de versionado distintos que a veces se confunden: **versiones de DSL** y **versiones del analizador de sintaxis**.

**DSL1 vs DSL2** se refiere a formas fundamentalmente diferentes de escribir pipelines de Nextflow.
DSL1 fue la sintaxis original donde los procesos se conectaban implícitamente a través de canales.
DSL2, introducido en Nextflow 20.07, agregó características de modularidad: la capacidad de importar procesos y workflows desde otros archivos, bloques `workflow` explícitos y salidas de proceso nombradas.
DSL1 fue deprecado en Nextflow 22.03 y eliminado en 22.12.
Todo el código moderno de Nextflow usa DSL2.

**Analizador de sintaxis v1 vs v2** se refiere a diferentes analizadores que ambos funcionan con código DSL2.
El analizador v1 es el original, más permisivo.
El analizador v2 es más estricto y habilita nuevas características del lenguaje como tipado estático (entradas y salidas tipadas) y directivas de salida a nivel de workflow.
El analizador v2 también proporciona mejores mensajes de error y detecta más errores en tiempo de análisis en lugar de en tiempo de ejecución.
El analizador v2 se convertirá en el predeterminado en Nextflow 26.04.

En resumen: DSL2 es el lenguaje que escribes; la versión del analizador de sintaxis determina qué tan estrictamente se interpreta ese lenguaje y qué características avanzadas están disponibles.

### Verificar y configurar la versión de Nextflow

Puedes verificar qué versión de Nextflow está instalada en tu sistema usando el comando `nextflow --version`.

Para más información sobre cómo actualizar tu versión de Nextflow, consulta la documentación de referencia sobre [Actualización de Nextflow](https://www.nextflow.io/docs/latest/updating-nextflow.html).

### Habilitar el analizador de sintaxis v2

Para **habilitar** el analizador de sintaxis v2 para tu sesión actual, ejecuta el siguiente comando en tu terminal:

```bash
export NXF_SYNTAX_PARSER=v2
```

Para hacer esto permanente (mientras v2 se convierte en el predeterminado en Nextflow 26.04), agrega el comando export a tu perfil de shell (`~/.bashrc`, `~/.zshrc`, etc.):

```bash
echo 'export NXF_SYNTAX_PARSER=v2' >> ~/.bashrc
source ~/.bashrc
```

Ten en cuenta que la variable de entorno `NXF_SYNTAX_PARSER=v2` es un requisito temporal.
A partir de Nextflow 26.04 en adelante, el analizador v2 se convertirá en el predeterminado y esta configuración ya no será necesaria.

### Deshabilitar el analizador de sintaxis v2

Para **deshabilitar** el analizador de sintaxis v2 para tu sesión actual, ejecuta el siguiente comando en tu terminal:

```bash
export NXF_SYNTAX_PARSER=v1
```

<!-- Will it be possible to disable it in versions after 26.04? -->

### Migrar código existente

Para obtener orientación sobre la migración de código existente para cumplir con versiones más recientes de Nextflow, consulta las [Notas de Migración](https://www.nextflow.io/docs/latest/migrations/index.html) en la documentación de referencia.

Estos dos artículos son particularmente útiles para migrar a la versión más reciente:

- [Migración a salidas de workflow](https://www.nextflow.io/docs/latest/tutorials/workflow-outputs.html)
- [Migración a tipos estáticos](https://www.nextflow.io/docs/latest/tutorials/static-types.html)

Ambas características se cubren como parte de la capacitación para principiantes a partir de la versión 3.0 de los materiales de capacitación.

Dependiendo de la generación de código de Nextflow que pretendas migrar, es posible que puedas hacer la mayor parte usando el linter de Nextflow con el comando `nextflow lint -format`.
Consulta la referencia CLI para [`lint`](https://www.nextflow.io/docs/latest/reference/cli.html#lint) para más detalles.

Esperamos que esto sea útil.
Si necesitas ayuda, comunícate en Slack o en el foro.

---

<div markdown class="homepage_logos">

![Seqera](../assets/img/seqera_logo.png#only-light)

![Seqera](../assets/img/seqera_logo_dark.png#only-dark)

</div>
