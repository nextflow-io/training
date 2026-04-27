---

# Entorno de Desarrollo

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Traducción asistida por IA - [más información y sugerencias](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

Los Entornos de Desarrollo Integrados (IDEs) modernos pueden transformar radicalmente su experiencia de desarrollo con Nextflow. Esta misión secundaria se enfoca específicamente en aprovechar VS Code y su extensión de Nextflow para escribir código más rápido, detectar errores temprano y navegar workflows complejos de manera eficiente.

!!! note "Esto no es un tutorial tradicional"

    A diferencia de otros módulos de capacitación, esta guía está organizada como una colección de consejos rápidos, sugerencias y ejemplos prácticos en lugar de un tutorial paso a paso. Cada sección puede explorarse de forma independiente según sus intereses y necesidades de desarrollo actuales. Siéntase libre de saltar entre secciones y enfocarse en las características que sean más útiles para su desarrollo de workflows.

## Lo que debe saber primero

Esta guía asume que ha completado el curso de capacitación [Hello Nextflow](../hello_nextflow/) y que está familiarizado con los conceptos fundamentales de Nextflow, incluyendo:

- **Estructura básica del workflow**: Comprensión de procesos, workflows y cómo se conectan entre sí
- **Operaciones con canales**: Creación de canales, paso de datos entre procesos y uso de operadores básicos
- **Módulos y organización**: Creación de módulos reutilizables y uso de sentencias include
- **Conceptos básicos de configuración**: Uso de `nextflow.config` para parámetros, directivas de proceso y perfiles

## Lo que aprenderá aquí

Esta guía se enfoca en las **características de productividad del IDE** que lo convertirán en un desarrollador de Nextflow más eficiente:

- **Resaltado de sintaxis avanzado**: Comprensión de lo que VS Code le muestra sobre la estructura de su código
- **Autocompletado inteligente**: Aprovechamiento de sugerencias contextuales para escribir código más rápido
- **Detección de errores y diagnósticos**: Identificación de errores de sintaxis antes de ejecutar su workflow
- **Navegación de código**: Desplazamiento rápido entre procesos, módulos y definiciones
- **Formato y organización**: Mantenimiento de un estilo de código consistente y legible
- **Desarrollo asistido por IA** (opcional): Uso de herramientas modernas de IA integradas con su IDE

!!! info "¿Por qué las características del IDE ahora?"

    Probablemente ya ha estado usando VS Code durante el curso [Hello Nextflow](../hello_nextflow/), pero mantuvimos el enfoque en aprender los fundamentos de Nextflow en lugar de las características del IDE. Ahora que está familiarizado con los conceptos básicos de Nextflow como procesos, workflows, canales y módulos, está listo para aprovechar las sofisticadas características del IDE que lo convertirán en un desarrollador más eficiente.

    Piense en esto como "subir de nivel" su entorno de desarrollo: el mismo editor que ha estado usando tiene capacidades mucho más poderosas que se vuelven verdaderamente valiosas una vez que comprende para qué le están ayudando.

---

## 0. Configuración y Calentamiento

Configuremos un espacio de trabajo específicamente para explorar las características del IDE:

```bash title="Navigate to the IDE features directory"
cd side-quests/ide_features
```

Abra este directorio en VS Code:

```bash title="Open VS Code in current directory"
code .
```

El directorio `ide_features` contiene workflows de ejemplo que demuestran varias características del IDE:

```bash title="Show directory structure"
tree .
```

```console title="Project structure"
tree .
.
├── basic_workflow.nf
├── complex_workflow.nf
├── data
│   ├── sample_001.fastq.gz
│   ├── sample_002.fastq.gz
│   ├── sample_003.fastq.gz
│   ├── sample_004.fastq.gz
│   ├── sample_005.fastq.gz
│   └── sample_data.csv
├── modules
│   ├── fastqc.nf
│   ├── star.nf
│   └── utils.nf
└── nextflow.config

3 directories, 12 files
```

!!! note "Sobre los archivos de ejemplo"

    - `basic_workflow.nf` es un workflow básico funcional que puede ejecutar y modificar
    - `complex_workflow.nf` está diseñado solo para ilustración y demostrar características de navegación; puede que no se ejecute correctamente, pero muestra una estructura realista de workflow con múltiples archivos

### Atajos de teclado

Algunas de las características de esta guía utilizan atajos de teclado opcionales. Es posible que esté accediendo a este material a través de GitHub Codespaces en el navegador, y en ese caso algunos atajos pueden no funcionar como se espera porque se usan para otras cosas en su sistema.

Si está ejecutando VS Code localmente, como probablemente lo hará cuando esté escribiendo workflows, los atajos funcionarán como se describe.

Si usa una Mac, algunos (no todos) los atajos de teclado usarán "cmd" en lugar de "ctrl", y lo indicaremos en el texto como `Ctrl/Cmd`.

### 0.1. Instalación de la extensión de Nextflow

!!! note "¿Ya usa Devcontainers?"

    Si está trabajando en **GitHub Codespaces** o usando un **devcontainer local**, la extensión de Nextflow probablemente ya está instalada y configurada para usted. Puede omitir los pasos de instalación manual a continuación y proceder directamente a explorar las características de la extensión.

Para instalar la extensión manualmente:

1. Abra VS Code
2. Vaya a la vista de Extensiones haciendo clic en el ícono de extensiones a la izquierda: ![ícono de extensiones](img/extensions_icon.png) (atajo `Ctrl/Cmd+Shift+X` si está ejecutando VSCode localmente)
3. Busque "Nextflow"
4. Instale la extensión oficial de Nextflow

![Instalar la extensión de Nextflow](img/install_extension.png)

### 0.2. Diseño del espacio de trabajo

Como ha estado usando VS Code durante Hello Nextflow, ya está familiarizado con los conceptos básicos. A continuación se explica cómo organizar su espacio de trabajo de manera eficiente para esta sesión:

- **Área del editor**: Para ver y editar archivos. Puede dividirla en múltiples paneles para comparar archivos lado a lado.
- **Explorador de archivos** (![ícono del explorador de archivos](img/files_icon.png)) (`Ctrl/Cmd+Shift+E`): Los archivos y carpetas locales en su sistema. Manténgalo abierto a la izquierda para navegar entre archivos.
- **Terminal integrada** (`Ctrl+Shift+` acento grave, tanto en Windows como en MacOS): Una terminal para interactuar con la computadora en la parte inferior. Úsela para ejecutar Nextflow u otros comandos.
- **Panel de problemas** (`Ctrl+Shift+M`): VS Code mostrará aquí los errores y problemas que detecte. Esto es útil para identificar problemas de un vistazo.

Puede arrastrar paneles o ocultarlos (`Ctrl/Cmd+B` para alternar la barra lateral) para personalizar su diseño mientras trabajamos en los ejemplos.

### Conclusión

Tiene VS Code configurado con la extensión de Nextflow y comprende el diseño del espacio de trabajo para un desarrollo eficiente.

### ¿Qué sigue?

Aprenda cómo el resaltado de sintaxis le ayuda a comprender la estructura del código Nextflow de un vistazo.

---

## 1. Resaltado de Sintaxis y Estructura del Código

Ahora que su espacio de trabajo está configurado, exploremos cómo el resaltado de sintaxis de VS Code le ayuda a leer y escribir código Nextflow de manera más efectiva.

### 1.1. Elementos de sintaxis de Nextflow

Abra `basic_workflow.nf` para ver el resaltado de sintaxis en acción:

![Presentación de sintaxis](img/syntax_showcase.png)

Observe cómo VS Code resalta:

- **Palabras clave** (`process`, `workflow`, `input`, `output`, `script`) en colores distintos
- **Literales de cadena** y **parámetros** con diferentes estilos
- **Comentarios** en un color atenuado
- **Variables** y **llamadas a funciones** con el énfasis apropiado
- **Bloques de código** con guías de sangría adecuadas

!!! note "Colores dependientes del tema"

    Los colores específicos que vea dependerán de su tema de VS Code (modo oscuro/claro), la configuración de colores y las personalizaciones que haya realizado. Lo importante es que los diferentes elementos de sintaxis se distingan visualmente entre sí, lo que facilita la comprensión de la estructura del código independientemente del esquema de colores elegido.

### 1.2. Comprensión de la estructura del código

El resaltado de sintaxis le ayuda a identificar rápidamente:

- **Límites de procesos**: Distinción clara entre diferentes procesos
- **Bloques de entrada/salida**: Fácil identificación de las definiciones de flujo de datos
- **Bloques de script**: Los comandos reales que se están ejecutando
- **Operaciones con canales**: Pasos de transformación de datos
- **Directivas de configuración**: Configuraciones específicas del proceso

Esta organización visual se vuelve invaluable cuando se trabaja con workflows complejos que contienen múltiples procesos y flujos de datos intrincados.

### Conclusión

Comprende cómo el resaltado de sintaxis de VS Code le ayuda a leer la estructura del código Nextflow e identificar diferentes elementos del lenguaje para un desarrollo más rápido.

### ¿Qué sigue?

Aprenda cómo el autocompletado inteligente acelera la escritura de código con sugerencias contextuales.

---

## 2. Autocompletado Inteligente

Las características de autocompletado de VS Code le ayudan a escribir código más rápido y con menos errores al sugerir opciones apropiadas según el contexto.

### 2.1. Sugerencias contextuales

Las opciones de autocompletado varían según dónde se encuentre en su código:

#### Operaciones con canales

Abra `basic_workflow.nf` nuevamente e intente escribir `channel.` en el bloque del workflow:

![Autocompletado de canales](img/autocomplete_channel.png)

Verá sugerencias para:

- `fromPath()` - Crear canal a partir de rutas de archivos
- `fromFilePairs()` - Crear canal a partir de archivos emparejados
- `of()` - Crear canal a partir de valores
- `fromSRA()` - Crear canal a partir de accesiones SRA
- Y muchos más...

Esto le ayuda a encontrar rápidamente el factory de canal correcto sin necesidad de recordar los nombres exactos de los métodos.

También puede descubrir los operadores disponibles para aplicar a los canales. Por ejemplo, escriba `FASTQC.out.html.` para ver las operaciones disponibles:

![Autocompletado de operaciones de canal](img/autocomplete_operators.png)

#### Directivas de proceso

Dentro de un bloque de script de proceso, escriba `task.` para ver las propiedades de tiempo de ejecución disponibles:

![Autocompletado de propiedades de tarea](img/autocomplete_task.png)

#### Configuración

Abra nextflow.config y escriba `process.` en cualquier lugar para ver las directivas de proceso disponibles:

![Autocompletado de configuración](img/autocomplete_config.png)

Verá sugerencias para:

- `executor`
- `memory`
- `cpus`

Esto ahorra tiempo al configurar procesos y funciona en diferentes ámbitos de configuración. Por ejemplo, intente escribir `docker.` para ver las opciones de configuración específicas de Docker.

### Conclusión

Puede usar el autocompletado inteligente de VS Code para descubrir operaciones de canal disponibles, directivas de proceso y opciones de configuración sin necesidad de memorizar la sintaxis.

### ¿Qué sigue?

Aprenda cómo la detección de errores en tiempo real le ayuda a identificar problemas antes de ejecutar su workflow, simplemente leyendo el código.

## 3. Detección de Errores y Diagnósticos

La detección de errores en tiempo real de VS Code le ayuda a identificar problemas antes de ejecutar su workflow.

### 3.1. Detección de errores de sintaxis

Creemos un error deliberado para ver la detección en acción. Abra `basic_workflow.nf` y cambie el nombre del proceso de `FASTQC` a `FASTQ` (o cualquier otro nombre no válido). VS Code resaltará inmediatamente el error en el bloque del workflow con un subrayado ondulado rojo:

![Subrayado de error](img/error_underline.png)

### 3.2. Panel de problemas

Más allá del resaltado individual de errores, VS Code proporciona un panel de Problemas centralizado que agrega todos los errores, advertencias y mensajes informativos de su espacio de trabajo. Ábralo con `Ctrl/Cmd+Shift+M` y use el ícono de filtro para mostrar solo los errores relevantes al archivo actual:

![Filtrar el panel de problemas](img/active_file.png)

Haga clic en cualquier problema para ir directamente a la línea problemática:

![Panel de problemas](img/problems_panel.png)

Corrija el error cambiando el nombre del proceso de vuelta a `FASTQC`.

### 3.3. Patrones de error comunes

Los errores comunes en la sintaxis de Nextflow incluyen:

- **Corchetes faltantes**: `{` o `}` sin coincidencia
- **Bloques incompletos**: Secciones requeridas faltantes en los procesos
- **Sintaxis no válida**: DSL de Nextflow mal formado
- **Errores tipográficos en palabras clave**: Directivas de proceso mal escritas
- **Incompatibilidades de canales**: Incompatibilidades de tipos

El servidor de lenguaje de Nextflow resalta estos problemas en el panel de Problemas. Puede revisarlos temprano para evitar errores de sintaxis al ejecutar un pipeline.

### Conclusión

Puede usar la detección de errores de VS Code y el panel de Problemas para identificar errores de sintaxis y problemas antes de ejecutar su workflow, ahorrando tiempo y evitando frustraciones.

### ¿Qué sigue?

Aprenda cómo navegar eficientemente entre procesos, módulos y definiciones en workflows complejos.

---

## 4. Navegación de Código y Gestión de Símbolos

La navegación eficiente es crucial cuando se trabaja con workflows complejos que abarcan múltiples archivos. Para comprender esto, reemplace la definición del proceso en `basic_workflow.nf` con una importación del módulo que le hemos proporcionado:

=== "Después"

    ```groovy title="basic_workflow.nf" linenums="3"
    include { FASTQC } from './modules/fastqc.nf'
    ```

=== "Antes"

    ```groovy title="basic_workflow.nf" linenums="3"
    process FASTQC {
        tag "${sample_id}"
        publishDir "${params.output_dir}/fastqc", mode: 'copy'

        input:
        tuple val(sample_id), path(reads)

        output:
        tuple val(sample_id), path("*.html"), emit: html
        tuple val(sample_id), path("*.zip"), emit: zip

        script:
        def args = task.ext.args ?: ''
        """
        fastqc \\
            ${args} \\
            --threads ${task.cpus} \\
            ${reads}
        """
    }
    ```

### 4.1. Ir a la definición

Si pasa el cursor sobre el nombre de un proceso como `FASTQC`, verá una ventana emergente con la interfaz del módulo (entradas y salidas):

![Ir a la definición](img/syntax.png)

Esta característica es particularmente valiosa al crear workflows, ya que le permite comprender la interfaz del módulo sin abrir el archivo del módulo directamente.

Puede navegar rápidamente a cualquier definición de proceso, módulo o variable usando **Ctrl/Cmd+clic**. Pase el cursor sobre el enlace al archivo del módulo en la parte superior del script y siga el enlace como se sugiere:

![Seguir enlace](img/follow_link.png)

Lo mismo funciona para los nombres de procesos. Vuelva a `basic_workflow.nf` e inténtelo con el nombre del proceso `FASTQC` en el bloque del workflow. Esto lo lleva directamente al nombre del proceso (que es el mismo que el archivo del módulo en este ejemplo, pero podría estar a mitad de un archivo mucho más grande).

Para volver a donde estaba, use **Alt+←** (o **Ctrl+-** en Mac). Esta es una forma poderosa de explorar el código sin perder su lugar.

Ahora exploremos la navegación en un workflow más complejo usando `complex_workflow.nf` (el archivo solo de ilustración mencionado anteriormente). Este workflow contiene múltiples procesos definidos en archivos de módulo separados, así como algunos en línea. Aunque las estructuras complejas de múltiples archivos pueden ser difíciles de navegar manualmente, la capacidad de saltar a las definiciones hace que la exploración sea mucho más manejable.

1. Abra `complex_workflow.nf`
2. Navegue a las definiciones de módulos
3. Use **Alt+←** (o **Ctrl+-**) para navegar hacia atrás
4. Navegue al nombre del proceso `FASTQC` en el bloque del workflow. Esto lo lleva directamente al nombre del proceso (que es el mismo que el archivo del módulo en este ejemplo, pero podría estar a mitad de un archivo mucho más grande).
5. Navegue hacia atrás nuevamente
6. Navegue al proceso `TRIM_GALORE` en el bloque del workflow. Este está definido en línea, por lo que no lo llevará a un archivo separado, pero aún le mostrará la definición del proceso y podrá navegar de vuelta a donde estaba.

### 4.2. Navegación de símbolos

Con `complex_workflow.nf` aún abierto, puede obtener una descripción general de todos los símbolos en el archivo escribiendo `@` en la barra de búsqueda en la parte superior de VSCode (el atajo de teclado es `Ctrl/Cmd+Shift+O`, pero puede que no funcione en Codespaces). Esto abre el panel de navegación de símbolos, que lista todos los símbolos en el archivo actual:

![Navegación de símbolos](img/symbols.png)

Esto muestra:

- Todas las definiciones de procesos
- Definiciones de workflows (hay dos workflows definidos en este archivo)
- Definiciones de funciones

Comience a escribir para filtrar los resultados.

### 4.3. Buscar todas las referencias

Comprender dónde se usa un proceso o variable en toda su base de código puede ser muy útil. Por ejemplo, si desea encontrar todas las referencias al proceso `FASTQC`, comience navegando a su definición. Puede hacerlo abriendo `modules/fastqc.nf` directamente, o usando la característica de navegación rápida de VS Code con `Ctrl/Cmd+clic` como hicimos anteriormente. Una vez en la definición del proceso, haga clic derecho en el nombre del proceso `FASTQC` y seleccione "Find All References" del menú contextual para ver todas las instancias donde se usa.

![Buscar referencias](img/references.png)

Esta característica muestra todas las instancias donde se hace referencia a `FASTQC` dentro de su espacio de trabajo, incluyendo su uso en los dos workflows distintos. Esta información es crucial para evaluar el impacto potencial de las modificaciones al proceso `FASTQC`.

### 4.4. Panel de esquema

El panel de Esquema, ubicado en la barra lateral del Explorador (haga clic en ![ícono del Explorador](img/files_icon.png)), proporciona una descripción general conveniente de todos los símbolos en su archivo actual. Esta característica le permite navegar y gestionar rápidamente la estructura de su código mostrando funciones, variables y otros elementos clave en una vista jerárquica.

![Panel de esquema](img/outline.png)

Use el panel de Esquema para navegar rápidamente a diferentes partes de su código sin usar el explorador de archivos.

### 4.5. Visualización del DAG

La extensión de Nextflow para VS Code puede visualizar su workflow como un Grafo Acíclico Dirigido (DAG). Esto le ayuda a comprender el flujo de datos y las dependencias entre procesos. Abra `complex_workflow.nf` y haga clic en el botón "Preview DAG" sobre `workflow {` (el segundo bloque `workflow` en este archivo):

![Vista previa del DAG](img/dag_preview.png)

Este es solo el workflow de 'entrada', pero también puede previsualizar el DAG para los workflows internos haciendo clic en el botón "Preview DAG" sobre el workflow `RNASEQ_PIPELINE {` más arriba:

![Vista previa del DAG del workflow interno](img/dag_preview_inner.png)

Para este workflow, puede usar los nodos del DAG para navegar a las definiciones de proceso correspondientes en el código. Haga clic en un nodo y lo llevará a la definición del proceso relevante en el editor. Especialmente cuando un workflow crece a un tamaño grande, esto puede ayudarle realmente a navegar por el código y comprender cómo están conectados los procesos.

### Conclusión

Puede navegar workflows complejos de manera eficiente usando ir a la definición, búsqueda de símbolos, buscar referencias y visualización del DAG para comprender la estructura del código y las dependencias.

### ¿Qué sigue?

Aprenda cómo trabajar eficazmente con múltiples archivos interconectados en proyectos Nextflow más grandes.

## 5. Trabajo con Múltiples Archivos

El desarrollo real con Nextflow implica trabajar con múltiples archivos interconectados. Exploremos cómo VS Code le ayuda a gestionar proyectos complejos de manera eficiente.

### 5.1. Navegación rápida de archivos

Con `complex_workflow.nf` abierto, notará que importa varios módulos. Practiquemos la navegación rápida entre ellos.

Presione **Ctrl+P** (o **Cmd+P**) y comience a escribir "fast":

VS Code le mostrará los archivos coincidentes. Seleccione `modules/fastqc.nf` para ir allí instantáneamente. Esto es mucho más rápido que hacer clic en el explorador de archivos cuando sabe aproximadamente qué archivo está buscando.

Pruebe esto con otros patrones:

- Escriba "star" para encontrar el archivo del módulo de alineamiento STAR (`star.nf`)
- Escriba "utils" para encontrar el archivo de funciones utilitarias (`utils.nf`)
- Escriba "config" para ir a los archivos de configuración (`nextflow.config`)

### 5.2. Editor dividido para desarrollo con múltiples archivos

Cuando trabaja con módulos, a menudo necesita ver tanto el workflow principal como las definiciones de módulos simultáneamente. Configuremos esto:

1. Abra `complex_workflow.nf`
2. Abra `modules/fastqc.nf` en una nueva pestaña
3. Haga clic derecho en la pestaña `modules/fastqc.nf` y seleccione "Split Right"
4. Ahora puede ver ambos archivos lado a lado

![Editor dividido](img/split_editor.png)

Esto es invaluable cuando:

- Verifica las interfaces de módulos mientras escribe llamadas al workflow y la vista previa no es suficiente
- Compara procesos similares en diferentes módulos
- Depura el flujo de datos entre el workflow y los módulos

### 5.3. Búsqueda en todo el proyecto

A veces necesita encontrar dónde se usan patrones específicos en todo su proyecto. Presione `Ctrl/Cmd+Shift+F` para abrir el panel de búsqueda.

Intente buscar `publishDir` en todo el espacio de trabajo:

![Búsqueda en el proyecto](img/project_search.png)

Esto le muestra cada archivo que usa directorios de publicación, ayudándole a:

- Comprender los patrones de organización de salidas
- Encontrar ejemplos de directivas específicas
- Garantizar la consistencia entre módulos

### Conclusión

Puede gestionar proyectos complejos con múltiples archivos usando navegación rápida de archivos, editores divididos y búsqueda en todo el proyecto para trabajar eficientemente entre workflows y módulos.

### ¿Qué sigue?

Aprenda cómo las características de formato y mantenimiento de código mantienen sus workflows organizados y legibles.

---

## 6. Formato y Mantenimiento del Código

El formato adecuado del código es esencial no solo por estética, sino también para mejorar la legibilidad, la comprensión y la facilidad de actualización de workflows complejos.

### 6.1. Formato automático en acción

Abra `basic_workflow.nf` y desordene deliberadamente el formato:

- Elimine algo de sangría: Seleccione todo el documento y presione `shift+tab` muchas veces para eliminar tantas sangrías como sea posible.
- Agregue espacios adicionales en lugares aleatorios: en la sentencia `channel.fromPath`, agregue 30 espacios después del `(`.
- Rompa algunas líneas de manera incómoda: Agregue una nueva línea entre el operador `.view {` y la cadena `Processing sample:` pero no agregue una nueva línea correspondiente antes del paréntesis de cierre `}`.

Ahora presione `Shift+Alt+F` (o `Shift+Option+F` en MacOS) para formatear automáticamente:

VS Code inmediatamente:

- Corrige la sangría para mostrar claramente la estructura del proceso
- Alinea elementos similares de manera consistente
- Elimina espacios en blanco innecesarios
- Mantiene saltos de línea legibles

Tenga en cuenta que el formato automático puede no resolver todos los problemas de estilo de código. El servidor de lenguaje de Nextflow intenta mantener su código ordenado, pero también respeta sus preferencias personales en ciertas áreas. Por ejemplo, si elimina la sangría dentro del bloque `script` de un proceso, el formateador lo dejará como está, ya que podría preferir intencionalmente ese estilo.

Actualmente, no existe una aplicación estricta de estilo para Nextflow, por lo que el servidor de lenguaje ofrece cierta flexibilidad. Sin embargo, aplicará consistentemente reglas de formato alrededor de las definiciones de métodos y funciones para mantener la claridad.

### 6.2. Características de organización del código

#### Comentado rápido

Seleccione un bloque de código en su workflow y presione **Ctrl+/** (o **Cmd+/**) para comentarlo:

```groovy
// workflow {
//     ch_input = channel.fromPath(params.input)
//         .splitCsv(header: true)
//         .map { row -> [row.sample_id, file(row.fastq_path)] }
//
//     FASTQC(ch_input)
// }
```

Esto es perfecto para:

- Deshabilitar temporalmente partes de workflows durante el desarrollo
- Agregar comentarios explicativos a operaciones de canal complejas
- Documentar secciones del workflow

Use **Ctrl+/** (o **Cmd+/**) nuevamente para descomentar el código.

#### Plegado de código para una vista general

En `complex_workflow.nf`, observe las pequeñas flechas junto a las definiciones de procesos. Haga clic en ellas para plegar (colapsar) los procesos:

![Plegado de código](img/code_folding.png)

Esto le da una vista general de alto nivel de la estructura de su workflow sin perderse en los detalles de implementación.

#### Coincidencia de corchetes

Coloque el cursor junto a cualquier corchete `{` o `}` y VS Code resalta el corchete coincidente. Use **Ctrl+Shift+\\** (o **Cmd+Shift+\\**) para saltar entre corchetes coincidentes.

Esto es crucial para:

- Comprender los límites de los procesos
- Encontrar corchetes faltantes o adicionales
- Navegar estructuras de workflow anidadas

#### Selección y edición multilínea

Para editar múltiples líneas simultáneamente, VS Code ofrece potentes capacidades de múltiples cursores:

- **Selección multilínea**: Mantenga presionado **Ctrl+Alt** (o **Cmd+Option** en MacOS) y use las teclas de flecha para seleccionar múltiples líneas
- **Sangría multilínea**: Seleccione múltiples líneas y use **Tab** para agregar sangría o **Shift+Tab** para quitarla en bloques completos

Esto es particularmente útil para:

- Aplicar sangría consistente a bloques de procesos completos
- Agregar comentarios a múltiples líneas a la vez
- Editar definiciones de parámetros similares en múltiples procesos

### Conclusión

Puede mantener un código limpio y legible usando formato automático, características de comentado, plegado de código, coincidencia de corchetes y edición multilínea para organizar workflows complejos de manera eficiente.

### ¿Qué sigue?

Aprenda cómo VS Code se integra con su flujo de trabajo de desarrollo más amplio más allá de la simple edición de código.

---

## 7. Integración con el Flujo de Trabajo de Desarrollo

VS Code se integra bien con su flujo de trabajo de desarrollo más allá de la simple edición de código.

### 7.1. Integración con control de versiones

!!! note "Codespaces e integración con Git"

    Si está trabajando en **GitHub Codespaces**, es posible que algunas características de integración con Git no funcionen como se espera, particularmente los atajos de teclado para el Control de código fuente. También es posible que haya rechazado abrir el directorio como repositorio Git durante la configuración inicial, lo cual está bien para propósitos de capacitación.

Si su proyecto es un repositorio git (como este), VS Code muestra:

- Archivos modificados con indicadores de color
- Estado de Git en la barra de estado
- Vistas de diferencias en línea
- Capacidades de confirmación y envío

Abra el panel de Control de código fuente usando el botón de control de código fuente (![ícono de control de código fuente](img/source_control_icon.png)) (`Ctrl+Shift+G` o `Cmd+Shift+G` si está trabajando con VSCode localmente) para ver los cambios de git y confirmar directamente en el editor.

![Panel de control de código fuente](img/source_control.png)

### 7.2. Ejecución e inspección de workflows

Ejecutemos un workflow y luego inspeccionemos los resultados. En la terminal integrada (`Ctrl+Shift+` acento grave, tanto en Windows como en MacOS), ejecute el workflow básico:

```bash title="Run the basic workflow"
nextflow run basic_workflow.nf --input data/sample_data.csv --output_dir results
```

Mientras el workflow se ejecuta, verá la salida en tiempo real en la terminal. Después de completarse, puede usar VS Code para inspeccionar los resultados sin salir del editor:

1. **Navegar a los directorios de trabajo**: Use el explorador de archivos o la terminal para explorar `.nextflow/work`
2. **Abrir archivos de registro**: Haga clic en las rutas de archivos de registro en la salida de la terminal para abrirlos directamente en VS Code
3. **Inspeccionar salidas**: Explore los directorios de resultados publicados en el explorador de archivos
4. **Ver informes de ejecución**: Abra informes HTML directamente en VS Code o en su navegador

Esto mantiene todo en un solo lugar en lugar de cambiar entre múltiples aplicaciones.

### Conclusión

Puede integrar VS Code con el control de versiones y la ejecución de workflows para gestionar todo su proceso de desarrollo desde una única interfaz.

### ¿Qué sigue?

Vea cómo todas estas características del IDE funcionan juntas en su flujo de trabajo de desarrollo diario.

---

## 8. Resumen y notas rápidas

Aquí hay algunas notas rápidas sobre cada una de las características del IDE discutidas anteriormente:

### 8.1. Comenzar una nueva característica

1. **Apertura rápida de archivos** (`Ctrl+P` o `Cmd+P`) para encontrar módulos existentes relevantes
2. **Editor dividido** para ver procesos similares lado a lado
3. **Navegación de símbolos** (`Ctrl+Shift+O` o `Cmd+Shift+O`) para comprender la estructura del archivo
4. **Autocompletado** para escribir nuevo código rápidamente

### 8.2. Depuración de problemas

1. **Panel de problemas** (`Ctrl+Shift+M` o `Cmd+Shift+M`) para ver todos los errores a la vez
2. **Ir a la definición** (`Ctrl+clic` o `Cmd+clic`) para comprender las interfaces de los procesos
3. **Buscar todas las referencias** para ver cómo se usan los procesos
4. **Búsqueda en todo el proyecto** para encontrar patrones o problemas similares

### 8.3. Refactorización y mejora

1. **Búsqueda en todo el proyecto** (`Ctrl+Shift+F` o `Cmd+Shift+F`) para encontrar patrones
2. **Formato automático** (`Shift+Alt+F` o `Shift+Option+F`) para mantener la consistencia
3. **Plegado de código** para enfocarse en la estructura
4. **Integración con Git** para rastrear cambios

---

## Resumen

Ha completado un recorrido rápido por las características del IDE de VS Code para el desarrollo con Nextflow. Estas herramientas lo harán significativamente más productivo al:

- **Reducir errores** mediante la verificación de sintaxis en tiempo real
- **Acelerar el desarrollo** con autocompletado inteligente
- **Mejorar la navegación** en workflows complejos con múltiples archivos
- **Mantener la calidad** mediante un formato consistente
- **Mejorar la comprensión** mediante resaltado avanzado y visualización de la estructura

No esperamos que recuerde todo, pero ahora que sabe que estas características existen, podrá encontrarlas cuando las necesite. A medida que continúe desarrollando workflows con Nextflow, estas características del IDE se volverán algo natural, permitiéndole enfocarse en escribir código de alta calidad en lugar de lidiar con la sintaxis y la estructura.

### ¿Qué sigue?

Aplique estas habilidades del IDE mientras trabaja en otros módulos de capacitación, por ejemplo:

- **[nf-test](nf-test.md)**: Cree suites de pruebas completas para sus workflows
- **[Hello nf-core](../../hello_nf-core/)**: Construya pipelines de calidad de producción con estándares de la comunidad

El verdadero poder de estas características del IDE emerge cuando trabaja en proyectos más grandes y complejos. Comience a incorporarlas en su flujo de trabajo gradualmente: en pocas sesiones, se volverán algo natural y transformarán la manera en que aborda el desarrollo con Nextflow.

Desde detectar errores antes de que lo ralenticen hasta navegar bases de código complejas con facilidad, estas herramientas lo convertirán en un desarrollador más seguro y eficiente.

¡Feliz programación!
