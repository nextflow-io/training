# Entorno de Desarrollo

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Traducción asistida por IA - [más información y sugerencias](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

Los Entornos de Desarrollo Integrados (IDEs) modernos pueden transformar dramáticamente su experiencia de desarrollo con Nextflow. Esta misión secundaria se enfoca específicamente en aprovechar VS Code y su extensión de Nextflow para escribir código más rápido, detectar errores temprano y navegar workflows complejos de manera eficiente.

!!! note "Este no es un tutorial tradicional"

    A diferencia de otros módulos de entrenamiento, esta guía está organizada como una colección de sugerencias rápidas, consejos y ejemplos prácticos en lugar de un tutorial paso a paso. Cada sección puede explorarse de forma independiente según sus intereses y necesidades actuales de desarrollo. Siéntase libre de moverse entre secciones y enfocarse en las características que serán más útiles inmediatamente para el desarrollo de su workflow.

## Qué debería saber primero

Esta guía asume que ha completado el curso de entrenamiento [Hello Nextflow](../hello_nextflow/) y se siente cómodo con conceptos fundamentales de Nextflow incluyendo:

- **Estructura básica de workflow**: Comprensión de procesos, workflows y cómo se conectan entre sí
- **Operaciones de canal**: Creación de canales, paso de datos entre procesos y uso de operadores básicos
- **Módulos y organización**: Creación de módulos reutilizables y uso de declaraciones include
- **Conceptos básicos de configuración**: Uso de `nextflow.config` para parámetros, directivas de proceso y perfiles

## Qué aprenderá aquí

Esta guía se enfoca en **características de productividad del IDE** que lo convertirán en un desarrollador de Nextflow más eficiente:

- **Resaltado de sintaxis avanzado**: Comprender qué está mostrando VS Code sobre la estructura de su código
- **Auto-completado inteligente**: Aprovechar sugerencias contextuales para escribir código más rápido
- **Detección de errores y diagnósticos**: Detectar errores de sintaxis antes de ejecutar su workflow
- **Navegación de código**: Moverse rápidamente entre procesos, módulos y definiciones
- **Formato y organización**: Mantener un estilo de código consistente y legible
- **Desarrollo asistido por IA** (opcional): Usar herramientas de IA modernas integradas con su IDE

!!! info "¿Por qué características del IDE ahora?"

    Es probable que ya haya estado usando VS Code durante el curso [Hello Nextflow](../hello_nextflow/), pero mantuvimos el enfoque en aprender los fundamentos de Nextflow en lugar de características del IDE. Ahora que se siente cómodo con conceptos básicos de Nextflow como procesos, workflows, canales y módulos, está listo para aprovechar las características sofisticadas del IDE que lo harán un desarrollador más eficiente.

    Piense en esto como "subir de nivel" su entorno de desarrollo - el mismo editor que ha estado usando tiene capacidades mucho más poderosas que se vuelven verdaderamente valiosas una vez que comprende en qué lo están ayudando.

---

## 0. Configuración y Preparación

Configuremos un espacio de trabajo específicamente para explorar características del IDE:

```bash title="Navegue al directorio de características del IDE"
cd side-quests/ide_features
```

Abra este directorio en VS Code:

```bash title="Abra VS Code en el directorio actual"
code .
```

El directorio `ide_features` contiene workflows de ejemplo que demuestran varias características del IDE:

```bash title="Mostrar estructura del directorio"
tree .
```

```console title="Estructura del proyecto"
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

!!! note "Acerca de los Archivos de Ejemplo"

    - `basic_workflow.nf` es un workflow básico funcional que puede ejecutar y modificar
    - `complex_workflow.nf` está diseñado solo para ilustración para demostrar características de navegación - puede que no se ejecute exitosamente pero muestra una estructura de workflow multi-archivo realista

### Atajos de Teclado

Algunas de las características en esta guía usarán atajos de teclado opcionales. Es muy posible que esté accediendo a este material vía GitHub Codespaces en el navegador, y en este caso a veces los atajos no funcionarán como se espera porque se usan para otras cosas en su sistema.

Si está ejecutando VS Code localmente, como probablemente lo hará cuando realmente esté escribiendo workflows, los atajos funcionarán como se describe.

Si está usando Mac, algunos (no todos) atajos de teclado usarán "cmd" en lugar de "ctrl", y lo indicaremos en el texto como `Ctrl/Cmd`.

### 0.1. Instalación de la Extensión de Nextflow

!!! note "¿Ya Usa Devcontainers?"

    Si está trabajando en **GitHub Codespaces** o usando un **devcontainer local**, la extensión de Nextflow probablemente ya está instalada y configurada para usted. Puede omitir los pasos de instalación manual a continuación y proceder directamente a explorar las características de la extensión.

Para instalar la extensión manualmente:

1. Abra VS Code
2. Vaya a la vista de Extensiones haciendo clic en el ícono de extensiones a la izquierda: ![ícono de extensiones](img/extensions_icon.png) (atajo `Ctrl/Cmd+Shift+X` si está ejecutando VS Code localmente)
3. Busque "Nextflow"
4. Instale la extensión oficial de Nextflow

![Instalar Extensión de Nextflow](img/install_extension.png)

### 0.2. Diseño del Espacio de Trabajo

Dado que ha estado usando VS Code durante Hello Nextflow, ya está familiarizado con lo básico. Aquí está cómo organizar su espacio de trabajo eficientemente para esta sesión:

- **Área del Editor**: Para ver y editar archivos. Puede dividir esto en múltiples paneles para comparar archivos lado a lado.
- **Explorador de Archivos** haga clic (![ícono del explorador de archivos](img/files_icon.png)) (`Ctrl/Cmd+Shift+E`): Los archivos y carpetas locales en su sistema. Mantenga esto abierto a la izquierda para navegar entre archivos
- **Terminal Integrada** (`Ctrl+Shift+` backtick tanto para Windows como MacOS): Una terminal para interactuar con la computadora en la parte inferior. Use esto para ejecutar Nextflow u otros comandos.
- **Panel de Problemas** (`Ctrl+Shift+M`): VS Code mostrará cualquier error y problema que detecte aquí. Esto es útil para resaltar problemas de un vistazo.

Puede arrastrar paneles o ocultarlos (`Ctrl/Cmd+B` para alternar la barra lateral) para personalizar su diseño mientras trabajamos con los ejemplos.

### Conclusión

Tiene VS Code configurado con la extensión de Nextflow y comprende el diseño del espacio de trabajo para desarrollo eficiente.

### ¿Qué sigue?

Aprenda cómo el resaltado de sintaxis le ayuda a comprender la estructura del código Nextflow de un vistazo.

---

## 1. Resaltado de Sintaxis y Estructura del Código

Ahora que su espacio de trabajo está configurado, exploremos cómo el resaltado de sintaxis de VS Code le ayuda a leer y escribir código Nextflow más efectivamente.

### 1.1. Elementos de Sintaxis de Nextflow

Abra `basic_workflow.nf` para ver el resaltado de sintaxis en acción:

![Demostración de Sintaxis](img/syntax_showcase.png)

Note cómo VS Code resalta:

- **Palabras clave** (`process`, `workflow`, `input`, `output`, `script`) en colores distintos
- **Literales de cadena** y **parámetros** con estilo diferente
- **Comentarios** en un color atenuado
- **Variables** y **llamadas a funciones** con énfasis apropiado
- **Bloques de código** con guías de indentación apropiadas

!!! note "Colores Dependientes del Tema"

    Los colores específicos que vea dependerán de su tema de VS Code (modo oscuro/claro), configuración de colores y cualquier personalización que haya hecho. Lo importante es que diferentes elementos de sintaxis se distingan visualmente entre sí, haciendo la estructura del código más fácil de entender independientemente de su esquema de color elegido.

### 1.2. Comprensión de la Estructura del Código

El resaltado de sintaxis le ayuda a identificar rápidamente:

- **Límites de proceso**: Distinción clara entre diferentes procesos
- **Bloques de entrada/salida**: Fácil de identificar definiciones de flujo de datos
- **Bloques de script**: Los comandos reales que se están ejecutando
- **Operaciones de canal**: Pasos de transformación de datos
- **Directivas de configuración**: Configuraciones específicas del proceso

Esta organización visual se vuelve invaluable cuando se trabaja con workflows complejos que contienen múltiples procesos y flujos de datos intrincados.

### Conclusión

Comprende cómo el resaltado de sintaxis de VS Code le ayuda a leer la estructura del código Nextflow e identificar diferentes elementos del lenguaje para un desarrollo más rápido.

### ¿Qué sigue?

Aprenda cómo el auto-completado inteligente acelera la escritura de código con sugerencias contextuales.

---

## 2. Auto-completado Inteligente

Las características de auto-completado de VS Code le ayudan a escribir código más rápido y con menos errores sugiriendo opciones apropiadas según el contexto.

### 2.1. Sugerencias Contextuales

Las opciones de auto-completado varían dependiendo de dónde esté en su código:

#### Operaciones de Canal

Abra `basic_workflow.nf` nuevamente e intente escribir `channel.` en el bloque workflow:

![Auto-completado de canal](img/autocomplete_channel.png)

Verá sugerencias para:

- `fromPath()` - Crear canal desde rutas de archivo
- `fromFilePairs()` - Crear canal desde archivos emparejados
- `of()` - Crear canal desde valores
- `fromSRA()` - Crear canal desde accesos SRA
- Y muchos más...

Esto le ayuda a encontrar rápidamente la fábrica de canal correcta sin necesidad de recordar nombres exactos de métodos.

También puede descubrir los operadores disponibles para aplicar a canales. Por ejemplo, escriba `FASTQC.out.html.` para ver operaciones disponibles:

![Auto-completado de operaciones de canal](img/autocomplete_operators.png)

#### Directivas de Proceso

Dentro de un bloque script de proceso, escriba `task.` para ver propiedades de tiempo de ejecución disponibles:

![Auto-completado de propiedades de tarea](img/autocomplete_task.png)

#### Configuración

Abra nextflow.config y escriba `process.` en cualquier lugar para ver directivas de proceso disponibles:

![Auto-completado de configuración](img/autocomplete_config.png)

Verá sugerencias para:

- `executor`
- `memory`
- `cpus`

Esto ahorra tiempo al configurar procesos y funciona a través de diferentes ámbitos de configuración. Por ejemplo, intente escribir `docker.` para ver opciones de configuración específicas de Docker.

### Conclusión

Puede usar el auto-completado inteligente de VS Code para descubrir operaciones de canal disponibles, directivas de proceso y opciones de configuración sin memorizar la sintaxis.

### ¿Qué sigue?

Aprenda cómo la detección de errores en tiempo real le ayuda a detectar problemas antes de ejecutar su workflow, simplemente leyendo el código.

## 3. Detección de Errores y Diagnósticos

La detección de errores en tiempo real de VS Code le ayuda a detectar problemas antes de ejecutar su workflow.

### 3.1. Detección de Errores de Sintaxis

Creemos un error deliberado para ver la detección en acción. Abra `basic_workflow.nf` y cambie el nombre del proceso de `FASTQC` a `FASTQ` (o cualquier otro nombre inválido). VS Code inmediatamente resaltará el error en el bloque workflow con un subrayado ondulado rojo:

![Subrayado de error](img/error_underline.png)

### 3.2. Panel de Problemas

Más allá del resaltado de errores individuales, VS Code proporciona un panel de Problemas centralizado que agrega todos los errores, advertencias y mensajes de información en su espacio de trabajo. Ábralo con `Ctrl/Cmd+Shift+M` y use el ícono de filtro para mostrar solo errores relevantes al archivo actual:

![Filtrar el panel de problemas](img/active_file.png)

Haga clic en cualquier problema para saltar directamente a la línea problemática

![Panel de Problemas](img/problems_panel.png)

Corrija el error cambiando el nombre del proceso de vuelta a `FASTQC`.

### 3.3. Patrones Comunes de Error

Los errores comunes en la sintaxis de Nextflow incluyen:

- **Corchetes faltantes**: `{` o `}` sin coincidencia
- **Bloques incompletos**: Secciones requeridas faltantes en procesos
- **Sintaxis inválida**: DSL de Nextflow mal formado
- **Errores tipográficos en palabras clave**: Directivas de proceso mal escritas
- **Desajustes de canal**: Incompatibilidades de tipo

El servidor de lenguaje de Nextflow resalta estos problemas en el panel de Problemas. Puede revisar estos temprano para evitar errores de sintaxis al ejecutar un pipeline.

### Conclusión

Puede usar la detección de errores de VS Code y el panel de Problemas para detectar errores de sintaxis y problemas antes de ejecutar su workflow, ahorrando tiempo y previniendo frustraciones.

### ¿Qué sigue?

Aprenda cómo navegar eficientemente entre procesos, módulos y definiciones en workflows complejos.

---

## 4. Navegación de Código y Gestión de Símbolos

La navegación eficiente es crucial cuando se trabaja con workflows complejos que abarcan múltiples archivos. Para entender esto, reemplace la definición del proceso en `basic_workflow.nf` con una importación para el módulo que le hemos proporcionado:

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

### 4.1. Ir a Definición

Si pasa el mouse sobre un nombre de proceso como `FASTQC`, verá una ventana emergente con la interfaz del módulo (entradas y salidas):

![Ir a definición](img/syntax.png)

Esta característica es particularmente valiosa cuando se crean workflows, ya que le permite comprender la interfaz del módulo sin abrir el archivo del módulo directamente.

Puede navegar rápidamente a cualquier definición de proceso, módulo o variable usando **Ctrl/Cmd-clic**. Pase el mouse sobre el enlace al archivo del módulo en la parte superior del script y siga el enlace como se sugiere:

![Seguir enlace](img/follow_link.png)

Lo mismo funciona para nombres de proceso. Regrese a `basic_workflow.nf` e intente esto en el nombre del proceso `FASTQC` en el bloque workflow. Esto lo enlaza directamente al nombre del proceso (que es el mismo que el archivo del módulo en este ejemplo, pero podría estar a mitad de camino en un archivo mucho más grande).

Para regresar a donde estaba, use **Alt+←** (o **Ctrl+-** en Mac). Esta es una forma poderosa de explorar código sin perder su lugar.

Ahora exploremos la navegación en un workflow más complejo usando `complex_workflow.nf` (el archivo solo para ilustración mencionado anteriormente). Este workflow contiene múltiples procesos definidos en archivos de módulos separados, así como algunos en línea. Aunque las estructuras complejas multi-archivo pueden ser desafiantes de navegar manualmente, la capacidad de saltar a definiciones hace la exploración mucho más manejable.

1. Abra `complex_workflow.nf`
2. Navegue a definiciones de módulos
3. Use **Alt+←** (o **Ctrl+-**) para navegar hacia atrás
4. Navegue al nombre del proceso `FASTQC` en el bloque workflow. Esto lo enlaza directamente al nombre del proceso (que es el mismo que el archivo del módulo en este ejemplo, pero podría estar a mitad de camino en un archivo mucho más grande).
5. Navegue hacia atrás nuevamente
6. Navegue al proceso `TRIM_GALORE` en el bloque workflow. Este está definido en línea, así que no lo llevará a un archivo separado, pero aún le mostrará la definición del proceso, y aún puede navegar de regreso a donde estaba.

### 4.2. Navegación de Símbolos

Con `complex_workflow.nf` aún abierto, puede obtener una visión general de todos los símbolos en el archivo escribiendo `@` en la barra de búsqueda en la parte superior de VS Code (el atajo de teclado es `Ctrl/Cmd+Shift+O`, pero puede no funcionar en Codespaces). Esto abre el panel de navegación de símbolos, que lista todos los símbolos en el archivo actual:

![Navegación de símbolos](img/symbols.png)

Esto muestra:

- Todas las definiciones de proceso
- Definiciones de workflow (hay dos workflows definidos en este archivo)
- Definiciones de función

Comience a escribir para filtrar resultados.

### 4.3. Encontrar Todas las Referencias

Comprender dónde se usa un proceso o variable en toda su base de código puede ser muy útil. Por ejemplo, si desea encontrar todas las referencias al proceso `FASTQC`, comience navegando a su definición. Puede hacer esto abriendo `modules/fastqc.nf` directamente, o usando la característica de navegación rápida de VS Code con `Ctrl/Cmd-clic` como hicimos arriba. Una vez en la definición del proceso, haga clic derecho en el nombre del proceso `FASTQC` y seleccione "Find All References" del menú contextual para ver todas las instancias donde se usa.

![Encontrar referencias](img/references.png)

Esta característica muestra todas las instancias donde se hace referencia a `FASTQC` dentro de su espacio de trabajo, incluyendo su uso en los dos workflows distintos. Esta visión es crucial para evaluar el impacto potencial de modificaciones al proceso `FASTQC`.

### 4.4. Panel de Esquema

El panel de Esquema, ubicado en la barra lateral del Explorador (haga clic en ![ícono del Explorador](img/files_icon.png)), proporciona una visión general conveniente de todos los símbolos en su archivo actual. Esta característica le permite navegar y gestionar rápidamente la estructura de su código mostrando funciones, variables y otros elementos clave en una vista jerárquica.

![Panel de esquema](img/outline.png)

Use el panel de Esquema para navegar rápidamente a diferentes partes de su código sin usar el navegador de archivos.

### 4.5. Visualización DAG

La extensión de Nextflow de VS Code puede visualizar su workflow como un Grafo Acíclico Dirigido (DAG). Esto le ayuda a comprender el flujo de datos y las dependencias entre procesos. Abra `complex_workflow.nf` y haga clic en el botón "Preview DAG" arriba de `workflow {` (el segundo bloque `workflow` en este archivo):

![Vista previa DAG](img/dag_preview.png)

Este es solo el workflow de 'entrada', pero también puede previsualizar el DAG para los workflows internos haciendo clic en el botón "Preview DAG" arriba del workflow `RNASEQ_PIPELINE {` más arriba:

![Vista previa DAG workflow interno](img/dag_preview_inner.png)

Para este workflow, puede usar los nodos en el DAG para navegar a las definiciones de proceso correspondientes en el código. Haga clic en un nodo y lo llevará a la definición de proceso relevante en el editor. Particularmente cuando un workflow crece a un tamaño grande, esto puede realmente ayudarlo a navegar por el código y comprender cómo los procesos están conectados.

### Conclusión

Puede navegar workflows complejos eficientemente usando ir a definición, búsqueda de símbolos, encontrar referencias y visualización DAG para comprender la estructura del código y las dependencias.

### ¿Qué sigue?

Aprenda cómo trabajar efectivamente a través de múltiples archivos interconectados en proyectos Nextflow más grandes.

## 5. Trabajo a Través de Múltiples Archivos

El desarrollo real de Nextflow implica trabajar con múltiples archivos interconectados. Exploremos cómo VS Code le ayuda a gestionar proyectos complejos eficientemente.

### 5.1. Navegación Rápida de Archivos

Con `complex_workflow.nf` abierto, notará que importa varios módulos. Practiquemos la navegación rápida entre ellos.

Presione **Ctrl+P** (o **Cmd+P**) y comience a escribir "fast":

VS Code le mostrará archivos coincidentes. Seleccione `modules/fastqc.nf` para saltar allí instantáneamente. Esto es mucho más rápido que hacer clic a través del explorador de archivos cuando sabe aproximadamente qué archivo está buscando.

Intente esto con otros patrones:

- Escriba "star" para encontrar el archivo del módulo de alineamiento STAR (`star.nf`)
- Escriba "utils" para encontrar el archivo de funciones de utilidad (`utils.nf`)
- Escriba "config" para saltar a archivos de configuración (`nextflow.config`)

### 5.2. Editor Dividido para Desarrollo Multi-archivo

Cuando trabaja con módulos, a menudo necesita ver tanto el workflow principal como las definiciones de módulos simultáneamente. Configuremos esto:

1. Abra `complex_workflow.nf`
2. Abra `modules/fastqc.nf` en una nueva pestaña
3. Haga clic derecho en la pestaña `modules/fastqc.nf` y seleccione "Split Right"
4. Ahora puede ver ambos archivos lado a lado

![Editor dividido](img/split_editor.png)

Esto es invaluable cuando:

- Verifica interfaces de módulos mientras escribe llamadas de workflow, y la vista previa no es suficiente
- Compara procesos similares a través de diferentes módulos
- Depura flujo de datos entre workflow y módulos

### 5.3. Búsqueda en Todo el Proyecto

A veces necesita encontrar dónde se usan patrones específicos en todo su proyecto. Presione `Ctrl/Cmd+Shift+F` para abrir el panel de búsqueda.

Intente buscar `publishDir` en todo el espacio de trabajo:

![Búsqueda de proyecto](img/project_search.png)

Esto le muestra cada archivo que usa directorios de publicación, ayudándolo a:

- Comprender patrones de organización de salida
- Encontrar ejemplos de directivas específicas
- Asegurar consistencia a través de módulos

### Conclusión

Puede gestionar proyectos complejos multi-archivo usando navegación rápida de archivos, editores divididos y búsqueda en todo el proyecto para trabajar eficientemente a través de workflows y módulos.

### ¿Qué sigue?

Aprenda cómo las características de formato de código y mantenimiento mantienen sus workflows organizados y legibles.

---

## 6. Formato de Código y Mantenimiento

El formato apropiado del código es esencial no solo para la estética sino también para mejorar la legibilidad, comprensión y la facilidad de actualizar workflows complejos.

### 6.1. Formato Automático en Acción

Abra `basic_workflow.nf` y desordene deliberadamente el formato:

- Elimine algo de indentación: Resalte todo el documento y presione `shift+tab` muchas veces para eliminar tantas indentaciones como sea posible.
- Agregue espacios extra en lugares aleatorios: en la declaración `channel.fromPath`, agregue 30 espacios después del `(`.
- Rompa algunas líneas de manera incómoda: Agregue una nueva línea entre el operador `.view {` y la cadena `Processing sample:` pero no agregue una nueva línea correspondiente antes del paréntesis de cierre `}`.

Ahora presione `Shift+Alt+F` (o `Shift+Option+F` en MacOS) para auto-formatear:

VS Code inmediatamente:

- Corrige la indentación para mostrar la estructura del proceso claramente
- Alinea elementos similares consistentemente
- Elimina espacios en blanco innecesarios
- Mantiene saltos de línea legibles

Note que el formato automático puede no resolver todos los problemas de estilo de código. El servidor de lenguaje de Nextflow intenta mantener su código ordenado, pero también respeta sus preferencias personales en ciertas áreas. Por ejemplo, si elimina indentación dentro del bloque `script` de un proceso, el formateador lo dejará como está, ya que podría preferir intencionalmente ese estilo.

Actualmente, no hay una aplicación de estilo estricta para Nextflow, así que el servidor de lenguaje ofrece cierta flexibilidad. Sin embargo, aplicará consistentemente reglas de formato alrededor de definiciones de métodos y funciones para mantener claridad.

### 6.2. Características de Organización de Código

#### Comentar Rápidamente

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
- Documentar secciones de workflow

Use **Ctrl+/** (o **Cmd+/**) nuevamente para descomentar el código.

#### Plegado de Código para Visión General

En `complex_workflow.nf`, note las pequeñas flechas junto a las definiciones de proceso. Haga clic en ellas para plegar (colapsar) procesos:

![Plegado de código](img/code_folding.png)

Esto le da una visión general de alto nivel de su estructura de workflow sin perderse en detalles de implementación.

#### Coincidencia de Corchetes

Coloque su cursor junto a cualquier corchete `{` o `}` y VS Code resalta el corchete coincidente. Use **Ctrl+Shift+\\** (o **Cmd+Shift+\\**) para saltar entre corchetes coincidentes.

Esto es crucial para:

- Comprender límites de proceso
- Encontrar corchetes faltantes o extra
- Navegar estructuras de workflow anidadas

#### Selección y Edición Multi-línea

Para editar múltiples líneas simultáneamente, VS Code ofrece capacidades poderosas de multi-cursor:

- **Selección multi-línea**: Mantenga **Ctrl+Alt** (o **Cmd+Option** para MacOS) y use las teclas de flecha para seleccionar múltiples líneas
- **Indentación multi-línea**: Seleccione múltiples líneas y use **Tab** para indentar o **Shift+Tab** para des-indentar bloques completos

Esto es particularmente útil para:

- Indentar bloques de proceso completos consistentemente
- Agregar comentarios a múltiples líneas a la vez
- Editar definiciones de parámetros similares a través de múltiples procesos

### Conclusión

Puede mantener código limpio y legible usando formato automático, características de comentarios, plegado de código, coincidencia de corchetes y edición multi-línea para organizar workflows complejos eficientemente.

### ¿Qué sigue?

Aprenda cómo VS Code se integra con su flujo de trabajo de desarrollo más amplio más allá de solo editar código.

---

## 7. Integración del Flujo de Trabajo de Desarrollo

VS Code se integra bien con su flujo de trabajo de desarrollo más allá de solo editar código.

### 7.1. Integración de Control de Versiones

!!! note "Codespaces e Integración Git"

    Si está trabajando en **GitHub Codespaces**, algunas características de integración Git pueden no funcionar como se espera, particularmente atajos de teclado para Control de Fuente. También pudo haber declinado abrir el directorio como un repositorio Git durante la configuración inicial, lo cual está bien para propósitos de entrenamiento.

Si su proyecto es un repositorio git (como este lo es), VS Code muestra:

- Archivos modificados con indicadores de color
- Estado de Git en la barra de estado
- Vistas de diferencias en línea
- Capacidades de commit y push

Abra el panel de Control de Fuente usando el botón de control de fuente (![ícono de Control de fuente](img/source_control_icon.png)) (`Ctrl+Shift+G` o `Cmd+Shift+G` si está trabajando con VS Code localmente) para ver cambios de git y realizar commits directamente en el editor.

![Panel de Control de Fuente](img/source_control.png)

### 7.2. Ejecución e Inspección de Workflows

Ejecutemos un workflow y luego inspeccionemos los resultados. En la terminal integrada (`Ctrl+Shift+` backtick tanto en Windows como MacOS), ejecute el workflow básico:

```bash title="Ejecute el workflow básico"
nextflow run basic_workflow.nf --input data/sample_data.csv --output_dir results
```

Mientras el workflow se ejecuta, verá salida en tiempo real en la terminal. Después de completarse, puede usar VS Code para inspeccionar resultados sin salir de su editor:

1. **Navegue a directorios de trabajo**: Use el explorador de archivos o terminal para navegar `.nextflow/work`
2. **Abra archivos de registro**: Haga clic en rutas de archivos de registro en la salida de la terminal para abrirlos directamente en VS Code
3. **Inspeccione salidas**: Navegue por directorios de resultados publicados en el explorador de archivos
4. **Vea reportes de ejecución**: Abra reportes HTML directamente en VS Code o su navegador

Esto mantiene todo en un solo lugar en lugar de cambiar entre múltiples aplicaciones.

### Conclusión

Puede integrar VS Code con control de versiones y ejecución de workflow para gestionar todo su proceso de desarrollo desde una sola interfaz.

### ¿Qué sigue?

Vea cómo todas estas características del IDE funcionan juntas en su flujo de trabajo de desarrollo diario.

---

## 8. Resumen y Notas Rápidas

Aquí hay algunas notas rápidas sobre cada una de las características del IDE discutidas arriba:

### 8.1. Iniciando una Nueva Característica

1. **Apertura rápida de archivos** (`Ctrl+P` o `Cmd+P`) para encontrar módulos existentes relevantes
2. **Editor dividido** para ver procesos similares lado a lado
3. **Navegación de símbolos** (`Ctrl+Shift+O` o `Cmd+Shift+O`) para comprender estructura de archivo
4. **Auto-completado** para escribir nuevo código rápidamente

### 8.2. Depuración de Problemas

1. **Panel de problemas** (`Ctrl+Shift+M` o `Cmd+Shift+M`) para ver todos los errores a la vez
2. **Ir a definición** (`Ctrl-clic` o `Cmd-clic`) para comprender interfaces de proceso
3. **Encontrar todas las referencias** para ver cómo se usan los procesos
4. **Búsqueda en todo el proyecto** para encontrar patrones o problemas similares

### 8.3. Refactorización y Mejora

1. **Búsqueda en todo el proyecto** (`Ctrl+Shift+F` o `Cmd+Shift+F`) para encontrar patrones
2. **Auto-formato** (`Shift+Alt+F` o `Shift+Option+F`) para mantener consistencia
3. **Plegado de código** para enfocarse en estructura
4. **Integración Git** para rastrear cambios

---

## Resumen

Ahora ha tenido un recorrido rápido de las características del IDE de VS Code para desarrollo de Nextflow. Estas herramientas lo harán significativamente más productivo al:

- **Reducir errores** a través de verificación de sintaxis en tiempo real
- **Acelerar el desarrollo** con auto-completado inteligente
- **Mejorar la navegación** en workflows complejos multi-archivo
- **Mantener calidad** a través de formato consistente
- **Mejorar la comprensión** a través de resaltado avanzado y visualización de estructura

No esperamos que recuerde todo, pero ahora que sabe que estas características existen podrá encontrarlas cuando las necesite. A medida que continúe desarrollando workflows de Nextflow, estas características del IDE se volverán una segunda naturaleza, permitiéndole enfocarse en escribir código de alta calidad en lugar de luchar con sintaxis y estructura.

### ¿Qué sigue?

Aplique estas habilidades del IDE mientras trabaja en otros módulos de entrenamiento, por ejemplo:

- **[nf-test](nf-test.md)**: Cree suites de pruebas completas para sus workflows
- **[Hello nf-core](../../hello_nf-core/)**: Construya pipelines de calidad de producción con estándares de la comunidad

El verdadero poder de estas características del IDE emerge a medida que trabaja en proyectos más grandes y complejos. Comience a incorporarlas en su flujo de trabajo gradualmente—dentro de unas pocas sesiones, se volverán una segunda naturaleza y transformarán cómo aborda el desarrollo de Nextflow.

Desde detectar errores antes de que lo retrasen hasta navegar bases de código complejas con facilidad, estas herramientas lo harán un desarrollador más confiado y eficiente.

¡Feliz codificación!
