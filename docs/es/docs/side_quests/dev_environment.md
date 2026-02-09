# Entorno de Desarrollo

Los Entornos de Desarrollo Integrados (IDEs) modernos pueden transformar dramáticamente tu experiencia de desarrollo con Nextflow. Esta misión secundaria se enfoca específicamente en aprovechar VS Code y su extensión de Nextflow para escribir código más rápido, detectar errores temprano y navegar flujos de trabajo complejos de manera eficiente.

!!! note "Este no es un tutorial tradicional"

    A diferencia de otros módulos de capacitación, esta guía está organizada como una colección de sugerencias rápidas, consejos y ejemplos prácticos en lugar de un tutorial paso a paso. Cada sección puede explorarse de forma independiente según tus intereses y necesidades actuales de desarrollo. Siéntete libre de saltar entre secciones y enfocarte en las características que serán más útiles inmediatamente para el desarrollo de tu workflow.

## Qué deberías saber primero

Esta guía asume que has completado el curso de capacitación [Hello Nextflow](../hello_nextflow/) y te sientes cómodo con los conceptos fundamentales de Nextflow, incluyendo:

- **Estructura básica de workflow**: Comprender procesos, workflows y cómo se conectan entre sí
- **Operaciones de canal**: Crear canales, pasar datos entre procesos y usar operadores básicos
- **Módulos y organización**: Crear módulos reutilizables y usar declaraciones include
- **Fundamentos de configuración**: Usar `nextflow.config` para parámetros, directivas de proceso y perfiles

## Qué aprenderás aquí

Esta guía se enfoca en **características de productividad del IDE** que te convertirán en un desarrollador de Nextflow más eficiente:

- **Resaltado de sintaxis avanzado**: Entender qué te muestra VS Code sobre la estructura de tu código
- **Autocompletado inteligente**: Aprovechar sugerencias contextuales para escribir código más rápido
- **Detección de errores y diagnósticos**: Detectar errores de sintaxis antes de ejecutar tu workflow
- **Navegación de código**: Moverte rápidamente entre procesos, módulos y definiciones
- **Formato y organización**: Mantener un estilo de código consistente y legible
- **Desarrollo asistido por IA** (opcional): Usar herramientas de IA modernas integradas con tu IDE

!!! info "¿Por qué características del IDE ahora?"

    Probablemente ya has estado usando VS Code durante el curso [Hello Nextflow](../hello_nextflow/), pero mantuvimos el enfoque en aprender los fundamentos de Nextflow en lugar de las características del IDE. Ahora que te sientes cómodo con conceptos básicos de Nextflow como procesos, workflows, canales y módulos, estás listo para aprovechar las sofisticadas características del IDE que te convertirán en un desarrollador más eficiente.

    Piensa en esto como "subir de nivel" tu entorno de desarrollo - el mismo editor que has estado usando tiene capacidades mucho más poderosas que se vuelven verdaderamente valiosas una vez que entiendes en qué te están ayudando.

---

## 0. Configuración y Calentamiento

Configuremos un espacio de trabajo específicamente para explorar las características del IDE:

```bash title="Navigate to the IDE features directory"
cd side-quests/ide_features
```

Abre este directorio en VS Code:

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

!!! note "Acerca de los Archivos de Ejemplo"

    - `basic_workflow.nf` es un workflow básico funcional que puedes ejecutar y modificar
    - `complex_workflow.nf` está diseñado solo para ilustración para demostrar características de navegación - puede que no se ejecute exitosamente pero muestra una estructura de workflow multi-archivo realista

### Atajos de Teclado

Algunas de las características en esta guía usarán atajos de teclado opcionales. Es posible que estés accediendo a este material a través de GitHub Codespaces en el navegador, y en este caso a veces los atajos no funcionarán como se espera porque se usan para otras cosas en tu sistema.

Si estás ejecutando VS Code localmente, como probablemente lo harás cuando realmente estés escribiendo workflows, los atajos funcionarán como se describe.

Si estás usando una Mac, algunos (no todos) atajos de teclado usarán "cmd" en lugar de "ctrl", y lo indicaremos en el texto como `Ctrl/Cmd`.

### 0.1. Instalación de la Extensión de Nextflow

!!! note "¿Ya Usas Devcontainers?"

    Si estás trabajando en **GitHub Codespaces** o usando un **devcontainer local**, la extensión de Nextflow probablemente ya está instalada y configurada para ti. Puedes omitir los pasos de instalación manual a continuación y proceder directamente a explorar las características de la extensión.

Para instalar la extensión manualmente:

1. Abre VS Code
2. Ve a la vista de Extensiones haciendo clic en el ícono de extensiones a la izquierda: ![ícono de extensiones](img/extensions_icon.png) (atajo `Ctrl/Cmd+Shift+X` si estás ejecutando VSCode localmente)
3. Busca "Nextflow"
4. Instala la extensión oficial de Nextflow

![Instalar Extensión de Nextflow](img/install_extension.png)

### 0.2. Diseño del Espacio de Trabajo

Como has estado usando VS Code durante Hello Nextflow, ya estás familiarizado con los conceptos básicos. Aquí te mostramos cómo organizar tu espacio de trabajo eficientemente para esta sesión:

- **Área del Editor**: Para ver y editar archivos. Puedes dividir esto en múltiples paneles para comparar archivos lado a lado.
- **Explorador de Archivos** clic (![ícono del explorador de archivos](img/files_icon.png)) (`Ctrl/Cmd+Shift+E`): Los archivos y carpetas locales en tu sistema. Mantén esto abierto a la izquierda para navegar entre archivos
- **Terminal Integrada** (`Ctrl+Shift+` acento grave tanto para Windows como MacOS): Una terminal para interactuar con la computadora en la parte inferior. Usa esto para ejecutar Nextflow u otros comandos.
- **Panel de Problemas** (`Ctrl+Shift+M`): VS Code mostrará aquí cualquier error y problema que detecte. Esto es útil para resaltar problemas de un vistazo.

Puedes arrastrar paneles o ocultarlos (`Ctrl/Cmd+B` para alternar la barra lateral) para personalizar tu diseño mientras trabajamos con los ejemplos.

### Conclusión

Tienes VS Code configurado con la extensión de Nextflow y entiendes el diseño del espacio de trabajo para un desarrollo eficiente.

### ¿Qué sigue?

Aprende cómo el resaltado de sintaxis te ayuda a entender la estructura del código Nextflow de un vistazo.

---

## 1. Resaltado de Sintaxis y Estructura del Código

Ahora que tu espacio de trabajo está configurado, exploremos cómo el resaltado de sintaxis de VS Code te ayuda a leer y escribir código Nextflow de manera más efectiva.

### 1.1. Elementos de Sintaxis de Nextflow

Abre `basic_workflow.nf` para ver el resaltado de sintaxis en acción:

![Demostración de Sintaxis](img/syntax_showcase.png)

Observa cómo VS Code resalta:

- **Palabras clave** (`process`, `workflow`, `input`, `output`, `script`) en colores distintos
- **Literales de cadena** y **parámetros** con diferentes estilos
- **Comentarios** en un color atenuado
- **Variables** y **llamadas a funciones** con énfasis apropiado
- **Bloques de código** con guías de indentación adecuadas

!!! note "Colores Dependientes del Tema"

    Los colores específicos que ves dependerán de tu tema de VS Code (modo oscuro/claro), configuración de colores y cualquier personalización que hayas hecho. Lo importante es que diferentes elementos de sintaxis se distingan visualmente entre sí, haciendo que la estructura del código sea más fácil de entender independientemente del esquema de colores que hayas elegido.

### 1.2. Entendiendo la Estructura del Código

El resaltado de sintaxis te ayuda a identificar rápidamente:

- **Límites de proceso**: Distinción clara entre diferentes procesos
- **Bloques de entrada/salida**: Fácil de detectar definiciones de flujo de datos
- **Bloques de script**: Los comandos reales que se están ejecutando
- **Operaciones de canal**: Pasos de transformación de datos
- **Directivas de configuración**: Configuraciones específicas del proceso

Esta organización visual se vuelve invaluable cuando trabajas con workflows complejos que contienen múltiples procesos y flujos de datos intrincados.

### Conclusión

Entiendes cómo el resaltado de sintaxis de VS Code te ayuda a leer la estructura del código Nextflow e identificar diferentes elementos del lenguaje para un desarrollo más rápido.

### ¿Qué sigue?

Aprende cómo el autocompletado inteligente acelera la escritura de código con sugerencias contextuales.

---

## 2. Autocompletado Inteligente

Las características de autocompletado de VS Code te ayudan a escribir código más rápido y con menos errores al sugerir opciones apropiadas según el contexto.

### 2.1. Sugerencias Contextuales

Las opciones de autocompletado varían dependiendo de dónde estés en tu código:

#### Operaciones de Canal

Abre `basic_workflow.nf` nuevamente e intenta escribir `channel.` en el bloque workflow:

![Autocompletado de canal](img/autocomplete_channel.png)

Verás sugerencias para:

- `fromPath()` - Crear canal desde rutas de archivo
- `fromFilePairs()` - Crear canal desde archivos emparejados
- `of()` - Crear canal desde valores
- `fromSRA()` - Crear canal desde accesiones SRA
- Y muchos más...

Esto te ayuda a encontrar rápidamente la factory de canal correcta para usar sin necesidad de recordar nombres exactos de métodos.

También puedes descubrir los operadores disponibles para aplicar a canales. Por ejemplo, escribe `FASTQC.out.html.` para ver las operaciones disponibles:

![Autocompletado de operaciones de canal](img/autocomplete_operators.png)

#### Directivas de Proceso

Dentro de un bloque de script de proceso, escribe `task.` para ver las propiedades de tiempo de ejecución disponibles:

![Autocompletado de propiedades de tarea](img/autocomplete_task.png)

#### Configuración

Abre nextflow.config y escribe `process.` en cualquier lugar para ver las directivas de proceso disponibles:

![Autocompletado de configuración](img/autocomplete_config.png)

Verás sugerencias para:

- `executor`
- `memory`
- `cpus`

Esto ahorra tiempo al configurar procesos y funciona en diferentes ámbitos de configuración. Por ejemplo, intenta escribir `docker.` para ver opciones de configuración específicas de Docker.

### Conclusión

Puedes usar el autocompletado inteligente de VS Code para descubrir operaciones de canal disponibles, directivas de proceso y opciones de configuración sin memorizar la sintaxis.

### ¿Qué sigue?

Aprende cómo la detección de errores en tiempo real te ayuda a detectar problemas antes de ejecutar tu workflow, simplemente leyendo el código.

## 3. Detección de Errores y Diagnósticos

La detección de errores en tiempo real de VS Code te ayuda a detectar problemas antes de ejecutar tu workflow.

### 3.1. Detección de Errores de Sintaxis

Creemos un error deliberado para ver la detección en acción. Abre `basic_workflow.nf` y cambia el nombre del proceso de `FASTQC` a `FASTQ` (o cualquier otro nombre inválido). VS Code inmediatamente resaltará el error en el bloque workflow con un subrayado ondulado rojo:

![Subrayado de error](img/error_underline.png)

### 3.2. Panel de Problemas

Más allá del resaltado de errores individuales, VS Code proporciona un Panel de Problemas centralizado que agrega todos los errores, advertencias y mensajes de información en tu espacio de trabajo. Ábrelo con `Ctrl/Cmd+Shift+M` y usa el ícono de filtro para mostrar solo los errores relevantes al archivo actual:

![Filtrar el panel de problemas](img/active_file.png)

Haz clic en cualquier problema para saltar directamente a la línea problemática

![Panel de Problemas](img/problems_panel.png)

Corrige el error cambiando el nombre del proceso de vuelta a `FASTQC`.

### 3.3. Patrones de Error Comunes

Los errores comunes en la sintaxis de Nextflow incluyen:

- **Corchetes faltantes**: `{` o `}` sin coincidencia
- **Bloques incompletos**: Secciones requeridas faltantes en procesos
- **Sintaxis inválida**: DSL de Nextflow mal formado
- **Errores tipográficos en palabras clave**: Directivas de proceso mal escritas
- **Desajustes de canal**: Incompatibilidades de tipo

El servidor de lenguaje de Nextflow resalta estos problemas en el Panel de Problemas. Puedes revisarlos temprano para evitar errores de sintaxis al ejecutar un pipeline.

### Conclusión

Puedes usar la detección de errores de VS Code y el Panel de Problemas para detectar errores de sintaxis y problemas antes de ejecutar tu workflow, ahorrando tiempo y evitando frustraciones.

### ¿Qué sigue?

Aprende cómo navegar eficientemente entre procesos, módulos y definiciones en workflows complejos.

---

## 4. Navegación de Código y Gestión de Símbolos

La navegación eficiente es crucial cuando trabajas con workflows complejos que abarcan múltiples archivos. Para entender esto, reemplaza la definición del proceso en `basic_workflow.nf` con una importación del módulo que te hemos proporcionado:

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

Si pasas el mouse sobre un nombre de proceso como `FASTQC`, verás una ventana emergente con la interfaz del módulo (entradas y salidas):

![Ir a definición](img/syntax.png)

Esta característica es particularmente valiosa al crear workflows, ya que te permite entender la interfaz del módulo sin abrir el archivo del módulo directamente.

Puedes navegar rápidamente a cualquier definición de proceso, módulo o variable usando **Ctrl/Cmd-clic**. Pasa el mouse sobre el enlace al archivo del módulo en la parte superior del script y sigue el enlace como se sugiere:

![Seguir enlace](img/follow_link.png)

Lo mismo funciona para nombres de proceso. Regresa a `basic_workflow.nf` e intenta esto en el nombre del proceso `FASTQC` en el bloque workflow. Esto te enlaza directamente al nombre del proceso (que es el mismo que el archivo del módulo en este ejemplo, pero podría estar a mitad de camino en un archivo mucho más grande).

Para regresar a donde estabas, usa **Alt+←** (o **Ctrl+-** en Mac). Esta es una forma poderosa de explorar código sin perder tu lugar.

Ahora exploremos la navegación en un workflow más complejo usando `complex_workflow.nf` (el archivo solo para ilustración mencionado anteriormente). Este workflow contiene múltiples procesos definidos en archivos de módulo separados, así como algunos en línea. Aunque las estructuras complejas de múltiples archivos pueden ser desafiantes de navegar manualmente, la capacidad de saltar a definiciones hace que la exploración sea mucho más manejable.

1. Abre `complex_workflow.nf`
2. Navega a definiciones de módulo
3. Usa **Alt+←** (o **Ctrl+-**) para navegar de regreso
4. Navega al nombre del proceso `FASTQC` en el bloque workflow. Esto te enlaza directamente al nombre del proceso (que es el mismo que el archivo del módulo en este ejemplo, pero podría estar a mitad de camino en un archivo mucho más grande).
5. Navega de regreso nuevamente
6. Navega al proceso `TRIM_GALORE` en el bloque workflow. Este está definido en línea, por lo que no te llevará a un archivo separado, pero aún te mostrará la definición del proceso, y aún puedes navegar de regreso a donde estabas.

### 4.2. Navegación de Símbolos

Con `complex_workflow.nf` aún abierto, puedes obtener una vista general de todos los símbolos en el archivo escribiendo `@` en la barra de búsqueda en la parte superior de VSCode (el atajo de teclado es `Ctrl/Cmd+Shift+O`, pero puede que no funcione en Codespaces). Esto abre el panel de navegación de símbolos, que lista todos los símbolos en el archivo actual:

![Navegación de símbolos](img/symbols.png)

Esto muestra:

- Todas las definiciones de proceso
- Definiciones de workflow (hay dos workflows definidos en este archivo)
- Definiciones de función

Comienza a escribir para filtrar resultados.

### 4.3. Encontrar Todas las Referencias

Entender dónde se usa un proceso o variable en todo tu código base puede ser muy útil. Por ejemplo, si quieres encontrar todas las referencias al proceso `FASTQC`, comienza navegando a su definición. Puedes hacer esto abriendo `modules/fastqc.nf` directamente, o usando la característica de navegación rápida de VS Code con `Ctrl/Cmd-clic` como hicimos arriba. Una vez en la definición del proceso, haz clic derecho en el nombre del proceso `FASTQC` y selecciona "Find All References" del menú contextual para ver todas las instancias donde se usa.

![Encontrar referencias](img/references.png)

Esta característica muestra todas las instancias donde se hace referencia a `FASTQC` dentro de tu espacio de trabajo, incluyendo su uso en los dos workflows distintos. Esta información es crucial para evaluar el impacto potencial de modificaciones al proceso `FASTQC`.

### 4.4. Panel de Esquema

El panel de Esquema, ubicado en la barra lateral del Explorador (haz clic en ![Ícono del Explorador](img/files_icon.png)), proporciona una vista general conveniente de todos los símbolos en tu archivo actual. Esta característica te permite navegar rápidamente y gestionar la estructura de tu código mostrando funciones, variables y otros elementos clave en una vista jerárquica.

![Panel de esquema](img/outline.png)

Usa el panel de Esquema para navegar rápidamente a diferentes partes de tu código sin usar el navegador de archivos.

### 4.5. Visualización DAG

La extensión de Nextflow de VS Code puede visualizar tu workflow como un Grafo Acíclico Dirigido (DAG). Esto te ayuda a entender el flujo de datos y las dependencias entre procesos. Abre `complex_workflow.nf` y haz clic en el botón "Preview DAG" arriba de `workflow {` (el segundo bloque `workflow` en este archivo):

![Vista previa DAG](img/dag_preview.png)

Este es solo el workflow de 'entrada', pero también puedes previsualizar el DAG para los workflows internos haciendo clic en el botón "Preview DAG" arriba del workflow `RNASEQ_PIPELINE {` más arriba:

![Vista previa DAG workflow interno](img/dag_preview_inner.png)

Para este workflow, puedes usar los nodos en el DAG para navegar a las definiciones de proceso correspondientes en el código. Haz clic en un nodo y te llevará a la definición de proceso relevante en el editor. Particularmente cuando un workflow crece a un tamaño grande, esto realmente puede ayudarte a navegar por el código y entender cómo están conectados los procesos.

### Conclusión

Puedes navegar workflows complejos eficientemente usando ir-a-definición, búsqueda de símbolos, encontrar referencias y visualización DAG para entender la estructura del código y las dependencias.

### ¿Qué sigue?

Aprende cómo trabajar efectivamente con múltiples archivos interconectados en proyectos Nextflow más grandes.

## 5. Trabajando con Múltiples Archivos

El desarrollo real de Nextflow implica trabajar con múltiples archivos interconectados. Exploremos cómo VS Code te ayuda a gestionar proyectos complejos eficientemente.

### 5.1. Navegación Rápida de Archivos

Con `complex_workflow.nf` abierto, notarás que importa varios módulos. Practiquemos la navegación rápida entre ellos.

Presiona **Ctrl+P** (o **Cmd+P**) y comienza a escribir "fast":

VS Code te mostrará archivos coincidentes. Selecciona `modules/fastqc.nf` para saltar allí instantáneamente. Esto es mucho más rápido que hacer clic a través del explorador de archivos cuando sabes aproximadamente qué archivo estás buscando.

Prueba esto con otros patrones:

- Escribe "star" para encontrar el archivo del módulo de alineamiento STAR (`star.nf`)
- Escribe "utils" para encontrar el archivo de funciones de utilidad (`utils.nf`)
- Escribe "config" para saltar a archivos de configuración (`nextflow.config`)

### 5.2. Editor Dividido para Desarrollo Multi-archivo

Cuando trabajas con módulos, a menudo necesitas ver tanto el workflow principal como las definiciones de módulo simultáneamente. Configuremos esto:

1. Abre `complex_workflow.nf`
2. Abre `modules/fastqc.nf` en una nueva pestaña
3. Haz clic derecho en la pestaña `modules/fastqc.nf` y selecciona "Split Right"
4. Ahora puedes ver ambos archivos lado a lado

![Editor dividido](img/split_editor.png)

Esto es invaluable cuando:

- Verificas interfaces de módulo mientras escribes llamadas de workflow, y la vista previa no es suficiente
- Comparas procesos similares en diferentes módulos
- Depuras flujo de datos entre workflow y módulos

### 5.3. Búsqueda en Todo el Proyecto

A veces necesitas encontrar dónde se usan patrones específicos en todo tu proyecto. Presiona `Ctrl/Cmd+Shift+F` para abrir el panel de búsqueda.

Intenta buscar `publishDir` en todo el espacio de trabajo:

![Búsqueda en proyecto](img/project_search.png)

Esto te muestra cada archivo que usa directorios de publicación, ayudándote a:

- Entender patrones de organización de salida
- Encontrar ejemplos de directivas específicas
- Asegurar consistencia entre módulos

### Conclusión

Puedes gestionar proyectos complejos de múltiples archivos usando navegación rápida de archivos, editores divididos y búsqueda en todo el proyecto para trabajar eficientemente con workflows y módulos.

### ¿Qué sigue?

Aprende cómo las características de formato de código y mantenimiento mantienen tus workflows organizados y legibles.

---

## 6. Formato de Código y Mantenimiento

El formato adecuado del código es esencial no solo para la estética sino también para mejorar la legibilidad, comprensión y facilidad de actualización de workflows complejos.

### 6.1. Formato Automático en Acción

Abre `basic_workflow.nf` y desordena deliberadamente el formato:

- Elimina algo de indentación: Resalta todo el documento y presiona `shift+tab` muchas veces para eliminar tantas indentaciones como sea posible.
- Agrega espacios extra en lugares aleatorios: en la declaración `channel.fromPath`, agrega 30 espacios después del `(`.
- Rompe algunas líneas de manera incómoda: Agrega una nueva línea entre el operador `.view {` y la cadena `Processing sample:` pero no agregues una nueva línea correspondiente antes del paréntesis de cierre `}`.

Ahora presiona `Shift+Alt+F` (o `Shift+Option+F` en MacOS) para auto-formatear:

VS Code inmediatamente:

- Corrige la indentación para mostrar la estructura del proceso claramente
- Alinea elementos similares consistentemente
- Elimina espacios en blanco innecesarios
- Mantiene saltos de línea legibles

Ten en cuenta que el formato automático puede no resolver todos los problemas de estilo de código. El servidor de lenguaje de Nextflow busca mantener tu código ordenado, pero también respeta tus preferencias personales en ciertas áreas. Por ejemplo, si eliminas la indentación dentro del bloque `script` de un proceso, el formateador lo dejará como está, ya que podrías preferir intencionalmente ese estilo.

Actualmente, no hay una aplicación estricta de estilo para Nextflow, por lo que el servidor de lenguaje ofrece cierta flexibilidad. Sin embargo, aplicará consistentemente reglas de formato alrededor de definiciones de métodos y funciones para mantener la claridad.

### 6.2. Características de Organización de Código

#### Comentarios Rápidos

Selecciona un bloque de código en tu workflow y presiona **Ctrl+/** (o **Cmd+/**) para comentarlo:

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

Usa **Ctrl+/** (o **Cmd+/**) nuevamente para descomentar el código.

#### Plegado de Código para Vista General

En `complex_workflow.nf`, observa las pequeñas flechas junto a las definiciones de proceso. Haz clic en ellas para plegar (colapsar) procesos:

![Plegado de código](img/code_folding.png)

Esto te da una vista general de alto nivel de la estructura de tu workflow sin perderte en detalles de implementación.

#### Coincidencia de Corchetes

Coloca tu cursor junto a cualquier corchete `{` o `}` y VS Code resalta el corchete coincidente. Usa **Ctrl+Shift+\\** (o **Cmd+Shift+\\**) para saltar entre corchetes coincidentes.

Esto es crucial para:

- Entender límites de proceso
- Encontrar corchetes faltantes o extra
- Navegar estructuras de workflow anidadas

#### Selección y Edición Multi-línea

Para editar múltiples líneas simultáneamente, VS Code ofrece poderosas capacidades de multi-cursor:

- **Selección multi-línea**: Mantén **Ctrl+Alt** (o **Cmd+Option** para MacOS) y usa las teclas de flecha para seleccionar múltiples líneas
- **Indentación multi-línea**: Selecciona múltiples líneas y usa **Tab** para indentar o **Shift+Tab** para reducir indentación de bloques enteros

Esto es particularmente útil para:

- Indentar bloques de proceso enteros consistentemente
- Agregar comentarios a múltiples líneas a la vez
- Editar definiciones de parámetros similares en múltiples procesos

### Conclusión

Puedes mantener código limpio y legible usando formato automático, características de comentarios, plegado de código, coincidencia de corchetes y edición multi-línea para organizar workflows complejos eficientemente.

### ¿Qué sigue?

Aprende cómo VS Code se integra con tu flujo de trabajo de desarrollo más amplio más allá de solo editar código.

---

## 7. Integración del Flujo de Trabajo de Desarrollo

VS Code se integra bien con tu flujo de trabajo de desarrollo más allá de solo editar código.

### 7.1. Integración de Control de Versiones

!!! note "Codespaces e Integración Git"

    Si estás trabajando en **GitHub Codespaces**, algunas características de integración Git pueden no funcionar como se espera, particularmente atajos de teclado para Control de Código Fuente. También puedes haber declinado abrir el directorio como un repositorio Git durante la configuración inicial, lo cual está bien para propósitos de capacitación.

Si tu proyecto es un repositorio git (como este lo es), VS Code muestra:

- Archivos modificados con indicadores de color
- Estado de Git en la barra de estado
- Vistas de diferencias en línea
- Capacidades de commit y push

Abre el panel de Control de Código Fuente usando el botón de control de código fuente (![Ícono de control de código fuente](img/source_control_icon.png)) (`Ctrl+Shift+G` o `Cmd+Shift+G` si estás trabajando con VSCode localmente) para ver cambios de git y hacer commits directamente en el editor.

![Panel de Control de Código Fuente](img/source_control.png)

### 7.2. Ejecución e Inspección de Workflows

Ejecutemos un workflow y luego inspeccionemos los resultados. En la terminal integrada (`Ctrl+Shift+` acento grave tanto en Windows como MacOS), ejecuta el workflow básico:

```bash title="Run the basic workflow"
nextflow run basic_workflow.nf --input data/sample_data.csv --output_dir results
```

Mientras el workflow se ejecuta, verás salida en tiempo real en la terminal. Después de completarse, puedes usar VS Code para inspeccionar resultados sin salir de tu editor:

1. **Navegar a directorios de trabajo**: Usa el explorador de archivos o terminal para navegar `.nextflow/work`
2. **Abrir archivos de log**: Haz clic en rutas de archivos de log en la salida de terminal para abrirlos directamente en VS Code
3. **Inspeccionar salidas**: Navega directorios de resultados publicados en el explorador de archivos
4. **Ver reportes de ejecución**: Abre reportes HTML directamente en VS Code o tu navegador

Esto mantiene todo en un solo lugar en lugar de cambiar entre múltiples aplicaciones.

### Conclusión

Puedes integrar VS Code con control de versiones y ejecución de workflow para gestionar todo tu proceso de desarrollo desde una sola interfaz.

### ¿Qué sigue?

Ve cómo todas estas características del IDE trabajan juntas en tu flujo de trabajo de desarrollo diario.

---

## 8. Recapitulación y notas rápidas

Aquí hay algunas notas rápidas sobre cada una de las características del IDE discutidas arriba:

### 8.1. Iniciando una Nueva Característica

1. **Apertura rápida de archivo** (`Ctrl+P` o `Cmd+P`) para encontrar módulos existentes relevantes
2. **Editor dividido** para ver procesos similares lado a lado
3. **Navegación de símbolos** (`Ctrl+Shift+O` o `Cmd+Shift+O`) para entender la estructura del archivo
4. **Autocompletado** para escribir nuevo código rápidamente

### 8.2. Depuración de Problemas

1. **Panel de problemas** (`Ctrl+Shift+M` o `Cmd+Shift+M`) para ver todos los errores a la vez
2. **Ir a definición** (`Ctrl-clic` o `Cmd-clic`) para entender interfaces de proceso
3. **Encontrar todas las referencias** para ver cómo se usan los procesos
4. **Búsqueda en todo el proyecto** para encontrar patrones o problemas similares

### 8.3. Refactorización y Mejora

1. **Búsqueda en todo el proyecto** (`Ctrl+Shift+F` o `Cmd+Shift+F`) para encontrar patrones
2. **Auto-formato** (`Shift+Alt+F` o `Shift+Option+F`) para mantener consistencia
3. **Plegado de código** para enfocarse en la estructura
4. **Integración Git** para rastrear cambios

---

## Resumen

Ahora has tenido un recorrido rápido de las características del IDE de VS Code para desarrollo de Nextflow. Estas herramientas te harán significativamente más productivo al:

- **Reducir errores** a través de verificación de sintaxis en tiempo real
- **Acelerar el desarrollo** con autocompletado inteligente
- **Mejorar la navegación** en workflows complejos de múltiples archivos
- **Mantener la calidad** a través de formato consistente
- **Mejorar la comprensión** a través de resaltado avanzado y visualización de estructura

No esperamos que recuerdes todo, pero ahora que sabes que estas características existen podrás encontrarlas cuando las necesites. A medida que continúes desarrollando workflows de Nextflow, estas características del IDE se volverán una segunda naturaleza, permitiéndote enfocarte en escribir código de alta calidad en lugar de luchar con sintaxis y estructura.

### ¿Qué sigue?

Aplica estas habilidades del IDE mientras trabajas en otros módulos de capacitación, por ejemplo:

- **[nf-test](nf-test.md)**: Crea suites de prueba completas para tus workflows
- **[Hello nf-core](../../hello_nf-core/)**: Construye pipelines de calidad de producción con estándares de la comunidad

El verdadero poder de estas características del IDE emerge a medida que trabajas en proyectos más grandes y complejos. Comienza a incorporarlas en tu flujo de trabajo gradualmente—en unas pocas sesiones, se volverán una segunda naturaleza y transformarán cómo abordas el desarrollo de Nextflow.

Desde detectar errores antes de que te ralenticen hasta navegar bases de código complejas con facilidad, estas herramientas te convertirán en un desarrollador más confiado y eficiente.

¡Feliz codificación!
