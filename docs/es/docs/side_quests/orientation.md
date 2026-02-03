# Orientación

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Traducción asistida por IA - [más información y sugerencias](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

El entorno de GitHub Codespaces contiene todo el software, código y datos necesarios para completar este curso de entrenamiento, por lo que no necesita instalar nada usted mismo.
Sin embargo, necesita una cuenta (gratuita) para iniciar sesión, y debería tomarse unos minutos para familiarizarse con la interfaz.

Si aún no lo ha hecho, por favor siga [este enlace](../../envsetup/) antes de continuar.

## Materiales proporcionados

A lo largo de este curso de entrenamiento, trabajaremos en el directorio `side-quests/`.
Este directorio contiene todos los archivos de código, datos de prueba y archivos accesorios que necesitará.

Siéntase libre de explorar los contenidos de este directorio; la forma más fácil de hacerlo es usar el explorador de archivos en el lado izquierdo del espacio de trabajo de GitHub Codespaces.
Alternativamente, puede usar el comando `tree`.
A lo largo del curso, usamos la salida de `tree` para representar la estructura y contenidos del directorio de forma legible, a veces con modificaciones menores para mayor claridad.

Aquí generamos una tabla de contenidos hasta el segundo nivel:

```bash
tree . -L 2
```

Si ejecuta esto dentro de `side-quests`, debería ver la siguiente salida:

```console title="Contenidos del directorio"
.
├── metadata
├── nf-core
├── nf-test
├── solutions
├── splitting_and_grouping
└── workflows_of_workflows
```

**Esto es un resumen de lo que debería saber para comenzar:**

- **Cada directorio corresponde a una misión secundaria individual.**
  Sus contenidos se detallan en la página de la misión secundaria correspondiente.

- **El directorio `solutions`** contiene los scripts de workflow y/o módulo completados que resultan de ejecutar varios pasos de cada misión secundaria.
  Están destinados a ser usados como referencia para verificar su trabajo y solucionar cualquier problema.

!!!tip "Consejo"

    Si por cualquier razón sale de este directorio, siempre puede ejecutar este comando para regresar a él:

    ```bash
    cd /workspaces/training/side-quests
    ```

Ahora, para comenzar el curso, haga clic en la flecha en la esquina inferior derecha de esta página.
