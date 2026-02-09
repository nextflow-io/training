# Orientación

El entorno de GitHub Codespaces contiene todo el software, código y datos necesarios para trabajar en este curso de capacitación, por lo que no necesitas instalar nada por tu cuenta.
Sin embargo, sí necesitas una cuenta (gratuita) para iniciar sesión, y deberías tomarte unos minutos para familiarizarte con la interfaz.

Si aún no lo has hecho, por favor sigue [este enlace](../../envsetup/) antes de continuar.

## Materiales proporcionados

A lo largo de este curso de capacitación, trabajaremos en el directorio `side-quests/`.
Este directorio contiene todos los archivos de código, datos de prueba y archivos accesorios que necesitarás.

Siéntete libre de explorar el contenido de este directorio; la forma más fácil de hacerlo es usar el explorador de archivos en el lado izquierdo del espacio de trabajo de GitHub Codespaces.
Alternativamente, puedes usar el comando `tree`.
A lo largo del curso, usamos la salida de `tree` para representar la estructura y contenido del directorio de forma legible, a veces con modificaciones menores para mayor claridad.

Aquí generamos una tabla de contenidos hasta el segundo nivel:

```bash
tree . -L 2
```

Si ejecutas esto dentro de `side-quests`, deberías ver la siguiente salida:

```console title="Directory contents"
.
├── metadata
├── nf-core
├── nf-test
├── solutions
├── splitting_and_grouping
└── workflows_of_workflows
```

**Aquí hay un resumen de lo que debes saber para comenzar:**

- **Cada directorio corresponde a una side quest individual.**
  Sus contenidos se detallan en la página de la side quest correspondiente.

- **El directorio `solutions`** contiene los scripts de workflow y/o módulos completados que resultan de ejecutar varios pasos de cada side quest.
  Están destinados a ser usados como referencia para verificar tu trabajo y solucionar cualquier problema.

!!!tip

    Si por alguna razón sales de este directorio, siempre puedes ejecutar este comando para regresar a él:

    ```bash
    cd /workspaces/training/side-quests
    ```

Ahora, para comenzar el curso, haz clic en la flecha en la esquina inferior derecha de esta página.
