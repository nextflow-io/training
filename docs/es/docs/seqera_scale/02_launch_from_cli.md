# Parte 2: Ejecutar pipelines desde la línea de comandos

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Traducción asistida por IA - [más información y sugerencias](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

En la [Parte 1](./01_run_with_seqera.md), ejecutó nf-core/rnaseq desde la interfaz web de Seqera.
Ahora haremos lo mismo desde la línea de comandos usando el CLI `tw`, y agregaremos un nuevo pipeline a su workspace.

---

## 1. Ejecutar pipelines desde la línea de comandos

En la vista de ejecución, haga clic en la pestaña **Command line**.
Verá el comando exacto `nextflow run` que Platform construyó y envió en su nombre — el mismo tipo de comando que ha estado ejecutando manualmente en el curso Use nf-core.

Platform no reemplaza a Nextflow; lo orquesta.
Todo lo que puede hacer a través de la interfaz web, también puede hacerlo desde una terminal usando el CLI `tw`, la herramienta de línea de comandos para interactuar con la API de Platform.
Esto es útil para automatizar ejecuciones desde scripts o pipelines de CI/CD.

Haremos esto ahora desde el mismo codespace que usó para los cursos anteriores.

### 1.1. Instalar el CLI tw

Ejecute los siguientes comandos en la terminal de su Codespace para descargar e instalar el binario `tw`:

```bash
curl -fsSL https://github.com/seqeralabs/tower-cli/releases/latest/download/tw-linux-x86_64 -o tw
chmod +x tw
sudo mv tw /usr/local/bin/
```

Verifique la instalación:

```bash
tw --version
```

??? success "Salida del comando"

    ```console
    Tower CLI version 0.40.0 (build 26db579)
    ```

El CLI `tw` está instalado y listo para configurarse.

### 1.2. Obtener un token de acceso

El CLI `tw` se autentica con Seqera usando un token de acceso personal.

1. En la interfaz web de Seqera, haga clic en su avatar en la esquina superior derecha y seleccione **Your tokens**.
2. Haga clic en **Add token**, asígnele un nombre (por ejemplo, `training`) y haga clic en **Add**.
3. Copie el valor del token — solo se mostrará una vez.
   Si no lo guarda de inmediato, deberá generar otro.

### 1.3. Configurar el CLI

Por conveniencia, vamos a configurar un archivo de configuración que contenga el
token de acceso que acaba de generar y el identificador del workspace.

Abra el archivo `.seqera_config` en este directorio en el editor y establezca las dos variables:

- **`TOWER_ACCESS_TOKEN`**: el token que generó en la sección 1.2
- **`TOWER_WORKSPACE_ID`**: el ID numérico de su workspace (la columna `ID` en `tw workspaces list`, que ejecutará en la sección 1.4)

Una vez que los valores estén completos, cargue la configuración:

```bash
source .seqera_config
```

Verifique la conexión:

```bash
tw info
```

??? success "Salida del comando"

    ```console
        Details
    -------------------------+-----------------------------
     Tower API endpoint      | https://api.cloud.seqera.io
     Tower API version       | 1.150.0
     Tower version           | 26.1.0-cycle54
     CLI version             | 0.30.0 (fde9dec)
     CLI minimum API version | 1.148.0
     Authenticated user      | <your-name>

    System health status
    ---------------------------------------+----
     Remote API server connection check    | OK
     Tower API version check               | OK
     Authentication API credential's token | OK
    ```

El CLI `tw` ahora está autenticado y conectado a su cuenta de Seqera.
Ejecute `source .seqera_config` al inicio de cada sesión de Codespace para recargar la configuración.

!!! tip "Consejo"

    Si su workspace no tiene un entorno de cómputo principal configurado, puede agregar `export TOWER_COMPUTE_ENV=<compute-env-name>` a su archivo de configuración para establecer uno predeterminado.
    Cualquier valor de configuración puede ser reemplazado en la línea de comandos pasando el flag explícitamente (por ejemplo, `--compute-env other-env`).
    Consulte la [referencia del CLI tw](https://docs.seqera.io/platform/latest/cli/reference) para ver la lista completa de opciones y variables de entorno.

### 1.4. Explorar su workspace desde el CLI

Liste los workspaces a los que tiene acceso:

```bash
tw workspaces list
```

??? success "Salida del comando"

    ```console
    Available workspaces:
    ID              | Name            | Full Name                | Visibility
    --------------- | --------------- | ------------------------ | ----------
    <workspace-id>  | my-workspace    | my-org/my-workspace      | PRIVATE
    ```

Vea las ejecuciones en su workspace, incluyendo la ejecución de nf-core/rnaseq que acaba de lanzar:

```bash
tw runs list
```

??? success "Salida del comando"

    ```console
    Pipeline runs for my-org/my-workspace workspace:
    ID        | Status   | Name          | Pipeline               | Run name
    --------- | -------- | ------------- | ---------------------- | --------
    <run-id>  | RUNNING  | ...           | nf-core/rnaseq         | happy_curie
    ```

La misma ejecución que está monitoreando en la interfaz web es visible aquí.

!!! note "Nota"

    Dado que `TOWER_WORKSPACE_ID` está configurado en `.seqera_config`, puede omitir `--workspace` en todos los comandos `tw`.
    Sin la configuración, tendría que pasarlo explícitamente:

    ```bash
    tw runs list --workspace <org>/<workspace>
    ```

Todo lo que es visible en la interfaz web es accesible desde el CLI.

### 1.5. Ejecutar nf-core/rnaseq desde el CLI

El pipeline que agregó a su workspace en la [Parte 1](./01_run_with_seqera.md) está disponible por nombre en el CLI.
Ejecútelo con el perfil `test`:

```bash
tw launch nf-core-rnaseq --profile test
```

??? success "Salida del comando"

    ```console
    Launching pipeline nf-core-rnaseq
    Run name: focused_einstein
    https://cloud.seqera.io/orgs/my-org/workspaces/my-workspace/watch/<run-id>
    ```

Abra el enlace en su navegador y confirme que la ejecución aparece en el panel **Runs**.

Una vez que pueda verla en ejecución, habrá confirmado que el CLI y la interfaz web son dos vistas del mismo workspace.

!!! note "Nota"

    También puede pasar una URL completa de GitHub directamente a `tw launch` sin agregar el pipeline a un workspace primero.
    Sin embargo, agregar el pipeline explícitamente antes de ejecutarlo es generalmente mejor: guarda la configuración del pipeline para ejecuciones futuras, lo hace disponible por nombre y lo hace visible para todos los miembros del workspace en el Launchpad.

    Es posible agregar un pipeline a un workspace directamente desde la línea de comandos usando `tw`.
    La siguiente sección muestra cómo hacerlo con el pipeline nf-core/demo.

### Conclusión

Sabe cómo autenticar el CLI `tw`, inspeccionar su workspace y ejecutar un pipeline guardado desde la terminal.

### ¿Qué sigue?

Agregue un nuevo pipeline a su workspace desde la línea de comandos y ejecútelo.

---

## 2. Agregar un nuevo pipeline y ejecutarlo

Cualquier pipeline de Nextflow en GitHub puede agregarse a su workspace con `tw pipelines add`, siempre que tenga un punto de entrada `main.nf` y un `nextflow.config` en su raíz.
nf-core/demo es un buen ejemplo para practicar: ya lo ejecutó en el curso Use nf-core, así que sabe qué hace y qué esperar.

### 2.1. Agregar nf-core/demo a su workspace

Ejecute el siguiente comando para registrar el pipeline en su workspace:

```bash
tw pipelines add \
    --name nf-core-demo \
    https://github.com/nf-core/demo
```

??? success "Salida del comando"

    ```console
    New pipeline 'nf-core-demo' added at my-org/my-workspace workspace
    ```

El pipeline ahora está registrado y aparecerá en el Launchpad.

### 2.2. Verificar que aparece en el Launchpad

Liste los pipelines en su workspace para confirmar que fue agregado:

```bash
tw pipelines list
```

??? success "Salida del comando"

    ```console
    Pipelines at my-org/my-workspace:
    ID   | Name            | Repository
    ---- | --------------- | -----------------------------------------
    ...  | nf-core-demo    | https://github.com/nf-core/demo
    ...  | nf-core-rnaseq  | ...
    ```

Abra su workspace en el navegador y haga clic en **Launchpad** para confirmar que nf-core/demo ahora aparece junto a nf-core/rnaseq.

!!! tip "Consejo"

    También puede agregar pipelines a través de la interfaz web: en la barra lateral izquierda, haga clic en **Launchpad**, luego en **Add pipeline**, y complete el formulario correspondiente.

Haga clic en el botón **Launch** en la entrada de nf-core/demo para abrir su formulario de ejecución.
Verá que los parámetros `input` y `outdir` están resaltados en rojo — son campos obligatorios sin valores predeterminados, porque `tw pipelines add` registra solo el código fuente del pipeline sin preconfigurar ningún parámetro.
Las siguientes dos secciones explican cómo proporcionar esos valores: primero a través del formulario web y luego desde la línea de comandos.

### 2.3. Ejecutar nf-core/demo desde la interfaz web

Con el formulario de ejecución abierto, complete los dos parámetros obligatorios.

Para `input`, ingrese la URL del samplesheet de prueba del perfil test de nf-core/demo.
Puede encontrarla en `conf/test.config` dentro del repositorio del pipeline, que examinó en el curso Use nf-core:

```
https://raw.githubusercontent.com/nf-core/test-datasets/viralrecon/samplesheet/samplesheet_test_illumina_amplicon.csv
```

Para `outdir`, ingrese una ruta de almacenamiento en la nube donde el pipeline pueda escribir sus resultados.
Use el bucket configurado para su workspace, con un subdirectorio para mantener las ejecuciones organizadas:

```
s3://my-bucket/demo-results
```

Una vez que ambos campos estén completos, haga clic en el botón azul **Launch**.

La ejecución aparece en el panel **Runs** y debería completarse en unos minutos con el conjunto de datos de prueba.
Haga clic en la ejecución para explorar la tabla de tareas y los reportes de ejecución.

### 2.4. Ejecutar nf-core/demo desde el CLI

A diferencia de `nextflow run`, el comando `tw launch` no acepta flags de parámetros individuales como `--input` o `--outdir`.
Los parámetros deben proporcionarse a través de un archivo en formato YAML o JSON, pasado con `--params-file`.
Esto fomenta la reproducibilidad: un archivo de parámetros guardado documenta exactamente qué valores se usaron en una ejecución, facilitando repetir o compartir una configuración de ejecución.

Cree un archivo de parámetros en su directorio de trabajo:

```bash
touch params.yaml
```

Ábralo en el editor y agregue la ruta de salida:

```yaml title="params.yaml"
outdir: "s3://my-bucket/demo-results-cli"
```

Ahora puede ejecutar el pipeline usando el perfil `test` (que proporciona el samplesheet de `input`) y el archivo de parámetros (que proporciona `outdir`):

```bash
tw launch nf-core-demo --profile test --params-file params.yaml
```

??? success "Salida del comando"

    ```console
    Launching pipeline nf-core-demo
    Run name: focused_feynman
    https://cloud.seqera.io/orgs/my-org/workspaces/my-workspace/watch/<run-id>
    ```

Abra el enlace para confirmar que la ejecución aparece en el panel **Runs**.

!!! tip "Consejo"

    Puede incluir el archivo de parámetros durante el paso de configuración inicial si desea establecer algunos valores predeterminados, así como algunas propiedades adicionales para coincidir con lo que hicimos antes a través del formulario web:

    ```bash
    tw pipelines add \
      --name nf-core-demo-gg \
      --description "Demo pipeline with defaults for testing" \
      --params-file params.yaml \
      --profile test \
      https://github.com/nf-core/demo
    ```

### Conclusión

Sabe cómo agregar cualquier pipeline de Nextflow alojado en GitHub a su workspace y ejecutarlo, tanto desde la interfaz web completando los parámetros manualmente, como desde el CLI `tw` combinando un perfil con un archivo de parámetros.

---

## Resumen

En esta parte aprendió a:

- Autenticar el CLI `tw` y ejecutar un pipeline guardado desde la terminal
- Agregar un nuevo pipeline desde GitHub usando el CLI y verificar que aparece en el Launchpad
- Ejecutar un pipeline desde la interfaz web de Seqera completando los parámetros obligatorios manualmente
- Ejecutar un pipeline desde el CLI usando un perfil de Nextflow y un archivo de parámetros
