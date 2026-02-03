# Parte 1: Más Contenedores

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Traducción asistida por IA - [más información y sugerencias](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

[TODO]

---

## 1. Cómo encontrar o crear imágenes de contenedores

Algunos desarrolladores de software proporcionan imágenes de contenedores para su software que están disponibles en registros de contenedores como Docker Hub, pero muchos no lo hacen.
En esta sección opcional, le mostraremos dos formas de obtener una imagen de contenedor para las herramientas que desea usar en sus pipelines de Nextflow: usando Seqera Containers y construyendo la imagen de contenedor usted mismo.

Obtendrá/construirá una imagen de contenedor para el paquete pip `quote`, que se usará en el ejercicio al final de esta sección.

### 1.1. Obtener una imagen de contenedor de Seqera Containers

Seqera Containers es un servicio gratuito que construye imágenes de contenedores para herramientas instalables mediante pip y conda (incluido bioconda).
Navegue a [Seqera Containers](https://www.seqera.io/containers/) y busque el paquete pip `quote`.

![Seqera Containers](img/seqera-containers-1.png)

Haga clic en "+Add" y luego en "Get Container" para solicitar una imagen de contenedor para el paquete pip `quote`.

![Seqera Containers](img/seqera-containers-2.png)

Si esta es la primera vez que se construye un contenedor comunitario para esta versión del paquete, puede tomar unos minutos completarse.
Haga clic para copiar el URI (por ejemplo, `community.wave.seqera.io/library/pip_quote:ae07804021465ee9`) de la imagen de contenedor que se creó para usted.

Ahora puede usar la imagen de contenedor para ejecutar el comando `quote` y obtener una frase aleatoria de Grace Hopper.

```bash
docker run --rm community.wave.seqera.io/library/pip_quote:ae07804021465ee9 quote "Grace Hopper"
```

Salida:

```console title="Salida"
Humans are allergic to change. They love to say, 'We've always done it
this way.' I try to fight that. That's why I have a clock on my wall
that runs counter-clockwise.
```

### 1.2. Construir la imagen de contenedor usted mismo

Usemos algunos detalles de construcción del sitio web de Seqera Containers para construir nosotros mismos la imagen de contenedor para el paquete pip `quote`.
Regrese al sitio web de Seqera Containers y haga clic en el botón "Build Details".

El primer elemento que veremos es el `Dockerfile`, un tipo de archivo de script que contiene todos los comandos necesarios para construir la imagen de contenedor.
Hemos agregado algunos comentarios explicativos al Dockerfile a continuación para ayudarle a entender qué hace cada parte.

```Dockerfile title="Dockerfile"
# Comenzar desde la imagen base de docker de micromamba
FROM mambaorg/micromamba:1.5.10-noble
# Copiar el archivo conda.yml dentro del contenedor
COPY --chown=$MAMBA_USER:$MAMBA_USER conda.yml /tmp/conda.yml
# Instalar varias utilidades para que Nextflow las use y los paquetes en el archivo conda.yml
RUN micromamba install -y -n base -f /tmp/conda.yml \
    && micromamba install -y -n base conda-forge::procps-ng \
    && micromamba env export --name base --explicit > environment.lock \
    && echo ">> CONDA_LOCK_START" \
    && cat environment.lock \
    && echo "<< CONDA_LOCK_END" \
    && micromamba clean -a -y
# Ejecutar el contenedor como el usuario root
USER root
# Establecer la variable de entorno PATH para incluir el directorio de instalación de micromamba
ENV PATH="$MAMBA_ROOT_PREFIX/bin:$PATH"
```

El segundo elemento que veremos es el archivo `conda.yml`, que contiene la lista de paquetes que deben instalarse en la imagen de contenedor.

```conda.yml title="conda.yml"
channels:
- conda-forge
- bioconda
dependencies:
- pip
- pip:
  - quote==3.0.0 #
```

Copie el contenido de estos archivos en los stubs ubicados en el directorio `containers/build`, luego ejecute el siguiente comando para construir la imagen de contenedor usted mismo.

!!! note "Nota"

    Usamos la bandera `-t quote:latest` para etiquetar la imagen de contenedor con el nombre `quote` y la etiqueta `latest`.
    Podremos usar esta etiqueta para referirnos a la imagen de contenedor al ejecutarla en este sistema.

```bash
docker build -t quote:latest containers/build
```

Después de que haya terminado de construirse, puede ejecutar la imagen de contenedor que acaba de construir.

```bash
docker run --rm quote:latest quote "Margaret Oakley Dayhoff"
```

### Conclusión

Ha aprendido dos formas diferentes de obtener una imagen de contenedor para una herramienta que desea usar en sus pipelines de Nextflow: usando Seqera Containers y construyendo la imagen de contenedor usted mismo.

### ¿Qué sigue?

Tiene todo lo que necesita para continuar con el [siguiente capítulo](./04_hello_genomics.md) de esta serie de entrenamiento.
También puede continuar con un ejercicio opcional para obtener citas sobre pioneros de la computación/biología usando el contenedor `quote` y mostrarlas usando el contenedor `cowsay`.

---

## 2. Hacer que la vaca cite científicos famosos

Esta sección contiene algunos ejercicios adicionales, para practicar lo que ha aprendido hasta ahora.
Hacer estos ejercicios _no es obligatorio_ para comprender las partes posteriores del entrenamiento, pero proporcionan una forma divertida de reforzar su aprendizaje al descubrir cómo hacer que la vaca cite científicos famosos.

```console title="cowsay-output-Grace-Hopper.txt"
  _________________________________________________
 /                                                 \
| Humans are allergic to change. They love to       |
| say, 'We've always done it this way.' I try to fi |
| ght that. That's why I have a clock on my wall th |
| at runs counter-clockwise.                        |
| -Grace Hopper                                     |
 \                                                 /
  =================================================
                                                 \
                                                  \
                                                    ^__^
                                                    (oo)\_______
                                                    (__)\       )\/\
                                                        ||----w |
                                                        ||     ||
```

### 2.1. Modificar el script `hello-containers.nf` para usar un proceso getQuote

Tenemos una lista de pioneros de la computación y la biología en el archivo `containers/data/pioneers.csv`.
A nivel general, para completar este ejercicio necesitará:

- Modificar el `params.input_file` predeterminado para que apunte al archivo `pioneers.csv`.
- Crear un proceso `getQuote` que use el contenedor `quote` para obtener una cita para cada entrada.
- Conectar la salida del proceso `getQuote` al proceso `cowsay` para mostrar la cita.

Para la imagen de contenedor `quote`, puede usar la que construyó usted mismo en el ejercicio adicional anterior o usar la que obtuvo de Seqera Containers.

!!! tip "Consejo"

    Una buena elección para el bloque `script` de su proceso getQuote podría ser:
        ```groovy
        script:
            def safe_author = author.tokenize(' ').join('-')
            """
            quote "$author" > quote-${safe_author}.txt
            echo "-${author}" >> quote-${safe_author}.txt
            """
        ```

Puede encontrar una solución a este ejercicio en `containers/solutions/hello-containers-4.1.nf`.

### 2.2. Modificar su pipeline de Nextflow para permitir que se ejecute en los modos `quote` y `sayHello`.

Agregue algo de lógica de ramificación a su pipeline para permitir que acepte entradas destinadas tanto para `quote` como para `sayHello`.
Aquí hay un ejemplo de cómo usar una sentencia `if` en un workflow de Nextflow:

```groovy title="hello-containers.nf"
workflow {
    if (params.quote) {
        ...
    }
    else {
        ...
    }
    cowSay(text_ch)
}
```

!!! tip "Consejo"

    Puede usar `new_ch = processName.out` para asignar un nombre al canal de salida de un proceso.

Puede encontrar una solución a este ejercicio en `containers/solutions/hello-containers-4.2.nf`.

### Conclusión

¡Sabe cómo usar contenedores en Nextflow para ejecutar procesos y cómo construir algo de lógica de ramificación en sus pipelines!

### ¿Qué sigue?

¡Celebre, tome un descanso para estirarse y beba algo de agua!

Cuando esté listo, continúe con la Parte 3 de esta serie de entrenamiento para aprender cómo aplicar lo que ha aprendido hasta ahora a un caso de uso de análisis de datos más realista.
