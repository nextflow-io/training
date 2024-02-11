---
description: How to set up a development environment to run Nextflow
---

# Configuración del entorno

Hay dos formas principales de comenzar con el curso de capacitación Nextflow de Seqera.

El primero es instalar los requisitos [localmente](#local-installation), lo cual resulta más cómodo si ya estás familiarizado con Git y Docker, o si trabajas sin conexión.

La segunda es usar [Gitpod](#gitpod), que resulta más sencillo si no conoces Nextflow en profundidad, ya que esta plataforma contiene todos los programas y datos necesarios para la ejecución. Simplemente haz clic en el enlace e inicia sesión con tu cuenta de GitHub para comenzar el tutorial:

[![Abrir en GitPod](https://img.shields.io/badge/Gitpod-%20Open%20in%20Gitpod-908a85?logo=gitpod)](https://gitpod.io/#https:// github.com/nextflow-io/entrenamiento)

## Instalación local

Nextflow se puede usar en cualquier sistema compatible con POSIX (Linux, macOS, subsistema de Windows para Linux, etc.).

#### Requisitos

- Bash
- [Java 11 (o posterior, hasta 18)](https://www.oracle.com/technetwork/java/javase/downloads/index.html)
- Git
- [Docker](https://docs.docker.com/get-docker/)

#### Requisitos opcionales para este tutorial

- [Singularity](https://github.com/sylabs/singularity) 2.5.x (o posterior)
- [Conda](https://conda.io/) 4.5 (o posterior)
- [Graphviz](http://www.graphviz.org/)
- [CLI de AWS](https://aws.amazon.com/cli/)
- Un entorno de AWS Batch configurado

### Descargar Nextflow

Introduce este comando en tu terminal:

```bash
wget -qO- https://get.nextflow.io | bash
```

O, si prefieres curl en lugar de wget:

```bash
curl -s https://get.nextflow.io | bash
```

Luego asegúrate de que el binario descargado sea ejecutable:

```bash
chmod +x nextflow
```

Y coloca el ejecutable `nextflow` en el `$PATH` (por ejemplo, `/usr/local/bin` o `/bin/`)

### Docker

Asegúrese de tener Docker Desktop ejecutándose en tu máquina. Descargue Docker en el siguiente [enlace](https://docs.docker.com/get-docker/).

### Material

Puede ver el material de capacitación aquí: <https://training.nextflow.io/>

Para descargar el material usa el comando:

```bash
git clone https://github.com/nextflow-io/training.git
```

Ejecute ahora `cd` en el directorio `nf-training`.

### Comprobando su instalación

Verifica la instalación correcta de `nextflow` ejecutando el siguiente comando:

```bash
nextflow info
```

Esto debería mostrar la versión actual, el sistema y el tiempo de ejecución.

##Gitpod

Gitpod te proporciona un entorno de desarrollo de Nextflow preconfigurado.

#### Requisitos

- Una cuenta de GitHub
- Un navegador web (Google Chrome o Firefox)
- Conexión a Internet

### Inicio rápido de Gitpod

Para ejecutar Gitpod:

- Haga clic en la siguiente URL: <https://gitpod.io/#https://github.com/nextflow-io/training>
    - Esta es la URL de nuestro repositorio de GitHub, con el prefijo `https://gitpod.io/#`
- Inicie sesión en su cuenta de GitHub (y permita la autorización).

Una vez que haya iniciado sesión, Gitpod debería cargarse (omite la compilación previa si se te solicita).

### Explora tu IDE de Gitpod

Ahora debes ver algo similar a lo siguiente:

![Bienvenida Gitpod](img/gitpod.welcome.png)

- **La barra lateral** te permite personalizar tu entorno Gitpod y realizar tareas básicas (copiar, pegar, abrir archivos, buscar, git, etc.). Haz clic en el botón Explorador para ver qué archivos hay en este repositorio.
- **La terminal** te permite ejecutar todos los programas del repositorio. Por ejemplo, tanto `nextflow` como `docker` están instalados y se pueden ejecutar.
- **La ventana principal** te permite ver y editar archivos. Al hacer clic en un archivo en el explorador, se abrirá dentro de la ventana principal. También deberías ver el navegador del material de formación de nextflow (<https://training.nextflow.io/>).

Para probar que el entorno funciona correctamente, escriba lo siguiente en la terminal:

```bash
nextflow info
```

Esto debería generar la versión de Nextflow y la información del tiempo de ejecución:

```
Version: 22.10.4 build 5836
Created: 09-12-2022 09:58 UTC
System: Linux 5.15.0-47-generic
Runtime: Groovy 3.0.13 on OpenJDK 64-Bit Server VM 17.0.3-internal+0-adhoc..src
Encoding: UTF-8 (UTF-8)
```

### Recursos de Gitpod

- Gitpod le brinda 500 créditos gratuitos al mes, lo que equivale a 50 horas de tiempo de ejecución de entorno gratuito utilizando el espacio de trabajo estándar (hasta 4 núcleos, 8 GB de RAM y 30 GB de almacenamiento).
- También hay una opción de espacio de trabajo que le brinda hasta 8 núcleos, 16 GB de RAM y 50 GB de almacenamiento. Sin embargo, el gran espacio de trabajo usará tus créditos más rápido y tendrás menos horas de acceso en este espacio.
- Gitpod expirará después de 30 minutos de inactividad y guardará sus cambios hasta por 2 semanas (vea la siguiente sección para reabrir una sesión agotada).

Consulte [gitpod.io](https://www.gitpod.io) para obtener más detalles.

### Reapertura de una sesión de Gitpod

Puedes reabrir un entorno de trabajo desde <https://gitpod.io/workspaces>. Busca tu entorno anterior en la lista, luego selecciona los puntos suspensivos (icono de tres puntos) y selecciona Abrir.

Si has guardado la URL para su entorno Gitpod anterior, simplemente puedes abrirla en tu navegador.

Alternativamente, puedes iniciar un nuevo espacio de trabajo en la URL de Gitpod: <https://gitpod.io/#https://github.com/nextflow-io/training>

Si has perdido tu entorno, puedes encontrar los scripts principales utilizados en este tutorial en el directorio `nf-training`.

### Guardar archivos de Gitpod en tu máquina local

Para guardar cualquier archivo desde el panel del explorador, haz clic derecho en el archivo y selecciona `Descargar`.

### Material

Puedes acceder al curso de capacitación en tu navegador desde <https://training.nextflow.io/>

## Seleccionar una versión de Nextflow

De forma predeterminada, Nextflow extraerá la última versión estable en el entorno.

Sin embargo, Nextflow evoluciona constantemente a medida que realizamos mejoras y solucionamos errores.

Los últimos lanzamientos se pueden ver en [GitHub](https://github.com/nextflow-io/nextflow).

Si deseas utilizar una versión específica de Nextflow, puedes configurar la variable `NXF_VER` como se muestra a continuación:

```bash
export NXF_VER=22.04.5
```

!!! Note

    Este taller tutorial requiere `NXF_VER=22.04.0`, o posterior. Esta versión utilizará DSL2 por defecto.

Ejecute `nextflow -version` nuevamente para confirmar que el cambio ha tenido efecto.
