# Devcontainers Locales

Si tienes una instalación local de Docker o estás dispuesto/a a instalar una, la forma más fácil de trabajar localmente con estos materiales es usar la función de devcontainer de Visual Studio Code. Este enfoque proporciona todas las herramientas y dependencias necesarias sin requerir instalación manual.

## Requisitos

Para usar la configuración de devcontainer local, necesitarás:

- [Visual Studio Code](https://code.visualstudio.com/)
- Una instalación local de Docker, por ejemplo:
  - [Docker Desktop](https://docs.docker.com/get-docker/) (para Windows/macOS)
  - [Docker Engine](https://docs.docker.com/engine/install/) (para Linux)
  - [Colima](https://github.com/abiosoft/colima) (alternativa para macOS)
- [Docker Buildx](https://docs.docker.com/build/concepts/overview/#install-buildx) (incluido en Docker Desktop, pero puede necesitar instalación separada con otras configuraciones de Docker)
- [Extensión Dev Containers](https://marketplace.visualstudio.com/items?itemName=ms-vscode-remote.remote-containers) para VS Code

Tu instalación de Docker debe estar ejecutándose antes de intentar abrir el devcontainer.

Para verificar que Docker buildx esté disponible, ejecuta:

```bash
docker buildx version
```

Si este comando falla, necesitarás instalar la extensión buildx antes de continuar.

## Instrucciones de Configuración

Sigue estos pasos para configurar tu entorno local usando devcontainers de VS Code:

### Instala la extensión "Dev Containers" en VS Code

- Abre VS Code
- Ve a Extensiones (Ctrl+Shift+X o Cmd+Shift+X en macOS)
- Busca "Dev Containers"
- Haz clic en "Install"

![Instalando la extensión Dev Containers en VS Code](img/install_extension.png)

### Clona el repositorio:

```bash
git clone https://github.com/nextflow-io/training.git
cd training
```

### Abre el repositorio en VS Code:

- Inicia VS Code
- Selecciona **File -> Open Folder** del menú
- Navega y selecciona la carpeta del repositorio de capacitación que acabas de clonar
- Haz clic en **Open**

### Reabrir en Contenedor

Si VS Code te solicita "Reopen in Container", haz clic en ello. Alternativamente:

- Presiona F1 (o Ctrl+Shift+P / Cmd+Shift+P en macOS)
- Escribe "Dev Containers: Reopen in Container"
- **Importante**: Cuando se te solicite seleccionar una configuración, elige la configuración de devcontainer **local-dev**

![Solicitud de Reopen in Container](img/reopen_prompt.png)

![Seleccionando configuración local](img/select_local_config.png)

Espera a que el contenedor se construya. Esto puede tomar algunos minutos la primera vez, ya que descarga y configura todos los componentes necesarios.

Una vez que el contenedor esté construido y ejecutándose, tendrás un entorno completamente configurado con todas las herramientas necesarias instaladas, incluyendo:

- Java
- Nextflow
- Docker
- Git
- Y todas las demás dependencias requeridas para la capacitación

![VS Code con devcontainer ejecutándose](img/running_container.png)

## Beneficios de Usar Devcontainers

Usar el enfoque de devcontainer ofrece varias ventajas:

- **Consistencia**: Asegura un entorno de desarrollo consistente en diferentes máquinas
- **Simplicidad**: Todas las dependencias están preinstaladas y configuradas
- **Aislamiento**: El entorno de desarrollo está aislado de tu sistema local
- **Reproducibilidad**: Todos los que usan el devcontainer obtienen la misma configuración
- **Sin instalación manual**: No es necesario instalar manualmente Java, Nextflow y otras herramientas

## Verificando tu Entorno

Una vez que tu devcontainer esté ejecutándose, puedes verificar que todo esté configurado correctamente ejecutando:

```bash
nextflow info
```

Esto debería mostrar la versión de Nextflow y la información de ejecución, confirmando que tu entorno está configurado correctamente.

## Solución de Problemas

Si encuentras problemas con la configuración del devcontainer:

1. Asegúrate de que tu instalación de Docker (Docker Desktop, Colima, Docker Engine, etc.) esté ejecutándose antes de abrir el devcontainer
2. Verifica que hayas seleccionado la configuración **local-dev** cuando se te solicite
3. Verifica que Docker buildx esté instalado y funcionando ejecutando `docker buildx version`
4. Si el contenedor no se construye, intenta reconstruirlo ejecutando el comando "Dev Containers: Rebuild Container"
5. Para problemas persistentes, consulta la [guía de solución de problemas de VS Code Dev Containers](https://code.visualstudio.com/docs/devcontainers/troubleshooting)
