# Parte 1: Lanzar pipelines desde la interfaz web

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Traducción asistida por IA - [más información y sugerencias](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

En esta parte del curso Scale with Seqera, configurará el acceso a Seqera Platform y lanzará un pipeline a escala de producción desde la interfaz web.

Asegúrese de que su directorio de trabajo esté configurado en `seqera-scale/` según las instrucciones de la página de [Primeros pasos](./00_orientation.md).

---

## 1. Primeros pasos con Seqera

Seqera proporciona una plataforma integral para lanzar, monitorear y gestionar pipelines de Nextflow.
Esta sección le guía a través del proceso de registro y orientación antes de ejecutar su primer pipeline.

### 1.1. Regístrese para obtener una cuenta gratuita

Vaya a [cloud.seqera.io](https://cloud.seqera.io) y cree una cuenta gratuita.
Puede registrarse usando su dirección de correo electrónico, GitHub o credenciales de Google.

Una cuenta gratuita le ofrece:

- **Espacio de trabajo personal**: su propio espacio para agregar pipelines, configurar entornos de cómputo y gestionar ejecuciones
- **Acceso al Community Showcase**: una colección curada de pipelines de nf-core y de la comunidad con configuraciones predefinidas y datos de ejemplo

Consulte la [documentación de Seqera](https://docs.seqera.io) para obtener una descripción completa de los niveles de cuenta y las funciones disponibles.

### 1.2. Explore el Community Showcase

Antes de lanzar sus propios pipelines, tómese unos minutos para explorar el Community Showcase.
Le ofrece una vista previa realista de cómo se ve la plataforma con pipelines y datos reales.

1. Inicie sesión en [cloud.seqera.io](https://cloud.seqera.io).
2. En la barra lateral izquierda, haga clic en **Showcase**.
3. Explore los pipelines disponibles — reconocerá varios pipelines de nf-core del curso Use nf-core.
4. Haga clic en un pipeline para ver su configuración y ajustes de lanzamiento.
5. Haga clic en **Runs** para explorar historiales de ejecuciones de ejemplo, incluyendo detalles a nivel de tarea e informes de ejecuciones anteriores.

Esta es una vista de solo lectura, pero le muestra cómo funciona la interfaz antes de ejecutar cualquier cosa por su cuenta.

### 1.3. Acceda a un espacio de trabajo con cómputo

Para lanzar pipelines se requiere un espacio de trabajo con un entorno de cómputo configurado.

Seqera admite dos formas de proporcionar cómputo:

- **Conecte su propia infraestructura**: AWS, Azure, Google Cloud y planificadores HPC (SLURM, LSF, PBS, entre otros).
  Consulte la [documentación de entornos de cómputo](https://docs.seqera.io) para las guías de configuración.
- **Seqera Compute**: un servicio administrado que proporciona entornos de cómputo preaprovisionados en AWS, con un costo, sin necesidad de configurar una cuenta en la nube.
  Puede activarlo directamente desde la configuración de su espacio de trabajo.

**Capacitación en grupo:**
Si está participando en una sesión de capacitación grupal, es posible que haya sido agregado a una organización y espacio de trabajo que ya tiene cómputo configurado.
Su instructor le proporcionará el nombre de la organización, el nombre del espacio de trabajo y cualquier otro detalle que necesite.

**Trabajo independiente:**
Si está siguiendo esta capacitación por su cuenta, deberá configurar un entorno de cómputo en su espacio de trabajo personal usando una de las opciones anteriores.
Los créditos gratuitos para probar Seqera Compute están [disponibles bajo solicitud](https://seqera.io/platform/compute/).

!!! note "Nota"

    El resto de este curso asume que tiene acceso a un espacio de trabajo con un entorno de cómputo configurado.
    Si está en una sesión de capacitación grupal, su instructor confirmará qué espacio de trabajo y entorno de cómputo utilizar.

### Conclusión

Tiene una cuenta de Seqera, ha explorado el Community Showcase y puede acceder a un espacio de trabajo con cómputo.

### ¿Qué sigue?

Lanzar un pipeline de RNA-seq a escala de producción desde la interfaz web de Seqera Cloud.

---

## 2. Lanzar nf-core/rnaseq desde la interfaz web

Como se vio en Use nf-core, el pipeline nf-core/rnaseq es un pipeline curado por la comunidad para el análisis de datos de secuenciación de RNA en bulk.

En esta sección, agregará el pipeline a su espacio de trabajo, lanzará una ejecución y monitoreará su progreso.

### 2.1. Agregar el pipeline a su espacio de trabajo

Convenientemente, nf-core/rnaseq forma parte de una colección curada de pipelines que se pueden agregar a su espacio de trabajo en pocos clics a través del servicio Seqera Pipelines.

_Más adelante en este curso le mostraremos cómo agregar sus propios pipelines._

1. Navegue a [**Seqera Pipelines**](https://seqera.io/pipelines) para explorar la colección de la comunidad.
2. Busque `rnaseq` y seleccione **nf-core/rnaseq**.
3. Haga clic en **Launch Pipeline** o desplácese hasta el final de la página hasta la sección **Launch Pipeline**.
4. Asegúrese de haber iniciado sesión y seleccione los valores apropiados en los menús desplegables **Organizations**, **Workspace** y **Compute Environment**.
   **Consejo para grupos:** Si está usando un espacio de trabajo compartido, agregue un identificador único (como su nombre de usuario) al nombre del pipeline.
5. Haga clic en **Add pipeline to your Seqera account**.

Aparecerá un cuadro con el mensaje: **Pipeline added: View Pipeline**.
Al hacer clic en el enlace, accederá a la entrada del pipeline en su launchpad.

El pipeline ahora aparece en el panel **Launchpad** de su espacio de trabajo y está listo para lanzarse.

### 2.2. Lanzar el pipeline

Haga clic en el botón **Launch** del pipeline, ya sea en el panel **Launchpad** o en la página de detalles del pipeline.
Esto abre la interfaz de configuración.

El pipeline ya está configurado con el perfil `test`, por lo que los datos de entrada, el directorio de salida y la referencia del genoma están completados previamente.
Por ahora puede ignorar el resto de los parámetros y la configuración avanzada.

Haga clic en el botón azul **Launch** para iniciar la ejecución.

### 2.3. Monitorear la ejecución

Después de lanzar, será llevado al panel **Runs** de su pipeline.

La vista de ejecución muestra:

- **Status**: estado actual de la ejecución (submitted, running, succeeded, failed)
- **Command line**: el comando exacto `nextflow run` que la plataforma construyó y envió
- **Parameters**: todos los valores de parámetros utilizados en esta ejecución
- **Tasks**: una tabla de cada llamada a proceso, con estado, duración y uso de recursos

Haga clic en cualquier fila de tarea para inspeccionar sus detalles de ejecución, incluyendo:

- El script `.command.sh` que se ejecutó
- Registros de stdout y stderr
- Métricas de CPU, memoria y E/S

La pestaña **Reports** mostrará un informe de MultiQC una vez que la ejecución se complete, agregando métricas de control de calidad de todas las muestras.

Esto tomará un tiempo en ejecutarse, así que continuaremos por ahora y volveremos más tarde para revisar las salidas y demás.

### Conclusión

Sabe cómo agregar un pipeline a un espacio de trabajo de Seqera, configurar y lanzar una ejecución, y monitorear la ejecución a escala.

### ¿Qué sigue?

Continúe con la [Parte 2](./02_launch_from_cli.md), donde aprenderá cómo hacer todo esto desde la línea de comandos usando el CLI `tw`.

---

## Resumen

En esta parte aprendió a:

- Registrarse en Seqera y explorar el Community Showcase
- Agregar un pipeline del catálogo curado, lanzar una ejecución a escala de producción y monitorear su ejecución
