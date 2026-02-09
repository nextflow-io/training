# GitHub Codespaces

GitHub Codespaces es una plataforma basada en la web que nos permite proporcionar un entorno preconfigurado para capacitación, respaldado por máquinas virtuales en la nube.
La plataforma es operada por Github (que es propiedad de Microsoft), y es accesible de forma gratuita (con cuotas de uso) para cualquier persona con una cuenta de Github.

!!! warning "Advertencia"

    Las cuentas asociadas a organizaciones pueden estar sujetas a ciertas restricciones adicionales.
    Si ese es tu caso, es posible que necesites usar una cuenta personal independiente, o usar una instalación local en su lugar.

## Crear una cuenta de GitHub

Puedes crear una cuenta gratuita de GitHub desde la [página principal de GitHub](https://github.com/).

## Iniciar tu GitHub Codespace

Una vez que hayas iniciado sesión en GitHub, abre este enlace en tu navegador para abrir el entorno de capacitación de Nextflow: <https://codespaces.new/nextflow-io/training?quickstart=1&ref=master>

Alternativamente, puedes hacer clic en el botón que se muestra a continuación, que se repite en cada curso de capacitación (típicamente en la página de Orientación).

[![Open in GitHub Codespaces](https://github.com/codespaces/badge.svg)](https://codespaces.new/nextflow-io/training?quickstart=1&ref=master)

Deberías ver una página donde puedes crear un nuevo GitHub Codespace:

![Create a GitHub Codespace](img/codespaces_create.png)

### Configuración

Para uso general, no deberías necesitar configurar nada.
A menos que se especifique lo contrario en el curso que estás comenzando, simplemente puedes hacer clic en el botón principal para continuar.

Sin embargo, es posible personalizar el entorno haciendo clic en el botón "Change options".

??? info "Opciones de configuración"

    Si haces clic en el botón "Change options", tendrás la opción de personalizar lo siguiente:

    #### Branch

    Esto te permite seleccionar una versión diferente de los materiales de capacitación.
    La rama `master` generalmente contiene correcciones de errores y materiales que han sido desarrollados y aprobados recientemente pero que aún no han sido publicados en el sitio web.
    Otras ramas contienen trabajo en progreso que puede no ser completamente funcional.

    #### Machine type

    Esto te permite personalizar la máquina virtual que usarás para trabajar en la capacitación.

    Usar una máquina con más núcleos te permite aprovechar mejor la capacidad de Nextflow para paralelizar la ejecución del workflow.
    Sin embargo, consumirá tu asignación de cuota gratuita más rápido, por lo que no recomendamos cambiar esta configuración a menos que se indique en las instrucciones del curso que planeas tomar.

    Consulta 'Cuotas de GitHub Codespaces' más abajo para más detalles sobre las cuotas.

### Tiempo de inicio

Abrir un nuevo entorno de GitHub Codespaces por primera vez puede tomar varios minutos, porque el sistema tiene que configurar tu máquina virtual, así que no te preocupes si hay un tiempo de espera.
Sin embargo, no debería tomar más de cinco minutos.

## Navegar por la interfaz de capacitación

Una vez que tu GitHub Codespaces se haya cargado, deberías ver algo similar a lo siguiente (que puede abrirse en modo claro dependiendo de las preferencias de tu cuenta):

![GitHub Codespaces welcome](img/codespaces_welcome.png)

Esta es la interfaz del IDE VSCode, una aplicación popular de desarrollo de código que recomendamos usar para el desarrollo con Nextflow.

- **El editor principal** es donde se abrirán el código de Nextflow y otros archivos de texto. Aquí es donde editarás el código. Cuando abras el codespace, esto te mostrará una vista previa del archivo `README.md`.
- **La terminal** debajo del editor principal te permite ejecutar comandos. Aquí es donde ejecutarás todas las líneas de comando dadas en las instrucciones del curso.
- **La barra lateral** te permite personalizar tu entorno y realizar tareas básicas (copiar, pegar, abrir archivos, buscar, git, etc.). Por defecto está abierta en el explorador de archivos, que te permite navegar por el contenido del repositorio. Hacer clic en un archivo en el explorador lo abrirá dentro de la ventana del editor principal.

Puedes ajustar las proporciones relativas de los paneles de la ventana como desees.

<!-- TODO (future) Link to development best practices side quest? -->

## Otras notas sobre el uso de GitHub Codespaces

### Reanudar una sesión

Una vez que hayas creado un entorno, puedes reanudarlo o reiniciarlo fácilmente y continuar desde donde lo dejaste.
Tu entorno expirará después de 30 minutos de inactividad y guardará tus cambios hasta por 2 semanas.

Puedes reabrir un entorno desde <https://github.com/codespaces/>.
Se listarán los entornos anteriores.
Haz clic en una sesión para reanudarla.

![List GitHub Codespace sessions](img/codespaces_list.png)

Si has guardado la URL de tu entorno anterior de GitHub Codespaces, simplemente puedes abrirla en tu navegador.
Alternativamente, haz clic en el mismo botón que usaste para crearlo en primer lugar:

[![Open in GitHub Codespaces](https://github.com/codespaces/badge.svg)](https://codespaces.new/nextflow-io/training?quickstart=1&ref=master)

Deberías ver la sesión anterior, la opción predeterminada es reanudarla:

![Resume a GitHub Codespace](img/codespaces_resume.png)

### Guardar archivos en tu máquina local

Para guardar cualquier archivo desde el panel del explorador, haz clic derecho en el archivo y selecciona `Download`.

### Administrar las cuotas de GitHub Codespaces

GitHub Codespaces te da hasta 15 GB-mes de almacenamiento por mes, y 120 horas-núcleo por mes.
Esto es equivalente a alrededor de 60 horas de tiempo de ejecución del entorno predeterminado usando el espacio de trabajo estándar (2 núcleos, 8 GB de RAM y 32 GB de almacenamiento).

Puedes crearlos con más recursos (ver explicación arriba), pero esto consumirá tu uso gratuito más rápido y tendrás menos horas de acceso a este espacio.
Por ejemplo, si seleccionas una máquina de 4 núcleos en lugar del predeterminado de 2 núcleos, tu cuota se agotará en la mitad del tiempo.

Opcionalmente, puedes comprar acceso a más recursos.

Para más información, consulta la documentación de GitHub:
[About billing for GitHub Codespaces](https://docs.github.com/en/billing/managing-billing-for-your-products/managing-billing-for-github-codespaces/about-billing-for-github-codespaces)
