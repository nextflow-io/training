# Parte 3: Perfilado y optimización de recursos

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Traducción asistida por IA - [más información y sugerencias](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

ESTO ES UN MARCADOR DE POSICIÓN

!!!note "Nota"

    Este módulo de entrenamiento está en proceso de rediseño.

---

POR HACER

### 1.1. Ejecutar el flujo de trabajo para generar un reporte de utilización de recursos

Para que Nextflow genere el reporte automáticamente, simplemente agregue `-with-report <nombre_archivo>.html` a su línea de comando.

```bash
nextflow run main.nf -profile my_laptop -with-report report-config-1.html
```

El reporte es un archivo html, que puede descargar y abrir en su navegador. También puede hacer clic derecho sobre él en el explorador de archivos de la izquierda y hacer clic en `Show preview` para visualizarlo en VS Code.

Tómese unos minutos para revisar el reporte y ver si puede identificar algunas oportunidades para ajustar los recursos.
Asegúrese de hacer clic en las pestañas que muestran los resultados de utilización como un porcentaje de lo que fue asignado.
Hay [documentación](https://www.nextflow.io/docs/latest/reports.html) disponible que describe todas las características.

<!-- TODO: insert images -->

Una observación es que `GATK_JOINTGENOTYPING` parece tener mucha demanda de CPU, lo cual tiene sentido ya que realiza muchos cálculos complejos.
Así que podríamos intentar aumentar eso y ver si reduce el tiempo de ejecución.

Sin embargo, parece que nos hemos excedido con las asignaciones de memoria; todos los procesos solo están usando una fracción de lo que les estamos dando.
Deberíamos reducir eso y ahorrar algunos recursos.

### 1.2. Ajustar las asignaciones de recursos para un proceso específico

Podemos especificar asignaciones de recursos para un proceso dado usando el selector de proceso `withName`.
La sintaxis se ve así cuando está por sí sola en un bloque process:

```groovy title="Sintaxis"
process {
    withName: 'GATK_JOINTGENOTYPING' {
        cpus = 4
    }
}
```

Agreguemos eso al bloque process existente en el archivo `nextflow.config`.

```groovy title="nextflow.config" linenums="11"
process {
    // valores predeterminados para todos los procesos
    cpus = 2
    memory = 2.GB
    // asignaciones para un proceso específico
    withName: 'GATK_JOINTGENOTYPING' {
        cpus = 4
    }
}
```

Con eso especificado, la configuración predeterminada se aplicará a todos los procesos **excepto** el proceso `GATK_JOINTGENOTYPING`, que es un caso especial que obtiene mucho más CPU.
Esperemos que eso tenga un efecto.

### 1.3. Ejecutar nuevamente con la configuración modificada

Ejecutemos el flujo de trabajo nuevamente con la configuración modificada y con la bandera de reporte activada, pero note que estamos dando al reporte un nombre diferente para poder diferenciarlos.

```bash
nextflow run main.nf -profile my_laptop -with-report report-config-2.html
```

Una vez más, probablemente no notará una diferencia sustancial en el tiempo de ejecución, porque esta es una carga de trabajo tan pequeña y las herramientas pasan más tiempo en tareas auxiliares que en realizar el trabajo 'real'.

Sin embargo, el segundo reporte muestra que nuestra utilización de recursos está más equilibrada ahora.

<!-- **TODO: screenshots?** -->

Como puede ver, este enfoque es útil cuando sus procesos tienen diferentes requisitos de recursos. Le permite dimensionar correctamente las asignaciones de recursos que configura para cada proceso basándose en datos reales, no en conjeturas.

!!!note "Nota"

    Esto es solo una pequeña muestra de lo que puede hacer para optimizar su uso de recursos.
    Nextflow tiene una [lógica de reintento dinámico](https://www.nextflow.io/docs/latest/process.html#dynamic-task-resources) realmente útil incorporada para reintentar trabajos que fallan debido a limitaciones de recursos.
    Además, Seqera Platform ofrece herramientas impulsadas por IA para optimizar sus asignaciones de recursos automáticamente también.

    Cubriremos ambos enfoques en una próxima parte de este curso de entrenamiento.

Dicho esto, puede haber algunas restricciones sobre lo que puede (o debe) asignar dependiendo del ejecutor de cómputo y la infraestructura de cómputo que esté usando. Por ejemplo, su clúster puede requerir que se mantenga dentro de ciertos límites que no se aplican cuando está ejecutando en otro lugar.
