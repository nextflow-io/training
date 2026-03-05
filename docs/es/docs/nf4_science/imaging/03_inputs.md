# Parte 3: Organización de entradas

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Traducción asistida por IA - [más información y sugerencias](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

En la Parte 2, ejecutamos molkart con múltiples parámetros en la línea de comandos.
Ahora aprenderemos dos enfoques mejores para gestionar entradas: **archivos de parámetros** y **hojas de muestras**.

## 1. Uso de archivos de parámetros

### 1.1. El problema con líneas de comandos largas

Recordemos nuestro comando de la Parte 2:

```bash
nextflow run ./molkart \
  --input 'data/samplesheet.csv' \
  --mindagap_tilesize 90 \
  --mindagap_boxsize 7 \
  --mindagap_loopnum 100 \
  --clahe_pyramid_tile 368 \
  --segmentation_method "cellpose" \
  --outdir results
```

Esto funciona, pero es difícil de reproducir, compartir o modificar.
¿Qué pasa si necesita ejecutar el mismo análisis nuevamente el próximo mes?
¿Qué pasa si un colaborador quiere usar exactamente su configuración?

### 1.2. Solución: Use un archivo de parámetros

Cree un archivo llamado `params.yaml`:

```yaml title="params.yaml"
input: "data/samplesheet.csv"
outdir: "results"
mindagap_tilesize: 90
mindagap_boxsize: 7
mindagap_loopnum: 100
clahe_pyramid_tile: 368
segmentation_method: "cellpose"
```

Ahora su comando se convierte en:

```bash
nextflow run ./molkart -params-file params.yaml -resume
```

¡Eso es todo! El archivo de parámetros documenta su configuración exacta y facilita la reejecución o el compartir.

### 1.3. Anulación de parámetros

Todavía puede anular parámetros específicos desde la línea de comandos:

```bash
nextflow run ./molkart -params-file params.yaml --segmentation_method "stardist" --outdir stardist_results -resume
```

La línea anterior cambia el `segmentation_method` a `stardist` y el nombre de `--outdir` a `stardist_results` en lugar de los parámetros en el archivo `params.yaml`.
Además, puede ver que la bandera `-resume` nos permitió reutilizar los resultados de preprocesamiento de la ejecución anterior, ahorrando tiempo.
Puede usar este patrón para probar rápidamente diferentes variaciones del pipeline.

### Conclusión

Los archivos de parámetros hacen que sus análisis sean reproducibles y fáciles de compartir.
Úselos para cualquier trabajo de análisis real.

### ¿Qué sigue?

Aprenda cómo las hojas de muestras organizan información sobre múltiples muestras.

---

## 2. El patrón de hoja de muestras

### 2.1. ¿Qué es una hoja de muestras?

Una hoja de muestras es un archivo CSV que describe sus muestras de entrada.
Cada fila es una muestra, y las columnas especifican los archivos y metadatos para esa muestra.

Veamos la hoja de muestras que hemos estado usando:

```bash
cat data/samplesheet.csv
```

```csv title="samplesheet.csv"
sample,nuclear_image,spot_table,membrane_image
mem_only,https://raw.githubusercontent.com/nf-core/test-datasets/molkart/test_data/input_data/nuclear.tiff,https://raw.githubusercontent.com/nf-core/test-datasets/molkart/test_data/input_data/spots.txt,https://raw.githubusercontent.com/nf-core/test-datasets/molkart/test_data/input_data/membrane.tiff
```

Las columnas son:

- `sample`: Identificador único de la muestra
- `nuclear_image`: Imagen de tinción nuclear (TIFF)
- `spot_table`: Puntos de transcripción (TXT)
- `membrane_image`: Imagen de tinción de membrana (TIFF, opcional)

### 2.2. Rutas de archivos

Las hojas de muestras aceptan múltiples tipos de rutas:

- **URLs**: Nextflow descarga automáticamente (como se muestra arriba)
- **Rutas locales**: `data/nuclear.tiff` o `/absolute/path/to/nuclear.tiff`
- **Almacenamiento en la nube**: `s3://bucket/nuclear.tiff`, `gs://bucket/nuclear.tiff`, `az://container/nuclear.tiff`

Puede mezclar tipos de rutas en la misma hoja de muestras.

### 2.3. Creación de su propia hoja de muestras

Primero, descarguemos los archivos de datos de prueba localmente:

```bash
cd /workspaces/training/nf4-science/imaging/data
curl -O https://raw.githubusercontent.com/nf-core/test-datasets/molkart/test_data/input_data/nuclear.tiff
curl -O https://raw.githubusercontent.com/nf-core/test-datasets/molkart/test_data/input_data/spots.txt
curl -O https://raw.githubusercontent.com/nf-core/test-datasets/molkart/test_data/input_data/membrane.tiff
cd ..
```

Ahora modifiquemos la hoja de muestras para hacer referencia a estos archivos locales:

```csv title="samplesheet.csv"
sample,nuclear_image,spot_table,membrane_image
mem_only,data/nuclear.tiff,data/spots.txt,data/membrane.tiff
```

!!! warning "Advertencia"

    Observe que las rutas en la hoja de muestras son relativas a donde **ejecuta** Nextflow, no a donde se encuentra la hoja de muestras.

Finalmente, ejecutemos nf-core/molkart una vez más con la hoja de muestras con rutas de archivos locales:

`nextflow run ./molkart -params-file params.yaml -resume`

Como puede ver, Nextflow ejecuta esta ejecución de manera similar a cuando los archivos se descargaron desde Github. Esta es una de las grandes características de Nextflow, organiza los datos adecuadamente para usted, independientemente de dónde se encuentren.

### Conclusión

Las hojas de muestras organizan conjuntos de datos con múltiples muestras de una manera que le permite definir explícitamente sus metadatos junto con las rutas de archivos.
La mayoría de los pipelines de nf-core usan este patrón.

### ¿Qué sigue?

Ahora que hemos cubierto las entradas, exploremos cómo configurar pipelines de Nextflow para diferentes entornos de cómputo.
