# Part 3: Organitzant les entrades

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Traducció assistida per IA - [més informació i suggeriments](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

A la Part 2, vam executar molkart amb múltiples paràmetres a la línia de comandes.
Ara aprendrem dues millors aproximacions per gestionar les entrades: **fitxers de paràmetres** i **fulls de mostres**.

## 1. Utilitzant fitxers de paràmetres

### 1.1. El problema amb les línies de comandes llargues

Recordeu la nostra comanda de la Part 2:

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

Això funciona, però és difícil de reproduir, compartir o modificar.
Què passa si necessiteu executar la mateixa anàlisi de nou el mes que ve?
Què passa si un col·laborador vol utilitzar exactament la vostra configuració?

### 1.2. Solució: Utilitzeu un fitxer de paràmetres

Creeu un fitxer anomenat `params.yaml`:

```yaml title="params.yaml"
input: "data/samplesheet.csv"
outdir: "results"
mindagap_tilesize: 90
mindagap_boxsize: 7
mindagap_loopnum: 100
clahe_pyramid_tile: 368
segmentation_method: "cellpose"
```

Ara la vostra comanda es converteix en:

```bash
nextflow run ./molkart -params-file params.yaml -resume
```

Això és tot! El fitxer de paràmetres documenta la vostra configuració exacta i facilita la reexecució o compartició.

### 1.3. Sobreescrivint paràmetres

Encara podeu sobreescriure paràmetres específics des de la línia de comandes:

```bash
nextflow run ./molkart -params-file params.yaml --segmentation_method "stardist" --outdir stardist_results -resume
```

La línia anterior canvia el `segmentation_method` a `stardist` i el nom de `--outdir` a `stardist_results` en lloc dels paràmetres del fitxer `params.yaml`.
A més, podeu veure que la bandera `-resume` ens va permetre reutilitzar els resultats de preprocessament de l'execució anterior, estalviant temps.
Podeu utilitzar aquest patró per provar ràpidament diferents variacions del pipeline.

### Conclusió

Els fitxers de paràmetres fan que les vostres anàlisis siguin reproduïbles i fàcils de compartir.
Utilitzeu-los per a qualsevol treball d'anàlisi real.

### Què segueix?

Apreneu com els fulls de mostres organitzen la informació sobre múltiples mostres.

---

## 2. El patró del full de mostres

### 2.1. Què és un full de mostres?

Un full de mostres és un fitxer CSV que descriu les vostres mostres d'entrada.
Cada fila és una mostra, i les columnes especifiquen els fitxers i les metadades per a aquesta mostra.

Vegem el full de mostres que hem estat utilitzant:

```bash
cat data/samplesheet.csv
```

```csv title="samplesheet.csv"
sample,nuclear_image,spot_table,membrane_image
mem_only,https://raw.githubusercontent.com/nf-core/test-datasets/molkart/test_data/input_data/nuclear.tiff,https://raw.githubusercontent.com/nf-core/test-datasets/molkart/test_data/input_data/spots.txt,https://raw.githubusercontent.com/nf-core/test-datasets/molkart/test_data/input_data/membrane.tiff
```

Les columnes són:

- `sample`: Identificador únic de la mostra
- `nuclear_image`: Imatge de tinció nuclear (TIFF)
- `spot_table`: Punts de transcripció (TXT)
- `membrane_image`: Imatge de tinció de membrana (TIFF, opcional)

### 2.2. Rutes de fitxers

Els fulls de mostres accepten múltiples tipus de rutes:

- **URLs**: Nextflow les descarrega automàticament (com es mostra a dalt)
- **Rutes locals**: `data/nuclear.tiff` o `/absolute/path/to/nuclear.tiff`
- **Emmagatzematge al núvol**: `s3://bucket/nuclear.tiff`, `gs://bucket/nuclear.tiff`, `az://container/nuclear.tiff`

Podeu barrejar tipus de rutes al mateix full de mostres.

### 2.3. Creant el vostre propi full de mostres

Primer, descarreguem els fitxers de dades de prova localment:

```bash
cd /workspaces/training/nf4-science/imaging/data
curl -O https://raw.githubusercontent.com/nf-core/test-datasets/molkart/test_data/input_data/nuclear.tiff
curl -O https://raw.githubusercontent.com/nf-core/test-datasets/molkart/test_data/input_data/spots.txt
curl -O https://raw.githubusercontent.com/nf-core/test-datasets/molkart/test_data/input_data/membrane.tiff
cd ..
```

Ara modifiquem el full de mostres per fer referència a aquests fitxers locals:

```csv title="samplesheet.csv"
sample,nuclear_image,spot_table,membrane_image
mem_only,data/nuclear.tiff,data/spots.txt,data/membrane.tiff
```

!!! Warning "Advertència"

    Noteu que les rutes al full de mostres són relatives a on **executeu** Nextflow, no a on es troba el full de mostres.

Finalment, executem nf-core/molkart una vegada més amb el full de mostres amb rutes de fitxers locals:

`nextflow run ./molkart -params-file params.yaml -resume`

Com podeu veure, Nextflow executa aquesta execució de manera similar a quan els fitxers es van descarregar de Github. Aquesta és una de les grans característiques de Nextflow, prepara les dades adequadament per a vosaltres, independentment d'on estiguin ubicades.

### Conclusió

Els fulls de mostres organitzen conjunts de dades de múltiples mostres d'una manera que us permet definir explícitament les vostres metadades juntament amb les rutes dels fitxers.
La majoria de pipelines nf-core utilitzen aquest patró.

### Què segueix?

Ara que hem cobert les entrades, explorem com configurar els pipelines Nextflow per a diferents entorns informàtics.
