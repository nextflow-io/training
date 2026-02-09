# Parte 3: Organizzare gli input

Nella Parte 2, abbiamo eseguito molkart con più parametri dalla riga di comando.
Ora impareremo due approcci migliori per gestire gli input: **file di parametri** e **samplesheet**.

## 1. Usare i file di parametri

### 1.1. Il problema con le righe di comando lunghe

Ricordiamo il nostro comando dalla Parte 2:

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

Funziona, ma è difficile da riprodurre, condividere o modificare.
Cosa succede se dovete eseguire la stessa analisi di nuovo il mese prossimo?
Cosa succede se un collaboratore vuole usare esattamente le vostre impostazioni?

### 1.2. Soluzione: Usare un file di parametri

Creiamo un file chiamato `params.yaml`:

```yaml title="params.yaml"
input: "data/samplesheet.csv"
outdir: "results"
mindagap_tilesize: 90
mindagap_boxsize: 7
mindagap_loopnum: 100
clahe_pyramid_tile: 368
segmentation_method: "cellpose"
```

Ora il vostro comando diventa:

```bash
nextflow run ./molkart -params-file params.yaml -resume
```

Ecco fatto! Il file di parametri documenta la vostra configurazione esatta e rende facile rieseguire o condividere l'analisi.

### 1.3. Sovrascrivere i parametri

Potete comunque sovrascrivere parametri specifici dalla riga di comando:

```bash
nextflow run ./molkart -params-file params.yaml --segmentation_method "stardist" --outdir stardist_results -resume
```

La riga sopra cambia il `segmentation_method` in `stardist` e il nome di `--outdir` in `stardist_results` invece dei parametri nel file `params.yaml`.
Inoltre, potete vedere che il flag `-resume` ci ha permesso di riutilizzare i risultati di pre-elaborazione dell'esecuzione precedente, risparmiando tempo.
Potete usare questo schema per testare rapidamente diverse varianti della pipeline.

### Takeaway

I file di parametri rendono le vostre analisi riproducibili e facili da condividere.
Usateli per qualsiasi lavoro di analisi reale.

### Cosa c'è dopo?

Impariamo come le samplesheet organizzano le informazioni su più campioni.

---

## 2. Il pattern delle samplesheet

### 2.1. Cos'è una samplesheet?

Una samplesheet è un file CSV che descrive i vostri campioni di input.
Ogni riga è un campione, e le colonne specificano i file e i metadati per quel campione.

Diamo un'occhiata alla samplesheet che abbiamo usato:

```bash
cat data/samplesheet.csv
```

```csv title="samplesheet.csv"
sample,nuclear_image,spot_table,membrane_image
mem_only,https://raw.githubusercontent.com/nf-core/test-datasets/molkart/test_data/input_data/nuclear.tiff,https://raw.githubusercontent.com/nf-core/test-datasets/molkart/test_data/input_data/spots.txt,https://raw.githubusercontent.com/nf-core/test-datasets/molkart/test_data/input_data/membrane.tiff
```

Le colonne sono:

- `sample`: Identificatore univoco del campione
- `nuclear_image`: Immagine di colorazione nucleare (TIFF)
- `spot_table`: Spot di trascritti (TXT)
- `membrane_image`: Immagine di colorazione della membrana (TIFF, opzionale)

### 2.2. Percorsi dei file

Le samplesheet accettano più tipi di percorso:

- **URL**: Nextflow scarica automaticamente (come mostrato sopra)
- **Percorsi locali**: `data/nuclear.tiff` o `/absolute/path/to/nuclear.tiff`
- **Cloud storage**: `s3://bucket/nuclear.tiff`, `gs://bucket/nuclear.tiff`, `az://container/nuclear.tiff`

Potete mescolare tipi di percorso nella stessa samplesheet.

### 2.3. Creare la propria samplesheet

Prima, scarichiamo i file di dati di test localmente:

```bash
cd /workspaces/training/nf4-science/imaging/data
curl -O https://raw.githubusercontent.com/nf-core/test-datasets/molkart/test_data/input_data/nuclear.tiff
curl -O https://raw.githubusercontent.com/nf-core/test-datasets/molkart/test_data/input_data/spots.txt
curl -O https://raw.githubusercontent.com/nf-core/test-datasets/molkart/test_data/input_data/membrane.tiff
cd ..
```

Ora modifichiamo la samplesheet per fare riferimento a questi file locali:

```csv title="samplesheet.csv"
sample,nuclear_image,spot_table,membrane_image
mem_only,data/nuclear.tiff,data/spots.txt,data/membrane.tiff
```

!!! warning "Attenzione"

    Notate che i percorsi nella samplesheet sono relativi a dove **eseguite** Nextflow, non a dove si trova la samplesheet.

Infine, eseguiamo nf-core/molkart ancora una volta con la samplesheet con percorsi di file locali:

`nextflow run ./molkart -params-file params.yaml -resume`

Come potete vedere, Nextflow esegue questa esecuzione in modo simile a quando i file venivano scaricati da Github. Questa è una delle grandi caratteristiche di Nextflow: prepara i dati correttamente per voi, indipendentemente da dove si trovano.

### Takeaway

Le samplesheet organizzano dataset multi-campione in un modo che vi permette di definire esplicitamente i vostri metadati insieme ai percorsi dei file.
La maggior parte delle pipeline nf-core usa questo schema.

### Cosa c'è dopo?

Ora che abbiamo trattato gli input, esploriamo come configurare le pipeline Nextflow per diversi ambienti di calcolo.
