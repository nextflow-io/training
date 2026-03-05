# Parte 3: Organizzazione degli input

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Traduzione assistita da IA - [scopri di più e suggerisci miglioramenti](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

Nella Parte 2, abbiamo eseguito molkart con più parametri dalla riga di comando.
Ora impareremo due approcci migliori per gestire gli input: **file di parametri** e **samplesheet**.

## 1. Utilizzo dei file di parametri

### 1.1. Il problema delle righe di comando lunghe

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

Questo funziona, ma è difficile da riprodurre, condividere o modificare.
Cosa succede se è necessario eseguire nuovamente la stessa analisi il prossimo mese?
Cosa succede se un collaboratore vuole utilizzare le vostre stesse impostazioni?

### 1.2. Soluzione: Utilizzare un file di parametri

Creare un file chiamato `params.yaml`:

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

Ecco fatto! Il file di parametri documenta la vostra configurazione esatta e rende facile rieseguire o condividere.

### 1.3. Sovrascrittura dei parametri

È ancora possibile sovrascrivere parametri specifici dalla riga di comando:

```bash
nextflow run ./molkart -params-file params.yaml --segmentation_method "stardist" --outdir stardist_results -resume
```

La riga sopra cambia il `segmentation_method` in `stardist` e il nome di `--outdir` in `stardist_results` invece dei parametri nel file `params.yaml`.
Inoltre, potete vedere che il flag `-resume` ci ha permesso di riutilizzare i risultati di pre-elaborazione dall'esecuzione precedente, risparmiando tempo.
Può utilizzare questo schema per testare rapidamente diverse variazioni della pipeline.

### Takeaway

I file di parametri rendono le vostre analisi riproducibili e facili da condividere.
Li utilizzi per qualsiasi lavoro di analisi reale.

### Prossimi passi

Scopra come i samplesheet organizzano le informazioni su più campioni.

---

## 2. Il pattern samplesheet

### 2.1. Cos'è un samplesheet?

Un samplesheet è un file CSV che descrive i vostri campioni di input.
Ogni riga è un campione e le colonne specificano i file e i metadati per quel campione.

Esaminiamo il samplesheet che abbiamo utilizzato:

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
- `spot_table`: Punti di trascrizione (TXT)
- `membrane_image`: Immagine di colorazione della membrana (TIFF, opzionale)

### 2.2. Percorsi dei file

I samplesheet accettano diversi tipi di percorso:

- **URL**: Nextflow scarica automaticamente (come mostrato sopra)
- **Percorsi locali**: `data/nuclear.tiff` o `/absolute/path/to/nuclear.tiff`
- **Cloud storage**: `s3://bucket/nuclear.tiff`, `gs://bucket/nuclear.tiff`, `az://container/nuclear.tiff`

Può combinare tipi di percorso nello stesso samplesheet.

### 2.3. Creazione del proprio samplesheet

Prima, scarichiamo i file di dati di test localmente:

```bash
cd /workspaces/training/nf4-science/imaging/data
curl -O https://raw.githubusercontent.com/nf-core/test-datasets/molkart/test_data/input_data/nuclear.tiff
curl -O https://raw.githubusercontent.com/nf-core/test-datasets/molkart/test_data/input_data/spots.txt
curl -O https://raw.githubusercontent.com/nf-core/test-datasets/molkart/test_data/input_data/membrane.tiff
cd ..
```

Ora modifichiamo il samplesheet per fare riferimento a questi file locali:

```csv title="samplesheet.csv"
sample,nuclear_image,spot_table,membrane_image
mem_only,data/nuclear.tiff,data/spots.txt,data/membrane.tiff
```

!!! warning "Avviso"

    Noti che i percorsi nel samplesheet sono relativi a dove **esegue** Nextflow, non a dove si trova il samplesheet.

Infine, eseguiamo nf-core/molkart un'altra volta con il samplesheet con percorsi di file locali:

`nextflow run ./molkart -params-file params.yaml -resume`

Come potete vedere, Nextflow esegue questa esecuzione in modo simile a quando i file sono stati scaricati da GitHub. Questa è una delle grandi caratteristiche di Nextflow: prepara i dati correttamente per voi, indipendentemente da dove si trovano.

### Takeaway

I samplesheet organizzano set di dati multi-campione in modo da permettervi di definire esplicitamente i vostri metadati insieme ai percorsi dei file.
La maggior parte delle pipeline nf-core utilizza questo pattern.

### Prossimi passi

Ora che abbiamo trattato gli input, esploriamo come configurare le pipeline Nextflow per diversi ambienti di calcolo.
