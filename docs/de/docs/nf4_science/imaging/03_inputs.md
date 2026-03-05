# Teil 3: Eingaben organisieren

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } KI-gestützte Übersetzung - [mehr erfahren & Verbesserungen vorschlagen](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

In Teil 2 haben wir molkart mit mehreren Parametern über die Kommandozeile ausgeführt.
Jetzt lernen wir zwei bessere Ansätze zum Verwalten von Eingaben kennen: **Parameterdateien** und **Samplesheets**.

## 1. Parameterdateien verwenden

### 1.1. Das Problem mit langen Befehlszeilen

Erinnere dich an unseren Befehl aus Teil 2:

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

Das funktioniert, ist aber schwer zu reproduzieren, zu teilen oder zu ändern.
Was, wenn du dieselbe Analyse nächsten Monat erneut ausführen musst?
Was, wenn ein Kollege genau deine Einstellungen verwenden möchte?

### 1.2. Lösung: Verwende eine Parameterdatei

Erstelle eine Datei namens `params.yaml`:

```yaml title="params.yaml"
input: "data/samplesheet.csv"
outdir: "results"
mindagap_tilesize: 90
mindagap_boxsize: 7
mindagap_loopnum: 100
clahe_pyramid_tile: 368
segmentation_method: "cellpose"
```

Jetzt wird dein Befehl zu:

```bash
nextflow run ./molkart -params-file params.yaml -resume
```

Das war's! Die Parameterdatei dokumentiert deine exakte Konfiguration und macht es einfach, sie erneut auszuführen oder zu teilen.

### 1.3. Parameter überschreiben

Du kannst immer noch spezifische Parameter über die Kommandozeile überschreiben:

```bash
nextflow run ./molkart -params-file params.yaml --segmentation_method "stardist" --outdir stardist_results -resume
```

Die obige Zeile ändert die `segmentation_method` zu `stardist` und den Namen von `--outdir` zu `stardist_results` anstelle der Parameter in der `params.yaml`-Datei.
Außerdem kannst du sehen, dass das `-resume`-Flag es uns ermöglicht hat, die Vorverarbeitungsergebnisse aus dem vorherigen Lauf wiederzuverwenden, was Zeit spart.
Du kannst dieses Muster nutzen, um schnell verschiedene Variationen der Pipeline zu testen.

### Fazit

Parameterdateien machen deine Analysen reproduzierbar und einfach zu teilen.
Verwende sie für jede echte Analysearbeit.

### Wie geht es weiter?

Erfahre, wie Samplesheets Informationen über mehrere Proben organisieren.

---

## 2. Das Samplesheet-Muster

### 2.1. Was ist ein Samplesheet?

Ein Samplesheet ist eine CSV-Datei, die deine Eingabeproben beschreibt.
Jede Zeile ist eine Probe, und die Spalten spezifizieren die Dateien und Metadaten für diese Probe.

Schauen wir uns das Samplesheet an, das wir verwendet haben:

```bash
cat data/samplesheet.csv
```

```csv title="samplesheet.csv"
sample,nuclear_image,spot_table,membrane_image
mem_only,https://raw.githubusercontent.com/nf-core/test-datasets/molkart/test_data/input_data/nuclear.tiff,https://raw.githubusercontent.com/nf-core/test-datasets/molkart/test_data/input_data/spots.txt,https://raw.githubusercontent.com/nf-core/test-datasets/molkart/test_data/input_data/membrane.tiff
```

Die Spalten sind:

- `sample`: Eindeutige Proben-Kennung
- `nuclear_image`: Kernfärbungsbild (TIFF)
- `spot_table`: Transkript-Spots (TXT)
- `membrane_image`: Membranfärbungsbild (TIFF, optional)

### 2.2. Dateipfade

Samplesheets akzeptieren mehrere Pfadtypen:

- **URLs**: Nextflow lädt automatisch herunter (wie oben gezeigt)
- **Lokale Pfade**: `data/nuclear.tiff` oder `/absolute/path/to/nuclear.tiff`
- **Cloud-Speicher**: `s3://bucket/nuclear.tiff`, `gs://bucket/nuclear.tiff`, `az://container/nuclear.tiff`

Du kannst Pfadtypen im selben Samplesheet mischen.

### 2.3. Dein eigenes Samplesheet erstellen

Lass uns zunächst die Testdaten lokal herunterladen:

```bash
cd /workspaces/training/nf4-science/imaging/data
curl -O https://raw.githubusercontent.com/nf-core/test-datasets/molkart/test_data/input_data/nuclear.tiff
curl -O https://raw.githubusercontent.com/nf-core/test-datasets/molkart/test_data/input_data/spots.txt
curl -O https://raw.githubusercontent.com/nf-core/test-datasets/molkart/test_data/input_data/membrane.tiff
cd ..
```

Jetzt ändern wir das Samplesheet, um auf diese lokalen Dateien zu verweisen:

```csv title="samplesheet.csv"
sample,nuclear_image,spot_table,membrane_image
mem_only,data/nuclear.tiff,data/spots.txt,data/membrane.tiff
```

!!! warning "Warnung"

    Beachte, dass die Pfade im Samplesheet relativ zu dem Ort sind, an dem du Nextflow **ausführst**, nicht wo sich das Samplesheet befindet.

Lass uns schließlich nf-core/molkart noch einmal mit dem Samplesheet mit lokalen Dateipfaden ausführen:

`nextflow run ./molkart -params-file params.yaml -resume`

Wie du sehen kannst, führt Nextflow diesen Lauf ähnlich aus wie bei den Dateien, die von Github heruntergeladen wurden. Dies ist eine der großartigen Funktionen von Nextflow: Es bereitet die Daten für dich korrekt vor, unabhängig davon, wo sie sich befinden.

### Fazit

Samplesheets organisieren Multi-Proben-Datensätze auf eine Weise, die es dir ermöglicht, deine Metadaten zusammen mit den Dateipfaden explizit zu definieren.
Die meisten nf-core-Pipelines verwenden dieses Muster.

### Wie geht es weiter?

Nachdem wir Eingaben behandelt haben, lass uns untersuchen, wie man Nextflow-Pipelines für verschiedene Rechenumgebungen konfiguriert.
