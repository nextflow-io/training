# Partie 3 : Organisation des entrées

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Traduction assistée par IA - [en savoir plus et suggérer des améliorations](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

Dans la Partie 2, nous avons exécuté molkart avec plusieurs paramètres sur la ligne de commande.
Nous allons maintenant apprendre deux meilleures approches pour gérer les entrées : les **fichiers de paramètres** et les **feuilles d'échantillons**.

## 1. Utilisation des fichiers de paramètres

### 1.1. Le problème des longues lignes de commande

Rappelons notre commande de la Partie 2 :

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

Cela fonctionne, mais c'est difficile à reproduire, partager ou modifier.
Que faire si vous devez exécuter la même analyse le mois prochain ?
Que faire si un collaborateur souhaite utiliser exactement vos paramètres ?

### 1.2. Solution : Utiliser un fichier de paramètres

Créez un fichier appelé `params.yaml` :

```yaml title="params.yaml"
input: "data/samplesheet.csv"
outdir: "results"
mindagap_tilesize: 90
mindagap_boxsize: 7
mindagap_loopnum: 100
clahe_pyramid_tile: 368
segmentation_method: "cellpose"
```

Maintenant votre commande devient :

```bash
nextflow run ./molkart -params-file params.yaml -resume
```

C'est tout ! Le fichier de paramètres documente votre configuration exacte et facilite la réexécution ou le partage.

### 1.3. Remplacement des paramètres

Vous pouvez toujours remplacer des paramètres spécifiques depuis la ligne de commande :

```bash
nextflow run ./molkart -params-file params.yaml --segmentation_method "stardist" --outdir stardist_results -resume
```

La ligne ci-dessus change le `segmentation_method` en `stardist` et le nom du `--outdir` en `stardist_results` au lieu des paramètres dans le fichier `params.yaml`.
De plus, vous pouvez voir que le flag `-resume` nous a permis de réutiliser les résultats de prétraitement de l'exécution précédente, ce qui économise du temps.
Vous pouvez utiliser ce modèle pour tester rapidement différentes variations du pipeline.

### Point clé

Les fichiers de paramètres rendent vos analyses reproductibles et faciles à partager.
Utilisez-les pour tout travail d'analyse réel.

### Et ensuite ?

Apprenez comment les feuilles d'échantillons organisent les informations sur plusieurs échantillons.

---

## 2. Le modèle de feuille d'échantillons

### 2.1. Qu'est-ce qu'une feuille d'échantillons ?

Une feuille d'échantillons est un fichier CSV qui décrit vos échantillons d'entrée.
Chaque ligne est un échantillon, et les colonnes spécifient les fichiers et les métadonnées pour cet échantillon.

Regardons la feuille d'échantillons que nous avons utilisée :

```bash
cat data/samplesheet.csv
```

```csv title="samplesheet.csv"
sample,nuclear_image,spot_table,membrane_image
mem_only,https://raw.githubusercontent.com/nf-core/test-datasets/molkart/test_data/input_data/nuclear.tiff,https://raw.githubusercontent.com/nf-core/test-datasets/molkart/test_data/input_data/spots.txt,https://raw.githubusercontent.com/nf-core/test-datasets/molkart/test_data/input_data/membrane.tiff
```

Les colonnes sont :

- `sample` : Identifiant unique de l'échantillon
- `nuclear_image` : Image de coloration nucléaire (TIFF)
- `spot_table` : Points de transcription (TXT)
- `membrane_image` : Image de coloration membranaire (TIFF, optionnel)

### 2.2. Chemins de fichiers

Les feuilles d'échantillons acceptent plusieurs types de chemins :

- **URLs** : Nextflow télécharge automatiquement (comme montré ci-dessus)
- **Chemins locaux** : `data/nuclear.tiff` ou `/absolute/path/to/nuclear.tiff`
- **Stockage cloud** : `s3://bucket/nuclear.tiff`, `gs://bucket/nuclear.tiff`, `az://container/nuclear.tiff`

Vous pouvez mélanger les types de chemins dans la même feuille d'échantillons.

### 2.3. Création de votre propre feuille d'échantillons

Tout d'abord, téléchargeons les fichiers de données de test localement :

```bash
cd /workspaces/training/nf4-science/imaging/data
curl -O https://raw.githubusercontent.com/nf-core/test-datasets/molkart/test_data/input_data/nuclear.tiff
curl -O https://raw.githubusercontent.com/nf-core/test-datasets/molkart/test_data/input_data/spots.txt
curl -O https://raw.githubusercontent.com/nf-core/test-datasets/molkart/test_data/input_data/membrane.tiff
cd ..
```

Maintenant, modifions la feuille d'échantillons pour référencer ces fichiers locaux :

```csv title="samplesheet.csv"
sample,nuclear_image,spot_table,membrane_image
mem_only,data/nuclear.tiff,data/spots.txt,data/membrane.tiff
```

!!! warning "Avertissement"

    Notez que les chemins dans la feuille d'échantillons sont relatifs à l'endroit où vous **exécutez** Nextflow, pas à l'endroit où se trouve la feuille d'échantillons.

Enfin, exécutons nf-core/molkart une fois de plus avec la feuille d'échantillons contenant des chemins de fichiers locaux :

`nextflow run ./molkart -params-file params.yaml -resume`

Comme vous pouvez le voir, Nextflow exécute cette exécution de manière similaire à lorsque les fichiers étaient téléchargés depuis Github. C'est l'une des grandes fonctionnalités de Nextflow : il organise les données correctement pour vous, quel que soit l'endroit où elles se trouvent.

### Point clé

Les feuilles d'échantillons organisent les ensembles de données multi-échantillons d'une manière qui vous permet de définir explicitement vos métadonnées avec les chemins de fichiers.
La plupart des pipelines nf-core utilisent ce modèle.

### Et ensuite ?

Maintenant que nous avons couvert les entrées, explorons comment configurer les pipelines Nextflow pour différents environnements informatiques.
