---
title: Modèles de script essentiels
description: Apprendre des techniques de programmation avancées dans Nextflow.
weight: 1200
---

# Modèles de script essentiels

Cette leçon vous guidera à travers les modèles de programmation essentiels qui sont fondamentaux pour créer des workflows Nextflow efficaces. Nous couvrirons des modèles pour gérer les entrées de données, transformer les valeurs, contrôler la logique de workflow, allouer dynamiquement les ressources, et plus encore.

## Objectifs d'apprentissage

- Comprendre les différences entre les paradigmes de flux de données et de scripting
- Appliquer les closures, les opérateurs ternaires et autres techniques Groovy
- Maîtriser les techniques de manipulation des métadonnées et d'extraction d'informations à partir de fichiers
- Utiliser les expressions régulières et le traitement de chaînes pour analyser les noms de fichiers
- Créer des fonctions réutilisables pour une logique complexe
- Implémenter une allocation dynamique des ressources et des stratégies de retry
- Ajouter une logique conditionnelle pour contrôler l'exécution du workflow
- Écrire un code robuste en utilisant les opérateurs de navigation sûre et Elvis
- Valider les entrées avec des messages d'erreur clairs
- Utiliser des gestionnaires d'événements pour gérer la fin du workflow

## Prérequis

- Compréhension de base des workflows Nextflow
- Familiarité avec la syntaxe DSL2
- Connaissance de base des channels et des processus
- Environnement de développement configuré avec Nextflow v23.04.0 ou supérieur
- Accès à Docker (ou Conda) pour les conteneurs logiciels

## Pour commencer

Ce tutoriel suppose que vous avez une connaissance de base de la syntaxe Nextflow et des opérations sur les channels.

Pour commencer, naviguez vers le dossier `side-quests/essential_scripting_patterns` :

```bash
cd side-quests/essential_scripting_patterns
```

Le dépôt contient plusieurs fichiers :

- `main.nf` : Le workflow principal
- `modules/fastp.nf` : Module pour le traitement de la qualité avec FASTP
- `modules/generate_report.nf` : Module pour générer des rapports
- `modules/trimgalore.nf` : Module pour un programme alternatif de trimming
- `data/samples.csv` : Un CSV avec des informations sur les échantillons
- `data/sequences/*.fastq` : Fichiers fastq d'exemple

Examinons les fichiers :

```bash
cat main.nf
```

Vous devriez voir un workflow minimal qui utilise des inclusions de modules :

```groovy title="main.nf" linenums="1"
include { FASTP } from './modules/fastp.nf'
include { GENERATE_REPORT } from './modules/generate_report.nf'

workflow {
    ch_samples = channel.fromPath("./data/samples.csv")
        .splitCsv(header: true)
        .map { row ->
            tuple(
                [id: row.sample_id],
                file(row.file_path)
            )
        }

    ch_fastp = FASTP(ch_samples)
    GENERATE_REPORT(ch_samples)
}
```

Les modules sont également simples :

```bash
cat modules/fastp.nf
```

```groovy title="modules/fastp.nf" linenums="1"
process FASTP {
    container 'community.wave.seqera.io/library/fastp:0.24.0--62c97b06e8447690'

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("${meta.id}_trimmed.fastq.gz"), emit: fastq
    tuple val(meta), path("${meta.id}.fastp.json"), emit: json

    script:
    """
    fastp --in1 ${reads} --out1 ${meta.id}_trimmed.fastq.gz --json ${meta.id}.fastp.json --html ${meta.id}.fastp.html
    """
}
```

```bash
cat modules/generate_report.nf
```

```groovy title="modules/generate_report.nf" linenums="1"
process GENERATE_REPORT {
    container 'community.wave.seqera.io/library/bash:5.0'

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("${meta.id}_report.txt"), emit: report

    script:
    """
    echo "Sample ID: ${meta.id}" > ${meta.id}_report.txt
    echo "File processed: ${reads}" >> ${meta.id}_report.txt
    echo "Process command: fastp --in1 ${reads} --out1 ${meta.id}_trimmed.fastq.gz" >> ${meta.id}_report.txt
    """
}
```

Et enfin, examinons les données d'échantillon :

```bash
cat data/samples.csv
```

```csv title="data/samples.csv"
sample_id,organism,tissue_type,sequencing_depth,quality_score,file_path
SAMPLE_001,human,liver,30000000,38.5,./data/sequences/SAMPLE_001_S1_L001_R1_001.fastq
SAMPLE_002,mouse,brain,25000000,35.2,./data/sequences/SAMPLE_002_S2_L001_R1_001.fastq
SAMPLE_003,human,kidney,45000000,42.1,./data/sequences/SAMPLE_003_S3_L001_R1_001.fastq
```

Exécutons d'abord le workflow de base pour nous assurer qu'il fonctionne :

```bash
nextflow run main.nf
```

Vous devriez voir le workflow s'exécuter et se terminer avec succès.

Maintenant, commençons à améliorer le workflow avec des modèles de script avancés !

---

## 1. Comprendre le flux de données vs le scripting

Nextflow utilise **deux paradigmes de programmation distincts** :

1. **Flux de données** : L'orchestration des channels via des opérateurs (`.map`, `.filter`, `.branch`)
2. **Scripting** : Code Groovy exécuté dans des closures ou des blocs de script dans des processus

Les deux sont cruciaux, mais fonctionnent différemment, donc cette première section clarifiera les différences.

### 1.1. Closures et opérateurs ternaires

Un concept fondamental dans Nextflow est la **closure** - un bloc de code qui peut être passé comme un objet et exécuté plus tard. Les closures sont essentielles pour les opérateurs de channel (`map`, `filter`, etc.).

Prenons notre opérateur `.map` actuel et améliorons-le pour extraire plus de métadonnées :

```groovy
.map { row ->
    tuple(
        [id: row.sample_id],
        file(row.file_path)
    )
}
```

Le bloc `{ row -> ... }` est une closure qui reçoit un argument `row`.

Maintenant, améliorons-le pour extraire plus de métadonnées du CSV :

=== "After"

    ```groovy title="main.nf" linenums="7" hl_lines="2-8"
        .map { row ->
            def sample_meta = [
                id: row.sample_id.toLowerCase(),
                organism: row.organism,
                tissue: row.tissue_type.replaceAll('_', ' ').toLowerCase(),
                depth: row.sequencing_depth.toInteger(),
                quality: row.quality_score.toDouble()
            ]
            tuple(sample_meta, file(row.file_path))
        }
    ```

=== "Before"

    ```groovy title="main.nf" linenums="7" hl_lines="2-5"
        .map { row ->
            tuple(
                [id: row.sample_id],
                file(row.file_path)
            )
        }
    ```

Maintenant, modifions `main.nf` pour utiliser ces métadonnées supplémentaires. D'abord, changeons la façon dont nous générons le rapport. Nous allons étiqueter les échantillons de haute qualité comme étant de "haute priorité" en utilisant un **opérateur ternaire**.

=== "After"

    ```groovy title="main.nf" linenums="7" hl_lines="9-10"
        .map { row ->
            def sample_meta = [
                id: row.sample_id.toLowerCase(),
                organism: row.organism,
                tissue: row.tissue_type.replaceAll('_', ' ').toLowerCase(),
                depth: row.sequencing_depth.toInteger(),
                quality: row.quality_score.toDouble()
            ]
            def priority = sample_meta.quality > 40 ? 'high' : 'normal'
            tuple(sample_meta + [priority: priority], file(row.file_path))
        }
    ```

=== "Before"

    ```groovy title="main.nf" linenums="7" hl_lines="9"
        .map { row ->
            def sample_meta = [
                id: row.sample_id.toLowerCase(),
                organism: row.organism,
                tissue: row.tissue_type.replaceAll('_', ' ').toLowerCase(),
                depth: row.sequencing_depth.toInteger(),
                quality: row.quality_score.toDouble()
            ]
            tuple(sample_meta, file(row.file_path))
        }
    ```

L'opérateur ternaire `condition ? valeur_si_vrai : valeur_si_faux` est une façon concise d'écrire une déclaration if/else. Si `sample_meta.quality > 40`, alors `priority` sera `'high'` ; sinon, ce sera `'normal'`.

Modifions également le processus `GENERATE_REPORT` pour inclure les métadonnées supplémentaires :

=== "After"

    ```groovy title="modules/generate_report.nf" linenums="11" hl_lines="2-7"
    script:
    """
    echo "Sample ID: ${meta.id}" > ${meta.id}_report.txt
    echo "Organism: ${meta.organism}" >> ${meta.id}_report.txt
    echo "Tissue: ${meta.tissue}" >> ${meta.id}_report.txt
    echo "Sequencing depth: ${meta.depth}" >> ${meta.id}_report.txt
    echo "Quality score: ${meta.quality}" >> ${meta.id}_report.txt
    echo "Priority: ${meta.priority}" >> ${meta.id}_report.txt
    echo "File processed: ${reads}" >> ${meta.id}_report.txt
    echo "Process command: fastp --in1 ${reads} --out1 ${meta.id}_trimmed.fastq.gz" >> ${meta.id}_report.txt
    """
    ```

=== "Before"

    ```groovy title="modules/generate_report.nf" linenums="11" hl_lines="2-4"
    script:
    """
    echo "Sample ID: ${meta.id}" > ${meta.id}_report.txt
    echo "File processed: ${reads}" >> ${meta.id}_report.txt
    echo "Process command: fastp --in1 ${reads} --out1 ${meta.id}_trimmed.fastq.gz" >> ${meta.id}_report.txt
    """
    ```

Maintenant, exécutons le workflow mis à jour :

```bash
nextflow run main.nf
```

??? success "Sortie de la commande"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `main.nf` [jolly_fourier] DSL2 - revision: d3e76a7fce

    executor >  local (6)
    [a3/6b9e80] process > FASTP (sample_003)           [100%] 3 of 3 ✔
    [29/54f4b6] process > GENERATE_REPORT (sample_001) [100%] 3 of 3 ✔
    ```

Nous pouvons vérifier le résultat en examinant l'un des fichiers de rapport :

```console
cat work/29/54f4b6b0eb90fed9e3a673b8e47629/sample_001_report.txt
```

Vous devriez maintenant voir toutes les métadonnées incluses :

```
Sample ID: sample_001
Organism: human
Tissue: liver
Sequencing depth: 30000000
Quality score: 38.5
Priority: normal
File processed: SAMPLE_001_S1_L001_R1_001.fastq
Process command: fastp --in1 SAMPLE_001_S1_L001_R1_001.fastq --out1 sample_001_trimmed.fastq.gz
```

### 1.2. La collection 'meta' vs 'reads'

Dans Nextflow, les channels, tuples et collections (comme les maps et les listes) sont des structures de données fondamentales. Il existe une différence importante entre les **collections dans le flux de données** (channels/opérateurs) et les **collections dans les blocs de script**.

Pour démontrer cela, modifions notre workflow pour extraire des métadonnées des noms de fichiers FASTQ en utilisant des expressions régulières. Il est courant que les noms de fichiers FASTQ suivent une convention comme : `SAMPLE_001_S1_L001_R1_001.fastq`.

Dans ce format :

- `SAMPLE_001` : ID de l'échantillon
- `S1` : Numéro d'échantillon
- `L001` : Numéro de lane
- `R1` : Numéro de read
- `001` : Fragment ou chunk

Extrayons ces métadonnées du nom de fichier FASTQ :

=== "After"

    ```groovy title="main.nf" linenums="7" hl_lines="12-20"
        .map { row ->
            def sample_meta = [
                id: row.sample_id.toLowerCase(),
                organism: row.organism,
                tissue: row.tissue_type.replaceAll('_', ' ').toLowerCase(),
                depth: row.sequencing_depth.toInteger(),
                quality: row.quality_score.toDouble()
            ]
            def priority = sample_meta.quality > 40 ? 'high' : 'normal'
            def fastq_path = file(row.file_path)

            // Extract metadata from filename using regex
            def m = (fastq_path.name =~ /^(.+)_S(\d+)_L(\d{3})_(R[12])_(\d{3})\.fastq(?:\.gz)?$/)
            def file_meta = m ? [
                sample_num: m[0][2].toInteger(),
                lane: m[0][3],
                read: m[0][4],
                chunk: m[0][5]
            ] : [:]

            tuple(sample_meta + file_meta + [priority: priority], fastq_path)
        }
    ```

=== "Before"

    ```groovy title="main.nf" linenums="7" hl_lines="10"
        .map { row ->
            def sample_meta = [
                id: row.sample_id.toLowerCase(),
                organism: row.organism,
                tissue: row.tissue_type.replaceAll('_', ' ').toLowerCase(),
                depth: row.sequencing_depth.toInteger(),
                quality: row.quality_score.toDouble()
            ]
            def priority = sample_meta.quality > 40 ? 'high' : 'normal'
            tuple(sample_meta + [priority: priority], file(row.file_path))
        }
    ```

Il y a beaucoup à déballer ici !

1. `fastq_path.name` obtient le nom du fichier (sans le chemin)
2. L'opérateur `=~` est pour la correspondance de modèles regex
3. Le modèle `/^(.+)_S(\d+)_L(\d{3})_(R[12])_(\d{3})\.fastq(?:\.gz)?$/` capture les composants
4. `m ? ... : ...` est un opérateur ternaire qui gère le cas où le nom de fichier ne correspond pas au modèle
5. `m[0][2]` accède au deuxième groupe de capture (les indices commencent à 1 pour les groupes)
6. La fonction `toInteger()` convertit la chaîne capturée en nombre entier
7. `sample_meta + file_meta + [priority: priority]` fusionne les trois maps en un seul

Modifions également le processus `GENERATE_REPORT` pour inclure ces nouvelles métadonnées :

=== "After"

    ```groovy title="modules/generate_report.nf" linenums="11" hl_lines="7-10"
    script:
    """
    echo "Sample ID: ${meta.id}" > ${meta.id}_report.txt
    echo "Organism: ${meta.organism}" >> ${meta.id}_report.txt
    echo "Tissue: ${meta.tissue}" >> ${meta.id}_report.txt
    echo "Sequencing depth: ${meta.depth}" >> ${meta.id}_report.txt
    echo "Quality score: ${meta.quality}" >> ${meta.id}_report.txt
    echo "Sample number: ${meta.sample_num}" >> ${meta.id}_report.txt
    echo "Lane: ${meta.lane}" >> ${meta.id}_report.txt
    echo "Read: ${meta.read}" >> ${meta.id}_report.txt
    echo "Chunk: ${meta.chunk}" >> ${meta.id}_report.txt
    echo "Priority: ${meta.priority}" >> ${meta.id}_report.txt
    echo "File processed: ${reads}" >> ${meta.id}_report.txt
    echo "Process command: fastp --in1 ${reads} --out1 ${meta.id}_trimmed.fastq.gz" >> ${meta.id}_report.txt
    """
    ```

=== "Before"

    ```groovy title="modules/generate_report.nf" linenums="11" hl_lines="7"
    script:
    """
    echo "Sample ID: ${meta.id}" > ${meta.id}_report.txt
    echo "Organism: ${meta.organism}" >> ${meta.id}_report.txt
    echo "Tissue: ${meta.tissue}" >> ${meta.id}_report.txt
    echo "Sequencing depth: ${meta.depth}" >> ${meta.id}_report.txt
    echo "Quality score: ${meta.quality}" >> ${meta.id}_report.txt
    echo "Priority: ${meta.priority}" >> ${meta.id}_report.txt
    echo "File processed: ${reads}" >> ${meta.id}_report.txt
    echo "Process command: fastp --in1 ${reads} --out1 ${meta.id}_trimmed.fastq.gz" >> ${meta.id}_report.txt
    """
    ```

Exécutons le workflow à nouveau :

```bash
nextflow run main.nf
```

??? success "Sortie de la commande"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `main.nf` [jolly_hawking] DSL2 - revision: cd0a5b0d29

    executor >  local (6)
    [b3/1cb89f] process > FASTP (sample_001)           [100%] 3 of 3 ✔
    [e6/c2f254] process > GENERATE_REPORT (sample_001) [100%] 3 of 3 ✔
    ```

Vérifions l'un des fichiers de rapport mis à jour :

```console
cat work/e6/c2f2542ec23ee6e7f0aa9c66c12e30/sample_001_report.txt
```

Vous devriez maintenant voir les métadonnées supplémentaires extraites du nom de fichier :

```
Sample ID: sample_001
Organism: human
Tissue: liver
Sequencing depth: 30000000
Quality score: 38.5
Sample number: 1
Lane: 001
Read: R1
Chunk: 001
Priority: normal
File processed: SAMPLE_001_S1_L001_R1_001.fastq
Process command: fastp --in1 SAMPLE_001_S1_L001_R1_001.fastq --out1 sample_001_trimmed.fastq.gz
```

### 1.3. Opérations sur les collections dans les closures vs manipulation des channels

Un point souvent source de confusion dans Nextflow concerne les opérations sur les collections à l'intérieur des closures. Les méthodes comme `collect()` fonctionnent différemment selon le contexte :

1. Dans les **channels Nextflow** : `channel.collect()` est un opérateur qui rassemble tous les éléments d'un channel dans une seule liste
2. Dans les **listes et maps Groovy** : `list.collect {...}` applique une fonction à chaque élément et retourne une nouvelle liste

Modifions le processus FASTP pour illustrer cette différence :

```groovy title="modules/fastp.nf" linenums="11" hl_lines="3-5"
script:
"""
# Démonstration d'opérations de collection dans un script
def options = ["--in1", reads, "--out1", "${meta.id}_trimmed.fastq.gz", "--json", "${meta.id}.fastp.json", "--html", "${meta.id}.fastp.html"]
fastp ${options.join(' ')}
"""
```

Cet exemple utilise `options.join(' ')` pour joindre les éléments de la liste en une seule chaîne avec des espaces entre eux.

Cependant, cela générera une erreur car nous essayons d'exécuter du code Groovy à l'intérieur d'un bloc de script bash. Changeons-le pour déplacer la logique de collection en dehors du bloc de script :

=== "After"

    ```groovy title="modules/fastp.nf" linenums="11" hl_lines="1-3"
    def options = ["--in1", reads, "--out1", "${meta.id}_trimmed.fastq.gz", "--json", "${meta.id}.fastp.json", "--html", "${meta.id}.fastp.html"]
    def cmd = "fastp ${options.join(' ')}"

    script:
    """
    $cmd
    """
    ```

=== "Before"

    ```groovy title="modules/fastp.nf" linenums="11" hl_lines="2-3"
    script:
    """
    fastp --in1 ${reads} --out1 ${meta.id}_trimmed.fastq.gz --json ${meta.id}.fastp.json --html ${meta.id}.fastp.html
    """
    ```

Exécutez le workflow à nouveau et cela devrait fonctionner ! Cela démontre l'utilisation du scripting Groovy pour la manipulation de collections avant de passer la commande au bloc de script.

```bash
nextflow run main.nf
```

??? success "Sortie de la commande"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `main.nf` [compassionate_shaw] DSL2 - revision: 3471dc57d9

    executor >  local (6)
    [32/8e94af] process > FASTP (sample_001)           [100%] 3 of 3 ✔
    [6e/e3e56d] process > GENERATE_REPORT (sample_001) [100%] 3 of 3 ✔
    ```

### Takeaway

Dans cette section, nous avons appris les différences entre le **flux de données** (opérations de channel) et le **scripting** (code à l'intérieur de closures et de blocs de script). Nous avons utilisé :

- Des **closures** pour extraire et transformer les métadonnées à partir des données d'entrée
- Des **opérateurs ternaires** (`condition ? valeur_si_vrai : valeur_si_faux`) pour une logique conditionnelle compacte
- Des **expressions régulières** avec l'opérateur `=~` pour extraire des composants à partir des noms de fichiers
- Des **manipulations de collections** comme `join()` pour créer des chaînes de commandes

Ces techniques sont fondamentales pour écrire des workflows Nextflow propres, efficaces et faciles à maintenir.

---

## 2. Traitement des chaînes pour gérer les noms de fichiers et les métadonnées

Le traitement des chaînes est une tâche courante dans les workflows bioinformatiques, surtout lors de l'extraction de métadonnées à partir des noms de fichiers ou de la génération dynamique de scripts. Nextflow dispose de puissantes capacités de manipulation de chaînes qui sont cruciales pour des workflows robustes.

### 2.1. Expressions régulières pour analyser les noms de fichiers

Nous avons déjà utilisé des expressions régulières pour extraire des métadonnées à partir des noms de fichiers FASTQ. Comprenons mieux comment fonctionne l'opérateur de correspondance de motifs `=~`.

Lorsque vous écrivez `x =~ /pattern/` :

1. Cela crée un objet `java.util.regex.Matcher`
2. S'il est évalué dans un contexte booléen, il vérifie s'il y a **une correspondance quelconque**
3. S'il est assigné à une variable, vous pouvez accéder aux **groupes de capture**

Voyons comment nous pourrions compliquer notre regex fastq pour gérer plus de variations :

=== "After"

    ```groovy title="main.nf" linenums="14" hl_lines="2"
            // Extract metadata from filename using regex
            def m = (fastq_path.name =~ /^(.+)_S(\d+)_L(\d{3})_(R[12])_(\d{3})\.fastq(?:\.gz)?$/)
            def file_meta = m ? [
                sample_num: m[0][2].toInteger(),
                lane: m[0][3],
                read: m[0][4],
                chunk: m[0][5]
            ] : [:]
    ```

=== "Before"

    ```groovy title="main.nf" linenums="14" hl_lines="2"
            // Extract metadata from filename using regex
            def m = (fastq_path.name =~ /^(.+)_S(\d+)_L(\d{3})_(R[12])_(\d{3})\.fastq$/)
            def file_meta = m ? [
                sample_num: m[0][2].toInteger(),
                lane: m[0][3],
                read: m[0][4],
                chunk: m[0][5]
            ] : [:]
    ```

Maintenant, `/^(.+)_S(\d+)_L(\d{3})_(R[12])_(\d{3})\.fastq(?:\.gz)?$/` correspondra à la fois à `.fastq` et à `.fastq.gz`.

!!! note "Groupes de capture dans les regex" - Les parenthèses `(...)` définissent des "groupes de capture" - `m[0]` est la correspondance complète - `m[0][1]`, `m[0][2]`, etc. sont les groupes de capture (en commençant par 1)

Si vous aviez des noms de fichiers avec des conventions différentes, vous pourriez utiliser l'opérateur OR (`|`) dans votre regex :

```groovy title="exemple de regex"
def m = (filename =~ /^(.+)_S(\d+)_L(\d{3})_(R[12])_(\d{3})\.fastq(?:\.gz)?$/ |
                      /^(.+)_(\d{6})_([ACGT]+)_L(\d{3})_(R[12])\.fastq(?:\.gz)?$/)
```

### 2.2. Génération dynamique de scripts

Une autre application puissante du traitement des chaînes est la génération dynamique de scripts bash basés sur les métadonnées ou les entrées. C'est particulièrement utile pour la logique conditionnelle à l'intérieur des processus.

Modifions le processus `GENERATE_REPORT` pour générer différents types de rapports en fonction de la priorité :

=== "After"

    ```groovy title="modules/generate_report.nf" linenums="10" hl_lines="3-11"
    script:
    // Generate different report content based on priority
    def report_content = meta.priority == 'high' ?
    """
    echo "ÉCHANTILLON HAUTE PRIORITÉ" > ${meta.id}_report.txt
    echo "===============================================" >> ${meta.id}_report.txt
    """ :
    """
    echo "Échantillon standard" > ${meta.id}_report.txt
    echo "---------------------------------------------" >> ${meta.id}_report.txt
    """

    """
    ${report_content}
    echo "Sample ID: ${meta.id}" >> ${meta.id}_report.txt
    echo "Organism: ${meta.organism}" >> ${meta.id}_report.txt
    echo "Tissue: ${meta.tissue}" >> ${meta.id}_report.txt
    echo "Sequencing depth: ${meta.depth}" >> ${meta.id}_report.txt
    echo "Quality score: ${meta.quality}" >> ${meta.id}_report.txt
    echo "Sample number: ${meta.sample_num}" >> ${meta.id}_report.txt
    echo "Lane: ${meta.lane}" >> ${meta.id}_report.txt
    echo "Read: ${meta.read}" >> ${meta.id}_report.txt
    echo "Chunk: ${meta.chunk}" >> ${meta.id}_report.txt
    echo "Priority: ${meta.priority}" >> ${meta.id}_report.txt
    echo "File processed: ${reads}" >> ${meta.id}_report.txt
    echo "Process command: fastp --in1 ${reads} --out1 ${meta.id}_trimmed.fastq.gz" >> ${meta.id}_report.txt
    """
    ```

=== "Before"

    ```groovy title="modules/generate_report.nf" linenums="11" hl_lines="2-13"
    script:
    """
    echo "Sample ID: ${meta.id}" > ${meta.id}_report.txt
    echo "Organism: ${meta.organism}" >> ${meta.id}_report.txt
    echo "Tissue: ${meta.tissue}" >> ${meta.id}_report.txt
    echo "Sequencing depth: ${meta.depth}" >> ${meta.id}_report.txt
    echo "Quality score: ${meta.quality}" >> ${meta.id}_report.txt
    echo "Sample number: ${meta.sample_num}" >> ${meta.id}_report.txt
    echo "Lane: ${meta.lane}" >> ${meta.id}_report.txt
    echo "Read: ${meta.read}" >> ${meta.id}_report.txt
    echo "Chunk: ${meta.chunk}" >> ${meta.id}_report.txt
    echo "Priority: ${meta.priority}" >> ${meta.id}_report.txt
    echo "File processed: ${reads}" >> ${meta.id}_report.txt
    echo "Process command: fastp --in1 ${reads} --out1 ${meta.id}_trimmed.fastq.gz" >> ${meta.id}_report.txt
    """
    ```

Exécutons le workflow pour voir le résultat :

```bash
nextflow run main.nf
```

??? success "Sortie de la commande"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `main.nf` [condescending_northcutt] DSL2 - revision: 1a3c16a96f

    executor >  local (6)
    [95/e633a2] process > FASTP (sample_001)           [100%] 3 of 3 ✔
    [a8/e4c214] process > GENERATE_REPORT (sample_001) [100%] 3 of 3 ✔
    ```

Examinons les fichiers de rapport générés - cherchons spécifiquement celui de haute priorité :

```console
find work -name "sample_003_report.txt" -exec cat {} \;
```

Vous devriez voir l'en-tête spécial :

```
ÉCHANTILLON HAUTE PRIORITÉ
===============================================
Sample ID: sample_003
Organism: human
...
```

### 2.3. Interpolation de variables : quand Nextflow évalue vs quand bash évalue

Un point subtil mais crucial à comprendre est de savoir quand l'interpolation de variables se produit :

1. `${var}` - Interpolé par Nextflow pendant la compilation du script
2. `\${var}` - Échappé, passé littéralement à bash en tant que `${var}` (pour les variables d'environnement bash)

Nous pouvons voir cela en action en mettant à jour le processus `GENERATE_REPORT` pour utiliser une variable d'environnement :

=== "After"

    ```groovy title="modules/generate_report.nf" linenums="11" hl_lines="15"
    script:
    // Generate different report content based on priority
    def report_content = meta.priority == 'high' ?
    """
    echo "ÉCHANTILLON HAUTE PRIORITÉ" > ${meta.id}_report.txt
    echo "===============================================" >> ${meta.id}_report.txt
    """ :
    """
    echo "Échantillon standard" > ${meta.id}_report.txt
    echo "---------------------------------------------" >> ${meta.id}_report.txt
    """

    """
    ${report_content}
    echo "Sample ID: ${meta.id}" >> ${meta.id}_report.txt
    echo "Processed by: \${USER}" >> ${meta.id}_report.txt
    echo "Organism: ${meta.organism}" >> ${meta.id}_report.txt
    echo "Tissue: ${meta.tissue}" >> ${meta.id}_report.txt
    echo "Sequencing depth: ${meta.depth}" >> ${meta.id}_report.txt
    echo "Quality score: ${meta.quality}" >> ${meta.id}_report.txt
    echo "Sample number: ${meta.sample_num}" >> ${meta.id}_report.txt
    echo "Lane: ${meta.lane}" >> ${meta.id}_report.txt
    echo "Read: ${meta.read}" >> ${meta.id}_report.txt
    echo "Chunk: ${meta.chunk}" >> ${meta.id}_report.txt
    echo "Priority: ${meta.priority}" >> ${meta.id}_report.txt
    echo "File processed: ${reads}" >> ${meta.id}_report.txt
    echo "Process command: fastp --in1 ${reads} --out1 ${meta.id}_trimmed.fastq.gz" >> ${meta.id}_report.txt
    """
    ```

=== "Before"

    ```groovy title="modules/generate_report.nf" linenums="11"
    script:
    // Generate different report content based on priority
    def report_content = meta.priority == 'high' ?
    """
    echo "ÉCHANTILLON HAUTE PRIORITÉ" > ${meta.id}_report.txt
    echo "===============================================" >> ${meta.id}_report.txt
    """ :
    """
    echo "Échantillon standard" > ${meta.id}_report.txt
    echo "---------------------------------------------" >> ${meta.id}_report.txt
    """

    """
    ${report_content}
    echo "Sample ID: ${meta.id}" >> ${meta.id}_report.txt
    echo "Organism: ${meta.organism}" >> ${meta.id}_report.txt
    echo "Tissue: ${meta.tissue}" >> ${meta.id}_report.txt
    echo "Sequencing depth: ${meta.depth}" >> ${meta.id}_report.txt
    echo "Quality score: ${meta.quality}" >> ${meta.id}_report.txt
    echo "Sample number: ${meta.sample_num}" >> ${meta.id}_report.txt
    echo "Lane: ${meta.lane}" >> ${meta.id}_report.txt
    echo "Read: ${meta.read}" >> ${meta.id}_report.txt
    echo "Chunk: ${meta.chunk}" >> ${meta.id}_report.txt
    echo "Priority: ${meta.priority}" >> ${meta.id}_report.txt
    echo "File processed: ${reads}" >> ${meta.id}_report.txt
    echo "Process command: fastp --in1 ${reads} --out1 ${meta.id}_trimmed.fastq.gz" >> ${meta.id}_report.txt
    """
    ```

Si nous exécutons le workflow tel quel, il échouera :

```bash
nextflow run main.nf
```

??? failure "Sortie de la commande"

    ```console
    Error modules/generate_report.nf:15:27: `USER` is not defined
    │  15 |     echo "Processed by: ${USER}" >> ${meta.id}_report.txt
    ╰     |                           ^^^^

    ERROR ~ Script compilation failed
    ```

Le problème est que Nextflow essaie d'interpréter `${USER}` comme une variable Nextflow, mais il n'y a pas de variable nommée `USER`. Nous devons échapper le `$` pour qu'il soit transmis à bash, qui a bien une variable d'environnement `USER` :

```groovy title="modules/generate_report.nf" linenums="15" hl_lines="1"
    echo "Processed by: \${USER}" >> ${meta.id}_report.txt
```

Ça marche maintenant ! La barre oblique inverse (`\`) dit à Nextflow "n'interprète pas ceci, passe-le à Bash."

### Takeaway

Dans cette section, vous avez appris des techniques de **traitement de chaînes** :

- **Expressions régulières pour l'analyse de fichiers** : Utilisation de l'opérateur `=~` et des modèles regex (`~/pattern/`) pour extraire des métadonnées à partir de conventions de nommage de fichiers complexes
- **Génération dynamique de scripts** : Utilisation de la logique conditionnelle (if/else, opérateurs ternaires) pour générer différentes chaînes de script basées sur les caractéristiques d'entrée
- **Interpolation de variables** : Comprendre quand Nextflow interprète les chaînes vs quand le shell le fait
  - `${var}` - Variables Nextflow (interpolées par Nextflow au moment de la compilation du workflow)
  - `\${var}` - Variables d'environnement shell (échappées, transmises à bash au moment de l'exécution)
  - `\$(cmd)` - Substitution de commandes shell (échappée, exécutée par bash au moment de l'exécution)

Ces modèles de traitement et de génération de chaînes sont essentiels pour gérer les divers formats de fichiers et conventions de nommage que vous rencontrerez dans les workflows bioinformatiques du monde réel.

---

## 3. Création de fonctions réutilisables

La logique de workflow complexe intégrée dans les opérateurs de canaux ou les définitions de processus réduit la lisibilité et la maintenabilité. Les **fonctions** vous permettent d'extraire cette logique dans des composants nommés et réutilisables.

Notre opération map est devenue longue et complexe. Extrayons-la dans une fonction réutilisable en utilisant le mot-clé `def`.

Pour illustrer à quoi cela ressemble avec notre workflow existant, faites la modification ci-dessous, en utilisant `def` pour définir une fonction réutilisable appelée `separateMetadata` :

=== "After"

    ```groovy title="main.nf" linenums="1" hl_lines="4-24 29"
    include { FASTP } from './modules/fastp.nf'
    include { GENERATE_REPORT } from './modules/generate_report.nf'

    def separateMetadata(row) {
        def sample_meta = [
            id: row.sample_id.toLowerCase(),
            organism: row.organism,
            tissue: row.tissue_type.replaceAll('_', ' ').toLowerCase(),
            depth: row.sequencing_depth.toInteger(),
            quality: row.quality_score.toDouble()
        ]
        def fastq_path = file(row.file_path)

        def m = (fastq_path.name =~ /^(.+)_S(\d+)_L(\d{3})_(R[12])_(\d{3})\.fastq(?:\.gz)?$/)
        def file_meta = m ? [
            sample_num: m[0][2].toInteger(),
            lane: m[0][3],
            read: m[0][4],
            chunk: m[0][5]
        ] : [:]

        def priority = sample_meta.quality > 40 ? 'high' : 'normal'
        return tuple(sample_meta + file_meta + [priority: priority], fastq_path)
    }

    workflow {
        ch_samples = channel.fromPath("./data/samples.csv")
            .splitCsv(header: true)
            .map{ row -> separateMetadata(row) }

        ch_fastp = FASTP(ch_samples)
        GENERATE_REPORT(ch_samples)
    }
    ```

=== "Before"

    ```groovy title="main.nf" linenums="1" hl_lines="7-27"
    include { FASTP } from './modules/fastp.nf'
    include { GENERATE_REPORT } from './modules/generate_report.nf'

    workflow {
        ch_samples = channel.fromPath("./data/samples.csv")
            .splitCsv(header: true)
            .map { row ->
                def sample_meta = [
                    id: row.sample_id.toLowerCase(),
                    organism: row.organism,
                    tissue: row.tissue_type.replaceAll('_', ' ').toLowerCase(),
                    depth: row.sequencing_depth.toInteger(),
                    quality: row.quality_score.toDouble()
                ]
                def fastq_path = file(row.file_path)

                def m = (fastq_path.name =~ /^(.+)_S(\d+)_L(\d{3})_(R[12])_(\d{3})\.fastq(?:\.gz)?$/)
                def file_meta = m ? [
                    sample_num: m[0][2].toInteger(),
                    lane: m[0][3],
                    read: m[0][4],
                    chunk: m[0][5]
                ] : [:]

                def priority = sample_meta.quality > 40 ? 'high' : 'normal'
                return tuple(sample_meta + file_meta + [priority: priority], fastq_path)
            }

        ch_fastp = FASTP(ch_samples)
        GENERATE_REPORT(ch_samples)
    }
    ```

En extrayant cette logique dans une fonction, nous avons réduit la logique réelle du workflow à quelque chose de beaucoup plus propre :

```groovy title="workflow minimal"
    ch_samples = channel.fromPath("./data/samples.csv")
        .splitCsv(header: true)
        .map{ row -> separateMetadata(row) }

    ch_fastp = FASTP(ch_samples)
    GENERATE_REPORT(ch_samples)
```

Cela rend la logique du workflow beaucoup plus facile à lire et à comprendre d'un coup d'œil. La fonction `separateMetadata` encapsule toute la logique complexe pour analyser et enrichir les métadonnées, la rendant réutilisable et testable.

Exécutez le workflow pour vous assurer qu'il fonctionne toujours :

```bash
nextflow run main.nf
```

??? success "Sortie de la commande"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `main.nf` [admiring_panini] DSL2 - revision: 8cc832e32f

    executor >  local (6)
    [8c/2e3f91] process > FASTP (3)           [100%] 3 of 3 ✔
    [7a/1b4c92] process > GENERATE_REPORT (3) [100%] 3 of 3 ✔
    ```

La sortie devrait montrer que les deux processus se terminent avec succès. Le workflow est maintenant beaucoup plus propre et plus facile à maintenir, avec toute la logique complexe de traitement des métadonnées encapsulée dans la fonction `separateMetadata`.

### Takeaway

Dans cette section, vous avez appris la **création de fonctions** :

- **Définition de fonctions avec `def`** : Le mot-clé pour créer des fonctions nommées (comme `def` en Python ou `function` en JavaScript)
- **Portée des fonctions** : Les fonctions définies au niveau du script sont accessibles dans tout votre workflow Nextflow
- **Valeurs de retour** : Les fonctions retournent automatiquement la dernière expression, ou utilisent un `return` explicite
- **Code plus propre** : Extraire une logique complexe dans des fonctions est une pratique fondamentale d'ingénierie logicielle dans n'importe quel langage

Ensuite, nous explorerons comment utiliser les closures dans les directives de processus pour une allocation dynamique des ressources.

---

## 4. Directives de ressources dynamiques avec closures

Jusqu'à présent, nous avons utilisé des scripts dans le bloc `script` des processus. Mais les **closures** (introduites dans la section 1.1) sont également incroyablement utiles dans les directives de processus, en particulier pour l'allocation dynamique des ressources. Ajoutons des directives de ressources à notre processus FASTP qui s'adaptent en fonction des caractéristiques de l'échantillon.

### 4.1. Allocation de ressources spécifique aux échantillons

Actuellement, notre processus FASTP utilise des ressources par défaut. Rendons-le plus intelligent en allouant plus de CPUs pour les échantillons de haute profondeur. Modifiez `modules/fastp.nf` pour inclure une directive `cpus` dynamique et une directive `memory` statique :

=== "After"

    ```groovy title="modules/fastp.nf" linenums="1" hl_lines="4-5"
    process FASTP {
        container 'community.wave.seqera.io/library/fastp:0.24.0--62c97b06e8447690'

        cpus { meta.depth > 40000000 ? 2 : 1 }
        memory 2.GB

        input:
        tuple val(meta), path(reads)
    ```

=== "Before"

    ```groovy title="modules/fastp.nf" linenums="1"
    process FASTP {
        container 'community.wave.seqera.io/library/fastp:0.24.0--62c97b06e8447690'

        input:
        tuple val(meta), path(reads)
    ```

La closure `{ meta.depth > 40000000 ? 2 : 1 }` utilise l'**opérateur ternaire** (couvert dans la section 1.1) et est évaluée pour chaque tâche, permettant une allocation de ressources par échantillon. Les échantillons à haute profondeur (>40M lectures) obtiennent 2 CPUs, tandis que les autres obtiennent 1 CPU.

!!! note "Accès aux variables d'entrée dans les directives"

    La closure peut accéder à n'importe quelle variable d'entrée (comme `meta` ici) parce que Nextflow évalue ces closures dans le contexte de chaque exécution de tâche.

Exécutez à nouveau le workflow avec l'option `-ansi-log false` pour qu'il soit plus facile de voir les hachages de tâches.

```bash
nextflow run main.nf -ansi-log false
```

??? success "Sortie de la commande"

    ```console
    N E X T F L O W  ~  version 25.10.2
    Launching `main.nf` [fervent_albattani] DSL2 - revision: fa8f249759
    [bd/ff3d41] Submitted process > FASTP (2)
    [a4/a3aab2] Submitted process > FASTP (1)
    [48/6db0c9] Submitted process > FASTP (3)
    [ec/83439d] Submitted process > GENERATE_REPORT (3)
    [bd/15d7cc] Submitted process > GENERATE_REPORT (2)
    [42/699357] Submitted process > GENERATE_REPORT (1)
    ```

Vous pouvez vérifier la commande `docker` exacte qui a été exécutée pour voir l'allocation de CPU pour n'importe quelle tâche donnée :

```console title="Vérification de la commande docker"
cat work/48/6db0c9e9d8aa65e4bb4936cd3bd59e/.command.run | grep "docker run"
```

Vous devriez voir quelque chose comme :

```bash title="commande docker"
    docker run -i --cpu-shares 4096 --memory 2048m -e "NXF_TASK_WORKDIR" -v /workspaces/training/side-quests/essential_scripting_patterns:/workspaces/training/side-quests/essential_scripting_patterns -w "$NXF_TASK_WORKDIR" --name $NXF_BOXID community.wave.seqera.io/library/fastp:0.24.0--62c97b06e8447690 /bin/bash -ue /workspaces/training/side-quests/essential_scripting_patterns/work/48/6db0c9e9d8aa65e4bb4936cd3bd59e/.command.sh
```

Dans cet exemple, nous avons choisi un exemple qui a demandé 2 CPUs (`--cpu-shares 2048`), car c'était un échantillon de haute profondeur, mais vous devriez voir différentes allocations de CPU en fonction de la profondeur de l'échantillon. Essayez cela pour les autres tâches aussi.

### 4.2. Stratégies de réessai

Un autre modèle puissant consiste à utiliser `task.attempt` pour des stratégies de réessai. Pour montrer pourquoi c'est utile, nous allons commencer par réduire l'allocation de mémoire à FASTP à moins que ce dont il a besoin. Changez la directive `memory` dans `modules/fastp.nf` à `1.GB` :

=== "After"

    ```groovy title="modules/fastp.nf" linenums="1" hl_lines="5"
    process FASTP {
        container 'community.wave.seqera.io/library/fastp:0.24.0--62c97b06e8447690'

        cpus { meta.depth > 40000000 ? 4 : 2 }
        memory 1.GB

        input:
        tuple val(meta), path(reads)
    ```

=== "Before"

    ```groovy title="modules/fastp.nf" linenums="1" hl_lines="5"
    process FASTP {
        container 'community.wave.seqera.io/library/fastp:0.24.0--62c97b06e8447690'

        cpus { meta.depth > 40000000 ? 4 : 2 }
        memory 2.GB

        input:
        tuple val(meta), path(reads)
    ```

... et exécutez à nouveau le workflow :

```bash
nextflow run main.nf
```

??? failure "Sortie de la commande"

    ```console hl_lines="2 11"
    Command exit status:
      137

    Command output:
      (empty)

    Command error:
      Detecting adapter sequence for read1...
      No adapter detected for read1

      .command.sh: line 7:   101 Killed                  fastp --in1 SAMPLE_002_S2_L001_R1_001.fastq --out1 sample_002_trimmed.fastq.gz --json sample_002.fastp.json --html sample_002.fastp.html --thread 2
    ```

Cela indique que le processus a été tué pour avoir dépassé les limites de mémoire.

C'est un scénario très courant dans les workflows du monde réel - parfois vous ne savez tout simplement pas combien de mémoire une tâche aura besoin jusqu'à ce que vous l'exécutiez.

Pour rendre notre workflow plus robuste, nous pouvons mettre en œuvre une stratégie de réessai qui augmente l'allocation de mémoire à chaque tentative, en utilisant encore une fois une closure Groovy. Modifiez la directive `memory` pour multiplier la mémoire de base par `task.attempt`, et ajoutez les directives `errorStrategy 'retry'` et `maxRetries 2` :

=== "After"

    ```groovy title="modules/fastp.nf" linenums="1" hl_lines="5-7"
    process FASTP {
        container 'community.wave.seqera.io/library/fastp:0.24.0--62c97b06e8447690'

        cpus { meta.depth > 40000000 ? 4 : 2 }
        memory { 1.GB * task.attempt }
        errorStrategy 'retry'
        maxRetries 2

        input:
        tuple val(meta), path(reads)
    ```

=== "Before"

    ```groovy title="modules/fastp.nf" linenums="1" hl_lines="5"
    process FASTP {
        container 'community.wave.seqera.io/library/fastp:0.24.0--62c97b06e8447690'

        cpus { meta.depth > 40000000 ? 4 : 2 }
        memory 2.GB

        input:
        tuple val(meta), path(reads)
    ```

Maintenant, si le processus échoue en raison d'une mémoire insuffisante, Nextflow réessaiera avec plus de mémoire :

- Première tentative : 1 GB (task.attempt = 1)
- Deuxième tentative : 2.GB (task.attempt = 2)

... et ainsi de suite, jusqu'à la limite `maxRetries`.

### Takeaway

Les directives dynamiques avec closures vous permettent de :

- Allouer des ressources en fonction des caractéristiques d'entrée
- Mettre en œuvre des stratégies de réessai automatiques avec des ressources croissantes
- Combiner plusieurs facteurs (métadonnées, numéro de tentative, priorités)
- Utiliser une logique conditionnelle pour des calculs de ressources complexes

Cela rend vos workflows à la fois plus efficaces (pas de surallocation) et plus robustes (réessai automatique avec plus de ressources).

---

## 5. Logique conditionnelle et contrôle des processus

Précédemment, nous avons utilisé `.map()` avec du scripting pour transformer des données de channels. Maintenant, nous allons utiliser une logique conditionnelle pour contrôler quels processus s'exécutent en fonction des données — essentiel pour des workflows flexibles s'adaptant à différents types d'échantillons.

Les [opérateurs de flux de données](https://www.nextflow.io/docs/latest/reference/operator.html) de Nextflow prennent des closures évaluées à l'exécution, permettant une logique conditionnelle pour diriger les décisions de workflow basées sur le contenu des channels.

### 5.1. Routage avec `.branch()`

Par exemple, imaginons que nos échantillons de séquençage doivent être découpés avec FASTP uniquement s'il s'agit d'échantillons humains avec une couverture supérieure à un certain seuil. Les échantillons de souris ou les échantillons à faible couverture devraient être exécutés avec Trimgalore à la place (c'est un exemple artificiel, mais il illustre le point).

Nous avons fourni un processus Trimgalore simple dans `modules/trimgalore.nf`, jetez-y un œil si vous voulez, mais les détails ne sont pas importants pour cet exercice. Le point clé est que nous voulons acheminer les échantillons en fonction de leurs métadonnées.

Incluez le nouveau module depuis `modules/trimgalore.nf` :

=== "After"

    ```groovy title="main.nf" linenums="1" hl_lines="2"
    include { FASTP } from './modules/fastp.nf'
    include { TRIMGALORE } from './modules/trimgalore.nf'
    ```

=== "Before"

    ```groovy title="main.nf" linenums="1"
    include { FASTP } from './modules/fastp.nf'
    ```

... puis modifiez votre workflow `main.nf` pour diviser les échantillons en fonction de leurs métadonnées et les acheminer via le processus de découpage approprié, comme ceci :

=== "After"

    ```groovy title="main.nf" linenums="28" hl_lines="5-12"
        ch_samples = channel.fromPath("./data/samples.csv")
            .splitCsv(header: true)
            .map { row -> separateMetadata(row) }

        trim_branches = ch_samples
            .branch { meta, reads ->
                fastp: meta.organism == 'human' && meta.depth >= 30000000
                trimgalore: true
            }

        ch_fastp = FASTP(trim_branches.fastp)
        ch_trimgalore = TRIMGALORE(trim_branches.trimgalore)
        GENERATE_REPORT(ch_samples)
    ```

=== "Before"

    ```groovy title="main.nf" linenums="28" hl_lines="5"
        ch_samples = channel.fromPath("./data/samples.csv")
            .splitCsv(header: true)
            .map { row -> separateMetadata(row) }

        ch_fastp = FASTP(ch_samples)
        GENERATE_REPORT(ch_samples)
    ```

Exécutez ce workflow modifié :

```bash
nextflow run main.nf
```

??? success "Sortie de la commande"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `main.nf` [adoring_galileo] DSL2 - revision: c9e83aaef1

    executor >  local (6)
    [1d/0747ac] process > FASTP (2)           [100%] 2 of 2 ✔
    [cc/c44caf] process > TRIMGALORE (1)      [100%] 1 of 1 ✔
    [34/bd5a9f] process > GENERATE_REPORT (1) [100%] 3 of 3 ✔
    ```

Ici, nous avons utilisé de petites mais puissantes expressions conditionnelles à l'intérieur de l'opérateur `.branch{}` pour acheminer les échantillons en fonction de leurs métadonnées. Les échantillons humains avec une couverture élevée passent par `FASTP`, tandis que tous les autres échantillons passent par `TRIMGALORE`.

### 5.2. Utiliser `.filter()` avec Truthiness

Un autre modèle puissant pour contrôler l'exécution du workflow est l'opérateur `.filter()`, qui utilise une closure pour déterminer quels éléments doivent continuer dans le pipeline. À l'intérieur de la closure de filtre, vous écrirez des **expressions booléennes** qui décident quels éléments passent.

Nextflow (comme beaucoup de langages dynamiques) a un concept de **"truthiness"** qui détermine quelles valeurs sont évaluées à `true` ou `false` dans des contextes booléens :

- **Truthy** : valeurs non nulles, chaînes non vides, nombres non nuls, collections non vides
- **Falsy** : `null`, chaînes vides `""`, zéro `0`, collections vides `[]` ou `[:]`, `false`

Cela signifie que `meta.id` seul (sans un explicite `!= null`) vérifie si l'ID existe et n'est pas vide. Utilisons cela pour filtrer les échantillons qui ne répondent pas à nos exigences de qualité.

Ajoutez ce qui suit avant l'opération de branchement :

=== "After"

    ```groovy title="main.nf" linenums="28" hl_lines="5-11"
        ch_samples = channel.fromPath("./data/samples.csv")
            .splitCsv(header: true)
            .map { row -> separateMetadata(row) }

        // Filtrer les échantillons invalides ou de faible qualité
        ch_valid_samples = ch_samples
            .filter { meta, reads ->
                meta.id && meta.organism && meta.depth >= 25000000
            }

        trim_branches = ch_valid_samples
            .branch { meta, reads ->
                fastp: meta.organism == 'human' && meta.depth >= 30000000
                trimgalore: true
            }
    ```

=== "Before"

    ```groovy title="main.nf" linenums="28" hl_lines="5"
        ch_samples = channel.fromPath("./data/samples.csv")
            .splitCsv(header: true)
            .map { row -> separateMetadata(row) }

        trim_branches = ch_samples
            .branch { meta, reads ->
                fastp: meta.organism == 'human' && meta.depth >= 30000000
                trimgalore: true
            }
    ```

Exécutez à nouveau le workflow :

```bash
nextflow run main.nf
```

??? success "Sortie de la commande"

    ```console
    N E X T F L O W  ~  version 25.10.2
    Launching `main.nf` [lonely_williams] DSL2 - revision: d0b3f121ec
    [94/b48eac] Submitted process > FASTP (2)
    [2c/d2b28f] Submitted process > GENERATE_REPORT (2)
    [65/2e3be4] Submitted process > GENERATE_REPORT (1)
    [94/b48eac] NOTE: Process `FASTP (2)` terminated with an error exit status (137) -- Execution is retried (1)
    [3e/0d8664] Submitted process > TRIMGALORE (1)
    [6a/9137b0] Submitted process > FASTP (1)
    [6a/9137b0] NOTE: Process `FASTP (1)` terminated with an error exit status (137) -- Execution is retried (1)
    [83/577ac0] Submitted process > GENERATE_REPORT (3)
    [a2/5117de] Re-submitted process > FASTP (1)
    [1f/a1a4ca] Re-submitted process > FASTP (2)
    ```

Parce que nous avons choisi un filtre qui exclut certains échantillons, moins de tâches ont été exécutées.

L'expression de filtre `meta.id && meta.organism && meta.depth >= 25000000` combine truthiness avec des comparaisons explicites :

- `meta.id && meta.organism` vérifie que les deux champs existent et ne sont pas vides (en utilisant truthiness)
- `meta.depth >= 25000000` assure une profondeur de séquençage suffisante avec une comparaison explicite

!!! note "Truthiness en pratique"

    L'expression `meta.id && meta.organism` est plus concise que d'écrire :
    ```groovy
    meta.id != null && meta.id != '' && meta.organism != null && meta.organism != ''
    ```

    Cela rend la logique de filtrage beaucoup plus propre et plus facile à lire.

### Takeaway

Dans cette section, vous avez appris à utiliser une logique conditionnelle pour contrôler l'exécution du workflow à l'aide des interfaces de closure des opérateurs Nextflow comme `.branch{}` et `.filter{}`, en exploitant truthiness pour écrire des expressions conditionnelles concises.

Notre pipeline achemine maintenant intelligemment les échantillons via des processus appropriés, mais les workflows de production doivent gérer gracieusement les données invalides. Rendons notre workflow robuste contre les valeurs manquantes ou nulles.

---

## 6. Opérateurs de navigation sûre et Elvis

Notre fonction `separateMetadata` suppose actuellement que tous les champs CSV sont présents et valides. Mais que se passe-t-il avec des données incomplètes ? Découvrons-le.

### 6.1. Le problème : Accéder à des propriétés qui n'existent pas

Disons que nous voulons ajouter un support pour les informations optionnelles d'exécution de séquençage. Dans certains laboratoires, les échantillons pourraient avoir un champ supplémentaire pour l'ID de l'exécution de séquençage ou le numéro de lot, mais notre CSV actuel n'a pas cette colonne. Essayons d'y accéder quand même.

Modifiez la fonction `separateMetadata` pour inclure un champ run_id :

=== "After"

    ```groovy title="main.nf" linenums="5" hl_lines="9"
    def separateMetadata(row) {
        def sample_meta = [
            id: row.sample_id.toLowerCase(),
            organism: row.organism,
            tissue: row.tissue_type.replaceAll('_', ' ').toLowerCase(),
            depth: row.sequencing_depth.toInteger(),
            quality: row.quality_score.toDouble()
        ]
        def run_id = row.run_id.toUpperCase()
    ```

=== "Before"

    ```groovy title="main.nf" linenums="5"
    def separateMetadata(row) {
        def sample_meta = [
            id: row.sample_id.toLowerCase(),
            organism: row.organism,
            tissue: row.tissue_type.replaceAll('_', ' ').toLowerCase(),
            depth: row.sequencing_depth.toInteger(),
            quality: row.quality_score.toDouble()
        ]
    ```

Maintenant, exécutez le workflow :

```bash
nextflow run main.nf
```

??? failure "Sortie de la commande"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `main.nf` [trusting_torvalds] DSL2 - revision: b56fbfbce2

    ERROR ~ Cannot invoke method toUpperCase() on null object

    -- Check script 'main.nf' at line: 13 or see '.nextflow.log' file for more details
    ```

Cela plante avec une NullPointerException.

Le problème est que `row.run_id` retourne `null` car la colonne `run_id` n'existe pas dans notre CSV. Quand nous essayons d'appeler `.toUpperCase()` sur `null`, cela plante. C'est là que l'opérateur de navigation sûre entre en jeu.

### 6.2. Opérateur de navigation sûre (`?.`)

L'opérateur de navigation sûre (`?.`) retourne `null` au lieu de lancer une exception quand il est appelé sur une valeur `null`. Si l'objet avant `?.` est `null`, l'expression entière est évaluée à `null` sans exécuter la méthode.

Mettez à jour la fonction pour utiliser la navigation sûre :

=== "After"

    ```groovy title="main.nf" linenums="4" hl_lines="9"
    def separateMetadata(row) {
        def sample_meta = [
            id: row.sample_id.toLowerCase(),
            organism: row.organism,
            tissue: row.tissue_type.replaceAll('_', ' ').toLowerCase(),
            depth: row.sequencing_depth.toInteger(),
            quality: row.quality_score.toDouble()
        ]
        def run_id = row.run_id?.toUpperCase()
    ```

=== "Before"

    ```groovy title="main.nf" linenums="4" hl_lines="9"
    def separateMetadata(row) {
        def sample_meta = [
            id: row.sample_id.toLowerCase(),
            organism: row.organism,
            tissue: row.tissue_type.replaceAll('_', ' ').toLowerCase(),
            depth: row.sequencing_depth.toInteger(),
            quality: row.quality_score.toDouble()
        ]
        def run_id = row.run_id.toUpperCase()
    ```

Exécutez à nouveau :

```bash
nextflow run main.nf
```

??? success "Sortie de la commande"

    ```console
    <!-- TODO: output -->
    ```

Pas de crash ! Le workflow gère maintenant gracieusement le champ manquant. Quand `row.run_id` est `null`, l'opérateur `?.` empêche l'appel `.toUpperCase()`, et `run_id` devient `null` au lieu de provoquer une exception.

### 6.3. Opérateur Elvis (`?:`) pour les valeurs par défaut

L'opérateur Elvis (`?:`) fournit des valeurs par défaut lorsque le côté gauche est "falsy" (comme expliqué précédemment). Il est nommé d'après Elvis Presley car `?:` ressemble à sa célèbre coiffure et ses yeux lorsqu'on le regarde de côté !

Maintenant que nous utilisons la navigation sûre, `run_id` sera `null` pour les échantillons sans ce champ. Utilisons l'opérateur Elvis pour fournir une valeur par défaut et l'ajouter à notre map `sample_meta` :

=== "After"

    ```groovy title="main.nf" linenums="5" hl_lines="9-10"
    def separateMetadata(row) {
        def sample_meta = [
            id: row.sample_id.toLowerCase(),
            organism: row.organism,
            tissue: row.tissue_type.replaceAll('_', ' ').toLowerCase(),
            depth: row.sequencing_depth.toInteger(),
            quality: row.quality_score.toDouble()
        ]
        def run_id = row.run_id?.toUpperCase() ?: 'NON_SPECIFIE'
        sample_meta.run = run_id
    ```

=== "Before"

    ```groovy title="main.nf" linenums="5" hl_lines="9"
    def separateMetadata(row) {
        def sample_meta = [
            id: row.sample_id.toLowerCase(),
            organism: row.organism,
            tissue: row.tissue_type.replaceAll('_', ' ').toLowerCase(),
            depth: row.sequencing_depth.toInteger(),
            quality: row.quality_score.toDouble()
        ]
        def run_id = row.run_id?.toUpperCase()
    ```

Ajoutez également un opérateur `view()` dans le workflow pour voir les résultats :

=== "After"

    ```groovy title="main.nf" linenums="30" hl_lines="4"
        ch_samples = channel.fromPath("./data/samples.csv")
            .splitCsv(header: true)
            .map{ row -> separateMetadata(row) }
            .view()
    ```

=== "Before"

    ```groovy title="main.nf" linenums="30"
        ch_samples = channel.fromPath("./data/samples.csv")
            .splitCsv(header: true)
            .map{ row -> separateMetadata(row) }
    ```

et exécutez le workflow :

```bash
nextflow run main.nf
```

??? success "Sortie de la commande"

    ```console
    [[id:sample_001, organism:human, tissue:liver, depth:30000000, quality:38.5, run:NON_SPECIFIE, sample_num:1, lane:001, read:R1, chunk:001, priority:normal], /workspaces/training/side-quests/essential_scripting_patterns/data/sequences/SAMPLE_001_S1_L001_R1_001.fastq]
    [[id:sample_002, organism:mouse, tissue:brain, depth:25000000, quality:35.2, run:NON_SPECIFIE, sample_num:2, lane:001, read:R1, chunk:001, priority:normal], /workspaces/training/side-quests/essential_scripting_patterns/data/sequences/SAMPLE_002_S2_L001_R1_001.fastq]
    [[id:sample_003, organism:human, tissue:kidney, depth:45000000, quality:42.1, run:NON_SPECIFIE, sample_num:3, lane:001, read:R1, chunk:001, priority:high], /workspaces/training/side-quests/essential_scripting_patterns/data/sequences/SAMPLE_003_S3_L001_R1_001.fastq]
    ```

Parfait ! Maintenant tous les échantillons ont un champ `run` avec soit leur ID d'exécution réel (en majuscules), soit la valeur par défaut 'NON_SPECIFIE'. La combinaison de `?.` et `?:` fournit à la fois la sécurité (pas de crashes) et des valeurs par défaut sensées.

Supprimez l'opérateur `.view()` maintenant que nous avons confirmé qu'il fonctionne.

!!! tip "Combinaison de la navigation sûre et d'Elvis"

    Le modèle `value?.method() ?: 'default'` est courant dans les workflows de production :

    - `value?.method()` - Appelle la méthode en toute sécurité, retourne `null` si `value` est `null`
    - `?: 'default'` - Fournit une valeur de repli si le résultat est `null`

    Ce modèle gère gracieusement les données manquantes/incomplètes.

Utilisez ces opérateurs de manière cohérente dans les fonctions, les closures d'opérateurs (`.map{}`, `.filter{}`), les scripts de processus et les fichiers de configuration. Ils empêchent les crashes lors de la manipulation de données du monde réel.

### Takeaway

- **Navigation sûre (`?.`)** : Empêche les crashes sur des valeurs nulles - retourne null au lieu de lancer une exception
- **Opérateur Elvis (`?:`)** : Fournit des valeurs par défaut - `value ?: 'default'`
- **Combinaison** : `value?.method() ?: 'default'` est le modèle courant

Ces opérateurs rendent les workflows résistants aux données incomplètes - essentiel pour le travail dans le monde réel.

---

## 7. Validation avec `error()` et `log.warn`

Parfois, vous devez arrêter le workflow immédiatement si les paramètres d'entrée sont invalides. Dans Nextflow, vous pouvez utiliser des fonctions intégrées comme `error()` et `log.warn`, ainsi que des constructions de programmation standard comme les instructions `if` et la logique booléenne, pour implémenter une logique de validation. Ajoutons une validation à notre workflow.

Créez une fonction de validation avant votre bloc workflow, appelez-la depuis le workflow, et changez la création du channel pour utiliser un paramètre pour le chemin du fichier CSV. Si le paramètre est manquant ou si le fichier n'existe pas, appelez `error()` pour arrêter l'exécution avec un message clair.

=== "After"

    ```groovy title="main.nf" linenums="1" hl_lines="5-20 23-24"
    include { FASTP } from './modules/fastp.nf'
    include { TRIMGALORE } from './modules/trimgalore.nf'
    include { GENERATE_REPORT } from './modules/generate_report.nf'

    def validateInputs() {
        // Vérification que le paramètre d'entrée est fourni
        if (!params.input) {
            error("Chemin du fichier CSV d'entrée non fourni. Veuillez spécifier --input <fichier.csv>")
        }

        // Vérification que le fichier CSV existe
        if (!file(params.input).exists()) {
            error("Fichier CSV d'entrée non trouvé : ${params.input}")
        }
    }
    ...
    workflow {
        validateInputs()
        ch_samples = channel.fromPath(params.input)
    ```

=== "Before"

    ```groovy title="main.nf" linenums="1"
    include { FASTP } from './modules/fastp.nf'
    include { TRIMGALORE } from './modules/trimgalore.nf'
    include { GENERATE_REPORT } from './modules/generate_report.nf'

    ...
    workflow {
        ch_samples = channel.fromPath("./data/samples.csv")
    ```

Maintenant, essayez d'exécuter sans le fichier CSV :

```bash
nextflow run main.nf
```

??? failure "Sortie de la commande"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `main.nf` [confident_coulomb] DSL2 - revision: 07059399ed

    WARN: Access to undefined parameter `input` -- Initialise it to a default value eg. `params.input = some_value`
    Chemin du fichier CSV d'entrée non fourni. Veuillez spécifier --input <fichier.csv>
    ```

Le workflow s'arrête immédiatement avec un message d'erreur clair au lieu d'échouer mystérieusement plus tard

Maintenant, exécutez-le avec un fichier inexistant :

```bash
nextflow run main.nf --input ./data/nonexistent.csv
```

??? failure "Sortie de la commande"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `main.nf` [cranky_gates] DSL2 - revision: 26839ae3eb

    Fichier CSV d'entrée non trouvé : ./data/nonexistent.csv
    ```

Enfin, exécutez-le avec le bon fichier :

```bash
nextflow run main.nf --input ./data/samples.csv
```

??? success "Sortie de la commande"

    ```console
    <!-- TODO: output -->
    ```

Cette fois, il s'exécute avec succès.

Vous pouvez également ajouter une validation dans la fonction `separateMetadata`. Utilisons le non-fatal `log.warn` pour émettre des avertissements pour les échantillons avec une faible profondeur de séquençage, mais permettons quand même au workflow de continuer :

=== "After"

    ```groovy title="main.nf" linenums="1" hl_lines="3-6"
        def priority = sample_meta.quality > 40 ? 'high' : 'normal'

        // Validation des données pour s'assurer qu'elles sont cohérentes
        if (sample_meta.depth < 30000000) {
            log.warn "Faible profondeur de séquençage pour ${sample_meta.id} : ${sample_meta.depth}"
        }

        return tuple(sample_meta + file_meta + [priority: priority], fastq_path)
    }
    ```

=== "Before"

    ```groovy title="main.nf" linenums="1"
        def priority = sample_meta.quality > 40 ? 'high' : 'normal'

        return tuple(sample_meta + file_meta + [priority: priority], fastq_path)
    }
    ```

Exécutez à nouveau le workflow avec le CSV original :

```bash
nextflow run main.nf --input ./data/samples.csv
```

??? warning "Sortie de la commande"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `main.nf` [awesome_goldwasser] DSL2 - revision: a31662a7c1

    executor >  local (5)
    [ce/df5eeb] process > FASTP (2)           [100%] 2 of 2 ✔
    [-        ] process > TRIMGALORE          -
    [d1/7d2b4b] process > GENERATE_REPORT (3) [100%] 3 of 3 ✔
    WARN: Faible profondeur de séquençage pour sample_002 : 25000000
    ```

Nous voyons un avertissement concernant une faible profondeur de séquençage pour l'un des échantillons.

### Takeaway

- **`error()`** : Arrête le workflow immédiatement avec un message clair
- **`log.warn`** : Émet des avertissements sans arrêter le workflow
- **Validation précoce** : Vérifiez les entrées avant le traitement pour échouer rapidement avec des erreurs utiles
- **Fonctions de validation** : Créez une logique de validation réutilisable qui peut être appelée au début du workflow

Une validation appropriée rend les workflows plus robustes et conviviaux en détectant les problèmes tôt avec des messages d'erreur clairs.

---

## 8. Gestionnaires d'événements de workflow

Jusqu'à présent, nous avons écrit du code dans nos scripts de workflow et définitions de processus. Mais il y a une autre fonctionnalité importante que vous devriez connaître : les gestionnaires d'événements de workflow.

Les gestionnaires d'événements sont des closures qui s'exécutent à des moments spécifiques du cycle de vie de votre workflow. Ils sont parfaits pour ajouter de la journalisation, des notifications ou des opérations de nettoyage. Ces gestionnaires doivent être définis dans votre script de workflow à côté de votre définition de workflow.

### 8.1. Le gestionnaire `onComplete`

Le gestionnaire d'événements le plus couramment utilisé est `onComplete`, qui s'exécute lorsque votre workflow se termine (qu'il ait réussi ou échoué). Ajoutons-en un pour résumer les résultats de notre pipeline.

Ajoutez le gestionnaire d'événements à votre fichier `main.nf`, à l'intérieur de votre définition de workflow :

=== "After"

    ```groovy title="main.nf" linenums="66" hl_lines="5-16"
        ch_fastp = FASTP(trim_branches.fastp)
        ch_trimgalore = TRIMGALORE(trim_branches.trimgalore)
        GENERATE_REPORT(ch_samples)

        workflow.onComplete = {
            println ""
            println "Résumé d'exécution du pipeline :"
            println "=========================="
            println "Terminé à : ${workflow.complete}"
            println "Durée     : ${workflow.duration}"
            println "Succès    : ${workflow.success}"
            println "workDir   : ${workflow.workDir}"
            println "statut sortie : ${workflow.exitStatus}"
            println ""
        }
    }
    ```

=== "Before"

    ```groovy title="main.nf" linenums="66" hl_lines="4"
        ch_fastp = FASTP(trim_branches.fastp)
        ch_trimgalore = TRIMGALORE(trim_branches.trimgalore)
        GENERATE_REPORT(ch_samples)
    }
    ```

Cette closure s'exécute lorsque le workflow se termine. À l'intérieur, vous avez accès à l'objet `workflow` qui fournit des propriétés utiles sur l'exécution.

Exécutez votre workflow et vous verrez ce résumé apparaître à la fin !

```bash
nextflow run main.nf --input ./data/samples.csv -ansi-log false
```

??? success "Sortie de la commande"

    ```console
    N E X T F L O W  ~  version 25.10.2
    Launching `main.nf` [marvelous_boltzmann] DSL2 - revision: a31662a7c1
    WARN: Faible profondeur de séquençage pour sample_002 : 25000000
    [9b/d48e40] Submitted process > FASTP (2)
    [6a/73867a] Submitted process > GENERATE_REPORT (2)
    [79/ad0ac5] Submitted process > GENERATE_REPORT (1)
    [f3/bda6cb] Submitted process > FASTP (1)
    [34/d5b52f] Submitted process > GENERATE_REPORT (3)

    Résumé d'exécution du pipeline :
    ==========================
    Terminé à : 2025-10-10T12:14:24.885384+01:00
    Durée     : 2.9s
    Succès    : true
    workDir   : /workspaces/training/side-quests/essential_scripting_patterns/work
    statut sortie : 0
    ```

Rendons-le plus utile en ajoutant une logique conditionnelle :

=== "After"

    ```groovy title="main.nf" linenums="66" hl_lines="5-22"
        ch_fastp = FASTP(trim_branches.fastp)
        ch_trimgalore = TRIMGALORE(trim_branches.trimgalore)
        GENERATE_REPORT(ch_samples)

        workflow.onComplete = {
            println ""
            println "Résumé d'exécution du pipeline :"
            println "=========================="
            println "Terminé à : ${workflow.complete}"
            println "Durée     : ${workflow.duration}"
            println "Succès    : ${workflow.success}"
            println "workDir   : ${workflow.workDir}"
            println "statut sortie : ${workflow.exitStatus}"
            println ""

            if (workflow.success) {
                println "✅ Pipeline terminé avec succès !"
            } else {
                println "❌ Pipeline échoué !"
                println "Erreur : ${workflow.errorMessage}"
            }
        }
    }
    ```

=== "Before"

    ```groovy title="main.nf" linenums="66" hl_lines="5-16"
        ch_fastp = FASTP(trim_branches.fastp)
        ch_trimgalore = TRIMGALORE(trim_branches.trimgalore)
        GENERATE_REPORT(ch_samples)

        workflow.onComplete = {
            println ""
            println "Résumé d'exécution du pipeline :"
            println "=========================="
            println "Terminé à : ${workflow.complete}"
            println "Durée     : ${workflow.duration}"
            println "Succès    : ${workflow.success}"
            println "workDir   : ${workflow.workDir}"
            println "statut sortie : ${workflow.exitStatus}"
            println ""
        }
    }
    ```

Maintenant, nous obtenons un résumé encore plus informatif, y compris un message de succès/échec et le répertoire de sortie si spécifié :

<!-- TODO: add run command -->

??? success "Sortie de la commande"

    ```console
    N E X T F L O W  ~  version 25.10.2
    Launching `main.nf` [boring_linnaeus] DSL2 - revision: a31662a7c1
    WARN: Faible profondeur de séquençage pour sample_002 : 25000000
    [e5/242efc] Submitted process > FASTP (2)
    [3b/74047c] Submitted process > GENERATE_REPORT (3)
    [8a/7a57e6] Submitted process > GENERATE_REPORT (1)
    [a8/b1a31f] Submitted process > GENERATE_REPORT (2)
    [40/648429] Submitted process > FASTP (1)

    Résumé d'exécution du pipeline :
    ==========================
    Terminé à : 2025-10-10T12:16:00.522569+01:00
    Durée     : 3.6s
    Succès    : true
    workDir   : /workspaces/training/side-quests/essential_scripting_patterns/work
    statut sortie : 0

    ✅ Pipeline terminé avec succès !
    ```

Vous pouvez également écrire le résumé dans un fichier en utilisant des opérations de fichier :

```groovy title="main.nf - Écriture du résumé dans un fichier"
workflow {
    // ... votre code de workflow ...

    workflow.onComplete = {
        def summary = """
        Résumé d'exécution du pipeline
        ===========================
        Terminé : ${workflow.complete}
        Durée   : ${workflow.duration}
        Succès  : ${workflow.success}
        Commande : ${workflow.commandLine}
        """

        println summary

        // Écrire dans un fichier journal
        def log_file = file("${workflow.launchDir}/pipeline_summary.txt")
        log_file.text = summary
    }
}
```

### 8.2. Le gestionnaire `onError`

En plus de `onComplete`, il y a un autre gestionnaire d'événements que vous pouvez utiliser : `onError`, qui s'exécute uniquement si le workflow échoue :

```groovy title="main.nf - gestionnaire onError"
workflow {
    // ... votre code de workflow ...

    workflow.onError = {
        println "="* 50
        println "Échec de l'exécution du pipeline !"
        println "Message d'erreur : ${workflow.errorMessage}"
        println "="* 50

        // Écrire un journal d'erreur détaillé
        def error_file = file("${workflow.launchDir}/error.log")
        error_file.text = """
        Rapport d'erreur du workflow
        =====================
        Heure : ${new Date()}
        Erreur : ${workflow.errorMessage}
        Rapport d'erreur : ${workflow.errorReport ?: 'Pas de rapport détaillé disponible'}
        """

        println "Détails de l'erreur écrits dans : ${error_file}"
    }
}
```

Vous pouvez utiliser plusieurs gestionnaires ensemble dans votre script de workflow :

```groovy title="main.nf - Gestionnaires combinés"
workflow {
    // ... votre code de workflow ...

    workflow.onError = {
        println "Workflow échoué : ${workflow.errorMessage}"
    }

    workflow.onComplete = {
        def duration_mins = workflow.duration.toMinutes().round(2)
        def status = workflow.success ? "SUCCÈS ✅" : "ÉCHEC ❌"

        println """
        Pipeline terminé : ${status}
        Durée : ${duration_mins} minutes
        """
    }
}
```

### Takeaway

Dans cette section, vous avez appris :

- **Closures de gestionnaires d'événements** : Closures dans votre script de workflow qui s'exécutent à différents moments du cycle de vie
- **Gestionnaire `onComplete`** : Pour les résumés d'exécution et les rapports de résultats
- **Gestionnaire `onError`** : Pour la gestion des erreurs et l'enregistrement des échecs
- **Propriétés de l'objet workflow** : Accéder à `workflow.success`, `workflow.duration`, `workflow.errorMessage`, etc.

Les gestionnaires d'événements montrent comment vous pouvez utiliser toute la puissance du langage Nextflow dans vos scripts de workflow pour ajouter des capacités sophistiquées de journalisation et de notification.

---

## Résumé

Félicitations, vous avez réussi !

Tout au long de cette quête secondaire, vous avez construit un pipeline de traitement d'échantillons complet qui a évolué depuis la gestion basique des métadonnées jusqu'à un workflow sophistiqué et prêt pour la production.
Chaque section s'est appuyée sur la précédente, démontrant comment les constructions de programmation transforment des workflows simples en systèmes puissants de traitement de données, avec les avantages suivants :

- **Code plus clair** : Comprendre le flux de données vs scripting vous aide à écrire des workflows plus organisés
- **Gestion robuste** : Les opérateurs de navigation sûre et Elvis rendent les workflows résistants aux données manquantes
- **Traitement flexible** : La logique conditionnelle permet à vos workflows de traiter différents types d'échantillons de manière appropriée
- **Ressources adaptatives** : Les directives dynamiques optimisent l'utilisation des ressources en fonction des caractéristiques d'entrée

Cette progression reflète l'évolution dans le monde réel des pipelines bioinformatiques, des prototypes de recherche qui gèrent quelques échantillons aux systèmes de production qui traitent des milliers d'échantillons à travers les laboratoires et les institutions.
Chaque défi que vous avez résolu et chaque modèle que vous avez appris reflète des problèmes réels auxquels les développeurs sont confrontés lors de la mise à l'échelle des workflows Nextflow.

L'application de ces modèles dans votre propre travail vous permettra de construire des workflows robustes et prêts pour la production.

### Modèles clés

1.  **Flux de données vs Scripting :** Vous avez appris à distinguer entre les opérations de flux de données (orchestration de channels) et le scripting (code qui manipule les données), y compris les différences cruciales entre les opérations sur différents types comme `collect` sur Channel vs List.

    - Flux de données : orchestration de channels

    ```groovy
    channel.fromPath('*.fastq').splitCsv(header: true)
    ```

    - Scripting : traitement de données sur des collections

    ```groovy
    sample_data.collect { it.toUpperCase() }
    ```

2.  **Traitement avancé des chaînes** : Vous avez maîtrisé les expressions régulières pour l'analyse des noms de fichiers, la génération dynamique de scripts dans les processus, et l'interpolation de variables (Nextflow vs Bash vs Shell).

    - Correspondance de motifs

    ```groovy
    filename =~ ~/^(\w+)_(\w+)_(\d+)\.fastq$/
    ```

    - Fonction avec retour conditionnel

    ```groovy
    def parseSample(filename) {
        def matcher = filename =~ pattern
        return matcher ? [valid: true, data: matcher[0]] : [valid: false]
    }
    ```

    - Collection de fichiers vers arguments de commande (dans un bloc de script de processus)

    ```groovy
    script:
    def file_args = input_files.collect { file -> "--input ${file}" }.join(' ')
    """
    analysis_tool ${file_args} --output results.txt
    """
    ```

3.  **Création de fonctions réutilisables** : Vous avez appris à extraire une logique complexe dans des fonctions nommées qui peuvent être appelées depuis des opérateurs de channels, rendant les workflows plus lisibles et maintenables.

    - Définir une fonction nommée

    ```groovy
    def separateMetadata(row) {
        def sample_meta = [ /* code caché pour brièveté */ ]
        def fastq_path = file(row.file_path)
        def m = (fastq_path.name =~ /^(.+)_S(\d+)_L(\d{3})_(R[12])_(\d{3})\.fastq(?:\.gz)?$/)
        def file_meta = m ? [ /* code caché pour brièveté */ ] : [:]
        def priority = sample_meta.quality > 40 ? 'high' : 'normal'

        return tuple(sample_meta + file_meta + [priority: priority], fastq_path)
    }
    ```

    - Appeler la fonction nommée dans un workflow

    ```groovy
    workflow {
        ch_samples = channel.fromPath("./data/samples.csv")
            .splitCsv(header: true)
            .map{ row -> separateMetadata(row) }

        ch_fastp = FASTP(ch_samples)
    }
    ```

4.  **Directives de ressources dynamiques avec closures** : Vous avez exploré l'utilisation de closures dans les directives de processus pour une allocation de ressources adaptative basée sur les caractéristiques d'entrée.

    - Closures nommées et composition

    ```groovy
    def enrichData = normalizeId >> addQualityCategory >> addFlags
    def processor = generalFunction.curry(fixedParam)
    ```

    - Closures avec accès à la portée

    ```groovy
    def collectStats = { data -> stats.count++; return data }
    ```

5.  **Logique conditionnelle et contrôle de processus** : Vous avez ajouté un routage intelligent en utilisant les opérateurs `.branch()` et `.filter()`, en exploitant truthiness pour des expressions conditionnelles concises.

    - Utilisez `.branch()` pour router les données à travers différentes branches de workflow

    ```groovy
    trim_branches = ch_samples
    .branch { meta, reads ->
        fastp: meta.organism == 'human' && meta.depth >= 30000000
        trimgalore: true
    }

    ch_fastp = FASTP(trim_branches.fastp)
    ch_trimgalore = TRIMGALORE(trim_branches.trimgalore)
    ```

    - Évaluation booléenne avec Groovy Truth

    ```groovy
    if (sample.files) println "A des fichiers"
    ```

    - Utilisez `filter()` pour sous-ensemble des données avec 'truthiness'

    ```groovy
    ch_valid_samples = ch_samples
        .filter { meta, reads ->
            meta.id && meta.organism && meta.depth >= 25000000
        }
    ```

6.  **Opérateurs de navigation sûre et Elvis** : Vous avez rendu le pipeline robuste contre les données manquantes en utilisant `?.` pour un accès sûr aux propriétés et `?:` pour fournir des valeurs par défaut.

    ```groovy
    def id = data?.sample?.id ?: 'inconnu'
    ```

7.  **Validation avec error() et log.warn** : Vous avez appris à valider les entrées tôt et à échouer rapidement avec des messages d'erreur clairs.

    ```groovy
    try {
        def errors = validateSample(sample)
        if (errors) throw new RuntimeException("Invalide : ${errors.join(', ')}")
    } catch (Exception e) {
        println "Erreur : ${e.message}"
    }
    ```

8.  **Gestionnaires d'événements de configuration** : Vous avez appris à utiliser des gestionnaires d'événements de workflow (`onComplete` et `onError`) pour la journalisation, les notifications et la gestion du cycle de vie.

    - Utilisation de `onComplete` pour journaliser et notifier

    ```groovy
    workflow.onComplete = {
        println "Succès    : ${workflow.success}"
        println "statut sortie : ${workflow.exitStatus}"

        if (workflow.success) {
            println "✅ Pipeline terminé avec succès !"
        } else {
            println "❌ Pipeline échoué !"
            println "Erreur : ${workflow.errorMessage}"
        }
    }
    ```

    - Utilisation de `onError` pour prendre des mesures spécifiquement en cas d'échec

    ```groovy
    workflow.onError = {
        // Écrire un journal d'erreur détaillé
        def error_file = file("${workflow.launchDir}/error.log")
        error_file.text = """
        Heure : ${new Date()}
        Erreur : ${workflow.errorMessage}
        Rapport d'erreur : ${workflow.errorReport ?: 'Pas de rapport détaillé disponible'}
        """

        println "Détails de l'erreur écrits dans : ${error_file}"
    }
    ```

### Ressources supplémentaires

- [Référence du langage Nextflow](https://nextflow.io/docs/latest/reference/syntax.html)
- [Opérateurs Nextflow](https://www.nextflow.io/docs/latest/operator.html)
- [Syntaxe des scripts Nextflow](https://www.nextflow.io/docs/latest/script.html)
- [Bibliothèque standard Nextflow](https://nextflow.io/docs/latest/reference/stdlib.html)

Assurez-vous de consulter ces ressources lorsque vous aurez besoin d'explorer des fonctionnalités plus avancées.

Vous bénéficierez de la pratique et de l'expansion de vos compétences pour :

- Écrire des workflows plus propres avec une séparation appropriée entre flux de données et scripting
- Maîtriser l'interpolation de variables pour éviter les pièges courants avec les variables Nextflow, Bash et shell
- Utiliser des directives de ressources dynamiques pour des workflows efficaces et adaptatifs
- Transformer des collections de fichiers en arguments de ligne de commande correctement formatés
- Gérer différentes conventions de nommage de fichiers et formats d'entrée avec élégance en utilisant regex et traitement de chaînes
- Construire un code réutilisable et maintenable en utilisant des modèles avancés de closure et de programmation fonctionnelle
- Traiter et organiser des ensembles de données complexes en utilisant des opérations de collection
- Ajouter validation, gestion d'erreurs et journalisation pour rendre vos workflows prêts pour la production
- Implémenter la gestion du cycle de vie du workflow avec des gestionnaires d'événements

---

## Et ensuite ?

Retournez au [menu des quêtes secondaires](./index.md) ou cliquez sur le bouton en bas à droite de la page pour passer au sujet suivant de la liste.
