# Partie 1 : Appel de variants par échantillon

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Traduction assistée par IA - [en savoir plus et suggérer des améliorations](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

Dans la première partie de ce cours, nous vous montrons comment construire un pipeline simple d'appel de variants qui applique l'appel de variants GATK à des échantillons de séquençage individuels.

### Vue d'ensemble de la méthode

L'appel de variants est une méthode d'analyse génomique qui vise à identifier les variations dans une séquence génomique par rapport à un génome de référence.
Ici, nous allons utiliser des outils et des méthodes conçus pour appeler les variants courts, _c'est-à-dire_ les SNP et les indels.

![Pipeline GATK](img/gatk-pipeline.png)

Un pipeline complet d'appel de variants implique généralement de nombreuses étapes, y compris le mapping sur la référence (parfois appelé alignement du génome) et le filtrage et la priorisation des variants.
Pour simplifier, dans cette partie du cours, nous allons nous concentrer uniquement sur la partie appel de variants.

### Jeu de données

Nous fournissons les données et ressources suivantes :

- **Un génome de référence** constitué d'une petite région du chromosome 20 humain (de hg19/b37) et ses fichiers accessoires (index et dictionnaire de séquence).
- **Trois échantillons de séquençage du génome entier** correspondant à un trio familial (mère, père et fils), qui ont été réduits à une petite tranche de données sur le chromosome 20 pour maintenir la taille des fichiers petite.
  Il s'agit de données de séquençage Illumina à lectures courtes qui ont déjà été mappées sur le génome de référence, fournies au format [BAM](https://samtools.github.io/hts-specs/SAMv1.pdf) (Binary Alignment Map, une version compressée de SAM, Sequence Alignment Map).
- **Une liste d'intervalles génomiques**, c'est-à-dire des coordonnées sur le génome où nos échantillons ont des données appropriées pour appeler des variants, fournie au format BED.

### Workflow

Dans cette partie du cours, nous allons développer un workflow qui effectue les opérations suivantes :

1. Générer un fichier d'index pour chaque fichier BAM d'entrée en utilisant [Samtools](https://www.htslib.org/)
2. Exécuter GATK HaplotypeCaller sur chaque fichier BAM d'entrée pour générer des appels de variants par échantillon au format VCF (Variant Call Format)

<figure class="excalidraw">
--8<-- "docs/en/docs/nf4_science/genomics/img/hello-gatk-1.svg"
</figure>

!!! note

    Les fichiers d'index sont une caractéristique courante des formats de fichiers bioinformatiques ; ils contiennent des informations sur la structure du fichier principal qui permettent à des outils comme GATK d'accéder à un sous-ensemble des données sans avoir à lire l'intégralité du fichier.
    C'est important en raison de la taille que peuvent atteindre ces fichiers.

---

## 0. Échauffement : Tester les commandes Samtools et GATK de manière interactive

Nous voulons d'abord essayer les commandes manuellement avant de tenter de les intégrer dans un workflow.
Les outils dont nous avons besoin (Samtools et GATK) ne sont pas installés dans l'environnement GitHub Codespaces, nous allons donc les utiliser via des conteneurs (voir [Hello Containers](../../hello_nextflow/05_hello_containers.md)).

!!! note

     Assurez-vous d'être dans le répertoire `nf4-science/genomics` de sorte que la dernière partie du chemin affiché lorsque vous tapez `pwd` soit `genomics`.

### 0.1. Indexer un fichier BAM d'entrée avec Samtools

Nous allons télécharger un conteneur Samtools, le démarrer de manière interactive et exécuter la commande `samtools index` sur l'un des fichiers BAM.

#### 0.1.1. Télécharger le conteneur Samtools

```bash
docker pull community.wave.seqera.io/library/samtools:1.20--b5dfbd93de237464
```

<!--
??? success "Sortie de la commande"

    ```console

    ```
-->

#### 0.1.2. Démarrer le conteneur Samtools de manière interactive

```bash
docker run -it -v ./data:/data community.wave.seqera.io/library/samtools:1.20--b5dfbd93de237464
```

<!--
??? success "Sortie de la commande"

    ```console

    ```
-->

#### 0.1.3. Exécuter la commande d'indexation

La [documentation Samtools](https://www.htslib.org/doc/samtools-index.html) nous donne la ligne de commande à exécuter pour indexer un fichier BAM.

Nous devons seulement fournir le fichier d'entrée ; l'outil générera automatiquement un nom pour la sortie en ajoutant `.bai` au nom du fichier d'entrée.

```bash
samtools index /data/bam/reads_mother.bam
```

Cela devrait se terminer immédiatement, et vous devriez maintenant voir un fichier appelé `reads_mother.bam.bai` dans le même répertoire que le fichier BAM d'entrée original.

??? abstract "Contenu du répertoire"

    ```console
    data/bam/
    ├── reads_father.bam
    ├── reads_mother.bam
    ├── reads_mother.bam.bai
    └── reads_son.bam
    ```

#### 0.1.4. Quitter le conteneur Samtools

```bash
exit
```

### 0.2. Appeler des variants avec GATK HaplotypeCaller

Nous allons télécharger un conteneur GATK, le démarrer de manière interactive et exécuter la commande `gatk HaplotypeCaller` sur le fichier BAM que nous venons d'indexer.

#### 0.2.1. Télécharger le conteneur GATK

```bash
docker pull community.wave.seqera.io/library/gatk4:4.5.0.0--730ee8817e436867
```

<!--
??? success "Sortie de la commande"

    ```console

    ```
-->

#### 0.2.2. Démarrer le conteneur GATK de manière interactive

```bash
docker run -it -v ./data:/data community.wave.seqera.io/library/gatk4:4.5.0.0--730ee8817e436867
```

<!--
??? success "Sortie de la commande"

    ```console

    ```
-->

#### 0.2.3. Exécuter la commande d'appel de variants

La [documentation GATK](https://gatk.broadinstitute.org/hc/en-us/articles/21905025322523-HaplotypeCaller) nous donne la ligne de commande à exécuter pour effectuer l'appel de variants sur un fichier BAM.

Nous devons fournir le fichier BAM d'entrée (`-I`) ainsi que le génome de référence (`-R`), un nom pour le fichier de sortie (`-O`) et une liste d'intervalles génomiques à analyser (`-L`).

Cependant, nous n'avons pas besoin de spécifier le chemin vers le fichier d'index ; l'outil le recherchera automatiquement dans le même répertoire, en se basant sur la convention établie de nommage et de co-localisation.
Il en va de même pour les fichiers accessoires du génome de référence (fichiers d'index et de dictionnaire de séquence, `*.fai` et `*.dict`).

```bash
gatk HaplotypeCaller \
        -R /data/ref/ref.fasta \
        -I /data/bam/reads_mother.bam \
        -O reads_mother.vcf \
        -L /data/ref/intervals.bed
```

<!--
??? success "Sortie de la commande"

    ```console

    ```
-->

Le fichier de sortie `reads_mother.vcf` est créé dans votre répertoire de travail dans le conteneur, vous ne le verrez donc pas dans l'explorateur de fichiers de VS Code à moins que vous ne changiez le chemin du fichier de sortie.
Cependant, c'est un petit fichier de test, vous pouvez donc utiliser `cat` pour l'ouvrir et afficher son contenu.
Si vous faites défiler jusqu'au début du fichier, vous trouverez un en-tête composé de nombreuses lignes de métadonnées, suivi d'une liste d'appels de variants, un par ligne.

```console title="reads_mother.vcf" linenums="26"
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	reads_mother
20_10037292_10066351	3480	.	C	CT	503.03	.	AC=2;AF=1.00;AN=2;DP=23;ExcessHet=0.0000;FS=0.000;MLEAC=2;MLEAF=1.00;MQ=60.00;QD=27.95;SOR=1.179	GT:AD:DP:GQ:PL	1/1:0,18:18:54:517,54,0
20_10037292_10066351	3520	.	AT	A	609.03	.	AC=2;AF=1.00;AN=2;DP=18;ExcessHet=0.0000;FS=0.000;MLEAC=2;MLEAF=1.00;MQ=60.00;QD=33.83;SOR=0.693	GT:AD:DP:GQ:PL	1/1:0,18:18:54:623,54,0
20_10037292_10066351	3529	.	T	A	155.64	.	AC=1;AF=0.500;AN=2;BaseQRankSum=-0.544;DP=21;ExcessHet=0.0000;FS=1.871;MLEAC=1;MLEAF=0.500;MQ=60.00;MQRankSum=0.000;QD=7.78;ReadPosRankSum=-1.158;SOR=1.034	GT:AD:DP:GQ:PL	0/1:12,8:20:99:163,0,328
```

Chaque ligne décrit un variant possible identifié dans les données de séquençage de l'échantillon. Pour obtenir des conseils sur l'interprétation du format VCF, consultez [cet article utile](https://www.ebi.ac.uk/training/online/courses/human-genetic-variation-introduction/variant-identification-and-analysis/understanding-vcf-format/).

Le fichier VCF de sortie est accompagné d'un fichier d'index appelé `reads_mother.vcf.idx` qui a été automatiquement créé par GATK.
Il a la même fonction que le fichier d'index BAM, permettre aux outils de rechercher et récupérer des sous-ensembles de données sans charger l'intégralité du fichier.

#### 0.2.4. Quitter le conteneur GATK

```bash
exit
```

### À retenir

Vous savez comment tester les commandes d'indexation Samtools et d'appel de variants GATK dans leurs conteneurs respectifs.

### Et ensuite ?

Apprenez à intégrer ces mêmes commandes dans un workflow en deux étapes qui utilise des conteneurs pour exécuter le travail.

---

## 1. Écrire un workflow à une seule étape qui exécute Samtools index sur un fichier BAM

Nous vous fournissons un fichier de workflow, `genomics-1.nf`, qui décrit les parties principales du workflow.
Il n'est pas fonctionnel ; son objectif est simplement de servir de squelette que vous utiliserez pour écrire le workflow réel.

### 1.1. Définir le processus d'indexation

Commençons par écrire un processus, que nous appellerons `SAMTOOLS_INDEX`, décrivant l'opération d'indexation.

```groovy title="genomics-1.nf" linenums="9"
/*
 * Générer le fichier d'index BAM
 */
process SAMTOOLS_INDEX {

    container 'community.wave.seqera.io/library/samtools:1.20--b5dfbd93de237464'

    input:
    path input_bam

    output:
    path "${input_bam}.bai"

    script:
    """
    samtools index '$input_bam'
    """
}
```

Vous devriez reconnaître tous les éléments de ce que vous avez appris dans la Partie 1 et la Partie 2 de cette série de formation.

Ce processus va nécessiter que nous lui passions un chemin de fichier via l'entrée `input_bam`, configurons donc cela ensuite.

### 1.2. Ajouter une déclaration de paramètre d'entrée

En haut du fichier, sous la section `Pipeline parameters`, nous déclarons un paramètre CLI appelé `reads_bam` et lui donnons une valeur par défaut.
De cette façon, nous pouvons être paresseux et ne pas spécifier l'entrée lorsque nous tapons la commande pour lancer le pipeline (à des fins de développement).

```groovy title="genomics-1.nf" linenums="3"
/*
 * Paramètres du pipeline
 */
params {
    // Entrée principale
    reads_bam: Path = "${projectDir}/data/bam/reads_mother.bam"
}
```

Maintenant nous avons un processus prêt, ainsi qu'un paramètre pour lui donner une entrée sur laquelle s'exécuter, câblons donc ces éléments ensemble.

!!! note

    `${projectDir}` est une variable Nextflow intégrée qui pointe vers le répertoire où se trouve le script de workflow Nextflow actuel (`genomics-1.nf`).

    Cela facilite le référencement des fichiers, des répertoires de données et d'autres ressources inclus dans le dépôt du workflow sans coder en dur des chemins absolus.

### 1.3. Ajouter un bloc workflow pour exécuter SAMTOOLS_INDEX

Dans le bloc `workflow`, nous devons configurer un **canal** pour alimenter l'entrée du processus `SAMTOOLS_INDEX` ; ensuite nous pouvons appeler le processus lui-même pour qu'il s'exécute sur le contenu de ce canal.

```groovy title="genomics-1.nf" linenums="24"
workflow {

    main:
    // Créer le canal d'entrée (fichier unique via paramètre CLI)
    reads_ch = channel.fromPath(params.reads_bam)

    // Créer le fichier d'index pour le fichier BAM d'entrée
    SAMTOOLS_INDEX(reads_ch)

    publish:
    bam_index = SAMTOOLS_INDEX.out
}
```

Le bloc workflow a deux sections :

- `main:` contient les opérations de canal et les appels de processus
- `publish:` déclare quelles sorties doivent être publiées, en les assignant à des cibles nommées

Vous remarquerez que nous utilisons le même factory de canal `.fromPath` que nous avons utilisé dans [Hello Channels](../../hello_nextflow/02_hello_channels.md).
En effet, nous faisons quelque chose de très similaire.
La différence est que nous disons à Nextflow de simplement charger le chemin du fichier lui-même dans le canal en tant qu'élément d'entrée, plutôt que de lire son contenu.

### 1.4. Ajouter un bloc output pour définir où les résultats sont publiés

Après le bloc workflow, nous ajoutons un bloc `output` qui spécifie où publier les sorties du workflow.

```groovy title="genomics-1.nf" linenums="37"
output {
    bam_index {
        path '.'
    }
}
```

Chaque cible nommée de la section `publish:` (comme `bam_index`) obtient son propre bloc où vous pouvez configurer le chemin de sortie relatif au répertoire de sortie de base.

!!! note

    Même si les fichiers de données que nous utilisons ici sont très petits, en génomique ils peuvent devenir très volumineux.
    Par défaut, Nextflow crée des liens symboliques vers les fichiers de sortie dans le répertoire de publication, ce qui évite les copies de fichiers inutiles.
    Vous pouvez modifier ce comportement en utilisant l'option `mode` (par exemple, `mode 'copy'`) pour créer des copies réelles à la place.
    Sachez que les liens symboliques seront rompus lorsque vous nettoierez votre répertoire `work`, donc pour les workflows de production vous voudrez peut-être utiliser `mode 'copy'`.

### 1.5. Configurer le répertoire de sortie

Le répertoire de sortie de base est défini via l'option de configuration `outputDir`. Ajoutez-la à `nextflow.config` :

=== "Après"

    ```groovy title="nextflow.config" hl_lines="2"
    docker.enabled = true
    outputDir = 'results_genomics'
    ```

=== "Avant"

    ```groovy title="nextflow.config"
    docker.enabled = true
    ```

### 1.6. Exécuter le workflow pour vérifier que l'étape d'indexation fonctionne

Exécutons le workflow ! Pour rappel, nous n'avons pas besoin de spécifier une entrée dans la ligne de commande car nous avons configuré une valeur par défaut pour l'entrée lors de la déclaration du paramètre d'entrée.

```bash
nextflow run genomics-1.nf
```

??? success "Sortie de la commande"

    ```console
    N E X T F L O W   ~  version 25.10.2

    ┃ Launching `genomics-1.nf` [reverent_sinoussi] DSL2 - revision: 41d43ad7fe

    executor >  local (1)
    [2a/e69536] SAMTOOLS_INDEX (1) | 1 of 1 ✔
    ```

Vous pouvez vérifier que le fichier d'index a été généré correctement en regardant dans le répertoire de travail ou dans le répertoire des résultats.

??? abstract "Contenu du répertoire de travail"

    ```console
    work/2a/e695367b2f60df09cf826b07192dc3
    ├── reads_mother.bam -> /workspaces/training/nf4-science/genomics/data/bam/reads_mother.bam
    └── reads_mother.bam.bai
    ```

??? abstract "Contenu du répertoire des résultats"

    ```console
    results_genomics/
    └── reads_mother.bam.bai -> */87/908bba*/reads_mother.bam.bai
    ```

Le voilà !

### À retenir

Vous savez comment intégrer un outil de génomique dans un workflow Nextflow à une seule étape et le faire fonctionner en utilisant un conteneur.

### Et ensuite ?

Ajouter une deuxième étape qui consomme la sortie de la première.

---

## 2. Ajouter un deuxième processus pour exécuter GATK HaplotypeCaller sur le fichier BAM indexé

Maintenant que nous avons un index pour notre fichier d'entrée, nous pouvons passer à la configuration de l'étape d'appel de variants, qui est la partie intéressante du workflow.

### 2.1. Définir le processus d'appel de variants

Écrivons un processus, que nous appellerons `GATK_HAPLOTYPECALLER`, décrivant l'opération d'appel de variants.

```groovy title="genomics-1.nf" linenums="44"
/*
 * Appeler des variants avec GATK HaplotypeCaller
 */
process GATK_HAPLOTYPECALLER {

    container "community.wave.seqera.io/library/gatk4:4.5.0.0--730ee8817e436867"

    input:
    path input_bam
    path input_bam_index
    path ref_fasta
    path ref_index
    path ref_dict
    path interval_list

    output:
    path "${input_bam}.vcf"     , emit: vcf
    path "${input_bam}.vcf.idx" , emit: idx

    script:
    """
    gatk HaplotypeCaller \
        -R ${ref_fasta} \
        -I ${input_bam} \
        -O ${input_bam}.vcf \
        -L ${interval_list}
    """
}
```

Vous remarquerez que nous avons introduit une nouvelle syntaxe ici (`emit:`) pour nommer de manière unique chacun de nos canaux de sortie, et les raisons de cela deviendront claires bientôt.

Cette commande prend pas mal plus d'entrées, car GATK a besoin de plus d'informations pour effectuer l'analyse par rapport à un simple travail d'indexation.
Mais vous noterez qu'il y a encore plus d'entrées définies dans le bloc d'entrées qu'il n'y en a listées dans la commande GATK. Pourquoi donc ?

!!! note

    GATK sait rechercher le fichier d'index BAM et les fichiers accessoires du génome de référence car il connaît les conventions entourant ces fichiers.
    Cependant, Nextflow est conçu pour être indépendant du domaine et ne sait rien des exigences des formats de fichiers bioinformatiques.

Nous devons dire explicitement à Nextflow qu'il doit placer ces fichiers dans le répertoire de travail au moment de l'exécution ; sinon il ne le fera pas, et GATK lancera (à juste titre) une erreur concernant les fichiers d'index manquants.

De même, nous devons lister explicitement le fichier d'index du VCF de sortie (le fichier `"${input_bam}.vcf.idx"`) pour que Nextflow sache garder la trace de ce fichier au cas où il serait nécessaire dans les étapes suivantes.

### 2.2. Ajouter des définitions pour les entrées accessoires

Puisque notre nouveau processus attend qu'un certain nombre de fichiers supplémentaires soient fournis, nous configurons des paramètres CLI pour eux sous la section `Pipeline parameters`, ainsi que quelques valeurs par défaut (pour les mêmes raisons qu'avant).

```groovy title="genomics-1.nf" linenums="8"
    // Fichiers accessoires
    reference: Path = "${projectDir}/data/ref/ref.fasta"
    reference_index: Path = "${projectDir}/data/ref/ref.fasta.fai"
    reference_dict: Path = "${projectDir}/data/ref/ref.dict"
    intervals: Path = "${projectDir}/data/ref/intervals.bed"
```

### 2.3. Créer des variables pour contenir les chemins des fichiers accessoires

Bien que les entrées de données principales soient diffusées dynamiquement à travers les canaux, il existe deux approches pour gérer les fichiers accessoires. L'approche recommandée consiste à créer des canaux explicites, ce qui rend le flux de données plus clair et plus cohérent. Alternativement, la fonction file() pour créer des variables peut être utilisée pour les cas plus simples, en particulier lorsque vous devez référencer le même fichier dans plusieurs processus - bien que sachez que cela crée toujours des canaux implicitement. <!-- TODO: Clarify: is this still necessary with typed inputs? -->

Ajoutez ceci au bloc workflow (après la création de `reads_ch`, à l'intérieur de la section `main:`) :

```groovy title="genomics-1.nf" linenums="79"
    // Charger les chemins de fichiers pour les fichiers accessoires (référence et intervalles)
    ref_file        = file(params.reference)
    ref_index_file  = file(params.reference_index)
    ref_dict_file   = file(params.reference_dict)
    intervals_file  = file(params.intervals)
```

Cela rendra les chemins des fichiers accessoires disponibles pour être fournis en entrée à tous les processus qui en ont besoin.

### 2.4. Ajouter un appel au bloc workflow pour exécuter GATK_HAPLOTYPECALLER

Maintenant que nous avons configuré notre deuxième processus et que toutes les entrées et fichiers accessoires sont prêts et disponibles, nous pouvons ajouter un appel au processus `GATK_HAPLOTYPECALLER` dans le corps du workflow.

```groovy title="genomics-1.nf" linenums="88"
    // Appeler les variants depuis le fichier BAM indexé
    GATK_HAPLOTYPECALLER(
        reads_ch,
        SAMTOOLS_INDEX.out,
        ref_file,
        ref_index_file,
        ref_dict_file,
        intervals_file
    )
```

Vous devriez reconnaître la syntaxe `*.out` de la Partie 1 de cette série de formation ; nous disons à Nextflow de prendre le canal de sortie de `SAMTOOLS_INDEX` et de le brancher dans l'appel du processus `GATK_HAPLOTYPECALLER`.

!!! note

    Vous remarquerez que les entrées sont fournies dans le même ordre exact dans l'appel au processus que celui dans lequel elles sont listées dans le bloc d'entrée du processus.
    Dans Nextflow, les entrées sont positionnelles, ce qui signifie que vous _devez_ suivre le même ordre ; et bien sûr il doit y avoir le même nombre d'éléments.

### 2.5. Mettre à jour la section publish et le bloc output

Nous devons mettre à jour la section `publish:` pour inclure les sorties VCF, et ajouter les cibles correspondantes dans le bloc `output`.

```groovy title="genomics-1.nf" linenums="99"
    publish:
    bam_index = SAMTOOLS_INDEX.out
    vcf = GATK_HAPLOTYPECALLER.out.vcf
    vcf_idx = GATK_HAPLOTYPECALLER.out.idx
}

output {
    bam_index {
        path '.'
    }
    vcf {
        path '.'
    }
    vcf_idx {
        path '.'
    }
}
```

### 2.6. Exécuter le workflow pour vérifier que l'étape d'appel de variants fonctionne

Exécutons le workflow étendu avec `-resume` pour ne pas avoir à exécuter à nouveau l'étape d'indexation.

```bash
nextflow run genomics-1.nf -resume
```

??? success "Sortie de la commande"

    ```console
    N E X T F L O W   ~  version 25.10.2

    ┃ Launching `genomics-1.nf` [grave_volta] DSL2 - revision: 4790abc96a

    executor >  local (1)
    [2a/e69536] SAMTOOLS_INDEX (1)       | 1 of 1, cached: 1 ✔
    [53/e18e98] GATK_HAPLOTYPECALLER (1) | 1 of 1 ✔
    ```

Maintenant si nous regardons la sortie de la console, nous voyons les deux processus listés.

Le premier processus a été ignoré grâce à la mise en cache, comme prévu, tandis que le deuxième processus a été exécuté puisqu'il est tout nouveau.

Vous trouverez les fichiers de sortie dans le répertoire des résultats (sous forme de liens symboliques vers le répertoire de travail).

??? abstract "Contenu du répertoire"

    ```console
    results_genomics/
    ├── reads_mother.bam.bai -> */87/908bba*/reads_mother.bam.bai
    ├── reads_mother.bam.vcf -> */cf/36f756*/reads_mother.bam.vcf
    └── reads_mother.bam.vcf.idx -> */cf/36f756*/reads_mother.bam.vcf.idx
    ```

Si vous ouvrez le fichier VCF, vous devriez voir le même contenu que dans le fichier que vous avez généré en exécutant la commande GATK directement dans le conteneur.

```console title="reads_mother.bam.vcf" linenums="26"
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	reads_mother
20_10037292_10066351	3480	.	C	CT	503.03	.	AC=2;AF=1.00;AN=2;DP=23;ExcessHet=0.0000;FS=0.000;MLEAC=2;MLEAF=1.00;MQ=60.00;QD=27.95;SOR=1.179	GT:AD:DP:GQ:PL	1/1:0,18:18:54:517,54,0
20_10037292_10066351	3520	.	AT	A	609.03	.	AC=2;AF=1.00;AN=2;DP=18;ExcessHet=0.0000;FS=0.000;MLEAC=2;MLEAF=1.00;MQ=60.00;QD=33.83;SOR=0.693	GT:AD:DP:GQ:PL	1/1:0,18:18:54:623,54,0
20_10037292_10066351	3529	.	T	A	155.64	.	AC=1;AF=0.500;AN=2;BaseQRankSum=-0.544;DP=21;ExcessHet=0.0000;FS=1.871;MLEAC=1;MLEAF=0.500;MQ=60.00;MQRankSum=0.000;QD=7.78;ReadPosRankSum=-1.158;SOR=1.034	GT:AD:DP:GQ:PL	0/1:12,8:20:99:163,0,328
```

C'est la sortie que nous voulons générer pour chaque échantillon dans notre étude.

### À retenir

Vous savez comment créer un workflow très basique à deux étapes qui effectue un véritable travail d'analyse et est capable de gérer les idiosyncrasies des formats de fichiers génomiques comme les fichiers accessoires.

### Et ensuite ?

Faire en sorte que le workflow gère plusieurs échantillons en masse.

---

## 3. Adapter le workflow pour qu'il s'exécute sur un lot d'échantillons

C'est bien d'avoir un workflow qui peut automatiser le traitement d'un seul échantillon, mais que se passe-t-il si vous avez 1000 échantillons ?
Devez-vous écrire un script bash qui boucle sur tous vos échantillons ?

Non, Dieu merci ! Il suffit de faire une petite modification du code et Nextflow gérera cela pour vous aussi.

### 3.1. Transformer la déclaration du paramètre d'entrée en un tableau listant les trois échantillons

Transformons ce chemin de fichier par défaut dans la déclaration du fichier BAM d'entrée en un tableau listant les chemins de fichiers pour nos trois échantillons de test, sous la section `Pipeline parameters`.

=== "Après"

    ```groovy title="genomics-1.nf" linenums="7"
    // Entrée principale (tableau de trois échantillons)
        reads_bam = [
            "${projectDir}/data/bam/reads_mother.bam",
            "${projectDir}/data/bam/reads_father.bam",
            "${projectDir}/data/bam/reads_son.bam"
        ]
    ```

=== "Avant"

    ```groovy title="genomics-1.nf" linenums="7"
        // Entrée principale
        reads_bam: Path = "${projectDir}/data/bam/reads_mother.bam"
    ```

!!! note

    Lors de l'utilisation de déclarations de paramètres typés (comme `reads_bam: Path`), vous ne pouvez pas assigner une valeur de tableau.
    Pour les tableaux, omettez l'annotation de type.

Et c'est en fait tout ce que nous devons faire, car le factory de canal que nous utilisons dans le corps du workflow (`.fromPath`) est tout aussi heureux d'accepter plusieurs chemins de fichiers à charger dans le canal d'entrée qu'il l'était d'en charger un seul.

!!! note

    Normalement, vous ne voudriez pas coder en dur la liste des échantillons dans votre fichier de workflow, mais nous le faisons ici pour garder les choses simples.
    Nous présenterons des façons plus élégantes de gérer les entrées plus tard dans cette série de formation.

### 3.2. Exécuter le workflow pour vérifier qu'il s'exécute sur les trois échantillons

Essayons maintenant d'exécuter le workflow maintenant que la plomberie est configurée pour s'exécuter sur les trois échantillons de test.

```bash
nextflow run genomics-1.nf -resume
```

Chose amusante : cela _pourrait fonctionner_, OU cela _pourrait échouer_. Par exemple, voici une exécution qui a réussi :

??? success "Sortie de la commande"

    ```console
    N E X T F L O W   ~  version 25.10.2

    ┃ Launching `genomics-1.nf` [peaceful_yalow] DSL2 - revision: a256d113ad

    executor >  local (6)
    [4f/7071b0] SAMTOOLS_INDEX (3)       | 3 of 3, cached: 1 ✔
    [7a/89bc43] GATK_HAPLOTYPECALLER (2) | 3 of 3, cached: 1 ✔
    ```

Si votre exécution de workflow a réussi, exécutez-la à nouveau jusqu'à ce que vous obteniez une erreur comme celle-ci :

??? failure "Sortie de la commande"

    ```console
    N E X T F L O W   ~  version 25.10.2

    ┃ Launching `genomics-1.nf` [loving_pasteur] DSL2 - revision: d2a8e63076

    executor >  local (4)
    [01/eea165] SAMTOOLS_INDEX (2)       | 3 of 3, cached: 1 ✔
    [a5/fa9fd0] GATK_HAPLOTYPECALLER (3) | 1 of 3, cached: 1
    ERROR ~ Error executing process > 'GATK_HAPLOTYPECALLER (2)'

    Caused by:
      Process `GATK_HAPLOTYPECALLER (2)` terminated with an error exit status (2)

    Command executed:

      gatk HaplotypeCaller         -R ref.fasta         -I reads_father.bam         -O reads_father.bam.vcf         -L intervals.bed

    Command exit status:
      2

    Command error:
      ...
      A USER ERROR has occurred: Traversal by intervals was requested but some input files are not indexed.
      ...
    ```

Si vous regardez la sortie d'erreur de la commande GATK, il y aura une ligne comme celle-ci :

```console
A USER ERROR has occurred: Traversal by intervals was requested but some input files are not indexed.
```

Eh bien, c'est étrange, considérant que nous avons explicitement indexé les fichiers BAM dans la première étape du workflow. Pourrait-il y avoir un problème avec la plomberie ?

#### 3.2.1. Vérifier les répertoires de travail pour les appels pertinents

Jetons un coup d'œil à l'intérieur du répertoire de travail pour l'appel de processus `GATK_HAPLOTYPECALLER` échoué listé dans la sortie de la console.

??? abstract "Contenu du répertoire"

    ```console
    work/a5/fa9fd0994b6beede5fb9ea073596c2
    ├── intervals.bed -> /workspaces/training/nf4-science/genomics/data/ref/intervals.bed
    ├── reads_father.bam.bai -> /workspaces/training/nf4-science/genomics/work/01/eea16597bd6e810fb4cf89e60f8c2d/reads_father.bam.bai
    ├── reads_son.bam -> /workspaces/training/nf4-science/genomics/data/bam/reads_son.bam
    ├── reads_son.bam.vcf
    ├── reads_son.bam.vcf.idx
    ├── ref.dict -> /workspaces/training/nf4-science/genomics/data/ref/ref.dict
    ├── ref.fasta -> /workspaces/training/nf4-science/genomics/data/ref/ref.fasta
    └── ref.fasta.fai -> /workspaces/training/nf4-science/genomics/data/ref/ref.fasta.fai
    ```

Portez une attention particulière aux noms du fichier BAM et de l'index BAM qui sont listés dans ce répertoire : `reads_son.bam` et `reads_father.bam.bai`.

Qu'est-ce que c'est que ça ? Nextflow a placé un fichier d'index dans le répertoire de travail de cet appel de processus, mais c'est le mauvais. Comment cela a-t-il pu se produire ?

#### 3.2.2. Utiliser l'opérateur [view()](https://www.nextflow.io/docs/latest/reference/operator.html#view) pour inspecter le contenu des canaux

Ajoutez ces deux lignes dans le corps du workflow avant l'appel du processus `GATK_HAPLOTYPER` :

```groovy title="genomics-1.nf" linenums="84"
    // diagnostics temporaires
    reads_ch.view()
    SAMTOOLS_INDEX.out.view()
```

Puis exécutez à nouveau la commande du workflow.

```bash
nextflow run genomics-1.nf
```

Une fois encore, cela peut réussir ou échouer. Voici une exécution réussie :

??? success "Sortie de la commande"

    ```console
    N E X T F L O W   ~  version 25.10.2

    ┃ Launching `genomics-1.nf` [fervent_pasteur] DSL2 - revision: a256d113ad

    /workspaces/training/nf4-science/genomics/data/bam/reads_mother.bam
    /workspaces/training/nf4-science/genomics/data/bam/reads_father.bam
    /workspaces/training/nf4-science/genomics/data/bam/reads_son.bam
    executor >  local (6)
    [4f/7071b0] SAMTOOLS_INDEX (3)       | 3 of 3 ✔
    /workspaces/training/nf4-science/genomics/work/b4/45a376f0e724be1dc626a6807f73d8/reads_mother.bam.bai
    /workspaces/training/nf4-science/genomics/work/4f/7071b082b45dd85b1c9b6b3b32cb69/reads_father.bam.bai
    /workspaces/training/nf4-science/genomics/work/3c/331645a9e20e67edae10da5ba17c7b/reads_son.bam.bai
    [a2/dbd8d5] GATK_HAPLOTYPECALLER (3) | 3 of 3 ✔
    ```

Et voici une échouée :

??? failure "Sortie de la commande"

    ```console
    N E X T F L O W   ~  version 25.10.2

    ┃ Launching `genomics-1.nf` [angry_hamilton] DSL2 - revision: a256d113ad

    /workspaces/training/nf4-science/genomics/data/bam/reads_mother.bam
    /workspaces/training/nf4-science/genomics/data/bam/reads_father.bam
    /workspaces/training/nf4-science/genomics/data/bam/reads_son.bam
    executor >  local (6)
    [4f/7071b0] SAMTOOLS_INDEX (3)       | 3 of 3 ✔
    /workspaces/training/nf4-science/genomics/work/4f/7071b082b45dd85b1c9b6b3b32cb69/reads_father.bam.bai
    /workspaces/training/nf4-science/genomics/work/3c/331645a9e20e67edae10da5ba17c7b/reads_son.bam.bai
    /workspaces/training/nf4-science/genomics/work/b4/45a376f0e724be1dc626a6807f73d8/reads_mother.bam.bai
    [a3/cf3a89] GATK_HAPLOTYPECALLER (3) | 1 of 3
    ERROR ~ Error executing process > 'GATK_HAPLOTYPECALLER (3)'
    ...
    ```

Vous devrez peut-être l'exécuter plusieurs fois pour qu'elle échoue à nouveau.
Cette erreur ne se reproduira pas systématiquement car elle dépend d'une certaine variabilité dans les temps d'exécution des appels de processus individuels.

Voici à quoi ressemble la sortie des deux appels `.view()` que nous avons ajoutés pour une exécution échouée :

```console
/workspaces/training/nf4-science/genomics/data/bam/reads_mother.bam
/workspaces/training/nf4-science/genomics/data/bam/reads_father.bam
/workspaces/training/nf4-science/genomics/data/bam/reads_son.bam
/workspaces/training/nf4-science/genomics/work/9c/53492e3518447b75363e1cd951be4b/reads_father.bam.bai
/workspaces/training/nf4-science/genomics/work/cc/37894fffdf6cc84c3b0b47f9b536b7/reads_son.bam.bai
/workspaces/training/nf4-science/genomics/work/4d/dff681a3d137ba7d9866e3d9307bd0/reads_mother.bam.bai
```

Les trois premières lignes correspondent au canal d'entrée et les secondes, au canal de sortie.
Vous pouvez voir que les fichiers BAM et les fichiers d'index pour les trois échantillons ne sont pas listés dans le même ordre !

!!! note

    Lorsque vous appelez un processus Nextflow sur un canal contenant plusieurs éléments, Nextflow essaiera de paralléliser l'exécution autant que possible, et collectera les sorties dans l'ordre dans lequel elles deviennent disponibles.
    La conséquence est que les sorties correspondantes peuvent être collectées dans un ordre différent de celui dans lequel les entrées originales ont été fournies.

Tel qu'il est écrit actuellement, notre script de workflow suppose que les fichiers d'index sortiront de l'étape d'indexation listés dans le même ordre mère/père/fils que les entrées ont été données.
Mais ce n'est pas garanti d'être le cas, c'est pourquoi parfois (mais pas toujours) les mauvais fichiers sont appariés dans la deuxième étape.

Pour corriger cela, nous devons nous assurer que les fichiers BAM et leurs fichiers d'index voyagent ensemble à travers les canaux.

!!! tip

    Les instructions `view()` dans le code du workflow ne font rien, donc ce n'est pas un problème de les laisser.
    Cependant, elles encombrent votre sortie de console, nous recommandons donc de les supprimer lorsque vous avez fini de déboguer le problème.

### 3.3. Changer la sortie du processus SAMTOOLS_INDEX en un tuple qui garde le fichier d'entrée et son index ensemble

La façon la plus simple d'assurer qu'un fichier BAM et son index restent étroitement associés est de les emballer ensemble dans un tuple sortant de la tâche d'index.

!!! note

    Un **tuple** est une liste ordonnée finie d'éléments qui est couramment utilisée pour retourner plusieurs valeurs d'une fonction. Les tuples sont particulièrement utiles pour passer plusieurs entrées ou sorties entre des processus tout en préservant leur association et leur ordre.

Tout d'abord, changeons la sortie du processus `SAMTOOLS_INDEX` pour inclure le fichier BAM dans sa déclaration de sortie.

=== "Après"

    ```groovy title="genomics-1.nf" linenums="20"
        output:
        tuple path(input_bam), path("${input_bam}.bai")
    ```

=== "Avant"

    ```groovy title="genomics-1.nf" linenums="20"
        output:
        path "${input_bam}.bai"
    ```

De cette façon, chaque fichier d'index sera étroitement couplé avec son fichier BAM original, et la sortie globale de l'étape d'indexation sera un seul canal contenant des paires de fichiers.

### 3.4. Changer l'entrée du processus GATK_HAPLOTYPECALLER en un tuple

Puisque nous avons changé la « forme » de la sortie du premier processus dans le workflow, nous devons mettre à jour la définition d'entrée du deuxième processus pour correspondre.

Spécifiquement, là où nous déclarions auparavant deux chemins d'entrée séparés dans le bloc d'entrée du processus `GATK_HAPLOTYPECALLER`, nous déclarons maintenant une seule entrée correspondant à la structure du tuple émis par `SAMTOOLS_INDEX`.

=== "Après"

    ```groovy title="genomics-1.nf" linenums="51"
        input:
        tuple path(input_bam), path(input_bam_index)
    ```

=== "Avant"

    ```groovy title="genomics-1.nf" linenums="51"
        input:
        path input_bam
        path input_bam_index
    ```

Bien sûr, puisque nous avons maintenant changé la forme des entrées que `GATK_HAPLOTYPECALLER` attend, nous devons mettre à jour l'appel du processus en conséquence dans le corps du workflow.

### 3.5. Mettre à jour l'appel à GATK_HAPLOTYPECALLER dans le bloc workflow

Nous n'avons plus besoin de fournir le `reads_ch` original au processus `GATK_HAPLOTYPECALLER`, puisque le fichier BAM est maintenant empaqueté dans le canal de sortie par `SAMTOOLS_INDEX`.

En conséquence, nous pouvons simplement supprimer cette ligne.

=== "Après"

    ```groovy title="genomics-1.nf" linenums="88"
    GATK_HAPLOTYPECALLER(
        SAMTOOLS_INDEX.out,
    ```

=== "Avant"

    ```groovy title="genomics-1.nf" linenums="88"
    GATK_HAPLOTYPECALLER(
        reads_ch,
        SAMTOOLS_INDEX.out,
    ```

C'est tout le recâblage nécessaire pour résoudre le problème de non-correspondance des index.

### 3.6. Mettre à jour la section publish et le bloc output pour le tuple

Puisque `SAMTOOLS_INDEX.out` est maintenant un tuple contenant à la fois le BAM et son index, les deux fichiers seront publiés ensemble.
Nous renommons la cible de `bam_index` à `indexed_bam` pour refléter qu'elle contient maintenant les deux fichiers.

=== "Après"

    ```groovy title="genomics-1.nf" hl_lines="2"
        publish:
        indexed_bam = SAMTOOLS_INDEX.out
    ```

=== "Avant"

    ```groovy title="genomics-1.nf"
        publish:
        bam_index = SAMTOOLS_INDEX.out
    ```

Nous devons également mettre à jour le bloc output pour utiliser le nouveau nom de cible :

=== "Après"

    ```groovy title="genomics-1.nf" hl_lines="2"
    output {
        indexed_bam {
            path '.'
        }
    ```

=== "Avant"

    ```groovy title="genomics-1.nf"
    output {
        bam_index {
            path '.'
        }
    ```

### 3.7. Exécuter le workflow pour vérifier qu'il fonctionne correctement sur les trois échantillons à chaque fois

Bien sûr, la preuve est dans le pudding, exécutons donc le workflow à nouveau quelques fois pour nous assurer que cela fonctionnera de manière fiable à l'avenir.

```bash
nextflow run genomics-1.nf
```

Cette fois (et à chaque fois) tout devrait fonctionner correctement :

??? success "Sortie de la commande"

    ```console
    N E X T F L O W   ~  version 25.10.2

    ┃ Launching `genomics-1.nf` [special_goldstine] DSL2 - revision: 4cbbf6ea3e

    executor >  local (6)
    [d6/10c2c4] SAMTOOLS_INDEX (1)       | 3 of 3 ✔
    [88/1783aa] GATK_HAPLOTYPECALLER (2) | 3 of 3 ✔
    ```

Le répertoire des résultats contient maintenant à la fois les fichiers BAM et BAI pour chaque échantillon (du tuple), ainsi que les sorties VCF :

??? abstract "Contenu du répertoire des résultats"

    ```console
    results_genomics/
    ├── reads_father.bam -> */60/e2614c*/reads_father.bam
    ├── reads_father.bam.bai -> */60/e2614c*/reads_father.bam.bai
    ├── reads_father.bam.vcf -> */b8/91b3c8*/reads_father.bam.vcf
    ├── reads_father.bam.vcf.idx -> */b8/91b3c8*/reads_father.bam.vcf.idx
    ├── reads_mother.bam -> */3e/fededc*/reads_mother.bam
    ├── reads_mother.bam.bai -> */3e/fededc*/reads_mother.bam.bai
    ├── reads_mother.bam.vcf -> */32/5ca037*/reads_mother.bam.vcf
    ├── reads_mother.bam.vcf.idx -> */32/5ca037*/reads_mother.bam.vcf.idx
    ├── reads_son.bam -> */3c/36d1c2*/reads_son.bam
    ├── reads_son.bam.bai -> */3c/36d1c2*/reads_son.bam.bai
    ├── reads_son.bam.vcf -> */d7/a6b046*/reads_son.bam.vcf
    └── reads_son.bam.vcf.idx -> */d7/a6b046*/reads_son.bam.vcf.idx
    ```

Si vous le souhaitez, vous pouvez utiliser `.view()` à nouveau pour jeter un coup d'œil au contenu du canal de sortie de `SAMTOOLS_INDEX` :

```groovy title="genomics-1.nf" linenums="92"
SAMTOOLS_INDEX.out.view()
```

Vous verrez que le canal contient les trois tuples attendus (chemins de fichiers tronqués pour la lisibilité).

```console title="Sortie"
[*/60/e2614c*/reads_father.bam, */60/e2614c*/reads_father.bam.bai]
[*/3e/fededc*/reads_mother.bam, */3e/fededc*/reads_mother.bam.bai]
[*/3c/36d1c2*/reads_son.bam, */3c/36d1c2*/reads_son.bam.bai]
```

Ce sera beaucoup plus sûr, à l'avenir.

### À retenir

Vous savez comment faire en sorte que votre workflow s'exécute sur plusieurs échantillons (de manière indépendante).

### Et ensuite ?

Faciliter la gestion des échantillons en masse.

---

## 4. Faire en sorte que le workflow accepte un fichier texte contenant un lot de fichiers d'entrée

Une façon très courante de fournir plusieurs fichiers de données d'entrée à un workflow est de le faire avec un fichier texte contenant les chemins de fichiers.
Il peut être aussi simple qu'un fichier texte listant un chemin de fichier par ligne et rien d'autre, ou le fichier peut contenir des métadonnées supplémentaires, auquel cas il est souvent appelé samplesheet.

Ici, nous allons vous montrer comment faire le cas simple.

### 4.1. Examiner le fichier texte fourni listant les chemins de fichiers d'entrée

Nous avons déjà créé un fichier texte listant les chemins de fichiers d'entrée, appelé `sample_bams.txt`, que vous pouvez trouver dans le répertoire `data/`.

```txt title="sample_bams.txt"
/workspaces/training/nf4-science/genomics/data/bam/reads_mother.bam
/workspaces/training/nf4-science/genomics/data/bam/reads_father.bam
/workspaces/training/nf4-science/genomics/data/bam/reads_son.bam
```

Comme vous pouvez le voir, nous avons listé un chemin de fichier par ligne, et ce sont des chemins absolus.

!!! note

    Les fichiers que nous utilisons ici sont simplement sur le système de fichiers local de votre GitHub Codespaces, mais nous pourrions également pointer vers des fichiers dans le stockage cloud.

### 4.2. Mettre à jour la valeur par défaut du paramètre

Changeons la valeur par défaut de notre paramètre d'entrée `reads_bam` pour pointer vers le fichier `sample_bams.txt`.

=== "Après"

    ```groovy title="genomics-1.nf" linenums="7"
        // Entrée principale (fichier de fichiers d'entrée, un par ligne)
        reads_bam: Path = "${projectDir}/data/sample_bams.txt"
    ```

=== "Avant"

    ```groovy title="genomics-1.nf" linenums="7"
    // Entrée principale (tableau de trois échantillons)
        reads_bam = [
            "${projectDir}/data/bam/reads_mother.bam",
            "${projectDir}/data/bam/reads_father.bam",
            "${projectDir}/data/bam/reads_son.bam"
        ]
    ```

De cette façon, nous pouvons continuer à être paresseux, mais la liste des fichiers ne vit plus dans le code du workflow lui-même, ce qui est un grand pas dans la bonne direction.

### 4.3. Mettre à jour le factory de canal pour lire les lignes d'un fichier

Actuellement, notre factory de canal d'entrée traite tous les fichiers que nous lui donnons comme les entrées de données que nous voulons alimenter au processus d'indexation.
Puisque nous lui donnons maintenant un fichier qui liste les chemins de fichiers d'entrée, nous devons changer son comportement pour analyser le fichier et traiter les chemins de fichiers qu'il contient comme les entrées de données.

Heureusement, nous pouvons faire cela très simplement, juste en ajoutant l'opérateur [`.splitText()`](https://www.nextflow.io/docs/latest/reference/operator.html#operator-splittext) à l'étape de construction du canal.

=== "Après"

    ```groovy title="genomics-1.nf" linenums="68"
        // Créer le canal d'entrée à partir d'un fichier texte listant les chemins de fichiers d'entrée
        reads_ch = channel.fromPath(params.reads_bam).splitText()
    ```

=== "Avant"

    ```groovy title="genomics-1.nf" linenums="68"
        // Créer le canal d'entrée (fichier unique via paramètre CLI)
        reads_ch = channel.fromPath(params.reads_bam)
    ```

!!! tip

    C'est une autre excellente opportunité d'utiliser l'opérateur `.view()` pour regarder à quoi ressemble le contenu du canal avant et après l'application d'un opérateur.

### 4.4. Exécuter le workflow pour vérifier qu'il fonctionne correctement

Exécutons le workflow une dernière fois. Cela devrait produire le même résultat qu'avant, n'est-ce pas ?

```bash
nextflow run genomics-1.nf -resume
```

??? success "Sortie de la commande"

    ```console
    N E X T F L O W   ~  version 25.10.2

    ┃ Launching `genomics-1.nf` [sick_albattani] DSL2 - revision: 46d84642f6

    [18/23b4bb] SAMTOOLS_INDEX (1)       | 3 of 3, cached: 3 ✔
    [12/f727bb] GATK_HAPLOTYPECALLER (3) | 3 of 3, cached: 3 ✔
    ```

Oui ! En fait, Nextflow détecte correctement que les appels de processus sont exactement les mêmes, et ne se donne même pas la peine de tout ré-exécuter, puisque nous exécutions avec `-resume`.

Et c'est tout ! Notre simple workflow d'appel de variants a toutes les fonctionnalités de base que nous voulions.

### À retenir

Vous savez comment créer un workflow linéaire multi-étapes pour indexer un fichier BAM et appliquer l'appel de variants par échantillon en utilisant GATK.

Plus généralement, vous avez appris à utiliser des composants et une logique Nextflow essentiels pour construire un pipeline de génomique simple qui fait du véritable travail, en tenant compte des idiosyncrasies des formats de fichiers génomiques et des exigences des outils.

### Et ensuite ?

Célébrez votre succès et prenez une pause extra longue !

Dans la prochaine partie de ce cours, vous apprendrez à utiliser quelques fonctionnalités Nextflow supplémentaires (y compris plus d'opérateurs de canal) pour appliquer l'appel de variants conjoint aux données.
