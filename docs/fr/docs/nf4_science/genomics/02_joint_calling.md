# Partie 2 : Appel conjoint de variants sur une cohorte

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Traduction assistée par IA - [en savoir plus et suggérer des améliorations](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

Dans la première partie de ce cours, vous avez construit un pipeline d'appel de variants complètement linéaire qui traitait les données de chaque échantillon indépendamment des autres.
Cependant, dans un cas d'usage génomique réel, vous aurez généralement besoin d'examiner les appels de variants de plusieurs échantillons ensemble.

Dans cette deuxième partie, nous vous montrons comment utiliser les canaux et les opérateurs de canaux pour implémenter l'appel conjoint de variants avec GATK, en nous basant sur le pipeline de la Partie 1.

### Aperçu de la méthode

La méthode d'appel de variants GATK que nous avons utilisée dans la première partie de ce cours générait simplement des appels de variants par échantillon.
C'est bien si vous voulez seulement examiner les variants de chaque échantillon isolément, mais cela fournit des informations limitées.
Il est souvent plus intéressant d'examiner comment les appels de variants diffèrent entre plusieurs échantillons, et pour ce faire, GATK offre une méthode alternative appelée appel conjoint de variants, que nous démontrons ici.

L'appel conjoint de variants consiste à générer un type spécial de sortie de variants appelée GVCF (pour Genomic VCF) pour chaque échantillon, puis à combiner les données GVCF de tous les échantillons et enfin, à exécuter une analyse statistique de « génotypage conjoint ».

![Analyse conjointe](img/joint-calling.png)

Ce qui est spécial avec le GVCF d'un échantillon, c'est qu'il contient des enregistrements résumant les statistiques de données de séquence sur toutes les positions dans la zone ciblée du génome, et pas seulement les positions où le programme a trouvé des preuves de variation.
Ceci est essentiel pour le calcul du génotypage conjoint ([lecture complémentaire](https://gatk.broadinstitute.org/hc/en-us/articles/360035890431-The-logic-of-joint-calling-for-germline-short-variants)).

Le GVCF est produit par GATK HaplotypeCaller, le même outil que nous avons utilisé dans la Partie 1, avec un paramètre supplémentaire (`-ERC GVCF`).
La combinaison des GVCFs est effectuée avec GATK GenomicsDBImport, qui combine les appels par échantillon dans un magasin de données (analogue à une base de données), puis l'analyse de « génotypage conjoint » proprement dite est effectuée avec GATK GenotypeGVCFs.

### Workflow

Donc, pour récapituler, dans cette partie du cours, nous allons développer un workflow qui fait ce qui suit :

<figure class="excalidraw">
--8<-- "docs/en/docs/nf4_science/genomics/img/hello-gatk-2.svg"
</figure>

1. Générer un fichier d'index pour chaque fichier BAM d'entrée en utilisant Samtools
2. Exécuter GATK HaplotypeCaller sur chaque fichier BAM d'entrée pour générer un GVCF d'appels de variants génomiques par échantillon
3. Collecter tous les GVCFs et les combiner dans un magasin de données GenomicsDB
4. Exécuter le génotypage conjoint sur le magasin de données GVCF combiné pour produire un VCF au niveau de la cohorte

Nous appliquerons ceci au même jeu de données que dans la Partie 1.

---

## 0. Échauffement : Exécuter Samtools et GATK directement

Tout comme précédemment, nous voulons essayer les commandes manuellement avant de tenter de les encapsuler dans un workflow.

!!! note

     Assurez-vous d'être dans le bon répertoire de travail :
     `cd /workspaces/training/nf4-science/genomics`

### 0.1. Indexer un fichier BAM d'entrée avec Samtools

Cette première étape est la même que dans la Partie 1, donc elle devrait être très familière, mais cette fois nous devons le faire pour les trois échantillons.

!!! note

    Nous avons techniquement déjà généré des fichiers d'index pour les trois échantillons via notre pipeline, nous pourrions donc aller les récupérer dans le répertoire de résultats. Cependant, il est plus propre de simplement refaire cela manuellement, et cela ne prendra qu'une minute.

#### 0.1.1. Démarrer le conteneur Samtools en mode interactif

```bash
docker run -it -v ./data:/data community.wave.seqera.io/library/samtools:1.20--b5dfbd93de237464
```

<!--
??? success "Sortie de la commande"

    ```console

    ```
-->

#### 0.1.2. Exécuter la commande d'indexation pour les trois échantillons

```bash
samtools index /data/bam/reads_mother.bam
samtools index /data/bam/reads_father.bam
samtools index /data/bam/reads_son.bam
```

Tout comme précédemment, cela devrait produire les fichiers d'index dans le même répertoire que les fichiers BAM correspondants.

??? abstract "Contenu du répertoire"

    ```console
    data/bam/
    ├── reads_father.bam
    ├── reads_father.bam.bai
    ├── reads_mother.bam
    ├── reads_mother.bam.bai
    ├── reads_son.bam
    └── reads_son.bam.bai
    ```

Maintenant que nous avons des fichiers d'index pour les trois échantillons, nous pouvons procéder à la génération des GVCFs pour chacun d'eux.

#### 0.1.3. Quitter le conteneur Samtools

```bash
exit
```

### 0.2. Appeler les variants avec GATK HaplotypeCaller en mode GVCF

Cette deuxième étape est très similaire à ce que nous avons fait dans la Partie 1 : Hello Genomics, mais nous allons maintenant exécuter GATK en « mode GVCF ».

#### 0.2.1. Démarrer le conteneur GATK en mode interactif

```bash
docker run -it -v ./data:/data community.wave.seqera.io/library/gatk4:4.5.0.0--730ee8817e436867
```

<!--
??? success "Sortie de la commande"

    ```console

    ```
-->

#### 0.2.2. Exécuter la commande d'appel de variants avec l'option GVCF

Afin de produire un VCF génomique (GVCF), nous ajoutons l'option `-ERC GVCF` à la commande de base, ce qui active le mode GVCF de HaplotypeCaller.

Nous modifions également l'extension du fichier de sortie de `.vcf` à `.g.vcf`.
Ce n'est techniquement pas une exigence, mais c'est une convention fortement recommandée.

```bash
gatk HaplotypeCaller \
        -R /data/ref/ref.fasta \
        -I /data/bam/reads_mother.bam \
        -O reads_mother.g.vcf \
        -L /data/ref/intervals.bed \
        -ERC GVCF
```

<!--
??? success "Sortie de la commande"

    ```console

    ```
-->

Ceci crée le fichier de sortie GVCF `reads_mother.g.vcf` dans le répertoire de travail actuel dans le conteneur.

Si vous utilisez `cat` pour visualiser le contenu, vous verrez qu'il est beaucoup plus long que le VCF équivalent que nous avons généré dans la Partie 1. Vous ne pouvez même pas faire défiler jusqu'au début du fichier, et la plupart des lignes semblent assez différentes de ce que nous avons vu dans le VCF de la Partie 1.

```console title="Sortie" linenums="1674"
20_10037292_10066351    14714   .       T       <NON_REF>       .       .       END=14718       GT:DP:GQ:MIN_DP:PL       0/0:37:99:37:0,99,1192
20_10037292_10066351    14719   .       T       <NON_REF>       .       .       END=14719       GT:DP:GQ:MIN_DP:PL       0/0:36:82:36:0,82,1087
20_10037292_10066351    14720   .       T       <NON_REF>       .       .       END=14737       GT:DP:GQ:MIN_DP:PL       0/0:42:99:37:0,100,1160
```

Celles-ci représentent des régions non-variantes où l'appeleur de variants n'a trouvé aucune preuve de variation, il a donc capturé certaines statistiques décrivant son niveau de confiance dans l'absence de variation. Cela permet de distinguer deux cas très différents : (1) il y a des données de bonne qualité montrant que l'échantillon est homozygote-référence, et (2) il n'y a pas assez de bonnes données disponibles pour faire une détermination dans un sens ou dans l'autre.

Dans un GVCF, il y a généralement beaucoup de telles lignes non-variantes, avec un plus petit nombre d'enregistrements de variants parsemés parmi elles. Essayez d'exécuter `head -176` sur le GVCF pour charger seulement les 176 premières lignes du fichier et trouver un appel de variant réel.

```console title="Sortie" linenums="174"
20_10037292_10066351    3479    .       T       <NON_REF>       .       .       END=3479        GT:DP:GQ:MIN_DP:PL       0/0:34:36:34:0,36,906
20_10037292_10066351    3480    .       C       CT,<NON_REF>    503.03  .       DP=23;ExcessHet=0.0000;MLEAC=2,0;MLEAF=1.00,0.00;RAW_MQandDP=82800,23    GT:AD:DP:GQ:PL:SB       1/1:0,18,0:18:54:517,54,0,517,54,517:0,0,7,11
20_10037292_10066351    3481    .       T       <NON_REF>       .       .       END=3481        GT:DP:GQ:MIN_DP:PL       0/0:21:51:21:0,51,765
```

La deuxième ligne montre le premier enregistrement de variant dans le fichier, qui correspond au premier variant dans le fichier VCF que nous avons examiné dans la Partie 1.

Tout comme le VCF original, le fichier de sortie GVCF est également accompagné d'un fichier d'index, appelé `reads_mother.g.vcf.idx`.

#### 0.2.3. Répéter le processus sur les deux autres échantillons

Afin de tester l'étape de génotypage conjoint, nous avons besoin de GVCFs pour les trois échantillons, générons-les donc manuellement maintenant.

```bash
gatk HaplotypeCaller \
        -R /data/ref/ref.fasta \
        -I /data/bam/reads_father.bam \
        -O reads_father.g.vcf \
        -L /data/ref/intervals.bed \
        -ERC GVCF
```

<!--
??? success "Sortie de la commande"

    ```console

    ```
-->

```bash
gatk HaplotypeCaller \
        -R /data/ref/ref.fasta \
        -I /data/bam/reads_son.bam \
        -O reads_son.g.vcf \
        -L /data/ref/intervals.bed \
        -ERC GVCF
```

<!--
??? success "Sortie de la commande"

    ```console

    ```
-->

Une fois terminé, vous devriez avoir trois fichiers se terminant par `.g.vcf` dans votre répertoire actuel (un par échantillon) et leurs fichiers d'index respectifs se terminant par `.g.vcf.idx`.

### 0.3. Exécuter le génotypage conjoint

Maintenant que nous avons tous les GVCFs, nous pouvons enfin essayer l'approche de génotypage conjoint pour générer des appels de variants pour une cohorte d'échantillons.
Pour rappel, c'est une méthode en deux étapes qui consiste à combiner les données de tous les GVCFs dans un magasin de données, puis à exécuter l'analyse de génotypage conjoint proprement dite pour générer le VCF final d'appels de variants conjoints.

#### 0.3.1. Combiner tous les GVCFs par échantillon

Cette première étape utilise un autre outil GATK, appelé GenomicsDBImport, pour combiner les données de tous les GVCFs dans un magasin de données GenomicsDB.

```bash
gatk GenomicsDBImport \
    -V reads_mother.g.vcf \
    -V reads_father.g.vcf \
    -V reads_son.g.vcf \
    -L /data/ref/intervals.bed \
    --genomicsdb-workspace-path family_trio_gdb
```

<!--
??? success "Sortie de la commande"

    ```console

    ```
-->

La sortie de cette étape est effectivement un répertoire contenant un ensemble de répertoires supplémentaires imbriqués contenant les données de variants combinées sous la forme de plusieurs fichiers différents.
Vous pouvez explorer mais vous verrez rapidement que ce format de magasin de données n'est pas destiné à être lu directement par des humains.

!!! note

    GATK inclut des outils qui permettent d'inspecter et d'extraire les données d'appel de variants du magasin de données selon les besoins.

#### 0.3.2. Exécuter l'analyse de génotypage conjoint proprement dite

Cette deuxième étape utilise encore un autre outil GATK, appelé GenotypeGVCFs, pour recalculer les statistiques de variants et les génotypes individuels à la lumière des données disponibles sur tous les échantillons de la cohorte.

```bash
gatk GenotypeGVCFs \
    -R /data/ref/ref.fasta \
    -V gendb://family_trio_gdb \
    -O family_trio.vcf
```

<!--
??? success "Sortie de la commande"

    ```console

    ```
-->

Ceci crée le fichier de sortie VCF `family_trio.vcf` dans le répertoire de travail actuel dans le conteneur.
C'est un autre fichier raisonnablement petit, vous pouvez donc utiliser `cat` pour visualiser son contenu et faire défiler pour trouver les premières lignes de variants.

```console title="family_trio.vcf" linenums="40"
#CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO    FORMAT  reads_father    reads_mother    reads_son
20_10037292_10066351    3480    .       C       CT      1625.89 .       AC=5;AF=0.833;AN=6;BaseQRankSum=0.220;DP=85;ExcessHet=0.0000;FS=2.476;MLEAC=5;MLEAF=0.833;MQ=60.00;MQRankSum=0.00;QD=21.68;ReadPosRankSum=-1.147e+00;SOR=0.487    GT:AD:DP:GQ:PL  0/1:15,16:31:99:367,0,375       1/1:0,18:18:54:517,54,0 1/1:0,26:26:78:756,78,0
20_10037292_10066351    3520    .       AT      A       1678.89 .       AC=5;AF=0.833;AN=6;BaseQRankSum=1.03;DP=80;ExcessHet=0.0000;FS=2.290;MLEAC=5;MLEAF=0.833;MQ=60.00;MQRankSum=0.00;QD=22.39;ReadPosRankSum=0.701;SOR=0.730 GT:AD:DP:GQ:PL   0/1:18,13:31:99:296,0,424       1/1:0,18:18:54:623,54,0 1/1:0,26:26:78:774,78,0
20_10037292_10066351    3529    .       T       A       154.29  .       AC=1;AF=0.167;AN=6;BaseQRankSum=-5.440e-01;DP=104;ExcessHet=0.0000;FS=1.871;MLEAC=1;MLEAF=0.167;MQ=60.00;MQRankSum=0.00;QD=7.71;ReadPosRankSum=-1.158e+00;SOR=1.034       GT:AD:DP:GQ:PL  0/0:44,0:44:99:0,112,1347       0/1:12,8:20:99:163,0,328        0/0:39,0:39:99:0,105,1194
```

Cela ressemble plus au VCF original que nous avons généré dans la Partie 1, sauf que cette fois nous avons des informations au niveau du génotype pour les trois échantillons.
Les trois dernières colonnes du fichier sont les blocs de génotypes pour les échantillons, listés par ordre alphabétique.

Si nous examinons les génotypes appelés pour notre trio familial de test pour le tout premier variant, nous voyons que le père est hétérozygote-variant (`0/1`), et la mère et le fils sont tous deux homozygotes-variants (`1/1`).

C'est finalement l'information que nous cherchons à extraire du jeu de données ! Alors encapsulons tout cela dans un workflow Nextflow afin de pouvoir le faire à grande échelle.

#### 0.3.3. Quitter le conteneur GATK

```bash
exit
```

### À retenir

Vous savez comment exécuter les commandes individuelles impliquées dans l'appel conjoint de variants dans le terminal pour vérifier qu'elles produiront les informations souhaitées.

### Et ensuite ?

Encapsuler ces commandes dans un pipeline réel.

---

## 1. Modifier l'étape d'appel de variants par échantillon pour produire un GVCF

La bonne nouvelle est que nous n'avons pas besoin de tout recommencer, puisque nous avons déjà écrit un workflow qui fait une partie de ce travail dans la Partie 1.
Cependant, ce pipeline produit des fichiers VCF, alors que maintenant nous voulons des fichiers GVCF afin de faire le génotypage conjoint.
Nous devons donc commencer par activer le mode d'appel de variants GVCF et mettre à jour l'extension du fichier de sortie.

!!! note

    Pour plus de commodité, nous allons travailler avec une nouvelle copie du workflow GATK tel qu'il se présente à la fin de la Partie 1, mais sous un nom différent : `genomics-2.nf`.

### 1.1. Indiquer à HaplotypeCaller d'émettre un GVCF et mettre à jour l'extension de sortie

Ouvrons le fichier `genomics-2.nf` dans l'éditeur de code.
Il devrait être très familier, mais n'hésitez pas à l'exécuter si vous voulez vous assurer qu'il fonctionne comme prévu.

Nous allons commencer par apporter deux modifications :

- Ajouter le paramètre `-ERC GVCF` à la commande GATK HaplotypeCaller ;
- Mettre à jour le chemin du fichier de sortie pour utiliser l'extension `.g.vcf` correspondante, selon la convention GATK.

Assurez-vous d'ajouter une barre oblique inverse (`\`) à la fin de la ligne précédente lorsque vous ajoutez `-ERC GVCF`.

=== "Après"

    ```groovy title="genomics-2.nf" linenums="56" hl_lines="4 6"
        """
        gatk HaplotypeCaller \
            -R ${ref_fasta} \
            -I ${input_bam} \
            -O ${input_bam}.g.vcf \
            -L ${interval_list} \
            -ERC GVCF
        """
    ```

=== "Avant"

    ```groovy title="genomics-2.nf" linenums="56" hl_lines="4"
        """
        gatk HaplotypeCaller \
            -R ${ref_fasta} \
            -I ${input_bam} \
            -O ${input_bam}.vcf \
            -L ${interval_list}
        """
    ```

Et c'est tout ce qu'il faut pour faire passer HaplotypeCaller à la génération de GVCFs au lieu de VCFs, n'est-ce pas ?

### 1.2. Exécuter le pipeline pour vérifier que vous pouvez générer des GVCFs

La commande d'exécution Nextflow est la même qu'avant, sauf pour le nom du fichier de workflow lui-même.
Assurez-vous de le mettre à jour de manière appropriée.

```bash
nextflow run genomics-2.nf
```

??? success "Sortie de la commande"

    ```console
    N E X T F L O W   ~  version 25.10.2

    ┃ Launching `genomics-2.nf` [crazy_venter] DSL2 - revision: a2d6f6f09f

    executor >  local (6)
    [f1/8d8486] SAMTOOLS_INDEX (1)       | 3 of 3 ✔
    [72/3249ca] GATK_HAPLOTYPECALLER (3) | 0 of 3
    ERROR ~ Error executing process > 'GATK_HAPLOTYPECALLER (2)'

    Caused by:
      Missing output file(s) `reads_son.bam.vcf` expected by process `GATK_HAPLOTYPECALLER (2)`

    Command executed:

      gatk HaplotypeCaller         -R ref.fasta         -I reads_son.bam         -O reads_son.bam.g.vcf         -L intervals.bed         -ERC GVCF
    ```

Et la sortie est... toute rouge ! Oh non.

La commande qui a été exécutée est correcte, nous avions donc raison de penser que c'était suffisant pour changer le comportement de l'outil GATK.
Mais regardez cette ligne sur le fichier de sortie manquant. Remarquez-vous quelque chose ?

C'est exact, nous avons oublié de dire à Nextflow d'attendre un nouveau nom de fichier. Oups.

### 1.3. Mettre à jour l'extension du fichier de sortie dans le bloc de sorties du processus également

Parce qu'il ne suffit pas de simplement changer l'extension du fichier dans la commande de l'outil elle-même, vous devez également dire à Nextflow que le nom du fichier de sortie attendu a changé.

=== "Après"

    ```groovy title="genomics-2.nf" linenums="50" hl_lines="2 3"
        output:
        path "${input_bam}.g.vcf"     , emit: vcf
        path "${input_bam}.g.vcf.idx" , emit: idx
    ```

=== "Avant"

    ```groovy title="genomics-2.nf" linenums="50" hl_lines="2 3"
        output:
        path "${input_bam}.vcf"     , emit: vcf
        path "${input_bam}.vcf.idx" , emit: idx
    ```

### 1.4. Mettre à jour les cibles de publication pour les nouvelles sorties GVCF

Puisque nous produisons maintenant des GVCFs au lieu de VCFs, nous devrions mettre à jour la section `publish:` du workflow pour utiliser des noms plus descriptifs.
Nous organiserons également les fichiers GVCF dans leur propre sous-répertoire pour plus de clarté.

=== "Après"

    ```groovy title="genomics-2.nf" linenums="88" hl_lines="3 4"
        publish:
        indexed_bam = SAMTOOLS_INDEX.out
        gvcf = GATK_HAPLOTYPECALLER.out.vcf
        gvcf_idx = GATK_HAPLOTYPECALLER.out.idx
    ```

=== "Avant"

    ```groovy title="genomics-2.nf" linenums="88"
        publish:
        indexed_bam = SAMTOOLS_INDEX.out
        vcf = GATK_HAPLOTYPECALLER.out.vcf
        vcf_idx = GATK_HAPLOTYPECALLER.out.idx
    ```

### 1.5. Mettre à jour le bloc output pour la nouvelle structure de répertoires

Nous devons également mettre à jour le bloc `output` pour placer les fichiers GVCF dans un sous-répertoire `gvcf`.

=== "Après"

    ```groovy title="genomics-2.nf" linenums="94" hl_lines="3 5 6 8 9"
    output {
        indexed_bam {
            path 'indexed_bam'
        }
        gvcf {
            path 'gvcf'
        }
        gvcf_idx {
            path 'gvcf'
        }
    }
    ```

=== "Avant"

    ```groovy title="genomics-2.nf" linenums="94"
    output {
        indexed_bam {
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

### 1.6. Exécuter à nouveau le pipeline

Exécutons-le avec `-resume` cette fois.

```bash
nextflow run genomics-2.nf -resume
```

??? success "Sortie de la commande"

    ```console
    N E X T F L O W   ~  version 25.10.2

    ┃ Launching `genomics-2.nf` [nostalgic_franklin] DSL2 - revision: f2c0a93c6a

    executor >  local (3)
    [cc/fbc705] SAMTOOLS_INDEX (3)       | 3 of 3, cached: 3 ✔
    [27/0d7eb9] GATK_HAPLOTYPECALLER (2) | 3 of 3 ✔
    ```

Cette fois, cela fonctionne.

La sortie Nextflow elle-même ne semble pas différente (comparée à une exécution réussie en mode VCF normal), mais maintenant nous pouvons trouver les fichiers `.g.vcf` et leurs fichiers d'index respectifs, pour les trois échantillons, organisés dans des sous-répertoires.

??? abstract "Contenu du répertoire (liens symboliques raccourcis)"

    ```console
    results_genomics/
    ├── gvcf/
    │   ├── reads_father.bam.g.vcf -> */27/0d7eb9*/reads_father.bam.g.vcf
    │   ├── reads_father.bam.g.vcf.idx -> */27/0d7eb9*/reads_father.bam.g.vcf.idx
    │   ├── reads_mother.bam.g.vcf -> */e4/4ed55e*/reads_mother.bam.g.vcf
    │   ├── reads_mother.bam.g.vcf.idx -> */e4/4ed55e*/reads_mother.bam.g.vcf.idx
    │   ├── reads_son.bam.g.vcf -> */08/e95962*/reads_son.bam.g.vcf
    │   └── reads_son.bam.g.vcf.idx -> */08/e95962*/reads_son.bam.g.vcf.idx
    └── indexed_bam/
        ├── reads_father.bam -> */9a/c7a873*/reads_father.bam
        ├── reads_father.bam.bai -> */9a/c7a873*/reads_father.bam.bai
        ├── reads_mother.bam -> */f1/8d8486*/reads_mother.bam
        ├── reads_mother.bam.bai -> */f1/8d8486*/reads_mother.bam.bai
        ├── reads_son.bam -> */cc/fbc705*/reads_son.bam
        └── reads_son.bam.bai -> */cc/fbc705*/reads_son.bam.bai
    ```

Si vous ouvrez l'un des fichiers GVCF et le parcourez, vous pouvez vérifier que GATK HaplotypeCaller a produit des fichiers GVCF comme demandé.

### À retenir

Bon, celle-ci était minime en termes d'apprentissage Nextflow...
Mais c'était une belle occasion de réitérer l'importance du bloc de sortie du processus !

### Et ensuite ?

Apprendre à collecter le contenu d'un canal et à le transmettre au processus suivant comme une entrée unique.

---

## 2. Collecter et combiner les données GVCF de tous les échantillons

Nous devons maintenant combiner les données de tous les GVCFs par échantillon dans une forme qui supporte l'analyse de génotypage conjoint que nous voulons faire.

### 2.1. Définir le processus qui combinera les GVCFs

Pour rappel de ce que nous avons fait plus tôt dans la section d'échauffement, combiner les GVCFs est un travail pour l'outil GATK GenomicsDBImport, qui produira un magasin de données dans le format dit GenomicsDB.

Écrivons un nouveau processus pour définir comment cela va fonctionner, basé sur la commande que nous avons utilisée plus tôt dans la section d'échauffement.

```groovy title="genomics-2.nf" linenums="66"
/*
 * Combiner les GVCFs dans un magasin de données GenomicsDB
 */
process GATK_GENOMICSDB {

    container "community.wave.seqera.io/library/gatk4:4.5.0.0--730ee8817e436867"

    input:
    path all_gvcfs
    path all_idxs
    path interval_list
    val cohort_name

    output:
    path "${cohort_name}_gdb"

    script:
    """
    gatk GenomicsDBImport \
        -V ${all_gvcfs} \
        -L ${interval_list} \
        --genomicsdb-workspace-path ${cohort_name}_gdb
    """
}
```

Qu'en pensez-vous, cela semble-t-il raisonnable ?

Branchons-le et voyons ce qui se passe.

### 2.2. Ajouter un paramètre `cohort_name` avec une valeur par défaut

Nous devons fournir un nom arbitraire pour la cohorte.
Plus tard dans la série de formations, vous apprendrez comment utiliser les métadonnées d'échantillons pour ce genre de choses, mais pour l'instant, nous déclarons simplement un paramètre CLI en utilisant `params` et lui donnons une valeur par défaut pour plus de commodité.

```groovy title="genomics-2.nf" linenums="16"
    // Nom de base pour le fichier de sortie final
    cohort_name: String = "family_trio"
```

### 2.3. Rassembler les sorties de GATK_HAPLOTYPECALLER pour tous les échantillons

Si nous devions simplement brancher le canal de sortie du processus `GATK_HAPLOTYPECALLER` tel quel, Nextflow appellerait le processus sur chaque GVCF d'échantillon séparément.
Cependant, nous voulons regrouper les trois GVCFs (et leurs fichiers d'index) de telle manière que Nextflow les transmette tous ensemble à un seul appel de processus.

Bonne nouvelle : nous pouvons le faire en utilisant l'opérateur de canal `collect()`. Ajoutons les lignes suivantes au corps du `workflow`, juste après l'appel à GATK_HAPLOTYPECALLER :

```groovy title="genomics-2.nf" linenums="118"
// Collecter les sorties d'appel de variants pour tous les échantillons
all_gvcfs_ch = GATK_HAPLOTYPECALLER.out.vcf.collect()
all_idxs_ch = GATK_HAPLOTYPECALLER.out.idx.collect()
```

Cela semble-t-il un peu compliqué ? Décomposons cela et traduisons-le en langage clair.

1. Nous prenons le canal de sortie du processus `GATK_HAPLOTYPECALLER`, référencé en utilisant la propriété `.out`.
2. Chaque « élément » sortant du canal est une paire de fichiers : le GVCF et son fichier d'index, dans cet ordre parce que c'est l'ordre dans lequel ils sont listés dans le bloc de sortie du processus. Heureusement, parce que dans la dernière session nous avons nommé les sorties de ce processus (en utilisant `emit:`), nous pouvons sélectionner les GVCFs d'un côté en ajoutant `.vcf` et les fichiers d'index de l'autre en ajoutant `.idx` après la propriété `.out`. Si nous n'avions pas nommé ces sorties, nous aurions dû y faire référence par `.out[0]` et `.out[1]`, respectivement.
3. Nous ajoutons l'opérateur de canal `collect()` pour regrouper tous les fichiers GVCF ensemble dans un seul élément d'un nouveau canal appelé `all_gvcfs_ch`, et faisons de même avec les fichiers d'index pour former le nouveau canal appelé `all_idxs_ch`.

!!! tip

    Si vous avez du mal à visualiser exactement ce qui se passe ici, rappelez-vous que vous pouvez utiliser l'opérateur `view()` pour inspecter le contenu des canaux avant et après l'application des opérateurs de canaux.

Les canaux `all_gvcfs_ch` et `all_idxs_ch` résultants sont ce que nous allons brancher dans le processus `GATK_GENOMICSDB` que nous venons d'écrire.

!!! note

    Au cas où vous vous poseriez la question, nous collectons les GVCFs et leurs fichiers d'index séparément parce que la commande GATK GenomicsDBImport ne veut voir que les chemins des fichiers GVCF. Heureusement, puisque Nextflow préparera tous les fichiers ensemble pour l'exécution, nous n'avons pas à nous soucier de l'ordre des fichiers comme nous l'avons fait pour les BAMs et leur index dans la Partie 1.

### 2.4. Ajouter un appel au bloc workflow pour exécuter GATK_GENOMICSDB

Nous avons un processus, et nous avons des canaux d'entrée. Nous devons juste ajouter l'appel de processus.

```groovy title="genomics-2.nf" linenums="122"
    // Combiner les GVCFs dans un magasin de données GenomicsDB
    GATK_GENOMICSDB(
        all_gvcfs_ch,
        all_idxs_ch,
        intervals_file,
        params.cohort_name
    )
```

Ok, tout est branché.

### 2.5. Exécuter le workflow

Voyons si cela fonctionne.

```bash
nextflow run genomics-2.nf -resume
```

??? failure "Sortie de la commande"

    ```console
    N E X T F L O W   ~  version 25.10.2

    ┃ Launching `genomics-2.nf` [disturbed_bell] DSL2 - revision: 57942246cc

    executor >  local (1)
    [f1/8d8486] SAMTOOLS_INDEX (1)       | 3 of 3, cached: 3 ✔
    [e4/4ed55e] GATK_HAPLOTYPECALLER (2) | 3 of 3, cached: 3 ✔
    [51/d350ea] GATK_GENOMICSDB          | 0 of 1
    ERROR ~ Error executing process > 'GATK_GENOMICSDB'

    Caused by:
      Process `GATK_GENOMICSDB` terminated with an error exit status (1)

    Command executed:

      gatk GenomicsDBImport         -V reads_son.bam.g.vcf reads_father.bam.g.vcf reads_mother.bam.g.vcf         -L intervals.bed         --genomicsdb-workspace-path family_trio_gdb
    ```

Il s'exécute assez rapidement, puisque nous exécutons avec `-resume`, mais il échoue !

Ah. Du côté positif, nous voyons que Nextflow a récupéré le processus `GATK_GENOMICSDB`, et l'a spécifiquement appelé une seule fois.
Cela suggère que l'approche `collect()` a fonctionné, dans une certaine mesure.
Mais, et c'est un gros mais, l'appel de processus a échoué.

Lorsque nous fouillons dans la sortie de la console ci-dessus, nous pouvons voir que la commande exécutée n'est pas correcte.

Pouvez-vous repérer l'erreur ?
Regardez cette partie : `-V reads_father.bam.g.vcf reads_son.bam.g.vcf reads_mother.bam.g.vcf`

Nous avons donné à `gatk GenomicsDBImport` plusieurs fichiers GVCF pour un seul argument `-V`, mais l'outil attend un argument `-V` séparé pour chaque fichier GVCF.

Pour rappel, voici la commande que nous avons exécutée dans le conteneur :

```bash
gatk GenomicsDBImport \
    -V reads_mother.g.vcf \
    -V reads_father.g.vcf \
    -V reads_son.g.vcf \
    -L /data/ref/intervals.bed \
    --genomicsdb-workspace-path family_trio_gdb
```

Cela signifie donc que nous devons d'une manière ou d'une autre transformer notre ensemble de fichiers GVCF en une chaîne de commande correctement formatée.

### 2.6. Construire une ligne de commande avec un argument `-V` séparé pour chaque GVCF d'entrée

C'est là que le fait que Nextflow soit basé sur Groovy s'avère pratique, car cela va nous permettre d'utiliser des manipulations de chaînes assez simples pour construire la chaîne de commande nécessaire.

Spécifiquement, en utilisant cette syntaxe : `all_gvcfs.collect { gvcf -> "-V ${gvcf}" }.join(' ')`

Encore une fois, décomposons cela en ses composants.

1. D'abord, nous prenons le contenu du canal d'entrée `all_gvcfs` et appliquons `.collect()` dessus (tout comme précédemment).
2. Cela nous permet de passer chaque chemin de fichier GVCF individuel dans l'ensemble à la **closure**, `{ gvcf -> "-V ${gvcf}" }`, où `gvcf` fait référence à ce chemin de fichier GVCF.
   La closure est une mini-fonction que nous utilisons pour préfixer `-V ` au chemin de fichier, sous la forme de `"-V ${gvcf}"`.
3. Ensuite, nous utilisons `.join(' ')` pour concaténer les trois chaînes avec un seul espace comme séparateur.

Avec un exemple concret, cela ressemble à ceci :

1. Nous avons trois fichiers :

   `[A.ext, B.ext, C.ext]`

2. La closure modifie chacun pour créer les chaînes :

   `"-V A.ext", "-V B.ext", "-V C.ext"`

3. L'opération `.join(' ')` génère la chaîne finale :

   `"-V A.ext -V B.ext -V C.ext"`

Une fois que nous avons cette chaîne, nous pouvons l'assigner à une variable locale, `gvcfs_line`, définie avec le mot-clé `def` :

`def gvcfs_line = all_gvcfs.collect { gvcf -> "-V ${gvcf}" }.join(' ')`

Ok, nous avons donc notre manipulation de chaîne. Où la mettons-nous ?

Nous voulons que cela aille quelque part dans la définition du processus, car nous voulons le faire _après_ avoir canalisé les chemins de fichiers GVCF dans le processus.
C'est parce que Nextflow doit les voir comme des chemins de fichiers afin de préparer les fichiers eux-mêmes correctement pour l'exécution.

Mais _où_ dans le processus pouvons-nous ajouter cela ?

Fait amusant : vous pouvez ajouter du code arbitraire après `script:` et avant les `"""` !

Parfait, ajoutons notre ligne de manipulation de chaîne là alors, et mettons à jour la commande `gatk GenomicsDBImport` pour utiliser la chaîne concaténée qu'elle produit.

=== "Après"

    ```groovy title="genomics-2.nf" linenums="87"  hl_lines="2 5"
        script:
        def gvcfs_line = all_gvcfs.collect { gvcf -> "-V ${gvcf}" }.join(' ')
        """
        gatk GenomicsDBImport \
            ${gvcfs_line} \
            -L ${interval_list} \
            --genomicsdb-workspace-path ${cohort_name}_gdb
        """
    ```

=== "Avant"

    ```groovy title="genomics-2.nf" linenums="87"  hl_lines="4"
        script:
        """
        gatk GenomicsDBImport \
            -V ${all_gvcfs} \
            -L ${interval_list} \
            --genomicsdb-workspace-path ${cohort_name}_gdb
        """
    ```

Cela devrait être tout ce qui est nécessaire pour fournir les entrées à `gatk GenomicsDBImport` correctement.

!!! tip

    Lorsque vous mettez à jour la commande `gatk GenomicsDBImport`, assurez-vous de supprimer le préfixe `-V ` lorsque vous remplacez par la variable `${gvcfs_line}`.

### 2.7. Exécuter le workflow pour vérifier qu'il génère la sortie GenomicsDB comme prévu

D'accord, voyons si cela a résolu le problème.

```bash
nextflow run genomics-2.nf -resume
```

??? success "Sortie de la commande"

    ```console
    N E X T F L O W   ~  version 25.10.2

    ┃ Launching `genomics-2.nf` [peaceful_gates] DSL2 - revision: ca0bf847ed

    executor >  local (1)
    [cc/fbc705] SAMTOOLS_INDEX (3)       | 3 of 3, cached: 3 ✔
    [27/0d7eb9] GATK_HAPLOTYPECALLER (2) | 3 of 3, cached: 3 ✔
    [76/d13861] GATK_GENOMICSDB          | 1 of 1 ✔
    ```

Aha ! Il semble que cela fonctionne maintenant.

Les deux premières étapes ont été ignorées avec succès, et la troisième étape a fonctionné comme un charme cette fois.
Le magasin de données GenomicsDB est créé dans le répertoire de travail mais n'est pas publié dans les résultats, puisque c'est juste un format intermédiaire que nous utiliserons pour le génotypage conjoint.

Au fait, nous n'avons rien eu à faire de spécial pour gérer le fait que la sortie soit un répertoire au lieu d'un seul fichier.

### À retenir

Vous savez maintenant comment collecter les sorties d'un canal et les regrouper comme une entrée unique pour un autre processus.
Vous savez également comment construire une ligne de commande pour fournir des entrées à un outil donné avec la syntaxe appropriée.

### Et ensuite ?

Apprendre à ajouter une deuxième commande au même processus.

---

## 3. Exécuter l'étape de génotypage conjoint dans le même processus

Maintenant que nous avons les appels de variants génomiques combinés, nous pouvons exécuter l'outil de génotypage conjoint, qui produira la sortie finale qui nous intéresse vraiment : le VCF d'appels de variants au niveau de la cohorte.

Pour des raisons logistiques, nous décidons d'inclure le génotypage conjoint dans le même processus.

### 3.1. Renommer le processus de GATK_GENOMICSDB à GATK_JOINTGENOTYPING

Puisque le processus exécutera plus d'un outil, nous changeons son nom pour faire référence à l'opération globale plutôt qu'à un seul nom d'outil.

=== "Après"

    ```groovy title="genomics-2.nf"
    /*
     * Combiner les GVCFs dans un magasin de données GenomicsDB et exécuter le génotypage conjoint pour produire des appels au niveau de la cohorte
     */
    process GATK_JOINTGENOTYPING {
    ```

=== "Avant"

    ```groovy title="genomics-2.nf"
    /*
     * Combiner les GVCFs dans un magasin de données GenomicsDB
     */
    process GATK_GENOMICSDB {
    ```

N'oubliez pas de garder vos noms de processus aussi descriptifs que possible, pour maximiser la lisibilité pour vos collègues — et votre futur vous !

### 3.2. Ajouter la commande de génotypage conjoint au processus GATK_JOINTGENOTYPING

Ajoutez simplement la deuxième commande après la première dans la section script.

=== "Après"

    ```groovy title="genomics-2.nf" linenums="89"  hl_lines="6-10"
        """
        gatk GenomicsDBImport \
            ${gvcfs_line} \
            -L ${interval_list} \
            --genomicsdb-workspace-path ${cohort_name}_gdb

        gatk GenotypeGVCFs \
            -R ${ref_fasta} \
            -V gendb://${cohort_name}_gdb \
            -L ${interval_list} \
            -O ${cohort_name}.joint.vcf
        """
    ```

=== "Avant"

    ```groovy title="genomics-2.nf" linenums="89"
        """
        gatk GenomicsDBImport \
            ${gvcfs_line} \
            -L ${interval_list} \
            --genomicsdb-workspace-path ${cohort_name}_gdb
        """
    ```

Les deux commandes seront exécutées en série, de la même manière qu'elles le seraient si nous les exécutions manuellement dans le terminal.

### 3.3. Ajouter les fichiers du génome de référence aux définitions d'entrée du processus GATK_JOINTGENOTYPING

La deuxième commande nécessite les fichiers du génome de référence, nous devons donc les ajouter aux entrées du processus.

=== "Après"

    ```groovy title="genomics-2.nf" linenums="78"  hl_lines="6-8"
        input:
        path all_gvcfs
        path all_idxs
        path interval_list
        val cohort_name
        path ref_fasta
        path ref_index
        path ref_dict
    ```

=== "Avant"

    ```groovy title="genomics-2.nf" linenums="78"
        input:
        path all_gvcfs
        path all_idxs
        path interval_list
        val cohort_name
    ```

Il peut sembler ennuyeux de taper tout cela, mais rappelez-vous, vous ne les tapez qu'une fois, et ensuite vous pouvez exécuter le workflow un million de fois. Ça en vaut la peine ?

### 3.4. Mettre à jour la définition de sortie du processus pour émettre le VCF d'appels de variants au niveau de la cohorte

Nous ne nous soucions pas vraiment de sauvegarder le magasin de données GenomicsDB, qui n'est qu'un format intermédiaire qui existe seulement pour des raisons logistiques, nous pouvons donc simplement le supprimer du bloc de sortie si nous le souhaitons.

La sortie qui nous intéresse vraiment est le VCF produit par la commande de génotypage conjoint.

=== "Après"

    ```groovy title="genomics-2.nf" linenums="87" hl_lines="2 3"
        output:
        path "${cohort_name}.joint.vcf"     , emit: vcf
        path "${cohort_name}.joint.vcf.idx" , emit: idx
    ```

=== "Avant"

    ```groovy title="genomics-2.nf" linenums="87"
        output:
        path "${cohort_name}_gdb"
    ```

Nous avons presque fini !

### 3.5. Mettre à jour l'appel de processus de GATK_GENOMICSDB à GATK_JOINTGENOTYPING

N'oublions pas de renommer l'appel de processus dans le corps du workflow de GATK_GENOMICSDB à GATK_JOINTGENOTYPING. Et pendant que nous y sommes, nous devrions également ajouter les fichiers du génome de référence comme entrées, puisque nous devons les fournir à l'outil de génotypage conjoint.

=== "Après"

    ```groovy title="genomics-2.nf" linenums="126"
    // Combiner les GVCFs dans un magasin de données GenomicsDB et appliquer le génotypage conjoint
    GATK_JOINTGENOTYPING(
        all_gvcfs_ch,
        all_idxs_ch,
        intervals_file,
        params.cohort_name,
        ref_file,
        ref_index_file,
        ref_dict_file
    )
    ```

=== "Avant"

    ```groovy title="genomics-2.nf" linenums="126"
    // Combiner les GVCFs dans un magasin de données GenomicsDB
    GATK_GENOMICSDB(
        all_gvcfs_ch,
        all_idxs_ch,
        intervals_file,
        params.cohort_name
    )
    ```

Maintenant le processus est complètement branché.

### 3.6. Ajouter le VCF conjoint à la section publish

Nous devons publier les sorties VCF conjointes du nouveau processus.
Ajoutez ces lignes à la section `publish:` du workflow :

```groovy title="genomics-2.nf" linenums="145"
    joint_vcf = GATK_JOINTGENOTYPING.out.vcf
    joint_vcf_idx = GATK_JOINTGENOTYPING.out.idx
```

### 3.7. Ajouter les cibles VCF conjointes au bloc output

Enfin, ajoutez des cibles de sortie pour les fichiers VCF conjoints.
Nous les placerons à la racine du répertoire de résultats puisque c'est la sortie finale.

```groovy title="genomics-2.nf" linenums="157"
    joint_vcf {
        path '.'
    }
    joint_vcf_idx {
        path '.'
    }
```

Maintenant tout devrait être complètement branché.

### 3.8. Exécuter le workflow

Enfin, nous pouvons exécuter le workflow modifié...

```bash
nextflow run genomics-2.nf -resume
```

??? success "Sortie de la commande"

    ```console
    N E X T F L O W   ~  version 25.10.2

    ┃ Launching `genomics-2.nf` [crazy_marconi] DSL2 - revision: 5da9afc841

    executor >  local (1)
    [9a/c7a873] SAMTOOLS_INDEX (2)       | 3 of 3, cached: 3 ✔
    [e4/4ed55e] GATK_HAPLOTYPECALLER (2) | 3 of 3, cached: 3 ✔
    [a6/7cc8ed] GATK_JOINTGENOTYPING     | 1 of 1 ✔
    ```

Et ça fonctionne !

Vous trouverez le fichier de sortie final, `family_trio.joint.vcf` (et son index de fichier), dans le répertoire de résultats.

??? abstract "Contenu du répertoire (liens symboliques raccourcis)"

    ```console
    results_genomics/
    ├── family_trio.joint.vcf -> */a6/7cc8ed*/family_trio.joint.vcf
    ├── family_trio.joint.vcf.idx -> */a6/7cc8ed*/family_trio.joint.vcf.idx
    ├── gvcf/
    │   ├── reads_father.bam.g.vcf -> */27/0d7eb9*/reads_father.bam.g.vcf
    │   ├── reads_father.bam.g.vcf.idx -> */27/0d7eb9*/reads_father.bam.g.vcf.idx
    │   ├── reads_mother.bam.g.vcf -> */e4/4ed55e*/reads_mother.bam.g.vcf
    │   ├── reads_mother.bam.g.vcf.idx -> */e4/4ed55e*/reads_mother.bam.g.vcf.idx
    │   ├── reads_son.bam.g.vcf -> */08/e95962*/reads_son.bam.g.vcf
    │   └── reads_son.bam.g.vcf.idx -> */08/e95962*/reads_son.bam.g.vcf.idx
    └── indexed_bam/
        ├── reads_father.bam -> */9a/c7a873*/reads_father.bam
        ├── reads_father.bam.bai -> */9a/c7a873*/reads_father.bam.bai
        ├── reads_mother.bam -> */f1/8d8486*/reads_mother.bam
        ├── reads_mother.bam.bai -> */f1/8d8486*/reads_mother.bam.bai
        ├── reads_son.bam -> */cc/fbc705*/reads_son.bam
        └── reads_son.bam.bai -> */cc/fbc705*/reads_son.bam.bai
    ```

Si vous êtes du genre sceptique (ou curieux·se), vous pouvez cliquer sur le fichier VCF conjoint pour l'ouvrir et vérifier que le workflow a généré les mêmes appels de variants que vous avez obtenus en exécutant les outils manuellement au début de cette section.

```console title="family_trio.joint.vcf" linenums="40"
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	reads_father	reads_mother	reads_son
20_10037292_10066351	3480	.	C	CT	1625.89	.	AC=5;AF=0.833;AN=6;BaseQRankSum=0.220;DP=85;ExcessHet=0.0000;FS=2.476;MLEAC=5;MLEAF=0.833;MQ=60.00;MQRankSum=0.00;QD=21.68;ReadPosRankSum=-1.147e+00;SOR=0.487	GT:AD:DP:GQ:PL	0/1:15,16:31:99:367,0,375	1/1:0,18:18:54:517,54,0	1/1:0,26:26:78:756,78,0
20_10037292_10066351	3520	.	AT	A	1678.89	.	AC=5;AF=0.833;AN=6;BaseQRankSum=1.03;DP=80;ExcessHet=0.0000;FS=2.290;MLEAC=5;MLEAF=0.833;MQ=60.00;MQRankSum=0.00;QD=22.39;ReadPosRankSum=0.701;SOR=0.730	GT:AD:DP:GQ:PL	0/1:18,13:31:99:296,0,424	1/1:0,18:18:54:623,54,0	1/1:0,26:26:78:774,78,0
20_10037292_10066351	3529	.	T	A	154.29	.	AC=1;AF=0.167;AN=6;BaseQRankSum=-5.440e-01;DP=104;ExcessHet=0.0000;FS=1.871;MLEAC=1;MLEAF=0.167;MQ=60.00;MQRankSum=0.00;QD=7.71;ReadPosRankSum=-1.158e+00;SOR=1.034	GT:AD:DP:GQ:PL	0/0:44,0:44:99:0,112,1347	0/1:12,8:20:99:163,0,328	0/0:39,0:39:99:0,105,1194
```

Vous avez maintenant un workflow d'appel conjoint de variants automatisé et entièrement reproductible !

!!! note

    Gardez à l'esprit que les fichiers de données que nous vous avons fournis ne couvrent qu'une toute petite portion du chromosome 20.
    La taille réelle d'un ensemble d'appels de variants se compterait en millions de variants.
    C'est pourquoi nous n'utilisons que de minuscules sous-ensembles de données à des fins de formation !

### À retenir

Vous savez comment utiliser certains opérateurs courants ainsi que des closures Groovy pour contrôler le flux de données dans votre workflow.

### Et ensuite ?

Célébrez votre succès et prenez une pause bien méritée.

Dans la prochaine partie de ce cours, vous apprendrez à modulariser votre workflow en extrayant les définitions de processus dans des modules réutilisables.
