# Partie 3 : Variant calling joint sur une cohorte

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Traduction assistée par IA - [en savoir plus et suggérer des améliorations](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

Dans la Partie 2, vous avez construit un pipeline de variant calling par échantillon qui traitait les données de chaque échantillon de manière indépendante.
Nous allons maintenant l'étendre pour implémenter le variant calling joint, comme décrit dans la [Partie 1](01_method.md).

## Exercice

Dans cette partie du cours, nous allons étendre le workflow pour effectuer les opérations suivantes :

<figure class="excalidraw">
--8<-- "docs/en/docs/nf4_science/genomics/img/hello-gatk-2.svg"
</figure>

1. Générer un fichier d'index pour chaque fichier BAM d'entrée en utilisant Samtools
2. Exécuter GATK HaplotypeCaller sur chaque fichier BAM d'entrée pour générer un GVCF des appels de variants génomiques par échantillon
3. Collecter tous les GVCF et les combiner dans un data store GenomicsDB
4. Exécuter le génotypage joint sur le data store GVCF combiné pour produire un VCF au niveau de la cohorte

Cette partie se construit directement sur le workflow produit par la Partie 2.

??? info "Comment commencer à partir de cette section"

    Cette section du cours suppose que vous avez terminé la [Partie 2 : Variant calling par échantillon](./02_per_sample_variant_calling.md) et que vous disposez d'un pipeline `genomics.nf` fonctionnel.

    Si vous n'avez pas terminé la Partie 2 ou si vous souhaitez repartir de zéro pour cette partie, vous pouvez utiliser la solution de la Partie 2 comme point de départ.
    Exécutez ces commandes depuis le répertoire `nf4-science/genomics/` :

    ```bash
    cp solutions/part2/genomics-2.nf genomics.nf
    cp solutions/part2/nextflow.config .
    cp solutions/part2/modules/* modules/
    ```

    Cela vous donne un workflow complet de variant calling par échantillon.
    Vous pouvez tester qu'il s'exécute avec succès en lançant la commande suivante :

    ```bash
    nextflow run genomics.nf -profile test
    ```

## Plan de la leçon

Nous avons décomposé cela en deux étapes :

1. **Modifier l'étape de variant calling par échantillon pour produire un GVCF.**
   Cela couvre la mise à jour des commandes et des sorties des processus.
2. **Ajouter une étape de génotypage joint qui combine et génotype les GVCF par échantillon.**
   Cela introduit l'opérateur `collect()`, les closures Groovy pour la construction de ligne de commande et les processus multi-commandes.

!!! note

     Assurez-vous d'être dans le bon répertoire de travail :
     `cd /workspaces/training/nf4-science/genomics`

---

## 1. Modifier l'étape de variant calling par échantillon pour produire un GVCF

Le pipeline de la Partie 2 produit des fichiers VCF, mais le variant calling joint nécessite des fichiers GVCF.
Nous devons activer le mode de variant calling GVCF et mettre à jour l'extension du fichier de sortie.

Rappelez-vous la commande de variant calling GVCF de la [Partie 1](01_method.md) :

```bash
gatk HaplotypeCaller \
        -R /data/ref/ref.fasta \
        -I /data/bam/reads_mother.bam \
        -O reads_mother.g.vcf \
        -L /data/ref/intervals.bed \
        -ERC GVCF
```

Par rapport à la commande HaplotypeCaller de base que nous avons encapsulée dans la Partie 2, les différences sont le paramètre `-ERC GVCF` et l'extension de sortie `.g.vcf`.

### 1.1. Indiquer à HaplotypeCaller d'émettre un GVCF et mettre à jour l'extension de sortie

Ouvrez le fichier de module `modules/gatk_haplotypecaller.nf` pour effectuer deux modifications :

- Ajoutez le paramètre `-ERC GVCF` à la commande GATK HaplotypeCaller ;
- Mettez à jour le chemin du fichier de sortie pour utiliser l'extension `.g.vcf` correspondante, selon la convention GATK.

Assurez-vous d'ajouter un antislash (`\`) à la fin de la ligne précédente lorsque vous ajoutez `-ERC GVCF`.

=== "Après"

    ```groovy title="modules/gatk_haplotypecaller.nf" linenums="22" hl_lines="5 7"
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

    ```groovy title="modules/gatk_haplotypecaller.nf" linenums="22" hl_lines="5"
        """
        gatk HaplotypeCaller \
            -R ${ref_fasta} \
            -I ${input_bam} \
            -O ${input_bam}.vcf \
            -L ${interval_list}
        """
    ```

Nous devons également mettre à jour le bloc de sortie pour correspondre à la nouvelle extension de fichier.
Puisque nous avons changé la sortie de la commande de `.vcf` à `.g.vcf`, le bloc `output:` du processus doit refléter le même changement.

### 1.2. Mettre à jour l'extension du fichier de sortie dans le bloc des sorties du processus

=== "Après"

    ```groovy title="modules/gatk_haplotypecaller.nf" linenums="17" hl_lines="2 3"
        output:
        path "${input_bam}.g.vcf"     , emit: vcf
        path "${input_bam}.g.vcf.idx" , emit: idx
    ```

=== "Avant"

    ```groovy title="modules/gatk_haplotypecaller.nf" linenums="17" hl_lines="2 3"
        output:
        path "${input_bam}.vcf"     , emit: vcf
        path "${input_bam}.vcf.idx" , emit: idx
    ```

Nous devons également mettre à jour la configuration de publication et de sortie du workflow pour refléter les nouvelles sorties GVCF.

### 1.3. Mettre à jour les cibles de publication pour les nouvelles sorties GVCF

Puisque nous produisons maintenant des GVCF au lieu de VCF, nous devons mettre à jour la section `publish:` du workflow pour utiliser des noms plus descriptifs.
Nous allons également organiser les fichiers GVCF dans leur propre sous-répertoire pour plus de clarté.

=== "Après"

    ```groovy title="genomics.nf" linenums="47" hl_lines="3 4"
        publish:
        indexed_bam = SAMTOOLS_INDEX.out
        gvcf = GATK_HAPLOTYPECALLER.out.vcf
        gvcf_idx = GATK_HAPLOTYPECALLER.out.idx
    ```

=== "Avant"

    ```groovy title="genomics.nf" linenums="47"
        publish:
        indexed_bam = SAMTOOLS_INDEX.out
        vcf = GATK_HAPLOTYPECALLER.out.vcf
        vcf_idx = GATK_HAPLOTYPECALLER.out.idx
    ```

Maintenant, mettez à jour le bloc de sortie pour correspondre.

### 1.4. Mettre à jour le bloc de sortie pour la nouvelle structure de répertoires

Nous devons également mettre à jour le bloc `output` pour placer les fichiers GVCF dans un sous-répertoire `gvcf`.

=== "Après"

    ```groovy title="genomics.nf" linenums="53" hl_lines="3 5 6 8 9"
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

    ```groovy title="genomics.nf" linenums="53"
    output {
        indexed_bam {
            path 'bam'
        }
        vcf {
            path 'vcf'
        }
        vcf_idx {
            path 'vcf'
        }
    }
    ```

Avec le module, les cibles de publication et le bloc de sortie tous mis à jour, nous pouvons tester les modifications.

### 1.5. Exécuter le pipeline

Exécutez le workflow pour vérifier que les modifications fonctionnent.

```bash
nextflow run genomics.nf
```

??? success "Sortie de la commande"

    ```console
    N E X T F L O W   ~  version 25.10.2

    ┃ Launching `genomics.nf` [nostalgic_franklin] DSL2 - revision: f2c0a93c6a

    executor >  local (6)
    [cc/fbc705] SAMTOOLS_INDEX (3)       | 3 of 3 ✔
    [27/0d7eb9] GATK_HAPLOTYPECALLER (2) | 3 of 3 ✔
    ```

La sortie Nextflow semble identique à celle d'avant, mais les fichiers `.g.vcf` et leurs fichiers d'index sont maintenant organisés dans des sous-répertoires.

??? abstract "Contenu du répertoire (liens symboliques raccourcis)"

    ```console
    results/
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

Si vous ouvrez l'un des fichiers GVCF et que vous le parcourez, vous pouvez vérifier que GATK HaplotypeCaller a produit des fichiers GVCF comme demandé.

### À retenir

Lorsque vous modifiez le nom de fichier de sortie d'une commande d'outil, le bloc `output:` du processus et la configuration de publication/sortie doivent être mis à jour en conséquence.

### Et ensuite ?

Apprenez à collecter le contenu d'un canal et à les transmettre au processus suivant comme une seule entrée.

---

## 2. Ajouter une étape de génotypage joint

Nous devons maintenant collecter les GVCF par échantillon, les combiner dans un data store GenomicsDB et exécuter le génotypage joint pour produire un VCF au niveau de la cohorte.
Comme expliqué dans la [Partie 1](01_method.md), il s'agit d'une opération à deux outils : GenomicsDBImport combine les GVCF, puis GenotypeGVCFs produit les appels de variants finaux.
Nous allons encapsuler les deux outils dans un seul processus appelé `GATK_JOINTGENOTYPING`.

Rappelez-vous les deux commandes de la [Partie 1](01_method.md) :

```bash
gatk GenomicsDBImport \
    -V reads_mother.g.vcf \
    -V reads_father.g.vcf \
    -V reads_son.g.vcf \
    -L /data/ref/intervals.bed \
    --genomicsdb-workspace-path family_trio_gdb
```

```bash
gatk GenotypeGVCFs \
    -R /data/ref/ref.fasta \
    -V gendb://family_trio_gdb \
    -O family_trio.vcf
```

La première commande prend les GVCF par échantillon et un fichier d'intervalles, et produit un data store GenomicsDB.
La seconde prend ce data store, un génome de référence, et produit le VCF final au niveau de la cohorte.
L'URI du conteneur est la même que pour HaplotypeCaller : `community.wave.seqera.io/library/gatk4:4.5.0.0--730ee8817e436867`.

### 2.1. Configurer les entrées

Le processus de génotypage joint nécessite deux types d'entrées que nous n'avons pas encore : un nom de cohorte arbitraire et les sorties GVCF collectées de tous les échantillons regroupées ensemble.

#### 2.1.1. Ajouter un paramètre `cohort_name`

Nous devons fournir un nom arbitraire pour la cohorte.
Plus tard dans la série de formations, vous apprendrez à utiliser les métadonnées d'échantillon pour ce genre de chose, mais pour l'instant nous déclarons simplement un paramètre CLI en utilisant `params` et nous lui donnons une valeur par défaut pour plus de commodité.

=== "Après"

    ```groovy title="genomics.nf" linenums="14" hl_lines="3-4"
        intervals: Path = "${projectDir}/data/ref/intervals.bed"

        // Nom de base pour le fichier de sortie final
        cohort_name: String = "family_trio"
    }
    ```

=== "Avant"

    ```groovy title="genomics.nf" linenums="14"
        intervals: Path = "${projectDir}/data/ref/intervals.bed"
    }
    ```

#### 2.1.2. Rassembler les sorties de HaplotypeCaller entre les échantillons

Si nous devions connecter directement le canal de sortie de `GATK_HAPLOTYPECALLER` au nouveau processus, Nextflow appellerait le processus sur chaque GVCF d'échantillon séparément.
Nous voulons regrouper tous les trois GVCF (et leurs fichiers d'index) pour que Nextflow les transmette tous ensemble à un seul appel de processus.

Nous pouvons faire cela en utilisant l'opérateur de canal `collect()`.
Ajoutez les lignes suivantes au corps du `workflow`, juste après l'appel à GATK_HAPLOTYPECALLER :

=== "Après"

    ```groovy title="genomics.nf" hl_lines="4-6"
            intervals_file
        )

        // Collecter les sorties de variant calling entre les échantillons
        all_gvcfs_ch = GATK_HAPLOTYPECALLER.out.vcf.collect()
        all_idxs_ch = GATK_HAPLOTYPECALLER.out.idx.collect()
    ```

=== "Avant"

    ```groovy title="genomics.nf"
            intervals_file
        )
    ```

Décortiquons cela :

1. Nous prenons le canal de sortie de `GATK_HAPLOTYPECALLER` en utilisant la propriété `.out`.
2. Parce que nous avons nommé les sorties en utilisant `emit:` dans la section 1, nous pouvons sélectionner les GVCF avec `.vcf` et les fichiers d'index avec `.idx`. Sans sorties nommées, nous devrions utiliser `.out[0]` et `.out[1]`.
3. L'opérateur `collect()` regroupe tous les fichiers en un seul élément, donc `all_gvcfs_ch` contient les trois GVCF ensemble, et `all_idxs_ch` contient les trois fichiers d'index ensemble.

Nous pouvons collecter les GVCF et leurs fichiers d'index séparément (au lieu de les garder ensemble dans des tuples) parce que Nextflow mettra en place tous les fichiers d'entrée ensemble pour l'exécution, donc les fichiers d'index seront présents aux côtés des GVCF.

!!! tip

    Vous pouvez utiliser l'opérateur `view()` pour inspecter le contenu des canaux avant et après l'application des opérateurs de canal.

### 2.2. Écrire le processus de génotypage joint et l'appeler dans le workflow

En suivant le même modèle que nous avons utilisé dans la Partie 2, nous allons écrire la définition du processus dans un fichier de module, l'importer dans le workflow et l'appeler sur les entrées que nous venons de préparer.

#### 2.2.1. Construire une chaîne pour donner à chaque GVCF un argument `-V`

Avant de commencer à remplir la définition du processus, il y a une chose à régler.
La commande GenomicsDBImport attend un argument `-V` séparé pour chaque fichier GVCF, comme ceci :

```bash
gatk GenomicsDBImport \
    -V reads_mother.bam.g.vcf \
    -V reads_father.bam.g.vcf \
    -V reads_son.bam.g.vcf \
    ...
```

Si nous devions écrire `-V ${all_gvcfs_ch}`, Nextflow concaténerait simplement les noms de fichiers et cette partie de la commande ressemblerait à ceci :

```groovy
-V reads_mother.bam.g.vcf reads_father.bam.g.vcf reads_son.bam.g.vcf
```

Mais nous avons besoin que la chaîne ressemble à ceci :

```groovy
-V reads_mother.bam.g.vcf -V reads_father.bam.g.vcf -V reads_son.bam.g.vcf
```

Surtout, nous devons construire cette chaîne de manière dynamique à partir de tous les fichiers présents dans le canal collecté.
Nextflow (via Groovy) fournit un moyen concis de le faire :

```groovy
def gvcfs_line = all_gvcfs.collect { gvcf -> "-V ${gvcf}" }.join(' ')
```

Décortiquons cela :

1. `all_gvcfs.collect { gvcf -> "-V ${gvcf}" }` itère sur chaque chemin de fichier et préfixe `-V ` devant lui, produisant `["-V A.g.vcf", "-V B.g.vcf", "-V C.g.vcf"]`.
2. `.join(' ')` les concatène avec des espaces : `"-V A.g.vcf -V B.g.vcf -V C.g.vcf"`.
3. Le résultat est assigné à une variable locale `gvcfs_line` (définie avec `def`), que nous pouvons interpoler dans le template de commande.

Cette ligne va dans le bloc `script:` du processus, avant le template de commande.
Vous pouvez placer du code Groovy arbitraire entre `script:` et le `"""` d'ouverture du template de commande.

Ensuite, vous pourrez vous référer à toute cette chaîne comme `gvcfs_line` dans le bloc `script:` du processus.

#### 2.2.2. Compléter le module pour le processus de génotypage joint

Nous pouvons maintenant aborder l'écriture du processus complet.

Ouvrez `modules/gatk_jointgenotyping.nf` et examinez le squelette de la définition du processus.

Continuez et complétez la définition du processus en utilisant les informations fournies ci-dessus, puis vérifiez votre travail par rapport à la solution dans l'onglet "Après" ci-dessous.

=== "Avant"

    ```groovy title="modules/gatk_jointgenotyping.nf" linenums="1"
    #!/usr/bin/env nextflow

    /*
     * Combiner les GVCF dans un data store GenomicsDB et exécuter le génotypage joint pour produire des appels au niveau de la cohorte
     */
    process GATK_JOINTGENOTYPING {

        container

        input:

        output:

        script:
        """

        """
    }
    ```

=== "Après"

    ```groovy title="modules/gatk_jointgenotyping.nf" linenums="1" hl_lines="8 11-17 20-21 24-25 29-33"
    #!/usr/bin/env nextflow

    /*
     * Combiner les GVCF dans un data store GenomicsDB et exécuter le génotypage joint pour produire des appels au niveau de la cohorte
     */
    process GATK_JOINTGENOTYPING {

        container "community.wave.seqera.io/library/gatk4:4.5.0.0--730ee8817e436867"

        input:
        path all_gvcfs
        path all_idxs
        path interval_list
        val cohort_name
        path ref_fasta
        path ref_index
        path ref_dict

        output:
        path "${cohort_name}.joint.vcf"     , emit: vcf
        path "${cohort_name}.joint.vcf.idx" , emit: idx

        script:
        def gvcfs_line = all_gvcfs.collect { gvcf -> "-V ${gvcf}" }.join(' ')
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
    }
    ```

Il y a plusieurs choses à souligner ici.

Comme précédemment, plusieurs entrées sont listées même si les commandes ne les référencent pas directement : `all_idxs`, `ref_index` et `ref_dict`.
Les lister garantit que Nextflow met en place ces fichiers dans le répertoire de travail aux côtés des fichiers qui apparaissent dans les commandes, ce que GATK s'attend à trouver selon les conventions de nommage.

La variable `gvcfs_line` utilise la closure Groovy décrite ci-dessus pour construire les arguments `-V` pour GenomicsDBImport.

Ce processus exécute deux commandes en série, comme vous le feriez dans le terminal.
GenomicsDBImport combine les GVCF par échantillon dans un data store, puis GenotypeGVCFs lit ce data store et produit le VCF final au niveau de la cohorte.
Le data store GenomicsDB (`${cohort_name}_gdb`) est un artefact intermédiaire utilisé uniquement au sein du processus ; il n'apparaît pas dans le bloc de sortie.

Une fois que vous avez terminé cela, le processus est prêt à être utilisé.
Pour l'utiliser dans le workflow, vous devrez importer le module et ajouter un appel de processus.

#### 2.2.3. Importer le module

Ajoutez l'instruction d'importation à `genomics.nf`, sous les instructions d'importation existantes :

=== "Après"

    ```groovy title="genomics.nf" linenums="21" hl_lines="3"
    include { SAMTOOLS_INDEX } from './modules/samtools_index.nf'
    include { GATK_HAPLOTYPECALLER } from './modules/gatk_haplotypecaller.nf'
    include { GATK_JOINTGENOTYPING } from './modules/gatk_jointgenotyping.nf'
    ```

=== "Avant"

    ```groovy title="genomics.nf" linenums="21"
    include { SAMTOOLS_INDEX } from './modules/samtools_index.nf'
    include { GATK_HAPLOTYPECALLER } from './modules/gatk_haplotypecaller.nf'
    ```

Le processus est maintenant disponible dans la portée du workflow.

#### 2.2.4. Ajouter l'appel de processus

Ajoutez l'appel à `GATK_JOINTGENOTYPING` dans le corps du workflow, après les lignes `collect()` :

=== "Après"

    ```groovy title="genomics.nf" hl_lines="3-12"
        all_idxs_ch = GATK_HAPLOTYPECALLER.out.idx.collect()

        // Combiner les GVCF dans un data store GenomicsDB et appliquer le génotypage joint
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

    ```groovy title="genomics.nf"
        all_idxs_ch = GATK_HAPLOTYPECALLER.out.idx.collect()
    ```

Le processus est maintenant entièrement câblé.
Ensuite, nous configurons la façon dont les sorties sont publiées.

### 2.3. Configurer la gestion des sorties

Nous devons publier les sorties du VCF joint.
Ajoutez des cibles de publication et des entrées de bloc de sortie pour les résultats du génotypage joint.

#### 2.3.1. Ajouter des cibles de publication pour le VCF joint

Ajoutez le VCF joint et son index à la section `publish:` du workflow :

=== "Après"

    ```groovy title="genomics.nf" hl_lines="5-6"
        publish:
        indexed_bam = SAMTOOLS_INDEX.out
        gvcf = GATK_HAPLOTYPECALLER.out.vcf
        gvcf_idx = GATK_HAPLOTYPECALLER.out.idx
        joint_vcf = GATK_JOINTGENOTYPING.out.vcf
        joint_vcf_idx = GATK_JOINTGENOTYPING.out.idx
    ```

=== "Avant"

    ```groovy title="genomics.nf"
        publish:
        indexed_bam = SAMTOOLS_INDEX.out
        gvcf = GATK_HAPLOTYPECALLER.out.vcf
        gvcf_idx = GATK_HAPLOTYPECALLER.out.idx
    ```

Maintenant, mettez à jour le bloc de sortie pour correspondre.

#### 2.3.2. Ajouter des entrées de bloc de sortie pour le VCF joint

Ajoutez des entrées pour les fichiers VCF joints.
Nous allons les placer à la racine du répertoire de résultats puisqu'il s'agit de la sortie finale.

=== "Après"

    ```groovy title="genomics.nf" hl_lines="11-16"
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
        joint_vcf {
            path '.'
        }
        joint_vcf_idx {
            path '.'
        }
    }
    ```

=== "Avant"

    ```groovy title="genomics.nf"
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

Avec le processus, les cibles de publication et le bloc de sortie tous en place, nous pouvons tester le workflow complet.

### 2.4. Exécuter le workflow

Exécutez le workflow pour vérifier que tout fonctionne.

```bash
nextflow run genomics.nf -resume
```

??? success "Sortie de la commande"

    ```console
    N E X T F L O W   ~  version 25.10.2

    ┃ Launching `genomics.nf` [crazy_marconi] DSL2 - revision: 5da9afc841

    executor >  local (1)
    [9a/c7a873] SAMTOOLS_INDEX (2)       | 3 of 3, cached: 3 ✔
    [e4/4ed55e] GATK_HAPLOTYPECALLER (2) | 3 of 3, cached: 3 ✔
    [a6/7cc8ed] GATK_JOINTGENOTYPING     | 1 of 1 ✔
    ```

Les deux premières étapes sont mises en cache depuis l'exécution précédente, et la nouvelle étape `GATK_JOINTGENOTYPING` s'exécute une fois sur les entrées collectées des trois échantillons.
Le fichier de sortie final, `family_trio.joint.vcf` (et son index), se trouve dans le répertoire de résultats.

??? abstract "Contenu du répertoire (liens symboliques raccourcis)"

    ```console
    results/
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

Si vous ouvrez le fichier VCF joint, vous pouvez vérifier que le workflow a produit les appels de variants attendus.

```console title="family_trio.joint.vcf" linenums="40"
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	reads_father	reads_mother	reads_son
20_10037292_10066351	3480	.	C	CT	1625.89	.	AC=5;AF=0.833;AN=6;BaseQRankSum=0.220;DP=85;ExcessHet=0.0000;FS=2.476;MLEAC=5;MLEAF=0.833;MQ=60.00;MQRankSum=0.00;QD=21.68;ReadPosRankSum=-1.147e+00;SOR=0.487	GT:AD:DP:GQ:PL	0/1:15,16:31:99:367,0,375	1/1:0,18:18:54:517,54,0	1/1:0,26:26:78:756,78,0
20_10037292_10066351	3520	.	AT	A	1678.89	.	AC=5;AF=0.833;AN=6;BaseQRankSum=1.03;DP=80;ExcessHet=0.0000;FS=2.290;MLEAC=5;MLEAF=0.833;MQ=60.00;MQRankSum=0.00;QD=22.39;ReadPosRankSum=0.701;SOR=0.730	GT:AD:DP:GQ:PL	0/1:18,13:31:99:296,0,424	1/1:0,18:18:54:623,54,0	1/1:0,26:26:78:774,78,0
20_10037292_10066351	3529	.	T	A	154.29	.	AC=1;AF=0.167;AN=6;BaseQRankSum=-5.440e-01;DP=104;ExcessHet=0.0000;FS=1.871;MLEAC=1;MLEAF=0.167;MQ=60.00;MQRankSum=0.00;QD=7.71;ReadPosRankSum=-1.158e+00;SOR=1.034	GT:AD:DP:GQ:PL	0/0:44,0:44:99:0,112,1347	0/1:12,8:20:99:163,0,328	0/0:39,0:39:99:0,105,1194
```

Vous disposez maintenant d'un workflow de variant calling joint automatisé et entièrement reproductible !

!!! note

    Gardez à l'esprit que les fichiers de données que nous vous avons fournis ne couvrent qu'une toute petite portion du chromosome 20.
    La taille réelle d'un ensemble d'appels de variants se compterait en millions de variants.
    C'est pourquoi nous n'utilisons que de minuscules sous-ensembles de données à des fins de formation !

### À retenir

Vous savez comment collecter les sorties d'un canal et les regrouper comme une seule entrée pour un autre processus.
Vous savez également comment construire une ligne de commande en utilisant des closures Groovy, et comment exécuter plusieurs commandes dans un seul processus.

### Et ensuite ?

Félicitations ! Vous avez terminé le cours Nextflow pour la génomique.

Rendez-vous au [résumé final du cours](./next_steps.md) pour revoir ce que vous avez appris et découvrir la suite.
