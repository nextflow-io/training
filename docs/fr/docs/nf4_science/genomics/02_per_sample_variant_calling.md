# Partie 2 : Appel de variants par échantillon

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Traduction assistée par IA - [en savoir plus et suggérer des améliorations](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

Dans la Partie 1, vous avez testé les commandes Samtools et GATK manuellement dans leurs conteneurs respectifs.
Nous allons maintenant encapsuler ces mêmes commandes dans un workflow Nextflow.

## Objectif

Dans cette partie du cours, nous allons développer un workflow qui effectue les opérations suivantes :

1. Générer un fichier d'index pour chaque fichier BAM d'entrée en utilisant [Samtools](https://www.htslib.org/)
2. Exécuter GATK HaplotypeCaller sur chaque fichier BAM d'entrée pour générer des appels de variants par échantillon au format VCF (Variant Call Format)

<figure class="excalidraw">
--8<-- "docs/en/docs/nf4_science/genomics/img/hello-gatk-1.svg"
</figure>

Cela reproduit les étapes de la Partie 1, où vous avez exécuté ces commandes manuellement dans leurs conteneurs.

Comme point de départ, nous vous fournissons un fichier de workflow, `genomics.nf`, qui décrit les principales parties du workflow, ainsi que deux fichiers de modules, samtools_index.nf et gatk_haplotypecaller.nf, qui décrivent la structure des modules.
Ces fichiers ne sont pas fonctionnels ; leur objectif est simplement de servir de squelettes que vous compléterez avec les parties intéressantes du code.

## Plan de la leçon

Afin de rendre le processus de développement plus pédagogique, nous avons divisé cela en quatre étapes :

1. **Écrire un workflow à une seule étape qui exécute Samtools index sur un fichier BAM.**
   Cela couvre la création d'un module, son importation et son appel dans un workflow.
2. **Ajouter un second processus pour exécuter GATK HaplotypeCaller sur le fichier BAM indexé.**
   Cela introduit le chaînage des sorties de processus vers les entrées et la gestion des fichiers accessoires.
3. **Adapter le workflow pour qu'il s'exécute sur un lot d'échantillons.**
   Cela couvre l'exécution parallèle et introduit les tuples pour garder les fichiers associés ensemble.
4. **Faire en sorte que le workflow accepte un fichier texte contenant un lot de fichiers d'entrée.**
   Cela démontre un modèle courant pour fournir des entrées en masse.

Chaque étape se concentre sur un aspect spécifique du développement de workflow.

---

## 1. Écrire un workflow à une seule étape qui exécute Samtools index sur un fichier BAM

Cette première étape se concentre sur les bases : charger un fichier BAM et générer un index pour celui-ci.

Rappelez-vous la commande `samtools index` de la [Partie 1](01_method.md) :

```bash
samtools index '<input_bam>'
```

La commande prend un fichier BAM en entrée et produit un fichier d'index `.bai` à côté de celui-ci.
L'URI du conteneur était `community.wave.seqera.io/library/samtools:1.20--b5dfbd93de237464`.

Nous allons prendre ces informations et les encapsuler dans Nextflow en trois étapes :

1. Configurer l'entrée
2. Écrire le processus d'indexation et l'appeler dans le workflow
3. Configurer la gestion de la sortie

### 1.1. Configurer l'entrée

Nous devons déclarer un paramètre d'entrée, créer un profil de test pour fournir une valeur par défaut pratique, et créer un canal d'entrée.

#### 1.1.1. Ajouter une déclaration de paramètre d'entrée

Dans le fichier de workflow principal `genomics.nf`, sous la section `Pipeline parameters`, déclarez un paramètre CLI appelé `reads_bam`.

=== "Après"

    ```groovy title="genomics.nf" linenums="5" hl_lines="4-7"
    /*
     * Pipeline parameters
     */
    params {
        // Entrée principale
        reads_bam: Path
    }
    ```

=== "Avant"

    ```groovy title="genomics.nf" linenums="5"
    /*
     * Pipeline parameters
     */

    // Entrée principale
    ```

Cela configure le paramètre CLI, mais nous ne voulons pas taper le chemin du fichier à chaque fois que nous exécutons le workflow pendant le développement.
Il existe plusieurs options pour fournir une valeur par défaut ; ici nous utilisons un profil de test.

#### 1.1.2. Créer un profil de test avec une valeur par défaut dans `nextflow.config`

Un profil de test fournit des valeurs par défaut pratiques pour essayer un workflow sans spécifier d'entrées sur la ligne de commande.
C'est une convention courante dans l'écosystème Nextflow (voir [Hello Config](../../hello_nextflow/06_hello_config.md) pour plus de détails).

Ajoutez un bloc `profiles` à `nextflow.config` avec un profil `test` qui définit le paramètre `reads_bam` sur l'un des fichiers BAM de test.

=== "Après"

    ```groovy title="nextflow.config" linenums="1" hl_lines="3-7"
    docker.enabled = true

    profiles {
        test {
            params.reads_bam = "${projectDir}/data/bam/reads_mother.bam"
        }
    }
    ```

=== "Avant"

    ```groovy title="nextflow.config" linenums="1"
    docker.enabled = true
    ```

Ici, nous utilisons `${projectDir}`, une variable Nextflow intégrée qui pointe vers le répertoire où se trouve le script de workflow.
Cela facilite la référence aux fichiers de données et autres ressources sans coder en dur des chemins absolus.

#### 1.1.3. Configurer le canal d'entrée

Dans le bloc workflow, créez un canal d'entrée à partir de la valeur du paramètre en utilisant la fabrique de canaux `.fromPath` (comme utilisée dans [Hello Channels](../../hello_nextflow/02_hello_channels.md)).

=== "Après"

    ```groovy title="genomics.nf" linenums="13" hl_lines="4-5"
    workflow {

        main:
        // Créer le canal d'entrée (fichier unique via paramètre CLI)
        reads_ch = channel.fromPath(params.reads_bam)
    ```

=== "Avant"

    ```groovy title="genomics.nf" linenums="13"
    workflow {

        main:
        // Créer le canal d'entrée
    ```

Maintenant nous devons créer le processus pour exécuter l'indexation sur cette entrée.

### 1.2. Écrire le processus d'indexation et l'appeler dans le workflow

Nous devons écrire la définition du processus dans le fichier de module, l'importer dans le workflow en utilisant une instruction include, et l'appeler sur l'entrée.

#### 1.2.1. Compléter le module pour le processus d'indexation

Ouvrez `modules/samtools_index.nf` et examinez le squelette de la définition du processus.
Vous devriez reconnaître les principaux éléments structurels ; sinon, envisagez de lire [Hello Nextflow](../../hello_nextflow/01_hello_world.md) pour un rappel.

Complétez la définition du processus par vous-même en utilisant les informations fournies ci-dessus, puis vérifiez votre travail par rapport à la solution dans l'onglet "Après" ci-dessous.

=== "Avant"

    ```groovy title="modules/samtools_index.nf" linenums="1"
    #!/usr/bin/env nextflow

    /*
     * Générer un fichier d'index BAM
     */
    process SAMTOOLS_INDEX {

        container

        input:

        output:

        script:
        """

        """
    }
    ```

=== "Après"

    ```groovy title="modules/samtools_index.nf" linenums="1" hl_lines="8 11 14 18"
    #!/usr/bin/env nextflow

    /*
     * Générer un fichier d'index BAM
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

Une fois que vous avez terminé cela, le processus est complet.
Pour l'utiliser dans le workflow, vous devrez importer le module et ajouter un appel de processus.

#### 1.2.2. Inclure le module

Dans `genomics.nf`, ajoutez une instruction `include` pour rendre le processus disponible au workflow :

=== "Après"

    ```groovy title="genomics.nf" linenums="3" hl_lines="2"
    // Instructions INCLUDE de modules
    include { SAMTOOLS_INDEX } from './modules/samtools_index.nf'
    ```

=== "Avant"

    ```groovy title="genomics.nf" linenums="3"
    // Instructions INCLUDE de modules
    ```

Le processus est maintenant disponible dans la portée du workflow.

#### 1.2.3. Appeler le processus d'indexation sur l'entrée

Maintenant, ajoutons un appel à `SAMTOOLS_INDEX` dans le bloc workflow, en passant le canal d'entrée comme argument.

=== "Après"

    ```groovy title="genomics.nf" linenums="14" hl_lines="7-8"
    workflow {

        main:
        // Créer le canal d'entrée (fichier unique via paramètre CLI)
        reads_ch = channel.fromPath(params.reads_bam)

        // Créer le fichier d'index pour le fichier BAM d'entrée
        SAMTOOLS_INDEX(reads_ch)
    ```

=== "Avant"

    ```groovy title="genomics.nf" linenums="14"
    workflow {

        main:
        // Créer le canal d'entrée (fichier unique via paramètre CLI)
        reads_ch = channel.fromPath(params.reads_bam)

        // Appeler les processus
    ```

Le workflow charge maintenant l'entrée et exécute le processus d'indexation dessus.
Ensuite, nous devons configurer comment la sortie est publiée.

### 1.3. Configurer la gestion de la sortie

Nous devons déclarer quelles sorties de processus publier et spécifier où elles doivent aller.

#### 1.3.1. Déclarer une sortie dans la section `publish:`

La section `publish:` à l'intérieur du bloc workflow déclare quelles sorties de processus doivent être publiées.
Assignez la sortie de `SAMTOOLS_INDEX` à une cible nommée appelée `bam_index`.

=== "Après"

    ```groovy title="genomics.nf" linenums="22" hl_lines="2"
        publish:
        bam_index = SAMTOOLS_INDEX.out
    }
    ```

=== "Avant"

    ```groovy title="genomics.nf" linenums="22"
        publish:
        // Déclarer les sorties à publier
    }
    ```

Maintenant nous devons dire à Nextflow où mettre la sortie publiée.

#### 1.3.2. Configurer la cible de sortie dans le bloc `output {}`

Le bloc `output {}` se situe en dehors du workflow et spécifie où chaque cible nommée est publiée.
Ajoutons une cible pour `bam_index` qui publie dans un sous-répertoire `bam/`.

=== "Après"

    ```groovy title="genomics.nf" linenums="26" hl_lines="2-4"
    output {
        bam_index {
            path 'bam'
        }
    }
    ```

=== "Avant"

    ```groovy title="genomics.nf" linenums="26"
    output {
        // Configurer les cibles de publication
    }
    ```

!!! note "Note"

    Par défaut, Nextflow publie les fichiers de sortie sous forme de liens symboliques, ce qui évite une duplication inutile.
    Même si les fichiers de données que nous utilisons ici sont très petits, en génomique ils peuvent devenir très volumineux.
    Les liens symboliques seront rompus lorsque vous nettoierez votre répertoire `work`, donc pour les workflows de production vous voudrez peut-être remplacer le mode de publication par défaut par `'copy'`.

### 1.4. Exécuter le workflow

À ce stade, nous avons un workflow d'indexation en une étape qui devrait être entièrement fonctionnel. Testons qu'il fonctionne !

Nous pouvons l'exécuter avec `-profile test` pour utiliser la valeur par défaut configurée dans le profil de test et éviter d'avoir à écrire le chemin sur la ligne de commande.

```bash
nextflow run genomics.nf -profile test
```

??? success "Sortie de la commande"

    ```console
    N E X T F L O W   ~  version 25.10.2

    ┃ Launching `genomics.nf` [reverent_sinoussi] DSL2 - revision: 41d43ad7fe

    executor >  local (1)
    [2a/e69536] SAMTOOLS_INDEX (1) | 1 of 1 ✔
    ```

Vous pouvez vérifier que le fichier d'index a été généré correctement en regardant dans le répertoire de travail ou dans le répertoire de résultats.

??? abstract "Contenu du répertoire de travail"

    ```console
    work/2a/e695367b2f60df09cf826b07192dc3
    ├── reads_mother.bam -> /workspaces/training/nf4-science/genomics/data/bam/reads_mother.bam
    └── reads_mother.bam.bai
    ```

??? abstract "Contenu du répertoire de résultats"

    ```console
    results/
    └── bam/
        └── reads_mother.bam.bai -> ...
    ```

Le voilà !

### À retenir

Vous savez comment créer un module contenant un processus, l'importer dans un workflow, l'appeler avec un canal d'entrée et publier les résultats.

### Et ensuite ?

Ajouter une seconde étape qui prend la sortie du processus d'indexation et l'utilise pour exécuter l'appel de variants.

---

## 2. Ajouter un second processus pour exécuter GATK HaplotypeCaller sur le fichier BAM indexé

Maintenant que nous avons un index pour notre fichier d'entrée, nous pouvons passer à la configuration de l'étape d'appel de variants.

Rappelez-vous la commande `gatk HaplotypeCaller` de la [Partie 1](01_method.md) :

```bash
gatk HaplotypeCaller \
        -R /data/ref/ref.fasta \
        -I /data/bam/reads_mother.bam \
        -O reads_mother.vcf \
        -L /data/ref/intervals.bed
```

La commande prend un fichier BAM (`-I`), un génome de référence (`-R`) et un fichier d'intervalles (`-L`), et produit un fichier VCF (`-O`) ainsi que son index.
L'outil s'attend également à ce que l'index BAM, l'index de référence et le dictionnaire de référence soient co-localisés avec leurs fichiers respectifs.
L'URI du conteneur était `community.wave.seqera.io/library/gatk4:4.5.0.0--730ee8817e436867`.

Nous suivons les trois mêmes étapes qu'auparavant :

1. Configurer les entrées
2. Écrire le processus d'appel de variants et l'appeler dans le workflow
3. Configurer la gestion de la sortie

### 2.1. Configurer les entrées

L'étape d'appel de variants nécessite plusieurs fichiers d'entrée supplémentaires.
Nous devons déclarer des paramètres pour eux, ajouter des valeurs par défaut au profil de test et créer des variables pour les charger.

#### 2.1.1. Ajouter des déclarations de paramètres pour les entrées accessoires

Puisque notre nouveau processus attend plusieurs fichiers supplémentaires à fournir, ajoutez des déclarations de paramètres pour eux dans `genomics.nf` sous la section `Pipeline parameters` :

=== "Après"

    ```groovy title="genomics.nf" linenums="9" hl_lines="5-9"
    params {
        // Entrée principale
        reads_bam: Path

        // Fichiers accessoires
        reference: Path
        reference_index: Path
        reference_dict: Path
        intervals: Path
    }
    ```

=== "Avant"

    ```groovy title="genomics.nf" linenums="9"
    params {
        // Entrée principale
        reads_bam: Path
    }
    ```

Comme auparavant, nous fournissons des valeurs par défaut via le profil de test plutôt qu'en ligne.

#### 2.1.2. Ajouter les valeurs par défaut des fichiers accessoires au profil de test

Tout comme nous l'avons fait pour `reads_bam` dans la section 1.1.2, ajoutez des valeurs par défaut pour les fichiers accessoires au profil de test dans `nextflow.config` :

=== "Après"

    ```groovy title="nextflow.config" linenums="4" hl_lines="3-6"
    test {
        params.reads_bam = "${projectDir}/data/bam/reads_mother.bam"
        params.reference = "${projectDir}/data/ref/ref.fasta"
        params.reference_index = "${projectDir}/data/ref/ref.fasta.fai"
        params.reference_dict = "${projectDir}/data/ref/ref.dict"
        params.intervals = "${projectDir}/data/ref/intervals.bed"
    }
    ```

=== "Avant"

    ```groovy title="nextflow.config" linenums="4"
    test {
        params.reads_bam = "${projectDir}/data/bam/reads_mother.bam"
    }
    ```

Maintenant nous devons créer des variables qui chargent ces chemins de fichiers pour une utilisation dans le workflow.

#### 2.1.3. Créer des variables pour les fichiers accessoires

Ajoutez des variables pour les chemins des fichiers accessoires à l'intérieur du bloc workflow :

=== "Après"

    ```groovy title="genomics.nf" linenums="21" hl_lines="7-11"
    workflow {

        main:
        // Créer le canal d'entrée (fichier unique via paramètre CLI)
        reads_ch = channel.fromPath(params.reads_bam)

        // Charger les chemins de fichiers pour les fichiers accessoires (référence et intervalles)
        ref_file        = file(params.reference)
        ref_index_file  = file(params.reference_index)
        ref_dict_file   = file(params.reference_dict)
        intervals_file  = file(params.intervals)

        // Créer le fichier d'index pour le fichier BAM d'entrée
        SAMTOOLS_INDEX(reads_ch)
    ```

=== "Avant"

    ```groovy title="genomics.nf" linenums="21"
    workflow {

        main:
        // Créer le canal d'entrée (fichier unique via paramètre CLI)
        reads_ch = channel.fromPath(params.reads_bam)

        // Créer le fichier d'index pour le fichier BAM d'entrée
        SAMTOOLS_INDEX(reads_ch)
    ```

La syntaxe `file()` indique explicitement à Nextflow de gérer ces entrées comme des chemins de fichiers.
Vous pouvez en apprendre plus à ce sujet dans la Quête secondaire [Working with files](../../side_quests/working_with_files.md).

### 2.2. Écrire le processus d'appel de variants et l'appeler dans le workflow

Nous devons écrire la définition du processus dans le fichier de module, l'importer dans le workflow en utilisant une instruction include, et l'appeler sur les lectures d'entrée plus la sortie de l'étape d'indexation et les fichiers accessoires.

#### 2.2.1. Compléter le module pour le processus d'appel de variants

Ouvrez `modules/gatk_haplotypecaller.nf` et examinez le squelette de la définition du processus.

Complétez la définition du processus par vous-même en utilisant les informations fournies ci-dessus, puis vérifiez votre travail par rapport à la solution dans l'onglet "Après" ci-dessous.

=== "Avant"

    ```groovy title="modules/gatk_haplotypecaller.nf" linenums="1"
    #!/usr/bin/env nextflow

    /*
     * Appeler les variants avec GATK HaplotypeCaller
     */
    process GATK_HAPLOTYPECALLER {

        container

        input:

        output:

        script:
        """

        """
    }
    ```

=== "Après"

    ```groovy title="modules/gatk_haplotypecaller.nf" linenums="1" hl_lines="8 11-16 19-20 24-28"
    #!/usr/bin/env nextflow

    /*
     * Appeler les variants avec GATK HaplotypeCaller
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

Vous remarquerez que ce processus a plus d'entrées que la commande GATK elle-même n'en nécessite.
GATK sait chercher le fichier d'index BAM et les fichiers accessoires du génome de référence en fonction de conventions de nommage, mais Nextflow est indépendant du domaine et ne connaît pas ces conventions.
Nous devons les lister explicitement pour que Nextflow les place dans le répertoire de travail au moment de l'exécution ; sinon GATK générera une erreur concernant des fichiers manquants.

De même, nous listons explicitement le fichier d'index du VCF de sortie (`"${input_bam}.vcf.idx"`) pour que Nextflow en garde la trace pour les étapes suivantes.
Nous utilisons la syntaxe `emit:` pour assigner un nom à chaque canal de sortie, ce qui deviendra utile lorsque nous connecterons les sorties dans le bloc publish.

Une fois que vous avez terminé cela, le processus est complet.
Pour l'utiliser dans le workflow, vous devrez importer le module et ajouter un appel de processus.

#### 2.2.2. Importer le nouveau module

Mettez à jour `genomics.nf` pour importer le nouveau module :

=== "Après"

    ```groovy title="genomics.nf" linenums="3" hl_lines="3"
    // Instructions INCLUDE de modules
    include { SAMTOOLS_INDEX } from './modules/samtools_index.nf'
    include { GATK_HAPLOTYPECALLER } from './modules/gatk_haplotypecaller.nf'
    ```

=== "Avant"

    ```groovy title="genomics.nf" linenums="3"
    // Instructions INCLUDE de modules
    include { SAMTOOLS_INDEX } from './modules/samtools_index.nf'
    ```

Le processus est maintenant disponible dans la portée du workflow.

#### 2.2.3. Ajouter l'appel de processus

Ajoutez l'appel de processus dans le corps du workflow, sous `main:` :

=== "Après"

    ```groovy title="genomics.nf" linenums="33" hl_lines="4-12"
        // Créer le fichier d'index pour le fichier BAM d'entrée
        SAMTOOLS_INDEX(reads_ch)

        // Appeler les variants à partir du fichier BAM indexé
        GATK_HAPLOTYPECALLER(
            reads_ch,
            SAMTOOLS_INDEX.out,
            ref_file,
            ref_index_file,
            ref_dict_file,
            intervals_file
        )
    ```

=== "Avant"

    ```groovy title="genomics.nf" linenums="33"
        // Créer le fichier d'index pour le fichier BAM d'entrée
        SAMTOOLS_INDEX(reads_ch)
    ```

Vous devriez reconnaître la syntaxe `*.out` de la série de formation Hello Nextflow ; nous disons à Nextflow de prendre le canal émis par `SAMTOOLS_INDEX` et de le brancher dans l'appel du processus `GATK_HAPLOTYPECALLER`.

!!! note "Note"

    Notez que les entrées sont fournies dans exactement le même ordre dans l'appel au processus que dans le bloc input du processus.
    Dans Nextflow, les entrées sont positionnelles, ce qui signifie que vous _devez_ suivre le même ordre ; et bien sûr il doit y avoir le même nombre d'éléments.

### 2.3. Configurer la gestion de la sortie

Nous devons ajouter les nouvelles sorties à la déclaration publish et configurer où elles vont.

#### 2.3.1. Ajouter des cibles de publication pour les sorties d'appel de variants

Ajoutez les sorties VCF et index à la section `publish:` :

=== "Après"

    ```groovy title="genomics.nf" linenums="45" hl_lines="3-4"
        publish:
        bam_index = SAMTOOLS_INDEX.out
        vcf = GATK_HAPLOTYPECALLER.out.vcf
        vcf_idx = GATK_HAPLOTYPECALLER.out.idx
    }
    ```

=== "Avant"

    ```groovy title="genomics.nf" linenums="45"
        publish:
        bam_index = SAMTOOLS_INDEX.out
    }
    ```

Maintenant nous devons dire à Nextflow où mettre les nouvelles sorties.

#### 2.3.2. Configurer les nouvelles cibles de sortie

Ajoutez des entrées pour les cibles `vcf` et `vcf_idx` dans le bloc `output {}`, en publiant les deux dans un sous-répertoire `vcf/` :

=== "Après"

    ```groovy title="genomics.nf" linenums="51" hl_lines="5-10"
    output {
        bam_index {
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

=== "Avant"

    ```groovy title="genomics.nf" linenums="49"
    output {
        bam_index {
            path 'bam'
        }
    }
    ```

Le VCF et son index sont publiés comme des cibles séparées qui vont toutes deux dans le sous-répertoire `vcf/`.

### 2.4. Exécuter le workflow

Exécutez le workflow étendu, en ajoutant `-resume` cette fois pour ne pas avoir à réexécuter l'étape d'indexation.

```bash
nextflow run genomics.nf -profile test -resume
```

??? success "Sortie de la commande"

    ```console
    N E X T F L O W   ~  version 25.10.2

    ┃ Launching `genomics.nf` [grave_volta] DSL2 - revision: 4790abc96a

    executor >  local (1)
    [2a/e69536] SAMTOOLS_INDEX (1)       | 1 of 1, cached: 1 ✔
    [53/e18e98] GATK_HAPLOTYPECALLER (1) | 1 of 1 ✔
    ```

Maintenant si nous regardons la sortie console, nous voyons les deux processus listés.

Le premier processus a été ignoré grâce au cache, comme prévu, alors que le second processus a été exécuté puisqu'il est tout nouveau.

Vous trouverez les fichiers de sortie dans le répertoire de résultats (sous forme de liens symboliques vers le répertoire de travail).

??? abstract "Contenu du répertoire"

    ```console
    results/
    ├── bam/
    │   └── reads_mother.bam.bai -> ...
    └── vcf/
        ├── reads_mother.bam.vcf -> ...
        └── reads_mother.bam.vcf.idx -> ...
    ```

Si vous ouvrez le fichier VCF, vous devriez voir le même contenu que dans le fichier que vous avez généré en exécutant la commande GATK directement dans le conteneur.

??? abstract "Contenu du fichier"

    ```console title="reads_mother.bam.vcf" linenums="26"
    #CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	reads_mother
    20_10037292_10066351	3480	.	C	CT	503.03	.	AC=2;AF=1.00;AN=2;DP=23;ExcessHet=0.0000;FS=0.000;MLEAC=2;MLEAF=1.00;MQ=60.00;QD=27.95;SOR=1.179	GT:AD:DP:GQ:PL	1/1:0,18:18:54:517,54,0
    20_10037292_10066351	3520	.	AT	A	609.03	.	AC=2;AF=1.00;AN=2;DP=18;ExcessHet=0.0000;FS=0.000;MLEAC=2;MLEAF=1.00;MQ=60.00;QD=33.83;SOR=0.693	GT:AD:DP:GQ:PL	1/1:0,18:18:54:623,54,0
    20_10037292_10066351	3529	.	T	A	155.64	.	AC=1;AF=0.500;AN=2;BaseQRankSum=-0.544;DP=21;ExcessHet=0.0000;FS=1.871;MLEAC=1;MLEAF=0.500;MQ=60.00;MQRankSum=0.000;QD=7.78;ReadPosRankSum=-1.158;SOR=1.034	GT:AD:DP:GQ:PL	0/1:12,8:20:99:163,0,328
    ```

C'est la sortie qui nous intéresse de générer pour chaque échantillon dans notre étude.

### À retenir

Vous savez comment créer un workflow modulaire en deux étapes qui effectue un véritable travail d'analyse et est capable de gérer les particularités des formats de fichiers génomiques comme les fichiers accessoires.

### Et ensuite ?

Faire en sorte que le workflow gère plusieurs échantillons en masse.

---

## 3. Adapter le workflow pour qu'il s'exécute sur un lot d'échantillons

C'est bien d'avoir un workflow qui peut automatiser le traitement sur un seul échantillon, mais que faire si vous avez 1000 échantillons ?
Devez-vous écrire un script bash qui boucle sur tous vos échantillons ?

Non, heureusement ! Faites simplement une petite modification du code et Nextflow gérera cela pour vous aussi.

### 3.1. Mettre à jour l'entrée pour lister trois échantillons

Pour exécuter sur plusieurs échantillons, mettez à jour le profil de test pour fournir un tableau de chemins de fichiers au lieu d'un seul.
C'est une façon rapide de tester l'exécution multi-échantillons ; dans la prochaine étape nous passerons à une approche plus évolutive utilisant un fichier d'entrées.

D'abord, commentez l'annotation de type dans la déclaration du paramètre, puisque les tableaux ne peuvent pas utiliser de déclarations typées :

=== "Après"

    ```groovy title="genomics.nf" linenums="10" hl_lines="1-2"
        // Entrée principale (tableau de trois échantillons)
        reads_bam //: Path
    ```

=== "Avant"

    ```groovy title="genomics.nf" linenums="10"
        // Entrée principale
        reads_bam: Path
    ```

Ensuite, mettez à jour le profil de test pour lister les trois échantillons :

=== "Après"

    ```groovy title="nextflow.config" linenums="4" hl_lines="2-6"
    test {
        params.reads_bam = [
            "${projectDir}/data/bam/reads_mother.bam",
            "${projectDir}/data/bam/reads_father.bam",
            "${projectDir}/data/bam/reads_son.bam"
        ]
        params.reference = "${projectDir}/data/ref/ref.fasta"
        params.reference_index = "${projectDir}/data/ref/ref.fasta.fai"
        params.reference_dict = "${projectDir}/data/ref/ref.dict"
        params.intervals = "${projectDir}/data/ref/intervals.bed"
    }
    ```

=== "Avant"

    ```groovy title="nextflow.config" linenums="4"
    test {
        params.reads_bam = "${projectDir}/data/bam/reads_mother.bam"
        params.reference = "${projectDir}/data/ref/ref.fasta"
        params.reference_index = "${projectDir}/data/ref/ref.fasta.fai"
        params.reference_dict = "${projectDir}/data/ref/ref.dict"
        params.intervals = "${projectDir}/data/ref/intervals.bed"
    }
    ```

La fabrique de canaux dans le corps du workflow (`.fromPath`) accepte plusieurs chemins de fichiers tout aussi bien qu'un seul, donc aucun autre changement n'est nécessaire.

### 3.2. Exécuter le workflow

Essayez d'exécuter le workflow maintenant que la plomberie est configurée pour s'exécuter sur les trois échantillons de test.

```bash
nextflow run genomics.nf -profile test -resume
```

Chose amusante : cela _pourrait fonctionner_, OU cela _pourrait échouer_. Par exemple, voici une exécution qui a réussi :

??? success "Sortie de la commande"

    ```console
    N E X T F L O W   ~  version 25.10.2

    ┃ Launching `genomics.nf` [peaceful_yalow] DSL2 - revision: a256d113ad

    executor >  local (6)
    [4f/7071b0] SAMTOOLS_INDEX (3)       | 3 of 3, cached: 1 ✔
    [7a/89bc43] GATK_HAPLOTYPECALLER (2) | 3 of 3, cached: 1 ✔
    ```

Si votre exécution de workflow a réussi, exécutez-la à nouveau jusqu'à obtenir une erreur comme celle-ci :

??? failure "Sortie de la commande"

    ```console
    N E X T F L O W   ~  version 25.10.2

    ┃ Launching `genomics.nf` [loving_pasteur] DSL2 - revision: d2a8e63076

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

Eh bien, c'est bizarre, étant donné que nous avons explicitement indexé les fichiers BAM dans la première étape du workflow. Pourrait-il y avoir un problème avec la plomberie ?

### 3.3. Dépanner le problème

Nous allons inspecter les répertoires de travail et utiliser l'opérateur `view()` pour comprendre ce qui s'est mal passé.

#### 3.3.1. Vérifier les répertoires de travail pour les appels pertinents

Jetez un œil à l'intérieur du répertoire de travail pour l'appel de processus `GATK_HAPLOTYPECALLER` échoué listé dans la sortie console.

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

Quoi ? Nextflow a placé un fichier d'index dans le répertoire de travail de cet appel de processus, mais c'est le mauvais. Comment cela a-t-il pu se produire ?

#### 3.3.2. Utiliser l'[opérateur view()](https://www.nextflow.io/docs/latest/reference/operator.html#view) pour inspecter le contenu des canaux

Ajoutez ces deux lignes dans le corps du workflow avant l'appel du processus `GATK_HAPLOTYPECALLER` pour visualiser le contenu du canal :

=== "Après"

    ```groovy title="genomics.nf" hl_lines="3-5"
        SAMTOOLS_INDEX(reads_ch)

        // diagnostics temporaires
        reads_ch.view()
        SAMTOOLS_INDEX.out.view()

        // Appeler les variants à partir du fichier BAM indexé
        GATK_HAPLOTYPECALLER(
    ```

=== "Avant"

    ```groovy title="genomics.nf"
        SAMTOOLS_INDEX(reads_ch)

        // Appeler les variants à partir du fichier BAM indexé
        GATK_HAPLOTYPECALLER(
    ```

Ensuite, exécutez à nouveau la commande du workflow.

```bash
nextflow run genomics.nf -profile test
```

Une fois de plus, cela peut réussir ou échouer. Voici à quoi ressemble la sortie des deux appels `.view()` pour une exécution échouée :

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

!!! note "Note"

    Lorsque vous appelez un processus Nextflow sur un canal contenant plusieurs éléments, Nextflow essaiera de paralléliser l'exécution autant que possible, et collectera les sorties dans l'ordre où elles deviennent disponibles.
    La conséquence est que les sorties correspondantes peuvent être collectées dans un ordre différent de celui dans lequel les entrées originales ont été fournies.

Tel qu'écrit actuellement, notre script de workflow suppose que les fichiers d'index sortiront de l'étape d'indexation listés dans le même ordre mère/père/fils que les entrées ont été données.
Mais cela n'est pas garanti, c'est pourquoi parfois (mais pas toujours) les mauvais fichiers sont appariés dans la seconde étape.

Pour corriger cela, nous devons nous assurer que les fichiers BAM et leurs fichiers d'index voyagent ensemble à travers les canaux.

!!! tip "Astuce"

    Les instructions `view()` dans le code du workflow ne font rien, donc ce n'est pas un problème de les laisser.
    Cependant elles encombrent votre sortie console, donc nous recommandons de les supprimer lorsque vous avez fini de dépanner le problème.

### 3.4. Mettre à jour le workflow pour gérer correctement les fichiers d'index

La correction consiste à regrouper chaque fichier BAM avec son index dans un tuple, puis à mettre à jour le processus en aval et la plomberie du workflow pour correspondre.

#### 3.4.1. Changer la sortie du module SAMTOOLS_INDEX en tuple

La façon la plus simple de s'assurer qu'un fichier BAM et son index restent étroitement associés est de les empaqueter ensemble dans un tuple sortant de la tâche d'index.

!!! note "Note"

    Un **tuple** est une liste finie et ordonnée d'éléments qui est couramment utilisée pour retourner plusieurs valeurs d'une fonction. Les tuples sont particulièrement utiles pour passer plusieurs entrées ou sorties entre processus tout en préservant leur association et leur ordre.

Mettez à jour la sortie dans `modules/samtools_index.nf` pour inclure le fichier BAM :

=== "Après"

    ```groovy title="modules/samtools_index.nf" linenums="14" hl_lines="2"
        output:
        tuple path(input_bam), path("${input_bam}.bai")
    ```

=== "Avant"

    ```groovy title="modules/samtools_index.nf" linenums="14"
        output:
        path "${input_bam}.bai"
    ```

De cette façon, chaque fichier d'index sera étroitement couplé avec son fichier BAM d'origine, et la sortie globale de l'étape d'indexation sera un seul canal contenant des paires de fichiers.

#### 3.4.2. Changer l'entrée du module GATK_HAPLOTYPECALLER pour accepter un tuple

Puisque nous avons changé la 'forme' de la sortie du premier processus, nous devons mettre à jour la définition d'entrée du second processus pour correspondre.

Mettez à jour `modules/gatk_haplotypecaller.nf` :

=== "Après"

    ```groovy title="modules/gatk_haplotypecaller.nf" linenums="11" hl_lines="2"
        input:
        tuple path(input_bam), path(input_bam_index)
    ```

=== "Avant"

    ```groovy title="modules/gatk_haplotypecaller.nf" linenums="11"
        input:
        path input_bam
        path input_bam_index
    ```

Maintenant nous devons mettre à jour le workflow pour refléter la nouvelle structure de tuple dans l'appel de processus et les cibles de publication.

#### 3.4.3. Mettre à jour l'appel à GATK_HAPLOTYPECALLER dans le workflow

Nous n'avons plus besoin de fournir le `reads_ch` original au processus `GATK_HAPLOTYPECALLER`, puisque le fichier BAM est maintenant regroupé dans le canal de sortie par `SAMTOOLS_INDEX`.

Mettez à jour l'appel dans `genomics.nf` :

=== "Après"

    ```groovy title="genomics.nf" linenums="42" hl_lines="2"
        GATK_HAPLOTYPECALLER(
            SAMTOOLS_INDEX.out,
    ```

=== "Avant"

    ```groovy title="genomics.nf" linenums="42"
        GATK_HAPLOTYPECALLER(
            reads_ch,
            SAMTOOLS_INDEX.out,
    ```

Enfin, nous devons mettre à jour les cibles de publication pour refléter la nouvelle structure de sortie.

#### 3.4.4. Mettre à jour la cible de publication pour la sortie BAM indexé

Puisque la sortie de SAMTOOLS_INDEX est maintenant un tuple contenant à la fois le fichier BAM et son index, renommez la cible de publication de `bam_index` à `indexed_bam` pour mieux refléter son contenu :

=== "Après"

    ```groovy title="genomics.nf" linenums="46" hl_lines="2 8"
        publish:
        indexed_bam = SAMTOOLS_INDEX.out
        vcf = GATK_HAPLOTYPECALLER.out.vcf
        vcf_idx = GATK_HAPLOTYPECALLER.out.idx
    }

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

=== "Avant"

    ```groovy title="genomics.nf" linenums="46"
        publish:
        bam_index = SAMTOOLS_INDEX.out
        vcf = GATK_HAPLOTYPECALLER.out.vcf
        vcf_idx = GATK_HAPLOTYPECALLER.out.idx
    }

    output {
        bam_index {
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

Avec ces changements, le BAM et son index sont garantis de voyager ensemble, donc l'appariement sera toujours correct.

### 3.5. Exécuter le workflow corrigé

Exécutez à nouveau le workflow pour vous assurer que cela fonctionnera de manière fiable à l'avenir.

```bash
nextflow run genomics.nf -profile test
```

Cette fois (et à chaque fois) tout devrait s'exécuter correctement :

??? success "Sortie de la commande"

    ```console
    N E X T F L O W   ~  version 25.10.2

    ┃ Launching `genomics.nf` [special_goldstine] DSL2 - revision: 4cbbf6ea3e

    executor >  local (6)
    [d6/10c2c4] SAMTOOLS_INDEX (1)       | 3 of 3 ✔
    [88/1783aa] GATK_HAPLOTYPECALLER (2) | 3 of 3 ✔
    ```

Le répertoire de résultats contient maintenant à la fois les fichiers BAM et BAI pour chaque échantillon (du tuple), ainsi que les sorties VCF :

??? abstract "Contenu du répertoire de résultats"

    ```console
    results/
    ├── bam/
    │   ├── reads_father.bam -> ...
    │   ├── reads_father.bam.bai -> ...
    │   ├── reads_mother.bam -> ...
    │   ├── reads_mother.bam.bai -> ...
    │   ├── reads_son.bam -> ...
    │   └── reads_son.bam.bai -> ...
    └── vcf/
        ├── reads_father.bam.vcf -> ...
        ├── reads_father.bam.vcf.idx -> ...
        ├── reads_mother.bam.vcf -> ...
        ├── reads_mother.bam.vcf.idx -> ...
        ├── reads_son.bam.vcf -> ...
        └── reads_son.bam.vcf.idx -> ...
    ```

En regroupant les fichiers associés dans des tuples, nous avons assuré que les bons fichiers voyagent toujours ensemble à travers le workflow.
Le workflow traite maintenant n'importe quel nombre d'échantillons de manière fiable, mais les lister individuellement dans la configuration n'est pas très évolutif.
Dans la prochaine étape, nous passerons à la lecture des entrées à partir d'un fichier.

### À retenir

Vous savez comment faire en sorte que votre workflow s'exécute sur plusieurs échantillons (indépendamment).

### Et ensuite ?

Faciliter la gestion des échantillons en masse.

---

## 4. Faire en sorte que le workflow accepte un fichier texte contenant un lot de fichiers d'entrée

Une façon très courante de fournir plusieurs fichiers de données d'entrée à un workflow est de le faire avec un fichier texte contenant les chemins de fichiers.
Il peut être aussi simple qu'un fichier texte listant un chemin de fichier par ligne et rien d'autre, ou le fichier peut contenir des métadonnées supplémentaires, auquel cas il est souvent appelé une feuille d'échantillons (samplesheet).

Ici nous allons vous montrer comment faire le cas simple.

### 4.1. Examiner le fichier texte fourni listant les chemins de fichiers d'entrée

Nous avons déjà créé un fichier texte listant les chemins de fichiers d'entrée, appelé `sample_bams.txt`, que vous pouvez trouver dans le répertoire `data/`.

```txt title="sample_bams.txt"
/workspaces/training/nf4-science/genomics/data/bam/reads_mother.bam
/workspaces/training/nf4-science/genomics/data/bam/reads_father.bam
/workspaces/training/nf4-science/genomics/data/bam/reads_son.bam
```

Comme vous pouvez le voir, nous avons listé un chemin de fichier par ligne, et ce sont des chemins absolus.

!!! note "Note"

    Les fichiers que nous utilisons ici sont simplement sur le système de fichiers local de votre GitHub Codespaces, mais nous pourrions également pointer vers des fichiers dans le stockage cloud.
    Si vous n'utilisez pas l'environnement Codespaces fourni, vous devrez peut-être adapter les chemins de fichiers pour correspondre à votre configuration locale.

### 4.2. Mettre à jour le paramètre et le profil de test

Changez le paramètre `reads_bam` pour pointer vers le fichier `sample_bams.txt` au lieu de lister les échantillons individuels.

Restaurez l'annotation de type dans le bloc params (puisque c'est à nouveau un seul chemin) :

=== "Après"

    ```groovy title="genomics.nf" linenums="10" hl_lines="1-2"
        // Entrée principale (fichier de fichiers d'entrée, un par ligne)
        reads_bam: Path
    ```

=== "Avant"

    ```groovy title="genomics.nf" linenums="10"
        // Entrée principale (tableau de trois échantillons)
        reads_bam
    ```

Ensuite, mettez à jour le profil de test pour pointer vers le fichier texte :

=== "Après"

    ```groovy title="nextflow.config" linenums="4" hl_lines="2"
    test {
        params.reads_bam = "${projectDir}/data/sample_bams.txt"
        params.reference = "${projectDir}/data/ref/ref.fasta"
        params.reference_index = "${projectDir}/data/ref/ref.fasta.fai"
        params.reference_dict = "${projectDir}/data/ref/ref.dict"
        params.intervals = "${projectDir}/data/ref/intervals.bed"
    }
    ```

=== "Avant"

    ```groovy title="nextflow.config" linenums="4"
    test {
        params.reads_bam = [
            "${projectDir}/data/bam/reads_mother.bam",
            "${projectDir}/data/bam/reads_father.bam",
            "${projectDir}/data/bam/reads_son.bam"
        ]
        params.reference = "${projectDir}/data/ref/ref.fasta"
        params.reference_index = "${projectDir}/data/ref/ref.fasta.fai"
        params.reference_dict = "${projectDir}/data/ref/ref.dict"
        params.intervals = "${projectDir}/data/ref/intervals.bed"
    }
    ```

La liste des fichiers ne vit plus du tout dans le code, ce qui est un grand pas dans la bonne direction.

### 4.3. Mettre à jour la fabrique de canaux pour lire les lignes d'un fichier

Actuellement, notre fabrique de canaux d'entrée traite tous les fichiers que nous lui donnons comme les entrées de données que nous voulons alimenter au processus d'indexation.
Puisque nous lui donnons maintenant un fichier qui liste les chemins de fichiers d'entrée, nous devons changer son comportement pour analyser le fichier et traiter les chemins de fichiers qu'il contient comme les entrées de données.

Nous pouvons faire cela en utilisant le même modèle que nous avons utilisé dans la [Partie 2 de Hello Nextflow](../../hello_nextflow/02_hello_channels.md#42-use-the-splitcsv-operator-to-parse-the-file) : appliquer l'opérateur [`splitCsv()`](https://nextflow.io/docs/latest/reference/operator.html#splitcsv) pour analyser le fichier, puis une opération `map` pour sélectionner le premier champ de chaque ligne.

=== "Après"

    ```groovy title="genomics.nf" linenums="24" hl_lines="1-4"
        // Créer le canal d'entrée à partir d'un fichier CSV listant les chemins de fichiers d'entrée
        reads_ch = Channel.fromPath(params.reads_bam)
                .splitCsv()
                .map { line -> file(line[0]) }
    ```

=== "Avant"

    ```groovy title="genomics.nf" linenums="24"
        // Créer le canal d'entrée (fichier unique via paramètre CLI)
        reads_ch = channel.fromPath(params.reads_bam)
    ```

Techniquement nous pourrions faire cela plus simplement en utilisant l'opérateur [`.splitText()`](https://www.nextflow.io/docs/latest/reference/operator.html#operator-splittext), puisque notre fichier d'entrée ne contient actuellement que des chemins de fichiers.
Cependant, en utilisant l'opérateur plus polyvalent `splitCsv` (complété par `map`), nous pouvons rendre notre workflow pérenne au cas où nous déciderions d'ajouter des métadonnées au fichier contenant les chemins de fichiers.

!!! tip "Astuce"

    Si vous n'êtes pas sûr·e de comprendre ce que font les opérateurs ici, c'est une autre excellente occasion d'utiliser l'opérateur `.view()` pour voir à quoi ressemble le contenu du canal avant et après leur application.

### 4.4. Exécuter le workflow

Exécutez le workflow une dernière fois. Cela devrait produire le même résultat qu'auparavant, n'est-ce pas ?

```bash
nextflow run genomics.nf -profile test -resume
```

??? success "Sortie de la commande"

    ```console
    N E X T F L O W   ~  version 25.10.2

    ┃ Launching `genomics.nf` [sick_albattani] DSL2 - revision: 46d84642f6

    [18/23b4bb] SAMTOOLS_INDEX (1)       | 3 of 3, cached: 3 ✔
    [12/f727bb] GATK_HAPLOTYPECALLER (3) | 3 of 3, cached: 3 ✔
    ```

Oui ! En fait, Nextflow détecte correctement que les appels de processus sont exactement les mêmes, et ne se donne même pas la peine de tout réexécuter, puisque nous exécutions avec `-resume`.

Et c'est tout ! Notre simple workflow d'appel de variants a toutes les fonctionnalités de base que nous voulions.

### À retenir

Vous savez comment créer un workflow modulaire en plusieurs étapes pour indexer un fichier BAM et appliquer l'appel de variants par échantillon en utilisant GATK.

Plus généralement, vous avez appris à utiliser les composants et la logique essentiels de Nextflow pour construire un pipeline génomique simple qui effectue un véritable travail, en tenant compte des particularités des formats de fichiers génomiques et des exigences des outils.

### Et ensuite ?

Célébrez votre succès et prenez une pause extra longue !

Dans la prochaine partie de ce cours, vous apprendrez comment transformer ce simple workflow d'appel de variants par échantillon pour appliquer l'appel de variants conjoint aux données.
