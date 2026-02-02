# Partie 3 : Déplacer le code dans des modules

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Traduction assistée par IA - [en savoir plus et suggérer des améliorations](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

Dans la première partie de cette formation, vous avez construit un pipeline d'appel de variants complètement linéaire qui traitait les données de chaque échantillon indépendamment des autres.

Dans la deuxième partie, nous vous avons montré comment utiliser les canaux et les opérateurs de canaux pour implémenter l'appel de variants joint avec GATK, en vous appuyant sur le pipeline de la Partie 1.

Dans cette partie, nous allons vous montrer comment convertir le code de ce workflow en modules. Pour suivre cette partie de la formation, vous devriez avoir terminé la Partie 1 et la Partie 2, ainsi que [Hello Modules](../../../hello_nextflow/hello_modules.md), qui couvre les bases des modules.

---

## 0. Échauffement

Lorsque nous avons commencé à développer notre workflow, nous avons tout mis dans un seul fichier de code.
Il est maintenant temps de s'attaquer à la **modularisation** de notre code, _c'est-à-dire_ extraire les définitions de processus dans des modules.

Nous allons commencer avec le même workflow que dans la Partie 2, que nous vous avons fourni dans le fichier `genomics-3.nf`.

!!! note "Note"

     Assurez-vous d'être dans le bon répertoire de travail :
     `cd /workspaces/training/nf4-science/genomics`

Exécutez le workflow pour vérifier le point de départ :

```bash
nextflow run genomics-3.nf -resume
```

```console title="Sortie"
 N E X T F L O W   ~  version 25.10.2

Launching `genomics-3.nf` [serene_borg] DSL2 - revision: 0cbebb67a1

executor >  local (7)
[6f/83ee72] SAMTOOLS_INDEX (3)       | 3 of 3 ✔
[53/b9d342] GATK_HAPLOTYPECALLER (1) | 3 of 3 ✔
[0c/fa6d15] GATK_JOINTGENOTYPING     | 1 of 1 ✔
```

Il y aura maintenant un répertoire `work` et un répertoire `results_genomics` à l'intérieur de votre répertoire de projet.

### À retenir

Vous êtes prêt·e à commencer la modularisation de votre workflow.

### Et ensuite ?

Déplacer les processus du workflow Genomics dans des modules.

---

## 1. Déplacer les processus dans des modules

Comme vous l'avez appris dans [Hello Modules](../../../hello_nextflow/hello_modules.md), vous pouvez créer un module simplement en copiant la définition du processus dans son propre fichier, dans n'importe quel répertoire, et vous pouvez nommer ce fichier comme vous le souhaitez.

Pour des raisons qui deviendront claires plus tard (en particulier lorsque nous aborderons les tests), dans cette formation nous suivrons la convention de nommer le fichier `main.nf`, et de le placer dans une structure de répertoires nommée d'après la boîte à outils et la commande.

### 1.1. Créer un module pour le processus `SAMTOOLS_INDEX`

Dans le cas du processus `SAMTOOLS_INDEX`, 'samtools' est la boîte à outils et 'index' est la commande. Nous allons donc créer une structure de répertoires `modules/samtools/index` et placer la définition du processus `SAMTOOLS_INDEX` dans le fichier `main.nf` à l'intérieur de ce répertoire.

```bash
mkdir -p modules/samtools/index
touch modules/samtools/index/main.nf
```

Ouvrez le fichier `main.nf` et copiez-y la définition du processus `SAMTOOLS_INDEX`.

```groovy title="modules/samtools/index/main.nf" linenums="1"
/*
 * Générer le fichier d'index BAM
 */
process SAMTOOLS_INDEX {

    container 'community.wave.seqera.io/library/samtools:1.20--b5dfbd93de237464'

    input:
    path input_bam

    output:
    tuple path(input_bam), path("${input_bam}.bai")

    script:
    """
    samtools index '$input_bam'
    """
}
```

Ensuite, supprimez la définition du processus `SAMTOOLS_INDEX` de `genomics-3.nf`, et ajoutez une déclaration d'import pour le module avant la définition du processus suivant, comme ceci :

=== "Après"

    ```groovy title="genomics-3.nf" linenums="1" hl_lines="1 2"
    // Inclure les modules
    include { SAMTOOLS_INDEX } from './modules/samtools/index/main.nf'

    /*
     * Appeler les variants avec GATK HaplotypeCaller
     */
    process GATK_HAPLOTYPECALLER {
    ```

=== "Avant"

    ```groovy title="genomics-3.nf" linenums="1" hl_lines="1"
    /*
     * Appeler les variants avec GATK HaplotypeCaller
     */
    process GATK_HAPLOTYPECALLER {
    ```

Vous pouvez maintenant exécuter à nouveau le workflow, et il devrait fonctionner de la même manière qu'auparavant. Si vous fournissez l'option `-resume`, aucune nouvelle tâche ne devrait même avoir besoin d'être exécutée :

```bash
nextflow run genomics-3.nf -resume
```

??? success "Sortie de la commande"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `genomics-3.nf` [sleepy_snyder] DSL2 - revision: aa68d06c43

    [0f/71b55e] SAMTOOLS_INDEX (1)       | 3 of 3, cached: 3 ✔
    [f1/18971b] GATK_HAPLOTYPECALLER (3) | 3 of 3, cached: 3 ✔
    [0c/fa6d15] GATK_JOINTGENOTYPING     | 1 of 1, cached: 1 ✔
    ```

### 1.2. Créer des modules pour les processus `GATK_HAPLOTYPECALLER` et `GATK_JOINTGENOTYPING`

Répétez les mêmes étapes pour les processus restants.
Pour chaque processus :

1. Créez la structure de répertoires (`modules/gatk/haplotypecaller/` et `modules/gatk/jointgenotyping/`)
2. Créez un fichier `main.nf` contenant la définition du processus
3. Supprimez la définition du processus de `genomics-3.nf`
4. Ajoutez une déclaration d'import pour le module

Une fois terminé, vérifiez que la structure de votre répertoire modules est correcte en exécutant :

```bash
tree modules/
```

??? abstract "Contenu du répertoire"

    ```console
    modules/
    ├── gatk
    │   ├── haplotypecaller
    │   │   └── main.nf
    │   └── jointgenotyping
    │       └── main.nf
    └── samtools
        └── index
            └── main.nf

    5 directories, 3 files
    ```

Vous devriez également avoir quelque chose comme ceci dans le fichier workflow principal, après la section des paramètres :

```
include { SAMTOOLS_INDEX } from './modules/samtools/index/main.nf'
include { GATK_HAPLOTYPECALLER } from './modules/gatk/haplotypecaller/main.nf'
include { GATK_JOINTGENOTYPING } from './modules/gatk/jointgenotyping/main.nf'

workflow {
```

### À retenir

Vous vous êtes exercé·e à modulariser un workflow, avec le workflow génomique comme exemple.

### Et ensuite ?

Tester le workflow modularisé.

---

## 2. Tester le workflow modularisé

Exécutez le workflow modularisé pour vérifier que tout fonctionne toujours.

```bash
nextflow run genomics-3.nf -resume
```

```console title="Sortie"
 N E X T F L O W   ~  version 25.10.2

Launching `genomics-3.nf` [astonishing_venter] DSL2 - revision: ca27264c13

[6f/83ee72] SAMTOOLS_INDEX (3)       | 3 of 3, cached: 3 ✔
[53/b9d342] GATK_HAPLOTYPECALLER (3) | 3 of 3, cached: 3 ✔
[0c/fa6d15] GATK_JOINTGENOTYPING     | 1 of 1, cached: 1 ✔
```

Tout fonctionne toujours, y compris la capacité de reprise du pipeline.
Les résultats continuent d'être publiés dans le répertoire `results_genomics`.

```console title="Contenu du répertoire"
results_genomics/
├── family_trio.joint.vcf
├── family_trio.joint.vcf.idx
├── gvcf
│   ├── reads_father.bam.g.vcf
│   ├── reads_father.bam.g.vcf.idx
│   ├── reads_mother.bam.g.vcf
│   ├── reads_mother.bam.g.vcf.idx
│   ├── reads_son.bam.g.vcf
│   └── reads_son.bam.g.vcf.idx
└── indexed_bam
    ├── reads_father.bam
    ├── reads_father.bam.bai
    ├── reads_mother.bam
    ├── reads_mother.bam.bai
    ├── reads_son.bam
    └── reads_son.bam.bai
```

### À retenir

Vous avez modularisé un workflow et vérifié qu'il fonctionne toujours de la même manière qu'auparavant.

### Et ensuite ?

Examiner ce que vous avez appris et anticiper les tests.

---

## 3. Résumé

Vous avez modularisé le workflow, et rien n'a changé dans le fonctionnement du pipeline.
C'est intentionnel : vous avez restructuré le code sans impacter sa fonction.

Les modules contiennent uniquement la logique des processus, les rendant propres et réutilisables.
Le script principal contrôle ce qui est publié et où, tandis que les modules restent concentrés sur leur tâche de calcul.

Vous avez posé les bases de choses qui rendront votre code plus facile à maintenir.
Par exemple, vous pouvez maintenant ajouter des tests à votre pipeline en utilisant le framework nf-test.
C'est ce que nous examinerons dans la prochaine partie de cette formation.
