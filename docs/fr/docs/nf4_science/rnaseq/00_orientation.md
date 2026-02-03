# Orientation

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Traduction assistée par IA - [en savoir plus et suggérer des améliorations](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

L'environnement de formation contient tous les logiciels, le code et les données nécessaires pour suivre ce cours de formation, vous n'avez donc rien à installer vous-même.
Cependant, vous avez besoin d'un compte (gratuit) pour vous connecter, et vous devriez prendre quelques minutes pour vous familiariser avec l'interface.

Si vous ne l'avez pas encore fait, veuillez suivre le mini-cours [Configuration de l'environnement](../../envsetup/) avant d'aller plus loin.

## Matériel fourni

Tout au long de ce cours de formation, nous travaillerons dans le répertoire `nf4-science/rnaseq/`, dans lequel vous devez vous placer lorsque vous ouvrez l'espace de travail de formation.
Ce répertoire contient tous les fichiers de code, les données de test et les fichiers accessoires dont vous aurez besoin.

N'hésitez pas à explorer le contenu de ce répertoire ; la façon la plus simple de le faire est d'utiliser l'explorateur de fichiers sur le côté gauche de l'espace de travail de formation dans l'interface VSCode.
Alternativement, vous pouvez utiliser la commande `tree`.
Tout au long du cours, nous utilisons la sortie de `tree` pour représenter la structure et le contenu des répertoires sous une forme lisible, parfois avec des modifications mineures pour plus de clarté.

Ici, nous générons une table des matières jusqu'au deuxième niveau :

```bash
tree . -L 3
```

??? success "Contenu du répertoire"

    ```console
    rnaseq
    ├── data
    │   ├── genome.fa
    │   ├── paired-end.csv
    │   ├── reads
    │   │   ├── ENCSR000COQ1_1.fastq.gz
    │   │   ├── ENCSR000COQ1_2.fastq.gz
    │   │   ├── ENCSR000COQ2_1.fastq.gz
    │   │   ├── ENCSR000COQ2_2.fastq.gz
    │   │   ├── ENCSR000COR1_1.fastq.gz
    │   │   ├── ENCSR000COR1_2.fastq.gz
    │   │   ├── ENCSR000COR2_1.fastq.gz
    │   │   ├── ENCSR000COR2_2.fastq.gz
    │   │   ├── ENCSR000CPO1_1.fastq.gz
    │   │   ├── ENCSR000CPO1_2.fastq.gz
    │   │   ├── ENCSR000CPO2_1.fastq.gz
    │   │   └── ENCSR000CPO2_2.fastq.gz
    │   └── single-end.csv
    ├── nextflow.config
    ├── rnaseq.nf
    └── solutions
        ├── modules
        │   ├── fastqc.nf
        │   ├── fastqc_pe.nf
        │   ├── hisat2_align.nf
        │   ├── hisat2_align_pe.nf
        │   ├── multiqc.nf
        │   ├── trim_galore.nf
        │   └── trim_galore_pe.nf
        ├── rnaseq-2.1.nf
        ├── rnaseq-2.2.nf
        ├── rnaseq-2.3.nf
        ├── rnaseq-3.1.nf
        ├── rnaseq-3.2.nf
        └── rnaseq_pe-3.3.nf
    ```

!!!note

    Ne vous inquiétez pas si cela semble beaucoup ; nous passerons en revue les éléments pertinents à chaque étape du cours.
    Ceci est simplement destiné à vous donner un aperçu.

**Voici un résumé de ce que vous devez savoir pour commencer :**

- **Le fichier `rnaseq.nf`** est l'ébauche du script de workflow que nous allons développer.

- **Le fichier `nextflow.config`** est un fichier de configuration qui définit les propriétés minimales de l'environnement. Vous pouvez l'ignorer pour le moment.

- **Le répertoire `data`** contient les données d'entrée et les ressources associées :

  - _Un génome de référence_ appelé `genome.fa` constitué d'une petite région du chromosome 20 humain (de hg19/b37).
  - _Des données RNAseq_ qui ont été extraites d'une petite région pour réduire la taille des fichiers, dans le répertoire `reads/`.
  - _Des fichiers CSV_ listant les identifiants et les chemins des fichiers de données d'exemple, pour un traitement par lots.

- **Le répertoire `solutions`** contient les scripts de workflow et les modules complets qui résultent de chaque étape du cours.
  Ils sont destinés à être utilisés comme référence pour vérifier votre travail et résoudre tout problème.
  Le numéro dans le nom du fichier correspond à l'étape de la partie pertinente du cours.

!!!tip

    Si pour une raison quelconque vous sortez de ce répertoire, vous pouvez toujours exécuter cette commande pour y retourner :

    ```bash
    cd /workspaces/training/nf4-science/rnaseq
    ```

Maintenant, pour commencer le cours, cliquez sur la flèche dans le coin inférieur droit de cette page.
