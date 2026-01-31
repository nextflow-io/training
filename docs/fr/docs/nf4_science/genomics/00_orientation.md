# Orientation

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Traduction assistée par IA - [en savoir plus et suggérer des améliorations](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

L'environnement de formation contient tous les logiciels, le code et les données nécessaires pour suivre cette formation, vous n'avez donc rien à installer vous-même.
Cependant, vous avez besoin d'un compte (gratuit) pour vous connecter, et vous devriez prendre quelques minutes pour vous familiariser avec l'interface.

Si vous ne l'avez pas encore fait, veuillez suivre [ce lien](../../../envsetup/) avant de continuer.

## Matériel fourni

Tout au long de cette formation, nous travaillerons dans le répertoire `nf4-science/genomics/`, dans lequel vous devez vous déplacer lorsque vous ouvrez l'espace de travail de formation.
Ce répertoire contient tous les fichiers de code, les données de test et les fichiers accessoires dont vous aurez besoin.

N'hésitez pas à explorer le contenu de ce répertoire ; le moyen le plus simple est d'utiliser l'explorateur de fichiers sur le côté gauche de l'espace de travail de formation dans l'interface VSCode.
Vous pouvez également utiliser la commande `tree`.
Tout au long de la formation, nous utilisons la sortie de `tree` pour représenter la structure et le contenu des répertoires sous une forme lisible, parfois avec des modifications mineures pour plus de clarté.

Ici, nous générons une table des matières jusqu'au deuxième niveau :

```bash
tree . -L 2
```

Si vous exécutez cette commande à l'intérieur de `nf4-science/genomics`, vous devriez voir la sortie suivante :

```console title="Contenu du répertoire"

.
├── data
│   ├── bam
│   ├── ref
│   ├── sample_bams.txt
│   └── samplesheet.csv
├── genomics-1.nf
├── genomics-2.nf
├── genomics-3.nf
├── genomics-4.nf
├── nextflow.config
└── solutions
    ├── modules
    ├── nf-test.config
    └── tests

6 directories, 8 files

```

!!!note

    Ne vous inquiétez pas si cela semble beaucoup ; nous passerons en revue les éléments pertinents à chaque étape de la formation.
    Ceci est juste destiné à vous donner un aperçu.

**Voici un résumé de ce que vous devez savoir pour commencer :**

- **Les fichiers `.nf`** sont des scripts de workflow nommés en fonction de la partie de la formation dans laquelle ils sont utilisés.

- **Le fichier `nextflow.config`** est un fichier de configuration qui définit des propriétés minimales de l'environnement.
  Vous pouvez l'ignorer pour l'instant.

- **Le répertoire `data`** contient les données d'entrée et les ressources associées, décrites plus tard dans la formation.

- **Le répertoire `solutions`** contient les fichiers de modules et les configurations de tests qui résultent des parties 3 et 4 de la formation.
  Ils sont destinés à être utilisés comme référence pour vérifier votre travail et résoudre tout problème.

!!!tip

    Si pour une raison quelconque vous sortez de ce répertoire, vous pouvez toujours exécuter cette commande pour y revenir :

    ```bash
    cd /workspaces/training/nf4-science/genomics
    ```

Maintenant, pour commencer la formation, cliquez sur la flèche dans le coin inférieur droit de cette page.
