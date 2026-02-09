# Orientation

L'environnement GitHub Codespaces contient tous les logiciels, le code et les données nécessaires pour suivre cette formation, vous n'avez donc rien à installer vous-même.
Cependant, vous avez besoin d'un compte (gratuit) pour vous connecter, et vous devriez prendre quelques minutes pour vous familiariser avec l'interface.

Si vous ne l'avez pas encore fait, veuillez suivre [ce lien](../../envsetup/) avant d'aller plus loin.

## Matériel fourni

Tout au long de cette formation, nous travaillerons dans le répertoire `side-quests/`.
Ce répertoire contient tous les fichiers de code, les données de test et les fichiers accessoires dont vous aurez besoin.

N'hésitez pas à explorer le contenu de ce répertoire ; le moyen le plus simple de le faire est d'utiliser l'explorateur de fichiers sur le côté gauche de l'espace de travail GitHub Codespaces.
Alternativement, vous pouvez utiliser la commande `tree`.
Tout au long du cours, nous utilisons la sortie de `tree` pour représenter la structure et le contenu des répertoires sous une forme lisible, parfois avec des modifications mineures pour plus de clarté.

Ici, nous générons une table des matières jusqu'au deuxième niveau :

```bash
tree . -L 2
```

Si vous exécutez cette commande dans `side-quests`, vous devriez voir la sortie suivante :

```console title="Directory contents"
.
├── metadata
├── nf-core
├── nf-test
├── solutions
├── splitting_and_grouping
└── workflows_of_workflows
```

**Voici un résumé de ce que vous devez savoir pour commencer :**

- **Chaque répertoire correspond à une quête secondaire individuelle.**
  Leur contenu est détaillé sur la page de la quête secondaire correspondante.

- **Le répertoire `solutions`** contient les scripts de workflow et/ou de module complétés qui résultent de l'exécution des différentes étapes de chaque quête secondaire.
  Ils sont destinés à être utilisés comme référence pour vérifier votre travail et résoudre d'éventuels problèmes.

!!!tip

    Si pour une raison quelconque vous sortez de ce répertoire, vous pouvez toujours exécuter cette commande pour y revenir :

    ```bash
    cd /workspaces/training/side-quests
    ```

Maintenant, pour commencer le cours, cliquez sur la flèche dans le coin inférieur droit de cette page.
