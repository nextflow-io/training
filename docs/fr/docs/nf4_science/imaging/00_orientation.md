# Orientation

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Traduction assistée par IA - [en savoir plus et suggérer des améliorations](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

Cette orientation suppose que vous avez déjà ouvert l'environnement de formation en cliquant sur le bouton "Open in GitHub Codespaces".
Si ce n'est pas le cas, veuillez le faire maintenant, idéalement dans une seconde fenêtre ou un second onglet de navigateur afin de pouvoir vous référer à ces instructions.

[![Open in GitHub Codespaces](https://github.com/codespaces/badge.svg)](https://github.com/codespaces/new?hide_repo_select=true&ref=master&repo=290790519&skip_quickstart=true&machine=premiumLinux&devcontainer_path=.devcontainer%2Fdevcontainer.json)

!!!warning "Exigence de taille de machine"

    Assurez-vous de sélectionner une **machine à 8 cœurs** lors de la création de votre Codespace pour cette formation. Les workflows de bio-imagerie nécessitent des ressources de calcul supplémentaires.

## GitHub Codespaces

L'environnement GitHub Codespaces contient tous les logiciels, le code et les données nécessaires pour suivre cette formation, vous n'avez donc rien à installer vous-même.
Cependant, vous avez besoin d'un compte GitHub (gratuit) pour vous connecter, et si vous n'êtes pas familier avec l'interface, vous devriez prendre quelques minutes pour vous familiariser en complétant le mini-cours [Orientation GitHub Codespaces](../../envsetup/index.md).

## Pré-téléchargement des images Docker

Une fois que vous avez ouvert votre Codespace, pré-téléchargeons toutes les images Docker dont nous aurons besoin pour cette formation.
Cela fera gagner du temps plus tard et garantira une exécution fluide des workflows.

Ouvrez un nouvel onglet de terminal et exécutez la commande suivante :

```bash
nextflow run nf-core/molkart -profile docker,test -stub -resume --outdir results
```

Cette commande téléchargera toutes les images Docker nécessaires en arrière-plan.
Vous pouvez continuer avec le reste de l'orientation pendant que cela s'exécute.

!!!tip "Astuce"

    Le flag `-stub` permet au pipeline de s'exécuter rapidement sans traiter de vraies données, ce qui est parfait pour télécharger les images. Vous pouvez surveiller la progression dans l'onglet du terminal.

## Répertoire de travail

Tout au long de cette formation, nous travaillerons dans le répertoire `nf4-science/imaging/`.

Changez de répertoire maintenant en exécutant cette commande dans le terminal :

```bash
cd nf4-science/imaging/
```

!!!tip "Astuce"

    Si pour une raison quelconque vous sortez de ce répertoire, vous pouvez toujours utiliser le chemin complet pour y revenir, en supposant que vous travaillez dans l'environnement de formation GitHub Codespaces :

    ```bash
    cd /workspaces/training/nf4-science/imaging
    ```

**Maintenant, pour commencer la formation, cliquez sur la flèche dans le coin inférieur droit de cette page.**
