---
title: Versions de Nextflow
description: Comprendre et gérer l'évolution des versions de syntaxe de Nextflow
hide:
  - toc
  - footer
---

## Version de syntaxe Nextflow actuellement prise en charge et exigences

Depuis la version 3.0 du portail de formation, tous nos cours de formation sont basés sur la version 25.10.2 de Nextflow, sauf indication contraire sur la page d'index du cours (à l'exception des supports obsolètes ou archivés qui peuvent ne pas inclure de mention de version).

Étant donné que les cours utilisent désormais des entrées typées au niveau du workflow ainsi que des directives de sortie au niveau du workflow, ils nécessitent l'utilisation de l'analyseur syntaxique V2.
Si vous prévoyez d'utiliser l'environnement que nous fournissons via [Github Codespaces](../envsetup/01_setup.md) ou les [devcontainers locaux](../envsetup/03_devcontainer.md), vous n'avez rien à faire sauf indication contraire dans les instructions du cours.
Cependant, si vous prévoyez de suivre les formations dans votre propre environnement ([Installation manuelle](../envsetup/02_local.md)), vous devrez vous assurer d'utiliser Nextflow version 25.10.2 ou ultérieure avec l'analyseur syntaxique v2 activé.

## Versions antérieures des supports de formation

Nos supports de formation sont versionnés depuis février 2025.

Vous pouvez accéder aux versions antérieures des supports de formation qui fonctionnent avec les versions de Nextflow **antérieures à 25.10.2** via le menu déroulant en haut de chaque page qui affiche la version numérotée des supports de formation.
Lorsque vous sélectionnez une version antérieure des supports de formation, les liens vers l'environnement de formation spécifieront automatiquement la version correspondante de l'environnement.

## Autres informations sur les versions de syntaxe Nextflow

Nextflow a deux concepts de versionnement distincts qui sont parfois confondus : les **versions DSL** et les **versions de l'analyseur syntaxique**.

**DSL1 vs DSL2** fait référence à des façons fondamentalement différentes d'écrire des pipelines Nextflow.
DSL1 était la syntaxe originale où les processus étaient implicitement connectés via des canaux.
DSL2, introduit dans Nextflow 20.07, a ajouté des fonctionnalités de modularité : la possibilité d'importer des processus et des workflows depuis d'autres fichiers, des blocs `workflow` explicites et des sorties de processus nommées.
DSL1 a été déprécié dans Nextflow 22.03 et supprimé dans la version 22.12.
Tout le code Nextflow moderne utilise DSL2.

**Analyseur syntaxique v1 vs v2** fait référence à différents analyseurs qui fonctionnent tous deux avec le code DSL2.
L'analyseur v1 est l'analyseur original, plus permissif.
L'analyseur v2 est plus strict et active de nouvelles fonctionnalités du langage telles que le typage statique (entrées et sorties typées) et les directives de sortie au niveau du workflow.
L'analyseur v2 fournit également de meilleurs messages d'erreur et détecte davantage d'erreurs au moment de l'analyse plutôt qu'à l'exécution.
L'analyseur v2 deviendra la valeur par défaut dans Nextflow 26.04.

En résumé : DSL2 est le langage que vous écrivez ; la version de l'analyseur syntaxique détermine avec quelle rigueur ce langage est interprété et quelles fonctionnalités avancées sont disponibles.

### Vérifier et définir la version de Nextflow

Vous pouvez vérifier quelle version de Nextflow est installée sur votre système en utilisant la commande `nextflow --version`.

Pour plus d'informations sur la mise à jour de votre version de Nextflow, veuillez consulter la documentation de référence sur la [Mise à jour de Nextflow](https://www.nextflow.io/docs/latest/updating-nextflow.html).

### Activer l'analyseur syntaxique v2

Pour **activer** l'analyseur syntaxique v2 pour votre session actuelle, exécutez la commande suivante dans votre terminal :

```bash
export NXF_SYNTAX_PARSER=v2
```

Pour rendre cela permanent (en attendant que v2 devienne la valeur par défaut dans Nextflow 26.04), ajoutez la commande export à votre profil shell (`~/.bashrc`, `~/.zshrc`, etc.) :

```bash
echo 'export NXF_SYNTAX_PARSER=v2' >> ~/.bashrc
source ~/.bashrc
```

Notez que la variable d'environnement `NXF_SYNTAX_PARSER=v2` est une exigence temporaire.
À partir de Nextflow 26.04, l'analyseur v2 deviendra la valeur par défaut et ce paramètre ne sera plus nécessaire.

### Désactiver l'analyseur syntaxique v2

Pour **désactiver** l'analyseur syntaxique v2 pour votre session actuelle, exécutez la commande suivante dans votre terminal :

```bash
export NXF_SYNTAX_PARSER=v1
```

<!-- Will it be possible to disable it in versions after 26.04? -->

### Migration du code existant

Pour des conseils concernant la migration du code existant afin de se conformer aux versions plus récentes de Nextflow, veuillez consulter les [Notes de migration](https://www.nextflow.io/docs/latest/migrations/index.html) dans la documentation de référence.

Ces deux articles sont particulièrement utiles pour migrer vers la version la plus récente :

- [Migration vers les sorties de workflow](https://www.nextflow.io/docs/latest/tutorials/workflow-outputs.html)
- [Migration vers les types statiques](https://www.nextflow.io/docs/latest/tutorials/static-types.html)

Ces deux fonctionnalités sont couvertes dans le cadre de la formation pour débutant·es à partir de la version 3.0 des supports de formation.

Selon la génération du code Nextflow que vous avez l'intention de migrer, vous pourrez peut-être en faire la majeure partie avec le linter Nextflow en utilisant la commande `nextflow lint -format`.
Consultez la référence CLI pour [`lint`](https://www.nextflow.io/docs/latest/reference/cli.html#lint) pour plus de détails.

Nous espérons que cela vous sera utile.
Si vous avez besoin d'aide, contactez-nous sur Slack ou sur le forum.

---

<div markdown class="homepage_logos">

![Seqera](../assets/img/seqera_logo.png#only-light)

![Seqera](../assets/img/seqera_logo_dark.png#only-dark)

</div>
