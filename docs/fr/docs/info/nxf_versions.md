---
title: Versions de Nextflow
description: Comprendre et gérer l'évolution des versions de syntaxe de Nextflow
hide:
  - toc
  - footer
---

## Version de syntaxe Nextflow actuellement supportée et exigences

À partir de la version 3.0 du portail de formation, tous nos cours de formation sont basés sur la version 25.10.2 de Nextflow, sauf indication contraire sur la page d'index du cours (à l'exception des matériaux obsolètes ou archivés qui peuvent ne pas inclure d'avis de version).

Étant donné que les cours utilisent maintenant des entrées typées au niveau du workflow ainsi que des directives de sortie au niveau du workflow, ils nécessitent l'utilisation du parser de syntaxe V2.
Si vous prévoyez d'utiliser l'environnement que nous fournissons via [Github Codespaces](../envsetup/01_setup.md) ou [devcontainers locaux](../envsetup/03_devcontainer.md), vous n'avez rien à faire sauf indication spécifique dans les instructions du cours.
Cependant, si vous prévoyez de travailler sur les formations dans votre propre environnement ([Installation manuelle](../envsetup/02_local.md)), vous devrez vous assurer d'utiliser Nextflow version 25.10.2 ou ultérieure avec le parser de syntaxe v2 activé.

## Versions antérieures des matériaux de formation

Nos matériaux de formation sont versionnés depuis février 2025.

Vous pouvez accéder aux versions antérieures des matériaux de formation qui fonctionnent avec les versions de Nextflow **antérieures à 25.10.2** via l'élément de menu déroulant en haut de chaque page qui affiche la version numérotée des matériaux de formation.
Lorsque vous sélectionnez une version antérieure des matériaux de formation, les liens vers l'environnement de formation spécifieront automatiquement la version correspondante de l'environnement.

## Autres informations sur les versions de syntaxe Nextflow

Nextflow a deux concepts de versionnement distincts qui sont parfois confondus : **les versions DSL** et **les versions du parser de syntaxe**.

**DSL1 vs DSL2** fait référence à des façons fondamentalement différentes d'écrire des pipelines Nextflow.
DSL1 était la syntaxe originale où les processes étaient implicitement connectés via des channels.
DSL2, introduit dans Nextflow 20.07, a ajouté des fonctionnalités de modularité : la possibilité d'importer des processes et des workflows depuis d'autres fichiers, des blocs `workflow` explicites et des sorties de process nommées.
DSL1 a été déprécié dans Nextflow 22.03 et supprimé dans 22.12.
Tout le code Nextflow moderne utilise DSL2.

**Parser de syntaxe v1 vs v2** fait référence à différents parsers qui fonctionnent tous deux avec du code DSL2.
Le parser v1 est l'original, plus permissif.
Le parser v2 est plus strict et permet de nouvelles fonctionnalités du langage telles que le typage statique (entrées et sorties typées) et les directives de sortie au niveau du workflow.
Le parser v2 fournit également de meilleurs messages d'erreur et détecte plus d'erreurs au moment de l'analyse plutôt qu'à l'exécution.
Le parser v2 deviendra la valeur par défaut dans Nextflow 26.04.

En résumé : DSL2 est le langage que vous écrivez ; la version du parser de syntaxe détermine à quel point ce langage est interprété strictement et quelles fonctionnalités avancées sont disponibles.

### Vérifier et définir la version de Nextflow

Vous pouvez vérifier quelle version de Nextflow est installée sur votre système en utilisant la commande `nextflow --version`.

Pour plus d'informations sur la mise à jour de votre version de Nextflow, veuillez consulter la documentation de référence sur [Updating Nextflow](https://www.nextflow.io/docs/latest/updating-nextflow.html).

### Activer le parser de syntaxe v2

Pour **activer** le parser de syntaxe v2 pour votre session actuelle, exécutez la commande suivante dans votre terminal :

```bash
export NXF_SYNTAX_PARSER=v2
```

Pour rendre cela permanent (en attendant que v2 devienne la valeur par défaut dans Nextflow 26.04), ajoutez la commande export à votre profil shell (`~/.bashrc`, `~/.zshrc`, etc.) :

```bash
echo 'export NXF_SYNTAX_PARSER=v2' >> ~/.bashrc
source ~/.bashrc
```

Notez que la variable d'environnement `NXF_SYNTAX_PARSER=v2` est une exigence temporaire.
À partir de Nextflow 26.04, le parser v2 deviendra la valeur par défaut et ce paramètre ne sera plus nécessaire.

### Désactiver le parser de syntaxe v2

Pour **désactiver** le parser de syntaxe v2 pour votre session actuelle, exécutez la commande suivante dans votre terminal :

```bash
export NXF_SYNTAX_PARSER=v1
```

<!-- Will it be possible to disable it in versions after 26.04? -->

### Migrer du code existant

Pour des conseils concernant la migration de code existant pour se conformer aux versions plus récentes de Nextflow, veuillez consulter les [Migration Notes](https://www.nextflow.io/docs/latest/migrations/index.html) dans la documentation de référence.

Ces deux articles sont particulièrement utiles pour migrer vers la version la plus récente :

- [Migrating to workflow outputs](https://www.nextflow.io/docs/latest/tutorials/workflow-outputs.html)
- [Migrating to static types](https://www.nextflow.io/docs/latest/tutorials/static-types.html)

Ces deux fonctionnalités sont couvertes dans le cadre de la formation pour débutants à partir de la version 3.0 des matériaux de formation.

Selon la génération de code Nextflow que vous souhaitez migrer, vous pourrez peut-être effectuer la majeure partie du travail avec le linter Nextflow en utilisant la commande `nextflow lint -format`.
Consultez la référence CLI pour [`lint`](https://www.nextflow.io/docs/latest/reference/cli.html#lint) pour plus de détails.

Nous espérons que cela vous sera utile.
Si vous avez besoin d'aide, contactez-nous sur Slack ou sur le forum.

---

<div markdown class="homepage_logos">

![Seqera](../assets/img/seqera_logo.png#only-light)

![Seqera](../assets/img/seqera_logo_dark.png#only-dark)

</div>
