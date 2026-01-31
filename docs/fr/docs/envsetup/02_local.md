# Installation manuelle

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Traduction assistée par IA - [en savoir plus et suggérer des améliorations](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

Il est possible d'installer tout ce dont vous avez besoin pour suivre la formation dans votre propre environnement local manuellement.

Nous avons documenté ici comment procéder sur les systèmes compatibles POSIX standard (en supposant une machine personnelle comme un ordinateur portable).
Gardez à l'esprit que certains détails peuvent différer selon votre système spécifique.

!!! tip "Astuce"

    Avant de continuer, avez-vous envisagé d'utiliser l'[approche Devcontainers](03_devcontainer.md) ?
    Elle fournit tous les outils et dépendances nécessaires sans nécessiter d'installation manuelle.

## Exigences logicielles générales

Nextflow peut être utilisé sur tout système compatible POSIX (Linux, macOS, Windows Subsystem for Linux, etc.) avec Java installé.
Nos cours de formation ont quelques exigences supplémentaires.

Au total, vous devrez avoir les logiciels suivants installés :

- Bash ou shell équivalent
- [Java 11 (ou plus récent, jusqu'à 21)](https://www.oracle.com/technetwork/java/javase/downloads/index.html)
- [Git](https://git-scm.com/)
- [Docker](https://docs.docker.com/get-docker/)
- [Conda](https://conda.io/) 4.5 (ou plus récent)
- [VSCode](https://code.visualstudio.com) avec l'[extension Nextflow](https://www.nextflow.io/docs/latest/developer-env.html#devenv-nextflow)

L'application VSCode est techniquement optionnelle mais nous vous recommandons fortement de l'utiliser pour suivre les cours ainsi que pour votre travail de développement Nextflow en général.

Le manuel de documentation Nextflow fournit des instructions pour installer ces dépendances sous [Configuration de l'environnement](https://www.nextflow.io/docs/latest/developer-env.html).

## Nextflow et outils nf-core

Vous devrez installer Nextflow lui-même, ainsi que les outils nf-core, comme détaillé dans les articles liés ci-dessous :

- [Installation de Nextflow](https://www.nextflow.io/docs/latest/install.html)
- [Outils nf-core](https://nf-co.re/docs/nf-core-tools/installation)

Nous recommandons d'utiliser l'option d'auto-installation pour Nextflow et l'option PyPI pour les outils nf-core.

!!! warning "Compatibilité des versions"

    <!-- Any update to this content needs to be copied to the home page -->
    **Depuis janvier 2026, tous nos cours de formation Nextflow nécessitent Nextflow version 25.10.2 ou ultérieure, avec la syntaxe stricte v2 activée, sauf indication contraire.**

    Pour plus d'informations sur les exigences de version et la syntaxe stricte v2, veuillez consulter le guide [Versions de Nextflow](../info/nxf_versions.md).

    Les anciennes versions du matériel de formation correspondant aux syntaxes antérieures sont disponibles via le sélecteur de version dans la barre de menu de cette page web.

## Supports de formation

Le moyen le plus simple de télécharger les supports de formation est de cloner l'intégralité du dépôt en utilisant cette commande :

```bash
git clone https://github.com/nextflow-io/training.git
```

Chaque cours a son propre répertoire.
Pour suivre un cours, ouvrez une fenêtre de terminal (idéalement depuis l'application VSCode) et faites `cd` dans le répertoire approprié.

Vous pouvez ensuite suivre les instructions du cours fournies sur le site web.
