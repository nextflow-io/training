---
title: Options d'environnement
description: Options pour configurer votre environnement pour les formations Nextflow
hide:
  - toc
  - footer
---

# Options d'environnement

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Traduction assistée par IA - [en savoir plus et suggérer des améliorations](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

Nous visons à fournir un environnement cohérent et minutieusement testé qui permet aux apprenant·es de se concentrer sur l'apprentissage de Nextflow sans avoir à consacrer du temps et des efforts à la gestion des logiciels.
À cette fin, nous avons développé un environnement conteneurisé qui contient tous les logiciels, fichiers de code et données d'exemple nécessaires pour suivre tous nos cours.

Cet environnement conteneurisé peut être exécuté directement sur GitHub Codespaces ou localement dans VS Code avec l'extension Devcontainers.

<div class="grid cards" markdown>

- :material-cloud-outline:{ .lg .middle } **GitHub Codespaces**

  ***

  GitHub Codespaces est un service web qui nous permet de fournir un environnement pré-construit pour la formation, avec tous les outils et données inclus, soutenu par des machines virtuelles dans le cloud. Il est accessible gratuitement à toute personne possédant un compte GitHub.

  [Utiliser GitHub Codespaces :material-arrow-right:](01_setup.md){ .md-button .md-button--primary .mt-1 }

- :material-laptop:{ .lg .middle } **Devcontainers locaux**

  ***

  VS Code avec Devcontainers fournit un environnement de développement conteneurisé exécuté localement avec tous les outils de formation préconfigurés. Il offre le même environnement pré-construit que Codespaces mais s'exécute entièrement sur votre matériel local.

  [Utiliser Devcontainers localement :material-arrow-right:](03_devcontainer.md){ .md-button .md-button--primary .mt-1 }

</div>

## Instructions pour l'installation manuelle

Si aucune des options ci-dessus ne convient à vos besoins, vous pouvez reproduire cet environnement sur votre propre système local en installant manuellement les dépendances logicielles et en clonant le dépôt de formation.

[Installation manuelle :material-arrow-right:](02_local.md){ .md-button .md-button--primary .mt-1 }

---

!!! info "Dépréciation de Gitpod"

    La formation Nextflow utilisait [Gitpod](https://gitpod.io) jusqu'en février 2025.
    Cependant, les créateurs de Gitpod ont décidé de retirer la fonctionnalité gratuite au profit du système [Gitpod Flex](https://www.gitpod.io/blog/introducing-gitpod-flex).
    Pour cette raison, nous sommes passés à l'utilisation de GitHub Codespaces, qui offre également un environnement de développement en un clic sans configuration préalable.

    Selon le moment où vous vous êtes inscrit·e à Gitpod et le moment exact où ils retireront le service, vous pourriez encore être en mesure de lancer la formation dans leur ancien IDE cloud, bien que nous ne puissions pas garantir un accès fiable à l'avenir :
    [Ouvrir dans Gitpod](https://gitpod.io/#https://github.com/nextflow-io/training).
