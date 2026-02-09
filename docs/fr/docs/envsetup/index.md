# Options d'environnement

Nous visons à fournir un environnement cohérent et rigoureusement testé qui permet aux apprenant·es de se concentrer sur l'apprentissage de Nextflow sans avoir à consacrer du temps et des efforts à la gestion des logiciels.
À cette fin, nous avons développé un environnement conteneurisé qui contient tous les logiciels nécessaires, les fichiers de code et les données d'exemple pour suivre tous nos cours.

Cet environnement conteneurisé peut être exécuté directement sur Github Codespaces ou localement dans VS Code avec l'extension Devcontainers.

<div class="grid cards" markdown>

-   :material-cloud-outline:{ .lg .middle } __Github Codespaces__

    ---

    GitHub Codespaces est un service web qui nous permet de fournir un environnement pré-configuré pour la formation, avec tous les outils et données inclus, soutenu par des machines virtuelles dans le cloud. Il est accessible gratuitement à toute personne disposant d'un compte Github.

    [Utiliser Github Codespaces :material-arrow-right:](01_setup.md){ .md-button .md-button--primary .mt-1 }

-   :material-laptop:{ .lg .middle } __Devcontainers locaux__

    ---

    VS Code avec Devcontainers fournit un environnement de développement conteneurisé exécuté localement avec tous les outils de formation pré-configurés. Il offre le même environnement pré-configuré que Codespaces mais s'exécute entièrement sur votre matériel local.

    [Utiliser Devcontainers localement :material-arrow-right:](03_devcontainer.md){ .md-button .md-button--primary .mt-1 }

</div>

## Instructions pour l'installation manuelle

Si aucune des options ci-dessus ne convient à vos besoins, vous pouvez reproduire cet environnement sur votre propre système local en installant les dépendances logicielles manuellement et en clonant le dépôt de formation.

[Installation manuelle :material-arrow-right:](02_local.md){ .md-button .md-button--primary .mt-1 }

---

!!! info "Abandon de Gitpod"

    Nextflow Training utilisait [Gitpod](https://gitpod.io) jusqu'en février 2025.
    Cependant, les créateurs de Gitpod ont décidé de retirer la fonctionnalité gratuite au profit du système [Gitpod Flex](https://www.gitpod.io/blog/introducing-gitpod-flex).
    Pour cette raison, nous sommes passés à l'utilisation de GitHub Codespaces, qui offre également un environnement de développement en un clic sans configuration préalable.

    Selon la date à laquelle vous vous êtes inscrit·e à Gitpod et le moment exact où ils retireront le service, vous pourrez peut-être encore lancer la formation dans leur ancien IDE cloud, bien que nous ne puissions garantir un accès fiable à l'avenir :
    [Ouvrir dans Gitpod](https://gitpod.io/#https://github.com/nextflow-io/training).
