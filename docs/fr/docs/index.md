---
title: Accueil
description: Bienvenue sur le portail de formation de la communauté Nextflow !
hide:
  - toc
  - footer
---

# Formation Nextflow

<div class="grid cards" markdown>

-   :material-book-open-variant:{ .lg .middle } __Cours en libre accès__

    ---

    **Bienvenue sur le portail de formation de la communauté Nextflow !**

    Suivez les cours ci-dessous à votre propre rythme, dans notre environnement web ou le vôtre.
    Chaque cours est pratique, avec des exercices orientés vers des objectifs que vous pouvez réaliser de manière autonome.

    [Explorer les cours :material-arrow-down:](#catalog-of-nextflow-training-courses){ .md-button .md-button--primary .mt-1 }

-   :material-account-group-outline:{ .lg .middle } __Événements de formation__

    ---

    **Vous cherchez quelque chose au-delà du libre accès ?**

    Retrouvez des événements de formation structurés, des conseils pour organiser vos propres formations, ainsi que notre licence open-source et notre politique de contribution.

    [Voir les événements de formation :material-arrow-right:](training_events.md){ .md-button .md-button--secondary .mt-1 }

</div>

!!! note "Traduction assistée par IA"

    Cette traduction a été créée à l'aide de l'intelligence artificielle et révisée par des traducteurs humains.
    Nous apprécions vos commentaires et suggestions d'amélioration.
    Consultez notre [guide de traduction](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md) pour plus d'informations.

## Catalog of Nextflow training courses

<div class="grid cards" markdown>

-   :material-account:{ .lg .middle } __Pour les utilisateur·trices__

    ---

    ### :material-play-circle:{.nextflow-primary} Exécuter des pipelines {.mt-1}

    Apprenez à exécuter des pipelines existants sans écrire de code.

    ??? courses "**Nextflow Run :** Exécuter des pipelines avec Nextflow"

        Une introduction rapide à l'exécution de pipelines Nextflow qui ne nécessite pas de comprendre le code. Couvre le lancement de pipelines, la récupération des sorties, l'utilisation de conteneurs et la configuration de l'exécution à un niveau de base.

        [Voir la formation :material-arrow-right:](nextflow_run/index.md){ .md-button .md-button--secondary }

    ??? courses "**Use nf-core :** Trouver et exécuter des pipelines créés par la communauté"

        Une introduction rapide à la recherche, l'exécution et la configuration de pipelines issus du projet communautaire nf-core, en commençant par un pipeline de démonstration minimal puis en passant à un pipeline d'analyse à l'échelle de la production.

        [Voir la formation :material-arrow-right:](nfcore_use/index.md){ .md-button .md-button--secondary }

    ??? courses "**Scale with Seqera :** Lancer et surveiller des pipelines à grande échelle"

        Une introduction pratique au lancement et à la surveillance de pipelines Nextflow avec Seqera Platform, depuis l'interface web et la ligne de commande.

        [Voir la formation :material-arrow-right:](seqera_scale/index.md){ .md-button .md-button--secondary }

    ---

    ### :material-tune:{.nextflow-primary} Gérer l'exécution {.mt-1}

    Apprenez à gérer efficacement l'exécution des pipelines.

    ??? courses "**Configure Execution :** Configurer les ressources, les nouvelles tentatives et les profils d'exécution"

        Une introduction pratique à la configuration de l'exécution de pipelines Nextflow : adaptation à différents environnements de calcul, contrôle des allocations de ressources et des nouvelles tentatives, et basculement entre des profils de configuration prédéfinis.

        [Voir la formation :material-arrow-right:](config_exec/index.md){ .md-button .md-button--secondary }

    !!! info compact "D'autres sujets à venir"

        L'optimisation des performances, l'exécution HPC/cloud et bien d'autres sujets sont prévus pour cette section.
        Votez pour les prochains sujets à couvrir dans notre [court sondage](https://seqera.typeform.com/to/JCs91e8v).

-   :material-code-tags:{ .lg .middle } __Pour les développeur·ses__

    ---

    ### :material-wrench:{.nextflow-primary} Écrire des pipelines {.mt-1}

    Apprenez à développer vos propres pipelines Nextflow.

    ??? courses "**Hello Nextflow :** Développer vos propres pipelines de zéro"

        Ce cours couvre les composants essentiels du langage Nextflow avec suffisamment de détails pour permettre le développement de pipelines simples mais entièrement fonctionnels, ainsi que les éléments clés de la conception, du développement et des bonnes pratiques de configuration des pipelines.

        [Voir la formation :material-arrow-right:](hello_nextflow/index.md){ .md-button .md-button--secondary }

    ??? courses "**Build with nf-core :** Utiliser les outils et les règles nf-core"

        Pour les développeur·ses Nextflow qui souhaitent apprendre à développer des pipelines conformes à [nf-core](https://nf-co.re/).
        Le cours couvre la structure des pipelines nf-core avec suffisamment de détails pour permettre le développement de pipelines simples mais entièrement fonctionnels qui tirent parti du template nf-core et des meilleures pratiques de développement, ainsi que l'utilisation de modules nf-core existants.

        [Voir la formation :material-arrow-right:](nfcore_build/index.md){ .md-button .md-button--secondary }

    ??? catalog "**Side Quests :** Plongez dans des sujets Nextflow avancés"

        Une collection de mini-cours autonomes destinés aux développeur·ses Nextflow qui souhaitent élargir leur champ de compétences et/ou approfondir leurs connaissances sur des sujets particuliers.
        Ils sont présentés de manière linéaire mais peuvent être suivis dans n'importe quel ordre (voir les dépendances dans l'aperçu de chaque mini-cours).

        [Parcourir les Quêtes secondaires :material-arrow-right:](side_quests/index.md){ .md-button .md-button--secondary }

    ---

    ### :material-microscope:{.nextflow-primary} Nextflow pour la science {.mt-1}

    Apprenez à développer des pipelines Nextflow pour des applications scientifiques spécifiques.

    ??? courses "**Genomics :** Développer un pipeline de détection de variants"

        Un cours pour les chercheur·euses qui souhaitent apprendre à développer leurs propres pipelines de génomique, en utilisant un cas d'usage de détection de variants pour illustrer les patterns essentiels de développement Nextflow.

        [Voir la formation :material-arrow-right:](nf4_science/genomics/index.md){ .md-button .md-button--secondary }

    ??? courses "**RNAseq :** Développer un pipeline de traitement RNAseq en masse"

        Un cours pour les chercheur·euses qui souhaitent apprendre à développer leurs propres pipelines RNAseq, en utilisant un cas d'usage de traitement RNAseq en masse pour illustrer les patterns essentiels de développement Nextflow.

        [Voir la formation :material-arrow-right:](nf4_science/rnaseq/index.md){ .md-button .md-button--secondary }

    ??? courses "**Bioimaging :** Exécuter et configurer des pipelines d'imagerie"

        Un cours pour les chercheur·euses qui souhaitent apprendre à exécuter et configurer des pipelines de bio-imagerie, en utilisant nf-core/molkart pour illustrer les patterns essentiels d'utilisation de Nextflow.

        [Voir la formation :material-arrow-right:](nf4_science/imaging/index.md){ .md-button .md-button--secondary }

</div>

## Configuration & Aide

<div class="grid cards mb-4" markdown>

-   :material-cog-outline:{ .lg .middle } __Environnement de formation__

    ---

    Options pour configurer votre environnement pour les formations Nextflow.

    [Voir les environnements de formation :material-arrow-right:](envsetup/index.md){ .md-button .md-button--secondary }

-   :material-tag-outline:{ .lg .middle } __Versions de Nextflow__

    ---

    Comprendre et gérer l'évolution des versions de syntaxe de Nextflow.

    [Vérifier les exigences de version :material-arrow-right:](info/nxf_versions.md){ .md-button .md-button--secondary }

-   :material-file-code-outline:{ .lg .middle } __Le pipeline Hello__

    ---

    Récapitulatif de ce que fait le pipeline Hello et de sa structure.

    [Lire le récapitulatif :material-arrow-right:](info/hello_pipeline.md){ .md-button .md-button--secondary }

-   :material-lifebuoy:{ .lg .middle } __Obtenir de l'aide__

    ---

    Ressources utiles lorsque vous rencontrez un problème avec la formation Nextflow.

    [Trouver de l'aide :material-arrow-right:](help.md){ .md-button .md-button--secondary }

</div>
