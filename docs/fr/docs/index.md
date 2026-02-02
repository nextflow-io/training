---
title: Accueil
description: Bienvenue sur le portail de formation de la communauté Nextflow !
hide:
  - toc
  - footer
---

# Formation Nextflow

<div class="grid cards" markdown>

-   :material-book-open-variant:{ .lg .middle } __Cours en libre-service__

    ---

    **Bienvenue sur le portail de formation de la communauté Nextflow !**

    Les cours de formation listés ci-dessous sont conçus pour être utilisés en libre-service.
    Vous pouvez les suivre à votre rythme, à tout moment, soit dans l'environnement web que nous fournissons via GitHub Codespaces, soit dans votre propre environnement.

    [Explorer les cours :material-arrow-right:](#catalogue-des-cours-de-formation-nextflow){ .md-button .md-button--primary .mt-1 }

-   :material-information-outline:{ .lg .middle } __Informations complémentaires__

    ---

    ??? warning "Compatibilité des versions"

        <!-- Any update to this content needs to be copied to the local installation page -->
        **Depuis janvier 2026, tous nos cours de formation Nextflow nécessitent Nextflow version 25.10.2 ou ultérieure, avec la syntaxe stricte activée, sauf indication contraire.**

        Pour plus d'informations sur les exigences de version et la syntaxe stricte, veuillez consulter le [guide de migration de la documentation Nextflow](https://nextflow.io/docs/latest/strict-syntax.html).

        Les anciennes versions du matériel de formation correspondant aux syntaxes antérieures sont disponibles via le sélecteur de version dans la barre de menu de cette page web.

    ??? terminal "Options d'environnement"

        Nous fournissons un environnement de formation web où tout ce dont vous avez besoin pour suivre la formation est préinstallé, accessible via GitHub Codespaces (nécessite un compte GitHub gratuit).

        [![Ouvrir dans GitHub Codespaces](https://github.com/codespaces/badge.svg)](https://codespaces.new/nextflow-io/training?quickstart=1&ref=master)

        Si cela ne convient pas à vos besoins, veuillez consulter les autres [options d'environnement](./envsetup/index.md).

    ??? learning "Événements de formation"

        Si vous préférez suivre une formation Nextflow dans le cadre d'un événement structuré, de nombreuses opportunités s'offrent à vous. Nous vous recommandons de consulter les options suivantes :

        - **[Semaines de formation]()** organisées trimestriellement par l'équipe communautaire
        - **[Événements Seqera](https://seqera.io/events/)** incluant des événements de formation en personne organisés par Seqera (recherchez « Seqera Sessions » et « Nextflow Summit »)
        - **[Ambassadeurs Nextflow]()** qui organisent des événements pour leur communauté locale
        - **[Événements nf-core](https://nf-co.re/events)** incluant des hackathons communautaires

    ??? people "Informations pour les formateurs"

        Si vous êtes un·e formateur·trice qui organise vos propres formations, vous êtes invité·e à utiliser nos supports directement depuis le portail de formation à condition d'attribuer les crédits appropriés. Voir « Crédits et contributions » ci-dessous pour plus de détails.

        De plus, nous aimerions avoir de vos nouvelles sur la façon dont nous pourrions mieux soutenir vos efforts de formation ! Veuillez nous contacter à [community@seqera.io](mailto:community@seqera.io) ou sur le forum communautaire (voir la page [Aide](help.md)).

    ??? licensing "Licence open-source et politique de contribution"

        [![Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International (CC BY-NC-SA 4.0)](assets/img/cc_by-nc-sa.svg){ align=right }](https://creativecommons.org/licenses/by-nc-sa/4.0/)

        Ce matériel de formation est développé et maintenu par [Seqera](https://seqera.io) et publié sous licence open-source ([CC BY-NC-SA](https://creativecommons.org/licenses/by-nc-sa/4.0/)) au bénéfice de la communauté. Si vous souhaitez utiliser ce matériel d'une manière qui sort du cadre de la licence (notez les limitations sur l'utilisation commerciale et la redistribution), veuillez nous contacter à [community@seqera.io](mailto:community@seqera.io) pour discuter de votre demande.

        Nous accueillons les améliorations, corrections et rapports de bogues de la communauté. Chaque page possède une icône :material-file-edit-outline: en haut à droite qui renvoie au dépôt de code, où vous pouvez signaler des problèmes ou proposer des modifications au matériel source de formation via une pull request. Consultez le fichier `README.md` dans le dépôt pour plus de détails.

</div>

!!! note "Traduction assistée par IA"

    Cette traduction a été réalisée à l'aide de l'intelligence artificielle et révisée par des traducteurs humains.
    Vos commentaires et suggestions d'amélioration sont les bienvenus.
    Consultez notre [guide de traduction](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md) pour plus d'informations.

## Catalogue des cours de formation Nextflow

<div class="grid cards" markdown>

-   :material-walk:{ .lg .middle } __Parcours d'introduction__

    ---

    ### :material-compass:{.nextflow-primary} Nextflow pour les débutants {.mt-1}

    Cours indépendants du domaine destinés à ceux qui découvrent complètement Nextflow. Chaque cours se compose d'une série de modules de formation conçus pour aider les apprenant·es à développer progressivement leurs compétences.

    ??? courses "**Hello Nextflow :** Apprenez à développer vos propres pipelines"

        Ce cours couvre les composants essentiels du langage Nextflow avec suffisamment de détails pour permettre le développement de pipelines simples mais pleinement fonctionnels, ainsi que les éléments clés des pratiques de conception, de développement et de configuration de pipelines.

        [Commencer la formation Hello Nextflow :material-arrow-right:](hello_nextflow/index.md){ .md-button .md-button--secondary }

    ??? courses "**Nextflow Run :** Apprenez à exécuter des pipelines existants"

        Une introduction concise à l'exécution et à la configuration des pipelines Nextflow, basée sur le cours Hello Nextflow pour développeur·ses mais avec moins d'accent sur le code. Couvre l'exécution, les sorties, la structure de base du code et la configuration pour différents environnements de calcul.

        [Commencer la formation Nextflow Run :material-arrow-right:](nextflow_run/index.md){ .md-button .md-button--secondary }

    ---

    ### :material-microscope:{.nextflow-primary} Nextflow pour la science {.mt-1}

    Apprenez à appliquer les concepts et composants présentés dans « Hello Nextflow » à des cas d'utilisation scientifiques spécifiques.

    ??? courses "**Nextflow pour la génomique** (appel de variants)"

        Pour les chercheur·euses qui souhaitent apprendre à développer leurs propres pipelines génomiques. Ce cours utilise un cas d'utilisation d'appel de variants pour démontrer comment développer un pipeline génomique simple mais fonctionnel.

        [Commencer la formation Nextflow pour la génomique :material-arrow-right:](nf4_science/genomics/){ .md-button .md-button--secondary }

    ??? courses "**Nextflow pour RNAseq** (RNAseq en vrac)"

        Pour les chercheur·euses qui souhaitent apprendre à développer leurs propres pipelines RNAseq. Ce cours utilise un cas d'utilisation de traitement RNAseq en vrac pour démontrer comment développer un pipeline RNAseq simple mais fonctionnel.

        [Commencer la formation Nextflow pour RNAseq :material-arrow-right:](nf4_science/rnaseq/){ .md-button .md-button--secondary }

    ??? courses "**Nextflow pour l'imagerie** (omiques spatiales)"

        Pour les chercheur·euses en imagerie et omiques spatiales qui souhaitent apprendre à exécuter et personnaliser des pipelines d'analyse. Ce cours utilise le pipeline nf-core/molkart pour fournir un pipeline biologiquement pertinent et démontrer comment exécuter, configurer et gérer les entrées pour les flux de travail Nextflow.

        [Commencer la formation Nextflow pour l'imagerie :material-arrow-right:](nf4_science/imaging/){ .md-button .md-button--secondary }

-   :material-run:{ .lg .middle } __Parcours avancé__

    ---

    ### :material-bridge:{.nextflow-primary} De Nextflow à nf-core {.mt-1}

    Apprenez à utiliser le code et les bonnes pratiques du projet communautaire [nf-core](https://nf-co.re/).

    Ces cours vous aident à passer des fondamentaux de Nextflow aux bonnes pratiques nf-core.
    Comprenez comment et pourquoi la communauté nf-core construit des pipelines, et comment vous pouvez contribuer et réutiliser ces techniques.

    ??? courses "**Hello nf-core :** Démarrer avec nf-core"

        Pour les développeur·ses qui souhaitent apprendre à exécuter et développer des pipelines conformes à [nf-core](https://nf-co.re/). Ce cours couvre la structure des pipelines nf-core avec suffisamment de détails pour permettre le développement de pipelines simples mais pleinement fonctionnels qui suivent le modèle nf-core et les bonnes pratiques de développement, ainsi que l'utilisation des modules nf-core existants.

        [Commencer la formation Hello nf-core :material-arrow-right:](hello_nf-core/index.md){ .md-button .md-button--secondary }

    ---

    ### :material-rocket-launch:{.nextflow-primary} Formation Nextflow avancée {.mt-1}

    Apprenez des concepts et mécanismes avancés pour développer et déployer des pipelines Nextflow afin de répondre à des cas d'utilisation réels.

    ??? courses "**Quêtes annexes :** Approfondissements sur des sujets autonomes"

        Mini-cours autonomes destinés aux développeur·ses Nextflow qui souhaitent élargir leur éventail de compétences et/ou approfondir leurs connaissances sur des sujets particuliers. Ils sont présentés de manière linéaire mais peuvent être suivis dans n'importe quel ordre (voir les dépendances dans chaque aperçu de mini-cours).

        [Parcourir les quêtes annexes :material-arrow-right:](side_quests/){ .md-button .md-button--secondary }

    ??? courses "**Collections de formation :** Parcours d'apprentissage recommandés à travers les quêtes annexes"

        Les collections de formation combinent plusieurs quêtes annexes afin de fournir une expérience d'apprentissage complète autour d'un thème ou d'un cas d'utilisation particulier.

        [Parcourir les collections de formation :material-arrow-right:](training_collections/){ .md-button .md-button--secondary }

</div>

!!! info "Vous recherchez des supports de formation archivés ?"

    Les anciens supports de formation (Formation Fondamentale, Formation Avancée et autres cours expérimentaux) ont été retirés du portail de formation car ils sont incompatibles avec la syntaxe stricte de Nextflow 3.0.
    Si vous avez besoin d'accéder à ces supports, ils sont disponibles dans l'[historique git](https://github.com/nextflow-io/training) avant janvier 2026.

---

<div markdown class="homepage_logos">

![Seqera](assets/img/seqera_logo.png#only-light)

![Seqera](assets/img/seqera_logo_dark.png#only-dark)

</div>
