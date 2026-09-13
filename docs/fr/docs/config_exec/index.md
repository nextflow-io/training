---
title: Configure Execution
hide:
  - toc
page_type: index_page
index_type: course
additional_information:
  technical_requirements: true
  learning_objectives:
    - Changer de technologie de packaging logiciel entre Docker et Conda
    - Sélectionner une plateforme d'exécution et comprendre comment Nextflow adapte l'exécution des tâches
    - Contrôler les allocations de ressources de calcul et relancer automatiquement les tâches qui échouent
    - Définir et combiner des profils pour basculer entre des configurations prédéfinies
  audience_prerequisites:
    - "**Audience :** Ce cours est conçu pour les apprenant·es qui savent déjà lancer des pipelines Nextflow en local et souhaitent configurer l'exécution de manière plus approfondie."
    - "**Compétences :** Une certaine familiarité avec la ligne de commande est supposée."
    - "**Cours :** Vous devez avoir suivi [Nextflow Run](../nextflow_run/index.md) ou être à l'aise avec l'exécution d'un pipeline local via `nextflow run`."
---

# Configure Execution

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Traduction assistée par IA - [en savoir plus et suggérer des améliorations](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>


**Configure Execution est une introduction pratique à l'adaptation de l'exécution de pipelines Nextflow à différents environnements de calcul.**

À travers des exercices orientés objectifs, vous apprendrez à changer de technologie de packaging logiciel, à sélectionner une plateforme d'exécution, à contrôler les allocations de ressources de calcul et les relances, et à regrouper la configuration en profils commutables.

Vous repartirez avec les compétences et la confiance nécessaires pour configurer l'exécution de pipelines Nextflow comme un·e pro.

<!-- additional_information -->

## Aperçu du cours

Ce cours est pratique et s'appuie sur les compétences abordées dans [Nextflow Run](../nextflow_run/index.md).

Vous reprendrez le même pipeline multi-étapes de ce cours et adapterez progressivement sa configuration à différents environnements de calcul, puis regrouperez tout dans des profils entre lesquels vous pourrez basculer au moment de l'exécution.

### Plan de cours

| Chapitre du cours                                                                        | Résumé                                                                                      | Durée estimée |
| ---------------------------------------------------------------------------------------- | ------------------------------------------------------------------------------------------- | ------------- |
| [Partie 1 : S'adapter à votre environnement de calcul](./01_packaging_and_execution.md)  | Changer de technologie de packaging logiciel et sélectionner une plateforme d'exécution     | 20 min        |
| [Partie 2 : Gérer les ressources de calcul et les échecs](./02_resources_and_retries.md) | Contrôler les allocations de ressources et relancer automatiquement les tâches qui échouent | 15 min        |
| [Partie 3 : Utiliser des profils pour changer de configuration](./03_profiles.md)        | Définir et combiner des profils, et inspecter la configuration entièrement résolue          | 15 min        |

À la fin de ce cours, vous serez à l'aise pour configurer des pipelines Nextflow pour différents environnements de calcul et basculer entre eux sans difficulté.

Prêt·e à suivre le cours ?

[Commencer l'apprentissage :material-arrow-right:](00_orientation.md){ .md-button .md-button--primary }
