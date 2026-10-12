---
title: Nextflow Run
hide:
  - toc
page_type: index_page
index_type: course
additional_information:
  technical_requirements: true
  learning_objectives:
    - Lancer et gérer des pipelines Nextflow depuis la ligne de commande
    - Comprendre comment les canaux et les opérateurs permettent des workflows multi-entrées et multi-étapes efficaces
    - Utiliser des conteneurs pour gérer les dépendances logicielles et assurer la reproductibilité
    - Configurer l'exécution des pipelines et leurs sorties
    - Générer des rapports d'exécution, inspecter l'historique des exécutions passées et nettoyer les anciens répertoires de travail
    - Exécuter des pipelines directement depuis des dépôts distants tels que GitHub
  audience_prerequisites:
    - "**Public :** Cette formation est conçue pour les apprenant·es qui sont complètement nouveaux·elles avec Nextflow et souhaitent exécuter des pipelines existants."
    - "**Compétences :** Une certaine familiarité avec la ligne de commande, les concepts de base du scripting et les formats de fichiers courants est supposée."
    - "**Domaine :** Les exercices sont tous indépendants du domaine, donc aucune connaissance scientifique préalable n'est requise."
---

# Nextflow Run

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Traduction assistée par IA - [en savoir plus et suggérer des améliorations](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

**Nextflow Run est une introduction pratique à l'exécution de workflows d'analyse de données reproductibles et évolutifs.**

En travaillant sur une série d'exercices orientés vers des objectifs, vous apprendrez les fondamentaux du lancement et de la gestion des pipelines Nextflow, vous comprendrez comment les canaux et les opérateurs permettent le traitement parallèle de plusieurs entrées, et vous utiliserez des conteneurs pour gérer les dépendances logicielles.

Vous en retirerez les compétences et la confiance pour commencer à exécuter des workflows avec Nextflow.

<!-- additional_information -->

## Aperçu de la formation

Cette formation est pratique, avec des exercices orientés vers des objectifs structurés pour introduire les informations progressivement.

Vous exécuterez plusieurs versions d'un pipeline Nextflow qui traite des entrées de texte, en commençant par une version simple constituée d'une seule étape, et en progressant vers une version multi-étapes qui prend un fichier CSV d'entrées, exécute quelques étapes de transformation, et produit un seul fichier texte contenant de l'art ASCII généré par un outil conteneurisé.

Cette formation se concentre sur l'exécution de pipelines (nommée d'après la commande principale `nextflow run`).
Si vous cherchez une introduction au développement de pipelines Nextflow, consultez [Hello Nextflow](../hello_nextflow/index.md).

!!! note "Note"

    Vous cherchez la version précédente de cette formation ? Elle est remplacée par la version sur cette page, mais reste consultable dans la [version 3.6.1](https://training.nextflow.io/3.6.1/nextflow_run/) du site de formation.

### Plan de la formation

| Chapitre de la formation                                                  | Résumé                                                                                                                         | Durée estimée |
| ------------------------------------------------------------------------- | ------------------------------------------------------------------------------------------------------------------------------ | ------------- |
| [Partie 1 : Exécuter Nextflow](./01_run_nextflow.md)                      | Lancer et gérer des pipelines Nextflow, et comprendre les mécanismes essentiels des workflows                                  | 25 mins       |
| [Partie 2 : Configurer le pipeline](./02_configure_pipeline.md)           | Configurer l'exécution du pipeline et ses sorties à l'aide de `nextflow.config`                                                | 20 mins       |
| [Partie 3 : Gérer les exécutions de workflows](./03_manage_executions.md) | Générer des rapports d'exécution, inspecter l'historique des exécutions passées et nettoyer les anciens répertoires de travail | 10 mins       |
| [Partie 4 : Exécuter des pipelines distants](./04_remote_repositories.md) | Exécuter un pipeline directement depuis GitHub et le fixer à une révision spécifique                                           | 10 mins       |

À la fin de cette formation, vous serez bien préparé·e pour aborder les prochaines étapes de votre parcours pour exécuter des workflows reproductibles pour vos besoins de calcul scientifique.

Prêt·e à suivre la formation ?

[Commencer l'apprentissage :material-arrow-right:](00_orientation.md){ .md-button .md-button--primary }
