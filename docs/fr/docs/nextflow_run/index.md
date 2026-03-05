---
title: Nextflow Run
hide:
  - toc
page_type: index_page
index_type: course
additional_information:
  technical_requirements: true
  learning_objectives:
    - Lancer et gérer l'exécution de workflows Nextflow
    - Trouver et interpréter les sorties (résultats) et les fichiers de log
    - Reconnaître les composants principaux de Nextflow dans un workflow multi-étapes simple
    - Configurer l'exécution d'un pipeline pour fonctionner sur des plateformes de calcul courantes, y compris HPC et cloud
    - Résumer les bonnes pratiques pour la reproductibilité, la portabilité et la réutilisation du code qui rendent les pipelines FAIR, y compris la modularité du code et les conteneurs logiciels
  audience_prerequisites:
    - "**Public :** Cette formation est conçue pour les apprenant·es qui sont complètement nouveaux·elles avec Nextflow et souhaitent exécuter des pipelines existants."
    - "**Compétences :** Une certaine familiarité avec la ligne de commande, les concepts de base du scripting et les formats de fichiers courants est supposée."
    - "**Domaine :** Les exercices sont tous indépendants du domaine, donc aucune connaissance scientifique préalable n'est requise."
---

# Nextflow Run

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Traduction assistée par IA - [en savoir plus et suggérer des améliorations](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

**Nextflow Run est une introduction pratique à l'exécution de workflows d'analyse de données reproductibles et évolutifs.**

En travaillant sur des exemples pratiques et des exercices guidés, vous apprendrez les fondamentaux de l'utilisation de Nextflow, notamment comment exécuter des pipelines, gérer les fichiers et les dépendances logicielles, paralléliser l'exécution sans effort, et exécuter des workflows dans différents environnements de calcul.

Vous en retirerez les compétences et la confiance pour commencer à exécuter des workflows avec Nextflow.

<!-- additional_information -->

## Aperçu de la formation

### Ce que vous allez faire

Cette formation est pratique, avec des exercices orientés vers des objectifs structurés pour introduire les informations progressivement.

Vous exécuterez plusieurs versions d'un pipeline Nextflow qui traite des entrées de texte.
Vous commencerez par une version simple constituée d'une seule étape, et progresserez éventuellement vers une version multi-étapes qui prend un fichier CSV de données textuelles tabulaires, exécute quelques étapes de transformation, et produit un seul fichier texte contenant une image ASCII d'un personnage disant le texte transformé.

Cette formation se concentre sur l'exécution de pipelines (nommée d'après la commande principale `nextflow run`).
Si vous cherchez une introduction au développement de pipelines Nextflow, consultez [Hello Nextflow](../hello_nextflow/index.md).

### Plan de la formation

Nous avons divisé cela en trois parties qui se concentreront chacune sur des aspects spécifiques de l'exécution et de la gestion des pipelines écrits en Nextflow.

| Chapitre de la formation                                     | Résumé                                                                                                                              | Durée estimée |
| ------------------------------------------------------------ | ----------------------------------------------------------------------------------------------------------------------------------- | ------------- |
| [Partie 1 : Exécuter les opérations de base](./01_basics.md) | Lancer et gérer l'exécution d'un workflow simple                                                                                    | 30 mins       |
| [Partie 2 : Exécuter de vrais pipelines](./02_pipeline.md)   | Traiter des entrées complexes, exécuter des workflows multi-étapes, utiliser des conteneurs et paralléliser l'exécution sans effort | 60 mins       |
| [Partie 3 : Configuration d'exécution](./03_config.md)       | Personnaliser le comportement du pipeline et optimiser l'utilisation dans différents environnements de calcul                       | 60 mins       |

À la fin de cette formation, vous serez bien préparé·e pour aborder les prochaines étapes de votre parcours pour exécuter des workflows reproductibles pour vos besoins de calcul scientifique.

Prêt·e à suivre la formation ?

[Commencer l'apprentissage :material-arrow-right:](00_orientation.md){ .md-button .md-button--primary }
