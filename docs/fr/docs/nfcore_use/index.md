---
title: Utiliser nf-core
hide:
  - toc
page_type: index_page
index_type: course
additional_information:
  technical_requirements: true
  learning_objectives:
    - Trouver, récupérer et exécuter des pipelines communautaires nf-core
    - Configurer l'exécution d'un pipeline à l'aide de paramètres et de fichiers de configuration
    - Comprendre comment les pipelines nf-core valident les paramètres et les données d'entrée
    - Exécuter un pipeline en production (nf-core/rnaseq) et remplacer ses allocations de ressources par défaut
  audience_prerequisites:
    - "**Public :** Ce cours est destiné aux apprenant·es qui savent déjà exécuter des pipelines Nextflow en local et qui découvrent nf-core, et souhaitent exécuter des pipelines communautaires existants."
    - "**Compétences :** Une certaine familiarité avec la ligne de commande, les concepts de base du scripting et les formats de fichiers courants est supposée."
    - "**Cours :** Avoir complété [Nextflow Run](../nextflow_run/index.md) ou être à l'aise avec l'exécution d'un pipeline local avec `nextflow run`."
    - "**Domaine :** Les exercices utilisent des pipelines bioinformatiques, mais aucune connaissance scientifique préalable n'est requise."
---

# Utiliser nf-core

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Traduction assistée par IA - [en savoir plus et suggérer des améliorations](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

**Utiliser nf-core est une introduction pratique à la recherche, l'exécution et la configuration de pipelines communautaires nf-core.**

En travaillant à travers des exemples concrets et des exercices guidés, vous apprendrez à trouver et récupérer des pipelines nf-core, à les exécuter à l'aide de leurs profils de test intégrés, et à personnaliser leur exécution via des paramètres et des fichiers de configuration.

Vous repartirez avec les compétences et la confiance nécessaires pour commencer à exécuter des pipelines nf-core pour vos propres analyses.

<!-- additional_information -->

## Aperçu du cours

Ce cours est pratique, avec des exercices orientés vers des objectifs, structurés pour introduire les informations progressivement.

Vous commencerez avec `nf-core/demo`, un pipeline minimal maintenu par le projet nf-core à des fins de formation, puis vous appliquerez ce que vous avez appris à `nf-core/rnaseq`, un pipeline de production largement utilisé pour l'analyse du séquençage RNA en masse.

Ce cours se concentre sur l'exécution de pipelines.
Si vous cherchez une introduction au développement de pipelines compatibles nf-core, consultez [Build with nf-core](../nfcore_build/index.md).

### Plan de cours

| Chapitre du cours                                                                | Résumé                                                                                                                     | Durée estimée |
| -------------------------------------------------------------------------------- | -------------------------------------------------------------------------------------------------------------------------- | ------------- |
| [Partie 1 : Exécuter un pipeline de démonstration](./01_run_demo.md)             | Trouver et récupérer un pipeline nf-core et l'exécuter à l'aide de son profil de test                                      | 20 min        |
| [Partie 2 : Configurer l'exécution du pipeline](./02_configure_execution.md)     | Définir des paramètres, comprendre la validation, et personnaliser l'allocation des ressources et les arguments des outils | 20 min        |
| [Partie 3 : Exécuter un pipeline en production](./03_run_production_pipeline.md) | Récupérer et exécuter nf-core/rnaseq, et remplacer ses allocations de ressources par défaut                                | 20 min        |

À la fin de ce cours, vous serez en mesure de tirer parti de la richesse des pipelines communautaires offerts par le projet nf-core.

Prêt·e à suivre le cours ?

[Commencer l'apprentissage :material-arrow-right:](00_orientation.md){ .md-button .md-button--primary }
