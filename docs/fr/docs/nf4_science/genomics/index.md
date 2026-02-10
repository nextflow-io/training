---
title: Nextflow for Genomics
hide:
  - toc
page_type: index_page
index_type: course
additional_information:
  technical_requirements: true
  learning_objectives:
    - Write a linear workflow to apply variant calling to a single sample
    - Handle accessory files such as index files and reference genome resources appropriately
    - Leverage Nextflow's dataflow paradigm to parallelize per-sample variant calling
    - Implement multi-sample joint calling using relevant channel operators
  audience_prerequisites:
    - "**Audience:** This course is designed for researchers in genomics and related fields who want to develop or customize data analysis pipelines."
    - "**Skills:** Some familiarity with the command line, basic scripting concepts, and common genomics file formats is assumed."
    - "**Prerequisites:** Foundational Nextflow concepts and tooling covered in [Hello Nextflow](../../hello_nextflow/)."
---

# Nextflow pour la génomique

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Traduction assistée par IA - [en savoir plus et suggérer des améliorations](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

**Un cours pratique appliquant Nextflow à un cas d'usage réel en génomique : l'appel de variants avec GATK.**

Ce cours s'appuie sur la formation pour débutant·es [Hello Nextflow](../../hello_nextflow/) et démontre comment utiliser Nextflow dans le contexte spécifique du domaine de la génomique.
Vous implémenterez un pipeline d'appel de variants avec [GATK](https://gatk.broadinstitute.org/) (Genome Analysis Toolkit), un logiciel largement utilisé pour l'analyse de données de séquençage à haut débit.

<!-- additional_information -->

## Aperçu du cours

Ce cours est pratique, avec des exercices orientés vers des objectifs et structurés pour introduire progressivement les informations.

Vous commencerez par exécuter les outils d'appel de variants manuellement dans le terminal pour comprendre la méthodologie, puis vous construirez progressivement un pipeline Nextflow qui automatise et met à l'échelle l'analyse.

### Plan de cours

Nous avons divisé ce cours en trois parties qui se concentrent chacune sur des aspects spécifiques de l'application de Nextflow à un cas d'usage en génomique.

| Chapitre du cours                                                                  | Résumé                                                                                                                             | Durée estimée |
| ---------------------------------------------------------------------------------- | ---------------------------------------------------------------------------------------------------------------------------------- | ------------- |
| [Partie 1 : Aperçu de la méthode](./01_method.md)                                  | Comprendre la méthodologie d'appel de variants et exécuter les outils manuellement                                                 | 30 min        |
| [Partie 2 : Appel de variants par échantillon](./02_per_sample_variant_calling.md) | Construire un pipeline qui indexe les fichiers BAM et appelle les variants, puis le mettre à l'échelle pour plusieurs échantillons | 60 min        |
| [Partie 3 : Appel conjoint sur une cohorte](./03_joint_calling.md)                 | Ajouter le génotypage conjoint multi-échantillons en utilisant des opérateurs de canaux pour agréger les sorties par échantillon   | 45 min        |

À la fin de ce cours, vous serez capable d'appliquer les concepts et outils fondamentaux de Nextflow à un cas d'usage typique en génomique.

Prêt·e à suivre le cours ?

[Commencer :material-arrow-right:](00_orientation.md){ .md-button .md-button--primary }

<!-- Clearfix for float -->
<div style="content: ''; clear: both; display: table;"></div>
