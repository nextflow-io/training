---
title: Nextflow for {DOMAIN}
hide:
  - toc
page_type: index_page
index_type: course
additional_information:
  technical_requirements: true
  learning_objectives:
    - Write a linear workflow to apply {METHOD} to a single sample
    - Handle accessory files such as {ACCESSORY_FILES} appropriately
    - Leverage Nextflow's dataflow paradigm to parallelize per-sample processing
    - Implement multi-sample aggregation using relevant channel operators
  audience_prerequisites:
    - "**Audience:** This course is designed for researchers in {DOMAIN} and related fields who want to develop or customize data analysis pipelines."
    - "**Skills:** Some familiarity with the command line, basic scripting concepts, and common {DOMAIN} file formats is assumed."
    - "**Prerequisites:** Foundational Nextflow concepts and tooling covered in [Hello Nextflow](../../hello_nextflow/)."
---

# Nextflow pour {DOMAIN}

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Traduction assistée par IA - [en savoir plus et suggérer des améliorations](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

**Un cours pratique appliquant Nextflow à un cas d'usage réel en {DOMAIN} : {METHOD_SHORT_DESCRIPTION}.**

Ce cours s'appuie sur la formation pour débutant·es [Hello Nextflow](../../hello_nextflow/) et démontre comment utiliser Nextflow dans le contexte spécifique du domaine {DOMAIN}.
Vous implémenterez un pipeline {METHOD} avec [{TOOL_A}]({TOOL_A_URL}) et [{TOOL_B}]({TOOL_B_URL}).

<!-- additional_information -->

## Aperçu du cours

Ce cours est pratique, avec des exercices orientés vers des objectifs et structurés pour introduire les informations progressivement.

Vous commencerez par exécuter les outils d'analyse manuellement dans le terminal pour comprendre la méthodologie, puis vous construirez progressivement un pipeline Nextflow qui automatise et met à l'échelle l'analyse.

### Plan de cours

Nous avons divisé ce cours en trois parties qui se concentrent chacune sur des aspects spécifiques de l'application de Nextflow à un cas d'usage en {DOMAIN}.

| Chapitre du cours                                         | Résumé                                                                                                                  | Durée estimée |
| --------------------------------------------------------- | ----------------------------------------------------------------------------------------------------------------------- | ------------- |
| [Partie 1 : Aperçu de la méthode](./01_method.md)         | Comprendre la méthodologie {METHOD} et exécuter les outils manuellement                                                 | 30 min        |
| [Partie 2 : Traitement mono-échantillon](./02_single_sample.md) | Construire un pipeline qui {PART2_SUMMARY}, puis le mettre à l'échelle pour plusieurs échantillons                     | 60 min        |
| [Partie 3 : Agrégation multi-échantillons](./03_multi_sample.md) | Ajouter l'agrégation multi-échantillons {AGGREGATION_SUMMARY} en utilisant des opérateurs de canaux pour agréger les sorties par échantillon | 45 min        |

À la fin de ce cours, vous serez capable d'appliquer les concepts et outils fondamentaux de Nextflow à un cas d'usage typique en {DOMAIN}.

Prêt·e à suivre le cours ?

[Commencer :material-arrow-right:](00_orientation.md){ .md-button .md-button--primary }

<!-- Clearfix for float -->
<div style="content: ''; clear: both; display: table;"></div>
