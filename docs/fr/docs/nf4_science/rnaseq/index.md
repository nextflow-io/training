---
title: ARNseq
hide:
  - toc
page_type: index_page
index_type: course
additional_information:
  technical_requirements: true
  learning_objectives:
    - Écrire un workflow linéaire pour appliquer des méthodes de base de traitement et de QC d'ARNseq
    - Gérer de manière appropriée les fichiers spécifiques au domaine tels que les FASTQ et les ressources de génome de référence
    - Gérer les données de séquençage single-end et paired-end
    - Exploiter le paradigme de flux de données de Nextflow pour paralléliser le traitement d'ARNseq par échantillon
    - Agréger les rapports de QC à travers plusieurs étapes et échantillons en utilisant les opérateurs de canal pertinents
  audience_prerequisites:
    - "**Public :** Ce cours est destiné aux chercheur·euses en transcriptomique et dans les domaines connexes qui souhaitent développer ou personnaliser des pipelines d'analyse de données."
    - "**Compétences :** Une certaine familiarité avec la ligne de commande, les concepts de base du scripting et les formats de fichiers ARNseq courants est supposée."
    - "**Prérequis :** Concepts et outils fondamentaux de Nextflow couverts dans [Hello Nextflow](../../hello_nextflow/)."
---

# Nextflow pour l'ARNseq

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Traduction assistée par IA - [en savoir plus et suggérer des améliorations](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

**Un cours pratique appliquant Nextflow à un cas d'usage réel de transcriptomique : traitement d'ARNseq en bulk avec Trim Galore, HISAT2 et FastQC.**

Ce cours s'appuie sur la formation débutant [Hello Nextflow](../../hello_nextflow/) et démontre comment utiliser Nextflow dans le contexte spécifique de l'analyse d'ARNseq en bulk.
Vous implémenterez un pipeline de traitement qui élimine les séquences adaptatrices, aligne les lectures sur un génome de référence et effectue un contrôle qualité (QC) à plusieurs étapes.

<!-- additional_information -->

## Aperçu du cours

Ce cours est pratique, avec des exercices orientés objectifs structurés pour introduire progressivement l'information.

Vous commencerez par exécuter les outils de traitement manuellement dans le terminal pour comprendre la méthodologie, puis vous construirez progressivement un pipeline Nextflow qui automatise et met à l'échelle l'analyse.

### Plan de cours

Nous avons divisé ce cours en trois parties qui se concentrent chacune sur des aspects spécifiques de l'application de Nextflow à un cas d'usage d'ARNseq.

| Chapitre du cours                                                               | Résumé                                                                                                                            | Durée estimée |
| ------------------------------------------------------------------------------- | --------------------------------------------------------------------------------------------------------------------------------- | ------------- |
| [Partie 1 : Aperçu de la méthode](./01_method.md)                               | Comprendre la méthodologie de traitement d'ARNseq et exécuter les outils manuellement                                             | 30 min        |
| [Partie 2 : Implémentation mono-échantillon](./02_single-sample.md)             | Construire un pipeline qui traite, aligne et contrôle la qualité d'un seul échantillon, puis le mettre à l'échelle pour plusieurs | 60 min        |
| [Partie 3 : Implémentation multi-échantillons paired-end](./03_multi-sample.md) | Étendre le pipeline pour gérer les données paired-end et agréger les rapports de QC à travers les échantillons                    | 45 min        |

À la fin de ce cours, vous serez capable d'appliquer les concepts et outils fondamentaux de Nextflow à un cas d'usage typique d'ARNseq.

Prêt·e à suivre le cours ?

[Commencer :material-arrow-right:](00_orientation.md){ .md-button .md-button--primary }

<!-- Clearfix for float -->
<div style="content: ''; clear: both; display: table;"></div>
