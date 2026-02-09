---
title: Nextflow pour la Génomique
hide:
  - toc
page_type: index_page
index_type: course
additional_information:
  technical_requirements: true
  learning_objectives:
    - Écrire un workflow linéaire pour appliquer l'appel de variants à un seul échantillon
    - Gérer correctement les fichiers accessoires tels que les fichiers d'index et les ressources du génome de référence
    - Exploiter le paradigme de flux de données de Nextflow pour paralléliser l'appel de variants par échantillon
    - Implémenter l'appel de variants multi-échantillons en utilisant les opérateurs de canaux appropriés
  audience_prerequisites:
    - "**Public :** Cette formation est destinée aux chercheur·euses en génomique et domaines connexes qui souhaitent développer ou personnaliser des pipelines d'analyse de données."
    - "**Compétences :** Une certaine familiarité avec la ligne de commande, les concepts de base du scripting et les formats de fichiers courants en génomique est supposée."
    - "**Prérequis :** Concepts et outils fondamentaux de Nextflow couverts dans [Hello Nextflow](../../hello_nextflow/)."
---

# Nextflow pour la Génomique

**Une formation pratique appliquant Nextflow à un cas d'usage réel en génomique : l'appel de variants avec GATK.**

Cette formation s'appuie sur la formation pour débutants [Hello Nextflow](../../hello_nextflow/) et démontre comment utiliser Nextflow dans le contexte spécifique du domaine de la génomique.
Vous allez implémenter un pipeline d'appel de variants avec [GATK](https://gatk.broadinstitute.org/) (Genome Analysis Toolkit), un ensemble logiciel largement utilisé pour l'analyse de données de séquençage à haut débit.

<!-- additional_information -->

## Aperçu de la formation

Cette formation est pratique, avec des exercices orientés objectifs structurés pour introduire l'information progressivement.

Vous commencerez par exécuter les outils d'appel de variants manuellement dans le terminal pour comprendre la méthodologie, puis vous construirez progressivement un pipeline Nextflow qui automatise et met à l'échelle l'analyse.

### Plan de cours

Nous avons divisé cette formation en trois parties qui se concentrent chacune sur des aspects spécifiques de l'application de Nextflow à un cas d'usage en génomique.

| Chapitre de la formation                                                           | Résumé                                                                                                                           | Durée estimée |
| ---------------------------------------------------------------------------------- | -------------------------------------------------------------------------------------------------------------------------------- | ------------- |
| [Partie 1 : Aperçu de la méthode](./01_method.md)                                  | Comprendre la méthodologie d'appel de variants et exécuter les outils manuellement                                               | 30 min        |
| [Partie 2 : Appel de variants par échantillon](./02_per_sample_variant_calling.md) | Construire un pipeline qui indexe les fichiers BAM et appelle les variants, puis mettre à l'échelle sur plusieurs échantillons   | 60 min        |
| [Partie 3 : Appel conjoint sur une cohorte](./03_joint_calling.md)                 | Ajouter le génotypage conjoint multi-échantillons en utilisant les opérateurs de canaux pour agréger les sorties par échantillon | 45 min        |

À la fin de cette formation, vous serez capable d'appliquer les concepts et outils fondamentaux de Nextflow à un cas d'usage typique en génomique.

Prêt·e à suivre la formation ?

[Commencer :material-arrow-right:](00_orientation.md){ .md-button .md-button--primary }

<!-- Clearfix for float -->
<div style="content: ''; clear: both; display: table;"></div>
