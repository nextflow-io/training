# Nextflow pour la génomique

**Un cours pratique appliquant Nextflow à un cas d'usage réel en génomique : l'appel de variants avec GATK.**

Ce cours s'appuie sur la formation pour débutant·es [Hello Nextflow](../../hello_nextflow/) et démontre comment utiliser Nextflow dans le contexte spécifique du domaine de la génomique.
Vous implémenterez un pipeline d'appel de variants avec [GATK](https://gatk.broadinstitute.org/) (Genome Analysis Toolkit), un logiciel largement utilisé pour l'analyse de données de séquençage à haut débit.

<!-- additional_information -->

## Aperçu du cours

Ce cours est pratique, avec des exercices orientés vers des objectifs et structurés pour introduire l'information progressivement.

Vous commencerez par exécuter les outils d'appel de variants manuellement dans le terminal pour comprendre la méthodologie, puis vous construirez progressivement un pipeline Nextflow qui automatise et met à l'échelle l'analyse.

### Plan de cours

Nous avons divisé ce cours en trois parties qui se concentrent chacune sur des aspects spécifiques de l'application de Nextflow à un cas d'usage en génomique.

| Chapitre du cours                                                        | Résumé                                                                                                                                | Durée estimée |
| ------------------------------------------------------------------------ | ------------------------------------------------------------------------------------------------------------------------------------- | ------------- |
| [Partie 1 : Aperçu de la méthode](./01_method.md)                        | Comprendre la méthodologie d'appel de variants et exécuter les outils manuellement                                                   | 30 min        |
| [Partie 2 : Appel de variants par échantillon](./02_per_sample_variant_calling.md) | Construire un pipeline qui indexe les fichiers BAM et appelle les variants, puis le mettre à l'échelle pour plusieurs échantillons | 60 min        |
| [Partie 3 : Appel conjoint sur une cohorte](./03_joint_calling.md)      | Ajouter le génotypage conjoint multi-échantillons en utilisant des opérateurs de canaux pour agréger les sorties par échantillon    | 45 min        |

À la fin de ce cours, vous serez capable d'appliquer les concepts et outils fondamentaux de Nextflow à un cas d'usage typique en génomique.

Prêt·e à suivre le cours ?

[Commencer :material-arrow-right:](00_orientation.md){ .md-button .md-button--primary }

<!-- Clearfix for float -->
<div style="content: ''; clear: both; display: table;"></div>
