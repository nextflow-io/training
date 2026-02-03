# Nextflow pour la Génomique

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Traduction assistée par IA - [en savoir plus et suggérer des améliorations](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

Cette formation est destinée aux chercheur·euses en génomique et domaines connexes qui souhaitent développer ou personnaliser des pipelines d'analyse de données.
Elle s'appuie sur la formation pour débutants [Hello Nextflow](../../hello_nextflow/) et démontre comment utiliser Nextflow dans le contexte spécifique du domaine de la génomique.

Plus précisément, ce cours démontre comment implémenter un pipeline simple d'appel de variants avec [GATK](https://gatk.broadinstitute.org/) (Genome Analysis Toolkit), un ensemble logiciel largement utilisé pour l'analyse de données de séquençage à haut débit.

Commençons ! Cliquez sur le bouton « Open in GitHub Codespaces » ci-dessous pour lancer l'environnement de formation (de préférence dans un onglet séparé), puis poursuivez la lecture pendant le chargement.

[![Open in GitHub Codespaces](https://github.com/codespaces/badge.svg)](https://codespaces.new/nextflow-io/training?quickstart=1&ref=master)

## Objectifs pédagogiques

En suivant ce cours, vous apprendrez comment appliquer les concepts et outils fondamentaux de Nextflow à un cas d'usage typique en génomique.

À la fin de cet atelier, vous serez capable de :

- Écrire un workflow linéaire pour appliquer l'appel de variants à un seul échantillon
- Gérer correctement les fichiers accessoires tels que les fichiers d'index et les ressources du génome de référence
- Exploiter le paradigme de flux de données de Nextflow pour paralléliser l'appel de variants par échantillon
- Implémenter l'appel de variants multi-échantillons en utilisant les opérateurs de canaux appropriés
- Implémenter des tests par étape et de bout en bout qui gèrent correctement les particularités spécifiques à la génomique

<!-- TODO for future expansion: add metadata/samplesheet handling -->

## Prérequis

Le cours suppose une familiarité minimale avec les éléments suivants :

- Outils et formats de fichiers couramment utilisés dans ce domaine scientifique
- Expérience avec la ligne de commande
- Concepts et outils fondamentaux de Nextflow couverts dans la formation pour débutants [Hello Nextflow](../../hello_nextflow/)

Pour les exigences techniques et la configuration de l'environnement, consultez le mini-cours [Configuration de l'Environnement](../../envsetup/).
