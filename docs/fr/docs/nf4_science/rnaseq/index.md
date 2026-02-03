---
title: Nextflow pour l'ARNseq
hide:
  - toc
---

# Nextflow pour l'ARNseq

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Traduction assistée par IA - [en savoir plus et suggérer des améliorations](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

Ce cours de formation est destiné aux chercheur·euses en transcriptomique et dans les domaines connexes qui souhaitent développer ou personnaliser des pipelines d'analyse de données.
Il s'appuie sur la formation débutant [Hello Nextflow](../../hello_nextflow/) et démontre comment utiliser Nextflow dans le contexte spécifique de l'analyse d'ARNseq en bulk.

Plus précisément, ce cours démontre comment implémenter un pipeline simple de traitement d'ARNseq en bulk pour éliminer les séquences adaptatrices, aligner les lectures sur un génome de référence et effectuer un contrôle qualité (QC) à plusieurs étapes.

Commençons ! Cliquez sur le bouton « Open in GitHub Codespaces » ci-dessous pour lancer l'environnement de formation (de préférence dans un onglet séparé), puis poursuivez la lecture pendant le chargement.

[![Open in GitHub Codespaces](https://github.com/codespaces/badge.svg)](https://codespaces.new/nextflow-io/training?quickstart=1&ref=master)

## Objectifs d'apprentissage

En suivant ce cours, vous apprendrez à appliquer les concepts et outils fondamentaux de Nextflow à un cas d'usage typique d'ARNseq.

À la fin de cet atelier, vous serez capable de :

- Écrire un workflow linéaire pour appliquer des méthodes de base de traitement et de QC d'ARNseq
- Gérer de manière appropriée les fichiers spécifiques au domaine tels que les FASTQ et les ressources de génome de référence
- Gérer les données de séquençage single-end et paired-end
- Exploiter le paradigme de flux de données de Nextflow pour paralléliser le traitement d'ARNseq par échantillon
- Agréger les rapports de QC à travers plusieurs étapes et échantillons en utilisant les opérateurs de canal pertinents

<!-- TODO
- Configure pipeline execution and manage and optimize resource allocations
- Implement per-step and end-to-end pipeline tests that handle RNAseq-specific idiosyncrasies appropriately
-->
<!-- TODO for future expansion: add metadata/samplesheet handling -->

## Prérequis

Ce cours suppose une familiarité minimale avec les éléments suivants :

- Outils et formats de fichiers couramment utilisés dans ce domaine scientifique
- Expérience avec la ligne de commande
- Concepts et outils fondamentaux de Nextflow couverts dans la formation débutant [Hello Nextflow](../../hello_nextflow/).

Pour les exigences techniques et la configuration de l'environnement, consultez le mini-cours [Configuration de l'environnement](../../envsetup/).
