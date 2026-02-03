# Nextflow run pour l'imagerie

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Traduction assistée par IA - [en savoir plus et suggérer des améliorations](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

Ce cours de formation s'adresse aux chercheur·euses en imagerie et biologie spatiale qui souhaitent exécuter et personnaliser des pipelines d'analyse de données.
Il enseigne les concepts fondamentaux de Nextflow liés à l'exécution, à l'organisation et à la configuration de workflows en utilisant [nf-core/molkart](https://nf-co.re/molkart), un pipeline pour le traitement de données de transcriptomique spatiale Molecular Cartography.
Les compétences que vous apprendrez ici sont transférables à n'importe quel pipeline Nextflow ou nf-core.

Commençons ! Cliquez sur le bouton « Open in GitHub Codespaces » ci-dessous pour lancer l'environnement de formation (de préférence dans un onglet séparé), puis continuez votre lecture pendant qu'il se charge.

[![Open in GitHub Codespaces](https://github.com/codespaces/badge.svg)](https://codespaces.new/nextflow-io/training?quickstart=1&ref=master)

## Objectifs d'apprentissage

En suivant ce cours, vous apprendrez à appliquer les concepts et outils fondamentaux de Nextflow à l'exécution de pipelines d'analyse d'imagerie.

À la fin de cet atelier, vous serez capable de :

- Lancer un workflow Nextflow localement et surveiller son exécution
- Trouver et interpréter les sorties (résultats) et fichiers de log générés par Nextflow
- Exécuter un pipeline nf-core avec des données de test et des entrées personnalisées
- Configurer l'exécution du pipeline en utilisant des profils et des fichiers de paramètres
- Gérer les entrées en utilisant des samplesheets et des paramètres en ligne de commande

## Public et prérequis

Ce cours suppose une familiarité minimale avec les éléments suivants :

- Expérience avec la ligne de commande
- Familiarité de base avec les formats de fichiers d'imagerie (images TIFF, données tabulaires)

Pour les exigences techniques et la configuration de l'environnement, consultez le mini-cours [Configuration de l'environnement](../../envsetup/).
