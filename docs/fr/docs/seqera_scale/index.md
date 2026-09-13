---
title: Mise à l'échelle avec Seqera
hide:
  - toc
page_type: index_page
index_type: course
additional_information:
  technical_requirements: true
  learning_objectives:
    - S'inscrire sur Seqera Platform et explorer la vitrine communautaire
    - Ajouter un pipeline à un espace de travail et le lancer depuis l'interface web
    - S'authentifier et lancer des pipelines depuis la ligne de commande avec le CLI `tw`
    - Enregistrer un pipeline hébergé sur GitHub et le lancer des deux façons
  audience_prerequisites:
    - "**Public :** Ce cours est conçu pour les apprenant·es qui souhaitent exécuter des pipelines Nextflow à grande échelle avec Seqera Platform."
    - "**Compétences :** Une familiarité avec l'exécution de pipelines nf-core depuis la ligne de commande est supposée."
    - "**Cours :** Avoir complété [Nextflow Run](../nextflow_run/index.md) et [Use nf-core](../nfcore_use/index.md), ou être à l'aise avec l'exécution de pipelines locaux et `nf-core/rnaseq`."
---

# Mise à l'échelle avec Seqera

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Traduction assistée par IA - [en savoir plus et suggérer des améliorations](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

**Mise à l'échelle avec Seqera est une introduction pratique au lancement et à la surveillance de pipelines Nextflow avec Seqera Platform.**

À travers des exemples concrets, vous configurerez l'accès à Seqera Platform, lancerez un pipeline à l'échelle de production depuis l'interface web et depuis la ligne de commande, et ajouterez un nouveau pipeline à votre espace de travail.

Vous repartirez avec les compétences et la confiance nécessaires pour exécuter et surveiller vos propres pipelines sur Seqera Platform.

<!-- additional_information -->

## Aperçu du cours

Ce cours est pratique et s'appuie sur les pipelines que vous avez déjà exécutés dans [Use nf-core](../nfcore_use/index.md).

Vous commencerez par vous inscrire sur Seqera Platform et lancer `nf-core/rnaseq`, un pipeline à l'échelle de production, depuis l'interface web.
Vous passerez ensuite à l'outil en ligne de commande `tw` pour faire de même depuis un terminal, et enfin vous enregistrerez un nouveau pipeline, `nf-core/demo`, et le lancerez des deux façons.

### Plan du cours

| Chapitre du cours                                                          | Résumé                                                                                                                    | Durée estimée |
| -------------------------------------------------------------------------- | ------------------------------------------------------------------------------------------------------------------------- | ------------- |
| [Partie 1 : Lancer des pipelines depuis l'interface web](./01_run_with_seqera.md) | Configurer l'accès à Seqera Platform et lancer un pipeline à l'échelle de production depuis l'interface web | 20 min        |
| [Partie 2 : Lancer des pipelines depuis la ligne de commande](./02_launch_from_cli.md) | Authentifier le CLI `tw`, lancer un pipeline enregistré et en enregistrer un nouveau depuis le CLI | 25 min        |

À la fin de ce cours, vous serez à l'aise pour lancer et surveiller des pipelines Nextflow sur Seqera Platform, que vous préfériez travailler depuis l'interface web ou depuis la ligne de commande.

Prêt·e à suivre le cours ?

[Commencer l'apprentissage :material-arrow-right:](00_orientation.md){ .md-button .md-button--primary }
