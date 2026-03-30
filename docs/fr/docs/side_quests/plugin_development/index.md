---
title: Développement de plugins
hide:
  - toc
---

# Développement de plugins

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Traduction assistée par IA - [en savoir plus et suggérer des améliorations](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

Le système de plugins de Nextflow vous permet d'étendre le langage avec des fonctions personnalisées, des hooks de surveillance, des backends d'exécution, et bien plus encore.
Les plugins permettent à la communauté d'ajouter des fonctionnalités à Nextflow sans modifier son cœur, ce qui les rend idéaux pour partager des fonctionnalités réutilisables entre pipelines.

Au cours de cette formation, vous apprendrez à utiliser des plugins existants et, si vous le souhaitez, à créer les vôtres.

## Public & prérequis

La partie 1 couvre l'utilisation de plugins existants et s'adresse à tous les utilisateur·trices de Nextflow.
Les parties 2 à 6 couvrent la création de vos propres plugins et impliquent du code Groovy ainsi que des outils de build.
Aucune expérience préalable en Java ou en Groovy n'est requise.

**Prérequis**

- Un compte GitHub OU une installation locale telle que décrite [ici](../../envsetup/02_local).
- Avoir suivi le cours [Hello Nextflow](../../hello_nextflow/index.md) ou équivalent.
- Java 21 ou version ultérieure (inclus dans l'environnement de formation ; uniquement nécessaire pour les parties 2 à 6).

**Répertoire de travail :** `side-quests/plugin_development`

## Objectifs pédagogiques

À la fin de cette formation, vous serez en mesure de :

**Utiliser des plugins (Partie 1) :**

- Installer et configurer des plugins existants dans vos workflows
- Importer et utiliser des fonctions de plugins

**Développer des plugins (Parties 2 à 6) :**

- Créer un nouveau projet de plugin à l'aide du générateur de projets intégré de Nextflow
- Implémenter des fonctions personnalisées appelables depuis des workflows
- Compiler, tester et installer votre plugin localement
- Surveiller les événements du workflow (par exemple, la fin d'une tâche, le démarrage/l'arrêt d'un pipeline) pour des journaux ou des notifications personnalisés
- Ajouter des options de configuration pour rendre les plugins personnalisables
- Distribuer votre plugin

## Plan de cours

#### Partie 1 : Les bases des plugins

Utiliser des plugins existants dans un workflow Nextflow et configurer leur comportement.

#### Partie 2 : Créer un projet de plugin

Générer un nouveau projet de plugin et examiner sa structure.

#### Partie 3 : Fonctions personnalisées

Implémenter des fonctions personnalisées, compiler votre plugin et l'exécuter dans un workflow.

#### Partie 4 : Tests

Écrire et exécuter des tests unitaires à l'aide du framework Spock.

#### Partie 5 : Surveillance du workflow

Réagir à des événements tels que la fin d'une tâche pour construire un compteur de tâches.

#### Partie 6 : Configuration & Distribution

Lire des paramètres depuis `nextflow.config` pour rendre votre plugin personnalisable, puis apprendre à le partager.

Prêt·e à suivre le cours ?

[Commencer l'apprentissage :material-arrow-right:](01_plugin_basics.md){ .md-button .md-button--primary }
