---
title: Hello Nextflow
hide:
  - toc
page_type: index_page
index_type: course
additional_information:
  technical_requirements: true
  learning_objectives:
    - Lancer et gérer l'exécution des flux de travail Nextflow
    - Trouver et interpréter les sorties (résultats) et fichiers de log générés par Nextflow
    - Résoudre les problèmes de base
    - Construire un flux de travail simple à plusieurs étapes à partir des composants essentiels de Nextflow
    - Distinguer les types essentiels de fabriques de canaux et d'opérateurs et les utiliser efficacement dans un flux de travail simple
    - Configurer l'exécution du pipeline pour s'exécuter sur des plateformes de calcul courantes, y compris HPC et cloud
    - Appliquer les bonnes pratiques de reproductibilité, portabilité et réutilisation du code qui rendent les pipelines FAIR, y compris la modularité du code et les conteneurs logiciels
  audience_prerequisites:
    - "**Public :** Ce cours est conçu pour les apprenant·es qui découvrent complètement Nextflow et souhaitent développer leurs propres pipelines."
    - "**Compétences :** Une certaine familiarité avec la ligne de commande, les concepts de base des scripts et les formats de fichiers courants est supposée."
    - "**Domaine :** Les exercices sont tous indépendants du domaine, aucune connaissance scientifique préalable n'est donc requise."
  videos_playlist: https://www.youtube.com/playlist?list=PLPZ8WHdZGxmXiHf8B26oB_fTfoKQdhlik
---

# Hello Nextflow

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Traduction assistée par IA - [en savoir plus et suggérer des améliorations](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

**Hello Nextflow est une introduction pratique à la construction de flux de travail d'analyse de données reproductibles et évolutifs.**

En travaillant sur des exemples pratiques et des exercices guidés, vous apprendrez les fondamentaux du développement de pipelines avec Nextflow, notamment comment définir des processus, les connecter dans des pipelines, gérer les fichiers et les dépendances logicielles, paralléliser l'exécution sans effort et exécuter des flux de travail dans différents environnements de calcul.

Vous repartirez avec les compétences et la confiance nécessaires pour commencer à développer et exécuter vos propres flux de travail avec Nextflow.

<!-- additional_information -->

## Aperçu du cours

Ce cours est conçu pour être pratique, avec des exercices orientés vers des objectifs structurés pour introduire l'information progressivement.

Vous développerez un pipeline Nextflow simple qui prend des entrées de texte, exécute quelques étapes de transformation et produit un fichier texte unique contenant une image ASCII d'un personnage prononçant le texte transformé.

### Plan de cours

Afin de ne pas vous submerger de concepts et de code, nous avons divisé cela en six parties qui se concentreront chacune sur des aspects spécifiques du développement de pipelines avec Nextflow.

| Chapitre du cours                                       | Résumé                                                                                                                   | Durée estimée |
| ------------------------------------------------------- | ------------------------------------------------------------------------------------------------------------------------ | ------------- |
| [Partie 1 : Hello World](./01_hello_world.md)           | Composants de base et principes impliqués dans l'assemblage et l'exécution d'un flux de travail Nextflow                 | 30 min        |
| [Partie 2 : Hello Channels](./02_hello_channels.md)     | Utilisation des canaux et opérateurs pour traiter les entrées et paralléliser l'exécution sans effort                    | 45 min        |
| [Partie 3 : Hello Workflow](./03_hello_workflow.md)     | Utilisation des canaux pour enchaîner plusieurs étapes et gérer le transfert de données entre les étapes                 | 60 min        |
| [Partie 4 : Hello Modules](./04_hello_modules.md)       | Application des principes de modularité du code pour augmenter la réutilisabilité et diminuer la charge de maintenance   | 20 min        |
| [Partie 5 : Hello Containers](./05_hello_containers.md) | Utilisation des conteneurs comme mécanisme de gestion des dépendances logicielles et augmentation de la reproductibilité | 60 min        |
| [Partie 6 : Hello Config](./06_hello_config.md)         | Personnalisation du comportement du pipeline et optimisation de l'utilisation dans différents environnements de calcul   | 60 min        |

À la fin de ce cours, vous serez bien préparé·e pour aborder les prochaines étapes de votre parcours pour développer des flux de travail reproductibles pour vos besoins de calcul scientifique.

Prêt·e à suivre le cours ?

[Commencer :material-arrow-right:](00_orientation.md){ .md-button .md-button--primary }

<!-- Clearfix for float -->
<div style="content: ''; clear: both; display: table;"></div>
