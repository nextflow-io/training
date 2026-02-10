---
title: Hello nf-core
hide:
  - toc
page_type: index_page
index_type: course
additional_information:
  technical_requirements: true
  learning_objectives:
    - Récupérer, lancer et gérer l'exécution de pipelines nf-core
    - Décrire la structure du code et l'organisation de projet des pipelines nf-core
    - Créer un pipeline compatible nf-core de base à partir d'un modèle
    - Mettre à niveau un workflow Nextflow simple pour respecter les standards nf-core
    - Ajouter des modules nf-core à un pipeline compatible nf-core
    - Contribuer vos propres modules à nf-core
    - Valider les entrées et les paramètres en utilisant les outils nf-core
  audience_prerequisites:
    - "**Public :** Ce cours est conçu pour les apprenant·es qui connaissent déjà les bases de Nextflow et souhaitent apprendre à utiliser les ressources et les bonnes pratiques nf-core."
    - "**Compétences :** Une familiarité avec la ligne de commande, les concepts de base de script et les formats de fichiers courants est supposée."
    - "**Cours :** Vous devez avoir terminé le cours [Hello Nextflow](../hello_nextflow/index.md) ou équivalent."
    - "**Domaine :** Les exercices sont tous indépendants du domaine scientifique, donc aucune connaissance scientifique préalable n'est requise."
---

# Hello nf-core

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Traduction assistée par IA - [en savoir plus et suggérer des améliorations](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

**Hello nf-core est une introduction pratique à l'utilisation des ressources et des bonnes pratiques nf-core.**

![logo nf-core](./img/nf-core-logo.png#only-light)
![logo nf-core](./img/nf-core-logo-darkbg.png#only-dark)

En travaillant sur des exemples pratiques et des exercices guidés, vous apprendrez à utiliser et développer des modules et pipelines compatibles nf-core, et à utiliser efficacement les outils nf-core.

Vous acquerrez les compétences et la confiance nécessaires pour commencer à développer des pipelines selon les bonnes pratiques nf-core.

<!-- additional_information -->

## Aperçu du cours

Ce cours est conçu pour être pratique, avec des exercices orientés objectifs structurés pour introduire l'information progressivement.

Vous serez initié·e à [**nf-core**](https://nf-co.re/), un effort communautaire pour développer et maintenir un ensemble organisé de pipelines scientifiques construits avec Nextflow, ainsi que les outils et directives pertinents qui favorisent le développement ouvert, les tests et l'évaluation par les pairs ([Nat Biotechnol 38, 276–278 (2020)](https://www.nature.com/articles/s41587-020-0439-x), [Genome Biol 26, 228 (2025)](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-025-03673-9)).

Les pipelines développés par la communauté nf-core sont conçus pour être modulaires, évolutifs et portables, permettant aux chercheur·euses de les adapter et de les exécuter facilement avec leurs propres données et ressources de calcul.
Les directives de bonnes pratiques imposées par le projet garantissent en outre que les pipelines sont robustes, bien documentés et validés sur des ensembles de données réels.
Cela contribue à augmenter la fiabilité et la reproductibilité des analyses scientifiques et permet finalement aux chercheur·euses d'accélérer leurs découvertes scientifiques.

Nous ne couvrirons pas tout ce qu'il y a à savoir sur les pipelines nf-core dans ce cours, car nf-core englobe de nombreuses fonctionnalités et conventions développées par la communauté au fil des années.
Au lieu de cela, nous nous concentrerons sur les concepts essentiels qui vous aideront à démarrer et à comprendre comment fonctionne nf-core.

### Plan de cours

Nous avons divisé cela en cinq parties qui se concentreront chacune sur des aspects spécifiques de l'utilisation des ressources nf-core.

| Chapitre du cours                                                        | Résumé                                                                                                                                                                 | Durée estimée |
| ------------------------------------------------------------------------ | ---------------------------------------------------------------------------------------------------------------------------------------------------------------------- | ------------- |
| [Partie 1 : Exécuter un pipeline de démonstration](./01_run_demo.md)     | Exécuter un pipeline nf-core existant et examiner sa structure de code pour avoir une idée de ce qui distingue ces pipelines des workflows Nextflow de base            | 30 min        |
| [Partie 2 : Réécrire Hello pour nf-core](./02_rewrite_hello.md)          | Adapter un workflow existant à la structure du modèle nf-core, en partant du workflow simple produit dans le cours [Hello Nextflow](../hello_nextflow/index.md)        | 60 min        |
| [Partie 3 : Utiliser un module nf-core](./03_use_module.md)              | Explorer la bibliothèque de modules de la communauté et apprendre à intégrer des modules pré-construits et testés qui encapsulent des outils bioinformatiques courants | 30 min        |
| [Partie 4 : Créer un module nf-core](./04_make_module.md)                | Créer votre propre module de style nf-core en utilisant la structure spécifique, les conventions de nommage et les exigences de métadonnées établies par nf-core       | 30 min        |
| [Partie 5 : Ajouter la validation des entrées](./05_input_validation.md) | Implémenter la validation des entrées pour les paramètres de ligne de commande et les fichiers de données d'entrée en utilisant nf-schema                              | 30 min        |

À la fin de ce cours, vous serez capable de tirer parti de l'énorme richesse de ressources offertes par le projet nf-core.

Prêt·e à commencer le cours ?

[Commencer l'apprentissage :material-arrow-right:](00_orientation.md){ .md-button .md-button--primary }
