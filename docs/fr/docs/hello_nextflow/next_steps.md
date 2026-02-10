# Résumé du cours

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Traduction assistée par IA - [en savoir plus et suggérer des améliorations](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

Félicitations pour avoir terminé le cours de formation Hello Nextflow ! 🎉

<div class="video-wrapper">
  <iframe width="560" height="315" src="https://www.youtube.com/embed/wnmUAfDL01E?si=y8lAedhEHWaTV4zd&amp;list=PLPZ8WHdZGxmWKozQuzr27jyMGqp9kElVK&amp;cc_load_policy=1" title="YouTube video player" frameborder="0" allow="accelerometer; autoplay; clipboard-write; encrypted-media; gyroscope; picture-in-picture; web-share" referrerpolicy="strict-origin-when-cross-origin" allowfullscreen></iframe>
</div>

/// caption
:fontawesome-brands-youtube:{ .youtube } Voir la [playlist complète sur la chaîne YouTube Nextflow](https://youtube.com/playlist?list=PLPZ8WHdZGxmWKozQuzr27jyMGqp9kElVK&si=eF7cLR62goy-lc6n).

:green_book: Vous pouvez lire la [transcription de la vidéo](./transcripts/07_next_steps.md) en parallèle de la vidéo.
///

## Votre parcours

Vous avez commencé avec un workflow très basique qui exécutait une commande codée en dur.
Au cours des six parties, vous avez transformé ce workflow basique en un pipeline modulaire à plusieurs étapes qui exploite les fonctionnalités clés de Nextflow, notamment les canaux, les opérateurs, le support intégré des conteneurs et les options de configuration.

### Ce que vous avez construit

- La forme finale du workflow Hello prend en entrée un fichier CSV contenant des salutations textuelles.
- Les quatre étapes sont implémentées comme des processus Nextflow (`sayHello`, `convertToUpper`, `collectGreetings` et `cowpy`) stockés dans des fichiers de modules séparés.
- Les résultats sont publiés dans un répertoire appelé `results/`.
- La sortie finale du pipeline est un fichier texte brut contenant de l'art ASCII d'un personnage prononçant les salutations en majuscules.

<figure class="excalidraw">
--8<-- "docs/en/docs/hello_nextflow/img/hello_pipeline_complete.svg"
</figure>

1. **`sayHello` :** Écrit chaque salutation dans son propre fichier de sortie (_par ex._ « Hello-output.txt »)
2. **`convertToUpper` :** Convertit chaque salutation en majuscules (_par ex._ « HELLO »)
3. **`collectGreetings` :** Collecte toutes les salutations en majuscules dans un seul fichier de lot
4. **`cowpy` :** Génère de l'art ASCII en utilisant l'outil `cowpy`

La configuration du workflow permet de fournir des entrées et des paramètres de manière flexible et reproductible.

### Compétences acquises

Grâce à ce cours pratique, vous avez appris à :

- Décrire et utiliser les composants essentiels de Nextflow suffisants pour construire un workflow simple à plusieurs étapes
- Décrire des concepts de niveau avancé tels que les opérateurs et les fabriques de canaux
- Lancer un workflow Nextflow localement
- Trouver et interpréter les sorties (résultats) et les fichiers de log générés par Nextflow
- Résoudre les problèmes de base

Vous êtes maintenant équipé·e des connaissances fondamentales pour commencer à développer vos propres pipelines avec Nextflow.

## Prochaines étapes pour développer vos compétences

Voici nos 3 principales suggestions pour ce qu'il faut faire ensuite :

- Appliquer Nextflow à un cas d'utilisation d'analyse scientifique avec [Nextflow pour la science](../nf4_science/index.md)
- Démarrer avec nf-core avec [Hello nf-core](../hello_nf-core/index.md)
- Explorer des fonctionnalités Nextflow plus avancées avec les [Quêtes secondaires](../side_quests/index.md)

Enfin, nous vous recommandons de jeter un œil à la [**Plateforme Seqera**](https://seqera.io/), une plateforme basée sur le cloud développée par les créateurs de Nextflow qui facilite encore plus le lancement et la gestion de vos workflows, ainsi que la gestion de vos données et l'exécution d'analyses interactives dans n'importe quel environnement.

## Enquête de satisfaction

Avant de passer à la suite, veuillez prendre une minute pour remplir l'enquête du cours ! Vos commentaires nous aident à améliorer nos supports de formation pour tout le monde.

[Répondre à l'enquête :material-arrow-right:](survey.md){ .md-button .md-button--primary }
