# Résumé du cours

Félicitations pour avoir terminé le cours de formation Hello Nextflow ! 🎉

<div class="video-wrapper">
  <iframe width="560" height="315" src="https://www.youtube.com/embed/wnmUAfDL01E?si=y8lAedhEHWaTV4zd&amp;list=PLPZ8WHdZGxmWKozQuzr27jyMGqp9kElVK&amp;cc_load_policy=1" title="YouTube video player" frameborder="0" allow="accelerometer; autoplay; clipboard-write; encrypted-media; gyroscope; picture-in-picture; web-share" referrerpolicy="strict-origin-when-cross-origin" allowfullscreen></iframe>
</div>

/// caption
:fontawesome-brands-youtube:{ .youtube } Voir la [playlist complète sur la chaîne YouTube Nextflow](https://youtube.com/playlist?list=PLPZ8WHdZGxmWKozQuzr27jyMGqp9kElVK&si=eF7cLR62goy-lc6n).

:green_book: Vous pouvez lire la [transcription vidéo](./transcripts/07_next_steps.md) en parallèle de la vidéo.
///

## Votre parcours

Vous avez commencé avec un workflow très basique qui exécutait une commande codée en dur.
Au cours de six parties, vous avez transformé ce workflow basique en un pipeline modulaire multi-étapes qui met en œuvre des fonctionnalités clés de Nextflow, notamment les canaux, les opérateurs, la prise en charge intégrée des conteneurs et les options de configuration.

### Ce que vous avez construit

- La forme finale du workflow Hello prend en entrée un fichier CSV contenant des messages de bienvenue textuels.
- Les quatre étapes sont implémentées sous forme de processus Nextflow (`sayHello`, `convertToUpper`, `collectGreetings` et `cowpy`) stockés dans des fichiers de modules séparés.
- Les résultats sont publiés dans un répertoire appelé `results/`.
- La sortie finale du pipeline est un fichier texte brut contenant de l'art ASCII d'un personnage disant les messages de bienvenue en majuscules.

<figure class="excalidraw">
--8<-- "docs/en/docs/hello_nextflow/img/hello_pipeline_complete.svg"
</figure>

1. **`sayHello` :** Écrit chaque message de bienvenue dans son propre fichier de sortie (_par exemple_ "Hello-output.txt")
2. **`convertToUpper` :** Convertit chaque message de bienvenue en majuscules (_par exemple_ "HELLO")
3. **`collectGreetings` :** Collecte tous les messages de bienvenue en majuscules dans un seul fichier batch
4. **`cowpy` :** Génère de l'art ASCII en utilisant l'outil `cowpy`

La configuration du workflow permet de fournir des entrées et des paramètres de manière flexible et reproductible.

### Compétences acquises

Grâce à ce cours pratique, vous avez appris à :

- Décrire et utiliser les composants principaux de Nextflow suffisants pour construire un workflow simple multi-étapes
- Décrire des concepts de niveau suivant tels que les opérateurs et les fabriques de canaux
- Lancer un workflow Nextflow localement
- Trouver et interpréter les sorties (résultats) et les fichiers journaux générés par Nextflow
- Résoudre des problèmes de base

Vous êtes maintenant équipé·e des connaissances fondamentales pour commencer à développer vos propres pipelines dans Nextflow.

## Prochaines étapes pour développer vos compétences

Voici nos 3 principales suggestions pour la suite :

- Appliquer Nextflow à un cas d'usage d'analyse scientifique avec [Nextflow for Science](../nf4_science/index.md)
- Démarrer avec nf-core avec [Hello nf-core](../hello_nf-core/index.md)
- Explorer des fonctionnalités plus avancées de Nextflow avec les [Quêtes secondaires](../side_quests/index.md)

Enfin, nous vous recommandons de jeter un œil à [**Seqera Platform**](https://seqera.io/), une plateforme cloud développée par les créateurs de Nextflow qui facilite encore davantage le lancement et la gestion de vos workflows, ainsi que la gestion de vos données et l'exécution d'analyses de manière interactive dans n'importe quel environnement.

## Questionnaire de satisfaction

Avant de continuer, veuillez prendre une minute pour compléter le questionnaire du cours ! Vos retours nous aident à améliorer nos supports de formation pour tout le monde.

[Répondre au questionnaire :material-arrow-right:](survey.md){ .md-button .md-button--primary }
