# Résumé de la formation

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Traduction assistée par IA - [en savoir plus et suggérer des améliorations](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

Félicitations pour avoir terminé la formation Nextflow Run !

<!-- placeholder for video -->

## Votre parcours

Vous avez commencé avec un workflow très basique, et avez appris à l'exécuter, trouver les sorties, et gérer son exécution.
Ensuite, vous avez travaillé à travers des versions de plus en plus complexes de ce workflow et avez appris à reconnaître les concepts et mécanismes essentiels qui alimentent les pipelines Nextflow, y compris les channels et les opérateurs, la modularisation du code, et les conteneurs.
Enfin, vous avez appris à personnaliser la configuration d'un pipeline pour s'adapter à vos préférences et à votre infrastructure de calcul.

### Ce que vous avez appris

Vous êtes maintenant capable de gérer l'exécution du pipeline Hello, de décrire comment il est structuré, et d'identifier les principales parties du code impliquées.

- La forme finale du workflow Hello prend en entrée un fichier CSV contenant des salutations textuelles.
- Les quatre étapes sont implémentées comme des processes Nextflow (`sayHello`, `convertToUpper`, `collectGreetings`, et `cowpy`) stockés dans des fichiers de module séparés.
- Les résultats sont publiés dans un répertoire appelé `results/`.
- La sortie finale du pipeline est un fichier texte brut contenant de l'art ASCII d'un personnage disant les salutations en majuscules.

<figure class="excalidraw">
--8<-- "docs/en/docs/hello_nextflow/img/hello_pipeline_complete.svg"
</figure>

1. **`sayHello` :** Écrit chaque salutation dans son propre fichier de sortie (par exemple « Hello-output.txt »)
2. **`convertToUpper` :** Convertit chaque salutation en majuscules (par exemple « HELLO »)
3. **`collectGreetings` :** Collecte toutes les salutations en majuscules dans un seul fichier de lot
4. **`cowpy` :** Génère de l'art ASCII en utilisant l'outil `cowpy`

La configuration du workflow prend en charge la fourniture d'entrées et de paramètres de manière flexible et reproductible.

### Compétences acquises

À travers cette formation pratique, vous avez appris à :

- Lancer un workflow Nextflow localement
- Trouver et interpréter les sorties (résultats) et les fichiers de log générés par Nextflow
- Reconnaître les composants principaux de Nextflow qui constituent un workflow multi-étapes simple
- Décrire des concepts avancés tels que les opérateurs et les fabriques de channel
- Configurer des pipelines pour différents environnements de calcul

Vous êtes maintenant équipé·e des connaissances fondamentales pour commencer à intégrer des pipelines Nextflow existants dans votre propre travail.

## Prochaines étapes pour développer vos compétences

Voici nos principales suggestions pour la suite :

- Ne vous contentez pas d'exécuter Nextflow, écrivez-le ! Devenez un·e développeur·se Nextflow avec [Hello Nextflow](../hello_nextflow/index.md)
- Appliquez Nextflow à un cas d'utilisation d'analyse scientifique avec [Nextflow for Science](../nf4_science/index.md)
- Commencez avec nf-core avec [Hello nf-core](../hello_nf-core/index.md)
- Apprenez les techniques de dépannage avec le [Debugging Side Quest](../side_quests/debugging.md)

Enfin, nous vous recommandons de jeter un œil à [**Seqera Platform**](https://seqera.io/), une plateforme basée sur le cloud développée par les créateurs de Nextflow qui facilite encore plus le lancement et la gestion de vos workflows, ainsi que la gestion de vos données et l'exécution d'analyses de manière interactive dans n'importe quel environnement.

## Obtenir de l'aide

Pour les ressources d'aide et le support communautaire, consultez la [page d'aide](../help.md).

## Enquête de satisfaction

Avant de passer à autre chose, veuillez prendre une minute pour compléter l'enquête de la formation ! Vos retours nous aident à améliorer nos supports de formation pour tout le monde.

[Répondre à l'enquête :material-arrow-right:](survey.md){ .md-button .md-button--primary }
