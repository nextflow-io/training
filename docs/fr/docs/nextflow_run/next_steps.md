# Résumé du cours

Félicitations pour avoir terminé le cours de formation Nextflow Run ! 🎉

<!-- placeholder for video -->

## Votre parcours

Vous avez commencé avec un workflow très basique, et avez appris à l'exécuter, à trouver les sorties et à gérer son exécution.
Ensuite, vous avez progressé à travers des versions de plus en plus complexes de ce workflow et avez appris à reconnaître les concepts et mécanismes essentiels qui alimentent les pipelines Nextflow, notamment les canaux et opérateurs, la modularisation du code et les conteneurs.
Enfin, vous avez appris à personnaliser la configuration d'un pipeline pour l'adapter à vos préférences et à votre infrastructure informatique.

### Ce que vous avez appris

Vous êtes maintenant capable de gérer l'exécution du pipeline Hello, de décrire comment il est structuré et d'identifier les principaux éléments de code impliqués.

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

- Lancer un workflow Nextflow localement
- Trouver et interpréter les sorties (résultats) et les fichiers journaux générés par Nextflow
- Reconnaître les composants Nextflow essentiels qui constituent un workflow simple en plusieurs étapes
- Décrire des concepts avancés tels que les opérateurs et les fabriques de canaux
- Configurer des pipelines pour différents environnements informatiques

Vous êtes maintenant équipé·e des connaissances fondamentales pour commencer à intégrer des pipelines Nextflow existants dans votre propre travail.

## Prochaines étapes pour développer vos compétences

Voici nos principales suggestions pour la suite :

- Ne vous contentez pas d'exécuter Nextflow, écrivez-le ! Devenez un·e développeur·se Nextflow avec [Hello Nextflow](../hello_nextflow/index.md)
- Appliquez Nextflow à un cas d'usage d'analyse scientifique avec [Nextflow for Science](../nf4_science/index.md)
- Démarrez avec nf-core grâce à [Hello nf-core](../hello_nf-core/index.md)
- Apprenez des techniques de dépannage avec la [Quête secondaire Débogage](../side_quests/debugging.md)

Enfin, nous vous recommandons de découvrir [**Seqera Platform**](https://seqera.io/), une plateforme cloud développée par les créateurs de Nextflow qui facilite encore davantage le lancement et la gestion de vos workflows, ainsi que la gestion de vos données et l'exécution d'analyses de manière interactive dans n'importe quel environnement.

## Obtenir de l'aide

Pour les ressources d'aide et le support communautaire, consultez la [page d'aide](../help.md).

## Questionnaire de satisfaction

Avant de continuer, veuillez prendre une minute pour compléter le questionnaire du cours ! Vos retours nous aident à améliorer nos supports de formation pour tout le monde.

[Répondre au questionnaire :material-arrow-right:](survey.md){ .md-button .md-button--primary }
