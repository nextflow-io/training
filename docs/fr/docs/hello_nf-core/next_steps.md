# Résumé du cours

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Traduction assistée par IA - [en savoir plus et suggérer des améliorations](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

Félicitations pour avoir terminé le cours de formation Hello nf-core ! 🎉

<!-- placeholder for video -->

## Votre parcours

Vous avez commencé par apprendre à récupérer et exécuter un pipeline de démonstration, puis vous avez abordé la conversion d'un simple workflow Nextflow en pipeline nf-core.
Vous avez appris à créer une structure de pipeline à l'aide d'un modèle et à greffer le pipeline existant sur cette structure.
Ensuite, vous avez progressivement affiné le pipeline en remplaçant l'un des modules locaux par un module nf-core, transformé un autre module local pour qu'il respecte les normes nf-core, et ajouté la validation des entrées.

### Ce que vous avez construit

Votre pipeline final `core-hello` possède maintenant :

- **Une structure standardisée** utilisant le modèle nf-core avec des répertoires organisés pour les workflows, les subworkflows, les modules et la configuration
- **Des modules communautaires** du dépôt nf-core (`cat/cat`) aux côtés de vos modules personnalisés
- **Une validation complète** qui vérifie à la fois les paramètres et les données d'entrée avant l'exécution du pipeline
- **Une configuration professionnelle** avec des profils pour différents environnements d'exécution
- **Une documentation complète** et des métadonnées suivant les conventions nf-core

### Compétences clés acquises

Grâce à ce cours pratique, vous avez appris à :

1. **Naviguer et comprendre** la structure d'un pipeline nf-core en explorant un pipeline existant
2. **Restructurer les workflows** pour les rendre composables et les adapter au modèle nf-core
3. **Trouver et intégrer** des modules pré-construits du dépôt communautaire
4. **Créer des modules personnalisés** en suivant les normes nf-core pour le nommage, la structure et les métadonnées
5. **Mettre en œuvre la validation** en utilisant nf-schema pour détecter les erreurs tôt avec des retours clairs

Vous êtes maintenant équipé·e des connaissances fondamentales pour construire des pipelines nf-core prêts pour la production qui suivent les meilleures pratiques de la communauté.

## Prochaines étapes pour développer vos compétences

Voici nos 3 principales suggestions pour la suite :

- Appliquez Nextflow à un cas d'usage d'analyse scientifique avec [Nextflow for Science](../nf4_science/index.md)
- Explorez des fonctionnalités Nextflow plus avancées avec les [Side Quests](../side_quests/index.md)
- Impliquez-vous en [rejoignant la communauté nf-core](https://nf-co.re/join).

Enfin, nous vous recommandons de jeter un œil à [**Seqera Platform**](https://seqera.io/), une plateforme cloud développée par les créateurs de Nextflow qui facilite encore davantage le lancement et la gestion de vos workflows, ainsi que la gestion de vos données et l'exécution d'analyses interactives dans n'importe quel environnement.

## Questionnaire de satisfaction

Avant de continuer, veuillez prendre une minute pour compléter le questionnaire du cours ! Vos retours nous aident à améliorer nos supports de formation pour tous.

[Répondre au questionnaire :material-arrow-right:](survey.md){ .md-button .md-button--primary }
