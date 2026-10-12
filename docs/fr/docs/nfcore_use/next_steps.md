# Résumé du cours

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Traduction assistée par IA - [en savoir plus et suggérer des améliorations](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

Félicitations pour avoir terminé le cours de formation Use nf-core ! 🎉

<!-- placeholder for video -->

## Votre parcours

Vous avez commencé par trouver et récupérer le pipeline `nf-core/demo`, puis vous avez appris à l'exécuter en utilisant son profil de test et à examiner ses sorties.
Ensuite, vous avez configuré son exécution via les paramètres du pipeline et les fichiers de configuration, et vous avez vu comment les pipelines nf-core valident les paramètres et les données d'entrée.
Enfin, vous avez appliqué ces mêmes compétences à `nf-core/rnaseq`, un pipeline à l'échelle de la production, et vous avez appris à remplacer ses allocations de ressources par défaut pour les adapter au matériel dont vous disposez.

### Ce que vous avez appris

Vous êtes maintenant en mesure de trouver, récupérer, exécuter et configurer des pipelines nf-core.

- Les pipelines nf-core sont récupérés avec `nextflow pull` et suivent une organisation du code standardisée.
- Chaque pipeline nf-core est livré avec un profil `test` pour une validation rapide sur un petit jeu de données.
- Les paramètres du pipeline (définis via `--param_name` ou `-params-file`) et la configuration (définie via `-c`) servent des objectifs différents : les entrées et les options d'analyse d'un côté, la logistique d'exécution comme l'allocation des ressources de l'autre.
- Les pipelines nf-core valident automatiquement les paramètres et les fichiers d'entrée, détectant les erreurs avant qu'aucun travail ne soit effectué.
- Les ressources par défaut sont assignées via des labels (`process_low`, `process_medium`, `process_high`) définis dans `conf/base.config`, que vous pouvez remplacer avec un fichier de configuration personnalisé.

### Compétences acquises

Au cours de cette formation pratique, vous avez appris à :

- Trouver un pipeline nf-core sur le site nf-co.re et récupérer son code source
- Exécuter un pipeline en utilisant son profil de test intégré et examiner ses sorties
- Obtenir de l'aide, définir des paramètres et comprendre la validation des paramètres et des entrées
- Personnaliser l'allocation des ressources et les arguments des outils via des fichiers de configuration
- Récupérer et exécuter un pipeline à l'échelle de la production, et remplacer ses labels de ressources par défaut

Vous êtes maintenant équipé·e des connaissances fondamentales pour commencer à exécuter des pipelines nf-core pour vos propres analyses.

## Prochaines étapes pour développer vos compétences

Voici nos meilleures suggestions pour la suite :

- Lancez et surveillez ces pipelines à grande échelle avec [Scale with Seqera](../seqera_scale/index.md)
- Ne vous contentez pas d'exécuter des pipelines nf-core, développez-en ! Apprenez les bonnes pratiques nf-core avec [Build with nf-core](../nfcore_build/index.md)
- Vous êtes nouveau·elle sur Nextflow ? Commencez par [Nextflow Run](../nextflow_run/index.md)
- Appliquez Nextflow à un cas d'usage d'analyse scientifique avec [Nextflow for Science](../nf4_science/index.md)
- Explorez des fonctionnalités Nextflow plus avancées avec les [Quêtes secondaires](../side_quests/index.md)

## Obtenir de l'aide

Pour les ressources d'aide et le soutien de la communauté, consultez la [page d'aide](../help.md).

## Enquête de satisfaction

Avant de continuer, veuillez prendre un moment pour compléter l'enquête du cours ! Vos retours nous aident à améliorer nos supports de formation pour tout le monde.

[Répondre à l'enquête :material-arrow-right:](survey.md){ .md-button .md-button--primary }
