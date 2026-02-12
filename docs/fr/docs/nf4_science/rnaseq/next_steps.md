# Résumé du cours

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Traduction assistée par IA - [en savoir plus et suggérer des améliorations](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

Félicitations pour avoir terminé le cours de formation Nextflow pour RNAseq !

## Votre parcours

Vous avez commencé par exécuter manuellement des outils de traitement RNAseq dans le terminal pour comprendre la méthodologie.
Ensuite, vous avez construit un pipeline Nextflow pour un seul échantillon afin d'automatiser le processus, vous l'avez étendu pour gérer plusieurs échantillons en parallèle, et vous l'avez enrichi pour traiter des données paired-end et agréger des rapports QC sur plusieurs échantillons.

### Ce que vous avez construit

- Un pipeline de traitement RNAseq qui prend des fichiers FASTQ en entrée et produit des lectures nettoyées, des alignements et des rapports QC agrégés en sortie.
- Des processus pour le nettoyage (Trim Galore), l'alignement (HISAT2), le contrôle qualité (FastQC) et l'agrégation de rapports (MultiQC) stockés dans des fichiers de modules séparés.
- Le pipeline parallélise automatiquement le traitement des échantillons d'entrée en utilisant le paradigme de flux de données de Nextflow.
- Le pipeline final gère les données de séquençage paired-end.

### Compétences acquises

Grâce à ce cours pratique, vous avez appris à :

- Écrire un workflow linéaire pour appliquer des méthodes de base de traitement et de QC RNAseq
- Gérer de manière appropriée des fichiers spécifiques au domaine tels que les FASTQ et les ressources de génome de référence
- Gérer les données de séquençage single-end et paired-end
- Exploiter le paradigme de flux de données de Nextflow pour paralléliser le traitement RNAseq par échantillon
- Agréger des rapports QC sur plusieurs étapes et échantillons en utilisant les opérateurs de canaux appropriés

Vous êtes maintenant équipé·e pour commencer à appliquer Nextflow aux workflows d'analyse RNAseq dans votre propre travail.

## Prochaines étapes pour développer vos compétences

Voici nos principales suggestions pour la suite :

- Appliquer Nextflow à d'autres cas d'usage d'analyse scientifique avec [Nextflow for Science](../index.md)
- Démarrer avec nf-core avec [Hello nf-core](../../hello_nf-core/index.md)
- Explorer des fonctionnalités Nextflow plus avancées avec les [Quêtes secondaires](../../side_quests/index.md)

Enfin, nous vous recommandons de découvrir [**Seqera Platform**](https://seqera.io/), une plateforme cloud développée par les créateurs de Nextflow qui facilite encore davantage le lancement et la gestion de vos workflows, ainsi que la gestion de vos données et l'exécution d'analyses de manière interactive dans n'importe quel environnement.

## Obtenir de l'aide

Pour les ressources d'aide et le support de la communauté, consultez la [page d'aide](../../help.md).

## Enquête de satisfaction

Avant de continuer, prenez une minute pour répondre à l'enquête du cours ! Vos retours nous aident à améliorer nos supports de formation pour tout le monde.

[Répondre à l'enquête :material-arrow-right:](survey.md){ .md-button .md-button--primary }
