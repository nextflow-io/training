# Résumé du cours

Félicitations pour avoir terminé le cours de formation Nextflow pour la Génomique ! 🎉

## Votre parcours

Vous avez commencé par exécuter manuellement des outils d'appel de variants dans le terminal pour comprendre la méthodologie.
Ensuite, vous avez construit un pipeline Nextflow pour un échantillon unique afin d'automatiser le processus, l'avez étendu pour gérer plusieurs échantillons en parallèle, et avez ajouté le génotypage conjoint multi-échantillons en utilisant des opérateurs de canaux.

### Ce que vous avez construit

- Un pipeline d'appel de variants qui prend des fichiers BAM en entrée et produit des fichiers VCF appelés conjointement en sortie.
- Trois processus (`SAMTOOLS_INDEX`, `GATK_HAPLOTYPECALLER` et `GATK_JOINTGENOTYPING`) stockés dans des fichiers de modules séparés.
- Le pipeline s'adapte automatiquement à n'importe quel nombre d'échantillons en entrée grâce au paradigme de flux de données de Nextflow.
- Les résultats sont publiés dans un répertoire appelé `results/`.

### Compétences acquises

À travers ce cours pratique, vous avez appris à :

- Écrire un workflow linéaire pour appliquer l'appel de variants à un échantillon unique
- Gérer correctement les fichiers accessoires tels que les fichiers d'index et les ressources du génome de référence
- Exploiter le paradigme de flux de données de Nextflow pour paralléliser l'appel de variants par échantillon
- Implémenter l'appel conjoint multi-échantillons en utilisant les opérateurs de canaux appropriés
  Vous êtes maintenant équipé·e pour commencer à appliquer Nextflow aux workflows d'analyse génomique dans votre propre travail.

## Prochaines étapes pour développer vos compétences

Voici nos principales suggestions pour la suite :

- Appliquer Nextflow à d'autres cas d'usage d'analyse scientifique avec [Nextflow for Science](../index.md)
- Commencer avec nf-core grâce à [Hello nf-core](../../hello_nf-core/index.md)
- Explorer des fonctionnalités plus avancées de Nextflow avec les [Quêtes secondaires](../../side_quests/index.md)

Enfin, nous vous recommandons de découvrir [**Seqera Platform**](https://seqera.io/), une plateforme cloud développée par les créateurs de Nextflow qui facilite encore davantage le lancement et la gestion de vos workflows, ainsi que la gestion de vos données et l'exécution d'analyses de manière interactive dans n'importe quel environnement.

## Obtenir de l'aide

Pour les ressources d'aide et le soutien de la communauté, consultez la [page d'aide](../../help.md).

## Enquête de satisfaction

Avant de continuer, veuillez prendre une minute pour répondre à l'enquête sur le cours ! Vos retours nous aident à améliorer nos supports de formation pour tous.

[Répondre à l'enquête :material-arrow-right:](survey.md){ .md-button .md-button--primary }
