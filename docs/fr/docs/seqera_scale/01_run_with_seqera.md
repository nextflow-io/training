# Partie 1 : Lancer des pipelines depuis l'interface web

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Traduction assistée par IA - [en savoir plus et suggérer des améliorations](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

Dans cette partie du cours Scale with Seqera, vous allez configurer l'accès à Seqera Platform et lancer un pipeline à l'échelle de la production depuis l'interface web.

Assurez-vous que votre répertoire de travail est défini sur `seqera-scale/` comme indiqué sur la page [Premiers pas](./00_orientation.md).

---

## 1. Premiers pas avec Seqera

Seqera fournit une plateforme complète pour lancer, surveiller et gérer des pipelines Nextflow.
Cette section vous guide à travers l'inscription et la prise en main avant d'exécuter votre premier pipeline.

### 1.1. Créer un compte gratuit

Rendez-vous sur [cloud.seqera.io](https://cloud.seqera.io) et créez un compte gratuit.
Vous pouvez vous inscrire avec votre adresse e-mail, GitHub ou vos identifiants Google.

Un compte gratuit vous donne accès à :

- **Espace de travail personnel** : votre propre espace pour ajouter des pipelines, configurer des environnements de calcul et gérer les exécutions
- **Accès au Community Showcase** : une collection sélectionnée de pipelines nf-core et communautaires avec des paramètres préconfigurés et des données d'exemple

Consultez la [documentation Seqera](https://docs.seqera.io) pour un aperçu complet des niveaux de compte et des fonctionnalités disponibles.

### 1.2. Explorer le Community Showcase

Avant de lancer vos propres pipelines, prenez quelques minutes pour explorer le Community Showcase.
Il vous donne un aperçu réaliste de ce à quoi ressemble la plateforme avec de vrais pipelines et de vraies données.

1. Connectez-vous sur [cloud.seqera.io](https://cloud.seqera.io).
2. Dans la barre latérale gauche, cliquez sur **Showcase**.
3. Parcourez les pipelines disponibles — vous reconnaîtrez plusieurs pipelines nf-core du cours Use nf-core.
4. Cliquez sur un pipeline pour afficher sa configuration et ses paramètres de lancement.
5. Cliquez sur **Runs** pour explorer les historiques d'exécution, y compris les détails au niveau des tâches et les rapports des exécutions précédentes.

Il s'agit d'une vue en lecture seule, mais elle vous montre comment fonctionne l'interface avant que vous n'exécutiez quoi que ce soit vous-même.

### 1.3. Accéder à un espace de travail avec des ressources de calcul

Le lancement de pipelines nécessite un espace de travail avec un environnement de calcul configuré.

Seqera prend en charge deux façons de fournir des ressources de calcul :

- **Connecter votre propre infrastructure** : AWS, Azure, Google Cloud et les ordonnanceurs HPC (SLURM, LSF, PBS et autres).
  Consultez la [documentation sur les environnements de calcul](https://docs.seqera.io) pour les guides de configuration.
- **Seqera Compute** : un service géré qui fournit des environnements de calcul pré-provisionnés sur AWS, moyennant des frais, sans configuration de compte cloud requise.
  Vous pouvez l'activer directement depuis les paramètres de votre espace de travail.

**Formation en groupe :**
Si vous participez à une session de formation en groupe, vous avez peut-être été ajouté·e à une organisation et à un espace de travail dont les ressources de calcul sont déjà configurées.
Votre formateur·trice vous communiquera le nom de l'organisation, le nom de l'espace de travail et tout autre détail dont vous avez besoin.

**Travail en autonomie :**
Si vous suivez cette formation par vous-même, vous devrez configurer un environnement de calcul dans votre espace de travail personnel en utilisant l'une des options ci-dessus.
Des crédits gratuits pour essayer Seqera Compute sont [disponibles sur demande](https://seqera.io/platform/compute/).

!!! note "Note"

    La suite de ce cours suppose que vous avez accès à un espace de travail avec un environnement de calcul configuré.
    Si vous participez à une session de formation en groupe, votre formateur·trice vous confirmera quel espace de travail et quel environnement de calcul utiliser.

### À retenir

Vous disposez d'un compte Seqera, vous avez exploré le Community Showcase et vous êtes en mesure d'accéder à un espace de travail avec des ressources de calcul.

### Et ensuite ?

Lancer un pipeline RNA-seq à l'échelle de la production depuis l'interface web de Seqera Cloud.

---

## 2. Lancer nf-core/rnaseq depuis l'interface web

Comme abordé dans Use nf-core, le pipeline nf-core/rnaseq est un pipeline développé par la communauté pour l'analyse de données de séquençage RNA en masse.

Dans cette section, vous allez ajouter le pipeline à votre espace de travail, lancer une exécution et surveiller son déroulement.

### 2.1. Ajouter le pipeline à votre espace de travail

nf-core/rnaseq fait partie d'une collection sélectionnée de pipelines qui peuvent être ajoutés à votre espace de travail en quelques clics via le service Seqera Pipelines.

_Nous vous montrerons comment ajouter vos propres pipelines plus loin dans ce cours._

1. Rendez-vous sur [**Seqera Pipelines**](https://seqera.io/pipelines) pour parcourir la collection communautaire.
2. Recherchez `rnaseq` et sélectionnez **nf-core/rnaseq**.
3. Cliquez sur **Launch Pipeline** ou faites défiler jusqu'en bas de la page jusqu'à la section **Launch Pipeline**.
4. Assurez-vous d'être connecté·e et sélectionnez les valeurs appropriées dans les menus déroulants **Organizations**, **Workspace** et **Compute Environment**.
   **Astuce pour les groupes :** Si vous utilisez un espace de travail partagé, ajoutez un identifiant unique (comme votre nom d'utilisateur·trice) au nom du pipeline.
5. Cliquez sur **Add pipeline to your Seqera account**

Une boîte apparaîtra avec le message : **Pipeline added: View Pipeline**.
Cliquer sur le lien vous amènera à l'entrée du pipeline dans votre launchpad.

Le pipeline est maintenant répertorié dans le panneau **Launchpad** de votre espace de travail et est prêt à être lancé.

### 2.2. Lancer le pipeline

Cliquez sur le bouton **Launch** du pipeline, soit dans le panneau **Launchpad**, soit sur la page de détails du pipeline.
Cela ouvre l'interface de configuration.

Le pipeline est déjà configuré avec le profil `test`, de sorte que les données d'entrée, le répertoire de sortie et la référence du génome sont pré-remplis.
Vous pouvez ignorer le reste des paramètres et des paramètres avancés pour l'instant.

Cliquez sur le bouton bleu **Launch** pour démarrer effectivement l'exécution.

### 2.3. Surveiller l'exécution

Après le lancement, vous serez redirigé·e vers le panneau **Runs** de votre pipeline.

La vue d'exécution affiche :

- **Status** : l'état actuel de l'exécution (soumise, en cours, réussie, échouée)
- **Command line** : la commande `nextflow run` exacte que la plateforme a construite et soumise
- **Parameters** : toutes les valeurs de paramètres utilisées pour cette exécution
- **Tasks** : un tableau de chaque appel de processus, avec le statut, la durée et l'utilisation des ressources

Cliquez sur n'importe quelle ligne de tâche pour inspecter ses détails d'exécution, notamment :

- Le script `.command.sh` qui a été exécuté
- Les journaux stdout et stderr
- Les métriques CPU, mémoire et E/S

L'onglet **Reports** affichera un rapport MultiQC une fois l'exécution terminée, agrégeant les métriques de contrôle qualité pour tous les échantillons.

Cela prendra un certain temps, nous allons donc continuer pour l'instant et revenir plus tard pour examiner les sorties, etc.

### À retenir

Vous savez comment ajouter un pipeline à un espace de travail Seqera, configurer et lancer une exécution, et surveiller l'exécution à grande échelle.

### Et ensuite ?

Passez à la [Partie 2](./02_launch_from_cli.md), où vous apprendrez à effectuer tout cela depuis la ligne de commande en utilisant le CLI `tw`.

---

## Résumé

Dans cette partie, vous avez appris à :

- Créer un compte Seqera et explorer le Community Showcase
- Ajouter un pipeline depuis le catalogue sélectionné, lancer une exécution à l'échelle de la production et surveiller son déroulement
