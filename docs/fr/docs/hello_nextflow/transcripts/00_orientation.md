# Orientation - Transcription Vidéo

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Traduction assistée par IA - [en savoir plus et suggérer des améliorations](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

<div class="video-wrapper">
  <iframe width="560" height="315" src="https://www.youtube.com/embed/G3CV-FcV-rc?si=nyLvwhrSB2m1NPc5&amp;list=PLPZ8WHdZGxmXiHf8B26oB_fTfoKQdhlik" title="YouTube video player" frameborder="0" allow="accelerometer; autoplay; clipboard-write; encrypted-media; gyroscope; picture-in-picture; web-share" referrerpolicy="strict-origin-when-cross-origin" allowfullscreen></iframe>
</div>

!!!note "Note importante"

    Cette page ne présente que la transcription. Pour les instructions complètes étape par étape, retournez au [matériel de formation](../00_orientation.md).

## Bienvenue

Bonjour, bienvenue à Hello Nextflow. Je m'appelle Phil Ewels. Je suis Chef de Produit pour l'Open Source chez Seqera, et je suis ravi d'être ici aujourd'hui pour vous guider à travers ce premier cours de formation Nextflow.

Nous allons passer en revue les bases de Nextflow, en expliquant comment écrire et exécuter des pipelines et les configurer.

Et vous allez construire votre propre pipeline simple à plusieurs étapes. Nous couvrirons la terminologie comme les opérateurs et les fabriques de canaux, et à la fin du cours, vous serez prêt·e à commencer à construire vos propres pipelines bioinformatiques.

Si vous avez des questions, n'hésitez pas à nous contacter sur community.seqera.io. Nous avons une communauté Nextflow très active, il y a une section dédiée à la formation, alors faites-nous savoir où vous êtes bloqué·e et quelqu'un pourra vous aider.

Très bien. Commençons.

## Site Web de Formation

Tout le matériel de formation pour les cours Nextflow se trouve sur training.nextflow.io. Vous pouvez y accéder dans votre navigateur web. Alors lancez-le maintenant et nous pouvons y jeter un œil.

Je vais utiliser la version 2.1.1. Nous publions de petites mises à jour et corrections de temps en temps, donc ne vous inquiétez pas si c'est légèrement différent, mais si le matériel a trop dérivé, vous pouvez toujours utiliser ce sélecteur de version en haut pour choisir la version exacte des matériaux que je vais présenter.

Si vous préférez le mode clair, vous pouvez changer le thème du site web ici.

Voyez les traductions ici, bien qu'au moment de l'enregistrement, c'est vraiment uniquement en anglais qui couvre ce nouveau matériel.

Et voyez également tout le code source du site web de formation et tout ce avec quoi nous allons travailler sur GitHub.

La page d'accueil ici liste tous les différents cours de matériel de formation que nous avons. Donc si je défile vers le bas, nous verrons Nextflow pour les nouveaux venus avec le cours Hello Nextflow que nous allons faire ici. Vous pouvez voir tous les autres cours que nous avons également, qui fonctionnent de manière similaire.

## Configuration de l'Environnement

Je vais en fait commencer par utiliser ce premier en haut, qui est commun à tous les cours de formation, et concerne spécifiquement la configuration de notre environnement.

Je clique dessus, cela m'amène à cette section, et nous pouvons voir les instructions pour le développement en local. Si vous voulez utiliser votre propre ordinateur portable avec votre propre copie de VS Code et vos propres installations logicielles, ou ce que nous attendons de la plupart des gens, qui est d'utiliser quelque chose appelé GitHub Codespaces.

Codespaces est un service fourni par GitHub où ils exécutent un serveur web dans le cloud, auquel vous pouvez vous connecter. Ce serveur a VS code installé, où vous pouvez l'exécuter dans votre navigateur web, ou si vous préférez, le connecter à votre installation locale de VS code. Tous les calculs, tous les fichiers, toutes les modifications se font à distance, ce qui signifie que tous les logiciels dont vous avez besoin sont pré-installés et sont les mêmes pour tout le monde.

## Créer un GitHub Codespace

Pour créer le codespace avec tout ce dont nous avons besoin, cherchez les boutons dans la documentation, qui disent "Open in GitHub Codespaces". Je vais cliquer dessus maintenant, l'ouvrir dans un nouvel onglet. Et je suis présenté avec cette page web. Maintenant vous pouvez voir que c'est pré-configuré avec nextflow-io training.

Je peux simplement cliquer sur créer un nouveau codespace. Mais en fait nous recommandons d'utiliser une machine légèrement plus grande pour la formation Nextflow avec quatre CPUs au lieu de deux. Vous pouvez changer quelle version du matériel il utilise. Donc cela utilise par défaut 2.1.1 parce que c'est la version de la documentation d'où j'ai suivi le lien. Mais je pourrais aussi le définir à une branche spécifique du dépôt si je le veux.

Maintenant je vais cliquer sur créer un codespace. Et il va commencer à configurer l'environnement pour moi.

## Création du Codespace

Maintenant, la première fois que vous faites cela, cela va prendre pas mal de temps, donc c'est le bon moment pour aller prendre une tasse de thé. Installez-vous confortablement, discutez avec la personne assise à côté de vous.

Si vous êtes intéressé·e, vous pouvez cliquer sur building codespace ici en bas pour voir les logs de la configuration. Et vous pouvez voir ici qu'il télécharge une image Docker avec tout ce dont j'ai besoin et configure l'environnement.

Maintenant, vous devez seulement attendre comme ça la première fois que vous créez un codespace. Si vous allez sur github.com/codespaces ici, vous verrez tous les différents Codespaces que vous avez ouverts. Voici celui que je viens de créer. La prochaine fois que vous faites cela, vous pouvez aller ici et vous pouvez sélectionner le codespace précédent et y revenir directement. Et c'est un processus beaucoup, beaucoup plus rapide pour réchauffer cet environnement existant. Cela conservera également toutes les modifications que vous avez apportées à VS Code et aux fichiers, donc vous ne perdrez pas votre progression si vous partez et revenez.

Vous pouvez cliquer sur les trois points ici pour effectuer d'autres actions. Par exemple, si vous l'avez configuré avec deux CPUs et maintenant vous en voulez quatre, vous pouvez changer le type de machine. Ou si vous voulez recommencer à zéro et frais, vous pouvez supprimer le codespace.

## Introduction à VS Code

D'accord, Codespaces a fini de configurer mon environnement et je suis maintenant présenté avec VS Code dans le navigateur web.

Si vous êtes habitué·e à VS code, cela vous semblera très familier. Si vous ne l'avez pas utilisé auparavant, c'est assez simple. Il y a quelques parties différentes de la page dont vous devez être conscient.

Ici sur la gauche, nous avons la barre latérale. Vous pouvez voir l'Explorateur configuré avec tous les différents fichiers dans le dépôt GitHub du dépôt de formation.

Sur ces boutons en bas à gauche, peuvent être différents outils. Dans la barre latérale. Je peux rechercher tous les fichiers dans tout le projet. Je peux travailler avec Git, je peux travailler avec GitHub, toutes sortes de choses comme ça.

En haut ici se trouve le menu principal. L'explorateur de fichiers est celui que nous aurons le plus ici, et vous pouvez faire un clic droit sur n'importe lequel de ces fichiers et faire les choses normales auxquelles vous vous attendez. Vous devrez peut-être cliquer à travers quelques avertissements comme celui-ci où il coupe copie et vous pouvez télécharger sur votre machine locale également.

Quand le codespace se charge, il nous donne un aperçu du fichier markdown dans cette zone principale ici. C'est exactement le même que celui qui s'affiche sur github.com. Je peux fermer ça et si je double-clique sur ce fichier Readme, vous verrez qu'il l'ouvre sous forme de code dans l'éditeur de code et tout comme avec n'importe quel autre fichier, nous pouvons éditer ce code directement.

Enfin en bas ici, nous avons la fenêtre de terminal. Je regardais les logs pendant la construction, donc c'est ce qu'elle affiche actuellement. Je peux aussi appuyer sur ce bouton plus pour démarrer une nouvelle session de terminal. Cela ne s'exécute pas sur ma machine. Rappelez-vous, cela s'exécute dans le cloud, et si je fais tree trois à une profondeur de deux, vous verrez tous les mêmes fichiers ici, qui étaient sur la gauche.

## Afficher uniquement les fichiers "hello-nextflow"

Ce dépôt GitHub contient tous les différents ensembles de formation, pas seulement celui que nous faisons. Donc si vous voulez, vous pouvez vous concentrer uniquement sur le dossier Hello Nextflow. Une façon de nettoyer un peu cela est d'aller dans le menu fichier puis ajouter un dossier à l'espace de travail.

Nous cliquons dessus, allons dans training. Hello nextflow, et cliquons sur ajouter. Cela actualisera votre écran. Et puis dans l'Explorateur, nous avons maintenant deux espaces de travail différents, celui que nous avions avant pour training et un avec juste Hello Nextflow.

Si vous voulez, vous pouvez faire un clic droit sur training et cliquer sur supprimer le dossier de l'espace de travail pour le retirer complètement de la barre latérale.

Maintenant nous avons juste les fichiers pour ce cours de formation particulier dans la barre latérale. Je peux masquer cet avertissement et maintenant je peux faire la même chose dans le terminal ici et faire CD pour changer de répertoire. Hello, Nextflow. Et encore, nous avons les mêmes fichiers ici, qui sont dans la barre latérale.

## Hello Nextflow : fichiers

En regardant ces fichiers pour le cours Hello Nextflow.

Nous avons un ensemble de fichiers .nf, qui sont pour Nextflow, et il y a un de ces fichiers pour chacun des chapitres du cours de formation. Nous allons travailler à travers ces fichiers et les modifier dans les exercices.

Nous avons également un fichier nextflow.config, qui a juste des paramètres de configuration de base pour exécuter Nextflow dans cet environnement, dont vous n'avez pas vraiment besoin de vous soucier à ce stade. Un fichier greetings.csv, que nous utiliserons pour traiter les données, qui sera introduit dans la prochaine partie de ce cours, et un fichier test-params.json, qui sera utilisé dans la partie six et que vous pouvez ignorer pour l'instant.

Ces fichiers Nextflow sont juste le début de chaque exercice. Si vous voulez voir à quoi ils devraient ressembler une fois terminés, vous pouvez aller dans un répertoire solutions et il y a les réponses pour chaque partie du cours de formation, donc vous pouvez voir une version fonctionnelle de ce vers quoi vous visez.

## Ouvrir un Terminal

Si à un moment donné vous fermez le terminal et ne vous rappelez plus comment y revenir, ne vous inquiétez pas. Ces boutons en haut à droite ouvrent et ferment différents panneaux dans l'espace de travail. Donc cliquez sur celui-ci pour le panneau du bas et il réapparaîtra. Et assurez-vous simplement que vous avez sélectionné terminal ici. Vous pouvez également cliquer sur ce bouton ici, la flèche sur le côté droit d'un terminal pour le mettre en plein écran.

Vous me verrez faire ça assez souvent parce que j'ai VS Code zoomé pour que vous puissiez lire le texte. Selon la taille de votre écran, vous pourriez avoir besoin ou non de faire cela. Il en va de même pour minimiser le panneau latéral.

Très bien. C'est suffisant pour l'environnement. Je pense que nous sommes prêts à commencer. Rejoignez-moi dans la prochaine vidéo pour le chapitre un.

[Transcription vidéo suivante :octicons-arrow-right-24:](01_hello_world.md)
