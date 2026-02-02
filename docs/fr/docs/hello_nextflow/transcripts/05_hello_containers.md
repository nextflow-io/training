# Partie 5 : Hello Containers - Transcription

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Traduction assistée par IA - [en savoir plus et suggérer des améliorations](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

<div class="video-wrapper">
  <iframe width="560" height="315" src="https://www.youtube.com/embed/5PyOWjKnNmg?si=QinuAnFwFj-Z8CrO&amp;list=PLPZ8WHdZGxmXiHf8B26oB_fTfoKQdhlik" title="YouTube video player" frameborder="0" allow="accelerometer; autoplay; clipboard-write; encrypted-media; gyroscope; picture-in-picture; web-share" referrerpolicy="strict-origin-when-cross-origin" allowfullscreen></iframe>
</div>

!!!note "Notes importantes"

    Cette page affiche uniquement la transcription. Pour les instructions détaillées étape par étape, retournez au [matériel de formation](../05_hello_containers.md).

    Les numéros de section affichés dans la transcription sont fournis à titre indicatif uniquement et peuvent ne pas inclure tous les numéros de section des supports.

## Bienvenue

Bonjour, bienvenue dans la Partie Cinq de la formation Hello Nextflow.

Ce chapitre s'intitule Hello Containers. Nous allons parler de la façon dont Nextflow s'intègre avec des outils tels que Docker et Singularity pour utiliser des conteneurs logiciels afin de fournir des logiciels aux utilisateur·trices de votre pipeline.

Cela signifie que lorsque les personnes exécutent votre pipeline, elles n'ont pas besoin d'aller installer tous les différents outils elles-mêmes. Nextflow le fera pour elles.

Les conteneurs sont une technologie extrêmement puissante et cruciale pour la reproductibilité et la facilité d'utilisation. Nous allons commencer par faire une brève introduction aux conteneurs eux-mêmes, exécuter quelques commandes docker manuellement, puis nous prendrons ces mêmes conteneurs et les intégrerons dans notre pipeline Nextflow.

D'accord. Commençons.

Comme précédemment, commençons par charger le matériel de formation. Allez sur training.nextflow.io. Hello Nextflow, Chapitre Cinq, Hello Containers.

Je vais accéder à mon environnement Codespaces et sur la gauche nous voyons hello containers point nf.

Comme avant, c'est le même script avec lequel nous avons terminé le chapitre quatre précédent, donc il devrait vous sembler familier.

Nous avons nos paramètres de ligne de commande pour spécifier le fichier d'entrée et le nom du lot. Nous incluons nos trois modules, et nous avons notre workflow où nous exécutons les trois processus.

## 0. Échauffement : Exécuter hello-containers.nf

N'hésitez pas à exécuter ce workflow à nouveau et vérifiez qu'il produit les sorties que vous attendez. Pour l'instant, je vais en fait le fermer et plonger dans le terminal.

## 1. Utiliser un conteneur 'manuellement'

Pour commencer ce chapitre, nous allons faire un petit récapitulatif sur la technologie des conteneurs. Si vous êtes très habitué·e à docker ou singularity ou d'autres technologies de conteneurs, alors considérez ceci comme un rappel, ou n'hésitez pas à le sauter complètement.

Nextflow supporte de nombreux types différents de technologies de conteneurs. Cela inclut Docker, Singularity, Podman, Shifter, Charliecloud, et plus encore.

Dans cette formation, nous allons nous concentrer sur Docker. Il est pré-installé dans les code spaces et est l'une des technologies de conteneurs les plus populaires, en particulier si vous développez sur votre propre ordinateur ou votre propre portable.

Si vous travaillez dans un environnement académique sur un HPC partagé, vous pourriez constater que Singularity est disponible et non Docker. Ce n'est pas grave. Tous les concepts sont exactement les mêmes. Quelques commandes manuelles sont différentes, mais si vous comprenez Docker, vous comprendrez aussi singularity.

En fait, Singularity est également installé dans l'environnement Code Spaces. Donc si vous le souhaitez, vous pouvez essayer de faire les mêmes tâches en utilisant Singularity au lieu de Docker.

D'accord, donc qu'est-ce que la technologie de conteneurs ? L'idée derrière Docker est qu'il peut récupérer une image depuis une source distante. La télécharger sur votre machine locale puis créer un conteneur basé sur cette image.

Ce conteneur en cours d'exécution est un peu comme une machine virtuelle fonctionnant sur votre ordinateur. Il est isolé de votre environnement, et il vient préemballé avec un système d'exploitation et un ensemble de logiciels disponibles.

## 1.1. Télécharger l'image du conteneur

La syntaxe dont nous avons besoin pour récupérer une image préexistante est "docker pull". Donc je vais taper ça dans mon terminal, mais maintenant nous avons besoin d'une image avec laquelle jouer.

Vous pouvez construire des images vous-même. Vous pouvez les trouver sur des registres publics comme Docker Hub ou quay.io. Mais une très bonne façon d'obtenir des images rapidement est d'utiliser Seqera Containers.

C'est un service communautaire gratuit que nous avons construit en 2024, que vous pouvez utiliser sans connexion ni rien.

Si vous allez sur seqera.io/containers ou cliquez sur containers en haut ici, vous êtes présenté·e avec une interface de recherche et vous pouvez taper le nom de n'importe quel outil disponible dans Conda ou sur le Python Package Index.

Par défaut, il recherche dans les canaux Bioconda et Conda Forge, mais vous pouvez préfixer n'importe quel canal Conda. Je suis ici si vous voulez.

Pour un peu de plaisir, utilisons cowpy. Je vais taper cowpy. Cela me donne des résultats de Python Package Index et Conda Forge. Je vais cliquer dessus pour l'ajouter à mon conteneur. Je pourrais ajouter plusieurs paquets ici si je le voulais. Sélectionner Docker, sélectionner linux/amd64, et cliquer sur Get Container.

Cela construit l'image pour moi à la demande si elle n'a pas déjà été créée, et me donne une URL que je peux copier.

Si vous êtes intéressé·e, vous pouvez cliquer sur view Build Details, et cela vous amène à une page qui montre le fichier d'environnement conda qui a été utilisé et le journal de construction complet pour la construction, ainsi que les résultats de l'analyse de sécurité.

Si je retourne à mes code spaces, je peux maintenant coller ce nom de conteneur et appuyer sur entrée.

Docker télécharge maintenant toutes les différentes couches dans cette image de conteneur, et nous dit maintenant que cette image est disponible à l'utilisation.

## Télécharger une image Singularity

Si vous utilisez singularity, le processus est fondamentalement le même. Nous sélectionnons nos paquets d'image, sélectionnons cowpy. Maintenant nous choisissons Singularity au lieu de Docker et cliquons sur Get Container. Cela nous donne une URL d'image utilisant oras://. Ou si vous préférez, vous pouvez utiliser https:// en cochant cette case. Copiez cette URL. Maintenant allez dans Code Spaces. Nous avons en fait Apptainer installé dans cet espace, qui est le même que Singularity, mais ils sont aliasés l'un à l'autre. Donc je vais faire apptainer pull et ensuite je vais l'appeler cowpy sif, mais vous pouvez l'appeler comme vous voulez. Collez l'URL. Et cela va télécharger cette image pour moi.

Je pourrais faire ls -lh et voir cowpy.sif

Singularity est différent de Docker, en ce que singularity stocke toutes les images dans des fichiers plats, alors que Docker a un registre où il conserve toutes les couches séparément sur votre machine hôte, et il a un démon en cours d'exécution pour garder trace de tout cela.

## 1.2. Utiliser le conteneur pour exécuter cowpy en tant que commande unique

D'accord, revenons à Docker. Nous pouvons maintenant essayer d'exécuter cette image que nous avons créée en faisant docker run.

Je vais faire dash dash rm, qui fait simplement une exécution unique de l'image. Et je vais coller l'URL de l'image. Et puis finalement, vous terminez ceci avec une commande que vous voulez exécuter.

L'image que nous avons générée avait cowpy installé, donc essayons cowpy.

Voilà. Elle a exécuté notre commande. Je n'ai pas cowpy installé localement. Vous pouvez voir que si j'essaie de l'exécuter, ça n'existe pas. Cependant, dans cette commande, je l'ai exécuté en utilisant Docker et il a correctement généré cette sortie.

## 1.3. Utiliser le conteneur pour exécuter cowpy de manière interactive

Nous pouvons aller plus loin si nous le souhaitons et démarrer un conteneur de manière interactive et regarder à l'intérieur. Encore une fois, je fais "docker run dash dash rm". Maintenant je vais faire dash it, qui indique à Docker que nous voulons un terminal interactif. Je fais à nouveau l'URL de l'image, et cette fois, au lieu de faire cowpy, je vais faire bin bash parce que la commande que nous voulons exécuter est bash.

Cela nous amène dans ce conteneur en cours d'exécution et vous pouvez voir que l'invite a changé maintenant.

Si je fais LS slash vous pouvez voir que les répertoires ici sont différents.

Si j'ouvre un deuxième terminal ici sur le côté droit, qui s'exécute simplement dans GitHub Code Spaces et je fais LS slash, vous voyez que nous avons des répertoires comme workspaces et temp, alors qu'ici dans Docker c'est différent.

Donc cet environnement est complètement séparé dans Docker et isolé de mon environnement hôte. C'est une bonne chose, parce que cela isole l'exécution de cette commande dans l'image Docker et la garde reproductible entre différentes personnes sur différents systèmes hôtes.

Si vous voulez utiliser des données de votre système hôte dans l'image Docker, vous devez explicitement monter cela dans le conteneur.

Nous allons faire cela dans une seconde.

## 1.3.2. Exécuter la ou les commandes de l'outil souhaité

D'abord cependant, voyons si nous pouvons exécuter cowpy. Encore une fois, la commande est disponible maintenant directement sur la ligne de commande, et nous pouvons commencer à faire des choses plus complexes et passer des arguments. Hello containers et au lieu de la vache, faisons le pingouin tux. Voyons ce que nous avons d'autre.

Faisons cheese. Merveilleux. Que diriez-vous de Dragon et Cow ? Plutôt bien.

## 1.3.3. Quitter le conteneur

D'accord. Je ne peux pas faire beaucoup plus parce que je n'ai pas de données dans ce conteneur. Donc sortons de cette image en cours d'exécution et voyons si nous pouvons monter des données dans le conteneur. Je peux faire ça en faisant control D ou en tapant exit. D'accord, je suis maintenant de retour dans mon GitHub code space régulier.

## 1.3.4. Monter des données dans le conteneur

Pour monter des données dans le conteneur Docker, je dois utiliser dash V. Donc je vais prendre ma commande docker précédente, retourner au début faire dash v. Je vais faire "." pour le répertoire de travail local actuel, puis deux points pour dire où cela devrait être monté dans le répertoire hôte et faire slash data. Donc cela monte ce répertoire particulier dans le conteneur à slash data.

Maintenant si je fais LS slash nous pouvons voir que nous avons un nouveau répertoire appelé data, et si je fais LS data, vous pouvez voir tous les fichiers que nous avons dans la barre latérale ici. Fantastique.

## 1.3.5. Utiliser les données montées

Maintenant nous pouvons commencer à utiliser certains des fichiers qui sont sur le système hôte dans l'image Docker. Donc je peux dire cat data greetings csv. Si vous vous souvenez, c'est notre fichier CSV avec nos différentes salutations d'avant, et je peux le piper vers cowpy. Fantastique. Maintenant nous progressons.

D'accord. C'est assez pour exécuter Docker de manière interactive. J'espère que vous avez maintenant une idée de ce qu'est Docker approximativement et comment l'utiliser à la fois pour exécuter une commande de manière unique, et aussi pour utiliser une image de manière interactive. Si vous utilisez singularity. Les commandes sont toutes très similaires sauf que vous faites des choses comme apptainer exec ou apptainer run, ou singularity exec ou singularity run.

## 2. Utiliser des conteneurs dans Nextflow

Ensuite, nous allons retourner à notre workflow Nextflow et voir comment utiliser cette technologie dans le pipeline Nextflow.

Fermons le terminal et ouvrons à nouveau Hello Containers.

## 2.1. Écrire un module cowpy

Pour rester avec notre exemple cowpy, créons un nouveau processus dans notre workflow, qui utilise cowpy. Allons dans modules, créons un nouveau fichier et appelons-le cowpy nf. Je vais maintenant tricher un peu et copier le code du matériel de formation et appuyer sur sauvegarder. Et regardons.

Donc c'est un processus simple. J'espère que maintenant vous comprenez à quoi ressemblent les éléments constitutifs d'un processus. Nous avons notre publishDir à nouveau, allant vers results. Nous avons deux entrées, un fichier d'entrée et une chaîne appelée character. Nous avons une sortie cowpy input file, et nous avons un script qui ressemble exactement à ce que nous avons exécuté manuellement dans notre image docker il y a une seconde : cat pour imprimer un fichier, le piper vers cowpy, dire quel type de caractère cowpy nous voulons utiliser, et sortir cela vers le fichier de sortie, que nous passons comme sortie ici.

## 2.2. Ajouter cowpy au workflow

D'accord, retournons à notre workflow, importons ce nouveau processus. Donc cowpy depuis modules cowpy nf. Créons un nouveau paramètre pour que nous puissions spécifier quel caractère nous voulions. Disons Turkey par défaut. Et puis appelons ce nouveau processus à la fin du workflow,

cowpy. Et utilisons la sortie ici de Collect Greetings. Donc collect greetings out, out file ici. Et puis nous avons besoin d'un second argument, qui est le nouveau paramètre que nous venons de créer. params dot character.

## 2.2.4. Exécuter le workflow pour vérifier qu'il fonctionne

D'accord, voyons si notre nouveau processus fonctionne. Nextflow run hello containers. Cela devrait exécuter ces trois premiers processus puis essayer d'exécuter cowpy à la fin.

Nous avons une erreur. Ce qu'il dit ici, cowpy a eu une erreur et il a eu un statut de sortie 127 et effectivement, commande sh cowpy commande non trouvée.

Nous n'avons pas dit à Nextflow que nous avons une image Docker disponible pour cowpy, donc il a essayé de l'exécuter sur notre système hôte et nous n'avons pas cowpy installé sur notre système hôte, donc cela a déclenché une erreur.

## 2.3. Utiliser un conteneur pour l'exécuter

Donc ce que nous devons faire, c'est que nous devons dire à Nextflow que nous avons un conteneur disponible. Allons à notre processus cowpy et nous allons ajouter une nouvelle directive en haut du processus appelée container.

Nous trouvons ensuite notre image, copions l'URL, et mettons cela dans une chaîne.

Ceci n'est pas suffisant en soi parce qu'un pipeline Nextflow peut avoir plusieurs façons de spécifier des logiciels. Je pourrais aussi faire conda conda-forge cowpy, par exemple. Et Nextflow a besoin de savoir laquelle de ces technologies vous voulez utiliser.

## 2.3.2. Activer l'utilisation de Docker via le fichier nextflow.config

Donc afin d'exécuter avec Docker activé, nous allons prendre un peu d'avance et utiliser le fichier de configuration Nextflow, qui est quelque chose que nous allons couvrir plus en détail dans le prochain chapitre. Vous pouvez voir dans ce répertoire que nous avons un fichier appelé Nextflow Config, et ici vous avez déjà docker.enabled False.

Nous allons changer cela en True pour activer Docker, et ensuite nous pouvons essayer d'exécuter le workflow à nouveau.

## 2.3.3. Exécuter le workflow avec Docker activé

Nextflow run hello containers nf et cette fois cowpy s'est exécuté avec succès. Regardons dans Results. cowpy collected test et voilà notre Turkey. Merveilleux.

Donc en arrière-plan là, Nextflow savait qu'il avait un conteneur disponible pour ce processus.

Il a récupéré l'image et il a exécuté les commandes pour nous.

## 2.3.4. Inspecter comment Nextflow a lancé la tâche conteneurisée

Si vous êtes curieux·se, nous pouvons en fait voir exactement ce qu'il a fait en regardant dans le répertoire work. Si je fais code work, puis le hash et ensuite command run, qui si vous vous souvenez est le fichier réel qui est exécuté pour cette tâche, nous pouvons entrer et nous pouvons chercher une fonction appelée NXF launch. Et ici vous pouvez voir la commande docker exacte que Nextflow a utilisée, qui ressemble beaucoup à ce que nous faisions manuellement dans le terminal plus tôt. Docker run. Lier ce répertoire hôte dans le conteneur, puis spécifier l'URL du conteneur.

Donc il n'y a pas de magie ici. C'est juste que Nextflow fait automatiquement le gros du travail pour vous d'une manière qui signifie que vous pouvez facilement spécifier des conteneurs dans votre pipeline, qui sont ensuite facilement disponibles pour quiconque d'autre utilise qui exécute votre workflow. Et ces personnes n'ont plus à penser à gérer les logiciels pour exécuter votre pipeline d'analyse.

Très, très simple, très pratique, et aussi vraiment reproductible. Bon dans l'ensemble.

D'accord, excellent travail. C'est la fin du Chapitre Cinq. Rejoignez-nous dans la prochaine vidéo pour la partie six, qui est la dernière partie de cette formation Hello Nextflow, où nous parlerons de la configuration Nextflow plus en détail.

À la prochaine vidéo.

[Transcription de la vidéo suivante :octicons-arrow-right-24:](06_hello_config.md)
