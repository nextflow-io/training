# Partie 6 : Hello Config - Transcription

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Traduction assistée par IA - [en savoir plus et suggérer des améliorations](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

<div class="video-wrapper">
  <iframe width="560" height="315" src="https://www.youtube.com/embed/IuDO2HeKvXk?si=tnXTi6mRkITY0zW_&amp;list=PLPZ8WHdZGxmXiHf8B26oB_fTfoKQdhlik" title="YouTube video player" frameborder="0" allow="accelerometer; autoplay; clipboard-write; encrypted-media; gyroscope; picture-in-picture; web-share" referrerpolicy="strict-origin-when-cross-origin" allowfullscreen></iframe>
</div>

!!!note "Notes importantes"

    Cette page présente uniquement la transcription. Pour des instructions complètes étape par étape, retournez au [matériel de formation](../06_hello_config.md).

    Les numéros de section affichés dans la transcription sont fournis à titre indicatif uniquement et peuvent ne pas inclure tous les numéros de section dans les supports.

## Bienvenue

Bonjour, bienvenue dans la partie six de la formation Hello Nextflow.

Ce chapitre s'intitule Hello Config, et c'est la dernière partie de notre formation.

Dans ce chapitre, nous allons parler de la configuration Nextflow. La configuration Nextflow est vraiment puissante. Elle nous permet d'exécuter le même pipeline sur plusieurs infrastructures de calcul différentes avec différentes méthodes de provisionnement logiciel et différentes options dans le pipeline lui-même.

Cela signifie que vous pouvez prendre des pipelines Nextflow construits par d'autres personnes et les exécuter sur votre système, même s'ils ont pu être conçus pour une infrastructure entièrement différente. Cette capacité à configurer Nextflow rend les workflows vraiment portables et partageables.

Dans ce chapitre, nous utiliserons le workflow que nous avons construit dans les parties précédentes, mais nous n'allons pas du tout modifier le code du workflow. Nous allons simplement examiner notre fichier de configuration Nextflow et voir comment la modification de la configuration altère la façon dont Nextflow s'exécute.

D'accord, commençons.

Comme auparavant, commençons par aller sur training.nextflow.io. Allez à gauche sur Hello Nextflow et chapitre six, Hello config. Je vais maintenant aller dans mon environnement GitHub Codespaces et vérifier le script que nous allons utiliser.

## 0. Échauffement : Vérifier que Docker est activé et exécuter le workflow Hello Config

Celui-ci s'appelle Hello Config, et il repart d'où nous étions auparavant. Il ressemble donc exactement au même avec nos trois paramètres. Greetings pour le fichier CSV, batch pour le nom du lot de sortie et character pour le nom cowpy. Nous avons nos quatre imports des différents processus, puis nous avons un workflow où nous les enchaînons.

Je vais effectivement fermer ce fichier maintenant car nous ne toucherons pas du tout au fichier Nextflow dans ce chapitre. Nous allons travailler uniquement dans le fichier de configuration. Si je regarde dans le fichier nextflow.config que nous avons brièvement examiné dans le chapitre cinq précédent, nous pouvons voir que nous avons une seule instruction ici : docker enabled equals true, qui indique à Nextflow d'utiliser Docker lorsqu'il exécute ce workflow.

J'utilise nextflow.config à la racine du pipeline ici, qui est chargé automatiquement lorsque j'exécute Nextflow. Mais rappelez-vous, Nextflow peut charger des fichiers de configuration depuis plusieurs endroits.

Si je vérifie avec la documentation Nextflow, je vais à Configuration, vous pouvez voir une liste de ces endroits et une priorité dans laquelle ils se chargent.

D'accord. Vérifions que notre workflow s'exécute comme nous l'attendons. Je fais apparaître un terminal. Je fais nextflow run hello-config et j'appuie sur entrée. Nous devrions avoir ces quatre processus en cours d'exécution, se terminant par une commande cowpy. Effectivement, cela a bien fonctionné. Docker était activé, il a téléchargé Docker et a exécuté cowpy pour moi, comme à la fin du chapitre cinq.

## 1. Déterminer quelle technologie d'empaquetage logiciel utiliser

D'accord. Disons que je travaille sur un HPC et que je n'ai pas Docker installé. La meilleure chose à faire dans ce scénario serait d'utiliser Singularity ou Apptainer. Si je devais faire cela, j'irais dans le module cowpy et changerais ce conteneur pour utiliser l'image singularity comme je l'ai montré dans le chapitre précédent, avec un oras://, que vous pouvez également obtenir depuis Seqera Containers.

J'irais ensuite dans nextflow.config, mettrais docker enabled à false et ferais singularity enabled equals true. Ou, si j'utilise Apptainer, apptainer enabled equals true et cela fonctionnerait.

Nextflow prend également en charge d'autres technologies en plus des conteneurs, quelque chose que vous connaissez peut-être est conda. Ici, nous pouvons faire conda enabled equals true et mettre Docker à false. conda n'utilise pas la même directive container. Au lieu de cela, nous pouvons en ajouter une nouvelle ici appelée conda. Nous spécifions ensuite le paquet conda que nous voulons utiliser. C'est une bonne pratique d'être aussi spécifique que possible pour essayer de rendre le pipeline aussi reproductible que possible. Je vais donc spécifier le canal conda, conda-forge, puis cowpy, et la version exacte, qui était 1.1.5.

Je pourrais aussi simplement écrire cowpy si je voulais, mais cela pourrait se résoudre en une version différente de cowpy lors de différentes exécutions du pipeline.

La chose intéressante à ce sujet est que je n'ai pas du tout touché la directive docker. Cette image Docker est toujours là. Je fournis juste deux alternatives maintenant, et celles-ci peuvent être activées ou désactivées en utilisant uniquement un fichier de configuration.

## 1.3. Exécuter le workflow pour vérifier qu'il peut utiliser Conda

Conda est maintenant activé, alors essayons-le.

Excellent. Il s'exécute et vous pouvez voir qu'il y a un message de Nextflow ici disant que Nextflow crée un environnement conda pour moi, et il utilise cet emplacement de cache.

En arrière-plan, Nextflow exécute des commandes "conda create" pour moi pour créer un nouvel environnement conda isolé avec juste les paquets que je veux, puis installe et récupère ces paquets conda afin qu'il puisse exécuter le processus.

Vous pouvez voir que cela a pris un peu de temps parce qu'il créait l'environnement et installait le logiciel pour la première fois. Cependant, il a mis en cache cet environnement, donc si j'exécute à nouveau la même commande Nextflow, cela devrait être beaucoup plus rapide car il réutilisera le même environnement conda.

L'une des choses intéressantes à ce sujet est que ces directives peuvent être spécifiées au niveau du processus, pas seulement pour le workflow entier. Donc si vous le souhaitez, vous pouvez mélanger et assortir quelle technologie est utilisée pour différents processus.

## 2. Allouer des ressources de calcul avec des directives de processus

Le fichier de configuration Nextflow peut faire beaucoup plus que simplement l'empaquetage logiciel. Nous pouvons également indiquer à Nextflow comment exécuter réellement les étapes du pipeline. Un exemple est de dire à un système hôte quelles ressources doivent être mises à disposition pour chaque tâche en cours d'exécution.

Par défaut, Nextflow ne donne pas grand-chose. Il donne un seul CPU et seulement deux gigaoctets de mémoire à chaque processus.

C'est probablement quelque chose que nous voudrions changer, afin que les processus qui prennent beaucoup de temps à s'exécuter puissent avoir plus de ressources et s'exécuter plus rapidement, mais il peut être difficile de savoir quoi allouer à un processus. Nextflow a quelques astuces intéressantes pour vous aider avec cela.

## 2.1. Exécuter le workflow pour générer un rapport d'utilisation des ressources

Exécutons à nouveau le workflow. Cette fois, je vais ajouter un argument supplémentaire, qui est -with-report. C'est une option Nextflow de base, donc c'est un seul trait d'union. Et puis le nom de fichier que je veux. Dans ce cas, je vais l'appeler report-config-one.html.

Je vais exécuter à nouveau le workflow. Il va s'exécuter exactement comme avant, mais il va me donner un rapport d'aide supplémentaire, que vous pouvez voir apparaître maintenant ici dans la barre latérale.

Je vais faire un clic droit sur ce fichier, cliquer sur télécharger, ce qui le télécharge de GitHub Codespaces vers mon système local, afin que je puisse facilement le visualiser dans le navigateur web ici.

Ce rapport peut être généré pour n'importe quelle exécution Nextflow, et il contient beaucoup d'informations. Il commence en haut avec des métadonnées sur la commande utilisée, quand le workflow a été exécuté, combien de temps cela a pris, mais en faisant défiler vers le bas, nous obtenons des informations plus détaillées sur les ressources qui ont été utilisées par chaque étape du pipeline.

Parce que chaque processus s'exécute plusieurs fois pour différentes tâches, nous avons une boîte à moustaches montrant la variation des ressources que nous avons utilisées pour chaque processus.

Si je fais défiler un peu plus bas, je vois des informations similaires sur la mémoire utilisée et la durée des tâches. Également la lecture-écriture disque.

Vous pouvez imaginer que pour un grand pipeline avec des tâches de longue durée, cela peut être très informatif sur la façon de peaufiner la configuration des ressources que vous demandez afin de ne pas en demander trop, mais aussi de fournir suffisamment pour qu'il s'exécute rapidement.

Si je continue à faire défiler le rapport, nous voyons également un tableau de tâches, qui nous montre des informations détaillées sur chaque tâche unique qui a été exécutée dans le workflow. Cela inclut des informations telles que le script résolu qui a été exécuté.

D'accord, revenons à notre fichier de configuration. Nous avons vu que nous n'avions vraiment pas besoin de beaucoup pour notre workflow, alors disons à Nextflow que nous n'avons besoin que d'un gigaoctet de mémoire pour chaque processus dans le workflow.

Maintenant, lorsque nous le définissons ainsi au niveau du processus, cela s'applique à chaque processus du pipeline.

## 2.3. Définir des allocations de ressources pour un processus individuel

Pour les besoins de l'argument, supposons que cowpy fait vraiment beaucoup de calculs lourds et qu'il a besoin de plus de ressources que les autres tâches. Nous pouvons définir un bloc de configuration supplémentaire ici, qui s'applique uniquement à ce processus en utilisant withName cowpy.

C'est ce qu'on appelle un sélecteur de configuration, et nous pouvons définir différents motifs ici pour correspondre à différents processus. Par exemple, je pourrais faire cow star. Je suis ensuite cela avec des accolades et donnons-lui deux gigaoctets de mémoire au lieu d'un et disons deux CPUs.

Maintenant Nextflow donnera à chaque processus du workflow un gigaoctet à part cette requête, qui est plus spécifique. Donc elle l'écrase. Et juste pour les processus qui s'appellent cowpy, ils obtiendront deux gigs de mémoire et deux CPUs.

Notez que Nextflow est intelligent concernant l'utilisation des ressources. Donc si vous commencez à mettre ces nombres à des valeurs plus élevées, vous verrez que Nextflow commence à mettre en file d'attente les soumissions de tâches l'une après l'autre, plutôt que de toutes les exécuter en parallèle, afin de ne pas sur-demander les ressources disponibles.

## 2.4. Exécuter le workflow avec la configuration modifiée

Essayons d'exécuter à nouveau un workflow et sauvegardons un nouveau rapport cette fois.

D'accord, nous pouvons télécharger ce fichier et y jeter un œil.

Oui, sans surprise, il ressemble exactement au même parce que c'est un workflow factice, qui ne fait rien de réel. Mais vous pouvez imaginer comment cette approche itérative de définition de limites et de réalisation de workflows réels avec ce type de rapports vous permet de faire une approche basée sur les preuves pour définir une configuration appropriée et vraiment tirer le meilleur parti des ressources de calcul dont vous disposez.

Vous pouvez commencer à être vraiment intelligent à ce sujet. Nextflow a une capacité intégrée pour réessayer les échecs, et vous pouvez en profiter dans votre fichier de configuration en utilisant une closure comme ceci et en définissant dynamiquement les ressources qui sont mises à disposition. Donc ici, j'ai dit à Nextflow de multiplier ces deux gigaoctets par la tentative de réessai. Donc le deuxième réessai obtiendra quatre gigs, le troisième réessai obtiendra six gigs et ainsi de suite. Cela va un peu au-delà de la portée de cette formation, mais si vous êtes intéressé·e, consultez la documentation Nextflow, qui a une belle section sur la logique de réessai dynamique.

## 2.5. Ajouter des limites de ressources

Maintenant, une chose que vous pourriez remarquer à ce sujet est que ce genre de chose peut rendre assez facile de dépasser accidentellement les ressources disponibles sur votre système. Si vous demandez plus de ressources que disponibles, Nextflow lancera une erreur concernant votre configuration et arrêtera l'exécution. Pour éviter cela, vous pouvez utiliser quelque chose appelé limites de ressources.

Sous la portée du processus, dans notre workflow, nous pouvons définir des limites de ressources comme ceci, qui prend un tableau, et nous pouvons spécifier la mémoire maximale, les CPUs et le temps qui sont disponibles sur ce système.

Définir des valeurs élevées ici n'augmente pas la quantité de ressources demandées. Nous allons toujours utiliser un gigaoctet dans nos requêtes, mais cela signifie que si l'une de ces requêtes atteint 750, elle atteindra ce plafond et rien de plus ne sera demandé, ce qui signifie que Nextflow continuera à s'exécuter et ne plantera pas en raison de ressources indisponibles.

C'est donc une bonne protection à utiliser, surtout si vous utilisez une logique dynamique avec votre allocation de ressources.

L'autre situation où c'est vraiment utile est si vous utilisez des pipelines publics et non contrôlés par vous. Ils peuvent venir avec des valeurs par défaut de configuration, et Nextflow prendra automatiquement la bonne approche en plafonnant toutes les demandes de ressources pour s'exécuter sur votre système.

D'accord, super. Nous avons parlé du logiciel. Nous avons parlé de l'allocation de ressources, et nous avons décrit différentes portées de configuration, à la fois pour tous les processus et pour des processus spécifiques.

## 3. Utiliser un fichier de paramètres pour stocker les paramètres du workflow

D'accord, ensuite nous allons tourner notre attention vers les paramètres. Nous pouvons définir des paramètres dans le fichier de configuration tout comme nous l'avons fait auparavant dans le script Nextflow. Donc params.greeting equals hello ou utiliser la portée params et définir foo equals bar.

Et c'est génial pour définir des valeurs par défaut pour votre workflow. Cependant, lorsque vous exécutez des pipelines, il peut être agréable de spécifier les paramètres dans un fichier JSON ou YAML.

Utiliser un fichier comme celui-ci est bien mieux que de spécifier des options en ligne de commande avec --. Car lorsque vous exécutez un workflow, vous devrez peut-être spécifier de nombreux paramètres et il peut être fastidieux de tous les écrire sur une seule CLI et sujet aux erreurs. De plus, il est peu probable que vous vous souveniez de tous les paramètres que vous avez utilisés, donc si vous les codez dans un fichier, il est plus facile de lancer à nouveau le workflow, en utilisant les mêmes paramètres à l'avenir.

Nous avons un exemple de fichier ici appelé test-params, et vous pouvez voir que cela spécifie les trois paramètres que nous avons dans notre workflow avec trois valeurs différentes. Personnellement, je trouve YAML plus facile à écrire que JSON. Donc juste pour démontrer que cela fonctionne, je vais créer un nouveau fichier appelé test.yaml et copier ceux-ci, me débarrasser des guillemets et sauvegarder.

Ces fichiers JSON et YAML peuvent être plus faciles à écrire car leur syntaxe est plus familière. Mais notez que ceux-ci sont uniquement pour les paramètres et ils ne prennent que la syntaxe clé-valeur comme ceci.

## 3.1. Exécuter le workflow en utilisant un fichier de paramètres

Essayons-le. Je fais la même commande qu'avant. Je me débarrasse du rapport et je vais faire -params-file test-params.yaml.

Non, c'est une option Nextflow de base, donc c'est un seul trait d'union.

D'accord. Il a exécuté le workflow et il a utilisé les paramètres dans ce fichier YAML au lieu que je les spécifie tous en ligne de commande. Cela peut sembler excessif juste pour cet exemple simple, mais vous pouvez imaginer que si vous avez 10 ou 20 paramètres différents, cela peut être pénible de les taper manuellement, et c'est juste beaucoup plus facile à éditer dans un éditeur de code et à conserver pour des raisons de reproductibilité.

## 3. Déterminer quel(s) exécuteur(s) doit/doivent être utilisé(s) pour effectuer le travail

D'accord. Nous avons parlé d'empaquetage logiciel avec Docker et conda. Nous avons parlé des exigences de ressources de processus avec les CPUs et la mémoire. Et nous avons parlé un peu de la façon de spécifier les paramètres lors de l'exécution des workflows.

Les dernières parties de la configuration concernent vraiment l'exécution, l'infrastructure de calcul sous-jacente elle-même, et c'est le vrai joyau de la couronne de Nextflow : que nous pouvons exécuter ces mêmes workflows sur plusieurs infrastructures de calcul différentes.

Je vais effectivement basculer vers le matériel de formation écrit pendant une seconde. Dans cette partie de la formation, nous pouvons voir quelques exemples différents de la façon dont différents exécuteurs, dans ce cas, des ordonnanceurs HPC, définissent les exigences de ressources nécessaires pour soumettre un travail.

Donc pour Slurm, vous avez ces en-têtes SBATCH, qui définissent --mem et le nombre de CPU. Si vous utilisez PBS, vous avez des en-têtes différents, et si vous utilisez Grid Engine, vous avez des en-têtes différents encore.

Vous pouvez imaginer que c'est encore plus différent si vous voulez exécuter sur le cloud, que ce soit AWS Batch, Google Cloud, Azure, ou plus.

Chacune de ces infrastructures de calcul sous-jacentes est appelée un exécuteur et Nextflow sait comment parler à tous ces différents exécuteurs afin de soumettre des tâches avec la syntaxe correcte.

La bonne nouvelle est que vous n'avez pas besoin de connaître cela. Tout ce que vous avez à faire est de dire à Nextflow quel exécuteur utiliser.

## 3.1. Cibler un backend différent

Nous retournons à notre fichier de configuration et au processus, nous faisons executor, et je vais taper local.

Local est en fait la valeur par défaut, si vous ne spécifiez aucun autre exécuteur, local est ce qui sera utilisé, et cela signifie simplement votre système hôte, où que vous ayez lancé Nextflow.

Je pourrais spécifier à la place, slurm. Et cela soumettrait des tâches Slurm, ou je pourrais dire awsbatch, et cela soumettrait des tâches à AWS Batch.

Vous avez besoin de configuration supplémentaire dans certains cas, par exemple, l'exécution sur le cloud nécessitera certaines informations d'identification, mais vraiment c'est le cœur de cela, et cela peut être aussi simple qu'une ou deux lignes de configuration pour exécuter votre workflow dans un environnement de calcul complètement différent.

Même si nous exécutons sur un système simple dans Codespaces, je peux toujours jouer avec cela un peu et prétendre que nous exécutons sur Slurm. Si je lance ensuite à nouveau le workflow, nextflow run hello-config, cela échouera car il ne pourra pas soumettre de tâches à Slurm. Mais nous pouvons toujours aller dans les répertoires de travail et voir ce que Nextflow a fait. Donc si nous allons dans ce répertoire de travail et regardons .command.run, vous pouvez voir en haut de ce fichier, nous avons maintenant ces lignes d'en-tête sbatch, qui ont essayé de spécifier les ressources nécessaires pour la tâche Slurm.

## 4. Utiliser des profils pour sélectionner des configurations prédéfinies

D'accord, nous y sommes presque. La dernière partie de ce chapitre concerne les profils de configuration. Si vous exécutez votre pipeline sur plusieurs systèmes différents, il pourrait être ennuyeux d'avoir tous ces fichiers de configuration Nextflow différents, que vous devez spécifier à chaque fois.

Au lieu de cela, vous pouvez encoder des groupes de configuration dans votre fichier de configuration Nextflow, et activer ou désactiver ces groupes en utilisant un indicateur de profil. Voyons à quoi cela ressemble.

## 4.1. Créer des profils pour basculer entre le développement local et l'exécution sur HPC

Nous allons créer deux profils dans notre exemple ici, un pour mon ordinateur portable et un pour un système HPC plus lourd. Je vais tricher un peu et simplement copier le code du matériel de formation et le mettre ici.

Nous avons une nouvelle portée appelée profiles, puis nous avons un nom pour chaque profil, qui peut être n'importe quoi. Et à l'intérieur de cela, nous avons une configuration, qui ressemble exactement à la configuration de niveau supérieur que nous avons déjà écrite. Donc encore, nous avons la portée process, la portée docker.

Sur le profil appelé my_laptop, je dis d'exécuter en utilisant l'exécuteur local, donc sur mon système hôte et d'utiliser Docker.

Sur le profil university_hpc ici, je dis d'utiliser Slurm pour soumettre des tâches, d'utiliser conda au lieu de Docker, et je spécifie différentes limites de ressources, qui peuvent correspondre à la taille du système de nœuds sur le HPC que j'utilise.

Par défaut, aucune de cette configuration ne sera utilisée lorsque j'exécute Nextflow, je dois spécifier que je veux utiliser l'un de ces profils.

## 4.2. Exécuter le workflow avec un profil

Faisons nextflow run hello-config. Et je vais faire -profile, trait d'union simple parce que c'est une option Nextflow de base. Et puis le nom que je lui ai donné, qui est my_laptop. Nextflow devrait maintenant utiliser le bloc de configuration qui a été spécifié dans ce profil de configuration, et l'appliquer lorsqu'il exécute Nextflow. Si je voulais utiliser l'autre bloc de configuration, je n'ai qu'à changer ce nom de profil. Beaucoup plus facile à retenir. Beaucoup plus facile à utiliser.

## 4.3. Créer un profil de test

Notez, les profils peuvent avoir n'importe quel type de configuration, donc cela ne doit pas être lié à votre environnement d'exécution. Par exemple, créons un nouveau profil ici, qui a un ensemble de paramètres. Nous pouvons changer cela en tux et changer en my profile, et maintenant lorsque nous faisons profile test, cela va spécifier ces paramètres, qui écraseront les paramètres qui sont spécifiés au niveau supérieur du workflow.

Lorsque vous exécutez Nextflow, vous pouvez enchaîner plusieurs profils et ils seront appliqués en séquence.

## 4.4. Exécuter le workflow localement avec le profil de test

Donc je peux prendre la commande précédente et faire virgule test. Cela appliquera d'abord la configuration my_laptop, puis il appliquera la configuration test. S'il y a un chevauchement, alors le profil à droite écrasera toute configuration dans les profils précédents. Si j'appuie sur entrée, voyons ce qui se passe.

D'accord, nous avons un nouveau fichier de résultats ici. Vous pouvez voir le My Profile, que j'ai spécifié comme l'une des options. Et nous pouvons également voir cowpy, my profile, et effectivement, il y a tux. Donc cela a fonctionné.

## Récapitulatif

D'accord ! Incroyable. C'est tout. Vous avez réussi jusqu'à la fin du cours. Vous obtenez un peu de confettis de célébration. Bravo d'avoir terminé ce chapitre.

[Transcription vidéo suivante :octicons-arrow-right-24:](07_next_steps.md)
