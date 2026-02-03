# Partie 1 : Hello World - Transcription

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Traduction assistée par IA - [en savoir plus et suggérer des améliorations](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

<div class="video-wrapper">
  <iframe width="560" height="315" src="https://www.youtube.com/embed/8X2hHI-9vms?si=F0t9LFYLjAWoyRXj&amp;list=PLPZ8WHdZGxmXiHf8B26oB_fTfoKQdhlik" title="YouTube video player" frameborder="0" allow="accelerometer; autoplay; clipboard-write; encrypted-media; gyroscope; picture-in-picture; web-share" referrerpolicy="strict-origin-when-cross-origin" allowfullscreen></iframe>
</div>

!!!note "Notes importantes"

    Cette page montre uniquement la transcription. Pour des instructions complètes étape par étape, retournez au [matériel de formation](../01_hello_world.md).

    Les numéros de section affichés dans la transcription sont fournis à titre indicatif uniquement et peuvent ne pas inclure tous les numéros de section du matériel.

## Bienvenue

Bonjour, bienvenue dans le Chapitre Un de Hello Nextflow.

Dans cette première partie d'un cours en six parties, nous allons voir les bases de Nextflow. Nous allons commencer par exécuter quelques commandes dans un terminal, puis nous prendrons ces commandes Bash et verrons comment les intégrer dans un script Nextflow.

Nous essaierons d'exécuter ce premier pipeline Nextflow, verrons ce que Nextflow fait, où il s'exécute, quels fichiers il crée, et quel est le but de ces fichiers.

Très bien, commençons.

## training.nextflow.io

Tout d'abord, allez sur training.nextflow.io. Comme précédemment, tout le matériel est écrit ici, et je vais le parcourir étape par étape. Je vais montrer mon écran pendant que je fais les étapes de la formation, mais tout ce que je dis se trouve dans le matériel de formation, donc vous pouvez le suivre à votre propre rythme, et vous pouvez tout trouver écrit là.

Cette vidéo a également des sous-titres activés, alors n'hésitez pas à les afficher et à suivre exactement ce que je dis au fur et à mesure.

D'accord, allons à Hello Nextflow. C'est le cours que nous allons faire aujourd'hui, et nous avons déjà fait l'orientation dans la première vidéo, donc nous allons directement passer à la partie un. Hello World.

Très bien, je vais maintenant quitter ce matériel de formation et passer à mon environnement Code Spaces. C'est ce que nous avons configuré dans la première vidéo. J'espère que vous avez quelque chose de très similaire à cela dans votre propre système. J'utilise VS Code et je regarde le matériel de formation et j'ai changé de répertoire dans le répertoire hello Nextflow.

## 0. Échauffement : Exécuter Hello World directement

D'accord. Commençons par quelques bases, qui devraient être familières à tout le monde. Je vais commencer par écrire une commande très basique dans le terminal. Ici en bas, je vais dire 'echo Hello World!"' appuyer sur entrée et, sans surprise, le terminal fait ce que je lui demande et renvoie cette chaîne. Hello world.

D'accord, ensuite je vais appuyer sur la flèche du haut pour récupérer cette commande et la modifier un peu. Cette fois, redirigeons cette sortie vers un fichier. Je vais plutôt l'écrire dans output.txt et appuyer sur entrée, rien dans le terminal cette fois parce que la sortie n'est pas allée au terminal. Elle est allée dans ce fichier.

Je peux ensuite lire ce fichier en faisant 'cat output.txt', appuyez sur tab ici pour auto-compléter le nom du fichier et voilà. Le fichier est là.

Je peux également voir ce fichier dans la barre latérale dans l'explorateur de fichiers dans VS Code. Je peux double-cliquer dessus et l'ouvrir ici. Si vous voulez l'ouvrir dans VS Code sans rien cliquer, vous pouvez également faire "code" puis "output.txt" et cela fait la même chose.

Parfait. C'est la première étape. Très simple.

## 1. Examiner le script de démarrage du workflow Hello World

D'accord. Nous allons maintenant faire exactement la même chose, mais dans Nextflow, au lieu de directement dans le terminal.

Nous allons utiliser le premier exemple de script pour commencer, ce fichier s'appelle Hello World. Je peux faire "ls" pour le voir dans un terminal, et je suis sur Mac, donc je peux faire commande + clic pour ouvrir ce fichier, ou j'aurais pu simplement double-cliquer dans la barre latérale ici.

Il y a quelques choses que nous pouvons voir dans ce fichier. Tout en haut, il y a une déclaration avec un dièse disant que c'est un fichier Nextflow et qu'il peut être exécuté. Il y a quelques commentaires ici, juste des commentaires de code réguliers en gris clair, qui n'affectent pas l'exécution, et aident simplement à lire le script.

Et ensuite il y a deux structures principales. Il y a un process ici et un workflow.

Les processes dans Nextflow sont les étapes du pipeline. Ce sont les parties qui font réellement la logique et effectuent le traitement.

Le workflow ensuite en bas relie ces processes ensemble et gouverne la logique du workflow, comment tout se connecte les uns aux autres.

Nous allons commencer par regarder un process. Nous reviendrons au workflow dans un moment.

## 1.2 La définition du process

Donc chaque process commence par un mot-clé process. Il a un nom et puis il a des accolades et tout ce qui se trouve dans ces accolades est ce process unique.

Un process doit avoir une section script, et contenu ici est un extrait bash dans une chaîne multi-lignes, qui est la partie du code qui est réellement exécutée dans l'environnement de calcul.

Nous avons également une déclaration output ici, qui dit à Nextflow quels fichiers sont censés être créés par le script. Notez que l'output ici a un mot-clé path, qui dit à Nextflow que c'est un fichier, pas une valeur, ou une chaîne.

Dans le bloc script, c'est juste une déclaration bash régulière, et c'est exactement la même que ce que nous avons écrit dans le terminal. Nous faisons echo de hello world vers un fichier appelé output.txt. Ce output.txt est ensuite récupéré par la définition output. La définition output ne fait rien en réalité. Elle dit juste à Nextflow quoi s'attendre, et si ce fichier n'était pas créé, Nextflow lancerait une erreur.

Notez que cet exemple n'est pas excellent parce que nous avons codé en dur le nom du fichier ici, output.txt et output.txt. Si l'un de ceux-ci était changé, cela causerait une erreur dans notre workflow.

Il y a une meilleure façon de faire cela avec des variables, que nous couvrirons dans une minute.

## 1.3 La définition du workflow

D'accord. En descendant au workflow, nous pouvons voir que nous avons un commentaire et ensuite nous exécutons le process appelé sayHello. C'est le même mot-clé qui est ici en haut. C'est à peu près aussi simple qu'un workflow peut être. Nous appelons juste un seul process sans entrée variable, donc nous ne le connectons à rien d'autre. Dans la partie suivante de ce cours, nous parlerons de comment rendre cela plus puissant en utilisant des entrées variables et en connectant des choses avec des channels.

## 2. Exécuter le workflow

D'accord, c'est tout ce dont nous avons besoin. Voyons si nous pouvons l'exécuter et voir ce qui se passe. Je vais juste effacer le terminal et ensuite je vais faire "nextflow run", et je vais appeler le nom du fichier, qui est hello-world.nf. C'est tout ce dont nous avons besoin pour exécuter un pipeline Nextflow. Ce pipeline ne prend pas d'entrée, donc nous n'avons pas besoin d'autres arguments.

Appuyons sur entrée et voyons ce qui se passe.

D'accord. J'espère que vous devriez avoir une sortie qui ressemble à ceci. Nous avons quelques informations nous disant que Nextflow s'est exécuté et quelle version il utilisait. Il nous dit quel script a été lancé et il nous donne un nom généré aléatoirement pour cette exécution de workflow particulière. Dans ce cas, le mien s'appelait "gloomy_crick".

La partie la plus importante de ceci cependant, c'est qu'il nous dit quelles étapes se sont exécutées dans le pipeline. Vous pouvez voir que notre process appelé sayHello s'est exécuté, et il s'est exécuté une fois et il était complet à cent pour cent.

Cette partie ici est le hash pour cette tâche de workflow particulière. Chaque process s'exécute une ou plusieurs fois, et chacune de ces exécutions est appelée une tâche.

## 2.2. Trouver la sortie et les logs dans le répertoire work

Chaque tâche obtient son propre répertoire isolé où elle s'exécute, donc elle est séparée du reste de l'exécution du workflow. Ce hash correspond à la structure de fichiers dans le répertoire work. Si je fais "tree work", nous pouvons voir a0, puis une version plus longue d'un hash court, puis notre fichier output.txt. Vous pouvez également le voir dans la barre latérale.

Vous pouvez voir dans la barre latérale qu'il y a quelques fichiers supplémentaires ici. La raison pour laquelle ceux-ci n'apparaissent pas dans le terminal est qu'ils sont des fichiers cachés, ils commencent par un point. Et en effet, si je fais "tree -a" pour all, et "work", nous pouvons les voir ici.

Ces fichiers dot sont présents dans chaque répertoire work que Nextflow crée, et chacun a une tâche légèrement différente. Premièrement .command.begin inclut juste quelques instructions pour Nextflow qui configure la tâche avant qu'elle ne s'exécute. .command.run sont les instructions réelles exécutées par Nextflow lui-même. Ensuite .command.sh est probablement le plus intéressant. C'est le script qui a été résolu à partir de notre bloc script du process.

Si je l'ouvre, vous pouvez voir que nous avons notre "echo Hello World" vers le fichier output.txt. C'est exactement la même chose que notre process dans ce cas, mais si nous avons des variables dans notre code Nextflow, chaque tâche aura un .command.sh différent, et vous pouvez voir comment ces variables ont été résolues.

Les autres fichiers ont à voir avec comment la tâche s'est exécutée. Donc .command.err, .log et .out sont l'erreur standard, la sortie standard et les deux combinés. Et .exitcode dit à Nextflow comment cette tâche s'est exécutée avec quel code de sortie, si elle a réussi ou non.

Enfin, nous avons notre fichier output.txt et effectivement, "Hello World", c'est ce que nous attendions et c'est ce qui a été créé.

D'accord, parfait. C'était votre toute première exécution Nextflow. Félicitations. C'est vraiment aussi simple que cela.

Ensuite, nous allons voir comment faire cela un peu plus facilement pour ne pas avoir à éditer le code chaque fois que nous voulons faire un changement dans la façon dont le pipeline s'exécute.

## 3. Gérer les exécutions de workflow

Cette structure de répertoires est excellente pour garder toutes les tâches séparées et tout organisé, mais bien sûr, ce n'est pas très pratique pour trouver vos fichiers de sortie. Vous ne voulez pas fouiller dans des tas de répertoires imbriqués pour essayer de trouver les résultats de votre pipeline.

## 3.1. Publier les sorties

La bonne nouvelle est que vous n'êtes pas censé le faire. Les répertoires work sont vraiment juste pour que Nextflow les utilise lui-même. Donc ce que nous allons faire, c'est que nous allons utiliser une fonction de Nextflow appelée "publishDir".

Nous retournons à notre workflow, allons au process. Nous pouvons ajouter une nouvelle déclaration ici appelée une directive. C'est ainsi que Nextflow appelle ces choses en haut des processes qui augmentent le fonctionnement, et celle que nous allons utiliser s'appelle publishDir.

Vous pouvez voir que j'ai commencé à taper ici et l'extension Nextflow pour VS Code m'a suggéré la directive, donc je peux juste appuyer sur entrée.

D'accord, je vais suivre cela avec un répertoire appelé "results" et nous allons lui dire de copier les fichiers de sortie là. Donc je vais dire mode copy. Parfait. Je vais sauvegarder et exécutons le workflow à nouveau.

nextflow run hello-world.nf

Il s'exécute exactement de la même manière. Bien que notez que nous avons un hash légèrement différent cette fois. Nextflow utilisera un hash différent à chaque fois que vous exécutez le workflow. Et nous avons un ensemble différent de répertoires work en conséquence. Les zones, l'une s'appelle EB à la place, mais vous pouvez voir que tous les fichiers sont les mêmes. Cependant, ce qui est nouveau cette fois, c'est que nous avons également un répertoire appelé "results".

Dans "results" ici, nous avons notre fichier de sortie. C'est ce que nous avons dit à Nextflow de faire. Nous avons dit, enregistrez les fichiers de résultats dans un répertoire appelé "results" et copiez-les là. Et donc c'est maintenant beaucoup plus facile à trouver. C'est juste là à côté de l'endroit où nous avons lancé le workflow et tous les différents fichiers peuvent y être organisés comme nous le souhaitons, indépendamment de où ou comment Nextflow a exécuté l'exécution réelle.

Notez que publishDir peut gérer les liens symboliques, ce qui est bien si vous travaillez sur un système de fichiers partagé et que vous voulez économiser de l'espace. Et aussi vous n'avez pas à définir tous les fichiers qui sont créés par un process comme une sortie.

Nextflow copiera uniquement les choses qui sont définies dans ce bloc output. Donc si vous avez des fichiers intermédiaires créés par l'étape, qui ne sont pas nécessaires en aval de ce process, vous ne les définissez simplement pas dans output et ils n'apparaîtront pas dans le publishDir. Donc c'est une façon de garder vos fichiers de sortie d'un pipeline propres et de supprimer facilement les fichiers intermédiaires une fois le travail terminé.

Une petite note ici. Il y a une nouvelle syntaxe Nextflow qui arrive appelée workflow output definitions, qui remplacera éventuellement publishDir. Cela nous donne un moyen de définir toutes les sorties d'un workflow au niveau du pipeline dans le bloc workflow. Ceci est décrit dans les docs Nextflow si vous voulez l'essayer. Mais pour l'instant, publishDir sera encore là pendant un moment, donc nous l'avons toujours dans la formation pour 2025.

## 3.2. Relancer un workflow avec -resume

D'accord. J'ai mentionné que le répertoire work ici a maintenant deux ensembles de résultats avec un hash différent de chaque fois que nous exécutons le workflow. C'est bien. Cependant, parfois nous ne voulons pas recalculer les étapes à chaque fois si nous n'en avons pas besoin.

Peut-être que vous êtes en train de construire itérativement votre workflow et vous ajoutez des étapes et vous voulez que les premières étapes réutilisent simplement les versions mises en cache. Ou peut-être que quelque chose s'est mal passé sur votre système de calcul à mi-chemin de votre workflow et vous voulez qu'il continue à partir de là où il s'est arrêté, mais saute les étapes qu'il avait déjà terminées.

Nextflow a une fonctionnalité intégrée pour cela appelée resume. Essayons-le. Donc tout d'abord, je vais juste regarder le répertoire work pour que nous puissions nous rappeler ce qui était là.

Et ensuite je vais faire "nextflow run hello-world.nf" et je vais ajouter une seule commande ici, "-resume".

Notez, un seul tiret, c'est vraiment important. Je vais l'exécuter et la sortie va ressembler essentiellement exactement à la même chose, avec quelques petites différences.

Notez ici qu'il dit "cached" en gris. Cela signifie que Nextflow n'a pas exécuté la tâche. Cette fois, il a trouvé quelque chose qui correspondait à ce que nous demandions et il a réutilisé ces sorties directement plutôt que de réexécuter l'étape.

Et effectivement, si vous regardez le hash ici, vous pouvez voir que cela correspond au hash existant que nous avions d'une exécution précédente.

## 3.3. Supprimer les anciens répertoires work

D'accord. Mais si vous développez de manière itérative, vous allez accumuler beaucoup de ces fichiers de workflow. Cela peut être un problème si vous manquez d'espace.

Nextflow peut nous aider à nettoyer ces répertoires work avec quelques commandes d'aide. Si je fais "nextflow log". Cela me donnera une liste de toutes les différentes exécutions de workflow que j'ai faites dans ce répertoire, et elles ont les noms d'exécution ici. Vous pouvez voir le gloomy quick, qui était le premier que nous avons exécuté, puis ces deux nouveaux.

Nous pouvons maintenant prendre ce nom et les utiliser avec la commande "nextflow clean". Je peux spécifier un seul nom d'exécution. Ou encore mieux, je peux dire à Nextflow de supprimer tout ce qui se trouve avant un seul nom de workflow avec "-before", et je vais mettre "stupefied_shaw". C'était ma dernière exécution, "-n".

La commande "-n" a dit à Nextflow de le faire comme un dry run sans réellement supprimer quoi que ce soit pour de vrai, et il nous dit lesquels des répertoires hash auraient été supprimés. Effectivement, c'est juste celui de la première exécution. Les deux secondes exécutions utilisent le même répertoire hash.

Je vais l'exécuter à nouveau, mais maintenant au lieu de "-n" pour dry run, je vais faire "-f" pour force et il a supprimé ce répertoire hash. Maintenant si je fais "tree work", nous pouvons voir, nous avons juste ce fichier de sortie restant.

Parfait. Donc nous avons réussi à nettoyer beaucoup d'espace disque là.

Quelques points à noter lors de la suppression de répertoires work, si vous faites des liens symboliques vers votre répertoire results, ces sources de liens symboliques seront maintenant supprimées et vos résultats seront perdus pour toujours. C'est pourquoi utiliser le mode copy est une chose plus sûre à faire, et généralement ce que nous recommandons.

Deuxièmement, la fonctionnalité resume de Nextflow repose sur ces répertoires work. Donc si vous les supprimez et que vous exécutez Nextflow à nouveau, la fonctionnalité resume ne fonctionnera plus. Donc c'est à vous de garder une trace de quelles choses vous pourriez avoir besoin ou non, et de supprimer uniquement les choses lorsque vous êtes sûr qu'il est sûr de le faire.

L'autre chose que nous pouvons faire est que nous pouvons simplement supprimer l'ensemble du répertoire work si nous avons terminé notre exécution de workflow et que nous sommes sûrs que nous n'en avons plus besoin.

Donc je peux faire "rm -r work". Je sais qu'il n'y avait rien d'important là-dedans. J'ai mes résultats qui m'intéressent dans le répertoire results où nous les avons copiés. Et donc il était sûr de supprimer le répertoire work. C'est à vous de décider laquelle de ces approches vous utilisez.

## 4. Utiliser une entrée variable passée sur la ligne de commande

D'accord, et ensuite ? J'ai mentionné que nous avions codé en dur certaines des valeurs dans notre script de workflow ici, le fichier output.txt, et qu'il pourrait y avoir une meilleure façon de faire cela.

Commençons par cela. Ce que nous allons faire, c'est trois choses. Nous allons ajouter une nouvelle entrée au process. Nous allons dire au script du process comment utiliser cette entrée, puis nous allons le câbler dans le workflow afin que nous puissions l'utiliser dynamiquement avec un flag de ligne de commande lors de l'exécution de Nextflow.

Donc tout d'abord. Ajoutons un bloc input ici. Tout comme output. C'est une nouvelle section pour le process, et je vais dire, "val greeting".

Notez ici, je dis "val", ce qui dit que c'est une variable, pas un path.

Je peux ensuite descendre dans le script et ensuite je peux retirer ce texte codé en dur ici et faire $greeting. Cela fonctionne comme n'importe quel autre langage de programmation. Nous définissons une variable ici et nous la référençons dans ce bloc script. Lorsque Nextflow exécute ce process, la variable sera interpolée. Et lorsque nous irons regarder ce fichier .command.sh, nous verrons la chaîne réelle codée en dur ici à la place.

## 4.1.3. Configurer un paramètre CLI et le fournir comme entrée à l'appel du process

D'accord, mais où fournissons-nous la variable ? Ensuite, nous descendons à la section workflow, et vous pouvez voir que l'extension ici dit que nous attendons maintenant une entrée, et elle m'a donné un avertissement.

Maintenant, la chose la plus simple que nous pourrions faire est simplement de la coder en dur. Je pourrais écrire "Hello World" et fournir cette entrée de chaîne au process. Mais encore une fois, cela ne résoudrait pas vraiment de problèmes. Nous devrions toujours retourner et éditer le code du pipeline chaque fois que nous voudrions changer quelque chose, ce qui n'est pas bon.

La bonne nouvelle est que Nextflow a un système intégré pour gérer les arguments de ligne de commande appelé paramètres. Donc à la place, je peux utiliser l'une de ces variables spéciales appelées params et je peux l'appeler comme je veux, mais je vais dire greeting pour qu'elle corresponde à la logique du workflow.

Sauvegardons et voyons ce que nous pouvons faire avec cela.

Donc si je retourne au terminal. Donc nous faisons "nextflow run hello-world.nf". Juste comme avant, mais la différence clé est que nous faisons --greeting

Notez, il y a deux tirets ici parce que c'est un paramètre. Quand nous avons repris le workflow avant, c'était un seul tiret. C'est parce que resume est une option Nextflow de base, et ceci est un paramètre qui est spécifique à notre pipeline.

Ne confondez pas les deux. C'est facile de faire cela. Si vous faisiez --resume au lieu d'un seul tiret, alors ce serait "params.resume", ce qui ne ferait rien. De même, si vous faisiez un seul tiret ici, Nextflow ne le reconnaîtrait pas comme un argument clé.

Donc c'est --greeting, qui correspond à params.greeting.

Je peux maintenant suivre cela avec n'importe quel texte que je veux. Donc je suis en Suède en ce moment, donc je vais dire, "Hej världen".

Donc exécutons-le, voyons ce qui se passe, moment de vérité.

D'accord, donc vous pouvez voir que le process s'est à nouveau exécuté, juste comme avant, sayHello avec une seule exécution.

Cela aura écrasé le fichier qui était dans le répertoire publishDir "results". Donc soyez prudent lorsque vous réexécutez les fichiers car les choses dans le répertoire publié seront écrasées.

Je peux maintenant faire "code results/output.txt", et effectivement, notre sortie a été mise à jour et dit maintenant "Hej världen".

## 4.2. Utiliser des valeurs par défaut pour les paramètres de ligne de commande

D'accord, c'est super. Mais le problème maintenant est que notre workflow repose sur nous définissant toujours ce paramètre, et c'est bien d'avoir des valeurs par défaut sensées pour que les choses s'exécutent de manière sensée pour votre workflow à moins que vous ne remplaciez les valeurs par défaut.

Donc la façon dont nous faisons cela est en définissant une valeur par défaut pour le paramètre dans notre script de workflow.

Donc si je retourne à mon fichier hello-world.nf, je peux aller dans le script juste au-dessus du workflow, taper "params.greeting" et le définir comme n'importe quelle autre variable. Donc mettons une chaîne ici et disons "Holà mundo!"

Maintenant ce paramètre a une valeur par défaut définie, qui sera utilisée ici, ou nous pouvons toujours la remplacer sur la ligne de commande avec --greeting, exactement comme nous l'avons fait avant.

Donc vérifions que ça fonctionne. "nextflow run hello-world.nf"

Pas d'arguments de ligne de commande cette fois, et vérifions si cela a fait la bonne chose.

"code results/output.txt". Et voilà. Nous avons notre valeur par défaut.

D'accord, essayons à nouveau, vérifions juste que je ne vous raconte pas de mensonges. Exécutons-le à nouveau, mais faisons --greeting, et utilisons l'exemple du matériel de formation, disons "Konnichiwa!"

Réexécute le workflow, et effectivement, notre fichier de sortie en haut vient d'être mis à jour avec la nouvelle valeur que nous avons fournie sur la ligne de commande.

Parfait. C'est un aspect vraiment central pour écrire n'importe quel workflow Nextflow. Définir des valeurs par défaut sensées dans votre code de pipeline, mais le rendre très facile à configurer pour l'utilisateur final en ayant des arguments de ligne de commande sur le terminal.

Notez que l'utilisateur final peut remplacer la configuration à plusieurs endroits différents. Vous pouvez avoir un fichier de configuration dans votre répertoire personnel, qui est appliqué à chaque exécution Nextflow que vous faites. Vous pouvez avoir un fichier de configuration dans un répertoire de lancement. Vous pouvez avoir un fichier de configuration dans un répertoire de pipeline. Tous ces différents emplacements de configuration sont chargés dans un ordre spécifique, qui est décrit dans les docs Nextflow.

D'accord, c'est la fin de la section un. Nous avons eu notre tout premier script de workflow dans Nextflow avec un process et un workflow. Nous avons examiné les entrées, les sorties, les scripts et la publication, et comment câbler les paramètres et un canal d'entrée dans notre process.

Félicitations, votre premier pas vers l'écriture de code Nextflow est terminé.

Faites une petite pause et je vous verrai de retour dans quelques minutes pour le chapitre deux.

[Transcription de la vidéo suivante :octicons-arrow-right-24:](02_hello_channels.md)
