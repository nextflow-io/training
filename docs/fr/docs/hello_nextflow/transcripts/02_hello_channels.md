# Partie 2 : Hello Channels - Transcription

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Traduction assistée par IA - [en savoir plus et suggérer des améliorations](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

<div class="video-wrapper">
  <iframe width="560" height="315" src="https://www.youtube.com/embed/lJ41WMMm44M?si=xCItHLiOQWqoqBB9&amp;list=PLPZ8WHdZGxmXiHf8B26oB_fTfoKQdhlik" title="YouTube video player" frameborder="0" allow="accelerometer; autoplay; clipboard-write; encrypted-media; gyroscope; picture-in-picture; web-share" referrerpolicy="strict-origin-when-cross-origin" allowfullscreen></iframe>
</div>

!!!note "Notes importantes"

    Cette page affiche uniquement la transcription. Pour des instructions complètes étape par étape, retournez au [matériel de formation](../02_hello_channels.md).

    Les numéros de section affichés dans la transcription sont fournis à titre indicatif uniquement et peuvent ne pas inclure tous les numéros de section dans les matériels.

## Bienvenue

Bonjour, bienvenue dans la partie deux de Hello Nextflow.

Ce chapitre s'appelle Hello Channels. Nous allons parler de cette partie fondamentale de Nextflow.

Les channels sont les éléments qui connectent les différentes étapes de votre pipeline, la façon dont vos données et votre logique circulent à travers votre workflow.

D'accord, plongeons-nous dedans.

Commençons par aller sur training.nextflow.io

Hello Nextflow dans la barre latérale et cliquons sur la partie deux, Hello Channels.

Tout le matériel est écrit ici, vous pouvez donc suivre à votre propre rythme et rattraper tout ce que vous auriez pu manquer.

Une fois que vous avez ouvert le site web, vous pouvez charger Codespaces et nous continuerons là où nous en étions à la fin du dernier chapitre.

## 0. Échauffement : Exécuter hello-channels.nf

Pour ce chapitre, nous allons éditer un fichier différent. Celui-ci s'appelle Hello Channels, vous pouvez donc le trouver dans la barre latérale, double-cliquez dessus pour l'ouvrir.

Maintenant, si vous venez juste du chapitre un, ce fichier vous semblera très familier. Le point de départ ici est essentiellement là où nous avons terminé le chapitre un, avec notre process appelé sayHello, notre entrée, sortie, notre publishDir et notre params.greeting, et notre workflow simple.

Nous commençons avec un nouveau fichier, c'est donc un terrain de jeu équitable pour tout le monde, mais vous pouvez continuer avec votre fichier précédent si vous préférez.

Notez que j'ai également supprimé tous les fichiers .nextflow\* et les répertoires work ici, juste pour avoir un point de départ propre. Peu importe si vous le faites ou non, c'est à vous de décider.

D'accord. Commençons par vérifier que ce pipeline fonctionne toujours comme nous l'attendons. Je vais ouvrir le terminal ici.

Faites "nextflow run hello-channels.nf" et appuyez sur entrée.

Il va exécuter ce petit workflow, exécute notre étape sayHello, génère un répertoire work avec ce hash, et voici notre répertoire results et voilà notre fichier de sortie, exactement comme nous l'attendions de notre params.greeting par défaut.

C'est donc excellent. Exactement pareil que le chapitre un, fonctionnant comme nous l'attendons.

## 1. Fournir des entrées variables via un channel explicitement

Dans le chapitre un, vous utilisiez déjà des channels, vous ne le réalisiez simplement pas. Lorsque nous avons spécifié une chaîne de caractères ici, Nextflow a automatiquement créé un channel autour de cette chaîne pour nous, simplement parce qu'il savait que nous appelions un process, donc nous avions besoin d'un channel d'entrée.

La première chose que nous allons faire est de le rendre explicite en écrivant réellement le channel lui-même.

## 1.1. Créer un channel d'entrée

Je vais donc aller dans le workflow ici en bas du script, et je vais dire greeting_ch. C'est une convention que nous utilisons souvent dans le code Nextflow d'avoir un underscore ch à la fin d'un nom de variable lorsqu'il s'agit d'un channel, juste pour qu'il soit facile d'identifier que c'est un channel, mais vous n'êtes pas obligé de le faire. Égale channel of Hello Channels.

Ce que nous venons d'utiliser s'appelle une "Channel Factory" dans le langage Nextflow. C'est cette chose ici, nous définissons cette variable à un nouveau channel, et cette channel factory ici crée un channel pour nous d'une manière particulière.

Il existe une poignée de channel factories différentes que Nextflow possède, pour créer des channels à partir de différents types d'entrées. Dot of est la plus simple, et prend simplement toutes les chaînes de caractères que nous lui donnons.

Notez que lorsque je survole ces mots dans VS Code, l'extension Nextflow me donne une fenêtre contextuelle expliquant ce que fait cette syntaxe, et il y a également un texte lire la suite en bas de cette fenêtre contextuelle.

Si je clique dessus, cela ouvrira la documentation Nextflow dans un nouvel onglet et m'amènera directement à la documentation pour cette chose spécifique. Dans ce cas pour channel.of.

## 1.2. Ajouter le channel comme entrée à l'appel du process

Notez que l'extension nous donne également un avertissement, disant que nous avons créé un nouveau channel ici, mais qu'il n'est utilisé par rien.

Alors, réglons ça. Je vais prendre le nouveau nom de channel et je vais remplacer ce params.greeting par notre nouveau channel.

Notez que nous n'utilisons plus le drapeau de ligne de commande --greeting maintenant, params.greeting n'est plus utilisé, nous revenons à coder en dur cette chaîne. Ce n'est pas grave. J'essaie juste de garder les choses simples. Nous reviendrons plus tard et utiliserons les params à nouveau.

## 1.3. Exécuter à nouveau la commande workflow

D'accord, vérifions juste que cela fonctionne. Ouvrons le terminal et notons à nouveau. Nextflow run hello channels. Vérifiez output.txt, et le voilà.

Excellent, un peu un exemple ennuyeux, faisant exactement la même chose que nous avons fait avant, mais maintenant au moins la logique est un peu plus claire. Nous sommes explicites sur l'écriture d'un nouveau channel.

Nous avons effectivement juste écrit plus de code pour faire la même chose. Mais cela commencera à avoir plus de sens alors que nous devenons un peu plus compliqués avec la façon dont nous créons nos channels.

## 2. Modifier le workflow pour s'exécuter sur plusieurs valeurs d'entrée

D'accord, rendons cela un peu plus intéressant. Il est très rare que vous souhaitiez exécuter un pipeline Nextflow sur une seule entrée, alors donnons-lui plusieurs entrées.

## 2.1. Charger plusieurs salutations dans le channel d'entrée

Depuis la documentation ici, je vais copier ces différentes chaînes de caractères, trois d'entre elles. Hello, Bonjour, Olà. Oh, j'espère que Copilot en suggère quelques autres. Alors appuyons sur tab pour entrer celles-ci.

La documentation Nextflow ici nous dit que nous pouvons donner plusieurs valeurs à cet opérateur, donc cela devrait fonctionner, mais essayons-le et voyons ce qui se passe.

## 2.1.2. Exécuter la commande et regarder la sortie du journal

Eh bien. Oui et non. Voyons voir. Il dit que cinq des cinq tâches ont été exécutées ici, mais il ne nous montre qu'un seul hash, ce qui est un peu étrange. Ce n'est pas grave. Tout est comme prévu ici. Par défaut, Nextflow utilise un type spécial de sortie vers le terminal appelé codes de contrôle ANSI, ce qui signifie qu'il écrase certaines lignes pour donner une vue compressée agréable de tous les différents processes qui sont exécutés.

Cela a beaucoup plus de sens lorsque vous avez des workflows plus importants et exécutez des centaines ou des milliers d'échantillons différents. Vous pouvez simplement générer tellement de sortie sur le terminal, c'est impossible à regarder, alors que cette vue de mise à jour vous donne une progression en temps réel.

## 2.1.3. Exécuter à nouveau la commande avec l'option -ansi-log false

Si vous le souhaitez, vous pouvez l'exécuter à nouveau, et cette fois je vais utiliser un argument de base Nextflow supplémentaire avec un seul trait d'union disant, "-ansi-log false". Cela utilise la version précédente de la sortie du journal Nextflow. Et ici vous pouvez voir tous les processes individuels qui ont été lancés.

C'est à vous de décider si vous faites cela ou non. La sortie de Nextflow est exactement la même dans les deux cas.

## 2.2. Assurer que les noms de fichiers de sortie seront uniques

D'accord, jetons un coup d'œil aux fichiers de sortie, puis nous irons dans results. Mais il n'y a qu'un seul fichier de sortie. Que s'est-il passé ? Nous avons vu que le process avait été exécuté plusieurs fois. Nous pouvons aller dans le répertoire work et voir tous les différents hashs, toutes les tâches ont été exécutées correctement. Mais si vous vous souvenez dans notre process ici, nous sauvegardons tout dans un fichier output.txt et ensuite nous publions cela dans ce répertoire.

Donc le même fichier a été créé cinq fois, puis il a été écrasé cinq fois. Et nous avons juste quelle que soit la tâche qui s'est exécutée en dernier.

## 2.2.1. Construire un nom de fichier de sortie dynamique

La façon dont nous résolvons cela est d'utiliser un nom de fichier de sortie dynamique. Ici nous avons déjà une variable appelée greeting dans le process, nous pouvons donc l'utiliser dans le nom du fichier de sortie. Je copie cela et je fais $greeting-output.txt.

Je vais entourer cela de guillemets, juste pour que bash ne soit pas confus par des espaces qui pourraient se glisser ici. Et ensuite je vais prendre le même nom de fichier et mettre à jour la sortie ici.

C'est vraiment important que la sortie corresponde à cela, parce que sinon, ce fichier ne sera pas trouvé et Nextflow plantera.

Je vais faire une autre modification vraiment importante, qui est que je vais changer ces guillemets simples pour des guillemets doubles. Notez que la couleur du code a changé lorsque j'ai fait cela. Cette variable n'est développée que si nous utilisons des guillemets doubles. Si j'utilise des guillemets simples ici, elle est utilisée comme valeur littérale, et j'obtiendrais un seul fichier appelé $greeting-output, ce qui n'est pas ce que je veux.

## 2.2.2. Exécuter le workflow

Alors remettons les guillemets doubles et essayons.

Je vais juste nettoyer mon répertoire avant de commencer, pour qu'il soit facile de voir les nouveaux fichiers. Je vais supprimer tout ce qui s'appelle .nextflow, work, et results.

Et je vais exécuter cette commande Nextflow à nouveau et voyons quels fichiers sont créés. Donc il exécute les cinq processes là. Si vous regardiez très attentivement, vous auriez pu voir cette ligne se mettre à jour pendant l'exécution.

Et maintenant nous pouvons aller dans le répertoire results, et effectivement, nous avons cinq sorties différentes, et elles sont toutes préfixées avec la différente salutation.

Si j'ouvre chacune d'elles, nous verrons qu'elles contiennent chacune la salutation correspondante. Fantastique. C'est ce que nous voulons.

## 3. Utiliser un opérateur pour transformer le contenu d'un channel

D'accord, donc maintenant nous savons ce que sont les channels et nous savons ce que sont les channel factories. Qu'en est-il des opérateurs ? C'est un autre terme pour une partie du langage Nextflow, qui est une série de fonctions qui nous permettent d'opérer sur les channels pour leur faire certaines choses. Nextflow, vient avec une suite d'opérateurs, qui nous permettent de manipuler les channels de diverses manières différentes.

## 3.1. Fournir un tableau de valeurs comme entrée au channel

Travaillons à travers ceci avec un exemple. Disons que nous voulons prendre ces chaînes d'entrée, mais au lieu de les mettre directement dans une channel factory, nous voulons les définir comme un tableau.

## 3.1.1. Configurer la variable d'entrée

Je vais donc prendre celles-ci et faire cela comme une nouvelle ligne au-dessus et dire, greetings, array.

Voilà. Je vais prendre cette variable de tableau et la mettre dans le channel.of, et appuyer sur sauvegarder.

## 3.1.3. Exécuter le workflow

Maintenant, voyons ce qui se passe. Retournons à mon terminal. Je vais juste nettoyer tous ces fichiers temporaires à nouveau. Et exécutons le workflow.

Pas bon. D'accord. Il a planté. Ce n'est pas grave. Je m'attendais à ce qu'il plante cette fois. Déboguer ce qui ne va pas lorsqu'un workflow Nextflow échoue est une partie clé d'être un·e développeur·se Nextflow. Cela arrivera beaucoup et il est important de comprendre ce que dit le message d'erreur et comment y faire face.

Les messages d'erreur Nextflow sont en fait assez structurés. Il nous dit quel process s'est mal passé. Il nous donne un message d'erreur pour une raison. Il dit quelle était la commande qu'il a essayé d'exécuter dans cette tâche particulière, quel était le statut de sortie, quelle était la sortie et où était le répertoire work de cette tâche.

Notez que je peux option-cliquer dessus dans VS Code et cela l'ouvre dans une barre latérale pour que je puisse y aller directement et voir tous ces fichiers cachés, dont nous avons parlé dans le chapitre précédent, y compris le fichier .command.sh. Vous pouvez voir que c'est le même que les commandes qui ont été exécutées ici.

En regardant ce fichier, nous pouvons avoir une idée de ce qui aurait pu mal tourner ici au lieu d'exécuter une seule tâche pour chaque élément du tableau comme il l'a fait la dernière fois, il a juste fourni le tableau entier en une seule fois comme une chaîne. Nous devons donc déballer ce tableau en valeurs individuelles avant de le passer dans le channel. Retournons en arrière et voyons si nous pouvons faire cela en utilisant un opérateur.

## 3.2. Utiliser un opérateur pour transformer le contenu du channel

Dans ce cas, nous n'allons pas changer le tableau avant de le passer dans le channel. Nous allons ajuster le channel pour qu'il se comporte de la manière que nous attendons. Nous allons faire cela en utilisant l'opérateur flatten peut faire dot commencer à taper et nous pouvons voir que l'extension VS Code commence à suggérer tous les différents opérateurs que nous avons disponibles.

## 3.2.1. Ajouter l'opérateur flatten()

Et je vais sélectionner flatten. Notez que l'espace blanc n'a pas d'importance dans ce contexte pour Nextflow. Vous pouvez donc mettre ces opérateurs sur une nouvelle ligne si vous le souhaitez. Je peux donc faire descendre cela ici et l'indenter pour qu'il se situe sous ".of" et vous verrez que les gens enchaînent souvent beaucoup d'opérateurs comme ceci sur un channel et l'indentent de cette façon pour qu'il soit plus facile à lire.

Vous pouvez également voir, comme avant, je peux survoler cela et lire ce que fait l'opérateur flatten, et aussi suivre un lien vers la documentation si je le souhaite.

Donc cet opérateur prend ce channel, qui a un seul tableau dedans, et sépare les valeurs du tableau.

## 3.2.2. Ajouter view() pour inspecter le contenu du channel

Nous pouvons jeter un coup d'œil dans les channels en utilisant l'opérateur spécial view, et je vais en ajouter quelques-uns ici. C'est un peu comme utiliser des instructions print dans d'autres langages. Je vais donc faire dot view et ensuite je vais utiliser ces accolades.

Ceci s'appelle une closure. Cela donne essentiellement du code supplémentaire à l'opérateur view, qu'il exécutera sur chaque élément dans le channel. Dans ce cas, je vais dire greeting before flatten. Greeting.

Je définis une variable ici, qui n'est que dans la portée de cette closure. Donc cette variable n'est utilisée qu'ici et je pourrais l'appeler comme je le veux. Peu importe vraiment. J'utilise juste greeting pour le rendre facile à lire.

Dans certains pipelines Nextflow, vous pourriez voir des gens utiliser une variable implicite spéciale appelée "$it". Comme ça. C'est une variable spéciale dans le code Nextflow, qui est un raccourci pour que vous n'ayez pas à faire la petite définition d'une variable. Cependant, au fil du temps, nous pensons que ce n'est pas très clair pour les personnes qui sont nouvelles à Nextflow, et nous décourageons l'utilisation de "$it" maintenant.

Je vais donc m'en tenir au comportement précédent de greeting et l'utiliser comme ceci parce que c'est plus explicite et c'est plus clair sur ce qui se passe.

Je vais ensuite copier cette ligne et faire exactement la même chose à nouveau après les arguments flatten. L'opérateur view est un peu spécial parce qu'il fait quelque chose sur les éléments, mais il continue également simplement à les passer au prochain opérateur pour que nous puissions l'enchaîner au milieu d'une chaîne d'opérations comme ceci, et il imprimera l'état là et continuera. J'espère donc que cela nous montrera à quoi ressemble le channel avant et après l'opérateur flatten.

## 3.2.3. Exécuter le workflow

Essayons-le. Nettoyer. Nettoyer tout dans l'espace de travail. Exécuter à nouveau le pipeline.

D'accord, nous pouvons voir qu'il a exécuté nos cinq processes. Encore une fois, il n'a pas planté avec une erreur, c'est donc définitivement bon. Et maintenant nous avons le before flatten et effectivement nous avons notre tableau et nous avons after flatten, imprimé cinq fois une fois pour chaque élément du tableau. C'est exactement ce que nous espérions. C'est donc vraiment une bonne nouvelle. Et cela correspond exactement à ce que nous attendrions du code.

Nous n'avons plus besoin de ces instructions de débogage, je peux donc soit les commenter soit les supprimer. Je vais les supprimer juste pour garder mon code propre et net. D'accord, super. Cet exemple fonctionne maintenant bien et nous pouvons commencer à voir comment les channels peuvent faire une logique un peu plus compliquée.

## 4. Utiliser un opérateur pour analyser les valeurs d'entrée à partir d'un fichier CSV

Maintenant nous allons essayer de faire cela en utilisant un fichier avec une série d'entrées à la place. C'est une façon très courante d'écrire des pipelines Nextflow en utilisant une feuille d'échantillons ou un CSV de métadonnées.

## 4.1. Modifier le script pour attendre un fichier CSV comme source de salutations

Si je vais dans la barre latérale, vous pouvez voir greetings.csv dans le dépôt d'exemples, et c'est un fichier CSV très, très simple qui contient juste trois lignes avec trois salutations différentes. Voyons si nous pouvons utiliser ce fichier CSV dans notre workflow.

Je vais maintenant revenir à l'utilisation de params comme nous l'avons fait dans le chapitre un, pour que nous puissions avoir une entrée en ligne de commande.

Je vais supprimer ce tableau greetings.

## 4.1.1. Changer le paramètre d'entrée pour pointer vers le fichier CSV

Je vais définir params greeting au nom de fichier, qui est greetings.csv, et je vais utiliser cette variable spéciale pour générer le channel. Je vais mettre ça là-dedans, et les erreurs disparaissent. Rappelez-vous que cela définit cette variable par défaut maintenant. Donc si j'exécute le pipeline sans aucun argument, il utilisera greetings.csv, mais je pourrais faire --greeting pour écraser cette variable si je le voulais.

## 4.1.2. Changer pour une channel factory conçue pour gérer un fichier

D'accord, nous passons maintenant un fichier plutôt qu'une chaîne ou un tableau de chaînes, nous avons donc probablement besoin d'une channel factory différente.

Nous allons nous débarrasser de "of" que nous avons utilisé jusqu'à présent, et utiliser à la place .fromPath. Cela fait exactement ce que cela semble. Il crée un channel avec des chemins au lieu de valeurs, en utilisant un nom de fichier de chaîne ou un glob. Je vais également supprimer l'opérateur flatten car nous n'en avons plus besoin maintenant que nous passons un fichier.

## 4.1.3. Exécuter le workflow

Je vais appuyer sur sauvegarder, ouvrir le terminal, exécuter le workflow, et ensuite voir ce qui se passe.

D'accord. Il a planté à nouveau. Ne vous inquiétez pas. Je m'attendais à celui-ci également. Jetons un coup d'œil au message d'erreur et voyons si nous pouvons comprendre ce qui ne va pas. Ici nous pouvons voir la commande exécutée, et un peu comme avant où nous avions tout le tableau imprimé. Maintenant nous avons le chemin de fichier qui est écho dans la commande, plutôt que de parcourir le contenu du fichier.

## 4.2. Utiliser l'opérateur splitCsv() pour analyser le fichier

Donc pour utiliser le contenu du fichier à la place, nous avons besoin d'un autre opérateur. L'opérateur que nous allons utiliser pour celui-ci s'appelle splitCsv. Cela a du sens, parce que c'est un fichier CSV que nous chargeons.

## 4.2.1. Appliquer splitCsv() au channel

Ok, donc splitCsv. Fermer la parenthèse. Nous n'avons besoin d'aucun argument ici. Et encore, je vais utiliser quelques opérateurs view pour donner un aperçu de ce qui se passe ici.

.view csv after splitCsv. Before split Csv.

## 4.2.2. Exécuter à nouveau le workflow

D'accord, essayons d'exécuter cela et voyons ce qui se passe.

D'accord, nous avons un peu plus de sortie cette fois, mais cela a quand même échoué. Nous pouvons regarder les instructions view, et ici vous pouvez voir before split CSV, et nous avons un chemin de fichier comme nous l'avons vu dans le message d'erreur précédent. After split CSV, nous avons maintenant trois valeurs correspondant aux trois lignes dans le fichier CSV.

Cependant, vous pouvez voir que chacune de ces valeurs est entourée de crochets. Donc chacune d'elles était un tableau en soi, et cela nous a donné la même zone que nous avions avant où il essaie d'écho un tableau plutôt qu'une simple chaîne.

Si nous pensons à un fichier CSV, cela a du sens. Typiquement, un fichier CSV aura des lignes et des colonnes, donc split CSV fait un tableau en deux dimensions. La première dimension du tableau est chaque ligne, et ensuite il y a une deuxième dimension, qui est chaque colonne pour chaque ligne.

Donc ici nous n'avons qu'une seule valeur sur chaque ligne, donc nous avons une seule colonne, donc nous avons un tableau à un élément pour chaque ligne du fichier.

Ce n'est pas grave. Nous avons juste besoin d'un autre opérateur pour effondrer ce tableau pour chaque ligne du fichier CSV analysé. Nettoyons ceci. Débarrassons-nous d'un terminal et voyons ce que nous pouvons faire.

## 4.3. Utiliser l'opérateur map() pour extraire les salutations

Maintenant nous pourrions utiliser l'opérateur flatten à nouveau, que nous avons utilisé avant. Nous avons vu comment cela peut effondrer un tableau en une série de valeurs, ce qui fonctionnerait très bien ici. Mais je vais utiliser l'opportunité pour démontrer un autre opérateur, qui est très courant dans les workflows appelé l'opérateur map.

## 4.3.1. Appliquer map() au channel

Je vais faire dot map et je vais faire item item[0].

Si vous écrivez beaucoup d'autres langages de code, vous pourriez déjà être familier avec l'opérateur map. Il prend un itérable, tel qu'un tableau ou un channel, et il fait une opération sur chaque valeur de celui-ci.

Ici nous disons que nous devrions définir une variable appelée item dans la portée de cette closure, et ensuite nous voulons retourner, juste la première valeur dans ce tableau. Donc item index zéro.

Cela aplatit effectivement le tableau. Vous pouvez voir comment nous pourrions étendre ceci pour être plus complexe, cependant : si notre fichier CSV avait six colonnes, mais nous ne sommes intéressés que par la quatrième colonne, nous pourrions accéder à un index spécifique ici. Ou faire tout autre type d'opération sur la valeur avant de la passer au traitement en aval.

Donc l'opérateur map est extrêmement flexible et très puissant pour modifier les channels en vol. Mettons une autre instruction view juste pour que nous puissions voir ce qu'il fait dans notre exécution. Peut juger cette ligne et la déplacer vers le bas. Et after map.

## 4.3.2. Exécuter le workflow une dernière fois

Ouvrons le terminal et essayons d'exécuter le workflow.

D'accord, pas d'erreurs cette fois. C'est bon signe. Nous pouvons maintenant parcourir toutes ces différentes sorties des instructions view. Before split CSV, nous avions un seul chemin. After split CSV, nous avions les tableaux à valeur unique, et ensuite after map, nous avons juste les valeurs sans aucune syntaxe de tableau. Allons dans le répertoire results, et voici nos fichiers de sortie se comportant exactement comme nous le voulions.

Il y a un petit bonus ici. Vous pouvez réellement voir que les opérateurs view sont légèrement mélangés dans l'ordre dans lequel ils ont fait la sortie. C'est parce que Nextflow fait de la parallélisation de ces différentes tâches. Donc après avoir divisé le CSV, il y a trois éléments dans ce channel, et il gère le traitement de ces trois éléments en parallèle automatiquement. Cela signifie que l'ordre des sorties est stochastique et peut varier. Dans ce cas, il s'est simplement trouvé que certains des opérateurs view ont retourné après que l'étape suivante ait été terminée, et donc c'est venu dans cet ordre.

Si j'exécute le même workflow à nouveau. Alors effectivement, c'est venu dans un ordre différent et cette fois nous avons les split CSVs et les maps dans l'ordre que nous attendrions.

Alors gardez juste à l'esprit, vous ne pouvez pas compter sur l'ordre des sorties d'une tâche de process parce que Nextflow gère cette parallélisation pour vous automatiquement. Nextflow fait cela pour vous avec sa logique de flux de données, et c'est la vraie puissance de Nextflow.

D'accord, c'est probablement l'un des chapitres les plus importants de toute la formation. Une fois que vous comprenez les channels, les channel factories et les opérateurs, vous commencez à exploiter la force de Nextflow et ce qui le rend unique en tant que langage de programmation. Cette fonctionnalité permet à Nextflow de paralléliser tous vos workflows pour vous et de générer une logique de workflow extrêmement complexe avec une syntaxe très propre et un modèle de flux de données push. Cela peut être un concept un peu étrange au début, mais une fois que vous vous habituez à écrire du code comme ceci, cela deviendra rapidement naturel et avant que vous le sachiez, vous écrirez des workflows fantastiques.

Faites une pause, une tasse de thé, promenez-vous et passons au chapitre trois, où nous commençons à étendre ces concepts dans des workflows plus complexes. À la prochaine vidéo.

[Transcription de la vidéo suivante :octicons-arrow-right-24:](03_hello_workflow.md)
