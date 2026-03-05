# Partie 5 : Hello Containers - Transcription vidéo

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Traduction assistée par IA - [en savoir plus et suggérer des améliorations](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

<div class="video-wrapper">
  <iframe width="560" height="315" src="https://www.youtube.com/embed/Xqr--bKEN9U?si=QinuAnFwFj-Z8CrO&amp;list=PLPZ8WHdZGxmWKozQuzr27jyMGqp9kElVK&amp;cc_load_policy=1" title="YouTube video player" frameborder="0" allow="accelerometer; autoplay; clipboard-write; encrypted-media; gyroscope; picture-in-picture; web-share" referrerpolicy="strict-origin-when-cross-origin" allowfullscreen></iframe>
</div>

!!!note "Notes importantes"

    Cette page présente uniquement la transcription. Pour les instructions complètes étape par étape, retournez au [matériel de cours](../05_hello_containers.md).

    Les numéros de section affichés dans la transcription sont fournis à titre indicatif uniquement et peuvent ne pas inclure tous les numéros de section du matériel.

## Bienvenue et contexte

Bonjour et bienvenue dans Hello Nextflow. Voici la partie cinq intitulée Hello Containers. Dans cette partie du cours, nous allons parler de la façon d'encapsuler les exigences logicielles d'un pipeline afin que les personnes qui exécutent le pipeline n'aient pas à se soucier de l'installation des logiciels.

Si vous travaillez en bioinformatique depuis aussi longtemps que moi, vous vous souvenez peut-être de ce que j'appelle souvent la mauvaise époque, où lorsque vous vouliez exécuter le pipeline de quelqu'un d'autre ou reproduire son travail, vous passiez des heures ou des jours à essayer d'installer tous les différents outils logiciels qu'ils utilisaient, dans les mêmes versions, en essayant de les compiler sur votre machine, et c'était un cauchemar. C'était vraiment difficile.

Si vous travailliez sur un HPC, vous avez peut-être utilisé des modules d'environnement où les administrateurs système essayaient d'installer des logiciels pour vous, ce qui était correct, mais toujours imparfait.

Mais maintenant nous avons de meilleures façons de faire cela. Nextflow a un support intégré pour différentes technologies de conteneurs logiciels. Docker est la plus courante. C'est celle que nous allons utiliser aujourd'hui. Elle fonctionne bien dans Codespaces. Elle fonctionne bien sur votre ordinateur local et elle fonctionne bien dans le cloud.

Mais aussi Singularity ou Apptainer, qui sont très courants sur les systèmes HPC et fonctionnent effectivement exactement de la même manière. Ou Podman, Shifter, il y a plein d'autres qui sont tous très similaires.

Le seul autre, qui est un peu similaire mais pas tout à fait, que Nextflow supporte est Conda. Et Nextflow peut gérer les environnements Conda pour vous sur une base par processus, ce qui est bien mieux que de gérer vos propres environnements Conda. Et encore une fois, peut être livré avec un pipeline.

Nous allons commencer ce chapitre en parlant un peu des technologies de conteneurs et de Docker et de leur fonctionnement. Et nous allons faire la première moitié manuellement dans Docker afin que vous compreniez ce qui se passe sous le capot et comment cela fonctionne. Parce que c'est vraiment important de comprendre ce que Nextflow fait et comment comprendre ce que votre workflow fait lorsqu'il est exécuté.

Alors. Passons à notre Codespaces. Maintenant j'ai tout nettoyé à nouveau, mais si nous allons dans Hello Containers, vous devriez voir que tous nos scripts et tout sont là, identiques à la fin du chapitre sur les modules. Nous avons donc nos différents modules ici, que j'ai créés dans le répertoire modules.

Ils sont toujours là. Ils doivent être là pour que cela puisse s'exécuter. Et le workflow et la sortie sont tous identiques sauf que nous avons changé le chemin de publication de sortie vers Hello Containers, afin que vos fichiers se retrouvent dans ce répertoire.

Nous pouvons exécuter cela maintenant pour vérifier que ça fonctionne si vous voulez, ou nous pouvons continuer avec le terminal.

## 1. Utiliser un conteneur 'manuellement'

Nous allons utiliser Docker pour gérer nos conteneurs, et je peux vérifier qu'il est installé sur mon Codespaces en faisant "docker -v", ce qui me montre la version qui est installée et tout, et que ça fonctionne correctement.

Maintenant les conteneurs et Docker ont deux concepts qui sont vraiment importants. L'un s'appelle image, et l'autre s'appelle conteneur. L'image est l'instantané, si vous voulez, de tout le système de fichiers que vous utiliserez, et le conteneur est l'environnement d'exécution. Vous créez donc un conteneur en utilisant une image.

Une fois que vous êtes dans ce conteneur, il fonctionne généralement comme un système d'exploitation complet. Il est coupé du monde extérieur. Il est séparé de tout le reste, et c'est une bonne chose. C'est ainsi que nous obtenons une si bonne reproductibilité avec Nextflow.

Parce que pour les tâches exécutées à l'intérieur d'un conteneur, elles ne sont pas contaminées par des fichiers de configuration sur votre système local. Aucune autre influence externe, elles s'exécutent dans leur propre petit bac à sable. Les fichiers sont ensuite produits de manière très, très reproductible parce que vous utilisez les mêmes bibliothèques sous-jacentes, toutes les mêmes dépendances, exactement le même logiciel pour chaque personne s'exécutant sur chaque environnement informatique différent. Ce qui franchement je pense est fantastique et incroyable que cela fonctionne. Et même, même à ce jour, cela me stupéfie encore que ce soit possible.

## 1.1. Récupérer l'image du conteneur

Nous allons donc essayer d'utiliser des images Docker et Docker, lorsque vous l'exécutez sur votre système, a un registre docker sur votre ordinateur, ou dans ce cas, dans le code space, qui garde une trace de toutes les différentes images qui ont été téléchargées et utilisées dans le passé, et des différentes couches dont elles sont constituées.

Nous pouvons voir quelles images nous avons localement avec Docker en faisant "docker image ls". Et dans ce cas, vous pouvez voir qu'il y a un tas d'images Docker ici, qui ont toutes à voir avec la configuration de ce Codespaces. Tout à voir avec les conteneurs de développement et autres. Donc vous n'avez pas besoin de trop vous en soucier, mais au fur et à mesure que nous ajoutons plus d'images et les téléchargeons, au fur et à mesure que ce cours avance, vous pouvez vérifier cette liste et vous verrez que le registre local garde une trace de toutes ces choses que nous avons récupérées.

Mais nous allons en récupérer une nouvelle en faisant "docker pull". Et cela dit à Docker de récupérer une nouvelle image depuis le web.

Nous mettons ensuite l'URI pour ce conteneur. Maintenant cela pourrait être une image docker que vous avez construite localement et ensuite poussée sur internet. Cela pourrait être une image que quelqu'un d'autre a faite. Il y a beaucoup, beaucoup, beaucoup de façons différentes de faire des images Docker, mais sans doute l'une des façons les plus simples est d'externaliser cela, et de laisser quelqu'un d'autre le faire pour vous.

Et ce que nous allons utiliser dans ce tutoriel est un service de Seqera appelé Seqera Containers.

Maintenant, Seqera Containers est totalement gratuit, et il utilise un logiciel open source que nous développons appelé Wave, qui a été construit pour gérer les conteneurs de manière complémentaire à Nextflow. Et il gère de nombreux cas d'usage courants auxquels nous sommes confrontés avec Nextflow.

Il est très courant que le logiciel dont nous avons besoin soit empaqueté dans Conda, dans les canaux Bioconda ou conda-forge ou d'autres canaux plus spécifiques au domaine. Et Wave et Seqera Containers sont vraiment bons pour construire des images à partir de cela.

Je peux donc aller sur cette interface web et nous allons jouer avec le package appelé "cowpy". Donc je tape le nom du package que je veux. Il recherche, il l'a trouvé sur l'index de packages Python, donc je peux utiliser cela. Ou si j'attends un peu plus longtemps, il recherche bioconda et conda-forge. Et vous pouvez voir, je peux spécifier n'importe quel canal conda ici. Donc si vous voulez trouver un canal Nvidia ou autre chose, cela devrait fonctionner aussi.

Et ensuite je peux spécifier si je veux qu'il construise une image docker pour moi ou une image singularity et aussi quelle architecture CPU je veux utiliser. Donc amd64 ou arm64.

Et une fois que les résultats bioconda sont listés, je peux maintenant voir toutes les différentes versions qui sont disponibles également. Je vais mettre ça dedans. Et maintenant je pourrais continuer à chercher et obtenir plus de packages de Conda si je veux et composer ce conteneur comme je veux, mais je veux juste celui-là. Donc je vais cliquer sur Get Container.

Maintenant, quelqu'un d'autre a déjà demandé le même conteneur auparavant et il est retourné depuis un registre, donc nous l'obtenons juste immédiatement. Mais si personne d'autre n'avait jamais demandé ce package logiciel ou cette combinaison de packages logiciels, Wave et Seqera Containers le construiraient à la volée pour nous.

Nous pouvons copier cette URL et nous pouvons également voir les détails de construction. Et cela nous montre ce que le service a fait en arrière-plan. Il a créé un fichier d'environnement conda. Un fichier docker, et ensuite c'est ça, exécutant le processus de construction docker. Il a également exécuté une analyse, une analyse de sécurité, donc vous pouvez voir toutes les CVE. Et il vous dit quand cela a été créé.

Wave et Seqera Containers peuvent faire beaucoup plus que cela, mais c'est un cas d'usage simple, qui est le plus courant. Et je dois dire que ces images sont hébergées pendant au moins cinq ans. Vous pouvez donc intégrer ces URL dans vos pipelines et savoir qu'elles ne vont pas disparaître de sitôt.

J'ai donc mon URL pour mon image docker pour cowpy.

Je peux maintenant faire "docker pull" cette URL, et cela va récupérer toutes les différentes couches et télécharger cette image pour qu'elle soit disponible pour moi localement.

## 1.2. Utiliser le conteneur pour exécuter cowpy comme une commande unique

D'accord, maintenant essayons de l'utiliser réellement. Donc maintenant je vais maintenant utiliser une commande "docker run" au lieu de "docker pull", et je vais utiliser le flag "--rm", qui dit juste à Docker d'arrêter ce conteneur une fois qu'il a fini de faire ce que je lui ai demandé. Et ensuite je mets l'identifiant pour le conteneur, qui est juste un URI.

Et ensuite à la fin, je spécifie la commande que je veux que Docker exécute à l'intérieur du conteneur généré à partir de cette image. Je vais juste dire cowpy, qui est le nom de l'outil qui est installé depuis Conda Forge, qui est disponible à l'intérieur de l'image.

Je vais appuyer sur entrée et voilà. Nous avons exécuté cowpy sur un système. Nous avons une petite vache qui nous donne des informations.

Maintenant notez que cowpy n'est pas installé sur mon système local. Donc si je l'exécute juste sans tout le truc Docker, il dit, commande non trouvée. Donc cela a récupéré une image. Il a créé un conteneur en utilisant Docker, et ensuite il est allé dans ce conteneur et a exécuté cette commande pour nous et nous a donné la sortie dans notre terminal. Très, très cool.

## 1.3. Utiliser le conteneur pour exécuter cowpy de manière interactive

D'accord, nous allons aller un peu plus loin maintenant et exécuter ce conteneur de manière interactive et jeter un œil à l'intérieur, pour que nous puissions voir ce qui se passe à l'intérieur du conteneur.

Donc si je remonte et je prends ma commande run et je vais me débarrasser de cowpy à la fin là, parce que je ne veux pas vraiment exécuter cowpy. Je veux exécuter un terminal Bash.

Et ensuite je vais revenir ici et je vais faire "-it", qui signifie Interactive et Terminal ou TTY, et je vais appuyer sur entrée.

Et maintenant vous pouvez voir l'invite, la partie avant que je tape, a changé. C'était l'invite Codespaces où il disait le répertoire, et maintenant il dit base et roots et tmp. Donc je suis maintenant à l'intérieur du conteneur, et si je fais "ls", vous verrez que les fichiers que je vois dans ce répertoire sont différents des fichiers que j'ai dans mon espace de travail.

Et en fait, je ne peux voir aucun des fichiers de mon espace de travail codespaces local ou de mon disque local à l'intérieur du conteneur Docker. Le runtime du conteneur docker est complètement isolé et il ne peut pas écrire ou lire des fichiers depuis un système de fichiers hôte à l'extérieur.

Je peux, cependant, voir le logiciel qui est installé à l'intérieur du conteneur et l'exécuter. Donc je peux exécuter cowpy et nous pouvons voir un peu plus sur comment utiliser cowpy. Ici je peux faire "cowpy 'Hello World'" et cela lui dit de mettre ma citation à l'intérieur d'une petite bulle de dialogue. Et vous pouvez également exécuter différents types de vaches, donc ça n'a pas besoin d'être une vache. Vous pouvez faire un "-c". Et je suis en Suède, donc je vais choisir un élan. Très bien. Je lui ai donné des bois.

Et il y a tout un tas de différents que vous pouvez essayer, que vous pouvez voir décrits dans les documents de formation.

## 1.3.4. Monter des données dans le conteneur

D'accord. Ce serait bien si nous pouvions exécuter cowpy sur les fichiers de notre système de fichiers.

Bien sûr, ce n'est pas très utile d'avoir juste le conteneur et aucun accès à quoi que ce soit du tout. Cela pourrait être sûr et reproductible, mais ce n'est pas très utile.

Alors comment faisons-nous cela ? Je vais sortir de ce conteneur Docker en tapant exit, et vous pouvez voir que l'invite nous dit que nous sommes maintenant de retour dans notre Codespaces normal.

Et je vais exécuter la même commande à nouveau. Mais cette fois je vais ajouter quelques flags supplémentaires ici. Et l'important est "-v", qui signifie monter un volume, qui est comme essentiellement une partie, partie d'un espace disque.

Le "-v" prend deux parties : il y a comme une chaîne et ensuite deux points et une chaîne. Et la première partie est le système de fichiers local, qui devrait être monté dans le conteneur. Et ensuite la deuxième partie est où cela devrait se retrouver à l'intérieur du conteneur.

Maintenant je veux juste charger tout mon système de fichiers local ici. Donc "." est le répertoire de travail actuel. Donc je vais juste faire "." et ensuite ":", et ensuite nous allons mettre cela dans un nouveau répertoire à l'intérieur du conteneur appelé "my_project". Cela pourrait vraiment s'appeler n'importe quoi.

Et ensuite je vais exécuter à nouveau.

Dans le répertoire de travail où je suis déposé, qui est /tmp, les fichiers ne sont pas là. Mais si je fais "ls my_project", nous y voilà : tous les mêmes fichiers que nous avions localement sur notre Codespaces sont maintenant disponibles à l'intérieur du conteneur à ce chemin.

C'est un accès en lecture et écriture donc je peux créer de nouveaux fichiers dans ce répertoire et ils apparaîtront sur mon système de fichiers hôte. Ce répertoire particulier, se comporte alors exactement comme si j'étais à l'extérieur du conteneur donc je peux maintenant lire et écrire et faire des choses.

## 1.3.5. Utiliser les données montées

D'accord, prouvons juste que nous pouvons faire cela. Je fais "cat /my_project/data/greetings.csv". Si vous vous souvenez, le contenu de ce fichier ressemble à ceci. Je peux maintenant rediriger cela vers cowpy et la vache imprimera les différentes sorties de ce fichier dans sa petite bulle de dialogue, ce qui est plutôt amusant.

Vous pouvez donc voir, nous pouvons maintenant utiliser le logiciel dans le conteneur pour interagir avec les fichiers sur notre système hôte.

D'accord, revenons en arrière et nous allons continuer avec le reste du matériel de formation.

## 2. Utiliser des conteneurs dans Nextflow

C'est donc vraiment cool d'utiliser des conteneurs. J'espère que cela a du sens. Et vous pouvez voir la valeur de ces conteneurs et pourquoi c'est utile pour exécuter des logiciels d'analyse.

Mais comment faisons-nous tout ce même processus à l'intérieur de Nextflow ? Nous ne voulons pas exécuter plein de commandes Docker nous-mêmes. Nous voulons juste laisser Nextflow gérer tout cela pour nous.

Alors travaillons sur cela. Nous allons ajouter un nouveau processus à notre pipeline, pour exécuter cowpy. D'accord, créons donc un nouveau module pour notre nouveau processus. Donc allons dans modules, appelons-le cowPy.nf, et ensuite je vais copier le code du matériel de formation ici.

Mais vous pouvez voir que le processus est très simple. Il ressemble beaucoup à ceux que nous avons faits jusqu'à présent, nous avons un bloc input avec un path, qui est notre fichier d'entrée, et aussi une value ici donc ce sera un caractère, donc nous pourrions utiliser un élan à nouveau si nous voulons.

Et ensuite une output, qui est un seul fichier ici, un path et ensuite un script. Et nous faisons la même chose que nous avons faite de manière interactive à l'intérieur du conteneur : nous faisons "cat" pour lire le fichier d'entrée. Nous redirigeons ce contenu vers cowpy. Nous choisissons un caractère spécifique basé sur cette entrée, nous écrivons dans un fichier de sortie appelé cowpy input file, qui est ensuite émis vers la sortie.

Parfait. Incluons cela. Donc include \{ cowpy \} from "./modules/cowpy.nf", l'ai-je appelé cowpy ? Oui.

Et ensuite appelons notre nouveau processus ici en bas dans le bloc main du workflow. Donc exécutons cowpy. Et nous allons prendre notre nouveau processus cowpy et nous allons dire collectGreetings.out.

Et ensuite si vous vous souvenez, il y avait deux sorties pour ce module. Une appelée outfile et une appelée report. L'extension VS Code nous suggère automatiquement celles-ci et nous voulons le .outfile.

Vous pouvez toujours aller dans ce processus ici. Soit vous survolez dessus et il devrait vous montrer rapidement quelles étaient les sorties. Et nous pouvons également faire commande clic dessus et il ouvrira le fichier module si vous voulez voir plus en détail.

Donc voilà. C'est le outfile là, et c'est le path. Donc ce sera maintenant le fichier d'entrée pour notre processus cowpy. Fantastique.

Maintenant si vous vous souvenez, un processus cowpy a deux entrées. Nous avions également le canal value pour le caractère. Donc nous pouvons ajouter "params.character" ici. J'aurais pu coder cela en dur si je voulais, mais faisons-en une option CLI pour que nous puissions faire dash, dash character.

Bien. Je dois maintenant définir le paramètre d'entrée que nous venons d'appeler et lui donner une valeur par défaut. Donc character, String. Et j'aime l'élan, donc je vais le définir sur moose par défaut.

Bien, essayons de l'exécuter. Donc si je fais Nextflow run hello containers, nous verrons ce qui se passe.

J'aurais pu utiliser dash resume si j'avais les anciens répertoires de travail qui traînaient. Et encore une fois, ces premiers processus auraient été mis en cache et ça aurait été un peu plus rapide, mais ça devrait être essentiellement la même chose.

Maintenant nous pouvons voir tout de suite qu'il a lancé une erreur quand il est arrivé à notre nouveau processus, il nous dit ici qu'il y a eu une erreur en exécutant le processus cowpy et il s'est terminé avec un statut de sortie 127. C'est la commande qu'il a essayé d'exécuter. Elle semble correcte, elle ressemble à ce que nous attendions. Elle prend ce nom de fichier de sortie, qui semble correct, elle l'exécute avec un caractère moose et essaie de sauvegarder.

Mais vous pouvez voir l'erreur de commande ici dit que la commande cowpy n'est pas trouvée. Et cela a du sens parce que nous n'avons pas encore dit à Nextflow d'utiliser un conteneur. Nous lui avons juste donné la commande cowpy. Et comme je l'ai dit avant, cowpy n'est pas installé sur notre système local. Donc quand il a essayé de l'exécuter, ça a échoué.

## 2.3.1. Spécifier un conteneur pour cowpy

Nous devons dire à Nextflow qu'il y a un conteneur disponible et qu'il peut l'utiliser. Alors comment faisons-nous cela ?

Si nous allons dans notre module ici, nous allons ajouter une nouvelle déclaration en haut appelée "container". Et nous allons ensuite définir cela sur une chaîne.

Maintenant, si vous vous souvenez, dans Seqera Containers, je peux copier cette URL et je la dépose juste entre guillemets ici.

Maintenant revenons et essayons de l'exécuter à nouveau.

Voyons si ça marche cette fois.

Malheureusement, ça échoue exactement de la même manière, même si maintenant nous avons défini un conteneur pour que le processus s'exécute. Donc pour utiliser notre image docker, nous devons dire à Nextflow d'activer l'utilisation de Docker lorsque nous exécutons le workflow.

Et nous allons faire cela en créant un nouveau fichier de configuration. Donc je vais dire touch nextflow.config.

C'est un nom de fichier spécial où s'il est dans le répertoire de travail pendant que je lance le pipeline, il sera chargé automatiquement. Donc si je vais dans ce fichier Nextflow dot config, vous pouvez voir qu'il existe déjà en fait, ce que j'avais oublié. Et nous avons docker.enabled ici déjà, mais il est défini sur false, ce qui est la valeur par défaut.

Donc si je change cela pour égale True à la place, docker.enabled. Et il y a des documents de référence pour toutes ces portées de configuration dans les documents Nextflow. Et aussi vous pouvez voir que lorsque je survole avec une extension VS Code, elle récupère les documents spécifiques à cela et me dit ce que cela signifie et comment le définir.

Donc maintenant nous l'avons défini sur true, et si j'exécute Nextflow à nouveau, Nextflow saura maintenant récupérer cette image docker pour nous si nous ne l'avons pas déjà localement, et ensuite exécuter ce processus avec cet environnement de conteneur.

Et donc nous pouvons voir qu'il s'est exécuté avec succès et nous avons une petite coche à côté de cowpy. Fantastique. Si je remonte et regarde dans le répertoire results, le fichier n'est pas encore là. Et c'est parce que nous devons encore publier ce fichier de sortie tout comme tous les autres.

Donc nous allons dans le bloc published dans le workflow, disons mycowpy égale cowpy.out.

Et ensuite ici en bas dans le bloc output, mycowpy, accolades path. Oups. Hello containers. Mode, copy.

Si j'exécute à nouveau maintenant, ça devrait s'exécuter exactement de la même manière. J'aurais pu utiliser dash resume et j'oublie à chaque fois. Et ensuite je remonte et maintenant nous avons un nouveau fichier créé appelé cowpy-COLLECTED, et voilà mon élan disant BONJOUR, HELLO, HOLà Fantastique.

Maintenant bien sûr je pourrais aussi passer maintenant "--character". Quelles sont les différentes options ? Je pense qu'il y a une dinde ? Donc je peux utiliser character Turkey. Ça va s'exécuter exactement de la même manière. J'ai raté une autre occasion d'utiliser dash resume, et maintenant si nous chargeons notre fichier et maintenant nous avons une dinde. Fantastique.

## 2.3.4. Inspecter comment Nextflow a lancé la tâche conteneurisée

D'accord. Dernière petite chose. Exécutons juste rapidement cette commande à nouveau, resume cette fois, et jetons un rapide coup d'œil dans le répertoire work pour voir ce que Nextflow fait sous le capot pour faire fonctionner tout cela pour nous.

Cette fois c'est super rapide, allons dans ce répertoire work, cd work/. Maintenant si vous vous souvenez nous avons un tas de fichiers point ici et celui qui nous intéresse dans ce cas est celui que j'ai dit que nous n'avons presque jamais besoin de regarder, appelé .command.run.

Si je fais code dot command run, ça va l'ouvrir dans l'éditeur. Et je peux chercher dans ce fichier et si je descends je devrais voir Docker run. Et vous pouvez voir que Nextflow fait la commande docker run pour nous, quand Docker est activé dans une config. Il a tout un tas de flags différents et de choses ici, mais vous pouvez voir le flag "-v" que nous avons utilisé nous-mêmes quand nous exécutions. Et vous pouvez voir qu'il monte le répertoire de l'espace de travail local dans le conteneur, pour que le conteneur puisse accéder à nos fichiers d'entrée et sauvegarder les sorties. Et ensuite à la fin, il exécute également .command.sh, qui est le script généré, qui a la commande cowpy dedans.

Et donc vous pouvez voir que Nextflow prend la logique du workflow, qui est le truc qui nous intéresse vraiment, qui est spécifique à notre analyse, et il fait tout le truc intelligent en coulisses pour faire fonctionner Docker sur notre système.

Et il fait cela d'une manière vraiment portable pour qu'un utilisateur final du pipeline puisse changer la technologie qu'il utilise : Docker, Singularity, Apptainer, Conda. Cela n'a pas vraiment d'importance pour la logique du pipeline, mais Nextflow gérera tous les besoins d'infrastructure sous-jacents, pour qu'il s'exécute n'importe où.

Et c'est vraiment le super pouvoir de Nextflow. C'est la reproductibilité et la portabilité. Et avec Nextflow vous pouvez réellement partager votre workflow et d'autres personnes peuvent l'exécuter sur leurs systèmes et ça fonctionnera juste.

C'est une chose vraiment, vraiment difficile à faire, et maintenant vous savez comment le faire aussi avec vos workflows.

D'accord, c'est tout pour ce chapitre. Si vous descendez à la fin du cours, vous trouverez un quiz à nouveau sur les conteneurs. J'espère que tout cela avait du sens. C'est une façon vraiment cool de travailler avec l'analyse. Et si vous êtes nouveau dans les conteneurs, j'espère que je vous ai convaincu que c'est la voie à suivre, et vous ne regarderez jamais en arrière.

Mais avec cela, faites peut-être une petite pause, et vous me rejoignez dans quelques minutes pour passer en revue la partie finale six de Hello Nextflow, qui concerne entièrement la configuration.

Merci beaucoup.
