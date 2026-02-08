# Partie 5 : Hello Containers - Transcription vidéo

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Traduction assistée par IA - [en savoir plus et suggérer des améliorations](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

<div class="video-wrapper">
  <iframe width="560" height="315" src="https://www.youtube.com/embed/Xqr--bKEN9U?si=QinuAnFwFj-Z8CrO&amp;list=PLPZ8WHdZGxmWKozQuzr27jyMGqp9kElVK&amp;cc_load_policy=1" title="YouTube video player" frameborder="0" allow="accelerometer; autoplay; clipboard-write; encrypted-media; gyroscope; picture-in-picture; web-share" referrerpolicy="strict-origin-when-cross-origin" allowfullscreen></iframe>
</div>

!!! note "Notes importantes"

    Cette page présente uniquement la transcription. Pour les instructions complètes étape par étape, retournez au [matériel de cours](../05_hello_containers.md).

    Les numéros de section indiqués dans la transcription sont fournis à titre indicatif uniquement et peuvent ne pas inclure tous les numéros de section du matériel.

## Bienvenue et contexte

Bonjour et bienvenue pour cette nouvelle partie de Hello Nextflow. Il s'agit de la partie cinq intitulée Hello Containers. Dans cette partie du cours, nous allons parler de comment encapsuler les exigences logicielles d'un pipeline afin que les personnes qui l'exécutent n'aient pas à se soucier de l'installation des logiciels.

Si vous travaillez en bio-informatique depuis aussi longtemps que moi, vous vous souvenez peut-être de ce que j'appelle souvent « le mauvais vieux temps », où quand vous vouliez exécuter le pipeline de quelqu'un d'autre ou reproduire son travail, vous passiez des heures ou des jours à essayer d'installer tous les différents outils logiciels qu'ils utilisaient, dans les mêmes versions, en essayant de les compiler sur votre machine, et c'était un cauchemar. C'était vraiment difficile.

Si vous travailliez sur un HPC, vous avez peut-être utilisé des modules d'environnement où les administrateurs système essayaient d'installer des logiciels pour vous, ce qui était correct, mais toujours imparfait.

Mais maintenant nous avons de meilleures façons de faire. Nextflow a un support intégré pour différentes technologies de conteneurs logiciels. Docker est la plus courante. C'est celle que nous allons utiliser aujourd'hui. Elle fonctionne bien dans Codespaces. Elle fonctionne bien sur votre ordinateur local et elle fonctionne bien dans le cloud.

Mais aussi Singularity ou Apptainer, qui sont très courants sur les systèmes HPC et fonctionnent effectivement exactement de la même manière. Ou Podman, Shifter, il y a plein d'autres qui sont tous très similaires.

La seule autre technologie, qui est un peu similaire mais pas tout à fait, que Nextflow supporte est Conda. Et Nextflow peut gérer des environnements Conda pour vous par processus, ce qui est bien mieux que de gérer vos propres environnements Conda. Et là encore, peut être livré avec un pipeline.

Nous allons commencer ce chapitre en parlant un peu des technologies de conteneurs et de Docker et de leur fonctionnement. Et nous allons faire la première moitié manuellement dans Docker pour que vous compreniez ce qui se passe sous le capot et comment cela fonctionne. Parce que c'est vraiment important de comprendre ce que fait Nextflow et comment comprendre ce que fait votre workflow quand il est exécuté.

Alors. Passons dans notre Codespaces. J'ai tout nettoyé à nouveau, mais si nous allons dans Hello Containers, vous devriez voir que tous nos scripts et tout sont là, identiques à la fin du chapitre sur les modules. Donc nous avons nos différents modules ici, que j'ai créés dans le répertoire modules.

Ils sont toujours là. Ils doivent être là pour que ça puisse s'exécuter. Et le workflow et la sortie sont tous identiques sauf que nous avons changé le chemin de publication de sortie vers Hello Containers, pour que vos fichiers se retrouvent dans ce répertoire.

Nous pouvons l'exécuter maintenant pour vérifier que ça marche si vous voulez, ou nous pouvons continuer avec le terminal.

## 1. Utiliser un conteneur « manuellement »

Nous allons utiliser Docker pour gérer nos conteneurs, et je peux vérifier qu'il est installé sur mon Codespaces en faisant « docker -v », ce qui me montre la version installée et tout, et que ça fonctionne correctement.

Maintenant les conteneurs et Docker ont deux concepts qui sont vraiment importants. L'un s'appelle image, et l'autre s'appelle conteneur. L'image est l'instantané, si vous voulez, de tout le système de fichiers que vous allez utiliser, et le conteneur est l'environnement d'exécution. Donc vous créez un conteneur en utilisant une image.

Une fois que vous êtes dans ce conteneur, il fonctionne généralement comme un système d'exploitation complet. Il est coupé du monde extérieur. Il est séparé de tout le reste, et c'est une bonne chose. C'est comme ça que nous obtenons une si bonne reproductibilité avec Nextflow.

Parce que pour les tâches exécutées à l'intérieur d'un conteneur, elles ne sont pas contaminées par des fichiers de configuration sur votre système local. Aucune autre influence externe, elles s'exécutent dans leur propre petit bac à sable. Les fichiers sont ensuite produits d'une manière très, très reproductible parce que vous utilisez les mêmes bibliothèques sous-jacentes, toutes les mêmes dépendances, exactement le même logiciel pour chaque personne s'exécutant sur chaque environnement informatique différent. Ce qui franchement je pense est fantastique et incroyable que ça marche. Et même, encore aujourd'hui ça me stupéfie toujours que ce soit possible.

## 1.1. Télécharger l'image du conteneur

Donc nous allons essayer d'utiliser des images Docker et Docker, quand vous l'exécutez sur votre système, a un registre Docker sur votre ordinateur, ou dans ce cas, dans le codespace, qui garde trace de toutes les différentes images qui ont été téléchargées et utilisées dans le passé, et les différentes couches dont elles sont composées.

Nous pouvons voir quelles images nous avons localement avec Docker en faisant « docker image ls ». Et dans ce cas vous pouvez voir qu'il y a plein d'images Docker ici, qui ont toutes à voir avec la configuration de ce Codespaces. Tout à voir avec les dev containers et autres. Donc vous n'avez pas besoin de trop vous en soucier, mais au fur et à mesure que nous ajoutons plus d'images et les téléchargeons, au fil de ce cours, vous pouvez vérifier cette liste et vous verrez que le registre local garde trace de toutes ces choses que nous avons téléchargées.

Mais nous allons en récupérer une nouvelle en faisant « docker pull ». Et ça dit à Docker d'aller chercher une nouvelle image depuis le web.

Nous mettons ensuite l'URI pour ce conteneur. Maintenant ça pourrait être une image Docker que vous avez construite localement et ensuite poussée sur internet. Ça pourrait être une image que quelqu'un d'autre a faite. Il y a beaucoup, beaucoup, beaucoup de façons différentes de faire des images Docker, mais sans doute l'une des façons les plus simples est d'externaliser ça, et de faire en sorte que quelqu'un d'autre le fasse pour vous.

Et ce que nous allons utiliser dans ce tutoriel est un service de Seqera appelé Seqera Containers.

Maintenant, Seqera Containers est totalement gratuit, et il utilise un logiciel open source que nous développons appelé Wave, qui a été construit pour gérer les conteneurs de manière complémentaire à Nextflow. Et il gère beaucoup des cas d'usage courants que nous rencontrons avec Nextflow.

Il est très courant que le logiciel dont nous avons besoin soit empaqueté dans Conda, dans les canaux Bioconda, ou conda-forge ou d'autres canaux plus spécifiques à un domaine. Et Wave et Seqera Containers sont vraiment bons pour construire des images à partir de ça.

Donc je peux aller sur cette interface web et nous allons jouer avec le package appelé « cowpy ». Donc je tape le nom du package que je veux. Il cherche, il l'a trouvé sur l'index de packages Python, donc je peux utiliser ça. Ou si j'attends un peu plus longtemps, il cherche dans bioconda et conda-forge. Et vous pouvez voir, je peux spécifier n'importe quel canal Conda ici. Donc si vous voulez trouver un canal Nvidia ou autre chose, ça devrait aussi fonctionner.

Et ensuite je peux spécifier si je veux qu'il construise une image Docker pour moi ou une image Singularity et aussi quelle architecture de CPU je veux utiliser. Donc amd64 ou arm64.

Et une fois que les résultats de bioconda sont listés, je peux maintenant voir toutes les différentes versions qui sont aussi disponibles. Je vais ajouter ça. Et maintenant je pourrais continuer à chercher et obtenir plus de packages de Conda si je veux et composer ce conteneur comme je le souhaite, mais je veux juste celui-là. Donc je vais cliquer sur Get Container.

Maintenant, quelqu'un d'autre a déjà demandé le même conteneur avant et il est retourné depuis un registre, donc nous l'obtenons immédiatement. Mais si personne d'autre n'avait jamais demandé ce package logiciel ou cette combinaison de packages logiciels, Wave et Seqera Containers le construiraient à la volée pour nous.

Nous pouvons copier cette URL et nous pouvons aussi voir les détails de construction. Et cela nous montre ce que le service a fait en arrière-plan. Il a créé un fichier d'environnement Conda. Un fichier Docker, et ensuite c'est ça, qui exécute le processus de construction Docker. Il a aussi exécuté un scan, un scan de sécurité, donc vous pouvez voir tous les CVE. Et il vous dit quand cela a été créé.

Wave et Seqera Containers peuvent faire beaucoup plus que ça, mais c'est un cas d'usage simple, qui est le plus courant. Et je dois dire que ces images sont hébergées pendant au moins cinq ans. Donc vous pouvez intégrer ces URLs dans vos pipelines et savoir qu'elles ne vont pas disparaître de sitôt.

Donc j'ai mon URL pour mon image Docker pour cowpy.

Je peux maintenant faire « docker pull » cette URL, et ça va récupérer toutes les différentes couches et télécharger cette image pour qu'elle soit disponible pour moi localement.

## 1.2. Utiliser le conteneur pour exécuter cowpy comme commande ponctuelle

D'accord, maintenant essayons de l'utiliser réellement. Donc maintenant je vais utiliser une commande « docker run » au lieu de « docker pull », et je vais utiliser le flag « --rm », qui dit juste à Docker de fermer ce conteneur une fois qu'il a fini de faire ce que je lui ai demandé. Et ensuite je mets l'identifiant pour le conteneur, qui est juste un URI.

Et ensuite à la fin, je spécifie la commande que je veux que Docker exécute à l'intérieur du conteneur généré à partir de cette image. Je vais juste dire cowpy, qui est le nom de l'outil qui est installé depuis Conda Forge, qui est disponible à l'intérieur de l'image.

Je vais appuyer sur entrée et voilà. Nous avons exécuté cowpy sur un système. Nous avons une petite vache qui nous donne des informations.

Maintenant notez que cowpy n'est pas installé sur mon système local. Donc si je l'exécute juste sans tout le truc Docker, ça dit, commande non trouvée. Donc ça a téléchargé une image. Ça a créé un conteneur en utilisant Docker, et ensuite c'est allé dans ce conteneur et a exécuté cette commande pour nous et nous a donné la sortie dans notre terminal. Très, très cool.

## 1.3. Utiliser le conteneur pour exécuter cowpy de manière interactive

D'accord, nous allons aller un peu plus loin maintenant et exécuter ce conteneur de manière interactive et jeter un coup d'œil, pour que nous puissions voir ce qui se passe à l'intérieur du conteneur.

Donc si je remonte et que je reprends ma commande run et je vais me débarrasser de cowpy à la fin là, parce que je ne veux pas vraiment exécuter cowpy. Je veux exécuter un terminal Bash.

Et ensuite je vais revenir ici et je vais faire « -it », qui signifie Interactive et Terminal ou TTY, et je vais appuyer sur entrée.

Et maintenant vous pouvez voir que le prompt, la partie avant que je tape, a changé. C'était le prompt Codespaces où il indiquait le répertoire, et maintenant il dit base et root et tmp. Donc je suis maintenant à l'intérieur du conteneur, et si je fais « ls », vous verrez que les fichiers que je vois dans ce répertoire sont différents des fichiers que j'ai dans mon workspace.

Et en fait, je ne peux voir aucun des fichiers de mon workspace codespaces local ou de mon disque local à l'intérieur du conteneur Docker. Le runtime du conteneur Docker, est complètement isolé et il ne peut pas écrire ou lire des fichiers depuis un système de fichiers hôte extérieur.

Je peux, cependant, voir le logiciel qui est installé à l'intérieur du conteneur et l'exécuter. Donc je peux exécuter cowpy et nous pouvons voir un peu plus sur comment utiliser cowpy. Ici je peux faire « cowpy 'Hello World' » et ça lui dit de mettre ma citation à l'intérieur d'une petite bulle de dialogue. Et vous pouvez aussi exécuter différents types de vaches, donc ça ne doit pas être une vache. Vous pouvez faire un « -c ». Et je suis en Suède, donc je vais choisir un élan. Très bien. Je lui ai donné des bois.

Et il y en a tout un tas de différents avec lesquels vous pouvez jouer, que vous pouvez voir décrits dans la documentation de formation.

## 1.3.4. Monter des données dans le conteneur

D'accord. Ce serait bien si nous pouvions exécuter cowpy sur les fichiers de notre système de fichiers.

Bien sûr, ce n'est pas très utile d'avoir juste le conteneur et aucun accès à quoi que ce soit du tout. Ça peut être sûr et reproductible, mais ce n'est pas très utile.

Donc comment faisons-nous ça ? Je vais sortir de ce conteneur Docker en tapant exit, et vous pouvez voir que le prompt nous dit que nous sommes maintenant de retour dans notre Codespaces régulier.

Et je vais exécuter la même commande à nouveau. Mais cette fois je vais ajouter des flags supplémentaires ici. Et celui qui est important est « -v », qui signifie monter un volume, qui est comme fondamentalement une partie, partie d'un espace disque.

Le « -v » prend deux parties : il y a comme une chaîne et ensuite deux points et une chaîne. Et la première partie est le système de fichiers local, qui devrait être monté dans le conteneur. Et ensuite la deuxième partie est où ça devrait se retrouver à l'intérieur du conteneur.

Maintenant je veux juste charger tout mon système de fichiers local ici. Donc « . » est le répertoire de travail actuel. Donc je vais juste faire « . » et ensuite « : », et ensuite nous allons mettre ça dans un nouveau répertoire à l'intérieur du conteneur appelé « my_project ». Ça pourrait vraiment s'appeler n'importe comment.

Et ensuite je vais exécuter à nouveau.

Dans le répertoire de travail où je suis placé, qui est /tmp, les fichiers ne sont pas là. Mais si je fais « ls my_project », nous y voilà : tous les mêmes fichiers que nous avions localement sur notre Codespaces sont maintenant disponibles à l'intérieur du conteneur à ce chemin.

C'est un accès en lecture et écriture donc je peux créer de nouveaux fichiers dans ce répertoire et ils apparaîtront sur mon système de fichiers hôte. Ce répertoire particulier, se comporte alors exactement comme si j'étais à l'extérieur du conteneur donc je peux maintenant lire et écrire et faire des choses.

## 1.3.5. Utiliser les données montées

D'accord, prouvons juste que nous pouvons faire ça. Je fais « cat /my_project/data/greetings.csv ». Si vous vous souvenez, le contenu de ce fichier ressemble à ça. Je peux maintenant rediriger ça vers cowpy et la vache va afficher les différentes sorties de ce fichier dans sa petite bulle de dialogue, ce qui est plutôt amusant.

Donc vous pouvez voir, nous pouvons maintenant utiliser le logiciel dans le conteneur pour interagir avec les fichiers sur notre système hôte.

D'accord, sortons et nous continuerons avec le reste du matériel de formation.

## 2. Utiliser des conteneurs dans Nextflow

Donc c'est vraiment cool d'utiliser des conteneurs. J'espère que ça a du sens. Et vous pouvez voir la valeur de ces conteneurs et pourquoi c'est utile pour exécuter des logiciels d'analyse.

Mais comment faisons-nous tout ce même processus à l'intérieur de Nextflow ? Nous ne voulons pas exécuter plein de commandes Docker nous-mêmes. Nous voulons juste laisser Nextflow gérer tout ça pour nous.

Alors travaillons là-dessus. Nous allons ajouter un nouveau processus à notre pipeline, pour exécuter cowpy. D'accord, alors créons un nouveau module pour notre nouveau processus. Donc allons dans modules, appelons-le cowPy.nf, et ensuite je vais copier le code du matériel de formation ici.

Mais vous pouvez voir que le processus est très simple. Il ressemble beaucoup à ceux que nous avons faits jusqu'à présent, nous avons un bloc input avec un path, qui est notre fichier d'entrée, et aussi une value ici pour que ce soit un caractère, donc nous pourrions utiliser un élan à nouveau si nous voulons.

Et ensuite une output, qui est un seul fichier ici, un path et ensuite un script. Et nous faisons la même chose que nous avons faite de manière interactive à l'intérieur du conteneur : nous faisons « cat » pour lire le fichier d'entrée. Nous redirigeons ce contenu vers cowpy. Nous choisissons un caractère spécifique basé sur cette entrée, nous écrivons dans un fichier de sortie appelé cowpy input file, qui est ensuite émis vers la sortie.

Très bien. Incluons ça. Donc include \{ cowpy \} from « ./modules/cowpy.nf », l'ai-je appelé cowpy ? Oui.

Et ensuite appelons notre nouveau processus en bas ici dans le bloc main du workflow. Donc exécutons cowpy. Et nous allons prendre notre nouveau processus cowpy et nous allons dire collectGreetings.out.

Et ensuite si vous vous souvenez, il y avait deux sorties pour ce module. Une appelée outfile et une appelée report. L'extension VS Code nous suggère automatiquement celles-ci et nous voulons .outfile.

Vous pouvez toujours aller dans ce processus ici. Vous survolez dessus et ça devrait vous montrer rapidement quelles étaient les sorties. Et nous pouvons aussi faire commande clic dessus et ça ouvrira le fichier de module si vous voulez voir plus en détail.

Donc voilà. C'est l'outfile là, et c'est le path. Donc ça sera maintenant le fichier d'entrée pour notre processus cowpy. Fantastique.

Maintenant si vous vous souvenez, un processus cowpy a deux entrées. Nous avions aussi le canal de valeur pour le caractère. Donc nous pouvons ajouter « params.character » ici. J'aurais pu coder ça en dur si je voulais, mais faisons-en une option CLI pour que nous puissions faire dash, dash character.

Bien. Je dois maintenant définir le paramètre d'entrée que nous venons d'appeler et lui donner une valeur par défaut. Donc character, String. Et j'aime l'élan, donc je vais le définir à moose par défaut.

Bien, essayons de l'exécuter. Donc si je fais Nextflow run hello containers, nous verrons ce qui se passe.

J'aurais pu utiliser dash resume si j'avais les anciens répertoires de travail qui traînaient. Et encore, ces premiers processus auraient été mis en cache et ça aurait été un peu plus rapide, mais ça devrait être fondamentalement pareil.

Maintenant nous pouvons voir tout de suite qu'il a lancé une erreur quand il est arrivé à notre nouveau processus, il nous dit ici qu'il y a eu une erreur en exécutant le processus cowpy et il s'est terminé avec un statut de sortie 127. C'est la commande qu'il a essayé d'exécuter. Ça a l'air correct, ça ressemble à ce que nous attendions. Ça prend ce nom de fichier de sortie, qui a l'air à peu près correct, ça l'exécute avec un caractère moose et essaie de le sauvegarder.

Mais vous pouvez voir que l'erreur de commande ici dit que la commande cowpy n'a pas été trouvée. Et ça a du sens parce que nous n'avons pas vraiment dit à Nextflow d'utiliser un conteneur encore. Nous lui avons juste donné la commande cowpy. Et comme je l'ai dit avant, cowpy n'est pas installé sur notre système local. Donc quand il a essayé de l'exécuter, ça a échoué.

## 2.3.1. Spécifier un conteneur pour cowpy

Nous devons dire à Nextflow qu'il y a un conteneur disponible et qu'il peut l'utiliser. Donc comment faisons-nous ça ?

Si nous allons dans notre module ici, nous allons ajouter une nouvelle déclaration en haut appelée « container ». Et nous allons ensuite définir ça à une chaîne.

Maintenant, si vous vous souvenez, dans Seqera Containers, je peux copier cette URL et je la dépose juste dans des guillemets ici.

Maintenant revenons et essayons de l'exécuter à nouveau.

Voyons si ça marche cette fois.

Malheureusement, ça échoue exactement de la même manière, même si maintenant nous avons défini un conteneur pour que le processus s'exécute. Donc pour utiliser notre image Docker, nous devons dire à Nextflow d'activer l'utilisation de Docker quand nous exécutons le workflow.

Et nous allons faire ça en créant un nouveau fichier de configuration. Donc je vais dire touch nextflow.config.

C'est un nom de fichier spécial où s'il est dans le répertoire de travail pendant que je lance le pipeline, il sera chargé automatiquement. Donc si je vais dans ce fichier nextflow.config, vous pouvez voir qu'il existe déjà en fait, ce que j'avais oublié. Et nous avons docker.enabled dedans déjà, mais il est défini à false, qui est la valeur par défaut.

Donc si je change ça à égale True à la place, docker.enabled. Et il y a de la documentation de référence pour tous ces scopes de configuration dans la documentation Nextflow. Et aussi vous pouvez voir que quand je survole avec l'extension VS Code, ça charge la documentation spécifique à ça et me dit ce que ça signifie et comment le définir.

Donc maintenant nous l'avons défini à true, et si j'exécute Nextflow à nouveau, Nextflow saura maintenant aller chercher cette image Docker pour nous si nous ne l'avons pas déjà localement, et ensuite exécuter ce processus avec cet environnement de conteneur.

Et donc nous pouvons voir qu'il s'est exécuté avec succès et nous avons une petite coche à côté d'un cowpy. Fantastique. Si je monte et regarde dans le répertoire results, le fichier n'est pas encore là. Et c'est parce que nous devons encore publier ce fichier de sortie exactement comme tous les autres.

Donc nous allons dans le bloc published dans le workflow, disons mycowpy égale cowpy.out.

Et ensuite en bas ici dans le bloc output, mycowpy, accolades path. Oups. Hello containers. Mode, copy.

Si j'exécute à nouveau maintenant, ça devrait s'exécuter exactement de la même manière. J'aurais pu utiliser dash resume et j'oublie à chaque fois. Et ensuite je monte et maintenant nous avons un nouveau fichier créé appelé cowpy-COLLECTED, et voilà mon élan disant BONJOUR, HELLO, HOLÀ Fantastique.

Maintenant bien sûr je pourrais aussi passer maintenant « --character ». Quelles sont les différentes options ? Je pense qu'il y a une dinde ? Donc je peux utiliser character Turkey. Ça va s'exécuter exactement de la même manière. J'ai manqué une autre occasion d'utiliser dash resume, et maintenant si nous chargeons notre fichier et maintenant nous avons une dinde. Fantastique.

## 2.3.4. Inspecter comment Nextflow a lancé la tâche conteneurisée

D'accord. Dernière petite chose. Exécutons rapidement cette commande à nouveau, resume cette fois, et jetons un coup d'œil rapide dans le répertoire de travail pour voir ce que fait Nextflow sous le capot pour faire tout ça fonctionner pour nous.

Cette fois c'est super rapide, allons dans ce répertoire de travail, cd work/. Maintenant si vous vous souvenez nous avons plein de fichiers point ici et celui qui nous intéresse dans ce cas est celui que j'ai dit que nous n'avons presque jamais besoin de regarder, appelé .command.run.

Si je fais code dot command run, ça va l'ouvrir dans l'éditeur. Et je peux chercher dans ce fichier et si je descends je devrais voir Docker run. Et vous pouvez voir que Nextflow fait la commande docker run pour nous, quand Docker est activé dans une configuration. Il a tout un tas de flags et de choses différents ici, mais vous pouvez voir le flag « -v » que nous avons utilisé nous-mêmes quand nous exécutions. Et vous pouvez voir qu'il monte le répertoire de workspace local dans le conteneur, pour que le conteneur puisse accéder à nos fichiers d'entrée et sauvegarder les sorties. Et ensuite à la fin, il exécute aussi .command.sh, qui est le script généré, qui contient la commande cowpy.

Et donc vous pouvez voir que Nextflow prend la logique du workflow, qui est le truc qui nous importe vraiment, qui est spécifique à notre analyse, et il fait tout le truc intelligent en coulisses pour faire fonctionner Docker sur notre système.

Et il fait ça d'une manière vraiment portable pour qu'un·e utilisateur·trice final·e du pipeline puisse changer la technologie qu'il·elle utilise : Docker, Singularity, Apptainer, Conda. Ça n'a pas vraiment d'importance pour la logique du pipeline, mais Nextflow gérera tous les besoins d'infrastructure sous-jacents, pour qu'il s'exécute n'importe où.

Et c'est vraiment le superpouvoir de Nextflow. C'est la reproductibilité et la portabilité. Et avec Nextflow vous pouvez réellement partager votre workflow et d'autres personnes peuvent l'exécuter sur leurs systèmes et ça fonctionnera tout simplement.

C'est une chose vraiment, vraiment difficile à faire, et maintenant vous savez comment le faire aussi avec vos workflows.

D'accord, c'est tout pour ce chapitre. Si vous descendez à la fin du cours, vous trouverez un quiz encore sur les conteneurs. J'espère que tout avait du sens. C'est une façon vraiment cool de travailler avec l'analyse. Et si vous débutez avec les conteneurs, j'espère vous avoir convaincu·e que c'est la façon de faire, et vous ne regarderez jamais en arrière.

Mais avec ça, prenez une petite pause peut-être, et vous me rejoindrez dans quelques minutes pour parcourir la partie finale six de Hello Nextflow, qui concerne entièrement la configuration.

Merci beaucoup.
