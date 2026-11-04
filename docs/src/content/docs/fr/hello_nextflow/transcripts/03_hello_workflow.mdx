# Partie 3 : Hello Workflow - Transcription vidéo

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Traduction assistée par IA - [en savoir plus et suggérer des améliorations](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

<div class="video-wrapper">
  <iframe width="560" height="315" src="https://www.youtube.com/embed/_aO56V3iXGI?si=Irl9nAQniDyICp2b&amp;list=PLPZ8WHdZGxmWKozQuzr27jyMGqp9kElVK&amp;cc_load_policy=1" title="YouTube video player" frameborder="0" allow="accelerometer; autoplay; clipboard-write; encrypted-media; gyroscope; picture-in-picture; web-share" referrerpolicy="strict-origin-when-cross-origin" allowfullscreen></iframe>
</div>

!!!note "Notes importantes"

    Cette page présente uniquement la transcription. Pour des instructions détaillées, retournez au [matériel du cours](../03_hello_workflow.md).

    Les numéros de section affichés dans la transcription sont fournis à titre indicatif uniquement et peuvent ne pas inclure tous les numéros de section des supports.

## Bienvenue et récapitulatif

Bonjour et bienvenue dans la partie trois de Hello Nextflow. Cette partie s'appelle Hello Workflow, et c'est dans cette partie du cours que nous allons vraiment justifier le nom pipeline ou workflow.

Nous allons prendre notre script de pipeline simple jusqu'à présent avec son unique processus, et nous allons commencer à ajouter des processus supplémentaires pour voir comment Nextflow gère cette orchestration et le flux de données à travers le pipeline.

Retournons à nos code spaces. Vous verrez que j'ai supprimé tous mes répertoires .nextflow\* et les répertoires work pour garder l'espace propre. Ne vous inquiétez pas si vous avez encore ces fichiers des parties précédentes du cours.

Nous allons travailler à partir d'un fichier appelé hello-workflow.nf. Comme précédemment, cela représente essentiellement le script que nous avons construit jusqu'à présent, et nous donne un point de départ propre. Et encore une fois, dans la sortie, nous pouvons voir que le chemin est maintenant hello_workflow. Donc les fichiers publiés devraient aller dans un sous-répertoire différent dans votre dossier results.

Pour récapituler où nous en sommes jusqu'à présent, nous avons un seul processus ici, avec une entrée greeting, une sortie greeting file. Et puis le simple script Bash, qui fait juste une commande echo vers un fichier.

Nous avons une seule entrée de workflow, le bloc params ici, où nous disons qu'il attend un chemin, et la valeur par défaut est data/greetings.csv, qui est ce fichier ici.

Ensuite, dans le workflow lui-même, nous avons un bloc main. Nous créons un canal. Nous analysons le CSV en lignes puis prenons le premier élément de chaque tableau, et nous passons ce canal dans ce processus, qui génère ensuite trois tâches, et nous publions depuis le workflow les sorties de ce processus.

Et enfin, dans le bloc output, nous disons à Nextflow de publier ces fichiers depuis ce canal vers le répertoire appelé hello_workflow. Et de copier ces fichiers plutôt que de créer des liens symboliques.

## 1. Ajouter une deuxième étape au workflow

D'accord, dans cette partie nous allons ajouter un deuxième processus à notre workflow. Nous allons prendre les sorties du processus sayHello, et les traiter dans une deuxième étape, qui va convertir toutes les lettres dans ces fichiers en majuscules, convertToUppercase.

C'est juste un exemple simple, c'est encore du traitement de chaînes simple, mais cela vous montre comment nous pouvons prendre la logique dans le workflow.

Nous allons utiliser une commande bash appelée "tr" pour cela, qui signifie translate. C'est une commande Unix qui existe depuis toujours. Si vous ne la connaissez pas, je ne vous en veux pas. Je ne pense pas l'avoir jamais utilisée avant la formation, mais vous pouvez l'essayer très rapidement dans le terminal. Si je fais "echo 'hello world'" et que je pipe vers 'tr' et ensuite entre guillemets vous dites plage de caractères, donc A à Z, minuscules, et ensuite vous voulez faire A à Z majuscules. Et ça dit simplement, traduire ces lettres en ces lettres.

Et quand j'appuie sur entrée, vous pouvez voir que tout est maintenant en majuscules. Très pratique si vous aimez crier sur les gens.

Donc c'est un style très simple de commande bash que nous allons utiliser dans notre deuxième processus.

## 1.2. Écrire l'étape de mise en majuscules comme processus Nextflow

Donc si je retourne à mon script, je vais un peu tricher et simplement copier le code de la documentation de formation. Mais vous pouvez voir exactement ce qui se passe.

Nous avons un nouveau processus ici. Celui-ci nous l'avons appelé convertToUpper, mais nous pourrions l'appeler comme nous voulons.

Nous avons une seule entrée path, comme nous l'avons fait précédemment. Ce n'est pas un canal de valeur, c'est un canal de chemin. Et puis une seule sortie.

Dans le bloc script, nous faisons "cat" sur le fichier d'entrée. Et nous pouvons mettre cela entre accolades si nous voulons. et qui prend cette variable. Et nous exécutons cette même commande bash dans le pipe et nous écrivons les résultats dans un fichier avec ce nom de fichier, et c'est capturé par le chemin de sortie.

Nous devons maintenant faire quelque chose avec ce nouveau processus. Donc nous allons descendre au workflow où nous construisons la logique différente d'un workflow, et après ce premier processus, nous allons exécuter notre deuxième processus. Donc convertToUpper est le nom du processus ici.

Il prend une entrée donc nous ne pouvons pas juste l'appeler tout seul. Nous voulons traiter la sortie du premier processus. Donc tout comme nous l'avons fait avec ceci, sayHello out où nous publions ces résultats. Nous voulons utiliser ces mêmes résultats ici comme entrée, donc nous pouvons les copier et les mettre là.

Nous voulons le processus sayHello ".out", et Nextflow sait que cela signifie un simple enregistrement de sortie unique ici, qui est ce fichier. Donc cela sera ensuite passé comme entrée à un deuxième processus.

## 1.5. Configurer la publication de sortie du workflow

D'accord. Et finalement, pour que nous sauvegardions réellement les résultats de ce deuxième processus, nous devons aussi les publier depuis le workflow, puis les définir dans le bloc output, même syntaxe qu'avant. Donc nous pouvons copier ceci et dire second outputs, ou comme vous voulez l'appeler.

Prendre le nom du processus qui nous intéresse, convertToUpper out, et ensuite ici dans le bloc output. Ajouter ceci et nous pourrions faire les mêmes attributs ici. Donc nous voulons aussi ces fichiers dans le sous-répertoire Hello Workflow, et nous voulons aussi les copier.

Super. Essayons de l'exécuter. Donc si j'ouvre le terminal et je fais "nextflow run hello-workflow.nf", et nous verrons ce qu'il fait. Voyons si ça a l'air différent des parties précédentes.

Donc il lance Nextflow. Dans la documentation, il dit de faire ceci avec "-resume", mais j'ai supprimé tout mon répertoire work, donc ça n'aurait fait aucune différence ici. Mais si vous l'aviez, alors cela fonctionnera aussi bien.

Et ça ressemble presque exactement à la même chose. Mais vous pouvez voir maintenant qu'il y a une deuxième ligne de sortie ici, où vous pouvez voir le nom du deuxième processus que nous venons d'ajouter. Et effectivement, vous pouvez voir qu'il s'est exécuté trois fois avec succès.

Excellent. Si j'avais mes répertoires work précédents et que j'avais fait ceci avec "-resume", ceux-ci auraient été mis en cache juste la première étape dans le pipeline. Parce que ces sorties étaient exactement les mêmes, donc Nextflow aurait su les réutiliser à nouveau.

Et donc vous pouvez voir comment vous pouvez utiliser -resume pour construire itérativement votre workflow, étape par étape, si vous en avez besoin.

D'accord, jetons un œil dans le répertoire results ici et voyons si ça a fonctionné. Nous pouvons voir que nous avons quelques fichiers supplémentaires ici. Nous avons nos fichiers originaux comme nous les avions avant du premier processus. Et effectivement, nous avons nos fichiers upper et les lettres sont toutes en majuscules, donc ça a fonctionné. C'est vraiment agréable à voir.

Il est également intéressant de regarder à l'intérieur de ces répertoires work. Comme avant le hash ici correspond aux répertoires work. Donc si je regarde dans "ls work", puis développe cela, nous verrons les différents fichiers ici.

Nous voyons le fichier de sortie du premier processus, qui a été importé ici comme entrée. Et nous pouvons voir le nouveau fichier de sortie qui a été généré.

Maintenant si je fais ceci avec "-la" pour lister et montrer tous les fichiers, nous verrons quelques choses de plus. Premièrement, vous verrez que ce fichier est en fait un lien symbolique vers le premier processus. C'est fondamentalement toujours un lien symbolique si possible, pour économiser l'espace disque. Nous ne publions pas les fichiers ici et il référence simplement ce fichier d'une première tâche dans une deuxième tâche afin que tout soit encapsulé dans ce répertoire de travail unique, et sûr et isolé de tout le reste.

Et cela doit être là parce que si nous regardons le fichier .command.sh, donc si je fais "cat work/b8/56\*", vous pouvez voir que les parties de fichier ici sont relatives, donc il fait cat de ce fichier d'entrée, qui a été lié symboliquement dans le même répertoire de travail.

Donc c'est ainsi que chaque répertoire work ressemblera. Quand vous le regardez dans Nextflow, vous aurez tous les fichiers d'entrée là mis en scène dans ce répertoire work. Et ensuite vous aurez aussi tous les fichiers de sortie qui ont été créés. Donc c'est super. Cela ressemble à ce que nous attendions.

## 2.1. Définir la commande de collecte et la tester dans le terminal

D'accord, revenons à notre workflow. Quelle est la prochaine étape que nous voulons faire ?

Nous avons maintenant deux processus et ils prennent ce fichier CSV unique, l'analysent et le divisent. Et ensuite nous avons trois tâches pour chacun de ces processus et Nextflow gère la parallélisation de tout cela, donc tout s'exécute côte à côte quand c'est possible.

Cette façon de diviser le travail pour exécuter les choses en parallèle est très courante. Et l'inverse de cela est ensuite de tout rassembler. Donc c'est ce que nous allons faire avec notre processus final dans le workflow, nous en aurons un troisième ici, qui prend ces trois sorties différentes et les combine toutes dans un seul fichier.

Nous pouvons faire cela assez simplement dans un terminal, juste pour avoir une idée de ce à quoi cela ressemblera.

Si je vais dans le dossier results. Donc, "cd results/hello_workflow/", et nous avons tous les fichiers UPPER ici. Je peux juste utiliser "cat", que nous utilisons pour imprimer le contenu de ce fichier, et vous pouvez donner plusieurs fichiers à "cat" et il les lira l'un après l'autre.

Donc je peux dire "UPPER-\*", ce qui me donne la même liste de trois noms de fichiers avec l'expansion Bash. Et je peux dire combined.txt. Je pense que dans la documentation, elle liste les noms de fichiers exacts, mais ça fait la même chose.

Maintenant, si j'utilise "cat combined.txt", nous pouvons voir que nous avons le contenu des trois fichiers.

Donc c'est fondamentalement tout ce que ce processus va faire, nous allons essayer de lui donner tous les différents fichiers de sortie d'un processus précédent dans une seule tâche de processus, et ensuite nous allons les "cat" ensemble et sauvegarder le fichier de sortie.

## 2.2. Créer un nouveau processus pour faire l'étape de collecte

D'accord, donc ajoutons notre nouveau processus. Je vais coller ceci depuis les supports de formation, et vous pouvez voir que ça nous a laissé un peu d'exercice pour le lecteur avec ces points d'interrogation. Mais vous pouvez voir le schéma général du processus est fondamentalement ce que nous venons de faire dans le terminal, où nous faisons "cat" d'un ensemble de fichiers d'entrée et l'écrivons dans un fichier de sortie ici appelé collected, et ensuite la sortie attend ce chemin unique à nouveau.

Donc nous avons besoin d'une sorte d'entrée ici et ce sera un ensemble de chemins. Donc encore une fois, nous définissons un canal d'entrée path et appelons-le input_files. Maintenant, cela nous a précédemment donné un chemin unique ici, mais un chemin peut aussi avoir plusieurs fichiers ici, même si c'est toujours une déclaration unique.

Je vais copier cela ici parce que nous voulons "cat" ces fichiers. Et vous pourriez penser que nous avons des problèmes ici avec l'impression d'un tableau ou des choses comme ça, mais Nextflow est généralement assez sensé quand il s'agit de cela. Et s'il reçoit un canal avec plusieurs fichiers dedans comme ceci, il va les mettre tous ensemble avec des séparateurs d'espace. Donc cela nous donnera la syntaxe correcte.

C'est super. Donc maintenant câblons notre nouveau processus. Je descends au workflow. Je vais dire combine the outputs, le nouveau nom de processus, et exactement comme avant. Je vais prendre ce processus précédent, convertToUpper et faire ".out".

Super. Essayons-le et voyons si ça fonctionne dans le terminal. Si je remonte juste de quelques répertoires et puis réexécute la commande Nextflow, et nous verrons ce qui se passe.

Donc le workflow s'est lancé et maintenant vous pouvez voir que nous avons trois noms de processus différents, ce qui est super. Les deux premiers ressemblent aux mêmes qu'avant, et le troisième nouveau s'exécute, ce qui est bien.

Cependant, il y a quelque chose d'un peu étrange ici. Nous voulions combiner ces fichiers de sortie dans un seul fichier, et pourtant ce processus que nous pouvons voir s'est exécuté trois fois, pas une fois.

Effectivement, si nous allons dans l'un de ces répertoires work. Et faisons "cat work/" "collected", alors nous verrons. Il n'y a qu'un seul mot ici, pas trois.

Et donc ce qui s'est passé est que Nextflow a continué cette parallélisation exactement comme il l'a fait dans les étapes précédentes. Et ce processus nous a donné un canal avec trois éléments, et ces trois éléments de canal ont été passés à notre processus en aval, qui a généré trois tâches de processus.

Il a essentiellement essayé de collecter trois fois séparément et à chaque fois il n'avait qu'un seul fichier, donc il a juste fait cat d'un fichier unique vers une sortie, et en fait, nous pouvons voir cela dans le fichier .command.sh également.

Si je fais .command.sh, nous pouvons voir qu'il n'a qu'un seul nom de fichier ici et seulement un seul fichier a été mis en scène dans ce répertoire de travail.

## 2.3. Ajouter l'étape de collecte au workflow

Donc d'une manière ou d'une autre nous devons dire à Nextflow de rassembler toutes ces sorties d'un processus précédent et de les donner à ce processus en aval comme un seul élément de canal, plutôt que trois.

Nous faisons cela avec un opérateur de canal appelé _collect_.

C'est un opérateur super utile, que vous verrez dans les pipelines Nextflow tout le temps. C'est un canal ici, ce canal de sortie, exactement comme celui que nous avons créé en haut. Et donc nous pouvons y ajouter des opérateurs de canal exactement comme nous l'avons fait avant. Nous pouvons juste faire point, et ensuite dans ce cas, collect, parenthèses.

Et c'est tout ce dont nous avons besoin. Cela va ensuite manipuler ce canal avant qu'il ne soit passé dans ce processus.

Si vous voulez voir ce qui lui arrive, nous pouvons aussi le voir ici. Donc ici, cela n'est pas lié à l'exécution de ce processus du tout, donc je pourrais le mettre à n'importe quel moment après l'exécution de ce processus. Mais nous prenons le même canal de sortie, et nous le regardons avec .view, et ensuite nous le regardons à nouveau avec .collect.view.

Et quand nous exécutons ceci, cela nous montrera les deux structures différentes de ce canal, avant et après collect. Donc essayons cela maintenant. D'accord, je viens juste de dézoomer un peu parce que certaines sorties sont assez longues, mais si j'exécute le pipeline, nous verrons si ça fonctionne.

J'espère qu'un troisième processus s'exécutera juste une fois, parce qu'il collecte les sorties et effectivement, vous pouvez voir collectGreetings comme un sur un. Donc ça a exécuté juste une tâche.

Et ensuite si nous regardons les déclarations view, nous avons trois déclarations view pour les trois éléments de avant, avec un chemin de fichier dans chacun.

Et ensuite après cette déclaration collect, cela s'est juste déclenché une fois parce qu'il y a un seul élément dans ce canal. Et maintenant nous avons cette liste de trois chemins de fichiers différents.

C'est exactement ce que nous espérions. Et vous pouvez voir, espérons-le, c'est fondamentalement l'inverse de cet opérateur "map" que nous avons fait pour aller des tableaux CSV en éléments de canal séparés. Maintenant nous prenons des éléments de canal séparés et les remettons dans un seul tableau.

Super, nous pouvons nettoyer ces déclarations view. Nous n'en avons plus besoin. Nous pouvons passer à l'étape suivante.

Avant d'aller plus loin, et avant que j'oublie, je vais ajouter une nouvelle déclaration publish ici. Third output. Vous pouvez appeler cela quelque chose de plus sémantique et descriptif dans votre workflow. Et ensuite je vais ajouter cela dans le bloc output à nouveau et dire path 'hello_workflow' mode 'copy'. Juste pour que le fichier de sortie généré par ce processus soit sauvegardé dans notre dossier results ici.

Juste pour vérifier rapidement que ça fonctionne. Devrait être un peu plus propre maintenant parce que nous n'avons pas ces déclarations view. Et, nous verrons si nous obtenons notre nouveau fichier de sortie ici. Une tâche sur une s'est exécutée, nous avons un nouveau fichier appelé collected, et maintenant nous avons les trois mots. Fantastique. Et ensuite ?

## 3. Passer des paramètres supplémentaires à un processus

D'accord. Ensuite nous allons regarder la gestion de plusieurs entrées dans un seul processus. Jusqu'à présent vous pouvez voir que tous nos processus prennent juste une chose comme entrée. Ils ont tous une seule ligne sous leur input.

Nous allons démontrer cela en permettant à Nextflow de spécifier un identifiant de lot différent afin que peut-être vous exécutiez ce workflow plusieurs fois et vous puissiez lui donner un ID de lot différent à chaque fois.

Je vais simplement ajouter une deuxième ligne dans l'entrée ici pour collectGreetings. Et je vais l'appeler "val", parce que c'est une chaîne. Maintenant c'est une valeur, pas un chemin, et je vais l'appeler "batch_name".

Ensuite je vais éditer le script ici pour utiliser cette variable, et je vais essayer de la mettre au même endroit que dans le matériel de formation. Donc je la mets au milieu de ce chemin de fichier COLLECTED-$\{batch_name\}-output.

Pas tout à fait fini. Rappelez-vous que nous devons dire à Nextflow ce que les noms de fichiers de sortie vont être. Donc nous devons aussi faire la même chose ici : COLLECTED-$\{batch_name\}-output.txt".

Fantastique. Nextflow obtient maintenant une deuxième entrée de variable et il l'interpole dans le script et la sortie.

Une dernière chose, nous devons maintenant trouver où cela est appelé, et nous devons passer la deuxième entrée au processus. C'est comme toute autre entrée dans une fonction dans n'importe quel autre langage.

Tout comme nous l'avons fait plus tôt dans la formation, je vais utiliser le "params" spécial ici, et nous allons l'appeler "params.batch" pour que nous puissions avoir une option CLI -- batch. Et maintenant vous pouvez voir que notre processus ici a deux entrées séparées juste séparées par des virgules, qui sont passées.

Il est vraiment important d'avoir le bon ordre, donc l'ordre des arguments ici pour channel puis le param doit correspondre. Le canal et le batch name là. C'est juste une correspondance positionnelle.

D'accord. Je peux exécuter ce pipeline maintenant tout de suite avec --batch, mais faisons d'abord la bonne chose et définissons-le dans l'entrée ici dans Params. Donc je vais l'ajouter à batch et ensuite nous allons dire que c'est une chaîne et donnons-lui une valeur par défaut. Donc appelons-le simplement batch. D'accord ? Maintenant essayons d'exécuter le workflow.

--batch Trio. Je pense que ça dit dans le matériel de formation, mais nous pourrions utiliser n'importe quelle chaîne que nous voulons là. Et espérons que nous verrons ce fichier de sortie de résultats apparaître ici.

Et effectivement, COLLECTED-trio-output - cela a bien fonctionné. Il a renommé notre fichier. Et vous pouvez imaginer maintenant que c'est utile parce que si j'exécute cela à nouveau avec un nom de lot différent, comme replicate_two, alors ça va nous donner un nom de lot différent ici.

Et et ça ne va alors pas écraser les fichiers de sortie dans ce cas. Donc c'est bien.

## 4. Ajouter une sortie à l'étape du collecteur

D'accord, donc nous avons maintenant plusieurs entrées à notre processus ici. Mais que se passe-t-il si nous voulons créer plusieurs sorties ? Notre exemple ici alors est que nous allons créer un rapport pour ce processus, disant simplement combien de fichiers ont été collectés.

Et nous ferons cela avec une commande echo ici. Donc nous pouvons dire echo. There were, je vais copier ceci depuis le matériel de formation, donc vous n'avez pas à me regarder le taper.

There were $\{count_greetings\} greetings in this batch, et sauvegarder cela dans un nouveau fichier maintenant appelé $\{batch_name\}, donc même variable, nous pouvons la réutiliser autant de fois que nous voulons, report.txt.

## 4.1.1. Compter le nombre de messages de bienvenue collectés

Nous devons réellement calculer cela d'une manière ou d'une autre. Nous pourrions faire cette logique dans le script Bash si nous voulions, en utilisant la logique Bash. Cependant, nous pouvons aussi juste faire du scripting directement dans le code Nextflow, tant que c'est dans le bloc script dans le processus et au-dessus de la section entre guillemets.

Tout ici ne sera pas inclus dans le script final rendu, et sera juste exécuté par Nextflow quand il rend une tâche.

Donc ici nous faisons juste de la logique. Nous créons une nouvelle variable appelée count_greetings. Nous prenons le canal input files ici, et nous appelons .size() dessus.

D'accord, cette fonction va me donner un nombre ici dans cette variable, et maintenant notre avertissement a disparu parce que cette variable est définie.

D'accord, donc nous créons ce deuxième fichier dans le répertoire work, mais nous devons dire à Nextflow de s'attendre à ce qu'il soit publié comme sortie de ce processus. Donc nous faisons cela par exactement la même syntaxe que nous avons faite pour le premier fichier.

Nous disons path parce que c'est, encore une fois, nous pourrions publier une variable ici si nous voulions avec "val", mais nous allons dire "path". Et ensuite le nom de fichier attendu. Remarquez qu'il n'est pas surligné ici. C'est parce que j'ai utilisé des guillemets simples. Je dois utiliser des guillemets doubles.

## 4.1.2. Émettre le fichier de rapport et nommer les sorties

D'accord, c'est super. Et nous pourrions maintenant commencer à accéder à ces sorties ici juste comme je l'ai fait ici. Mais maintenant c'est un tableau d'objets différents, donc je pourrais faire collectGreetings.out[0] pour obtenir le premier, ou un pour obtenir le second, qui est notre nouveau rapport.

Mais je n'aime pas vraiment faire cela beaucoup parce que c'est assez facile de se tromper dans le comptage des index. Et vous restez assis là à compter les lignes beaucoup et vous ajoutez une nouvelle sortie et soudainement tout casse. Donc

c'est beaucoup plus agréable de tout référencer par nom à la place. Et nous pouvons faire cela avec une clé spéciale ici appelée "emit".

Donc nous pouvons appeler cela comme nous voulons. Appelons cela emit outfile, et emit reports. Si vous définissez ceux-ci et vous pouvez le faire sur un ou plusieurs, c'est à vous de décider. Maintenant je peux descendre ici et à la place je peux aller dot out dot reports et juste l'appeler par nom, ce qui est beaucoup plus facile pour comprendre votre code quand vous le lisez, et c'est plus sûr face aux changements dans le code.

J'ai ajouté le .out.report ici, mais en fait j'ai besoin d'avoir deux sorties différentes qui sont publiées. Donc je vais renommer comme quelque chose de plus intéressant comme collected et report et est-ce que c'est ce que je l'ai appelé ? Je l'ai appelé out file, désolé. Donc ce nom emit ici outfile et report. parce que nous publions deux canaux de sortie différents et donc nous devons référencer les deux dans le bloc publish.

Ensuite nous devons aussi définir ceux-ci dans le bloc output. Donc j'ai renommé cela collected, et encore une fois, pour reports, un peu verbeux ici, mais c'est vraiment utile quand vous venez lire un nouveau workflow, de voir toutes les différentes sorties ici, tous les différents canaux listés côte à côte, et il y a des façons de rendre cela moins verbeux, que nous aborderons plus tard.

D'accord, essayons-le et exécutons notre workflow et voyons ce qui se passe.

Espérons que maintenant il devrait s'exécuter fondamentalement de la même façon qu'avant. Et nous allons obtenir un nouveau fichier de sortie ici appelé replicate_two, report. Et voilà. Il s'est ouvert et il dit qu'il y a trois messages de bienvenue dans le lot, ce qui est ce que nous attendions, donc c'est parfait.

Si je vais dans le répertoire work ici juste pour vous prouver qu'il a été exécuté dans le code Nextflow plutôt que dans le script bash, je peux aller à cat work/ command.sh, et vous verrez ici qu'il fait juste echo de cette chaîne directement. There were three greetings in this batch, et donc cette variable a été interpolée par Nextflow. Elle a été calculée dans le bloc script avant qu'il n'écrive le fichier .command.sh. Donc le calcul de variable résultant est fondamentalement codé en dur dans celui-ci avant qu'il ne soit exécuté sur votre environnement de calcul dans ce cas.

Et donc vous pouvez voir cette séparation entre le script. Bloquez ici et tout ce qui est au-dessus. J'espère que ça a du sens.

## À retenir et quiz

D'accord, c'est la fin de cette partie de Hello Nextflow. Donc comme avant, allez jeter un œil au quiz. Faites-le sur la page web ou dans la CLI, parcourez certaines des questions et vérifiez simplement que vous avez compris une partie du matériel que nous avons couvert. Voyez s'il y a quelque chose là qui met en évidence quelque chose que vous pourriez ne pas avoir compris. Pas trop de questions. Agréable et facile à faire. Ou vous pouvez le faire sur la page web ici aussi.

Et faites une petite pause, une petite promenade et revenez et rejoignez-nous dans la partie quatre de Hello Nextflow, où nous parlerons des modules. Merci beaucoup.
