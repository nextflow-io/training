# Partie 2 : Hello Channels - Transcription vidéo

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Traduction assistée par IA - [en savoir plus et suggérer des améliorations](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

<div class="video-wrapper">
  <iframe width="560" height="315" src="https://www.youtube.com/embed/yDR66fzAMOg?si=xCItHLiOQWqoqBB9&amp;list=PLPZ8WHdZGxmWKozQuzr27jyMGqp9kElVK&amp;cc_load_policy=1" title="YouTube video player" frameborder="0" allow="accelerometer; autoplay; clipboard-write; encrypted-media; gyroscope; picture-in-picture; web-share" referrerpolicy="strict-origin-when-cross-origin" allowfullscreen></iframe>
</div>

!!! note "Notes importantes"

    Cette page présente uniquement la transcription. Pour les instructions complètes étape par étape, retournez au [matériel de cours](../02_hello_channels.md).

    Les numéros de section indiqués dans la transcription sont fournis à titre indicatif uniquement et peuvent ne pas inclure tous les numéros de section du matériel.

## Bienvenue

Bonjour et bienvenue dans la Partie 2 de Hello Nextflow. Ce chapitre s'appelle Hello Channels.

Les canaux sont comme la colle de votre pipeline Nextflow. Ce sont les éléments qui maintiennent ensemble tous les différents processus, que Nextflow utilise pour faire circuler toutes les informations et orchestrer votre workflow.

Il y a un autre aspect des canaux qui concerne les opérateurs. Ce sont essentiellement des fonctions que nous pouvons utiliser sur les canaux pour modifier leur contenu. Plongeons dans VS Code et voyons où nous en sommes.

J'ai beaucoup zoomé sur VS Code, donc pour garder les choses propres et organisées, j'ai supprimé tous les fichiers _.nextflow\*_ et le répertoire _work/_, ainsi que _results/_ et tout ce qui concerne le Chapitre Un. Je repars simplement à zéro ici. Mais ne vous inquiétez pas trop pour ça. Si vous ne voulez pas, vous pouvez laisser ces fichiers. Ils ne causeront aucun problème.

Nous allons commencer par travailler sur _hello-channels.nf_ pour ce chapitre, et si j'ouvre ce fichier, il devrait ressembler beaucoup au fichier sur lequel nous travaillions précédemment. Il se peut que certaines parties soient à des endroits différents du script, mais tout devrait être essentiellement identique.

Une chose qui est différente, c'est que le chemin dans le bloc output ici est maintenant _hello_channels_ pour cette partie, ce qui signifie que les fichiers de résultats seront stockés dans un sous-répertoire différent dans vos résultats si vous l'avez encore. Ce devrait donc être un endroit propre et bien rangé pour commencer sans être confus au sujet des sorties.

Bon, rappelons-nous rapidement ce que fait ce script quand nous exécutons ce workflow. Nous faisons _"nextflow run hello-channels.nf"_. Nous pouvons faire _"--input myinput"_, et quand nous exécutons ceci, il va utiliser ce paramètre, params.input, qui a été passé comme variable pour le processus sayHello ici, qui va dans greeting et est sauvegardé dans output.txt. Et nous pouvons voir ça dans le fichier de résultats. Parfait.

## 1. Fournir des entrées variables via un canal explicitement

C'est bien. Mais c'est assez simpliste. Nous avons une variable dans ce paramètre, qui entre dans un processus qui s'exécute une fois, et cela ne fait pas vraiment passer à l'échelle. Et nous ne pouvons pas lui donner beaucoup de fichiers différents à créer ici. Nous ne pouvons pas lui donner beaucoup de salutations différentes. Nous n'en avons qu'une.

En réalité, Nextflow concerne la mise à l'échelle de votre analyse. Vous voulez donc probablement qu'il fasse plus d'une chose. Et nous faisons cela avec des _canaux_.

Les canaux sont un concept un peu unique pour de nombreuses personnes qui découvrent Nextflow. Cela vient de ce genre de concepts de programmation fonctionnelle, et il peut falloir un peu de temps pour bien comprendre, mais une fois que ça fait tilt, ils débloquent vraiment la puissance de Nextflow et c'est essentiel à la façon dont vous écrivez vos workflows.

## 1.1. Créer un canal d'entrée

Commençons par prendre ce script et faire en sorte qu'il utilise un _canal_ au lieu d'un simple _param_.

Nous descendons au workflow, qui est l'endroit où toute notre logique de workflow consiste à assembler les choses. Et je vais aller ici et je vais créer un nouveau canal.

Créer un nouveau canal.

Et je vais l'appeler "_greeting_ch"_. C'est une convention de faire "_\_ch"_ comme ceci, juste pour que vous puissiez vous rappeler que cette variable est un canal. Mais vous pouvez l'appeler comme vous voulez.

Et ensuite je vais dire égal, et je vais faire _"channel.of"._

Channel est comme l'espace de noms pour tout ce qui concerne les canaux. Avec un "c" minuscule si vous avez déjà utilisé Nextflow. Et le _".of"_ est quelque chose qu'on appelle une fabrique de canaux (Channel factory), qui est essentiellement un moyen de créer un canal.

Il existe de nombreuses fabriques de canaux différentes. Si je fais juste "." ici, vous pouvez voir que VS Code en suggère plein, mais _".of"_ est la plus simple et prend simplement une entrée ici.

Donc je peux mettre des parenthèses et je vais dire _"Hello Channels!"_.

Parfait. J'ai un canal. Fantastique. Je peux sauvegarder, je pourrais l'exécuter à nouveau, mais rien d'intéressant ne va se passer. VS Code m'a mis une ligne d'avertissement orange ici et m'a dit que c'est configuré : vous avez créé ceci, mais vous ne l'avez jamais vraiment utilisé pour quoi que ce soit. Ce canal n'est pas consommé.

D'accord, alors comment l'utilisons-nous ? Très simple. Je vais prendre ceci, le copier, et je vais supprimer _params.input_ et je vais mettre _"greeting_ch"_ ici à la place. Nous allons donc passer ce canal comme entrée à sayHello.

Notez que j'ai codé en dur cette chaîne pour l'instant. C'est un peu un retour en arrière après notre joli param que nous avons utilisé à la fin du dernier chapitre, mais ça simplifie juste les choses pour commencer afin que vous puissiez voir la logique.

D'accord, je vais aller dans mon terminal et je vais exécuter le workflow à nouveau. Sans aucun _"--input"_ cette fois, et ça va s'exécuter et ça va utiliser ce canal que nous avons créé et j'espère que nous devrions avoir un fichier ici dans _results/hello_channels/_ et il dit maintenant "Hello Channels!". Fantastique. C'est donc ce que nous espérions de notre canal ici. Parfait.

## 1.4. Utiliser view() pour inspecter le contenu du canal

Une autre chose à ajouter ici, juste une introduction rapide à une autre fonction que nous pouvons utiliser sur les canaux appelée "_.view"_.

C'est analogue à la commande _print_ en Python ou dans d'autres langages auxquels vous pourriez être habitué·e, et elle affiche simplement le contenu de ce canal dans le terminal quand nous l'exécutons.

Donc faire "_.view"_, et ensuite si je réexécute le workflow à nouveau, il devrait afficher dans le terminal quel est le contenu de ce canal, au moment où nous l'avons créé.

Effectivement, vous pouvez voir qu'il est affiché dans le terminal ici. _"Hello Channels!"_.

Notez que vous pouvez diviser ces choses sur plusieurs lignes si vous le souhaitez, et en fait, le formateur automatique de Nextflow essaiera de le faire pour vous. Les espaces blancs ne sont pas vraiment importants ici, donc vous pouvez enchaîner ces choses les unes après les autres.

## 2. Modifier le workflow pour s'exécuter sur plusieurs valeurs d'entrée

D'accord, donc notre canal a une chose dedans, ce qui est bien, mais c'est essentiellement la même chose qu'avant. Alors rendons-le un peu plus compliqué. Ajoutons quelques éléments supplémentaires dans notre canal.

La fabrique de canaux "_.of()"_ peut prendre plusieurs éléments, alors écrivons-en quelques autres. Nous allons faire _Hello, Bonjour, Hej_. Et ensuite nous pouvons exécuter ce workflow à nouveau et nous verrons ce qui se passe.

Devrait s'exécuter à nouveau. Et nous avons maintenant affiché _"Hello", "Bonjour"_ et _"Hej"_ dans le terminal avec notre instruction view. Fantastique.

## 2.1.2. Exécuter la commande et examiner la sortie du journal

Vous pourriez penser que nous avons terminé à ce stade. Mais en fait, il y a un piège ici qui va nous faire trébucher. Si nous regardons notre fichier de sortie ici, vous pouvez voir qu'il contient _"Hello"_, mais il n'a aucune des autres sorties. En fait, c'est juste celui-ci.

Si nous exécutons ce workflow plusieurs fois, nous pourrions même voir que parfois il a _"Bonjour"_, parfois il a _"Hej"_. C'est un peu aléatoire.

Si nous regardons le terminal, nous pouvons voir qu'il s'est exécuté trois fois et nous pouvons voir les différentes sorties view. Mais si je vais dans le répertoire work, je peux faire _"cat work"_. Mettre ce hash et développer ça et _output.txt_. Vous pouvez voir que ce fichier dans le répertoire work est différent du répertoire results, et celui-ci est _"Hej"_. Il y a donc quelque chose qui ne fonctionne pas tout à fait correctement ici.

Et le point clé est que nous avons trois tâches qui se sont exécutées. La sortie de Nextflow essaie de résumer cela au fur et à mesure du traitement, pour qu'elle ne prenne pas complètement le contrôle de tout votre terminal, et ce logging ANSI utilise des codes d'échappement ANSI, a essentiellement écrasé les autres tâches. Il vous montre donc juste la dernière qui s'est trouvée être mise à jour.

## 2.1.3. Exécuter la commande à nouveau avec l'option -ansi-log false

Il y a quelques choses que nous pouvons faire pour vraiment mieux comprendre cela. Nous pouvons regarder dans le répertoire work lui-même et vous pouvez voir tous les différents répertoires work là, mais c'est un peu confus car ils seront mélangés avec différentes exécutions Nextflow.

Ou nous pouvons dire à Nextflow de ne pas utiliser les codes d'échappement ANSI.

Donc si j'exécute la commande à nouveau, mais cette fois je dis _"-ansi-log false"_ pour le désactiver, je pourrais aussi utiliser les variables d'environnement _$NO_COLOR_ ou _"$NXF_ANSI_LOG=false"_. Ensuite, il utilise le style de journalisation Nextflow plus à l'ancienne sans aucun de ces codes d'échappement. Il affiche simplement directement dans un terminal sans mises à jour astucieuses.

Et maintenant nous pouvons voir tous ces trois processus qui se sont exécutés. Et chacun d'eux a son propre hash de tâche. Et si nous allons dans ces répertoires work, nous verrons les trois différentes salutations que nous avons spécifiées.

Donc ça a un peu plus de sens maintenant. J'espère que vous comprenez que Nextflow faisait cela, il était juste un peu malin avec ce qu'il vous montrait dans le terminal avec ces répertoires work.

Cependant, cela a résolu un problème avec les répertoires work, mais ça n'a pas résolu le problème avec le fichier de sortie. Nous avons toujours juste un fichier de sortie qui dit _"Hello"_.

## 2.2. S'assurer que les noms de fichiers de sortie seront uniques

Maintenant pour comprendre ceci, nous devons revenir à notre script de workflow. Nous générons notre canal ici, nous le passons à notre processus, et si nous regardons le processus, nous écrivons la salutation dans un fichier appelé _"output.txt"_ et nous repassons ce fichier de sortie au bloc output ici, en le publiant.

Cependant, chaque fois que ce processus s'exécute, ces trois tâches différentes génèrent toutes un fichier appelé _"output.txt"_, tous ces fichiers de sortie sont publiés dans le répertoire results, et ils s'écrasent tous les uns les autres. Donc quel que soit le fichier de résultat que vous obtenez là, c'est juste le dernier qui a été généré, mais qui a écrasé tous les autres. Ce n'est pas vraiment ce que nous voulons.

## 2.2.1. Construire un nom de fichier de sortie dynamique

Il existe différentes façons de gérer cela, mais la plus simple pour l'instant est simplement de créer différents noms de fichiers uniques. Donc chaque fois que la tâche s'exécute avec une salutation différente, elle générera un fichier de sortie différent, qui ne sera plus en conflit lors de la publication. Et ensuite nous obtiendrons trois fichiers de sortie uniques.

Nous faisons cela exactement de la même manière. Nous pouvons utiliser cette variable n'importe où dans le bloc script et nous pouvons l'utiliser plusieurs fois.

Donc je peux la coller ici, _"$\{greeting\}\_output.txt"_, et ensuite je dois aussi la coller ici parce que nous ne créons plus un fichier appelé _output.txt_. Donc si je ne mets pas à jour ceci, Nextflow plantera avec une erreur disant qu'il attendait un fichier qui n'a jamais été généré.

Donc je dois faire la même chose là et je dois utiliser des guillemets doubles, pas des guillemets simples, pour que cette variable soit comprise.

D'accord, essayons et voyons si ça a marché. Nous allons exécuter le workflow à nouveau. J'espère qu'il nous montrera les trois tâches différentes dans les trois répertoires work différents. Et effectivement, vous pouvez voir dans le dossier results ici à gauche, nous avons maintenant trois fichiers différents avec trois noms de fichiers différents et chacun avec le contenu différent que nous attendons. Donc les fichiers ne s'écrasent plus les uns les autres, et tout est là comme nous l'attendons.

C'est une configuration un peu triviale que nous avons parcourue ici, mais elle souligne certains des concepts clés que vous devez comprendre sur le fonctionnement de la publication de fichiers, et certaines des choses dans lesquelles vous pourriez tomber comme des pièges. J'espère donc que vous pourrez éviter cela dans vos propres workflows.

Il vaut également la peine de noter que ce que nous avons fait ici est un peu impraticable dans des situations de la vie réelle. Nous avons pris des données d'entrée et nous utilisons ces données, mais nous nommons aussi le fichier d'après ces données, ce que vous ne pouvez généralement pas faire.

Donc dans de vrais pipelines Nextflow plus matures, vous passerez souvent un objet meta avec toutes les métadonnées associées à un échantillon donné. Vous pouvez ensuite créer des noms de fichiers dynamiques basés sur cela, ce qui est beaucoup plus pratique.

Si vous êtes intéressé·e par comment faire cela avec les meilleures pratiques, il y a une quête secondaire sur _training.nextflow.io_, qui concerne spécifiquement les métadonnées et les meta maps, donc vous pouvez creuser là pour plus de détails.

## 3. Fournir plusieurs entrées via un tableau

D'accord. Ensuite, nous allons explorer un peu comment les canaux sont structurés et comment ils diffèrent d'autres types de structures de données dans le langage de codage. Et je vais réfléchir un peu à comment je pourrais potentiellement utiliser un tableau, qui pourrait être un concept familier si vous venez d'autres langages.

Puis-je utiliser un tableau dans un canal ? Essayons. Je vais créer un tableau, et j'ai copié ceci des docs, _"greetings_array"_ et _"Hello", "Bonjour"_ et _"Holà"_. Et ensuite je vais mettre ça ici au lieu de mes chaînes codées en dur. Donc je vais dire "channel.of" _"greetings_array"_, en passant ce tableau dans un canal. Essayons.

Ouvrir le terminal, et exécuter le pipeline.

D'accord. Vous pouvez voir que l'instruction view ici a bien affiché notre tableau comme prévu, mais ensuite tout ce texte rouge, ou il ne sera pas rouge si vous avez toujours _"-ansi-log"_ désactivé, mais tout ce texte rouge nous dit que quelque chose s'est mal passé.

Nous n'avons plus une belle coche verte ici. Nous avons une croix rouge, et si je rends juste cela un peu plus large pour que ce soit plus facile à lire, Nextflow nous dit ce qui s'est mal passé.

Donc décomposons ceci section par section. Il dit que l'erreur a été causée par, et ensuite la raison de l'erreur, qui est des fichiers de sortie manquants. Donc essentiellement ce bloc output disait que ce fichier devrait être créé et il ne l'a pas été. Ensuite, il dit que c'est la commande qui a été exécutée. Donc c'est essentiellement le contenu de ce fichier _.command.sh_. C'est à quoi il ressemblait après que toutes ces variables aient été insérées.

Et vous pouvez voir ici que notre commande echo n'a en fait été exécutée qu'une seule fois et elle a utilisé le tableau entier, mais dans une représentation en chaîne, ce qui n'est pas vraiment ce que nous voulions.

Et ensuite la commande s'est terminée comme ça, et c'était le répertoire work où nous pouvons aller voir les fichiers pour mieux comprendre.

D'accord. Donc ce qui s'est passé c'est que Nextflow a juste passé ce tableau entier comme un seul élément de canal au processus, ce qui signifie que le processus ne s'est exécuté qu'une seule fois. Il avait une tâche et il n'a pas utilisé les données dans une structure que nous attendions.

## 3.2. Utiliser un opérateur pour transformer le contenu du canal

Donc nous devons faire quelque chose à ce canal d'abord, avant qu'il puisse être utilisé. Et cela prépare le terrain pour l'utilisation d'opérateurs, qui sont des fonctions spéciales que nous pouvons utiliser sur les canaux pour manipuler le contenu des canaux.

Dans ce cas, nous allons utiliser quelque chose appelé _flatten_. Que nous passons à la fin du canal ici. Donc nous créons le canal et ensuite nous exécutons _flatten_. Et encore, si nous passons la souris dessus, il nous montre la documentation pour cette commande directement dans VS Code, ce qui est très utile. Vous pouvez aussi trouver toutes ces docs sur le site web de Nextflow, dans la documentation.

Je pourrais juste exécuter ce code maintenant et voir si ça marche, mais c'est aussi une bonne occasion d'introduire comment faire du code dynamique dans les opérateurs et dans le code Nextflow, qui s'appellent des closures.

Donc je vais rajouter une commande view ici avant d'exécuter _flatten_. Et ici celle-ci a ces accolades, qui sont la closure dynamique. Et il y a juste du code arbitraire à l'intérieur ici qui sera exécuté, dans le contexte d'un opérateur view.

Ici, cela dit prendre la salutation, qui est l'entrée de l'opérateur view, et c'est ici. Je pourrais appeler ceci comme je veux, je pourrais appeler ceci _"foo"_ et je dois juste y faire référence comme _"foo"_ plus tard. Et ensuite je dis avec ceci, retourner ceci.

Et ensuite définir le retour d'une chaîne qui dit avant le flatten pour une variable. Très simple.

Je vais maintenant en ajouter un autre exactement pareil, mais je vais dire après _flatten_.

Donc ce que cela fait, parce que cela s'exécute en séquence, vous allez voir à quoi ressemble le canal avant d'exécuter _flatten_, et ensuite à nouveau après avoir exécuté _flatten_.

Et ensuite ce canal greeting est toujours créé, donc il va toujours être passé au processus. Et j'espère que maintenant le workflow va s'exécuter. Essayons.

Parfait. Donc tout d'abord, le pipeline n'a pas planté cette fois. Nous avons eu trois processus qui se sont exécutés correctement et nous avons une petite coche. Et ensuite nous pouvons voir que nos instructions view ont fonctionné.

Nous avons avant _flatten_, qui est ce tableau que nous avons vu avant de l'échec, et ensuite nous avons trois fois le après _flatten_ a été appelé où nous avons _"Hello", "Bonjour"_, et tous ces trois éléments séparés dans le tableau, qui sont maintenant comme nous l'espérions, trois éléments séparés dans le canal.

Et vous pouvez voir que l'opérateur _view_ a été exécuté trois fois. Et c'est parce que ce canal après _flatten_ a maintenant trois éléments. Et donc l'opérateur est appelé trois fois.

Très rapidement, je mentionnerais juste que quand je créais des fabriques de canaux avant, j'ai fait _"."_, et ensuite nous avons vu qu'il y avait beaucoup de façons différentes de créer des canaux, et l'une d'elles s'appelle "_fromList"_. Et c'est en fait spécifiquement conçu pour faire cette même opération. Donc nous aurions pu juste faire fromList greetings away, et ça marchera. C'est une syntaxe légèrement plus propre et plus agréable. Mais pour les besoins de cette démonstration, nous voulions le faire étape par étape pour que vous puissiez voir comment le canal est manipulé et comment différents opérateurs peuvent changer ce qui est dans le contenu d'un canal.

## 4. Lire les valeurs d'entrée depuis un fichier CSV

D'accord, comment pouvons-nous rendre cela un peu plus réaliste ? Vous n'allez probablement pas vouloir créer beaucoup de code dans votre pipeline Nextflow avec des tableaux codés en dur. Vous allez probablement vouloir prendre les données de l'extérieur quand vous lancez, et ces données seront presque certainement dans des fichiers.

Donc la prochaine chose que nous allons faire, c'est que nous allons répliquer cela, mais au lieu de prendre les données d'un seul paramètre CLI ou d'une chaîne ou d'un tableau codé en dur, nous allons les prendre d'un fichier.

Alors débarrassons-nous de notre greetings away. Et maintenant nous allons changer cette fabrique de canaux à nouveau. Je viens de dire qu'il y en avait plein parmi lesquelles choisir et il y en a une appelée _".fromPath"_. Et je vais lui dire de, dans ce cas, prendre _params.input_, ce qui revient à notre entrée que nous utilisions plus tôt.

Maintenant ce paramètre n'est pas vraiment prêt à être utilisé encore. Nous disons toujours que c'est une chaîne et elle est codée en dur ici avec une valeur par défaut, mais nous pourrions écraser cette chaîne. Nous voulons maintenant que ce soit un fichier à la place. Donc le type est différent. Ce n'est plus une _String_. C'est un _Path_.

Et ensuite nous pouvons définir la valeur par défaut si nous voulons, à nouveau, vers un Path. Et si je regarde dans explorer à gauche, vous pouvez voir dans ce dépôt, dans ce répertoire de travail, j'ai un répertoire appelé data. J'ai un fichier dedans appelé _"greetings.csv"._

Donc je peux juste définir la valeur par défaut ici à _"data/greetings.csv"_. Maintenant, quand j'exécute ce pipeline à nouveau sans aucune option de ligne de commande, il utilisera cette valeur par défaut. Il sait que c'est un path, donc il sait qu'il devrait gérer cela comme un path et non une chaîne.

Et ensuite il va passer cela dans une fabrique de canaux à partir de ce _params.input_ et créer notre canal, qui va ensuite être utilisé dans ce processus appelé _sayHello_. Essayons.

D'accord. Échec. Ne vous inquiétez pas. C'était attendu. Et si vous suivez le matériel de formation, vous verrez que c'était attendu là aussi. Voyons ce qui se passe ici.

Il a essayé d'exécuter le pipeline. Il a essayé d'exécuter le processus, et il a eu une erreur assez similaire à celle que nous avons vue avant.

Ici il dit : nous avons essayé d'exécuter _echo_, mais au lieu d'afficher le contenu de ce fichier CSV, il a juste affiché le chemin. Et vous pouvez voir que c'est le chemin absolu complet ici vers ce fichier CSV.

Et ensuite bien sûr, parce qu'il a essayé d'écrire cela dans ce chemin vraiment compliqué, il ne savait pas vraiment quoi faire. Et c'était en dehors de la portée du répertoire work du processus.

J'ai mentionné au début que Nextflow encapsule chaque tâche exécutée dans un répertoire work spécial. Et si vous essayez d'écrire vers des données qui sont en dehors de ce répertoire work, Nextflow vous arrêtera par mesure de sécurité. Et c'est ce qui s'est passé ici. Nous avons essayé d'écrire vers un chemin absolu et Nextflow a échoué et nous a empêchés.

## 4.2. Utiliser l'opérateur splitCsv() pour analyser le fichier

D'accord, jetons un œil à ce canal et voyons à quoi il ressemble. Nous pouvons faire _".view"_, et j'ai copié ceci du site web. Donc _.view_, et nous avons une closure dynamique ici et nous disons un nom de variable "_csv"_ comme entrée. Donc c'est le contenu du canal, et nous disons avant splitCsv, et c'est à quoi il ressemble.

Si je l'exécute à nouveau, il échouera toujours, mais il nous montrera ce qu'il y a dans ce canal. Ce n'est pas particulièrement excitant. C'est cette variable _path_. Donc vous pouvez voir que c'est juste une chaîne ici parce qu'elle est affichée dans un terminal, mais c'est un objet _path_, qui contient les informations et les métadonnées sur ce fichier.

Nous ne voulons pas passer les métadonnées du fichier à l'entrée. Nous voulons passer le contenu de ce fichier. Si nous regardons le fichier _greetings.csv_, vous pouvez voir ici qu'il a ces différentes variables ici. _Hello, Bonjour, Holà_ à nouveau. Et ce sont vraiment ces choses que nous voulons passer à notre processus, pas juste le fichier lui-même comme un seul objet.

Donc nous devons analyser ce fichier CSV. Nous devons le déballer, accéder au contenu du fichier CSV, et ensuite passer le contenu dans le canal au processus.

Comme vous pouvez probablement le dire d'après le message de log, nous voulons utiliser le _splitCsv_, qui est un autre opérateur, un autre opérateur de canal. Donc si je fais "_dot" "s"_, et ensuite vous pouvez voir qu'il est auto-suggéré. Oups, _splitCsv_ et des parenthèses.

Et ensuite après _splitCsv_, je vais mettre une autre instruction _view_ juste pour que nous puissions voir à quoi ça ressemble après. Exécutons le pipeline et voyons ce que nous avons.

D'accord. Il a toujours échoué, mais d'une nouvelle façon passionnante, ce qui est un progrès.

Cette fois encore, nous avons un problème avec notre script, qui a été rendu. Maintenant. Nous n'avons plus le chemin final, mais nous avons un tableau de variables, qui ressemble beaucoup à l'erreur que nous avons eue plus tôt quand nous passions un tableau comme entrée fixe.

Avec notre journalisation de l'opérateur view, nous pouvons voir qu'avant _splitCsv_ c'était le chemin. Et effectivement, après _splitCsv_, nous avons trois sorties différentes et chacune de ces sorties ressemble énormément à chacune des lignes du fichier _greetings.csv_, ce qui est logique.

Donc ce qui s'est passé ici, c'est que Nextflow a analysé ce fichier CSV, nous a donné trois objets, un tableau pour chaque ligne du fichier CSV. Donc ensuite trois fois nous avons passé un tableau de variables au canal au lieu d'une seule valeur de chaîne.

D'accord, donc la dernière fois que nous avons eu ce problème, nous avons utilisé _flatten_. Essayons juste très rapidement. Essayons flatten et voyons ce qui se passe.

Je peux appeler ces variables comme je veux. Donc je vais l'appeler _myarray_ parce que ce n'est plus vraiment un CSV. Essayons d'exécuter à nouveau et voyons ce qui se passe avec _flatten_.

Donc cette fois nous allons exécuter, nous avons analysé le CSV en trois objets tableau, et ensuite nous l'avons aplati. Et cette fois ça a réussi. Et le pipeline Nextflow s'est exécuté. Cependant, vous pouvez voir que _flatten_ y va vraiment à fond et aplatit tout. Et donc nous obtenons trois entrées de tableau indépendantes pour chaque ligne. Et donc il a exécuté le processus trois fois pour chaque ligne du CSV. Et maintenant nous avons tout un tas de fichiers de résultats, et 123, 456, et toutes sortes de choses, pas seulement cette première colonne du CSV, qui est ce que nous voulions vraiment.

## 4.3. Utiliser l'opérateur map() pour extraire les salutations

Donc comment accédons-nous juste à la première colonne ? Si flatten est trop simpliste ici, nous avons besoin d'un opérateur plus complexe où nous pouvons réellement personnaliser et lui dire ce que nous voulons du CSV.

Pour faire cela, nous allons utiliser _map_. Essentiellement _map_ dit juste, exécuter du code, une fonction sur chaque élément que je reçois et faire une sorte de transformation dessus. Et parce que c'est si flexible, vous le verrez apparaître dans le code Nextflow tout le temps.

En soi, il ne fait rien. Donc nous ne voulons pas de parenthèses ordinaires, nous voulons une closure ici et nous devons lui dire quoi faire. Donc je vais dire _"row"_, parce qu'on reçoit des lignes du CSV, donc c'est un nom de variable logique. C'est l'entrée. Et je veux retourner juste le premier élément de ce tableau.

Les tableaux dans Nextflow sont basés sur zéro, donc nous allons dire juste le premier élément, qui est row zéro. Si nous voulions la deuxième colonne, ce serait un ou la troisième colonne serait deux, et ainsi de suite. Nous pouvons retourner ce que nous voulons ici, mais je vais retourner juste la première valeur.

Et maintenant, nous pouvons exécuter le pipeline à nouveau et voir s'il fait ce que nous attendons.

Effectivement, après _splitCsv_ nous avons nos tableaux, et ensuite après le _map,_ nous avons nos belles chaînes propres, juste _"Hello", "Bonjour"_ et _"Holà"_. Et le pipeline fait maintenant ce que nous voulons qu'il fasse. Fantastique.

Donc nous pouvons nous débarrasser de toutes ces commandes view maintenant. Nous n'en avons plus besoin.

## Récapitulatif

Nous avons terminé notre genre de débogage et c'est le code avec lequel nous nous retrouvons. Prenant notre paramètre CLI appelé _input_, qui est classé comme un _Path_. Nextflow trouve le chemin, le charge, et comprend le fichier CSV. Retourne toutes les différentes lignes. Et ensuite nous mappons juste le premier élément de cette ligne dans le canal qui nous donne en quelque sorte le contenu du canal, qui est passé au processus.

Et le processus s'exécute sur chaque élément du canal, qui sont trois. Et il exécute le processus trois fois, lui donnant trois tâches. Et ces résultats sont ensuite publiés depuis le workflow, récupérés par la sortie du processus. Publiés depuis un workflow et sauvegardés dans le bloc output vers un sous-répertoire appelé _"hello_channels"_.

Plutôt cool. Nous arrivons maintenant à quelque chose qui ressemble plus étroitement à un pipeline Nextflow de la vie réelle que vous pourriez exécuter pour une vraie analyse.

## À retenir

D'accord. J'espère que vous avez maintenant une idée de ce que sont les canaux et opérateurs Nextflow et comment les opérateurs fonctionnent sur les canaux et comment vous pouvez les créer.

Les canaux, comme je l'ai dit au début de cette vidéo, sont la colle de Nextflow. Et vous pouvez voir ici que nous pouvons prendre différentes entrées et les manipuler et prendre ces données et ensuite les passer dans la logique de workflow en aval.

Et ce bloc workflow ici est vraiment l'endroit où vous construisez toute cette parallélisation et toute la logique intelligente, et expliquez à Nextflow comment construire votre DAG de workflow, et comment orchestrer votre pipeline.

Les canaux ne sont pas le concept le plus facile à comprendre. Donc prenez une pause, réfléchissez un peu à cela, peut-être relisez le matériel, et assurez-vous vraiment que vous avez bien compris ces concepts parce que c'est la clé de votre compréhension de Nextflow et plus vous comprenez les canaux et les différents opérateurs de canaux et les différentes fabriques de canaux, plus vous vous amuserez à écrire Nextflow et plus vos pipelines seront puissants.

Ce n'est pas la même chose que la programmation régulière en Python ou dans d'autres langages. Nous n'utilisons pas d'instructions _if_ ici, c'est de la programmation de flux fonctionnelle utilisant des canaux et des opérateurs. Donc c'est un peu différent, mais c'est aussi super puissant.

C'est la fin de ce chapitre. Allez faire une petite pause et je vous verrai dans la prochaine vidéo pour la partie trois où nous allons parcourir Hello Workflow, et parler un peu plus des workflows.

Tout comme le chapitre précédent, il y a quelques questions de quiz en bas de la page web ici, donc vous pouvez les parcourir rapidement et vous assurer que vous comprenez toutes les différentes parties du matériel que nous venons de faire. Et à part cela, je vous verrai dans la prochaine vidéo. Merci beaucoup.

D'accord.
