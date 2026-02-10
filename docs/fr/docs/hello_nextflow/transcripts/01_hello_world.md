# Partie 1 : Hello World - Transcription vidéo

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Traduction assistée par IA - [en savoir plus et suggérer des améliorations](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

<div class="video-wrapper">
  <iframe width="560" height="315" src="https://www.youtube.com/embed/tOukLxWCHiA?si=F0t9LFYLjAWoyRXj&amp;list=PLPZ8WHdZGxmWKozQuzr27jyMGqp9kElVK&amp;cc_load_policy=1" title="YouTube video player" frameborder="0" allow="accelerometer; autoplay; clipboard-write; encrypted-media; gyroscope; picture-in-picture; web-share" referrerpolicy="strict-origin-when-cross-origin" allowfullscreen></iframe>
</div>

!!!note "Notes importantes"

    Cette page présente uniquement la transcription. Pour les instructions complètes étape par étape, retournez au [matériel de formation](../01_hello_world.md).

    Les numéros de section indiqués dans la transcription sont fournis à titre indicatif uniquement et peuvent ne pas inclure tous les numéros de section présents dans le matériel.

## Bienvenue

Bonjour et bienvenue.

Vous êtes maintenant dans la Partie 1 de la formation « Hello Nextflow » intitulée « Hello World ». Dans ce chapitre, nous allons commencer à développer une compréhension des bases de Nextflow.

Vous devriez maintenant être configuré·e dans Codespaces ou un environnement équivalent avec VS Code en cours d'exécution, et vous devriez avoir votre dossier Hello Nextflow dans l'espace de travail dans l'Explorateur avec tous ces différents fichiers ici.

Nous allons commencer par faire des choses très basiques dans le terminal en utilisant Bash, puis nous verrons si nous pouvons faire les mêmes choses dans Nextflow pour que vous ayez une idée de ce à quoi ressemble la syntaxe.

## 0. Échauffement

Commençons vraiment simplement. Commençons juste avec « echo », pour afficher quelque chose dans un terminal. « Hello World ». J'appuie sur entrée et cela s'affiche dans un terminal. Hello World. J'espère que ce n'est une surprise pour personne qui regarde cette formation.

D'accord, faisons quelque chose avec ça. Plutôt que de simplement l'afficher dans le terminal, écrivons-le dans un fichier. Je vais appuyer sur la flèche haut de mon clavier, ce qui parcourt l'historique Bash, donc cela me donne ma dernière commande, et je vais ajouter à la fin un petit symbole supérieur à, qui redirige la sortie de cette commande vers un fichier, et je vais l'appeler output.txt.

J'appuie à nouveau sur entrée pour exécuter cette commande, rien dans le terminal cette fois, mais nous pouvons voir sur le côté gauche qu'un nouveau fichier est apparu ici, appelé output.txt.

Nous pouvons visualiser cela dans un terminal avec quelque chose comme cat. Donc cat output.txt et effectivement il affiche « Hello World ». Nous pouvons également double-cliquer dessus et cela l'ouvre dans l'éditeur de code dans VS Code.

## 1.1. Examiner le code

Très bien. Je vous avais dit que c'était simple. Et ensuite ? Essayons de reprendre ce processus et de le refaire, mais cette fois, faisons-le dans Nextflow.

Comme je l'ai dit, tous les différents chapitres de cette formation commencent par un script et celui-ci s'appelle Hello World. Je vais donc trouver Hello World. Il s'affiche en aperçu si je clique une fois dessus, je vais double-cliquer pour l'ouvrir dans l'éditeur ici. Et je vais rapidement me débarrasser du terminal.

Maintenant, c'est un script très simple, aussi simple que possible. Il ne fait que 22 lignes, et il fait essentiellement la même chose. En fait, une partie de cela devrait vous sembler familière. C'est ce que nous venons de taper. Nous pouvons voir notre commande bash redirigeant vers un fichier là.

D'accord. Quoi d'autre ? Également, dans ce fichier, nous pouvons commencer à voir certains des concepts fondamentaux de Nextflow. Nous avons un **process** en rouge ici et un **workflow**. Ce sont des mots-clés spéciaux et une terminologie spéciale dans Nextflow.

## 1.1.1. La définition du process

Différents **process** au sein d'un workflow encapsulent différentes unités logiques de votre workflow. Chaque **process** fait une chose.

Lorsque nous l'exécutons, il génère une tâche ou plusieurs tâches, qui sont des étapes d'exécution réelles d'un pipeline. Tous les **process** sont ensuite orchestrés dans un bloc **workflow**, que nous voyons en bas, et dans ce cas, exécute simplement ce seul **process**.

Le nom du **process** suit ce mot-clé ici, et cela peut être essentiellement n'importe quoi. Et ensuite, le contenu du **process** se trouve entre ces accolades.

Il n'y a vraiment qu'une seule exigence pour un **process**, qui est qu'il inclue une sorte de bloc script ou exec. C'est dans les triples guillemets ici, et c'est le script bash qui est écrit dans le répertoire de travail lorsque nous exécutons le pipeline et c'est la chose qui s'exécute réellement sur votre ordinateur ou serveur.

C'est généralement du bash, mais vous pouvez également mettre un shebang différent ici en haut, et cela pourrait être un script Python ou un script R. Peu importe. Tout ce qui se trouve dans ce script sera exécuté.

Il y a une autre chose que nous avons ajoutée dans ce **process** ici, qui est la déclaration output. Cela indique à Nextflow que ce **process** attend un fichier de sortie appelé output.txt. Il indique que c'est un path, donc il doit être traité comme un fichier, pas disons, si c'était val, cela indiquerait que c'est comme une variable ou une valeur.

Notez que cela ne crée pas ce fichier. Cela ne le génère pas réellement. C'est fait par le script ici en bas. Cela indique simplement à Nextflow de s'attendre à un fichier de sortie avec ce nom de fichier.

## 1.1.2. La définition du workflow

D'accord. Et ensuite en bas nous avons un **workflow** ici, et encore une fois, nous avons une déclaration. Celui-ci s'appelle Main. C'est l'équivalent workflow d'un bloc script, si vous voulez. C'est la partie du workflow qui fait quelque chose. Et dans ce cas, nous disons, appeler le **process** appelé sayHello.

Normalement, bien sûr, votre pipeline aura l'air beaucoup plus complexe que cela. Vous aurez probablement plus d'un **process**, et vous utiliserez des canaux pour orchestrer le flux de données entre eux. Nous allons aborder cela dans les prochaines parties de cette formation, mais pour l'instant, c'est suffisant. C'est un pipeline valide, qui devrait fonctionner.

Je peux même cliquer sur preview DAG ici dans VS Code. Le DAG est une représentation d'une structure de flux de données dans le pipeline, et nous pouvons le voir rendu sur le côté sous forme de diagramme mermaid. Dans ce cas, c'est très simple. Il y a une boîte, qui est le workflow et un **process**, qui s'appelle sayHello, mais cela pourrait devenir plus intéressant au fur et à mesure.

## 1.2. Exécuter le workflow

D'accord, essayons d'exécuter ce workflow et voyons ce qui se passe.

Je vais faire réapparaître le terminal en bas, effacer la sortie, et je vais taper Nextflow Run. Et ensuite je vais juste taper le nom du script, qui est hello-world.nf. Et je vais appuyer sur entrée.

D'accord, il y a des informations standard en haut, qui nous indiquent que Nextflow s'est exécuté et quelle version était en cours d'exécution et quel était le nom du script et tout.

Et vraiment la chose importante que nous recherchons ici est _ici_, qui est un résumé des différentes tâches qui ont été exécutées.

Si le vôtre ressemble à ceci avec une petite coche verte, alors bravo. Vous venez d'exécuter votre premier pipeline. Fantastique.

Il nous indique ici le nom du **process**, qui s'est exécuté, qui s'appelait Say Hello, et il nous a dit qu'il s'est exécuté une fois et qu'il a réussi. Cela se met à jour au fur et à mesure, donc lorsque vous exécutez un pipeline plus important, vous verrez la progression représentée ici. Mais parce que celui-ci est si petit, il s'exécute essentiellement immédiatement.

## 1.2.2. Trouver la sortie et les logs dans le répertoire work

Maintenant, lorsque vous exécutez un pipeline Nextflow, chacun de ces **process** est assemblé, et chaque **process**, comme je l'ai dit avant, peut générer des tâches une ou plusieurs. Donc dans ce cas, nous avions une seule tâche de ce **process**. Il s'est juste exécuté une fois et cela a été fait sous ce _hash_ de tâche.

Nextflow ne traite pas directement les fichiers dans votre répertoire de travail, il crée un dossier spécial appelé work. Et si je fais « ls », nous verrons qu'il est apparu ici : _work_, et à l'intérieur il y a des sous-répertoires pour chaque tâche qui s'exécute. Et cela correspond à ce hash. Donc vous pouvez voir si je vais à « ls work/c4 », et ensuite c'est tronqué, mais il commence par 203, et c'est le répertoire de travail, qui a été créé par ce **process** lorsque nous avons exécuté le pipeline. Et vous pouvez le voir sur le côté également.

Lorsque je liste ces fichiers, vous pouvez voir que le fichier output.txt a été généré. Vous pouvez le voir ici également. Et il y a un tas de fichiers cachés, qui n'apparaissent pas avec mon « ls » normal.

Si je clique sur output.txt, effectivement, nous avons notre sortie. Fantastique. Donc le pipeline a fonctionné.

Cela peut sembler beaucoup de code standard pour exécuter ce qui était essentiellement un script bash d'une ligne, mais cela aura plus de sens à mesure que nos **process** deviendront plus compliqués. Et ce répertoire work avec Nextflow et ces fichiers, qui sont créés, est vraiment l'épine dorsale de ce qui rend Nextflow si puissant.

Chaque tâche, chaque élément d'un pipeline est isolé de toutes les autres tâches. C'est reproductible. Ils n'entrent pas en conflit les uns avec les autres, et tout peut s'exécuter en parallèle. C'est en fait une très belle façon une fois que vous vous y habituez à cause de cette isolation que vous pouvez entrer et voir exactement ce qui s'est passé pour une seule tâche et déboguer.

Jetons un coup d'œil rapide à ces autres fichiers dans le répertoire work. De haut en bas, nous avons un fichier appelé _.command.begin_. Celui-ci est vide. C'est juste ce qu'on appelle un fichier sentinelle, créé par Nextflow disant, d'accord, je commence la tâche. Rien d'intéressant là.

Ensuite il y a _.command.error_, _.command.log_ et _.command.out_. Ce sont toutes des sorties de la commande bash ou de ce script qui s'est exécuté. C'est la sortie d'erreur standard. C'est la sortie standard, et c'est les deux combinées telles qu'elles sont sorties. Donc vous obtenez l'ordre logique.

D'accord, ceux-ci étaient tous vides pour cela également, donc pas très intéressant, mais les choses deviennent plus intéressantes lorsque vous arrivez à _.command.run_.

C'est généralement un script très long. Et c'est ce que Nextflow exécute réellement. Si vous entrez ici, vous commencerez à voir toute la logique interne de Nextflow et voir ce qu'il fait et comment il exécute votre **process**. Cela dépendra de l'endroit où vous exécutez, si nous exécutons localement ou si nous le soumettons comme un job à SLURM, auquel cas nous aurons des en-têtes SLURM en haut. Toutes ces différentes configurations.

Généralement, vous n'avez pas vraiment besoin de regarder dans ce fichier, cependant. Il est autogénéré par Nextflow et il n'y a rien de vraiment particulièrement unique à votre pipeline, qui s'y trouve. Mais c'est vraiment le cœur de ce qui s'exécute.

Le suivant est beaucoup plus intéressant. _.command.sh_ est le script généré, qui provient de votre **process**, et ici vous pouvez voir que Nextflow a ajouté l'en-tête Bash, puis il a exécuté notre commande, qui était dans notre bloc script.

Et c'est tout ce que fait le fichier _.command.run_, il exécute simplement ce fichier _.command.sh_.

C'est vraiment utile, c'est celui que vous finissez généralement par regarder le plus lorsque vous essayez de déboguer quelque chose et de vérifier que la logique de votre pipeline Nextflow fait ce que vous attendez qu'il fasse.

Enfin, nous avons un fichier appelé _.exitcode_, et cela capture simplement le code de sortie d'une tâche, qui dans ce cas a réussi. Donc le code de sortie était zéro.

Si quelque chose ne va pas, vous manquez de mémoire ou autre chose et cela échoue, alors c'est très utile pour comprendre ce qui n'a pas fonctionné.

## 1.3. Exécuter à nouveau le workflow

Une dernière chose à comprendre à propos des répertoires work est que si je continue à exécuter ce pipeline de manière répétée, donc si je fais _« nextflow run hello-world.nf »_, cela va faire exactement la même chose, mais cette fois il aura un nouvel id de tâche. Vous pouvez voir que ce hash ici est différent, et maintenant si je regarde dans work, il y a deux répertoires hash. Et ceux-ci sont, encore une fois, séparés l'un de l'autre.

Donc chaque fois que vous exécutez un workflow Nextflow, à moins que vous n'utilisiez le resume, qui utilise le cache, nous y reviendrons plus tard, il va réexécuter ces **process** dans de nouveaux répertoires work, qui sont séparés les uns des autres. Vous n'aurez aucune collision de noms de fichiers, vous n'aurez aucun problème comme ça. Tout est isolé et propre.

Et si nous entrons dans ce répertoire, vous pouvez voir tous les mêmes fichiers et le même _output.txt_, qui a été recréé à partir de zéro.

## 2. Publier les sorties

D'accord, c'est génial pour Nextflow lui-même, pendant qu'il exécute votre pipeline pour que toutes les choses soient séparées les unes des autres et propres et puissent être gérées.

Mais ce n'est pas super utile si vous êtes une personne essayant d'explorer vos résultats. Vous ne voulez pas vraiment fouiller dans des milliers et des milliers de répertoires work différents en essayant de trouver vos fichiers de résultats. Et vous n'êtes pas vraiment censé·e le faire. Les répertoires work ne sont pas censés être l'état final de l'endroit où vos fichiers sont créés.

Nous faisons cela en publiant nos fichiers.

## 2.1.1. Déclarer la sortie du process sayHello

Donc si je retourne à notre script, nous allons travailler dans notre bloc **workflow** ici. Nous allons lui dire quels fichiers attendre, quels fichiers nous intéressent, puis nous allons créer un nouveau bloc en dessous appelé le bloc output.

C'est la nouvelle syntaxe, qui est venue avec l'analyseur syntaxique et sera par défaut dans la version 26.04 de Nextflow. Donc si vous avez utilisé Nextflow un peu avant, c'est l'une des choses qui est nouvelle.

Donc nous avons le bloc main, et ensuite je vais dire publish et je vais dire à Nextflow quoi attendre de la publication. Nous allons l'appeler _first_output_, et nous allons l'appeler, _sayHello.out_.

J'ai accidentellement fait une faute de frappe là, mais c'est une bonne occasion de souligner également certaines des fonctionnalités de l'extension Nextflow VS Code. Vous pouvez voir que tout de suite elle m'a donné une petite ligne ondulée rouge en dessous disant que quelque chose ne va pas. Et si je survole, elle va me dire que cette variable n'est pas définie. Je ne sais pas ce que c'est.

C'est assez évident dans ce cas, j'ai fait une faute de frappe. Je voulais taper, sayHello, et ensuite la ligne ondulée disparaît.

Maintenant c'est violet. L'analyseur syntaxique Nextflow sait que c'est un **process** et quand je survole, il me donne une représentation réduite de ce à quoi ressemble ce **process**. Donc je peux voir très rapidement d'un coup d'œil qu'il ne prend aucune entrée et qu'il nous donne cette sortie. Donc travailler dans VS Code avec cette extension vous donne beaucoup d'informations contextuelles pendant que vous écrivez du code.

Notez que nous pouvons faire référence à la sortie de ce **process** avec la syntaxe _.out_. Et pour le moment nous pouvons appeler cela comme nous voulons, c'est juste un nom de variable arbitraire.

## 2.1.2. Ajouter un bloc output: au script

Là où cela devient important, c'est lorsque nous faisons notre nouveau bloc ici, et c'est en dessous du bloc **workflow** maintenant, nous ne sommes plus à l'intérieur du **workflow**. Accolades à nouveau. Et c'est là que nous disons simplement à Nextflow où mettre tous les fichiers, qui sont créés par le workflow.

Maintenant je vais prendre ce nom de variable, que j'ai créé ici, et je vais le mettre là et mettre des accolades pour cela. Et je vais dire à Nextflow d'utiliser un path. Oups. Path, entre guillemets. Et je vais utiliser un point. Cela dit simplement à Nextflow de mettre le fichier à la racine du répertoire results. Donc pas de sous-répertoires ou quoi que ce soit.

Essayons d'exécuter à nouveau notre workflow. Si je fais _« nextflow run hello-world.nf »_, alors j'espère que cela devrait avoir l'air essentiellement exactement pareil. Rien n'a vraiment changé avec Nextflow ici. Il exécute les mêmes choses. Il les fait juste dans des répertoires work à nouveau.

Mais maintenant si je fais _« ls results/ »_, vous verrez qu'il y a un nouveau répertoire ici qui a été créé appelé results, qui est le répertoire de base par défaut pour la publication du workflow. Et dedans il y a un fichier appelé _output.txt_.

Si je fais _« ls -l results »_, vous verrez que c'est en fait un lien symbolique vers le répertoire work. Donc ce n'est pas un vrai fichier, il est lié au répertoire work et il a collecté tous les fichiers là pour nous.

## 2.2. Définir un emplacement personnalisé

« Results » est le nom par défaut pour ce chemin. Si je réexécute le workflow, et cette fois je fais _dash_ trait d'union simple, c'est parce que c'est une option Nextflow de base. _« -Output-dir **my** results »_. Je pourrais aussi juste faire _« -o »_ pour faire court. Alors cela va définir un répertoire de base différent pour l'endroit où les fichiers sont stockés et encore une fois, ici en haut dans _myresults/_, maintenant nous avons un _output.txt_.

C'est génial, mais nous ne voulons probablement pas tous les fichiers juste à la racine. Nous voulons une certaine organisation, donc nous pouvons également créer un sous-répertoire ici appelé comme nous voulons. Disons _« path 'hello_world' »_, et je réexécute simplement cela. _« nextflow run hello-world.nf »_. Cela devrait aller dans le répertoire results dans un sous-répertoire et effectivement, maintenant sous results ici en haut nous avons _hello_world/_ et nous avons _output.txt_.

Chose importante à noter, l'ancien fichier _output.txt_ est toujours là. Le répertoire results n'est pas effacé lorsque vous faites cela. Seuls les nouveaux fichiers sont copiés dedans. Ils écraseront les fichiers qui sont déjà là s'ils ont le même nom de fichier, mais ils n'effaceront pas les anciens. Donc vous devez être un peu prudent·e lorsque vous réexécutez des pipelines. Si vous ne voulez pas qu'ils soient par-dessus les fichiers qui sont déjà là. Assurez-vous d'utiliser un répertoire vide.

## 2.3. Définir le mode de publication sur copy

D'accord, j'ai mentionné que ces fichiers sont des liens symboliques, donc si je fais _« ls -l results/hello_world/ »_, vous pouvez voir qu'il fait un lien symbolique vers le répertoire work. C'est généralement une bonne chose si vous travaillez sur quelque chose comme HPC, et ce sont des fichiers vraiment énormes et vous ne voulez pas les dupliquer, parce que cela signifie que les fichiers ne sont stockés qu'une seule fois sur le système de fichiers.

Cependant, cela signifie que si vous supprimez le répertoire work : si je fais _« rm -r work »_ et efface tous ces fichiers intermédiaires qui ont été créés. Maintenant, si j'essaie de lire ce fichier _« results/hello_world/ »_. Il va pointer comme un lien symbolique vers un fichier qui n'existe plus et les données sont perdues pour toujours et sont irrécupérables, ce qui n'est peut-être pas génial.

Donc généralement nous, je dis que c'est une bonne pratique de copier les fichiers au lieu de faire des liens symboliques si vous le pouvez, parce que c'est plus sûr. Soyez juste conscient·e que cela utilisera deux fois plus d'espace disque à moins que vous ne supprimiez ces répertoires work.

Pour faire cela avec le bloc output, je vais aller au first output ici. J'ai défini le path avant et maintenant je vais définir le mode et vous pouvez voir pendant que je tape, l'extension VS code suggère des choses qu'elle connaît, c'est une directive output ici. Et je vais dire copy. J'appuie sur enregistrer.

Réexécutons le workflow. Il va créer les fichiers à nouveau, nouveau répertoire work.

Maintenant, si je vais à _« ls -l results/hello_world/ »_ vous pouvez voir que c'est un vrai fichier et ce n'est plus un lien symbolique, et Nextflow l'a copié. Bon à savoir. Donc path et mode sont des choses que vous vous retrouverez à écrire assez souvent.

Maintenant, bien sûr, c'est très simple. Nous allons rendre cela plus complexe et puissant au fur et à mesure, et vous verrez comment rendre ces choses dynamiques et pas trop verbeuses.

## 2.4. Note sur les directives publishDir au niveau du process

Maintenant, j'ai dit au début que c'est une forme de syntaxe assez nouvelle. Elle n'est disponible que dans les dernières versions de Nextflow au moment où j'enregistre ceci, et cela s'appelle Workflow Outputs.

Si vous utilisez cela, c'est génial. Cela déverrouille beaucoup d'autres fonctionnalités intéressantes dans Nextflow, telles que Nextflow Lineage pour aider à suivre l'héritage de ces fichiers au fur et à mesure qu'ils sont créés, et bientôt sera la valeur par défaut dans 26.04. Et à une date ultérieure dans le futur, ce sera la seule façon d'écrire vos workflows.

Cependant, comme nous sommes dans cette phase de transition en ce moment, vous pourriez bien voir des pipelines dans la nature, qui utilisent quelque chose appelé publishDir, qui est l'ancienne façon de le faire, et cela est défini non pas au niveau du workflow et de la sortie, mais cela est défini au niveau du **process**.

Et cette déclaration dit essentiellement la même chose. Elle dit, publier les fichiers de résultats dans un répertoire appelé results, et utiliser un mode copy. Donc vous pouvez voir que la syntaxe est très similaire. Mais lorsque vous écrivez de nouveaux pipelines maintenant, essayez de ne pas utiliser cette directive publishDir, même si vous la voyez, dans les résultats d'IA ou dans la documentation ou d'autres pipelines, parce que c'est l'ancienne façon de le faire.

En 2026, nous devrions tous utiliser les workflow outputs.

Tout cela est documenté, si vous faites cela et que vous avez utilisé Nextflow avant, vous pouvez aller sur la documentation Nextflow ici, nextflow.io/docs/. Et si je descends jusqu'à tutorials, il y a un tutoriel appelé _Migrating to Workflow Outputs_.

Il est vraiment bon. Il passe en revue toute la syntaxe, comment elle est équivalente à l'ancienne syntaxe, pourquoi nous l'avons changée, et, avoir une chronologie et tout. Et il passe en revue tous les différents scénarios avec des tonnes et des tonnes d'exemples. Donc vous pouvez facilement convertir le code Nextflow existant vers la nouvelle syntaxe.

## 3.1. Modifier le process sayHello pour attendre une entrée variable

D'accord, donc nous avons notre script simple, qui exécute un **process**, crée un fichier, dit à Nextflow que c'est une sortie, puis nous disons à Nextflow où enregistrer ce fichier. C'est un bon début.

Mais ce serait plus intéressant si tout n'était pas codé en dur. Donc ensuite, réfléchissons à comment dire à Nextflow que ce **process** peut prendre une entrée variable, qui est quelque chose que nous pouvons contrôler au moment de l'exécution lorsque nous lançons un workflow.

Nous devons faire quelques choses différentes pour que cela se produise.

Premièrement, nous devons dire à ce **process** qu'il peut accepter une variable d'entrée et nous tapons _input_ ici comme un nouveau bloc de déclaration. Et nous allons appeler cela _« val greeting »_.

La partie val est l'équivalent d'un path en bas. Elle dit à Nextflow que c'est une variable, comme une chaîne de caractères dans ce cas. Et si vous survolez à nouveau, elle vous dit de l'extension ce que cela signifie.

Ensuite, nous allons dire à Nextflow quoi faire avec cela. Ce n'est pas suffisant de simplement dire qu'il y a une variable. Vous devez dire dans le script comment utiliser cette variable. Et donc je vais me débarrasser de cette chaîne codée en dur ici, et je vais mettre une variable.

Je vais rapidement le faire sans accolades juste pour vous montrer que c'est autorisé, et c'est l'ancienne façon de le faire. Mais maintenant avec la nouvelle syntaxe, nous recommandons vraiment de le mettre entre accolades comme ceci, et cela rend vraiment clair que cela est interpolé par Nextflow ici.

Génial. Donc _« input greeting »_ va dans _$\{greeting\}._ La dernière chose est que nous devons dire à Nextflow au niveau du workflow que ce **process** prend maintenant une entrée. Et pour faire cela, nous allons essentiellement lui donner une variable.

## 3.2. Configurer un paramètre en ligne de commande pour capturer l'entrée utilisateur·trice

Nous pourrions le coder en dur à nouveau, comme Hello World, et cela fonctionnerait bien, mais évidemment cela ne nous donne pas vraiment d'avantage. Nous voulions pouvoir configurer cela au moment de l'exécution, donc nous voulons pouvoir le faire sur la CLI, lorsque vous lancez Nextflow.

Et la façon dont nous faisons cela est un concept Nextflow spécial appelé _params_. Nous allons appeler cela _params.input_.

Ce que cela fait, c'est qu'il expose cette variable input sur la CLI et c'est là que nous utilisons un double tiret lorsque nous lançons Nextflow.

Je peux appeler cela comme je veux, je peux l'appeler _hello, greeting_. Peu importe. Tout ce que je fais là sera exposé comme une option CLI lorsque nous lançons un pipeline. Et c'est un vrai tour de magie de Nextflow parce que cela signifie que vous pouvez construire votre script de workflow très rapidement avec ces paramètres, et vous construisez essentiellement une CLI personnalisée pour votre pipeline, ce qui facilite vraiment la personnalisation de différentes options à la volée lorsque vous lancez.

Donc. Essayons-le. Retournons à notre terminal. Nous avons notre commande _« nextflow run »_ ici. Et maintenant je vais faire _« --input »_, qui correspond au _« params.input »_ que nous avons vu avant. Je pense que dans la documentation c'est en français. Geraldine aime parler français. Je vais le faire en suédois parce que je vis en Suède. donc je vais dire, « _Hej Världen_ » et appuyer sur entrée.

On peut utiliser des guillemets simples ou doubles, cela affecte juste comment Bash l'interprète.

Il exécute le pipeline Nextflow exactement de la même manière. Vous pouvez voir que le répertoire de travail et tout est le même. Mais maintenant si je vais jusqu'à _« results/hello_world/output »_. Nous pouvons voir notre joli suédois ici à la place.

Donc nous avons passé dynamiquement une entrée d'une CLI à un paramètre. Nous avons passé cela comme une entrée à un **process** et le **process** a interprété cela et l'a mis dans un bloc script, qui a ensuite changé dynamiquement la sortie de ce résultat de script. Plutôt cool.

Une logique assez complexe avec très peu de syntaxe ici. Et vous pouvez espérons-le voir comment cela commence maintenant à évoluer. Et c'est comme ça que nous construisons vraiment la logique et la personnalisabilité de nos pipelines dans le script Nextflow.

## 3.4. Utiliser des valeurs par défaut pour les paramètres de ligne de commande

D'accord, c'est génial. Le problème maintenant cependant est que chaque fois que j'exécute ce pipeline, je dois faire dash, input pour qu'il s'exécute.

Si j'essaie d'exécuter sans ce paramètre, maintenant Nextflow va générer une erreur disant qu'il avait besoin de ce paramètre et qu'il n'a pas été défini. et donc il ne savait pas quoi faire.

C'est une nouvelle chose cool, au fait. Dans le passé, Nextflow aurait juste exécuté avec une chaîne vide, et vous auriez eu toutes sortes d'erreurs bizarres, qui auraient été difficiles à comprendre. Mais dans le nouvel analyseur syntaxique Nextflow, il est un peu plus prudent et il vous le dit tout de suite.

Donc nous ne voulons pas toujours spécifier chaque option. C'est une bonne pratique de spécifier des valeurs par défaut sensées. Alors comment faisons-nous cela dans notre script ?

Vous remarquerez que lorsque nous avons écrit cela, nous avons juste mis _params.input_ directement là où nous l'utilisons. Donc la solution évidente est que nous définissons une valeur par défaut, et nous faisons cela en haut du script ici dans un bloc params spécial dans le workflow. C'est dans le script workflow ici.

Encore une fois, une nouvelle syntaxe ici, donc faites attention. C'est vraiment cool. Nous avons le nom du paramètre, qui sera attendu ici.

Et ensuite après ce caractère deux-points, nous définissons un type de la variable. Vous n'êtes pas obligé·e de faire cela, vous pouvez juste le laisser vide, mais c'est vraiment bien. Cela dit à Nextflow que nous attendons une chaîne de caractères et de la traiter comme telle.

Si nous voulons un nombre à la place, par exemple, nous pourrions écrire float, et cela dirait que nous voulons un nombre à virgule flottante. Et si nous essayons d'exécuter avec cela, alors il générera une erreur. Si nous lui donnons une chaîne, qui n'est pas un float. Et il le passera également comme tel. Comme si nous faisons string, alors il sait que c'est une chaîne. Et même s'il a des zéros en tête et est entièrement numérique, il le passera quand même comme une vraie chaîne.

Donc cette sécurité de type est une fonctionnalité très nouvelle de Nextflow, mais vraiment puissante pour rendre votre code plus sûr à écrire et à exécuter.

Ensuite après cela nous avons un symbole égal puis la valeur par défaut ici. Nextflow a été écrit à Barcelone à l'origine, donc il semble approprié que nous ayons un peu d'espagnol ici, _« Holà mundo ! »_ comme valeur par défaut.

Bon je vais enregistrer ce script, retourner, exécuter le script à nouveau sans _--input_. Et cette fois il devrait s'exécuter et il créera notre nouveau fichier dans _results_. Et dans ce fichier maintenant il dit _« Holà mundo ! »_.

Ce n'est qu'une valeur par défaut cependant, donc cela ne signifie pas que nous ne pouvons toujours pas faire la même chose qu'avant. Si je retourne et trouve mon ancien script ici, _« Hej Världen »_, parce que je fais _--input_ sur la ligne de commande, cela écrasera cette valeur par défaut et utilisera cela à nouveau dans le fichier output.txt.

Donc cela dans le script n'est que la valeur par défaut que je définis.

Au fur et à mesure que nous construisons notre workflow pour qu'il soit plus complexe et inclue plus de paramètres, ce bloc params en haut du script commencera à tous les collecter en un seul endroit.

Et vous vous retrouvez avec cette symétrie assez agréable dans votre script, où vous avez effectivement toutes vos entrées de workflow ici et vos sorties de workflow en bas. Et il est très clair quelle est l'interface de votre workflow avec le monde extérieur. Donc vous pouvez récupérer un nouveau pipeline très rapidement avec la nouvelle syntaxe et comprendre comment l'utiliser.

Une dernière chose cool. Nous n'avons pas à définir une valeur par défaut avec cela. Si nous faisons params input mais ne définissons pas de valeur par défaut, alors cela dit à Nextflow que ce paramètre est requis, et encore une fois, le pipeline échouera à s'exécuter sans lui, mais il vous donnera un message d'erreur plus utile plutôt que quelque chose à propos d'être null.

Donc il dit que nous attendons que son input soit requis, mais il n'a pas été spécifié sur la ligne de commande. Très bien.

D'accord, donc j'espère que maintenant il est clair comment configurer votre pipeline Nextflow avec des entrées et des paramètres variables, comment définir la valeur par défaut, définir les types, cela pourrait être un booléen vrai faux flag ou un entier ou différents types ici. Comment les passer dans votre workflow, où cela passe, puis interpole dans votre **process**. Et ensuite vous savez également comment personnaliser ceux-ci sur la ligne de commande lorsque vous lancez Nextflow. Cela commence à avoir l'air plus intéressant que notre simple commande bash.

## 4. Gérer les exécutions de workflow

D'accord. Et ensuite ? Pour la dernière partie de ce chapitre, nous allons parler un peu de comment gérer toutes les différentes exécutions de workflow. Si vous regardez dans ma barre latérale ici et l'Explorateur sous work, vous verrez que j'ai exécuté un tas de pipelines différents et ces répertoires work deviennent assez longs, il y en a beaucoup.

Et l'autre chose est, comme je l'ai dit avant, chaque fois que je réexécute ce pipeline, il crée un nouveau jeu de répertoires work, et il réexécute tous les **process** à partir de zéro, ce qui est une bonne chose. C'est le comportement prévu. C'est reproductible et il régénère tout à neuf. Mais évidemment, si vous exécutez des **process** de très longue durée, c'est ennuyeux de toujours devoir démarrer votre pipeline depuis le début s'il a planté à mi-chemin, ou si vous changez quelque chose à la fin du pipeline.

## 4.1. Relancer un workflow avec -resume

Heureusement, Nextflow est vraiment bon pour savoir ce qui a été précédemment exécuté et ce qui est disponible, et pour réutiliser ces anciens résultats c'est très simple. Nous ajoutons juste un nouveau flag à la fin de la commande _« -resume »_.

Maintenant, notez qu'il y a deux tirets sur input parce que c'est le paramètre. Il n'y a qu'un seul tiret sur resume parce que c'est une option Nextflow de base.

Cela fait trébucher les gens tout le temps, même si vous utilisez Nextflow depuis longtemps. Donc rappelez-vous toujours un ou deux tirets. Cela dépend si c'est une option Nextflow de base.

D'accord, donc maintenant je fais _-resume_ et j'exécute exactement le même workflow à nouveau. Et cette fois cela devrait avoir l'air à peu près exactement pareil avec une différence clé.

Dans la sortie ici, vous pouvez voir que les résultats ont été mis en cache. Et en fait, ce hash de tâche ici est exactement le même que l'exécution précédente, et il a juste réutilisé ce répertoire work dans son intégralité. Les entrées et les sorties et le script étaient tous non modifiés. Et donc il prend juste ce fichier de là et s'il y a des étapes en aval dans le **process**, il les passerait à l'étape suivante dans le pipeline.

Donc il exécute toujours le pipeline entier du début à la fin, mais il utilise des résultats mis en cache pour chacune de ces tâches, là où il le peut.

Maintenant, lorsque vous faites _-resume_, il reprend juste la dernière exécution de pipeline dans votre répertoire de travail, quelle qu'elle soit. Mais vous pouvez en fait reprendre à partir de n'importe quelle exécution précédente que vous avez faite là. Et nous en avons fait pas mal maintenant.

## 4.2. Inspecter le journal des exécutions passées

Pour les regarder toutes, nous pouvons faire _« nextflow log »_ au lieu de _« nextflow run »_, et cela nous donnera une belle sortie montrant toutes ces différentes... J'ai besoin de rendre mon écran un peu plus petit pour que nous puissions le voir, toutes ces différentes exécutions quand nous les avons faites, l'id de session, la commande et tout.

Et nous pouvons regarder ici et nous pouvons prendre le nom d'exécution de n'importe laquelle de celles-ci et ensuite reprendre une de ces spécifiques. Donc je peux retourner et je peux reprendre celle appelée _hungry_ekeblad_. Et je mets juste cela après le _resume_.

Si vous êtes curieux·se, au fait, tous ces adjectifs et noms de scientifiques sont dans le code source Nextflow. C'est un très bon moyen d'obtenir votre toute première pull request à Nextflow en allant le trouver et en ajoutant votre scientifique préféré·e.

Et de toute façon, donc j'ai fait cela et il est retourné et il a regardé les résultats mis en cache de cette exécution de workflow, a réalisé qu'il pouvait toujours les réutiliser, et il l'a fait. Donc j'ai obtenu les résultats mis en cache à nouveau.

## 4.3. Supprimer les anciens répertoires work

C'est génial. Qu'en est-il si je veux nettoyer ces répertoires work ? Il y en a des tonnes ici. Il y a des tonnes de fichiers. Peut-être que je sais avec certitude que je veux reprendre à partir des deux dernières exécutions de pipeline, mais je ne me soucie pas de toutes celles d'avant.

Alors je peux en choisir une ici et je peux utiliser une autre commande Nextflow, qui est _« nextflow clean »_, et je peux faire _« nextflow clean »_, je vais faire _« -before »_, et le nom d'exécution particulier, qui dans ce cas était _reverent_pike_ et je vais faire _« -n »_, qui dit à Nextflow de juste faire une exécution à sec. Donc il me dit juste ce qu'il supprimera. Sans réellement faire quoi que ce soit, donc il supprimerait ces répertoires work.

Cela semble sensé. Donc je vais faire la même commande à nouveau, mais au lieu de _« -n »_ je vais faire _« -f »_ pour réellement faire le nettoyage. Et cette fois il a réellement supprimé tous ces répertoires. Et si j'entre et regarde les répertoires work, c'est maintenant beaucoup plus léger. Fantastique.

Donc c'est comme ça qu'on nettoie tous vos répertoires work locaux d'une manière assez sûre sans complètement détruire le cache. Donc vous pouvez toujours reprendre si vous voulez.

Si jamais vous oubliez ce que sont ces flags pour chaque commande Nextflow vous pouvez faire _« nextflow help »_, puis le nom de la commande. Donc si je fais _« nextflow help clean »_, vous pouvez voir toutes les différentes options : _-after, -before, -but_, toutes les différentes façons de configurer ce comportement de nettoyage. Plutôt cool.

## À retenir

D'accord, c'est la fin de la partie un de Hello Nextflow. C'est un début assez intense pour la formation, mais j'espère que maintenant vous avez une assez bonne compréhension de ce à quoi ressemble un script Nextflow ; avec différentes parties clés, les **process**, les **workflow**, les sorties, et les paramètres. Vous savez comment les configurer avec des remplacements de base depuis la ligne de commande, comment créer un bloc input dynamique avec un script dynamique et vous savez comment gérer toutes vos exécutions de charge de travail : voir ce que vous avez déjà exécuté, reprendre, nettoyer. Il y a beaucoup de choses. Vous avez parcouru un long chemin. Donc si vous voulez faire une pause et faire une petite promenade et prendre une tasse de thé, c'est probablement le bon moment maintenant. Vous l'avez mérité.

À partir d'ici, nous construisons essentiellement sur cette fondation. Comment pouvons-nous rendre cela plus complexe, plus puissant ? Comment pouvons-nous le rendre plus flexible ? Faire les choses que nous voulons faire notre analyse à grande échelle.

## Quiz

Maintenant si vous descendez jusqu'à la partie un, hello world, sur la page web vous verrez un petit quiz et c'est quelque chose de nouveau que nous avons fait pour cette version de la formation Nextflow. Et vous pouvez parcourir et vous tester pour vérifier que vous avez compris tout le matériel que nous avons fait dans ce chapitre.

Cela ne nous est pas envoyé ou quoi que ce soit, c'est juste stocké dans votre navigateur. Donc nous ne savons pas quelles sont vos réponses, mais c'est juste une petite auto-vérification pour vous assurer que vous n'avez rien manqué ou mal compris quoi que ce soit. Et vous pouvez l'essayer autant de fois que vous le souhaitez.

Si vous êtes comme moi, peut-être que vous voulez rester dans le terminal dans votre instance VS Code, auquel cas vous pouvez taper la commande _quiz_ puis juste lui dire sur quel chapitre vous êtes. Donc nous faisons _« Hello World »_, et ensuite vous pouvez faire exactement les mêmes questions de quiz, qui sont dans le navigateur web, mais juste dans votre terminal.

Cool. D'accord. J'espère que vous appréciez cela. Amusez-vous un peu et, nous vous verrons dans le prochain chapitre dans juste une minute pour parler de tous les canaux Nextflow.
