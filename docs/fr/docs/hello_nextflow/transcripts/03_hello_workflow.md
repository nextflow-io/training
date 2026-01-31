# Partie 3 : Hello Workflow - Transcription

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Traduction assistée par IA - [en savoir plus et suggérer des améliorations](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

<div class="video-wrapper">
  <iframe width="560" height="315" src="https://www.youtube.com/embed/zJP7cUYPEbA?si=Irl9nAQniDyICp2b&amp;list=PLPZ8WHdZGxmXiHf8B26oB_fTfoKQdhlik" title="YouTube video player" frameborder="0" allow="accelerometer; autoplay; clipboard-write; encrypted-media; gyroscope; picture-in-picture; web-share" referrerpolicy="strict-origin-when-cross-origin" allowfullscreen></iframe>
</div>

!!!note "Notes importantes"

    Cette page présente uniquement la transcription. Pour les instructions complètes étape par étape, retournez au [matériel de formation](../03_hello_workflow.md).

    Les numéros de section indiqués dans la transcription sont fournis à titre indicatif uniquement et peuvent ne pas inclure tous les numéros de section du matériel.

## Bienvenue

Bonjour, bienvenue dans la partie trois de la formation « Hello Nextflow ».

Ce chapitre s'intitule « Hello Workflow ».

Dans le chapitre deux, nous avons construit un workflow simple d'un seul process, mais en réalité, les pipelines sont utiles car ils peuvent enchaîner plusieurs étapes d'analyse ensemble.

Dans ce chapitre, nous allons prendre cet exemple initial et l'étendre pour le rendre un peu plus réaliste.

Nous allons ajouter quelques étapes supplémentaires et nous allons voir comment utiliser les channels pour connecter ces étapes.

Nous allons examiner plusieurs tâches, qui peuvent se regrouper en un seul process, et nous allons examiner les processes qui peuvent avoir plusieurs entrées et plusieurs sorties.

D'accord, commençons.

Alors, commençons. Comme auparavant. Allons sur training.nextflow.io. Hello Nextflow, chapitre trois. Hello Workflow. Et ouvrons notre espace de travail. J'ai nettoyé tous mes fichiers de travail des chapitres précédents et je vais ouvrir Hello Workflow.

Maintenant, c'est le même fichier sur lequel nous avons travaillé jusqu'à présent, donc cela devrait être familier. Nous avons notre process say hello. Nous avons notre params.greeting avec son fichier greetings CSV, et nous avons notre workflow en bas, qui charge ce fichier CSV, crée le channel et le transmet à notre process.

## 0. Échauffement : Exécuter hello-workflow.nf

Si vous voulez, nous pouvons essayer cela et vérifier qu'il fonctionne comme prévu. Ouvrez un terminal pour nextflow run hello workflow nf et appuyez sur entrée.

D'accord, excellent. Nos trois processes s'exécutent. Nous avons notre répertoire results avec nos trois sorties. Bonjour. Hello. Holà. Alors fermons ces fichiers, fermons le terminal, revenons au script.

## 1. Ajouter une deuxième étape au workflow

D'accord. Pour notre exemple, nous restons basiques et nous essayons de rester indépendants du domaine. Notre deuxième process va donc simplement manipuler ces chaînes de caractères, ces mots, de manière simple. Nous allons utiliser la commande Unix translate pour prendre ces fichiers et les mettre tous en majuscules. Nous faisons cela avec la commande « tr ».

## 1.1. Définir la commande de mise en majuscules et la tester dans le terminal

Nous pouvons essayer cela simplement dans le terminal bash, et voir si cela fonctionne. Donc vous faites echo, Hello World, puis vous passez cela avec le caractère pipe à tr, et nous lui donnons un motif de reconnaissance, a à z et ce en quoi il doit traduire. A à Z en majuscules.

C'est très simple car cela traite littéralement les caractères A à Z. Donc cela ne fonctionnera pas sur quoi que ce soit d'accentué ou similaire. Mais pour les besoins de l'exemple, vous devriez comprendre l'idée.

Je vais appuyer sur entrée et cela affiche dans le terminal, HELLO WORLD en majuscules. Et tout comme avant, nous pourrions rediriger cela vers un fichier si nous le voulions. Outfile.

D'accord. Nettoyons cela.

## 1.1. Écrire l'étape de mise en majuscules en tant que process Nextflow

Retournons à notre script et écrivons un nouveau process pour gérer cette commande bash. Je vais copier le process précédent, le coller en dessous, et l'appeler convert to upper. Pour majuscules. Je vais utiliser le même publishDir results, mais je vais faire quelques modifications ici. Au lieu de prendre un val, je vais prendre un path input file, et je vais avoir un préfixe ici upper, pour que nos fichiers de sortie n'écrasent pas la sortie. Et je vais utiliser le nom de variable de l'entrée. Et ensuite je vais modifier le script ici en bas, et à la place je vais utiliser cat sur le fichier d'entrée et tout comme nous l'avons fait dans Bash TR, a-z, upper input file .txt. D'accord, cliquons sur enregistrer.

## 1.2. Ajouter un appel au nouveau process dans le bloc workflow

Maintenant, si je descends, nous devons réellement appeler ce process. Ajouter simplement le process dans le script ne suffit pas. Nous devons indiquer à Nextflow que nous devons exécuter ce process et où le faire.

Donc je vais écrire ici, convert to upper et

d'accord, nous obtenons une erreur ici disant qu'il attend un argument. Effectivement, nous devons passer quelque chose à ce process pour qu'il ait effectivement quelque chose à faire.

## 1.3. Passer la sortie du premier process au deuxième process

Ce que nous allons faire, c'est prendre la sortie de ce process. Donc je prends le nom, say hello, et quand je fais dot out.

Pour un exemple simple comme celui-ci, où nous avons un process qui n'a qu'une seule sortie et nous passons cela à un nouveau process, donc il a une entrée, cela devrait être tout ce dont nous avons besoin. Donc je vais cliquer sur enregistrer, ouvrir le terminal, et essayons d'exécuter cela à nouveau.

## 1.4. Exécuter à nouveau le workflow

Maintenant, je n'ai pas nettoyé mon répertoire work de la dernière fois que j'ai exécuté ce workflow. Je vais l'exécuter à nouveau et je vais utiliser cela comme une opportunité pour montrer comment fonctionne la mise en cache partielle. Donc si je fais simple tiret resume. Avec un peu de chance, il devrait réutiliser les sorties de ce premier process, qui étaient exactement les mêmes que la dernière fois que j'ai exécuté. Mais maintenant nous avons un nouveau process ici qui n'a jamais été exécuté auparavant, qui s'exécute à partir de zéro. Et effectivement, vous pouvez voir que le premier process a utilisé les sorties en cache, et la deuxième sortie a exécuté trois sur trois. Vous pouvez également voir que nous avons maintenant nos deux processes ici, notre premier process, say hello, exécuté trois fois, et notre deuxième process convert to upper exécuté trois fois.

Si j'exécute cela à nouveau, pour rappel, avec -ansi-log false, nous devrions voir que six tâches de process différentes s'exécutent, trois pour chacune d'elles. Donc cela fait exactement ce que nous espérions. Le premier process s'exécute trois fois, passant ces sorties à un deuxième process, qui s'exécute ensuite trois fois.

Alors jetons un coup d'œil dans le répertoire work et voyons comment Nextflow gère ces entrées de fichiers. Si je prends ce répertoire de hachage ici du deuxième process, nous pouvons utiliser à nouveau la commande tree avec -a juste pour regarder ces fichiers. Vous pouvez voir ici que nous avons notre fichier d'entrée, qui est le fichier Bonjour-output.txt, et c'est en fait un lien symbolique. C'est ce que cette flèche nous montre, et elle pointe vers le fichier dans le répertoire work précédent.

Cela a du sens. Nextflow gère l'exécution de chaque tâche dans son propre répertoire encapsulé, donc il est complètement autonome. Cependant, il doit fournir les fichiers des étapes précédentes comme entrée. Plutôt que de sortir du répertoire work pour obtenir ces fichiers, Nextflow les organise dans le répertoire work.

Si nous avons un système de fichiers partagé comme ici, il le fait en utilisant un lien symbolique pour ne pas utiliser d'espace de fichier supplémentaire. Si nous utilisons le stockage cloud avec des buckets dans différents emplacements, il récupérerait ces fichiers et les copierait réellement dans le répertoire work.

Regardons le fichier command sh. Si je fais code work, command sh, vous pouvez voir, effectivement, il accède à ce fichier depuis le répertoire local. Donc tout est très autonome et propre.

Nous pouvons également vérifier le répertoire results et nous assurer que ces fichiers ont été correctement sortis. Et effectivement, dans results, nous pouvons voir tous les fichiers de sortie du premier process et tous les fichiers de sortie du second. Et ils sont tous en majuscules comme nous l'espérions.

C'est là que la puissance de Nextflow commence à briller. Avec un code très minimal, Nextflow a géré l'exécution en parallèle de ces tâches avec une encapsulation propre dans des répertoires work séparés et l'organisation des fichiers d'entrée et de sortie et la publication des fichiers, tout automatiquement pour nous, directement prêt à l'emploi. Donc vous pouvez voir comment, à mesure que nous augmentons la complexité de nos workflows d'analyse, cette fonctionnalité est vraiment, vraiment précieuse.

## 2. Ajouter une troisième étape pour collecter toutes les salutations

D'accord. Ces étapes étaient un pour un. Nous avions une sortie du premier process allant à une entrée pour le deuxième process. Ensuite, nous allons parler de la façon de collecter ces différentes sorties dans une seule tâche de process, ce qui est encore une fois, une chose très courante à faire. Alors ouvrons rapidement le terminal et faisons un test à sec de cela.

## 2.1. Définir la commande de collecte et la tester dans le terminal

Je vais tricher et copier l'exemple de code bash du matériel de formation et simplement appuyer sur entrée.

Ce que nous pouvons voir ici, c'est que nous avons exécuté cette commande echo trois fois vers trois fichiers de sortie différents, que je peux voir ici. Et ensuite nous avons utilisé la commande cat pour afficher la sortie de chacun de ces trois fichiers différents, et rediriger cela vers un seul fichier collecté.

Et si je fais « cat COLLECTED-output », vous pouvez voir qu'il contient le contenu de ces trois fichiers différents, maintenant dans un seul fichier.

## 2.2. Créer un nouveau process pour effectuer l'étape de collecte

Alors voyons si nous pouvons reproduire la même chose dans notre pipeline Nextflow.

Remontons et créons un troisième process. Je vais copier celui précédent, et cette fois je vais l'appeler Collect Greetings.

Dans le terminal bash, nous l'avons appelé collected output txt. Donc je vais dire le même path output ici. Et je vais faire la redirection ici, donc c'est sauvegardé de la même manière.

D'accord. Nous devons changer ce qui se passe au début de cette commande, et nous devons réfléchir à ce qu'est le fichier d'entrée ici. En fait, ce process va prendre plusieurs fichiers d'entrée. Je vais garder path et je vais changer cela en une nouvelle variable appelée input files, au pluriel.

Je vais ensuite à nouveau, les cat comme nous l'avons fait dans notre script bash. Et je vais utiliser la variable ici.

Maintenant, vous pourriez penser que cela ne fonctionnerait pas. Nous avons vu précédemment des échecs où un tableau de chaînes ou un tableau de chemins a été passé à un process et cela a causé une erreur. Mais en fait, ici Nextflow va gérer cela automatiquement pour nous de la bonne manière. Il va prendre plusieurs fichiers d'entrée différents, et il va simplement afficher les différents chemins de fichiers ici.

Bien sûr, cela aide que la commande cat puisse prendre une série de noms de fichiers comme cela. Si j'utilisais une commande différente qui nécessitait un argument avant chaque chemin de fichier ou quelque chose comme ça, nous devrions avoir un peu plus de code ici et de logique pour pouvoir gérer l'itération de ces chemins de fichiers. Mais dans ce cas, cela devrait simplement fonctionner.

## 2.3. Ajouter l'étape de collecte au workflow

D'accord, descendons au workflow et ajoutons notre nouveau process. Collect greetings. Et encore une fois, prenons la sortie de convert to upper out. Enregistrons cela.

Essayons. nextflow run hello workflow.

D'accord, le workflow s'est exécuté, mais quelque chose est un peu étrange ici. Nous avons trois exécutions de la première étape, ce que nous attendons. Trois tâches pour la deuxième, mais nous avons également trois tâches à la fin alors que nous nous attendions à n'avoir qu'une seule tâche ici fusionnant toutes les sorties.

Si nous allons dans notre répertoire results. Nous voyons également que la sortie collectée n'a qu'une seule valeur plutôt que les trois. C'est parce que ce fichier de sortie a été écrasé trois fois avec trois valeurs différentes.

Cela a du sens car nous avons passé une sortie à une entrée ici de la même manière que nous l'avons fait à l'étape précédente.

## 2.4. Utiliser un opérateur pour collecter les salutations en une seule entrée

Nous avons donc besoin d'un opérateur ici pour prendre ce channel avec trois éléments et les réduire à un seul élément, pour que ce process final ne s'exécute qu'une seule fois.

Pour ce faire, nous allons utiliser l'opérateur collect. Je peux faire cela directement dans le workflow. Je peux faire .out et enchaîner un opérateur ici à la fin .collect.

Cliquez sur enregistrer. Et puis pour les besoins de cette formation, je vais également faire quelques opérateurs view comme nous l'avons fait avant, pour que nous puissions jeter un œil à ce channel avant et après avoir utilisé l'opérateur collect, pour que nous puissions comprendre ce qui se passe.

Je vais prendre ce channel, supprimer le collect et dot view greetings, puis je vais dupliquer cette ligne, ajouter l'opérateur collect. Et changer cela en after.

C'est séparé de l'endroit où nous appelons cela, mais ce n'est pas grave car nous utilisons les mêmes appels d'opérateurs sur le même channel de sortie.

D'accord, enregistrons et essayons-le dans le terminal. Je vais faire nextflow run. Hello, workflow. Réexécutons notre script.

D'accord. Cela a l'air mieux. Comme avant, nous pouvons voir que les deux premiers processes s'exécutent trois fois et maintenant notre process final ne s'est exécuté qu'une seule fois.

Si nous regardons ce qui a été affiché par l'opérateur view, ici en bas, nous avons dit before collect, qui est cette sortie ici, et c'est affiché trois fois. Et vous pouvez voir qu'il y a un seul chemin pour chacun d'eux. Et ensuite after collect, vous pouvez voir que nous avons ce tableau de trois chemins. Donc c'est comme prévu.

D'accord, vérifions le fichier results et voyons si c'est ce que nous attendons cette fois. Effectivement, il y a maintenant trois lignes dans le fichier - cela a réussi à concaténer ces trois sorties en un seul fichier de sortie. Fantastique.

D'accord, je vais nettoyer et passons à l'étape suivante. Et je vais supprimer ces instructions view juste pour garder les choses propres.

## 3. Passer plus d'une entrée à un process afin de nommer le fichier de sortie final de manière unique

D'accord. Jusqu'à présent, tous nos processes n'ont pris qu'une seule entrée. Nous allons maintenant faire un exercice où nous ajoutons plus d'une entrée à un process pour voir comment cela fonctionne. Pour ce faire, nous allons utiliser cet exemple collect greetings.

Chaque fois que j'ai exécuté le workflow, il a écrasé ce fichier dans le répertoire results, ce qui peut ne pas être ce que nous voulons.

## 3.1. Modifier le process de collecte pour accepter un nom défini par l'utilisateur pour le fichier de sortie

Donc pour cet exemple, nous allons passer un paramètre supplémentaire pour pouvoir personnaliser le nom du fichier de sortie.

Ajouter une deuxième entrée à un process est très simple. J'ajoute simplement une deuxième ligne dans le bloc input. Cette fois, ce sera une value, plutôt qu'un path, car nous voulons passer une chaîne et je vais l'appeler batch underscore name.

Je peux maintenant utiliser cette variable dans le bloc script, et je vais dire collected dash dollar batch name.

J'utilise des accolades ici autour du nom de la variable. C'est juste pour la séparer du reste de la chaîne, et ce n'est probablement pas nécessaire dans ce cas, mais je pense que cela rend la lecture plus facile.

D'accord. Enfin, n'oubliez pas de mettre à jour le path de sortie car maintenant le nom du fichier a changé, donc je vais faire la même chose et mettre le batch name dans la sortie du path comme prévu.

## 3.2. Ajouter un paramètre de ligne de commande batch

Nous devons maintenant passer un batch name de quelque part, et je vais créer un deuxième paramètre pour ce faire afin que nous puissions le faire en ligne de commande lorsque nous exécutons le workflow.

Donc je vais faire params batch name, et par défaut, appelons cela test batch. Maintenant je peux utiliser cette variable de paramètre spéciale en bas, là où nous appelons le process.

Et effectivement VS Code nous dit qu'il n'y a pas assez d'arguments pour ce process maintenant, et qu'il attend une deuxième entrée.

Simplement faire virgule et passer notre nouvelle variable et l'erreur disparaît.

Notez que l'ordre des entrées ici est vraiment important. La première entrée du process était le path, et la deuxième entrée est le name. Si je change l'ordre ici, je dois également changer l'ordre lorsque j'appelle le process. Sinon. Nextflow passera le mauvais channel à la mauvaise entrée.

## 3.3. Exécuter le workflow

D'accord, essayons et voyons si cela fonctionne. Faisons « nextflow run hello- workflow ». D'accord, il s'est exécuté comme avant. Regardons dans le répertoire results.

Effectivement, notre nom de fichier ici s'appelle maintenant « collected test batch output txt ». Fantastique.

Et maintenant voyons si nous pouvons écraser cela en exécutant à nouveau. Cette fois, je vais faire --batch_name pour correspondre à ce nom de variable de paramètre spécial ici. Et je vais l'appeler demo output.

Exécutons à nouveau le workflow et nous verrons si quelque chose se passe.

D'accord, nous avons maintenant un collected demo output .txt. Et parce que ce nom de fichier est différent de celui-là, il ne l'a pas écrasé. Les deux sont maintenant présents dans le répertoire results.

## 4. Ajouter une sortie à l'étape de collecte

D'accord, donc là nous avons montré comment donner plusieurs entrées à un process, mais qu'en est-il de plusieurs sorties ? Pour cet exemple, nous allons calculer le nombre de salutations qui sont traitées et sortir cela comme une sortie secondaire pour cette étape collect greeting.

## 4.1. Modifier le process pour compter et sortir le nombre de salutations

Nous allons faire une petite astuce ici. Les processes Nextflow ont ce bloc script avec une chaîne multi-lignes, et cela est passé comme sortie bash au dot command dot sh. Mais nous pouvons en fait écrire n'importe quel code personnalisé au-dessus de cela, et cela sera exécuté dans le cadre d'une tâche mais pas inclus dans le script bash.

Une des fonctions intégrées dans la syntaxe Nextflow s'appelle size. Donc je vais prendre l'entrée path, et je vais dire count underscore greetings, juste pour définir un nom de variable. Je vais prendre les input files et je vais appeler « size » dessus.

Cette fonction comptera la taille de ce channel d'entrée et l'attribuera à une variable.

Nous pouvons maintenant retourner cette variable dans le cadre du bloc output. Donc nous disons, val, car c'est une value, pas un fichier. Et count greetings.

Maintenant cela suffit en soi, et nous pourrions maintenant accéder à ces différentes sorties de ce process. Cependant, nous devrions y accéder de manière positionnelle. Donc en utilisant une clé d'index comme zéro et un.

Pour rendre un peu plus facile l'accès aux sorties, nous pouvons les nommer et nous faisons cela en utilisant une instruction emit.

Donc nous faisons, virgule emit out file ou quel que soit le nom que je veux donner à cela. Et je fais ici emit count. C'est essentiellement juste un décorateur, qui nous aide à écrire un code légèrement plus propre pour que nous puissions facilement référencer les sorties spécifiques plus tard dans le bloc workflow.

## 4.2. Rapporter la sortie à la fin du workflow

D'accord. Si je descends au bloc workflow, je peux maintenant prendre les sorties de collect greetings, faire collect greetings, dot out, et nous pouvons voir que nos deux sorties nommées sont suggérées ici par l'extension VS Code. Très pratique.

Donc je vais faire dot count pour obtenir la valeur count que nous venons de créer, et je vais faire view, pour qu'elle s'affiche dans la ligne de commande. Donc nous pouvons la voir lorsque nous exécutons le workflow.

Écrivons quelque chose dans la closure ici juste pour que ce soit un peu plus joli. num greetings, there were greetings greetings.

Et nous ne nous soucions pas vraiment de l'autre sortie car nous ne l'utilisons pas comme entrée pour d'autres processes. Mais vous pouvez voir comment nous pourrions facilement passer cela comme entrée à un autre process si nous le voulions, en aval.

## 4.3. Exécuter le workflow

Nous allons cliquer sur enregistrer. Jetons un coup d'œil au terminal et essayons-le.

D'accord, fantastique. Nous y voilà. There are three greetings. C'est exactement ça.

D'accord, excellent. C'est la fin de ce chapitre. Nous avons tout terminé pour être arrivés jusqu'ici. Vous commencez maintenant à construire un workflow assez réaliste, où nous sommes capables de gérer les entrées et les sorties et la logique dans notre workflow.

À mesure que ces fichiers de workflow deviennent plus longs, ils commencent à devenir un peu lourds. Donc dans le prochain chapitre, nous examinerons comment nous pouvons modulariser le code Nextflow dans des fichiers séparés pour qu'il soit plus facile de trouver et de maintenir le code dans le workflow.

Rejoignez-nous dans la prochaine vidéo pour le chapitre quatre. Hello Modules.

[Transcription de la vidéo suivante :octicons-arrow-right-24:](04_hello_modules.md)
