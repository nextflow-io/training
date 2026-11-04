# Partie 4 : Hello Modules - Transcription vidéo

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Traduction assistée par IA - [en savoir plus et suggérer des améliorations](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

<div class="video-wrapper">
  <iframe width="560" height="315" src="https://www.youtube.com/embed/43Ot-f0iOME?si=0AWnXB7xqHAzJdJV&amp;list=PLPZ8WHdZGxmWKozQuzr27jyMGqp9kElVK&amp;cc_load_policy=1" title="YouTube video player" frameborder="0" allow="accelerometer; autoplay; clipboard-write; encrypted-media; gyroscope; picture-in-picture; web-share" referrerpolicy="strict-origin-when-cross-origin" allowfullscreen></iframe>
</div>

!!!note "Notes importantes"

    Cette page présente uniquement la transcription. Pour les instructions complètes étape par étape, retournez au [matériel de cours](../04_hello_modules.md).

    Les numéros de section indiqués dans la transcription sont fournis à titre indicatif uniquement et peuvent ne pas inclure tous les numéros de section du matériel.

## Bienvenue

Bonjour et bienvenue dans la partie quatre de Hello Nextflow. Cette section est entièrement consacrée aux modules, et c'est une section assez courte du cours. Nous n'allons pas vraiment écrire beaucoup de code, il s'agit plutôt de la façon dont nous organisons le code dans notre pipeline.

Jusqu'à présent, nous avons tout mis dans un seul fichier, ce qui est bien, et c'est d'ailleurs comme ça que nous construisions les pipelines Nextflow autrefois.

Mais à mesure que ce pipeline grandit, le script devient de plus en plus long et de plus en plus difficile à parcourir, à maintenir, et cela signifie aussi que nous ne pouvons pas vraiment partager le code.

Les modules Nextflow nous permettent d'extraire les processus de ce script principal et ensuite de les importer. Cela signifie que le code est plus facile à parcourir et cela signifie aussi que nous pouvons partager ce code de module entre différents pipelines.

Ce petit diagramme sur la page principale de la documentation illustre bien le concept. Au lieu d'un énorme script, nous allons inclure ces fichiers de modules séparés, provenant de différents scripts de modules, et tout va être intégré dans le workflow, mais ça va toujours s'exécuter exactement de la même manière.

Alors passons à GitHub Codespaces et jetons un petit coup d'œil. Comme précédemment, j'ai un peu nettoyé mon espace de travail ici. J'ai supprimé les anciens répertoires Nextflow et le répertoire work et ainsi de suite. Mais ce n'est pas grave si vous avez encore ces fichiers.

Je vais commencer à travailler dans le fichier hello modules, qui est essentiellement là où nous l'avons laissé à la fin du chapitre précédent. Nous avons nos trois processus ici. Nous avons quelques params, le bloc workflow, où nous exécutons ces trois processus et les assemblons avec des canaux. Ensuite, nous publions les canaux de sortie et nous avons le bloc output indiquant comment publier ces fichiers.

## 1. Créer un répertoire pour stocker les modules

Maintenant, comme je l'ai dit, nous n'allons pas vraiment écrire ou modifier beaucoup de code. Nous allons simplement déplacer le code que nous avons déjà. Les fichiers de modules Nextflow contiennent généralement un seul processus, et par convention nous les conservons normalement dans un répertoire appelé modules. Mais vous pouvez l'appeler comme vous voulez. Mais je vais garder un répertoire modules dans mon dépôt ici, et ensuite je vais créer un fichier pour chaque processus. Donc je vais dire nouveau fichier, sayHello.nf.

## 2. Créer un module pour sayHello()

Maintenant je vais prendre mon processus et je vais simplement sélectionner ce code, le couper du fichier principal hello modules et le coller ici.

Évidemment, cela ne fait rien par lui-même. Notre script principal a toujours besoin de ce processus, donc nous devons le réintégrer d'une manière ou d'une autre. Et nous faisons cela avec l'instruction include.

Donc je tape include et des accolades, et ensuite je prends le nom du processus. Et je dis from, et ensuite je donne un chemin de fichier relatif. Donc ça dit, commence par ./ parce que c'est relatif à l'endroit où ce script est enregistré. Donc c'est modules sayHello.nf.

Remarquez que l'extension VS Code est assez utile ici. Elle nous indique si elle peut trouver ce fichier et si elle peut trouver un processus, que je nomme. Si je mets délibérément une faute de frappe ici, elle me donne une erreur tout de suite et elle me dira qu'elle ne peut pas trouver ce processus que j'essaie d'importer. Donc gardez un œil sur les erreurs que vous trouvez.

Et c'est vraiment tout. Nous avons toujours notre processus ici. Il n'y a pas de modifications nécessaires en bas ici. Le processus a le même nom et il est exécuté exactement de la même manière. C'est juste que le code réel du processus est maintenant dans un fichier séparé.

Nous pouvons exécuter à nouveau le workflow Nextflow, il va fonctionner exactement de la même manière. Et c'est essentiellement le reste de ce chapitre du cours qui consiste simplement à déplacer ces trois processus dans leurs propres fichiers.

Alors faisons-le maintenant. Je vais rapidement créer un nouveau fichier de module pour le deuxième processus : convertToUpper.nf. Je vais couper ce code, le coller là. Et ensuite je vais inclure celui-là. Voyons, super.

Et ensuite je vais créer un nouveau fichier pour collectGreetings.nf. Couper ça.

Beaucoup de couper, couper et copier et coller.

Et maintenant notre script de workflow principal semble soudainement beaucoup, beaucoup plus court, beaucoup plus abordable et beaucoup plus facile à lire.

Et vous pouvez voir comment le projet commence maintenant à se construire avec nos différents fichiers. Nous pouvons plonger dans les détails aux endroits que nous voulons. Naviguer pour trouver des étapes spécifiques dans le pipeline beaucoup plus facilement, et obtenir un aperçu de ce que fait le pipeline rapidement.

## Naviguer dans les modules avec VS Code

Maintenant, bien sûr, l'inconvénient de faire cela est que si vous avez un gros pipeline, vous aurez beaucoup de fichiers de modules et ils pourraient être organisés dans plusieurs sous-répertoires ou toutes sortes de choses. Maintenant, encore une fois, un petit conseil ici. L'extension VS Code est plutôt bonne pour naviguer dans votre base de code et aussi pour vous informer sur le code là-bas.

Vous pouvez voir que VS Code comprend ce qu'est ce processus et me donne un petit aperçu lorsque je survole, donc je peux voir sans avoir à aller chercher le code source, quelles sont les entrées et les sorties, ce qui est généralement la chose la plus importante quand je l'utilise dans un workflow.

Et aussi si je maintiens command, je suis sur Mac, et je clique sur le nom du processus, il ouvre le fichier directement tout de suite. Le fait apparaître. Donc je peux simplement y aller directement sans même penser à ce que sont les chemins de fichiers réels. Et ça fonctionne partout, je peux faire ça aussi, là où les processus sont appelés. Donc ça rend vraiment rapide.

## 4.4. Exécuter le workflow

D'accord, vérifions simplement que le pipeline s'exécute toujours comme nous l'attendons. Donc j'ouvre le terminal. Faisons "nextflow run hello modules", et voyons s'il s'exécute sans problèmes.

J'espère que tout l'intérêt de ceci est que le pipeline est essentiellement inchangé, donc vous ne devriez pas vraiment voir de changements par rapport à quand nous l'avons exécuté avant. La sortie ici semble exactement la même, et vous pouvez voir notre répertoire results avec tous les mêmes fichiers, donc c'est super. Pas de changement, c'est bien.

## Une note sur nf-core/modules

Juste avant de conclure, je veux rapidement aborder la puissance de la collaboration en ce qui concerne les modules. Ces fichiers sont dans mon dépôt, donc il n'est pas évident tout de suite comment nous pourrions collaborer dessus. Et il existe de nombreuses façons différentes de le faire, mais probablement l'exemple le plus grand et le plus connu de cela est nf-core.

Si je vais sur le site web nf-core, je vais dans ressources, et modules. Vous pouvez voir que nf-core a une énorme bibliothèque de modules, près de juste en dessous de 1700 modules quand je regarde ceci. Et donc je peux taper le nom de n'importe lequel de mes outils préférés, aller voir si quelqu'un d'autre a déjà écrit un module pour ça, et voir ce processus de module pré-écrit ici avec toutes les entrées, les sorties, les conteneurs logiciels, toutes ces informations, et vous pouvez voir sur le côté ici, combien de pipelines nf-core différents utilisent tous ce processus partagé unique.

C'est un exemple un peu extrême, mais vous pouvez voir que c'est vraiment réutiliser ce code. Et si je clique pour aller à la source GitHub de ceci, c'est exactement la même chose que ce que nous faisons. C'est juste un gros processus dans un fichier.

Maintenant du côté nf-core, nous faisons quelques astuces pour pouvoir partager ces fichiers et les intégrer dans différents dépôts. Et si vous voulez en savoir plus à ce sujet, allez consulter le cours que nous avons sur l'utilisation et la construction avec nf-core spécifiquement. Mais je voulais juste vous donner une idée de à quel point ce concept de réutilisation de code peut être puissant.

## Conclusion

Bon, c'est tout pour les modules. Je vous avais dit que c'était une section courte du cours. Consultez le quiz, assurez-vous de bien comprendre et assurez-vous que tout fonctionne toujours correctement. Et je vous retrouve dans la prochaine vidéo, qui est entièrement consacrée aux conteneurs logiciels. Merci beaucoup.

I.
