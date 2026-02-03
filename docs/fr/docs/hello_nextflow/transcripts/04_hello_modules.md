# Partie 4 : Hello Modules - Transcription

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Traduction assistée par IA - [en savoir plus et suggérer des améliorations](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

<div class="video-wrapper">
  <iframe width="560" height="315" src="https://www.youtube.com/embed/Xxp_menS0E8?si=0AWnXB7xqHAzJdJV&amp;list=PLPZ8WHdZGxmXiHf8B26oB_fTfoKQdhlik" title="YouTube video player" frameborder="0" allow="accelerometer; autoplay; clipboard-write; encrypted-media; gyroscope; picture-in-picture; web-share" referrerpolicy="strict-origin-when-cross-origin" allowfullscreen></iframe>
</div>

!!!note "Notes importantes"

    Cette page présente uniquement la transcription. Pour des instructions complètes étape par étape, retournez au [matériel de formation](../04_hello_modules.md).

    Les numéros de section affichés dans la transcription sont fournis à titre indicatif uniquement et peuvent ne pas inclure tous les numéros de section du matériel.

## Bienvenue

Bonjour, bienvenue dans la Partie Quatre de la formation Hello Nextflow.

Ce chapitre s'intitule Hello Modules, et nous allons parler de la façon de modulariser le code Nextflow. Ce que nous allons faire, c'est prendre notre script de workflow unique et le diviser en fichiers séparés.

Cela rend le code plus facile à naviguer et à maintenir au fur et à mesure que votre workflow grandit, et permet également de partager des modules entre pipelines, de sorte que si vous avez plusieurs pipelines utilisant le même outil, vous n'avez besoin d'écrire ce process qu'une seule fois.

Un exemple classique de ceci est le dépôt de modules nf-core, qui contient des milliers d'outils différents dans des modules prêts à l'emploi, que vous pouvez installer et utiliser dans votre workflow.

Nextflow peut également fonctionner avec des sub workflows, qui sont comme des modules, mais ils ont plusieurs processes. Cela dépasse le cadre de cette formation, mais cela fonctionne essentiellement de la même manière.

D'accord. Jetons un coup d'œil.

Comme d'habitude, commencez par aller sur training.nextflow.io.

Allez dans « Hello Nextflow » dans la barre latérale, et nous faisons la partie quatre : « Hello Modules ».

Je vais maintenant passer à mon environnement GitHub Code Spaces et examiner le fichier « hello-modules ».

Comme précédemment, nous commençons à partir du point final du chapitre précédent, donc ce script devrait vous sembler familier. Nous avons nos trois processes, say hello, convert to upper et collect greetings, et dans un workflow simple, qui exécute ces trois commandes et émet un message à la fin. Nous avons deux paramètres appelés greeting et batch, qui spécifie le nom utilisé pour le fichier de sortie collecté à la fin.

## 0. Échauffement : Exécuter hello-modules.nf

Nous pouvons vérifier que ce workflow fonctionne toujours comme prévu en faisant nextflow run hello, modules.

Parfait. Il a exécuté trois tâches avec chacun de ces processes, une tâche collectée, et il nous a dit qu'il y a trois salutations dans ce lot. Si nous allons dans results, nous avons nos différents fichiers de sortie ici, y compris le test collecté, la sortie batch.

## 1. Créer un répertoire pour stocker les modules

Bien. Faisons un peu de modularisation.

C'est généralement une bonne idée de mettre les modules dans un sous-répertoire de votre dépôt de pipeline, juste pour garder les choses organisées. Vous pouvez appeler cela comme vous voulez, mais par convention, nous l'appelons habituellement modules.

Alors allons-y, ouvrons un terminal et faisons make the modules. Vous pouvez le voir apparaître dans la barre latérale et VS Code ici.

## 2. Créer un module pour sayHello()

Je vais ensuite créer un nouveau fichier pour mon premier module. Vous pouvez faire « touch » ou « code » ou vous pouvez le faire dans la barre latérale, cela n'a vraiment pas d'importance. Donc je vais faire code modules et je vais le nommer d'après le process. Donc sayHello.nf . NF est une extension de fichier traditionnelle pour les fichiers Nextflow.

Je vais sauvegarder ici et nous pouvons voir notre nouveau fichier de module apparaître.

## 2.2. Déplacer le code du process sayHello vers le fichier module

Bien, ensuite je vais prendre le code du module depuis le workflow. Je vais aussi prendre le hash bang ici et le copier en premier pour qu'il soit clairement un fichier Nextflow. Et puis je vais prendre ce process et je vais le couper. Donc je vais le retirer de mon script de workflow principal et je vais le coller dans ce nouveau module.

C'est tout le contenu que ce fichier module va contenir. Juste un seul process, pas de workflow, pas de logique, juste un process seul.

Je peux maintenant fermer ce fichier.

## 2.3. Ajouter une déclaration d'import avant le bloc workflow

Maintenant mon workflow manque ce premier process, donc nous devons le ramener en l'important. La syntaxe pour cela est très similaire à d'autres langages de programmation, donc elle peut sembler familière. Nous faisons include accolades, le nom du process, say hello, puis from le chemin du fichier modules, say hello, nf. Fantastique.

Quelques astuces ici. L'extension VS Code est intelligente à ce sujet. Elle reconnaît ce chemin de fichier et vous pouvez passer la souris dessus et faire follow link. Ou je suis sur Mac, je peux faire option click et cela ouvre ce fichier. Donc nous pouvons rapidement y accéder.

Ce nom de process est maintenant utilisé par le workflow en bas ici, et nous pouvons faire la même chose ici. Il nous montre un peu d'information sur ce process, et encore une fois, je peux maintenir option, cliquer dessus, et cela l'ouvrira dans l'éditeur.

C'est donc un moyen très rapide lorsque vous avez beaucoup de fichiers pour vos différents processes de naviguer rapidement dans votre base de code dans VS Code.

D'accord. C'est essentiellement tout pour ce chapitre. Maintenant nous faisons juste la même chose à nouveau pour les autres processes.

## 3. Modulariser le process convertToUpper()

Donc créons un nouveau fichier ici. Appelons-le Convert to upper nf. Encore une fois, copions le hash bang. Et puis coupons le process.

Copions le nom du process là, include une nouvelle instruction include avec le nouveau nom de process.

## 4. Modulariser le process collectGreetings()

Et puis faisons la même chose pour le troisième process. Nouveau fichier, connect. Greetings,

faisons le hash bang. Coupons le process, collons le process, et faisons une nouvelle instruction include.

Maintenant vous pouvez voir ici j'ai un soulignement d'erreur ici disant invalid include source. Et c'est en fait une véritable erreur que j'ai faite parce que j'allais un peu trop vite. Si vous regardez attentivement, vous pouvez voir que j'ai manqué le T dans convert to upper

Donc VS Code m'a très utilement dit que j'avais fait une erreur là. Si je corrige ce nom de fichier, l'erreur disparaît. C'est un bon exemple de pourquoi la vérification des erreurs dans VS Code est si utile pour écrire du code Nextflow. Sinon je ne l'aurais pas remarqué et je ne l'aurais découvert que bien plus tard lorsque j'aurais essayé d'exécuter le workflow.

Notre script de pipeline principal est maintenant beaucoup plus simple. Il n'a aucun process dedans, nous avons juste trois instructions include et notre workflow. Nous n'avons changé aucune logique du workflow. Nous n'avons changé aucun code de process, donc j'espère qu'il devrait fonctionner exactement de la même manière.

## 4.4. Exécuter le workflow pour vérifier qu'il fait la même chose qu'avant

Vérifions. Je vais ouvrir un terminal et je vais exécuter exactement la même commande qu'avant.

Effectivement, il a exécuté nos processes, say hello, convert to upper collect greetings, et nous a donné trois salutations à nouveau.

Donc nous avons déplacé notre code, mais nous n'avons rien changé à la façon dont le workflow s'exécute et il est complètement inchangé. La seule différence est que nous avons maintenant un code plus propre, plus facile à maintenir et plus facile à partager avec les autres.

Et c'est tout. C'était un chapitre court. C'est un concept simple, mais il est très puissant et essentiel à la façon dont nous écrivons des workflows Nextflow plus complexes. Il est donc important que vous le compreniez et que vous preniez l'habitude de l'utiliser.

Dans le chapitre suivant, nous allons changer un peu de rythme et arrêter de penser autant à la syntaxe d'écriture du code Nextflow, et réfléchir un peu à la façon dont nous utilisons les logiciels dans les processes eux-mêmes. Rejoignez-nous dans la partie cinq pour Hello Containers.

[Transcription vidéo suivante :octicons-arrow-right-24:](05_hello_containers.md)
