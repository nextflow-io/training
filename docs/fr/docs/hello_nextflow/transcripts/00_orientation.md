# Orientation - Transcription vidéo

<div class="video-wrapper">
  <iframe width="560" height="315" src="https://www.youtube.com/embed/PIjOdFaYwWA?si=y8lAedhEHWaTV4zd&amp;list=PLPZ8WHdZGxmWKozQuzr27jyMGqp9kElVK&amp;cc_load_policy=1" title="YouTube video player" frameborder="0" allow="accelerometer; autoplay; clipboard-write; encrypted-media; gyroscope; picture-in-picture; web-share" referrerpolicy="strict-origin-when-cross-origin" allowfullscreen></iframe>
</div>

!!!note "Note importante"

    Cette page présente uniquement la transcription. Pour les instructions complètes étape par étape, retournez au [contenu du cours](../00_orientation.md).

## Bienvenue

Bonjour et bienvenue dans Hello Nextflow. Je m'appelle Phil Ewels. Je suis Product Manager pour les logiciels open source chez Seqera, l'entreprise derrière Nextflow.

Ce cours est une introduction pratique à la construction de workflows avec Nextflow. Il est conçu pour les personnes complètement nouvelles à Nextflow qui souhaitent développer leurs propres pipelines.

Les exemples sont tous basés sur du traitement de texte simple, afin que vous puissiez vous concentrer sur les concepts Nextflow sans avoir besoin d'expertise métier, juste une certaine familiarité avec la ligne de commande.

Nous allons parcourir les bases de Nextflow : écrire des processus, les connecter dans des workflows multi-étapes, gérer les dépendances logicielles avec des conteneurs, et configurer des pipelines pour différents environnements informatiques. À la fin, vous aurez construit un pipeline fonctionnel à partir de zéro.

Ce cours se concentre sur le _développement_ de pipelines. Si vous souhaitez simplement _exécuter_ des pipelines existants sans trop plonger dans le code, nous avons un cours plus court « Nextflow Run » qui pourrait mieux vous convenir.

Une fois que vous aurez acquis les bases ici, nous avons également des cours de suivi qui appliquent ces concepts à de véritables analyses scientifiques. Nous vous apprendrons comment utiliser les pipelines et les bonnes pratiques de la communauté nf-core.

Si vous êtes bloqué·e, rendez-vous sur community.seqera.io. Il y a un forum communautaire actif avec une section dédiée uniquement aux questions de formation. Vous pouvez l'utiliser à tout moment, cependant, nous organisons également des semaines de formation trimestrielles avec des personnes disponibles spécifiquement pour vous aider. Donc si vous suivez la formation pendant l'une de ces semaines, n'hésitez vraiment pas à demander de l'aide.

Vous pouvez également essayer de demander de l'aide à Seqera AI. C'est excellent pour expliquer le code Nextflow et vous aider au débogage.

Lorsque vous serez prêt·e à exécuter Nextflow à grande échelle, Seqera Platform est le meilleur endroit pour le faire. Il fonctionne sur votre infrastructure sans aucun verrouillage fournisseur, avec tout, du lancement de pipeline à la surveillance en temps réel, en passant par les environnements d'analyse interactifs. Mais pour l'instant, concentrons-nous simplement sur les fondamentaux.

Bien, commençons.

## training.nextflow.io

D'accord. La première chose à noter est que tous les cours de formation sur training.nextflow.io sont très interactifs. L'idée est que vous suiviez le matériel de formation et mes instructions, et nous parcourons le matériel de formation ensemble. Vous aurez donc besoin de deux choses : vous aurez besoin de votre ordinateur portable et vous aurez besoin de ce site web ouvert. Et c'est à peu près tout.

Voici donc la page d'accueil telle qu'elle apparaît aujourd'hui lorsque j'enregistre ceci. Vous pouvez voir qu'il y a un aperçu des différentes choses, le contexte, et les différents cours que nous avons, dont la liste s'allonge constamment.

Nextflow for newcomers est là où nous sommes. Il y a deux cours ici, Nextflow Run, qui est un cours différent, et Hello Nextflow, qui est celui qui nous intéresse.

Et vous pouvez également voir tous les différents cours dans la barre latérale. Je peux passer à Hello Nextflow, et nous pouvons voir tous les différents chapitres que nous allons parcourir ensemble.

Il y a quelques autres choses importantes à noter ici. Premièrement, le matériel de formation est versionné, donc vous pouvez voir ici en haut. Il indique 3.0 latest, qui au moment où j'enregistre est la dernière version stable. Cela changera avec le temps. Nous publions de nouveaux cours et nous mettons à jour le matériel au fil du temps. Donc si c'est 3.1 ou 3.2, ne vous inquiétez pas trop. Si c'est 4.0, alors il y a probablement une nouvelle vidéo, et vous devriez peut-être aller la trouver car il y aura probablement des mises à jour importantes.

Un autre menu déroulant en haut est celui-ci, pour la langue. C'est tout nouveau pour la version 3.0. Nous avons pris le matériel précédemment traduit, qui avait été fait par des humains, à la main, et nous l'avons passé dans un LLM et mis en place toute cette nouvelle infrastructure pour maintenir différentes traductions du matériel de formation en utilisant la traduction LLM.

Nous avons donc maintenant toutes ces fantastiques traductions ici. Donc si vous voulez suivre en coréen, vous pouvez charger tout le site web en coréen. Et suivre là-bas. Pareil pour toutes ces autres langues, hindi et allemand et ainsi de suite. Je vais suivre en anglais. C'est la langue principale dans laquelle nous écrivons le matériel.

Quelques autres boutons si vous aimez avoir le mode clair. Au lieu de ce mode, vous pouvez suivre le site web en mode clair en haut ici.

Et puis aussi tout ce que nous regardons est dans un seul dépôt GitHub, qui est open source, appelé nextflow-io/training. Et si vous cliquez sur ce bouton à tout moment, cela vous amènera au dépôt GitHub. Nous y reviendrons dans une minute.

## Configuration de GitHub Codespaces

D'accord, maintenant que vous avez ceci ouvert dans l'onglet du navigateur. Allons sur Hello Nextflow et cliquons. Vous pouvez voir sur la page d'introduction, elle nous indique quelques-unes des exigences, l'aperçu, et le plan de cours de ce que nous allons couvrir approximativement, puis nous allons plonger dans la mise en route.

Il existe différentes façons de faire ce tutoriel interactif. Si vous êtes à l'aise, vous êtes bienvenu·e à faire cela localement sur votre propre ordinateur avec votre propre installation de Nextflow. Si nous cliquons sur Options d'environnement, vous pouvez voir qu'il y a plus de détails sur la façon de faire cela soit en utilisant des Devcontainers locaux, soit vous pouvez également simplement installer tous les logiciels localement, avec une installation manuelle.

Nous travaillons à faire fonctionner cela correctement avec Seqera Studios, donc c'est une autre option. Mais la plus courante en ce moment est d'utiliser GitHub Codespaces.

Codespaces configure un environnement sandbox sur un serveur distant géré par GitHub. Et c'est gratuit pour une certaine quantité d'utilisation, ce qui est généralement suffisant pour la formation. Et il vous configurera une instance VS Code, un IDE où vous pouvez accéder à tous les fichiers du dépôt, exécuter Nextflow et tout. Et nous avons préconfiguré Codespaces pour vous. Il a donc tout ce dont vous avez besoin.

La beauté de ceci est qu'il suffit d'un clic pour configurer un Codespace. C'est pareil pour tout le monde, et nous savons que vous avez déjà tous les prérequis installés, donc c'est rapide et efficace.

Donc la première chose à faire est d'aller sur « Getting Started ». Cherchez ce bouton, qui dit _Open in Codespaces_. Je vais faire cmd \+ clic dessus pour l'ouvrir dans un nouvel onglet, et cela nous amène à GitHub.

Voici à quoi cela ressemble. Nous pouvons voir que nous avons défini toutes les options ici pour vous. Si vous le souhaitez, vous pouvez cliquer sur change options. Quelques choses que vous pouvez faire ici. Vous pouvez donner une instance machine plus grande, par exemple, si vous constatez qu'elle plante parce qu'elle manque de mémoire ou quoi que ce soit de ce genre. Ou définir des versions spécifiques du matériel de formation. Mais généralement, vous pouvez simplement suivre ce que nous avons configuré ici et vous pouvez le voir. Dans ce cas, il utilise la version 3.0.

Je vais donc cliquer sur create new Codespace. Et cela m'amène dedans.

Notez également qu'il indique no Codespace to resume là. Si j'ai précédemment créé un Codespace, cliquer à nouveau sur ce bouton sur le matériel de formation m'amènera à la même page et listera tous les Codespaces que j'ai déjà en cours d'exécution. Ensuite, vous pouvez simplement y retourner directement et continuer là où vous vous étiez arrêté·e. Donc peu importe si vous avez fermé votre ordinateur portable.

Ils s'arrêtent automatiquement après quelques minutes d'inactivité, mais ce n'est pas un problème. Vous pouvez simplement les redémarrer.

Une fois que vous démarrez un nouveau Codespace, il va rester sur cette page comme ceci et il va charger pendant un bon moment. C'est donc le bon moment pour faire une petite pause. Peut-être que vous avez oublié d'aller aux toilettes ou vous voulez une tasse de thé avant de commencer ? Allez-y maintenant pendant que vous attendez, car cela va tourner là pendant un moment.

Rapidement pendant que nous attendons qu'il charge, je vais également aller sur github.com/codespaces et juste montrer que c'est la page d'aperçu où vous pouvez voir tous les différents Codespaces que vous avez actuellement en cours d'exécution.

Vous pouvez voir que j'en ai un ici pour nextflow-io/training. Aucun changement, car je n'ai encore rien fait dedans. La quantité de ressources qu'il utilise, et vous pouvez voir qu'en ce moment il se configure. Je peux aller ici, cliquer sur ce petit menu déroulant et cliquer sur delete. Donc si vous configurez accidentellement plusieurs Codespaces et que vous n'en utilisez pas certains, vous pouvez supprimer les anciens et nettoyer.

Enfin, une autre façon d'accéder à ceci. Si nous allons au dépôt GitHub. Et cela fonctionne pour n'importe quel dépôt GitHub. Cliquez sur code. Vous pouvez avoir des commandes pour cloner le dépôt localement. Et il y a un onglet appelé Codespaces. Et encore une fois, vous pouvez en créer un nouveau, et vous pouvez voir ceux qui sont déjà en cours d'exécution.

Donc encore une fois, si vous oubliez comment vous avez créé votre Codespace, vous pouvez toujours y revenir de cette façon.

## L'interface VS Code

D'accord, la construction est terminée et elle commence maintenant à charger GitHub Codespaces. Cela ne prend pas toujours aussi longtemps, donc ne vous inquiétez pas. C'est juste la première fois que vous créez le Codespace. Si vous retournez dans un qui existe déjà, c'est beaucoup plus rapide.

Ne soyez pas trop impatient·e si c'est la première fois, ce n'est pas encore terminé, même si cela commence à nous donner une interface.

Mais pendant que nous attendons que les dernières choses soient configurées, je vais juste vous faire découvrir l'interface au cas où vous ne seriez pas très familier·ère avec VS Code.

Premièrement, il y a la barre latérale de chat pour les trucs IA, dont nous n'avons pas besoin. Je vais donc fermer cela, m'en débarrasser et libérer de l'espace.

Sur la gauche, nous avons l'explorateur de fichiers qui nous montre tous les fichiers dans le dépôt Git, qui est l'espace de travail que nous avons créé. Notez que ce ne sont pas des fichiers locaux. Tout cela est sur le serveur distant où nous travaillons. Vous pouvez glisser-déposer des fichiers locaux et des choses, mais pour la plupart, nous n'allons pas penser à cela aujourd'hui. Nous travaillons simplement purement à distance.

Il y a d'autres outils dans cette barre latérale, par exemple, la recherche. Vous pouvez donc rechercher tous les fichiers dans un dépôt en une seule fois. Et si nous faisions du travail de développement sur le dépôt de formation, nous pourrions faire l'intégration avec le contrôle de source avec Git et le débogage et d'autres choses.

D'autres choses sont, il y a une fenêtre principale d'édition de code ici, qui vient de charger un aperçu du readme, qui est pour le matériel de formation. Donc dans ce cas, il visualise du markdown, mais normalement ce sera un éditeur de code.

Et puis en dessous, nous avons le terminal, qui est où nous allons exécuter toutes nos commandes et interagir directement avec Nextflow.

Tout dans le Codespace est préinstallé, donc la commande Nextflow est déjà là et ainsi de suite.

D'accord. Quand vous arrivez jusque-là, cela devrait être à peu près terminé. Vous pouvez voir maintenant qu'il a téléchargé le serveur de langage Nextflow et il a configuré quelques extensions pour nous dans VS code, y compris l'extension Nextflow, qui va être utile. Je peux donc fermer cela et je peux fermer le README.md.

Et maintenant vous pouvez voir que j'ai un peu plus sur le côté gauche. Je suis un peu zoomé ici, mais si je dézoome vous pouvez voir qu'un des boutons dit Nextflow avec l'icône Nextflow. Et cela a des trucs sympas ici pour explorer le projet et des choses, auxquelles nous reviendrons plus tard.

D'accord. Au cas où vous perdriez l'un de ces panneaux, ces boutons en haut à droite sont vraiment utiles et ils montrent et cachent simplement les choses. Donc cela montre et cache l'Explorateur, montre et cache le terminal en bas. Et ainsi de suite.

Je vais utiliser ceux-ci assez souvent parce que je suis très zoomé, donc j'essaie de vous aider à voir tout le texte sur mon écran, et donc c'est utile de pouvoir passer en plein écran avec le terminal puis le cacher quand nous regardons le code. Mais la plupart du temps, vous pouvez simplement avoir tout ce truc ouvert en même temps.

D'accord, quoi d'autre à regarder ? Pas grand-chose de plus. Notez que Nextflow, comme je le dis, est installé. Je peux donc taper « nextflow -version » et cela devrait apparaître en disant quelle version nous avons installée.

Il y a d'autres trucs installés ici aussi. À la fin de chaque chapitre, nous avons un ensemble de questions de quiz, par exemple, sur le site web. Et vous pouvez également les faire dans le terminal si vous le souhaitez en tapant quiz.

Il y a d'autres raccourcis clavier que je vais utiliser, juste au cas où vous seriez curieux·euse. Par exemple, juste là j'ai appuyé sur cmd\+K sur mon Mac et cela a effacé le terminal, pour se débarrasser de toute la sortie précédente. C'est donc agréable pour garder les choses propres. Si vous me voyez faire, c'est comme ça que je le fais.

Et aussi si vous êtes nouveau·elle au terminal, rappelez-vous que vous pouvez utiliser tab pour auto-compléter, ce que je vais faire beaucoup pour auto-compléter les chemins.

Je peux donc voir sur la gauche ici qu'il y a un dossier appelé Hello Nextflow, qui est ce que nous allons parcourir. Si je fais « ls » pour lister les fichiers, je peux faire « hel », appuyer sur tab, auto-complète. Et donc c'est un moyen très rapide de compléter les chemins.

## Ouvrir uniquement le dossier Hello Nextflow

D'accord. C'est génial. Il y a beaucoup de choses dans ce dépôt cependant.

Il y a tous les fichiers pour générer le site web, et il y a plusieurs cours différents ici, et vous pouvez le faire depuis cette racine et simplement cliquer dans le dossier « Hello Nextflow ». Mais c'est agréable de se concentrer purement sur ceci.

Vous pouvez définir ceci comme votre espace de travail avec un tas de clics ici et là et définir un répertoire de projet et des trucs. Mais le moyen le plus simple est de taper code, qui est la commande CLI pour lancer VS Code, puis « hello-nextflow ».

Cela ouvrira un nouvel onglet de navigateur et vous pouvez fermer l'ancien. Et cela ressemble exactement à la même chose. Mais maintenant vous pouvez voir que nous sommes dans ce sous-répertoire et tous les autres fichiers sont invisibles, et nous avons une configuration plus propre.

Vous pouvez voir ici qu'également le répertoire de travail actuel est maintenant dans le dossier Hello Nextflow. Donc propre et net. Nous n'avons pas besoin de nous inquiéter d'être au mauvais endroit. D'accord.

## Nouvelle syntaxe Nextflow pour 2026

Il y a une chose spéciale que je dois mentionner à ce stade. En ce moment, au début de 2026, nous commençons à introduire différentes fonctionnalités dans Nextflow, et l'une des grandes nouvelles est un nouveau parseur de syntaxe de langage à l'intérieur de Nextflow.

Fondamentalement, le moteur qui lit vos fichiers Nextflow et comprend cela, pour l'exécution. Il y a quelques changements dans la syntaxe, et il est vraiment important que vous utilisiez Nextflow avec le parseur de syntaxe correct activé.

Vous avez besoin de deux choses pour cela. Vous avez besoin d'une version à jour de Nextflow et vous devez vous assurer qu'elle est activée.

Si je fais « nextflow -version » à nouveau, vous verrez que le Codespaces fonctionne avec 25.10.2 et 25.10 est la version minimale pour pouvoir utiliser ce truc.

Si vous utilisez 26.04, qui pour moi n'est pas encore sorti, mais le sera bientôt. Alors cela exécutera le nouveau parseur de syntaxe par défaut, et vous n'avez rien d'autre à faire.

Mais si vous exécutez 25.10, vous devez activer le parseur de syntaxe strict, comme on l'appelle, ou parseur de syntaxe v2.

Cela se fait avec une variable d'environnement. Elle est déjà définie dans les Codespaces, donc vous n'avez rien à faire. Mais si vous exécutez localement, vous devez définir cela, et je peux vérifier cela en faisant « echo $NXF_SYNTAX_PARSER », et cela devrait être défini sur v2.

Donc si vous exécutez localement, faites simplement « export NXF_SYNTAX_PARSER=v2 ». Simple comme ça. Mais n'oubliez pas de faire cela. Sinon vous allez voir des divergences bizarres et des erreurs au fur et à mesure.

Si vous avez le moindre doute sur tout ce truc autour de la version Nextflow et du parseur de syntaxe, premièrement, rappelez-vous, vous n'avez pas besoin de vous inquiéter si vous êtes dans Codespaces. Tout devrait être configuré correctement. Mais deuxièmement, si vous allez sur le matériel de formation Nextflow, si vous descendez, parlez des exigences de version, il y a un lien ici qui vous amène à la page d'aide autour des versions d'exploration, et cela passe en revue tout cela en détail.

Cela vaut la peine de lire ceci si vous avez un moment. Car cela aide à clarifier ce que sont certains des différents termes, que vous pourriez entendre lorsque vous commencez à utiliser Nextflow. Des choses comme DSL1, DSL2, parseur de syntaxe un, parseur de syntaxe deux, et ainsi de suite. Donc cela vaut la peine de juste vérifier cela et cela répète une partie de ce que je viens de dire.

C'est également vraiment utile si vous avez précédemment écrit du code Nextflow et que vous revenez pour un rappel. Cela vous indique certaines des choses qui changent et vous renvoie vers des parties de la documentation Nextflow, qui vous indiquent comment mettre à jour votre code Nextflow.

## Fichiers du cours

D'accord. Dernière chose pour nous familiariser, c'est juste voir les fichiers, qui sont dans ce répertoire. Vous pouvez soit regarder dans la barre latérale, soit souvent dans le matériel de formation, nous utilisons la commande tree, -L, qui est le nombre de niveaux à regarder. Nous dirons deux, et si je mets ceci en plein écran, vous verrez que cela reflète exactement ce que nous voyons dans la barre latérale là, mais cela exclut les fichiers cachés, qui commencent par un point.

Donc les fichiers \*.nf, cela signifie Nextflow. Ce sont donc les fichiers de script Nextflow, et il y a un fichier de démarrage ici pour chacun des différents chapitres du matériel de formation, que nous ouvrirons et explorerons puis éditerons.

Nous changerons ces fichiers au fur et à mesure, et donc à la fin de chaque chapitre, les fichiers devraient ressembler à peu près au début du chapitre pour le suivant. Mais nous vous donnons ces différents fichiers afin que vous puissiez toujours recommencer à zéro et ne pas trop vous inquiéter de gâcher la syntaxe.

Si vous avez besoin de comparer à quelque chose qui devrait définitivement fonctionner. Vous pouvez vérifier dans le dossier solutions, et c'est comme un état final pour chacun des chapitres, donc vous pouvez comparer ce que vous avez écrit avec ce qui est là.

Il y a un répertoire data. Cela a juste un fichier greetings.csv, que nous utiliserons comme exemple, données d'entrée dans une partie du cours, et des choses comme un fichier de configuration et quelques paramètres, que nous décrirons plus tard dans le cours.

## Conclusion

D'accord. Donc maintenant, j'espère que tout fonctionne. Votre écran ressemble au mien et vous comprenez comment accéder à tout et ce que sont tous les différents fichiers.

Si vous faites défiler vers le bas de la page sur la mise en route, petite case à cocher vous devriez dire que je comprends ce que je fais. Mon environnement est opérationnel et vous êtes prêt·e, votre répertoire de travail est correctement défini sur le dossier « Hello Nextflow ».

Si vous avez coché tout cela et qu'ils sont verts. Nous pouvons continuer vers la vidéo suivante et le chapitre suivant, qui est la partie un. Hello World. À tout de suite.
