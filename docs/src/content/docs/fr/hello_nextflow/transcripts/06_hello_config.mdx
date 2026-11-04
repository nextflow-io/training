# Partie 6 : Hello Config - Transcription de la vidéo

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Traduction assistée par IA - [en savoir plus et suggérer des améliorations](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

<div class="video-wrapper">
  <iframe width="560" height="315" src="https://www.youtube.com/embed/FcZTi6mRkITY0zW_&amp;list=PLPZ8WHdZGxmWKozQuzr27jyMGqp9kElVK&amp;cc_load_policy=1" title="YouTube video player" frameborder="0" allow="accelerometer; autoplay; clipboard-write; encrypted-media; gyroscope; picture-in-picture; web-share" referrerpolicy="strict-origin-when-cross-origin" allowfullscreen></iframe>
</div>

!!!note "Notes importantes"

    Cette page présente uniquement la transcription. Pour les instructions complètes étape par étape, consultez le [matériel de cours](../06_hello_config.md).

    Les numéros de section indiqués dans la transcription sont fournis à titre indicatif uniquement et peuvent ne pas inclure tous les numéros de section du matériel.

## Bienvenue

Bonjour, et bienvenue dans la Partie Six de Hello Nextflow. Cette section est entièrement consacrée aux configurations, et c'est la dernière partie de ce cours.

Nextflow excelle particulièrement dans deux domaines : la reproductibilité et la portabilité. C'est avec les configurations que nous voyons vraiment briller le second de ces atouts. La capacité de configurer un pipeline Nextflow pour qu'il s'exécute de différentes manières et fonctionne sur différents systèmes, sans avoir à modifier le code sous-jacent du pipeline.

Ce super pouvoir permet aux pipelines Nextflow d'être réutilisés par d'autres personnes dans différents endroits, ou sur différentes infrastructures auxquelles vous pourriez avoir accès vous-même.

Cela signifie que vous pouvez développer le code du pipeline sur votre ordinateur portable, le déployer sur le cloud, l'exécuter sur votre HPC, et c'est le même code de pipeline et il s'exécute partout.

Dans cette section, nous allons aborder quelques sujets. Nous commencerons par la manière dont Nextflow gère les fichiers de configuration, d'où il les charge, comment vous les écrivez et comment vous les structurez, et cette séparation entre le pipeline lui-même et ce qui devrait aller dans un fichier de configuration.

Ensuite, nous passerons à certains cas d'usage courants tels que la modification de l'emplacement de stockage des fichiers de sortie, et également comment faire fonctionner le pipeline sur différentes infrastructures, en utilisant différents types d'empaquetage logiciel ou en soumettant des tâches à différentes infrastructures.

## Hiérarchies de fichiers de configuration

Très bien, commençons. En ce qui concerne le chargement des fichiers de configuration, Nextflow peut les extraire de nombreux endroits différents, ce qui est une bonne chose et peut aussi être un peu risqué car parfois il peut être un peu difficile de savoir d'où il obtient un fichier de configuration et dans quel ordre il charge les choses.

Je vous recommande donc vraiment de cliquer sur ce lien ici, qui nous amène à la documentation Nextflow. Et sur cette page de configuration, elle liste les emplacements clés d'où la configuration est chargée, et surtout, l'ordre dans lequel ces éléments sont chargés.

Vous pouvez voir que vous pouvez placer un fichier de configuration dans votre répertoire home Nextflow, qui est généralement ".nextflow" dans votre répertoire home. Et ce fichier sera toujours chargé par chaque exécution de Nextflow sur votre système.

Le prochain endroit où chercher est un fichier à la racine de votre dépôt ou répertoire de pipeline appelé "nextflow.config".

Ensuite, un autre fichier appelé "nextflow.config", mais cette fois dans le répertoire à partir duquel vous lancez Nextflow : le répertoire de lancement.

Enfin, vous pouvez fournir des chemins de fichiers de configuration sur la ligne de commande avec un argument "-c", et vous pouvez le faire plusieurs fois. Et ils sont appliqués dans l'ordre que vous spécifiez.

Vous pouvez fournir des fichiers de configuration dans tous ces emplacements si vous le souhaitez, et ils seront chargés de manière itérative, chacun écrasant le précédent uniquement dans les scopes de configuration où ils entrent en conflit.

C'est un système vraiment puissant car il signifie que vous pouvez définir des valeurs par défaut sensées et ensuite devenir de plus en plus spécifique à mesure que vous vous concentrez sur cette configuration.

## 0. Échauffement : exécuter hello-config.nf

D'accord, fermons cela et passons à notre Codespaces pour commencer. Comme précédemment, j'ai fait le ménage ici, j'ai supprimé mes répertoires de résultats précédents, mes répertoires Nextflow et work, etc. Ne vous inquiétez pas si vous avez encore ces fichiers qui traînent. C'est juste parce que j'ai beaucoup zoomé et donc les choses deviennent vite désordonnées autrement.

Nous allons travailler avec hello-config.nf, le dernier fichier de notre répertoire, et cela devrait faire suite à ce que nous avons fait dans la section précédente.

Nous avons donc nos quatre processus différents, qui sont inclus à partir de fichiers de modules. Nous avons nos paramètres de pipeline, notre bloc workflow où nous appelons les différents processus et assemblons les canaux, publions les canaux de sortie, puis le bloc output en bas où nous définissons où ces fichiers doivent être stockés et comment ils doivent être copiés.

Nous avons également déjà un fichier "nextflow.config" du dernier chapitre, où nous activons Docker, et nous allons développer ce fichier aujourd'hui.

Comme précédemment, nous avons modifié le chemin de sortie dans ce script principal pour hello config, juste pour qu'il n'entre pas en conflit avec les résultats précédents que vous avez générés.

D'accord, vérifions simplement rapidement que tout fonctionne toujours comme prévu. Ouvrons un terminal et faisons nextflow run hello-config.nf. Nextflow se charge. Devrait exécuter nos quatre processus différents. Générer du bel art ASCII en utilisant cowpy puis sauvegarder nos résultats dans nos fichiers de résultats dans ce répertoire.

Je peux jeter un coup d'œil rapide ici juste pour m'assurer que ces fichiers ressemblent à ce que nous attendons, et effectivement, voici notre dinde géante. Parfait.

## 1.1. Déplacer les valeurs par défaut vers nextflow.config

Maintenant, la première chose que nous allons faire est de commencer à déplacer certaines choses de notre script vers notre fichier de configuration.

Et ce qui nous intéresse, ce sont principalement les paramètres à ce stade. Nous voulons prendre les valeurs par défaut et les mettre dans le fichier de configuration, pour qu'il soit plus clair quelles sont les valeurs par défaut et qu'il soit plus facile pour les gens de les écraser.

Je vais juste prendre ce bloc params ici du script et le mettre dans le fichier de configuration. Et nous devons être un peu prudent·e·s ici, car actuellement la syntaxe est légèrement différente entre la configuration et les scripts. Le fichier de configuration ne peut pas prendre de déclarations de types car nous ne définissons pas vraiment ces paramètres, nous les référençons simplement. Je vais donc me débarrasser de ceux-ci.

Mais sinon c'est pratiquement la même chose. Nous avons un bloc params puis nous avons nos différents paramètres d'entrée, paramètre batch, paramètre character.

Je peux maintenant retourner à mon script et je n'ai plus besoin de définir ces valeurs par défaut car ces valeurs sont maintenant dans mon fichier Nextflow config.

Cependant, je laisse les noms de paramètres et leurs types, pour que Nextflow connaisse cette information et puisse toujours faire toute la sécurité de type et tout.

D'accord. Nous sauvegardons ces fichiers et vérifions rapidement que tout fonctionne toujours de la même manière qu'avant. Il ne devrait y avoir aucun changement ici. Nous avons conservé les mêmes valeurs. Nous avons juste déplacé l'endroit où elles ont été définies.

Parfait.

## 1.2. Utiliser un fichier de configuration spécifique à l'exécution

Maintenant, jusqu'à présent, nous avons lancé Nextflow depuis le même répertoire où nous avons notre script de pipeline. Donc notre répertoire de lancement et notre répertoire de pipeline sont en quelque sorte la même chose.

Pour montrer comment nous pouvons avoir différents fichiers de configuration avec différents répertoires de lancement, nous allons créer un nouveau sous-répertoire maintenant.

Je vais donc faire mkdir, et nous allons l'appeler tux-run.

Ensuite je vais faire cd, changer de répertoire vers tux-run. Et notez que notre répertoire de travail n'est maintenant plus dans le même répertoire que les scripts de pipeline.

D'accord, créons un nouveau fichier "nextflow.config". Donc touch Nextflow config, et ouvrons-le simplement dans VS Code. Vous pouvez voir aussi dans la barre latérale ici que nous sommes maintenant dans ce sous-répertoire.

Maintenant nous pouvons prendre le même bloc params que nous avions dans le nextflow.config de niveau supérieur, copier ceci et maintenant nous pouvons changer ces valeurs.

Premièrement, les données sont maintenant un chemin relatif différent car nous sommes dans un sous-répertoire, nous devons donc mettre à jour cela. Ensuite nous allons changer batch en experiment, et nous allons changer le character de Turkey à tux.

Maintenant cliquons sur sauvegarder là, et essayons. Tout comme avec data, je dois maintenant dire ../ pour arriver au script. Donc c'est Hello config. Et j'appuie sur entrée.

Le code du pipeline n'a pas changé du tout, mais maintenant nous allons avoir deux ensembles de configurations qui se chargent, et le fichier de configuration du répertoire de lancement devrait écraser les valeurs par défaut, qui ont été définies dans le nextflow.config du pipeline, et nous devrions obtenir différents ensembles de résultats.

Effectivement, dans notre répertoire ici, dans tux-run, vous pouvez voir que nous avons un répertoire dot Nextflow et un répertoire work et c'est parce que ceux-ci sont toujours créés dans votre répertoire de lancement. Donc ce sont différents des répertoires work et results que nous avions des exécutions précédentes.

Maintenant, si je regarde dans results, nous pouvons voir notre collected et voici notre petit personnage tux. Donc vous pouvez voir que ces paramètres ont été correctement interprétés.

## 1.3. Utiliser un fichier de paramètres

D'accord. Avant, quand je parlais des différents fichiers de configuration qui pouvaient être chargés, j'ai omis un autre endroit d'où nous pouvons obtenir la configuration.

Vous pouvez l'obtenir à partir d'une ligne de commande comme nous l'avons vu avec les noms de paramètres précédés de deux tirets, mais nous pouvons également fournir un fichier YAML ou JSON, juste de paramètres.

Le fichier de configuration peut avoir tous les différents types de scopes, mais ces fichiers ne contiennent que des paramètres, et c'est une manière conviviale de fournir de nombreux paramètres en une fois, et peut-être une façon un peu plus reproductible car vous les écrivez dans un fichier, il est donc facile de les récupérer ultérieurement.

Alors retournons à notre terminal et juste avant d'oublier, assurons-nous de remonter d'un répertoire, donc je ne suis plus dans le sous-répertoire, et je vais regarder le fichier YAML que nous avons ici appelé test-params.yaml.

Donc si je fais juste code test-params.yaml, vous pouvez voir que c'est juste un fichier YAML normal. Rien de spécial à ce sujet. Avec les clés étant nos noms de paramètres, avec le formatage YAML donc deux-points ici, puis une valeur.

Notez que ce n'est pas du code Nextflow, donc nous ne pouvons pas mettre des choses comme des variables ici. Ce sont juste des valeurs statiques.

Également parce que JSON est en fait analysé comme YAML, nous pouvons aussi avoir un fichier test-params.json, qui semble très similaire. C'est juste un format de données différent.

Donc nous avons deux fichiers de test différents ici et nous avons des variables légèrement différentes.

D'accord, alors comment les donnons-nous à Nextflow ? C'est très simple. Nous faisons Nextflow run hello config, comme avant. Et au lieu de "-c" pour fichier de configuration, ou charger ces noms de fichiers par défaut, nous faisons -params-file. Tiret simple car c'est une option Nextflow de base.

Ensuite nous passons le chemin de ce fichier. Donc je vais faire "-params-file test-params.yaml", et nous verrons si ceux-ci sont correctement chargés.

D'accord. Ça a fonctionné. Rappelons-nous simplement ce qu'il y avait dans ce fichier YAML. Donc le batch était défini sur YAML, donc c'est ce qu'il devrait s'appeler, et il devrait avoir un stégosaure. Alors allons-y et regardons dans results. Et nous avons COLLECTED-yaml. Alors voyons si nous avons un Stégosaure. Fantastique, un Stégosaure portant un chapeau. C'est ce que nous aimons.

Donc cela a très bien fonctionné, et c'est exactement pareil avec le fichier JSON. Nous échangeons juste l'extension du fichier ici et Nextflow sait comment lire cela.

Et dans ce cas, nous devrions avoir un batch appelé JSON et nous devrions avoir une tortue. Jetons un coup d'œil. Merveilleux. Un de mes outils CLI préférés.

## 2.1. Personnaliser le répertoire de sortie avec -output-dir

D'accord, donc cela concernait principalement les entrées du pipeline et la modification des paramètres. Qu'en est-il des sorties ?

Maintenant, bien que nous ayons modifié les sous-répertoires en utilisant params, vous avez peut-être remarqué que tous nos fichiers vont toujours dans results.

Nous pouvons changer ce répertoire de base où tous les fichiers sont publiés avec un flag de ligne de commande appelé -output-dir. Donc si je fais Nextflow run hello config, puis je fais -output-dir, et nous allons l'appeler "custom-outdir-cli". Je ne sais pas taper. Juste pour que nous nous souvenions d'où viennent ces fichiers.

C'est une option Nextflow de base et c'est une option très nouvelle. Cela n'a été ajouté que récemment, et c'est l'une des choses que nous pouvons faire avec le nouveau parseur de langage et tout.

C'est un peu long à taper. Vous pouvez aussi simplement l'appeler "-o" si vous voulez. Donc si je reviens juste en arrière. Je peux juste raccourcir ça en "-o", ce qui est un peu plus simple.

D'accord. Nous lançons cela. Nous n'avons rien changé dans notre pipeline ni même dans notre configuration à ce stade, et cela devrait espérons-le sauvegarder tous nos résultats dans un répertoire de niveau supérieur différent. Et vous pouvez imaginer que vous pouvez définir cela sur pratiquement n'importe quel chemin que vous voulez.

Il vient d'arriver en haut. Nous avons un custom-outdir-cli, et tous les fichiers y sont organisés exactement de la même manière, avec leurs mêmes sous-répertoires et noms de fichiers. C'est donc un moyen vraiment facile de changer simplement où le pipeline publie ses résultats, sans trop réfléchir à la façon dont ces résultats sont organisés.

## 2.1.2. Supprimer les chemins codés en dur du bloc output

Si je regarde dans ce répertoire, nous pouvons voir que nous avons toujours un sous-répertoire appelé Hello Config, ce qui semble un peu redondant maintenant.

Alors chargeons simplement notre script à nouveau et nous pouvons maintenant supprimer ce sous-répertoire du bloc output en bas. Car nous n'en avons plus vraiment besoin. Donc nous pouvons juste faire cela maintenant, supprimer cela d'ici. Et ensuite si c'est juste ça, vous pouvez soit supprimer cela complètement, soit le laisser comme une chaîne vide. Je vais le laisser comme une chaîne vide pour l'instant, parce que nous allons revenir et mettre différentes choses à sa place à l'avenir. Mais si vous ne vous souciez pas des sous-répertoires, il est plus propre de supprimer complètement la déclaration path là.

D'accord, sauvegardons. Essayons-le rapidement à nouveau. Je vais en fait supprimer mon répertoire "custom-outdir-cli" pour que nous ne soyons pas confus par des fichiers existants là. Car rappelez-vous, quand vous publiez des choses, cela ne supprime pas les fichiers qui étaient déjà là. Cela ajoute simplement de nouveaux fichiers. Lançons à nouveau cette commande, custom-outdir-cli.

Et maintenant si vous faites "ls custom-outdir-cli", il n'y a plus de répertoire là appelé Hello Config.

## 2.2.1. Définir outputDir dans le fichier de configuration

D'accord, le flag de ligne de commande ici, "-o" ou "-output-dir" est bien. Mais comment définir des valeurs par défaut pour cela dans la configuration ? Comment faisons-nous cela ?

J'ouvre le fichier "nextflow.config", ferme tout le reste et me débarrasse de ça. Nous pouvons ajouter une nouvelle option de configuration ici, que j'ai juste copiée du site web du matériel de formation, et elle s'appelle outputDir.

Elle n'est sous aucun scope. Elle n'est pas sous params ou quoi que ce soit. Elle est au niveau supérieur, et nous pouvons la définir sur une chaîne. Maintenant une chose simple à faire est de la changer simplement en n'importe quoi d'autre que results comme chaîne codée en dur. Mais parce que c'est dans un fichier de configuration Nextflow, nous pouvons être un peu malins ici et également inclure des variables.

Et vous pouvez voir ici que nous avons inclus une variable params, params.batch, qui fait partie de cette chaîne. Cela signifie que nous pouvons réutiliser des variables qui proviennent d'autres endroits. Et dans ce cas, si nous faisons --batch, lorsque nous exécutons le Pipeline Nextflow, nous allons obtenir un sous-répertoire dans notre chemin personnalisé basé sur le nom du batch.

D'accord, alors essayons cela et jetons juste un coup d'œil rapide pour voir comment, comment les résultats apparaissent. Donc si je fais Nextflow run hello config et --batch my_run. Rappelons-nous à quoi ressemblait la configuration. Donc c'est custom-outdir-config.

Tree custom-outdir-config. Et vous pouvez voir que le batch s'appelait my_run. Et puis nous avons ce sous-répertoire appelé my_run. Donc ce chemin de fichier dynamique a fonctionné.

Et en plus de cela, il n'est plus allé dans un répertoire results par défaut, et je n'ai pas eu à spécifier quoi que ce soit sur la ligne de commande pour changer le répertoire de base. Donc nous avons réussi à réinitialiser la valeur par défaut pour le outputDir par défaut.

## 2.2.2. Sous-répertoires avec les noms de batch et de processus

D'accord, allons un peu plus loin. C'est une variable dynamique dans le fichier de configuration. Qu'en est-il du script lui-même ? Maintenant, jusqu'à présent nous avions ces chemins ici et ceux-ci peuvent également être dynamiques. Donc au lieu de simplement coder en dur quelque chose, nous pouvons mettre des accolades et mettre quelque chose de dynamique.

Par exemple, nous avons nos processus appelés sayHello. Nous pourrions faire sayHello.name, qui est un attribut du processus, ce qui est un peu ennuyeux car c'est juste "sayHello" dans ce cas. Mais c'est variable.

Donc cela vous donne une idée. Donc nous pouvons mettre cela ici et dire convertToUpper.name, collectGreetings.name, collectGreetings.name à nouveau, et cowpy.

Maintenant quand nous exécutons, le répertoire de base va toujours être custom-outdir-config. Et il va être dans un sous-répertoire appelé params.batch, mais les sous-répertoires en dessous devraient être organisés par nom de processus.

Essayons simplement cela et voyons si ça fonctionne. Donc je vais supprimer le répertoire précédent pour que nous ne soyons pas confus, et utiliser exactement la même commande Nextflow Run.

Cela devrait s'exécuter de la même manière. J'aurais pu utiliser dash resume sur tous ces éléments pour le rendre un peu plus rapide et utiliser les résultats calculés précédemment. Maintenant, si je fais tree custom-outdir-config, vous pouvez voir qu'il n'est pas dans results, il est dans notre répertoire de base avec le nom du batch. Et vous pouvez voir que tous les résultats sont maintenant organisés dans des sous-répertoires nommés d'après le processus. Donc nous avons deux endroits différents où nous définissons des chemins de sortie dynamiques ici.

D'accord. Dernière chose, rajoutons ces dossiers intermédiaires, que nous avions avant car ils étaient plutôt sympas. Intermediates.

Et nous pouvons également réfléchir un peu à ce params.batch, peut-être qu'en tant que développeur·se de pipeline j'aimais vraiment avoir cela dans le sous-répertoire, mais si les utilisateur·trices finaux·ales du pipeline définissent "-o" ou -output-dir sur le CLI, cela écrase complètement cette déclaration entière, et nous perdons ce sous-répertoire.

Donc ce que nous pouvons faire, c'est sortir ce chemin dynamique du outputDir config, qui serait écrasé, et le mettre dans le chemin output, qui n'est pas écrasé.

Donc nous pouvons faire params.batch slash intermediates slash sayHello.name, et faire tout cela dans une chaîne entre guillemets doubles, donc c'est interpolé par Nextflow.

On peut maintenant copier, oups. Copier ceux-ci vers les autres processus. N'oubliez pas de tous les mettre entre guillemets. Et supprimer intermediates de ces sorties particulières.

D'accord ? Cela a l'air légèrement plus complexe maintenant, mais vous pouvez voir que nous commençons vraiment à construire une belle structure de répertoire de sortie organisée dans notre code.

Et ce qui est vraiment bien, c'est que cette complexité supplémentaire dans le code ne passe pas au CLI. Donc nous pouvons exécuter notre commande avec -output-dir et quelles que soient les variables batch, en pensant juste à comment exécuter le pipeline et pas vraiment en pensant trop à ce qui est dans le code. Et nos fichiers de sortie vont être construits vraiment bien d'une manière très bien organisée, ce qui est agréable pour les personnes utilisant le pipeline en gros.

Parfait. En écrivant ceci, je réalise que j'ai fait une erreur. Voyons si quelqu'un m'a pris en flagrant délit ici. Nous avons collectGreetings.name, donc quelque chose a mal tourné. Et ouais, effectivement, j'ai accidentellement oublié de mettre ceux-ci dans des accolades.

Donc rappelez-vous, faites attention quand vous écrivez votre code et assurez-vous de dire à Nextflow ce qui est une variable et ce qui est juste une chaîne. Car il fera exactement ce que vous lui dites de faire. Et rien de plus. Comme tous les bons ordinateurs. D'accord, cela devrait le corriger.

## 2.3. Définir le mode de publication au niveau du workflow

Il y a une partie de ce script que je n'aime toujours pas, c'est le fait que nous écrivons mode copy encore et encore, et s'il y a une chose que nous n'aimons pas, c'est nous répéter.

Donc nous pouvons nettoyer un peu cela en prenant cela et en le déplaçant dans la configuration. Et en fait, nous pouvons le définir pour tout le pipeline en une seule fois. Donc nous n'avons pas à le dire plusieurs fois.

Nous allons à notre fichier de configuration et nous avons un nouveau scope ici appelé workflow. Et nous pouvons soit faire des accolades, soit utiliser la notation par points. Cela ne fait aucune différence. Le site web du matériel de formation utilise la notation par points. Je peux dire output et nous pouvons mélanger et assortir, donc mode equals copy. Parfait.

Et maintenant nous pouvons revenir ici et supprimer ceux-ci. Maintenant nous pourrions les laisser en place. La configuration écrase essentiellement ce qui est écrit ici, mais comme nous l'avons dans la configuration au niveau du pipeline, et que ces deux fichiers sont livrés ensemble, il n'y a vraiment aucune raison de le faire deux fois.

D'accord. Juste une vérification de bon sens, parce qu'apparemment nous faisons des erreurs. Lançons cela à nouveau et vérifions simplement que nous utilisons correctement le mode copy pour publier les fichiers. Donc nous allons exécuter le script à nouveau et cette fois nous avons mis les résultats dans un répertoire appelé config-output-mode, voyons à quoi ressemblent les fichiers là-dedans.

Et puis si je fais "ls -l" pour regarder batch, et nous pouvons regarder cowpy, par exemple. Et nous devrions voir, ouais, que c'est un fichier approprié ici, qui n'est pas un lien symbolique, donc cet attribut de configuration a été appliqué correctement.

## 3. Sélectionner une technologie d'empaquetage logiciel

D'accord. Jusqu'à présent, nous nous sommes concentré·e·s sur les entrées et les sorties, les fichiers avec lesquels le workflow s'exécute. Mais qu'en est-il de l'infrastructure ? J'ai dit au début que Nextflow vous permet d'exécuter le même pipeline sur différentes configurations informatiques. Alors à quoi cela ressemble-t-il ?

Pour montrer cela, nous allons passer de l'utilisation de Docker pour exécuter cowpy, et à la place nous utiliserons Conda pour faire la même chose.

Je peux faire cela très simplement. Si je vais à code, "nextflow.config". Si vous vous souvenez en haut, nous avons défini docker.enabled plus tôt, dans le dernier chapitre pour que nous puissions utiliser le conteneur avec cowpy dedans.

Je vais dire à Nextflow de ne pas utiliser Docker. Définir cela sur false. Et je vais dire Conda enabled equals true. Donc dire à Nextflow, veuillez utiliser Conda.

Maintenant, activer Conda n'est pas suffisant en soi. Tout comme nous l'avons fait avec Docker, nous devons dire à Nextflow où il peut obtenir le logiciel dont il a besoin.

Donc si nous allons dans les modules ici. Et ouvrons le script cowpy. Nous pouvons voir que nous avons une déclaration container en haut. Et le conteneur est utilisé par Docker, mais aussi Singularity, Apptainer, et beaucoup d'autres outils logiciels.

Mais il ne peut pas être utilisé pour Conda, donc nous avons une déclaration séparée appelée "conda", et nous pourrions simplement écrire "cowpy". Et cela laissera à la résolution de paquets conda le soin de déterminer la meilleure façon de résoudre cela, selon votre environnement conda local.

Ou c'est une bonne pratique de faire ce que le site web du matériel de formation dit de faire, qui est de définir un canal conda spécifique avec sa notation à doubles deux-points, et définir certainement une version spécifique du logiciel pour que chaque personne qui exécute le pipeline obtienne la même version.

Notez que les conteneurs sont un peu supérieurs à cet égard, car lorsque vous installez quelque chose avec Conda, il va toujours déterminer toutes les dépendances pour ce paquet, et elles peuvent changer avec le temps. C'est ce qu'on appelle la dérive des dépendances.

Donc les conteneurs, cependant, verrouillent toute la pile de dépendances logicielles jusqu'en bas, donc vous pouvez être un peu plus confiant·e que A, ça va fonctionner, et B, ce sera reproductible.

Donc si vous êtes en mesure d'utiliser Docker ou Singularity ou Apptainer, je recommanderais certainement cela.

Maintenant ce qui est bien à propos de cela, c'est que le fichier de module, qui est écrit par le développeur·se de pipeline, a maintenant à la fois Container et Conda, et donc nous disons à la personne qui exécute ce pipeline, peu nous importe quelle solution d'empaquetage logiciel vous utilisez. Cela fonctionnera à la fois avec Docker et avec Conda, et voici où obtenir le logiciel dans les deux cas.

Nous pouvons ouvrir le terminal et essayons cela. Donc Nextflow run hello config --batch conda. Et la première fois que cela s'exécute avec conda, ça va être un peu lent quand ça arrive à ce processus particulier, parce qu'il doit exécuter "conda install".

Et il crée un environnement conda spécial juste pour ce seul processus. Donc il n'utilise pas mon environnement conda global, que j'ai sur mon terminal. Il en crée un juste pour ce seul processus. C'est bien parce que cela évite des choses comme les conflits de dépendances entre différents processus de votre workflow. Si vos processus ont des outils qui nécessitent différentes versions de Python ou des choses comme ça, ce n'est pas grave car ils utilisent différents environnements conda.

Nextflow met en cache ces environnements conda localement, vous pouvez voir qu'il vous dit où se trouve ce chemin, c'est dans le répertoire work ici. Et donc la prochaine fois que j'exécuterai ce script avec Conda, ce sera beaucoup plus rapide car il trouvera cet environnement conda existant et le réutilisera simplement. Mais la première fois que nous le faisons, il doit aller le chercher, le résoudre, télécharger toutes les dépendances et tout configurer.

D'accord, super, ça a fonctionné. Nous pouvons juste nous rappeler ce que le pipeline est actuellement configuré pour utiliser. Si nous regardons dans le fichier de configuration, c'était "custom-outdir-config" pour moi en ce moment. Voyons si je vais à ce répertoire de base. Et j'ai fait --batch conda. Voici notre sous-répertoire conda. Donc ça a fonctionné et voici notre sortie cowpy.

Donc il a récupéré cowpy, l'a installé sur mon système local en utilisant conda, et a exécuté le processus. Et ce qui est génial, c'est qu'en tant qu'utilisateur·trice final·e, je n'ai pas eu à penser du tout à la gestion du logiciel là. Nextflow l'a juste arrangé pour moi. J'ai dit, j'ai besoin d'utiliser conda sur ce système. Le développeur·se du pipeline a dit quels paquets j'avais besoin. Et Nextflow a fait le reste. Très puissant.

Notez que vous pouvez en fait utiliser un mélange de différentes technologies. Donc je peux activer Docker pour des processus spécifiques, et conda pour d'autres processus, ou dire que certains processus devraient juste utiliser le logiciel local que j'avais installé. C'est assez inhabituel, mais c'est possible, et dans certains cas, par exemple, si vous utilisez certains logiciels qui pourraient être difficiles à empaqueter dans Docker, vous avez une échappatoire.

## 4. Sélectionner une plateforme d'exécution

Donc c'est l'empaquetage logiciel. L'autre partie de la portabilité vers d'autres systèmes est l'endroit où les tâches s'exécutent réellement. Pour le moment, je m'exécute sur essentiellement mon ordinateur portable ou dans ce Codespaces, qui est un seul ordinateur. Il n'y a rien de compliqué. Nextflow est un peu malin dans la parallélisation des tâches du mieux qu'il peut, mais c'est tout sur un système.

Maintenant, si vous exécutez sur un HPC, vous avez probablement une sorte de planificateur de tâches comme SLURM ou PBS ou quelque chose, et vous soumettez des tâches à ce planificateur et il répartit toutes les tâches sur différents nœuds de calcul.

Une autre façon d'exécuter est sur le cloud. Donc peut-être que vous utilisez AWS Batch, ou Azure Cloud, ou Google. Et ceux-ci fonctionnent tous dans un système similaire où vous avez un planificateur et vous soumettez des tâches et elles sont soumises à différents endroits pour être calculées.

Maintenant dans un passé lointain quand j'ai commencé à faire de la bioinformatique, tous les logiciels de tout le monde pour exécuter des analyses étaient très liés à leur infrastructure informatique, ce qui rendait presque impossible la réplication.

Mais avec cette séparation de configuration dans Nextflow, et avec la capacité de Nextflow à interagir avec de très nombreux backends d'infrastructure de calcul différents, c'est très simple de prendre notre pipeline sans modifier du tout le code du pipeline et de simplement changer cela.

## 4.1. Cibler un backend différent

Donc si nous allons à notre fichier "nextflow.config", et nous pouvons maintenant mettre une configuration au niveau du processus. Donc si je mets en haut le scope process et je peux définir l'executor, et ici il est défini sur local, qui est la valeur par défaut.

Remarquez que parce que c'est au niveau du processus, nous pouvons cibler des choses vers différents processus. Et donc vous pouvez en fait configurer des executors pour être spécifiques à un processus et avoir une exécution hybride, où certaines tâches pourraient s'exécuter localement, là où la tâche Nextflow est exécutée. Certaines sont soumises à différents HPC et certaines pourraient être soumises au cloud. Vous pouvez être aussi malin·e que vous le souhaitez.

Maintenant, il est très difficile de faire une démo de cela dans un environnement de formation comme celui-ci parce que je n'ai pas de HPC auquel soumettre. Mais ce que je peux faire, c'est si je tape slurm, nous pouvons tricher un peu et vous pouvez avoir une idée de cela.

Et c'est vraiment le plus intéressant pour les personnes qui ont l'habitude d'exécuter sur SLURM et qui savent à quoi ressemblent les en-têtes SLURM. Mais si je fais Nextflow run, hello config. Cela va échouer parce qu'il va essayer de soumettre des tâches à un cluster qui n'existe pas. Donc nous aurons une sorte d'erreur à propos de sbatch qui n'est pas disponible.

Ouais, écrit. C'est l'outil. C'est l'outil CLI que vous utilisez pour soumettre des tâches à un cluster slurm. Mais ce que nous pouvons faire, c'est aller regarder dans notre répertoire work ici en faisant command clic, ouvrir ce répertoire et regarder le .command.run. Et vous pouvez voir en haut du fichier .command.run, nous avons nos en-têtes sbatch, disant à un cluster SLURM théorique comment gérer cette soumission de tâche.

Et donc vous pouvez voir que Nextflow est malin, il fait toutes les bonnes choses. C'est juste que nous n'avions pas de cluster auquel soumettre.

## 5. Contrôler les allocations de ressources de calcul

Qu'est-ce qui est encore différent entre différentes infrastructures informatiques ? Une autre chose est la quantité de ressources disponibles que vous avez, et en fait, dans de nombreux environnements de calcul, c'est une exigence que vous devez spécifier combien de CPUs et quelle quantité de mémoire une tâche nécessite.

Encore une fois, Nextflow abstrait cela pour nous, de sorte qu'il n'est plus spécifique à un seul type d'environnement de calcul, et nous pouvons taper dans le scope au niveau du processus ici. CPUs equals un, memory equals deux gigaoctets. Notre pipeline n'est pas très exigeant, donc cela devrait aller.

Maintenant, j'ai juste deviné ces chiffres ici, mais comment savez-vous quelle est une quantité sensée de ressources à utiliser ? C'est un travail assez difficile d'aller fouiller dans tous ces différents processus d'un gros pipeline de nombreux échantillons et de comprendre quelle était l'utilisation des ressources.

Donc une bonne approche pour cela est de définir ces valeurs sur des nombres élevés pour commencer, juste pour que votre pipeline s'exécute sans aucune erreur, puis demander à Nextflow de générer un rapport d'utilisation pour vous.

C'est super facile à faire, donc je vais retourner à un terminal. Oh, je dois me rappeler de remettre ça sur local pour que mon pipeline s'exécute réellement. Et je vais dire Nextflow run, et je vais utiliser un flag de ligne de commande -with-report.

Et je peux laisser cela vide et il donnera un nom de fichier par défaut, mais je vais lui donner un nom de fichier spécifique pour que cela soit sauvegardé à un endroit spécifique.

Appuyons sur Entrée, et le pipeline s'exécute exactement comme d'habitude, mais quand il se termine, il va générer un beau rapport HTML pour moi.

Donc dans la barre latérale ici, j'ai ce fichier HTML. Si j'exécutais cela localement, je l'ouvrirais simplement. Parce que je suis dans Codespaces, je vais faire un clic droit dessus et cliquer sur télécharger, ce qui va le télécharger sur mon ordinateur local. Et je peux simplement l'ouvrir facilement dans le navigateur web.

Nextflow peut générer un rapport comme celui-ci pour n'importe quel pipeline et il contient des informations vraiment intéressantes. Donc c'est une bonne pratique de toujours sauvegarder ces choses. Il nous dit quand nous avons exécuté, où nous avons exécuté, si c'est réussi ou non, quels paramètres ont été utilisés, quelle était la commande CLI, des choses comme ça.

Et il y a aussi ces graphiques sur l'utilisation des ressources. Donc il nous dit quel pourcentage d'appels CPU ont été utilisés pour chaque processus sous forme de diagramme en boîte ici, parce qu'il y a de nombreuses tâches pour chaque processus, donc nous pouvons voir la distribution.

Vous pouvez voir nos processus ici, cowpy et collectGreetings n'avaient qu'une seule tâche, donc c'est juste une seule ligne. Et nous avons à la fois CPU et mémoire et durée de tâche, et ils étaient très rapides.

Si vous utilisez Seqera Platform, soit dit en passant, vous obtenez les mêmes graphiques intégrés dans l'interface Platform sans avoir à faire quoi que ce soit. Donc vous obtenez toujours cette information à portée de main.

D'accord, donc nous pouvons utiliser ce rapport et sur une vraie exécution, et avoir une idée de combien de CPUs et combien de mémoire sont utilisés par notre pipeline et revenir et mettre ces valeurs dans notre fichier de configuration, pour que la prochaine fois peut-être nous ne demandions pas autant. Et nous pouvons être un peu plus économes.

Maintenant vous pouvez devenir vraiment malin·e dans la configuration des fichiers de configuration de pipeline. Et encore une fois, si vous utilisez Seqera Platform, recherchez un petit bouton qui ressemble à une ampoule. Car si vous cliquez dessus, il générera un fichier de configuration hautement optimisé, qui est spécifiquement adapté à vos données, votre exécution et votre pipeline. Pour l'exécuter de la manière la plus efficace possible.

Mais pour l'instant, je vais dire qu'en fait le nombre par défaut de CPUs que Nextflow donnait était bien et mais nous n'avons besoin que d'un gigaoctet de mémoire.

## 5.3. Définir les allocations de ressources pour un processus spécifique

Maintenant, dans la vraie vie, il est assez inhabituel que tous les processus de votre pipeline aient besoin des mêmes exigences. Vous pourriez avoir quelque chose comme MultiQC comme outil de reporting, qui nécessite très peu en termes de ressources et s'exécute assez rapidement.

Et puis peut-être que vous avez quelque chose qui indexe un génome de référence ou qui fait un alignement ou qui fait un autre travail. Peu importe ce que c'est, qui prend beaucoup de ressources. Et donc pour ces différentes soumissions de tâches à un planificateur, vous voulez donner différentes quantités de ressources.

Sous ce scope process, nous pouvons définir une configuration, qui cible des processus spécifiques de différentes manières.

Ici nous utilisons withName, nous pouvons aussi utiliser des labels, et ceux-ci peuvent utiliser un motif pour cibler un ou plusieurs processus. Ici nous disons simplement que tous les processus qui ont un nom cowpy sont définis sur deux gigaoctets de mémoire et deux CPUs, et parce que c'est un sélecteur plus spécifique que le processus de niveau supérieur, ceci est écrasé dans ces cas, donc vous pouvez construire un beau fichier de configuration ici, qui adapte vraiment tous vos différents processus dans votre pipeline pour les rendre vraiment efficaces.

## 5.5. Ajouter des limites de ressources

Maintenant en tant que développeur·se de pipeline, je connais probablement assez bien les outils, et je veux que tout s'exécute aussi rapidement et aussi bien que possible. Donc il se peut que je mette des chiffres assez élevés pour certains d'entre eux parce que je sais que ça s'exécutera beaucoup plus rapidement si je donne à cowpy 20 CPUs au lieu de deux.

C'est bien jusqu'à ce que vous alliez exécuter sur votre ordinateur portable ou sur les tests d'Intégration Continue GitHub Actions, ou un autre système, qui n'a peut-être pas 20 CPUs disponibles.

Maintenant quand vous essayez d'exécuter le pipeline, il va planter parce que Nextflow dira, je ne peux pas soumettre cette tâche nulle part. Je n'ai pas les ressources disponibles.

Maintenant pour éviter ce plantage dur, nous pouvons ajouter un peu plus de configuration, qui est spécifique à notre système maintenant, appelée limites de ressources. Et cela ressemble à ceci. C'est sous le scope process à nouveau.

Et limites de ressources, vous pouvez spécifier essentiellement le plafond de ce que vous avez disponible. C'est une map ici, et vous pouvez, dans cette map, vous pouvez définir la mémoire, les CPUs et le temps.

Maintenant ce qui se passe, c'est que lorsque Nextflow soumet une tâche d'un processus, il regarde ce qui est demandé et il fait essentiellement un minimum entre cela et cela. Donc si nous avons demandé 20 CPUs, mais seulement quatre sont disponibles, il demandera quatre. Le pipeline ne plante pas et il utilise aussi près que possible de ce qui a été conçu par le développeur·se du pipeline.

## 6. Utiliser des profils pour basculer entre des configurations prédéfinies

D'accord. J'ai dit que les limites de ressources ici pourraient être spécifiques au système, et peut-être que j'ai un fichier Nextflow config dans mon pipeline, et je sais que les gens vont l'utiliser dans une gamme d'endroits différents. Maintenant, au lieu de forcer tout le monde à créer son propre fichier Nextflow config à chaque fois, ce que je peux faire, c'est regrouper différents préréglages de configuration ensemble dans des profils de configuration.

Je vais descendre un peu ici et puis juste après params, parce que l'ordre du fichier de configuration ici est important, le fichier de configuration est chargé séquentiellement, donc je vais mettre ces profils après tout le reste pour qu'ils écrasent les params précédemment définis. Et je vais coller ces profils du matériel de formation.

Donc il y a un nouveau scope de niveau supérieur appelé profiles. Nous pouvons avoir des noms arbitraires ici. Donc nous avons my_laptop et univ_hpc. Et ici nous pouvons voir que nous définissons les mêmes paramètres de configuration qu'avant. Maintenant juste dans un profil. Donc nous avons un executor local pour exécuter sur my_laptop et je soumets à un cluster SLURM sur le HPC.

J'utilise Docker localement, conda sur le HPC, et le système HPC a des limites de ressources beaucoup plus élevées.

Maintenant je peux exécuter le pipeline avec l'option CLI -profile, dire quel profil je veux utiliser. Donc je vais utiliser my_laptop, et Nextflow appliquera toute la configuration dans ce scope de profil. Donc je peux essayer cela maintenant. C'est la même commande qu'avant. Nextflow run hello config, et je fais dash profile, tiret simple car c'est l'option Nextflow de base, dash profile my_laptop.

Il va maintenant appliquer par lot toute cette option de configuration. Oh, et vous pouvez voir, j'ai dit avant que cela pourrait arriver que l'exigence du processus, il a demandé quatre CPUs et je n'en ai que deux sur cette instance Codespaces.

Donc c'est une bonne opportunité juste pour essayer les limites de ressources de processus, et dire que je n'ai que deux CPUs sur my_laptop, ou dans ce Codespaces. Maintenant si nous l'exécutons à nouveau, cela devrait plafonner cette exigence à deux et espérons-le le pipeline s'exécutera. Parfait.

## 6.2. Créer un profil de paramètres de test

Notez que ces profils n'ont pas à avoir uniquement une configuration concernant leur infrastructure. Vous pouvez avoir des regroupements de toute configuration ici, y compris des paramètres.

Donc une autre chose que vous verrez très souvent dans les pipelines des gens est un profil de test, qui inclut des paramètres, que vous soumettriez normalement par utilisateur·trice. Mais ici nous avons, essentiellement différentes valeurs par défaut sensées pour quand je veux exécuter des cas de test.

Et c'est génial parce que je n'ai pas nécessairement à aller spécifier toutes ces choses, qui pourraient être des paramètres requis. Sinon je peux juste dire dash profile test et ça fonctionnera directement.

Maintenant quelque chose à noter, c'est que les profils peuvent aussi être combinés plus d'un. Donc je peux faire profile my_laptop ici, et ensuite aussi ajouter test. Je ne fais pas profile deux fois. Je fais juste une liste séparée par des virgules ici sans espaces. Et il va appliquer ces profils dans l'ordre. Donc il prendra la configuration du profil my_laptop, puis il appliquera la configuration de test par-dessus.

Vraiment pratique et vous pouvez voir comment vous pouvez configurer beaucoup de groupes par défaut sensés ici pour faciliter l'exécution de votre pipeline.

## 6.3. Utiliser nextflow config pour voir la configuration résolue

Espérons-le, je vous ai convaincu·e que la résolution de configuration Nextflow est puissante, mais je ne vous en voudrais pas si vous avez un peu les yeux qui se croisent à ce stade après que j'ai dit environ 20 façons différentes de fournir la configuration et donné toutes ces différentes couches comme une pelure d'oignon.

Donc si jamais vous vous sentez incertain·e de ce qu'est la configuration résolue finale pour Nextflow, sachez qu'il y a une commande appelée "nextflow config", et nous pouvons l'exécuter et elle nous dira quelle est la configuration résolue à notre emplacement actuel.

Donc quand je l'exécute ici, il trouve le fichier "nextflow.config" dans le répertoire de travail actuel, et il traite toute la configuration différente, et il me donne la sortie résolue.

Notez que le fichier de configuration Nextflow peut également prendre l'option CLI profile. Donc si je lui dis de résoudre dans les profils my_laptop et test, et vous pouvez voir qu'il a également appliqué les limites de ressources ici de l'option de configuration my_laptop et également défini les params, qui étaient dans le test.

Donc c'est une façon agréable d'explorer simplement comment fonctionne la résolution de configuration, si vous avez le moindre doute.

## Conclusion

D'accord, c'est tout. C'est la configuration Nextflow en bref. Vous pouvez faire beaucoup de choses avec la configuration. C'est vraiment puissant. Mais ce sont la plupart des cas d'usage courants que vous vous retrouverez à faire, et ces concepts s'appliquent à toutes les différentes options.

Félicitez-vous car c'est la fin du cours de formation Hello Nextflow. Vous êtes espérons-le maintenant confiant·e à la fois pour écrire votre propre pipeline Nextflow à partir de zéro, le configurer et l'exécuter, et vous connaissez tous les tenants et aboutissants et les choses à surveiller.

Il y a un quiz de plus que vous pouvez essayer sur la page de formation sur la configuration. Alors descendez et essayez-le et assurez-vous que vous avez compris toutes ces parties sur la configuration.

Et rejoignez-nous dans la dernière vidéo juste pour une rapide conclusion sur certaines des prochaines étapes qui pourraient être bonnes à faire après ce cours de formation.

Merci d'être resté·e·s avec nous. Bravo et je vous verrai dans la prochaine vidéo.
