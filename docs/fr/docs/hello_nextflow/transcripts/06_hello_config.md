# Partie 6 : Hello Config - Transcription vidéo

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Traduction assistée par IA - [en savoir plus et suggérer des améliorations](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

<div class="video-wrapper">
  <iframe width="560" height="315" src="https://www.youtube.com/embed/FcZTiE25TeA?si=tnXTi6mRkITY0zW_&amp;list=PLPZ8WHdZGxmWKozQuzr27jyMGqp9kElVK&amp;cc_load_policy=1" title="YouTube video player" frameborder="0" allow="accelerometer; autoplay; clipboard-write; encrypted-media; gyroscope; picture-in-picture; web-share" referrerpolicy="strict-origin-when-cross-origin" allowfullscreen></iframe>
</div>

!!!note "Notes importantes"

    Cette page présente uniquement la transcription. Pour les instructions détaillées étape par étape, retournez au [matériel de formation](../06_hello_config.md).

    Les numéros de section indiqués dans la transcription sont fournis à titre indicatif uniquement et peuvent ne pas inclure tous les numéros de section présents dans le matériel.

## Bienvenue

Bonjour et bienvenue dans la Partie Six de Hello Nextflow. Cette section est entièrement consacrée aux configurations, et c'est la dernière partie de ce cours.

Nextflow excelle particulièrement dans deux domaines : la reproductibilité et la portabilité. Les configurations sont l'endroit où nous voyons vraiment briller la seconde de ces qualités. La capacité de configurer un pipeline Nextflow pour qu'il s'exécute de différentes manières et fonctionne sur différents systèmes, sans avoir à modifier le code sous-jacent du pipeline.

Ce super-pouvoir permet aux pipelines Nextflow d'être réutilisés par d'autres personnes dans différents endroits, ou sur différentes infrastructures auxquelles vous pourriez avoir accès vous-même.

Cela signifie que vous pouvez développer le code du pipeline sur votre ordinateur portable, le pousser vers le cloud, l'exécuter sur votre HPC, et c'est le même code de pipeline qui s'exécute partout.

Dans cette section, nous allons aborder quelques sujets. Nous commencerons par la façon dont Nextflow gère les fichiers de configuration, d'où il les charge, comment vous les écrivez et comment vous les structurez, et cette séparation entre le pipeline lui-même et ce qui devrait aller dans un fichier de configuration.

Ensuite, nous passerons à quelques cas d'usage courants tels que la modification de l'emplacement de stockage des fichiers de sortie, et aussi comment faire fonctionner le pipeline sur différentes infrastructures, en utilisant différents types d'empaquetage logiciel ou en soumettant des tâches à différentes infrastructures.

## Hiérarchies des fichiers de configuration

D'accord, commençons. En ce qui concerne le chargement des fichiers de configuration, Nextflow peut extraire des informations de nombreux endroits différents, ce qui est une bonne chose et peut aussi être une chose légèrement risquée car parfois il peut être un peu difficile de savoir d'où il obtient un fichier de configuration et dans quel ordre il charge les choses.

Je vous recommande donc vraiment de cliquer sur ce lien ici, qui nous amène à la documentation Nextflow. Et sur cette page de configuration, elle liste les principaux endroits d'où la configuration est chargée, et surtout, l'ordre dans lequel ces éléments sont chargés.

Vous pouvez voir que vous pouvez placer un fichier de configuration dans votre répertoire home Nextflow, qui est généralement ".nextflow" dans votre répertoire home. Et ce fichier sera toujours chargé par chaque exécution Nextflow sur votre système.

L'endroit suivant à examiner est un fichier à la racine de votre dépôt ou répertoire de pipeline appelé "nextflow.config".

Ensuite, un autre fichier appelé "nextflow.config", mais cette fois dans le répertoire à partir duquel vous lancez Nextflow : le répertoire de lancement.

Enfin, vous pouvez fournir des chemins de fichiers de configuration sur la ligne de commande avec un argument "-c", et vous pouvez le faire plusieurs fois. Et ils sont appliqués dans l'ordre que vous spécifiez.

Vous pouvez fournir des fichiers de configuration dans tous ces emplacements si vous le souhaitez, et ils seront chargés de manière itérative, chacun écrasant le précédent uniquement dans les portées de configuration où ils entrent en conflit.

C'est un système vraiment puissant car cela signifie que vous pouvez définir des valeurs par défaut sensées, puis devenir progressivement de plus en plus spécifique à mesure que vous affinez cette configuration.

## 0. Échauffement : Exécuter hello-config.nf

D'accord, fermons ceci et passons à notre Codespaces pour commencer. Comme précédemment, j'ai nettoyé ici, j'ai supprimé mes répertoires de résultats précédents, mes répertoires Nextflow et work, etc. Ne vous inquiétez pas si vous avez encore ces fichiers qui traînent. C'est juste parce que je suis très zoomé et donc les choses deviennent désordonnées très rapidement sinon.

Nous allons travailler avec hello-config.nf, le dernier fichier de notre répertoire, et cela devrait faire suite à ce que nous avons laissé dans la section précédente.

Nous avons donc nos quatre processus différents, qui sont inclus à partir de fichiers de modules. Nous avons nos paramètres de pipeline, notre bloc workflow où nous appelons les différents processus et assemblons les canaux, publions les canaux de sortie, puis le bloc output en bas où nous définissons où ces fichiers doivent être stockés et comment ils doivent être copiés.

Nous avons également déjà un fichier "nextflow.config" du dernier chapitre, où nous activons Docker, et nous allons enrichir ce fichier aujourd'hui.

Comme précédemment, nous avons changé le chemin de sortie dans ce script principal en hello config, juste pour qu'il n'entre pas en conflit avec les résultats précédents que vous avez générés.

D'accord, vérifions rapidement que tout fonctionne toujours comme prévu. Ouvrons un terminal et faisons nextflow run hello-config.nf. Nextflow se charge. Devrait exécuter nos quatre processus différents. Générer de jolis dessins ASCII en utilisant cowpy puis sauvegarder nos résultats dans nos fichiers de résultats dans ce répertoire.

Je peux jeter un coup d'œil rapide ici juste pour m'assurer que ces fichiers ressemblent à ce que nous attendons, et effectivement, voici notre dinde géante. Parfait.

## 1.1. Déplacer les valeurs par défaut vers nextflow.config

Maintenant, la première chose que nous allons faire est de commencer à déplacer certaines choses de notre script vers notre fichier de configuration.

Et ce qui nous intéresse, ce sont principalement les paramètres à ce stade. Nous voulons prendre les valeurs par défaut et les mettre dans le fichier de configuration, pour qu'il soit plus clair quelles sont les valeurs par défaut et qu'il soit plus facile pour les gens de les écraser.

Je vais juste prendre ce bloc params ici du script et le mettre dans le fichier de configuration. Et nous devons être un peu prudent·e·s ici, car pour le moment la syntaxe est légèrement différente entre la configuration et les scripts. Le fichier de configuration ne peut pas prendre de déclarations de type car nous ne définissons pas vraiment ces paramètres, nous les référençons simplement. Je vais donc me débarrasser de ceux-ci.

Mais sinon, c'est très similaire. Nous avons un bloc params puis nous avons nos différents paramètres d'entrée, paramètre batch, paramètre character.

Je peux maintenant retourner à mon script et je n'ai plus besoin de définir ces valeurs par défaut car ces valeurs sont maintenant dans mon fichier Nextflow config.

Cependant, je laisse les noms de paramètres et leurs types, pour que Nextflow connaisse ces informations et puisse toujours faire toute la vérification de type et tout le reste.

D'accord. Nous sauvegardons ces fichiers et vérifions rapidement que tout fonctionne toujours de la même manière qu'avant. Il ne devrait y avoir aucun changement ici. Nous avons conservé les mêmes valeurs. Nous avons juste déplacé l'endroit où elles ont été définies.

Parfait.

## 1.2. Utiliser un fichier de configuration spécifique à l'exécution

Maintenant, jusqu'à présent, nous avons lancé Nextflow depuis le même répertoire où nous avons notre script de pipeline. Donc notre répertoire de lancement et notre répertoire de pipeline sont en quelque sorte la même chose.

Pour montrer comment nous pouvons avoir différents fichiers de configuration avec différents répertoires de lancement, nous allons créer un nouveau sous-répertoire maintenant.

Je vais donc dire mkdir, et nous allons l'appeler tux-run.

Puis je vais faire cd, changer de répertoire vers tux-run. Et notez que notre répertoire de travail n'est maintenant plus dans le même répertoire que les scripts de pipeline.

D'accord, créons un nouveau fichier "nextflow.config". Donc touch nextflow config, et ouvrons-le simplement dans VS Code. Vous pouvez voir aussi dans la barre latérale ici que nous sommes maintenant dans ce sous-répertoire.

Maintenant, nous pouvons prendre le même bloc params que nous avions dans le nextflow.config de niveau supérieur, le copier ici et maintenant nous pouvons changer ces valeurs.

Premièrement, les données sont maintenant un chemin relatif différent car nous sommes dans un sous-répertoire, nous devons donc mettre à jour cela. Et puis nous allons changer batch en experiment, et nous allons changer le character de Turkey à tux.

Maintenant, cliquons sur sauvegarder là, et essayons. Tout comme avec data, je dois maintenant dire ../ pour accéder au script. Donc c'est Hello config. Et j'appuie sur entrée.

Le code du pipeline n'a pas changé du tout, mais maintenant nous allons avoir deux ensembles de configuration qui se chargent, et le fichier de configuration du répertoire de lancement devrait écraser les valeurs par défaut, qui ont été définies dans le nextflow.config du pipeline, et nous devrions obtenir différents ensembles de résultats.

Effectivement, dans notre répertoire ici, dans tux-run, vous pouvez voir que nous avons un répertoire dot Nextflow et un répertoire work et c'est parce que ceux-ci sont toujours créés dans votre répertoire de lancement. Donc ceux-ci sont différents des répertoires work et results que nous avions des exécutions précédentes.

Maintenant, si je regarde dans results, nous pouvons voir notre collected et voici notre petit personnage tux. Vous pouvez donc voir que ces paramètres ont été correctement interprétés.

## 1.3. Utiliser un fichier de paramètres

D'accord. Avant, quand je parlais des différents fichiers de configuration qui pouvaient être chargés, j'ai omis un autre endroit d'où nous pouvons obtenir la configuration.

Vous pouvez l'obtenir à partir d'une ligne de commande comme nous l'avons vu avec les noms de paramètres précédés de deux tirets, mais nous pouvons également fournir un fichier YAML ou JSON, juste de paramètres.

Le fichier de configuration peut avoir tous les différents types de portées, mais ces fichiers ne contiennent que des paramètres, et c'est une façon conviviale de fournir de nombreux paramètres à la fois, et peut-être une façon un peu plus reproductible car vous les écrivez dans un fichier, il est donc facile de les récupérer ultérieurement.

Retournons donc à notre terminal et juste avant d'oublier, assurons-nous de remonter d'un répertoire, donc je ne suis plus dans le sous-répertoire, et je vais regarder le fichier YAML que nous avons ici appelé test-params.yaml.

Donc si je fais juste code test-params.yaml, vous pouvez voir que c'est juste un fichier YAML ordinaire. Rien de spécial. Avec les clés étant nos noms de paramètres, avec le formatage YAML donc deux-points ici, puis une valeur.

Notez que ce n'est pas du code Nextflow, nous ne pouvons donc pas mettre des choses comme des variables ici. Ce sont juste des valeurs statiques.

Aussi, parce que JSON s'analyse en fait comme YAML, nous pouvons également avoir un fichier test-params.json, qui ressemble beaucoup. C'est juste un format de données différent.

Nous avons donc deux fichiers de test différents ici et nous avons des variables légèrement différentes.

D'accord, alors comment les donnons-nous à Nextflow ? C'est très simple. Nous faisons Nextflow run hello config, comme avant. Et au lieu de "-c" pour fichier de configuration, ou de charger ces noms de fichiers par défaut, nous faisons -params-file. Tiret simple car c'est une option Nextflow de base.

Puis nous passons le chemin de ce fichier. Je vais donc faire "-params-file test-params.yaml", et nous verrons si ceux-ci sont correctement chargés.

D'accord. Ça a fonctionné. Rappelons-nous ce qui était dans ce fichier YAML. Donc le batch était défini sur YAML, c'est donc ainsi qu'il devrait être appelé, et il devrait avoir un stegosaurus. Allons donc voir dans results. Et nous avons COLLECTED-yaml. Voyons donc si nous avons un Stegosaurus. Fantastique, un Stegosaurus portant un chapeau. C'est ce que nous aimons.

Donc cela a très bien fonctionné, et c'est exactement la même chose avec le fichier JSON. Nous changeons juste l'extension du fichier ici et Nextflow sait comment le lire.

Et dans ce cas, nous devrions avoir un batch appelé JSON et nous devrions avoir une tortue. Jetons un coup d'œil. Merveilleux. Un de mes outils CLI préférés.

## 2.1. Personnaliser le répertoire de sortie avec -output-dir

D'accord, donc cela concernait principalement les entrées du pipeline et le changement de paramètres. Qu'en est-il des sorties ?

Maintenant, bien que nous ayons changé les sous-répertoires en utilisant params, vous avez peut-être remarqué que tous nos fichiers vont toujours dans results.

Nous pouvons changer ce répertoire de base où tous les fichiers sont publiés avec un indicateur de ligne de commande appelé -output-dir. Donc si je fais Nextflow run hello config, puis je fais -output-dir, et nous allons l'appeler "custom-outdir-cli". Je ne sais pas taper. Juste pour nous rappeler d'où viennent ces fichiers.

C'est une option Nextflow de base et c'est une option très récente. Cela n'a été ajouté que récemment, et c'est l'une des choses que nous pouvons faire avec le nouveau parseur de langage et tout.

C'est un peu long à taper. Vous pouvez aussi simplement l'appeler "-o" si vous voulez. Donc si je reviens en arrière. Je peux juste raccourcir cela à "-o", ce qui est un peu plus simple.

D'accord. Nous exécutons cela. Nous n'avons rien changé dans notre pipeline ou même dans notre configuration à ce stade, et cela devrait espérons-le sauvegarder tous nos résultats dans un répertoire de niveau supérieur différent. Et vous pouvez imaginer que vous pouvez définir cela sur pratiquement n'importe quel chemin que vous voulez.

Il vient juste d'arriver en haut. Nous avons un custom-outdir-cli, et tous les fichiers y sont organisés exactement de la même manière, avec leurs mêmes sous-répertoires et noms de fichiers. C'est donc un moyen vraiment facile de changer simplement où le pipeline publie ses résultats, sans trop réfléchir à la façon dont ces résultats sont organisés.

## 2.1.2. Supprimer les chemins codés en dur du bloc output

Si je regarde dans ce répertoire, nous pouvons voir que nous avons toujours un sous-répertoire appelé Hello Config, ce qui semble un peu redondant maintenant.

Chargeons donc à nouveau notre script et nous pouvons maintenant supprimer ce sous-répertoire du bloc output en bas. Car nous n'en avons plus vraiment besoin. Nous pouvons donc le faire maintenant, supprimer cela d'ici. Et puis si c'est juste cela, vous pouvez soit supprimer cela complètement, soit le laisser comme une chaîne vide. Je vais le laisser comme une chaîne vide pour l'instant, car nous allons revenir et mettre différentes choses à sa place à l'avenir. Mais si vous ne vous souciez pas des sous-répertoires, il est plus propre de supprimer complètement la déclaration path là.

D'accord, sauvegardons. Essayons rapidement à nouveau. Je vais en fait supprimer mon répertoire "custom-outdir-cli" pour ne pas être confus·e·s par des fichiers existants là. Car rappelez-vous, quand vous publiez des choses, cela ne supprime pas les fichiers qui étaient déjà là. Cela ajoute juste de nouveaux fichiers. Exécutons à nouveau cette commande, custom-outdir-cli.

Et maintenant si vous faites "ls custom-outdir-cli", il n'y a plus de répertoire là appelé Hello Config.

## 2.2.1. Définir outputDir dans le fichier de configuration

D'accord, l'indicateur de ligne de commande ici, "-o" ou "-output-dir" est bien. Mais qu'en est-il de définir des valeurs par défaut pour cela dans la configuration ? Comment faisons-nous cela ?

J'ouvre le fichier "nextflow.config", ferme tout le reste et me débarrasse de cela. Nous pouvons ajouter une nouvelle option de configuration ici, que j'ai juste copiée du site web du matériel de formation, et elle s'appelle outputDir.

Elle n'est sous aucune portée. Elle n'est pas sous params ou quoi que ce soit. Elle est au niveau supérieur, et nous pouvons la définir sur une chaîne. Maintenant, une chose simple à faire est juste de la changer en n'importe quoi d'autre que results comme une chaîne codée en dur. Mais parce que c'est dans un fichier de configuration Nextflow, nous pouvons être un peu malins ici et aussi inclure des variables.

Et vous pouvez voir ici que nous avons inclus une variable params, params.batch, qui fait partie de cette chaîne. Cela signifie que nous pouvons réutiliser des variables qui viennent d'autres endroits. Et dans ce cas, si nous faisons --batch, quand nous exécutons le pipeline Nextflow, nous allons obtenir un sous-répertoire dans notre chemin personnalisé basé sur le nom du batch.

D'accord, essayons cela et jetons juste un coup d'œil rapide pour voir comment les résultats se présentent. Donc si je fais Nextflow run hello config et --batch my_run. Rappelons-nous à quoi ressemblait la configuration. Donc c'est custom-outdir-config.

Tree custom-outdir-config. Et vous pouvez voir que le batch s'appelait my_run. Et puis nous avons ce sous-répertoire appelé my_run. Donc ce chemin de fichier dynamique a fonctionné.

Et non seulement cela, il n'est plus allé dans un répertoire results par défaut, et je n'ai pas eu à spécifier quoi que ce soit sur la ligne de commande pour changer le répertoire de base. Nous avons donc réussi à réinitialiser la valeur par défaut pour le outputDir par défaut.

## 2.2.2. Sous-répertoires avec les noms de batch et de processus

D'accord, allons un peu plus loin. C'est une variable dynamique dans le fichier de configuration. Qu'en est-il du script lui-même ? Maintenant, jusqu'à présent, nous avions ces chemins ici et ceux-ci peuvent aussi être dynamiques. Donc au lieu de simplement coder en dur quelque chose, nous pouvons mettre des accolades et mettre quelque chose de dynamique.

Par exemple, nous avons nos processus appelés sayHello. Nous pourrions faire sayHello.name, qui est un attribut du processus, ce qui est un peu ennuyeux car c'est juste "sayHello" dans ce cas. Mais c'est variable.

Cela vous donne donc une idée. Nous pouvons donc mettre cela ici et dire convertToUpper.name, collectGreetings.name, collectGreetings.name à nouveau, et cowpy.

Maintenant, quand nous exécutons, le répertoire de base va toujours être custom-outdir-config. Et il va être dans un sous-répertoire appelé params.batch, mais les sous-répertoires en dessous devraient être organisés par nom de processus.

Essayons simplement cela et voyons si cela fonctionne. Je vais donc supprimer le répertoire précédent pour ne pas être confus·e·s, et utiliser exactement la même commande Nextflow Run.

Cela devrait s'exécuter de la même manière. Je pourrais utiliser dash resume sur tous ceux-ci pour que ce soit un peu plus rapide et utiliser les résultats précédemment calculés. Maintenant, si je fais tree custom-outdir-config, vous pouvez voir que ce n'est pas dans results, c'est dans notre répertoire de base avec le nom du batch. Et vous pouvez voir que tous les résultats sont maintenant organisés dans des sous-répertoires nommés d'après le processus. Nous avons donc deux endroits différents où nous définissons des chemins de sortie dynamiques ici.

D'accord. Dernière chose, rajoutons ces dossiers intermédiaires, que nous avions avant car ils étaient plutôt sympas. Intermediates.

Et nous pouvons aussi réfléchir un peu à ce params.batch, peut-être qu'en tant que développeur·se de pipeline, j'aimais vraiment avoir cela dans le sous-répertoire, mais si les utilisateur·trices finaux·ales du pipeline définissent "-o" ou -output-dir sur la CLI, cela écrase complètement cette déclaration entière, et nous perdons ce sous-répertoire.

Ce que nous pouvons faire, c'est prendre ce chemin dynamique hors de la configuration outputDir, qui serait écrasée, et le mettre dans le chemin de sortie, qui n'est pas écrasé.

Nous pouvons donc faire params.batch slash intermediates slash sayHello.name, et faire tout cela dans une chaîne entre guillemets doubles, pour qu'elle soit interpolée par Nextflow.

Je peux maintenant copier, oups. Copier ceux-ci vers les autres processus. N'oubliez pas de tous les mettre entre guillemets. Et supprimer intermediates de ces sorties particulières.

D'accord ? Cela semble légèrement plus complexe maintenant, mais vous pouvez voir que nous commençons vraiment à construire une belle structure de répertoire de sortie organisée dans notre code.

Et ce qui est vraiment bien, c'est que cette complexité supplémentaire dans le code ne passe pas à la CLI. Nous pouvons donc exécuter notre commande avec -output-dir et quelles que soient les variables batch, en pensant juste à comment exécuter le pipeline et sans vraiment trop penser à ce qui est dans le code. Et nos fichiers de sortie vont être construits de manière vraiment agréable d'une manière très bien organisée, ce qui est agréable pour les personnes utilisant le pipeline en gros.

Parfait. En écrivant ceci, je réalise que j'ai fait une erreur. Voyons si quelqu'un m'a attrapé ici. Nous avons collectGreetings.name, donc quelque chose s'est mal passé. Et oui, effectivement, j'ai accidentellement oublié de mettre ceux-ci entre accolades.

Rappelez-vous donc, soyez prudent·e·s quand vous écrivez votre code et assurez-vous de dire à Nextflow ce qui est une variable et ce qui est juste une chaîne. Car il fera exactement ce que vous lui dites de faire. Et rien de plus. Comme tous les bons ordinateurs. D'accord, cela devrait le corriger.

## 2.3. Définir le mode de publication au niveau du workflow

Il y a une partie de ce script que je n'aime toujours pas, c'est le fait que nous écrivons mode copy encore et encore, et s'il y a une chose que nous n'aimons pas, c'est nous répéter.

Nous pouvons donc nettoyer cela un peu en prenant ceci et en le déplaçant dans la configuration. Et en fait, nous pouvons le définir pour l'ensemble du pipeline en une seule fois. Nous n'avons donc pas à le dire plusieurs fois.

Nous allons dans notre fichier de configuration et nous avons une nouvelle portée ici appelée workflow. Et nous pouvons soit faire des accolades, soit utiliser la notation par points. Cela ne fait aucune différence. Le site web du matériel de formation utilise la notation par points. Je peux dire output et nous pouvons mélanger et assortir, donc mode equals copy. Parfait.

Et maintenant nous pouvons revenir ici et supprimer ceux-ci. Nous pourrions les laisser en place. La configuration écrase essentiellement ce qui est écrit ici, mais comme nous l'avons dans la configuration au niveau du pipeline, et que ces deux fichiers sont livrés ensemble, il n'y a vraiment aucune raison de le faire deux fois.

D'accord. Juste une vérification de bon sens, car apparemment nous faisons des erreurs. Exécutons cela à nouveau et vérifions simplement que nous utilisons correctement le mode copy pour publier les fichiers. Nous allons donc exécuter le script à nouveau et cette fois nous avons mis les résultats dans un répertoire appelé config-output-mode, voyons à quoi ressemblent les fichiers là-dedans.

Et puis si je fais "ls -l" pour regarder batch, et nous pouvons regarder cowpy, par exemple. Et nous devrions voir, oui, que c'est un fichier approprié ici, qui n'est pas un lien symbolique, donc cet attribut de configuration a été appliqué correctement.

## 3. Sélectionner une technologie d'empaquetage logiciel

D'accord. Jusqu'à présent, nous nous sommes concentré·e·s sur les entrées et les sorties, les fichiers avec lesquels le workflow s'exécute. Mais qu'en est-il de l'infrastructure ? J'ai dit au début que Nextflow vous permet d'exécuter le même pipeline sur différentes configurations informatiques. Alors à quoi cela ressemble-t-il ?

Pour montrer cela, nous allons passer de l'utilisation de Docker pour exécuter cowpy, et à la place nous utiliserons Conda pour faire la même chose.

Je peux faire cela très simplement. Si je vais dans code, "nextflow.config". Si vous vous souvenez en haut, nous avons défini docker.enabled plus tôt, et dans le dernier chapitre pour pouvoir utiliser le conteneur avec cowpy dedans.

Je vais dire à Nextflow de ne pas utiliser Docker. Définir cela sur false. Et je vais dire Conda enabled equals true. Donc dire à Nextflow, s'il vous plaît utilisez Conda.

Maintenant, activer simplement Conda n'est pas suffisant en soi. Tout comme nous l'avons fait avec Docker, nous devons dire à Nextflow où il peut obtenir le logiciel dont il a besoin.

Donc si nous allons dans les modules ici. Et ouvrons le script cowpy. Nous pouvons voir que nous avons une déclaration container en haut. Et le conteneur est utilisé par Docker, mais aussi Singularity, Apptainer, et beaucoup d'autres outils logiciels.

Mais il ne peut pas être utilisé pour Conda, nous avons donc une déclaration séparée appelée "conda", et nous pourrions simplement écrire "cowpy". Et cela laissera à la résolution de paquets conda le soin de déterminer la meilleure façon de résoudre cela, selon votre environnement conda local.

Ou c'est une bonne pratique de faire ce que le site web du matériel de formation dit de faire, qui est de définir un canal conda spécifique avec sa notation à double deux-points, et définitivement définir une version spécifique du logiciel pour que chaque personne qui exécute le pipeline obtienne la même version.

Notez que les conteneurs sont un peu supérieurs à cet égard, car lorsque vous installez quelque chose avec Conda, il va toujours déterminer toutes les dépendances pour ce paquet, et elles peuvent changer au fil du temps. C'est ce qu'on appelle la dérive des dépendances.

Les conteneurs, cependant, verrouillent toute la pile de dépendances logicielles jusqu'en bas, vous pouvez donc être un peu plus confiant·e que A, cela va fonctionner, et B, ce sera reproductible.

Donc si vous êtes en mesure d'utiliser Docker ou Singularity ou Apptainer, je recommanderais définitivement cela.

Maintenant, ce qui est bien avec ceci, c'est que le fichier de module, qui est écrit par le développeur·se du pipeline, a maintenant à la fois Container et Conda, et nous disons donc à la personne qui exécute ce pipeline, peu importe quelle solution d'empaquetage logiciel vous utilisez. Cela fonctionnera à la fois avec Docker et avec Conda, et voici où obtenir le logiciel dans les deux cas.

Nous pouvons ouvrir le terminal et essayons cela. Donc Nextflow run hello config --batch conda. Et la première fois que cela s'exécute avec conda, cela va être un peu lent quand il arrive à ce processus particulier, car il doit exécuter "conda install".

Et il crée un environnement conda spécial juste pour ce processus. Il n'utilise donc pas mon environnement conda global, que j'ai sur mon terminal. Il en crée un juste pour ce processus. C'est bien car cela évite des choses comme les conflits de dépendances entre différents processus dans votre workflow. Si vos processus ont des outils qui ont besoin de différentes versions de Python ou des choses comme ça, ce n'est pas grave car ils utilisent différents environnements conda.

Nextflow met en cache ces environnements conda localement, vous pouvez voir qu'il vous dit où se trouve ce chemin, c'est dans le répertoire work ici. Et donc la prochaine fois que j'exécute ce script avec Conda, ce sera beaucoup plus rapide car il trouvera cet environnement conda existant et le réutilisera simplement. Mais la première fois que nous le faisons, il doit aller le chercher, le résoudre, télécharger toutes les dépendances et tout configurer.

D'accord, super, ça a fonctionné. Nous pouvons juste nous rappeler ce que le pipeline est actuellement configuré pour utiliser. Si nous regardons dans le fichier de configuration, c'était "custom-outdir-config" en ce moment pour moi. Voyons si je vais dans ce répertoire de base. Et j'ai fait --batch conda. Voici notre sous-répertoire conda. Donc cela a fonctionné et voici notre sortie cowpy.

Il a donc récupéré cowpy, l'a installé sur mon système local en utilisant conda, et a exécuté le processus. Et ce qui est génial, c'est qu'en tant qu'utilisateur·trice final·e, je n'ai pas eu à penser du tout à la gestion des logiciels là. Nextflow l'a juste arrangé pour moi. J'ai dit, j'ai besoin d'utiliser conda sur ce système. Le développeur·se du pipeline a dit quels paquets j'avais besoin. Et Nextflow a fait le reste. Très puissant.

Notez que vous pouvez en fait utiliser un mélange de différentes technologies. Je peux donc activer Docker pour des processus spécifiques, et conda pour d'autres processus, ou dire que certains processus devraient juste utiliser quel que soit le logiciel local que j'avais installé. C'est assez inhabituel, mais c'est possible, et dans certains cas, par exemple, si vous utilisez certains logiciels qui pourraient être difficiles à empaqueter dans Docker, vous avez une échappatoire.

## 4. Sélectionner une plateforme d'exécution

Donc c'est l'empaquetage logiciel. L'autre partie de la portabilité vers d'autres systèmes est l'endroit où les tâches s'exécutent réellement. Pour le moment, je m'exécute sur essentiellement mon ordinateur portable ou dans ce Codespaces, qui est un seul ordinateur. Il n'y a rien de sophistiqué. Nextflow est un peu malin pour paralléliser les tâches du mieux qu'il peut, mais tout est sur un seul système.

Maintenant, si vous exécutez sur un HPC, vous avez probablement une sorte de planificateur de tâches tel que SLURM ou PBS ou quelque chose, et vous soumettez des tâches à ce planificateur et il distribue toutes les tâches à différents nœuds de calcul.

Une autre façon d'exécuter est sur le cloud. Peut-être que vous utilisez AWS Batch, ou Azure Cloud, ou Google. Et tous ceux-ci fonctionnent dans un système similaire où vous avez un planificateur et vous soumettez des tâches et elles sont soumises à différents endroits pour être calculées.

Maintenant, dans un passé lointain quand j'ai commencé à faire de la bioinformatique, tous les logiciels de tout le monde pour exécuter des analyses étaient très liés à leur infrastructure informatique, ce qui rendait presque impossible la réplication.

Mais avec cette séparation de configuration dans Nextflow, et avec la capacité de Nextflow à interagir avec de très nombreux backends d'infrastructure de calcul différents, il est très simple de prendre notre pipeline sans modifier du tout le code du pipeline et de simplement changer cela.

## 4.1. Cibler un backend différent

Donc si nous allons dans notre fichier "nextflow.config", et nous pouvons maintenant mettre une configuration au niveau du processus. Donc si je mets en haut la portée process et je peux définir l'executor, et ici il est défini sur local, qui est la valeur par défaut.

Notez que parce que c'est au niveau du processus, nous pouvons cibler des choses vers différents processus. Et vous pouvez donc en fait configurer des executors pour être spécifiques au processus et avoir une exécution hybride, où certaines tâches pourraient s'exécuter localement, où que la tâche Nextflow soit exécutée. Certaines sont soumises à différents HPC et certaines pourraient être soumises au cloud. Vous pouvez être aussi malin·e que vous le souhaitez.

Maintenant, il est très difficile de démontrer cela dans un environnement de formation comme celui-ci car je n'ai pas de HPC auquel soumettre. Mais ce que je peux faire, c'est si je tape slurm, nous pouvons tricher un peu et vous pouvez avoir une idée de cela.

Et c'est vraiment plus intéressant pour les personnes qui ont l'habitude d'exécuter sur SLURM et connaissent à quoi ressemblent les en-têtes SLURM. Mais si je fais Nextflow run, hello config. Cela va échouer car cela va essayer de soumettre des tâches à un cluster qui n'existe pas. Nous obtiendrons donc une sorte d'erreur à propos de sbatch qui n'est pas disponible.

Oui, écrit. C'est l'outil. C'est l'outil CLI que vous utilisez pour soumettre des tâches à un cluster slurm. Mais ce que nous pouvons faire, c'est aller regarder dans notre répertoire work ici en cliquant avec commande, ouvrir ce répertoire et regarder le .command.run. Et vous pouvez voir en haut du fichier .command.run, nous avons nos en-têtes sbatch, disant à un cluster SLURM théorique comment gérer cette soumission de tâche.

Et vous pouvez donc voir que Nextflow est malin, il fait toutes les bonnes choses. C'est juste que nous n'avions pas de cluster auquel soumettre.

## 5. Contrôler les allocations de ressources de calcul

Qu'est-ce qui est différent d'autre entre différentes infrastructures informatiques ? Une autre chose est la quantité de ressources disponibles que vous avez, et en fait, dans de nombreux environnements de calcul, c'est une exigence que vous devez spécifier combien de CPUs et combien de mémoire une tâche a besoin.

Encore une fois, Nextflow abstrait cela pour nous, de sorte que ce n'est plus spécifique à un seul type d'environnement de calcul, et nous pouvons taper dans la portée au niveau du processus ici. CPUs equals one, memory equals two gigabytes. Notre pipeline n'est pas très exigeant, donc cela devrait aller.

Maintenant, j'ai juste deviné ces nombres ici, mais comment savez-vous quelle est une quantité sensée de ressources à utiliser ? C'est un travail assez difficile d'aller fouiller dans tous ces différents processus d'un grand pipeline de nombreux échantillons et de comprendre quelle était l'utilisation des ressources.

Une bonne approche pour cela est donc de définir ces valeurs à des nombres élevés pour commencer, juste pour que votre pipeline s'exécute sans erreurs, puis de demander à Nextflow de générer un rapport d'utilisation pour vous.

C'est super facile à faire, donc je vais retourner à un terminal. Oh, je dois me rappeler de remettre cela sur local pour que mon pipeline s'exécute réellement. Et je vais dire Nextflow run, et je vais utiliser un indicateur de ligne de commande -with-report.

Et je peux laisser cela vide et il donnera un nom de fichier par défaut, mais je vais lui donner un nom de fichier spécifique pour que cela soit sauvegardé à un endroit spécifique.

Appuyez sur Entrée, et le pipeline s'exécute exactement comme d'habitude, mais quand il se termine, il va générer un joli rapport HTML pour moi.

Donc dans la barre latérale ici, j'ai ce fichier HTML. Si j'exécutais cela localement, je l'ouvrirais simplement. Comme je suis dans Codespaces, je vais faire un clic droit dessus et cliquer sur télécharger, ce qui va le télécharger sur mon ordinateur local. Et je peux juste facilement l'ouvrir dans le navigateur web.

Nextflow peut générer un rapport comme celui-ci pour n'importe quel pipeline et il contient des informations vraiment intéressantes. C'est donc une bonne pratique de toujours sauvegarder ces choses. Il nous dit quand nous avons exécuté, où nous avons exécuté, si c'était réussi ou non, quels paramètres ont été utilisés, quelle était la commande CLI, des choses comme ça.

Et il y a aussi ces graphiques sur l'utilisation des ressources. Il nous dit donc quel pourcentage d'appels CPU ont été utilisés pour chaque processus sous forme de boîte à moustaches ici, car il y a de nombreuses tâches pour chaque processus, nous pouvons donc voir la distribution.

Vous pouvez voir nos processus ici, cowpy et collectGreetings n'avaient qu'une seule tâche, donc c'est juste une seule ligne. Et nous avons à la fois CPU et mémoire et durée de la tâche, et ils étaient très rapides.

Si vous utilisez Seqera Platform, soit dit en passant, vous obtenez les mêmes graphiques intégrés dans l'interface Platform sans avoir à faire quoi que ce soit. Vous obtenez donc toujours ces informations à portée de main.

D'accord, nous pouvons donc utiliser ce rapport et sur une vraie exécution, et avoir une idée de combien de CPUs et combien de mémoire sont utilisés par notre pipeline et revenir et remettre ces valeurs dans notre fichier de configuration, pour que la prochaine fois peut-être nous ne demandions pas autant. Et nous pouvons être un peu plus économes.

Maintenant, vous pouvez devenir vraiment malin·e dans la configuration des fichiers de configuration de pipeline. Et encore une fois, si vous utilisez Seqera Platform, recherchez un petit bouton qui ressemble à une ampoule. Car si vous cliquez dessus, il générera un fichier de configuration hautement optimisé, qui est adapté spécifiquement à vos données, votre exécution et votre pipeline. Pour l'exécuter de la manière la plus efficace possible.

Mais pour l'instant, je vais dire qu'en fait le nombre par défaut de CPUs que Nextflow donnait était bien et nous n'avons besoin que d'un gigaoctet de mémoire.

## 5.3. Définir les allocations de ressources pour un processus spécifique

Maintenant, dans la vraie vie, il est assez inhabituel que tous les processus de votre pipeline aient besoin des mêmes exigences. Vous pourriez avoir quelque chose comme MultiQC comme outil de rapport, qui a besoin de très peu en termes de ressources et s'exécute assez rapidement.

Et puis peut-être que vous avez quelque chose qui indexe un génome de référence ou fait un alignement ou fait un autre travail. Peu importe ce que c'est, qui prend beaucoup de ressources. Et donc pour ces différentes soumissions de tâches à un planificateur, vous voulez donner différentes quantités de ressources.

Sous cette portée process, nous pouvons définir une configuration, qui cible des processus spécifiques de différentes manières.

Ici, nous utilisons withName, nous pouvons aussi utiliser des labels, et ceux-ci peuvent utiliser un motif pour cibler un ou plusieurs processus. Ici, nous disons simplement que tous les processus qui ont un nom cowpy sont définis sur deux gigaoctets de mémoire et deux CPUs, et parce que c'est un sélecteur plus spécifique que celui au niveau supérieur du processus, celui-ci est écrasé dans ces cas, vous pouvez donc construire un joli fichier de configuration ici, qui adapte vraiment tous vos différents processus dans votre pipeline pour les rendre vraiment efficaces.

## 5.5. Ajouter des limites de ressources

Maintenant, en tant que développeur·se de pipeline, je connais probablement assez bien les outils, et je veux que tout s'exécute aussi rapidement et aussi bien que possible. Il se peut donc que je mette des nombres assez élevés pour certains d'entre eux car je sais que cela s'exécutera beaucoup plus rapidement si je donne à cowpy 20 CPUs au lieu de deux.

C'est bien jusqu'à ce que vous alliez exécuter sur votre ordinateur portable ou sur GitHub Actions Continuous Integration test, ou un autre système, qui n'a peut-être pas 20 CPUs disponibles.

Maintenant, quand vous essayez d'exécuter le pipeline, il va planter car Nextflow dira, je ne peux pas soumettre cette tâche nulle part. Je n'ai pas les ressources disponibles.

Maintenant, pour éviter ce plantage brutal, nous pouvons ajouter un peu plus de configuration, qui est spécifique à notre système maintenant, appelée limites de ressources. Et cela ressemble à ceci. C'est sous la portée process à nouveau.

Et limites de ressources, vous pouvez spécifier essentiellement le plafond de ce que vous avez disponible. C'est une map ici, et vous pouvez, dans cette map, vous pouvez définir la mémoire, les CPUs et le temps.

Maintenant, ce qui se passe, c'est que lorsque Nextflow soumet une tâche à partir d'un processus, il regarde ce qui est demandé et il fait essentiellement un minimum entre cela et cela. Donc si nous avons demandé 20 CPUs, mais seulement quatre sont disponibles, il demandera quatre. Le pipeline ne plante pas et il utilise aussi proche que possible de ce qui a été conçu par le développeur·se du pipeline.

## 6. Utiliser des profils pour basculer entre des configurations prédéfinies

D'accord. J'ai dit que les limites de ressources ici pourraient être spécifiques au système, et peut-être que j'ai un fichier Nextflow config dans mon pipeline, et je sais que les gens vont utiliser cela dans une gamme d'endroits différents. Maintenant, au lieu de forcer tout le monde à créer son propre fichier Nextflow config à chaque fois, ce que je peux faire, c'est regrouper différents préréglages de configuration ensemble dans des profils de configuration.

Je vais descendre un peu ici et juste après params, car l'ordre du fichier de configuration ici est important, le fichier de configuration est chargé séquentiellement, je vais donc mettre ces profils après tout le reste pour qu'ils écrasent les params précédemment définis. Et je vais coller ces profils du matériel de formation.

Il y a donc une nouvelle portée de niveau supérieur appelée profiles. Nous pouvons avoir des noms arbitraires ici. Nous avons donc my_laptop et univ_hpc. Et ici nous pouvons voir que nous définissons les mêmes paramètres de configuration que nous avions avant. Maintenant juste dans un profil. Nous avons donc un executor local pour exécuter sur my_laptop et je soumets à un cluster SLURM sur le HPC.

J'utilise Docker localement, conda sur le HPC, et le système HPC a des limites de ressources beaucoup plus élevées.

Maintenant, je peux exécuter le pipeline avec l'option CLI -profile, dire quel profil je veux utiliser. Je vais donc utiliser my_laptop, et Nextflow appliquera toute la configuration dans cette portée de profil. Je peux donc essayer cela maintenant. C'est la même commande qu'avant. Nextflow run hello config, et je fais dash profile, tiret simple car c'est l'option Nextflow de base, dash profile my_laptop.

Il va maintenant appliquer par lot toute cette option de configuration. Oh, et vous pouvez voir, j'ai dit avant que cela pourrait arriver que l'exigence du processus, il a demandé quatre CPUs et je n'en ai que deux sur cette instance Codespaces.

C'est donc une bonne opportunité juste pour essayer les limites de ressources du processus, et dire que je n'ai que deux CPUs sur my_laptop, ou dans ce Codespaces. Maintenant, si nous l'exécutons à nouveau, il devrait plafonner cette exigence à deux et espérons-le le pipeline s'exécutera. Parfait.

## 6.2. Créer un profil de paramètres de test

Notez que ces profils ne doivent pas seulement avoir une configuration sur leur infrastructure. Vous pouvez avoir des regroupements de n'importe quelle configuration ici, y compris des paramètres.

Une autre chose que vous verrez très souvent dans les pipelines des gens est un profil test, qui inclut des paramètres, que vous soumettriez normalement par utilisateur·trice. Mais ici nous avons, essentiellement différentes valeurs par défaut sensées pour quand je veux exécuter des cas de test.

Et c'est génial car je n'ai pas nécessairement à aller spécifier toutes ces choses, qui pourraient être des paramètres requis. Sinon, je peux juste dire dash profile test et cela s'exécutera juste prêt à l'emploi.

Maintenant, quelque chose à noter est que les profils peuvent aussi être combinés plus d'un. Je peux donc faire profile my_laptop ici, puis aussi ajouter test. Je ne fais pas profile deux fois. Je fais juste une liste séparée par des virgules ici sans espaces. Et il va appliquer ces profils dans l'ordre. Il prendra donc la configuration du profil my_laptop, puis il appliquera la configuration test par-dessus.

Vraiment pratique et vous pouvez voir comment vous pouvez configurer beaucoup de groupes par défaut sensés ici pour faciliter l'exécution de votre pipeline.

## 6.3. Utiliser nextflow config pour voir la configuration résolue

J'espère vous avoir convaincu·e·s que la résolution de configuration Nextflow est puissante, mais je ne vous en voudrais pas si vous avez un peu les yeux qui tournent à ce stade après que j'ai dit environ 20 façons différentes de fournir la configuration et donné toutes ces différentes couches comme une pelure d'oignon.

Donc si jamais vous vous sentez incertain·e·s sur ce qu'est la configuration finale résolue pour Nextflow, sachez qu'il y a une commande appelée "nextflow config", et nous pouvons l'exécuter et elle nous dira quelle est la configuration résolue à notre emplacement actuel.

Donc quand je l'exécute ici, il trouve le fichier "nextflow.config" dans le répertoire de travail actuel, et il traite toute la configuration différente, et il me donne la sortie résolue.

Notez que le fichier de configuration Nextflow peut aussi prendre l'option CLI profile. Donc si je lui dis de résoudre dans les profils my_laptop et test, et vous pouvez voir qu'il a aussi appliqué les limites de ressources ici de l'option de configuration my_laptop et aussi défini les params, qui étaient dans le test.

C'est donc une bonne façon juste d'explorer comment fonctionne la résolution de configuration, si vous avez le moindre doute.

## Conclusion

D'accord, c'est tout. C'est la configuration Nextflow en bref. Vous pouvez faire beaucoup de choses avec la configuration. C'est vraiment puissant. Mais ce sont la plupart des cas d'usage courants que vous vous trouverez à faire, et ces concepts s'appliquent à toutes les différentes options.

Félicitez-vous car c'est la fin du cours de formation Hello Nextflow. Vous êtes espérons-le maintenant confiant·e·s à la fois pour écrire votre propre pipeline Nextflow à partir de zéro, le configurer et l'exécuter, et vous connaissez tous les tenants et aboutissants et les choses à surveiller.

Il y a un dernier quiz que vous pouvez essayer sur la page de formation sur la configuration. Descendez donc et essayez-le et assurez-vous d'avoir compris toutes ces parties sur la configuration.

Et rejoignez-nous dans la dernière vidéo juste pour une conclusion rapide sur certaines des prochaines étapes qui pourraient être bonnes à faire après ce cours de formation.

Merci d'être resté·e·s avec nous. Bravo et je vous verrai dans la prochaine vidéo.
