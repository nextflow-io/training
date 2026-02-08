# Prochaines étapes - Transcription vidéo

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Traduction assistée par IA - [en savoir plus et suggérer des améliorations](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

<div class="video-wrapper">
  <iframe width="560" height="315" src="https://www.youtube.com/embed/wnmUAfDL01E?si=Lp8hS8RdaMwbp5j5&amp;list=PLPZ8WHdZGxmWKozQuzr27jyMGqp9kElVK&amp;cc_load_policy=1" title="YouTube video player" frameborder="0" allow="accelerometer; autoplay; clipboard-write; encrypted-media; gyroscope; picture-in-picture; web-share" referrerpolicy="strict-origin-when-cross-origin" allowfullscreen></iframe>
</div>

!!!note "Notes importantes"

    Cette page affiche uniquement la transcription. Pour les instructions complètes étape par étape, retournez au [contenu du cours](../next_steps.md).

## Bienvenue

​

Félicitations, vous avez réussi !

Vous êtes arrivé·e au bout et avez terminé le cours de formation Hello Nextflow. Nous espérons sincèrement que vous l'avez apprécié. Merci d'être resté·e avec nous tout au long du parcours, et nous apprécions vraiment le temps et les efforts que vous avez consacrés à l'apprentissage de Nextflow. Nous espérons vraiment que cela sera utile pour votre travail.

## Autres cours sur training.nextflow.io

N'oubliez pas de revenir régulièrement sur training.nextflow.io. Nous ajoutons de nouveaux cours courts tout le temps et nous actualisons également une grande partie du matériel qui est déjà présent. Ce cours de formation Hello Nextflow sera donc mis à jour au fil du temps.

Ceci est particulièrement important car nous mettons à jour la syntaxe dans Nextflow, et 2026 verra l'arrivée de nombreuses nouvelles fonctionnalités, donc ce cours aura un aspect et une sensation un peu différents la prochaine fois que nous le ferons en 2027.

Plus précisément, je tiens à mentionner la page « Nextflow for Science ». Ce sont des cours courts et ils sont conçus pour faire suite à ce cours Hello Nextflow. Ils montrent comment utiliser Nextflow avec différents cas d'usage spécifiques, qu'il s'agisse de génomique, de RNAseq ou de toutes sortes de choses différentes. Nous essayons d'ajouter constamment davantage de cas d'usage scientifiques.

Il y a aussi les Quêtes secondaires. Lorsque nous développons un cours comme Hello Nextflow, il y a tellement de sujets que nous pourrions couvrir qu'il est difficile de garder tout dans le périmètre. Donc, si un sujet particulier nous semble intéressant pour les participant·es et mérite d'être approfondi, nous le plaçons dans une Quête secondaire.

Allez jeter un œil et si certains sujets pourraient être pertinents pour votre travail, comme nf-test ou faire différentes choses avec les métadonnées, et les modèles de script courants, consultez les Quêtes secondaires pour voir si cela pourrait être utile d'en apprendre davantage.

Il y a aussi le cours sur nf-core. Espérons que vous connaissez le projet à ce stade, mais si ce n'est pas le cas, allez le découvrir. Il existe près de 150 pipelines différents pour différents types d'analyse et différents types de données, il est donc tout à fait possible qu'il existe un pipeline prêt à l'emploi pour le type d'analyse de données dont vous avez besoin.

Il est important de noter qu'il existe également des composants dans nf-core, près de 1700 modules différents, différents processus et encapsuleurs pour des outils. Et avec les outils fournis avec nf-core, vous pouvez les combiner et construire votre propre pipeline comme des briques Lego. Bien plus rapidement et de manière plus reproductible.

## Seqera Platform

À mesure que vous intensifiez votre utilisation de Nextflow, découvrez Seqera Platform, c'est la meilleure façon d'exécuter Nextflow. Vous pouvez l'exécuter sur votre propre infrastructure, que ce soit HPC ou AWS, Azure, Google Cloud, Oracle et plus encore. Vous pouvez également utiliser notre propre Seqera Compute si vous ne souhaitez gérer aucune infrastructure de calcul.

Seqera Platform simplifie vraiment la configuration de ces infrastructures cloud complexes avec des fonctionnalités comme Batch Forge, qui crée l'environnement pour vous. Et cela aide également vraiment avec l'observabilité, la journalisation d'audit et la conformité.

Cela permet objectivement de faire fonctionner les pipelines de manière moins coûteuse et plus rapide avec des technologies comme Fusion, qui optimisent l'accès disque et les transferts de données. Et il y a aussi l'optimisation de pipeline pour s'assurer que la configuration de vos pipelines est aussi finement réglée que possible.

Il existe des fonctionnalités totalement différentes en plus de l'exécution de pipelines. Nous avons Studios où vous pouvez exécuter des analyses interactives et créer des environnements à partir de n'importe quelle image docker personnalisée que vous créez. Et Data Explorer, qui vous aide à explorer vos différents systèmes de fichiers où qu'ils se trouvent.

Il existe un niveau gratuit pour Seqera Platform, vous pouvez donc utiliser pratiquement toutes ces fonctionnalités gratuitement dès maintenant. Et nous vous offrirons même cent dollars de crédit de calcul gratuit avec Seqera Compute si vous vous inscrivez avec votre adresse e-mail professionnelle. Enfin, il existe un programme académique, donc si vous travaillez dans une université, consultez la page de tarification, trouvez le formulaire là-bas et faites-le nous savoir, et nous vous mettrons à niveau vers Cloud Pro gratuitement.

## Aide communautaire et événements

Très bien. Pour aller de l'avant. Si vous avez besoin d'aide avec Nextflow, consultez community.seqera.io. C'est vraiment actif et nous espérons vous y voir pour discuter de vos différents problèmes et cas d'usage, et peut-être maintenant vous pourrez même aider d'autres personnes.

Nous organisons également de nombreux événements. Nous avons des événements communautaires provenant de nf-core et Nextflow. Nous avons un hackathon nf-core en ligne et distribué en mars, nous avons eu plus d'un millier de participant·es l'année dernière avec des sites partout dans le monde. Rejoignez-nous donc si vous le pouvez.

Et nous avons également des événements Nextflow Summit, un à Boston, puis un événement à Barcelone et en ligne. Des présentations fantastiques où vous pouvez entendre parler de personnes utilisant Nextflow de manières vraiment massives, folles et passionnantes. Et il y a aussi des hackathons associés à ceux-ci et des formations en présentiel.

## Podcast et blog Nextflow

Si vous souhaitez rester au courant de ce qui se passe dans l'écosystème Nextflow, consultez seqera.io/blog.

Il existe une section pour Nextflow où vous pouvez lire des articles de blog communautaires de personnes travaillant dans la communauté, ainsi que des articles de blog de Seqera sur les mises à jour de Nextflow et des autres outils que nous générons.

J'aimerais également faire un peu de promotion pour mon projet personnel, qui est le Nextflow Podcast. Écoutez-le sur Spotify, Apple Music ou YouTube. Nous publions de nouveaux épisodes périodiquement où je discute avec d'autres personnes, soit travaillant avec Nextflow ou des technologies associées, soit des personnes de la communauté. Et nous faisons de véritables plongées techniques approfondies sur le fonctionnement des choses et ce que les gens font. Donc si vous êtes intéressé·e, écoutez-les. Ils sont vraiment amusants.

## Remerciements

D'accord, j'aimerais faire une série de remerciements. L'équipe de formation de Seqera est responsable de ce matériel. Je suis assis devant une caméra, mais en réalité, tout le travail difficile a été fait par ces autres personnes. Mention spéciale pour Geraldine, qui a écrit et actualisé ce matériel de formation pour le cours Hello Nextflow et d'autres. Et aussi Jon, qui a vraiment aidé, en particulier avec la mise à jour de la syntaxe pour la nouvelle syntaxe Nextflow et également en écrivant de nombreux cours lui-même. D'autres membres de l'équipe de développement scientifique tels que Rike, Rob, Florian et beaucoup d'autres ont eu une contribution énorme au matériel avec lequel nous avons travaillé.

J'aimerais également remercier les personnes de la communauté. Les nouvelles traductions, par exemple, qui sont très récentes, ont été fortement influencées par les membres du programme d'ambassadeur·trices et ailleurs. Et vraiment, la nature open source du matériel de formation signifie que nous recevons des pull requests et des issues assez fréquemment, ce qui nous aide vraiment.

## Enquête

Maintenant que vous avez terminé, si vous ne l'avez pas encore fait, veuillez rapidement remplir l'enquête de retour d'expérience. Elle se trouve sur le site web training.nextflow.io juste en dessous de la section Hello Nextflow.

Il n'y a que cinq questions. C'est vraiment, vraiment rapide, mais cela nous permet de suivre à peu près combien de personnes suivent la formation et vous pouvez également nous dire comment améliorer le matériel de formation. Nous vérifions vraiment toutes les réponses, donc nous apprécions vraiment vos commentaires.

## Au revoir

Encore une fois, merci beaucoup de nous avoir rejoint pour ce cours et pour ce voyage. Déposez une issue ou une Pull Request GitHub si vous avez repéré quelque chose dans le matériel de formation qui pourrait être amélioré selon vous. Et j'espère vraiment vous voir dans un autre cours de formation Nextflow, ou peut-être lors d'un hackathon ou d'un événement. Merci encore.​
