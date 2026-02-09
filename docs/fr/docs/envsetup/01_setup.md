# GitHub Codespaces

GitHub Codespaces est une plateforme web qui nous permet de fournir un environnement préconfiguré pour la formation, soutenu par des machines virtuelles dans le cloud.
La plateforme est exploitée par Github (qui appartient à Microsoft), et est accessible gratuitement (avec des quotas d'utilisation) à toute personne disposant d'un compte Github.

!!! warning "Avertissement"

    Les comptes rattachés à des organisations peuvent être soumis à certaines restrictions supplémentaires.
    Si c'est votre cas, vous devrez peut-être utiliser un compte personnel indépendant, ou utiliser une installation locale à la place.

## Créer un compte GitHub

Vous pouvez créer un compte GitHub gratuit depuis la [page d'accueil GitHub](https://github.com/).

## Lancer votre GitHub Codespace

Une fois connecté·e à GitHub, ouvrez ce lien dans votre navigateur pour ouvrir l'environnement de formation Nextflow : <https://codespaces.new/nextflow-io/training?quickstart=1&ref=master>

Alternativement, vous pouvez cliquer sur le bouton affiché ci-dessous, qui est répété dans chaque cours de formation (généralement sur la page d'orientation).

[![Open in GitHub Codespaces](https://github.com/codespaces/badge.svg)](https://codespaces.new/nextflow-io/training?quickstart=1&ref=master)

Une page devrait s'afficher où vous pouvez créer un nouveau GitHub Codespace :

![Create a GitHub Codespace](img/codespaces_create.png)

### Configuration

Pour une utilisation générale, vous ne devriez pas avoir besoin de configurer quoi que ce soit.
Sauf indication contraire dans le cours que vous commencez, vous pouvez simplement cliquer sur le bouton principal pour continuer.

Cependant, il est possible de personnaliser l'environnement en cliquant sur le bouton « Change options ».

??? info "Options de configuration"

    Si vous cliquez sur le bouton « Change options », vous aurez la possibilité de personnaliser les éléments suivants :

    #### Branch

    Cela vous permet de sélectionner une version différente des supports de formation.
    La branche `master` contient généralement des corrections de bugs et des supports qui ont été récemment développés et approuvés mais qui n'ont pas encore été publiés sur le site web.
    Les autres branches contiennent des travaux en cours qui peuvent ne pas être entièrement fonctionnels.

    #### Machine type

    Cela vous permet de personnaliser la machine virtuelle que vous utiliserez pour suivre la formation.

    L'utilisation d'une machine avec plus de cœurs vous permet de mieux profiter de la capacité de Nextflow à paralléliser l'exécution du workflow.
    Cependant, cela consommera votre quota gratuit plus rapidement, nous ne recommandons donc pas de modifier ce paramètre sauf si cela est conseillé dans les instructions du cours que vous prévoyez de suivre.

    Voir « Quotas GitHub Codespaces » ci-dessous pour plus de détails sur les quotas.

### Temps de démarrage

L'ouverture d'un nouvel environnement GitHub Codespaces pour la première fois peut prendre plusieurs minutes, car le système doit configurer votre machine virtuelle, donc ne vous inquiétez pas s'il y a un temps d'attente.
Cependant, cela ne devrait pas prendre plus de cinq minutes.

## Naviguer dans l'interface de formation

Une fois votre GitHub Codespaces chargé, vous devriez voir quelque chose de similaire à ce qui suit (qui peut s'ouvrir en mode clair selon les préférences de votre compte) :

![GitHub Codespaces welcome](img/codespaces_welcome.png)

Il s'agit de l'interface de l'IDE VSCode, une application de développement de code populaire que nous recommandons d'utiliser pour le développement Nextflow.

- **L'éditeur principal** est l'endroit où le code Nextflow et les autres fichiers texte s'ouvriront. C'est ici que vous éditerez le code. Lorsque vous ouvrez le codespace, cela vous montrera un aperçu du fichier `README.md`.
- **Le terminal** sous l'éditeur principal vous permet d'exécuter des commandes. C'est ici que vous exécuterez toutes les lignes de commande données dans les instructions du cours.
- **La barre latérale** vous permet de personnaliser votre environnement et d'effectuer des tâches de base (copier, coller, ouvrir des fichiers, rechercher, git, etc.). Par défaut, elle est ouverte sur l'explorateur de fichiers, qui vous permet de parcourir le contenu du dépôt. Cliquer sur un fichier dans l'explorateur l'ouvrira dans la fenêtre de l'éditeur principal.

Vous pouvez ajuster les proportions relatives des panneaux de fenêtre comme vous le souhaitez.

<!-- TODO (future) Link to development best practices side quest? -->

## Autres remarques sur l'utilisation de GitHub Codespaces

### Reprendre une session

Une fois que vous avez créé un environnement, vous pouvez facilement le reprendre ou le redémarrer et continuer là où vous vous étiez arrêté·e.
Votre environnement expirera après 30 minutes d'inactivité et sauvegardera vos modifications pendant jusqu'à 2 semaines.

Vous pouvez rouvrir un environnement depuis <https://github.com/codespaces/>.
Les environnements précédents seront listés.
Cliquez sur une session pour la reprendre.

![List GitHub Codespace sessions](img/codespaces_list.png)

Si vous avez sauvegardé l'URL de votre environnement GitHub Codespaces précédent, vous pouvez simplement l'ouvrir dans votre navigateur.
Alternativement, cliquez sur le même bouton que vous avez utilisé pour le créer en premier lieu :

[![Open in GitHub Codespaces](https://github.com/codespaces/badge.svg)](https://codespaces.new/nextflow-io/training?quickstart=1&ref=master)

Vous devriez voir la session précédente, l'option par défaut est de la reprendre :

![Resume a GitHub Codespace](img/codespaces_resume.png)

### Sauvegarder des fichiers sur votre machine locale

Pour sauvegarder n'importe quel fichier depuis le panneau de l'explorateur, faites un clic droit sur le fichier et sélectionnez `Download`.

### Gérer les quotas GitHub Codespaces

GitHub Codespaces vous donne jusqu'à 15 Go-mois de stockage par mois, et 120 heures-cœur par mois.
Cela équivaut à environ 60 heures d'exécution de l'environnement par défaut utilisant l'espace de travail standard (2 cœurs, 8 Go de RAM et 32 Go de stockage).

Vous pouvez les créer avec plus de ressources (voir l'explication ci-dessus), mais cela consommera votre utilisation gratuite plus rapidement et vous aurez moins d'heures d'accès à cet espace.
Par exemple, si vous sélectionnez une machine à 4 cœurs au lieu des 2 cœurs par défaut, votre quota s'épuisera deux fois plus vite.

En option, vous pouvez acheter l'accès à plus de ressources.

Pour plus d'informations, consultez la documentation GitHub :
[About billing for GitHub Codespaces](https://docs.github.com/en/billing/managing-billing-for-your-products/managing-billing-for-github-codespaces/about-billing-for-github-codespaces)
