# Partie 3 : Profilage et optimisation des ressources

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Traduction assistée par IA - [en savoir plus et suggérer des améliorations](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

CECI EST UN ESPACE RÉSERVÉ

!!!note "Note"

    Ce module de formation est en cours de refonte.

---

TODO

### 1.1. Exécuter le workflow pour générer un rapport d'utilisation des ressources

Pour que Nextflow génère le rapport automatiquement, ajoutez simplement `-with-report <nom_fichier>.html` à votre ligne de commande.

```bash
nextflow run main.nf -profile my_laptop -with-report report-config-1.html
```

Le rapport est un fichier html, que vous pouvez télécharger et ouvrir dans votre navigateur. Vous pouvez également faire un clic droit dessus dans l'explorateur de fichiers à gauche et cliquer sur `Show preview` afin de le visualiser dans VS Code.

Prenez quelques minutes pour parcourir le rapport et voir si vous pouvez identifier des opportunités d'ajustement des ressources.
Assurez-vous de cliquer sur les onglets qui montrent les résultats d'utilisation en pourcentage de ce qui a été alloué.
Il existe une [documentation](https://www.nextflow.io/docs/latest/reports.html) décrivant toutes les fonctionnalités disponibles.

<!-- TODO: insert images -->

Une observation est que le processus `GATK_JOINTGENOTYPING` semble être très gourmand en CPU, ce qui est logique puisqu'il effectue de nombreux calculs complexes.
Nous pourrions donc essayer d'augmenter cette allocation et voir si cela réduit le temps d'exécution.

Cependant, nous semblons avoir surestimé les allocations de mémoire ; tous les processus n'utilisent qu'une fraction de ce que nous leur donnons.
Nous devrions réduire cela et économiser des ressources.

### 1.2. Ajuster les allocations de ressources pour un processus spécifique

Nous pouvons spécifier les allocations de ressources pour un processus donné en utilisant le sélecteur de processus `withName`.
La syntaxe ressemble à ceci lorsqu'elle est seule dans un bloc process :

```groovy title="Syntaxe"
process {
    withName: 'GATK_JOINTGENOTYPING' {
        cpus = 4
    }
}
```

Ajoutons cela au bloc process existant dans le fichier `nextflow.config`.

```groovy title="nextflow.config" linenums="11"
process {
    // valeurs par défaut pour tous les processus
    cpus = 2
    memory = 2.GB
    // allocations pour un processus spécifique
    withName: 'GATK_JOINTGENOTYPING' {
        cpus = 4
    }
}
```

Avec cette spécification, les paramètres par défaut s'appliqueront à tous les processus **sauf** au processus `GATK_JOINTGENOTYPING`, qui est un cas particulier recevant beaucoup plus de CPU.
Espérons que cela aura un effet.

### 1.3. Exécuter à nouveau avec la configuration modifiée

Exécutons à nouveau le workflow avec la configuration modifiée et avec l'indicateur de rapport activé, mais notez que nous donnons un nom différent au rapport afin de pouvoir les différencier.

```bash
nextflow run main.nf -profile my_laptop -with-report report-config-2.html
```

Une fois de plus, vous ne remarquerez probablement pas de différence substantielle dans le temps d'exécution, car il s'agit d'une charge de travail si petite et les outils passent plus de temps dans des tâches auxiliaires que dans l'exécution du « vrai » travail.

Cependant, le deuxième rapport montre que notre utilisation des ressources est maintenant plus équilibrée.

<!-- **TODO: screenshots?** -->

Comme vous pouvez le voir, cette approche est utile lorsque vos processus ont des besoins en ressources différents. Elle vous permet d'ajuster précisément les allocations de ressources que vous configurez pour chaque processus en fonction de données réelles, et non de suppositions.

!!!note "Note"

    Ceci n'est qu'un petit aperçu de ce que vous pouvez faire pour optimiser votre utilisation des ressources.
    Nextflow lui-même dispose d'une [logique de nouvelle tentative dynamique](https://www.nextflow.io/docs/latest/process.html#dynamic-task-resources) vraiment intéressante intégrée pour réessayer les tâches qui échouent en raison de limitations de ressources.
    De plus, Seqera Platform offre des outils pilotés par l'IA pour optimiser automatiquement vos allocations de ressources également.

    Nous aborderons ces deux approches dans une prochaine partie de ce cours de formation.

Cela étant dit, il peut y avoir certaines contraintes sur ce que vous pouvez (ou devez) allouer selon l'exécuteur de calcul et l'infrastructure de calcul que vous utilisez. Par exemple, votre cluster peut vous obliger à rester dans certaines limites qui ne s'appliquent pas lorsque vous exécutez ailleurs.
