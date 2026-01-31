# Partie 4 : Configuration

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Traduction assistée par IA - [en savoir plus et suggérer des améliorations](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

Dans les Parties 1-3, nous avons appris à exécuter Nextflow, à lancer un pipeline nf-core et à gérer les entrées avec des fichiers de paramètres et des samplesheets.
Nous allons maintenant explorer comment configurer les pipelines pour différents environnements informatiques en utilisant des **fichiers de configuration** et des **profils**.

## Objectifs d'apprentissage

À la fin de cette partie, vous serez capable de :

- Comprendre comment Nextflow résout la configuration à partir de plusieurs sources
- Utiliser les profils intégrés nf-core pour les conteneurs et les tests
- Créer des profils personnalisés pour différents environnements informatiques
- Personnaliser les demandes de ressources en utilisant les étiquettes de processus
- Gérer les limites de ressources dans des environnements contraints
- Inspecter la configuration résolue avec `nextflow config`

---

## 1. Comprendre la configuration Nextflow

### 1.1. Qu'est-ce qu'un fichier de configuration ?

Nextflow utilise des fichiers de configuration pour séparer la **logique du workflow** (quoi faire) des **paramètres d'exécution** (comment et où le faire).

Les fichiers de configuration contrôlent :

- Les moteurs de conteneurs (Docker, Singularity, Conda)
- Les ressources de calcul (CPUs, mémoire, temps)
- Les plateformes d'exécution (local, HPC, cloud)
- Les paramètres du pipeline

### 1.2. Priorité de configuration

Nextflow charge la configuration à partir de plusieurs sources, les sources ultérieures remplaçant les précédentes :

1. **Configuration du pipeline** : `nextflow.config` dans le dépôt du pipeline
2. **Configuration du répertoire** : `nextflow.config` dans votre répertoire de travail actuel
3. **Configuration utilisateur** : `~/.nextflow/config`
4. **Ligne de commande** : Paramètres et options passés directement

Cette approche en couches vous permet de conserver les valeurs par défaut dans le pipeline, de les remplacer par des paramètres spécifiques à l'utilisateur et d'effectuer des ajustements rapides en ligne de commande.

### 1.3. Notre configuration actuelle

Examinons la configuration que nous avons utilisée :

```groovy title="nextflow.config"
docker.enabled = true
process {
    resourceLimits = [
        cpus: 2,
        memory: '7.GB',
    ]
}

```

Commentons ou modifions la ligne `docker.enabled = true` de la Partie 2, et voyons comment nous pouvons obtenir le même résultat en utilisant un profil dans molkart à la place.

---

## 2. Utilisation des profils

### 2.1. Qu'est-ce que les profils ?

Les profils sont des ensembles nommés de configuration qui peuvent être activés avec le drapeau `-profile` via la commande `nextflow run`.
Ils facilitent le passage d'un scénario de calcul à l'autre sans modifier les fichiers de configuration.

Tous les pipelines nf-core sont fournis avec plusieurs profils par défaut que nous pouvons utiliser.

### 2.2. Inspection des profils intégrés

Inspectons-les dans le fichier `molkart/nextflow.config` associé au code du pipeline :

```bash
code molkart/nextflow.config
```

Recherchez le bloc `profiles` :

```groovy title="molkart/nextflow.config (extrait)"
profiles {
    docker {
        docker.enabled          = true
        singularity.enabled     = false
        conda.enabled           = false
    }
    singularity {
        singularity.enabled     = true
        singularity.autoMounts  = true
        docker.enabled          = false
    }
    conda {
        conda.enabled           = true
        docker.enabled          = false
        conda.channels          = ['conda-forge', 'bioconda']
    }
}
```

Profils de conteneurs courants :

- `docker` : Utilise des conteneurs Docker (le plus courant pour le développement local)
- `singularity` : Utilise Singularity/Apptainer (courant sur HPC)
- `conda` : Utilise des environnements Conda
- `apptainer` : Utilise des conteneurs Apptainer

### 2.3. Ré-exécution avec des profils au lieu de nextflow.config

Maintenant que nous avons désactivé la configuration docker dans notre fichier `nextflow.config` local et que nous comprenons les profils, ré-exécutons le pipeline en utilisant le drapeau `-profile`.

Précédemment dans la Partie 3, nous avons créé un fichier `params.yaml` avec nos paramètres personnalisés.
Nous pouvons maintenant le combiner avec le profil Docker intégré :

```bash
nextflow run ./molkart \
  -profile docker \
  -params-file params.yaml \
  -resume
```

Décomposons ce que fait chaque drapeau :

- `-profile docker` : Active le profil Docker du fichier `nextflow.config` de molkart, qui définit `docker.enabled = true`
- `-params-file params.yaml` : Charge tous les paramètres du pipeline à partir de notre fichier YAML
- `-resume` : Réutilise les résultats en cache des exécutions précédentes

Comme nous utilisons `-resume`, Nextflow vérifiera si quelque chose a changé depuis la dernière exécution.
Si les paramètres, les entrées et le code sont identiques, toutes les tâches seront récupérées du cache et le pipeline se terminera presque instantanément.

```console title="Sortie (extrait)"
executor >  local (12)
...
[1a/2b3c4d] NFCORE_MOLKART:MOLKART:MINDAGAP_MINDAGAP (mem_only)   [100%] 2 of 2, cached: 2 ✔
[5e/6f7g8h] NFCORE_MOLKART:MOLKART:CLAHE (mem_only)               [100%] 2 of 2, cached: 2 ✔
...
-[nf-core/molkart] Pipeline completed successfully-
```

Remarquez que tous les processus affichent `cached: 2` ou `cached: 1` - rien n'a été ré-exécuté !

### 2.4. Profils de test

Les profils de test fournissent des moyens rapides de spécifier des paramètres d'entrée et des fichiers de données par défaut pour vous permettre de vérifier que le pipeline fonctionne.
Les pipelines nf-core incluront toujours au moins deux profils de test :

- `test` : Petit jeu de données avec des paramètres rapides pour des tests rapides
- `test_full` : Test plus complet avec des données plus volumineuses

Examinons de plus près le profil `test` dans molkart qui est inclus à l'aide de la directive `includeConfig` :

```groovy title="molkart/nextflow.config (extrait)"
profiles {
  ...
    test      { includeConfig 'conf/test.config'      }
}
```

Cela signifie que chaque fois que nous exécutons le pipeline avec `-profile test`, Nextflow chargera la configuration depuis `conf/test.config`.

```groovy title="molkart/conf/test.config (extrait)"
params {
    config_profile_name        = 'Test profile'
    config_profile_description = 'Minimal test dataset to check pipeline function'

    input = 'https://raw.githubusercontent.com/nf-core/test-datasets/molkart/test_data/samplesheets/samplesheet_membrane.csv'
    mindagap_tilesize = 90
    mindagap_boxsize = 7
    mindagap_loopnum = 100
    clahe_pyramid_tile = 368
    segmentation_method = "mesmer,cellpose,stardist"
}

process {
    resourceLimits = [
        cpus: 4,
        memory: '15.GB',
        time: '1.h'
    ]
}
```

Remarquez que ce profil contient les mêmes paramètres que ceux que nous avons utilisés dans notre fichier `params.yaml` précédemment.

Vous pouvez activer plusieurs profils en les séparant par des virgules.
Utilisons cela pour tester notre pipeline sans avoir besoin de notre fichier de paramètres :

```bash
nextflow run ./molkart -profile docker,test --outdir results -resume
```

Cela combine :

- `docker` : Active les conteneurs Docker
- `test` : Utilise le jeu de données et les paramètres de test

Les profils sont appliqués de gauche à droite, donc les profils ultérieurs remplacent les précédents s'ils définissent les mêmes valeurs.

### À retenir

Les pipelines nf-core sont fournis avec des profils intégrés pour les conteneurs, les tests et les environnements spéciaux.
Vous pouvez combiner plusieurs profils pour construire la configuration dont vous avez besoin.

### Et ensuite ?

Apprenez à créer vos propres profils personnalisés pour différents environnements informatiques.

---

## 3. Création de profils personnalisés

### 3.1. Créer des profils pour basculer entre le développement local et l'exécution sur HPC

Créons des profils personnalisés pour deux scénarios :

1. Développement local avec Docker
2. HPC universitaire avec planificateur Slurm et Singularity

Ajoutez ce qui suit à votre `nextflow.config` :

```groovy title="nextflow.config"
profiles {
    local_dev {
        docker.enabled = true
        process.executor = 'local'
    }

    hpc_cluster {
        singularity.enabled = true
        process.executor = 'slurm'
        process.queue = 'standard_queue'
        singularity.cacheDir = '/shared/containers'
    }
}
```

Vous pouvez maintenant basculer facilement entre les environnements :

```bash
# Pour le développement local
nextflow run ./molkart -profile local_dev --input data/samplesheet.csv --outdir results

# Pour HPC (lorsque disponible)
nextflow run ./molkart -profile hpc_cluster --input data/samplesheet.csv --outdir results
```

!!! Note

    Nous ne pouvons pas tester le profil HPC dans cet environnement de formation car nous n'avons pas accès à un planificateur Slurm.
    Mais cela montre comment vous le configureriez pour une utilisation réelle.

### 3.2. Utiliser `nextflow config` pour inspecter la configuration

La commande `nextflow config` affiche la configuration entièrement résolue sans exécuter le pipeline.

Afficher la configuration par défaut :

```bash
nextflow config ./molkart
```

Afficher la configuration avec un profil spécifique :

```bash
nextflow config -profile local_dev ./molkart
```

Ceci est extrêmement utile pour :

- Déboguer les problèmes de configuration
- Comprendre quelles valeurs seront réellement utilisées
- Vérifier comment plusieurs profils interagissent

### À retenir

Les profils personnalisés vous permettent de basculer entre différents environnements informatiques avec un seul drapeau en ligne de commande.
Utilisez `nextflow config` pour inspecter la configuration résolue avant l'exécution.

### Et ensuite ?

Apprenez à personnaliser les demandes de ressources pour des processus individuels en utilisant le système d'étiquettes de processus de nf-core.

---

## 4. Personnalisation des demandes de ressources

### 4.1. Comprendre les étiquettes de processus dans les pipelines nf-core

Pour plus de simplicité, les pipelines nf-core utilisent des [**étiquettes de processus**](https://www.nextflow.io/docs/latest/reference/process.html#process-label) pour standardiser l'allocation de ressources dans tous les pipelines.
Chaque processus est étiqueté avec une étiquette comme `process_low`, `process_medium` ou `process_high` pour décrire respectivement des besoins en ressources de calcul faibles, moyens ou élevés.
Ces étiquettes sont converties en demandes de ressources spécifiques dans l'un des fichiers de configuration situés dans le répertoire `conf/` du pipeline.

```groovy title="molkart/conf/base.config (extrait)"
process {
    cpus   = { 1      * task.attempt }
    memory = { 6.GB   * task.attempt }
    time   = { 4.h    * task.attempt }

    errorStrategy = { task.exitStatus in ((130..145) + 104) ? 'retry' : 'finish' }
    maxRetries    = 1
    maxErrors     = '-1'

    withLabel:process_single {
        cpus   = { 1                   }
        memory = { 6.GB * task.attempt }
        time   = { 4.h  * task.attempt }
    }
    withLabel:process_low {
        cpus   = { 2     * task.attempt }
        memory = { 12.GB * task.attempt }
        time   = { 4.h   * task.attempt }
    }
    withLabel:process_medium {
        cpus   = { 6     * task.attempt }
        memory = { 36.GB * task.attempt }
        time   = { 8.h   * task.attempt }
    }
    withLabel:process_high {
        cpus   = { 12    * task.attempt }
        memory = { 72.GB * task.attempt }
        time   = { 16.h  * task.attempt }
    }
}
```

Remarquez le multiplicateur `task.attempt` - cela permet aux nouvelles tentatives de tâches ultérieures de demander plus de ressources, si le pipeline est configuré avec `process.maxRetries > 1`.

### 4.2. Remplacement des ressources pour des processus spécifiques

Pour un contrôle précis, ciblez des processus individuels par nom :

```groovy title="nextflow.config"
process {
    withName: 'NFCORE_MOLKART:MOLKART:CELLPOSE' {
        cpus   = 16
        memory = 32.GB
    }
}
```

Si nous essayons d'exécuter ce pipeline avec le remplacement ci-dessus, le processus `CELLPOSE` demandera 16 CPUs et 32 GB de mémoire au lieu des valeurs par défaut définies par son étiquette.
Cela provoquera l'échec du pipeline dans notre environnement actuel car nous n'avons pas autant de RAM disponible.
Nous apprendrons comment prévenir ces types d'échecs dans la section suivante.

!!! Tip "Astuce"

    Pour trouver les noms de processus, consultez la sortie d'exécution du pipeline ou vérifiez `.nextflow.log`.
    Les noms de processus suivent le modèle `WORKFLOW:SUBWORKFLOW:PROCESS`.

### À retenir

Les pipelines nf-core utilisent des étiquettes de processus pour standardiser l'allocation de ressources.
Vous pouvez remplacer les ressources par étiquette (affecte plusieurs processus) ou par nom (affecte un processus spécifique).

### Et ensuite ?

Apprenez à gérer les limites de ressources dans des environnements contraints comme GitHub Codespaces.

---

## 5. Gestion des ressources dans des environnements contraints

### 5.1. Le problème des limites de ressources

Si nous essayions d'exécuter molkart avec un processus demandant 16 CPUs et 32 GB de mémoire (comme montré dans la section 4.2), cela échouerait dans notre environnement actuel car nous n'avons pas autant de ressources disponibles.
Dans un environnement de cluster avec des nœuds plus grands, de telles demandes seraient soumises au planificateur.

Dans des environnements contraints comme GitHub Codespaces, sans limites, Nextflow refuserait d'exécuter des processus qui dépassent les ressources disponibles.

### 5.2. Définition des limites de ressources

La directive `resourceLimits` plafonne les demandes de ressources aux valeurs spécifiées :

```groovy title="nextflow.config"
process {
    resourceLimits = [ cpus: 2, memory: 7.GB ]
}
```

Cela indique à Nextflow : « Si un processus demande plus de 2 CPUs ou 7 GB de mémoire, limitez-le plutôt à ces limites. »

### 5.3. Ajout de limites de ressources aux profils personnalisés

Mettez à jour vos profils personnalisés pour inclure des limites appropriées :

```groovy title="nextflow.config"
profiles {
    local_dev {
        docker.enabled = true
        process.executor = 'local'
        process.resourceLimits = [
            cpus: 2,
            memory: 7.GB
        ]
    }

    hpc_cluster {
        singularity.enabled = true
        process.executor = 'slurm'
        process.queue = 'batch'
        process.resourceLimits = [
            cpus: 32,
            memory: 128.GB,
            time: 24.h
        ]
    }
}
```

!!! Warning "Avertissement"

    Définir des limites de ressources trop basses peut entraîner l'échec ou le ralentissement des processus.
    Le pipeline peut avoir besoin d'utiliser des algorithmes moins gourmands en mémoire ou de traiter les données par plus petits morceaux.

### À retenir

Utilisez `resourceLimits` pour exécuter des pipelines dans des environnements à ressources limitées en plafonnant les demandes de ressources des processus.
Différents profils peuvent avoir des limites différentes adaptées à leur environnement.

### Et ensuite ?

Vous avez terminé la formation de base Nextflow pour la bio-imagerie !

---

## Conclusion

Vous comprenez maintenant comment configurer les pipelines Nextflow pour différents environnements informatiques.

Compétences clés que vous avez acquises :

- **Priorité de configuration** : Comment Nextflow résout les paramètres à partir de plusieurs sources
- **Profils nf-core** : Utilisation des profils intégrés pour les conteneurs, les tests et les utilitaires
- **Profils personnalisés** : Création de vos propres profils pour différents environnements
- **Étiquettes de processus** : Compréhension et remplacement des demandes de ressources par étiquette
- **Limites de ressources** : Gestion d'environnements contraints avec `resourceLimits`
- **Inspection de configuration** : Utilisation de `nextflow config` pour déboguer et vérifier les paramètres

Ces compétences de configuration sont transférables à tout pipeline Nextflow et vous aideront à exécuter des workflows efficacement sur des machines locales, des clusters HPC et des plateformes cloud.

### Et ensuite ?

Félicitations pour avoir terminé le cours Nextflow pour la bio-imagerie !

Prochaines étapes :

- Remplissez le questionnaire du cours pour fournir des commentaires
- Consultez [Hello Nextflow](../hello_nextflow/index.md) pour en savoir plus sur le développement de workflows
- Explorez [Hello nf-core](../hello_nf-core/index.md) pour approfondir les outils nf-core
- Parcourez d'autres cours dans les [collections de formation](../training_collections/index.md)
