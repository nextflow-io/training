# Partie 2 : Exécuter nf-core/molkart

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Traduction assistée par IA - [en savoir plus et suggérer des améliorations](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

Dans la Partie 1, nous avons exécuté un workflow simple Hello World pour comprendre les bases de l'exécution de Nextflow.
Nous allons maintenant exécuter un véritable pipeline de bio-imagerie : **nf-core/molkart**.

Ce pipeline traite les données de transcriptomique spatiale Molecular Cartography de Resolve Bioscience.
Cependant, les modèles Nextflow que vous apprendrez ici s'appliquent à n'importe quel pipeline nf-core ou workflow de production.

## 1. Comprendre les pipelines nf-core

Avant d'exécuter le pipeline, comprenons ce qu'est nf-core et pourquoi c'est important pour l'exécution de workflows.

### 1.1. Qu'est-ce que nf-core ?

[nf-core](https://nf-co.re/) est une collection pilotée par la communauté de pipelines Nextflow de haute qualité.
Tous les pipelines nf-core suivent la même structure et les mêmes conventions, ce qui signifie qu'une fois que vous avez appris à en exécuter un, vous pouvez tous les exécuter.

Caractéristiques clés des pipelines nf-core :

- **Structure standardisée** : Tous les pipelines ont des noms de paramètres et des modèles d'utilisation cohérents
- **Données de test intégrées** : Chaque pipeline inclut des profils de test pour une validation rapide
- **Documentation complète** : Instructions d'utilisation détaillées et descriptions des paramètres
- **Contrôle qualité** : Rapports de QC automatisés utilisant MultiQC
- **Support des conteneurs** : Conteneurs pré-construits pour la reproductibilité

!!! tip "Vous voulez en savoir plus sur nf-core ?"

    Pour une introduction approfondie au développement de pipelines nf-core, consultez le cours de formation [Hello nf-core](../../hello_nf-core/index.md).
    Il couvre comment créer et personnaliser des pipelines nf-core à partir de zéro.

### 1.2. Le pipeline molkart

![Pipeline nf-core/molkart](img/molkart.png)

Le pipeline [nf-core/molkart](https://nf-co.re/molkart) traite les données d'imagerie de transcriptomique spatiale à travers plusieurs étapes :

1. **Prétraitement d'image** : Remplissage de motif de grille et amélioration optionnelle du contraste
2. **Segmentation cellulaire** : Plusieurs options d'algorithmes (Cellpose, Mesmer, ilastik, Stardist)
3. **Attribution de spots** : Attribuer les spots de transcrit aux cellules segmentées
4. **Contrôle qualité** : Générer des rapports de QC complets

Les sorties clés sont :

- Tables de comptage cellule par transcrit
- Masques de segmentation
- Rapport de contrôle qualité MultiQC

---

## 2. Exécuter molkart avec des données de test

Avant de commencer, clonons le dépôt molkart localement afin de pouvoir inspecter son code :

```bash
cd /workspaces/training/nf4-science/imaging
git clone --branch 1.2.0 --depth 1 https://github.com/nf-core/molkart
```

Cela crée un répertoire `molkart/` contenant le code source complet du pipeline.

!!! note "Pourquoi clonons-nous localement ?"

    Typiquement, vous exécuteriez les pipelines nf-core directement depuis GitHub en utilisant `nextflow run nf-core/molkart -r 1.2.0`.
    Nextflow télécharge automatiquement la version du pipeline demandée pour vous dans `$HOME/.nextflow/assets/nf-core/molkart` et l'exécute à partir de là.
    Cependant, pour cette formation, nous clonons le pipeline dans un répertoire local différent afin de pouvoir inspecter plus facilement le code.

### 2.1. Comprendre les exigences en matière de conteneurs

Avant d'exécuter le pipeline complet, apprenons pourquoi les conteneurs sont essentiels pour les pipelines nf-core.

Essayons d'exécuter le pipeline en utilisant l'ensemble de données de test et les paramètres de la configuration de test molkart :

```bash
nextflow run ./molkart \
  --input 'data/samplesheet.csv' \
  --mindagap_tilesize 90 \
  --mindagap_boxsize 7 \
  --mindagap_loopnum 100 \
  --clahe_pyramid_tile 368 \
  --segmentation_method "mesmer,cellpose,stardist" \
  --outdir results
```

Décomposons ces paramètres :

- `--input` : Chemin vers la feuille d'échantillons contenant les métadonnées des échantillons
- `--mindagap_tilesize`, `--mindagap_boxsize`, `--mindagap_loopnum` : Paramètres pour le remplissage de motif de grille
- `--clahe_pyramid_tile` : Taille du noyau pour l'amélioration du contraste
- `--segmentation_method` : Quel(s) algorithme(s) utiliser pour la segmentation cellulaire
- `--outdir` : Où enregistrer les résultats

!!! Warning "Cette commande échouera - c'est intentionnel !"

    Nous exécutons délibérément cela sans conteneurs pour démontrer pourquoi ils sont nécessaires.

Après quelques instants, vous verrez une erreur comme celle-ci :

??? failure "Sortie de la commande"

    ```console
    ERROR ~ Error executing process > 'NFCORE_MOLKART:MOLKART:MINDAGAP_DUPLICATEFINDER (mem_only)'

    Caused by:
      Process `NFCORE_MOLKART:MOLKART:MINDAGAP_DUPLICATEFINDER (mem_only)` terminated with an error exit status (127)

    Command executed:

      duplicate_finder.py \
          spots.txt \
          90

    Command exit status:
      127

    Command error:
      .command.sh: line 3: duplicate_finder.py: command not found
    ```

**Que se passe-t-il ici ?**

L'erreur `command not found` (statut de sortie 127) signifie que Nextflow a essayé d'exécuter `duplicate_finder.py` mais n'a pas pu le trouver sur votre système.
C'est parce que :

1. Le pipeline attend que des logiciels bioinformatiques spécialisés soient installés
2. Ces outils (comme `duplicate_finder.py`, `apply_clahe.dask.py`, etc.) ne font pas partie des distributions Linux standard
3. Sans conteneurs, Nextflow essaie d'exécuter les commandes directement sur votre machine locale

**D'où ces outils sont-ils censés provenir ?**

Inspectons l'un des modules de processus pour voir comment il déclare ses exigences logicielles.

Ouvrez le module de prétraitement CLAHE :

```bash
code molkart/modules/local/clahe/main.nf
```

Regardez la ligne 5 - vous verrez :

```groovy
container 'ghcr.io/schapirolabor/molkart-local:v0.0.4'
```

Cette ligne indique à Nextflow : « Pour exécuter ce processus, utilisez l'image Docker `ghcr.io/schapirolabor/molkart-local:v0.0.4`, qui contient tous les logiciels requis. »

Chaque processus déclare quelle image de conteneur fournit ses outils requis.
Cependant, Nextflow n'utilise ces conteneurs que si vous lui dites de le faire !

**La solution : Activer Docker dans la configuration**

### 2.2. Configurer Docker et lancer le pipeline

Pour activer Docker, nous devons changer `docker.enabled` de `false` à `true` dans le fichier `nextflow.config`.

Ouvrez le fichier de configuration :

```bash
code nextflow.config
```

Changez `docker.enabled = false` en `docker.enabled = true` :

```groovy
docker.enabled = true
process {
    resourceLimits = [
        cpus: 2,
        memory: '7.GB',
    ]
}
```

Maintenant, exécutez à nouveau le pipeline avec la même commande :

```bash
nextflow run ./molkart \
  --input 'data/samplesheet.csv' \
  --mindagap_tilesize 90 \
  --mindagap_boxsize 7 \
  --mindagap_loopnum 100 \
  --clahe_pyramid_tile 368 \
  --segmentation_method "cellpose,mesmer,stardist" \
  --outdir results
```

Cette fois, Nextflow va :

1. Lire le paramètre `docker.enabled = true` de la configuration
2. Récupérer les images Docker requises (première fois seulement)
3. Exécuter chaque processus à l'intérieur de son conteneur spécifié
4. S'exécuter avec succès car tous les outils sont disponibles à l'intérieur des conteneurs

!!! Tip "Pourquoi les conteneurs sont importants"

    La plupart des pipelines nf-core **nécessitent** la conteneurisation (Docker, Singularity, Podman, etc.) car :

    - Ils utilisent des logiciels bioinformatiques spécialisés non disponibles dans les environnements standard
    - Les conteneurs garantissent la reproductibilité - exactement les mêmes versions de logiciels s'exécutent partout
    - Vous n'avez pas besoin d'installer manuellement des dizaines d'outils et leurs dépendances

    Pour plus de détails sur les conteneurs dans Nextflow, consultez [Hello Containers](../../hello_nextflow/05_hello_containers.md) de la formation Hello Nextflow.

### 2.3. Surveiller l'exécution

Pendant l'exécution du pipeline, vous verrez une sortie similaire à ceci :

??? success "Sortie de la commande"

    ```console
    Nextflow 25.04.8 is available - Please consider updating your version to it

    N E X T F L O W   ~  version 25.04.3

    Launching `https://github.com/nf-core/molkart` [soggy_kalam] DSL2 - revision: 5e54b29cb3 [dev]


    ------------------------------------------------------
                                            ,--./,-.
            ___     __   __   __   ___     /,-._.--~'
      |\ | |__  __ /  ` /  \ |__) |__         }  {
      | \| |       \__, \__/ |  \ |___     \`-._,-`-,
                                            `._,._,'
      nf-core/molkart 1.2.0dev
    ------------------------------------------------------
    Segmentation methods and options
      segmentation_method       : mesmer,cellpose,stardist

    Image preprocessing
      mindagap_boxsize          : 7
      mindagap_loopnum          : 100
      clahe_kernel              : 25
      mindagap_tilesize         : 90
      clahe_pyramid_tile        : 368

    Input/output options
      input                     : https://raw.githubusercontent.com/nf-core/test-datasets/molkart/test_data/samplesheets/samplesheet_membrane.csv
      outdir                    : results

    Institutional config options
      config_profile_name       : Test profile
      config_profile_description: Minimal test dataset to check pipeline function

    Generic options
      trace_report_suffix       : 2025-10-18_22-22-21

    Core Nextflow options
      revision                  : dev
      runName                   : soggy_kalam
      containerEngine           : docker
      launchDir                 : /workspaces/training/nf4-science/imaging
      workDir                   : /workspaces/training/nf4-science/imaging/work
      projectDir                : /workspaces/.nextflow/assets/nf-core/molkart
      userName                  : root
      profile                   : docker,test
      configFiles               :

    !! Only displaying parameters that differ from the pipeline defaults !!
    ------------------------------------------------------
    * The pipeline
        https://doi.org/10.5281/zenodo.10650748

    * The nf-core framework
        https://doi.org/10.1038/s41587-020-0439-x

    * Software dependencies
        https://github.com/nf-core/molkart/blob/master/CITATIONS.md

    executor >  local (22)
    [c1/da5009] NFCORE_MOLKART:MOLKART:MINDAGAP_MINDAGAP (mem_only)        [100%] 2 of 2 ✔
    [73/8f5e8a] NFCORE_MOLKART:MOLKART:CLAHE (mem_only)                    [100%] 2 of 2 ✔
    [ec/8f84d5] NFCORE_MOLKART:MOLKART:CREATE_STACK (mem_only)             [100%] 1 of 1 ✔
    [a2/99349b] NFCORE_MOLKART:MOLKART:MINDAGAP_DUPLICATEFINDER (mem_only) [100%] 1 of 1 ✔
    [95/c9b4b1] NFCORE_MOLKART:MOLKART:DEEPCELL_MESMER (mem_only)          [100%] 1 of 1 ✔
    [d4/1ebd1e] NFCORE_MOLKART:MOLKART:STARDIST (mem_only)                 [100%] 1 of 1 ✔
    [3e/3c0736] NFCORE_MOLKART:MOLKART:CELLPOSE (mem_only)                 [100%] 1 of 1 ✔
    [a0/415c6a] NFCORE_MOLKART:MOLKART:MASKFILTER (mem_only)               [100%] 3 of 3 ✔
    [14/a830c9] NFCORE_MOLKART:MOLKART:SPOT2CELL (mem_only)                [100%] 3 of 3 ✔
    [b5/391836] NFCORE_MOLKART:MOLKART:CREATE_ANNDATA (mem_only)           [100%] 3 of 3 ✔
    [77/aed558] NFCORE_MOLKART:MOLKART:MOLKARTQC (mem_only)                [100%] 3 of 3 ✔
    [e6/b81475] NFCORE_MOLKART:MOLKART:MULTIQC                             [100%] 1 of 1 ✔
    -[nf-core/molkart] Pipeline completed successfully-
    Completed at: 19-Oct-2025 22:23:01
    Duration    : 2m 52s
    CPU hours   : 0.1
    Succeeded   : 22
    ```

Remarquez comment cette sortie est plus détaillée que notre exemple Hello World en raison des conventions nf-core que le pipeline suit :

- Le pipeline affiche sa version et son logo
- Les paramètres de configuration sont affichés
- Plusieurs processus s'exécutent en parallèle (indiqué par plusieurs lignes de processus)
- Les noms de processus incluent le chemin complet du module (par ex., `NFCORE_MOLKART:MOLKART:MINDAGAP_MINDAGAP`)

### 2.4. Comprendre l'exécution des processus

La ligne executor `executor > local (22)` vous indique :

- **executor** : Quel environnement de calcul est utilisé (`local` = votre machine)
- **(22)** : Nombre total de tâches lancées

Chaque ligne de processus montre :

- **Hash** (`[1a/2b3c4d]`) : Identifiant du répertoire de travail (comme avant)
- **Nom du processus** : Chemin complet du module et nom du processus
- **Identifiant d'entrée** : Nom de l'échantillon entre parenthèses
- **Progression** : Pourcentage terminé et comptage (par ex., `1 of 1 ✔`)

### À retenir

Vous savez comment lancer un pipeline nf-core avec des données de test et interpréter sa sortie d'exécution.

### Et maintenant ?

Apprenez où trouver les résultats et comment les interpréter.

---

## 3. Trouver et examiner les sorties

Lorsque le pipeline se termine avec succès, vous verrez un message de complétion et un résumé d'exécution.

### 3.1. Localiser le répertoire des résultats

Par défaut, les pipelines nf-core écrivent les sorties dans un répertoire spécifié par le paramètre `outdir`, que nous avons défini sur `results/`.

Listez le contenu :

```bash
tree results/
```

Vous devriez voir plusieurs sous-répertoires :

```console title="results/"
results/
├── anndata/
├── clahe/
├── mindagap/
├── molkartqc/
├── multiqc/
├── pipeline_info/
├── segmentation/
├── spot2cell/
└── stack/
```

Chaque sous-répertoire contient des sorties d'une étape spécifique du pipeline :

- **mindagap/** : Images remplies par grille de l'étape de prétraitement MindaGap
- **clahe/** : Images avec contraste amélioré du prétraitement CLAHE
- **stack/** : Piles d'images multi-canaux créées pour la segmentation
- **segmentation/** : Résultats de segmentation de différents algorithmes (cellpose/, mesmer/, stardist/, filtered_masks/)
- **spot2cell/** : Tables de comptage cellule par transcrit
- **anndata/** : Objets AnnData contenant des matrices cellule par transcrit et des coordonnées spatiales
- **molkartqc/** : Métriques de contrôle qualité pour l'attribution de spots
- **multiqc/** : Rapport complet de contrôle qualité
- **pipeline_info/** : Rapports d'exécution et journaux

### 3.2. Examiner le rapport MultiQC

Le rapport MultiQC est un fichier HTML complet qui agrège les métriques de qualité de toutes les étapes du pipeline.

Ouvrez le rapport dans l'explorateur de fichiers puis cliquez sur le bouton « Show Preview » pour le voir rendu directement dans VS Code.

Le rapport inclut :

- Statistiques générales pour tous les échantillons
- Métriques de prétraitement
- Métriques de qualité de segmentation
- Nombre de cellules et de spots détectés

!!! Tip

    Les rapports MultiQC sont généralement inclus dans tous les pipelines nf-core.
    Ils fournissent toujours un aperçu de haut niveau de l'exécution du pipeline et de la qualité des données.

### 3.3. Examiner les tables cellule par transcrit

La sortie scientifique la plus importante est la table de comptage cellule par transcrit.
Cela vous indique combien de chaque transcrit a été détecté dans chaque cellule.

Naviguez vers le répertoire spot2cell :

```bash
ls results/spot2cell/
```

Vous trouverez des fichiers comme :

- `cellxgene_mem_only_cellpose.csv` : Table cellule par transcrit utilisant la segmentation Cellpose
- `cellxgene_mem_only_mesmer.csv` : Table cellule par transcrit utilisant la segmentation Mesmer
- `cellxgene_mem_only_stardist.csv` : Table cellule par transcrit utilisant la segmentation Stardist

Nous n'avons exécuté qu'un seul échantillon dans cet ensemble de données de test, mais dans une expérience réelle, nous aurions ces tables pour chaque échantillon.
Remarquez comment Nextflow est capable de traiter plusieurs méthodes de segmentation en parallèle, facilitant la comparaison des résultats.

### 3.4. Voir les rapports d'exécution

Nextflow génère automatiquement plusieurs rapports d'exécution.

Vérifiez le répertoire pipeline_info :

```bash
ls results/pipeline_info/
```

Fichiers clés :

- **execution_report.html** : Chronologie et visualisation de l'utilisation des ressources
- **execution_timeline.html** : Diagramme de Gantt de l'exécution des processus
- **execution_trace.txt** : Métriques détaillées d'exécution des tâches
- **pipeline_dag.html** : Graphe acyclique dirigé montrant la structure du workflow

Ouvrez le rapport d'exécution pour voir l'utilisation des ressources :

```bash
code results/pipeline_info/execution_report.html
```

Cela montre :

- Combien de temps chaque processus a pris
- Utilisation du CPU et de la mémoire
- Quelles tâches ont été mises en cache ou exécutées

!!! Tip

    Ces rapports sont incroyablement utiles pour optimiser l'allocation des ressources et dépanner les problèmes de performance.

### À retenir

Vous savez comment localiser les sorties du pipeline, examiner les rapports de contrôle qualité et accéder aux métriques d'exécution.

### Et maintenant ?

Apprenez le répertoire de travail et comment Nextflow gère les fichiers intermédiaires.

---

## 4. Explorer le répertoire de travail

Tout comme avec notre exemple Hello World, tout le travail réel se passe dans le répertoire `work/`.

### 4.1. Comprendre la structure du répertoire de travail

Le répertoire de travail contient un sous-répertoire pour chaque tâche qui a été exécutée.
Pour ce pipeline avec 12 tâches, il y aura 12 sous-répertoires de travail.

Listez le répertoire de travail :

```bash
ls -d work/*/*/ | head -5
```

Cela montre les 5 premiers répertoires de tâches.

### 4.2. Inspecter un répertoire de tâche

Prenez l'un des hashs de processus de segmentation de la sortie console (par ex., `[3m/4n5o6p]`) et regardez à l'intérieur :

```bash
ls -la work/3m/4n5o6p*/
```

Vous verrez :

- **Fichiers .command.\*** : Scripts et journaux d'exécution Nextflow (comme avant)
- **Fichiers d'entrée préparés** : Liens symboliques vers les fichiers d'entrée réels
- **Fichiers de sortie** : Masques de segmentation, résultats intermédiaires, etc.

La différence clé par rapport à Hello World :

- Les pipelines réels préparent de gros fichiers d'entrée (images, données de référence)
- Les fichiers de sortie peuvent être assez volumineux (masques de segmentation, images traitées)
- Plusieurs fichiers d'entrée et de sortie par tâche

!!! Tip

    Si un processus échoue, vous pouvez naviguer vers son répertoire de travail, examiner `.command.err` pour les messages d'erreur, et même réexécuter `.command.sh` manuellement pour déboguer le problème.

### 4.3. Nettoyage du répertoire de travail

Le répertoire de travail peut devenir assez volumineux sur plusieurs exécutions de pipeline.
Comme nous l'avons appris dans la Partie 1, vous pouvez utiliser `nextflow clean` pour supprimer les répertoires de travail des anciennes exécutions.

Cependant, pour les pipelines nf-core avec de gros fichiers intermédiaires, il est particulièrement important de nettoyer régulièrement.

### À retenir

Vous comprenez comment les pipelines nf-core organisent leurs répertoires de travail et comment inspecter des tâches individuelles pour le débogage.

### Et maintenant ?

Apprenez le cache Nextflow et comment reprendre les exécutions de pipeline échouées.

---

## 5. Reprendre une exécution de pipeline

L'une des fonctionnalités les plus puissantes de Nextflow est la capacité de reprendre un pipeline à partir du point d'échec.

### 5.1. Le mécanisme de cache

Lorsque vous exécutez un pipeline avec `-resume`, Nextflow :

1. Vérifie le cache pour chaque tâche
2. Si les entrées, le code et les paramètres sont identiques, réutilise le résultat en cache
3. Réexécute uniquement les tâches qui ont changé ou échoué

Ceci est essentiel pour les pipelines de longue durée où des échecs peuvent survenir tard dans l'exécution.

### 5.2. Essayer resume avec molkart

Exécutez à nouveau la même commande, mais ajoutez `-resume` :

```bash
nextflow run ./molkart \
  --input 'data/samplesheet.csv' \
  --mindagap_tilesize 90 \
  --mindagap_boxsize 7 \
  --mindagap_loopnum 100 \
  --clahe_pyramid_tile 368 \
  --segmentation_method "cellpose" \
  --outdir results \
  -resume
```

Vous devriez voir une sortie comme : <!-- TODO: full output -->

```console
executor >  local (0)
[1a/2b3c4d] NFCORE_MOLKART:MOLKART:MINDAGAP_MINDAGAP (mem_only)        [100%] 2 of 2, cached: 2 ✔
[5e/6f7g8h] NFCORE_MOLKART:MOLKART:CLAHE (mem_only)                    [100%] 2 of 2, cached: 2 ✔
[7f/8g9h0i] NFCORE_MOLKART:MOLKART:CREATE_STACK (mem_only)             [100%] 1 of 1, cached: 1 ✔
[9h/0i1j2k] NFCORE_MOLKART:MOLKART:MINDAGAP_DUPLICATEFINDER (mem_only) [100%] 1 of 1, cached: 1 ✔
[2k/3l4m5n] NFCORE_MOLKART:MOLKART:CELLPOSE (mem_only)                 [100%] 1 of 1, cached: 1 ✔
...
```

Remarquez `cached: 2` ou `cached: 1` pour chaque processus - rien n'a été réexécuté !

### 5.3. Quand resume est utile

Resume est particulièrement utile quand :

- Un pipeline échoue en raison de limites de ressources (mémoire insuffisante, limite de temps dépassée)
- Vous devez modifier des processus en aval sans réexécuter les étapes en amont
- Votre connexion réseau tombe pendant le téléchargement de données
- Vous voulez ajouter des sorties supplémentaires sans refaire le calcul

!!! Warning

    Resume ne fonctionne que si vous n'avez pas modifié les données d'entrée, le code du pipeline ou les paramètres.
    Si vous modifiez l'un de ces éléments, Nextflow réexécutera correctement les tâches affectées.

### À retenir

Vous savez comment utiliser `-resume` pour réexécuter efficacement les pipelines sans répéter les tâches réussies.

### Et maintenant ?

Maintenant que vous pouvez exécuter nf-core/molkart avec des données de test, vous êtes prêt·e à apprendre comment le configurer pour vos propres ensembles de données.
