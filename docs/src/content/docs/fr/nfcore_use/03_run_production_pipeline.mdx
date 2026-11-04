# Partie 3 : Exécuter un pipeline de production

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Traduction assistée par IA - [en savoir plus et suggérer des améliorations](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

Dans la [Partie 2](./02_configure_execution.md), vous avez appris à définir des paramètres et à personnaliser la configuration pour nf-core/demo.
Nous allons maintenant appliquer ce que vous avez appris à un vrai pipeline de production, nf-core/rnaseq.

---

## 1. Télécharger et exécuter nf-core/rnaseq

Jusqu'à présent, nous avons utilisé `nf-core/demo`, qui est un pipeline minimal conçu pour la formation.
Nous allons maintenant télécharger un vrai pipeline de production et l'exécuter avec son profil de test.

Le pipeline `nf-core/rnaseq` effectue les étapes principales de l'analyse RNA-seq en masse : contrôle qualité, suppression des adaptateurs, alignement des lectures et quantification au niveau des gènes.
C'est probablement le pipeline nf-core le plus utilisé à ce jour.

### 1.1. Télécharger le pipeline

Exécutez la commande suivante pour le télécharger.

```bash
nextflow pull nf-core/rnaseq
```

??? success "Sortie de la commande"

    ```console
    Checking nf-core/rnaseq ...
     downloaded from https://github.com/nf-core/rnaseq.git - revision: 1f03b53ef7 [master]
    ```

Le pipeline est maintenant mis en cache localement et prêt à être exécuté.

### 1.2. Exécuter le profil de test

Exécutez-le avec le profil de test et Docker :

```bash
nextflow run nf-core/rnaseq -profile test,docker --outdir rnaseq-results
```

??? failure "Sortie de la commande"

    ```console
     N E X T F L O W   ~  version 26.04.4

    Launching `https://github.com/nf-core/rnaseq` [suspicious_dijkstra] revision: 1f03b53ef7 [master]

    ------------------------------------------------------
                                            ,--./,-.
            ___     __   __   __   ___     /,-._.--~'
      |\ | |__  __ /  ` /  \ |__) |__         }  {
      | \| |       \__, \__/ |  \ |___     \`-._,-`-,
                                            `._,._,'
      nf-core/rnaseq 3.26.0
    ------------------------------------------------------
    ...
    [-        ] NFCORE_RNASEQ:RNASEQ:FASTQ_QC_TRIM_FILTER_SETSTRANDEDNESS:FQ_LINT  -
    Plus 47 more processes waiting for tasks…

    Execution cancelled -- Finishing pending tasks before exit
    -[nf-core/rnaseq] Pipeline completed with errors-
    ERROR ~ Error executing process > 'NFCORE_RNASEQ:RNASEQ:FASTQ_QC_TRIM_FILTER_SETSTRANDEDNESS:FQ_LINT (WT_REP2)'

    Caused by:
      Process requirement exceeds available memory -- req: 12 GB; avail: 7.7 GB


    Command executed:

      fq lint \
          --disable-validator P001 \
          SRR6357072_1.fastq.gz SRR6357072_2.fastq.gz > WT_REP2.fq_lint.txt

    Command exit status:
      -

    Command output:
      (empty)

    Work dir:
      /workspaces/training/nfcore-use/work/xx/xxxxxxxxxxxxxxxxxxxxxx

    Container:
      quay.io/biocontainers/fq:0.12.0--h9ee0642_0

    Tip: you can replicate the issue by changing to the process work dir and entering the command `bash .command.run`

     -- Check '.nextflow.log' file for details
    ERROR ~ Pipeline failed. Please refer to troubleshooting docs for common issues: https://nf-co.re/docs/running/troubleshooting
    ```

La ligne clé dans cette erreur est :

```console
Process requirement exceeds available memory -- req: 12 GB; avail: 7.7 GB
```

La machine Codespaces par défaut dispose de 8 Go de RAM, ce qui correspond également à la valeur par défaut typique pour Docker Desktop.
Le pipeline demande 12 Go pour le processus `FQ_LINT` — plus que ce que la machine peut fournir.

Ces 12 Go proviennent du label de ressources `process_low` défini dans `conf/base.config` :

```groovy title="conf/base.config"
withLabel:process_low {
    cpus   = { 2     * task.attempt }
    memory = { 12.GB * task.attempt }
    time   = { 4.h   * task.attempt }
}
```

Une option serait d'utiliser un type de machine plus puissant, mais à des fins de test, nous souhaitons pouvoir exécuter le pipeline sur n'importe quel matériel disponible.
La meilleure approche consiste à remplacer les valeurs de ressources par défaut dans un fichier de configuration personnalisé.

### 1.3. Réexécuter avec une configuration personnalisée

Nous vous fournissons un fichier de configuration personnalisé qui remplace les valeurs de ressources par défaut basées sur les labels.

??? full-code "laptop.config"

    ```groovy title="laptop.config"
    process {
        withLabel: 'process_low' {
            cpus   = 2
            memory = 6.GB
        }
        withLabel: 'process_medium' {
            cpus   = 4
            memory = 6.GB
        }
        withLabel: 'process_high' {
            cpus   = 6
            memory = 6.GB
        }
        withLabel: 'process_high_memory' {
            memory = 6.GB
        }
    }
    ```

La [Partie 2](./02_configure_execution.md) a présenté `withName:` pour cibler un seul processus par son nom.
Ici, nous utilisons `withLabel:` pour cibler en une seule fois tous les processus qui partagent un label.

Ce fichier est déjà présent dans votre répertoire de travail.
Passez-le avec `-c` pour appliquer les valeurs personnalisées :

```bash
nextflow run nf-core/rnaseq -profile test,docker -c laptop.config --outdir rnaseq-results
```

??? success "Sortie de la commande (lancement du pipeline)"

    ```console
     N E X T F L O W   ~  version 26.04.4

    Launching `https://github.com/nf-core/rnaseq` [romantic_faraday] revision: 1f03b53ef7 [master]

    ------------------------------------------------------
                                            ,--./,-.
            ___     __   __   __   ___     /,-._.--~'
      |\ | |__  __ /  ` /  \ |__) |__         }  {
      | \| |       \__, \__/ |  \ |___     \`-._,-`-,
                                            `._,._,'
      nf-core/rnaseq 3.26.0
    ------------------------------------------------------
    ...
    executor >  local (7)
    [xx/xxxxxx] NFCORE_RNASEQ:RNASEQ:FQ_LINT (RAP1_IAA_30M_REP1)   | 3 of 5
    [xx/xxxxxx] NFCORE_RNASEQ:RNASEQ:FASTQC (RAP1_IAA_30M_REP1)    | 2 of 5
    ...
    ```

Le pipeline est maintenant en cours d'exécution, et vous pouvez observer les tâches se terminer une par une.
Sur ce jeu de données de test minimal, l'exécution se terminera en 15 à 20 minutes, en exécutant plus de 200 tâches au total.

Les expériences RNA-seq réelles impliquent généralement des dizaines d'échantillons et s'exécutent pendant des heures, voire des jours.
Nextflow prend en charge les ordonnanceurs HPC (SLURM, PBS, LSF) et les plateformes cloud (AWS, Google Cloud, Azure), qui peuvent réduire considérablement le temps d'exécution en distribuant le travail sur de nombreux nœuds.
La mise en place de ces environnements, cependant, ajoute une complexité significative.

La plateforme Seqera (développée par les créateurs de Nextflow) fournit une interface web pour lancer des pipelines Nextflow sur une infrastructure HPC ou cloud (la vôtre ou une gérée pour vous), avec des capacités de gestion des calculs et des données qui simplifient le processus d'exécution de pipelines à grande échelle.

!!! tip "Astuce"

    Les chercheur·euses académiques peuvent accéder à la plateforme Seqera gratuitement via le [programme académique Seqera](https://seqera.io/academic-program/).

### À retenir

Vous avez téléchargé `nf-core/rnaseq`, découvert le fonctionnement des labels de ressources nf-core, et appris à les remplacer avec un fichier de configuration personnalisé.
Plus important encore, vous avez compris pourquoi l'exécution locale est un point de départ plutôt qu'une destination pour une analyse à grande échelle réelle.

### Et ensuite ?

Vous avez couvert les fondamentaux de l'exécution des pipelines nf-core.
Consultez les [Prochaines étapes](next_steps.md) pour savoir comment continuer.

---

## Résumé

Dans cette partie, vous avez appris à :

- Télécharger et exécuter un pipeline à l'échelle de la production (nf-core/rnaseq), et remplacer ses labels de ressources par défaut
