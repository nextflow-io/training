# Partie 4 : Ajout de tests

Dans la première partie de ce cours, vous avez construit un pipeline d'appel de variants complètement linéaire qui traitait les données de chaque échantillon indépendamment des autres.

Dans la deuxième partie, nous vous avons montré comment utiliser les canaux et les opérateurs de canaux pour implémenter l'appel de variants conjoint avec GATK.

Dans la troisième partie, nous avons modularisé le pipeline.

Dans cette partie de la formation, nous allons vous montrer comment utiliser [**nf-test**](https://www.nf-test.com/), un framework de test qui s'intègre bien avec Nextflow et facilite l'ajout de tests au niveau des modules et au niveau du workflow à votre pipeline. Pour suivre cette partie de la formation, vous devriez avoir complété la Partie 1, la Partie 2 et la Partie 3, ainsi que la [quête annexe nf-test](../../side_quests/nf-test.md), qui couvre les bases de nf-test et pourquoi les tests sont importants.

---

## 0. Échauffement

!!! note

    Assurez-vous d'être dans le répertoire de travail correct :
    `cd /workspaces/training/nf4-science/genomics`

Si vous avez suivi les parties précédentes de ce cours de formation, vous devriez avoir une version fonctionnelle du pipeline de génomique avec la structure de répertoires de modules appropriée.

??? abstract "Contenu du répertoire"

    ```console
    modules/
    ├── gatk
    │   ├── haplotypecaller
    │   │   └── main.nf
    │   └── jointgenotyping
    │       └── main.nf
    └── samtools
        └── index
            └── main.nf
    ```

Ce répertoire modules se trouve dans le répertoire `solutions` si vous en avez besoin.

Nous allons commencer avec le même workflow que dans la Partie 3, que nous vous avons fourni dans le fichier `genomics-4.nf`. Exactement comme pour la [quête annexe nf-test](../../side_quests/nf-test.md), nous allons ajouter plusieurs types de tests différents aux trois processus de ce pipeline, ainsi qu'un test au niveau du workflow.

### 0.1. Vérifier que le workflow s'exécute

Avant de commencer à ajouter des tests, assurez-vous que le workflow s'exécute comme prévu.

```bash
nextflow run genomics-4.nf -resume
```

Cela devrait vous sembler très familier maintenant si vous avez suivi ce cours de formation depuis le début.

??? success "Sortie de la commande"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `genomics-4.nf` [gloomy_poincare] DSL2 - revision: 43203316e0

    executor >  local (7)
    [18/89dfa4] SAMTOOLS_INDEX (1)       | 3 of 3 ✔
    [30/b2522b] GATK_HAPLOTYPECALLER (2) | 3 of 3 ✔
    [a8/d2c189] GATK_JOINTGENOTYPING     | 1 of 1 ✔
    ```

Comme précédemment, il y aura maintenant un répertoire `work` et un répertoire `results_genomics` dans votre répertoire de projet. Nous utiliserons en fait ces résultats plus tard dans nos tests. Mais à partir de maintenant, nous allons utiliser le package `nf-test` pour tester le pipeline.

### 0.2. Initialiser `nf-test`

Comme pour la [quête annexe nf-test](../../side_quests/nf-test.md), nous devons initialiser le package `nf-test`.

```bash
nf-test init
```

??? success "Sortie de la commande"

    ```console
    🚀 nf-test 0.9.3
    https://www.nf-test.com
    (c) 2021 - 2024 Lukas Forer and Sebastian Schoenherr

    Project configured. Configuration is stored in nf-test.config
    ```

??? abstract "Contenu de nf-test.config"

    ```groovy title="nf-test.config"
    config {

        testsDir "tests"
        workDir ".nf-test"
        configFile "tests/nextflow.config"
        profile ""

    }
    ```

Cela crée également un répertoire `tests` contenant une ébauche de fichier de configuration.

### À retenir

Nous sommes maintenant prêts à commencer à écrire des tests pour notre pipeline de génomique.

### Et maintenant ?

Écrire des tests de base qui évaluent si les appels de processus ont réussi et produit les sorties correctes.

---

## 1. Tester un processus pour le succès et la correspondance des sorties

Nous commencerons par tester le processus `SAMTOOLS_INDEX`, qui crée des fichiers d'index pour les fichiers BAM afin de permettre un accès aléatoire efficace. C'est un bon premier cas de test car :

1. Il a une seule entrée bien définie (un fichier BAM)
2. Il produit une sortie prévisible (un fichier d'index BAI)
3. La sortie devrait être identique pour des entrées identiques

### 1.1. Générer une ébauche de fichier de test

Tout d'abord, générez une ébauche de fichier de test :

```bash
nf-test generate process modules/samtools/index/main.nf
```

??? success "Sortie de la commande"

    ```console
    Load source file '/workspaces/training/nf4-science/genomics/modules/samtools/index/main.nf'
    Wrote process test file '/workspaces/training/nf4-science/genomics/tests/modules/samtools/index/main.nf.test'

    SUCCESS: Generated 1 test files.
    ```

Cela crée un fichier dans le même répertoire que `main.nf`.
Vous pouvez naviguer vers le répertoire dans l'explorateur de fichiers et ouvrir le fichier, qui devrait contenir le code suivant :

```groovy title="tests/modules/samtools/index/main.nf.test" linenums="1"
nextflow_process {

    name "Test Process SAMTOOLS_INDEX"
    script "modules/samtools/index/main.nf"
    process "SAMTOOLS_INDEX"

    test("Should run without failures") {

        when {
            params {
                // define parameters here. Example:
                // outdir = "tests/results"
            }
            process {
                """
                // define inputs of the process here. Example:
                // input[0] = file("test-file.txt")
                """
            }
        }

        then {
            assert process.success
            assert snapshot(process.out).match()
        }

    }

}
```

Les assertions de départ devraient être familières de la [quête annexe nf-test](../../side_quests/nf-test.md) :

- `assert process.success` indique que nous attendons que le processus s'exécute avec succès et se termine sans échec.
- `snapshot(process.out).match()` indique que nous attendons que le résultat de l'exécution soit identique au résultat obtenu lors d'une exécution précédente (le cas échéant).
  Nous en discutons plus en détail plus tard.

En utilisant ceci comme point de départ, nous devons ajouter les bonnes entrées de test pour le processus samtools index, et tous les paramètres le cas échéant.

### 1.2. Déplacer le fichier de test et mettre à jour le chemin du script

Avant de commencer à remplir le test, nous devons déplacer le fichier vers son emplacement définitif. Une partie de la raison pour laquelle nous avons ajouté un répertoire pour chaque module est que nous pouvons maintenant inclure les tests dans un répertoire `tests` colocalisé avec le fichier `main.nf` de chaque module. Créez ce répertoire et déplacez le fichier de test là-bas.

```bash
mkdir -p modules/samtools/index/tests
mv tests/modules/samtools/index/main.nf.test modules/samtools/index/tests/
```

Nous pouvons maintenant simplifier la section `script` du fichier de test en un chemin relatif :

=== "Après"

    ```groovy title="modules/samtools/index/tests/main.nf.test" linenums="3" hl_lines="2"
    name "Test Process SAMTOOLS_INDEX"
    script "../main.nf"
    process "SAMTOOLS_INDEX"
    ```

=== "Avant"

    ```groovy title="modules/samtools/index/tests/main.nf.test" linenums="3" hl_lines="2"
    name "Test Process SAMTOOLS_INDEX"
    script "modules/samtools/index/main.nf"
    process "SAMTOOLS_INDEX"
    ```

Cela indique au test où trouver le fichier `main.nf` du module, sans avoir à spécifier le chemin complet.

### 1.3. Fournir des entrées de test pour SAMTOOLS_INDEX

Le fichier d'ébauche inclut un espace réservé que nous devons remplacer par une entrée de test réelle, appropriée à l'entrée de `samtools index`. L'entrée appropriée est un fichier BAM, que nous avons disponible dans le répertoire `data/bam`.

=== "Après"

    ```groovy title="modules/samtools/index/tests/main.nf.test" linenums="14"
    process {
        """
        input[0] = file("${projectDir}/data/bam/reads_son.bam")
        """
    }
    ```

=== "Avant"

    ```groovy title="modules/samtools/index/tests/main.nf.test" linenums="14"
    process {
        """
        // define inputs of the process here. Example:
        // input[0] = file("test-file.txt")
        """
    }
    ```

### 1.4. Nommer le test en fonction de la fonctionnalité

Comme nous l'avons appris précédemment, c'est une bonne pratique de renommer le test avec quelque chose qui a du sens dans le contexte du test.

=== "Après"

    ```groovy title="modules/samtools/index/tests/main.nf.test" linenums="7"
    test("Should index reads_son.bam correctly") {
    ```

    Cela prend une chaîne arbitraire, nous pourrions donc y mettre ce que nous voulons.
    Ici, nous choisissons de faire référence au nom du fichier et à son format.

=== "Avant"

    ```groovy title="modules/samtools/index/tests/main.nf.test" linenums="7"
    test("Should run without failures") {
    ```

### 1.5. Exécuter le test et examiner la sortie

Exécutez le test :

```bash
nf-test test modules/samtools/index/tests/main.nf.test
```

??? success "Sortie de la commande"

    ```console
    🚀 nf-test 0.9.3
    https://www.nf-test.com
    (c) 2021 - 2024 Lukas Forer and Sebastian Schoenherr


    Test Process SAMTOOLS_INDEX

      Test [625e39ee] 'Should index reads_son.bam correctly' PASSED (7.717s)
      Snapshots:
        1 created [Should index reads_son.bam correctly]


    Snapshot Summary:
      1 created

    SUCCESS: Executed 1 tests in 7.727s
    ```

Comme nous l'avons appris précédemment, cela a vérifié l'assertion de base sur le succès du processus et créé un fichier snapshot basé sur la sortie du processus. Nous pouvons voir le contenu du fichier snapshot dans le fichier `tests/modules/samtools/index/tests/main.nf.test.snap` :

```json title="modules/samtools/index/tests/main.nf.test.snap" linenums="1"
{
  "Should index reads_son.bam correctly": {
    "content": [
      {
        "0": [
          [
            "reads_son.bam:md5,af5956d9388ba017944bef276b71d809",
            "reads_son.bam.bai:md5,a2ca7b84998218ee77eff14af8eb8ca2"
          ]
        ]
      }
    ],
    "meta": {
      "nf-test": "0.9.3",
      "nextflow": "25.10.2"
    },
    "timestamp": "2026-01-27T15:09:48.394063389"
  }
}
```

Nous pouvons également exécuter à nouveau le test et voir qu'il réussit, car la sortie est identique au snapshot :

```bash
nf-test test modules/samtools/index/tests/main.nf.test
```

??? success "Sortie de la commande"

    ```console
    🚀 nf-test 0.9.3
    https://www.nf-test.com
    (c) 2021 - 2024 Lukas Forer and Sebastian Schoenherr


    Test Process SAMTOOLS_INDEX

      Test [625e39ee] 'Should index reads_son.bam correctly' PASSED (7.938s)


    SUCCESS: Executed 1 tests in 7.987s
    ```

### 1.6. Ajouter plus de tests à `SAMTOOLS_INDEX`

Parfois, il est utile de tester une gamme de différents fichiers d'entrée pour s'assurer que nous testons une variété de problèmes potentiels. Ajoutez des tests pour les fichiers BAM de la mère et du père dans le trio de nos données de test.

```groovy
    test("Should index reads_mother.bam correctly") {

        when {
            params {
                // define parameters here
            }
            process {
                """
                input[0] = file("${projectDir}/data/bam/reads_mother.bam")
                """
            }
        }

        then {
            assert process.success
            assert snapshot(process.out).match()
        }

    }

    test("Should index reads_father.bam correctly") {

        when {
            params {
                // define parameters here
            }
            process {
                """
                input[0] = file("${projectDir}/data/bam/reads_father.bam")
                """
            }
        }

        then {
            assert process.success
            assert snapshot(process.out).match()
        }

    }
```

Ensuite, vous pouvez exécuter à nouveau le test :

```bash
nf-test test modules/samtools/index/tests/main.nf.test
```

??? success "Sortie de la commande"

    ```console
    🚀 nf-test 0.9.3
    https://www.nf-test.com
    (c) 2021 - 2024 Lukas Forer and Sebastian Schoenherr


    Test Process SAMTOOLS_INDEX

      Test [625e39ee] 'Should index reads_son.bam correctly' PASSED (7.185s)
      Test [a8b28f36] 'Should index reads_mother.bam correctly' PASSED (6.576s)
      Test [c15852a1] 'Should index reads_father.bam correctly' PASSED (6.31s)
      Snapshots:
        2 created [Should index reads_father.bam correctly, Should index reads_mother.bam correctly]


    Snapshot Summary:
      2 created

    SUCCESS: Executed 3 tests in 20.117s
    ```

Notez l'avertissement, faisant référence à l'effet du paramètre `--update-snapshot`.

!!! note

    Ici, nous utilisons des données de test que nous avons utilisées précédemment pour démontrer les sorties scientifiques du pipeline.
    Si nous avions prévu d'utiliser ces tests dans un environnement de production, nous aurions généré des entrées plus petites à des fins de test.

    En général, il est important de garder les tests unitaires aussi légers que possible en utilisant les plus petits morceaux de données nécessaires et suffisants pour évaluer la fonctionnalité du processus, sinon le temps d'exécution total peut s'additionner de manière assez sérieuse.
    Une suite de tests qui prend trop de temps à s'exécuter régulièrement est une suite de tests susceptible d'être ignorée dans l'intérêt de l'opportunité.

### À retenir

Vous avez écrit votre premier test de module pour un processus de génomique, vérifiant que `SAMTOOLS_INDEX` crée correctement des fichiers d'index pour différents fichiers BAM. La suite de tests garantit que :

1. Le processus s'exécute avec succès
2. Les fichiers d'index sont créés
3. Les sorties sont cohérentes entre les exécutions
4. Le processus fonctionne pour tous les fichiers BAM d'échantillons

### Et maintenant ?

Apprendre à écrire des tests pour d'autres processus dans notre workflow de génomique, en utilisant la méthode setup pour gérer les processus chaînés. Nous évaluerons également si les sorties, en particulier nos fichiers VCF, contiennent les appels de variants attendus.

---

## 2. Ajouter des tests à un processus chaîné et tester le contenu

Pour tester `GATK_HAPLOTYPECALLER`, nous devons fournir au processus la sortie de `SAMTOOLS_INDEX` comme entrée. Nous pourrions le faire en exécutant `SAMTOOLS_INDEX`, en récupérant ses sorties et en les stockant avec les données de test pour le workflow. C'est en fait l'approche recommandée pour un pipeline soigné, mais nf-test fournit une approche alternative, utilisant la méthode `setup`.

Avec la méthode setup, nous pouvons déclencher le processus `SAMTOOLS_INDEX` dans le cadre de la configuration du test, puis utiliser sa sortie comme entrée pour `GATK_HAPLOTYPECALLER`. Cela a un coût : nous allons devoir exécuter le processus `SAMTOOLS_INDEX` à chaque fois que nous exécutons le test pour `GATK_HAPLOTYPECALLER`. Cependant, peut-être que nous développons encore le workflow et ne voulons pas pré-générer des données de test que nous pourrions devoir modifier plus tard. Le processus `SAMTOOLS_INDEX` est également très rapide, donc peut-être que les avantages de pré-générer et de stocker ses sorties sont négligeables. Voici comment fonctionne la méthode setup.

### 2.1. Générer et placer le fichier de test

Comme précédemment, nous générons d'abord l'ébauche du fichier :

```bash
nf-test generate process modules/gatk/haplotypecaller/main.nf
```

??? success "Sortie de la commande"

    ```console
    Load source file '/workspaces/training/nf4-science/genomics/modules/gatk/haplotypecaller/main.nf'
    Wrote process test file '/workspaces/training/nf4-science/genomics/tests/modules/gatk/haplotypecaller/main.nf.test'

    SUCCESS: Generated 1 test files.
    ```

Cela produit l'ébauche de test suivante :

```groovy title="tests/modules/gatk/haplotypecaller/main.nf.test" linenums="1"
nextflow_process {

    name "Test Process GATK_HAPLOTYPECALLER"
    script "modules/gatk/haplotypecaller/main.nf"
    process "GATK_HAPLOTYPECALLER"

    test("Should run without failures") {

        when {
            params {
                // define parameters here. Example:
                // outdir = "tests/results"
            }
            process {
                """
                // define inputs of the process here. Example:
                // input[0] = file("test-file.txt")
                """
            }
        }

        then {
            assert process.success
            assert snapshot(process.out).match()
        }

    }

}
```

### 2.2. Déplacer le fichier de test et mettre à jour le chemin du script

Nous créons un répertoire pour le fichier de test colocalisé avec le fichier `main.nf` du module :

```bash
mkdir -p modules/gatk/haplotypecaller/tests
```

Et nous déplaçons le fichier d'ébauche de test là-bas :

```bash
mv tests/modules/gatk/haplotypecaller/main.nf.test modules/gatk/haplotypecaller/tests/
```

Enfin, n'oubliez pas de mettre à jour le chemin du script :

=== "Après"

    ```groovy title="modules/gatk/haplotypecaller/tests/main.nf.test" linenums="3" hl_lines="2"
        name "Test Process GATK_HAPLOTYPECALLER"
        script "../main.nf"
        process "GATK_HAPLOTYPECALLER"
    ```

=== "Avant"

    ```groovy title="modules/gatk/haplotypecaller/tests/main.nf.test" linenums="3" hl_lines="2"
        name "Test Process GATK_HAPLOTYPECALLER"
        script "modules/gatk/haplotypecaller/main.nf"
        process "GATK_HAPLOTYPECALLER"
    ```

### 2.3. Fournir des entrées en utilisant la méthode setup

Nous insérons un bloc `setup` avant le bloc `when`, où nous pouvons déclencher une exécution du processus `SAMTOOLS_INDEX` sur l'un de nos fichiers d'entrée d'origine. N'oubliez pas non plus de changer le nom du test en quelque chose de significatif.

=== "Après"

    ```groovy title="modules/gatk/haplotypecaller/tests/main.nf.test" linenums="7" hl_lines="1-12"
        test("Should call son's haplotype correctly") {

            setup {
                run("SAMTOOLS_INDEX") {
                    script "../../../samtools/index/main.nf"
                    process {
                        """
                        input[0] =  file("${projectDir}/data/bam/reads_son.bam")
                        """
                    }
                }
            }

            when {
    ```

=== "Avant"

    ```groovy title="modules/gatk/haplotypecaller/tests/main.nf.test" linenums="7"  hl_lines="1"
    test("Should run without failures") {

        when {
    ```

Ensuite, nous pouvons faire référence à la sortie de ce processus dans le bloc `when` où nous spécifions les entrées de test :

```groovy title="modules/gatk/haplotypecaller/tests/main.nf.test" linenums="20"
        when {
            params {
                // define parameters here
            }
            process {
                """
                input[0] = SAMTOOLS_INDEX.out
                input[1] = file("${projectDir}/data/ref/ref.fasta")
                input[2] = file("${projectDir}/data/ref/ref.fasta.fai")
                input[3] = file("${projectDir}/data/ref/ref.dict")
                input[4] = file("${projectDir}/data/ref/intervals.bed")
                """
            }
        }
```

Effectuez cette modification et exécutez à nouveau le test :

```bash
nf-test test modules/gatk/haplotypecaller/tests/main.nf.test
```

??? success "Sortie de la commande"

    ```console
    🚀 nf-test 0.9.3
    https://www.nf-test.com
    (c) 2021 - 2024 Lukas Forer and Sebastian Schoenherr


    Test Process GATK_HAPLOTYPECALLER

      Test [c5156c2b] 'Should call son's haplotype correctly' PASSED (40.53s)
      Snapshots:
        1 created [Should call son's haplotype correctly]


    Snapshot Summary:
      1 created

    SUCCESS: Executed 1 tests in 40.555s
    ```

Cela produit également un fichier snapshot comme précédemment.

### 2.4. Exécuter à nouveau et observer l'échec

Fait intéressant, si vous exécutez exactement la même commande à nouveau, cette fois le test échouera.

```bash
nf-test test modules/gatk/haplotypecaller/tests/main.nf.test
```

??? failure "Sortie de la commande"

    ```console
    🚀 nf-test 0.9.3
    https://www.nf-test.com
    (c) 2021 - 2024 Lukas Forer and Sebastian Schoenherr


    Test Process GATK_HAPLOTYPECALLER

      Test [c5156c2b] 'Should call son's haplotype correctly' FAILED (40.123s)

      java.lang.RuntimeException: Different Snapshot:
      [                                                                                           [
          {                                                                                           {
              "0": [                                                                                      "0": [
                  "reads_son.bam.g.vcf:md5,069316cdd4328542ffc6ae247b1dac39"                         |                 "reads_son.bam.g.vcf:md5,005f1a13ee39f11b0fc9bea094850eac"
              ],                                                                                          ],
              "1": [                                                                                      "1": [
                  "reads_son.bam.g.vcf.idx:md5,dc36c18f2afdc546f41e68b2687e9334"                     |                 "reads_son.bam.g.vcf.idx:md5,dbad4b76a4b90c158ffc9c9740764242"
              ],                                                                                          ],
              "idx": [                                                                                    "idx": [
                  "reads_son.bam.g.vcf.idx:md5,dc36c18f2afdc546f41e68b2687e9334"                     |                 "reads_son.bam.g.vcf.idx:md5,dbad4b76a4b90c158ffc9c9740764242"
              ],                                                                                          ],
              "vcf": [                                                                                    "vcf": [
                  "reads_son.bam.g.vcf:md5,069316cdd4328542ffc6ae247b1dac39"                         |                 "reads_son.bam.g.vcf:md5,005f1a13ee39f11b0fc9bea094850eac"
              ]                                                                                           ]
          }                                                                                           }
      ]                                                                                           ]

      Nextflow stdout:

      Nextflow stderr:


        Obsolete snapshots can only be checked if all tests of a file are executed successful.


    FAILURE: Executed 1 tests in 40.156s (1 failed)
    ```

Le message d'erreur vous indique qu'il y avait des différences entre les snapshots pour les deux exécutions ; plus précisément, les valeurs md5sum sont différentes pour les fichiers VCF.

Pourquoi ? Pour faire court, l'outil HaplotypeCaller inclut un horodatage dans l'en-tête VCF qui est différent à chaque fois (par définition).
En conséquence, nous ne pouvons pas simplement nous attendre à ce que les fichiers aient des md5sum identiques même s'ils ont un contenu identique en termes d'appels de variants eux-mêmes.

Comment gérons-nous cela ?

### 2.5. Utiliser une méthode d'assertion de contenu pour vérifier un variant spécifique

Une façon de résoudre le problème est d'utiliser un [type d'assertion différent](https://nf-co.re/docs/contributing/tutorials/nf-test_assertions).
Dans ce cas, nous allons vérifier un contenu spécifique au lieu d'affirmer l'identité.
Plus exactement, nous allons faire lire les lignes du fichier VCF par l'outil et vérifier l'existence de lignes spécifiques.

En pratique, nous remplaçons la deuxième assertion dans le bloc `then` comme suit :

=== "Après"

    ```groovy title="modules/gatk/haplotypecaller/tests/main.nf.test" linenums="35" hl_lines="3 4"
            then {
                assert process.success
                assert path(process.out[0][0]).readLines().contains('#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	reads_son')
                assert path(process.out[0][0]).readLines().contains('20_10037292_10066351	3277	.	G	<NON_REF>	.	.	END=3282	GT:DP:GQ:MIN_DP:PL	0/0:25:72:24:0,72,719')
            }
    ```

=== "Avant"

    ```groovy title="modules/gatk/haplotypecaller/tests/main.nf.test" linenums="35" hl_lines="3"
    then {
        assert process.success
        assert snapshot(process.out).match()
    }
    ```

Ici, nous lisons le contenu complet du fichier de sortie VCF et recherchons une correspondance de contenu, ce qui est acceptable sur un petit fichier de test, mais vous ne voudriez pas faire cela sur un fichier plus volumineux.
Vous pourriez plutôt choisir de lire des lignes spécifiques.

Cette approche nécessite de choisir plus soigneusement ce que nous voulons utiliser comme « signal » à tester.
Le bon côté, c'est qu'elle peut être utilisée pour tester avec une grande précision si un outil d'analyse peut systématiquement identifier des caractéristiques « difficiles » (comme des variants rares) au fur et à mesure de son développement.

### 2.6. Exécuter à nouveau et observer le succès

Une fois que nous avons modifié le test de cette façon, nous pouvons exécuter le test plusieurs fois, et il réussira systématiquement.

```bash
nf-test test modules/gatk/haplotypecaller/tests/main.nf.test
```

??? success "Sortie de la commande"

    ```console
    🚀 nf-test 0.9.3
    https://www.nf-test.com
    (c) 2021 - 2024 Lukas Forer and Sebastian Schoenherr


    Test Process GATK_HAPLOTYPECALLER

      Test [c5156c2b] 'Should call son's haplotype correctly' PASSED (40.53s)


    SUCCESS: Executed 1 tests in 40.555s
    ```

### 2.7. Ajouter plus de tests

Ajoutez des tests similaires pour les échantillons de la mère et du père :

```groovy title="modules/gatk/haplotypecaller/tests/main.nf.test" linenums="43"
    test("Should call mother's haplotype correctly") {

        setup {
            run("SAMTOOLS_INDEX") {
                script "../../../samtools/index/main.nf"
                process {
                    """
                    input[0] =  file("${projectDir}/data/bam/reads_mother.bam")
                    """
                }
            }
        }

        when {
            params {
                // define parameters here
            }
            process {
                """
                input[0] = SAMTOOLS_INDEX.out
                input[1] = file("${projectDir}/data/ref/ref.fasta")
                input[2] = file("${projectDir}/data/ref/ref.fasta.fai")
                input[3] = file("${projectDir}/data/ref/ref.dict")
                input[4] = file("${projectDir}/data/ref/intervals.bed")
                """
            }
        }

        then {
            assert process.success
            assert path(process.out[0][0]).readLines().contains('#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	reads_mother')
            assert path(process.out[0][0]).readLines().contains('20_10037292_10066351	3277	.	G	<NON_REF>	.	.	END=3278	GT:DP:GQ:MIN_DP:PL	0/0:38:99:37:0,102,1530')
        }
    }

    test("Should call father's haplotype correctly") {

        setup {
            run("SAMTOOLS_INDEX") {
                script "../../../samtools/index/main.nf"
                process {
                    """
                    input[0] =  file("${projectDir}/data/bam/reads_father.bam")
                    """
                }
            }
        }

        when {
            params {
                // define parameters here
            }
            process {
                """
                input[0] = SAMTOOLS_INDEX.out
                input[1] = file("${projectDir}/data/ref/ref.fasta")
                input[2] = file("${projectDir}/data/ref/ref.fasta.fai")
                input[3] = file("${projectDir}/data/ref/ref.dict")
                input[4] = file("${projectDir}/data/ref/intervals.bed")
                """
            }
        }

        then {
            assert process.success
            assert path(process.out[0][0]).readLines().contains('#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	reads_father')
            assert path(process.out[0][0]).readLines().contains('20_10037292_10066351	3277	.	G	<NON_REF>	.	.	END=3281	GT:DP:GQ:MIN_DP:PL	0/0:44:99:42:0,120,1800')
        }
    }
```

### 2.8. Exécuter la commande de test

```bash
nf-test test modules/gatk/haplotypecaller/tests/main.nf.test
```

??? success "Sortie de la commande"

    ```console
    🚀 nf-test 0.9.3
    https://www.nf-test.com
    (c) 2021 - 2024 Lukas Forer and Sebastian Schoenherr


    Test Process GATK_HAPLOTYPECALLER

      Test [c5156c2b] 'Should call son's haplotype correctly' PASSED (40.53s)
      Test [10de94a8] 'Should call mother's haplotype correctly' PASSED (41.47s)
      Test [c0386fc7] 'Should call father's haplotype correctly' PASSED (45.556s)


    SUCCESS: Executed 3 tests in 127.586s
    ```

Cela complète le plan de test de base pour cette deuxième étape dans le pipeline. Passons au troisième et dernier test au niveau du module !

### À retenir

Vous avez appris comment :

1. Tester des processus qui dépendent des sorties d'autres processus
2. Vérifier des variants génomiques spécifiques dans les fichiers de sortie VCF
3. Gérer les sorties non déterministes en vérifiant un contenu spécifique
4. Tester l'appel de variants sur plusieurs échantillons

### Et maintenant ?

Apprendre à écrire des tests qui utilisent des données de test pré-générées pour l'étape de génotypage conjoint.

---

## 3. Utiliser des données de test pré-générées

Pour l'étape de génotypage conjoint, nous utiliserons une approche différente - l'utilisation de données de test pré-générées. Ceci est souvent préférable pour :

1. Les processus complexes avec plusieurs dépendances
2. Les processus qui prennent beaucoup de temps à s'exécuter
3. Les processus qui font partie d'un pipeline stable de production

### 3.1. Générer des données de test

Inspectez les résultats que nous avons générés au début de cette section :

```bash
tree results_genomics/
```

```console title="Contenu du répertoire de résultats"
results_genomics/
├── family_trio.joint.vcf
├── family_trio.joint.vcf.idx
├── gvcf
│   ├── reads_father.bam.g.vcf -> /workspaces/training/nf4-science/genomics/work/30/b2522b83c63baff8c3cf75704512a2/reads_father.bam.g.vcf
│   ├── reads_father.bam.g.vcf.idx -> /workspaces/training/nf4-science/genomics/work/30/b2522b83c63baff8c3cf75704512a2/reads_father.bam.g.vcf.idx
│   ├── reads_mother.bam.g.vcf -> /workspaces/training/nf4-science/genomics/work/f6/be2efa58e625d08cf8d0da1d0e9f09/reads_mother.bam.g.vcf
│   ├── reads_mother.bam.g.vcf.idx -> /workspaces/training/nf4-science/genomics/work/f6/be2efa58e625d08cf8d0da1d0e9f09/reads_mother.bam.g.vcf.idx
│   ├── reads_son.bam.g.vcf -> /workspaces/training/nf4-science/genomics/work/fe/2f22d56aa16ed45f8bc419312894f6/reads_son.bam.g.vcf
│   └── reads_son.bam.g.vcf.idx -> /workspaces/training/nf4-science/genomics/work/fe/2f22d56aa16ed45f8bc419312894f6/reads_son.bam.g.vcf.idx
└── indexed_bam
    ├── reads_father.bam -> /workspaces/training/nf4-science/genomics/work/42/a3bf19dbfaf1f3672b16a5d5e6a8be/reads_father.bam
    ├── reads_father.bam.bai -> /workspaces/training/nf4-science/genomics/work/cf/289c2d264f496d60a69e3e9ba6463e/reads_father.bam.bai
    ├── reads_mother.bam -> /workspaces/training/nf4-science/genomics/work/af/f31a6ade82cc0cf853c4f61c8bc473/reads_mother.bam
    ├── reads_mother.bam.bai -> /workspaces/training/nf4-science/genomics/work/18/89dfa40a3def17e45421e54431a126/reads_mother.bam.bai
    ├── reads_son.bam -> /workspaces/training/nf4-science/genomics/work/9f/9615dd553d6f13d8bec4f006ac395f/reads_son.bam
    └── reads_son.bam.bai -> /workspaces/training/nf4-science/genomics/work/4d/cb384a97db5687cc9daab002017c7c/reads_son.bam.bai

2 directories, 14 files
```

L'étape de génotypage conjoint nécessite les fichiers VCF produits par les étapes de l'appeleur d'haplotypes comme entrées, ainsi que les indices. Copions donc les résultats que nous avons dans le répertoire de tests du module `jointgenotyping`.

```bash
mkdir -p modules/gatk/jointgenotyping/tests/inputs/
cp results_genomics/gvcf/*.g.vcf results_genomics/gvcf/*.g.vcf.idx modules/gatk/jointgenotyping/tests/inputs/
```

Nous pouvons maintenant utiliser ces fichiers comme entrées pour le test que nous allons écrire pour l'étape de génotypage conjoint.

### 3.2. Générer l'ébauche du fichier de test

Comme précédemment, nous générons d'abord l'ébauche du fichier :

```bash
nf-test generate process modules/gatk/jointgenotyping/main.nf
```

??? success "Sortie de la commande"

    ```console
    Load source file '/workspaces/training/nf4-science/genomics/modules/gatk/jointgenotyping/main.nf'
    Wrote process test file '/workspaces/training/nf4-science/genomics/tests/modules/gatk/jointgenotyping/main.nf.test'

    SUCCESS: Generated 1 test files.
    ```

Cela produit l'ébauche de test suivante :

```groovy title="tests/modules/gatk/jointgenotyping/main.nf.test" linenums="1"
nextflow_process {

    name "Test Process GATK_JOINTGENOTYPING"
    script "modules/gatk/jointgenotyping/main.nf"
    process "GATK_JOINTGENOTYPING"

    test("Should run without failures") {

        when {
            params {
                // define parameters here. Example:
                // outdir = "tests/results"
            }
            process {
                """
                // define inputs of the process here. Example:
                // input[0] = file("test-file.txt")
                """
            }
        }

        then {
            assert process.success
            assert snapshot(process.out).match()
        }

    }

}
```

### 3.3. Déplacer le fichier de test et mettre à jour le chemin du script

Cette fois, nous avons déjà un répertoire pour les tests colocalisé avec le fichier `main.nf` du module, nous pouvons donc déplacer le fichier d'ébauche de test là-bas :

```bash
mv tests/modules/gatk/jointgenotyping/main.nf.test modules/gatk/jointgenotyping/tests/
```

Et n'oubliez pas de mettre à jour le chemin du script :

=== "Après"

    ```groovy title="modules/gatk/jointgenotyping/tests/main.nf.test" linenums="3" hl_lines="2"
    name "Test Process GATK_JOINTGENOTYPING"
    script "../main.nf"
    process "GATK_JOINTGENOTYPING"
    ```

=== "Avant"

    ```groovy title="modules/gatk/jointgenotyping/tests/main.nf.test" linenums="3" hl_lines="2"
    name "Test Process GATK_JOINTGENOTYPING"
    script "modules/gatk/jointgenotyping/main.nf"
    process "GATK_JOINTGENOTYPING"
    ```

### 3.4. Fournir les entrées

Remplissez les entrées en fonction des définitions d'entrée du processus et renommez le test en conséquence :

```groovy title="modules/gatk/jointgenotyping/tests/main.nf.test" linenums="7"
    test("Should call trio's joint genotype correctly") {

        when {
            params {
                // define parameters here
            }
            process {
                """
                input[0] = [
                    file("${projectDir}/modules/gatk/jointgenotyping/tests/inputs/reads_father.bam.g.vcf"),
                    file("${projectDir}/modules/gatk/jointgenotyping/tests/inputs/reads_mother.bam.g.vcf"),
                    file("${projectDir}/modules/gatk/jointgenotyping/tests/inputs/reads_son.bam.g.vcf")
                ]
                input[1] = [
                    file("${projectDir}/modules/gatk/jointgenotyping/tests/inputs/reads_father.bam.g.vcf.idx"),
                    file("${projectDir}/modules/gatk/jointgenotyping/tests/inputs/reads_mother.bam.g.vcf.idx"),
                    file("${projectDir}/modules/gatk/jointgenotyping/tests/inputs/reads_son.bam.g.vcf.idx")
                ]
                input[2] = file("${projectDir}/data/ref/intervals.bed")
                input[3] = "family_trio"
                input[4] = file("${projectDir}/data/ref/ref.fasta")
                input[5] = file("${projectDir}/data/ref/ref.fasta.fai")
                input[6] = file("${projectDir}/data/ref/ref.dict")
                """
            }
        }
```

### 3.5. Utiliser des assertions de contenu

La sortie de l'étape de génotypage conjoint est un autre fichier VCF, nous allons donc utiliser à nouveau une assertion de contenu.

```groovy title="modules/gatk/jointgenotyping/tests/main.nf.test" linenums="25"
    then {
        assert process.success
        assert path(process.out[0][0]).readLines().contains('#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	reads_father	reads_mother	reads_son')
        assert path(process.out[0][0]).readLines().contains('20_10037292_10066351	3480	.	C	CT	1625.89	.	AC=5;AF=0.833;AN=6;BaseQRankSum=0.220;DP=85;ExcessHet=0.0000;FS=2.476;MLEAC=5;MLEAF=0.833;MQ=60.00;MQRankSum=0.00;QD=21.68;ReadPosRankSum=-1.147e+00;SOR=0.487	GT:AD:DP:GQ:PL	0/1:15,16:31:99:367,0,375	1/1:0,18:18:54:517,54,0	1/1:0,26:26:78:756,78,0')
    }
```

En vérifiant le contenu d'un variant spécifique dans le fichier de sortie, ce test vérifie que :

1. Le processus de génotypage conjoint s'exécute avec succès
2. Le VCF de sortie contient les trois échantillons dans le bon ordre
3. Un variant spécifique est appelé correctement avec :
   - Des génotypes précis pour chaque échantillon (0/1 pour le père, 1/1 pour la mère et le fils)
   - Des profondeurs de lecture et des qualités de génotype correctes
   - Des statistiques au niveau de la population comme la fréquence allélique (AF=0.833)

Nous n'avons pas fait de snapshot du fichier entier, mais en vérifiant un variant spécifique, nous pouvons être confiants que le processus de génotypage conjoint fonctionne comme prévu.

### 3.6. Exécuter le test

```bash
nf-test test modules/gatk/jointgenotyping/tests/main.nf.test
```

??? success "Sortie de la commande"

    ```console
    🚀 nf-test 0.9.3
    https://www.nf-test.com
    (c) 2021 - 2024 Lukas Forer and Sebastian Schoenherr


    Test Process GATK_JOINTGENOTYPING

      Test [ac2067de] 'Should call trio's joint genotype correctly' PASSED (53.827s)


    SUCCESS: Executed 1 tests in 53.837s
    ```

Le test réussit, vérifiant que notre processus de génotypage conjoint :

1. Combine correctement les VCF d'échantillons individuels
2. Effectue l'appel de variants conjoint
3. Produit un VCF multi-échantillons avec des appels de génotype cohérents entre les exécutions

### À retenir

Vous savez comment :

- Utiliser des résultats précédemment générés comme entrées pour les tests
- Écrire des tests en utilisant des données de test pré-générées

### Et maintenant ?

Ajouter un test au niveau du workflow pour vérifier que l'ensemble du pipeline d'appel de variants fonctionne de bout en bout.

---

## 4. Ajouter un test au niveau du workflow

Nous allons maintenant tester le pipeline complet d'appel de variants, des fichiers BAM aux génotypes conjoints. Cela vérifie que :

1. Tous les processus fonctionnent correctement ensemble
2. Les données circulent correctement entre les étapes
3. Les appels de variants finaux sont cohérents

### 4.1. Générer le test de workflow

Générez un fichier de test pour le pipeline complet :

```bash
nf-test generate pipeline genomics-4.nf
```

??? success "Sortie de la commande"

    ```console
    Load source file '/workspaces/training/nf4-science/genomics/genomics-4.nf'
    Wrote pipeline test file '/workspaces/training/nf4-science/genomics/tests/genomics-4.nf.test'

    SUCCESS: Generated 1 test files.
    ```

Cela crée une ébauche de test de base :

```groovy title="tests/genomics-4.nf.test" linenums="1"
nextflow_pipeline {

    name "Test Workflow genomics-4.nf"
    script "genomics-4.nf"

    test("Should run without failures") {

        when {
            params {
                // define parameters here. Example:
                // outdir = "tests/results"
            }
        }

        then {
            assert workflow.success
        }

    }

}
```

Corrigez simplement le nom en quelque chose de significatif (vous verrez pourquoi c'est utile sous peu).

=== "Après"

    ```groovy title="tests/genomics-4.nf.test" linenums="6" hl_lines="1"
        test("Should run the pipeline without failures") {
    ```

=== "Avant"

    ```groovy title="tests/genomics-4.nf.test" linenums="6" hl_lines="1"
        test("Should run without failures") {
    ```

!!! note

    Dans ce cas, le fichier de test peut rester là où `nf-test` l'a créé.

### 4.2. Spécifier les paramètres d'entrée

Nous devons encore spécifier les entrées, ce qui se fait de manière légèrement différente au niveau du workflow par rapport aux tests au niveau du module.
Il existe plusieurs façons de le faire, notamment en spécifiant un profil.
Cependant, une façon plus simple est de configurer un bloc `params {}` dans le fichier `nextflow.config` que `nf-test init` a créé à l'origine dans le répertoire `tests`.

```groovy title="tests/nextflow.config" linenums="1"
/*
========================================================================================
    Fichier de configuration Nextflow pour l'exécution des tests
========================================================================================
*/

// Répertoire de sortie pour les sorties du workflow
outputDir = 'results_genomics'

/*
 * Paramètres du pipeline
 */

params {
    // Entrée principale (fichier de fichiers d'entrée, un par ligne)
    reads_bam = "${projectDir}/data/sample_bams.txt"

    // Fichiers accessoires
    reference = "${projectDir}/data/ref/ref.fasta"
    reference_index = "${projectDir}/data/ref/ref.fasta.fai"
    reference_dict = "${projectDir}/data/ref/ref.dict"
    intervals = "${projectDir}/data/ref/intervals.bed"

    // Nom de base pour le fichier de sortie final
    cohort_name = "family_trio"
}
```

Lorsque nous exécuterons le test, `nf-test` récupérera ce fichier de configuration et importera les entrées en conséquence.

### 4.3. Exécuter le test de workflow

```bash
nf-test test tests/genomics-4.nf.test
```

??? success "Sortie de la commande"

    ```console
    🚀 nf-test 0.9.3
    https://www.nf-test.com
    (c) 2021 - 2024 Lukas Forer and Sebastian Schoenherr


    Test Workflow genomics-4.nf

      Test [1b4c6936] 'Should run the pipeline without failures' PASSED (171.019s)


    SUCCESS: Executed 1 tests in 171.056s
    ```

Le test réussit, confirmant que notre pipeline complet d'appel de variants :

1. Traite avec succès tous les échantillons
2. Enchaîne correctement toutes les étapes

### 4.4. Exécuter TOUS les tests

nf-test a encore un tour dans son sac. Nous pouvons exécuter tous les tests en une seule fois ! Modifiez le fichier `nf-test.config` pour que nf-test recherche dans tous les répertoires les fichiers nf-test. Vous pouvez le faire en modifiant le paramètre `testsDir` :

=== "Après"

    ```groovy title="nf-test.config" linenums="1" hl_lines="3"
    config {

        testsDir "."
        workDir ".nf-test"
        configFile "tests/nextflow.config"
        profile ""

    }
    ```

=== "Avant"

    ```groovy title="nf-test.config" linenums="1" hl_lines="3"
    config {

        testsDir "tests"
        workDir ".nf-test"
        configFile "tests/nextflow.config"
        profile ""

    }
    ```

Maintenant, nous pouvons simplement exécuter nf-test et il exécutera _chaque test unique_ dans notre dépôt :

```bash
nf-test test
```

??? success "Sortie de la commande"

    ```console
    🚀 nf-test 0.9.3
    https://www.nf-test.com
    (c) 2021 - 2024 Lukas Forer and Sebastian Schoenherr


    Test Process GATK_HAPLOTYPECALLER

      Test [c5156c2b] 'Should call son's haplotype correctly' PASSED (39.947s)
      Test [10de94a8] 'Should call mother's haplotype correctly' PASSED (43.17s)
      Test [c0386fc7] 'Should call father's haplotype correctly' PASSED (44.244s)

    Test Process GATK_JOINTGENOTYPING

      Test [ac2067de] 'Should call trio's joint genotype correctly' PASSED (61.129s)

    Test Process SAMTOOLS_INDEX

      Test [625e39ee] 'Should index reads_son.bam correctly' PASSED (8.671s)
      Test [a8b28f36] 'Should index reads_mother.bam correctly' PASSED (8.518s)
      Test [c15852a1] 'Should index reads_father.bam correctly' PASSED (5.378s)

    Test Workflow genomics-4.nf

      Test [1b4c6936] 'Should run the pipeline without failures' PASSED (169.714s)


    SUCCESS: Executed 8 tests in 380.801s
    ```

8 tests en 1 commande ! Nous avons passé beaucoup de temps à configurer de nombreux tests, mais quand il s'agit de les exécuter, c'était très rapide et facile. Vous pouvez voir à quel point c'est utile lors de la maintenance d'un grand pipeline, qui pourrait inclure des centaines d'éléments différents. Nous passons du temps à écrire des tests une fois pour pouvoir gagner du temps en les exécutant plusieurs fois.

De plus, nous pouvons automatiser cela ! Imaginez que des tests s'exécutent à chaque fois que vous ou un collègue essayez d'ajouter du nouveau code. C'est ainsi que nous nous assurons que nos pipelines maintiennent un niveau de qualité élevé.

## À retenir

Vous savez maintenant comment écrire et exécuter plusieurs types de tests pour votre pipeline de génomique en utilisant nf-test. Ce framework de test aide à garantir que votre workflow d'appel de variants produit des résultats cohérents et fiables dans différents environnements et au fur et à mesure que vous apportez des modifications au code.

Vous avez appris à tester des composants critiques tels que :

- Le processus `SAMTOOLS_INDEX` qui prépare les fichiers BAM pour l'appel de variants
- Le processus `GATK_HAPLOTYPECALLER` qui identifie les variants dans des échantillons individuels
- Le processus `GATK_JOINTGENOTYPING` qui combine les appels de variants à travers une cohorte

Vous avez également mis en œuvre différentes stratégies de test spécifiques aux données de génomique :

- Vérifier que les fichiers VCF contiennent les appels de variants attendus malgré des éléments non déterministes comme les horodatages
- Tester avec un jeu de données de trio familial pour assurer une identification appropriée des variants à travers des échantillons apparentés
- Vérifier des coordonnées génomiques spécifiques et des informations de variants dans vos fichiers de sortie

Ces compétences en tests sont essentielles pour développer des pipelines de bioinformatique robustes qui peuvent traiter de manière fiable des données génomiques et produire des appels de variants précis. Au fur et à mesure que vous continuez à travailler avec Nextflow pour l'analyse génomique, cette base de tests vous aidera à maintenir un code de haute qualité qui produit des résultats scientifiques dignes de confiance.
