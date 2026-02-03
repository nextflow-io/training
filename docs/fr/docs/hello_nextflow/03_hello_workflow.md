# Partie 3 : Hello Workflow

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Traduction assistée par IA - [en savoir plus et suggérer des améliorations](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

<!--
<div class="video-wrapper">
  <iframe width="560" height="315" src="https://www.youtube.com/embed/zJP7cUYPEbA?si=Irl9nAQniDyICp2b&amp;list=PLPZ8WHdZGxmXiHf8B26oB_fTfoKQdhlik" title="YouTube video player" frameborder="0" allow="accelerometer; autoplay; clipboard-write; encrypted-media; gyroscope; picture-in-picture; web-share" referrerpolicy="strict-origin-when-cross-origin" allowfullscreen></iframe>
</div>

/// caption
:fontawesome-brands-youtube:{ .youtube } Voir [la playlist complète](https://www.youtube.com/playlist?list=PLPZ8WHdZGxmXiHf8B26oB_fTfoKQdhlik) sur la chaîne YouTube de Nextflow.

:green_book: La transcription de la vidéo est disponible [ici](./transcripts/03_hello_workflow.md).
///
-->

La plupart des workflows réels impliquent plus d'une étape.
Dans ce module de formation, vous apprendrez à connecter des processus ensemble dans un workflow multi-étapes.

Cela vous enseignera la méthode Nextflow pour accomplir les tâches suivantes :

1. Faire circuler les données d'un processus à l'autre
2. Collecter les sorties de plusieurs appels de processus dans un seul appel de processus
3. Passer plus d'une entrée à un processus
4. Gérer plusieurs sorties provenant d'un processus

Pour illustrer, nous continuerons à construire sur l'exemple Hello World indépendant du domaine des Parties 1 et 2.
Cette fois, nous allons apporter les modifications suivantes à notre workflow pour mieux refléter la façon dont les gens construisent de véritables workflows :

1. Ajouter une deuxième étape qui convertit le message de bienvenue en majuscules.
2. Ajouter une troisième étape qui collecte tous les messages transformés et les écrit dans un seul fichier.
3. Ajouter un paramètre pour nommer le fichier de sortie final et le passer comme entrée secondaire à l'étape de collecte.
4. Faire en sorte que l'étape de collecte rapporte également une statistique simple sur ce qui a été traité.

??? info "Comment commencer à partir de cette section"

    Cette section du cours suppose que vous avez complété les Parties 1-2 du cours [Hello Nextflow](./index.md), mais si vous êtes à l'aise avec les bases couvertes dans ces sections, vous pouvez commencer ici sans rien faire de spécial.

---

## 0. Échauffement : Exécuter `hello-workflow.nf`

Nous allons utiliser le script de workflow `hello-workflow.nf` comme point de départ.
Il est équivalent au script produit en travaillant à travers la Partie 2 de ce cours de formation, sauf que nous avons supprimé les instructions `view()` et changé la destination de sortie :

```groovy title="hello-workflow.nf" linenums="37" hl_lines="3"
output {
    first_output {
        path 'hello_workflow'
        mode 'copy'
    }
}
```

Juste pour s'assurer que tout fonctionne, exécutez le script une fois avant d'effectuer des modifications :

```bash
nextflow run hello-workflow.nf
```

??? success "Sortie de la commande"

    ```console
     N E X T F L O W   ~  version 25.10.2

    Launching `hello-workflow.nf` [admiring_lamarr] DSL2 - revision: 4d4053520d

    executor >  local (3)
    [b1/5826b5] process > sayHello (2) [100%] 3 of 3 ✔
    ```

Comme précédemment, vous trouverez les fichiers de sortie à l'emplacement spécifié dans le bloc `output`.
Pour ce chapitre, c'est sous `results/hello_workflow/`.

??? abstract "Contenu du répertoire"

    ```console
    results/hello_workflow
    ├── Bonjour-output.txt
    ├── Hello-output.txt
    └── Holà-output.txt
    ```

Si cela a fonctionné pour vous, vous êtes prêt·e à apprendre comment assembler un workflow multi-étapes.

---

## 1. Ajouter une deuxième étape au workflow

Nous allons ajouter une étape pour convertir chaque message de bienvenue en majuscules.

<figure class="excalidraw">
--8<-- "docs/en/docs/hello_nextflow/img/hello-multistep.svg"
</figure>

Pour ce faire, nous devons accomplir trois choses :

- Définir la commande que nous allons utiliser pour effectuer la conversion en majuscules.
- Écrire un nouveau processus qui encapsule la commande de mise en majuscules.
- Appeler le nouveau processus dans le bloc workflow et le configurer pour prendre la sortie du processus `sayHello()` comme entrée.

### 1.1. Définir la commande de mise en majuscules et la tester dans le terminal

Pour effectuer la conversion des messages de bienvenue en majuscules, nous allons utiliser un outil UNIX classique appelé `tr` pour « remplacement de texte », avec la syntaxe suivante :

```bash title="Syntaxe"
tr '[a-z]' '[A-Z]'
```

C'est un one-liner de remplacement de texte très naïf qui ne tient pas compte des lettres accentuées, donc par exemple « Holà » deviendra « HOLà », mais il fera un travail suffisant pour démontrer les concepts Nextflow et c'est ce qui compte.

Pour le tester, nous pouvons exécuter la commande `echo 'Hello World'` et rediriger sa sortie vers la commande `tr` :

```bash
echo 'Hello World' | tr '[a-z]' '[A-Z]' > UPPER-output.txt
```

La sortie est un fichier texte appelé `UPPER-output.txt` qui contient la version en majuscules de la chaîne `Hello World`.

??? abstract "Contenu du fichier"

    ```console title="UPPER-output.txt"
    HELLO WORLD
    ```

C'est essentiellement ce que nous allons essayer de faire avec notre workflow.

### 1.2. Écrire l'étape de mise en majuscules comme un processus Nextflow

Nous pouvons modeler notre nouveau processus sur le premier, puisque nous voulons utiliser tous les mêmes composants.

Ajoutez la définition de processus suivante au script de workflow, juste en dessous de la première :

```groovy title="hello-workflow.nf" linenums="20"
/*
 * Utilise un outil de remplacement de texte pour convertir la salutation en majuscules
 */
process convertToUpper {

    input:
    path input_file

    output:
    path "UPPER-${input_file}"

    script:
    """
    cat '$input_file' | tr '[a-z]' '[A-Z]' > 'UPPER-${input_file}'
    """
}
```

Dans celui-ci, nous composons le deuxième nom de fichier de sortie basé sur le nom du fichier d'entrée, de manière similaire à ce que nous avions fait initialement pour la sortie du premier processus.

### 1.3. Ajouter un appel au nouveau processus dans le bloc workflow

Maintenant, nous devons dire à Nextflow d'appeler effectivement le processus que nous venons de définir.

Dans le bloc workflow, effectuez la modification de code suivante :

=== "Après"

    ```groovy title="hello-workflow.nf" linenums="44" hl_lines="10-11"
    workflow {

        main:
        // créer un canal pour les entrées depuis un fichier CSV
        greeting_ch = channel.fromPath(params.input)
                            .splitCsv()
                            .map { line -> line[0] }
        // émettre une salutation
        sayHello(greeting_ch)
        // convertir la salutation en majuscules
        convertToUpper()

        publish:
        first_output = sayHello.out
    }
    ```

=== "Avant"

    ```groovy title="hello-workflow.nf" linenums="44"
    workflow {

        main:
        // créer un canal pour les entrées depuis un fichier CSV
        greeting_ch = channel.fromPath(params.input)
                            .splitCsv()
                            .map { line -> line[0] }
        // émettre une salutation
        sayHello(greeting_ch)

        publish:
        first_output = sayHello.out
    }
    ```

Ceci n'est pas encore fonctionnel car nous n'avons pas spécifié ce qui doit être passé en entrée au processus `convertToUpper()`.

### 1.4. Passer la sortie du premier processus au deuxième processus

Maintenant, nous devons faire en sorte que la sortie du processus `sayHello()` circule vers le processus `convertToUpper()`.

De manière pratique, Nextflow empaquette automatiquement la sortie d'un processus dans un canal appelé `<process>.out`.
Ainsi, la sortie du processus `sayHello` est un canal appelé `sayHello.out`, que nous pouvons brancher directement dans l'appel à `convertToUpper()`.

Dans le bloc workflow, effectuez la modification de code suivante :

=== "Après"

    ```groovy title="hello-workflow.nf" linenums="53" hl_lines="2"
        // convertir la salutation en majuscules
        convertToUpper(sayHello.out)
    ```

=== "Avant"

    ```groovy title="hello-workflow.nf" linenums="53" hl_lines="2"
        // convertir la salutation en majuscules
        convertToUpper()
    ```

Pour un cas simple comme celui-ci (une sortie vers une entrée), c'est tout ce que nous devons faire pour connecter deux processus !

### 1.5. Configurer la publication des sorties du workflow

Enfin, mettons à jour les sorties du workflow pour publier également les résultats du deuxième processus.

#### 1.5.1. Mettre à jour la section `publish:` du bloc `workflow`

Dans le bloc `workflow`, effectuez la modification de code suivante :

=== "Après"

    ```groovy title="hello-workflow.nf" linenums="56" hl_lines="3"
        publish:
        first_output = sayHello.out
        uppercased = convertToUpper.out
    }
    ```

=== "Avant"

    ```groovy title="hello-workflow.nf" linenums="56"
        publish:
        first_output = sayHello.out
    }
    ```

La logique est la même que précédemment.

#### 1.5.2. Mettre à jour le bloc `output`

Dans le bloc `output`, effectuez la modification de code suivante :

=== "Après"

    ```groovy title="hello-workflow.nf" linenums="61" hl_lines="6-9"
    output {
        first_output {
            path 'hello_workflow'
            mode 'copy'
        }
        uppercased {
            path 'hello_workflow'
            mode 'copy'
        }
    }
    ```

=== "Avant"

    ```groovy title="hello-workflow.nf" linenums="61"
    output {
        first_output {
            path 'hello_workflow'
            mode 'copy'
        }
    }
    ```

Encore une fois, la logique est la même qu'avant.

Cela vous montre que vous pouvez contrôler les paramètres de sortie à un niveau très granulaire, pour chaque sortie individuelle.
N'hésitez pas à essayer de changer les chemins ou le mode de publication pour l'un des processus pour voir ce qui se passe.

Bien sûr, cela signifie que nous répétons certaines informations ici, ce qui pourrait devenir gênant si nous voulions mettre à jour l'emplacement pour toutes les sorties de la même manière.
Plus tard dans le cours, vous apprendrez comment configurer ces paramètres pour plusieurs sorties de manière structurée.

### 1.6. Exécuter le workflow avec `-resume`

Testons ceci en utilisant le flag `-resume`, puisque nous avons déjà exécuté la première étape du workflow avec succès.

```bash
nextflow run hello-workflow.nf -resume
```

??? success "Sortie de la commande"

    ```console
     N E X T F L O W   ~  version 25.10.2

    Launching `hello-workflow.nf` [high_cantor] DSL2 - revision: d746983511

    executor >  local (3)
    [ab/816321] process > sayHello (3)       [100%] 3 of 3, cached: 3 ✔
    [e0/ecf81b] process > convertToUpper (3) [100%] 3 of 3 ✔
    ```

Il y a maintenant une ligne supplémentaire dans la sortie de la console qui correspond au nouveau processus que nous venons d'ajouter.

Vous trouverez les sorties dans le répertoire `results/hello_workflow` comme défini dans le bloc `output`.

??? abstract "Contenu du répertoire"

    ```console
    results/hello_workflow/
    ├── Bonjour-output.txt
    ├── Hello-output.txt
    ├── Holà-output.txt
    ├── UPPER-Bonjour-output.txt
    ├── UPPER-Hello-output.txt
    └── UPPER-Holà-output.txt
    ```

C'est pratique ! Mais cela vaut quand même la peine de jeter un œil à l'intérieur du répertoire de travail d'un des appels au deuxième processus.

??? abstract "Contenu du répertoire"

    ```console
    work/e0/ecf81b4cacc648b9b994218d5b29d7/
    ├── Holà-output.txt -> /workspaces/training/hello-nextflow/work/ab/81632178cd37e9e815959278808819/Holà-output.txt
    └── UPPER-Holà-output.txt
    ```

Remarquez qu'il y a deux fichiers `*-output` : la sortie du premier processus ainsi que la sortie du deuxième.

La sortie du premier processus est là parce que Nextflow l'a **stagée** là afin d'avoir tout ce qui est nécessaire pour l'exécution dans le même sous-répertoire.

Cependant, c'est en fait un lien symbolique pointant vers le fichier original dans le sous-répertoire de l'appel du premier processus.
Par défaut, lors de l'exécution sur une seule machine comme nous le faisons ici, Nextflow utilise des liens symboliques plutôt que des copies pour stager les fichiers d'entrée et intermédiaires.

Maintenant, avant de continuer, réfléchissez à comment tout ce que nous avons fait est de connecter la sortie de `sayHello` à l'entrée de `convertToUpper` et les deux processus ont pu être exécutés en série.
Nextflow a fait le travail difficile de gérer les fichiers d'entrée et de sortie individuels et de les passer entre les deux commandes pour nous.

C'est l'une des raisons pour lesquelles les canaux Nextflow sont si puissants : ils s'occupent du travail fastidieux impliqué dans la connexion des étapes du workflow ensemble.

### À retenir

Vous savez comment enchaîner des processus ensemble en fournissant la sortie d'une étape comme entrée à l'étape suivante.

### Et ensuite ?

Apprenez à collecter les sorties des appels de processus par lots et à les alimenter dans un seul processus.

---

## 2. Ajouter une troisième étape pour collecter tous les messages de bienvenue

Lorsque nous utilisons un processus pour appliquer une transformation à chacun des éléments d'un canal, comme nous le faisons ici avec les multiples messages de bienvenue, nous voulons parfois collecter les éléments du canal de sortie de ce processus et les alimenter dans un autre processus qui effectue une sorte d'analyse ou de sommation.

Pour illustrer, nous allons ajouter une nouvelle étape à notre pipeline qui collecte tous les messages de bienvenue en majuscules produits par le processus `convertToUpper` et les écrit dans un seul fichier.

<figure class="excalidraw">
--8<-- "docs/en/docs/hello_nextflow/img/hello-collect.svg"
</figure>

Sans vouloir gâcher la surprise, cela va impliquer un opérateur très utile.

### 2.1. Définir la commande de collecte et la tester dans le terminal

L'étape de collecte que nous voulons ajouter à notre workflow utilisera la commande `cat` pour concaténer plusieurs messages de bienvenue en majuscules dans un seul fichier.

Exécutons la commande seule dans le terminal pour vérifier qu'elle fonctionne comme prévu, comme nous l'avons fait précédemment.

Exécutez ce qui suit dans votre terminal :

```bash
echo 'Hello' | tr '[a-z]' '[A-Z]' > UPPER-Hello-output.txt
echo 'Bonjour' | tr '[a-z]' '[A-Z]' > UPPER-Bonjour-output.txt
echo 'Holà' | tr '[a-z]' '[A-Z]' > UPPER-Holà-output.txt
cat UPPER-Hello-output.txt UPPER-Bonjour-output.txt UPPER-Holà-output.txt > COLLECTED-output.txt
```

La sortie est un fichier texte appelé `COLLECTED-output.txt` qui contient les versions en majuscules des messages de bienvenue originaux.

??? abstract "Contenu du fichier"

    ```console title="COLLECTED-output.txt"
    HELLO
    BONJOUR
    HOLà
    ```

C'est le résultat que nous voulons obtenir avec notre workflow.

### 2.2. Créer un nouveau processus pour effectuer l'étape de collecte

Créons un nouveau processus et appelons-le `collectGreetings()`.
Nous pouvons commencer à l'écrire en nous basant sur ce que nous avons vu auparavant.

#### 2.2.1. Écrire les parties « évidentes » du processus

Ajoutez la définition de processus suivante au script de workflow :

```groovy title="hello-workflow.nf" linenums="37"
/*
 * Collecter les salutations en majuscules dans un seul fichier de sortie
 */
process collectGreetings {

    input:
    ???

    output:
    path "COLLECTED-output.txt"

    script:
    """
    cat ??? > 'COLLECTED-output.txt'
    """
}
```

C'est ce que nous pouvons écrire avec confiance en nous basant sur ce que vous avez appris jusqu'à présent.
Mais ce n'est pas fonctionnel !
Cela omet la définition des entrées et la première moitié de la commande de script parce que nous devons déterminer comment les écrire.

#### 2.2.2. Définir les entrées de `collectGreetings()`

Nous devons collecter les messages de bienvenue de tous les appels au processus `convertToUpper()`.
Que savons-nous pouvoir obtenir de l'étape précédente du workflow ?

Le canal sorti par `convertToUpper()` contiendra les chemins vers les fichiers individuels contenant les messages de bienvenue en majuscules.
Cela correspond à un emplacement d'entrée ; appelons-le `input_files` pour simplifier.

Dans le bloc du processus, effectuez la modification de code suivante :

=== "Après"

    ```groovy title="hello-workflow.nf" linenums="42" hl_lines="2"
          input:
          path input_files
    ```

=== "Avant"

    ```groovy title="hello-workflow.nf" linenums="42" hl_lines="2"
          input:
          ???
    ```

Notez que nous utilisons le préfixe `path` même si nous nous attendons à ce que cela contienne plusieurs fichiers.

#### 2.2.3. Composer la commande de concaténation

C'est là que les choses pourraient devenir un peu délicates, car nous devons être capables de gérer un nombre arbitraire de fichiers d'entrée.
Plus précisément, nous ne pouvons pas écrire la commande à l'avance, nous devons donc dire à Nextflow comment la composer à l'exécution en fonction des entrées qui circulent dans le processus.

En d'autres termes, si nous avons un canal d'entrée contenant l'élément `[file1.txt, file2.txt, file3.txt]`, nous avons besoin que Nextflow transforme cela en `cat file1.txt file2.txt file3.txt`.

Heureusement, Nextflow est tout à fait capable de faire cela pour nous si nous écrivons simplement `cat ${input_files}` dans la commande de script.

Dans le bloc du processus, effectuez la modification de code suivante :

=== "Après"

    ```groovy title="hello-workflow.nf" linenums="54" hl_lines="3"
        script:
        """
        cat ${input_files} > 'COLLECTED-output.txt'
        """
    ```

=== "Avant"

    ```groovy title="hello-workflow.nf" linenums="54"
        script:
        """
        cat ??? > 'COLLECTED-output.txt'
        """
    ```

En théorie, cela devrait gérer n'importe quel nombre arbitraire de fichiers d'entrée.

!!! tip "Astuce"

    Certains outils en ligne de commande nécessitent de fournir un argument (comme `-input`) pour chaque fichier d'entrée.
    Dans ce cas, nous devrions faire un peu de travail supplémentaire pour composer la commande.
    Vous pouvez voir un exemple de ceci dans le cours de formation [Nextflow for Genomics](../../nf4_science/genomics/).

### 2.3. Ajouter l'étape de collecte au workflow

Maintenant, nous devrions juste avoir besoin d'appeler le processus de collecte sur la sortie de l'étape de mise en majuscules.

#### 2.3.1. Connecter les appels de processus

Dans le bloc workflow, effectuez la modification de code suivante :

=== "Après"

    ```groovy title="hello-workflow.nf" linenums="75" hl_lines="4 5"
        // convertir la salutation en majuscules
        convertToUpper(sayHello.out)

        // collecter toutes les salutations dans un seul fichier
        collectGreetings(convertToUpper.out)
    }
    ```

=== "Avant"

    ```groovy title="hello-workflow.nf" linenums="75"
        // convertir la salutation en majuscules
        convertToUpper(sayHello.out)
    }
    ```

Cela connecte la sortie de `convertToUpper()` à l'entrée de `collectGreetings()`.

#### 2.3.2. Exécuter le workflow avec `-resume`

Essayons.

```bash
nextflow run hello-workflow.nf -resume
```

??? success "Sortie de la commande"

    ```console hl_lines="8"
    N E X T F L O W   ~  version 25.10.2

    Launching `hello-workflow.nf` [mad_gilbert] DSL2 - revision: 6acfd5e28d

    executor >  local (3)
    [79/33b2f0] sayHello (2)         | 3 of 3, cached: 3 ✔
    [99/79394f] convertToUpper (3)   | 3 of 3, cached: 3 ✔
    [47/50fe4a] collectGreetings (1) | 3 of 3 ✔
    ```

Il s'exécute avec succès, y compris la troisième étape.

Cependant, regardez le nombre d'appels pour `collectGreetings()` sur la dernière ligne.
Nous n'en attendions qu'un, mais il y en a trois.

Maintenant, regardez le contenu du fichier de sortie final.

??? abstract "Contenu du fichier"

    ```console title="results/COLLECTED-output.txt"
    Holà
    ```

Oh non. L'étape de collecte a été exécutée individuellement sur chaque message de bienvenue, ce qui N'EST PAS ce que nous voulions.

Nous devons faire quelque chose pour dire explicitement à Nextflow que nous voulons que cette troisième étape s'exécute sur tous les éléments du canal sorti par `convertToUpper()`.

### 2.4. Utiliser un opérateur pour collecter les messages de bienvenue dans une seule entrée

Oui, encore une fois la réponse à notre problème est un opérateur.

Plus précisément, nous allons utiliser l'opérateur bien nommé [`collect()`](https://www.nextflow.io/docs/latest/reference/operator.html#collect).

#### 2.4.1. Ajouter l'opérateur `collect()`

Cette fois, cela va sembler un peu différent parce que nous n'ajoutons pas l'opérateur dans le contexte d'une channel factory ; nous l'ajoutons à un canal de sortie.

Nous prenons le `convertToUpper.out` et ajoutons l'opérateur `collect()`, ce qui nous donne `convertToUpper.out.collect()`.
Nous pouvons brancher cela directement dans l'appel du processus `collectGreetings()`.

Dans le bloc workflow, effectuez la modification de code suivante :

=== "Après"

    ```groovy title="hello-workflow.nf" linenums="73" hl_lines="2"
        // collecter toutes les salutations dans un seul fichier
        collectGreetings(convertToUpper.out.collect())
    }
    ```

=== "Avant"

    ```groovy title="hello-workflow.nf" linenums="73" hl_lines="2"
        // collecter toutes les salutations dans un seul fichier
        collectGreetings(convertToUpper.out)
    }
    ```

#### 2.4.2. Ajouter quelques instructions `view()`

Incluons également quelques instructions `view()` pour visualiser les états avant et après du contenu du canal.

=== "Après"

    ```groovy title="hello-workflow.nf" linenums="73" hl_lines="4-6"
        // collecter toutes les salutations dans un seul fichier
        collectGreetings(convertToUpper.out.collect())

        // optional view statements
        convertToUpper.out.view { contents -> "Before collect: $contents" }
        convertToUpper.out.collect().view { contents -> "After collect: $contents" }
    }
    ```

=== "Avant"

    ```groovy title="hello-workflow.nf" linenums="73"
        // collecter toutes les salutations dans un seul fichier
        collectGreetings(convertToUpper.out.collect())
    }
    ```

Les instructions `view()` peuvent aller où vous voulez ; nous les avons mises juste après l'appel pour la lisibilité.

#### 2.4.3. Exécuter à nouveau le workflow avec `-resume`

Essayons :

```bash
nextflow run hello-workflow.nf -resume
```

??? success "Sortie de la commande"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `hello-workflow.nf` [soggy_franklin] DSL2 - revision: bc8e1b2726

    [d6/cdf466] sayHello (1)       | 3 of 3, cached: 3 ✔
    [99/79394f] convertToUpper (2) | 3 of 3, cached: 3 ✔
    [1e/83586c] collectGreetings   | 1 of 1 ✔
    Before collect: /workspaces/training/hello-nextflow/work/b3/d52708edba8b864024589285cb3445/UPPER-Bonjour-output.txt
    Before collect: /workspaces/training/hello-nextflow/work/99/79394f549e3040dfc2440f69ede1fc/UPPER-Hello-output.txt
    Before collect: /workspaces/training/hello-nextflow/work/aa/56bfe7cf00239dc5badc1d04b60ac4/UPPER-Holà-output.txt
    After collect: [/workspaces/training/hello-nextflow/work/b3/d52708edba8b864024589285cb3445/UPPER-Bonjour-output.txt, /workspaces/training/hello-nextflow/work/99/79394f549e3040dfc2440f69ede1fc/UPPER-Hello-output.txt, /workspaces/training/hello-nextflow/work/aa/56bfe7cf00239dc5badc1d04b60ac4/UPPER-Holà-output.txt]
    ```

Il s'exécute avec succès, bien que la sortie du journal puisse sembler un peu plus désordonnée que ceci (nous l'avons nettoyée pour la lisibilité).

Cette fois, la troisième étape n'a été appelée qu'une seule fois !

En regardant la sortie des instructions `view()`, nous voyons ce qui suit :

- Trois instructions `Before collect:`, une pour chaque message de bienvenue : à ce stade, les chemins de fichiers sont des éléments individuels dans le canal.
- Une seule instruction `After collect:` : les trois chemins de fichiers sont maintenant empaquetés dans un seul élément.

Regardez le contenu du fichier de sortie final.

??? abstract "Contenu du fichier"

    ```console title="results/COLLECTED-output.txt"
    BONJOUR
    HELLO
    HOLà
    ```

Cette fois, nous avons les trois messages de bienvenue dans le fichier de sortie final. Succès !

!!! note "Note"

    Si vous exécutez ceci plusieurs fois sans `-resume`, vous verrez que l'ordre des messages de bienvenue change d'une exécution à l'autre.
    Cela vous montre que l'ordre dans lequel les éléments circulent à travers les appels de processus n'est pas garanti d'être cohérent.

#### 2.4.4. Supprimer les instructions `view()` pour la lisibilité

Avant de passer à la section suivante, nous vous recommandons de supprimer les instructions `view()` pour éviter d'encombrer la sortie de la console.

=== "Après"

    ```groovy title="hello-workflow.nf" linenums="73"
        // collecter toutes les salutations dans un seul fichier
        collectGreetings(convertToUpper.out.collect())
    ```

=== "Avant"

    ```groovy title="hello-workflow.nf" linenums="73" hl_lines="4-6"
        // collecter toutes les salutations dans un seul fichier
        collectGreetings(convertToUpper.out.collect())

        // optional view statements
        convertToUpper.out.view { contents -> "Before collect: $contents" }
        convertToUpper.out.collect().view { contents -> "After collect: $contents" }
    ```

C'est essentiellement l'opération inverse du point 2.4.2.

### À retenir

Vous savez comment collecter les sorties d'un lot d'appels de processus et les alimenter dans une étape d'analyse ou de sommation conjointe.

### Et ensuite ?

Apprenez à passer plus d'une entrée à un processus.

---

## 3. Passer plus d'une entrée à un processus

Nous voulons pouvoir nommer le fichier de sortie final avec un nom spécifique afin de traiter des lots ultérieurs de messages de bienvenue sans écraser les résultats finaux.

Pour ce faire, nous allons apporter les améliorations suivantes au workflow :

- Modifier le processus collecteur pour accepter un nom défini par l'utilisateur pour le fichier de sortie
- Ajouter un paramètre de ligne de commande au workflow et le passer au processus collecteur

### 3.1. Modifier le processus collecteur

Nous allons devoir déclarer l'entrée supplémentaire et l'intégrer dans le nom du fichier de sortie.

#### 3.1.1. Déclarer l'entrée supplémentaire

Bonne nouvelle : nous pouvons déclarer autant de variables d'entrée que nous voulons dans la définition du processus.
Appelons celle-ci `batch_name`.

Dans le bloc du processus, effectuez la modification de code suivante :

=== "Après"

    ```groovy title="hello-workflow.nf" linenums="42" hl_lines="3"
        input:
        path input_files
        val batch_name
    ```

=== "Avant"

    ```groovy title="hello-workflow.nf" linenums="42"
        input:
        path input_files
    ```

Vous pouvez configurer vos processus pour attendre autant d'entrées que vous le souhaitez.
Pour le moment, ce sont toutes des entrées requises ; vous _devez_ fournir une valeur pour que le workflow fonctionne.

Vous apprendrez comment gérer les entrées requises vs. optionnelles plus tard dans votre parcours Nextflow.

#### 3.1.2. Utiliser la variable `batch_name` dans le nom du fichier de sortie

Nous pouvons insérer la variable dans le nom du fichier de sortie de la même manière que nous avons composé des noms de fichiers dynamiques auparavant.

Dans le bloc du processus, effectuez la modification de code suivante :

=== "Après"

    ```groovy title="hello-workflow.nf" linenums="46" hl_lines="2 6"
        output:
        path "COLLECTED-${batch_name}-output.txt"

        script:
        """
        cat ${input_files} > 'COLLECTED-${batch_name}-output.txt'
        """
    ```

=== "Avant"

    ```groovy title="hello-workflow.nf" linenums="46" hl_lines="2 6"
        output:
        path "COLLECTED-output.txt"

        script:
        """
        cat ${input_files} > 'COLLECTED-output.txt'
        """
    ```

Cela configure le processus pour utiliser la valeur `batch_name` afin de générer un nom de fichier spécifique pour la sortie finale du workflow.

### 3.2. Ajouter un paramètre de ligne de commande `batch`

Maintenant, nous avons besoin d'un moyen de fournir la valeur pour `batch_name` et de l'alimenter à l'appel du processus.

#### 3.2.1. Utiliser `params` pour configurer le paramètre

Vous savez déjà comment utiliser le système `params` pour déclarer des paramètres CLI.
Utilisons cela pour déclarer un paramètre `batch` (avec une valeur par défaut parce que nous sommes paresseux).

Dans la section des paramètres du pipeline, effectuez les modifications de code suivantes :

=== "Après"

    ```groovy title="hello-workflow.nf" linenums="55" hl_lines="6"
    /*
     * Paramètres du pipeline
     */
    params {
        input: Path = 'data/greetings.csv'
        batch: String = 'batch'
    }
    ```

=== "Avant"

    ```groovy title="hello-workflow.nf" linenums="55"
    /*
     * Paramètres du pipeline
     */
    params {
        input: Path = 'data/greetings.csv'
    }
    ```

Tout comme nous l'avons démontré pour `--input`, vous pouvez remplacer cette valeur par défaut en spécifiant une valeur avec `--batch` sur la ligne de commande.

#### 3.2.2. Passer le paramètre `batch` au processus

Pour fournir la valeur du paramètre au processus, nous devons l'ajouter dans l'appel du processus.

Dans le bloc workflow, effectuez la modification de code suivante :

=== "Après"

    ```groovy title="hello-workflow.nf" linenums="74" hl_lines="2"
        // collecter toutes les salutations dans un seul fichier
        collectGreetings(convertToUpper.out.collect(), params.batch)
    ```

=== "Avant"

    ```groovy title="hello-workflow.nf" linenums="74" hl_lines="2"
        // collecter toutes les salutations dans un seul fichier
        collectGreetings(convertToUpper.out.collect())
    ```

Vous voyez que pour fournir plusieurs entrées à un processus, vous les listez simplement dans les parenthèses de l'appel, séparées par des virgules.

!!! warning "Avertissement"

    Vous DEVEZ fournir les entrées au processus dans le MÊME ORDRE EXACT qu'elles sont listées dans le bloc de définition des entrées du processus.

### 3.3. Exécuter le workflow

Essayons d'exécuter ceci avec un nom de lot sur la ligne de commande.

```bash
nextflow run hello-workflow.nf -resume --batch trio
```

??? success "Sortie de la commande"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `hello-workflow.nf` [confident_rutherford] DSL2 - revision: bc58af409c

    executor >  local (1)
    [79/33b2f0] sayHello (2)       | 3 of 3, cached: 3 ✔
    [99/79394f] convertToUpper (2) | 3 of 3, cached: 3 ✔
    [b5/f19efe] collectGreetings   | 1 of 1 ✔
    ```

Il s'exécute avec succès et produit la sortie désirée :

??? abstract "Contenu du fichier"

    ```console title="results/COLLECTED-trio-output.txt"
    HELLO
    BONJOUR
    HOLà
    ```

Maintenant, tant que nous spécifions le paramètre de manière appropriée, les exécutions ultérieures sur d'autres lots d'entrées n'écraseront pas les résultats précédents.

### À retenir

Vous savez comment passer plus d'une entrée à un processus.

### Et ensuite ?

Apprenez à émettre plusieurs sorties et à les gérer de manière pratique.

---

## 4. Ajouter une sortie à l'étape collecteur

Jusqu'à présent, nous avons utilisé des processus qui ne produisaient qu'une seule sortie chacun.
Nous avons pu accéder à leurs sorties respectives de manière très pratique en utilisant la syntaxe `<process>.out`, que nous avons utilisée à la fois dans le contexte du passage d'une sortie au processus suivant (par exemple `convertToUpper(sayHello.out)`) et dans le contexte de la section `publish:` (par exemple `first_output = sayHello.out`).

Que se passe-t-il quand un processus en produit plus d'une ?
Comment gérons-nous les sorties multiples ?
Pouvons-nous sélectionner et utiliser une sortie spécifique ?

Toutes d'excellentes questions, et la réponse courte est oui, nous le pouvons !

Les sorties multiples seront empaquetées dans des canaux séparés.
Nous pouvons soit choisir de donner des noms à ces canaux de sortie, ce qui facilite la référence individuelle plus tard, soit nous y référer par index.

Explorons avec un exemple.

À des fins de démonstration, disons que nous voulons compter le nombre de messages de bienvenue qui sont collectés pour un lot donné d'entrées et le rapporter dans un fichier.

### 4.1. Modifier le processus pour compter et sortir le nombre de messages de bienvenue

Cela nécessitera deux changements clés à la définition du processus : nous avons besoin d'un moyen de compter les messages de bienvenue et d'écrire un fichier de rapport, puis nous devons ajouter ce fichier de rapport au bloc `output` du processus.

#### 4.1.1. Compter le nombre de messages de bienvenue collectés

De manière pratique, Nextflow nous permet d'ajouter du code arbitraire dans le bloc `script:` de la définition du processus, ce qui s'avère vraiment pratique pour faire des choses comme ceci.

Cela signifie que nous pouvons utiliser la fonction intégrée `size()` de Nextflow pour obtenir le nombre de fichiers dans le tableau `input_files`, et écrire le résultat dans un fichier avec une commande `echo`.

Dans le bloc du processus `collectGreetings`, effectuez les modifications de code suivantes :

=== "Après"

    ```groovy title="hello-workflow.nf" linenums="55" hl_lines="2 5"
        script:
        count_greetings = input_files.size()
        """
        cat ${input_files} > 'COLLECTED-${batch_name}-output.txt'
        echo 'There were ${count_greetings} greetings in this batch.' > '${batch_name}-report.txt'
        """
    ```

=== "Avant"

    ```groovy title="hello-workflow.nf" linenums="55"
        script:
        """
        cat ${input_files} > 'COLLECTED-${batch_name}-output.txt'
        """
    ```

La variable `count_greetings` sera calculée à l'exécution.

#### 4.1.2. Émettre le fichier de rapport et nommer les sorties

En principe, tout ce que nous devons faire est d'ajouter le fichier de rapport au bloc `output:`.

Cependant, pendant que nous y sommes, nous allons également ajouter quelques balises `emit:` à nos déclarations de sortie. Celles-ci nous permettront de sélectionner les sorties par nom au lieu d'avoir à utiliser des indices positionnels.

Dans le bloc du processus, effectuez la modification de code suivante :

=== "Après"

    ```groovy title="hello-workflow.nf" linenums="46" hl_lines="2 3"
        output:
        path "COLLECTED-${batch_name}-output.txt", emit: outfile
        path "${batch_name}-report.txt", emit: report
    ```

=== "Avant"

    ```groovy title="hello-workflow.nf" linenums="46"
        output:
        path "COLLECTED-${batch_name}-output.txt"
    ```

Les balises `emit:` sont optionnelles, et nous aurions pu ajouter une balise à une seule des sorties.
Mais comme dit le dicton, pourquoi pas les deux ?

!!! tip "Astuce"

    Si vous ne nommez pas les sorties d'un processus en utilisant `emit:`, vous pouvez toujours y accéder individuellement en utilisant leur index respectif (basé sur zéro).
    Par exemple, vous utiliseriez `<process>.out[0]` pour obtenir la première sortie, `<process>.out[1]` pour obtenir la deuxième sortie, et ainsi de suite.

    Nous préférons nommer les sorties parce que sinon, il est trop facile de saisir le mauvais index par erreur, surtout quand le processus produit beaucoup de sorties.

### 4.2. Mettre à jour les sorties du workflow

Maintenant que nous avons deux sorties provenant du processus `collectGreetings`, la sortie `collectGreetings.out` contient deux canaux :

- `collectGreetings.out.outfile` contient le fichier de sortie final
- `collectGreetings.out.report` contient le fichier de rapport

Nous devons mettre à jour les sorties du workflow en conséquence.

#### 4.2.1. Mettre à jour la section `publish:`

Dans le bloc `workflow`, effectuez la modification de code suivante :

=== "Après"

    ```groovy title="hello-workflow.nf" linenums="80" hl_lines="4 5"
        publish:
        first_output = sayHello.out
        uppercased = convertToUpper.out
        collected = collectGreetings.out.outfile
        batch_report = collectGreetings.out.report
    ```

=== "Avant"

    ```groovy title="hello-workflow.nf" linenums="80" hl_lines="4"
        publish:
        first_output = sayHello.out
        uppercased = convertToUpper.out
        collected = collectGreetings.out
    ```

Comme vous pouvez le voir, se référer à des sorties de processus spécifiques est maintenant trivial.
Quand nous ajouterons une étape de plus à notre pipeline dans la Partie 5 (Conteneurs), nous pourrons facilement nous référer à `collectGreetings.out.outfile` et le passer au nouveau processus (spoiler : le nouveau processus s'appelle `cowpy`).

Mais pour l'instant, finissons de mettre à jour les sorties au niveau du workflow.

#### 4.2.2. Mettre à jour le bloc `output`

Dans le bloc `output`, effectuez la modification de code suivante :

=== "Après"

    ```groovy title="hello-workflow.nf" linenums="86" hl_lines="14-17"
    output {
        first_output {
            path 'hello_workflow'
            mode 'copy'
        }
        uppercased {
            path 'hello_workflow'
            mode 'copy'
        }
        collected {
            path 'hello_workflow'
            mode 'copy'
        }
        batch_report {
            path 'hello_workflow'
            mode 'copy'
        }
    }
    ```

=== "Avant"

    ```groovy title="hello-workflow.nf" linenums="80"
    output {
        first_output {
            path 'hello_workflow'
            mode 'copy'
        }
        uppercased {
            path 'hello_workflow'
            mode 'copy'
        }
        collected {
            path 'hello_workflow'
            mode 'copy'
        }
    }
    ```

Nous n'avons pas besoin de mettre à jour la définition de sortie `collected` puisque ce nom n'a pas changé.
Nous devons juste ajouter la nouvelle sortie.

### 4.3. Exécuter le workflow

Essayons d'exécuter ceci avec le lot actuel de messages de bienvenue.

```bash
nextflow run hello-workflow.nf -resume --batch trio
```

??? success "Sortie de la commande"

    ```console
     N E X T F L O W   ~  version 25.10.2

    Launching `hello-workflow.nf` [ecstatic_wilson] DSL2 - revision: c80285f8c8

    executor >  local (1)
    [c5/4c6ca9] sayHello (3)       [100%] 3 of 3, cached: 3 ✔
    [0e/6cbc59] convertToUpper (3) [100%] 3 of 3, cached: 3 ✔
    [02/61ead2] collectGreetings   [100%] 1 of 1 ✔
    ```

Si vous regardez dans le répertoire `results/hello_workflow/`, vous trouverez le nouveau fichier de rapport, `trio-report.txt`.
Ouvrez-le pour vérifier que le workflow a correctement rapporté le nombre de messages de bienvenue qui ont été traités.

??? abstract "Contenu du fichier"

    ```txt title="trio-report.txt"
    There were 3 greetings in this batch.
    ```

N'hésitez pas à ajouter plus de messages de bienvenue au CSV et tester ce qui se passe.

### À retenir

Vous savez comment faire en sorte qu'un processus émette plusieurs sorties nommées et comment les gérer de manière appropriée au niveau du workflow.

Plus généralement, vous comprenez les principes clés impliqués dans la connexion des processus ensemble de manières courantes.

### Et ensuite ?

Prenez une pause extra longue, vous l'avez bien méritée.

Quand vous êtes prêt·e, passez à la [**Partie 4 : Hello Modules**](./04_hello_modules.md) pour apprendre comment modulariser votre code pour une meilleure maintenabilité et efficacité du code.

---

## Quiz

<quiz>
Comment accédez-vous à la sortie d'un processus dans le bloc workflow ?
- [ ] `process.output`
- [ ] `output.processName`
- [x] `processName.out`
- [ ] `get(processName)`

En savoir plus : [1.4. Passer la sortie du premier processus au deuxième processus](#14-passer-la-sortie-du-premier-processus-au-deuxieme-processus)
</quiz>

<quiz>
Qu'est-ce qui détermine l'ordre d'exécution des processus dans Nextflow ?
- [ ] L'ordre dans lequel les processus sont écrits dans le bloc workflow
- [ ] L'ordre alphabétique par nom de processus
- [x] Les dépendances de données entre les processus
- [ ] Ordre aléatoire pour l'exécution parallèle

En savoir plus : [1.4. Passer la sortie du premier processus au deuxième processus](#14-passer-la-sortie-du-premier-processus-au-deuxieme-processus)
</quiz>

<quiz>
Quel opérateur doit remplacer `???` pour rassembler toutes les sorties dans une seule liste pour le processus en aval ?

```groovy hl_lines="4"
workflow {
    greetings_ch = Channel.of('Hello', 'Bonjour', 'Hola')
    SAYHELLO(greetings_ch)
    GATHER_ALL(SAYHELLO.out.???)
}
```

- [ ] `flatten()`
- [x] `collect()`
- [ ] `mix()`
- [ ] `join()`

En savoir plus : [2.4. Utiliser un opérateur pour collecter les messages de bienvenue dans une seule entrée](#24-utiliser-un-operateur-pour-collecter-les-messages-de-bienvenue-dans-une-seule-entree)
</quiz>

<quiz>
Quand devriez-vous utiliser l'opérateur `collect()` ?
- [ ] Quand vous voulez traiter les éléments en parallèle
- [ ] Quand vous devez filtrer le contenu d'un canal
- [x] Quand un processus en aval a besoin de tous les éléments d'un processus en amont
- [ ] Quand vous voulez diviser les données entre plusieurs processus

En savoir plus : [2.4. Utiliser un opérateur pour collecter les messages de bienvenue dans une seule entrée](#24-utiliser-un-operateur-pour-collecter-les-messages-de-bienvenue-dans-une-seule-entree)
</quiz>

<quiz>
Comment accédez-vous à une sortie nommée d'un processus ?
- [ ] `processName.outputName`
- [ ] `processName.get(outputName)`
- [x] `processName.out.outputName`
- [ ] `output.processName.outputName`

En savoir plus : [4.1.2. Émettre le fichier de rapport et nommer les sorties](#412-emettre-le-fichier-de-rapport-et-nommer-les-sorties)
</quiz>

<quiz>
Quelle est la syntaxe correcte pour nommer une sortie dans un processus ?
- [ ] `name: outputName`
- [ ] `output: outputName`
- [x] `emit: outputName`
- [ ] `label: outputName`

En savoir plus : [4.1.2. Émettre le fichier de rapport et nommer les sorties](#412-emettre-le-fichier-de-rapport-et-nommer-les-sorties)
</quiz>

<quiz>
Lors de la fourniture de plusieurs entrées à un processus, qu'est-ce qui doit être vrai ?
- [ ] Toutes les entrées doivent être du même type
- [ ] Les entrées doivent être fournies dans l'ordre alphabétique
- [x] L'ordre des entrées doit correspondre à l'ordre défini dans le bloc d'entrée
- [ ] Seules deux entrées peuvent être fournies à la fois

En savoir plus : [3. Passer plus d'une entrée à un processus](#3-passer-plus-dune-entree-a-un-processus)
</quiz>
