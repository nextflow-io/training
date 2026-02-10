# Environnement de Développement

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Traduction assistée par IA - [en savoir plus et suggérer des améliorations](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

Les environnements de développement intégrés (IDE) modernes peuvent transformer radicalement votre expérience de développement Nextflow. Cette quête secondaire se concentre spécifiquement sur l'exploitation de VS Code et de son extension Nextflow pour écrire du code plus rapidement, détecter les erreurs tôt et naviguer efficacement dans des workflows complexes.

!!! note "Ce n'est pas un tutoriel traditionnel"

    Contrairement aux autres modules de formation, ce guide est organisé comme une collection de conseils rapides, d'astuces et d'exemples pratiques plutôt qu'un tutoriel étape par étape. Chaque section peut être explorée indépendamment en fonction de vos intérêts et de vos besoins actuels de développement. N'hésitez pas à naviguer librement et à vous concentrer sur les fonctionnalités qui seront les plus immédiatement utiles pour le développement de votre workflow.

## Ce que vous devriez savoir d'abord

Ce guide suppose que vous avez terminé le cours de formation [Hello Nextflow](../hello_nextflow/) et que vous êtes à l'aise avec les concepts fondamentaux de Nextflow, notamment :

- **Structure de workflow de base** : Comprendre les processus, les workflows et comment ils se connectent ensemble
- **Opérations sur les canaux** : Créer des canaux, transmettre des données entre processus et utiliser les opérateurs de base
- **Modules et organisation** : Créer des modules réutilisables et utiliser les instructions include
- **Bases de la configuration** : Utiliser `nextflow.config` pour les paramètres, les directives de processus et les profils

## Ce que vous apprendrez ici

Ce guide se concentre sur les **fonctionnalités de productivité de l'IDE** qui feront de vous un·e développeur·se Nextflow plus efficace :

- **Coloration syntaxique avancée** : Comprendre ce que VS Code vous montre sur la structure de votre code
- **Auto-complétion intelligente** : Exploiter les suggestions contextuelles pour écrire du code plus rapidement
- **Détection d'erreurs et diagnostics** : Détecter les erreurs de syntaxe avant d'exécuter votre workflow
- **Navigation dans le code** : Se déplacer rapidement entre processus, modules et définitions
- **Formatage et organisation** : Maintenir un style de code cohérent et lisible
- **Développement assisté par IA** (optionnel) : Utiliser des outils IA modernes intégrés à votre IDE

!!! info "Pourquoi les fonctionnalités IDE maintenant ?"

    Vous avez probablement déjà utilisé VS Code pendant le cours [Hello Nextflow](../hello_nextflow/), mais nous avons maintenu l'accent sur l'apprentissage des fondamentaux de Nextflow plutôt que sur les fonctionnalités de l'IDE. Maintenant que vous êtes à l'aise avec les concepts de base de Nextflow comme les processus, workflows, canaux et modules, vous êtes prêt·e à exploiter les fonctionnalités sophistiquées de l'IDE qui feront de vous un développeur·se plus efficace.

    Considérez cela comme une « montée en compétence » de votre environnement de développement - l'éditeur que vous utilisez déjà possède des capacités bien plus puissantes qui deviennent vraiment précieuses une fois que vous comprenez ce qu'elles vous apportent.

---

## 0. Configuration et Préparation

Configurons un espace de travail spécifiquement pour explorer les fonctionnalités de l'IDE :

```bash title="Naviguer vers le répertoire des fonctionnalités IDE"
cd side-quests/ide_features
```

Ouvrez ce répertoire dans VS Code :

```bash title="Ouvrir VS Code dans le répertoire courant"
code .
```

Le répertoire `ide_features` contient des exemples de workflows qui démontrent diverses fonctionnalités de l'IDE :

```bash title="Afficher la structure du répertoire"
tree .
```

```console title="Structure du projet"
tree .
.
├── basic_workflow.nf
├── complex_workflow.nf
├── data
│   ├── sample_001.fastq.gz
│   ├── sample_002.fastq.gz
│   ├── sample_003.fastq.gz
│   ├── sample_004.fastq.gz
│   ├── sample_005.fastq.gz
│   └── sample_data.csv
├── modules
│   ├── fastqc.nf
│   ├── star.nf
│   └── utils.nf
└── nextflow.config

3 directories, 12 files
```

!!! note "À propos des fichiers d'exemple"

    - `basic_workflow.nf` est un workflow de base fonctionnel que vous pouvez exécuter et modifier
    - `complex_workflow.nf` est conçu uniquement à titre d'illustration pour démontrer les fonctionnalités de navigation - il peut ne pas s'exécuter correctement mais montre une structure de workflow multi-fichiers réaliste

### Raccourcis Clavier

Certaines des fonctionnalités de ce guide utiliseront des raccourcis clavier optionnels. Vous accédez peut-être à ce matériel via GitHub Codespaces dans un navigateur, et dans ce cas, les raccourcis ne fonctionneront parfois pas comme prévu car ils sont utilisés pour d'autres fonctions dans votre système.

Si vous exécutez VS Code localement, comme vous le ferez probablement lorsque vous écrirez réellement des workflows, les raccourcis fonctionneront comme décrit.

Si vous utilisez un Mac, certains (pas tous) raccourcis clavier utiliseront « cmd » au lieu de « ctrl », et nous l'indiquerons dans le texte comme `Ctrl/Cmd`.

### 0.1. Installation de l'Extension Nextflow

!!! note "Vous utilisez déjà des Devcontainers ?"

    Si vous travaillez dans **GitHub Codespaces** ou utilisez un **devcontainer local**, l'extension Nextflow est probablement déjà installée et configurée pour vous. Vous pouvez ignorer les étapes d'installation manuelle ci-dessous et passer directement à l'exploration des fonctionnalités de l'extension.

Pour installer l'extension manuellement :

1. Ouvrez VS Code
2. Accédez à la vue Extensions en cliquant sur l'icône des extensions à gauche : ![icône des extensions](img/extensions_icon.png) (raccourci `Ctrl/Cmd+Shift+X` si vous exécutez VSCode localement)
3. Recherchez « Nextflow »
4. Installez l'extension officielle Nextflow

![Installer l'extension Nextflow](img/install_extension.png)

### 0.2. Disposition de l'Espace de Travail

Puisque vous avez utilisé VS Code tout au long de Hello Nextflow, vous connaissez déjà les bases. Voici comment organiser efficacement votre espace de travail pour cette session :

- **Zone d'Édition** : Pour visualiser et éditer les fichiers. Vous pouvez diviser cette zone en plusieurs panneaux pour comparer des fichiers côte à côte.
- **Explorateur de Fichiers** cliquez (![icône de l'explorateur de fichiers](img/files_icon.png)) (`Ctrl/Cmd+Shift+E`) : Les fichiers et dossiers locaux sur votre système. Gardez-le ouvert à gauche pour naviguer entre les fichiers
- **Terminal Intégré** (`Ctrl+Shift+` backtick pour Windows et MacOS) : Un terminal pour interagir avec l'ordinateur en bas. Utilisez-le pour exécuter Nextflow ou d'autres commandes.
- **Panneau Problèmes** (`Ctrl+Shift+M`) : VS Code affichera ici toutes les erreurs et problèmes qu'il détecte. Ceci est utile pour mettre en évidence les problèmes d'un coup d'œil.

Vous pouvez faire glisser les panneaux ou les masquer (`Ctrl/Cmd+B` pour basculer la barre latérale) afin de personnaliser votre disposition pendant que nous travaillons sur les exemples.

### Points Clés

Vous avez VS Code configuré avec l'extension Nextflow et comprenez la disposition de l'espace de travail pour un développement efficace.

### Et Ensuite ?

Découvrez comment la coloration syntaxique vous aide à comprendre la structure du code Nextflow d'un seul coup d'œil.

---

## 1. Coloration Syntaxique et Structure du Code

Maintenant que votre espace de travail est configuré, explorons comment la coloration syntaxique de VS Code vous aide à lire et écrire du code Nextflow plus efficacement.

### 1.1. Éléments de Syntaxe Nextflow

Ouvrez `basic_workflow.nf` pour voir la coloration syntaxique en action :

![Présentation de la syntaxe](img/syntax_showcase.png)

Remarquez comment VS Code met en évidence :

- **Les mots-clés** (`process`, `workflow`, `input`, `output`, `script`) dans des couleurs distinctes
- **Les chaînes littérales** et **les paramètres** avec un style différent
- **Les commentaires** dans une couleur atténuée
- **Les variables** et **les appels de fonction** avec l'emphase appropriée
- **Les blocs de code** avec des guides d'indentation appropriés

!!! note "Couleurs Dépendantes du Thème"

    Les couleurs spécifiques que vous voyez dépendront de votre thème VS Code (mode sombre/clair), des paramètres de couleur et de toute personnalisation que vous avez effectuée. L'important est que les différents éléments de syntaxe soient visuellement distingués les uns des autres, rendant la structure du code plus facile à comprendre quel que soit le schéma de couleurs choisi.

### 1.2. Comprendre la Structure du Code

La coloration syntaxique vous aide à identifier rapidement :

- **Les limites des processus** : Distinction claire entre différents processus
- **Les blocs d'entrée/sortie** : Facile à repérer les définitions de flux de données
- **Les blocs script** : Les commandes réellement exécutées
- **Les opérations sur les canaux** : Les étapes de transformation de données
- **Les directives de configuration** : Les paramètres spécifiques aux processus

Cette organisation visuelle devient inestimable lorsque vous travaillez avec des workflows complexes contenant plusieurs processus et des flux de données complexes.

### Points Clés

Vous comprenez comment la coloration syntaxique de VS Code vous aide à lire la structure du code Nextflow et à identifier différents éléments du langage pour un développement plus rapide.

### Et Ensuite ?

Découvrez comment l'auto-complétion intelligente accélère l'écriture de code avec des suggestions contextuelles.

---

## 2. Auto-complétion Intelligente

Les fonctionnalités d'auto-complétion de VS Code vous aident à écrire du code plus rapidement et avec moins d'erreurs en suggérant des options appropriées en fonction du contexte.

### 2.1. Suggestions Contextuelles

Les options d'auto-complétion varient en fonction de l'endroit où vous vous trouvez dans votre code :

#### Opérations sur les Canaux

Ouvrez à nouveau `basic_workflow.nf` et essayez de taper `channel.` dans le bloc workflow :

![Auto-complétion des canaux](img/autocomplete_channel.png)

Vous verrez des suggestions pour :

- `fromPath()` - Créer un canal à partir de chemins de fichiers
- `fromFilePairs()` - Créer un canal à partir de fichiers appariés
- `of()` - Créer un canal à partir de valeurs
- `fromSRA()` - Créer un canal à partir d'accessions SRA
- Et bien d'autres...

Cela vous aide à trouver rapidement la bonne fabrique de canaux à utiliser sans avoir besoin de mémoriser les noms exacts des méthodes.

Vous pouvez également découvrir les opérateurs disponibles à appliquer aux canaux. Par exemple, tapez `FASTQC.out.html.` pour voir les opérations disponibles :

![Auto-complétion des opérateurs de canaux](img/autocomplete_operators.png)

#### Directives de Processus

À l'intérieur d'un bloc script de processus, tapez `task.` pour voir les propriétés d'exécution disponibles :

![Auto-complétion des propriétés task](img/autocomplete_task.png)

#### Configuration

Ouvrez nextflow.config et tapez `process.` n'importe où pour voir les directives de processus disponibles :

![Auto-complétion de configuration](img/autocomplete_config.png)

Vous verrez des suggestions pour :

- `executor`
- `memory`
- `cpus`

Cela fait gagner du temps lors de la configuration des processus et fonctionne dans différents contextes de configuration. Par exemple, essayez de taper `docker.` pour voir les options de configuration spécifiques à Docker.

### Points Clés

Vous pouvez utiliser l'auto-complétion intelligente de VS Code pour découvrir les opérations de canaux disponibles, les directives de processus et les options de configuration sans mémoriser la syntaxe.

### Et Ensuite ?

Découvrez comment la détection d'erreurs en temps réel vous aide à détecter les problèmes avant d'exécuter votre workflow, simplement en lisant le code.

## 3. Détection d'Erreurs et Diagnostics

La détection d'erreurs en temps réel de VS Code vous aide à détecter les problèmes avant d'exécuter votre workflow.

### 3.1. Détection d'Erreurs de Syntaxe

Créons une erreur délibérée pour voir la détection en action. Ouvrez `basic_workflow.nf` et changez le nom du processus de `FASTQC` à `FASTQ` (ou tout autre nom invalide). VS Code mettra immédiatement en évidence l'erreur dans le bloc workflow avec un soulignement ondulé rouge :

![Soulignement d'erreur](img/error_underline.png)

### 3.2. Panneau Problèmes

Au-delà de la mise en évidence individuelle des erreurs, VS Code fournit un panneau Problèmes centralisé qui agrège toutes les erreurs, avertissements et messages d'information dans votre espace de travail. Ouvrez-le avec `Ctrl/Cmd+Shift+M` et utilisez l'icône de filtre pour afficher uniquement les erreurs pertinentes pour le fichier actuel :

![Filtrer le panneau problèmes](img/active_file.png)

Cliquez sur n'importe quel problème pour accéder directement à la ligne problématique

![Panneau Problèmes](img/problems_panel.png)

Corrigez l'erreur en changeant le nom du processus pour revenir à `FASTQC`.

### 3.3. Modèles d'Erreurs Courants

Les erreurs courantes dans la syntaxe Nextflow incluent :

- **Parenthèses manquantes** : `{` ou `}` non appariés
- **Blocs incomplets** : Sections requises manquantes dans les processus
- **Syntaxe invalide** : DSL Nextflow mal formé
- **Fautes de frappe dans les mots-clés** : Directives de processus mal orthographiées
- **Incompatibilités de canaux** : Incompatibilités de types

Le serveur de langage Nextflow met en évidence ces problèmes dans le panneau Problèmes. Vous pouvez les vérifier tôt pour éviter les erreurs de syntaxe lors de l'exécution d'un pipeline.

### Points Clés

Vous pouvez utiliser la détection d'erreurs de VS Code et le panneau Problèmes pour détecter les erreurs de syntaxe et les problèmes avant d'exécuter votre workflow, économisant du temps et évitant la frustration.

### Et Ensuite ?

Découvrez comment naviguer efficacement entre processus, modules et définitions dans des workflows complexes.

---

## 4. Navigation dans le Code et Gestion des Symboles

Une navigation efficace est cruciale lorsque vous travaillez avec des workflows complexes répartis sur plusieurs fichiers. Pour comprendre cela, remplacez la définition du processus dans `basic_workflow.nf` par une importation du module que nous vous avons fourni :

=== "Après"

    ```groovy title="basic_workflow.nf" linenums="3"
    include { FASTQC } from './modules/fastqc.nf'
    ```

=== "Avant"

    ```groovy title="basic_workflow.nf" linenums="3"
    process FASTQC {
        tag "${sample_id}"
        publishDir "${params.output_dir}/fastqc", mode: 'copy'

        input:
        tuple val(sample_id), path(reads)

        output:
        tuple val(sample_id), path("*.html"), emit: html
        tuple val(sample_id), path("*.zip"), emit: zip

        script:
        def args = task.ext.args ?: ''
        """
        fastqc \\
            ${args} \\
            --threads ${task.cpus} \\
            ${reads}
        """
    }
    ```

### 4.1. Aller à la Définition

Si vous passez la souris sur un nom de processus comme `FASTQC`, vous verrez une fenêtre contextuelle avec l'interface du module (entrées et sorties) :

![Aller à la définition](img/syntax.png)

Cette fonctionnalité est particulièrement précieuse lors de la création de workflows, car elle vous permet de comprendre l'interface du module sans ouvrir directement le fichier du module.

Vous pouvez naviguer rapidement vers n'importe quelle définition de processus, module ou variable en utilisant **Ctrl/Cmd-clic**. Passez la souris sur le lien vers le fichier du module en haut du script, et suivez le lien comme suggéré :

![Suivre le lien](img/follow_link.png)

La même chose fonctionne pour les noms de processus. Retournez à `basic_workflow.nf` et essayez ceci sur le nom du processus `FASTQC` dans le bloc workflow. Cela vous lie directement au nom du processus (qui est le même que le fichier du module dans cet exemple, mais pourrait être à mi-chemin d'un fichier beaucoup plus grand).

Pour revenir à l'endroit où vous étiez, utilisez **Alt+←** (ou **Ctrl+-** sur Mac). C'est une façon puissante d'explorer le code sans perdre votre position.

Explorons maintenant la navigation dans un workflow plus complexe en utilisant `complex_workflow.nf` (le fichier d'illustration uniquement mentionné précédemment). Ce workflow contient plusieurs processus définis dans des fichiers de module séparés, ainsi que certains en ligne. Bien que les structures multi-fichiers complexes puissent être difficiles à naviguer manuellement, la capacité à accéder aux définitions rend l'exploration beaucoup plus gérable.

1. Ouvrez `complex_workflow.nf`
2. Naviguez vers les définitions de modules
3. Utilisez **Alt+←** (ou **Ctrl+-**) pour revenir en arrière
4. Naviguez vers le nom du processus `FASTQC` dans le bloc workflow. Cela vous lie directement au nom du processus (qui est le même que le fichier du module dans cet exemple, mais pourrait être à mi-chemin d'un fichier beaucoup plus grand).
5. Naviguez à nouveau en arrière
6. Naviguez vers le processus `TRIM_GALORE` dans le bloc workflow. Celui-ci est défini en ligne, il ne vous mènera donc pas vers un fichier séparé, mais il vous montrera quand même la définition du processus, et vous pourrez toujours revenir à l'endroit où vous étiez.

### 4.2. Navigation par Symboles

Avec `complex_workflow.nf` toujours ouvert, vous pouvez obtenir un aperçu de tous les symboles du fichier en tapant `@` dans la barre de recherche en haut de VSCode (le raccourci clavier est `Ctrl/Cmd+Shift+O`, mais peut ne pas fonctionner dans Codespaces). Cela ouvre le panneau de navigation par symboles, qui liste tous les symboles du fichier actuel :

![Navigation par symboles](img/symbols.png)

Cela affiche :

- Toutes les définitions de processus
- Les définitions de workflow (il y a deux workflows définis dans ce fichier)
- Les définitions de fonctions

Commencez à taper pour filtrer les résultats.

### 4.3. Rechercher Toutes les Références

Comprendre où un processus ou une variable est utilisé dans votre base de code peut être très utile. Par exemple, si vous voulez trouver toutes les références au processus `FASTQC`, commencez par naviguer vers sa définition. Vous pouvez le faire en ouvrant `modules/fastqc.nf` directement, ou en utilisant la fonctionnalité de navigation rapide de VS Code avec `Ctrl/Cmd-clic` comme nous l'avons fait ci-dessus. Une fois à la définition du processus, faites un clic droit sur le nom du processus `FASTQC` et sélectionnez « Find All References » dans le menu contextuel pour voir toutes les instances où il est utilisé.

![Rechercher les références](img/references.png)

Cette fonctionnalité affiche toutes les instances où `FASTQC` est référencé dans votre espace de travail, y compris son utilisation dans les deux workflows distincts. Cette information est cruciale pour évaluer l'impact potentiel des modifications du processus `FASTQC`.

### 4.4. Panneau Plan

Le panneau Plan, situé dans la barre latérale de l'Explorateur (cliquez ![Icône Explorateur](img/files_icon.png)), fournit un aperçu pratique de tous les symboles de votre fichier actuel. Cette fonctionnalité vous permet de naviguer rapidement et de gérer la structure de votre code en affichant les fonctions, variables et autres éléments clés dans une vue hiérarchique.

![Panneau Plan](img/outline.png)

Utilisez le panneau Plan pour naviguer rapidement vers différentes parties de votre code sans utiliser l'explorateur de fichiers.

### 4.5. Visualisation DAG

L'extension Nextflow de VS Code peut visualiser votre workflow sous forme de graphe acyclique dirigé (DAG). Cela vous aide à comprendre le flux de données et les dépendances entre processus. Ouvrez `complex_workflow.nf` et cliquez sur le bouton « Preview DAG » au-dessus de `workflow {` (le deuxième bloc `workflow` dans ce fichier) :

![Aperçu DAG](img/dag_preview.png)

Il ne s'agit que du workflow « d'entrée », mais vous pouvez également prévisualiser le DAG pour les workflows internes en cliquant sur le bouton « Preview DAG » au-dessus du workflow `RNASEQ_PIPELINE {` plus haut :

![Aperçu DAG workflow interne](img/dag_preview_inner.png)

Pour ce workflow, vous pouvez utiliser les nœuds du DAG pour naviguer vers les définitions de processus correspondantes dans le code. Cliquez sur un nœud, et il vous mènera à la définition de processus pertinente dans l'éditeur. Particulièrement lorsqu'un workflow atteint une grande taille, cela peut vraiment vous aider à naviguer dans le code et à comprendre comment les processus sont connectés.

### Points Clés

Vous pouvez naviguer efficacement dans des workflows complexes en utilisant l'aller-à-la-définition, la recherche de symboles, la recherche de références et la visualisation DAG pour comprendre la structure du code et les dépendances.

### Et Ensuite ?

Découvrez comment travailler efficacement sur plusieurs fichiers interconnectés dans de plus grands projets Nextflow.

## 5. Travailler sur Plusieurs Fichiers

Le développement réel de Nextflow implique de travailler avec plusieurs fichiers interconnectés. Explorons comment VS Code vous aide à gérer efficacement des projets complexes.

### 5.1. Navigation Rapide entre Fichiers

Avec `complex_workflow.nf` ouvert, vous remarquerez qu'il importe plusieurs modules. Pratiquons la navigation rapide entre eux.

Appuyez sur **Ctrl+P** (ou **Cmd+P**) et commencez à taper « fast » :

VS Code vous montrera les fichiers correspondants. Sélectionnez `modules/fastqc.nf` pour y accéder instantanément. C'est beaucoup plus rapide que de cliquer dans l'explorateur de fichiers lorsque vous savez à peu près quel fichier vous recherchez.

Essayez ceci avec d'autres modèles :

- Tapez « star » pour trouver le fichier du module d'alignement STAR (`star.nf`)
- Tapez « utils » pour trouver le fichier des fonctions utilitaires (`utils.nf`)
- Tapez « config » pour accéder aux fichiers de configuration (`nextflow.config`)

### 5.2. Éditeur Divisé pour le Développement Multi-fichiers

Lorsque vous travaillez avec des modules, vous devez souvent voir à la fois le workflow principal et les définitions de modules simultanément. Configurons cela :

1. Ouvrez `complex_workflow.nf`
2. Ouvrez `modules/fastqc.nf` dans un nouvel onglet
3. Faites un clic droit sur l'onglet `modules/fastqc.nf` et sélectionnez « Split Right »
4. Maintenant vous pouvez voir les deux fichiers côte à côte

![Éditeur divisé](img/split_editor.png)

Ceci est inestimable lorsque :

- Vous vérifiez les interfaces de modules pendant l'écriture d'appels de workflow, et l'aperçu ne suffit pas
- Vous comparez des processus similaires dans différents modules
- Vous déboguez le flux de données entre workflow et modules

### 5.3. Recherche à l'Échelle du Projet

Parfois, vous devez trouver où des modèles spécifiques sont utilisés dans l'ensemble de votre projet. Appuyez sur `Ctrl/Cmd+Shift+F` pour ouvrir le panneau de recherche.

Essayez de rechercher `publishDir` dans l'espace de travail :

![Recherche de projet](img/project_search.png)

Cela vous montre tous les fichiers qui utilisent des répertoires de publication, vous aidant à :

- Comprendre les modèles d'organisation de sortie
- Trouver des exemples de directives spécifiques
- Assurer la cohérence entre les modules

### Points Clés

Vous pouvez gérer des projets multi-fichiers complexes en utilisant la navigation rapide entre fichiers, les éditeurs divisés et la recherche à l'échelle du projet pour travailler efficacement sur les workflows et les modules.

### Et Ensuite ?

Découvrez comment les fonctionnalités de formatage et de maintenance du code maintiennent vos workflows organisés et lisibles.

---

## 6. Formatage et Maintenance du Code

Un formatage correct du code est essentiel non seulement pour l'esthétique mais aussi pour améliorer la lisibilité, la compréhension et la facilité de mise à jour de workflows complexes.

### 6.1. Formatage Automatique en Action

Ouvrez `basic_workflow.nf` et dérangez délibérément le formatage :

- Supprimez quelques indentations : Mettez en surbrillance l'ensemble du document et appuyez plusieurs fois sur `shift+tab` pour supprimer autant d'indentations que possible.
- Ajoutez des espaces supplémentaires à des endroits aléatoires : dans l'instruction `channel.fromPath`, ajoutez 30 espaces après le `(`.
- Cassez maladroitement certaines lignes : Ajoutez une nouvelle ligne entre l'opérateur `.view {` et la chaîne `Processing sample:` mais n'ajoutez pas de nouvelle ligne correspondante avant la parenthèse fermante `}`.

Maintenant, appuyez sur `Shift+Alt+F` (ou `Shift+Option+F` sur MacOS) pour auto-formater :

VS Code immédiatement :

- Corrige l'indentation pour montrer clairement la structure du processus
- Aligne les éléments similaires de manière cohérente
- Supprime les espaces inutiles
- Maintient des sauts de ligne lisibles

Notez que le formatage automatique peut ne pas résoudre tous les problèmes de style de code. Le serveur de langage Nextflow vise à garder votre code propre, mais il respecte également vos préférences personnelles dans certains domaines. Par exemple, si vous supprimez l'indentation à l'intérieur du bloc `script` d'un processus, le formateur la laissera telle quelle, car vous pourriez intentionnellement préférer ce style.

Actuellement, il n'y a pas d'application stricte de style pour Nextflow, donc le serveur de langage offre une certaine flexibilité. Cependant, il appliquera systématiquement les règles de formatage autour des définitions de méthodes et de fonctions pour maintenir la clarté.

### 6.2. Fonctionnalités d'Organisation du Code

#### Commentaires Rapides

Sélectionnez un bloc de code dans votre workflow et appuyez sur **Ctrl+/** (ou **Cmd+/**) pour le commenter :

```groovy
// workflow {
//     ch_input = channel.fromPath(params.input)
//         .splitCsv(header: true)
//         .map { row -> [row.sample_id, file(row.fastq_path)] }
//
//     FASTQC(ch_input)
// }
```

C'est parfait pour :

- Désactiver temporairement des parties de workflows pendant le développement
- Ajouter des commentaires explicatifs à des opérations de canaux complexes
- Documenter des sections de workflow

Utilisez à nouveau **Ctrl+/** (ou **Cmd+/**) pour décommenter le code.

#### Pliage de Code pour Aperçu

Dans `complex_workflow.nf`, remarquez les petites flèches à côté des définitions de processus. Cliquez dessus pour plier (réduire) les processus :

![Pliage de code](img/code_folding.png)

Cela vous donne un aperçu de haut niveau de la structure de votre workflow sans vous perdre dans les détails d'implémentation.

#### Correspondance des Parenthèses

Placez votre curseur à côté de n'importe quelle parenthèse `{` ou `}` et VS Code met en évidence la parenthèse correspondante. Utilisez **Ctrl+Shift+\\** (ou **Cmd+Shift+\\**) pour sauter entre parenthèses correspondantes.

Ceci est crucial pour :

- Comprendre les limites des processus
- Trouver les parenthèses manquantes ou supplémentaires
- Naviguer dans les structures de workflow imbriquées

#### Sélection et Édition Multi-lignes

Pour éditer plusieurs lignes simultanément, VS Code offre de puissantes capacités multi-curseurs :

- **Sélection multi-lignes** : Maintenez **Ctrl+Alt** (ou **Cmd+Option** pour MacOS) et utilisez les touches fléchées pour sélectionner plusieurs lignes
- **Indentation multi-lignes** : Sélectionnez plusieurs lignes et utilisez **Tab** pour indenter ou **Shift+Tab** pour désindenter des blocs entiers

Ceci est particulièrement utile pour :

- Indenter des blocs de processus entiers de manière cohérente
- Ajouter des commentaires à plusieurs lignes à la fois
- Éditer des définitions de paramètres similaires dans plusieurs processus

### Points Clés

Vous pouvez maintenir un code propre et lisible en utilisant le formatage automatique, les fonctionnalités de commentaires, le pliage de code, la correspondance des parenthèses et l'édition multi-lignes pour organiser efficacement des workflows complexes.

### Et Ensuite ?

Découvrez comment VS Code s'intègre à votre workflow de développement plus large au-delà de la simple édition de code.

---

## 7. Intégration du Workflow de Développement

VS Code s'intègre bien avec votre workflow de développement au-delà de la simple édition de code.

### 7.1. Intégration du Contrôle de Version

!!! note "Codespaces et Intégration Git"

    Si vous travaillez dans **GitHub Codespaces**, certaines fonctionnalités d'intégration Git peuvent ne pas fonctionner comme prévu, en particulier les raccourcis clavier pour le contrôle de source. Vous avez peut-être également refusé d'ouvrir le répertoire en tant que dépôt Git lors de la configuration initiale, ce qui est acceptable à des fins de formation.

Si votre projet est un dépôt git (comme c'est le cas), VS Code affiche :

- Les fichiers modifiés avec des indicateurs colorés
- Le statut Git dans la barre d'état
- Des vues de différences en ligne
- Des capacités de commit et de push

Ouvrez le panneau Contrôle de Source en utilisant le bouton de contrôle de source (![Icône de contrôle de source](img/source_control_icon.png)) (`Ctrl+Shift+G` ou `Cmd+Shift+G` si vous travaillez avec VSCode localement) pour voir les modifications git et valider les commits directement dans l'éditeur.

![Panneau Contrôle de Source](img/source_control.png)

### 7.2. Exécution et Inspection des Workflows

Exécutons un workflow puis inspectons les résultats. Dans le terminal intégré (`Ctrl+Shift+` backtick sous Windows et MacOS), exécutez le workflow de base :

```bash title="Exécuter le workflow de base"
nextflow run basic_workflow.nf --input data/sample_data.csv --output_dir results
```

Pendant l'exécution du workflow, vous verrez la sortie en temps réel dans le terminal. Après l'achèvement, vous pouvez utiliser VS Code pour inspecter les résultats sans quitter votre éditeur :

1. **Naviguer vers les répertoires de travail** : Utilisez l'explorateur de fichiers ou le terminal pour parcourir `.nextflow/work`
2. **Ouvrir les fichiers journaux** : Cliquez sur les chemins de fichiers journaux dans la sortie du terminal pour les ouvrir directement dans VS Code
3. **Inspecter les sorties** : Parcourez les répertoires de résultats publiés dans l'explorateur de fichiers
4. **Voir les rapports d'exécution** : Ouvrez les rapports HTML directement dans VS Code ou votre navigateur

Cela garde tout au même endroit plutôt que de basculer entre plusieurs applications.

### Points Clés

Vous pouvez intégrer VS Code avec le contrôle de version et l'exécution de workflow pour gérer l'ensemble de votre processus de développement depuis une interface unique.

### Et Ensuite ?

Voyez comment toutes ces fonctionnalités IDE fonctionnent ensemble dans votre workflow de développement quotidien.

---

## 8. Récapitulatif et Notes Rapides

Voici quelques notes rapides sur chacune des fonctionnalités IDE discutées ci-dessus :

### 8.1. Commencer une Nouvelle Fonctionnalité

1. **Ouverture rapide de fichier** (`Ctrl+P` ou `Cmd+P`) pour trouver les modules existants pertinents
2. **Éditeur divisé** pour voir des processus similaires côte à côte
3. **Navigation par symboles** (`Ctrl+Shift+O` ou `Cmd+Shift+O`) pour comprendre la structure du fichier
4. **Auto-complétion** pour écrire rapidement du nouveau code

### 8.2. Débogage des Problèmes

1. **Panneau Problèmes** (`Ctrl+Shift+M` ou `Cmd+Shift+M`) pour voir toutes les erreurs en une fois
2. **Aller à la définition** (`Ctrl-clic` ou `Cmd-clic`) pour comprendre les interfaces de processus
3. **Rechercher toutes les références** pour voir comment les processus sont utilisés
4. **Recherche à l'échelle du projet** pour trouver des modèles ou des problèmes similaires

### 8.3. Refactorisation et Amélioration

1. **Recherche à l'échelle du projet** (`Ctrl+Shift+F` ou `Cmd+Shift+F`) pour trouver des modèles
2. **Auto-formatage** (`Shift+Alt+F` ou `Shift+Option+F`) pour maintenir la cohérence
3. **Pliage de code** pour se concentrer sur la structure
4. **Intégration Git** pour suivre les modifications

---

## Résumé

Vous avez maintenant eu un tour d'horizon rapide des fonctionnalités IDE de VS Code pour le développement Nextflow. Ces outils vous rendront nettement plus productif en :

- **Réduisant les erreurs** grâce à la vérification syntaxique en temps réel
- **Accélérant le développement** avec l'auto-complétion intelligente
- **Améliorant la navigation** dans des workflows multi-fichiers complexes
- **Maintenant la qualité** grâce à un formatage cohérent
- **Améliorant la compréhension** grâce à la coloration avancée et à la visualisation de structure

Nous ne nous attendons pas à ce que vous vous souveniez de tout, mais maintenant que vous savez que ces fonctionnalités existent, vous pourrez les trouver quand vous en aurez besoin. À mesure que vous continuerez à développer des workflows Nextflow, ces fonctionnalités IDE deviendront une seconde nature, vous permettant de vous concentrer sur l'écriture de code de haute qualité plutôt que de vous battre avec la syntaxe et la structure.

### Et Ensuite ?

Appliquez ces compétences IDE en travaillant sur d'autres modules de formation, par exemple :

- **[nf-test](nf-test.md)** : Créez des suites de tests complètes pour vos workflows
- **[Hello nf-core](../../hello_nf-core/)** : Construisez des pipelines de qualité production avec les standards de la communauté

Le véritable pouvoir de ces fonctionnalités IDE émerge à mesure que vous travaillez sur des projets plus grands et plus complexes. Commencez à les incorporer progressivement dans votre workflow - en quelques sessions, elles deviendront une seconde nature et transformeront votre approche du développement Nextflow.

De la détection d'erreurs avant qu'elles ne vous ralentissent à la navigation dans des bases de code complexes avec facilité, ces outils feront de vous un·e développeur·se plus confiant·e et efficace.

Bon codage !
