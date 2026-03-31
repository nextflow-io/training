# Partie 1 : Les bases des plugins

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Traduction assistée par IA - [en savoir plus et suggérer des améliorations](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

Dans cette section, vous apprendrez comment les plugins étendent Nextflow, puis vous essaierez trois plugins différents pour les voir en action.

---

## 1. Comment fonctionnent les plugins

Les plugins étendent Nextflow via plusieurs types d'extensions :

| Type d'extension      | Ce qu'il fait                                              | Exemple                      |
| --------------------- | ---------------------------------------------------------- | ---------------------------- |
| Fonctions             | Ajoute des fonctions personnalisées appelables depuis les workflows | `samplesheetToList()`        |
| Moniteurs de workflow | Réagit aux événements comme la complétion d'une tâche      | Journalisation personnalisée, alertes Slack |
| Executors             | Ajoute des backends d'exécution de tâches                  | AWS Batch, Kubernetes        |
| Systèmes de fichiers  | Ajoute des backends de stockage                            | S3, Azure Blob               |

Les fonctions et les moniteurs de workflow (appelés "trace observers" dans l'API Nextflow) sont les types les plus courants pour les auteur·trices de plugins.
Les executors et les systèmes de fichiers sont généralement créés par des éditeurs de plateformes.

Les exercices suivants vous montrent des plugins de fonctions et un plugin observateur, afin que vous puissiez voir les deux types en action.

---

## 2. Utiliser des plugins de fonctions

Les plugins de fonctions ajoutent des fonctions appelables que vous importez dans vos workflows.
Vous en essaierez deux : nf-hello (un exemple simple) et nf-schema (un plugin réel largement utilisé).
Les deux exercices modifient le même pipeline `hello.nf`, afin que vous puissiez voir comment les plugins améliorent un workflow existant.

### 2.1. nf-hello : remplacer du code écrit à la main

Le plugin [nf-hello](https://github.com/nextflow-io/nf-hello) fournit une fonction `randomString` qui génère des chaînes de caractères aléatoires.
Le pipeline définit déjà sa propre version inline de cette fonction, que vous remplacerez par celle du plugin.

#### 2.1.1. Voir le point de départ

Examinez le pipeline :

```bash
cat hello.nf
```

```groovy title="Output"
#!/usr/bin/env nextflow

params.input = 'greetings.csv'

/**
 * Génère une chaîne alphanumérique aléatoire
 */
def randomString(int length) {
    def chars = ('a'..'z') + ('A'..'Z') + ('0'..'9')
    def random = new Random()
    return (1..length).collect { chars[random.nextInt(chars.size())] }.join()
}

process SAY_HELLO {
    input:
        val greeting
    output:
        stdout
    script:
    """
    echo '$greeting'
    """
}

workflow {
    greeting_ch = channel.fromPath(params.input)
        .splitCsv(header: true)
        .map { row -> "${row.greeting}_${randomString(8)}" }
    SAY_HELLO(greeting_ch)
    SAY_HELLO.out.view { result -> "Output: ${result.trim()}" }
}
```

Le pipeline définit sa propre fonction `randomString` inline, puis l'utilise pour ajouter un identifiant aléatoire à chaque message de salutation.

Exécutez-le :

```bash
nextflow run hello.nf
```

```console title="Output"
Output: Hello_aBcDeFgH
Output: Bonjour_xYzWvUtS
Output: Holà_qRsPdMnK
Output: Ciao_jLhGfEcB
Output: Hallo_tNwOiAuR
```

L'ordre de votre sortie et les chaînes aléatoires seront différents, et si vous exécutez le script à nouveau, vous obtiendrez un ensemble différent de salutations aléatoires.

#### 2.1.2. Configurer le plugin

Remplacez la fonction inline par une fonction provenant du plugin. Ajoutez ceci à votre `nextflow.config` :

```groovy title="nextflow.config"
// Configuration pour les exercices de développement de plugins
plugins {
    id 'nf-hello@0.5.0'
}
```

Les plugins sont déclarés dans `nextflow.config` en utilisant le bloc `plugins {}`.
Nextflow les télécharge automatiquement depuis le [Nextflow Plugin Registry](https://registry.nextflow.io/), un dépôt central de plugins communautaires et officiels.

#### 2.1.3. Utiliser la fonction du plugin

Remplacez la fonction `randomString` inline par la version du plugin :

=== "Après"

    ```groovy title="hello.nf" hl_lines="3"
    #!/usr/bin/env nextflow

    include { randomString } from 'plugin/nf-hello'

    params.input = 'greetings.csv'

    process SAY_HELLO {
        input:
            val greeting
        output:
            stdout
        script:
        """
        echo '$greeting'
        """
    }

    workflow {
        greeting_ch = channel.fromPath(params.input)
            .splitCsv(header: true)
            .map { row -> "${row.greeting}_${randomString(8)}" }
        SAY_HELLO(greeting_ch)
        SAY_HELLO.out.view { result -> "Output: ${result.trim()}" }
    }
    ```

=== "Avant"

    ```groovy title="hello.nf" hl_lines="5-12"
    #!/usr/bin/env nextflow

    params.input = 'greetings.csv'

    /**
     * Génère une chaîne alphanumérique aléatoire
     */
    def randomString(int length) {
        def chars = ('a'..'z') + ('A'..'Z') + ('0'..'9')
        def random = new Random()
        return (1..length).collect { chars[random.nextInt(chars.size())] }.join()
    }

    process SAY_HELLO {
        input:
            val greeting
        output:
            stdout
        script:
        """
        echo '$greeting'
        """
    }

    workflow {
        greeting_ch = channel.fromPath(params.input)
            .splitCsv(header: true)
            .map { row -> "${row.greeting}_${randomString(8)}" }
        SAY_HELLO(greeting_ch)
        SAY_HELLO.out.view { result -> "Output: ${result.trim()}" }
    }
    ```

L'instruction `include` importe `randomString` depuis une bibliothèque éprouvée, testée et maintenue par un plus grand nombre de contributeur·trices qui peuvent détecter et corriger les bugs.
Au lieu que chaque pipeline maintienne sa propre copie de la fonction, chaque pipeline qui utilise le plugin bénéficie de la même implémentation validée.
Cela réduit la duplication de code et la charge de maintenance qui en découle.
La syntaxe `#!groovy include { function } from 'plugin/plugin-id'` est le même `include` utilisé pour les modules Nextflow, avec un préfixe `plugin/`.
Vous pouvez consulter le [code source de `randomString`](https://github.com/nextflow-io/nf-hello/blob/e67bddebfa589c7ae51f41bf780c92068dc09e93/plugins/nf-hello/src/main/nextflow/hello/HelloExtension.groovy#L110) dans le dépôt nf-hello sur GitHub.

#### 2.1.4. Exécuter le pipeline

```bash
nextflow run hello.nf
```

```console title="Output"
Pipeline is starting! 🚀
Output: Hello_yqvtclcc
Output: Bonjour_vwwpyzcs
Output: Holà_wrghmgab
Output: Ciao_noniajuy
Output: Hallo_tvrtuxtp
Pipeline complete! 👋
```

(Vos chaînes aléatoires seront différentes.)

La sortie contient toujours des suffixes aléatoires, mais maintenant `randomString` provient du plugin nf-hello plutôt que du code inline.
Les messages "Pipeline is starting !" et "Pipeline complete !" sont nouveaux.
Ils proviennent du composant observateur du plugin, que vous explorerez dans la Partie 5.

Nextflow télécharge les plugins automatiquement la première fois qu'ils sont utilisés, de sorte que tout pipeline qui déclare `nf-hello@0.5.0` obtient exactement la même fonction `randomString` testée sans copier de code entre les projets.

Vous avez maintenant vu les trois étapes pour utiliser un plugin de fonctions : le déclarer dans `nextflow.config`, importer la fonction avec `include`, et l'appeler dans votre workflow.
L'exercice suivant applique ces mêmes étapes à un plugin réel.

### 2.2. nf-schema : analyse CSV validée

Le plugin [nf-schema](https://github.com/nextflow-io/nf-schema) est l'un des plugins Nextflow les plus largement utilisés.
Il fournit `samplesheetToList`, une fonction qui analyse les fichiers CSV/TSV en utilisant un schéma JSON définissant les colonnes et les types attendus.

Le pipeline lit actuellement `greetings.csv` en utilisant `splitCsv` et un `map` manuel, mais nf-schema peut remplacer cela par une analyse validée et pilotée par schéma.
Un fichier de schéma JSON (`greetings_schema.json`) est déjà fourni dans le répertoire de l'exercice.

??? info "Qu'est-ce qu'un schéma ?"

    Un schéma est une description formelle de ce à quoi ressemblent des données valides.
    Il définit des éléments tels que les colonnes attendues, le type de chaque valeur (chaîne de caractères, nombre, etc.) et les champs obligatoires.

    Considérez-le comme un contrat : si les données d'entrée ne correspondent pas au schéma, l'outil peut détecter le problème tôt plutôt que de laisser celui-ci provoquer des erreurs confuses plus tard dans le pipeline.

#### 2.2.1. Examiner le schéma

```bash
cat greetings_schema.json
```

```json title="Output"
{
  "$schema": "https://json-schema.org/draft/2020-12/schema",
  "type": "array",
  "items": {
    "type": "object",
    "properties": {
      "greeting": {
        "type": "string",
        "description": "The greeting text"
      },
      "language": {
        "type": "string",
        "description": "The language of the greeting"
      }
    },
    "required": ["greeting"]
  }
}
```

Le schéma définit deux colonnes (`greeting` et `language`) et marque `greeting` comme obligatoire.
Si quelqu'un fournit un CSV sans la colonne `greeting`, nf-schema détecte l'erreur avant l'exécution du pipeline.

#### 2.2.2. Ajouter nf-schema à la configuration

Mettez à jour `nextflow.config` pour inclure les deux plugins :

=== "Après"

    ```groovy title="nextflow.config" hl_lines="3"
    plugins {
        id 'nf-hello@0.5.0'
        id 'nf-schema@2.6.1'
    }
    ```

=== "Avant"

    ```groovy title="nextflow.config"
    plugins {
        id 'nf-hello@0.5.0'
    }
    ```

#### 2.2.3. Mettre à jour hello.nf pour utiliser samplesheetToList

Remplacez l'entrée `splitCsv` par `samplesheetToList` :

=== "Après"

    ```groovy title="hello.nf" hl_lines="4 20 21 22"
    #!/usr/bin/env nextflow

    include { randomString } from 'plugin/nf-hello'
    include { samplesheetToList } from 'plugin/nf-schema'

    params.input = 'greetings.csv'

    process SAY_HELLO {
        input:
            val greeting
        output:
            stdout
        script:
        """
        echo '$greeting'
        """
    }

    workflow {
        def samplesheet_list = samplesheetToList(params.input, 'greetings_schema.json')
        greeting_ch = Channel.fromList(samplesheet_list)
            .map { row -> "${row[0]}_${randomString(8)}" }
        SAY_HELLO(greeting_ch)
        SAY_HELLO.out.view { result -> "Output: ${result.trim()}" }
    }
    ```

=== "Avant"

    ```groovy title="hello.nf" hl_lines="19 20 21"
    #!/usr/bin/env nextflow

    include { randomString } from 'plugin/nf-hello'

    params.input = 'greetings.csv'

    process SAY_HELLO {
        input:
            val greeting
        output:
            stdout
        script:
        """
        echo '$greeting'
        """
    }

    workflow {
        greeting_ch = channel.fromPath(params.input)
            .splitCsv(header: true)
            .map { row -> "${row.greeting}_${randomString(8)}" }
        SAY_HELLO(greeting_ch)
        SAY_HELLO.out.view { result -> "Output: ${result.trim()}" }
    }
    ```

Le code d'analyse personnalisé avec `splitCsv` et `map` est remplacé par `samplesheetToList`, une fonction éprouvée et testée qui valide également la samplesheet par rapport au schéma avant l'exécution du pipeline.
Cela réduit la charge de maintenance de la logique d'analyse écrite à la main tout en améliorant l'expérience des utilisateur·trices du pipeline, qui reçoivent des messages d'erreur clairs lorsque leurs données d'entrée ne correspondent pas au format attendu.
Chaque ligne devient une liste de valeurs dans l'ordre des colonnes, donc `row[0]` est le message de salutation et `row[1]` est la langue.

#### 2.2.4. Exécuter le pipeline

```bash
nextflow run hello.nf
```

```console title="Output"
Pipeline is starting! 🚀
Output: Hello_diozjdwm
Output: Bonjour_speathmm
Output: Holà_dllxnzap
Output: Ciao_wzueddzc
Output: Hallo_hsxwrjbh
Pipeline complete! 👋
```

(Vos chaînes aléatoires seront différentes.)

La sortie est identique, mais maintenant le schéma valide la structure du CSV avant l'exécution du pipeline.
Dans des pipelines réels avec des samplesheets complexes et de nombreuses colonnes, ce type de validation prévient des erreurs que `splitCsv` + `map` manuels laisseraient passer.

#### 2.2.5. Voir la validation en action

Pour voir ce que la validation par schéma détecte, essayez d'introduire des erreurs dans `greetings.csv`.

Renommez la colonne obligatoire `greeting` en `message` :

```csv title="greetings.csv" hl_lines="1"
message,language
Hello,English
Bonjour,French
Holà,Spanish
Ciao,Italian
Hallo,German
```

Exécutez le pipeline :

```bash
nextflow run hello.nf
```

```console title="Output"
ERROR ~ Validation of samplesheet failed!

The following errors have been detected in greetings.csv:

-> Entry 1: Missing required field(s): greeting
-> Entry 2: Missing required field(s): greeting
-> Entry 3: Missing required field(s): greeting
-> Entry 4: Missing required field(s): greeting
-> Entry 5: Missing required field(s): greeting
```

Le pipeline refuse de s'exécuter car le schéma exige une colonne `greeting` et ne peut pas en trouver une.

Maintenant, restaurez la colonne obligatoire mais renommez la colonne optionnelle `language` en `lang` :

```csv title="greetings.csv" hl_lines="1"
greeting,lang
Hello,English
Bonjour,French
Holà,Spanish
Ciao,Italian
Hallo,German
```

```bash
nextflow run hello.nf
```

Cette fois, le pipeline s'exécute, mais affiche un avertissement :

```console title="Output (partial)"
WARN: Found the following unidentified headers in greetings.csv:
	- lang
```

Les colonnes obligatoires provoquent des erreurs bloquantes ; les colonnes optionnelles provoquent des avertissements.
C'est le type de retour précoce qui fait gagner du temps de débogage dans des pipelines réels avec des dizaines de colonnes.

#### 2.2.6. Configurer le comportement de validation

L'avertissement concernant `lang` est utile, mais vous pouvez contrôler sa sévérité via la configuration.
Les plugins peuvent inclure leurs propres portées de configuration qui contrôlent leur comportement.
Le plugin nf-schema inclut la portée de configuration `validation` ; en modifiant les paramètres ici, vous pouvez changer le comportement de nf-schema.

Ajoutez un bloc `validation` à `nextflow.config` pour que les en-têtes non reconnus provoquent une erreur plutôt qu'un avertissement :

=== "Après"

    ```groovy title="nextflow.config" hl_lines="6-10"
    plugins {
        id 'nf-hello@0.5.0'
        id 'nf-schema@2.6.1'
    }

    validation {
        logging {
            unrecognisedHeaders = "error"
        }
    }
    ```

=== "Avant"

    ```groovy title="nextflow.config"
    plugins {
        id 'nf-hello@0.5.0'
        id 'nf-schema@2.6.1'
    }
    ```

Exécutez à nouveau le pipeline avec la colonne `lang` toujours en place :

```bash
nextflow run hello.nf
```

```console title="Output (partial)"
Found the following unidentified headers in greetings.csv:
	- lang
 -- Check script 'hello.nf' at line: 20 or see '.nextflow.log' file for more details
```

Le pipeline échoue maintenant au lieu d'afficher un avertissement.
Le code du pipeline n'a pas changé ; seule la configuration a été modifiée.

Restaurez `greetings.csv` à son état d'origine et supprimez le bloc `validation` avant de continuer :

```csv title="greetings.csv"
greeting,language
Hello,English
Bonjour,French
Holà,Spanish
Ciao,Italian
Hallo,German
```

```groovy title="nextflow.config"
plugins {
    id 'nf-hello@0.5.0'
    id 'nf-schema@2.6.1'
}
```

nf-hello et nf-schema sont tous deux des plugins de fonctions : ils fournissent des fonctions que vous importez avec `include` et appelez dans votre code de workflow.
L'exercice suivant présente un type de plugin différent qui fonctionne sans aucune instruction `include`.

---

## 3. Utiliser un plugin observateur : nf-co2footprint

Tous les plugins ne fournissent pas des fonctions à importer.
Le plugin [nf-co2footprint](https://github.com/nextflow-io/nf-co2footprint) utilise un **trace observer** pour surveiller l'utilisation des ressources de votre pipeline et estimer son empreinte carbone.
Vous n'avez pas besoin de modifier le code du pipeline ; il suffit de l'ajouter à la configuration.

### 3.1. Ajouter nf-co2footprint à la configuration

Mettez à jour `nextflow.config` :

=== "Après"

    ```groovy title="nextflow.config" hl_lines="4"
    plugins {
        id 'nf-hello@0.5.0'
        id 'nf-schema@2.6.1'
        id 'nf-co2footprint@1.2.0'
    }
    ```

=== "Avant"

    ```groovy title="nextflow.config"
    plugins {
        id 'nf-hello@0.5.0'
        id 'nf-schema@2.6.1'
    }
    ```

### 3.2. Exécuter le pipeline

```bash
nextflow run hello.nf
```

Le plugin produit plusieurs messages INFO et WARN pendant l'exécution.
Ceux-ci sont normaux pour un petit exemple s'exécutant sur une machine locale :

```console title="Output (partial)"
nf-co2footprint plugin  ~  version 1.2.0
WARN - [nf-co2footprint] Target zone null not found. Attempting to retrieve carbon intensity for fallback zone GLOBAL.
INFO - [nf-co2footprint] Using fallback carbon intensity from GLOBAL from CI table: 480.0 gCO₂eq/kWh.
WARN - [nf-co2footprint] Executor 'null' not mapped.
WARN - [nf-co2footprint] Fallback to: `machineType = null`, `pue = 1.0`. ...
...
WARN - [nf-co2footprint] No CPU model detected. Using default CPU power draw value (11.41 W).
WARN - [nf-co2footprint] 🔁 Requested memory is null for task 2. Using maximum consumed memory/`peak_rss` (0 GB) for CO₂e footprint computation.
```

Les avertissements concernant la zone, l'executor, le modèle de CPU et la mémoire apparaissent parce que le plugin ne peut pas détecter tous les détails matériels d'un environnement de formation local.
Dans un environnement de production (par exemple, un cluster HPC ou le cloud), ces valeurs seraient disponibles et les estimations plus précises.

À la fin, recherchez une ligne comme :

```console title="Output (partial)"
🌱 The workflow run used 126.76 uWh of electricity, resulting in the release of 60.84 ug of CO₂ equivalents into the atmosphere.
```

(Vos chiffres seront différents.)

### 3.3. Consulter le rapport

Le plugin génère des fichiers de sortie dans votre répertoire de travail :

```bash
ls co2footprint_*
```

```console title="Output"
co2footprint_report_<timestamp>.html
co2footprint_summary_<timestamp>.txt
co2footprint_trace_<timestamp>.txt
```

Examinez le résumé :

```bash
cat co2footprint_summary_*.txt
```

```console title="Output"
Total CO₂e footprint measures of this workflow run (including cached tasks):
  CO₂e emissions: 60.84 ug
  Energy consumption: 126.76 uWh
  CO₂e emissions (market): -

Which equals:
  - 3.48E-7 km travelled by car
  - It takes one tree 0.17s to sequester the equivalent amount of CO₂ from the atmosphere
  - 1.22E-7 % of a flight from Paris to London
```

(Vos chiffres seront différents.)

La première section affiche les chiffres bruts d'énergie et d'émissions.
La section "Which equals" met ces chiffres en perspective en les convertissant en équivalents familiers.
Le résumé inclut également une section listant les options de configuration du plugin et une citation de l'article de recherche [Green Algorithms](https://doi.org/10.1002/advs.202100707) sur lequel la méthode de calcul est basée.

### 3.4. Configurer le plugin

L'avertissement "Target zone null" de la section 3.2 est apparu parce qu'aucun emplacement n'était configuré pour le plugin.
Le plugin nf-co2footprint définit une portée de configuration `co2footprint` où vous pouvez définir votre emplacement géographique.

Ajoutez un bloc `co2footprint` à `nextflow.config` :

=== "Après"

    ```groovy title="nextflow.config" hl_lines="7-9"
    plugins {
        id 'nf-hello@0.5.0'
        id 'nf-schema@2.6.1'
        id 'nf-co2footprint@1.2.0'
    }

    co2footprint {
        location = 'GB'
    }
    ```

=== "Avant"

    ```groovy title="nextflow.config"
    plugins {
        id 'nf-hello@0.5.0'
        id 'nf-schema@2.6.1'
        id 'nf-co2footprint@1.2.0'
    }
    ```

!!! tip "Astuce"

    Utilisez le code de votre propre pays si vous préférez (par exemple, `'US'`, `'DE'`, `'FR'`).

Exécutez le pipeline :

```bash
nextflow run hello.nf
```

```console title="Output (partial)"
INFO - [nf-co2footprint] Using fallback carbon intensity from GB from CI table: 163.92 gCO₂eq/kWh.
```

L'avertissement de zone a disparu.
Le plugin utilise maintenant l'intensité carbone spécifique à la Grande-Bretagne (163,92 gCO₂eq/kWh) au lieu de la valeur de repli mondiale (480,0 gCO₂eq/kWh).

!!! note "Note"

    Vous pouvez également voir un message `WARN: Unrecognized config option 'co2footprint.location'`.
    Il s'agit d'un message cosmétique qui peut être ignoré en toute sécurité ; le plugin lit toujours la valeur correctement.

Dans la Partie 6, vous créerez une portée de configuration pour votre propre plugin.

Ce plugin fonctionne entièrement via le mécanisme observateur, en se connectant aux événements du cycle de vie du workflow pour collecter des métriques de ressources et générer son rapport lorsque le pipeline se termine.

Vous avez maintenant essayé des plugins de fonctions (importés avec `include`) et un plugin observateur (activé via la configuration seule).
Ce sont les deux types d'extensions les plus courants, mais comme le montre le tableau de la section 1, les plugins peuvent également ajouter des executors et des systèmes de fichiers.

---

## 4. Découvrir les plugins

Le [Nextflow Plugin Registry](https://registry.nextflow.io/) est le hub central pour trouver les plugins disponibles.

![La page du plugin nf-hello sur registry.nextflow.io](img/plugin-registry-nf-hello.png)

Chaque page de plugin affiche sa description, les versions disponibles, les instructions d'installation et des liens vers la documentation.

---

## 5. Se préparer au développement de plugins

Les sections suivantes (Parties 2 à 6) utilisent un fichier de pipeline séparé, `greet.nf`, qui dépend de nf-schema mais pas de nf-hello ni de nf-co2footprint.

Mettez à jour `nextflow.config` pour ne conserver que nf-schema :

```groovy title="nextflow.config"
// Configuration pour les exercices de développement de plugins
plugins {
    id 'nf-schema@2.6.1'
}
```

Supprimez les fichiers de sortie de co2footprint :

```bash
rm -f co2footprint_*
```

Le fichier `hello.nf` conserve votre travail de la Partie 1 pour référence ; à partir de maintenant, vous travaillerez avec `greet.nf`.

---

## À retenir

Vous avez utilisé trois plugins différents :

- **nf-hello** : Un plugin de fonctions fournissant `randomString`, importé avec `include`
- **nf-schema** : Un plugin de fonctions fournissant `samplesheetToList` pour l'analyse CSV validée par schéma
- **nf-co2footprint** : Un plugin observateur qui surveille automatiquement l'utilisation des ressources, sans `include` nécessaire

Modèles clés :

- Les plugins sont déclarés dans `nextflow.config` avec `#!groovy plugins { id 'plugin-name@version' }`
- Les plugins de fonctions nécessitent `#!groovy include { function } from 'plugin/plugin-id'`
- Les plugins observateurs fonctionnent automatiquement une fois déclarés dans la configuration
- Les plugins peuvent définir des portées de configuration (par exemple, `#!groovy validation {}`, `#!groovy co2footprint {}`) pour personnaliser leur comportement
- Le [Nextflow Plugin Registry](https://registry.nextflow.io/) liste les plugins disponibles

---

## Et ensuite ?

Les sections suivantes vous montrent comment créer votre propre plugin.
Si le développement de plugins ne vous intéresse pas, vous pouvez vous arrêter ici ou passer directement au [Résumé](summary.md).

[Continuer vers la Partie 2 :material-arrow-right:](02_create_project.md){ .md-button .md-button--primary }
