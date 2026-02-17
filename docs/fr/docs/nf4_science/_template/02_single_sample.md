# Partie 2 : Traitement d'un échantillon unique

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Traduction assistée par IA - [en savoir plus et suggérer des améliorations](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

Dans la Partie 1, vous avez testé les commandes {TOOL_A} et {TOOL_B} manuellement dans leurs conteneurs respectifs.
Nous allons maintenant encapsuler ces mêmes commandes dans un workflow Nextflow.

## Objectif

Dans cette partie du cours, nous allons développer un workflow qui effectue les opérations suivantes :

1. {PROCESS_1_DESCRIPTION}
2. {PROCESS_2_DESCRIPTION}

<figure class="excalidraw">
--8<-- "docs/en/docs/nf4_science/{DOMAIN_DIR}/img/{DIAGRAM_1}.svg"
</figure>

Cela reproduit les étapes de la Partie 1, où vous avez exécuté ces commandes manuellement dans leurs conteneurs.

Comme point de départ, nous vous fournissons un fichier de workflow, `{DOMAIN_DIR}.nf`, qui décrit les principales parties du workflow, ainsi que deux fichiers de modules, {TOOL_A_MODULE}.nf et {TOOL_B_MODULE}.nf, qui décrivent la structure des modules.
Ces fichiers ne sont pas fonctionnels ; leur objectif est simplement de servir de squelettes que vous compléterez avec les parties intéressantes du code.

## Plan de la leçon

Afin de rendre le processus de développement plus pédagogique, nous l'avons décomposé en {N} étapes :

1. **Écrire un workflow à une seule étape qui exécute {TOOL_A_ACTION}.**
   Cela couvre la création d'un module, son importation et son appel dans un workflow.
2. **Ajouter un second processus pour exécuter {TOOL_B_ACTION}.**
   Cela introduit le chaînage des sorties de processus vers les entrées et la gestion des fichiers accessoires.
3. **Adapter le workflow pour qu'il s'exécute sur un lot d'échantillons.**
   Cela couvre l'exécution parallèle et introduit les tuples pour maintenir les fichiers associés ensemble.
4. **Faire en sorte que le workflow accepte une feuille d'échantillons contenant un lot de fichiers d'entrée.**
   Cela démontre un modèle courant pour fournir des entrées en masse.

Chaque étape se concentre sur un aspect spécifique du développement de workflow.

---

## 1. Écrire un workflow à une seule étape qui exécute {TOOL_A_ACTION}

Cette première étape se concentre sur les bases : charger {PRIMARY_INPUT_TYPE} et {TOOL_A_OUTPUT_DESCRIPTION}.

Rappelez-vous la commande `{TOOL_A_COMMAND_NAME}` de la [Partie 1](01_method.md) :

```bash
{TOOL_A_COMMAND_SUMMARY}
```

La commande prend {INPUT_DESCRIPTION} et produit {OUTPUT_DESCRIPTION}.
L'URI du conteneur était `{TOOL_A_CONTAINER_URI}`.

Nous allons prendre ces informations et les encapsuler dans Nextflow en trois étapes :

1. Configurer l'entrée
2. Écrire le processus et l'appeler dans le workflow
3. Configurer la gestion de la sortie

### 1.1. Configurer l'entrée

Nous devons déclarer un paramètre d'entrée, créer un profil de test pour fournir une valeur par défaut pratique, et créer un canal d'entrée.

#### 1.1.1. Ajouter une déclaration de paramètre d'entrée

Dans le fichier de workflow principal `{DOMAIN_DIR}.nf`, sous la section `Pipeline parameters`, déclarez un paramètre CLI appelé `{PRIMARY_PARAM_NAME}`.

=== "Après"

    ```groovy title="{DOMAIN_DIR}.nf" linenums="5" hl_lines="4-7"
    /*
     * Pipeline parameters
     */
    params {
        // Primary input
        {PRIMARY_PARAM_NAME}: Path
    }
    ```

=== "Avant"

    ```groovy title="{DOMAIN_DIR}.nf" linenums="5"
    /*
     * Pipeline parameters
     */

    // Primary input
    ```

Cela configure le paramètre CLI, mais nous ne voulons pas saisir le chemin du fichier à chaque fois que nous exécutons le workflow pendant le développement.
Il existe plusieurs options pour fournir une valeur par défaut ; ici, nous utilisons un profil de test.

#### 1.1.2. Créer un profil de test avec une valeur par défaut dans `nextflow.config`

Un profil de test fournit des valeurs par défaut pratiques pour essayer un workflow sans spécifier d'entrées sur la ligne de commande.
C'est une convention courante dans l'écosystème Nextflow (voir [Hello Config](../../hello_nextflow/06_hello_config.md) pour plus de détails).

Ajoutez un bloc `profiles` à `nextflow.config` avec un profil `test` qui définit le paramètre `{PRIMARY_PARAM_NAME}` sur l'un des {PRIMARY_INPUT_TYPE}s de test.

=== "Après"

    ```groovy title="nextflow.config" linenums="1" hl_lines="3-7"
    docker.enabled = true

    profiles {
        test {
            params.{PRIMARY_PARAM_NAME} = "${projectDir}/data/{TEST_INPUT_PATH}"
        }
    }
    ```

=== "Avant"

    ```groovy title="nextflow.config" linenums="1"
    docker.enabled = true
    ```

Ici, nous utilisons `${projectDir}`, une variable Nextflow intégrée qui pointe vers le répertoire où se trouve le script de workflow.
Cela facilite la référence aux fichiers de données et autres ressources sans coder en dur des chemins absolus.

#### 1.1.3. Configurer le canal d'entrée

{INPUT_CHANNEL_INSTRUCTIONS}

### 1.2. Écrire le module {TOOL_A_NAME}

{MODULE_INSTRUCTIONS_COVERING: container, input, output, script}

### 1.3. Importer et appeler le module dans le workflow

{IMPORT_AND_CALL_INSTRUCTIONS}

### 1.4. Exécuter le workflow

À ce stade, nous avons un workflow en une étape qui devrait être entièrement fonctionnel.

Nous pouvons l'exécuter avec `-profile test` pour utiliser la valeur par défaut configurée dans le profil de test et éviter d'avoir à écrire le chemin sur la ligne de commande.

```bash
nextflow run {DOMAIN_DIR}.nf -profile test
```

??? success "Sortie de la commande"

    ```console
    N E X T F L O W   ~  version 25.10.2

    ┃ Launching `{DOMAIN_DIR}.nf` [{RUN_NAME}] DSL2 - revision: {HASH}

    executor >  local (1)
    [{HASH}] {PROCESS_A_NAME} (1) | 1 of 1 ✔
    ```

{VERIFY_OUTPUT_INSTRUCTIONS}

### À retenir

Vous savez comment créer un module contenant un processus, l'importer dans un workflow, l'appeler avec un canal d'entrée et publier les résultats.

### Et ensuite ?

Ajoutez un second processus pour chaîner des étapes d'analyse supplémentaires.

---

## 2. Ajouter un second processus pour exécuter {TOOL_B_ACTION}

{DESCRIBE_WHAT_THIS_STEP_ADDS}

Rappelez-vous la commande `{TOOL_B_COMMAND_NAME}` de la [Partie 1](01_method.md) :

```bash
{TOOL_B_COMMAND_SUMMARY}
```

### 2.1. Écrire le module {TOOL_B_NAME}

{MODULE_INSTRUCTIONS}

### 2.2. Importer et appeler le module dans le workflow

{IMPORT_AND_CALL_INSTRUCTIONS_CHAINING_OUTPUTS}

### 2.3. Exécuter le workflow

```bash
nextflow run {DOMAIN_DIR}.nf -profile test
```

??? success "Sortie de la commande"

    ```console
    {EXPECTED_OUTPUT}
    ```

{VERIFY_OUTPUT_INSTRUCTIONS}

### À retenir

Vous savez comment chaîner les sorties de processus vers les entrées et gérer les fichiers accessoires dans le workflow.

### Et ensuite ?

Adaptez le workflow pour traiter plusieurs échantillons en parallèle.

---

## 3. Adapter le workflow pour qu'il s'exécute sur un lot d'échantillons

Jusqu'à présent, le workflow traite un seul échantillon.
Pour gérer plusieurs échantillons, nous devons modifier la façon dont les entrées sont fournies et exploiter le paradigme de flux de données de Nextflow pour paralléliser l'exécution.

### 3.1. {BATCH_INPUT_INSTRUCTIONS}

{INSTRUCTIONS_FOR_SWITCHING_TO_SAMPLESHEET_OR_BATCH_INPUT}

### 3.2. Exécuter le workflow sur plusieurs échantillons

```bash
nextflow run {DOMAIN_DIR}.nf -profile test
```

??? success "Sortie de la commande"

    ```console
    {EXPECTED_OUTPUT_MULTIPLE_SAMPLES}
    ```

### À retenir

Vous savez comment exploiter le paradigme de flux de données de Nextflow pour paralléliser le traitement par échantillon sur plusieurs échantillons d'entrée.

### Et ensuite ?

Dans la [Partie 3](03_multi_sample.md), vous ajouterez une agrégation multi-échantillons pour combiner les résultats de tous les échantillons.
