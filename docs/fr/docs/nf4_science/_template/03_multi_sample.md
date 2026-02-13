# Partie 3 : Agrégation multi-échantillons

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Traduction assistée par IA - [en savoir plus et suggérer des améliorations](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

Dans la Partie 2, vous avez construit un pipeline de traitement par échantillon qui gérait chaque échantillon de manière indépendante.
Nous allons maintenant l'étendre pour implémenter {AGGREGATION_METHOD} multi-échantillons, comme présenté dans la [Partie 1](01_method.md).

## Objectif

Dans cette partie du cours, nous allons étendre le workflow pour effectuer les opérations suivantes :

<figure class="excalidraw">
--8<-- "docs/en/docs/nf4_science/{DOMAIN_DIR}/img/{DIAGRAM_2}.svg"
</figure>

1. {PER_SAMPLE_STEP_1}
2. {PER_SAMPLE_STEP_2}
3. {AGGREGATION_STEP_1}
4. {AGGREGATION_STEP_2}

Cette partie s'appuie directement sur le workflow produit dans la Partie 2.

??? info "Comment commencer à partir de cette section"

    Cette section du cours suppose que vous avez terminé la [Partie 2 : Traitement d'un échantillon unique](./02_single_sample.md) et que vous disposez d'un pipeline `{DOMAIN_DIR}.nf` fonctionnel.

    Si vous n'avez pas terminé la Partie 2 ou souhaitez repartir de zéro pour cette partie, vous pouvez utiliser la solution de la Partie 2 comme point de départ.
    Exécutez ces commandes depuis le répertoire `nf4-science/{DOMAIN_DIR}/` :

    ```bash
    cp solutions/part2/{DOMAIN_DIR}-2.nf {DOMAIN_DIR}.nf
    cp solutions/part2/nextflow.config .
    cp solutions/part2/modules/* modules/
    ```

    Cela vous donne un workflow complet de traitement d'échantillon unique.
    Vous pouvez vérifier qu'il s'exécute correctement en lançant la commande suivante :

    ```bash
    nextflow run {DOMAIN_DIR}.nf -profile test
    ```

## Plan de la leçon

Nous avons divisé cette partie en deux étapes :

1. **{MODIFICATION_STEP_SUMMARY}.**
   Cela couvre la mise à jour des commandes et sorties des processus.
2. **{AGGREGATION_STEP_SUMMARY}.**
   Cela introduit l'opérateur `collect()` {AND_OTHER_CONCEPTS}.

!!! note "Note"

     Assurez-vous d'être dans le bon répertoire de travail :
     `cd /workspaces/training/nf4-science/{DOMAIN_DIR}`

---

## 1. {MODIFICATION_STEP_TITLE}

{DESCRIPTION_OF_MODIFICATION_TO_EXISTING_PROCESSES}

Rappelez-vous la commande modifiée de la [Partie 1](01_method.md) :

```bash
{MODIFIED_COMMAND}
```

{EXPLAIN_DIFFERENCES_FROM_PART_2}

### 1.1. {MODIFICATION_SUBSTEP}

{INSTRUCTIONS_WITH_BEFORE_AFTER_TABS}

=== "Après"

    ```groovy title="modules/{MODULE_FILE}.nf" linenums="{N}" hl_lines="{LINES}"
    {UPDATED_CODE}
    ```

=== "Avant"

    ```groovy title="modules/{MODULE_FILE}.nf" linenums="{N}" hl_lines="{LINES}"
    {ORIGINAL_CODE}
    ```

### 1.2. Exécuter le workflow pour vérifier la modification

```bash
nextflow run {DOMAIN_DIR}.nf -profile test
```

??? success "Sortie de la commande"

    ```console
    {EXPECTED_OUTPUT}
    ```

{VERIFY_OUTPUT_INSTRUCTIONS}

### À retenir

Vous savez comment modifier les commandes et sorties des processus pour adapter le comportement du workflow.

### Et ensuite ?

Ajoutez l'étape d'agrégation multi-échantillons.

---

## 2. {AGGREGATION_STEP_TITLE}

{DESCRIPTION_OF_AGGREGATION}

### 2.1. Écrire le module d'agrégation

{MODULE_INSTRUCTIONS}

### 2.2. Collecter les sorties par échantillon et les transmettre au processus d'agrégation

{INSTRUCTIONS_USING_COLLECT_OPERATOR}

### 2.3. Exécuter le workflow complet

```bash
nextflow run {DOMAIN_DIR}.nf -profile test
```

??? success "Sortie de la commande"

    ```console
    {EXPECTED_OUTPUT}
    ```

{VERIFY_OUTPUT_INSTRUCTIONS}

---

### À retenir

Vous disposez d'un pipeline complet qui traite les échantillons individuellement et agrège les résultats de tous les échantillons.
Vous savez comment utiliser des opérateurs de canaux comme `collect()` pour agréger les sorties par échantillon pour une analyse multi-échantillons.

### Et ensuite ?

Félicitations pour avoir terminé ce cours ! Rendez-vous au [résumé du cours](next_steps.md) pour revoir ce que vous avez appris et explorer les prochaines étapes.
