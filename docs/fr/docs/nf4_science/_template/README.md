# Modèle de cours nf4_science

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Traduction assistée par IA - [en savoir plus et suggérer des améliorations](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

Ceci est un modèle générique pour créer un nouveau cours spécifique à un domaine sous `nf4_science/`.
Il est basé sur les modèles communs trouvés dans les cours Genomics et RNAseq.

## Comment l'utiliser

1. Copiez ce répertoire et le répertoire de scripts correspondant pour créer votre cours :

   ```bash
   # Documentation
   cp -r docs/en/docs/nf4_science/_template docs/en/docs/nf4_science/{your_domain}

   # Scripts
   cp -r nf4-science/_template nf4-science/{your_domain}
   ```

2. Renommez le script de workflow :

   ```bash
   mv nf4-science/{your_domain}/_template.nf nf4-science/{your_domain}/{your_domain}.nf
   ```

3. Renommez les fichiers de modules pour correspondre à vos outils :

   ```bash
   mv nf4-science/{your_domain}/modules/_tool_a.nf nf4-science/{your_domain}/modules/{tool_a}.nf
   mv nf4-science/{your_domain}/modules/_tool_b.nf nf4-science/{your_domain}/modules/{tool_b}.nf
   ```

4. Recherchez toutes les valeurs `{PLACEHOLDER}` dans la documentation et les scripts et remplacez-les par votre contenu.

5. Ajoutez à la navigation de `mkdocs.yml` (voir l'extrait ci-dessous).

6. Remplissez le répertoire `data/` avec des jeux de données de test.

7. Développez le répertoire `solutions/` avec du code fonctionnel pour chaque partie.

## Extrait de navigation mkdocs.yml

```yaml
- { Your Domain }:
    - nf4_science/{your_domain}/index.md
    - nf4_science/{your_domain}/00_orientation.md
    - nf4_science/{your_domain}/01_method.md
    - nf4_science/{your_domain}/02_single_sample.md
    - nf4_science/{your_domain}/03_multi_sample.md
    - nf4_science/{your_domain}/survey.md
    - nf4_science/{your_domain}/next_steps.md
```

## Référence des placeholders

### Placeholders au niveau du cours

| Placeholder                  | Description                            | Exemple (Genomics)                                    |
| ---------------------------- | -------------------------------------- | ----------------------------------------------------- |
| `{DOMAIN}`                   | Nom du domaine (casse titre)           | Genomics                                              |
| `{DOMAIN_DIR}`               | Nom du répertoire (minuscules)         | genomics                                              |
| `{METHOD}`                   | Nom de la méthode d'analyse            | variant calling                                       |
| `{METHOD_SHORT_DESCRIPTION}` | Description de la méthode en une ligne | variant calling with GATK                             |
| `{ACCESSORY_FILES}`          | Types de fichiers accessoires utilisés | fichiers d'index et ressources du génome de référence |

### Placeholders des outils

| Placeholder                           | Description                         | Exemple (Genomics)                                  |
| ------------------------------------- | ----------------------------------- | --------------------------------------------------- |
| `{TOOL_A}` / `{TOOL_B}`               | Noms d'affichage des outils         | Samtools / GATK                                     |
| `{TOOL_A_URL}` / `{TOOL_B_URL}`       | URLs de documentation des outils    | https://www.htslib.org/                             |
| `{TOOL_A_CONTAINER_URI}`              | URI complète du conteneur           | community.wave.seqera.io/library/samtools:1.20--... |
| `{TOOL_A_MODULE}` / `{TOOL_B_MODULE}` | Noms de fichiers modules (sans .nf) | samtools_index / gatk_haplotypecaller               |
| `{TOOL_A_PROCESS_NAME}`               | Nom du processus en MAJUSCULES      | SAMTOOLS_INDEX                                      |
| `{TOOL_A_COMMAND}`                    | La commande shell à exécuter        | samtools index '<input_bam>'                        |

### Placeholders d'entrée/sortie

| Placeholder            | Description                             | Exemple (Genomics)   |
| ---------------------- | --------------------------------------- | -------------------- |
| `{PRIMARY_INPUT_TYPE}` | Type de fichier d'entrée principal      | fichier BAM          |
| `{PRIMARY_PARAM_NAME}` | Nom du paramètre Nextflow               | reads_bam            |
| `{TEST_INPUT_PATH}`    | Chemin vers l'entrée de test dans data/ | bam/reads_mother.bam |
| `{FINAL_OUTPUT_TYPE}`  | Type de sortie finale du pipeline       | VCFs joint-called    |

### Placeholders de contenu

| Placeholder             | Description                                      |
| ----------------------- | ------------------------------------------------ |
| `{PART2_SUMMARY}`       | Résumé en une ligne de la portée de la Partie 2  |
| `{AGGREGATION_SUMMARY}` | Résumé en une ligne de l'étape d'agrégation      |
| `{PROCESS_LIST}`        | Liste des noms de processus séparés par virgules |
| `{TYPEFORM_ID}`         | ID d'intégration Typeform pour le questionnaire  |

## Structure pédagogique

Le modèle suit un arc en trois parties :

1. **Partie 1 (01_method.md)** : Tests manuels dans des conteneurs Docker pour comprendre la méthodologie
2. **Partie 2 (02_single_sample.md)** : Encapsuler les commandes dans Nextflow ; échantillon unique, puis batch
3. **Partie 3 (03_multi_sample.md)** : Ajouter l'agrégation multi-échantillons en utilisant les opérateurs de canaux

### Conventions clés

- Utiliser des **blocs de code à onglets Avant/Après** avec `hl_lines` pour montrer les modifications de code
- Utiliser `??? success "Command output"` pour la sortie attendue repliable
- Utiliser `??? abstract "Directory contents"` pour les arborescences de répertoires repliables
- Terminer chaque section majeure par les sous-sections **À retenir** et **Et ensuite ?**
- Utiliser des règles horizontales `---` pour séparer les sections numérotées majeures
- Fournir des fichiers squelettes (workflow + modules) à compléter par les apprenant·es
- Organiser les solutions par partie (`solutions/part2/`, `solutions/part3/`)

## Structure des répertoires

```
docs/en/docs/nf4_science/{domain}/
├── index.md                    # Course overview with frontmatter
├── 00_orientation.md           # Environment setup
├── 01_method.md                # Manual testing in containers
├── 02_single_sample.md         # Single-sample Nextflow implementation
├── 03_multi_sample.md          # Multi-sample aggregation
├── survey.md                   # Typeform feedback survey
├── next_steps.md               # Course summary and suggestions
└── img/                        # Diagrams (.excalidraw.svg, .png)

nf4-science/{domain}/
├── {domain}.nf                 # Skeleton workflow file
├── nextflow.config             # Minimal config (docker.enabled = true)
├── data/                       # Test datasets and resources
│   └── samplesheet.csv         # Sample metadata
├── modules/                    # Skeleton module files
│   ├── {tool_a}.nf
│   └── {tool_b}.nf
└── solutions/
    ├── part2/                  # Complete Part 2 solution
    │   ├── {domain}-2.nf
    │   ├── nextflow.config
    │   └── modules/
    └── part3/                  # Complete Part 3 solution
        ├── {domain}-3.nf
        ├── nextflow.config
        └── modules/
```
