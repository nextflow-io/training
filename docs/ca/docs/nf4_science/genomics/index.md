---
title: Genòmica
hide:
  - toc
page_type: index_page
index_type: course
additional_information:
  technical_requirements: true
  learning_objectives:
    - Escriure un workflow lineal per aplicar la detecció de variants a una sola mostra
    - Gestionar fitxers accessoris com ara fitxers d'índex i recursos del genoma de referència de manera apropiada
    - Aprofitar el paradigma de flux de dades de Nextflow per paral·lelitzar la detecció de variants per mostra
    - Implementar la detecció conjunta de múltiples mostres utilitzant operadors de canal rellevants
  audience_prerequisites:
    - "**Audiència:** Aquest curs està dissenyat per a investigadors en genòmica i camps relacionats que volen desenvolupar o personalitzar pipelines d'anàlisi de dades."
    - "**Habilitats:** S'assumeix certa familiaritat amb la línia de comandes, conceptes bàsics de scripting i formats de fitxers genòmics comuns."
    - "**Prerequisits:** Conceptes i eines fonamentals de Nextflow coberts a [Hello Nextflow](../../hello_nextflow/)."
---

# Nextflow per a Genòmica

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Traducció assistida per IA - [més informació i suggeriments](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

**Un curs pràctic que aplica Nextflow a un cas d'ús real de genòmica: detecció de variants amb GATK.**

Aquest curs es basa en la formació per a principiants [Hello Nextflow](../../hello_nextflow/) i demostra com utilitzar Nextflow en el context específic del domini de la genòmica.
Implementareu un pipeline de detecció de variants amb [GATK](https://gatk.broadinstitute.org/) (Genome Analysis Toolkit), un paquet de programari àmpliament utilitzat per analitzar dades de seqüenciació d'alt rendiment.

<!-- additional_information -->

## Visió general del curs

Aquest curs és pràctic, amb exercicis orientats a objectius estructurats per introduir informació gradualment.

Començareu executant les eines de detecció de variants manualment al terminal per entendre la metodologia, i després construireu progressivament un pipeline de Nextflow que automatitza i escala l'anàlisi.

### Pla de lliçons

Hem dividit això en tres parts que se centren cadascuna en aspectes específics de l'aplicació de Nextflow a un cas d'ús de genòmica.

| Capítol del curs                                                         | Resum                                                                                                                  | Durada estimada |
| ------------------------------------------------------------------------ | ---------------------------------------------------------------------------------------------------------------------- | --------------- |
| [Part 1: Visió general del mètode](./01_method.md)                      | Comprendre la metodologia de detecció de variants i executar les eines manualment                                      | 30 min          |
| [Part 2: Detecció de variants per mostra](./02_per_sample_variant_calling.md) | Construir un pipeline que indexa fitxers BAM i detecta variants, i després escalar a múltiples mostres                | 60 min          |
| [Part 3: Detecció conjunta en una cohort](./03_joint_calling.md)        | Afegir genotipat conjunt de múltiples mostres utilitzant operadors de canal per agregar sortides per mostra           | 45 min          |

Al final d'aquest curs, podreu aplicar conceptes i eines fonamentals de Nextflow a un cas d'ús típic de genòmica.

Preparat per fer el curs?

[Comença :material-arrow-right:](00_orientation.md){ .md-button .md-button--primary }

<!-- Clearfix for float -->
<div style="content: ''; clear: both; display: table;"></div>
