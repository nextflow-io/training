---
title: RNAseq
hide:
  - toc
page_type: index_page
index_type: course
additional_information:
  technical_requirements: true
  learning_objectives:
    - Escriure un workflow lineal per aplicar mètodes bàsics de processament i control de qualitat d'RNAseq
    - Gestionar fitxers específics del domini com FASTQ i recursos de genoma de referència de manera apropiada
    - Gestionar dades de seqüenciació single-end i paired-end
    - Aprofitar el paradigma de flux de dades de Nextflow per paral·lelitzar el processament d'RNAseq per mostra
    - Agregar informes de control de qualitat a través de múltiples passos i mostres utilitzant operadors de canal rellevants
  audience_prerequisites:
    - "**Audiència:** Aquest curs està dissenyat per a investigadors en transcriptòmica i camps relacionats que volen desenvolupar o personalitzar pipelines d'anàlisi de dades."
    - "**Habilitats:** S'assumeix certa familiaritat amb la línia de comandes, conceptes bàsics de scripting i formats comuns de fitxers d'RNAseq."
    - "**Prerequisits:** Conceptes fonamentals de Nextflow i eines cobertes a [Hello Nextflow](../../hello_nextflow/)."
---

# Nextflow per a RNAseq

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Traducció assistida per IA - [més informació i suggeriments](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

**Un curs pràctic que aplica Nextflow a un cas d'ús real de transcriptòmica: processament d'RNAseq bulk amb Trim Galore, HISAT2 i FastQC.**

Aquest curs es basa en la formació per a principiants [Hello Nextflow](../../hello_nextflow/) i demostra com utilitzar Nextflow en el context específic de l'anàlisi d'RNAseq bulk.
Implementareu un pipeline de processament que retalla seqüències adaptadores, alinea lectures a un genoma de referència i realitza control de qualitat (QC) en diverses etapes.

<!-- additional_information -->

## Visió general del curs

Aquest curs és pràctic, amb exercicis orientats a objectius estructurats per introduir informació gradualment.

Començareu executant les eines de processament manualment al terminal per entendre la metodologia, i després construireu progressivament un pipeline de Nextflow que automatitza i escala l'anàlisi.

### Pla de lliçons

Hem dividit això en tres parts que se centren cadascuna en aspectes específics de l'aplicació de Nextflow a un cas d'ús d'RNAseq.

| Capítol del curs                                                      | Resum                                                                                                                  | Durada estimada |
| --------------------------------------------------------------------- | ---------------------------------------------------------------------------------------------------------------------- | --------------- |
| [Part 1: Visió general del mètode](./01_method.md)                    | Comprendre la metodologia de processament d'RNAseq i executar les eines manualment                                     | 30 min          |
| [Part 2: Implementació d'una sola mostra](./02_single-sample.md)      | Construir un pipeline que retalla, alinea i fa QC d'una sola mostra, i després escalar per gestionar múltiples mostres | 60 min          |
| [Part 3: Implementació multi-mostra paired-end](./03_multi-sample.md) | Estendre el pipeline per gestionar dades paired-end i agregar informes de control de qualitat a través de les mostres  | 45 min          |

Al final d'aquest curs, podreu aplicar conceptes fonamentals de Nextflow i eines a un cas d'ús típic d'RNAseq.

Preparat per fer el curs?

[Comença :material-arrow-right:](00_orientation.md){ .md-button .md-button--primary }

<!-- Clearfix for float -->
<div style="content: ''; clear: both; display: table;"></div>
