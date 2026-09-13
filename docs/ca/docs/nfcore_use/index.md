---
title: Utilitza nf-core
hide:
  - toc
page_type: index_page
index_type: course
additional_information:
  technical_requirements: true
  learning_objectives:
    - Trobar, obtenir i executar pipelines de la comunitat nf-core
    - Configurar l'execució de pipelines mitjançant paràmetres i fitxers de configuració
    - Entendre com els pipelines nf-core validen els paràmetres i les dades d'entrada
    - Executar un pipeline a escala de producció (nf-core/rnaseq) i sobreescriure les seves assignacions de recursos per defecte
  audience_prerequisites:
    - "**Audiència:** Aquest curs està dissenyat per a persones que ja saben executar pipelines locals de Nextflow i que són noves a nf-core, i volen executar pipelines de la comunitat existents."
    - "**Habilitats:** S'assumeix una certa familiaritat amb la línia de comandes, conceptes bàsics de scripting i formats de fitxer comuns."
    - "**Cursos:** Cal haver completat [Nextflow Run](../nextflow_run/index.md) o estar còmode executant un pipeline local amb `nextflow run`."
    - "**Àmbit:** Els exercicis utilitzen pipelines de bioinformàtica, però no es requereix cap coneixement científic previ del domini."
---

# Utilitza nf-core

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Traducció assistida per IA - [més informació i suggeriments](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

**Utilitza nf-core és una introducció pràctica per trobar, executar i configurar pipelines de la comunitat nf-core.**

Treballant amb exemples pràctics i exercicis guiats, aprendreu a trobar i obtenir pipelines nf-core, executar-los utilitzant els seus perfils de prova integrats, i personalitzar-ne l'execució mitjançant paràmetres i fitxers de configuració.

Adquirireu les habilitats i la confiança necessàries per començar a executar pipelines nf-core per a les vostres pròpies anàlisis.

<!-- additional_information -->

## Visió general del curs

Aquest curs és pràctic, amb exercicis orientats a objectius estructurats per introduir la informació de manera gradual.

Començareu amb `nf-core/demo`, un pipeline mínim mantingut pel projecte nf-core amb finalitats de formació, i després aplicareu el que heu après a `nf-core/rnaseq`, un pipeline de producció àmpliament utilitzat per a l'anàlisi de seqüenciació d'RNA en bloc.

Aquest curs se centra en l'execució de pipelines.
Si busqueu una introducció al desenvolupament de pipelines compatibles amb nf-core, consulteu [Build with nf-core](../nfcore_build/index.md).

### Pla de lliçons

| Capítol del curs                                                     | Resum                                                                                                          | Durada estimada |
| -------------------------------------------------------------------- | -------------------------------------------------------------------------------------------------------------- | --------------- |
| [Part 1: Executar un pipeline de demostració](./01_run_demo.md)      | Trobar i obtenir un pipeline nf-core i executar-lo utilitzant el seu perfil de prova                           | 20 min          |
| [Part 2: Configurar l'execució del pipeline](./02_configure_execution.md) | Establir paràmetres, entendre la validació i personalitzar l'assignació de recursos i els arguments d'eines | 20 min          |
| [Part 3: Executar un pipeline de producció](./03_run_production_pipeline.md) | Obtenir i executar nf-core/rnaseq, i sobreescriure les seves assignacions de recursos per defecte        | 20 min          |

Al final d'aquest curs, podreu aprofitar la gran quantitat de pipelines de la comunitat que ofereix el projecte nf-core.

Esteu preparats per fer el curs?

[Comenceu a aprendre :material-arrow-right:](00_orientation.md){ .md-button .md-button--primary }
