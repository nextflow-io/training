---
title: Escala amb Seqera
hide:
  - toc
page_type: index_page
index_type: course
additional_information:
  technical_requirements: true
  learning_objectives:
    - Registrar-se a Seqera Platform i explorar el Community Showcase
    - Afegir un pipeline a un workspace i llançar-lo des de la interfície web
    - Autenticar-se i llançar pipelines des de la línia de comandes amb el `tw` CLI
    - Registrar un pipeline allotjat a GitHub i llançar-lo de les dues maneres
  audience_prerequisites:
    - "**Audiència:** Aquest curs està dissenyat per a persones que volen executar pipelines de Nextflow a escala utilitzant Seqera Platform."
    - "**Habilitats:** Es pressuposa familiaritat amb l'execució de pipelines nf-core des de la línia de comandes."
    - "**Cursos:** Cal haver completat [Nextflow Run](../nextflow_run/index.md) i [Use nf-core](../nfcore_use/index.md), o bé tenir experiència executant pipelines locals i `nf-core/rnaseq`."
---

# Escala amb Seqera

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Traducció assistida per IA - [més informació i suggeriments](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

**Scale with Seqera és una introducció pràctica per llançar i monitoritzar pipelines de Nextflow amb Seqera Platform.**

Treballant amb exemples pràctics, configurareu l'accés a Seqera Platform, llançareu un pipeline a escala de producció tant des de la interfície web com des de la línia de comandes, i afegireu un nou pipeline al vostre workspace.

Adquirireu les habilitats i la confiança necessàries per executar i monitoritzar els vostres propis pipelines a Seqera Platform.

<!-- additional_information -->

## Visió general del curs

Aquest curs és pràctic i es basa en els pipelines que ja heu executat a [Use nf-core](../nfcore_use/index.md).

Començareu registrant-vos a Seqera Platform i llançant `nf-core/rnaseq`, un pipeline a escala de producció, des de la interfície web.
Després passareu a l'eina de línia de comandes `tw` per fer el mateix des d'un terminal, i finalment registrareu un nou pipeline, `nf-core/demo`, i el llançareu de les dues maneres.

### Pla de lliçons

| Capítol del curs                                                                 | Resum                                                                                                     | Durada estimada |
| -------------------------------------------------------------------------------- | --------------------------------------------------------------------------------------------------------- | --------------- |
| [Part 1: Llançar pipelines des de la interfície web](./01_run_with_seqera.md)    | Configurar l'accés a Seqera Platform i llançar un pipeline a escala de producció des de la interfície web | 20 min          |
| [Part 2: Llançar pipelines des de la línia de comandes](./02_launch_from_cli.md) | Autenticar el `tw` CLI, llançar un pipeline desat i registrar un nou pipeline des del CLI                 | 25 min          |

Al final d'aquest curs, us sentireu còmodes llançant i monitoritzant pipelines de Nextflow a Seqera Platform, tant si preferiu treballar des de la interfície web com des de la línia de comandes.

Esteu preparats per començar el curs?

[Comenceu a aprendre :material-arrow-right:](00_orientation.md){ .md-button .md-button--primary }
