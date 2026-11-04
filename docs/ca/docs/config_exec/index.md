---
title: Configuració d'Execució
hide:
  - toc
page_type: index_page
index_type: course
additional_information:
  technical_requirements: true
  learning_objectives:
    - Canviar la tecnologia d'empaquetament de programari entre Docker i Conda
    - Seleccionar una plataforma d'execució i entendre com Nextflow adapta l'execució de tasques
    - Controlar l'assignació de recursos de còmput i reintentar automàticament les tasques que fallen
    - Definir i combinar perfils per canviar entre configuracions predefinides
  audience_prerequisites:
    - "**Audiència:** Aquest curs està dissenyat per a persones que ja saben com llançar pipelines de Nextflow localment i volen configurar l'execució amb més profunditat."
    - "**Habilitats:** S'assumeix una certa familiaritat amb la línia de comandes."
    - "**Cursos:** Cal haver completat [Nextflow Run](../nextflow_run/index.md) o estar còmode executant un pipeline local amb `nextflow run`."
---

# Configuració d'Execució

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Traducció assistida per IA - [més informació i suggeriments](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

**Configure Execution és una introducció pràctica a l'adaptació de l'execució de pipelines de Nextflow a diferents entorns de còmput.**

Mitjançant exercicis orientats a objectius, aprendreu a canviar la tecnologia d'empaquetament de programari, seleccionar una plataforma d'execució, controlar l'assignació de recursos de còmput i els reintents, i agrupar la configuració en perfils commutables.

Adquirireu les habilitats i la confiança per configurar l'execució de pipelines de Nextflow com un professional.

<!-- additional_information -->

## Visió general del curs

Aquest curs és pràctic i es basa en les habilitats cobertes a [Nextflow Run](../nextflow_run/index.md).

Partireu del mateix pipeline de múltiples passos d'aquell curs i adaptareu progressivament la seva configuració a diferents entorns de còmput, per després agrupar-ho tot en perfils entre els quals podreu canviar en temps d'execució.

### Pla de lliçons

| Capítol del curs                                                                 | Resum                                                                                       | Durada estimada |
| -------------------------------------------------------------------------------- | ------------------------------------------------------------------------------------------- | --------------- |
| [Part 1: Adaptar-se al vostre entorn de còmput](./01_packaging_and_execution.md) | Canviar la tecnologia d'empaquetament de programari i seleccionar una plataforma d'execució | 20 min          |
| [Part 2: Gestionar recursos de còmput i fallades](./02_resources_and_retries.md) | Controlar l'assignació de recursos i reintentar automàticament les tasques que fallen       | 15 min          |
| [Part 3: Utilitzar perfils per canviar de configuració](./03_profiles.md)        | Definir i combinar perfils, i inspeccionar la configuració completament resolta             | 15 min          |

Al final d'aquest curs, estareu còmodes configurant pipelines de Nextflow per a una varietat d'entorns de còmput i canviant entre ells amb el mínim d'inconvenients.

Esteu preparats per fer el curs?

[Comenceu a aprendre :material-arrow-right:](00_orientation.md){ .md-button .md-button--primary }
