---
title: Hello nf-core
hide:
  - toc
page_type: index_page
index_type: course
additional_information:
  technical_requirements: true
  learning_objectives:
    - Recuperar, executar i gestionar l'execució de pipelines nf-core
    - Descriure l'estructura del codi i l'organització del projecte de pipelines nf-core
    - Crear un pipeline bàsic compatible amb nf-core a partir d'una plantilla
    - Actualitzar un workflow de Nextflow senzill per complir els estàndards nf-core
    - Afegir mòduls nf-core a un pipeline compatible amb nf-core
    - Contribuir els teus propis mòduls a nf-core
    - Validar entrades i paràmetres utilitzant les eines nf-core
  audience_prerequisites:
    - "**Audiència:** Aquest curs està dissenyat per a estudiants que ja estan familiaritzats amb Nextflow bàsic i volen aprendre a utilitzar recursos i bones pràctiques d'nf-core."
    - "**Habilitats:** S'assumeix familiaritat amb la línia de comandes, conceptes bàsics de scripting i formats de fitxer comuns."
    - "**Cursos:** Cal haver completat el curs [Hello Nextflow](../hello_nextflow/index.md) o equivalent."
    - "**Àmbit:** Tots els exercicis són independents del domini, per tant no es requereix coneixement científic previ."
---

# Hello nf-core

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Traducció assistida per IA - [més informació i suggeriments](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

**Hello nf-core és una introducció pràctica a l'ús de recursos i bones pràctiques d'nf-core.**

![nf-core logo](./img/nf-core-logo.png#only-light)
![nf-core logo](./img/nf-core-logo-darkbg.png#only-dark)

Treballant amb exemples pràctics i exercicis guiats, aprendràs a utilitzar i desenvolupar mòduls i pipelines compatibles amb nf-core, i a utilitzar les eines nf-core de manera efectiva.

Adquiriràs les habilitats i la confiança per començar a desenvolupar pipelines seguint les bones pràctiques d'nf-core.

<!-- additional_information -->

## Visió general del curs

Aquest curs està dissenyat per ser pràctic, amb exercicis orientats a objectius estructurats per introduir informació gradualment.

Se t'introduirà a [**nf-core**](https://nf-co.re/), un esforç comunitari per desenvolupar i mantenir un conjunt seleccionat de pipelines científics construïts amb Nextflow, així com eines i directrius rellevants que promouen el desenvolupament obert, les proves i la revisió per parells ([Nat Biotechnol 38, 276–278 (2020)](https://www.nature.com/articles/s41587-020-0439-x), [Genome Biol 26, 228 (2025)](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-025-03673-9)).

Els pipelines desenvolupats per la comunitat nf-core estan dissenyats per ser modulars, escalables i portables, permetent als investigadors adaptar-los i executar-los fàcilment utilitzant les seves pròpies dades i recursos de càlcul.
Les directrius de bones pràctiques aplicades pel projecte asseguren encara més que els pipelines siguin robustos, ben documentats i validats amb conjunts de dades reals.
Això ajuda a augmentar la fiabilitat i la reproducibilitat de les anàlisis científiques i, en última instància, permet als investigadors accelerar els seus descobriments científics.

No cobrirem tot el que hi ha per saber sobre els pipelines nf-core en aquest curs, perquè nf-core engloba moltes característiques i convencions desenvolupades per la comunitat durant anys.
En canvi, ens centrarem en els conceptes essencials que t'ajudaran a començar i entendre com funciona nf-core.

### Pla de lliçons

Hem dividit això en cinc parts que se centraran cadascuna en aspectes específics de l'ús de recursos nf-core.

| Capítol del curs                                                | Resum                                                                                                                                                                   | Durada estimada |
| --------------------------------------------------------------- | ----------------------------------------------------------------------------------------------------------------------------------------------------------------------- | --------------- |
| [Part 1: Executar un pipeline de demostració](./01_run_demo.md) | Executar un pipeline nf-core existent i examinar la seva estructura de codi per tenir una idea del que fa aquests pipelines diferents dels workflows bàsics de Nextflow | 30 min          |
| [Part 2: Reescriure Hello per a nf-core](./02_rewrite_hello.md) | Adaptar un workflow existent a l'estructura de la plantilla nf-core, començant pel workflow senzill produït al curs [Hello Nextflow](../hello_nextflow/index.md)        | 60 min          |
| [Part 3: Utilitzar un mòdul nf-core](./03_use_module.md)        | Explorar la biblioteca de mòduls de la comunitat i aprendre a integrar mòduls preconstruïts i provats que encapsulen eines bioinformàtiques comunes                     | 30 min          |
| [Part 4: Crear un mòdul nf-core](./04_make_module.md)           | Crear el teu propi mòdul a l'estil nf-core utilitzant l'estructura específica, convencions de nomenclatura i requisits de metadades establerts per nf-core              | 30 min          |
| [Part 5: Afegir validació d'entrada](./05_input_validation.md)  | Implementar validació d'entrada tant per a paràmetres de línia de comandes com per a fitxers de dades d'entrada utilitzant nf-schema                                    | 30 min          |

Al final d'aquest curs, podràs aprofitar l'enorme riquesa de recursos que ofereix el projecte nf-core.

Preparat per fer el curs?

[Comença a aprendre :material-arrow-right:](00_orientation.md){ .md-button .md-button--primary }
