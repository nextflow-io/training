---
title: Nextflow Run
hide:
  - toc
page_type: index_page
index_type: course
additional_information:
  technical_requirements: true
  learning_objectives:
    - Llançar i gestionar l'execució de workflows de Nextflow
    - Trobar i interpretar sortides (resultats) i fitxers de registre
    - Reconèixer components bàsics de Nextflow en un workflow senzill de múltiples passos
    - Configurar l'execució de pipelines per executar-se en plataformes de computació comunes, incloent HPC i núvol
    - Resumir les millors pràctiques per a la reproducibilitat, portabilitat i reutilització de codi que fan que els pipelines siguin FAIR, incloent la modularitat del codi i els contenidors de programari
  audience_prerequisites:
    - "**Audiència:** Aquest curs està dissenyat per a estudiants que són completament nous a Nextflow i volen executar pipelines existents."
    - "**Habilitats:** S'assumeix certa familiaritat amb la línia de comandes, conceptes bàsics de scripting i formats de fitxer comuns."
    - "**Àmbit:** Tots els exercicis són independents del domini, per tant no es requereix coneixement científic previ."
---

# Nextflow Run

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Traducció assistida per IA - [més informació i suggeriments](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

**Nextflow Run és una introducció pràctica a l'execució de workflows d'anàlisi de dades reproduïbles i escalables.**

Treballant amb exemples pràctics i exercicis guiats, aprendràs els fonaments de l'ús de Nextflow, incloent com executar pipelines, gestionar fitxers i dependències de programari, paral·lelitzar l'execució sense esforç, i executar workflows en diferents entorns de computació.

Obtindràs les habilitats i la confiança per començar a executar workflows amb Nextflow.

<!-- additional_information -->

## Visió general del curs

### Què faràs

Aquest curs és pràctic, amb exercicis orientats a objectius estructurats per introduir informació gradualment.

Executaràs diverses versions d'un pipeline de Nextflow que processa entrades de text.
Començaràs amb una versió senzilla que consisteix en un sol pas, i eventualment progressaràs a una versió de múltiples passos que pren un fitxer CSV d'entrades de text tabulars, executa uns quants passos de transformació, i genera un únic fitxer de text que conté una imatge ASCII d'un personatge dient el text transformat.

Aquest curs se centra en l'execució de pipelines (anomenat així per la comanda bàsica `nextflow run`).
Si busques una introducció al desenvolupament de pipelines de Nextflow, consulta [Hello Nextflow](../hello_nextflow/index.md).

### Pla de lliçons

Hem dividit això en tres parts que se centraran cadascuna en aspectes específics de l'execució i gestió de pipelines escrits en Nextflow.

| Capítol del curs                                      | Resum                                                                                                                                       | Durada estimada |
| ----------------------------------------------------- | ------------------------------------------------------------------------------------------------------------------------------------------- | --------------- |
| [Part 1: Operacions bàsiques d'execució](./01_basics.md) | Llançar i gestionar l'execució d'un workflow senzill                                                                                        | 30 mins         |
| [Part 2: Executar pipelines reals](./02_pipeline.md)     | Processar entrades complexes, executar workflows de múltiples passos, utilitzar contenidors i paral·lelitzar l'execució sense esforç       | 60 mins         |
| [Part 3: Configuració d'execució](./03_config.md)        | Personalitzar el comportament del pipeline i optimitzar l'ús en diferents entorns computacionals                                           | 60 mins         |

Al final d'aquest curs, estaràs ben preparat/da per abordar els següents passos en el teu viatge per executar workflows reproduïbles per a les teves necessitats de computació científica.

Preparat/da per fer el curs?

[Comença a aprendre :material-arrow-right:](00_orientation.md){ .md-button .md-button--primary }
