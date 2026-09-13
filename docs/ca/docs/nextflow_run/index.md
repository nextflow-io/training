---
title: Nextflow Run
hide:
  - toc
page_type: index_page
index_type: course
additional_information:
  technical_requirements: true
  learning_objectives:
    - Llançar i gestionar pipelines de Nextflow des de la línia de comandes
    - Entendre com els canals i els operadors permeten workflows eficients amb múltiples entrades i múltiples passos
    - Utilitzar contenidors per gestionar les dependències de programari i garantir la reproducibilitat
    - Configurar l'execució de pipelines i les sortides
    - Generar informes d'execució, inspeccionar l'historial d'execucions anteriors i netejar directoris de treball antics
    - Executar pipelines directament des de repositoris remots com GitHub
  audience_prerequisites:
    - "**Audiència:** Aquest curs està dissenyat per a estudiants que són completament nous a Nextflow i volen executar pipelines existents."
    - "**Habilitats:** S'assumeix certa familiaritat amb la línia de comandes, conceptes bàsics de scripting i formats de fitxer comuns."
    - "**Àmbit:** Tots els exercicis són independents del domini, per tant no es requereix coneixement científic previ."
---

# Nextflow Run

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Traducció assistida per IA - [més informació i suggeriments](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

**Nextflow Run és una introducció pràctica a l'execució de workflows d'anàlisi de dades reproduïbles i escalables.**

Treballant amb una sèrie d'exercicis orientats a objectius, aprendràs els fonaments per llançar i gestionar pipelines de Nextflow, entendràs com els canals i els operadors permeten el processament paral·lel de múltiples entrades, i utilitzaràs contenidors per gestionar les dependències de programari.

Obtindràs les habilitats i la confiança per començar a executar workflows amb Nextflow.

<!-- additional_information -->

## Visió general del curs

Aquest curs és pràctic, amb exercicis orientats a objectius estructurats per introduir informació gradualment.

Executaràs diverses versions d'un pipeline de Nextflow que processa entrades de text, començant amb una versió senzilla d'un sol pas i progressant fins a una versió de múltiples passos que pren un fitxer CSV d'entrades, executa uns quants passos de transformació i genera un únic fitxer de text que conté art ASCII generat per una eina en contenidor.

Aquest curs se centra en l'execució de pipelines (anomenat així per la comanda bàsica `nextflow run`).
Si busques una introducció al desenvolupament de pipelines de Nextflow, consulta [Hello Nextflow](../hello_nextflow/index.md).

!!! note "Nota"

    Busques la versió anterior d'aquest curs? Ha estat substituïda per la versió d'aquesta pàgina, però encara es pot consultar a la [versió 3.6.1](https://training.nextflow.io/3.6.1/nextflow_run/) del lloc de formació.

### Pla de lliçons

| Capítol del curs                                                           | Resum                                                                                                               | Durada estimada |
| -------------------------------------------------------------------------- | ------------------------------------------------------------------------------------------------------------------- | --------------- |
| [Part 1: Executar Nextflow](./01_run_nextflow.md)                          | Llançar i gestionar pipelines de Nextflow, i entendre la mecànica essencial dels workflows                          | 25 min          |
| [Part 2: Configurar el pipeline](./02_configure_pipeline.md)               | Configurar l'execució del pipeline i les sortides mitjançant `nextflow.config`                                      | 20 min          |
| [Part 3: Gestionar les execucions del workflow](./03_manage_executions.md) | Generar informes d'execució, inspeccionar l'historial d'execucions anteriors i netejar directoris de treball antics | 10 min          |
| [Part 4: Executar pipelines remots](./04_remote_repositories.md)           | Executar un pipeline directament des de GitHub i fixar-lo a una revisió específica                                  | 10 min          |

Al final d'aquest curs, estaràs ben preparat/da per abordar els següents passos en el teu viatge per executar workflows reproduïbles per a les teves necessitats de computació científica.

Preparat/da per fer el curs?

[Comença a aprendre :material-arrow-right:](00_orientation.md){ .md-button .md-button--primary }
