---
title: Inici
description: Benvingut/da al portal de formació de la comunitat Nextflow!
hide:
  - toc
  - footer
---

# Formació en Nextflow

<div class="grid cards" markdown>

-   :material-book-open-variant:{ .lg .middle } __Cursos d'autoservei__

    ---

    **Benvingut/da al portal de formació de la comunitat Nextflow!**

    Seguiu els cursos següents al vostre ritme, en el nostre entorn web o en el vostre propi.
    Cada curs és pràctic, amb exercicis orientats a objectius que podeu completar de manera independent.

    [Exploreu els cursos :material-arrow-down:](#catalog-of-nextflow-training-courses){ .md-button .md-button--primary .mt-1 }

-   :material-account-group-outline:{ .lg .middle } __Esdeveniments de formació__

    ---

    **Busqueu alguna cosa més enllà de l'autoservei?**

    Trobeu esdeveniments de formació estructurats, orientació per organitzar les vostres pròpies formacions, i la nostra llicència de codi obert i política de contribució.

    [Vegeu els esdeveniments de formació :material-arrow-right:](training_events.md){ .md-button .md-button--secondary .mt-1 }

</div>

!!! note "Traducció assistida per IA"

    Aquesta traducció ha estat creada utilitzant intel·ligència artificial i revisada per traductors humans.
    Agraïm els vostres comentaris i suggeriments de millora.
    Consulteu la nostra [guia de traducció](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md) per a més informació.

## Catalog of Nextflow training courses

<div class="grid cards" markdown>

-   :material-account:{ .lg .middle } __Per a usuaris__

    ---

    ### :material-play-circle:{.nextflow-primary} Executar pipelines {.mt-1}

    Apreneu a executar pipelines existents sense escriure cap codi.

    ??? courses "**Nextflow Run:** Executeu pipelines amb Nextflow"

        Una introducció ràpida a l'execució de pipelines de Nextflow que no requereix comprendre el codi. Cobreix el llançament de pipelines, la recuperació de sortides, l'ús de contenidors i la configuració de l'execució a un nivell bàsic.

        [Veure la formació :material-arrow-right:](nextflow_run/index.md){ .md-button .md-button--secondary }

    ??? courses "**Use nf-core:** Trobeu i executeu pipelines curades per la comunitat"

        Una introducció ràpida a la cerca, execució i configuració de pipelines del projecte comunitari nf-core, començant amb un pipeline de demostració mínim i escalant fins a un pipeline d'anàlisi a escala de producció.

        [Veure la formació :material-arrow-right:](nfcore_use/index.md){ .md-button .md-button--secondary }

    ??? courses "**Scale with Seqera:** Llanceu i monitoritzeu pipelines a escala"

        Una introducció pràctica al llançament i la monitorització de pipelines de Nextflow amb Seqera Platform, tant des de la interfície web com des de la línia de comandes.

        [Veure la formació :material-arrow-right:](seqera_scale/index.md){ .md-button .md-button--secondary }

    ---

    ### :material-tune:{.nextflow-primary} Gestionar l'execució {.mt-1}

    Apreneu a gestionar l'execució de pipelines de manera eficaç.

    ??? courses "**Execution Config:** Configureu pipelines com un professional"

        Una introducció pràctica a la configuració de l'execució de pipelines de Nextflow: adaptació a diferents entorns de càlcul, control de l'assignació de recursos i reintents, i canvi entre perfils de configuració predefinits.

        [Veure la formació :material-arrow-right:](execution_config/index.md){ .md-button .md-button--secondary }

    !!! info compact "Més temes properament"

        L'ajust del rendiment, l'execució en HPC/cloud i més estan planificats per a aquesta secció.
        Voteu sobre què cobrir a continuació en la nostra [breu enquesta d'interès](https://seqera.typeform.com/to/JCs91e8v).

-   :material-code-tags:{ .lg .middle } __Per a desenvolupadors__

    ---

    ### :material-wrench:{.nextflow-primary} Escriure pipelines {.mt-1}

    Apreneu a desenvolupar els vostres propis pipelines de Nextflow.

    ??? courses "**Hello Nextflow:** Desenvolupeu els vostres propis pipelines des de zero"

        Aquest curs cobreix els components principals del llenguatge Nextflow amb prou detall per permetre el desenvolupament de pipelines senzills però completament funcionals, a més d'elements clau de disseny, desenvolupament i pràctiques de configuració de pipelines.

        [Veure la formació :material-arrow-right:](hello_nextflow/index.md){ .md-button .md-button--secondary }

    ??? courses "**Build with nf-core:** Utilitzeu les eines i les regles de nf-core"

        Per a desenvolupadors de Nextflow que vulguin aprendre a desenvolupar pipelines compatibles amb [nf-core](https://nf-co.re/).
        El curs cobreix l'estructura dels pipelines de nf-core amb prou detall per permetre el desenvolupament de pipelines senzills però completament funcionals que aprofiten la plantilla de nf-core i les millors pràctiques de desenvolupament, així com l'ús de mòduls de nf-core existents.

        [Veure la formació :material-arrow-right:](nfcore_build/index.md){ .md-button .md-button--secondary }

    ??? catalog "**Side Quests:** Endinseuvos en temes avançats de Nextflow"

        Una col·lecció de mini-cursos independents destinats a desenvolupadors de Nextflow que vulguin ampliar el seu ventall i/o aprofundir les seves habilitats en temes concrets.
        Es presenten de manera lineal però es poden fer en qualsevol ordre (vegeu les dependències a la descripció general de cada mini-curs).

        [Exploreu els Side Quests :material-arrow-right:](side_quests/index.md){ .md-button .md-button--secondary }

    ---

    ### :material-microscope:{.nextflow-primary} Nextflow per a la Ciència {.mt-1}

    Apreneu a desenvolupar pipelines de Nextflow per a aplicacions científiques específiques.

    ??? courses "**Genomics:** Desenvolupeu un pipeline de detecció de variants"

        Un curs per a investigadors que vulguin aprendre a desenvolupar els seus propis pipelines de genòmica, utilitzant un cas d'ús de detecció de variants per demostrar patrons essencials de desenvolupament amb Nextflow.

        [Veure la formació :material-arrow-right:](nf4_science/genomics/index.md){ .md-button .md-button--secondary }

    ??? courses "**RNAseq:** Desenvolupeu un pipeline de processament de RNAseq massiu"

        Un curs per a investigadors que vulguin aprendre a desenvolupar els seus propis pipelines de RNAseq, utilitzant un cas d'ús de processament de RNAseq massiu per demostrar patrons essencials de desenvolupament amb Nextflow.

        [Veure la formació :material-arrow-right:](nf4_science/rnaseq/index.md){ .md-button .md-button--secondary }

    ??? courses "**Bioimaging:** Executeu i configureu pipelines d'imatge"

        Un curs per a investigadors que vulguin aprendre a executar i configurar pipelines de bioimatge, utilitzant nf-core/molkart per demostrar patrons essencials d'ús de Nextflow.

        [Veure la formació :material-arrow-right:](nf4_science/imaging/index.md){ .md-button .md-button--secondary }

</div>

## Configuració i ajuda

<div class="grid cards mb-4" markdown>

-   :material-cog-outline:{ .lg .middle } __Entorn de formació__

    ---

    Opcions per configurar el vostre entorn per a les formacions de Nextflow.

    [Vegeu els entorns de formació :material-arrow-right:](envsetup/index.md){ .md-button .md-button--secondary }

-   :material-tag-outline:{ .lg .middle } __Versions de Nextflow__

    ---

    Comprensió i gestió de l'evolució de les versions de sintaxi de Nextflow.

    [Comproveu els requisits de versió :material-arrow-right:](info/nxf_versions.md){ .md-button .md-button--secondary }

-   :material-file-code-outline:{ .lg .middle } __El pipeline Hello__

    ---

    Resum de què fa el pipeline Hello i com està estructurat.

    [Llegiu el resum :material-arrow-right:](info/hello_pipeline.md){ .md-button .md-button--secondary }

-   :material-lifebuoy:{ .lg .middle } __Obtenir ajuda__

    ---

    Recursos útils quan teniu un problema amb la formació de Nextflow.

    [Trobeu ajuda :material-arrow-right:](help.md){ .md-button .md-button--secondary }

</div>
