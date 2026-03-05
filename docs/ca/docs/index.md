---
title: Home
description: Welcome to the Nextflow community training portal!
hide:
  - toc
  - footer
---

# Formació de Nextflow

<div class="grid cards" markdown>

-   :material-book-open-variant:{ .lg .middle } __Cursos d'autoservei__

    ---

    **Et donem la benvinguda al portal de formació de la comunitat Nextflow!**

    Els cursos de formació que es llisten a continuació estan dissenyats per ser utilitzats com a recurs d'autoservei.
    Pots treballar-los pel teu compte en qualsevol moment, ja sigui a l'entorn basat en web que proporcionem via Github Codespaces o al teu propi entorn.

    [Explora els cursos :material-arrow-right:](#catalog-of-nextflow-training-courses){ .md-button .md-button--primary .mt-1 }

-   :material-information-outline:{ .lg .middle } __Informació addicional__

    ---

    ??? warning "Compatibilitat de versions"

        <!-- Any update to this content needs to be copied to the local installation page -->
        **A partir de gener de 2026, tots els nostres cursos de formació de Nextflow requereixen la versió 25.10.2 o posterior de Nextflow, amb la sintaxi estricta activada, tret que s'indiqui el contrari.**

        Per a més informació sobre els requisits de versió i la sintaxi estricta, consulteu la [guia de migració de la documentació de Nextflow](https://nextflow.io/docs/latest/strict-syntax.html).

        Les versions antigues del material de formació corresponents a la sintaxi anterior estan disponibles mitjançant el selector de versions a la barra de menú d'aquesta pàgina web.

    ??? terminal "Opcions d'entorn"

        Proporcionem un entorn de formació basat en web on tot el que necessites per fer la formació està preinstal·lat, disponible a través de Github Codespaces (requereix un compte gratuït de GitHub).

        [![Open in GitHub Codespaces](https://github.com/codespaces/badge.svg)](https://codespaces.new/nextflow-io/training?quickstart=1&ref=master)

        Si això no s'ajusta a les teves necessitats, consulta les altres [Opcions d'entorn](./envsetup/index.md).

    ??? learning "Esdeveniments de formació"

        Si prefereixes fer la formació de Nextflow com a part d'un esdeveniment estructurat, hi ha moltes oportunitats per fer-ho. Recomanem consultar les següents opcions:

        - **[Training Weeks]()** organitzades trimestralment per l'equip de la Comunitat
        - **[Seqera Events](https://seqera.io/events/)** inclouen esdeveniments de formació presencials organitzats per Seqera (cerca 'Seqera Sessions' i 'Nextflow Summit')
        - **[Nextflow Ambassadors]()** organitzen esdeveniments per a la seva comunitat local
        - **[nf-core events](https://nf-co.re/events)** inclouen hackathons de la comunitat

    ??? people "Informació per a formadors"

        Si ets un instructor que organitza les seves pròpies formacions, ets benvingut/da a utilitzar els nostres materials directament des del portal de formació sempre que atribueixis el crèdit adequat. Consulta 'Crèdits i contribucions' a continuació per a més detalls.

        A més, ens encantaria saber de tu sobre com podríem donar millor suport als teus esforços de formació! Contacta'ns a [community@seqera.io](mailto:community@seqera.io) o al fòrum de la comunitat (consulta la pàgina d'[Ajuda](help.md)).

    ??? licensing "Llicència de codi obert i política de contribució"

        [![Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International (CC BY-NC-SA 4.0)](assets/img/cc_by-nc-sa.svg){ align=right }](https://creativecommons.org/licenses/by-nc-sa/4.0/)

        Aquest material de formació està desenvolupat i mantingut per [Seqera](https://seqera.io) i publicat sota una llicència de codi obert ([CC BY-NC-SA](https://creativecommons.org/licenses/by-nc-sa/4.0/)) en benefici de la comunitat. Si vols utilitzar aquest material d'una manera que quedi fora de l'abast de la llicència (tingues en compte les limitacions sobre l'ús comercial i la redistribució), contacta'ns a [community@seqera.io](mailto:community@seqera.io) per discutir la teva sol·licitud.

        Donem la benvinguda a millores, correccions i informes d'errors de la comunitat. Cada pàgina té una icona :material-file-edit-outline: a la part superior dreta de la pàgina que enllaça amb el repositori de codi, on pots informar de problemes o proposar canvis al material font de formació mitjançant una pull request. Consulta el `README.md` al repositori per a més detalls.

</div>

!!! note "Traducció assistida per IA"

    Aquesta traducció ha estat creada utilitzant intel·ligència artificial i revisada per traductors humans.
    Agraïm els vostres comentaris i suggeriments de millora.
    Consulteu la nostra [guia de traducció](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md) per a més informació.

## Catàleg de cursos de formació de Nextflow

<div class="grid cards" markdown>

-   :material-walk:{ .lg .middle } __Itinerari introductori__

    ---

    ### :material-compass:{.nextflow-primary} Nextflow per a Principiants {.mt-1}

    Cursos independents del domini destinats a aquells que són completament nous a Nextflow. Cada curs consisteix en una sèrie de mòduls de formació dissenyats per ajudar els estudiants a desenvolupar les seves habilitats progressivament.

    ??? courses "**Hello Nextflow:** Aprèn a desenvolupar els teus propis pipelines"

        Aquest curs cobreix els components principals del llenguatge Nextflow amb prou detall per permetre desenvolupar pipelines simples però completament funcionals, a més d'elements clau de disseny, desenvolupament i pràctiques de configuració de pipelines.

        [Comença la formació Hello Nextflow :material-arrow-right:](hello_nextflow/index.md){ .md-button .md-button--secondary }

    ??? courses "**Nextflow Run:** Aprèn a executar pipelines existents"

        Una introducció concisa a l'execució i configuració de pipelines de Nextflow, basada en el curs de desenvolupadors Hello Nextflow però amb menys èmfasi en el codi. Cobreix l'execució, les sortides, l'estructura bàsica del codi i la configuració per a diferents entorns de càlcul.

        [Comença la formació Nextflow Run :material-arrow-right:](nextflow_run/index.md){ .md-button .md-button--secondary }

    ---

    ### :material-microscope:{.nextflow-primary} Nextflow per a la Ciència {.mt-1}

    Aprèn a aplicar els conceptes i components presentats a 'Hello Nextflow' a casos d'ús científics específics.

    ??? courses "**Nextflow for Genomics** (detecció de variants)"

        Per a investigadors que vulguin aprendre a desenvolupar els seus propis pipelines de genòmica. El curs utilitza un cas d'ús de detecció de variants per demostrar com desenvolupar un pipeline de genòmica simple però funcional.

        [Comença la formació Nextflow for Genomics :material-arrow-right:](nf4_science/genomics/){ .md-button .md-button--secondary }

    ??? courses "**Nextflow for RNAseq** (RNAseq massiu)"

        Per a investigadors que vulguin aprendre a desenvolupar els seus propis pipelines de RNAseq. El curs utilitza un cas d'ús de processament de RNAseq massiu per demostrar com desenvolupar un pipeline de RNAseq simple però funcional.

        [Comença la formació Nextflow for RNAseq :material-arrow-right:](nf4_science/rnaseq/){ .md-button .md-button--secondary }

    ??? courses "**Nextflow for Imaging** (òmica espacial)"

        Per a investigadors en imatge i òmica espacial que vulguin aprendre a executar i personalitzar pipelines d'anàlisi. El curs utilitza el pipeline nf-core/molkart per proporcionar un pipeline biològicament rellevant que demostra com executar, configurar i gestionar entrades per a workflows de pipelines de Nextflow.

        [Comença la formació Nextflow for Imaging :material-arrow-right:](nf4_science/imaging/){ .md-button .md-button--secondary }

-   :material-run:{ .lg .middle } __Itinerari avançat__

    ---

    ### :material-bridge:{.nextflow-primary} De Nextflow a nf-core {.mt-1}

    Aprèn a utilitzar codi i bones pràctiques del projecte comunitari [nf-core](https://nf-co.re/).

    Aquests cursos t'ajuden a passar dels fonaments de Nextflow a les bones pràctiques d'nf-core.
    Comprèn com i per què la comunitat nf-core construeix pipelines, i com pots contribuir i reutilitzar aquestes tècniques.

    ??? courses "**Hello nf-core:** Comença amb nf-core"

        Per a desenvolupadors que vulguin aprendre a executar i desenvolupar pipelines compatibles amb [nf-core](https://nf-co.re/). El curs cobreix l'estructura dels pipelines nf-core amb prou detall per permetre desenvolupar pipelines simples però completament funcionals que segueixen la plantilla nf-core i les bones pràctiques de desenvolupament, així com utilitzar mòduls nf-core existents.

        [Comença la formació Hello nf-core :material-arrow-right:](hello_nf-core/index.md){ .md-button .md-button--secondary }

    ---

    ### :material-rocket-launch:{.nextflow-primary} Formació Avançada de Nextflow {.mt-1}

    Aprèn conceptes i mecanismes avançats per desenvolupar i desplegar pipelines de Nextflow per abordar casos d'ús del món real.

    ??? courses "**Side Quests:** Aprofundiments en temes independents"

        Mini-cursos independents destinats a desenvolupadors de Nextflow que vulguin ampliar el seu rang i/o aprofundir les seves habilitats en temes particulars. Es presenten linealment però es poden fer en qualsevol ordre (consulta les dependències a la visió general de cada mini-curs).

        [Explora els Side Quests :material-arrow-right:](side_quests/){ .md-button .md-button--secondary }

    ??? courses "**Training Collections:** Itineraris d'aprenentatge recomanats a través dels Side Quests"

        Les Training Collections combinen múltiples Side Quests per proporcionar una experiència d'aprenentatge completa al voltant d'un tema o cas d'ús particular.

        [Explora les Training Collections :material-arrow-right:](training_collections/){ .md-button .md-button--secondary }

</div>

!!! info "Busques materials de formació arxivats?"

    Els materials de formació antics (Fundamentals Training, Advanced Training i altres cursos experimentals) s'han eliminat del portal de formació ja que són incompatibles amb la sintaxi estricta de Nextflow 3.0.
    Si necessites accés a aquests materials, estan disponibles a l'[historial de git](https://github.com/nextflow-io/training) abans de gener de 2026.

---

<div markdown class="homepage_logos">

![Seqera](assets/img/seqera_logo.png#only-light)

![Seqera](assets/img/seqera_logo_dark.png#only-dark)

</div>
