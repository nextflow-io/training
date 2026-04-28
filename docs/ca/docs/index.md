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

    Els cursos de formació que es llisten a continuació estan dissenyats per ser utilitzats com a recurs d'autoservei.
    Podeu treballar-hi pel vostre compte en qualsevol moment, ja sigui en l'entorn basat en web que proporcionem mitjançant Github Codespaces o en el vostre propi entorn.

    [Exploreu els cursos :material-arrow-right:](#catalog-of-nextflow-training-courses){ .md-button .md-button--primary .mt-1 }

-   :material-information-outline:{ .lg .middle } __Informació addicional__

    ---

    ??? warning "Compatibilitat de versions"

        <!-- Qualsevol actualització d'aquest contingut s'ha de copiar a la pàgina d'instal·lació local -->
        **A partir del gener de 2026, tots els nostres cursos de formació en Nextflow requereixen la versió 25.10.2 de Nextflow o posterior, amb la sintaxi estricta activada, tret que s'indiqui el contrari.**

        Per obtenir més informació sobre els requisits de versió i la sintaxi estricta, consulteu la [guia de migració de la documentació de Nextflow](https://nextflow.io/docs/latest/strict-syntax.html).

        Les versions anteriors del material de formació corresponents a la sintaxi prèvia estan disponibles mitjançant el selector de versions a la barra de menú d'aquesta pàgina web.

    ??? terminal "Opcions d'entorn"

        Proporcionem un entorn de formació basat en web on tot el que necessiteu per fer la formació està preinstal·lat, disponible a través de Github Codespaces (requereix un compte de GitHub gratuït).

        [![Obriu a GitHub Codespaces](https://github.com/codespaces/badge.svg)](https://codespaces.new/nextflow-io/training?quickstart=1&ref=master)

        Si això no s'adapta a les vostres necessitats, consulteu les altres [opcions d'entorn](./envsetup/index.md).

    ??? learning "Esdeveniments de formació"

        Si preferiu fer la formació en Nextflow com a part d'un esdeveniment estructurat, hi ha moltes oportunitats per fer-ho. Us recomanem que consulteu les opcions següents:

        - **[Setmanes de formació]()** organitzades trimestralment per l'equip de la Comunitat
        - **[Esdeveniments de Seqera](https://seqera.io/events/)** inclouen esdeveniments de formació presencial organitzats per Seqera (cerqueu 'Seqera Sessions' i 'Nextflow Summit')
        - **[Ambaixadors de Nextflow]()** organitzen esdeveniments per a la seva comunitat local
        - **[Esdeveniments de nf-core](https://nf-co.re/events)** inclouen hackathons de la comunitat

    ??? people "Informació per a formadors"

        Si sou un instructor que organitza les vostres pròpies formacions, podeu utilitzar els nostres materials directament des del portal de formació sempre que atribuïu el crèdit corresponent. Vegeu 'Crèdits i contribucions' a continuació per obtenir més detalls.

        A més, ens encantaria saber com podríem donar-vos millor suport en els vostres esforços de formació. Poseu-vos en contacte amb nosaltres a [community@seqera.io](mailto:community@seqera.io) o al fòrum de la comunitat (vegeu la pàgina d'[Ajuda](help.md)).

    ??? licensing "Llicència de codi obert i política de contribució"

        [![Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International (CC BY-NC-SA 4.0)](assets/img/cc_by-nc-sa.svg){ align=right }](https://creativecommons.org/licenses/by-nc-sa/4.0/)

        Aquest material de formació és desenvolupat i mantingut per [Seqera](https://seqera.io) i publicat sota una llicència de codi obert ([CC BY-NC-SA](https://creativecommons.org/licenses/by-nc-sa/4.0/)) en benefici de la comunitat. Si voleu utilitzar aquest material d'una manera que quedi fora de l'àmbit de la llicència (tingueu en compte les limitacions sobre l'ús comercial i la redistribució), poseu-vos en contacte amb nosaltres a [community@seqera.io](mailto:community@seqera.io) per parlar de la vostra sol·licitud.

        Acceptem millores, correccions i informes d'errors de la comunitat. Cada pàgina té una icona :material-file-edit-outline: a la part superior dreta que enllaça amb el repositori de codi, on podeu informar de problemes o proposar canvis al material de formació font mitjançant una pull request. Vegeu el fitxer `README.md` al repositori per obtenir més detalls.

</div>

!!! note "Traducció assistida per IA"

    Aquesta traducció ha estat creada utilitzant intel·ligència artificial i revisada per traductors humans.
    Agraïm els vostres comentaris i suggeriments de millora.
    Consulteu la nostra [guia de traducció](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md) per a més informació.

## Catàleg de cursos de formació en Nextflow

<div class="grid cards" markdown>

-   :material-walk:{ .lg .middle } __Itinerari introductori__

    ---

    ### :material-compass:{.nextflow-primary} Nextflow per a nouvinguts {.mt-1}

    Cursos independents del domini destinats a aquells que són completament nous a Nextflow. Cada curs consisteix en una sèrie de mòduls de formació dissenyats per ajudar els estudiants a desenvolupar les seves habilitats de manera progressiva.

    ??? courses "**Hello Nextflow:** Apreneu a desenvolupar els vostres propis pipelines"

        Aquest curs cobreix els components principals del llenguatge Nextflow amb prou detall per permetre el desenvolupament de pipelines simples però completament funcionals, a més dels elements clau de disseny, desenvolupament i pràctiques de configuració de pipelines.

        [Comenceu la formació Hello Nextflow :material-arrow-right:](hello_nextflow/index.md){ .md-button .md-button--secondary }

    ??? courses "**Nextflow Run:** Apreneu a executar pipelines existents"

        Una introducció concisa a l'execució i configuració de pipelines Nextflow, basada en el curs per a desenvolupadors Hello Nextflow però amb menys èmfasi en el codi. Cobreix l'execució, les sortides, l'estructura bàsica del codi i la configuració per a diferents entorns de càlcul.

        [Comenceu la formació Nextflow Run :material-arrow-right:](nextflow_run/index.md){ .md-button .md-button--secondary }

    ---

    ### :material-microscope:{.nextflow-primary} Nextflow per a la ciència {.mt-1}

    Apreneu a aplicar els conceptes i components presentats a 'Hello Nextflow' a casos d'ús científics específics.

    ??? courses "**Nextflow per a Genòmica** (detecció de variants)"

        Per a investigadors que volen aprendre a desenvolupar els seus propis pipelines de genòmica. El curs utilitza un cas d'ús de detecció de variants per demostrar com desenvolupar un pipeline de genòmica simple però funcional.

        [Comenceu la formació Nextflow per a Genòmica :material-arrow-right:](nf4_science/genomics/index.md){ .md-button .md-button--secondary }

    ??? courses "**Nextflow per a RNAseq** (RNAseq massiu)"

        Per a investigadors que volen aprendre a desenvolupar els seus propis pipelines de RNAseq. El curs utilitza un cas d'ús de processament de RNAseq massiu per demostrar com desenvolupar un pipeline de RNAseq simple però funcional.

        [Comenceu la formació Nextflow per a RNAseq :material-arrow-right:](nf4_science/rnaseq/index.md){ .md-button .md-button--secondary }

    ??? courses "**Nextflow per a Imatge** (òmica espacial)"

        Per a investigadors en imatge i òmica espacial que volen aprendre a executar i personalitzar pipelines d'anàlisi. El curs utilitza el pipeline nf-core/molkart per proporcionar un pipeline biològicament rellevant que demostra com executar, configurar i gestionar les entrades per a workflows de pipelines Nextflow.

        [Comenceu la formació Nextflow per a Imatge :material-arrow-right:](nf4_science/imaging/index.md){ .md-button .md-button--secondary }

-   :material-run:{ .lg .middle } __Itinerari avançat__

    ---

    ### :material-bridge:{.nextflow-primary} De Nextflow a nf-core {.mt-1}

    Apreneu a utilitzar codi i bones pràctiques del projecte de la comunitat [nf-core](https://nf-co.re/).

    Aquests cursos us ajuden a passar dels fonaments de Nextflow a les bones pràctiques de nf-core.
    Enteneu com i per què la comunitat nf-core construeix pipelines, i com podeu contribuir-hi i reutilitzar aquestes tècniques.

    ??? courses "**Hello nf-core:** Primers passos amb nf-core"

        Per a desenvolupadors que volen aprendre a executar i desenvolupar pipelines compatibles amb [nf-core](https://nf-co.re/). El curs cobreix l'estructura dels pipelines nf-core amb prou detall per permetre el desenvolupament de pipelines simples però completament funcionals que segueixen la plantilla nf-core i les bones pràctiques de desenvolupament, així com l'ús de mòduls nf-core existents.

        [Comenceu la formació Hello nf-core :material-arrow-right:](hello_nf-core/index.md){ .md-button .md-button--secondary }

    ---

    ### :material-rocket-launch:{.nextflow-primary} Formació avançada en Nextflow {.mt-1}

    Apreneu conceptes i mecanismes avançats per desenvolupar i desplegar pipelines Nextflow per abordar casos d'ús del món real.

    ??? courses "**Side Quests:** Immersions en temes independents"

        Mini-cursos independents destinats a desenvolupadors de Nextflow que volen ampliar el seu ventall i/o aprofundir les seves habilitats en temes concrets. Es presenten de manera lineal però es poden fer en qualsevol ordre (vegeu les dependències a la visió general de cada mini-curs).

        [Exploreu els Side Quests :material-arrow-right:](side_quests/index.md){ .md-button .md-button--secondary }

    ??? courses "**Col·leccions de formació:** Itineraris d'aprenentatge recomanats a través dels Side Quests"

        Les col·leccions de formació combinen múltiples Side Quests per proporcionar una experiència d'aprenentatge completa al voltant d'un tema o cas d'ús particular.

        [Exploreu les col·leccions de formació :material-arrow-right:](training_collections/index.md){ .md-button .md-button--secondary }

</div>

!!! info "Busqueu materials de formació arxivats?"

    Els materials de formació antics (Formació Fonamental, Formació Avançada i altres cursos experimentals) s'han eliminat del portal de formació perquè són incompatibles amb la sintaxi estricta de Nextflow 3.0.
    Si necessiteu accedir a aquests materials, estan disponibles a l'[historial de git](https://github.com/nextflow-io/training) anterior al gener de 2026.

---

<div markdown class="homepage_logos">

![Seqera](assets/img/seqera_logo.png#only-light)

![Seqera](assets/img/seqera_logo_dark.png#only-dark)

</div>
