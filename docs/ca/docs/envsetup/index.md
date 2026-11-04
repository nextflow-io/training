---
title: Environment options
description: Options for setting up your environment for the Nextflow trainings
hide:
  - toc
  - footer
---

# Opcions d'entorn

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Traducció assistida per IA - [més informació i suggeriments](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

Volem proporcionar un entorn consistent i exhaustivament provat que permeti als estudiants centrar-se en aprendre Nextflow sense haver de dedicar temps i esforç a gestionar programari.
Amb aquesta finalitat, hem desenvolupat un entorn contenidoritzat que conté tot el programari necessari, els fitxers de codi i les dades d'exemple per treballar amb tots els nostres cursos.

Aquest entorn contenidoritzat es pot executar directament a Github Codespaces o localment a VS Code amb l'extensió Devcontainers.

<div class="grid cards" markdown>

-   :material-cloud-outline:{ .lg .middle } __Github Codespaces__

    ---

    GitHub Codespaces és un servei basat en web que ens permet proporcionar un entorn preconfigurat per a la formació, amb totes les eines i dades incloses, recolzat per màquines virtuals al núvol. És accessible gratuïtament per a qualsevol persona amb un compte de Github.

    [Utilitzar Github Codespaces:material-arrow-right:](01_setup.md){ .md-button .md-button--primary .mt-1 }

-   :material-laptop:{ .lg .middle } __Devcontainers locals__

    ---

    VS Code amb Devcontainers proporciona un entorn de desenvolupament contenidoritzat d'execució local amb totes les eines de formació preconfigurades. Ofereix el mateix entorn preconfigurat que Codespaces però executant-se completament al vostre maquinari local.

    [Utilitzar Devcontainers localment :material-arrow-right:](03_devcontainer.md){ .md-button .md-button--primary .mt-1 }

</div>

## Instruccions per a la instal·lació manual

Si cap de les opcions anteriors s'ajusta a les vostres necessitats, podeu replicar aquest entorn al vostre propi sistema local instal·lant les dependències de programari manualment i clonant el repositori de formació.

[Instal·lació manual :material-arrow-right:](02_local.md){ .md-button .md-button--primary .mt-1 }

---

!!! info "Desaprovació de Gitpod"

    La formació de Nextflow utilitzava [Gitpod](https://gitpod.io) fins al febrer de 2025.
    No obstant això, els creadors de Gitpod van decidir retirar la funcionalitat gratuïta en favor del sistema [Gitpod Flex](https://www.gitpod.io/blog/introducing-gitpod-flex).
    Per aquest motiu, hem canviat a l'ús de GitHub Codespaces, que també ofereix un entorn de desenvolupament d'un sol clic sense configuració prèvia.

    Depenent de quan us vau registrar a Gitpod i de quan exactament retirin el servei, és possible que encara pugueu iniciar la formació al seu antic IDE al núvol, tot i que no podem garantir un accés fiable en el futur:
    [Obrir a Gitpod](https://gitpod.io/#https://github.com/nextflow-io/training).
