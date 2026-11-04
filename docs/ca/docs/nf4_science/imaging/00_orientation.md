# Orientació

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Traducció assistida per IA - [més informació i suggeriments](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

Aquesta orientació assumeix que ja heu obert l'entorn de formació fent clic al botó "Open in GitHub Codespaces".
Si no ho heu fet, feu-ho ara, idealment en una segona finestra o pestanya del navegador perquè pugueu consultar aquestes instruccions.

[![Open in GitHub Codespaces](https://github.com/codespaces/badge.svg)](https://github.com/codespaces/new?hide_repo_select=true&ref=master&repo=290790519&skip_quickstart=true&machine=premiumLinux&devcontainer_path=.devcontainer%2Fdevcontainer.json)

!!!warning "Requisit de mida de màquina"

    Assegureu-vos de seleccionar una **màquina de 8 nuclis** quan creeu el vostre Codespace per a aquest curs de formació. Els workflows de bioimatge requereixen recursos de càlcul addicionals.

## GitHub Codespaces

L'entorn de GitHub Codespaces conté tot el programari, codi i dades necessaris per treballar en aquest curs de formació, així que no cal que instal·leu res vosaltres mateixos.
No obstant això, necessiteu un compte de GitHub (gratuït) per iniciar sessió, i si no esteu familiaritzats amb la interfície, hauríeu de dedicar uns minuts a familiaritzar-vos-hi completant el mini-curs d'[Orientació de GitHub Codespaces](../../envsetup/index.md).

## Pre-descàrrega d'imatges Docker

Un cop hàgiu obert el vostre Codespace, descarreguem totes les imatges Docker que necessitarem per a aquest curs de formació.
Això estalviarà temps més endavant i garantirà una execució fluida dels workflows.

Obriu una nova pestanya de terminal i executeu la comanda següent:

```bash
nextflow run nf-core/molkart -profile docker,test -stub -resume --outdir results
```

Aquesta comanda descarregarà totes les imatges Docker necessàries en segon pla.
Podeu continuar amb la resta de l'orientació mentre s'executa.

!!!tip

    La bandera `-stub` permet que el pipeline s'executi ràpidament sense processar dades reals, cosa que és perfecta per descarregar imatges. Podeu monitoritzar el progrés a la pestanya del terminal.

## Directori de treball

Al llarg d'aquest curs de formació, treballarem al directori `nf4-science/imaging/`.

Canvieu de directori ara executant aquesta comanda al terminal:

```bash
cd nf4-science/imaging/
```

!!!tip

    Si per qualsevol motiu sortiu d'aquest directori, sempre podeu utilitzar el camí complet per tornar-hi, assumint que esteu executant això dins de l'entorn de formació de GitHub Codespaces:

    ```bash
    cd /workspaces/training/nf4-science/imaging
    ```

**Ara, per començar el curs, feu clic a la fletxa de la cantonada inferior dreta d'aquesta pàgina.**
