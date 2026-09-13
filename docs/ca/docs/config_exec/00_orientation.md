# Primers passos

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Traducció assistida per IA - [més informació i suggeriments](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>


## Iniciar un entorn de formació

Per utilitzar l'entorn preconfigurado que proporcionem a GitHub Codespaces, feu clic al botó "Open in GitHub Codespaces" que trobareu a continuació. Per a altres opcions, consulteu [Opcions d'entorn](../envsetup/index.md).

Us recomanem obrir l'entorn de formació en una nova pestanya o finestra del navegador (feu clic amb el botó dret, ctrl-clic o cmd-clic segons el vostre equip) per poder continuar llegint mentre es carrega l'entorn.
Haureu de mantenir aquestes instruccions obertes en paral·lel per treballar el curs.

[![Open in GitHub Codespaces](https://github.com/codespaces/badge.svg)](https://codespaces.new/nextflow-io/training?quickstart=1&ref=master)

### Conceptes bàsics de l'entorn

Aquest entorn de formació conté tot el programari, el codi i les dades necessàries per treballar el curs de formació, de manera que no cal que instal·leu res vosaltres mateixos.

El codespace està configurat amb una interfície VSCode, que inclou un explorador de fitxers, un editor de codi i un terminal.
Totes les instruccions donades durant el curs (p. ex., "obriu el fitxer", "editeu el codi" o "executeu aquesta comanda") fan referència a aquestes tres parts de la interfície VSCode, llevat que s'indiqui el contrari.

Si esteu treballant aquest curs pel vostre compte, familiaritzeu-vos amb els [conceptes bàsics de l'entorn](../envsetup/01_setup.md) per obtenir més detalls.

### Requisits de versió

Aquest curs requereix Nextflow 25.10.2 o posterior, amb l'analitzador de sintaxi v2 activat (el valor per defecte a partir de la versió 25.10+).
Si esteu utilitzant un entorn local o personalitzat, assegureu-vos d'estar utilitzant la configuració correcta tal com es documenta [aquí](../info/nxf_versions.md).

## Prepareu-vos per treballar

Un cop el vostre codespace estigui en funcionament, hi ha dues coses a fer abans de començar: establir el directori de treball i fer una ullada als materials proporcionats.

### Establir el directori de treball

Per defecte, el codespace s'obre a l'arrel de tots els cursos de formació.
Per a aquest curs, canvieu al directori `config-exec/`:

```bash
cd config-exec/
```

A continuació, configureu VSCode perquè es centri en aquest directori, de manera que només apareguin els fitxers rellevants a la barra lateral de l'explorador de fitxers:

```bash
code .
```

!!! tip "Consell"

    Si per qualsevol motiu sortiu d'aquest directori (p. ex., el vostre codespace entra en repòs), sempre podeu utilitzar la ruta completa per tornar-hi, assumint que esteu executant-ho dins de l'entorn de formació de Github Codespaces:

    ```bash
    cd /workspaces/training/config-exec
    ```

### Explorar els materials proporcionats

Podeu explorar els materials del curs utilitzant l'explorador de fitxers de l'esquerra, o amb la comanda `tree`.
Executeu el següent des del terminal per veure l'estructura completa:

```bash
tree . -L 2
```

??? abstract "Contingut del directori"

    ```console
    .
    ├── data
    │   └── greetings.csv
    ├── main.nf
    ├── modules
    │   ├── collectGreetings.nf
    │   ├── convertToUpper.nf
    │   ├── cowpy.nf
    │   └── sayHello.nf
    └── nextflow.config
    ```

Els fitxers **`main.nf`** i **`modules/`** són el mateix pipeline de múltiples passos de [Nextflow Run](../nextflow_run/index.md), i el fitxer **`nextflow.config`** és la mateixa configuració que ja heu vist allà.
Ampliareu tots dos al llarg d'aquests exercicis.

El directori **`data/`** conté el fitxer d'entrada CSV que llegeix el pipeline.

## Llista de verificació de preparació

Creieu que esteu a punt per començar?

- [ ] Entenc l'objectiu d'aquest curs i els seus prerequisits
- [ ] El meu entorn està en funcionament
- [ ] He establert el meu directori de treball correctament

Si podeu marcar totes les caselles, esteu a punt per continuar.

**Per continuar a la [Part 1: Adaptació al vostre entorn de còmput](./01_packaging_and_execution.md), feu clic a la fletxa a la cantonada inferior dreta d'aquesta pàgina.**
