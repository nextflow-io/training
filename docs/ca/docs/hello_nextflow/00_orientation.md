# Primers passos

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Traducció assistida per IA - [més informació i suggeriments](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

<div class="video-wrapper">
  <iframe width="560" height="315" src="https://www.youtube.com/embed/PIjOdFaYwWA?si=gZxlXgkVxuLEzOsC" title="YouTube video player" frameborder="0" allow="accelerometer; autoplay; clipboard-write; encrypted-media; gyroscope; picture-in-picture; web-share" referrerpolicy="strict-origin-when-cross-origin" allowfullscreen></iframe>
</div>

/// caption
:fontawesome-brands-youtube:{ .youtube } Consulteu [la llista de reproducció completa](https://youtube.com/playlist?list=PLPZ8WHdZGxmWKozQuzr27jyMGqp9kElVK&si=eF7cLR62goy-lc6n) al canal de YouTube de Nextflow.

:green_book: La transcripció del vídeo està disponible [aquí](./transcripts/00_orientation.md).
///

!!! tip "Consell"

    Els vídeos de YouTube tenen superpoders!

    - :fontawesome-solid-closed-captioning: Subtítols d'alta qualitat (curats manualment). Activeu-los amb la icona :material-subtitles:
    - :material-bookmark: Capítols de vídeo a la línia de temps que corresponen als encapçalaments de la pàgina.

## Inicieu un entorn de formació

Per utilitzar l'entorn preconfigurat que proporcionem a GitHub Codespaces, feu clic al botó "Open in GitHub Codespaces" a continuació. Per a altres opcions, consulteu [Opcions d'entorn](../envsetup/index.md).

Recomanem obrir l'entorn de formació en una pestanya o finestra nova del navegador (utilitzeu clic dret, ctrl-clic o cmd-clic segons el vostre equip) perquè pugueu continuar llegint mentre es carrega l'entorn.
Haureu de mantenir aquestes instruccions obertes en paral·lel per treballar durant el curs.

[![Open in GitHub Codespaces](https://github.com/codespaces/badge.svg)](https://codespaces.new/nextflow-io/training?quickstart=1&ref=master)

### Conceptes bàsics de l'entorn

Aquest entorn de formació conté tot el programari, codi i dades necessaris per treballar durant el curs de formació, així que no cal que instal·leu res vosaltres mateixos.

El codespace està configurat amb una interfície VSCode, que inclou un explorador de fitxers, un editor de codi i un terminal de línia de comandes.
Totes les instruccions donades durant el curs (per exemple, 'obriu el fitxer', 'editeu el codi' o 'executeu aquesta comanda') es refereixen a aquestes tres parts de la interfície VSCode tret que s'especifiqui el contrari.

Si esteu treballant en aquest curs pel vostre compte, familiaritzeu-vos amb els [conceptes bàsics de l'entorn](../envsetup/01_setup.md) per a més detalls.

### Requisits de versió

Aquesta formació està dissenyada per a Nextflow 25.10.2 o posterior **amb l'analitzador de sintaxi v2 ACTIVAT**.
Si utilitzeu un entorn local o personalitzat, assegureu-vos que utilitzeu la configuració correcta tal com es documenta [aquí](../info/nxf_versions.md).

## Prepareu-vos per treballar

Un cop el vostre codespace estigui en funcionament, hi ha dues coses que heu de fer abans d'endinsar-vos en la formació: establir el vostre directori de treball per a aquest curs específic i donar un cop d'ull als materials proporcionats.

### Establiu el directori de treball

Per defecte, el codespace s'obre amb el directori de treball establert a l'arrel de tots els cursos de formació, però per a aquest curs, treballarem al directori `hello-nextflow/`.

Canvieu de directori ara executant aquesta comanda al terminal:

```bash
cd hello-nextflow/
```

Podeu configurar VSCode perquè se centri en aquest directori, de manera que només els fitxers rellevants es mostrin a la barra lateral de l'explorador de fitxers:

```bash
code .
```

!!! tip "Consell"

    Si per qualsevol motiu sortiu d'aquest directori (per exemple, el vostre codespace s'adorm), sempre podeu utilitzar la ruta completa per tornar-hi, assumint que esteu executant això dins de l'entorn de formació de Github Codespaces:

    ```bash
    cd /workspaces/training/hello-nextflow
    ```

Ara donem un cop d'ull als continguts.

### Exploreu els materials proporcionats

Podeu explorar els continguts d'aquest directori utilitzant l'explorador de fitxers al costat esquerre de l'espai de treball de formació.
Alternativament, podeu utilitzar la comanda `tree`.

Durant el curs, utilitzem la sortida de `tree` per representar l'estructura i els continguts del directori de forma llegible, de vegades amb modificacions menors per claredat.

Aquí generem una taula de continguts fins al segon nivell:

```bash
tree . -L 2
```

??? abstract "Contingut del directori"

    ```console
    .
    ├── data
    │   └── greetings.csv
    ├── hello-channels.nf
    ├── hello-config.nf
    ├── hello-containers.nf
    ├── hello-modules.nf
    ├── hello-workflow.nf
    ├── hello-world.nf
    ├── nextflow.config
    ├── solutions
    │   ├── 1-hello-world
    │   ├── 2-hello-channels
    │   ├── 3-hello-workflow
    │   ├── 4-hello-modules
    │   ├── 5-hello-containers
    │   └── 6-hello-config
    ├── test-params.json
    └── test-params.yaml
    ```

Feu clic a la caixa de color per expandir la secció i veure els seus continguts.
Utilitzem seccions desplegables com aquesta per incloure la sortida esperada de les comandes de manera concisa.

- **Els fitxers `.nf`** són scripts de workflow que es nomenen segons la part del curs en què s'utilitzen.

- **El fitxer `nextflow.config`** és un fitxer de configuració que estableix propietats mínimes de l'entorn.
  Podeu ignorar-lo de moment.

- **El fitxer `greetings.csv`** sota `data/` conté dades d'entrada que utilitzarem en la major part del curs. Es descriu a la Part 2 (Canals), quan l'introduïm per primera vegada.

- **Els fitxers `test-params.*`** són fitxers de configuració que utilitzarem a la Part 6 (Configuració). Podeu ignorar-los de moment.

- **El directori `solutions`** conté els scripts de workflow completats que resulten de cada pas del curs.
  Estan pensats per ser utilitzats com a referència per comprovar el vostre treball i resoldre qualsevol problema.

## Llista de verificació de preparació

Creieu que esteu preparats per començar?

- [ ] Entenc l'objectiu d'aquest curs i els seus prerequisits
- [ ] El meu entorn està en funcionament
- [ ] He establert el meu directori de treball adequadament

Si podeu marcar totes les caselles, esteu preparats per començar.

**Per continuar a la [Part 1: Hello World](./01_hello_world.md), feu clic a la fletxa a la cantonada inferior dreta d'aquesta pàgina.**
