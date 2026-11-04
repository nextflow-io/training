# Primers passos

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Traducció assistida per IA - [més informació i suggeriments](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

## Iniciar un entorn de formació

Per utilitzar l'entorn preconfigurat que proporcionem a GitHub Codespaces, feu clic al botó "Open in GitHub Codespaces" que hi ha a continuació. Per a altres opcions, consulteu [Opcions d'entorn](../../envsetup/index.md).

Recomanem obrir l'entorn de formació en una pestanya o finestra nova del navegador (utilitzeu clic dret, ctrl+clic o cmd+clic segons el vostre equip) perquè pugueu continuar llegint mentre es carrega l'entorn.
Haureu de mantenir aquestes instruccions obertes en paral·lel per treballar durant el curs.

[![Open in GitHub Codespaces](https://github.com/codespaces/badge.svg)](https://codespaces.new/nextflow-io/training?quickstart=1&ref=master)

### Conceptes bàsics de l'entorn

Aquest entorn de formació conté tot el programari, codi i dades necessaris per treballar durant el curs de formació, així que no cal que instal·leu res vosaltres mateixos.

El codespace està configurat amb una interfície VSCode, que inclou un explorador de fitxers, un editor de codi i un terminal de línia de comandes.
Totes les instruccions donades durant el curs (per exemple, "obriu el fitxer", "editeu el codi" o "executeu aquesta comanda") es refereixen a aquestes tres parts de la interfície VSCode tret que s'especifiqui el contrari.

Si esteu treballant en aquest curs per compte propi, familiaritzeu-vos amb els [conceptes bàsics de l'entorn](../../envsetup/01_setup.md) per a més detalls.

### Requisits de versió

Aquesta formació està dissenyada per a Nextflow 25.10.2 o posterior **amb l'analitzador de sintaxi v2 ACTIVAT**.
Si utilitzeu un entorn local o personalitzat, assegureu-vos que utilitzeu la configuració correcta tal com es documenta [aquí](../../info/nxf_versions.md).

## Prepareu-vos per treballar

Un cop el vostre codespace estigui en funcionament, hi ha dues coses que heu de fer abans d'endinsar-vos en la formació: establir el vostre directori de treball per a aquest curs específic i donar un cop d'ull als materials proporcionats.

### Establir el directori de treball

Per defecte, el codespace s'obre amb el directori de treball establert a l'arrel de tots els cursos de formació, però per a aquest curs, treballarem al directori `nf4-science/{DOMAIN_DIR}/`.

Canvieu de directori ara executant aquesta comanda al terminal:

```bash
cd nf4-science/{DOMAIN_DIR}/
```

Podeu configurar VSCode perquè se centri en aquest directori, de manera que només els fitxers rellevants es mostrin a la barra lateral de l'explorador de fitxers:

```bash
code .
```

!!! tip "Consell"

    Si per qualsevol motiu sortiu d'aquest directori (per exemple, el vostre codespace s'adorm), sempre podeu utilitzar la ruta completa per tornar-hi, assumint que esteu executant això dins de l'entorn de formació de Github Codespaces:

    ```bash
    cd /workspaces/training/nf4-science/{DOMAIN_DIR}
    ```

Ara donem un cop d'ull als continguts.

### Explorar els materials proporcionats

Podeu explorar els continguts d'aquest directori utilitzant l'explorador de fitxers a la part esquerra de l'espai de treball de formació.
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
    │   ├── {DATA_SUBDIRS}
    │   └── samplesheet.csv
    ├── {DOMAIN_DIR}.nf
    ├── modules
    │   ├── {TOOL_A_MODULE}.nf
    │   └── {TOOL_B_MODULE}.nf
    ├── nextflow.config
    └── solutions
        ├── modules
        ├── part2
        └── part3

    N directories, N files
    ```

Feu clic a la caixa de color per expandir la secció i veure els seus continguts.
Utilitzem seccions desplegables com aquesta per mostrar la sortida esperada de les comandes així com els continguts de directoris i fitxers de manera concisa.

- **El fitxer `{DOMAIN_DIR}.nf`** és un script de workflow que construireu durant el curs.

- **El directori `modules`** conté fitxers de mòdul esquelet que omplireu durant el curs.

- **El fitxer `nextflow.config`** és un fitxer de configuració que estableix propietats mínimes de l'entorn.
  Podeu ignorar-lo de moment.

- **El directori `data`** conté dades d'entrada i recursos relacionats, descrits més endavant en el curs.

- **El directori `solutions`** conté fitxers de mòdul completats i solucions específiques de cada part que poden servir com a punt de partida per a la següent part.
  Estan pensats per ser utilitzats com a referència per comprovar el vostre treball i resoldre qualsevol problema.

## Llista de verificació de preparació

Creieu que esteu preparats per començar?

- [ ] Entenc l'objectiu d'aquest curs i els seus prerequisits
- [ ] El meu entorn està en funcionament
- [ ] He establert el meu directori de treball adequadament

Si podeu marcar totes les caselles, esteu preparats per començar.

**Per continuar a [Part 1: Visió general del mètode](./01_method.md), feu clic a la fletxa de la cantonada inferior dreta d'aquesta pàgina.**
