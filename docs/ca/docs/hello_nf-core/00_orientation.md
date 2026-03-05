# Primers passos

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Traducció assistida per IA - [més informació i suggeriments](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

## Iniciar un entorn de formació

Per utilitzar l'entorn preconfigurat que proporcionem a GitHub Codespaces, feu clic al botó "Open in GitHub Codespaces" que apareix a continuació. Per a altres opcions, consulteu [Opcions d'entorn](../envsetup/index.md).

Recomanem obrir l'entorn de formació en una pestanya o finestra nova del navegador (utilitzeu clic dret, ctrl+clic o cmd+clic segons el vostre equip) perquè pugueu continuar llegint mentre es carrega l'entorn.
Haureu de mantenir aquestes instruccions obertes en paral·lel per treballar durant el curs.

[![Open in GitHub Codespaces](https://github.com/codespaces/badge.svg)](https://codespaces.new/nextflow-io/training?quickstart=1&ref=master)

### Conceptes bàsics de l'entorn

Aquest entorn de formació conté tot el programari, codi i dades necessaris per treballar durant el curs de formació, de manera que no cal que instal·leu res vosaltres mateixos.

El codespace està configurat amb una interfície VSCode, que inclou un explorador del sistema de fitxers, un editor de codi i un terminal de línia de comandes.
Totes les instruccions donades durant el curs (per exemple, "obriu el fitxer", "editeu el codi" o "executeu aquesta comanda") es refereixen a aquestes tres parts de la interfície VSCode tret que s'especifiqui el contrari.

Si esteu treballant en aquest curs pel vostre compte, familiaritzeu-vos amb els [conceptes bàsics de l'entorn](../envsetup/01_setup.md) per a més detalls.

### Requisits de versió

Aquesta formació està dissenyada per a **Nextflow 25.10.2** o posterior **amb l'analitzador de sintaxi v2 DESACTIVAT**.

#### Si utilitzeu el nostre entorn de formació:

HEU d'executar la comanda següent abans de continuar:

```bash
export NXF_SYNTAX_PARSER=v1
```

#### Si utilitzeu un entorn local o personalitzat:

Assegureu-vos que utilitzeu la configuració correcta tal com es documenta [aquí](../info/nxf_versions.md).

La formació també requereix **nf-core tools 3.4.1**.
Si utilitzeu una versió diferent de les eines nf-core, podeu tenir dificultats per seguir el curs.

Podeu comprovar quina versió està instal·lada al vostre entorn utilitzant la comanda `nf-core --version`.

## Prepareu-vos per treballar

Un cop el vostre codespace estigui en funcionament, hi ha dues coses que heu de fer abans d'endinsar-vos en la formació: establir el vostre directori de treball per a aquest curs específic i donar un cop d'ull als materials proporcionats.

### Establir el directori de treball

Per defecte, el codespace s'obre amb el directori de treball establert a l'arrel de tots els cursos de formació, però per a aquest curs, treballarem al directori `hello-nf-core/`.

Canvieu de directori ara executant aquesta comanda al terminal:

```bash
cd hello-nf-core/
```

!!! tip "Consell"

    Si per qualsevol motiu sortiu d'aquest directori (per exemple, el vostre codespace s'adorm), sempre podeu utilitzar el camí complet per tornar-hi, assumint que esteu executant això dins de l'entorn de formació de Github Codespaces:

    ```bash
    cd /workspaces/training/hello-nf-core
    ```

Ara donem un cop d'ull al contingut d'aquest directori.

### Explorar els materials proporcionats

Podeu explorar el contingut d'aquest directori utilitzant l'explorador de fitxers a la part esquerra de l'espai de treball de formació.
Alternativament, podeu utilitzar la comanda `tree`.

Durant el curs, utilitzem la sortida de `tree` per representar l'estructura i el contingut del directori de forma llegible, de vegades amb modificacions menors per claredat.

Aquí generem una taula de continguts fins al segon nivell:

```bash
tree . -L 2
```

??? abstract "Contingut del directori"

    ```console
    .
    ├── greetings.csv
    ├── original-hello
    │   ├── hello.nf
    │   ├── modules
    │   └── nextflow.config
    └── solutions
        ├── composable-hello
        ├── core-hello-part2
        ├── core-hello-part3
        ├── core-hello-part4
        ├── core-hello-part5
        └── core-hello-start
    ```

Feu clic a la caixa de color per expandir la secció i veure el seu contingut.
Utilitzem seccions desplegables com aquesta per incloure la sortida esperada de les comandes de manera concisa.

- **El fitxer `greetings.csv`** és un CSV que conté algunes dades columnars mínimes que utilitzem per a proves.

- **El directori `original-hello`** conté una còpia del codi font produït en treballar durant la sèrie completa de formació Hello Nextflow (amb Docker activat).

- **El directori `solutions`** conté els scripts de workflow completats que resulten de cada pas del curs.
  Estan pensats per ser utilitzats com a referència per comprovar el vostre treball i resoldre qualsevol problema.

## Llista de verificació de preparació

Creieu que esteu preparats per començar?

- [ ] Entenc l'objectiu d'aquest curs i els seus prerequisits
- [ ] El meu entorn està en funcionament
- [ ] M'he assegurat que l'analitzador de sintaxi està establert a **v1**
- [ ] He establert el meu directori de treball adequadament

Si podeu marcar totes les caselles, esteu a punt per començar.

**Per continuar a la Part 1, feu clic a la fletxa de la cantonada inferior dreta d'aquesta pàgina.**
