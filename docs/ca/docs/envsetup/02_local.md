# Instal·lació manual

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Traducció assistida per IA - [més informació i suggeriments](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

És possible instal·lar tot el que necessites per executar la formació al teu propi entorn local de manera manual.

Aquí hem documentat com fer-ho en sistemes estàndard compatibles amb POSIX (assumint una màquina personal com ara un portàtil).
Tingues en compte que alguns detalls poden ser diferents segons el teu sistema específic.

!!! tip "Consell"

    Abans de continuar, has considerat utilitzar l'[enfocament de Devcontainers](03_devcontainer.md)?
    Proporciona totes les eines i dependències necessàries sense requerir instal·lació manual.

## Requisits generals de programari

Nextflow es pot utilitzar en qualsevol sistema compatible amb POSIX (Linux, macOS, Windows Subsystem for Linux, etc.) amb Java instal·lat.
Els nostres cursos de formació tenen alguns requisits addicionals.

En total, necessitaràs tenir el següent programari instal·lat:

- Bash o shell equivalent
- [Java 11 (o posterior, fins a 21)](https://www.oracle.com/technetwork/java/javase/downloads/index.html)
- [Git](https://git-scm.com/)
- [Docker](https://docs.docker.com/get-docker/)
- [Conda](https://conda.io/) 4.5 (o posterior)
- [VSCode](https://code.visualstudio.com) amb l'[extensió de Nextflow](https://www.nextflow.io/docs/latest/developer-env.html#devenv-nextflow)

L'aplicació VSCode és tècnicament opcional però recomanem fermament que la utilitzis per treballar els cursos així com per al teu treball de desenvolupament amb Nextflow en general.

El manual de documentació de Nextflow proporciona instruccions per instal·lar aquestes dependències a [Environment setup](https://www.nextflow.io/docs/latest/developer-env.html).

## Nextflow i eines nf-core

Necessitaràs instal·lar Nextflow mateix, més les eines nf-core, tal com es detalla als articles enllaçats a continuació:

- [Instal·lació de Nextflow](https://www.nextflow.io/docs/latest/install.html)
- [Eines nf-core](https://nf-co.re/docs/nf-core-tools/installation)

Recomanem utilitzar l'opció d'autoinstal·lació per a Nextflow i l'opció PyPI per a les eines nf-core.

!!! warning "Advertència sobre compatibilitat de versions"

    <!-- Any update to this content needs to be copied to the home page -->
    **A partir de gener de 2026, tots els nostres cursos de formació de Nextflow requereixen la versió 25.10.2 de Nextflow o posterior, amb la sintaxi estricta v2 activada, llevat que s'indiqui el contrari.**

    Per a més informació sobre els requisits de versió i la sintaxi estricta v2, consulteu la guia de [versions de Nextflow](../info/nxf_versions.md).

    Les versions antigues del material de formació corresponents a la sintaxi anterior estan disponibles mitjançant el selector de versions a la barra de menú d'aquesta pàgina web.

## Materials de formació

La manera més fàcil de descarregar els materials de formació és clonar tot el repositori utilitzant aquesta comanda:

```bash
git clone https://github.com/nextflow-io/training.git
```

Cada curs té el seu propi directori.
Per treballar un curs, obre una finestra de terminal (idealment, des de dins de l'aplicació VSCode) i fes `cd` al directori corresponent.

Després pots seguir les instruccions del curs proporcionades al lloc web.
