---
title: Versions de Nextflow
description: Comprendre i gestionar l'evolució de les versions de sintaxi de Nextflow
hide:
  - toc
  - footer
---

## Versió de sintaxi de Nextflow i requisits actuals

A partir de la versió 3.0 del portal de formació, tots els nostres cursos de formació es basen en la versió 25.10.2 de Nextflow, tret que s'especifiqui el contrari a la pàgina d'índex del curs (excepte materials obsolets o arxivats que poden no incloure un avís de versió).

Com que els cursos ara utilitzen entrades tipades a nivell de workflow així com directives de sortida a nivell de workflow, requereixen l'ús de l'analitzador de sintaxi V2.
Si teniu previst utilitzar l'entorn que proporcionem a través de [Github Codespaces](../envsetup/01_setup.md) o [devcontainers locals](../envsetup/03_devcontainer.md), no cal que feu res tret que s'indiqui específicament a les instruccions del curs.
No obstant això, si teniu previst treballar amb les formacions al vostre propi entorn ([Instal·lació manual](../envsetup/02_local.md)), haureu d'assegurar-vos d'utilitzar Nextflow versió 25.10.2 o posterior amb l'analitzador de sintaxi v2 habilitat.

## Versions anteriors dels materials de formació

Els nostres materials de formació tenen versions des de febrer de 2025.

Podeu accedir a versions anteriors dels materials de formació que funcionen amb versions de Nextflow **anteriors a la 25.10.2** mitjançant el menú desplegable a la part superior de cada pàgina que mostra la versió numerada dels materials de formació.
Quan seleccioneu una versió anterior dels materials de formació, els enllaços a l'entorn de formació especificaran automàticament la versió corresponent de l'entorn.

## Altra informació sobre les versions de sintaxi de Nextflow

Nextflow té dos conceptes de versionat diferents que de vegades es confonen: **versions DSL** i **versions de l'analitzador de sintaxi**.

**DSL1 vs DSL2** es refereix a maneres fonamentalment diferents d'escriure pipelines de Nextflow.
DSL1 era la sintaxi original on els processos es connectaven implícitament a través de canals.
DSL2, introduït a Nextflow 20.07, va afegir funcionalitats de modularitat: la capacitat d'importar processos i workflows des d'altres fitxers, blocs `workflow` explícits i sortides de procés amb nom.
DSL1 va quedar obsolet a Nextflow 22.03 i es va eliminar a la versió 22.12.
Tot el codi modern de Nextflow utilitza DSL2.

**Analitzador de sintaxi v1 vs v2** es refereix a diferents analitzadors que tots dos funcionen amb codi DSL2.
L'analitzador v1 és l'original, més permissiu.
L'analitzador v2 és més estricte i habilita noves funcionalitats del llenguatge com el tipat estàtic (entrades i sortides tipades) i directives de sortida a nivell de workflow.
L'analitzador v2 també proporciona millors missatges d'error i detecta més errors en temps d'anàlisi en lloc de temps d'execució.
L'analitzador v2 es convertirà en el predeterminat a Nextflow 26.04.

En resum: DSL2 és el llenguatge que escriviu; la versió de l'analitzador de sintaxi determina com d'estrictament s'interpreta aquest llenguatge i quines funcionalitats avançades estan disponibles.

### Comprovar i establir la versió de Nextflow

Podeu comprovar quina versió de Nextflow està instal·lada al vostre sistema utilitzant la comanda `nextflow --version`.

Per a més informació sobre com actualitzar la vostra versió de Nextflow, consulteu la documentació de referència sobre [Actualitzar Nextflow](https://www.nextflow.io/docs/latest/updating-nextflow.html).

### Habilitar l'analitzador de sintaxi v2

Per **habilitar** l'analitzador de sintaxi v2 per a la vostra sessió actual, executeu la següent comanda al vostre terminal:

```bash
export NXF_SYNTAX_PARSER=v2
```

Per fer-ho permanent (pendent que v2 es converteixi en el predeterminat a Nextflow 26.04), afegiu la comanda export al vostre perfil de shell (`~/.bashrc`, `~/.zshrc`, etc.):

```bash
echo 'export NXF_SYNTAX_PARSER=v2' >> ~/.bashrc
source ~/.bashrc
```

Tingueu en compte que la variable d'entorn `NXF_SYNTAX_PARSER=v2` és un requisit temporal.
A partir de Nextflow 26.04, l'analitzador v2 es convertirà en el predeterminat i aquesta configuració ja no serà necessària.

### Deshabilitar l'analitzador de sintaxi v2

Per **deshabilitar** l'analitzador de sintaxi v2 per a la vostra sessió actual, executeu la següent comanda al vostre terminal:

```bash
export NXF_SYNTAX_PARSER=v1
```

<!-- Will it be possible to disable it in versions after 26.04? -->

### Migrar codi existent

Per a orientació sobre la migració de codi existent per complir amb versions més recents de Nextflow, consulteu les [Notes de migració](https://www.nextflow.io/docs/latest/migrations/index.html) a la documentació de referència.

Aquests dos articles són especialment útils per migrar a la versió més recent:

- [Migrar a sortides de workflow](https://www.nextflow.io/docs/latest/tutorials/workflow-outputs.html)
- [Migrar a tipus estàtics](https://www.nextflow.io/docs/latest/tutorials/static-types.html)

Ambdues funcionalitats es cobreixen com a part de la formació per a principiants a partir de la versió 3.0 dels materials de formació.

Depenent de la generació de codi Nextflow que tingueu intenció de migrar, és possible que pugueu fer la major part utilitzant el linter de Nextflow amb la comanda `nextflow lint -format`.
Consulteu la referència CLI per a [`lint`](https://www.nextflow.io/docs/latest/reference/cli.html#lint) per a més detalls.

Esperem que això us sigui útil.
Si necessiteu ajuda, contacteu-nos a Slack o al fòrum.

---

<div markdown class="homepage_logos">

![Seqera](../assets/img/seqera_logo.png#only-light)

![Seqera](../assets/img/seqera_logo_dark.png#only-dark)

</div>
