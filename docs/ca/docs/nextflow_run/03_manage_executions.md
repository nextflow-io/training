# Part 3: Gestionar les execucions del workflow

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Traducció assistida per IA - [més informació i suggeriments](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

A mesura que executeu i torneu a executar pipelines, acumuleu historial d'execucions i directoris `work/` antics.
A la [Part 1](./01_run_nextflow.md#23-use-resume-to-skip-completed-work) ja heu utilitzat `-resume` per ometre la feina que ja estava feta.
Aquí aprendreu a generar informes sobre una execució, inspeccionar l'historial d'execucions passades amb [`nextflow log`](https://nextflow.io/docs/latest/reference/cli.html#log), i eliminar directoris de treball antics que ja no necessiteu amb [`nextflow clean`](https://nextflow.io/docs/latest/reference/cli.html#clean).

---

## 1. Generar informes del pipeline

Nextflow pot generar diversos tipus d'informes sobre una execució, cadascun afegit amb el seu propi indicador `-with-*`: un informe d'execució (`-with-report`), una línia de temps d'execució (`-with-timeline`), un fitxer de traça de tasques (`-with-trace`), i un diagrama del workflow (`-with-dag`).
Aquí generarem els dos primers; consulteu [Execution reports](https://nextflow.io/docs/latest/reports.html) a la referència de Nextflow per als altres.

### 1.1. Generar un informe d'execució

Afegiu `-with-report` a qualsevol comanda `nextflow run` per generar un informe HTML un cop el pipeline hagi finalitzat:

```bash
nextflow run main.nf --input data/greetings.csv --character turkey -with-report
```

??? success "Sortida de la comanda"

    ```console hl_lines="6 7 8 9"
    N E X T F L O W   ~  version 26.04.4

    Launching `main.nf` [intergalactic_dalembert] revision: ce74f81996

    executor >  local (8)
    [34/23f10f] sayHello (2)       | 3 of 3 ✔
    [af/cab69d] convertToUpper (3) | 3 of 3 ✔
    [9e/d73afb] collectGreetings   | 1 of 1 ✔
    [3c/392db0] cowpy              | 1 of 1 ✔

    Outputs:

      /workspaces/training/nextflow-run/results

      first_output:
        - batch/intermediates/Bonjour-output.txt
        - batch/intermediates/Hello-output.txt
        - batch/intermediates/Hola-output.txt

      uppercased:
        - batch/intermediates/UPPER-Bonjour-output.txt
        - batch/intermediates/UPPER-Hello-output.txt
        - batch/intermediates/UPPER-Hola-output.txt

      collected: batch/intermediates/COLLECTED-batch-output.txt

      batch_report: batch/batch-report.txt

      cowpy_art: batch/cowpy-COLLECTED-batch-output.txt
    ```

Nextflow escriu l'informe en un fitxer anomenat `report-<timestamp>.html` al directori de treball.
Obriu-lo en un navegador per veure un resum de l'execució, una taula de cada tasca amb el seu estat i temps d'execució, i gràfics d'ús de recursos desglossat per procés.

La pestanya **Tasks** llista totes les tasques que ha executat el pipeline, amb el nom del procés, l'estat i l'ús de recursos:

![Taula de tasques de l'informe d'execució](img/execution_report_tasks.png)

L'informe és especialment útil quan un pipeline tarda més del previst o una tasca falla: la taula de tasques mostra exactament on s'ha invertit el temps i quines tasques han tingut èxit o han fallat.

### 1.2. Generar una línia de temps d'execució

Afegiu `-with-timeline` a una execució per obtenir una vista estil diagrama de Gantt de quan s'ha executat cada tasca:

```bash
nextflow run main.nf --input data/greetings.csv --character turkey -with-timeline
```

??? success "Sortida de la comanda"

    ```console hl_lines="6 7 8 9"
    N E X T F L O W   ~  version 26.04.4

    Launching `main.nf` [jolly_noyce] revision: ce74f81996

    executor >  local (8)
    [ad/e92ef3] sayHello (3)       | 3 of 3 ✔
    [2a/df8a8d] convertToUpper (2) | 3 of 3 ✔
    [be/7fb72a] collectGreetings   | 1 of 1 ✔
    [63/dc9bd6] cowpy              | 1 of 1 ✔

    Outputs:

      /workspaces/training/nextflow-run/results

      first_output:
        - batch/intermediates/Bonjour-output.txt
        - batch/intermediates/Hello-output.txt
        - batch/intermediates/Hola-output.txt

      uppercased:
        - batch/intermediates/UPPER-Bonjour-output.txt
        - batch/intermediates/UPPER-Hello-output.txt
        - batch/intermediates/UPPER-Hola-output.txt

      collected: batch/intermediates/COLLECTED-batch-output.txt

      batch_report: batch/batch-report.txt

      cowpy_art: batch/cowpy-COLLECTED-batch-output.txt
    ```

Nextflow escriu la línia de temps en un fitxer anomenat `timeline-<timestamp>.html`.
Obriu-lo en un navegador per veure una barra per a cada tasca, posicionada i dimensionada segons quan s'ha executat i quant de temps ha trigat:

![Línia de temps d'execució](img/execution_timeline.png)

La línia de temps fa visible d'un cop d'ull la forma d'expansió i contracció de la [Part 1](./01_run_nextflow.md#31-run-the-workflow): les tres tasques `sayHello` s'executen en paral·lel, després les tres tasques `convertToUpper`, i finalment `collectGreetings` i `cowpy` s'executen una darrere l'altra, ja que cadascuna depèn de tot el que hi ha abans.

### Conclusio

Sabeu com generar un informe d'execució HTML amb `-with-report` i una línia de temps d'execució amb `-with-timeline`, i on trobar els altres tipus d'informes que admet Nextflow.

### Què segueix?

Apreneu a inspeccionar l'historial d'execucions passades.

---

## 2. Inspeccionar el registre d'execucions passades

Tant si esteu desenvolupant un pipeline com si l'executeu en producció, en algun moment haureu de consultar informació sobre execucions passades.

### 2.1. El fitxer d'historial

Cada vegada que llanceu un workflow de Nextflow, s'escriu una línia en un fitxer de registre anomenat `history`, dins d'un directori ocult anomenat `.nextflow` al directori de treball actual.

??? abstract "Contingut del fitxer"

    ```txt title=".nextflow/history" linenums="1"
    2026-09-13 01:47:58	2.8s	determined_ramanujan	OK	ce74f81996b8ce999853fdbb25ba7969	03078950-5011-414f-992b-8fe1bd219ac2	nextflow run main.nf --input data/greetings.csv --character turkey
    2026-09-13 01:48:03	3s	prickly_cuvier	OK	ce74f81996b8ce999853fdbb25ba7969	5a9e8bf7-a7db-4c2c-b1d8-3bd4c7266528	nextflow run main.nf --input data/greetings.csv --character tux
    2026-09-13 01:48:09	2.3s	elegant_panini	OK	ce74f81996b8ce999853fdbb25ba7969	517ccc29-9db7-4d8d-8ea0-c8db6e9ac653	nextflow run main.nf --input data/greetings.csv --character stegosaurus
    2026-09-13 01:48:14	1.5s	deadly_lamport	OK	ce74f81996b8ce999853fdbb25ba7969	517ccc29-9db7-4d8d-8ea0-c8db6e9ac653	nextflow run main.nf --input data/greetings.csv --character stegosaurus -resume
    ```

Cada línia us proporciona la marca de temps, la durada, el nom de l'execució, l'estat, l'ID de revisió, l'ID de sessió i la línia de comanda completa d'una execució llançada des d'aquest directori.

Fixeu-vos en les dues últimes línies: són dues invocacions separades (una normal i una amb `-resume`) de la mateixa comanda exacta, i comparteixen el mateix ID de sessió.
L'ID de sessió només canvia quan llanceu una execució genuïnament nova; l'ús de `-resume` el manté, que és com Nextflow sap quina memòria cau reutilitzar.

### 2.2. Utilitzar `nextflow log` per a una vista més clara

Llegir el fitxer d'historial directament funciona, però `nextflow log` formata la mateixa informació amb una capçalera:

```bash
nextflow log
```

??? success "Sortida de la comanda"

    ```console linenums="1"
    TIMESTAMP          	DURATION	RUN NAME            	STATUS	REVISION ID	SESSION ID                          	COMMAND
    2026-09-13 01:47:58	2.8s    	determined_ramanujan	OK    	ce74f81996 	03078950-5011-414f-992b-8fe1bd219ac2	nextflow run main.nf --input data/greetings.csv --character turkey
    2026-09-13 01:48:03	3s      	prickly_cuvier      	OK    	ce74f81996 	5a9e8bf7-a7db-4c2c-b1d8-3bd4c7266528	nextflow run main.nf --input data/greetings.csv --character tux
    2026-09-13 01:48:09	2.3s    	elegant_panini      	OK    	ce74f81996 	517ccc29-9db7-4d8d-8ea0-c8db6e9ac653	nextflow run main.nf --input data/greetings.csv --character stegosaurus
    2026-09-13 01:48:14	1.5s    	deadly_lamport      	OK    	ce74f81996 	517ccc29-9db7-4d8d-8ea0-c8db6e9ac653	nextflow run main.nf --input data/greetings.csv --character stegosaurus -resume
    ```

Nextflow agrupa la informació de memòria cau que utilitza per a `-resume` sota `.nextflow/cache`, indexada per ID de sessió.
Per això, cercar el nom d'execució o l'ID de sessió correctes aquí és el primer pas sempre que necessiteu investigar o netejar una execució passada.

### Conclusio

Sabeu on Nextflow registra l'historial d'execucions passades i com inspeccionar-lo amb `nextflow log`.

### Què segueix?

Apreneu a eliminar directoris de treball antics que ja no necessiteu.

---

## 3. Eliminar directoris de treball antics

Cada execució deixa els seus directoris de tasques sota `work/`, fins i tot després d'haver copiat les sortides que us interessen a `results/`.
Executeu prou pipelines durant el desenvolupament i aquests subdirectoris s'acumulen, de manera que Nextflow proporciona `nextflow clean` per eliminar els que ja no necessiteu.

### 3.1. Determinar els criteris d'eliminació

`nextflow clean` admet diverses maneres de seleccionar què eliminar; consulteu la [documentació de referència](https://www.nextflow.io/docs/latest/reference/cli.html#clean) per a la llista completa.
Aquí eliminareu tot el que prové d'execucions anteriors a una execució determinada, utilitzant el seu nom d'execució.

Cerqueu l'execució més recent que voleu conservar amb `nextflow log`; en l'[exemple de la secció 2.2](#22-use-nextflow-log-for-a-friendlier-view) és `elegant_panini`, l'última execució normal abans de la que utilitza `-resume`.
El nom d'execució és la cadena de dues parts generada automàticament que apareix a la línia de consola `Launching (...)`, o a la columna `RUN NAME` de `nextflow log`.

### 3.2. Fer una execució de prova

Afegiu `-n` primer per comprovar què eliminaria una comanda determinada sense eliminar res realment:

```bash
nextflow clean -before elegant_panini -n
```

??? success "Sortida de la comanda"

    ```console
    Would remove /workspaces/training/nextflow-run/work/e5/9ec3fe24deb5d1ffa478dff5e5e1d4
    Would remove /workspaces/training/nextflow-run/work/75/36b1ece86b213c331852d0a3dc659b
    Would remove /workspaces/training/nextflow-run/work/be/8157808f28699259b5ff37e4bf9f72
    Would remove /workspaces/training/nextflow-run/work/f3/d828b7ba6b315eb08e0f96168119f5
    Would remove /workspaces/training/nextflow-run/work/94/7a2d7ef67d81b6e66430f300816698
    Would remove /workspaces/training/nextflow-run/work/eb/e37c5c8638c45916de3d9794d6d22f
    Would remove /workspaces/training/nextflow-run/work/3f/f56a705a6370f6c615aff22d3f240a
    Would remove /workspaces/training/nextflow-run/work/25/4922fa35ba2c6087777ddce20e9f2f
    Would remove /workspaces/training/nextflow-run/work/a4/2de7bcc2f8fc974261f7280c1480f7
    Would remove /workspaces/training/nextflow-run/work/48/8c49bded446d4fcb78c8125318c5e3
    Would remove /workspaces/training/nextflow-run/work/8e/c9834a9bdd7c21de2991d959201a9d
    Would remove /workspaces/training/nextflow-run/work/d5/a9055dc8c817c6c1f436494685955c
    Would remove /workspaces/training/nextflow-run/work/27/9473df9e55e35c8869d9571b480926
    Would remove /workspaces/training/nextflow-run/work/d3/ebc066e1d844232b9f5ee0a7d87464
    Would remove /workspaces/training/nextflow-run/work/91/17ef1815ed06b8f177d0135d8ef0dc
    Would remove /workspaces/training/nextflow-run/work/12/0443f9fcff27a552d4bc73aaedf24a
    ```

Són 16 directoris de tasques: les 8 tasques de l'execució `turkey` més les 8 de l'execució `tux`, exactament les que esperaríeu per a dues execucions completes d'aquest pipeline de quatre processos.
L'execució `elegant_panini` en si mateixa, i les tasques en memòria cau que l'execució amb `-resume` ha reutilitzat d'ella, es deixen intactes.

La vostra sortida llistarà noms de directoris diferents, i el nombre de línies que obteniu depèn de quantes execucions hàgiu fet; si no veieu cap línia, o bé el nom d'execució no coincideix amb cap del vostre registre, o bé no hi ha res a eliminar abans d'ell.

### 3.3. Procedir amb l'eliminació

Un cop la prova sembli correcta, torneu a executar la mateixa comanda amb `-f` en lloc de `-n`:

```bash
nextflow clean -before elegant_panini -f
```

??? success "Sortida de la comanda"

    ```console
    Removed /workspaces/training/nextflow-run/work/e5/9ec3fe24deb5d1ffa478dff5e5e1d4
    Removed /workspaces/training/nextflow-run/work/75/36b1ece86b213c331852d0a3dc659b
    Removed /workspaces/training/nextflow-run/work/be/8157808f28699259b5ff37e4bf9f72
    Removed /workspaces/training/nextflow-run/work/f3/d828b7ba6b315eb08e0f96168119f5
    Removed /workspaces/training/nextflow-run/work/94/7a2d7ef67d81b6e66430f300816698
    Removed /workspaces/training/nextflow-run/work/eb/e37c5c8638c45916de3d9794d6d22f
    Removed /workspaces/training/nextflow-run/work/3f/f56a705a6370f6c615aff22d3f240a
    Removed /workspaces/training/nextflow-run/work/25/4922fa35ba2c6087777ddce20e9f2f
    Removed /workspaces/training/nextflow-run/work/a4/2de7bcc2f8fc974261f7280c1480f7
    Removed /workspaces/training/nextflow-run/work/48/8c49bded446d4fcb78c8125318c5e3
    Removed /workspaces/training/nextflow-run/work/8e/c9834a9bdd7c21de2991d959201a9d
    Removed /workspaces/training/nextflow-run/work/d5/a9055dc8c817c6c1f436494685955c
    Removed /workspaces/training/nextflow-run/work/27/9473df9e55e35c8869d9571b480926
    Removed /workspaces/training/nextflow-run/work/d3/ebc066e1d844232b9f5ee0a7d87464
    Removed /workspaces/training/nextflow-run/work/91/17ef1815ed06b8f177d0135d8ef0dc
    Removed /workspaces/training/nextflow-run/work/12/0443f9fcff27a552d4bc73aaedf24a
    ```

`nextflow clean` buida els directoris de tasques però deixa els directoris pare de dos caràcters (com `e5/`) al seu lloc.

!!! warning "Advertència"

    Eliminar directoris de treball d'execucions passades els elimina de la memòria cau de Nextflow i suprimeix qualsevol sortida emmagatzemada únicament allà.
    Això trenca la capacitat de Nextflow de reprendre l'execució sense tornar a executar els processos corresponents, de manera que només netegeu execucions de les quals esteu segurs que no haureu de reprendre.
    Per això també val la pena publicar tot el que us interessa a `results/` amb `mode 'copy'` en lloc de dependre del directori `work/` o d'un mode de publicació `symlink`.

### Conclusio

Sabeu com eliminar directoris de treball antics amb `nextflow clean`, i per què fer-ho implica perdre la capacitat de reprendre des d'aquestes execucions.

### Què segueix?

Apreneu a executar pipelines directament des de repositoris remots com GitHub a la [Part 4](./04_remote_repositories.md).

---

## Resum

En aquesta part heu après a:

- Generar un informe d'execució HTML amb `-with-report` i una línia de temps d'execució amb `-with-timeline`
- Inspeccionar l'historial d'execucions passades amb `nextflow log`
- Eliminar directoris de treball antics amb `nextflow clean`, i entendre la implicació sobre la capacitat de reprendre que comporta
