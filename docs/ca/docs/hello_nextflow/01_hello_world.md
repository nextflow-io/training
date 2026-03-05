# Part 1: Hello World

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Traducció assistida per IA - [més informació i suggeriments](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

<div class="video-wrapper">
  <iframe width="560" height="315" src="https://www.youtube.com/embed/tOukLxWCHiA?si=y8lAedhEHWaTV4zd&amp;list=PLPZ8WHdZGxmWKozQuzr27jyMGqp9kElVK&amp;cc_load_policy=1" title="YouTube video player" frameborder="0" allow="accelerometer; autoplay; clipboard-write; encrypted-media; gyroscope; picture-in-picture; web-share" referrerpolicy="strict-origin-when-cross-origin" allowfullscreen></iframe>
</div>

/// caption
:fontawesome-brands-youtube:{ .youtube } Consulteu [la llista de reproducció completa](https://youtube.com/playlist?list=PLPZ8WHdZGxmWKozQuzr27jyMGqp9kElVK&si=eF7cLR62goy-lc6n) al canal de YouTube de Nextflow.

:green_book: La transcripció del vídeo està disponible [aquí](./transcripts/01_hello_world.md).
///

En aquesta primera part del curs de formació Hello Nextflow, comencem amb un exemple molt bàsic de Hello World independent del domini, que anirem construint progressivament per demostrar l'ús de la lògica i els components fonamentals de Nextflow.

??? info "Què és un exemple Hello World?"

    Un "Hello World!" és un exemple minimalista que pretén demostrar la sintaxi i l'estructura bàsiques d'un llenguatge de programació o marc de programari.
    L'exemple normalment consisteix a imprimir la frase "Hello, World!" al dispositiu de sortida, com ara la consola o el terminal, o escriure-la en un fitxer.

---

## 0. Escalfament: Executeu un exemple Hello World directament

Demostrem això amb una comanda simple que executem directament al terminal, per mostrar què fa abans d'embolicar-la amb Nextflow.

!!! tip "Consell"

    Recordeu que ara hauríeu d'estar dins del directori `hello-nextflow/` tal com es descriu a la pàgina [Primers passos](00_orientation.md).

### 0.1. Feu que el terminal digui hola

Executeu la comanda següent al vostre terminal.

```bash
echo 'Hello World!'
```

??? success "Sortida de la comanda"

    ```console
    Hello World!
    ```

Això mostra el text 'Hello World' directament al terminal.

### 0.2. Escriviu la sortida a un fitxer

Executar pipelines implica principalment llegir dades de fitxers i escriure resultats a altres fitxers, així que modifiquem la comanda per escriure la sortida de text a un fitxer per fer l'exemple una mica més rellevant.

```bash
echo 'Hello World!' > output.txt
```

??? success "Sortida de la comanda"

    ```console

    ```

Això no mostra res al terminal.

### 0.3. Trobeu la sortida

El text 'Hello World' ara hauria d'estar al fitxer de sortida que hem especificat, anomenat `output.txt`.
Podeu obrir-lo a l'explorador de fitxers o des de la línia de comandes utilitzant la utilitat `cat`, per exemple.

??? abstract "Contingut del fitxer"

    ```console title="output.txt" linenums="1"
    Hello World!
    ```

Això és el que intentarem replicar amb el nostre primer workflow de Nextflow.

### Conclusió

Ara sabeu com executar una comanda simple al terminal que mostra text i, opcionalment, com fer que escrigui la sortida a un fitxer.

### Què segueix?

Descobriu com seria això escrit com un workflow de Nextflow.

---

## 1. Examineu l'script i executeu-lo

Us proporcionem un script de workflow totalment funcional encara que minimalista anomenat `hello-world.nf` que fa el mateix que abans (escriure 'Hello World!') però amb Nextflow.

Per començar, obrim l'script del workflow perquè pugueu tenir una idea de com està estructurat.
Després l'executarem i buscarem les seves sortides.

### 1.1. Examineu el codi

Trobareu l'script `hello-world.nf` al vostre directori actual, que hauria de ser `hello-nextflow`. Obriu-lo al panell de l'editor.

??? full-code "Fitxer de codi complet"

    ```groovy title="hello-world.nf" linenums="1"
    #!/usr/bin/env nextflow

    /*
    * Utilitza echo per imprimir 'Hello World!' a un fitxer
    */
    process sayHello {

        output:
        path 'output.txt'

        script:
        """
        echo 'Hello World!' > output.txt
        """
    }

    workflow {

        main:
        // emet una salutació
        sayHello()
    }
    ```

Un script de workflow de Nextflow normalment inclou una o més definicions de [**process**](https://nextflow.io/docs/latest/process.html) i el [**workflow**](https://nextflow.io/docs/latest/workflow.html) en si, més alguns blocs opcionals (no presents aquí) que introduirem més endavant.

Cada **process** descriu quines operacions ha de realitzar el pas corresponent del pipeline, mentre que el **workflow** descriu la lògica de flux de dades que connecta els diversos passos.

Primer examinarem més de prop el bloc **process**, després mirarem el bloc **workflow**.

#### 1.1.1. La definició de `process`

El primer bloc de codi descriu un **process**.

La definició del procés comença amb la paraula clau `process`, seguida del nom del procés i finalment el cos del procés delimitat per claus.
El cos del procés ha de contenir un bloc script que especifica la comanda a executar, que pot ser qualsevol cosa que poguéssiu executar en un terminal de línia de comandes.

```groovy title="hello-world.nf" linenums="3"
/*
* Utilitza echo per imprimir 'Hello World!' a un fitxer
*/
process sayHello {

    output:
    path 'output.txt'

    script:
    """
    echo 'Hello World!' > output.txt
    """
}
```

Aquí tenim un **process** anomenat `sayHello` que escriu la seva **sortida** a un fitxer anomenat `output.txt`.

<figure class="excalidraw">
--8<-- "docs/en/docs/hello_nextflow/img/hello_world.svg"
</figure>

Aquesta és una definició de procés molt mínima que només conté una definició d'`output` i l'`script` a executar.

La definició d'`output` inclou el qualificador `path`, que indica a Nextflow que això s'ha de gestionar com un camí (inclou tant camins de directori com fitxers).
Un altre qualificador comú és `val`.

És important destacar que la definició de sortida no _determina_ quina sortida es crearà.
Simplement _declara_ quina és la sortida esperada, perquè Nextflow la pugui buscar un cop l'execució estigui completa.
Això és necessari per verificar que la comanda s'ha executat correctament i per passar la sortida a processos posteriors si cal. La sortida produïda que no coincideixi amb el que es declara al bloc de sortida no es passarà als processos posteriors.

!!! warning "Advertència"

    Aquest exemple és fràgil perquè hem codificat el nom del fitxer de sortida en dos llocs separats (l'script i els blocs de sortida).
    Si canviem un però no l'altre, l'script fallarà.
    Més endavant, aprendreu maneres d'utilitzar variables per mitigar aquest problema.

En un pipeline del món real, un procés normalment conté blocs addicionals com directives i entrades, que introduirem d'aquí a poc.

#### 1.1.2. La definició de `workflow`

El segon bloc de codi descriu el **workflow** en si.
La definició del workflow comença amb la paraula clau `workflow`, seguida d'un nom opcional, després el cos del workflow delimitat per claus.

Aquí tenim un **workflow** que consisteix en un bloc `main:` (que diu 'aquest és el cos principal del workflow') que conté una crida al procés `sayHello`.

```groovy title="hello-world.nf" linenums="17"
workflow {

    main:
    // emet una salutació
    sayHello()
}
```

Aquesta és una definició de **workflow** molt mínima.
En un pipeline del món real, el workflow normalment conté múltiples crides a **processos** connectats per **canals**, i els processos esperen una o més **entrades** variables.

Aprendreu com afegir entrades variables més endavant en aquest mòdul de formació; i aprendreu com afegir més processos i connectar-los mitjançant canals a la Part 3 d'aquest curs.

!!! tip "Consell"

    Tècnicament la línia `main:` no és necessària per a workflows simples com aquest, així que podeu trobar workflows que no la tenen.
    Però la necessitarem per aprofitar les sortides a nivell de workflow, així que és millor incloure-la des del principi.

### 1.2. Executeu el workflow

Mirar codi no és tan divertit com executar-lo, així que provem-ho a la pràctica.

#### 1.2.1. Llanceu el workflow i superviseu l'execució

Al terminal, executeu la comanda següent:

```bash
nextflow run hello-world.nf
```

??? success "Sortida de la comanda"

    ```console hl_lines="7"
    N E X T F L O W   ~  version 25.10.2

    Launching `hello-world.nf` [goofy_torvalds] DSL2 - revision: c33d41f479

    executor >  local (1)
    [65/7be2fa] sayHello | 1 of 1 ✔
    ```

Si la sortida de la vostra consola s'assembla a això, felicitats, acabeu d'executar el vostre primer workflow de Nextflow!

La sortida més important aquí és l'última línia, que està ressaltada a la sortida anterior:

```console
[65/7be2fa] sayHello | 1 of 1 ✔
```

Això ens indica que el procés `sayHello` s'ha executat correctament una vegada (`1 of 1 ✔`).

És important destacar que aquesta línia també us indica on trobar la sortida de la crida al procés `sayHello`.
Mirem-ho ara.

#### 1.2.2. Trobeu la sortida i els registres al directori `work`

Quan executeu Nextflow per primera vegada en un directori determinat, crea un directori anomenat `work` on escriurà tots els fitxers (i qualsevol enllaç simbòlic) generats durant l'execució.

Dins del directori `work`, Nextflow organitza les sortides i els registres per crida de procés.
Per a cada crida de procés, Nextflow crea un subdirectori imbricat, anomenat amb un hash per fer-lo únic, on prepararà totes les entrades necessàries (utilitzant enllaços simbòlics per defecte), escriurà fitxers auxiliars i escriurà registres i qualsevol sortida del procés.

El camí a aquest subdirectori es mostra en forma truncada entre claudàtors a la sortida de la consola.
Mirant el que hem obtingut per a l'execució mostrada anteriorment, la línia de registre de la consola per al procés sayHello comença amb `[65/7be2fa]`. Això correspon al camí de directori següent: `work/65/7be2fad5e71e5f49998f795677fd68`

Vegem què hi ha allà.

??? abstract "Contingut del directori"

    ```console
    work
    └── 65
        └── 7be2fad5e71e5f49998f795677fd68
            ├── .command.begin
            ├── .command.err
            ├── .command.log
            ├── .command.out
            ├── .command.run
            ├── .command.sh
            ├── .exitcode
            └── output.txt
    ```

??? question "No veieu el mateix?"

    Els noms exactes dels subdirectoris seran diferents al vostre sistema.

    Si navegueu pel contingut del subdirectori de tasques a l'explorador de fitxers de VSCode, veureu tots els fitxers immediatament.
    No obstant això, els fitxers de registre estan configurats per ser invisibles al terminal, així que si voleu utilitzar `ls` o `tree` per veure'ls, haureu d'establir l'opció rellevant per mostrar fitxers invisibles.

    ```bash
    tree -a work
    ```

El primer que voleu mirar és la sortida real del workflow, és a dir, el fitxer `output.txt` produït pel procés `sayHello`.
Obriu-lo i trobareu la salutació `Hello World!`, que era l'objectiu del nostre workflow minimalista.

??? abstract "Contingut del fitxer"

    ```console title="output.txt"
    Hello World!
    ```

Ha funcionat!

És cert que pot semblar molt codi d'embolcall per a un resultat tan petit, però el valor de tot aquest codi d'embolcall serà més obvi quan comencem a llegir fitxers d'entrada i a encadenar múltiples passos.

Dit això, mirem també els altres fitxers d'aquest directori. Aquests són fitxers auxiliars i de registre produïts per Nextflow com a part de l'execució de la tasca.

- **`.command.begin`**: Metadades relacionades amb l'inici de l'execució de la crida al procés
- **`.command.err`**: Missatges d'error (`stderr`) emesos per la crida al procés
- **`.command.log`**: Sortida de registre completa emesa per la crida al procés
- **`.command.out`**: Sortida regular (`stdout`) de la crida al procés
- **`.command.run`**: Script complet executat per Nextflow per executar la crida al procés
- **`.command.sh`**: La comanda que realment va executar la crida al procés
- **`.exitcode`**: El codi de sortida resultant de la comanda

El fitxer `.command.sh` és especialment útil perquè us indica la comanda principal que Nextflow va executar, sense incloure tota la comptabilitat i la configuració de tasques/entorn.

??? abstract "Contingut del fitxer"

    ```console title=".command.sh"
    #!/bin/bash -ue
    echo 'Hello World!' > output.txt
    ```

Això coincideix amb el que vam executar abans manualment.

En aquest cas és molt senzill perquè la comanda del procés estava codificada, però més endavant al curs veureu comandes de procés que impliquen alguna interpolació de variables.
Això fa especialment valuós poder veure exactament com Nextflow va interpretar el codi i quina comanda es va produir quan esteu solucionant problemes d'una execució fallida.

### 1.3. Executeu el workflow de nou

Proveu de tornar a executar el workflow unes quantes vegades, després mireu els directoris de tasques sota `work/`.

??? abstract "Contingut del directori"

    ```console
    work
    ├── 0f
    │   └── 52b7e07b0e274a80843fca48ed21b8
    │       ├── .command.begin
    │       ├── .command.err
    │       ├── .command.log
    │       ├── .command.out
    │       ├── .command.run
    │       ├── .command.sh
    │       ├── .exitcode
    │       └── output.txt
    ├── 65
        └── 7be2fad5e71e5f49998f795677fd68
    │   │   ├── .command.begin
    │   │   ├── .command.err
    │   │   ├── .command.log
    │   │   ├── .command.out
    │   │   ├── .command.run
    │   │   ├── .command.sh
    │   │   ├── .exitcode
    │   │   └── output.txt
    │   └── e029f2e75305874a9ab263d21ebc2c
    │       ├── .command.begin
    │       ├── .command.err
    │       ├── .command.log
    │       ├── .command.out
    │       ├── .command.run
    │       ├── .command.sh
    │       ├── .exitcode
    │       └── output.txt
    ├── 6c
    │   └── d4fd787e0b01b3c82e85696c297500
    │       ├── .command.begin
    │       ├── .command.err
    │       ├── .command.log
    │       ├── .command.out
    │       ├── .command.run
    │       ├── .command.sh
    │       ├── .exitcode
    │       └── output.txt
    └── e8
        └── ab99fad46ade52905ec973ff39bb80
            ├── .command.begin
            ├── .command.err
            ├── .command.log
            ├── .command.out
            ├── .command.run
            ├── .command.sh
            ├── .exitcode
            └── output.txt
    ```

Veieu que s'ha creat un nou subdirectori amb un conjunt complet de fitxers de sortida i registre per a cada execució.
Això us mostra que executar el mateix workflow diverses vegades no sobreescriurà els resultats d'execucions anteriors.

### Conclusió

Sabeu com desxifrar un script simple de Nextflow, executar-lo i trobar la sortida i els fitxers de registre rellevants al directori work.

### Què segueix?

Apreneu com publicar les sortides del workflow a una ubicació més convenient.

---

## 2. Publiqueu sortides

Com acabeu d'aprendre, la sortida produïda pel nostre pipeline està enterrada en un directori de treball diverses capes de profunditat.
Això es fa a propòsit; Nextflow té el control d'aquest directori i no se suposa que hi interactuem.
No obstant això, això fa que sigui incòmode recuperar sortides que ens importen.

Afortunadament, Nextflow proporciona una manera de publicar sortides a un directori designat utilitzant [definicions de sortida de workflow](https://nextflow.io/docs/latest/workflow.html#workflow-outputs).

### 2.1. Ús bàsic

Això implicarà dos nous fragments de codi:

1. Un bloc `publish:` dins del cos del `workflow`, declarant sortides de procés.
2. Un bloc `output` a l'script especificant opcions de sortida com el mode i la ubicació.

#### 2.1.1. Declareu la sortida del procés `sayHello`

Hem d'afegir un bloc `publish:` al cos del workflow (el mateix tipus d'element de codi que el bloc `main:`) i llistar la sortida del procés `sayHello()`.

Al fitxer d'script del workflow `hello-world.nf`, afegiu les línies de codi següents:

=== "Després"

    ```groovy title="hello-world.nf" linenums="17" hl_lines="7-8"
    workflow {

        main:
        // emet una salutació
        sayHello()

        publish:
        first_output = sayHello.out
    }
    ```

=== "Abans"

    ```groovy title="hello-world.nf" linenums="17"
    workflow {

        main:
        // emet una salutació
        sayHello()
    }
    ```

Veieu que podem referir-nos a la sortida del procés simplement fent `sayHello().out`, i assignar-li un nom arbitrari, `first_output`.

#### 2.1.2. Afegiu un bloc `output:` a l'script

Ara només cal afegir el bloc `output:` on s'especificarà el camí del directori de sortida. Tingueu en compte que aquest nou bloc se situa **fora** i **sota** el bloc `workflow` dins de l'script.

Al fitxer d'script del workflow `hello-world.nf`, afegiu les línies de codi següents:

=== "Després"

    ```groovy title="hello-world.nf" linenums="17" hl_lines="11-15"
    workflow {

        main:
        // emet una salutació
        sayHello()

        publish:
        first_output = sayHello.out
    }

    output {
        first_output {
            path '.'
        }
    }
    ```

=== "Abans"

    ```groovy title="hello-world.nf" linenums="17"
    workflow {

        main:
        // emet una salutació
        sayHello()

        publish:
        first_output = sayHello.out
    }
    ```

Podem utilitzar això per assignar camins específics a qualsevol sortida de procés declarada al bloc `workflow`.
Més endavant, aprendreu maneres de generar estructures de directoris de sortida sofisticades, però de moment, només estem codificant un camí mínim per simplicitat.

#### 2.1.3. Executeu el workflow

Ara executeu l'script del workflow modificat:

```bash
nextflow run hello-world.nf
```

??? success "Sortida de la comanda"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `hello-world.nf` [jovial_mayer] DSL2 - revision: 35bd3425e5

    executor >  local (1)
    [9f/48ef97] sayHello | 1 of 1 ✔
    ```

La sortida del terminal hauria de semblar familiar. Externament, no ha canviat res.

No obstant això, comproveu el vostre explorador de fitxers: aquesta vegada, Nextflow ha creat un nou directori anomenat `results/`.

??? abstract "Contingut del directori"

    ```console hl_lines="10-11 22"
    .
    ├── greetings.csv
    ├── hello-channels.nf
    ├── hello-config.nf
    ├── hello-containers.nf
    ├── hello-modules.nf
    ├── hello-workflow.nf
    ├── hello-world.nf
    ├── nextflow.config
    ├── results
    │   └── output.txt -> /workspaces/training/hello-nextflow/work/9f/48ef97f110b0dbd83635d7cbe288d2/output.txt
    ├── solutions
    │   ├── 1-hello-world
    │   ├── 2-hello-channels
    │   ├── 3-hello-workflow
    │   ├── 4-hello-modules
    │   ├── 5-hello-containers
    │   └── 6-hello-config
    ├── test-params.json
    └── work
        ├── 65
        └── 9f
    ```

Dins del directori `results`, trobem un enllaç simbòlic a l'`output.txt` produït al directori work per la comanda que acabem d'executar.

Això ens permet recuperar fàcilment fitxers de sortida sense haver de buscar pel subdirectori work.

### 2.2. Establiu una ubicació personalitzada

Tenir una ubicació per defecte és genial, però potser voleu personalitzar on es desen els resultats i com s'organitzen.

Per exemple, potser voleu organitzar les vostres sortides en subdirectoris.
La manera més senzilla de fer-ho és assignar un camí de sortida específic per sortida.

#### 2.2.1. Modifiqueu el camí de sortida

Una vegada més, modificar el comportament de publicació per a una sortida específica és realment senzill.
Per establir una ubicació personalitzada, només cal editar el `path` en conseqüència:

=== "Després"

    ```groovy title="hello-world.nf" linenums="27" hl_lines="3"
    output {
        first_output {
            path 'hello_world'
        }
    }
    ```

=== "Abans"

    ```groovy title="hello-world.nf" linenums="27" hl_lines="3"
    output {
        first_output {
            path '.'
        }
    }
    ```

Com que això s'estableix al nivell de la sortida individual, podeu especificar diferents ubicacions i subdirectoris per adaptar-se a les vostres necessitats.

#### 2.2.2. Executeu el workflow de nou

Provem-ho.

```bash
nextflow run hello-world.nf
```

??? success "Sortida de la comanda"

    ```console
     N E X T F L O W   ~  version 25.10.2

    Launching `hello-world.nf` [tiny_shaw] DSL2 - revision: 757723adc1

    executor >  local (1)
    [8c/79499c] process > sayHello [100%] 1 of 1 ✔
    ```

Aquesta vegada el resultat s'escriu sota el subdirectori especificat.

??? abstract "Contingut del directori"

    ```console hl_lines="2-3"
    results/
    ├── hello_world
    │   └── output.txt -> /workspaces/training/hello-nextflow/work/8c/79499c2e506b79e2e01acb808d9d12/output.txt
    └── output.txt -> /workspaces/training/hello-nextflow/work/65/f56f2cd75df1352e106fcdd084b97b/output.txt
    ```

Veieu que el resultat de l'execució anterior encara hi és.

<figure class="excalidraw">
--8<-- "docs/en/docs/hello_nextflow/img/hello_world_output.svg"
</figure>

Podeu utilitzar tants nivells d'imbricació com vulgueu.
També és possible utilitzar el nom del procés o altres variables per anomenar els directoris utilitzats per organitzar resultats, i és possible canviar el nom per defecte del directori de sortida de nivell superior (que està controlat per la bandera CLI `-o` o la variable de configuració `outputDir`).
Cobrirem aquestes opcions més endavant a la formació.

### 2.3. Establiu el mode de publicació a còpia

Per defecte, les sortides es publiquen com a enllaços simbòlics des del directori `work`.
Això significa que només hi ha un únic fitxer al sistema de fitxers.

Això és genial quan esteu tractant amb fitxers molt grans, dels quals no voleu emmagatzemar múltiples còpies.
No obstant això, si suprimiu el directori work en algun moment (cobrirem les operacions de neteja en breu), perdreu l'accés al fitxer.
Així que heu de tenir un pla per desar còpies de qualsevol fitxer important en un lloc segur.

Una opció fàcil és canviar el mode de publicació a còpia per a les sortides que us importen.

#### 2.3.1. Afegiu la directiva mode

Aquesta part és realment senzilla.
Només cal afegir `mode 'copy'` a la definició de sortida a nivell de workflow rellevant:

=== "Després"

    ```groovy title="hello-world.nf" linenums="27" hl_lines="4"
    output {
        first_output {
            path 'hello_world'
            mode 'copy'
        }
    }
    ```

=== "Abans"

    ```groovy title="hello-world.nf" linenums="27"
    output {
        first_output {
            path 'hello_world'
        }
    }
    ```

Això estableix el mode de publicació per a aquesta sortida específica.

#### 2.3.2. Executeu el workflow de nou

Provem-ho.

```bash
nextflow run hello-world.nf
```

??? success "Sortida de la comanda"

    ```console
     N E X T F L O W   ~  version 25.10.2

    Launching `hello-world.nf` [tiny_shaw] DSL2 - revision: 757723adc1

    executor >  local (1)
    [df/521638] process > sayHello [100%] 1 of 1 ✔
    ```

Aquesta vegada, si mireu els resultats, el fitxer és una còpia adequada en lloc de només un enllaç simbòlic.

??? abstract "Contingut del directori"

    ```console hl_lines="3"
    results/
    ├── hello_world
    │   └── output.txt
    └── output.txt -> /workspaces/training/hello-nextflow/work/65/f56f2cd75df1352e106fcdd084b97b/output.txt
    ```

Com que això també s'estableix al nivell de la sortida individual, us permet establir el mode de publicació de manera granular.
Això serà especialment útil més endavant quan passem a pipelines de múltiples passos, on potser voleu copiar només les sortides finals i deixar les sortides intermèdies com a enllaços simbòlics, per exemple.

Com s'ha indicat anteriorment, hi ha altres opcions més sofisticades per controlar com es publiquen les sortides.
Us mostrarem com utilitzar-les a temps en el vostre viatge amb Nextflow.

### 2.4. Nota sobre les directives `publishDir` a nivell de procés

Fins fa molt poc, la manera establerta de publicar sortides era fer-ho al nivell de cada procés individual utilitzant una directiva `publishDir`.

Per aconseguir el que acabem de fer per a les sortides del procés `sayHello`, hauríem afegit en canvi la línia següent a la definició del procés:

```groovy title="hello-world.nf" linenums="6" hl_lines="3"
process sayHello {

    publishDir 'results/hello_world', mode: 'copy'

    output:
    path 'output.txt'

    script:
    """
    echo 'Hello World!' > output.txt
    """
}
```

Encara trobareu aquest patró de codi per tot arreu en pipelines i mòduls de procés de Nextflow més antics, així que és important ser-ne conscient.
No obstant això, no recomanem utilitzar-lo en cap treball nou, ja que eventualment no es permetrà en futures versions del llenguatge Nextflow.

### Conclusió

Sabeu com publicar sortides de workflow a una ubicació més convenient.

### Què segueix?

Apreneu a proporcionar una entrada variable mitjançant un paràmetre de línia de comandes i utilitzar valors per defecte de manera efectiva.

---

## 3. Utilitzeu una entrada variable passada per línia de comandes

En el seu estat actual, el nostre workflow utilitza una salutació codificada a la comanda del procés.
Volem afegir una mica de flexibilitat utilitzant una variable d'entrada, perquè puguem canviar més fàcilment la salutació en temps d'execució.

Això requereix que fem tres conjunts de canvis al nostre script:

1. Canviar el procés per esperar una entrada variable
2. Configurar un paràmetre de línia de comandes per capturar l'entrada de l'usuari
3. Passar l'entrada al procés al cos del workflow

Fem aquests canvis un a la vegada.

### 3.1. Canvieu el procés `sayHello` per esperar una entrada variable

Hem d'editar la definició del procés per (1) acceptar una variable d'entrada i (2) utilitzar aquesta variable a la línia de comandes.

#### 3.1.1. Afegiu un bloc d'entrada a la definició del procés

Primer, adaptem la definició del procés per acceptar una entrada anomenada `greeting`.

Al bloc del procés, feu el canvi de codi següent:

=== "Després"

    ```groovy title="hello-world.nf" linenums="6" hl_lines="3-4"
    process sayHello {

        input:
        val greeting

        output:
        path 'output.txt'
    ```

=== "Abans"

    ```groovy title="hello-world.nf" linenums="6"
    process sayHello {

        output:
        path 'output.txt'
    ```

La variable `greeting` està prefixada per `val` per indicar a Nextflow que és un valor (no un camí).

#### 3.1.2. Editeu la comanda del procés per utilitzar la variable d'entrada

Ara intercanviem el valor codificat original pel valor de la variable d'entrada que esperem rebre.

Al bloc del procés, feu el canvi de codi següent:

=== "Després"

    ```groovy title="hello-world.nf" linenums="14" hl_lines="3"
    script:
    """
    echo '${greeting}' > output.txt
    """
    ```

=== "Abans"

    ```groovy title="hello-world.nf" linenums="14" hl_lines="3"
    script:
    """
    echo 'Hello World!' > output.txt
    """
    ```

El símbol `$` i les claus (`{ }`) indiquen a Nextflow que això és un nom de variable que cal substituir pel valor d'entrada real (=interpolat).

!!! tip "Consell"

    Les claus (`{ }`) eren tècnicament opcionals en versions anteriors de Nextflow, així que podeu veure workflows més antics on això s'escriu com `echo '$greeting' > output.txt`.

Ara que el procés `sayHello()` està preparat per acceptar una entrada variable, necessitem una manera de proporcionar un valor d'entrada a la crida del procés al nivell del workflow.

### 3.2. Configureu un paràmetre de línia de comandes per capturar l'entrada de l'usuari

Podríem simplement codificar una entrada directament fent la crida al procés `sayHello('Hello World!')`.
No obstant això, quan estem fent un treball real amb el nostre workflow, voldrem poder controlar les seves entrades des de la línia de comandes, perquè puguem fer alguna cosa així:

<figure class="excalidraw">
--8<-- "docs/en/docs/hello_nextflow/img/hello_world_input.svg"
</figure>

Afortunadament, Nextflow té un sistema de paràmetres de workflow integrat anomenat [`params`](https://nextflow.io/docs/latest/config.html#params) que facilita declarar i utilitzar paràmetres CLI.

La sintaxi general és declarar `params.<nom_parametre>` per indicar a Nextflow que esperi un paràmetre `--<nom_parametre>` a la línia de comandes.

Aquí, volem crear un paràmetre anomenat `--input`, així que hem de declarar `params.input` en algun lloc del workflow.
En principi podem escriure-ho a qualsevol lloc; però com que el volem donar a la crida del procés `sayHello()`, podem connectar-lo allà directament escrivint `sayHello(params.input)`.

Al bloc del workflow, feu el canvi de codi següent:

=== "Després"

    ```groovy title="hello-world.nf" linenums="23" hl_lines="2"
    // emet una salutació
    sayHello(params.input)
    ```

=== "Abans"

    ```groovy title="hello-world.nf" linenums="23" hl_lines="2"
    // emet una salutació
    sayHello()
    ```

Això indica a Nextflow que executi el procés `sayHello` amb el valor proporcionat mitjançant el paràmetre `--input`.

En efecte, hem aconseguit els passos (2) i (3) descrits al principi de la secció d'una sola vegada.

### 3.3. Executeu la comanda del workflow

Executem-ho!

```bash
nextflow run hello-world.nf --input 'Bonjour le monde!'
```

??? success "Sortida de la comanda"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `hello-world.nf` [elated_lavoisier] DSL2 - revision: 7c031b42ea

    executor >  local (1)
    [4b/654319] sayHello | 1 of 1 ✔
    ```

Si heu fet totes aquestes edicions correctament, hauríeu d'obtenir una altra execució exitosa.

Assegureu-vos d'obrir el fitxer de sortida per comprovar que ara teniu la nova versió de la salutació.

??? abstract "Contingut del fitxer"

    ```console title="results/hello_world/output.txt"
    Bonjour le monde!
    ```

Et voilà!

Tingueu en compte com la nova execució ha sobreescrit el fitxer de sortida publicat al directori `results`.
No obstant això, els resultats de les execucions anteriors encara es conserven als directoris de tasques sota `work`.

!!! tip "Consell"

    Podeu distingir fàcilment els paràmetres a nivell de Nextflow dels paràmetres a nivell de pipeline.

    - Els paràmetres que s'apliquen a un pipeline sempre porten un guió doble (`--`).
    - Els paràmetres que modifiquen una configuració de Nextflow, _p. ex._ la funció `-resume` que vam utilitzar anteriorment, porten un guió simple (`-`).

### 3.4. Utilitzeu valors per defecte per als paràmetres de línia de comandes

D'acord, això era convenient, però en molts casos, té sentit proporcionar un valor per defecte per a un paràmetre determinat perquè no hàgiu d'especificar-lo per a cada execució.

#### 3.4.1. Establiu un valor per defecte per al paràmetre CLI

Donem al paràmetre `input` un valor per defecte declarant-lo abans de la definició del workflow.

```groovy title="hello-world.nf" linenums="20"
/*
 * Paràmetres del pipeline
 */
params {
    input: String = 'Holà mundo!'
}
```

Com veieu, podem especificar el tipus d'entrada que el workflow espera (Nextflow 25.10.2 i posteriors).
La sintaxi és `nom: Tipus = valor_per_defecte`.
Els tipus suportats inclouen `String`, `Integer`, `Float`, `Boolean` i `Path`.

!!! info "Info"

    En workflows més antics, podeu veure que tot aquest bloc `params` s'escriu simplement com `input = 'Holà mundo!'`.

A mesura que afegiu més paràmetres al vostre pipeline, hauríeu d'afegir-los tots a aquest bloc, tant si necessiteu donar-los un valor per defecte com si no.
Això facilitarà trobar tots els paràmetres configurables d'una ullada.

#### 3.4.2. Executeu el workflow de nou sense especificar el paràmetre

Ara que teniu un valor per defecte establert, podeu executar el workflow de nou sense haver d'especificar un valor a la línia de comandes.

```bash
nextflow run hello-world.nf
```

??? success "Sortida de la comanda"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `hello-world.nf` [determined_edison] DSL2 - revision: 3539118582

    executor >  local (1)
    [72/394147] sayHello | 1 of 1 ✔
    ```

??? question "Si no ha funcionat"

    Si això ha fallat amb un error que s'assembla a això:

    ```
    ERROR ~ Script compilation error
    - file : /workspaces/training/hello-nextflow/solutions/1-hello-world/hello-world-3.nf
    - cause: you tried to assign a value to the class 'java.lang.String'
    @ line 24, column 12.
          input: String = 'Holà mundo!'
                  ^

    1 error


    -- Check '.nextflow.log' file for details
    ```

    Llavors probablement esteu utilitzant l'analitzador de llenguatge Nextflow v1 més antic.
    Això es va esmentar al principi del curs, però potser us ho vau perdre.
    Consulteu el material d'ajuda sobre [versions de Nextflow](../info/nxf_versions.md).

    En resum, si esteu utilitzant Nextflow `25.10` llavors heu d'habilitar l'analitzador de llenguatge v2:

    ```bash
    export NXF_SYNTAX_PARSER=v2
    ```

La sortida estarà al mateix lloc que anteriorment, però el contingut hauria d'estar actualitzat amb el nou text.

??? abstract "Contingut del fitxer"

    ```console title="results/hello_world/output.txt"
    Holà mundo!
    ```

Nextflow va utilitzar el valor per defecte del paràmetre de salutació per crear la sortida.

#### 3.4.3. Sobreescriviu el valor per defecte

Si proporcioneu el paràmetre a la línia de comandes, el valor CLI sobreescriurà el valor per defecte.

Proveu-ho:

```bash
nextflow run hello-world.nf --input 'Konnichiwa!'
```

??? success "Sortida de la comanda"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `hello-world.nf` [elegant_faraday] DSL2 - revision: 3539118582

    executor >  local (1)
    [6f/a12a91] sayHello | 1 of 1 ✔
    ```

Una vegada més, hauríeu de trobar la sortida actualitzada corresponent al vostre directori de resultats.

??? abstract "Contingut del fitxer"

    ```console title="results/hello_world/output.txt"
    Konnichiwa!
    ```

!!! note "Nota"

    A Nextflow, hi ha múltiples llocs on podeu especificar valors per als paràmetres.
    Si el mateix paràmetre s'estableix a valors diferents en múltiples llocs, Nextflow determinarà quin valor utilitzar en funció de l'ordre de precedència que es descriu [aquí](https://www.nextflow.io/docs/latest/config.html).

    Cobrirem això amb més detall a la Part 6 (Configuració).

### Conclusió

Sabeu com utilitzar una entrada variable simple proporcionada en temps d'execució mitjançant un paràmetre de línia de comandes, així com configurar, utilitzar i sobreescriure valors per defecte.

### Què segueix?

Apreneu com gestionar execucions de manera més convenient.

---

## 4. Gestioneu execucions de workflow

Saber com llançar workflows i recuperar sortides és genial, però aviat trobareu que hi ha alguns altres aspectes de la gestió de workflows que us facilitaran la vida, especialment si esteu desenvolupant els vostres propis workflows.

Aquí us mostrem com utilitzar la funció [`-resume`](https://nextflow.io/docs/latest/cache-and-resume.html) per quan necessiteu tornar a llançar el mateix workflow, com inspeccionar el registre d'execucions passades amb [`nextflow log`](https://nextflow.io/docs/latest/reference/cli.html#log), i com suprimir directoris work més antics amb [`nextflow clean`](https://nextflow.io/docs/latest/reference/cli.html#clean).

### 4.1. Torneu a llançar un workflow amb `-resume`

De vegades, voldreu tornar a executar un pipeline que ja heu llançat anteriorment sense refer cap pas que ja s'hagi completat correctament.

Nextflow té una opció anomenada [`-resume`](https://nextflow.io/docs/latest/cache-and-resume.html) que us permet fer això.
Específicament, en aquest mode, qualsevol procés que ja s'hagi executat amb exactament el mateix codi, configuració i entrades s'ometrà.
Això significa que Nextflow només executarà processos que hàgiu afegit o modificat des de l'última execució, o als quals proporcioneu noves configuracions o entrades.

Hi ha dos avantatges clau de fer això:

- Si esteu enmig del desenvolupament del vostre pipeline, podeu iterar més ràpidament ja que només heu d'executar el(s) procés(os) en què esteu treballant activament per provar els vostres canvis.
- Si esteu executant un pipeline en producció i alguna cosa va malament, en molts casos podeu solucionar el problema i tornar a llançar el pipeline, i es reprendrà l'execució des del punt de fallada, cosa que us pot estalviar molt de temps i càlcul.

Per utilitzar-lo, simplement afegiu `-resume` a la vostra comanda i executeu-la:

```bash
nextflow run hello-world.nf -resume
```

??? success "Sortida de la comanda"

    ```console hl_lines="5"
    N E X T F L O W   ~  version 25.10.2

    Launching `hello-world.nf` [golden_cantor] DSL2 - revision: 35bd3425e5

    [62/49a1f8] sayHello | 1 of 1, cached: 1 ✔
    ```

La sortida de la consola hauria de semblar familiar, però hi ha una cosa que és una mica diferent en comparació amb abans.

Busqueu la part `cached:` que s'ha afegit a la línia d'estat del procés (línia 5), que significa que Nextflow ha reconegut que ja ha fet aquesta feina i simplement ha reutilitzat el resultat de l'execució exitosa anterior.

També podeu veure que el hash del subdirectori work és el mateix que a l'execució anterior.
Nextflow literalment us està assenyalant l'execució anterior i dient "Ja ho vaig fer allà".

!!! tip "Consell"

    Quan torneu a executar un pipeline amb `resume`, Nextflow no sobreescriu cap fitxer publicat fora del directori work per cap execució que s'hagi executat correctament anteriorment.

### 4.2. Inspeccioneu el registre d'execucions passades

Tant si esteu desenvolupant un nou pipeline com si esteu executant pipelines en producció, en algun moment probablement necessitareu buscar informació sobre execucions passades.
Així és com fer-ho.

Cada vegada que llanceu un workflow de nextflow, s'escriu una línia a un fitxer de registre anomenat `history`, sota un directori ocult anomenat `.nextflow` al directori de treball actual.

??? abstract "Contingut del fitxer"

    ```txt title=".nextflow/history" linenums="1"
    2025-07-04 19:27:09	1.8s	wise_watson	OK	3539118582ccde68dde471cc2c66295c	a02c9c46-c3c7-4085-9139-d1b9b5b194c8	nextflow run 1-hello.nf --input 'Hello World'
    2025-07-04 19:27:20	2.9s	spontaneous_blackwell	OK	3539118582ccde68dde471cc2c66295c	59a5db23-d83c-4c02-a54e-37ddb73a337e	nextflow run 1-hello.nf --input Bonjour
    2025-07-04 19:27:31	1.8s	gigantic_yonath	OK	3539118582ccde68dde471cc2c66295c	5acaa83a-6ad6-4509-bebc-cb25d5d7ddd0	nextflow run 1-hello.nf --input 'Dobry den'
    2025-07-04 19:27:45	2.4s	backstabbing_swartz	OK	3539118582ccde68dde471cc2c66295c	5f4b3269-5b53-404a-956c-cac915fbb74e	nextflow run 1-hello.nf --input Konnichiwa
    2025-07-04 19:27:57	2.1s	goofy_wilson	OK	3539118582ccde68dde471cc2c66295c	5f4b3269-5b53-404a-956c-cac915fbb74e	nextflow run 1-hello.nf --input Konnichiwa -resume
    ```

Aquest fitxer us proporciona la marca de temps, el nom d'execució, l'estat, l'ID de revisió, l'ID de sessió i la línia de comandes completa per a cada execució de Nextflow que s'ha llançat des del directori de treball actual.

Una manera més convenient d'accedir a aquesta informació és utilitzar la comanda `nextflow log`.

```bash
nextflow log
```

??? success "Sortida de la comanda"

    ```console linenums="1"
    TIMESTAMP               DURATION        RUN NAME                STATUS  REVISION ID     SESSION ID                              COMMAND
    2025-07-04 19:27:09     1.8s            wise_watson             OK       3539118582     a02c9c46-c3c7-4085-9139-d1b9b5b194c8    nextflow run 1-hello.nf --input 'Hello World'
    2025-07-04 19:27:20     2.9s            spontaneous_blackwell   OK       3539118582     59a5db23-d83c-4c02-a54e-37ddb73a337e    nextflow run 1-hello.nf --input Bonjour
    2025-07-04 19:27:31     1.8s            gigantic_yonath         OK       3539118582     5acaa83a-6ad6-4509-bebc-cb25d5d7ddd0    nextflow run 1-hello.nf --input 'Dobry den'
    2025-07-04 19:27:45     2.4s            backstabbing_swartz     OK       3539118582     5f4b3269-5b53-404a-956c-cac915fbb74e    nextflow run 1-hello.nf --input Konnichiwa
    2025-07-04 19:27:57     2.1s            goofy_wilson            OK       3539118582     5f4b3269-5b53-404a-956c-cac915fbb74e    nextflow run 1-hello.nf --input Konnichiwa -resume
    ```

Això mostrarà el contingut del fitxer de registre al terminal, augmentat amb una línia de capçalera.

Notareu que l'ID de sessió canvia cada vegada que executeu una nova comanda `nextflow run`, EXCEPTE si esteu utilitzant l'opció `-resume`.
En aquest cas, l'ID de sessió es manté igual.

Nextflow utilitza l'ID de sessió per agrupar informació de memòria cau d'execució sota el directori `cache`, també ubicat sota `.nextflow`.

### 4.3. Suprimiu directoris work més antics

Durant el procés de desenvolupament, normalment executareu el vostre esborrany de pipeline un gran nombre de vegades, cosa que pot conduir a una acumulació de molts fitxers en molts subdirectoris.

Afortunadament, Nextflow inclou una subcomanda `clean` útil que pot suprimir automàticament els subdirectoris work per a execucions passades que ja no us importen.

#### 4.3.1. Determineu els criteris de supressió

Hi ha múltiples [opcions](https://www.nextflow.io/docs/latest/reference/cli.html#clean) per determinar què suprimir.

Aquí us mostrem un exemple que suprimeix tots els subdirectoris d'execucions anteriors a una execució determinada, especificada utilitzant el seu nom d'execució.

Busqueu l'execució exitosa més recent on no vau utilitzar `-resume`; en el nostre cas el nom d'execució era `golden_cantor`.

El nom d'execució és la cadena de dues parts generada per la màquina que es mostra entre claudàtors a la línia de sortida de consola `Launching (...)`.
També podeu utilitzar el registre de Nextflow per buscar una execució en funció de la seva marca de temps i/o línia de comandes.

#### 4.3.2. Feu una execució de prova

Primer utilitzem la bandera d'execució de prova `-n` per comprovar què se suprimirà donada la comanda:

```bash
nextflow clean -before golden_cantor -n
```

??? success "Sortida de la comanda"

    ```console
    Would remove /workspaces/training/hello-nextflow/work/a3/7be2fad5e71e5f49998f795677fd68
    ```

La vostra sortida tindrà noms de directori de tasques diferents i pot tenir un nombre diferent de línies, però hauria de semblar similar a l'exemple.

Si no veieu cap línia de sortida, o bé no heu proporcionat un nom d'execució vàlid o no hi ha execucions passades per suprimir. Assegureu-vos de canviar `golden_cantor` a la comanda d'exemple per qualsevol que sigui el nom d'execució més recent corresponent al vostre registre.

#### 4.3.3. Procediu amb la supressió

Si la sortida sembla com s'esperava i voleu procedir amb la supressió, torneu a executar la comanda amb la bandera `-f` en lloc de `-n`:

```bash
nextflow clean -before golden_cantor -f
```

??? success "Sortida de la comanda"

    ```console
    Removed /workspaces/training/hello-nextflow/work/a3/7be2fad5e71e5f49998f795677fd68
    ```

La sortida hauria de ser similar a abans, però ara dient 'Removed' en lloc de 'Would remove'.
Tingueu en compte que això no elimina els subdirectoris de dos caràcters (com `a3/` anteriorment) però sí que buida el seu contingut.

!!! Warning "Advertència"

    Suprimir subdirectoris work d'execucions passades els elimina de la memòria cau de Nextflow i suprimeix qualsevol sortida que s'hagi emmagatzemat en aquests directoris.
    Això significa que trenca la capacitat de Nextflow de reprendre l'execució sense tornar a executar els processos corresponents.

    Sou responsables de desar qualsevol sortida que us importi o en què planegeu confiar! Aquesta és la raó principal per la qual preferim utilitzar el mode `copy` en lloc del mode `symlink` per a la directiva `publish`.

### Conclusió

Sabeu com publicar sortides a un directori específic, tornar a llançar un pipeline sense repetir passos que ja s'han executat de manera idèntica, i utilitzar la comanda `nextflow clean` per netejar directoris work antics.

Més en general, sabeu com interpretar un workflow simple de Nextflow, gestionar la seva execució i recuperar sortides.

### Què segueix?

Feu una petita pausa, us l'heu guanyat!

Quan estigueu preparats, passeu a [**Part 2: Hello Channels**](./02_hello_channels.md) per aprendre com utilitzar canals per alimentar entrades al vostre workflow, cosa que us permetrà aprofitar el paral·lelisme de flux de dades integrat de Nextflow i altres funcions potents.

---

## Qüestionari

<quiz>
Quins són els components mínims requerits d'un procés de Nextflow?
- [ ] Només blocs d'entrada i sortida
- [x] Blocs de sortida i script
- [ ] Blocs d'entrada, sortida i script
- [ ] Només un bloc script

Més informació: [1.1.1. La definició de procés](#111-the-process-definition)
</quiz>

<quiz>
Quin és el propòsit del bloc de sortida en un procés?
- [ ] Imprimir resultats a la consola
- [ ] Desar fitxers al directori work
- [x] Declarar sortides esperades del procés
- [ ] Definir variables d'entorn

Més informació: [1.1.1. La definició de procés](#111-the-process-definition)
</quiz>

<quiz>
Quina comanda s'utilitza per executar un workflow de Nextflow?
- [ ] `nextflow start`
- [ ] `nextflow execute`
- [x] `nextflow run`
- [ ] `nextflow launch`
</quiz>

<quiz>
Mirant el directori work d'una tasca, quin fitxer conté la comanda real que es va executar?

```
work/a3/7be2fa.../
├── .command.begin
├── .command.err
├── .command.log
├── .command.out
├── .command.run
├── .command.sh
├── .exitcode
└── output.txt
```

- [ ] `.command.run`
- [x] `.command.sh`
- [ ] `.command.log`
- [ ] `.command.out`

Més informació: [1.2.2. Trobeu la sortida i els registres al directori `work`](#122-find-the-output-and-logs-in-the-work-directory)
</quiz>

<quiz>
Què fa la bandera `-resume`?
- [ ] Reinicia el workflow des del principi
- [ ] Pausa el workflow
- [x] Omet processos que ja s'han completat correctament
- [ ] Crea una còpia de seguretat del workflow

Més informació: [4.1. Torneu a llançar un workflow amb `-resume`](#41-re-launch-a-workflow-with--resume)
</quiz>

<quiz>
Quin és el mode per defecte per publicar sortides de workflow?
- [ ] Copiar fitxers al directori de sortida
- [x] Crear enllaços simbòlics al directori de sortida
- [ ] Moure fitxers al directori de sortida
- [ ] Comprimir fitxers al directori de sortida

Més informació: [2.3. Establiu el mode de publicació a còpia](#23-set-the-publish-mode-to-copy)
</quiz>

<quiz>
Com passeu un valor de paràmetre a un workflow de Nextflow des de la línia de comandes?
- [ ] `-parameter value`
- [ ] `--parameter:value`
- [x] `--parameter value`
- [ ] `-p parameter=value`

Més informació: [3.2. Configureu un paràmetre de línia de comandes per capturar l'entrada de l'usuari](#32-set-up-a-command-line-parameter-to-capture-user-input)
</quiz>

<quiz>
Com es fa referència a una variable dins d'un bloc script de Nextflow?
- [ ] Utilitzeu la sintaxi `%variable%`
- [x] Utilitzeu la sintaxi `#!groovy ${variable}`
- [ ] Utilitzeu la sintaxi `{{variable}}`
- [ ] Utilitzeu la sintaxi `[variable]`
</quiz>
