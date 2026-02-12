# Part 1: Executar operacions bàsiques

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Traducció assistida per IA - [més informació i suggeriments](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

En aquesta primera part del curs de formació Nextflow Run, comencem el tema amb un exemple molt bàsic de Hello World independent del domini, que utilitzarem per demostrar operacions essencials i assenyalar els components de codi Nextflow corresponents.

??? info "Què és un exemple Hello World?"

    Un "Hello World!" és un exemple minimalista que pretén demostrar la sintaxi bàsica i l'estructura d'un llenguatge de programació o marc de programari.
    L'exemple normalment consisteix a imprimir la frase "Hello, World!" al dispositiu de sortida, com ara la consola o el terminal, o escriure-la en un fitxer.

---

## 1. Executar un Hello World directament

Demostrem aquest concepte amb una comanda simple que executem directament al terminal, per mostrar què fa abans d'embolicar-la en Nextflow.

!!! tip "Consell"

    Recordeu que ara hauríeu d'estar dins del directori `nextflow-run/` tal com es descriu a la pàgina [Primers passos](00_orientation.md).

### 1.1. Fer que el terminal digui hola

Executeu la comanda següent al vostre terminal.

```bash
echo 'Hello World!'
```

??? success "Sortida de la comanda"

    ```console
    Hello World!
    ```

Això mostra el text 'Hello World' directament al terminal.

### 1.2. Escriure la sortida a un fitxer

Executar pipelines implica principalment llegir dades de fitxers i escriure resultats a altres fitxers, així que modifiquem la comanda per escriure la sortida de text a un fitxer per fer l'exemple una mica més rellevant.

```bash
echo 'Hello World!' > output.txt
```

??? success "Sortida de la comanda"

    ```console

    ```

Això no mostra res al terminal.

### 1.3. Trobar la sortida

El text 'Hello World' ara hauria d'estar al fitxer de sortida que hem especificat, anomenat `output.txt`.
Podeu obrir-lo a l'explorador de fitxers o des de la línia de comandes utilitzant la utilitat `cat`, per exemple.

??? abstract "Contingut del fitxer"

    ```console title="output.txt" linenums="1"
    Hello World!
    ```

Això és el que intentarem replicar amb el nostre primer workflow de Nextflow.

### Conclusió

Ara sabeu com executar una comanda simple al terminal que mostra algun text, i opcionalment, com fer que escrigui la sortida a un fitxer.

### Què segueix?

Descobriu què cal fer per executar un workflow de Nextflow que aconsegueixi el mateix resultat.

---

## 2. Executar el workflow

Us proporcionem un script de workflow anomenat `1-hello.nf` que pren una salutació d'entrada mitjançant un argument de línia de comandes anomenat `--input` i produeix un fitxer de text que conté aquesta salutació.

No mirarem el codi encara; primer vegem com és executar-lo.

### 2.1. Llançar el workflow i monitoritzar l'execució

Al terminal, executeu la comanda següent.

```bash
nextflow run 1-hello.nf --input 'Hello World!'
```

??? success "Sortida de la comanda"

    ```console hl_lines="6"
    N E X T F L O W   ~  version 25.10.2

    Launching `1-hello.nf` [goofy_torvalds] DSL2 - revision: c33d41f479

    executor >  local (1)
    [a3/7be2fa] sayHello | 1 of 1 ✔
    ```

Si la sortida de la vostra consola s'assembla a això, felicitats, acabeu d'executar el vostre primer workflow de Nextflow!

??? question "Si no ha funcionat"

    Si ha fallat amb un error que s'assembla a això:

    ```
    Parameter `input` was specified on the command line or params file but is not declared in the script or config

    -- Check script '1-hello.nf' at line: 23 or see '.nextflow.log' file for more details
    ```

    Probablement esteu utilitzant l'analitzador de llenguatge v1 de Nextflow més antic.
    Això es va esmentar al començament del curs, però potser us ho heu perdut.
    Consulteu el material d'ajuda sobre [versions de Nextflow](../info/nxf_versions.md).

    En resum, si utilitzeu Nextflow `25.10` necessiteu habilitar l'analitzador de llenguatge v2:

    ```bash
    export NXF_SYNTAX_PARSER=v2
    ```

La sortida més important aquí és l'última línia, que està ressaltada a la sortida anterior:

```console
[a3/7be2fa] sayHello | 1 of 1 ✔
```

Això ens diu que el procés `sayHello` s'ha executat amb èxit una vegada (`1 of 1 ✔`).

Això és genial, però potser us esteu preguntant: on és la sortida?

### 2.2. Trobar el fitxer de sortida al directori `results`

Aquest workflow està configurat per publicar la seva sortida a un directori de resultats.
Si mireu el vostre directori actual, veureu que quan heu executat el workflow, Nextflow ha creat un directori nou anomenat `results`, així com un subdirectori anomenat `1-hello` sota aquest, que conté un fitxer anomenat `output.txt`.

```console title="results/"
results
└── 1-hello
    └── output.txt
```

Obriu el fitxer; el contingut hauria de coincidir amb la cadena que heu especificat a la línia de comandes.

```console title="results/1-hello/output.txt" linenums="1"
Hello World!
```

Genial, el nostre workflow ha fet el que havia de fer!

### 2.3. Desar els resultats a un directori diferent

Per defecte, Nextflow desarà les sortides del pipeline a un directori anomenat `results` al vostre camí actual.
Per canviar on es publiquen els vostres fitxers, utilitzeu la bandera CLI `-output-dir` (o `-o` de forma abreujada)

!!! danger "Perill"

    Tingueu en compte que `--input` té dos guions i `-output-dir` en té un!
    Això és perquè `--input` és un _paràmetre_ del pipeline i `-output-dir` és una bandera CLI bàsica de Nextflow.
    Més sobre això més endavant.

```bash
nextflow run 1-hello.nf --input 'Hello World!' -output-dir hello_results
```

??? success "Sortida de la comanda"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `1-hello.nf` [hungry_celsius] DSL2 - revision: f048d6ea78

    executor >  local (1)
    [a3/1e1535] sayHello [100%] 1 of 1 ✔
    ```

Hauríeu de veure que les vostres sortides ara es publiquen a un directori anomenat `hello_results` en lloc de `results`:

```console title="hello_results/"
hello_results
└── 1-hello
    └── output.txt
```

Els fitxers dins d'aquest directori són exactament els mateixos que abans, només és diferent el directori de nivell superior.
Tanmateix, tingueu en compte que en tots dos casos el resultat 'publicat' és una còpia (o en alguns casos un enllaç simbòlic) de la sortida real produïda per Nextflow quan ha executat el workflow.

Així que ara, mirarem sota el capó per veure on Nextflow ha executat realment el treball.

!!! Warning "Advertència"

    No tots els workflows estaran configurats per publicar sortides a un directori de resultats, i/o els noms i l'estructura dels directoris poden ser diferents.
    Una mica més endavant en aquesta secció, us mostrarem com esbrinar on s'especifica aquest comportament.

### 2.4. Trobar la sortida original i els registres al directori `work/`

Quan executeu un workflow, Nextflow crea un 'directori de tasca' distint per a cada invocació de cada procés del workflow (=cada pas del pipeline).
Per a cadascun, prepararà les entrades necessàries, executarà la(les) instrucció(ns) rellevant(s) i escriurà sortides i fitxers de registre dins d'aquest directori, que es nomena automàticament utilitzant un hash per fer-lo únic.

Tots aquests directoris de tasca viuran sota un directori anomenat `work` dins del vostre directori actual (on esteu executant la comanda).

Això pot sonar confús, així que vegem com es veu a la pràctica.

Tornant a la sortida de consola del workflow que hem executat abans, teníem aquesta línia:

```console
[a3/1e1535] sayHello [100%] 1 of 1 ✔
```

Veieu com la línia comença amb `[a3/1e1535]`?
Aquesta és una forma truncada del camí del directori de tasca per a aquesta crida de procés, i us indica on trobar la sortida de la crida del procés `sayHello` dins del camí del directori `work/`.

Podeu trobar el camí complet escrivint la comanda següent (substituint `a3/1e1535` pel que veieu al vostre propi terminal) i prement la tecla tab per autocompletar el camí o afegint un asterisc:

```bash
ls work/a3/1e1535*
```

Això hauria de donar el camí complet del directori: `work/a3/1e153543b0a7f9d2c4735ddb4ab231`

Vegem què hi ha dins.

??? abstract "Contingut del directori"

    ```console
    work
    ├── a3
    │   └── 1e153543b0a7f9d2c4735ddb4ab231
    │       ├── .command.begin
    │       ├── .command.err
    │       ├── .command.log
    │       ├── .command.out
    │       ├── .command.run
    │       ├── .command.sh
    │       ├── .exitcode
    │       └── output.txt
    └── a4
        └── aa3694b8808bdcc1135ef4a1187a4d
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

    Si navegueu pel contingut del subdirectori de tasca a l'explorador de fitxers de VSCode, veureu tots els fitxers immediatament.
    Tanmateix, els fitxers de registre estan configurats per ser invisibles al terminal, així que si voleu utilitzar `ls` o `tree` per veure'ls, haureu d'establir l'opció rellevant per mostrar fitxers invisibles.

    ```bash
    tree -a work
    ```

Hi ha dos conjunts de directoris a `work/`, de les dues execucions diferents del pipeline que hem fet.
Cada execució de tasca obté el seu propi directori aïllat per treballar.
En aquest cas el pipeline ha fet el mateix les dues vegades, així que el contingut de cada directori de tasca és idèntic

Hauríeu de reconèixer immediatament el fitxer `output.txt`, que de fet és la sortida original del procés `sayHello` que es va publicar al directori `results`.
Si l'obriu, trobareu la salutació `Hello World!` de nou.

```console title="work/a3/1e153543b0a7f9d2c4735ddb4ab231/output.txt"
Hello World!
```

I què passa amb tots aquests altres fitxers?

Aquests són els fitxers auxiliars i de registre que Nextflow va escriure com a part de l'execució de la tasca:

- **`.command.begin`**: Fitxer sentinella creat tan aviat com es llança la tasca.
- **`.command.err`**: Missatges d'error (`stderr`) emesos per la crida del procés
- **`.command.log`**: Sortida de registre completa emesa per la crida del procés
- **`.command.out`**: Sortida regular (`stdout`) per la crida del procés
- **`.command.run`**: Script complet executat per Nextflow per executar la crida del procés
- **`.command.sh`**: La comanda que realment va executar la crida del procés
- **`.exitcode`**: El codi de sortida resultant de la comanda

El fitxer `.command.sh` és especialment útil perquè us mostra la comanda principal que Nextflow va executar, sense incloure tota la comptabilitat i configuració de tasca/entorn.

```console title="work/a3/1e153543b0a7f9d2c4735ddb4ab231/.command.sh"
#!/bin/bash -ue
echo 'Hello World!' > output.txt

```

Així que això confirma que el workflow va compondre la mateixa comanda que vam executar directament a la línia de comandes abans.

Quan alguna cosa va malament i necessiteu solucionar problemes sobre què ha passat, pot ser útil mirar l'script `command.sh` per comprovar exactament quina comanda va compondre Nextflow basant-se en les instruccions del workflow, la interpolació de variables, etc.

### 2.5. Tornar a executar el workflow amb diferents salutacions

Proveu de tornar a executar el workflow diverses vegades amb diferents valors per a l'argument `--input`, després mireu els directoris de tasca.

??? abstract "Contingut del directori"

    ```console
    work/
    ├── 09
    │   └── 5ea8665939daf6f04724286c9b3c8a
    │       ├── .command.begin
    │       ├── .command.err
    │       ├── .command.log
    │       ├── .command.out
    │       ├── .command.run
    │       ├── .command.sh
    │       ├── .exitcode
    │       └── output.txt
    ├── 92
    │   └── ceb95e05d87621c92a399da9bd2067
    │       ├── .command.begin
    │       ├── .command.err
    │       ├── .command.log
    │       ├── .command.out
    │       ├── .command.run
    │       ├── .command.sh
    │       ├── .exitcode
    │       └── output.txt
    ├── 93
    │   └── 6708dbc20c7efdc6769cbe477061ec
    │       ├── .command.begin
    │       ├── .command.err
    │       ├── .command.log
    │       ├── .command.out
    │       ├── .command.run
    │       ├── .command.sh
    │       ├── .exitcode
    │       └── output.txt
    ├── a3
    │   └── 1e153543b0a7f9d2c4735ddb4ab231
    │       ├── .command.begin
    │       ├── .command.err
    │       ├── .command.log
    │       ├── .command.out
    │       ├── .command.run
    │       ├── .command.sh
    │       ├── .exitcode
    │       └── output.txt
    └── a4
        └── aa3694b8808bdcc1135ef4a1187a4d
            ├── .command.begin
            ├── .command.err
            ├── .command.log
            ├── .command.out
            ├── .command.run
            ├── .command.sh
            ├── .exitcode
            └── output.txt
    ```

Veieu que s'ha creat un subdirectori nou amb un conjunt complet de fitxers de sortida i registre per a cada execució.

En canvi, si mireu el directori `results`, encara només hi ha un conjunt de resultats, i el contingut del fitxer de sortida correspon al que heu executat per última vegada.

??? abstract "Contingut del directori"

    ```console title="results/"
    results
    └── 1-hello
        └── output.txt
    ```

Això us mostra que els resultats publicats se sobreescriuran per execucions posteriors, mentre que els directoris de tasca sota `work/` es preserven.

### Conclusió

Sabeu com executar un script simple de Nextflow, monitoritzar la seva execució i trobar les seves sortides.

### Què segueix?

Apreneu a llegir un script bàsic de Nextflow i identifiqueu com els seus components es relacionen amb la seva funcionalitat.

---

## 3. Examinar l'script inicial del workflow Hello World

El que hem fet allà bàsicament era tractar l'script del workflow com una caixa negra.
Ara que hem vist què fa, obrim la caixa i mirem a dins.

El nostre objectiu aquí no és memoritzar la sintaxi del codi Nextflow, sinó formar alguna intuïció bàsica sobre quins són els components principals i com estan organitzats.

### 3.1. Examinar l'estructura general del codi

Trobareu l'script `1-hello.nf` al vostre directori actual, que hauria de ser `nextflow-run`. Obriu-lo al panell de l'editor.

??? full-code "Fitxer de codi complet"

    ```groovy title="1-hello.nf" linenums="1"
    #!/usr/bin/env nextflow

    /*
    * Use echo to print 'Hello World!' to a file
    */
    process sayHello {

        input:
        val greeting

        output:
        path 'output.txt'

        script:
        """
        echo '${greeting}' > output.txt
        """
    }

    /*
    * Pipeline parameters
    */
    params {
        input: String
    }

    workflow {

        main:
        // emit a greeting
        sayHello(params.input)

        publish:
        first_output = sayHello.out
    }

    output {
        first_output {
            path '1-hello'
            mode 'copy'
        }
    }
    ```

Un script de workflow de Nextflow normalment inclou una o més definicions de **process**, el **workflow** en si, i alguns blocs opcionals com **params** i **output**.

Cada **process** descriu quina(es) operació(ns) hauria d'acomplir el pas corresponent del pipeline, mentre que el **workflow** descriu la lògica de flux de dades que connecta els diversos passos.

Vegem més de prop primer el bloc **process**, després mirarem el bloc **workflow**.

### 3.2. La definició del `process`

El primer bloc de codi descriu un [**process**](https://nextflow.io/docs/latest/process.html).
La definició del procés comença amb la paraula clau `process`, seguida del nom del procés i finalment el cos del procés delimitat per claus.
El cos del procés ha de contenir un bloc script que especifica la comanda a executar, que pot ser qualsevol cosa que poguéssiu executar en un terminal de línia de comandes.

```groovy title="1-hello.nf" linenums="3"
/*
* Use echo to print a greeting to a file
*/
process sayHello {

    input:
    val greeting

    output:
    path 'output.txt'

    script:
    """
    echo '${greeting}' > output.txt
    """
}
```

Aquí tenim un **process** anomenat `sayHello` que pren una variable d'**input** anomenada `greeting` i escriu la seva **output** a un fitxer anomenat `output.txt`.

<figure class="excalidraw">
--8<-- "docs/en/docs/nextflow_run/img/sayhello_with_input.svg"
</figure>

Aquesta és una definició de procés molt mínima que només conté una definició d'`input`, una definició d'`output` i l'`script` a executar.

La definició d'`input` inclou el qualificador `val`, que indica a Nextflow que esperi un valor d'algun tipus (pot ser una cadena, un número, el que sigui).

La definició d'`output` inclou el qualificador `path`, que indica a Nextflow que això s'ha de gestionar com un camí (inclou tant camins de directori com fitxers).

### 3.3. La definició del `workflow`

El segon bloc de codi descriu el [**workflow**](https://nextflow.io/docs/latest/workflow.html) en si.
La definició del workflow comença amb la paraula clau `workflow`, seguida d'un nom opcional, després el cos del workflow delimitat per claus.

Aquí tenim un **workflow** que consisteix en un bloc `main:` i un bloc `publish:`.
El bloc `main:` és el cos principal del workflow i el bloc `publish:` llista les sortides que s'han de publicar al directori `results`.

```groovy title="1-hello.nf" linenums="27"
workflow {

    main:
    // emet una salutacio
    sayHello(params.input)

    publish:
    first_output = sayHello.out
}
```

En aquest cas el bloc `main:` conté una crida al procés `sayHello` i li dóna una entrada anomenada `params.input` per utilitzar com a salutació.

Com discutirem amb més detall d'aquí a un moment, `params.input` conté el valor que vam donar al paràmetre `--input` a la nostra línia de comandes.

El bloc `publish:` llista la sortida de la crida del procés `sayHello()`, a la qual es refereix com `sayHello.out` i dóna el nom `first_output` (això pot ser qualsevol cosa que l'autor del workflow vulgui).

Aquesta és una definició de **workflow** molt mínima.
En un pipeline del món real, el workflow normalment conté múltiples crides a **processes** connectats per **channels**, i pot haver-hi valors per defecte configurats per a les entrades variables.

Entrarem en això a la Part 2 del curs.
Per ara, vegem més de prop com el nostre workflow està gestionant entrades i sortides.

### 3.4. El sistema `params` de paràmetres de línia de comandes

El `params.input` que proporcionem a la crida del procés `sayHello()` és un fragment de codi Nextflow interessant i val la pena dedicar-hi un minut extra.

Com s'ha esmentat anteriorment, així és com passem el valor del paràmetre de línia de comandes `--input` a la crida del procés `sayHello()`.
De fet, simplement declarar `params.someParameterName` és suficient per donar al workflow un paràmetre anomenat `--someParameterName` des de la línia de comandes.

Aquí hem formalitzat aquesta declaració de paràmetre configurant un bloc `params` que especifica el tipus d'entrada que el workflow espera (Nextflow 25.10.2 i posteriors).

```groovy title="1-hello.nf" linenums="20"
/*
 * Pipeline parameters
 */
params {
    input: String
}
```

Els tipus suportats inclouen `String`, `Integer`, `Float`, `Boolean` i `Path`.
Per aprendre més, consulteu [Workflow parameters](https://nextflow.io/docs/latest/config.html#workflow-parameters) a la documentació de referència de Nextflow.

!!! tip "Consell"

    Recordeu que els paràmetres de _workflow_ declarats utilitzant el sistema `params` sempre prenen dos guions a la línia de comandes (`--`).
    Això els distingeix de les banderes CLI de _nivell Nextflow_, que només prenen un guió (`-`).

### 3.5. La directiva `publish`

A l'altre extrem del workflow, ja hem donat una ullada al bloc `publish:`.
Aquesta és una meitat del sistema de gestió de sortides; l'altra meitat és el bloc `output` situat a sota.

```groovy title="1-hello.nf" linenums="37"
output {
    first_output {
        path '1-hello'
        mode 'copy'
    }
}
```

Això especifica que la sortida `first_output` llistada al bloc `publish:` s'ha de copiar a un subdirectori anomenat `1-hello` sota el directori de sortida `results` per defecte.

La línia `mode 'copy'` sobreescriu el comportament per defecte del sistema, que és fer un enllaç simbòlic (o symlink) al fitxer original al directori `work/` en lloc d'una còpia pròpia.

Hi ha més opcions que les mostrades aquí per controlar el comportament de publicació; en cobrirem algunes més endavant.
També veureu que quan un workflow genera múltiples sortides, cadascuna es llista d'aquesta manera al bloc `output`.

Per aprendre més, consulteu [Publishing outputs](https://nextflow.io/docs/latest/workflow.html#publishing-outputs) a la documentació de referència de Nextflow.

??? info "Sintaxi antiga per publicar sortides utilitzant `publishDir`"

    Fins fa molt poc, la manera establerta de publicar sortides era fer-ho al nivell de cada procés individual utilitzant una directiva `publishDir`.

    Encara trobareu aquest patró de codi per tot arreu en pipelines de Nextflow més antics i mòduls de procés, així que és important ser-ne conscient.

    En lloc de tenir un bloc `publish:` al workflow i un bloc `output` al nivell superior, veuríeu una línia `publishDir` a la definició del procés `sayHello`:

    ```groovy title="Syntax example" linenums="1" hl_lines="3"
    process sayHello {

        publishDir 'results/1-hello', mode: 'copy'

        output:
        path 'output.txt'

        script:
        """
        echo 'Hello World!' > output.txt
        """
    }
    ```

    Tanmateix, no recomanem utilitzar això en cap treball nou ja que eventualment no es permetrà en versions futures del llenguatge Nextflow.

### Conclusió

Ara sabeu com està estructurat un workflow simple de Nextflow, i com els components bàsics es relacionen amb la seva funcionalitat.

### Què segueix?

Apreneu a gestionar les vostres execucions de workflow de manera convenient.

---

## 4. Gestionar execucions de workflow

Saber com llançar workflows i recuperar sortides és genial, però ràpidament trobareu que hi ha alguns altres aspectes de la gestió de workflows que us facilitaran la vida.

Aquí us mostrem com aprofitar la funció `resume` per quan necessiteu tornar a llançar el mateix workflow, com inspeccionar els registres d'execució amb `nextflow log`, i com eliminar directoris de treball més antics amb `nextflow clean`.

### 4.1. Tornar a llançar un workflow amb `-resume`

De vegades, voldreu tornar a executar un pipeline que ja heu llançat prèviament sense refer cap treball que ja s'hagi completat amb èxit.

Nextflow té una opció anomenada `-resume` que us permet fer això.
Específicament, en aquest mode, qualsevol procés que ja s'hagi executat amb exactament el mateix codi, configuració i entrades s'ometrà.
Això significa que Nextflow només executarà processos que hàgiu afegit o modificat des de l'última execució, o als quals estigueu proporcionant noves configuracions o entrades.

Hi ha dos avantatges clau de fer això:

- Si esteu desenvolupant un pipeline, podeu iterar més ràpidament ja que només heu d'executar el(s) procés(os) en què esteu treballant activament per provar els vostres canvis.
- Si esteu executant un pipeline en producció i alguna cosa va malament, en molts casos podeu solucionar el problema i tornar a llançar el pipeline, i es reprendrà l'execució des del punt de fallada, cosa que us pot estalviar molt de temps i càlcul.

Per utilitzar-lo, simplement afegiu `-resume` a la vostra comanda i executeu-la:

```bash
nextflow run 1-hello.nf --input 'Hello World!' -resume
```

??? success "Sortida de la comanda"

    ```console linenums="1"
    N E X T F L O W   ~  version 25.10.2

    Launching `1-hello.nf` [tiny_noyce] DSL2 - revision: c33d41f479

    [a3/7be2fa] sayHello | 1 of 1, cached: 1 ✔
    ```

La sortida de consola hauria de semblar familiar, però hi ha una cosa que és una mica diferent en comparació amb abans.

Busqueu la part `cached:` que s'ha afegit a la línia d'estat del procés (línia 5), que significa que Nextflow ha reconegut que ja ha fet aquest treball i simplement ha reutilitzat el resultat de l'execució anterior amb èxit.

També podeu veure que el hash del subdirectori de treball és el mateix que a l'execució anterior.
Nextflow literalment us està assenyalant l'execució anterior i dient "Ja ho vaig fer allà."

!!! tip "Consell"

    Quan torneu a executar un pipeline amb `resume`, Nextflow no sobreescriu cap fitxer publicat fora del directori de treball per cap execució que s'hagi executat amb èxit prèviament.

    Per aprendre més, consulteu [Cache and resume](https://nextflow.io/docs/latest/cache-and-resume.html) a la documentació de referència de Nextflow.

### 4.2. Inspeccionar el registre d'execucions passades

Cada vegada que llanceu un workflow de nextflow, s'escriu una línia a un fitxer de registre anomenat `history`, sota un directori ocult anomenat `.nextflow` al directori de treball actual.

??? abstract "Contingut del fitxer"

    ```txt title=".nextflow/history" linenums="1"
    2025-07-04 19:27:09	1.8s	wise_watson	OK	3539118582ccde68dde471cc2c66295c	a02c9c46-c3c7-4085-9139-d1b9b5b194c8	nextflow run 1-hello.nf --input 'Hello World'
    2025-07-04 19:27:20	2.9s	spontaneous_blackwell	OK	3539118582ccde68dde471cc2c66295c	59a5db23-d83c-4c02-a54e-37ddb73a337e	nextflow run 1-hello.nf --input Bonjour
    2025-07-04 19:27:31	1.8s	gigantic_yonath	OK	3539118582ccde68dde471cc2c66295c	5acaa83a-6ad6-4509-bebc-cb25d5d7ddd0	nextflow run 1-hello.nf --input 'Dobry den'
    2025-07-04 19:27:45	2.4s	backstabbing_swartz	OK	3539118582ccde68dde471cc2c66295c	5f4b3269-5b53-404a-956c-cac915fbb74e	nextflow run 1-hello.nf --input Konnichiwa
    2025-07-04 19:27:57	2.1s	goofy_wilson	OK	3539118582ccde68dde471cc2c66295c	5f4b3269-5b53-404a-956c-cac915fbb74e	nextflow run 1-hello.nf --input Konnichiwa -resume
    ```

Aquest fitxer us dóna la marca de temps, el nom d'execució, l'estat, l'ID de revisió, l'ID de sessió i la línia de comandes completa per a cada execució de Nextflow que s'ha llançat des del directori de treball actual.

Una manera més convenient d'accedir a aquesta informació és utilitzar la comanda [`nextflow log`](https://nextflow.io/docs/latest/reference/cli.html#log).

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

### 4.3. Eliminar directoris de treball més antics

Si executeu molts pipelines, podeu acabar acumulant molts fitxers a través de molts subdirectoris.
Com que els subdirectoris es nomenen aleatòriament, és difícil dir pels seus noms quines són execucions més antigues vs. més recents.

Afortunadament Nextflow inclou una comanda útil anomenada [`nextflow clean`](https://www.nextflow.io/docs/latest/reference/cli.html#clean) que pot eliminar automàticament els subdirectoris de treball per a execucions passades que ja no us importen.

#### 4.3.1. Determinar criteris d'eliminació

Hi ha múltiples opcions per determinar què eliminar, que podeu explorar a la documentació enllaçada anteriorment.
Aquí us mostrem un exemple que elimina tots els subdirectoris d'execucions anteriors a una execució donada, especificada utilitzant el seu nom d'execució.

Busqueu l'execució amb èxit més recent on no vau utilitzar `-resume`; en el nostre cas el nom d'execució era `backstabbing_swartz`.

El nom d'execució és la cadena de dues parts generada per màquina que es mostra entre claudàtors a la línia de sortida de consola `Launching (...)`.
També podeu utilitzar el registre de Nextflow per buscar una execució basant-vos en la seva marca de temps i/o línia de comandes.

#### 4.3.2. Fer una execució de prova

Primer utilitzem la bandera d'execució de prova `-n` per comprovar què s'eliminarà donada la comanda:

```bash
nextflow clean -before backstabbing_swartz -n
```

??? success "Sortida de la comanda"

    ```console
    Would remove /workspaces/training/hello-nextflow/work/eb/1a5de36637b475afd88fca7f79e024
    Would remove /workspaces/training/hello-nextflow/work/6b/19b0e002ea13486d3a0344c336c1d0
    Would remove /workspaces/training/hello-nextflow/work/45/9a6dd7ab771f93003d040956282883
    ```

La vostra sortida tindrà noms de directori de tasca diferents i pot tenir un nombre diferent de línies, però hauria de semblar similar a l'exemple.

Si no veieu cap línia de sortida, o no heu proporcionat un nom d'execució vàlid o no hi ha execucions passades per eliminar. Assegureu-vos de canviar `backstabbing_swartz` a la comanda d'exemple pel que sigui el nom d'execució més recent corresponent al vostre registre.

#### 4.3.3. Procedir amb l'eliminació

Si la sortida sembla com s'esperava i voleu procedir amb l'eliminació, torneu a executar la comanda amb la bandera `-f` en lloc de `-n`:

```bash
nextflow clean -before backstabbing_swartz -f
```

??? success "Sortida de la comanda"

    ```console
    Removed /workspaces/training/hello-nextflow/work/eb/1a5de36637b475afd88fca7f79e024
    Removed /workspaces/training/hello-nextflow/work/6b/19b0e002ea13486d3a0344c336c1d0
    Removed /workspaces/training/hello-nextflow/work/45/9a6dd7ab771f93003d040956282883
    ```

La sortida hauria de ser similar a abans, però ara dient 'Removed' en lloc de 'Would remove'.
Tingueu en compte que això no elimina els subdirectoris de dos caràcters (com `eb/` anteriorment) però sí que buida el seu contingut.

!!! Warning "Advertència"

    Eliminar subdirectoris de treball d'execucions passades els elimina de la memòria cau de Nextflow i elimina qualsevol sortida que s'hagi emmagatzemat en aquests directoris.
    Això significa que trenca la capacitat de Nextflow de reprendre l'execució sense tornar a executar els processos corresponents.

    Sou responsables de desar qualsevol sortida que us importi! Aquesta és la raó principal per la qual preferim utilitzar el mode `copy` en lloc del mode `symlink` per a la directiva `publish`.

### Conclusió

Sabeu com tornar a llançar un pipeline sense repetir passos que ja s'han executat de manera idèntica, inspeccionar el registre d'execució i utilitzar la comanda `nextflow clean` per netejar directoris de treball antics.

### Què segueix?

Feu una petita pausa! Acabeu d'absorbir els blocs de construcció de la sintaxi de Nextflow i les instruccions d'ús bàsiques.

A la següent secció d'aquesta formació, mirarem quatre versions successivament més realistes del pipeline Hello World que demostraran com Nextflow us permet processar múltiples entrades de manera eficient, executar workflows compostos de múltiples passos connectats entre si, aprofitar components de codi modulars i utilitzar contenidors per a una major reproductibilitat i portabilitat.

---

## Qüestionari

<quiz>
A la línia de sortida de consola `[a3/7be2fa] SAYHELLO | 1 of 1 ✔`, què representa `[a3/7be2fa]`?
- [ ] El número de versió del procés
- [ ] Un identificador d'execució únic
- [x] El camí truncat al directori de treball de la tasca
- [ ] El checksum del fitxer de sortida

Aprèn més: [2.3. Trobar la sortida original i els registres al directori `work/`](#23-find-the-original-output-and-logs-in-the-work-directory)
</quiz>

<quiz>
Quin és el propòsit del fitxer `.command.sh` en un directori de tasca?
- [ ] Emmagatzema la configuració de la tasca
- [x] Mostra la comanda real que va ser executada pel procés
- [ ] Conté missatges d'error de tasques fallides
- [ ] Llista fitxers d'entrada preparats per a la tasca

Aprèn més: [2.3. Trobar la sortida original i els registres al directori `work/`](#23-find-the-original-output-and-logs-in-the-work-directory)
</quiz>

<quiz>
Què passa amb els resultats publicats quan torneu a executar un workflow sense `-resume`?
- [ ] Es preserven en directoris separats amb marca de temps
- [x] Se sobreescriuen per la nova execució
- [ ] Nextflow impedeix la sobreescriptura i falla
- [ ] Es fa una còpia de seguretat automàticament

Aprèn més: [2.4. Tornar a executar el workflow amb diferents salutacions](#24-re-run-the-workflow-with-different-greetings)
</quiz>

<quiz>
Què indica aquesta sortida de consola?

```console
[skipped  ] process > sayHello (1) [100%] 1 of 1, cached: 1 ✔
```

- [ ] La tasca ha fallat i s'ha omès
- [ ] La tasca està esperant en una cua
- [x] Nextflow ha reutilitzat resultats d'una execució idèntica anterior
- [ ] La tasca s'ha cancel·lat manualment

Aprèn més: [4.1. Tornar a llançar un workflow amb `-resume`](#41-re-launch-a-workflow-with--resume)
</quiz>

<quiz>
On emmagatzema Nextflow l'historial d'execució que mostra la comanda `nextflow log`?
- [ ] Al directori results
- [ ] Al directori work
- [x] Al fitxer `.nextflow/history`
- [ ] A `nextflow.config`

Aprèn més: [4.2. Inspeccionar el registre d'execucions passades](#42-inspect-the-log-of-past-executions)
</quiz>

<quiz>
Quin és el propòsit del bloc `params` en un fitxer de workflow?
- [ ] Definir requisits de recursos del procés
- [ ] Configurar l'executor
- [x] Declarar i tipificar paràmetres d'entrada del workflow
- [ ] Especificar opcions de publicació de sortida

Aprèn més: [3.4. El sistema params de paràmetres de línia de comandes](#34-the-params-system-of-command-line-parameters)
</quiz>

<quiz>
Al bloc `output` del workflow, què fa `mode 'copy'`?
- [ ] Crea una còpia de seguretat del directori de treball
- [x] Fa una còpia completa dels fitxers en lloc d'enllaços simbòlics
- [ ] Copia l'script del workflow als resultats
- [ ] Habilita la còpia incremental de fitxers

Aprèn més: [3.5. La directiva publish](#35-the-publish-directive)
</quiz>

<quiz>
Quina és la bandera recomanada per utilitzar amb la comanda `nextflow clean` abans d'eliminar fitxers realment?
- [x] `-n` (execució de prova) per previsualitzar què s'eliminaria
- [ ] `-v` (verbose) per veure sortida detallada
- [ ] `-a` (all) per seleccionar tots els directoris
- [ ] `-q` (quiet) per suprimir advertències

Aprèn més: [4.3. Eliminar directoris de treball més antics](#43-delete-older-work-directories)
</quiz>
