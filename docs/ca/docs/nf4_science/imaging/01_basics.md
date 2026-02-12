# Part 1: Executar operacions bàsiques

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Traducció assistida per IA - [més informació i suggeriments](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

En aquesta primera part del curs de formació de Nextflow per a Bioimatge, utilitzarem un exemple molt bàsic de Hello World independent del domini per demostrar operacions essencials i assenyalar els components de codi de Nextflow corresponents.

## 1. Executar el workflow

Us proporcionem un script de workflow anomenat `hello-world.nf` que pren una entrada mitjançant un argument de línia de comandes anomenat `--greeting` i produeix un fitxer de text que conté aquesta salutació.
Encara no mirarem el codi; primer vegem com és executar-lo.

### 1.1. Llançar el workflow i monitoritzar l'execució

Al terminal, executeu la comanda següent:

```bash
nextflow run hello-world.nf --greeting 'Hello World!'
```

La sortida de la consola hauria de semblar-se a això:

```console title="Output" linenums="1"
 N E X T F L O W   ~  version 25.04.3

Launching `hello-world.nf` [goofy_torvalds] DSL2 - revision: c33d41f479

executor >  local (1)
[a3/7be2fa] sayHello | 1 of 1 ✔
```

Felicitats, acabeu d'executar el vostre primer workflow de Nextflow!

La sortida més important aquí és l'última línia (línia 6):

```console title="Output" linenums="6"
[a3/7be2fa] sayHello | 1 of 1 ✔
```

Això ens indica que el procés `sayHello` s'ha executat correctament una vegada (`1 of 1 ✔`).

Això és genial, però potser us esteu preguntant: on és la sortida?

### 1.2. Trobar el fitxer de sortida al directori `results`

Aquest workflow està configurat per publicar la seva sortida a un directori anomenat `results`.
Si mireu el vostre directori actual, veureu que quan heu executat el workflow, Nextflow ha creat un nou directori anomenat `results`, que conté un fitxer anomenat `output.txt`.

```console title="results/" linenums="1"
results
└── output.txt
```

Obriu el fitxer; el contingut hauria de coincidir amb la salutació que heu especificat a la línia de comandes.

<details>
  <summary>Contingut del fitxer</summary>

```console title="results/output.txt" linenums="1"
Hello World!
```

</details>

Això és genial, el nostre workflow ha fet el que havia de fer!

Tanmateix, tingueu en compte que el resultat 'publicat' és una còpia (o en alguns casos un enllaç simbòlic) de la sortida real produïda per Nextflow quan ha executat el workflow.

Així que ara, mirarem sota el capó per veure on Nextflow ha executat realment el treball.

!!! warning "Advertència"

    No tots els workflows estaran configurats per publicar sortides a un directori de resultats, i/o el nom del directori pot ser diferent.
    Una mica més endavant en aquesta secció, us mostrarem com esbrinar on s'especifica aquest comportament.

### 1.3. Trobar la sortida original i els registres al directori `work/`

Quan executeu un workflow, Nextflow crea un 'directori de tasca' distint per a cada invocació de cada procés del workflow (=cada pas del pipeline).
Per a cadascun, prepararà les entrades necessàries, executarà les instruccions rellevants i escriurà sortides i fitxers de registre dins d'aquest únic directori, que es nomena automàticament utilitzant un hash per tal de fer-lo únic.

Tots aquests directoris de tasca es trobaran sota un directori anomenat `work` dins del vostre directori actual (on esteu executant la comanda).

Això pot sonar confús, així que vegem com es veu a la pràctica.

Tornant a la sortida de consola del workflow que hem executat anteriorment, teníem aquesta línia:

```console title="Excerpt of command output" linenums="6"
[a3/7be2fa] sayHello | 1 of 1 ✔
```

Veieu com la línia comença amb `[a3/7be2fa]`?
Aquesta és una forma truncada del camí del directori de tasca per a aquesta crida de procés, i us indica on trobar la sortida de la crida del procés `sayHello` dins del camí del directori `work/`.

Podeu trobar el camí complet escrivint la comanda següent (substituint `a3/7be2fa` pel que veieu al vostre propi terminal) i prement la tecla de tabulació per autocompletar el camí o afegint un asterisc:

```bash
tree work/a3/7be2fa*
```

Això hauria de produir el camí complet del directori: `work/a3/7be2fa7be2fad5e71e5f49998f795677fd68`

Vegem què hi ha allà.

!!! Tip "Consell"

    Si navegueu pel contingut del subdirectori de tasca a l'explorador de fitxers de VSCode, veureu tots els fitxers immediatament.
    Tanmateix, els fitxers de registre estan configurats per ser invisibles al terminal, així que si voleu utilitzar `ls` o `tree` per veure'ls, haureu d'establir l'opció rellevant per mostrar fitxers invisibles.

    ```bash
    tree -a work
    ```

Els noms exactes dels subdirectoris seran diferents al vostre sistema.

<details>
  <summary>Contingut del directori</summary>

```console title="work/"
work
└── a3
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

</details>

Hauríeu de reconèixer immediatament el fitxer `output.txt`, que de fet és la sortida original del procés `sayHello` que es va publicar al directori `results`.
Si l'obriu, trobareu la salutació `Hello World!` de nou.

<details>
  <summary>Contingut del fitxer output.txt</summary>

```console title="work/a3/7be2fa7be2fad5e71e5f49998f795677fd68/output.txt" linenums="1"
Hello World!
```

</details>

Així que, què passa amb tots aquests altres fitxers?

Aquests són els fitxers auxiliars i de registre que Nextflow va escriure com a part de l'execució de la tasca:

- **`.command.begin`**: Fitxer sentinella creat tan aviat com es llança la tasca.
- **`.command.err`**: Missatges d'error (`stderr`) emesos per la crida del procés
- **`.command.log`**: Sortida de registre completa emesa per la crida del procés
- **`.command.out`**: Sortida regular (`stdout`) per la crida del procés
- **`.command.run`**: Script complet executat per Nextflow per executar la crida del procés
- **`.command.sh`**: La comanda que realment va executar la crida del procés
- **`.exitcode`**: El codi de sortida resultant de la comanda

El fitxer `.command.sh` és especialment útil perquè us mostra la comanda principal que Nextflow va executar sense incloure tota la comptabilitat i la configuració de tasca/entorn.

<details>
  <summary>Contingut del fitxer</summary>

```console title="work/a3/7be2fa7be2fad5e71e5f49998f795677fd68/.command.sh" linenums="1"
#!/bin/bash -ue
echo 'Hello World!' > output.txt

```

</details>

!!! Tip "Consell"

    Quan alguna cosa va malament i necessiteu solucionar problemes sobre què ha passat, pot ser útil mirar l'script `command.sh` per comprovar exactament quina comanda va compondre Nextflow basant-se en les instruccions del workflow, la interpolació de variables, etc.

### 1.4. Exercici opcional: tornar a executar amb diferents salutacions

Proveu de tornar a executar el workflow diverses vegades amb diferents valors per a l'argument `--greeting`, després mireu tant el contingut del directori `results/` com els directoris de tasca.

Observeu com les sortides i els registres dels directoris de tasca aïllats es preserven, mentre que el contingut del directori `results` es sobreescriu per la sortida d'execucions posteriors.

### Conclusió

Sabeu com executar un script simple de Nextflow, monitoritzar la seva execució i trobar les seves sortides.

### Què segueix?

Apreneu a llegir un script bàsic de Nextflow i identificar com els seus components es relacionen amb la seva funcionalitat.

---

## 2. Examinar l'script inicial del workflow Hello World

El que hem fet allà bàsicament ha estat tractar l'script del workflow com una caixa negra.
Ara que hem vist què fa, obrim la caixa i mirem a dins.

_L'objectiu aquí no és memoritzar la sintaxi del codi de Nextflow, sinó formar una intuïció bàsica de quins són els components principals i com estan organitzats._

### 2.1. Examinar l'estructura general del codi

Obrim l'script `hello-world.nf` al panell de l'editor.

<details>
  <summary>Codi</summary>

```groovy title="hello-world.nf" linenums="1"
#!/usr/bin/env nextflow

/*
 * Use echo to print a greeting to a file
 */
process sayHello {

    publishDir 'results', mode: 'copy'

    input:
    val greeting

    output:
    path 'output.txt'

    script:
    """
    echo '$greeting' > output.txt
    """
}

workflow {

    // emet una salutacio
    sayHello(params.greeting)
}
```

</details>

Un script de Nextflow implica dos tipus principals de components bàsics: un o més **processos**, i el **workflow** en si mateix.
Cada **procés** descriu quines operacions ha d'acomplir el pas corresponent del pipeline, mentre que el **workflow** descriu la lògica de flux de dades que connecta els diversos passos.

Vegem més de prop primer el bloc de **procés**, després mirarem el bloc de **workflow**.

### 2.2. La definició de `process`

El primer bloc de codi descriu un **procés**.
La definició del procés comença amb la paraula clau `process`, seguida del nom del procés i finalment el cos del procés delimitat per claus.
El cos del procés ha de contenir un bloc de script que especifica la comanda a executar, que pot ser qualsevol cosa que poguéssiu executar en un terminal de línia de comandes.

Aquí tenim un **procés** anomenat `sayHello` que pren una variable d'**entrada** anomenada `greeting` i escriu la seva **sortida** a un fitxer anomenat `output.txt`.

<details>
  <summary>Codi</summary>

```groovy title="hello-world.nf" linenums="3"
/*
 * Use echo to print a greeting to a file
 */
process sayHello {

    publishDir 'results', mode: 'copy'

    input:
    val greeting

    output:
    path 'output.txt'

    script:
    """
    echo '$greeting' > output.txt
    """
}
```

</details>

Aquesta és una definició de procés molt mínima que només conté una definició d'`input`, una definició d'`output` i l'`script` a executar.

La definició d'`input` inclou el qualificador `val`, que indica a Nextflow que esperi un valor d'algun tipus (pot ser una cadena de text, un número, el que sigui).

La definició d'`output` inclou el qualificador `path`, que indica a Nextflow que això s'ha de gestionar com un camí (inclou tant camins de directori com fitxers).

!!! Tip "Consell"

    La definició de sortida no _determina_ quina sortida es crearà.
    Simplement _declara_ on trobar els fitxers de sortida esperats, perquè Nextflow pugui buscar-los una vegada l'execució estigui completa.

    Això és necessari per verificar que la comanda s'ha executat correctament i per passar la sortida a processos posteriors si cal.
    La sortida produïda que no coincideixi amb el que es declara al bloc de sortida no es passarà als processos posteriors.

En un pipeline del món real, un procés normalment conté informació addicional com ara directives de procés, que introduirem d'aquí a una estona.

### 2.3. La definició de `workflow`

El segon bloc de codi descriu el **workflow** en si mateix.
La definició del workflow comença amb la paraula clau `workflow`, seguida d'un nom opcional, després el cos del workflow delimitat per claus.

Aquí tenim un **workflow** que consisteix en una crida al procés `sayHello`, que pren una entrada, `params.greeting`, que conté el valor que hem donat al paràmetre `--greeting`.

```groovy title="hello-world.nf" linenums="22"
workflow {

    // emet una salutacio
    sayHello(params.greeting)
}
```

Aquesta és una definició de **workflow** molt mínima.
En un pipeline del món real, el workflow normalment conté múltiples crides a **processos** connectats per **canals**, i pot haver-hi valors per defecte configurats per a les entrades de variables.

Veurem això en acció quan executem nf-core/molkart a la Part 2 del curs.

### 2.4. El sistema `params` de paràmetres de línia de comandes

El `params.greeting` que proporcionem a la crida del procés `sayHello()` és un fragment de codi de Nextflow interessant i val la pena dedicar-hi un minut extra.

Com s'ha esmentat anteriorment, així és com passem el valor del paràmetre de línia de comandes `--greeting` a la crida del procés `sayHello()`.
De fet, simplement declarant `params.someParameterName` ens permetrà donar al workflow un paràmetre anomenat `--someParameterName` des de la línia de comandes.

!!! Tip "Consell"

    Aquests paràmetres de workflow declarats utilitzant el sistema `params` sempre prenen dos guions (`--`).
    Això els distingeix dels paràmetres de nivell de Nextflow, que només prenen un guió (`-`).

### Conclusió

Ara sabeu com està estructurat un workflow simple de Nextflow, i com els components bàsics es relacionen amb la seva funcionalitat.

### Què segueix?

Apreneu a gestionar les vostres execucions de workflow de manera convenient.

---

## 3. Gestionar execucions de workflow

Saber com llançar workflows i recuperar sortides és genial, però aviat trobareu que hi ha alguns altres aspectes de la gestió de workflows que us facilitaran la vida.

Aquí us mostrem com aprofitar la funció `resume` per quan necessiteu tornar a llançar el mateix workflow, com inspeccionar els registres d'execució amb `nextflow log`, i com eliminar directoris de treball antics amb `nextflow clean`.

### 3.1. Tornar a llançar un workflow amb `-resume`

De vegades, voldreu tornar a executar un pipeline que ja heu llançat anteriorment sense refer cap treball que ja s'hagi completat correctament.

Nextflow té una opció anomenada `-resume` que us permet fer això.
Específicament, en aquest mode, qualsevol procés que ja s'hagi executat amb exactament el mateix codi, configuració i entrades es saltarà.
Això significa que Nextflow només executarà processos que hàgiu afegit o modificat des de l'última execució, o als quals estigueu proporcionant noves configuracions o entrades.

Hi ha dos avantatges clau de fer això:

- Si esteu en procés de desenvolupar un pipeline, podeu iterar més ràpidament ja que només heu d'executar els processos en els quals esteu treballant activament per provar els vostres canvis.
- Si esteu executant un pipeline en producció i alguna cosa va malament, en molts casos podeu solucionar el problema i tornar a llançar el pipeline, i es reprendrà l'execució des del punt de fallada, cosa que us pot estalviar molt de temps i càlcul.

Per utilitzar-lo, simplement afegiu `-resume` a la vostra comanda i executeu-la:

```bash
nextflow run hello-world.nf --greeting 'Hello World!' -resume
```

??? success "Sortida de la comanda"

    ```console
    N E X T F L O W   ~  version 25.04.3

    Launching `hello-world.nf` [tiny_noyce] DSL2 - revision: c33d41f479

    [a3/7be2fa] process > sayHello [100%] 1 of 1, cached: 1 ✔
    ```

Busqueu la part `cached:` que s'ha afegit a la línia d'estat del procés (línia 5), que significa que Nextflow ha reconegut que ja ha fet aquest treball i simplement ha reutilitzat el resultat de l'execució anterior correcta.

També podeu veure que el hash del subdirectori de treball és el mateix que a l'execució anterior.
Nextflow literalment us està assenyalant l'execució anterior i dient "Ja ho vaig fer allà."

!!! Tip "Consell"

    Quan torneu a executar un pipeline amb `resume`, Nextflow no sobreescriu cap fitxer escrit a un directori `publishDir` per cap crida de procés que s'hagi executat anteriorment amb èxit.

### 3.2. Inspeccionar el registre d'execucions passades

Cada vegada que llanceu un workflow de nextflow, s'escriu una línia a un fitxer de registre anomenat `history`, sota un directori ocult anomenat `.nextflow` al directori de treball actual.

Una manera més convenient d'accedir a aquesta informació és utilitzar la comanda `nextflow log`.

```bash
nextflow log
```

Això mostrarà el contingut del fitxer de registre al terminal, mostrant-vos la marca de temps, el nom de l'execució, l'estat i la línia de comandes completa per a cada execució de Nextflow que s'hagi llançat des del directori de treball actual.

### 3.3. Eliminar directoris de treball antics

Durant el procés de desenvolupament, normalment executareu els vostres esborranys de pipelines un gran nombre de vegades, cosa que pot portar a una acumulació de molts fitxers a través de molts subdirectoris.
Com que els subdirectoris es nomenen aleatòriament, és difícil dir pels seus noms quines són execucions més antigues vs. més recents.

Nextflow inclou una subcomanda `clean` convenient que pot eliminar automàticament els subdirectoris de treball per a execucions passades que ja no us importen, amb diverses [opcions](https://www.nextflow.io/docs/latest/reference/cli.html#clean) per controlar què s'eliminarà.

Podeu utilitzar el registre de Nextflow per buscar una execució basant-vos en la seva marca de temps i/o línia de comandes, després utilitzar `nextflow clean -before <run_name> -f` per eliminar directoris de treball d'execucions anteriors.

!!! Warning "Advertència"

    Eliminar subdirectoris de treball d'execucions passades els elimina de la memòria cau de Nextflow i elimina qualsevol sortida que s'hagi emmagatzemat en aquests directoris.
    Això significa que trenca la capacitat de Nextflow de reprendre l'execució sense tornar a executar els processos corresponents.

    Sou responsables de desar qualsevol sortida que us importi o en la qual tingueu previst confiar! Si esteu utilitzant la directiva `publishDir` per a aquest propòsit, assegureu-vos d'utilitzar el mode `copy`, no el mode `symlink`.

### Conclusió

Sabeu com tornar a llançar un pipeline sense repetir passos que ja s'han executat d'una manera idèntica, inspeccionar el registre d'execució i utilitzar la comanda `nextflow clean` per netejar directoris de treball antics.

### Què segueix?

Ara que enteneu les operacions bàsiques de Nextflow, esteu preparats per executar un pipeline real de bioimatge amb nf-core/molkart.
