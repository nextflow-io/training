# Part 2: Hello Channels

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Traducció assistida per IA - [més informació i suggeriments](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

<div class="video-wrapper">
  <iframe width="560" height="315" src="https://www.youtube.com/embed/yDR66fzAMOg?si=y8lAedhEHWaTV4zd&amp;list=PLPZ8WHdZGxmWKozQuzr27jyMGqp9kElVK&amp;cc_load_policy=1" title="YouTube video player" frameborder="0" allow="accelerometer; autoplay; clipboard-write; encrypted-media; gyroscope; picture-in-picture; web-share" referrerpolicy="strict-origin-when-cross-origin" allowfullscreen></iframe>
</div>

/// caption
:fontawesome-brands-youtube:{ .youtube } Consulteu [la llista de reproducció completa](https://youtube.com/playlist?list=PLPZ8WHdZGxmWKozQuzr27jyMGqp9kElVK&si=eF7cLR62goy-lc6n) al canal de YouTube de Nextflow.

:green_book: La transcripció del vídeo està disponible [aquí](./transcripts/02_hello_channels.md).
///

A la Part 1 d'aquest curs (Hello World), us vam mostrar com proporcionar una entrada variable a un procés proporcionant l'entrada directament a la crida del procés: `sayHello(params.input)`.
Aquest va ser un enfocament deliberadament simplificat.
A la pràctica, aquest enfocament té limitacions importants; concretament, només funciona per a casos molt simples on només volem executar el procés una vegada, amb un únic valor.
En la majoria de casos d'ús realistes de workflows, volem processar múltiples valors (dades experimentals per a múltiples mostres, per exemple), així que necessitem una manera més sofisticada de gestionar les entrades.

Per això existeixen els [**canals**](https://nextflow.io/docs/latest/channel.html) de Nextflow.
Els canals són cues dissenyades per gestionar entrades de manera eficient i transportar-les d'un pas a un altre en workflows de múltiples passos, alhora que proporcionen paral·lelisme integrat i molts beneficis addicionals.

En aquesta part del curs, aprendreu com utilitzar un canal per gestionar múltiples entrades de diverses fonts diferents.
També aprendreu a utilitzar [**operadors**](https://nextflow.io/docs/latest/reference/operator.html) per transformar el contingut dels canals segons sigui necessari.

??? info "Com començar des d'aquesta secció"

    Aquesta secció del curs assumeix que heu completat la Part 1 del curs [Hello Nextflow](./index.md), però si us sentiu còmodes amb els conceptes bàsics tractats en aquella secció, podeu començar des d'aquí sense fer res especial.

---

## 0. Escalfament: Executeu `hello-channels.nf`

Utilitzarem l'script de workflow `hello-channels.nf` com a punt de partida.
És equivalent a l'script produït en treballar la Part 1 d'aquest curs de formació, excepte que hem canviat la destinació de sortida:

```groovy title="hello-channels.nf" linenums="37" hl_lines="3"
output {
    first_output {
        path 'hello_channels'
        mode 'copy'
    }
}
```

Només per assegurar-nos que tot funciona, executeu l'script una vegada abans de fer cap canvi:

```bash
nextflow run hello-channels.nf --input 'Hello Channels!'
```

??? success "Sortida de la comanda"

    ```console
     N E X T F L O W   ~  version 25.10.2

    Launching `hello-channels.nf` [wise_jennings] DSL2 - revision: b24f4902d6

    executor >  local (1)
    [6f/824bc1] process > sayHello [100%] 1 of 1 ✔
    ```

Com abans, trobareu el fitxer de sortida anomenat `output.txt` al directori `results/hello_channels` (tal com s'especifica al bloc `output` de l'script de workflow, mostrat més amunt).

??? abstract "Contingut del directori"

    ```console title="results/hello_channels" hl_lines="2-3"
    results
    ├── hello_channels
    │   └── output.txt
    ├── hello_world
    │   └── output.txt
    └── output.txt -> /workspaces/training/hello-nextflow/work/8c/79499c11beea6e9d43605141f2817f/output.txt
    ```

??? abstract "Contingut del fitxer"

    ```console title="results/hello_channels/output.txt"
    Hello Channels!
    ```

Si això us ha funcionat, esteu preparats per aprendre sobre els canals.

---

## 1. Proporcioneu entrades variables mitjançant un canal explícitament

Crearem un **canal** per passar l'entrada variable al procés `sayHello()` en lloc de confiar en la gestió implícita, que té certes limitacions.

### 1.1. Creeu un canal d'entrada

Hi ha una varietat de [**constructors de canals**](https://nextflow.io/docs/latest/reference/channel.html) que podem utilitzar per configurar un canal.
Per mantenir les coses simples de moment, utilitzarem el constructor de canals més bàsic, anomenat [`channel.of`](https://nextflow.io/docs/latest/reference/channel.html#of), que crearà un canal que conté un únic valor.
Funcionalment això serà similar a com ho teníem configurat abans, però en lloc de fer que Nextflow creï un canal implícitament, ara ho estem fent explícitament.

Aquesta és la línia de codi que utilitzarem:

```console title="Syntax"
greeting_ch = channel.of('Hello Channels!')
```

Això crea un canal anomenat `greeting_ch` utilitzant el constructor de canals `channel.of()`, que configura un canal de cua simple, i carrega la cadena `'Hello Channels!'` per utilitzar-la com a valor de salutació.

<figure class="excalidraw">
--8<-- "docs/en/docs/hello_nextflow/img/hello-pipeline-channel.svg"
</figure>

!!! note "Nota"

    Estem tornant temporalment a cadenes codificades en lloc d'utilitzar un paràmetre de CLI per motius de llegibilitat. Tornarem a utilitzar paràmetres de CLI un cop hàgim cobert què està passant al nivell del canal.

Al bloc workflow, afegiu el codi del constructor de canals:

=== "Després"

    ```groovy title="hello-channels.nf" linenums="27" hl_lines="4 5"
    workflow {

        main:
        // crea un canal per a les entrades
        greeting_ch = channel.of('Hello Channels!')
        // emet una salutacio
        sayHello(params.input)

        publish:
        first_output = sayHello.out
    }
    ```

=== "Abans"

    ```groovy title="hello-channels.nf" linenums="27"
    workflow {

        main:
        // emet una salutacio
        sayHello(params.input)

        publish:
        first_output = sayHello.out
    }
    ```

Això encara no és funcional ja que encara no hem canviat l'entrada a la crida del procés.

### 1.2. Afegiu el canal com a entrada a la crida del procés

Ara necessitem connectar el nostre canal recentment creat a la crida del procés `sayHello()`, substituint el paràmetre de CLI que estàvem proporcionant directament abans.

Al bloc workflow, feu el següent canvi de codi:

=== "Després"

    ```groovy title="hello-channels.nf" linenums="27" hl_lines="7"
    workflow {

        main:
        // crea un canal per a les entrades
        greeting_ch = channel.of('Hello Channels!')
        // emet una salutacio
        sayHello(greeting_ch)

        publish:
        first_output = sayHello.out
    }
    ```

=== "Abans"

    ```groovy title="hello-channels.nf" linenums="27" hl_lines="7"
    workflow {

        main:
        // crea un canal per a les entrades
        greeting_ch = channel.of('Hello Channels!')
        // emet una salutacio
        sayHello(params.input)

        publish:
        first_output = sayHello.out
    }
    ```

Això indica a Nextflow que executi el procés `sayHello` amb el contingut del canal `greeting_ch`.

Ara el nostre workflow és adequadament funcional; és l'equivalent explícit d'escriure `sayHello('Hello Channels!')`.

### 1.3. Executeu el workflow

Executem-ho!

```bash
nextflow run hello-channels.nf
```

??? success "Sortida de la comanda"

    ```console
     N E X T F L O W   ~  version 25.10.2

    Launching `hello-channels.nf` [fabulous_crick] DSL2 - revision: 23e20f76e8

    executor >  local (1)
    [c0/4f1872] process > sayHello (1) [100%] 1 of 1 ✔
    ```

Si heu fet les dues edicions correctament, hauríeu d'obtenir una execució exitosa.
Podeu comprovar el directori de resultats per assegurar-vos que el resultat és encara el mateix que abans.

??? abstract "Contingut del fitxer"

    ```console title="results/hello_channels/output.txt"
    Hello Channels!
    ```

Així que hem augmentat la flexibilitat del nostre workflow mentre aconseguim el mateix resultat final.
Això pot semblar que estem escrivint més codi sense cap benefici tangible, però el valor quedarà clar tan aviat com comencem a gestionar més entrades.

Com a previsualització d'això, vegem una cosa més abans de continuar: un petit però convenient benefici d'utilitzar un canal explícit per gestionar l'entrada de dades.

### 1.4. Utilitzeu `view()` per inspeccionar el contingut del canal

Els canals de Nextflow estan construïts d'una manera que ens permet operar sobre el seu contingut utilitzant operadors, que tractarem en detall més endavant en aquest capítol.

De moment, només us mostrarem com utilitzar un operador súper simple anomenat [`view()`](https://www.nextflow.io/docs/latest/reference/operator.html#view) per inspeccionar el contingut d'un canal.
Podeu pensar en `view()` com una eina de depuració, com una instrucció `print()` en Python, o el seu equivalent en altres llenguatges.

Afegiu aquesta petita línia al bloc workflow:

=== "Després"

    ```groovy title="hello-channels.nf" linenums="27" hl_lines="7"
    workflow {

        main:
        // crea un canal per a les entrades
        greeting_ch = channel.of('Hello Channels!')
                             .view()
        // emet una salutacio
        sayHello(greeting_ch)

        publish:
        first_output = sayHello.out
    }
    ```

=== "Abans"

    ```groovy title="hello-channels.nf" linenums="27"
    workflow {

        main:
        // crea un canal per a les entrades
        greeting_ch = channel.of('Hello Channels!')
        // emet una salutacio
        sayHello(greeting_ch)

        publish:
        first_output = sayHello.out
    }
    ```

La quantitat exacta d'espais no importa sempre que sigui un múltiple de 4; només estem intentant alinear l'inici de la instrucció `.view()` amb la part `.of()` de la construcció del canal.

Ara executeu el workflow de nou:

```bash
nextflow run hello-channels.nf
```

??? success "Sortida de la comanda"

    ```console hl_lines="7"
     N E X T F L O W   ~  version 25.10.2

    Launching `hello-channels.nf` [scruffy_shaw] DSL2 - revision: 2ede41e14a

    executor >  local (1)
    [ef/f7e40a] sayHello (1) [100%] 1 of 1 ✔
    Hello Channels!
    ```

Com podeu veure, això mostra el contingut del canal a la consola.
Aquí només tenim un element, però quan comencem a carregar múltiples valors al canal a la següent secció, veureu que està configurat per mostrar un element per línia.

### Conclusió

Sabeu com utilitzar un constructor de canals bàsic per proporcionar una entrada a un procés.

### Què segueix?

Apreneu com utilitzar canals per fer que el workflow iterin sobre múltiples valors d'entrada.

---

## 2. Modifiqueu el workflow per executar-se amb múltiples valors d'entrada

Els workflows normalment s'executen amb lots d'entrades que estan destinades a ser processades en bloc, així que volem actualitzar el workflow per acceptar múltiples valors d'entrada.

### 2.1. Carregueu múltiples salutacions al canal d'entrada

Convenientment, el constructor de canals `channel.of()` que hem estat utilitzant està molt content d'acceptar més d'un valor, així que no necessitem modificar-lo en absolut.
Simplement podem carregar múltiples valors al canal.

Fem-los `'Hello'`, `'Bonjour'` i `'Holà'`.

#### 2.1.1. Afegiu més salutacions

Al bloc workflow, feu el següent canvi de codi:

=== "Després"

    ```groovy title="hello-channels.nf" linenums="30" hl_lines="2"
    // crea un canal per a les entrades
    greeting_ch = channel.of('Hello','Bonjour','Holà')
                         .view()
    ```

=== "Abans"

    ```groovy title="hello-channels.nf" linenums="30" hl_lines="2"
    // crea un canal per a les entrades
    greeting_ch = channel.of('Hello Channels')
                         .view()
    ```

La documentació ens diu que això hauria de funcionar. Pot ser realment tan simple?

#### 2.1.2. Executeu la comanda i mireu la sortida del registre

Provem-ho.

```bash
nextflow run hello-channels.nf
```

??? success "Sortida de la comanda"

    ```console hl_lines="6"
     N E X T F L O W   ~  version 25.10.2

    Launching `hello-channels.nf` [amazing_crick] DSL2 - revision: 59a9a5888a

    executor >  local (3)
    [f4/c9962c] process > sayHello (1) [100%] 3 of 3 ✔
    Hello
    Bonjour
    Holà
    ```

Certament sembla que s'ha executat correctament.
El monitor d'execució mostra que es van fer `3 de 3` crides al procés `sayHello`, i veiem les tres salutacions enumerades per la instrucció `view()`, una per línia com es va prometre.

No obstant això, encara només hi ha una sortida al directori de resultats:

??? abstract "Contingut del directori"

    ```console title="results/hello_channels" hl_lines="3"
    results
    ├── hello_channels
    │   └── output.txt
    ├── hello_world
    │   └── output.txt
    └── output.txt -> /workspaces/training/hello-nextflow/work/8c/79499c11beea6e9d43605141f2817f/output.txt
    ```

??? abstract "Contingut del fitxer"

    ```console title="results/hello_channels/output.txt"
    Holà
    ```

Hauríeu de veure una de les tres salutacions allà, encara que la que heu obtingut pot ser diferent del que es mostra aquí.
Podeu pensar per què podria ser això?

<figure class="excalidraw">
--8<-- "docs/en/docs/hello_nextflow/img/hello-pipeline-channel-multi.svg"
</figure>

_Al diagrama, el canal es representa en verd, i l'ordre dels elements es representa com boles en una canonada: el primer carregat està a la dreta, després el segon al mig, després el tercer està a l'esquerra._

Mirant enrere al monitor d'execució, només ens va donar un camí de subdirectori (`f4/c9962c`).
Donem-hi una ullada.

??? abstract "Contingut del directori"

    ```console hl_lines="9"
    work/f4/c9962ce91ef87480babcb86b2b9042/
    ├── .command.begin
    ├── .command.err
    ├── .command.log
    ├── .command.out
    ├── .command.run
    ├── .command.sh
    ├── .exitcode
    └── output.txt
    ```

??? abstract "Contingut del fitxer"

    ```console title="work/f4/c9962ce91ef87480babcb86b2b9042/output.txt"
    Hello
    ```

Aquesta ni tan sols és la mateixa salutació que vam obtenir al directori de resultats! Què està passant?

En aquest punt, hem de dir-vos que per defecte, el sistema de registre ANSI escriu el registre de múltiples crides al mateix procés a la mateixa línia.
Així que l'estat de les tres crides al procés sayHello() estan aterrant al mateix lloc.

Afortunadament, podem desactivar aquest comportament per veure la llista completa de crides de procés.

#### 2.1.3. Executeu la comanda de nou amb l'opció `-ansi-log false`

Per expandir el registre per mostrar una línia per crida de procés, afegiu `-ansi-log false` a la comanda.

```bash
nextflow run hello-channels.nf -ansi-log false
```

??? success "Sortida de la comanda"

    ```console
     N E X T F L O W  ~  version 25.10.2
    Launching `hello-channels.nf` [desperate_monod] DSL2 - revision: 59a9a5888a
    Hello
    Bonjour
    Holà
    [23/871c7e] Submitted process > sayHello (2)
    [7f/21e2c2] Submitted process > sayHello (1)
    [f4/ea10a6] Submitted process > sayHello (3)
    ```

Aquesta vegada veiem les tres execucions de procés i els seus subdirectoris de treball associats llistats a la sortida.

Això és molt millor, almenys per a un workflow simple.
Per a un workflow complex, o un gran nombre d'entrades, tenir la llista completa mostrada al terminal seria una mica aclaparador.
Per això `-ansi-log false` no és el comportament per defecte.

!!! tip "Consell"

    La manera com es reporta l'estat és una mica diferent entre els dos modes de registre.
    En el mode condensat, Nextflow informa si les crides es van completar amb èxit o no.
    En aquest mode expandit, només informa que es van enviar.

De totes maneres, ara que tenim els subdirectoris de cada crida de procés, podem buscar els seus registres i sortides.

??? abstract "Contingut del directori"

    ```console
    work/23/871c7ec3642a898ecd5e6090d21300/
    ├── .command.begin
    ├── .command.err
    ├── .command.log
    ├── .command.out
    ├── .command.run
    ├── .command.sh
    ├── .exitcode
    └── output.txt
    ```

    ```console
    work/7f/21e2c2f3cc8833ef3858b236e5575c/
    ├── .command.begin
    ├── .command.err
    ├── .command.log
    ├── .command.out
    ├── .command.run
    ├── .command.sh
    ├── .exitcode
    └── output.txt
    ```

    ```console
    work/f4/ea10a680d5687596d3eaa3fcf69272/
    ├── .command.begin
    ├── .command.err
    ├── .command.log
    ├── .command.out
    ├── .command.run
    ├── .command.sh
    ├── .exitcode
    └── output.txt
    ```

??? abstract "Contingut del fitxer"

    ```txt title="work/23/871c7ec3642a898ecd5e6090d21300/output.txt"
    Bonjour
    ```

    ```txt title="work/7f/21e2c2f3cc8833ef3858b236e5575c/output.txt"
    Hello
    ```

    ```txt title="work/f4/ea10a680d5687596d3eaa3fcf69272/output.txt"
    Holà
    ```

Això mostra que els tres processos s'han executat amb èxit (visca).

Dit això, encara tenim el problema que només hi ha un fitxer de sortida al directori de resultats.

Potser recordeu que vam codificar el nom del fitxer de sortida per al procés `sayHello`, així que les tres crides van produir un fitxer anomenat `output.txt`.

<figure class="excalidraw">
--8<-- "docs/en/docs/hello_nextflow/img/hello-channels-task-dirs.svg"
</figure>

Mentre els fitxers de sortida es quedin als subdirectoris de treball, aïllats dels altres processos, això està bé.
Però quan es publiquen al mateix directori de resultats, el que es va copiar primer és sobreescrit pel següent, i així successivament.

### 2.2. Assegureu-vos que els noms dels fitxers de sortida seran únics

Podem continuar publicant totes les sortides al mateix directori de resultats, però necessitem assegurar-nos que tindran noms únics.
Específicament, necessitem modificar el primer procés per generar un nom de fitxer dinàmicament perquè els noms de fitxer finals siguin únics.

Així que com fem que els noms de fitxer siguin únics?
Una manera comuna de fer-ho és utilitzar alguna peça única de metadades de les entrades (rebudes del canal d'entrada) com a part del nom del fitxer de sortida.
Aquí, per conveniència, només utilitzarem la salutació mateixa ja que és només una cadena curta, i la posarem davant del nom base del fitxer de sortida.

<figure class="excalidraw">
--8<-- "docs/en/docs/hello_nextflow/img/hello-pipeline-channel-multi-unique.svg"
</figure>

#### 2.2.1. Construïu un nom de fitxer de sortida dinàmic

Al bloc de procés, feu els següents canvis de codi:

=== "Després"

    ```groovy title="hello-channels.nf" linenums="6" hl_lines="7 11"
    process sayHello {

        input:
        val greeting

        output:
        path "${greeting}-output.txt"

        script:
        """
        echo '${greeting}' > '${greeting}-output.txt'
        """
    }
    ```

=== "Abans"

    ```groovy title="hello-channels.nf" linenums="6" hl_lines="7 11"
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

Assegureu-vos de substituir `output.txt` tant a la definició de sortida com al bloc de comandes `script:`.

!!! tip "Consell"

    A la definició de sortida, HAURIEU d'utilitzar cometes dobles al voltant de l'expressió del nom de fitxer de sortida (NO cometes simples), altrament fallarà.

Això hauria de produir un nom de fitxer de sortida únic cada vegada que es crida el procés, de manera que es pugui distingir de les sortides d'altres crides al mateix procés al directori de sortida.

#### 2.2.2. Executeu el workflow

Executem-ho. Tingueu en compte que hem tornat a executar amb la configuració de registre ANSI per defecte.

```bash
nextflow run hello-channels.nf
```

??? success "Sortida de la comanda"

    ```console
     N E X T F L O W   ~  version 25.10.2

    Launching `hello-channels.nf` [sharp_minsky] DSL2 - revision: 16a291febe

    executor >  local (3)
    [e8/33ee64] sayHello (2) [100%] 3 of 3 ✔
    Hello
    Bonjour
    Holà
    ```

Tornant a la vista de resum, la sortida es resumeix en una línia de nou.
Doneu una ullada al directori `results` per veure si totes les salutacions de sortida hi són.

??? abstract "Contingut del directori"

    ```console
    results/hello_channels/
    ├── Bonjour-output.txt
    ├── Hello-output.txt
    ├── Holà-output.txt
    └── output.txt
    ```

Sí! I cadascun té el contingut esperat.

??? abstract "Contingut del fitxer"

    ```console title="Bonjour-output.txt"
    Bonjour
    ```

    ```console title="Hello-output.txt"
    Hello
    ```

    ```console title="Holà-output.txt"
    Holà
    ```

Èxit! Ara podem afegir tantes salutacions com vulguem sense preocupar-nos que els fitxers de sortida siguin sobreescrits.

!!! tip "Consell"

    A la pràctica, anomenar fitxers basant-se en les dades d'entrada mateixes és gairebé sempre impràctic.
    La millor manera de generar noms de fitxer dinàmics és passar metadades a un procés juntament amb els fitxers d'entrada.
    Les metadades normalment es proporcionen mitjançant un 'full de mostres' o equivalents.
    Aprendreu com fer-ho més endavant a la vostra formació de Nextflow (vegeu [Missió secundària de metadades](../side_quests/metadata.md)).

### Conclusió

Sabeu com alimentar múltiples elements d'entrada a través d'un canal.

### Què segueix?

Apreneu a utilitzar un operador per transformar el contingut d'un canal.

---

## 3. Proporcioneu múltiples entrades mitjançant un array

Acabem de mostrar-vos com gestionar múltiples elements d'entrada que estaven codificats directament al constructor de canals.
Què passa si volíem proporcionar aquestes múltiples entrades d'una manera diferent?

Per exemple, imagineu que configurem una variable d'entrada que conté un array d'elements com aquest:

`greetings_array = ['Hello','Bonjour','Holà']`

Podem carregar això al nostre canal de sortida i esperar que funcioni?

<figure class="excalidraw">
--8<-- "docs/en/docs/hello_nextflow/img/hello-pipeline-multi-inputs-array.svg"
</figure>

Descobrim-ho.

### 3.1. Proporcioneu un array de valors com a entrada al canal

El sentit comú suggereix que hauríem de poder simplement passar un array de valors en lloc d'un únic valor.
Provem-ho; necessitarem configurar la variable d'entrada i carregar-la al constructor de canals.

#### 3.1.1. Configureu la variable d'entrada

Prenem la variable `greetings_array` que acabem d'imaginar i fem-la realitat afegint-la al bloc workflow:

=== "Després"

    ```groovy title="hello-channels.nf" linenums="27" hl_lines="4 5"
    workflow {

        main:
        // declara un array de salutacions d'entrada
        greetings_array = ['Hello','Bonjour','Holà']
        // crea un canal per a les entrades
        greeting_ch = channel.of('Hello','Bonjour','Holà')
                             .view()
        // emet una salutacio
        sayHello(greeting_ch)

        publish:
        first_output = sayHello.out
    }
    ```

=== "Abans"

    ```groovy title="hello-channels.nf" linenums="27"
    workflow {

        main:
        // crea un canal per a les entrades
        greeting_ch = channel.of('Hello','Bonjour','Holà')
                             .view()
        // emet una salutacio
        sayHello(greeting_ch)

        publish:
        first_output = sayHello.out
    }
    ```

Això encara no és funcional, només hem afegit una declaració per a l'array.

#### 3.1.2. Establiu l'array de salutacions com a entrada al constructor de canals

Ara substituirem els valors `'Hello','Bonjour','Holà'` actualment codificats al constructor de canals amb el `greetings_array` que acabem de crear.

Al bloc workflow, feu el següent canvi:

=== "Després"

    ```groovy title="hello-channels.nf" linenums="27" hl_lines="7"
    workflow {

        main:
        // declara un array de salutacions d'entrada
        greetings_array = ['Hello','Bonjour','Holà']
        // crea un canal per a les entrades
        greeting_ch = channel.of(greetings_array)
                             .view()
        // emet una salutacio
        sayHello(greeting_ch)

        publish:
        first_output = sayHello.out
    }
    ```

=== "Abans"

    ```groovy title="hello-channels.nf" linenums="27" hl_lines="7"
    workflow {

        main:
        // declara un array de salutacions d'entrada
        greetings_array = ['Hello','Bonjour','Holà']
        // crea un canal per a les entrades
        greeting_ch = channel.of('Hello','Bonjour','Holà')
                             .view()
        // emet una salutacio
        sayHello(greeting_ch)

        publish:
        first_output = sayHello.out
    }
    ```

Això hauria de ser funcional ara.

#### 3.1.3. Executeu el workflow

Provem d'executar-ho:

```bash
nextflow run hello-channels.nf
```

??? failure "Sortida de la comanda"

    ```console hl_lines="7 11 16"
    N E X T F L O W   ~  version 25.10.2

    Launching `hello-channels.nf` [friendly_koch] DSL2 - revision: 97256837a7

    executor >  local (1)
    [a8/1f6ead] sayHello (1) | 0 of 1
    [Hello, Bonjour, Holà]
    ERROR ~ Error executing process > 'sayHello (1)'

    Caused by:
      Missing output file(s) `[Hello, Bonjour, Holà]-output.txt` expected by process `sayHello (1)`


    Command executed:

      echo '[Hello, Bonjour, Holà]' > '[Hello, Bonjour, Holà]-output.txt'

    Command exit status:
      0

    Command output:
      (empty)

    Work dir:
      /workspaces/training/hello-nextflow/work/a8/1f6ead5f3fa30a3c508e2e7cf83ffb

    Tip: you can replicate the issue by changing to the process work dir and entering the command `bash .command.run`

    -- Check '.nextflow.log' file for details
    ```

Oh no! Hi ha un error!

Mireu la sortida de `view()` i els missatges d'error.

Sembla que Nextflow va intentar executar una única crida de procés, utilitzant `[Hello, Bonjour, Holà]` com un únic valor de cadena, en lloc d'utilitzar les tres cadenes de l'array com a valors separats.

<figure class="excalidraw">
--8<-- "docs/en/docs/hello_nextflow/img/hello-channels-array-fail.svg"
</figure>

Així que és l''embalatge' el que està causant el problema.
Com fem que Nextflow desempaqueti l'array i carregui les cadenes individuals al canal?

### 3.2. Utilitzeu un operador per transformar el contingut del canal

Aquí és on entren en joc els [**operadors**](https://nextflow.io/docs/latest/reference/operator.html).
Ja heu utilitzat l'operador `.view()`, que només mira què hi ha allà.
Ara mirarem operadors que ens permeten actuar sobre el contingut d'un canal.

Si fullegeu la [llista d'operadors](https://nextflow.io/docs/latest/reference/operator.html) a la documentació de Nextflow, trobareu [`flatten()`](https://nextflow.io/docs/latest/reference/operator.html#flatten), que fa exactament el que necessitem: desempaquetar el contingut d'un array i emetre'ls com a elements individuals.

#### 3.2.1. Afegiu l'operador `flatten()`

Per aplicar l'operador `flatten()` al nostre canal d'entrada, l'afegim a la declaració del constructor de canals.

Al bloc workflow, feu el següent canvi de codi:

=== "Després"

    ```groovy title="hello-channels.nf" linenums="27" hl_lines="9"
    workflow {

        main:
        // declara un array de salutacions d'entrada
        greetings_array = ['Hello','Bonjour','Holà']
        // crea un canal per a les entrades
        greeting_ch = channel.of(greetings_array)
                             .view()
                             .flatten()
        // emet una salutacio
        sayHello(greeting_ch)

        publish:
        first_output = sayHello.out
    }
    ```

=== "Abans"

    ```groovy title="hello-channels.nf" linenums="27"
    workflow {

        main:
        // declara un array de salutacions d'entrada
        greetings_array = ['Hello','Bonjour','Holà']
        // crea un canal per a les entrades
        greeting_ch = channel.of(greetings_array)
                             .view()
        // emet una salutacio
        sayHello(greeting_ch)

        publish:
        first_output = sayHello.out
    }
    ```

Aquí hem afegit l'operador a la següent línia per llegibilitat, però podeu afegir operadors a la mateixa línia que el constructor de canals si ho preferiu, així:
`greeting_ch = channel.of(greetings_array).view().flatten()`

<figure class="excalidraw">
--8<-- "docs/en/docs/hello_nextflow/img/hello-channels-array-success.svg"
</figure>

#### 3.2.2. Refineu les instruccions `view()`

Podríem executar això immediatament per provar si funciona, però mentre hi som, refinarem com inspeccionem el contingut del canal.

Volem poder contrastar com es veu el contingut abans i després que s'apliqui l'operador `flatten()`, així que n'afegirem un segon, I afegirem una mica de codi per fer que s'etiquetin més clarament a la sortida.

Al bloc workflow, feu el següent canvi de codi:

=== "Després"

    ```groovy title="hello-channels.nf" linenums="27" hl_lines="8-10"
    workflow {

        main:
        // declara un array de salutacions d'entrada
        greetings_array = ['Hello','Bonjour','Holà']
        // crea un canal per a les entrades
        greeting_ch = channel.of(greetings_array)
                             .view { greeting -> "Before flatten: $greeting" }
                             .flatten()
                             .view { greeting -> "After flatten: $greeting" }
        // emet una salutacio
        sayHello(greeting_ch)

        publish:
        first_output = sayHello.out
    }
    ```

=== "Abans"

    ```groovy title="hello-channels.nf" linenums="27" hl_lines="8-9"
    workflow {

        main:
        // declara un array de salutacions d'entrada
        greetings_array = ['Hello','Bonjour','Holà']
        // crea un canal per a les entrades
        greeting_ch = channel.of(greetings_array)
                             .view()
                             .flatten()
        // emet una salutacio
        sayHello(greeting_ch)

        publish:
        first_output = sayHello.out
    }
    ```

Veieu que hem afegit una segona instrucció `.view`, i per a cadascuna d'elles, hem substituït els parèntesis buits (`()`) amb claus que contenen codi, com ara `{ greeting -> "Abans de flatten: $greeting" }`.

Aquestes s'anomenen _closures_. El codi que contenen s'executarà per a cada element del canal.
Definim una variable temporal per al valor intern, aquí anomenada `greeting` (però podria ser qualsevol nom arbitrari), que només s'utilitza dins de l'àmbit d'aquesta closure.

En aquest exemple, `$greeting` representa cada element individual carregat al canal.
Això resultarà en una sortida de consola ben etiquetada.

!!! info "Info"

    En alguns pipelines podeu veure una variable especial anomenada `$it` utilitzada dins de closures d'operadors.
    Aquesta és una variable _implícita_ que permet un accés abreujat a la variable interna,
    sense necessitat de definir-la amb un `->`.

    Preferim ser explícits per ajudar a la claredat del codi, per tant la sintaxi `$it` està desaconsellada i s'eliminarà gradualment del llenguatge Nextflow.

#### 3.2.3. Executeu el workflow

Finalment, podeu provar d'executar el workflow de nou!

```bash
nextflow run hello-channels.nf
```

??? success "Sortida de la comanda"

    ```console hl_lines="7-10"
     N E X T F L O W   ~  version 25.10.2

    Launching `hello-channels.nf` [sleepy_gutenberg] DSL2 - revision: 1db4f760ee

    executor >  local (3)
    [b1/6a1e15] sayHello (2) [100%] 3 of 3 ✔
    Before flatten: [Hello, Bonjour, Holà]
    After flatten: Hello
    After flatten: Bonjour
    After flatten: Holà
    ```

Aquesta vegada funciona I ens dóna la visió addicional de com es veu el contingut del canal abans i després d'executar l'operador `flatten()`.

- Una única instrucció `Abans de flatten:` perquè en aquell moment el canal conté un element, l'array original.
- Tres instruccions separades `Després de flatten:`, una per a cada salutació, que ara són elements individuals al canal.

Importantment, això significa que cada element ara pot ser processat separadament pel workflow.

!!! tip "Consell"

    És tècnicament possible aconseguir els mateixos resultats utilitzant un constructor de canals diferent, [`channel.fromList`](https://nextflow.io/docs/latest/reference/channel.html#fromlist), que inclou un pas de mapatge implícit en la seva operació.
    Aquí vam triar no utilitzar-lo per demostrar l'ús d'un operador en un cas d'ús simple.

### Conclusió

Sabeu com utilitzar un operador com `flatten()` per transformar el contingut d'un canal, i com utilitzar l'operador `view()` per inspeccionar el contingut del canal abans i després d'aplicar un operador.

### Què segueix?

Apreneu com fer que el workflow prengui un fitxer com a font dels seus valors d'entrada.

---

## 4. Llegiu valors d'entrada d'un fitxer CSV

Realisticament, gairebé mai començarem des d'un array de valors.
Molt probablement, tindrem un o més fitxers que contenen les dades que necessiten ser processades, en algun tipus de format estructurat.

Hem preparat un fitxer CSV anomenat `greetings.csv` que conté diverses salutacions d'entrada, imitant el tipus de dades en columnes que potser voldríeu processar en una anàlisi de dades real, emmagatzemat sota `data/`.
(Els números no són significatius, només hi són amb finalitats il·lustratives.)

```csv title="data/greetings.csv" linenums="1"
Hello,English,123
Bonjour,French,456
Holà,Spanish,789
```

La nostra següent tasca és adaptar el nostre workflow per llegir els valors d'aquest fitxer.

<figure class="excalidraw">
--8<-- "docs/en/docs/hello_nextflow/img/hello-pipeline-multi-inputs-csv.svg"
</figure>

Vegem com podem fer que això passi.

### 4.1. Modifiqueu l'script per esperar un fitxer CSV com a font de salutacions

Per començar, necessitarem fer dos canvis clau a l'script:

- Canviar el paràmetre d'entrada per apuntar al fitxer CSV
- Canviar el constructor de canals a un dissenyat per gestionar un fitxer

#### 4.1.1. Canvieu el paràmetre d'entrada per apuntar al fitxer CSV

Recordeu el paràmetre `params.input` que vam configurar a la Part 1?
L'actualitzarem per apuntar al fitxer CSV que conté les nostres salutacions.

Feu la següent edició a la declaració del paràmetre:

=== "Després"

    ```groovy title="hello-channels.nf" linenums="20" hl_lines="5"
    /*
    * Parametres del pipeline
    */
    params {
        input: Path = 'data/greetings.csv'
    }
    ```

=== "Abans"

    ```groovy title="hello-channels.nf" linenums="20" hl_lines="5"
    /*
     * Parametres del pipeline
     */
    input: String = 'Holà mundo!'
    ```

Això assumeix que el fitxer està col·locat amb el codi del workflow.
Aprendreu com tractar amb altres ubicacions de dades més endavant al vostre viatge amb Nextflow.

#### 4.1.2. Canvieu a un constructor de canals dissenyat per gestionar un fitxer

Ja que ara volem utilitzar un fitxer en lloc de cadenes simples com a entrada, no podem utilitzar el constructor de canals `channel.of()` d'abans.
Necessitem canviar a utilitzar un nou constructor de canals, [`channel.fromPath()`](https://nextflow.io/docs/latest/reference/channel.html#frompath), que té alguna funcionalitat integrada per gestionar camins de fitxers.

Al bloc workflow, feu el següent canvi de codi:

=== "Després"

    ```groovy title="hello-channels.nf" linenums="27" hl_lines="4-8"
    workflow {

        main:
        // crea un canal per a entrades des d'un fitxer CSV
        greeting_ch = channel.fromPath(params.input)
                             .view { greeting -> "Before flatten: $greeting" }
                             // .flatten()
                             // .view { greeting -> "Després de flatten: $greeting" }
        // emet una salutacio
        sayHello(greeting_ch)

        publish:
        first_output = sayHello.out
    }
    ```

=== "Abans"

    ```groovy title="hello-channels.nf" linenums="27" hl_lines="4-8"
    workflow {

        main:
        // declara un array de salutacions d'entrada
        greetings_array = ['Hello','Bonjour','Holà']
        // crea un canal per a les entrades
        greeting_ch = channel.of(greetings_array)
                             .view { greeting -> "Before flatten: $greeting" }
                             .flatten()
                             .view { greeting -> "After flatten: $greeting" }
        // emet una salutacio
        sayHello(greeting_ch)

        publish:
        first_output = sayHello.out
    }
    ```

Notareu que hem canviat l'entrada del canal de nou a `param.input`, i hem eliminat la declaració `greetings_array` ja que ja no la necessitarem.
També hem comentat el `flatten()` i la segona instrucció `view()`.

#### 4.1.3. Executeu el workflow

Provem d'executar el workflow amb el nou constructor de canals i el fitxer d'entrada.

```bash
nextflow run hello-channels.nf
```

??? failure "Sortida de la comanda"

    ```console hl_lines="5 6 9 14"
     N E X T F L O W   ~  version 25.10.2

    Launching `hello-channels.nf` [peaceful_poisson] DSL2 - revision: a286c08ad5

    [-        ] sayHello [  0%] 0 of 1
    Before flatten: /workspaces/training/hello-nextflow/data/greetings.csv
    ERROR ~ Error executing process > 'sayHello (1)'

    Caused by:
      File `/workspaces/training/hello-nextflow/data/greetings.csv-output.txt` is outside the scope of the process work directory: /workspaces/training/hello-nextflow/work/30/e610cb4ea5ae8693f456ac3329c92f


    Command executed:

      echo '/workspaces/training/hello-nextflow/data/greetings.csv' > '/workspaces/training/hello-nextflow/data/greetings.csv-output.txt'

    Command exit status:
      -

    Command output:
      (empty)

    Work dir:
      /workspaces/training/hello-nextflow/work/30/e610cb4ea5ae8693f456ac3329c92f

    Tip: when you have fixed the problem you can continue the execution adding the option `-resume` to the run command line

    -- Check '.nextflow.log' file for details
    ```

Oh no, no funciona. Doneu una ullada a l'inici de la sortida de consola i al missatge d'error.
La part `Command executed:` és especialment útil aquí.

Això pot semblar una mica familiar.
Sembla que Nextflow va intentar executar una única crida de procés utilitzant el camí del fitxer mateix com un valor de cadena.
Així que ha resolt el camí del fitxer correctament, però no ha analitzat realment el seu contingut, que és el que volíem.

Com fem que Nextflow obri el fitxer i carregui el seu contingut al canal?

Sembla que necessitem un altre [operador](https://nextflow.io/docs/latest/reference/operator.html)!

### 4.2. Utilitzeu l'operador `splitCsv()` per analitzar el fitxer

Mirant de nou la llista d'operadors, trobem [`splitCsv()`](https://nextflow.io/docs/latest/reference/operator.html#splitcsv), que està dissenyat per analitzar i dividir text amb format CSV.

#### 4.2.1. Apliqueu `splitCsv()` al canal

Per aplicar l'operador, l'afegim a la línia del constructor de canals com abans.

Al bloc workflow, feu el següent canvi de codi per substituir `flatten()` amb `splitCsv()` (sense comentar):

=== "Després"

    ```groovy title="hello-channels.nf" linenums="27" hl_lines="6-8"
    workflow {

        main:
        // crea un canal per a entrades des d'un fitxer CSV
        greeting_ch = channel.fromPath(params.input)
                             .view { csv -> "Before splitCsv: $csv" }
                             .splitCsv()
                             .view { csv -> "After splitCsv: $csv" }
        // emet una salutacio
        sayHello(greeting_ch)

        publish:
        first_output = sayHello.out
    }
    ```

=== "Abans"

    ```groovy title="hello-channels.nf" linenums="27" hl_lines="6-8"
    workflow {

        main:
        // crea un canal per a entrades des d'un fitxer CSV
        greeting_ch = channel.fromPath(params.input)
                             .view { greeting -> "Before flatten: $greeting" }
                             // .flatten()
                             // .view { greeting -> "Després de flatten: $greeting" }
        // emet una salutacio
        sayHello(greeting_ch)

        publish:
        first_output = sayHello.out
    }
    ```

Com podeu veure, també hem actualitzat les instruccions `view()` abans/després.
Tècnicament podríem haver utilitzat el mateix nom de variable (`greeting`) però l'hem actualitzat a alguna cosa més apropiada (`csv`) per fer el codi més llegible per altres.

#### 4.2.2. Executeu el workflow de nou

Provem d'executar el workflow amb la lògica d'anàlisi CSV afegida.

```bash
nextflow run hello-channels.nf
```

??? failure "Sortida de la comanda"

    ```console hl_lines="7-11 14 19"
     N E X T F L O W   ~  version 25.10.2

    Launching `hello-channels.nf` [insane_fermat] DSL2 - revision: 8e62fcbeb1

    executor >  local (3)
    [24/76da2f] sayHello (2) [  0%] 0 of 3 ✘
    Before splitCsv: /workspaces/training/hello-nextflow/data/greetings.csv
    After splitCsv: [Hello, English, 123]
    After splitCsv: [Bonjour, French, 456]
    After splitCsv: [Holà, Spanish, 789]
    ERROR ~ Error executing process > 'sayHello (2)'

    Caused by:
      Missing output file(s) `[Bonjour, French, 456]-output.txt` expected by process `sayHello (2)`


    Command executed:

      echo '[Bonjour, French, 456]' > '[Bonjour, French, 456]-output.txt'

    Command exit status:
      0

    Command output:
      (empty)

    Work dir:
      /workspaces/training/hello-nextflow/work/24/76da2fcc4876b61632749f99e26a50

    Tip: you can try to figure out what's wrong by changing to the process work dir and showing the script file named `.command.sh`

    -- Check '.nextflow.log' file for details
    ```

Curiosament, això també falla, però amb un error diferent.
Aquesta vegada Nextflow ha analitzat el contingut del fitxer (visca!) però ha carregat cada fila com un array, i cada array és un element al canal.

<figure class="excalidraw">
--8<-- "docs/en/docs/hello_nextflow/img/hello-channels-split-fail.svg"
</figure>

Necessitem dir-li que només prengui la primera columna de cada fila.
Així que com desempaquetem això?

Hem utilitzat prèviament `flatten()` per desempaquetar el contingut d'un canal, però això no funcionaria aquí perquè flatten desempaqueta _tot_ (sentiu-vos lliures de provar-ho si voleu veure-ho per vosaltres mateixos).

En lloc d'això, utilitzarem un altre operador anomenat `map()` que és realment útil i apareix molt en pipelines de Nextflow.

### 4.3. Utilitzeu l'operador `map()` per extreure les salutacions

L'operador [`map()`](https://nextflow.io/docs/latest/reference/operator.html#map) és una petita eina molt pràctica que ens permet fer tot tipus de mapatges al contingut d'un canal.

En aquest cas, l'utilitzarem per extreure aquell element que volem de cada fila al nostre fitxer de dades.
Aquesta és la sintaxi:

```groovy title="Syntax"
.map { row -> row[0] }
```

Això significa 'per a cada fila al canal, pren el 0è (primer) element que conté'.

Així que apliquem això al nostre anàlisi CSV.

#### 4.3.1. Apliqueu `map()` al canal

Al bloc workflow, feu el següent canvi de codi:

=== "Després"

    ```groovy title="hello-channels.nf" linenums="27" hl_lines="9 10"
    workflow {

        main:
        // crea un canal per a entrades des d'un fitxer CSV
        greeting_ch = channel.fromPath(params.input)
                             .view { csv -> "Before splitCsv: $csv" }
                             .splitCsv()
                             .view { csv -> "After splitCsv: $csv" }
                             .map { item -> item[0] }
                             .view { csv -> "After map: $csv" }
        // emet una salutacio
        sayHello(greeting_ch)

        publish:
        first_output = sayHello.out
    }
    ```

=== "Abans"

    ```groovy title="hello-channels.nf" linenums="27"
    workflow {

        main:
        // crea un canal per a entrades des d'un fitxer CSV
        greeting_ch = channel.fromPath(params.input)
                             .view { csv -> "Before splitCsv: $csv" }
                             .splitCsv()
                             .view { csv -> "After splitCsv: $csv" }
        // emet una salutacio
        sayHello(greeting_ch)

        publish:
        first_output = sayHello.out
    }
    ```

Veieu que hem afegit una altra crida `view()` per confirmar que l'operador fa el que esperem.

#### 4.3.2. Executeu el workflow

Executem això una vegada més:

```bash
nextflow run hello-channels.nf
```

??? success "Sortida de la comanda"

    ```console
     N E X T F L O W   ~  version 25.10.2

    Launching `hello-channels.nf` [focused_volhard] DSL2 - revision: de435e45be

    executor >  local (3)
    [54/6eebe3] sayHello (3) [100%] 3 of 3 ✔
    Before splitCsv: /workspaces/training/hello-nextflow/data/greetings.csv
    After splitCsv: [Hello, English, 123]
    After splitCsv: [Bonjour, French, 456]
    After splitCsv: [Holà, Spanish, 789]
    After map: Hello
    After map: Bonjour
    After map: Holà
    ```

Aquesta vegada hauria d'executar-se sense error.

Mirant la sortida de les instruccions `view()`, veieu el següent:

- Una única instrucció `Abans de splitCsv:`: en aquell moment el canal conté un element, el camí del fitxer original.
- Tres instruccions separades `Després de splitCsv:`: una per a cada salutació, però cadascuna està continguda dins d'un array que correspon a aquella línia del fitxer.
- Tres instruccions separades `Després de map:`: una per a cada salutació, que ara són elements individuals al canal.

_Tingueu en compte que les línies poden aparèixer en un ordre diferent a la vostra sortida._

<figure class="excalidraw">
--8<-- "docs/en/docs/hello_nextflow/img/hello-channels-split-and-map.svg"
</figure>

També podeu mirar els fitxers de sortida per verificar que cada salutació va ser correctament extreta i processada a través del workflow.

Hem aconseguit el mateix resultat que abans, però ara tenim molta més flexibilitat per afegir més elements al canal de salutacions que volem processar modificant un fitxer d'entrada, sense modificar cap codi.
Aprendreu enfocaments més sofisticats per gestionar entrades complexes en una formació posterior.

### Conclusió

Sabeu com utilitzar el constructor de canals `.fromPath()` i els operadors `splitCsv()` i `map()` per llegir un fitxer de valors d'entrada i gestionar-los adequadament.

Més generalment, teniu una comprensió bàsica de com Nextflow utilitza **canals** per gestionar entrades a processos i **operadors** per transformar el seu contingut.
També heu vist com els canals gestionen l'execució paral·lela implícitament.

<figure class="excalidraw">
--8<-- "docs/en/docs/hello_nextflow/img/hello-channels-parallel.svg"
</figure>

### Què segueix?

Feu una gran pausa, heu treballat dur en aquesta!

Quan estigueu preparats, passeu a [**Part 3: Hello Workflow**](./03_hello_workflow.md) per aprendre com afegir més passos i connectar-los junts en un workflow adequat.

---

## Qüestionari

<quiz>
Què és un canal a Nextflow?
- [ ] Una especificació de camí de fitxer
- [ ] Una definició de procés
- [x] Una estructura tipus cua per passar dades entre processos
- [ ] Una configuració

Més informació: [1.1. Creeu un canal d'entrada](#11-create-an-input-channel)
</quiz>

<quiz>
Què mostrarà aquest codi?

```groovy
channel.of('Hello', 'Bonjour', 'Hola')
    .view()
```

- [ ] `['Hello', 'Bonjour', 'Hola']` (una única llista)
- [x] Cada element en una línia separada: `Hello`, `Bonjour`, `Hola`
- [ ] Res (els canals no imprimeixen per defecte)
- [ ] Un error (sintaxi no vàlida)

Més informació: [1.1. Creeu un canal d'entrada](#11-create-an-input-channel)
</quiz>

<quiz>
Quan un canal conté múltiples valors, com gestiona Nextflow l'execució del procés?
- [ ] El procés s'executa una vegada amb tots els valors
- [x] El procés s'executa una vegada per a cada valor al canal
- [ ] El procés s'executa només amb el primer valor
- [ ] El procés s'executa només amb l'últim valor

Més informació: [2. Modifiqueu el workflow per executar-se amb múltiples valors d'entrada](#2-modify-the-workflow-to-run-on-multiple-input-values)
</quiz>

<quiz>
Què fa l'operador `flatten()`?
- [ ] Combina múltiples canals en un
- [ ] Ordena els elements del canal
- [x] Desempaqueta arrays en elements individuals
- [ ] Elimina elements duplicats

Més informació: [3.2.1. Afegiu l'operador `flatten()`](#321-add-the-flatten-operator)
</quiz>

<quiz>
Quin és el propòsit de l'operador `view()`?
- [ ] Filtrar el contingut del canal
- [ ] Transformar els elements del canal
- [x] Inspeccionar i depurar el contingut del canal
- [ ] Desar el contingut del canal a un fitxer

Més informació: [1.4. Utilitzeu `view()` per inspeccionar el contingut del canal](#14-use-view-to-inspect-the-channel-contents)
</quiz>

<quiz>
Què fa `splitCsv()`?
- [ ] Crea un fitxer CSV a partir del contingut del canal
- [ ] Divideix una cadena per comes
- [x] Analitza un fitxer CSV en arrays que representen cada fila
- [ ] Fusiona múltiples fitxers CSV

Més informació: [4.2. Utilitzeu l'operador `splitCsv()` per analitzar el fitxer](#42-use-the-splitcsv-operator-to-parse-the-file)
</quiz>

<quiz>
Quin és el propòsit de l'operador `map()`?
- [ ] Filtrar elements d'un canal
- [ ] Combinar múltiples canals
- [x] Transformar cada element en un canal
- [ ] Comptar elements en un canal

Més informació: [4.3. Utilitzeu l'operador `map()` per extreure les salutacions](#43-use-the-map-operator-to-extract-the-greetings)
</quiz>

<quiz>
Per què és important utilitzar noms de fitxer de sortida dinàmics quan es processen múltiples entrades?
- [ ] Per millorar el rendiment
- [ ] Per reduir l'espai de disc
- [x] Per evitar que els fitxers de sortida se sobreescriguin entre ells
- [ ] Per habilitar la funcionalitat de resume

Més informació: [2.2. Assegureu-vos que els noms dels fitxers de sortida seran únics](#22-ensure-the-output-file-names-will-be-unique)
</quiz>
