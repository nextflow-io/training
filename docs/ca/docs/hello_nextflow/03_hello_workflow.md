# Part 3: Hello Workflow

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Traducció assistida per IA - [més informació i suggeriments](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

<div class="video-wrapper">
  <iframe width="560" height="315" src="https://www.youtube.com/embed/_aO56V3iXGI?si=y8lAedhEHWaTV4zd&amp;list=PLPZ8WHdZGxmWKozQuzr27jyMGqp9kElVK&amp;cc_load_policy=1" title="YouTube video player" frameborder="0" allow="accelerometer; autoplay; clipboard-write; encrypted-media; gyroscope; picture-in-picture; web-share" referrerpolicy="strict-origin-when-cross-origin" allowfullscreen></iframe>
</div>

/// caption
:fontawesome-brands-youtube:{ .youtube } Consulteu [la llista de reproducció completa](https://youtube.com/playlist?list=PLPZ8WHdZGxmWKozQuzr27jyMGqp9kElVK&si=eF7cLR62goy-lc6n) al canal de YouTube de Nextflow.

:green_book: La transcripció del vídeo està disponible [aquí](./transcripts/03_hello_workflow.md).
///

La majoria de workflows del món real impliquen més d'un pas.
En aquest mòdul de formació, aprendreu com connectar processos entre si en un workflow de múltiples passos.

Això us ensenyarà la manera Nextflow d'aconseguir el següent:

1. Fer que les dades flueixin d'un procés al següent
2. Recollir sortides de múltiples crides de procés en una única crida de procés
3. Passar paràmetres addicionals a un procés
4. Gestionar múltiples sortides que surten d'un procés

Per demostrar-ho, continuarem construint sobre l'exemple Hello World independent del domini de les Parts 1 i 2.
Aquesta vegada, farem els següents canvis al nostre workflow per reflectir millor com la gent construeix workflows reals:

1. Afegir un segon pas que converteixi la salutació a majúscules.
2. Afegir un tercer pas que reculli totes les salutacions transformades i les escrigui en un únic fitxer.
3. Afegir un paràmetre per nomenar el fitxer de sortida final i passar-lo com a entrada secundària al pas de recollida.
4. Fer que el pas de recollida també informi d'una estadística simple sobre el que s'ha processat.

??? info "Com començar des d'aquesta secció"

    Aquesta secció del curs assumeix que heu completat les Parts 1-2 del curs [Hello Nextflow](./index.md), però si us sentiu còmodes amb els conceptes bàsics coberts en aquestes seccions, podeu començar des d'aquí sense fer res especial.

---

## 0. Escalfament: Executeu `hello-workflow.nf`

Utilitzarem l'script de workflow `hello-workflow.nf` com a punt de partida.
És equivalent a l'script produït en treballar la Part 2 d'aquest curs de formació, excepte que hem eliminat les instruccions `view()` i hem canviat la destinació de sortida:

```groovy title="hello-workflow.nf" linenums="37" hl_lines="3"
output {
    first_output {
        path 'hello_workflow'
        mode 'copy'
    }
}
```

Aquest diagrama resumeix l'operació actual del workflow.
Hauria de semblar familiar, excepte que ara estem mostrant explícitament que les sortides del procés s'empaqueten en un canal, igual que les entrades.
Utilitzarem aquest canal de sortida en un moment.

<figure class="excalidraw">
--8<-- "docs/en/docs/hello_nextflow/img/hello-workflow-channels.svg"
</figure>

Només per assegurar-nos que tot funciona, executeu l'script una vegada abans de fer cap canvi:

```bash
nextflow run hello-workflow.nf
```

??? success "Sortida de la comanda"

    ```console
     N E X T F L O W   ~  version 25.10.2

    Launching `hello-workflow.nf` [admiring_lamarr] DSL2 - revision: 4d4053520d

    executor >  local (3)
    [b1/5826b5] process > sayHello (2) [100%] 3 of 3 ✔
    ```

Com abans, trobareu els fitxers de sortida a la ubicació especificada al bloc `output`.
Per a aquest capítol, és sota `results/hello_workflow/`.

??? abstract "Contingut del directori"

    ```console
    results/hello_workflow
    ├── Bonjour-output.txt
    ├── Hello-output.txt
    └── Holà-output.txt
    ```

Si això us ha funcionat, esteu preparats per aprendre com muntar un workflow de múltiples passos.

---

## 1. Afegiu un segon pas al workflow

Afegirem un pas per convertir cada salutació a majúscules.

<figure class="excalidraw">
--8<-- "docs/en/docs/hello_nextflow/img/hello-multistep.svg"
</figure>

Per fer-ho, necessitem fer tres coses:

- Definir la comanda que utilitzarem per fer la conversió a majúscules.
- Escriure un nou procés que encapsuli la comanda de conversió a majúscules.
- Cridar el nou procés al bloc workflow i configurar-lo per prendre la sortida del procés `sayHello()` com a entrada.

### 1.1. Definiu la comanda de conversió a majúscules i proveu-la al terminal

Per fer la conversió de les salutacions a majúscules, utilitzarem una eina UNIX clàssica anomenada `tr` per 'text replacement', amb la següent sintaxi:

```bash title="Syntax"
tr '[a-z]' '[A-Z]'
```

Aquesta és una substitució de text molt naïf que no té en compte les lletres accentuades, així que per exemple 'Holà' es convertirà en 'HOLà', però farà una feina prou bona per demostrar els conceptes de Nextflow i això és el que importa.

Per provar-ho, podem executar la comanda `echo 'Hello World'` i canalitzar la seva sortida a la comanda `tr`:

```bash
echo 'Hello World' | tr '[a-z]' '[A-Z]' > UPPER-output.txt
```

La sortida és un fitxer de text anomenat `UPPER-output.txt` que conté la versió en majúscules de la cadena `Hello World`.

??? abstract "Contingut del fitxer"

    ```console title="UPPER-output.txt"
    HELLO WORLD
    ```

Això és bàsicament el que intentarem fer amb el nostre workflow.

### 1.2. Escriviu el pas de conversió a majúscules com a procés Nextflow

Podem modelar el nostre nou procés sobre el primer, ja que volem utilitzar tots els mateixos components.

Afegiu la següent definició de procés a l'script de workflow, just sota el primer:

```groovy title="hello-workflow.nf" linenums="20"
/*
 * Utilitza una eina de substitució de text per convertir la salutació a majúscules
 */
process convertToUpper {

    input:
    path input_file

    output:
    path "UPPER-${input_file}"

    script:
    """
    cat '${input_file}' | tr '[a-z]' '[A-Z]' > 'UPPER-${input_file}'
    """
}
```

En aquest, compondrem el segon nom de fitxer de sortida basat en el nom del fitxer d'entrada, de manera similar al que vam fer originalment per a la sortida del primer procés.

### 1.3. Afegiu una crida al nou procés al bloc workflow

Ara necessitem dir a Nextflow que realment cridi el procés que acabem de definir.

Al bloc workflow, feu el següent canvi de codi:

=== "Després"

    ```groovy title="hello-workflow.nf" linenums="44" hl_lines="10-11"
    workflow {

        main:
        // crea un canal per a entrades des d'un fitxer CSV
        greeting_ch = channel.fromPath(params.input)
                            .splitCsv()
                            .map { line -> line[0] }
        // emet una salutació
        sayHello(greeting_ch)
        // converteix la salutació a majúscules
        convertToUpper()

        publish:
        first_output = sayHello.out
    }
    ```

=== "Abans"

    ```groovy title="hello-workflow.nf" linenums="44"
    workflow {

        main:
        // crea un canal per a entrades des d'un fitxer CSV
        greeting_ch = channel.fromPath(params.input)
                            .splitCsv()
                            .map { line -> line[0] }
        // emet una salutació
        sayHello(greeting_ch)

        publish:
        first_output = sayHello.out
    }
    ```

Això encara no és funcional perquè no hem especificat què hauria de ser l'entrada del procés `convertToUpper()`.

### 1.4. Passeu la sortida del primer procés al segon procés

Ara necessitem fer que la sortida del procés `sayHello()` flueixi cap al procés `convertToUpper()`.

Convenientment, Nextflow empaqueta automàticament la sortida d'un procés en un canal, com es mostra al diagrama de la secció d'escalfament.
Podem referir-nos al canal de sortida d'un procés com `<process>.out`.

Així doncs, la sortida del procés `sayHello` és un canal anomenat `sayHello.out`, que podem connectar directament a la crida a `convertToUpper()`.

<figure class="excalidraw">
--8<-- "docs/en/docs/hello_nextflow/img/hello-multistep-connector.svg"
</figure>

Al bloc workflow, feu el següent canvi de codi:

=== "Després"

    ```groovy title="hello-workflow.nf" linenums="53" hl_lines="2"
        // converteix la salutació a majúscules
        convertToUpper(sayHello.out)
    ```

=== "Abans"

    ```groovy title="hello-workflow.nf" linenums="53" hl_lines="2"
        // converteix la salutació a majúscules
        convertToUpper()
    ```

Per a un cas simple com aquest (una sortida a una entrada), això és tot el que necessitem fer per connectar dos processos!

### 1.5. Configureu la publicació de sortida del workflow

Finalment, actualitzem les sortides del workflow per publicar també els resultats del segon procés.

#### 1.5.1. Actualitzeu la secció `publish:` del bloc `workflow`

Al bloc `workflow`, feu el següent canvi de codi:

=== "Després"

    ```groovy title="hello-workflow.nf" linenums="56" hl_lines="3"
        publish:
        first_output = sayHello.out
        uppercased = convertToUpper.out
    }
    ```

=== "Abans"

    ```groovy title="hello-workflow.nf" linenums="56"
        publish:
        first_output = sayHello.out
    }
    ```

La lògica és la mateixa que abans.

#### 1.5.2. Actualitzeu el bloc `output`

Al bloc `output`, feu el següent canvi de codi:

=== "Després"

    ```groovy title="hello-workflow.nf" linenums="61" hl_lines="6-9"
    output {
        first_output {
            path 'hello_workflow'
            mode 'copy'
        }
        uppercased {
            path 'hello_workflow'
            mode 'copy'
        }
    }
    ```

=== "Abans"

    ```groovy title="hello-workflow.nf" linenums="61"
    output {
        first_output {
            path 'hello_workflow'
            mode 'copy'
        }
    }
    ```

Una vegada més, la lògica és la mateixa que abans.

Això us mostra que podeu controlar la configuració de sortida a un nivell molt granular, per a cada sortida individual.
Sentiu-vos lliures de provar canviar els camins o el mode de publicació per a un dels processos per veure què passa.

Per descomptat, això significa que estem repetint alguna informació aquí, cosa que podria resultar inconvenient si volguéssim actualitzar la ubicació per a totes les sortides de la mateixa manera.
Més endavant al curs, aprendreu com configurar aquests paràmetres per a múltiples sortides de manera estructurada.

### 1.6. Executeu el workflow amb `-resume`

Provem això utilitzant la bandera `-resume`, ja que ja hem executat el primer pas del workflow amb èxit.

```bash
nextflow run hello-workflow.nf -resume
```

??? success "Sortida de la comanda"

    ```console
     N E X T F L O W   ~  version 25.10.2

    Launching `hello-workflow.nf` [high_cantor] DSL2 - revision: d746983511

    executor >  local (3)
    [ab/816321] process > sayHello (3)       [100%] 3 of 3, cached: 3 ✔
    [e0/ecf81b] process > convertToUpper (3) [100%] 3 of 3 ✔
    ```

Ara hi ha una línia extra a la sortida de la consola que correspon al nou procés que acabem d'afegir.

Trobareu les sortides al directori `results/hello_workflow` tal com s'ha establert al bloc `output`.

??? abstract "Contingut del directori"

    ```console
    results/hello_workflow/
    ├── Bonjour-output.txt
    ├── Hello-output.txt
    ├── Holà-output.txt
    ├── UPPER-Bonjour-output.txt
    ├── UPPER-Hello-output.txt
    └── UPPER-Holà-output.txt
    ```

Això és convenient! Però encara val la pena donar una ullada dins del directori de treball d'una de les crides al segon procés.

??? abstract "Contingut del directori"

    ```console
    work/e0/ecf81b4cacc648b9b994218d5b29d7/
    ├── Holà-output.txt -> /workspaces/training/hello-nextflow/work/ab/81632178cd37e9e815959278808819/Holà-output.txt
    └── UPPER-Holà-output.txt
    ```

Observeu que hi ha dos fitxers `*-output`: la sortida del primer procés així com la sortida del segon.

La sortida del primer procés hi és perquè Nextflow la va **preparar** allà per tenir tot el necessari per a l'execució dins del mateix subdirectori.

No obstant això, en realitat és un enllaç simbòlic que apunta al fitxer original al subdirectori de la primera crida de procés.
Per defecte, quan s'executa en una sola màquina com estem fent aquí, Nextflow utilitza enllaços simbòlics en lloc de còpies per preparar fitxers d'entrada i intermedis.

Ara, abans de continuar, penseu en com tot el que vam fer va ser connectar la sortida de `sayHello` a l'entrada de `convertToUpper` i els dos processos es van poder executar en sèrie.
Nextflow va fer la feina dura de gestionar fitxers d'entrada i sortida individuals i passar-los entre les dues comandes per nosaltres.

Aquesta és una de les raons per les quals els canals de Nextflow són tan potents: s'encarreguen de la feina rutinària implicada en connectar passos de workflow junts.

### Conclusió

Sabeu com encadenar processos junts proporcionant la sortida d'un pas com a entrada al següent pas.

### Què segueix?

Apreneu com recollir sortides de crides de procés per lots i alimentar-les en un únic procés.

---

## 2. Afegiu un tercer pas per recollir totes les salutacions

Quan utilitzem un procés per aplicar una transformació a cadascun dels elements d'un canal, com estem fent aquí amb les múltiples salutacions, de vegades volem recollir elements del canal de sortida d'aquest procés i alimentar-los en un altre procés que realitza algun tipus d'anàlisi o sumatori.

Per demostrar-ho, afegirem un nou pas al nostre pipeline que recull totes les salutacions en majúscules produïdes pel procés `convertToUpper` i les escriu en un únic fitxer.

<figure class="excalidraw">
--8<-- "docs/en/docs/hello_nextflow/img/hello-collect.svg"
</figure>

Per no espatllar la sorpresa, però això implicarà un operador molt útil.

### 2.1. Definiu la comanda de recollida i proveu-la al terminal

El pas de recollida que volem afegir al nostre workflow utilitzarà la comanda `cat` per concatenar múltiples salutacions en majúscules en un únic fitxer.

Executem la comanda per si sola al terminal per verificar que funciona com s'espera, igual que hem fet abans.

Executeu el següent al vostre terminal:

```bash
echo 'Hello' | tr '[a-z]' '[A-Z]' > UPPER-Hello-output.txt
echo 'Bonjour' | tr '[a-z]' '[A-Z]' > UPPER-Bonjour-output.txt
echo 'Holà' | tr '[a-z]' '[A-Z]' > UPPER-Holà-output.txt
cat UPPER-Hello-output.txt UPPER-Bonjour-output.txt UPPER-Holà-output.txt > COLLECTED-output.txt
```

La sortida és un fitxer de text anomenat `COLLECTED-output.txt` que conté les versions en majúscules de les salutacions originals.

??? abstract "Contingut del fitxer"

    ```console title="COLLECTED-output.txt"
    HELLO
    BONJOUR
    HOLà
    ```

Aquest és el resultat que volem aconseguir amb el nostre workflow.

### 2.2. Creeu un nou procés per fer el pas de recollida

Creem un nou procés i l'anomenem `collectGreetings()`.
Podem començar a escriure'l basat en el que hem vist abans.

#### 2.2.1. Escriviu les parts 'òbvies' del procés

Afegiu la següent definició de procés a l'script de workflow:

```groovy title="hello-workflow.nf" linenums="37"
/*
 * Recull salutacions en majúscules en un únic fitxer de sortida
 */
process collectGreetings {

    input:
    ???

    output:
    path "COLLECTED-output.txt"

    script:
    """
    cat ??? > 'COLLECTED-output.txt'
    """
}
```

Això és el que podem escriure amb confiança basat en el que heu après fins ara.
Però això no és funcional!
Deixa de banda la(es) definició(ns) d'entrada i la primera meitat de la comanda script perquè necessitem esbrinar com escriure això.

#### 2.2.2. Definiu les entrades a `collectGreetings()`

Necessitem recollir les salutacions de totes les crides al procés `convertToUpper()`.
Què sabem que podem obtenir del pas anterior al workflow?

El canal de sortida de `convertToUpper()` contindrà els camins als fitxers individuals que contenen les salutacions en majúscules.
Això equival a una ranura d'entrada; anomenem-la `input_files` per simplicitat.

Al bloc de procés, feu el següent canvi de codi:

=== "Després"

    ```groovy title="hello-workflow.nf" linenums="42" hl_lines="2"
          input:
          path input_files
    ```

=== "Abans"

    ```groovy title="hello-workflow.nf" linenums="42" hl_lines="2"
          input:
          ???
    ```

Observeu que utilitzem el prefix `path` encara que esperem que això contingui múltiples fitxers.

#### 2.2.3. Compondreu la comanda de concatenació

Aquí és on les coses podrien complicar-se una mica, perquè necessitem poder gestionar un nombre arbitrari de fitxers d'entrada.
Específicament, no podem escriure la comanda per endavant, així que necessitem dir a Nextflow com compondre-la en temps d'execució basat en quines entrades flueixen cap al procés.

En altres paraules, si tenim un canal d'entrada que conté l'element `[file1.txt, file2.txt, file3.txt]`, necessitem que Nextflow ho converteixi en `cat file1.txt file2.txt file3.txt`.

Afortunadament, Nextflow està molt content de fer-ho per nosaltres si simplement escrivim `cat ${input_files}` a la comanda script.

Al bloc de procés, feu el següent canvi de codi:

=== "Després"

    ```groovy title="hello-workflow.nf" linenums="54" hl_lines="3"
        script:
        """
        cat ${input_files} > 'COLLECTED-output.txt'
        """
    ```

=== "Abans"

    ```groovy title="hello-workflow.nf" linenums="54"
        script:
        """
        cat ??? > 'COLLECTED-output.txt'
        """
    ```

En teoria això hauria de gestionar qualsevol nombre arbitrari de fitxers d'entrada.

!!! tip "Consell"

    Algunes eines de línia de comandes requereixen proporcionar un argument (com `-input`) per a cada fitxer d'entrada.
    En aquest cas, hauríem de fer una mica de feina extra per compondre la comanda.
    Podeu veure un exemple d'això al curs de formació [Nextflow for Genomics](../../nf4_science/genomics/).

### 2.3. Afegiu el pas de recollida al workflow

Ara només hauríem de necessitar cridar el procés de recollida sobre la sortida del pas de conversió a majúscules.
Això també és un canal, anomenat `convertToUpper.out`.

<figure class="excalidraw">
--8<-- "docs/en/docs/hello_nextflow/img/hello-collect-connector.svg"
</figure>

#### 2.3.1. Connecteu les crides de procés

Al bloc workflow, feu el següent canvi de codi:

=== "Després"

    ```groovy title="hello-workflow.nf" linenums="75" hl_lines="4 5"
        // converteix la salutació a majúscules
        convertToUpper(sayHello.out)

        // recull totes les salutacions en un fitxer
        collectGreetings(convertToUpper.out)
    }
    ```

=== "Abans"

    ```groovy title="hello-workflow.nf" linenums="75"
        // converteix la salutació a majúscules
        convertToUpper(sayHello.out)
    }
    ```

Això connecta la sortida de `convertToUpper()` a l'entrada de `collectGreetings()`.

#### 2.3.2. Executeu el workflow amb `-resume`

Provem-ho.

```bash
nextflow run hello-workflow.nf -resume
```

??? success "Sortida de la comanda"

    ```console hl_lines="8"
    N E X T F L O W   ~  version 25.10.2

    Launching `hello-workflow.nf` [mad_gilbert] DSL2 - revision: 6acfd5e28d

    executor >  local (3)
    [79/33b2f0] sayHello (2)         | 3 of 3, cached: 3 ✔
    [99/79394f] convertToUpper (3)   | 3 of 3, cached: 3 ✔
    [47/50fe4a] collectGreetings (1) | 3 of 3 ✔
    ```

S'executa amb èxit, incloent el tercer pas.

No obstant això, mireu el nombre de crides per a `collectGreetings()` a l'última línia.
Només esperàvem una, però n'hi ha tres.

Ara doneu una ullada al contingut del fitxer de sortida final.

??? abstract "Contingut del fitxer"

    ```console title="results/COLLECTED-output.txt"
    Holà
    ```

Oh no. El pas de recollida es va executar individualment en cada salutació, cosa que NO és el que volíem.

<figure class="excalidraw">
--8<-- "docs/en/docs/hello_nextflow/img/hello-collect-no-operator.svg"
</figure>

Necessitem fer alguna cosa per dir a Nextflow explícitament que volem que aquest tercer pas s'executi sobre tots els elements del canal de sortida de `convertToUpper()`.

### 2.4. Utilitzeu un operador per recollir les salutacions en una única entrada

Sí, una vegada més la resposta al nostre problema és un operador.

Específicament, utilitzarem l'operador amb el nom apropiat [`collect()`](https://nextflow.io/docs/latest/reference/operator.html#collect).

#### 2.4.1. Afegiu l'operador `collect()`

Aquesta vegada semblarà una mica diferent perquè no estem afegint l'operador en el context d'una factoria de canals; l'estem afegint a un canal de sortida.

Prenem el `convertToUpper.out` i afegim l'operador `collect()`, cosa que ens dóna `convertToUpper.out.collect()`.
Podem connectar això directament a la crida del procés `collectGreetings()`.

Al bloc workflow, feu el següent canvi de codi:

=== "Després"

    ```groovy title="hello-workflow.nf" linenums="73" hl_lines="2"
        // recull totes les salutacions en un fitxer
        collectGreetings(convertToUpper.out.collect())
    }
    ```

=== "Abans"

    ```groovy title="hello-workflow.nf" linenums="73" hl_lines="2"
        // recull totes les salutacions en un fitxer
        collectGreetings(convertToUpper.out)
    }
    ```

#### 2.4.2. Afegiu algunes instruccions `view()`

També inclourem un parell d'instruccions `view()` per visualitzar els estats abans i després del contingut del canal.

=== "Després"

    ```groovy title="hello-workflow.nf" linenums="73" hl_lines="4-6"
        // recull totes les salutacions en un fitxer
        collectGreetings(convertToUpper.out.collect())

        // instruccions view opcionals
        convertToUpper.out.view { contents -> "Before collect: $contents" }
        convertToUpper.out.collect().view { contents -> "After collect: $contents" }
    }
    ```

=== "Abans"

    ```groovy title="hello-workflow.nf" linenums="73"
        // recull totes les salutacions en un fitxer
        collectGreetings(convertToUpper.out.collect())
    }
    ```

Les instruccions `view()` poden anar on vulgueu; les hem posat just després de la crida per llegibilitat.

#### 2.4.3. Executeu el workflow de nou amb `-resume`

Provem-ho:

```bash
nextflow run hello-workflow.nf -resume
```

??? success "Sortida de la comanda"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `hello-workflow.nf` [soggy_franklin] DSL2 - revision: bc8e1b2726

    [d6/cdf466] sayHello (1)       | 3 of 3, cached: 3 ✔
    [99/79394f] convertToUpper (2) | 3 of 3, cached: 3 ✔
    [1e/83586c] collectGreetings   | 1 of 1 ✔
    Before collect: /workspaces/training/hello-nextflow/work/b3/d52708edba8b864024589285cb3445/UPPER-Bonjour-output.txt
    Before collect: /workspaces/training/hello-nextflow/work/99/79394f549e3040dfc2440f69ede1fc/UPPER-Hello-output.txt
    Before collect: /workspaces/training/hello-nextflow/work/aa/56bfe7cf00239dc5badc1d04b60ac4/UPPER-Holà-output.txt
    After collect: [/workspaces/training/hello-nextflow/work/b3/d52708edba8b864024589285cb3445/UPPER-Bonjour-output.txt, /workspaces/training/hello-nextflow/work/99/79394f549e3040dfc2440f69ede1fc/UPPER-Hello-output.txt, /workspaces/training/hello-nextflow/work/aa/56bfe7cf00239dc5badc1d04b60ac4/UPPER-Holà-output.txt]
    ```

S'executa amb èxit, encara que la sortida del registre pot semblar una mica més desordenada que això (l'hem netejat per llegibilitat).

Aquesta vegada el tercer pas només es va cridar una vegada!
Mirant la sortida de les instruccions `view()`, veiem el següent:

- Tres instruccions `Abans de collect:`, una per a cada salutació: en aquest punt els camins de fitxer són elements individuals al canal.
- Una única instrucció `Després de collect:`: els tres camins de fitxer ara estan empaquetats en un únic element.

Podem resumir això amb el següent diagrama:

<figure class="excalidraw">
--8<-- "docs/en/docs/hello_nextflow/img/hello-collect-WITH-operator.svg"
</figure>

Finalment, podeu donar una ullada al contingut del fitxer de sortida per satisfer-vos que tot ha funcionat correctament.

??? abstract "Contingut del fitxer"

    ```console title="results/COLLECTED-output.txt"
    BONJOUR
    HELLO
    HOLà
    ```

Aquesta vegada tenim les tres salutacions al fitxer de sortida final. Èxit!

!!! note "Nota"

    Si executeu això diverses vegades sense `-resume`, veureu que l'ordre de les salutacions canvia d'una execució a l'altra.
    Això us mostra que l'ordre en què els elements flueixen a través de les crides de procés no està garantit que sigui consistent.

#### 2.4.4. Elimineu les instruccions `view()` per llegibilitat

Abans de passar a la següent secció, us recomanem que elimineu les instruccions `view()` per evitar embolicar la sortida de la consola.

=== "Després"

    ```groovy title="hello-workflow.nf" linenums="73"
        // recull totes les salutacions en un fitxer
        collectGreetings(convertToUpper.out.collect())
    ```

=== "Abans"

    ```groovy title="hello-workflow.nf" linenums="73" hl_lines="4-6"
        // recull totes les salutacions en un fitxer
        collectGreetings(convertToUpper.out.collect())

        // instruccions view opcionals
        convertToUpper.out.view { contents -> "Before collect: $contents" }
        convertToUpper.out.collect().view { contents -> "After collect: $contents" }
    ```

Això és bàsicament l'operació inversa del punt 2.4.2.

### Conclusió

Sabeu com recollir sortides d'un lot de crides de procés i alimentar-les en un pas d'anàlisi o sumatori conjunt.

Per recapitular, això és el que heu construït fins ara:

<figure class="excalidraw">
--8<-- "docs/en/docs/hello_nextflow/img/hello-collect.svg"
</figure>

### Què segueix?

Apreneu com passar més d'una entrada a un procés.

---

## 3. Passeu paràmetres addicionals a un procés

Volem poder nomenar el fitxer de sortida final alguna cosa específica per poder processar lots posteriors de salutacions sense sobreescriure els resultats finals.

Per fer-ho, farem els següents refinaments al workflow:

- Modificar el procés col·lector per acceptar un nom definit per l'usuari per al fitxer de sortida (`batch_name`)
- Afegir un paràmetre de línia de comandes al workflow (`--batch`) i passar-lo al procés col·lector

<figure class="excalidraw">
--8<-- "docs/en/docs/hello_nextflow/img/hello-collect-batch.svg"
</figure>

### 3.1. Modifiqueu el procés col·lector

Necessitarem declarar l'entrada addicional i integrar-la al nom del fitxer de sortida.

#### 3.1.1. Declareu l'entrada addicional

Bones notícies: podem declarar tantes variables d'entrada com vulguem a la definició del procés.
Anomenem aquesta `batch_name`.

Al bloc de procés, feu el següent canvi de codi:

=== "Després"

    ```groovy title="hello-workflow.nf" linenums="42" hl_lines="3"
        input:
        path input_files
        val batch_name
    ```

=== "Abans"

    ```groovy title="hello-workflow.nf" linenums="42"
        input:
        path input_files
    ```

Podeu configurar els vostres processos per esperar tantes entrades com vulgueu.
Ara mateix, totes estan configurades per ser entrades requerides; _heu_ de proporcionar un valor perquè el workflow funcioni.

Aprendreu com gestionar entrades requerides vs. opcionals més endavant al vostre viatge amb Nextflow.

#### 3.1.2. Utilitzeu la variable `batch_name` al nom del fitxer de sortida

Podem inserir la variable al nom del fitxer de sortida de la mateixa manera que hem compost noms de fitxer dinàmics abans.

Al bloc de procés, feu el següent canvi de codi:

=== "Després"

    ```groovy title="hello-workflow.nf" linenums="46" hl_lines="2 6"
        output:
        path "COLLECTED-${batch_name}-output.txt"

        script:
        """
        cat ${input_files} > 'COLLECTED-${batch_name}-output.txt'
        """
    ```

=== "Abans"

    ```groovy title="hello-workflow.nf" linenums="46" hl_lines="2 6"
        output:
        path "COLLECTED-output.txt"

        script:
        """
        cat ${input_files} > 'COLLECTED-output.txt'
        """
    ```

Això configura el procés per utilitzar el valor `batch_name` per generar un nom de fitxer específic per a la sortida final del workflow.

### 3.2. Afegiu un paràmetre de línia de comandes `batch`

Ara necessitem una manera de subministrar el valor per a `batch_name` i alimentar-lo a la crida del procés.

#### 3.2.1. Utilitzeu `params` per configurar el paràmetre

Ja sabeu com utilitzar el sistema `params` per declarar paràmetres CLI.
Utilitzem això per declarar un paràmetre `batch` (amb un valor per defecte perquè som mandrosos).

A la secció de paràmetres del pipeline, feu els següents canvis de codi:

=== "Després"

    ```groovy title="hello-workflow.nf" linenums="55" hl_lines="6"
    /*
     * Paràmetres del pipeline
     */
    params {
        input: Path = 'data/greetings.csv'
        batch: String = 'batch'
    }
    ```

=== "Abans"

    ```groovy title="hello-workflow.nf" linenums="55"
    /*
     * Paràmetres del pipeline
     */
    params {
        input: Path = 'data/greetings.csv'
    }
    ```

Igual que vam demostrar per a `--input`, podeu sobreescriure aquest valor per defecte especificant un valor amb `--batch` a la línia de comandes.

#### 3.2.2. Passeu el paràmetre `batch` al procés

Per proporcionar el valor del paràmetre al procés, necessitem afegir-lo a la crida del procés.

Al bloc workflow, feu el següent canvi de codi:

=== "Després"

    ```groovy title="hello-workflow.nf" linenums="74" hl_lines="2"
        // recull totes les salutacions en un fitxer
        collectGreetings(convertToUpper.out.collect(), params.batch)
    ```

=== "Abans"

    ```groovy title="hello-workflow.nf" linenums="74" hl_lines="2"
        // recull totes les salutacions en un fitxer
        collectGreetings(convertToUpper.out.collect())
    ```

Veieu que per proporcionar múltiples entrades a un procés, simplement les llisteu als parèntesis de la crida, separades per comes.

!!! warning "Advertència"

    HEU de proporcionar les entrades al procés en el MATEIX ORDRE EXACTE en què estan llistades al bloc de definició d'entrada del procés.

### 3.3. Executeu el workflow

Provem d'executar això amb un nom de lot a la línia de comandes.

```bash
nextflow run hello-workflow.nf -resume --batch trio
```

??? success "Sortida de la comanda"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `hello-workflow.nf` [confident_rutherford] DSL2 - revision: bc58af409c

    executor >  local (1)
    [79/33b2f0] sayHello (2)       | 3 of 3, cached: 3 ✔
    [99/79394f] convertToUpper (2) | 3 of 3, cached: 3 ✔
    [b5/f19efe] collectGreetings   | 1 of 1 ✔
    ```

S'executa amb èxit i produeix la sortida desitjada:

??? abstract "Contingut del fitxer"

    ```console title="results/COLLECTED-trio-output.txt"
    HELLO
    BONJOUR
    HOLà
    ```

Ara, sempre que especifiquem el paràmetre adequadament, les execucions posteriors sobre altres lots d'entrades no destruiran els resultats anteriors.

### Conclusió

Sabeu com passar més d'una entrada a un procés.

### Què segueix?

Apreneu com emetre múltiples sortides i gestionar-les convenientment.

---

## 4. Afegiu una sortida al pas col·lector

Fins ara hem estat utilitzant processos que només produïen una sortida cadascun.
Vam poder accedir a les seves respectives sortides molt convenientment utilitzant la sintaxi `<process>.out`, que vam utilitzar tant en el context de passar una sortida al següent procés (p. ex. `convertToUpper(sayHello.out)`) com en el context de la secció `publish:` (p. ex. `first_output = sayHello.out`).

Què passa quan un procés produeix més d'una?
Com gestionem les múltiples sortides?
Podem seleccionar i utilitzar una sortida específica?

Totes són excel·lents preguntes, i la resposta curta és que sí que podem!

Les múltiples sortides s'empaquetaran en canals separats.
Podem optar per donar noms a aquests canals de sortida, cosa que facilita referir-nos-hi individualment més endavant, o podem referir-nos-hi per índex.

Per a fins de demostració, diguem que volem comptar el nombre de salutacions que s'estan recollint per a un lot donat d'entrades i informar-ho en un fitxer.

### 4.1. Modifiqueu el procés per comptar i emetre el nombre de salutacions

Això requerirà dos canvis clau a la definició del procés: necessitem una manera de comptar les salutacions i escriure un fitxer d'informe, després necessitem afegir aquest fitxer d'informe al bloc `output` del procés.

#### 4.1.1. Compteu el nombre de salutacions recollides

Convenientment, Nextflow ens permet afegir codi arbitrari al bloc `script:` de la definició del procés, cosa que és molt útil per fer coses com aquesta.

Això significa que podem utilitzar la funció integrada de Nextflow `size()` per obtenir el nombre de fitxers a l'array `input_files`, i escriure el resultat a un fitxer amb una comanda `echo`.

Al bloc de procés `collectGreetings`, feu els següents canvis de codi:

=== "Després"

    ```groovy title="hello-workflow.nf" linenums="55" hl_lines="2 5"
        script:
        count_greetings = input_files.size()
        """
        cat ${input_files} > 'COLLECTED-${batch_name}-output.txt'
        echo 'There were ${count_greetings} greetings in this batch.' > '${batch_name}-report.txt'
        """
    ```

=== "Abans"

    ```groovy title="hello-workflow.nf" linenums="55"
        script:
        """
        cat ${input_files} > 'COLLECTED-${batch_name}-output.txt'
        """
    ```

La variable `count_greetings` es calcularà en temps d'execució.

#### 4.1.2. Emeteu el fitxer d'informe i anomeneu les sortides

En principi tot el que necessitem fer és afegir el fitxer d'informe al bloc `output:`.

No obstant això, mentre hi som, també afegirem algunes etiquetes `emit:` a les nostres declaracions de sortida. Aquestes ens permetran seleccionar les sortides per nom en lloc d'haver d'utilitzar índexs posicionals.

Al bloc de procés, feu el següent canvi de codi:

=== "Després"

    ```groovy title="hello-workflow.nf" linenums="46" hl_lines="2 3"
        output:
        path "COLLECTED-${batch_name}-output.txt", emit: outfile
        path "${batch_name}-report.txt", emit: report
    ```

=== "Abans"

    ```groovy title="hello-workflow.nf" linenums="46"
        output:
        path "COLLECTED-${batch_name}-output.txt"
    ```

Les etiquetes `emit:` són opcionals, i podríem haver afegit una etiqueta només a una de les sortides.
Però com diu el refrany, per què no ambdues?

!!! tip "Consell"

    Si no anomeneu les sortides d'un procés utilitzant `emit:`, encara podeu accedir-hi individualment utilitzant el seu respectiu índex (basat en zero).
    Per exemple, utilitzaríeu `<process>.out[0]` per obtenir la primera sortida, `<process>.out[1]` per obtenir la segona sortida, i així successivament.

    Preferim anomenar les sortides perquè altrament és massa fàcil agafar l'índex equivocat per error, especialment quan el procés produeix moltes sortides.

### 4.2. Actualitzeu les sortides del workflow

Ara que tenim dues sortides sortint del procés `collectGreetings`, la sortida `collectGreetings.out` conté dos canals:

- `collectGreetings.out.outfile` conté el fitxer de sortida final
- `collectGreetings.out.report` conté el fitxer d'informe

<figure class="excalidraw">
--8<-- "docs/en/docs/hello_nextflow/img/hello-collect-report.svg"
</figure>

Necessitem actualitzar les sortides del workflow en conseqüència.

#### 4.2.1. Actualitzeu la secció `publish:`

Al `bloc workflow`, feu el següent canvi de codi:

=== "Després"

    ```groovy title="hello-workflow.nf" linenums="80" hl_lines="4 5"
        publish:
        first_output = sayHello.out
        uppercased = convertToUpper.out
        collected = collectGreetings.out.outfile
        batch_report = collectGreetings.out.report
    ```

=== "Abans"

    ```groovy title="hello-workflow.nf" linenums="80" hl_lines="4"
        publish:
        first_output = sayHello.out
        uppercased = convertToUpper.out
        collected = collectGreetings.out
    ```

Com podeu veure, referir-se a sortides de procés específiques ara és trivial.
Quan anem a afegir un pas més al nostre pipeline a la Part 5 (Contenidors), podrem referir-nos fàcilment a `collectGreetings.out.outfile` i passar-lo al nou procés (spoiler: el nou procés s'anomena `cowpy`).

Però per ara, acabem d'actualitzar les sortides a nivell de workflow.

#### 4.2.2. Actualitzeu el bloc `output`

Al bloc `output`, feu el següent canvi de codi:

=== "Després"

    ```groovy title="hello-workflow.nf" linenums="86" hl_lines="14-17"
    output {
        first_output {
            path 'hello_workflow'
            mode 'copy'
        }
        uppercased {
            path 'hello_workflow'
            mode 'copy'
        }
        collected {
            path 'hello_workflow'
            mode 'copy'
        }
        batch_report {
            path 'hello_workflow'
            mode 'copy'
        }
    }
    ```

=== "Abans"

    ```groovy title="hello-workflow.nf" linenums="80"
    output {
        first_output {
            path 'hello_workflow'
            mode 'copy'
        }
        uppercased {
            path 'hello_workflow'
            mode 'copy'
        }
        collected {
            path 'hello_workflow'
            mode 'copy'
        }
    }
    ```

No necessitem actualitzar la definició de sortida `collected` ja que aquest nom no ha canviat.
Només necessitem afegir la nova sortida.

### 4.3. Executeu el workflow

Provem d'executar això amb el lot actual de salutacions.

```bash
nextflow run hello-workflow.nf -resume --batch trio
```

??? success "Sortida de la comanda"

    ```console
     N E X T F L O W   ~  version 25.10.2

    Launching `hello-workflow.nf` [ecstatic_wilson] DSL2 - revision: c80285f8c8

    executor >  local (1)
    [c5/4c6ca9] sayHello (3)       [100%] 3 of 3, cached: 3 ✔
    [0e/6cbc59] convertToUpper (3) [100%] 3 of 3, cached: 3 ✔
    [02/61ead2] collectGreetings   [100%] 1 of 1 ✔
    ```

Si mireu al directori `results/hello_workflow/`, trobareu el nou fitxer d'informe, `trio-report.txt`.
Obriu-lo per verificar que el workflow ha informat correctament del recompte de salutacions que es van processar.

??? abstract "Contingut del fitxer"

    ```txt title="trio-report.txt"
    There were 3 greetings in this batch.
    ```

<figure class="excalidraw">
--8<-- "docs/en/docs/hello_nextflow/img/hello-collect-4-way.svg"
</figure>

Sentiu-vos lliures d'afegir més salutacions al CSV i provar què passa.

### Conclusió

Sabeu com fer que un procés emeti múltiples sortides amb nom i com gestionar-les adequadament a nivell de workflow.

Més generalment, enteneu els principis clau implicats en connectar processos junts de maneres comunes.

### Què segueix?

Preneu-vos un descans extra llarg, us l'heu guanyat.

Quan estigueu preparats, passeu a [**Part 4: Hello Modules**](./04_hello_modules.md) per aprendre com modularitzar el vostre codi per a una millor mantenibilitat i eficiència del codi.

---

## Qüestionari

<quiz>
Com accediu a la sortida d'un procés al bloc workflow?
- [ ] `process.output`
- [ ] `output.processName`
- [x] `processName.out`
- [ ] `get(processName)`

Més informació: [1.4. Passeu la sortida del primer procés al segon procés](#14-pass-the-output-of-the-first-process-to-the-second-process)
</quiz>

<quiz>
Què determina l'ordre d'execució dels processos a Nextflow?
- [ ] L'ordre en què els processos estan escrits al bloc workflow
- [ ] Ordre alfabètic per nom de procés
- [x] Dependències de dades entre processos
- [ ] Ordre aleatori per a execució paral·lela

Més informació: [1.4. Passeu la sortida del primer procés al segon procés](#14-pass-the-output-of-the-first-process-to-the-second-process)
</quiz>

<quiz>
Quin operador hauria de substituir `???` per recollir totes les sortides en una única llista per al procés posterior?

```groovy hl_lines="4"
workflow {
    greetings_ch = Channel.of('Hello', 'Bonjour', 'Hola')
    SAYHELLO(greetings_ch)
    GATHER_ALL(SAYHELLO.out.???)
}
```

- [ ] `flatten()`
- [x] `collect()`
- [ ] `mix()`
- [ ] `join()`

Més informació: [2.4. Utilitzeu un operador per recollir les salutacions en una única entrada](#24-use-an-operator-to-collect-the-greetings-into-a-single-input)
</quiz>

<quiz>
Quan hauríeu d'utilitzar l'operador `collect()`?
- [ ] Quan voleu processar elements en paral·lel
- [ ] Quan necessiteu filtrar el contingut del canal
- [x] Quan un procés posterior necessita tots els elements d'un procés anterior
- [ ] Quan voleu dividir dades entre múltiples processos

Més informació: [2.4. Utilitzeu un operador per recollir les salutacions en una única entrada](#24-use-an-operator-to-collect-the-greetings-into-a-single-input)
</quiz>

<quiz>
Com accediu a una sortida amb nom d'un procés?
- [ ] `processName.outputName`
- [ ] `processName.get(outputName)`
- [x] `processName.out.outputName`
- [ ] `output.processName.outputName`

Més informació: [4.1.2. Emeteu el fitxer d'informe i anomeneu les sortides](#412-emit-the-report-file-and-name-outputs)
</quiz>

<quiz>
Quina és la sintaxi correcta per anomenar una sortida en un procés?
- [ ] `name: outputName`
- [ ] `output: outputName`
- [x] `emit: outputName`
- [ ] `label: outputName`

Més informació: [4.1.2. Emeteu el fitxer d'informe i anomeneu les sortides](#412-emit-the-report-file-and-name-outputs)
</quiz>

<quiz>
Quan proporcioneu múltiples entrades a un procés, què ha de ser cert?
- [ ] Totes les entrades han de ser del mateix tipus
- [ ] Les entrades s'han de proporcionar en ordre alfabètic
- [x] L'ordre de les entrades ha de coincidir amb l'ordre definit al bloc d'entrada
- [ ] Només es poden proporcionar dues entrades alhora

Més informació: [3. Passeu més d'una entrada a un procés](#3-pass-more-than-one-input-to-a-process)
</quiz>
