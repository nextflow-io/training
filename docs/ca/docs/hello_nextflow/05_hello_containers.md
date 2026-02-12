# Part 5: Hello Containers

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Traducció assistida per IA - [més informació i suggeriments](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

<div class="video-wrapper">
  <iframe width="560" height="315" src="https://www.youtube.com/embed/Xqr--bKEN9U?si=y8lAedhEHWaTV4zd&amp;list=PLPZ8WHdZGxmWKozQuzr27jyMGqp9kElVK&amp;cc_load_policy=1" title="YouTube video player" frameborder="0" allow="accelerometer; autoplay; clipboard-write; encrypted-media; gyroscope; picture-in-picture; web-share" referrerpolicy="strict-origin-when-cross-origin" allowfullscreen></iframe>
</div>

/// caption
:fontawesome-brands-youtube:{ .youtube } Consulteu [la llista de reproducció completa](https://youtube.com/playlist?list=PLPZ8WHdZGxmWKozQuzr27jyMGqp9kElVK&si=eF7cLR62goy-lc6n) al canal de YouTube de Nextflow.

:green_book: La transcripció del vídeo està disponible [aquí](./transcripts/05_hello_containers.md).
///

A les Parts 1-4 d'aquest curs de formació, heu après a utilitzar els blocs de construcció bàsics de Nextflow per muntar un workflow senzill capaç de processar text, paral·lelitzar l'execució si hi havia múltiples entrades i recollir els resultats per a un processament posterior.

Tanmateix, estàveu limitats a eines UNIX bàsiques disponibles al vostre entorn.
Les tasques del món real sovint requereixen diverses eines i paquets que no s'inclouen per defecte.
Normalment, hauríeu d'instal·lar aquestes eines, gestionar les seves dependències i resoldre qualsevol conflicte.

Tot això és molt tediós i molest, així que us mostrarem com utilitzar **contenidors** per resoldre aquest problema de manera molt més convenient.

Un **contenidor** és una unitat de programari lleugera, autònoma i executable creada a partir d'una **imatge** de contenidor que inclou tot el necessari per executar una aplicació, incloent codi, biblioteques del sistema i configuració.
Com us podeu imaginar, això serà molt útil per fer els vostres pipelines més reproduïbles.

Tingueu en compte que ensenyarem això utilitzant [Docker](https://www.docker.com/get-started/), però recordeu que Nextflow també suporta [diverses altres tecnologies de contenidors](https://nextflow.io/docs/latest/container.html).

??? info "Com començar des d'aquesta secció"

    Aquesta secció del curs assumeix que heu completat les Parts 1-4 del curs [Hello Nextflow](./index.md) i teniu un pipeline complet i funcional.

    Si esteu començant el curs des d'aquest punt, haureu de copiar el directori `modules` des de les solucions:

    ```bash
    cp -r solutions/4-hello-modules/modules .
    ```

---

## 0. Escalfament: Executeu `hello-containers.nf`

Utilitzarem l'script de workflow `hello-containers.nf` com a punt de partida.
És equivalent a l'script produït en treballar la Part 4 d'aquest curs de formació, excepte que hem canviat les destinacions de sortida:

```groovy title="hello-containers.nf" linenums="37" hl_lines="3 7 11 15"
output {
    first_output {
        path 'hello_containers'
        mode 'copy'
    }
    uppercased {
        path 'hello_containers'
        mode 'copy'
    }
    collected {
        path 'hello_containers'
        mode 'copy'
    }
    batch_report {
        path 'hello_containers'
        mode 'copy'
    }
}
```

Només per assegurar-nos que tot funciona, executeu l'script una vegada abans de fer cap canvi:

```bash
nextflow run hello-containers.nf
```

??? success "Sortida de la comanda"

    ```console
     N E X T F L O W   ~  version 25.10.2

    Launching `hello-containers.nf` [nice_escher] DSL2 - revision: d5dfdc9872

    executor > local (7)
    [5a/ec1fa1] sayHello (2) [100%] 3 of 3 ✔
    [30/32b5b8] convertToUpper (3) [100%] 3 of 3 ✔
    [d3/be01bc] collectGreetings [100%] 1 of 1 ✔

    ```

Com abans, trobareu els fitxers de sortida al directori especificat al bloc `output` (`results/hello_containers/`).

??? abstract "Contingut del directori"

    ```console
    results/hello_containers/
    ├── Bonjour-output.txt
    ├── COLLECTED-batch-output.txt
    ├── Hello-output.txt
    ├── Holà-output.txt
    ├── batch-report.txt
    ├── UPPER-Bonjour-output.txt
    ├── UPPER-Hello-output.txt
    └── UPPER-Holà-output.txt
    ```

Si això us ha funcionat, esteu preparats per aprendre a utilitzar contenidors.

---

## 1. Utilitzeu un contenidor 'manualment'

El que volem fer és afegir un pas al nostre workflow que utilitzarà un contenidor per a l'execució.

Tanmateix, primer repassarem alguns conceptes i operacions bàsiques per consolidar la vostra comprensió del que són els contenidors abans de començar a utilitzar-los a Nextflow.

### 1.1. Descarregueu la imatge del contenidor

Per utilitzar un contenidor, normalment descarregueu o _pull_ una imatge de contenidor d'un registre de contenidors, i després executeu la imatge del contenidor per crear una instància de contenidor.

La sintaxi general és la següent:

```bash title="Syntax"
docker pull '<container>'
```

La part `docker pull` és la instrucció al sistema de contenidors per descarregar una imatge de contenidor d'un repositori.

La part `'<container>'` és l'adreça URI de la imatge del contenidor.

Com a exemple, descarreguem una imatge de contenidor que conté [cowpy](https://github.com/jeffbuttars/cowpy), una implementació en python d'una eina anomenada `cowsay` que genera art ASCII per mostrar entrades de text arbitràries d'una manera divertida.

```txt title="Example"
 ________________________
< Are we having fun yet? >
 ------------------------
    \                                  ___-------___
     \                             _-~~             ~~-_
      \                         _-~                    /~-_
             /^\__/^\         /~  \                   /    \
           /|  O|| O|        /      \_______________/        \
          | |___||__|      /       /                \          \
          |          \    /      /                    \          \
          |   (_______) /______/                        \_________ \
          |         / /         \                      /            \
           \         \^\\         \                  /               \     /
             \         ||           \______________/      _-_       //\__//
               \       ||------_-~~-_ ------------- \ --/~   ~\    || __/
                 ~-----||====/~     |==================|       |/~~~~~
                  (_(__/  ./     /                    \_\      \.
                         (_(___/                         \_____)_)
```

Hi ha diversos repositoris on podeu trobar contenidors publicats.
Hem utilitzat el servei [Seqera Containers](https://seqera.io/containers/) per generar aquesta imatge de contenidor Docker a partir del paquet Conda `cowpy`: `'community.wave.seqera.io/library/cowpy:1.1.5--3db457ae1977a273'`.

Executeu la comanda pull completa:

```bash
docker pull 'community.wave.seqera.io/library/cowpy:1.1.5--3db457ae1977a273'
```

??? success "Sortida de la comanda"

    ```console
    1.1.5--3db457ae1977a273: Pulling from library/cowpy
    dafa2b0c44d2: Pull complete
    dec6b097362e: Pull complete
    f88da01cff0b: Pull complete
    4f4fb700ef54: Pull complete
    92dc97a3ef36: Pull complete
    403f74b0f85e: Pull complete
    10b8c00c10a5: Pull complete
    17dc7ea432cc: Pull complete
    bb36d6c3110d: Pull complete
    0ea1a16bbe82: Pull complete
    030a47592a0a: Pull complete
    c23bdb422167: Pull complete
    e1686ff32a11: Pull complete
    Digest: sha256:1ebc0043e8cafa61203bf42d29fd05bd14e7b4298e5e8cf986504c15f5aa4160
    Status: Downloaded newer image for community.wave.seqera.io/library/cowpy:1.1.5--3db457ae1977a273
    community.wave.seqera.io/library/cowpy:1.1.5--3db457ae1977a273
    ```

Si mai no heu descarregat la imatge abans, això pot trigar un minut a completar-se.
Un cop fet, teniu una còpia local de la imatge del contenidor.

### 1.2. Utilitzeu el contenidor per executar `cowpy` com a comanda única

Una manera molt comuna en què la gent utilitza contenidors és executar-los directament, _és a dir_, de manera no interactiva.
Això és genial per executar comandes úniques.

La sintaxi general és la següent:

```bash title="Syntax"
docker run --rm '<container>' [tool command]
```

La part `docker run --rm '<container>'` és la instrucció al sistema de contenidors per crear una instància de contenidor a partir d'una imatge de contenidor i executar una comanda en ella.
La bandera `--rm` indica al sistema que tanqui la instància del contenidor després que la comanda s'hagi completat.

La sintaxi `[tool command]` depèn de l'eina que esteu utilitzant i de com està configurat el contenidor.
Comencem simplement amb `cowpy`.

Completament muntada, la comanda d'execució del contenidor es veu així; endavant, executeu-la.

```bash
docker run --rm 'community.wave.seqera.io/library/cowpy:1.1.5--3db457ae1977a273' cowpy
```

??? success "Sortida de la comanda"

    ```console
    ______________________________________________________
    < Cowacter, eyes:default, tongue:False, thoughts:False >
    ------------------------------------------------------
        \   ^__^
          \  (oo)\_______
            (__)\       )\/\
              ||----w |
              ||     ||
    ```

El sistema ha creat el contenidor, ha executat la comanda `cowpy` amb els seus paràmetres, ha enviat la sortida a la consola i finalment, ha tancat la instància del contenidor.

### 1.3. Utilitzeu el contenidor per executar `cowpy` interactivament

També podeu executar un contenidor interactivament, cosa que us dóna un prompt de shell dins del contenidor i us permet jugar amb la comanda.

#### 1.3.1. Inicieu el contenidor

Per executar interactivament, només afegim `-it` a la comanda `docker run`.
Opcionalment, podem especificar el shell que volem utilitzar dins del contenidor afegint _p. ex._ `/bin/bash` a la comanda.

```bash
docker run --rm -it 'community.wave.seqera.io/library/cowpy:1.1.5--3db457ae1977a273' /bin/bash
```

Noteu que el vostre prompt canvia a alguna cosa com `(base) root@b645838b3314:/tmp#`, que indica que ara esteu dins del contenidor.

Podeu verificar-ho executant `ls /` per llistar el contingut del directori des de l'arrel del sistema de fitxers:

```bash
ls /
```

??? abstract "Sortida de la comanda"

    ```console
    bin  boot  dev  etc  home  lib  lib64  media  mnt  opt  proc  root  run  sbin  srv  sys  tmp  usr  var
    ```

Utilitzem `ls` aquí en lloc de `tree` perquè la utilitat `tree` no està disponible en aquest contenidor.
Podeu veure que el sistema de fitxers dins del contenidor és diferent del sistema de fitxers del vostre sistema amfitrió.

Una limitació del que acabem de fer és que el contenidor està completament aïllat del sistema amfitrió per defecte.
Això significa que el contenidor no pot accedir a cap fitxer del sistema amfitrió tret que li permeteu fer-ho explícitament.

Us mostrarem com fer-ho en un minut.

#### 1.3.2. Executeu la comanda de l'eina desitjada

Ara que esteu dins del contenidor, podeu executar la comanda `cowpy` directament i donar-li alguns paràmetres.
Per exemple, la documentació de l'eina diu que podem canviar el personatge ('cowacter') amb `-c`.

```bash
cowpy "Hello Containers" -c tux
```

??? success "Sortida de la comanda"

    ```console
    __________________
    < Hello Containers >
    ------------------
      \
        \
            .--.
          |o_o |
          |:_/ |
          //   \ \
        (|     | )
        /'\_   _/`\
        \___)=(___/
    ```

Ara la sortida mostra el pingüí de Linux, Tux, en lloc de la vaca per defecte, perquè hem especificat el paràmetre `-c tux`.

Com que esteu dins del contenidor, podeu executar la comanda `cowpy` tantes vegades com vulgueu, variant els paràmetres d'entrada, sense haver de preocupar-vos per les comandes de Docker.

!!! Tip "Consell"

    Utilitzeu la bandera '-c' per triar un personatge diferent, incloent:
    `beavis`, `cheese`, `daemon`, `dragonandcow`, `ghostbusters`, `kitty`, `moose`, `milk`, `stegosaurus`, `turkey`, `turtle`, `tux`

Això és genial. El que seria encara més genial és si poguéssim alimentar el nostre `greetings.csv` com a entrada en això.
Però com que no tenim accés al sistema de fitxers, no podem.

Solucionem això.

#### 1.3.3. Sortiu del contenidor

Per sortir del contenidor, podeu escriure `exit` al prompt o utilitzar la drecera de teclat ++ctrl+d++.

```bash
exit
```

El vostre prompt ara hauria de tornar al que era abans d'iniciar el contenidor.

#### 1.3.4. Munteu dades al contenidor

Com s'ha indicat anteriorment, el contenidor està aïllat del sistema amfitrió per defecte.

Per permetre que el contenidor accedeixi al sistema de fitxers de l'amfitrió, podeu **muntar** un **volum** del sistema amfitrió al contenidor utilitzant la següent sintaxi:

```bash title="Syntax"
-v <outside_path>:<inside_path>
```

En el nostre cas `<outside_path>` serà el directori de treball actual, així que podem utilitzar simplement un punt (`.`), i `<inside_path>` és només un àlies que inventem; anomenem-lo `/my_project` (el camí intern ha de ser absolut).

Per muntar un volum, reemplacem els camins i afegim l'argument de muntatge de volum a la comanda docker run de la següent manera:

```bash
docker run --rm -it -v .:/my_project 'community.wave.seqera.io/library/cowpy:1.1.5--3db457ae1977a273' /bin/bash
```

Això munta el directori de treball actual com un volum que serà accessible sota `/my_project` dins del contenidor.

Podeu comprovar que funciona llistant el contingut de `/my_project`:

```bash
ls /my_project
```

??? success "Sortida de la comanda"

    ```console
    data               hello-config.nf      hello-modules.nf   hello-world.nf  nextflow.config  solutions         work
    hello-channels.nf  hello-containers.nf  hello-workflow.nf  modules         results          test-params.json
    ```

Ara podeu veure el contingut del directori de treball des de dins del contenidor, incloent el fitxer `greetings.csv` sota `data/`.

Això ha establert efectivament un túnel a través de la paret del contenidor que podeu utilitzar per accedir a aquesta part del vostre sistema de fitxers.

#### 1.3.5. Utilitzeu les dades muntades

Ara que hem muntat el directori de treball al contenidor, podem utilitzar la comanda `cowpy` per mostrar el contingut del fitxer `greetings.csv`.

Per fer això, utilitzarem `cat /my_project/data/greetings.csv | ` per canalitzar el contingut del fitxer CSV a la comanda `cowpy`.

```bash
cat /my_project/data/greetings.csv | cowpy -c turkey
```

??? success "Sortida de la comanda"

    ```console title="data/greetings.csv"
     ____________________
    / Hello,English,123  \
    | Bonjour,French,456 |
    \ Holà,Spanish,789   /
    --------------------
      \                                  ,+*^^*+___+++_
      \                           ,*^^^^              )
        \                       _+*                     ^**+_
        \                    +^       _ _++*+_+++_,         )
                  _+^^*+_    (     ,+*^ ^          \+_        )
                {       )  (    ,(    ,_+--+--,      ^)      ^\
                { (\@)    } f   ,(  ,+-^ __*_*_  ^^\_   ^\       )
              {:;-/    (_+*-+^^^^^+*+*<_ _++_)_    )    )      /
              ( /  (    (        ,___    ^*+_+* )   <    <      \
              U _/     )    *--<  ) ^\-----++__)   )    )       )
                (      )  _(^)^^))  )  )\^^^^^))^*+/    /       /
              (      /  (_))_^)) )  )  ))^^^^^))^^^)__/     +^^
            (     ,/    (^))^))  )  ) ))^^^^^^^))^^)       _)
              *+__+*       (_))^)  ) ) ))^^^^^^))^^^^^)____*^
              \             \_)^)_)) ))^^^^^^^^^^))^^^^)
              (_             ^\__^^^^^^^^^^^^))^^^^^^^)
                ^\___            ^\__^^^^^^))^^^^^^^^)\\
                      ^^^^^\uuu/^^\uuu/^^^^\^\^\^\^\^\^\^\
                        ___) >____) >___   ^\_\_\_\_\_\_\)
                        ^^^//\\_^^//\\_^       ^(\_\_\_\)
                          ^^^ ^^ ^^^ ^
    ```

Això produeix l'art ASCII desitjat d'un gall dindi recitant les nostres salutacions d'exemple!
Excepte que aquí el gall dindi està repetint les files completes en lloc de només les salutacions.
Ja sabem que el nostre workflow de Nextflow farà una feina millor!

Sentiu-vos lliures de jugar amb aquesta comanda.
Quan hàgiu acabat, sortiu del contenidor com abans:

```bash
exit
```

Us trobareu de nou al vostre shell normal.

### Conclusió

Sabeu com descarregar un contenidor i executar-lo ja sigui com a comanda única o interactivament. També sabeu com fer les vostres dades accessibles des de dins del vostre contenidor, cosa que us permet provar qualsevol eina que us interessi amb dades reals sense haver d'instal·lar cap programari al vostre sistema.

### Què segueix?

Apreneu a utilitzar contenidors per a l'execució de processos de Nextflow.

---

## 2. Utilitzeu contenidors a Nextflow

Nextflow té suport integrat per executar processos dins de contenidors per permetre-us executar eines que no teniu instal·lades al vostre entorn de càlcul.
Això significa que podeu utilitzar qualsevol imatge de contenidor que vulgueu per executar els vostres processos, i Nextflow s'encarregarà de descarregar la imatge, muntar les dades i executar el procés dins d'ella.

Per demostrar això, afegirem un pas `cowpy` al pipeline que hem estat desenvolupant, després del pas `collectGreetings`.

<figure class="excalidraw">
--8<-- "docs/en/docs/hello_nextflow/img/hello-pipeline-cowpy.svg"
</figure>

### 2.1. Escriviu un mòdul `cowpy`

Primer, creem el mòdul de procés `cowpy`.

#### 2.1.1. Creeu un fitxer buit per al nou mòdul

Creeu un fitxer buit per al mòdul anomenat `cowpy.nf`.

```bash
touch modules/cowpy.nf
```

Això ens dóna un lloc on posar el codi del procés.

#### 2.1.2. Copieu el codi del procés `cowpy` al fitxer del mòdul

Podem modelar el nostre procés `cowpy` en els altres processos que hem escrit anteriorment.

```groovy title="modules/cowpy.nf" linenums="1"
#!/usr/bin/env nextflow

// Genera art ASCII amb cowpy
process cowpy {

    input:
    path input_file
    val character

    output:
    path "cowpy-${input_file}"

    script:
    """
    cat ${input_file} | cowpy -c "${character}" > cowpy-${input_file}
    """

}
```

El procés espera un `input_file` que conté les salutacions així com un valor `character`.

La sortida serà un nou fitxer de text que conté l'art ASCII generat per l'eina `cowpy`.

### 2.2. Afegiu cowpy al workflow

Ara necessitem importar el mòdul i cridar el procés.

#### 2.2.1. Importeu el procés `cowpy` a `hello-containers.nf`

Inseriu la declaració d'importació sobre el bloc de workflow i ompliu-la adequadament.

=== "Després"

    ```groovy title="hello-containers.nf" linenums="3" hl_lines="5"
    // Inclou mòduls
    include { sayHello } from './modules/sayHello.nf'
    include { convertToUpper } from './modules/convertToUpper.nf'
    include { collectGreetings } from './modules/collectGreetings.nf'
    include { cowpy } from './modules/cowpy.nf'
    ```

=== "Abans"

    ```groovy title="hello-containers.nf" linenums="3"
    // Inclou mòduls
    include { sayHello } from './modules/sayHello.nf'
    include { convertToUpper } from './modules/convertToUpper.nf'
    include { collectGreetings } from './modules/collectGreetings.nf'
    ```

Ara el mòdul `cowpy` està disponible per utilitzar al workflow.

#### 2.2.2. Afegiu una crida al procés `cowpy` al workflow

Connectem el procés `cowpy()` a la sortida del procés `collectGreetings()`, que com recordareu produeix dues sortides:

- `collectGreetings.out.outfile` conté el fitxer de sortida <--_el que volem_
- `collectGreetings.out.report` conté el fitxer d'informe amb el recompte de salutacions per lot

Al bloc de workflow, feu el següent canvi de codi:

=== "Després"

    ```groovy title="hello-containers.nf" linenums="19" hl_lines="12-13"
        main:
        // crea un canal per a entrades des d'un fitxer CSV
        greeting_ch = channel.fromPath(params.input)
                            .splitCsv()
                            .map { line -> line[0] }
        // emet una salutacio
        sayHello(greeting_ch)
        // converteix la salutacio a majuscules
        convertToUpper(sayHello.out)
        // recull totes les salutacions en un fitxer
        collectGreetings(convertToUpper.out.collect(), params.batch)
        // genera art ASCII de les salutacions amb cowpy
        cowpy(collectGreetings.out.outfile, params.character)
    ```

=== "Abans"

    ```groovy title="hello-containers.nf" linenums="19"
        main:
        // crea un canal per a entrades des d'un fitxer CSV
        greeting_ch = channel.fromPath(params.input)
                            .splitCsv()
                            .map { line -> line[0] }
        // emet una salutacio
        sayHello(greeting_ch)
        // converteix la salutacio a majuscules
        convertToUpper(sayHello.out)
        // recull totes les salutacions en un fitxer
        collectGreetings(convertToUpper.out.collect(), params.batch)
    ```

Noteu que hem declarat un nou paràmetre CLI, `params.character`, per especificar quin personatge volem que digui les salutacions.

#### 2.2.3. Afegiu el paràmetre `character` al bloc `params`

Això és tècnicament opcional però és la pràctica recomanada i és una oportunitat per establir un valor per defecte per al personatge mentre hi som.

=== "Després"

    ```groovy title="hello-containers.nf" linenums="9" hl_lines="7"
    /*
    * Paràmetres del pipeline
    */
    params {
        input: Path = 'data/greetings.csv'
        batch: String = 'batch'
        character: String = 'turkey'
    }
    ```

=== "Abans"

    ```groovy title="hello-containers.nf" linenums="9"
    /*
    * Paràmetres del pipeline
    */
    params {
        input: Path = 'data/greetings.csv'
        batch: String = 'batch'
    }
    ```

Ara podem ser mandrosos i saltar-nos escriure el paràmetre de personatge a les nostres línies de comandes.

#### 2.2.4. Actualitzeu les sortides del workflow

Necessitem actualitzar les sortides del workflow per publicar la sortida del procés `cowpy`.

##### 2.2.4.1. Actualitzeu la secció `publish:`

Al `bloc de workflow`, feu el següent canvi de codi:

=== "Després"

    ```groovy title="hello-containers.nf" linenums="34" hl_lines="6"
        publish:
        first_output = sayHello.out
        uppercased = convertToUpper.out
        collected = collectGreetings.out.outfile
        batch_report = collectGreetings.out.report
        cowpy_art = cowpy.out
    ```

=== "Abans"

    ```groovy title="hello-containers.nf" linenums="34"
        publish:
        first_output = sayHello.out
        uppercased = convertToUpper.out
        collected = collectGreetings.out.outfile
        batch_report = collectGreetings.out.report
    ```

El procés `cowpy` només produeix una sortida així que podem referir-nos-hi de la manera habitual afegint `.out`.

Però per ara, acabem d'actualitzar les sortides a nivell de workflow.

##### 2.2.4.2. Actualitzeu el bloc `output`

Necessitem afegir la sortida final `cowpy_art` al bloc `output`. Mentre hi som, editem també les destinacions de publicació ja que ara el nostre pipeline està complet i sabem quines sortides realment ens importen.

Al bloc `output`, feu els següents canvis de codi:

=== "Després"

    ```groovy title="hello-containers.nf" linenums="42" hl_lines="3 7 11 15 18-21"
    output {
        first_output {
            path 'hello_containers/intermediates'
            mode 'copy'
        }
        uppercased {
            path 'hello_containers/intermediates'
            mode 'copy'
        }
        collected {
            path 'hello_containers/intermediates'
            mode 'copy'
        }
        batch_report {
            path 'hello_containers'
            mode 'copy'
        }
        cowpy_art {
            path 'hello_containers'
            mode 'copy'
        }
    }
    ```

=== "Abans"

    ```groovy title="hello-containers.nf" linenums="42" hl_lines="3 7 11 15"
    output {
        first_output {
            path 'hello_containers'
            mode 'copy'
        }
        uppercased {
            path 'hello_containers'
            mode 'copy'
        }
        collected {
            path 'hello_containers'
            mode 'copy'
        }
        batch_report {
            path 'hello_containers'
            mode 'copy'
        }
    }
    ```

Ara les sortides publicades estaran una mica més organitzades.

#### 2.2.5. Executeu el workflow

Només per recapitular, això és el que estem buscant:

<figure class="excalidraw">
--8<-- "docs/en/docs/hello_nextflow/img/hello_pipeline_complete.svg"
</figure>

Creieu que funcionarà?

Esborrem les sortides publicades anteriors per tenir una pissarra neta, i executem el workflow amb la bandera `-resume`.

```bash
rm -r hello_containers/
nextflow run hello-containers.nf -resume
```

??? failure "Sortida de la comanda (editada per claredat)"

    ```console hl_lines="10 13 20-21 26-27"
     N E X T F L O W   ~  version 25.10.2

    Launching `hello-containers.nf` [lonely_woese] DSL2 - revision: abf1dccf7f

    executor >  local (1)
    [c9/f5c686] sayHello (3)       [100%] 3 of 3, cached: 3 ✔
    [ef/3135a8] convertToUpper (3) [100%] 3 of 3, cached: 3 ✔
    [7f/f435e3] collectGreetings   [100%] 1 of 1, cached: 1 ✔
    [9b/02e776] cowpy              [  0%] 0 of 1 ✘
    ERROR ~ Error executing process > 'cowpy'

    Caused by:
      Process `cowpy` terminated with an error exit status (127)


    Command executed:

      cat COLLECTED-batch-output.txt | cowpy -c "turkey" > cowpy-COLLECTED-batch-output.txt

    Command exit status:
      127

    Command output:
      (empty)

    Command error:
      .command.sh: line 2: cowpy: command not found

    Work dir:
      /workspaces/training/hello-nextflow/work/9b/02e7761db848f82db3c3e59ff3a9b6

    Tip: when you have fixed the problem you can continue the execution adding the option `-resume` to the run command line

    -- Check '.nextflow.log' file for details
    ERROR ~ Cannot access first() element from an empty List

    -- Check '.nextflow.log' file for details
    ```

Oh no, hi ha un error!
El codi d'error donat per `error exit status (127)` significa que l'executable que hem demanat no s'ha trobat.

Això té sentit, ja que estem cridant l'eina `cowpy` però encara no hem especificat un contenidor (ups).

### 2.3. Utilitzeu un contenidor per executar el procés `cowpy`

Necessitem especificar un contenidor i dir a Nextflow que l'utilitzi per al procés `cowpy()`.

#### 2.3.1. Especifiqueu un contenidor per a `cowpy`

Podem utilitzar la mateixa imatge que estàvem utilitzant directament a la primera secció d'aquest tutorial.

Editeu el mòdul `cowpy.nf` per afegir la directiva `container` a la definició del procés de la següent manera:

=== "Després"

    ```groovy title="modules/cowpy.nf" linenums="4" hl_lines="3"
    process cowpy {

        container 'community.wave.seqera.io/library/cowpy:1.1.5--3db457ae1977a273'

        input:
        path input_file
        val character

        output:
        path "cowpy-${input_file}"

        script:
        """
        cat ${input_file} | cowpy -c "${character}" > cowpy-${input_file}
        """
    }
    ```

=== "Abans"

    ```groovy title="modules/cowpy.nf" linenums="4"
    process cowpy {

        input:
        path input_file
        val character

        output:
        path "cowpy-${input_file}"

        script:
        """
        cat ${input_file} | cowpy -c "${character}" > cowpy-${input_file}
        """
    }
    ```

Això indica a Nextflow que _si l'ús de Docker està habilitat_, hauria d'utilitzar la imatge de contenidor especificada aquí per executar el procés.

#### 2.3.2. Habiliteu l'ús de Docker mitjançant el fitxer `nextflow.config`

Noteu que hem dit _'si l'ús de Docker està habilitat'_. Per defecte, no ho està, així que necessitem dir a Nextflow que se li permet utilitzar Docker.
Per això, anticiparem lleugerament el tema de la següent i última part d'aquest curs (Part 6), que cobreix la configuració.

Una de les principals maneres que Nextflow ofereix per configurar l'execució del workflow és utilitzar un fitxer `nextflow.config`.
Quan aquest fitxer està present al directori actual, Nextflow el carregarà automàticament i aplicarà qualsevol configuració que contingui.

Hem proporcionat un fitxer `nextflow.config` amb una sola línia de codi que deshabilita explícitament Docker: `docker.enabled = false`.

Ara, canviem això a `true` per habilitar Docker:

=== "Després"

    ```console title="nextflow.config" linenums="1" hl_lines="1"
    docker.enabled = true
    ```

=== "Abans"

    ```console title="nextflow.config" linenums="1" hl_lines="1"
    docker.enabled = false
    ```

!!! tip "Consell"

    És possible habilitar l'execució de Docker des de la línia de comandes, per execució, utilitzant el paràmetre `-with-docker <container>`.
    Tanmateix, això només ens permet especificar un contenidor per a tot el workflow, mentre que l'enfocament que acabem de mostrar-vos ens permet especificar un contenidor diferent per procés.
    Això és millor per a la modularitat, el manteniment del codi i la reproduïbilitat.

#### 2.3.3. Executeu el workflow amb Docker habilitat

Executeu el workflow amb la bandera `-resume`:

```bash
nextflow run hello-containers.nf -resume
```

??? success "Sortida de la comanda"

    ```console
     N E X T F L O W   ~  version 25.10.2

    Launching `hello-containers.nf` [drunk_perlman] DSL2 - revision: abf1dccf7f

    executor >  local (1)
    [c9/f5c686] sayHello (3)       [100%] 3 of 3, cached: 3 ✔
    [ef/3135a8] convertToUpper (3) [100%] 3 of 3, cached: 3 ✔
    [7f/f435e3] collectGreetings   [100%] 1 of 1, cached: 1 ✔
    [98/656c6c] cowpy              [100%] 1 of 1 ✔
    ```

Aquesta vegada sí que funciona!
Com sempre podeu trobar les sortides del workflow al directori de resultats corresponent, encara que aquesta vegada estan una mica més ordenades, amb només l'informe i la sortida final al nivell superior, i tots els fitxers intermedis apartats en un subdirectori.

??? abstract "Contingut del directori"

    ```console
    results/hello_containers/
    ├── cowpy-COLLECTED-batch-output.txt
    ├── intermediates
    │   ├── Bonjour-output.txt
    │   ├── COLLECTED-batch-output.txt
    │   ├── Hello-output.txt
    │   ├── Holà-output.txt
    │   ├── UPPER-Bonjour-output.txt
    │   ├── UPPER-Hello-output.txt
    │   └── UPPER-Holà-output.txt
    └── batch-report.txt
    ```

La sortida d'art ASCII final està al directori `results/hello_containers/`, amb el nom `cowpy-COLLECTED-batch-output.txt`.

??? abstract "Contingut del fitxer"

    ```console title="results/hello_containers/cowpy-COLLECTED-batch-output.txt"
    _________
    / HOLà    \
    | HELLO   |
    \ BONJOUR /
    ---------
      \                                  ,+*^^*+___+++_
      \                           ,*^^^^              )
        \                       _+*                     ^**+_
        \                    +^       _ _++*+_+++_,         )
                  _+^^*+_    (     ,+*^ ^          \+_        )
                {       )  (    ,(    ,_+--+--,      ^)      ^\
                { (\@)    } f   ,(  ,+-^ __*_*_  ^^\_   ^\       )
              {:;-/    (_+*-+^^^^^+*+*<_ _++_)_    )    )      /
              ( /  (    (        ,___    ^*+_+* )   <    <      \
              U _/     )    *--<  ) ^\-----++__)   )    )       )
                (      )  _(^)^^))  )  )\^^^^^))^*+/    /       /
              (      /  (_))_^)) )  )  ))^^^^^))^^^)__/     +^^
            (     ,/    (^))^))  )  ) ))^^^^^^^))^^)       _)
              *+__+*       (_))^)  ) ) ))^^^^^^))^^^^^)____*^
              \             \_)^)_)) ))^^^^^^^^^^))^^^^)
              (_             ^\__^^^^^^^^^^^^))^^^^^^^)
                ^\___            ^\__^^^^^^))^^^^^^^^)\\
                      ^^^^^\uuu/^^\uuu/^^^^\^\^\^\^\^\^\^\
                        ___) >____) >___   ^\_\_\_\_\_\_\)
                        ^^^//\\_^^//\\_^       ^(\_\_\_\)
                          ^^^ ^^ ^^^ ^
    ```

I aquí està, el nostre bonic gall dindi dient les salutacions com desitjàvem.

#### 2.3.4. Inspeccioneu com Nextflow va llançar la tasca contenidoritzada

Com a coda final d'aquesta secció, donem una ullada al subdirectori de treball per a una de les crides del procés `cowpy` per obtenir una mica més d'informació sobre com Nextflow treballa amb contenidors sota el capó.

Comproveu la sortida de la vostra comanda `nextflow run` per trobar el camí al subdirectori de treball per al procés `cowpy`.
Mirant el que hem obtingut per a l'execució mostrada anteriorment, la línia de registre de consola per al procés `cowpy` comença amb `[98/656c6c]`.
Això correspon al següent camí de directori truncat: `work/98/656c6c`.

En aquest directori, trobareu el fitxer `.command.run` que conté totes les comandes que Nextflow va executar en nom vostre en el curs d'executar el pipeline.

??? abstract "Contingut del fitxer"

    ```console title="work/98/656c6c90cce1667c094d880f4b6dcc/.command.run"
    #!/bin/bash
    ### ---
    ### name: 'cowpy'
    ### container: 'community.wave.seqera.io/library/cowpy:1.1.5--3db457ae1977a273'
    ### outputs:
    ### - 'cowpy-COLLECTED-batch-output.txt'
    ### ...
    set -e
    set -u
    NXF_DEBUG=${NXF_DEBUG:=0}; [[ $NXF_DEBUG > 1 ]] && set -x
    NXF_ENTRY=${1:-nxf_main}


    nxf_sleep() {
      sleep $1 2>/dev/null || sleep 1;
    }

    nxf_date() {
        local ts=$(date +%s%3N);
        if [[ ${#ts} == 10 ]]; then echo ${ts}000
        elif [[ $ts == *%3N ]]; then echo ${ts/\%3N/000}
        elif [[ $ts == *3N ]]; then echo ${ts/3N/000}
        elif [[ ${#ts} == 13 ]]; then echo $ts
        else echo "Unexpected timestamp value: $ts"; exit 1
        fi
    }

    nxf_env() {
        echo '============= task environment ============='
        env | sort | sed "s/\(.*\)AWS\(.*\)=\(.\{6\}\).*/\1AWS\2=\3xxxxxxxxxxxxx/"
        echo '============= task output =================='
    }

    nxf_kill() {
        declare -a children
        while read P PP;do
            children[$PP]+=" $P"
        done < <(ps -e -o pid= -o ppid=)

        kill_all() {
            [[ $1 != $$ ]] && kill $1 2>/dev/null || true
            for i in ${children[$1]:=}; do kill_all $i; done
        }

        kill_all $1
    }

    nxf_mktemp() {
        local base=${1:-/tmp}
        mkdir -p "$base"
        if [[ $(uname) = Darwin ]]; then mktemp -d $base/nxf.XXXXXXXXXX
        else TMPDIR="$base" mktemp -d -t nxf.XXXXXXXXXX
        fi
    }

    nxf_fs_copy() {
      local source=$1
      local target=$2
      local basedir=$(dirname $1)
      mkdir -p $target/$basedir
      cp -fRL $source $target/$basedir
    }

    nxf_fs_move() {
      local source=$1
      local target=$2
      local basedir=$(dirname $1)
      mkdir -p $target/$basedir
      mv -f $source $target/$basedir
    }

    nxf_fs_rsync() {
      rsync -rRl $1 $2
    }

    nxf_fs_rclone() {
      rclone copyto $1 $2/$1
    }

    nxf_fs_fcp() {
      fcp $1 $2/$1
    }

    on_exit() {
        local last_err=$?
        local exit_status=${nxf_main_ret:=0}
        [[ ${exit_status} -eq 0 && ${nxf_unstage_ret:=0} -ne 0 ]] && exit_status=${nxf_unstage_ret:=0}
        [[ ${exit_status} -eq 0 && ${last_err} -ne 0 ]] && exit_status=${last_err}
        printf -- $exit_status > /workspaces/training/hello-nextflow/work/98/656c6c90cce1667c094d880f4b6dcc/.exitcode
        set +u
        docker rm $NXF_BOXID &>/dev/null || true
        exit $exit_status
    }

    on_term() {
        set +e
        docker stop $NXF_BOXID
    }

    nxf_launch() {
        docker run -i --cpu-shares 1024 -e "NXF_TASK_WORKDIR" -v /workspaces/training/hello-nextflow/work:/workspaces/training/hello-nextflow/work -w "$NXF_TASK_WORKDIR" --name $NXF_BOXID community.wave.seqera.io/library/cowpy:1.1.5--3db457ae1977a273 /bin/bash -ue /workspaces/training/hello-nextflow/work/98/656c6c90cce1667c094d880f4b6dcc/.command.sh
    }

    nxf_stage() {
        true
        # stage input files
        rm -f COLLECTED-batch-output.txt
        ln -s /workspaces/training/hello-nextflow/work/7f/f435e3f2cf95979b5f3d7647ae6696/COLLECTED-batch-output.txt COLLECTED-batch-output.txt
    }

    nxf_unstage_outputs() {
        true
    }

    nxf_unstage_controls() {
        true
    }

    nxf_unstage() {
        if [[ ${nxf_main_ret:=0} == 0 ]]; then
            (set -e -o pipefail; (nxf_unstage_outputs | tee -a .command.out) 3>&1 1>&2 2>&3 | tee -a .command.err)
            nxf_unstage_ret=$?
        fi
        nxf_unstage_controls
    }

    nxf_main() {
        trap on_exit EXIT
        trap on_term TERM INT USR2
        trap '' USR1

        [[ "${NXF_CHDIR:-}" ]] && cd "$NXF_CHDIR"
        export NXF_BOXID="nxf-$(dd bs=18 count=1 if=/dev/urandom 2>/dev/null | base64 | tr +/ 0A | tr -d '\r\n')"
        NXF_SCRATCH=''
        [[ $NXF_DEBUG > 0 ]] && nxf_env
        touch /workspaces/training/hello-nextflow/work/98/656c6c90cce1667c094d880f4b6dcc/.command.begin
        set +u
        set -u
        [[ $NXF_SCRATCH ]] && cd $NXF_SCRATCH
        export NXF_TASK_WORKDIR="$PWD"
        nxf_stage

        set +e
        (set -o pipefail; (nxf_launch | tee .command.out) 3>&1 1>&2 2>&3 | tee .command.err) &
        pid=$!
        wait $pid || nxf_main_ret=$?
        nxf_unstage
    }

    $NXF_ENTRY

    ```

Si cerqueu `nxf_launch` en aquest fitxer, hauríeu de veure alguna cosa com això:

```console
nxf_launch() {
    docker run -i --cpu-shares 1024 -e "NXF_TASK_WORKDIR" -v /workspaces/training/hello-nextflow/work:/workspaces/training/hello-nextflow/work -w "$NXF_TASK_WORKDIR" --name $NXF_BOXID community.wave.seqera.io/library/cowpy:1.1.5--3db457ae1977a273 /bin/bash -ue /workspaces/training/hello-nextflow/work/98/656c6c90cce1667c094d880f4b6dcc/.command.sh
}
```

Com podeu veure, Nextflow està utilitzant la comanda `docker run` per llançar la crida del procés.
També munta el subdirectori de treball corresponent al contenidor, estableix el directori de treball dins del contenidor en conseqüència, i executa el nostre script bash amb plantilla al fitxer `.command.sh`.

Tota la feina dura que havíem de fer manualment a la primera secció? Nextflow ho fa per nosaltres entre bastidors!

```txt
 _______________________
< Hurray for robots...! >
 -----------------------
                                   ,-----.
                                   |     |
                                ,--|     |-.
                         __,----|  |     | |
                       ,;::     |  `_____' |
                       `._______|    i^i   |
                                `----| |---'| .
                           ,-------._| |== ||//
                           |       |_|P`.  /'/
                           `-------' 'Y Y/'/'
                                     .==\ /_\
   ^__^                             /   /'|  `i
   (oo)\_______                   /'   /  |   |
   (__)\       )\/\             /'    /   |   `i
       ||----w |           ___,;`----'.___L_,-'`\__
       ||     ||          i_____;----\.____i""\____\
```

### Conclusió

Sabeu com utilitzar contenidors a Nextflow per executar processos.

### Què segueix?

Feu una pausa!

Quan estigueu preparats, passeu a [**Part 6: Hello Config**](./06_hello_config.md) per aprendre a configurar l'execució del vostre pipeline per adaptar-se a la vostra infraestructura així com gestionar la configuració d'entrades i paràmetres.

És l'última part, i després haureu acabat aquest curs!

---

## Qüestionari

<quiz>
Què és un contenidor?
- [ ] Un tipus de màquina virtual
- [ ] Un format de compressió de fitxers
- [x] Una unitat executable lleugera i autònoma que inclou tot el necessari per executar una aplicació
- [ ] Un protocol de xarxa
</quiz>

<quiz>
Quina és la diferència entre una imatge de contenidor i una instància de contenidor?
- [ ] Són la mateixa cosa
- [x] Una imatge és una plantilla; una instància és un contenidor en execució creat a partir d'aquesta imatge
- [ ] Una instància és una plantilla; una imatge és un contenidor en execució
- [ ] Les imatges són per a Docker; les instàncies són per a Singularity
</quiz>

<quiz>
Què fa la bandera `-v` en una comanda `docker run`?
- [ ] Habilita la sortida detallada
- [ ] Valida el contenidor
- [x] Munta un volum del sistema amfitrió al contenidor
- [ ] Especifica la versió del contenidor

Més informació: [1.3.4. Munteu dades al contenidor](#134-mount-data-into-the-container)
</quiz>

<quiz>
Per què necessiteu muntar volums quan utilitzeu contenidors?
- [ ] Per millorar el rendiment del contenidor
- [ ] Per estalviar espai de disc
- [x] Perquè els contenidors estan aïllats del sistema de fitxers de l'amfitrió per defecte
- [ ] Per habilitar la xarxa

Més informació: [1.3.4. Munteu dades al contenidor](#134-mount-data-into-the-container)
</quiz>

<quiz>
Com especifiqueu un contenidor per a un procés de Nextflow?
- [ ] `docker 'container-uri'`
- [ ] `image 'container-uri'`
- [x] `container 'container-uri'`
- [ ] `use 'container-uri'`

Més informació: [2.3.1. Especifiqueu un contenidor per a cowpy](#231-specify-a-container-for-cowpy)
</quiz>

<quiz>
Quina configuració de `nextflow.config` habilita Docker per al vostre workflow?
- [ ] `#!groovy process.docker = true`
- [x] `#!groovy docker.enabled = true`
- [ ] `#!groovy container.engine = 'docker'`
- [ ] `#!groovy docker.activate = true`

Més informació: [2.3.2. Habiliteu l'ús de Docker mitjançant el fitxer `nextflow.config`](#232-enable-use-of-docker-via-the-nextflowconfig-file)
</quiz>

<quiz>
Què gestiona automàticament Nextflow quan executa un procés en un contenidor? (Seleccioneu totes les que corresponguin)
- [x] Descarregar la imatge del contenidor si cal
- [x] Muntar el directori de treball
- [x] Executar l'script del procés dins del contenidor
- [x] Netejar la instància del contenidor després de l'execució

Més informació: [2.3.4. Inspeccioneu com Nextflow va llançar la tasca contenidoritzada](#234-inspect-how-nextflow-launched-the-containerized-task)
</quiz>
