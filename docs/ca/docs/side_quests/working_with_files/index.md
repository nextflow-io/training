# Processament de fitxers d'entrada

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Traducció assistida per IA - [més informació i suggeriments](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

Els workflows d'anàlisi científica sovint impliquen el processament d'un gran nombre de fitxers.
Nextflow proporciona eines potents per gestionar fitxers de manera eficient, ajudant-vos a organitzar i processar les vostres dades amb un mínim de codi.

### Objectius d'aprenentatge

En aquesta missió secundària, explorarem com Nextflow gestiona els fitxers, des d'operacions bàsiques fins a tècniques més avançades per treballar amb col·leccions de fitxers.
Aprendreu a extreure metadades dels noms de fitxer, que és un requisit habitual en els pipelines d'anàlisi científica.

Al final d'aquesta missió secundària, sereu capaços de:

- Crear objectes Path a partir de cadenes de rutes de fitxer utilitzant el mètode `file()` de Nextflow
- Accedir als atributs dels fitxers com ara el nom, l'extensió i el directori pare
- Gestionar fitxers locals i remots de manera transparent mitjançant URIs
- Utilitzar canals per automatitzar la gestió de fitxers amb `channel.fromPath()` i `channel.fromFilePairs()`
- Extreure i estructurar metadades dels noms de fitxer mitjançant la manipulació de cadenes
- Agrupar fitxers relacionats utilitzant coincidència de patrons i expressions glob
- Integrar operacions de fitxers en processos Nextflow amb una gestió correcta de l'entrada
- Organitzar les sortides dels processos mitjançant estructures de directoris basades en metadades

Aquestes habilitats us ajudaran a construir workflows que puguin gestionar diferents tipus d'entrades de fitxers amb gran flexibilitat.

### Prerequisits

Abans d'abordar aquesta missió secundària, hauríeu de:

- Haver completat el tutorial [Hello Nextflow](../../hello_nextflow/) o un curs equivalent per a principiants.
- Estar còmodes amb els conceptes i mecanismes bàsics de Nextflow (processos, canals, operadors)

<!-- He eliminat el suggeriment de fer primer la missió secundària de metamaps perquè funciona de manera més natural després -->

---

## 0. Primers passos

#### Obriu l'espai de treball de formació

Si encara no ho heu fet, assegureu-vos d'obrir l'entorn de formació tal com es descriu a la [Configuració de l'entorn](../envsetup/index.md).

[![Open in GitHub Codespaces](https://github.com/codespaces/badge.svg)](https://codespaces.new/nextflow-io/training?quickstart=1&ref=master)

#### Moveu-vos al directori del projecte

Anem al directori on es troben els fitxers d'aquest tutorial.

```bash
cd side-quests/working_with_files
```

Podeu configurar VSCode perquè es centri en aquest directori:

```bash
code .
```

#### Reviseu els materials

Trobareu un fitxer de workflow senzill anomenat `main.nf`, un directori `modules` que conté dos fitxers de mòdul, i un directori `data` que conté alguns fitxers de dades d'exemple.

??? abstract "Contingut del directori"

    ```console
    .
    ├── data
    │   ├── patientA_rep1_normal_R1_001.fastq.gz
    │   ├── patientA_rep1_normal_R2_001.fastq.gz
    │   ├── patientA_rep1_tumor_R1_001.fastq.gz
    │   ├── patientA_rep1_tumor_R2_001.fastq.gz
    │   ├── patientA_rep2_normal_R1_001.fastq.gz
    │   ├── patientA_rep2_normal_R2_001.fastq.gz
    │   ├── patientA_rep2_tumor_R1_001.fastq.gz
    │   ├── patientA_rep2_tumor_R2_001.fastq.gz
    │   ├── patientB_rep1_normal_R1_001.fastq.gz
    │   ├── patientB_rep1_normal_R2_001.fastq.gz
    │   ├── patientB_rep1_tumor_R1_001.fastq.gz
    │   ├── patientB_rep1_tumor_R2_001.fastq.gz
    │   ├── patientC_rep1_normal_R1_001.fastq.gz
    │   ├── patientC_rep1_normal_R2_001.fastq.gz
    │   ├── patientC_rep1_tumor_R1_001.fastq.gz
    │   └── patientC_rep1_tumor_R2_001.fastq.gz
    ├── main.nf
    └── modules
        ├── analyze_reads.nf
        └── count_lines.nf
    ```

Aquest directori conté dades de seqüenciació paired-end de tres pacients (A, B, C).

Per a cada pacient, tenim mostres de tipus `tumor` (que normalment provenen de biòpsies tumorals) o `normal` (extretes de teixit sa o sang).
Si no esteu familiaritzats amb l'anàlisi del càncer, simplement sapigueu que això correspon a un model experimental que utilitza mostres tumor/normal aparellades per realitzar anàlisis contrastives.

Específicament per al pacient A, tenim dos conjunts de rèpliques tècniques (repeticions).

Els fitxers de dades de seqüenciació s'anomenen amb la convenció típica `_R1_` i `_R2_` per al que es coneix com a 'lectures directes' i 'lectures inverses'.

_No us preocupeu si no esteu familiaritzats amb aquest disseny experimental, no és crític per entendre aquest tutorial._

#### Reviseu l'assignació

El vostre repte és escriure un workflow Nextflow que:

1. **Carregui** els fitxers d'entrada utilitzant els mètodes de gestió de fitxers de Nextflow
2. **Extregui** metadades (ID del pacient, rèplica, tipus de mostra) de l'estructura del nom de fitxer
3. **Agrupe** els fitxers aparellats (R1/R2) junts utilitzant `channel.fromFilePairs()`
4. **Processi** els fitxers amb un mòdul d'anàlisi proporcionat
5. **Organitzi** les sortides en una estructura de directoris basada en les metadades extretes

#### Llista de verificació de preparació

Creieu que esteu preparats per començar?

- [ ] Entenc l'objectiu d'aquest curs i els seus prerequisits
- [ ] El meu espai de treball està en funcionament
- [ ] He configurat el meu directori de treball correctament
- [ ] Entenc l'assignació

Si podeu marcar totes les caselles, esteu a punt per començar.

---

## 1. Operacions bàsiques amb fitxers

### 1.1. Identificar el tipus d'un objecte amb `.class`

Feu una ullada al fitxer de workflow `main.nf`:

```groovy title="main.nf" linenums="1"
#!/usr/bin/env nextflow

workflow {

    // Crea un objecte Path a partir d'una ruta de cadena
    myFile = 'data/patientA_rep1_normal_R1_001.fastq.gz'

    println "${myFile} is of class ${myFile.class}"
}
```

Aquest és un mini-workflow (sense cap procés) que fa referència a una única ruta de fitxer en el seu workflow, i després l'imprimeix a la consola, juntament amb la seva classe.

??? info "Què és `.class`?"

    A Nextflow, `.class` ens indica quin tipus d'objecte estem treballant. És com preguntar "quin tipus de cosa és aquesta?" per esbrinar si és una cadena, un número, un fitxer o alguna altra cosa.
    Això ens ajudarà a il·lustrar la diferència entre una cadena simple i un objecte Path a les seccions següents.

Executem el workflow:

```bash
nextflow run main.nf
```

??? success "Sortida de la comanda"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `main.nf` [romantic_chandrasekhar] DSL2 - revision: 5a4a89bc3a

    data/patientA_rep1_normal_R1_001.fastq.gz is of class java.lang.String
    ```

Com podeu veure, Nextflow ha imprès la ruta de cadena exactament tal com l'hem escrit.

Això és simplement sortida de text; Nextflow encara no ha fet res especial amb ella.
També hem confirmat que, pel que fa a Nextflow, això és només una cadena (de classe `java.lang.String`).
Té sentit, ja que encara no li hem dit a Nextflow que correspon a un fitxer.

### 1.2. Crear un objecte Path amb file()

Podem dir a Nextflow com gestionar els fitxers creant [objectes Path](https://www.nextflow.io/docs/latest/reference/stdlib-types.html#path) a partir de cadenes de rutes.

En el nostre workflow, podem convertir la cadena de ruta `data/patientA_rep1_normal_R1_001.fastq.gz` en un objecte Path utilitzant el mètode `file()`, que proporciona accés a les propietats i operacions del fitxer.

Editeu el `main.nf` per embolcallar la cadena amb `file()` de la manera següent:

=== "Després"

    ```groovy title="main.nf" linenums="5" hl_lines="2"
        // Crea un objecte Path a partir d'una ruta de cadena
        myFile = file('data/patientA_rep1_normal_R1_001.fastq.gz')

        println "${myFile} is of class ${myFile.class}"
    ```

=== "Abans"

    ```groovy title="main.nf" linenums="5" hl_lines="2"
        // Crea un objecte Path a partir d'una ruta de cadena
        myFile = 'data/patientA_rep1_normal_R1_001.fastq.gz'

        println "${myFile} is of class ${myFile.class}"
    ```

Ara executeu el workflow de nou:

```bash
nextflow run main.nf
```

??? success "Sortida de la comanda"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `main.nf` [kickass_coulomb] DSL2 - revision: 5af44b1b59

    /workspaces/training/side-quests/working_with_files/data/patientA_rep1_normal_R1_001.fastq.gz is of class class sun.nio.fs.UnixPath
    ```

Aquesta vegada, veieu la ruta absoluta completa en lloc de la ruta relativa que hem proporcionat com a entrada.

Nextflow ha convertit la nostra cadena en un objecte Path i l'ha resolt a la ubicació real del fitxer al sistema.
La ruta del fitxer ara serà absoluta, com ara `/workspaces/training/side-quests/working_with_files/data/patientA_rep1_normal_R1_001.fastq.gz`.

Observeu també que la classe de l'objecte Path és `sun.nio.fs.UnixPath`: aquesta és la manera de Nextflow de representar fitxers locals.
Com veurem més endavant, els fitxers remots tindran noms de classe diferents (com ara `nextflow.file.http.XPath` per a fitxers HTTP), però tots funcionen exactament de la mateixa manera i es poden utilitzar de manera idèntica en els vostres workflows.

!!! tip "Consell"

    **La diferència clau:**

    - **Cadena de ruta**: Simplement text que Nextflow tracta com a caràcters
    - **Objecte Path**: Una referència de fitxer intel·ligent amb la qual Nextflow pot treballar

    Penseu-hi d'aquesta manera: una cadena de ruta és com escriure una adreça en paper, mentre que un objecte Path és com tenir l'adreça carregada en un dispositiu GPS que sap com navegar fins allà i us pot donar detalls sobre el trajecte.

### 1.3. Accedir als atributs dels fitxers

Per què és útil això? Doncs bé, ara que Nextflow entén que `myFile` és un objecte Path i no simplement una cadena, podem accedir als diversos atributs de l'objecte Path.

Actualitzem el nostre workflow per imprimir els atributs de fitxer integrats:

=== "Després"

    ```groovy title="main.nf" linenums="5" hl_lines="4-9"
        // Crea un objecte Path a partir d'una ruta de cadena
        myFile = file('data/patientA_rep1_normal_R1_001.fastq.gz')

        // Imprimeix els atributs del fitxer
        println "File object class: ${myFile.class}"
        println "File name: ${myFile.name}"
        println "Simple name: ${myFile.simpleName}"
        println "Extension: ${myFile.extension}"
        println "Parent directory: ${myFile.parent}"
    ```

=== "Abans"

    ```groovy title="main.nf" linenums="5" hl_lines="4"
        // Crea un objecte Path a partir d'una ruta de cadena
        myFile = file('data/patientA_rep1_normal_R1_001.fastq.gz')

        println "${myFile} is of class ${myFile.class}"
    ```

Executeu el workflow:

```bash
nextflow run main.nf
```

??? success "Sortida de la comanda"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `main.nf` [ecstatic_ampere] DSL2 - revision: f3fa3dcb48

    File object class: sun.nio.fs.UnixPath
    File name: patientA_rep1_normal_R1_001.fastq.gz
    Simple name: patientA_rep1_normal_R1_001
    Extension: gz
    Parent directory: /workspaces/training/side-quests/working_with_files/data
    ```

Podeu veure els diversos atributs del fitxer impresos a la consola.

### 1.4. Passar el fitxer a un procés

La diferència entre cadenes i objectes Path es torna crítica quan comenceu a construir workflows reals amb processos.
Fins ara hem verificat que Nextflow ara tracta el nostre fitxer d'entrada com a fitxer, però vegem si realment podem executar alguna cosa sobre aquest fitxer en un procés.

#### 1.4.1. Importar el procés i examinar el codi

Us proporcionem un mòdul de procés pre-escrit anomenat `COUNT_LINES` que pren un fitxer d'entrada i compta quantes línies conté.

Per utilitzar el procés en el workflow, simplement heu d'afegir una declaració include abans del bloc workflow:

=== "Després"

    ```groovy title="main.nf" linenums="1" hl_lines="3"
    #!/usr/bin/env nextflow

    include { COUNT_LINES } from './modules/count_lines.nf'

    workflow {
    ```

=== "Abans"

    ```groovy title="main.nf" linenums="1"
    #!/usr/bin/env nextflow

    workflow {
    ```

Podeu obrir el fitxer del mòdul per examinar el seu codi:

```groovy title="modules/count_lines.nf" linenums="1"
#!/usr/bin/env nextflow

process COUNT_LINES {
    debug true

    input:
    path input_file

    script:
    """
    set -o pipefail
    echo "Processing file: $input_file"
    gzip -dc $input_file | wc -l
    """
}
```

Com podeu veure, és un script bastant senzill que descomprimeix el fitxer i compta quantes línies conté.

??? info "Què fa `debug true`?"

    La directiva `debug true` en la definició del procés fa que Nextflow imprimeixi la sortida del vostre script (com el recompte de línies "40") directament al registre d'execució.
    Sense això, només veuríeu l'estat d'execució del procés però no la sortida real del vostre script.

    Per obtenir més informació sobre la depuració de processos Nextflow, consulteu la missió secundària [Debugging Nextflow Workflows](debugging.md).

#### 1.4.2. Afegir una crida a `COUNT_LINES`

Ara que el procés està disponible per al workflow, podem afegir una crida al procés `COUNT_LINES` per executar-lo sobre el fitxer d'entrada.

Feu les edicions següents al workflow:

=== "Després"

    ```groovy title="main.nf" linenums="7" hl_lines="11-12"
        // Crea un objecte Path a partir d'una ruta de cadena
        myFile = file('data/patientA_rep1_normal_R1_001.fastq.gz')

        // Imprimeix els atributs del fitxer
        println "File object class: ${myFile.class}"
        println "File name: ${myFile.name}"
        println "Simple name: ${myFile.simpleName}"
        println "Extension: ${myFile.extension}"
        println "Parent directory: ${myFile.parent}"

        // Compta les línies del fitxer
        COUNT_LINES(myFile)
    ```

=== "Abans"

    ```groovy title="main.nf" linenums="7" hl_lines="4-9"
        // Crea un objecte Path a partir d'una ruta de cadena
        myFile = file('data/patientA_rep1_normal_R1_001.fastq.gz')

        // Imprimeix els atributs del fitxer
        println "File object class: ${myFile.class}"
        println "File name: ${myFile.name}"
        println "Simple name: ${myFile.simpleName}"
        println "Extension: ${myFile.extension}"
        println "Parent directory: ${myFile.parent}"
    ```

I ara executeu el workflow:

```bash
nextflow run main.nf
```

??? success "Sortida de la comanda"

    ```console
     N E X T F L O W   ~  version 25.10.2

    Launching `main.nf` [cheeky_hypatia] DSL2 - revision: 281d13c414

    File object class: class sun.nio.fs.UnixPath
    File name: patientA_rep1_normal_R1_001.fastq.gz
    Simple name: patientA_rep1_normal_R1_001
    Extension: gz
    Parent directory: /workspaces/training/side-quests/working_with_files/data
    executor >  local (1)
    [e9/341c05] COUNT_LINES [100%] 1 of 1 ✔
    Processing file: /workspaces/training/side-quests/working_with_files/data/patientA_rep1_normal_R1_001.fastq.gz
    40
    ```

Això mostra que som capaços d'operar sobre el fitxer correctament dins d'un procés.

Concretament, Nextflow ha dut a terme les operacions següents amb èxit:

- Ha preparat (staged) el fitxer al directori de treball
- Ha descomprimit el fitxer .gz
- Ha comptat les línies (40 línies en aquest cas)
- Ha completat sense errors

La clau d'aquesta operació fluida és que estem dient explícitament a Nextflow que la nostra entrada és un fitxer i s'ha de tractar com a tal.

### 1.5. Solucionar errors bàsics d'entrada de fitxers

Això sovint confon els nouvinguts a Nextflow, així que dediquem uns minuts a veure què passa quan ho feu malament.

Hi ha dos llocs principals on podeu equivocar-vos en la gestió de fitxers: al nivell del workflow i al nivell del procés.

#### 1.5.1. Error a nivell de workflow

Vegem què passa si tornem a tractar el fitxer com una cadena quan especifiquem l'entrada en el bloc workflow.

Feu les edicions següents al workflow, assegurant-vos de comentar les declaracions d'impressió específiques de ruta:

=== "Després"

    ```groovy title="main.nf" linenums="7" hl_lines="2 6-11"
        // Crea un objecte Path a partir d'una ruta de cadena
        myFile = 'data/patientA_rep1_normal_R1_001.fastq.gz'

        // Imprimeix els atributs del fitxer
        println "File object class: ${myFile.class}"
        /*
        println "File name: ${myFile.name}"
        println "Simple name: ${myFile.simpleName}"
        println "Extension: ${myFile.extension}"
        println "Parent directory: ${myFile.parent}"
        */

        // Compta les línies del fitxer
        COUNT_LINES(myFile)
    ```

=== "Abans"

    ```groovy title="main.nf" linenums="7" hl_lines="4-9"
        // Crea un objecte Path a partir d'una ruta de cadena
        myFile = file('data/patientA_rep1_normal_R1_001.fastq.gz')

        // Imprimeix els atributs del fitxer
        println "File object class: ${myFile.class}"
        println "File name: ${myFile.name}"
        println "Simple name: ${myFile.simpleName}"
        println "Extension: ${myFile.extension}"
        println "Parent directory: ${myFile.parent}"

        // Compta les línies del fitxer
        COUNT_LINES(myFile)
    ```

I ara executeu el workflow:

```bash
nextflow run main.nf
```

??? failure "Sortida de la comanda"

    ```console
     N E X T F L O W   ~  version 25.10.2

    Launching `main.nf` [friendly_goodall] DSL2 - revision: ae50609b20

    [-        ] COUNT_LINES -
    ERROR ~ Error executing process > 'COUNT_LINES'

    Caused by:
      Not a valid path value: 'data/patientA_rep1_normal_R1_001.fastq.gz'



    Tip: view the complete command output by changing to the process work dir and entering the command `cat .command.out`

    -- Check '.nextflow.log' file for details
    ```

Aquesta és la part important:

```console
Not a valid path value: 'data/patientA_rep1_normal_R1_001.fastq.gz'
```

Quan especifiqueu una entrada `path`, Nextflow valida que esteu passant referències de fitxers reals, no simplement cadenes.
Aquest error us indica que `'data/patientA_rep1_normal_R1_001.fastq.gz'` no és un valor de ruta vàlid perquè és una cadena, no un objecte Path.

Nextflow ha detectat immediatament el problema i s'ha aturat abans fins i tot d'iniciar el procés.

#### 1.5.2. Error a nivell de procés

L'altre lloc on podríem oblidar especificar que volem que Nextflow tracti l'entrada com a fitxer és en la definició del procés.

!!! warning "Manteniu l'error de workflow de 1.5.1"

    Perquè aquesta prova funcioni correctament, manteniu el workflow en el seu estat trencat (utilitzant una cadena simple en lloc de `file()`).
    Quan es combina amb `val` en el procés, això produeix l'error que es mostra a continuació.

Feu la edició següent al mòdul:

=== "Després"

    ```groovy title="modules/count_lines.nf" linenums="3" hl_lines="5"
    process COUNT_LINES {
        debug true

        input:
        val input_file
    ```

=== "Abans"

    ```groovy title="modules/count_lines.nf" linenums="3" hl_lines="5"
    process COUNT_LINES {
        debug true

        input:
        path input_file
    ```

I ara executeu el workflow de nou:

```bash
nextflow run main.nf
```

??? failure "Sortida de la comanda"

    ```console
     N E X T F L O W   ~  version 25.10.2

    Launching `main.nf` [soggy_golick] DSL2 - revision: ae50609b20

    executor >  local (1)
    [b3/b3023c] COUNT_LINES [  0%] 0 of 1 ✘
    ERROR ~ Error executing process > 'COUNT_LINES'

    Caused by:
      Process `COUNT_LINES` terminated with an error exit status (1)


    Command executed:

      set -o pipefail
      echo "Processing file: data/patientA_rep1_normal_R1_001.fastq.gz"
      gzip -dc data/patientA_rep1_normal_R1_001.fastq.gz | wc -l

    Command exit status:
      1

    Command output:
      Processing file: data/patientA_rep1_normal_R1_001.fastq.gz
      0

    Command error:
      Processing file: data/patientA_rep1_normal_R1_001.fastq.gz
      gzip: data/patientA_rep1_normal_R1_001.fastq.gz: No such file or directory
      0

    Work dir:
      /workspaces/training/side-quests/working_with_files/work/b3/b3023cb2ccb986851301d8e369e79f

    Tip: when you have fixed the problem you can continue the execution adding the option `-resume` to the run command line

    -- Check '.nextflow.log' file for details
    ```

Això mostra molts detalls sobre l'error perquè el procés està configurat per mostrar informació de depuració, tal com s'ha indicat anteriorment.

Aquestes són les seccions més rellevants:

```console
Command executed:

  set -o pipefail
  echo "Processing file: data/patientA_rep1_normal_R1_001.fastq.gz"
  gzip -dc data/patientA_rep1_normal_R1_001.fastq.gz | wc -l
```

```console
Command error:
  Processing file: data/patientA_rep1_normal_R1_001.fastq.gz
  gzip: data/patientA_rep1_normal_R1_001.fastq.gz: No such file or directory
  0
```

Això indica que el sistema no ha pogut trobar el fitxer; però si busqueu la ruta, hi ha un fitxer amb aquest nom en aquella ubicació.

Quan hem executat això, Nextflow ha passat el valor de cadena a l'script, però no ha _preparat_ (staged) el fitxer real al directori de treball.
Així que el procés ha intentat utilitzar la cadena relativa, `data/patientA_rep1_normal_R1_001.fastq.gz`, però aquest fitxer no existeix dins del directori de treball del procés.

En conjunt, aquests dos exemples us mostren la importància de dir a Nextflow si una entrada s'ha de gestionar com a fitxer.

!!! note "Nota"

    Assegureu-vos de tornar enrere i corregir tots dos errors intencionals abans de continuar a la secció següent.

### Conclusió

- Cadenes de ruta vs objectes Path: Les cadenes són simplement text, els objectes Path són referències de fitxer intel·ligents
- El mètode `file()` converteix una cadena de ruta en un objecte Path amb el qual Nextflow pot treballar
- Podeu accedir a les propietats dels fitxers com ara `name`, `simpleName`, `extension` i `parent` [utilitzant els atributs de fitxer](https://www.nextflow.io/docs/latest/working-with-files.html#getting-file-attributes)
- L'ús d'objectes Path en lloc de cadenes permet a Nextflow gestionar correctament els fitxers en el vostre workflow
- Resultats de l'entrada del procés: La gestió correcta de fitxers requereix objectes Path, no cadenes, per garantir que els fitxers es preparin (staged) i siguin accessibles per als processos.

---

## 2. Ús de fitxers remots

Una de les característiques clau de Nextflow és la capacitat de canviar de manera transparent entre fitxers locals (a la mateixa màquina) i fitxers remots accessibles per internet.

Si ho feu correctament, mai hauríeu de canviar la lògica del vostre workflow per adaptar-vos a fitxers que provenen de diferents ubicacions.
Tot el que heu de fer per utilitzar un fitxer remot és especificar el prefix adequat a la ruta del fitxer quan el proporcioneu al workflow.

Per exemple, `/path/to/data` no té cap prefix, indicant que és una ruta de fitxer local 'normal', mentre que `s3://path/to/data` inclou el prefix `s3://`, indicant que es troba a l'emmagatzematge d'objectes S3 d'Amazon.

Molts protocols diferents estan suportats:

- HTTP(S)/FTP (http://, https://, ftp://)
- Amazon S3 (s3://)
- Azure Blob Storage (az://)
- Google Cloud Storage (gs://)

Per utilitzar qualsevol d'aquests, simplement especifiqueu el prefix rellevant a la cadena, que tècnicament s'anomena Identificador Uniforme de Recursos (URI) en lloc d'una ruta de fitxer.
Nextflow gestionarà l'autenticació i la preparació (staging) dels fitxers al lloc correcte, descarregant o pujant i totes les altres operacions de fitxers que esperaríeu.

La fortalesa clau d'aquest sistema és que ens permet canviar entre entorns sense canviar cap lògica del pipeline.
Per exemple, podeu desenvolupar amb un conjunt de proves petit i local abans de canviar a un conjunt de proves a escala completa ubicat en emmagatzematge remot simplement canviant l'URI.

### 2.1. Utilitzar un fitxer d'internet

Provem-ho canviant la ruta local que proporcionem al nostre workflow per una ruta HTTPS que apunta a una còpia de les mateixes dades emmagatzemades a Github.

!!! warning "Advertència"

    Això només funcionarà si teniu una connexió a internet activa.

Obriu `main.nf` de nou i canvieu la ruta d'entrada de la manera següent:

=== "Després"

    ```groovy title="main.nf" linenums="2" hl_lines="2"
        // Utilitza un fitxer remot d'internet
        myFile = file('https://raw.github.com/nextflow-io/training/master/side-quests/working_with_files/data/patientA_rep1_normal_R1_001.fastq.gz')

        // Imprimeix els atributs del fitxer
        println "File object class: ${myFile.class}"
        println "File name: ${myFile.name}"
        println "Simple name: ${myFile.simpleName}"
        println "Extension: ${myFile.extension}"
        println "Parent directory: ${myFile.parent}"
    ```

=== "Abans"

    ```groovy title="main.nf" linenums="2" hl_lines="2"
        // Crea un objecte Path a partir d'una ruta de cadena
        myFile = file('data/patientA_rep1_normal_R1_001.fastq.gz')

        // Imprimeix els atributs del fitxer
        println "File object class: ${myFile.class}"
        println "File name: ${myFile.name}"
        println "Simple name: ${myFile.simpleName}"
        println "Extension: ${myFile.extension}"
        println "Parent directory: ${myFile.parent}"
    ```

Executem el workflow:

```bash
nextflow run main.nf
```

??? success "Sortida de la comanda"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `main.nf` [insane_swartz] DSL2 - revision: fff18abe6d

    File object class: class nextflow.file.http.XPath
    File name: patientA_rep1_normal_R1_001.fastq.gz
    Simple name: patientA_rep1_normal_R1_001
    Extension: gz
    Parent directory: /nextflow-io/training/master/side-quests/working_with_files/data
    executor >  local (1)
    [8a/2ab7ca] COUNT_LINES [100%] 1 of 1 ✔
    Processing file: patientA_rep1_normal_R1_001.fastq.gz
    40
    ```

Funciona! Podeu veure que molt poc ha canviat.

L'única diferència a la sortida de la consola és que la classe de l'objecte Path ara és `nextflow.file.http.XPath`, mentre que per a la ruta local la classe era `sun.nio.fs.UnixPath`.
No cal que recordeu aquestes classes; simplement les mencionem per demostrar que Nextflow identifica i gestiona les diferents ubicacions de manera adequada.

Entre bastidors, Nextflow ha descarregat el fitxer a un directori de preparació (staging) ubicat dins del directori de treball.
Aquest fitxer preparat es pot tractar com un fitxer local i es pot enllaçar simbòlicament al directori del procés rellevant.

Podeu verificar que això ha passat aquí mirant el contingut del directori de treball ubicat al valor hash del procés.

??? abstract "Contingut del directori work"

    Si el hash del procés era `8a/2ab7ca`, podríeu explorar el directori de treball:

    ```console
    $ ls -la work/8a/2ab7ca*/
    total 16
    drwxr-xr-x  6 user  staff   192 Jan 28 10:00 .
    drwxr-xr-x  3 user  staff    96 Jan 28 10:00 ..
    -rw-r--r--  1 user  staff     0 Jan 28 10:00 .command.begin
    -rw-r--r--  1 user  staff   127 Jan 28 10:00 .command.sh
    lrwxr-xr-x  1 user  staff    89 Jan 28 10:00 patientA_rep1_normal_R1_001.fastq.gz -> /path/to/work/stage/.../patientA_rep1_normal_R1_001.fastq.gz
    ```

    L'enllaç simbòlic apunta a una còpia preparada del fitxer remot que Nextflow ha descarregat automàticament.

Tingueu en compte que per a fitxers més grans, el pas de descàrrega trigarà un temps addicional en comparació amb l'execució sobre fitxers locals.
No obstant això, Nextflow comprova si ja té una còpia preparada per evitar descàrregues innecessàries.
Així que si torneu a executar sobre el mateix fitxer i no heu eliminat el fitxer preparat, Nextflow utilitzarà la còpia preparada.

Això mostra com és de fàcil canviar entre dades locals i remotes utilitzant Nextflow, que és una característica clau de Nextflow.

!!! note "Nota"

    L'única excepció important a aquest principi és que no podeu utilitzar patrons glob o rutes de directori amb HTTPS perquè HTTPS no pot llistar múltiples fitxers, de manera que heu d'especificar URLs de fitxer exactes.
    No obstant això, altres protocols d'emmagatzematge com l'emmagatzematge blob (`s3://`, `az://`, `gs://`) poden utilitzar tant globs com rutes de directori.

    Aquí teniu com podríeu utilitzar patrons glob amb emmagatzematge al núvol:

    ```groovy title="Cloud storage examples (not runnable in this environment)"
    // S3 amb patrons glob - coincidiria amb múltiples fitxers
    ch_s3_files = channel.fromPath('s3://my-bucket/data/*.fastq.gz')

    // Azure Blob Storage amb patrons glob
    ch_azure_files = channel.fromPath('az://container/data/patient*_R{1,2}.fastq.gz')

    // Google Cloud Storage amb patrons glob
    ch_gcs_files = channel.fromPath('gs://bucket/data/sample_*.fastq.gz')
    ```

    Us mostrarem com treballar amb globs a la pràctica a la secció següent.

### 2.2. Tornar al fitxer local

Tornarem a utilitzar els nostres fitxers d'exemple locals per a la resta d'aquesta missió secundària, així que tornem a canviar l'entrada del workflow al fitxer original:

=== "Després"

    ```groovy title="main.nf" linenums="2" hl_lines="2"
        // Crea un objecte Path a partir d'una ruta de cadena
        myFile = file('data/patientA_rep1_normal_R1_001.fastq.gz')

        // Imprimeix els atributs del fitxer
        println "File object class: ${myFile.class}"
        println "File name: ${myFile.name}"
        println "Simple name: ${myFile.simpleName}"
        println "Extension: ${myFile.extension}"
        println "Parent directory: ${myFile.parent}"
    ```

=== "Abans"

    ```groovy title="main.nf" linenums="2" hl_lines="2"
        // Crea un objecte Path a partir d'una ruta de cadena
        myFile = file('https://raw.github.com/nextflow-io/training/master/side-quests/working_with_files/data/patientA_rep1_normal_R1_001.fastq.gz')

        // Imprimeix els atributs del fitxer
        println "File object class: ${myFile.class}"
        println "File name: ${myFile.name}"
        println "Simple name: ${myFile.simpleName}"
        println "Extension: ${myFile.extension}"
        println "Parent directory: ${myFile.parent}"
    ```

### Conclusió

- Les dades remotes s'accedeixen mitjançant un URI (HTTP, FTP, S3, Azure, Google Cloud)
- Nextflow descarregarà i prepararà (staged) automàticament les dades al lloc correcte, sempre que aquestes rutes s'estiguin passant a processos
- No escriviu lògica per descarregar o pujar fitxers remots!
- Els fitxers locals i remots produeixen tipus d'objectes diferents però funcionen de manera idèntica
- **Important**: HTTP/HTTPS només funciona amb fitxers individuals (sense patrons glob)
- L'emmagatzematge al núvol (S3, Azure, GCS) admet tant fitxers individuals com patrons glob
- Podeu canviar de manera transparent entre fonts de dades locals i remotes sense canviar la lògica del codi (sempre que el protocol admeti les operacions requerides)

---

## 3. Ús de la fàbrica de canals `fromPath()`

Fins ara hem estat treballant amb un sol fitxer a la vegada, però a Nextflow, normalment voldrem crear un canal d'entrada amb múltiples fitxers d'entrada per processar.

Una manera ingènua de fer-ho seria combinar el mètode `file()` amb [`channel.of()`](https://www.nextflow.io/docs/latest/reference/channel.html#of) d'aquesta manera:

```groovy title="Syntax example"
ch_files = channel.of([file('data/patientA_rep1_normal_R1_001.fastq.gz')],
                      [file('data/patientA_rep1_normal_R1_001.fastq.gz')])
```

Això funciona, però és poc elegant.

!!! tip "Quan utilitzar `file()` vs `channel.fromPath()`"

    - Utilitzeu `file()` quan necessiteu un únic objecte Path per a la manipulació directa (comprovar si un fitxer existeix, llegir els seus atributs o passar-lo a una única invocació de procés)
    - Utilitzeu `channel.fromPath()` quan necessiteu un canal que pugui contenir múltiples fitxers, especialment amb patrons glob, o quan els fitxers fluiran a través de múltiples processos

Aquí és on entra [`channel.fromPath()`](https://www.nextflow.io/docs/latest/reference/channel.html#frompath): una fàbrica de canals convenient que agrupa tota la funcionalitat que necessitem per generar un canal a partir d'una o més cadenes de fitxers estàtiques, així com patrons glob.

### 3.1. Afegir la fàbrica de canals

Actualitzem el nostre workflow per utilitzar `channel.fromPath`.

=== "Després"

    ```groovy title="main.nf" linenums="7" hl_lines="1-3"
        // Carrega fitxers amb channel.fromPath
        ch_files = channel.fromPath('data/patientA_rep1_normal_R1_001.fastq.gz')
        ch_files.view { myFile -> "Found file: $myFile" }

        // Imprimeix els atributs del fitxer
        /* Comentem-los per ara, hi tornarem!
        println "File object class: ${myFile.class}"
        println "File name: ${myFile.name}"
        println "Simple name: ${myFile.simpleName}"
        println "Extension: ${myFile.extension}"
        println "Parent directory: ${myFile.parent}"
        */

        // Compta les línies del fitxer
        // COUNT_LINES(myFile)
    ```

=== "Abans"

    ```groovy title="main.nf" linenums="7" hl_lines="1-2"
        // Crea un objecte Path a partir d'una ruta de cadena
        myFile = file('data/patientA_rep1_normal_R1_001.fastq.gz')

        // Imprimeix els atributs del fitxer
        println "File object class: ${myFile.class}"
        println "File name: ${myFile.name}"
        println "Simple name: ${myFile.simpleName}"
        println "Extension: ${myFile.extension}"
        println "Parent directory: ${myFile.parent}"

        // Compta les línies del fitxer
        COUNT_LINES(myFile)
    ```

També hem comentat el codi que imprimeix els atributs per ara, i hem afegit una declaració `.view` per imprimir només el nom del fitxer.

Executeu el workflow:

```bash
nextflow run main.nf
```

??? success "Sortida de la comanda"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `main.nf` [grave_meucci] DSL2 - revision: b09964a583

    Found file: /workspaces/training/side-quests/working_with_files/data/patientA_rep1_normal_R1_001.fastq.gz
    ```

Com podeu veure, la ruta del fitxer s'està carregant com un objecte de tipus `Path` al canal.
Això és similar al que hauria fet `file()`, excepte que ara tenim un canal al qual podem afegir més fitxers si volem.

L'ús de `channel.fromPath()` és una manera convenient de crear un nou canal poblat per una llista de fitxers.

### 3.2. Veure els atributs dels fitxers al canal

En el nostre primer intent d'utilitzar la fàbrica de canals, hem simplificat el codi i simplement hem imprès el nom del fitxer.

Tornem a imprimir els atributs complets del fitxer:

=== "Després"

    ```groovy title="main.nf" linenums="7" hl_lines="3-9 12"
        // Carrega fitxers amb channel.fromPath
        ch_files = channel.fromPath('data/patientA_rep1_normal_R1_001.fastq.gz')
        ch_files.view { myFile ->
            println "File object class: ${myFile.class}"
            println "File name: ${myFile.name}"
            println "Simple name: ${myFile.simpleName}"
            println "Extension: ${myFile.extension}"
            println "Parent directory: ${myFile.parent}"
        }

        // Compta les línies del fitxer
        COUNT_LINES(ch_files)
    ```

=== "Abans"

    ```groovy title="main.nf" linenums="7" hl_lines="3"
        // Carrega fitxers amb channel.fromPath
        ch_files = channel.fromPath('data/patientA_rep1_normal_R1_001.fastq.gz')
        ch_files.view { myFile -> "Found file: $myFile" }

        // Compta les línies del fitxer
        // COUNT_LINES(ch_files)
    ```

També estem reactivant la crida al procés `COUNT_LINES` per verificar que el processament de fitxers encara funciona correctament amb el nostre enfocament basat en canals.

Executeu el workflow:

```bash
nextflow run main.nf
```

??? success "Sortida de la comanda"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `main.nf` [furious_swanson] DSL2 - revision: c35c34950d

    executor >  local (1)
    [9d/6701a6] COUNT_LINES (1) [100%] 1 of 1 ✔
    File object class: sun.nio.fs.UnixPath
    File name: patientA_rep1_normal_R1_001.fastq.gz
    Simple name: patientA_rep1_normal_R1_001
    Extension: gz
    Parent directory: /workspaces/training/side-quests/working_with_files/data
    Processing file: patientA_rep1_normal_R1_001.fastq.gz
    40
    ```

I aquí ho teniu, els mateixos resultats que abans però ara tenim el fitxer en un canal, de manera que podem afegir-ne més.

### 3.3. Utilitzar un glob per coincidir amb múltiples fitxers

Hi ha diverses maneres de carregar més fitxers al canal.
Aquí us mostrarem com utilitzar patrons glob, que són una manera convenient de coincidir i recuperar noms de fitxers i directoris basant-se en caràcters comodí.
El procés de coincidència d'aquests patrons s'anomena "globbing" o "expansió de noms de fitxer".

!!! note "Nota"

    Com s'ha indicat anteriorment, Nextflow admet globbing per gestionar fitxers d'entrada i sortida en la majoria dels casos, excepte amb rutes de fitxer HTTPS perquè HTTPS no pot llistar múltiples fitxers.

Suposem que volem recuperar tots dos fitxers d'un parell de fitxers associats a un pacient determinat, `patientA`:

```console
patientA_rep1_normal_R1_001.fastq.gz
patientA_rep1_normal_R2_001.fastq.gz
```

Com que l'única diferència entre els noms de fitxer és el número de rèplica, _és a dir_, el número després de `R`, podem utilitzar el caràcter comodí `*` per substituir el número de la manera següent:

```console
patientA_rep1_normal_R*_001.fastq.gz
```

Aquest és el patró glob que necessitem.

Ara tot el que hem de fer és actualitzar la ruta del fitxer a la fàbrica de canals per utilitzar aquest patró glob de la manera següent:

=== "Després"

    ```groovy title="main.nf" linenums="7"
      // Carrega fitxers amb channel.fromPath
        ch_files = channel.fromPath('data/patientA_rep1_normal_R*_001.fastq.gz')
    ```

=== "Abans"

    ```groovy title="main.nf" linenums="7"
      // Carrega fitxers amb channel.fromPath
        ch_files = channel.fromPath('data/patientA_rep1_normal_R1_001.fastq.gz')
    ```

Nextflow reconeixerà automàticament que aquest és un patró glob i el gestionarà adequadament.

Executeu el workflow per provar-ho:

```bash
nextflow run main.nf
```

??? success "Sortida de la comanda"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `main.nf` [boring_sammet] DSL2 - revision: d2aa789c9a

    executor >  local (2)
    [3c/a65de5] COUNT_LINES (2) [100%] 2 of 2 ✔
    File object class: class sun.nio.fs.UnixPath
    File name: patientA_rep1_normal_R1_001.fastq.gz
    Simple name: patientA_rep1_normal_R1_001
    Extension: gz
    Parent directory: /workspaces/training/side-quests/working_with_files/data
    File object class: class sun.nio.fs.UnixPath
    File name: patientA_rep1_normal_R2_001.fastq.gz
    Simple name: patientA_rep1_normal_R2_001
    Extension: gz
    Parent directory: /workspaces/training/side-quests/working_with_files/data
    Processing file: patientA_rep1_normal_R1_001.fastq.gz
    40

    Processing file: patientA_rep1_normal_R2_001.fastq.gz
    40
    ```

Com podeu veure, ara tenim dos objectes Path al nostre canal, la qual cosa mostra que Nextflow ha fet l'expansió de noms de fitxer correctament, i ha carregat i processat tots dos fitxers com s'esperava.

Utilitzant aquest mètode, podem recuperar tants o tan pocs fitxers com vulguem simplement canviant el patró glob. Si el féssim més ampli, per exemple substituint totes les parts variables dels noms de fitxer per `*` (_p. ex._ `data/patient*_rep*_*_R*_001.fastq.gz`) podríem agafar tots els fitxers d'exemple del directori `data`.

### Conclusió

- `channel.fromPath()` crea un canal amb fitxers que coincideixen amb un patró
- Cada fitxer s'emet com un element separat al canal
- Podem utilitzar un patró glob per coincidir amb múltiples fitxers
- Els fitxers es converteixen automàticament en objectes Path amb tots els atributs
- El mètode `.view()` permet inspeccionar el contingut del canal

---

## 4. Extracció de metadades bàsiques dels noms de fitxer

En la majoria de dominis científics, és molt comú tenir metadades codificades en els noms dels fitxers que contenen les dades.
Per exemple, en bioinformàtica, els fitxers que contenen dades de seqüenciació sovint s'anomenen d'una manera que codifica informació sobre la mostra, la condició, la rèplica i el número de lectura.

Si els noms de fitxer es construeixen d'acord amb una convenció consistent, podeu extreure aquestes metadades d'una manera estandarditzada i utilitzar-les en el curs de la vostra anàlisi.
Això és un gran 'si', per descomptat, i hauríeu de ser molt cautelosos sempre que us baseu en l'estructura del nom de fitxer; però la realitat és que aquest enfocament s'utilitza molt àmpliament, així que vegem com es fa a Nextflow.

En el cas de les nostres dades d'exemple, sabem que els noms de fitxer inclouen metadades estructurades de manera consistent.
Per exemple, el nom de fitxer `patientA_rep1_normal_R2_001` codifica el següent:

- ID del pacient: `patientA`
- ID de la rèplica: `rep1`
- tipus de mostra: `normal` (en contraposició a `tumor`)
- conjunt de lectures: `R1` (en contraposició a `R2`)

Modificarem el nostre workflow per recuperar aquesta informació en tres passos:

1. Recuperar el `simpleName` del fitxer, que inclou les metadades
2. Separar les metadades utilitzant un mètode anomenat `tokenize()`
3. Utilitzar un map per organitzar les metadades

!!! warning "Advertència"

    Mai hauríeu de codificar informació sensible als noms de fitxer, com ara noms de pacients o altres característiques identificatives, ja que això pot comprometre la privacitat del pacient o altres restriccions de seguretat rellevants.

### 4.1. Recuperar el `simpleName`

El `simpleName` és un atribut de fitxer que correspon al nom del fitxer sense la seva ruta ni extensió.

Feu les edicions següents al workflow:

=== "Després"

    ```groovy title="main.nf" linenums="7" hl_lines="3-6"
        // Carrega fitxers amb channel.fromPath
        ch_files = channel.fromPath('data/patientA_rep1_normal_R*_001.fastq.gz')
        ch_files.map { myFile ->
            [ myFile.simpleName, myFile ]
        }
        .view()
    ```

=== "Abans"

    ```groovy title="main.nf" linenums="7" hl_lines="3-9"
        // Carrega fitxers amb channel.fromPath
        ch_files = channel.fromPath('data/patientA_rep1_normal_R*_001.fastq.gz')
        ch_files.view { myFile ->
            println "File object class: ${myFile.class}"
            println "File name: ${myFile.name}"
            println "Simple name: ${myFile.simpleName}"
            println "Extension: ${myFile.extension}"
            println "Parent directory: ${myFile.parent}"
        }
    ```

Això recupera el `simpleName` i l'associa amb l'objecte de fitxer complet utilitzant una operació `map()`.

Executeu el workflow per provar que funciona:

```bash
nextflow run main.nf
```

??? success "Sortida de la comanda"

    ```console hl_lines="7-8"
    N E X T F L O W   ~  version 25.10.2

    Launching `main.nf` [suspicious_mahavira] DSL2 - revision: ae8edc4e48

    executor >  local (2)
    [e9/55774b] COUNT_LINES (2) [100%] 2 of 2 ✔
    [patientA_rep1_normal_R2_001, /workspaces/training/side-quests/working_with_files/data/patientA_rep1_normal_R2_001.fastq.gz]
    [patientA_rep1_normal_R1_001, /workspaces/training/side-quests/working_with_files/data/patientA_rep1_normal_R1_001.fastq.gz]
    Processing file: patientA_rep1_normal_R1_001.fastq.gz
    40

    Processing file: patientA_rep1_normal_R2_001.fastq.gz
    40
    ```

Cada element del canal és ara una tupla que conté el `simpleName` i l'objecte de fitxer original.

### 4.2. Extreure les metadades del `simpleName`

En aquest punt, les metadades que volem estan incrustades al `simpleName`, però no podem accedir als elements individuals directament.
Així que hem de dividir el `simpleName` en els seus components.
Afortunadament, aquests components estan simplement separats per guions baixos al nom de fitxer original, de manera que podem aplicar un mètode comú de Nextflow anomenat `tokenize()` que és perfecte per a aquesta tasca.

Feu les edicions següents al workflow:

=== "Després"

    ```groovy title="main.nf" linenums="7" hl_lines="4"
        // Carrega fitxers amb channel.fromPath
        ch_files = channel.fromPath('data/patientA_rep1_normal_R*_001.fastq.gz')
        ch_files.map { myFile ->
            [ myFile.simpleName.tokenize('_'), myFile ]
        }
    ```

=== "Abans"

    ```groovy title="main.nf" linenums="7" hl_lines="4"
        // Carrega fitxers amb channel.fromPath
        ch_files = channel.fromPath('data/patientA_rep1_normal_R*_001.fastq.gz')
        ch_files.map { myFile ->
            [ myFile.simpleName, myFile ]
        }
    ```

El mètode `tokenize()` dividirà la cadena `simpleName` allà on trobi guions baixos, i retornarà una llista que conté les subcadenes.

Executeu el workflow:

```bash
nextflow run main.nf
```

??? success "Sortida de la comanda"

    ```console hl_lines="7-8"
    N E X T F L O W   ~  version 25.10.2

    Launching `main.nf` [gigantic_gauss] DSL2 - revision: a39baabb57

    executor >  local (2)
    [e7/da2f4b] COUNT_LINES (2) [100%] 2 of 2 ✔
    [[patientA, rep1, normal, R2, 001], /workspaces/training/side-quests/working_with_files/data/patientA_rep1_normal_R2_001.fastq.gz]
    [[patientA, rep1, normal, R1, 001], /workspaces/training/side-quests/working_with_files/data/patientA_rep1_normal_R1_001.fastq.gz]
    Processing file: patientA_rep1_normal_R2_001.fastq.gz
    40

    Processing file: patientA_rep1_normal_R1_001.fastq.gz
    40
    ```

Ara la tupla per a cada element del nostre canal conté la llista de metadades (_p. ex._ `[patientA, rep1, normal, R1, 001]`) i l'objecte de fitxer original.

Molt bé!
Hem desglossat la informació del nostre pacient d'una única cadena en una llista de cadenes.
Ara podem gestionar cada part de la informació del pacient per separat.

### 4.3. Utilitzar un map per organitzar les metadades

Les nostres metadades ara mateix són simplement una llista plana.
És prou fàcil d'utilitzar però difícil de llegir.

```console
[patientA, rep1, normal, R1, 001]
```

Quin és l'element a l'índex 3? Podeu dir-ho sense referir-vos a l'explicació original de l'estructura de metadades?

Aquesta és una gran oportunitat per utilitzar un magatzem clau-valor, on cada element té un conjunt de claus i els seus valors associats, de manera que podeu referir-vos fàcilment a cada clau per obtenir el valor corresponent.

En el nostre exemple, això significa passar d'aquesta organització:

```groovy
data = [patientA, 1, normal, R1]

println data[3]
```

A aquesta:

```groovy
data = [id: patientA, replicate: 1, type: normal, readNum: 1]

println data.readNum
```

A Nextflow, això s'anomena [map](https://nextflow.io/docs/latest/script.html#maps).

Convertim ara la nostra llista plana en un map.
Feu les edicions següents al workflow:

=== "Després"

    ```groovy title="main.nf" linenums="7" hl_lines="4-13"
        // Carrega fitxers amb channel.fromPath
        ch_files = channel.fromPath('data/patientA_rep1_normal_R*_001.fastq.gz')
        ch_files.map { myFile ->
            def (patient, replicate, type, readNum) = myFile.simpleName.tokenize('_')
            [
              [
                id: patient,
                replicate: replicate.replace('rep', ''),
                type: type,
                readNum: readNum.replace('R', ''),
              ],
              myFile
            ]
        }
    ```

=== "Abans"

    ```groovy title="main.nf" linenums="7" hl_lines="4"
        // Carrega fitxers amb channel.fromPath
        ch_files = channel.fromPath('data/patientA_rep1_normal_R*_001.fastq.gz')
        ch_files.map { myFile ->
            [ myFile.simpleName.tokenize('_'), myFile ]
        }
    ```

Els canvis clau aquí són:

- **Assignació per desestructuració**: `def (patient, replicate, type, readNum) = ...` extreu els valors tokenitzats en variables amb nom en una sola línia
- **Sintaxi literal de map**: `[id: patient, replicate: ...]` crea un map on cada clau (com `id`) s'associa amb un valor (com `patient`)
- **Estructura niada**: La llista exterior `[..., myFile]` aparella el map de metadades amb l'objecte de fitxer original

També hem simplificat un parell de les cadenes de metadades utilitzant un mètode de substitució de cadenes anomenat `replace()` per eliminar alguns caràcters que són innecessaris (_p. ex._ `replicate.replace('rep', '')` per conservar només el número dels IDs de rèplica).

Executem el workflow de nou:

```bash
nextflow run main.nf
```

??? success "Sortida de la comanda"

    ```console hl_lines="7-8"
    N E X T F L O W   ~  version 25.10.2

    Launching `main.nf` [infallible_swartz] DSL2 - revision: 7f4e68c0cb

    executor >  local (2)
    [1b/e7fb27] COUNT_LINES (1) [100%] 2 of 2 ✔
    [[id:patientA, replicate:1, type:normal, readNum:2], /workspaces/training/side-quests/working_with_files/data/patientA_rep1_normal_R2_001.fastq.gz]
    [[id:patientA, replicate:1, type:normal, readNum:1], /workspaces/training/side-quests/working_with_files/data/patientA_rep1_normal_R1_001.fastq.gz]
    Processing file: patientA_rep1_normal_R2_001.fastq.gz
    40

    Processing file: patientA_rep1_normal_R1_001.fastq.gz
    40
    ```

Ara les metadades estan clarament etiquetades (_p. ex._ `[id:patientA, replicate:1, type:normal, readNum:2]`) de manera que és molt més fàcil saber què és cada cosa.

També serà molt més fàcil fer ús real dels elements de metadades al workflow, i farà que el nostre codi sigui més fàcil de llegir i mantenir.

### Conclusió

- Podem gestionar els noms de fitxer a Nextflow amb el poder d'un llenguatge de programació complet
- Podem tractar els noms de fitxer com a cadenes per extreure informació rellevant
- L'ús de mètodes com `tokenize()` i `replace()` ens permet manipular cadenes al nom de fitxer
- L'operació `.map()` transforma els elements del canal preservant l'estructura
- Les metadades estructurades (maps) fan el codi més llegible i mantenible que les llistes posicionals

A continuació, veurem com gestionar fitxers de dades aparellats.

---

## 5. Gestió de fitxers de dades aparellats

Molts dissenys experimentals produeixen fitxers de dades aparellats que es beneficien de ser gestionats d'una manera explícitament aparellada.
Per exemple, en bioinformàtica, les dades de seqüenciació sovint es generen en forma de lectures aparellades, és a dir, cadenes de seqüència que provenen del mateix fragment d'ADN (sovint anomenades 'directes' i 'inverses' perquè es llegeixen des d'extrems oposats).

Aquest és el cas de les nostres dades d'exemple, on R1 i R2 fan referència als dos conjunts de lectures.

```console
data/patientA_rep1_normal_R1_001.fastq.gz
data/patientA_rep1_normal_R2_001.fastq.gz
```

Nextflow proporciona una fàbrica de canals especialitzada per treballar amb fitxers aparellats com aquest anomenada `channel.fromFilePairs()`, que agrupa automàticament els fitxers basant-se en un patró de nomenclatura compartit. Això us permet associar els fitxers aparellats de manera més estreta amb menys esforç.

Modificarem el nostre workflow per aprofitar-ho.
Caldrà fer-ho en dos passos:

1. Canviar la fàbrica de canals a `channel.fromFilePairs()`
2. Extreure i organitzar les metadades

### 5.1. Canviar la fàbrica de canals a `channel.fromFilePairs()`

Per utilitzar `channel.fromFilePairs`, hem d'especificar el patró que Nextflow hauria d'utilitzar per identificar els dos membres d'un parell.

Tornant a les nostres dades d'exemple, podem formalitzar el patró de nomenclatura de la manera següent:

```console
data/patientA_rep1_normal_R{1,2}_001.fastq.gz
```

Això és similar al patró glob que hem utilitzat anteriorment, excepte que aquest enumera específicament les subcadenes (ja sigui `1` o `2` just després de la R) que identifiquen els dos membres del parell.

Actualitzem el workflow `main.nf` en conseqüència:

=== "Després"

    ```groovy title="main.nf" linenums="7" hl_lines="1-2"
        // Carrega fitxers amb channel.fromFilePairs
        ch_files = channel.fromFilePairs('data/patientA_rep1_normal_R{1,2}_001.fastq.gz')
        /* Comentem el map per ara, hi tornarem!
        ch_files.map { myFile ->
            def (sample, replicate, type, readNum) = myFile.simpleName.tokenize('_')
            [
                [
                    id: sample,
                    replicate: replicate.replace('rep', ''),
                    type: type,
                    readNum: readNum,
                ],
                myFile
            ]
        }
        */
        .view()
    ```

=== "Abans"

    ```groovy title="main.nf" linenums="7" hl_lines="1-2"
        // Carrega fitxers amb channel.fromFilePairs
        ch_files = channel.fromPath('data/patientA_rep1_normal_R*_001.fastq.gz')
        ch_files.map { myFile ->
            def (sample, replicate, type, readNum) = myFile.simpleName.tokenize('_')
            [
                [
                    id: sample,
                    replicate: replicate.replace('rep', ''),
                    type: type,
                    readNum: readNum,
                ],
                myFile
            ]
        }
        .view()
    ```

Hem canviat la fàbrica de canals i adaptat el patró de coincidència de fitxers, i de passada hem comentat l'operació map.
L'afegirem de nou més endavant, amb algunes modificacions.

Executeu el workflow per provar-ho:

```bash
nextflow run main.nf
```

??? failure "Sortida de la comanda"

    ```console hl_lines="7-8"
     N E X T F L O W   ~  version 25.10.2

    Launching `main.nf` [angry_koch] DSL2 - revision: 44fdf66105

    [-        ] COUNT_LINES -
    [-        ] COUNT_LINES -
    [patientA_rep1_normal_R, [/workspaces/training/side-quests/working_with_files/data/patientA_rep1_normal_R1_001.fastq.gz, /workspaces/training/side-quests/working_with_files/data/patientA_rep1_normal_R2_001.fastq.gz]]
    ERROR ~ Error executing process > 'COUNT_LINES (1)'

    Caused by:
      Not a valid path value: 'patientA_rep1_normal_R'



    Tip: when you have fixed the problem you can continue the execution adding the option `-resume` to the run command line

    -- Check '.nextflow.log' file for details
    ```

Vaja, aquesta vegada l'execució ha fallat!

La part rellevant del missatge d'error és aquí:

```console
Not a valid path value: 'patientA_rep1_normal_R'
```

Això és perquè hem canviat la fàbrica de canals.
Fins ara, el canal d'entrada original només contenia les rutes dels fitxers.
Tota la manipulació de metadades que hem estat fent no ha afectat realment el contingut del canal.

Ara que estem utilitzant la fàbrica de canals `.fromFilePairs`, el contingut del canal resultant és diferent.
Veiem només un element de canal, compost per una tupla que conté dos elements: la part del `simpleName` compartida pels dos fitxers, que serveix com a identificador, i una tupla que conté els dos objectes de fitxer, en el format `id, [ file1, file2 ]`.

Això és fantàstic, perquè Nextflow ha fet la feina difícil d'extreure el nom del pacient examinant el prefix compartit i utilitzant-lo com a identificador del pacient.

No obstant això, trenca el nostre workflow actual.
Si volguéssim continuar executant `COUNT_LINES` de la mateixa manera sense canviar el procés, hauríem d'aplicar una operació de mapping per extreure les rutes dels fitxers.
Però no ho farem, perquè el nostre objectiu final és utilitzar un procés diferent, `ANALYZE_READS`, que gestiona els parells de fitxers adequadament.

Així que simplement comentem (o eliminem) la crida a `COUNT_LINES` i continuem.

=== "Després"

    ```groovy title="main.nf" linenums="26" hl_lines="2"
        // Compta les línies del fitxer
        // COUNT_LINES(ch_files)
    ```

=== "Abans"

    ```groovy title="main.nf" linenums="26" hl_lines="2"
        // Compta les línies del fitxer
        COUNT_LINES(ch_files)
    ```

També podeu comentar o eliminar la declaració include de `COUNT_LINES`, però això no tindrà cap efecte funcional.

Ara executem el workflow de nou:

```bash
nextflow run main.nf
```

??? success "Sortida de la comanda"

    ```console hl_lines="5"
     N E X T F L O W   ~  version 25.10.2

    Launching `main.nf` [fabulous_davinci] DSL2 - revision: 22b53268dc

    [patientA_rep1_normal_R, [/workspaces/training/side-quests/working_with_files/data/patientA_rep1_normal_R1_001.fastq.gz, /workspaces/training/side-quests/working_with_files/data/patientA_rep1_normal_R2_001.fastq.gz]]
    ```

Ara el workflow s'executa amb èxit!

No obstant això, encara hem d'extreure la resta de les metadades del camp `id`.

### 5.2. Extreure i organitzar metadades dels parells de fitxers

La nostra operació `map` anterior no funcionarà perquè no coincideix amb l'estructura de dades, però podem modificar-la perquè funcioni.

Ja tenim accés a l'identificador real del pacient a la cadena que `fromFilePairs()` ha utilitzat com a identificador, de manera que podem utilitzar-la per extreure les metadades sense obtenir el `simpleName` de l'objecte Path com hem fet abans.

Descomenteu l'operació map al workflow i feu les edicions següents:

=== "Després"

    ```groovy title="main.nf" linenums="7" hl_lines="3-4 9 11 13"
        // Carrega fitxers amb channel.fromFilePairs
        ch_files = channel.fromFilePairs('data/patientA_rep1_normal_R{1,2}_001.fastq.gz')
        ch_files.map { id, files ->
            def (sample, replicate, type) = id.tokenize('_')
            [
                [
                    id: sample,
                    replicate: replicate.replace('rep', ''),
                    type: type
                ],
                files
            ]
        }
        .view()
    ```

=== "Abans"

    ```groovy title="main.nf" linenums="7" hl_lines="3-5 11 13"
        // Carrega fitxers amb channel.fromFilePairs
        ch_files = channel.fromFilePairs('data/patientA_rep1_normal_R{1,2}_001.fastq.gz')
        /* Comentem el map per ara, hi tornarem!
        ch_files.map { myFile ->
            def (sample, replicate, type, readNum) = myFile.simpleName.tokenize('_')
            [
                [
                    id: sample,
                    replicate: replicate.replace('rep', ''),
                    type: type,
                    readNum: readNum,
                ],
                myFile
            ]
        }
        */
        .view()
    ```

Aquesta vegada el map comença des de `id, files` en lloc de simplement `myFile`, i `tokenize()` s'aplica a `id` en lloc de a `myFile.simpleName`.

Observeu també que hem eliminat `readNum` de la línia `tokenize()`; qualsevol subcadena que no nombrem específicament (començant des de l'esquerra) es descartarà silenciosament.
Podem fer-ho perquè els fitxers aparellats ara estan estretament associats, de manera que ja no necessitem `readNum` al map de metadades.

Executem el workflow:

```bash
nextflow run main.nf
```

??? success "Sortida de la comanda"

    ```console

    N E X T F L O W   ~  version 25.10.2

    Launching `main.nf` [prickly_stonebraker] DSL2 - revision: f62ab10a3f

    [[id:patientA, replicate:1, type:normal], [/workspaces/training/side-quests/working_with_files/data/patientA_rep1_normal_R1_001.fastq.gz, /workspaces/training/side-quests/working_with_files/data/patientA_rep1_normal_R2_001.fastq.gz]]
    ```

I aquí ho teniu: tenim el map de metadades (`[id:patientA, replicate:1, type:normal]`) a la primera posició de la tupla de sortida, seguit de la tupla de fitxers aparellats, tal com s'havia previst.

Per descomptat, això només recollirà i processarà aquell parell específic de fitxers.
Si voleu experimentar amb el processament de múltiples parells, podeu provar d'afegir comodins al patró d'entrada i veure què passa.
Per exemple, proveu d'utilitzar `data/patientA_rep1_*_R{1,2}_001.fastq.gz`

### Conclusió

- [`channel.fromFilePairs()` troba i aparella automàticament fitxers relacionats](https://www.nextflow.io/docs/latest/reference/channel.html#fromfilepairs)
- Això simplifica la gestió de lectures paired-end al vostre pipeline
- Els fitxers aparellats es poden agrupar com a tuples `[id, [file1, file2]]`
- L'extracció de metadades es pot fer a partir de l'ID del parell de fitxers en lloc de fitxers individuals

---

## 6. Ús d'operacions de fitxers en processos

Ara posem tot això junt en un procés senzill per reforçar com utilitzar operacions de fitxers dins d'un procés Nextflow.

Us proporcionem un mòdul de procés pre-escrit anomenat `ANALYZE_READS` que pren una tupla de metadades i un parell de fitxers d'entrada i els analitza.
Podríem imaginar que això fa alineament de seqüències, o detecció de variants o qualsevol altre pas que tingui sentit per a aquest tipus de dades.

Comencem.

### 6.1. Importar el procés i examinar el codi

Per utilitzar aquest procés al workflow, simplement hem d'afegir una declaració include de mòdul abans del bloc workflow.

Feu la edició següent al workflow:

=== "Després"

    ```groovy title="main.nf" linenums="1" hl_lines="3"
    #!/usr/bin/env nextflow

    include { ANALYZE_READS } from './modules/analyze_reads.nf'

    workflow {
    ```

=== "Abans"

    ```groovy title="main.nf" linenums="1"
    #!/usr/bin/env nextflow

    workflow {
    ```

Podeu obrir el fitxer del mòdul per examinar el seu codi:

```groovy title="modules/analyze_reads.nf - process example" linenums="1"
#!/usr/bin/env nextflow

process ANALYZE_READS {
    tag { meta.id }

    publishDir { "results/${meta.id}" }, mode: 'copy'

    input:
    tuple val(meta), path(files)

    output:
    tuple val(meta.id), path("${meta.id}_stats.txt")

    script:
    """
    echo "Sample metadata: ${meta.id}" > ${meta.id}_stats.txt
    echo "Replicate: ${meta.replicate}" >> ${meta.id}_stats.txt
    echo "Type: ${meta.type}" >> ${meta.id}_stats.txt
    echo "Read 1: ${files[0]}" >> ${meta.id}_stats.txt
    echo "Read 2: ${files[1]}" >> ${meta.id}_stats.txt
    echo "Read 1 size: \$(gunzip -dc ${files[0]} | wc -l | awk '{print \$1/4}') reads" >> ${meta.id}_stats.txt
    echo "Read 2 size: \$(gunzip -dc ${files[1]} | wc -l | awk '{print \$1/4}') reads" >> ${meta.id}_stats.txt
    """
}
```

!!! note "Nota"

    Les directives `tag` i `publishDir` utilitzen la sintaxi de closure (`{ ... }`) en lloc de la interpolació de cadenes (`"${...}"`).
    Això és perquè aquestes directives fan referència a variables d'entrada (`meta`) que no estan disponibles fins al moment d'execució.
    La sintaxi de closure difereix l'avaluació fins que el procés s'executa realment.

!!! note "Nota"

    Anomenem el nostre map de metadades `meta` per convenció.
    Per a una immersió més profunda en els meta maps, consulteu la missió secundària [Metadata and meta maps](../metadata/).

### 6.2. Cridar el procés al workflow

Ara que el procés està disponible per al workflow, podem afegir una crida al procés `ANALYZE_READS` per executar-lo.

Per executar-lo sobre les nostres dades d'exemple, haurem de fer dues coses:

1. Donar un nom al canal remapat
2. Afegir una crida al procés

#### 6.2.1. Donar nom al canal d'entrada remapat

Anteriorment hem aplicat les manipulacions de mapping directament al canal d'entrada.
Per tal d'alimentar el contingut remapat al procés `ANALYZE_READS` (i fer-ho d'una manera clara i fàcil de llegir) volem crear un nou canal anomenat `ch_samples`.

Podem fer-ho utilitzant l'operador [`set`](https://www.nextflow.io/docs/latest/reference/operator.html#set).

Al workflow principal, substituïu l'operador `.view()` per `.set { ch_samples }`, i afegiu una línia per verificar que podem fer referència al canal pel seu nom.

=== "Després"

    ```groovy title="main.nf" linenums="7" hl_lines="14 16-17"
        // Carrega fitxers amb channel.fromFilePairs
        ch_files = channel.fromFilePairs('data/patientA_rep1_normal_R{1,2}_001.fastq.gz')
        ch_files.map { id,  files ->
           def (sample, replicate, type, readNum) = id.tokenize('_')
           [
               [
                   id: sample,
                   replicate: replicate.replace('rep', ''),
                   type: type
               ],
               files
           ]
        }
            .set { ch_samples }

        // Temporal: inspeccionem ch_samples
        ch_samples.view()
    ```

=== "Abans"

    ```groovy title="main.nf" linenums="7" hl_lines="14"
        // Carrega fitxers amb channel.fromFilePairs
        ch_files = channel.fromFilePairs('data/patientA_rep1_normal_R{1,2}_001.fastq.gz')
        ch_files.map { id,  files ->
           def (sample, replicate, type, readNum) = id.tokenize('_')
           [
               [
                   id: sample,
                   replicate: replicate.replace('rep', ''),
                   type: type
               ],
               files
           ]
        }
        .view()
    }
    ```

Executem-ho:

```bash
nextflow run main.nf
```

??? success "Sortida de la comanda"

    ```console
     N E X T F L O W   ~  version 25.10.2

    Launching `main.nf` [goofy_kirch] DSL2 - revision: 3313283e42

    [[id:patientA, replicate:1, type:normal], [/workspaces/training/side-quests/working_with_files/data/patientA_rep1_normal_R1_001.fastq.gz, /workspaces/training/side-quests/working_with_files/data/patientA_rep1_normal_R2_001.fastq.gz]]
    ```

Això confirma que ara podem fer referència al canal pel seu nom.

#### 6.2.2. Cridar el procés sobre les dades

Ara cridem realment el procés `ANALYZE_READS` sobre el canal `ch_samples`.

Al workflow principal, feu els canvis de codi següents:

=== "Després"

    ```groovy title="main.nf" linenums="23"
        // Executa l'anàlisi
        ANALYZE_READS(ch_samples)
    ```

=== "Abans"

    ```groovy title="main.nf" linenums="23"
        // Temporal: inspeccionem ch_samples
        ch_samples.view()
    ```

Executem-ho:

```bash
nextflow run main.nf
```

??? success "Sortida de la comanda"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `./main.nf` [shrivelled_cori] DSL2 - revision: b546a31769

    executor >  local (1)
    [b5/110360] process > ANALYZE_READS (patientA) [100%] 1 of 1 ✔
    ```

Aquest procés està configurat per publicar les seves sortides a un directori `results`, així que feu una ullada allà.

??? abstract "Contingut de directori i fitxer"

    ```console
    results
    └── patientA
        └── patientA_stats.txt
    ```

    ```txt title="patientA_stats.txt"
    Sample metadata: patientA
    Replicate: 1
    Type: normal
    Read 1: patientA_rep1_normal_R1_001.fastq.gz
    Read 2: patientA_rep1_normal_R2_001.fastq.gz
    Read 1 size: 10 reads
    Read 2 size: 10 reads
    ```

El procés ha pres les nostres entrades i ha creat un nou fitxer que conté les metadades del pacient, tal com estava dissenyat.
Excel·lent!

### 6.3. Incloure molts més pacients

Per descomptat, això simplement processa un únic parell de fitxers per a un únic pacient, que no és exactament el tipus d'alt rendiment que espereu obtenir amb Nextflow.
Probablement voldreu processar moltes més dades alhora.

Recordeu que `channel.fromPath()` accepta un _glob_ com a entrada, la qual cosa significa que pot acceptar qualsevol nombre de fitxers que coincideixin amb el patró.
Per tant, si volem incloure tots els pacients, simplement podem modificar la cadena d'entrada per incloure'n més, tal com s'ha indicat de passada anteriorment.

Suposem que volem ser tan inclusius com sigui possible.
Feu les edicions següents al workflow:

=== "Després"

    ```groovy title="main.nf" linenums="7" hl_lines="2"
        // Carrega fitxers amb channel.fromFilePairs
        ch_files = channel.fromFilePairs('data/*_R{1,2}_001.fastq.gz')
    ```

=== "Abans"

    ```groovy title="main.nf" linenums="7" hl_lines="2"
        // Carrega fitxers amb channel.fromFilePairs
        ch_files = channel.fromFilePairs('data/patientA_rep1_normal_R{1,2}_001.fastq.gz')
    ```

Executeu el pipeline de nou:

```bash
nextflow run main.nf
```

??? success "Sortida de la comanda"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `./main.nf` [big_stonebraker] DSL2 - revision: f7f9b8a76c

    executor >  local (8)
    [d5/441891] process > ANALYZE_READS (patientC) [100%] 8 of 8 ✔
    ```

El directori de resultats ara hauria de contenir resultats per a totes les dades disponibles.

??? abstract "Contingut del directori"

    ```console
    results
    ├── patientA
    │   └── patientA_stats.txt
    ├── patientB
    │   └── patientB_stats.txt
    └── patientC
        └── patientC_stats.txt
    ```

Èxit! Hem analitzat tots els pacients d'una sola vegada! Oi?

Potser no.
Si mireu més de prop, tenim un problema: tenim dues rèpliques per al pacientA, però només un fitxer de sortida!
Estem sobreescrivint el fitxer de sortida cada vegada.

### 6.4. Fer únics els fitxers publicats

Com que tenim accés a les metadades del pacient, podem utilitzar-les per fer únics els fitxers publicats incloent metadades diferenciadores, ja sigui a l'estructura de directoris o als propis noms de fitxer.

Feu el canvi següent al workflow:

=== "Després"

    ```groovy title="modules/analyze_reads.nf" linenums="6"
        publishDir { "results/${meta.type}/${meta.id}/${meta.replicate}" }, mode: 'copy'
    ```

=== "Abans"

    ```groovy title="modules/analyze_reads.nf" linenums="6"
        publishDir { "results/${meta.id}" }, mode: 'copy'
    ```

Aquí mostrem l'opció d'utilitzar nivells de directori addicionals per tenir en compte els tipus de mostra i les rèpliques, però podríeu experimentar fent-ho també a nivell del nom de fitxer.

Ara executeu el pipeline una vegada més, però assegureu-vos d'eliminar primer el directori de resultats per tenir un espai de treball net:

```bash
rm -r results
nextflow run main.nf
```

??? success "Sortida de la comanda"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `./main.nf` [insane_swartz] DSL2 - revision: fff18abe6d

    executor >  local (8)
    [e3/449081] process > ANALYZE_READS (patientC) [100%] 8 of 8 ✔
    ```

Comproveu el directori de resultats ara:

??? abstract "Contingut del directori"

    ```console
    results/
    ├── normal
    │   ├── patientA
    │   │   ├── 1
    │   │   │   └── patientA_stats.txt
    │   │   └── 2
    │   │       └── patientA_stats.txt
    │   ├── patientB
    │   │   └── 1
    │   │       └── patientB_stats.txt
    │   └── patientC
    │       └── 1
    │           └── patientC_stats.txt
    └── tumor
        ├── patientA
        │   ├── 1
        │   │   └── patientA_stats.txt
        │   └── 2
        │       └── patientA_stats.txt
        ├── patientB
        │   └── 1
        │       └── patientB_stats.txt
        └── patientC
            └── 1
                └── patientC_stats.txt
    ```

I aquí ho teniu, totes les nostres metadades, ben organitzades. Això és un èxit!

Hi ha molt més que podeu fer un cop teniu les vostres metadades carregades en un map com aquest:

1. Crear directoris de sortida organitzats basats en els atributs del pacient
2. Prendre decisions en els processos basades en les propietats del pacient
3. Dividir, unir i recombinar dades basant-se en valors de metadades

Aquest patró de mantenir les metadades explícites i adjuntes a les dades (en lloc de codificades als noms de fitxer) és una pràctica recomanada fonamental a Nextflow que permet construir workflows d'anàlisi robustos i mantenibles.
Podeu aprendre més sobre això a la missió secundària [Metadata and meta maps](../metadata/).

### Conclusió

- La directiva `publishDir` pot organitzar les sortides basant-se en valors de metadades
- Les metadades en tuples permeten una organització estructurada dels resultats
- Aquest enfocament crea workflows mantenibles amb una clara procedència de les dades
- Els processos poden prendre tuples de metadades i fitxers com a entrada
- La directiva `tag` proporciona identificació del procés als registres d'execució
- L'estructura del workflow separa la creació de canals de l'execució de processos

---

## Resum

En aquesta missió secundària, heu après a treballar amb fitxers a Nextflow, des d'operacions bàsiques fins a tècniques més avançades per gestionar col·leccions de fitxers.

L'aplicació d'aquestes tècniques en el vostre propi treball us permetrà construir workflows més eficients i mantenibles, especialment quan treballeu amb un gran nombre de fitxers amb convencions de nomenclatura complexes.

### Patrons clau

1.  **Operacions bàsiques amb fitxers:** Hem creat objectes Path amb `file()` i hem accedit als atributs dels fitxers com ara el nom, l'extensió i el directori pare, aprenent la diferència entre cadenes i objectes Path.

    - Crear un objecte Path amb `file()`

    ```groovy
    myFile = file('path/to/file.txt')
    ```

    - Obtenir atributs del fitxer

    ```groovy
    println myFile.name       // file.txt
    println myFile.baseName   // file
    println myFile.extension  // txt
    println myFile.parent     // path/to
    ```

2.  **Ús de fitxers remots**: Hem après a canviar de manera transparent entre fitxers locals i remots utilitzant URIs, demostrant la capacitat de Nextflow per gestionar fitxers de diverses fonts sense canviar la lògica del workflow.

    - Fitxer local

    ```groovy
    myFile = file('path/to/file.txt')
    ```

    - FTP

    ```groovy
    myFile = file('ftp://path/to/file.txt')
    ```

    - HTTPS

    ```groovy
    myFile = file('https://path/to/file.txt')
    ```

    - Amazon S3

    ```groovy
    myFile = file('s3://path/to/file.txt')
    ```

    - Azure Blob Storage

    ```groovy
    myFile = file('az://path/to/file.txt')
    ```

    - Google Cloud Storage

    ```groovy
    myFile = file('gs://path/to/file.txt')
    ```

3.  **Càrrega de fitxers amb la fàbrica de canals `fromPath()`:** Hem creat canals a partir de patrons de fitxers amb `channel.fromPath()` i hem vist els seus atributs de fitxer, incloent els tipus d'objectes.

    - Crear un canal a partir d'un patró de fitxer

    ```groovy
     ch_files = channel.fromPath('data/*.fastq.gz')
    ```

    - Obtenir atributs del fitxer

    ```groovy
     ch_files.view { myFile ->
        println "File object class: ${myFile.class}"
        println "File name: ${myFile.name}"
        println "Simple name: ${myFile.simpleName}"
        println "Extension: ${myFile.extension}"
        println "Parent directory: ${myFile.parent}"
    }
    ```

4.  **Extracció de metadades del pacient dels noms de fitxer:** Hem utilitzat `tokenize()` i `replace()` per extreure i estructurar metadades dels noms de fitxer, convertint-les en maps organitzats.

    ```groovy
    def name = file.name.tokenize('_')
    def patientId = name[0]
    def replicate = name[1].replace('rep', '')
    def type = name[2]
    def readNum = name[3].replace('R', '')
    ```

5.  **Simplificació amb channel.fromFilePairs:** Hem utilitzat `channel.fromFilePairs()` per aparellar automàticament fitxers relacionats i extreure metadades dels IDs dels parells de fitxers.

    ```groovy
    ch_pairs = channel.fromFilePairs('data/*_R{1,2}_001.fastq.gz')
    ```

6.  **Ús d'operacions de fitxers en processos:** Hem integrat operacions de fitxers en processos Nextflow amb una gestió correcta de l'entrada, utilitzant `publishDir` per organitzar les sortides basant-se en metadades.

    - Associar un meta map amb les entrades del procés

    ```groovy
    ch_files = channel.fromFilePairs('data/patientA_rep1_normal_R{1,2}_001.fastq.gz')
    ch_files.map { id,  files ->
        def (sample, replicate, type, readNum) = id.tokenize('_')
        [
            [
                id: sample,
                replicate: replicate.replace('rep', ''),
                type: type
            ],
             files
        ]
    }
        .set { ch_samples }

    ANALYZE_READS(ch_samples)
    ```

    - Organitzar les sortides basant-se en metadades

    ```groovy
    publishDir { "results/${meta.type}/${meta.id}/${meta.replicate}" }, mode: 'copy'
    ```

### Recursos addicionals

- [Documentació de Nextflow: Treballar amb fitxers](https://www.nextflow.io/docs/latest/working-with-files.html)
- [channel.fromPath](https://www.nextflow.io/docs/latest/reference/channel.html#frompath)
- [channel.fromFilePairs](https://www.nextflow.io/docs/latest/reference/channel.html#fromfilepairs)

---

## Què segueix?

Torneu al [menú de missions secundàries](../) o feu clic al botó a la part inferior dreta de la pàgina per continuar amb el tema següent de la llista.
