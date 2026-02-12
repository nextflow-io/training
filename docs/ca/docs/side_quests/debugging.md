# Depuració de Workflows

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Traducció assistida per IA - [més informació i suggeriments](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

La depuració és una habilitat crítica que us pot estalviar hores de frustració i ajudar-vos a convertir-vos en un desenvolupador de Nextflow més efectiu. Al llarg de la vostra carrera, especialment quan esteu començant, trobareu errors mentre construïu i manteniu els vostres workflows. Aprendre enfocaments sistemàtics de depuració us ajudarà a identificar i resoldre problemes ràpidament.

### Objectius d'aprenentatge

En aquesta missió secundària, explorarem **tècniques de depuració sistemàtiques** per a workflows de Nextflow:

- **Depuració d'errors de sintaxi**: Utilitzar funcionalitats de l'IDE i missatges d'error de Nextflow de manera efectiva
- **Depuració de canals**: Diagnosticar problemes de flux de dades i problemes d'estructura de canals
- **Depuració de processos**: Investigar fallades d'execució i problemes de recursos
- **Eines de depuració integrades**: Aprofitar el mode de previsualització de Nextflow, execució stub i directoris de treball
- **Enfocaments sistemàtics**: Una metodologia de quatre fases per a una depuració eficient

Al final, tindreu una metodologia de depuració robusta que transforma missatges d'error frustrants en fulls de ruta clars cap a solucions.

### Prerequisits

Abans d'emprendre aquesta missió secundària, hauríeu de:

- Haver completat el tutorial [Hello Nextflow](../hello_nextflow/README.md) o un curs equivalent per a principiants.
- Sentir-vos còmodes utilitzant conceptes i mecanismes bàsics de Nextflow (processos, canals, operadors)

**Opcional:** Recomanem completar primer la missió secundària [IDE Features for Nextflow Development](./ide_features.md).
Aquesta cobreix de manera exhaustiva les funcionalitats de l'IDE que donen suport a la depuració (ressaltat de sintaxi, detecció d'errors, etc.), que utilitzarem intensament aquí.

---

## 0. Primers passos

#### Obrir el codespace de formació

Si encara no ho heu fet, assegureu-vos d'obrir l'entorn de formació tal com es descriu a [Configuració de l'entorn](../envsetup/index.md).

[![Open in GitHub Codespaces](https://github.com/codespaces/badge.svg)](https://codespaces.new/nextflow-io/training?quickstart=1&ref=master)

#### Moure's al directori del projecte

Anem al directori on es troben els fitxers per a aquest tutorial.

```bash
cd side-quests/debugging
```

Podeu configurar VSCode perquè se centri en aquest directori:

```bash
code .
```

#### Revisar els materials

Trobareu un conjunt de workflows d'exemple amb diversos tipus d'errors que utilitzarem per a la pràctica:

??? abstract "Contingut del directori"

    ```console
    .
    ├── bad_bash_var.nf
    ├── bad_channel_shape.nf
    ├── bad_channel_shape_viewed_debug.nf
    ├── bad_channel_shape_viewed.nf
    ├── bad_number_inputs.nf
    ├── badpractice_syntax.nf
    ├── bad_resources.nf
    ├── bad_syntax.nf
    ├── buggy_workflow.nf
    ├── data
    │   ├── sample_001.fastq.gz
    │   ├── sample_002.fastq.gz
    │   ├── sample_003.fastq.gz
    │   ├── sample_004.fastq.gz
    │   ├── sample_005.fastq.gz
    │   └── sample_data.csv
    ├── exhausted.nf
    ├── invalid_process.nf
    ├── missing_output.nf
    ├── missing_software.nf
    ├── missing_software_with_stub.nf
    ├── nextflow.config
    └── no_such_var.nf
    ```

Aquests fitxers representen escenaris de depuració comuns que trobareu en el desenvolupament real.

#### Revisar l'assignació

El vostre repte és executar cada workflow, identificar l'error o errors, i corregir-los.

Per a cada workflow amb errors:

1. **Executar el workflow** i observar l'error
2. **Analitzar el missatge d'error**: què us està dient Nextflow?
3. **Localitzar el problema** al codi utilitzant les pistes proporcionades
4. **Corregir l'error** i verificar que la vostra solució funciona
5. **Restablir el fitxer** abans de passar a la següent secció (utilitzeu `git checkout <filename>`)

Els exercicis progressen des d'errors de sintaxi simples fins a problemes d'execució més subtils.
Les solucions es discuteixen en línia, però intenteu resoldre cada un vosaltres mateixos abans de llegir més endavant.

#### Llista de verificació de preparació

Creieu que esteu preparats per començar?

- [ ] Entenc l'objectiu d'aquest curs i els seus prerequisits
- [ ] El meu codespace està en funcionament
- [ ] He establert el meu directori de treball adequadament
- [ ] Entenc l'assignació

Si podeu marcar totes les caselles, esteu preparats per començar.

---

## 1. Errors de Sintaxi

Els errors de sintaxi són el tipus d'error més comú que trobareu quan escriviu codi Nextflow. Es produeixen quan el codi no es conforma a les regles de sintaxi esperades del DSL de Nextflow. Aquests errors impedeixen que el vostre workflow s'executi del tot, així que és important aprendre a identificar-los i corregir-los ràpidament.

### 1.1. Claus que falten

Un dels errors de sintaxi més comuns, i de vegades un dels més complexos de depurar, són les **claus que falten o no coincideixen**.

Comencem amb un exemple pràctic.

#### Executar el pipeline

```bash
nextflow run bad_syntax.nf
```

??? failure "Sortida de la comanda"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `bad_syntax.nf` [stupefied_bhabha] DSL2 - revision: ca6327fad2

    Error bad_syntax.nf:24:1: Unexpected input: '<EOF>'

    ERROR ~ Script compilation failed

     -- Check '.nextflow.log' file for details
    ```

**Elements clau dels missatges d'error de sintaxi:**

- **Fitxer i ubicació**: Mostra quin fitxer i línia/columna contenen l'error (`bad_syntax.nf:24:1`)
- **Descripció de l'error**: Explica què ha trobat l'analitzador que no esperava (`Unexpected input: '<EOF>'`)
- **Indicador EOF**: El missatge `<EOF>` (End Of File) indica que l'analitzador ha arribat al final del fitxer mentre encara esperava més contingut - un signe clàssic de claus no tancades

#### Comprovar el codi

Ara, examinem `bad_syntax.nf` per entendre què està causant l'error:

```groovy title="bad_syntax.nf" hl_lines="14" linenums="1"
#!/usr/bin/env nextflow

process PROCESS_FILES {
    input:
    val sample_name

    output:
    path "${sample_name}_output.txt"

    script:
    """
    echo "Processing ${sample_name}" > ${sample_name}_output.txt
    """
// Missing closing brace for the process

workflow {

    // Create input channel
    input_ch = channel.of('sample1', 'sample2', 'sample3')

    // Call the process with the input channel
    PROCESS_FILES(input_ch)
}
```

Per a l'objectiu d'aquest exemple hem deixat un comentari perquè vegeu on és l'error. L'extensió de Nextflow per a VSCode també us hauria de donar algunes pistes sobre què pot estar malament, posant la clau no coincident en vermell i ressaltant el final prematur del fitxer:

![Bad syntax](img/bad_syntax.png)

**Estratègia de depuració per a errors de claus:**

1. Utilitzar la coincidència de claus de VS Code (col·locar el cursor al costat d'una clau)
2. Comprovar el panell de Problemes per a missatges relacionats amb claus
3. Assegurar-se que cada `{` d'obertura té una `}` de tancament corresponent

#### Corregir el codi

Substituïu el comentari amb la clau de tancament que falta:

=== "Després"

    ```groovy title="bad_syntax.nf" hl_lines="14" linenums="1"
    #!/usr/bin/env nextflow

    process PROCESS_FILES {
        input:
        val sample_name

        output:
        path "${sample_name}_output.txt"

        script:
        """
        echo "Processing ${sample_name}" > ${sample_name}_output.txt
        """
    }  // Afegir la clau de tancament que falta

    workflow {

        // Create input channel
        input_ch = channel.of('sample1', 'sample2', 'sample3')

        // Call the process with the input channel
        PROCESS_FILES(input_ch)
    }
    ```

=== "Abans"

    ```groovy title="bad_syntax.nf" hl_lines="14" linenums="1"
    #!/usr/bin/env nextflow

    process PROCESS_FILES {
        input:
        val sample_name

        output:
        path "${sample_name}_output.txt"

        script:
        """
        echo "Processing ${sample_name}" > ${sample_name}_output.txt
        """
    // Missing closing brace for the process

    workflow {

        // Create input channel
        input_ch = channel.of('sample1', 'sample2', 'sample3')

        // Call the process with the input channel
        PROCESS_FILES(input_ch)
    }
    ```

#### Executar el pipeline

Ara executeu el workflow de nou per confirmar que funciona:

```bash
nextflow run bad_syntax.nf
```

??? success "Sortida de la comanda"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `bad_syntax.nf` [insane_faggin] DSL2 - revision: 961938ee2b

    executor >  local (3)
    [48/cd7f54] PROCESS_FILES (1) | 3 of 3 ✔
    ```

### 1.2. Utilitzar paraules clau o directives de procés incorrectes

Un altre error de sintaxi comú és una **definició de procés invàlida**. Això pot passar si oblideu definir blocs requerits o utilitzeu directives incorrectes a la definició del procés.

#### Executar el pipeline

```bash
nextflow run invalid_process.nf
```

??? failure "Sortida de la comanda"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `invalid_process.nf` [nasty_jepsen] DSL2 - revision: da9758d614

    Error invalid_process.nf:3:1: Invalid process definition -- check for missing or out-of-order section labels
    │   3 | process PROCESS_FILES {
    │     | ^^^^^^^^^^^^^^^^^^^^^^^
    │   4 |     inputs:
    │   5 |     val sample_name
    │   6 |
    ╰   7 |     output:

    ERROR ~ Script compilation failed

     -- Check '.nextflow.log' file for details
    ```

#### Comprovar el codi

L'error indica una "Invalid process definition" i mostra el context al voltant del problema. Mirant les línies 3-7, podem veure `inputs:` a la línia 4, que és el problema. Examinem `invalid_process.nf`:

```groovy title="invalid_process.nf" hl_lines="4" linenums="1"
#!/usr/bin/env nextflow

process PROCESS_FILES {
    inputs:  // ERROR: Hauria de ser 'input' no 'inputs'
    val sample_name

    output:
    path "${sample_name}_output.txt"

    script:
    """
    echo "Processing ${sample_name}" > ${sample_name}_output.txt
    """
}

workflow {

    // Create input channel
    input_ch = channel.of('sample1', 'sample2', 'sample3')

    // Call the process with the input channel
    PROCESS_FILES(input_ch)
}
```

Mirant la línia 4 al context de l'error, podem detectar el problema: estem utilitzant `inputs` en lloc de la directiva correcta `input`. L'extensió de Nextflow per a VSCode també marcarà això:

![Invalid process message](img/invalid_process_message.png)

#### Corregir el codi

Substituïu la paraula clau incorrecta per la correcta consultant [la documentació](https://www.nextflow.io/docs/latest/process.html#):

=== "Després"

    ```groovy title="invalid_process.nf" hl_lines="4" linenums="1"
    #!/usr/bin/env nextflow

    process PROCESS_FILES {
        input:  // Corregit: Canviat 'inputs' a 'input'
        val sample_name

        output:
        path "${sample_name}_output.txt"

        script:
        """
        echo "Processing ${sample_name}" > ${sample_name}_output.txt
        """
    }

    workflow {

        // Create input channel
        input_ch = channel.of('sample1', 'sample2', 'sample3')

        // Call the process with the input channel
        PROCESS_FILES(input_ch)
    }
    ```

=== "Abans"

    ```groovy title="invalid_process.nf" hl_lines="4" linenums="1"
    #!/usr/bin/env nextflow

    process PROCESS_FILES {
        inputs:  // ERROR: Hauria de ser 'input' no 'inputs'
        val sample_name

        output:
        path "${sample_name}_output.txt"

        script:
        """
        echo "Processing ${sample_name}" > ${sample_name}_output.txt
        """
    }

    workflow {

        // Create input channel
        input_ch = channel.of('sample1', 'sample2', 'sample3')

        // Call the process with the input channel
        PROCESS_FILES(input_ch)
    }
    ```

#### Executar el pipeline

Ara executeu el workflow de nou per confirmar que funciona:

```bash
nextflow run invalid_process.nf
```

??? success "Sortida de la comanda"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `invalid_process.nf` [silly_fermi] DSL2 - revision: 961938ee2b

    executor >  local (3)
    [b7/76cd9d] PROCESS_FILES (2) | 3 of 3 ✔
    ```

### 1.3. Utilitzar noms de variables incorrectes

Els noms de variables que utilitzeu als vostres blocs de script han de ser vàlids, derivats ja sigui d'entrades o de codi groovy inserit abans de l'script. Però quan esteu gestionant complexitat al començament del desenvolupament del pipeline, és fàcil cometre errors en la denominació de variables, i Nextflow us ho farà saber ràpidament.

#### Executar el pipeline

```bash
nextflow run no_such_var.nf
```

??? failure "Sortida de la comanda"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `no_such_var.nf` [gloomy_meninsky] DSL2 - revision: 0c4d3bc28c

    Error no_such_var.nf:17:39: `undefined_var` is not defined
    │  17 |     echo "Using undefined variable: ${undefined_var}" >> ${output_pref
    ╰     |                                       ^^^^^^^^^^^^^

    ERROR ~ Script compilation failed

     -- Check '.nextflow.log' file for details
    ```

L'error es detecta en temps de compilació i apunta directament a la variable no definida a la línia 17, amb un accent circumflex indicant exactament on és el problema.

#### Comprovar el codi

Examinem `no_such_var.nf`:

```groovy title="no_such_var.nf" hl_lines="17" linenums="1"
#!/usr/bin/env nextflow

process PROCESS_FILES {
    input:
    val sample_name

    output:
    path "${sample_name}_processed.txt"

    script:
    // Definir variables en codi Groovy abans de l'script
    def output_prefix = "${sample_name}_processed"
    def timestamp = new Date().format("yyyy-MM-dd")

    """
    echo "Processing ${sample_name} on ${timestamp}" > ${output_prefix}.txt
    echo "Using undefined variable: ${undefined_var}" >> ${output_prefix}.txt  // ERROR: undefined_var no està definida
    """
}

workflow {
    input_ch = channel.of('sample1', 'sample2', 'sample3')
    PROCESS_FILES(input_ch)
}
```

El missatge d'error indica que la variable no es reconeix a la plantilla de l'script, i allà ho teniu - hauríeu de poder veure `${undefined_var}` utilitzada al bloc de script, però no definida en cap altre lloc.

#### Corregir el codi

Si obteniu un error 'No such variable', podeu corregir-lo definint la variable (corregint noms de variables d'entrada o editant codi groovy abans de l'script), o eliminant-la del bloc de script si no és necessària:

=== "Després"

    ```groovy title="no_such_var.nf" hl_lines="15-17" linenums="1"
    #!/usr/bin/env nextflow

    process PROCESS_FILES {
        input:
        val sample_name

        output:
        path "${sample_name}_output.txt"

        script:
        // Definir variables en codi Groovy abans de l'script
        def output_prefix = "${sample_name}_processed"
        def timestamp = new Date().format("yyyy-MM-dd")

        """
        echo "Processing ${sample_name} on ${timestamp}" > ${output_prefix}.txt
        """  // Eliminada la línia amb undefined_var
    }

    workflow {
        input_ch = channel.of('sample1', 'sample2', 'sample3')
        PROCESS_FILES(input_ch)
    }
    ```

=== "Abans"

    ```groovy title="no_such_var.nf" hl_lines="17" linenums="1"
    #!/usr/bin/env nextflow

    process PROCESS_FILES {
        input:
        val sample_name

        output:
        path "${sample_name}_output.txt"

        script:
        // Definir variables en codi Groovy abans de l'script
        def output_prefix = "${sample_name}_processed"
        def timestamp = new Date().format("yyyy-MM-dd")

        """
        echo "Processing ${sample_name} on ${timestamp}" > ${output_prefix}.txt
        echo "Using undefined variable: ${undefined_var}" >> ${output_prefix}.txt  // ERROR: undefined_var no està definida
        """
    }

    workflow {
        input_ch = channel.of('sample1', 'sample2', 'sample3')
        PROCESS_FILES(input_ch)
    }
    ```

#### Executar el pipeline

Ara executeu el workflow de nou per confirmar que funciona:

```bash
nextflow run no_such_var.nf
```

??? success "Sortida de la comanda"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `no_such_var.nf` [suspicious_venter] DSL2 - revision: 6ba490f7c5

    executor >  local (3)
    [21/237300] PROCESS_FILES (2) | 3 of 3 ✔
    ```

### 1.4. Mal ús de variables Bash

Quan comenceu amb Nextflow, pot ser difícil entendre la diferència entre variables de Nextflow (Groovy) i Bash. Això pot generar una altra forma de l'error de variable incorrecta que apareix quan s'intenten utilitzar variables al contingut Bash del bloc de script.

#### Executar el pipeline

```bash
nextflow run bad_bash_var.nf
```

??? failure "Sortida de la comanda"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `bad_bash_var.nf` [infallible_mandelbrot] DSL2 - revision: 0853c11080

    Error bad_bash_var.nf:13:42: `prefix` is not defined
    │  13 |     echo "Processing ${sample_name}" > ${prefix}.txt
    ╰     |                                          ^^^^^^

    ERROR ~ Script compilation failed

     -- Check '.nextflow.log' file for details
    ```

#### Comprovar el codi

L'error apunta a la línia 13 on s'utilitza `${prefix}`. Examinem `bad_bash_var.nf` per veure què està causant el problema:

```groovy title="bad_bash_var.nf" hl_lines="13" linenums="1"
#!/usr/bin/env nextflow

process PROCESS_FILES {
    input:
    val sample_name

    output:
    path "${sample_name}_output.txt"

    script:
    """
    prefix="${sample_name}_output"
    echo "Processing ${sample_name}" > ${prefix}.txt  # ERROR: ${prefix} és sintaxi Groovy, no Bash
    """
}
```

En aquest exemple, estem definint la variable `prefix` en Bash, però en un procés Nextflow la sintaxi `$` que hem utilitzat per referir-nos-hi (`${prefix}`) s'interpreta com una variable Groovy, no Bash. La variable no existeix al context Groovy, així que obtenim un error 'no such variable'.

#### Corregir el codi

Si voleu utilitzar una variable Bash, heu d'escapar el signe de dòlar així:

=== "Després"

    ```groovy title="bad_bash_var.nf" hl_lines="13" linenums="1"
    #!/usr/bin/env nextflow

    process PROCESS_FILES {
        input:
        val sample_name

        output:
        path "${sample_name}_output.txt"

        script:
        """
        prefix="${sample_name}_output"
        echo "Processing ${sample_name}" > \${prefix}.txt  # Corregit: Escapat el signe de dòlar
        """
    }

    workflow {
        input_ch = channel.of('sample1', 'sample2', 'sample3')
        PROCESS_FILES(input_ch)
    }
    ```

=== "Abans"

    ```groovy title="bad_bash_var.nf" hl_lines="13" linenums="1"
    #!/usr/bin/env nextflow

    process PROCESS_FILES {
        input:
        val sample_name

        output:
        path "${sample_name}_output.txt"

        script:
        """
        prefix="${sample_name}_output"
        echo "Processing ${sample_name}" > ${prefix}.txt  # ERROR: ${prefix} és sintaxi Groovy, no Bash
        """
    }
    ```

Això indica a Nextflow que interpreti això com una variable Bash.

#### Executar el pipeline

Ara executeu el workflow de nou per confirmar que funciona:

```bash
nextflow run bad_bash_var.nf
```

??? success "Sortida de la comanda"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `bad_bash_var.nf` [naughty_franklin] DSL2 - revision: 58c1c83709

    executor >  local (3)
    [4e/560285] PROCESS_FILES (2) | 3 of 3 ✔
    ```

!!! tip "Variables Groovy vs Bash"

    Per a manipulacions de variables simples com concatenació de cadenes o operacions de prefix/sufix, normalment és més llegible utilitzar variables Groovy a la secció script en lloc de variables Bash al bloc script:

    ```groovy linenums="1"
    script:
    def output_prefix = "${sample_name}_processed"
    def output_file = "${output_prefix}.txt"
    """
    echo "Processing ${sample_name}" > ${output_file}
    """
    ```

    Aquest enfocament evita la necessitat d'escapar signes de dòlar i fa que el codi sigui més fàcil de llegir i mantenir.

### 1.5. Sentències fora del bloc Workflow

L'extensió de Nextflow per a VSCode ressalta problemes amb l'estructura del codi que causaran errors. Un exemple comú és definir canals fora del bloc `workflow {}` - això ara s'aplica com un error de sintaxi.

#### Executar el pipeline

```bash
nextflow run badpractice_syntax.nf
```

??? failure "Sortida de la comanda"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `badpractice_syntax.nf` [intergalactic_colden] DSL2 - revision: 5e4b291bde

    Error badpractice_syntax.nf:3:1: Statements cannot be mixed with script declarations -- move statements into a process or workflow
    │   3 | input_ch = channel.of('sample1', 'sample2', 'sample3')
    ╰     | ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

    ERROR ~ Script compilation failed

     -- Check '.nextflow.log' file for details
    ```

El missatge d'error indica clarament el problema: les sentències (com definicions de canals) no es poden barrejar amb declaracions de script fora d'un bloc workflow o process.

#### Comprovar el codi

Examinem `badpractice_syntax.nf` per veure què està causant l'error:

```groovy title="badpractice_syntax.nf" hl_lines="3" linenums="1"
#!/usr/bin/env nextflow

input_ch = channel.of('sample1', 'sample2', 'sample3')  // ERROR: Canal definit fora del workflow

process PROCESS_FILES {
    input:
    val sample_name

    output:
    path "${sample_name}_processed.txt"

    script:
    // Definir variables en codi Groovy abans de l'script
    def output_prefix = "${sample_name}_processed"
    def timestamp = new Date().format("yyyy-MM-dd")

    """
    echo "Processing ${sample_name} on ${timestamp}" > ${output_prefix}.txt
    """
}

workflow {
    PROCESS_FILES(input_ch)
}
```

L'extensió VSCode també ressaltarà la variable `input_ch` com definida fora del bloc workflow:

![Non-lethal syntax error](img/nonlethal.png)

#### Corregir el codi

Moveu la definició del canal dins del bloc workflow:

=== "Després"

    ```groovy title="badpractice_syntax.nf" hl_lines="21" linenums="1"
    #!/usr/bin/env nextflow

    process PROCESS_FILES {
        input:
        val sample_name

        output:
        path "${sample_name}_processed.txt"

        script:
        // Definir variables en codi Groovy abans de l'script
        def output_prefix = "${sample_name}_processed"
        def timestamp = new Date().format("yyyy-MM-dd")

        """
        echo "Processing ${sample_name} on ${timestamp}" > ${output_prefix}.txt
        """
    }

    workflow {
        input_ch = channel.of('sample1', 'sample2', 'sample3')  // Mogut dins del bloc workflow
        PROCESS_FILES(input_ch)
    }
    ```

=== "Abans"

    ```groovy title="badpractice_syntax.nf" hl_lines="3" linenums="1"
    #!/usr/bin/env nextflow

    input_ch = channel.of('sample1', 'sample2', 'sample3')  // ERROR: Canal definit fora del workflow

    process PROCESS_FILES {
        input:
        val sample_name

        output:
        path "${sample_name}_processed.txt"

        script:
        // Definir variables en codi Groovy abans de l'script
        def output_prefix = "${sample_name}_processed"
        def timestamp = new Date().format("yyyy-MM-dd")

        """
        echo "Processing ${sample_name} on ${timestamp}" > ${output_prefix}.txt
        """
    }

    workflow {
        PROCESS_FILES(input_ch)
    }
    ```

#### Executar el pipeline

Executeu el workflow de nou per confirmar que la correcció funciona:

```bash
nextflow run badpractice_syntax.nf
```

??? success "Sortida de la comanda"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `badpractice_syntax.nf` [naughty_ochoa] DSL2 - revision: 5e4b291bde

    executor >  local (3)
    [6a/84a608] PROCESS_FILES (2) | 3 of 3 ✔
    ```

Mantingueu els vostres canals d'entrada definits dins del bloc workflow, i en general seguiu qualsevol altra recomanació que faci l'extensió.

### Conclusió

Podeu identificar i corregir errors de sintaxi sistemàticament utilitzant missatges d'error de Nextflow i indicadors visuals de l'IDE. Els errors de sintaxi comuns inclouen claus que falten, paraules clau de procés incorrectes, variables no definides i ús inadequat de variables Bash vs. Nextflow. L'extensió VSCode ajuda a detectar molts d'aquests abans de l'execució. Amb aquestes habilitats de depuració de sintaxi al vostre conjunt d'eines, podreu resoldre ràpidament els errors de sintaxi de Nextflow més comuns i passar a abordar problemes d'execució més complexos.

### Què segueix?

Apreneu a depurar errors d'estructura de canals més complexos que es produeixen fins i tot quan la sintaxi és correcta.

---

## 2. Errors d'Estructura de Canals

Els errors d'estructura de canals són més subtils que els errors de sintaxi perquè el codi és sintàcticament correcte, però les formes de les dades no coincideixen amb el que els processos esperen. Nextflow intentarà executar el pipeline, però pot trobar que el nombre d'entrades no coincideix amb el que espera i fallar. Aquests errors normalment només apareixen en temps d'execució i requereixen una comprensió de les dades que flueixen pel vostre workflow.

!!! tip "Depuració de Canals amb `.view()`"

    Al llarg d'aquesta secció, recordeu que podeu utilitzar l'operador `.view()` per inspeccionar el contingut del canal en qualsevol punt del vostre workflow. Aquesta és una de les eines de depuració més potents per entendre problemes d'estructura de canals. Explorarem aquesta tècnica en detall a la secció 2.4, però sentiu-vos lliures d'utilitzar-la mentre treballeu amb els exemples.

    ```groovy
    my_channel.view()  // Mostra què està fluint pel canal
    ```

### 2.1. Nombre incorrecte de canals d'entrada

Aquest error es produeix quan passeu un nombre diferent de canals del que un procés espera.

#### Executar el pipeline

```bash
nextflow run bad_number_inputs.nf
```

??? failure "Sortida de la comanda"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `bad_number_inputs.nf` [happy_swartz] DSL2 - revision: d83e58dcd3

    Error bad_number_inputs.nf:23:5: Incorrect number of call arguments, expected 1 but received 2
    │  23 |     PROCESS_FILES(samples_ch, files_ch)
    ╰     |     ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

    ERROR ~ Script compilation failed

     -- Check '.nextflow.log' file for details
    ```

#### Comprovar el codi

El missatge d'error indica clarament que la crida esperava 1 argument però va rebre 2, i apunta a la línia 23. Examinem `bad_number_inputs.nf`:

```groovy title="bad_number_inputs.nf" hl_lines="5 23" linenums="1"
#!/usr/bin/env nextflow

process PROCESS_FILES {
    input:
        val sample_name  // El procés espera només 1 entrada

    output:
        path "${sample_name}_output.txt"

    script:
    """
    echo "Processing ${sample_name}" > ${sample_name}_output.txt
    """
}

workflow {

    // Crear dos canals separats
    samples_ch = channel.of('sample1', 'sample2', 'sample3')
    files_ch = channel.of('file1.txt', 'file2.txt', 'file3.txt')

    // ERROR: Passant 2 canals però el procés només espera 1
    PROCESS_FILES(samples_ch, files_ch)
}
```

Hauríeu de veure la crida `PROCESS_FILES` no coincident, subministrant múltiples canals d'entrada quan el procés només en defineix un. L'extensió VSCode també subratlla la crida del procés en vermell, i proporciona un missatge de diagnòstic quan passeu el ratolí per sobre:

![Incorrect number of args message](img/incorrect_num_args.png)

#### Corregir el codi

Per a aquest exemple específic, el procés espera un sol canal i no requereix el segon canal, així que podem corregir-ho passant només el canal `samples_ch`:

=== "Després"

    ```groovy title="bad_number_inputs.nf" hl_lines="23" linenums="1"
    #!/usr/bin/env nextflow

    process PROCESS_FILES {
        input:
            val sample_name  // El procés espera només 1 entrada

        output:
            path "${sample_name}_output.txt"

        script:
        """
        echo "Processing ${sample_name}" > ${sample_name}_output.txt
        """
    }

    workflow {

        // Crear dos canals separats
        samples_ch = channel.of('sample1', 'sample2', 'sample3')
        files_ch = channel.of('file1.txt', 'file2.txt', 'file3.txt')

        // Corregit: Passar només el canal que el procés espera
        PROCESS_FILES(samples_ch)
    }
    ```

=== "Abans"

    ```groovy title="bad_number_inputs.nf" hl_lines="5 23" linenums="1"
    #!/usr/bin/env nextflow

    process PROCESS_FILES {
        input:
            val sample_name  // El procés espera només 1 entrada

        output:
            path "${sample_name}_output.txt"

        script:
        """
        echo "Processing ${sample_name}" > ${sample_name}_output.txt
        """
    }

    workflow {

        // Crear dos canals separats
        samples_ch = channel.of('sample1', 'sample2', 'sample3')
        files_ch = channel.of('file1.txt', 'file2.txt', 'file3.txt')

        // ERROR: Passant 2 canals però el procés només espera 1
        PROCESS_FILES(samples_ch, files_ch)
    }
    ```

#### Executar el pipeline

```bash
nextflow run bad_number_inputs.nf
```

??? success "Sortida de la comanda"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `bad_number_inputs.nf` [big_euler] DSL2 - revision: e302bd87be

    executor >  local (3)
    [48/497f7b] PROCESS_FILES (3) | 3 of 3 ✔
    ```

Més comunament que aquest exemple, podríeu afegir entrades addicionals a un procés i oblidar actualitzar la crida del workflow en conseqüència, cosa que pot portar a aquest tipus d'error. Afortunadament, aquest és un dels errors més fàcils d'entendre i corregir, ja que el missatge d'error és força clar sobre la discrepància.

### 2.2. Esgotament de Canal (El Procés s'Executa Menys Vegades del que s'Esperava)

Alguns errors d'estructura de canals són molt més subtils i no produeixen cap error. Probablement el més comú d'aquests reflecteix un repte que els nous usuaris de Nextflow afronten en entendre que els canals de cua poden esgotar-se i quedar-se sense elements, cosa que significa que el workflow acaba prematurament.

#### Executar el pipeline

```bash
nextflow run exhausted.nf
```

??? success "Sortida de la comanda"

```console title="Exhausted channel output"
 N E X T F L O W   ~  version 25.10.2

Launching `exhausted.nf` [extravagant_gauss] DSL2 - revision: 08cff7ba2a

executor >  local (1)
[bd/f61fff] PROCESS_FILES (1) [100%] 1 of 1 ✔
```

Aquest workflow es completa sense errors, però només processa una sola mostra!

#### Comprovar el codi

Examinem `exhausted.nf` per veure si això és correcte:

```groovy title="exhausted.nf" hl_lines="23 24" linenums="1"
#!/usr/bin/env nextflow

process PROCESS_FILES {
    input:
    val reference
    val sample_name

    output:
    path "${output_prefix}.txt"

    script:
    // Definir variables en codi Groovy abans de l'script
    output_prefix = "${reference}_${sample_name}"
    def timestamp = new Date().format("yyyy-MM-dd")

    """
    echo "Processing ${sample_name} on ${timestamp}" > ${output_prefix}.txt
    """
}

workflow {

    reference_ch = channel.of('baseline_reference')
    input_ch = channel.of('sample1', 'sample2', 'sample3')

    PROCESS_FILES(reference_ch, input_ch)
}
```

El procés només s'executa una vegada en lloc de tres perquè el canal `reference_ch` és un canal de cua que s'esgota després de la primera execució del procés. Quan un canal s'esgota, tot el procés s'atura, fins i tot si altres canals encara tenen elements.

Aquest és un patró comú on teniu un sol fitxer de referència que necessita ser reutilitzat a través de múltiples mostres. La solució és convertir el canal de referència en un canal de valor que pot ser reutilitzat indefinidament.

#### Corregir el codi

Hi ha un parell de maneres d'abordar això depenent de quants fitxers estan afectats.

**Opció 1**: Teniu un sol fitxer de referència que esteu reutilitzant molt. Podeu simplement crear un tipus de canal de valor, que pot ser utilitzat una i altra vegada. Hi ha tres maneres de fer-ho:

**1a** Utilitzar `channel.value()`:

```groovy title="exhausted.nf (fixed - Option 1a)" hl_lines="2" linenums="21"
workflow {
    reference_ch = channel.value('baseline_reference')  // El canal de valor pot ser reutilitzat
    input_ch = channel.of('sample1', 'sample2', 'sample3')

    PROCESS_FILES(reference_ch, input_ch)
}
```

**1b** Utilitzar l'[operador](https://www.nextflow.io/docs/latest/reference/operator.html#first) `first()`:

```groovy title="exhausted.nf (fixed - Option 1b)" hl_lines="2" linenums="21"
workflow {
    reference_ch = channel.of('baseline_reference').first()  // Convertir a canal de valor
    input_ch = channel.of('sample1', 'sample2', 'sample3')

    PROCESS_FILES(reference_ch, input_ch)
}
```

**1c.** Utilitzar l'[operador](https://www.nextflow.io/docs/latest/reference/operator.html#collect) `collect()`:

```groovy title="exhausted.nf (fixed - Option 1c)" hl_lines="2" linenums="21"
workflow {
    reference_ch = channel.of('baseline_reference').collect()  // Convertir a canal de valor
    input_ch = channel.of('sample1', 'sample2', 'sample3')

    PROCESS_FILES(reference_ch, input_ch)
}
```

**Opció 2**: En escenaris més complexos, potser on teniu múltiples fitxers de referència per a totes les mostres al canal de mostres, podeu utilitzar l'operador `combine` per crear un nou canal que combini els dos canals en tuples:

```groovy title="exhausted.nf (fixed - Option 2)" hl_lines="4" linenums="21"
workflow {
    reference_ch = channel.of('baseline_reference','other_reference')
    input_ch = channel.of('sample1', 'sample2', 'sample3')
    combined_ch = reference_ch.combine(input_ch)  // Crea producte cartesià

    PROCESS_FILES(combined_ch)
}
```

L'operador `.combine()` genera un producte cartesià dels dos canals, així que cada element a `reference_ch` s'aparellarà amb cada element a `input_ch`. Això permet que el procés s'executi per a cada mostra mentre encara utilitza la referència.

Això requereix que l'entrada del procés sigui ajustada. Al nostre exemple, l'inici de la definició del procés hauria de ser ajustat de la següent manera:

```groovy title="exhausted.nf (fixed - Option 2)" hl_lines="5" linenums="1"
#!/usr/bin/env nextflow

process PROCESS_FILES {
    input:
        tuple val(reference), val(sample_name)
```

Aquest enfocament pot no ser adequat en totes les situacions.

#### Executar el pipeline

Proveu una de les correccions anteriors i executeu el workflow de nou:

```bash
nextflow run exhausted.nf
```

??? success "Sortida de la comanda"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `exhausted.nf` [maniac_leavitt] DSL2 - revision: f372a56a7d

    executor >  local (3)
    [80/0779e9] PROCESS_FILES (3) | 3 of 3 ✔
    ```

Ara hauríeu de veure totes tres mostres sent processades en lloc de només una.

### 2.3. Estructura de Contingut de Canal Incorrecta

Quan els workflows arriben a un cert nivell de complexitat, pot ser una mica difícil fer un seguiment de les estructures internes de cada canal, i la gent comunament genera discrepàncies entre el que el procés espera i el que el canal realment conté. Això és més subtil que el problema que vam discutir anteriorment, on el nombre de canals era incorrecte. En aquest cas, podeu tenir el nombre correcte de canals d'entrada, però l'estructura interna d'un o més d'aquests canals no coincideix amb el que el procés espera.

#### Executar el pipeline

```bash
nextflow run bad_channel_shape.nf
```

??? failure "Sortida de la comanda"

    ```console
    Launching `bad_channel_shape.nf` [hopeful_pare] DSL2 - revision: ffd66071a1

    executor >  local (3)
    executor >  local (3)
    [3f/c2dcb3] PROCESS_FILES (3) [  0%] 0 of 3 ✘
    ERROR ~ Error executing process > 'PROCESS_FILES (1)'

    Caused by:
      Missing output file(s) `[sample1, file1.txt]_output.txt` expected by process `PROCESS_FILES (1)`


    Command executed:

      echo "Processing [sample1, file1.txt]" > [sample1, file1.txt]_output.txt

    Command exit status:
      0

    Command output:
      (empty)

    Work dir:
      /workspaces/training/side-quests/debugging/work/d6/1fb69d1d93300bbc9d42f1875b981e

    Tip: when you have fixed the problem you can continue the execution adding the option `-resume` to the run command line

    -- Check '.nextflow.log' file for details
    ```

#### Comprovar el codi

Els claudàtors al missatge d'error proporcionen la pista aquí - el procés està tractant la tupla com un sol valor, cosa que no és el que volem. Examinem `bad_channel_shape.nf`:

```groovy title="bad_channel_shape.nf" hl_lines="5 20-22" linenums="1"
#!/usr/bin/env nextflow

process PROCESS_FILES {
    input:
        val sample_name  // Espera un sol valor, rep una tupla

    output:
        path "${sample_name}_output.txt"

    script:
    """
    echo "Processing ${sample_name}" > ${sample_name}_output.txt
    """
}

workflow {

    // El canal emet tuples, però el procés espera valors únics
    input_ch = channel.of(
      ['sample1', 'file1.txt'],
      ['sample2', 'file2.txt'],
      ['sample3', 'file3.txt']
    )
    PROCESS_FILES(input_ch)
}
```

Podeu veure que estem generant un canal compost de tuples: `['sample1', 'file1.txt']`, però el procés espera un sol valor, `val sample_name`. La comanda executada mostra que el procés està intentant crear un fitxer anomenat `[sample3, file3.txt]_output.txt`, que no és la sortida prevista.

#### Corregir el codi

Per corregir això, si el procés requereix ambdues entrades podríem ajustar el procés per acceptar una tupla:

=== "Opció 1: Acceptar tupla al procés"

    === "Després"

        ```groovy title="bad_channel_shape.nf" hl_lines="5"  linenums="1"
        #!/usr/bin/env nextflow

        process PROCESS_FILES {
            input:
                tuple val(sample_name), val(file_name)  // Corregit: Acceptar tupla

            output:
                path "${sample_name}_output.txt"

            script:
            """
            echo "Processing ${sample_name}" > ${sample_name}_output.txt
            """
        }

        workflow {

            // El canal emet tuples, però el procés espera valors únics
            input_ch = channel.of(
              ['sample1', 'file1.txt'],
              ['sample2', 'file2.txt'],
              ['sample3', 'file3.txt']
            )
            PROCESS_FILES(input_ch)
        }
        ```

    === "Abans"

        ```groovy title="bad_channel_shape.nf" hl_lines="5" linenums="1"
        #!/usr/bin/env nextflow

        process PROCESS_FILES {
            input:
                val sample_name  // Espera un sol valor, rep una tupla

            output:
                path "${sample_name}_output.txt"

            script:
            """
            echo "Processing ${sample_name}" > ${sample_name}_output.txt
            """
        }

        workflow {

            // El canal emet tuples, però el procés espera valors únics
            input_ch = channel.of(
              ['sample1', 'file1.txt'],
              ['sample2', 'file2.txt'],
              ['sample3', 'file3.txt']
            )
            PROCESS_FILES(input_ch)
        }
        ```

=== "Opció 2: Extreure el primer element"

    === "Després"

        ```groovy title="bad_channel_shape.nf" hl_lines="9" linenums="16"
        workflow {

            // El canal emet tuples, però el procés espera valors únics
            input_ch = channel.of(
              ['sample1', 'file1.txt'],
              ['sample2', 'file2.txt'],
              ['sample3', 'file3.txt']
            )
            PROCESS_FILES(input_ch.map { it[0] })  // Corregit: Extreure el primer element
        }
        ```

    === "Abans"

        ```groovy title="bad_channel_shape.nf" hl_lines="9" linenums="16"
        workflow {

            // El canal emet tuples, però el procés espera valors únics
            input_ch = channel.of(
              ['sample1', 'file1.txt'],
              ['sample2', 'file2.txt'],
              ['sample3', 'file3.txt']
            )
            PROCESS_FILES(input_ch)
        }
        ```

#### Executar el pipeline

Trieu una de les solucions i torneu a executar el workflow:

```bash
nextflow run bad_channel_shape.nf
```

??? success "Sortida de la comanda"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `bad_channel_shape.nf` [clever_thompson] DSL2 - revision: 8cbcae3746

    executor >  local (3)
    [bb/80a958] PROCESS_FILES (2) | 3 of 3 ✔
    ```

### 2.4. Tècniques de Depuració de Canals

#### Utilitzar `.view()` per a la Inspecció de Canals

L'eina de depuració més potent per a canals és l'operador `.view()`. Amb `.view()`, podeu entendre la forma dels vostres canals en totes les etapes per ajudar amb la depuració.

#### Executar el pipeline

Executeu `bad_channel_shape_viewed.nf` per veure això en acció:

```bash
nextflow run bad_channel_shape_viewed.nf
```

??? success "Sortida de la comanda"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `bad_channel_shape_viewed.nf` [maniac_poisson] DSL2 - revision: b4f24dc9da

    executor >  local (3)
    [c0/db76b3] PROCESS_FILES (3) [100%] 3 of 3 ✔
    Channel content: [sample1, file1.txt]
    Channel content: [sample2, file2.txt]
    Channel content: [sample3, file3.txt]
    After mapping: sample1
    After mapping: sample2
    After mapping: sample3
    ```

#### Comprovar el codi

Examinem `bad_channel_shape_viewed.nf` per veure com s'utilitza `.view()`:

```groovy title="bad_channel_shape_viewed.nf" linenums="16" hl_lines="9 11"
workflow {

    // El canal emet tuples, però el procés espera valors únics
    input_ch = channel.of(
      ['sample1', 'file1.txt'],
      ['sample2', 'file2.txt'],
      ['sample3', 'file3.txt']
    )
    .view { "Channel content: $it" }  // Debug: Mostrar contingut original del canal
    .map { tuple -> tuple[0] }        // Transform: Extreure el primer element
    .view { "After mapping: $it" }    // Debug: Mostrar contingut transformat del canal

    PROCESS_FILES(input_ch)
}
```

#### Corregir el codi

Per estalviar-vos d'utilitzar operacions `.view()` excessivament en el futur per entendre el contingut del canal, és aconsellable afegir alguns comentaris per ajudar:

```groovy title="bad_channel_shape_viewed.nf (with comments)" linenums="16" hl_lines="8 9"
workflow {

    // El canal emet tuples, però el procés espera valors únics
    input_ch = channel.of(
            ['sample1', 'file1.txt'],
            ['sample2', 'file2.txt'],
            ['sample3', 'file3.txt'],
        ) // [sample_name, file_name]
        .map { tuple -> tuple[0] } // sample_name

    PROCESS_FILES(input_ch)
}
```

Això es tornarà més important a mesura que els vostres workflows creixin en complexitat i l'estructura del canal es torni més opaca.

#### Executar el pipeline

```bash
nextflow run bad_channel_shape_viewed.nf
```

??? success "Sortida de la comanda"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `bad_channel_shape_viewed.nf` [marvelous_koch] DSL2 - revision: 03e79cdbad

    executor >  local (3)
    [ff/d67cec] PROCESS_FILES (2) | 3 of 3 ✔
    Channel content: [sample1, file1.txt]
    Channel content: [sample2, file2.txt]
    Channel content: [sample3, file3.txt]
    After mapping: sample1
    After mapping: sample2
    After mapping: sample3
    ```

### Conclusió

Molts errors d'estructura de canals es poden crear amb sintaxi Nextflow vàlida. Podeu depurar errors d'estructura de canals entenent el flux de dades, utilitzant operadors `.view()` per a la inspecció, i reconeixent patrons de missatges d'error com claudàtors indicant estructures de tupla inesperades.

### Què segueix?

Apreneu sobre errors creats per definicions de processos.

---

## 3. Errors d'Estructura de Processos

La majoria dels errors que trobareu relacionats amb processos estaran relacionats amb errors que heu comès en formar la comanda, o amb problemes relacionats amb el programari subjacent. Dit això, de manera similar als problemes de canals anteriors, podeu cometre errors a la definició del procés que no qualifiquen com a errors de sintaxi, però que causaran errors en temps d'execució.

### 3.1. Fitxers de Sortida que Falten

Un error comú quan s'escriuen processos és fer alguna cosa que genera una discrepància entre el que el procés espera i el que es genera.

#### Executar el pipeline

```bash
nextflow run missing_output.nf
```

??? failure "Sortida de la comanda"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `missing_output.nf` [zen_stone] DSL2 - revision: 37ff61f926

    executor >  local (3)
    executor >  local (3)
    [fd/2642e9] process > PROCESS_FILES (2) [ 66%] 2 of 3, failed: 2
    ERROR ~ Error executing process > 'PROCESS_FILES (3)'

    Caused by:
      Missing output file(s) `sample3.txt` expected by process `PROCESS_FILES (3)`


    Command executed:

      echo "Processing sample3" > sample3_output.txt

    Command exit status:
      0

    Command output:
      (empty)

    Work dir:
      /workspaces/training/side-quests/debugging/work/02/9604d49fb8200a74d737c72a6c98ed

    Tip: when you have fixed the problem you can continue the execution adding the option `-resume` to the run command line

    -- Check '.nextflow.log' file for details
    ```

#### Comprovar el codi

El missatge d'error indica que el procés esperava produir un fitxer de sortida anomenat `sample3.txt`, però l'script realment crea `sample3_output.txt`. Examinem la definició del procés a `missing_output.nf`:

```groovy title="missing_output.nf" linenums="3" hl_lines="6 10"
process PROCESS_FILES {
    input:
    val sample_name

    output:
    path "${sample_name}.txt"  // Espera: sample3.txt

    script:
    """
    echo "Processing ${sample_name}" > ${sample_name}_output.txt  // Crea: sample3_output.txt
    """
}
```

Hauríeu de veure que hi ha una discrepància entre el nom del fitxer de sortida al bloc `output:`, i el que s'utilitza a l'script. Aquesta discrepància fa que el procés falli. Si trobeu aquest tipus d'error, torneu enrere i comproveu que les sortides coincideixin entre la vostra definició de procés i el vostre bloc de sortida.

Si el problema encara no està clar, comproveu el directori de treball mateix per identificar els fitxers de sortida reals creats:

```bash
❯ ls -h work/02/9604d49fb8200a74d737c72a6c98ed
sample3_output.txt
```

Per a aquest exemple això ens destacaria que s'està incorporant un sufix `_output` al nom del fitxer de sortida, contrari a la nostra definició `output:`.

#### Corregir el codi

Corregiu la discrepància fent que el nom del fitxer de sortida sigui consistent:

=== "Després"

    ```groovy title="missing_output.nf" hl_lines="6 10" linenums="3"
    process PROCESS_FILES {
        input:
        val sample_name

        output:
        path "${sample_name}_output.txt"  // Corregit: Coincideix amb la sortida de l'script

        script:
        """
        echo "Processing ${sample_name}" > ${sample_name}_output.txt
        """
    }
    ```

=== "Abans"

    ```groovy title="missing_output.nf" hl_lines="6 10" linenums="3"
    process PROCESS_FILES {
        input:
        val sample_name

        output:
        path "${sample_name}.txt"  // Espera: sample3.txt

        script:
        """
        echo "Processing ${sample_name}" > ${sample_name}_output.txt  // Crea: sample3_output.txt
        """
    }
    ```

#### Executar el pipeline

```bash
nextflow run missing_output.nf
```

??? success "Sortida de la comanda"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `missing_output.nf` [elated_hamilton] DSL2 - revision: 961938ee2b

    executor >  local (3)
    [16/1c437c] PROCESS_FILES (3) | 3 of 3 ✔
    ```

### 3.2. Programari que falta

Una altra classe d'errors es produeix a causa d'errors en l'aprovisionament de programari. `missing_software.nf` és un workflow sintàcticament vàlid, però depèn d'algun programari extern per proporcionar la comanda `cowpy` que utilitza.

#### Executar el pipeline

```bash
nextflow run missing_software.nf
```

??? failure "Sortida de la comanda"

    ```console hl_lines="12 18"
    ERROR ~ Error executing process > 'PROCESS_FILES (3)'

    Caused by:
      Process `PROCESS_FILES (3)` terminated with an error exit status (127)


    Command executed:

      cowpy sample3 > sample3_output.txt

    Command exit status:
      127

    Command output:
      (empty)

    Command error:
      .command.sh: line 2: cowpy: command not found

    Work dir:
      /workspaces/training/side-quests/debugging/work/82/42a5bfb60c9c6ee63ebdbc2d51aa6e

    Tip: you can try to figure out what's wrong by changing to the process work directory and showing the script file named `.command.sh`

    -- Check '.nextflow.log' file for details
    ```

El procés no té accés a la comanda que estem especificant. De vegades això és perquè un script està present al directori `bin` del workflow, però no s'ha fet executable. Altres vegades és perquè el programari no està instal·lat al contenidor o entorn on s'està executant el workflow.

#### Comprovar el codi

Vigileu amb aquest codi de sortida `127` - us diu exactament el problema. Examinem `missing_software.nf`:

```groovy title="missing_software.nf" linenums="3" hl_lines="3"
process PROCESS_FILES {

    container 'community.wave.seqera.io/library/cowpy:1.1.5--3db457ae1977a273'

    input:
    val sample_name

    output:
    path "${sample_name}_output.txt"

    script:
    """
    cowpy ${sample_name} > ${sample_name}_output.txt
    """
}
```

#### Corregir el codi

Hem estat una mica poc sincers aquí, i en realitat no hi ha res malament amb el codi. Només necessitem especificar la configuració necessària per executar el procés de manera que tingui accés a la comanda en qüestió. En aquest cas el procés té una definició de contenidor, així que tot el que necessitem fer és executar el workflow amb Docker habilitat.

#### Executar el pipeline

Hem configurat un perfil Docker per a vosaltres a `nextflow.config`, així que podeu executar el workflow amb:

```bash
nextflow run missing_software.nf -profile docker
```

??? success "Sortida de la comanda"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `missing_software.nf` [awesome_stonebraker] DSL2 - revision: 0296d12839

    executor >  local (3)
    [38/ab20d1] PROCESS_FILES (1) | 3 of 3 ✔
    ```

!!! note "Nota"

    Per aprendre més sobre com Nextflow utilitza contenidors, vegeu [Hello Nextflow](../hello_nextflow/05_hello_containers.md)

### 3.3. Mala configuració de recursos

En ús de producció, estareu configurant recursos als vostres processos. Per exemple `memory` defineix la quantitat màxima de memòria disponible per al vostre procés, i si el procés excedeix això, el vostre planificador normalment matarà el procés i retornarà un codi de sortida de `137`. No podem demostrar això aquí perquè estem utilitzant l'executor `local`, però podem mostrar alguna cosa similar amb `time`.

#### Executar el pipeline

`bad_resources.nf` té configuració de procés amb un límit poc realista de temps d'1 mil·lisegon:

```bash
nextflow run bad_resources.nf -profile docker
```

??? failure "Sortida de la comanda"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `bad_resources.nf` [disturbed_elion] DSL2 - revision: 27d2066e86

    executor >  local (3)
    [c0/ded8e1] PROCESS_FILES (3) | 0 of 3 ✘
    ERROR ~ Error executing process > 'PROCESS_FILES (2)'

    Caused by:
      Process exceeded running time limit (1ms)

    Command executed:

      cowpy sample2 > sample2_output.txt

    Command exit status:
      -

    Command output:
      (empty)

    Work dir:
      /workspaces/training/side-quests/debugging/work/53/f0a4cc56d6b3dc2a6754ff326f1349

    Container:
      community.wave.seqera.io/library/cowpy:1.1.5--3db457ae1977a273

    Tip: you can replicate the issue by changing to the process work dir and entering the command `bash .command.run`

     -- Check '.nextflow.log' file for details
    ```

#### Comprovar el codi

Examinem `bad_resources.nf`:

```groovy title="bad_resources.nf" linenums="3" hl_lines="3"
process PROCESS_FILES {

    time '1 ms'  // ERROR: Límit de temps poc realista

    input:
    val sample_name

    output:
    path "${sample_name}_output.txt"

    script:
    """
    sleep 1  // Triga 1 segon, però el límit de temps és 1ms
    cowpy ${sample_name} > ${sample_name}_output.txt
    """
}
```

Sabem que el procés trigarà més d'un segon (hem afegit un sleep allà per assegurar-nos), però el procés està configurat per esgotar el temps després d'1 mil·lisegon. Algú ha estat una mica poc realista amb la seva configuració!

#### Corregir el codi

Augmenteu el límit de temps a un valor realista:

=== "Després"

    ```groovy title="bad_resources.nf" hl_lines="3" linenums="3"
    process PROCESS_FILES {

        time '100 s'  // Corregit: Límit de temps realista

        input:
        val sample_name

        output:
        path "${sample_name}_output.txt"

        script:
        """
        sleep 1
        cowpy ${sample_name} > ${sample_name}_output.txt
        """
    }
    ```

=== "Abans"

    ```groovy title="bad_resources.nf" hl_lines="3" linenums="3"
    process PROCESS_FILES {

        time '1 ms'  // ERROR: Límit de temps poc realista

        input:
        val sample_name

        output:
        path "${sample_name}_output.txt"

        script:
        """
        sleep 1  // Triga 1 segon, però el límit de temps és 1ms
        cowpy ${sample_name} > ${sample_name}_output.txt
        """
    }
    ```

#### Executar el pipeline

```bash
nextflow run bad_resources.nf -profile docker
```

??? success "Sortida de la comanda"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `bad_resources.nf` [friendly_mcclintock] DSL2 - revision: 381567d2c1

    executor >  local (3)
    [c2/9b4c41] PROCESS_FILES (3) | 3 of 3 ✔
    ```

Si us assegureu de llegir els vostres missatges d'error, fallades com aquesta no us haurien de desconcertar durant massa temps. Però assegureu-vos d'entendre els requisits de recursos de les comandes que esteu executant perquè pugueu configurar les vostres directives de recursos adequadament.

### 3.4. Tècniques de Depuració de Processos

Quan els processos fallen o es comporten de manera inesperada, necessiteu tècniques sistemàtiques per investigar què va anar malament. El directori de treball conté tota la informació que necessiteu per depurar l'execució del procés.

#### Utilitzar la Inspecció del Directori de Treball

L'eina de depuració més potent per a processos és examinar el directori de treball. Quan un procés falla, Nextflow crea un directori de treball per a aquesta execució específica del procés que conté tots els fitxers necessaris per entendre què va passar.

#### Executar el pipeline

Utilitzem l'exemple `missing_output.nf` d'abans per demostrar la inspecció del directori de treball (regenereu una discrepància de denominació de sortida si cal):

```bash
nextflow run missing_output.nf
```

??? failure "Sortida de la comanda"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `missing_output.nf` [irreverent_payne] DSL2 - revision: 3d5117f7e2

    executor >  local (3)
    [5d/d544a4] PROCESS_FILES (2) | 0 of 3 ✘
    ERROR ~ Error executing process > 'PROCESS_FILES (1)'

    Caused by:
      Missing output file(s) `sample1.txt` expected by process `PROCESS_FILES (1)`

    Command executed:

      echo "Processing sample1" > sample1_output.txt

    Command exit status:
      0

    Command output:
      (empty)

    Work dir:
      /workspaces/training/side-quests/debugging/work/1e/2011154d0b0f001cd383d7364b5244

    Tip: you can replicate the issue by changing to the process work dir and entering the command `bash .command.run`

     -- Check '.nextflow.log' file for details
    ```

#### Comprovar el directori de treball

Quan obteniu aquest error, el directori de treball conté tota la informació de depuració. Trobeu el camí del directori de treball del missatge d'error i examineu els seus continguts:

```bash
# Trobar el directori de treball del missatge d'error
ls work/02/9604d49fb8200a74d737c72a6c98ed/
```

Llavors podeu examinar els fitxers clau:

##### Comprovar l'Script de Comanda

El fitxer `.command.sh` mostra exactament quina comanda es va executar:

```bash
# Veure la comanda executada
cat work/02/9604d49fb8200a74d737c72a6c98ed/.command.sh
```

Això revela:

- **Substitució de variables**: Si les variables de Nextflow es van expandir correctament
- **Camins de fitxers**: Si els fitxers d'entrada es van localitzar correctament
- **Estructura de comanda**: Si la sintaxi de l'script és correcta

Problemes comuns a buscar:

- **Cometes que falten**: Variables que contenen espais necessiten cometes adequades
- **Camins de fitxers incorrectes**: Fitxers d'entrada que no existeixen o estan en ubicacions incorrectes
- **Noms de variables incorrectes**: Errors tipogràfics en referències de variables
- **Configuració d'entorn que falta**: Comandes que depenen d'entorns específics

##### Comprovar la Sortida d'Error

El fitxer `.command.err` conté els missatges d'error reals:

```bash
# Veure la sortida d'error
cat work/02/9604d49fb8200a74d737c72a6c98ed/.command.err
```

Aquest fitxer mostrarà:

- **Codis de sortida**: 127 (comanda no trobada), 137 (matat), etc.
- **Errors de permisos**: Problemes d'accés a fitxers
- **Errors de programari**: Missatges d'error específics de l'aplicació
- **Errors de recursos**: Límit de memòria/temps excedit

##### Comprovar la Sortida Estàndard

El fitxer `.command.out` mostra què va produir la vostra comanda:

```bash
# Veure la sortida estàndard
cat work/02/9604d49fb8200a74d737c72a6c98ed/.command.out
```

Això ajuda a verificar:

- **Sortida esperada**: Si la comanda va produir els resultats correctes
- **Execució parcial**: Si la comanda va començar però va fallar a mig camí
- **Informació de depuració**: Qualsevol sortida de diagnòstic del vostre script

##### Comprovar el Codi de Sortida

El fitxer `.exitcode` conté el codi de sortida del procés:

```bash
# Veure el codi de sortida
cat work/*/*/.exitcode
```

Codis de sortida comuns i els seus significats:

- **Codi de sortida 127**: Comanda no trobada - comprovar la instal·lació del programari
- **Codi de sortida 137**: Procés matat - comprovar límits de memòria/temps

##### Comprovar l'Existència de Fitxers

Quan els processos fallen a causa de fitxers de sortida que falten, comproveu quins fitxers es van crear realment:

```bash
# Llistar tots els fitxers al directori de treball
ls -la work/02/9604d49fb8200a74d737c72a6c98ed/
```

Això ajuda a identificar:

- **Discrepàncies de denominació de fitxers**: Fitxers de sortida amb noms diferents dels esperats
- **Problemes de permisos**: Fitxers que no es van poder crear
- **Problemes de camí**: Fitxers creats en directoris incorrectes

Al nostre exemple anterior, això ens va confirmar que mentre el nostre `sample3.txt` esperat no estava present, `sample3_output.txt` sí que ho estava:

```bash
❯ ls -h work/02/9604d49fb8200a74d737c72a6c98ed
sample3_output.txt
```

### Conclusió

La depuració de processos requereix examinar directoris de treball per entendre què va anar malament. Els fitxers clau inclouen `.command.sh` (l'script executat), `.command.err` (missatges d'error), i `.command.out` (sortida estàndard). Codis de sortida com 127 (comanda no trobada) i 137 (procés matat) proporcionen pistes de diagnòstic immediates sobre el tipus de fallada.

### Què segueix?

Apreneu sobre les eines de depuració integrades de Nextflow i enfocaments sistemàtics per a la resolució de problemes.

---

## 4. Eines de Depuració Integrades i Tècniques Avançades

Nextflow proporciona diverses eines integrades potents per a la depuració i l'anàlisi de l'execució del workflow. Aquestes eines us ajuden a entendre què va anar malament, on va anar malament, i com solucionar-ho eficientment.

### 4.1. Sortida de Procés en Temps Real

De vegades necessiteu veure què està passant dins dels processos en execució. Podeu habilitar la sortida de procés en temps real, que us mostra exactament què està fent cada tasca mentre s'executa.

#### Executar el pipeline

`bad_channel_shape_viewed.nf` dels nostres exemples anteriors va imprimir contingut de canal utilitzant `.view()`, però també podem utilitzar la directiva `debug` per fer eco de variables des de dins del procés mateix, cosa que demostrem a `bad_channel_shape_viewed_debug.nf`. Executeu el workflow:

```bash
nextflow run bad_channel_shape_viewed_debug.nf
```

??? success "Sortida de la comanda"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `bad_channel_shape_viewed_debug.nf` [agitated_crick] DSL2 - revision: ea3676d9ec

    executor >  local (3)
    [c6/2dac51] process > PROCESS_FILES (3) [100%] 3 of 3 ✔
    Channel content: [sample1, file1.txt]
    Channel content: [sample2, file2.txt]
    Channel content: [sample3, file3.txt]
    After mapping: sample1
    After mapping: sample2
    After mapping: sample3
    Sample name inside process is sample2

    Sample name inside process is sample1

    Sample name inside process is sample3
    ```

#### Comprovar el codi

Examinem `bad_channel_shape_viewed_debug.nf` per veure com funciona la directiva `debug`:

```groovy title="bad_channel_shape_viewed_debug.nf" linenums="3" hl_lines="2"
process PROCESS_FILES {
    debug true  // Habilitar sortida en temps real

    input:
    val sample_name

    output:
    path "${sample_name}_output.txt"

    script:
    """
    echo "Sample name inside process is ${sample_name}"
    echo "Processing ${sample_name}" > ${sample_name}_output.txt
    """
}
```

La directiva `debug` pot ser una manera ràpida i convenient d'entendre l'entorn d'un procés.

### 4.2. Mode de Previsualització

De vegades voleu detectar problemes abans que s'executin els processos. Nextflow proporciona una bandera per a aquest tipus de depuració proactiva: `-preview`.

#### Executar el pipeline

El mode de previsualització us permet provar la lògica del workflow sense executar comandes. Això pot ser força útil per comprovar ràpidament l'estructura del vostre workflow i assegurar-vos que els processos estan connectats correctament sense executar cap comanda real.

!!! note "Nota"

    Si vau corregir `bad_syntax.nf` anteriorment, reintroduïu l'error de sintaxi eliminant la clau de tancament després del bloc script abans d'executar aquesta comanda.

Executeu aquesta comanda:

```bash
nextflow run bad_syntax.nf -preview
```

??? failure "Sortida de la comanda"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `bad_syntax.nf` [magical_mercator] DSL2 - revision: 550b9a8873

    Error bad_syntax.nf:24:1: Unexpected input: '<EOF>'

    ERROR ~ Script compilation failed

     -- Check '.nextflow.log' file for details
    ```

El mode de previsualització és particularment útil per detectar errors de sintaxi aviat sense executar cap procés. Valida l'estructura del workflow i les connexions de processos abans de l'execució.

### 4.3. Execució Stub per a Proves de Lògica

De vegades els errors són difícils de depurar perquè les comandes triguen massa temps, requereixen programari especial, o fallen per raons complexes. L'execució stub us permet provar la lògica del workflow sense executar les comandes reals.

#### Executar el pipeline

Quan esteu desenvolupant un procés Nextflow, podeu utilitzar la directiva `stub` per definir comandes 'fictícies' que generen sortides de la forma correcta sense executar la comanda real. Aquest enfocament és particularment valuós quan voleu verificar que la vostra lògica de workflow és correcta abans de tractar amb les complexitats del programari real.

Per exemple, recordeu el nostre `missing_software.nf` d'abans? El que tenia programari que faltava que impedia que el workflow s'executés fins que vam afegir `-profile docker`? `missing_software_with_stub.nf` és un workflow molt similar. Si l'executem de la mateixa manera, generarem el mateix error:

```bash
nextflow run missing_software_with_stub.nf
```

??? failure "Sortida de la comanda"

    ```console hl_lines="12 18"
    ERROR ~ Error executing process > 'PROCESS_FILES (3)'

    Caused by:
      Process `PROCESS_FILES (3)` terminated with an error exit status (127)


    Command executed:

      cowpy sample3 > sample3_output.txt

    Command exit status:
      127

    Command output:
      (empty)

    Command error:
      .command.sh: line 2: cowpy: command not found

    Work dir:
      /workspaces/training/side-quests/debugging/work/82/42a5bfb60c9c6ee63ebdbc2d51aa6e

    Tip: you can try to figure out what's wrong by changing to the process work directory and showing the script file named `.command.sh`

    -- Check '.nextflow.log' file for details
    ```

No obstant això, aquest workflow no produirà errors si l'executem amb `-stub-run`, fins i tot sense el perfil `docker`:

```bash
nextflow run missing_software_with_stub.nf -stub-run
```

??? success "Sortida de la comanda"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `missing_software_with_stub.nf` [astonishing_shockley] DSL2 - revision: f1f4f05d7d

    executor >  local (3)
    [b5/2517a3] PROCESS_FILES (3) | 3 of 3 ✔
    ```

#### Comprovar el codi

Examinem `missing_software_with_stub.nf`:

```groovy title="missing_software.nf (with stub)" hl_lines="16-19" linenums="3"
process PROCESS_FILES {

    container 'community.wave.seqera.io/library/cowpy:1.1.5--3db457ae1977a273'

    input:
    val sample_name

    output:
    path "${sample_name}_output.txt"

    script:
    """
    cowpy ${sample_name} > ${sample_name}_output.txt
    """

    stub:
    """
    touch ${sample_name}_output.txt
    """
}
```

Relatiu a `missing_software.nf`, aquest procés té una directiva `stub:` especificant una comanda a utilitzar en lloc de la especificada a `script:`, en el cas que Nextflow s'executi en mode stub.

La comanda `touch` que estem utilitzant aquí no depèn de cap programari o entrades apropiades, i s'executarà en totes les situacions, permetent-nos depurar la lògica del workflow sense preocupar-nos pels interns del procés.

**L'execució stub ajuda a depurar:**

- Estructura de canals i flux de dades
- Connexions de processos i dependències
- Propagació de paràmetres
- Lògica de workflow sense dependències de programari

### 4.4. Enfocament Sistemàtic de Depuració

Ara que heu après tècniques individuals de depuració - des de fitxers de traça i directoris de treball fins al mode de previsualització, execució stub i monitorització de recursos - liguem-les juntes en una metodologia sistemàtica. Tenir un enfocament estructurat us impedeix sentir-vos aclaparats per errors complexos i assegura que no us perdeu pistes importants.

Aquesta metodologia combina totes les eines que hem cobert en un workflow eficient:

**Mètode de Depuració de Quatre Fases:**

**Fase 1: Resolució d'Errors de Sintaxi (5 minuts)**

1. Comprovar subratllats vermells a VSCode o al vostre IDE
2. Executar `nextflow run workflow.nf -preview` per identificar problemes de sintaxi
3. Corregir tots els errors de sintaxi (claus que falten, comes finals, etc.)
4. Assegurar-se que el workflow s'analitza correctament abans de continuar

**Fase 2: Avaluació Ràpida (5 minuts)**

1. Llegir els missatges d'error d'execució amb cura
2. Comprovar si és un error d'execució, lògica o recursos
3. Utilitzar el mode de previsualització per provar la lògica bàsica del workflow

**Fase 3: Investigació Detallada (15-30 minuts)**

1. Trobar el directori de treball de la tasca fallida
2. Examinar fitxers de registre
3. Afegir operadors `.view()` per inspeccionar canals
4. Utilitzar `-stub-run` per provar la lògica del workflow sense execució

**Fase 4: Corregir i Validar (15 minuts)**

1. Fer correccions mínimes dirigides
2. Provar amb resume: `nextflow run workflow.nf -resume`
3. Verificar l'execució completa del workflow

!!! tip "Utilitzar Resume per a una Depuració Eficient"

    Un cop heu identificat un problema, necessiteu una manera eficient de provar les vostres correccions sense perdre temps reexecutant parts reeixides del vostre workflow. La funcionalitat `-resume` de Nextflow és invaluable per a la depuració.

    Haureu trobat `-resume` si heu treballat amb [Hello Nextflow](../hello_nextflow/), i és important que en feu bon ús quan depureu per estalviar-vos esperar mentre s'executen els processos abans del vostre procés problemàtic.

    **Estratègia de depuració amb resume:**

    1. Executar el workflow fins a la fallada
    2. Examinar el directori de treball per a la tasca fallida
    3. Corregir el problema específic
    4. Reprendre per provar només la correcció
    5. Repetir fins que el workflow es completi

#### Perfil de Configuració de Depuració

Per fer aquest enfocament sistemàtic encara més eficient, podeu crear una configuració de depuració dedicada que habiliti automàticament totes les eines que necessiteu:

```groovy title="nextflow.config (debug profile)" linenums="1"
profiles {
    debug {
        process {
            debug = true
            cleanup = false

            // Recursos conservadors per a depuració
            maxForks = 1
            memory = '2.GB'
            cpus = 1
        }
    }
}
```

Llavors podeu executar el pipeline amb aquest perfil habilitat:

```bash
nextflow run workflow.nf -profile debug
```

Aquest perfil habilita la sortida en temps real, preserva els directoris de treball, i limita la paral·lelització per a una depuració més fàcil.

### 4.5. Exercici Pràctic de Depuració

Ara és el moment de posar en pràctica l'enfocament sistemàtic de depuració. El workflow `buggy_workflow.nf` conté diversos errors comuns que representen els tipus de problemes que trobareu en el desenvolupament real.

!!! exercise "Exercici"

    Utilitzeu l'enfocament sistemàtic de depuració per identificar i corregir tots els errors a `buggy_workflow.nf`. Aquest workflow intenta processar dades de mostra d'un fitxer CSV però conté múltiples errors intencionats que representen escenaris de depuració comuns.

    Comenceu executant el workflow per veure el primer error:

    ```bash
    nextflow run buggy_workflow.nf
    ```

    ??? failure "Sortida de la comanda"

        ```console
        N E X T F L O W   ~  version 25.10.2

        Launching `buggy_workflow.nf` [wise_ramanujan] DSL2 - revision: d51a8e83fd

        ERROR ~ Range [11, 12) out of bounds for length 11

         -- Check '.nextflow.log' file for details
        ```

        Aquest error críptic indica un problema d'anàlisi al voltant de la línia 11-12 al bloc `params{}`. L'analitzador v2 detecta problemes estructurals aviat.

    Apliqueu el mètode de depuració de quatre fases que heu après:

    **Fase 1: Resolució d'Errors de Sintaxi**
    - Comprovar subratllats vermells a VSCode o al vostre IDE
    - Executar `nextflow run workflow.nf -preview` per identificar problemes de sintaxi
    - Corregir tots els errors de sintaxi (claus que falten, comes finals, etc.)
    - Assegurar-se que el workflow s'analitza correctament abans de continuar

    **Fase 2: Avaluació Ràpida**
    - Llegir els missatges d'error d'execució amb cura
    - Identificar si els errors són d'execució, lògica o relacionats amb recursos
    - Utilitzar el mode `-preview` per provar la lògica bàsica del workflow

    **Fase 3: Investigació Detallada**
    - Examinar directoris de treball per a tasques fallides
    - Afegir operadors `.view()` per inspeccionar canals
    - Comprovar fitxers de registre als directoris de treball
    - Utilitzar `-stub-run` per provar la lògica del workflow sense execució

    **Fase 4: Corregir i Validar**
    - Fer correccions dirigides
    - Utilitzar `-resume` per provar correccions eficientment
    - Verificar l'execució completa del workflow

    **Eines de Depuració a la Vostra Disposició:**
    ```bash
    # Mode de previsualització per a comprovació de sintaxi
    nextflow run buggy_workflow.nf -preview

    # Perfil de depuració per a sortida detallada
    nextflow run buggy_workflow.nf -profile debug

    # Execució stub per a proves de lògica
    nextflow run buggy_workflow.nf -stub-run

    # Resume després de correccions
    nextflow run buggy_workflow.nf -resume
    ```

    ??? solution "Solució"
        El `buggy_workflow.nf` conté 9 o 10 errors distints (depenent de com es comptin) cobrint totes les categories principals de depuració. Aquí hi ha un desglossament sistemàtic de cada error i com corregir-lo

        Comencem amb aquests errors de sintaxi:

        **Error 1: Error de Sintaxi - Coma Final**
        ```groovy linenums="21"
        output:
            path "${sample_id}_result.txt",  // ERROR: Coma final
        ```
        **Correcció:** Eliminar la coma final
        ```groovy linenums="21"
        output:
            path "${sample_id}_result.txt"
        ```

        **Error 2: Error de Sintaxi - Clau de Tancament que Falta**
        ```groovy linenums="24"
        script:
        """
        echo "Processing: ${sample}"
        cat ${input_file} > ${sample}_result.txt
        """
        // ERROR: Clau de tancament que falta per al procés processFiles
        ```
        **Correcció:** Afegir la clau de tancament que falta
        ```groovy linenums="29"
        """
        echo "Processing: ${sample_id}"
        cat ${input_file} > ${sample_id}_result.txt
        """
        }  // Afegir clau de tancament que falta
        ```

        **Error 3: Error de Nom de Variable**
        ```groovy linenums="26"
        echo "Processing: ${sample}"     // ERROR: hauria de ser sample_id
        cat ${input_file} > ${sample}_result.txt  // ERROR: hauria de ser sample_id
        ```
        **Correcció:** Utilitzar el nom de variable d'entrada correcte
        ```groovy linenums="26"
        echo "Processing: ${sample_id}"
        cat ${input_file} > ${sample_id}_result.txt
        ```

        **Error 4: Error de Variable No Definida**
        ```groovy linenums="87"
        heavy_ch = heavyProcess(sample_ids)  // ERROR: sample_ids no definit
        ```
        **Correcció:** Utilitzar el canal correcte i extreure IDs de mostra
        ```groovy linenums="87"
        heavy_ch = heavyProcess(input_ch)
        ```

        En aquest punt el workflow s'executarà, però encara obtindrem errors (p. ex. `Path value cannot be null` a `processFiles`), causats per una mala estructura de canal.

        **Error 5: Error d'Estructura de Canal - Sortida de Map Incorrecta**
        ```groovy linenums="83"
        .map { row -> row.sample_id }  // ERROR: processFiles espera tupla
        ```
        **Correcció:** Retornar l'estructura de tupla que processFiles espera
        ```groovy linenums="83"
        .map { row -> [row.sample_id, file(row.fastq_path)] }
        ```

        Però això trencarà la nostra crida per executar `heavyProcess()` anterior, així que necessitarem utilitzar un map per passar només els IDs de mostra a aquest procés:

        **Error 6: Mala estructura de canal per a heavyProcess**
        ```groovy linenums="87"
        heavy_ch = heavyProcess(input_ch)  // ERROR: input_ch ara té 2 elements per emissió- heavyProcess només necessita 1 (el primer)
        ```
        **Correcció:** Utilitzar el canal correcte i extreure IDs de mostra
        ```groovy linenums="87"
        heavy_ch = heavyProcess(input_ch.map{it[0]})
        ```

        Ara arribem una mica més lluny però rebem un error sobre `No such variable: i`, perquè no vam escapar una variable Bash.

        **Error 7: Error d'Escapament de Variable Bash**
        ```groovy linenums="48"
        echo "Heavy computation $i for ${sample_id}"  // ERROR: $i no escapat
        ```
        **Correcció:** Escapar la variable bash
        ```groovy linenums="48"
        echo "Heavy computation \${i} for ${sample_id}"
        ```

        Ara obtenim `Process exceeded running time limit (1ms)`, així que corregim el límit de temps d'execució per al procés rellevant:

        **Error 8: Error de Configuració de Recursos**
        ```groovy linenums="36"
        time '1 ms'  // ERROR: Límit de temps poc realista
        ```
        **Correcció:** Augmentar a un límit de temps realista
        ```groovy linenums="36"
        time '100 s'
        ```

        A continuació tenim un error `Missing output file(s)` per resoldre:

        **Error 9: Discrepància de Nom de Fitxer de Sortida**
        ```groovy linenums="49"
        done > ${sample_id}.txt  // ERROR: Nom de fitxer incorrecte, hauria de coincidir amb la declaració de sortida
        ```
        **Correcció:** Coincidir amb la declaració de sortida
        ```groovy linenums="49"
        done > ${sample_id}_heavy.txt
        ```

        Els dos primers processos van executar-se, però no el tercer.

        **Error 10: Discrepància de Nom de Fitxer de Sortida**
        ```groovy linenums="88"
        file_ch = channel.fromPath("*.txt") // Error: intentant prendre entrada del pwd en lloc d'un procés
        handleFiles(file_ch)
        ```
        **Correcció:** Prendre la sortida del procés anterior
        ```groovy linenums="88"
        handleFiles(heavyProcess.out)
        ```

        Amb això, tot el workflow hauria d'executar-se.

        **Workflow Complet Corregit:**
        ```groovy linenums="1"
        #!/usr/bin/env nextflow

        /*
        * Workflow amb errors per a exercicis de depuració
        * Aquest workflow conté diversos errors intencionats amb finalitats d'aprenentatge
        */

        params{
            // Paràmetres amb validació que falta
            input: Path = 'data/sample_data.csv'
            output: String = 'results'
        }

        /*
        * Procés amb discrepància entrada/sortida
        */
        process processFiles {
            publishDir "${params.output}/processed", mode: 'copy'

            input:
                tuple val(sample_id), path(input_file)

            output:
                path "${sample_id}_result.txt"

            script:
            """
            echo "Processing: ${sample_id}"
            cat ${input_file} > ${sample_id}_result.txt
            """
        }

        /*
        * Procés amb problemes de recursos
        */
        process heavyProcess {
            publishDir "${params.output}/heavy", mode: 'copy'

            time '100 s'

            input:
                val sample_id

            output:
                path "${sample_id}_heavy.txt"

            script:
            """
            # Simular computació pesada
            for i in {1..1000000}; do
                echo "Heavy computation \$i for ${sample_id}"
            done > ${sample_id}_heavy.txt
            """
        }

        /*
        * Procés amb problemes de gestió de fitxers
        */
        process handleFiles {
            publishDir "${params.output}/files", mode: 'copy'

            input:
                path input_file

            output:
                path "processed_${input_file}"

            script:
            """
            if [ -f "${input_file}" ]; then
                cp ${input_file} processed_${input_file}
            fi
            """
        }

        /*
        * Workflow principal amb problemes de canals
        */
        workflow {

            // Canal amb ús incorrecte
            input_ch = channel
                .fromPath(params.input)
                .splitCsv(header: true)
                .map { row -> [row.sample_id, file(row.fastq_path)] }

            processed_ch = processFiles(input_ch)

            heavy_ch = heavyProcess(input_ch.map{it[0]})

            handleFiles(heavyProcess.out)
        }
        ```

**Categories d'Errors Cobertes:**

- **Errors de sintaxi**: Claus que falten, comes finals, variables no definides
- **Errors d'estructura de canals**: Formes de dades incorrectes, canals no definits
- **Errors de processos**: Discrepàncies de fitxers de sortida, escapament de variables
- **Errors de recursos**: Límits de temps poc realistes

**Lliçons Clau de Depuració:**

1. **Llegir els missatges d'error amb cura** - sovint apunten directament al problema
2. **Utilitzar enfocaments sistemàtics** - corregir un error a la vegada i provar amb `-resume`
3. **Entendre el flux de dades** - els errors d'estructura de canals són sovint els més subtils
4. **Comprovar directoris de treball** - quan els processos fallen, els registres us diuen exactament què va anar malament

---

## Resum

En aquesta missió secundària, heu après un conjunt de tècniques sistemàtiques per depurar workflows de Nextflow.
Aplicar aquestes tècniques al vostre propi treball us permetrà passar menys temps lluitant amb el vostre ordinador, resoldre problemes més ràpidament i protegir-vos de problemes futurs.

### Patrons clau

**1. Com identificar i corregir errors de sintaxi**:

- Interpretar missatges d'error de Nextflow i localitzar problemes
- Errors de sintaxi comuns: claus que falten, paraules clau incorrectes, variables no definides
- Distingir entre variables de Nextflow (Groovy) i Bash
- Utilitzar funcionalitats de l'extensió VS Code per a detecció primerenca d'errors

```groovy
// Clau que falta - buscar subratllats vermells a l'IDE
process FOO {
    script:
    """
    echo "hello"
    """
// } <-- falta!

// Paraula clau incorrecta
inputs:  // Hauria de ser 'input:'

// Variable no definida - escapar amb barra invertida per a variables Bash
echo "${undefined_var}"      // Variable Nextflow (error si no està definida)
echo "\${bash_var}"          // Variable Bash (escapada)
```

**2. Com depurar problemes d'estructura de canals**:

- Entendre problemes de cardinalitat i esgotament de canals
- Depurar discrepàncies d'estructura de contingut de canals
- Utilitzar operadors `.view()` per a inspecció de canals
- Reconèixer patrons de missatges d'error com claudàtors indicant estructures de tupla inesperades

```groovy
// Inspeccionar contingut de canal
my_channel.view { "Content: $it" }

// Convertir cua a canal de valor (prevé l'esgotament)
reference_ch = channel.value('ref.fa')
// o
reference_ch = channel.of('ref.fa').first()
```

**3. Com resoldre problemes d'execució de processos**:

- Diagnosticar errors de fitxers de sortida que falten
- Entendre codis de sortida (127 per a programari que falta, 137 per a problemes de memòria)
- Investigar directoris de treball i fitxers de comandes
- Configurar recursos adequadament

```bash
# Comprovar què es va executar realment
cat work/ab/cdef12/.command.sh

# Comprovar sortida d'error
cat work/ab/cdef12/.command.err

# Codi de sortida 127 = comanda no trobada
# Codi de sortida 137 = matat (límit de memòria/temps)
```

**4. Com utilitzar les eines de depuració integrades de Nextflow**:

- Aprofitar el mode de previsualització i depuració en temps real
- Implementar execució stub per a proves de lògica
- Aplicar resume per a cicles de depuració eficients
- Seguir una metodologia sistemàtica de depuració de quatre fases

!!! tip "Referència Ràpida de Depuració"

    **Errors de sintaxi?** → Comprovar advertències de VSCode, executar `nextflow run workflow.nf -preview`

    **Problemes de canals?** → Utilitzar `.view()` per inspeccionar contingut: `my_channel.view()`

    **Fallades de processos?** → Comprovar fitxers del directori de treball:

    - `.command.sh` - l'script executat
    - `.command.err` - missatges d'error
    - `.exitcode` - estat de sortida (127 = comanda no trobada, 137 = matat)

    **Comportament misteriós?** → Executar amb `-stub-run` per provar la lògica del workflow

    **Heu fet correccions?** → Utilitzar `-resume` per estalviar temps provant: `nextflow run workflow.nf -resume`

---

### Recursos addicionals

- [Guia de resolució de problemes de Nextflow](https://www.nextflow.io/docs/latest/troubleshooting.html): Documentació oficial de resolució de problemes
- [Entendre els canals de Nextflow](https://www.nextflow.io/docs/latest/channel.html): Immersió profunda en tipus i comportament de canals
- [Referència de directives de procés](https://www.nextflow.io/docs/latest/process.html#directives): Totes les opcions de configuració de procés disponibles
- [nf-test](https://www.nf-test.com/): Marc de proves per a pipelines de Nextflow
- [Comunitat Slack de Nextflow](https://www.nextflow.io/slack-invite.html): Obtenir ajuda de la comunitat

Per a workflows de producció, considereu:

- Configurar [Seqera Platform](https://seqera.io/platform/) per a monitorització i depuració a escala
- Utilitzar [contenidors Wave](https://seqera.io/wave/) per a entorns de programari reproduïbles

**Recordeu:** La depuració efectiva és una habilitat que millora amb la pràctica. La metodologia sistemàtica i el conjunt d'eines exhaustiu que heu adquirit aquí us serviran bé al llarg del vostre viatge de desenvolupament amb Nextflow.

---

## Què segueix?

Torneu al [menú de Missions Secundàries](./index.md) o feu clic al botó a la part inferior dreta de la pàgina per passar al següent tema de la llista.
