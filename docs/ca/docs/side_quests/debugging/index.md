# Depuració de Workflows

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Traducció assistida per IA - [més informació i suggeriments](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

La depuració és una habilitat crítica que us pot estalviar hores de frustració i ajudar-vos a ser un desenvolupador de Nextflow més eficient. Al llarg de la vostra carrera, especialment quan esteu començant, trobareu errors mentre construïu i manteniu els vostres workflows. Aprendre enfocaments sistemàtics de depuració us ajudarà a identificar i resoldre problemes ràpidament.

### Objectius d'aprenentatge

En aquesta missió secundària, explorarem **tècniques sistemàtiques de depuració** per a workflows de Nextflow:

- **Depuració d'errors de sintaxi**: Ús eficaç de les funcions de l'IDE i dels missatges d'error de Nextflow
- **Depuració de canals**: Diagnosi de problemes de flux de dades i problemes d'estructura de canals
- **Depuració de processos**: Investigació de fallades d'execució i problemes de recursos
- **Eines de depuració integrades**: Aprofitament del mode de previsualització de Nextflow, l'execució stub i els directoris de treball
- **Enfocaments sistemàtics**: Una metodologia de quatre fases per a una depuració eficient

Al final, tindreu una metodologia de depuració robusta que transforma els missatges d'error frustrants en fulls de ruta clars cap a les solucions.

### Prerequisits

Abans d'emprendre aquesta missió secundària, hauríeu de:

- Haver completat el tutorial [Hello Nextflow](../hello_nextflow/README.md) o un curs equivalent per a principiants.
- Estar còmodes amb els conceptes i mecanismes bàsics de Nextflow (processos, canals, operadors)

**Opcional:** Recomanem completar primer la missió secundària [IDE Features for Nextflow Development](../dev_environment/).
Aquesta cobreix les funcions de l'IDE que donen suport a la depuració (ressaltat de sintaxi, detecció d'errors, etc.), que farem servir àmpliament aquí.

---

## 0. Primers passos

#### Obriu l'espai de treball de formació

Si encara no ho heu fet, assegureu-vos d'obrir l'entorn de formació tal com es descriu a la [Configuració de l'entorn](../envsetup/index.md).

[![Open in GitHub Codespaces](https://github.com/codespaces/badge.svg)](https://codespaces.new/nextflow-io/training?quickstart=1&ref=master)

#### Moveu-vos al directori del projecte

Movem-nos al directori on es troben els fitxers d'aquest tutorial.

```bash
cd side-quests/debugging
```

Podeu configurar VSCode perquè es centri en aquest directori:

```bash
code .
```

#### Reviseu els materials

Trobareu un conjunt de workflows d'exemple amb diversos tipus d'errors que farem servir per practicar:

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

Aquests fitxers representen escenaris de depuració habituals que trobareu en el desenvolupament real.

#### Reviseu l'assignació

El vostre repte és executar cada workflow, identificar els errors i corregir-los.

Per a cada workflow amb errors:

1. **Executeu el workflow** i observeu l'error
2. **Analitzeu el missatge d'error**: què us diu Nextflow?
3. **Localitzeu el problema** al codi usant les pistes proporcionades
4. **Corregiu l'error** i verifiqueu que la vostra solució funciona
5. **Restabliu el fitxer** abans de passar a la secció següent (useu `git checkout <filename>`)

Els exercicis progressen des d'errors de sintaxi simples fins a problemes d'execució més subtils.
Les solucions es discuteixen en línia, però intenteu resoldre cada una vosaltres mateixos abans de llegir endavant.

#### Llista de comprovació de preparació

Creieu que esteu preparats per submergir-vos?

- [ ] Entenc l'objectiu d'aquest curs i els seus prerequisits
- [ ] El meu espai de treball està en funcionament
- [ ] He configurat el meu directori de treball adequadament
- [ ] Entenc l'assignació

Si podeu marcar totes les caselles, esteu a punt.

---

## 1. Errors de Sintaxi

Els errors de sintaxi són el tipus d'error més comú que trobareu quan escriviu codi Nextflow. Es produeixen quan el codi no s'ajusta a les regles de sintaxi esperades del DSL de Nextflow. Aquests errors impedeixen que el vostre workflow s'executi, de manera que és important aprendre a identificar-los i corregir-los ràpidament.

### 1.1. Claus que falten

Un dels errors de sintaxi més comuns, i de vegades un dels més complexos de depurar, és el de **claus que falten o no coincideixen**.

Comencem amb un exemple pràctic.

#### Executeu el pipeline

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

- **Fitxer i ubicació**: Mostra quin fitxer i quina línia/columna contenen l'error (`bad_syntax.nf:24:1`)
- **Descripció de l'error**: Explica què ha trobat l'analitzador que no esperava (`Unexpected input: '<EOF>'`)
- **Indicador EOF**: El missatge `<EOF>` (End Of File) indica que l'analitzador ha arribat al final del fitxer mentre encara esperava més contingut - un signe clàssic de claus no tancades

#### Comproveu el codi

Ara, examinem `bad_syntax.nf` per entendre què causa l'error:

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
// Falta la clau de tancament del procés

workflow {

    // Crea el canal d'entrada
    input_ch = channel.of('sample1', 'sample2', 'sample3')

    // Crida el procés amb el canal d'entrada
    PROCESS_FILES(input_ch)
}
```

Per a aquest exemple hem deixat un comentari per mostrar-vos on és l'error. L'extensió VSCode de Nextflow també us hauria de donar algunes pistes sobre el que podria estar malament, posant la clau no coincident en vermell i ressaltant el final prematur del fitxer:

![Bad syntax](img/bad_syntax.png)

**Estratègia de depuració per a errors de claus:**

1. Useu la coincidència de claus de VS Code (col·loqueu el cursor al costat d'una clau)
2. Comproveu el panell de Problemes per a missatges relacionats amb claus
3. Assegureu-vos que cada `{` d'obertura té un `}` de tancament corresponent

#### Corregiu el codi

Substituïu el comentari per la clau de tancament que falta:

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
    }  // Afegeix la clau de tancament que faltava

    workflow {

        // Crea el canal d'entrada
        input_ch = channel.of('sample1', 'sample2', 'sample3')

        // Crida el procés amb el canal d'entrada
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
    // Falta la clau de tancament del procés

    workflow {

        // Crea el canal d'entrada
        input_ch = channel.of('sample1', 'sample2', 'sample3')

        // Crida el procés amb el canal d'entrada
        PROCESS_FILES(input_ch)
    }
    ```

#### Executeu el pipeline

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

### 1.2. Ús de paraules clau o directives de procés incorrectes

Un altre error de sintaxi comú és una **definició de procés no vàlida**. Això pot passar si oblideu definir blocs requerits o useu directives incorrectes a la definició del procés.

#### Executeu el pipeline

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

#### Comproveu el codi

El missatge d'error indica una "definició de procés no vàlida" i mostra el context al voltant del problema. Mirant les línies 3-7, podem veure `inputs:` a la línia 4, que és el problema. Examinem `invalid_process.nf`:

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

    // Crea el canal d'entrada
    input_ch = channel.of('sample1', 'sample2', 'sample3')

    // Crida el procés amb el canal d'entrada
    PROCESS_FILES(input_ch)
}
```

Mirant la línia 4 en el context de l'error, podem detectar el problema: estem usant `inputs` en lloc de la directiva correcta `input`. L'extensió VSCode de Nextflow també ho marcarà:

![Invalid process message](img/invalid_process_message.png)

#### Corregiu el codi

Substituïu la paraula clau incorrecta per la correcta consultant [la documentació](https://www.nextflow.io/docs/latest/process.html#):

=== "Després"

    ```groovy title="invalid_process.nf" hl_lines="4" linenums="1"
    #!/usr/bin/env nextflow

    process PROCESS_FILES {
        input:  // Corregit: Canviat 'inputs' per 'input'
        val sample_name

        output:
        path "${sample_name}_output.txt"

        script:
        """
        echo "Processing ${sample_name}" > ${sample_name}_output.txt
        """
    }

    workflow {

        // Crea el canal d'entrada
        input_ch = channel.of('sample1', 'sample2', 'sample3')

        // Crida el procés amb el canal d'entrada
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

        // Crea el canal d'entrada
        input_ch = channel.of('sample1', 'sample2', 'sample3')

        // Crida el procés amb el canal d'entrada
        PROCESS_FILES(input_ch)
    }
    ```

#### Executeu el pipeline

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

### 1.3. Ús de noms de variables incorrectes

Els noms de variables que useu als vostres blocs de script han de ser vàlids, derivats bé de les entrades o bé de codi Groovy inserit abans del script. Però quan esteu gestionant complexitat al principi del desenvolupament del pipeline, és fàcil cometre errors en la nomenclatura de variables, i Nextflow us ho farà saber ràpidament.

#### Executeu el pipeline

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

L'error es detecta en temps de compilació i apunta directament a la variable no definida a la línia 17, amb un accent circumflex que indica exactament on és el problema.

#### Comproveu el codi

Examinem `no_such_var.nf`:

```groovy title="no_such_var.nf" hl_lines="17" linenums="1"
#!/usr/bin/env nextflow

process PROCESS_FILES {
    input:
    val sample_name

    output:
    path "${sample_name}_processed.txt"

    script:
    // Defineix variables en codi Groovy abans del script
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

El missatge d'error indica que la variable no es reconeix a la plantilla del script, i aquí ho teniu: hauríeu de poder veure `${undefined_var}` usat al bloc del script, però no definit en cap altre lloc.

#### Corregiu el codi

Si obteniu un error 'No such variable', podeu corregir-lo bé definint la variable (corregint els noms de variables d'entrada o editant el codi Groovy abans del script), o bé eliminant-la del bloc del script si no és necessària:

=== "Després"

    ```groovy title="no_such_var.nf" hl_lines="15-17" linenums="1"
    #!/usr/bin/env nextflow

    process PROCESS_FILES {
        input:
        val sample_name

        output:
        path "${sample_name}_output.txt"

        script:
        // Defineix variables en codi Groovy abans del script
        def output_prefix = "${sample_name}_processed"
        def timestamp = new Date().format("yyyy-MM-dd")

        """
        echo "Processing ${sample_name} on ${timestamp}" > ${output_prefix}.txt
        """  // S'ha eliminat la línia amb undefined_var
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
        // Defineix variables en codi Groovy abans del script
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

#### Executeu el pipeline

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

### 1.4. Ús incorrecte de variables Bash

Quan es comença amb Nextflow, pot ser difícil entendre la diferència entre les variables de Nextflow (Groovy) i les de Bash. Això pot generar una altra forma de l'error de variable incorrecta que apareix quan s'intenten usar variables al contingut Bash del bloc del script.

#### Executeu el pipeline

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

#### Comproveu el codi

L'error apunta a la línia 13 on s'usa `${prefix}`. Examinem `bad_bash_var.nf` per veure què causa el problema:

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

En aquest exemple, estem definint la variable `prefix` en Bash, però en un procés Nextflow la sintaxi `$` que hem usat per referir-nos-hi (`${prefix}`) s'interpreta com una variable Groovy, no Bash. La variable no existeix en el context Groovy, de manera que obtenim un error 'no such variable'.

#### Corregiu el codi

Si voleu usar una variable Bash, heu d'escapar el signe de dòlar d'aquesta manera:

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
        echo "Processing ${sample_name}" > \${prefix}.txt  # Corregit: S'ha escapat el signe de dòlar
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

Això indica a Nextflow que ho interpreti com una variable Bash.

#### Executeu el pipeline

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

    Per a manipulacions de variables simples com la concatenació de strings o operacions de prefix/sufix, normalment és més llegible usar variables Groovy a la secció del script en lloc de variables Bash al bloc del script:

    ```groovy linenums="1"
    script:
    def output_prefix = "${sample_name}_processed"
    def output_file = "${output_prefix}.txt"
    """
    echo "Processing ${sample_name}" > ${output_file}
    """
    ```

    Aquest enfocament evita la necessitat d'escapar els signes de dòlar i fa que el codi sigui més fàcil de llegir i mantenir.

### 1.5. Instruccions fora del bloc Workflow

L'extensió VSCode de Nextflow ressalta problemes amb l'estructura del codi que causaran errors. Un exemple comú és definir canals fora del bloc `workflow {}` - ara s'aplica com un error de sintaxi.

#### Executeu el pipeline

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

El missatge d'error indica clarament el problema: les instruccions (com les definicions de canals) no es poden barrejar amb les declaracions del script fora d'un bloc workflow o process.

#### Comproveu el codi

Examinem `badpractice_syntax.nf` per veure què causa l'error:

```groovy title="badpractice_syntax.nf" hl_lines="3" linenums="1"
#!/usr/bin/env nextflow

input_ch = channel.of('sample1', 'sample2', 'sample3')  // ERROR: Canal definit fora del workflow

process PROCESS_FILES {
    input:
    val sample_name

    output:
    path "${sample_name}_processed.txt"

    script:
    // Defineix variables en codi Groovy abans del script
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

L'extensió VSCode també ressaltarà la variable `input_ch` com a definida fora del bloc workflow:

![Non-lethal syntax error](img/nonlethal.png)

#### Corregiu el codi

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
        // Defineix variables en codi Groovy abans del script
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
        // Defineix variables en codi Groovy abans del script
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

#### Executeu el pipeline

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

Manteniu els vostres canals d'entrada definits dins del bloc workflow, i en general seguiu qualsevol altra recomanació que faci l'extensió.

### Conclusió

Podeu identificar i corregir sistemàticament errors de sintaxi usant els missatges d'error de Nextflow i els indicadors visuals de l'IDE. Els errors de sintaxi comuns inclouen claus que falten, paraules clau de procés incorrectes, variables no definides i ús incorrecte de variables Bash vs. Nextflow. L'extensió VSCode ajuda a detectar molts d'aquests errors abans de l'execució. Amb aquestes habilitats de depuració de sintaxi al vostre repertori, podreu resoldre ràpidament els errors de sintaxi de Nextflow més comuns i passar a abordar problemes d'execució més complexos.

### Què segueix?

Apreneu a depurar errors d'estructura de canal més complexos que es produeixen fins i tot quan la sintaxi és correcta.

---

## 2. Errors d'Estructura de Canal

Els errors d'estructura de canal són més subtils que els errors de sintaxi perquè el codi és sintàcticament correcte, però les formes de les dades no coincideixen amb el que esperen els processos. Nextflow intentarà executar el pipeline, però pot trobar que el nombre d'entrades no coincideix amb el que espera i fallar. Aquests errors normalment només apareixen en temps d'execució i requereixen una comprensió del flux de dades a través del vostre workflow.

!!! tip "Depuració de canals amb `.view()`"

    Al llarg d'aquesta secció, recordeu que podeu usar l'operador `.view()` per inspeccionar el contingut del canal en qualsevol punt del vostre workflow. Aquesta és una de les eines de depuració més potents per entendre els problemes d'estructura de canal. Explorarem aquesta tècnica en detall a la secció 2.4, però podeu usar-la lliurement mentre treballeu amb els exemples.

    ```groovy
    my_channel.view()  // Mostra el que flueix pel canal
    ```

### 2.1. Nombre incorrecte de canals d'entrada

Aquest error es produeix quan passeu un nombre diferent de canals del que espera un procés.

#### Executeu el pipeline

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

#### Comproveu el codi

El missatge d'error indica clarament que la crida esperava 1 argument però en va rebre 2, i apunta a la línia 23. Examinem `bad_number_inputs.nf`:

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

    // Crea dos canals separats
    samples_ch = channel.of('sample1', 'sample2', 'sample3')
    files_ch = channel.of('file1.txt', 'file2.txt', 'file3.txt')

    // ERROR: Es passen 2 canals però el procés espera només 1
    PROCESS_FILES(samples_ch, files_ch)
}
```

Hauríeu de veure la crida `PROCESS_FILES` no coincident, que subministra múltiples canals d'entrada quan el procés només en defineix un. L'extensió VSCode també subratllará la crida del procés en vermell i proporcionarà un missatge de diagnòstic quan hi passeu el ratolí per sobre:

![Incorrect number of args message](img/incorrect_num_args.png)

#### Corregiu el codi

Per a aquest exemple específic, el procés espera un únic canal i no requereix el segon canal, de manera que podem corregir-ho passant només el canal `samples_ch`:

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

        // Crea dos canals separats
        samples_ch = channel.of('sample1', 'sample2', 'sample3')
        files_ch = channel.of('file1.txt', 'file2.txt', 'file3.txt')

        // Corregit: Passa només el canal que espera el procés
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

        // Crea dos canals separats
        samples_ch = channel.of('sample1', 'sample2', 'sample3')
        files_ch = channel.of('file1.txt', 'file2.txt', 'file3.txt')

        // ERROR: Es passen 2 canals però el procés espera només 1
        PROCESS_FILES(samples_ch, files_ch)
    }
    ```

#### Executeu el pipeline

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

Més habitualment que en aquest exemple, podríeu afegir entrades addicionals a un procés i oblidar actualitzar la crida del workflow en conseqüència, cosa que pot portar a aquest tipus d'error. Afortunadament, aquest és un dels errors més fàcils d'entendre i corregir, ja que el missatge d'error és força clar sobre la discrepància.

### 2.2. Esgotament del canal (el procés s'executa menys vegades de les esperades)

Alguns errors d'estructura de canal són molt més subtils i no produeixen cap error. Probablement el més comú d'aquests reflecteix un repte que els nous usuaris de Nextflow afronten en entendre que els queue channels es poden esgotar i quedar-se sense elements, cosa que significa que el workflow acaba prematurament.

#### Executeu el pipeline

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

Aquest workflow es completa sense error, però només processa una única mostra!

#### Comproveu el codi

Examinem `exhausted.nf` per veure si és correcte:

```groovy title="exhausted.nf" hl_lines="23 24" linenums="1"
#!/usr/bin/env nextflow

process PROCESS_FILES {
    input:
    val reference
    val sample_name

    output:
    path "${output_prefix}.txt"

    script:
    // Defineix variables en codi Groovy abans del script
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

El procés només s'executa una vegada en lloc de tres perquè el canal `reference_ch` és un queue channel que s'esgota després de la primera execució del procés. Quan un canal s'esgota, tot el procés s'atura, fins i tot si altres canals encara tenen elements.

Aquest és un patró comú on teniu un únic fitxer de referència que cal reutilitzar en múltiples mostres. La solució és convertir el canal de referència en un value channel que es pugui reutilitzar indefinidament.

#### Corregiu el codi

Hi ha un parell de maneres d'abordar això depenent de quants fitxers es veuen afectats.

**Opció 1**: Teniu un únic fitxer de referència que reutilitzeu molt. Podeu simplement crear un value channel, que es pot usar una i altra vegada. Hi ha tres maneres de fer-ho:

**1a** Useu `channel.value()`:

```groovy title="exhausted.nf (fixed - Option 1a)" hl_lines="2" linenums="21"
workflow {
    reference_ch = channel.value('baseline_reference')  // El value channel es pot reutilitzar
    input_ch = channel.of('sample1', 'sample2', 'sample3')

    PROCESS_FILES(reference_ch, input_ch)
}
```

**1b** Useu l'[operador](https://www.nextflow.io/docs/latest/reference/operator.html#first) `first()`:

```groovy title="exhausted.nf (fixed - Option 1b)" hl_lines="2" linenums="21"
workflow {
    reference_ch = channel.of('baseline_reference').first()  // Converteix a value channel
    input_ch = channel.of('sample1', 'sample2', 'sample3')

    PROCESS_FILES(reference_ch, input_ch)
}
```

**1c.** Useu l'[operador](https://www.nextflow.io/docs/latest/reference/operator.html#collect) `collect()`:

```groovy title="exhausted.nf (fixed - Option 1c)" hl_lines="2" linenums="21"
workflow {
    reference_ch = channel.of('baseline_reference').collect()  // Converteix a value channel
    input_ch = channel.of('sample1', 'sample2', 'sample3')

    PROCESS_FILES(reference_ch, input_ch)
}
```

**Opció 2**: En escenaris més complexos, potser on teniu múltiples fitxers de referència per a totes les mostres del canal de mostres, podeu usar l'operador `combine` per crear un nou canal que combini els dos canals en tuples:

```groovy title="exhausted.nf (fixed - Option 2)" hl_lines="4" linenums="21"
workflow {
    reference_ch = channel.of('baseline_reference','other_reference')
    input_ch = channel.of('sample1', 'sample2', 'sample3')
    combined_ch = reference_ch.combine(input_ch)  // Crea el producte cartesià

    PROCESS_FILES(combined_ch)
}
```

L'operador `.combine()` genera un producte cartesià dels dos canals, de manera que cada element de `reference_ch` s'aparellarà amb cada element de `input_ch`. Això permet que el procés s'executi per a cada mostra mentre segueix usant la referència.

Això requereix que s'ajusti l'entrada del procés. En el nostre exemple, l'inici de la definició del procés s'hauria d'ajustar de la manera següent:

```groovy title="exhausted.nf (fixed - Option 2)" hl_lines="5" linenums="1"
#!/usr/bin/env nextflow

process PROCESS_FILES {
    input:
        tuple val(reference), val(sample_name)
```

Aquest enfocament pot no ser adequat en totes les situacions.

#### Executeu el pipeline

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

Ara hauríeu de veure les tres mostres processades en lloc d'una sola.

### 2.3. Estructura de contingut de canal incorrecta

Quan els workflows arriben a un cert nivell de complexitat, pot ser una mica difícil fer un seguiment de les estructures internes de cada canal, i la gent sovint genera discrepàncies entre el que espera el procés i el que conté realment el canal. Això és més subtil que el problema que hem discutit anteriorment, on el nombre de canals era incorrecte. En aquest cas, podeu tenir el nombre correcte de canals d'entrada, però l'estructura interna d'un o més d'aquests canals no coincideix amb el que espera el procés.

#### Executeu el pipeline

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

#### Comproveu el codi

Els claudàtors al missatge d'error proporcionen la pista aquí - el procés tracta la tupla com un valor únic, que no és el que volem. Examinem `bad_channel_shape.nf`:

```groovy title="bad_channel_shape.nf" hl_lines="5 20-22" linenums="1"
#!/usr/bin/env nextflow

process PROCESS_FILES {
    input:
        val sample_name  // Espera un valor únic, rep una tupla

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

Podeu veure que estem generant un canal compost de tuples: `['sample1', 'file1.txt']`, però el procés espera un valor únic, `val sample_name`. La comanda executada mostra que el procés intenta crear un fitxer anomenat `[sample3, file3.txt]_output.txt`, que no és la sortida prevista.

#### Corregiu el codi

Per corregir-ho, si el procés requereix ambdues entrades podríem ajustar el procés per acceptar una tupla:

=== "Opció 1: Acceptar tupla al procés"

    === "Després"

        ```groovy title="bad_channel_shape.nf" hl_lines="5"  linenums="1"
        #!/usr/bin/env nextflow

        process PROCESS_FILES {
            input:
                tuple val(sample_name), val(file_name)  // Corregit: Accepta tupla

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
                val sample_name  // Espera un valor únic, rep una tupla

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
            PROCESS_FILES(input_ch.map { it[0] })  // Corregit: Extreu el primer element
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

#### Executeu el pipeline

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

### 2.4. Tècniques de depuració de canals

#### Ús de `.view()` per a la inspecció de canals

L'eina de depuració més potent per als canals és l'operador `.view()`. Amb `.view()`, podeu entendre la forma dels vostres canals en totes les etapes per ajudar amb la depuració.

#### Executeu el pipeline

Executeu `bad_channel_shape_viewed.nf` per veure-ho en acció:

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

#### Comproveu el codi

Examinem `bad_channel_shape_viewed.nf` per veure com s'usa `.view()`:

```groovy title="bad_channel_shape_viewed.nf" linenums="16" hl_lines="9 11"
workflow {

    // El canal emet tuples, però el procés espera valors únics
    input_ch = channel.of(
      ['sample1', 'file1.txt'],
      ['sample2', 'file2.txt'],
      ['sample3', 'file3.txt']
    )
    .view { "Channel content: $it" }  // Depuració: Mostra el contingut original del canal
    .map { tuple -> tuple[0] }        // Transformació: Extreu el primer element
    .view { "After mapping: $it" }    // Depuració: Mostra el contingut del canal transformat

    PROCESS_FILES(input_ch)
}
```

#### Corregiu el codi

Per estalviar-vos d'usar operacions `.view()` excessivament en el futur per entendre el contingut del canal, és aconsellable afegir alguns comentaris per ajudar:

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

Això serà més important a mesura que els vostres workflows creixin en complexitat i l'estructura del canal es torni més opaca.

#### Executeu el pipeline

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

Molts errors d'estructura de canal es poden crear amb sintaxi Nextflow vàlida. Podeu depurar errors d'estructura de canal entenent el flux de dades, usant operadors `.view()` per a la inspecció i reconeixent patrons de missatges d'error com els claudàtors que indiquen estructures de tupla inesperades.

### Què segueix?

Apreneu sobre els errors creats per les definicions de processos.

---

## 3. Errors d'Estructura de Procés

La majoria dels errors que trobareu relacionats amb els processos estaran relacionats amb errors que heu comès en formar la comanda, o amb problemes relacionats amb el programari subjacent. Dit això, de manera similar als problemes de canal anteriors, podeu cometre errors a la definició del procés que no qualifiquen com a errors de sintaxi, però que causaran errors en temps d'execució.

### 3.1. Fitxers de sortida que falten

Un error comú en escriure processos és fer alguna cosa que genera una discrepància entre el que espera el procés i el que es genera.

#### Executeu el pipeline

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

#### Comproveu el codi

El missatge d'error indica que el procés esperava produir un fitxer de sortida anomenat `sample3.txt`, però el script en realitat crea `sample3_output.txt`. Examinem la definició del procés a `missing_output.nf`:

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

Hauríeu de veure que hi ha una discrepància entre el nom del fitxer de sortida al bloc `output:` i el que s'usa al script. Aquesta discrepància fa que el procés falli. Si trobeu aquest tipus d'error, torneu enrere i comproveu que les sortides coincideixen entre la definició del vostre procés i el vostre bloc de sortida.

Si el problema encara no és clar, comproveu el directori de treball per identificar els fitxers de sortida reals creats:

```bash
❯ ls -h work/02/9604d49fb8200a74d737c72a6c98ed
sample3_output.txt
```

Per a aquest exemple, això ens destacaria que s'incorpora un sufix `_output` al nom del fitxer de sortida, contràriament a la nostra definició `output:`.

#### Corregiu el codi

Corregiu la discrepància fent que el nom del fitxer de sortida sigui consistent:

=== "Després"

    ```groovy title="missing_output.nf" hl_lines="6 10" linenums="3"
    process PROCESS_FILES {
        input:
        val sample_name

        output:
        path "${sample_name}_output.txt"  // Corregit: Coincideix amb la sortida del script

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

#### Executeu el pipeline

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

Una altra classe d'errors es produeix a causa d'errors en el proveïment de programari. `missing_software.nf` és un workflow sintàcticament vàlid, però depèn d'algun programari extern per proporcionar la comanda `cowpy` que usa.

#### Executeu el pipeline

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

El procés no té accés a la comanda que estem especificant. De vegades és perquè un script és present al directori `bin` del workflow, però no s'ha fet executable. Altres vegades és perquè el programari no està instal·lat al contenidor o entorn on s'executa el workflow.

#### Comproveu el codi

Fixeu-vos en aquell codi de sortida `127` - us indica exactament el problema. Examinem `missing_software.nf`:

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

#### Corregiu el codi

Hem estat una mica poc sincers aquí, i en realitat no hi ha res malament amb el codi. Simplement hem d'especificar la configuració necessària per executar el procés de manera que tingui accés a la comanda en qüestió. En aquest cas el procés té una definició de contenidor, de manera que tot el que hem de fer és executar el workflow amb Docker habilitat.

#### Executeu el pipeline

Hem configurat un perfil Docker per a vosaltres a `nextflow.config`, de manera que podeu executar el workflow amb:

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

    Per aprendre més sobre com Nextflow usa els contenidors, vegeu [Hello Nextflow](../hello_nextflow/05_hello_containers.md)

### 3.3. Configuració de recursos incorrecta

En ús en producció, estareu configurant recursos als vostres processos. Per exemple, `memory` defineix la quantitat màxima de memòria disponible per al vostre procés, i si el procés la supera, el vostre planificador normalment matar el procés i retornarà un codi de sortida de `137`. No podem demostrar-ho aquí perquè estem usant l'executor `local`, però podem mostrar alguna cosa similar amb `time`.

#### Executeu el pipeline

`bad_resources.nf` té una configuració de procés amb un límit de temps poc realista d'1 mil·lisegon:

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

#### Comproveu el codi

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

Sabem que el procés trigarà més d'un segon (hem afegit un sleep per assegurar-nos-en), però el procés està configurat per expirar després d'1 mil·lisegon. Algú ha estat una mica poc realista amb la seva configuració!

#### Corregiu el codi

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

#### Executeu el pipeline

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

Si us assegureu de llegir els vostres missatges d'error, fallades com aquesta no us haurien de desconcertar durant massa temps. Però assegureu-vos d'entendre els requisits de recursos de les comandes que esteu executant per poder configurar les vostres directives de recursos adequadament.

### 3.4. Tècniques de depuració de processos

Quan els processos fallen o es comporten de manera inesperada, necessiteu tècniques sistemàtiques per investigar el que ha anat malament. El directori de treball conté tota la informació que necessiteu per depurar l'execució del procés.

#### Ús de la inspecció del directori de treball

L'eina de depuració més potent per als processos és examinar el directori de treball. Quan un procés falla, Nextflow crea un directori de treball per a aquesta execució específica del procés que conté tots els fitxers necessaris per entendre el que ha passat.

#### Executeu el pipeline

Usem l'exemple `missing_output.nf` d'abans per demostrar la inspecció del directori de treball (torneu a generar una discrepància de nomenclatura de sortida si cal):

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

#### Comproveu el directori de treball

Quan obteniu aquest error, el directori de treball conté tota la informació de depuració. Trobeu la ruta del directori de treball del missatge d'error i examineu el seu contingut:

```bash
# Trobeu el directori de treball del missatge d'error
ls work/02/9604d49fb8200a74d737c72a6c98ed/
```

A continuació podeu examinar els fitxers clau:

##### Comproveu l'script de la comanda

El fitxer `.command.sh` mostra exactament quina comanda s'ha executat:

```bash
# Visualitzeu la comanda executada
cat work/02/9604d49fb8200a74d737c72a6c98ed/.command.sh
```

Això revela:

- **Substitució de variables**: Si les variables de Nextflow s'han expandit correctament
- **Rutes de fitxers**: Si els fitxers d'entrada s'han localitzat correctament
- **Estructura de la comanda**: Si la sintaxi del script és correcta

Problemes comuns a buscar:

- **Cometes que falten**: Les variables que contenen espais necessiten cometes adequades
- **Rutes de fitxers incorrectes**: Fitxers d'entrada que no existeixen o estan en ubicacions incorrectes
- **Noms de variables incorrectes**: Errors tipogràfics en les referències de variables
- **Configuració d'entorn que falta**: Comandes que depenen d'entorns específics

##### Comproveu la sortida d'error

El fitxer `.command.err` conté els missatges d'error reals:

```bash
# Visualitzeu la sortida d'error
cat work/02/9604d49fb8200a74d737c72a6c98ed/.command.err
```

Aquest fitxer mostrarà:

- **Codis de sortida**: 127 (comanda no trobada), 137 (eliminat), etc.
- **Errors de permisos**: Problemes d'accés a fitxers
- **Errors de programari**: Missatges d'error específics de l'aplicació
- **Errors de recursos**: Memòria/límit de temps superat

##### Comproveu la sortida estàndard

El fitxer `.command.out` mostra el que ha produït la vostra comanda:

```bash
# Visualitzeu la sortida estàndard
cat work/02/9604d49fb8200a74d737c72a6c98ed/.command.out
```

Això ajuda a verificar:

- **Sortida esperada**: Si la comanda ha produït els resultats correctes
- **Execució parcial**: Si la comanda ha començat però ha fallat a mig camí
- **Informació de depuració**: Qualsevol sortida de diagnòstic del vostre script

##### Comproveu el codi de sortida

El fitxer `.exitcode` conté el codi de sortida del procés:

```bash
# Visualitzeu el codi de sortida
cat work/*/*/.exitcode
```

Codis de sortida comuns i els seus significats:

- **Codi de sortida 127**: Comanda no trobada - comproveu la instal·lació del programari
- **Codi de sortida 137**: Procés eliminat - comproveu els límits de memòria/temps

##### Comproveu l'existència de fitxers

Quan els processos fallen a causa de fitxers de sortida que falten, comproveu quins fitxers s'han creat realment:

```bash
# Llisteu tots els fitxers al directori de treball
ls -la work/02/9604d49fb8200a74d737c72a6c98ed/
```

Això ajuda a identificar:

- **Discrepàncies en els noms de fitxers**: Fitxers de sortida amb noms diferents dels esperats
- **Problemes de permisos**: Fitxers que no s'han pogut crear
- **Problemes de ruta**: Fitxers creats en directoris incorrectes

En el nostre exemple anterior, això ens va confirmar que mentre el nostre `sample3.txt` esperat no era present, `sample3_output.txt` sí que ho era:

```bash
❯ ls -h work/02/9604d49fb8200a74d737c72a6c98ed
sample3_output.txt
```

### Conclusió

La depuració de processos requereix examinar els directoris de treball per entendre el que ha anat malament. Els fitxers clau inclouen `.command.sh` (l'script executat), `.command.err` (missatges d'error) i `.command.out` (sortida estàndard). Els codis de sortida com 127 (comanda no trobada) i 137 (procés eliminat) proporcionen pistes de diagnòstic immediates sobre el tipus de fallada.

### Què segueix?

Apreneu sobre les eines de depuració integrades de Nextflow i els enfocaments sistemàtics per a la resolució de problemes.

---

## 4. Eines de Depuració Integrades i Tècniques Avançades

Nextflow proporciona diverses eines integrades potents per depurar i analitzar l'execució del workflow. Aquestes eines us ajuden a entendre el que ha anat malament, on ha anat malament i com corregir-ho eficientment.

### 4.1. Sortida de procés en temps real

De vegades necessiteu veure el que passa dins dels processos en execució. Podeu habilitar la sortida de procés en temps real, que us mostra exactament el que fa cada tasca mentre s'executa.

#### Executeu el pipeline

`bad_channel_shape_viewed.nf` dels nostres exemples anteriors imprimia el contingut del canal usant `.view()`, però també podem usar la directiva `debug` per fer eco de variables des de dins del procés mateix, que demostrem a `bad_channel_shape_viewed_debug.nf`. Executeu el workflow:

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

#### Comproveu el codi

Examinem `bad_channel_shape_viewed_debug.nf` per veure com funciona la directiva `debug`:

```groovy title="bad_channel_shape_viewed_debug.nf" linenums="3" hl_lines="2"
process PROCESS_FILES {
    debug true  // Habilita la sortida en temps real

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

### 4.2. Mode de previsualització

De vegades voleu detectar problemes abans que s'executi cap procés. Nextflow proporciona un indicador per a aquest tipus de depuració proactiva: `-preview`.

#### Executeu el pipeline

El mode de previsualització us permet provar la lògica del workflow sense executar comandes. Això pot ser molt útil per comprovar ràpidament l'estructura del vostre workflow i assegurar-vos que els processos estan connectats correctament sense executar cap comanda real.

!!! note "Nota"

    Si heu corregit `bad_syntax.nf` anteriorment, torneu a introduir l'error de sintaxi eliminant la clau de tancament després del bloc del script abans d'executar aquesta comanda.

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

El mode de previsualització és particularment útil per detectar errors de sintaxi aviat sense executar cap procés. Valida l'estructura del workflow i les connexions dels processos abans de l'execució.

### 4.3. Execució stub per a proves de lògica

De vegades els errors són difícils de depurar perquè les comandes triguen massa, requereixen programari especial o fallen per raons complexes. L'execució stub us permet provar la lògica del workflow sense executar les comandes reals.

#### Executeu el pipeline

Quan esteu desenvolupant un procés Nextflow, podeu usar la directiva `stub` per definir comandes 'fictícies' que generin sortides de la forma correcta sense executar la comanda real. Aquest enfocament és particularment valuós quan voleu verificar que la lògica del vostre workflow és correcta abans d'afrontar les complexitats del programari real.

Per exemple, recordeu el nostre `missing_software.nf` d'abans? El que tenia programari que faltava i impedia que el workflow s'executés fins que afegíem `-profile docker`? `missing_software_with_stub.nf` és un workflow molt similar. Si l'executem de la mateixa manera, generarem el mateix error:

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

Tanmateix, aquest workflow no produirà errors si l'executem amb `-stub-run`, fins i tot sense el perfil `docker`:

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

#### Comproveu el codi

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

En relació amb `missing_software.nf`, aquest procés té una directiva `stub:` que especifica una comanda que s'ha d'usar en lloc de la especificada a `script:`, en el cas que Nextflow s'executi en mode stub.

La comanda `touch` que estem usant aquí no depèn de cap programari ni entrades adequades, i s'executarà en totes les situacions, permetent-nos depurar la lògica del workflow sense preocupar-nos pels internals del procés.

**L'execució stub ajuda a depurar:**

- Estructura del canal i flux de dades
- Connexions i dependències de processos
- Propagació de paràmetres
- Lògica del workflow sense dependències de programari

### 4.4. Enfocament sistemàtic de depuració

Ara que heu après tècniques de depuració individuals - des de fitxers de traça i directoris de treball fins al mode de previsualització, l'execució stub i la monitorització de recursos - liguem-les totes en una metodologia sistemàtica. Tenir un enfocament estructurat us impedeix sentir-vos aclaparats per errors complexos i assegura que no us perdeu pistes importants.

Aquesta metodologia combina totes les eines que hem cobert en un workflow eficient:

**Mètode de depuració de quatre fases:**

**Fase 1: Resolució d'errors de sintaxi (5 minuts)**

1. Comproveu si hi ha subratllats vermells a VSCode o al vostre IDE
2. Executeu `nextflow run workflow.nf -preview` per identificar problemes de sintaxi
3. Corregiu tots els errors de sintaxi (claus que falten, comes finals, etc.)
4. Assegureu-vos que el workflow s'analitza correctament abans de continuar

**Fase 2: Avaluació ràpida (5 minuts)**

1. Llegiu els missatges d'error d'execució amb atenció
2. Comproveu si és un error d'execució, de lògica o de recursos
3. Useu el mode de previsualització per provar la lògica bàsica del workflow

**Fase 3: Investigació detallada (15-30 minuts)**

1. Trobeu el directori de treball de la tasca fallida
2. Examineu els fitxers de registre
3. Afegiu operadors `.view()` per inspeccionar els canals
4. Useu `-stub-run` per provar la lògica del workflow sense execució

**Fase 4: Correcció i validació (15 minuts)**

1. Feu correccions mínimes i específiques
2. Proveu amb resume: `nextflow run workflow.nf -resume`
3. Verifiqueu l'execució completa del workflow

!!! tip "Ús de Resume per a una depuració eficient"

    Un cop heu identificat un problema, necessiteu una manera eficient de provar les vostres correccions sense perdre temps tornant a executar les parts reeixides del vostre workflow. La funcionalitat `-resume` de Nextflow és inestimable per a la depuració.

    Haureu trobat `-resume` si heu treballat a través de [Hello Nextflow](../hello_nextflow/), i és important que en feu bon ús quan depureu per estalviar-vos esperar mentre s'executen els processos anteriors al vostre procés problemàtic.

    **Estratègia de depuració amb resume:**

    1. Executeu el workflow fins a la fallada
    2. Examineu el directori de treball de la tasca fallida
    3. Corregiu el problema específic
    4. Repreneu per provar només la correcció
    5. Repetiu fins que el workflow es completi

#### Perfil de configuració de depuració

Per fer aquest enfocament sistemàtic encara més eficient, podeu crear una configuració de depuració dedicada que habiliti automàticament totes les eines que necessiteu:

```groovy title="nextflow.config (debug profile)" linenums="1"
profiles {
    debug {
        process {
            debug = true
            cleanup = false

            // Recursos conservadors per a la depuració
            maxForks = 1
            memory = '2.GB'
            cpus = 1
        }
    }
}
```

A continuació podeu executar el pipeline amb aquest perfil habilitat:

```bash
nextflow run workflow.nf -profile debug
```

Aquest perfil habilita la sortida en temps real, preserva els directoris de treball i limita la paral·lelització per facilitar la depuració.

### 4.5. Exercici pràctic de depuració

Ara és el moment de posar en pràctica l'enfocament sistemàtic de depuració. El workflow `buggy_workflow.nf` conté diversos errors comuns que representen els tipus de problemes que trobareu en el desenvolupament real.

!!! exercise "Exercici"

    Useu l'enfocament sistemàtic de depuració per identificar i corregir tots els errors a `buggy_workflow.nf`. Aquest workflow intenta processar dades de mostra d'un fitxer CSV però conté múltiples errors intencionals que representen escenaris de depuració comuns.

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

        Aquest error críptic indica un problema d'anàlisi al voltant de les línies 11-12 al bloc `params{}`. L'analitzador v2 detecta problemes estructurals aviat.

    Apliqueu el mètode de depuració de quatre fases que heu après:

    **Fase 1: Resolució d'errors de sintaxi**
    - Comproveu si hi ha subratllats vermells a VSCode o al vostre IDE
    - Executeu `nextflow run workflow.nf -preview` per identificar problemes de sintaxi
    - Corregiu tots els errors de sintaxi (claus que falten, comes finals, etc.)
    - Assegureu-vos que el workflow s'analitza correctament abans de continuar

    **Fase 2: Avaluació ràpida**
    - Llegiu els missatges d'error d'execució amb atenció
    - Identifiqueu si els errors són d'execució, de lògica o de recursos
    - Useu el mode `-preview` per provar la lògica bàsica del workflow

    **Fase 3: Investigació detallada**
    - Examineu els directoris de treball de les tasques fallides
    - Afegiu operadors `.view()` per inspeccionar els canals
    - Comproveu els fitxers de registre als directoris de treball
    - Useu `-stub-run` per provar la lògica del workflow sense execució

    **Fase 4: Correcció i validació**
    - Feu correccions específiques
    - Useu `-resume` per provar les correccions eficientment
    - Verifiqueu l'execució completa del workflow

    **Eines de depuració a la vostra disposició:**
    ```bash
    # Mode de previsualització per a la comprovació de sintaxi
    nextflow run buggy_workflow.nf -preview

    # Perfil de depuració per a sortida detallada
    nextflow run buggy_workflow.nf -profile debug

    # Execució stub per a proves de lògica
    nextflow run buggy_workflow.nf -stub-run

    # Resume després de les correccions
    nextflow run buggy_workflow.nf -resume
    ```

    ??? solution "Solució"
        El `buggy_workflow.nf` conté 9 o 10 errors diferents (depenent de com es comptin) que cobreixen totes les categories principals de depuració. Aquí teniu un desglossament sistemàtic de cada error i com corregir-lo.

        Comencem amb els errors de sintaxi:

        **Error 1: Error de sintaxi - Coma final**
        ```groovy linenums="21"
        output:
            path "${sample_id}_result.txt",  // ERROR: Coma final
        ```
        **Correcció:** Elimineu la coma final
        ```groovy linenums="21"
        output:
            path "${sample_id}_result.txt"
        ```

        **Error 2: Error de sintaxi - Clau de tancament que falta**
        ```groovy linenums="24"
        script:
        """
        echo "Processing: ${sample}"
        cat ${input_file} > ${sample}_result.txt
        """
        // ERROR: Falta la clau de tancament del procés processFiles
        ```
        **Correcció:** Afegiu la clau de tancament que falta
        ```groovy linenums="29"
        """
        echo "Processing: ${sample_id}"
        cat ${input_file} > ${sample_id}_result.txt
        """
        }  // Afegiu la clau de tancament que falta
        ```

        **Error 3: Error de nom de variable**
        ```groovy linenums="26"
        echo "Processing: ${sample}"     // ERROR: hauria de ser sample_id
        cat ${input_file} > ${sample}_result.txt  // ERROR: hauria de ser sample_id
        ```
        **Correcció:** Useu el nom de variable d'entrada correcte
        ```groovy linenums="26"
        echo "Processing: ${sample_id}"
        cat ${input_file} > ${sample_id}_result.txt
        ```

        **Error 4: Error de variable no definida**
        ```groovy linenums="87"
        heavy_ch = heavyProcess(sample_ids)  // ERROR: sample_ids no definida
        ```
        **Correcció:** Useu el canal correcte i extraieu els IDs de mostra
        ```groovy linenums="87"
        heavy_ch = heavyProcess(input_ch)
        ```

        En aquest punt el workflow s'executarà, però seguirem obtenint errors (p. ex. `Path value cannot be null` a `processFiles`), causats per una estructura de canal incorrecta.

        **Error 5: Error d'estructura de canal - Sortida de map incorrecta**
        ```groovy linenums="83"
        .map { row -> row.sample_id }  // ERROR: processFiles espera una tupla
        ```
        **Correcció:** Retorneu l'estructura de tupla que espera processFiles
        ```groovy linenums="83"
        .map { row -> [row.sample_id, file(row.fastq_path)] }
        ```

        Però això trencarà la nostra crida per executar `heavyProcess()` anterior, de manera que haurem d'usar un map per passar només els IDs de mostra a aquell procés:

        **Error 6: Estructura de canal incorrecta per a heavyProcess**
        ```groovy linenums="87"
        heavy_ch = heavyProcess(input_ch)  // ERROR: input_ch ara té 2 elements per emissió - heavyProcess només necessita 1 (el primer)
        ```
        **Correcció:** Useu el canal correcte i extraieu els IDs de mostra
        ```groovy linenums="87"
        heavy_ch = heavyProcess(input_ch.map{it[0]})
        ```

        Ara avancem una mica més però rebem un error sobre `No such variable: i`, perquè no hem escapat una variable Bash.

        **Error 7: Error d'escapament de variable Bash**
        ```groovy linenums="48"
        echo "Heavy computation $i for ${sample_id}"  // ERROR: $i no escapada
        ```
        **Correcció:** Escapeu la variable bash
        ```groovy linenums="48"
        echo "Heavy computation \${i} for ${sample_id}"
        ```

        Ara obtenim `Process exceeded running time limit (1ms)`, de manera que corregim el límit de temps d'execució per al procés rellevant:

        **Error 8: Error de configuració de recursos**
        ```groovy linenums="36"
        time '1 ms'  // ERROR: Límit de temps poc realista
        ```
        **Correcció:** Augmenteu a un límit de temps realista
        ```groovy linenums="36"
        time '100 s'
        ```

        A continuació tenim un error `Missing output file(s)` a resoldre:

        **Error 9: Discrepància en el nom del fitxer de sortida**
        ```groovy linenums="49"
        done > ${sample_id}.txt  // ERROR: Nom de fitxer incorrecte, hauria de coincidir amb la declaració de sortida
        ```
        **Correcció:** Feu coincidir la declaració de sortida
        ```groovy linenums="49"
        done > ${sample_id}_heavy.txt
        ```

        Els dos primers processos s'han executat, però no el tercer.

        **Error 10: Discrepància en el nom del fitxer de sortida**
        ```groovy linenums="88"
        file_ch = channel.fromPath("*.txt") // Error: intentant prendre l'entrada del pwd en lloc d'un procés
        handleFiles(file_ch)
        ```
        **Correcció:** Preneu la sortida del procés anterior
        ```groovy linenums="88"
        handleFiles(heavyProcess.out)
        ```

        Amb això, tot el workflow hauria d'executar-se.

        **Workflow corregit complet:**
        ```groovy linenums="1"
        #!/usr/bin/env nextflow

        /*
        * Workflow amb errors per a exercicis de depuració
        * Aquest workflow conté diversos errors intencionals amb finalitats d'aprenentatge
        */

        params{
            // Paràmetres sense validació
            input: Path = 'data/sample_data.csv'
            output: String = 'results'
        }

        /*
        * Procés amb discrepància d'entrada/sortida
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
            # Simula un càlcul pesat
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
        * Workflow principal amb problemes de canal
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

**Categories d'errors cobertes:**

- **Errors de sintaxi**: Claus que falten, comes finals, variables no definides
- **Errors d'estructura de canal**: Formes de dades incorrectes, canals no definits
- **Errors de procés**: Discrepàncies en noms de fitxers de sortida, escapament de variables
- **Errors de recursos**: Límits de temps poc realistes

**Lliçons clau de depuració:**

1. **Llegiu els missatges d'error amb atenció** - sovint apunten directament al problema
2. **Useu enfocaments sistemàtics** - corregiu un error a la vegada i proveu amb `-resume`
3. **Enteneu el flux de dades** - els errors d'estructura de canal sovint són els més subtils
4. **Comproveu els directoris de treball** - quan els processos fallen, els registres us diuen exactament el que ha anat malament

---

## Resum

En aquesta missió secundària, heu après un conjunt de tècniques sistemàtiques per depurar workflows de Nextflow.
Aplicar aquestes tècniques en el vostre propi treball us permetrà passar menys temps lluitant amb l'ordinador, resoldre problemes més ràpidament i protegir-vos de problemes futurs.

### Patrons clau

**1. Com identificar i corregir errors de sintaxi**:

- Interpretació dels missatges d'error de Nextflow i localització de problemes
- Errors de sintaxi comuns: claus que falten, paraules clau incorrectes, variables no definides
- Distinció entre variables de Nextflow (Groovy) i Bash
- Ús de les funcions de l'extensió VS Code per a la detecció d'errors aviat

```groovy
// Clau que falta - busqueu subratllats vermells a l'IDE
process FOO {
    script:
    """
    echo "hello"
    """
// } <-- falta!

// Paraula clau incorrecta
inputs:  // Hauria de ser 'input:'

// Variable no definida - escapeu amb barra invertida per a variables Bash
echo "${undefined_var}"      // Variable Nextflow (error si no està definida)
echo "\${bash_var}"          // Variable Bash (escapada)
```

**2. Com depurar problemes d'estructura de canal**:

- Comprensió de la cardinalitat del canal i els problemes d'esgotament
- Depuració de discrepàncies en l'estructura del contingut del canal
- Ús d'operadors `.view()` per a la inspecció del canal
- Reconeixement de patrons d'error com els claudàtors a la sortida

```groovy
// Inspeccioneu el contingut del canal
my_channel.view { "Content: $it" }

// Convertiu queue channel a value channel (evita l'esgotament)
reference_ch = channel.value('ref.fa')
// o
reference_ch = channel.of('ref.fa').first()
```

**3. Com solucionar problemes d'execució de processos**:

- Diagnosi d'errors de fitxers de sortida que falten
- Comprensió dels codis de sortida (127 per a programari que falta, 137 per a problemes de memòria)
- Investigació de directoris de treball i fitxers de comandes
- Configuració adequada dels recursos

```bash
# Comproveu el que s'ha executat realment
cat work/ab/cdef12/.command.sh

# Comproveu la sortida d'error
cat work/ab/cdef12/.command.err

# Codi de sortida 127 = comanda no trobada
# Codi de sortida 137 = eliminat (límit de memòria/temps)
```

**4. Com usar les eines de depuració integrades de Nextflow**:

- Aprofitament del mode de previsualització i la depuració en temps real
- Implementació de l'execució stub per a proves de lògica
- Aplicació de resume per a cicles de depuració eficients
- Seguiment d'una metodologia sistemàtica de depuració de quatre fases

!!! tip "Referència ràpida de depuració"

    **Errors de sintaxi?** → Comproveu els avisos de VSCode, executeu `nextflow run workflow.nf -preview`

    **Problemes de canal?** → Useu `.view()` per inspeccionar el contingut: `my_channel.view()`

    **Fallades de procés?** → Comproveu els fitxers del directori de treball:

    - `.command.sh` - l'script executat
    - `.command.err` - missatges d'error
    - `.exitcode` - estat de sortida (127 = comanda no trobada, 137 = eliminat)

    **Comportament misteriós?** → Executeu amb `-stub-run` per provar la lògica del workflow

    **Heu fet correccions?** → Useu `-resume` per estalviar temps en les proves: `nextflow run workflow.nf -resume`

---

### Recursos addicionals

- [Guia de resolució de problemes de Nextflow](https://www.nextflow.io/docs/latest/troubleshooting.html): Documentació oficial de resolució de problemes
- [Comprensió dels canals de Nextflow](https://www.nextflow.io/docs/latest/channel.html): Immersió profunda en els tipus de canal i el comportament
- [Referència de directives de procés](https://www.nextflow.io/docs/latest/process.html#directives): Totes les opcions de configuració de procés disponibles
- [nf-test](https://www.nf-test.com/): Marc de proves per a pipelines de Nextflow
- [Comunitat de Slack de Nextflow](https://www.nextflow.io/slack-invite.html): Obteniu ajuda de la comunitat

Per a workflows en producció, considereu:

- Configurar [Seqera Platform](https://seqera.io/platform/) per a la monitorització i depuració a escala
- Usar [contenidors Wave](https://seqera.io/wave/) per a entorns de programari reproduïbles

**Recordeu:** La depuració eficaç és una habilitat que millora amb la pràctica. La metodologia sistemàtica i el conjunt d'eines complet que heu adquirit aquí us serviran bé al llarg del vostre camí de desenvolupament amb Nextflow.

---

## Què segueix?

Torneu al [menú de missions secundàries](../) o feu clic al botó a la part inferior dreta de la pàgina per passar al tema següent de la llista.
