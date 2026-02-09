# Depuração de Fluxos de Trabalho

A depuração é uma habilidade crítica que pode economizar horas de frustração e ajudá-lo a se tornar um desenvolvedor Nextflow mais eficaz. Ao longo de sua carreira, especialmente quando você está começando, você encontrará bugs ao construir e manter seus fluxos de trabalho. Aprender abordagens sistemáticas de depuração ajudará você a identificar e resolver problemas rapidamente.

### Objetivos de aprendizado

Nesta missão secundária, exploraremos **técnicas sistemáticas de depuração** para fluxos de trabalho Nextflow:

- **Depuração de erros de sintaxe**: Usando recursos de IDE e mensagens de erro do Nextflow de forma eficaz
- **Depuração de canais**: Diagnosticando problemas de fluxo de dados e problemas de estrutura de canais
- **Depuração de processos**: Investigando falhas de execução e problemas de recursos
- **Ferramentas de depuração integradas**: Aproveitando o modo de visualização, execução stub e diretórios de trabalho do Nextflow
- **Abordagens sistemáticas**: Uma metodologia de quatro fases para depuração eficiente

Ao final, você terá uma metodologia robusta de depuração que transforma mensagens de erro frustrantes em roteiros claros para soluções.

### Pré-requisitos

Antes de assumir esta missão secundária, você deve:

- Ter completado o tutorial [Hello Nextflow](../hello_nextflow/README.md) ou curso equivalente para iniciantes.
- Estar confortável usando conceitos e mecanismos básicos do Nextflow (processos, canais, operadores)

**Opcional:** Recomendamos completar a missão secundária [Recursos de IDE para Desenvolvimento Nextflow](./ide_features.md) primeiro.
Ela cobre uma cobertura abrangente de recursos de IDE que suportam depuração (destaque de sintaxe, detecção de erros, etc.), que usaremos bastante aqui.

---

## 0. Primeiros passos

#### Abra o codespace de treinamento

Se você ainda não o fez, certifique-se de abrir o ambiente de treinamento conforme descrito em [Configuração do Ambiente](../envsetup/index.md).

[![Abrir no GitHub Codespaces](https://github.com/codespaces/badge.svg)](https://codespaces.new/nextflow-io/training?quickstart=1&ref=master)

#### Mova-se para o diretório do projeto

Vamos mover para o diretório onde os arquivos para este tutorial estão localizados.

```bash
cd side-quests/debugging
```

Você pode configurar o VSCode para focar neste diretório:

```bash
code .
```

#### Revise os materiais

Você encontrará um conjunto de fluxos de trabalho de exemplo com vários tipos de bugs que usaremos para prática:

??? abstract "Conteúdo do diretório"

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

Esses arquivos representam cenários comuns de depuração que você encontrará no desenvolvimento do mundo real.

#### Revise a tarefa

Seu desafio é executar cada fluxo de trabalho, identificar o(s) erro(s) e corrigi-los.

Para cada fluxo de trabalho com bugs:

1. **Execute o fluxo de trabalho** e observe o erro
2. **Analise a mensagem de erro**: o que o Nextflow está dizendo?
3. **Localize o problema** no código usando as pistas fornecidas
4. **Corrija o bug** e verifique se sua solução funciona
5. **Redefina o arquivo** antes de passar para a próxima seção (use `git checkout <filename>`)

Os exercícios progridem de erros de sintaxe simples para problemas de tempo de execução mais sutis.
As soluções são discutidas inline, mas tente resolver cada uma você mesmo antes de ler adiante.

#### Lista de verificação de prontidão

Acha que está pronto para mergulhar?

- [ ] Entendo o objetivo deste curso e seus pré-requisitos
- [ ] Meu codespace está funcionando
- [ ] Defini meu diretório de trabalho apropriadamente
- [ ] Entendo a tarefa

Se você pode marcar todas as caixas, está pronto para começar.

---

## 1. Erros de Sintaxe

Erros de sintaxe são o tipo mais comum de erro que você encontrará ao escrever código Nextflow. Eles ocorrem quando o código não está em conformidade com as regras de sintaxe esperadas do DSL Nextflow. Esses erros impedem que seu fluxo de trabalho seja executado, então é importante aprender a identificá-los e corrigi-los rapidamente.

### 1.1. Chaves ausentes

Um dos erros de sintaxe mais comuns, e às vezes um dos mais complexos de depurar, são **colchetes ausentes ou incompatíveis**.

Vamos começar com um exemplo prático.

#### Execute o pipeline

```bash
nextflow run bad_syntax.nf
```

??? failure "Saída do comando"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `bad_syntax.nf` [stupefied_bhabha] DSL2 - revision: ca6327fad2

    Error bad_syntax.nf:24:1: Unexpected input: '<EOF>'

    ERROR ~ Script compilation failed

     -- Check '.nextflow.log' file for details
    ```

**Elementos-chave das mensagens de erro de sintaxe:**

- **Arquivo e localização**: Mostra qual arquivo e linha/coluna contêm o erro (`bad_syntax.nf:24:1`)
- **Descrição do erro**: Explica o que o analisador encontrou que não esperava (`Unexpected input: '<EOF>'`)
- **Indicador EOF**: A mensagem `<EOF>` (End Of File) indica que o analisador chegou ao final do arquivo enquanto ainda esperava mais conteúdo - um sinal clássico de chaves não fechadas

#### Verifique o código

Agora, vamos examinar `bad_syntax.nf` para entender o que está causando o erro:

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

Para o propósito deste exemplo, deixamos um comentário para você mostrar onde está o erro. A extensão Nextflow VSCode também deve estar dando algumas dicas sobre o que pode estar errado, colocando a chave incompatível em vermelho e destacando o fim prematuro do arquivo:

![Sintaxe ruim](img/bad_syntax.png)

**Estratégia de depuração para erros de colchetes:**

1. Use a correspondência de colchetes do VS Code (coloque o cursor ao lado de um colchete)
2. Verifique o painel Problemas para mensagens relacionadas a colchetes
3. Certifique-se de que cada `{` de abertura tenha um `}` de fechamento correspondente

#### Corrija o código

Substitua o comentário pela chave de fechamento ausente:

=== "Depois"

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
    }  // Add the missing closing brace

    workflow {

        // Create input channel
        input_ch = channel.of('sample1', 'sample2', 'sample3')

        // Call the process with the input channel
        PROCESS_FILES(input_ch)
    }
    ```

=== "Antes"

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

#### Execute o pipeline

Agora execute o fluxo de trabalho novamente para confirmar que funciona:

```bash
nextflow run bad_syntax.nf
```

??? success "Saída do comando"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `bad_syntax.nf` [insane_faggin] DSL2 - revision: 961938ee2b

    executor >  local (3)
    [48/cd7f54] PROCESS_FILES (1) | 3 of 3 ✔
    ```

### 1.2. Usando palavras-chave ou diretivas de processo incorretas

Outro erro de sintaxe comum é uma **definição de processo inválida**. Isso pode acontecer se você esquecer de definir blocos obrigatórios ou usar diretivas incorretas na definição do processo.

#### Execute o pipeline

```bash
nextflow run invalid_process.nf
```

??? failure "Saída do comando"

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

#### Verifique o código

O erro indica uma "Definição de processo inválida" e mostra o contexto em torno do problema. Olhando para as linhas 3-7, podemos ver `inputs:` na linha 4, que é o problema. Vamos examinar `invalid_process.nf`:

```groovy title="invalid_process.nf" hl_lines="4" linenums="1"
#!/usr/bin/env nextflow

process PROCESS_FILES {
    inputs:  // ERROR: Should be 'input' not 'inputs'
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

Olhando para a linha 4 no contexto do erro, podemos identificar o problema: estamos usando `inputs` em vez da diretiva correta `input`. A extensão Nextflow VSCode também sinalizará isso:

![Mensagem de processo inválido](img/invalid_process_message.png)

#### Corrija o código

Substitua a palavra-chave incorreta pela correta referenciando [a documentação](https://www.nextflow.io/docs/latest/process.html#):

=== "Depois"

    ```groovy title="invalid_process.nf" hl_lines="4" linenums="1"
    #!/usr/bin/env nextflow

    process PROCESS_FILES {
        input:  // Fixed: Changed 'inputs' to 'input'
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

=== "Antes"

    ```groovy title="invalid_process.nf" hl_lines="4" linenums="1"
    #!/usr/bin/env nextflow

    process PROCESS_FILES {
        inputs:  // ERROR: Should be 'input' not 'inputs'
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

#### Execute o pipeline

Agora execute o fluxo de trabalho novamente para confirmar que funciona:

```bash
nextflow run invalid_process.nf
```

??? success "Saída do comando"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `invalid_process.nf` [silly_fermi] DSL2 - revision: 961938ee2b

    executor >  local (3)
    [b7/76cd9d] PROCESS_FILES (2) | 3 of 3 ✔
    ```

### 1.3. Usando nomes de variáveis ruins

Os nomes de variáveis que você usa em seus blocos de script devem ser válidos, derivados de entradas ou de código groovy inserido antes do script. Mas quando você está lidando com complexidade no início do desenvolvimento do pipeline, é fácil cometer erros na nomeação de variáveis, e o Nextflow avisará rapidamente.

#### Execute o pipeline

```bash
nextflow run no_such_var.nf
```

??? failure "Saída do comando"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `no_such_var.nf` [gloomy_meninsky] DSL2 - revision: 0c4d3bc28c

    Error no_such_var.nf:17:39: `undefined_var` is not defined
    │  17 |     echo "Using undefined variable: ${undefined_var}" >> ${output_pref
    ╰     |                                       ^^^^^^^^^^^^^

    ERROR ~ Script compilation failed

     -- Check '.nextflow.log' file for details
    ```

O erro é capturado em tempo de compilação e aponta diretamente para a variável indefinida na linha 17, com um cursor indicando exatamente onde está o problema.

#### Verifique o código

Vamos examinar `no_such_var.nf`:

```groovy title="no_such_var.nf" hl_lines="17" linenums="1"
#!/usr/bin/env nextflow

process PROCESS_FILES {
    input:
    val sample_name

    output:
    path "${sample_name}_processed.txt"

    script:
    // Define variables in Groovy code before the script
    def output_prefix = "${sample_name}_processed"
    def timestamp = new Date().format("yyyy-MM-dd")

    """
    echo "Processing ${sample_name} on ${timestamp}" > ${output_prefix}.txt
    echo "Using undefined variable: ${undefined_var}" >> ${output_prefix}.txt  // ERROR: undefined_var not defined
    """
}

workflow {
    input_ch = channel.of('sample1', 'sample2', 'sample3')
    PROCESS_FILES(input_ch)
}
```

A mensagem de erro indica que a variável não é reconhecida no template do script, e lá está - você deve ser capaz de ver `${undefined_var}` usado no bloco de script, mas não definido em outro lugar.

#### Corrija o código

Se você receber um erro 'No such variable', pode corrigi-lo definindo a variável (corrigindo nomes de variáveis de entrada ou editando código groovy antes do script), ou removendo-a do bloco de script se não for necessária:

=== "Depois"

    ```groovy title="no_such_var.nf" hl_lines="15-17" linenums="1"
    #!/usr/bin/env nextflow

    process PROCESS_FILES {
        input:
        val sample_name

        output:
        path "${sample_name}_output.txt"

        script:
        // Define variables in Groovy code before the script
        def output_prefix = "${sample_name}_processed"
        def timestamp = new Date().format("yyyy-MM-dd")

        """
        echo "Processing ${sample_name} on ${timestamp}" > ${output_prefix}.txt
        """  // Removed the line with undefined_var
    }

    workflow {
        input_ch = channel.of('sample1', 'sample2', 'sample3')
        PROCESS_FILES(input_ch)
    }
    ```

=== "Antes"

    ```groovy title="no_such_var.nf" hl_lines="17" linenums="1"
    #!/usr/bin/env nextflow

    process PROCESS_FILES {
        input:
        val sample_name

        output:
        path "${sample_name}_output.txt"

        script:
        // Define variables in Groovy code before the script
        def output_prefix = "${sample_name}_processed"
        def timestamp = new Date().format("yyyy-MM-dd")

        """
        echo "Processing ${sample_name} on ${timestamp}" > ${output_prefix}.txt
        echo "Using undefined variable: ${undefined_var}" >> ${output_prefix}.txt  // ERROR: undefined_var not defined
        """
    }

    workflow {
        input_ch = channel.of('sample1', 'sample2', 'sample3')
        PROCESS_FILES(input_ch)
    }
    ```

#### Execute o pipeline

Agora execute o fluxo de trabalho novamente para confirmar que funciona:

```bash
nextflow run no_such_var.nf
```

??? success "Saída do comando"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `no_such_var.nf` [suspicious_venter] DSL2 - revision: 6ba490f7c5

    executor >  local (3)
    [21/237300] PROCESS_FILES (2) | 3 of 3 ✔
    ```

### 1.4. Uso incorreto de variáveis Bash

Começando no Nextflow, pode ser difícil entender a diferença entre variáveis Nextflow (Groovy) e Bash. Isso pode gerar outra forma do erro de variável ruim que aparece ao tentar usar variáveis no conteúdo Bash do bloco de script.

#### Execute o pipeline

```bash
nextflow run bad_bash_var.nf
```

??? failure "Saída do comando"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `bad_bash_var.nf` [infallible_mandelbrot] DSL2 - revision: 0853c11080

    Error bad_bash_var.nf:13:42: `prefix` is not defined
    │  13 |     echo "Processing ${sample_name}" > ${prefix}.txt
    ╰     |                                          ^^^^^^

    ERROR ~ Script compilation failed

     -- Check '.nextflow.log' file for details
    ```

#### Verifique o código

O erro aponta para a linha 13 onde `${prefix}` é usado. Vamos examinar `bad_bash_var.nf` para ver o que está causando o problema:

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
    echo "Processing ${sample_name}" > ${prefix}.txt  # ERROR: ${prefix} is Groovy syntax, not Bash
    """
}
```

Neste exemplo, estamos definindo a variável `prefix` em Bash, mas em um processo Nextflow a sintaxe `$` que usamos para referenciá-la (`${prefix}`) é interpretada como uma variável Groovy, não Bash. A variável não existe no contexto Groovy, então recebemos um erro 'no such variable'.

#### Corrija o código

Se você quiser usar uma variável Bash, deve escapar o cifrão assim:

=== "Depois"

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
        echo "Processing ${sample_name}" > \${prefix}.txt  # Fixed: Escaped the dollar sign
        """
    }

    workflow {
        input_ch = channel.of('sample1', 'sample2', 'sample3')
        PROCESS_FILES(input_ch)
    }
    ```

=== "Antes"

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
        echo "Processing ${sample_name}" > ${prefix}.txt  # ERROR: ${prefix} is Groovy syntax, not Bash
        """
    }
    ```

Isso diz ao Nextflow para interpretar isso como uma variável Bash.

#### Execute o pipeline

Agora execute o fluxo de trabalho novamente para confirmar que funciona:

```bash
nextflow run bad_bash_var.nf
```

??? success "Saída do comando"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `bad_bash_var.nf` [naughty_franklin] DSL2 - revision: 58c1c83709

    executor >  local (3)
    [4e/560285] PROCESS_FILES (2) | 3 of 3 ✔
    ```

!!! tip "Variáveis Groovy vs Bash"

    Para manipulações simples de variáveis como concatenação de strings ou operações de prefixo/sufixo, geralmente é mais legível usar variáveis Groovy na seção script em vez de variáveis Bash no bloco de script:

    ```groovy linenums="1"
    script:
    def output_prefix = "${sample_name}_processed"
    def output_file = "${output_prefix}.txt"
    """
    echo "Processing ${sample_name}" > ${output_file}
    """
    ```

    Esta abordagem evita a necessidade de escapar cifrões e torna o código mais fácil de ler e manter.

### 1.5. Declarações Fora do Bloco Workflow

A extensão Nextflow VSCode destaca problemas com a estrutura do código que causarão erros. Um exemplo comum é definir canais fora do bloco `workflow {}` - isso agora é aplicado como um erro de sintaxe.

#### Execute o pipeline

```bash
nextflow run badpractice_syntax.nf
```

??? failure "Saída do comando"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `badpractice_syntax.nf` [intergalactic_colden] DSL2 - revision: 5e4b291bde

    Error badpractice_syntax.nf:3:1: Statements cannot be mixed with script declarations -- move statements into a process or workflow
    │   3 | input_ch = channel.of('sample1', 'sample2', 'sample3')
    ╰     | ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

    ERROR ~ Script compilation failed

     -- Check '.nextflow.log' file for details
    ```

A mensagem de erro indica claramente o problema: declarações (como definições de canais) não podem ser misturadas com declarações de script fora de um bloco workflow ou process.

#### Verifique o código

Vamos examinar `badpractice_syntax.nf` para ver o que está causando o erro:

```groovy title="badpractice_syntax.nf" hl_lines="3" linenums="1"
#!/usr/bin/env nextflow

input_ch = channel.of('sample1', 'sample2', 'sample3')  // ERROR: Channel defined outside workflow

process PROCESS_FILES {
    input:
    val sample_name

    output:
    path "${sample_name}_processed.txt"

    script:
    // Define variables in Groovy code before the script
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

A extensão VSCode também destacará a variável `input_ch` como sendo definida fora do bloco workflow:

![Erro de sintaxe não letal](img/nonlethal.png)

#### Corrija o código

Mova a definição do canal para dentro do bloco workflow:

=== "Depois"

    ```groovy title="badpractice_syntax.nf" hl_lines="21" linenums="1"
    #!/usr/bin/env nextflow

    process PROCESS_FILES {
        input:
        val sample_name

        output:
        path "${sample_name}_processed.txt"

        script:
        // Define variables in Groovy code before the script
        def output_prefix = "${sample_name}_processed"
        def timestamp = new Date().format("yyyy-MM-dd")

        """
        echo "Processing ${sample_name} on ${timestamp}" > ${output_prefix}.txt
        """
    }

    workflow {
        input_ch = channel.of('sample1', 'sample2', 'sample3')  // Moved inside workflow block
        PROCESS_FILES(input_ch)
    }
    ```

=== "Antes"

    ```groovy title="badpractice_syntax.nf" hl_lines="3" linenums="1"
    #!/usr/bin/env nextflow

    input_ch = channel.of('sample1', 'sample2', 'sample3')  // ERROR: Channel defined outside workflow

    process PROCESS_FILES {
        input:
        val sample_name

        output:
        path "${sample_name}_processed.txt"

        script:
        // Define variables in Groovy code before the script
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

#### Execute o pipeline

Execute o fluxo de trabalho novamente para confirmar que a correção funciona:

```bash
nextflow run badpractice_syntax.nf
```

??? success "Saída do comando"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `badpractice_syntax.nf` [naughty_ochoa] DSL2 - revision: 5e4b291bde

    executor >  local (3)
    [6a/84a608] PROCESS_FILES (2) | 3 of 3 ✔
    ```

Mantenha seus canais de entrada definidos dentro do bloco workflow e, em geral, siga quaisquer outras recomendações que a extensão fizer.

### Conclusão

Você pode identificar e corrigir erros de sintaxe sistematicamente usando mensagens de erro do Nextflow e indicadores visuais da IDE. Erros de sintaxe comuns incluem chaves ausentes, palavras-chave de processo incorretas, variáveis indefinidas e uso inadequado de variáveis Bash vs. Nextflow. A extensão VSCode ajuda a capturar muitos desses antes do tempo de execução. Com essas habilidades de depuração de sintaxe em seu kit de ferramentas, você será capaz de resolver rapidamente os erros de sintaxe Nextflow mais comuns e passar para lidar com problemas de tempo de execução mais complexos.

### O que vem a seguir?

Aprenda a depurar erros de estrutura de canais mais complexos que ocorrem mesmo quando a sintaxe está correta.

---

## 2. Erros de Estrutura de Canais

Erros de estrutura de canais são mais sutis do que erros de sintaxe porque o código está sintaticamente correto, mas as formas dos dados não correspondem ao que os processos esperam. O Nextflow tentará executar o pipeline, mas pode descobrir que o número de entradas não corresponde ao que espera e falhar. Esses erros normalmente só aparecem em tempo de execução e requerem uma compreensão dos dados fluindo através do seu fluxo de trabalho.

!!! tip "Depurando Canais com `.view()`"

    Ao longo desta seção, lembre-se de que você pode usar o operador `.view()` para inspecionar o conteúdo do canal em qualquer ponto do seu fluxo de trabalho. Esta é uma das ferramentas de depuração mais poderosas para entender problemas de estrutura de canais. Exploraremos esta técnica em detalhes na seção 2.4, mas sinta-se à vontade para usá-la enquanto trabalha nos exemplos.

    ```groovy
    my_channel.view()  // Shows what's flowing through the channel
    ```

### 2.1. Número Errado de Canais de Entrada

Este erro ocorre quando você passa um número diferente de canais do que um processo espera.

#### Execute o pipeline

```bash
nextflow run bad_number_inputs.nf
```

??? failure "Saída do comando"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `bad_number_inputs.nf` [happy_swartz] DSL2 - revision: d83e58dcd3

    Error bad_number_inputs.nf:23:5: Incorrect number of call arguments, expected 1 but received 2
    │  23 |     PROCESS_FILES(samples_ch, files_ch)
    ╰     |     ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

    ERROR ~ Script compilation failed

     -- Check '.nextflow.log' file for details
    ```

#### Verifique o código

A mensagem de erro afirma claramente que a chamada esperava 1 argumento mas recebeu 2, e aponta para a linha 23. Vamos examinar `bad_number_inputs.nf`:

```groovy title="bad_number_inputs.nf" hl_lines="5 23" linenums="1"
#!/usr/bin/env nextflow

process PROCESS_FILES {
    input:
        val sample_name  // Process expects only 1 input

    output:
        path "${sample_name}_output.txt"

    script:
    """
    echo "Processing ${sample_name}" > ${sample_name}_output.txt
    """
}

workflow {

    // Create two separate channels
    samples_ch = channel.of('sample1', 'sample2', 'sample3')
    files_ch = channel.of('file1.txt', 'file2.txt', 'file3.txt')

    // ERROR: Passing 2 channels but process expects only 1
    PROCESS_FILES(samples_ch, files_ch)
}
```

Você deve ver a chamada `PROCESS_FILES` incompatível, fornecendo múltiplos canais de entrada quando o processo define apenas um. A extensão VSCode também sublinhará a chamada do processo em vermelho e fornecerá uma mensagem de diagnóstico quando você passar o mouse:

![Mensagem de número incorreto de argumentos](img/incorrect_num_args.png)

#### Corrija o código

Para este exemplo específico, o processo espera um único canal e não requer o segundo canal, então podemos corrigi-lo passando apenas o canal `samples_ch`:

=== "Depois"

    ```groovy title="bad_number_inputs.nf" hl_lines="23" linenums="1"
    #!/usr/bin/env nextflow

    process PROCESS_FILES {
        input:
            val sample_name  // Process expects only 1 input

        output:
            path "${sample_name}_output.txt"

        script:
        """
        echo "Processing ${sample_name}" > ${sample_name}_output.txt
        """
    }

    workflow {

        // Create two separate channels
        samples_ch = channel.of('sample1', 'sample2', 'sample3')
        files_ch = channel.of('file1.txt', 'file2.txt', 'file3.txt')

        // Fixed: Pass only the channel the process expects
        PROCESS_FILES(samples_ch)
    }
    ```

=== "Antes"

    ```groovy title="bad_number_inputs.nf" hl_lines="5 23" linenums="1"
    #!/usr/bin/env nextflow

    process PROCESS_FILES {
        input:
            val sample_name  // Process expects only 1 input

        output:
            path "${sample_name}_output.txt"

        script:
        """
        echo "Processing ${sample_name}" > ${sample_name}_output.txt
        """
    }

    workflow {

        // Create two separate channels
        samples_ch = channel.of('sample1', 'sample2', 'sample3')
        files_ch = channel.of('file1.txt', 'file2.txt', 'file3.txt')

        // ERROR: Passing 2 channels but process expects only 1
        PROCESS_FILES(samples_ch, files_ch)
    }
    ```

#### Execute o pipeline

```bash
nextflow run bad_number_inputs.nf
```

??? success "Saída do comando"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `bad_number_inputs.nf` [big_euler] DSL2 - revision: e302bd87be

    executor >  local (3)
    [48/497f7b] PROCESS_FILES (3) | 3 of 3 ✔
    ```

Mais comumente do que este exemplo, você pode adicionar entradas adicionais a um processo e esquecer de atualizar a chamada do fluxo de trabalho de acordo, o que pode levar a este tipo de erro. Felizmente, este é um dos erros mais fáceis de entender e corrigir, pois a mensagem de erro é bastante clara sobre a incompatibilidade.

### 2.2. Esgotamento de Canal (Processo Executa Menos Vezes do que o Esperado)

Alguns erros de estrutura de canais são muito mais sutis e não produzem erros. Provavelmente o mais comum deles reflete um desafio que novos usuários do Nextflow enfrentam ao entender que canais de fila podem ser esgotados e ficar sem itens, significando que o fluxo de trabalho termina prematuramente.

#### Execute o pipeline

```bash
nextflow run exhausted.nf
```

??? success "Saída do comando"

```console title="Exhausted channel output"
 N E X T F L O W   ~  version 25.10.2

Launching `exhausted.nf` [extravagant_gauss] DSL2 - revision: 08cff7ba2a

executor >  local (1)
[bd/f61fff] PROCESS_FILES (1) [100%] 1 of 1 ✔
```

Este fluxo de trabalho é concluído sem erro, mas processa apenas uma única amostra!

#### Verifique o código

Vamos examinar `exhausted.nf` para ver se está certo:

```groovy title="exhausted.nf" hl_lines="23 24" linenums="1"
#!/usr/bin/env nextflow

process PROCESS_FILES {
    input:
    val reference
    val sample_name

    output:
    path "${output_prefix}.txt"

    script:
    // Define variables in Groovy code before the script
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

O processo executa apenas uma vez em vez de três vezes porque o canal `reference_ch` é um canal de fila que se esgota após a primeira execução do processo. Quando um canal se esgota, o processo inteiro para, mesmo que outros canais ainda tenham itens.

Este é um padrão comum onde você tem um único arquivo de referência que precisa ser reutilizado em várias amostras. A solução é converter o canal de referência em um canal de valor que pode ser reutilizado indefinidamente.

#### Corrija o código

Existem algumas maneiras de abordar isso dependendo de quantos arquivos são afetados.

**Opção 1**: Você tem um único arquivo de referência que está reutilizando muito. Você pode simplesmente criar um tipo de canal de valor, que pode ser usado repetidamente. Existem três maneiras de fazer isso:

**1a** Use `channel.value()`:

```groovy title="exhausted.nf (fixed - Option 1a)" hl_lines="2" linenums="21"
workflow {
    reference_ch = channel.value('baseline_reference')  // Value channel can be reused
    input_ch = channel.of('sample1', 'sample2', 'sample3')

    PROCESS_FILES(reference_ch, input_ch)
}
```

**1b** Use o operador `first()` [operator](https://www.nextflow.io/docs/latest/reference/operator.html#first):

```groovy title="exhausted.nf (fixed - Option 1b)" hl_lines="2" linenums="21"
workflow {
    reference_ch = channel.of('baseline_reference').first()  // Convert to value channel
    input_ch = channel.of('sample1', 'sample2', 'sample3')

    PROCESS_FILES(reference_ch, input_ch)
}
```

**1c.** Use o operador `collect()` [operator](https://www.nextflow.io/docs/latest/reference/operator.html#collect):

```groovy title="exhausted.nf (fixed - Option 1c)" hl_lines="2" linenums="21"
workflow {
    reference_ch = channel.of('baseline_reference').collect()  // Convert to value channel
    input_ch = channel.of('sample1', 'sample2', 'sample3')

    PROCESS_FILES(reference_ch, input_ch)
}
```

**Opção 2**: Em cenários mais complexos, talvez onde você tenha múltiplos arquivos de referência para todas as amostras no canal de amostras, você pode usar o operador `combine` para criar um novo canal que combina os dois canais em tuplas:

```groovy title="exhausted.nf (fixed - Option 2)" hl_lines="4" linenums="21"
workflow {
    reference_ch = channel.of('baseline_reference','other_reference')
    input_ch = channel.of('sample1', 'sample2', 'sample3')
    combined_ch = reference_ch.combine(input_ch)  // Creates cartesian product

    PROCESS_FILES(combined_ch)
}
```

O operador `.combine()` gera um produto cartesiano dos dois canais, então cada item em `reference_ch` será pareado com cada item em `input_ch`. Isso permite que o processo execute para cada amostra enquanto ainda usa a referência.

Isso requer que a entrada do processo seja ajustada. Em nosso exemplo, o início da definição do processo precisaria ser ajustado da seguinte forma:

```groovy title="exhausted.nf (fixed - Option 2)" hl_lines="5" linenums="1"
#!/usr/bin/env nextflow

process PROCESS_FILES {
    input:
        tuple val(reference), val(sample_name)
```

Esta abordagem pode não ser adequada em todas as situações.

#### Execute o pipeline

Tente uma das correções acima e execute o fluxo de trabalho novamente:

```bash
nextflow run exhausted.nf
```

??? success "Saída do comando"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `exhausted.nf` [maniac_leavitt] DSL2 - revision: f372a56a7d

    executor >  local (3)
    [80/0779e9] PROCESS_FILES (3) | 3 of 3 ✔
    ```

Você deve agora ver todas as três amostras sendo processadas em vez de apenas uma.

### 2.3. Estrutura de Conteúdo de Canal Errada

Quando os fluxos de trabalho atingem um certo nível de complexidade, pode ser um pouco difícil acompanhar as estruturas internas de cada canal, e as pessoas comumente geram incompatibilidades entre o que o processo espera e o que o canal realmente contém. Isso é mais sutil do que o problema que discutimos anteriormente, onde o número de canais estava incorreto. Neste caso, você pode ter o número correto de canais de entrada, mas a estrutura interna de um ou mais desses canais não corresponde ao que o processo espera.

#### Execute o pipeline

```bash
nextflow run bad_channel_shape.nf
```

??? failure "Saída do comando"

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

#### Verifique o código

Os colchetes na mensagem de erro fornecem a pista aqui - o processo está tratando a tupla como um único valor, o que não é o que queremos. Vamos examinar `bad_channel_shape.nf`:

```groovy title="bad_channel_shape.nf" hl_lines="5 20-22" linenums="1"
#!/usr/bin/env nextflow

process PROCESS_FILES {
    input:
        val sample_name  // Expects single value, gets tuple

    output:
        path "${sample_name}_output.txt"

    script:
    """
    echo "Processing ${sample_name}" > ${sample_name}_output.txt
    """
}

workflow {

    // Channel emits tuples, but process expects single values
    input_ch = channel.of(
      ['sample1', 'file1.txt'],
      ['sample2', 'file2.txt'],
      ['sample3', 'file3.txt']
    )
    PROCESS_FILES(input_ch)
}
```

Você pode ver que estamos gerando um canal composto de tuplas: `['sample1', 'file1.txt']`, mas o processo espera um único valor, `val sample_name`. O comando executado mostra que o processo está tentando criar um arquivo chamado `[sample3, file3.txt]_output.txt`, que não é a saída pretendida.

#### Corrija o código

Para corrigir isso, se o processo requer ambas as entradas, poderíamos ajustar o processo para aceitar uma tupla:

=== "Opção 1: Aceitar tupla no processo"

    === "Depois"

        ```groovy title="bad_channel_shape.nf" hl_lines="5"  linenums="1"
        #!/usr/bin/env nextflow

        process PROCESS_FILES {
            input:
                tuple val(sample_name), val(file_name)  // Fixed: Accept tuple

            output:
                path "${sample_name}_output.txt"

            script:
            """
            echo "Processing ${sample_name}" > ${sample_name}_output.txt
            """
        }

        workflow {

            // Channel emits tuples, but process expects single values
            input_ch = channel.of(
              ['sample1', 'file1.txt'],
              ['sample2', 'file2.txt'],
              ['sample3', 'file3.txt']
            )
            PROCESS_FILES(input_ch)
        }
        ```

    === "Antes"

        ```groovy title="bad_channel_shape.nf" hl_lines="5" linenums="1"
        #!/usr/bin/env nextflow

        process PROCESS_FILES {
            input:
                val sample_name  // Expects single value, gets tuple

            output:
                path "${sample_name}_output.txt"

            script:
            """
            echo "Processing ${sample_name}" > ${sample_name}_output.txt
            """
        }

        workflow {

            // Channel emits tuples, but process expects single values
            input_ch = channel.of(
              ['sample1', 'file1.txt'],
              ['sample2', 'file2.txt'],
              ['sample3', 'file3.txt']
            )
            PROCESS_FILES(input_ch)
        }
        ```

=== "Opção 2: Extrair primeiro elemento"

    === "Depois"

        ```groovy title="bad_channel_shape.nf" hl_lines="9" linenums="16"
        workflow {

            // Channel emits tuples, but process expects single values
            input_ch = channel.of(
              ['sample1', 'file1.txt'],
              ['sample2', 'file2.txt'],
              ['sample3', 'file3.txt']
            )
            PROCESS_FILES(input_ch.map { it[0] })  // Fixed: Extract first element
        }
        ```

    === "Antes"

        ```groovy title="bad_channel_shape.nf" hl_lines="9" linenums="16"
        workflow {

            // Channel emits tuples, but process expects single values
            input_ch = channel.of(
              ['sample1', 'file1.txt'],
              ['sample2', 'file2.txt'],
              ['sample3', 'file3.txt']
            )
            PROCESS_FILES(input_ch)
        }
        ```

#### Execute o pipeline

Escolha uma das soluções e execute novamente o fluxo de trabalho:

```bash
nextflow run bad_channel_shape.nf
```

??? success "Saída do comando"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `bad_channel_shape.nf` [clever_thompson] DSL2 - revision: 8cbcae3746

    executor >  local (3)
    [bb/80a958] PROCESS_FILES (2) | 3 of 3 ✔
    ```

### 2.4. Técnicas de Depuração de Canais

#### Usando `.view()` para Inspeção de Canais

A ferramenta de depuração mais poderosa para canais é o operador `.view()`. Com `.view()`, você pode entender a forma de seus canais em todos os estágios para ajudar na depuração.

#### Execute o pipeline

Execute `bad_channel_shape_viewed.nf` para ver isso em ação:

```bash
nextflow run bad_channel_shape_viewed.nf
```

??? success "Saída do comando"

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

#### Verifique o código

Vamos examinar `bad_channel_shape_viewed.nf` para ver como `.view()` é usado:

```groovy title="bad_channel_shape_viewed.nf" linenums="16" hl_lines="9 11"
workflow {

    // Channel emits tuples, but process expects single values
    input_ch = channel.of(
      ['sample1', 'file1.txt'],
      ['sample2', 'file2.txt'],
      ['sample3', 'file3.txt']
    )
    .view { "Channel content: $it" }  // Debug: Show original channel content
    .map { tuple -> tuple[0] }        // Transform: Extract first element
    .view { "After mapping: $it" }    // Debug: Show transformed channel content

    PROCESS_FILES(input_ch)
}
```

#### Corrija o código

Para evitar que você use operações `.view()` excessivamente no futuro para entender o conteúdo do canal, é aconselhável adicionar alguns comentários para ajudar:

```groovy title="bad_channel_shape_viewed.nf (with comments)" linenums="16" hl_lines="8 9"
workflow {

    // Channel emits tuples, but process expects single values
    input_ch = channel.of(
            ['sample1', 'file1.txt'],
            ['sample2', 'file2.txt'],
            ['sample3', 'file3.txt'],
        ) // [sample_name, file_name]
        .map { tuple -> tuple[0] } // sample_name

    PROCESS_FILES(input_ch)
}
```

Isso se tornará mais importante à medida que seus fluxos de trabalho crescem em complexidade e a estrutura do canal se torna mais opaca.

#### Execute o pipeline

```bash
nextflow run bad_channel_shape_viewed.nf
```

??? success "Saída do comando"

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

### Conclusão

Muitos erros de estrutura de canais podem ser criados com sintaxe Nextflow válida. Você pode depurar erros de estrutura de canais entendendo o fluxo de dados, usando operadores `.view()` para inspeção e reconhecendo padrões de mensagens de erro como colchetes indicando estruturas de tuplas inesperadas.

### O que vem a seguir?

Aprenda sobre erros criados por definições de processos.

---

## 3. Erros de Estrutura de Processos

A maioria dos erros que você encontrará relacionados a processos estará relacionada a erros que você cometeu ao formar o comando, ou a problemas relacionados ao software subjacente. Dito isso, semelhante aos problemas de canais acima, você pode cometer erros na definição do processo que não se qualificam como erros de sintaxe, mas que causarão erros em tempo de execução.

### 3.1. Arquivos de Saída Ausentes

Um erro comum ao escrever processos é fazer algo que gera uma incompatibilidade entre o que o processo espera e o que é gerado.

#### Execute o pipeline

```bash
nextflow run missing_output.nf
```

??? failure "Saída do comando"

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

#### Verifique o código

A mensagem de erro indica que o processo esperava produzir um arquivo de saída chamado `sample3.txt`, mas o script realmente cria `sample3_output.txt`. Vamos examinar a definição do processo em `missing_output.nf`:

```groovy title="missing_output.nf" linenums="3" hl_lines="6 10"
process PROCESS_FILES {
    input:
    val sample_name

    output:
    path "${sample_name}.txt"  // Expects: sample3.txt

    script:
    """
    echo "Processing ${sample_name}" > ${sample_name}_output.txt  // Creates: sample3_output.txt
    """
}
```

Você deve ver que há uma incompatibilidade entre o nome do arquivo de saída no bloco `output:` e o usado no script. Essa incompatibilidade faz com que o processo falhe. Se você encontrar esse tipo de erro, volte e verifique se as saídas correspondem entre sua definição de processo e seu bloco de saída.

Se o problema ainda não estiver claro, verifique o próprio diretório de trabalho para identificar os arquivos de saída reais criados:

```bash
❯ ls -h work/02/9604d49fb8200a74d737c72a6c98ed
sample3_output.txt
```

Para este exemplo, isso destacaria para nós que um sufixo `_output` está sendo incorporado ao nome do arquivo de saída, contrário à nossa definição `output:`.

#### Corrija o código

Corrija a incompatibilidade tornando o nome do arquivo de saída consistente:

=== "Depois"

    ```groovy title="missing_output.nf" hl_lines="6 10" linenums="3"
    process PROCESS_FILES {
        input:
        val sample_name

        output:
        path "${sample_name}_output.txt"  // Fixed: Match the script output

        script:
        """
        echo "Processing ${sample_name}" > ${sample_name}_output.txt
        """
    }
    ```

=== "Antes"

    ```groovy title="missing_output.nf" hl_lines="6 10" linenums="3"
    process PROCESS_FILES {
        input:
        val sample_name

        output:
        path "${sample_name}.txt"  // Expects: sample3.txt

        script:
        """
        echo "Processing ${sample_name}" > ${sample_name}_output.txt  // Creates: sample3_output.txt
        """
    }
    ```

#### Execute o pipeline

```bash
nextflow run missing_output.nf
```

??? success "Saída do comando"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `missing_output.nf` [elated_hamilton] DSL2 - revision: 961938ee2b

    executor >  local (3)
    [16/1c437c] PROCESS_FILES (3) | 3 of 3 ✔
    ```

### 3.2. Software ausente

Outra classe de erros ocorre devido a erros no provisionamento de software. `missing_software.nf` é um fluxo de trabalho sintaticamente válido, mas depende de algum software externo para fornecer o comando `cowpy` que usa.

#### Execute o pipeline

```bash
nextflow run missing_software.nf
```

??? failure "Saída do comando"

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

O processo não tem acesso ao comando que estamos especificando. Às vezes isso ocorre porque um script está presente no diretório `bin` do fluxo de trabalho, mas não foi tornado executável. Outras vezes é porque o software não está instalado no contêiner ou ambiente onde o fluxo de trabalho está sendo executado.

#### Verifique o código

Fique atento a esse código de saída `127` - ele diz exatamente qual é o problema. Vamos examinar `missing_software.nf`:

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

#### Corrija o código

Fomos um pouco desonestos aqui, e na verdade não há nada de errado com o código. Só precisamos especificar a configuração necessária para executar o processo de forma que ele tenha acesso ao comando em questão. Neste caso, o processo tem uma definição de contêiner, então tudo o que precisamos fazer é executar o fluxo de trabalho com Docker habilitado.

#### Execute o pipeline

Configuramos um perfil Docker para você em `nextflow.config`, então você pode executar o fluxo de trabalho com:

```bash
nextflow run missing_software.nf -profile docker
```

??? success "Saída do comando"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `missing_software.nf` [awesome_stonebraker] DSL2 - revision: 0296d12839

    executor >  local (3)
    [38/ab20d1] PROCESS_FILES (1) | 3 of 3 ✔
    ```

!!! note

    Para saber mais sobre como o Nextflow usa contêineres, veja [Hello Nextflow](../hello_nextflow/05_hello_containers.md)

### 3.3. Configuração de recursos ruim

No uso em produção, você estará configurando recursos em seus processos. Por exemplo, `memory` define a quantidade máxima de memória disponível para seu processo, e se o processo exceder isso, seu agendador normalmente matará o processo e retornará um código de saída de `137`. Não podemos demonstrar isso aqui porque estamos usando o executor `local`, mas podemos mostrar algo semelhante com `time`.

#### Execute o pipeline

`bad_resources.nf` tem configuração de processo com um limite de tempo irrealista de 1 milissegundo:

```bash
nextflow run bad_resources.nf -profile docker
```

??? failure "Saída do comando"

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

#### Verifique o código

Vamos examinar `bad_resources.nf`:

```groovy title="bad_resources.nf" linenums="3" hl_lines="3"
process PROCESS_FILES {

    time '1 ms'  // ERROR: Unrealistic time limit

    input:
    val sample_name

    output:
    path "${sample_name}_output.txt"

    script:
    """
    sleep 1  // Takes 1 second, but time limit is 1ms
    cowpy ${sample_name} > ${sample_name}_output.txt
    """
}
```

Sabemos que o processo levará mais de um segundo (adicionamos um sleep lá para garantir), mas o processo está configurado para expirar após 1 milissegundo. Alguém foi um pouco irrealista com sua configuração!

#### Corrija o código

Aumente o limite de tempo para um valor realista:

=== "Depois"

    ```groovy title="bad_resources.nf" hl_lines="3" linenums="3"
    process PROCESS_FILES {

        time '100 s'  // Fixed: Realistic time limit

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

=== "Antes"

    ```groovy title="bad_resources.nf" hl_lines="3" linenums="3"
    process PROCESS_FILES {

        time '1 ms'  // ERROR: Unrealistic time limit

        input:
        val sample_name

        output:
        path "${sample_name}_output.txt"

        script:
        """
        sleep 1  // Takes 1 second, but time limit is 1ms
        cowpy ${sample_name} > ${sample_name}_output.txt
        """
    }
    ```

#### Execute o pipeline

```bash
nextflow run bad_resources.nf -profile docker
```

??? success "Saída do comando"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `bad_resources.nf` [friendly_mcclintock] DSL2 - revision: 381567d2c1

    executor >  local (3)
    [c2/9b4c41] PROCESS_FILES (3) | 3 of 3 ✔
    ```

Se você se certificar de ler suas mensagens de erro, falhas como esta não devem confundi-lo por muito tempo. Mas certifique-se de entender os requisitos de recursos dos comandos que está executando para que possa configurar suas diretivas de recursos apropriadamente.

### 3.4. Técnicas de Depuração de Processos

Quando os processos falham ou se comportam inesperadamente, você precisa de técnicas sistemáticas para investigar o que deu errado. O diretório de trabalho contém todas as informações necessárias para depurar a execução do processo.

#### Usando Inspeção do Diretório de Trabalho

A ferramenta de depuração mais poderosa para processos é examinar o diretório de trabalho. Quando um processo falha, o Nextflow cria um diretório de trabalho para essa execução específica do processo contendo todos os arquivos necessários para entender o que aconteceu.

#### Execute o pipeline

Vamos usar o exemplo `missing_output.nf` de antes para demonstrar a inspeção do diretório de trabalho (regenere uma incompatibilidade de nomenclatura de saída se precisar):

```bash
nextflow run missing_output.nf
```

??? failure "Saída do comando"

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

#### Verifique o diretório de trabalho

Quando você recebe este erro, o diretório de trabalho contém todas as informações de depuração. Encontre o caminho do diretório de trabalho na mensagem de erro e examine seu conteúdo:

```bash
# Find the work directory from the error message
ls work/02/9604d49fb8200a74d737c72a6c98ed/
```

Você pode então examinar os arquivos-chave:

##### Verifique o Script de Comando

O arquivo `.command.sh` mostra exatamente qual comando foi executado:

```bash
# View the executed command
cat work/02/9604d49fb8200a74d737c72a6c98ed/.command.sh
```

Isso revela:

- **Substituição de variáveis**: Se as variáveis Nextflow foram expandidas corretamente
- **Caminhos de arquivo**: Se os arquivos de entrada foram localizados corretamente
- **Estrutura do comando**: Se a sintaxe do script está correta

Problemas comuns a procurar:

- **Aspas ausentes**: Variáveis contendo espaços precisam de aspas adequadas
- **Caminhos de arquivo errados**: Arquivos de entrada que não existem ou estão em locais errados
- **Nomes de variáveis incorretos**: Erros de digitação em referências de variáveis
- **Configuração de ambiente ausente**: Comandos que dependem de ambientes específicos

##### Verifique a Saída de Erro

O arquivo `.command.err` contém as mensagens de erro reais:

```bash
# View error output
cat work/02/9604d49fb8200a74d737c72a6c98ed/.command.err
```

Este arquivo mostrará:

- **Códigos de saída**: 127 (comando não encontrado), 137 (morto), etc.
- **Erros de permissão**: Problemas de acesso a arquivos
- **Erros de software**: Mensagens de erro específicas da aplicação
- **Erros de recursos**: Limite de memória/tempo excedido

##### Verifique a Saída Padrão

O arquivo `.command.out` mostra o que seu comando produziu:

```bash
# View standard output
cat work/02/9604d49fb8200a74d737c72a6c98ed/.command.out
```

Isso ajuda a verificar:

- **Saída esperada**: Se o comando produziu os resultados corretos
- **Execução parcial**: Se o comando começou mas falhou no meio do caminho
- **Informações de depuração**: Qualquer saída de diagnóstico do seu script

##### Verifique o Código de Saída

O arquivo `.exitcode` contém o código de saída do processo:

```bash
# View exit code
cat work/*/*/.exitcode
```

Códigos de saída comuns e seus significados:

- **Código de saída 127**: Comando não encontrado - verifique a instalação do software
- **Código de saída 137**: Processo morto - verifique limites de memória/tempo

##### Verifique a Existência de Arquivos

Quando os processos falham devido a arquivos de saída ausentes, verifique quais arquivos foram realmente criados:

```bash
# List all files in the work directory
ls -la work/02/9604d49fb8200a74d737c72a6c98ed/
```

Isso ajuda a identificar:

- **Incompatibilidades de nomenclatura de arquivos**: Arquivos de saída com nomes diferentes do esperado
- **Problemas de permissão**: Arquivos que não puderam ser criados
- **Problemas de caminho**: Arquivos criados em diretórios errados

Em nosso exemplo anterior, isso confirmou para nós que enquanto nosso `sample3.txt` esperado não estava presente, `sample3_output.txt` estava:

```bash
❯ ls -h work/02/9604d49fb8200a74d737c72a6c98ed
sample3_output.txt
```

### Conclusão

A depuração de processos requer examinar diretórios de trabalho para entender o que deu errado. Os arquivos-chave incluem `.command.sh` (o script executado), `.command.err` (mensagens de erro) e `.command.out` (saída padrão). Códigos de saída como 127 (comando não encontrado) e 137 (processo morto) fornecem pistas de diagnóstico imediatas sobre o tipo de falha.

### O que vem a seguir?

Aprenda sobre as ferramentas de depuração integradas do Nextflow e abordagens sistemáticas para solução de problemas.

---

## 4. Ferramentas de Depuração Integradas e Técnicas Avançadas

O Nextflow fornece várias ferramentas integradas poderosas para depuração e análise de execução de fluxo de trabalho. Essas ferramentas ajudam você a entender o que deu errado, onde deu errado e como corrigi-lo de forma eficiente.

### 4.1. Saída de Processo em Tempo Real

Às vezes você precisa ver o que está acontecendo dentro dos processos em execução. Você pode habilitar a saída de processo em tempo real, que mostra exatamente o que cada tarefa está fazendo enquanto executa.

#### Execute o pipeline

`bad_channel_shape_viewed.nf` de nossos exemplos anteriores imprimiu o conteúdo do canal usando `.view()`, mas também podemos usar a diretiva `debug` para ecoar variáveis de dentro do próprio processo, o que demonstramos em `bad_channel_shape_viewed_debug.nf`. Execute o fluxo de trabalho:

```bash
nextflow run bad_channel_shape_viewed_debug.nf
```

??? success "Saída do comando"

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

#### Verifique o código

Vamos examinar `bad_channel_shape_viewed_debug.nf` para ver como a diretiva `debug` funciona:

```groovy title="bad_channel_shape_viewed_debug.nf" linenums="3" hl_lines="2"
process PROCESS_FILES {
    debug true  // Enable real-time output

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

A diretiva `debug` pode ser uma maneira rápida e conveniente de entender o ambiente de um processo.

### 4.2. Modo de Visualização

Às vezes você quer detectar problemas antes que qualquer processo seja executado. O Nextflow fornece uma flag para este tipo de depuração proativa: `-preview`.

#### Execute o pipeline

O modo de visualização permite testar a lógica do fluxo de trabalho sem executar comandos. Isso pode ser bastante útil para verificar rapidamente a estrutura do seu fluxo de trabalho e garantir que os processos estejam conectados corretamente sem executar nenhum comando real.

!!! note

    Se você corrigiu `bad_syntax.nf` anteriormente, reintroduza o erro de sintaxe removendo a chave de fechamento após o bloco de script antes de executar este comando.

Execute este comando:

```bash
nextflow run bad_syntax.nf -preview
```

??? failure "Saída do comando"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `bad_syntax.nf` [magical_mercator] DSL2 - revision: 550b9a8873

    Error bad_syntax.nf:24:1: Unexpected input: '<EOF>'

    ERROR ~ Script compilation failed

     -- Check '.nextflow.log' file for details
    ```

O modo de visualização é particularmente útil para detectar erros de sintaxe cedo sem executar nenhum processo. Ele valida a estrutura do fluxo de trabalho e as conexões de processos antes da execução.

### 4.3. Execução Stub para Teste de Lógica

Às vezes os erros são difíceis de depurar porque os comandos demoram muito, requerem software especial ou falham por razões complexas. A execução stub permite testar a lógica do fluxo de trabalho sem executar os comandos reais.

#### Execute o pipeline

Quando você está desenvolvendo um processo Nextflow, pode usar a diretiva `stub` para definir comandos 'fictícios' que geram saídas da forma correta sem executar o comando real. Esta abordagem é particularmente valiosa quando você quer verificar se a lógica do seu fluxo de trabalho está correta antes de lidar com as complexidades do software real.

Por exemplo, lembra do nosso `missing_software.nf` de antes? Aquele onde tínhamos software ausente que impediu a execução do fluxo de trabalho até adicionarmos `-profile docker`? `missing_software_with_stub.nf` é um fluxo de trabalho muito semelhante. Se o executarmos da mesma forma, geraremos o mesmo erro:

```bash
nextflow run missing_software_with_stub.nf
```

??? failure "Saída do comando"

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

No entanto, este fluxo de trabalho não produzirá erros se o executarmos com `-stub-run`, mesmo sem o perfil `docker`:

```bash
nextflow run missing_software_with_stub.nf -stub-run
```

??? success "Saída do comando"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `missing_software_with_stub.nf` [astonishing_shockley] DSL2 - revision: f1f4f05d7d

    executor >  local (3)
    [b5/2517a3] PROCESS_FILES (3) | 3 of 3 ✔
    ```

#### Verifique o código

Vamos examinar `missing_software_with_stub.nf`:

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

Relativo a `missing_software.nf`, este processo tem uma diretiva `stub:` especificando um comando a ser usado em vez do especificado em `script:`, no caso de o Nextflow ser executado em modo stub.

O comando `touch` que estamos usando aqui não depende de nenhum software ou entradas apropriadas, e será executado em todas as situações, permitindo-nos depurar a lógica do fluxo de trabalho sem nos preocuparmos com os internos do processo.

**A execução stub ajuda a depurar:**

- Estrutura de canais e fluxo de dados
- Conexões e dependências de processos
- Propagação de parâmetros
- Lógica de fluxo de trabalho sem dependências de software

### 4.4. Abordagem Sistemática de Depuração

Agora que você aprendeu técnicas individuais de depuração - de arquivos de rastreamento e diretórios de trabalho ao modo de visualização, execução stub e monitoramento de recursos - vamos uni-las em uma metodologia sistemática. Ter uma abordagem estruturada evita que você fique sobrecarregado por erros complexos e garante que você não perca pistas importantes.

Esta metodologia combina todas as ferramentas que cobrimos em um fluxo de trabalho eficiente:

**Método de Depuração em Quatro Fases:**

**Fase 1: Resolução de Erros de Sintaxe (5 minutos)**

1. Verifique sublinhados vermelhos no VSCode ou sua IDE
2. Execute `nextflow run workflow.nf -preview` para identificar problemas de sintaxe
3. Corrija todos os erros de sintaxe (chaves ausentes, vírgulas finais, etc.)
4. Certifique-se de que o fluxo de trabalho seja analisado com sucesso antes de prosseguir

**Fase 2: Avaliação Rápida (5 minutos)**

1. Leia as mensagens de erro de tempo de execução cuidadosamente
2. Verifique se é um erro de tempo de execução, lógica ou recursos
3. Use o modo de visualização para testar a lógica básica do fluxo de trabalho

**Fase 3: Investigação Detalhada (15-30 minutos)**

1. Encontre o diretório de trabalho da tarefa que falhou
2. Examine arquivos de log
3. Adicione operadores `.view()` para inspecionar canais
4. Use `-stub-run` para testar a lógica do fluxo de trabalho sem execução

**Fase 4: Corrigir e Validar (15 minutos)**

1. Faça correções mínimas direcionadas
2. Teste com resume: `nextflow run workflow.nf -resume`
3. Verifique a execução completa do fluxo de trabalho

!!! tip "Usando Resume para Depuração Eficiente"

    Uma vez que você identificou um problema, precisa de uma maneira eficiente de testar suas correções sem desperdiçar tempo reexecutando partes bem-sucedidas do seu fluxo de trabalho. A funcionalidade `-resume` do Nextflow é inestimável para depuração.

    Você terá encontrado `-resume` se trabalhou em [Hello Nextflow](../hello_nextflow/), e é importante que você faça bom uso dele ao depurar para economizar tempo de espera enquanto os processos antes do seu processo problemático são executados.

    **Estratégia de depuração com resume:**

    1. Execute o fluxo de trabalho até a falha
    2. Examine o diretório de trabalho da tarefa que falhou
    3. Corrija o problema específico
    4. Retome para testar apenas a correção
    5. Repita até que o fluxo de trabalho seja concluído

#### Perfil de Configuração de Depuração

Para tornar esta abordagem sistemática ainda mais eficiente, você pode criar uma configuração de depuração dedicada que habilita automaticamente todas as ferramentas necessárias:

```groovy title="nextflow.config (debug profile)" linenums="1"
profiles {
    debug {
        process {
            debug = true
            cleanup = false

            // Conservative resources for debugging
            maxForks = 1
            memory = '2.GB'
            cpus = 1
        }
    }
}
```

Então você pode executar o pipeline com este perfil habilitado:

```bash
nextflow run workflow.nf -profile debug
```

Este perfil habilita saída em tempo real, preserva diretórios de trabalho e limita a paralelização para depuração mais fácil.

### 4.5. Exercício Prático de Depuração

Agora é hora de colocar a abordagem sistemática de depuração em prática. O fluxo de trabalho `buggy_workflow.nf` contém vários erros comuns que representam os tipos de problemas que você encontrará no desenvolvimento do mundo real.

!!! exercise

    Use a abordagem sistemática de depuração para identificar e corrigir todos os erros em `buggy_workflow.nf`. Este fluxo de trabalho tenta processar dados de amostra de um arquivo CSV, mas contém múltiplos bugs intencionais representando cenários comuns de depuração.

    Comece executando o fluxo de trabalho para ver o primeiro erro:

    ```bash
    nextflow run buggy_workflow.nf
    ```

    ??? failure "Saída do comando"

        ```console
        N E X T F L O W   ~  version 25.10.2

        Launching `buggy_workflow.nf` [wise_ramanujan] DSL2 - revision: d51a8e83fd

        ERROR ~ Range [11, 12) out of bounds for length 11

         -- Check '.nextflow.log' file for details
        ```

        Este erro críptico indica um problema de análise em torno da linha 11-12 no bloco `params{}`. O analisador v2 captura problemas estruturais cedo.

    Aplique o método de depuração em quatro fases que você aprendeu:

    **Fase 1: Resolução de Erros de Sintaxe**
    - Verifique sublinhados vermelhos no VSCode ou sua IDE
    - Execute `nextflow run workflow.nf -preview` para identificar problemas de sintaxe
    - Corrija todos os erros de sintaxe (chaves ausentes, vírgulas finais, etc.)
    - Certifique-se de que o fluxo de trabalho seja analisado com sucesso antes de prosseguir

    **Fase 2: Avaliação Rápida**
    - Leia as mensagens de erro de tempo de execução cuidadosamente
    - Identifique se os erros são de tempo de execução, lógica ou relacionados a recursos
    - Use o modo `-preview` para testar a lógica básica do fluxo de trabalho

    **Fase 3: Investigação Detalhada**
    - Examine diretórios de trabalho para tarefas que falharam
    - Adicione operadores `.view()` para inspecionar canais
    - Verifique arquivos de log em diretórios de trabalho
    - Use `-stub-run` para testar a lógica do fluxo de trabalho sem execução

    **Fase 4: Corrigir e Validar**
    - Faça correções direcionadas
    - Use `-resume` para testar correções de forma eficiente
    - Verifique a execução completa do fluxo de trabalho

    **Ferramentas de Depuração à Sua Disposição:**
    ```bash
    # Preview mode for syntax checking
    nextflow run buggy_workflow.nf -preview

    # Debug profile for detailed output
    nextflow run buggy_workflow.nf -profile debug

    # Stub running for logic testing
    nextflow run buggy_workflow.nf -stub-run

    # Resume after fixes
    nextflow run buggy_workflow.nf -resume
    ```

    ??? solution
        O `buggy_workflow.nf` contém 9 ou 10 erros distintos (dependendo de como você conta) cobrindo todas as principais categorias de depuração. Aqui está uma análise sistemática de cada erro e como corrigi-lo

        Vamos começar com esses erros de sintaxe:

        **Erro 1: Erro de Sintaxe - Vírgula Final**
        ```groovy linenums="21"
        output:
            path "${sample_id}_result.txt",  // ERROR: Trailing comma
        ```
        **Correção:** Remova a vírgula final
        ```groovy linenums="21"
        output:
            path "${sample_id}_result.txt"
        ```

        **Erro 2: Erro de Sintaxe - Chave de Fechamento Ausente**
        ```groovy linenums="24"
        script:
        """
        echo "Processing: ${sample}"
        cat ${input_file} > ${sample}_result.txt
        """
        // ERROR: Missing closing brace for processFiles process
        ```
        **Correção:** Adicione a chave de fechamento ausente
        ```groovy linenums="29"
        """
        echo "Processing: ${sample_id}"
        cat ${input_file} > ${sample_id}_result.txt
        """
        }  // Add missing closing brace
        ```

        **Erro 3: Erro de Nome de Variável**
        ```groovy linenums="26"
        echo "Processing: ${sample}"     // ERROR: should be sample_id
        cat ${input_file} > ${sample}_result.txt  // ERROR: should be sample_id
        ```
        **Correção:** Use o nome de variável de entrada correto
        ```groovy linenums="26"
        echo "Processing: ${sample_id}"
        cat ${input_file} > ${sample_id}_result.txt
        ```

        **Erro 4: Erro de Variável Indefinida**
        ```groovy linenums="87"
        heavy_ch = heavyProcess(sample_ids)  // ERROR: sample_ids undefined
        ```
        **Correção:** Use o canal correto e extraia IDs de amostra
        ```groovy linenums="87"
        heavy_ch = heavyProcess(input_ch)
        ```

        Neste ponto o fluxo de trabalho será executado, mas ainda estaremos recebendo erros (por exemplo, `Path value cannot be null` em `processFiles`), causados por estrutura de canal ruim.

        **Erro 5: Erro de Estrutura de Canal - Saída de Map Errada**
        ```groovy linenums="83"
        .map { row -> row.sample_id }  // ERROR: processFiles expects tuple
        ```
        **Correção:** Retorne a estrutura de tupla que processFiles espera
        ```groovy linenums="83"
        .map { row -> [row.sample_id, file(row.fastq_path)] }
        ```

        Mas isso quebrará nossa correção para executar `heavyProcess()` acima, então precisaremos usar um map para passar apenas os IDs de amostra para esse processo:

        **Erro 6: Estrutura de canal ruim para heavyProcess**
        ```groovy linenums="87"
        heavy_ch = heavyProcess(input_ch)  // ERROR: input_ch now has 2 elements per emission- heavyProcess only needs 1 (the first)
        ```
        **Correção:** Use o canal correto e extraia IDs de amostra
        ```groovy linenums="87"
        heavy_ch = heavyProcess(input_ch.map{it[0]})
        ```

        Agora chegamos um pouco mais longe mas recebemos um erro sobre `No such variable: i`, porque não escapamos uma variável Bash.

        **Erro 7: Erro de Escape de Variável Bash**
        ```groovy linenums="48"
        echo "Heavy computation $i for ${sample_id}"  // ERROR: $i not escaped
        ```
        **Correção:** Escape a variável bash
        ```groovy linenums="48"
        echo "Heavy computation \${i} for ${sample_id}"
        ```

        Agora recebemos `Process exceeded running time limit (1ms)`, então corrigimos o limite de tempo de execução para o processo relevante:

        **Erro 8: Erro de Configuração de Recursos**
        ```groovy linenums="36"
        time '1 ms'  // ERROR: Unrealistic time limit
        ```
        **Correção:** Aumente para um limite de tempo realista
        ```groovy linenums="36"
        time '100 s'
        ```

        Em seguida, temos um erro `Missing output file(s)` para resolver:

        **Erro 9: Incompatibilidade de Nome de Arquivo de Saída**
        ```groovy linenums="49"
        done > ${sample_id}.txt  // ERROR: Wrong filename, should match output declaration
        ```
        **Correção:** Corresponda à declaração de saída
        ```groovy linenums="49"
        done > ${sample_id}_heavy.txt
        ```

        Os dois primeiros processos foram executados, mas não o terceiro.

        **Erro 10: Incompatibilidade de Nome de Arquivo de Saída**
        ```groovy linenums="88"
        file_ch = channel.fromPath("*.txt") // Error: attempting to take input from the pwd rather than a process
        handleFiles(file_ch)
        ```
        **Correção:** Pegue a saída do processo anterior
        ```groovy linenums="88"
        handleFiles(heavyProcess.out)
        ```

        Com isso, todo o fluxo de trabalho deve ser executado.

        **Fluxo de Trabalho Corrigido Completo:**
        ```groovy linenums="1"
        #!/usr/bin/env nextflow

        /*
        * Buggy workflow for debugging exercises
        * This workflow contains several intentional bugs for learning purposes
        */

        params{
            // Parameters with missing validation
            input: Path = 'data/sample_data.csv'
            output: String = 'results'
        }

        /*
        * Process with input/output mismatch
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
        * Process with resource issues
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
            # Simulate heavy computation
            for i in {1..1000000}; do
                echo "Heavy computation \$i for ${sample_id}"
            done > ${sample_id}_heavy.txt
            """
        }

        /*
        * Process with file handling issues
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
        * Main workflow with channel issues
        */
        workflow {

            // Channel with incorrect usage
            input_ch = channel
                .fromPath(params.input)
                .splitCsv(header: true)
                .map { row -> [row.sample_id, file(row.fastq_path)] }

            processed_ch = processFiles(input_ch)

            heavy_ch = heavyProcess(input_ch.map{it[0]})

            handleFiles(heavyProcess.out)
        }
        ```

**Categorias de Erros Cobertas:**

- **Erros de sintaxe**: Chaves ausentes, vírgulas finais, variáveis indefinidas
- **Erros de estrutura de canais**: Formas de dados erradas, canais indefinidos
- **Erros de processos**: Incompatibilidades de arquivos de saída, escape de variáveis
- **Erros de recursos**: Limites de tempo irrealistas

**Lições-Chave de Depuração:**

1. **Leia as mensagens de erro cuidadosamente** - elas frequentemente apontam diretamente para o problema
2. **Use abordagens sistemáticas** - corrija um erro de cada vez e teste com `-resume`
3. **Entenda o fluxo de dados** - erros de estrutura de canais são frequentemente os mais sutis
4. **Verifique diretórios de trabalho** - quando os processos falham, os logs dizem exatamente o que deu errado

---

## Resumo

Nesta missão secundária, você aprendeu um conjunto de técnicas sistemáticas para depurar fluxos de trabalho Nextflow.
Aplicar essas técnicas em seu próprio trabalho permitirá que você gaste menos tempo lutando com seu computador, resolva problemas mais rapidamente e se proteja de problemas futuros.

### Padrões-chave

**1. Como identificar e corrigir erros de sintaxe**:

- Interpretando mensagens de erro do Nextflow e localizando problemas
- Erros de sintaxe comuns: chaves ausentes, palavras-chave incorretas, variáveis indefinidas
- Distinguindo entre variáveis Nextflow (Groovy) e Bash
- Usando recursos de extensão do VS Code para detecção precoce de erros

```groovy
// Missing brace - look for red underlines in IDE
process FOO {
    script:
    """
    echo "hello"
    """
// } <-- missing!

// Wrong keyword
inputs:  // Should be 'input:'

// Undefined variable - escape with backslash for Bash variables
echo "${undefined_var}"      // Nextflow variable (error if not defined)
echo "\${bash_var}"          // Bash variable (escaped)
```

**2. Como depurar problemas de estrutura de canais**:

- Entendendo problemas de cardinalidade e esgotamento de canais
- Depurando incompatibilidades de estrutura de conteúdo de canais
- Usando operadores `.view()` para inspeção de canais
- Reconhecendo padrões de mensagens de erro como colchetes indicando estruturas de tuplas inesperadas

```groovy
// Inspect channel content
my_channel.view { "Content: $it" }

// Convert queue to value channel (prevents exhaustion)
reference_ch = channel.value('ref.fa')
// or
reference_ch = channel.of('ref.fa').first()
```

**3. Como solucionar problemas de execução de processos**:

- Diagnosticando erros de arquivos de saída ausentes
- Entendendo códigos de saída (127 para software ausente, 137 para problemas de memória)
- Investigando diretórios de trabalho e arquivos de comando
- Configurando recursos apropriadamente

```bash
# Check what was actually executed
cat work/ab/cdef12/.command.sh

# Check error output
cat work/ab/cdef12/.command.err

# Exit code 127 = command not found
# Exit code 137 = killed (memory/time limit)
```

**4. Como usar as ferramentas de depuração integradas do Nextflow**:

- Aproveitando o modo de visualização e depuração em tempo real
- Implementando execução stub para teste de lógica
- Aplicando resume para ciclos de depuração eficientes
- Seguindo uma metodologia sistemática de depuração em quatro fases

!!! tip "Referência Rápida de Depuração"

    **Erros de sintaxe?** → Verifique avisos do VSCode, execute `nextflow run workflow.nf -preview`

    **Problemas de canais?** → Use `.view()` para inspecionar conteúdo: `my_channel.view()`

    **Falhas de processos?** → Verifique arquivos do diretório de trabalho:

    - `.command.sh` - o script executado
    - `.command.err` - mensagens de erro
    - `.exitcode` - status de saída (127 = comando não encontrado, 137 = morto)

    **Comportamento misterioso?** → Execute com `-stub-run` para testar a lógica do fluxo de trabalho

    **Fez correções?** → Use `-resume` para economizar tempo testando: `nextflow run workflow.nf -resume`

---

### Recursos adicionais

- [Guia de solução de problemas do Nextflow](https://www.nextflow.io/docs/latest/troubleshooting.html): Documentação oficial de solução de problemas
- [Entendendo canais Nextflow](https://www.nextflow.io/docs/latest/channel.html): Mergulho profundo em tipos e comportamento de canais
- [Referência de diretivas de processos](https://www.nextflow.io/docs/latest/process.html#directives): Todas as opções de configuração de processos disponíveis
- [nf-test](https://www.nf-test.com/): Framework de testes para pipelines Nextflow
- [Comunidade Nextflow no Slack](https://www.nextflow.io/slack-invite.html): Obtenha ajuda da comunidade

Para fluxos de trabalho de produção, considere:

- Configurar [Seqera Platform](https://seqera.io/platform/) para monitoramento e depuração em escala
- Usar [contêineres Wave](https://seqera.io/wave/) para ambientes de software reproduzíveis

**Lembre-se:** Depuração eficaz é uma habilidade que melhora com a prática. A metodologia sistemática e o kit de ferramentas abrangente que você adquiriu aqui servirão bem ao longo de sua jornada de desenvolvimento Nextflow.

---

## O que vem a seguir?

Retorne ao [menu de Missões Secundárias](./index.md) ou clique no botão no canto inferior direito da página para passar para o próximo tópico da lista.
