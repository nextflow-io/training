# Depurando Workflows

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Tradução assistida por IA - [saiba mais e sugira melhorias](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

Depuração é uma habilidade crítica que pode economizar horas de frustração e ajudá-lo a se tornar um desenvolvedor Nextflow mais eficaz. Ao longo da sua carreira, especialmente quando você está começando, você encontrará bugs ao construir e manter seus workflows. Aprender abordagens sistemáticas de depuração ajudará você a identificar e resolver problemas rapidamente.

### Objetivos de aprendizagem

Nesta missão lateral, exploraremos **técnicas sistemáticas de depuração** para workflows Nextflow:

- **Depuração de erros de sintaxe**: Usando recursos de IDE e mensagens de erro do Nextflow efetivamente
- **Depuração de canais**: Diagnosticando problemas de fluxo de dados e problemas de estrutura de canais
- **Depuração de processos**: Investigando falhas de execução e problemas de recursos
- **Ferramentas de depuração integradas**: Aproveitando o modo de visualização, execução stub e diretórios de trabalho do Nextflow
- **Abordagens sistemáticas**: Uma metodologia de quatro fases para depuração eficiente

Ao final, você terá uma metodologia robusta de depuração que transforma mensagens de erro frustrantes em roteiros claros para soluções.

### Pré-requisitos

Antes de assumir esta missão lateral, você deve:

- Ter completado o tutorial [Hello Nextflow](../hello_nextflow/README.md) ou curso equivalente para iniciantes.
- Estar confortável usando conceitos e mecanismos básicos do Nextflow (processos, canais, operadores)

**Opcional:** Recomendamos completar a missão lateral [Recursos de IDE para Desenvolvimento Nextflow](./ide_features.md) primeiro.
Ela cobre cobertura abrangente de recursos de IDE que suportam depuração (destaque de sintaxe, detecção de erros, etc.), que usaremos bastante aqui.

---

## 0. Começando

#### Abra o codespace de treinamento

Se você ainda não o fez, certifique-se de abrir o ambiente de treinamento conforme descrito na [Configuração do Ambiente](../envsetup/index.md).

[![Abrir no GitHub Codespaces](https://github.com/codespaces/badge.svg)](https://codespaces.new/nextflow-io/training?quickstart=1&ref=master)

#### Mova para o diretório do projeto

Vamos mover para o diretório onde os arquivos para este tutorial estão localizados.

```bash
cd side-quests/debugging
```

Você pode configurar o VSCode para focar neste diretório:

```bash
code .
```

#### Revise os materiais

Você encontrará um conjunto de workflows de exemplo com vários tipos de bugs que usaremos para prática:

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

Esses arquivos representam cenários comuns de depuração que você encontrará no desenvolvimento real.

#### Revise a tarefa

Seu desafio é executar cada workflow, identificar o(s) erro(s) e corrigi-los.

Para cada workflow com bugs:

1. **Execute o workflow** e observe o erro
2. **Analise a mensagem de erro**: o que o Nextflow está dizendo?
3. **Localize o problema** no código usando as pistas fornecidas
4. **Corrija o bug** e verifique se sua solução funciona
5. **Restaure o arquivo** antes de passar para a próxima seção (use `git checkout <filename>`)

Os exercícios progridem de erros de sintaxe simples para problemas de runtime mais sutis.
Soluções são discutidas inline, mas tente resolver cada um você mesmo antes de ler adiante.

#### Lista de verificação de prontidão

Acha que está pronto para mergulhar?

- [ ] Eu entendo o objetivo deste curso e seus pré-requisitos
- [ ] Meu codespace está funcionando
- [ ] Configurei meu diretório de trabalho adequadamente
- [ ] Eu entendo a tarefa

Se você pode marcar todas as caixas, está pronto para começar.

---

## 1. Erros de Sintaxe

Erros de sintaxe são o tipo mais comum de erro que você encontrará ao escrever código Nextflow. Eles ocorrem quando o código não está conforme as regras de sintaxe esperadas do Nextflow DSL. Esses erros impedem que seu workflow seja executado, então é importante aprender a identificá-los e corrigi-los rapidamente.

### 1.1. Chaves ausentes

Um dos erros de sintaxe mais comuns, e às vezes um dos mais complexos de depurar, é **chaves ausentes ou incompatíveis**.

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
- **Descrição do erro**: Explica o que o parser encontrou que não esperava (`Unexpected input: '<EOF>'`)
- **Indicador EOF**: A mensagem `<EOF>` (End Of File) indica que o parser alcançou o fim do arquivo enquanto ainda esperava mais conteúdo - um sinal clássico de chaves não fechadas

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
// Chave de fechamento ausente para o processo

workflow {

    // Cria canal de entrada
    input_ch = channel.of('sample1', 'sample2', 'sample3')

    // Chama o processo com o canal de entrada
    PROCESS_FILES(input_ch)
}
```

Para o propósito deste exemplo, deixamos um comentário para você mostrar onde está o erro. A extensão Nextflow do VSCode também deve estar dando algumas dicas sobre o que pode estar errado, colocando a chave incompatível em vermelho e destacando o fim prematuro do arquivo:

![Bad syntax](img/bad_syntax.png)

**Estratégia de depuração para erros de chaves:**

1. Use a correspondência de chaves do VS Code (coloque o cursor ao lado de uma chave)
2. Verifique o painel Problemas para mensagens relacionadas a chaves
3. Certifique-se de que cada `{` de abertura tem um `}` de fechamento correspondente

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
    }  // Adiciona a chave de fechamento ausente

    workflow {

        // Cria canal de entrada
        input_ch = channel.of('sample1', 'sample2', 'sample3')

        // Chama o processo com o canal de entrada
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
    // Chave de fechamento ausente para o processo

    workflow {

        // Cria canal de entrada
        input_ch = channel.of('sample1', 'sample2', 'sample3')

        // Chama o processo com o canal de entrada
        PROCESS_FILES(input_ch)
    }
    ```

#### Execute o pipeline

Agora execute o workflow novamente para confirmar que funciona:

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

    // Cria canal de entrada
    input_ch = channel.of('sample1', 'sample2', 'sample3')

    // Chama o processo com o canal de entrada
    PROCESS_FILES(input_ch)
}
```

Olhando para a linha 4 no contexto do erro, podemos identificar o problema: estamos usando `inputs` em vez da diretiva correta `input`. A extensão Nextflow do VSCode também sinalizará isso:

![Invalid process message](img/invalid_process_message.png)

#### Corrija o código

Substitua a palavra-chave incorreta pela correta referenciando [a documentação](https://www.nextflow.io/docs/latest/process.html#):

=== "Depois"

    ```groovy title="invalid_process.nf" hl_lines="4" linenums="1"
    #!/usr/bin/env nextflow

    process PROCESS_FILES {
        input:  // Corrigido: Mudou 'inputs' para 'input'
        val sample_name

        output:
        path "${sample_name}_output.txt"

        script:
        """
        echo "Processing ${sample_name}" > ${sample_name}_output.txt
        """
    }

    workflow {

        // Cria canal de entrada
        input_ch = channel.of('sample1', 'sample2', 'sample3')

        // Chama o processo com o canal de entrada
        PROCESS_FILES(input_ch)
    }
    ```

=== "Antes"

    ```groovy title="invalid_process.nf" hl_lines="4" linenums="1"
    #!/usr/bin/env nextflow

    process PROCESS_FILES {
        inputs:  // ERRO: Deveria ser 'input' não 'inputs'
        val sample_name

        output:
        path "${sample_name}_output.txt"

        script:
        """
        echo "Processing ${sample_name}" > ${sample_name}_output.txt
        """
    }

    workflow {

        // Cria canal de entrada
        input_ch = channel.of('sample1', 'sample2', 'sample3')

        // Chama o processo com o canal de entrada
        PROCESS_FILES(input_ch)
    }
    ```

#### Execute o pipeline

Agora execute o workflow novamente para confirmar que funciona:

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

Os nomes de variáveis que você usa em seus blocos de script devem ser válidos, derivados de entradas ou de código groovy inserido antes do script. Mas quando você está lidando com complexidade no início do desenvolvimento do pipeline, é fácil cometer erros na nomeação de variáveis, e o Nextflow vai avisar rapidamente.

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

O erro é capturado no momento da compilação e aponta diretamente para a variável indefinida na linha 17, com um cursor indicando exatamente onde está o problema.

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
    // Define variáveis no código Groovy antes do script
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

A mensagem de erro indica que a variável não é reconhecida no template do script, e lá você vê - `${undefined_var}` usado no bloco de script, mas não definido em outro lugar.

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
        // Define variáveis no código Groovy antes do script
        def output_prefix = "${sample_name}_processed"
        def timestamp = new Date().format("yyyy-MM-dd")

        """
        echo "Processing ${sample_name} on ${timestamp}" > ${output_prefix}.txt
        """  // Removida a linha com undefined_var
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
        // Define variáveis no código Groovy antes do script
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

Agora execute o workflow novamente para confirmar que funciona:

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

### 1.4. Mau uso de variáveis Bash

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

Agora execute o workflow novamente para confirmar que funciona:

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

    Para manipulações simples de variáveis como concatenação de strings ou operações de prefixo/sufixo, geralmente é mais legível usar variáveis Groovy na seção de script em vez de variáveis Bash no bloco de script:

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

A extensão Nextflow do VSCode destaca problemas com a estrutura do código que causarão erros. Um exemplo comum é definir canais fora do bloco `workflow {}` - isso agora é imposto como um erro de sintaxe.

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

A mensagem de erro indica claramente o problema: declarações (como definições de canal) não podem ser misturadas com declarações de script fora de um bloco workflow ou process.

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
    // Define variáveis no código Groovy antes do script
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

![Non-lethal syntax error](img/nonlethal.png)

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
        // Define variáveis no código Groovy antes do script
        def output_prefix = "${sample_name}_processed"
        def timestamp = new Date().format("yyyy-MM-dd")

        """
        echo "Processing ${sample_name} on ${timestamp}" > ${output_prefix}.txt
        """
    }

    workflow {
        input_ch = channel.of('sample1', 'sample2', 'sample3')  // Movido para dentro do bloco workflow
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
        // Define variáveis no código Groovy antes do script
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

Execute o workflow novamente para confirmar que a correção funciona:

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

Mantenha seus canais de entrada definidos dentro do bloco workflow e, em geral, siga quaisquer outras recomendações que a extensão faça.

### Conclusão

Você pode identificar e corrigir erros de sintaxe sistematicamente usando mensagens de erro do Nextflow e indicadores visuais da IDE. Erros de sintaxe comuns incluem chaves ausentes, palavras-chave de processo incorretas, variáveis indefinidas e uso inadequado de variáveis Bash vs. Nextflow. A extensão VSCode ajuda a capturar muitos desses erros antes do runtime. Com essas habilidades de depuração de sintaxe em seu toolkit, você será capaz de resolver rapidamente os erros de sintaxe mais comuns do Nextflow e passar para lidar com problemas de runtime mais complexos.

### O que vem a seguir?

Aprenda a depurar erros de estrutura de canal mais complexos que ocorrem mesmo quando a sintaxe está correta.

---

## 2. Erros de Estrutura de Canal

Erros de estrutura de canal são mais sutis do que erros de sintaxe porque o código está sintaticamente correto, mas as formas dos dados não correspondem ao que os processos esperam. O Nextflow tentará executar o pipeline, mas pode descobrir que o número de entradas não corresponde ao que espera e falhar. Esses erros geralmente aparecem apenas no runtime e requerem compreensão dos dados fluindo através do seu workflow.

!!! tip "Depurando Canais com `.view()`"

    Ao longo desta seção, lembre-se de que você pode usar o operador `.view()` para inspecionar o conteúdo do canal em qualquer ponto do seu workflow. Esta é uma das ferramentas de depuração mais poderosas para entender problemas de estrutura de canal. Exploraremos essa técnica em detalhes na seção 2.4, mas sinta-se livre para usá-la enquanto trabalha nos exemplos.

    ```groovy
    my_channel.view()  // Mostra o que está fluindo através do canal
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
        val sample_name  // Processo espera apenas 1 entrada

    output:
        path "${sample_name}_output.txt"

    script:
    """
    echo "Processing ${sample_name}" > ${sample_name}_output.txt
    """
}

workflow {

    // Cria dois canais separados
    samples_ch = channel.of('sample1', 'sample2', 'sample3')
    files_ch = channel.of('file1.txt', 'file2.txt', 'file3.txt')

    // ERROR: Passing 2 channels but process expects only 1
    PROCESS_FILES(samples_ch, files_ch)
}
```

Você deve ver a chamada `PROCESS_FILES` incompatível, fornecendo múltiplos canais de entrada quando o processo define apenas um. A extensão VSCode também sublinhará a chamada do processo em vermelho e fornecerá uma mensagem de diagnóstico quando você passar o mouse:

![Incorrect number of args message](img/incorrect_num_args.png)

#### Corrija o código

Para este exemplo específico, o processo espera um único canal e não requer o segundo canal, então podemos corrigi-lo passando apenas o canal `samples_ch`:

=== "Depois"

    ```groovy title="bad_number_inputs.nf" hl_lines="23" linenums="1"
    #!/usr/bin/env nextflow

    process PROCESS_FILES {
        input:
            val sample_name  // Processo espera apenas 1 entrada

        output:
            path "${sample_name}_output.txt"

        script:
        """
        echo "Processing ${sample_name}" > ${sample_name}_output.txt
        """
    }

    workflow {

        // Cria dois canais separados
        samples_ch = channel.of('sample1', 'sample2', 'sample3')
        files_ch = channel.of('file1.txt', 'file2.txt', 'file3.txt')

        // Corrigido: Passa apenas o canal que o processo espera
        PROCESS_FILES(samples_ch)
    }
    ```

=== "Antes"

    ```groovy title="bad_number_inputs.nf" hl_lines="5 23" linenums="1"
    #!/usr/bin/env nextflow

    process PROCESS_FILES {
        input:
            val sample_name  // Processo espera apenas 1 entrada

        output:
            path "${sample_name}_output.txt"

        script:
        """
        echo "Processing ${sample_name}" > ${sample_name}_output.txt
        """
    }

    workflow {

        // Cria dois canais separados
        samples_ch = channel.of('sample1', 'sample2', 'sample3')
        files_ch = channel.of('file1.txt', 'file2.txt', 'file3.txt')

        // ERRO: Passando 2 canais mas o processo espera apenas 1
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

Mais comumente do que neste exemplo, você pode adicionar entradas adicionais a um processo e esquecer de atualizar a chamada do workflow de acordo, o que pode levar a este tipo de erro. Felizmente, este é um dos erros mais fáceis de entender e corrigir, já que a mensagem de erro é bastante clara sobre a incompatibilidade.

### 2.2. Exaustão de Canal (Processo Executa Menos Vezes do que o Esperado)

Alguns erros de estrutura de canal são muito mais sutis e não produzem erros. Provavelmente o mais comum destes reflete um desafio que novos usuários Nextflow enfrentam ao entender que canais de fila podem ser esgotados e ficar sem itens, significando que o workflow termina prematuramente.

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

Este workflow completa sem erro, mas processa apenas uma única amostra!

#### Verifique o código

Vamos examinar `exhausted.nf` para ver se está correto:

```groovy title="exhausted.nf" hl_lines="23 24" linenums="1"
#!/usr/bin/env nextflow

process PROCESS_FILES {
    input:
    val reference
    val sample_name

    output:
    path "${output_prefix}.txt"

    script:
    // Define variáveis no código Groovy antes do script
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

O processo executa apenas uma vez em vez de três vezes porque o canal `reference_ch` é um canal de fila que é esgotado após a primeira execução do processo. Quando um canal é esgotado, o processo inteiro para, mesmo que outros canais ainda tenham itens.

Este é um padrão comum onde você tem um único arquivo de referência que precisa ser reutilizado em múltiplas amostras. A solução é converter o canal de referência em um canal de valor que pode ser reutilizado indefinidamente.

#### Corrija o código

Existem algumas maneiras de resolver isso dependendo de quantos arquivos são afetados.

**Opção 1**: Você tem um único arquivo de referência que está reutilizando muito. Você pode simplesmente criar um tipo de canal de valor, que pode ser usado repetidamente. Existem três maneiras de fazer isso:

**1a** Use `channel.value()`:

```groovy title="exhausted.nf (corrigido - Opção 1a)" hl_lines="2" linenums="21"
workflow {
    reference_ch = channel.value('baseline_reference')  // Canal de valor pode ser reutilizado
    input_ch = channel.of('sample1', 'sample2', 'sample3')

    PROCESS_FILES(reference_ch, input_ch)
}
```

**1b** Use o operador `first()` [operator](https://www.nextflow.io/docs/latest/reference/operator.html#first):

```groovy title="exhausted.nf (corrigido - Opção 1b)" hl_lines="2" linenums="21"
workflow {
    reference_ch = channel.of('baseline_reference').first()  // Converte para canal de valor
    input_ch = channel.of('sample1', 'sample2', 'sample3')

    PROCESS_FILES(reference_ch, input_ch)
}
```

**1c.** Use o operador `collect()` [operator](https://www.nextflow.io/docs/latest/reference/operator.html#collect):

```groovy title="exhausted.nf (corrigido - Opção 1c)" hl_lines="2" linenums="21"
workflow {
    reference_ch = channel.of('baseline_reference').collect()  // Converte para canal de valor
    input_ch = channel.of('sample1', 'sample2', 'sample3')

    PROCESS_FILES(reference_ch, input_ch)
}
```

**Opção 2**: Em cenários mais complexos, talvez onde você tenha múltiplos arquivos de referência para todas as amostras no canal de amostras, você pode usar o operador `combine` para criar um novo canal que combine os dois canais em tuplas:

```groovy title="exhausted.nf (corrigido - Opção 2)" hl_lines="4" linenums="21"
workflow {
    reference_ch = channel.of('baseline_reference','other_reference')
    input_ch = channel.of('sample1', 'sample2', 'sample3')
    combined_ch = reference_ch.combine(input_ch)  // Cria produto cartesiano

    PROCESS_FILES(combined_ch)
}
```

O operador `.combine()` gera um produto cartesiano dos dois canais, então cada item em `reference_ch` será pareado com cada item em `input_ch`. Isso permite que o processo execute para cada amostra enquanto ainda usa a referência.

Isso requer que a entrada do processo seja ajustada. No nosso exemplo, o início da definição do processo precisaria ser ajustado da seguinte forma:

```groovy title="exhausted.nf (corrigido - Opção 2)" hl_lines="5" linenums="1"
#!/usr/bin/env nextflow

process PROCESS_FILES {
    input:
        tuple val(reference), val(sample_name)
```

Esta abordagem pode não ser adequada em todas as situações.

#### Execute o pipeline

Tente uma das correções acima e execute o workflow novamente:

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

Quando os workflows atingem um certo nível de complexidade, pode ser um pouco difícil acompanhar as estruturas internas de cada canal, e as pessoas comumente geram incompatibilidades entre o que o processo espera e o que o canal realmente contém. Isso é mais sutil do que o problema que discutimos anteriormente, onde o número de canais estava incorreto. Neste caso, você pode ter o número correto de canais de entrada, mas a estrutura interna de um ou mais desses canais não corresponde ao que o processo espera.

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
        val sample_name  // Espera valor único, recebe tupla

    output:
        path "${sample_name}_output.txt"

    script:
    """
    echo "Processing ${sample_name}" > ${sample_name}_output.txt
    """
}

workflow {

    // Canal emite tuplas, mas o processo espera valores únicos
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
                tuple val(sample_name), val(file_name)  // Corrigido: Aceita tupla

            output:
                path "${sample_name}_output.txt"

            script:
            """
            echo "Processing ${sample_name}" > ${sample_name}_output.txt
            """
        }

        workflow {

        // Canal emite tuplas, mas o processo espera valores únicos
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
                val sample_name  // Espera valor único, recebe tupla

            output:
                path "${sample_name}_output.txt"

            script:
            """
            echo "Processing ${sample_name}" > ${sample_name}_output.txt
            """
        }

        workflow {

        // Canal emite tuplas, mas o processo espera valores únicos
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

        // Canal emite tuplas, mas o processo espera valores únicos
            input_ch = channel.of(
              ['sample1', 'file1.txt'],
              ['sample2', 'file2.txt'],
              ['sample3', 'file3.txt']
            )
            PROCESS_FILES(input_ch.map { it[0] })  // Corrigido: Extrai primeiro elemento
        }
        ```

    === "Antes"

        ```groovy title="bad_channel_shape.nf" hl_lines="9" linenums="16"
        workflow {

        // Canal emite tuplas, mas o processo espera valores únicos
            input_ch = channel.of(
              ['sample1', 'file1.txt'],
              ['sample2', 'file2.txt'],
              ['sample3', 'file3.txt']
            )
            PROCESS_FILES(input_ch)
        }
        ```

#### Execute o pipeline

Escolha uma das soluções e execute novamente o workflow:

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

### 2.4. Técnicas de Depuração de Canal

#### Usando `.view()` para Inspeção de Canal

A ferramenta de depuração mais poderosa para canais é o operador `.view()`. Com `.view()`, você pode entender a forma dos seus canais em todos os estágios para ajudar na depuração.

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

    // Canal emite tuplas, mas o processo espera valores únicos
    input_ch = channel.of(
      ['sample1', 'file1.txt'],
      ['sample2', 'file2.txt'],
      ['sample3', 'file3.txt']
    )
    .view { "Channel content: $it" }  // Debug: Mostra conteúdo original do canal
    .map { tuple -> tuple[0] }        // Transforma: Extrai primeiro elemento
    .view { "After mapping: $it" }    // Debug: Mostra conteúdo do canal transformado

    PROCESS_FILES(input_ch)
}
```

#### Corrija o código

Para evitar usar operações `.view()` excessivamente no futuro para entender o conteúdo do canal, é aconselhável adicionar alguns comentários para ajudar:

```groovy title="bad_channel_shape_viewed.nf (com comentários)" linenums="16" hl_lines="8 9"
workflow {

    // Canal emite tuplas, mas o processo espera valores únicos
    input_ch = channel.of(
            ['sample1', 'file1.txt'],
            ['sample2', 'file2.txt'],
            ['sample3', 'file3.txt'],
        ) // [sample_name, file_name]
        .map { tuple -> tuple[0] } // sample_name

    PROCESS_FILES(input_ch)
}
```

Isso se tornará mais importante à medida que seus workflows crescem em complexidade e a estrutura do canal se torna mais opaca.

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

Muitos erros de estrutura de canal podem ser criados com sintaxe Nextflow válida. Você pode depurar erros de estrutura de canal compreendendo o fluxo de dados, usando operadores `.view()` para inspeção e reconhecendo padrões de mensagem de erro como colchetes indicando estruturas de tupla inesperadas.

### O que vem a seguir?

Aprenda sobre erros criados por definições de processo.

---

## 3. Erros de Estrutura de Processo

A maioria dos erros que você encontrará relacionados a processos estará relacionada a erros que você cometeu na formação do comando, ou a problemas relacionados ao software subjacente. Dito isso, similarmente aos problemas de canal acima, você pode cometer erros na definição do processo que não se qualificam como erros de sintaxe, mas que causarão erros no runtime.

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
    path "${sample_name}.txt"  // Espera: sample3.txt

    script:
    """
    echo "Processing ${sample_name}" > ${sample_name}_output.txt  // Cria: sample3_output.txt
    """
}
```

Você deve ver que há uma incompatibilidade entre o nome do arquivo de saída no bloco `output:` e o usado no script. Esta incompatibilidade faz com que o processo falhe. Se você encontrar esse tipo de erro, volte e verifique se as saídas correspondem entre sua definição de processo e seu bloco de saída.

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
        path "${sample_name}_output.txt"  // Corrigido: Corresponde à saída do script

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
        path "${sample_name}.txt"  // Espera: sample3.txt

        script:
        """
        echo "Processing ${sample_name}" > ${sample_name}_output.txt  // Cria: sample3_output.txt
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

Outra classe de erros ocorre devido a erros no provisionamento de software. `missing_software.nf` é um workflow sintaticamente válido, mas depende de algum software externo para fornecer o comando `cowpy` que ele usa.

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

O processo não tem acesso ao comando que estamos especificando. Às vezes isso ocorre porque um script está presente no diretório `bin` do workflow, mas não foi tornado executável. Outras vezes é porque o software não está instalado no contêiner ou ambiente onde o workflow está sendo executado.

#### Verifique o código

Fique atento ao código de saída `127` - ele diz exatamente qual é o problema. Vamos examinar `missing_software.nf`:

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

Fomos um pouco desonestos aqui, e na verdade não há nada de errado com o código. Só precisamos especificar a configuração necessária para executar o processo de forma que ele tenha acesso ao comando em questão. Neste caso, o processo tem uma definição de contêiner, então tudo o que precisamos fazer é executar o workflow com Docker habilitado.

#### Execute o pipeline

Configuramos um perfil Docker para você em `nextflow.config`, então você pode executar o workflow com:

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

    Para aprender mais sobre como o Nextflow usa contêineres, veja [Hello Nextflow](../hello_nextflow/05_hello_containers.md)

### 3.3. Má configuração de recursos

No uso em produção, você estará configurando recursos em seus processos. Por exemplo, `memory` define a quantidade máxima de memória disponível para seu processo, e se o processo exceder isso, seu agendador normalmente matará o processo e retornará um código de saída de `137`. Não podemos demonstrar isso aqui porque estamos usando o executor `local`, mas podemos mostrar algo similar com `time`.

#### Execute o pipeline

`bad_resources.nf` tem configuração de processo com um limite não realista de tempo de 1 milissegundo:

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

```groovy

```
