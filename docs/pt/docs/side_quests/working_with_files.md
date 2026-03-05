# Processamento de entrada de arquivos

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Tradução assistida por IA - [saiba mais e sugira melhorias](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

Fluxos de trabalho de análise científica frequentemente envolvem o processamento de grandes números de arquivos.
Nextflow fornece ferramentas poderosas para lidar com arquivos de forma eficiente, ajudando você a organizar e processar seus dados com o mínimo de código.

### Objetivos de aprendizado

Nesta side quest, vamos explorar como Nextflow lida com arquivos, desde operações básicas até técnicas mais avançadas para trabalhar com coleções de arquivos.
Você aprenderá como extrair metadados de nomes de arquivos, que é um requisito comum em pipelines de análise científica.

Ao final desta side quest, você será capaz de:

- Criar objetos Path a partir de strings de caminhos de arquivos usando o método `file()` do Nextflow
- Acessar atributos de arquivo como nome, extensão e diretório pai
- Lidar com arquivos locais e remotos de forma transparente usando URIs
- Usar canais para automatizar o manuseio de arquivos com `channel.fromPath()` e `channel.fromFilePairs()`
- Extrair e estruturar metadados de nomes de arquivos usando manipulação de strings
- Agrupar arquivos relacionados usando correspondência de padrões e expressões glob
- Integrar operações de arquivo em processos Nextflow com manuseio de entrada apropriado
- Organizar saídas de processos usando estruturas de diretórios baseadas em metadados

Essas habilidades ajudarão você a construir fluxos de trabalho que podem lidar com diferentes tipos de entradas de arquivo com grande flexibilidade.

### Pré-requisitos

Antes de enfrentar esta side quest, você deve:

- Ter completado o tutorial [Hello Nextflow](../../hello_nextflow/) ou curso equivalente para iniciantes.
- Estar confortável usando conceitos e mecanismos básicos do Nextflow (processos, canais, operadores)

<!-- I removed the suggestion to do the metamaps SQ first because that works more naturally after -->

---

## 0. Comece

#### Abra o codespace de treinamento

Se você ainda não fez isso, certifique-se de abrir o ambiente de treinamento conforme descrito em [Configuração do Ambiente](../envsetup/index.md).

[![Open in GitHub Codespaces](https://github.com/codespaces/badge.svg)](https://codespaces.new/nextflow-io/training?quickstart=1&ref=master)

#### Entre no diretório do projeto

Vamos entrar no diretório onde os arquivos para este tutorial estão localizados.

```bash
cd side-quests/working_with_files
```

Você pode configurar o VSCode para focar neste diretório:

```bash
code .
```

#### Revise os materiais

Você encontrará um arquivo de fluxo de trabalho simples chamado `main.nf`, um diretório `modules` contendo dois arquivos de módulo, e um diretório `data` contendo alguns arquivos de dados de exemplo.

??? abstract "Conteúdo do diretório"

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

Este diretório contém dados de sequenciamento paired-end de três pacientes (A, B, C).

Para cada paciente, temos amostras que são do tipo `tumor` (tipicamente originando de biópsias de tumor) ou `normal` (retiradas de tecido saudável ou sangue).
Se você não está familiarizado com análise de câncer, saiba apenas que isso corresponde a um modelo experimental que usa amostras pareadas tumor/normal para realizar análises contrastivas.

Para o paciente A especificamente, temos dois conjuntos de replicatas técnicas (repetições).

Os arquivos de dados de sequenciamento são nomeados com uma convenção típica `_R1_` e `_R2_` para o que são conhecidas como 'reads forward' e 'reads reverse'.

_Não se preocupe se você não está familiarizado com este desenho experimental, não é crítico para entender este tutorial._

#### Revise a tarefa

Seu desafio é escrever um fluxo de trabalho Nextflow que vai:

1. **Carregar** arquivos de entrada usando os métodos de manuseio de arquivos do Nextflow
2. **Extrair** metadados (ID do paciente, replicata, tipo de amostra) da estrutura do nome do arquivo
3. **Agrupar** arquivos pareados (R1/R2) juntos usando `channel.fromFilePairs()`
4. **Processar** os arquivos com um módulo de análise fornecido
5. **Organizar** saídas em uma estrutura de diretórios baseada nos metadados extraídos

#### Lista de verificação de prontidão

Acha que está pronto para mergulhar?

- [ ] Eu entendo o objetivo deste curso e seus pré-requisitos
- [ ] Meu codespace está funcionando
- [ ] Eu configurei meu diretório de trabalho apropriadamente
- [ ] Eu entendo a tarefa

Se você pode marcar todas as caixas, está pronto para começar.

---

## 1. Operações básicas de arquivo

### 1.1. Identifique o tipo de um objeto com `.class`

Dê uma olhada no arquivo de fluxo de trabalho `main.nf`:

```groovy title="main.nf" linenums="1"
#!/usr/bin/env nextflow

workflow {

    // Cria um objeto Path a partir de uma string de caminho
    myFile = 'data/patientA_rep1_normal_R1_001.fastq.gz'

    println "${myFile} is of class ${myFile.class}"
}
```

Este é um mini-fluxo de trabalho (sem nenhum processo) que se refere a um único caminho de arquivo em seu fluxo de trabalho, depois o imprime no console, junto com sua classe.

??? info "O que é `.class`?"

    No Nextflow, `.class` nos diz com que tipo de objeto estamos trabalhando. É como perguntar "que tipo de coisa é isso?" para descobrir se é uma string, um número, um arquivo ou outra coisa.
    Isso vai nos ajudar a ilustrar a diferença entre uma string simples e um objeto Path nas próximas seções.

Vamos executar o fluxo de trabalho:

```bash
nextflow run main.nf
```

??? success "Saída do comando"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `main.nf` [romantic_chandrasekhar] DSL2 - revision: 5a4a89bc3a

    data/patientA_rep1_normal_R1_001.fastq.gz is of class java.lang.String
    ```

Como você pode ver, Nextflow imprimiu a string de caminho exatamente como nós escrevemos.

Isso é apenas saída de texto; Nextflow ainda não fez nada especial com isso.
Também confirmamos que para o Nextflow, isso é apenas uma string (da classe `java.lang.String`).
Isso faz sentido, já que ainda não dissemos ao Nextflow que corresponde a um arquivo.

### 1.2. Crie um objeto Path com file()

Podemos dizer ao Nextflow como lidar com arquivos criando [objetos Path](https://www.nextflow.io/docs/latest/reference/stdlib-types.html#path) a partir de strings de caminho.

Em nosso fluxo de trabalho, podemos converter a string de caminho `data/patientA_rep1_normal_R1_001.fastq.gz` em um objeto Path usando o método `file()`, que fornece acesso a propriedades e operações de arquivo.

Edite o `main.nf` para envolver a string com `file()` da seguinte forma:

=== "Depois"

    ```groovy title="main.nf" linenums="5" hl_lines="2"
        // Cria um objeto Path a partir de uma string de caminho
        myFile = file('data/patientA_rep1_normal_R1_001.fastq.gz')

        println "${myFile} is of class ${myFile.class}"
    ```

=== "Antes"

    ```groovy title="main.nf" linenums="5" hl_lines="2"
        // Cria um objeto Path a partir de uma string de caminho
        myFile = 'data/patientA_rep1_normal_R1_001.fastq.gz'

        println "${myFile} is of class ${myFile.class}"
    ```

Agora execute o fluxo de trabalho novamente:

```bash
nextflow run main.nf
```

??? success "Saída do comando"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `main.nf` [kickass_coulomb] DSL2 - revision: 5af44b1b59

    /workspaces/training/side-quests/working_with_files/data/patientA_rep1_normal_R1_001.fastq.gz is of class class sun.nio.fs.UnixPath
    ```

Desta vez, você vê o caminho absoluto completo em vez do caminho relativo que fornecemos como entrada.

Nextflow converteu nossa string em um objeto Path e resolveu para a localização real do arquivo no sistema.
O caminho do arquivo agora será absoluto, como em `/workspaces/training/side-quests/working_with_files/data/patientA_rep1_normal_R1_001.fastq.gz`.

Note também que a classe do objeto Path é `sun.nio.fs.UnixPath`: esta é a forma do Nextflow de representar arquivos locais.
Como veremos mais tarde, arquivos remotos terão nomes de classe diferentes (como `nextflow.file.http.XPath` para arquivos HTTP), mas todos funcionam exatamente da mesma forma e podem ser usados identicamente em seus fluxos de trabalho.

!!! tip

    **A diferença chave:**

    - **String de caminho**: Apenas texto que Nextflow trata como caracteres
    - **Objeto Path**: Uma referência de arquivo inteligente com a qual Nextflow pode trabalhar

    Pense nisso assim: uma string de caminho é como escrever um endereço no papel, enquanto um objeto Path é como ter o endereço carregado em um dispositivo GPS que sabe como navegar até lá e pode te dizer detalhes sobre a jornada.

### 1.3. Acesse atributos de arquivo

Por que isso é útil? Bem, agora que Nextflow entende que `myFile` é um objeto Path e não apenas uma string, podemos acessar os vários atributos do objeto Path.

Vamos atualizar nosso fluxo de trabalho para imprimir os atributos de arquivo embutidos:

=== "Depois"

    ```groovy title="main.nf" linenums="5" hl_lines="4-9"
        // Cria um objeto Path a partir de uma string de caminho
        myFile = file('data/patientA_rep1_normal_R1_001.fastq.gz')

        // Imprime atributos do arquivo
        println "File object class: ${myFile.class}"
        println "File name: ${myFile.name}"
        println "Simple name: ${myFile.simpleName}"
        println "Extension: ${myFile.extension}"
        println "Parent directory: ${myFile.parent}"
    ```

=== "Antes"

    ```groovy title="main.nf" linenums="5" hl_lines="4"
        // Cria um objeto Path a partir de uma string de caminho
        myFile = file('data/patientA_rep1_normal_R1_001.fastq.gz')

        println "${myFile} is of class ${myFile.class}"
    ```

Execute o fluxo de trabalho:

```bash
nextflow run main.nf
```

??? success "Saída do comando"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `main.nf` [ecstatic_ampere] DSL2 - revision: f3fa3dcb48

    File object class: sun.nio.fs.UnixPath
    File name: patientA_rep1_normal_R1_001.fastq.gz
    Simple name: patientA_rep1_normal_R1_001
    Extension: gz
    Parent directory: /workspaces/training/side-quests/working_with_files/data
    ```

Você vê os vários atributos de arquivo impressos no console acima.

### 1.4. Forneça o arquivo para um processo

A diferença entre strings e objetos Path se torna crítica quando você começa a construir fluxos de trabalho reais com processos.
Até agora verificamos que Nextflow está agora tratando nosso arquivo de entrada como um arquivo, mas vamos ver se realmente podemos executar algo nesse arquivo em um processo.

#### 1.4.1. Importe o processo e examine o código

Fornecemos um módulo de processo pré-escrito chamado `COUNT_LINES` que recebe uma entrada de arquivo e conta quantas linhas há nele.

Para usar o processo no fluxo de trabalho, você só precisa adicionar uma declaração include antes do bloco workflow:

=== "Depois"

    ```groovy title="main.nf" linenums="1" hl_lines="3"
    #!/usr/bin/env nextflow

    include { COUNT_LINES } from './modules/count_lines.nf'

    workflow {
    ```

=== "Antes"

    ```groovy title="main.nf" linenums="1"
    #!/usr/bin/env nextflow

    workflow {
    ```

Você pode abrir o arquivo de módulo para examinar seu código:

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

Como você pode ver, é um script bem direto que descompacta o arquivo e conta quantas linhas ele contém.

??? info "O que faz `debug true`?"

    A diretiva `debug true` na definição do processo faz com que Nextflow imprima a saída do seu script (como a contagem de linhas "40") diretamente no log de execução.
    Sem isso, você veria apenas o status de execução do processo, mas não a saída real do seu script.

    Para mais informações sobre depuração de processos Nextflow, veja a side quest [Debugging Nextflow Workflows](debugging.md).

#### 1.4.2. Adicione uma chamada para `COUNT_LINES`

Agora que o processo está disponível para o fluxo de trabalho, podemos adicionar uma chamada ao processo `COUNT_LINES` para executá-lo no arquivo de entrada.

Faça as seguintes edições no fluxo de trabalho:

=== "Depois"

    ```groovy title="main.nf" linenums="7" hl_lines="11-12"
        // Cria um objeto Path a partir de uma string de caminho
        myFile = file('data/patientA_rep1_normal_R1_001.fastq.gz')

        // Imprime atributos do arquivo
        println "File object class: ${myFile.class}"
        println "File name: ${myFile.name}"
        println "Simple name: ${myFile.simpleName}"
        println "Extension: ${myFile.extension}"
        println "Parent directory: ${myFile.parent}"

        // Conta as linhas no arquivo
        COUNT_LINES(myFile)
    ```

=== "Antes"

    ```groovy title="main.nf" linenums="7" hl_lines="4-9"
        // Cria um objeto Path a partir de uma string de caminho
        myFile = file('data/patientA_rep1_normal_R1_001.fastq.gz')

        // Imprime atributos do arquivo
        println "File object class: ${myFile.class}"
        println "File name: ${myFile.name}"
        println "Simple name: ${myFile.simpleName}"
        println "Extension: ${myFile.extension}"
        println "Parent directory: ${myFile.parent}"
    ```

E agora execute o fluxo de trabalho:

```bash
nextflow run main.nf
```

??? success "Saída do comando"

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

Isso mostra que somos capazes de operar no arquivo apropriadamente dentro de um processo.

Especificamente, Nextflow realizou as seguintes operações com sucesso:

- Staged o arquivo no diretório de trabalho
- Descomprimiu o arquivo .gz
- Contou as linhas (40 linhas neste caso)
- Completou sem erro

A chave para esta operação suave é que estamos explicitamente dizendo ao Nextflow que nossa entrada é um arquivo e deve ser tratada como tal.

### 1.5. Solucione erros básicos de entrada de arquivo

Isso frequentemente confunde recém-chegados ao Nextflow, então vamos dedicar alguns minutos para ver o que acontece quando você faz isso errado.

Há dois lugares principais onde você pode errar no manuseio de arquivos: no nível do fluxo de trabalho e no nível do processo.

#### 1.5.1. Erro no nível do fluxo de trabalho

Vamos ver o que acontece se voltarmos a tratar o arquivo como uma string quando especificamos a entrada no bloco workflow.

Faça as seguintes edições no fluxo de trabalho, certificando-se de comentar as declarações print específicas de caminho:

=== "Depois"

    ```groovy title="main.nf" linenums="7" hl_lines="2 6-11"
        // Cria um objeto Path a partir de uma string de caminho
        myFile = 'data/patientA_rep1_normal_R1_001.fastq.gz'

        // Imprime atributos do arquivo
        println "File object class: ${myFile.class}"
        /*
        println "File name: ${myFile.name}"
        println "Simple name: ${myFile.simpleName}"
        println "Extension: ${myFile.extension}"
        println "Parent directory: ${myFile.parent}"
        */

        // Conta as linhas no arquivo
        COUNT_LINES(myFile)
    ```

=== "Antes"

    ```groovy title="main.nf" linenums="7" hl_lines="4-9"
        // Cria um objeto Path a partir de uma string de caminho
        myFile = file('data/patientA_rep1_normal_R1_001.fastq.gz')

        // Imprime atributos do arquivo
        println "File object class: ${myFile.class}"
        println "File name: ${myFile.name}"
        println "Simple name: ${myFile.simpleName}"
        println "Extension: ${myFile.extension}"
        println "Parent directory: ${myFile.parent}"

        // Conta as linhas no arquivo
        COUNT_LINES(myFile)
    ```

E agora execute o fluxo de trabalho:

```bash
nextflow run main.nf
```

??? failure "Saída do comando"

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

Esta é a parte importante:

```console
Not a valid path value: 'data/patientA_rep1_normal_R1_001.fastq.gz'
```

Quando você especifica uma entrada `path`, Nextflow valida que você está passando referências de arquivo reais, não apenas strings.
Este erro está dizendo que `'data/patientA_rep1_normal_R1_001.fastq.gz'` não é um valor de caminho válido porque é uma string, não um objeto Path.

Nextflow detectou imediatamente o problema e parou antes mesmo de iniciar o processo.

#### 1.5.2. Erro no nível do processo

O outro lugar onde podemos esquecer de especificar que queremos que Nextflow trate a entrada como um arquivo é na definição do processo.

!!! warning "Mantenha o erro do fluxo de trabalho de 1.5.1"

    Para este teste funcionar corretamente, mantenha o fluxo de trabalho em seu estado quebrado (usando uma string simples em vez de `file()`).
    Quando combinado com `val` no processo, isso produz o erro mostrado abaixo.

Faça a seguinte edição no módulo:

=== "Depois"

    ```groovy title="modules/count_lines.nf" linenums="3" hl_lines="5"
    process COUNT_LINES {
        debug true

        input:
        val input_file
    ```

=== "Antes"

    ```groovy title="modules/count_lines.nf" linenums="3" hl_lines="5"
    process COUNT_LINES {
        debug true

        input:
        path input_file
    ```

E agora execute o fluxo de trabalho novamente:

```bash
nextflow run main.nf
```

??? failure "Saída do comando"

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

Isso mostra muitos detalhes sobre o erro porque o processo está configurado para produzir informações de depuração, como notado acima.

Estas são as seções mais relevantes:

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

Isso diz que o sistema não conseguiu encontrar o arquivo; no entanto, se você procurar o caminho, há um arquivo com esse nome naquele local.

Quando executamos isso, Nextflow passou o valor da string para o script, mas não _staged_ o arquivo real no diretório de trabalho.
Então o processo tentou usar a string relativa, `data/patientA_rep1_normal_R1_001.fastq.gz`, mas esse arquivo não existe dentro do diretório de trabalho do processo.

Tomados em conjunto, esses dois exemplos mostram como é importante dizer ao Nextflow se uma entrada deve ser tratada como um arquivo.

!!! note

    Certifique-se de voltar e corrigir ambos os erros intencionais antes de continuar para a próxima seção.

### Conclusão

- Strings de caminho vs objetos Path: Strings são apenas texto, objetos Path são referências de arquivo inteligentes
- O método `file()` converte uma string de caminho em um objeto Path com o qual Nextflow pode trabalhar
- Você pode acessar propriedades de arquivo como `name`, `simpleName`, `extension` e `parent` [usando atributos de arquivo](https://www.nextflow.io/docs/latest/working-with-files.html#getting-file-attributes)
- Usar objetos Path em vez de strings permite ao Nextflow gerenciar arquivos adequadamente em seu fluxo de trabalho
- Resultados de Entrada de Processo: O manuseio adequado de arquivos requer objetos Path, não strings, para garantir que os arquivos sejam corretamente staged e acessíveis para uso pelos processos.

---

## 2. Usando arquivos remotos

Uma das principais características do Nextflow é a capacidade de alternar perfeitamente entre arquivos locais (na mesma máquina) e arquivos remotos acessíveis pela internet.

Se você estiver fazendo certo, nunca precisará mudar a lógica do seu fluxo de trabalho para acomodar arquivos vindos de diferentes locais.
Tudo o que você precisa fazer para usar um arquivo remoto é especificar o prefixo apropriado no caminho do arquivo quando você está fornecendo-o ao fluxo de trabalho.

Por exemplo, `/path/to/data` não tem prefixo, indicando que é um caminho de arquivo local 'normal', enquanto `s3://path/to/data` inclui o prefixo `s3://`, indicando que está localizado no armazenamento de objetos S3 da Amazon.

Muitos protocolos diferentes são suportados:

- HTTP(S)/FTP (http://, https://, ftp://)
- Amazon S3 (s3://)
- Azure Blob Storage (az://)
- Google Cloud Storage (gs://)

Para usar qualquer um destes, simplesmente especifique o prefixo relevante na string, que então é tecnicamente chamada de Uniform Resource Identifier (URI) em vez de caminho de arquivo.
Nextflow cuidará da autenticação e do staging dos arquivos no lugar certo, fazendo download ou upload e todas as outras operações de arquivo que você esperaria.

A força chave deste sistema é que ele nos permite alternar entre ambientes sem mudar nenhuma lógica do pipeline.
Por exemplo, você pode desenvolver com um conjunto de teste pequeno e local antes de mudar para um conjunto de teste em escala completa localizado em armazenamento remoto simplesmente mudando a URI.

### 2.1. Use um arquivo da internet

Vamos testar isso mudando o caminho local que estamos fornecendo ao nosso fluxo de trabalho com um caminho HTTPS apontando para uma cópia dos mesmos dados que está armazenada no Github.

!!! warning

    Isso só funcionará se você tiver uma conexão ativa com a internet.

Abra `main.nf` novamente e mude o caminho de entrada da seguinte forma:

=== "Depois"

    ```groovy title="main.nf" linenums="2" hl_lines="2"
        // Usando um arquivo remoto da internet
        myFile = file('https://raw.github.com/nextflow-io/training/master/side-quests/working_with_files/data/patientA_rep1_normal_R1_001.fastq.gz')

        // Imprime atributos do arquivo
        println "File object class: ${myFile.class}"
        println "File name: ${myFile.name}"
        println "Simple name: ${myFile.simpleName}"
        println "Extension: ${myFile.extension}"
        println "Parent directory: ${myFile.parent}"
    ```

=== "Antes"

    ```groovy title="main.nf" linenums="2" hl_lines="2"
        // Cria um objeto Path a partir de uma string de caminho
        myFile = file('data/patientA_rep1_normal_R1_001.fastq.gz')

        // Imprime atributos do arquivo
        println "File object class: ${myFile.class}"
        println "File name: ${myFile.name}"
        println "Simple name: ${myFile.simpleName}"
        println "Extension: ${myFile.extension}"
        println "Parent directory: ${myFile.parent}"
    ```

Vamos executar o fluxo de trabalho:

```bash
nextflow run main.nf
```

??? success "Saída do comando"

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

Funciona! Você pode ver que muito pouco mudou.

A única diferença na saída do console é que a classe do objeto path agora é `nextflow.file.http.XPath`, enquanto para o caminho local a classe era `sun.nio.fs.UnixPath`.
Você não precisa lembrar essas classes; mencionamos isso apenas para demonstrar que Nextflow identifica e lida com as diferentes localizações apropriadamente.

Nos bastidores, Nextflow baixou o arquivo para um diretório de staging localizado dentro do diretório work.
Esse arquivo staged pode então ser tratado como um arquivo local e criar um link simbólico no diretório do processo relevante.

Você pode verificar que isso aconteceu aqui olhando o conteúdo do diretório de trabalho localizado no valor de hash do processo.

??? abstract "Conteúdo do diretório work"

    Se o hash do processo fosse `8a/2ab7ca`, você poderia explorar o diretório work:

    ```console
    $ ls -la work/8a/2ab7ca*/
    total 16
    drwxr-xr-x  6 user  staff   192 Jan 28 10:00 .
    drwxr-xr-x  3 user  staff    96 Jan 28 10:00 ..
    -rw-r--r--  1 user  staff     0 Jan 28 10:00 .command.begin
    -rw-r--r--  1 user  staff   127 Jan 28 10:00 .command.sh
    lrwxr-xr-x  1 user  staff    89 Jan 28 10:00 patientA_rep1_normal_R1_001.fastq.gz -> /path/to/work/stage/.../patientA_rep1_normal_R1_001.fastq.gz
    ```

    O link simbólico aponta para uma cópia staged do arquivo remoto que Nextflow baixou automaticamente.

Note que para arquivos maiores, o passo de download levará algum tempo extra comparado à execução em arquivos locais.
No entanto, Nextflow verifica se já tem uma cópia staged para evitar downloads desnecessários.
Então se você executar novamente no mesmo arquivo e não tiver deletado o arquivo staged, Nextflow usará a cópia staged.

Isso mostra como é fácil alternar entre dados locais e remotos usando Nextflow, que é uma característica chave do Nextflow.

!!! note

    A única exceção importante a este princípio é que você não pode usar padrões glob ou caminhos de diretório com HTTPS porque HTTPS não pode listar múltiplos arquivos, então você deve especificar URLs exatas de arquivo.
    No entanto, outros protocolos de armazenamento como blob storage (`s3://`, `az://`, `gs://`) podem usar tanto globs quanto caminhos de diretório.

    Aqui está como você poderia usar padrões glob com armazenamento em nuvem:

    ```groovy title="Exemplos de armazenamento em nuvem (não executáveis neste ambiente)"
    // S3 com padrões glob - corresponderia a múltiplos arquivos
    ch_s3_files = channel.fromPath('s3://my-bucket/data/*.fastq.gz')

    // Azure Blob Storage com padrões glob
    ch_azure_files = channel.fromPath('az://container/data/patient*_R{1,2}.fastq.gz')

    // Google Cloud Storage com padrões glob
    ch_gcs_files = channel.fromPath('gs://bucket/data/sample_*.fastq.gz')
    ```

    Vamos mostrar como trabalhar com globs na prática na próxima seção.

### 2.2. Volte para o arquivo local

Vamos voltar a usar nossos arquivos de exemplo locais para o resto desta side quest, então vamos mudar a entrada do fluxo de trabalho de volta para o arquivo original:

=== "Depois"

    ```groovy title="main.nf" linenums="2" hl_lines="2"
        // Cria um objeto Path a partir de uma string de caminho
        myFile = file('data/patientA_rep1_normal_R1_001.fastq.gz')

        // Imprime atributos do arquivo
        println "File object class: ${myFile.class}"
        println "File name: ${myFile.name}"
        println "Simple name: ${myFile.simpleName}"
        println "Extension: ${myFile.extension}"
        println "Parent directory: ${myFile.parent}"
    ```

=== "Antes"

    ```groovy title="main.nf" linenums="2" hl_lines="2"
        // Cria um objeto Path a partir de uma string de caminho
        myFile = file('https://raw.github.com/nextflow-io/training/master/side-quests/working_with_files/data/patientA_rep1_normal_R1_001.fastq.gz')

        // Imprime atributos do arquivo
        println "File object class: ${myFile.class}"
        println "File name: ${myFile.name}"
        println "Simple name: ${myFile.simpleName}"
        println "Extension: ${myFile.extension}"
        println "Parent directory: ${myFile.parent}"
    ```

### Conclusão

- Dados remotos são acessados usando uma URI (HTTP, FTP, S3, Azure, Google Cloud)
- Nextflow automaticamente fará download e stage dos dados no lugar certo, desde que esses caminhos sejam alimentados em processos
- Não escreva lógica para baixar ou fazer upload de arquivos remotos!
- Arquivos locais e remotos produzem tipos de objetos diferentes mas funcionam identicamente
- **Importante**: HTTP/HTTPS só funcionam com arquivos únicos (sem padrões glob)
- Armazenamento em nuvem (S3, Azure, GCS) suporta tanto arquivos únicos quanto padrões glob
- Você pode alternar perfeitamente entre fontes de dados locais e remotas sem mudar a lógica do código (desde que o protocolo suporte suas operações requeridas)

---

## 3. Usando o channel factory `fromPath()`

Até agora temos trabalhado com um arquivo por vez, mas no Nextflow, tipicamente vamos querer criar um canal de entrada com múltiplos arquivos de entrada para processar.

Uma forma ingênua de fazer isso seria combinar o método `file()` com [`channel.of()`](https://www.nextflow.io/docs/latest/reference/channel.html#of) assim:

```groovy title="Exemplo de sintaxe"
ch_files = channel.of([file('data/patientA_rep1_normal_R1_001.fastq.gz')],
                      [file('data/patientA_rep1_normal_R1_001.fastq.gz')])
```

Isso funciona, mas é desajeitado.

!!! tip "Quando usar `file()` vs `channel.fromPath()`"

    - Use `file()` quando você precisar de um único objeto Path para manipulação direta (verificar se um arquivo existe, ler seus atributos, ou passar para uma única invocação de processo)
    - Use `channel.fromPath()` quando você precisar de um canal que pode conter múltiplos arquivos, especialmente com padrões glob, ou quando arquivos fluirão através de múltiplos processos

É aqui que [`channel.fromPath()`](https://www.nextflow.io/docs/latest/reference/channel.html#frompath) entra: um channel factory conveniente que agrupa toda a funcionalidade que precisamos para gerar um canal a partir de uma ou mais strings de arquivo estáticas, bem como padrões glob.

### 3.1. Adicione o channel factory

Vamos atualizar nosso fluxo de trabalho para usar `channel.fromPath`.

=== "Depois"

    ```groovy title="main.nf" linenums="7" hl_lines="1-3"
        // Carrega arquivos com channel.fromPath
        ch_files = channel.fromPath('data/patientA_rep1_normal_R1_001.fastq.gz')
        ch_files.view { myFile -> "Found file: $myFile" }

        // Imprime atributos do arquivo
        /* Comente isso por enquanto, voltaremos a eles!
        println "File object class: ${myFile.class}"
        println "File name: ${myFile.name}"
        println "Simple name: ${myFile.simpleName}"
        println "Extension: ${myFile.extension}"
        println "Parent directory: ${myFile.parent}"
        */

        // Conta as linhas no arquivo
        // COUNT_LINES(myFile)
    ```

=== "Antes"

    ```groovy title="main.nf" linenums="7" hl_lines="1-2"
        // Cria um objeto Path a partir de uma string de caminho
        myFile = file('data/patientA_rep1_normal_R1_001.fastq.gz')

        // Imprime atributos do arquivo
        println "File object class: ${myFile.class}"
        println "File name: ${myFile.name}"
        println "Simple name: ${myFile.simpleName}"
        println "Extension: ${myFile.extension}"
        println "Parent directory: ${myFile.parent}"

        // Conta as linhas no arquivo
        COUNT_LINES(myFile)
    ```

Também comentamos o código que imprime os atributos por enquanto, e adicionamos uma declaração `.view` para imprimir apenas o nome do arquivo.

Execute o fluxo de trabalho:

```bash
nextflow run main.nf
```

??? success "Saída do comando"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `main.nf` [grave_meucci] DSL2 - revision: b09964a583

    Found file: /workspaces/training/side-quests/working_with_files/data/patientA_rep1_normal_R1_001.fastq.gz
    ```

Como você pode ver, o caminho do arquivo está sendo carregado como um objeto do tipo `Path` no canal.
Isso é semelhante ao que `file()` teria feito, exceto que agora temos um canal no qual podemos carregar mais arquivos se quisermos.

Usar `channel.fromPath()` é uma forma conveniente de criar um novo canal populado por uma lista de arquivos.

### 3.2. Visualize atributos de arquivos no canal

Em nossa primeira passagem usando o channel factory, simplificamos o código e apenas imprimimos o nome do arquivo.

Vamos voltar a imprimir os atributos completos do arquivo:

=== "Depois"

    ```groovy title="main.nf" linenums="7" hl_lines="3-9 12"
        // Carrega arquivos com channel.fromPath
        ch_files = channel.fromPath('data/patientA_rep1_normal_R1_001.fastq.gz')
        ch_files.view { myFile ->
            println "File object class: ${myFile.class}"
            println "File name: ${myFile.name}"
            println "Simple name: ${myFile.simpleName}"
            println "Extension: ${myFile.extension}"
            println "Parent directory: ${myFile.parent}"
        }

        // Conta as linhas no arquivo
        COUNT_LINES(ch_files)
    ```

=== "Antes"

    ```groovy title="main.nf" linenums="7" hl_lines="3"
        // Carrega arquivos com channel.fromPath
        ch_files = channel.fromPath('data/patientA_rep1_normal_R1_001.fastq.gz')
        ch_files.view { myFile -> "Found file: $myFile" }

        // Conta as linhas no arquivo
        // COUNT_LINES(ch_files)
    ```

Também estamos reabilitando a chamada do processo `COUNT_LINES` para verificar que o processamento de arquivo ainda funciona corretamente com nossa abordagem baseada em canal.

Execute o fluxo de trabalho:

```bash
nextflow run main.nf
```

??? success "Saída do comando"

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

E aí está, mesmos resultados de antes, mas agora temos o arquivo em um canal, então podemos adicionar mais.

### 3.3. Usando um glob para corresponder múltiplos arquivos

Há várias formas que poderíamos carregar mais arquivos no canal.
Aqui vamos mostrar como usar padrões glob, que são uma forma conveniente de corresponder e recuperar nomes de arquivos e diretórios baseados em caracteres curinga.
O processo de corresponder esses padrões é chamado de "globbing" ou "expansão de nome de arquivo".

!!! note

    Como notado anteriormente, Nextflow suporta globbing para gerenciar arquivos de entrada e saída na maioria dos casos, exceto com caminhos de arquivo HTTPS porque HTTPS não pode listar múltiplos arquivos.

Digamos que queremos recuperar ambos os arquivos em um par de arquivos associados com um dado paciente, `patientA`:

```console
patientA_rep1_normal_R1_001.fastq.gz
patientA_rep1_normal_R2_001.fastq.gz
```

Já que a única diferença entre os nomes de arquivo é o número de replicata, _i.e._ o número após `R`, podemos usar o caractere curinga `*` para substituir o número da seguinte forma:

```console
patientA_rep1_normal_R*_001.fastq.gz
```

Esse é o padrão glob que precisamos.

Agora tudo o que precisamos fazer é atualizar o caminho do arquivo no channel factory para usar esse padrão glob da seguinte forma:

=== "Depois"

    ```groovy title="main.nf" linenums="7"
      // Carrega arquivos com channel.fromPath
        ch_files = channel.fromPath('data/patientA_rep1_normal_R*_001.fastq.gz')
    ```

=== "Antes"

    ```groovy title="main.nf" linenums="7"
      // Carrega arquivos com channel.fromPath
        ch_files = channel.fromPath('data/patientA_rep1_normal_R1_001.fastq.gz')
    ```

Nextflow reconhecerá automaticamente que este é um padrão glob e o tratará apropriadamente.

Execute o fluxo de trabalho para testar isso:

```bash
nextflow run main.nf
```

??? success "Saída do comando"

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

Como você pode ver, agora temos dois objetos Path em nosso canal, o que mostra que Nextflow fez a expansão de nome de arquivo corretamente, e carregou e processou ambos os arquivos como esperado.

Usando este método, podemos recuperar tantos ou tão poucos arquivos quanto quisermos apenas mudando o padrão glob. Se o tornássemos mais generoso, por exemplo substituindo todas as partes variáveis dos nomes de arquivo por `*` (_e.g._ `data/patient*_rep*_*_R*_001.fastq.gz`) poderíamos pegar todos os arquivos de exemplo no diretório `data`.

### Conclusão

- `channel.fromPath()` cria um canal com arquivos correspondendo a um padrão
- Cada arquivo é emitido como um elemento separado no canal
- Podemos usar um padrão glob para corresponder múltiplos arquivos
- Arquivos são automaticamente convertidos em objetos Path com atributos completos
- O método `.view()` permite inspeção do conteúdo do canal

---

## 4. Extraindo metadados básicos de nomes de arquivos

Na maioria dos domínios científicos, é muito comum ter metadados codificados nos nomes dos arquivos que contêm os dados.
Por exemplo, em bioinformática, arquivos contendo dados de sequenciamento são frequentemente nomeados de uma forma que codifica informações sobre a amostra, condição, replicata e número de read.

Se os nomes de arquivo são construídos de acordo com uma convenção consistente, você pode extrair esses metadados de maneira padronizada e usá-los no curso de sua análise.
Isso é um grande 'se', é claro, e você deve ser muito cauteloso sempre que depender da estrutura do nome do arquivo; mas a realidade é que esta abordagem é muito amplamente usada, então vamos ver como é feito no Nextflow.

No caso de nossos dados de exemplo, sabemos que os nomes de arquivo incluem metadados estruturados de forma consistente.
Por exemplo, o nome de arquivo `patientA_rep1_normal_R2_001` codifica o seguinte:

- ID do paciente: `patientA`
- ID da replicata: `rep1`
- tipo de amostra: `normal` (em oposição a `tumor`)
- conjunto de read: `R1` (em oposição a `R2`)

Vamos modificar nosso fluxo de trabalho para recuperar esta informação em três passos:

1. Recuperar o `simpleName` do arquivo, que inclui os metadados
2. Separar os metadados usando um método chamado `tokenize()`
3. Usar um mapa para organizar os metadados

!!! warning

    Você nunca deve codificar informações sensíveis em nomes de arquivos, como nomes de pacientes ou outras características identificadoras, pois isso pode comprometer a privacidade do paciente ou outras restrições de segurança relevantes.

### 4.1. Recupere o `simpleName`

O `simpleName` é um atributo de arquivo que corresponde ao nome do arquivo desprovido de seu caminho e extensão.

Faça as seguintes edições no fluxo de trabalho:

=== "Depois"

    ```groovy title="main.nf" linenums="7" hl_lines="3-6"
        // Carrega arquivos com channel.fromPath
        ch_files = channel.fromPath('data/patientA_rep1_normal_R*_001.fastq.gz')
        ch_files.map { myFile ->
            [ myFile.simpleName, myFile ]
        }
        .view()
    ```

=== "Antes"

    ```groovy title="main.nf" linenums="7" hl_lines="3-9"
        // Carrega arquivos com channel.fromPath
        ch_files = channel.fromPath('data/patientA_rep1_normal_R*_001.fastq.gz')
        ch_files.view { myFile ->
            println "File object class: ${myFile.class}"
            println "File name: ${myFile.name}"
            println "Simple name: ${myFile.simpleName}"
            println "Extension: ${myFile.extension}"
            println "Parent directory: ${myFile.parent}"
        }
    ```

Isso recupera o `simpleName` e o associa com o objeto de arquivo completo usando uma operação `map()`.

Execute o fluxo de trabalho para testar que funciona:

```bash
nextflow run main.nf
```

??? success "Saída do comando"

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

Cada elemento no canal agora é uma tupla contendo o `simpleName` e o objeto de arquivo original.

### 4.2. Extraia os metadados do `simplename`

Neste ponto, os metadados que queremos estão incorporados no `simplename`, mas não podemos acessar itens individuais diretamente.
Então precisamos dividir o `simplename` em seus componentes.
Felizmente, esses componentes são simplesmente separados por sublinhados no nome de arquivo original, então podemos aplicar um método comum do Nextflow chamado `tokenize()` que é perfeito para esta tarefa.

Faça as seguintes edições no fluxo de trabalho:

=== "Depois"

    ```groovy title="main.nf" linenums="7" hl_lines="4"
        // Carrega arquivos com channel.fromPath
        ch_files = channel.fromPath('data/patientA_rep1_normal_R*_001.fastq.gz')
        ch_files.map { myFile ->
            [ myFile.simpleName.tokenize('_'), myFile ]
        }
    ```

=== "Antes"

    ```groovy title="main.nf" linenums="7" hl_lines="4"
        // Carrega arquivos com channel.fromPath
        ch_files = channel.fromPath('data/patientA_rep1_normal_R*_001.fastq.gz')
        ch_files.map { myFile ->
            [ myFile.simpleName, myFile ]
        }
    ```

O método `tokenize()` dividirá a string `simpleName` onde quer que encontre sublinhados, e retornará uma lista contendo as substrings.

Execute o fluxo de trabalho:

```bash
nextflow run main.nf
```

??? success "Saída do comando"

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

Agora a tupla para cada elemento em nosso canal contém a lista de metadados (_e.g._ `[patientA, rep1, normal, R1, 001]`) e o objeto de arquivo original.

Isso é ótimo!
Quebramos nossa informação de paciente de uma única string em uma lista de strings.
Agora podemos lidar com cada parte da informação do paciente separadamente.

### 4.3. Use um mapa para organizar os metadados

Nossos metadados são apenas uma lista plana no momento.
É fácil o suficiente de usar mas difícil de ler.

```console
[patientA, rep1, normal, R1, 001]
```

O que é o item no índice 3? Você pode dizer sem se referir de volta à explicação original da estrutura de metadados?

Esta é uma grande oportunidade para usar um armazenamento de chave-valor, onde cada item tem um conjunto de chaves e seus valores associados, então você pode facilmente se referir a cada chave para obter o valor correspondente.

Em nosso exemplo, isso significa ir desta organização:

```groovy
data = [patientA, 1, normal, R1]

println data[3]
```

Para esta:

```groovy
data = [id: patientA, replicate: 1, type: normal, readNum: 1]

println data.readNum
```

No Nextflow, isso é chamado de [mapa](https://nextflow.io/docs/latest/script.html#maps).

Vamos converter nossa lista plana em um mapa agora.
Faça as seguintes edições no fluxo de trabalho:

=== "Depois"

    ```groovy title="main.nf" linenums="7" hl_lines="4-13"
        // Carrega arquivos com channel.fromPath
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

=== "Antes"

    ```groovy title="main.nf" linenums="7" hl_lines="4"
        // Carrega arquivos com channel.fromPath
        ch_files = channel.fromPath('data/patientA_rep1_normal_R*_001.fastq.gz')
        ch_files.map { myFile ->
            [ myFile.simpleName.tokenize('_'), myFile ]
        }
    ```

As mudanças chave aqui são:

- **Atribuição desestruturante**: `def (patient, replicate, type, readNum) = ...` extrai os valores tokenizados em variáveis nomeadas em uma linha
- **Sintaxe literal de mapa**: `[id: patient, replicate: ...]` cria um mapa onde cada chave (como `id`) é associada com um valor (como `patient`)
- **Estrutura aninhada**: A lista externa `[..., myFile]` emparelha o mapa de metadados com o objeto de arquivo original

Também simplificamos algumas das strings de metadados usando um método de substituição de string chamado `replace()` para remover alguns caracteres que são desnecessários (_e.g._ `replicate.replace('rep', '')` para manter apenas o número dos IDs de replicata).

Vamos executar o fluxo de trabalho novamente:

```bash
nextflow run main.nf
```

??? success "Saída do comando"

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

Agora os metadados estão bem rotulados (_e.g._ `[id:patientA, replicate:1, type:normal, readNum:2]`) então é muito mais fácil dizer o que é o quê.

Também será muito mais fácil realmente fazer uso de elementos de metadados no fluxo de trabalho, e tornará nosso código mais fácil de ler e mais sustentável.

### Conclusão

- Podemos lidar com nomes de arquivos no Nextflow com o poder de uma linguagem de programação completa
- Podemos tratar os nomes de arquivos como strings para extrair informação relevante
- Uso de métodos como `tokenize()` e `replace()` permite manipular strings no nome do arquivo
- A operação `.map()` transforma elementos de canal enquanto preserva a estrutura
- Metadados estruturados (mapas) tornam o código mais legível e sustentável do que listas posicionais

A seguir, veremos como lidar com arquivos de dados pareados.

---

## 5. Lidando com arquivos de dados pareados

Muitos desenhos experimentais produzem arquivos de dados pareados que se beneficiam de serem tratados de forma explicitamente pareada.
Por exemplo, em bioinformática, dados de sequenciamento são frequentemente gerados na forma de reads pareados, significando strings de sequência que se originam do mesmo fragmento de DNA (frequentemente chamados de 'forward' e 'reverse' porque são lidos de extremidades opostas).

Esse é o caso de nossos dados de exemplo, onde R1 e R2 se referem aos dois conjuntos de reads.

```console
data/patientA_rep1_normal_R1_001.fastq.gz
data/patientA_rep1_normal_R2_001.fastq.gz
```

Nextflow fornece um channel factory especializado para trabalhar com arquivos pareados como este chamado `channel.fromFilePairs()`, que automaticamente agrupa arquivos baseado em um padrão de nomeação compartilhado. Isso permite que você associe os arquivos pareados mais estreitamente com menos esforço.

Vamos modificar nosso fluxo de trabalho para aproveitar isso.
Vai levar dois passos:

1. Mudar o channel factory para `channel.fromFilePairs()`
2. Extrair e mapear os metadados

### 5.1. Mude o channel factory para `channel.fromFilePairs()`

Para usar `channel.fromFilePairs`, precisamos especificar o padrão que Nextflow deve usar para identificar os dois membros em um par.

Voltando aos nossos dados de exemplo, podemos formalizar o padrão de nomeação da seguinte forma:

```console
data/patientA_rep1_normal_R{1,2}_001.fastq.gz
```

Isso é semelhante ao padrão glob que usamos anteriormente, exceto que isso enumera especificamente as substrings (seja `1` ou `2` vindo logo após o R) que identificam os dois membros do par.

Vamos atualizar o fluxo de trabalho `main.nf` de acordo:

=== "Depois"

    ```groovy title="main.nf" linenums="7" hl_lines="1-2"
        // Carrega arquivos com channel.fromFilePairs
        ch_files = channel.fromFilePairs('data/patientA_rep1_normal_R{1,2}_001.fastq.gz')
        /* Comente o mapeamento por enquanto, voltaremos a ele!
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

=== "Antes"

    ```groovy title="main.nf" linenums="7" hl_lines="1-2"
        // Carrega arquivos com channel.fromFilePairs
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

Mudamos o channel factory e adaptamos o padrão de correspondência de arquivos, e enquanto fazíamos isso, comentamos a operação map.
Adicionaremos isso de volta mais tarde, com algumas modificações.

Execute o fluxo de trabalho para testá-lo:

```bash
nextflow run main.nf
```

??? failure "Saída do comando"

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

Uh-oh, desta vez a execução falhou!

A parte relevante da mensagem de erro está aqui:

```console
Not a valid path value: 'patientA_rep1_normal_R'
```

Isso é porque mudamos o channel factory.
Até agora, o canal de entrada original continha apenas os caminhos de arquivo.
Toda a manipulação de metadados que temos feito não afetou realmente o conteúdo do canal.

Agora que estamos usando o channel factory `.fromFilePairs`, o conteúdo do canal resultante é diferente.
Vemos apenas um elemento de canal, composto de uma tupla contendo dois itens: a parte do `simpleName` compartilhada pelos dois arquivos, que serve como um identificador, e uma tupla contendo os dois objetos de arquivo, no formato `id, [ file1, file2 ]`.

Isso é ótimo, porque Nextflow fez o trabalho difícil de extrair o nome do paciente examinando o prefixo compartilhado e usando-o como um identificador de paciente.

No entanto, isso quebra nosso fluxo de trabalho atual.
Se quiséssemos ainda executar `COUNT_LINES` da mesma forma sem mudar o processo, teríamos que aplicar uma operação de mapeamento para extrair os caminhos de arquivo.
Mas não vamos fazer isso, porque nosso objetivo final é usar um processo diferente, `ANALYZE_READS`, que lida com pares de arquivos apropriadamente.

Então vamos simplesmente comentar (ou deletar) a chamada para `COUNT_LINES` e seguir em frente.

=== "Depois"

    ```groovy title="main.nf" linenums="26" hl_lines="2"
        // Conta as linhas no arquivo
        // COUNT_LINES(ch_files)
    ```

=== "Antes"

    ```groovy title="main.nf" linenums="26" hl_lines="2"
        // Conta as linhas no arquivo
        COUNT_LINES(ch_files)
    ```

Você também pode comentar ou deletar a declaração include de `COUNT_LINES`, mas isso não terá efeito funcional.

Agora vamos executar o fluxo de trabalho novamente:

```bash
nextflow run main.nf
```

??? success "Saída do comando"

    ```console hl_lines="5"
     N E X T F L O W   ~  version 25.10.2

    Launching `main.nf` [fabulous_davinci] DSL2 - revision: 22b53268dc

    [patientA_rep1_normal_R, [/workspaces/training/side-quests/working_with_files/data/patientA_rep1_normal_R1_001.fastq.gz, /workspaces/training/side-quests/working_with_files/data/patientA_rep1_normal_R2_001.fastq.gz]]
    ```

Yay, desta vez o fluxo de trabalho tem sucesso!

No entanto, ainda precisamos extrair o resto dos metadados do campo `id`.

### 5.2. Extraia e organize metadados de pares de arquivos

Nossa operação `map` de antes não funcionará porque não corresponde à estrutura de dados, mas podemos modificá-la para funcionar.

Já temos acesso ao identificador de paciente real na string que `fromFilePairs()` usou como identificador, então podemos usar isso para extrair os metadados sem obter o `simpleName` do objeto Path como fizemos antes.

Descomente a operação map no fluxo de trabalho e faça as seguintes edições:

=== "Depois"

    ```groovy title="main.nf" linenums="7" hl_lines="3-4 9 11 13"
        // Carrega arquivos com channel.fromFilePairs
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

=== "Antes"

    ```groovy title="main.nf" linenums="7" hl_lines="3-5 11 13"
        // Carrega arquivos com channel.fromFilePairs
        ch_files = channel.fromFilePairs('data/patientA_rep1_normal_R{1,2}_001.fastq.gz')
        /* Comente o mapeamento por enquanto, voltaremos a ele!
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

Desta vez o map começa de `id, files` em vez de apenas `myFile`, e `tokenize()` é aplicado a `id` em vez de a `myFile.simpleName`.

Note também que descartamos `readNum` da linha `tokenize()`; quaisquer substrings que não nomeamos especificamente (começando da esquerda) serão silenciosamente descartadas.
Podemos fazer isso porque os arquivos pareados agora estão estreitamente associados, então não precisamos mais de `readNum` no mapa de metadados.

Vamos executar o fluxo de trabalho:

```bash
nextflow run main.nf
```

??? success "Saída do comando"

    ```console

    N E X T F L O W   ~  version 25.10.2

    Launching `main.nf` [prickly_stonebraker] DSL2 - revision: f62ab10a3f

    [[id:patientA, replicate:1, type:normal], [/workspaces/training/side-quests/working_with_files/data/patientA_rep1_normal_R1_001.fastq.gz, /workspaces/training/side-quests/working_with_files/data/patientA_rep1_normal_R2_001.fastq.gz]]
    ```

E aí está: temos o mapa de metadados (`[id:patientA, replicate:1, type:normal]`) na primeira posição da tupla de saída, seguido pela tupla de arquivos pareados, como pretendido.

É claro, isso só pegará e processará aquele par específico de arquivos.
Se você quiser experimentar processar múltiplos pares, pode tentar adicionar curingas no padrão de entrada e ver o que acontece.
Por exemplo, tente usar `data/patientA_rep1_*_R{1,2}_001.fastq.gz`

### Conclusão

- [`channel.fromFilePairs()` automaticamente encontra e emparelha arquivos relacionados](https://www.nextflow.io/docs/latest/reference/channel.html#fromfilepairs)
- Isso simplifica o manuseio de reads paired-end em seu pipeline
- Arquivos pareados podem ser agrupados como tuplas `[id, [file1, file2]]`
- A extração de metadados pode ser feita a partir do ID de arquivo pareado em vez de arquivos individuais

---

## 6. Usando operações de arquivo em processos

Agora vamos juntar tudo isso em um processo simples para reforçar como usar operações de arquivo dentro de um processo Nextflow.

Fornecemos um módulo de processo pré-escrito chamado `ANALYZE_READS` que recebe uma tupla de metadados e um par de arquivos de entrada e os analisa.
Poderíamos imaginar que isso está fazendo alinhamento de sequências, ou chamada de variantes ou qualquer outro passo que faça sentido para este tipo de dado.

Vamos começar.

### 6.1. Importe o processo e examine o código

Para usar este processo no fluxo de trabalho, apenas precisamos adicionar uma declaração include de módulo antes do bloco workflow.

Faça a seguinte edição no fluxo de trabalho:

=== "Depois"

    ```groovy title="main.nf" linenums="1" hl_lines="3"
    #!/usr/bin/env nextflow

    include { ANALYZE_READS } from './modules/analyze_reads.nf'

    workflow {
    ```

=== "Antes"

    ```groovy title="main.nf" linenums="1"
    #!/usr/bin/env nextflow

    workflow {
    ```

Você pode abrir o arquivo de módulo para examinar seu código:

```groovy title="modules/analyze_reads.nf - exemplo de processo" linenums="1"
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

!!! note

    As diretivas `tag` e `publishDir` usam sintaxe de closure (`{ ... }`) em vez de interpolação de string (`"${...}"`).
    Isso é porque essas diretivas referenciam variáveis de entrada (`meta`) que não estão disponíveis até o tempo de execução.
    A sintaxe de closure adia a avaliação até que o processo realmente execute.

!!! note

    Estamos chamando nosso mapa de metadados de `meta` por convenção.
    Para um mergulho mais profundo em meta maps, veja a side quest [Metadata and meta maps](./metadata.md).

### 6.2. Chame o processo no fluxo de trabalho

Agora que o processo está disponível para o fluxo de trabalho, podemos adicionar uma chamada ao processo `ANALYZE_READS` para executá-lo.

Para executá-lo em nossos dados de exemplo, precisaremos fazer duas coisas:

1. Dar um nome ao canal remapeado
2. Adicionar uma chamada ao processo

#### 6.2.1. Nomeie o canal de entrada remapeado

Anteriormente aplicamos as manipulações de mapeamento diretamente ao canal de entrada.
Para alimentar os conteúdos remapeados ao processo `ANALYZE_READS` (e fazer isso de uma forma que seja clara e fácil de ler) queremos criar um novo canal chamado `ch_samples`.

Podemos fazer isso usando o operador [`set`](https://www.nextflow.io/docs/latest/reference/operator.html#set).

No fluxo de trabalho principal, substitua o operador `.view()` por `.set { ch_samples }`, e adicione uma linha testando que podemos nos referir ao canal por nome.

=== "Depois"

    ```groovy title="main.nf" linenums="7" hl_lines="14 16-17"
        // Carrega arquivos com channel.fromFilePairs
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

        // Temporário: espie dentro de ch_samples
        ch_samples.view()
    ```

=== "Antes"

    ```groovy title="main.nf" linenums="7" hl_lines="14"
        // Carrega arquivos com channel.fromFilePairs
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

Vamos executar isso:

```bash
nextflow run main.nf
```

??? success "Saída do comando"

    ```console
     N E X T F L O W   ~  version 25.10.2

    Launching `main.nf` [goofy_kirch] DSL2 - revision: 3313283e42

    [[id:patientA, replicate:1, type:normal], [/workspaces/training/side-quests/working_with_files/data/patientA_rep1_normal_R1_001.fastq.gz, /workspaces/training/side-quests/working_with_files/data/patientA_rep1_normal_R2_001.fastq.gz]]
    ```

Isso confirma que agora podemos nos referir ao canal por nome.

#### 6.2.2. Chame o processo nos dados

Agora vamos realmente chamar o processo `ANALYZE_READS` no canal `ch_samples`.

No fluxo de trabalho principal, faça as seguintes mudanças de código:

=== "Depois"

    ```groovy title="main.nf" linenums="23"
        // Executa a análise
        ANALYZE_READS(ch_samples)
    ```

=== "Antes"

    ```groovy title="main.nf" linenums="23"
        // Temporário: espie dentro de ch_samples
        ch_samples.view()
    ```

Vamos executar isso:

```bash
nextflow run main.nf
```

??? success "Saída do comando"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `./main.nf` [shrivelled_cori] DSL2 - revision: b546a31769

    executor >  local (1)
    [b5/110360] process > ANALYZE_READS (patientA) [100%] 1 of 1 ✔
    ```

Este processo está configurado para publicar suas saídas em um diretório `results`, então dê uma olhada lá.

??? abstract "Conteúdo do diretório e arquivo"

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

O processo pegou nossas entradas e criou um novo arquivo contendo os metadados do paciente, como projetado.
Esplêndido!

### 6.3. Inclua muitos mais pacientes

É claro, isso está apenas processando um único par de arquivos para um único paciente, o que não é exatamente o tipo de alto desempenho que você espera obter com Nextflow.
Você provavelmente vai querer processar muito mais dados de uma vez.

Lembre-se de que `channel.fromPath()` aceita um _glob_ como entrada, o que significa que pode aceitar qualquer número de arquivos que correspondam ao padrão.
Portanto, se quisermos incluir todos os pacientes, podemos simplesmente modificar a string de entrada para incluir mais pacientes, como notado de passagem anteriormente.

Vamos fingir que queremos ser o mais ganancioso possível.
Faça as seguintes edições no fluxo de trabalho:

=== "Depois"

    ```groovy title="main.nf" linenums="7" hl_lines="2"
        // Carrega arquivos com channel.fromFilePairs
        ch_files = channel.fromFilePairs('data/*_R{1,2}_001.fastq.gz')
    ```

=== "Antes"

    ```groovy title="main.nf" linenums="7" hl_lines="2"
        // Carrega arquivos com channel.fromFilePairs
        ch_files = channel.fromFilePairs('data/patientA_rep1_normal_R{1,2}_001.fastq.gz')
    ```

Execute o pipeline novamente:

```bash
nextflow run main.nf
```

??? success "Saída do comando"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `./main.nf` [big_stonebraker] DSL2 - revision: f7f9b8a76c

    executor >  local (8)
    [d5/441891] process > ANALYZE_READS (patientC) [100%] 8 of 8 ✔
    ```

O diretório results agora deve conter resultados para todos os dados disponíveis.

??? abstract "Conteúdo do diretório"

    ```console
    results
    ├── patientA
    │   └── patientA_stats.txt
    ├── patientB
    │   └── patientB_stats.txt
    └── patientC
        └── patientC_stats.txt
    ```

Sucesso! Analisamos todos os pacientes de uma vez! Certo?

Talvez não.
Se você olhar mais de perto, temos um problema: temos duas replicatas para patientA, mas apenas um arquivo de saída!
Estamos sobrescrevendo o arquivo de saída cada vez.

### 6.4. Torne os arquivos publicados únicos

Já que temos acesso aos metadados do paciente, podemos usá-los para tornar os arquivos publicados únicos incluindo metadados diferenciadores, seja na estrutura de diretórios ou nos próprios nomes de arquivo.

Faça a seguinte mudança no fluxo de trabalho:

=== "Depois"

    ```groovy title="modules/analyze_reads.nf" linenums="6"
        publishDir { "results/${meta.type}/${meta.id}/${meta.replicate}" }, mode: 'copy'
    ```

=== "Antes"

    ```groovy title="modules/analyze_reads.nf" linenums="6"
        publishDir { "results/${meta.id}" }, mode: 'copy'
    ```

Aqui mostramos a opção de usar níveis de diretório adicionais para considerar tipos de amostra e replicatas, mas você poderia experimentar fazer isso no nível de nome de arquivo também.

Agora execute o pipeline mais uma vez, mas certifique-se de remover o diretório results primeiro para dar a si mesmo um espaço de trabalho limpo:

```bash
rm -r results
nextflow run main.nf
```

??? success "Saída do comando"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `./main.nf` [insane_swartz] DSL2 - revision: fff18abe6d

    executor >  local (8)
    [e3/449081] process > ANALYZE_READS (patientC) [100%] 8 of 8 ✔
    ```

Verifique o diretório results agora:

??? abstract "Conteúdo do diretório"

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

E aí está, todos os nossos metadados, bem organizados. Isso é sucesso!

Há muito mais que você pode fazer uma vez que tenha seus metadados carregados em um mapa como este:

1. Criar diretórios de saída organizados baseados em atributos do paciente
2. Tomar decisões em processos baseadas em propriedades do paciente
3. Dividir, juntar e recombinar dados baseado em valores de metadados

Este padrão de manter metadados explícitos e anexados aos dados (em vez de codificados em nomes de arquivo) é uma prática recomendada central no Nextflow que permite construir fluxos de trabalho de análise robustos e mantíveis.
Você pode aprender mais sobre isso na side quest [Metadata and meta maps](./metadata.md).

### Conclusão

- A diretiva `publishDir` pode organizar saídas baseadas em valores de metadados
- Metadados em tuplas permitem organização estruturada de resultados
- Esta abordagem cria fluxos de trabalho mantíveis com proveniência de dados clara
- Processos podem receber tuplas de metadados e arquivos como entrada
- A diretiva `tag` fornece identificação de processo nos logs de execução
- A estrutura do fluxo de trabalho separa criação de canal de execução de processo

---

## Resumo

Nesta side quest, você aprendeu como trabalhar com arquivos no Nextflow, desde operações básicas até técnicas mais avançadas para lidar com coleções de arquivos.

Aplicar essas técnicas em seu próprio trabalho permitirá que você construa fluxos de trabalho mais eficientes e mantíveis, especialmente ao trabalhar com grandes números de arquivos com convenções de nomeação complexas.

### Padrões-chave

1.  **Operações Básicas de Arquivo:** Criamos objetos Path com `file()` e acessamos atributos de arquivo como nome, extensão e diretório pai, aprendendo a diferença entre strings e objetos Path.

    - Crie um objeto Path com `file()`

    ```groovy
    myFile = file('path/to/file.txt')
    ```

    - Obtenha atributos de arquivo

    ```groovy
    println myFile.name       // file.txt
    println myFile.baseName   // file
    println myFile.extension  // txt
    println myFile.parent     // path/to
    ```

2.  **Usando Arquivos Remotos**: Aprendemos como alternar transparentemente entre arquivos locais e remotos usando URIs, demonstrando a capacidade do Nextflow de lidar com arquivos de várias fontes sem mudar a lógica do fluxo de trabalho.

    - Arquivo local

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

3.  **Carregando arquivos usando o channel factory `fromPath()`:** Criamos canais a partir de padrões de arquivo com `channel.fromPath()` e visualizamos seus atributos de arquivo, incluindo tipos de objetos.

    - Crie um canal a partir de um padrão de arquivo

    ```groovy
     ch_files = channel.fromPath('data/*.fastq.gz')
    ```

    - Obtenha atributos de arquivo

    ```groovy
     ch_files.view { myFile ->
        println "File object class: ${myFile.class}"
        println "File name: ${myFile.name}"
        println "Simple name: ${myFile.simpleName}"
        println "Extension: ${myFile.extension}"
        println "Parent directory: ${myFile.parent}"
    }
    ```

4.  **Extraindo Metadados de Paciente de Nomes de Arquivo:** Usamos `tokenize()` e `replace()` para extrair e estruturar metadados de nomes de arquivo, convertendo-os em mapas organizados.

    ```groovy
    def name = file.name.tokenize('_')
    def patientId = name[0]
    def replicate = name[1].replace('rep', '')
    def type = name[2]
    def readNum = name[3].replace('R', '')
    ```

5.  **Simplificando com channel.fromFilePairs:** Usamos `channel.fromFilePairs()` para automaticamente emparelhar arquivos relacionados e extrair metadados de IDs de arquivo pareado.

    ```groovy
    ch_pairs = channel.fromFilePairs('data/*_R{1,2}_001.fastq.gz')
    ```

6.  **Usando Operações de Arquivo em Processos:** Integramos operações de arquivo em processos Nextflow com manuseio apropriado de entrada, usando `publishDir` para organizar saídas baseadas em metadados.

    - Associe um meta map com as entradas do processo

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

    - Organize saídas baseadas em metadados

    ```groovy
    publishDir { "results/${meta.type}/${meta.id}/${meta.replicate}" }, mode: 'copy'
    ```

### Recursos adicionais

- [Documentação Nextflow: Working with Files](https://www.nextflow.io/docs/latest/working-with-files.html)
- [channel.fromPath](https://www.nextflow.io/docs/latest/reference/channel.html#frompath)
- [channel.fromFilePairs](https://www.nextflow.io/docs/latest/reference/channel.html#fromfilepairs)

---

## O que vem a seguir?

Retorne ao [menu de Side Quests](./index.md) ou clique no botão no canto inferior direito da página para avançar para o próximo tópico da lista.
