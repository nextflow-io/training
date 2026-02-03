# Parte 3: Hello Workflow

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Tradução assistida por IA - [saiba mais e sugira melhorias](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

<!--
<div class="video-wrapper">
  <iframe width="560" height="315" src="https://www.youtube.com/embed/zJP7cUYPEbA?si=Irl9nAQniDyICp2b&amp;list=PLPZ8WHdZGxmXiHf8B26oB_fTfoKQdhlik" title="YouTube video player" frameborder="0" allow="accelerometer; autoplay; clipboard-write; encrypted-media; gyroscope; picture-in-picture; web-share" referrerpolicy="strict-origin-when-cross-origin" allowfullscreen></iframe>
</div>

/// caption
:fontawesome-brands-youtube:{ .youtube } Veja [a playlist completa](https://www.youtube.com/playlist?list=PLPZ8WHdZGxmXiHf8B26oB_fTfoKQdhlik) no canal do Nextflow no YouTube.

:green_book: A transcrição do vídeo está disponível [aqui](./transcripts/03_hello_workflow.md).
///
-->

A maioria dos fluxos de trabalho do mundo real envolve mais de uma etapa.
Neste módulo de treinamento, você aprenderá como conectar processos em um fluxo de trabalho de múltiplas etapas.

Isso ensinará a você a maneira Nextflow de realizar o seguinte:

1. Fazer os dados fluírem de um processo para o próximo
2. Coletar saídas de múltiplas chamadas de processo em uma única chamada de processo
3. Passar parâmetros adicionais para um processo
4. Gerenciar múltiplas saídas vindas de um processo

Para demonstrar, continuaremos construindo sobre o exemplo Hello World agnóstico de domínio das Partes 1 e 2.
Desta vez, faremos as seguintes alterações em nosso fluxo de trabalho para refletir melhor como as pessoas constroem fluxos de trabalho reais:

1. Adicionar uma segunda etapa que converte a saudação para maiúsculas.
2. Adicionar uma terceira etapa que coleta todas as saudações transformadas e as escreve em um único arquivo.
3. Adicionar um parâmetro para nomear o arquivo de saída final e passá-lo como uma entrada secundária para a etapa de coleta.
4. Fazer a etapa de coleta também relatar uma estatística simples sobre o que foi processado.

??? info "Como começar a partir desta seção"

    Esta seção do curso pressupõe que você completou as Partes 1-2 do curso [Hello Nextflow](./index.md), mas se você está confortável com os conceitos básicos cobertos nessas seções, pode começar a partir daqui sem fazer nada especial.

---

## 0. Aquecimento: Execute `hello-workflow.nf`

Vamos usar o script de fluxo de trabalho `hello-workflow.nf` como ponto de partida.
Ele é equivalente ao script produzido ao trabalhar na Parte 2 deste curso de treinamento, exceto que removemos as instruções `view()` e alteramos o destino de saída:

```groovy title="hello-workflow.nf" linenums="37" hl_lines="3"
output {
    first_output {
        path 'hello_workflow'
        mode 'copy'
    }
}
```

Este diagrama resume a operação atual do fluxo de trabalho.
Deve parecer familiar, exceto que agora estamos mostrando explicitamente que as saídas do processo são empacotadas em um canal, assim como as entradas eram.
Vamos dar um bom uso a esse canal de saída em um minuto.

<figure class="excalidraw">
--8<-- "docs/en/docs/hello_nextflow/img/hello-workflow-channels.svg"
</figure>

Apenas para ter certeza de que tudo está funcionando, execute o script uma vez antes de fazer quaisquer alterações:

```bash
nextflow run hello-workflow.nf
```

??? success "Saída do comando"

    ```console
     N E X T F L O W   ~  version 25.10.2

    Launching `hello-workflow.nf` [admiring_lamarr] DSL2 - revision: 4d4053520d

    executor >  local (3)
    [b1/5826b5] process > sayHello (2) [100%] 3 of 3 ✔
    ```

Como anteriormente, você encontrará os arquivos de saída no local especificado no bloco `output`.
Para este capítulo, está em `results/hello_workflow/`.

??? abstract "Conteúdo do diretório"

    ```console
    results/hello_workflow
    ├── Bonjour-output.txt
    ├── Hello-output.txt
    └── Holà-output.txt
    ```

Se isso funcionou para você, você está pronto para aprender como montar um fluxo de trabalho de múltiplas etapas.

---

## 1. Adicione uma segunda etapa ao fluxo de trabalho

Vamos adicionar uma etapa para converter cada saudação para maiúsculas.

<figure class="excalidraw">
--8<-- "docs/en/docs/hello_nextflow/img/hello-multistep.svg"
</figure>

Para isso, precisamos fazer três coisas:

- Definir o comando que vamos usar para fazer a conversão para maiúsculas.
- Escrever um novo processo que envolve o comando de conversão para maiúsculas.
- Chamar o novo processo no bloco de fluxo de trabalho e configurá-lo para receber a saída do processo `sayHello()` como entrada.

### 1.1. Defina o comando de conversão para maiúsculas e teste-o no terminal

Para fazer a conversão das saudações para maiúsculas, vamos usar uma ferramenta UNIX clássica chamada `tr` para 'text replacement', com a seguinte sintaxe:

```bash title="Sintaxe"
tr '[a-z]' '[A-Z]'
```

Esta é uma substituição de texto muito simples que não leva em conta letras acentuadas, então por exemplo 'Holà' se tornará 'HOLà', mas fará um trabalho bom o suficiente para demonstrar os conceitos do Nextflow e isso é o que importa.

Para testá-lo, podemos executar o comando `echo 'Hello World'` e direcionar sua saída para o comando `tr`:

```bash
echo 'Hello World' | tr '[a-z]' '[A-Z]' > UPPER-output.txt
```

A saída é um arquivo de texto chamado `UPPER-output.txt` que contém a versão em maiúsculas da string `Hello World`.

??? abstract "Conteúdo do arquivo"

    ```console title="UPPER-output.txt"
    HELLO WORLD
    ```

Isso é basicamente o que vamos tentar fazer com nosso fluxo de trabalho.

### 1.2. Escreva a etapa de conversão para maiúsculas como um processo Nextflow

Podemos modelar nosso novo processo no primeiro, já que queremos usar todos os mesmos componentes.

Adicione a seguinte definição de processo ao script de fluxo de trabalho, logo abaixo da primeira:

```groovy title="hello-workflow.nf" linenums="20"
/*
 * Use uma ferramenta de substituição de texto para converter a saudação para maiúsculas
 */
process convertToUpper {

    input:
    path input_file

    output:
    path "UPPER-${input_file}"

    script:
    """
    cat '$input_file' | tr '[a-z]' '[A-Z]' > 'UPPER-${input_file}'
    """
}
```

Neste, compomos o segundo nome de arquivo de saída baseado no nome do arquivo de entrada, de forma similar ao que fizemos originalmente para a saída do primeiro processo.

### 1.3. Adicione uma chamada ao novo processo no bloco de fluxo de trabalho

Agora precisamos dizer ao Nextflow para realmente chamar o processo que acabamos de definir.

No bloco de fluxo de trabalho, faça a seguinte alteração de código:

=== "Depois"

    ```groovy title="hello-workflow.nf" linenums="44" hl_lines="10-11"
    workflow {

        main:
        // cria um canal para entradas de um arquivo CSV
        greeting_ch = channel.fromPath(params.input)
                            .splitCsv()
                            .map { line -> line[0] }
        // emite uma saudação
        sayHello(greeting_ch)
        // converte a saudação para maiúsculas
        convertToUpper()

        publish:
        first_output = sayHello.out
    }
    ```

=== "Antes"

    ```groovy title="hello-workflow.nf" linenums="44"
    workflow {

        main:
        // cria um canal para entradas de um arquivo CSV
        greeting_ch = channel.fromPath(params.input)
                            .splitCsv()
                            .map { line -> line[0] }
        // emite uma saudação
        sayHello(greeting_ch)

        publish:
        first_output = sayHello.out
    }
    ```

Isso ainda não é funcional porque não especificamos o que deve ser entrada para o processo `convertToUpper()`.

### 1.4. Passe a saída do primeiro processo para o segundo processo

Agora precisamos fazer a saída do processo `sayHello()` fluir para o processo `convertToUpper()`.

Convenientemente, o Nextflow empacota automaticamente a saída de um processo em um canal, como mostrado no diagrama na seção de aquecimento.
Podemos nos referir ao canal de saída de um processo como `<process>.out`.

Então a saída do processo `sayHello` é um canal chamado `sayHello.out`, que podemos conectar diretamente na chamada para `convertToUpper()`.

<figure class="excalidraw">
--8<-- "docs/en/docs/hello_nextflow/img/hello-multistep-connector.svg"
</figure>

No bloco de fluxo de trabalho, faça a seguinte alteração de código:

=== "Depois"

    ```groovy title="hello-workflow.nf" linenums="53" hl_lines="2"
        // converte a saudação para maiúsculas
        convertToUpper(sayHello.out)
    ```

=== "Antes"

    ```groovy title="hello-workflow.nf" linenums="53" hl_lines="2"
        // converte a saudação para maiúsculas
        convertToUpper()
    ```

Para um caso simples como este (uma saída para uma entrada), isso é tudo que precisamos fazer para conectar dois processos!

### 1.5. Configure a publicação da saída do fluxo de trabalho

Finalmente, vamos atualizar as saídas do fluxo de trabalho para publicar os resultados do segundo processo também.

#### 1.5.1. Atualize a seção `publish:` do bloco `workflow`

No bloco `workflow`, faça a seguinte alteração de código:

=== "Depois"

    ```groovy title="hello-workflow.nf" linenums="56" hl_lines="3"
        publish:
        first_output = sayHello.out
        uppercased = convertToUpper.out
    }
    ```

=== "Antes"

    ```groovy title="hello-workflow.nf" linenums="56"
        publish:
        first_output = sayHello.out
    }
    ```

A lógica é a mesma de antes.

#### 1.5.2. Atualize o bloco `output`

No bloco `output`, faça a seguinte alteração de código:

=== "Depois"

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

=== "Antes"

    ```groovy title="hello-workflow.nf" linenums="61"
    output {
        first_output {
            path 'hello_workflow'
            mode 'copy'
        }
    }
    ```

Mais uma vez, a lógica é a mesma de antes.

Isso mostra que você pode controlar as configurações de saída em um nível muito granular, para cada saída individual.
Sinta-se à vontade para tentar mudar os caminhos ou o modo de publicação para um dos processos para ver o que acontece.

Claro, isso significa que estamos repetindo algumas informações aqui, o que pode se tornar inconveniente se quisermos atualizar a localização para todas as saídas da mesma forma.
Mais tarde no curso, você aprenderá como configurar essas definições para múltiplas saídas de forma estruturada.

### 1.6. Execute o fluxo de trabalho com `-resume`

Vamos testar isso usando a flag `-resume`, já que já executamos a primeira etapa do fluxo de trabalho com sucesso.

```bash
nextflow run hello-workflow.nf -resume
```

??? success "Saída do comando"

    ```console
     N E X T F L O W   ~  version 25.10.2

    Launching `hello-workflow.nf` [high_cantor] DSL2 - revision: d746983511

    executor >  local (3)
    [ab/816321] process > sayHello (3)       [100%] 3 of 3, cached: 3 ✔
    [e0/ecf81b] process > convertToUpper (3) [100%] 3 of 3 ✔
    ```

Agora há uma linha extra na saída do console que corresponde ao novo processo que acabamos de adicionar.

Você encontrará as saídas no diretório `results/hello_workflow` conforme definido no bloco `output`.

??? abstract "Conteúdo do diretório"

    ```console
    results/hello_workflow/
    ├── Bonjour-output.txt
    ├── Hello-output.txt
    ├── Holà-output.txt
    ├── UPPER-Bonjour-output.txt
    ├── UPPER-Hello-output.txt
    └── UPPER-Holà-output.txt
    ```

Isso é conveniente! Mas ainda vale a pena dar uma olhada dentro do diretório de trabalho de uma das chamadas ao segundo processo.

??? abstract "Conteúdo do diretório"

    ```console
    work/e0/ecf81b4cacc648b9b994218d5b29d7/
    ├── Holà-output.txt -> /workspaces/training/hello-nextflow/work/ab/81632178cd37e9e815959278808819/Holà-output.txt
    └── UPPER-Holà-output.txt
    ```

Note que há dois arquivos `*-output`: a saída do primeiro processo assim como a saída do segundo.

A saída do primeiro processo está lá porque o Nextflow a **preparou** (staged) ali para ter tudo o que é necessário para execução dentro do mesmo subdiretório.

No entanto, é na verdade um link simbólico apontando para o arquivo original no subdiretório da primeira chamada de processo.
Por padrão, ao executar em uma única máquina como estamos fazendo aqui, o Nextflow usa links simbólicos em vez de cópias para preparar arquivos de entrada e intermediários.

Agora, antes de prosseguir, pense em como tudo o que fizemos foi conectar a saída de `sayHello` à entrada de `convertToUpper` e os dois processos puderam ser executados em série.
O Nextflow fez o trabalho pesado de gerenciar arquivos individuais de entrada e saída e passá-los entre os dois comandos para nós.

Esta é uma das razões pelas quais os canais do Nextflow são tão poderosos: eles cuidam do trabalho tedioso envolvido em conectar as etapas do fluxo de trabalho.

### Conclusão

Você sabe como encadear processos conectando a saída de uma etapa como entrada para a próxima etapa.

### O que vem a seguir?

Aprenda como coletar saídas de chamadas de processo em lote e alimentá-las em um único processo.

---

## 2. Adicione uma terceira etapa para coletar todas as saudações

Quando usamos um processo para aplicar uma transformação a cada um dos elementos em um canal, como estamos fazendo aqui para as múltiplas saudações, às vezes queremos coletar elementos do canal de saída desse processo e alimentá-los em outro processo que realiza algum tipo de análise ou soma.

Para demonstrar, adicionaremos uma nova etapa ao nosso pipeline que coleta todas as saudações em maiúsculas produzidas pelo processo `convertToUpper` e as escreve em um único arquivo.

<figure class="excalidraw">
--8<-- "docs/en/docs/hello_nextflow/img/hello-collect.svg"
</figure>

Sem estragar a surpresa, mas isso vai envolver um operador muito útil.

### 2.1. Defina o comando de coleta e teste-o no terminal

A etapa de coleta que queremos adicionar ao nosso fluxo de trabalho usará o comando `cat` para concatenar múltiplas saudações em maiúsculas em um único arquivo.

Vamos executar o comando sozinho no terminal para verificar que funciona conforme esperado, assim como fizemos anteriormente.

Execute o seguinte em seu terminal:

```bash
echo 'Hello' | tr '[a-z]' '[A-Z]' > UPPER-Hello-output.txt
echo 'Bonjour' | tr '[a-z]' '[A-Z]' > UPPER-Bonjour-output.txt
echo 'Holà' | tr '[a-z]' '[A-Z]' > UPPER-Holà-output.txt
cat UPPER-Hello-output.txt UPPER-Bonjour-output.txt UPPER-Holà-output.txt > COLLECTED-output.txt
```

A saída é um arquivo de texto chamado `COLLECTED-output.txt` que contém as versões em maiúsculas das saudações originais.

??? abstract "Conteúdo do arquivo"

    ```console title="COLLECTED-output.txt"
    HELLO
    BONJOUR
    HOLà
    ```

Esse é o resultado que queremos alcançar com nosso fluxo de trabalho.

### 2.2. Crie um novo processo para fazer a etapa de coleta

Vamos criar um novo processo e chamá-lo de `collectGreetings()`.
Podemos começar a escrevê-lo com base no que já vimos antes.

#### 2.2.1. Escreva as partes 'óbvias' do processo

Adicione a seguinte definição de processo ao script de fluxo de trabalho:

```groovy title="hello-workflow.nf" linenums="37"
/*
 * Coleta saudações em maiúsculas em um único arquivo de saída
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

Isso é o que podemos escrever com confiança com base no que você aprendeu até agora.
Mas isso não é funcional!
Deixa de fora a(s) definição(ões) de entrada e a primeira metade do comando do script porque precisamos descobrir como escrever isso.

#### 2.2.2. Defina as entradas para `collectGreetings()`

Precisamos coletar as saudações de todas as chamadas ao processo `convertToUpper()`.
O que sabemos que podemos obter da etapa anterior no fluxo de trabalho?

O canal de saída de `convertToUpper()` conterá os caminhos para os arquivos individuais contendo as saudações em maiúsculas.
Isso equivale a um slot de entrada; vamos chamá-lo de `input_files` por simplicidade.

No bloco de processo, faça a seguinte alteração de código:

=== "Depois"

    ```groovy title="hello-workflow.nf" linenums="42" hl_lines="2"
          input:
          path input_files
    ```

=== "Antes"

    ```groovy title="hello-workflow.nf" linenums="42" hl_lines="2"
          input:
          ???
    ```

Note que usamos o prefixo `path` mesmo esperando que isso contenha múltiplos arquivos.

#### 2.2.3. Componha o comando de concatenação

É aqui que as coisas podem ficar um pouco complicadas, porque precisamos ser capazes de lidar com um número arbitrário de arquivos de entrada.
Especificamente, não podemos escrever o comando antecipadamente, então precisamos dizer ao Nextflow como compô-lo em tempo de execução com base nas entradas que fluem para o processo.

Em outras palavras, se tivermos um canal de entrada contendo o elemento `[file1.txt, file2.txt, file3.txt]`, precisamos que o Nextflow transforme isso em `cat file1.txt file2.txt file3.txt`.

Felizmente, o Nextflow fica feliz em fazer isso por nós se simplesmente escrevermos `cat ${input_files}` no comando do script.

No bloco de processo, faça a seguinte alteração de código:

=== "Depois"

    ```groovy title="hello-workflow.nf" linenums="54" hl_lines="3"
        script:
        """
        cat ${input_files} > 'COLLECTED-output.txt'
        """
    ```

=== "Antes"

    ```groovy title="hello-workflow.nf" linenums="54"
        script:
        """
        cat ??? > 'COLLECTED-output.txt'
        """
    ```

Em teoria, isso deve lidar com qualquer número arbitrário de arquivos de entrada.

!!! tip

    Algumas ferramentas de linha de comando exigem fornecer um argumento (como `-input`) para cada arquivo de entrada.
    Nesse caso, teríamos que fazer um pouco de trabalho extra para compor o comando.
    Você pode ver um exemplo disso no curso de treinamento [Nextflow for Genomics](../../nf4_science/genomics/).

### 2.3. Adicione a etapa de coleta ao fluxo de trabalho

Agora devemos apenas precisar chamar o processo de coleta na saída da etapa de conversão para maiúsculas.
Esse também é um canal, chamado `convertToUpper.out`.

<figure class="excalidraw">
--8<-- "docs/en/docs/hello_nextflow/img/hello-collect-connector.svg"
</figure>

#### 2.3.1. Conecte as chamadas de processo

No bloco de fluxo de trabalho, faça a seguinte alteração de código:

=== "Depois"

    ```groovy title="hello-workflow.nf" linenums="75" hl_lines="4 5"
        // converte a saudação para maiúsculas
        convertToUpper(sayHello.out)

        // coleta todas as saudações em um arquivo
        collectGreetings(convertToUpper.out)
    }
    ```

=== "Antes"

    ```groovy title="hello-workflow.nf" linenums="75"
        // converte a saudação para maiúsculas
        convertToUpper(sayHello.out)
    }
    ```

Isso conecta a saída de `convertToUpper()` à entrada de `collectGreetings()`.

#### 2.3.2. Execute o fluxo de trabalho com `-resume`

Vamos tentar.

```bash
nextflow run hello-workflow.nf -resume
```

??? success "Saída do comando"

    ```console hl_lines="8"
    N E X T F L O W   ~  version 25.10.2

    Launching `hello-workflow.nf` [mad_gilbert] DSL2 - revision: 6acfd5e28d

    executor >  local (3)
    [79/33b2f0] sayHello (2)         | 3 of 3, cached: 3 ✔
    [99/79394f] convertToUpper (3)   | 3 of 3, cached: 3 ✔
    [47/50fe4a] collectGreetings (1) | 3 of 3 ✔
    ```

Ele é executado com sucesso, incluindo a terceira etapa.

No entanto, olhe o número de chamadas para `collectGreetings()` na última linha.
Estávamos esperando apenas uma, mas há três.

Agora dê uma olhada no conteúdo do arquivo de saída final.

??? abstract "Conteúdo do arquivo"

    ```console title="results/COLLECTED-output.txt"
    Holà
    ```

Oh não. A etapa de coleta foi executada individualmente em cada saudação, o que NÃO é o que queríamos.

<figure class="excalidraw">
--8<-- "docs/en/docs/hello_nextflow/img/hello-collect-no-operator.svg"
</figure>

Precisamos fazer algo para dizer ao Nextflow explicitamente que queremos que a terceira etapa seja executada em todos os elementos no canal de saída de `convertToUpper()`.

### 2.4. Use um operador para coletar as saudações em uma única entrada

Sim, mais uma vez a resposta para nosso problema é um operador.

Especificamente, vamos usar o operador apropriadamente chamado [`collect()`](https://nextflow.io/docs/latest/reference/operator.html#collect).

#### 2.4.1. Adicione o operador `collect()`

Desta vez vai parecer um pouco diferente porque não estamos adicionando o operador no contexto de uma fábrica de canal; estamos adicionando-o a um canal de saída.

Pegamos o `convertToUpper.out` e acrescentamos o operador `collect()`, o que nos dá `convertToUpper.out.collect()`.
Podemos conectar isso diretamente na chamada do processo `collectGreetings()`.

No bloco de fluxo de trabalho, faça a seguinte alteração de código:

=== "Depois"

    ```groovy title="hello-workflow.nf" linenums="73" hl_lines="2"
        // coleta todas as saudações em um arquivo
        collectGreetings(convertToUpper.out.collect())
    }
    ```

=== "Antes"

    ```groovy title="hello-workflow.nf" linenums="73" hl_lines="2"
        // coleta todas as saudações em um arquivo
        collectGreetings(convertToUpper.out)
    }
    ```

#### 2.4.2. Adicione algumas instruções `view()`

Vamos também incluir algumas instruções `view()` para visualizar os estados antes e depois do conteúdo do canal.

=== "Depois"

    ```groovy title="hello-workflow.nf" linenums="73" hl_lines="4-6"
        // coleta todas as saudações em um arquivo
        collectGreetings(convertToUpper.out.collect())

        // instruções view opcionais
        convertToUpper.out.view { contents -> "Antes do collect: $contents" }
        convertToUpper.out.collect().view { contents -> "Depois do collect: $contents" }
    }
    ```

=== "Antes"

    ```groovy title="hello-workflow.nf" linenums="73"
        // coleta todas as saudações em um arquivo
        collectGreetings(convertToUpper.out.collect())
    }
    ```

As instruções `view()` podem ir em qualquer lugar que você quiser; as colocamos logo após a chamada para legibilidade.

#### 2.4.3. Execute o fluxo de trabalho novamente com `-resume`

Vamos tentar:

```bash
nextflow run hello-workflow.nf -resume
```

??? success "Saída do comando"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `hello-workflow.nf` [soggy_franklin] DSL2 - revision: bc8e1b2726

    [d6/cdf466] sayHello (1)       | 3 of 3, cached: 3 ✔
    [99/79394f] convertToUpper (2) | 3 of 3, cached: 3 ✔
    [1e/83586c] collectGreetings   | 1 of 1 ✔
    Antes do collect: /workspaces/training/hello-nextflow/work/b3/d52708edba8b864024589285cb3445/UPPER-Bonjour-output.txt
    Antes do collect: /workspaces/training/hello-nextflow/work/99/79394f549e3040dfc2440f69ede1fc/UPPER-Hello-output.txt
    Antes do collect: /workspaces/training/hello-nextflow/work/aa/56bfe7cf00239dc5badc1d04b60ac4/UPPER-Holà-output.txt
    Depois do collect: [/workspaces/training/hello-nextflow/work/b3/d52708edba8b864024589285cb3445/UPPER-Bonjour-output.txt, /workspaces/training/hello-nextflow/work/99/79394f549e3040dfc2440f69ede1fc/UPPER-Hello-output.txt, /workspaces/training/hello-nextflow/work/aa/56bfe7cf00239dc5badc1d04b60ac4/UPPER-Holà-output.txt]
    ```

Ele é executado com sucesso, embora a saída de log possa parecer um pouco mais bagunçada do que isso (nós a limpamos para legibilidade).

Desta vez a terceira etapa foi chamada apenas uma vez!
Olhando para a saída das instruções `view()`, vemos o seguinte:

- Três instruções `Antes do collect:`, uma para cada saudação: nesse ponto os caminhos dos arquivos são itens individuais no canal.
- Uma única instrução `Depois do collect:`: os três caminhos de arquivos agora estão empacotados em um único elemento.

Podemos resumir isso com o seguinte diagrama:

<figure class="excalidraw">
--8<-- "docs/en/docs/hello_nextflow/img/hello-collect-WITH-operator.svg"
</figure>

Finalmente, você pode dar uma olhada no conteúdo do arquivo de saída para se convencer de que tudo funcionou corretamente.

??? abstract "Conteúdo do arquivo"

    ```console title="results/COLLECTED-output.txt"
    BONJOUR
    HELLO
    HOLà
    ```

Desta vez temos todas as três saudações no arquivo de saída final. Sucesso!

!!! note

    Se você executar isso várias vezes sem `-resume`, verá que a ordem das saudações muda de uma execução para outra.
    Isso mostra que a ordem na qual os elementos fluem através das chamadas de processo não é garantida ser consistente.

#### 2.4.4. Remova as instruções `view()` para legibilidade

Antes de passar para a próxima seção, recomendamos que você delete as instruções `view()` para evitar poluir a saída do console.

=== "Depois"

    ```groovy title="hello-workflow.nf" linenums="73"
        // coleta todas as saudações em um arquivo
        collectGreetings(convertToUpper.out.collect())
    ```

=== "Antes"

    ```groovy title="hello-workflow.nf" linenums="73" hl_lines="4-6"
        // coleta todas as saudações em um arquivo
        collectGreetings(convertToUpper.out.collect())

        // instruções view opcionais
        convertToUpper.out.view { contents -> "Antes do collect: $contents" }
        convertToUpper.out.collect().view { contents -> "Depois do collect: $contents" }
    ```

Esta é basicamente a operação inversa do ponto 2.4.2.

### Conclusão

Você sabe como coletar saídas de um lote de chamadas de processo e alimentá-las em uma etapa de análise conjunta ou soma.

Para recapitular, isto é o que você construiu até agora:

<figure class="excalidraw">
--8<-- "docs/en/docs/hello_nextflow/img/hello-collect.svg"
</figure>

### O que vem a seguir?

Aprenda como passar mais de uma entrada para um processo.

---

## 3. Passe parâmetros adicionais para um processo

Queremos ser capazes de nomear o arquivo de saída final com algo específico para processar lotes subsequentes de saudações sem sobrescrever os resultados finais.

Para isso, faremos os seguintes refinamentos no fluxo de trabalho:

- Modificar o processo coletor para aceitar um nome definido pelo usuário para o arquivo de saída (`batch_name`)
- Adicionar um parâmetro de linha de comando ao fluxo de trabalho (`--batch`) e passá-lo ao processo coletor

<figure class="excalidraw">
--8<-- "docs/en/docs/hello_nextflow/img/hello-collect-batch.svg"
</figure>

### 3.1. Modifique o processo coletor

Vamos precisar declarar a entrada adicional e integrá-la ao nome do arquivo de saída.

#### 3.1.1. Declare a entrada adicional

Boas notícias: podemos declarar quantas variáveis de entrada quisermos na definição do processo.
Vamos chamar esta de `batch_name`.

No bloco de processo, faça a seguinte alteração de código:

=== "Depois"

    ```groovy title="hello-workflow.nf" linenums="42" hl_lines="3"
        input:
        path input_files
        val batch_name
    ```

=== "Antes"

    ```groovy title="hello-workflow.nf" linenums="42"
        input:
        path input_files
    ```

Você pode configurar seus processos para esperar quantas entradas quiser.
Agora, todas estas estão configuradas para serem entradas obrigatórias; você _deve_ fornecer um valor para o fluxo de trabalho funcionar.

Você aprenderá como gerenciar entradas obrigatórias versus opcionais mais tarde em sua jornada Nextflow.

#### 3.1.2. Use a variável `batch_name` no nome do arquivo de saída

Podemos inserir a variável no nome do arquivo de saída da mesma forma que compusemos nomes de arquivos dinâmicos antes.

No bloco de processo, faça a seguinte alteração de código:

=== "Depois"

    ```groovy title="hello-workflow.nf" linenums="46" hl_lines="2 6"
        output:
        path "COLLECTED-${batch_name}-output.txt"

        script:
        """
        cat ${input_files} > 'COLLECTED-${batch_name}-output.txt'
        """
    ```

=== "Antes"

    ```groovy title="hello-workflow.nf" linenums="46" hl_lines="2 6"
        output:
        path "COLLECTED-output.txt"

        script:
        """
        cat ${input_files} > 'COLLECTED-output.txt'
        """
    ```

Isso configura o processo para usar o valor `batch_name` para gerar um nome de arquivo específico para a saída final do fluxo de trabalho.

### 3.2. Adicione um parâmetro de linha de comando `batch`

Agora precisamos de uma forma de fornecer o valor para `batch_name` e alimentá-lo à chamada do processo.

#### 3.2.1. Use `params` para configurar o parâmetro

Você já sabe como usar o sistema `params` para declarar parâmetros CLI.
Vamos usar isso para declarar um parâmetro `batch` (com um valor padrão porque somos preguiçosos).

Na seção de parâmetros do pipeline, faça as seguintes alterações de código:

=== "Depois"

    ```groovy title="hello-workflow.nf" linenums="55" hl_lines="6"
    /*
     * Parâmetros do pipeline
     */
    params {
        input: Path = 'data/greetings.csv'
        batch: String = 'batch'
    }
    ```

=== "Antes"

    ```groovy title="hello-workflow.nf" linenums="55"
    /*
     * Parâmetros do pipeline
     */
    params {
        input: Path = 'data/greetings.csv'
    }
    ```

Assim como demonstramos para `--input`, você pode sobrescrever esse valor padrão especificando um valor com `--batch` na linha de comando.

#### 3.2.2. Passe o parâmetro `batch` para o processo

Para fornecer o valor do parâmetro ao processo, precisamos adicioná-lo na chamada do processo.

No bloco de fluxo de trabalho, faça a seguinte alteração de código:

=== "Depois"

    ```groovy title="hello-workflow.nf" linenums="74" hl_lines="2"
        // coleta todas as saudações em um arquivo
        collectGreetings(convertToUpper.out.collect(), params.batch)
    ```

=== "Antes"

    ```groovy title="hello-workflow.nf" linenums="74" hl_lines="2"
        // coleta todas as saudações em um arquivo
        collectGreetings(convertToUpper.out.collect())
    ```

Você vê que para fornecer múltiplas entradas a um processo, você simplesmente as lista nos parênteses da chamada, separadas por vírgulas.

!!! warning

    Você DEVE fornecer as entradas ao processo na MESMA ORDEM EXATA em que estão listadas no bloco de definição de entrada do processo.

### 3.3. Execute o fluxo de trabalho

Vamos tentar executar isso com um nome de lote na linha de comando.

```bash
nextflow run hello-workflow.nf -resume --batch trio
```

??? success "Saída do comando"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `hello-workflow.nf` [confident_rutherford] DSL2 - revision: bc58af409c

    executor >  local (1)
    [79/33b2f0] sayHello (2)       | 3 of 3, cached: 3 ✔
    [99/79394f] convertToUpper (2) | 3 of 3, cached: 3 ✔
    [b5/f19efe] collectGreetings   | 1 of 1 ✔
    ```

Ele é executado com sucesso e produz a saída desejada:

??? abstract "Conteúdo do arquivo"

    ```console title="results/COLLECTED-trio-output.txt"
    HELLO
    BONJOUR
    HOLà
    ```

Agora, contanto que especifiquemos o parâmetro apropriadamente, execuções subsequentes em outros lotes de entradas não destruirão os resultados anteriores.

### Conclusão

Você sabe como passar mais de uma entrada para um processo.

### O que vem a seguir?

Aprenda como emitir múltiplas saídas e manuseá-las convenientemente.

---

## 4. Adicione uma saída à etapa coletora

Até agora estivemos usando processos que produziam apenas uma saída cada.
Conseguimos acessar suas respectivas saídas muito convenientemente usando a sintaxe `<process>.out`, que usamos tanto no contexto de passar uma saída para o próximo processo (por exemplo, `convertToUpper(sayHello.out)`) quanto no contexto da seção `publish:` (por exemplo, `first_output = sayHello.out`).

O que acontece quando um processo produz mais de uma?
Como lidamos com as múltiplas saídas?
Podemos selecionar e usar uma saída específica?

Todas excelentes perguntas, e a resposta curta é sim, podemos!

Múltiplas saídas serão empacotadas em canais separados.
Podemos escolher dar nomes a esses canais de saída, o que torna fácil referenciá-los individualmente mais tarde, ou podemos referenciá-los por índice.

Para fins de demonstração, digamos que queremos contar o número de saudações que estão sendo coletadas para um determinado lote de entradas e relatá-lo em um arquivo.

### 4.1. Modifique o processo para contar e gerar o número de saudações

Isso exigirá duas mudanças-chave na definição do processo: precisamos de uma forma de contar as saudações e escrever um arquivo de relatório, então precisamos adicionar esse arquivo de relatório ao bloco `output` do processo.

#### 4.1.1. Conte o número de saudações coletadas

Convenientemente, o Nextflow nos permite adicionar código arbitrário no bloco `script:` da definição do processo, o que é muito útil para fazer coisas como essa.

Isso significa que podemos usar a função integrada `size()` do Nextflow para obter o número de arquivos no array `input_files`, e escrever o resultado em arquivo com um comando `echo`.

No bloco de processo `collectGreetings`, faça as seguintes alterações de código:

=== "Depois"

    ```groovy title="hello-workflow.nf" linenums="55" hl_lines="2 5"
        script:
        count_greetings = input_files.size()
        """
        cat ${input_files} > 'COLLECTED-${batch_name}-output.txt'
        echo 'Havia ${count_greetings} saudações neste lote.' > '${batch_name}-report.txt'
        """
    ```

=== "Antes"

    ```groovy title="hello-workflow.nf" linenums="55"
        script:
        """
        cat ${input_files} > 'COLLECTED-${batch_name}-output.txt'
        """
    ```

A variável `count_greetings` será computada em tempo de execução.

#### 4.1.2. Emita o arquivo de relatório e nomeie as saídas

Em princípio, tudo o que precisamos fazer é adicionar o arquivo de relatório ao bloco `output:`.

No entanto, enquanto estamos fazendo isso, também vamos adicionar algumas tags `emit:` às nossas declarações de saída. Estas nos permitirão selecionar as saídas por nome em vez de ter que usar índices posicionais.

No bloco de processo, faça a seguinte alteração de código:

=== "Depois"

    ```groovy title="hello-workflow.nf" linenums="46" hl_lines="2 3"
        output:
        path "COLLECTED-${batch_name}-output.txt", emit: outfile
        path "${batch_name}-report.txt", emit: report
    ```

=== "Antes"

    ```groovy title="hello-workflow.nf" linenums="46"
        output:
        path "COLLECTED-${batch_name}-output.txt"
    ```

As tags `emit:` são opcionais, e poderíamos ter adicionado uma tag a apenas uma das saídas.
Mas como dizem, por que não ambos?

!!! tip

    Se você não nomear as saídas de um processo usando `emit:`, ainda pode acessá-las individualmente usando seu respectivo índice (baseado em zero).
    Por exemplo, você usaria `<process>.out[0]` para obter a primeira saída, `<process>.out[1]` para obter a segunda saída, e assim por diante.

    Preferimos nomear saídas porque caso contrário, é muito fácil pegar o índice errado por erro, especialmente quando o processo produz muitas saídas.

### 4.2. Atualize as saídas do fluxo de trabalho

Agora que temos duas saídas saindo do processo `collectGreetings`, a saída `collectGreetings.out` contém dois canais:

- `collectGreetings.out.outfile` contém o arquivo de saída final
- `collectGreetings.out.report` contém o arquivo de relatório

<figure class="excalidraw">
--8<-- "docs/en/docs/hello_nextflow/img/hello-collect-report.svg"
</figure>

Precisamos atualizar as saídas do fluxo de trabalho adequadamente.

#### 4.2.1. Atualize a seção `publish:`

No bloco `workflow`, faça a seguinte alteração de código:

=== "Depois"

    ```groovy title="hello-workflow.nf" linenums="80" hl_lines="4 5"
        publish:
        first_output = sayHello.out
        uppercased = convertToUpper.out
        collected = collectGreetings.out.outfile
        batch_report = collectGreetings.out.report
    ```

=== "Antes"

    ```groovy title="hello-workflow.nf" linenums="80" hl_lines="4"
        publish:
        first_output = sayHello.out
        uppercased = convertToUpper.out
        collected = collectGreetings.out
    ```

Como você pode ver, referir-se a saídas específicas de processos agora é trivial.
Quando formos adicionar mais um passo ao nosso pipeline na Parte 5 (Contêineres), poderemos facilmente nos referir a `collectGreetings.out.outfile` e passá-lo ao novo processo (spoiler: o novo processo se chama `cowpy`).

Mas por enquanto, vamos terminar de atualizar as saídas no nível do fluxo de trabalho.

#### 4.2.2. Atualize o bloco `output`

No bloco `output`, faça a seguinte alteração de código:

=== "Depois"

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

=== "Antes"

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

Não precisamos atualizar a definição de saída `collected` já que esse nome não mudou.
Só precisamos adicionar a nova saída.

### 4.3. Execute o fluxo de trabalho

Vamos tentar executar isso com o lote atual de saudações.

```bash
nextflow run hello-workflow.nf -resume --batch trio
```

??? success "Saída do comando"

    ```console
     N E X T F L O W   ~  version 25.10.2

    Launching `hello-workflow.nf` [ecstatic_wilson] DSL2 - revision: c80285f8c8

    executor >  local (1)
    [c5/4c6ca9] sayHello (3)       [100%] 3 of 3, cached: 3 ✔
    [0e/6cbc59] convertToUpper (3) [100%] 3 of 3, cached: 3 ✔
    [02/61ead2] collectGreetings   [100%] 1 of 1 ✔
    ```

Se você olhar no diretório `results/hello_workflow/`, encontrará o novo arquivo de relatório, `trio-report.txt`.
Abra-o para verificar que o fluxo de trabalho relatou corretamente a contagem de saudações que foram processadas.

??? abstract "Conteúdo do arquivo"

    ```txt title="trio-report.txt"
    Havia 3 saudações neste lote.
    ```

<figure class="excalidraw">
--8<-- "docs/en/docs/hello_nextflow/img/hello-collect-4-way.svg"
</figure>

Sinta-se à vontade para adicionar mais saudações ao CSV e testar o que acontece.

### Conclusão

Você sabe como fazer um processo emitir múltiplas saídas nomeadas e como manuseá-las apropriadamente no nível do fluxo de trabalho.

De forma mais geral, você entende os princípios-chave envolvidos em conectar processos de formas comuns.

### O que vem a seguir?

Faça uma pausa extra longa, você a merece.

Quando estiver pronto, passe para [**Parte 4: Hello Modules**](./04_hello_modules.md) para aprender como modularizar seu código para melhor manutenibilidade e eficiência de código.

---

## Quiz

<quiz>
Como você acessa a saída de um processo no bloco de fluxo de trabalho?
- [ ] `process.output`
- [ ] `output.processName`
- [x] `processName.out`
- [ ] `get(processName)`

Saiba mais: [1.4. Passe a saída do primeiro processo para o segundo processo](#14-passe-a-saída-do-primeiro-processo-para-o-segundo-processo)
</quiz>

<quiz>
O que determina a ordem de execução de processos no Nextflow?
- [ ] A ordem em que os processos são escritos no bloco de fluxo de trabalho
- [ ] Ordem alfabética pelo nome do processo
- [x] Dependências de dados entre processos
- [ ] Ordem aleatória para execução paralela

Saiba mais: [1.4. Passe a saída do primeiro processo para o segundo processo](#14-passe-a-saída-do-primeiro-processo-para-o-segundo-processo)
</quiz>

<quiz>
Qual operador deve substituir `???` para reunir todas as saídas em uma única lista para o processo downstream?

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

Saiba mais: [2.4. Use um operador para coletar as saudações em uma única entrada](#24-use-um-operador-para-coletar-as-saudações-em-uma-única-entrada)
</quiz>

<quiz>
Quando você deve usar o operador `collect()`?
- [ ] Quando você quer processar itens em paralelo
- [ ] Quando você precisa filtrar conteúdo do canal
- [x] Quando um processo downstream precisa de todos os itens de um processo upstream
- [ ] Quando você quer dividir dados entre múltiplos processos

Saiba mais: [2.4. Use um operador para coletar as saudações em uma única entrada](#24-use-um-operador-para-coletar-as-saudações-em-uma-única-entrada)
</quiz>

<quiz>
Como você acessa uma saída nomeada de um processo?
- [ ] `processName.outputName`
- [ ] `processName.get(outputName)`
- [x] `processName.out.outputName`
- [ ] `output.processName.outputName`

Saiba mais: [4.1.2. Emita o arquivo de relatório e nomeie as saídas](#412-emita-o-arquivo-de-relatório-e-nomeie-as-saídas)
</quiz>

<quiz>
Qual é a sintaxe correta para nomear uma saída em um processo?
- [ ] `name: outputName`
- [ ] `output: outputName`
- [x] `emit: outputName`
- [ ] `label: outputName`

Saiba mais: [4.1.2. Emita o arquivo de relatório e nomeie as saídas](#412-emita-o-arquivo-de-relatório-e-nomeie-as-saídas)
</quiz>

<quiz>
Ao fornecer múltiplas entradas a um processo, o que deve ser verdadeiro?
- [ ] Todas as entradas devem ser do mesmo tipo
- [ ] As entradas devem ser fornecidas em ordem alfabética
- [x] A ordem das entradas deve corresponder à ordem definida no bloco de entrada
- [ ] Apenas duas entradas podem ser fornecidas por vez

Saiba mais: [3. Passe parâmetros adicionais para um processo](#3-passe-parâmetros-adicionais-para-um-processo)
</quiz>
