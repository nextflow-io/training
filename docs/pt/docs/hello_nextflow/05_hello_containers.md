# Parte 5: Hello Containers

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Tradução assistida por IA - [saiba mais e sugira melhorias](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

<div class="video-wrapper">
  <iframe width="560" height="315" src="https://www.youtube.com/embed/Xqr--bKEN9U?si=y8lAedhEHWaTV4zd&amp;list=PLPZ8WHdZGxmWKozQuzr27jyMGqp9kElVK&amp;cc_load_policy=1" title="YouTube video player" frameborder="0" allow="accelerometer; autoplay; clipboard-write; encrypted-media; gyroscope; picture-in-picture; web-share" referrerpolicy="strict-origin-when-cross-origin" allowfullscreen></iframe>
</div>

/// caption
:fontawesome-brands-youtube:{ .youtube } Veja [a playlist completa](https://youtube.com/playlist?list=PLPZ8WHdZGxmWKozQuzr27jyMGqp9kElVK&si=eF7cLR62goy-lc6n) no canal do Nextflow no YouTube.

:green_book: A transcrição do vídeo está disponível [aqui](./transcripts/05_hello_containers.md).
///

Nas Partes 1-4 deste curso de treinamento, você aprendeu a usar os blocos de construção básicos do Nextflow para montar um fluxo de trabalho simples capaz de processar texto, paralelizar a execução se houvesse múltiplas entradas e coletar os resultados para processamento adicional.

No entanto, você estava limitado às ferramentas básicas do UNIX disponíveis no seu ambiente.
Tarefas do mundo real frequentemente requerem várias ferramentas e pacotes não incluídos por padrão.
Normalmente, você precisaria instalar essas ferramentas, gerenciar suas dependências e resolver quaisquer conflitos.

Tudo isso é muito tedioso e irritante, então vamos mostrar como usar **contêineres** para resolver esse problema de forma muito mais conveniente.

Um **contêiner** é uma unidade de software leve, autônoma e executável criada a partir de uma **imagem** de contêiner que inclui tudo o que é necessário para executar uma aplicação, incluindo código, bibliotecas do sistema e configurações.
Como você pode imaginar, isso será muito útil para tornar seus pipelines mais reproduzíveis.

Note que ensinaremos isso usando [Docker](https://www.docker.com/get-started/), mas tenha em mente que o Nextflow suporta [várias outras tecnologias de contêiner](https://nextflow.io/docs/latest/container.html) também.

??? info "Como começar a partir desta seção"

    Esta seção do curso assume que você completou as Partes 1-4 do curso [Hello Nextflow](./index.md) e tem um pipeline completo e funcional.

    Se você está começando o curso a partir deste ponto, precisará copiar o diretório `modules` das soluções:

    ```bash
    cp -r solutions/4-hello-modules/modules .
    ```

---

## 0. Aquecimento: Execute `hello-containers.nf`

Vamos usar o script de fluxo de trabalho `hello-containers.nf` como ponto de partida.
Ele é equivalente ao script produzido ao trabalhar na Parte 4 deste curso de treinamento, exceto que mudamos os destinos de saída:

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

Só para garantir que tudo está funcionando, execute o script uma vez antes de fazer qualquer alteração:

```bash
nextflow run hello-containers.nf
```

??? success "Saída do comando"

    ```console
     N E X T F L O W   ~  version 25.10.2

    Launching `hello-containers.nf` [nice_escher] DSL2 - revision: d5dfdc9872

    executor > local (7)
    [5a/ec1fa1] sayHello (2) [100%] 3 of 3 ✔
    [30/32b5b8] convertToUpper (3) [100%] 3 of 3 ✔
    [d3/be01bc] collectGreetings [100%] 1 of 1 ✔

    ```

Como anteriormente, você encontrará os arquivos de saída no diretório especificado no bloco `output` (`results/hello_containers/`).

??? abstract "Conteúdo do diretório"

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

Se isso funcionou para você, você está pronto para aprender como usar contêineres.

---

## 1. Use um contêiner 'manualmente'

O que queremos fazer é adicionar uma etapa ao nosso fluxo de trabalho que usará um contêiner para execução.

No entanto, primeiro vamos revisar alguns conceitos e operações básicas para solidificar seu entendimento do que são contêineres antes de começarmos a usá-los no Nextflow.

### 1.1. Baixe a imagem do contêiner

Para usar um contêiner, você geralmente baixa ou _puxa_ uma imagem de contêiner de um registro de contêineres, e então executa a imagem do contêiner para criar uma instância de contêiner.

A sintaxe geral é a seguinte:

```bash title="Syntax"
docker pull '<container>'
```

A parte `docker pull` é a instrução para o sistema de contêiner puxar uma imagem de contêiner de um repositório.

A parte `'<container>'` é o endereço URI da imagem do contêiner.

Como exemplo, vamos puxar uma imagem de contêiner que contém [cowpy](https://github.com/jeffbuttars/cowpy), uma implementação em Python de uma ferramenta chamada `cowsay` que gera arte ASCII para exibir entradas de texto arbitrárias de forma divertida.

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

Existem vários repositórios onde você pode encontrar contêineres publicados.
Usamos o serviço [Seqera Containers](https://seqera.io/containers/) para gerar esta imagem de contêiner Docker a partir do pacote Conda `cowpy`: `'community.wave.seqera.io/library/cowpy:1.1.5--3db457ae1977a273'`.

Execute o comando pull completo:

```bash
docker pull 'community.wave.seqera.io/library/cowpy:1.1.5--3db457ae1977a273'
```

??? success "Saída do comando"

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

Se você nunca baixou a imagem antes, isso pode levar um minuto para completar.
Uma vez concluído, você tem uma cópia local da imagem do contêiner.

### 1.2. Use o contêiner para executar `cowpy` como um comando único

Uma maneira muito comum de as pessoas usarem contêineres é executá-los diretamente, _ou seja_, de forma não interativa.
Isso é ótimo para executar comandos únicos.

A sintaxe geral é a seguinte:

```bash title="Syntax"
docker run --rm '<container>' [tool command]
```

A parte `docker run --rm '<container>'` é a instrução para o sistema de contêiner criar uma instância de contêiner a partir de uma imagem de contêiner e executar um comando nela.
A flag `--rm` diz ao sistema para desligar a instância do contêiner após o comando ter sido completado.

A sintaxe `[comando da ferramenta]` depende da ferramenta que você está usando e de como o contêiner está configurado.
Vamos começar apenas com `cowpy`.

Totalmente montado, o comando de execução do contêiner fica assim; vá em frente e execute-o.

```bash
docker run --rm 'community.wave.seqera.io/library/cowpy:1.1.5--3db457ae1977a273' cowpy
```

??? success "Saída do comando"

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

O sistema criou o contêiner, executou o comando `cowpy` com seus parâmetros, enviou a saída para o console e finalmente, desligou a instância do contêiner.

### 1.3. Use o contêiner para executar `cowpy` interativamente

Você também pode executar um contêiner interativamente, o que lhe dá um prompt de shell dentro do contêiner e permite que você brinque com o comando.

#### 1.3.1. Inicie o contêiner

Para executar interativamente, apenas adicionamos `-it` ao comando `docker run`.
Opcionalmente, podemos especificar o shell que queremos usar dentro do contêiner adicionando _por exemplo_ `/bin/bash` ao comando.

```bash
docker run --rm -it 'community.wave.seqera.io/library/cowpy:1.1.5--3db457ae1977a273' /bin/bash
```

Note que seu prompt muda para algo como `(base) root@b645838b3314:/tmp#`, o que indica que você está agora dentro do contêiner.

Você pode verificar isso executando `ls /` para listar o conteúdo do diretório a partir da raiz do sistema de arquivos:

```bash
ls /
```

??? abstract "Saída do comando"

    ```console
    bin  boot  dev  etc  home  lib  lib64  media  mnt  opt  proc  root  run  sbin  srv  sys  tmp  usr  var
    ```

Usamos `ls` aqui em vez de `tree` porque o utilitário `tree` não está disponível neste contêiner.
Você pode ver que o sistema de arquivos dentro do contêiner é diferente do sistema de arquivos no seu sistema host.

Uma limitação do que acabamos de fazer é que o contêiner está completamente isolado do sistema host por padrão.
Isso significa que o contêiner não pode acessar nenhum arquivo no sistema host a menos que você explicitamente permita que ele faça isso.

Vamos mostrar como fazer isso em um minuto.

#### 1.3.2. Execute o(s) comando(s) da ferramenta desejada

Agora que você está dentro do contêiner, pode executar o comando `cowpy` diretamente e dar alguns parâmetros a ele.
Por exemplo, a documentação da ferramenta diz que podemos mudar o personagem ('cowacter') com `-c`.

```bash
cowpy "Hello Containers" -c tux
```

??? success "Saída do comando"

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

Agora a saída mostra o pinguim do Linux, Tux, em vez da vaca padrão, porque especificamos o parâmetro `-c tux`.

Como você está dentro do contêiner, pode executar o comando `cowpy` quantas vezes quiser, variando os parâmetros de entrada, sem ter que se preocupar com comandos Docker.

!!! Tip "Dica"

    Use a flag '-c' para escolher um personagem diferente, incluindo:
    `beavis`, `cheese`, `daemon`, `dragonandcow`, `ghostbusters`, `kitty`, `moose`, `milk`, `stegosaurus`, `turkey`, `turtle`, `tux`

Isso é legal. O que seria ainda mais legal é se pudéssemos alimentar nosso `greetings.csv` como entrada nisso.
Mas como não temos acesso ao sistema de arquivos, não podemos.

Vamos consertar isso.

#### 1.3.3. Saia do contêiner

Para sair do contêiner, você pode digitar `exit` no prompt ou usar o atalho de teclado ++ctrl+d++.

```bash
exit
```

Seu prompt agora deve estar de volta ao que era antes de você iniciar o contêiner.

#### 1.3.4. Monte dados no contêiner

Como observado anteriormente, o contêiner está isolado do sistema host por padrão.

Para permitir que o contêiner acesse o sistema de arquivos do host, você pode **montar** um **volume** do sistema host no contêiner usando a seguinte sintaxe:

```bash title="Syntax"
-v <outside_path>:<inside_path>
```

No nosso caso, `<caminho_externo>` será o diretório de trabalho atual, então podemos apenas usar um ponto (`.`), e `<caminho_interno>` é apenas um alias que inventamos; vamos chamá-lo de `/my_project` (o caminho interno deve ser absoluto).

Para montar um volume, substituímos os caminhos e adicionamos o argumento de montagem de volume ao comando docker run da seguinte forma:

```bash
docker run --rm -it -v .:/my_project 'community.wave.seqera.io/library/cowpy:1.1.5--3db457ae1977a273' /bin/bash
```

Isso monta o diretório de trabalho atual como um volume que será acessível em `/my_project` dentro do contêiner.

Você pode verificar que funciona listando o conteúdo de `/my_project`:

```bash
ls /my_project
```

??? success "Saída do comando"

    ```console
    data               hello-config.nf      hello-modules.nf   hello-world.nf  nextflow.config  solutions         work
    hello-channels.nf  hello-containers.nf  hello-workflow.nf  modules         results          test-params.json
    ```

Agora você pode ver o conteúdo do diretório de trabalho de dentro do contêiner, incluindo o arquivo `greetings.csv` em `data/`.

Isso efetivamente estabeleceu um túnel através da parede do contêiner que você pode usar para acessar essa parte do seu sistema de arquivos.

#### 1.3.5. Use os dados montados

Agora que montamos o diretório de trabalho no contêiner, podemos usar o comando `cowpy` para exibir o conteúdo do arquivo `greetings.csv`.

Para fazer isso, usaremos `cat /my_project/data/greetings.csv | ` para canalizar o conteúdo do arquivo CSV para o comando `cowpy`.

```bash
cat /my_project/data/greetings.csv | cowpy -c turkey
```

??? success "Saída do comando"

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

Isso produz a arte ASCII desejada de um peru recitando nossas saudações de exemplo!
Exceto que aqui o peru está repetindo as linhas completas em vez de apenas as saudações.
Já sabemos que nosso fluxo de trabalho Nextflow fará um trabalho melhor!

Sinta-se à vontade para brincar com este comando.
Quando terminar, saia do contêiner como anteriormente:

```bash
exit
```

Você se encontrará de volta ao seu shell normal.

### Conclusão

Você sabe como puxar um contêiner e executá-lo como um comando único ou interativamente. Você também sabe como tornar seus dados acessíveis de dentro do seu contêiner, o que permite que você experimente qualquer ferramenta em que esteja interessado com dados reais sem ter que instalar nenhum software no seu sistema.

### O que vem a seguir?

Aprenda como usar contêineres para a execução de processos Nextflow.

---

## 2. Use contêineres no Nextflow

O Nextflow tem suporte integrado para executar processos dentro de contêineres para permitir que você execute ferramentas que não tem instaladas no seu ambiente de computação.
Isso significa que você pode usar qualquer imagem de contêiner que desejar para executar seus processos, e o Nextflow cuidará de puxar a imagem, montar os dados e executar o processo dentro dela.

Para demonstrar isso, vamos adicionar uma etapa `cowpy` ao pipeline que estivemos desenvolvendo, após a etapa `collectGreetings`.

<figure class="excalidraw">
--8<-- "docs/en/docs/hello_nextflow/img/hello-pipeline-cowpy.svg"
</figure>

### 2.1. Escreva um módulo `cowpy`

Primeiro, vamos criar o módulo de processo `cowpy`.

#### 2.1.1. Crie um arquivo stub para o novo módulo

Crie um arquivo vazio para o módulo chamado `cowpy.nf`.

```bash
touch modules/cowpy.nf
```

Isso nos dá um lugar para colocar o código do processo.

#### 2.1.2. Copie o código do processo `cowpy` no arquivo do módulo

Podemos modelar nosso processo `cowpy` nos outros processos que escrevemos anteriormente.

```groovy title="modules/cowpy.nf" linenums="1"
#!/usr/bin/env nextflow

// Gera arte ASCII com cowpy
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

O processo espera um `input_file` contendo as saudações, bem como um valor `character`.

A saída será um novo arquivo de texto contendo a arte ASCII gerada pela ferramenta `cowpy`.

### 2.2. Adicione cowpy ao fluxo de trabalho

Agora precisamos importar o módulo e chamar o processo.

#### 2.2.1. Importe o processo `cowpy` em `hello-containers.nf`

Insira a declaração de importação acima do bloco workflow e preencha-a apropriadamente.

=== "Depois"

    ```groovy title="hello-containers.nf" linenums="3" hl_lines="5"
    // Inclui módulos
    include { sayHello } from './modules/sayHello.nf'
    include { convertToUpper } from './modules/convertToUpper.nf'
    include { collectGreetings } from './modules/collectGreetings.nf'
    include { cowpy } from './modules/cowpy.nf'
    ```

=== "Antes"

    ```groovy title="hello-containers.nf" linenums="3"
    // Inclui módulos
    include { sayHello } from './modules/sayHello.nf'
    include { convertToUpper } from './modules/convertToUpper.nf'
    include { collectGreetings } from './modules/collectGreetings.nf'
    ```

Agora o módulo `cowpy` está disponível para uso no fluxo de trabalho.

#### 2.2.2. Adicione uma chamada ao processo `cowpy` no fluxo de trabalho

Vamos conectar o processo `cowpy()` à saída do processo `collectGreetings()`, que como você deve se lembrar produz duas saídas:

- `collectGreetings.out.outfile` contém o arquivo de saída <--_o que queremos_
- `collectGreetings.out.report` contém o arquivo de relatório com a contagem de saudações por lote

No bloco workflow, faça a seguinte alteração de código:

=== "Depois"

    ```groovy title="hello-containers.nf" linenums="19" hl_lines="12-13"
        main:
        // cria um canal para entradas de um arquivo CSV
        greeting_ch = channel.fromPath(params.input)
                            .splitCsv()
                            .map { line -> line[0] }
        // emite uma saudação
        sayHello(greeting_ch)
        // converte a saudação para maiúsculas
        convertToUpper(sayHello.out)
        // coleta todas as saudações em um arquivo
        collectGreetings(convertToUpper.out.collect(), params.batch)
        // gera arte ASCII das saudações com cowpy
        cowpy(collectGreetings.out.outfile, params.character)
    ```

=== "Antes"

    ```groovy title="hello-containers.nf" linenums="19"
        main:
        // cria um canal para entradas de um arquivo CSV
        greeting_ch = channel.fromPath(params.input)
                            .splitCsv()
                            .map { line -> line[0] }
        // emite uma saudação
        sayHello(greeting_ch)
        // converte a saudação para maiúsculas
        convertToUpper(sayHello.out)
        // coleta todas as saudações em um arquivo
        collectGreetings(convertToUpper.out.collect(), params.batch)
    ```

Note que declaramos um novo parâmetro CLI, `params.character`, para especificar qual personagem queremos que diga as saudações.

#### 2.2.3. Adicione o parâmetro `character` ao bloco `params`

Isso é tecnicamente opcional, mas é a prática recomendada e é uma oportunidade de definir um valor padrão para o personagem enquanto estamos nisso.

=== "Depois"

    ```groovy title="hello-containers.nf" linenums="9" hl_lines="7"
    /*
    * Parâmetros do pipeline
    */
    params {
        input: Path = 'data/greetings.csv'
        batch: String = 'batch'
        character: String = 'turkey'
    }
    ```

=== "Antes"

    ```groovy title="hello-containers.nf" linenums="9"
    /*
    * Parâmetros do pipeline
    */
    params {
        input: Path = 'data/greetings.csv'
        batch: String = 'batch'
    }
    ```

Agora podemos ser preguiçosos e pular a digitação do parâmetro character em nossas linhas de comando.

#### 2.2.4. Atualize as saídas do fluxo de trabalho

Precisamos atualizar as saídas do fluxo de trabalho para publicar a saída do processo `cowpy`.

##### 2.2.4.1. Atualize a seção `publish:`

No `bloco workflow`, faça a seguinte alteração de código:

=== "Depois"

    ```groovy title="hello-containers.nf" linenums="34" hl_lines="6"
        publish:
        first_output = sayHello.out
        uppercased = convertToUpper.out
        collected = collectGreetings.out.outfile
        batch_report = collectGreetings.out.report
        cowpy_art = cowpy.out
    ```

=== "Antes"

    ```groovy title="hello-containers.nf" linenums="34"
        publish:
        first_output = sayHello.out
        uppercased = convertToUpper.out
        collected = collectGreetings.out.outfile
        batch_report = collectGreetings.out.report
    ```

O processo `cowpy` produz apenas uma saída, então podemos nos referir a ela da maneira usual adicionando `.out`.

Mas por enquanto, vamos terminar de atualizar as saídas no nível do fluxo de trabalho.

##### 2.2.4.2. Atualize o bloco `output`

Precisamos adicionar a saída final `cowpy_art` ao bloco `output`. Enquanto estamos nisso, vamos também editar os destinos de publicação, já que agora nosso pipeline está completo e sabemos quais saídas realmente nos importam.

No bloco `output`, faça as seguintes alterações de código:

=== "Depois"

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

=== "Antes"

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

Agora as saídas publicadas estarão um pouco mais organizadas.

#### 2.2.5. Execute o fluxo de trabalho

Só para recapitular, isso é o que estamos buscando:

<figure class="excalidraw">
--8<-- "docs/en/docs/hello_nextflow/img/hello_pipeline_complete.svg"
</figure>

Você acha que vai funcionar?

Vamos deletar as saídas publicadas anteriores para ter uma tela limpa, e executar o fluxo de trabalho com a flag `-resume`.

```bash
rm -r hello_containers/
nextflow run hello-containers.nf -resume
```

??? failure "Saída do comando (editada para clareza)"

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

Oh não, há um erro!
O código de erro dado por `error exit status (127)` significa que o executável que pedimos não foi encontrado.

Isso faz sentido, já que estamos chamando a ferramenta `cowpy` mas ainda não especificamos um contêiner (ops).

### 2.3. Use um contêiner para executar o processo `cowpy`

Precisamos especificar um contêiner e dizer ao Nextflow para usá-lo para o processo `cowpy()`.

#### 2.3.1. Especifique um contêiner para `cowpy`

Podemos usar a mesma imagem que estávamos usando diretamente na primeira seção deste tutorial.

Edite o módulo `cowpy.nf` para adicionar a diretiva `container` à definição do processo da seguinte forma:

=== "Depois"

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

=== "Antes"

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

Isso diz ao Nextflow que _se o uso do Docker estiver habilitado_, ele deve usar a imagem de contêiner especificada aqui para executar o processo.

#### 2.3.2. Habilite o uso do Docker via arquivo `nextflow.config`

Note que dissemos _'se o uso do Docker estiver habilitado'_. Por padrão, não está, então precisamos dizer ao Nextflow que é permitido usar Docker.
Para isso, vamos antecipar ligeiramente o tópico da próxima e última parte deste curso (Parte 6), que cobre configuração.

Uma das principais maneiras que o Nextflow oferece para configurar a execução do fluxo de trabalho é usar um arquivo `nextflow.config`.
Quando tal arquivo está presente no diretório atual, o Nextflow o carregará automaticamente e aplicará qualquer configuração que ele contenha.

Fornecemos um arquivo `nextflow.config` com uma única linha de código que explicitamente desabilita o Docker: `docker.enabled = false`.

Agora, vamos mudar isso para `true` para habilitar o Docker:

=== "Depois"

    ```console title="nextflow.config" linenums="1" hl_lines="1"
    docker.enabled = true
    ```

=== "Antes"

    ```console title="nextflow.config" linenums="1" hl_lines="1"
    docker.enabled = false
    ```

!!! tip "Dica"

    É possível habilitar a execução do Docker a partir da linha de comando, por execução, usando o parâmetro `-with-docker <container>`.
    No entanto, isso só nos permite especificar um contêiner para todo o fluxo de trabalho, enquanto a abordagem que acabamos de mostrar permite especificar um contêiner diferente por processo.
    Isso é melhor para modularidade, manutenção de código e reprodutibilidade.

#### 2.3.3. Execute o fluxo de trabalho com Docker habilitado

Execute o fluxo de trabalho com a flag `-resume`:

```bash
nextflow run hello-containers.nf -resume
```

??? success "Saída do comando"

    ```console
     N E X T F L O W   ~  version 25.10.2

    Launching `hello-containers.nf` [drunk_perlman] DSL2 - revision: abf1dccf7f

    executor >  local (1)
    [c9/f5c686] sayHello (3)       [100%] 3 of 3, cached: 3 ✔
    [ef/3135a8] convertToUpper (3) [100%] 3 of 3, cached: 3 ✔
    [7f/f435e3] collectGreetings   [100%] 1 of 1, cached: 1 ✔
    [98/656c6c] cowpy              [100%] 1 of 1 ✔
    ```

Desta vez realmente funciona!
Como de costume, você pode encontrar as saídas do fluxo de trabalho no diretório de resultados correspondente, embora desta vez elas estejam um pouco mais organizadas, com apenas o relatório e a saída final no nível superior, e todos os arquivos intermediários empurrados para fora do caminho em um subdiretório.

??? abstract "Conteúdo do diretório"

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

A saída final de arte ASCII está no diretório `results/hello_containers/`, sob o nome `cowpy-COLLECTED-batch-output.txt`.

??? abstract "Conteúdo do arquivo"

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

E aí está, nosso lindo peru dizendo as saudações como desejado.

#### 2.3.4. Inspecione como o Nextflow lançou a tarefa containerizada

Como uma coda final para esta seção, vamos dar uma olhada no subdiretório de trabalho para uma das chamadas do processo `cowpy` para obter um pouco mais de insight sobre como o Nextflow trabalha com contêineres nos bastidores.

Verifique a saída do seu comando `nextflow run` para encontrar o caminho para o subdiretório de trabalho do processo `cowpy`.
Olhando o que obtivemos para a execução mostrada acima, a linha de log do console para o processo `cowpy` começa com `[98/656c6c]`.
Isso corresponde ao seguinte caminho de diretório truncado: `work/98/656c6c`.

Nesse diretório, você encontrará o arquivo `.command.run` que contém todos os comandos que o Nextflow executou em seu nome no curso da execução do pipeline.

??? abstract "Conteúdo do arquivo"

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

Se você procurar por `nxf_launch` neste arquivo, você deve ver algo assim:

```console
nxf_launch() {
    docker run -i --cpu-shares 1024 -e "NXF_TASK_WORKDIR" -v /workspaces/training/hello-nextflow/work:/workspaces/training/hello-nextflow/work -w "$NXF_TASK_WORKDIR" --name $NXF_BOXID community.wave.seqera.io/library/cowpy:1.1.5--3db457ae1977a273 /bin/bash -ue /workspaces/training/hello-nextflow/work/98/656c6c90cce1667c094d880f4b6dcc/.command.sh
}
```

Como você pode ver, o Nextflow está usando o comando `docker run` para lançar a chamada do processo.
Ele também monta o subdiretório de trabalho correspondente no contêiner, define o diretório de trabalho dentro do contêiner adequadamente e executa nosso script bash modelado no arquivo `.command.sh`.

Todo o trabalho duro que tivemos que fazer manualmente na primeira seção? O Nextflow faz isso por nós nos bastidores!

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

### Conclusão

Você sabe como usar contêineres no Nextflow para executar processos.

### O que vem a seguir?

Faça uma pausa!

Quando estiver pronto, prossiga para [**Parte 6: Hello Config**](./06_hello_config.md) para aprender como configurar a execução do seu pipeline para se adequar à sua infraestrutura, bem como gerenciar a configuração de entradas e parâmetros.

É a última parte, e então você terá concluído este curso!

---

## Quiz

<quiz>
O que é um contêiner?
- [ ] Um tipo de máquina virtual
- [ ] Um formato de compressão de arquivo
- [x] Uma unidade executável leve e autônoma que inclui tudo o que é necessário para executar uma aplicação
- [ ] Um protocolo de rede
</quiz>

<quiz>
Qual é a diferença entre uma imagem de contêiner e uma instância de contêiner?
- [ ] São a mesma coisa
- [x] Uma imagem é um modelo; uma instância é um contêiner em execução criado a partir dessa imagem
- [ ] Uma instância é um modelo; uma imagem é um contêiner em execução
- [ ] Imagens são para Docker; instâncias são para Singularity
</quiz>

<quiz>
O que a flag `-v` faz em um comando `docker run`?
- [ ] Habilita saída detalhada
- [ ] Valida o contêiner
- [x] Monta um volume do sistema host no contêiner
- [ ] Especifica a versão do contêiner

Saiba mais: [1.3.4. Monte dados no contêiner](#134-mount-data-into-the-container)
</quiz>

<quiz>
Por que você precisa montar volumes ao usar contêineres?
- [ ] Para melhorar o desempenho do contêiner
- [ ] Para economizar espaço em disco
- [x] Porque os contêineres estão isolados do sistema de arquivos do host por padrão
- [ ] Para habilitar rede

Saiba mais: [1.3.4. Monte dados no contêiner](#134-mount-data-into-the-container)
</quiz>

<quiz>
Como você especifica um contêiner para um processo Nextflow?
- [ ] `docker 'container-uri'`
- [ ] `image 'container-uri'`
- [x] `container 'container-uri'`
- [ ] `use 'container-uri'`

Saiba mais: [2.3.1. Especifique um contêiner para cowpy](#231-specify-a-container-for-cowpy)
</quiz>

<quiz>
Qual configuração do `nextflow.config` habilita o Docker para seu fluxo de trabalho?
- [ ] `#!groovy process.docker = true`
- [x] `#!groovy docker.enabled = true`
- [ ] `#!groovy container.engine = 'docker'`
- [ ] `#!groovy docker.activate = true`

Saiba mais: [2.3.2. Habilite o uso do Docker via arquivo `nextflow.config`](#232-enable-use-of-docker-via-the-nextflowconfig-file)
</quiz>

<quiz>
O que o Nextflow gerencia automaticamente ao executar um processo em um contêiner? (Selecione todas as opções aplicáveis)
- [x] Puxar a imagem do contêiner se necessário
- [x] Montar o diretório de trabalho
- [x] Executar o script do processo dentro do contêiner
- [x] Limpar a instância do contêiner após a execução

Saiba mais: [2.3.4. Inspecione como o Nextflow lançou a tarefa containerizada](#234-inspect-how-nextflow-launched-the-containerized-task)
</quiz>
