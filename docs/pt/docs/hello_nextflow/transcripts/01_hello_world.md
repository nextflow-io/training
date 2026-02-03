# Parte 1: Olá Mundo - Transcrição

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Tradução assistida por IA - [saiba mais e sugira melhorias](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

<div class="video-wrapper">
  <iframe width="560" height="315" src="https://www.youtube.com/embed/8X2hHI-9vms?si=F0t9LFYLjAWoyRXj&amp;list=PLPZ8WHdZGxmXiHf8B26oB_fTfoKQdhlik" title="YouTube video player" frameborder="0" allow="accelerometer; autoplay; clipboard-write; encrypted-media; gyroscope; picture-in-picture; web-share" referrerpolicy="strict-origin-when-cross-origin" allowfullscreen></iframe>
</div>

!!!note "Notas importantes"

    Esta página mostra apenas a transcrição. Para instruções passo a passo completas, retorne ao [material do curso](../01_hello_world.md).

    Os números das seções mostrados na transcrição são fornecidos apenas para fins indicativos e podem não incluir todos os números de seção dos materiais.

## Boas-vindas

Olá, bem-vindo ao Capítulo Um de Olá Nextflow.

Nesta primeira parte de um curso de seis partes, vamos entrar nos conceitos mais básicos do Nextflow. Vamos começar executando alguns comandos em um terminal e então vamos pegar esses comandos Bash e ver como transformá-los em um script Nextflow.

Vamos tentar executar esse primeiro pipeline Nextflow, ver o que o Nextflow faz, onde ele executa, quais arquivos ele cria e qual é o propósito desses arquivos.

Tudo bem, vamos começar.

## training.nextflow.io

Primeiro de tudo, vá para training.nextflow.io. Assim como antes, todo o material está escrito aqui, e vou trabalhar através dele passo a passo. Vou mostrar minha tela enquanto faço as etapas do treinamento, mas tudo o que estou dizendo está no material de treinamento, então você pode seguir no seu próprio ritmo, e pode encontrar tudo escrito lá.

Este vídeo também tem legendas habilitadas, então sinta-se à vontade para ativá-las e acompanhar exatamente o que estou dizendo enquanto falo.

Ok, vamos para Olá Nextflow. Esse é o curso que vamos fazer hoje, e já fizemos a orientação no primeiro vídeo, então vamos direto para a parte um: Olá Mundo.

Ok, vou sair deste material de treinamento agora e entrar no meu ambiente Codespaces. Isso é o que configuramos no primeiro vídeo. Espero que você tenha algo muito similar a isso em seu próprio sistema. Estou usando VS Code e estou olhando o material de treinamento e mudei de diretório para o diretório hello Nextflow.

## 0. Aquecimento: Execute Olá Mundo diretamente

Ok. Vamos começar com alguns conceitos básicos, que esperançosamente serão familiares para todos. Vou começar apenas escrevendo um comando muito básico no terminal. Aqui embaixo vou dizer 'echo Hello World!"' pressiono enter e, sem surpresas, o terminal faz o que peço e retorna aquela string. Hello world.

Ok, então vou pressionar para cima para obter aquele comando e editá-lo um pouco mais. Vamos desta vez redirecionar aquela saída para um arquivo. Vou escrevê-la em vez disso para output.txt e pressionar enter, nada no terminal desta vez porque a saída não veio para o terminal. Ela foi para aquele arquivo.

Posso então ler aquele arquivo fazendo 'cat output.txt' pressiono tab aqui para expandir automaticamente o nome do arquivo e pronto. O arquivo está lá.

Também posso ver aquele arquivo na barra lateral no explorador de arquivos no VS Code. Posso clicar duas vezes nele e abri-lo aqui. Se você quiser abri-lo no VS Code sem clicar em nada, também pode fazer "code" e então "output.txt" e faz a mesma coisa.

Ótimo. Esse é o primeiro passo. Muito simples.

## 1. Examine o script inicial do fluxo de trabalho Olá Mundo

Ok. Agora vamos fazer exatamente a mesma coisa, mas no Nextflow, em vez de diretamente no terminal.

Vamos usar o primeiro script de exemplo para começar, este arquivo é chamado Hello World. Posso fazer "ls" para visualizá-lo em um terminal, e estou no Mac, então posso fazer command click para abrir aquele arquivo, ou poderia ter apenas clicado duas vezes na barra lateral aqui.

Há algumas coisas que podemos ver neste arquivo. Logo no topo, há uma declaração com hash dizendo que este é um arquivo Nextflow e é assim que pode ser executado. Há alguns comentários aqui, apenas comentários regulares de código em cinza claro, que não afetam a execução e apenas nos ajudam a ler o script.

E então há duas estruturas principais. Há um processo aqui e um fluxo de trabalho.

Processos no Nextflow são as etapas do pipeline. São as partes que realmente fazem a lógica e fazem o processamento.

O fluxo de trabalho então na parte inferior conecta esses processos juntos e governa a lógica do fluxo de trabalho, como tudo se conecta um ao outro.

Vamos começar olhando um processo. Voltaremos ao fluxo de trabalho em um momento.

## 1.2 A definição do processo

Então todo processo começa com uma palavra-chave process. Tem um nome e então tem algumas chaves e tudo dentro dessas chaves é aquele único processo.

Um processo deve ter uma seção script, e contida aqui está um trecho bash em uma string de múltiplas linhas, que é a parte do código que é realmente executada no ambiente de computação.

Também temos uma declaração de saída aqui, que diz ao Nextflow quais arquivos são esperados serem criados pelo script. Note que a saída aqui tem uma palavra-chave path, que diz ao Nextflow que isso é um arquivo, não um valor ou uma string.

Dentro do bloco script, isso é apenas uma declaração bash regular, e é exatamente a mesma que escrevemos no terminal. Estamos ecoando hello world para um arquivo chamado output.txt. Este output.txt é então capturado pela definição de saída. A definição de saída não está realmente fazendo nada. Está apenas dizendo ao Nextflow o que esperar, e se este arquivo não fosse criado, o Nextflow lançaria um erro.

Note que este exemplo não é muito bom porque codificamos fixamente o nome do arquivo aqui, output.txt e output.txt. Se qualquer um destes fosse alterado, isso causaria um erro em nosso fluxo de trabalho.

Há uma maneira melhor de fazer isso com variáveis, que vamos cobrir em um minuto.

## 1.3 A definição do fluxo de trabalho

Ok. Descendo para o fluxo de trabalho, podemos ver que temos um comentário e então executamos o processo chamado sayHello. Esta é a mesma palavra-chave que está aqui em cima. Isso é tão simples quanto um fluxo de trabalho pode ser. Estamos apenas chamando um único processo sem entrada variável, então não estamos conectando-o a mais nada. Na parte posterior deste curso, vamos falar sobre como tornar isso mais poderoso usando entradas variáveis e conectando coisas com canais.

## 2. Execute o fluxo de trabalho

Ok, isso é tudo que precisamos. Vamos ver se podemos executá-lo e ver o que acontece. Vou apenas limpar o terminal e então vou fazer "nextflow run", e vou chamar o nome do arquivo, que é hello-world.nf. Isso é tudo que precisamos para executar um pipeline Nextflow. Este pipeline não recebe nenhuma entrada, então não precisamos de nenhum outro argumento.

Vamos pressionar enter e ver o que acontece.

Ok. Esperançosamente você deve ter alguma saída que se parece com isso. Temos alguns pedaços de informação nos dizendo que o Nextflow executou e qual versão estava usando. Nos diz qual script foi lançado e nos dá um nome gerado aleatoriamente para esta execução particular do fluxo de trabalho. Neste caso, o meu foi chamado "gloomy_crick".

A parte mais importante disso, porém, é que nos diz quais etapas executaram no pipeline. Você pode ver que nosso processo chamado sayHello executou, e executou uma vez e estava cem por cento completo.

Esta parte aqui é o hash para aquela tarefa particular do fluxo de trabalho. Cada processo executa uma ou mais vezes, e cada uma dessas execuções é chamada de tarefa.

## 2.2. Encontre a saída e os logs no diretório work

Cada tarefa obtém seu próprio diretório isolado onde executa, então está separada do resto da execução do fluxo de trabalho. Este hash corresponde à estrutura de arquivos dentro do diretório work. Se eu fizer "tree work", podemos ver a0, e então uma versão mais longa de um hash curto, e então nosso arquivo output.txt. Você também pode vê-lo na barra lateral.

Você pode ver na barra lateral que há alguns arquivos adicionais aqui. A razão pela qual estes não apareceram no terminal é porque são arquivos ocultos, eles começam com um ponto. E de fato, se eu fizer "tree -a" para todos, e "work", podemos vê-los aqui.

Estes arquivos com ponto estão presentes em cada único diretório work que o Nextflow cria, e cada um tem uma tarefa ligeiramente diferente. Primeiro .command.begin apenas inclui algumas instruções para o Nextflow que configura a tarefa antes de executar. .command.run são as instruções reais executadas pelo próprio Nextflow. Então .command.sh é provavelmente o mais interessante. Este é o script que foi resolvido do nosso bloco script do processo.

Se eu abri-lo, você pode ver que temos nosso "echo Hello World" para o arquivo output.txt. Isso é exatamente o mesmo que nosso processo neste caso, mas se tivermos quaisquer variáveis dentro do nosso código Nextflow, cada tarefa terá um .command.sh diferente, e você pode ver como essas variáveis foram resolvidas.

Os outros arquivos são sobre como a tarefa executou. Então .command.err, .log e .out são o erro padrão, saída padrão e os dois combinados. E .exitcode diz ao Nextflow como esta tarefa executou com qual código de saída, se foi bem-sucedida ou não.

Finalmente, temos nosso arquivo output.txt e com certeza, "Hello World" isso é o que estávamos esperando e isso é o que foi criado.

Ok, ótimo. Essa foi sua primeira execução Nextflow. Parabéns. É realmente tão simples assim.

A seguir, vamos ver como fazer isso de uma maneira um pouco mais conveniente para que não tenhamos que editar o código toda vez que quisermos fazer uma mudança em como o pipeline executa.

## 3. Gerencie as execuções do fluxo de trabalho

Esta estrutura de diretório é ótima para manter todas as tarefas separadas e tudo organizado, mas claro, não é muito conveniente para encontrar seus arquivos de saída. Você não quer estar cavando através de muitos diretórios aninhados tentando encontrar os resultados do seu pipeline.

## 3.1. Publique saídas

A boa notícia é que você não deveria. Os diretórios work são realmente apenas para o próprio Nextflow usar. Então o que vamos fazer é usar uma função do Nextflow chamada "publishDir".

Voltamos ao nosso fluxo de trabalho, vamos ao processo. Podemos adicionar uma nova declaração aqui chamada diretiva. Isso é o que o Nextflow chama essas coisas no topo dos processos que aumentam como a funcionalidade funciona, e a que vamos usar é chamada publishDir.

Você pode ver que comecei a digitar aqui e a extensão Nextflow para VS Code sugeriu a diretiva para mim, então posso apenas pressionar enter.

Ok. Vou seguir isso com um diretório chamado "results" e vamos dizer para copiar os arquivos de saída lá. Então vou dizer mode copy. Ótimo. Vou salvar e vamos executar o fluxo de trabalho novamente.

nextflow run hello-world.nf

Executa exatamente da mesma forma. Embora note que temos um hash ligeiramente diferente desta vez. O Nextflow usará um hash diferente toda vez que você executar o fluxo de trabalho. E temos um conjunto diferente de diretórios work como resultado. Áreas, um chamado EB em vez disso, mas você pode ver que todos os arquivos são os mesmos. No entanto, o que é novo desta vez é que também temos um diretório chamado "results".

Dentro de "results" aqui temos nosso arquivo de saída. Isso é o que dissemos ao Nextflow para fazer. Dissemos, salve os arquivos de resultados em um diretório chamado "results" e copie-os lá. E então isso agora é muito mais fácil de encontrar. Está apenas lá ao lado de onde lançamos o fluxo de trabalho e todos os diferentes arquivos podem ser organizados lá como quisermos, independentemente de onde ou como o Nextflow executou a execução real.

Note que publishDir pode lidar com links simbólicos, o que é bom se você está trabalhando em um sistema de arquivos compartilhado e quer economizar espaço. E também você não tem que definir todos os arquivos que são criados por um processo como uma saída.

O Nextflow só copiará as coisas que estão definidas neste bloco output. Então, se você tem arquivos intermediários criados pela etapa, que não são necessários a jusante deste processo, você simplesmente não os define na saída e eles não aparecerão no publishDir. Então esta é uma maneira de manter seus arquivos de saída de um pipeline limpos e facilmente deletar arquivos intermediários uma vez que o local de trabalho tenha terminado.

Uma nota rápida aqui. Há uma nova sintaxe Nextflow chegando chamada definições de saída de fluxo de trabalho, que eventualmente substituirá publishDir. Isso nos dá uma maneira de definir todas as saídas de um fluxo de trabalho em nível de pipeline no bloco workflow. Isso é descrito na documentação do Nextflow se você quiser experimentar. Mas por enquanto, publishDir ainda estará por aí por um tempo, então ainda temos isso no treinamento para 2025.

## 3.2. Relance um fluxo de trabalho com -resume

Ok. Mencionei que o diretório work aqui agora tem dois conjuntos de resultados com um hash diferente de cada vez que executamos o fluxo de trabalho. Isso é bom. No entanto, às vezes não queremos recomputar etapas toda vez se não precisarmos.

Talvez você esteja construindo iterativamente seu fluxo de trabalho e esteja adicionando etapas e queira que as primeiras etapas apenas reutilizem as versões em cache. Ou talvez algo deu errado em seu sistema de computação no meio do seu fluxo de trabalho e você quer que ele continue de onde parou, mas pule as etapas que já tinha completado.

O Nextflow tem funcionalidade embutida para isso chamada resume. Vamos experimentar. Então primeiro, vou apenas dar uma olhada no diretório work para que possamos lembrar o que estava lá.

E então vou fazer "nextflow run hello-world.nf" e vou adicionar um único comando aqui, "-resume".

Note, um único traço, isso é realmente importante. Vou executá-lo e a saída vai parecer basicamente exatamente a mesma, com algumas pequenas diferenças.

Note aqui que diz "cached" em cinza. Isso significa que o Nextflow não executou a tarefa. Desta vez ele encontrou algo que correspondia ao que eram requisitos e reutilizou essas saídas diretamente em vez de reexecutar a etapa.

E com certeza, se você olhar o hash aqui, você pode ver que isso corresponde ao hash existente que tínhamos de uma execução anterior.

## 3.3. Delete diretórios work mais antigos

Ok. Mas se você está desenvolvendo iterativamente, vai acumular muitos desses arquivos de fluxo de trabalho. Isso pode ser um problema se você pode estar com pouco espaço.

O Nextflow pode nos ajudar a limpar esses diretórios work com alguns comandos auxiliares. Se eu fizer "nextflow log". Isso me dará uma lista de todas as diferentes execuções de fluxo de trabalho que fiz neste diretório, e elas têm os nomes de execução aqui. Você pode ver o gloomy quick, que foi o primeiro que executamos, e então esses dois novos.

Agora podemos pegar esse nome e usá-los com o comando "nextflow clean". Posso especificar um único nome de execução. Ou ainda melhor, posso dizer ao Nextflow para deletar tudo de antes de um único nome de fluxo de trabalho com "-before", e vou colocar "stupefied_shaw". Essa foi minha execução mais recente, "-n".

O comando "-n" disse ao Nextflow para fazer isso como uma execução de teste sem realmente deletar nada de verdade, e nos diz quais dos diretórios hash teriam sido removidos. Com certeza, é apenas aquele da primeira execução. Ambas as segundas execuções usam o mesmo diretório hash.

Vou executá-lo novamente, mas agora em vez de "-n" para execução de teste, vou fazer "-f" para forçar e ele removeu aquele diretório hash. Agora se eu fizer "tree work", podemos ver, temos apenas este arquivo de saída restante.

Ótimo. Então conseguimos limpar um monte de espaço em disco ali.

Algumas coisas a notar ao deletar diretórios work, se você fizer links simbólicos para seu diretório de resultados, essas fontes de links simbólicos agora serão deletadas e seus resultados se perderão para sempre. Então é por isso que usar o modo copy é uma coisa mais segura de fazer, e geralmente o que recomendamos.

Em segundo lugar, a funcionalidade resume do Nextflow depende desses diretórios work. Então, se você deletá-los e executar o Nextflow novamente, a funcionalidade resume não funcionará mais. Então cabe a você acompanhar quais coisas você pode precisar ou não precisar, e apenas deletar coisas quando você tem certeza de que é seguro fazê-lo.

A outra coisa que podemos fazer é apenas deletar o diretório work inteiro se terminarmos nossa execução de fluxo de trabalho e tivermos certeza de que não precisamos mais dele.

Então posso fazer "rm -r work". Eu sei que não havia nada importante lá. Tenho meus resultados com os quais me importo no diretório de resultados onde os copiamos. E então foi seguro deletar o diretório work. Cabe a você qual dessas abordagens usar.

## 4. Use uma entrada variável passada na linha de comando

Ok, qual é o próximo passo? Mencionei que tínhamos codificado fixamente alguns dos valores em nosso script de fluxo de trabalho aqui, o arquivo output.txt, e que pode haver uma maneira melhor de fazer isso.

Vamos começar com isso. O que vamos fazer são três coisas. Vamos adicionar uma nova entrada ao processo. Vamos dizer ao script do processo como usar essa entrada, e então vamos conectá-la no fluxo de trabalho para que possamos usá-la dinamicamente com uma flag de linha de comando ao executar o Nextflow.

Então, primeiro de tudo. Vamos adicionar um bloco input aqui. Assim como output. Esta é uma nova seção para o processo, e vou dizer, "val greeting".

Note aqui, estou dizendo "val", que diz que isso é uma variável, não um path.

Posso então descer no script e então posso tirar este texto codificado fixamente aqui e fazer $greeting. Isso funciona como qualquer outra linguagem de programação. Estamos definindo uma variável aqui e estamos referenciando-a dentro deste bloco script. Quando o Nextflow executa este processo, a variável será interpolada. E quando formos olhar aquele arquivo .command.sh, veremos a string codificada fixamente real aqui em vez disso.

## 4.1.3. Configure um parâmetro CLI e forneça-o como entrada para a chamada do processo

Ok, mas onde fornecemos a variável? A seguir descemos para a seção workflow, e você pode ver que a extensão aqui está dizendo, agora esperamos uma entrada, e me deu um aviso.

Agora, a coisa mais simples que poderíamos fazer é apenas codificá-la fixamente. Eu poderia escrever "Hello World" e fornecer aquela entrada de string para o processo. Mas novamente, isso realmente não resolveria nenhum problema. Ainda teríamos que voltar e editar o código do pipeline toda vez que quiséssemos mudar algo, o que não é bom.

A boa notícia é que o Nextflow tem um sistema embutido para lidar com argumentos de linha de comando chamado parâmetros. Então, em vez disso, posso usar uma dessas variáveis especiais chamada params e posso chamá-la do que eu quiser, mas vou dizer greeting para que corresponda à lógica do fluxo de trabalho.

Salvar e vamos ver o que podemos fazer com isso.

Então, se eu voltar ao terminal. Então fazemos "nextflow run hello-world.nf". Assim como antes, mas a diferença chave é que fazemos --greeting

Note, há dois traços aqui porque este é um parâmetro. Quando retomamos o fluxo de trabalho antes, era um único traço. Isso é porque resume é uma opção central do Nextflow, e este é um parâmetro que é específico para nosso pipeline.

Não misture os dois. É fácil fazer isso. Se você fizesse --resume em vez de apenas um traço, então isso seria "params.resume", que não faria nada. Da mesma forma, se você fizesse um único traço aqui, o Nextflow não reconheceria como um argumento chave.

Então é --greeting, que corresponde a params.greeting.

Agora posso seguir isso com qualquer texto que eu quiser. Então estou na Suécia no momento, então vou dizer, "Hej världen".

Então vamos executar, ver o que acontece, momento da verdade.

Ok, então você pode ver que o processo executou novamente, assim como antes, sayHello com uma única execução.

Isso terá sobrescrito o arquivo que estava no diretório publishDir "results". E então tenha cuidado quando estiver reexecutando os arquivos porque coisas no publishDir serão sobrescritas.

Agora posso fazer "code results/output.txt", e com certeza, nossa saída foi atualizada e agora diz "Hej världen".

## 4.2. Use valores padrão para parâmetros de linha de comando

Ok, isso é ótimo. Mas o problema agora é que nosso fluxo de trabalho depende de sempre definirmos este parâmetro, e é bom ter padrões sensatos para que as coisas executem de uma maneira sensata para seu fluxo de trabalho, a menos que você sobrescreva os padrões.

Então a maneira que fazemos isso é definindo um valor padrão para o parâmetro em nosso script de fluxo de trabalho.

Então, se eu voltar ao meu arquivo hello-world.nf, posso ir no script logo acima de workflow, digitar "params.greeting" e defini-lo como qualquer outra variável. Então vamos colocar uma string aqui e vamos dizer "Holà mundo!"

Agora este parâmetro tem um padrão definido, que será usado aqui, ou ainda podemos sobrescrevê-lo na linha de comando com --greeting, assim como fizemos antes.

Então vamos verificar se funciona. "nextflow run hello-world.nf"

Sem argumentos de linha de comando desta vez, e verificar se fez a coisa certa.

"code results/output.txt". E lá está. Obtivemos nosso padrão.

Ok, vamos tentar novamente, apenas verificar que não estou contando nenhuma mentira. Vamos executar novamente, mas fazer --greeting, e usar o exemplo do material de treinamento, vamos dizer "Konnichiwa!"

Reexecuta o fluxo de trabalho, e com certeza, nosso arquivo de saída no topo foi apenas atualizado com o novo valor que fornecemos na linha de comando.

Ótimo. Este é um aspecto realmente central para escrever qualquer fluxo de trabalho Nextflow. Definir padrões sensatos em seu código de pipeline, mas tornar muito fácil de configurar para o usuário final tendo argumentos de linha de comando no terminal.

Note que o usuário final pode sobrescrever a configuração em múltiplos lugares diferentes. Você pode ter um arquivo de configuração em seu diretório home, que é aplicado a cada única execução Nextflow que você faz. Você pode ter um arquivo de configuração em um diretório de lançamento. Você pode ter um arquivo de configuração em um diretório de pipeline. Todos esses diferentes locais de configuração são carregados em uma ordem específica, que é descrita na documentação do Nextflow.

Ok, esse é o fim da seção um. Tivemos nosso primeiro script de fluxo de trabalho no Nextflow com um processo e um fluxo de trabalho. Olhamos entradas, saídas, scripts e publicação, e como conectar parâmetros e um canal de entrada em nosso processo.

Parabéns, seu primeiro passo em direção a escrever código Nextflow está completo.

Faça uma pequena pausa e te vejo de volta em alguns minutos para o capítulo dois.

[Próxima transcrição de vídeo :octicons-arrow-right-24:](02_hello_channels.md)
