# Parte 1: Hello World - Transcrição do Vídeo

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Tradução assistida por IA - [saiba mais e sugira melhorias](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

<div class="video-wrapper">
  <iframe width="560" height="315" src="https://www.youtube.com/embed/tOukLxWCHiA?si=F0t9LFYLjAWoyRXj&amp;list=PLPZ8WHdZGxmWKozQuzr27jyMGqp9kElVK&amp;cc_load_policy=1" title="YouTube video player" frameborder="0" allow="accelerometer; autoplay; clipboard-write; encrypted-media; gyroscope; picture-in-picture; web-share" referrerpolicy="strict-origin-when-cross-origin" allowfullscreen></iframe>
</div>

!!!note "Notas importantes"

    Esta página mostra apenas a transcrição. Para instruções passo a passo completas, retorne ao [material do curso](../01_hello_world.md).

    Os números das seções mostrados na transcrição são fornecidos apenas para fins indicativos e podem não incluir todos os números de seção dos materiais.

## Boas-vindas

Olá, e bem-vindo de volta.

Você está agora na Parte Um do curso "Hello Nextflow" chamada "Hello World". Neste capítulo, vamos começar a construir algum entendimento dos conceitos mais básicos do Nextflow.

Então, esperamos que você esteja agora configurado no Codespaces ou em algum lugar equivalente com o VS Code rodando, e você tem sua pasta Hello Nextflow no workspace no Explorer com todos esses arquivos diferentes aqui.

Vamos começar fazendo apenas algumas coisas muito básicas no terminal usando Bash, e então vamos ver se conseguimos fazer as mesmas coisas dentro do Nextflow para que você tenha uma noção de como a sintaxe se parece.

## 0. Aquecimento

Então vamos começar muito, muito simples. Vamos começar apenas com "echo", para imprimir algo em um terminal. "Hello World". Eu pressiono enter e isso vai para um terminal. Hello World. Esperamos que isso não seja uma surpresa para ninguém assistindo este curso.

Ok, vamos fazer algo com isso. Em vez de apenas imprimi-lo no terminal, vamos escrever em um arquivo. Vou pressionar a seta para cima no meu teclado, que percorre o histórico do Bash, então me dá meu último comando, e vou adicionar no final dele ali, um pequeno símbolo de maior que, que redireciona a saída deste comando para um arquivo, e vou chamá-lo de output.txt.

Enter novamente, para executar esse comando, nada no terminal desta vez, mas podemos ver no lado esquerdo que o novo arquivo apareceu aqui, chamado output.txt.

Podemos visualizá-lo em um terminal com algo como cat. Então cat output.txt e com certeza diz "Hello World". Também podemos dar um clique duplo nele e ele abre no editor de código no VS Code.

## 1.1. Examine o código

Tudo bem. Eu disse que era simples. O que vem a seguir? Vamos tentar pegar este processo e fazer novamente, mas desta vez, vamos fazê-lo dentro do Nextflow.

Como eu disse, todos os diferentes capítulos neste curso começam com um script e este se chama Hello World. Então vou encontrar Hello World. Ele mostra uma prévia se eu clicar uma vez nele, vou dar um clique duplo para abri-lo no editor aqui. E vou rapidamente me livrar do terminal.

Agora este é um script muito, muito simples, tão simples quanto possível. Tem apenas 22 linhas de comprimento, e faz basicamente a mesma coisa. Na verdade. Parte disso deve parecer familiar. É o que acabamos de digitar. Podemos ver nosso comando bash redirecionando para um arquivo ali.

Ok. O que mais? Também, neste arquivo, podemos começar a ver alguns dos conceitos centrais do Nextflow. Temos um processo em vermelho aqui e um fluxo de trabalho. Estas são palavras-chave especiais e terminologia especial no Nextflow.

## 1.1.1. A definição do processo

Diferentes processos dentro de um fluxo de trabalho envolvem diferentes unidades lógicas do seu fluxo de trabalho. Cada processo faz uma coisa.

Quando executamos, ele gera uma tarefa ou múltiplas tarefas, que são uma etapa real de execução de um pipeline. Todos os processos são então orquestrados dentro de um bloco de fluxo de trabalho, que vemos na parte inferior, e neste caso, apenas executa aquele único processo.

O nome do processo segue esta palavra-chave aqui, e isso pode ser basicamente qualquer coisa. E então o conteúdo do processo está dentro dessas chaves.

Há apenas realmente um requisito para o processo, que é que ele inclua algum tipo de bloco script ou exec. Está nas aspas triplas aqui, e este é o script bash que é escrito no diretório de trabalho quando executamos o pipeline e há uma coisa que realmente roda no seu computador ou servidor.

Isso é bash tipicamente, mas você também pode colocar uma shebang diferente aqui no topo, e poderia ser um script Python ou um script R. Não importa. O que quer que esteja neste script será executado.

Há uma outra coisa que adicionamos neste processo aqui, que é a declaração de saída. Isso diz ao Nextflow que este processo está esperando um arquivo de saída chamado output.txt. Diz que é um path, então deve ser tratado como um arquivo, não diga, se isso fosse val, diria que é como uma variável ou valor.

Note que isso não está criando este arquivo. Não está realmente gerando ele. Isso é feito pelo script aqui embaixo. Está apenas dizendo ao Nextflow para esperar um arquivo de saída com este nome de arquivo.

## 1.1.2. A definição do fluxo de trabalho

Ok. E então na parte inferior temos um fluxo de trabalho aqui, e novamente, temos uma declaração. Este se chama Main. Este é o equivalente do fluxo de trabalho de um bloco de script, se quiser. É a parte do fluxo de trabalho que faz algo. E neste caso, estamos dizendo, chame o processo chamado sayHello.

Normalmente, é claro, seu pipeline parecerá muito mais complexo do que isso. Você provavelmente terá mais de um processo, e usará canais para orquestrar o fluxo de dados entre eles. Vamos chegar nisso nas próximas partes deste curso, mas por enquanto, isso é suficiente. Este é um pipeline válido, que deve funcionar.

Posso até clicar em preview DAG aqui no VS Code. O DAG ou DAG é uma representação de uma estrutura de fluxo de dados no pipeline, e podemos vê-lo renderizado no lado como um diagrama mermaid. Neste caso é muito, muito simples. Há uma caixa, que é o fluxo de trabalho e um processo, que se chama sayHello, mas isso pode parecer mais interessante conforme avançamos.

## 1.2. Execute o fluxo de trabalho

Ok, vamos tentar executar este fluxo de trabalho e ver o que acontece.

Vou trazer o terminal de volta na parte inferior, limpar a saída, e vou digitar Nextflow Run. E então vou apenas digitar o nome do script, que é hello-world.nf. E vou pressionar enter.

Ok, tem algumas coisas padrão no topo, que nos dizem que o Nextflow executou e qual versão estava rodando e qual era o nome do script e tudo mais.

E realmente a coisa importante que estamos procurando aqui é _aqui_, que é um resumo das diferentes tarefas que foram executadas.

Se o seu parecer com isso com um pequeno marcador verde, então parabéns. Você acabou de executar seu primeiro pipeline. Fantástico.

Diz aqui o nome do processo que executou, que se chamava Say Hello, e nos disse que executou uma vez e que foi bem-sucedido. Isso atualiza conforme você avança, então quando você está executando um pipeline maior, você verá o progresso representado aqui. Mas porque isso é tão pequeno, roda basicamente imediatamente.

## 1.2.2. Encontre a saída e os logs no diretório work

Agora quando você executa um pipeline Nextflow, cada um desses processos é costurado junto, e cada processo, como eu disse antes, pode gerar tarefas uma ou múltiplas. Então neste caso, tivemos uma única tarefa deste processo. Ele apenas executou uma vez e isso foi feito sob este _hash_ de tarefa.

O Nextflow não lida com os arquivos no seu diretório de trabalho diretamente, ele cria uma pasta especial chamada work. E se eu fizer "ls", veremos que apareceu aqui: _work_, e dentro há subdiretórios para cada tarefa única que executa. E isso corresponde a este hash. Então você pode ver se eu for para "ls work/c4", e então está truncado, mas começa 203, e esse é o diretório de trabalho, que foi criado por este processo quando executamos o pipeline. E você pode vê-lo no lado também.

Quando listo esses arquivos, você pode ver que o arquivo output.txt foi gerado. Você pode vê-lo aqui também. E há um monte de arquivos ocultos, que não estão aparecendo com meu "ls" regular.

Se eu clicar em output.txt, com certeza, temos nossa saída. Fantástico. Então o pipeline funcionou.

Pode parecer muito boilerplate para executar o que era essencialmente um script bash de uma linha, mas fará mais sentido conforme nossos processos ficarem mais complicados. E este diretório work com o Nextflow e esses arquivos, que são criados é realmente a espinha dorsal do que torna o Nextflow tão poderoso.

Cada tarefa, cada elemento de um pipeline está isolado de todas as outras tarefas. É reproduzível. Eles não entram em conflito uns com os outros, e tudo pode executar em paralelo. Na verdade é um jeito muito bom quando você se acostuma porque por causa desse isolamento você pode entrar e ver exatamente o que aconteceu para uma única tarefa e depurar.

Vamos dar uma olhada rápida nesses outros arquivos no diretório work. De cima para baixo, temos um arquivo chamado _.command.begin_. Isso está vazio. É apenas o que é chamado de arquivo sentinela, criado pelo Nextflow dizendo, ok, estou começando a tarefa. Nada de interessante ali.

Então há _.command.error_, _.command.log_ e _.command.out_. Estes são todos saídas do comando bash ou deste script que executou. Este é standard error. Este é standard out, e este é os dois combinados conforme saíram. Então você obtém a ordem lógica.

Ok, esses estavam todos vazios para isso também, então não muito interessantes, mas as coisas ficam mais interessantes quando você chega ao _.command.run_.

Este é tipicamente um script muito longo. E isto é o que o Nextflow realmente executa. Se você entrar aqui, começará a ver toda a lógica interna do Nextflow e ver o que ele está fazendo e como está executando seu processo. Isso dependerá de onde você está executando, se estamos executando localmente ou submetendo como um job para SLURM, nesse caso teremos cabeçalhos SLURM no topo. Todas essas diferentes configurações.

Geralmente, você realmente não precisa nunca olhar neste arquivo. É gerado automaticamente pelo Nextflow e não há nada realmente particularmente único para seu pipeline, que está nele. Mas esse é realmente o núcleo do que está rodando.

O próximo é muito mais interessante. _.command.sh_ é o script gerado, que veio do seu processo, e aqui você pode ver que o Nextflow adicionou o cabeçalho Bash, e então executou nosso comando, que estava no nosso bloco de script.

E é tudo que o arquivo _.command.run_ faz é apenas executar este arquivo _.command.sh_.

Este é realmente útil, que é o que você geralmente acaba olhando mais quando está tentando depurar algo e verificar que a lógica do seu pipeline Nextflow está fazendo o que você espera que faça.

Finalmente, temos um arquivo chamado _.exitcode_, e isso apenas captura o código de saída de uma tarefa, que neste caso foi bem-sucedido. Então o código de saída foi zero.

Se algo der errado, você ficar sem memória ou outra coisa e falhar, então isso é muito útil para entender o que deu errado.

## 1.3. Execute o fluxo de trabalho novamente

Mais uma coisa para entender sobre diretórios work é que se eu continuar executando este pipeline repetidamente, então se eu _"nextflow run hello-world.nf"_, vai fazer exatamente a mesma coisa, mas desta vez terá um novo id de tarefa. Você pode ver que este hash aqui é diferente, e agora se eu olhar em work, há dois diretórios de hash. E estes são, novamente, separados um do outro.

Então toda vez que você executa um fluxo de trabalho Nextflow, a menos que você use o resume, que usa o cache, vamos tocar nisso mais tarde, vai reexecutar esses processos em novos diretórios work, que são separados uns dos outros. Você não terá nenhuma colisão de nome de arquivo, você não terá nenhum problema assim. Tudo está isolado e limpo.

E se entrarmos neste diretório, você pode ver todos os mesmos arquivos e o mesmo _output.txt_, que foi recriado do zero.

## 2. Publique saídas

Ok, isso é ótimo para o Nextflow em si, enquanto está executando seu pipeline para que todas as coisas estejam separadas umas das outras e limpas e possam ser gerenciadas.

Mas não é super útil se você é uma pessoa tentando explorar seus resultados. Você realmente não quer estar cavucando através de milhares e milhares de diferentes diretórios work tentando encontrar seus arquivos de resultado. E você realmente não deveria. Os diretórios work não são destinados a ser o estado final de onde seus arquivos são criados.

Fazemos isso publicando nossos arquivos.

## 2.1.1. Declare a saída do processo sayHello

Então se eu voltar ao nosso script, vamos trabalhar no nosso bloco de fluxo de trabalho aqui. Vamos dizer quais arquivos esperar, quais arquivos nos importamos, e então vamos criar um novo bloco embaixo chamado bloco de saída.

Esta é a nova sintaxe, que veio com o analisador sintático e é padrão na versão 26.04 do Nextflow. Então se você usou o Nextflow um pouco antes, esta é uma das coisas que é nova.

Então temos o bloco main, e a seguir vou dizer publish e vou dizer ao Nextflow o que esperar da publicação. Vamos chamar de _first_output_, e vamos chamar de _sayHello.out_.

Eu acidentalmente cometi um erro de digitação ali, mas esta é uma boa oportunidade para também apontar alguns dos recursos da extensão VS Code do Nextflow. Você pode ver que imediatamente me deu uma pequena linha ondulada vermelha embaixo disso dizendo que algo está errado. E se eu passar o mouse sobre ele, vai me dizer que esta variável não está definida. Não sei o que é.

É bastante óbvio neste caso, cometi um erro de digitação. Eu queria digitar sayHello, e então a linha ondulada vai embora.

Agora está roxo. O analisador sintático do Nextflow sabe que este é um processo e quando passo o mouse sobre ele, me dá uma representação reduzida de como este processo se parece. Então posso ver muito rapidamente de relance que não recebe nenhuma entrada e nos dá esta saída. Então trabalhar no VS Code com esta extensão dá muita informação contextual conforme você está escrevendo código.

Note que podemos nos referir à saída deste processo com a sintaxe _.out_. E no momento podemos chamar isso do que quisermos, é apenas um nome de variável arbitrário.

## 2.1.2. Adicione um bloco output: ao script

Onde se torna importante é quando fazemos nosso novo bloco aqui, e este está abaixo do bloco de fluxo de trabalho agora, não estamos mais dentro do workflow. Chaves novamente. E é aqui que apenas dizemos ao Nextflow onde colocar todos os arquivos, que são criados pelo fluxo de trabalho.

Agora vou pegar este nome de variável, que criei aqui, e vou colocar ali e colocar algumas chaves para isso. E vou dizer ao Nextflow para usar um path. Ops. Path, entre aspas. E vou usar ponto. Isso apenas diz ao Nextflow para colocar o arquivo na raiz do diretório de resultados. Então não há subdiretórios ou nada.

Vamos tentar executar nosso fluxo de trabalho novamente. Se eu fizer _"nextflow run hello-world.nf"_, então esperançosamente deve parecer basicamente exatamente o mesmo. Nada realmente mudou com o Nextflow aqui. Está executando as mesmas coisas. Está apenas fazendo-as em diretórios work novamente.

Mas agora se eu fizer _"ls results/"_, você verá que há um novo diretório aqui que foi criado chamado results, que é o diretório base padrão para publicação de fluxo de trabalho. E lá dentro há um arquivo chamado _output.txt_.

Se eu fizer _"ls -l results"_, você verá que isso é na verdade um link simbólico para o diretório work. Então este não é um arquivo real, está linkado ao diretório work e coletou todos os arquivos ali para nós.

## 2.2. Defina uma localização personalizada

"Results" é o nome padrão para este caminho. Se eu executar o fluxo de trabalho novamente, e desta vez eu faço _traço_ hífen simples, isso é, porque é uma opção central do Nextflow. _" -output-dir **meus** resultados"._ Também poderia apenas fazer _"-o"_ de forma abreviada. Então vai definir um diretório base diferente para onde os arquivos são armazenados e mais uma vez, aqui em _myresults/_, agora temos um _output.txt_.

Isso é ótimo, mas provavelmente não queremos todos os arquivos apenas na raiz. Queremos alguma organização, então também podemos criar um subdiretório aqui chamado do que quisermos. Vamos dizer _"path 'hello_world'"_, e eu apenas executo isso novamente. _"nextflow run hello-world.nf"_. Deve ir para o diretório results em um subdiretório e com certeza, agora sob results aqui no topo temos _hello_world/_ e temos _output.txt_.

Coisa importante a notar, o antigo arquivo _output.txt_ ainda está lá. O diretório results não é limpo quando você faz isso. Apenas novos arquivos são copiados para lá. Eles vão sobrescrever arquivos que já estão lá se tiverem o mesmo nome de arquivo, mas não vão limpar os antigos. Então você precisa ser um pouco cuidadoso sobre quando reexecuta pipelines. Se você não quer que eles estejam em cima dos arquivos que já estão lá. Certifique-se de usar um diretório em branco vazio.

## 2.3. Defina o modo de publicação para copiar

Ok, mencionei que esses arquivos são links simbólicos, então se eu fizer _"ls -l results/hello_world/"_, você pode ver que está fazendo link simbólico para o diretório work. Isso é geralmente uma coisa boa se você está trabalhando em algo como HPC, e esses são arquivos realmente enormes e você não quer duplicá-los, porque significa que os arquivos são armazenados apenas uma vez no sistema de arquivos.

No entanto, isso significa que se você deletar o diretório work: se eu fizer _"rm -r work"_ e limpar todos aqueles arquivos intermediários que foram criados. Agora, se eu tentar ler este arquivo _"results/hello_world/"_. Vai estar apontando como um link simbólico para um arquivo que não existe mais e os dados se foram para sempre e são irrecuperáveis, o que talvez não seja ótimo.

Então geralmente nós, eu digo que é boa prática copiar os arquivos em vez de fazer link simbólico se você puder, porque é mais seguro. Apenas esteja ciente de que vai usar o dobro de espaço em disco a menos que você delete esses diretórios work.

Para fazer isso com o bloco de saída, vou para a primeira saída aqui. Defini o path antes e agora vou definir o mode e você pode ver conforme eu digito, a extensão VS code está, sugerindo coisas que sabe que é uma diretiva de saída aqui. E vou dizer copy. Eu aperto salvar.

Vamos reexecutar o fluxo de trabalho. Vai criar os arquivos novamente, novo diretório work.

Agora, se eu for para _"ls -l results/hello_world/"_ você pode ver que este é um arquivo real e não é mais um link simbólico, e o Nextflow copiou isso. Bom saber. Então path e mode são coisas que você vai se ver escrevendo bastante.

Agora, é claro, isso é muito simples. Vamos tornar isso mais complexo e poderoso conforme avançamos, e você verá como tornar essas coisas dinâmicas e não muito verbosas.

## 2.4. Nota sobre diretivas publishDir no nível de processo

Agora, eu disse quando começamos sobre isso, que esta é uma forma de sintaxe bastante nova. Está disponível apenas nas versões mais recentes do Nextflow enquanto gravo isso, e é chamada de Workflow Outputs.

Se você usar isso, é ótimo. Desbloqueia muitos outros recursos legais dentro do Nextflow, como, Nextflow Lineage para ajudar a rastrear a herança desses arquivos conforme são criados, e em breve será o padrão em 26.04. E em uma data posterior no futuro, esta será a única maneira de escrever seus fluxos de trabalho.

No entanto, como estamos nesta fase de transição agora, você pode muito bem ver pipelines em uso, que você usa algo chamado publishDir, que é a maneira antiga de fazer isso, e isso é definido não no nível de fluxo de trabalho e saída, mas é definido no nível de processo.

E esta declaração diz basicamente a mesma coisa. Diz, publique os arquivos de resultados em um diretório chamado results, e use um modo de cópia. Então você pode ver que a sintaxe é muito similar. Mas quando você está escrevendo novos pipelines agora, tente não usar essa diretiva publishDir, mesmo se você ver, em resultados de IA ou em documentação ou outros pipelines, porque essa é a maneira antiga de fazer isso.

Em 2026 todos devemos estar usando workflow outputs.

Isso está tudo documentado, se você está fazendo isso e já usou Nextflow antes, você pode ir para a documentação do Nextflow aqui, nextflow.io/docs/. E se eu rolar para baixo até tutoriais, há um tutorial chamado _Migrating to Workflow Outputs_.

É realmente bom. Passa por toda a sintaxe, como é equivalente à sintaxe antiga, por que mudamos, e, tem uma linha do tempo e tudo. E passa por todos os diferentes cenários com montes e montes de exemplos. Então você pode facilmente converter código Nextflow existente para a nova sintaxe.

## 3.1. Altere o processo sayHello para esperar uma entrada variável

Ok, então temos nosso script simples, que está executando um processo, criando um arquivo, dizendo ao Nextflow que é uma saída, e então estamos dizendo ao Nextflow onde salvar esse arquivo. Isso é um bom começo.

Mas seria mais interessante se não estivesse tudo codificado. Então a seguir, vamos pensar sobre como dizer ao Nextflow que este processo pode receber uma entrada variável, que é algo que podemos controlar em tempo de execução quando lançamos um fluxo de trabalho.

Precisamos fazer algumas coisas diferentes para fazer isso acontecer.

Primeiramente, precisamos dizer a este processo que ele pode aceitar uma variável de entrada e digitamos _input_ aqui como um novo bloco de declaração. E vamos chamar isso de _"val greeting"_.

A parte val é o equivalente de um path aqui embaixo. Diz ao Nextflow que isso é uma variável, como uma string neste caso. E se você passar o mouse sobre ela novamente, diz a você da extensão o que isso significa.

A seguir vamos dizer ao Nextflow o que fazer com isso. Não é, suficiente apenas dizer que há uma variável. Você tem que dizer no script como usar essa variável. E então vou me livrar desta string codificada aqui, e vou colocar uma variável.

Vou rapidamente fazê-lo sem chaves apenas para mostrar que isso é, permitido, e este é o estilo antigo de fazê-lo. Mas agora com a nova sintaxe, realmente recomendamos colocá-la dentro de chaves assim, e deixa realmente claro que isso está sendo interpolado pelo Nextflow aqui.

Ótimo. Então _"input greeting"_ vai para _$\{greeting\}._ Última coisa é precisamos dizer ao Nextflow no nível de fluxo de trabalho que este processo agora recebe uma entrada. E para fazer isso, vamos basicamente dar a ele uma variável.

## 3.2. Configure um parâmetro de linha de comando para capturar entrada do usuário

Poderíamos codificá-lo novamente, como Hello World, e isso funcionaria bem, mas obviamente não nos dá realmente nenhuma vantagem. Queríamos ser capazes de configurar isso em tempo de execução, então queremos ser capazes de fazê-lo na CLI, quando você lança o Nextflow.

E a maneira que fazemos isso é um conceito especial do Nextflow chamado _params_. Vamos chamar isso de _params.input_.

O que isso faz é expor essa variável de entrada na CLI e é aí que usamos um traço duplo quando lançamos o Nextflow.

Posso chamar isso do que eu quiser, posso chamar de _hello, greeting_. Não importa. O que quer que eu faça ali será exposto como uma opção de CLI quando lançamos um pipeline. E este é um verdadeiro truque mágico pelo Nextflow porque significa que você pode construir seu script de fluxo de trabalho muito rapidamente com esses parâmetros, e você está essencialmente construindo uma CLI customizada para seu pipeline, tornando realmente fácil customizar diferentes opções de forma dinâmica quando você lança.

Então. Vamos tentar. Voltar ao nosso terminal. Temos nosso comando _"nextflow run"_ aqui. E agora vou fazer _"--input"_, que corresponde ao _"params.input"_ que vimos antes. Acho que na documentação está em francês. Geraldine gosta de falar francês. Vou fazer em sueco porque vivo na Suécia. então vou dizer, "_Hej Världen_" e apertar enter.

Pode usar aspas simples ou duplas, apenas afeta como o Bash interpreta.

Executa o pipeline Nextflow exatamente da mesma maneira. Você pode ver o diretório de trabalho e tudo é o mesmo. Mas agora se eu subir para _"results/hello_world/output"_. Podemos ver nosso bom sueco aqui.

Então passamos dinamicamente uma entrada de uma CLI para um parâmetro. Passamos isso como uma entrada para o processo e o processo interpretou isso e colocou em um bloco de script, que então mudou dinamicamente a saída daquele resultado de script. Muito legal.

Lógica bastante complexa com muito, pouca sintaxe aqui. E você pode espero que veja como isso agora começa a escalar. E é assim que realmente construímos a lógica e a customizabilidade dos nossos pipelines no script Nextflow.

## 3.4. Use valores padrão para parâmetros de linha de comando

Ok, isso é ótimo. O problema agora é, toda vez que executo este pipeline, preciso fazer traço, input para que ele execute.

Se eu tentar executar sem este parâmetro, agora o Nextflow vai lançar um erro dizendo que precisava deste parâmetro e não foi definido. e então não sabia o que fazer.

Isso é uma coisa nova legal, a propósito. No passado, o Nextflow teria apenas executado com uma string vazia, e você teria tido todo tipo de erros estranhos, que teriam sido difíceis de entender. Mas no novo analisador sintático do Nextflow, é um pouco mais cuidadoso e diz logo de cara.

Então não queremos sempre especificar todas as opções. É boa prática especificar padrões sensatos. Então como fazemos isso no nosso script?

Você notará que quando escrevemos isso, apenas colocamos _params.input_ direto onde estamos usando. Então a solução óbvia é definirmos um padrão, e fazemos isso no topo do script aqui em um bloco params especial no fluxo de trabalho. Isto está no script do fluxo de trabalho aqui.

Novamente, alguma sintaxe nova aqui, então preste atenção. Isso é realmente coisa legal. Temos o nome do parâmetro, que será esperado aqui.

E então depois deste caractere de dois pontos, estamos definindo um tipo da variável. Você não precisa fazer isso, você pode apenas deixar em branco, mas é realmente bom. Diz ao Nextflow que estamos esperando uma string e tratá-la como tal.

Se quiséssemos um número em vez disso, por exemplo, poderíamos escrever float, e isso diria que queremos um número de ponto flutuante. E se tentarmos executar com isso, então vai lançar um erro. Se dermos uma string, que não é um float. E também vai passá-la como tal. Como se fizermos string, então sabe que é uma string. E mesmo se tiver zeros à esquerda e for todo numérico, ainda vai passá-la como uma string real.

Então essa segurança de tipo é um recurso muito novo do Nextflow, mas realmente poderoso para tornar seu código mais seguro de escrever e executar.

Então depois disso temos um símbolo de igual e então o valor padrão aqui. O Nextflow foi escrito em Barcelona originalmente, então parece apropriado que tenhamos algum, espanhol aqui, _"Holà mundo!"_ como padrão.

Certo, vou salvar esse script, voltar, executar o script novamente sem _--input_. E desta vez deve executar e vai criar nosso novo arquivo em _results_. E neste arquivo agora diz _"Holà mundo!"_.

Este é apenas um padrão no entanto, então não significa que não possamos ainda fazer a mesma coisa de antes. Se eu voltar e encontrar meu script antigo aqui, _"Hej Världen"_, porque faço _--input_ na linha de comando, isso vai sobrescrever aquele padrão e usar aquilo novamente no arquivo output.txt.

Então isso no script é apenas o valor padrão que estou definindo.

Conforme construímos nosso fluxo de trabalho para ser mais complexo e incluir mais parâmetros, este bloco params no topo do script começará a coletá-los todos em um lugar.

E você acaba com essa simetria bastante agradável no seu script, onde você efetivamente tem todas as suas entradas de fluxo de trabalho aqui e suas saídas de fluxo de trabalho na parte inferior. E está muito claro qual é a interface do seu fluxo de trabalho para o mundo externo. Então você pode pegar um novo pipeline muito rapidamente com a nova sintaxe e entender como usá-lo.

Uma última coisa legal. Não temos que definir um valor padrão com isso. Se fizermos params input mas não definirmos um valor padrão, então diz ao Nextflow que este parâmetro é obrigatório, e novamente, o pipeline falhará ao executar sem ele, mas vai dar uma mensagem de erro mais útil em vez de algo sobre ser nulo.

Então diz que estamos esperando que sua entrada seja obrigatória, mas não foi especificada na linha de comando. Muito bom.

Ok, então esperançosamente agora está claro sobre como configurar seu pipeline Nextflow com entradas e parâmetros variáveis, como definir o padrão, definir, os tipos, poderia ser um flag Booleano verdadeiro falso ou um inteiro ou tipos diferentes aqui. Como passá-los para seu fluxo de trabalho, onde vai através, e então interpola em seu processo. E então você também sabe como customizar esses na linha de comando quando lança o Nextflow. Isso está começando a parecer mais interessante do que nosso simples comando bash.

## 4. Gerencie execuções de fluxo de trabalho

Ok. O que vem a seguir? Para a parte final deste capítulo, vamos falar um pouco sobre como gerenciar todas as diferentes execuções de fluxo de trabalho. Se você olhar na minha barra lateral aqui e no Explorer embaixo de work, você verá que executei um monte de pipelines diferentes e esses diretórios work estão ficando bastante longos, há muitos deles.

E a outra coisa é, como eu disse antes, toda vez que reexecuto este pipeline, está criando um novo conjunto de diretórios work, e está reexecutando todos os processos do zero, o que é uma coisa boa. Esse é o comportamento pretendido. É reproduzível e está regenerando tudo fresco. Mas obviamente, se você está executando processos de longa duração, é chato sempre ter que começar seu pipeline do início se ele caiu no meio, ou se você muda algo no final do pipeline.

## 4.1. Relance um fluxo de trabalho com -resume

Felizmente, o Nextflow é realmente, bom em saber o que foi previamente executado e o que está disponível, e reusar aqueles resultados antigos é muito, simples. Apenas adicionamos uma, nova flag no final do comando _"-resume"_.

Agora, note que há dois hífens em input porque é o parâmetro. Há apenas um hífen em resume porque é uma opção central do Nextflow.

Confunde as pessoas o tempo todo, mesmo se você tem usado Nextflow por muito tempo. Então sempre lembre um ou dois hífens. Depende se é uma opção central do Nextflow.

Ok, então agora faço _-resume_ e executo exatamente o mesmo fluxo de trabalho novamente. E desta vez deve parecer praticamente exatamente o mesmo com uma diferença chave.

Na saída aqui, você pode ver que os resultados foram cacheados. E na verdade, este hash de tarefa aqui é exatamente o mesmo da execução anterior, e apenas reutilizou aquele diretório work em sua totalidade. As entradas e as saídas e o script foram todos não modificados. E então apenas pega aquele arquivo daquilo e se há etapas posteriores no processo, passaria eles para a próxima etapa no pipeline.

Então ainda está executando o pipeline inteiro do início ao fim, mas está usando resultados cacheados para cada uma dessas tarefas, onde pode.

Agora, quando você faz _-resume_, apenas retoma a última execução de pipeline no seu diretório de trabalho, qualquer que fosse. Mas você pode realmente retomar de qualquer execução anterior que você fez lá. E fizemos bastante agora.

## 4.2. Inspecione o log de execuções passadas

Para olhar todas elas, podemos fazer _"nextflow log"_ em vez de _"nextflow run"_, e isso vai nos dar uma saída bacana mostrando todas essas diferentes.. Preciso fazer minha tela um pouco menor para que possamos ver, todas essas diferentes execuções quando fizemos, o id da sessão, o comando e tudo.

E podemos olhar aqui e podemos pegar o nome de execução de qualquer uma dessas e então retomar uma daquelas específicas. Então posso voltar e posso retomar aquela chamada _hungry_ekeblad_. E apenas coloco isso depois do _resume_.

Se você está curioso, a propósito, todos esses adjetivos e nomes de cientistas estão no código-fonte do Nextflow. É um jeito realmente bom de obter seu primeiro pull request para o Nextflow indo e encontrando e adicionando seu cientista favorito.

E de qualquer forma, então fiz isso e voltou e olhou os resultados cacheados desta execução de fluxo de trabalho, percebeu que ainda poderia reutilizá-los, e fez. Então obtive os resultados cacheados novamente.

## 4.3. Delete diretórios work mais antigos

Isso é ótimo. E se eu quiser limpar esses diretórios work? Há montes deles aqui. Há montes de arquivos. Talvez eu saiba com certeza que quero retomar das últimas execuções de pipeline, mas não me importo com todas as anteriores a isso.

Então posso escolher uma aqui e posso usar outro comando Nextflow, que é _"nextflow clean"_, e posso fazer _"nextflow clean"_, vou fazer _"-before"_, e o nome de execução particular, que neste caso era _reverent_pike_ e vou fazer _"-n"_, que diz ao Nextflow apenas para fazer uma execução seca. Então apenas me diz o que vai deletar. Sem realmente fazer nada, então removeria esses diretórios work.

Isso parece sensato. Então vou fazer o mesmo comando novamente, mas em vez de _"-n"_ vou fazer _"-f"_ para realmente fazer a limpeza. E desta vez realmente removeu todos esses diretórios. E se eu entrar e olhar os diretórios work, agora está parecendo muito mais leve. Fantástico.

Então é assim que limpar todos os seus diretórios work locais de uma maneira bastante segura sem completamente destruir o cache. Então você ainda pode retomar se quiser.

Se você esquecer o que essas flags são para cada comando Nextflow você pode fazer _"nextflow help"_, e então o nome do comando. Então se eu fizer _"nextflow help clean"_, você pode ver todas as diferentes opções: _-after, -before, -but_, todas diferentes maneiras de configurar esse comportamento de limpeza. Muito legal.

## Conclusão

Ok, esse é o fim da parte um do Hello Nextflow. É um começo bastante intenso para o curso, mas esperançosamente agora você tem um entendimento bastante bom de como um script Nextflow se parece; com diferentes partes chave, os processos, os fluxos de trabalho, as saídas e os parâmetros. Você sabe como configurá-los com sobrescritas básicas da linha de comando, como fazer um bloco de entrada dinâmica com um script dinâmico e você sabe como gerenciar todas as suas execuções de carga de trabalho: vendo o que você já executou, retomando, limpando. Há um monte de coisas. Você percorreu um longo caminho. Então se você quiser fazer uma pausa e dar uma caminhada rápida e uma xícara de chá, agora é provavelmente um bom momento. Você mereceu.

Daqui em diante, estamos basicamente construindo sobre essa fundação. Como podemos tornar isso mais complexo, mais poderoso? Como podemos torná-lo mais flexível? Fazer as coisas que queremos fazer nossa análise em escala.

## Quiz

Agora se você rolar para baixo até a parte um, hello world, na página web você verá um pequeno quiz e isso é algo novo que fizemos para esta versão do treinamento Nextflow. E você pode passar e testar a si mesmo para verificar que você entendeu todo o material que fizemos neste capítulo.

Isso não é enviado para nós ou nada, é apenas armazenado no seu navegador. Então não sabemos quais são suas respostas, mas é apenas uma pequena verificação para ter certeza de que você não perdeu nada ou entendeu mal alguma coisa. E você pode tentar quantas vezes quiser.

Se você é como eu, talvez queira ficar no terminal na sua instância de VS Code, nesse caso você pode digitar o comando _quiz_ e então apenas dizer em qual, capítulo você está. Então fazemos _"Hello World"_, e então você pode fazer exatamente o mesmo, perguntas de quiz, que estão no navegador web, mas apenas no seu terminal.

Legal. Ok. Espero que goste disso. Divirta-se um pouco e, nos vemos no próximo capítulo em apenas um minuto para falar tudo sobre canais Nextflow.
