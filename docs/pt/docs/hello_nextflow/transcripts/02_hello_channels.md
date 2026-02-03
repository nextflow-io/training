# Parte 2: Hello Channels - Transcrição

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Tradução assistida por IA - [saiba mais e sugira melhorias](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

<div class="video-wrapper">
  <iframe width="560" height="315" src="https://www.youtube.com/embed/lJ41WMMm44M?si=xCItHLiOQWqoqBB9&amp;list=PLPZ8WHdZGxmXiHf8B26oB_fTfoKQdhlik" title="YouTube video player" frameborder="0" allow="accelerometer; autoplay; clipboard-write; encrypted-media; gyroscope; picture-in-picture; web-share" referrerpolicy="strict-origin-when-cross-origin" allowfullscreen></iframe>
</div>

!!!note "Notas importantes"

    Esta página mostra apenas a transcrição. Para instruções passo a passo completas, retorne ao [material do curso](../02_hello_channels.md).

    Os números de seção mostrados na transcrição são fornecidos apenas para fins indicativos e podem não incluir todos os números de seção dos materiais.

## Boas-vindas

Olá, bem-vindo à parte dois do Hello Nextflow.

Este capítulo é chamado Hello Channels. Vamos falar tudo sobre esta parte fundamental do Nextflow.

Canais são as coisas que conectam as diferentes etapas do seu pipeline, a maneira como seus dados e lógica fluem através do seu fluxo de trabalho.

Ok, vamos começar.

Vamos começar indo para training.nextflow.io

Hello Nextflow na barra lateral e clicando na parte dois. Hello Channels.

Todo o material está escrito aqui, então você pode seguir no seu próprio ritmo e revisar qualquer coisa que possa ter perdido.

Depois de abrir o site, você pode carregar o Codespaces e continuaremos de onde paramos no final do último capítulo.

## 0. Aquecimento: Execute hello-channels.nf

Para este capítulo, vamos editar um arquivo diferente. Este se chama Hello Channels, então você pode encontrá-lo na barra lateral, clique duas vezes para abrir.

Agora, se você acabou de vir do capítulo um, este arquivo parecerá muito familiar. O ponto de partida aqui é basicamente onde terminamos o capítulo um, com nosso processo chamado sayHello, nossa entrada, saída, nosso publishDir e nosso params.greeting, e nosso fluxo de trabalho simples.

Estamos começando com um novo arquivo, então é um terreno nivelado para todos, mas você pode continuar com seu arquivo anterior se preferir.

Observe, eu também deletei todos os arquivos .nextflow\* e os diretórios work aqui, apenas para que seja um ponto de partida limpo. Não importa se você faz isso ou não, fica a seu critério.

Ok. Vamos começar verificando se este pipeline ainda funciona como esperamos. Vou abrir o terminal aqui.

Faço "nextflow run hello-channels.nf" e aperto enter.

Vai executar aquele pequeno fluxo de trabalho, executa nossa etapa sayHello, gera um diretório work com aquele hash, e aqui está nossa pasta results e lá está nosso arquivo de saída, exatamente como esperávamos do nosso params.greeting padrão.

Então isso é ótimo. Exatamente o mesmo que o capítulo um, funcionando como esperávamos.

## 1. Forneça entradas variáveis através de um canal explicitamente

No capítulo um, você na verdade já estava usando canais, você só não percebeu. Quando especificamos uma string aqui, o Nextflow automaticamente criou um canal em torno dessa string para nós, apenas porque sabia que estávamos chamando um processo, então precisávamos de um canal de entrada.

A primeira coisa que vamos fazer é torná-lo explícito, escrevendo realmente o canal em si.

## 1.1. Crie um canal de entrada

Então vou ir para o fluxo de trabalho aqui no final do script, e vou dizer greeting_ch. Esta é uma convenção que frequentemente usamos no código Nextflow de ter um underscore ch no final de um nome de variável quando é um canal, apenas para que seja fácil identificar que é um canal, mas você não precisa fazer isso. Igual a channel of Hello Channels.

O que acabamos de usar é algo chamado "Channel Factory" na linguagem do Nextflow. É essa coisa aqui, estamos definindo esta variável para um novo canal, e esta fábrica de canais aqui está criando um canal para nós de uma maneira particular.

Existem um punhado de diferentes fábricas de canais que o Nextflow tem, para criar canais de diferentes tipos de entradas. Dot of é a mais simples, e apenas recebe quaisquer strings que fornecemos.

Observe que quando passo o mouse sobre essas palavras no VS Code, a extensão Nextflow está me dando um popup explicando o que esta sintaxe faz, e também há um texto de ler mais na parte inferior dessa janela popup.

Se eu clicar nisso, vai abrir os documentos do Nextflow. Em uma nova aba e me levar direto para a documentação deste item específico. Neste caso para channel.of.

## 1.2. Adicione o canal como entrada para a chamada do processo

Observe que a extensão também está nos dando um aviso, dizendo que criamos um novo canal aqui, mas ele não está sendo usado por nada.

Então, vamos corrigir isso. Vou pegar o novo nome do canal e vou substituir este params.greeting pelo nosso novo canal.

Observe que não estamos mais usando o sinalizador de linha de comando --greeting agora, params.greeting não está sendo usado, estamos voltando a codificar esta string diretamente. Tudo bem. Estou apenas tentando manter as coisas simples. Voltaremos mais tarde e usaremos os params novamente.

## 1.3. Execute o comando workflow novamente

Ok, vamos apenas verificar se isso funciona. Abra o terminal e observe novamente. Nextflow run hello channels. Verifique output.txt, e lá está.

Ótimo, um exemplo meio chato, fazendo exatamente a mesma coisa que fizemos antes, mas agora pelo menos a lógica está um pouco mais clara. Estamos sendo explícitos sobre escrever um novo canal.

Efetivamente acabamos de escrever mais código para fazer a mesma coisa. Mas isso começará a fazer mais sentido conforme nos tornarmos um pouco mais complicados com a forma como criamos nossos canais.

## 2. Modifique o fluxo de trabalho para executar em múltiplos valores de entrada

Ok, vamos tornar isso um pouco mais interessante. É muito raro que você queira executar um pipeline Nextflow em uma única entrada, então vamos dar várias entradas a ele.

## 2.1. Carregue múltiplas saudações no canal de entrada

Da documentação aqui. Vou copiar essas diferentes strings, três delas. Hello, Bonjour, Olà. Oh, espere. O Copilot está sugerindo algumas outras. Então vamos apertar tab e digitar essas.

A documentação do Nextflow aqui nos diz que podemos dar múltiplos valores para este operador, então deve funcionar, mas vamos tentar e ver o que acontece.

## 2.1.2. Execute o comando e observe a saída do log

Bem. Sim e não. Vamos ver. Diz que cinco de cinco tarefas foram executadas aqui, mas mostra apenas um hash, o que é um pouco estranho. Tudo bem. Tudo está como esperado aqui. Por padrão. O Nextflow usa um tipo especial de saída para um terminal chamado códigos de controle ANSI, o que significa que sobrescreve certas linhas para dar uma visão comprimida agradável de todos os diferentes processos que estão sendo executados.

Isso faz muito mais sentido quando você tem fluxos de trabalho maiores e está executando centenas ou milhares de diferentes amostras. Você pode simplesmente gerar tanta saída no terminal que é impossível de ver, enquanto esta visão de atualização lhe dá um progresso em tempo real para você.

## 2.1.3. Execute o comando novamente com a opção -ansi-log false

Se você quiser, pode executá-lo novamente, e desta vez vou usar um argumento central adicional do Nextflow com um único hífen dizendo, "-ansi-log false". Isso usa a versão anterior da saída de log do Nextflow. E aqui você pode ver todos os processos individuais que foram lançados.

Fica a seu critério fazer isso ou não. A saída do Nextflow é exatamente a mesma em ambos os casos.

## 2.2. Certifique-se de que os nomes dos arquivos de saída serão únicos

Ok, vamos dar uma olhada nos arquivos de saída, então iremos para results. Mas há apenas um único arquivo de saída. O que aconteceu? Vimos que o processo havia sido executado várias vezes. Podemos ir ao diretório work e ver todos os diferentes hashes, todas as tarefas foram executadas corretamente. Mas se você se lembrar no nosso processo aqui, estamos salvando tudo em um arquivo output.txt e depois publicando isso para este diretório.

Então o mesmo arquivo foi criado cinco vezes, e então foi sobrescrito cinco vezes. E temos apenas qualquer tarefa que tenha sido executada por último.

## 2.2.1. Construa um nome de arquivo de saída dinâmico

A maneira como corrigimos isso é usando um nome de arquivo de saída dinâmico. Aqui já temos uma variável chamada greeting dentro do processo, então podemos usar isso no nome do arquivo de saída. Copio isso e faço $greeting-output.txt.

Vou cercar isso com aspas, apenas para que o bash não fique confuso com quaisquer espaços que possam aparecer aqui. E então vou pegar o mesmo nome de arquivo e atualizar a saída aqui.

É realmente importante que a saída corresponda a isso, porque caso contrário, este arquivo não será encontrado e o Nextflow vai travar.

Vou fazer mais uma edição realmente importante, que é mudar essas aspas simples por aspas duplas. Observe que a cor do código mudou quando fiz isso. Esta variável só é expandida se usarmos aspas duplas. Se eu usar aspas simples aqui, é usado como um valor literal, e eu teria um único arquivo chamado $greeting-output, o que não é o que eu quero.

## 2.2.2. Execute o fluxo de trabalho

Então vamos colocar as aspas duplas de volta e tentar.

Vou apenas limpar meu diretório antes de começar, para que seja fácil ver os novos arquivos. Vou deletar qualquer coisa chamada .nextflow, work e results.

E vou executar aquele comando Nextflow novamente e vamos ver quais arquivos são criados. Então ele executa os cinco processos ali. Se você estivesse assistindo muito atentamente, pode ter visto aquela linha atualizar enquanto estava executando.

E agora podemos ir ao diretório results, e com certeza, temos cinco saídas diferentes, e todas elas são prefixadas com a saudação diferente.

Se eu abrir cada uma delas, veremos que cada uma contém a saudação correspondente. Fantástico. É isso que queremos.

## 3. Use um operador para transformar o conteúdo de um canal

Ok, então agora sabemos o que são canais e sabemos o que são fábricas de canais. E quanto aos operadores? Este é outro termo para parte da linguagem Nextflow, que é uma série de funções que nos permitem operar em canais para fazer certas coisas com eles. Nextflow, vem com um conjunto de operadores, que nos permitem manipular canais de várias maneiras diferentes.

## 3.1. Forneça um array de valores como entrada para o canal

Vamos trabalhar nisso com um exemplo. Digamos que queremos pegar essas strings de entrada, mas em vez de apenas colocá-las diretamente em uma fábrica de canais, queremos defini-las como um array.

## 3.1.1. Configure a variável de entrada

Então vou pegar isso e fazer isso como uma nova linha acima e dizer, greetings, array.

Pronto. Vou pegar aquela variável array e colocá-la no channel.of, e salvar.

## 3.1.3. Execute o fluxo de trabalho

Agora, vamos ver o que acontece. Volto ao meu terminal. Vou apenas limpar todos aqueles arquivos temporários novamente. E vamos executar o fluxo de trabalho.

Nada bom. Ok. Quebrou. Tudo bem. Eu estava esperando que quebrasse desta vez. Depurar o que dá errado quando um fluxo de trabalho Nextflow falha é uma parte fundamental de ser um desenvolvedor Nextflow. Isso vai acontecer muito e é importante entender o que a mensagem de erro diz e como lidar com isso.

As mensagens de erro do Nextflow são na verdade bastante estruturadas. Nos diz qual processo deu errado. Nos dá uma mensagem de erro por uma razão. Diz qual foi o comando que tentou executar dentro daquela tarefa particular, qual foi o status de saída, qual foi a saída e onde estava o diretório work daquela tarefa.

Observe que posso clicar com option nisto no VS Code e ele abre em uma barra lateral para que eu possa ir direto lá e visualizar todos esses arquivos ocultos, sobre os quais falamos no capítulo anterior, incluindo o arquivo .command.sh. Você pode ver que este é o mesmo que os comandos que foi executado aqui.

Ao olhar para este arquivo, podemos ter uma ideia do que pode ter dado errado aqui em vez de executar uma única tarefa para cada elemento no array como fez da última vez, ele apenas forneceu o array inteiro de uma vez como uma string. Então precisamos desempacotar aquele array em valores individuais antes de passá-lo para o canal. Vamos voltar e ver se podemos fazer isso usando um operador.

## 3.2. Use um operador para transformar o conteúdo do canal

Neste caso, não vamos alterar o array antes de passá-lo para o canal. Vamos ajustar o canal para que ele se comporte da maneira que esperamos. Vamos fazer isso usando o operador flatten pode fazer dot começar a digitar e podemos ver que a extensão do VS Code começa a sugerir todos os diferentes operadores que temos disponíveis.

## 3.2.1. Adicione o operador flatten()

E vou selecionar flatten. Observe que o espaço em branco não importa neste contexto para o Nextflow. Então você pode colocar esses operadores em uma nova linha se quiser. Então posso soltar isso aqui e indentar para que fique embaixo de ".of" e você verá que as pessoas frequentemente encadeiam muitos operadores assim em um canal e indentam desta forma para que seja mais fácil de ler.

Você também pode ver, como antes, posso passar o mouse sobre isso e ler o que o operador flatten está fazendo, e também seguir um link para a documentação se eu quiser.

Então este operador está pegando este canal, que tem um único array dentro dele, e separando os valores do array.

## 3.2.2. Adicione view() para inspecionar o conteúdo do canal

Podemos espiar dentro dos canais usando o operador especial view, e vou adicionar alguns deles aqui. Isso é um pouco como usar instruções print em outras linguagens. Então vou fazer dot view e então vou usar esses colchetes ondulados.

Isso é chamado de closure. Isso basicamente fornece código adicional para o operador view, que ele executará em cada item dentro do canal. Neste caso, vou dizer greeting before flatten. Greeting.

Estou definindo uma variável aqui, que está apenas dentro do escopo deste closure. Então esta variável é usada apenas aqui e eu poderia chamá-la como quisesse. Não importa muito. Estou apenas usando greeting para tornar mais fácil de ler.

Em alguns pipelines Nextflow, você pode ver pessoas usando uma variável implícita especial chamada "$it". Assim. Esta é uma variável especial dentro do código Nextflow, que é um atalho para que você não precise fazer a pequena definição de uma variável. No entanto, ao longo do tempo estamos pensando, isso não é muito claro para pessoas que são novas no Nextflow, e desencorajamos o uso de "$it" agora.

Então vou ficar com o comportamento anterior de greeting e usá-lo assim porque isso é mais explícito e é mais claro sobre o que está acontecendo.

Vou então copiar esta linha e fazer exatamente a mesma coisa novamente após os argumentos flatten. O operador view é um pouco especial porque faz algo nos elementos, mas também apenas continua passando-os para o próximo operador, então podemos encadeá-lo no meio de uma cadeia de operações assim, e ele imprimirá o status lá e continuará. Então, esperançosamente, isso nos mostrará como o canal se parece antes e depois do operador flatten.

## 3.2.3. Execute o fluxo de trabalho

Vamos tentar. Limpar. Limpar tudo no espaço de trabalho. Execute o pipeline novamente.

Ok, então podemos ver que executou nossos cinco processos. Novamente, não travou com um erro, então isso é definitivamente bom. E agora temos o before flatten e com certeza temos nosso array e temos after flatten, impresso cinco vezes uma vez para cada elemento do array. Isso é exatamente o que esperávamos. Então isso é realmente uma boa notícia. E isso se encaixa exatamente com o que esperaríamos do código.

Não precisamos mais dessas instruções de depuração, então posso comentá-las ou deletá-las. Vou deletá-las apenas para manter meu código limpo e organizado. Ok, ótimo. Este exemplo agora está funcionando bem e podemos começar a ver como os canais podem fazer uma lógica um pouco mais complicada.

## 4. Use um operador para analisar valores de entrada de um arquivo CSV

Agora vamos tentar fazer isso usando um arquivo com uma série de entradas. Esta é uma maneira muito comum de escrever pipelines Nextflow usando uma planilha de amostras ou um CSV de metadados.

## 4.1. Modifique o script para esperar um arquivo CSV como fonte de saudações

Se eu for à barra lateral, você pode ver greetings.csv no repositório de exemplo, e este é um arquivo CSV muito, muito simples que apenas contém três linhas com três saudações diferentes. Vamos ver se podemos usar este arquivo CSV dentro do nosso fluxo de trabalho.

Agora vou voltar a usar params como fizemos no capítulo um, para que possamos ter uma entrada de linha de comando.

Vou deletar este array greetings.

## 4.1.1. Altere o parâmetro de entrada para apontar para o arquivo CSV

Vou definir params greeting para o nome do arquivo, que é greetings.csv, e vou usar esta variável especial para gerar o canal. Vou colocar isso aí, e os erros desaparecem. Lembre-se de que isso está definindo esta variável por padrão agora. Então, se eu executar o pipeline sem nenhum argumento, ele usará greetings.csv, mas eu poderia fazer --greeting para sobrescrever esta variável se quisesse.

## 4.1.2. Mude para uma fábrica de canais projetada para lidar com um arquivo

Ok, estamos passando um arquivo agora em vez de uma string ou um array de strings, então provavelmente precisamos de uma fábrica de canais diferente.

Vamos nos livrar de "of" que estamos usando até agora, e em vez disso usar .fromPath. Isso faz exatamente o que parece. Cria um canal com caminhos em vez de valores, usando um nome de arquivo ou glob string. Também vou remover o operador flatten, pois não precisamos mais disso, agora que estamos passando um arquivo.

## 4.1.3. Execute o fluxo de trabalho

Vou salvar, abrir o terminal, executar o fluxo de trabalho e então ver o que acontece.

Ok. Travou novamente. Não se preocupe. Eu estava esperando este também. Vamos dar uma olhada na mensagem de erro e ver se conseguimos descobrir o que está dando errado. Aqui podemos ver o comando executado, e um pouco como antes onde tínhamos todo o array impresso. Agora temos o caminho do arquivo sendo ecoado no comando, em vez de passar pelo conteúdo do arquivo.

## 4.2. Use o operador splitCsv() para analisar o arquivo

Então, para usar o conteúdo do arquivo em vez disso, precisamos de outro operador. O operador que vamos usar para este é chamado splitCsv. Faz sentido, porque é um arquivo CSV que estamos carregando.

## 4.2.1. Aplique splitCsv() ao canal

Ok, então splitCsv. Fecha parêntese. Não precisamos de nenhum argumento aqui. E novamente, vou usar alguns operadores view para dar uma visão do que está acontecendo aqui.

.view csv after splitCsv. Before split Csv.s

## 4.2.2. Execute o fluxo de trabalho novamente

Ok, vamos tentar executar isso e ver o que acontece.

Ok, temos um pouco mais de saída desta vez, mas ainda falhou. Podemos olhar as instruções view, e aqui você pode ver before split CSV, e temos um caminho de arquivo como vimos na mensagem de erro anterior. After split CSV, agora temos três valores correspondentes às três linhas no arquivo CSV.

No entanto, você pode ver que cada um desses valores está cercado por colchetes. Então cada um deles era um array por si só, e isso nos deu a mesma área que tínhamos antes, onde está tentando ecoar um array em vez de apenas uma única string.

Se pensarmos sobre um arquivo CSV, isso meio que faz sentido. Normalmente, um arquivo CSV terá linhas e colunas, então split CSV faz array bidimensional. A primeira dimensão do array é cada linha, e então há uma segunda dimensão, que é cada coluna para cada linha.

Então aqui temos apenas um único valor em cada linha, então temos uma única coluna, então temos um array de um elemento para cada linha do arquivo.

Tudo bem. Precisamos apenas de outro operador para colapsar aquele array para cada linha do arquivo CSV analisado. Vamos limpar isso. Livrar-se de um terminal e ver o que podemos fazer.

## 4.3. Use o operador map() para extrair as saudações

Agora poderíamos usar o operador flatten novamente, que usamos antes. Vimos como ele pode colapsar um array em uma série de valores, o que funcionaria muito bem aqui. Mas vou usar a oportunidade para demonstrar outro operador, que é muito comum dentro de fluxos de trabalho chamado operador map.

## 4.3.1. Aplique map() ao canal

Vou fazer dot map e vou fazer item item[0].

Se você escreve muito código em outras linguagens, pode já estar familiarizado com o operador map. Ele pega um iterável, como um array ou um canal, e faz alguma operação em cada valor dele.

Aqui estamos dizendo que devemos definir uma variável chamada item dentro do escopo deste closure, e então queremos retornar, apenas o primeiro valor naquele array. Então item índice zero.

Isso está efetivamente achatando o array. Você pode ver como poderíamos estender isso para ser mais complexo, no entanto: se nosso arquivo CSV tivesse seis colunas, mas estivéssemos apenas interessados na quarta coluna, poderíamos acessar um índice específico aqui. Ou fazer qualquer outro tipo de operação no valor antes de passá-lo para processamento downstream.

Então o operador map é extremamente flexível e muito poderoso para modificar canais em tempo de execução. Vamos colocar outra instrução view apenas para que possamos ver o que está fazendo em nossa execução. Pode adjudicar aquela linha e movê-la para baixo. E after map.

## 4.3.2. Execute o fluxo de trabalho mais uma vez

Vamos abrir o terminal e tentar executar o fluxo de trabalho.

Ok, sem erros desta vez. Isso é um bom sinal. Agora podemos passar por todas essas diferentes saídas das instruções view. Before split CSV, tínhamos um único caminho. After split CSV, tínhamos os arrays de valor único, e então after map, temos apenas os valores sem nenhuma sintaxe de array. Vamos até o diretório results, e aqui estão nossos arquivos de saída se comportando exatamente como queríamos.

Há um pequeno bônus aqui. Você pode realmente ver que os operadores view estão ligeiramente misturados na ordem em que fizeram a saída. Isso é porque o Nextflow está fazendo paralelização dessas diferentes tarefas. Então, depois que dividiu o CSV, há três elementos neste canal, e está lidando com o processamento desses três elementos em paralelo automaticamente. Isso significa que a ordem das saídas é estocástica e pode variar. Neste caso, apenas aconteceu que alguns dos operadores view retornaram depois que a etapa subsequente havia sido concluída, e então veio nesta ordem.

Se eu executar o mesmo fluxo de trabalho novamente. Então com certeza, veio em uma ordem diferente e desta vez temos os split CSVs e os maps na ordem que esperaríamos.

Então apenas tenha em mente, você não pode confiar na ordem das saídas de uma tarefa de processo porque o Nextflow está lidando com essa paralelização para você automaticamente. O Nextflow faz isso para você com sua lógica de fluxo de dados, e esse é o verdadeiro poder do Nextflow.

Ok, este é provavelmente um dos capítulos mais importantes de todo o treinamento. Uma vez que você entende canais, fábricas de canais e operadores, você começa a entrar na força do Nextflow e no que o torna único como linguagem de programação. Esta funcionalidade permite ao Nextflow paralelizar todos os seus fluxos de trabalho para você e gerar lógica de fluxo de trabalho extremamente complexa com uma sintaxe muito limpa e um modelo de fluxo de dados push. Pode ser um conceito um pouco estranho no início, mas uma vez que você se acostuma a escrever código assim, rapidamente se sentirá natural e antes que você perceba, estará escrevendo fluxos de trabalho fantásticos.

Faça uma pausa, uma xícara de chá, dê uma volta e vamos para o capítulo três, onde começamos a estender esses conceitos para fluxos de trabalho mais complexos. Vejo você no próximo vídeo.

[Próxima transcrição de vídeo :octicons-arrow-right-24:](03_hello_workflow.md)
