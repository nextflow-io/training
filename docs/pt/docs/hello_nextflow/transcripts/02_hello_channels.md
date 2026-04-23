# Parte 2: Hello Channels - Transcrição do Vídeo

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Tradução assistida por IA - [saiba mais e sugira melhorias](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

<div class="video-wrapper">
  <iframe width="560" height="315" src="https://www.youtube.com/embed/yDR66fzAMOg?si=xCItHLiOQWqoqBB9&amp;list=PLPZ8WHdZGxmWKozQuzr27jyMGqp9kElVK&amp;cc_load_policy=1" title="YouTube video player" frameborder="0" allow="accelerometer; autoplay; clipboard-write; encrypted-media; gyroscope; picture-in-picture; web-share" referrerpolicy="strict-origin-when-cross-origin" allowfullscreen></iframe>
</div>

!!!note "Notas importantes"

    Esta página mostra apenas a transcrição. Para instruções passo a passo completas, retorne ao [material do curso](../02_hello_channels.md).

    Os números de seção mostrados na transcrição são fornecidos apenas para fins indicativos e podem não incluir todos os números de seção dos materiais.

## Boas-vindas

Olá e bem-vindo de volta à Parte 2 do Hello Nextflow. Este capítulo se chama Hello Channels.

Canais são como a cola no seu pipeline Nextflow. São as partes que mantêm todos os diferentes processos juntos, que o Nextflow usa para passar todas as informações e orquestrar seu fluxo de trabalho.

Há outra parte dos canais que são os operadores. Estes são basicamente funções que podemos usar em canais para modificar o conteúdo. Vamos mergulhar no VS Code e ver onde estamos.

Estou bem ampliado neste VS Code, então para manter as coisas limpas e organizadas, removi todos os arquivos _.nextflow\*_ e o diretório _work/_ e o _results/_ e tudo do Capítulo Um. E estou apenas começando do zero aqui. Mas não se preocupe muito com isso. Se você não quiser, pode deixar esses arquivos por aí. Eles não causarão problemas.

Vamos começar trabalhando em _hello-channels.nf_ para este capítulo, e se eu abrir isso, deve parecer muito similar ao arquivo em que estávamos trabalhando anteriormente. Pode ser que diferentes partes estejam em diferentes partes do script, mas tudo deve ser basicamente o mesmo.

Uma coisa que é diferente é que o caminho no bloco de saída aqui agora é _hello_channels_ para esta parte, o que significa que os arquivos de resultado serão armazenados em um subdiretório diferente em seus resultados se você ainda tiver isso lá. Então deve ser um lugar limpo e agradável para começar sem ficar confuso sobre as saídas.

Ok, então vamos lembrar rapidamente o que este script faz quando executamos este fluxo de trabalho. Fazemos _"nextflow run hello-channels.nf"_. Podemos fazer _"--input myinput"_, e quando executamos isso, ele vai usar este parâmetro, params.input, que foi passado como a variável para o processo sayHello aqui em cima, que vai para greeting e é salvo em output.txt. E podemos ver isso no arquivo de resultados. Ótimo.

## 1. Fornecer entradas variáveis via um canal explicitamente

Isso é legal. Mas é, é bem simplista. Temos uma variável neste parâmetro, que vai para um processo que executa uma vez, e não escala muito bem. E não podemos dar a ele muitos arquivos diferentes para criar aqui. Não podemos dar a ele muitas saudações diferentes. Temos apenas uma.

Na realidade, o Nextflow é tudo sobre escalar sua análise. Então você provavelmente quer que ele faça mais de uma coisa. E fazemos isso com _canais_.

Canais são um conceito um pouco único para muitas pessoas que estão começando com Nextflow. Vem desses conceitos de programação funcional, e pode levar um pouco de tempo para entender, mas uma vez que você entende, eles realmente desbloqueiam o poder do Nextflow e são fundamentais para como você escreve seus fluxos de trabalho.

## 1.1. Criar um canal de entrada

Vamos começar pegando este script e fazendo-o usar um _canal_ em vez de apenas um _parâmetro_.

Vamos para o fluxo de trabalho, que é onde toda a nossa lógica de fluxo de trabalho está sobre juntar as coisas. E vou entrar aqui e vou criar um novo canal.

Criar um novo canal.

E vou chamá-lo de "_greeting_ch"_. Esta é a convenção de fazer "_\_ch"_ assim, apenas para que você possa lembrar que esta variável é um canal. Mas você pode chamá-la do que quiser.

E então vou dizer igual, e vou fazer _"channel.of"._

Channel é como o namespace para tudo relacionado a canais. "c" minúsculo se você já usou Nextflow antes. E o _".of"_ é algo chamado de fábrica de canais, que é basicamente uma maneira de criar um canal.

Existem muitas fábricas de canais diferentes. Se eu fizer apenas "." aqui, você pode ver que o VS Code está sugerindo várias delas, mas _".of"_ é a mais simples e apenas recebe uma entrada aqui.

Então posso fazer alguns parênteses e vou dizer _"Hello Channels!"_.

Ótimo. Tenho um canal. Fantástico. Posso salvar, poderia executá-lo novamente, mas nada interessante vai acontecer. O VS Code me deu uma linha de aviso laranja aqui embaixo e me disse que isso está configurado: você criou isso, mas nunca realmente usou para nada. Este canal não está sendo consumido.

Ok, então como usamos isso? Muito simples. Vou pegar isso, copiar, e vou deletar _params.input_ e vou colocar _"greeting_ch"_ aqui em vez disso. Então vamos passar este canal como entrada para sayHello.

Note que eu codifiquei esta string por enquanto. Isso é um pouco de um passo para trás depois do nosso bom parâmetro que usamos no final do último capítulo, mas apenas mantém as coisas simples para começar para que você possa ver a lógica.

Ok, vou entrar no meu terminal e vou executar o fluxo de trabalho novamente. Sem nenhum _"--input"_ desta vez, e vai executar e vai usar aquele canal que criamos e esperançosamente devemos ter um arquivo aqui em _results/hello_channels/_ e agora diz "Hello Channels!". Fantástico. Então é isso que estamos esperando do nosso canal aqui. Ótimo.

## 1.4. Usar view() para inspecionar o conteúdo do canal

Mais uma coisa para adicionar aqui, apenas uma rápida introdução a outra função que podemos usar em canais chamada "_.view"_.

Isso é análogo ao comando _print_ em Python ou outras linguagens que você pode estar acostumado, e apenas despeja o conteúdo deste canal no terminal quando o executamos.

Então faça "_.view"_, e então se eu executar o fluxo de trabalho novamente, ele deve imprimir no terminal qual é o conteúdo daquele canal, no momento em que o criamos.

Com certeza, você pode ver que foi impresso no terminal aqui. _"Hello Channels!"_.

Note que você pode quebrar essas coisas em linhas se quiser, e de fato, o formatador automático do Nextflow tentará fazer isso para você. Espaço em branco não é realmente importante aqui, então você pode encadear essas coisas uma após a outra.

## 2. Modificar o fluxo de trabalho para executar em múltiplos valores de entrada

Ok, então nosso canal tem uma coisa que é legal, mas é basicamente o mesmo de antes. Então vamos torná-lo um pouco mais complicado. Vamos adicionar mais algumas coisas ao nosso canal.

A fábrica de canais "_.of()"_ pode receber múltiplos itens, então vamos escrever mais alguns. Faremos _Hello, Bonjour, Hej_. E então podemos executar este fluxo de trabalho novamente e veremos o que acontece.

Deve executar novamente. E imprimimos agora. _"Hello", "Bonjour"_ e _"Hej"_ no terminal com nossa instrução view. Fantástico.

## 2.1.2. Executar o comando e olhar a saída do log

Você pode pensar que terminamos neste ponto. Mas na verdade há uma pegadinha aqui, que vai nos atrapalhar. Se olharmos nosso arquivo de saída aqui. Você pode ver que tem _"Hello"_ dentro, mas não tem nenhuma das outras saídas. Na verdade, é apenas este.

Se executarmos este fluxo de trabalho várias vezes, podemos até ver que às vezes tem _"Bonjour"_, às vezes tem _"Hej"_. É um pouco aleatório.

Se olharmos o terminal, podemos ver que executou três vezes e podemos ver as diferentes saídas do view. Mas se eu for para o diretório work, posso fazer _"cat work"_. Colocar este hash e expandir isso e _output.txt_. Você pode ver que este arquivo no diretório work é diferente do diretório results, e este é _"Hej"._ Então há algo não funcionando direito aqui.

E a chave é que temos três tarefas que executaram. A saída do Nextflow tenta resumir isso conforme o processamento continua, para que não tome completamente todo o seu terminal, e esse log ANSI usa códigos de escape ANSI, basicamente sobrescreveu as outras tarefas. Então apenas mostra a última que aconteceu de ser atualizada.

## 2.1.3. Executar o comando novamente com a opção -ansi-log false

Há algumas coisas que podemos fazer para realmente entender isso um pouco melhor. Podemos olhar no próprio diretório work e você pode ver todos os diferentes diretórios work lá, mas isso é um pouco confuso porque estará misturado com diferentes execuções do Nextflow.

Ou podemos dizer ao Nextflow para não usar os códigos de escape ANSI.

Então se eu executar o comando novamente, mas desta vez eu digo _"-ansi-log false"_ para desligá-lo, eu também poderia usar as variáveis de ambiente _$NO_COLOR_ ou _"$NXF_ANSI_LOG=false"_. Então ele usa o tipo de estilo mais antigo de log do Nextflow sem nenhum desses códigos de escape. Apenas imprime diretamente no terminal sem atualizações inteligentes.

E agora podemos ver todos esses três processos que executaram. E cada um deles seu próprio hash de tarefa. E se formos nesses diretórios work, veremos as três diferentes saudações que especificamos.

Então isso faz um pouco mais de sentido agora. Esperançosamente você entende que o Nextflow estava fazendo isso, estava apenas sendo um pouco inteligente com o que mostrava no terminal com aqueles diretórios work.

No entanto, isso corrigiu um problema com os diretórios work, mas não corrigiu um problema com o arquivo de saída. Ainda temos apenas um arquivo de saída que diz _"Hello"_.

## 2.2. Garantir que os nomes dos arquivos de saída sejam únicos

Agora para entender isso, precisamos voltar ao nosso script de fluxo de trabalho. Estamos gerando nosso canal aqui, estamos passando-o para nosso processo, e se olharmos o processo, estamos escrevendo a saudação em um arquivo chamado _"output.txt"_ e passando aquele arquivo de saída de volta para o bloco de saída aqui embaixo, publicando-o.

No entanto, cada três vezes que este processo executa essas três tarefas diferentes. Todas elas geram um arquivo chamado _"output.txt"_, todos esses arquivos de saída são publicados no diretório results, e todos eles sobrescrevem uns aos outros. Então qualquer arquivo de resultado que você obtenha lá é apenas o último que foi gerado, mas destruiu todos os outros. Isso não é realmente o que queremos.

## 2.2.1. Construir um nome de arquivo de saída dinâmico

Existem diferentes maneiras de lidar com isso, mas a mais simples por enquanto é apenas criar nomes de arquivo únicos diferentes. Então cada vez que a tarefa executa com uma saudação diferente, ela gerará um arquivo de saída diferente, que não entrará mais em conflito quando publicado. E então teremos três arquivos de saída únicos.

Fazemos isso exatamente da mesma maneira. Podemos usar esta variável em qualquer lugar dentro do bloco script e podemos usá-la múltiplas vezes.

Então posso colar aqui, _"$\{greeting\}\_output.txt"_, e então também preciso colar aqui em cima porque não estamos mais criando um arquivo chamado _output.txt_. Então se eu não atualizar isso, o Nextflow vai travar com um erro dizendo que esperava um arquivo, que nunca foi gerado.

Então preciso fazer o mesmo lá e preciso usar aspas duplas, não aspas simples, para que esta variável seja entendida.

Ok, vamos tentar e ver se funcionou. Vamos executar o fluxo de trabalho novamente. Esperançosamente mostrará as três tarefas diferentes dentro dos três diretórios work diferentes. E com certeza, você pode ver na pasta results aqui em cima à esquerda. Agora temos três arquivos diferentes com três nomes de arquivo diferentes e cada um com o conteúdo diferente que esperamos. Então os arquivos não estão mais destruindo uns aos outros, e tudo está lá como esperamos.

Esta é uma configuração um pouco trivial pela qual passamos aqui, mas ressalta alguns dos conceitos-chave que você precisa entender sobre como a publicação de arquivos funciona, e algumas das coisas em que você pode cair como armadilhas. Então esperançosamente você pode evitar isso em seus próprios fluxos de trabalho.

Vale notar também que o que fizemos aqui é um pouco impraticável em situações da vida real. Pegamos alguns dados de entrada e estamos usando esses dados, mas também estamos nomeando o arquivo com base nesses dados, o que você geralmente não pode fazer.

Então em pipelines Nextflow mais maduros e reais, você frequentemente passará um objeto meta com todos os metadados associados a uma determinada amostra. Você pode então criar nomes de arquivo dinâmicos baseados nisso, o que é muito mais prático.

Se você está interessado em como fazer isso com as melhores práticas, há uma missão secundária em _training.nextflow.io_, que é toda sobre especificamente metadados e mapas meta, então você pode se aprofundar lá para mais detalhes.

## 3. Fornecer múltiplas entradas via um array

Ok. A seguir vamos explorar um pouco sobre como os canais são estruturados e como eles diferem de outros tipos de estruturas de dados na linguagem de codificação. E vou pensar um pouco sobre como eu poderia potencialmente usar um array, que pode ser um conceito familiar se você veio de outras linguagens.

Posso usar um array em um canal? Vamos tentar. Vou criar um array, e copiei isso da documentação, _"greetings_array"_ e _"Hello", "Bonjour"_ e _"Holà"_. E então vou colocar isso aqui em vez das minhas strings codificadas. Então vou dizer "channel.of" _"greetings_array"_, passando este array para um canal. Vamos tentar.

Abrir o terminal, e executar o pipeline.

Ok. Você pode ver que a instrução view aqui imprimiu nosso array como esperado, mas então todo esse texto vermelho, ou não será vermelho se você ainda tiver _"-ansi-log"_ desligado, mas todo esse texto vermelho está nos dizendo que algo deu errado.

Não temos mais um marcador verde aqui. Temos uma cruz vermelha, e se eu apenas tornar isso um pouco mais largo para que seja mais fácil de ler, o Nextflow está nos dizendo o que deu errado.

Então vamos quebrar isso seção por seção. Diz que o erro foi causado por, e então a razão do erro, que são arquivos de saída ausentes. Então basicamente aquele bloco de saída disse que este arquivo deveria ser criado e não foi. Em seguida diz que este é o comando que foi executado. Então isso é basicamente o conteúdo daquele arquivo _.command.sh_. É assim que parecia depois que todas aquelas variáveis foram colocadas.

E você pode ver aqui nosso comando echo na verdade só foi executado uma vez e usou o array inteiro, mas em uma representação de string, o que não é realmente o que queríamos.

E então o comando saiu assim, e aquele era o diretório work onde podemos ir e ver os arquivos para entender um pouco mais.

Ok. Então o que aconteceu foi. O Nextflow apenas passou este array inteiro como um único elemento de canal para o processo, o que significou que o processo só executou uma vez. Teve uma tarefa e não usou os dados em uma estrutura que esperávamos.

## 3.2. Usar um operador para transformar o conteúdo do canal

Então precisamos fazer algo com este canal primeiro, antes que possa ser usado. E isso está preparando o cenário para usar operadores, que são funções especiais que podemos usar em canais para manipular o conteúdo do canal.

Neste caso, vamos usar algo chamado _flatten_. Que passamos no final do canal aqui. Então criamos o canal e então executamos _flatten_. E novamente, se passarmos o mouse sobre ele, ele nos mostra a documentação para este comando imediatamente no VS Code, o que é muito útil. Você também pode encontrar toda essa documentação no site do Nextflow, a documentação.

Eu poderia apenas executar este código agora e ver se funciona, mas também é uma boa oportunidade para introduzir como fazer código dinâmico dentro de operadores e dentro do código Nextflow, que são chamados de closures.

Então vou adicionar de volta um comando view aqui antes de executarmos _flatten_. E aqui este tem essas chaves onduladas, que é o closure dinâmico. E há apenas algum código arbitrário dentro aqui que será executado, dentro do contexto de um operador view.

Aqui, isso está dizendo pegue a saudação, que é a entrada do operador view, e isso está aqui. Eu poderia chamar isso do que eu quisesse, poderia chamar isso de _"foo"_ e só preciso me referir a ele como _"foo"_ depois. E então digo com isso, retorne isso.

E então defina retornando uma string que diz antes do flatten para uma variável. muito simples.

Agora vou adicionar outro exatamente igual, mas vou dizer depois de _flatten_.

Então o que isso faz, porque isso executa em sequência, você vai ver como o canal se parece antes de executarmos _flatten_, e então novamente depois de executarmos _flatten_.

E então este canal greeting ainda é criado, então ainda vai ser passado para o processo. E esperançosamente agora o fluxo de trabalho vai executar. Vamos tentar.

Ótimo. Então primeiro de tudo é que o pipeline não travou desta vez. Tivemos três processos que executaram corretamente e temos um pequeno marcador de verificação. E então podemos ver que nossas instruções view funcionaram.

Temos antes de _flatten_, que é aquele array que vimos antes da falha, e então temos três vezes o depois de _flatten_ foi chamado onde temos _"Hello", "Bonjour"_, e todos aqueles outros três elementos separados no array, que agora são como esperávamos, três elementos separados no canal.

E você pode ver que o operador _view_ foi executado três vezes. E isso é porque este canal depois de _flatten_ agora tem três elementos. E então o operador é chamado três vezes.

Muito rapidamente, eu apenas mencionaria que quando estava criando fábricas de canais antes, fiz _"."_, e então vimos que havia muitas maneiras diferentes de criar canais, e uma delas é chamada "_fromList"_. E isso é na verdade especificamente projetado para fazer esta mesma operação. Então poderíamos ter apenas feito from list greetings away, e isso funcionará. É uma sintaxe um pouco mais limpa e agradável. Mas para os propósitos desta demonstração, queríamos torná-la um pouco mais passo a passo para que você pudesse ver como o canal está sendo manipulado e como diferentes operadores podem mudar o que está no conteúdo de um canal.

## 4. Ler valores de entrada de um arquivo CSV

Ok, como podemos tornar isso um pouco mais realista? Você provavelmente não vai querer estar criando muito código no seu pipeline Nextflow com arrays codificados. Você provavelmente vai querer pegar os dados de fora quando você lançar, e esses dados quase certamente estarão em arquivos.

Então a próxima coisa que vamos fazer é replicar isso, mas em vez de pegar os dados de um único parâmetro CLI ou de uma string ou array codificado, vamos pegá-los de um arquivo.

Então vamos nos livrar do nosso greetings away. E agora vamos mudar esta fábrica de canais novamente. Eu acabei de dizer que havia um monte para escolher e há uma chamada _".fromPath"_. E vou dizer a ela para, neste caso, pegar _params.input_, que está voltando para nossa entrada que estávamos usando anteriormente.

Agora esse parâmetro não está realmente pronto para ser usado ainda. Ainda estamos dizendo que é uma string e está codificado aqui com um padrão, mas poderíamos sobrescrever essa string. Agora queremos que isso seja um arquivo em vez disso. Então o tipo é diferente. Não é mais uma _String_. É um _Path_.

E então podemos definir o padrão se quisermos, novamente, para um Path. E se eu olhar no explorador à esquerda, você pode ver neste repositório, neste diretório de trabalho, tenho um diretório chamado data. Tenho um arquivo lá chamado _"greetings.csv"._

Então posso apenas definir o padrão aqui para _"data/greetings.csv"_. Agora, quando eu executar este pipeline novamente sem nenhuma opção de linha de comando, ele usará este valor padrão. Ele sabe que é um caminho, então sabe que deve lidar com isso como um caminho e não uma string.

E então vai passar isso para uma fábrica de canais deste _params.input_ e criar nosso canal, que então vai ser usado neste processo chamado _sayHello_. Vamos tentar.

Ok. Falhou. Não se preocupe. Isso era esperado. E se você está seguindo o material de treinamento, verá que era esperado lá também. Vamos ver o que está acontecendo aqui.

Tentou executar o pipeline. Tentou executar o processo, e obteve um erro bem similar ao que vimos antes.

Aqui diz: tentamos executar _echo_, mas em vez de ecoar o conteúdo deste arquivo CSV, apenas ecoou o caminho. E você pode ver que é o caminho absoluto completo aqui para este arquivo CSV.

E então com certeza, porque tentou escrever isso para este caminho realmente complicado, não sabia realmente o que fazer. E estava fora do escopo do diretório work do processo.

Mencionei no início que o Nextflow encapsula cada tarefa executada dentro de um diretório work especial. E se você tentar escrever para dados, que estão fora daquele diretório work, o Nextflow vai impedi-lo como uma precaução de segurança. E é isso que aconteceu aqui. Tentamos escrever para um caminho absoluto e o Nextflow falhou e nos impediu.

## 4.2. Usar o operador splitCsv() para analisar o arquivo

Ok, vamos dar uma olhada neste canal e ver como ele se parece. Podemos fazer _".view"_, e copiei isso do site. Então _.view_, e temos um closure dinâmico aqui e dizemos um nome de variável "_csv"_ como entrada. Então esse é o conteúdo do canal, e dizemos antes de splitCsv, e é assim que se parece.

Se eu executar novamente, ainda vai falhar, mas vai nos mostrar o que está dentro deste canal. Não é particularmente emocionante. É aquela variável _path_. Então você pode ver que é apenas uma string aqui porque está sendo impressa no terminal, mas é um objeto _path_, que contém a informação e metadados sobre este arquivo.

Não queremos passar os metadados do arquivo para a entrada. Queremos passar o conteúdo daquele arquivo. Se olharmos o arquivo _greetings.csv_, você pode ver aqui que tem essas diferentes variáveis aqui. _Hello, Bonjour, Holà_ novamente. E essas são as coisas que realmente queremos estar passando para nosso processo, não apenas o arquivo em si como um único objeto.

Então precisamos analisar este arquivo CSV. Precisamos desempacotá-lo, chegar ao conteúdo do arquivo CSV, e então passar o conteúdo dentro do canal para o processo.

Como você provavelmente pode dizer pela mensagem de log, queremos usar o _splitCsv_, que é outro operador, outro operador de canal. Então se eu fizer "_dot" "s"_, e então você pode ver que foi auto sugerido. Ops, _splitCsv_ e alguns parênteses.

E então depois de _splitCsv_, vou colocar outra instrução _view_ apenas para que possamos ver como fica depois. Vamos executar o pipeline e ver o que temos.

Ok. Ainda falhou, mas de uma maneira nova e emocionante, o que é progresso.

Desta vez novamente, temos algum problema com nosso script, que foi renderizado. Agora. Não temos mais o caminho final, mas temos um array de variáveis, que se parece muito com o erro que tivemos anteriormente quando estávamos passando um array como uma entrada fixa.

Com nosso log do operador view, podemos ver antes de _splitCsv_ era o caminho. E com certeza, depois de _splitCsv_, temos três saídas diferentes e cada uma dessas saídas se parece muito com cada uma das linhas do arquivo _greetings.csv_, o que faz sentido.

Então o que aconteceu aqui é que o Nextflow analisou este arquivo CSV nos deu três objetos, um array para cada linha do arquivo CSV. Então três vezes passamos um array de variáveis para o canal em vez de um único valor de string.

Ok, então da última vez que tivemos este problema, usamos _flatten_. Vamos apenas muito rapidamente. Tentar flatten e ver o que acontece.

Posso chamar essas variáveis do que quiser. Então vou chamá-la de _myarray_ porque não é mais realmente um CSV. Vamos tentar executar novamente e ver o que acontece com _flatten_.

Então desta vez vamos executar, analisamos o CSV em três objetos de array, e então achatamos. E desta vez passou. E o pipeline Nextflow executou. No entanto você pode ver que _flatten_ realmente vai com tudo e achata tudo. E então obtemos três entradas de array independentes para cada linha. E então executou o processo três vezes cada linha de um CSV. E agora temos um monte de arquivos de resultados, e 123, 456, e todos os tipos de coisas, não apenas aquela primeira coluna do CSV, que é o que realmente queríamos.

## 4.3. Usar o operador map() para extrair as saudações

Então como chegamos apenas à primeira coluna? Se flatten é muito simplista aqui, precisamos de um operador mais complexo onde possamos realmente personalizar e dizer o que queremos do CSV.

Para fazer isso, vamos usar _map_. Basicamente _map_ apenas diz, execute algum código, alguma função sobre cada elemento que eu receber e faça algum tipo de transformação nele. E porque é tão flexível, você verá aparecer no código Nextflow o tempo todo.

Por si só, não faz nada. Então não queremos parênteses regulares, queremos um closure aqui e precisamos dizer o que fazer. Então vou dizer _"row"_, porque está sendo dado linhas do CSV, então é um nome de variável lógico. É a entrada. E quero retornar apenas o primeiro elemento daquele array.

Arrays no Nextflow são baseados em zero, então vamos dizer apenas o primeiro elemento, que é linha zero. Se quiséssemos a segunda coluna, eu poderia ser um ou a terceira coluna ser dois, e assim por diante. Podemos retornar o que quisermos aqui, mas vou retornar apenas o primeiro valor.

E agora, podemos executar o pipeline novamente e ver se faz o que esperamos.

Com certeza, depois de _splitCsv_ temos nossos arrays, e então depois do _map_, temos nossas strings limpas e agradáveis, apenas _"Hello", "Bonjour"_ e _"Holà"_. E o pipeline agora está fazendo o que queremos. Fantástico.

Então podemos nos livrar de todos esses comandos view agora. Não precisamos mais deles.

## Recapitulação

Terminamos nossa depuração e este é o código com que terminamos. Pegando nosso parâmetro CLI chamado _input_, que é classificado como um _Path_. O Nextflow encontra o caminho, carrega-o e entende o arquivo CSV. Retorna todas as diferentes linhas. E então mapeamos apenas o primeiro elemento daquela linha no canal que meio que nos dá o conteúdo do canal, que é passado para o processo.

E o processo executa sobre cada elemento no canal, que são três. E executa o processo três vezes, dando-lhe três tarefas. E esses resultados são então publicados do fluxo de trabalho, capturados pela saída do processo. Publicados de um fluxo de trabalho e salvos no bloco de saída para um subdiretório chamado _"hello_channels"_.

Muito legal. Estamos chegando agora a algo que se assemelha mais a um pipeline Nextflow da vida real que você pode executar para alguma análise real.

## Conclusão

Ok. Esperançosamente você está agora tendo uma noção do que são canais e operadores do Nextflow e como os operadores trabalham em canais e como você pode criá-los.

Canais, como eu disse no início deste vídeo, são a cola do Nextflow. E você pode ver aqui que podemos pegar diferentes entradas e manipulá-las e pegar esses dados e então passá-los para a lógica de fluxo de trabalho downstream.

E este bloco de fluxo de trabalho aqui é realmente onde você constrói toda aquela paralelização e toda a lógica inteligente, e explica ao Nextflow como construir seu DAG de fluxo de trabalho, e como orquestrar seu pipeline.

Canais não são o conceito mais fácil de entender. Então faça uma pausa, pense um pouco sobre isso, talvez leia o material novamente, e realmente certifique-se de que você entendeu esses conceitos porque isso é fundamental para seu entendimento do Nextflow e quanto melhor você entender canais e os diferentes operadores de canal e as diferentes fábricas de canais. Mais divertido você terá escrevendo Nextflow e mais poderosos seus pipelines serão.

Isso não é o mesmo que programação regular em Python ou outras linguagens. Não estamos usando instruções _if_ aqui, isso é programação de fluxo funcional usando canais e operadores. Então é um pouco diferente, mas também é super poderoso.

Esse é o fim deste capítulo. Vá e faça uma pausa rápida e te vejo no próximo vídeo para a parte três onde vamos passar por Hello Workflow, e falar um pouco mais sobre os fluxos de trabalho.

Assim como o capítulo anterior, há algumas questões de quiz na parte inferior da página web aqui, então você pode dar uma passada rápida por elas e certificar-se de que entende todas as diferentes partes do material que acabamos de fazer. E além disso, te vejo no próximo vídeo. Muito obrigado.

Ok.
