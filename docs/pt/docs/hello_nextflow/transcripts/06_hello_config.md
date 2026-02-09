# Parte 6: Hello Config - Transcrição do Vídeo

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Tradução assistida por IA - [saiba mais e sugira melhorias](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

<div class="video-wrapper">
  <iframe width="560" height="315" src="https://www.youtube.com/embed/FcZTiE25TeA?si=tnXTi6mRkITY0zW_&amp;list=PLPZ8WHdZGxmWKozQuzr27jyMGqp9kElVK&amp;cc_load_policy=1" title="YouTube video player" frameborder="0" allow="accelerometer; autoplay; clipboard-write; encrypted-media; gyroscope; picture-in-picture; web-share" referrerpolicy="strict-origin-when-cross-origin" allowfullscreen></iframe>
</div>

!!!note "Notas importantes"

    Esta página mostra apenas a transcrição. Para instruções passo a passo completas, retorne ao [material do curso](../06_hello_config.md).

    Os números de seção mostrados na transcrição são fornecidos apenas para fins indicativos e podem não incluir todos os números de seção dos materiais.

## Boas-vindas

Olá, e bem-vindo de volta à Parte Seis de Hello Nextflow. Esta seção é toda sobre configurações, e é a última parte deste curso.

Nextflow é particularmente bom em duas coisas: reprodutibilidade e portabilidade. Configurações é onde vemos a segunda dessas características realmente brilhar. A capacidade de configurar um pipeline Nextflow para executar de diferentes maneiras e funcionar em diferentes sistemas, sem ter que editar o código subjacente do pipeline.

Este superpoder permite que pipelines Nextflow sejam reutilizados por outras pessoas em lugares diferentes, ou através de diferentes infraestruturas às quais você mesmo pode ter acesso.

Isso significa que você pode desenvolver código de pipeline no seu laptop, enviá-lo para a nuvem, executá-lo no seu HPC, e é o mesmo código de pipeline e ele roda em todos os lugares.

Nesta seção, vamos passar por alguns tópicos. Começaremos com como o Nextflow lida com arquivos de configuração, de onde ele os carrega, e como você os escreve e como os estrutura, e essa separação entre o pipeline em si e o que deve ir em um arquivo de configuração.

Então passaremos para alguns casos de uso comuns, como mudar onde os arquivos de saída são armazenados, e também como fazer o pipeline funcionar em diferentes infraestruturas, tanto usando diferentes tipos de empacotamento de software quanto submetendo jobs para diferentes infraestruturas.

## Hierarquias de arquivos de configuração

Ok, vamos começar. Quando se trata de carregar arquivos de configuração, o Nextflow pode buscar de muitos lugares diferentes, o que é uma coisa boa e também pode ser uma coisa um pouco arriscada porque às vezes pode ser um pouco difícil saber de onde ele está obtendo um arquivo de configuração e em que ordem ele carrega as coisas.

Então eu realmente recomendo que você clique neste link aqui, que nos leva à documentação do Nextflow. E nesta página de configuração, ela lista os principais lugares de onde a configuração é carregada, e importante, a ordem em que essas coisas são carregadas.

Então você pode ver, você pode colocar um arquivo de configuração no seu diretório home do Nextflow, que é tipicamente ".nextflow" no seu diretório home. E esse arquivo sempre será carregado por toda execução do Nextflow no seu sistema.

O próximo lugar a procurar é um arquivo na raiz do seu repositório ou diretório de pipeline chamado "nextflow.config".

Depois disso, outro arquivo chamado "nextflow.config", mas desta vez no diretório de onde você está lançando o Nextflow: o diretório de lançamento.

Finalmente, você pode fornecer caminhos de arquivos de configuração na linha de comando com um argumento "-c", e você pode fazer isso várias vezes. E eles são aplicados na ordem em que você os especifica.

Você pode fornecer arquivos de configuração em todos esses locais se quiser, e eles serão carregados iterativamente, cada um sobrescrevendo o anterior apenas nos escopos de configuração onde eles conflitam.

Este é um sistema realmente poderoso porque significa que você pode definir padrões sensatos e então ficar gradualmente mais e mais específico à medida que você se concentra naquela configuração.

## 0. Aquecimento: Execute hello-config.nf

Ok, vamos fechar isso e pular para nossos Codespaces e começar. Como antes, eu limpei aqui, removi meus diretórios de resultados anteriores, meu Nextflow e meus diretórios de trabalho e assim por diante. Não se preocupe se você ainda tem esses arquivos por aí. É só porque estou muito ampliado e então as coisas ficam bagunçadas muito rapidamente caso contrário.

Vamos trabalhar com hello-config.nf, o último arquivo em nosso diretório, e isso deve continuar de onde paramos na seção anterior.

Então temos nossos quatro processos diferentes, que estão sendo incluídos de arquivos de módulo. Temos nossos parâmetros de pipeline, nosso bloco de workflow onde estamos chamando os diferentes processos e conectando os canais, publicando os canais de saída, e então o bloco de output na parte inferior onde definimos onde esses arquivos devem ser armazenados e como eles devem ser copiados.

Também já temos um arquivo "nextflow.config" do último capítulo, onde habilitamos o Docker, e vamos construir sobre este arquivo hoje.

Como antes, mudamos o caminho de saída neste script principal para hello config, apenas para que não conflite com resultados anteriores que você gerou.

Ok, vamos apenas verificar rapidamente se tudo ainda está funcionando como esperamos. Abrir um terminal e fazemos nextflow run hello-config.nf. O Nextflow carrega. Deve executar nossos quatro processos diferentes. Gerar alguma arte ascii legal usando cowpy e então salvar nossos resultados em nossos arquivos de resultados naquele diretório.

Posso dar uma olhada rápida aqui apenas para ter certeza de que esses arquivos parecem como esperamos, e com certeza, lá está nosso Peru gigante. Ótimo.

## 1.1. Mova valores padrão para nextflow.config

Agora a primeira coisa que vamos fazer é começar a mover algumas coisas do nosso script para nosso arquivo de configuração.

E o que nos importa são principalmente os parâmetros neste estágio. Queremos levar os valores padrão para o arquivo de configuração, para que fique mais claro quais são os padrões e seja mais fácil para as pessoas sobrescrevê-los.

Vou apenas pegar este bloco params aqui do script e colocá-lo no arquivo de configuração. E precisamos ter um pouco de cuidado aqui, porque agora a sintaxe é ligeiramente diferente entre configuração e scripts. O arquivo de configuração não pode aceitar declarações de tipo porque não estamos realmente definindo esses params, estamos apenas referenciando-os. Então vou me livrar deles.

Mas fora isso é muito parecido. Temos um bloco params e então temos nossos diferentes parâmetros de entrada, parâmetro batch, parâmetro character.

Agora posso voltar ao meu script e não preciso mais definir esses padrões porque esses valores agora estão no meu arquivo Nextflow config.

No entanto, deixo os nomes dos parâmetros e seus tipos, para que o Nextflow saiba essa informação e ainda possa fazer toda a segurança de tipo e tudo mais.

Ok. Salvamos esses arquivos e verificamos rapidamente se tudo ainda funciona da mesma forma que antes. Não deve haver mudanças aqui. Mantivemos os valores iguais. Apenas movemos onde eles foram definidos.

Ótimo.

## 1.2. Use um arquivo de configuração específico da execução

Agora, até agora estivemos lançando o Nextflow do mesmo diretório onde temos nosso script de pipeline. Então nosso diretório de lançamento e nosso diretório de pipeline são meio que a mesma coisa.

Para mostrar como podemos ter diferentes arquivos de configuração com diferentes diretórios de lançamento, vamos criar um novo subdiretório agora.

Então vou dizer mkdir, e vamos chamá-lo de tux-run.

E então vou fazer cd, mudar de diretório para tux-run. E note que agora estamos em nosso diretório de trabalho não está mais no mesmo diretório que os scripts de pipeline.

Ok, vamos criar um novo arquivo "nextflow.config". Então touch Nextflow config, e vamos apenas abri-lo no VS Code. Você pode ver também na barra lateral aqui que agora estamos neste subdiretório.

Agora podemos pegar o mesmo bloco params que tínhamos no nextflow.config de nível superior, copiar isso e agora podemos mudar esses valores.

Primeiro, os dados agora são um caminho relativo diferente porque estamos em um subdiretório, então precisamos atualizar isso. E então vamos mudar batch para experiment, e vamos mudar o character de Turkey para tux.

Agora clique em salvar aí, e vamos testar. Assim como com data, preciso agora dizer ../ para chegar ao script. Então é Hello config. E pressiono enter.

O código do pipeline não mudou nada, mas agora vamos ter dois conjuntos de configuração carregando, e o arquivo de configuração do diretório de lançamento deve sobrescrever os padrões, que foram definidos no nextflow.config do pipeline, e devemos obter diferentes conjuntos de resultados.

Com certeza, dentro do nosso diretório aqui, dentro de tux-run, você pode ver que temos um diretório dot Nextflow e um diretório work e isso é porque estes são criados sempre no seu diretório de lançamento. Então estes são diferentes dos de work e results que tínhamos de execuções anteriores.

Agora, se eu olhar em results, podemos ver nosso collected e lá está nosso pequeno personagem tux. Então você pode ver que esses parâmetros foram interpretados adequadamente.

## 1.3. Use um arquivo de parâmetros

Ok. Antes quando eu estava falando sobre os diferentes arquivos de configuração que poderiam ser carregados, deixei de fora um outro lugar de onde podemos obter configuração.

Você pode obtê-la da linha de comando como vimos com dash dash nomes de parâmetros, mas também podemos fornecer um arquivo YAML ou JSON, apenas de params.

O arquivo de configuração pode ter todos os diferentes tipos de escopos, mas esses arquivos são apenas parâmetros, e é uma maneira amigável ao usuário de fornecer muitos parâmetros de uma vez, e talvez uma maneira um pouco mais reproduzível porque você os escreve em arquivo, então é fácil obtê-los em um estágio posterior.

Então vamos voltar ao nosso terminal e antes que esqueçamos, certifique-se de que voltamos um diretório, então não estou mais no subdiretório, e vou olhar o arquivo YAML que temos aqui chamado test-params.yaml.

Então se eu apenas fizer code test-params.yaml, você pode ver que este é apenas um arquivo YAML regular. Nada de especial sobre ele. Com as chaves sendo nossos nomes de parâmetros, com a formatação YAML então dois pontos aqui, e então um valor.

Note que este não é código Nextflow, então não podemos colocar coisas como variáveis aqui. Estes são apenas valores estáticos.

Também porque JSON na verdade analisa como YAML, também podemos ter um arquivo test-params.json, que parece muito similar. É apenas um formato de dados diferente.

Então temos dois arquivos de teste diferentes aqui e temos variáveis ligeiramente diferentes.

Ok, então como damos estes ao Nextflow? É muito simples. Fazemos Nextflow run hello config, como antes. E em vez de "-c" para arquivo de configuração, ou carregar aqueles nomes de arquivo padrão, fazemos -params-file. Hífen único porque é uma opção central do Nextflow.

E então passamos o caminho para aquele arquivo. Então vou fazer "-params-file test-params.yaml", e veremos se esses são carregados adequadamente.

Ok. Executou. Vamos apenas nos lembrar do que estava neste arquivo YAML. Então o batch foi definido como YAML, então é assim que deve ser chamado, e deve ter um stegosaurus. Então vamos subir e olhar em results. E temos COLLECTED-yaml. Então vamos ver se temos um Stegosaurus. Fantástico, um Stegosaurus usando um chapéu. É isso que gostamos.

Então isso funcionou muito bem, e é exatamente o mesmo com o arquivo JSON. Apenas trocamos a extensão do arquivo aqui e o Nextflow sabe como ler isso.

E neste caso, devemos ter um batch chamado JSON e devemos ter uma tartaruga. Vamos dar uma olhada. Maravilhoso. Uma das minhas ferramentas CLI favoritas.

## 2.1. Personalize o diretório de saída com -output-dir

Ok, então isso tem sido principalmente pensar sobre entradas para o pipeline e mudar parâmetros. E quanto às saídas?

Agora, embora tenhamos estado mudando os subdiretórios usando params, você pode ter notado que todos os nossos arquivos ainda estão indo para results.

Podemos mudar aquele diretório base para o qual todos os arquivos são publicados com uma flag de linha de comando chamada -output-dir. Então se eu fizer Nextflow run hello config, e então faço -output-dir, e vamos chamá-lo de "custom-outdir-cli". Não consigo digitar. Apenas para lembrarmos de onde esses arquivos vieram.

Esta é uma opção central do Nextflow e é muito nova. Isso foi adicionado apenas recentemente, e esta é uma das coisas que podemos fazer com o novo analisador de linguagem e tudo mais.

É um bocado para digitar. Você também pode apenas chamá-lo de "-o" se quiser. Então se eu apenas voltar. Posso apenas encurtar isso para "-o", que é um pouco mais simples.

Ok. Executamos isso. Não mudamos nada em nosso pipeline ou mesmo em nossa configuração neste ponto, e deve esperançosamente salvar todos os nossos resultados em um diretório de nível superior diferente. E você pode imaginar que pode definir isso para basicamente qualquer caminho que quiser.

Acabou de chegar no topo. Temos um custom-outdir-cli, e todos os arquivos estão organizados lá exatamente da mesma maneira, com seus mesmos subdiretórios e nomes de arquivo. Então esta é uma maneira realmente fácil de apenas mudar onde o pipeline publica seus resultados, sem pensar muito sobre como esses resultados são organizados.

## 2.1.2. Remova caminhos codificados do bloco output

Se eu olhar dentro deste diretório, podemos ver que ainda temos um subdiretório chamado Hello Config, que parece um pouco redundante agora.

Então vamos apenas carregar nosso script novamente e agora podemos remover aquele subdiretório do bloco output na parte inferior. Porque não precisamos mais dele realmente. Então podemos apenas fazer isso agora, deletar isso daqui. E então se é apenas isso, você pode ou deletar isso completamente ou deixar como uma string vazia. Vou deixar como uma string vazia por enquanto, porque vamos voltar e colocar algumas coisas diferentes em seu lugar no futuro. Mas se você não se importa com subdiretórios, é mais limpo apenas remover completamente a declaração path lá.

Ok, vamos salvar. Apenas tentar rapidamente de novo. Na verdade vou remover meu diretório "custom-outdir-cli" para não ficarmos confusos com quaisquer arquivos existentes lá. Porque lembre-se, quando você publica coisas, não remove os arquivos que estavam lá antes. Apenas adiciona novos. Vamos executar aquele comando novamente, custom-outdir-cli.

E agora se você fizer "ls custom-outdir-cli", não há mais diretório lá chamado Hello Config.

## 2.2.1. Defina outputDir no arquivo de configuração

Ok, a flag de linha de comando aqui, "-o" ou "-output-dir" é boa. Mas e quanto a definir padrões para isso em config? Como fazemos isso?

Abro o arquivo "nextflow.config", fecho todo o resto e me livro disso. Podemos adicionar uma nova opção de configuração aqui, que eu apenas copiei do site de material de treinamento, e é chamada outputDir.

Não está sob nenhum escopo. Não está sob params ou nada. É de nível superior, e podemos definir isso como uma string. Agora uma coisa simples a fazer é apenas mudá-la para qualquer coisa diferente de results como uma string codificada. Mas porque isso está em um arquivo de configuração Nextflow, podemos ser um pouco espertos aqui e também incluir variáveis.

E você pode ver aqui que incluímos uma variável params, params.batch, que é parte desta string. Isso significa que podemos reutilizar variáveis que estão vindo de outros lugares. E neste caso, se fizermos --batch, quando executarmos o Pipeline Nextflow, vamos obter um subdiretório em nosso caminho personalizado baseado em qual era o nome do batch.

Ok, então vamos testar isso e apenas dar uma olhada rápida para ver como, como os resultados parecem. Então se eu fizer Nextflow run hello config e --batch my_run. Vamos nos lembrar de como a configuração parecia. Então é custom-outdir-config.

Tree custom-outdir-config. E você pode ver que o batch foi chamado my_run. E então temos aquele subdiretório chamado my_run. Então aquele caminho de arquivo dinâmico funcionou.

E não apenas isso, não foi mais para um diretório results padrão, e eu não tive que especificar nada na linha de comando para mudar o diretório base. Então redefinimos com sucesso o valor padrão para o outputDir padrão.

## 2.2.2. Subdiretórios com nomes de batch e processo

Ok, vamos levar isso um pouco mais longe. Essa é uma variável dinâmica dentro do arquivo de configuração. E quanto ao script em si? Agora, até agora tivemos esses caminhos aqui e estes também podem ser dinâmicos. Então em vez de apenas codificar algo, podemos colocar alguns colchetes ondulados e colocar algo dinâmico.

Então por exemplo, temos nossos processos chamados sayHello. Poderíamos fazer sayHello.name, que é um atributo do processo, que é meio chato porque é apenas "sayHello" neste caso. Mas é variável.

Então isso te dá uma ideia. Então podemos colocar isso aqui e dizer convertToUpper.name, collectGreetings.name, collectGreetings.name novamente, e cowpy.

Agora quando executarmos, o diretório base ainda vai ser custom-outdir-config. E vai estar em um subdiretório chamado params.batch, mas os subdiretórios sob isso devem ser organizados por nome de processo.

Vamos apenas testar isso e ver se funciona. Então vou remover o diretório anterior para não ficarmos confusos, e apenas usar exatamente o mesmo comando Nextflow Run.

Deve executar da mesma maneira. Eu poderia estar usando dash resume em todos estes para torná-lo um pouco mais rápido e usar os resultados previamente calculados. Agora, se eu fizer tree custom-outdir-config, você pode ver que não está em results, está em nosso diretório base com o nome do batch. E você pode ver que todos os resultados agora estão organizados dentro de subdiretórios nomeados após o processo. Então temos dois lugares diferentes onde estamos definindo caminhos de saída dinâmicos aqui.

Ok. Última coisa, vamos adicionar de volta aquelas pastas intermediárias, que tínhamos antes porque eram meio legais. Intermediates.

E também podemos pensar um pouco sobre este params.batch, talvez como desenvolvedor de pipeline eu realmente gostei de ter isso no subdiretório, mas se usuários finais do pipeline estão definindo "-o" ou -output-dir no CLI, está completamente sobrescrevendo toda esta declaração, e perdemos aquele subdiretório.

Então o que podemos fazer é podemos tirar aquele caminho dinâmico do outputDir config, que seria esmagado, e colocá-lo no caminho de output, que não é esmagado.

Então podemos fazer params.batch barra intermediates barra sayHello.name, e fazer tudo isso em uma string entre aspas duplas, então é interpolado pelo Nextflow.

Posso agora copiar, ops. Copiar estes para os outros processos. Lembre-se de colocá-los todos entre aspas. E remover intermediates destas saídas particulares.

Ok? Está parecendo ligeiramente mais complexo agora, mas você pode ver que estamos realmente começando a construir uma estrutura de diretório de saída bem organizada em nosso código.

E o que é realmente legal é que essa complexidade extra no código não passa para o CLI. Então podemos executar nosso comando com -output-dir e quaisquer variáveis batch, apenas pensando sobre como executar o pipeline e não realmente pensando muito sobre o que está no código. E nossos arquivos de saída vão ser construídos muito bem de uma maneira muito bem organizada, o que é legal para pessoas usando o pipeline basicamente.

Ótimo. Enquanto escrevo isso, percebo que cometi um erro. Veja se alguém me pegou aqui. Temos collectGreetings.name, então algo deu errado. E sim, com certeza, eu acidentalmente esqueci de colocar estes em colchetes ondulados.

Então lembre-se, seja cuidadoso quando estiver escrevendo seu código e certifique-se de que você diz ao Nextflow o que é uma variável e o que é apenas uma string. Porque ele vai fazer exatamente o que você mandar fazer. E nada mais. Como todos os bons computadores. Ok, isso deve consertar.

## 2.3. Defina o modo de publicação no nível do workflow

Há um pedaço deste script, que eu ainda não amo, que é o fato de que estamos escrevendo mode copy de novo e de novo, e se há uma coisa que não gostamos, é nos repetir.

Então podemos limpar isso um pouco pegando isso e movendo para a configuração. E de fato, podemos definir isso para o pipeline inteiro de uma vez. Então não temos que dizer isso várias vezes.

Vamos para nosso arquivo de configuração e temos um novo escopo aqui chamado workflow. E podemos ou fazer colchetes ondulados ou podemos fazer notação de ponto. Não faz diferença. O site de material de treinamento usa notação de ponto. Posso dizer output e podemos misturar e combinar, então mode equals copy. Ótimo.

E agora podemos voltar aqui e deletar estes. Agora poderíamos deixá-los no lugar. A configuração está basicamente sobrescrevendo o que está escrito aqui, mas como temos isso na configuração de nível de pipeline, e esses dois arquivos são enviados juntos, não há razão para realmente fazer isso duas vezes.

Ok. Apenas uma verificação de sanidade, porque aparentemente cometemos erros. Vamos executar isso novamente e apenas verificar que estamos usando corretamente o modo copy para publicar arquivos. Então vamos executar o script novamente e desta vez colocamos os resultados em um diretório chamado config-output-mode, ver como os arquivos parecem lá.

E então se eu fizer "ls -l" para olhar batch, e podemos olhar cowpy, por exemplo. E devemos ver, sim, que este é um arquivo adequado aqui, que não é um link simbólico, então aquele atributo de configuração foi aplicado adequadamente.

## 3. Selecione uma tecnologia de empacotamento de software

Ok. Até agora estivemos focando nas entradas e saídas, os arquivos com os quais o fluxo de trabalho está executando. Mas e quanto à infraestrutura? Eu disse no início que o Nextflow permite que você execute o mesmo pipeline em diferentes configurações de computação. Então como isso parece?

Para mostrar isso, vamos mudar de usar Docker para executar cowpy, e em vez disso usaremos Conda para fazer a mesma coisa.

Posso fazer isso muito simplesmente. Se eu for para code, "nextflow.config". Se você se lembra no topo, definimos docker.enabled anteriormente, e no último capítulo para que pudéssemos usar o contêiner com cowpy dentro.

Vou dizer ao Nextflow para não usar Docker. Definir isso como false. E vou dizer Conda enabled equals true. Então dizer ao Nextflow, por favor use Conda.

Agora apenas habilitar Conda não é suficiente por si só. Assim como fizemos com Docker, temos que dizer ao Nextflow onde ele pode obter o software que precisa.

Então se pularmos para os módulos aqui. E abrir o script cowpy. Podemos ver que temos uma declaração container no topo. E o contêiner é usado pelo Docker, mas também Singularity, Apptainer, e muitas das outras ferramentas de software.

Mas não pode ser usado para Conda, então temos uma declaração separada chamada "conda", e poderíamos apenas escrever "cowpy". E isso deixará para a resolução de pacotes conda descobrir a melhor maneira de resolver isso de acordo com seu ambiente conda local.

Ou é uma boa prática fazer o que o site de material de treinamento diz para fazer, que é definir um canal conda específico com sua notação de dois pontos duplos, e definitivamente definir uma versão específica do software para que cada pessoa que execute o pipeline obtenha a mesma versão.

Note que contêineres são um pouco superiores neste aspecto, porque quando você instala algo com Conda, ainda vai descobrir todas as dependências para aquele pacote, e elas podem mudar ao longo do tempo. Chamado de desvio de dependência.

Então contêineres, no entanto, travam toda a pilha de dependências de software até o fim, então você pode estar um pouco mais confiante de que A, vai funcionar, e B, será reproduzível.

Então se você for capaz de usar Docker ou Singularity ou Apptainer, eu definitivamente recomendaria isso.

Agora o que é legal sobre isso é que o arquivo de módulo, que é escrito pelo desenvolvedor do pipeline, agora tem tanto Container quanto Conda, e então estamos dizendo à pessoa que está executando este pipeline, não nos importamos qual solução de empacotamento de software você usa. Vai funcionar tanto com Docker quanto com Conda, e é aqui que obter o software em ambos os casos.

Podemos abrir o terminal e vamos tentar isso. Então Nextflow run hello config --batch conda. E a primeira vez que isso executa com conda, vai ser um pouco lento quando chegar naquele processo particular, porque tem que executar "conda install".

E está criando um ambiente conda especial apenas para este único processo. Então não está usando meu ambiente conda global, que tenho no meu terminal. Está criando um apenas para aquele único processo. Isso é bom porque evita coisas como conflitos de dependência entre diferentes processos em seu fluxo de trabalho. Se seus processos têm ferramentas que precisam de diferentes versões de Python ou coisas assim, tudo bem porque estão usando diferentes ambientes conda.

O Nextflow armazena em cache esses ambientes conda localmente, você pode ver que ele te diz onde está aquele caminho, está no diretório work aqui. E então da próxima vez que eu executar este script com Conda, será muito mais rápido porque vai encontrar aquele ambiente conda existente e apenas reutilizá-lo. Mas a primeira vez que fazemos isso, tem que ir e buscá-lo, resolvê-lo, baixar todas as dependências, e configurar tudo.

Ok, ótimo, executou. Podemos apenas nos lembrar do que o pipeline está atualmente configurado para usar. Se olharmos no arquivo de configuração, era "custom-outdir-config" agora para mim. Veja se eu subir para aquele diretório base. E eu fiz --batch conda. Lá está nosso subdiretório conda. Então funcionou e lá está nossa saída cowpy.

Então buscou cowpy, instalou no meu sistema local usando conda, e executou o processo. E o que é ótimo é que, como aquele usuário final, eu não tive que pensar nada sobre qualquer gerenciamento de software lá. O Nextflow apenas resolveu para mim. Eu disse, preciso usar conda neste sistema. O desenvolvedor do pipeline disse quais pacotes eu precisava. E o Nextflow fez o resto. Muito poderoso.

Note que você pode na verdade usar uma mistura de diferentes tecnologias. Então posso habilitar Docker para processos específicos, e conda para outros processos, ou dizer que alguns processos devem apenas usar qualquer software local que eu tinha instalado. Isso é bem incomum, mas é possível, e em alguns casos, por exemplo, se você está usando certo software que pode ser difícil de empacotar no Docker, você tem uma saída.

## 4. Selecione uma plataforma de execução

Então isso é empacotamento de software. A outra parte da portabilidade para outros sistemas é onde os jobs realmente executam. No momento, estou executando basicamente no meu laptop ou neste Codespaces, que é um único computador. Não há nada sofisticado. O Nextflow está sendo um pouco esperto sobre paralelizar os jobs da melhor forma que pode, mas está tudo em um sistema.

Agora, se você está executando em um HPC, você provavelmente tem algum tipo de agendador de jobs como SLURM ou PBS ou algo assim, e você submeterá jobs para aquele agendador e ele distribuirá todos os jobs para diferentes nós de computação.

Outra maneira de executar é na nuvem. Então talvez você esteja usando AWS Batch, ou Azure Cloud, ou Google. E todos estes funcionam em um sistema similar onde você tem um agendador e você submete jobs e eles são submetidos para diferentes lugares para serem computados.

Agora no passado distante quando comecei a fazer bioinformática, o software de todos para executar análise era muito ligado à sua infraestrutura computacional, o que tornava quase impossível replicar.

Mas com esta separação de configuração no Nextflow, e com a capacidade do Nextflow de interagir com muitos backends de infraestrutura de computação diferentes, é muito simples pegar nosso pipeline sem modificar o código do pipeline de forma alguma e apenas trocar isso.

## 4.1. Mirando um backend diferente

Então se formos ao nosso arquivo "nextflow.config", e agora podemos colocar alguma configuração de nível de processo. Então se eu colocar no topo escopo process e posso definir o executor, e aqui está definido como local, que é o padrão.

Note que porque isso é nível de processo, podemos mirar coisas para diferentes processos. E então você pode na verdade configurar executores para serem específicos de processo e ter uma execução híbrida, onde alguns jobs podem executar localmente, onde quer que o job Nextflow esteja sendo executado. Alguns são submetidos para diferentes HPC e alguns podem ser submetidos para a nuvem. Você pode ser tão esperto quanto quiser.

Agora, é muito difícil demonstrar isso em um ambiente de treinamento como este porque não tenho um HPC para submeter. Mas posso fazer é se eu digitar slurm, podemos trapacear um pouco e você pode ter uma sensação disso.

E isso é realmente mais interessante para pessoas que estão acostumadas a executar no SLURM e sabem como os cabeçalhos SLURM parecem. Mas se eu fizer Nextflow run, hello config. Vai falhar porque vai tentar submeter jobs para um cluster que não existe. Então vamos obter algum tipo de erro sobre sbatch não estar disponível.

Sim, escrito. Essa é a ferramenta. Essa é a ferramenta CLI que você usa para submeter jobs para um cluster slurm. Mas o que podemos fazer é podemos ir e olhar em nosso diretório work aqui por comando, clicar, abrir aquele diretório e olhar o .command.run. E você pode ver no topo do arquivo .command.run, temos nossos cabeçalhos sbatch, dizendo a um cluster SLURM teórico como lidar com esta submissão de job.

E então você pode ver que o Nextflow está sendo esperto, está fazendo todas as coisas certas. É só que não tínhamos um cluster para submeter.

## 5. Controle alocações de recursos computacionais

O que mais é diferente entre diferentes infraestruturas de computação? Outra coisa é quanto de recursos disponíveis você tem, e de fato, em muitos ambientes de computação, é um requisito que você tenha que especificar quantos CPUs e quanta memória um job precisa.

Novamente, o Nextflow abstrai isso para nós, de modo que não é mais específico para um único tipo de ambiente de computação, e podemos digitar no escopo de nível de processo aqui. CPUs equals one, memory equals two gigabytes. Nosso pipeline não é muito exigente, então isso deve estar bem.

Agora, eu apenas adivinhei esses números aqui, mas como você sabe qual é uma quantidade sensata de recursos para usar? É um trabalho bastante difícil ir e cavar através de todos esses diferentes processos de um grande pipeline de muitas amostras e entender qual foi a utilização de recursos.

Então uma boa abordagem para isso é definir esses valores para números altos para começar, apenas para que seu pipeline execute sem quaisquer erros, e então pedir ao Nextflow para gerar um relatório de uso para você.

Isso é super fácil de fazer, então vou voltar a um terminal. Oh, preciso lembrar de definir isso de volta para local para que meu pipeline realmente execute. E vou dizer Nextflow run, e vou usar uma flag de linha de comando -with-report.

E posso deixar isso em branco e ele dará um nome de arquivo padrão, mas vou dar um nome de arquivo específico para que isso seja salvo em um lugar específico.

Pressiono Enter, e o pipeline executa exatamente como normal, mas quando terminar, vai gerar um relatório HTML legal para mim.

Então na barra lateral aqui, tenho este arquivo HTML. Se eu estivesse executando isso localmente, apenas abriria. Estou, porque estou em Codespaces, vou clicar com o botão direito nisso e clicar em download, o que vai baixá-lo para meu computador local. E posso apenas facilmente abri-lo no navegador web.

O Nextflow pode gerar um relatório como este para qualquer pipeline e tem algumas informações realmente legais. Então é uma boa prática sempre salvar essas coisas. Nos diz quando executamos, onde executamos, se foi bem-sucedido ou não, quais parâmetros foram usados, qual foi o comando CLI, coisas assim.

E também há esses gráficos sobre uso de recursos. Então nos diz qual porcentagem de chamadas de CPU foram usadas para cada processo como um gráfico de caixa aqui, porque há muitas tarefas para cada processo, então podemos ver a distribuição.

Você pode ver nossos processos aqui, cowpy e collectGreetings apenas tiveram uma única tarefa, então é apenas uma única linha. E temos tanto CPU quanto memória e duração do job, e foram muito rápidos.

Se você está usando Seqera Platform, a propósito, você obtém os mesmos gráficos integrados na interface da Platform sem ter que fazer nada. Então você sempre tem essa informação na ponta dos dedos.

Ok, então podemos usar este relatório e em uma execução real, e ter uma sensação de quantos CPUs e quanta memória está sendo usada pelo nosso pipeline e voltar e colocar aqueles valores de volta em nosso arquivo de configuração, para que da próxima vez talvez não solicitemos tanto. E podemos ser um pouco mais enxutos.

Agora você pode ficar realmente esperto sobre configurar arquivos de configuração de pipeline. E novamente, se você está usando Seqera Platform, procure por um pequeno botão que parece uma lâmpada. Porque se você clicar nisso, vai gerar um arquivo de configuração altamente otimizado, que é adaptado especificamente para seus dados, sua execução e seu pipeline. Para executá-lo da maneira mais eficiente possível.

Mas por enquanto, vou dizer que na verdade o número padrão de CPUs que o Nextflow estava dando estava bom e mas só preciso de um gigabyte de memória.

## 5.3. Defina alocações de recursos para um processo específico

Agora, na vida real, é bastante incomum que todos os processos em seu pipeline vão precisar dos mesmos requisitos. Você pode ter algo como MultiQC como uma ferramenta de relatório, que precisa muito pouco em termos de recursos e executa bastante rapidamente.

E então talvez você tenha algo que está indexando um genoma de referência ou fazendo algum alinhamento ou fazendo algum outro job. Não importa o que é, que leva muitos recursos. E então para essas diferentes submissões de job para um agendador, você quer dar diferentes quantidades de recursos.

Sob este escopo process, podemos definir uma configuração, que mira processos específicos de diferentes maneiras.

Aqui estamos usando withName, também podemos usar labels, e estes podem usar um padrão para mirar um ou múltiplos processos. Aqui estamos apenas dizendo quaisquer processos que tenham um nome cowpy definido para dois gigabytes de memória e dois CPUs, e porque este é um seletor mais específico do que o de nível superior process, este é sobrescrito nestes casos, então você pode construir um arquivo de configuração legal aqui, que realmente adapta todos os seus diferentes processos em seu pipeline para torná-los realmente eficientes.

## 5.5. Adicione limites de recursos

Agora como desenvolvedor de pipeline, eu provavelmente conheço as ferramentas muito bem, e quero que tudo execute o mais rápido e o melhor possível. Então pode ser que eu coloque números bem altos para alguns destes porque sei que vai executar muito mais rápido se eu der ao cowpy 20 CPUs em vez de dois.

Isso está bem até você ir executar no seu laptop ou no GitHub Actions Continuous Integration test, ou algum outro sistema, que talvez não tenha 20 CPUs disponíveis.

Agora quando você tenta executar o pipeline, vai travar porque o Nextflow vai dizer, não posso submeter este job em lugar nenhum. Não tenho os recursos disponíveis.

Agora para evitar aquela falha dura, podemos adicionar um pouco mais de configuração, que é específica para nosso sistema agora, chamada limites de recursos. E isso parece assim. Está sob o escopo process novamente.

E limites de recursos, você pode especificar basicamente o teto do que você tem disponível. É um mapa aqui, e você pode, dentro deste mapa, você pode definir a memória, os CPUs, e o tempo.

Agora o que acontece é quando o Nextflow submete uma tarefa de um processo, ele olha o que é solicitado e ele basicamente apenas faz um mínimo entre isso e aquilo. Então se solicitamos 20 CPUs, mas apenas quatro estão disponíveis, vai solicitar quatro. O pipeline não trava e usa o mais próximo possível do que foi projetado pelo desenvolvedor do pipeline.

## 6. Use perfis para alternar entre configurações predefinidas

Ok. Eu disse que os limites de recursos aqui podem ser específicos do sistema, e talvez eu tenha um arquivo Nextflow config no meu pipeline, e sei que as pessoas vão estar usando isso em uma variedade de lugares diferentes. Agora, em vez de forçar todos a criar seu próprio arquivo Nextflow config toda vez, o que posso fazer é posso agrupar diferentes predefinições de configuração juntas em perfis de configuração.

Vou rolar um pouco para baixo aqui e então logo após params, porque a ordem do arquivo de configuração aqui é importante, o arquivo de configuração é carregado sequencialmente, então vou colocar esses perfis depois de tudo mais para que sobrescreva os params previamente definidos. E vou colar esses perfis do material de treinamento.

Então há um novo escopo de nível superior chamado profiles. Podemos ter nomes arbitrários aqui. Então temos my_laptop e univ_hpc. E aqui podemos ver que estamos definindo os outros mesmos parâmetros de configuração que estávamos antes. Agora dentro apenas de um perfil. Então temos um executor local para executar no my_laptop e estou submetendo para um cluster SLURM no HPC.

Estou usando Docker localmente, conda no HPC, e o sistema HPC tem limites de recursos muito mais altos.

Agora posso executar o pipeline com a opção CLI -profile, dizer qual perfil quero usar. Então vou usar my_laptop, e o Nextflow vai aplicar toda a configuração dentro daquele escopo de perfil. Então posso tentar isso agora. É o mesmo comando de antes. Nextflow run hello config, e faço dash profile, traço único porque é a opção central do Nextflow, dash profile my_laptop.

Agora vai aplicar em lote toda aquela opção de configuração. Oh, e você pode ver, eu disse antes que isso pode acontecer que o requisito do processo, pediu quatro CPUs e eu só tenho dois nesta instância Codespaces.

Então esta é uma boa oportunidade apenas para testar os limites de recursos do processo, e dizer que só tenho dois CPUs no my_laptop, ou neste Codespaces. Agora se executarmos novamente, deve limitar aquele requisito a dois e esperançosamente o pipeline vai executar. Ótimo.

## 6.2. Crie um perfil de parâmetros de teste

Note que esses perfis não têm que ter apenas configuração sobre sua infraestrutura. Você pode ter agrupamentos de qualquer configuração aqui, incluindo parâmetros.

Então outra coisa que você verá muito frequentemente nos pipelines das pessoas é um perfil de teste, que inclui parâmetros, que você normalmente submeteria por usuário. Mas aqui temos, basicamente diferentes padrões sensatos para quando quero executar casos de teste.

E isso é ótimo porque não tenho que necessariamente ir e especificar todas essas coisas, que podem ser parâmetros obrigatórios. Caso contrário posso apenas dizer dash profile test e vai apenas executar pronto para uso.

Agora algo a notar é que perfis também podem ser combinados mais de um. Então posso fazer profile my_laptop aqui, e então também adicionar test. Não faço profile duas vezes. Apenas faço uma lista separada por vírgulas aqui sem espaços. E vai aplicar esses perfis em ordem. Então vai pegar a configuração do perfil my_laptop, e então vai aplicar a configuração test por cima.

Realmente conveniente e você pode ver como você pode configurar muitos grupos padrão sensatos aqui para tornar fácil executar seu pipeline.

## 6.3. Use nextflow config para ver a configuração resolvida

Esperançosamente, eu te convenci de que a resolução de configuração do Nextflow é poderosa, mas eu não te culparia se você estiver ficando um pouco vesgo neste ponto depois que eu disse cerca de 20 maneiras diferentes de fornecer configuração e dar todas essas diferentes camadas como uma casca de cebola.

Então se alguma vez você estiver se sentindo inseguro sobre qual é a configuração resolvida final para o Nextflow, saiba que há um comando chamado "nextflow config", e podemos executar isso e ele vai nos dizer qual é a configuração resolvida em nossa localização atual.

Então quando eu executo aqui, ele encontra o arquivo "nextflow.config" no diretório de trabalho atual, e processa toda a configuração diferente, e me dá a saída resolvida.

Note que o arquivo de configuração Nextflow também pode aceitar a opção CLI profile. Então se eu disser para resolver nos perfis my_laptop e test, e você pode ver que também aplicou os limites de recursos aqui da opção de configuração my_laptop e também definiu os params, que estavam no test.

Então esta é uma maneira legal apenas de explorar como a resolução de configuração está funcionando, se você estiver de alguma forma inseguro.

## Encerramento

Ok, é isso. Isso é configuração Nextflow em poucas palavras. Você pode fazer muita coisa com configuração. É realmente poderoso. Mas estes são a maioria dos casos de uso comuns que você vai se encontrar fazendo, e esses conceitos se aplicam a todas as diferentes opções.

Dê um tapinha nas costas porque este é o fim do curso de treinamento Hello Nextflow. Você está esperançosamente agora confiante tanto em escrever seu próprio pipeline Nextflow do zero, configurá-lo e executá-lo, e você conhece todos os detalhes e as coisas para ficar atento.

Há mais um quiz que você pode tentar na página de treinamento de configuração. Então desça e tente isso e certifique-se de que entendeu todas essas partes sobre a configuração.

E junte-se a nós no último vídeo apenas para um rápido encerramento sobre alguns dos próximos passos que podem ser bons de fazer após este curso de treinamento.

Obrigado por ficar conosco. Muito bem e te vejo no próximo vídeo.
