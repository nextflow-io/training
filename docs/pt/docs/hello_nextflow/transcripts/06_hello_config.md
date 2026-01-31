# Parte 6: Hello Config - Transcrição

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Tradução assistida por IA - [saiba mais e sugira melhorias](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

<div class="video-wrapper">
  <iframe width="560" height="315" src="https://www.youtube.com/embed/IuDO2HeKvXk?si=tnXTi6mRkITY0zW_&amp;list=PLPZ8WHdZGxmXiHf8B26oB_fTfoKQdhlik" title="YouTube video player" frameborder="0" allow="accelerometer; autoplay; clipboard-write; encrypted-media; gyroscope; picture-in-picture; web-share" referrerpolicy="strict-origin-when-cross-origin" allowfullscreen></iframe>
</div>

!!!note "Notas importantes"

    Esta página mostra apenas a transcrição. Para instruções completas passo a passo, retorne ao [material do curso](../06_hello_config.md).

    Os números das seções mostrados na transcrição são fornecidos apenas para fins indicativos e podem não incluir todos os números de seção nos materiais.

## Boas-vindas

Olá, bem-vindo à parte seis do curso de treinamento Hello Nextflow.

Este capítulo se chama Hello Config, e é a parte final do nosso curso de treinamento.

Neste capítulo, vamos falar sobre configuração do Nextflow. A configuração do Nextflow é muito poderosa. Ela nos permite executar o mesmo pipeline em múltiplas infraestruturas computacionais diferentes com diferentes provisionamentos de software e diferentes opções no próprio pipeline.

Isso significa que você pode pegar pipelines Nextflow construídos por outras pessoas e executá-los no seu sistema, mesmo que eles tenham sido construídos para uma infraestrutura completamente diferente. Esta capacidade de configurar o Nextflow torna os fluxos de trabalho verdadeiramente portáteis e compartilháveis.

Neste capítulo, usaremos o fluxo de trabalho que construímos nas partes anteriores, mas não vamos editar o código do fluxo de trabalho. Vamos apenas olhar nosso arquivo de configuração do Nextflow e ver como alterar a configuração modifica a forma como o Nextflow é executado.

Ok, vamos começar.

Assim como antes, vamos começar indo para training.nextflow.io. Vá para o lado esquerdo em Hello Nextflow e capítulo seis. Hello config. Agora vou entrar no meu ambiente GitHub Codespaces e verificar o script que usaremos.

## 0. Aquecimento: Verificar que o Docker está habilitado e executar o fluxo de trabalho Hello Config

Este se chama Hello Config, e está começando de onde estávamos antes. Então parecendo exatamente o mesmo com nossos três parâmetros. Greetings para o arquivo CSV, batch para o nome do lote de saída e character para o nome do cowpy. Temos nossas quatro importações dos diferentes processos, e então temos um fluxo de trabalho onde os encadeamos.

Na verdade, vou fechar este arquivo agora porque não vamos tocar no arquivo Nextflow em absoluto neste capítulo. Vamos trabalhar puramente dentro do arquivo de configuração. Se eu olhar para o arquivo Nextflow dot config que vimos brevemente no capítulo cinco anterior, podemos ver que temos uma única instrução aqui: Docker enabled equals true, que está dizendo ao Nextflow para usar Docker quando executar este fluxo de trabalho.

Estou usando Nextflow dot config na raiz do pipeline aqui, que é carregado automaticamente quando executo o Nextflow. Mas lembre-se, o Nextflow pode carregar arquivos de configuração de múltiplos lugares.

Se eu verificar com Nextflow docs vou para Configuration, você pode ver uma lista desses lugares e uma prioridade na qual eles são carregados.

Ok. Vamos verificar se nosso fluxo de trabalho está executando como esperamos. Abrir um terminal. Fazer Nextflow. Run. Hello, config. E pressionar enter. Devemos ter esses quatro processos executando, terminando com um comando cowpy. Com certeza, isso funcionou corretamente. Eu tinha o Docker habilitado, puxou o Docker e executou o cowpy para mim, exatamente como fez no final do capítulo cinco.

## 1. Determinar qual tecnologia de empacotamento de software usar

Ok. Digamos que estou executando em um HPC e não tenho o Docker instalado. A melhor coisa a fazer neste cenário seria usar Singularity ou Apptainer. Se eu fosse fazer isso, iria para o módulo cowpy e mudaria este contêiner para usar a imagem singularity como mostrei no capítulo anterior, com um oras://, que você também pode obter do Seqera Containers.

Eu então iria para Nextflow dot config definir Docker enabled para false e fazer singularity enabled equals true. Ou, se usar Apptainer, apptainer enabled equals true e isso funcionaria.

O Nextflow suporta outras tecnologias também, além de contêineres, algo com que você pode estar familiarizado é conda. Aqui podemos fazer conda enabled equals true e definir Docker para false. conda não usa a mesma diretiva container. Em vez disso, podemos adicionar uma nova aqui chamada conda. Então especificamos o pacote conda que queremos usar. É uma boa prática ser o mais específico possível para tentar tornar o pipeline o mais reproduzível possível. Então vou especificar o canal conda, conda-forge, e então cowpy, e a versão exata, que era 1.1.5.

Eu também poderia apenas escrever cowpy se quisesse, mas isso poderia resolver para uma versão diferente do cowpy em diferentes execuções do pipeline.

O legal sobre isso é que não toquei na diretiva docker de forma alguma. Esta imagem Docker ainda está lá. Estou apenas fornecendo duas alternativas agora, e estas podem ser ligadas ou desligadas usando apenas um arquivo de configuração.

## 1.3. Executar o fluxo de trabalho para verificar que ele pode usar Conda

Conda está agora habilitado, então vamos experimentar.

Ótimo. Está executando e você pode ver que há uma mensagem do Nextflow aqui dizendo que o Nextflow está criando um ambiente conda para mim, e está usando este local de cache.

Nos bastidores, o Nextflow está executando comandos "conda create" para mim para criar um novo ambiente conda isolado com apenas os pacotes que eu quero, e então instalando e buscando esses pacotes conda para que possa executar o processo.

Você pode ver que levou um pouco de tempo porque estava criando o ambiente e instalando o software pela primeira vez. No entanto, ele armazenou em cache este ambiente, então se eu executar o mesmo comando Nextflow novamente, deve ser muito mais rápido porque ele reutilizará o mesmo ambiente conda.

Uma das coisas legais sobre isso é que essas diretivas podem ser especificadas no nível do processo, não apenas para o fluxo de trabalho inteiro. Então, se você quiser, pode misturar e combinar qual tecnologia é usada para diferentes processos.

## 2. Alocar recursos computacionais com diretivas de processo

O arquivo de configuração do Nextflow pode fazer muito mais do que apenas empacotamento de software. Também podemos dizer ao Nextflow como realmente executar os passos no pipeline. Um exemplo é dizer a um sistema host quais recursos devem ser disponibilizados para cada tarefa em execução.

Por padrão, o Nextflow não dá muito. Ele dá uma única CPU e apenas dois gigabytes de memória para cada processo.

Isso é provavelmente algo que gostaríamos de mudar, para que processos que demoram muito para executar possam ter mais recursos e executar mais rapidamente, mas pode ser difícil saber o que alocar para um processo. O Nextflow tem alguns truques legais para ajudá-lo com isso.

## 2.1. Executar o fluxo de trabalho para gerar um relatório de utilização de recursos

Vamos executar o fluxo de trabalho novamente. Desta vez, vou adicionar um argumento adicional, que é dash with reports. É uma opção central do Nextflow, então é um único hífen. E então qualquer nome de arquivo que eu goste. Neste caso, vou chamá-lo de report config one html.

Vou executar o fluxo de trabalho novamente. Vai executar exatamente como antes, mas vai me dar um relatório auxiliar adicional, que você pode ver que apareceu aqui na barra lateral.

Vou clicar com o botão direito neste arquivo, clicar em download, que o baixa do GitHub Codespaces para o meu sistema local, para que eu possa visualizá-lo facilmente no navegador da web aqui em cima.

Este relatório pode ser gerado para qualquer execução do Nextflow, e tem muita informação. Começa no topo com alguns metadados sobre qual comando foi usado, quando o fluxo de trabalho foi executado, quanto tempo levou, mas à medida que você rola para baixo, obtemos informações mais detalhadas sobre os recursos que foram usados por cada passo no pipeline.

Como cada processo é executado várias vezes para diferentes tarefas, temos um gráfico de caixa mostrando a variação dos recursos que usamos para cada processo.

Se eu rolar um pouco mais para baixo, vejo informações similares sobre memória usada e duração do job. Também leitura e escrita de disco.

Você pode imaginar que para um pipeline grande com tarefas de longa execução, isso pode ser muito informativo sobre como ajustar finamente a configuração dos recursos que você está solicitando para que você não solicite em excesso, mas também para que você possa fornecer o suficiente para que execute rapidamente.

Se eu continuar rolando o relatório para baixo, também vemos uma tabela de tarefas, que nos mostra informações detalhadas sobre cada tarefa individual que foi executada no fluxo de trabalho. Isso inclui informações como o script resolvido, que foi executado.

Ok, vamos voltar ao nosso arquivo de configuração. Vimos que realmente não precisávamos de muito para nosso fluxo de trabalho, então vamos dizer ao Nextflow que só precisamos de um gigabyte de memória para cada processo no fluxo de trabalho.

Agora, quando definimos assim no nível de processo, isso é aplicado a cada processo individual no pipeline.

## 2.3. Definir alocações de recursos para um processo individual

Para fins de argumento, vamos fingir que cowpy está realmente fazendo muito trabalho pesado e precisa de mais recursos do que as outras tarefas. Podemos definir um bloco extra de configuração aqui, que se aplica apenas a esse processo usando, with name cowpy.

Isso é chamado de seletor de configuração, e podemos definir diferentes padrões aqui para corresponder a diferentes processos. Por exemplo, eu poderia fazer cow star. Então eu sigo isso com algumas chaves e vamos dar dois gigabytes de memória em vez de um e vamos dizer duas CPUs.

Agora o Nextflow estará dando a cada processo no fluxo de trabalho um gigabyte, exceto por esta solicitação, que é mais específica. Então ela sobrescreve. E apenas para quaisquer processos que são chamados cowpy, receberão dois gigs de memória e duas CPUs.

Note que o Nextflow é inteligente sobre a utilização de recursos. Então, se você começar a colocar esses números em valores mais altos, verá que o Nextflow começa a enfileirar submissões de jobs uma após a outra, em vez de executar todas em paralelo, para que não solicite em excesso os recursos que estão disponíveis.

## 2.4. Executar o fluxo de trabalho com a configuração modificada

Vamos tentar executar o fluxo de trabalho novamente e vamos salvar um novo relatório desta vez.

Ok, podemos baixar este arquivo e dar uma olhada.

Sim, sem surpresas, parece basicamente exatamente o mesmo porque este é um fluxo de trabalho fictício, que não está fazendo nada real. Mas você pode imaginar como esta abordagem iterativa de definir limites e fazer fluxos de trabalho da vida real com este tipo de relatório permite que você faça uma abordagem baseada em evidências para definir configuração apropriada e realmente aproveitar ao máximo os recursos computacionais que você tem disponíveis.

Você pode começar a ser realmente inteligente sobre isso. O Nextflow tem uma capacidade embutida para tentar novamente falhas, e você pode aproveitar em seu arquivo de configuração usando um closure como este e definindo dinamicamente os recursos que são disponibilizados. Então aqui eu disse ao Nextflow para multiplicar aqueles dois gigabytes pela tentativa de retry. Então a segunda tentativa receberá quatro gigs, a terceira tentativa receberá seis gigs e assim por diante. Isso está um pouco além do escopo deste curso de treinamento, mas se você estiver interessado, confira a documentação do Nextflow, que tem uma seção legal sobre lógica de retry dinâmica.

## 2.5. Adicionar limites de recursos

Agora, uma coisa que você pode notar sobre isso é que esse tipo de coisa pode tornar muito fácil acidentalmente ir além dos recursos disponíveis no seu sistema. Se você solicitar mais recursos do que estão disponíveis, o Nextflow lançará um erro sobre sua configuração e interromperá a execução. Para evitar isso, você pode usar algo chamado limites de recursos.

No escopo de processo, em nosso fluxo de trabalho, podemos definir limites de recursos assim, que recebe um array, e podemos especificar a memória máxima, CPUs e tempo que estão disponíveis neste sistema.

Definir valores altos aqui não aumenta a quantidade de recursos que são solicitados. Ainda vamos usar um gigabyte em nossas solicitações, mas significa que se qualquer uma dessas solicitações chegar a 750, elas atingirão esse teto e nada mais do que isso será solicitado, o que significa que o Nextflow continuará a executar e não travará por causa de recursos indisponíveis.

Então, esta é uma boa proteção para usar, especialmente se você estiver usando lógica dinâmica com sua alocação de recursos.

A outra situação onde isso é realmente útil é se você estiver usando pipelines que são públicos e não controlados por você. Eles podem vir com padrões de configuração, e o Nextflow automaticamente tomará a abordagem correta de limitar quaisquer solicitações de recursos para executar em seu sistema.

Ok, ótimo. Falamos sobre software. Falamos sobre alocação de recursos, e descrevemos diferentes escopos de configuração, tanto para todos os processos quanto para processos específicos.

## 3. Usar um arquivo de parâmetros para armazenar parâmetros do fluxo de trabalho

Ok, a seguir vamos voltar nossa atenção para parâmetros. Podemos definir parâmetros no arquivo de configuração assim como fizemos antes no script Nextflow. Então params dot greeting equals hello ou ou usar escopo params e definir foo equals bar.

E isso é ótimo para definir padrões para seu fluxo de trabalho. No entanto, quando você está executando pipelines, pode ser legal especificar parâmetros em um arquivo JSON ou YAML.

Usar um arquivo como este é muito melhor do que especificar opções de linha de comando com dash dash. Pois quando você executa um fluxo de trabalho, você pode ter que especificar muitos parâmetros e pode ser tedioso escrevê-los todos em uma única CLI e propenso a erros. Além disso, é improvável que você se lembre de todos os parâmetros que usou, então se você codificar isso em um arquivo, é mais fácil lançar o fluxo de trabalho novamente, usando os mesmos parâmetros no futuro.

Temos um arquivo de exemplo aqui chamado test params, e você pode ver que isso especifica os três parâmetros que temos em nosso fluxo de trabalho com três valores diferentes. Pessoalmente, acho YAML mais fácil de escrever do que JSON. Então, apenas para demonstrar que funciona, vou criar um novo arquivo chamado Test yaml e copiar estes, me livrar das aspas. E salvar.

Esses arquivos JSON e YAML podem ser mais fáceis de escrever, pois são sintaxe mais familiar. Mas note que estes são apenas para parâmetros e eles só aceitam sintaxe de chave-valor assim.

## 3.1. Executar o fluxo de trabalho usando um arquivo de parâmetros

Vamos experimentar. Fazer o mesmo comando de antes. Livrar-se do relatório e vou fazer dash params file test params yaml.

Não, esta é uma opção central do Nextflow, então é um único hífen.

Ok. Executou o fluxo de trabalho e usou os parâmetros naquele arquivo YAML em vez de eu especificá-los todos na linha de comando. Pode parecer exagero apenas para este exemplo simples, mas você pode imaginar se você tem 10 ou 20 parâmetros diferentes, pode ser um incômodo digitar manualmente, e isso é apenas muito mais fácil de editar em um editor de código e manter para fins de reprodutibilidade.

## 3. Determinar qual(is) executor(es) deve(m) ser usado(s) para fazer o trabalho

Ok. Falamos sobre empacotamento de software com Docker e conda. Falamos sobre requisitos de recursos de processo com CPUs e memória. E falamos um pouco sobre como especificar parâmetros ao executar fluxos de trabalho.

As partes finais da configuração realmente são a execução, a infraestrutura computacional subjacente em si, e esta é a verdadeira joia da coroa do Nextflow: que podemos executar esses mesmos fluxos de trabalho em múltiplas infraestruturas computacionais diferentes.

Na verdade, vou mudar para o material de treinamento escrito por um segundo. Nesta parte do treinamento, podemos ver alguns exemplos diferentes de como diferentes executores, neste caso, escalonadores HPC, definem os requisitos de recursos necessários para submeter um job.

Então para Slurm, você tem esses cabeçalhos SBATCH, que definem dash dash mem e o número de CPU. Se você está usando PBS, você tem cabeçalhos diferentes, e se você usa Grid Engine, você tem cabeçalhos diferentes novamente.

Você pode imaginar que é ainda mais diferente se você quiser executar na nuvem, seja AWS batch, Google Cloud, Azure, ou mais.

Cada uma dessas infraestruturas computacionais subjacentes é chamada de executor e o Nextflow sabe como falar com todos esses diferentes executores para submeter jobs com a sintaxe correta.

A boa notícia é que você não precisa saber sobre isso. Tudo que você tem que fazer é dizer ao Nextflow, qual executor usar.

## 3.1. Direcionando para um backend diferente

Voltamos ao nosso arquivo de configuração e ao processo fazemos executor, e vou digitar local.

Local é na verdade o padrão, se você não especificar nenhum outro executor, local é o que será usado, e isso apenas significa seu sistema host, onde quer que você tenha lançado o Nextflow,

Eu poderia especificar em vez disso, Slurm. E isso submeteria jobs Slurm, ou eu poderia dizer AWS batch, e isso submeteria jobs para AWS batch.

Você precisa de alguma configuração adicional em alguns casos, por exemplo, executar na nuvem precisará de certas credenciais, mas realmente este é o núcleo disso, e pode ser tão simples quanto uma ou duas linhas de configuração para executar seu fluxo de trabalho em um ambiente computacional completamente diferente.

Mesmo que estejamos executando em um sistema simples dentro do Codespaces, ainda posso brincar com isso um pouco e fingir que estamos executando no Slurm. Se eu então lançar o fluxo de trabalho novamente, Nextflow run, hello config. Vai falhar porque não será capaz de submeter jobs para o Slurm. Mas ainda podemos ir para os diretórios de trabalho e ver o que o Nextflow fez. Então, se formos para este diretório de trabalho e olharmos para Command Run. Você pode ver no topo deste arquivo, agora temos essas linhas de cabeçalho sbatch, que tentaram especificar os recursos necessários para o job Slurm.

## 4. Usar perfis para selecionar configurações predefinidas

Ok, estamos quase lá. A parte final deste capítulo é falar sobre perfis de configuração. Se você está executando seu pipeline em vários sistemas diferentes, pode ser irritante ter todos esses arquivos de configuração Nextflow diferentes, que você precisa especificar toda vez.

Em vez disso, você pode codificar agrupamentos de configuração dentro do seu arquivo Nextflow config, e ligar e desligar esses grupos usando uma flag de perfil. Vamos ver como isso se parece.

## 4.1. Criar perfis para alternar entre desenvolvimento local e execução em HPC

Vamos criar dois perfis em nosso exemplo aqui, um para o meu laptop e um para um sistema HPC mais pesado. Vou trapacear um pouco e apenas copiar o código do material de treinamento e colocá-lo aqui.

Temos um novo escopo chamado profiles, e então temos um nome para cada perfil, que pode ser qualquer coisa. E dentro disso temos configuração, que parece exatamente o mesmo que a configuração de nível superior que já escrevemos. Então, novamente, temos escopo process. Escopo Docker.

No perfil chamado my laptop. Estou dizendo para executar usando o executor local, então no meu sistema host e usar Docker.

No perfil university HPC aqui estou dizendo para usar Slurm para submeter jobs, usar conda em vez de Docker, e estou especificando diferentes limites de recursos, que podem corresponder ao tamanho do sistema de nós no HPC que estou usando.

Por padrão, nenhuma dessas configurações será usada quando eu executar o Nextflow, tenho que especificar que quero usar um desses perfis.

## 4.2. Executar o fluxo de trabalho com um perfil

Vamos fazer nextflow run hello config. E vou fazer dash profile, único hífen porque é uma opção central do Nextflow. E então o nome que dei a ele, que é my laptop. O Nextflow deve agora usar o bloco de configuração que foi especificado dentro daquele perfil de configuração, e aplicá-lo quando executar o Nextflow. Se eu quisesse usar o outro bloco de configuração, só tenho que mudar esse nome de perfil. Muito mais fácil de lembrar. Muito mais fácil de usar.

## 4.3. Criar um perfil de teste

Note, os perfis podem ter qualquer tipo de configuração, então não precisa estar relacionado ao seu ambiente de execução. Por exemplo, vamos criar um novo perfil aqui, que tem um conjunto de parâmetros. Podemos mudar isso para tux e mudar para my profile, e agora quando fazemos profile test, vai especificar esses parâmetros, que sobrescreverão os parâmetros que são especificados no nível superior do fluxo de trabalho.

Quando você executa o Nextflow, você pode encadear múltiplos perfis e eles serão aplicados em sequência.

## 4.4. Executar o fluxo de trabalho localmente com o perfil de teste

Então posso pegar o comando anterior e fazer vírgula test. Isso aplicará o, my laptop config primeiro, e então aplicará o test config. Se houver alguma sobreposição, então o perfil à direita sobrescreverá qualquer configuração em perfis anteriores. Se eu pressionar enter, vamos ver o que acontece.

Ok, temos um novo arquivo de resultados aqui. Você pode ver o My Profile, que eu especifiquei como uma das opções. E também podemos ver cowpy, my profile, e com certeza, está o tux. Então funcionou.

## Conclusão

Ok! Incrível. É isso. Você chegou ao final do curso. Você ganha um pouco de confete de celebração. Parabéns por terminar este capítulo.

[Próxima transcrição do vídeo :octicons-arrow-right-24:](07_next_steps.md)
