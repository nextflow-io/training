# Parte 5: Hello Containers - Transcrição do Vídeo

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Tradução assistida por IA - [saiba mais e sugira melhorias](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

<div class="video-wrapper">
  <iframe width="560" height="315" src="https://www.youtube.com/embed/Xqr--bKEN9U?si=QinuAnFwFj-Z8CrO&amp;list=PLPZ8WHdZGxmWKozQuzr27jyMGqp9kElVK&amp;cc_load_policy=1" title="YouTube video player" frameborder="0" allow="accelerometer; autoplay; clipboard-write; encrypted-media; gyroscope; picture-in-picture; web-share" referrerpolicy="strict-origin-when-cross-origin" allowfullscreen></iframe>
</div>

!!!note "Notas importantes"

    Esta página mostra apenas a transcrição. Para instruções completas passo a passo, retorne ao [material do curso](../05_hello_containers.md).

    Os números de seção mostrados na transcrição são fornecidos apenas para fins indicativos e podem não incluir todos os números de seção dos materiais.

## Boas-vindas e contexto

Olá, e bem-vindo de volta ao Hello Nextflow. Esta é a parte cinco chamada Hello Containers. E nesta parte do curso, vamos falar tudo sobre como encapsular os requisitos de software para um pipeline, para que as pessoas que executam o pipeline não precisem pensar em instalar o software.

Se você trabalha em bioinformática há tanto tempo quanto eu, você pode se lembrar do que eu costumo chamar de "os velhos tempos ruins", onde quando você queria executar o pipeline de outra pessoa ou replicar o trabalho dela, você passava horas ou dias tentando instalar todas as diferentes ferramentas de software que elas usaram, nas mesmas versões, tentando compilá-las na sua máquina, e era um pesadelo. Era realmente difícil.

Se você estava executando em um HPC, você pode ter usado módulos de ambiente onde os administradores de sistema tentavam instalar software para você, o que era aceitável, mas ainda imperfeito.

Mas agora temos maneiras melhores de fazer isso. O Nextflow tem suporte integrado para diferentes tecnologias de contêineres de software. Docker é a mais comum. É essa que vamos usar hoje. Funciona bem no Codespaces. Funciona bem no seu computador local e funciona bem na nuvem.

Mas também Singularity ou Apptainer, que são muito comuns em sistemas HPC e efetivamente funcionam exatamente da mesma maneira. Ou Podman, Shifter, há várias outras que são todas muito similares.

A única extra, que é meio similar mas não exatamente, que o Nextflow suporta, é o Conda. E o Nextflow pode gerenciar ambientes Conda para você em uma base por processo, o que é muito melhor do que fazer seus próprios ambientes Conda. E novamente, pode ser enviado com um pipeline.

Vamos começar este capítulo falando um pouco sobre tecnologias de contêineres e Docker e como eles funcionam. E vamos fazer a primeira metade apenas manualmente no Docker para que você entenda o que está acontecendo por baixo dos panos e como isso funciona. Porque isso é realmente importante para entender o que o Nextflow está fazendo e como entender o que seu fluxo de trabalho está fazendo quando está sendo executado.

Então. Vamos pular para nosso Codespaces. Agora eu limpei tudo novamente, mas se pularmos para Hello Containers, você deve ver que todos os nossos scripts e tudo estão lá, os mesmos do final do capítulo de módulos. Então temos nossos diferentes módulos aqui, que eu criei no diretório modules.

Eles ainda estão lá. Eles precisam estar lá para que possa executar. E o fluxo de trabalho e a saída são todos iguais, exceto que mudamos o caminho de publicação da saída para Hello Containers, para que seus arquivos acabem nesse diretório.

Podemos executar isso agora para verificar se funciona, se você quiser, ou podemos seguir em frente com o terminal.

## 1. Usar um contêiner 'manualmente'

Vamos usar o Docker para gerenciar nossos contêineres, e posso verificar se ele está instalado no meu Codespaces fazendo "docker -v", que me mostra a versão que está instalada e tudo, e que está funcionando corretamente.

Agora contêineres e Docker têm dois conceitos que são realmente importantes. Um é chamado de imagem, e outro é chamado de contêiner. A imagem é o instantâneo, por assim dizer, de todo o sistema de arquivos que você estará usando, e o contêiner é o ambiente em execução. Então você cria um contêiner usando uma imagem.

Uma vez que você está dentro desse contêiner, ele normalmente funciona como um sistema operacional completo. Ele é isolado do mundo externo. Ele é separado de todo o resto, e isso é uma coisa boa. É assim que conseguimos uma reprodutibilidade tão boa com o Nextflow.

Porque para tarefas executadas dentro de um contêiner, elas não são afetadas por quaisquer arquivos de configuração no seu sistema local. Quaisquer outras influências externas, elas são executadas em sua própria pequena caixa de areia. Os arquivos são então produzidos de uma maneira muito, muito reproduzível porque você está usando as mesmas bibliotecas subjacentes, todas as mesmas dependências, exatamente o mesmo software para cada pessoa executando em cada ambiente de computação diferente. O que francamente eu acho fantástico e incrível que funcione. E mesmo, mesmo até hoje ainda meio que me deixa impressionado que isso seja possível.

## 1.1. Baixar a imagem do contêiner

Então vamos experimentar usar algumas imagens Docker e Docker, quando você o executa no seu sistema, tem um registro docker no seu computador, ou neste caso, no code space, que mantém o controle de todas as diferentes imagens que foram baixadas e usadas no passado, e as diferentes camadas das quais elas são construídas.

Podemos ver quais imagens temos localmente com Docker fazendo "docker image ls". E neste caso você pode ver que há várias imagens Docker aqui, que todas têm a ver com a configuração deste Codespaces. Todas têm a ver com contêineres de desenvolvimento e coisas assim. Então você não precisa se preocupar muito com elas, mas à medida que adicionamos mais imagens e as baixamos, conforme este curso avança, você pode verificar essa lista e verá que o registro local está mantendo o controle de todas essas coisas que baixamos.

Mas vamos pegar uma nova fazendo "docker pull". E isso diz ao Docker para buscar uma nova imagem da web.

Em seguida, colocamos o URI para esse contêiner. Agora isso poderia ser uma imagem docker que você construiu localmente e depois enviou para a internet. Poderia ser uma imagem que outra pessoa fez. Há muitas, muitas, muitas maneiras diferentes de fazer imagens Docker, mas sem dúvida uma das maneiras mais simples é terceirizar isso, e fazer com que outra pessoa faça isso para você.

E o que vamos usar neste tutorial é um serviço da Seqera chamado Seqera Containers.

Agora, o Seqera Containers é totalmente gratuito, e usa um software de código aberto que desenvolvemos chamado Wave, que foi construído para gerenciar contêineres de uma maneira complementar ao Nextflow. E ele lida com muitos dos casos de uso comuns que nos encontramos lidando, com Nextflow.

É muito comum que o software que precisamos esteja empacotado no Conda, nos canais Bioconda, ou conda-forge ou outros canais mais específicos de domínio. E Wave e Seqera Containers são muito bons em construir imagens a partir disso.

Então posso ir a esta interface web e vamos mexer com o pacote chamado "cowpy". Então eu digito o nome do pacote que quero. Ele procura, encontrou no índice de pacotes Python, então posso usar isso. Ou se eu esperar um pouco mais, ele está procurando bioconda e conda-forge. E você pode ver, posso especificar qualquer canal conda aqui. Então, se você quiser encontrar um canal Nvidia ou qualquer outra coisa, isso também deve funcionar.

E então posso especificar se quero que ele construa uma imagem docker para mim ou uma imagem singularity e também qual arquitetura de CPU eu quero usar. Então amd64 ou arm64.

E uma vez que os resultados do bioconda estão listados, posso agora ver todas as diferentes versões que estão disponíveis também. Vou colocar isso. E agora eu poderia continuar procurando e obter mais pacotes do Conda se eu quiser e compor este contêiner como eu quiser, mas só quero esse. Então vou clicar em Get Container.

Agora, outra pessoa já solicitou o mesmo contêiner antes e ele foi retornado de um registro, então nós apenas o obtemos imediatamente. Mas se ninguém mais tivesse pedido este pacote de software ou esta combinação de pacotes de software, Wave e Seqera Containers o construiriam em tempo real para nós.

Podemos copiar esta URL e também podemos ver os detalhes da compilação. E isso nos mostra o que o serviço fez no backend. Ele criou um arquivo de ambiente conda. Um arquivo docker, e então isso é ele, executando o processo de construção docker. Ele também executou uma varredura, uma varredura de segurança, então você pode ver quaisquer CVEs. E ele diz quando isso foi criado.

Wave e Seqera Containers podem fazer muito mais do que isso, mas este é meio que um caso de uso simples, que é o mais comum. E devo dizer que essas imagens são hospedadas por pelo menos cinco anos. Então você pode construir essas URLs em seus pipelines e saber que elas não vão desaparecer tão cedo.

Então eu tenho minha URL para minha imagem docker para cowpy.

Posso agora fazer "docker pull" dessa URL, e ela vai buscar todas as diferentes camadas e baixar esta imagem para que ela esteja disponível para mim localmente.

## 1.2. Usar o contêiner para executar cowpy como um comando único

Ok, agora vamos tentar realmente usá-lo. Então agora vou usar um comando "docker run" em vez de "docker pull", e vou usar a flag "--rm", que apenas diz ao Docker para desligar este contêiner uma vez que ele terminou de fazer o que pedi a ele. E então coloco o identificador para o contêiner, que é apenas um URI.

E então no final, especifico o comando que quero que o Docker execute dentro do contêiner gerado a partir desta imagem. Vou apenas dizer cowpy, que é o nome da ferramenta que está instalada do Conda Forge, que está disponível dentro da imagem.

Vou apertar enter e pronto. Executamos cowpy em um sistema. Temos uma vaquinha nos dando algumas informações.

Agora note que cowpy não está instalado no meu sistema local. Então se eu executá-lo apenas sem todo o material do Docker, ele diz, comando não encontrado. Então isso baixou uma imagem. Ele criou um contêiner usando Docker, e então entrou nesse contêiner e executou este comando para nós e nos deu a saída de volta ao nosso terminal. Muito, muito legal.

## 1.3. Usar o contêiner para executar cowpy interativamente

Ok, vamos dar um passo adiante agora e executar este contêiner interativamente e dar uma olhada por dentro, para que possamos ver o que está acontecendo dentro do contêiner.

Então se eu voltar e pegar meu comando run e vou me livrar do cowpy no final, porque na verdade não quero executar o cowpy. Quero executar um terminal Bash.

E então vou voltar aqui e vou fazer "-it", que significa Interativo e Terminal ou TTY, e vou pressionar enter.

E agora você pode ver o prompt, a parte antes de eu digitar, mudou. Este era o prompt do Codespaces onde dizia o diretório, e agora diz base e root e tmp. Então estou agora dentro do contêiner, e se eu fizer "ls", você verá que os arquivos que vejo neste diretório são diferentes dos arquivos que tenho no meu workspace.

E na verdade, não posso ver nenhum dos arquivos do meu workspace local do codespaces ou do meu drive local dentro do contêiner Docker. O tempo de execução do contêiner docker está completamente isolado e não pode escrever ou ler quaisquer arquivos de um sistema de arquivos host externo.

Posso, no entanto, ver o software que está instalado dentro do contêiner e executá-lo. Então posso executar cowpy e podemos ver um pouco mais sobre como usar cowpy. Aqui posso fazer "cowpy 'Hello World'" e isso diz, diz a ele para realmente colocar minha citação dentro de um pequeno balão de fala. E você também pode executar diferentes tipos de vacas, então não precisa ser uma vaca. Você pode fazer um "-c". E estou na Suécia, então vou escolher um alce. Muito legal. Dei a ele alguns chifres.

E há um monte de diferentes que você pode brincar, que você pode ver descritos nos documentos de treinamento.

## 1.3.4. Montar dados no contêiner

Ok. Seria bom se pudéssemos executar cowpy nos arquivos em nosso sistema de arquivos.

Claro, não é super útil apenas ter o contêiner e nenhum acesso a nada. Pode ser seguro e reproduzível, mas não é muito útil.

Então como fazemos isso? Vou sair deste contêiner Docker digitando exit, e você pode ver que o prompt nos diz que agora estamos de volta em nosso Codespaces regular novamente.

E vou executar o mesmo comando novamente. Mas desta vez vou adicionar algumas flags adicionais de volta aqui. E a importante é "-v", que significa montar um volume, que é basicamente uma parte, parte de um espaço em disco.

O "-v" pega duas partes: há como uma string e depois dois pontos e uma string. E a primeira parte é o sistema de arquivos local, que deve ser montado no contêiner. E então a segunda parte é onde isso deve acabar dentro do contêiner.

Agora eu só quero carregar todo o meu sistema de arquivos local aqui. Então "." é o diretório de trabalho atual. Então vou apenas fazer "." e depois ":", e então vamos colocar isso em um novo diretório dentro do contêiner chamado "my_project". Isso poderia realmente ser chamado de qualquer coisa.

E então vou executar novamente.

No diretório de trabalho onde sou despejado, que é /tmp, os arquivos não estão lá. Mas se eu fizer "ls my_project", lá está: todos os mesmos arquivos que tínhamos localmente em nosso Codespaces agora estão disponíveis dentro do contêiner nesse caminho.

Este é acesso de leitura e gravação, então posso criar novos arquivos neste diretório e eles aparecerão no meu sistema de arquivos host. Este diretório particular, então se comporta exatamente como se eu estivesse fora do contêiner, então posso agora ler e escrever e fazer coisas.

## 1.3.5. Usar os dados montados

Ok, vamos apenas provar que podemos fazer isso. Eu faço "cat /my_project/data/greetings.csv". Se você se lembra, o conteúdo deste arquivo parece com isso. Posso agora canalizar isso para cowpy e a vaca vai imprimir as diferentes saídas desse arquivo em seu pequeno balão de fala, o que é meio divertido.

Então você pode ver, agora podemos usar o software no contêiner para interagir com os arquivos em nosso sistema host.

Ok, vamos sair e vamos continuar com o resto do material de treinamento.

## 2. Usar contêineres no Nextflow

Então isso é muito legal usando contêineres. Espero que faça sentido. E você pode ver o valor desses contêineres e por que isso é útil para executar software de análise.

Mas como fazemos todo esse mesmo processo dentro do Nextflow? Não queremos estar executando cargas de comandos Docker nós mesmos. Queremos apenas deixar o Nextflow lidar com tudo isso para nós.

Então vamos trabalhar através disso. Vamos adicionar um novo processo ao nosso pipeline, para executar cowpy. Ok, então vamos criar um novo módulo para nosso novo processo. Então vá para modules, vamos chamá-lo de cowPy.nf, e então vou copiar o código do material de treinamento aqui.

Mas você pode ver que o processo é muito simples. Parece muito com os que fizemos até agora, temos um bloco de entrada com um caminho, que é nosso arquivo de entrada, e também um valor aqui para que este seja um caractere, então poderíamos usar um alce novamente se quisermos.

E então uma saída, que é um único arquivo aqui, um caminho e depois um script. E estamos fazendo a mesma coisa que fizemos interativamente dentro do contêiner: estamos fazendo "cat" para ler o arquivo de entrada. Estamos canalizando esse conteúdo para cowpy. Estamos escolhendo um caractere específico com base nessa entrada, estamos escrevendo para um arquivo de saída chamado cowpy input file, que é então ecoado para a saída.

Ótimo. Vamos incluir isso. Então include \{ cowpy \} from "./modules/cowpy.nf", eu o chamei de cowpy? Sim.

E então vamos chamar nosso novo processo aqui embaixo no bloco principal do fluxo de trabalho. Então vamos executar cowpy. E vamos pegar nosso novo processo cowpy e vamos dizer collectGreetings.out.

E então se você se lembra, havia duas saídas para este módulo. Uma chamada outfile e uma chamada report. A extensão VS Code está auto-sugerindo isso para nós e queremos o .outfile.

Você sempre pode entrar neste processo aqui. Você pode passar o mouse sobre ele e ele deve mostrar rapidamente quais eram as saídas. E também podemos fazer command click nele e ele abrirá o arquivo do módulo se você quiser ver com mais detalhes.

Então aqui vamos nós. Esse é o outfile lá, e esse é o caminho. Então isso agora será o arquivo de entrada para nosso processo cowpy. Fantástico.

Agora, se você se lembra, um processo cowpy tem duas entradas. Também tínhamos o canal de valor para o caractere. Então podemos adicionar "params.character" aqui. Eu poderia ter codificado isso se quisesse, mas vamos torná-lo uma opção CLI para que possamos fazer dash, dash character.

Certo. Agora preciso definir o parâmetro de entrada que acabamos de chamar e dar a ele um padrão. Então character, String. E eu gosto do alce, então vou defini-lo para moose por padrão.

Certo, vamos tentar executá-lo. Então se eu fizer Nextflow run hello containers, veremos o que acontece.

Eu poderia ter usado dash resume se tivesse os diretórios de trabalho antigos circulando por aí. E novamente, esses primeiros processos teriam sido, cacheados e teria sido um pouco mais rápido, mas deveria ser basicamente o mesmo.

Agora podemos ver imediatamente que ele lançou um erro quando chegou ao nosso novo processo, está nos dizendo aqui que houve um erro ao executar o processo cowpy e ele saiu com um status de saída 127. Este é o comando que ele tentou executar. Parece certo, parece como esperávamos. Está pegando esse nome de arquivo de saída, que parece certo, está executando com um caractere de alce e tentando salvar.

Mas você pode ver o erro do comando aqui está dizendo que o comando cowpy não foi encontrado. E isso faz sentido porque ainda não dissemos ao Nextflow para usar um contêiner. Acabamos de dar a ele o comando cowpy. E como eu disse antes, cowpy não está instalado em nosso sistema local. Então quando ele tentou executá-lo, falhou.

## 2.3.1. Especificar um contêiner para cowpy

Precisamos dizer ao Nextflow que há um contêiner disponível e ele pode usá-lo. Então como fazemos isso?

Se entrarmos em nosso módulo aqui, vamos adicionar uma nova declaração no topo chamada "container". E então vamos definir isso para uma string.

Agora, se você se lembra, lá no Seqera Containers, posso copiar essa URL e apenas soltá-la entre aspas aqui.

Agora volte e tente executá-lo novamente.

Deixe-me ver se funciona desta vez.

Infelizmente, falha exatamente da mesma maneira, mesmo agora que definimos um contêiner para o processo executar. Então, para usar nossa imagem docker, precisamos dizer ao Nextflow para habilitar o uso do Docker quando executarmos o fluxo de trabalho.

E vamos fazer isso criando um novo arquivo de configuração. Então vou dizer touch nextflow.config.

Este é um nome de arquivo especial onde se ele estiver no diretório de trabalho enquanto lanço o pipeline, ele será carregado automaticamente. Então se eu entrar neste arquivo nextflow.config, você pode ver que ele já existe, o que eu havia esquecido. E temos docker.enabled aqui já, mas está definido como falso, que é o padrão.

Então se eu mudar isso para equals True em vez disso, docker.enabled. E há documentos de referência para todos esses escopos de configuração nos documentos do Nextflow. E também você pode ver que quando passo o mouse sobre com a extensão VS Code, ela puxa os documentos específicos para isso e me diz o que significa e como configurá-lo.

Então agora definimos como verdadeiro, e se eu executar o Nextflow novamente, o Nextflow agora saberá buscar aquela imagem docker para nós se ainda não a tivermos localmente, e então executar esse processo com aquele ambiente de contêiner.

E então podemos ver que ele foi executado com sucesso e temos um pequeno tick ao lado do cowpy. Fantástico. Se eu subir e olhar no diretório de resultados, o arquivo ainda não está lá. E isso é porque ainda precisamos publicar este arquivo de saída da mesma forma que todos os outros.

Então vamos para o bloco published dentro do fluxo de trabalho, dizer mycowpy equals cowpy.out.

E então aqui embaixo no bloco de saída, mycowpy, chaves path. Ops. Hello containers. Mode, copy.

Se eu executar novamente agora, deve executar exatamente da mesma maneira. Eu poderia ter usado dash resume e eu esqueço todas as vezes. E então eu subo e agora temos um novo arquivo criado chamado cowpy-COLLECTED, e lá está meu alce dizendo BONJOUR, HELLO, HOLÀ Fantástico.

Agora é claro que eu também poderia passar agora "--character". Quais são as diferentes opções? Acho que há um Peru? Então posso usar character Turkey. Vai executar exatamente da mesma maneira. Perdi outra oportunidade de usar dash resume, e agora se carregarmos nosso arquivo e agora temos um Peru. Fantástico.

## 2.3.4. Inspecionar como o Nextflow lançou a tarefa containerizada

Ok. Última coisinha. Vamos apenas executar rapidamente este comando novamente, retomar desta vez, e dar uma olhada rápida no diretório de trabalho para ver o que é que o Nextflow está fazendo por baixo dos panos para fazer tudo isso funcionar para nós.

Desta vez é super rápido, vamos entrar neste diretório de trabalho, cd work/. Agora, se você se lembra, temos vários arquivos ponto aqui e o que nos interessa neste caso é aquele que eu disse que quase nunca precisamos olhar, chamado .command.run.

Se eu fizer code dot command run, vai abri-lo no editor. E posso pesquisar neste arquivo e se eu rolar para baixo devo ver Docker run. E você pode ver que o Nextflow está fazendo o comando docker run para nós, quando o Docker está habilitado em uma configuração. Tem um monte de diferentes flags e coisas aqui, mas você pode ver a flag "-v" que usamos nós mesmos quando estávamos executando. E você pode ver que está montando o diretório workspace local no contêiner, para que o contêiner possa acessar nossos arquivos de entrada e salvar as saídas. E então no final, ele também está executando .command.sh, que é o script gerado, que tem o comando cowpy dentro.

E então você pode ver que o Nextflow está pegando a lógica do fluxo de trabalho, que é o material que realmente nos importa, que é específico para nossa análise, e está fazendo todo o material inteligente dos bastidores para fazer o Docker funcionar em nosso sistema.

E está fazendo isso de uma maneira realmente portátil para que um usuário final do pipeline possa trocar a tecnologia que está usando: Docker, Singularity, Apptainer, Conda. Isso realmente não importa para a lógica do pipeline, mas o Nextflow lidará com todas as necessidades de infraestrutura subjacentes, para que execute em qualquer lugar.

E esse é realmente o superpoder do Nextflow. É reprodutibilidade e portabilidade. E com o Nextflow você pode realmente compartilhar seu fluxo de trabalho e outras pessoas podem executá-lo em seus sistemas e vai simplesmente funcionar.

Essa é uma coisa realmente, realmente difícil de fazer, e agora você sabe como fazer isso também com seus fluxos de trabalho.

Ok, é isso para este capítulo. Se você descer até o final do curso, encontrará um quiz novamente sobre alguns contêineres. Espero que tudo tenha feito sentido. É uma maneira realmente legal de trabalhar com análise. E se você é novo em contêineres, espero ter te convencido de que é o caminho a seguir, e você nunca olhará para trás.

Mas com isso, faça uma pequena pausa talvez, e você me encontre em alguns minutos para passar pela parte final seis do Hello Nextflow, que é toda sobre configuração.

Muito obrigado.
