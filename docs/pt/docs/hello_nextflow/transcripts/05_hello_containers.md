# Parte 5: Hello Containers - Transcrição do Vídeo

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Tradução assistida por IA - [saiba mais e sugira melhorias](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

<div class="video-wrapper">
  <iframe width="560" height="315" src="https://www.youtube.com/embed/Xqr--bKEN9U?si=QinuAnFwFj-Z8CrO&amp;list=PLPZ8WHdZGxmWKozQuzr27jyMGqp9kElVK&amp;cc_load_policy=1" title="YouTube video player" frameborder="0" allow="accelerometer; autoplay; clipboard-write; encrypted-media; gyroscope; picture-in-picture; web-share" referrerpolicy="strict-origin-when-cross-origin" allowfullscreen></iframe>
</div>

!!!note "Notas importantes"

    Esta página mostra apenas a transcrição. Para instruções completas passo a passo, retorne ao [material do curso](../05_hello_containers.md).

    Os números de seção mostrados na transcrição são fornecidos apenas para fins indicativos e podem não incluir todos os números de seção dos materiais.

## Boas-vindas e contexto

Olá, e bem-vindo de volta ao Hello Nextflow. Esta é a parte cinco chamada Hello Containers. E nesta parte do curso, vamos falar tudo sobre como encapsular os requisitos de software para um pipeline, para que as pessoas executando o pipeline não precisem pensar em instalar o software.

Se você trabalha em bioinformática há tanto tempo quanto eu, você pode se lembrar do que eu costumo chamar de "os velhos tempos ruins", onde quando você queria executar o pipeline de outra pessoa ou replicar seu trabalho, você passava horas ou dias tentando instalar todas as diferentes ferramentas de software que eles usaram, nas mesmas versões, tentando compilá-las na sua máquina, e era um pesadelo. Era realmente difícil.

Se você está executando em um HPC, você pode ter usado módulos de ambiente onde os administradores de sistema tentavam instalar software para você, o que era, ok, mas ainda imperfeito.

Mas agora temos maneiras melhores de fazer isso. O Nextflow tem suporte integrado para diferentes tecnologias de contêineres de software. Docker é a mais comum. É a que vamos usar hoje. Funciona bem no Codespaces. Funciona bem no seu computador local e funciona bem na nuvem.

Mas também Singularity ou Apptainer, que são muito comuns em sistemas HPC e efetivamente funcionam exatamente da mesma maneira. Ou Podman, Shifter, há um monte de outros que são todos muito similares.

O único extra, que é meio similar mas não exatamente, que o Nextflow suporta é o Conda. E o Nextflow pode gerenciar ambientes Conda para você por processo, o que é muito melhor do que fazer seus próprios ambientes Conda. E novamente, pode ser distribuído com um pipeline.

Vamos começar este capítulo falando um pouco sobre tecnologias de contêineres e Docker e como eles funcionam. E vamos fazer a primeira metade apenas manualmente no Docker para que você entenda o que está acontecendo por baixo dos panos e como isso funciona. Porque isso é realmente importante para entender o que o Nextflow está fazendo e como entender o que seu fluxo de trabalho está fazendo quando está sendo executado.

Então. Vamos pular para nosso Codespaces. Agora eu limpei tudo novamente, mas se pularmos para Hello Containers, você deve ver que todos os nossos scripts e tudo estão lá, iguais ao final do capítulo de módulos. Então temos nossos diferentes módulos aqui, que eu criei no diretório modules.

Eles ainda estão lá. Eles precisam estar lá para que possa executar. E o fluxo de trabalho e a saída são todos iguais, exceto que mudamos o caminho de publicação de saída para Hello Containers, para que seus arquivos acabem naquele diretório.

Podemos executar isso agora para verificar se funciona, se você quiser, ou podemos continuar com o terminal.

## 1. Usar um contêiner 'manualmente'

Vamos usar o Docker para gerenciar nossos contêineres, e posso verificar se está instalado no meu Codespaces fazendo "docker -v", que me mostra a versão que está instalada e tudo, e que está funcionando corretamente.

Agora contêineres e Docker têm dois conceitos que são realmente importantes. Um é chamado de imagem, e um é chamado de contêiner. A imagem é o snapshot, se você quiser, de todo o sistema de arquivos que você estará usando, e o contêiner é o ambiente em execução. Então você cria um contêiner usando uma imagem.

Uma vez que você está naquele contêiner, ele normalmente funciona como um sistema operacional completo. Ele é isolado do mundo exterior. Ele é separado de tudo mais, e isso é uma coisa boa. É assim que obtemos tão boa reprodutibilidade com o Nextflow.

Porque para tarefas executadas dentro de um contêiner, elas não são contaminadas por nenhum arquivo de configuração no seu sistema local. Quaisquer outras influências externas, elas executam em sua própria pequena caixa de areia. Os arquivos são então produzidos de uma maneira muito, muito reprodutível porque você está usando as mesmas bibliotecas subjacentes, todas as mesmas dependências, exatamente o mesmo software para cada pessoa executando em cada ambiente de computação diferente. O que francamente eu acho fantástico e incrível que funcione. E mesmo, mesmo até hoje ainda meio que me surpreende que isso seja possível.

## 1.1. Baixar a imagem do contêiner

Então vamos experimentar usar algumas imagens Docker e Docker, quando você o executa no seu sistema, tem um registro docker no seu computador, ou neste caso, no code space, que mantém o controle de todas as diferentes imagens que foram baixadas e usadas no passado, e as diferentes camadas das quais elas são construídas.

Podemos ver quais imagens temos localmente com Docker fazendo "docker image ls". E neste caso você pode ver que há um monte de imagens Docker aqui, que todas têm a ver com configurar este Codespaces. Todas têm a ver com dev containers e coisas assim. Então você não precisa se preocupar muito com elas, mas à medida que adicionamos mais imagens e as baixamos, conforme este curso avança, você pode verificar essa lista e verá que o registro local está mantendo o controle de todas essas coisas que baixamos.

Mas vamos pegar uma nova fazendo "docker pull". E isso diz ao Docker para buscar uma nova imagem da web.

Então colocamos o URI para aquele contêiner. Agora isso poderia ser uma imagem docker que você construiu localmente e então enviou para a internet. Poderia ser uma imagem que outra pessoa fez. Há muitas, muitas, muitas maneiras diferentes de fazer imagens Docker, mas sem dúvida uma das maneiras mais simples é terceirizar isso, e conseguir que outra pessoa faça isso para você.

E o que vamos usar neste tutorial é um serviço da Seqera chamado Seqera Containers.

Agora, Seqera Containers é totalmente gratuito, e usa um software de código aberto que desenvolvemos chamado Wave, que foi construído para gerenciar contêineres de uma maneira complementar ao Nextflow. E ele lida com muitos dos casos de uso comuns que nos encontramos lidando, com Nextflow.

É muito comum que o software que precisamos esteja empacotado no Conda, nos canais Bioconda ou conda-forge ou outros canais mais específicos de domínio. E Wave e Seqera Containers são realmente bons em construir imagens a partir disso.

Então posso ir a esta interface web e vamos mexer com o pacote chamado "cowpy". Então eu digito o nome do pacote que quero. Ele pesquisa, encontrou no índice de pacotes Python, então posso usar isso. Ou se eu esperar um pouco mais, está pesquisando bioconda e conda-forge. E você pode ver, posso especificar qualquer canal conda aqui. Então se você quiser encontrar um canal Nvidia ou qualquer outra coisa, isso também deve funcionar.

E então posso especificar se quero que ele construa uma imagem docker para mim ou uma imagem singularity e também qual arquitetura de CPU eu quero usar. Então amd64 ou arm64.

E uma vez que os resultados do bioconda são listados, agora posso ver todas as diferentes versões que estão disponíveis também. Vou colocar isso. E agora eu poderia continuar pesquisando e obter mais pacotes do Conda se eu quiser e compor este contêiner como eu quiser, mas eu só quero aquele. Então vou clicar em Get Container.

Agora, outra pessoa já solicitou o mesmo contêiner antes e ele é retornado de um registro, então nós apenas o obtemos imediatamente. Mas se ninguém mais tivesse pedido este pacote de software ou esta combinação de pacotes de software, Wave e Seqera Containers o construiriam em tempo real para nós.

Podemos copiar esta URL e também podemos ver os detalhes da construção. E isso nos mostra o que o serviço fez no backend. Ele criou um arquivo de ambiente conda. Um arquivo docker, e então é isso, executando o processo de construção docker. Ele também executou uma varredura, uma varredura de segurança, então você pode ver quaisquer CVEs. E ele te diz quando isso foi criado.

Wave e Seqera Containers podem fazer muito mais do que isso, mas este é meio que um caso de uso simples, que é o mais comum. E devo dizer que essas imagens são hospedadas por pelo menos cinco anos. Então você pode construir essas URLs em seus pipelines e saber que elas não vão desaparecer tão cedo.

Então eu tenho minha URL para minha imagem docker para cowpy.

Agora posso fazer "docker pull" daquela URL, e ela vai buscar todas as diferentes camadas e baixar esta imagem para que esteja disponível para mim localmente.

## 1.2. Usar o contêiner para executar cowpy como um comando único

Ok, agora vamos tentar realmente usá-lo. Então agora vou usar um comando "docker run" em vez de "docker pull", e vou usar a flag "--rm", que apenas diz ao Docker para desligar este contêiner uma vez que terminou de fazer o que pedi. E então coloco o identificador para o contêiner, que é apenas um URI.

E então no final, especifico o comando que quero que o Docker execute dentro do contêiner gerado a partir desta imagem. Vou apenas dizer cowpy, que é o nome da ferramenta que está instalada do Conda Forge, que está disponível dentro da imagem.

Vou apertar enter e pronto. Executamos cowpy em um sistema. Temos uma pequena vaca nos dando algumas informações.

Agora note que cowpy não está instalado no meu sistema local. Então se eu executá-lo apenas sem todo o material do Docker, ele diz, comando não encontrado. Então isso baixou uma imagem. Ele criou um contêiner usando Docker, e então entrou naquele contêiner e executou este comando para nós e nos deu a saída de volta ao nosso terminal. Muito, muito legal.

## 1.3. Usar o contêiner para executar cowpy interativamente

Ok, vamos ir um passo além agora e executar este contêiner interativamente e dar uma olhada, para que possamos ver o que está acontecendo dentro do contêiner.

Então se eu voltar e pegar meu comando run e vou me livrar do cowpy no final, porque eu realmente não quero executar cowpy. Eu quero executar um terminal Bash.

E então vou voltar aqui e vou fazer "-it", que significa Interactive e Terminal ou TTY, e vou apertar enter.

E agora você pode ver o prompt, a parte antes de eu digitar, mudou. Este era o prompt do Codespaces onde dizia o diretório, e agora diz base e roots e tmp. Então agora estou dentro do contêiner, e se eu fizer "ls", você verá que os arquivos que vejo neste diretório são diferentes dos arquivos que tenho no meu workspace.

E de fato, não posso ver nenhum dos arquivos do meu workspace local do codespaces ou do meu drive local dentro do contêiner Docker. O runtime do contêiner docker está completamente isolado e não pode escrever ou ler nenhum arquivo de um sistema de arquivos host externo.

Eu posso, no entanto, ver o software que está instalado dentro do contêiner e executá-lo. Então posso executar cowpy e podemos ver um pouco mais sobre como usar cowpy. Aqui posso fazer "cowpy 'Hello World'" e isso diz a ele para realmente colocar minha citação dentro de um pequeno balão de fala. E você também pode executar diferentes tipos de vacas, então não precisa ser uma vaca. Você pode fazer um "-c". E estou na Suécia, então vou escolher um alce. Muito legal. Dei a ele alguns chifres.

E há um monte de diferentes que você pode brincar, que você pode ver descritos nos documentos de treinamento.

## 1.3.4. Montar dados no contêiner

Ok. Seria bom se pudéssemos executar cowpy nos arquivos em nosso sistema de arquivos.

Claro, não é muito útil apenas ter o contêiner e nenhum acesso a nada. Pode ser seguro e reprodutível, mas não é muito útil.

Então como fazemos isso? Vou sair deste contêiner Docker digitando exit, e você pode ver que o prompt nos diz que agora estamos de volta ao nosso Codespaces regular novamente.

E vou executar o mesmo comando novamente. Mas desta vez vou adicionar algumas flags adicionais aqui. E a importante é "-v", que significa montar um volume, que é basicamente uma parte, parte de um espaço em disco.

O "-v" recebe duas partes: há como uma string e então dois pontos e uma string. E a primeira parte é o sistema de arquivos local, que deve ser montado no contêiner. E então a segunda parte é onde isso deve acabar dentro do contêiner.

Agora eu só quero carregar todo o meu sistema de arquivos local aqui. Então "." é o diretório de trabalho atual. Então vou apenas fazer "." e então ":", e então vamos colocar isso em um novo diretório dentro do contêiner chamado "my_project". Isso poderia realmente ser chamado de qualquer coisa.

E então vou executar novamente.

No diretório de trabalho onde sou colocado, que é /tmp, os arquivos não estão lá. Mas se eu fizer "ls my_project", lá está: todos os mesmos arquivos que tínhamos localmente em nosso Codespaces agora estão disponíveis dentro do contêiner naquele caminho.

Este é acesso de leitura e escrita, então posso criar novos arquivos neste diretório e eles aparecerão no meu sistema de arquivos host. Este diretório em particular, então se comporta exatamente como se eu estivesse fora do contêiner, então agora posso ler e escrever e fazer coisas.

## 1.3.5. Usar os dados montados

Ok, vamos apenas provar que podemos fazer isso. Faço "cat /my_project/data/greetings.csv". Se você se lembra, o conteúdo deste arquivo se parece com isso. Agora posso canalizar isso para cowpy e a vaca vai imprimir as diferentes saídas daquele arquivo em seu pequeno balão de fala, o que é meio divertido.

Então você pode ver, agora podemos usar o software no contêiner para interagir com os arquivos em nosso sistema host.

Ok, vamos voltar e vamos continuar com o resto do material de treinamento.

## 2. Usar contêineres no Nextflow

Então isso é muito legal usar contêineres. Espero que faça sentido. E você pode ver o valor desses contêineres e por que isso é útil para executar software de análise.

Mas como fazemos todo esse mesmo processo dentro do Nextflow? Não queremos estar executando um monte de comandos Docker nós mesmos. Queremos apenas deixar o Nextflow lidar com tudo isso para nós.

Então vamos trabalhar nisso. Vamos adicionar um novo processo ao nosso pipeline, para executar cowpy. Ok, então vamos criar um novo módulo para nosso novo processo. Então vá em modules, vamos chamá-lo de cowPy.nf, e então vou copiar o código do material de treinamento aqui.

Mas você pode ver que o processo é muito simples. Parece muito com os que fizemos até agora, temos um bloco de entrada com um path, que é nosso arquivo de entrada, e também um value aqui para que este seja um caractere, então poderíamos usar um alce novamente se quisermos.

E então uma saída, que é um único arquivo aqui, um path e então um script. E estamos fazendo a mesma coisa que fizemos interativamente dentro do contêiner: estamos fazendo "cat" para ler o arquivo de entrada. Estamos canalizando aquele conteúdo para cowpy. Estamos escolhendo um caractere específico baseado naquela entrada, estamos escrevendo para um arquivo de saída chamado cowpy input file, que é então ecoado para a saída.

Ótimo. Vamos incluir isso. Então include \{ cowpy \} from "./modules/cowpy.nf", eu chamei de cowpy? Sim.

E então vamos chamar nosso novo processo aqui embaixo no bloco main do fluxo de trabalho. Então vamos executar cowpy. E vamos pegar nosso novo processo cowpy e vamos dizer collectGreetings.out.

E então se você se lembra, havia duas saídas para este módulo. Uma chamada outfile e uma chamada report. A extensão VS Code está auto-sugerindo essas para nós e queremos o .outfile.

Você sempre pode pular para este processo aqui. Você pode passar o mouse sobre ele e deve mostrar rapidamente quais eram as saídas. E também podemos fazer command click nele e ele abrirá o arquivo do módulo se você quiser ver com mais detalhes.

Então aqui vamos nós. Esse é o outfile lá, e esse é o path. Então isso agora será o arquivo de entrada para nosso processo cowpy. Fantástico.

Agora se você se lembra, um processo cowpy tem duas entradas. Também tínhamos o canal de valor para o caractere. Então podemos adicionar "params.character" aqui. Eu poderia ter codificado isso se quisesse, mas vamos torná-lo uma opção CLI para que possamos fazer dash, dash character.

Certo. Agora preciso definir o parâmetro de entrada que acabamos de chamar e dar a ele um padrão. Então character, String. E eu gosto do alce, então vou configurá-lo para moose por padrão.

Certo, vamos tentar executá-lo. Então se eu fizer Nextflow run hello containers, vamos ver o que acontece.

Eu poderia ter usado dash resume se tivesse os diretórios de trabalho antigos por aí. E novamente, esses primeiros processos teriam sido armazenados em cache e teria sido um pouco mais rápido, mas deve ser basicamente o mesmo.

Agora podemos ver imediatamente que lançou um erro quando chegou ao nosso novo processo, está nos dizendo aqui que houve um erro executando o processo cowpy e ele saiu com um status de saída 127. Este é o comando que tentou executar. Parece certo, parece como esperávamos. Está pegando aquele nome de arquivo de saída, que parece certo, está executando com um caractere moose e tentando salvar.

Mas você pode ver o erro de comando aqui está dizendo que o comando cowpy não foi encontrado. E isso faz sentido porque ainda não dissemos ao Nextflow para usar um contêiner. Apenas demos a ele o comando cowpy. E como eu disse antes, cowpy não está instalado em nosso sistema local. Então quando tentou executá-lo, falhou.

## 2.3.1. Especificar um contêiner para cowpy

Precisamos dizer ao Nextflow que há um contêiner disponível e ele pode usá-lo. Então como fazemos isso?

Se entrarmos em nosso módulo aqui, vamos adicionar uma nova declaração no topo chamada "container". E vamos então configurá-la para uma string.

Agora, se você se lembra, lá no Seqera Containers, posso copiar aquela URL e apenas solto isso entre aspas aqui.

Agora volte e tente executá-lo novamente.

Deixe-me ver se funciona desta vez.

Infelizmente, falha exatamente da mesma maneira, mesmo que agora tenhamos definido um contêiner para o processo executar. Então, para usar nossa imagem docker, precisamos dizer ao Nextflow para habilitar o uso do Docker quando executamos o fluxo de trabalho.

E vamos fazer isso criando um novo arquivo de configuração. Então vou dizer touch nextflow.config.

Este é um nome de arquivo especial onde se estiver no diretório de trabalho enquanto eu lanço o pipeline, ele será carregado automaticamente. Então se eu entrar neste arquivo Nextflow dot config, você pode ver que ele realmente já existe, o que eu tinha esquecido. E temos docker.enabled aqui já, mas está configurado para false, que é o padrão.

Então se eu mudar isso para equals True em vez disso, docker.enabled. E há documentos de referência para todos esses escopos de configuração nos documentos do Nextflow. E também você pode ver que quando passo o mouse com uma extensão VS Code, ele puxa os documentos específicos para isso e me diz o que significa e como configurá-lo.

Então agora configuramos para true, e se eu executar o Nextflow novamente, o Nextflow agora saberá buscar aquela imagem docker para nós se ainda não a tivermos localmente, e então executar aquele processo com aquele ambiente de contêiner.

E então podemos ver que ele executou com sucesso e temos um pequeno tick ao lado de cowpy. Fantástico. Se eu subir e olhar no diretório de resultados, o arquivo ainda não está lá. E isso é porque ainda precisamos publicar este arquivo de saída assim como todos os outros.

Então vamos ao bloco published dentro do fluxo de trabalho, dizer mycowpy equals cowpy.out.

E então aqui embaixo no bloco de saída, mycowpy, chaves path. Ops. Hello containers. Mode, copy.

Se eu executar novamente agora, deve executar exatamente da mesma maneira. Eu poderia ter usado dash resume e esqueço toda vez. E então eu subo e agora temos um novo arquivo criado chamado cowpy-COLLECTED, e lá está meu alce dizendo BONJOUR, HELLO, HOLÀ Fantástico.

Agora, claro, eu também poderia passar agora "--character". Quais são as diferentes opções? Acho que há um Peru? Então posso usar character Turkey. Vai executar exatamente da mesma maneira. Perdi outra oportunidade de usar dash resume, e agora se carregarmos nosso arquivo e agora temos um Peru. Fantástico.

## 2.3.4. Inspecionar como o Nextflow lançou a tarefa containerizada

Ok. Última coisinha. Vamos apenas executar rapidamente este comando novamente, resume desta vez, e dar uma olhada rápida no diretório de trabalho para ver o que é que o Nextflow está fazendo por baixo dos panos para fazer tudo isso funcionar para nós.

Desta vez é super rápido, vamos entrar neste diretório de trabalho, cd work/. Agora se você se lembra temos um monte de arquivos ponto aqui e o que nos interessa neste caso é o que eu disse que quase nunca precisamos olhar, chamado .command.run.

Se eu fizer code dot command run, vai abri-lo no editor. E posso pesquisar neste arquivo e se eu rolar para baixo devo ver Docker run. E você pode ver que o Nextflow está fazendo o comando docker run para nós, quando o Docker está habilitado em uma configuração. Tem um monte de diferentes flags e coisas aqui, mas você pode ver a flag "-v" que usamos nós mesmos quando estávamos executando. E você pode ver que está montando o diretório de workspace local no contêiner, para que o contêiner possa acessar nossos arquivos de entrada e salvar as saídas. E então no final, também está executando .command.sh, que é o script gerado, que tem o comando cowpy dentro.

E então você pode ver que o Nextflow está pegando a lógica do fluxo de trabalho, que é o material com o qual realmente nos importamos, que é específico para nossa análise, e está fazendo todo o material inteligente nos bastidores para fazer o Docker funcionar em nosso sistema.

E está fazendo isso de uma maneira realmente portável para que um usuário final do pipeline possa trocar a tecnologia que está usando: Docker, Singularity, Apptainer, Conda. Isso realmente não importa para a lógica do pipeline, mas o Nextflow vai lidar com todas as necessidades de infraestrutura subjacentes, para que execute em qualquer lugar.

E esse é realmente o superpoder do Nextflow. É reprodutibilidade e portabilidade. E com o Nextflow você pode realmente compartilhar seu fluxo de trabalho e outras pessoas podem executá-lo em seus sistemas e vai simplesmente funcionar.

Isso é uma coisa realmente, realmente difícil de fazer, e agora você sabe como fazer isso também com seus fluxos de trabalho.

Ok, é isso para este capítulo. Se você descer até o final do curso, você encontrará um quiz novamente sobre alguns contêineres. Espero que tudo tenha feito sentido. É uma maneira realmente legal de trabalhar com análise. E se você é novo em contêineres, espero ter te convencido de que é o caminho a seguir, e você nunca vai olhar para trás.

Mas com isso, faça uma pequena pausa talvez, e você se junta a mim em alguns minutos para passar pela parte final seis do Hello Nextflow, que é tudo sobre configuração.

Muito obrigado.
