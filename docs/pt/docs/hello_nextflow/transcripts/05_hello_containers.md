# Parte 5: Olá Contêineres - Transcrição

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Tradução assistida por IA - [saiba mais e sugira melhorias](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

<div class="video-wrapper">
  <iframe width="560" height="315" src="https://www.youtube.com/embed/5PyOWjKnNmg?si=QinuAnFwFj-Z8CrO&amp;list=PLPZ8WHdZGxmXiHf8B26oB_fTfoKQdhlik" title="YouTube video player" frameborder="0" allow="accelerometer; autoplay; clipboard-write; encrypted-media; gyroscope; picture-in-picture; web-share" referrerpolicy="strict-origin-when-cross-origin" allowfullscreen></iframe>
</div>

!!!note "Notas importantes"

    Esta página mostra apenas a transcrição. Para instruções passo a passo completas, retorne ao [material do curso](../05_hello_containers.md).

    Os números de seção mostrados na transcrição são fornecidos apenas para fins indicativos e podem não incluir todos os números de seção nos materiais.

## Boas-vindas

Olá, bem-vindo à Parte Cinco do curso de treinamento Hello Nextflow.

Este capítulo se chama Olá Contêineres. Vamos falar sobre como o Nextflow se integra com ferramentas como Docker e Singularity para usar contêineres de software para provisionar software aos usuários do seu pipeline.

Isso significa que quando as pessoas executarem seu pipeline, elas não precisarão instalar todas as diferentes ferramentas por conta própria. O Nextflow fará isso por elas.

Os contêineres são uma tecnologia extremamente poderosa e crucial para reprodutibilidade e facilidade de uso. Vamos começar fazendo uma breve introdução aos próprios contêineres, executando alguns comandos docker manualmente, e então pegaremos esses mesmos contêineres e os colocaremos no nosso pipeline Nextflow.

Ok. Vamos começar.

Então, assim como antes, vamos começar carregando o material de treinamento. Vá para training.nextflow.io. Hello Nextflow, Capítulo Cinco, Olá Contêineres.

Vou entrar no meu ambiente Codespaces e aqui à esquerda vemos hello containers ponto nf.

Assim como antes, este é o mesmo script com o qual terminamos o capítulo quatro anterior, então deve parecer familiar.

Temos nossos parâmetros de linha de comando para especificar o arquivo de entrada e o nome do lote. Estamos incluindo nossos três módulos, e temos nosso fluxo de trabalho onde executamos os três processos.

## 0. Aquecimento: Execute hello-containers.nf

Sinta-se à vontade para executar este fluxo de trabalho novamente e verificar se ele está produzindo as saídas que você espera. Por enquanto, na verdade vou fechá-lo e mergulhar no terminal.

## 1. Use um contêiner 'manualmente'

Para começar este capítulo, vamos fazer uma pequena recapitulação sobre tecnologia de contêineres. Se você está muito acostumado com docker ou singularity ou outras tecnologias de contêiner, então trate isso como uma atualização, ou sinta-se à vontade para pular completamente.

O Nextflow suporta muitos tipos diferentes de tecnologias de contêiner. Isso inclui Docker, Singularity, Podman, Shifter, Charliecloud, e mais.

Neste treinamento, vamos focar no Docker. Ele vem pré-instalado nos code spaces e é uma das tecnologias de contêiner mais populares, especialmente se você está desenvolvendo no seu próprio computador ou laptop.

Se você está trabalhando em um ambiente acadêmico em um HPC compartilhado, você pode descobrir que o Singularity está disponível e não o Docker. Tudo bem. Todos os conceitos são exatamente os mesmos. Alguns dos comandos manuais são diferentes, mas se você entende Docker, você também entenderá singularity.

Na verdade, Singularity também está instalado no ambiente Code Spaces. Então, se quiser, você pode tentar fazer as mesmas tarefas usando Singularity em vez de Docker.

Ok, então o que é tecnologia de contêiner? A ideia por trás do Docker é que ele pode buscar uma imagem de uma fonte remota. Puxá-la para sua máquina local e então criar um contêiner baseado naquela imagem.

Este contêiner em execução é um pouco como uma máquina virtual rodando no seu computador. Ele está isolado do seu ambiente, e vem pré-empacotado com um sistema operacional e um conjunto de software disponível.

## 1.1. Faça o pull da imagem do contêiner

A sintaxe que precisamos para buscar uma imagem pré-existente é "docker pull". Então vou digitar isso no meu terminal, mas agora precisamos de uma imagem para brincar.

Você pode construir imagens você mesmo. Você pode encontrá-las em registros públicos como Docker Hub ou quay.io. Mas uma maneira realmente boa de obter imagens rapidamente é usando Seqera Containers.

Este é um serviço comunitário gratuito que construímos em 2024, que você pode usar sem login ou qualquer coisa.

Se você for para seqera.io/containers ou clicar em containers no topo aqui, você é apresentado a uma interface de busca e pode digitar o nome de qualquer ferramenta disponível no Conda ou no Python Package Index.

Por padrão, ele busca nos canais Bioconda e Conda Forge, mas você pode prefixar qualquer canal Conda. Estou aqui se você quiser.

Para nos divertir um pouco, vamos usar cowpy. Vou digitar cowpy. Me dá resultados do Python Package Index e Conda Forge. Vou clicar nisso para adicioná-lo ao meu contêiner. Eu poderia adicionar múltiplos pacotes aqui se quisesse. Seleciono Docker, seleciono linux/amd64, e clico em Get Container.

Isso constrói a imagem para mim sob demanda se ela ainda não foi criada, e me dá uma URL que posso copiar.

Se você estiver interessado, pode clicar em view Build Details, e isso te leva a uma página que mostra o arquivo de ambiente conda que foi usado e o log de construção completo para a construção, junto com os resultados da varredura de segurança.

Se eu voltar para o meu code spaces, agora posso colar este nome de contêiner e apertar enter.

O Docker agora baixa todas as diferentes camadas dentro desta imagem de contêiner, e agora nos diz que esta imagem está disponível para uso.

## Fazendo pull de uma imagem Singularity

Se você está usando singularity, o processo é basicamente o mesmo. Selecionamos nossos pacotes de imagem, selecionamos cowpy. Agora escolhemos Singularity em vez de Docker e clicamos em Get Container. Isso nos dá uma URL de imagem usando oras://. Ou se você preferir, pode usar https:// marcando aquela caixa. Copie aquela URL. Agora vá para Code Spaces. Na verdade temos Apptainer instalado neste espaço, que é o mesmo que Singularity, mas eles são aliases um do outro. Então vou fazer apptainer pull e então vou chamá-lo de cowpy sif, mas você pode chamá-lo do que quiser. Cole a URL. E isso vai baixar aquela imagem para mim.

Eu poderia fazer ls -lh e ver cowpy.sif

Singularity é diferente do Docker, pois singularity armazena todas as imagens em arquivos planos, enquanto Docker tem um registro onde mantém todas as camadas separadamente na sua máquina host, e tem um daemon em execução para manter controle de tudo isso.

## 1.2. Use o contêiner para executar cowpy como um comando único

Ok, vamos voltar ao Docker. Agora podemos tentar executar esta imagem que criamos fazendo docker run.

Vou fazer dash dash rm, que apenas faz uma execução única da imagem. E vou colar a URL da imagem. E então finalmente, você finaliza isso com um comando que quer executar.

A imagem que geramos tinha cowpy instalado, então vamos tentar cowpy.

Pronto. Executou nosso comando. Não tenho cowpy instalado localmente. Você pode ver que se eu tento executá-lo, ele não existe. No entanto, neste comando, eu o executei usando Docker e ele gerou corretamente esta saída.

## 1.3. Use o contêiner para executar cowpy interativamente

Podemos ir além disso se quisermos e iniciar um contêiner interativamente e olhar por dentro. Novamente, faço "docker run dash dash rm". Agora vou fazer dash it, que diz ao Docker que queremos um terminal interativo. Faço a URL da imagem novamente, e desta vez, em vez de fazer cowpy, vou fazer bin bash porque o comando que queremos executar é bash.

Isso nos leva para dentro deste contêiner em execução e você pode ver que o prompt mudou agora.

Se eu fizer LS slash você pode ver que os diretórios aqui são diferentes.

Se eu abrir um segundo terminal aqui do lado direito, que está apenas executando no GitHub Code Spaces e faço LS slash, você vê que temos diretórios como workspaces e temp, enquanto aqui no Docker é diferente.

Então este ambiente está completamente separado dentro do Docker e isolado do meu ambiente host. Isso é uma coisa boa, porque isso isola a execução deste comando na imagem Docker e mantém reprodutível entre diferentes pessoas em diferentes sistemas host.

Se você quiser usar dados do seu sistema host dentro da imagem Docker, você tem que explicitamente montar isso no contêiner.

Vamos fazer isso em um segundo.

## 1.3.2. Execute o(s) comando(s) da ferramenta desejada

Primeiro, porém, vamos ver se conseguimos executar cowpy. Lá novamente, o comando está disponível agora diretamente na linha de comando, e podemos começar a fazer coisas mais complexas e passar argumentos. Hello containers e em vez da vaca, vamos fazer o pinguim tux. Vamos ver o que mais temos.

Vamos fazer cheese. Maravilhoso. Que tal Dragon e Cow? Muito bom.

## 1.3.3. Saia do contêiner

Ok. Não posso fazer muito mais porque não tenho nenhum dado neste contêiner. Então vamos sair desta imagem em execução e ver se conseguimos montar alguns dados no contêiner. Posso fazer isso fazendo control D ou digitando exit. Ok, agora estou de volta no meu code space regular do GitHub.

## 1.3.4. Monte dados no contêiner

Para montar alguns dados no contêiner Docker, preciso usar dash V. Então vou pegar meu comando docker anterior, voltar ao início fazer dash v. Vou fazer "." para o diretório de trabalho local atual, e então dois pontos para dizer onde isso deve ser montado no diretório host e faço slash data. Então isso está montando este diretório particular no contêiner em slash data.

Agora se eu fizer LS slash podemos ver que temos um novo diretório chamado data, e se eu fizer LS data, você pode ver todos os arquivos que temos na barra lateral aqui. Fantástico.

## 1.3.5. Use os dados montados

Agora podemos começar a usar alguns dos arquivos que estão no sistema host dentro da imagem Docker. Então posso dizer cat data greetings csv. Se você se lembra, este é nosso arquivo CSV com nossas diferentes saudações de antes, e posso direcionar isso para cowpy. Fantástico. Agora estamos chegando a algum lugar.

Ok. Isso é suficiente para executar Docker interativamente. Esperançosamente você agora tem uma noção do que é Docker e como usá-lo tanto para executar um comando de forma única, quanto para usar uma imagem interativamente. Se você está usando singularity. Os comandos são todos muito similares exceto que você faz coisas como apptainer exec ou apptainer run, ou singularity exec ou singularity run.

## 2. Use contêineres no Nextflow

A seguir vamos voltar ao nosso fluxo de trabalho Nextflow e ver como usar esta tecnologia dentro do pipeline Nextflow.

Vamos fechar o terminal e abrir Hello Containers novamente.

## 2.1. Escreva um módulo cowpy

Para ficar com nosso exemplo cowpy, vamos criar um novo processo no nosso fluxo de trabalho, que usa cowpy. Vamos até módulos, criar um novo arquivo e chamá-lo de cowpy nf. Agora vou trapacear um pouco e copiar o código do material de treinamento e apertar salvar. E vamos dar uma olhada.

Então este é um processo simples. Esperançosamente agora você entende como são os blocos de construção de um processo. Temos nosso publishDir novamente, indo para results. Temos duas entradas, um arquivo de entrada e uma string chamada character. Temos uma saída cowpy input file, e temos um script que parece exatamente o mesmo que o que executamos manualmente dentro da nossa imagem docker há um segundo atrás: cat para imprimir um arquivo, direcionando para cowpy, dizendo qual tipo de caractere cowpy queremos usar, e direcionando isso para o arquivo de saída, que passamos como a saída aqui.

## 2.2. Adicione cowpy ao fluxo de trabalho

Ok, vamos voltar ao nosso fluxo de trabalho, importar este novo processo. Então cowpy de modules cowpy nf. Vamos criar um novo parâmetro para que possamos especificar qual caractere queríamos. Vamos dizer Turkey por padrão. E então vamos chamar este novo processo no final do fluxo de trabalho,

cowpy. E vamos usar a saída aqui de Collect Greetings. Então collect greetings out, out file aqui. E então precisamos de um segundo argumento, que são os novos params que acabamos de fazer. params ponto character.

## 2.2.4. Execute o fluxo de trabalho para verificar se funciona

Ok, vamos ver se nosso novo processo funciona. Nextflow run hello containers. Isso deve executar aqueles primeiros três processos e então tentar executar cowpy no final.

Tivemos um erro. O que está dizendo aqui, cowpy teve um erro e teve um status de saída 127 e com certeza, comando sh cowpy comando não encontrado.

Não dissemos ao Nextflow que temos uma imagem Docker disponível para cowpy, então ele tentou executá-lo no nosso sistema host e não temos cowpy instalado no nosso sistema host, então disparou um erro.

## 2.3. Use um contêiner para executá-lo

Então o que precisamos fazer é precisamos dizer ao Nextflow que temos um contêiner disponível. Vamos ao nosso processo cowpy e vamos adicionar uma nova diretiva no topo do processo chamada container.

Então encontramos nossa imagem, copiamos a URL, e colocamos isso em uma string.

Isso não é suficiente por si só porque um pipeline Nextflow pode ter várias maneiras de especificar software. Eu poderia também fazer conda conda-forge cowpy, por exemplo. E o Nextflow precisa saber qual dessas tecnologias você quer usar.

## 2.3.2. Habilite o uso de Docker via o arquivo nextflow.config

Então para executar com Docker habilitado, vamos nos adiantar um pouco e usar o arquivo de configuração Nextflow, que é algo que vamos cobrir com mais detalhes no próximo capítulo. Você pode ver neste diretório que temos um arquivo chamado Nextflow Config, e aqui você já tem docker.enabled False.

Vamos mudar isso para True para habilitar Docker, e então podemos tentar executar o fluxo de trabalho novamente.

## 2.3.3. Execute o fluxo de trabalho com Docker habilitado

Nextflow run hello containers nf e desta vez cowpy executou com sucesso. Vamos olhar em Results. cowpy collected test e lá está nosso Turkey. Maravilhoso.

Então nos bastidores ali, o Nextflow sabia que tinha um contêiner disponível para aquele processo.

Ele buscou a imagem e executou os comandos para nós.

## 2.3.4. Inspecione como o Nextflow lançou a tarefa conteinerizada

Se você está curioso, podemos realmente ver exatamente o que ele fez olhando no diretório work. Se eu fizer code work, e então o hash e então command run, que se você se lembra é o arquivo real que é executado para aquela tarefa, podemos entrar e podemos procurar por uma função chamada NXF launch. E aqui você pode ver o comando docker exato que o Nextflow usou, que parece muito com o que estávamos fazendo manualmente no terminal anteriormente. Docker run. Vinculando este diretório host no contêiner, e então especificando a URL do contêiner.

Então não há mágica aqui. É apenas que o Nextflow está automaticamente fazendo o trabalho pesado para você de uma forma que significa que você pode facilmente especificar contêineres no seu pipeline, que estão então prontamente disponíveis para qualquer outra pessoa usar que execute seu fluxo de trabalho. E essas pessoas não precisam mais pensar sobre gerenciar software para executar seu pipeline de análise.

Muito, muito simples, muito conveniente, e também realmente reprodutível. Bom em todos os aspectos.

Ok, ótimo trabalho. Esse é o fim do Capítulo Cinco. Junte-se a nós no próximo vídeo para a parte seis, que é a parte final deste treinamento Hello Nextflow, onde falaremos sobre configuração Nextflow com mais detalhes.

Vejo você no próximo vídeo.

[Transcrição do próximo vídeo :octicons-arrow-right-24:](06_hello_config.md)
