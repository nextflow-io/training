# Orientação - Transcrição do Vídeo

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Tradução assistida por IA - [saiba mais e sugira melhorias](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

<div class="video-wrapper">
  <iframe width="560" height="315" src="https://www.youtube.com/embed/G3CV-FcV-rc?si=nyLvwhrSB2m1NPc5&amp;list=PLPZ8WHdZGxmXiHf8B26oB_fTfoKQdhlik" title="YouTube video player" frameborder="0" allow="accelerometer; autoplay; clipboard-write; encrypted-media; gyroscope; picture-in-picture; web-share" referrerpolicy="strict-origin-when-cross-origin" allowfullscreen></iframe>
</div>

!!!note "Nota importante"

    Esta página mostra apenas a transcrição. Para instruções passo a passo completas, retorne ao [material do curso](../00_orientation.md).

## Boas-vindas

Olá, bem-vindo ao Hello Nextflow. Meu nome é Phil Ewels. Sou Gerente de Produto para Open Source na Seqera, e estou muito feliz em estar aqui hoje para guiá-lo por este primeiro curso de treinamento em Nextflow.

Vamos passar pelos fundamentos do Nextflow, explicando como escrever e executar pipelines e configurá-los.

E você vai construir seu próprio pipeline simples de múltiplas etapas. Vamos cobrir terminologia como operadores e fábricas de canais, e ao final do curso, você estará pronto para começar a construir seus próprios pipelines de bioinformática.

Se você tiver alguma dúvida, entre em contato em community.seqera.io. Temos uma comunidade Nextflow muito ativa, há uma seção dedicada ao treinamento, então apenas nos avise onde você está com dificuldades e alguém poderá ajudar.

Certo. Vamos começar.

## Site de Treinamento

Todo o material de treinamento para os cursos de Nextflow está em training.nextflow.io. Você pode acessá-lo no seu navegador. Então abra isso agora e podemos dar uma olhada.

Estarei executando isso com a versão 2.1.1. Fazemos pequenas atualizações e correções aqui e ali, então não se preocupe se estiver um pouco diferente, mas se o material tiver mudado muito, você sempre pode usar este seletor de versão no topo para escolher a versão exata dos materiais que vou estar abordando.

Se você prefere o modo claro, pode alterar o tema do site aqui.

Veja as traduções aqui, embora no momento da gravação, realmente só o inglês cobre este novo material.

E também veja todo o código-fonte do site de treinamento e tudo com o que trabalharemos no GitHub.

A página inicial aqui lista todos os diferentes cursos de material de treinamento que temos. Então eu rolo para baixo, veremos Nextflow para iniciantes com o curso Hello Nextflow que faremos aqui. Você pode ver todos os outros cursos que também temos, que funcionam de maneira similar.

## Configuração do Ambiente

Na verdade, vou começar usando este primeiro aqui no topo, que é comum para todos os cursos de treinamento, e é especificamente sobre configurar nosso ambiente.

Eu clico, me leva a esta seção, e podemos ver instruções para desenvolver localmente. Se você quiser usar seu próprio laptop com sua própria cópia do VS Code e suas próprias instalações de software, ou o que esperamos que a maioria das pessoas faça, que é usar algo chamado GitHub Codespaces.

Codespaces é um serviço fornecido pelo GitHub onde eles executam um servidor web na nuvem, ao qual você pode se conectar. Esse servidor tem o VS Code instalado, onde você pode executá-lo no seu navegador, ou se preferir, conectá-lo à sua instalação local do VS Code. Toda a computação, todos os arquivos, toda a edição acontecem remotamente, o que significa que todo o software que você precisa vem pré-instalado e é o mesmo para todos.

## Criando um GitHub Codespace

Para criar o codespace com tudo que precisamos, procure pelos botões no material da documentação, que dizem "Abrir no GitHub Codespaces". Vou clicar nisso agora, abrir em uma nova aba. E me é apresentada esta página web. Agora você pode ver que está pré-configurado para definir com nextflow-io training.

Posso apenas clicar em criar novo codespace. Mas na verdade recomendamos que usemos uma máquina um pouco maior para o treinamento Nextflow com quatro CPUs em vez de duas. Você pode alterar qual versão do material ele usa. Então isso está padronizado para 2.1.1 porque essa é a versão dos documentos de onde segui o link. Mas eu também poderia defini-lo para um branch específico do repositório se eu quiser.

Agora vou clicar em criar codespace. E ele vai começar a configurar o ambiente para mim.

## Criação do Codespace

Agora, a primeira vez que você fizer isso, vai levar bastante tempo, então agora é um bom momento para ir tomar uma xícara de chá. Fique confortável, converse com a pessoa ao seu lado.

Se você estiver interessado, pode clicar em building codespace aqui embaixo para ver os logs da configuração. E você pode ver aqui que está baixando uma imagem Docker com tudo que preciso e configurando o ambiente.

Agora, você só tem que esperar assim na primeira vez que criar um codespace. Se você for para github.com/codespaces aqui, verá todos os diferentes Codespaces que você tem abertos. Aqui está o que acabei de criar. Na próxima vez que fizer isso, você pode vir aqui e pode selecionar o codespace anterior e simplesmente voltar direto para ele. E é um processo muito, muito mais rápido para reativar aquele ambiente existente. Isso também manterá todas as alterações que você fez no VS Code e nos arquivos, então você não perderá seu progresso se sair e voltar.

Você pode clicar nos três pontos aqui para realizar outras ações. Por exemplo, se você o configurou com duas CPUs e agora quer quatro, pode alterar o tipo de máquina. Ou se você quiser começar do zero e fresco, pode deletar o codespace.

## Introdução ao VS Code

Ok, Codespaces terminou de configurar meu ambiente e agora me é apresentado o VS Code no navegador.

Se você está acostumado com o VS Code. Isso vai parecer muito familiar se você não o usou antes, é bastante simples. Há algumas partes diferentes da página que você precisa conhecer.

Aqui à esquerda, temos a barra lateral. Você pode ver o Explorer configurado com todos os diferentes arquivos no repositório GitHub do repositório de treinamento.

Nestes botões na parte esquerda, podem estar diferentes ferramentas. Na barra lateral. Posso pesquisar todos os arquivos em todo o projeto. Posso trabalhar com Git, posso trabalhar com GitHub, todos os diferentes tipos de coisas assim.

No topo aqui está o menu principal. O explorador de arquivos é o que teremos aberto na maioria das vezes aqui, e você pode clicar com o botão direito em qualquer um desses arquivos e fazer as coisas normais que você esperaria. Você pode precisar clicar através de alguns avisos como este onde ele como cortar copiar e você pode baixar para sua máquina local também.

Quando o codespace carrega, ele nos dá uma prévia do arquivo markdown nesta área principal aqui. Este é exatamente o mesmo que o que renderiza em github.com. Posso fechar isso e se eu clicar duas vezes naquele arquivo Readme, você verá que ele abre como código no editor de código e assim como com qualquer outro arquivo, podemos editar este código diretamente.

Finalmente, na parte inferior aqui, temos a janela do terminal. Eu estava olhando os logs enquanto ele construía, então isso é o que a coisa atual está mostrando. Também posso pressionar este botão de mais para iniciar uma nova sessão de terminal. Isso não está rodando na minha máquina. Lembre-se, isso está rodando na nuvem, e se eu fizer tree três até profundidade de dois, você verá todos os mesmos arquivos aqui, que estavam à esquerda.

## Mostrando apenas os arquivos "hello-nextflow"

Este repositório GitHub contém todos os diferentes conjuntos de treinamento, não apenas o que estamos fazendo. Então, se você quiser, pode focar apenas na pasta Hello Nextflow. Uma maneira de limpar isso um pouco é ir ao menu arquivo e então adicionar pasta ao workspace.

Clicamos nisso vamos para training. Hello nextflow, e clicamos em adicionar. Vai atualizar sua tela. E então no Explorer, agora temos dois workspaces diferentes, o que tínhamos antes para training e um com apenas Hello Nextflow.

Se você quiser, pode clicar com o botão direito em training e clicar em remover pasta do workspace para se livrar completamente dela da barra lateral.

Agora temos apenas os arquivos para este curso de treinamento específico na lateral. Posso esconder aquele aviso e agora posso fazer a mesma coisa no terminal aqui e fazer CD para mudar de diretório. Hello, Nextflow. E novamente, temos os mesmos arquivos aqui, que estão na barra lateral.

## Hello Nextflow: arquivos

Olhando esses arquivos para o curso Hello Nextflow.

Temos um monte de arquivos .nf, que são para Nextflow, e há um desses arquivos para cada um dos capítulos do curso de treinamento. Trabalharemos através desses arquivos e os modificaremos nos exercícios.

Também temos um arquivo nextflow.config, que apenas tem configurações básicas de config para executar Nextflow neste ambiente, com as quais você realmente não precisa se preocupar neste ponto. Um arquivo greetings.csv, que usaremos para processar dados, que será introduzido na próxima parte deste curso, e um arquivo test-params.json, que será usado na parte seis e você pode ignorar por enquanto.

Esses arquivos Nextflow são apenas o início de cada exercício. Se você quiser ver como eles devem ficar quando terminados, você pode ir para um diretório solutions e lá estão as respostas para cada parte do curso de treinamento, então você pode ver uma versão funcional do que você está almejando.

## Abrindo um terminal

Se a qualquer momento você fechar o terminal e não conseguir lembrar como voltar, não se preocupe com isso. Estes botões no topo, à direita, abrem e fecham diferentes painéis no workspace. Então clique neste para o painel inferior e ele reaparecerá. E apenas certifique-se de que você tem terminal selecionado aqui. Você também pode clicar neste botão aqui, a seta no lado direito de um terminal para deixá-lo em tela cheia.

Você me verá fazendo isso bastante porque tenho o VS Code ampliado para que você possa ler o texto. Dependendo do tamanho da sua tela, você pode ou não precisar fazer isso. O mesmo vale para minimizar o painel lateral.

Certo. Isso é suficiente para o ambiente. Acho que estamos prontos para começar. Junte-se a mim de volta no próximo vídeo para o capítulo um.

[Próxima transcrição de vídeo :octicons-arrow-right-24:](01_hello_world.md)
