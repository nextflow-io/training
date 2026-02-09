# Orientação - Transcrição do Vídeo

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Tradução assistida por IA - [saiba mais e sugira melhorias](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

<div class="video-wrapper">
  <iframe width="560" height="315" src="https://www.youtube.com/embed/PIjOdFaYwWA?si=y8lAedhEHWaTV4zd&amp;list=PLPZ8WHdZGxmWKozQuzr27jyMGqp9kElVK&amp;cc_load_policy=1" title="YouTube video player" frameborder="0" allow="accelerometer; autoplay; clipboard-write; encrypted-media; gyroscope; picture-in-picture; web-share" referrerpolicy="strict-origin-when-cross-origin" allowfullscreen></iframe>
</div>

!!!note "Nota importante"

    Esta página mostra apenas a transcrição. Para instruções completas passo a passo, retorne ao [material do curso](../00_orientation.md).

## Boas-vindas

Olá, e bem-vindo ao Hello Nextflow. Meu nome é Phil Ewels. Sou Gerente de Produto para Software de Código Aberto na Seqera, a empresa por trás do Nextflow.

Este curso é uma introdução prática à construção de fluxos de trabalho com Nextflow. Ele foi projetado para pessoas que são completamente novas no Nextflow e querem desenvolver seus próprios pipelines.

Os exemplos são todos de processamento de texto simples, para que você possa se concentrar nos conceitos do Nextflow sem precisar de conhecimento especializado do domínio, apenas alguma familiaridade com a linha de comando.

Vamos passar pelos fundamentos do Nextflow: escrever processos, conectá-los em fluxos de trabalho de múltiplas etapas, gerenciar dependências de software com contêineres e configurar pipelines para diferentes ambientes de computação. Ao final, você terá construído um pipeline funcional do zero.

Este curso se concentra em _desenvolver_ pipelines. Se você só quer _executar_ pipelines existentes sem mergulhar muito no código, temos um curso mais curto "Nextflow Run" que pode ser mais adequado para você.

Depois que você dominar os fundamentos aqui, também temos cursos de acompanhamento que aplicam esses conceitos à análise científica real. Vamos ensinar você a usar os pipelines e as melhores práticas da comunidade nf-core.

Se você ficar preso, vá para community.seqera.io. Há um fórum da comunidade ativo lá com uma seção dedicada apenas para perguntas sobre treinamento. Você pode usá-lo a qualquer momento, no entanto, também realizamos semanas de treinamento trimestrais com pessoas disponíveis especificamente para ajudar. Então, se você estiver fazendo o treinamento durante uma dessas, definitivamente não seja tímido e peça ajuda.

Você também pode tentar pedir ajuda ao Seqera AI. É ótimo para explicar código Nextflow e ajudá-lo com depuração.

Quando você estiver pronto para executar Nextflow em escala, Seqera Platform é o melhor lugar para fazer isso. Ele é executado em sua infraestrutura sem nenhum bloqueio de fornecedor, com tudo, desde lançamento de pipeline até monitoramento em tempo real, até ambientes de análise interativos. Mas por enquanto, vamos apenas focar nos fundamentos.

Certo, vamos começar.

## training.nextflow.io

Ok. A primeira coisa a notar é que todos os cursos de treinamento em training.nextflow.io são muito interativos. A ideia é que você siga o material de treinamento e minhas instruções, e passamos pelo material de treinamento juntos. Então você vai precisar de duas coisas: você vai precisar do seu laptop e você vai precisar deste site aberto. E é basicamente isso.

Então esta é a página inicial como ela aparece hoje quando eu gravo isso. Você pode ver que há uma visão geral das diferentes coisas, o histórico e os diferentes cursos que temos, que é uma lista que está crescendo o tempo todo.

Nextflow para iniciantes é onde estamos. Há dois cursos aqui, Nextflow Run, que é um curso diferente, e o Hello Nextflow, que é o que nos interessa.

E você também pode ver todos os diferentes cursos na barra lateral. Posso pular para Hello Nextflow, e podemos ver todos os diferentes capítulos que vamos trabalhar juntos.

Há algumas outras coisas importantes a notar aqui. Primeiro, o material de treinamento é versionado, então você pode ver aqui em cima. Diz 3.0 latest, que no momento em que estou gravando é a versão estável mais recente. Isso vai mudar com o tempo. Lançamos novos cursos e atualizamos o material ao longo do tempo. Então, se for 3.1 ou 3.2, não se preocupe muito. Se for 4.0, então provavelmente há um novo vídeo, e você deveria talvez ir e encontrá-lo porque provavelmente haverá atualizações significativas.

Outro menu suspenso no topo é este de idioma. Agora isso é novidade para a versão 3.0. Pegamos o material previamente traduzido, que foi feito por humanos, manualmente, e passamos isso para um LLM e configuramos toda essa nova infraestrutura para manter diferentes traduções do material de treinamento usando tradução LLM.

Então agora temos todas essas traduções fantásticas aqui. Então, se você quiser ouvir em coreano, pode carregar o site inteiro em coreano. E seguir por lá. O mesmo para todos esses outros idiomas, hindi e alemão e assim por diante. Vou seguir em inglês. Esse é o idioma principal no qual escrevemos o material.

Alguns outros botões se você gosta de ter o modo claro. Em vez do modo escuro, você pode seguir o site no modo claro aqui em cima.

E então também tudo que olhamos está em um único repositório GitHub, que é de código aberto, chamado nextflow-io/training. E se você clicar neste botão a qualquer momento, ele vai para o repositório GitHub. Voltaremos a isso em um minuto.

## Configurando GitHub Codespaces

Ok, então agora você tem isso aberto na aba do navegador. Vamos para Hello Nextflow e clicar. Você pode ver na página de introdução, ela nos diz alguns dos requisitos, a visão geral e o plano de aula do que vamos cobrir aproximadamente, e então vamos mergulhar em começar.

Existem diferentes maneiras de fazer este tutorial interativo. Se você estiver confortável, você é bem-vindo para fazer isso localmente em seu próprio computador com sua própria instalação do Nextflow. Se clicarmos em Opções de Ambiente, você pode ver que há mais detalhes sobre como fazer isso usando Devcontainers locais ou você também pode simplesmente instalar todo o software localmente, com instalação manual.

Estamos trabalhando para fazer isso funcionar bem com Seqera Studios, então essa é outra opção. Mas a mais comum agora é usar GitHub Codespaces.

Codespaces configura um ambiente sandbox em um servidor remoto executado pelo GitHub. E é gratuito para uma certa quantidade de uso, que geralmente é suficiente para treinamento. E ele vai configurar você com uma instância VS Code, um IDE onde você pode acessar todos os arquivos do repositório, executar Nextflow e tudo. E pré-configuramos Codespaces para você. Então ele tem tudo que você precisa.

A beleza disso é que é apenas um clique para configurar um Codespace. É o mesmo para todos, e sabemos que você já tem todos os pré-requisitos instalados, então é rápido e fácil.

Então a primeira coisa a fazer é ir para "Começando". Procure este botão, que diz _Abrir em Codespaces_. Vou fazer command \+ clique nele para abrir em uma nova aba, e isso nos leva ao GitHub.

É assim que parece. Podemos ver, definimos todas as opções aqui para você. Se você quiser, pode clicar em alterar opções. Algumas coisas que você pode fazer aqui. Você pode dar uma máquina de instância maior, por exemplo, se você achar que ela trava porque fica sem memória ou algo assim. Ou definir versões específicas do material de treinamento. Mas geralmente você pode apenas ir com o que configuramos aqui e você pode ver. Neste caso, está usando a versão 3.0.

Então vou clicar em criar novo Codespace. E isso me leva para dentro.

Note também, diz nenhum Codespace para retomar lá. Se eu criei anteriormente um Codespace, clicar nesse botão novamente no material de treinamento vai me levar para a mesma página e vai listar todos os Codespaces que já tenho em execução. Então você pode simplesmente voltar direto para eles e continuar de onde parou. Então não importa se você fechou seu laptop.

Eles se desligam automaticamente após alguns minutos de inatividade, mas não há problema. Você pode simplesmente reiniciá-los.

Depois de iniciar um novo Codespace, ele vai ficar nesta página assim e vai carregar por um bom tempo. Então agora é um bom momento para fazer uma pausa rápida. Talvez você tenha esquecido de ir ao banheiro ou queira uma xícara de chá antes de começarmos? Vá agora enquanto você está esperando por isso, porque vai ficar girando lá por um tempo.

Rapidamente enquanto esperamos carregar, também vou para github.com/codespaces e apenas mostrar esta é a página de visão geral onde você pode ver todos os diferentes Codespaces que você tem em execução atualmente.

Você pode ver que tenho um aqui para nextflow-io/training. Sem alterações, porque ainda não fiz nada nele. A quantidade de recursos que está usando, e você pode ver no momento está configurando. Posso ir aqui, clicar neste pequeno menu suspenso e clicar em excluir. Então, se você acidentalmente configurar vários Codespaces e não estiver usando alguns, pode excluir os antigos e limpar.

Finalmente, mais uma maneira de entrar nisso. Se formos ao repositório GitHub. E isso funciona para qualquer repositório GitHub. Clique em code. Você pode ter comandos para clonar o repositório localmente. E há uma aba chamada Codespaces. E novamente, você pode criar um novo, e pode ver os que já estão em execução.

Então novamente, se você esquecer como criou seu Codespace, você sempre pode voltar a ele desta forma.

## A interface do VS Code

Ok, a construção terminou e agora está começando a carregar o GitHub Codespaces. Nem sempre leva tanto tempo, então não se preocupe. É apenas a primeira vez que você cria o Codespace. Se você voltar para um que já existe, é muito mais rápido.

Não seja muito impaciente se esta é a primeira vez, ainda não terminou, mesmo que esteja começando a nos dar uma interface.

Mas enquanto esperamos as coisas finais serem configuradas, vou apenas levá-lo pela interface caso você não esteja familiarizado com o VS Code.

Primeiro, há a barra lateral de chat para coisas de IA, que não precisamos. Então vou fechar isso, me livrar disso e liberar algum espaço.

À esquerda, temos o explorador de arquivos que nos mostra todos os arquivos no repositório Git, que é o espaço de trabalho que criamos. Note, estes não são arquivos locais. Isso tudo está no servidor remoto onde estamos trabalhando. Você pode arrastar e soltar arquivos locais e coisas assim, mas na maior parte, não vamos pensar nisso hoje. Estamos apenas trabalhando puramente remotamente.

Existem outras ferramentas nesta barra lateral, por exemplo, pesquisa. Então você pode pesquisar todos os arquivos em um repositório de uma vez. E se estivéssemos fazendo trabalho de desenvolvimento no repositório de treinamento, poderíamos fazer integração com controle de origem com Git e depuração e outras coisas.

Outras coisas são, há uma janela principal de edição de código aqui em cima, que acabou de carregar uma prévia do readme, que é para o material de treinamento. Então, neste caso, está visualizando markdown, mas normalmente isso será um editor de código.

E então abaixo disso temos o terminal, que é onde vamos executar todos os nossos comandos e interagir diretamente com o Nextflow.

Tudo no Codespace está pré-instalado, então o comando Nextflow já está lá e assim por diante.

Ok. Quando você chegar até aqui, deve estar quase pronto. Você pode ver agora que baixou o servidor de linguagem Nextflow e configurou algumas extensões para nós no VS code, incluindo a extensão Nextflow, que vai ser útil. Então posso fechar isso e posso fechar o README.md.

E agora você pode ver que tenho mais alguns no lado esquerdo. Estou um pouco ampliado aqui, mas se eu diminuir o zoom você pode ver que um dos botões diz Nextflow com o ícone do Nextflow. E isso tem algumas coisas legais aqui para explorar o projeto e coisas assim, que voltaremos mais tarde.

Ok. Caso você perca algum desses painéis, esses botões no canto superior direito são realmente úteis e eles apenas mostram e escondem coisas. Então isso mostra e esconde o Explorer, mostra e esconde o terminal na parte inferior. E assim por diante.

Vou usar bastante isso porque estou muito ampliado, então tento ajudá-lo a ver todo o texto na minha tela, e então é útil poder ir em tela cheia com o terminal e depois escondê-lo quando estivermos olhando o código. Mas na maior parte do tempo você pode simplesmente ter todas essas coisas abertas ao mesmo tempo.

Ok, o que mais olhar? Não muito mais. Note que Nextflow, como eu disse, está instalado. Então posso digitar "nextflow -version" e deve aparecer dizendo qual versão temos instalada.

Há algumas outras coisas instaladas aqui também. No final de cada capítulo, temos um conjunto de perguntas de quiz, por exemplo, no site. E você também pode fazer isso no terminal se quiser digitando quiz.

Há alguns outros atalhos de teclado que vou usar, só para você saber. Por exemplo, agora mesmo pressionei cmd\+K no meu Mac e isso limpou o terminal, para me livrar de toda a saída anterior. Então isso é bom para manter as coisas limpas. Se você me ver fazendo isso, é assim que estou fazendo.

E também se você é novo no terminal, lembre-se que você pode usar tab para auto completar, o que vou fazer muito para auto completar caminhos.

Então posso ver no lado esquerdo aqui há uma pasta chamada Hello Nextflow, que é o que vamos trabalhar. Se eu fizer "ls" para listar arquivos, posso fazer "hel", apertar tab, auto completa. E então esta é uma maneira muito rápida de completar caminhos.

## Abrindo apenas a pasta Hello Nextflow

Ok. Isso é ótimo. Há muita coisa neste repositório, no entanto.

Há todos os arquivos para gerar o site, e há vários cursos diferentes aqui, e você pode fazer isso a partir desta raiz e apenas clicar na pasta "Hello Nextflow". Mas é bom realmente focar apenas nisso.

Você pode definir isso como seu espaço de trabalho com um monte de cliques por aqui e definindo um diretório de projeto e coisas assim. Mas a maneira mais fácil é digitar code, que é o comando CLI para iniciar o VS Code, e então "hello-nextflow".

Isso vai abrir uma nova aba do navegador e você pode fechar a antiga. E parece exatamente o mesmo. Mas agora você pode ver que estamos neste subdiretório e todos os outros arquivos são invisíveis, e temos uma configuração mais limpa.

Você pode ver aqui que também o diretório de trabalho atual agora está dentro da pasta Hello Nextflow. Então limpo e organizado. Não precisamos nos preocupar em estar no lugar errado. Ok.

## Nova Sintaxe do Nextflow para 2026

Há uma coisa especial que preciso mencionar neste ponto. Agora, no início de 2026, estamos começando a trazer recursos diferentes para o Nextflow, e um dos grandes novos é um novo analisador de sintaxe de linguagem dentro do Nextflow.

Basicamente o motor que lê seus arquivos Nextflow e entende isso, para tempo de execução. Há algumas mudanças na sintaxe, e é realmente importante que você use o Nextflow com o analisador de sintaxe correto habilitado.

Duas coisas que você precisa para isso. Você precisa de uma versão atualizada do Nextflow e você precisa ter certeza de que está habilitado.

Se eu fizer "nextflow -version" novamente, você verá que o Codespaces está executando com 25.10.2 e 25.10 é a versão mínima para poder usar essas coisas.

Se você estiver usando 26.04, que para mim ainda não saiu, mas vai sair em breve. Então isso estará executando o novo analisador de sintaxe por padrão, e você não precisa fazer mais nada.

Mas se você estiver executando 25.10, você precisa habilitar o analisador de sintaxe estrito, como é chamado, ou analisador de sintaxe v2.

Isso é feito com uma variável de ambiente. Já está definido no Codespaces, então você não precisa fazer nada. Mas se você estiver executando localmente, você precisa definir isso, e posso verificar isso fazendo "echo $NXF_SYNTAX_PARSER", e deve estar definido como v2.

Então, se você estiver executando localmente, apenas faça "export NXF_SYNTAX_PARSER=v2". Simples assim. Mas lembre-se de fazer isso. Porque caso contrário você vai ver algumas discrepâncias estranhas e erros conforme avançamos.

Se você tiver alguma dúvida sobre qualquer uma dessas coisas sobre versão do Nextflow e analisador de sintaxe, primeiro, lembre-se, você não precisa se preocupar se estiver no Codespaces. Tudo deve estar configurado corretamente. Mas segundo, se você for ao material de treinamento do Nextflow, se você descer, falar sobre requisitos de versão, há um link aqui que leva você para a página de ajuda sobre explorar versões, e isso meio que passa por tudo em detalhes.

Vale a pena ler isso se você tiver um momento. Porque ajuda a esclarecer quais são alguns dos diferentes termos, que você pode ouvir quando começar a usar o Nextflow. Coisas como DSL1, DSL2, analisador de sintaxe um, analisador de sintaxe dois, e assim por diante. Então vale a pena apenas dar uma olhada nisso e isso repete um pouco do que acabei de dizer.

Também é realmente útil se você escreveu código Nextflow anteriormente e está voltando para uma atualização. Ele diz algumas das coisas que mudam e vincula você a partes da documentação do Nextflow, que diz como atualizar seu código Nextflow.

## Arquivos do curso

Ok. Última coisa para nos familiarizarmos é apenas ver os arquivos, que estão neste diretório. Você pode olhar na barra lateral ou frequentemente no material de treinamento, usamos o comando tree, -L, que é o número de níveis para olhar. Vamos dizer dois, e se eu fizer isso em tela cheia, você verá que isso basicamente espelha exatamente o que vemos na barra lateral lá, mas exclui arquivos ocultos, que começam com um ponto.

Então os arquivos \*.nf, significa Nextflow. Então estes são os arquivos de script Nextflow, e há um arquivo inicial aqui para cada um dos diferentes capítulos do material de treinamento, que vamos abrir e explorar e depois editar.

Vamos mudar esses arquivos conforme avançamos, e então no final de cada capítulo, os arquivos devem estar parecendo praticamente iguais ao início do capítulo para o próximo. Mas damos a você esses arquivos diferentes para que você possa sempre começar do zero e não se preocupar muito em bagunçar a sintaxe.

Se você precisar comparar com algo que definitivamente deveria funcionar. Você pode verificar na pasta solutions, e este é como um estado final para cada um dos capítulos, então você pode comparar o que você escreveu com o que está lá.

Há um diretório data. Este tem apenas um arquivo greetings.csv, que usaremos como exemplo, dados de entrada em alguns dos cursos, e coisas como um arquivo de configuração e alguns parâmetros, que descreveremos mais tarde no curso.

## Conclusão

Ok, então agora esperançosamente tudo está funcionando. Sua tela parece igual à minha e você entende como acessar tudo e quais são todos os diferentes arquivos.

Se você rolar para baixo até o final da página em começando, pequena caixa de seleção você deve dizer que eu entendo o que estou fazendo. Meu ambiente está funcionando e você está definido, você está no diretório de trabalho corretamente para a pasta "Hello Nextflow".

Se você marcou todos esses e eles parecem verdes. Podemos continuar para o próximo vídeo e o próximo capítulo, que é a parte um. Hello World. Vejo você em um momento.
