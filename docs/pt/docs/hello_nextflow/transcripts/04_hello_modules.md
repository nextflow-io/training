# Parte 4: Olá Módulos - Transcrição

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Tradução assistida por IA - [saiba mais e sugira melhorias](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

<div class="video-wrapper">
  <iframe width="560" height="315" src="https://www.youtube.com/embed/Xxp_menS0E8?si=0AWnXB7xqHAzJdJV&amp;list=PLPZ8WHdZGxmXiHf8B26oB_fTfoKQdhlik" title="YouTube video player" frameborder="0" allow="accelerometer; autoplay; clipboard-write; encrypted-media; gyroscope; picture-in-picture; web-share" referrerpolicy="strict-origin-when-cross-origin" allowfullscreen></iframe>
</div>

!!!note "Notas importantes"

    Esta página mostra apenas a transcrição. Para instruções completas passo a passo, retorne ao [material do curso](../04_hello_modules.md).

    Os números de seção mostrados na transcrição são fornecidos apenas para fins indicativos e podem não incluir todos os números de seção dos materiais.

## Boas-vindas

Olá, bem-vindo à Parte Quatro do curso de treinamento Olá Nextflow.

Este capítulo se chama Olá Módulos, e vamos falar sobre como modularizar código Nextflow. O que vamos fazer é pegar nosso script de fluxo de trabalho único e dividi-lo em arquivos separados.

Isso torna o código mais fácil de navegar e manter conforme seu fluxo de trabalho cresce, e também torna possível compartilhar módulos entre pipelines, de modo que se você tiver múltiplos pipelines usando a mesma ferramenta, você só precisa escrever esse processo uma vez.

Um exemplo clássico disso é o repositório de módulos nf-core, que possui milhares de diferentes ferramentas em módulos prontos para uso, que você pode instalar e usar em seu fluxo de trabalho.

O Nextflow também pode trabalhar com sub fluxos de trabalho, que são como módulos, mas têm múltiplos processos. Isso está fora do escopo deste treinamento, mas funciona basicamente da mesma forma.

Certo. Vamos dar uma olhada.

Como sempre, comece indo para training.nextflow.io.

Vá para "Olá Nextflow" na barra lateral, e estamos fazendo a parte quatro: "Olá Módulos".

Agora vou pular para meu ambiente GitHub Code Spaces e dar uma olhada no arquivo "hello-modules".

Assim como antes, estamos começando no ponto final do capítulo anterior, então este script deve parecer familiar. Temos nossos três processos, say hello, convert to upper e collect greetings, e em um fluxo de trabalho simples, que executa esses três comandos e emite uma mensagem no final. Temos dois parâmetros chamados greeting e batch, que especifica o nome, que é usado para o arquivo de saída coletada no final.

## 0. Aquecimento: Execute hello-modules.nf

Podemos verificar que este fluxo de trabalho ainda funciona como esperamos fazendo nextflow run hello, modules.

Ótimo. Ele executou três tarefas com cada um desses processos, uma tarefa coletada, e nos disse que há três saudações neste lote. Se entrarmos em results, temos nossos diferentes arquivos de saída aqui, incluindo a saída do lote de teste coletado.

## 1. Crie um diretório para armazenar módulos

Certo. Vamos fazer alguma modularização.

É geralmente uma boa ideia colocar módulos em uma subpasta no repositório do seu pipeline, apenas para manter as coisas organizadas. Você pode chamar isso do que quiser, mas por convenção geralmente chamamos de modules.

Então vamos em frente, vamos para um terminal e fazemos make the modules. Você pode vê-lo aparecer na barra lateral e VS Code aqui.

## 2. Crie um módulo para sayHello()

Vou então criar um novo arquivo para meu primeiro módulo. Você pode fazer "touch" ou "code" ou você pode fazer isso na barra lateral, realmente não importa. Então vou fazer code modules e vou nomeá-lo após o processo. Então sayHello.nf. NF é uma extensão de arquivo tradicional para arquivos Nextflow.

Vou apertar salvar aqui e podemos ver nosso novo arquivo de módulo aparecer.

## 2.2. Mova o código do processo sayHello para o arquivo do módulo

Certo, a seguir vou pegar o código do módulo do fluxo de trabalho. Também vou pegar o hash bang aqui e copiar isso primeiro para que seja claramente um arquivo Nextflow. E então vou pegar este processo e vou cortar. Então vou removê-lo do meu script de fluxo de trabalho principal e vou colá-lo neste novo módulo.

Esse é todo o conteúdo que este arquivo de módulo vai conter. Apenas um único processo, sem fluxo de trabalho, sem lógica, apenas um processo sozinho.

Agora posso fechar este arquivo.

## 2.3. Adicione uma declaração de importação antes do bloco workflow

Agora meu fluxo de trabalho está faltando esse primeiro processo, então precisamos trazê-lo de volta importando-o. A sintaxe para isso é muito similar a outras linguagens de programação, então pode parecer familiar. Fazemos include chaves, o nome do processo, say hello, e então from o caminho do arquivo modules, say hello, nf. Fantástico.

Alguns truques aqui. A extensão VS Code é inteligente sobre isso. Ela reconhece este caminho de arquivo e você pode passar o mouse sobre ele e fazer follow link. Ou estou no Mac, posso fazer option click e ele abre este arquivo. Então podemos pular rapidamente para ele.

Este nome de processo agora está sendo usado pelo fluxo de trabalho aqui embaixo, e podemos fazer a mesma coisa aqui. Ele nos mostra um pouco de informação sobre esse processo, e novamente, posso segurar option, clicar nele, e ele vai abri-lo no editor.

Então é uma maneira realmente rápida quando você tem muitos arquivos para seus diferentes processos de navegar rapidamente pela sua base de código no VS Code.

Ok. Isso é basicamente tudo para este capítulo. Agora apenas fazemos a mesma coisa novamente para os outros processos.

## 3. Modularize o processo convertToUpper()

Então vamos criar um novo arquivo aqui. Chamá-lo de Convert to upper nf. Novamente, copiar o hash bang. E então cortar o processo.

Copiar o nome do processo ali, incluir uma nova declaração include com o novo nome do processo.

## 4. Modularize o processo collectGreetings()

E então fazer o mesmo para o terceiro processo. Novo arquivo, connect. Greetings,

fazer o hash bang. Cortar o processo, colar o processo, e fazer uma nova declaração include.

Agora você pode ver aqui que tenho um sublinhado de erro aqui dizendo invalid include source. E este é na verdade um erro genuíno que cometi porque estava me movendo um pouco rápido demais. Se você olhar com atenção, pode ver que perdi o T em convert to upper

Então o VS Code muito utilmente me disse que cometi um erro ali. Se eu corrigir esse nome de arquivo, o erro desaparece. É um bom exemplo de por que a verificação de erros dentro do VS Code é tão útil para escrever código Nextflow. Caso contrário, eu não teria percebido isso e só teria descoberto muito mais tarde quando tentasse executar o fluxo de trabalho.

Nosso script principal do pipeline agora está parecendo muito mais simples. Ele não tem nenhum processo, apenas temos três declarações include e nosso fluxo de trabalho. Não mudamos nenhuma lógica do fluxo de trabalho. Não mudamos nenhum código de processo, então esperançosamente deve funcionar exatamente da mesma forma.

## 4.4. Execute o fluxo de trabalho para verificar que ele faz a mesma coisa de antes

Vamos verificar. Vou abrir um terminal e vou executar exatamente o mesmo comando de antes.

Com certeza, ele executou nossos processos, say hello, convert to upper collect greetings, e nos deu três saudações novamente.

Então movemos nosso código, mas não mudamos nada sobre como o fluxo de trabalho executa e está completamente inalterado. A única diferença é que agora temos código mais limpo, mais fácil de manter e mais fácil de compartilhar com outros.

E é isso. Foi um capítulo curto. É um conceito simples, mas é muito poderoso e fundamental para como escrevemos fluxos de trabalho Nextflow mais complexos. Então é importante que você entenda isso e adquira o hábito de usá-lo.

No próximo capítulo, vamos ter uma mudança de ritmo e parar de pensar tanto sobre a sintaxe de escrever código Nextflow, e pensar um pouco sobre como usamos software nos processos em si. Junte-se a nós na parte cinco para Olá Contêineres.

[Próxima transcrição de vídeo :octicons-arrow-right-24:](05_hello_containers.md)
