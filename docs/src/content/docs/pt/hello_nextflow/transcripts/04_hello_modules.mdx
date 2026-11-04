# Parte 4: Hello Modules - Transcrição do Vídeo

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Tradução assistida por IA - [saiba mais e sugira melhorias](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

<div class="video-wrapper">
  <iframe width="560" height="315" src="https://www.youtube.com/embed/43Ot-f0iOME?si=0AWnXB7xqHAzJdJV&amp;list=PLPZ8WHdZGxmWKozQuzr27jyMGqp9kElVK&amp;cc_load_policy=1" title="YouTube video player" frameborder="0" allow="accelerometer; autoplay; clipboard-write; encrypted-media; gyroscope; picture-in-picture; web-share" referrerpolicy="strict-origin-when-cross-origin" allowfullscreen></iframe>
</div>

!!!note "Notas importantes"

    Esta página mostra apenas a transcrição. Para instruções completas passo a passo, retorne ao [material do curso](../04_hello_modules.md).

    Os números das seções mostrados na transcrição são fornecidos apenas para fins indicativos e podem não incluir todos os números de seção dos materiais.

## Boas-vindas

Olá, e bem-vindo de volta à parte quatro do Hello Nextflow. Esta seção é toda sobre módulos, e é uma seção bem curta do curso. Nós não vamos realmente escrever muito código, é mais sobre como organizamos o código em nosso pipeline.

Até agora, nós colocamos tudo em um único arquivo, o que é bom, e na verdade é assim que costumávamos construir pipelines Nextflow nos velhos tempos.

Mas conforme esse pipeline cresce, o script fica cada vez mais longo e mais difícil de navegar, manter, e também significa que não podemos realmente compartilhar nenhum código.

Os módulos Nextflow nos permitem extrair processos daquele script principal e então importá-los. Isso significa que o código é mais fácil de navegar e também significa que podemos compartilhar esse código de módulo entre diferentes pipelines.

Este pequeno diagrama na página principal da documentação mostra bem o conceito. Em vez de um script enorme, vamos incluir esses arquivos de módulo separados, de diferentes scripts de módulo, e tudo vai ser puxado para o fluxo de trabalho, mas ainda vai rodar exatamente da mesma forma.

Então vamos entrar no GitHub Codespaces e dar uma olhada. Como antes, limpei um pouco meu espaço de trabalho aqui. Removi os antigos diretórios Nextflow e o diretório work e assim por diante. Mas não importa se você ainda tem esses arquivos por aí.

Vou começar a trabalhar no arquivo hello modules, que está basicamente onde deixamos no final de um capítulo anterior. Temos nossos três processos aqui. Temos alguns parâmetros, o bloco workflow, onde estamos executando esses três processos e conectando-os com canais. Depois publicamos os canais de saída e temos o bloco output dizendo como publicar esses arquivos.

## 1. Criar um diretório para armazenar módulos

Agora, como eu disse, não vamos realmente escrever ou editar muito código. Vamos apenas mover o código que já temos. Os arquivos de módulo Nextflow normalmente têm um único processo neles, e por convenção normalmente os mantemos em um diretório chamado modules. Mas você pode chamar isso do que quiser. Mas vou manter um diretório modules no meu repositório aqui, e então vou criar um arquivo para cada processo. Então vou dizer novo arquivo, sayHello.nf.

## 2. Criar um módulo para sayHello()

Agora vou pegar meu processo e simplesmente vou selecionar este código, cortá-lo do arquivo principal hello modules e colá-lo aqui.

Obviamente isso não faz nada por si só. Nosso script principal ainda precisa desse processo, então precisamos trazê-lo de volta de alguma forma. E fazemos isso com a declaração include.

Então eu digito include e algumas chaves, e então pego o nome do processo. E digo from, e então dou um caminho de arquivo relativo. Então diz, começa com ./ porque é relativo de onde este script está salvo. Então é modules sayHello.nf.

Observe que a extensão do VS Code é bastante útil aqui. Ela nos diz, se consegue encontrar este arquivo e se consegue encontrar um processo, que estou nomeando. Se eu deliberadamente colocar um erro de digitação aqui, ele me dá um erro imediatamente e me dirá que não consegue encontrar este processo que estou tentando importar. Então fique de olho em quaisquer erros que você encontrar.

E isso é realmente tudo. Ainda temos nosso processo aqui. Não são necessárias mudanças aqui embaixo. O processo tem o mesmo nome e é executado exatamente da mesma forma. É só que o código real do processo agora está em um arquivo separado.

Podemos executar o fluxo de trabalho Nextflow novamente, vai funcionar exatamente da mesma forma. E isso é basicamente o resto deste capítulo do curso é apenas mover esses três processos para seus próprios arquivos.

Então vamos fazer isso agora. Vou rapidamente criar um novo arquivo de módulo para o segundo processo: convertToUpper.nf. Vou cortar esse código, colar lá. E então vou incluir esse. vamos lá, ótimo.

E então vou criar um novo arquivo para collectGreetings.nf. Cortar isso.

Muito cortar, cortar e copiar e colar.

E agora nosso script de fluxo de trabalho principal está de repente parecendo muito, muito mais curto, muito mais acessível e muito mais fácil de ler.

E você pode ver como o projeto agora começa a se construir com nossos diferentes arquivos. Podemos mergulhar nos detalhes nos lugares que queremos. Navegar pelo nosso caminho para encontrar etapas específicas no pipeline muito mais facilmente, e obter uma visão geral do que o pipeline está fazendo rapidamente.

## Navegando módulos com VS Code

Agora, claro, a desvantagem de fazer isso é que se você tem um pipeline grande, terá muitos arquivos de módulo e eles poderiam estar organizados em múltiplos subdiretórios ou todo tipo de coisa. Agora, novamente, uma pequena dica aqui. A extensão do VS Code é muito boa em navegar sua base de código para você e também em te contar sobre o código lá.

Você pode ver que o VS Code entende o que é este processo e me dá uma pequena visão geral dele quando passo o mouse, então posso ver sem ter que sair e encontrar o código fonte, quais são as entradas e as saídas, que é tipicamente a coisa mais importante quando estou usando em um fluxo de trabalho.

E também se eu segurar command, estou em um Mac, e clicar no nome do processo, ele abre o arquivo diretamente imediatamente. Puxa ele. Então posso simplesmente pular direto para lá sem nem pensar sobre quais são os caminhos reais do arquivo. E isso funciona em qualquer lugar, posso fazer isso também, onde quer que processos estejam sendo chamados. Então isso torna realmente rápido.

## 4.4. Executar o fluxo de trabalho

Ok, vamos apenas verificar que o pipeline ainda executa como esperamos. Então abra o terminal. Vamos fazer "nextflow run hello modules", e ver se executa sem problemas.

Esperançosamente, todo o ponto disso é que o pipeline está basicamente inalterado, então você não deveria, realmente ver nenhuma mudança de quando executamos antes. A saída aqui parece exatamente a mesma, e você pode ver nosso diretório results com todos os mesmos arquivos, então isso é ótimo. Nenhuma mudança é bom.

## Uma nota sobre nf-core/modules

Só antes de terminarmos, quero tocar rapidamente no poder da colaboração quando se trata de módulos. Esses arquivos estão no meu repositório, então não é óbvio de imediato como podemos colaborar neles. E há muitas formas diferentes que você pode fazer isso, mas provavelmente o maior e mais conhecido exemplo disso é o nf-core.

Se eu for ao site do nf-core, vou em resources, e modules. Você pode ver que o nf-core tem uma enorme biblioteca de módulos, quase 1700 módulos quando vejo isso. E então posso digitar o nome de qualquer uma das minhas ferramentas favoritas, ir e encontrar se alguém já escreveu um módulo para ela, e ver este processo de módulo pré-escrito aqui com todas as entradas, as saídas, os contêineres de software, toda essa informação, e você pode ver no lado aqui, quantos diferentes pipelines nf-core estão todos usando este único processo compartilhado.

Este é um exemplo meio extremo, mas você pode ver que isso está realmente reutilizando este código. E se eu clicar através para o código fonte no GitHub disso, é exatamente o mesmo que estamos fazendo. É apenas um grande processo em um arquivo.

Agora no lado do nf-core, fazemos alguns truques para poder compartilhar esses arquivos e trazê-los para diferentes repositórios. E se você quer saber mais sobre isso, vá conferir o curso que temos sobre usar e construir com nf-core especificamente. Mas eu queria apenas dar uma ideia de quão poderoso este conceito de reutilização de código pode ser.

## Encerramento

Certo, isso é tudo para módulos. Eu disse que era uma seção curta do curso. Confira o quiz, certifique-se de que entendeu e certifique-se de que tudo ainda está funcionando corretamente. E te vejo de volta no próximo vídeo, que é todo sobre contêineres de software. Muito obrigado.

I.
