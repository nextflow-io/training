# Parte 3: Olá Workflow - Transcrição

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Tradução assistida por IA - [saiba mais e sugira melhorias](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

<div class="video-wrapper">
  <iframe width="560" height="315" src="https://www.youtube.com/embed/zJP7cUYPEbA?si=Irl9nAQniDyICp2b&amp;list=PLPZ8WHdZGxmXiHf8B26oB_fTfoKQdhlik" title="YouTube video player" frameborder="0" allow="accelerometer; autoplay; clipboard-write; encrypted-media; gyroscope; picture-in-picture; web-share" referrerpolicy="strict-origin-when-cross-origin" allowfullscreen></iframe>
</div>

!!!note "Notas importantes"

    Esta página mostra apenas a transcrição. Para instruções passo a passo completas, retorne ao [material do curso](../03_hello_workflow.md).

    Os números de seção mostrados na transcrição são fornecidos apenas para fins indicativos e podem não incluir todos os números de seção nos materiais.

## Boas-vindas

Olá, bem-vindo à parte três do curso de treinamento "Olá Nextflow".

Este capítulo é chamado "Olá Workflow".

No capítulo dois, construímos um fluxo de trabalho simples de um processo, mas na realidade, pipelines são úteis porque podem encadear múltiplas etapas de análise juntas.

Neste capítulo, vamos pegar aquele exemplo inicial e estendê-lo para ser um pouco mais realista.

Vamos adicionar algumas etapas adicionais e veremos como usamos canais para conectar essas etapas.

Vamos olhar para múltiplas tarefas, que podem colapsar em um único processo e vamos olhar para processos que podem ter múltiplas entradas e múltiplas saídas.

Ok, vamos começar.

Então vamos começar. Como antes. Vamos para training.nextflow.io. Olá Nextflow, capítulo três. Olá Workflow. E vamos abrir nosso workspace. Eu limpei todos os meus arquivos de trabalho dos capítulos anteriores e vou abrir Olá Workflow.

Agora este é o mesmo arquivo com o qual temos trabalhado até agora, então isso deve parecer familiar. Temos nosso processo say hello. Temos nosso params.greeting com seu arquivo greetings CSV, e temos nosso fluxo de trabalho na parte inferior, que carrega aquele arquivo CSV, cria o canal e o passa para nosso processo.

## 0. Aquecimento: Execute hello-workflow.nf

Se você quiser, podemos experimentar isso e verificar se está funcionando como esperamos. Abra um terminal para nextflow run hello workflow nf e clique em enter.

Ok, ótimo. Nossos três processos executaram. Temos nosso diretório results com nossas três saídas. Bonjour. Hello. Holà. Então vamos fechar esses arquivos, fechar o terminal, voltar ao script.

## 1. Adicione um segundo passo ao fluxo de trabalho

Ok. Para nosso exemplo, estamos mantendo o básico e tentando ser agnósticos de domínio. Então nosso segundo processo vai apenas manipular essas strings, essas palavras, de uma maneira simples. Vamos usar o comando Unix translate para pegar esses arquivos e torná-los todos em maiúsculas. Fazemos isso com o comando "tr".

## 1.1. Defina o comando de conversão para maiúsculas e teste-o no terminal

Podemos tentar isso apenas no terminal bash, e ver se funciona. Então você faz echo, Hello World, e então passa isso com o caractere pipe para tr, e damos a ele um padrão de reconhecimento, a to z e para o que deve traduzir. A to Z em maiúsculas.

Isso é muito simples porque está literalmente fazendo os caracteres A a Z. Então não funcionará em nada que seja acentuado ou algo assim. Mas para os propósitos do exemplo, você deve entender a ideia.

Vou pressionar enter e ele imprime no terminal, HELLO WORLD em maiúsculas. E assim como antes, poderíamos redirecionar isso para um arquivo se quiséssemos. Outfile.

Ok. Vamos limpar isso.

## 1.1. Escreva a etapa de conversão para maiúsculas como um processo Nextflow

Vamos voltar ao nosso script e escrever um novo processo para lidar com este comando bash. Vou copiar o processo anterior, colá-lo abaixo, e chamá-lo de convert to upper. Para maiúsculas. Vou usar o mesmo publishDir results, mas vou fazer algumas mudanças aqui. Em vez de receber um val, vou receber um path input file, e vou ter um prefixo aqui upper, para que nossos arquivos de saída não sobrescrevam a saída. E vou usar o nome da variável da entrada. E então vou mudar o script aqui embaixo, e em vez disso vou usar cat no arquivo de entrada e assim como fizemos em Bash TR, a-z, upper input file .txt. Ok, vamos clicar em salvar.

## 1.2. Adicione uma chamada ao novo processo no bloco workflow

Agora se eu rolar para baixo, precisamos realmente chamar este processo. Apenas adicionar o processo em um script não é suficiente. Temos que dizer ao Nextflow que precisamos executar este processo e onde fazer isso.

Então vou fazer aqui, convert to upper e

ok, estamos recebendo um erro aqui dizendo que espera um argumento. Com certeza, precisamos passar algo para este processo para que ele realmente tenha algo para fazer.

## 1.3. Passe a saída do primeiro processo para o segundo processo

O que vamos fazer é pegar a saída deste processo. Então eu pego o nome, say hello, e quando eu faço dot out.

Para um exemplo simples como este, onde temos um processo que tem apenas uma saída e estamos passando isso para um novo processo, então ele tem uma entrada, isso deve ser tudo que precisamos. Então vou clicar em salvar, abrir o terminal, e vamos tentar executar isso novamente.

## 1.4. Execute o fluxo de trabalho novamente

Agora, eu não limpei meu diretório work da última vez que executei este fluxo de trabalho. Vou executá-lo novamente e vou usar isso como uma oportunidade para mostrar como o cache parcial funciona. Então se eu fizer traço único resume. Esperançosamente deve reutilizar as saídas daquele primeiro processo, que eram exatamente as mesmas da última vez que executei. Mas agora temos um novo processo aqui que não foi executado antes, que executa do zero. E com certeza, você pode ver que o primeiro processo usou as saídas do cache, e a segunda saída executou três de três. Você também pode ver que temos ambos os nossos processos aqui agora, nosso primeiro processo, say hello, executou três vezes, e nosso segundo processo convert to upper executou três vezes.

Se eu executar isso novamente, como lembrete, com -ansi-log false, devemos ver que seis tarefas de processo diferentes executam três para cada uma delas. Então isso está fazendo exatamente o que esperávamos. O primeiro processo está executando três vezes, passando essas saídas para um segundo processo, que está então executando três vezes.

Então vamos dar uma olhada dentro do diretório work e ver como o Nextflow está manipulando essas entradas de arquivo. Se eu pegar este diretório hash aqui do segundo processo, podemos usar um comando tree novamente com -a apenas para olhar esses arquivos. Você pode ver aqui que temos nosso arquivo de entrada, que é o arquivo Bonjour-output.txt, e isso é na verdade um symlink. É isso que esta seta está nos mostrando, e está apontando para o arquivo no diretório work anterior.

Isso faz sentido. O Nextflow manipula a execução de cada tarefa em seu próprio diretório encapsulado, então é completamente auto-contido. No entanto, ele precisa fornecer os arquivos de etapas anteriores como entrada. Em vez de alcançar fora do diretório work para obter esses arquivos, o Nextflow os prepara no diretório work.

Se tivermos um sistema de arquivos compartilhado como aqui, ele faz isso usando um symlink para que não use nenhum espaço de arquivo adicional. Se usarmos armazenamento em nuvem com buckets em diferentes localizações, ele buscaria esses arquivos e realmente os copiaria para o diretório work.

Vamos dar uma olhada no arquivo command sh. Se eu fizer code work, command sh, você pode ver, com certeza, ele está acessando aquele arquivo do diretório local. Então tudo é muito auto-contido e limpo.

Também podemos verificar o diretório results e ter certeza de que esses arquivos foram emitidos adequadamente. E com certeza, em results, podemos ver todos os arquivos de saída do primeiro processo e todos os arquivos de saída do segundo. E eles estão todos em maiúsculas como esperávamos.

É aqui que o poder do Nextflow começa a brilhar. Com algum código muito mínimo, o Nextflow lidou com a execução em paralelo dessas tarefas com encapsulamento limpo dentro de diretórios work separados e preparação de arquivos de entrada e saída e publicação de arquivos, tudo automaticamente para nós, logo de cara. Então você pode ver como, à medida que escalamos essa complexidade de nossos fluxos de trabalho de análise, essa funcionalidade é realmente, realmente valiosa.

## 2. Adicione um terceiro passo para coletar todas as saudações

Ok. Essas etapas foram um-para-um. Tivemos uma saída do primeiro processo indo para uma entrada para o segundo processo. Em seguida, vamos falar sobre como coletar essas diferentes saídas em uma única tarefa de processo, que é novamente, uma coisa muito comum de fazer. Então vamos rapidamente abrir o terminal e fazer uma execução de teste disso.

## 2.1. Defina o comando de coleta e teste-o no terminal

Vou trapacear e copiar o código bash de exemplo do material de treinamento e apenas pressionar enter.

O que podemos ver aqui é que executamos este comando echo três vezes para três arquivos de saída diferentes, que posso ver aqui. E então usamos o comando cat para imprimir a saída de cada um desses três arquivos diferentes, e redirecionar isso para um único arquivo coletado.

E se eu fizer "cat COLLECTED-output", você pode ver que tem o conteúdo daqueles três arquivos diferentes, agora em um único arquivo.

## 2.2. Crie um novo processo para fazer a etapa de coleta

Então vamos ver se podemos replicar a mesma coisa dentro do nosso pipeline Nextflow.

Vamos rolar para cima e criar um terceiro processo. Vou copiar este anterior, e desta vez vou chamá-lo de Collect Greetings.

No terminal bash, chamamos de collected output txt. Então vou dizer o mesmo path output aqui. E vou fazer o redirecionamento aqui, então é salvo da mesma maneira.

Ok. Precisamos mudar o que acontece no início daquele comando, e precisamos pensar sobre o que é o arquivo de entrada aqui. Na verdade, este processo vai receber múltiplos arquivos de entrada. Vou manter path e vou mudar isso para uma nova variável chamada input files, plural.

Vou então novamente, dar cat neles como fizemos em nosso script bash. E vou usar a variável aqui.

Agora, você pode pensar que isso não funcionaria. Vimos anteriormente falhas onde um array de strings ou um array de paths foi passado para um processo e isso causou um erro. Mas na verdade, aqui o Nextflow vai lidar com isso automaticamente para nós da maneira certa. Ele vai pegar vários arquivos de entrada diferentes, e vai apenas imprimir os diferentes caminhos de arquivo aqui.

Claro que ajuda que o comando cat pode receber uma série de nomes de arquivo assim. Se eu estivesse usando um comando diferente que requer um argumento antes de cada caminho de arquivo ou algo assim, teríamos que ter um pouco mais de código aqui e lógica para poder lidar com a iteração desses caminhos de arquivo. Mas neste caso, deve apenas funcionar.

## 2.3. Adicione a etapa de coleta ao fluxo de trabalho

Ok, vamos descer para o fluxo de trabalho e adicionar nosso novo processo. Collect greetings. E novamente, vamos pegar a saída de convert to upper out. Vamos salvar isso.

Vamos tentar. nextflow run hello workflow.

Ok, o fluxo de trabalho executou, mas algo está um pouco estranho aqui. Temos três execuções da primeira etapa, que esperamos. Três tarefas para a segunda, mas também temos três tarefas no final quando esperávamos ter apenas uma única tarefa aqui mesclando todas as saídas.

Se formos ao nosso diretório results. Também vemos que o output coletado tem apenas um único valor em vez de todos os três. Isso é porque aquele arquivo de saída foi sobrescrito três vezes com três valores diferentes.

Isso faz sentido porque passamos uma saída para uma entrada aqui da mesma maneira que fizemos na etapa anterior.

## 2.4. Use um operador para coletar as saudações em uma única entrada

Então precisamos de um operador aqui para pegar este canal com três elementos e colapsá-los em um único elemento, para que aquele processo final execute apenas uma vez.

Para fazer isso, vamos usar o operador collect. Posso fazer isso diretamente dentro do fluxo de trabalho. Posso fazer .out e encadear em um operador aqui no final .collect.

Pressione salvar. E então para os propósitos deste treinamento, também vou fazer alguns operadores view como fizemos antes, para que possamos dar uma olhada neste canal antes e depois de usarmos o operador collect, para que possamos entender o que está acontecendo.

Vou pegar este canal, remover o collect e dot view greetings, e então vou duplicar esta linha, adicionar o operador collect. E mudar isso para after.

Isso é separado de onde estamos chamando isso, mas tudo bem porque estamos usando as mesmas chamadas de operador no mesmo canal de saída.

Ok, vamos salvar e vamos experimentar no terminal. Vou fazer nextflow run. Hello, workflow. Reexecutar nosso script.

Ok. Isso está parecendo melhor. Como antes, podemos ver os primeiros dois processos executarem três vezes e agora nosso processo final executou apenas uma vez.

Se olharmos o que foi impresso pelo operador view, aqui embaixo, dissemos before collect, que é esta saída aqui, e isso é impresso três vezes. E você pode ver que há um único path para cada um desses. E então after collect, você pode ver que temos este array de três paths. Então isso é como esperamos.

Ok, vamos verificar o arquivo results e ver se é o que esperamos desta vez. Com certeza, agora há três linhas no arquivo - isso concatenou com sucesso essas três saídas em um único arquivo de saída. Fantástico.

Ok, vou limpar e vamos para o próximo passo. E vou deletar essas declarações view apenas para manter as coisas limpas.

## 3. Passe mais de uma entrada para um processo a fim de nomear o arquivo de saída final de forma única

Ok. Até agora, todos os nossos processos receberam apenas uma única entrada. Agora vamos fazer um exercício onde adicionamos mais de uma entrada a um processo para ver como isso funciona. Para fazer isso, vamos usar este exemplo collect greetings.

Cada vez que executei o fluxo de trabalho, ele sobrescreveu aquele arquivo no diretório results, o que pode não ser o que queremos.

## 3.1. Modifique o processo coletor para aceitar um nome definido pelo usuário para o arquivo de saída

Então para este exemplo, vamos passar um parâmetro adicional para que possamos personalizar o nome do arquivo de saída.

Adicionar uma segunda entrada a um processo é muito simples. Apenas adiciono uma segunda linha no bloco input. Desta vez vai ser um value, em vez de um path, porque queremos passar uma string e vou chamá-lo de batch underscore name.

Agora posso usar esta variável no bloco script, e vou dizer collected dash dollar batch name.

Estou usando chaves aqui em torno do nome da variável. Isso é apenas para mantê-lo separado do resto da string, e provavelmente não é necessário neste caso, mas acho que torna mais fácil de ler.

Ok. Finalmente, lembre-se de atualizar o path de saída porque agora o nome do arquivo mudou, então vou fazer a mesma coisa e colocar o batch name na saída do path como esperado.

## 3.2. Adicione um parâmetro de linha de comando batch

Agora precisamos passar um batch name de algum lugar, e vou criar um segundo parâmetro para fazer isso para que possamos fazer isso na linha de comando quando executarmos o fluxo de trabalho.

Então vou fazer params batch name, e por padrão, vamos chamar isso de test batch. Agora posso usar esta variável de parâmetros especiais abaixo, onde chamamos o processo.

E com certeza o VS Code está nos dizendo que não há argumentos suficientes para este processo agora, e que ele espera uma segunda entrada.

Simplesmente faço vírgula e passo nossa nova variável e o erro desaparece.

Note que a ordem das entradas aqui é realmente importante. A primeira entrada do processo foi o path, e a segunda entrada é o nome. Se eu mudar a ordem aqui, também devo mudar a ordem quando chamo o processo. Caso contrário. Em seguida, vamos passar o canal errado para a entrada errada.

## 3.3. Execute o fluxo de trabalho

Ok, vamos tentar e ver se funciona. Vamos fazer "nextflow run hello- workflow. Ok, executou como antes. Vamos dar uma olhada no diretório results.

Com certeza, nosso nome de arquivo aqui agora é chamado "collected test batch output txt". Fantástico.

E agora vamos ver se podemos sobrescrever isso executando novamente. Desta vez vou fazer --batch_name para combinar com aquele nome de variável de parâmetro especial aqui. E vou chamá-lo de demo output.

Execute o fluxo de trabalho novamente e veremos se algo acontece.

Ok, agora temos um collected demo output .txt. E porque este nome de arquivo é diferente daquele, ele não o sobrescreveu. Ambos agora estão presentes no diretório results.

## 4. Adicione uma saída à etapa coletora

Ok, então mostramos dar múltiplas entradas a um processo, mas e quanto a múltiplas saídas? Para este exemplo, vamos calcular o número de saudações que são processadas e emitir isso como uma saída secundária para esta etapa collect greeting.

## 4.1. Modifique o processo para contar e emitir o número de saudações

Vamos fazer um truque aqui. Os processos Nextflow têm este bloco script com uma string multi-linha, e isso é passado como saída bash para o ponto comando ponto sh. Mas na verdade podemos escrever qualquer código personalizado acima disso, e isso será executado como parte de uma tarefa mas não incluído dentro do script bash.

Uma das funções integradas na sintaxe Nextflow é chamada size. Então vou pegar a entrada path, e vou dizer count underscore greetings, apenas para definir um nome de variável. Vou pegar os input files e vou chamar "size" nele.

Esta função vai contar o tamanho deste canal de entrada e atribuir a uma variável.

Agora podemos retornar aquela variável como parte do bloco output. Então dizemos, val, porque é valor, não um arquivo. E count greetings.

Agora isso é suficiente por si só, e agora poderíamos acessar essas diferentes saídas deste processo. No entanto, teríamos que acessá-las de maneira posicional. Então usando uma chave de índice como zero e um.

Para tornar um pouco mais fácil obter as saídas, podemos nomeá-las e fazemos isso usando uma declaração emit.

Então fazemos vírgula emit out file ou como eu quiser chamar isso. E faço aqui emit count. Isso é basicamente apenas um decorador, que apenas nos ajuda a escrever código um pouco mais limpo para que possamos facilmente referenciar as saídas específicas mais tarde no bloco workflow.

## 4.2. Relate a saída no final do fluxo de trabalho

Ok. Se eu rolar para baixo para o bloco workflow, agora posso pegar as saídas de collect greetings, fazer collect greetings, dot out, e podemos ver nossas duas saídas nomeadas são sugeridas aqui pela extensão VS Code. Muito útil.

Então vou fazer dot count para obter o valor de contagem que acabamos de criar, e vou fazer view, para que ele imprima na linha de comando. Então podemos vê-lo quando executarmos o fluxo de trabalho.

Vamos escrever algo no closure aqui apenas para torná-lo um pouco mais agradável. num greetings, havia greetings saudações.

E na verdade não nos importamos com a outra saída porque não estamos usando isso como entrada para quaisquer outros processos. Mas você pode ver como poderíamos facilmente passar isso como entrada para outro processo se quiséssemos, downstream.

## 4.3. Execute o fluxo de trabalho

Vamos clicar em salvar. Vamos dar uma olhada no terminal e tentar isso.

Ok, fantástico. Aqui vamos nós. Há três saudações. Isso está exatamente certo.

Ok, ótimo. Esse é o fim deste capítulo. Terminamos por ter chegado até aqui. Você está agora começando a construir um fluxo de trabalho bastante realista, onde somos capazes de lidar com entradas e saídas e lógica dentro do nosso fluxo de trabalho.

À medida que esses arquivos de fluxo de trabalho ficam mais longos, eles começam a se tornar um pouco difíceis de manejar. Então no próximo capítulo, vamos olhar para como podemos modularizar o código Nextflow em arquivos separados para que seja mais fácil encontrar e manter o código dentro do fluxo de trabalho.

Junte-se a nós no próximo vídeo para o capítulo quatro. Olá Modules.

[Próxima transcrição de vídeo :octicons-arrow-right-24:](04_hello_modules.md)
