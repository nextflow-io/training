# Parte 3: Hello Workflow - Transcrição do Vídeo

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Tradução assistida por IA - [saiba mais e sugira melhorias](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

<div class="video-wrapper">
  <iframe width="560" height="315" src="https://www.youtube.com/embed/_aO56V3iXGI?si=Irl9nAQniDyICp2b&amp;list=PLPZ8WHdZGxmWKozQuzr27jyMGqp9kElVK&amp;cc_load_policy=1" title="YouTube video player" frameborder="0" allow="accelerometer; autoplay; clipboard-write; encrypted-media; gyroscope; picture-in-picture; web-share" referrerpolicy="strict-origin-when-cross-origin" allowfullscreen></iframe>
</div>

!!!note "Notas importantes"

    Esta página mostra apenas a transcrição. Para instruções completas passo a passo, retorne ao [material do curso](../03_hello_workflow.md).

    Os números das seções mostrados na transcrição são fornecidos apenas para fins indicativos e podem não incluir todos os números de seção dos materiais.

## Boas-vindas e recapitulação

Olá, e bem-vindo de volta à parte três de Hello Nextflow. Esta parte se chama Hello Workflow, e é nesta parte do curso onde realmente começamos a justificar o nome pipeline ou fluxo de trabalho.

Vamos pegar nosso script de pipeline simples até agora com seu único processo, e vamos começar a adicionar processos adicionais e ver como Nextflow lida com essa orquestração e o fluxo de dados através do pipeline.

Vamos voltar aos nossos code spaces. Você verá que eu deletei todos os diretórios .nextflow\* e os diretórios work e tudo para tentar mantê-lo limpo. Não se preocupe se você ainda tem esses arquivos das partes anteriores do curso.

Vamos trabalhar a partir de um arquivo chamado hello-workflow.nf. Como antes, isso basicamente representa o script que construímos até este ponto, e nos dá um ponto de partida limpo. E novamente, lá embaixo na saída podemos ver que o caminho agora é hello_workflow. Então os arquivos publicados devem estar indo para um subdiretório diferente na sua pasta results.

Para recapitular onde estamos até agora, temos um único processo aqui, com uma entrada greeting, uma saída greeting file. E então o script Bash simples, que apenas faz um comando echo para um arquivo.

Temos uma única entrada de fluxo de trabalho, o bloco params aqui, onde dizemos que está esperando um caminho, e o padrão é data/greetings.csv, que é este arquivo aqui em cima.

Então no fluxo de trabalho em si, temos um bloco main. Estamos criando um canal. Estamos analisando o CSV em linhas e então pegando o primeiro elemento de cada array, e estamos passando esse canal para aquele processo, que está então gerando três tarefas, e estamos publicando a partir do fluxo de trabalho, as saídas daquele processo.

E então finalmente, no bloco output, estamos dizendo ao Nextflow para publicar esses arquivos deste canal para o diretório chamado hello_workflow. E para copiar esses arquivos em vez de criar soft links para eles.

## 1. Adicione uma segunda etapa ao fluxo de trabalho

Okay, nesta parte vamos adicionar um segundo processo ao nosso fluxo de trabalho. Vamos pegar as saídas do processo sayHello, e processá-las em uma segunda etapa, que vai converter todas as letras dentro daqueles arquivos convertToUppercase.

Este é apenas um exemplo bobo, é apenas algum processamento simples de string novamente, mas mostra como podemos pegar a lógica, dentro do fluxo de trabalho.

Vamos usar um comando bash chamado "tr" para isso, que é a abreviação de translate. É um comando Unix que existe há muito tempo. Se você não está familiarizado com ele, não te culpo. Acho que nunca usei antes do treinamento, mas você pode testá-lo muito rapidamente no terminal. Se eu fizer "echo 'hello world'" e então pipe para 'tr' e então entre aspas você diz intervalo de caracteres, então A até Z, minúsculo, e então você quer fazer A até Z maiúsculo. E ele apenas diz, traduza estas letras para estas letras.

E quando eu apertar enter, você pode ver que agora colocou tudo em maiúsculas. Muito bom se você gosta de gritar com as pessoas.

Então esse é um estilo muito simples de comando bash que vamos usar em nosso segundo processo.

## 1.2. Escreva a etapa de conversão para maiúsculas como um processo Nextflow

Então se eu voltar ao meu script, vou trapacear um pouco e apenas copiar o código dos, dos documentos do treinamento. Mas você pode ver exatamente o que está acontecendo.

Temos um novo processo aqui. Este nós chamamos de convertToUpper, mas poderíamos chamá-lo do que quiséssemos.

Temos uma única entrada path, como fizemos antes. Não é um canal de valor, é um canal de caminho. E então uma única saída.

No bloco script fazemos "cat" no arquivo de entrada. E podemos colocar isso entre chaves se quisermos. e que pega aquela variável. E executamos o mesmo comando bash no pipe e escrevemos os resultados para um arquivo com este nome de arquivo, e isso é capturado pelo output path.

Agora precisamos fazer algo com este novo processo. Então vamos descer até o fluxo de trabalho onde construímos a lógica diferente de um fluxo de trabalho, e após aquele primeiro processo, vamos executar nosso segundo processo. Então convertToUpper é o nome do processo aqui.

Ele recebe uma entrada, então não podemos apenas chamá-lo sozinho. Queremos processar a saída do primeiro processo. Então assim como fizemos com isso, sayHello out onde estamos publicando esses resultados. Queremos usar esses mesmos resultados aqui como entrada, então podemos copiá-los e colocá-los lá.

Queremos o processo sayHello ".out", e Nextflow sabe que isso significa um registro de saída único simples aqui, que é este arquivo. Então isso será passado como entrada para um segundo processo.

## 1.5. Configure a publicação de saída do fluxo de trabalho

Okay. E finalmente, para que realmente salvemos os resultados deste segundo processo, também precisamos publicá-los a partir do fluxo de trabalho, e então defini-los no bloco output, mesma sintaxe de antes. Então podemos copiar isso e dizer second outputs, ou o que você quiser chamar.

Pegar o nome do processo que estamos interessados, convertToUpper out, e então aqui embaixo no bloco output. Adicionar isso e poderíamos fazer os mesmos atributos aqui. Então também queremos esses arquivos no subdiretório Hello Workflow, e também queremos copiá-los.

Ótimo. Vamos tentar executá-lo. Então se eu abrir o terminal e fizer "nextflow run hello-workflow.nf", e veremos o que ele faz. Veja se parece diferente das partes anteriores.

Então ele lança o Nextflow. Nos documentos, diz para fazer isso com "-resume", mas eu deletei todo o meu diretório work, então não teria feito diferença aqui. Mas se você fez, então isso funcionará também.

E parece quase exatamente o mesmo. Mas você pode ver agora que há uma segunda linha de saída aqui, onde você pode ver o nome do segundo processo que acabamos de adicionar. E com certeza, você pode ver que executou três vezes com sucesso.

Brilhante. Se eu tivesse meus diretórios work anteriores por perto e tivesse feito isso com "-resume", estes teriam sido, em cache apenas a primeira etapa no pipeline. Porque aquelas saídas eram exatamente as mesmas, então o Nextflow saberia para reutilizá-las novamente.

E então você pode ver como pode usar -resume para construir iterativamente seu fluxo de trabalho, passo a passo, se precisar.

Okay, vamos dar uma olhada no diretório results aqui em cima e ver se funcionou. Podemos ver que temos alguns arquivos a mais aqui em cima. Temos nossos arquivos originais como fizemos antes do primeiro processo. E com certeza, temos nossos arquivos upper e as letras estão todas em maiúsculas, então funcionou. É realmente bom de ver.

Também é interessante apenas verificar dentro desses diretórios work. Como antes, o hash aqui corresponde aos diretórios work. Então se eu olhar em "ls work", e então expandir aquilo, veremos os diferentes arquivos aqui.

Vemos o arquivo de saída do primeiro processo, que foi puxado aqui como entrada. E podemos ver o novo arquivo de saída que foi gerado.

Agora se eu fizer isso com "-la" para listar e mostrar todos os arquivos, veremos mais algumas coisas. Primeiro, você verá que este arquivo é na verdade um soft link para o primeiro processo. Isso é basicamente sempre um soft link se puder ser, para economizar espaço de arquivo. Não estamos publicando os arquivos aqui e apenas referencia aquele arquivo de uma primeira tarefa para uma segunda tarefa de modo que tudo esteja encapsulado dentro daquele diretório de trabalho, e seguro e isolado de todo o resto.

E isso precisa estar lá porque se olharmos o arquivo .command.sh, então se eu fizer "cat work/b8/56\*", você pode ver que as partes de arquivo aqui são relativas, então está fazendo cat daquele arquivo de entrada, que foi soft linkado para o mesmo diretório de trabalho.

Então é assim que todo diretório work parecerá. Quando você olha no Nextflow, você terá todos os arquivos de entrada lá preparados naquele diretório work. E então você também terá quaisquer arquivos de saída que foram criados. Então isso é ótimo. Isso parece como esperávamos.

## 2.1. Defina o comando de coleta e teste-o no terminal

Okay, vamos voltar ao nosso fluxo de trabalho. Qual é a próxima etapa que queremos fazer?

Temos dois processos agora e eles estão pegando este arquivo CSV, analisando e dividindo. E então temos três tarefas para cada um desses processos e Nextflow lida com a paralelização de tudo isso, então tudo executa lado a lado quando possível.

Essa forma de dividir o trabalho para executar coisas em paralelo é muito comum. E o inverso disso é então reunir tudo de volta. Então é isso que vamos fazer com nosso processo final no fluxo de trabalho, teremos um terceiro aqui, que pega essas três saídas diferentes e combina todas elas em um único arquivo.

Podemos fazer isso de forma bem simples em um terminal, apenas para ter uma noção de como isso parecerá.

Se eu for para a pasta results. Então, "cd results/hello_workflow/", e temos todos os arquivos UPPER aqui. Posso apenas usar "cat", que usamos para imprimir o conteúdo daquele arquivo, e você pode dar múltiplos arquivos para "cat" e ele lerá um após o outro.

Então posso dizer "UPPER-\*", o que me dá a mesma lista de três nomes de arquivo com expansão Bash. E posso dizer combined.txt. Acho que nos documentos, lista os nomes exatos dos arquivos, mas está fazendo a mesma coisa.

Agora, se eu usar "cat combined.txt", podemos ver que temos o conteúdo dos arquivos dos três daqueles arquivos.

Então isso é basicamente tudo que este processo vai fazer é vamos tentar dar a ele todos os diferentes arquivos de saída de um processo anterior em uma única tarefa de processo, e então vamos fazer "cat" deles juntos e salvar o arquivo de saída.

## 2.2. Crie um novo processo para fazer a etapa de coleta

Okay, então vamos adicionar nosso novo processo. Vou colar isso dos materiais de treinamento, e você pode ver que nos deixou um pouco de exercício para o leitor aqui com essas interrogações. Mas você pode ver o esboço geral do processo é basicamente o que acabamos de fazer no terminal, onde estamos fazendo "cat" de um monte de arquivos de entrada e escrevendo para um arquivo de saída aqui chamado collected, e então a saída espera aquele caminho único novamente.

Então precisamos de algum tipo de entrada aqui e vão ser um conjunto de caminhos. Então novamente, definimos um canal de entrada path e vamos chamá-lo de input_files. Agora, isso anteriormente nos deu um único caminho aqui, mas um caminho também pode ter múltiplos arquivos aqui, mesmo que ainda seja uma única declaração.

Vou copiar isso aqui embaixo porque queremos fazer "cat" desses arquivos. E você pode pensar que temos alguns problemas aqui com imprimir um array ou coisas assim, mas Nextflow é geralmente bem sensato quando se trata disso. E se for dado um canal com múltiplos arquivos nele assim, ele colocará todos juntos com separadores de espaço. Então isso nos dará a sintaxe correta.

Isso é ótimo. Então agora vamos conectar nosso novo processo. Vou descer até o fluxo de trabalho. Vou dizer combine the outputs, o novo nome do processo, e exatamente como antes. Vou pegar este processo anterior, convertToUpper e fazer ".out".

Ótimo. Vamos testar e ver se funciona no terminal. Se eu apenas voltar alguns diretórios e então executar novamente o comando Nextflow, e veremos o que acontece.

Então o fluxo de trabalho foi lançado e agora você pode ver que temos três nomes de processos diferentes, o que é ótimo. Os dois primeiros parecem os mesmos de antes, e o terceiro novo executa, o que é bom.

No entanto, há algo um pouco estranho aqui. Queríamos combinar aqueles arquivos de saída em um único arquivo, e ainda assim este processo podemos ver que executou três vezes, não uma.

Com certeza, se formos em um desses diretórios work. E fazer "cat work/" "collected", então veremos. Há apenas uma única palavra aqui, não três.

E então o que aconteceu é que Nextflow continuou aquela paralelização exatamente como fez nas etapas anteriores. E este processo nos deu um canal com três elementos, e aqueles três elementos do canal foram passados para aquele nosso processo downstream, que gerou três tarefas de processo.

Basicamente tentou coletar três vezes separadas e cada vez tinha apenas um único arquivo, então apenas fez cat de arquivo único para uma saída, e de fato, podemos ver isso no arquivo .command.sh também.

Se eu fizer .command.sh, podemos ver que tem apenas um único nome de arquivo aqui e apenas um único arquivo foi preparado naquele diretório de trabalho.

## 2.3. Adicione a etapa de coleta ao fluxo de trabalho

Então de alguma forma precisamos dizer ao Nextflow para juntar todas aquelas saídas de um processo anterior e dá-las para este processo downstream como um único elemento de canal, em vez de três.

Fazemos isso com um operador de canal chamado _collect_.

Este é um operador super útil, que você verá em pipelines Nextflow o tempo todo. Este é um canal aqui, este canal de saída, exatamente como aquele que criamos lá em cima. E então podemos anexar operadores de canal a ele exatamente como fizemos antes. Podemos apenas fazer ponto, e então neste caso, collect, parênteses.

E isso é tudo que precisamos. Isso vai então manipular este canal antes de ser passado para este processo.

Se você quiser ver o que está acontecendo com ele, também podemos visualizá-lo aqui. Então aqui, isso não está relacionado a executar este processo de forma alguma, então eu poderia colocá-lo em qualquer ponto após executar aquele processo. Mas pegamos o mesmo canal de saída, e estamos olhando com .view, e então estamos olhando novamente com .collect.view.

E quando executarmos isso, mostrará as duas estruturas diferentes daquele canal, antes e depois de collect. Então vamos tentar isso agora. Okay, apenas dei um zoom out um pouco porque algumas das saídas são bem longas, mas se eu executar o pipeline, veremos se funciona.

Estou esperando que o terceiro processo execute apenas uma vez, porque está coletando as saídas e com certeza, você pode ver collectGreetings como um de um. Então isso executou apenas uma tarefa.

E então se olharmos as declarações view, temos três declarações view para os três elementos de antes, com um caminho de arquivo em cada um.

E então após aquela declaração collect, isso apenas foi acionado uma vez porque há um único elemento naquele canal. E agora temos esta lista de três caminhos de arquivo diferentes.

Isso é exatamente o que esperávamos. E você pode ver, espero, isso é basicamente o inverso daquele operador "map" que fizemos para ir de arrays CSV para elementos de canal separados. Agora estamos pegando elementos de canal separados e colocando de volta em um único array.

Ótimo, podemos limpar essas declarações view. Não precisamos mais delas. Podemos ir para a próxima etapa.

Antes de ir mais longe, e antes que eu esqueça, vou adicionar uma nova declaração publish aqui. Third output. Você pode chamar isso de algo mais semântico e descritivo em seu fluxo de trabalho. E então vou adicionar isso ao bloco output novamente e dizer path 'hello_workflow' mode 'copy'. Apenas para que o arquivo de saída gerado por este processo seja salvo em nossa pasta results aqui em cima.

Apenas para verificar rapidamente que funciona. Deve estar um pouco mais limpo agora porque não temos aquelas declarações view. E veremos se conseguimos nosso novo arquivo de saída aqui em cima. Um de, uma tarefa executou, conseguiu um novo arquivo chamado collected, e agora temos todas as três palavras. Fantástico. O que vem a seguir?

## 3. Passe parâmetros adicionais para um processo

Okay. A seguir vamos olhar para lidar com múltiplas entradas em um único processo. Até agora você pode ver que todos os nossos processos estão apenas pegando uma coisa como entrada. Todos têm uma única linha sob sua entrada.

Vamos demonstrar isso permitindo que Nextflow especifique um identificador de lote diferente para que talvez você execute este fluxo de trabalho várias vezes e possa dar um ID de lote diferente cada vez.

Vou simplesmente adicionar uma segunda linha na entrada aqui para collectGreetings. E vou chamá-lo de "val", porque isso é uma string. Agora é um valor, não um caminho, e vou chamá-lo de "batch_name".

Então vou editar o script aqui embaixo para usar esta variável, e vou tentar colocá-la no mesmo lugar que o material de treinamento. Então coloquei no meio deste caminho de arquivo COLLECTED-$\{batch_name\}-output.

Ainda não terminei. Lembre-se que temos que dizer ao Nextflow quais serão os nomes dos arquivos de saída. Então também temos que fazer a mesma coisa aqui em cima: COLLECTED-$\{batch_name\}-output.txt".

Fantástico. Nextflow agora está recebendo uma segunda entrada de variável e está interpolando isso no script e na saída.

Uma última coisa, agora temos que encontrar onde isso está sendo chamado, e temos que passar a segunda entrada para o processo. Isso é como qualquer outra entrada em uma função em qualquer outra linguagem.

Assim como fizemos anteriormente no treinamento, vou usar o "params" especial aqui, e vamos chamá-lo de "params.batch" para que possamos ter uma opção CLI --batch. E agora você pode ver que nosso processo aqui tem duas entradas separadas apenas separadas por vírgula, que estão sendo passadas.

É realmente importante acertar a ordem, então a ordem dos argumentos aqui para canal e então o param devem corresponder. O canal e o batch name lá. Isso é apenas correspondência posicional.

Okay. Posso executar este pipeline agora imediatamente com --batch, mas vamos primeiro fazer a coisa certa e defini-lo na entrada aqui em Params. Então vou adicioná-lo ao batch e então vamos dizer que é uma string e vamos dar um padrão. Então vamos apenas chamá-lo de batch. Okay? Agora vamos tentar executar o fluxo de trabalho.

--batch Trio. Acho que diz no material de treinamento, mas poderíamos usar qualquer string que quiséssemos lá. E esperançosamente veremos aquele arquivo de saída de results aparecer aqui.

E com certeza, COLLECTED-trio-output - isso funcionou corretamente. Renomeou nosso arquivo. E você pode imaginar agora isso é útil porque se eu executar isso novamente com um nome de lote diferente, como replicate_two, então vai nos dar um nome de lote diferente aqui em cima.

E e não vai então sobrescrever os arquivos de saída neste caso. Então isso é bom.

## 4. Adicione uma saída à etapa coletora

Okay, então agora temos múltiplas entradas para nosso processo aqui. Mas o que acontece se quisermos criar múltiplas saídas? Nosso exemplo aqui então é que vamos criar um relatório para este processo, apenas dizendo este é quantos arquivos foram coletados.

E faremos isso com um comando echo aqui. Então podemos dizer echo. There were, vou copiar isso de um material de treinamento, para que você não tenha que me assistir digitá-lo.

There were $\{count_greetings\} greetings in this batch, e salvar isso em um novo arquivo agora chamado $\{batch_name\}, então mesma variável, podemos reutilizar isso quantas vezes quisermos, report.txt.

## 4.1.1. Conte o número de saudações coletadas

Precisamos realmente calcular isso de alguma forma. Poderíamos fazer essa lógica no script Bash se quiséssemos, usando lógica Bash. No entanto, também podemos apenas fazer script diretamente dentro do código Nextflow, contanto que esteja dentro do bloco script no processo e acima da seção entre aspas.

Qualquer coisa aqui não será incluída no script final renderizado, e será apenas executada pelo Nextflow quando renderizar uma tarefa.

Então aqui estamos apenas fazendo alguma lógica. Estamos criando uma nova variável chamada count_greetings. Pegamos o canal de arquivos de entrada aqui, e estamos chamando .size() nele.

Okay, aquela função vai me dar um número aqui nesta variável, e agora nosso aviso foi embora porque esta variável está sendo definida.

Okay, então estamos criando aquele segundo arquivo no diretório work, mas precisamos dizer ao Nextflow para esperá-lo como uma saída publicada deste processo. Então fazemos isso pela exatamente mesma sintaxe que fizemos para o primeiro arquivo.

Dizemos path porque, novamente, poderíamos estar publicando uma variável aqui se quiséssemos com "val", mas vamos dizer "path". E então o nome de arquivo esperado. Note que não está destacado aqui. Isso porque usei aspas simples. Tenho que usar aspas duplas.

## 4.1.2. Emita o arquivo de relatório e nomeie as saídas

Okay, isso é ótimo. E agora poderíamos começar a acessar essas saídas aqui embaixo exatamente como fiz aqui. Mas agora é um array de objetos diferentes, então eu poderia fazer collectGreetings.out[0] para pegar o primeiro, ou um para pegar o segundo, que é nosso novo relatório.

Mas eu realmente não gosto muito de fazer isso porque é muito fácil bagunçar a contagem de índice. E você fica lá contando linhas muito e você adiciona uma nova saída e de repente tudo quebra. Então

é muito mais legal referenciar tudo por nome. E podemos fazer isso com uma chave especial aqui chamada "emit".

Então podemos chamar isso do que quisermos. Vamos chamar de emit outfile, e emit reports. Se você definir estes e pode fazer em um ou muitos, é com você. Agora posso descer aqui e em vez disso posso ir dot out dot reports e apenas chamá-lo por nome, o que é muito mais fácil de entender seu código quando você o lê, e é mais seguro para mudanças no código.

Adicionei o .out.report aqui, mas na verdade preciso ter duas saídas diferentes sendo publicadas. Então vou renomear como algo mais interessante como collected e report e é assim que chamei? Chamei de out file, desculpa. Então aquele nome emit aqui outfile e report. porque estamos publicando dois canais de saída diferentes e então precisamos referenciar ambos no bloco publish.

Então também precisamos definir estes no bloco output. Então renomeei aquele collected, e novamente, para reports, um pouco verboso aqui, mas é realmente útil quando você vem ler um novo fluxo de trabalho, ver todas as diferentes saídas aqui, todos os diferentes canais listados lado a lado, e há maneiras de tornar isso menos verboso, o que tocaremos mais tarde.

Okay, vamos tentar executá-lo e rodar nosso fluxo de trabalho e ver o que acontece.

Esperançosamente agora deve executar basicamente o mesmo que fez antes. E vamos conseguir um novo arquivo de saída aqui em cima chamado replicate_two, report. E pronto. Abriu e diz que há três saudações no lote, o que é o que esperávamos, então está perfeito.

Se eu entrar no diretório work aqui apenas para provar que foi executado no código Nextflow, em vez do script bash, posso ir para cat work/ command.sh, e você verá aqui que está apenas fazendo echo desta string diretamente. There were three greetings in this batch, e então aquela variável foi interpolada pelo Nextflow. Foi calculada no bloco script antes de escrever o arquivo .command.sh. Então o cálculo da variável resultante é basicamente codificado nisto antes de ser executado em seu ambiente de computação neste caso.

E então você pode ver aquela separação entre o script. Bloco aqui e qualquer coisa acima dele. Espero que faça sentido.

## Conclusão e questionário

Okay, esse é o fim desta parte de Hello Nextflow. Então como antes, vá e confira o questionário. Faça na página web ou no CLI, passe por algumas das perguntas e apenas verifique se você entendeu algum do material que cobrimos. Veja se há algo lá que destaca algo que você pode não ter entendido. Não muitas perguntas. Agradável e fácil de fazer. Ou você pode fazer na página web aqui embaixo também.

E faça uma pequena pausa, uma pequena caminhada e volte e junte-se a nós na parte quatro de Hello Nextflow, onde falaremos sobre módulos. Muito obrigado.
