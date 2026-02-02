# Resumo do curso

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Tradução assistida por IA - [saiba mais e sugira melhorias](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

Parabéns por concluir o curso de treinamento Nextflow Run! 🎉

<!-- placeholder for video -->

## Sua jornada

Você começou com um fluxo de trabalho muito básico e aprendeu a executá-lo, encontrar as saídas e gerenciar sua execução.
Então, você trabalhou através de versões cada vez mais complexas daquele fluxo de trabalho e aprendeu a reconhecer os conceitos e mecanismos essenciais que alimentam os pipelines Nextflow, incluindo canais e operadores, modularização de código e contêineres.
Finalmente, você aprendeu como personalizar a configuração de um pipeline para se adequar às suas preferências e sua infraestrutura computacional.

### O que você aprendeu

Você agora é capaz de gerenciar a execução do pipeline Hello, descrever como ele está estruturado e identificar as principais peças de código envolvidas.

- A forma final do fluxo de trabalho Hello recebe como entrada um arquivo CSV contendo saudações de texto.
- As quatro etapas são implementadas como processos Nextflow (`sayHello`, `convertToUpper`, `collectGreetings` e `cowpy`) armazenados em arquivos de módulo separados.
- Os resultados são publicados em um diretório chamado `results/`.
- A saída final do pipeline é um arquivo de texto simples contendo arte ASCII de um personagem dizendo as saudações em maiúsculas.

<figure class="excalidraw">
--8<-- "docs/en/docs/hello_nextflow/img/hello_pipeline_complete.svg"
</figure>

1. **`sayHello`:** Escreve cada saudação em seu próprio arquivo de saída (_ex._ "Hello-output.txt")
2. **`convertToUpper`:** Converte cada saudação para maiúsculas (_ex._ "HELLO")
3. **`collectGreetings`:** Coleta todas as saudações em maiúsculas em um único arquivo de lote
4. **`cowpy`:** Gera arte ASCII usando a ferramenta `cowpy`

A configuração do fluxo de trabalho suporta fornecer entradas e parâmetros de forma flexível e reproduzível.

### Habilidades adquiridas

Através deste curso prático, você aprendeu como:

- Iniciar um fluxo de trabalho Nextflow localmente
- Encontrar e interpretar saídas (resultados) e arquivos de log gerados pelo Nextflow
- Reconhecer os componentes principais do Nextflow que constituem um fluxo de trabalho simples de múltiplas etapas
- Descrever conceitos de próximos passos como operadores e fábricas de canal
- Configurar pipelines para diferentes ambientes de computação

Você agora está equipado com o conhecimento fundamental para começar a integrar pipelines Nextflow existentes em seu próprio trabalho.

## Próximos passos para construir suas habilidades

Aqui estão nossas principais sugestões do que fazer a seguir:

- Não apenas execute Nextflow, escreva! Torne-se um desenvolvedor Nextflow com [Hello Nextflow](../hello_nextflow/index.md)
- Aplique Nextflow a um caso de uso de análise científica com [Nextflow for Science](../nf4_science/index.md)
- Comece com nf-core com [Hello nf-core](../hello_nf-core/index.md)
- Aprenda técnicas de solução de problemas com a [Side Quest de Debugging](../side_quests/debugging.md)

Finalmente, recomendamos que você dê uma olhada na [**Seqera Platform**](https://seqera.io/), uma plataforma baseada em nuvem desenvolvida pelos criadores do Nextflow que torna ainda mais fácil lançar e gerenciar seus fluxos de trabalho, bem como gerenciar seus dados e executar análises interativamente em qualquer ambiente.

## Obtendo ajuda

Para recursos de ajuda e suporte da comunidade, veja a [página de Ajuda](../help.md).

## Pesquisa de feedback

Antes de seguir em frente, por favor tire um minuto para completar a pesquisa do curso! Seu feedback nos ajuda a melhorar nossos materiais de treinamento para todos.

[Fazer a pesquisa :material-arrow-right:](survey.md){ .md-button .md-button--primary }
