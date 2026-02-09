# Resumo do curso

Parabéns por concluir o curso de treinamento Nextflow Run! 🎉

<!-- placeholder for video -->

## Sua jornada

Você começou com um fluxo de trabalho muito básico e aprendeu a executá-lo, encontrar as saídas e gerenciar sua execução.
Em seguida, você trabalhou com versões cada vez mais complexas desse fluxo de trabalho e aprendeu a reconhecer os conceitos e mecanismos essenciais que impulsionam os pipelines Nextflow, incluindo canais e operadores, modularização de código e contêineres.
Por fim, você aprendeu a personalizar a configuração de um pipeline para se adequar às suas preferências e à sua infraestrutura computacional.

### O que você aprendeu

Agora você é capaz de gerenciar a execução do pipeline Hello, descrever como ele está estruturado e identificar as principais partes do código envolvidas.

- A forma final do fluxo de trabalho Hello recebe como entrada um arquivo CSV contendo saudações em texto.
- As quatro etapas são implementadas como processos Nextflow (`sayHello`, `convertToUpper`, `collectGreetings` e `cowpy`) armazenados em arquivos de módulo separados.
- Os resultados são publicados em um diretório chamado `results/`.
- A saída final do pipeline é um arquivo de texto simples contendo arte ASCII de um personagem dizendo as saudações em maiúsculas.

<figure class="excalidraw">
--8<-- "docs/en/docs/hello_nextflow/img/hello_pipeline_complete.svg"
</figure>

1. **`sayHello`:** Escreve cada saudação em seu próprio arquivo de saída (_ex._ "Hello-output.txt")
2. **`convertToUpper`:** Converte cada saudação para maiúsculas (_ex._ "HELLO")
3. **`collectGreetings`:** Coleta todas as saudações em maiúsculas em um único arquivo em lote
4. **`cowpy`:** Gera arte ASCII usando a ferramenta `cowpy`

A configuração do fluxo de trabalho suporta o fornecimento de entradas e parâmetros de forma flexível e reproduzível.

### Habilidades adquiridas

Através deste curso prático, você aprendeu como:

- Executar um fluxo de trabalho Nextflow localmente
- Encontrar e interpretar saídas (resultados) e arquivos de log gerados pelo Nextflow
- Reconhecer os componentes principais do Nextflow que constituem um fluxo de trabalho simples de múltiplas etapas
- Descrever conceitos de próximo nível, como operadores e fábricas de canais
- Configurar pipelines para diferentes ambientes computacionais

Agora você está equipado com o conhecimento fundamental para começar a integrar pipelines Nextflow existentes em seu próprio trabalho.

## Próximos passos para desenvolver suas habilidades

Aqui estão nossas principais sugestões para o que fazer a seguir:

- Não apenas execute Nextflow, escreva-o! Torne-se um desenvolvedor Nextflow com [Hello Nextflow](../hello_nextflow/index.md)
- Aplique Nextflow a um caso de uso de análise científica com [Nextflow for Science](../nf4_science/index.md)
- Comece com nf-core com [Hello nf-core](../hello_nf-core/index.md)
- Aprenda técnicas de solução de problemas com a [Missão Secundária de Depuração](../side_quests/debugging.md)

Por fim, recomendamos que você dê uma olhada na [**Seqera Platform**](https://seqera.io/), uma plataforma baseada em nuvem desenvolvida pelos criadores do Nextflow que torna ainda mais fácil executar e gerenciar seus fluxos de trabalho, além de gerenciar seus dados e executar análises interativamente em qualquer ambiente.

## Obtendo ajuda

Para recursos de ajuda e suporte da comunidade, consulte a [página de Ajuda](../help.md).

## Pesquisa de feedback

Antes de prosseguir, reserve um minuto para completar a pesquisa do curso! Seu feedback nos ajuda a melhorar nossos materiais de treinamento para todos.

[Responder à pesquisa :material-arrow-right:](survey.md){ .md-button .md-button--primary }
