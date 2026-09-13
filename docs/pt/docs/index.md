---
title: Início
description: Bem-vindo ao portal de treinamento da comunidade Nextflow!
hide:
  - toc
  - footer
---

# Treinamento Nextflow

<div class="grid cards" markdown>

-   :material-book-open-variant:{ .lg .middle } __Cursos autônomos__

    ---

    **Bem-vindo ao portal de treinamento da comunidade Nextflow!**

    Faça os cursos abaixo no seu próprio ritmo, em nosso ambiente baseado na web ou no seu próprio.
    Cada curso é prático, com exercícios orientados a objetivos que você pode completar de forma independente.

    [Explorar os cursos :material-arrow-down:](#catalog-of-nextflow-training-courses){ .md-button .md-button--primary .mt-1 }

-   :material-account-group-outline:{ .lg .middle } __Eventos de Treinamento__

    ---

    **Procurando algo além dos cursos autônomos?**

    Encontre eventos de treinamento estruturados, orientações para conduzir seus próprios treinamentos e nossa licença de código aberto e política de contribuição.

    [Ver eventos de treinamento :material-arrow-right:](training_events.md){ .md-button .md-button--secondary .mt-1 }

</div>

!!! note "Tradução assistida por IA"

    Esta tradução foi criada utilizando inteligência artificial e revisada por tradutores humanos.
    Agradecemos seu feedback e sugestões de melhorias.
    Consulte nosso [guia de tradução](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md) para mais informações.

## Catalog of Nextflow training courses

<div class="grid cards" markdown>

-   :material-account:{ .lg .middle } __Para usuários__

    ---

    ### :material-play-circle:{.nextflow-primary} Executar pipelines {.mt-1}

    Aprenda a executar pipelines existentes sem escrever nenhum código.

    ??? courses "**Nextflow Run:** Execute pipelines com Nextflow"

        Uma introdução rápida à execução de pipelines Nextflow que não requer entendimento de código. Abrange o lançamento de pipelines, recuperação de saídas, uso de contêineres e configuração de execução em nível básico.

        [Ver o treinamento :material-arrow-right:](nextflow_run/index.md){ .md-button .md-button--secondary }

    ??? courses "**Use nf-core:** Encontre e execute pipelines curados pela comunidade"

        Uma introdução rápida para encontrar, executar e configurar pipelines do projeto comunitário nf-core, começando com um pipeline de demonstração mínimo e escalando até um pipeline de análise em escala de produção.

        [Ver o treinamento :material-arrow-right:](nfcore_use/index.md){ .md-button .md-button--secondary }

    ??? courses "**Scale with Seqera:** Lance e monitore pipelines em escala"

        Uma introdução prática ao lançamento e monitoramento de pipelines Nextflow com a Seqera Platform, tanto pela interface web quanto pela linha de comando.

        [Ver o treinamento :material-arrow-right:](seqera_scale/index.md){ .md-button .md-button--secondary }

    ---

    ### :material-tune:{.nextflow-primary} Gerenciar execução {.mt-1}

    Aprenda a gerenciar a execução de pipelines de forma eficaz.

    ??? courses "**Execution Config:** Configure pipelines como um profissional"

        Uma introdução prática à configuração da execução de pipelines Nextflow: adaptação a diferentes ambientes de computação, controle de alocações de recursos e tentativas de reexecução, e alternância entre perfis de configuração predefinidos.

        [Ver o treinamento :material-arrow-right:](execution_config/index.md){ .md-button .md-button--secondary }

    !!! info compact "Mais tópicos em breve"

        Ajuste de desempenho, execução em HPC/nuvem e muito mais estão planejados para esta seção.
        Vote no que abordar a seguir em nossa [pesquisa rápida de interesse](https://seqera.typeform.com/to/JCs91e8v).

-   :material-code-tags:{ .lg .middle } __Para desenvolvedores__

    ---

    ### :material-wrench:{.nextflow-primary} Escrever pipelines {.mt-1}

    Aprenda a desenvolver seus próprios pipelines Nextflow.

    ??? courses "**Hello Nextflow:** Desenvolva seus próprios pipelines do zero"

        Este curso abrange os componentes principais da linguagem Nextflow com detalhes suficientes para permitir o desenvolvimento de pipelines simples, mas totalmente funcionais, além de elementos-chave de design de pipelines, práticas de desenvolvimento e configuração.

        [Ver o treinamento :material-arrow-right:](hello_nextflow/index.md){ .md-button .md-button--secondary }

    ??? courses "**Build with nf-core:** Use as ferramentas e regras do nf-core"

        Para desenvolvedores Nextflow que desejam aprender a desenvolver pipelines compatíveis com o [nf-core](https://nf-co.re/).
        O curso abrange a estrutura dos pipelines nf-core com detalhes suficientes para permitir o desenvolvimento de pipelines simples, mas totalmente funcionais, que aproveitam o template nf-core e as melhores práticas de desenvolvimento, além do uso de módulos nf-core existentes.

        [Ver o treinamento :material-arrow-right:](nfcore_build/index.md){ .md-button .md-button--secondary }

    ??? catalog "**Side Quests:** Mergulhe em tópicos avançados de Nextflow"

        Uma coleção de mini-cursos independentes destinados a desenvolvedores Nextflow que desejam ampliar seu repertório e/ou aprofundar suas habilidades em tópicos específicos.
        São apresentados de forma linear, mas podem ser feitos em qualquer ordem (veja as dependências na visão geral de cada mini-curso).

        [Explorar os Side Quests :material-arrow-right:](side_quests/index.md){ .md-button .md-button--secondary }

    ---

    ### :material-microscope:{.nextflow-primary} Nextflow para Ciência {.mt-1}

    Aprenda a desenvolver pipelines Nextflow para aplicações científicas específicas.

    ??? courses "**Genomics:** Desenvolva um pipeline de chamada de variantes"

        Um curso para pesquisadores que desejam aprender a desenvolver seus próprios pipelines de genômica, usando um caso de uso de chamada de variantes para demonstrar padrões essenciais de desenvolvimento em Nextflow.

        [Ver o treinamento :material-arrow-right:](nf4_science/genomics/index.md){ .md-button .md-button--secondary }

    ??? courses "**RNAseq:** Desenvolva um pipeline de processamento de RNAseq em massa"

        Um curso para pesquisadores que desejam aprender a desenvolver seus próprios pipelines de RNAseq, usando um caso de uso de processamento de RNAseq em massa para demonstrar padrões essenciais de desenvolvimento em Nextflow.

        [Ver o treinamento :material-arrow-right:](nf4_science/rnaseq/index.md){ .md-button .md-button--secondary }

    ??? courses "**Bioimaging:** Execute e configure pipelines de imagem"

        Um curso para pesquisadores que desejam aprender a executar e configurar pipelines de bioimagem, usando o nf-core/molkart para demonstrar padrões essenciais de uso do Nextflow.

        [Ver o treinamento :material-arrow-right:](nf4_science/imaging/index.md){ .md-button .md-button--secondary }

</div>

## Configuração e Ajuda

<div class="grid cards mb-4" markdown>

-   :material-cog-outline:{ .lg .middle } __Ambiente de Treinamento__

    ---

    Opções para configurar seu ambiente para os treinamentos de Nextflow.

    [Ver os ambientes de treinamento :material-arrow-right:](envsetup/index.md){ .md-button .md-button--secondary }

-   :material-tag-outline:{ .lg .middle } __Versões do Nextflow__

    ---

    Entendendo e gerenciando a evolução das versões de sintaxe do Nextflow.

    [Verificar requisitos de versão :material-arrow-right:](info/nxf_versions.md){ .md-button .md-button--secondary }

-   :material-file-code-outline:{ .lg .middle } __O pipeline Hello__

    ---

    Resumo do que o pipeline Hello faz e como ele é estruturado.

    [Ler o resumo :material-arrow-right:](info/hello_pipeline.md){ .md-button .md-button--secondary }

-   :material-lifebuoy:{ .lg .middle } __Obtendo ajuda__

    ---

    Recursos úteis quando você tiver algum problema com o treinamento de Nextflow.

    [Encontrar ajuda :material-arrow-right:](help.md){ .md-button .md-button--secondary }

</div>
