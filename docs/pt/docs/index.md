---
title: Início
description: Bem-vindo ao portal de treinamento da comunidade Nextflow!
hide:
  - toc
  - footer
---

# Treinamento Nextflow

<div class="grid cards" markdown>

-   :material-book-open-variant:{ .lg .middle } __Cursos de autoatendimento__

    ---

    **Bem-vindo ao portal de treinamento da comunidade Nextflow!**

    Os cursos de treinamento listados abaixo foram desenvolvidos para serem utilizados como um recurso de autoatendimento.
    Você pode realizá-los por conta própria a qualquer momento, seja no ambiente baseado na web que disponibilizamos via Github Codespaces ou no seu próprio ambiente.

    [Explorar os cursos :material-arrow-right:](#catalog-of-nextflow-training-courses){ .md-button .md-button--primary .mt-1 }

-   :material-information-outline:{ .lg .middle } __Informações adicionais__

    ---

    ??? warning "Compatibilidade de versões"

        <!-- Qualquer atualização neste conteúdo precisa ser copiada para a página de instalação local -->
        **A partir de janeiro de 2026, todos os nossos cursos de treinamento Nextflow exigem a versão 25.10.2 ou posterior do Nextflow, com sintaxe estrita ativada, salvo indicação em contrário.**

        Para mais informações sobre os requisitos de versão e sintaxe estrita, consulte o [guia de migração da documentação do Nextflow](https://nextflow.io/docs/latest/strict-syntax.html).

        Versões mais antigas do material de treinamento correspondentes à sintaxe anterior estão disponíveis através do seletor de versão na barra de menu desta página.

    ??? terminal "Opções de ambiente"

        Disponibilizamos um ambiente de treinamento baseado na web onde tudo o que você precisa para realizar o treinamento já está pré-instalado, disponível através do Github Codespaces (requer uma conta gratuita no GitHub).

        [![Abrir no GitHub Codespaces](https://github.com/codespaces/badge.svg)](https://codespaces.new/nextflow-io/training?quickstart=1&ref=master)

        Se isso não atender às suas necessidades, consulte as outras [opções de ambiente](./envsetup/index.md).

    ??? learning "Eventos de treinamento"

        Se você preferir realizar o treinamento Nextflow como parte de um evento estruturado, há muitas oportunidades para isso. Recomendamos verificar as seguintes opções:

        - **[Training Weeks]()** organizadas trimestralmente pela equipe da Comunidade
        - **[Seqera Events](https://seqera.io/events/)** incluem eventos de treinamento presenciais organizados pela Seqera (pesquise por 'Seqera Sessions' e 'Nextflow Summit')
        - **[Nextflow Ambassadors]()** organizam eventos para suas comunidades locais
        - **[nf-core events](https://nf-co.re/events)** incluem hackathons da comunidade

    ??? people "Informações para instrutores"

        Se você é um instrutor que conduz seus próprios treinamentos, é bem-vindo a utilizar nossos materiais diretamente do portal de treinamento, desde que atribua os devidos créditos. Consulte 'Créditos e contribuições' abaixo para mais detalhes.

        Além disso, adoraríamos ouvir de você sobre como poderíamos apoiar melhor seus esforços de treinamento! Entre em contato conosco em [community@seqera.io](mailto:community@seqera.io) ou no fórum da comunidade (consulte a página de [Ajuda](help.md)).

    ??? licensing "Licença de código aberto e política de contribuição"

        [![Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International (CC BY-NC-SA 4.0)](assets/img/cc_by-nc-sa.svg){ align=right }](https://creativecommons.org/licenses/by-nc-sa/4.0/)

        Este material de treinamento é desenvolvido e mantido pela [Seqera](https://seqera.io) e disponibilizado sob uma licença de código aberto ([CC BY-NC-SA](https://creativecommons.org/licenses/by-nc-sa/4.0/)) em benefício da comunidade. Se você deseja utilizar este material de uma forma que esteja fora do escopo da licença (observe as limitações sobre uso comercial e redistribuição), entre em contato conosco em [community@seqera.io](mailto:community@seqera.io) para discutir sua solicitação.

        Agradecemos melhorias, correções e relatórios de bugs da comunidade. Cada página possui um ícone :material-file-edit-outline: no canto superior direito que leva ao repositório de código, onde você pode reportar problemas ou propor alterações ao material de treinamento via pull request. Consulte o `README.md` no repositório para mais detalhes.

</div>

!!! note "Tradução assistida por IA"

    Esta tradução foi criada utilizando inteligência artificial e revisada por tradutores humanos.
    Agradecemos seu feedback e sugestões de melhorias.
    Consulte nosso [guia de tradução](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md) para mais informações.

## Catálogo de cursos de treinamento Nextflow

<div class="grid cards" markdown>

-   :material-walk:{ .lg .middle } __Trilha introdutória__

    ---

    ### :material-compass:{.nextflow-primary} Nextflow para Iniciantes {.mt-1}

    Cursos independentes de domínio destinados a quem é completamente novo no Nextflow. Cada curso consiste em uma série de módulos de treinamento projetados para ajudar os participantes a desenvolver suas habilidades progressivamente.

    ??? courses "**Hello Nextflow:** Aprenda a desenvolver seus próprios pipelines"

        Este curso aborda os componentes principais da linguagem Nextflow com detalhes suficientes para permitir o desenvolvimento de pipelines simples, mas totalmente funcionais, além de elementos-chave de design, desenvolvimento e práticas de configuração de pipelines.

        [Iniciar o treinamento Hello Nextflow :material-arrow-right:](hello_nextflow/index.md){ .md-button .md-button--secondary }

    ??? courses "**Nextflow Run:** Aprenda a executar pipelines existentes"

        Uma introdução concisa à execução e configuração de pipelines Nextflow, baseada no curso para desenvolvedores Hello Nextflow, mas com menos foco em código. Aborda execução, saídas, estrutura básica de código e configuração para diferentes ambientes de computação.

        [Iniciar o treinamento Nextflow Run :material-arrow-right:](nextflow_run/index.md){ .md-button .md-button--secondary }

    ---

    ### :material-microscope:{.nextflow-primary} Nextflow para Ciência {.mt-1}

    Aprenda a aplicar os conceitos e componentes apresentados no 'Hello Nextflow' a casos de uso científicos específicos.

    ??? courses "**Nextflow para Genômica** (chamada de variantes)"

        Para pesquisadores que desejam aprender a desenvolver seus próprios pipelines de genômica. O curso utiliza um caso de uso de chamada de variantes para demonstrar como desenvolver um pipeline de genômica simples, mas funcional.

        [Iniciar o treinamento Nextflow para Genômica :material-arrow-right:](nf4_science/genomics/index.md){ .md-button .md-button--secondary }

    ??? courses "**Nextflow para RNAseq** (RNAseq em massa)"

        Para pesquisadores que desejam aprender a desenvolver seus próprios pipelines de RNAseq. O curso utiliza um caso de uso de processamento de RNAseq em massa para demonstrar como desenvolver um pipeline de RNAseq simples, mas funcional.

        [Iniciar o treinamento Nextflow para RNAseq :material-arrow-right:](nf4_science/rnaseq/index.md){ .md-button .md-button--secondary }

    ??? courses "**Nextflow para Imagem** (omics espacial)"

        Para pesquisadores em imagem e omics espacial que desejam aprender a executar e personalizar pipelines de análise. O curso utiliza o pipeline nf-core/molkart para fornecer um pipeline biologicamente relevante e demonstrar como executar, configurar e gerenciar entradas para fluxos de trabalho de pipelines Nextflow.

        [Iniciar o treinamento Nextflow para Imagem :material-arrow-right:](nf4_science/imaging/index.md){ .md-button .md-button--secondary }

-   :material-run:{ .lg .middle } __Trilha avançada__

    ---

    ### :material-bridge:{.nextflow-primary} De Nextflow para nf-core {.mt-1}

    Aprenda a utilizar código e boas práticas do projeto comunitário [nf-core](https://nf-co.re/).

    Esses cursos ajudam você a ir dos fundamentos do Nextflow às boas práticas do nf-core.
    Entenda como e por que a comunidade nf-core desenvolve pipelines, e como você pode contribuir e reutilizar essas técnicas.

    ??? courses "**Hello nf-core:** Comece com o nf-core"

        Para desenvolvedores que desejam aprender a executar e desenvolver pipelines compatíveis com o [nf-core](https://nf-co.re/). O curso aborda a estrutura dos pipelines nf-core com detalhes suficientes para permitir o desenvolvimento de pipelines simples, mas totalmente funcionais, que seguem o template e as boas práticas de desenvolvimento do nf-core, além de utilizar módulos nf-core existentes.

        [Iniciar o treinamento Hello nf-core :material-arrow-right:](hello_nf-core/index.md){ .md-button .md-button--secondary }

    ---

    ### :material-rocket-launch:{.nextflow-primary} Treinamento Avançado em Nextflow {.mt-1}

    Aprenda conceitos e mecanismos avançados para desenvolver e implantar pipelines Nextflow para atender a casos de uso do mundo real.

    ??? courses "**Side Quests:** Mergulhos profundos em tópicos independentes"

        Mini-cursos independentes destinados a desenvolvedores Nextflow que desejam ampliar seu repertório e/ou aprofundar suas habilidades em tópicos específicos. São apresentados de forma linear, mas podem ser realizados em qualquer ordem (consulte as dependências na visão geral de cada mini-curso).

        [Explorar os Side Quests :material-arrow-right:](side_quests/index.md){ .md-button .md-button--secondary }

    ??? courses "**Training Collections:** Trilhas de aprendizado recomendadas pelos Side Quests"

        As Training Collections combinam múltiplos Side Quests para oferecer uma experiência de aprendizado abrangente em torno de um tema ou caso de uso específico.

        [Explorar as Training Collections :material-arrow-right:](training_collections/index.md){ .md-button .md-button--secondary }

</div>

!!! info "Procurando materiais de treinamento arquivados?"

    Materiais de treinamento mais antigos (Treinamento Fundamental, Treinamento Avançado e outros cursos experimentais) foram removidos do portal de treinamento por serem incompatíveis com a sintaxe estrita do Nextflow 3.0.
    Se você precisar acessar esses materiais, eles estão disponíveis no [histórico do git](https://github.com/nextflow-io/training) anterior a janeiro de 2026.

---

<div markdown class="homepage_logos">

![Seqera](assets/img/seqera_logo.png#only-light)

![Seqera](assets/img/seqera_logo_dark.png#only-dark)

</div>
