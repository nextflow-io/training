---
title: Início
description: Bem-vindo ao portal de treinamento da comunidade Nextflow!
hide:
  - toc
  - footer
---

# Treinamento Nextflow

<div class="grid cards" markdown>

-   :material-book-open-variant:{ .lg .middle } __Cursos de autoaprendizagem__

    ---

    **Bem-vindo ao portal de treinamento da comunidade Nextflow!**

    Os cursos de treinamento listados abaixo foram projetados para serem usados como um recurso de autoaprendizagem.
    Você pode trabalhar neles por conta própria a qualquer momento, seja no ambiente baseado na web que fornecemos via Github Codespaces ou em seu próprio ambiente.

    [Explore os cursos :material-arrow-right:](#catalog-of-nextflow-training-courses){ .md-button .md-button--primary .mt-1 }

-   :material-information-outline:{ .lg .middle } __Informações adicionais__

    ---

    ??? warning "Compatibilidade de versão"

        <!-- Qualquer atualização neste conteúdo precisa ser copiada para a página de instalação local -->
        **A partir de janeiro de 2026, todos os nossos cursos de treinamento Nextflow requerem a versão 25.10.2 ou posterior do Nextflow, com sintaxe estrita ativada, salvo indicação em contrário.**

        Para mais informações sobre requisitos de versão e sintaxe estrita, consulte o [guia de migração da documentação do Nextflow](https://nextflow.io/docs/latest/strict-syntax.html).

        Versões mais antigas do material de treinamento correspondentes à sintaxe anterior estão disponíveis através do seletor de versão na barra de menu desta página web.

    ??? terminal "Opções de ambiente"

        Fornecemos um ambiente de treinamento baseado na web onde tudo o que você precisa para fazer o treinamento está pré-instalado, disponível através do Github Codespaces (requer uma conta gratuita do GitHub).

        [![Abrir no GitHub Codespaces](https://github.com/codespaces/badge.svg)](https://codespaces.new/nextflow-io/training?quickstart=1&ref=master)

        Se isso não atender às suas necessidades, consulte as outras [Opções de ambiente](./envsetup/index.md).

    ??? learning "Eventos de treinamento"

        Se você preferir fazer o treinamento Nextflow como parte de um evento estruturado, há muitas oportunidades para isso. Recomendamos conferir as seguintes opções:

        - **[Training Weeks]()** organizadas trimestralmente pela equipe da Comunidade
        - **[Seqera Events](https://seqera.io/events/)** incluem eventos de treinamento presenciais organizados pela Seqera (procure por 'Seqera Sessions' e 'Nextflow Summit')
        - **[Nextflow Ambassadors]()** organizam eventos para sua comunidade local
        - **[nf-core events](https://nf-co.re/events)** incluem hackathons da comunidade

    ??? people "Informações para instrutores"

        Se você é um instrutor realizando seus próprios treinamentos, você é bem-vindo para usar nossos materiais diretamente do portal de treinamento, desde que atribua o crédito apropriado. Veja 'Créditos e contribuições' abaixo para detalhes.

        Além disso, adoraríamos ouvir de você sobre como poderíamos apoiar melhor seus esforços de treinamento! Entre em contato conosco em [community@seqera.io](mailto:community@seqera.io) ou no fórum da comunidade (veja a página de [Ajuda](help.md)).

    ??? licensing "Licença de código aberto e política de contribuição"

        [![Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International (CC BY-NC-SA 4.0)](assets/img/cc_by-nc-sa.svg){ align=right }](https://creativecommons.org/licenses/by-nc-sa/4.0/)

        Este material de treinamento é desenvolvido e mantido pela [Seqera](https://seqera.io) e lançado sob uma licença de código aberto ([CC BY-NC-SA](https://creativecommons.org/licenses/by-nc-sa/4.0/)) para o benefício da comunidade. Se você deseja usar este material de uma forma que esteja fora do escopo da licença (observe as limitações sobre uso comercial e redistribuição), entre em contato conosco em [community@seqera.io](mailto:community@seqera.io) para discutir sua solicitação.

        Recebemos melhorias, correções e relatórios de bugs da comunidade. Cada página tem um ícone :material-file-edit-outline: no canto superior direito da página que leva ao repositório de código, onde você pode relatar problemas ou propor mudanças no material de treinamento através de um pull request. Veja o `README.md` no repositório para mais detalhes.

</div>

!!! note "Tradução assistida por IA"

    Esta tradução foi criada utilizando inteligência artificial e revisada por tradutores humanos.
    Agradecemos seu feedback e sugestões de melhorias.
    Consulte nosso [guia de tradução](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md) para mais informações.

## Catálogo de cursos de treinamento Nextflow

<div class="grid cards" markdown>

-   :material-walk:{ .lg .middle } __Trilha introdutória__

    ---

    ### :material-compass:{.nextflow-primary} Nextflow for Newcomers {.mt-1}

    Cursos independentes de domínio destinados àqueles que são completamente novos no Nextflow. Cada curso consiste em uma série de módulos de treinamento projetados para ajudar os alunos a desenvolver suas habilidades progressivamente.

    ??? courses "**Hello Nextflow:** Aprenda a desenvolver seus próprios pipelines"

        Este curso cobre os componentes principais da linguagem Nextflow em detalhes suficientes para permitir o desenvolvimento de pipelines simples, mas totalmente funcionais, além de elementos-chave de práticas de design, desenvolvimento e configuração de pipelines.

        [Iniciar o treinamento Hello Nextflow :material-arrow-right:](hello_nextflow/index.md){ .md-button .md-button--secondary }

    ??? courses "**Nextflow Run:** Aprenda a executar pipelines existentes"

        Uma introdução concisa à execução e configuração de pipelines Nextflow, baseada no curso de desenvolvedor Hello Nextflow, mas com menos foco em código. Cobre execução, saídas, estrutura básica de código e configuração para diferentes ambientes de computação.

        [Iniciar o treinamento Nextflow Run :material-arrow-right:](nextflow_run/index.md){ .md-button .md-button--secondary }

    ---

    ### :material-microscope:{.nextflow-primary} Nextflow for Science {.mt-1}

    Aprenda a aplicar os conceitos e componentes apresentados em 'Hello Nextflow' a casos de uso científicos específicos.

    ??? courses "**Nextflow for Genomics** (chamada de variantes)"

        Para pesquisadores que desejam aprender como desenvolver seus próprios pipelines de genômica. O curso usa um caso de uso de chamada de variantes para demonstrar como desenvolver um pipeline de genômica simples, mas funcional.

        [Iniciar o treinamento Nextflow for Genomics :material-arrow-right:](nf4_science/genomics/){ .md-button .md-button--secondary }

    ??? courses "**Nextflow for RNAseq** (RNAseq em massa)"

        Para pesquisadores que desejam aprender como desenvolver seus próprios pipelines de RNAseq. O curso usa um caso de uso de processamento de RNAseq em massa para demonstrar como desenvolver um pipeline de RNAseq simples, mas funcional.

        [Iniciar o treinamento Nextflow for RNAseq :material-arrow-right:](nf4_science/rnaseq/){ .md-button .md-button--secondary }

    ??? courses "**Nextflow for Imaging** (ômica espacial)"

        Para pesquisadores em imagem e ômica espacial que desejam aprender como executar e personalizar pipelines de análise. O curso usa o pipeline nf-core/molkart para fornecer um pipeline biologicamente relevante que demonstra como executar, configurar e gerenciar entradas para fluxos de trabalho de pipelines Nextflow.

        [Iniciar o treinamento Nextflow for Imaging :material-arrow-right:](nf4_science/imaging/){ .md-button .md-button--secondary }

-   :material-run:{ .lg .middle } __Trilha avançada__

    ---

    ### :material-bridge:{.nextflow-primary} From Nextflow to nf-core {.mt-1}

    Aprenda a utilizar código e melhores práticas do projeto comunitário [nf-core](https://nf-co.re/).

    Esses cursos ajudam você a ir dos fundamentos do Nextflow às melhores práticas do nf-core.
    Entenda como e por que a comunidade nf-core constrói pipelines, e como você pode contribuir e reutilizar essas técnicas.

    ??? courses "**Hello nf-core:** Comece com o nf-core"

        Para desenvolvedores que desejam aprender a executar e desenvolver pipelines compatíveis com [nf-core](https://nf-co.re/). O curso cobre a estrutura dos pipelines nf-core em detalhes suficientes para permitir o desenvolvimento de pipelines simples, mas totalmente funcionais, que seguem o template nf-core e as melhores práticas de desenvolvimento, além de usar módulos nf-core existentes.

        [Iniciar o treinamento Hello nf-core :material-arrow-right:](hello_nf-core/index.md){ .md-button .md-button--secondary }

    ---

    ### :material-rocket-launch:{.nextflow-primary} Advanced Nextflow Training {.mt-1}

    Aprenda conceitos e mecanismos avançados para desenvolver e implantar pipelines Nextflow para abordar casos de uso do mundo real.

    ??? courses "**Side Quests:** Mergulhos profundos em tópicos independentes"

        Mini-cursos independentes destinados a desenvolvedores Nextflow que desejam ampliar seu alcance e/ou aprofundar suas habilidades em tópicos específicos. Eles são apresentados linearmente, mas podem ser feitos em qualquer ordem (veja as dependências na visão geral de cada mini-curso).

        [Navegar pelos Side Quests :material-arrow-right:](side_quests/){ .md-button .md-button--secondary }

    ??? courses "**Training Collections:** Caminhos de aprendizagem recomendados através dos Side Quests"

        As Training Collections combinam múltiplos Side Quests para fornecer uma experiência de aprendizagem abrangente em torno de um tema ou caso de uso específico.

        [Navegar pelas Training Collections :material-arrow-right:](training_collections/){ .md-button .md-button--secondary }

</div>

!!! info "Procurando materiais de treinamento arquivados?"

    Materiais de treinamento mais antigos (Fundamentals Training, Advanced Training e outros cursos experimentais) foram removidos do portal de treinamento por serem incompatíveis com a sintaxe estrita do Nextflow 3.0.
    Se você precisar acessar esses materiais, eles estão disponíveis no [histórico do git](https://github.com/nextflow-io/training) antes de janeiro de 2026.

---

<div markdown class="homepage_logos">

![Seqera](assets/img/seqera_logo.png#only-light)

![Seqera](assets/img/seqera_logo_dark.png#only-dark)

</div>
