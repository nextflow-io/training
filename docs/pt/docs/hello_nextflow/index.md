---
title: Hello Nextflow
hide:
    - toc
page_type: index_page
index_type: course
additional_information:
    technical_requirements: true
    learning_objectives:
        - Iniciar e gerenciar a execução de fluxos de trabalho Nextflow
        - Encontrar e interpretar saídas (resultados) e arquivos de log gerados pelo Nextflow
        - Solucionar problemas básicos
        - Construir um fluxo de trabalho simples de múltiplas etapas a partir de componentes principais do Nextflow
        - Distinguir entre tipos essenciais de fábricas de canais e operadores e utilizá-los efetivamente em um fluxo de trabalho simples
        - Configurar a execução de pipelines para rodar em plataformas de computação comuns, incluindo HPC e nuvem
        - Aplicar melhores práticas para reprodutibilidade, portabilidade e reutilização de código que tornam pipelines FAIR, incluindo modularidade de código e contêineres de software
    audience_prerequisites:
        - "**Público:** Este curso é projetado para alunos que são completamente novos no Nextflow e desejam desenvolver seus próprios pipelines."
        - "**Habilidades:** Assume-se alguma familiaridade com a linha de comando, conceitos básicos de script e formatos de arquivo comuns."
        - "**Domínio:** Os exercícios são todos agnósticos de domínio, portanto nenhum conhecimento científico prévio é necessário."
    videos_playlist: https://youtube.com/playlist?list=PLPZ8WHdZGxmWKozQuzr27jyMGqp9kElVK&si=eF7cLR62goy-lc6n
---

# Hello Nextflow

**Hello Nextflow é uma introdução prática à construção de fluxos de trabalho de análise de dados reprodutíveis e escaláveis.**

Trabalhando através de exemplos práticos e exercícios guiados, você aprenderá os fundamentos do desenvolvimento de pipelines com Nextflow, incluindo como definir processos, conectá-los em pipelines, gerenciar arquivos e dependências de software, paralelizar a execução sem esforço e executar fluxos de trabalho em diferentes ambientes computacionais.

Você sairá com as habilidades e confiança para começar a desenvolver e executar seus próprios fluxos de trabalho com Nextflow.

<!-- additional_information -->

## Visão geral do curso

Este curso é projetado para ser prático, com exercícios orientados a objetivos estruturados para introduzir informações gradualmente.

Você desenvolverá um pipeline Nextflow simples que recebe algumas entradas de texto, executa algumas etapas de transformação e gera um único arquivo de texto contendo uma imagem ASCII de um personagem dizendo o texto transformado.

### Plano de aula

Para evitar sobrecarregá-lo com conceitos e código, dividimos isso em seis partes que focarão em aspectos específicos do desenvolvimento de pipelines com Nextflow.

| Capítulo do curso                                    | Resumo                                                                                                                | Duração estimada |
| ---------------------------------------------------- | --------------------------------------------------------------------------------------------------------------------- | ---------------- |
| [Parte 1: Hello World](./01_hello_world.md)          | Componentes básicos e princípios envolvidos na montagem e execução de um fluxo de trabalho Nextflow                  | 30 min           |
| [Parte 2: Hello Channels](./02_hello_channels.md)    | Usando canais e operadores para processar entradas e paralelizar a execução sem esforço                               | 45 min           |
| [Parte 3: Hello Workflow](./03_hello_workflow.md)    | Usando canais para encadear múltiplas etapas e lidar com a transferência de dados entre etapas                        | 60 min           |
| [Parte 4: Hello Modules](./04_hello_modules.md)      | Aplicando princípios de modularidade de código para aumentar a reutilização e diminuir a carga de manutenção         | 20 min           |
| [Parte 5: Hello Containers](./05_hello_containers.md)| Usando contêineres como mecanismo para gerenciar dependências de software e aumentar a reprodutibilidade             | 60 min           |
| [Parte 6: Hello Config](./06_hello_config.md)        | Personalizando o comportamento do pipeline e otimizando o uso em diferentes ambientes computacionais                  | 60 min           |

Ao final deste curso, você estará bem preparado para enfrentar os próximos passos em sua jornada para desenvolver fluxos de trabalho reprodutíveis para suas necessidades de computação científica.

Pronto para fazer o curso?

[Começar :material-arrow-right:](00_orientation.md){ .md-button .md-button--primary }

<!-- Clearfix for float -->
<div style="content: ''; clear: both; display: table;"></div>
