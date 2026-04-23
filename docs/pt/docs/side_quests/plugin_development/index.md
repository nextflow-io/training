---
title: Desenvolvimento de Plugins
hide:
  - toc
---

# Desenvolvimento de Plugins

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Tradução assistida por IA - [saiba mais e sugira melhorias](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

O sistema de plugins do Nextflow permite que você estenda a linguagem com funções personalizadas, hooks de monitoramento, backends de execução e muito mais.
Os plugins permitem que a comunidade adicione funcionalidades ao Nextflow sem modificar seu núcleo, tornando-os ideais para compartilhar funcionalidades reutilizáveis entre pipelines.

Durante este treinamento, você aprenderá como usar plugins existentes e, opcionalmente, criar os seus próprios.

## Público & pré-requisitos

A Parte 1 aborda o uso de plugins existentes e é relevante para todos os usuários do Nextflow.
As Partes 2-6 abordam a criação dos seus próprios plugins e envolvem código Groovy e ferramentas de build.
Não é necessária experiência prévia com Java ou Groovy.

**Pré-requisitos**

- Uma conta no GitHub OU uma instalação local conforme descrito [aqui](../../envsetup/02_local).
- Ter concluído o curso [Hello Nextflow](../../hello_nextflow/index.md) ou equivalente.
- Java 21 ou superior (incluído no ambiente de treinamento; necessário apenas para as Partes 2-6).

**Diretório de trabalho:** `side-quests/plugin_development`

## Objetivos de aprendizado

Ao final deste treinamento, você será capaz de:

**Usando plugins (Parte 1):**

- Instalar e configurar plugins existentes nos seus fluxos de trabalho
- Importar e usar funções de plugins

**Desenvolvendo plugins (Partes 2-6):**

- Criar um novo projeto de plugin usando o gerador de projetos integrado do Nextflow
- Implementar funções personalizadas que podem ser chamadas a partir de fluxos de trabalho
- Compilar, testar e instalar seu plugin localmente
- Monitorar eventos do fluxo de trabalho (por exemplo, conclusão de tarefas, início/fim do pipeline) para logging personalizado ou notificações
- Adicionar opções de configuração para tornar os plugins personalizáveis
- Distribuir seu plugin

## Plano de aulas

#### Parte 1: Conceitos básicos de plugins

Use plugins existentes em um fluxo de trabalho Nextflow e configure seu comportamento.

#### Parte 2: Criar um projeto de plugin

Gere um novo projeto de plugin e examine sua estrutura.

#### Parte 3: Funções personalizadas

Implemente funções personalizadas, compile seu plugin e execute-o em um fluxo de trabalho.

#### Parte 4: Testes

Escreva e execute testes unitários usando o framework Spock.

#### Parte 5: Monitoramento de fluxo de trabalho

Responda a eventos como a conclusão de tarefas para construir um contador de tarefas.

#### Parte 6: Configuração & Distribuição

Leia configurações do `nextflow.config` para tornar seu plugin personalizável e, em seguida, aprenda como compartilhá-lo.

Pronto para começar o curso?

[Começar a aprender :material-arrow-right:](01_plugin_basics.md){ .md-button .md-button--primary }
