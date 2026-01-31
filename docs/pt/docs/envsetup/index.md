---
title: Opções de ambiente
description: Opções para configurar seu ambiente para os treinamentos Nextflow
hide:
  - toc
  - footer
---

# Opções de ambiente

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Tradução assistida por IA - [saiba mais e sugira melhorias](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

Nosso objetivo é fornecer um ambiente consistente e completamente testado que permita aos alunos se concentrarem em aprender Nextflow sem ter que gastar tempo e esforço gerenciando software.
Para isso, desenvolvemos um ambiente em contêiner que contém todo o software necessário, arquivos de código e dados de exemplo para trabalhar em todos os nossos cursos.

Este ambiente em contêiner pode ser executado diretamente no Github Codespaces ou localmente no VS Code com a extensão Devcontainers.

<div class="grid cards" markdown>

-   :material-cloud-outline:{ .lg .middle } __Github Codespaces__

    ---

    GitHub Codespaces é um serviço baseado na web que nos permite fornecer um ambiente pré-construído para treinamento, com todas as ferramentas e dados incluídos, apoiado por máquinas virtuais na nuvem. É acessível gratuitamente para qualquer pessoa com uma conta Github.

    [Usar Github Codespaces:material-arrow-right:](01_setup.md){ .md-button .md-button--primary .mt-1 }

-   :material-laptop:{ .lg .middle } __Devcontainers Locais__

    ---

    VS Code com Devcontainers fornece um ambiente de desenvolvimento em contêiner executado localmente com todas as ferramentas de treinamento pré-configuradas. Oferece o mesmo ambiente pré-construído do Codespaces, mas executando inteiramente em seu hardware local.

    [Usar Devcontainers localmente :material-arrow-right:](03_devcontainer.md){ .md-button .md-button--primary .mt-1 }

</div>

## Instruções para instalação manual

Se nenhuma das opções acima atender às suas necessidades, você pode replicar este ambiente em seu próprio sistema local instalando as dependências de software manualmente e clonando o repositório de treinamento.

[Instalação manual :material-arrow-right:](02_local.md){ .md-button .md-button--primary .mt-1 }

---

!!! info "Descontinuação do Gitpod"

    O Treinamento Nextflow costumava usar o [Gitpod](https://gitpod.io) até fevereiro de 2025.
    No entanto, os criadores do Gitpod decidiram encerrar a funcionalidade gratuita em favor do sistema [Gitpod Flex](https://www.gitpod.io/blog/introducing-gitpod-flex).
    Por esse motivo, mudamos para o uso do GitHub Codespaces, que também oferece um ambiente de desenvolvimento com um clique, sem configuração prévia.

    Dependendo de quando você se cadastrou no Gitpod e quando exatamente eles encerrarem o serviço, você ainda pode conseguir iniciar o treinamento em seu antigo IDE na nuvem, embora não possamos garantir acesso confiável daqui para frente:
    [Abrir no Gitpod](https://gitpod.io/#https://github.com/nextflow-io/training).
