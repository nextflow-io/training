# Instalação manual

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Tradução assistida por IA - [saiba mais e sugira melhorias](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

É possível instalar tudo o que você precisa para executar o treinamento no seu próprio ambiente local manualmente.

Aqui documentamos como fazer isso em sistemas compatíveis com POSIX padrão (assumindo uma máquina pessoal como um laptop).
Tenha em mente que alguns detalhes podem ser diferentes dependendo do seu sistema específico.

!!! tip "Dica"

    Antes de prosseguir, você considerou usar a [abordagem de Devcontainers](03_devcontainer.md)?
    Ela fornece todas as ferramentas e dependências necessárias sem exigir instalação manual.

## Requisitos gerais de software

Nextflow pode ser usado em qualquer sistema compatível com POSIX (Linux, macOS, Windows Subsystem for Linux, etc.) com Java instalado.
Nossos cursos de treinamento têm alguns requisitos adicionais.

No total, você precisará ter o seguinte software instalado:

- Bash ou shell equivalente
- [Java 11 (ou posterior, até 21)](https://www.oracle.com/technetwork/java/javase/downloads/index.html)
- [Git](https://git-scm.com/)
- [Docker](https://docs.docker.com/get-docker/)
- [Conda](https://conda.io/) 4.5 (ou posterior)
- [VSCode](https://code.visualstudio.com) com a [extensão Nextflow](https://www.nextflow.io/docs/latest/developer-env.html#devenv-nextflow)

A aplicação VSCode é tecnicamente opcional, mas recomendamos fortemente que você a use para trabalhar nos cursos, bem como para seu trabalho de desenvolvimento em Nextflow em geral.

O manual de documentação do Nextflow fornece instruções para instalar essas dependências em [Environment setup](https://www.nextflow.io/docs/latest/developer-env.html).

## Nextflow e ferramentas nf-core

Você precisará instalar o próprio Nextflow, mais as ferramentas nf-core, conforme detalhado nos artigos vinculados abaixo:

- [Instalação do Nextflow](https://www.nextflow.io/docs/latest/install.html)
- [Ferramentas nf-core](https://nf-co.re/docs/nf-core-tools/installation)

Recomendamos usar a opção de auto-instalação para Nextflow e a opção PyPI para ferramentas nf-core.

!!! warning "Aviso"

    <!-- Any update to this content needs to be copied to the home page -->
    **A partir de janeiro de 2026, todos os nossos cursos de treinamento em Nextflow requerem Nextflow versão 25.10.2 ou posterior, com sintaxe v2 estrita ativada, salvo indicação em contrário.**

    Para mais informações sobre requisitos de versão e sintaxe v2 estrita, consulte o guia [Versões do Nextflow](../info/nxf_versions.md).

    Versões mais antigas do material de treinamento correspondentes à sintaxe anterior estão disponíveis através do seletor de versão na barra de menu desta página web.

## Materiais de treinamento

A maneira mais fácil de baixar os materiais de treinamento é clonar o repositório inteiro usando este comando:

```bash
git clone https://github.com/nextflow-io/training.git
```

Cada curso tem seu próprio diretório.
Para trabalhar em um curso, abra uma janela de terminal (idealmente, de dentro da aplicação VSCode) e use `cd` para entrar no diretório relevante.

Você pode então seguir as instruções do curso fornecidas no site.
