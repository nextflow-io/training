# Devcontainers Locais

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Tradução assistida por IA - [saiba mais e sugira melhorias](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

Se você tem uma instalação local do Docker ou está disposto a instalar uma, a maneira mais fácil de trabalhar localmente com esses materiais é usar o recurso de devcontainer do Visual Studio Code. Essa abordagem fornece todas as ferramentas e dependências necessárias sem exigir instalação manual.

## Requisitos

Para usar a configuração de devcontainer local, você precisará de:

- [Visual Studio Code](https://code.visualstudio.com/)
- Uma instalação local do Docker, por exemplo:
  - [Docker Desktop](https://docs.docker.com/get-docker/) (para Windows/macOS)
  - [Docker Engine](https://docs.docker.com/engine/install/) (para Linux)
  - [Colima](https://github.com/abiosoft/colima) (alternativa para macOS)
- [Docker Buildx](https://docs.docker.com/build/concepts/overview/#install-buildx) (incluído no Docker Desktop, mas pode precisar de instalação separada com outras configurações do Docker)
- [Extensão Dev Containers](https://marketplace.visualstudio.com/items?itemName=ms-vscode-remote.remote-containers) para VS Code

Sua instalação do Docker deve estar em execução antes de você tentar abrir o devcontainer.

Para verificar se o Docker buildx está disponível, execute:

```bash
docker buildx version
```

Se este comando falhar, você precisará instalar a extensão buildx antes de prosseguir.

## Instruções de Configuração

Siga estas etapas para configurar seu ambiente local usando devcontainers do VS Code:

### Instalar a extensão "Dev Containers" no VS Code

- Abra o VS Code
- Vá para Extensões (Ctrl+Shift+X ou Cmd+Shift+X no macOS)
- Pesquise por "Dev Containers"
- Clique em "Install"

![Instalando a extensão Dev Containers no VS Code](img/install_extension.png)

### Clone o repositório:

```bash
git clone https://github.com/nextflow-io/training.git
cd training
```

### Abra o repositório no VS Code:

- Inicie o VS Code
- Selecione **File -> Open Folder** no menu
- Navegue até e selecione a pasta do repositório de treinamento que você acabou de clonar
- Clique em **Open**

### Reabrir no Contêiner

Se solicitado pelo VS Code para "Reopen in Container", clique nele. Alternativamente:

- Pressione F1 (ou Ctrl+Shift+P / Cmd+Shift+P no macOS)
- Digite "Dev Containers: Reopen in Container"
- **Importante**: Quando solicitado a selecionar uma configuração, escolha a configuração de devcontainer **local-dev**

![Prompt Reopen in Container](img/reopen_prompt.png)

![Selecionando configuração local](img/select_local_config.png)

Aguarde a construção do contêiner. Isso pode levar alguns minutos na primeira vez, pois ele baixa e configura todos os componentes necessários.

Uma vez que o contêiner esteja construído e em execução, você terá um ambiente totalmente configurado com todas as ferramentas necessárias instaladas, incluindo:

- Java
- Nextflow
- Docker
- Git
- E todas as outras dependências necessárias para o treinamento

![VS Code com devcontainer em execução](img/running_container.png)

## Benefícios de Usar Devcontainers

Usar a abordagem de devcontainer oferece várias vantagens:

- **Consistência**: Garante um ambiente de desenvolvimento consistente em diferentes máquinas
- **Simplicidade**: Todas as dependências são pré-instaladas e configuradas
- **Isolamento**: O ambiente de desenvolvimento é isolado do seu sistema local
- **Reprodutibilidade**: Todos que usam o devcontainer obtêm a mesma configuração
- **Sem instalação manual**: Não é necessário instalar manualmente Java, Nextflow e outras ferramentas

## Verificando Seu Ambiente

Uma vez que seu devcontainer esteja em execução, você pode verificar se tudo está configurado corretamente executando:

```bash
nextflow info
```

Isso deve exibir a versão do Nextflow e informações de tempo de execução, confirmando que seu ambiente está configurado corretamente.

## Solução de Problemas

Se você encontrar problemas com a configuração do devcontainer:

1. Certifique-se de que sua instalação do Docker (Docker Desktop, Colima, Docker Engine, etc.) está em execução antes de abrir o devcontainer
2. Verifique se você selecionou a configuração **local-dev** quando solicitado
3. Verifique se o Docker buildx está instalado e funcionando executando `docker buildx version`
4. Se o contêiner falhar ao construir, tente reconstruí-lo executando o comando "Dev Containers: Rebuild Container"
5. Para problemas persistentes, consulte o [guia de solução de problemas do VS Code Dev Containers](https://code.visualstudio.com/docs/devcontainers/troubleshooting)
