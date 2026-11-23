# Parte 1: Executar pipelines pela interface web

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Tradução assistida por IA - [saiba mais e sugira melhorias](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

Nesta parte do curso Scale with Seqera, você vai configurar o acesso à Seqera Platform e executar um pipeline em escala de produção pela interface web.

Certifique-se de que seu diretório de trabalho está definido como `seqera-scale/`, conforme instruído na página [Primeiros passos](./00_orientation.md).

---

## 1. Primeiros passos com Seqera

A Seqera oferece uma plataforma abrangente para executar, monitorar e gerenciar pipelines Nextflow.
Esta seção orienta você no processo de cadastro e familiarização antes de executar seu primeiro pipeline.

### 1.1. Crie uma conta gratuita

Acesse [cloud.seqera.io](https://cloud.seqera.io) e crie uma conta gratuita.
Você pode se cadastrar usando seu endereço de e-mail, GitHub ou credenciais do Google.

Uma conta gratuita oferece:

- **Workspace pessoal**: seu próprio espaço para adicionar pipelines, configurar ambientes de computação e gerenciar execuções
- **Acesso ao Community Showcase**: uma coleção curada de pipelines nf-core e da comunidade, com configurações pré-definidas e dados de exemplo para execução

Consulte a [documentação da Seqera](https://docs.seqera.io) para uma visão completa dos planos de conta e recursos disponíveis.

### 1.2. Explore o Community Showcase

Antes de executar seus próprios pipelines, reserve alguns minutos para explorar o Community Showcase.
Ele oferece uma prévia realista de como a Platform funciona com pipelines e dados reais.

1. Faça login em [cloud.seqera.io](https://cloud.seqera.io).
2. Na barra lateral esquerda, clique em **Showcase**.
3. Navegue pelos pipelines disponíveis — você vai reconhecer vários pipelines nf-core do curso Use nf-core.
4. Clique em um pipeline para ver sua configuração e as opções de execução.
5. Clique em **Runs** para explorar o histórico de execuções de exemplo, incluindo detalhes por tarefa e relatórios de execuções anteriores.

Esta é uma visualização somente leitura, mas mostra como a interface funciona antes de você executar qualquer coisa por conta própria.

### 1.3. Acesse um workspace com computação

Para executar pipelines, é necessário um workspace com um ambiente de computação configurado.

A Seqera oferece duas formas de fornecer computação:

- **Conecte sua própria infraestrutura**: AWS, Azure, Google Cloud e schedulers HPC (SLURM, LSF, PBS, entre outros).
  Consulte a [documentação de ambientes de computação](https://docs.seqera.io) para guias de configuração.
- **Seqera Compute**: um serviço gerenciado que fornece ambientes de computação pré-provisionados na AWS, mediante pagamento, sem necessidade de configurar uma conta na nuvem.
  Você pode ativá-lo diretamente nas configurações do seu workspace.

**Treinamento em grupo:**
Se você está participando de uma sessão de treinamento em grupo, pode ter sido adicionado a uma organização e workspace que já possui computação configurada.
Seu instrutor fornecerá o nome da organização, o nome do workspace e quaisquer outros detalhes necessários.

**Trabalhando de forma independente:**
Se você está seguindo este treinamento por conta própria, será necessário configurar um ambiente de computação no seu workspace pessoal usando uma das opções acima.
Créditos gratuitos para experimentar o Seqera Compute estão [disponíveis mediante solicitação](https://seqera.io/platform/compute/).

!!! note "Nota"

    O restante deste curso pressupõe que você tem acesso a um workspace com um ambiente de computação configurado.
    Se você estiver em uma sessão de treinamento em grupo, seu instrutor confirmará qual workspace e ambiente de computação utilizar.

### Conclusão

Você tem uma conta na Seqera, explorou o Community Showcase e consegue acessar um workspace com computação.

### O que vem a seguir?

Execute um pipeline RNA-seq em escala de produção pela interface web do Seqera Cloud.

---

## 2. Executar nf-core/rnaseq pela interface web

Conforme abordado no curso Use nf-core, o pipeline nf-core/rnaseq é um pipeline curado pela comunidade para análise de dados de sequenciamento de RNA em bulk.

Nesta seção, você vai adicionar o pipeline ao seu workspace, executar uma análise e monitorar sua execução.

### 2.1. Adicione o pipeline ao seu workspace

Convenientemente, o nf-core/rnaseq faz parte de uma coleção curada de pipelines que podem ser adicionados ao seu workspace com poucos cliques por meio do serviço Seqera Pipelines.

_Mostraremos como adicionar seus próprios pipelines mais adiante neste curso._

1. Acesse [**Seqera Pipelines**](https://seqera.io/pipelines) para navegar pela coleção da comunidade.
2. Pesquise por `rnaseq` e selecione **nf-core/rnaseq**.
3. Clique em **Launch Pipeline** ou role até o final da página até a seção **Launch Pipeline**.
4. Certifique-se de estar logado e selecione os valores apropriados nos menus suspensos **Organizations**, **Workspace** e **Compute Environment**.
   **Dica para grupos:** Se você estiver usando um workspace compartilhado, adicione um identificador único (como seu nome de usuário) ao nome do pipeline.
5. Clique em **Add pipeline to your Seqera account**

Uma caixa aparecerá com a mensagem: **Pipeline added: View Pipeline**.
Clicar no link levará você à entrada do pipeline no seu launchpad.

O pipeline agora está listado no painel **Launchpad** do seu workspace e está pronto para ser executado.

### 2.2. Execute o pipeline

Clique no botão **Launch** do pipeline, seja no painel **Launchpad** ou na página de detalhes do pipeline.
Isso abre a interface de configuração.

O pipeline já está configurado com o perfil `test`, portanto os dados de entrada, o diretório de saída e a referência do genoma já estão preenchidos.
Você pode ignorar os demais parâmetros e configurações avançadas por enquanto.

Clique no botão azul **Launch** para iniciar a execução de fato.

### 2.3. Monitore a execução

Após executar, você será direcionado ao painel **Runs** do seu pipeline.

A visualização da execução mostra:

- **Status**: estado atual da execução (submitted, running, succeeded, failed)
- **Command line**: o comando `nextflow run` exato que a Platform construiu e submeteu
- **Parameters**: todos os valores de parâmetros usados nesta execução
- **Tasks**: uma tabela com cada chamada de processo, com status, duração e uso de recursos

Clique em qualquer linha de tarefa para inspecionar seus detalhes de execução, incluindo:

- O script `.command.sh` que foi executado
- Logs de stdout e stderr
- Métricas de CPU, memória e I/O

A aba **Reports** exibirá um relatório MultiQC assim que a execução for concluída, agregando métricas de controle de qualidade de todas as amostras.

Isso levará algum tempo para ser executado, então vamos continuar por enquanto e voltar mais tarde para verificar as saídas e assim por diante.

### Conclusão

Você sabe como adicionar um pipeline a um workspace da Seqera, configurar e executar uma análise, e monitorar a execução em escala.

### O que vem a seguir?

Siga para a [Parte 2](./02_launch_from_cli.md), onde você aprenderá como fazer tudo isso pela linha de comando usando o CLI `tw`.

---

## Resumo

Nesta parte você aprendeu a:

- Criar uma conta na Seqera e explorar o Community Showcase
- Adicionar um pipeline do catálogo curado, executar uma análise em escala de produção e monitorar a execução
