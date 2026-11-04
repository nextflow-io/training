# Parte 2: Executar pipelines pela linha de comando

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Tradução assistida por IA - [saiba mais e sugira melhorias](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

Na [Parte 1](./01_run_with_seqera.md), você executou o nf-core/rnaseq pela interface web do Seqera.
Agora faremos o mesmo pela linha de comando usando o `tw` CLI, e adicionaremos um novo pipeline ao seu workspace.

---

## 1. Executar pipelines pela linha de comando

Na visualização de execução, clique na aba **Command line**.
Você verá o comando `nextflow run` exato que a Platform construiu e enviou em seu nome — o mesmo tipo de comando que você executou manualmente no curso Use nf-core.

A Platform não substitui o Nextflow; ela o orquestra.
Tudo o que você pode fazer pela interface web, também pode fazer em um terminal usando o `tw` CLI, a ferramenta de linha de comando para interagir com a API da Platform.
Isso é útil para automatizar execuções a partir de scripts ou pipelines de CI/CD.

Faremos isso agora a partir do mesmo codespace que você usou nos cursos anteriores.

### 1.1. Instalar o tw CLI

Execute os seguintes comandos no terminal do seu Codespace para baixar e instalar o binário `tw`:

```bash
curl -fsSL https://github.com/seqeralabs/tower-cli/releases/latest/download/tw-linux-x86_64 -o tw
chmod +x tw
sudo mv tw /usr/local/bin/
```

Verifique a instalação:

```bash
tw --version
```

??? success "Saída do comando"

    ```console
    Tower CLI version 0.40.0 (build 26db579)
    ```

O `tw` CLI está instalado e pronto para ser configurado.

### 1.2. Obter um token de acesso

O `tw` CLI se autentica no Seqera usando um token de acesso pessoal.

1. Na interface web do Seqera, clique no seu avatar no canto superior direito e selecione **Your tokens**.
2. Clique em **Add token**, dê um nome a ele (por exemplo, `training`) e clique em **Add**.
3. Copie o valor do token — ele será exibido apenas uma vez.
   Se você não salvá-lo em algum lugar imediatamente, precisará gerar outro.

### 1.3. Configurar o CLI

Por conveniência, vamos configurar um arquivo de configuração contendo o
token de acesso que você acabou de gerar e o identificador do workspace.

Abra o arquivo `.seqera_config` neste diretório no editor e defina as duas variáveis:

- **`TOWER_ACCESS_TOKEN`**: o token que você gerou na seção 1.2
- **`TOWER_WORKSPACE_ID`**: o ID numérico do seu workspace (a coluna `ID` em `tw workspaces list`, que você executa na seção 1.4)

Depois de preencher os valores, carregue a configuração:

```bash
source .seqera_config
```

Verifique a conexão:

```bash
tw info
```

??? success "Saída do comando"

    ```console
        Details
    -------------------------+-----------------------------
     Tower API endpoint      | https://api.cloud.seqera.io
     Tower API version       | 1.150.0
     Tower version           | 26.1.0-cycle54
     CLI version             | 0.30.0 (fde9dec)
     CLI minimum API version | 1.148.0
     Authenticated user      | <your-name>

    System health status
    ---------------------------------------+----
     Remote API server connection check    | OK
     Tower API version check               | OK
     Authentication API credential's token | OK
    ```

O `tw` CLI agora está autenticado e conectado à sua conta Seqera.
Execute `source .seqera_config` no início de cada sessão do Codespace para recarregar a configuração.

!!! tip "Dica"

    Se o seu workspace não tiver um ambiente de computação primário definido, você pode adicionar `export TOWER_COMPUTE_ENV=<compute-env-name>` ao seu arquivo de configuração para definir um padrão.
    Qualquer valor de configuração pode ser substituído na linha de comando passando a flag explicitamente (por exemplo, `--compute-env other-env`).
    Consulte a [referência do tw CLI](https://docs.seqera.io/platform/latest/cli/reference) para a lista completa de opções e variáveis de ambiente.

### 1.4. Explorar seu workspace pelo CLI

Liste os workspaces aos quais você tem acesso:

```bash
tw workspaces list
```

??? success "Saída do comando"

    ```console
    Available workspaces:
    ID              | Name            | Full Name                | Visibility
    --------------- | --------------- | ------------------------ | ----------
    <workspace-id>  | my-workspace    | my-org/my-workspace      | PRIVATE
    ```

Visualize as execuções no seu workspace, incluindo a execução do nf-core/rnaseq que você acabou de iniciar:

```bash
tw runs list
```

??? success "Saída do comando"

    ```console
    Pipeline runs for my-org/my-workspace workspace:
    ID        | Status   | Name          | Pipeline               | Run name
    --------- | -------- | ------------- | ---------------------- | --------
    <run-id>  | RUNNING  | ...           | nf-core/rnaseq         | happy_curie
    ```

A mesma execução que você está monitorando na interface web está visível aqui.

!!! note "Nota"

    Como `TOWER_WORKSPACE_ID` está definido em `.seqera_config`, você pode omitir `--workspace` de todos os comandos `tw`.
    Sem a configuração, você o passaria explicitamente:

    ```bash
    tw runs list --workspace <org>/<workspace>
    ```

Tudo o que está visível na interface web é acessível pelo CLI.

### 1.5. Executar nf-core/rnaseq pelo CLI

O pipeline que você adicionou ao seu workspace na [Parte 1](./01_run_with_seqera.md) está disponível pelo nome no CLI.
Execute-o com o perfil `test`:

```bash
tw launch nf-core-rnaseq --profile test
```

??? success "Saída do comando"

    ```console
    Launching pipeline nf-core-rnaseq
    Run name: focused_einstein
    https://cloud.seqera.io/orgs/my-org/workspaces/my-workspace/watch/<run-id>
    ```

Abra o link no seu navegador e confirme que a execução aparece no painel **Runs**.

Assim que você puder vê-la em execução, terá confirmado que o CLI e a interface web são duas visualizações do mesmo workspace.

!!! note "Nota"

    Você também pode passar uma URL completa do GitHub diretamente para `tw launch` sem adicionar o pipeline a um workspace primeiro.
    No entanto, adicionar o pipeline explicitamente antes de executá-lo é geralmente melhor: isso salva a configuração do pipeline para execuções futuras, torna-o disponível pelo nome e o torna visível para todos os membros do workspace no Launchpad.

    É possível adicionar um pipeline a um workspace diretamente pela linha de comando usando `tw`.
    A próxima seção mostra como fazer isso com o pipeline nf-core/demo.

### Conclusão

Você sabe como autenticar o `tw` CLI, inspecionar seu workspace e executar um pipeline salvo pelo terminal.

### O que vem a seguir?

Adicionar um novo pipeline ao seu workspace pela linha de comando e executá-lo.

---

## 2. Adicionar um novo pipeline e executá-lo

Qualquer pipeline Nextflow no GitHub pode ser adicionado ao seu workspace com `tw pipelines add`, desde que tenha um ponto de entrada `main.nf` e um `nextflow.config` na raiz.
O nf-core/demo é um bom exemplo para praticar: você já o executou no curso Use nf-core, então você sabe o que ele faz e o que esperar.

### 2.1. Adicionar nf-core/demo ao seu workspace

Execute o seguinte comando para registrar o pipeline no seu workspace:

```bash
tw pipelines add \
    --name nf-core-demo \
    https://github.com/nf-core/demo
```

??? success "Saída do comando"

    ```console
    New pipeline 'nf-core-demo' added at my-org/my-workspace workspace
    ```

O pipeline agora está registrado e aparecerá no Launchpad.

### 2.2. Verificar se ele aparece no Launchpad

Liste os pipelines no seu workspace para confirmar que foi adicionado:

```bash
tw pipelines list
```

??? success "Saída do comando"

    ```console
    Pipelines at my-org/my-workspace:
    ID   | Name            | Repository
    ---- | --------------- | -----------------------------------------
    ...  | nf-core-demo    | https://github.com/nf-core/demo
    ...  | nf-core-rnaseq  | ...
    ```

Abra seu workspace no navegador e clique em **Launchpad** para confirmar que o nf-core/demo agora aparece ao lado do nf-core/rnaseq.

!!! tip "Dica"

    Você também pode adicionar pipelines pela interface web: na barra lateral esquerda, clique em **Launchpad**, depois em **Add pipeline**, e preencha o formulário adequadamente.

Clique no botão **Launch** na entrada do nf-core/demo para abrir seu formulário de execução.
Você verá que os parâmetros `input` e `outdir` estão destacados em vermelho — são campos obrigatórios sem valores padrão, porque `tw pipelines add` registra apenas o código-fonte do pipeline sem pré-configurar nenhum parâmetro.
As próximas duas seções mostram como fornecer esses valores: primeiro pelo formulário web, depois pela linha de comando.

### 2.3. Executar nf-core/demo pela interface web

Com o formulário de execução aberto, preencha os dois parâmetros obrigatórios.

Para `input`, insira a URL do samplesheet de teste do perfil de teste do nf-core/demo.
Você pode encontrá-la em `conf/test.config` dentro do repositório do pipeline, que você examinou no curso Use nf-core:

```
https://raw.githubusercontent.com/nf-core/test-datasets/viralrecon/samplesheet/samplesheet_test_illumina_amplicon.csv
```

Para `outdir`, insira um caminho de armazenamento em nuvem onde o pipeline possa gravar seus resultados.
Use o bucket configurado para o seu workspace, com um subdiretório para manter as execuções organizadas:

```
s3://my-bucket/demo-results
```

Depois de preencher ambos os campos, clique no botão azul **Launch**.

A execução aparece no painel **Runs** e deve ser concluída em alguns minutos com o conjunto de dados de teste.
Clique na execução para explorar a tabela de tarefas e quaisquer relatórios de execução.

### 2.4. Executar nf-core/demo pelo CLI

Ao contrário do `nextflow run`, o comando `tw launch` não aceita flags de parâmetros individuais como `--input` ou `--outdir`.
Os parâmetros devem ser fornecidos por meio de um arquivo no formato YAML ou JSON, passado com `--params-file`.
Isso incentiva a reprodutibilidade: um arquivo de parâmetros salvo documenta exatamente quais valores foram usados em uma execução, facilitando a repetição ou o compartilhamento de uma configuração de execução.

Crie um arquivo de parâmetros no seu diretório de trabalho:

```bash
touch params.yaml
```

Abra-o no editor e adicione o caminho de saída:

```yaml title="params.yaml"
outdir: "s3://my-bucket/demo-results-cli"
```

Agora você pode executar o pipeline usando o perfil `test` (que fornece o samplesheet de `input`) e o arquivo de parâmetros (que fornece o `outdir`):

```bash
tw launch nf-core-demo --profile test --params-file params.yaml
```

??? success "Saída do comando"

    ```console
    Launching pipeline nf-core-demo
    Run name: focused_feynman
    https://cloud.seqera.io/orgs/my-org/workspaces/my-workspace/watch/<run-id>
    ```

Abra o link para confirmar que a execução aparece no painel **Runs**.

!!! tip "Dica"

    Você pode incluir o arquivo de parâmetros durante a etapa de configuração inicial se quiser definir alguns valores padrão, bem como algumas propriedades adicionais para corresponder ao que fizemos anteriormente pelo formulário web:

    ```bash
    tw pipelines add \
      --name nf-core-demo-gg \
      --description "Demo pipeline with defaults for testing" \
      --params-file params.yaml \
      --profile test \
      https://github.com/nf-core/demo
    ```

### Conclusão

Você sabe como adicionar qualquer pipeline Nextflow hospedado no GitHub ao seu workspace e executá-lo, tanto pela interface web preenchendo os parâmetros manualmente, quanto pelo `tw` CLI combinando um perfil com um arquivo de parâmetros.

---

## Resumo

Nesta parte você aprendeu a:

- Autenticar o `tw` CLI e executar um pipeline salvo pelo terminal
- Adicionar um novo pipeline do GitHub usando o CLI e verificar se ele aparece no Launchpad
- Executar um pipeline pela interface web do Seqera preenchendo os parâmetros obrigatórios manualmente
- Executar um pipeline pelo CLI usando um perfil Nextflow e um arquivo de parâmetros
