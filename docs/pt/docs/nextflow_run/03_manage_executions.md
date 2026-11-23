# Parte 3: Gerenciar execuções de fluxos de trabalho

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Tradução assistida por IA - [saiba mais e sugira melhorias](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

À medida que você executa e re-executa pipelines, você acumula histórico de execuções e diretórios `work/` antigos.
Na [Parte 1](./01_run_nextflow.md#23-use-resume-to-skip-completed-work) você já usou `-resume` para pular trabalhos que já haviam sido concluídos.
Aqui você aprenderá como gerar relatórios sobre uma execução, inspecionar o histórico de execuções passadas com [`nextflow log`](https://nextflow.io/docs/latest/reference/cli.html#log), e excluir diretórios de trabalho antigos que você não precisa mais com [`nextflow clean`](https://nextflow.io/docs/latest/reference/cli.html#clean).

---

## 1. Gerar relatórios do pipeline

O Nextflow pode gerar vários tipos de relatórios sobre uma execução, cada um adicionado com sua própria flag `-with-*`: um relatório de execução (`-with-report`), uma linha do tempo de execução (`-with-timeline`), um arquivo de rastreamento de tarefas (`-with-trace`), e um diagrama do fluxo de trabalho (`-with-dag`).
Vamos gerar os dois primeiros aqui; consulte [Execution reports](https://nextflow.io/docs/latest/reports.html) na referência do Nextflow para os demais.

### 1.1. Gerar um relatório de execução

Adicione `-with-report` a qualquer comando `nextflow run` para gerar um relatório HTML após a conclusão do pipeline:

```bash
nextflow run main.nf --input data/greetings.csv --character turkey -with-report
```

??? success "Saída do comando"

    ```console hl_lines="6 7 8 9"
    N E X T F L O W   ~  version 26.04.4

    Launching `main.nf` [intergalactic_dalembert] revision: ce74f81996

    executor >  local (8)
    [34/23f10f] sayHello (2)       | 3 of 3 ✔
    [af/cab69d] convertToUpper (3) | 3 of 3 ✔
    [9e/d73afb] collectGreetings   | 1 of 1 ✔
    [3c/392db0] cowpy              | 1 of 1 ✔

    Outputs:

      /workspaces/training/nextflow-run/results

      first_output:
        - batch/intermediates/Bonjour-output.txt
        - batch/intermediates/Hello-output.txt
        - batch/intermediates/Hola-output.txt

      uppercased:
        - batch/intermediates/UPPER-Bonjour-output.txt
        - batch/intermediates/UPPER-Hello-output.txt
        - batch/intermediates/UPPER-Hola-output.txt

      collected: batch/intermediates/COLLECTED-batch-output.txt

      batch_report: batch/batch-report.txt

      cowpy_art: batch/cowpy-COLLECTED-batch-output.txt
    ```

O Nextflow grava o relatório em um arquivo chamado `report-<timestamp>.html` no diretório de trabalho.
Abra-o em um navegador para ver um resumo da execução, uma tabela com cada tarefa e seu status e tempo de execução, e gráficos de uso de recursos divididos por processo.

A aba **Tasks** lista cada tarefa que o pipeline executou, com o nome do processo, status e uso de recursos:

![Tabela de tarefas do relatório de execução](img/execution_report_tasks.png)

O relatório é especialmente útil quando um pipeline demora mais do que o esperado ou uma tarefa falha: a tabela de tarefas mostra exatamente onde o tempo foi gasto e quais tarefas tiveram sucesso ou falharam.

### 1.2. Gerar uma linha do tempo de execução

Adicione `-with-timeline` a uma execução para obter uma visualização no estilo de gráfico de Gantt de quando cada tarefa foi executada:

```bash
nextflow run main.nf --input data/greetings.csv --character turkey -with-timeline
```

??? success "Saída do comando"

    ```console hl_lines="6 7 8 9"
    N E X T F L O W   ~  version 26.04.4

    Launching `main.nf` [jolly_noyce] revision: ce74f81996

    executor >  local (8)
    [ad/e92ef3] sayHello (3)       | 3 of 3 ✔
    [2a/df8a8d] convertToUpper (2) | 3 of 3 ✔
    [be/7fb72a] collectGreetings   | 1 of 1 ✔
    [63/dc9bd6] cowpy              | 1 of 1 ✔

    Outputs:

      /workspaces/training/nextflow-run/results

      first_output:
        - batch/intermediates/Bonjour-output.txt
        - batch/intermediates/Hello-output.txt
        - batch/intermediates/Hola-output.txt

      uppercased:
        - batch/intermediates/UPPER-Bonjour-output.txt
        - batch/intermediates/UPPER-Hello-output.txt
        - batch/intermediates/UPPER-Hola-output.txt

      collected: batch/intermediates/COLLECTED-batch-output.txt

      batch_report: batch/batch-report.txt

      cowpy_art: batch/cowpy-COLLECTED-batch-output.txt
    ```

O Nextflow grava a linha do tempo em um arquivo chamado `timeline-<timestamp>.html`.
Abra-o em um navegador para ver uma barra para cada tarefa, posicionada e dimensionada de acordo com quando ela foi executada e quanto tempo levou:

![Linha do tempo de execução](img/execution_timeline.png)

A linha do tempo torna visível de relance o formato de expansão-e-convergência da [Parte 1](./01_run_nextflow.md#31-run-the-workflow): as três tarefas `sayHello` são executadas em paralelo, depois as três tarefas `convertToUpper`, e então `collectGreetings` e `cowpy` são executadas uma após a outra, já que cada uma depende de tudo que veio antes.

### Conclusão

Você sabe como gerar um relatório de execução HTML com `-with-report` e uma linha do tempo de execução com `-with-timeline`, e onde procurar os outros tipos de relatório que o Nextflow suporta.

### O que vem a seguir?

Aprenda como inspecionar o histórico de execuções passadas.

---

## 2. Inspecionar o log de execuções passadas

Seja desenvolvendo um pipeline ou executando-o em produção, em algum momento você precisará consultar informações sobre execuções anteriores.

### 2.1. O arquivo de histórico

Toda vez que você inicia um fluxo de trabalho Nextflow, uma linha é gravada em um arquivo de log chamado `history`, dentro de um diretório oculto chamado `.nextflow` no diretório de trabalho atual.

??? abstract "Conteúdo do arquivo"

    ```txt title=".nextflow/history" linenums="1"
    2026-09-13 01:47:58	2.8s	determined_ramanujan	OK	ce74f81996b8ce999853fdbb25ba7969	03078950-5011-414f-992b-8fe1bd219ac2	nextflow run main.nf --input data/greetings.csv --character turkey
    2026-09-13 01:48:03	3s	prickly_cuvier	OK	ce74f81996b8ce999853fdbb25ba7969	5a9e8bf7-a7db-4c2c-b1d8-3bd4c7266528	nextflow run main.nf --input data/greetings.csv --character tux
    2026-09-13 01:48:09	2.3s	elegant_panini	OK	ce74f81996b8ce999853fdbb25ba7969	517ccc29-9db7-4d8d-8ea0-c8db6e9ac653	nextflow run main.nf --input data/greetings.csv --character stegosaurus
    2026-09-13 01:48:14	1.5s	deadly_lamport	OK	ce74f81996b8ce999853fdbb25ba7969	517ccc29-9db7-4d8d-8ea0-c8db6e9ac653	nextflow run main.nf --input data/greetings.csv --character stegosaurus -resume
    ```

Cada linha fornece o timestamp, duração, nome da execução, status, ID de revisão, ID de sessão e linha de comando completa de uma execução iniciada neste diretório.

Observe as duas últimas linhas: são duas invocações separadas (uma simples, outra com `-resume`) do mesmo comando exato, e elas compartilham o mesmo ID de sessão.
O ID de sessão só muda quando você inicia uma execução genuinamente nova; usar `-resume` o mantém, e é assim que o Nextflow sabe qual cache reutilizar.

### 2.2. Use `nextflow log` para uma visualização mais amigável

Ler o arquivo de histórico bruto funciona, mas `nextflow log` formata as mesmas informações com um cabeçalho:

```bash
nextflow log
```

??? success "Saída do comando"

    ```console linenums="1"
    TIMESTAMP          	DURATION	RUN NAME            	STATUS	REVISION ID	SESSION ID                          	COMMAND
    2026-09-13 01:47:58	2.8s    	determined_ramanujan	OK    	ce74f81996 	03078950-5011-414f-992b-8fe1bd219ac2	nextflow run main.nf --input data/greetings.csv --character turkey
    2026-09-13 01:48:03	3s      	prickly_cuvier      	OK    	ce74f81996 	5a9e8bf7-a7db-4c2c-b1d8-3bd4c7266528	nextflow run main.nf --input data/greetings.csv --character tux
    2026-09-13 01:48:09	2.3s    	elegant_panini      	OK    	ce74f81996 	517ccc29-9db7-4d8d-8ea0-c8db6e9ac653	nextflow run main.nf --input data/greetings.csv --character stegosaurus
    2026-09-13 01:48:14	1.5s    	deadly_lamport      	OK    	ce74f81996 	517ccc29-9db7-4d8d-8ea0-c8db6e9ac653	nextflow run main.nf --input data/greetings.csv --character stegosaurus -resume
    ```

O Nextflow agrupa as informações de cache que usa para `-resume` em `.nextflow/cache`, indexadas pelo ID de sessão.
É por isso que consultar o nome da execução ou o ID de sessão correto aqui é o primeiro passo sempre que você precisar investigar ou limpar uma execução passada.

### Conclusão

Você sabe onde o Nextflow registra o histórico de execuções passadas e como inspecioná-lo com `nextflow log`.

### O que vem a seguir?

Aprenda como remover diretórios de trabalho antigos que você não precisa mais.

---

## 3. Excluir diretórios de trabalho antigos

Cada execução deixa seus diretórios de tarefas em `work/`, mesmo depois de você ter copiado as saídas que importam para `results/`.
Execute pipelines suficientes durante o desenvolvimento e esses subdiretórios vão se acumulando, por isso o Nextflow fornece `nextflow clean` para remover os que você não precisa mais.

### 3.1. Determinar os critérios de exclusão

`nextflow clean` suporta várias formas de selecionar o que remover; consulte a [documentação de referência](https://www.nextflow.io/docs/latest/reference/cli.html#clean) para a lista completa.
Aqui você excluirá tudo de execuções anteriores a uma determinada execução, usando seu nome.

Consulte a execução mais recente que você deseja manter usando `nextflow log`; no [exemplo da seção 2.2](#22-use-nextflow-log-for-a-friendlier-view) é `elegant_panini`, a última execução simples antes da execução com `-resume`.
O nome da execução é a string de duas partes gerada automaticamente exibida na linha do console `Launching (...)`, ou na coluna `RUN NAME` do `nextflow log`.

### 3.2. Fazer uma execução de teste

Adicione `-n` primeiro para verificar o que um determinado comando excluiria sem realmente excluir nada:

```bash
nextflow clean -before elegant_panini -n
```

??? success "Saída do comando"

    ```console
    Would remove /workspaces/training/nextflow-run/work/e5/9ec3fe24deb5d1ffa478dff5e5e1d4
    Would remove /workspaces/training/nextflow-run/work/75/36b1ece86b213c331852d0a3dc659b
    Would remove /workspaces/training/nextflow-run/work/be/8157808f28699259b5ff37e4bf9f72
    Would remove /workspaces/training/nextflow-run/work/f3/d828b7ba6b315eb08e0f96168119f5
    Would remove /workspaces/training/nextflow-run/work/94/7a2d7ef67d81b6e66430f300816698
    Would remove /workspaces/training/nextflow-run/work/eb/e37c5c8638c45916de3d9794d6d22f
    Would remove /workspaces/training/nextflow-run/work/3f/f56a705a6370f6c615aff22d3f240a
    Would remove /workspaces/training/nextflow-run/work/25/4922fa35ba2c6087777ddce20e9f2f
    Would remove /workspaces/training/nextflow-run/work/a4/2de7bcc2f8fc974261f7280c1480f7
    Would remove /workspaces/training/nextflow-run/work/48/8c49bded446d4fcb78c8125318c5e3
    Would remove /workspaces/training/nextflow-run/work/8e/c9834a9bdd7c21de2991d959201a9d
    Would remove /workspaces/training/nextflow-run/work/d5/a9055dc8c817c6c1f436494685955c
    Would remove /workspaces/training/nextflow-run/work/27/9473df9e55e35c8869d9571b480926
    Would remove /workspaces/training/nextflow-run/work/d3/ebc066e1d844232b9f5ee0a7d87464
    Would remove /workspaces/training/nextflow-run/work/91/17ef1815ed06b8f177d0135d8ef0dc
    Would remove /workspaces/training/nextflow-run/work/12/0443f9fcff27a552d4bc73aaedf24a
    ```

São 16 diretórios de tarefas: as 8 tarefas da execução `turkey` mais as 8 da execução `tux`, exatamente o número esperado para duas execuções completas deste pipeline de quatro processos.
A própria execução `elegant_panini`, e as tarefas em cache que a execução com `-resume` reutilizou dela, são preservadas.

Sua saída listará nomes de diretórios diferentes, e o número de linhas depende de quantas execuções você realizou; se você não ver nenhuma linha, ou o nome da execução não corresponde a nenhum no seu log, ou não há nada para excluir antes dele.

### 3.3. Prosseguir com a exclusão

Quando a execução de teste parecer correta, execute o mesmo comando com `-f` em vez de `-n`:

```bash
nextflow clean -before elegant_panini -f
```

??? success "Saída do comando"

    ```console
    Removed /workspaces/training/nextflow-run/work/e5/9ec3fe24deb5d1ffa478dff5e5e1d4
    Removed /workspaces/training/nextflow-run/work/75/36b1ece86b213c331852d0a3dc659b
    Removed /workspaces/training/nextflow-run/work/be/8157808f28699259b5ff37e4bf9f72
    Removed /workspaces/training/nextflow-run/work/f3/d828b7ba6b315eb08e0f96168119f5
    Removed /workspaces/training/nextflow-run/work/94/7a2d7ef67d81b6e66430f300816698
    Removed /workspaces/training/nextflow-run/work/eb/e37c5c8638c45916de3d9794d6d22f
    Removed /workspaces/training/nextflow-run/work/3f/f56a705a6370f6c615aff22d3f240a
    Removed /workspaces/training/nextflow-run/work/25/4922fa35ba2c6087777ddce20e9f2f
    Removed /workspaces/training/nextflow-run/work/a4/2de7bcc2f8fc974261f7280c1480f7
    Removed /workspaces/training/nextflow-run/work/48/8c49bded446d4fcb78c8125318c5e3
    Removed /workspaces/training/nextflow-run/work/8e/c9834a9bdd7c21de2991d959201a9d
    Removed /workspaces/training/nextflow-run/work/d5/a9055dc8c817c6c1f436494685955c
    Removed /workspaces/training/nextflow-run/work/27/9473df9e55e35c8869d9571b480926
    Removed /workspaces/training/nextflow-run/work/d3/ebc066e1d844232b9f5ee0a7d87464
    Removed /workspaces/training/nextflow-run/work/91/17ef1815ed06b8f177d0135d8ef0dc
    Removed /workspaces/training/nextflow-run/work/12/0443f9fcff27a552d4bc73aaedf24a
    ```

`nextflow clean` esvazia os diretórios de tarefas, mas mantém os diretórios pai de dois caracteres (como `e5/`) no lugar.

!!! warning "Aviso"

    Excluir diretórios de trabalho de execuções passadas os remove do cache do Nextflow e apaga quaisquer saídas armazenadas apenas neles.
    Isso compromete a capacidade do Nextflow de retomar a execução sem re-executar os processos correspondentes, portanto, limpe apenas as execuções das quais você tem certeza de que não precisará retomar.
    É também por isso que vale a pena publicar tudo que importa em `results/` com `mode 'copy'` em vez de depender do diretório `work/` ou de um modo de publicação com `symlink`.

### Conclusão

Você sabe como remover diretórios de trabalho antigos com `nextflow clean`, e por que fazer isso implica abrir mão da capacidade de retomar a execução a partir dessas execuções.

### O que vem a seguir?

Aprenda como executar pipelines diretamente de repositórios remotos como o GitHub na [Parte 4](./04_remote_repositories.md).

---

## Resumo

Nesta parte você aprendeu a:

- Gerar um relatório de execução HTML com `-with-report` e uma linha do tempo de execução com `-with-timeline`
- Inspecionar o histórico de execuções passadas com `nextflow log`
- Remover diretórios de trabalho antigos com `nextflow clean`, e entender a troca que isso implica em relação ao resume
