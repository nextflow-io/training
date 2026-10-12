# Parte 2: Gerenciar recursos de computação e falhas

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Tradução assistida por IA - [saiba mais e sugira melhorias](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

Na [Parte 1](./01_packaging_and_execution.md), você adaptou onde e como as tarefas de um pipeline são executadas.
Aqui você vai adaptar quanto recurso computacional cada tarefa recebe e o que acontece quando uma tarefa falha apesar da sua melhor estimativa de alocação.

---

## 1. Controlar alocações de recursos de computação

Por padrão, o Nextflow aloca um único CPU para cada processo via diretiva `cpus`, e não impõe um limite de memória a menos que você defina um:

```groovy title="Built-in configuration"
process {
    cpus = 1
}
```

Você já sabe pelo [Nextflow Run](../nextflow_run/index.md) que a configuração deste pipeline define `memory` como 1 GB para todos os processos.
Mas como você sabe quais valores realmente usar nos seus próprios pipelines?

### 1.1. Gerar um relatório de utilização de recursos

Você já gerou um relatório de execução com `-with-report` no [Nextflow Run](../nextflow_run/02_configure_pipeline.md).
Esse mesmo relatório é como você descobre quanto de CPU e memória seus processos realmente precisam: execute o fluxo de trabalho com algumas alocações padrão, registre o uso real e ajuste a partir daí.

```bash
nextflow run main.nf -with-report report-config-1.html
```

O relatório é um arquivo HTML que você pode abrir no navegador.
Ele detalha o tempo de execução e a utilização de recursos por processo, incluindo qual porcentagem dos recursos alocados foi realmente utilizada.
Veja o que ele mostra para `cowpy` com os padrões atuais (1 CPU, 1 GB de memória):

| Métrica               | Valor  |
| --------------------- | ------ |
| Uso de CPU            | 116%   |
| Pico de memória usada | 6,4 MB |
| Memória alocada       | 1 GB   |

`cowpy` usa bem menos de 1% da sua alocação de 1 GB; o `%cpu` acima de 100% significa apenas que ele usa brevemente mais do que o equivalente a um CPU de processamento dentro do contêiner, em rajadas curtas.

Consulte [Reports](https://nextflow.io/docs/latest/reports.html) para a lista completa de funcionalidades disponíveis.

### 1.2. Definir alocações de recursos para um processo específico

O relatório acima mostra que `cowpy` está confortavelmente dentro da sua alocação atual, mas digamos que você queira dar mais margem de segurança mesmo assim, por exemplo porque espera entradas maiores em produção.
Você pode substituir os padrões para um único processo com `withName`.

=== "Depois"

    ```groovy title="nextflow.config" linenums="6" hl_lines="5-6"
    process {
        memory = 1.GB
        withName: 'cowpy' {
            conda = 'conda-forge::cowpy==1.1.5'
            memory = 2.GB
            cpus = 2
        }
    }
    ```

=== "Antes"

    ```groovy title="nextflow.config" linenums="6"
    process {
        memory = 1.GB
        withName: 'cowpy' {
            conda = 'conda-forge::cowpy==1.1.5'
        }
    }
    ```

Com isso configurado, cada processo solicita 1 GB de memória e um único CPU, exceto `cowpy`, que solicita 2 GB e 2 CPUs (além da configuração do `conda` da [Parte 1](./01_packaging_and_execution.md)).

!!! info "Info"

    Se sua máquina tiver poucos CPUs e você alocar um número alto por processo, as chamadas de tarefa podem ficar na fila esperando umas pelas outras, já que o Nextflow não solicitará mais CPUs do que os disponíveis.

Execute novamente com um nome de arquivo de relatório diferente, para que você possa comparar antes e depois.

```bash
nextflow run main.nf -with-report report-config-2.html
```

??? success "Saída do comando"

    ```console
    N E X T F L O W   ~  version 26.04.4

    Launching `main.nf` [voluminous_venter] revision: c3c85dec78

    executor >  local (8)
    [a1/0e96d4] sayHello (1)       | 3 of 3 ✔
    [a3/7173a3] convertToUpper (2) | 3 of 3 ✔
    [4f/a8ae3d] collectGreetings   | 1 of 1 ✔
    [91/3724f8] cowpy              | 1 of 1 ✔

    Outputs:

      /workspaces/training/config-exec/results

      first_output:
        - full_pipeline/intermediates/Bonjour-output.txt
        - full_pipeline/intermediates/Hola-output.txt
        - full_pipeline/intermediates/Hello-output.txt

      uppercased:
        - full_pipeline/intermediates/UPPER-Bonjour-output.txt
        - full_pipeline/intermediates/UPPER-Hello-output.txt
        - full_pipeline/intermediates/UPPER-Hola-output.txt

      collected: full_pipeline/intermediates/COLLECTED-batch-output.txt

      batch_report: full_pipeline/batch-report.txt

      cowpy_art: full_pipeline/cowpy-COLLECTED-batch-output.txt
    ```

Comparando os dois relatórios para `cowpy`:

| Métrica               | Antes (1 CPU, 1 GB) | Depois (2 CPUs, 2 GB) |
| --------------------- | ------------------- | --------------------- |
| Pico de memória usada | 6,4 MB              | 6,4 MB                |
| Uso de CPU            | 116%                | 118%                  |

Dobrar a alocação não alterou o uso real em nada, o que indica que o 1 GB / 1 CPU original já era generoso para essa carga de trabalho de exemplo.
Em um pipeline real processando dados não triviais, você esperaria que os números em si diferissem significativamente entre os processos, que é exatamente por isso que você faz o profiling antes de decidir o que alocar, em vez de adivinhar.

### 1.3. Adicionar limites de recursos

Dependendo da sua infraestrutura de computação, pode haver restrições rígidas sobre o que você pode solicitar, por exemplo um limite máximo em todo o cluster.
A diretiva `resourceLimits` permite que você defina esses limites:

```groovy title="Syntax example"
process {
    resourceLimits = [
        memory: 750.GB,
        cpus: 200,
        time: 30.d
    ]
}
```

O Nextflow traduz esses limites para o que o executor de destino espera.
Se um processo solicitar mais do que o limite, a solicitação é reduzida ao limite em vez de ser rejeitada.

!!! warning "Aviso"

    Isso não é algo que você pode executar no ambiente de treinamento, pois requer infraestrutura HPC para ter efeito.

??? info "Configurações de referência institucionais"

    O projeto nf-core mantém uma [coleção de arquivos de configuração](https://nf-co.re/configs/) compartilhados por instituições ao redor do mundo, cobrindo uma ampla variedade de executores HPC e de nuvem.
    São um ponto de partida útil, independentemente de sua própria instituição estar entre elas ou não.

### Conclusão

Você sabe como gerar um relatório de profiling para avaliar a utilização de recursos, substituir alocações de recursos para um processo específico e limitar alocações com `resourceLimits`.

### O que vem a seguir?

Aprenda como fazer um pipeline se recuperar automaticamente quando uma tarefa falha, independentemente de sua estimativa de alocação de recursos estar correta ou não.

---

## 2. Lidar com falhas de tarefas com retentativas

O profiling informa o que um processo precisa na maior parte do tempo, mas cargas de trabalho reais variam: uma alocação confortável para a maioria das entradas ainda pode ser insuficiente para uma entrada excepcionalmente grande, e estimativas simplesmente podem estar erradas.
Em vez de deixar uma única tarefa com falha derrubar toda a execução, o Nextflow pode repetir uma tarefa com falha automaticamente, opcionalmente dando a ela mais recursos a cada tentativa.

### 2.1. Repetir uma tarefa com falha automaticamente

Para ver isso em ação, defina deliberadamente a alocação de memória do `cowpy` abaixo do que ele realmente precisa: lembre-se da [seção 1.1](#11-generate-a-resource-utilization-report) que ele atinge um pico de cerca de 6,4 MB, então 6 MB deve ser um pouco menos do que o suficiente.

=== "Depois"

    ```groovy title="nextflow.config" linenums="6" hl_lines="5-7"
    process {
        memory = 1.GB
        withName: 'cowpy' {
            conda = 'conda-forge::cowpy==1.1.5'
            memory = 6.MB
            errorStrategy = 'retry'
            maxRetries = 2
        }
    }
    ```

=== "Antes"

    ```groovy title="nextflow.config" linenums="6"
    process {
        memory = 1.GB
        withName: 'cowpy' {
            conda = 'conda-forge::cowpy==1.1.5'
            memory = 2.GB
            cpus = 2
        }
    }
    ```

`errorStrategy` informa ao Nextflow o que fazer quando uma tarefa falha: `'retry'` resubmete a tarefa em vez de parar todo o pipeline.
`maxRetries` limita quantas tentativas extras ela recebe antes de o Nextflow desistir.

```bash
nextflow run main.nf
```

??? failure "Saída do comando (resumida)"

    ```console
    N E X T F L O W   ~  version 26.04.4

    Launching `main.nf` [desperate_brazil] revision: c3c85dec78

    executor >  local (10)
    [67/fe1f49] sayHello (1)       | 3 of 3 ✔
    [8a/f13335] convertToUpper (1) | 3 of 3 ✔
    [39/1b24ed] collectGreetings   | 1 of 1 ✔
    [7a/d5eb6f] cowpy              | 0 of 1, retries: 2 ✘
    [9d/b79eb3] NOTE: Process `cowpy` terminated with an error exit status (137) -- Execution is retried (1)
    [6d/1d9d84] NOTE: Process `cowpy` terminated with an error exit status (137) -- Execution is retried (2)
    ERROR ~ Error executing process > 'cowpy'

    Caused by:
      Process `cowpy` terminated with an error exit status (137)

    Command executed:
      cat COLLECTED-batch-output.txt | cowpy -c "turkey" > cowpy-COLLECTED-batch-output.txt

    Command exit status:
      137

    Command output:
      (empty)

    Command error:
      /usr/local/bin/_activate_current_env.sh: line 35:    14 Killed                  micromamba activate "${ENV_NAME:-base}"

    Work dir:
      /workspaces/training/config-exec/work/7a/d5eb6feeac0eed18d95d3da7a7aeb4

    Tip: when you have fixed the problem you can continue the execution adding the option `-resume` to the run command line

    -- Check '.nextflow.log' file for details
    ```

O código de saída 137 é o sinal padrão para um encerramento por falta de memória: o contêiner não tinha memória suficiente para executar `cowpy`.
O Nextflow repetiu a tarefa duas vezes, três tentativas no total, correspondendo a `maxRetries = 2`.
Como a alocação de memória nunca mudou entre as tentativas, cada tentativa encontrou o mesmo obstáculo; uma vez esgotadas as retentativas, o Nextflow relata a falha completa e para o pipeline, saindo com um status diferente de zero.

Repetir por si só não resolve nada se a causa subjacente não mudar entre as tentativas.

### 2.2. Aumentar recursos a cada retentativa

Dentro de uma diretiva de processo, `task.attempt` contém o número da tentativa atual, começando em 1.
Você pode usá-lo em uma closure para escalar uma alocação de recursos para cima a cada retentativa.

=== "Depois"

    ```groovy title="nextflow.config" linenums="6" hl_lines="5 7"
    process {
        memory = 1.GB
        withName: 'cowpy' {
            conda = 'conda-forge::cowpy==1.1.5'
            memory = { 6.MB * task.attempt }
            errorStrategy = 'retry'
            maxRetries = 3
        }
    }
    ```

=== "Antes"

    ```groovy title="nextflow.config" linenums="6"
    process {
        memory = 1.GB
        withName: 'cowpy' {
            conda = 'conda-forge::cowpy==1.1.5'
            memory = 6.MB
            errorStrategy = 'retry'
            maxRetries = 2
        }
    }
    ```

Execute o fluxo de trabalho novamente:

```bash
nextflow run main.nf
```

??? success "Saída do comando (resumida)"

    ```console
    N E X T F L O W   ~  version 26.04.4

    Launching `main.nf` [grave_joliot] revision: c3c85dec78

    executor >  local (9)
    [0f/211b8a] sayHello (2)       | 3 of 3 ✔
    [26/301ea2] convertToUpper (3) | 3 of 3 ✔
    [22/5a895b] collectGreetings   | 1 of 1 ✔
    [e1/beee86] cowpy              | 1 of 1, retries: 1 ✔
    [b6/7aed6a] NOTE: Process `cowpy` terminated with an error exit status (137) -- Execution is retried (1)

    Outputs:

      /workspaces/training/config-exec/results

      first_output:
        - full_pipeline/intermediates/Hola-output.txt
        - full_pipeline/intermediates/Hello-output.txt
        - full_pipeline/intermediates/Bonjour-output.txt

      uppercased:
        - full_pipeline/intermediates/UPPER-Hola-output.txt
        - full_pipeline/intermediates/UPPER-Hello-output.txt
        - full_pipeline/intermediates/UPPER-Bonjour-output.txt

      collected: full_pipeline/intermediates/COLLECTED-batch-output.txt

      batch_report: full_pipeline/batch-report.txt

      cowpy_art: full_pipeline/cowpy-COLLECTED-batch-output.txt
    ```

A primeira tentativa ainda falha com 6 MB, mas a retentativa é executada com 12 MB (`6.MB * 2`) e tem sucesso, e o pipeline é concluído com todas as saídas publicadas.

!!! warning "Aviso"

    A saída do console ainda inclui uma linha `NOTE:` relatando a primeira tentativa com falha, mesmo que o pipeline como um todo tenha sido bem-sucedido: o Nextflow registra cada retentativa individualmente, mas uma falha com retentativa não afeta o resultado geral.
    Verifique o resumo `Outputs:`, ou o status de saída do comando, para confirmar se a execução realmente foi bem-sucedida.

Consulte [Dynamic computing resources](https://nextflow.io/docs/latest/process.html#dynamic-task-resources) na documentação do Nextflow para padrões de retentativa mais avançados, incluindo escalonamento baseado no erro específico que ocorreu.

### Conclusão

Você sabe como fazer um pipeline repetir automaticamente tarefas com falha e como escalar alocações de recursos a cada retentativa usando `task.attempt`.

### O que vem a seguir?

Siga para a [Parte 3](./03_profiles.md), onde você aprenderá como agrupar configurações como essas em perfis intercambiáveis.

---

## Resumo

Nesta parte você aprendeu a:

- Gerar um relatório de profiling de recursos e definir alocações de recursos por processo
- Limitar solicitações de recursos com `resourceLimits`
- Repetir automaticamente uma tarefa com falha usando `errorStrategy` e `maxRetries`
- Escalar uma alocação de recursos para cima a cada retentativa usando `task.attempt`
