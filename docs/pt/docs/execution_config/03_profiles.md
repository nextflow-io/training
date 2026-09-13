# Parte 3: Use profiles para alternar configurações

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Tradução assistida por IA - [saiba mais e sugira melhorias](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

Ao longo da [Parte 1](./01_packaging_and_execution.md) e da [Parte 2](./02_resources_and_retries.md), você acumulou algumas opções de configuração: empacotamento de software, plataforma de execução e alocação de recursos.
Na prática, você frequentemente vai querer alternar entre conjuntos completos dessas opções dependendo de onde está executando, por exemplo, um laptop para desenvolvimento e um cluster HPC para produção.

O Nextflow permite que você configure qualquer número de [profiles](https://nextflow.io/docs/latest/config.html#profiles) descrevendo diferentes configurações, e selecione um (ou vários) em tempo de execução com uma única flag.

Você já usou um: o profile `test` do [Nextflow Run](../nextflow_run/index.md) substitui os parâmetros de entrada por um conjunto pequeno e bem definido.
Agora você vai criar seus próprios profiles de infraestrutura e combiná-los com ele.

---

## 1. Crie profiles para diferentes ambientes

### 1.1. Configure os profiles

Adicione dois profiles ao `nextflow.config`: um para executar em um laptop comum com Docker, e outro para um cluster HPC universitário com um agendador Slurm e Conda.

=== "Depois"

    ```groovy title="nextflow.config" linenums="35" hl_lines="10-19"
    profiles {
        test {
            params.input = 'data/greetings.csv'
            params.batch = 'test'
            params.character = 'tux'
        }
        conda {
            docker.enabled = false
            conda.enabled = true
        }
        my_laptop {
            process.executor = 'local'
            docker.enabled = true
        }
        univ_hpc {
            process.executor = 'slurm'
            conda.enabled = true
            process.resourceLimits = [
                memory: 750.GB,
                cpus: 200,
                time: 30.d
            ]
        }
    }
    ```

=== "Antes"

    ```groovy title="nextflow.config" linenums="35"
    profiles {
        test {
            params.input = 'data/greetings.csv'
            params.batch = 'test'
            params.character = 'tux'
        }
        conda {
            docker.enabled = false
            conda.enabled = true
        }
    }
    ```

O profile `univ_hpc` também define limites de recursos, pois isso é tipicamente exigido em infraestruturas HPC compartilhadas.

### 1.2. Execute o fluxo de trabalho com um profile

Selecione um profile em tempo de execução com `-profile`.

```bash
nextflow run main.nf -profile my_laptop
```

??? success "Saída do comando"

    ```console
    N E X T F L O W   ~  version 26.04.4

    Launching `main.nf` [nice_kimura] revision: c3c85dec78

    executor >  local (8)
    [df/b86b2f] sayHello (2)       | 3 of 3 ✔
    [b2/e2cb9c] convertToUpper (3) | 3 of 3 ✔
    [af/2ecff7] collectGreetings   | 1 of 1 ✔
    [c6/713c96] cowpy              | 1 of 1 ✔

    Outputs:

      /workspaces/training/execution-config/results

      first_output:
        - full_pipeline/intermediates/Bonjour-output.txt
        - full_pipeline/intermediates/Hello-output.txt
        - full_pipeline/intermediates/Hola-output.txt

      uppercased:
        - full_pipeline/intermediates/UPPER-Bonjour-output.txt
        - full_pipeline/intermediates/UPPER-Hello-output.txt
        - full_pipeline/intermediates/UPPER-Hola-output.txt

      collected: full_pipeline/intermediates/COLLECTED-batch-output.txt

      batch_report: full_pipeline/batch-report.txt

      cowpy_art: full_pipeline/cowpy-COLLECTED-batch-output.txt
    ```

!!! warning "Aviso"

    O profile `univ_hpc` não vai funcionar no ambiente de treinamento, pois não há um agendador Slurm disponível.

Se você encontrar outras configurações que sempre pertencem juntas, adicione-as ao profile correspondente.
Você também pode criar profiles adicionais para agrupar qualquer outra combinação que precisar.

### 1.3. Execute com múltiplos profiles

Os profiles não são mutuamente exclusivos.
Você pode ativar vários ao mesmo tempo com `-profile <profile1>,<profile2>`.
Combine `my_laptop` com o profile `test` que você já conhece do Nextflow Run.

```bash
nextflow run main.nf -profile my_laptop,test
```

??? success "Saída do comando"

    ```console
    N E X T F L O W   ~  version 26.04.4

    Launching `main.nf` [romantic_volhard] revision: c3c85dec78

    executor >  local (8)
    [22/81de4f] sayHello (3)       | 3 of 3 ✔
    [93/e16b04] convertToUpper (1) | 3 of 3 ✔
    [ca/3e6a91] collectGreetings   | 1 of 1 ✔
    [16/0cc7e3] cowpy              | 1 of 1 ✔

    Outputs:

      /workspaces/training/execution-config/results

      first_output:
        - full_pipeline/intermediates/Hola-output.txt
        - full_pipeline/intermediates/Hello-output.txt
        - full_pipeline/intermediates/Bonjour-output.txt

      uppercased:
        - full_pipeline/intermediates/UPPER-Hello-output.txt
        - full_pipeline/intermediates/UPPER-Bonjour-output.txt
        - full_pipeline/intermediates/UPPER-Hola-output.txt

      collected: full_pipeline/intermediates/COLLECTED-test-output.txt

      batch_report: full_pipeline/test-report.txt

      cowpy_art: full_pipeline/cowpy-COLLECTED-test-output.txt
    ```

Os nomes dos arquivos individuais refletem corretamente `batch = 'test'` do profile `test` (`COLLECTED-test-output.txt`, e assim por diante).

Se você combinar profiles que definem a mesma opção, o Nextflow resolve o conflito usando o valor lido por último, ou seja, aquele que aparece mais tarde no arquivo.
Se as configurações conflitantes vierem de fontes de configuração completamente diferentes, a [ordem de precedência](https://www.nextflow.io/docs/latest/config.html) padrão é aplicada.

### Conclusão

Você sabe como definir profiles que agrupam configurações específicas de infraestrutura, selecionar um em tempo de execução com `-profile`, combinar múltiplos profiles em uma única execução, e como o Nextflow resolve conflitos quando mais de um profile define a mesma opção.

### O que vem a seguir?

Aprenda como inspecionar a configuração completamente resolvida antes de executar qualquer coisa.

---

## 2. Inspecione a configuração resolvida

Você já usou `nextflow config -profile test` no [Nextflow Run](../nextflow_run/02_configure_pipeline.md) para verificar o que um único profile resolve.
Esse comando se torna especialmente útil quando você está combinando múltiplos profiles: como você acabou de ver, quando dois profiles definem a mesma opção, pode ser complicado determinar manualmente qual valor prevalece.
O comando `nextflow config` resolve tudo isso por você, sem executar o pipeline.

### 2.1. Resolva a configuração padrão

```bash
nextflow config
```

??? success "Saída do comando"

    ```groovy
    params {
       input = 'data/greetings.csv'
       batch = 'batch'
       character = 'turkey'
    }

    docker {
       enabled = true
    }

    process {
       memory = '1 GB'
       withName:cowpy {
          conda = 'conda-forge::cowpy==1.1.5'
          memory = '2 GB'
          cpus = 2
       }
    }
    ```

Isso é exatamente o que seria aplicado se você executasse o pipeline sem flags adicionais.

### 2.2. Resolva a configuração com profiles ativados

Adicione os mesmos profiles que você usaria em uma execução real.

```bash
nextflow config -profile my_laptop,test
```

??? success "Saída do comando"

    ```groovy
    params {
       input = 'data/greetings.csv'
       batch = 'test'
       character = 'tux'
    }

    docker {
       enabled = true
    }

    process {
       memory = '1 GB'
       withName:cowpy {
          conda = 'conda-forge::cowpy==1.1.5'
          memory = '2 GB'
          cpus = 2
       }
       executor = 'local'
    }
    ```

Comparar os dois confirma o que mudou: `params.batch`, `params.character` e `process.executor` refletem os profiles `my_laptop,test`.
Isso se torna especialmente valioso para pipelines com muitas camadas de configuração, onde determinar manualmente as configurações resolvidas seria tedioso e sujeito a erros.

### Conclusão

Você sabe como usar `nextflow config` para inspecionar a configuração completamente resolvida para qualquer combinação de profiles, antes de executar qualquer coisa.

### O que vem a seguir?

Você cobriu os fundamentos da configuração de pipelines Nextflow.
Consulte o [Resumo do curso](next_steps.md) para saber o que fazer a partir daqui.

---

## Resumo

Nesta parte você aprendeu a:

- Definir profiles que agrupam configurações específicas de infraestrutura
- Combinar múltiplos profiles em uma única execução e entender como os conflitos entre eles são resolvidos
- Usar `nextflow config` para inspecionar a configuração completamente resolvida
