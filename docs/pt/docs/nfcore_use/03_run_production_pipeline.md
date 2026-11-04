# Parte 3: Executar um pipeline de produção

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Tradução assistida por IA - [saiba mais e sugira melhorias](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

Na [Parte 2](./02_configure_execution.md), você aprendeu como definir parâmetros e personalizar a configuração para o nf-core/demo.
Agora vamos aplicar o que você aprendeu a um pipeline de produção real, o nf-core/rnaseq.

---

## 1. Baixar e executar o nf-core/rnaseq

Até agora usamos o `nf-core/demo`, que é um pipeline mínimo desenvolvido para treinamento.
Agora vamos baixar um pipeline de produção real e executá-lo com seu perfil de teste.

O pipeline `nf-core/rnaseq` realiza as etapas principais da análise de RNA sequenciamento em massa: controle de qualidade, remoção de adaptadores, alinhamento de leituras e quantificação em nível de genes.
É provavelmente o pipeline nf-core mais amplamente utilizado até hoje.

### 1.1. Baixar o pipeline

Execute o seguinte comando para baixá-lo.

```bash
nextflow pull nf-core/rnaseq
```

??? success "Saída do comando"

    ```console
    Checking nf-core/rnaseq ...
     downloaded from https://github.com/nf-core/rnaseq.git - revision: 1f03b53ef7 [master]
    ```

O pipeline agora está armazenado em cache localmente e pronto para ser executado.

### 1.2. Executar o perfil de teste

Execute-o com o perfil de teste e Docker:

```bash
nextflow run nf-core/rnaseq -profile test,docker --outdir rnaseq-results
```

??? failure "Saída do comando"

    ```console
     N E X T F L O W   ~  version 26.04.4

    Launching `https://github.com/nf-core/rnaseq` [suspicious_dijkstra] revision: 1f03b53ef7 [master]

    ------------------------------------------------------
                                            ,--./,-.
            ___     __   __   __   ___     /,-._.--~'
      |\ | |__  __ /  ` /  \ |__) |__         }  {
      | \| |       \__, \__/ |  \ |___     \`-._,-`-,
                                            `._,._,'
      nf-core/rnaseq 3.26.0
    ------------------------------------------------------
    ...
    [-        ] NFCORE_RNASEQ:RNASEQ:FASTQ_QC_TRIM_FILTER_SETSTRANDEDNESS:FQ_LINT  -
    Plus 47 more processes waiting for tasks…

    Execution cancelled -- Finishing pending tasks before exit
    -[nf-core/rnaseq] Pipeline completed with errors-
    ERROR ~ Error executing process > 'NFCORE_RNASEQ:RNASEQ:FASTQ_QC_TRIM_FILTER_SETSTRANDEDNESS:FQ_LINT (WT_REP2)'

    Caused by:
      Process requirement exceeds available memory -- req: 12 GB; avail: 7.7 GB


    Command executed:

      fq lint \
          --disable-validator P001 \
          SRR6357072_1.fastq.gz SRR6357072_2.fastq.gz > WT_REP2.fq_lint.txt

    Command exit status:
      -

    Command output:
      (empty)

    Work dir:
      /workspaces/training/nfcore-use/work/xx/xxxxxxxxxxxxxxxxxxxxxx

    Container:
      quay.io/biocontainers/fq:0.12.0--h9ee0642_0

    Tip: you can replicate the issue by changing to the process work dir and entering the command `bash .command.run`

     -- Check '.nextflow.log' file for details
    ERROR ~ Pipeline failed. Please refer to troubleshooting docs for common issues: https://nf-co.re/docs/running/troubleshooting
    ```

A linha principal nesse erro é:

```console
Process requirement exceeds available memory -- req: 12 GB; avail: 7.7 GB
```

A máquina padrão do Codespaces tem 8 GB de RAM, que também é o padrão típico para o Docker Desktop.
O pipeline está solicitando 12 GB para o processo `FQ_LINT` — mais do que a máquina pode fornecer.

Esses 12 GB vêm do rótulo de recurso `process_low` definido em `conf/base.config`:

```groovy title="conf/base.config"
withLabel:process_low {
    cpus   = { 2     * task.attempt }
    memory = { 12.GB * task.attempt }
    time   = { 4.h   * task.attempt }
}
```

Uma opção seria usar um tipo de máquina maior, mas para fins de teste queremos poder executar em qualquer hardware disponível.
A abordagem mais adequada é substituir os padrões de recursos em um arquivo de configuração personalizado.

### 1.3. Executar novamente com uma configuração personalizada

Fornecemos um arquivo de configuração personalizado que substitui os padrões de recursos baseados em rótulos.

??? full-code "laptop.config"

    ```groovy title="laptop.config"
    process {
        withLabel: 'process_low' {
            cpus   = 2
            memory = 6.GB
        }
        withLabel: 'process_medium' {
            cpus   = 4
            memory = 6.GB
        }
        withLabel: 'process_high' {
            cpus   = 6
            memory = 6.GB
        }
        withLabel: 'process_high_memory' {
            memory = 6.GB
        }
    }
    ```

A [Parte 2](./02_configure_execution.md) apresentou `withName:` para direcionar um único processo pelo nome.
Aqui usamos `withLabel:` para direcionar todos os processos que compartilham um rótulo de uma só vez.

Este arquivo já está presente no seu diretório de trabalho.
Passe-o com `-c` para aplicar as substituições:

```bash
nextflow run nf-core/rnaseq -profile test,docker -c laptop.config --outdir rnaseq-results
```

??? success "Saída do comando (pipeline iniciando)"

    ```console
     N E X T F L O W   ~  version 26.04.4

    Launching `https://github.com/nf-core/rnaseq` [romantic_faraday] revision: 1f03b53ef7 [master]

    ------------------------------------------------------
                                            ,--./,-.
            ___     __   __   __   ___     /,-._.--~'
      |\ | |__  __ /  ` /  \ |__) |__         }  {
      | \| |       \__, \__/ |  \ |___     \`-._,-`-,
                                            `._,._,'
      nf-core/rnaseq 3.26.0
    ------------------------------------------------------
    ...
    executor >  local (7)
    [xx/xxxxxx] NFCORE_RNASEQ:RNASEQ:FQ_LINT (RAP1_IAA_30M_REP1)   | 3 of 5
    [xx/xxxxxx] NFCORE_RNASEQ:RNASEQ:FASTQC (RAP1_IAA_30M_REP1)    | 2 of 5
    ...
    ```

O pipeline agora está em execução, e você pode acompanhar as tarefas sendo concluídas uma a uma.
Neste conjunto de dados de teste mínimo, ele será concluído em 15–20 minutos, executando mais de 200 tarefas no total.

Experimentos reais de RNA-seq normalmente envolvem dezenas de amostras e são executados por horas ou dias.
O Nextflow oferece suporte a agendadores HPC (SLURM, PBS, LSF) e plataformas de nuvem (AWS, Google Cloud, Azure), que podem reduzir drasticamente o tempo de execução ao distribuir o trabalho entre muitos nós.
Configurar esses ambientes, no entanto, adiciona uma complexidade significativa.

A plataforma Seqera (desenvolvida pelos criadores do Nextflow) fornece uma interface web para lançar pipelines Nextflow em infraestrutura HPC ou de nuvem (seja a sua própria ou uma gerenciada para você), com capacidades de gerenciamento de computação e dados que simplificam o processo de execução de pipelines em escala.

!!! tip "Dica"

    Pesquisadores acadêmicos podem acessar a Seqera Platform gratuitamente por meio do [programa acadêmico da Seqera](https://seqera.io/academic-program/).

### Conclusão

Você baixou o `nf-core/rnaseq`, viu como os rótulos de recursos do nf-core funcionam e aprendeu a substituí-los com um arquivo de configuração personalizado.
Mais importante, você viu por que a execução local é um ponto de partida, e não um destino, para análises em escala real.

### O que vem a seguir?

Você cobriu os fundamentos da execução de pipelines nf-core.
Consulte os [Próximos passos](next_steps.md) para saber para onde ir a partir daqui.

---

## Resumo

Nesta parte você aprendeu a:

- Baixar e executar um pipeline em escala de produção (nf-core/rnaseq) e substituir seus rótulos de recursos padrão
