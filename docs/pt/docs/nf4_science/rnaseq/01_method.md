# Parte 1: Visão geral do método e teste manual

Existem múltiplos métodos válidos para processar e analisar dados de RNAseq em bulk.
Para este curso, estamos seguindo o método descrito [aqui](https://www.bioinformatics.babraham.ac.uk/training/RNASeq_Course/Analysing%20RNA-Seq%20data%20Exercise.pdf) pelos Drs. Simon Andrews e Laura Biggins no [Babraham Institute](https://www.babraham.ac.uk/).

Nosso objetivo é desenvolver um fluxo de trabalho que implementa as seguintes etapas de processamento: executar controle de qualidade inicial nas leituras de uma amostra de RNAseq em bulk, remover sequências de adaptadores das leituras, alinhar as leituras a um genoma de referência e produzir um relatório abrangente de controle de qualidade (QC).

<figure class="excalidraw">
--8<-- "docs/nf4_science/rnaseq/img/preprocess.svg"
</figure>

- **FASTQC:** Realizar QC nos dados de leitura antes da remoção usando FastQC
- **TRIM_GALORE:** Remover sequências de adaptadores e realizar QC após a remoção usando Trim Galore (agrupa Cutadapt e FastQC)
- **HISAT2_ALIGN:** Alinhar leituras ao genoma de referência usando Hisat2
- **MULTIQC:** Gerar um relatório abrangente de QC usando MultiQC

No entanto, antes de mergulharmos na escrita de qualquer código de fluxo de trabalho, vamos testar os comandos manualmente em alguns dados de teste.
As ferramentas que precisamos não estão instaladas no ambiente GitHub Codespaces, então vamos usá-las via contêineres (veja [Hello Containers](../../hello_nextflow/05_hello_containers.md)).

!!! note "Nota"

     Certifique-se de estar no diretório `nf4-science/rnaseq`. A última parte do caminho mostrada quando você digita `pwd` deve ser `rnaseq`.

---

## 1. QC inicial e remoção de adaptadores

Vamos baixar uma imagem de contêiner que tem tanto `fastqc` quanto `trim_galore` instalados, iniciá-la interativamente e executar os comandos de remoção e QC em um dos arquivos de dados de exemplo.

### 1.1. Baixar o contêiner

```bash
docker pull community.wave.seqera.io/library/trim-galore:0.6.10--1bf8ca4e1967cd18
```

Isso fornece a seguinte saída no console enquanto o sistema baixa a imagem:

??? success "Saída do comando"

    ```console
    0.6.10--1bf8ca4e1967cd18: Pulling from library/trim-galore
    dafa2b0c44d2: Pull complete
    dec6b097362e: Pull complete
    f88da01cff0b: Pull complete
    4f4fb700ef54: Pull complete
    92dc97a3ef36: Pull complete
    403f74b0f85e: Pull complete
    10b8c00c10a5: Pull complete
    17dc7ea432cc: Pull complete
    bb36d6c3110d: Pull complete
    0ea1a16bbe82: Pull complete
    030a47592a0a: Pull complete
    32ec762be2d0: Pull complete
    d2cb90387285: Pull complete
    Digest: sha256:4f00e7b2a09f3c8d8a9ce955120e177152fb1e56f63a2a6e186088b1250d9907
    Status: Downloaded newer image for community.wave.seqera.io/library/trim-galore:0.6.10--1bf8ca4e1967cd18
    community.wave.seqera.io/library/trim-galore:0.6.10--1bf8ca4e1967cd18
    ```

### 1.2. Iniciar o contêiner interativamente

```bash
docker run -it -v ./data:/data community.wave.seqera.io/library/trim-galore:0.6.10--1bf8ca4e1967cd18
```

<!--
??? success "Saída do comando"

    ```console

    ```
-->

Seu prompt mudará para algo como `(base) root@b645838b3314:/tmp#`, o que indica que você está agora dentro do contêiner.

A parte `-v ./data:/data` do comando nos permitirá acessar o conteúdo do diretório `data/` de dentro do contêiner.

```bash
ls /data/reads
```

??? success "Saída do comando"

    ```console
    ENCSR000COQ1_1.fastq.gz  ENCSR000COQ2_2.fastq.gz  ENCSR000COR2_1.fastq.gz  ENCSR000CPO1_2.fastq.gz
    ENCSR000COQ1_2.fastq.gz  ENCSR000COR1_1.fastq.gz  ENCSR000COR2_2.fastq.gz  ENCSR000CPO2_1.fastq.gz
    ENCSR000COQ2_1.fastq.gz  ENCSR000COR1_2.fastq.gz  ENCSR000CPO1_1.fastq.gz  ENCSR000CPO2_2.fastq.gzO
    ```

### 1.3. Executar o primeiro comando `fastqc`

Vamos executar `fastqc` para coletar métricas de controle de qualidade nos dados de leitura.

```bash
fastqc /data/reads/ENCSR000COQ1_1.fastq.gz
```

??? success "Saída do comando"

    ```console
    application/gzip
    Started analysis of ENCSR000COQ1_1.fastq.gz
    Approx 5% complete for ENCSR000COQ1_1.fastq.gz
    Approx 10% complete for ENCSR000COQ1_1.fastq.gz
    Approx 15% complete for ENCSR000COQ1_1.fastq.gz
    Approx 20% complete for ENCSR000COQ1_1.fastq.gz
    Approx 25% complete for ENCSR000COQ1_1.fastq.gz
    Approx 30% complete for ENCSR000COQ1_1.fastq.gz
    Approx 35% complete for ENCSR000COQ1_1.fastq.gz
    Approx 40% complete for ENCSR000COQ1_1.fastq.gz
    Approx 45% complete for ENCSR000COQ1_1.fastq.gz
    Approx 50% complete for ENCSR000COQ1_1.fastq.gz
    Approx 55% complete for ENCSR000COQ1_1.fastq.gz
    Approx 60% complete for ENCSR000COQ1_1.fastq.gz
    Approx 65% complete for ENCSR000COQ1_1.fastq.gz
    Approx 70% complete for ENCSR000COQ1_1.fastq.gz
    Approx 75% complete for ENCSR000COQ1_1.fastq.gz
    Approx 80% complete for ENCSR000COQ1_1.fastq.gz
    Approx 85% complete for ENCSR000COQ1_1.fastq.gz
    Approx 90% complete for ENCSR000COQ1_1.fastq.gz
    Approx 95% complete for ENCSR000COQ1_1.fastq.gz
    Analysis complete for ENCSR000COQ1_1.fastq.gz
    ```

Isso deve executar muito rapidamente.
Você pode encontrar os arquivos de saída no mesmo diretório que os dados originais:

```bash
ls /data/reads/ENCSR000COQ1_1_fastqc*
```

<!-- switch to tree -->

```console title="Saída"
/data/reads/ENCSR000COQ1_1_fastqc.html  /data/reads/ENCSR000COQ1_1_fastqc.zip
```

### 1.4. Remover sequências de adaptadores com `trim_galore`

Agora vamos executar `trim_galore`, que agrupa Cutadapt e FastQC, para remover as sequências de adaptadores e coletar métricas de QC pós-remoção.

```bash
trim_galore --fastqc /data/reads/ENCSR000COQ1_1.fastq.gz
```

A flag `--fastqc` faz com que o comando execute automaticamente uma etapa de coleta de QC após a conclusão da remoção.

_A saída é muito detalhada, então o que segue está abreviado._

??? success "Saída do comando"

    ```console
    Multicore support not enabled. Proceeding with single-core trimming.
    Path to Cutadapt set as: 'cutadapt' (default)
    Cutadapt seems to be working fine (tested command 'cutadapt --version')
    Cutadapt version: 4.9
    single-core operation.
    igzip command line interface 2.31.0
    igzip detected. Using igzip for decompressing

    <...>

    Analysis complete for ENCSR000COQ1_1_trimmed.fq.gz
    ```

Você pode encontrar os arquivos de saída no diretório de trabalho:

```bash
ls ENCSR000COQ1_1*
```

```console title="Saída"
ENCSR000COQ1_1.fastq.gz_trimming_report.txt  ENCSR000COQ1_1_trimmed_fastqc.html
ENCSR000COQ1_1_trimmed.fq.gz                 ENCSR000COQ1_1_trimmed_fastqc.zip
```

### 1.5. Mover os arquivos de saída para o sistema de arquivos fora do contêiner

Qualquer coisa que permaneça dentro do contêiner ficará inacessível para trabalhos futuros, então vamos mover estes para um novo diretório.

```bash
mkdir /data/trimmed
mv ENCSR000COQ1_1* /data/trimmed
```

### 1.6. Sair do contêiner

```bash
exit
```

---

## 2. Alinhar as leituras ao genoma de referência

Vamos baixar uma imagem de contêiner que tem `hisat2` instalado, iniciá-la interativamente e executar o comando de alinhamento para alinhar os dados de RNAseq a um genoma de referência.

### 2.1. Baixar o contêiner `hisat2`

```bash
docker pull community.wave.seqera.io/library/hisat2_samtools:5e49f68a37dc010e
```

??? success "Saída do comando"

    ```console
    Unable to find image 'community.wave.seqera.io/library/hisat2_samtools:5e49f68a37dc010e' locally
    5e49f68a37dc010e: Pulling from library/hisat2_samtools
    dafa2b0c44d2: Already exists
    dec6b097362e: Already exists
    f88da01cff0b: Already exists
    4f4fb700ef54: Already exists
    92dc97a3ef36: Already exists
    403f74b0f85e: Already exists
    10b8c00c10a5: Already exists
    17dc7ea432cc: Already exists
    bb36d6c3110d: Already exists
    0ea1a16bbe82: Already exists
    030a47592a0a: Already exists
    e74ed5dd390b: Pull complete
    abfcf0185e51: Pull complete
    Digest: sha256:29d8e1a3172a2bdde7be813f7ebec22d331388194a7c0de872b4ccca4bed8f45
    Status: Downloaded newer image for community.wave.seqera.io/library/hisat2_samtools:5e49f68a37dc010e
    ```

### 2.2. Iniciar o contêiner `hisat2` interativamente

```bash
docker run -it -v ./data:/data community.wave.seqera.io/library/hisat2_samtools:5e49f68a37dc010e
```

O comando é o mesmo de antes, com o URI do contêiner relevante substituído.

### 2.3. Criar os arquivos de índice do genoma Hisat2

Hisat2 requer que a referência do genoma seja fornecida em um formato muito específico e não pode simplesmente consumir o arquivo FASTA `genome.fa` que fornecemos, então vamos aproveitar esta oportunidade para criar os recursos relevantes.

```bash
hisat2-build /data/genome.fa genome_index
```

A saída é muito detalhada, então o seguinte está abreviado:

<!-- TODO: switch to full output -->

??? success "Saída do comando"

    ```console
    Settings:
      Output files: "genome_index.*.ht2"
    <...>
    Total time for call to driver() for forward index: 00:00:16
    ```

Isso cria múltiplos arquivos de índice do genoma, que você pode encontrar no diretório de trabalho.

```bash
ls genome_index.*
```

```console title="Saída"
genome_index.1.ht2  genome_index.3.ht2  genome_index.5.ht2  genome_index.7.ht2
genome_index.2.ht2  genome_index.4.ht2  genome_index.6.ht2  genome_index.8.ht2
```

Vamos usá-los daqui a pouco, mas primeiro vamos gerar um tarball compactado com gzip com esses arquivos de índice do genoma; precisaremos deles mais tarde e gerar estes não é tipicamente algo que queremos fazer como parte de um fluxo de trabalho.

```bash
tar -czvf /data/genome_index.tar.gz genome_index.*
```

Isso armazena um tarball `genome_index.tar.gz` contendo os arquivos de índice do genoma no diretório `data/` no nosso sistema de arquivos, o que será útil na Parte 2 deste curso.

### 2.4. Executar o comando `hisat2`

Agora podemos executar o comando de alinhamento, que realiza a etapa de alinhamento com `hisat2` e então redireciona a saída para `samtools` para escrever a saída como um arquivo BAM.

A entrada de dados de leitura é o arquivo `/data/trimmed/ENCSR000COQ1_1_trimmed.fq.gz` que geramos com `trim_galore` na etapa anterior.

```bash
hisat2 -x genome_index -U /data/trimmed/ENCSR000COQ1_1_trimmed.fq.gz \
    --new-summary --summary-file ENCSR000COQ1_1_trimmed.hisat2.log | \
    samtools view -bS -o ENCSR000COQ1_1_trimmed.bam
```

??? success "Saída do comando"

    ```console
    HISAT2 summary stats:
            Total reads: 27816
                    Aligned 0 time: 1550 (5.57%)
                    Aligned 1 time: 25410 (91.35%)
                    Aligned >1 times: 856 (3.08%)
            Overall alignment rate: 94.43%
    ```

Isso executa quase instantaneamente porque é um arquivo de teste muito pequeno.
Em escala real, isso poderia levar muito mais tempo.

Mais uma vez, você pode encontrar os arquivos de saída no diretório de trabalho:

```bash
ls ENCSR000COQ1_1*
```

```console title="Saída"
ENCSR000COQ1_1_trimmed.bam  ENCSR000COQ1_1_trimmed.hisat2.log
```

### 2.5. Mover os arquivos de saída para o sistema de arquivos fora do contêiner

```bash
mkdir /data/aligned
mv ENCSR000COQ1_1* /data/aligned
```

### 2.6. Sair do contêiner

```bash
exit
```

---

## 3. Gerar um relatório abrangente de QC

Vamos baixar uma imagem de contêiner que tem `multiqc` instalado, iniciá-la interativamente e executar um comando de geração de relatório nos arquivos de relatório FastQC antes/depois.

### 3.1. Baixar o contêiner `multiqc`

```bash
docker pull community.wave.seqera.io/library/pip_multiqc:a3c26f6199d64b7c
```

??? success "Saída do comando"

    ```console
    ad8f247edb55897c: Pulling from library/pip_multiqc
    dafa2b0c44d2: Already exists
    dec6b097362e: Already exists
    f88da01cff0b: Already exists
    4f4fb700ef54: Already exists
    92dc97a3ef36: Already exists
    403f74b0f85e: Already exists
    10b8c00c10a5: Already exists
    17dc7ea432cc: Already exists
    bb36d6c3110d: Already exists
    0ea1a16bbe82: Already exists
    030a47592a0a: Already exists
    3f229294c69a: Pull complete
    5a5ad47fd84c: Pull complete
    Digest: sha256:0ebb1d9605395a7df49ad0eb366b21f46afd96a5090376b0d8941cf5294a895a
    Status: Downloaded newer image for community.wave.seqera.io/library/pip_multiqc:a3c26f6199d64b7c
    community.wave.seqera.io/library/pip_multiqc:a3c26f6199d64b7c
    ```

### 3.2. Iniciar o contêiner `multiqc` interativamente

```bash
docker run -it -v ./data:/data community.wave.seqera.io/library/pip_multiqc:a3c26f6199d64b7c
```

### 3.3. Executar o comando `multiqc`

```bash
multiqc /data/reads /data/trimmed /data/aligned -n ENCSR000COQ1_1_QC
```

??? success "Saída do comando"

    ```console

    /// MultiQC 🔍 v1.27.1

          file_search | Search path: /data/reads
          file_search | Search path: /data/trimmed
          file_search | Search path: /data/aligned
            searching | ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━ 100% 20/20
                hisat2 | Found 1 reports
              cutadapt | Found 1 reports
                fastqc | Found 1 reports
        write_results | Data        : ENCSR000COQ1_1_QC_data
        write_results | Report      : ENCSR000COQ1_1_QC.html
              multiqc | MultiQC complete
    ```

MultiQC é capaz de pesquisar em diretórios por relatórios de QC compatíveis e agregará tudo o que encontrar.

Aqui vemos que a ferramenta encontrou todos os três relatórios de QC que geramos: o QC inicial que fizemos com `fastqc`, o relatório pós-remoção do `cutadapt` (feito via `trim_galore`) e o QC pós-alinhamento produzido por `hisat2`.

Os arquivos de saída estão mais uma vez no diretório de trabalho:

```bash
ls ENCSR000COQ1_1_QC*
```

```console title="Saída"
ENCSR000COQ1_1_QC.html

ENCSR000COQ1_1_QC_data:
cutadapt_filtered_reads_plot.txt                     fastqc_top_overrepresented_sequences_table.txt
cutadapt_trimmed_sequences_plot_3_Counts.txt         hisat2_se_plot.txt
cutadapt_trimmed_sequences_plot_3_Obs_Exp.txt        multiqc.log
fastqc-status-check-heatmap.txt                      multiqc_citations.txt
fastqc_adapter_content_plot.txt                      multiqc_cutadapt.txt
fastqc_per_base_n_content_plot.txt                   multiqc_data.json
fastqc_per_base_sequence_quality_plot.txt            multiqc_fastqc.txt
fastqc_per_sequence_gc_content_plot_Counts.txt       multiqc_general_stats.txt
fastqc_per_sequence_gc_content_plot_Percentages.txt  multiqc_hisat2.txt
fastqc_per_sequence_quality_scores_plot.txt          multiqc_software_versions.txt
fastqc_sequence_counts_plot.txt                      multiqc_sources.txt
fastqc_sequence_duplication_levels_plot.txt
```

### 3.4. Mover os arquivos de saída para o sistema de arquivos fora do contêiner

```bash
mkdir /data/final_qc
mv ENCSR000COQ1_1_QC** /data/final_qc
```

### 3.5. Sair do contêiner

```bash
exit
```

---

### Conclusão

Você testou todos os comandos individuais interativamente nos contêineres relevantes.

### Qual é o próximo passo?

Aprenda como encapsular esses mesmos comandos em um fluxo de trabalho de múltiplas etapas que usa contêineres para executar o trabalho.
