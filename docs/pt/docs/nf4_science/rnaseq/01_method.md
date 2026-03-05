# Parte 1: Visão geral do método

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Tradução assistida por IA - [saiba mais e sugira melhorias](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

Existem múltiplos métodos válidos para processar e analisar dados de RNAseq em bulk.
Para este curso, estamos seguindo o método descrito [aqui](https://www.bioinformatics.babraham.ac.uk/training/RNASeq_Course/Analysing%20RNA-Seq%20data%20Exercise.pdf) pelos Drs. Simon Andrews e Laura Biggins no [Babraham Institute](https://www.babraham.ac.uk/).

Nosso objetivo é desenvolver um fluxo de trabalho que implementa as seguintes etapas de processamento: executar controle de qualidade inicial nas leituras de uma amostra de RNAseq em bulk, remover sequências de adaptadores das leituras, alinhar as leituras a um genoma de referência e produzir um relatório abrangente de controle de qualidade (QC).

<figure class="excalidraw">
--8<-- "docs/en/docs/nf4_science/rnaseq/img/rnaseq-method-03.svg"
</figure>

- **FASTQC:** Realizar QC nos dados de leitura antes da remoção usando FastQC
- **TRIM_GALORE:** Remover sequências de adaptadores e realizar QC após a remoção usando Trim Galore (agrupa Cutadapt e FastQC)
- **HISAT2_ALIGN:** Alinhar leituras ao genoma de referência usando Hisat2
- **MULTIQC:** Gerar um relatório abrangente de QC usando MultiQC

### Métodos

Vamos mostrar como aplicar essas etapas de processamento em duas fases.
Primeiro começaremos com **processamento de amostra única** que executa as ferramentas de QC, remoção e alinhamento em uma amostra.
Depois estenderemos para **processamento de múltiplas amostras** que executa as mesmas ferramentas em múltiplas amostras e gera um relatório agregado de controle de qualidade.

Antes de mergulharmos na escrita de qualquer código de fluxo de trabalho para qualquer abordagem, vamos testar os comandos manualmente em alguns dados de teste.

### Conjunto de dados

Fornecemos os seguintes dados e recursos relacionados:

- **Dados de RNAseq** (`reads/`): arquivos FASTQ de seis amostras, reduzidos a uma pequena região para manter os tamanhos de arquivo baixos. Cada amostra tem leituras paired-end (dois arquivos por amostra), embora comecemos trabalhando apenas com leituras single-end.
- **Um genoma de referência** (`genome.fa`): uma pequena região do cromossomo humano 20 (de hg19/b37).
- **Planilhas CSV** (`single-end.csv` e `paired-end.csv`): arquivos listando os IDs e caminhos dos arquivos de dados de exemplo.

### Software

As quatro ferramentas principais envolvidas são [FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/) para coleta de métricas de controle de qualidade, [Trim Galore](https://www.bioinformatics.babraham.ac.uk/projects/trim_galore/) para remoção de adaptadores (agrupa Cutadapt e FastQC para QC pós-remoção), [HISAT2](http://daehwankimlab.github.io/hisat2/) para alinhamento com splicing a um genoma de referência, e [MultiQC](https://multiqc.info/) para geração de relatório agregado de QC.

Essas ferramentas não estão instaladas no ambiente GitHub Codespaces, então vamos usá-las via contêineres obtidos através do serviço Seqera Containers (veja [Hello Containers](../../hello_nextflow/05_hello_containers.md)).

!!! tip "Dica"

     Certifique-se de estar no diretório `nf4-science/rnaseq`. A última parte do caminho mostrada quando você digita `pwd` deve ser `rnaseq`.

---

## 1. Processamento de amostra única

Nesta seção testamos os comandos que processam uma única amostra de RNAseq: controle de qualidade, remoção de adaptadores e alinhamento a um genoma de referência.
Estes são os comandos que vamos encapsular em um fluxo de trabalho Nextflow na Parte 2 deste curso.

1. Executar QC inicial em um arquivo FASTQ usando FastQC
2. Remover sequências de adaptadores e executar QC pós-remoção usando Trim Galore
3. Alinhar as leituras removidas ao genoma de referência usando HISAT2

<figure class="excalidraw">
--8<-- "docs/en/docs/nf4_science/rnaseq/img/rnaseq-method-02.svg"
</figure>

Começamos testando esses comandos em apenas uma amostra.

### 1.1. QC e remoção de adaptadores

Primeiro, queremos executar os comandos de QC e remoção em um dos arquivos de dados de exemplo.

#### 1.1.1. Baixar o contêiner

Vamos baixar uma imagem de contêiner que tem tanto `fastqc` quanto `trim_galore` instalados:

```bash
docker pull community.wave.seqera.io/library/trim-galore:0.6.10--1bf8ca4e1967cd18
```

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

Se você não baixou esta imagem antes, pode levar um minuto para completar.
Uma vez concluído, você tem uma cópia local da imagem do contêiner.

#### 1.1.2. Iniciar o contêiner interativamente

Para executar o contêiner interativamente, use `docker run` com as flags `-it`.
A opção `-v ./data:/data` monta nosso diretório local `data/` para que possamos acessar os arquivos de entrada de dentro do contêiner.

```bash
docker run -it -v ./data:/data community.wave.seqera.io/library/trim-galore:0.6.10--1bf8ca4e1967cd18
```

??? success "Saída do comando"

    ```console
    (base) root@b645838b3314:/tmp#
    ```

Seu prompt mudará para algo como `(base) root@b645838b3314:/tmp#`, o que indica que você está agora dentro do contêiner.

Verifique se você pode ver os arquivos de dados de sequência em `/data/reads`:

```bash
ls /data/reads
```

??? abstract "Conteúdo do diretório"

    ```console
    ENCSR000COQ1_1.fastq.gz  ENCSR000COQ2_2.fastq.gz  ENCSR000COR2_1.fastq.gz  ENCSR000CPO1_2.fastq.gz
    ENCSR000COQ1_2.fastq.gz  ENCSR000COR1_1.fastq.gz  ENCSR000COR2_2.fastq.gz  ENCSR000CPO2_1.fastq.gz
    ENCSR000COQ2_1.fastq.gz  ENCSR000COR1_2.fastq.gz  ENCSR000CPO1_1.fastq.gz  ENCSR000CPO2_2.fastq.gz
    ```

Com isso, você está pronto para testar seu primeiro comando.

#### 1.1.3. Executar o comando FastQC

O método referenciado acima nos fornece a linha de comando para executar QC em um único arquivo.
Precisamos apenas fornecer o arquivo de entrada; a ferramenta gerará automaticamente arquivos de saída no mesmo diretório que os dados originais.

Execute o comando `fastqc` em um arquivo de dados:

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

??? abstract "Conteúdo do diretório"

    ```console
    /data/reads/ENCSR000COQ1_1_fastqc.html  /data/reads/ENCSR000COQ1_1_fastqc.zip
    ```

Você deve ver um relatório HTML e um arquivo ZIP contendo as métricas de QC.
Isso completa o teste da primeira etapa.

#### 1.1.4. Remover sequências de adaptadores com Trim Galore

Agora vamos executar `trim_galore`, que agrupa Cutadapt e FastQC, para remover as sequências de adaptadores e coletar métricas de QC pós-remoção.
Como observado acima, o software está incluído no mesmo contêiner, então nenhuma mudança é necessária.

O comando é direto; simplesmente precisamos adicionar a flag `--fastqc` para executar automaticamente uma etapa de coleta de QC após a conclusão da remoção.

```bash
trim_galore --fastqc /data/reads/ENCSR000COQ1_1.fastq.gz
```

??? success "Saída do comando"

    ```console hl_lines="54 55 56 58 59 60"
    Multicore support not enabled. Proceeding with single-core trimming.
    Path to Cutadapt set as: 'cutadapt' (default)
    Cutadapt seems to be working fine (tested command 'cutadapt --version')
    Cutadapt version: 4.9
    single-core operation.
    igzip command line interface 2.31.0
    igzip detected. Using igzip for decompressing

    No quality encoding type selected. Assuming that the data provided uses Sanger encoded Phred scores (default)



    AUTO-DETECTING ADAPTER TYPE
    ===========================
    Attempting to auto-detect adapter type from the first 1 million sequences of the first file (>> /data/reads/ENCSR000COQ1_1.fastq.gz <<)

    Found perfect matches for the following adapter sequences:
    Adapter type	Count	Sequence	Sequences analysed	Percentage
    Illumina	9	AGATCGGAAGAGC	27816	0.03
    smallRNA	0	TGGAATTCTCGG	27816	0.00
    Nextera	0	CTGTCTCTTATA	27816	0.00
    Using Illumina adapter for trimming (count: 9). Second best hit was smallRNA (count: 0)

    Writing report to 'ENCSR000COQ1_1.fastq.gz_trimming_report.txt'

    SUMMARISING RUN PARAMETERS
    ==========================
    Input filename: /data/reads/ENCSR000COQ1_1.fastq.gz
    Trimming mode: single-end
    Trim Galore version: 0.6.10
    Cutadapt version: 4.9
    Number of cores used for trimming: 1
    Quality Phred score cutoff: 20
    Quality encoding type selected: ASCII+33
    Adapter sequence: 'AGATCGGAAGAGC' (Illumina TruSeq, Sanger iPCR; auto-detected)
    Maximum trimming error rate: 0.1 (default)
    Minimum required adapter overlap (stringency): 1 bp
    Minimum required sequence length before a sequence gets removed: 20 bp
    Running FastQC on the data once trimming has completed
    Output file(s) will be GZIP compressed

    Cutadapt seems to be fairly up-to-date (version 4.9). Setting -j 1
    Writing final adapter and quality trimmed output to ENCSR000COQ1_1_trimmed.fq.gz


      >>> Now performing quality (cutoff '-q 20') and adapter trimming in a single pass for the adapter sequence: 'AGATCGGAAGAGC' from file /data/reads/ENCSR000COQ1_1.fastq.gz <<<
    This is cutadapt 4.9 with Python 3.12.7
    Command line parameters: -j 1 -e 0.1 -q 20 -O 1 -a AGATCGGAAGAGC /data/reads/ENCSR000COQ1_1.fastq.gz
    Processing single-end reads on 1 core ...
    Finished in 0.373 s (13.399 µs/read; 4.48 M reads/minute).

    === Summary ===

    Total reads processed:                  27,816
    Reads with adapters:                     9,173 (33.0%)
    Reads written (passing filters):        27,816 (100.0%)

    Total basepairs processed:     2,114,016 bp
    Quality-trimmed:                       0 bp (0.0%)
    Total written (filtered):      2,100,697 bp (99.4%)

    === Adapter 1 ===

    Sequence: AGATCGGAAGAGC; Type: regular 3'; Length: 13; Trimmed: 9173 times

    Minimum overlap: 1
    No. of allowed errors:
    1-9 bp: 0; 10-13 bp: 1

    Bases preceding removed adapters:
      A: 27.4%
      C: 37.4%
      G: 20.9%
      T: 14.3%
      none/other: 0.0%

    Overview of removed sequences
    length	count	expect	max.err	error counts
    1	6229	6954.0	0	6229
    2	2221	1738.5	0	2221
    3	581	434.6	0	581
    4	88	108.7	0	88
    5	33	27.2	0	33
    6	2	6.8	0	2
    7	1	1.7	0	1
    9	1	0.1	0	1
    10	2	0.0	1	2
    12	1	0.0	1	0 1
    14	4	0.0	1	3 1
    16	1	0.0	1	1
    19	1	0.0	1	1
    22	1	0.0	1	1
    29	4	0.0	1	0 4
    33	3	0.0	1	3

    RUN STATISTICS FOR INPUT FILE: /data/reads/ENCSR000COQ1_1.fastq.gz
    =============================================
    27816 sequences processed in total
    Sequences removed because they became shorter than the length cutoff of 20 bp:	0 (0.0%)


      >>> Now running FastQC on the data <<<

    application/gzip
    Started analysis of ENCSR000COQ1_1_trimmed.fq.gz
    Approx 5% complete for ENCSR000COQ1_1_trimmed.fq.gz
    Approx 10% complete for ENCSR000COQ1_1_trimmed.fq.gz
    Approx 15% complete for ENCSR000COQ1_1_trimmed.fq.gz
    Approx 20% complete for ENCSR000COQ1_1_trimmed.fq.gz
    Approx 25% complete for ENCSR000COQ1_1_trimmed.fq.gz
    Approx 30% complete for ENCSR000COQ1_1_trimmed.fq.gz
    Approx 35% complete for ENCSR000COQ1_1_trimmed.fq.gz
    Approx 40% complete for ENCSR000COQ1_1_trimmed.fq.gz
    Approx 45% complete for ENCSR000COQ1_1_trimmed.fq.gz
    Approx 50% complete for ENCSR000COQ1_1_trimmed.fq.gz
    Approx 55% complete for ENCSR000COQ1_1_trimmed.fq.gz
    Approx 60% complete for ENCSR000COQ1_1_trimmed.fq.gz
    Approx 65% complete for ENCSR000COQ1_1_trimmed.fq.gz
    Approx 70% complete for ENCSR000COQ1_1_trimmed.fq.gz
    Approx 75% complete for ENCSR000COQ1_1_trimmed.fq.gz
    Approx 80% complete for ENCSR000COQ1_1_trimmed.fq.gz
    Approx 85% complete for ENCSR000COQ1_1_trimmed.fq.gz
    Approx 90% complete for ENCSR000COQ1_1_trimmed.fq.gz
    Approx 95% complete for ENCSR000COQ1_1_trimmed.fq.gz
    Analysis complete for ENCSR000COQ1_1_trimmed.fq.gz
    ```

A saída é muito detalhada, então destacamos as linhas mais relevantes no exemplo acima.
Você pode encontrar os arquivos de saída no diretório de trabalho:

```bash
ls ENCSR000COQ1_1*
```

??? abstract "Conteúdo do diretório"

    ```console
    ENCSR000COQ1_1.fastq.gz_trimming_report.txt  ENCSR000COQ1_1_trimmed_fastqc.html
    ENCSR000COQ1_1_trimmed.fq.gz                 ENCSR000COQ1_1_trimmed_fastqc.zip
    ```

Isso inclui as leituras removidas, relatório de remoção e arquivos de QC pós-remoção.

#### 1.1.5. Mover os arquivos de saída

Qualquer coisa que permaneça dentro do contêiner ficará inacessível para trabalhos futuros, então precisamos mover esses arquivos para um diretório no sistema de arquivos montado.

```bash
mkdir /data/trimmed
mv ENCSR000COQ1_1* /data/trimmed
```

??? abstract "Conteúdo do diretório"

    ```console
    /data/trimmed
    ├── ENCSR000COQ1_1.fastq.gz_trimming_report.txt
    ├── ENCSR000COQ1_1_trimmed.fq.gz
    ├── ENCSR000COQ1_1_trimmed_fastqc.html
    └── ENCSR000COQ1_1_trimmed_fastqc.zip
    ```

Os arquivos agora estão acessíveis no seu sistema de arquivos normal.

#### 1.1.6. Sair do contêiner

Para sair do contêiner, digite `exit`.

```bash
exit
```

Seu prompt deve voltar ao normal; isso completa o teste das duas primeiras etapas.

### 1.2. Alinhar leituras ao genoma de referência

Em seguida, queremos executar o comando de alinhamento para alinhar as leituras de RNAseq removidas a um genoma de referência.

#### 1.2.1. Baixar o contêiner

Vamos baixar uma imagem de contêiner que tem `hisat2` e `samtools` instalados:

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

Você notará que algumas camadas mostram `Already exists` porque são compartilhadas com a imagem do contêiner Trim Galore que baixamos anteriormente.
Como resultado, este download deve ser mais rápido que o primeiro.

#### 1.2.2. Iniciar o contêiner interativamente

Inicie o contêiner interativamente, usando a mesma abordagem de antes com o URI do contêiner relevante substituído.

```bash
docker run -it -v ./data:/data community.wave.seqera.io/library/hisat2_samtools:5e49f68a37dc010e
```

Seu prompt mudará novamente para indicar que você está dentro do contêiner.

#### 1.2.3. Criar os arquivos de índice do genoma

HISAT2 requer que a referência do genoma seja fornecida em um formato muito específico e não pode simplesmente consumir o arquivo FASTA `genome.fa` que fornecemos, então vamos aproveitar esta oportunidade para criar os recursos relevantes.

```bash
hisat2-build /data/genome.fa genome_index
```

??? success "Saída do comando"

    ```console hl_lines="1 2 218"
    Settings:
      Output files: "genome_index.*.ht2"
      Line rate: 6 (line is 64 bytes)
      Lines per side: 1 (side is 64 bytes)
      Offset rate: 4 (one in 16)
      FTable chars: 10
      Strings: unpacked
      Local offset rate: 3 (one in 8)
      Local fTable chars: 6
      Local sequence length: 57344
      Local sequence overlap between two consecutive indexes: 1024
      Endianness: little
      Actual local endianness: little
      Sanity checking: disabled
      Assertions: disabled
      Random seed: 0
      Sizeofs: void*:8, int:4, long:8, size_t:8
    Input files DNA, FASTA:
      /data/genome.fa
    Reading reference sizes
      Time reading reference sizes: 00:00:00
    Calculating joined length
    Writing header
    Reserving space for joined string
    Joining reference sequences
      Time to join reference sequences: 00:00:00
      Time to read SNPs and splice sites: 00:00:00
    Using parameters --bmax 6542727 --dcv 1024
      Doing ahead-of-time memory usage test
      Passed!  Constructing with these parameters: --bmax 6542727 --dcv 1024
    Constructing suffix-array element generator
    Building DifferenceCoverSample
      Building sPrime
      Building sPrimeOrder
      V-Sorting samples
      V-Sorting samples time: 00:00:01
      Allocating rank array
      Ranking v-sort output
      Ranking v-sort output time: 00:00:00
      Invoking Larsson-Sadakane on ranks
      Invoking Larsson-Sadakane on ranks time: 00:00:00
      Sanity-checking and returning
    Building samples
    Reserving space for 12 sample suffixes
    Generating random suffixes
    QSorting 12 sample offsets, eliminating duplicates
    QSorting sample offsets, eliminating duplicates time: 00:00:00
    Multikey QSorting 12 samples
      (Using difference cover)
      Multikey QSorting samples time: 00:00:00
    Calculating bucket sizes
    Splitting and merging
      Splitting and merging time: 00:00:00
    Split 1, merged 7; iterating...
    Splitting and merging
      Splitting and merging time: 00:00:00
    Avg bucket size: 4.98493e+06 (target: 6542726)
    Converting suffix-array elements to index image
    Allocating ftab, absorbFtab
    Entering GFM loop
    Getting block 1 of 7
      Reserving size (6542727) for bucket 1
      Calculating Z arrays for bucket 1
      Entering block accumulator loop for bucket 1:
      bucket 1: 10%
      bucket 1: 20%
      bucket 1: 30%
      bucket 1: 40%
      bucket 1: 50%
      bucket 1: 60%
      bucket 1: 70%
      bucket 1: 80%
      bucket 1: 90%
      bucket 1: 100%
      Sorting block of length 3540952 for bucket 1
      (Using difference cover)
      Sorting block time: 00:00:01
    Returning block of 3540953 for bucket 1
    Getting block 2 of 7
      Reserving size (6542727) for bucket 2
      Calculating Z arrays for bucket 2
      Entering block accumulator loop for bucket 2:
      bucket 2: 10%
      bucket 2: 20%
      bucket 2: 30%
      bucket 2: 40%
      bucket 2: 50%
      bucket 2: 60%
      bucket 2: 70%
      bucket 2: 80%
      bucket 2: 90%
      bucket 2: 100%
      Sorting block of length 6195795 for bucket 2
      (Using difference cover)
      Sorting block time: 00:00:01
    Returning block of 6195796 for bucket 2
    Getting block 3 of 7
      Reserving size (6542727) for bucket 3
      Calculating Z arrays for bucket 3
      Entering block accumulator loop for bucket 3:
      bucket 3: 10%
      bucket 3: 20%
      bucket 3: 30%
      bucket 3: 40%
      bucket 3: 50%
      bucket 3: 60%
      bucket 3: 70%
      bucket 3: 80%
      bucket 3: 90%
      bucket 3: 100%
      Sorting block of length 6199288 for bucket 3
      (Using difference cover)
      Sorting block time: 00:00:01
    Returning block of 6199289 for bucket 3
    Getting block 4 of 7
      Reserving size (6542727) for bucket 4
      Calculating Z arrays for bucket 4
      Entering block accumulator loop for bucket 4:
      bucket 4: 10%
      bucket 4: 20%
      bucket 4: 30%
      bucket 4: 40%
      bucket 4: 50%
      bucket 4: 60%
      bucket 4: 70%
      bucket 4: 80%
      bucket 4: 90%
      bucket 4: 100%
      Sorting block of length 6454986 for bucket 4
      (Using difference cover)
      Sorting block time: 00:00:00
    Returning block of 6454987 for bucket 4
    Getting block 5 of 7
      Reserving size (6542727) for bucket 5
      Calculating Z arrays for bucket 5
      Entering block accumulator loop for bucket 5:
      bucket 5: 10%
      bucket 5: 20%
      bucket 5: 30%
      bucket 5: 40%
      bucket 5: 50%
      bucket 5: 60%
      bucket 5: 70%
      bucket 5: 80%
      bucket 5: 90%
      bucket 5: 100%
      Sorting block of length 3493181 for bucket 5
      (Using difference cover)
      Sorting block time: 00:00:00
    Returning block of 3493182 for bucket 5
    Getting block 6 of 7
      Reserving size (6542727) for bucket 6
      Calculating Z arrays for bucket 6
      Entering block accumulator loop for bucket 6:
      bucket 6: 10%
      bucket 6: 20%
      bucket 6: 30%
      bucket 6: 40%
      bucket 6: 50%
      bucket 6: 60%
      bucket 6: 70%
      bucket 6: 80%
      bucket 6: 90%
      bucket 6: 100%
      Sorting block of length 5875908 for bucket 6
      (Using difference cover)
      Sorting block time: 00:00:00
    Returning block of 5875909 for bucket 6
    Getting block 7 of 7
      Reserving size (6542727) for bucket 7
      Calculating Z arrays for bucket 7
      Entering block accumulator loop for bucket 7:
      bucket 7: 10%
      bucket 7: 20%
      bucket 7: 30%
      bucket 7: 40%
      bucket 7: 50%
      bucket 7: 60%
      bucket 7: 70%
      bucket 7: 80%
      bucket 7: 90%
      bucket 7: 100%
      Sorting block of length 3134429 for bucket 7
      (Using difference cover)
      Sorting block time: 00:00:00
    Returning block of 3134430 for bucket 7
    Exited GFM loop
    fchr[A]: 0
    fchr[C]: 9094775
    fchr[G]: 17470759
    fchr[T]: 25839994
    fchr[$]: 34894545
    Exiting GFM::buildToDisk()
    Returning from initFromVector
    Wrote 15826295 bytes to primary GFM file: genome_index.1.ht2
    Wrote 8723644 bytes to secondary GFM file: genome_index.2.ht2
    Re-opening _in1 and _in2 as input streams
    Returning from GFM constructor
    Returning from initFromVector
    Wrote 15353415 bytes to primary GFM file: genome_index.5.ht2
    Wrote 8883598 bytes to secondary GFM file: genome_index.6.ht2
    Re-opening _in5 and _in5 as input streams
    Returning from HGFM constructor
    Headers:
        len: 34894545
        gbwtLen: 34894546
        nodes: 34894546
        sz: 8723637
        gbwtSz: 8723637
        lineRate: 6
        offRate: 4
        offMask: 0xfffffff0
        ftabChars: 10
        eftabLen: 0
        eftabSz: 0
        ftabLen: 1048577
        ftabSz: 4194308
        offsLen: 2180910
        offsSz: 8723640
        lineSz: 64
        sideSz: 64
        sideGbwtSz: 48
        sideGbwtLen: 192
        numSides: 181743
        numLines: 181743
        gbwtTotLen: 11631552
        gbwtTotSz: 11631552
        reverse: 0
        linearFM: Yes
    Total time for call to driver() for forward index: 00:00:12
    ```

A saída é muito detalhada, então destacamos algumas linhas relevantes no exemplo acima.

Isso cria múltiplos arquivos de índice do genoma, que você pode encontrar no diretório de trabalho.

```bash
ls genome_index.*
```

??? abstract "Conteúdo do diretório"

    ```console
    genome_index.1.ht2  genome_index.3.ht2  genome_index.5.ht2  genome_index.7.ht2
    genome_index.2.ht2  genome_index.4.ht2  genome_index.6.ht2  genome_index.8.ht2
    ```

Precisaremos desses arquivos mais tarde, e gerá-los não é tipicamente algo que queremos fazer como parte de um fluxo de trabalho, então vamos gerar um tarball compactado com gzip contendo os arquivos de índice do genoma que podemos facilmente passar conforme necessário.

```bash
tar -czvf /data/genome_index.tar.gz genome_index.*
```

??? success "Saída do comando"

    ```console
    genome_index.1.ht2
    genome_index.2.ht2
    genome_index.3.ht2
    genome_index.4.ht2
    genome_index.5.ht2
    genome_index.6.ht2
    genome_index.7.ht2
    genome_index.8.ht2
    ```

Vamos mover o tarball `genome_index.tar.gz` resultante contendo os arquivos de índice do genoma para o diretório `data/` no nosso sistema de arquivos em alguns minutos.
Isso será útil na Parte 2 deste curso.

#### 1.2.4. Executar o comando de alinhamento

Agora podemos executar o comando de alinhamento, que realiza a etapa de alinhamento com `hisat2` e então redireciona a saída para `samtools` para escrever a saída como um arquivo BAM.

A entrada de dados de leitura é o arquivo `/data/trimmed/ENCSR000COQ1_1_trimmed.fq.gz` que geramos com `trim_galore` na etapa anterior.

```bash
hisat2 -x genome_index -U /data/trimmed/ENCSR000COQ1_1_trimmed.fq.gz \
    --new-summary --summary-file ENCSR000COQ1_1_trimmed.hisat2.log | \
    samtools view -bS -o ENCSR000COQ1_1_trimmed.bam
```

??? success "Saída do comando"

    ```console hl_lines="6"
    HISAT2 summary stats:
    	Total reads: 27816
    		Aligned 0 time: 1550 (5.57%)
    		Aligned 1 time: 25410 (91.35%)
    		Aligned >1 times: 856 (3.08%)
    	Overall alignment rate: 94.43%
    ```

Isso executa quase instantaneamente porque é um arquivo de teste muito pequeno.
Em escala real, isso poderia levar muito mais tempo.

Mais uma vez você pode encontrar os arquivos de saída no diretório de trabalho:

```bash
ls ENCSR000COQ1_1*
```

??? abstract "Conteúdo do diretório"

    ```console title="Output"
    ENCSR000COQ1_1_trimmed.bam  ENCSR000COQ1_1_trimmed.hisat2.log
    ```

O alinhamento produziu um arquivo BAM e um arquivo de log com estatísticas de alinhamento.

#### 1.2.5. Mover os arquivos de saída

Como antes, mova os arquivos de saída para um diretório no sistema de arquivos montado para que permaneçam acessíveis após sairmos do contêiner.

```bash
mkdir /data/aligned
mv ENCSR000COQ1_1* /data/aligned
```

Com isso feito, temos tudo o que precisamos.

#### 1.2.6. Sair do contêiner

Para sair do contêiner, digite `exit`.

```bash
exit
```

Seu prompt deve voltar ao normal.
Isso conclui a execução de teste do processamento de amostra única.

!!! example "Escreva como um fluxo de trabalho!"

    Sinta-se à vontade para seguir para a [Parte 2](./02_single-sample.md) imediatamente se quiser começar a implementar esta análise como um fluxo de trabalho Nextflow.
    Você só precisará voltar para completar a segunda rodada de testes antes de seguir para a Parte 3.

---

## 2. Agregação de QC de múltiplas amostras

Os comandos que acabamos de testar processam uma amostra por vez.
Na prática, tipicamente precisamos processar muitas amostras e então agregar resultados de QC em todas elas para avaliar a qualidade do conjunto de dados geral.

[MultiQC](https://multiqc.info/) é uma ferramenta que pesquisa através de diretórios por relatórios de QC de muitas ferramentas comuns de bioinformática e os agrega em um único relatório HTML abrangente.
Ela pode reconhecer saída de FastQC, Cutadapt (via Trim Galore) e HISAT2, entre muitas outras.

Aqui processamos duas amostras adicionais através das mesmas ferramentas por amostra, depois usamos MultiQC para agregar relatórios de QC em todas as três amostras.
Estes são os comandos que vamos encapsular em um fluxo de trabalho Nextflow na Parte 3 deste curso.

1. Executar QC e remoção em amostras adicionais usando Trim Galore
2. Executar alinhamento em amostras adicionais usando HISAT2
3. Agregar todos os relatórios de QC em um relatório abrangente usando MultiQC

<figure class="excalidraw">
--8<-- "docs/en/docs/nf4_science/rnaseq/img/rnaseq-method-03.svg"
</figure>

### 2.1. QC e remoção de amostras adicionais

Os comandos de QC e remoção por amostra são idênticos ao que executamos na seção 1.1.
Já baixamos a imagem do contêiner, então podemos iniciá-la diretamente.

#### 2.1.1. Iniciar o contêiner

Já baixamos esta imagem de contêiner na seção 1.1, então podemos iniciá-la diretamente:

```bash
docker run -it -v ./data:/data community.wave.seqera.io/library/trim-galore:0.6.10--1bf8ca4e1967cd18
```

Seu prompt muda para indicar que você está dentro do contêiner.

#### 2.1.2. Executar QC e remoção em amostras adicionais

Execute FastQC e Trim Galore em mais duas amostras, uma após a outra.

```bash
fastqc /data/reads/ENCSR000COQ2_1.fastq.gz
fastqc /data/reads/ENCSR000COR1_1.fastq.gz
```

```bash
trim_galore --fastqc /data/reads/ENCSR000COQ2_1.fastq.gz
trim_galore --fastqc /data/reads/ENCSR000COR1_1.fastq.gz
```

Uma vez concluído, você deve ter arquivos de saída do Trim Galore para ambas as amostras no diretório de trabalho.

#### 2.1.3. Mover os arquivos de saída

Mova os arquivos de saída do Trim Galore para o mesmo diretório que usamos na seção 1.

```bash
mv ENCSR000COQ2_1* ENCSR000COR1_1* /data/trimmed
```

??? abstract "Conteúdo do diretório"

    ```console
    /data/trimmed
    ├── ENCSR000COQ1_1.fastq.gz_trimming_report.txt
    ├── ENCSR000COQ1_1_trimmed.fq.gz
    ├── ENCSR000COQ1_1_trimmed_fastqc.html
    ├── ENCSR000COQ1_1_trimmed_fastqc.zip
    ├── ENCSR000COQ2_1.fastq.gz_trimming_report.txt
    ├── ENCSR000COQ2_1_trimmed.fq.gz
    ├── ENCSR000COQ2_1_trimmed_fastqc.html
    ├── ENCSR000COQ2_1_trimmed_fastqc.zip
    ├── ENCSR000COR1_1.fastq.gz_trimming_report.txt
    ├── ENCSR000COR1_1_trimmed.fq.gz
    ├── ENCSR000COR1_1_trimmed_fastqc.html
    └── ENCSR000COR1_1_trimmed_fastqc.zip
    ```

Os arquivos agora estão acessíveis no seu sistema de arquivos normal.

#### 2.1.4. Sair do contêiner

Para sair do contêiner, digite `exit`.

```bash
exit
```

Seu prompt deve voltar ao normal.

### 2.2. Alinhar amostras adicionais

Os comandos de alinhamento são idênticos ao que executamos na seção 1.2.
Precisamos extrair o índice do genoma do tarball que salvamos anteriormente, já que os arquivos de índice originais foram criados dentro de um contêiner que não existe mais.

#### 2.2.1. Iniciar o contêiner

Já baixamos esta imagem de contêiner na seção 1.2, então podemos iniciá-la diretamente:

```bash
docker run -it -v ./data:/data community.wave.seqera.io/library/hisat2_samtools:5e49f68a37dc010e
```

Seu prompt muda para indicar que você está dentro do contêiner.

#### 2.2.2. Extrair o índice do genoma

Extraia os arquivos de índice do genoma do tarball que salvamos no sistema de arquivos montado:

```bash
tar -xzf /data/genome_index.tar.gz
```

Isso restaura os arquivos `genome_index.*` no diretório de trabalho.

#### 2.2.3. Executar alinhamento em amostras adicionais

Execute o alinhamento HISAT2 nas duas amostras recém-removidas, uma após a outra.

```bash
hisat2 -x genome_index -U /data/trimmed/ENCSR000COQ2_1_trimmed.fq.gz \
    --new-summary --summary-file ENCSR000COQ2_1_trimmed.hisat2.log | \
    samtools view -bS -o ENCSR000COQ2_1_trimmed.bam
```

??? success "Saída do comando"

    ```console hl_lines="6"
    HISAT2 summary stats:
    	Total reads: 18736
    		Aligned 0 time: 1531 (8.17%)
    		Aligned 1 time: 16726 (89.27%)
    		Aligned >1 times: 479 (2.56%)
    	Overall alignment rate: 91.83%
    ```

```bash
hisat2 -x genome_index -U /data/trimmed/ENCSR000COR1_1_trimmed.fq.gz \
    --new-summary --summary-file ENCSR000COR1_1_trimmed.hisat2.log | \
    samtools view -bS -o ENCSR000COR1_1_trimmed.bam
```

??? success "Saída do comando"

    ```console hl_lines="6"
    HISAT2 summary stats:
    	Total reads: 38056
    		Aligned 0 time: 2311 (6.07%)
    		Aligned 1 time: 33289 (87.47%)
    		Aligned >1 times: 2456 (6.45%)
    	Overall alignment rate: 93.93%
    ```

Uma vez concluído, você deve ter arquivos BAM e de log para ambas as amostras no diretório de trabalho.

#### 2.2.4. Mover os arquivos de saída

Mova os arquivos de saída do alinhamento para o mesmo diretório que usamos na seção 1.

```bash
mv ENCSR000COQ2_1* ENCSR000COR1_1* /data/aligned
```

??? abstract "Conteúdo do diretório"

    ```console
    /data/aligned
    ├── ENCSR000COQ1_1_trimmed.bam
    ├── ENCSR000COQ1_1_trimmed.hisat2.log
    ├── ENCSR000COQ2_1_trimmed.bam
    ├── ENCSR000COQ2_1_trimmed.hisat2.log
    ├── ENCSR000COR1_1_trimmed.bam
    └── ENCSR000COR1_1_trimmed.hisat2.log
    ```

Os arquivos agora estão acessíveis no seu sistema de arquivos normal.

#### 2.2.5. Sair do contêiner

Para sair do contêiner, digite `exit`.

```bash
exit
```

Seu prompt deve voltar ao normal.

### 2.3. Gerar um relatório abrangente de QC

Agora que temos saída de QC, remoção e alinhamento para três amostras, podemos usar MultiQC para agregá-las em um único relatório.
MultiQC pesquisa através de diretórios por relatórios de QC compatíveis e agrega tudo o que encontrar.

#### 2.3.1. Baixar o contêiner

Vamos baixar uma imagem de contêiner que tem `multiqc` instalado:

```bash
docker pull community.wave.seqera.io/library/pip_multiqc:a3c26f6199d64b7c
```

??? success "Saída do comando"

    ```console
    a3c26f6199d64b7c: Pulling from library/pip_multiqc
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
    2ed162b168e8: Pull complete
    ca06fe148f21: Pull complete
    Digest: sha256:af0e9de56896805aa2a065f7650362956f4213d99e95314f6fec472c6a3bf091
    Status: Downloaded newer image for community.wave.seqera.io/library/pip_multiqc:a3c26f6199d64b7c
    community.wave.seqera.io/library/pip_multiqc:a3c26f6199d64b7c
    ```

Você notará que algumas camadas mostram `Already exists` porque são compartilhadas com as imagens de contêiner que baixamos anteriormente.
Como resultado, este download deve ser mais rápido que os anteriores.

#### 2.3.2. Iniciar o contêiner interativamente

Inicie o contêiner interativamente com o diretório de dados montado, como antes.

```bash
docker run -it -v ./data:/data community.wave.seqera.io/library/pip_multiqc:a3c26f6199d64b7c
```

Seu prompt mudará para indicar que você está dentro do contêiner.

#### 2.3.3. Executar o comando MultiQC

Execute `multiqc`, apontando-o para os diretórios onde armazenamos arquivos de saída relacionados a QC para todas as três amostras.
A flag `-n` define o nome do relatório de saída.

```bash
multiqc /data/reads /data/trimmed /data/aligned -n all_samples_QC
```

??? success "Saída do comando"

    ```console hl_lines="8 9 10 11 12"

    /// MultiQC 🔍 v1.32

          file_search | Search path: /data/reads
          file_search | Search path: /data/trimmed
          file_search | Search path: /data/aligned
            searching | ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━ 100% 36/36
               hisat2 | Found 3 reports
             cutadapt | Found 3 reports
               fastqc | Found 3 reports
        write_results | Data        : all_samples_QC_data
        write_results | Report      : all_samples_QC.html
              multiqc | MultiQC complete
    ```

Aqui vemos que a ferramenta encontrou relatórios de QC para todas as três amostras: o QC inicial do `fastqc`, os relatórios pós-remoção do `cutadapt` (via `trim_galore`) e os resumos de alinhamento do `hisat2`.

Os arquivos de saída estão no diretório de trabalho:

```bash
ls all_samples_QC*
```

??? abstract "Conteúdo do diretório"

    ```console
    all_samples_QC.html

    all_samples_QC_data:
    cutadapt_filtered_reads_plot.txt                     multiqc.log
    cutadapt_trimmed_sequences_plot_3_Counts.txt         multiqc.parquet
    cutadapt_trimmed_sequences_plot_3_Obs_Exp.txt        multiqc_citations.txt
    fastqc-status-check-heatmap.txt                      multiqc_cutadapt.txt
    fastqc_adapter_content_plot.txt                      multiqc_data.json
    fastqc_overrepresented_sequences_plot.txt            multiqc_fastqc.txt
    fastqc_per_base_n_content_plot.txt                   multiqc_general_stats.txt
    fastqc_per_base_sequence_quality_plot.txt            multiqc_hisat2.txt
    fastqc_per_sequence_gc_content_plot_Counts.txt       multiqc_software_versions.txt
    fastqc_per_sequence_gc_content_plot_Percentages.txt  multiqc_sources.txt
    fastqc_per_sequence_quality_scores_plot.txt
    fastqc_sequence_counts_plot.txt
    fastqc_sequence_duplication_levels_plot.txt
    fastqc_top_overrepresented_sequences_table.txt
    hisat2_se_plot.txt
    llms-full.txt
    ```

A saída principal é o relatório `all_samples_QC.html`, acompanhado por um diretório de dados contendo as métricas subjacentes.

#### 2.3.4. Mover os arquivos de saída

Mova o relatório e seu diretório de dados para o sistema de arquivos montado.

```bash
mkdir /data/multiqc
mv all_samples_QC* /data/multiqc
```

Os arquivos agora estão acessíveis no seu sistema de arquivos normal.

#### 2.3.5. Sair do contêiner

Para sair do contêiner, digite `exit`.

```bash
exit
```

Seu prompt deve voltar ao normal.
Isso conclui o teste de todos os comandos de processamento de RNAseq.

---

### Conclusão

Você sabe como executar os comandos FastQC, Trim Galore, HISAT2 e MultiQC em seus respectivos contêineres, incluindo como processar múltiplas amostras e agregar relatórios de QC.

### O que vem a seguir?

Faça uma pausa, depois siga para a [Parte 2](./02_single-sample.md) para aprender como encapsular esses mesmos comandos em fluxos de trabalho que usam contêineres para executar o trabalho.
