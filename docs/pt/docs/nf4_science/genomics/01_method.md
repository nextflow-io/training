# Parte 1: Visão geral do método e testes manuais

Chamada de variantes é um método de análise genômica que visa identificar variações em uma sequência genômica em relação a um genoma de referência.
Aqui vamos usar ferramentas e métodos projetados para chamar variantes germinativas curtas, _i.e._ SNPs e indels, em dados de sequenciamento de genoma completo.

![Pipeline GATK](img/gatk-pipeline.png)

Um pipeline completo de chamada de variantes normalmente envolve muitas etapas, incluindo mapeamento para a referência (às vezes chamado de alinhamento do genoma) e filtragem e priorização de variantes.
Para simplificar, neste curso vamos focar apenas na parte de chamada de variantes.

### Métodos

Vamos mostrar duas maneiras de aplicar chamada de variantes a amostras de sequenciamento de genoma completo para identificar SNPs e indels germinativos.
Primeiro começaremos com uma **abordagem por amostra** simples que chama variantes independentemente de cada amostra.
Em seguida, mostraremos uma **abordagem de chamada conjunta** mais sofisticada que analisa múltiplas amostras juntas, produzindo resultados mais precisos e informativos.

Antes de mergulharmos na escrita de qualquer código de fluxo de trabalho para qualquer uma das abordagens, vamos testar os comandos manualmente em alguns dados de teste.

### Conjunto de dados

Fornecemos os seguintes dados e recursos relacionados:

- **Um genoma de referência** consistindo de uma pequena região do cromossomo 20 humano (de hg19/b37) e seus arquivos acessórios (índice e dicionário de sequência).
- **Três amostras de sequenciamento de genoma completo** correspondentes a um trio familiar (mãe, pai e filho), que foram reduzidas a uma pequena fatia de dados no cromossomo 20 para manter os tamanhos de arquivo pequenos.
  Estes são dados de sequenciamento Illumina de reads curtos que já foram mapeados para o genoma de referência, fornecidos em formato [BAM](https://samtools.github.io/hts-specs/SAMv1.pdf) (Binary Alignment Map, uma versão comprimida de SAM, Sequence Alignment Map).
- **Uma lista de intervalos genômicos**, i.e. coordenadas no genoma onde nossas amostras têm dados adequados para chamar variantes, fornecida em formato BED.

### Software

As duas principais ferramentas envolvidas são [Samtools](https://www.htslib.org/), um kit de ferramentas amplamente usado para manipular arquivos de alinhamento de sequência, e [GATK](https://gatk.broadinstitute.org/) (Genome Analysis Toolkit), um conjunto de ferramentas para descoberta de variantes desenvolvido no Broad Institute.

Essas ferramentas não estão instaladas no ambiente GitHub Codespaces, então vamos usá-las via contêineres (veja [Hello Containers](../../hello_nextflow/05_hello_containers.md)).

!!! note "Nota"

     Certifique-se de estar no diretório `nf4-science/genomics` para que a última parte do caminho mostrado quando você digita `pwd` seja `genomics`.

---

## 1. Chamada de variantes por amostra

A chamada de variantes por amostra processa cada amostra independentemente: o chamador de variantes examina os dados de sequenciamento de uma amostra por vez e identifica posições onde a amostra difere da referência.

Nesta seção testamos os dois comandos que compõem a abordagem de chamada de variantes por amostra: indexar um arquivo BAM com Samtools e chamar variantes com GATK HaplotypeCaller.
Estes são os comandos que vamos envolver em um fluxo de trabalho Nextflow na Parte 2 deste curso.

1. Gerar um arquivo de índice para um arquivo de entrada BAM usando [Samtools](https://www.htslib.org/)
2. Executar o GATK HaplotypeCaller no arquivo BAM indexado para gerar chamadas de variantes por amostra em VCF (Variant Call Format)

<figure class="excalidraw">
--8<-- "docs/en/docs/nf4_science/genomics/img/hello-gatk-1.svg"
</figure>

Começamos testando os dois comandos em apenas uma amostra.

### 1.1. Indexar um arquivo de entrada BAM com Samtools

Arquivos de índice são uma característica comum de formatos de arquivo de bioinformática; eles contêm informações sobre a estrutura do arquivo principal que permite que ferramentas como GATK acessem um subconjunto dos dados sem ter que ler todo o arquivo.
Isso é importante por causa de quão grandes esses arquivos podem ficar.

Arquivos BAM são frequentemente fornecidos sem um índice, então o primeiro passo em muitos fluxos de trabalho de análise é gerar um usando `samtools index`.

Vamos baixar um contêiner Samtools, iniciá-lo interativamente e executar o comando `samtools index` em um dos arquivos BAM.

#### 1.1.1. Baixar o contêiner Samtools

Execute o comando `docker pull` para baixar a imagem do contêiner Samtools:

```bash
docker pull community.wave.seqera.io/library/samtools:1.20--b5dfbd93de237464
```

??? success "Saída do comando"

    ```console
    1.20--b5dfbd93de237464: Pulling from library/samtools
    6360b3717211: Pull complete
    2ec3f7ad9b3c: Pull complete
    7716ca300600: Pull complete
    4f4fb700ef54: Pull complete
    8c61d418774c: Pull complete
    03dae77ff45c: Pull complete
    aab7f787139d: Pull complete
    4f4fb700ef54: Pull complete
    837d55536720: Pull complete
    897362c12ca7: Pull complete
    3893cbe24e91: Pull complete
    d1b61e94977b: Pull complete
    c72ff66fb90f: Pull complete
    0e0388f29b6d: Pull complete
    Digest: sha256:bbfc45b4f228975bde86cba95e303dd94ecf2fdacea5bfb2e2f34b0d7b141e41
    Status: Downloaded newer image for community.wave.seqera.io/library/samtools:1.20--b5dfbd93de237464
    community.wave.seqera.io/library/samtools:1.20--b5dfbd93de237464
    ```

Se você não baixou esta imagem antes, pode levar um minuto para completar.
Uma vez concluído, você tem uma cópia local da imagem do contêiner.

#### 1.1.2. Iniciar o contêiner Samtools interativamente

Para executar o contêiner interativamente, use `docker run` com as flags `-it`.
A opção `-v ./data:/data` monta o diretório local `data` dentro do contêiner para que as ferramentas possam acessar os arquivos de entrada.

```bash
docker run -it -v ./data:/data community.wave.seqera.io/library/samtools:1.20--b5dfbd93de237464
```

Seu prompt muda para algo como `(base) root@a1b2c3d4e5f6:/tmp#`, indicando que você está agora dentro do contêiner.
Os arquivos de dados são acessíveis em `/data`.

#### 1.1.3. Executar o comando de indexação

A [documentação do Samtools](https://www.htslib.org/doc/samtools-index.html) nos fornece a linha de comando para executar para indexar um arquivo BAM.

Precisamos apenas fornecer o arquivo de entrada; a ferramenta gerará automaticamente um nome para a saída anexando `.bai` ao nome do arquivo de entrada.

```bash
samtools index /data/bam/reads_mother.bam
```

??? abstract "Conteúdo do diretório"

    ```console
    data/bam/
    ├── reads_father.bam
    ├── reads_mother.bam
    ├── reads_mother.bam.bai
    └── reads_son.bam
    ```

Você deve agora ver um arquivo chamado `reads_mother.bam.bai` no mesmo diretório do arquivo de entrada BAM original.

#### 1.1.4. Sair do contêiner Samtools

Para sair do contêiner, digite `exit`.

```bash
exit
```

Seu prompt deve agora estar de volta ao que era antes de você iniciar o contêiner.

### 1.2. Chamar variantes com GATK HaplotypeCaller

Vamos baixar um contêiner GATK, iniciá-lo interativamente e executar o comando `gatk HaplotypeCaller` no arquivo BAM que acabamos de indexar.

#### 1.2.1. Baixar o contêiner GATK

Execute o comando `docker pull` para baixar a imagem do contêiner GATK:

```bash
docker pull community.wave.seqera.io/library/gatk4:4.5.0.0--730ee8817e436867
```

??? success "Saída do comando"

    Algumas camadas mostram `Already exists` porque são compartilhadas com a imagem do contêiner Samtools que baixamos anteriormente.

    ```console
    4.5.0.0--730ee8817e436867: Pulling from library/gatk4
    6360b3717211: Already exists
    2ec3f7ad9b3c: Already exists
    7716ca300600: Already exists
    4f4fb700ef54: Already exists
    8c61d418774c: Already exists
    03dae77ff45c: Already exists
    aab7f787139d: Already exists
    4f4fb700ef54: Already exists
    837d55536720: Already exists
    897362c12ca7: Already exists
    3893cbe24e91: Already exists
    d1b61e94977b: Already exists
    e5c558f54708: Pull complete
    087cce32d294: Pull complete
    Digest: sha256:e33413b9100f834fcc62fd5bc9edc1e881e820aafa606e09301eac2303d8724b
    Status: Downloaded newer image for community.wave.seqera.io/library/gatk4:4.5.0.0--730ee8817e436867
    community.wave.seqera.io/library/gatk4:4.5.0.0--730ee8817e436867
    ```

Isso deve ser mais rápido que o primeiro pull porque as duas imagens de contêiner compartilham a maioria de suas camadas.

#### 1.2.2. Iniciar o contêiner GATK interativamente

Inicie o contêiner GATK interativamente com o diretório de dados montado, assim como fizemos para Samtools.

```bash
docker run -it -v ./data:/data community.wave.seqera.io/library/gatk4:4.5.0.0--730ee8817e436867
```

Seu prompt muda para indicar que você está agora dentro do contêiner GATK.

#### 1.2.3. Executar o comando de chamada de variantes

A [documentação do GATK](https://gatk.broadinstitute.org/hc/en-us/articles/21905025322523-HaplotypeCaller) nos fornece a linha de comando para executar para realizar chamada de variantes em um arquivo BAM.

Precisamos fornecer o arquivo de entrada BAM (`-I`), bem como o genoma de referência (`-R`), um nome para o arquivo de saída (`-O`) e uma lista de intervalos genômicos a analisar (`-L`).

No entanto, não precisamos especificar o caminho para o arquivo de índice; a ferramenta procurará automaticamente por ele no mesmo diretório, com base na convenção estabelecida de nomenclatura e co-localização.
O mesmo se aplica aos arquivos acessórios do genoma de referência (arquivos de índice e dicionário de sequência, `*.fai` e `*.dict`).

```bash
gatk HaplotypeCaller \
        -R /data/ref/ref.fasta \
        -I /data/bam/reads_mother.bam \
        -O reads_mother.vcf \
        -L /data/ref/intervals.bed
```

??? success "Saída do comando"

    A ferramenta produz saída de log detalhada. As linhas destacadas confirmam a conclusão bem-sucedida.

    ```console hl_lines="37 51 56 57"
    Using GATK jar /opt/conda/share/gatk4-4.5.0.0-0/gatk-package-4.5.0.0-local.jar
    Running:
        java -Dsamjdk.use_async_io_read_samtools=false -Dsamjdk.use_async_io_write_samtools=true -Dsamjdk.use_async_io_write_tribble=false -Dsamjdk.compression_level=2 -jar /opt/conda/share/gatk4-4.5.0.0-0/gatk-package-4.5.0.0-local.jar HaplotypeCaller -R /data/ref/ref.fasta -I /data/bam/reads_mother.bam -O reads_mother.vcf -L /data/ref/intervals.bed
    00:27:50.687 INFO  NativeLibraryLoader - Loading libgkl_compression.so from jar:file:/opt/conda/share/gatk4-4.5.0.0-0/gatk-package-4.5.0.0-local.jar!/com/intel/gkl/native/libgkl_compression.so
    00:27:50.854 INFO  HaplotypeCaller - ------------------------------------------------------------
    00:27:50.858 INFO  HaplotypeCaller - The Genome Analysis Toolkit (GATK) v4.5.0.0
    00:27:50.858 INFO  HaplotypeCaller - For support and documentation go to https://software.broadinstitute.org/gatk/
    00:27:50.858 INFO  HaplotypeCaller - Executing as root@a1fe8ff42d07 on Linux v6.10.14-linuxkit amd64
    00:27:50.858 INFO  HaplotypeCaller - Java runtime: OpenJDK 64-Bit Server VM v17.0.11-internal+0-adhoc..src
    00:27:50.859 INFO  HaplotypeCaller - Start Date/Time: February 8, 2026 at 12:27:50 AM GMT
    00:27:50.859 INFO  HaplotypeCaller - ------------------------------------------------------------
    00:27:50.859 INFO  HaplotypeCaller - ------------------------------------------------------------
    00:27:50.861 INFO  HaplotypeCaller - HTSJDK Version: 4.1.0
    00:27:50.861 INFO  HaplotypeCaller - Picard Version: 3.1.1
    00:27:50.861 INFO  HaplotypeCaller - Built for Spark Version: 3.5.0
    00:27:50.862 INFO  HaplotypeCaller - HTSJDK Defaults.COMPRESSION_LEVEL : 2
    00:27:50.862 INFO  HaplotypeCaller - HTSJDK Defaults.USE_ASYNC_IO_READ_FOR_SAMTOOLS : false
    00:27:50.862 INFO  HaplotypeCaller - HTSJDK Defaults.USE_ASYNC_IO_WRITE_FOR_SAMTOOLS : true
    00:27:50.863 INFO  HaplotypeCaller - HTSJDK Defaults.USE_ASYNC_IO_WRITE_FOR_TRIBBLE : false
    00:27:50.864 INFO  HaplotypeCaller - Deflater: IntelDeflater
    00:27:50.864 INFO  HaplotypeCaller - Inflater: IntelInflater
    00:27:50.864 INFO  HaplotypeCaller - GCS max retries/reopens: 20
    00:27:50.864 INFO  HaplotypeCaller - Requester pays: disabled
    00:27:50.865 INFO  HaplotypeCaller - Initializing engine
    00:27:50.991 INFO  FeatureManager - Using codec BEDCodec to read file file:///data/ref/intervals.bed
    00:27:51.016 INFO  IntervalArgumentCollection - Processing 6369 bp from intervals
    00:27:51.029 INFO  HaplotypeCaller - Done initializing engine
    00:27:51.040 INFO  NativeLibraryLoader - Loading libgkl_utils.so from jar:file:/opt/conda/share/gatk4-4.5.0.0-0/gatk-package-4.5.0.0-local.jar!/com/intel/gkl/native/libgkl_utils.so
    00:27:51.042 INFO  NativeLibraryLoader - Loading libgkl_smithwaterman.so from jar:file:/opt/conda/share/gatk4-4.5.0.0-0/gatk-package-4.5.0.0-local.jar!/com/intel/gkl/native/libgkl_smithwaterman.so
    00:27:51.042 INFO  SmithWatermanAligner - Using AVX accelerated SmithWaterman implementation
    00:27:51.046 INFO  HaplotypeCallerEngine - Disabling physical phasing, which is supported only for reference-model confidence output
    00:27:51.063 INFO  NativeLibraryLoader - Loading libgkl_pairhmm_omp.so from jar:file:/opt/conda/share/gatk4-4.5.0.0-0/gatk-package-4.5.0.0-local.jar!/com/intel/gkl/native/libgkl_pairhmm_omp.so
    00:27:51.085 INFO  IntelPairHmm - Flush-to-zero (FTZ) is enabled when running PairHMM
    00:27:51.086 INFO  IntelPairHmm - Available threads: 10
    00:27:51.086 INFO  IntelPairHmm - Requested threads: 4
    00:27:51.086 INFO  PairHMM - Using the OpenMP multi-threaded AVX-accelerated native PairHMM implementation
    00:27:51.128 INFO  ProgressMeter - Starting traversal
    00:27:51.136 INFO  ProgressMeter -        Current Locus  Elapsed Minutes     Regions Processed   Regions/Minute
    00:27:51.882 WARN  InbreedingCoeff - InbreedingCoeff will not be calculated at position 20_10037292_10066351:3480 and possibly subsequent; at least 10 samples must have called genotypes
    00:27:52.969 INFO  HaplotypeCaller - 7 read(s) filtered by: MappingQualityReadFilter
    0 read(s) filtered by: MappingQualityAvailableReadFilter
    0 read(s) filtered by: MappedReadFilter
    0 read(s) filtered by: NotSecondaryAlignmentReadFilter
    0 read(s) filtered by: NotDuplicateReadFilter
    0 read(s) filtered by: PassesVendorQualityCheckReadFilter
    0 read(s) filtered by: NonZeroReferenceLengthAlignmentReadFilter
    0 read(s) filtered by: GoodCigarReadFilter
    0 read(s) filtered by: WellformedReadFilter
    7 total reads filtered out of 1867 reads processed
    00:27:52.971 INFO  ProgressMeter - 20_10037292_10066351:13499              0.0                    35           1145.7
    00:27:52.971 INFO  ProgressMeter - Traversal complete. Processed 35 total regions in 0.0 minutes.
    00:27:52.976 INFO  VectorLoglessPairHMM - Time spent in setup for JNI call : 0.003346916
    00:27:52.976 INFO  PairHMM - Total compute time in PairHMM computeLogLikelihoods() : 0.045731709
    00:27:52.977 INFO  SmithWatermanAligner - Total compute time in native Smith-Waterman : 0.02 sec
    00:27:52.981 INFO  HaplotypeCaller - Shutting down engine
    [February 8, 2026 at 12:27:52 AM GMT] org.broadinstitute.hellbender.tools.walkers.haplotypecaller.HaplotypeCaller done. Elapsed time: 0.04 minutes.
    Runtime.totalMemory()=203423744
    ```

O arquivo de saída `reads_mother.vcf` é criado dentro do seu diretório de trabalho no contêiner, então você não o verá no explorador de arquivos do VS Code a menos que você altere o caminho do arquivo de saída.
No entanto, é um arquivo de teste pequeno, então você pode usar `cat` para abri-lo e visualizar o conteúdo.
Se você rolar até o início do arquivo, encontrará um cabeçalho composto de muitas linhas de metadados, seguido por uma lista de chamadas de variantes, uma por linha.

??? abstract "Conteúdo do arquivo"

    ```console title="reads_mother.vcf" linenums="26"
    #CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	reads_mother
    20_10037292_10066351	3480	.	C	CT	503.03	.	AC=2;AF=1.00;AN=2;DP=23;ExcessHet=0.0000;FS=0.000;MLEAC=2;MLEAF=1.00;MQ=60.00;QD=27.95;SOR=1.179	GT:AD:DP:GQ:PL	1/1:0,18:18:54:517,54,0
    20_10037292_10066351	3520	.	AT	A	609.03	.	AC=2;AF=1.00;AN=2;DP=18;ExcessHet=0.0000;FS=0.000;MLEAC=2;MLEAF=1.00;MQ=60.00;QD=33.83;SOR=0.693	GT:AD:DP:GQ:PL	1/1:0,18:18:54:623,54,0
    20_10037292_10066351	3529	.	T	A	155.64	.	AC=1;AF=0.500;AN=2;BaseQRankSum=-0.544;DP=21;ExcessHet=0.0000;FS=1.871;MLEAC=1;MLEAF=0.500;MQ=60.00;MQRankSum=0.000;QD=7.78;ReadPosRankSum=-1.158;SOR=1.034	GT:AD:DP:GQ:PL	0/1:12,8:20:99:163,0,328
    ```

Cada linha descreve uma possível variante identificada nos dados de sequenciamento da amostra. Para orientação sobre como interpretar o formato VCF, veja [este artigo útil](https://www.ebi.ac.uk/training/online/courses/human-genetic-variation-introduction/variant-identification-and-analysis/understanding-vcf-format/).

O arquivo de saída VCF é acompanhado por um arquivo de índice chamado `reads_mother.vcf.idx` que foi criado automaticamente pelo GATK.
Ele tem a mesma função que o arquivo de índice BAM, permitir que ferramentas busquem e recuperem subconjuntos de dados sem carregar o arquivo inteiro.

#### 1.2.4. Sair do contêiner GATK

Para sair do contêiner, digite `exit`.

```bash
exit
```

Seu prompt deve estar de volta ao normal.
Isso conclui o teste de chamada de variantes por amostra.

---

## 2. Chamada conjunta em uma coorte

A abordagem de chamada de variantes que acabamos de usar gera chamadas de variantes por amostra.
Isso é adequado para analisar variantes de cada amostra isoladamente, mas fornece informações limitadas.
Muitas vezes é mais interessante ver como as chamadas de variantes diferem entre múltiplas amostras.
O GATK oferece um método alternativo chamado chamada de variantes conjunta para esse propósito.

A chamada de variantes conjunta envolve gerar um tipo especial de saída de variante chamado GVCF (para Genomic VCF) para cada amostra, então combinar os dados GVCF de todas as amostras e executar uma análise estatística de 'genotipagem conjunta'.

![Análise conjunta](img/joint-calling.png)

O que é especial sobre o GVCF de uma amostra é que ele contém registros resumindo estatísticas de dados de sequência sobre todas as posições na área alvo do genoma, não apenas as posições onde o programa encontrou evidência de variação.
Isso é crítico para o cálculo de genotipagem conjunta ([leitura adicional](https://gatk.broadinstitute.org/hc/en-us/articles/360035890431-The-logic-of-joint-calling-for-germline-short-variants)).

O GVCF é produzido pelo GATK HaplotypeCaller, a mesma ferramenta que acabamos de testar, com um parâmetro adicional (`-ERC GVCF`).
Combinar os GVCFs é feito com GATK GenomicsDBImport, que combina as chamadas por amostra em um armazenamento de dados (análogo a um banco de dados).
A análise de 'genotipagem conjunta' real é então feita com GATK GenotypeGVCFs.

Aqui testamos os comandos necessários para gerar GVCFs e executar genotipagem conjunta.
Estes são os comandos que vamos envolver em um fluxo de trabalho Nextflow na Parte 3 deste curso.

1. Gerar um arquivo de índice para cada arquivo de entrada BAM usando Samtools
2. Executar o GATK HaplotypeCaller em cada arquivo de entrada BAM para gerar um GVCF de chamadas de variantes genômicas por amostra
3. Coletar todos os GVCFs e combiná-los em um armazenamento de dados GenomicsDB
4. Executar genotipagem conjunta no armazenamento de dados GVCF combinado para produzir um VCF de nível de coorte

<figure class="excalidraw">
--8<-- "docs/en/docs/nf4_science/genomics/img/hello-gatk-2.svg"
</figure>

Agora precisamos testar todos esses comandos, começando com a indexação de todos os três arquivos BAM.

### 2.1. Indexar arquivos BAM para todas as três amostras

Na primeira seção acima, indexamos apenas um arquivo BAM.
Agora precisamos indexar todas as três amostras para que o GATK HaplotypeCaller possa processá-las.

#### 2.1.1. Iniciar o contêiner Samtools interativamente

Já baixamos a imagem do contêiner Samtools, então podemos iniciá-lo diretamente:

```bash
docker run -it -v ./data:/data community.wave.seqera.io/library/samtools:1.20--b5dfbd93de237464
```

Seu prompt muda para indicar que você está dentro do contêiner, com o diretório de dados montado como antes.

#### 2.1.2. Executar o comando de indexação nas três amostras

Execute o comando de indexação em cada um dos três arquivos BAM:

```bash
samtools index /data/bam/reads_mother.bam
samtools index /data/bam/reads_father.bam
samtools index /data/bam/reads_son.bam
```

??? abstract "Conteúdo do diretório"

    ```console
    data/bam/
    ├── reads_father.bam
    ├── reads_father.bam.bai
    ├── reads_mother.bam
    ├── reads_mother.bam.bai
    ├── reads_son.bam
    └── reads_son.bam.bai
    ```

Isso deve produzir os arquivos de índice no mesmo diretório dos arquivos BAM correspondentes.

#### 2.1.3. Sair do contêiner Samtools

Para sair do contêiner, digite `exit`.

```bash
exit
```

Seu prompt deve estar de volta ao normal.

### 2.2. Gerar GVCFs para todas as três amostras

Para executar a etapa de genotipagem conjunta, precisamos de GVCFs para todas as três amostras.

#### 2.2.1. Iniciar o contêiner GATK interativamente

Já baixamos a imagem do contêiner GATK anteriormente, então podemos iniciá-lo diretamente:

```bash
docker run -it -v ./data:/data community.wave.seqera.io/library/gatk4:4.5.0.0--730ee8817e436867
```

Seu prompt muda para indicar que você está dentro do contêiner GATK.

#### 2.2.2. Executar o comando de chamada de variantes com a opção GVCF

Para produzir um VCF genômico (GVCF), adicionamos a opção `-ERC GVCF` ao comando base, que ativa o modo GVCF do HaplotypeCaller.

Também alteramos a extensão do arquivo de saída de `.vcf` para `.g.vcf`.
Tecnicamente isso não é um requisito, mas é uma convenção fortemente recomendada.

```bash
gatk HaplotypeCaller \
        -R /data/ref/ref.fasta \
        -I /data/bam/reads_mother.bam \
        -O reads_mother.g.vcf \
        -L /data/ref/intervals.bed \
        -ERC GVCF
```

??? success "Saída do comando"

    ```console hl_lines="39 53 58 59"
    Using GATK jar /opt/conda/share/gatk4-4.5.0.0-0/gatk-package-4.5.0.0-local.jar
    Running:
        java -Dsamjdk.use_async_io_read_samtools=false -Dsamjdk.use_async_io_write_samtools=true -Dsamjdk.use_async_io_write_tribble=false -Dsamjdk.compression_level=2 -jar /opt/conda/share/gatk4-4.5.0.0-0/gatk-package-4.5.0.0-local.jar HaplotypeCaller -R /data/ref/ref.fasta -I /data/bam/reads_mother.bam -O reads_mother.g.vcf -L /data/ref/intervals.bed -ERC GVCF
    00:28:03.593 INFO  NativeLibraryLoader - Loading libgkl_compression.so from jar:file:/opt/conda/share/gatk4-4.5.0.0-0/gatk-package-4.5.0.0-local.jar!/com/intel/gkl/native/libgkl_compression.so
    00:28:03.765 INFO  HaplotypeCaller - ------------------------------------------------------------
    00:28:03.768 INFO  HaplotypeCaller - The Genome Analysis Toolkit (GATK) v4.5.0.0
    00:28:03.768 INFO  HaplotypeCaller - For support and documentation go to https://software.broadinstitute.org/gatk/
    00:28:03.768 INFO  HaplotypeCaller - Executing as root@8515e5a0598e on Linux v6.10.14-linuxkit amd64
    00:28:03.768 INFO  HaplotypeCaller - Java runtime: OpenJDK 64-Bit Server VM v17.0.11-internal+0-adhoc..src
    00:28:03.769 INFO  HaplotypeCaller - Start Date/Time: February 8, 2026 at 12:28:03 AM GMT
    00:28:03.769 INFO  HaplotypeCaller - ------------------------------------------------------------
    00:28:03.770 INFO  HaplotypeCaller - ------------------------------------------------------------
    00:28:03.772 INFO  HaplotypeCaller - HTSJDK Version: 4.1.0
    00:28:03.773 INFO  HaplotypeCaller - Picard Version: 3.1.1
    00:28:03.773 INFO  HaplotypeCaller - Built for Spark Version: 3.5.0
    00:28:03.773 INFO  HaplotypeCaller - HTSJDK Defaults.COMPRESSION_LEVEL : 2
    00:28:03.773 INFO  HaplotypeCaller - HTSJDK Defaults.USE_ASYNC_IO_READ_FOR_SAMTOOLS : false
    00:28:03.773 INFO  HaplotypeCaller - HTSJDK Defaults.USE_ASYNC_IO_WRITE_FOR_SAMTOOLS : true
    00:28:03.774 INFO  HaplotypeCaller - HTSJDK Defaults.USE_ASYNC_IO_WRITE_FOR_TRIBBLE : false
    00:28:03.774 INFO  HaplotypeCaller - Deflater: IntelDeflater
    00:28:03.774 INFO  HaplotypeCaller - Inflater: IntelInflater
    00:28:03.775 INFO  HaplotypeCaller - GCS max retries/reopens: 20
    00:28:03.775 INFO  HaplotypeCaller - Requester pays: disabled
    00:28:03.776 INFO  HaplotypeCaller - Initializing engine
    00:28:03.896 INFO  FeatureManager - Using codec BEDCodec to read file file:///data/ref/intervals.bed
    00:28:03.919 INFO  IntervalArgumentCollection - Processing 6369 bp from intervals
    00:28:03.934 INFO  HaplotypeCaller - Done initializing engine
    00:28:03.935 INFO  HaplotypeCallerEngine - Tool is in reference confidence mode and the annotation, the following changes will be made to any specified annotations: 'StrandBiasBySample' will be enabled. 'ChromosomeCounts', 'FisherStrand', 'StrandOddsRatio' and 'QualByDepth' annotations have been disabled
    00:28:03.943 INFO  NativeLibraryLoader - Loading libgkl_utils.so from jar:file:/opt/conda/share/gatk4-4.5.0.0-0/gatk-package-4.5.0.0-local.jar!/com/intel/gkl/native/libgkl_utils.so
    00:28:03.945 INFO  NativeLibraryLoader - Loading libgkl_smithwaterman.so from jar:file:/opt/conda/share/gatk4-4.5.0.0-0/gatk-package-4.5.0.0-local.jar!/com/intel/gkl/native/libgkl_smithwaterman.so
    00:28:03.946 INFO  SmithWatermanAligner - Using AVX accelerated SmithWaterman implementation
    00:28:03.955 INFO  HaplotypeCallerEngine - Standard Emitting and Calling confidence set to -0.0 for reference-model confidence output
    00:28:03.956 INFO  HaplotypeCallerEngine - All sites annotated with PLs forced to true for reference-model confidence output
    00:28:03.972 INFO  NativeLibraryLoader - Loading libgkl_pairhmm_omp.so from jar:file:/opt/conda/share/gatk4-4.5.0.0-0/gatk-package-4.5.0.0-local.jar!/com/intel/gkl/native/libgkl_pairhmm_omp.so
    00:28:03.993 INFO  IntelPairHmm - Flush-to-zero (FTZ) is enabled when running PairHMM
    00:28:03.994 INFO  IntelPairHmm - Available threads: 10
    00:28:03.994 INFO  IntelPairHmm - Requested threads: 4
    00:28:03.994 INFO  PairHMM - Using the OpenMP multi-threaded AVX-accelerated native PairHMM implementation
    00:28:04.044 INFO  ProgressMeter - Starting traversal
    00:28:04.070 INFO  ProgressMeter -        Current Locus  Elapsed Minutes     Regions Processed   Regions/Minute
    00:28:04.874 WARN  InbreedingCoeff - InbreedingCoeff will not be calculated at position 20_10037292_10066351:3480 and possibly subsequent; at least 10 samples must have called genotypes
    00:28:06.535 INFO  HaplotypeCaller - 7 read(s) filtered by: MappingQualityReadFilter
    0 read(s) filtered by: MappingQualityAvailableReadFilter
    0 read(s) filtered by: MappedReadFilter
    0 read(s) filtered by: NotSecondaryAlignmentReadFilter
    0 read(s) filtered by: NotDuplicateReadFilter
    0 read(s) filtered by: PassesVendorQualityCheckReadFilter
    0 read(s) filtered by: NonZeroReferenceLengthAlignmentReadFilter
    0 read(s) filtered by: GoodCigarReadFilter
    0 read(s) filtered by: WellformedReadFilter
    7 total reads filtered out of 1867 reads processed
    00:28:06.537 INFO  ProgressMeter - 20_10037292_10066351:13499              0.0                    35            851.6
    00:28:06.538 INFO  ProgressMeter - Traversal complete. Processed 35 total regions in 0.0 minutes.
    00:28:06.543 INFO  VectorLoglessPairHMM - Time spent in setup for JNI call : 0.003648749
    00:28:06.544 INFO  PairHMM - Total compute time in PairHMM computeLogLikelihoods() : 0.031498916
    00:28:06.544 INFO  SmithWatermanAligner - Total compute time in native Smith-Waterman : 0.02 sec
    00:28:06.547 INFO  HaplotypeCaller - Shutting down engine
    [February 8, 2026 at 12:28:06 AM GMT] org.broadinstitute.hellbender.tools.walkers.haplotypecaller.HaplotypeCaller done. Elapsed time: 0.05 minutes.
    Runtime.totalMemory()=281018368
    ```

Isso cria o arquivo de saída GVCF `reads_mother.g.vcf` no diretório de trabalho atual no contêiner.

Se você usar `cat` para visualizar o conteúdo, verá que é muito mais longo que o VCF equivalente que geramos na seção 1. Você nem consegue rolar até o início do arquivo, e a maioria das linhas parece bem diferente do que vimos no VCF.

??? abstract "Conteúdo do arquivo"

    ```console title="reads_mother.g.vcf" linenums="1674"
    20_10037292_10066351    14714   .       T       <NON_REF>       .       .       END=14718       GT:DP:GQ:MIN_DP:PL       0/0:37:99:37:0,99,1192
    20_10037292_10066351    14719   .       T       <NON_REF>       .       .       END=14719       GT:DP:GQ:MIN_DP:PL       0/0:36:82:36:0,82,1087
    20_10037292_10066351    14720   .       T       <NON_REF>       .       .       END=14737       GT:DP:GQ:MIN_DP:PL       0/0:42:99:37:0,100,1160
    ```

Estas representam regiões não-variantes onde o chamador de variantes não encontrou evidência de variação, então ele capturou algumas estatísticas descrevendo seu nível de confiança na ausência de variação.
Isso torna possível distinguir entre dois casos muito diferentes: (1) há dados de boa qualidade mostrando que a amostra é homozigota-referência, e (2) não há dados bons suficientes disponíveis para fazer uma determinação de qualquer forma.

Em um GVCF, normalmente há muitas dessas linhas não-variantes, com um número menor de registros de variantes espalhados entre elas.
Tente executar `head -176` no GVCF para carregar apenas as primeiras 176 linhas do arquivo para encontrar uma chamada de variante real.

??? abstract "Conteúdo do arquivo"

    ```console title="reads_mother.g.vcf" linenums="174"
    20_10037292_10066351    3479    .       T       <NON_REF>       .       .       END=3479        GT:DP:GQ:MIN_DP:PL       0/0:34:36:34:0,36,906
    20_10037292_10066351    3480    .       C       CT,<NON_REF>    503.03  .       DP=23;ExcessHet=0.0000;MLEAC=2,0;MLEAF=1.00,0.00;RAW_MQandDP=82800,23    GT:AD:DP:GQ:PL:SB       1/1:0,18,0:18:54:517,54,0,517,54,517:0,0,7,11
    20_10037292_10066351    3481    .       T       <NON_REF>       .       .       END=3481        GT:DP:GQ:MIN_DP:PL       0/0:21:51:21:0,51,765
    ```

A segunda linha mostra o primeiro registro de variante no arquivo, que corresponde à primeira variante no arquivo VCF que vimos anteriormente.

Assim como o VCF original, o arquivo de saída GVCF também é acompanhado por um arquivo de índice, chamado `reads_mother.g.vcf.idx`.

#### 2.2.3. Repetir o processo nas outras duas amostras

Gere GVCFs para as duas amostras restantes executando os comandos abaixo, um após o outro.

```bash
gatk HaplotypeCaller \
        -R /data/ref/ref.fasta \
        -I /data/bam/reads_father.bam \
        -O reads_father.g.vcf \
        -L /data/ref/intervals.bed \
        -ERC GVCF
```

```bash
gatk HaplotypeCaller \
        -R /data/ref/ref.fasta \
        -I /data/bam/reads_son.bam \
        -O reads_son.g.vcf \
        -L /data/ref/intervals.bed \
        -ERC GVCF
```

Uma vez concluído, você deve ter três arquivos terminando em `.g.vcf` no seu diretório atual (um por amostra) e seus respectivos arquivos de índice terminando em `.g.vcf.idx`.

Mas não saia do contêiner!
Vamos usar o mesmo contêiner na próxima etapa.

### 2.3. Executar genotipagem conjunta

Agora que temos todos os GVCFs, podemos testar a abordagem de genotipagem conjunta para gerar chamadas de variantes para uma coorte de amostras.
É um método de duas etapas que consiste em combinar os dados de todos os GVCFs em um armazenamento de dados, e então executar a análise de genotipagem conjunta propriamente dita para gerar o VCF final de variantes chamadas conjuntamente.

#### 2.3.1. Combinar todos os GVCFs por amostra

Esta primeira etapa usa outra ferramenta GATK, chamada GenomicsDBImport, para combinar os dados de todos os GVCFs em um armazenamento de dados GenomicsDB.

```bash
gatk GenomicsDBImport \
    -V reads_mother.g.vcf \
    -V reads_father.g.vcf \
    -V reads_son.g.vcf \
    -L /data/ref/intervals.bed \
    --genomicsdb-workspace-path family_trio_gdb
```

??? success "Saída do comando"

    ```console hl_lines="33 36 37 39 40"
    Using GATK jar /opt/conda/share/gatk4-4.5.0.0-0/gatk-package-4.5.0.0-local.jar
    Running:
        java -Dsamjdk.use_async_io_read_samtools=false -Dsamjdk.use_async_io_write_samtools=true -Dsamjdk.use_async_io_write_tribble=false -Dsamjdk.compression_level=2 -jar /opt/conda/share/gatk4-4.5.0.0-0/gatk-package-4.5.0.0-local.jar GenomicsDBImport -V reads_mother.g.vcf -V reads_father.g.vcf -V reads_son.g.vcf -L /data/ref/intervals.bed --genomicsdb-workspace-path family_trio_gdb
    00:28:20.772 INFO  NativeLibraryLoader - Loading libgkl_compression.so from jar:file:/opt/conda/share/gatk4-4.5.0.0-0/gatk-package-4.5.0.0-local.jar!/com/intel/gkl/native/libgkl_compression.so
    00:28:20.914 INFO  GenomicsDBImport - ------------------------------------------------------------
    00:28:20.917 INFO  GenomicsDBImport - The Genome Analysis Toolkit (GATK) v4.5.0.0
    00:28:20.917 INFO  GenomicsDBImport - For support and documentation go to https://software.broadinstitute.org/gatk/
    00:28:20.917 INFO  GenomicsDBImport - Executing as root@8515e5a0598e on Linux v6.10.14-linuxkit amd64
    00:28:20.917 INFO  GenomicsDBImport - Java runtime: OpenJDK 64-Bit Server VM v17.0.11-internal+0-adhoc..src
    00:28:20.918 INFO  GenomicsDBImport - Start Date/Time: February 8, 2026 at 12:28:20 AM GMT
    00:28:20.918 INFO  GenomicsDBImport - ------------------------------------------------------------
    00:28:20.918 INFO  GenomicsDBImport - ------------------------------------------------------------
    00:28:20.920 INFO  GenomicsDBImport - HTSJDK Version: 4.1.0
    00:28:20.921 INFO  GenomicsDBImport - Picard Version: 3.1.1
    00:28:20.921 INFO  GenomicsDBImport - Built for Spark Version: 3.5.0
    00:28:20.922 INFO  GenomicsDBImport - HTSJDK Defaults.COMPRESSION_LEVEL : 2
    00:28:20.923 INFO  GenomicsDBImport - HTSJDK Defaults.USE_ASYNC_IO_READ_FOR_SAMTOOLS : false
    00:28:20.923 INFO  GenomicsDBImport - HTSJDK Defaults.USE_ASYNC_IO_WRITE_FOR_SAMTOOLS : true
    00:28:20.923 INFO  GenomicsDBImport - HTSJDK Defaults.USE_ASYNC_IO_WRITE_FOR_TRIBBLE : false
    00:28:20.923 INFO  GenomicsDBImport - Deflater: IntelDeflater
    00:28:20.924 INFO  GenomicsDBImport - Inflater: IntelInflater
    00:28:20.924 INFO  GenomicsDBImport - GCS max retries/reopens: 20
    00:28:20.924 INFO  GenomicsDBImport - Requester pays: disabled
    00:28:20.925 INFO  GenomicsDBImport - Initializing engine
    00:28:21.144 INFO  FeatureManager - Using codec BEDCodec to read file file:///data/ref/intervals.bed
    00:28:21.152 INFO  IntervalArgumentCollection - Processing 6369 bp from intervals
    00:28:21.157 INFO  GenomicsDBImport - Done initializing engine
    00:28:21.287 INFO  GenomicsDBLibLoader - GenomicsDB native library version : 1.5.1-84e800e
    00:28:21.290 INFO  GenomicsDBImport - Vid Map JSON file will be written to /tmp/family_trio_gdb/vidmap.json
    00:28:21.290 INFO  GenomicsDBImport - Callset Map JSON file will be written to /tmp/family_trio_gdb/callset.json
    00:28:21.291 INFO  GenomicsDBImport - Complete VCF Header will be written to /tmp/family_trio_gdb/vcfheader.vcf
    00:28:21.291 INFO  GenomicsDBImport - Importing to workspace - /tmp/family_trio_gdb
    00:28:21.453 INFO  GenomicsDBImport - Importing batch 1 with 3 samples
    00:28:21.757 INFO  GenomicsDBImport - Importing batch 1 with 3 samples
    00:28:21.859 INFO  GenomicsDBImport - Importing batch 1 with 3 samples
    00:28:21.979 INFO  GenomicsDBImport - Done importing batch 1/1
    00:28:21.988 INFO  GenomicsDBImport - Import completed!
    00:28:21.988 INFO  GenomicsDBImport - Shutting down engine
    [February 8, 2026 at 12:28:21 AM GMT] org.broadinstitute.hellbender.tools.genomicsdb.GenomicsDBImport done. Elapsed time: 0.02 minutes.
    Runtime.totalMemory()=305135616
    ```

A saída desta etapa é efetivamente um diretório contendo um conjunto de diretórios ainda mais aninhados contendo os dados de variantes combinados na forma de múltiplos arquivos diferentes.
Você pode explorar, mas rapidamente verá que este formato de armazenamento de dados não é destinado a ser lido diretamente por humanos.

!!! note "Nota"

    O GATK inclui ferramentas que tornam possível inspecionar e extrair dados de chamada de variantes do armazenamento de dados conforme necessário.

#### 2.3.2. Executar a análise de genotipagem conjunta propriamente dita

Esta segunda etapa usa ainda outra ferramenta GATK, chamada GenotypeGVCFs, para recalcular estatísticas de variantes e genótipos individuais à luz dos dados disponíveis em todas as amostras da coorte.

```bash
gatk GenotypeGVCFs \
    -R /data/ref/ref.fasta \
    -V gendb://family_trio_gdb \
    -O family_trio.vcf
```

??? success "Saída do comando"

    ```console hl_lines="30 35 37 38"
    Using GATK jar /opt/conda/share/gatk4-4.5.0.0-0/gatk-package-4.5.0.0-local.jar
    Running:
        java -Dsamjdk.use_async_io_read_samtools=false -Dsamjdk.use_async_io_write_samtools=true -Dsamjdk.use_async_io_write_tribble=false -Dsamjdk.compression_level=2 -jar /opt/conda/share/gatk4-4.5.0.0-0/gatk-package-4.5.0.0-local.jar GenotypeGVCFs -R /data/ref/ref.fasta -V gendb://family_trio_gdb -O family_trio.vcf
    00:28:24.625 INFO  NativeLibraryLoader - Loading libgkl_compression.so from jar:file:/opt/conda/share/gatk4-4.5.0.0-0/gatk-package-4.5.0.0-local.jar!/com/intel/gkl/native/libgkl_compression.so
    00:28:24.798 INFO  GenotypeGVCFs - ------------------------------------------------------------
    00:28:24.801 INFO  GenotypeGVCFs - The Genome Analysis Toolkit (GATK) v4.5.0.0
    00:28:24.801 INFO  GenotypeGVCFs - For support and documentation go to https://software.broadinstitute.org/gatk/
    00:28:24.801 INFO  GenotypeGVCFs - Executing as root@8515e5a0598e on Linux v6.10.14-linuxkit amd64
    00:28:24.801 INFO  GenotypeGVCFs - Java runtime: OpenJDK 64-Bit Server VM v17.0.11-internal+0-adhoc..src
    00:28:24.802 INFO  GenotypeGVCFs - Start Date/Time: February 8, 2026 at 12:28:24 AM GMT
    00:28:24.802 INFO  GenotypeGVCFs - ------------------------------------------------------------
    00:28:24.802 INFO  GenotypeGVCFs - ------------------------------------------------------------
    00:28:24.804 INFO  GenotypeGVCFs - HTSJDK Version: 4.1.0
    00:28:24.804 INFO  GenotypeGVCFs - Picard Version: 3.1.1
    00:28:24.804 INFO  GenotypeGVCFs - Built for Spark Version: 3.5.0
    00:28:24.805 INFO  GenotypeGVCFs - HTSJDK Defaults.COMPRESSION_LEVEL : 2
    00:28:24.805 INFO  GenotypeGVCFs - HTSJDK Defaults.USE_ASYNC_IO_READ_FOR_SAMTOOLS : false
    00:28:24.806 INFO  GenotypeGVCFs - HTSJDK Defaults.USE_ASYNC_IO_WRITE_FOR_SAMTOOLS : true
    00:28:24.806 INFO  GenotypeGVCFs - HTSJDK Defaults.USE_ASYNC_IO_WRITE_FOR_TRIBBLE : false
    00:28:24.806 INFO  GenotypeGVCFs - Deflater: IntelDeflater
    00:28:24.806 INFO  GenotypeGVCFs - Inflater: IntelInflater
    00:28:24.807 INFO  GenotypeGVCFs - GCS max retries/reopens: 20
    00:28:24.807 INFO  GenotypeGVCFs - Requester pays: disabled
    00:28:24.808 INFO  GenotypeGVCFs - Initializing engine
    00:28:25.023 INFO  GenomicsDBLibLoader - GenomicsDB native library version : 1.5.1-84e800e
    00:28:25.081 INFO  NativeGenomicsDB - pid=162 tid=163 No valid combination operation found for INFO field InbreedingCoeff  - the field will NOT be part of INFO fields in the generated VCF records
    00:28:25.082 INFO  NativeGenomicsDB - pid=162 tid=163 No valid combination operation found for INFO field MLEAC  - the field will NOT be part of INFO fields in the generated VCF records
    00:28:25.082 INFO  NativeGenomicsDB - pid=162 tid=163 No valid combination operation found for INFO field MLEAF  - the field will NOT be part of INFO fields in the generated VCF records
    00:28:25.109 INFO  GenotypeGVCFs - Done initializing engine
    00:28:25.184 INFO  ProgressMeter - Starting traversal
    00:28:25.187 INFO  ProgressMeter -        Current Locus  Elapsed Minutes    Variants Processed  Variants/Minute
    00:28:25.446 WARN  InbreedingCoeff - InbreedingCoeff will not be calculated at position 20_10037292_10066351:3480 and possibly subsequent; at least 10 samples must have called genotypes
    GENOMICSDB_TIMER,GenomicsDB iterator next() timer,Wall-clock time(s),0.15034835899999904,Cpu time(s),0.1355218420000006
    00:28:26.189 INFO  ProgressMeter - 20_10037292_10066351:13953              0.0                  3390         202994.0
    00:28:26.190 INFO  ProgressMeter - Traversal complete. Processed 3390 total variants in 0.0 minutes.
    00:28:26.194 INFO  GenotypeGVCFs - Shutting down engine
    [February 8, 2026 at 12:28:26 AM GMT] org.broadinstitute.hellbender.tools.walkers.GenotypeGVCFs done. Elapsed time: 0.03 minutes.
    Runtime.totalMemory()=296747008
    ```

Isso cria o arquivo de saída VCF `family_trio.vcf` no diretório de trabalho atual no contêiner.
É outro arquivo razoavelmente pequeno, então você pode usar `cat` neste arquivo para visualizar seu conteúdo, e rolar para cima para encontrar as primeiras linhas de variante.

??? abstract "Conteúdo do arquivo"

    ```console title="family_trio.vcf" linenums="40"
    #CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO    FORMAT  reads_father    reads_mother    reads_son
    20_10037292_10066351    3480    .       C       CT      1625.89 .       AC=5;AF=0.833;AN=6;BaseQRankSum=0.220;DP=85;ExcessHet=0.0000;FS=2.476;MLEAC=5;MLEAF=0.833;MQ=60.00;MQRankSum=0.00;QD=21.68;ReadPosRankSum=-1.147e+00;SOR=0.487    GT:AD:DP:GQ:PL  0/1:15,16:31:99:367,0,375       1/1:0,18:18:54:517,54,0 1/1:0,26:26:78:756,78,0
    20_10037292_10066351    3520    .       AT      A       1678.89 .       AC=5;AF=0.833;AN=6;BaseQRankSum=1.03;DP=80;ExcessHet=0.0000;FS=2.290;MLEAC=5;MLEAF=0.833;MQ=60.00;MQRankSum=0.00;QD=22.39;ReadPosRankSum=0.701;SOR=0.730 GT:AD:DP:GQ:PL   0/1:18,13:31:99:296,0,424       1/1:0,18:18:54:623,54,0 1/1:0,26:26:78:774,78,0
    20_10037292_10066351    3529    .       T       A       154.29  .       AC=1;AF=0.167;AN=6;BaseQRankSum=-5.440e-01;DP=104;ExcessHet=0.0000;FS=1.871;MLEAC=1;MLEAF=0.167;MQ=60.00;MQRankSum=0.00;QD=7.71;ReadPosRankSum=-1.158e+00;SOR=1.034       GT:AD:DP:GQ:PL  0/0:44,0:44:99:0,112,1347       0/1:12,8:20:99:163,0,328        0/0:39,0:39:99:0,105,1194
    ```

Isso parece similar ao VCF que geramos anteriormente, exceto que desta vez temos informações de nível de genótipo para todas as três amostras.
As últimas três colunas no arquivo são os blocos de genótipo para as amostras, listadas em ordem alfabética.

Se olharmos para os genótipos chamados para nosso trio familiar de teste para a primeira variante, vemos que o pai é heterozigoto-variante (`0/1`), e a mãe e o filho são ambos homozigoto-variante (`1/1`).

Essa é, em última análise, a informação que estamos buscando extrair do conjunto de dados!

#### 2.3.3. Sair do contêiner GATK

Para sair do contêiner, digite `exit`.

```bash
exit
```

Seu prompt deve estar de volta ao normal.
Isso conclui o teste manual dos comandos de chamada de variantes.

---

### Conclusão

Você sabe como testar os comandos de indexação Samtools e chamada de variantes GATK em seus respectivos contêineres, incluindo como gerar GVCFs e executar genotipagem conjunta em múltiplas amostras.

### O que vem a seguir?

Aprenda como envolver esses mesmos comandos em fluxos de trabalho que usam contêineres para executar o trabalho.
