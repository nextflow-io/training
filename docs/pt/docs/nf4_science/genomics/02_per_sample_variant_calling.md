# Parte 2: Chamada de variantes por amostra

Na Parte 1, você testou os comandos do Samtools e do GATK manualmente em seus respectivos contêineres.
Agora vamos envolver esses mesmos comandos em um fluxo de trabalho Nextflow.

## Tarefa

Nesta parte do curso, vamos desenvolver um fluxo de trabalho que faz o seguinte:

1. Gerar um arquivo de índice para cada arquivo BAM de entrada usando [Samtools](https://www.htslib.org/)
2. Executar o GATK HaplotypeCaller em cada arquivo BAM de entrada para gerar chamadas de variantes por amostra em VCF (Variant Call Format)

<figure class="excalidraw">
--8<-- "docs/en/docs/nf4_science/genomics/img/hello-gatk-1.svg"
</figure>

Isso replica as etapas da Parte 1, onde você executou esses comandos manualmente em seus contêineres.

Como ponto de partida, fornecemos um arquivo de fluxo de trabalho, `genomics.nf`, que descreve as principais partes do fluxo de trabalho, além de dois arquivos de módulo, samtools_index.nf e gatk_haplotypecaller.nf, que descrevem a estrutura dos módulos.
Esses arquivos não são funcionais; seu propósito é apenas servir como estruturas para você preencher com as partes interessantes do código.

## Plano de aula

Para tornar o processo de desenvolvimento mais educativo, dividimos isso em quatro etapas:

1. **Escrever um fluxo de trabalho de único estágio que executa o Samtools index em um arquivo BAM.**
   Isso cobre a criação de um módulo, importá-lo e chamá-lo em um fluxo de trabalho.
2. **Adicionar um segundo processo para executar o GATK HaplotypeCaller no arquivo BAM indexado.**
   Isso introduz o encadeamento de saídas de processos para entradas e o tratamento de arquivos acessórios.
3. **Adaptar o fluxo de trabalho para executar em um lote de amostras.**
   Isso cobre a execução paralela e introduz tuplas para manter arquivos associados juntos.
4. **Fazer o fluxo de trabalho aceitar um arquivo de texto contendo um lote de arquivos de entrada.**
   Isso demonstra um padrão comum para fornecer entradas em massa.

Cada etapa se concentra em um aspecto específico do desenvolvimento do fluxo de trabalho.

---

## 1. Escrever um fluxo de trabalho de único estágio que executa o Samtools index em um arquivo BAM

Esta primeira etapa se concentra no básico: carregar um arquivo BAM e gerar um índice para ele.

Lembre-se do comando `samtools index` da [Parte 1](01_method.md):

```bash
samtools index '<input_bam>'
```

O comando recebe um arquivo BAM como entrada e produz um arquivo de índice `.bai` ao lado dele.
A URI do contêiner era `community.wave.seqera.io/library/samtools:1.20--b5dfbd93de237464`.

Vamos pegar essas informações e envolvê-las em Nextflow em três estágios:

1. Configurar a entrada
2. Escrever o processo de indexação e chamá-lo no fluxo de trabalho
3. Configurar o tratamento da saída

### 1.1. Configurar a entrada

Precisamos declarar um parâmetro de entrada, criar um perfil de teste para fornecer um valor padrão conveniente e criar um canal de entrada.

#### 1.1.1. Adicionar uma declaração de parâmetro de entrada

No arquivo principal do fluxo de trabalho `genomics.nf`, na seção `Pipeline parameters`, declare um parâmetro CLI chamado `reads_bam`.

=== "Depois"

    ```groovy title="genomics.nf" linenums="5" hl_lines="4-7"
    /*
     * Pipeline parameters
     */
    params {
        // Primary input
        reads_bam: Path
    }
    ```

=== "Antes"

    ```groovy title="genomics.nf" linenums="5"
    /*
     * Pipeline parameters
     */

    // Primary input
    ```

Isso configura o parâmetro CLI, mas não queremos digitar o caminho do arquivo toda vez que executamos o fluxo de trabalho durante o desenvolvimento.
Existem várias opções para fornecer um valor padrão; aqui usamos um perfil de teste.

#### 1.1.2. Criar um perfil de teste com um valor padrão em `nextflow.config`

Um perfil de teste fornece valores padrão convenientes para experimentar um fluxo de trabalho sem especificar entradas na linha de comando.
Esta é uma convenção comum no ecossistema Nextflow (veja [Hello Config](../../hello_nextflow/06_hello_config.md) para mais detalhes).

Adicione um bloco `profiles` ao `nextflow.config` com um perfil `test` que define o parâmetro `reads_bam` para um dos arquivos BAM de teste.

=== "Depois"

    ```groovy title="nextflow.config" linenums="1" hl_lines="3-7"
    docker.enabled = true

    profiles {
        test {
            params.reads_bam = "${projectDir}/data/bam/reads_mother.bam"
        }
    }
    ```

=== "Antes"

    ```groovy title="nextflow.config" linenums="1"
    docker.enabled = true
    ```

Aqui, estamos usando `${projectDir}`, uma variável integrada do Nextflow que aponta para o diretório onde o script do fluxo de trabalho está localizado.
Isso facilita a referência a arquivos de dados e outros recursos sem codificar caminhos absolutos.

#### 1.1.3. Configurar o canal de entrada

No bloco workflow, crie um canal de entrada a partir do valor do parâmetro usando a factory de canal `.fromPath` (como usado em [Hello Channels](../../hello_nextflow/02_hello_channels.md)).

=== "Depois"

    ```groovy title="genomics.nf" linenums="13" hl_lines="4-5"
    workflow {

        main:
        // Create input channel (single file via CLI parameter)
        reads_ch = channel.fromPath(params.reads_bam)
    ```

=== "Antes"

    ```groovy title="genomics.nf" linenums="13"
    workflow {

        main:
        // Create input channel
    ```

Agora precisamos criar o processo para executar a indexação nesta entrada.

### 1.2. Escrever o processo de indexação e chamá-lo no fluxo de trabalho

Precisamos escrever a definição do processo no arquivo do módulo, importá-lo para o fluxo de trabalho usando uma instrução include e chamá-lo na entrada.

#### 1.2.1. Preencher o módulo para o processo de indexação

Abra `modules/samtools_index.nf` e examine a estrutura da definição do processo.
Você deve reconhecer os principais elementos estruturais; caso contrário, considere ler [Hello Nextflow](../../hello_nextflow/01_hello_world.md) para relembrar.

Vá em frente e preencha a definição do processo por conta própria usando as informações fornecidas acima, depois verifique seu trabalho com a solução na aba "Depois" abaixo.

=== "Antes"

    ```groovy title="modules/samtools_index.nf" linenums="1"
    #!/usr/bin/env nextflow

    /*
     * Generate BAM index file
     */
    process SAMTOOLS_INDEX {

        container

        input:

        output:

        script:
        """

        """
    }
    ```

=== "Depois"

    ```groovy title="modules/samtools_index.nf" linenums="1" hl_lines="8 11 14 18"
    #!/usr/bin/env nextflow

    /*
     * Gera arquivo de índice BAM
     */
    process SAMTOOLS_INDEX {

        container 'community.wave.seqera.io/library/samtools:1.20--b5dfbd93de237464'

        input:
        path input_bam

        output:
        path "${input_bam}.bai"

        script:
        """
        samtools index '$input_bam'
        """
    }
    ```

Depois de concluir isso, o processo está completo.
Para usá-lo no fluxo de trabalho, você precisará importar o módulo e adicionar uma chamada de processo.

#### 1.2.2. Incluir o módulo

Em `genomics.nf`, adicione uma instrução `include` para tornar o processo disponível para o fluxo de trabalho:

=== "Depois"

    ```groovy title="genomics.nf" linenums="3" hl_lines="2"
    // Module INCLUDE statements
    include { SAMTOOLS_INDEX } from './modules/samtools_index.nf'
    ```

=== "Antes"

    ```groovy title="genomics.nf" linenums="3"
    // Module INCLUDE statements
    ```

O processo agora está disponível no escopo do fluxo de trabalho.

#### 1.2.3. Chamar o processo de indexação na entrada

Agora, vamos adicionar uma chamada para `SAMTOOLS_INDEX` no bloco workflow, passando o canal de entrada como argumento.

=== "Depois"

    ```groovy title="genomics.nf" linenums="14" hl_lines="7-8"
    workflow {

        main:
        // Create input channel (single file via CLI parameter)
        reads_ch = channel.fromPath(params.reads_bam)

        // Create index file for input BAM file
        SAMTOOLS_INDEX(reads_ch)
    ```

=== "Antes"

    ```groovy title="genomics.nf" linenums="14"
    workflow {

        main:
        // Create input channel (single file via CLI parameter)
        reads_ch = channel.fromPath(params.reads_bam)

        // Call processes
    ```

O fluxo de trabalho agora carrega a entrada e executa o processo de indexação nela.
Em seguida, precisamos configurar como a saída é publicada.

### 1.3. Configurar o tratamento da saída

Precisamos declarar quais saídas do processo publicar e especificar para onde elas devem ir.

#### 1.3.1. Declarar uma saída na seção `publish:`

A seção `publish:` dentro do bloco workflow declara quais saídas do processo devem ser publicadas.
Atribua a saída de `SAMTOOLS_INDEX` a um alvo nomeado chamado `bam_index`.

=== "Depois"

    ```groovy title="genomics.nf" linenums="22" hl_lines="2"
        publish:
        bam_index = SAMTOOLS_INDEX.out
    }
    ```

=== "Antes"

    ```groovy title="genomics.nf" linenums="22"
        publish:
        // Declare outputs to publish
    }
    ```

Agora precisamos dizer ao Nextflow onde colocar a saída publicada.

#### 1.3.2. Configurar o alvo de saída no bloco `output {}`

O bloco `output {}` fica fora do fluxo de trabalho e especifica onde cada alvo nomeado é publicado.
Vamos adicionar um alvo para `bam_index` que publica em um subdiretório `bam/`.

=== "Depois"

    ```groovy title="genomics.nf" linenums="26" hl_lines="2-4"
    output {
        bam_index {
            path 'bam'
        }
    }
    ```

=== "Antes"

    ```groovy title="genomics.nf" linenums="26"
    output {
        // Configure publish targets
    }
    ```

!!! note

    Por padrão, o Nextflow publica arquivos de saída como links simbólicos, o que evita duplicação desnecessária.
    Embora os arquivos de dados que estamos usando aqui sejam muito pequenos, em genômica eles podem ficar muito grandes.
    Os links simbólicos serão quebrados quando você limpar seu diretório `work`, então para fluxos de trabalho de produção você pode querer substituir o modo de publicação padrão para `'copy'`.

### 1.4. Executar o fluxo de trabalho

Neste ponto, temos um fluxo de trabalho de indexação de uma etapa que deve ser totalmente funcional. Vamos testar se funciona!

Podemos executá-lo com `-profile test` para usar o valor padrão configurado no perfil de teste e evitar ter que escrever o caminho na linha de comando.

```bash
nextflow run genomics.nf -profile test
```

??? success "Saída do comando"

    ```console
    N E X T F L O W   ~  version 25.10.2

    ┃ Launching `genomics.nf` [reverent_sinoussi] DSL2 - revision: 41d43ad7fe

    executor >  local (1)
    [2a/e69536] SAMTOOLS_INDEX (1) | 1 of 1 ✔
    ```

Você pode verificar se o arquivo de índice foi gerado corretamente olhando no diretório de trabalho ou no diretório de resultados.

??? abstract "Conteúdo do diretório de trabalho"

    ```console
    work/2a/e695367b2f60df09cf826b07192dc3
    ├── reads_mother.bam -> /workspaces/training/nf4-science/genomics/data/bam/reads_mother.bam
    └── reads_mother.bam.bai
    ```

??? abstract "Conteúdo do diretório de resultados"

    ```console
    results/
    └── bam/
        └── reads_mother.bam.bai -> ...
    ```

Aí está!

### Conclusão

Você sabe como criar um módulo contendo um processo, importá-lo para um fluxo de trabalho, chamá-lo com um canal de entrada e publicar os resultados.

### O que vem a seguir?

Adicionar uma segunda etapa que pega a saída do processo de indexação e a usa para executar a chamada de variantes.

---

## 2. Adicionar um segundo processo para executar o GATK HaplotypeCaller no arquivo BAM indexado

Agora que temos um índice para nosso arquivo de entrada, podemos passar para a configuração da etapa de chamada de variantes.

Lembre-se do comando `gatk HaplotypeCaller` da [Parte 1](01_method.md):

```bash
gatk HaplotypeCaller \
        -R /data/ref/ref.fasta \
        -I /data/bam/reads_mother.bam \
        -O reads_mother.vcf \
        -L /data/ref/intervals.bed
```

O comando recebe um arquivo BAM (`-I`), um genoma de referência (`-R`) e um arquivo de intervalos (`-L`), e produz um arquivo VCF (`-O`) junto com seu índice.
A ferramenta também espera que o índice do BAM, o índice da referência e o dicionário da referência estejam colocalizados com seus respectivos arquivos.
A URI do contêiner era `community.wave.seqera.io/library/gatk4:4.5.0.0--730ee8817e436867`.

Seguimos os mesmos três estágios de antes:

1. Configurar as entradas
2. Escrever o processo de chamada de variantes e chamá-lo no fluxo de trabalho
3. Configurar o tratamento da saída

### 2.1. Configurar as entradas

A etapa de chamada de variantes requer vários arquivos de entrada adicionais.
Precisamos declarar parâmetros para eles, adicionar valores padrão ao perfil de teste e criar variáveis para carregá-los.

#### 2.1.1. Adicionar declarações de parâmetros para entradas acessórias

Como nosso novo processo espera alguns arquivos adicionais para serem fornecidos, adicione declarações de parâmetros para eles em `genomics.nf` na seção `Pipeline parameters`:

=== "Depois"

    ```groovy title="genomics.nf" linenums="9" hl_lines="5-9"
    params {
        // Primary input
        reads_bam: Path

        // Accessory files
        reference: Path
        reference_index: Path
        reference_dict: Path
        intervals: Path
    }
    ```

=== "Antes"

    ```groovy title="genomics.nf" linenums="9"
    params {
        // Primary input
        reads_bam: Path
    }
    ```

Como antes, fornecemos valores padrão através do perfil de teste em vez de inline.

#### 2.1.2. Adicionar padrões de arquivos acessórios ao perfil de teste

Assim como fizemos para `reads_bam` na seção 1.1.2, adicione valores padrão para os arquivos acessórios ao perfil de teste em `nextflow.config`:

=== "Depois"

    ```groovy title="nextflow.config" linenums="4" hl_lines="3-6"
    test {
        params.reads_bam = "${projectDir}/data/bam/reads_mother.bam"
        params.reference = "${projectDir}/data/ref/ref.fasta"
        params.reference_index = "${projectDir}/data/ref/ref.fasta.fai"
        params.reference_dict = "${projectDir}/data/ref/ref.dict"
        params.intervals = "${projectDir}/data/ref/intervals.bed"
    }
    ```

=== "Antes"

    ```groovy title="nextflow.config" linenums="4"
    test {
        params.reads_bam = "${projectDir}/data/bam/reads_mother.bam"
    }
    ```

Agora precisamos criar variáveis que carreguem esses caminhos de arquivo para uso no fluxo de trabalho.

#### 2.1.3. Criar variáveis para os arquivos acessórios

Adicione variáveis para os caminhos dos arquivos acessórios dentro do bloco workflow:

=== "Depois"

    ```groovy title="genomics.nf" linenums="21" hl_lines="7-11"
    workflow {

        main:
        // Create input channel (single file via CLI parameter)
        reads_ch = channel.fromPath(params.reads_bam)

        // Load the file paths for the accessory files (reference and intervals)
        ref_file        = file(params.reference)
        ref_index_file  = file(params.reference_index)
        ref_dict_file   = file(params.reference_dict)
        intervals_file  = file(params.intervals)

        // Create index file for input BAM file
        SAMTOOLS_INDEX(reads_ch)
    ```

=== "Antes"

    ```groovy title="genomics.nf" linenums="21"
    workflow {

        main:
        // Create input channel (single file via CLI parameter)
        reads_ch = channel.fromPath(params.reads_bam)

        // Create index file for input BAM file
        SAMTOOLS_INDEX(reads_ch)
    ```

A sintaxe `file()` diz explicitamente ao Nextflow para tratar essas entradas como caminhos de arquivo.
Você pode aprender mais sobre isso na Side Quest [Trabalhando com arquivos](../../side_quests/working_with_files.md).

### 2.2. Escrever o processo de chamada de variantes e chamá-lo no fluxo de trabalho

Precisamos escrever a definição do processo no arquivo do módulo, importá-lo para o fluxo de trabalho usando uma instrução include e chamá-lo nas leituras de entrada mais a saída da etapa de indexação e os arquivos acessórios.

#### 2.2.1. Preencher o módulo para o processo de chamada de variantes

Abra `modules/gatk_haplotypecaller.nf` e examine a estrutura da definição do processo.

Vá em frente e preencha a definição do processo por conta própria usando as informações fornecidas acima, depois verifique seu trabalho com a solução na aba "Depois" abaixo.

=== "Antes"

    ```groovy title="modules/gatk_haplotypecaller.nf" linenums="1"
    #!/usr/bin/env nextflow

    /*
     * Call variants with GATK HaplotypeCaller
     */
    process GATK_HAPLOTYPECALLER {

        container

        input:

        output:

        script:
        """

        """
    }
    ```

=== "Depois"

    ```groovy title="modules/gatk_haplotypecaller.nf" linenums="1" hl_lines="8 11-16 19-20 24-28"
    #!/usr/bin/env nextflow

    /*
     * Chama variantes com GATK HaplotypeCaller
     */
    process GATK_HAPLOTYPECALLER {

        container "community.wave.seqera.io/library/gatk4:4.5.0.0--730ee8817e436867"

        input:
        path input_bam
        path input_bam_index
        path ref_fasta
        path ref_index
        path ref_dict
        path interval_list

        output:
        path "${input_bam}.vcf"     , emit: vcf
        path "${input_bam}.vcf.idx" , emit: idx

        script:
        """
        gatk HaplotypeCaller \
            -R ${ref_fasta} \
            -I ${input_bam} \
            -O ${input_bam}.vcf \
            -L ${interval_list}
        """
    }
    ```

Você notará que este processo tem mais entradas do que o comando GATK realmente requer.
O GATK sabe procurar o arquivo de índice do BAM e os arquivos acessórios do genoma de referência com base em convenções de nomenclatura, mas o Nextflow é agnóstico de domínio e não sabe sobre essas convenções.
Precisamos listá-los explicitamente para que o Nextflow os prepare no diretório de trabalho em tempo de execução; caso contrário, o GATK lançará um erro sobre arquivos ausentes.

Da mesma forma, listamos explicitamente o arquivo de índice da saída VCF (`"${input_bam}.vcf.idx"`) para que o Nextflow mantenha o controle dele para etapas subsequentes.
Usamos a sintaxe `emit:` para atribuir um nome a cada canal de saída, o que se tornará útil quando conectarmos as saídas ao bloco publish.

Depois de concluir isso, o processo está completo.
Para usá-lo no fluxo de trabalho, você precisará importar o módulo e adicionar uma chamada de processo.

#### 2.2.2. Importar o novo módulo

Atualize `genomics.nf` para importar o novo módulo:

=== "Depois"

    ```groovy title="genomics.nf" linenums="3" hl_lines="3"
    // Module INCLUDE statements
    include { SAMTOOLS_INDEX } from './modules/samtools_index.nf'
    include { GATK_HAPLOTYPECALLER } from './modules/gatk_haplotypecaller.nf'
    ```

=== "Antes"

    ```groovy title="genomics.nf" linenums="3"
    // Module INCLUDE statements
    include { SAMTOOLS_INDEX } from './modules/samtools_index.nf'
    ```

O processo agora está disponível no escopo do fluxo de trabalho.

#### 2.2.3. Adicionar a chamada do processo

Adicione a chamada do processo no corpo do fluxo de trabalho, sob `main:`:

=== "Depois"

    ```groovy title="genomics.nf" linenums="33" hl_lines="4-12"
        // Create index file for input BAM file
        SAMTOOLS_INDEX(reads_ch)

        // Call variants from the indexed BAM file
        GATK_HAPLOTYPECALLER(
            reads_ch,
            SAMTOOLS_INDEX.out,
            ref_file,
            ref_index_file,
            ref_dict_file,
            intervals_file
        )
    ```

=== "Antes"

    ```groovy title="genomics.nf" linenums="33"
        // Create index file for input BAM file
        SAMTOOLS_INDEX(reads_ch)
    ```

Você deve reconhecer a sintaxe `*.out` da série de treinamento Hello Nextflow; estamos dizendo ao Nextflow para pegar a saída do canal por `SAMTOOLS_INDEX` e conectá-la à chamada do processo `GATK_HAPLOTYPECALLER`.

!!! note

    Observe que as entradas são fornecidas exatamente na mesma ordem na chamada do processo como estão listadas no bloco de entrada do processo.
    No Nextflow, as entradas são posicionais, o que significa que você _deve_ seguir a mesma ordem; e, é claro, deve haver o mesmo número de elementos.

### 2.3. Configurar o tratamento da saída

Precisamos adicionar as novas saídas à declaração publish e configurar para onde elas vão.

#### 2.3.1. Adicionar alvos de publicação para as saídas de chamada de variantes

Adicione as saídas VCF e índice à seção `publish:`:

=== "Depois"

    ```groovy title="genomics.nf" linenums="45" hl_lines="3-4"
        publish:
        bam_index = SAMTOOLS_INDEX.out
        vcf = GATK_HAPLOTYPECALLER.out.vcf
        vcf_idx = GATK_HAPLOTYPECALLER.out.idx
    }
    ```

=== "Antes"

    ```groovy title="genomics.nf" linenums="45"
        publish:
        bam_index = SAMTOOLS_INDEX.out
    }
    ```

Agora precisamos dizer ao Nextflow onde colocar as novas saídas.

#### 2.3.2. Configurar os novos alvos de saída

Adicione entradas para os alvos `vcf` e `vcf_idx` no bloco `output {}`, publicando ambos em um subdiretório `vcf/`:

=== "Depois"

    ```groovy title="genomics.nf" linenums="51" hl_lines="5-10"
    output {
        bam_index {
            path 'bam'
        }
        vcf {
            path 'vcf'
        }
        vcf_idx {
            path 'vcf'
        }
    }
    ```

=== "Antes"

    ```groovy title="genomics.nf" linenums="49"
    output {
        bam_index {
            path 'bam'
        }
    }
    ```

O VCF e seu índice são publicados como alvos separados que ambos vão para o subdiretório `vcf/`.

### 2.4. Executar o fluxo de trabalho

Execute o fluxo de trabalho expandido, adicionando `-resume` desta vez para que não tenhamos que executar a etapa de indexação novamente.

```bash
nextflow run genomics.nf -profile test -resume
```

??? success "Saída do comando"

    ```console
    N E X T F L O W   ~  version 25.10.2

    ┃ Launching `genomics.nf` [grave_volta] DSL2 - revision: 4790abc96a

    executor >  local (1)
    [2a/e69536] SAMTOOLS_INDEX (1)       | 1 of 1, cached: 1 ✔
    [53/e18e98] GATK_HAPLOTYPECALLER (1) | 1 of 1 ✔
    ```

Agora, se olharmos para a saída do console, vemos os dois processos listados.

O primeiro processo foi pulado graças ao cache, como esperado, enquanto o segundo processo foi executado já que é totalmente novo.

Você encontrará os arquivos de saída no diretório de resultados (como links simbólicos para o diretório de trabalho).

??? abstract "Conteúdo do diretório"

    ```console
    results/
    ├── bam/
    │   └── reads_mother.bam.bai -> ...
    └── vcf/
        ├── reads_mother.bam.vcf -> ...
        └── reads_mother.bam.vcf.idx -> ...
    ```

Se você abrir o arquivo VCF, deverá ver o mesmo conteúdo do arquivo que você gerou ao executar o comando GATK diretamente no contêiner.

??? abstract "Conteúdo do arquivo"

    ```console title="reads_mother.bam.vcf" linenums="26"
    #CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	reads_mother
    20_10037292_10066351	3480	.	C	CT	503.03	.	AC=2;AF=1.00;AN=2;DP=23;ExcessHet=0.0000;FS=0.000;MLEAC=2;MLEAF=1.00;MQ=60.00;QD=27.95;SOR=1.179	GT:AD:DP:GQ:PL	1/1:0,18:18:54:517,54,0
    20_10037292_10066351	3520	.	AT	A	609.03	.	AC=2;AF=1.00;AN=2;DP=18;ExcessHet=0.0000;FS=0.000;MLEAC=2;MLEAF=1.00;MQ=60.00;QD=33.83;SOR=0.693	GT:AD:DP:GQ:PL	1/1:0,18:18:54:623,54,0
    20_10037292_10066351	3529	.	T	A	155.64	.	AC=1;AF=0.500;AN=2;BaseQRankSum=-0.544;DP=21;ExcessHet=0.0000;FS=1.871;MLEAC=1;MLEAF=0.500;MQ=60.00;MQRankSum=0.000;QD=7.78;ReadPosRankSum=-1.158;SOR=1.034	GT:AD:DP:GQ:PL	0/1:12,8:20:99:163,0,328
    ```

Esta é a saída que nos importa gerar para cada amostra em nosso estudo.

### Conclusão

Você sabe como fazer um fluxo de trabalho modular de duas etapas que faz trabalho de análise real e é capaz de lidar com as idiossincrasia dos formatos de arquivo de genômica, como os arquivos acessórios.

### O que vem a seguir?

Fazer o fluxo de trabalho lidar com várias amostras em massa.

---

## 3. Adaptar o fluxo de trabalho para executar em um lote de amostras

É muito bom ter um fluxo de trabalho que possa automatizar o processamento de uma única amostra, mas e se você tiver 1000 amostras?
Você precisa escrever um script bash que percorre todas as suas amostras?

Não, felizmente! Basta fazer um pequeno ajuste no código e o Nextflow cuidará disso para você também.

### 3.1. Atualizar a entrada para listar três amostras

Para executar em várias amostras, atualize o perfil de teste para fornecer um array de caminhos de arquivo em vez de um único.
Esta é uma maneira rápida de testar a execução de várias amostras; na próxima etapa, mudaremos para uma abordagem mais escalável usando um arquivo de entradas.

Primeiro, comente a anotação de tipo na declaração do parâmetro, já que arrays não podem usar declarações tipadas:

=== "Depois"

    ```groovy title="genomics.nf" linenums="10" hl_lines="1-2"
        // Primary input (array of three samples)
        reads_bam //: Path
    ```

=== "Antes"

    ```groovy title="genomics.nf" linenums="10"
        // Primary input
        reads_bam: Path
    ```

Em seguida, atualize o perfil de teste para listar todas as três amostras:

=== "Depois"

    ```groovy title="nextflow.config" linenums="4" hl_lines="2-6"
    test {
        params.reads_bam = [
            "${projectDir}/data/bam/reads_mother.bam",
            "${projectDir}/data/bam/reads_father.bam",
            "${projectDir}/data/bam/reads_son.bam"
        ]
        params.reference = "${projectDir}/data/ref/ref.fasta"
        params.reference_index = "${projectDir}/data/ref/ref.fasta.fai"
        params.reference_dict = "${projectDir}/data/ref/ref.dict"
        params.intervals = "${projectDir}/data/ref/intervals.bed"
    }
    ```

=== "Antes"

    ```groovy title="nextflow.config" linenums="4"
    test {
        params.reads_bam = "${projectDir}/data/bam/reads_mother.bam"
        params.reference = "${projectDir}/data/ref/ref.fasta"
        params.reference_index = "${projectDir}/data/ref/ref.fasta.fai"
        params.reference_dict = "${projectDir}/data/ref/ref.dict"
        params.intervals = "${projectDir}/data/ref/intervals.bed"
    }
    ```

A factory de canal no corpo do fluxo de trabalho (`.fromPath`) aceita vários caminhos de arquivo tão bem quanto um único, então nenhuma outra alteração é necessária.

### 3.2. Executar o fluxo de trabalho

Tente executar o fluxo de trabalho agora que a estrutura está configurada para executar em todas as três amostras de teste.

```bash
nextflow run genomics.nf -profile test -resume
```

Coisa engraçada: isso _pode funcionar_, OU _pode falhar_. Por exemplo, aqui está uma execução que teve sucesso:

??? success "Saída do comando"

    ```console
    N E X T F L O W   ~  version 25.10.2

    ┃ Launching `genomics.nf` [peaceful_yalow] DSL2 - revision: a256d113ad

    executor >  local (6)
    [4f/7071b0] SAMTOOLS_INDEX (3)       | 3 of 3, cached: 1 ✔
    [7a/89bc43] GATK_HAPLOTYPECALLER (2) | 3 of 3, cached: 1 ✔
    ```

Se a execução do seu fluxo de trabalho teve sucesso, execute-o novamente até obter um erro como este:

??? failure "Saída do comando"

    ```console
    N E X T F L O W   ~  version 25.10.2

    ┃ Launching `genomics.nf` [loving_pasteur] DSL2 - revision: d2a8e63076

    executor >  local (4)
    [01/eea165] SAMTOOLS_INDEX (2)       | 3 of 3, cached: 1 ✔
    [a5/fa9fd0] GATK_HAPLOTYPECALLER (3) | 1 of 3, cached: 1
    ERROR ~ Error executing process > 'GATK_HAPLOTYPECALLER (2)'

    Caused by:
      Process `GATK_HAPLOTYPECALLER (2)` terminated with an error exit status (2)

    Command executed:

      gatk HaplotypeCaller         -R ref.fasta         -I reads_father.bam         -O reads_father.bam.vcf         -L intervals.bed

    Command exit status:
      2

    Command error:
      ...
      A USER ERROR has occurred: Traversal by intervals was requested but some input files are not indexed.
      ...
    ```

Se você olhar para a saída de erro do comando GATK, haverá uma linha como esta:

```console
A USER ERROR has occurred: Traversal by intervals was requested but some input files are not indexed.
```

Bem, isso é estranho, considerando que indexamos explicitamente os arquivos BAM na primeira etapa do fluxo de trabalho. Poderia haver algo errado com a estrutura?

### 3.3. Solucionar o problema

Vamos inspecionar os diretórios de trabalho e usar o operador `view()` para descobrir o que deu errado.

#### 3.3.1. Verificar os diretórios de trabalho para as chamadas relevantes

Dê uma olhada dentro do diretório de trabalho para a chamada do processo `GATK_HAPLOTYPECALLER` falhada listada na saída do console.

??? abstract "Conteúdo do diretório"

    ```console
    work/a5/fa9fd0994b6beede5fb9ea073596c2
    ├── intervals.bed -> /workspaces/training/nf4-science/genomics/data/ref/intervals.bed
    ├── reads_father.bam.bai -> /workspaces/training/nf4-science/genomics/work/01/eea16597bd6e810fb4cf89e60f8c2d/reads_father.bam.bai
    ├── reads_son.bam -> /workspaces/training/nf4-science/genomics/data/bam/reads_son.bam
    ├── reads_son.bam.vcf
    ├── reads_son.bam.vcf.idx
    ├── ref.dict -> /workspaces/training/nf4-science/genomics/data/ref/ref.dict
    ├── ref.fasta -> /workspaces/training/nf4-science/genomics/data/ref/ref.fasta
    └── ref.fasta.fai -> /workspaces/training/nf4-science/genomics/data/ref/ref.fasta.fai
    ```

Preste atenção especial aos nomes do arquivo BAM e do índice BAM que estão listados neste diretório: `reads_son.bam` e `reads_father.bam.bai`.

O que diabos? O Nextflow preparou um arquivo de índice no diretório de trabalho desta chamada de processo, mas é o errado. Como isso pode ter acontecido?

#### 3.3.2. Usar o [operador view()](https://www.nextflow.io/docs/latest/reference/operator.html#view) para inspecionar o conteúdo do canal

Adicione estas duas linhas no corpo do fluxo de trabalho antes da chamada do processo `GATK_HAPLOTYPECALLER` para visualizar o conteúdo do canal:

=== "Depois"

    ```groovy title="genomics.nf" hl_lines="3-5"
        SAMTOOLS_INDEX(reads_ch)

        // temporary diagnostics
        reads_ch.view()
        SAMTOOLS_INDEX.out.view()

        // Call variants from the indexed BAM file
        GATK_HAPLOTYPECALLER(
    ```

=== "Antes"

    ```groovy title="genomics.nf"
        SAMTOOLS_INDEX(reads_ch)

        // Call variants from the indexed BAM file
        GATK_HAPLOTYPECALLER(
    ```

Em seguida, execute o comando do fluxo de trabalho novamente.

```bash
nextflow run genomics.nf -profile test
```

Mais uma vez, isso pode ter sucesso ou falhar. Aqui está como a saída das duas chamadas `.view()` se parece para uma execução falhada:

```console
/workspaces/training/nf4-science/genomics/data/bam/reads_mother.bam
/workspaces/training/nf4-science/genomics/data/bam/reads_father.bam
/workspaces/training/nf4-science/genomics/data/bam/reads_son.bam
/workspaces/training/nf4-science/genomics/work/9c/53492e3518447b75363e1cd951be4b/reads_father.bam.bai
/workspaces/training/nf4-science/genomics/work/cc/37894fffdf6cc84c3b0b47f9b536b7/reads_son.bam.bai
/workspaces/training/nf4-science/genomics/work/4d/dff681a3d137ba7d9866e3d9307bd0/reads_mother.bam.bai
```

As três primeiras linhas correspondem ao canal de entrada e a segunda, ao canal de saída.
Você pode ver que os arquivos BAM e os arquivos de índice para as três amostras não estão listados na mesma ordem!

!!! note

    Quando você chama um processo Nextflow em um canal contendo vários elementos, o Nextflow tentará paralelizar a execução o máximo possível e coletará saídas em qualquer ordem em que fiquem disponíveis.
    A consequência é que as saídas correspondentes podem ser coletadas em uma ordem diferente da ordem em que as entradas originais foram fornecidas.

Como está escrito atualmente, nosso script de fluxo de trabalho assume que os arquivos de índice sairão da etapa de indexação listados na mesma ordem mãe/pai/filho que as entradas foram fornecidas.
Mas isso não é garantido, e é por isso que às vezes (embora nem sempre) os arquivos errados são emparelhados na segunda etapa.

Para corrigir isso, precisamos garantir que os arquivos BAM e seus arquivos de índice viajem juntos pelos canais.

!!! tip

    As instruções `view()` no código do fluxo de trabalho não fazem nada, então não é um problema deixá-las.
    No entanto, elas irão desordenar sua saída do console, então recomendamos removê-las quando você terminar de solucionar o problema.

### 3.4. Atualizar o fluxo de trabalho para lidar com os arquivos de índice corretamente

A correção é empacotar cada arquivo BAM com seu índice em uma tupla, depois atualizar o processo downstream e a estrutura do fluxo de trabalho para corresponder.

#### 3.4.1. Alterar a saída do módulo SAMTOOLS_INDEX para uma tupla

A maneira mais simples de garantir que um arquivo BAM e seu índice permaneçam intimamente associados é empacotá-los juntos em uma tupla saindo da tarefa de índice.

!!! note

    Uma **tupla** é uma lista finita e ordenada de elementos que é comumente usada para retornar vários valores de uma função. As tuplas são particularmente úteis para passar várias entradas ou saídas entre processos, preservando sua associação e ordem.

Atualize a saída em `modules/samtools_index.nf` para incluir o arquivo BAM:

=== "Depois"

    ```groovy title="modules/samtools_index.nf" linenums="14" hl_lines="2"
        output:
        tuple path(input_bam), path("${input_bam}.bai")
    ```

=== "Antes"

    ```groovy title="modules/samtools_index.nf" linenums="14"
        output:
        path "${input_bam}.bai"
    ```

Dessa forma, cada arquivo de índice será fortemente acoplado ao seu arquivo BAM original, e a saída geral da etapa de indexação será um único canal contendo pares de arquivos.

#### 3.4.2. Alterar a entrada do módulo GATK_HAPLOTYPECALLER para aceitar uma tupla

Como mudamos a 'forma' da saída do primeiro processo, precisamos atualizar a definição de entrada do segundo processo para corresponder.

Atualize `modules/gatk_haplotypecaller.nf`:

=== "Depois"

    ```groovy title="modules/gatk_haplotypecaller.nf" linenums="11" hl_lines="2"
        input:
        tuple path(input_bam), path(input_bam_index)
    ```

=== "Antes"

    ```groovy title="modules/gatk_haplotypecaller.nf" linenums="11"
        input:
        path input_bam
        path input_bam_index
    ```

Agora precisamos atualizar o fluxo de trabalho para refletir a nova estrutura de tupla na chamada do processo e nos alvos de publicação.

#### 3.4.3. Atualizar a chamada para GATK_HAPLOTYPECALLER no fluxo de trabalho

Não precisamos mais fornecer o `reads_ch` original ao processo `GATK_HAPLOTYPECALLER`, já que o arquivo BAM agora está empacotado na saída do canal por `SAMTOOLS_INDEX`.

Atualize a chamada em `genomics.nf`:

=== "Depois"

    ```groovy title="genomics.nf" linenums="42" hl_lines="2"
        GATK_HAPLOTYPECALLER(
            SAMTOOLS_INDEX.out,
    ```

=== "Antes"

    ```groovy title="genomics.nf" linenums="42"
        GATK_HAPLOTYPECALLER(
            reads_ch,
            SAMTOOLS_INDEX.out,
    ```

Finalmente, precisamos atualizar os alvos de publicação para refletir a nova estrutura de saída.

#### 3.4.4. Atualizar o alvo de publicação para a saída do BAM indexado

Como a saída do SAMTOOLS_INDEX agora é uma tupla contendo tanto o arquivo BAM quanto seu índice, renomeie o alvo de publicação de `bam_index` para `indexed_bam` para melhor refletir seu conteúdo:

=== "Depois"

    ```groovy title="genomics.nf" linenums="46" hl_lines="2 8"
        publish:
        indexed_bam = SAMTOOLS_INDEX.out
        vcf = GATK_HAPLOTYPECALLER.out.vcf
        vcf_idx = GATK_HAPLOTYPECALLER.out.idx
    }

    output {
        indexed_bam {
            path 'bam'
        }
        vcf {
            path 'vcf'
        }
        vcf_idx {
            path 'vcf'
        }
    }
    ```

=== "Antes"

    ```groovy title="genomics.nf" linenums="46"
        publish:
        bam_index = SAMTOOLS_INDEX.out
        vcf = GATK_HAPLOTYPECALLER.out.vcf
        vcf_idx = GATK_HAPLOTYPECALLER.out.idx
    }

    output {
        bam_index {
            path 'bam'
        }
        vcf {
            path 'vcf'
        }
        vcf_idx {
            path 'vcf'
        }
    }
    ```

Com essas mudanças, o BAM e seu índice são garantidos de viajar juntos, então o emparelhamento sempre estará correto.

### 3.5. Executar o fluxo de trabalho corrigido

Execute o fluxo de trabalho novamente para garantir que isso funcionará de forma confiável daqui para frente.

```bash
nextflow run genomics.nf -profile test
```

Desta vez (e todas as vezes) tudo deve funcionar corretamente:

??? success "Saída do comando"

    ```console
    N E X T F L O W   ~  version 25.10.2

    ┃ Launching `genomics.nf` [special_goldstine] DSL2 - revision: 4cbbf6ea3e

    executor >  local (6)
    [d6/10c2c4] SAMTOOLS_INDEX (1)       | 3 of 3 ✔
    [88/1783aa] GATK_HAPLOTYPECALLER (2) | 3 of 3 ✔
    ```

O diretório de resultados agora contém tanto os arquivos BAM quanto BAI para cada amostra (da tupla), junto com as saídas VCF:

??? abstract "Conteúdo do diretório de resultados"

    ```console
    results/
    ├── bam/
    │   ├── reads_father.bam -> ...
    │   ├── reads_father.bam.bai -> ...
    │   ├── reads_mother.bam -> ...
    │   ├── reads_mother.bam.bai -> ...
    │   ├── reads_son.bam -> ...
    │   └── reads_son.bam.bai -> ...
    └── vcf/
        ├── reads_father.bam.vcf -> ...
        ├── reads_father.bam.vcf.idx -> ...
        ├── reads_mother.bam.vcf -> ...
        ├── reads_mother.bam.vcf.idx -> ...
        ├── reads_son.bam.vcf -> ...
        └── reads_son.bam.vcf.idx -> ...
    ```

Ao empacotar arquivos associados em tuplas, garantimos que os arquivos corretos sempre viajem juntos pelo fluxo de trabalho.
O fluxo de trabalho agora processa qualquer número de amostras de forma confiável, mas listá-las individualmente no config não é muito escalável.
Na próxima etapa, mudaremos para a leitura de entradas de um arquivo.

### Conclusão

Você sabe como fazer seu fluxo de trabalho executar em várias amostras (independentemente).

### O que vem a seguir?

Facilitar o tratamento de amostras em massa.

---

## 4. Fazer o fluxo de trabalho aceitar um arquivo de texto contendo um lote de arquivos de entrada

Uma maneira muito comum de fornecer vários arquivos de dados de entrada para um fluxo de trabalho é fazê-lo com um arquivo de texto contendo os caminhos dos arquivos.
Pode ser tão simples quanto um arquivo de texto listando um caminho de arquivo por linha e nada mais, ou o arquivo pode conter metadados adicionais, caso em que é frequentemente chamado de samplesheet.

Aqui vamos mostrar como fazer o caso simples.

### 4.1. Examinar o arquivo de texto fornecido listando os caminhos dos arquivos de entrada

Já fizemos um arquivo de texto listando os caminhos dos arquivos de entrada, chamado `sample_bams.txt`, que você pode encontrar no diretório `data/`.

```txt title="sample_bams.txt"
/workspaces/training/nf4-science/genomics/data/bam/reads_mother.bam
/workspaces/training/nf4-science/genomics/data/bam/reads_father.bam
/workspaces/training/nf4-science/genomics/data/bam/reads_son.bam
```

Como você pode ver, listamos um caminho de arquivo por linha, e eles são caminhos absolutos.

!!! note

    Os arquivos que estamos usando aqui estão apenas no sistema de arquivos local do seu GitHub Codespaces, mas também poderíamos apontar para arquivos no armazenamento em nuvem.
    Se você não estiver usando o ambiente Codespaces fornecido, pode ser necessário adaptar os caminhos dos arquivos para corresponder à sua configuração local.

### 4.2. Atualizar o parâmetro e o perfil de teste

Mude o parâmetro `reads_bam` para apontar para o arquivo `sample_bams.txt` em vez de listar amostras individuais.

Restaure a anotação de tipo no bloco params (já que é um único caminho novamente):

=== "Depois"

    ```groovy title="genomics.nf" linenums="10" hl_lines="1-2"
        // Primary input (file of input files, one per line)
        reads_bam: Path
    ```

=== "Antes"

    ```groovy title="genomics.nf" linenums="10"
        // Primary input (array of three samples)
        reads_bam
    ```

Em seguida, atualize o perfil de teste para apontar para o arquivo de texto:

=== "Depois"

    ```groovy title="nextflow.config" linenums="4" hl_lines="2"
    test {
        params.reads_bam = "${projectDir}/data/sample_bams.txt"
        params.reference = "${projectDir}/data/ref/ref.fasta"
        params.reference_index = "${projectDir}/data/ref/ref.fasta.fai"
        params.reference_dict = "${projectDir}/data/ref/ref.dict"
        params.intervals = "${projectDir}/data/ref/intervals.bed"
    }
    ```

=== "Antes"

    ```groovy title="nextflow.config" linenums="4"
    test {
        params.reads_bam = [
            "${projectDir}/data/bam/reads_mother.bam",
            "${projectDir}/data/bam/reads_father.bam",
            "${projectDir}/data/bam/reads_son.bam"
        ]
        params.reference = "${projectDir}/data/ref/ref.fasta"
        params.reference_index = "${projectDir}/data/ref/ref.fasta.fai"
        params.reference_dict = "${projectDir}/data/ref/ref.dict"
        params.intervals = "${projectDir}/data/ref/intervals.bed"
    }
    ```

A lista de arquivos não vive mais no código, o que é um grande passo na direção certa.

### 4.3. Atualizar a factory de canal para ler linhas de um arquivo

Atualmente, nossa factory de canal de entrada trata quaisquer arquivos que damos a ela como as entradas de dados que queremos alimentar no processo de indexação.
Como agora estamos dando a ela um arquivo que lista caminhos de arquivos de entrada, precisamos mudar seu comportamento para analisar o arquivo e tratar os caminhos de arquivo que ele contém como as entradas de dados.

Podemos fazer isso usando o mesmo padrão que usamos na [Parte 2 do Hello Nextflow](../../hello_nextflow/02_hello_channels.md#42-use-the-splitcsv-operator-to-parse-the-file): aplicando o operador [`splitCsv()`](https://nextflow.io/docs/latest/reference/operator.html#splitcsv) para analisar o arquivo, depois uma operação `map` para selecionar o primeiro campo de cada linha.

=== "Depois"

    ```groovy title="genomics.nf" linenums="24" hl_lines="1-4"
        // Create input channel from a CSV file listing input file paths
        reads_ch = Channel.fromPath(params.reads_bam)
                .splitCsv()
                .map { line -> file(line[0]) }
    ```

=== "Antes"

    ```groovy title="genomics.nf" linenums="24"
        // Create input channel (single file via CLI parameter)
        reads_ch = channel.fromPath(params.reads_bam)
    ```

Tecnicamente poderíamos fazer isso mais simplesmente usando o operador [`.splitText()`](https://www.nextflow.io/docs/latest/reference/operator.html#operator-splittext), já que nosso arquivo de entrada atualmente contém apenas caminhos de arquivo.
No entanto, ao usar o operador mais versátil `splitCsv` (suplementado por `map`), podemos tornar nosso fluxo de trabalho à prova de futuro caso decidamos adicionar metadados ao arquivo contendo caminhos de arquivo.

!!! tip

    Se você não tem confiança de que entende o que os operadores estão fazendo aqui, esta é outra ótima oportunidade para usar o operador `.view()` para ver como o conteúdo do canal se parece antes e depois de aplicá-los.

### 4.4. Executar o fluxo de trabalho

Execute o fluxo de trabalho mais uma vez. Isso deve produzir o mesmo resultado de antes, certo?

```bash
nextflow run genomics.nf -profile test -resume
```

??? success "Saída do comando"

    ```console
    N E X T F L O W   ~  version 25.10.2

    ┃ Launching `genomics.nf` [sick_albattani] DSL2 - revision: 46d84642f6

    [18/23b4bb] SAMTOOLS_INDEX (1)       | 3 of 3, cached: 3 ✔
    [12/f727bb] GATK_HAPLOTYPECALLER (3) | 3 of 3, cached: 3 ✔
    ```

Sim! Na verdade, o Nextflow detecta corretamente que as chamadas do processo são exatamente as mesmas e nem se preocupa em executar tudo novamente, já que estávamos executando com `-resume`.

E é isso! Nosso fluxo de trabalho simples de chamada de variantes tem todos os recursos básicos que queríamos.

### Conclusão

Você sabe como fazer um fluxo de trabalho modular de várias etapas para indexar um arquivo BAM e aplicar chamada de variantes por amostra usando GATK.

Mais geralmente, você aprendeu como usar componentes e lógica essenciais do Nextflow para construir um pipeline de genômica simples que faz trabalho real, levando em conta as idiossincrasia dos formatos de arquivo de genômica e requisitos de ferramentas.

### O que vem a seguir?

Celebre seu sucesso e faça uma pausa extra longa!

Na próxima parte deste curso, você aprenderá como transformar este fluxo de trabalho simples de chamada de variantes por amostra para aplicar chamada de variantes conjunta aos dados.
