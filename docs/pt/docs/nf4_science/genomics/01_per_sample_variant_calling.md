# Parte 1: Chamada de variantes por amostra

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Tradução assistida por IA - [saiba mais e sugira melhorias](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

Na primeira parte deste curso, mostramos como construir um pipeline simples de chamada de variantes que aplica a chamada de variantes do GATK a amostras individuais de sequenciamento.

### Visão geral do método

A chamada de variantes é um método de análise genômica que visa identificar variações em uma sequência genômica em relação a um genoma de referência.
Aqui vamos usar ferramentas e métodos projetados para chamar variantes curtas, _ou seja_, SNPs e indels.

![Pipeline GATK](img/gatk-pipeline.png)

Um pipeline completo de chamada de variantes normalmente envolve muitas etapas, incluindo mapeamento para a referência (às vezes chamado de alinhamento genômico) e filtragem e priorização de variantes.
Para simplificar, nesta parte do curso vamos focar apenas na parte de chamada de variantes.

### Conjunto de dados

Fornecemos os seguintes dados e recursos relacionados:

- **Um genoma de referência** consistindo de uma pequena região do cromossomo 20 humano (de hg19/b37) e seus arquivos acessórios (índice e dicionário de sequências).
- **Três amostras de sequenciamento de genoma completo** correspondentes a um trio familiar (mãe, pai e filho), que foram reduzidas a uma pequena fatia de dados no cromossomo 20 para manter os tamanhos de arquivo pequenos.
  Estes são dados de sequenciamento Illumina de leitura curta que já foram mapeados para o genoma de referência, fornecidos em formato [BAM](https://samtools.github.io/hts-specs/SAMv1.pdf) (Binary Alignment Map, uma versão comprimida de SAM, Sequence Alignment Map).
- **Uma lista de intervalos genômicos**, ou seja, coordenadas no genoma onde nossas amostras têm dados adequados para chamar variantes, fornecida em formato BED.

### Fluxo de trabalho

Nesta parte do curso, vamos desenvolver um fluxo de trabalho que faz o seguinte:

1. Gerar um arquivo de índice para cada arquivo BAM de entrada usando [Samtools](https://www.htslib.org/)
2. Executar o GATK HaplotypeCaller em cada arquivo BAM de entrada para gerar chamadas de variantes por amostra em VCF (Variant Call Format)

<figure class="excalidraw">
--8<-- "docs/en/docs/nf4_science/genomics/img/hello-gatk-1.svg"
</figure>

!!! note "Nota"

    Arquivos de índice são uma característica comum dos formatos de arquivo de bioinformática; eles contêm informações sobre a estrutura do arquivo principal que permite que ferramentas como GATK acessem um subconjunto dos dados sem ter que ler o arquivo inteiro.
    Isso é importante por causa de quão grandes esses arquivos podem ficar.

---

## 0. Aquecimento: Testar os comandos Samtools e GATK interativamente

Primeiro queremos experimentar os comandos manualmente antes de tentar envolvê-los em um fluxo de trabalho.
As ferramentas que precisamos (Samtools e GATK) não estão instaladas no ambiente GitHub Codespaces, então vamos usá-las via contêineres (veja [Hello Containers](../../hello_nextflow/05_hello_containers.md)).

!!! note "Nota"

     Certifique-se de estar no diretório `nf4-science/genomics` para que a última parte do caminho mostrado quando você digita `pwd` seja `genomics`.

### 0.1. Indexar um arquivo BAM de entrada com Samtools

Vamos baixar um contêiner Samtools, iniciá-lo interativamente e executar o comando `samtools index` em um dos arquivos BAM.

#### 0.1.1. Baixar o contêiner Samtools

```bash
docker pull community.wave.seqera.io/library/samtools:1.20--b5dfbd93de237464
```

<!--
??? success "Saída do comando"

    ```console

    ```
-->

#### 0.1.2. Iniciar o contêiner Samtools interativamente

```bash
docker run -it -v ./data:/data community.wave.seqera.io/library/samtools:1.20--b5dfbd93de237464
```

<!--
??? success "Saída do comando"

    ```console

    ```
-->

#### 0.1.3. Executar o comando de indexação

A [documentação do Samtools](https://www.htslib.org/doc/samtools-index.html) nos fornece a linha de comando para executar para indexar um arquivo BAM.

Só precisamos fornecer o arquivo de entrada; a ferramenta irá gerar automaticamente um nome para a saída adicionando `.bai` ao nome do arquivo de entrada.

```bash
samtools index /data/bam/reads_mother.bam
```

Isso deve ser concluído imediatamente, e você deve ver agora um arquivo chamado `reads_mother.bam.bai` no mesmo diretório que o arquivo BAM de entrada original.

??? abstract "Conteúdo do diretório"

    ```console
    data/bam/
    ├── reads_father.bam
    ├── reads_mother.bam
    ├── reads_mother.bam.bai
    └── reads_son.bam
    ```

#### 0.1.4. Sair do contêiner Samtools

```bash
exit
```

### 0.2. Chamar variantes com GATK HaplotypeCaller

Vamos baixar um contêiner GATK, iniciá-lo interativamente e executar o comando `gatk HaplotypeCaller` no arquivo BAM que acabamos de indexar.

#### 0.2.1. Baixar o contêiner GATK

```bash
docker pull community.wave.seqera.io/library/gatk4:4.5.0.0--730ee8817e436867
```

<!--
??? success "Saída do comando"

    ```console

    ```
-->

#### 0.2.2. Iniciar o contêiner GATK interativamente

```bash
docker run -it -v ./data:/data community.wave.seqera.io/library/gatk4:4.5.0.0--730ee8817e436867
```

<!--
??? success "Saída do comando"

    ```console

    ```
-->

#### 0.2.3. Executar o comando de chamada de variantes

A [documentação do GATK](https://gatk.broadinstitute.org/hc/en-us/articles/21905025322523-HaplotypeCaller) nos fornece a linha de comando para executar para realizar chamada de variantes em um arquivo BAM.

Precisamos fornecer o arquivo BAM de entrada (`-I`) assim como o genoma de referência (`-R`), um nome para o arquivo de saída (`-O`) e uma lista de intervalos genômicos para analisar (`-L`).

No entanto, não precisamos especificar o caminho para o arquivo de índice; a ferramenta irá procurá-lo automaticamente no mesmo diretório, baseado na convenção estabelecida de nomenclatura e colocalização.
O mesmo se aplica aos arquivos acessórios do genoma de referência (arquivos de índice e dicionário de sequências, `*.fai` e `*.dict`).

```bash
gatk HaplotypeCaller \
        -R /data/ref/ref.fasta \
        -I /data/bam/reads_mother.bam \
        -O reads_mother.vcf \
        -L /data/ref/intervals.bed
```

<!--
??? success "Saída do comando"

    ```console

    ```
-->

O arquivo de saída `reads_mother.vcf` é criado dentro do seu diretório de trabalho no contêiner, então você não o verá no explorador de arquivos do VS Code a menos que você mude o caminho do arquivo de saída.
No entanto, é um arquivo de teste pequeno, então você pode usar `cat` para abri-lo e visualizar o conteúdo.
Se você rolar até o início do arquivo, encontrará um cabeçalho composto de muitas linhas de metadados, seguido por uma lista de chamadas de variantes, uma por linha.

```console title="reads_mother.vcf" linenums="26"
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	reads_mother
20_10037292_10066351	3480	.	C	CT	503.03	.	AC=2;AF=1.00;AN=2;DP=23;ExcessHet=0.0000;FS=0.000;MLEAC=2;MLEAF=1.00;MQ=60.00;QD=27.95;SOR=1.179	GT:AD:DP:GQ:PL	1/1:0,18:18:54:517,54,0
20_10037292_10066351	3520	.	AT	A	609.03	.	AC=2;AF=1.00;AN=2;DP=18;ExcessHet=0.0000;FS=0.000;MLEAC=2;MLEAF=1.00;MQ=60.00;QD=33.83;SOR=0.693	GT:AD:DP:GQ:PL	1/1:0,18:18:54:623,54,0
20_10037292_10066351	3529	.	T	A	155.64	.	AC=1;AF=0.500;AN=2;BaseQRankSum=-0.544;DP=21;ExcessHet=0.0000;FS=1.871;MLEAC=1;MLEAF=0.500;MQ=60.00;MQRankSum=0.000;QD=7.78;ReadPosRankSum=-1.158;SOR=1.034	GT:AD:DP:GQ:PL	0/1:12,8:20:99:163,0,328
```

Cada linha descreve uma possível variante identificada nos dados de sequenciamento da amostra. Para orientação sobre como interpretar o formato VCF, veja [este artigo útil](https://www.ebi.ac.uk/training/online/courses/human-genetic-variation-introduction/variant-identification-and-analysis/understanding-vcf-format/).

O arquivo VCF de saída é acompanhado por um arquivo de índice chamado `reads_mother.vcf.idx` que foi criado automaticamente pelo GATK.
Ele tem a mesma função que o arquivo de índice BAM, permitir que ferramentas busquem e recuperem subconjuntos de dados sem carregar o arquivo inteiro.

#### 0.2.4. Sair do contêiner GATK

```bash
exit
```

### Conclusão

Você sabe como testar os comandos de indexação do Samtools e chamada de variantes do GATK em seus respectivos contêineres.

### O que vem a seguir?

Aprender como envolver esses mesmos comandos em um fluxo de trabalho de duas etapas que usa contêineres para executar o trabalho.

---

## 1. Escrever um fluxo de trabalho de etapa única que executa Samtools index em um arquivo BAM

Fornecemos um arquivo de fluxo de trabalho, `genomics-1.nf`, que delineia as partes principais do fluxo de trabalho.
Ele não é funcional; seu propósito é apenas servir como um esqueleto que você usará para escrever o fluxo de trabalho real.

### 1.1. Definir o processo de indexação

Vamos começar escrevendo um processo, que chamaremos de `SAMTOOLS_INDEX`, descrevendo a operação de indexação.

```groovy title="genomics-1.nf" linenums="9"
/*
 * Gerar arquivo de índice BAM
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

Você deve reconhecer todas as partes do que aprendeu na Parte 1 & Parte 2 desta série de treinamento.

Este processo vai exigir que passemos um caminho de arquivo via entrada `input_bam`, então vamos configurar isso em seguida.

### 1.2. Adicionar uma declaração de parâmetro de entrada

No topo do arquivo, na seção `Pipeline parameters`, declaramos um parâmetro CLI chamado `reads_bam` e damos a ele um valor padrão.
Dessa forma, podemos ser preguiçosos e não especificar a entrada quando digitarmos o comando para iniciar o pipeline (para fins de desenvolvimento).

```groovy title="genomics-1.nf" linenums="3"
/*
 * Parâmetros do pipeline
 */
params {
    // Entrada principal
    reads_bam: Path = "${projectDir}/data/bam/reads_mother.bam"
}
```

Agora temos um processo pronto, assim como um parâmetro para dar a ele uma entrada para executar, então vamos conectar essas coisas juntas.

!!! note "Nota"

    `${projectDir}` é uma variável embutida do Nextflow que aponta para o diretório onde o script de fluxo de trabalho Nextflow atual (`genomics-1.nf`) está localizado.

    Isso facilita a referência a arquivos, diretórios de dados e outros recursos incluídos no repositório do fluxo de trabalho sem codificar caminhos absolutos.

### 1.3. Adicionar bloco de fluxo de trabalho para executar SAMTOOLS_INDEX

No bloco `workflow`, precisamos configurar um **canal** para alimentar a entrada do processo `SAMTOOLS_INDEX`; então podemos chamar o próprio processo para executar no conteúdo desse canal.

```groovy title="genomics-1.nf" linenums="24"
workflow {

    main:
    // Criar canal de entrada (arquivo único via parâmetro CLI)
    reads_ch = channel.fromPath(params.reads_bam)

    // Criar arquivo de índice para arquivo BAM de entrada
    SAMTOOLS_INDEX(reads_ch)

    publish:
    bam_index = SAMTOOLS_INDEX.out
}
```

O bloco workflow tem duas seções:

- `main:` contém as operações de canal e chamadas de processo
- `publish:` declara quais saídas devem ser publicadas, atribuindo-as a alvos nomeados

Você notará que estamos usando a mesma factory de canal `.fromPath` que usamos em [Hello Channels](../../hello_nextflow/02_hello_channels.md).
De fato, estamos fazendo algo muito similar.
A diferença é que estamos dizendo ao Nextflow para apenas carregar o próprio caminho do arquivo no canal como um elemento de entrada, em vez de ler seu conteúdo.

### 1.4. Adicionar um bloco de saída para definir onde os resultados são publicados

Após o bloco workflow, adicionamos um bloco `output` que especifica onde publicar as saídas do fluxo de trabalho.

```groovy title="genomics-1.nf" linenums="37"
output {
    bam_index {
        path '.'
    }
}
```

Cada alvo nomeado da seção `publish:` (como `bam_index`) recebe seu próprio bloco onde você pode configurar o caminho de saída relativo ao diretório de saída base.

!!! note "Nota"

    Embora os arquivos de dados que estamos usando aqui sejam muito pequenos, em genômica eles podem ficar muito grandes.
    Por padrão, o Nextflow cria links simbólicos para os arquivos de saída no diretório de publicação, o que evita cópias de arquivos desnecessárias.
    Você pode mudar esse comportamento usando a opção `mode` (por exemplo, `mode 'copy'`) para criar cópias reais.
    Esteja ciente de que os links simbólicos quebrarão quando você limpar seu diretório `work`, então para fluxos de trabalho de produção você pode querer usar `mode 'copy'`.

### 1.5. Configurar o diretório de saída

O diretório de saída base é definido via opção de configuração `outputDir`. Adicione-o ao `nextflow.config`:

=== "Depois"

    ```groovy title="nextflow.config" hl_lines="2"
    docker.enabled = true
    outputDir = 'results_genomics'
    ```

=== "Antes"

    ```groovy title="nextflow.config"
    docker.enabled = true
    ```

### 1.6. Executar o fluxo de trabalho para verificar que a etapa de indexação funciona

Vamos executar o fluxo de trabalho! Como lembrete, não precisamos especificar uma entrada na linha de comando porque configuramos um valor padrão para a entrada quando declaramos o parâmetro de entrada.

```bash
nextflow run genomics-1.nf
```

??? success "Saída do comando"

    ```console
    N E X T F L O W   ~  version 25.10.2

    ┃ Launching `genomics-1.nf` [reverent_sinoussi] DSL2 - revision: 41d43ad7fe

    executor >  local (1)
    [2a/e69536] SAMTOOLS_INDEX (1) | 1 of 1 ✔
    ```

Você pode verificar que o arquivo de índice foi gerado corretamente olhando no diretório de trabalho ou no diretório de resultados.

??? abstract "Conteúdo do diretório de trabalho"

    ```console
    work/2a/e695367b2f60df09cf826b07192dc3
    ├── reads_mother.bam -> /workspaces/training/nf4-science/genomics/data/bam/reads_mother.bam
    └── reads_mother.bam.bai
    ```

??? abstract "Conteúdo do diretório de resultados"

    ```console
    results_genomics/
    └── reads_mother.bam.bai -> */87/908bba*/reads_mother.bam.bai
    ```

Lá está!

### Conclusão

Você sabe como envolver uma ferramenta de genômica em um fluxo de trabalho Nextflow de etapa única e fazê-la executar usando um contêiner.

### O que vem a seguir?

Adicionar uma segunda etapa que consome a saída da primeira.

---

## 2. Adicionar um segundo processo para executar GATK HaplotypeCaller no arquivo BAM indexado

Agora que temos um índice para nosso arquivo de entrada, podemos passar para configurar a etapa de chamada de variantes, que é a parte interessante do fluxo de trabalho.

### 2.1. Definir o processo de chamada de variantes

Vamos escrever um processo, que chamaremos de `GATK_HAPLOTYPECALLER`, descrevendo a operação de chamada de variantes.

```groovy title="genomics-1.nf" linenums="44"
/*
 * Chamar variantes com GATK HaplotypeCaller
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

Você notará que introduzimos uma nova sintaxe aqui (`emit:`) para nomear exclusivamente cada um dos nossos canais de saída, e as razões para isso ficarão claras em breve.

Este comando recebe muito mais entradas, porque o GATK precisa de mais informações para realizar a análise em comparação com um trabalho simples de indexação.
Mas você notará que existem ainda mais entradas definidas no bloco de entradas do que as listadas no comando GATK. Por quê?

!!! note "Nota"

    O GATK sabe procurar o arquivo de índice BAM e os arquivos acessórios do genoma de referência porque está ciente das convenções em torno desses arquivos.
    No entanto, o Nextflow é projetado para ser independente de domínio e não sabe nada sobre requisitos de formato de arquivo de bioinformática.

Precisamos dizer ao Nextflow explicitamente que ele tem que preparar esses arquivos no diretório de trabalho em tempo de execução; caso contrário, ele não fará isso, e o GATK (corretamente) lançará um erro sobre os arquivos de índice estarem faltando.

Da mesma forma, temos que listar o arquivo de índice do VCF de saída (o arquivo `"${input_bam}.vcf.idx"`) explicitamente para que o Nextflow saiba rastrear esse arquivo caso seja necessário em etapas subsequentes.

### 2.2. Adicionar definições para entradas acessórias

Como nosso novo processo espera um punhado de arquivos adicionais serem fornecidos, configuramos alguns parâmetros CLI para eles na seção `Pipeline parameters`, junto com alguns valores padrão (mesmas razões de antes).

```groovy title="genomics-1.nf" linenums="8"
    // Arquivos acessórios
    reference: Path = "${projectDir}/data/ref/ref.fasta"
    reference_index: Path = "${projectDir}/data/ref/ref.fasta.fai"
    reference_dict: Path = "${projectDir}/data/ref/ref.dict"
    intervals: Path = "${projectDir}/data/ref/intervals.bed"
```

### 2.3. Criar variáveis para conter os caminhos dos arquivos acessórios

Enquanto entradas de dados principais são transmitidas dinamicamente através de canais, existem duas abordagens para lidar com arquivos acessórios. A abordagem recomendada é criar canais explícitos, o que torna o fluxo de dados mais claro e consistente. Alternativamente, a função file() para criar variáveis pode ser usada para casos mais simples, particularmente quando você precisa referenciar o mesmo arquivo em múltiplos processos - embora esteja ciente de que isso ainda cria canais implicitamente. <!-- TODO: Clarificar: isso ainda é necessário com entradas tipadas? -->

Adicione isso ao bloco workflow (após a criação do `reads_ch`, dentro da seção `main:`):

```groovy title="genomics-1.nf" linenums="79"
    // Carregar os caminhos de arquivo para os arquivos acessórios (referência e intervalos)
    ref_file        = file(params.reference)
    ref_index_file  = file(params.reference_index)
    ref_dict_file   = file(params.reference_dict)
    intervals_file  = file(params.intervals)
```

Isso tornará os caminhos dos arquivos acessórios disponíveis para fornecer como entrada a quaisquer processos que precisem deles.

### 2.4. Adicionar uma chamada ao bloco de fluxo de trabalho para executar GATK_HAPLOTYPECALLER

Agora que temos nosso segundo processo configurado e todas as entradas e arquivos acessórios estão prontos e disponíveis, podemos adicionar uma chamada ao processo `GATK_HAPLOTYPECALLER` no corpo do fluxo de trabalho.

```groovy title="genomics-1.nf" linenums="88"
    // Chamar variantes do arquivo BAM indexado
    GATK_HAPLOTYPECALLER(
        reads_ch,
        SAMTOOLS_INDEX.out,
        ref_file,
        ref_index_file,
        ref_dict_file,
        intervals_file
    )
```

Você deve reconhecer a sintaxe `*.out` da Parte 1 desta série de treinamento; estamos dizendo ao Nextflow para pegar a saída do canal por `SAMTOOLS_INDEX` e conectar isso na chamada do processo `GATK_HAPLOTYPECALLER`.

!!! note "Nota"

    Você notará que as entradas são fornecidas exatamente na mesma ordem na chamada ao processo como estão listadas no bloco de entrada do processo.
    No Nextflow, as entradas são posicionais, significando que você _deve_ seguir a mesma ordem; e é claro que tem que haver o mesmo número de elementos.

### 2.5. Atualizar a seção publish e o bloco de saída

Precisamos atualizar a seção `publish:` para incluir as saídas VCF, e adicionar alvos correspondentes no bloco `output`.

```groovy title="genomics-1.nf" linenums="99"
    publish:
    bam_index = SAMTOOLS_INDEX.out
    vcf = GATK_HAPLOTYPECALLER.out.vcf
    vcf_idx = GATK_HAPLOTYPECALLER.out.idx
}

output {
    bam_index {
        path '.'
    }
    vcf {
        path '.'
    }
    vcf_idx {
        path '.'
    }
}
```

### 2.6. Executar o fluxo de trabalho para verificar que a etapa de chamada de variantes funciona

Vamos executar o fluxo de trabalho expandido com `-resume` para que não tenhamos que executar a etapa de indexação novamente.

```bash
nextflow run genomics-1.nf -resume
```

??? success "Saída do comando"

    ```console
    N E X T F L O W   ~  version 25.10.2

    ┃ Launching `genomics-1.nf` [grave_volta] DSL2 - revision: 4790abc96a

    executor >  local (1)
    [2a/e69536] SAMTOOLS_INDEX (1)       | 1 of 1, cached: 1 ✔
    [53/e18e98] GATK_HAPLOTYPECALLER (1) | 1 of 1 ✔
    ```

Agora, se olharmos para a saída do console, vemos os dois processos listados.

O primeiro processo foi pulado graças ao cache, como esperado, enquanto o segundo processo foi executado por ser novo.

Você encontrará os arquivos de saída no diretório de resultados (como links simbólicos para o diretório de trabalho).

??? abstract "Conteúdo do diretório"

    ```console
    results_genomics/
    ├── reads_mother.bam.bai -> */87/908bba*/reads_mother.bam.bai
    ├── reads_mother.bam.vcf -> */cf/36f756*/reads_mother.bam.vcf
    └── reads_mother.bam.vcf.idx -> */cf/36f756*/reads_mother.bam.vcf.idx
    ```

Se você abrir o arquivo VCF, deve ver o mesmo conteúdo do arquivo que você gerou executando o comando GATK diretamente no contêiner.

```console title="reads_mother.bam.vcf" linenums="26"
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	reads_mother
20_10037292_10066351	3480	.	C	CT	503.03	.	AC=2;AF=1.00;AN=2;DP=23;ExcessHet=0.0000;FS=0.000;MLEAC=2;MLEAF=1.00;MQ=60.00;QD=27.95;SOR=1.179	GT:AD:DP:GQ:PL	1/1:0,18:18:54:517,54,0
20_10037292_10066351	3520	.	AT	A	609.03	.	AC=2;AF=1.00;AN=2;DP=18;ExcessHet=0.0000;FS=0.000;MLEAC=2;MLEAF=1.00;MQ=60.00;QD=33.83;SOR=0.693	GT:AD:DP:GQ:PL	1/1:0,18:18:54:623,54,0
20_10037292_10066351	3529	.	T	A	155.64	.	AC=1;AF=0.500;AN=2;BaseQRankSum=-0.544;DP=21;ExcessHet=0.0000;FS=1.871;MLEAC=1;MLEAF=0.500;MQ=60.00;MQRankSum=0.000;QD=7.78;ReadPosRankSum=-1.158;SOR=1.034	GT:AD:DP:GQ:PL	0/1:12,8:20:99:163,0,328
```

Esta é a saída que nos importa gerar para cada amostra em nosso estudo.

### Conclusão

Você sabe como fazer um fluxo de trabalho básico de duas etapas que faz trabalho de análise real e é capaz de lidar com idiossincrasias de formato de arquivo de genômica como os arquivos acessórios.

### O que vem a seguir?

Fazer o fluxo de trabalho lidar com múltiplas amostras em lote.

---

## 3. Adaptar o fluxo de trabalho para executar em um lote de amostras

É muito bom ter um fluxo de trabalho que pode automatizar o processamento em uma única amostra, mas e se você tiver 1000 amostras?
Você precisa escrever um script bash que percorre todas as suas amostras?

Não, graças aos céus! Basta fazer um pequeno ajuste no código e o Nextflow cuidará disso para você também.

### 3.1. Transformar a declaração de parâmetro de entrada em um array listando as três amostras

Vamos transformar aquele caminho de arquivo padrão na declaração do arquivo BAM de entrada em um array listando caminhos de arquivo para nossas três amostras de teste, na seção `Pipeline parameters`.

=== "Depois"

    ```groovy title="genomics-1.nf" linenums="7"
    // Entrada principal (array de três amostras)
        reads_bam = [
            "${projectDir}/data/bam/reads_mother.bam",
            "${projectDir}/data/bam/reads_father.bam",
            "${projectDir}/data/bam/reads_son.bam"
        ]
    ```

=== "Antes"

    ```groovy title="genomics-1.nf" linenums="7"
        // Entrada principal
        reads_bam: Path = "${projectDir}/data/bam/reads_mother.bam"
    ```

!!! note "Nota"

    Ao usar declarações de parâmetros tipadas (como `reads_bam: Path`), você não pode atribuir um valor de array.
    Para arrays, omita a anotação de tipo.

E isso é realmente tudo que precisamos fazer, porque a factory de canal que usamos no corpo do fluxo de trabalho (`.fromPath`) está tão feliz em aceitar múltiplos caminhos de arquivo para carregar no canal de entrada quanto estava em carregar um único.

!!! note "Nota"

    Normalmente, você não iria querer codificar a lista de amostras no seu arquivo de fluxo de trabalho, mas estamos fazendo isso aqui para manter as coisas simples.
    Apresentaremos maneiras mais elegantes de lidar com entradas mais tarde nesta série de treinamento.

### 3.2. Executar o fluxo de trabalho para verificar que ele executa em todas as três amostras

Vamos tentar executar o fluxo de trabalho agora que o encanamento está configurado para executar em todas as três amostras de teste.

```bash
nextflow run genomics-1.nf -resume
```

Coisa engraçada: isso _pode funcionar_, OU _pode falhar_. Por exemplo, aqui está uma execução que teve sucesso:

??? success "Saída do comando"

    ```console
    N E X T F L O W   ~  version 25.10.2

    ┃ Launching `genomics-1.nf` [peaceful_yalow] DSL2 - revision: a256d113ad

    executor >  local (6)
    [4f/7071b0] SAMTOOLS_INDEX (3)       | 3 of 3, cached: 1 ✔
    [7a/89bc43] GATK_HAPLOTYPECALLER (2) | 3 of 3, cached: 1 ✔
    ```

Se sua execução de fluxo de trabalho teve sucesso, execute-a novamente até obter um erro como este:

??? failure "Saída do comando"

    ```console
    N E X T F L O W   ~  version 25.10.2

    ┃ Launching `genomics-1.nf` [loving_pasteur] DSL2 - revision: d2a8e63076

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

Bem, isso é estranho, considerando que indexamos explicitamente os arquivos BAM na primeira etapa do fluxo de trabalho. Poderia haver algo errado com o encanamento?

#### 3.2.1. Verificar os diretórios de trabalho para as chamadas relevantes

Vamos dar uma olhada dentro do diretório de trabalho para a chamada de processo `GATK_HAPLOTYPECALLER` com falha listada na saída do console.

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

Preste atenção particular aos nomes do arquivo BAM e do índice BAM que estão listados neste diretório: `reads_son.bam` e `reads_father.bam.bai`.

Mas o quê? O Nextflow preparou um arquivo de índice no diretório de trabalho desta chamada de processo, mas é o errado. Como isso pôde acontecer?

#### 3.2.2. Usar o [operador view()](https://www.nextflow.io/docs/latest/reference/operator.html#view) para inspecionar o conteúdo do canal

Adicione estas duas linhas no corpo do fluxo de trabalho antes da chamada do processo `GATK_HAPLOTYPER`:

```groovy title="genomics-1.nf" linenums="84"
    // diagnósticos temporários
    reads_ch.view()
    SAMTOOLS_INDEX.out.view()
```

Então execute o comando do fluxo de trabalho novamente.

```bash
nextflow run genomics-1.nf
```

Mais uma vez, isso pode ter sucesso ou falhar. Aqui está uma execução bem-sucedida:

??? success "Saída do comando"

    ```console
    N E X T F L O W   ~  version 25.10.2

    ┃ Launching `genomics-1.nf` [fervent_pasteur] DSL2 - revision: a256d113ad

    /workspaces/training/nf4-science/genomics/data/bam/reads_mother.bam
    /workspaces/training/nf4-science/genomics/data/bam/reads_father.bam
    /workspaces/training/nf4-science/genomics/data/bam/reads_son.bam
    executor >  local (6)
    [4f/7071b0] SAMTOOLS_INDEX (3)       | 3 of 3 ✔
    /workspaces/training/nf4-science/genomics/work/b4/45a376f0e724be1dc626a6807f73d8/reads_mother.bam.bai
    /workspaces/training/nf4-science/genomics/work/4f/7071b082b45dd85b1c9b6b3b32cb69/reads_father.bam.bai
    /workspaces/training/nf4-science/genomics/work/3c/331645a9e20e67edae10da5ba17c7b/reads_son.bam.bai
    [a2/dbd8d5] GATK_HAPLOTYPECALLER (3) | 3 of 3 ✔
    ```

E aqui está uma com falha:

??? failure "Saída do comando"

    ```console
    N E X T F L O W   ~  version 25.10.2

    ┃ Launching `genomics-1.nf` [angry_hamilton] DSL2 - revision: a256d113ad

    /workspaces/training/nf4-science/genomics/data/bam/reads_mother.bam
    /workspaces/training/nf4-science/genomics/data/bam/reads_father.bam
    /workspaces/training/nf4-science/genomics/data/bam/reads_son.bam
    executor >  local (6)
    [4f/7071b0] SAMTOOLS_INDEX (3)       | 3 of 3 ✔
    /workspaces/training/nf4-science/genomics/work/4f/7071b082b45dd85b1c9b6b3b32cb69/reads_father.bam.bai
    /workspaces/training/nf4-science/genomics/work/3c/331645a9e20e67edae10da5ba17c7b/reads_son.bam.bai
    /workspaces/training/nf4-science/genomics/work/b4/45a376f0e724be1dc626a6807f73d8/reads_mother.bam.bai
    [a3/cf3a89] GATK_HAPLOTYPECALLER (3) | 1 of 3
    ERROR ~ Error executing process > 'GATK_HAPLOTYPECALLER (3)'
    ...
    ```

Você pode precisar executá-lo várias vezes para que falhe novamente.
Este erro não se reproduzirá consistentemente porque depende de alguma variabilidade nos tempos de execução das chamadas de processo individuais.

Isso é o que a saída das duas chamadas `.view()` que adicionamos parece para uma execução com falha:

```console
/workspaces/training/nf4-science/genomics/data/bam/reads_mother.bam
/workspaces/training/nf4-science/genomics/data/bam/reads_father.bam
/workspaces/training/nf4-science/genomics/data/bam/reads_son.bam
/workspaces/training/nf4-science/genomics/work/9c/53492e3518447b75363e1cd951be4b/reads_father.bam.bai
/workspaces/training/nf4-science/genomics/work/cc/37894fffdf6cc84c3b0b47f9b536b7/reads_son.bam.bai
/workspaces/training/nf4-science/genomics/work/4d/dff681a3d137ba7d9866e3d9307bd0/reads_mother.bam.bai
```

As primeiras três linhas correspondem ao canal de entrada e a segunda, ao canal de saída.
Você pode ver que os arquivos BAM e arquivos de índice para as três amostras não estão listados na mesma ordem!

!!! note "Nota"

    Quando você chama um processo Nextflow em um canal contendo múltiplos elementos, o Nextflow tentará paralelizar a execução o máximo possível, e coletará saídas em qualquer ordem que elas se tornem disponíveis.
    A consequência é que as saídas correspondentes podem ser coletadas em uma ordem diferente da que as entradas originais foram alimentadas.

Como escrito atualmente, nosso script de fluxo de trabalho assume que os arquivos de índice sairão da etapa de indexação listados na mesma ordem mãe/pai/filho que as entradas foram dadas.
Mas isso não é garantido ser o caso, que é por que às vezes (embora nem sempre) os arquivos errados são emparelhados na segunda etapa.

Para corrigir isso, precisamos garantir que os arquivos BAM e seus arquivos de índice viajem juntos através dos canais.

!!! tip "Dica"

    As instruções `view()` no código do fluxo de trabalho não fazem nada, então não é um problema deixá-las.
    No entanto, elas vão poluir sua saída do console, então recomendamos removê-las quando você terminar de resolver o problema.

### 3.3. Mudar a saída do processo SAMTOOLS_INDEX para uma tupla que mantém o arquivo de entrada e seu índice juntos

A maneira mais simples de garantir que um arquivo BAM e seu índice permaneçam intimamente associados é empacotá-los juntos em uma tupla saindo da tarefa de índice.

!!! note "Nota"

    Uma **tupla** é uma lista finita e ordenada de elementos que é comumente usada para retornar múltiplos valores de uma função. Tuplas são particularmente úteis para passar múltiplas entradas ou saídas entre processos enquanto preservam sua associação e ordem.

Primeiro, vamos mudar a saída do processo `SAMTOOLS_INDEX` para incluir o arquivo BAM em sua declaração de saída.

=== "Depois"

    ```groovy title="genomics-1.nf" linenums="20"
        output:
        tuple path(input_bam), path("${input_bam}.bai")
    ```

=== "Antes"

    ```groovy title="genomics-1.nf" linenums="20"
        output:
        path "${input_bam}.bai"
    ```

Dessa forma, cada arquivo de índice será estreitamente acoplado com seu arquivo BAM original, e a saída geral da etapa de indexação será um único canal contendo pares de arquivos.

### 3.4. Mudar a entrada para o processo GATK_HAPLOTYPECALLER para ser uma tupla

Como mudamos a 'forma' da saída do primeiro processo no fluxo de trabalho, precisamos atualizar a definição de entrada do segundo processo para corresponder.

Especificamente, onde anteriormente declaramos dois caminhos de entrada separados no bloco de entrada do processo `GATK_HAPLOTYPECALLER`, agora declaramos uma única entrada correspondendo à estrutura da tupla emitida por `SAMTOOLS_INDEX`.

=== "Depois"

    ```groovy title="genomics-1.nf" linenums="51"
        input:
        tuple path(input_bam), path(input_bam_index)
    ```

=== "Antes"

    ```groovy title="genomics-1.nf" linenums="51"
        input:
        path input_bam
        path input_bam_index
    ```

É claro, já que mudamos a forma das entradas que `GATK_HAPLOTYPECALLER` espera, precisamos atualizar a chamada do processo de acordo no corpo do fluxo de trabalho.

### 3.5. Atualizar a chamada para GATK_HAPLOTYPECALLER no bloco de fluxo de trabalho

Não precisamos mais fornecer o `reads_ch` original ao processo `GATK_HAPLOTYPECALLER`, já que o arquivo BAM agora está empacotado na saída do canal por `SAMTOOLS_INDEX`.

Como resultado, podemos simplesmente deletar aquela linha.

=== "Depois"

    ```groovy title="genomics-1.nf" linenums="88"
    GATK_HAPLOTYPECALLER(
        SAMTOOLS_INDEX.out,
    ```

=== "Antes"

    ```groovy title="genomics-1.nf" linenums="88"
    GATK_HAPLOTYPECALLER(
        reads_ch,
        SAMTOOLS_INDEX.out,
    ```

Essa é toda a reconexão necessária para resolver o problema de incompatibilidade de índice.

### 3.6. Atualizar a seção publish e o bloco de saída para a tupla

Como `SAMTOOLS_INDEX.out` agora é uma tupla contendo tanto o BAM quanto seu índice, ambos os arquivos serão publicados juntos.
Renomeamos o alvo de `bam_index` para `indexed_bam` para refletir que agora contém ambos os arquivos.

=== "Depois"

    ```groovy title="genomics-1.nf" hl_lines="2"
        publish:
        indexed_bam = SAMTOOLS_INDEX.out
    ```

=== "Antes"

    ```groovy title="genomics-1.nf"
        publish:
        bam_index = SAMTOOLS_INDEX.out
    ```

Também precisamos atualizar o bloco output para usar o novo nome de alvo:

=== "Depois"

    ```groovy title="genomics-1.nf" hl_lines="2"
    output {
        indexed_bam {
            path '.'
        }
    ```

=== "Antes"

    ```groovy title="genomics-1.nf"
    output {
        bam_index {
            path '.'
        }
    ```

### 3.7. Executar o fluxo de trabalho para verificar que funciona corretamente em todas as três amostras sempre

É claro, a prova está no pudim, então vamos executar o fluxo de trabalho novamente algumas vezes para garantir que isso funcionará de forma confiável daqui para frente.

```bash
nextflow run genomics-1.nf
```

Desta vez (e todas as vezes) tudo deve executar corretamente:

??? success "Saída do comando"

    ```console
    N E X T F L O W   ~  version 25.10.2

    ┃ Launching `genomics-1.nf` [special_goldstine] DSL2 - revision: 4cbbf6ea3e

    executor >  local (6)
    [d6/10c2c4] SAMTOOLS_INDEX (1)       | 3 of 3 ✔
    [88/1783aa] GATK_HAPLOTYPECALLER (2) | 3 of 3 ✔
    ```

O diretório de resultados agora contém tanto arquivos BAM quanto BAI para cada amostra (da tupla), junto com as saídas VCF:

??? abstract "Conteúdo do diretório de resultados"

    ```console
    results_genomics/
    ├── reads_father.bam -> */60/e2614c*/reads_father.bam
    ├── reads_father.bam.bai -> */60/e2614c*/reads_father.bam.bai
    ├── reads_father.bam.vcf -> */b8/91b3c8*/reads_father.bam.vcf
    ├── reads_father.bam.vcf.idx -> */b8/91b3c8*/reads_father.bam.vcf.idx
    ├── reads_mother.bam -> */3e/fededc*/reads_mother.bam
    ├── reads_mother.bam.bai -> */3e/fededc*/reads_mother.bam.bai
    ├── reads_mother.bam.vcf -> */32/5ca037*/reads_mother.bam.vcf
    ├── reads_mother.bam.vcf.idx -> */32/5ca037*/reads_mother.bam.vcf.idx
    ├── reads_son.bam -> */3c/36d1c2*/reads_son.bam
    ├── reads_son.bam.bai -> */3c/36d1c2*/reads_son.bam.bai
    ├── reads_son.bam.vcf -> */d7/a6b046*/reads_son.bam.vcf
    └── reads_son.bam.vcf.idx -> */d7/a6b046*/reads_son.bam.vcf.idx
    ```

Se você quiser, pode usar `.view()` novamente para espiar como o conteúdo do canal de saída `SAMTOOLS_INDEX` se parece:

```groovy title="genomics-1.nf" linenums="92"
SAMTOOLS_INDEX.out.view()
```

Você verá que o canal contém as três tuplas esperadas (caminhos de arquivo truncados para legibilidade).

```console title="Saída"
[*/60/e2614c*/reads_father.bam, */60/e2614c*/reads_father.bam.bai]
[*/3e/fededc*/reads_mother.bam, */3e/fededc*/reads_mother.bam.bai]
[*/3c/36d1c2*/reads_son.bam, */3c/36d1c2*/reads_son.bam.bai]
```

Isso será muito mais seguro, daqui para frente.

### Conclusão

Você sabe como fazer seu fluxo de trabalho executar em múltiplas amostras (independentemente).

### O que vem a seguir?

Facilitar o manuseio de amostras em lote.

---

## 4. Fazer o fluxo de trabalho aceitar um arquivo de texto contendo um lote de arquivos de entrada

Uma maneira muito comum de fornecer múltiplos arquivos de dados de entrada para um fluxo de trabalho é fazê-lo com um arquivo de texto contendo os caminhos dos arquivos.
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

!!! note "Nota"

    Os arquivos que estamos usando aqui estão apenas no sistema de arquivos local do seu GitHub Codespaces, mas também poderíamos apontar para arquivos em armazenamento na nuvem.

### 4.2. Atualizar o padrão do parâmetro

Vamos mudar o valor padrão para nosso parâmetro de entrada `reads_bam` para apontar para o arquivo `sample_bams.txt`.

=== "Depois"

    ```groovy title="genomics-1.nf" linenums="7"
        // Entrada principal (arquivo de arquivos de entrada, um por linha)
        reads_bam: Path = "${projectDir}/data/sample_bams.txt"
    ```

=== "Antes"

    ```groovy title="genomics-1.nf" linenums="7"
    // Entrada principal (array de três amostras)
        reads_bam = [
            "${projectDir}/data/bam/reads_mother.bam",
            "${projectDir}/data/bam/reads_father.bam",
            "${projectDir}/data/bam/reads_son.bam"
        ]
    ```

Dessa forma podemos continuar a ser preguiçosos, mas a lista de arquivos não vive mais no próprio código do fluxo de trabalho, o que é um grande passo na direção certa.

### 4.3. Atualizar a factory de canal para ler linhas de um arquivo

Atualmente, nossa factory de canal de entrada trata quaisquer arquivos que demos a ela como as entradas de dados que queremos alimentar ao processo de indexação.
Como agora estamos dando a ela um arquivo que lista caminhos de arquivo de entrada, precisamos mudar seu comportamento para analisar o arquivo e tratar os caminhos de arquivo que ele contém como as entradas de dados.

Felizmente podemos fazer isso muito simplesmente, apenas adicionando o [operador `.splitText()`](https://www.nextflow.io/docs/latest/reference/operator.html#operator-splittext) à etapa de construção do canal.

=== "Depois"

    ```groovy title="genomics-1.nf" linenums="68"
        // Criar canal de entrada de um arquivo de texto listando caminhos de arquivo de entrada
        reads_ch = channel.fromPath(params.reads_bam).splitText()
    ```

=== "Antes"

    ```groovy title="genomics-1.nf" linenums="68"
        // Criar canal de entrada (arquivo único via parâmetro CLI)
        reads_ch = channel.fromPath(params.reads_bam)
    ```

!!! tip "Dica"

    Esta é outra ótima oportunidade para usar o operador `.view()` para ver como o conteúdo do canal se parece antes e depois de aplicar um operador.

### 4.4. Executar o fluxo de trabalho para verificar que funciona corretamente

Vamos executar o fluxo de trabalho mais uma vez. Isso deve produzir o mesmo resultado de antes, certo?

```bash
nextflow run genomics-1.nf -resume
```

??? success "Saída do comando"

    ```console
    N E X T F L O W   ~  version 25.10.2

    ┃ Launching `genomics-1.nf` [sick_albattani] DSL2 - revision: 46d84642f6

    [18/23b4bb] SAMTOOLS_INDEX (1)       | 3 of 3, cached: 3 ✔
    [12/f727bb] GATK_HAPLOTYPECALLER (3) | 3 of 3, cached: 3 ✔
    ```

Sim! Na verdade, o Nextflow detecta corretamente que as chamadas de processo são exatamente as mesmas, e nem se preocupa em executar tudo novamente, já que estávamos executando com `-resume`.

E é isso! Nosso fluxo de trabalho simples de chamada de variantes tem todas as características básicas que queríamos.

### Conclusão

Você sabe como fazer um fluxo de trabalho linear de múltiplas etapas para indexar um arquivo BAM e aplicar chamada de variantes por amostra usando GATK.

De forma mais geral, você aprendeu como usar componentes e lógica essenciais do Nextflow para construir um pipeline de genômica simples que faz trabalho real, levando em conta as idiossincrasias dos formatos de arquivo de genômica e requisitos de ferramentas.

### O que vem a seguir?

Celebre seu sucesso e faça uma pausa extra longa!

Na próxima parte deste curso, você aprenderá como usar alguns recursos adicionais do Nextflow (incluindo mais operadores de canal) para aplicar chamada conjunta de variantes aos dados.
