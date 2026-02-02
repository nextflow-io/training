# Parte 2: Chamada conjunta em uma coorte

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Tradução assistida por IA - [saiba mais e sugira melhorias](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

Na primeira parte deste curso, você construiu um pipeline de chamada de variantes que era completamente linear e processava os dados de cada amostra independentemente das outras.
No entanto, em um caso de uso real de genômica, você normalmente precisará examinar as chamadas de variantes de múltiplas amostras em conjunto.

Nesta segunda parte, mostramos como usar canais e operadores de canal para implementar a chamada conjunta de variantes com GATK, construindo sobre o pipeline da Parte 1.

### Visão geral do método

O método de chamada de variantes do GATK que usamos na primeira parte deste curso simplesmente gerou chamadas de variantes por amostra.
Isso é adequado se você quiser apenas examinar as variantes de cada amostra isoladamente, mas isso fornece informações limitadas.
Frequentemente é mais interessante examinar como as chamadas de variantes diferem entre múltiplas amostras, e para fazer isso, o GATK oferece um método alternativo chamado chamada conjunta de variantes, que demonstramos aqui.

A chamada conjunta de variantes envolve gerar um tipo especial de saída de variantes chamado GVCF (Genomic VCF) para cada amostra, depois combinar os dados GVCF de todas as amostras e, finalmente, executar uma análise estatística de 'genotipagem conjunta'.

![Análise conjunta](img/joint-calling.png)

O que é especial sobre o GVCF de uma amostra é que ele contém registros resumindo estatísticas de dados de sequência sobre todas as posições na área alvo do genoma, não apenas as posições onde o programa encontrou evidências de variação.
Isso é crítico para o cálculo de genotipagem conjunta ([leitura adicional](https://gatk.broadinstitute.org/hc/en-us/articles/360035890431-The-logic-of-joint-calling-for-germline-short-variants)).

O GVCF é produzido pelo GATK HaplotypeCaller, a mesma ferramenta que usamos na Parte 1, com um parâmetro adicional (`-ERC GVCF`).
A combinação dos GVCFs é feita com GATK GenomicsDBImport, que combina as chamadas por amostra em um armazenamento de dados (análogo a um banco de dados), então a análise de 'genotipagem conjunta' propriamente dita é feita com GATK GenotypeGVCFs.

### Fluxo de trabalho

Então, para recapitular, nesta parte do curso, vamos desenvolver um fluxo de trabalho que faz o seguinte:

<figure class="excalidraw">
--8<-- "docs/en/docs/nf4_science/genomics/img/hello-gatk-2.svg"
</figure>

1. Gerar um arquivo de índice para cada arquivo BAM de entrada usando Samtools
2. Executar o GATK HaplotypeCaller em cada arquivo BAM de entrada para gerar um GVCF de chamadas de variantes genômicas por amostra
3. Coletar todos os GVCFs e combiná-los em um armazenamento de dados GenomicsDB
4. Executar genotipagem conjunta nos dados GVCF combinados para produzir um VCF em nível de coorte

Vamos aplicar isso ao mesmo conjunto de dados da Parte 1.

---

## 0. Aquecimento: Execute Samtools e GATK diretamente

Assim como anteriormente, queremos testar os comandos manualmente antes de tentar envolvê-los em um fluxo de trabalho.

!!! note

     Certifique-se de estar no diretório de trabalho correto:
     `cd /workspaces/training/nf4-science/genomics`

### 0.1. Indexar um arquivo BAM de entrada com Samtools

Este primeiro passo é o mesmo da Parte 1, então deve parecer muito familiar, mas desta vez precisamos fazer isso para todas as três amostras.

!!! note

    Tecnicamente já geramos arquivos de índice para as três amostras através do nosso pipeline, então poderíamos ir buscá-los no diretório de resultados. No entanto, é mais limpo apenas refazer isso manualmente, e levará apenas um minuto.

#### 0.1.1. Iniciar o contêiner Samtools interativamente

```bash
docker run -it -v ./data:/data community.wave.seqera.io/library/samtools:1.20--b5dfbd93de237464
```

<!--
??? success "Saída do comando"

    ```console

    ```
-->

#### 0.1.2. Executar o comando de indexação para as três amostras

```bash
samtools index /data/bam/reads_mother.bam
samtools index /data/bam/reads_father.bam
samtools index /data/bam/reads_son.bam
```

Assim como anteriormente, isso deve produzir os arquivos de índice no mesmo diretório que os arquivos BAM correspondentes.

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

Agora que temos arquivos de índice para todas as três amostras, podemos prosseguir para gerar os GVCFs para cada uma delas.

#### 0.1.3. Sair do contêiner Samtools

```bash
exit
```

### 0.2. Chamar variantes com GATK HaplotypeCaller em modo GVCF

Este segundo passo é muito similar ao que fizemos na Parte 1: Hello Genomics, mas agora vamos executar o GATK em 'modo GVCF'.

#### 0.2.1. Iniciar o contêiner GATK interativamente

```bash
docker run -it -v ./data:/data community.wave.seqera.io/library/gatk4:4.5.0.0--730ee8817e436867
```

<!--
??? success "Saída do comando"

    ```console

    ```
-->

#### 0.2.2. Executar o comando de chamada de variantes com a opção GVCF

Para produzir um VCF genômico (GVCF), adicionamos a opção `-ERC GVCF` ao comando base, que ativa o modo GVCF do HaplotypeCaller.

Também mudamos a extensão do arquivo de saída de `.vcf` para `.g.vcf`.
Tecnicamente isso não é um requisito, mas é uma convenção fortemente recomendada.

```bash
gatk HaplotypeCaller \
        -R /data/ref/ref.fasta \
        -I /data/bam/reads_mother.bam \
        -O reads_mother.g.vcf \
        -L /data/ref/intervals.bed \
        -ERC GVCF
```

<!--
??? success "Saída do comando"

    ```console

    ```
-->

Isso cria o arquivo de saída GVCF `reads_mother.g.vcf` no diretório de trabalho atual no contêiner.

Se você usar `cat` para visualizar o conteúdo, verá que é muito mais longo do que o VCF equivalente que geramos na Parte 1. Você nem consegue rolar até o início do arquivo, e a maioria das linhas parece bem diferente do que vimos no VCF na Parte 1.

```console title="Saída" linenums="1674"
20_10037292_10066351    14714   .       T       <NON_REF>       .       .       END=14718       GT:DP:GQ:MIN_DP:PL       0/0:37:99:37:0,99,1192
20_10037292_10066351    14719   .       T       <NON_REF>       .       .       END=14719       GT:DP:GQ:MIN_DP:PL       0/0:36:82:36:0,82,1087
20_10037292_10066351    14720   .       T       <NON_REF>       .       .       END=14737       GT:DP:GQ:MIN_DP:PL       0/0:42:99:37:0,100,1160
```

Estas representam regiões não variantes onde o chamador de variantes não encontrou evidências de variação, então ele capturou algumas estatísticas descrevendo seu nível de confiança na ausência de variação. Isso torna possível distinguir entre dois casos muito diferentes: (1) há dados de boa qualidade mostrando que a amostra é homozigota-referência, e (2) não há dados bons suficientes disponíveis para fazer uma determinação de qualquer maneira.

Em um GVCF, normalmente há muitas dessas linhas não variantes, com um número menor de registros variantes espalhados entre elas. Tente executar `head -176` no GVCF para carregar apenas as primeiras 176 linhas do arquivo para encontrar uma chamada de variante real.

```console title="Saída" linenums="174"
20_10037292_10066351    3479    .       T       <NON_REF>       .       .       END=3479        GT:DP:GQ:MIN_DP:PL       0/0:34:36:34:0,36,906
20_10037292_10066351    3480    .       C       CT,<NON_REF>    503.03  .       DP=23;ExcessHet=0.0000;MLEAC=2,0;MLEAF=1.00,0.00;RAW_MQandDP=82800,23    GT:AD:DP:GQ:PL:SB       1/1:0,18,0:18:54:517,54,0,517,54,517:0,0,7,11
20_10037292_10066351    3481    .       T       <NON_REF>       .       .       END=3481        GT:DP:GQ:MIN_DP:PL       0/0:21:51:21:0,51,765
```

A segunda linha mostra o primeiro registro variante no arquivo, que corresponde à primeira variante no arquivo VCF que examinamos na Parte 1.

Assim como o VCF original, o arquivo de saída GVCF também é acompanhado por um arquivo de índice, chamado `reads_mother.g.vcf.idx`.

#### 0.2.3. Repetir o processo nas outras duas amostras

Para testar a etapa de genotipagem conjunta, precisamos de GVCFs para todas as três amostras, então vamos gerar esses manualmente agora.

```bash
gatk HaplotypeCaller \
        -R /data/ref/ref.fasta \
        -I /data/bam/reads_father.bam \
        -O reads_father.g.vcf \
        -L /data/ref/intervals.bed \
        -ERC GVCF
```

<!--
??? success "Saída do comando"

    ```console

    ```
-->

```bash
gatk HaplotypeCaller \
        -R /data/ref/ref.fasta \
        -I /data/bam/reads_son.bam \
        -O reads_son.g.vcf \
        -L /data/ref/intervals.bed \
        -ERC GVCF
```

<!--
??? success "Saída do comando"

    ```console

    ```
-->

Uma vez que isso seja concluído, você deve ter três arquivos terminando em `.g.vcf` no seu diretório atual (um por amostra) e seus respectivos arquivos de índice terminando em `.g.vcf.idx`.

### 0.3. Executar genotipagem conjunta

Agora que temos todos os GVCFs, podemos finalmente testar a abordagem de genotipagem conjunta para gerar chamadas de variantes para uma coorte de amostras.
Como lembrete, é um método de duas etapas que consiste em combinar os dados de todos os GVCFs em um armazenamento de dados, e então executar a análise de genotipagem conjunta propriamente dita para gerar o VCF final de variantes chamadas conjuntamente.

#### 0.3.1. Combinar todos os GVCFs por amostra

Esta primeira etapa usa outra ferramenta GATK, chamada GenomicsDBImport, para combinar os dados de todos os GVCFs em um armazenamento de dados GenomicsDB.

```bash
gatk GenomicsDBImport \
    -V reads_mother.g.vcf \
    -V reads_father.g.vcf \
    -V reads_son.g.vcf \
    -L /data/ref/intervals.bed \
    --genomicsdb-workspace-path family_trio_gdb
```

<!--
??? success "Saída do comando"

    ```console

    ```
-->

A saída desta etapa é efetivamente um diretório contendo um conjunto de diretórios ainda mais aninhados contendo os dados de variantes combinados na forma de múltiplos arquivos diferentes.
Você pode explorar, mas rapidamente verá que este formato de armazenamento de dados não é destinado a ser lido diretamente por humanos.

!!! note

    GATK inclui ferramentas que tornam possível inspecionar e extrair dados de chamadas de variantes do armazenamento de dados conforme necessário.

#### 0.3.2. Executar a análise de genotipagem conjunta propriamente dita

Esta segunda etapa usa ainda outra ferramenta GATK, chamada GenotypeGVCFs, para recalcular estatísticas de variantes e genótipos individuais à luz dos dados disponíveis em todas as amostras da coorte.

```bash
gatk GenotypeGVCFs \
    -R /data/ref/ref.fasta \
    -V gendb://family_trio_gdb \
    -O family_trio.vcf
```

<!--
??? success "Saída do comando"

    ```console

    ```
-->

Isso cria o arquivo de saída VCF `family_trio.vcf` no diretório de trabalho atual no contêiner.
É outro arquivo razoavelmente pequeno, então você pode usar `cat` neste arquivo para ver seu conteúdo, e rolar para cima para encontrar as primeiras linhas de variantes.

```console title="family_trio.vcf" linenums="40"
#CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO    FORMAT  reads_father    reads_mother    reads_son
20_10037292_10066351    3480    .       C       CT      1625.89 .       AC=5;AF=0.833;AN=6;BaseQRankSum=0.220;DP=85;ExcessHet=0.0000;FS=2.476;MLEAC=5;MLEAF=0.833;MQ=60.00;MQRankSum=0.00;QD=21.68;ReadPosRankSum=-1.147e+00;SOR=0.487    GT:AD:DP:GQ:PL  0/1:15,16:31:99:367,0,375       1/1:0,18:18:54:517,54,0 1/1:0,26:26:78:756,78,0
20_10037292_10066351    3520    .       AT      A       1678.89 .       AC=5;AF=0.833;AN=6;BaseQRankSum=1.03;DP=80;ExcessHet=0.0000;FS=2.290;MLEAC=5;MLEAF=0.833;MQ=60.00;MQRankSum=0.00;QD=22.39;ReadPosRankSum=0.701;SOR=0.730 GT:AD:DP:GQ:PL   0/1:18,13:31:99:296,0,424       1/1:0,18:18:54:623,54,0 1/1:0,26:26:78:774,78,0
20_10037292_10066351    3529    .       T       A       154.29  .       AC=1;AF=0.167;AN=6;BaseQRankSum=-5.440e-01;DP=104;ExcessHet=0.0000;FS=1.871;MLEAC=1;MLEAF=0.167;MQ=60.00;MQRankSum=0.00;QD=7.71;ReadPosRankSum=-1.158e+00;SOR=1.034       GT:AD:DP:GQ:PL  0/0:44,0:44:99:0,112,1347       0/1:12,8:20:99:163,0,328        0/0:39,0:39:99:0,105,1194
```

Isso se parece mais com o VCF original que geramos na Parte 1, exceto que desta vez temos informações em nível de genótipo para todas as três amostras.
As últimas três colunas no arquivo são os blocos de genótipos para as amostras, listadas em ordem alfabética.

Se olharmos para os genótipos chamados para nosso trio familiar de teste para a primeira variante, vemos que o pai é heterozigoto-variante (`0/1`), e a mãe e o filho são ambos homozigotos-variante (`1/1`).

Essa é, em última análise, a informação que estamos procurando extrair do conjunto de dados! Então vamos envolver tudo isso em um fluxo de trabalho Nextflow para que possamos fazer isso em escala.

#### 0.3.3. Sair do contêiner GATK

```bash
exit
```

### Conclusão

Você sabe como executar os comandos individuais envolvidos na chamada conjunta de variantes no terminal para verificar que eles produzirão as informações que você deseja.

### O que vem a seguir?

Envolver esses comandos em um pipeline real.

---

## 1. Modificar a etapa de chamada de variantes por amostra para produzir um GVCF

A boa notícia é que não precisamos começar tudo de novo, já que já escrevemos um fluxo de trabalho que faz parte desse trabalho na Parte 1.
No entanto, esse pipeline produz arquivos VCF, enquanto agora queremos arquivos GVCF para fazer a genotipagem conjunta.
Então precisamos começar ativando o modo de chamada de variantes GVCF e atualizando a extensão do arquivo de saída.

!!! note

    Por conveniência, vamos trabalhar com uma cópia nova do fluxo de trabalho GATK como está no final da Parte 1, mas sob um nome diferente: `genomics-2.nf`.

### 1.1. Dizer ao HaplotypeCaller para emitir um GVCF e atualizar a extensão de saída

Vamos abrir o arquivo `genomics-2.nf` no editor de código.
Ele deve parecer muito familiar, mas sinta-se à vontade para executá-lo se quiser se convencer de que ele funciona como esperado.

Vamos começar fazendo duas mudanças:

- Adicionar o parâmetro `-ERC GVCF` ao comando GATK HaplotypeCaller;
- Atualizar o caminho do arquivo de saída para usar a extensão correspondente `.g.vcf`, conforme convenção do GATK.

Certifique-se de adicionar uma barra invertida (`\`) no final da linha anterior quando adicionar `-ERC GVCF`.

=== "Depois"

    ```groovy title="genomics-2.nf" linenums="56" hl_lines="4 6"
        """
        gatk HaplotypeCaller \
            -R ${ref_fasta} \
            -I ${input_bam} \
            -O ${input_bam}.g.vcf \
            -L ${interval_list} \
            -ERC GVCF
        """
    ```

=== "Antes"

    ```groovy title="genomics-2.nf" linenums="56" hl_lines="4"
        """
        gatk HaplotypeCaller \
            -R ${ref_fasta} \
            -I ${input_bam} \
            -O ${input_bam}.vcf \
            -L ${interval_list}
        """
    ```

E isso é tudo que é necessário para fazer o HaplotypeCaller gerar GVCFs em vez de VCFs, certo?

### 1.2. Executar o pipeline para verificar que você pode gerar GVCFs

O comando de execução Nextflow é o mesmo de antes, exceto pelo próprio nome do arquivo de fluxo de trabalho.
Certifique-se de atualizar isso apropriadamente.

```bash
nextflow run genomics-2.nf
```

??? success "Saída do comando"

    ```console
    N E X T F L O W   ~  version 25.10.2

    ┃ Launching `genomics-2.nf` [crazy_venter] DSL2 - revision: a2d6f6f09f

    executor >  local (6)
    [f1/8d8486] SAMTOOLS_INDEX (1)       | 3 of 3 ✔
    [72/3249ca] GATK_HAPLOTYPECALLER (3) | 0 of 3
    ERROR ~ Error executing process > 'GATK_HAPLOTYPECALLER (2)'

    Caused by:
      Missing output file(s) `reads_son.bam.vcf` expected by process `GATK_HAPLOTYPECALLER (2)`

    Command executed:

      gatk HaplotypeCaller         -R ref.fasta         -I reads_son.bam         -O reads_son.bam.g.vcf         -L intervals.bed         -ERC GVCF
    ```

E a saída é... toda vermelha! Oh não.

O comando que foi executado está correto, então estávamos certos de que isso era suficiente para mudar o comportamento da ferramenta GATK.
Mas olhe aquela linha sobre o arquivo de saída ausente. Nota algo?

Isso mesmo, esquecemos de dizer ao Nextflow para esperar um novo nome de arquivo. Ops.

### 1.3. Atualizar a extensão do arquivo de saída no bloco de saídas do processo também

Porque não é suficiente apenas mudar a extensão do arquivo no próprio comando da ferramenta, você também precisa dizer ao Nextflow que o nome do arquivo de saída esperado mudou.

=== "Depois"

    ```groovy title="genomics-2.nf" linenums="50" hl_lines="2 3"
        output:
        path "${input_bam}.g.vcf"     , emit: vcf
        path "${input_bam}.g.vcf.idx" , emit: idx
    ```

=== "Antes"

    ```groovy title="genomics-2.nf" linenums="50" hl_lines="2 3"
        output:
        path "${input_bam}.vcf"     , emit: vcf
        path "${input_bam}.vcf.idx" , emit: idx
    ```

### 1.4. Atualizar os alvos de publicação para as novas saídas GVCF

Já que agora estamos produzindo GVCFs em vez de VCFs, devemos atualizar a seção `publish:` do fluxo de trabalho para usar nomes mais descritivos.
Também vamos organizar os arquivos GVCF em seu próprio subdiretório para maior clareza.

=== "Depois"

    ```groovy title="genomics-2.nf" linenums="88" hl_lines="3 4"
        publish:
        indexed_bam = SAMTOOLS_INDEX.out
        gvcf = GATK_HAPLOTYPECALLER.out.vcf
        gvcf_idx = GATK_HAPLOTYPECALLER.out.idx
    ```

=== "Antes"

    ```groovy title="genomics-2.nf" linenums="88"
        publish:
        indexed_bam = SAMTOOLS_INDEX.out
        vcf = GATK_HAPLOTYPECALLER.out.vcf
        vcf_idx = GATK_HAPLOTYPECALLER.out.idx
    ```

### 1.5. Atualizar o bloco de saída para a nova estrutura de diretório

Também precisamos atualizar o bloco `output` para colocar os arquivos GVCF em um subdiretório `gvcf`.

=== "Depois"

    ```groovy title="genomics-2.nf" linenums="94" hl_lines="3 5 6 8 9"
    output {
        indexed_bam {
            path 'indexed_bam'
        }
        gvcf {
            path 'gvcf'
        }
        gvcf_idx {
            path 'gvcf'
        }
    }
    ```

=== "Antes"

    ```groovy title="genomics-2.nf" linenums="94"
    output {
        indexed_bam {
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

### 1.6. Executar o pipeline novamente

Vamos executá-lo com `-resume` desta vez.

```bash
nextflow run genomics-2.nf -resume
```

??? success "Saída do comando"

    ```console
    N E X T F L O W   ~  version 25.10.2

    ┃ Launching `genomics-2.nf` [nostalgic_franklin] DSL2 - revision: f2c0a93c6a

    executor >  local (3)
    [cc/fbc705] SAMTOOLS_INDEX (3)       | 3 of 3, cached: 3 ✔
    [27/0d7eb9] GATK_HAPLOTYPECALLER (2) | 3 of 3 ✔
    ```

Desta vez funciona.

A própria saída do Nextflow não parece diferente (comparada a uma execução bem-sucedida no modo VCF normal), mas agora podemos encontrar os arquivos `.g.vcf` e seus respectivos arquivos de índice, para todas as três amostras, organizados em subdiretórios.

??? abstract "Conteúdo do diretório (symlinks encurtados)"

    ```console
    results_genomics/
    ├── gvcf/
    │   ├── reads_father.bam.g.vcf -> */27/0d7eb9*/reads_father.bam.g.vcf
    │   ├── reads_father.bam.g.vcf.idx -> */27/0d7eb9*/reads_father.bam.g.vcf.idx
    │   ├── reads_mother.bam.g.vcf -> */e4/4ed55e*/reads_mother.bam.g.vcf
    │   ├── reads_mother.bam.g.vcf.idx -> */e4/4ed55e*/reads_mother.bam.g.vcf.idx
    │   ├── reads_son.bam.g.vcf -> */08/e95962*/reads_son.bam.g.vcf
    │   └── reads_son.bam.g.vcf.idx -> */08/e95962*/reads_son.bam.g.vcf.idx
    └── indexed_bam/
        ├── reads_father.bam -> */9a/c7a873*/reads_father.bam
        ├── reads_father.bam.bai -> */9a/c7a873*/reads_father.bam.bai
        ├── reads_mother.bam -> */f1/8d8486*/reads_mother.bam
        ├── reads_mother.bam.bai -> */f1/8d8486*/reads_mother.bam.bai
        ├── reads_son.bam -> */cc/fbc705*/reads_son.bam
        └── reads_son.bam.bai -> */cc/fbc705*/reads_son.bam.bai
    ```

Se você abrir um dos arquivos GVCF e percorrer por ele, poderá verificar que o GATK HaplotypeCaller produziu arquivos GVCF conforme solicitado.

### Conclusão

Ok, esta foi mínima em termos de aprendizado de Nextflow...
Mas foi uma boa oportunidade para reiterar a importância do bloco de saída do processo!

### O que vem a seguir?

Aprender a coletar o conteúdo de um canal e passá-los para o próximo processo como uma única entrada.

---

## 2. Coletar e combinar os dados GVCF em todas as amostras

Agora precisamos combinar os dados de todos os GVCFs por amostra em uma forma que suporte a análise de genotipagem conjunta que queremos fazer.

### 2.1. Definir o processo que vai combinar os GVCFs

Como lembrete do que fizemos anteriormente na seção de aquecimento, combinar os GVCFs é um trabalho para a ferramenta GATK GenomicsDBImport, que produzirá um armazenamento de dados no chamado formato GenomicsDB.

Vamos escrever um novo processo para definir como isso vai funcionar, baseado no comando que usamos anteriormente na seção de aquecimento.

```groovy title="genomics-2.nf" linenums="66"
/*
 * Combinar GVCFs em armazenamento de dados GenomicsDB
 */
process GATK_GENOMICSDB {

    container "community.wave.seqera.io/library/gatk4:4.5.0.0--730ee8817e436867"

    input:
    path all_gvcfs
    path all_idxs
    path interval_list
    val cohort_name

    output:
    path "${cohort_name}_gdb"

    script:
    """
    gatk GenomicsDBImport \
        -V ${all_gvcfs} \
        -L ${interval_list} \
        --genomicsdb-workspace-path ${cohort_name}_gdb
    """
}
```

O que você acha, parece razoável?

Vamos conectá-lo e ver o que acontece.

### 2.2. Adicionar um parâmetro `cohort_name` com um valor padrão

Precisamos fornecer um nome arbitrário para a coorte.
Mais tarde na série de treinamento você aprenderá como usar metadados de amostra para esse tipo de coisa, mas por enquanto apenas declaramos um parâmetro CLI usando `params` e damos a ele um valor padrão por conveniência.

```groovy title="genomics-2.nf" linenums="16"
    // Nome base para o arquivo de saída final
    cohort_name: String = "family_trio"
```

### 2.3. Reunir as saídas de GATK_HAPLOTYPECALLER entre amostras

Se fôssemos apenas conectar o canal de saída do processo `GATK_HAPLOTYPECALLER` como está, o Nextflow chamaria o processo em cada GVCF de amostra separadamente.
No entanto, queremos agrupar todos os três GVCFs (e seus arquivos de índice) de forma que o Nextflow entregue todos eles juntos a uma única chamada de processo.

Boas notícias: podemos fazer isso usando o operador de canal `collect()`. Vamos adicionar as seguintes linhas ao corpo do `workflow`, logo após a chamada a GATK_HAPLOTYPECALLER:

```groovy title="genomics-2.nf" linenums="118"
// Coletar saídas de chamada de variantes entre amostras
all_gvcfs_ch = GATK_HAPLOTYPECALLER.out.vcf.collect()
all_idxs_ch = GATK_HAPLOTYPECALLER.out.idx.collect()
```

Isso parece um pouco complicado? Vamos quebrar isso e traduzir para linguagem simples.

1. Estamos pegando o canal de saída do processo `GATK_HAPLOTYPECALLER`, referenciado usando a propriedade `.out`.
2. Cada 'elemento' saindo do canal é um par de arquivos: o GVCF e seu arquivo de índice, nessa ordem porque essa é a ordem em que estão listados no bloco de saída do processo. Convenientemente, porque na última sessão nomeamos as saídas deste processo (usando `emit:`), podemos selecionar os GVCFs por um lado adicionando `.vcf` e os arquivos de índice por outro adicionando `.idx` após a propriedade `.out`. Se não tivéssemos nomeado essas saídas, teríamos que nos referir a elas por `.out[0]` e `.out[1]`, respectivamente.
3. Anexamos o operador de canal `collect()` para agrupar todos os arquivos GVCF juntos em um único elemento em um novo canal chamado `all_gvcfs_ch`, e fazemos o mesmo com os arquivos de índice para formar o novo canal chamado `all_idxs_ch`.

!!! tip

    Se você está tendo dificuldade em visualizar exatamente o que está acontecendo aqui, lembre-se de que você pode usar o operador `view()` para inspecionar o conteúdo dos canais antes e depois de aplicar operadores de canal.

Os canais resultantes `all_gvcfs_ch` e `all_idxs_ch` são o que vamos conectar ao processo `GATK_GENOMICSDB` que acabamos de escrever.

!!! note

    Caso você estivesse se perguntando, coletamos os GVCFs e seus arquivos de índice separadamente porque o comando GATK GenomicsDBImport só quer ver os caminhos dos arquivos GVCF. Felizmente, já que o Nextflow vai organizar todos os arquivos juntos para execução, não precisamos nos preocupar com a ordem dos arquivos como fizemos para BAMs e seus índices na Parte 1.

### 2.4. Adicionar uma chamada ao bloco de fluxo de trabalho para executar GATK_GENOMICSDB

Temos um processo, e temos canais de entrada. Só precisamos adicionar a chamada do processo.

```groovy title="genomics-2.nf" linenums="122"
    // Combinar GVCFs em um armazenamento de dados GenomicsDB
    GATK_GENOMICSDB(
        all_gvcfs_ch,
        all_idxs_ch,
        intervals_file,
        params.cohort_name
    )
```

Ok, tudo está conectado.

### 2.5. Executar o fluxo de trabalho

Vamos ver se isso funciona.

```bash
nextflow run genomics-2.nf -resume
```

??? failure "Saída do comando"

    ```console
    N E X T F L O W   ~  version 25.10.2

    ┃ Launching `genomics-2.nf` [disturbed_bell] DSL2 - revision: 57942246cc

    executor >  local (1)
    [f1/8d8486] SAMTOOLS_INDEX (1)       | 3 of 3, cached: 3 ✔
    [e4/4ed55e] GATK_HAPLOTYPECALLER (2) | 3 of 3, cached: 3 ✔
    [51/d350ea] GATK_GENOMICSDB          | 0 of 1
    ERROR ~ Error executing process > 'GATK_GENOMICSDB'

    Caused by:
      Process `GATK_GENOMICSDB` terminated with an error exit status (1)

    Command executed:

      gatk GenomicsDBImport         -V reads_son.bam.g.vcf reads_father.bam.g.vcf reads_mother.bam.g.vcf         -L intervals.bed         --genomicsdb-workspace-path family_trio_gdb
    ```

Ele executa razoavelmente rápido, já que estamos executando com `-resume`, mas falha!

Ah. Pelo lado positivo, vemos que o Nextflow pegou o processo `GATK_GENOMICSDB`, e especificamente o chamou apenas uma vez.
Isso sugere que a abordagem `collect()` funcionou, até certo ponto.
Mas, e é um grande mas, a chamada do processo falhou.

Quando investigamos a saída do console acima, podemos ver que o comando executado não está correto.

Você consegue identificar o erro?
Olhe para este pedaço: `-V reads_father.bam.g.vcf reads_son.bam.g.vcf reads_mother.bam.g.vcf`

Demos ao `gatk GenomicsDBImport` múltiplos arquivos GVCF para um único argumento `-V`, mas a ferramenta espera um argumento `-V` separado para cada arquivo GVCF.

Como lembrete, este foi o comando que executamos no contêiner:

```bash
gatk GenomicsDBImport \
    -V reads_mother.g.vcf \
    -V reads_father.g.vcf \
    -V reads_son.g.vcf \
    -L /data/ref/intervals.bed \
    --genomicsdb-workspace-path family_trio_gdb
```

Então isso significa que precisamos de alguma forma transformar nosso pacote de arquivos GVCF em uma string de comando formatada adequadamente.

### 2.6. Construir uma linha de comando com um argumento `-V` separado para cada GVCF de entrada

É aqui que o Nextflow ser baseado em Groovy é útil, porque vai nos permitir usar algumas manipulações de string bastante diretas para construir a string de comando necessária.

Especificamente, usando esta sintaxe: `all_gvcfs.collect { gvcf -> "-V ${gvcf}" }.join(' ')`

Mais uma vez, vamos quebrar isso em seus componentes.

1. Primeiro, pegamos o conteúdo do canal de entrada `all_gvcfs` e aplicamos `.collect()` nele (assim como antes).
2. Isso nos permite passar cada caminho de arquivo GVCF individual no pacote para o **closure**, `{ gvcf -> "-V ${gvcf}" }`, onde `gvcf` se refere a esse caminho de arquivo GVCF.
   O closure é uma mini-função que usamos para prefixar `-V ` ao caminho do arquivo, na forma de `"-V ${gvcf}"`.
3. Então usamos `.join(' ')` para concatenar todas as três strings com um único espaço como separador.

Com um exemplo concreto, fica assim:

1. Temos três arquivos:

   `[A.ext, B.ext, C.ext]`

2. O closure modifica cada um para criar as strings:

   `"-V A.ext", "-V B.ext", "-V C.ext"`

3. A operação `.join(' ')` gera a string final:

   `"-V A.ext -V B.ext -V C.ext"`

Uma vez que temos essa string, podemos atribuí-la a uma variável local, `gvcfs_line`, definida com a palavra-chave `def`:

`def gvcfs_line = all_gvcfs.collect { gvcf -> "-V ${gvcf}" }.join(' ')`

Ok, então temos nossa coisinha de manipulação de string. Onde a colocamos?

Queremos que isso vá dentro da definição do processo em algum lugar, porque queremos fazer isso _depois_ que canalizamos os caminhos dos arquivos GVCF para o processo.
Isso é porque o Nextflow deve vê-los como caminhos de arquivo para organizar os próprios arquivos corretamente para execução.

Mas _onde_ no processo podemos adicionar isso?

Fato curioso: você pode adicionar código arbitrário depois de `script:` e antes do `"""` !

Ótimo, vamos adicionar nossa linha de manipulação de string lá então, e atualizar o comando `gatk GenomicsDBImport` para usar a string concatenada que ela produz.

=== "Depois"

    ```groovy title="genomics-2.nf" linenums="87"  hl_lines="2 5"
        script:
        def gvcfs_line = all_gvcfs.collect { gvcf -> "-V ${gvcf}" }.join(' ')
        """
        gatk GenomicsDBImport \
            ${gvcfs_line} \
            -L ${interval_list} \
            --genomicsdb-workspace-path ${cohort_name}_gdb
        """
    ```

=== "Antes"

    ```groovy title="genomics-2.nf" linenums="87"  hl_lines="4"
        script:
        """
        gatk GenomicsDBImport \
            -V ${all_gvcfs} \
            -L ${interval_list} \
            --genomicsdb-workspace-path ${cohort_name}_gdb
        """
    ```

Isso deve ser tudo que é necessário para fornecer as entradas ao `gatk GenomicsDBImport` corretamente.

!!! tip

    Quando você atualizar o comando `gatk GenomicsDBImport`, certifique-se de remover o prefixo `-V ` quando você trocar pela variável `${gvcfs_line}`.

### 2.7. Executar o fluxo de trabalho para verificar que ele gera a saída GenomicsDB conforme esperado

Certo, vamos ver se isso resolveu o problema.

```bash
nextflow run genomics-2.nf -resume
```

??? success "Saída do comando"

    ```console
    N E X T F L O W   ~  version 25.10.2

    ┃ Launching `genomics-2.nf` [peaceful_gates] DSL2 - revision: ca0bf847ed

    executor >  local (1)
    [cc/fbc705] SAMTOOLS_INDEX (3)       | 3 of 3, cached: 3 ✔
    [27/0d7eb9] GATK_HAPLOTYPECALLER (2) | 3 of 3, cached: 3 ✔
    [76/d13861] GATK_GENOMICSDB          | 1 of 1 ✔
    ```

Aha! Parece estar funcionando agora.

As duas primeiras etapas foram puladas com sucesso, e a terceira etapa funcionou perfeitamente desta vez.
O armazenamento de dados GenomicsDB é criado no diretório de trabalho mas não publicado em resultados, já que é apenas um formato intermediário que usaremos para genotipagem conjunta.

A propósito, não tivemos que fazer nada especial para lidar com a saída sendo um diretório em vez de um único arquivo.

### Conclusão

Agora você sabe como coletar saídas de um canal e agrupá-las como uma única entrada para outro processo.
Você também sabe como construir uma linha de comando para fornecer entradas a uma determinada ferramenta com a sintaxe apropriada.

### O que vem a seguir?

Aprender como adicionar um segundo comando ao mesmo processo.

---

## 3. Executar a etapa de genotipagem conjunta como parte do mesmo processo

Agora que temos as chamadas de variantes genômicas combinadas, podemos executar a ferramenta de genotipagem conjunta, que produzirá a saída final que realmente nos interessa: o VCF de chamadas de variantes em nível de coorte.

Por razões logísticas, decidimos incluir a genotipagem conjunta dentro do mesmo processo.

### 3.1. Renomear o processo de GATK_GENOMICSDB para GATK_JOINTGENOTYPING

Como o processo estará executando mais de uma ferramenta, mudamos seu nome para referir-se à operação geral em vez de um único nome de ferramenta.

=== "Depois"

    ```groovy title="genomics-2.nf"
    /*
     * Combinar GVCFs em armazenamento de dados GenomicsDB e executar genotipagem conjunta para produzir chamadas em nível de coorte
     */
    process GATK_JOINTGENOTYPING {
    ```

=== "Antes"

    ```groovy title="genomics-2.nf"
    /*
     * Combinar GVCFs em armazenamento de dados GenomicsDB
     */
    process GATK_GENOMICSDB {
    ```

Lembre-se de manter seus nomes de processo o mais descritivos possível, para maximizar a legibilidade para seus colegas —e seu eu futuro!

### 3.2. Adicionar o comando de genotipagem conjunta ao processo GATK_JOINTGENOTYPING

Simplesmente adicione o segundo comando após o primeiro dentro da seção de script.

=== "Depois"

    ```groovy title="genomics-2.nf" linenums="89"  hl_lines="6-10"
        """
        gatk GenomicsDBImport \
            ${gvcfs_line} \
            -L ${interval_list} \
            --genomicsdb-workspace-path ${cohort_name}_gdb

        gatk GenotypeGVCFs \
            -R ${ref_fasta} \
            -V gendb://${cohort_name}_gdb \
            -L ${interval_list} \
            -O ${cohort_name}.joint.vcf
        """
    ```

=== "Antes"

    ```groovy title="genomics-2.nf" linenums="89"
        """
        gatk GenomicsDBImport \
            ${gvcfs_line} \
            -L ${interval_list} \
            --genomicsdb-workspace-path ${cohort_name}_gdb
        """
    ```

Os dois comandos serão executados em série, da mesma forma que seriam se os executássemos manualmente no terminal.

### 3.3. Adicionar os arquivos do genoma de referência às definições de entrada do processo GATK_JOINTGENOTYPING

O segundo comando requer os arquivos do genoma de referência, então precisamos adicioná-los às entradas do processo.

=== "Depois"

    ```groovy title="genomics-2.nf" linenums="78"  hl_lines="6-8"
        input:
        path all_gvcfs
        path all_idxs
        path interval_list
        val cohort_name
        path ref_fasta
        path ref_index
        path ref_dict
    ```

=== "Antes"

    ```groovy title="genomics-2.nf" linenums="78"
        input:
        path all_gvcfs
        path all_idxs
        path interval_list
        val cohort_name
    ```

Pode parecer irritante digitar isso tudo, mas lembre-se, você só digita uma vez, e então pode executar o fluxo de trabalho um milhão de vezes. Vale a pena?

### 3.4. Atualizar a definição de saída do processo para emitir o VCF de chamadas de variantes em nível de coorte

Realmente não nos importamos em salvar o armazenamento de dados GenomicsDB, que é apenas um formato intermediário que só existe por razões logísticas, então podemos simplesmente removê-lo do bloco de saída se quisermos.

A saída que realmente nos interessa é o VCF produzido pelo comando de genotipagem conjunta.

=== "Depois"

    ```groovy title="genomics-2.nf" linenums="87" hl_lines="2 3"
        output:
        path "${cohort_name}.joint.vcf"     , emit: vcf
        path "${cohort_name}.joint.vcf.idx" , emit: idx
    ```

=== "Antes"

    ```groovy title="genomics-2.nf" linenums="87"
        output:
        path "${cohort_name}_gdb"
    ```

Estamos quase terminando!

### 3.5. Atualizar a chamada do processo de GATK_GENOMICSDB para GATK_JOINTGENOTYPING

Não vamos esquecer de renomear a chamada do processo no corpo do fluxo de trabalho de GATK_GENOMICSDB para GATK_JOINTGENOTYPING. E já que estamos nisso, também devemos adicionar os arquivos do genoma de referência como entradas, já que precisamos fornecê-los à ferramenta de genotipagem conjunta.

=== "Depois"

    ```groovy title="genomics-2.nf" linenums="126"
    // Combinar GVCFs em um armazenamento de dados GenomicsDB e aplicar genotipagem conjunta
    GATK_JOINTGENOTYPING(
        all_gvcfs_ch,
        all_idxs_ch,
        intervals_file,
        params.cohort_name,
        ref_file,
        ref_index_file,
        ref_dict_file
    )
    ```

=== "Antes"

    ```groovy title="genomics-2.nf" linenums="126"
    // Combinar GVCFs em um armazenamento de dados GenomicsDB
    GATK_GENOMICSDB(
        all_gvcfs_ch,
        all_idxs_ch,
        intervals_file,
        params.cohort_name
    )
    ```

Agora o processo está completamente conectado.

### 3.6. Adicionar o VCF conjunto à seção de publicação

Precisamos publicar as saídas do VCF conjunto do novo processo.
Adicione estas linhas à seção `publish:` do fluxo de trabalho:

```groovy title="genomics-2.nf" linenums="145"
    joint_vcf = GATK_JOINTGENOTYPING.out.vcf
    joint_vcf_idx = GATK_JOINTGENOTYPING.out.idx
```

### 3.7. Adicionar os alvos do VCF conjunto ao bloco de saída

Finalmente, adicione alvos de saída para os arquivos VCF conjunto.
Vamos colocá-los na raiz do diretório de resultados já que esta é a saída final.

```groovy title="genomics-2.nf" linenums="157"
    joint_vcf {
        path '.'
    }
    joint_vcf_idx {
        path '.'
    }
```

Agora tudo deve estar completamente conectado.

### 3.8. Executar o fluxo de trabalho

Finalmente, podemos executar o fluxo de trabalho modificado...

```bash
nextflow run genomics-2.nf -resume
```

??? success "Saída do comando"

    ```console
    N E X T F L O W   ~  version 25.10.2

    ┃ Launching `genomics-2.nf` [crazy_marconi] DSL2 - revision: 5da9afc841

    executor >  local (1)
    [9a/c7a873] SAMTOOLS_INDEX (2)       | 3 of 3, cached: 3 ✔
    [e4/4ed55e] GATK_HAPLOTYPECALLER (2) | 3 of 3, cached: 3 ✔
    [a6/7cc8ed] GATK_JOINTGENOTYPING     | 1 of 1 ✔
    ```

E funciona!

Você encontrará o arquivo de saída final, `family_trio.joint.vcf` (e seu índice de arquivo), no diretório de resultados.

??? abstract "Conteúdo do diretório (symlinks encurtados)"

    ```console
    results_genomics/
    ├── family_trio.joint.vcf -> */a6/7cc8ed*/family_trio.joint.vcf
    ├── family_trio.joint.vcf.idx -> */a6/7cc8ed*/family_trio.joint.vcf.idx
    ├── gvcf/
    │   ├── reads_father.bam.g.vcf -> */27/0d7eb9*/reads_father.bam.g.vcf
    │   ├── reads_father.bam.g.vcf.idx -> */27/0d7eb9*/reads_father.bam.g.vcf.idx
    │   ├── reads_mother.bam.g.vcf -> */e4/4ed55e*/reads_mother.bam.g.vcf
    │   ├── reads_mother.bam.g.vcf.idx -> */e4/4ed55e*/reads_mother.bam.g.vcf.idx
    │   ├── reads_son.bam.g.vcf -> */08/e95962*/reads_son.bam.g.vcf
    │   └── reads_son.bam.g.vcf.idx -> */08/e95962*/reads_son.bam.g.vcf.idx
    └── indexed_bam/
        ├── reads_father.bam -> */9a/c7a873*/reads_father.bam
        ├── reads_father.bam.bai -> */9a/c7a873*/reads_father.bam.bai
        ├── reads_mother.bam -> */f1/8d8486*/reads_mother.bam
        ├── reads_mother.bam.bai -> */f1/8d8486*/reads_mother.bam.bai
        ├── reads_son.bam -> */cc/fbc705*/reads_son.bam
        └── reads_son.bam.bai -> */cc/fbc705*/reads_son.bam.bai
    ```

Se você é do tipo cético, pode clicar no arquivo VCF conjunto para abri-lo e verificar que o fluxo de trabalho gerou as mesmas chamadas de variantes que você obteve executando as ferramentas manualmente no início desta seção.

```console title="family_trio.joint.vcf" linenums="40"
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	reads_father	reads_mother	reads_son
20_10037292_10066351	3480	.	C	CT	1625.89	.	AC=5;AF=0.833;AN=6;BaseQRankSum=0.220;DP=85;ExcessHet=0.0000;FS=2.476;MLEAC=5;MLEAF=0.833;MQ=60.00;MQRankSum=0.00;QD=21.68;ReadPosRankSum=-1.147e+00;SOR=0.487	GT:AD:DP:GQ:PL	0/1:15,16:31:99:367,0,375	1/1:0,18:18:54:517,54,0	1/1:0,26:26:78:756,78,0
20_10037292_10066351	3520	.	AT	A	1678.89	.	AC=5;AF=0.833;AN=6;BaseQRankSum=1.03;DP=80;ExcessHet=0.0000;FS=2.290;MLEAC=5;MLEAF=0.833;MQ=60.00;MQRankSum=0.00;QD=22.39;ReadPosRankSum=0.701;SOR=0.730	GT:AD:DP:GQ:PL	0/1:18,13:31:99:296,0,424	1/1:0,18:18:54:623,54,0	1/1:0,26:26:78:774,78,0
20_10037292_10066351	3529	.	T	A	154.29	.	AC=1;AF=0.167;AN=6;BaseQRankSum=-5.440e-01;DP=104;ExcessHet=0.0000;FS=1.871;MLEAC=1;MLEAF=0.167;MQ=60.00;MQRankSum=0.00;QD=7.71;ReadPosRankSum=-1.158e+00;SOR=1.034	GT:AD:DP:GQ:PL	0/0:44,0:44:99:0,112,1347	0/1:12,8:20:99:163,0,328	0/0:39,0:39:99:0,105,1194
```

Agora você tem um fluxo de trabalho de chamada conjunta de variantes automatizado e totalmente reproduzível!

!!! note

    Tenha em mente que os arquivos de dados que fornecemos a você cobrem apenas uma pequena porção do cromossomo 20.
    O tamanho real de um conjunto de chamadas de variantes seria contado em milhões de variantes.
    É por isso que usamos apenas pequenos subconjuntos de dados para fins de treinamento!

### Conclusão

Você sabe como usar alguns operadores comuns, bem como closures Groovy, para controlar o fluxo de dados em seu fluxo de trabalho.

### O que vem a seguir?

Celebre seu sucesso e faça uma pausa bem merecida.

Na próxima parte deste curso, você aprenderá como modularizar seu fluxo de trabalho extraindo definições de processo em módulos reutilizáveis.
