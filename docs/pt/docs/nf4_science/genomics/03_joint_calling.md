# Parte 3: Chamada conjunta de variantes em uma coorte

Na Parte 2, você construiu um pipeline de chamada de variantes por amostra que processou os dados de cada amostra independentemente.
Agora vamos estendê-lo para implementar a chamada conjunta de variantes, conforme abordado na [Parte 1](01_method.md).

## Tarefa

Nesta parte do curso, vamos estender o fluxo de trabalho para fazer o seguinte:

<figure class="excalidraw">
--8<-- "docs/en/docs/nf4_science/genomics/img/hello-gatk-2.svg"
</figure>

1. Gerar um arquivo de índice para cada arquivo BAM de entrada usando o Samtools
2. Executar o GATK HaplotypeCaller em cada arquivo BAM de entrada para gerar um GVCF de chamadas de variantes genômicas por amostra
3. Coletar todos os GVCFs e combiná-los em um armazenamento de dados GenomicsDB
4. Executar a genotipagem conjunta no armazenamento de dados GVCF combinado para produzir um VCF em nível de coorte

Esta parte se baseia diretamente no fluxo de trabalho produzido pela Parte 2.

??? info "Como começar a partir desta seção"

    Esta seção do curso pressupõe que você completou a [Parte 2: Chamada de variantes por amostra](./02_per_sample_variant_calling.md) e tem um pipeline `genomics.nf` funcional.

    Se você não completou a Parte 2 ou deseja começar do zero para esta parte, pode usar a solução da Parte 2 como ponto de partida.
    Execute estes comandos dentro do diretório `nf4-science/genomics/`:

    ```bash
    cp solutions/part2/genomics-2.nf genomics.nf
    cp solutions/part2/nextflow.config .
    cp solutions/part2/modules/* modules/
    ```

    Isso lhe dá um fluxo de trabalho completo de chamada de variantes por amostra.
    Você pode testar se ele é executado com sucesso executando o seguinte comando:

    ```bash
    nextflow run genomics.nf -profile test
    ```

## Plano de aula

Dividimos isso em duas etapas:

1. **Modificar a etapa de chamada de variantes por amostra para produzir um GVCF.**
   Isso abrange a atualização de comandos e saídas do processo.
2. **Adicionar uma etapa de genotipagem conjunta que combina e genotipa os GVCFs por amostra.**
   Isso introduz o operador `collect()`, closures Groovy para construção de linha de comando e processos com múltiplos comandos.

!!! note

     Certifique-se de estar no diretório de trabalho correto:
     `cd /workspaces/training/nf4-science/genomics`

---

## 1. Modificar a etapa de chamada de variantes por amostra para produzir um GVCF

O pipeline da Parte 2 produz arquivos VCF, mas a chamada conjunta requer arquivos GVCF.
Precisamos ativar o modo de chamada de variantes GVCF e atualizar a extensão do arquivo de saída.

Relembre o comando de chamada de variantes GVCF da [Parte 1](01_method.md):

```bash
gatk HaplotypeCaller \
        -R /data/ref/ref.fasta \
        -I /data/bam/reads_mother.bam \
        -O reads_mother.g.vcf \
        -L /data/ref/intervals.bed \
        -ERC GVCF
```

Comparado ao comando base do HaplotypeCaller que encapsulamos na Parte 2, as diferenças são o parâmetro `-ERC GVCF` e a extensão de saída `.g.vcf`.

### 1.1. Informar ao HaplotypeCaller para emitir um GVCF e atualizar a extensão de saída

Abra o arquivo de módulo `modules/gatk_haplotypecaller.nf` para fazer duas alterações:

- Adicionar o parâmetro `-ERC GVCF` ao comando GATK HaplotypeCaller;
- Atualizar o caminho do arquivo de saída para usar a extensão `.g.vcf` correspondente, conforme convenção do GATK.

Certifique-se de adicionar uma barra invertida (`\`) no final da linha anterior quando adicionar `-ERC GVCF`.

=== "Depois"

    ```groovy title="modules/gatk_haplotypecaller.nf" linenums="22" hl_lines="5 7"
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

    ```groovy title="modules/gatk_haplotypecaller.nf" linenums="22" hl_lines="5"
        """
        gatk HaplotypeCaller \
            -R ${ref_fasta} \
            -I ${input_bam} \
            -O ${input_bam}.vcf \
            -L ${interval_list}
        """
    ```

Também precisamos atualizar o bloco de saída para corresponder à nova extensão de arquivo.
Como mudamos a saída do comando de `.vcf` para `.g.vcf`, o bloco `output:` do processo deve refletir a mesma mudança.

### 1.2. Atualizar a extensão do arquivo de saída no bloco de saídas do processo

=== "Depois"

    ```groovy title="modules/gatk_haplotypecaller.nf" linenums="17" hl_lines="2 3"
        output:
        path "${input_bam}.g.vcf"     , emit: vcf
        path "${input_bam}.g.vcf.idx" , emit: idx
    ```

=== "Antes"

    ```groovy title="modules/gatk_haplotypecaller.nf" linenums="17" hl_lines="2 3"
        output:
        path "${input_bam}.vcf"     , emit: vcf
        path "${input_bam}.vcf.idx" , emit: idx
    ```

Também precisamos atualizar a configuração de publicação e saída do fluxo de trabalho para refletir as novas saídas GVCF.

### 1.3. Atualizar os alvos de publicação para as novas saídas GVCF

Como agora estamos produzindo GVCFs em vez de VCFs, devemos atualizar a seção `publish:` do fluxo de trabalho para usar nomes mais descritivos.
Também organizaremos os arquivos GVCF em seu próprio subdiretório para maior clareza.

=== "Depois"

    ```groovy title="genomics.nf" linenums="47" hl_lines="3 4"
        publish:
        indexed_bam = SAMTOOLS_INDEX.out
        gvcf = GATK_HAPLOTYPECALLER.out.vcf
        gvcf_idx = GATK_HAPLOTYPECALLER.out.idx
    ```

=== "Antes"

    ```groovy title="genomics.nf" linenums="47"
        publish:
        indexed_bam = SAMTOOLS_INDEX.out
        vcf = GATK_HAPLOTYPECALLER.out.vcf
        vcf_idx = GATK_HAPLOTYPECALLER.out.idx
    ```

Agora atualize o bloco de saída para corresponder.

### 1.4. Atualizar o bloco de saída para a nova estrutura de diretórios

Também precisamos atualizar o bloco `output` para colocar os arquivos GVCF em um subdiretório `gvcf`.

=== "Depois"

    ```groovy title="genomics.nf" linenums="53" hl_lines="3 5 6 8 9"
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

    ```groovy title="genomics.nf" linenums="53"
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

Com o módulo, os alvos de publicação e o bloco de saída todos atualizados, podemos testar as alterações.

### 1.5. Executar o pipeline

Execute o fluxo de trabalho para verificar se as alterações funcionam.

```bash
nextflow run genomics.nf
```

??? success "Saída do comando"

    ```console
    N E X T F L O W   ~  version 25.10.2

    ┃ Launching `genomics.nf` [nostalgic_franklin] DSL2 - revision: f2c0a93c6a

    executor >  local (6)
    [cc/fbc705] SAMTOOLS_INDEX (3)       | 3 of 3 ✔
    [27/0d7eb9] GATK_HAPLOTYPECALLER (2) | 3 of 3 ✔
    ```

A saída do Nextflow parece a mesma de antes, mas os arquivos `.g.vcf` e seus arquivos de índice agora estão organizados em subdiretórios.

??? abstract "Conteúdo do diretório (links simbólicos encurtados)"

    ```console
    results/
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

Se você abrir um dos arquivos GVCF e percorrê-lo, poderá verificar que o GATK HaplotypeCaller produziu arquivos GVCF conforme solicitado.

### Conclusão

Quando você altera o nome do arquivo de saída de um comando de ferramenta, o bloco `output:` do processo e a configuração de publicação/saída devem ser atualizados para corresponder.

### O que vem a seguir?

Aprenda a coletar o conteúdo de um canal e passá-lo para o próximo processo como uma única entrada.

---

## 2. Adicionar uma etapa de genotipagem conjunta

Agora precisamos coletar os GVCFs por amostra, combiná-los em um armazenamento de dados GenomicsDB e executar a genotipagem conjunta para produzir um VCF em nível de coorte.
Conforme abordado na [Parte 1](01_method.md), esta é uma operação de duas ferramentas: GenomicsDBImport combina os GVCFs, depois GenotypeGVCFs produz as chamadas de variantes finais.
Vamos encapsular ambas as ferramentas em um único processo chamado `GATK_JOINTGENOTYPING`.

Relembre os dois comandos da [Parte 1](01_method.md):

```bash
gatk GenomicsDBImport \
    -V reads_mother.g.vcf \
    -V reads_father.g.vcf \
    -V reads_son.g.vcf \
    -L /data/ref/intervals.bed \
    --genomicsdb-workspace-path family_trio_gdb
```

```bash
gatk GenotypeGVCFs \
    -R /data/ref/ref.fasta \
    -V gendb://family_trio_gdb \
    -O family_trio.vcf
```

O primeiro comando recebe os GVCFs por amostra e um arquivo de intervalos, e produz um armazenamento de dados GenomicsDB.
O segundo recebe esse armazenamento de dados, um genoma de referência e produz o VCF final em nível de coorte.
O URI do contêiner é o mesmo do HaplotypeCaller: `community.wave.seqera.io/library/gatk4:4.5.0.0--730ee8817e436867`.

### 2.1. Configurar as entradas

O processo de genotipagem conjunta precisa de dois tipos de entradas que ainda não temos: um nome de coorte arbitrário e as saídas GVCF coletadas de todas as amostras agrupadas juntas.

#### 2.1.1. Adicionar um parâmetro `cohort_name`

Precisamos fornecer um nome arbitrário para a coorte.
Mais adiante na série de treinamento, você aprenderá como usar metadados de amostras para esse tipo de coisa, mas por enquanto apenas declaramos um parâmetro CLI usando `params` e damos a ele um valor padrão por conveniência.

=== "Depois"

    ```groovy title="genomics.nf" linenums="14" hl_lines="3-4"
        intervals: Path = "${projectDir}/data/ref/intervals.bed"

        // Base name for final output file
        cohort_name: String = "family_trio"
    }
    ```

=== "Antes"

    ```groovy title="genomics.nf" linenums="14"
        intervals: Path = "${projectDir}/data/ref/intervals.bed"
    }
    ```

#### 2.1.2. Coletar as saídas do HaplotypeCaller entre amostras

Se conectássemos o canal de saída de `GATK_HAPLOTYPECALLER` diretamente ao novo processo, o Nextflow chamaria o processo em cada GVCF de amostra separadamente.
Queremos agrupar todos os três GVCFs (e seus arquivos de índice) para que o Nextflow entregue todos eles juntos em uma única chamada de processo.

Podemos fazer isso usando o operador de canal `collect()`.
Adicione as seguintes linhas ao corpo do `workflow`, logo após a chamada a GATK_HAPLOTYPECALLER:

=== "Depois"

    ```groovy title="genomics.nf" hl_lines="4-6"
            intervals_file
        )

        // Collect variant calling outputs across samples
        all_gvcfs_ch = GATK_HAPLOTYPECALLER.out.vcf.collect()
        all_idxs_ch = GATK_HAPLOTYPECALLER.out.idx.collect()
    ```

=== "Antes"

    ```groovy title="genomics.nf"
            intervals_file
        )
    ```

Analisando isso:

1. Pegamos o canal de saída de `GATK_HAPLOTYPECALLER` usando a propriedade `.out`.
2. Como nomeamos as saídas usando `emit:` na seção 1, podemos selecionar os GVCFs com `.vcf` e os arquivos de índice com `.idx`. Sem saídas nomeadas, teríamos que usar `.out[0]` e `.out[1]`.
3. O operador `collect()` agrupa todos os arquivos em um único elemento, então `all_gvcfs_ch` contém todos os três GVCFs juntos, e `all_idxs_ch` contém todos os três arquivos de índice juntos.

Podemos coletar os GVCFs e seus arquivos de índice separadamente (em vez de mantê-los juntos em tuplas) porque o Nextflow colocará todos os arquivos de entrada juntos para execução, então os arquivos de índice estarão presentes junto aos GVCFs.

!!! tip

    Você pode usar o operador `view()` para inspecionar o conteúdo dos canais antes e depois de aplicar operadores de canal.

### 2.2. Escrever o processo de genotipagem conjunta e chamá-lo no fluxo de trabalho

Seguindo o mesmo padrão que usamos na Parte 2, escreveremos a definição do processo em um arquivo de módulo, importaremos para o fluxo de trabalho e o chamaremos nas entradas que acabamos de preparar.

#### 2.2.1. Construir uma string para dar a cada GVCF um argumento `-V`

Antes de começarmos a preencher a definição do processo, há uma coisa a resolver.
O comando GenomicsDBImport espera um argumento `-V` separado para cada arquivo GVCF, assim:

```bash
gatk GenomicsDBImport \
    -V reads_mother.bam.g.vcf \
    -V reads_father.bam.g.vcf \
    -V reads_son.bam.g.vcf \
    ...
```

Se escrevêssemos `-V ${all_gvcfs_ch}`, o Nextflow simplesmente concatenaria os nomes de arquivo e essa parte do comando ficaria assim:

```groovy
-V reads_mother.bam.g.vcf reads_father.bam.g.vcf reads_son.bam.g.vcf
```

Mas precisamos que a string fique assim:

```groovy
-V reads_mother.bam.g.vcf -V reads_father.bam.g.vcf -V reads_son.bam.g.vcf
```

Importante, precisamos construir essa string dinamicamente a partir de quaisquer arquivos que estejam no canal coletado.
O Nextflow (via Groovy) fornece uma maneira concisa de fazer isso:

```groovy
def gvcfs_line = all_gvcfs.collect { gvcf -> "-V ${gvcf}" }.join(' ')
```

Analisando isso:

1. `all_gvcfs.collect { gvcf -> "-V ${gvcf}" }` itera sobre cada caminho de arquivo e anexa `-V ` a ele, produzindo `["-V A.g.vcf", "-V B.g.vcf", "-V C.g.vcf"]`.
2. `.join(' ')` os concatena com espaços: `"-V A.g.vcf -V B.g.vcf -V C.g.vcf"`.
3. O resultado é atribuído a uma variável local `gvcfs_line` (definida com `def`), que podemos interpolar no template de comando.

Esta linha vai dentro do bloco `script:` do processo, antes do template de comando.
Você pode colocar código Groovy arbitrário entre `script:` e as `"""` de abertura do template de comando.

Então você poderá se referir a toda essa string como `gvcfs_line` no bloco `script:` do processo.

#### 2.2.2. Preencher o módulo para o processo de genotipagem conjunta

Agora podemos começar a escrever o processo completo.

Abra `modules/gatk_jointgenotyping.nf` e examine o esboço da definição do processo.

Vá em frente e preencha a definição do processo usando as informações fornecidas acima, depois verifique seu trabalho contra a solução na aba "Depois" abaixo.

=== "Antes"

    ```groovy title="modules/gatk_jointgenotyping.nf" linenums="1"
    #!/usr/bin/env nextflow

    /*
     * Combinar GVCFs em armazenamento de dados GenomicsDB e executar genotipagem conjunta para produzir chamadas em nível de coorte
     */
    process GATK_JOINTGENOTYPING {

        container

        input:

        output:

        script:
        """

        """
    }
    ```

=== "Depois"

    ```groovy title="modules/gatk_jointgenotyping.nf" linenums="1" hl_lines="8 11-17 20-21 24-25 29-33"
    #!/usr/bin/env nextflow

    /*
     * Combinar GVCFs em armazenamento de dados GenomicsDB e executar genotipagem conjunta para produzir chamadas em nível de coorte
     */
    process GATK_JOINTGENOTYPING {

        container "community.wave.seqera.io/library/gatk4:4.5.0.0--730ee8817e436867"

        input:
        path all_gvcfs
        path all_idxs
        path interval_list
        val cohort_name
        path ref_fasta
        path ref_index
        path ref_dict

        output:
        path "${cohort_name}.joint.vcf"     , emit: vcf
        path "${cohort_name}.joint.vcf.idx" , emit: idx

        script:
        def gvcfs_line = all_gvcfs.collect { gvcf -> "-V ${gvcf}" }.join(' ')
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
    }
    ```

Há várias coisas que vale a pena destacar aqui.

Como anteriormente, várias entradas são listadas mesmo que os comandos não as referenciem diretamente: `all_idxs`, `ref_index` e `ref_dict`.
Listá-las garante que o Nextflow coloque esses arquivos no diretório de trabalho junto aos arquivos que aparecem nos comandos, que o GATK espera encontrar com base em convenções de nomenclatura.

A variável `gvcfs_line` usa a closure Groovy descrita acima para construir os argumentos `-V` para GenomicsDBImport.

Este processo executa dois comandos em série, assim como você faria no terminal.
GenomicsDBImport combina os GVCFs por amostra em um armazenamento de dados, depois GenotypeGVCFs lê esse armazenamento de dados e produz o VCF final em nível de coorte.
O armazenamento de dados GenomicsDB (`${cohort_name}_gdb`) é um artefato intermediário usado apenas dentro do processo; ele não aparece no bloco de saída.

Uma vez que você tenha completado isso, o processo está pronto para uso.
Para usá-lo no fluxo de trabalho, você precisará importar o módulo e adicionar uma chamada de processo.

#### 2.2.3. Importar o módulo

Adicione a instrução de importação a `genomics.nf`, abaixo das instruções de importação existentes:

=== "Depois"

    ```groovy title="genomics.nf" linenums="21" hl_lines="3"
    include { SAMTOOLS_INDEX } from './modules/samtools_index.nf'
    include { GATK_HAPLOTYPECALLER } from './modules/gatk_haplotypecaller.nf'
    include { GATK_JOINTGENOTYPING } from './modules/gatk_jointgenotyping.nf'
    ```

=== "Antes"

    ```groovy title="genomics.nf" linenums="21"
    include { SAMTOOLS_INDEX } from './modules/samtools_index.nf'
    include { GATK_HAPLOTYPECALLER } from './modules/gatk_haplotypecaller.nf'
    ```

O processo agora está disponível no escopo do fluxo de trabalho.

#### 2.2.4. Adicionar a chamada do processo

Adicione a chamada a `GATK_JOINTGENOTYPING` no corpo do fluxo de trabalho, após as linhas `collect()`:

=== "Depois"

    ```groovy title="genomics.nf" hl_lines="3-12"
        all_idxs_ch = GATK_HAPLOTYPECALLER.out.idx.collect()

        // Combine GVCFs into a GenomicsDB data store and apply joint genotyping
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

    ```groovy title="genomics.nf"
        all_idxs_ch = GATK_HAPLOTYPECALLER.out.idx.collect()
    ```

O processo agora está totalmente conectado.
Em seguida, configuramos como as saídas são publicadas.

### 2.3. Configurar o tratamento de saída

Precisamos publicar as saídas VCF conjuntas.
Adicione alvos de publicação e entradas de bloco de saída para os resultados de genotipagem conjunta.

#### 2.3.1. Adicionar alvos de publicação para o VCF conjunto

Adicione o VCF conjunto e seu índice à seção `publish:` do fluxo de trabalho:

=== "Depois"

    ```groovy title="genomics.nf" hl_lines="5-6"
        publish:
        indexed_bam = SAMTOOLS_INDEX.out
        gvcf = GATK_HAPLOTYPECALLER.out.vcf
        gvcf_idx = GATK_HAPLOTYPECALLER.out.idx
        joint_vcf = GATK_JOINTGENOTYPING.out.vcf
        joint_vcf_idx = GATK_JOINTGENOTYPING.out.idx
    ```

=== "Antes"

    ```groovy title="genomics.nf"
        publish:
        indexed_bam = SAMTOOLS_INDEX.out
        gvcf = GATK_HAPLOTYPECALLER.out.vcf
        gvcf_idx = GATK_HAPLOTYPECALLER.out.idx
    ```

Agora atualize o bloco de saída para corresponder.

#### 2.3.2. Adicionar entradas de bloco de saída para o VCF conjunto

Adicione entradas para os arquivos VCF conjuntos.
Vamos colocá-los na raiz do diretório de resultados, pois esta é a saída final.

=== "Depois"

    ```groovy title="genomics.nf" hl_lines="11-16"
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
        joint_vcf {
            path '.'
        }
        joint_vcf_idx {
            path '.'
        }
    }
    ```

=== "Antes"

    ```groovy title="genomics.nf"
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

Com o processo, os alvos de publicação e o bloco de saída todos no lugar, podemos testar o fluxo de trabalho completo.

### 2.4. Executar o fluxo de trabalho

Execute o fluxo de trabalho para verificar se tudo funciona.

```bash
nextflow run genomics.nf -resume
```

??? success "Saída do comando"

    ```console
    N E X T F L O W   ~  version 25.10.2

    ┃ Launching `genomics.nf` [crazy_marconi] DSL2 - revision: 5da9afc841

    executor >  local (1)
    [9a/c7a873] SAMTOOLS_INDEX (2)       | 3 of 3, cached: 3 ✔
    [e4/4ed55e] GATK_HAPLOTYPECALLER (2) | 3 of 3, cached: 3 ✔
    [a6/7cc8ed] GATK_JOINTGENOTYPING     | 1 of 1 ✔
    ```

As duas primeiras etapas estão em cache da execução anterior, e a nova etapa `GATK_JOINTGENOTYPING` é executada uma vez nas entradas coletadas de todas as três amostras.
O arquivo de saída final, `family_trio.joint.vcf` (e seu índice), estão no diretório de resultados.

??? abstract "Conteúdo do diretório (links simbólicos encurtados)"

    ```console
    results/
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

Se você abrir o arquivo VCF conjunto, poderá verificar que o fluxo de trabalho produziu as chamadas de variantes esperadas.

```console title="family_trio.joint.vcf" linenums="40"
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	reads_father	reads_mother	reads_son
20_10037292_10066351	3480	.	C	CT	1625.89	.	AC=5;AF=0.833;AN=6;BaseQRankSum=0.220;DP=85;ExcessHet=0.0000;FS=2.476;MLEAC=5;MLEAF=0.833;MQ=60.00;MQRankSum=0.00;QD=21.68;ReadPosRankSum=-1.147e+00;SOR=0.487	GT:AD:DP:GQ:PL	0/1:15,16:31:99:367,0,375	1/1:0,18:18:54:517,54,0	1/1:0,26:26:78:756,78,0
20_10037292_10066351	3520	.	AT	A	1678.89	.	AC=5;AF=0.833;AN=6;BaseQRankSum=1.03;DP=80;ExcessHet=0.0000;FS=2.290;MLEAC=5;MLEAF=0.833;MQ=60.00;MQRankSum=0.00;QD=22.39;ReadPosRankSum=0.701;SOR=0.730	GT:AD:DP:GQ:PL	0/1:18,13:31:99:296,0,424	1/1:0,18:18:54:623,54,0	1/1:0,26:26:78:774,78,0
20_10037292_10066351	3529	.	T	A	154.29	.	AC=1;AF=0.167;AN=6;BaseQRankSum=-5.440e-01;DP=104;ExcessHet=0.0000;FS=1.871;MLEAC=1;MLEAF=0.167;MQ=60.00;MQRankSum=0.00;QD=7.71;ReadPosRankSum=-1.158e+00;SOR=1.034	GT:AD:DP:GQ:PL	0/0:44,0:44:99:0,112,1347	0/1:12,8:20:99:163,0,328	0/0:39,0:39:99:0,105,1194
```

Agora você tem um fluxo de trabalho de chamada conjunta de variantes automatizado e totalmente reproduzível!

!!! note

    Tenha em mente que os arquivos de dados que fornecemos cobrem apenas uma pequena porção do cromossomo 20.
    O tamanho real de um conjunto de chamadas de variantes seria contado em milhões de variantes.
    É por isso que usamos apenas pequenos subconjuntos de dados para fins de treinamento!

### Conclusão

Você sabe como coletar saídas de um canal e agrupá-las como uma única entrada para outro processo.
Você também sabe como construir uma linha de comando usando closures Groovy e como executar múltiplos comandos em um único processo.

### O que vem a seguir?

Dê um tapinha nas costas! Você completou o curso Nextflow for Genomics.

Vá para o [resumo final do curso](./next_steps.md) para revisar o que você aprendeu e descobrir o que vem a seguir.
