# Parte 3: Movendo código para módulos

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Tradução assistida por IA - [saiba mais e sugira melhorias](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

Na primeira parte deste curso, você construiu um pipeline de chamada de variantes que era completamente linear e processava os dados de cada amostra independentemente das outras.

Na segunda parte, mostramos como usar canais e operadores de canal para implementar chamada conjunta de variantes com GATK, construindo sobre o pipeline da Parte 1.

Nesta parte, mostraremos como converter o código desse fluxo de trabalho em módulos. Para acompanhar esta parte do treinamento, você deve ter completado a Parte 1 e a Parte 2, assim como [Hello Modules](../../../hello_nextflow/hello_modules.md), que cobre os fundamentos de módulos.

---

## 0. Aquecimento

Quando começamos a desenvolver nosso fluxo de trabalho, colocamos tudo em um único arquivo de código.
Agora é hora de abordar a **modularização** do nosso código, _ou seja_, extrair as definições de processo em módulos.

Vamos começar com o mesmo fluxo de trabalho da Parte 2, que fornecemos para você no arquivo `genomics-3.nf`.

!!! note "Nota"

     Certifique-se de estar no diretório de trabalho correto:
     `cd /workspaces/training/nf4-science/genomics`

Execute o fluxo de trabalho para verificar o ponto de partida:

```bash
nextflow run genomics-3.nf -resume
```

```console title="Saída"
 N E X T F L O W   ~  version 25.10.2

Launching `genomics-3.nf` [serene_borg] DSL2 - revision: 0cbebb67a1

executor >  local (7)
[6f/83ee72] SAMTOOLS_INDEX (3)       | 3 of 3 ✔
[53/b9d342] GATK_HAPLOTYPECALLER (1) | 3 of 3 ✔
[0c/fa6d15] GATK_JOINTGENOTYPING     | 1 of 1 ✔
```

Agora haverá um diretório `work` e um diretório `results_genomics` dentro do seu diretório do projeto.

### Conclusão

Você está pronto para começar a modularizar seu fluxo de trabalho.

### Qual é o próximo passo?

Mover os processos do fluxo de trabalho de Genômica para módulos.

---

## 1. Mover processos para módulos

Como você aprendeu em [Hello Modules](../../../hello_nextflow/hello_modules.md), você pode criar um módulo simplesmente copiando a definição do processo para seu próprio arquivo, em qualquer diretório, e pode nomear esse arquivo como quiser.

Por razões que ficarão claras mais tarde (em particular quando chegarmos aos testes), neste treinamento seguiremos a convenção de nomear o arquivo como `main.nf`, e colocá-lo em uma estrutura de diretórios nomeada após o kit de ferramentas e o comando.

### 1.1. Criar um módulo para o processo `SAMTOOLS_INDEX`

No caso do processo `SAMTOOLS_INDEX`, 'samtools' é o kit de ferramentas e 'index' é o comando. Então, criaremos uma estrutura de diretórios `modules/samtools/index` e colocaremos a definição do processo `SAMTOOLS_INDEX` no arquivo `main.nf` dentro desse diretório.

```bash
mkdir -p modules/samtools/index
touch modules/samtools/index/main.nf
```

Abra o arquivo `main.nf` e copie a definição do processo `SAMTOOLS_INDEX` para ele.

```groovy title="modules/samtools/index/main.nf" linenums="1"
/*
 * Gerar arquivo de índice BAM
 */
process SAMTOOLS_INDEX {

    container 'community.wave.seqera.io/library/samtools:1.20--b5dfbd93de237464'

    input:
    path input_bam

    output:
    tuple path(input_bam), path("${input_bam}.bai")

    script:
    """
    samtools index '$input_bam'
    """
}
```

Então, remova a definição do processo `SAMTOOLS_INDEX` de `genomics-3.nf`, e adicione uma declaração de importação para o módulo antes da próxima definição de processo, assim:

=== "Depois"

    ```groovy title="genomics-3.nf" linenums="1" hl_lines="1 2"
    // Incluir módulos
    include { SAMTOOLS_INDEX } from './modules/samtools/index/main.nf'

    /*
     * Chamar variantes com GATK HaplotypeCaller
     */
    process GATK_HAPLOTYPECALLER {
    ```

=== "Antes"

    ```groovy title="genomics-3.nf" linenums="1" hl_lines="1"
    /*
     * Chamar variantes com GATK HaplotypeCaller
     */
    process GATK_HAPLOTYPECALLER {
    ```

Agora você pode executar o fluxo de trabalho novamente, e ele deve funcionar da mesma forma que antes. Se você fornecer a flag `-resume`, nenhuma tarefa nova deve nem precisar ser executada:

```bash
nextflow run genomics-3.nf -resume
```

??? success "Saída do comando"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `genomics-3.nf` [sleepy_snyder] DSL2 - revision: aa68d06c43

    [0f/71b55e] SAMTOOLS_INDEX (1)       | 3 of 3, cached: 3 ✔
    [f1/18971b] GATK_HAPLOTYPECALLER (3) | 3 of 3, cached: 3 ✔
    [0c/fa6d15] GATK_JOINTGENOTYPING     | 1 of 1, cached: 1 ✔
    ```

### 1.2. Criar módulos para os processos `GATK_HAPLOTYPECALLER` e `GATK_JOINTGENOTYPING`

Repita os mesmos passos para os processos restantes.
Para cada processo:

1. Crie a estrutura de diretórios (`modules/gatk/haplotypecaller/` e `modules/gatk/jointgenotyping/`)
2. Crie um arquivo `main.nf` contendo a definição do processo
3. Remova a definição do processo de `genomics-3.nf`
4. Adicione uma declaração de importação para o módulo

Quando terminar, verifique se a estrutura de diretórios dos seus módulos está correta executando:

```bash
tree modules/
```

??? abstract "Conteúdo do diretório"

    ```console
    modules/
    ├── gatk
    │   ├── haplotypecaller
    │   │   └── main.nf
    │   └── jointgenotyping
    │       └── main.nf
    └── samtools
        └── index
            └── main.nf

    5 directories, 3 files
    ```

Você também deve ter algo assim no arquivo principal do fluxo de trabalho, após a seção de parâmetros:

```
include { SAMTOOLS_INDEX } from './modules/samtools/index/main.nf'
include { GATK_HAPLOTYPECALLER } from './modules/gatk/haplotypecaller/main.nf'
include { GATK_JOINTGENOTYPING } from './modules/gatk/jointgenotyping/main.nf'

workflow {
```

### Conclusão

Você praticou modularizar um fluxo de trabalho, com o fluxo de trabalho de genômica como exemplo.

### Qual é o próximo passo?

Testar o fluxo de trabalho modularizado.

---

## 2. Testar o fluxo de trabalho modularizado

Execute o fluxo de trabalho modularizado para verificar se tudo ainda funciona.

```bash
nextflow run genomics-3.nf -resume
```

```console title="Saída"
 N E X T F L O W   ~  version 25.10.2

Launching `genomics-3.nf` [astonishing_venter] DSL2 - revision: ca27264c13

[6f/83ee72] SAMTOOLS_INDEX (3)       | 3 of 3, cached: 3 ✔
[53/b9d342] GATK_HAPLOTYPECALLER (3) | 3 of 3, cached: 3 ✔
[0c/fa6d15] GATK_JOINTGENOTYPING     | 1 of 1, cached: 1 ✔
```

Tudo ainda funciona, incluindo a capacidade de retomada do pipeline.
Os resultados continuam sendo publicados no diretório `results_genomics`.

```console title="Conteúdo do diretório"
results_genomics/
├── family_trio.joint.vcf
├── family_trio.joint.vcf.idx
├── gvcf
│   ├── reads_father.bam.g.vcf
│   ├── reads_father.bam.g.vcf.idx
│   ├── reads_mother.bam.g.vcf
│   ├── reads_mother.bam.g.vcf.idx
│   ├── reads_son.bam.g.vcf
│   └── reads_son.bam.g.vcf.idx
└── indexed_bam
    ├── reads_father.bam
    ├── reads_father.bam.bai
    ├── reads_mother.bam
    ├── reads_mother.bam.bai
    ├── reads_son.bam
    └── reads_son.bam.bai
```

### Conclusão

Você modularizou um fluxo de trabalho e verificou que ele ainda funciona da mesma forma que antes.

### Qual é o próximo passo?

Revisar o que você aprendeu e olhar adiante para testes.

---

## 3. Resumo

Você modularizou o fluxo de trabalho, e nada mudou em como o pipeline funciona.
Isso é intencional: você reestruturou o código sem impactar sua função.

Os módulos contêm apenas a lógica do processo, tornando-os limpos e reutilizáveis.
O script principal controla o que é publicado e onde, enquanto os módulos permanecem focados em sua tarefa computacional.

Você estabeleceu uma base para coisas que tornarão seu código mais fácil de manter.
Por exemplo, agora você pode adicionar testes ao seu pipeline usando o framework nf-test.
É isso que veremos na próxima parte deste curso.
