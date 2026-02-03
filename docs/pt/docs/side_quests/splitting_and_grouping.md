# Divisão e Agrupamento

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Tradução assistida por IA - [saiba mais e sugira melhorias](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

O Nextflow fornece ferramentas poderosas para trabalhar com dados de forma flexível. Uma capacidade-chave é dividir dados em diferentes fluxos e depois agrupar itens relacionados de volta. Isso é especialmente valioso em fluxos de trabalho de bioinformática onde você precisa processar diferentes tipos de amostras separadamente antes de combinar resultados para análise.

Pense nisso como classificar correspondências: você separa cartas por destino, processa cada pilha de forma diferente e depois recombina itens indo para a mesma pessoa. O Nextflow usa operadores especiais para realizar isso com dados científicos. Essa abordagem também é comumente conhecida como o padrão **scatter/gather** em computação distribuída e fluxos de trabalho de bioinformática.

O sistema de canais do Nextflow está no centro dessa flexibilidade. Os canais conectam diferentes partes do seu fluxo de trabalho, permitindo que os dados fluam através da sua análise. Você pode criar múltiplos canais a partir de uma única fonte de dados, processar cada canal de forma diferente e depois mesclar canais de volta quando necessário. Essa abordagem permite que você projete fluxos de trabalho que naturalmente espelham os caminhos de ramificação e convergência de análises complexas de bioinformática.

### Objetivos de aprendizado

Nesta missão secundária, você aprenderá a dividir e agrupar dados usando os operadores de canal do Nextflow.
Começaremos com um arquivo CSV contendo informações de amostras e arquivos de dados associados, depois manipularemos e reorganizaremos esses dados.

Ao final desta missão secundária, você será capaz de separar e combinar fluxos de dados de forma eficaz, usando as seguintes técnicas:

- Ler dados de arquivos usando `splitCsv`
- Filtrar e transformar dados com `filter` e `map`
- Combinar dados relacionados usando `join` e `groupTuple`
- Criar combinações de dados com `combine` para processamento paralelo
- Otimizar a estrutura de dados usando `subMap` e estratégias de deduplicação
- Construir funções reutilizáveis com closures nomeados para ajudar a manipular estruturas de canal

Essas habilidades ajudarão você a construir fluxos de trabalho que podem lidar com múltiplos arquivos de entrada e diferentes tipos de dados de forma eficiente, mantendo uma estrutura de código limpa e de fácil manutenção.

### Pré-requisitos

Antes de começar esta missão secundária, você deve:

- Ter completado o tutorial [Hello Nextflow](../hello_nextflow/README.md) ou curso equivalente para iniciantes.
- Estar confortável usando conceitos e mecanismos básicos do Nextflow (processos, canais, operadores, trabalhando com arquivos, metadados)

**Opcional:** Recomendamos completar a missão secundária [Metadados em fluxos de trabalho](./metadata.md) primeiro.
Ela cobre os fundamentos de leitura de arquivos CSV com `splitCsv` e criação de mapas de metadados, que usaremos intensamente aqui.

---

## 0. Primeiros passos

#### Abrir o codespace de treinamento

Se você ainda não o fez, certifique-se de abrir o ambiente de treinamento conforme descrito em [Configuração do Ambiente](../envsetup/index.md).

[![Open in GitHub Codespaces](https://github.com/codespaces/badge.svg)](https://codespaces.new/nextflow-io/training?quickstart=1&ref=master)

#### Mover para o diretório do projeto

Vamos mover para o diretório onde os arquivos deste tutorial estão localizados.

```bash
cd side-quests/splitting_and_grouping
```

Você pode configurar o VSCode para focar neste diretório:

```bash
code .
```

#### Revisar os materiais

Você encontrará um arquivo de fluxo de trabalho principal e um diretório `data` contendo uma planilha de amostras chamada `samplesheet.csv`.

```console title="Conteúdo do diretório"
.
├── data
│   └── samplesheet.csv
└── main.nf
```

A planilha de amostras contém informações sobre amostras de diferentes pacientes, incluindo o ID do paciente, número de repetição da amostra, tipo (normal ou tumor) e caminhos para arquivos de dados hipotéticos (que não existem de fato, mas vamos fingir que existem).

```console title="samplesheet.csv"
id,repeat,type,bam
patientA,1,normal,patientA_rep1_normal.bam
patientA,1,tumor,patientA_rep1_tumor.bam
patientA,2,normal,patientA_rep2_normal.bam
patientA,2,tumor,patientA_rep2_tumor.bam
patientB,1,normal,patientB_rep1_normal.bam
patientB,1,tumor,patientB_rep1_tumor.bam
patientC,1,normal,patientC_rep1_normal.bam
patientC,1,tumor,patientC_rep1_tumor.bam
```

Esta planilha lista oito amostras de três pacientes (A, B, C).

Para cada paciente, temos amostras que são do tipo `tumor` (tipicamente originadas de biópsias de tumor) ou `normal` (coletadas de tecido saudável ou sangue).
Se você não está familiarizado com análise de câncer, saiba apenas que isso corresponde a um modelo experimental que usa amostras pareadas tumor/normal para realizar análises contrastivas.

Para o paciente A especificamente, temos dois conjuntos de réplicas técnicas (repetições).

!!! note "Nota"

    Não se preocupe se você não está familiarizado com este desenho experimental, não é crítico para entender este tutorial.

#### Revisar a tarefa

Seu desafio é escrever um fluxo de trabalho Nextflow que irá:

1. **Ler** dados de amostras de um arquivo CSV e estruturá-los com mapas de metadados
2. **Separar** amostras em diferentes canais com base no tipo (normal vs tumor)
3. **Unir** pares correspondentes tumor/normal por ID do paciente e número de réplica
4. **Distribuir** amostras através de intervalos genômicos para processamento paralelo
5. **Agrupar** amostras relacionadas de volta para análise downstream

Isso representa um padrão comum de bioinformática onde você precisa dividir dados para processamento independente, depois recombinar itens relacionados para análise comparativa.

#### Lista de verificação de prontidão

Acha que está pronto para começar?

- [ ] Entendo o objetivo deste curso e seus pré-requisitos
- [ ] Meu codespace está funcionando
- [ ] Defini meu diretório de trabalho apropriadamente
- [ ] Entendo a tarefa

Se você pode marcar todas as caixas, está pronto para começar.

---

## 1. Ler dados de amostras

### 1.1. Ler dados de amostras com `splitCsv` e criar mapas de metadados

Vamos começar lendo os dados de amostras com `splitCsv` e organizando-os no padrão de mapa de metadados. No `main.nf`, você verá que já começamos o fluxo de trabalho.

```groovy title="main.nf" linenums="1" hl_lines="2"
workflow {
    ch_samplesheet = channel.fromPath("./data/samplesheet.csv")
}
```

!!! note "Nota"

    Ao longo deste tutorial, usaremos o prefixo `ch_` para todas as variáveis de canal para indicar claramente que são canais Nextflow.

Se você completou a missão secundária [Metadados em fluxos de trabalho](./metadata.md), você reconhecerá esse padrão. Usaremos `splitCsv` para ler o CSV e imediatamente estruturar os dados com um mapa de metadados para separar metadados de caminhos de arquivo.

!!! info

    Encontraremos dois conceitos diferentes chamados `map` neste treinamento:

    - **Estrutura de dados**: O mapa Groovy (equivalente a dicionários/hashes em outras linguagens) que armazena pares chave-valor
    - **Operador de canal**: O operador `.map()` que transforma itens em um canal

    Esclareceremos qual deles queremos dizer no contexto, mas essa distinção é importante para entender ao trabalhar com Nextflow.

Aplique estas mudanças ao `main.nf`:

=== "Depois"

    ```groovy title="main.nf" linenums="2" hl_lines="2-6"
        ch_samples = channel.fromPath("./data/samplesheet.csv")
            .splitCsv(header: true)
            .map{ row ->
              [[id:row.id, repeat:row.repeat, type:row.type], row.bam]
            }
            .view()
    ```

=== "Antes"

    ```groovy title="main.nf" linenums="2" hl_lines="1"
        ch_samplesheet = channel.fromPath("./data/samplesheet.csv")
    ```

Isso combina a operação `splitCsv` (lendo o CSV com cabeçalhos) e a operação `map` (estruturando dados como tuplas `[meta, file]`) em uma etapa. Aplique essa mudança e execute o pipeline:

```bash
nextflow run main.nf
```

??? success "Saída do comando"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `main.nf` [deadly_mercator] DSL2 - revision: bd6b0224e9

    [[id:patientA, repeat:1, type:normal], patientA_rep1_normal.bam]
    [[id:patientA, repeat:1, type:tumor], patientA_rep1_tumor.bam]
    [[id:patientA, repeat:2, type:normal], patientA_rep2_normal.bam]
    [[id:patientA, repeat:2, type:tumor], patientA_rep2_tumor.bam]
    [[id:patientB, repeat:1, type:normal], patientB_rep1_normal.bam]
    [[id:patientB, repeat:1, type:tumor], patientB_rep1_tumor.bam]
    [[id:patientC, repeat:1, type:normal], patientC_rep1_normal.bam]
    [[id:patientC, repeat:1, type:tumor], patientC_rep1_tumor.bam]
    ```

Agora temos um canal onde cada item é uma tupla `[meta, file]` - metadados separados de caminhos de arquivo. Esta estrutura nos permite dividir e agrupar nossa carga de trabalho com base em campos de metadados.

---

## 2. Filtrar e transformar dados

### 2.1. Filtrar dados com `filter`

Podemos usar o [operador `filter`](https://www.nextflow.io/docs/latest/operator.html#filter) para filtrar os dados com base em uma condição. Digamos que queremos processar apenas amostras normais. Podemos fazer isso filtrando os dados com base no campo `type`. Vamos inserir isso antes do operador `view`.

=== "Depois"

    ```groovy title="main.nf" linenums="2" hl_lines="6"
        ch_samples = channel.fromPath("./data/samplesheet.csv")
            .splitCsv(header: true)
            .map{ row ->
              [[id:row.id, repeat:row.repeat, type:row.type], row.bam]
            }
            .filter { meta, file -> meta.type == 'normal' }
            .view()
    ```

=== "Antes"

    ```groovy title="main.nf" linenums="2"
        ch_samples = channel.fromPath("./data/samplesheet.csv")
            .splitCsv(header: true)
            .map{ row ->
              [[id:row.id, repeat:row.repeat, type:row.type], row.bam]
            }
            .view()
    ```

Execute o fluxo de trabalho novamente para ver o resultado filtrado:

```bash
nextflow run main.nf
```

??? success "Saída do comando"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `main.nf` [admiring_brown] DSL2 - revision: 194d61704d

    [[id:patientA, repeat:1, type:normal], patientA_rep1_normal.bam]
    [[id:patientA, repeat:2, type:normal], patientA_rep2_normal.bam]
    [[id:patientB, repeat:1, type:normal], patientB_rep1_normal.bam]
    [[id:patientC, repeat:1, type:normal], patientC_rep1_normal.bam]
    ```

Filtramos com sucesso os dados para incluir apenas amostras normais. Vamos recapitular como isso funciona.

O operador `filter` recebe um closure que é aplicado a cada elemento no canal. Se o closure retornar `true`, o elemento é incluído; se retornar `false`, o elemento é excluído.

No nosso caso, queremos manter apenas amostras onde `meta.type == 'normal'`. O closure usa a tupla `meta,file` para se referir a cada amostra, acessa o tipo de amostra com `meta.type` e verifica se é igual a `'normal'`.

Isso é realizado com o único closure que introduzimos acima:

```groovy title="main.nf" linenums="7"
    .filter { meta, file -> meta.type == 'normal' }
```

### 2.2. Criar canais filtrados separados

Atualmente estamos aplicando o filtro ao canal criado diretamente do CSV, mas queremos filtrar isso de mais de uma maneira, então vamos reescrever a lógica para criar um canal filtrado separado para amostras normais:

=== "Depois"

    ```groovy title="main.nf" linenums="2" hl_lines="6 8"
        ch_samples = channel.fromPath("./data/samplesheet.csv")
            .splitCsv(header: true)
            .map{ row ->
                [[id:row.id, repeat:row.repeat, type:row.type], row.bam]
            }
        ch_normal_samples = ch_samples
            .filter { meta, file -> meta.type == 'normal' }
        ch_normal_samples
            .view()
    ```

=== "Antes"

    ```groovy title="main.nf" linenums="2"
        ch_samples = channel.fromPath("./data/samplesheet.csv")
            .splitCsv(header: true)
            .map{ row ->
              [[id:row.id, repeat:row.repeat, type:row.type], row.bam]
            }
            .filter { meta, file -> meta.type == 'normal' }
            .view()
    ```

Execute o pipeline para ver os resultados:

```bash
nextflow run main.nf
```

??? success "Saída do comando"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `main.nf` [trusting_poisson] DSL2 - revision: 639186ee74

    [[id:patientA, repeat:1, type:normal], patientA_rep1_normal.bam]
    [[id:patientA, repeat:2, type:normal], patientA_rep2_normal.bam]
    [[id:patientB, repeat:1, type:normal], patientB_rep1_normal.bam]
    [[id:patientC, repeat:1, type:normal], patientC_rep1_normal.bam]
    ```

Filtramos com sucesso os dados e criamos um canal separado para amostras normais.

Vamos criar um canal filtrado para as amostras de tumor também:

=== "Depois"

    ```groovy title="main.nf" linenums="7" hl_lines="3-8"
        ch_normal_samples = ch_samples
            .filter { meta, file -> meta.type == 'normal' }
        ch_tumor_samples = ch_samples
            .filter { meta, file -> meta.type == 'tumor' }
        ch_normal_samples
            .view{'Normal sample: ' + it}
        ch_tumor_samples
            .view{'Tumor sample: ' + it}
    ```

=== "Antes"

    ```groovy title="main.nf" linenums="7" hl_lines="3 4"
        ch_normal_samples = ch_samples
            .filter { meta, file -> meta.type == 'normal' }
        ch_normal_samples
            .view()
    ```

```bash
nextflow run main.nf
```

??? success "Saída do comando"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `main.nf` [maniac_boltzmann] DSL2 - revision: 3636b6576b

    Tumor sample: [[id:patientA, repeat:1, type:tumor], patientA_rep1_tumor.bam]
    Tumor sample: [[id:patientA, repeat:2, type:tumor], patientA_rep2_tumor.bam]
    Normal sample: [[id:patientA, repeat:1, type:normal], patientA_rep1_normal.bam]
    Normal sample: [[id:patientA, repeat:2, type:normal], patientA_rep2_normal.bam]
    Normal sample: [[id:patientB, repeat:1, type:normal], patientB_rep1_normal.bam]
    Normal sample: [[id:patientC, repeat:1, type:normal], patientC_rep1_normal.bam]
    Tumor sample: [[id:patientB, repeat:1, type:tumor], patientB_rep1_tumor.bam]
    Tumor sample: [[id:patientC, repeat:1, type:tumor], patientC_rep1_tumor.bam]
    ```

Separamos as amostras normais e de tumor em dois canais diferentes e usamos um closure fornecido a `view()` para rotulá-las de forma diferente na saída: `ch_tumor_samples.view{'Tumor sample: ' + it}`.

### Conclusão

Nesta seção, você aprendeu:

- **Filtrar dados**: Como filtrar dados com `filter`
- **Dividir dados**: Como dividir dados em diferentes canais com base em uma condição
- **Visualizar dados**: Como usar `view` para imprimir os dados e rotular a saída de diferentes canais

Agora separamos as amostras normais e de tumor em dois canais diferentes. Em seguida, vamos unir as amostras normais e de tumor no campo `id`.

---

## 3. Unir canais por identificadores

Na seção anterior, separamos as amostras normais e de tumor em dois canais diferentes. Elas poderiam ser processadas independentemente usando processos ou fluxos de trabalho específicos com base em seu tipo. Mas o que acontece quando queremos comparar as amostras normais e de tumor do mesmo paciente? Neste ponto, precisamos uni-las de volta, certificando-nos de combinar as amostras com base em seu campo `id`.

O Nextflow inclui muitos métodos para combinar canais, mas neste caso o operador mais apropriado é [`join`](https://www.nextflow.io/docs/latest/operator.html#join). Se você está familiarizado com SQL, ele age como a operação `JOIN`, onde especificamos a chave para unir e o tipo de junção a ser executada.

### 3.1. Usar `map` e `join` para combinar com base no ID do paciente

Se verificarmos a documentação do [`join`](https://www.nextflow.io/docs/latest/operator.html#join), podemos ver que por padrão ele une dois canais com base no primeiro item em cada tupla.

#### 3.1.1. Verificar a estrutura de dados

Se você não tem a saída do console ainda disponível, vamos executar o pipeline para verificar nossa estrutura de dados e ver como precisamos modificá-la para unir no campo `id`.

```bash
nextflow run main.nf
```

??? success "Saída do comando"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `main.nf` [maniac_boltzmann] DSL2 - revision: 3636b6576b

    Tumor sample: [[id:patientA, repeat:1, type:tumor], patientA_rep1_tumor.bam]
    Tumor sample: [[id:patientA, repeat:2, type:tumor], patientA_rep2_tumor.bam]
    Normal sample: [[id:patientA, repeat:1, type:normal], patientA_rep1_normal.bam]
    Normal sample: [[id:patientA, repeat:2, type:normal], patientA_rep2_normal.bam]
    Normal sample: [[id:patientB, repeat:1, type:normal], patientB_rep1_normal.bam]
    Normal sample: [[id:patientC, repeat:1, type:normal], patientC_rep1_normal.bam]
    Tumor sample: [[id:patientB, repeat:1, type:tumor], patientB_rep1_tumor.bam]
    Tumor sample: [[id:patientC, repeat:1, type:tumor], patientC_rep1_tumor.bam]
    ```

Podemos ver que o campo `id` é o primeiro elemento em cada mapa de metadados. Para que `join` funcione, devemos isolar o campo `id` em cada tupla. Depois disso, podemos simplesmente usar o operador `join` para combinar os dois canais.

#### 3.1.2. Isolar o campo `id`

Para isolar o campo `id`, podemos usar o [operador `map`](https://www.nextflow.io/docs/latest/operator.html#map) para criar uma nova tupla com o campo `id` como o primeiro elemento.

=== "Depois"

    ```groovy title="main.nf" linenums="7" hl_lines="3 6"
        ch_normal_samples = ch_samples
            .filter { meta, file -> meta.type == 'normal' }
            .map { meta, file -> [meta.id, meta, file] }
        ch_tumor_samples = ch_samples
            .filter { meta, file -> meta.type == 'tumor' }
            .map { meta, file -> [meta.id, meta, file] }
        ch_normal_samples
            .view{'Normal sample: ' + it}
        ch_tumor_samples
            .view{'Tumor sample: ' + it}
    ```

=== "Antes"

    ```groovy title="main.nf" linenums="7"
        ch_normal_samples = ch_samples
            .filter { meta, file -> meta.type == 'normal' }
        ch_tumor_samples = ch_samples
            .filter { meta, file -> meta.type == 'tumor' }
        ch_normal_samples
            .view{'Normal sample: ' + it}
        ch_tumor_samples
            .view{'Tumor sample: ' + it}
    ```

```bash
nextflow run main.nf
```

??? success "Saída do comando"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `main.nf` [mad_lagrange] DSL2 - revision: 9940b3f23d

    Tumor sample: [patientA, [id:patientA, repeat:1, type:tumor], patientA_rep1_tumor.bam]
    Tumor sample: [patientA, [id:patientA, repeat:2, type:tumor], patientA_rep2_tumor.bam]
    Normal sample: [patientA, [id:patientA, repeat:1, type:normal], patientA_rep1_normal.bam]
    Normal sample: [patientA, [id:patientA, repeat:2, type:normal], patientA_rep2_normal.bam]
    Tumor sample: [patientB, [id:patientB, repeat:1, type:tumor], patientB_rep1_tumor.bam]
    Tumor sample: [patientC, [id:patientC, repeat:1, type:tumor], patientC_rep1_tumor.bam]
    Normal sample: [patientB, [id:patientB, repeat:1, type:normal], patientB_rep1_normal.bam]
    Normal sample: [patientC, [id:patientC, repeat:1, type:normal], patientC_rep1_normal.bam]
    ```

Pode ser sutil, mas você deve ser capaz de ver que o primeiro elemento em cada tupla é o campo `id`.

#### 3.1.3. Combinar os dois canais

Agora podemos usar o operador `join` para combinar os dois canais com base no campo `id`.

Mais uma vez, usaremos `view` para imprimir as saídas unidas.

=== "Depois"

    ```groovy title="main.nf" linenums="7" hl_lines="7-9"
        ch_normal_samples = ch_samples
            .filter { meta, file -> meta.type == 'normal' }
            .map { meta, file -> [meta.id, meta, file] }
        ch_tumor_samples = ch_samples
            .filter { meta, file -> meta.type == 'tumor' }
            .map { meta, file -> [meta.id, meta, file] }
        ch_joined_samples = ch_normal_samples
            .join(ch_tumor_samples)
        ch_joined_samples.view()
    ```

=== "Antes"

    ```groovy title="main.nf" linenums="7" hl_lines="7-10"
        ch_normal_samples = ch_samples
            .filter { meta, file -> meta.type == 'normal' }
            .map { meta, file -> [meta.id, meta, file] }
        ch_tumor_samples = ch_samples
            .filter { meta, file -> meta.type == 'tumor' }
            .map { meta, file -> [meta.id, meta, file] }
        ch_normal_samples
            .view{'Normal sample: ' + it}
        ch_tumor_samples
            .view{'Tumor sample: ' + it}
    ```

```bash
nextflow run main.nf
```

??? success "Saída do comando"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `main.nf` [soggy_wiles] DSL2 - revision: 3bc1979889

    [patientA, [id:patientA, repeat:1, type:normal], patientA_rep1_normal.bam, [id:patientA, repeat:1, type:tumor], patientA_rep1_tumor.bam]
    [patientA, [id:patientA, repeat:2, type:normal], patientA_rep2_normal.bam, [id:patientA, repeat:2, type:tumor], patientA_rep2_tumor.bam]
    [patientB, [id:patientB, repeat:1, type:normal], patientB_rep1_normal.bam, [id:patientB, repeat:1, type:tumor], patientB_rep1_tumor.bam]
    [patientC, [id:patientC, repeat:1, type:normal], patientC_rep1_normal.bam, [id:patientC, repeat:1, type:tumor], patientC_rep1_tumor.bam]
    ```

É um pouco difícil de ver porque é muito largo, mas você deve ser capaz de ver que as amostras foram unidas pelo campo `id`. Cada tupla agora tem o formato:

- `id`: O ID da amostra
- `normal_meta_map`: Os metadados da amostra normal incluindo tipo, réplica e caminho para o arquivo bam
- `normal_sample_file`: O arquivo da amostra normal
- `tumor_meta_map`: Os metadados da amostra de tumor incluindo tipo, réplica e caminho para o arquivo bam
- `tumor_sample`: A amostra de tumor incluindo tipo, réplica e caminho para o arquivo bam

!!! warning "Aviso"

    O operador `join` descartará quaisquer tuplas não correspondidas. Neste exemplo, garantimos que todas as amostras fossem correspondidas para tumor e normal, mas se isso não for verdade você deve usar o parâmetro `remainder: true` para manter as tuplas não correspondidas. Consulte a [documentação](https://www.nextflow.io/docs/latest/operator.html#join) para mais detalhes.

Então agora você sabe como usar `map` para isolar um campo em uma tupla, e como usar `join` para combinar tuplas com base no primeiro campo.
Com esse conhecimento, podemos combinar com sucesso canais com base em um campo compartilhado.

Em seguida, consideraremos a situação onde você quer unir em múltiplos campos.

### 3.2. Unir em múltiplos campos

Temos 2 réplicas para a sampleA, mas apenas 1 para sampleB e sampleC. Neste caso fomos capazes de uni-las efetivamente usando o campo `id`, mas o que aconteceria se estivessem fora de sincronia? Poderíamos misturar as amostras normais e de tumor de diferentes réplicas!

Para evitar isso, podemos unir em múltiplos campos. Na verdade, existem várias maneiras de fazer isso, mas vamos focar em criar uma nova chave de junção que inclui tanto o `id` da amostra quanto o número de `replicate`.

Vamos começar criando uma nova chave de junção. Podemos fazer isso da mesma forma que antes, usando o [operador `map`](https://www.nextflow.io/docs/latest/operator.html#map) para criar uma nova tupla com os campos `id` e `repeat` como o primeiro elemento.

=== "Depois"

    ```groovy title="main.nf" linenums="7" hl_lines="3 6"
        ch_normal_samples = ch_samples
            .filter { meta, file -> meta.type == 'normal' }
            .map { meta, file -> [[meta.id, meta.repeat], meta, file] }
        ch_tumor_samples = ch_samples
            .filter { meta, file -> meta.type == 'tumor' }
            .map { meta, file -> [[meta.id, meta.repeat], meta, file] }
    ```

=== "Antes"

    ```groovy title="main.nf" linenums="7" hl_lines="3 6"
        ch_normal_samples = ch_samples
            .filter { meta, file -> meta.type == 'normal' }
            .map { meta, file -> [meta.id, meta, file] }
        ch_tumor_samples = ch_samples
            .filter { meta, file -> meta.type == 'tumor' }
            .map { meta, file -> [meta.id, meta, file] }
    ```

Agora devemos ver a junção ocorrendo, mas usando tanto os campos `id` quanto `repeat`. Execute o fluxo de trabalho:

```bash
nextflow run main.nf
```

??? success "Saída do comando"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `main.nf` [prickly_wing] DSL2 - revision: 3bebf22dee

    [[patientA, 1], [id:patientA, repeat:1, type:normal], patientA_rep1_normal.bam, [id:patientA, repeat:1, type:tumor], patientA_rep1_tumor.bam]
    [[patientA, 2], [id:patientA, repeat:2, type:normal], patientA_rep2_normal.bam, [id:patientA, repeat:2, type:tumor], patientA_rep2_tumor.bam]
    [[patientB, 1], [id:patientB, repeat:1, type:normal], patientB_rep1_normal.bam, [id:patientB, repeat:1, type:tumor], patientB_rep1_tumor.bam]
    [[patientC, 1], [id:patientC, repeat:1, type:normal], patientC_rep1_normal.bam, [id:patientC, repeat:1, type:tumor], patientC_rep1_tumor.bam]
    ```

Note como temos uma tupla de dois elementos (campos `id` e `repeat`) como o primeiro elemento de cada resultado unido. Isso demonstra como itens complexos podem ser usados como uma chave de junção, permitindo correspondências bastante intrincadas entre amostras das mesmas condições.

Se você quiser explorar mais maneiras de unir em chaves diferentes, consulte a [documentação do operador join](https://www.nextflow.io/docs/latest/operator.html#join) para opções e exemplos adicionais.

### 3.3. Usar `subMap` para criar uma nova chave de junção

A abordagem anterior perde os nomes dos campos da nossa chave de junção - os campos `id` e `repeat` se tornam apenas uma lista de valores. Para reter os nomes dos campos para acesso posterior, podemos usar o [método `subMap`](<https://docs.groovy-lang.org/latest/html/groovy-jdk/java/util/Map.html#subMap(java.util.Collection)>).

O método `subMap` extrai apenas os pares chave-valor especificados de um mapa. Aqui extrairemos apenas os campos `id` e `repeat` para criar nossa chave de junção.

=== "Depois"

    ```groovy title="main.nf" linenums="7" hl_lines="3 6"
        ch_normal_samples = ch_samples
            .filter { meta, file -> meta.type == 'normal' }
            .map { meta, file -> [meta.subMap(['id', 'repeat']), meta, file] }
        ch_tumor_samples = ch_samples
            .filter { meta, file -> meta.type == 'tumor' }
            .map { meta, file -> [meta.subMap(['id', 'repeat']), meta, file] }
    ```

=== "Antes"

    ```groovy title="main.nf" linenums="7" hl_lines="3 6"
        ch_normal_samples = ch_samples
            .filter { meta, file -> meta.type == 'normal' }
            .map { meta, file -> [[meta.id, meta.repeat], meta, file] }
        ch_tumor_samples = ch_samples
            .filter { meta, file -> meta.type == 'tumor' }
            .map { meta, file -> [[meta.id, meta.repeat], meta, file] }
    ```

```bash
nextflow run main.nf
```

??? success "Saída do comando"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `main.nf` [reverent_wing] DSL2 - revision: 847016c3b7

    [[id:patientA, repeat:1], [id:patientA, repeat:1, type:normal], patientA_rep1_normal.bam, [id:patientA, repeat:1, type:tumor], patientA_rep1_tumor.bam]
    [[id:patientA, repeat:2], [id:patientA, repeat:2, type:normal], patientA_rep2_normal.bam, [id:patientA, repeat:2, type:tumor], patientA_rep2_tumor.bam]
    [[id:patientB, repeat:1], [id:patientB, repeat:1, type:normal], patientB_rep1_normal.bam, [id:patientB, repeat:1, type:tumor], patientB_rep1_tumor.bam]
    [[id:patientC, repeat:1], [id:patientC, repeat:1, type:normal], patientC_rep1_normal.bam, [id:patientC, repeat:1, type:tumor], patientC_rep1_tumor.bam]
    ```

Agora temos uma nova chave de junção que não apenas inclui os campos `id` e `repeat`, mas também retém os nomes dos campos para que possamos acessá-los mais tarde por nome, por exemplo `meta.id` e `meta.repeat`.

### 3.4. Usar um closure nomeado em map

Para evitar duplicação e reduzir erros, podemos usar um closure nomeado. Um closure nomeado nos permite criar uma função reutilizável que podemos chamar em vários lugares.

Para fazer isso, primeiro definimos o closure como uma nova variável:

=== "Depois"

    ```groovy title="main.nf" linenums="2" hl_lines="7"
        ch_samples = channel.fromPath("./data/samplesheet.csv")
            .splitCsv(header: true)
            .map{ row ->
                [[id:row.id, repeat:row.repeat, type:row.type], row.bam]
            }

        getSampleIdAndReplicate = { meta, bam -> [ meta.subMap(['id', 'repeat']), meta, file(bam) ] }

        ch_normal_samples = ch_samples
            .filter { meta, file -> meta.type == 'normal' }
    ```

=== "Antes"

    ```groovy title="main.nf" linenums="2"
        ch_samples = channel.fromPath("./data/samplesheet.csv")
            .splitCsv(header: true)
            .map{ row ->
                [[id:row.id, repeat:row.repeat, type:row.type], row.bam]
            }
        ch_normal_samples = ch_samples
            .filter { meta, file -> meta.type == 'normal' }
    ```

Definimos a transformação do map como uma variável nomeada que podemos reutilizar.

Note que também convertemos o caminho do arquivo para um objeto Path usando `file()` para que qualquer processo recebendo este canal possa lidar com o arquivo corretamente (para mais informações veja [Trabalhando com arquivos](./working_with_files.md)).

Vamos implementar o closure em nosso fluxo de trabalho:

=== "Depois"

    ```groovy title="main.nf" linenums="10" hl_lines="3 6"
        ch_normal_samples = ch_samples
            .filter { meta, file -> meta.type == 'normal' }
             .map ( getSampleIdAndReplicate )
        ch_tumor_samples = ch_samples
            .filter { meta, file -> meta.type == 'tumor' }
             .map ( getSampleIdAndReplicate )

    ```

=== "Antes"

    ```groovy title="main.nf" linenums="10" hl_lines="3 6"
        ch_normal_samples = ch_samples
            .filter { meta, file -> meta.type == 'normal' }
            .map { meta, file -> [meta.subMap(['id', 'repeat']), meta, file] }
        ch_tumor_samples = ch_samples
            .filter { meta, file -> meta.type == 'tumor' }
            .map { meta, file -> [meta.subMap(['id', 'repeat']), meta, file] }
    ```

!!! note "Nota"

    O operador `map` mudou de usar `{ }` para usar `( )` para passar o closure como um argumento. Isso ocorre porque o operador `map` espera um closure como argumento e `{ }` é usado para definir um closure anônimo. Ao chamar um closure nomeado, use a sintaxe `( )`.

Execute o fluxo de trabalho mais uma vez para verificar se tudo ainda está funcionando:

```bash
nextflow run main.nf
```

??? success "Saída do comando"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `main.nf` [angry_meninsky] DSL2 - revision: 2edc226b1d

    [[id:patientA, repeat:1], [id:patientA, repeat:1, type:normal], patientA_rep1_normal.bam, [id:patientA, repeat:1, type:tumor], patientA_rep1_tumor.bam]
    [[id:patientA, repeat:2], [id:patientA, repeat:2, type:normal], patientA_rep2_normal.bam, [id:patientA, repeat:2, type:tumor], patientA_rep2_tumor.bam]
    [[id:patientB, repeat:1], [id:patientB, repeat:1, type:normal], patientB_rep1_normal.bam, [id:patientB, repeat:1, type:tumor], patientB_rep1_tumor.bam]
    [[id:patientC, repeat:1], [id:patientC, repeat:1, type:normal], patientC_rep1_normal.bam, [id:patientC, repeat:1, type:tumor], patientC_rep1_tumor.bam]
    ```

Usar um closure nomeado nos permite reutilizar a mesma transformação em vários lugares, reduzindo o risco de erros e tornando o código mais legível e de fácil manutenção.

### 3.5. Reduzir duplicação de dados

Temos muitos dados duplicados em nosso fluxo de trabalho. Cada item nas amostras unidas repete os campos `id` e `repeat`. Como essa informação já está disponível na chave de agrupamento, podemos evitar essa redundância. Como lembrete, nossa estrutura de dados atual se parece com isso:

```groovy
[
  [
    "id": "sampleC",
    "repeat": "1",
  ],
  [
    "id": "sampleC",
    "repeat": "1",
    "type": "normal",
  ],
  "sampleC_rep1_normal.bam"
  [
    "id": "sampleC",
    "repeat": "1",
    "type": "tumor",
  ],
  "sampleC_rep1_tumor.bam"
]
```

Como os campos `id` e `repeat` estão disponíveis na chave de agrupamento, vamos removê-los do resto de cada item do canal para evitar duplicação. Podemos fazer isso usando o método `subMap` para criar um novo mapa com apenas o campo `type`. Essa abordagem nos permite manter todas as informações necessárias enquanto eliminamos redundância em nossa estrutura de dados.

=== "Depois"

    ```groovy title="main.nf" linenums="8" hl_lines="1"
        getSampleIdAndReplicate = { meta, bam -> [ meta.subMap(['id', 'repeat']), meta.subMap(['type']), file(bam) ] }
    ```

=== "Antes"

    ```groovy title="main.nf" linenums="8" hl_lines="1"
        getSampleIdAndReplicate = { meta, bam -> [ meta.subMap(['id', 'repeat']), meta, file(bam) ] }
    ```

Agora o closure retorna uma tupla onde o primeiro elemento contém os campos `id` e `repeat`, e o segundo elemento contém apenas o campo `type`. Isso elimina redundância armazenando as informações de `id` e `repeat` uma vez na chave de agrupamento, mantendo todas as informações necessárias.

Execute o fluxo de trabalho para ver como isso se parece:

```bash
nextflow run main.nf
```

??? success "Saída do comando"

    ```console
    [[id:patientA, repeat:1], [type:normal], /workspaces/training/side-quests/splitting_and_grouping/patientA_rep1_normal.bam, [type:tumor], /workspaces/training/side-quests/splitting_and_grouping/patientA_rep1_tumor.bam]
    [[id:patientA, repeat:2], [type:normal], /workspaces/training/side-quests/splitting_and_grouping/patientA_rep2_normal.bam, [type:tumor], /workspaces/training/side-quests/splitting_and_grouping/patientA_rep2_tumor.bam]
    [[id:patientB, repeat:1], [type:normal], /workspaces/training/side-quests/splitting_and_grouping/patientB_rep1_normal.bam, [type:tumor], /workspaces/training/side-quests/splitting_and_grouping/patientB_rep1_tumor.bam]
    [[id:patientC, repeat:1], [type:normal], /workspaces/training/side-quests/splitting_and_grouping/patientC_rep1_normal.bam, [type:tumor], /workspaces/training/side-quests/splitting_and_grouping/patientC_rep1_tumor.bam]
    ```

Podemos ver que só declaramos os campos `id` e `repeat` uma vez na chave de agrupamento e temos o campo `type` nos dados da amostra. Não perdemos nenhuma informação, mas conseguimos tornar o conteúdo do nosso canal mais sucinto.

### 3.6. Remover informações redundantes

Removemos informações duplicadas acima, mas ainda temos algumas outras informações redundantes em nossos canais.

No início, separamos as amostras normais e de tumor usando `filter`, depois as unimos com base nas chaves `id` e `repeat`. O operador `join` preserva a ordem em que as tuplas são mescladas, então no nosso caso, com amostras normais no lado esquerdo e amostras de tumor no lado direito, o canal resultante mantém esta estrutura: `id, <elementos normais>, <elementos tumor>`.

Como sabemos a posição de cada elemento em nosso canal, podemos simplificar ainda mais a estrutura removendo os metadados `[type:normal]` e `[type:tumor]`.

=== "Depois"

    ```groovy title="main.nf" linenums="8" hl_lines="1"
        getSampleIdAndReplicate = { meta, file -> [ meta.subMap(['id', 'repeat']), file ] }
    ```

=== "Antes"

    ```groovy title="main.nf" linenums="8" hl_lines="1"
        getSampleIdAndReplicate = { meta, file -> [ meta.subMap(['id', 'repeat']), meta.subMap(['type']), file ] }
    ```

Execute novamente para ver o resultado:

```bash
nextflow run main.nf
```

??? success "Saída do comando"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `main.nf` [confident_leavitt] DSL2 - revision: a2303895bd

    [[id:patientA, repeat:1], patientA_rep1_normal.bam, patientA_rep1_tumor.bam]
    [[id:patientA, repeat:2], patientA_rep2_normal.bam, patientA_rep2_tumor.bam]
    [[id:patientB, repeat:1], patientB_rep1_normal.bam, patientB_rep1_tumor.bam]
    [[id:patientC, repeat:1], patientC_rep1_normal.bam, patientC_rep1_tumor.bam]
    ```

### Conclusão

Nesta seção, você aprendeu:

- **Manipular Tuplas**: Como usar `map` para isolar um campo em uma tupla
- **Unir Tuplas**: Como usar `join` para combinar tuplas com base no primeiro campo
- **Criar Chaves de Junção**: Como usar `subMap` para criar uma nova chave de junção
- **Closures Nomeados**: Como usar um closure nomeado em map
- **Junção em Múltiplos Campos**: Como unir em múltiplos campos para correspondência mais precisa
- **Otimização da Estrutura de Dados**: Como simplificar a estrutura do canal removendo informações redundantes

Você agora tem um fluxo de trabalho que pode dividir uma planilha de amostras, filtrar as amostras normais e de tumor, uni-las por ID de amostra e número de réplica, depois imprimir os resultados.

Este é um padrão comum em fluxos de trabalho de bioinformática onde você precisa combinar amostras ou outros tipos de dados após processamento independente, então é uma habilidade útil. Em seguida, vamos olhar para repetir uma amostra várias vezes.

## 4. Distribuir amostras em intervalos

Um padrão-chave em fluxos de trabalho de bioinformática é distribuir análises através de regiões genômicas. Por exemplo, a chamada de variantes pode ser paralelizada dividindo o genoma em intervalos (como cromossomos ou regiões menores). Essa estratégia de paralelização melhora significativamente a eficiência do pipeline distribuindo a carga computacional entre múltiplos núcleos ou nós, reduzindo o tempo total de execução.

Na seção seguinte, demonstraremos como distribuir nossos dados de amostras através de múltiplos intervalos genômicos. Vamos parear cada amostra com cada intervalo, permitindo processamento paralelo de diferentes regiões genômicas. Isso multiplicará o tamanho do nosso conjunto de dados pelo número de intervalos, criando múltiplas unidades de análise independentes que podem ser reunidas mais tarde.

### 4.1. Distribuir amostras em intervalos usando `combine`

Vamos começar criando um canal de intervalos. Para manter a vida simples, usaremos apenas 3 intervalos que definiremos manualmente. Em um fluxo de trabalho real, você poderia lê-los de uma entrada de arquivo ou até criar um canal com vários arquivos de intervalo.

=== "Depois"

    ```groovy title="main.nf" linenums="17" hl_lines="2"
            .join(ch_tumor_samples)
        ch_intervals = channel.of('chr1', 'chr2', 'chr3')
    ```

=== "Antes"

    ```groovy title="main.nf" linenums="17" hl_lines="2"
            .join(ch_tumor_samples)
        ch_joined_samples.view()
    ```

Agora lembre-se, queremos repetir cada amostra para cada intervalo. Isso às vezes é chamado de produto cartesiano das amostras e intervalos. Podemos conseguir isso usando o [operador `combine`](https://www.nextflow.io/docs/latest/operator.html#combine). Isso pegará cada item do canal 1 e o repetirá para cada item no canal 2. Vamos adicionar um operador combine ao nosso fluxo de trabalho:

=== "Depois"

    ```groovy title="main.nf" linenums="18" hl_lines="3-5"
        ch_intervals = channel.of('chr1', 'chr2', 'chr3')

        ch_combined_samples = ch_joined_samples
            .combine(ch_intervals)
            .view()
    ```

=== "Antes"

    ```groovy title="main.nf" linenums="18"
        ch_intervals = channel.of('chr1', 'chr2', 'chr3')
    ```

Agora vamos executá-lo e ver o que acontece:

```bash
nextflow run main.nf
```

??? success "Saída do comando"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `main.nf` [mighty_tesla] DSL2 - revision: ae013ab70b

    [[id:patientA, repeat:1], patientA_rep1_normal.bam, patientA_rep1_tumor.bam, chr1]
    [[id:patientA, repeat:1], patientA_rep1_normal.bam, patientA_rep1_tumor.bam, chr2]
    [[id:patientA, repeat:1], patientA_rep1_normal.bam, patientA_rep1_tumor.bam, chr3]
    [[id:patientA, repeat:2], patientA_rep2_normal.bam, patientA_rep2_tumor.bam, chr1]
    [[id:patientA, repeat:2], patientA_rep2_normal.bam, patientA_rep2_tumor.bam, chr2]
    [[id:patientA, repeat:2], patientA_rep2_normal.bam, patientA_rep2_tumor.bam, chr3]
    [[id:patientB, repeat:1], patientB_rep1_normal.bam, patientB_rep1_tumor.bam, chr1]
    [[id:patientB, repeat:1], patientB_rep1_normal.bam, patientB_rep1_tumor.bam, chr2]
    [[id:patientB, repeat:1], patientB_rep1_normal.bam, patientB_rep1_tumor.bam, chr3]
    [[id:patientC, repeat:1], patientC_rep1_normal.bam, patientC_rep1_tumor.bam, chr1]
    [[id:patientC, repeat:1], patientC_rep1_normal.bam, patientC_rep1_tumor.bam, chr2]
    [[id:patientC, repeat:1], patientC_rep1_normal.bam, patientC_rep1_tumor.bam, chr3]
    ```

Sucesso! Repetimos cada amostra para cada intervalo em nossa lista de 3 intervalos. Efetivamente triplicamos o número de itens em nosso canal.

É um pouco difícil de ler, então na próxima seção vamos organizá-lo.

### 4.2. Organizar o canal

Podemos usar o operador `map` para organizar e refatorar nossos dados de amostras para que seja mais fácil de entender. Vamos mover a string de intervalos para o mapa de junção no primeiro elemento.

=== "Depois"

    ```groovy title="main.nf" linenums="20" hl_lines="3-9"
        ch_combined_samples = ch_joined_samples
            .combine(ch_intervals)
            .map { grouping_key, normal, tumor, interval ->
                [
                    grouping_key + [interval: interval],
                    normal,
                    tumor
                ]
            }
            .view()
    ```

=== "Antes"

    ```groovy title="main.nf" linenums="20"
        ch_combined_samples = ch_joined_samples
            .combine(ch_intervals)
            .view()
    ```

Vamos analisar o que essa operação map faz passo a passo.

Primeiro, usamos parâmetros nomeados para tornar o código mais legível. Ao usar os nomes `grouping_key`, `normal`, `tumor` e `interval`, podemos nos referir aos elementos na tupla por nome em vez de por índice:

```groovy
        .map { grouping_key, normal, tumor, interval ->
```

Em seguida, combinamos o `grouping_key` com o campo `interval`. O `grouping_key` é um mapa contendo os campos `id` e `repeat`. Criamos um novo mapa com o `interval` e os mesclamos usando a adição de mapas do Groovy (`+`):

```groovy
                grouping_key + [interval: interval],
```

Finalmente, retornamos isso como uma tupla com três elementos: o mapa de metadados combinado, o arquivo da amostra normal e o arquivo da amostra de tumor:

```groovy
            [
                grouping_key + [interval: interval],
                normal,
                tumor
            ]
```

Vamos executá-lo novamente e verificar o conteúdo do canal:

```bash
nextflow run main.nf
```

??? success "Saída do comando"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `main.nf` [sad_hawking] DSL2 - revision: 1f6f6250cd

    [[id:patientA, repeat:1, interval:chr1], patientA_rep1_normal.bam, patientA_rep1_tumor.bam]
    [[id:patientA, repeat:1, interval:chr2], patientA_rep1_normal.bam, patientA_rep1_tumor.bam]
    [[id:patientA, repeat:1, interval:chr3], patientA_rep1_normal.bam, patientA_rep1_tumor.bam]
    [[id:patientA, repeat:2, interval:chr1], patientA_rep2_normal.bam, patientA_rep2_tumor.bam]
    [[id:patientA, repeat:2, interval:chr2], patientA_rep2_normal.bam, patientA_rep2_tumor.bam]
    [[id:patientA, repeat:2, interval:chr3], patientA_rep2_normal.bam, patientA_rep2_tumor.bam]
    [[id:patientB, repeat:1, interval:chr1], patientB_rep1_normal.bam, patientB_rep1_tumor.bam]
    [[id:patientB, repeat:1, interval:chr2], patientB_rep1_normal.bam, patientB_rep1_tumor.bam]
    [[id:patientB, repeat:1, interval:chr3], patientB_rep1_normal.bam, patientB_rep1_tumor.bam]
    [[id:patientC, repeat:1, interval:chr1], patientC_rep1_normal.bam, patientC_rep1_tumor.bam]
    [[id:patientC, repeat:1, interval:chr2], patientC_rep1_normal.bam, patientC_rep1_tumor.bam]
    [[id:patientC, repeat:1, interval:chr3], patientC_rep1_normal.bam, patientC_rep1_tumor.bam]
    ```

Usar `map` para forçar seus dados na estrutura correta pode ser complicado, mas é crucial para manipulação eficaz de dados.

Agora temos cada amostra repetida em todos os intervalos genômicos, criando múltiplas unidades de análise independentes que podem ser processadas em paralelo. Mas e se quisermos trazer amostras relacionadas de volta? Na próxima seção, aprenderemos como agrupar amostras que compartilham atributos comuns.

### Conclusão

Nesta seção, você aprendeu:

- **Distribuir amostras em intervalos**: Como usar `combine` para repetir amostras em intervalos
- **Criar produtos cartesianos**: Como gerar todas as combinações de amostras e intervalos
- **Organizar estrutura de canal**: Como usar `map` para reestruturar dados para melhor legibilidade
- **Preparação para processamento paralelo**: Como configurar dados para análise distribuída

## 5. Agregar amostras usando `groupTuple`

Nas seções anteriores, aprendemos como dividir dados de um arquivo de entrada e filtrar por campos específicos (no nosso caso amostras normais e de tumor). Mas isso cobre apenas um tipo de junção. E se quisermos agrupar amostras por um atributo específico? Por exemplo, em vez de unir pares correspondentes normal-tumor, podemos querer processar todas as amostras de "sampleA" juntas independentemente de seu tipo. Esse padrão é comum em fluxos de trabalho de bioinformática onde você pode querer processar amostras relacionadas separadamente por razões de eficiência antes de comparar ou combinar os resultados no final.

O Nextflow inclui métodos integrados para fazer isso, o principal que vamos olhar é `groupTuple`.

Vamos começar agrupando todas as nossas amostras que têm os mesmos campos `id` e `interval`, isso seria típico de uma análise onde quiséssemos agrupar réplicas técnicas mas manter amostras significativamente diferentes separadas.

Para fazer isso, devemos separar nossas variáveis de agrupamento para que possamos usá-las isoladamente.

O primeiro passo é similar ao que fizemos na seção anterior. Devemos isolar nossa variável de agrupamento como o primeiro elemento da tupla. Lembre-se, nosso primeiro elemento é atualmente um mapa dos campos `id`, `repeat` e `interval`:

```groovy title="main.nf" linenums="1"
{
  "id": "sampleA",
  "repeat": "1",
  "interval": "chr1"
}
```

Podemos reutilizar o método `subMap` de antes para isolar nossos campos `id` e `interval` do mapa. Como antes, usaremos o operador `map` para aplicar o método `subMap` ao primeiro elemento da tupla para cada amostra.

=== "Depois"

    ```groovy title="main.nf" linenums="20" hl_lines="11-19"
        ch_combined_samples = ch_joined_samples
            .combine(ch_intervals)
            .map { grouping_key, normal, tumor, interval ->
                [
                    grouping_key + [interval: interval],
                    normal,
                    tumor
                ]
            }

        ch_grouped_samples = ch_combined_samples
            .map { grouping_key, normal, tumor ->
                [
                    grouping_key.subMap('id', 'interval'),
                    normal,
                    tumor
                ]
              }
              .view()
    ```

=== "Antes"

    ```groovy title="main.nf" linenums="20" hl_lines="10"
        ch_combined_samples = ch_joined_samples
            .combine(ch_intervals)
            .map { grouping_key, normal, tumor, interval ->
                [
                    grouping_key + [interval: interval],
                    normal,
                    tumor
                ]
            }
            .view()
    ```

Vamos executá-lo novamente e verificar o conteúdo do canal:

```bash
nextflow run main.nf
```

??? success "Saída do comando"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `main.nf` [hopeful_brenner] DSL2 - revision: 7f4f7fea76

    [[id:patientA, interval:chr1], patientA_rep1_normal.bam, patientA_rep1_tumor.bam]
    [[id:patientA, interval:chr2], patientA_rep1_normal.bam, patientA_rep1_tumor.bam]
    [[id:patientA, interval:chr3], patientA_rep1_normal.bam, patientA_rep1_tumor.bam]
    [[id:patientA, interval:chr1], patient
