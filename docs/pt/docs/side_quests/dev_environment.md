# Ambiente de Desenvolvimento

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Tradução assistida por IA - [saiba mais e sugira melhorias](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

Ambientes de Desenvolvimento Integrados (IDEs) modernos podem transformar dramaticamente sua experiência de desenvolvimento com Nextflow. Este side quest foca especificamente em aproveitar o VS Code e sua extensão Nextflow para escrever código mais rapidamente, detectar erros precocemente e navegar workflows complexos de forma eficiente.

!!! note "Este não é um tutorial tradicional"

    Ao contrário de outros módulos de treinamento, este guia é organizado como uma coleção de dicas rápidas, sugestões e exemplos práticos em vez de um tutorial passo a passo. Cada seção pode ser explorada independentemente com base em seus interesses e necessidades atuais de desenvolvimento. Sinta-se à vontade para pular entre seções e focar nos recursos que serão mais imediatamente úteis para o desenvolvimento do seu fluxo de trabalho.

## O que você deve saber primeiro

Este guia assume que você completou o curso de treinamento [Hello Nextflow](../hello_nextflow/) e está confortável com conceitos fundamentais do Nextflow, incluindo:

- **Estrutura básica de fluxo de trabalho**: Compreensão de processos, fluxos de trabalho e como eles se conectam
- **Operações de canal**: Criação de canais, passagem de dados entre processos e uso de operadores básicos
- **Módulos e organização**: Criação de módulos reutilizáveis e uso de declarações include
- **Fundamentos de configuração**: Uso do `nextflow.config` para parâmetros, diretivas de processo e perfis

## O que você aprenderá aqui

Este guia foca em **recursos de produtividade do IDE** que tornarão você um desenvolvedor Nextflow mais eficiente:

- **Destaque de sintaxe avançado**: Compreender o que o VS Code está mostrando sobre a estrutura do seu código
- **Autocompletar inteligente**: Aproveitar sugestões sensíveis ao contexto para escrever código mais rapidamente
- **Detecção de erros e diagnósticos**: Capturar erros de sintaxe antes de executar seu fluxo de trabalho
- **Navegação de código**: Mover-se rapidamente entre processos, módulos e definições
- **Formatação e organização**: Manter estilo de código consistente e legível
- **Desenvolvimento assistido por IA** (opcional): Usar ferramentas modernas de IA integradas ao seu IDE

!!! info "Por que recursos de IDE agora?"

    Você provavelmente já estava usando o VS Code durante o curso [Hello Nextflow](../hello_nextflow/), mas mantivemos o foco no aprendizado dos fundamentos do Nextflow em vez dos recursos do IDE. Agora que você está confortável com conceitos básicos do Nextflow como processos, fluxos de trabalho, canais e módulos, você está pronto para aproveitar os recursos sofisticados do IDE que o tornarão um desenvolvedor mais eficiente.

    Pense nisso como "subir de nível" seu ambiente de desenvolvimento - o mesmo editor que você tem usado possui capacidades muito mais poderosas que se tornam verdadeiramente valiosas uma vez que você entende com o que elas estão te ajudando.

---

## 0. Configuração e Aquecimento

Vamos configurar um espaço de trabalho especificamente para explorar recursos do IDE:

```bash title="Navegue para o diretório de recursos do IDE"
cd side-quests/ide_features
```

Abra este diretório no VS Code:

```bash title="Abra o VS Code no diretório atual"
code .
```

O diretório `ide_features` contém fluxos de trabalho de exemplo que demonstram vários recursos do IDE:

```bash title="Mostre a estrutura de diretórios"
tree .
```

```console title="Estrutura do projeto"
tree .
.
├── basic_workflow.nf
├── complex_workflow.nf
├── data
│   ├── sample_001.fastq.gz
│   ├── sample_002.fastq.gz
│   ├── sample_003.fastq.gz
│   ├── sample_004.fastq.gz
│   ├── sample_005.fastq.gz
│   └── sample_data.csv
├── modules
│   ├── fastqc.nf
│   ├── star.nf
│   └── utils.nf
└── nextflow.config

3 directories, 12 files
```

!!! note "Sobre os Arquivos de Exemplo"

    - `basic_workflow.nf` é um fluxo de trabalho básico funcional que você pode executar e modificar
    - `complex_workflow.nf` é projetado apenas para ilustração para demonstrar recursos de navegação - pode não executar com sucesso, mas mostra uma estrutura de fluxo de trabalho realista com múltiplos arquivos

### Atalhos de Teclado

Alguns dos recursos neste guia usarão atalhos de teclado opcionais. Você pode estar acessando este material via GitHub Codespaces no navegador, e neste caso às vezes os atalhos não funcionarão como esperado porque são usados para outras coisas no seu sistema.

Se você estiver executando o VS Code localmente, como provavelmente estará quando realmente escrever fluxos de trabalho, os atalhos funcionarão conforme descrito.

Se você estiver usando um Mac, alguns (não todos) atalhos de teclado usarão "cmd" em vez de "ctrl", e indicaremos isso no texto como `Ctrl/Cmd`.

### 0.1. Instalando a Extensão Nextflow

!!! note "Já Usando Devcontainers?"

    Se você estiver trabalhando no **GitHub Codespaces** ou usando um **devcontainer local**, a extensão Nextflow provavelmente já está instalada e configurada para você. Você pode pular as etapas de instalação manual abaixo e prosseguir diretamente para explorar os recursos da extensão.

Para instalar a extensão manualmente:

1. Abra o VS Code
2. Vá para a visualização de Extensões clicando no ícone de extensões à esquerda: ![ícone de extensões](img/extensions_icon.png) (atalho `Ctrl/Cmd+Shift+X` se você estiver executando o VSCode localmente)
3. Procure por "Nextflow"
4. Instale a extensão oficial do Nextflow

![Instalar Extensão Nextflow](img/install_extension.png)

### 0.2. Layout do Espaço de Trabalho

Como você tem usado o VS Code durante o Hello Nextflow, você já está familiarizado com o básico. Aqui está como organizar seu espaço de trabalho eficientemente para esta sessão:

- **Área do Editor**: Para visualizar e editar arquivos. Você pode dividir isso em vários painéis para comparar arquivos lado a lado.
- **Explorador de Arquivos** clique (![ícone do explorador de arquivos](img/files_icon.png)) (`Ctrl/Cmd+Shift+E`): Os arquivos e pastas locais no seu sistema. Mantenha isso aberto à esquerda para navegar entre arquivos
- **Terminal Integrado** (`Ctrl+Shift+` crase para Windows e MacOS): Um terminal para interagir com o computador na parte inferior. Use isso para executar Nextflow ou outros comandos.
- **Painel de Problemas** (`Ctrl+Shift+M`): O VS Code mostrará quaisquer erros e problemas que detectar aqui. Isso é útil para destacar questões rapidamente.

Você pode arrastar painéis ou ocultá-los (`Ctrl/Cmd+B` para alternar a barra lateral) para personalizar seu layout enquanto trabalhamos nos exemplos.

### Conclusão

Você tem o VS Code configurado com a extensão Nextflow e entende o layout do espaço de trabalho para desenvolvimento eficiente.

### Qual é o próximo passo?

Aprenda como o destaque de sintaxe ajuda você a entender a estrutura do código Nextflow rapidamente.

---

## 1. Destaque de Sintaxe e Estrutura de Código

Agora que seu espaço de trabalho está configurado, vamos explorar como o destaque de sintaxe do VS Code ajuda você a ler e escrever código Nextflow de forma mais eficaz.

### 1.1. Elementos de Sintaxe Nextflow

Abra `basic_workflow.nf` para ver o destaque de sintaxe em ação:

![Demonstração de Sintaxe](img/syntax_showcase.png)

Observe como o VS Code destaca:

- **Palavras-chave** (`process`, `workflow`, `input`, `output`, `script`) em cores distintas
- **Literais de string** e **parâmetros** com estilos diferentes
- **Comentários** em uma cor suave
- **Variáveis** e **chamadas de função** com ênfase apropriada
- **Blocos de código** com guias de indentação adequadas

!!! note "Cores Dependentes do Tema"

    As cores específicas que você vê dependerão do seu tema do VS Code (modo escuro/claro), configurações de cores e quaisquer personalizações que você tenha feito. O importante é que diferentes elementos de sintaxe sejam visualmente distinguidos uns dos outros, tornando a estrutura do código mais fácil de entender independentemente do esquema de cores escolhido.

### 1.2. Compreendendo a Estrutura de Código

O destaque de sintaxe ajuda você a identificar rapidamente:

- **Limites de processo**: Distinção clara entre diferentes processos
- **Blocos de entrada/saída**: Fácil de identificar definições de fluxo de dados
- **Blocos de script**: Os comandos reais sendo executados
- **Operações de canal**: Etapas de transformação de dados
- **Diretivas de configuração**: Configurações específicas do processo

Esta organização visual se torna inestimável ao trabalhar com fluxos de trabalho complexos contendo múltiplos processos e fluxos de dados intrincados.

### Conclusão

Você entende como o destaque de sintaxe do VS Code ajuda você a ler a estrutura do código Nextflow e identificar diferentes elementos da linguagem para um desenvolvimento mais rápido.

### Qual é o próximo passo?

Aprenda como o autocompletar inteligente acelera a escrita de código com sugestões sensíveis ao contexto.

---

## 2. Autocompletar Inteligente

Os recursos de autocompletar do VS Code ajudam você a escrever código mais rapidamente e com menos erros, sugerindo opções apropriadas baseadas no contexto.

### 2.1. Sugestões Sensíveis ao Contexto

As opções de autocompletar variam dependendo de onde você está no seu código:

#### Operações de Canal

Abra `basic_workflow.nf` novamente e tente digitar `channel.` no bloco workflow:

![Autocompletar de canal](img/autocomplete_channel.png)

Você verá sugestões para:

- `fromPath()` - Criar canal a partir de caminhos de arquivo
- `fromFilePairs()` - Criar canal a partir de arquivos pareados
- `of()` - Criar canal a partir de valores
- `fromSRA()` - Criar canal a partir de acessos SRA
- E muito mais...

Isso ajuda você a encontrar rapidamente a factory de canal certa para usar sem precisar lembrar nomes exatos de métodos.

Você também pode descobrir os operadores disponíveis para aplicar a canais. Por exemplo, digite `FASTQC.out.html.` para ver operações disponíveis:

![Autocompletar de operações de canal](img/autocomplete_operators.png)

#### Diretivas de Processo

Dentro de um bloco de script de processo, digite `task.` para ver propriedades de runtime disponíveis:

![Autocompletar de propriedades de task](img/autocomplete_task.png)

#### Configuração

Abra nextflow.config e digite `process.` em qualquer lugar para ver diretivas de processo disponíveis:

![Autocompletar de configuração](img/autocomplete_config.png)

Você verá sugestões para:

- `executor`
- `memory`
- `cpus`

Isso economiza tempo ao configurar processos e funciona em diferentes escopos de configuração. Por exemplo, tente digitar `docker.` para ver opções de configuração específicas do Docker.

### Conclusão

Você pode usar o autocompletar inteligente do VS Code para descobrir operações de canal disponíveis, diretivas de processo e opções de configuração sem memorizar sintaxe.

### Qual é o próximo passo?

Aprenda como a detecção de erros em tempo real ajuda você a capturar problemas antes de executar seu fluxo de trabalho, simplesmente lendo o código.

## 3. Detecção de Erros e Diagnósticos

A detecção de erros em tempo real do VS Code ajuda você a capturar problemas antes de executar seu fluxo de trabalho.

### 3.1. Detecção de Erros de Sintaxe

Vamos criar um erro deliberado para ver a detecção em ação. Abra `basic_workflow.nf` e mude o nome do processo de `FASTQC` para `FASTQ` (ou qualquer outro nome inválido). O VS Code imediatamente destacará o erro no bloco workflow com um sublinhado ondulado vermelho:

![Sublinhado de erro](img/error_underline.png)

### 3.2. Painel de Problemas

Além do destaque individual de erros, o VS Code fornece um Painel de Problemas centralizado que agrega todos os erros, avisos e mensagens de informação em todo o seu espaço de trabalho. Abra-o com `Ctrl/Cmd+Shift+M` e use o ícone de filtro para mostrar apenas erros relevantes ao arquivo atual:

![Filtrar o painel de problemas](img/active_file.png)

Clique em qualquer problema para ir diretamente à linha problemática

![Painel de Problemas](img/problems_panel.png)

Corrija o erro mudando o nome do processo de volta para `FASTQC`.

### 3.3. Padrões de Erros Comuns

Erros comuns na sintaxe Nextflow incluem:

- **Colchetes faltando**: `{` ou `}` não correspondidos
- **Blocos incompletos**: Seções obrigatórias ausentes em processos
- **Sintaxe inválida**: DSL Nextflow mal formado
- **Erros de digitação em palavras-chave**: Diretivas de processo escritas incorretamente
- **Incompatibilidades de canal**: Incompatibilidades de tipo

O language server do Nextflow destaca esses problemas no Painel de Problemas. Você pode verificá-los cedo para evitar erros de sintaxe ao executar um pipeline.

### Conclusão

Você pode usar a detecção de erros do VS Code e o Painel de Problemas para capturar erros de sintaxe e problemas antes de executar seu fluxo de trabalho, economizando tempo e evitando frustrações.

### Qual é o próximo passo?

Aprenda como navegar eficientemente entre processos, módulos e definições em fluxos de trabalho complexos.

---

## 4. Navegação de Código e Gerenciamento de Símbolos

Navegação eficiente é crucial ao trabalhar com fluxos de trabalho complexos que abrangem múltiplos arquivos. Para entender isso, substitua a definição do processo em `basic_workflow.nf` com uma importação do módulo que fornecemos:

=== "Depois"

    ```groovy title="basic_workflow.nf" linenums="3"
    include { FASTQC } from './modules/fastqc.nf'
    ```

=== "Antes"

    ```groovy title="basic_workflow.nf" linenums="3"
    process FASTQC {
        tag "${sample_id}"
        publishDir "${params.output_dir}/fastqc", mode: 'copy'

        input:
        tuple val(sample_id), path(reads)

        output:
        tuple val(sample_id), path("*.html"), emit: html
        tuple val(sample_id), path("*.zip"), emit: zip

        script:
        def args = task.ext.args ?: ''
        """
        fastqc \\
            ${args} \\
            --threads ${task.cpus} \\
            ${reads}
        """
    }
    ```

### 4.1. Ir para Definição

Se você passar o mouse sobre um nome de processo como `FASTQC`, você verá um popup com a interface do módulo (entradas e saídas):

![Ir para definição](img/syntax.png)

Este recurso é particularmente valioso ao criar fluxos de trabalho, pois permite que você entenda a interface do módulo sem abrir o arquivo do módulo diretamente.

Você pode navegar rapidamente para qualquer definição de processo, módulo ou variável usando **Ctrl/Cmd-click**. Passe o mouse sobre o link para o arquivo do módulo no topo do script e siga o link conforme sugerido:

![Seguir link](img/follow_link.png)

A mesma coisa funciona para nomes de processo. Volte para `basic_workflow.nf` e tente isso no nome do processo `FASTQC` no bloco workflow. Isso leva você diretamente ao nome do processo (que é o mesmo que o arquivo do módulo neste exemplo, mas poderia estar no meio de um arquivo muito maior).

Para voltar para onde você estava, use **Alt+←** (ou **Ctrl+-** no Mac). Esta é uma maneira poderosa de explorar código sem perder sua posição.

Agora vamos explorar a navegação em um fluxo de trabalho mais complexo usando `complex_workflow.nf` (o arquivo somente para ilustração mencionado anteriormente). Este fluxo de trabalho contém múltiplos processos definidos em arquivos de módulo separados, além de alguns inline. Embora estruturas complexas de múltiplos arquivos possam ser desafiadoras para navegar manualmente, a capacidade de pular para definições torna a exploração muito mais gerenciável.

1. Abra `complex_workflow.nf`
2. Navegue para definições de módulo
3. Use **Alt+←** (ou **Ctrl+-**) para navegar de volta
4. Navegue para o nome do processo `FASTQC` no bloco workflow. Isso leva você diretamente ao nome do processo (que é o mesmo que o arquivo do módulo neste exemplo, mas poderia estar no meio de um arquivo muito maior).
5. Navegue de volta novamente
6. Navegue para o processo `TRIM_GALORE` no bloco workflow. Este é definido inline, então não levará você a um arquivo separado, mas ainda mostrará a definição do processo, e você ainda pode navegar de volta para onde estava.

### 4.2. Navegação de Símbolos

Com `complex_workflow.nf` ainda aberto, você pode obter uma visão geral de todos os símbolos no arquivo digitando `@` na barra de pesquisa no topo do VSCode (o atalho de teclado é `Ctrl/Cmd+Shift+O`, mas pode não funcionar no Codespaces). Isso abre o painel de navegação de símbolos, que lista todos os símbolos no arquivo atual:

![Navegação de símbolos](img/symbols.png)

Isso mostra:

- Todas as definições de processo
- Definições de workflow (há dois workflows definidos neste arquivo)
- Definições de função

Comece a digitar para filtrar resultados.

### 4.3. Encontrar Todas as Referências

Entender onde um processo ou variável é usado em toda a sua base de código pode ser muito útil. Por exemplo, se você quiser encontrar todas as referências ao processo `FASTQC`, comece navegando até sua definição. Você pode fazer isso abrindo `modules/fastqc.nf` diretamente, ou usando o recurso de navegação rápida do VS Code com `Ctrl/Cmd-click` como fizemos acima. Uma vez na definição do processo, clique com o botão direito no nome do processo `FASTQC` e selecione "Find All References" no menu de contexto para ver todas as instâncias onde ele é usado.

![Encontrar referências](img/references.png)

Este recurso exibe todas as instâncias onde `FASTQC` é referenciado dentro do seu espaço de trabalho, incluindo seu uso nos dois workflows distintos. Esta percepção é crucial para avaliar o impacto potencial de modificações no processo `FASTQC`.

### 4.4. Painel Outline

O painel Outline, localizado na barra lateral Explorer (clique ![ícone do Explorer](img/files_icon.png)), fornece uma visão geral conveniente de todos os símbolos no seu arquivo atual. Este recurso permite que você navegue e gerencie rapidamente a estrutura do seu código exibindo funções, variáveis e outros elementos-chave em uma visualização hierárquica.

![Painel Outline](img/outline.png)

Use o painel Outline para navegar rapidamente para diferentes partes do seu código sem usar o navegador de arquivos.

### 4.5. Visualização DAG

A extensão Nextflow do VS Code pode visualizar seu fluxo de trabalho como um Grafo Acíclico Direcionado (DAG). Isso ajuda você a entender o fluxo de dados e dependências entre processos. Abra `complex_workflow.nf` e clique no botão "Preview DAG" acima de `workflow {` (o segundo bloco `workflow` neste arquivo):

![Visualização DAG](img/dag_preview.png)

Este é apenas o workflow de 'entrada', mas você também pode visualizar o DAG para os workflows internos clicando no botão "Preview DAG" acima do workflow `RNASEQ_PIPELINE {` mais acima:

![Visualização DAG workflow interno](img/dag_preview_inner.png)

Para este fluxo de trabalho, você pode usar os nós no DAG para navegar até as definições de processo correspondentes no código. Clique em um nó e ele levará você à definição de processo relevante no editor. Particularmente quando um fluxo de trabalho cresce para um tamanho grande, isso pode realmente ajudá-lo a navegar pelo código e entender como os processos estão conectados.

### Conclusão

Você pode navegar fluxos de trabalho complexos eficientemente usando ir-para-definição, busca de símbolos, encontrar referências e visualização DAG para entender a estrutura do código e dependências.

### Qual é o próximo passo?

Aprenda como trabalhar efetivamente com múltiplos arquivos interconectados em projetos Nextflow maiores.

## 5. Trabalhando com Múltiplos Arquivos

O desenvolvimento real com Nextflow envolve trabalhar com múltiplos arquivos interconectados. Vamos explorar como o VS Code ajuda você a gerenciar projetos complexos eficientemente.

### 5.1. Navegação Rápida de Arquivos

Com `complex_workflow.nf` aberto, você notará que ele importa vários módulos. Vamos praticar navegação rápida entre eles.

Pressione **Ctrl+P** (ou **Cmd+P**) e comece a digitar "fast":

O VS Code mostrará arquivos correspondentes. Selecione `modules/fastqc.nf` para ir lá instantaneamente. Isso é muito mais rápido do que clicar no explorador de arquivos quando você sabe aproximadamente qual arquivo está procurando.

Tente isso com outros padrões:

- Digite "star" para encontrar o arquivo do módulo de alinhamento STAR (`star.nf`)
- Digite "utils" para encontrar o arquivo de funções utilitárias (`utils.nf`)
- Digite "config" para ir para arquivos de configuração (`nextflow.config`)

### 5.2. Editor Dividido para Desenvolvimento Multi-arquivo

Ao trabalhar com módulos, você frequentemente precisa ver tanto o fluxo de trabalho principal quanto as definições de módulo simultaneamente. Vamos configurar isso:

1. Abra `complex_workflow.nf`
2. Abra `modules/fastqc.nf` em uma nova aba
3. Clique com o botão direito na aba `modules/fastqc.nf` e selecione "Split Right"
4. Agora você pode ver ambos os arquivos lado a lado

![Editor dividido](img/split_editor.png)

Isso é inestimável quando:

- Verificando interfaces de módulo ao escrever chamadas de workflow, e a visualização não é suficiente
- Comparando processos similares em diferentes módulos
- Depurando fluxo de dados entre workflow e módulos

### 5.3. Pesquisa em Todo o Projeto

Às vezes você precisa encontrar onde padrões específicos são usados em todo o seu projeto. Pressione `Ctrl/Cmd+Shift+F` para abrir o painel de pesquisa.

Tente pesquisar por `publishDir` em todo o espaço de trabalho:

![Pesquisa no projeto](img/project_search.png)

Isso mostra todo arquivo que usa diretórios de publicação, ajudando você a:

- Entender padrões de organização de saída
- Encontrar exemplos de diretivas específicas
- Garantir consistência entre módulos

### Conclusão

Você pode gerenciar projetos complexos de múltiplos arquivos usando navegação rápida de arquivos, editores divididos e pesquisa em todo o projeto para trabalhar eficientemente com fluxos de trabalho e módulos.

### Qual é o próximo passo?

Aprenda como recursos de formatação de código e manutenção mantêm seus fluxos de trabalho organizados e legíveis.

---

## 6. Formatação de Código e Manutenção

Formatação adequada de código é essencial não apenas para estética, mas também para melhorar a legibilidade, compreensão e facilidade de atualização de fluxos de trabalho complexos.

### 6.1. Formatação Automática em Ação

Abra `basic_workflow.nf` e deliberadamente bagunce a formatação:

- Remova alguma indentação: Destaque o documento inteiro e pressione `shift+tab` várias vezes para remover o máximo de indentações possível.
- Adicione espaços extras em lugares aleatórios: na declaração `channel.fromPath`, adicione 30 espaços após o `(`.
- Quebre algumas linhas de forma estranha: Adicione uma nova linha entre o operador `.view {` e a string `Processing sample:`, mas não adicione uma nova linha correspondente antes do parêntese de fechamento `}`.

Agora pressione `Shift+Alt+F` (ou `Shift+Option+F` no MacOS) para formatar automaticamente:

O VS Code imediatamente:

- Corrige a indentação para mostrar a estrutura do processo claramente
- Alinha elementos similares de forma consistente
- Remove espaços em branco desnecessários
- Mantém quebras de linha legíveis

Note que a formatação automática pode não resolver todos os problemas de estilo de código. O language server do Nextflow visa manter seu código organizado, mas também respeita suas preferências pessoais em certas áreas. Por exemplo, se você remover a indentação dentro do bloco `script` de um processo, o formatador deixará como está, já que você pode intencionalmente preferir esse estilo.

Atualmente, não há uma imposição de estilo rigorosa para Nextflow, então o language server oferece alguma flexibilidade. No entanto, ele aplicará consistentemente regras de formatação em torno de definições de métodos e funções para manter a clareza.

### 6.2. Recursos de Organização de Código

#### Comentar Rapidamente

Selecione um bloco de código no seu fluxo de trabalho e pressione **Ctrl+/** (ou **Cmd+/**) para comentá-lo:

```groovy
// workflow {
//     ch_input = channel.fromPath(params.input)
//         .splitCsv(header: true)
//         .map { row -> [row.sample_id, file(row.fastq_path)] }
//
//     FASTQC(ch_input)
// }
```

Isso é perfeito para:

- Desabilitar temporariamente partes de fluxos de trabalho durante o desenvolvimento
- Adicionar comentários explicativos a operações de canal complexas
- Documentar seções de fluxo de trabalho

Use **Ctrl+/** (ou **Cmd+/**) novamente para descomentar o código.

#### Dobramento de Código para Visão Geral

Em `complex_workflow.nf`, note as pequenas setas ao lado das definições de processo. Clique nelas para dobrar (colapsar) processos:

![Dobramento de código](img/code_folding.png)

Isso fornece uma visão geral de alto nível da estrutura do seu fluxo de trabalho sem se perder em detalhes de implementação.

#### Correspondência de Colchetes

Coloque seu cursor próximo a qualquer colchete `{` ou `}` e o VS Code destaca o colchete correspondente. Use **Ctrl+Shift+\\** (ou **Cmd+Shift+\\**) para pular entre colchetes correspondentes.

Isso é crucial para:

- Entender limites de processo
- Encontrar colchetes faltando ou extras
- Navegar estruturas de fluxo de trabalho aninhadas

#### Seleção e Edição Multi-linha

Para editar múltiplas linhas simultaneamente, o VS Code oferece recursos poderosos de múltiplos cursores:

- **Seleção multi-linha**: Segure **Ctrl+Alt** (ou **Cmd+Option** para MacOS) e use as setas do teclado para selecionar múltiplas linhas
- **Indentação multi-linha**: Selecione múltiplas linhas e use **Tab** para indentar ou **Shift+Tab** para remover indentação de blocos inteiros

Isso é particularmente útil para:

- Indentar blocos de processo inteiros consistentemente
- Adicionar comentários a múltiplas linhas de uma vez
- Editar definições de parâmetros similares em múltiplos processos

### Conclusão

Você pode manter código limpo e legível usando formatação automática, recursos de comentário, dobramento de código, correspondência de colchetes e edição multi-linha para organizar fluxos de trabalho complexos eficientemente.

### Qual é o próximo passo?

Aprenda como o VS Code se integra ao seu fluxo de trabalho de desenvolvimento mais amplo além de apenas editar código.

---

## 7. Integração com Fluxo de Trabalho de Desenvolvimento

O VS Code se integra bem ao seu fluxo de trabalho de desenvolvimento além de apenas editar código.

### 7.1. Integração com Controle de Versão

!!! note "Codespaces e Integração Git"

    Se você estiver trabalhando no **GitHub Codespaces**, alguns recursos de integração Git podem não funcionar como esperado, particularmente atalhos de teclado para Controle de Código. Você também pode ter recusado abrir o diretório como um repositório Git durante a configuração inicial, o que é adequado para fins de treinamento.

Se seu projeto é um repositório git (como este é), o VS Code mostra:

- Arquivos modificados com indicadores coloridos
- Status Git na barra de status
- Visualizações de diff inline
- Capacidades de commit e push

Abra o painel de Controle de Código usando o botão de controle de código (![Ícone de controle de código](img/source_control_icon.png)) (`Ctrl+Shift+G` ou `Cmd+Shift+G` se você estiver trabalhando com VSCode localmente) para ver alterações git e fazer stage de commits diretamente no editor.

![Painel de Controle de Código](img/source_control.png)

### 7.2. Executando e Inspecionando Fluxos de Trabalho

Vamos executar um fluxo de trabalho e então inspecionar os resultados. No terminal integrado (`Ctrl+Shift+` crase tanto no Windows quanto no MacOS), execute o fluxo de trabalho básico:

```bash title="Execute o fluxo de trabalho básico"
nextflow run basic_workflow.nf --input data/sample_data.csv --output_dir results
```

Enquanto o fluxo de trabalho executa, você verá saída em tempo real no terminal. Após a conclusão, você pode usar o VS Code para inspecionar resultados sem sair do seu editor:

1. **Navegar para diretórios de trabalho**: Use o explorador de arquivos ou terminal para navegar em `.nextflow/work`
2. **Abrir arquivos de log**: Clique em caminhos de arquivo de log na saída do terminal para abri-los diretamente no VS Code
3. **Inspecionar saídas**: Navegue pelos diretórios de resultados publicados no explorador de arquivos
4. **Visualizar relatórios de execução**: Abra relatórios HTML diretamente no VS Code ou no seu navegador

Isso mantém tudo em um só lugar em vez de alternar entre múltiplas aplicações.

### Conclusão

Você pode integrar o VS Code com controle de versão e execução de fluxo de trabalho para gerenciar todo o seu processo de desenvolvimento a partir de uma única interface.

### Qual é o próximo passo?

Veja como todos esses recursos do IDE trabalham juntos no seu fluxo de trabalho de desenvolvimento diário.

---

## 8. Recapitulação e notas rápidas

Aqui estão algumas notas rápidas sobre cada um dos recursos do IDE discutidos acima:

### 8.1. Iniciando um Novo Recurso

1. **Abertura rápida de arquivo** (`Ctrl+P` ou `Cmd+P`) para encontrar módulos existentes relevantes
2. **Editor dividido** para visualizar processos similares lado a lado
3. **Navegação de símbolos** (`Ctrl+Shift+O` ou `Cmd+Shift+O`) para entender a estrutura do arquivo
4. **Autocompletar** para escrever novo código rapidamente

### 8.2. Depurando Problemas

1. **Painel de problemas** (`Ctrl+Shift+M` ou `Cmd+Shift+M`) para ver todos os erros de uma vez
2. **Ir para definição** (`Ctrl-click` ou `Cmd-click`) para entender interfaces de processo
3. **Encontrar todas as referências** para ver como processos são usados
4. **Pesquisa em todo o projeto** para encontrar padrões ou problemas similares

### 8.3. Refatoração e Melhoria

1. **Pesquisa em todo o projeto** (`Ctrl+Shift+F` ou `Cmd+Shift+F`) para encontrar padrões
2. **Formatação automática** (`Shift+Alt+F` ou `Shift+Option+F`) para manter consistência
3. **Dobramento de código** para focar na estrutura
4. **Integração Git** para rastrear mudanças

---

## Resumo

Você agora teve um tour rápido dos recursos do IDE do VS Code para desenvolvimento Nextflow. Essas ferramentas tornarão você significativamente mais produtivo ao:

- **Reduzir erros** através de verificação de sintaxe em tempo real
- **Acelerar o desenvolvimento** com autocompletar inteligente
- **Melhorar a navegação** em fluxos de trabalho complexos de múltiplos arquivos
- **Manter a qualidade** através de formatação consistente
- **Melhorar a compreensão** através de destaque avançado e visualização de estrutura

Não esperamos que você lembre de tudo, mas agora você sabe que esses recursos existem e poderá encontrá-los quando precisar. À medida que continua desenvolvendo fluxos de trabalho Nextflow, esses recursos do IDE se tornarão uma segunda natureza, permitindo que você se concentre em escrever código de alta qualidade em vez de lutar com sintaxe e estrutura.

### Qual é o próximo passo?

Aplique essas habilidades de IDE enquanto trabalha em outros módulos de treinamento, por exemplo:

- **[nf-test](nf-test.md)**: Crie suítes de teste abrangentes para seus fluxos de trabalho
- **[Hello nf-core](../../hello_nf-core/)**: Construa pipelines de qualidade de produção com padrões da comunidade

O verdadeiro poder desses recursos do IDE emerge à medida que você trabalha em projetos maiores e mais complexos. Comece a incorporá-los ao seu fluxo de trabalho gradualmente—dentro de algumas sessões, eles se tornarão uma segunda natureza e transformarão como você aborda o desenvolvimento Nextflow.

Desde capturar erros antes que eles te desacelerem até navegar bases de código complexas com facilidade, essas ferramentas tornarão você um desenvolvedor mais confiante e eficiente.

Bom código!
