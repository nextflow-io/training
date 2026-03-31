---

# Ambiente de Desenvolvimento

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Tradução assistida por IA - [saiba mais e sugira melhorias](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

Ambientes de Desenvolvimento Integrados (IDEs) modernos podem transformar radicalmente sua experiência de desenvolvimento com Nextflow. Esta missão secundária foca especificamente em aproveitar o VS Code e sua extensão para Nextflow para escrever código mais rápido, detectar erros cedo e navegar por fluxos de trabalho complexos com eficiência.

!!! note "Isso não é um tutorial tradicional"

    Ao contrário de outros módulos de treinamento, este guia é organizado como uma coleção de dicas rápidas, sugestões e exemplos práticos, em vez de um tutorial passo a passo. Cada seção pode ser explorada de forma independente, de acordo com seus interesses e necessidades de desenvolvimento atuais. Sinta-se à vontade para pular entre as seções e focar nos recursos que serão mais imediatamente úteis para o desenvolvimento do seu fluxo de trabalho.

## O que você deve saber antes

Este guia assume que você concluiu o curso de treinamento [Hello Nextflow](../hello_nextflow/) e está confortável com os conceitos fundamentais do Nextflow, incluindo:

- **Estrutura básica do fluxo de trabalho**: Entender processos, fluxos de trabalho e como eles se conectam
- **Operações com canais**: Criar canais, passar dados entre processos e usar operadores básicos
- **Módulos e organização**: Criar módulos reutilizáveis e usar declarações include
- **Noções básicas de configuração**: Usar `nextflow.config` para parâmetros, diretivas de processo e perfis

## O que você aprenderá aqui

Este guia foca em **recursos de produtividade do IDE** que farão de você um desenvolvedor Nextflow mais eficiente:

- **Realce de sintaxe avançado**: Entender o que o VS Code está mostrando sobre a estrutura do seu código
- **Auto-completar inteligente**: Aproveitar sugestões contextuais para escrever código mais rápido
- **Detecção de erros e diagnósticos**: Identificar erros de sintaxe antes de executar seu fluxo de trabalho
- **Navegação no código**: Mover-se rapidamente entre processos, módulos e definições
- **Formatação e organização**: Manter um estilo de código consistente e legível
- **Desenvolvimento assistido por IA** (opcional): Usar ferramentas modernas de IA integradas ao seu IDE

!!! info "Por que recursos de IDE agora?"

    Você provavelmente já estava usando o VS Code durante o curso [Hello Nextflow](../hello_nextflow/), mas mantivemos o foco no aprendizado dos fundamentos do Nextflow em vez dos recursos do IDE. Agora que você está confortável com os conceitos básicos do Nextflow — como processos, fluxos de trabalho, canais e módulos — você está pronto para aproveitar os sofisticados recursos do IDE que farão de você um desenvolvedor mais eficiente.

    Pense nisso como "subir de nível" no seu ambiente de desenvolvimento — o mesmo editor que você tem usado possui capacidades muito mais poderosas que se tornam verdadeiramente valiosas quando você entende o que elas estão te ajudando a fazer.

---

## 0. Configuração e Aquecimento

Vamos configurar um espaço de trabalho especificamente para explorar os recursos do IDE:

```bash title="Navigate to the IDE features directory"
cd side-quests/ide_features
```

Abra este diretório no VS Code:

```bash title="Open VS Code in current directory"
code .
```

O diretório `ide_features` contém fluxos de trabalho de exemplo que demonstram vários recursos do IDE:

```bash title="Show directory structure"
tree .
```

```console title="Project structure"
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
    - `complex_workflow.nf` foi criado apenas para ilustração, para demonstrar recursos de navegação — ele pode não ser executado com sucesso, mas mostra uma estrutura realista de fluxo de trabalho com múltiplos arquivos

### Atalhos de Teclado

Alguns dos recursos neste guia usam atalhos de teclado opcionais. Se você estiver acessando este material via GitHub Codespaces no navegador, alguns atalhos podem não funcionar como esperado, pois são usados para outras funções no seu sistema.

Se você estiver executando o VS Code localmente, como provavelmente fará quando estiver escrevendo fluxos de trabalho de verdade, os atalhos funcionarão conforme descrito.

Se você estiver usando um Mac, alguns (não todos) atalhos de teclado usarão "cmd" em vez de "ctrl", e indicaremos isso no texto como `Ctrl/Cmd`.

### 0.1. Instalando a Extensão Nextflow

!!! note "Já Usando Devcontainers?"

    Se você estiver trabalhando no **GitHub Codespaces** ou usando um **devcontainer local**, a extensão Nextflow provavelmente já está instalada e configurada para você. Você pode pular as etapas de instalação manual abaixo e ir diretamente para explorar os recursos da extensão.

Para instalar a extensão manualmente:

1. Abra o VS Code
2. Vá para a visualização de Extensões clicando no ícone de extensões à esquerda: ![ícone de extensões](img/extensions_icon.png) (atalho `Ctrl/Cmd+Shift+X` se você estiver executando o VSCode localmente)
3. Pesquise por "Nextflow"
4. Instale a extensão oficial do Nextflow

![Instalar a Extensão Nextflow](img/install_extension.png)

### 0.2. Layout do Espaço de Trabalho

Como você já usou o VS Code durante o Hello Nextflow, você já está familiarizado com o básico. Veja como organizar seu espaço de trabalho de forma eficiente para esta sessão:

- **Área do Editor**: Para visualizar e editar arquivos. Você pode dividir em múltiplos painéis para comparar arquivos lado a lado.
- **Explorador de Arquivos** clique em (![ícone do explorador de arquivos](img/files_icon.png)) (`Ctrl/Cmd+Shift+E`): Os arquivos e pastas locais no seu sistema. Mantenha aberto à esquerda para navegar entre arquivos
- **Terminal Integrado** (`Ctrl+Shift+` acento grave para Windows e MacOS): Um terminal para interagir com o computador na parte inferior. Use para executar Nextflow ou outros comandos.
- **Painel de Problemas** (`Ctrl+Shift+M`): O VS Code mostrará aqui todos os erros e problemas detectados. Útil para identificar problemas rapidamente.

Você pode arrastar painéis ou ocultá-los (`Ctrl/Cmd+B` para alternar a barra lateral) para personalizar seu layout enquanto trabalhamos pelos exemplos.

### Conclusão

Você configurou o VS Code com a extensão Nextflow e entende o layout do espaço de trabalho para um desenvolvimento eficiente.

### O que vem a seguir?

Aprenda como o realce de sintaxe ajuda você a entender a estrutura do código Nextflow de relance.

---

## 1. Realce de Sintaxe e Estrutura do Código

Agora que seu espaço de trabalho está configurado, vamos explorar como o realce de sintaxe do VS Code ajuda você a ler e escrever código Nextflow com mais eficiência.

### 1.1. Elementos de Sintaxe do Nextflow

Abra `basic_workflow.nf` para ver o realce de sintaxe em ação:

![Demonstração de Sintaxe](img/syntax_showcase.png)

Observe como o VS Code destaca:

- **Palavras-chave** (`process`, `workflow`, `input`, `output`, `script`) em cores distintas
- **Literais de string** e **parâmetros** com estilos diferentes
- **Comentários** em uma cor discreta
- **Variáveis** e **chamadas de função** com ênfase apropriada
- **Blocos de código** com guias de indentação adequadas

!!! note "Cores Dependentes do Tema"

    As cores específicas que você verá dependerão do tema do VS Code (modo escuro/claro), configurações de cores e quaisquer personalizações que você tenha feito. O importante é que diferentes elementos de sintaxe sejam visualmente distinguidos uns dos outros, tornando a estrutura do código mais fácil de entender independentemente do esquema de cores escolhido.

### 1.2. Entendendo a Estrutura do Código

O realce de sintaxe ajuda você a identificar rapidamente:

- **Limites de processo**: Distinção clara entre diferentes processos
- **Blocos de entrada/saída**: Fácil de identificar definições de fluxo de dados
- **Blocos de script**: Os comandos reais sendo executados
- **Operações com canais**: Etapas de transformação de dados
- **Diretivas de configuração**: Configurações específicas de processo

Essa organização visual se torna inestimável ao trabalhar com fluxos de trabalho complexos contendo múltiplos processos e fluxos de dados intrincados.

### Conclusão

Você entende como o realce de sintaxe do VS Code ajuda a ler a estrutura do código Nextflow e identificar diferentes elementos da linguagem para um desenvolvimento mais rápido.

### O que vem a seguir?

Aprenda como o auto-completar inteligente acelera a escrita de código com sugestões contextuais.

---

## 2. Auto-completar Inteligente

Os recursos de auto-completar do VS Code ajudam você a escrever código mais rápido e com menos erros, sugerindo opções apropriadas com base no contexto.

### 2.1. Sugestões Contextuais

As opções de auto-completar variam dependendo de onde você está no seu código:

#### Operações com Canais

Abra `basic_workflow.nf` novamente e tente digitar `channel.` no bloco workflow:

![Auto-completar de canal](img/autocomplete_channel.png)

Você verá sugestões para:

- `fromPath()` - Criar canal a partir de caminhos de arquivo
- `fromFilePairs()` - Criar canal a partir de arquivos pareados
- `of()` - Criar canal a partir de valores
- `fromSRA()` - Criar canal a partir de acessos SRA
- E muito mais...

Isso ajuda você a encontrar rapidamente o factory de canal correto sem precisar lembrar os nomes exatos dos métodos.

Você também pode descobrir os operadores disponíveis para aplicar aos canais. Por exemplo, digite `FASTQC.out.html.` para ver as operações disponíveis:

![Auto-completar de operações de canal](img/autocomplete_operators.png)

#### Diretivas de Processo

Dentro de um bloco script de processo, digite `task.` para ver as propriedades de runtime disponíveis:

![Auto-completar de propriedades de task](img/autocomplete_task.png)

#### Configuração

Abra nextflow.config e digite `process.` em qualquer lugar para ver as diretivas de processo disponíveis:

![Auto-completar de configuração](img/autocomplete_config.png)

Você verá sugestões para:

- `executor`
- `memory`
- `cpus`

Isso economiza tempo ao configurar processos e funciona em diferentes escopos de configuração. Por exemplo, tente digitar `docker.` para ver as opções de configuração específicas do Docker.

### Conclusão

Você pode usar o auto-completar inteligente do VS Code para descobrir operações de canal disponíveis, diretivas de processo e opções de configuração sem memorizar a sintaxe.

### O que vem a seguir?

Aprenda como a detecção de erros em tempo real ajuda você a identificar problemas antes de executar seu fluxo de trabalho, simplesmente lendo o código.

## 3. Detecção de Erros e Diagnósticos

A detecção de erros em tempo real do VS Code ajuda você a identificar problemas antes de executar seu fluxo de trabalho.

### 3.1. Detecção de Erros de Sintaxe

Vamos criar um erro deliberado para ver a detecção em ação. Abra `basic_workflow.nf` e altere o nome do processo de `FASTQC` para `FASTQ` (ou qualquer outro nome inválido). O VS Code imediatamente destacará o erro no bloco workflow com um sublinhado vermelho ondulado:

![Sublinhado de erro](img/error_underline.png)

### 3.2. Painel de Problemas

Além do destaque individual de erros, o VS Code fornece um Painel de Problemas centralizado que agrega todos os erros, avisos e mensagens informativas em seu espaço de trabalho. Abra-o com `Ctrl/Cmd+Shift+M` e use o ícone de filtro para mostrar apenas os erros relevantes ao arquivo atual:

![Filtrar o painel de problemas](img/active_file.png)

Clique em qualquer problema para ir diretamente à linha problemática

![Painel de Problemas](img/problems_panel.png)

Corrija o erro alterando o nome do processo de volta para `FASTQC`.

### 3.3. Padrões Comuns de Erro

Erros comuns na sintaxe do Nextflow incluem:

- **Colchetes faltando**: `{` ou `}` sem correspondência
- **Blocos incompletos**: Seções obrigatórias ausentes em processos
- **Sintaxe inválida**: DSL Nextflow malformado
- **Erros de digitação em palavras-chave**: Diretivas de processo com erros ortográficos
- **Incompatibilidades de canal**: Incompatibilidades de tipo

O servidor de linguagem Nextflow destaca esses problemas no Painel de Problemas. Você pode verificá-los antecipadamente para evitar erros de sintaxe ao executar um pipeline.

### Conclusão

Você pode usar a detecção de erros e o Painel de Problemas do VS Code para identificar erros de sintaxe e problemas antes de executar seu fluxo de trabalho, economizando tempo e evitando frustrações.

### O que vem a seguir?

Aprenda como navegar eficientemente entre processos, módulos e definições em fluxos de trabalho complexos.

---

## 4. Navegação no Código e Gerenciamento de Símbolos

A navegação eficiente é crucial ao trabalhar com fluxos de trabalho complexos que abrangem múltiplos arquivos. Para entender isso, substitua a definição de processo em `basic_workflow.nf` por uma importação do módulo que fornecemos:

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

Se você passar o mouse sobre um nome de processo como `FASTQC`, verá um popup com a interface do módulo (entradas e saídas):

![Ir para definição](img/syntax.png)

Este recurso é particularmente valioso ao criar fluxos de trabalho, pois permite entender a interface do módulo sem abrir o arquivo do módulo diretamente.

Você pode navegar rapidamente para qualquer definição de processo, módulo ou variável usando **Ctrl/Cmd-clique**. Passe o mouse sobre o link para o arquivo do módulo no topo do script e siga o link conforme sugerido:

![Seguir link](img/follow_link.png)

O mesmo funciona para nomes de processo. Volte para `basic_workflow.nf` e tente isso no nome do processo `FASTQC` no bloco workflow. Isso leva você diretamente ao nome do processo (que é o mesmo que o arquivo do módulo neste exemplo, mas poderia estar no meio de um arquivo muito maior).

Para voltar ao ponto anterior, use **Alt+←** (ou **Ctrl+-** no Mac). Esta é uma forma poderosa de explorar o código sem perder seu lugar.

Agora vamos explorar a navegação em um fluxo de trabalho mais complexo usando `complex_workflow.nf` (o arquivo apenas para ilustração mencionado anteriormente). Este fluxo de trabalho contém múltiplos processos definidos em arquivos de módulo separados, além de alguns inline. Embora estruturas complexas com múltiplos arquivos possam ser desafiadoras de navegar manualmente, a capacidade de pular para definições torna a exploração muito mais gerenciável.

1. Abra `complex_workflow.nf`
2. Navegue para as definições de módulo
3. Use **Alt+←** (ou **Ctrl+-**) para navegar de volta
4. Navegue para o nome do processo `FASTQC` no bloco workflow. Isso leva você diretamente ao nome do processo (que é o mesmo que o arquivo do módulo neste exemplo, mas poderia estar no meio de um arquivo muito maior).
5. Navegue de volta novamente
6. Navegue para o processo `TRIM_GALORE` no bloco workflow. Este está definido inline, então não levará você a um arquivo separado, mas ainda mostrará a definição do processo, e você ainda pode navegar de volta ao ponto anterior.

### 4.2. Navegação por Símbolos

Com `complex_workflow.nf` ainda aberto, você pode obter uma visão geral de todos os símbolos no arquivo digitando `@` na barra de pesquisa no topo do VSCode (o atalho de teclado é `Ctrl/Cmd+Shift+O`, mas pode não funcionar no Codespaces). Isso abre o painel de navegação por símbolos, que lista todos os símbolos no arquivo atual:

![Navegação por símbolos](img/symbols.png)

Isso mostra:

- Todas as definições de processo
- Definições de workflow (há dois workflows definidos neste arquivo)
- Definições de função

Comece a digitar para filtrar os resultados.

### 4.3. Encontrar Todas as Referências

Entender onde um processo ou variável é usado em toda a sua base de código pode ser muito útil. Por exemplo, se você quiser encontrar todas as referências ao processo `FASTQC`, comece navegando para sua definição. Você pode fazer isso abrindo `modules/fastqc.nf` diretamente, ou usando o recurso de navegação rápida do VS Code com `Ctrl/Cmd-clique` como fizemos acima. Uma vez na definição do processo, clique com o botão direito no nome do processo `FASTQC` e selecione "Find All References" no menu de contexto para ver todas as instâncias onde ele é usado.

![Encontrar referências](img/references.png)

Este recurso exibe todas as instâncias onde `FASTQC` é referenciado em seu espaço de trabalho, incluindo seu uso nos dois fluxos de trabalho distintos. Essa informação é crucial para avaliar o impacto potencial de modificações no processo `FASTQC`.

### 4.4. Painel de Estrutura

O painel de Estrutura (Outline), localizado na barra lateral do Explorador (clique em ![ícone do Explorador](img/files_icon.png)), fornece uma visão geral conveniente de todos os símbolos no arquivo atual. Este recurso permite navegar e gerenciar rapidamente a estrutura do seu código, exibindo funções, variáveis e outros elementos-chave em uma visualização hierárquica.

![Painel de Estrutura](img/outline.png)

Use o painel de Estrutura para navegar rapidamente para diferentes partes do seu código sem usar o explorador de arquivos.

### 4.5. Visualização do DAG

A extensão Nextflow do VS Code pode visualizar seu fluxo de trabalho como um Grafo Acíclico Dirigido (DAG). Isso ajuda você a entender o fluxo de dados e as dependências entre processos. Abra `complex_workflow.nf` e clique no botão "Preview DAG" acima de `workflow {` (o segundo bloco `workflow` neste arquivo):

![Pré-visualização do DAG](img/dag_preview.png)

Este é apenas o workflow de 'entrada', mas você também pode pré-visualizar o DAG para os workflows internos clicando no botão "Preview DAG" acima do workflow `RNASEQ_PIPELINE {` mais acima:

![Pré-visualização do DAG do workflow interno](img/dag_preview_inner.png)

Para este fluxo de trabalho, você pode usar os nós no DAG para navegar para as definições de processo correspondentes no código. Clique em um nó e ele levará você à definição de processo relevante no editor. Especialmente quando um fluxo de trabalho cresce muito, isso pode realmente ajudar a navegar pelo código e entender como os processos estão conectados.

### Conclusão

Você pode navegar por fluxos de trabalho complexos com eficiência usando ir para definição, pesquisa de símbolos, encontrar referências e visualização do DAG para entender a estrutura do código e as dependências.

### O que vem a seguir?

Aprenda como trabalhar efetivamente com múltiplos arquivos interconectados em projetos Nextflow maiores.

## 5. Trabalhando com Múltiplos Arquivos

O desenvolvimento real com Nextflow envolve trabalhar com múltiplos arquivos interconectados. Vamos explorar como o VS Code ajuda você a gerenciar projetos complexos com eficiência.

### 5.1. Navegação Rápida entre Arquivos

Com `complex_workflow.nf` aberto, você notará que ele importa vários módulos. Vamos praticar a navegação rápida entre eles.

Pressione **Ctrl+P** (ou **Cmd+P**) e comece a digitar "fast":

O VS Code mostrará os arquivos correspondentes. Selecione `modules/fastqc.nf` para ir diretamente para lá. Isso é muito mais rápido do que clicar pelo explorador de arquivos quando você sabe aproximadamente qual arquivo está procurando.

Tente isso com outros padrões:

- Digite "star" para encontrar o arquivo do módulo de alinhamento STAR (`star.nf`)
- Digite "utils" para encontrar o arquivo de funções utilitárias (`utils.nf`)
- Digite "config" para ir para os arquivos de configuração (`nextflow.config`)

### 5.2. Editor Dividido para Desenvolvimento com Múltiplos Arquivos

Ao trabalhar com módulos, muitas vezes você precisa ver tanto o fluxo de trabalho principal quanto as definições de módulo simultaneamente. Vamos configurar isso:

1. Abra `complex_workflow.nf`
2. Abra `modules/fastqc.nf` em uma nova aba
3. Clique com o botão direito na aba `modules/fastqc.nf` e selecione "Split Right"
4. Agora você pode ver ambos os arquivos lado a lado

![Editor dividido](img/split_editor.png)

Isso é inestimável quando:

- Verificando interfaces de módulo ao escrever chamadas de workflow, e a pré-visualização não é suficiente
- Comparando processos semelhantes em diferentes módulos
- Depurando o fluxo de dados entre workflow e módulos

### 5.3. Pesquisa em Todo o Projeto

Às vezes você precisa encontrar onde padrões específicos são usados em todo o seu projeto. Pressione `Ctrl/Cmd+Shift+F` para abrir o painel de pesquisa.

Tente pesquisar por `publishDir` em todo o espaço de trabalho:

![Pesquisa no projeto](img/project_search.png)

Isso mostra todos os arquivos que usam diretórios de publicação, ajudando você a:

- Entender padrões de organização de saída
- Encontrar exemplos de diretivas específicas
- Garantir consistência entre módulos

### Conclusão

Você pode gerenciar projetos complexos com múltiplos arquivos usando navegação rápida entre arquivos, editores divididos e pesquisa em todo o projeto para trabalhar eficientemente entre fluxos de trabalho e módulos.

### O que vem a seguir?

Aprenda como os recursos de formatação e manutenção de código mantêm seus fluxos de trabalho organizados e legíveis.

---

## 6. Formatação e Manutenção do Código

A formatação adequada do código é essencial não apenas para a estética, mas também para melhorar a legibilidade, a compreensão e a facilidade de atualização de fluxos de trabalho complexos.

### 6.1. Formatação Automática em Ação

Abra `basic_workflow.nf` e deliberadamente bagunce a formatação:

- Remova alguma indentação: Selecione o documento inteiro e pressione `shift+tab` várias vezes para remover o máximo de indentações possível.
- Adicione espaços extras em lugares aleatórios: na instrução `channel.fromPath`, adicione 30 espaços após o `(`.
- Quebre algumas linhas de forma estranha: Adicione uma nova linha entre o operador `.view {` e a string `Processing sample:`, mas não adicione uma nova linha correspondente antes do parêntese de fechamento `}`.

Agora pressione `Shift+Alt+F` (ou `Shift+Option+F` no MacOS) para formatar automaticamente:

O VS Code imediatamente:

- Corrige a indentação para mostrar a estrutura do processo claramente
- Alinha elementos semelhantes de forma consistente
- Remove espaços em branco desnecessários
- Mantém quebras de linha legíveis

Observe que a formatação automática pode não resolver todos os problemas de estilo de código. O servidor de linguagem Nextflow visa manter seu código organizado, mas também respeita suas preferências pessoais em certas áreas. Por exemplo, se você remover a indentação dentro do bloco `script` de um processo, o formatador deixará como está, pois você pode intencionalmente preferir esse estilo.

Atualmente, não há aplicação estrita de estilo para Nextflow, então o servidor de linguagem oferece alguma flexibilidade. No entanto, ele aplicará consistentemente regras de formatação em torno de definições de métodos e funções para manter a clareza.

### 6.2. Recursos de Organização do Código

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
- Documentar seções do fluxo de trabalho

Use **Ctrl+/** (ou **Cmd+/**) novamente para descomentar o código.

#### Dobramento de Código para Visão Geral

Em `complex_workflow.nf`, observe as pequenas setas ao lado das definições de processo. Clique nelas para dobrar (recolher) os processos:

![Dobramento de código](img/code_folding.png)

Isso fornece uma visão geral de alto nível da estrutura do seu fluxo de trabalho sem se perder nos detalhes de implementação.

#### Correspondência de Colchetes

Posicione o cursor ao lado de qualquer colchete `{` ou `}` e o VS Code destaca o colchete correspondente. Use **Ctrl+Shift+\\** (ou **Cmd+Shift+\\**) para pular entre colchetes correspondentes.

Isso é crucial para:

- Entender os limites do processo
- Encontrar colchetes faltando ou extras
- Navegar por estruturas de workflow aninhadas

#### Seleção e Edição de Múltiplas Linhas

Para editar múltiplas linhas simultaneamente, o VS Code oferece poderosas capacidades de múltiplos cursores:

- **Seleção de múltiplas linhas**: Segure **Ctrl+Alt** (ou **Cmd+Option** no MacOS) e use as teclas de seta para selecionar múltiplas linhas
- **Indentação de múltiplas linhas**: Selecione múltiplas linhas e use **Tab** para indentar ou **Shift+Tab** para desindentar blocos inteiros

Isso é particularmente útil para:

- Indentar blocos de processo inteiros de forma consistente
- Adicionar comentários a múltiplas linhas de uma vez
- Editar definições de parâmetros semelhantes em múltiplos processos

### Conclusão

Você pode manter um código limpo e legível usando formatação automática, recursos de comentário, dobramento de código, correspondência de colchetes e edição de múltiplas linhas para organizar fluxos de trabalho complexos com eficiência.

### O que vem a seguir?

Aprenda como o VS Code se integra ao seu fluxo de trabalho de desenvolvimento mais amplo, além de apenas editar código.

---

## 7. Integração com o Fluxo de Trabalho de Desenvolvimento

O VS Code se integra bem ao seu fluxo de trabalho de desenvolvimento além de apenas editar código.

### 7.1. Integração com Controle de Versão

!!! note "Codespaces e Integração com Git"

    Se você estiver trabalhando no **GitHub Codespaces**, alguns recursos de integração com Git podem não funcionar como esperado, especialmente atalhos de teclado para Controle de Código-Fonte. Você também pode ter recusado abrir o diretório como um repositório Git durante a configuração inicial, o que é adequado para fins de treinamento.

Se o seu projeto for um repositório git (como este é), o VS Code mostra:

- Arquivos modificados com indicadores coloridos
- Status do Git na barra de status
- Visualizações de diff inline
- Capacidades de commit e push

Abra o painel de Controle de Código-Fonte usando o botão de controle de código-fonte (![ícone de controle de código-fonte](img/source_control_icon.png)) (`Ctrl+Shift+G` ou `Cmd+Shift+G` se você estiver trabalhando com VSCode localmente) para ver as alterações do git e fazer commits diretamente no editor.

![Painel de Controle de Código-Fonte](img/source_control.png)

### 7.2. Executando e Inspecionando Fluxos de Trabalho

Vamos executar um fluxo de trabalho e depois inspecionar os resultados. No terminal integrado (`Ctrl+Shift+` acento grave no Windows e MacOS), execute o fluxo de trabalho básico:

```bash title="Run the basic workflow"
nextflow run basic_workflow.nf --input data/sample_data.csv --output_dir results
```

Enquanto o fluxo de trabalho é executado, você verá a saída em tempo real no terminal. Após a conclusão, você pode usar o VS Code para inspecionar os resultados sem sair do editor:

1. **Navegar para os diretórios de trabalho**: Use o explorador de arquivos ou terminal para navegar por `.nextflow/work`
2. **Abrir arquivos de log**: Clique nos caminhos de arquivos de log na saída do terminal para abri-los diretamente no VS Code
3. **Inspecionar saídas**: Navegue pelos diretórios de resultados publicados no explorador de arquivos
4. **Visualizar relatórios de execução**: Abra relatórios HTML diretamente no VS Code ou no seu navegador

Isso mantém tudo em um só lugar, em vez de alternar entre múltiplos aplicativos.

### Conclusão

Você pode integrar o VS Code com controle de versão e execução de fluxo de trabalho para gerenciar todo o seu processo de desenvolvimento a partir de uma única interface.

### O que vem a seguir?

Veja como todos esses recursos do IDE funcionam juntos no seu fluxo de trabalho de desenvolvimento diário.

---

## 8. Recapitulação e Notas Rápidas

Aqui estão algumas notas rápidas sobre cada um dos recursos do IDE discutidos acima:

### 8.1. Iniciando um Novo Recurso

1. **Abertura rápida de arquivo** (`Ctrl+P` ou `Cmd+P`) para encontrar módulos existentes relevantes
2. **Editor dividido** para visualizar processos semelhantes lado a lado
3. **Navegação por símbolos** (`Ctrl+Shift+O` ou `Cmd+Shift+O`) para entender a estrutura do arquivo
4. **Auto-completar** para escrever novo código rapidamente

### 8.2. Depurando Problemas

1. **Painel de Problemas** (`Ctrl+Shift+M` ou `Cmd+Shift+M`) para ver todos os erros de uma vez
2. **Ir para definição** (`Ctrl-clique` ou `Cmd-clique`) para entender interfaces de processo
3. **Encontrar todas as referências** para ver como os processos são usados
4. **Pesquisa em todo o projeto** para encontrar padrões ou problemas semelhantes

### 8.3. Refatoração e Melhoria

1. **Pesquisa em todo o projeto** (`Ctrl+Shift+F` ou `Cmd+Shift+F`) para encontrar padrões
2. **Formatação automática** (`Shift+Alt+F` ou `Shift+Option+F`) para manter consistência
3. **Dobramento de código** para focar na estrutura
4. **Integração com Git** para rastrear alterações

---

## Resumo

Você agora fez um tour rápido pelos recursos do IDE do VS Code para desenvolvimento com Nextflow. Essas ferramentas farão de você um desenvolvedor significativamente mais produtivo ao:

- **Reduzir erros** por meio da verificação de sintaxe em tempo real
- **Acelerar o desenvolvimento** com auto-completar inteligente
- **Melhorar a navegação** em fluxos de trabalho complexos com múltiplos arquivos
- **Manter a qualidade** por meio de formatação consistente
- **Aprimorar a compreensão** por meio de realce avançado e visualização de estrutura

Não esperamos que você se lembre de tudo, mas agora que você sabe que esses recursos existem, poderá encontrá-los quando precisar. À medida que você continua desenvolvendo fluxos de trabalho Nextflow, esses recursos do IDE se tornarão naturais, permitindo que você se concentre em escrever código de alta qualidade em vez de lutar com sintaxe e estrutura.

### O que vem a seguir?

Aplique essas habilidades de IDE enquanto trabalha em outros módulos de treinamento, por exemplo:

- **[nf-test](nf-test.md)**: Crie suítes de teste abrangentes para seus fluxos de trabalho
- **[Hello nf-core](../../hello_nf-core/)**: Construa pipelines de qualidade de produção com padrões da comunidade

O verdadeiro poder desses recursos do IDE emerge à medida que você trabalha em projetos maiores e mais complexos. Comece a incorporá-los ao seu fluxo de trabalho gradualmente — em poucas sessões, eles se tornarão naturais e transformarão a forma como você aborda o desenvolvimento com Nextflow.

Desde detectar erros antes que eles te atrasem até navegar por bases de código complexas com facilidade, essas ferramentas farão de você um desenvolvedor mais confiante e eficiente.

Bom código!
