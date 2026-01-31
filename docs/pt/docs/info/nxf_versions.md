---
title: Versões do Nextflow
description: Entendendo e gerenciando a evolução das versões de sintaxe do Nextflow
hide:
  - toc
  - footer
---

## Versão de sintaxe do Nextflow atualmente suportada e requisitos

A partir da versão 3.0 do portal de treinamento, todos os nossos cursos de treinamento são baseados na versão 25.10.2 do Nextflow, a menos que especificado de outra forma na página de índice do curso (exceto materiais descontinuados ou arquivados que podem não incluir uma nota de versão).

Como os cursos agora usam entradas tipadas no nível do fluxo de trabalho, bem como diretivas de saída no nível do fluxo de trabalho, eles requerem o uso do analisador de sintaxe V2.
Se você planeja usar o ambiente que fornecemos através do [Github Codespaces](../envsetup/01_setup.md) ou [devcontainers locais](../envsetup/03_devcontainer.md), você não precisa fazer nada, a menos que seja especificamente indicado nas instruções do curso.
No entanto, se você está planejando trabalhar através dos treinamentos em seu próprio ambiente ([Instalação manual](../envsetup/02_local.md)), você precisará garantir o uso do Nextflow versão 25.10.2 ou posterior com o analisador de sintaxe v2 habilitado.

## Versões mais antigas dos materiais de treinamento

Nossos materiais de treinamento são versionados desde fevereiro de 2025.

Você pode acessar versões mais antigas dos materiais de treinamento que funcionam com versões do Nextflow **anteriores a 25.10.2** através do menu suspenso no topo de cada página que mostra a versão numerada dos materiais de treinamento.
Quando você seleciona uma versão mais antiga dos materiais de treinamento, os links para o ambiente de treinamento especificarão automaticamente a versão correspondente do ambiente.

## Outras informações sobre versões de sintaxe do Nextflow

O Nextflow tem dois conceitos de versionamento distintos que às vezes são confundidos: **versões DSL** e **versões do analisador de sintaxe**.

**DSL1 vs DSL2** refere-se a maneiras fundamentalmente diferentes de escrever pipelines Nextflow.
DSL1 era a sintaxe original onde os processos eram implicitamente conectados através de canais.
DSL2, introduzido no Nextflow 20.07, adicionou recursos de modularidade: a capacidade de importar processos e fluxos de trabalho de outros arquivos, blocos `workflow` explícitos e saídas de processos nomeadas.
DSL1 foi descontinuado no Nextflow 22.03 e removido em 22.12.
Todo código Nextflow moderno usa DSL2.

**Analisador de sintaxe v1 vs v2** refere-se a diferentes analisadores que ambos funcionam com código DSL2.
O analisador v1 é o original, mais permissivo.
O analisador v2 é mais rigoroso e habilita novos recursos de linguagem, como tipagem estática (entradas e saídas tipadas) e diretivas de saída no nível do fluxo de trabalho.
O analisador v2 também fornece melhores mensagens de erro e detecta mais erros no momento da análise, em vez de em tempo de execução.
O analisador v2 se tornará o padrão no Nextflow 26.04.

Em resumo: DSL2 é a linguagem que você escreve; a versão do analisador de sintaxe determina quão rigorosamente essa linguagem é interpretada e quais recursos avançados estão disponíveis.

### Verificando e definindo a versão do Nextflow

Você pode verificar qual versão do Nextflow está instalada em seu sistema usando o comando `nextflow --version`.

Para mais informações sobre como atualizar sua versão do Nextflow, consulte a documentação de referência sobre [Atualizando o Nextflow](https://www.nextflow.io/docs/latest/updating-nextflow.html).

### Habilitando o analisador de sintaxe v2

Para **habilitar** o analisador de sintaxe v2 para sua sessão atual, execute o seguinte comando em seu terminal:

```bash
export NXF_SYNTAX_PARSER=v2
```

Para tornar isso permanente (enquanto v2 não se torna o padrão no Nextflow 26.04), adicione o comando export ao perfil do seu shell (`~/.bashrc`, `~/.zshrc`, etc.):

```bash
echo 'export NXF_SYNTAX_PARSER=v2' >> ~/.bashrc
source ~/.bashrc
```

Note que a variável de ambiente `NXF_SYNTAX_PARSER=v2` é um requisito temporário.
A partir do Nextflow 26.04, o analisador v2 se tornará o padrão e essa configuração não será mais necessária.

### Desabilitando o analisador de sintaxe v2

Para **desabilitar** o analisador de sintaxe v2 para sua sessão atual, execute o seguinte comando em seu terminal:

```bash
export NXF_SYNTAX_PARSER=v1
```

<!-- Will it be possible to disable it in versions after 26.04? -->

### Migrando código existente

Para orientações sobre a migração de código existente para estar em conformidade com versões mais recentes do Nextflow, consulte as [Notas de Migração](https://www.nextflow.io/docs/latest/migrations/index.html) na documentação de referência.

Estes dois artigos são particularmente úteis para migrar para a versão mais recente:

- [Migrando para saídas de fluxo de trabalho](https://www.nextflow.io/docs/latest/tutorials/workflow-outputs.html)
- [Migrando para tipos estáticos](https://www.nextflow.io/docs/latest/tutorials/static-types.html)

Ambos os recursos são abordados como parte do treinamento para iniciantes a partir da versão 3.0 dos materiais de treinamento.

Dependendo da geração do código Nextflow que você pretende migrar, você pode conseguir fazer a maior parte do trabalho pelo linter do Nextflow usando o comando `nextflow lint -format`.
Consulte a referência CLI para [`lint`](https://www.nextflow.io/docs/latest/reference/cli.html#lint) para mais detalhes.

Esperamos que isso seja útil.
Se você precisar de ajuda, entre em contato no Slack ou no fórum.

---

<div markdown class="homepage_logos">

![Seqera](../assets/img/seqera_logo.png#only-light)

![Seqera](../assets/img/seqera_logo_dark.png#only-dark)

</div>
