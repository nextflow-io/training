# Parte 3: Perfil de recursos e otimização

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Tradução assistida por IA - [saiba mais e sugira melhorias](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

ESTE É UM PLACEHOLDER

!!!note "Nota"

    Este módulo de treinamento está sendo refeito.

---

TODO

### 1.1. Execute o fluxo de trabalho para gerar um relatório de utilização de recursos

Para que o Nextflow gere o relatório automaticamente, basta adicionar `-with-report <nome-do-arquivo>.html` à sua linha de comando.

```bash
nextflow run main.nf -profile my_laptop -with-report report-config-1.html
```

O relatório é um arquivo html, que você pode baixar e abrir no seu navegador. Você também pode clicar com o botão direito nele no explorador de arquivos à esquerda e clicar em `Show preview` para visualizá-lo no VS Code.

Reserve alguns minutos para analisar o relatório e ver se você consegue identificar algumas oportunidades para ajustar recursos.
Certifique-se de clicar nas abas que mostram os resultados de utilização como uma porcentagem do que foi alocado.
Há alguma [documentação](https://www.nextflow.io/docs/latest/reports.html) descrevendo todos os recursos disponíveis.

<!-- TODO: insert images -->

Uma observação é que o `GATK_JOINTGENOTYPING` parece estar muito faminto por CPU, o que faz sentido já que ele realiza muitos cálculos complexos.
Então poderíamos tentar aumentar isso e ver se reduz o tempo de execução.

No entanto, parece que exageramos nas alocações de memória; todos os processos estão usando apenas uma fração do que estamos dando a eles.
Devemos reduzir isso e economizar alguns recursos.

### 1.2. Ajuste as alocações de recursos para um processo específico

Podemos especificar alocações de recursos para um determinado processo usando o seletor de processo `withName`.
A sintaxe fica assim quando está sozinha em um bloco de processo:

```groovy title="Sintaxe"
process {
    withName: 'GATK_JOINTGENOTYPING' {
        cpus = 4
    }
}
```

Vamos adicionar isso ao bloco de processo existente no arquivo `nextflow.config`.

```groovy title="nextflow.config" linenums="11"
process {
    // padrões para todos os processos
    cpus = 2
    memory = 2.GB
    // alocações para um processo específico
    withName: 'GATK_JOINTGENOTYPING' {
        cpus = 4
    }
}
```

Com isso especificado, as configurações padrão se aplicarão a todos os processos **exceto** o processo `GATK_JOINTGENOTYPING`, que é um floco de neve especial que recebe muito mais CPU.
Esperamos que isso tenha um efeito.

### 1.3. Execute novamente com a configuração modificada

Vamos executar o fluxo de trabalho novamente com a configuração modificada e com o sinalizador de relatório ativado, mas observe que estamos dando ao relatório um nome diferente para que possamos diferenciá-los.

```bash
nextflow run main.nf -profile my_laptop -with-report report-config-2.html
```

Mais uma vez, você provavelmente não notará uma diferença substancial no tempo de execução, porque esta é uma carga de trabalho muito pequena e as ferramentas gastam mais tempo em tarefas auxiliares do que na execução do trabalho 'real'.

No entanto, o segundo relatório mostra que nossa utilização de recursos está mais equilibrada agora.

<!-- **TODO: screenshots?** -->

Como você pode ver, esta abordagem é útil quando seus processos têm requisitos de recursos diferentes. Ela permite que você dimensione adequadamente as alocações de recursos configuradas para cada processo com base em dados reais, não em suposições.

!!!note "Nota"

    Isto é apenas um pequeno aperitivo do que você pode fazer para otimizar o uso de recursos.
    O próprio Nextflow tem uma [lógica de repetição dinâmica](https://www.nextflow.io/docs/latest/process.html#dynamic-task-resources) realmente interessante integrada para repetir tarefas que falham devido a limitações de recursos.
    Além disso, a Seqera Platform oferece ferramentas orientadas por IA para otimizar suas alocações de recursos automaticamente também.

    Abordaremos ambas as abordagens em uma parte futura deste curso de treinamento.

Dito isso, pode haver algumas restrições sobre o que você pode (ou deve) alocar dependendo de qual executor de computação e infraestrutura de computação você está usando. Por exemplo, seu cluster pode exigir que você permaneça dentro de certos limites que não se aplicam quando você está executando em outro lugar.
