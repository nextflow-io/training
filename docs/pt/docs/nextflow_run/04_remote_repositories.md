# Parte 4: Executar pipelines remotos

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Tradução assistida por IA - [saiba mais e sugira melhorias](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

Até agora, você executou scripts de fluxo de trabalho armazenados localmente.
Na prática, muitas vezes você vai querer executar pipelines publicados em repositórios remotos, como o GitHub, sem precisar baixá-los você mesmo.

O Nextflow torna isso simples: você pode executar qualquer pipeline diretamente a partir de uma URL de repositório Git.

---

## 1. Executar um pipeline do GitHub

A sintaxe básica para executar um pipeline remoto é `nextflow run <repositório>`, onde `<repositório>` pode ser um caminho de repositório GitHub como `nextflow-io/hello`, uma URL completa, ou um caminho para GitLab, Bitbucket ou outro serviço de hospedagem Git.

### 1.1. Iniciar o pipeline

Execute o pipeline de demonstração oficial "hello" do Nextflow.
Este é um pipeline diferente e muito mais simples do que o que você vem executando neste curso: ele é anterior ao pipeline "Hello" usado ao longo deste treinamento, e apenas imprime uma saudação para cada uma de algumas linguagens predefinidas, portanto não espere a entrada CSV ou a arte ASCII a que você está acostumado.

```bash
nextflow run nextflow-io/hello
```

??? success "Saída do comando"

    ```console
    N E X T F L O W   ~  version 26.04.4

    Pulling nextflow-io/hello ...
     downloaded from https://github.com/nextflow-io/hello.git
    Launching `https://github.com/nextflow-io/hello` [sleepy_swanson] revision: 3c2cdc9823 [master]

    executor >  local (4)
    [ba/08236d] sayHello (4) | 4 of 4 ✔
    Ciao world!

    Hello world!

    Bonjour world!

    Hola world!
    ```

### 1.2. Encontrar onde o pipeline está em cache

Na primeira vez que você executa um pipeline remoto, o Nextflow o baixa e armazena em cache localmente.
Execuções subsequentes reutilizam a versão em cache, a menos que você solicite explicitamente uma atualização.

Por padrão, o Nextflow salva os pipelines baixados em `$NXF_HOME/assets`.
Para descobrir onde um pipeline específico foi salvo e quais revisões estão disponíveis, pergunte diretamente ao Nextflow:

```bash
nextflow info nextflow-io/hello
```

??? success "Saída do comando"

    ```console
     project name: nextflow-io/hello
     repository  : https://github.com/nextflow-io/hello
     local path  : /workspaces/.nextflow/assets/.repos/nextflow-io/hello
     main script : main.nf
     revisions   :
     > master (default)
       mybranch
       testing
       v1.1 [t]
       v1.2 [t]
       v1.3 [t]
    ```

    O Nextflow marca com `>` cada revisão que você já obteve localmente; as demais estão disponíveis, mas ainda não foram baixadas para uma cópia de trabalho.

Você também pode listar todos os pipelines que já baixou com `nextflow list`:

```bash
nextflow list
```

??? success "Saída do comando"

    ```console
    nextflow-io/hello
    ```

O curso [Use nf-core](../nfcore_use/01_run_demo.md#12-retrieve-the-pipeline-code) aborda esse mecanismo de cache com mais profundidade, incluindo como navegar pelo código-fonte de um pipeline baixado.

### Conclusão

Você sabe como executar um pipeline diretamente de um repositório GitHub sem baixá-lo você mesmo, e onde encontrá-lo localmente depois.

### O que vem a seguir?

Aprenda como fixar uma versão específica de um pipeline remoto para garantir reprodutibilidade.

---

## 2. Especificar uma versão para reprodutibilidade

Por padrão, o Nextflow executa a revisão mais recente do branch padrão.
Você pode fixar uma versão específica (tag), branch ou commit usando a flag `-r`.

### 2.1. Fixar uma revisão específica

```bash
nextflow run nextflow-io/hello -r v1.3
```

??? success "Saída do comando"

    ```console
    N E X T F L O W   ~  version 26.04.4

    Pulling nextflow-io/hello:v1.3 ...
     downloaded from https://github.com/nextflow-io/hello.git
    Launching `https://github.com/nextflow-io/hello` [sick_carson] revision: 2ce0b0e294 [v1.3]

    executor >  local (4)
    [61/e11f77] sayHello (4) | 4 of 4 ✔
    Ciao world!

    Bonjour world!

    Hello world!

    Hola world!
    ```

O Nextflow baixa essa revisão na primeira vez que você a solicita, daí as linhas `Pulling` e `downloaded from`; solicitar a mesma revisão novamente mais tarde vai direto para `Launching`.
Fixar uma revisão exata é essencial para a reprodutibilidade.
Isso garante que você e seus colaboradores executem exatamente o mesmo código de pipeline, independentemente do que tenha mudado no repositório desde então.

### 2.2. As revisões se aplicam apenas por invocação

Fixar uma revisão com `-r` afeta apenas a execução em que você a especifica: não altera o que uma execução posterior com `nextflow run` simples utilizará.
Tente executar o pipeline novamente sem `-r`:

```bash
nextflow run nextflow-io/hello
```

??? success "Saída do comando"

    ```console
    N E X T F L O W   ~  version 26.04.4

    Launching `https://github.com/nextflow-io/hello` [hungry_maxwell] revision: 3c2cdc9823 [master]

    executor >  local (4)
    [ba/08236d] sayHello (4) | 4 of 4 ✔
    Ciao world!

    Hello world!

    Bonjour world!

    Hola world!
    ```

Mesmo que a execução anterior tenha fixado explicitamente `v1.3`, esta execução volta direto ao branch padrão (`master`).
O Nextflow mantém uma cópia de trabalho local separada para cada revisão que você usou, que é o que os marcadores `>` no `nextflow info` mostram, mas nunca lembra qual foi a última que você executou.
Você pode descobrir o nome do branch padrão de um pipeline executando `nextflow info <pipeline>`; é o marcado como `(default)`.
A reprodutibilidade é inteiramente sua responsabilidade: sempre passe `-r` explicitamente quando isso for importante, em vez de assumir que uma revisão fixada em uma execução anterior ainda se aplica.

### Conclusão

Você sabe como fixar um pipeline remoto em uma versão, branch ou commit específico para execução reprodutível, e que a fixação se aplica apenas àquela invocação, não a execuções posteriores.

### O que vem a seguir?

Você cobriu os fundamentos de execução e gerenciamento de pipelines Nextflow.
Consulte o [Resumo do curso](next_steps.md) para saber o que fazer a partir daqui.

---

## Resumo

Nesta parte você aprendeu a:

- Executar um pipeline diretamente de um repositório GitHub sem baixá-lo
- Fixar um pipeline remoto em uma revisão específica para reprodutibilidade, e entender que a fixação se aplica apenas àquela invocação
