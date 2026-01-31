# Parte 4: Adicionando testes

Na primeira parte deste curso, você construiu um pipeline de chamada de variantes que era completamente linear e processava os dados de cada amostra independentemente das outras.

Na segunda parte, mostramos como usar canais e operadores de canal para implementar chamada conjunta de variantes com GATK.

Na terceira parte, modularizamos o pipeline.

Nesta parte do treinamento, vamos mostrar como usar [**nf-test**](https://www.nf-test.com/), um framework de testes que se integra bem com Nextflow e facilita a adição de testes tanto em nível de módulo quanto em nível de fluxo de trabalho ao seu pipeline. Para acompanhar esta parte do treinamento, você deve ter completado a Parte 1, Parte 2 e Parte 3, bem como a [missão secundária nf-test](../../side_quests/nf-test.md), que cobre os fundamentos do nf-test e por que testar é importante.

---

## 0. Aquecimento

!!! note "Nota"

    Certifique-se de estar no diretório de trabalho correto:
    `cd /workspaces/training/nf4-science/genomics`

Se você trabalhou nas partes anteriores deste curso de treinamento, você deve ter uma versão funcional do pipeline de genômica com a estrutura de diretórios de módulos apropriada.

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
    ```

Este diretório de módulos pode ser encontrado no diretório `solutions` se você precisar.

Vamos começar com o mesmo fluxo de trabalho da Parte 3, que fornecemos para você no arquivo `genomics-4.nf`. Exatamente como na [missão secundária nf-test](../../side_quests/nf-test.md), vamos adicionar alguns tipos diferentes de testes aos três processos neste pipeline, bem como um teste em nível de fluxo de trabalho.

### 0.1. Verificar se o fluxo de trabalho executa

Antes de começarmos a adicionar testes, certifique-se de que o fluxo de trabalho executa conforme esperado.

```bash
nextflow run genomics-4.nf -resume
```

Isso deve parecer muito familiar agora se você está acompanhando este curso de treinamento desde o início.

??? success "Saída do comando"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `genomics-4.nf` [gloomy_poincare] DSL2 - revision: 43203316e0

    executor >  local (7)
    [18/89dfa4] SAMTOOLS_INDEX (1)       | 3 of 3 ✔
    [30/b2522b] GATK_HAPLOTYPECALLER (2) | 3 of 3 ✔
    [a8/d2c189] GATK_JOINTGENOTYPING     | 1 of 1 ✔
    ```

Como anteriormente, agora haverá um diretório `work` e um diretório `results_genomics` dentro do seu diretório de projeto. Na verdade, usaremos esses resultados mais tarde em nossos testes. Mas a partir de agora vamos usar o pacote `nf-test` para testar o pipeline.

### 0.2. Inicializar `nf-test`

Como na [missão secundária nf-test](../../side_quests/nf-test.md), precisamos inicializar o pacote `nf-test`.

```bash
nf-test init
```

??? success "Saída do comando"

    ```console
    🚀 nf-test 0.9.3
    https://www.nf-test.com
    (c) 2021 - 2024 Lukas Forer and Sebastian Schoenherr

    Project configured. Configuration is stored in nf-test.config
    ```

??? abstract "Conteúdo do nf-test.config"

    ```groovy title="nf-test.config"
    config {

        testsDir "tests"
        workDir ".nf-test"
        configFile "tests/nextflow.config"
        profile ""

    }
    ```

Ele também cria um diretório `tests` contendo um esboço de arquivo de configuração.

### Resumo

Agora estamos prontos para começar a escrever testes para nosso pipeline de genômica.

### O que vem a seguir?

Escrever testes básicos que avaliem se as chamadas de processo foram bem-sucedidas e produziram as saídas corretas.

---

## 1. Testar um processo quanto a sucesso e saídas correspondentes

Começaremos testando o processo `SAMTOOLS_INDEX`, que cria arquivos de índice para arquivos BAM para permitir acesso aleatório eficiente. Este é um bom primeiro caso de teste porque:

1. Possui uma única entrada bem definida (um arquivo BAM)
2. Produz uma saída previsível (um arquivo de índice BAI)
3. A saída deve ser idêntica para entradas idênticas

### 1.1. Gerar um esboço de arquivo de teste

Primeiro, gere um esboço de arquivo de teste:

```bash
nf-test generate process modules/samtools/index/main.nf
```

??? success "Saída do comando"

    ```console
    Load source file '/workspaces/training/nf4-science/genomics/modules/samtools/index/main.nf'
    Wrote process test file '/workspaces/training/nf4-science/genomics/tests/modules/samtools/index/main.nf.test'

    SUCCESS: Generated 1 test files.
    ```

Isso cria um arquivo no mesmo diretório que `main.nf`.
Você pode navegar até o diretório no explorador de arquivos e abrir o arquivo, que deve conter o seguinte código:

```groovy title="tests/modules/samtools/index/main.nf.test" linenums="1"
nextflow_process {

    name "Test Process SAMTOOLS_INDEX"
    script "modules/samtools/index/main.nf"
    process "SAMTOOLS_INDEX"

    test("Should run without failures") {

        when {
            params {
                // define parameters here. Example:
                // outdir = "tests/results"
            }
            process {
                """
                // define inputs of the process here. Example:
                // input[0] = file("test-file.txt")
                """
            }
        }

        then {
            assert process.success
            assert snapshot(process.out).match()
        }

    }

}
```

As asserções iniciais devem ser familiares da [missão secundária nf-test](../../side_quests/nf-test.md):

- `assert process.success` declara que esperamos que o processo execute com sucesso e complete sem falhas.
- `snapshot(process.out).match()` declara que esperamos que o resultado da execução seja idêntico ao resultado obtido em uma execução anterior (se aplicável).
  Discutimos isso em mais detalhes posteriormente.

Usando isso como ponto de partida, precisamos adicionar as entradas de teste corretas para o processo de índice do samtools, e quaisquer parâmetros se aplicável.

### 1.2. Mover o arquivo de teste e atualizar o caminho do script

Antes de começarmos a trabalhar no preenchimento do teste, precisamos mover o arquivo para sua localização definitiva. Parte da razão pela qual adicionamos um diretório para cada módulo é que agora podemos enviar testes em um diretório `tests` co-localizado com o arquivo `main.nf` de cada módulo. Crie esse diretório e mova o arquivo de teste para lá.

```bash
mkdir -p modules/samtools/index/tests
mv tests/modules/samtools/index/main.nf.test modules/samtools/index/tests/
```

Agora podemos simplificar a seção `script` do arquivo de teste para um caminho relativo:

=== "Depois"

    ```groovy title="modules/samtools/index/tests/main.nf.test" linenums="3" hl_lines="2"
    name "Test Process SAMTOOLS_INDEX"
    script "../main.nf"
    process "SAMTOOLS_INDEX"
    ```

=== "Antes"

    ```groovy title="modules/samtools/index/tests/main.nf.test" linenums="3" hl_lines="2"
    name "Test Process SAMTOOLS_INDEX"
    script "modules/samtools/index/main.nf"
    process "SAMTOOLS_INDEX"
    ```

Isso informa ao teste onde encontrar o arquivo `main.nf` do módulo, sem ter que especificar o caminho completo.

### 1.3. Fornecer entradas de teste para SAMTOOLS_INDEX

O arquivo esboço inclui um placeholder que precisamos substituir por uma entrada de teste real, apropriada para a entrada do `samtools index`. A entrada apropriada é um arquivo BAM, que temos disponível no diretório `data/bam`.

=== "Depois"

    ```groovy title="modules/samtools/index/tests/main.nf.test" linenums="14"
    process {
        """
        input[0] = file("${projectDir}/data/bam/reads_son.bam")
        """
    }
    ```

=== "Antes"

    ```groovy title="modules/samtools/index/tests/main.nf.test" linenums="14"
    process {
        """
        // define inputs of the process here. Example:
        // input[0] = file("test-file.txt")
        """
    }
    ```

### 1.4. Nomear o teste baseado na funcionalidade

Como aprendemos antes, é uma boa prática renomear o teste para algo que faça sentido no contexto do teste.

=== "Depois"

    ```groovy title="modules/samtools/index/tests/main.nf.test" linenums="7"
    test("Should index reads_son.bam correctly") {
    ```

    Isso recebe uma string arbitrária, então poderíamos colocar qualquer coisa que quisermos.
    Aqui escolhemos nos referir ao nome do arquivo e seu formato.

=== "Antes"

    ```groovy title="modules/samtools/index/tests/main.nf.test" linenums="7"
    test("Should run without failures") {
    ```

### 1.5. Executar o teste e examinar a saída

Execute o teste:

```bash
nf-test test modules/samtools/index/tests/main.nf.test
```

??? success "Saída do comando"

    ```console
    🚀 nf-test 0.9.3
    https://www.nf-test.com
    (c) 2021 - 2024 Lukas Forer and Sebastian Schoenherr


    Test Process SAMTOOLS_INDEX

      Test [625e39ee] 'Should index reads_son.bam correctly' PASSED (7.717s)
      Snapshots:
        1 created [Should index reads_son.bam correctly]


    Snapshot Summary:
      1 created

    SUCCESS: Executed 1 tests in 7.727s
    ```

Como aprendemos anteriormente, isso verificou a asserção básica sobre o sucesso do processo e criou um arquivo de snapshot baseado na saída do processo. Podemos ver o conteúdo do arquivo de snapshot em `tests/modules/samtools/index/tests/main.nf.test.snap`:

```json title="modules/samtools/index/tests/main.nf.test.snap" linenums="1"
{
  "Should index reads_son.bam correctly": {
    "content": [
      {
        "0": [
          [
            "reads_son.bam:md5,af5956d9388ba017944bef276b71d809",
            "reads_son.bam.bai:md5,a2ca7b84998218ee77eff14af8eb8ca2"
          ]
        ]
      }
    ],
    "meta": {
      "nf-test": "0.9.3",
      "nextflow": "25.10.2"
    },
    "timestamp": "2026-01-27T15:09:48.394063389"
  }
}
```

Também podemos executar o teste novamente e ver que ele passa, porque a saída é idêntica ao snapshot:

```bash
nf-test test modules/samtools/index/tests/main.nf.test
```

??? success "Saída do comando"

    ```console
    🚀 nf-test 0.9.3
    https://www.nf-test.com
    (c) 2021 - 2024 Lukas Forer and Sebastian Schoenherr


    Test Process SAMTOOLS_INDEX

      Test [625e39ee] 'Should index reads_son.bam correctly' PASSED (7.938s)


    SUCCESS: Executed 1 tests in 7.987s
    ```

### 1.6. Adicionar mais testes ao SAMTOOLS_INDEX

Às vezes é útil testar uma variedade de diferentes arquivos de entrada para garantir que estamos testando uma variedade de possíveis problemas. Adicione testes para os arquivos BAM da mãe e do pai no trio de nossos dados de teste.

```groovy
    test("Should index reads_mother.bam correctly") {

        when {
            params {
                // define parameters here
            }
            process {
                """
                input[0] = file("${projectDir}/data/bam/reads_mother.bam")
                """
            }
        }

        then {
            assert process.success
            assert snapshot(process.out).match()
        }

    }

    test("Should index reads_father.bam correctly") {

        when {
            params {
                // define parameters here
            }
            process {
                """
                input[0] = file("${projectDir}/data/bam/reads_father.bam")
                """
            }
        }

        then {
            assert process.success
            assert snapshot(process.out).match()
        }

    }
```

Então você pode executar o teste novamente:

```bash
nf-test test modules/samtools/index/tests/main.nf.test
```

??? success "Saída do comando"

    ```console
    🚀 nf-test 0.9.3
    https://www.nf-test.com
    (c) 2021 - 2024 Lukas Forer and Sebastian Schoenherr


    Test Process SAMTOOLS_INDEX

      Test [625e39ee] 'Should index reads_son.bam correctly' PASSED (7.185s)
      Test [a8b28f36] 'Should index reads_mother.bam correctly' PASSED (6.576s)
      Test [c15852a1] 'Should index reads_father.bam correctly' PASSED (6.31s)
      Snapshots:
        2 created [Should index reads_father.bam correctly, Should index reads_mother.bam correctly]


    Snapshot Summary:
      2 created

    SUCCESS: Executed 3 tests in 20.117s
    ```

Observe o aviso, referindo-se ao efeito do parâmetro `--update-snapshot`.

!!! note "Nota"

    Aqui estamos usando dados de teste que usamos anteriormente para demonstrar as saídas científicas do pipeline.
    Se tivéssemos planejado operar esses testes em um ambiente de produção, teríamos gerado entradas menores para fins de teste.

    Em geral, é importante manter os testes unitários o mais leves possível usando os menores pedaços de dados necessários e suficientes para avaliar a funcionalidade do processo, caso contrário o tempo total de execução pode aumentar consideravelmente.
    Uma suíte de testes que leva muito tempo para executar regularmente é uma suíte de testes que provavelmente será pulada no interesse da conveniência.

### Resumo

Você escreveu seu primeiro teste de módulo para um processo de genômica, verificando que `SAMTOOLS_INDEX` cria corretamente arquivos de índice para diferentes arquivos BAM. A suíte de testes garante que:

1. O processo executa com sucesso
2. Arquivos de índice são criados
3. As saídas são consistentes entre execuções
4. O processo funciona para todos os arquivos BAM de amostra

### O que vem a seguir?

Aprenda como escrever testes para outros processos em nosso fluxo de trabalho de genômica, usando o método setup para lidar com processos encadeados. Também avaliaremos se as saídas, especificamente nossos arquivos VCF, contêm chamadas de variantes esperadas.

---

## 2. Adicionar testes a um processo encadeado e testar o conteúdo

Para testar `GATK_HAPLOTYPECALLER`, precisamos fornecer ao processo a saída de `SAMTOOLS_INDEX` como entrada. Poderíamos fazer isso executando `SAMTOOLS_INDEX`, recuperando suas saídas e armazenando-as com os dados de teste do fluxo de trabalho. Essa é na verdade a abordagem recomendada para um pipeline polido, mas nf-test fornece uma abordagem alternativa, usando o método `setup`.

Com o método setup, podemos acionar o processo `SAMTOOLS_INDEX` como parte da configuração do teste, e então usar sua saída como entrada para `GATK_HAPLOTYPECALLER`. Isso tem um custo: vamos ter que executar o processo `SAMTOOLS_INDEX` toda vez que executarmos o teste para `GATK_HAPLOTYPECALLER`. No entanto, talvez ainda estejamos desenvolvendo o fluxo de trabalho e não queiramos pré-gerar dados de teste que talvez tenhamos que mudar mais tarde. O processo `SAMTOOLS_INDEX` também é muito rápido, então talvez os benefícios de pré-gerar e armazenar suas saídas sejam negligenciáveis. Aqui está como o método setup funciona.

### 2.1. Gerar e colocar o arquivo de teste

Como anteriormente, primeiro geramos o esboço do arquivo:

```bash
nf-test generate process modules/gatk/haplotypecaller/main.nf
```

??? success "Saída do comando"

    ```console
    Load source file '/workspaces/training/nf4-science/genomics/modules/gatk/haplotypecaller/main.nf'
    Wrote process test file '/workspaces/training/nf4-science/genomics/tests/modules/gatk/haplotypecaller/main.nf.test'

    SUCCESS: Generated 1 test files.
    ```

Isso produz o seguinte esboço de teste:

```groovy title="tests/modules/gatk/haplotypecaller/main.nf.test" linenums="1"
nextflow_process {

    name "Test Process GATK_HAPLOTYPECALLER"
    script "modules/gatk/haplotypecaller/main.nf"
    process "GATK_HAPLOTYPECALLER"

    test("Should run without failures") {

        when {
            params {
                // define parameters here. Example:
                // outdir = "tests/results"
            }
            process {
                """
                // define inputs of the process here. Example:
                // input[0] = file("test-file.txt")
                """
            }
        }

        then {
            assert process.success
            assert snapshot(process.out).match()
        }

    }

}
```

### 2.2. Mover o arquivo de teste e atualizar o caminho do script

Criamos um diretório para o arquivo de teste co-localizado com o arquivo `main.nf` do módulo:

```bash
mkdir -p modules/gatk/haplotypecaller/tests
```

E movemos o arquivo esboço de teste para lá:

```bash
mv tests/modules/gatk/haplotypecaller/main.nf.test modules/gatk/haplotypecaller/tests/
```

Finalmente, não esqueça de atualizar o caminho do script:

=== "Depois"

    ```groovy title="modules/gatk/haplotypecaller/tests/main.nf.test" linenums="3" hl_lines="2"
        name "Test Process GATK_HAPLOTYPECALLER"
        script "../main.nf"
        process "GATK_HAPLOTYPECALLER"
    ```

=== "Antes"

    ```groovy title="modules/gatk/haplotypecaller/tests/main.nf.test" linenums="3" hl_lines="2"
        name "Test Process GATK_HAPLOTYPECALLER"
        script "modules/gatk/haplotypecaller/main.nf"
        process "GATK_HAPLOTYPECALLER"
    ```

### 2.3. Fornecer entradas usando o método setup

Inserimos um bloco `setup` antes do bloco `when`, onde podemos acionar uma execução do processo `SAMTOOLS_INDEX` em um de nossos arquivos de entrada originais. Além disso, lembre-se como antes de mudar o nome do teste para algo significativo.

=== "Depois"

    ```groovy title="modules/gatk/haplotypecaller/tests/main.nf.test" linenums="7" hl_lines="1-12"
        test("Should call son's haplotype correctly") {

            setup {
                run("SAMTOOLS_INDEX") {
                    script "../../../samtools/index/main.nf"
                    process {
                        """
                        input[0] =  file("${projectDir}/data/bam/reads_son.bam")
                        """
                    }
                }
            }

            when {
    ```

=== "Antes"

    ```groovy title="modules/gatk/haplotypecaller/tests/main.nf.test" linenums="7"  hl_lines="1"
    test("Should run without failures") {

        when {
    ```

Então podemos nos referir à saída desse processo no bloco `when` onde especificamos as entradas do teste:

```groovy title="modules/gatk/haplotypecaller/tests/main.nf.test" linenums="20"
        when {
            params {
                // define parameters here
            }
            process {
                """
                input[0] = SAMTOOLS_INDEX.out
                input[1] = file("${projectDir}/data/ref/ref.fasta")
                input[2] = file("${projectDir}/data/ref/ref.fasta.fai")
                input[3] = file("${projectDir}/data/ref/ref.dict")
                input[4] = file("${projectDir}/data/ref/intervals.bed")
                """
            }
        }
```

Faça essa mudança e execute o teste novamente:

```bash
nf-test test modules/gatk/haplotypecaller/tests/main.nf.test
```

??? success "Saída do comando"

    ```console
    🚀 nf-test 0.9.3
    https://www.nf-test.com
    (c) 2021 - 2024 Lukas Forer and Sebastian Schoenherr


    Test Process GATK_HAPLOTYPECALLER

      Test [c5156c2b] 'Should call son's haplotype correctly' PASSED (40.53s)
      Snapshots:
        1 created [Should call son's haplotype correctly]


    Snapshot Summary:
      1 created

    SUCCESS: Executed 1 tests in 40.555s
    ```

Ele também produz um arquivo de snapshot como anteriormente.

### 2.4. Executar novamente e observar a falha

Curiosamente, se você executar o mesmo comando novamente, desta vez o teste falhará.

```bash
nf-test test modules/gatk/haplotypecaller/tests/main.nf.test
```

??? failure "Saída do comando"

    ```console
    🚀 nf-test 0.9.3
    https://www.nf-test.com
    (c) 2021 - 2024 Lukas Forer and Sebastian Schoenherr


    Test Process GATK_HAPLOTYPECALLER

      Test [c5156c2b] 'Should call son's haplotype correctly' FAILED (40.123s)

      java.lang.RuntimeException: Different Snapshot:
      [                                                                                           [
          {                                                                                           {
              "0": [                                                                                      "0": [
                  "reads_son.bam.g.vcf:md5,069316cdd4328542ffc6ae247b1dac39"                         |                 "reads_son.bam.g.vcf:md5,005f1a13ee39f11b0fc9bea094850eac"
              ],                                                                                          ],
              "1": [                                                                                      "1": [
                  "reads_son.bam.g.vcf.idx:md5,dc36c18f2afdc546f41e68b2687e9334"                     |                 "reads_son.bam.g.vcf.idx:md5,dbad4b76a4b90c158ffc9c9740764242"
              ],                                                                                          ],
              "idx": [                                                                                    "idx": [
                  "reads_son.bam.g.vcf.idx:md5,dc36c18f2afdc546f41e68b2687e9334"                     |                 "reads_son.bam.g.vcf.idx:md5,dbad4b76a4b90c158ffc9c9740764242"
              ],                                                                                          ],
              "vcf": [                                                                                    "vcf": [
                  "reads_son.bam.g.vcf:md5,069316cdd4328542ffc6ae247b1dac39"                         |                 "reads_son.bam.g.vcf:md5,005f1a13ee39f11b0fc9bea094850eac"
              ]                                                                                           ]
          }                                                                                           }
      ]                                                                                           ]

      Nextflow stdout:

      Nextflow stderr:


        Obsolete snapshots can only be checked if all tests of a file are executed successful.


    FAILURE: Executed 1 tests in 40.156s (1 failed)
    ```

A mensagem de erro informa que houve diferenças entre os snapshots para as duas execuções; especificamente, os valores de md5sum são diferentes para os arquivos VCF.

Por quê? Para encurtar a história, a ferramenta HaplotypeCaller inclui um timestamp no cabeçalho VCF que é diferente toda vez (por definição).
Como resultado, não podemos simplesmente esperar que os arquivos tenham md5sums idênticos mesmo se tiverem conteúdo idêntico em termos das próprias chamadas de variantes.

Como lidamos com isso?

### 2.5. Usar um método de asserção de conteúdo para verificar uma variante específica

Uma maneira de resolver o problema é usar um [tipo diferente de asserção](https://nf-co.re/docs/contributing/tutorials/nf-test_assertions).
Neste caso, vamos verificar conteúdo específico em vez de afirmar identidade.
Mais exatamente, faremos a ferramenta ler as linhas do arquivo VCF e verificar a existência de linhas específicas.

Na prática, substituímos a segunda asserção no bloco `then` da seguinte forma:

=== "Depois"

    ```groovy title="modules/gatk/haplotypecaller/tests/main.nf.test" linenums="35" hl_lines="3 4"
            then {
                assert process.success
                assert path(process.out[0][0]).readLines().contains('#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	reads_son')
                assert path(process.out[0][0]).readLines().contains('20_10037292_10066351	3277	.	G	<NON_REF>	.	.	END=3282	GT:DP:GQ:MIN_DP:PL	0/0:25:72:24:0,72,719')
            }
    ```

=== "Antes"

    ```groovy title="modules/gatk/haplotypecaller/tests/main.nf.test" linenums="35" hl_lines="3"
    then {
        assert process.success
        assert snapshot(process.out).match()
    }
    ```

Aqui estamos lendo o conteúdo completo do arquivo de saída VCF e procurando uma correspondência de conteúdo, o que é aceitável fazer em um pequeno arquivo de teste, mas você não gostaria de fazer isso em um arquivo maior.
Você pode optar por ler linhas específicas.

Esta abordagem requer escolher com mais cuidado o que queremos usar como 'sinal' para testar.
Pelo lado positivo, ela pode ser usada para testar com grande precisão se uma ferramenta de análise pode identificar consistentemente características 'difíceis' (como variantes raras) à medida que passa por desenvolvimento adicional.

### 2.6. Executar novamente e observar sucesso

Uma vez que modificamos o teste desta forma, podemos executar o teste várias vezes, e ele passará consistentemente.

```bash
nf-test test modules/gatk/haplotypecaller/tests/main.nf.test
```

??? success "Saída do comando"

    ```console
    🚀 nf-test 0.9.3
    https://www.nf-test.com
    (c) 2021 - 2024 Lukas Forer and Sebastian Schoenherr


    Test Process GATK_HAPLOTYPECALLER

      Test [c5156c2b] 'Should call son's haplotype correctly' PASSED (40.53s)


    SUCCESS: Executed 1 tests in 40.555s
    ```

### 2.7. Adicionar mais testes

Adicione testes similares para as amostras da mãe e do pai:

```groovy title="modules/gatk/haplotypecaller/tests/main.nf.test" linenums="43"
    test("Should call mother's haplotype correctly") {

        setup {
            run("SAMTOOLS_INDEX") {
                script "../../../samtools/index/main.nf"
                process {
                    """
                    input[0] =  file("${projectDir}/data/bam/reads_mother.bam")
                    """
                }
            }
        }

        when {
            params {
                // define parameters here
            }
            process {
                """
                input[0] = SAMTOOLS_INDEX.out
                input[1] = file("${projectDir}/data/ref/ref.fasta")
                input[2] = file("${projectDir}/data/ref/ref.fasta.fai")
                input[3] = file("${projectDir}/data/ref/ref.dict")
                input[4] = file("${projectDir}/data/ref/intervals.bed")
                """
            }
        }

        then {
            assert process.success
            assert path(process.out[0][0]).readLines().contains('#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	reads_mother')
            assert path(process.out[0][0]).readLines().contains('20_10037292_10066351	3277	.	G	<NON_REF>	.	.	END=3278	GT:DP:GQ:MIN_DP:PL	0/0:38:99:37:0,102,1530')
        }
    }

    test("Should call father's haplotype correctly") {

        setup {
            run("SAMTOOLS_INDEX") {
                script "../../../samtools/index/main.nf"
                process {
                    """
                    input[0] =  file("${projectDir}/data/bam/reads_father.bam")
                    """
                }
            }
        }

        when {
            params {
                // define parameters here
            }
            process {
                """
                input[0] = SAMTOOLS_INDEX.out
                input[1] = file("${projectDir}/data/ref/ref.fasta")
                input[2] = file("${projectDir}/data/ref/ref.fasta.fai")
                input[3] = file("${projectDir}/data/ref/ref.dict")
                input[4] = file("${projectDir}/data/ref/intervals.bed")
                """
            }
        }

        then {
            assert process.success
            assert path(process.out[0][0]).readLines().contains('#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	reads_father')
            assert path(process.out[0][0]).readLines().contains('20_10037292_10066351	3277	.	G	<NON_REF>	.	.	END=3281	GT:DP:GQ:MIN_DP:PL	0/0:44:99:42:0,120,1800')
        }
    }
```

### 2.8. Executar o comando de teste

```bash
nf-test test modules/gatk/haplotypecaller/tests/main.nf.test
```

??? success "Saída do comando"

    ```console
    🚀 nf-test 0.9.3
    https://www.nf-test.com
    (c) 2021 - 2024 Lukas Forer and Sebastian Schoenherr


    Test Process GATK_HAPLOTYPECALLER

      Test [c5156c2b] 'Should call son's haplotype correctly' PASSED (40.53s)
      Test [10de94a8] 'Should call mother's haplotype correctly' PASSED (41.47s)
      Test [c0386fc7] 'Should call father's haplotype correctly' PASSED (45.556s)


    SUCCESS: Executed 3 tests in 127.586s
    ```

Isso completa o plano de testes básico para esta segunda etapa do pipeline. Vamos para o terceiro e último teste em nível de módulo!

### Resumo

Você aprendeu como:

1. Testar processos que dependem de saídas de outros processos
2. Verificar variantes genômicas específicas em arquivos de saída VCF
3. Lidar com saídas não determinísticas verificando conteúdo específico
4. Testar chamada de variantes em várias amostras

### O que vem a seguir?

Aprenda como escrever testes que usam dados de teste pré-gerados para a etapa de genotipagem conjunta.

---

## 3. Usar dados de teste pré-gerados

Para a etapa de genotipagem conjunta, usaremos uma abordagem diferente - usando dados de teste pré-gerados. Isso geralmente é preferível para:

1. Processos complexos com múltiplas dependências
2. Processos que levam muito tempo para executar
3. Processos que fazem parte de um pipeline estável e de produção

### 3.1. Gerar dados de teste

Inspecione os resultados que geramos no início desta seção:

```bash
tree results_genomics/
```

```console title="Conteúdo do diretório de resultados"
results_genomics/
├── family_trio.joint.vcf
├── family_trio.joint.vcf.idx
├── gvcf
│   ├── reads_father.bam.g.vcf -> /workspaces/training/nf4-science/genomics/work/30/b2522b83c63baff8c3cf75704512a2/reads_father.bam.g.vcf
│   ├── reads_father.bam.g.vcf.idx -> /workspaces/training/nf4-science/genomics/work/30/b2522b83c63baff8c3cf75704512a2/reads_father.bam.g.vcf.idx
│   ├── reads_mother.bam.g.vcf -> /workspaces/training/nf4-science/genomics/work/f6/be2efa58e625d08cf8d0da1d0e9f09/reads_mother.bam.g.vcf
│   ├── reads_mother.bam.g.vcf.idx -> /workspaces/training/nf4-science/genomics/work/f6/be2efa58e625d08cf8d0da1d0e9f09/reads_mother.bam.g.vcf.idx
│   ├── reads_son.bam.g.vcf -> /workspaces/training/nf4-science/genomics/work/fe/2f22d56aa16ed45f8bc419312894f6/reads_son.bam.g.vcf
│   └── reads_son.bam.g.vcf.idx -> /workspaces/training/nf4-science/genomics/work/fe/2f22d56aa16ed45f8bc419312894f6/reads_son.bam.g.vcf.idx
└── indexed_bam
    ├── reads_father.bam -> /workspaces/training/nf4-science/genomics/work/42/a3bf19dbfaf1f3672b16a5d5e6a8be/reads_father.bam
    ├── reads_father.bam.bai -> /workspaces/training/nf4-science/genomics/work/cf/289c2d264f496d60a69e3e9ba6463e/reads_father.bam.bai
    ├── reads_mother.bam -> /workspaces/training/nf4-science/genomics/work/af/f31a6ade82cc0cf853c4f61c8bc473/reads_mother.bam
    ├── reads_mother.bam.bai -> /workspaces/training/nf4-science/genomics/work/18/89dfa40a3def17e45421e54431a126/reads_mother.bam.bai
    ├── reads_son.bam -> /workspaces/training/nf4-science/genomics/work/9f/9615dd553d6f13d8bec4f006ac395f/reads_son.bam
    └── reads_son.bam.bai -> /workspaces/training/nf4-science/genomics/work/4d/cb384a97db5687cc9daab002017c7c/reads_son.bam.bai

2 directories, 14 files
```

A etapa de genotipagem conjunta precisa dos arquivos VCF produzidos pelas etapas do haplotype caller como entradas, junto com os índices. Então vamos copiar os resultados que temos para o diretório de testes do módulo `jointgenotyping`.

```bash
mkdir -p modules/gatk/jointgenotyping/tests/inputs/
cp results_genomics/gvcf/*.g.vcf results_genomics/gvcf/*.g.vcf.idx modules/gatk/jointgenotyping/tests/inputs/
```

Agora podemos usar esses arquivos como entradas para o teste que vamos escrever para a etapa de genotipagem conjunta.

### 3.2. Gerar o esboço do arquivo de teste

Como anteriormente, primeiro geramos o esboço do arquivo:

```bash
nf-test generate process modules/gatk/jointgenotyping/main.nf
```

??? success "Saída do comando"

    ```console
    Load source file '/workspaces/training/nf4-science/genomics/modules/gatk/jointgenotyping/main.nf'
    Wrote process test file '/workspaces/training/nf4-science/genomics/tests/modules/gatk/jointgenotyping/main.nf.test'

    SUCCESS: Generated 1 test files.
    ```

Isso produz o seguinte esboço de teste:

```groovy title="tests/modules/gatk/jointgenotyping/main.nf.test" linenums="1"
nextflow_process {

    name "Test Process GATK_JOINTGENOTYPING"
    script "modules/gatk/jointgenotyping/main.nf"
    process "GATK_JOINTGENOTYPING"

    test("Should run without failures") {

        when {
            params {
                // define parameters here. Example:
                // outdir = "tests/results"
            }
            process {
                """
                // define inputs of the process here. Example:
                // input[0] = file("test-file.txt")
                """
            }
        }

        then {
            assert process.success
            assert snapshot(process.out).match()
        }

    }

}
```

### 3.3. Mover o arquivo de teste e atualizar o caminho do script

Desta vez já temos um diretório para testes co-localizado com o arquivo `main.nf` do módulo, então podemos mover o arquivo esboço de teste para lá:

```bash
mv tests/modules/gatk/jointgenotyping/main.nf.test modules/gatk/jointgenotyping/tests/
```

E não esqueça de atualizar o caminho do script:

=== "Depois"

    ```groovy title="modules/gatk/jointgenotyping/tests/main.nf.test" linenums="3" hl_lines="2"
    name "Test Process GATK_JOINTGENOTYPING"
    script "../main.nf"
    process "GATK_JOINTGENOTYPING"
    ```

=== "Antes"

    ```groovy title="modules/gatk/jointgenotyping/tests/main.nf.test" linenums="3" hl_lines="2"
    name "Test Process GATK_JOINTGENOTYPING"
    script "modules/gatk/jointgenotyping/main.nf"
    process "GATK_JOINTGENOTYPING"
    ```

### 3.4. Fornecer entradas

Preencha as entradas com base nas definições de entrada do processo e renomeie o teste adequadamente:

```groovy title="modules/gatk/jointgenotyping/tests/main.nf.test" linenums="7"
    test("Should call trio's joint genotype correctly") {

        when {
            params {
                // define parameters here
            }
            process {
                """
                input[0] = [
                    file("${projectDir}/modules/gatk/jointgenotyping/tests/inputs/reads_father.bam.g.vcf"),
                    file("${projectDir}/modules/gatk/jointgenotyping/tests/inputs/reads_mother.bam.g.vcf"),
                    file("${projectDir}/modules/gatk/jointgenotyping/tests/inputs/reads_son.bam.g.vcf")
                ]
                input[1] = [
                    file("${projectDir}/modules/gatk/jointgenotyping/tests/inputs/reads_father.bam.g.vcf.idx"),
                    file("${projectDir}/modules/gatk/jointgenotyping/tests/inputs/reads_mother.bam.g.vcf.idx"),
                    file("${projectDir}/modules/gatk/jointgenotyping/tests/inputs/reads_son.bam.g.vcf.idx")
                ]
                input[2] = file("${projectDir}/data/ref/intervals.bed")
                input[3] = "family_trio"
                input[4] = file("${projectDir}/data/ref/ref.fasta")
                input[5] = file("${projectDir}/data/ref/ref.fasta.fai")
                input[6] = file("${projectDir}/data/ref/ref.dict")
                """
            }
        }
```

### 3.5. Usar asserções de conteúdo

A saída da etapa de genotipagem conjunta é outro arquivo VCF, então vamos usar uma asserção de conteúdo novamente.

```groovy title="modules/gatk/jointgenotyping/tests/main.nf.test" linenums="25"
    then {
        assert process.success
        assert path(process.out[0][0]).readLines().contains('#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	reads_father	reads_mother	reads_son')
        assert path(process.out[0][0]).readLines().contains('20_10037292_10066351	3480	.	C	CT	1625.89	.	AC=5;AF=0.833;AN=6;BaseQRankSum=0.220;DP=85;ExcessHet=0.0000;FS=2.476;MLEAC=5;MLEAF=0.833;MQ=60.00;MQRankSum=0.00;QD=21.68;ReadPosRankSum=-1.147e+00;SOR=0.487	GT:AD:DP:GQ:PL	0/1:15,16:31:99:367,0,375	1/1:0,18:18:54:517,54,0	1/1:0,26:26:78:756,78,0')
    }
```

Ao verificar o conteúdo de uma variante específica no arquivo de saída, este teste verifica que:

1. O processo de genotipagem conjunta executa com sucesso
2. O VCF de saída contém todas as três amostras na ordem correta
3. Uma variante específica é chamada corretamente com:
   - Genótipos precisos para cada amostra (0/1 para o pai, 1/1 para a mãe e o filho)
   - Profundidades de leitura e qualidades de genótipo corretas
   - Estatísticas em nível populacional como frequência alélica (AF=0.833)

Não fizemos snapshot do arquivo inteiro, mas ao verificar uma variante específica, podemos ter confiança de que o processo de genotipagem conjunta está funcionando conforme esperado.

### 3.6. Executar o teste

```bash
nf-test test modules/gatk/jointgenotyping/tests/main.nf.test
```

??? success "Saída do comando"

    ```console
    🚀 nf-test 0.9.3
    https://www.nf-test.com
    (c) 2021 - 2024 Lukas Forer and Sebastian Schoenherr


    Test Process GATK_JOINTGENOTYPING

      Test [ac2067de] 'Should call trio's joint genotype correctly' PASSED (53.827s)


    SUCCESS: Executed 1 tests in 53.837s
    ```

O teste passa, verificando que nosso processo de genotipagem conjunta corretamente:

1. Combina VCFs de amostras individuais
2. Realiza chamada conjunta de variantes
3. Produz um VCF multi-amostra com chamadas de genótipo consistentes entre execuções

### Resumo

Você sabe como:

- Usar resultados previamente gerados como entradas para testes
- Escrever testes usando dados de teste pré-gerados

### O que vem a seguir?

Adicionar um teste em nível de fluxo de trabalho para verificar se todo o pipeline de chamada de variantes funciona de ponta a ponta.

---

## 4. Adicionar um teste em nível de fluxo de trabalho

Agora testaremos o pipeline completo de chamada de variantes, de arquivos BAM a genótipos conjuntos. Isso verifica que:

1. Todos os processos trabalham juntos corretamente
2. Os dados fluem adequadamente entre as etapas
3. As chamadas de variantes finais são consistentes

### 4.1. Gerar o teste de fluxo de trabalho

Gere um arquivo de teste para o pipeline completo:

```bash
nf-test generate pipeline genomics-4.nf
```

??? success "Saída do comando"

    ```console
    Load source file '/workspaces/training/nf4-science/genomics/genomics-4.nf'
    Wrote pipeline test file '/workspaces/training/nf4-science/genomics/tests/genomics-4.nf.test'

    SUCCESS: Generated 1 test files.
    ```

Isso cria um esboço de teste básico:

```groovy title="tests/genomics-4.nf.test" linenums="1"
nextflow_pipeline {

    name "Test Workflow genomics-4.nf"
    script "genomics-4.nf"

    test("Should run without failures") {

        when {
            params {
                // define parameters here. Example:
                // outdir = "tests/results"
            }
        }

        then {
            assert workflow.success
        }

    }

}
```

Apenas corrija o nome para algo significativo (você verá por que isso é útil em breve).

=== "Depois"

    ```groovy title="tests/genomics-4.nf.test" linenums="6" hl_lines="1"
        test("Should run the pipeline without failures") {
    ```

=== "Antes"

    ```groovy title="tests/genomics-4.nf.test" linenums="6" hl_lines="1"
        test("Should run without failures") {
    ```

!!! note "Nota"

    Neste caso, o arquivo de teste pode permanecer onde `nf-test` o criou.

### 4.2. Especificar parâmetros de entrada

Ainda precisamos especificar entradas, o que é feito de forma ligeiramente diferente em testes em nível de fluxo de trabalho comparado a testes em nível de módulo.
Existem várias maneiras de fazer isso, incluindo especificando um perfil.
No entanto, uma maneira mais simples é configurar um bloco `params {}` no arquivo `nextflow.config` que `nf-test init` originalmente criou no diretório `tests`.

```groovy title="tests/nextflow.config" linenums="1"
/*
========================================================================================
    Nextflow config file for running tests
========================================================================================
*/

// Output directory for workflow outputs
outputDir = 'results_genomics'

/*
 * Pipeline parameters
 */

params {
    // Primary input (file of input files, one per line)
    reads_bam = "${projectDir}/data/sample_bams.txt"

    // Accessory files
    reference = "${projectDir}/data/ref/ref.fasta"
    reference_index = "${projectDir}/data/ref/ref.fasta.fai"
    reference_dict = "${projectDir}/data/ref/ref.dict"
    intervals = "${projectDir}/data/ref/intervals.bed"

    // Base name for final output file
    cohort_name = "family_trio"
}
```

Quando executarmos o teste, `nf-test` pegará este arquivo de configuração e puxará as entradas adequadamente.

### 4.3. Executar o teste de fluxo de trabalho

```bash
nf-test test tests/genomics-4.nf.test
```

??? success "Saída do comando"

    ```console
    🚀 nf-test 0.9.3
    https://www.nf-test.com
    (c) 2021 - 2024 Lukas Forer and Sebastian Schoenherr


    Test Workflow genomics-4.nf

      Test [1b4c6936] 'Should run the pipeline without failures' PASSED (171.019s)


    SUCCESS: Executed 1 tests in 171.056s
    ```

O teste passa, confirmando que nosso pipeline completo de chamada de variantes:

1. Processa com sucesso todas as amostras
2. Encadeia corretamente todas as etapas

### 4.4. Executar TODOS os testes

nf-test tem mais um truque na manga. Podemos executar todos os testes de uma vez! Modifique o arquivo `nf-test.config` para que nf-test procure em cada diretório por arquivos nf-test. Você pode fazer isso modificando o parâmetro `testsDir`:

=== "Depois"

    ```groovy title="nf-test.config" linenums="1" hl_lines="3"
    config {

        testsDir "."
        workDir ".nf-test"
        configFile "tests/nextflow.config"
        profile ""

    }
    ```

=== "Antes"

    ```groovy title="nf-test.config" linenums="1" hl_lines="3"
    config {

        testsDir "tests"
        workDir ".nf-test"
        configFile "tests/nextflow.config"
        profile ""

    }
    ```

Agora, podemos simplesmente executar nf-test e ele executará _cada teste_ em nosso repositório:

```bash
nf-test test
```

??? success "Saída do comando"

    ```console
    🚀 nf-test 0.9.3
    https://www.nf-test.com
    (c) 2021 - 2024 Lukas Forer and Sebastian Schoenherr


    Test Process GATK_HAPLOTYPECALLER

      Test [c5156c2b] 'Should call son's haplotype correctly' PASSED (39.947s)
      Test [10de94a8] 'Should call mother's haplotype correctly' PASSED (43.17s)
      Test [c0386fc7] 'Should call father's haplotype correctly' PASSED (44.244s)

    Test Process GATK_JOINTGENOTYPING

      Test [ac2067de] 'Should call trio's joint genotype correctly' PASSED (61.129s)

    Test Process SAMTOOLS_INDEX

      Test [625e39ee] 'Should index reads_son.bam correctly' PASSED (8.671s)
      Test [a8b28f36] 'Should index reads_mother.bam correctly' PASSED (8.518s)
      Test [c15852a1] 'Should index reads_father.bam correctly' PASSED (5.378s)

    Test Workflow genomics-4.nf

      Test [1b4c6936] 'Should run the pipeline without failures' PASSED (169.714s)


    SUCCESS: Executed 8 tests in 380.801s
    ```

8 testes em 1 comando! Gastamos muito tempo configurando muitos e muitos testes, mas quando chegou a hora de executá-los foi muito rápido e fácil. Você pode ver como isso é útil ao manter um pipeline grande, que pode incluir centenas de elementos diferentes. Gastamos tempo escrevendo testes uma vez para que possamos economizar tempo executando-os muitas vezes.

Além disso, podemos automatizar isso! Imagine testes executando toda vez que você ou um colega tenta adicionar novo código. É assim que garantimos que nossos pipelines mantêm um alto padrão.

## Resumo

Você agora sabe como escrever e executar vários tipos de testes para seu pipeline de genômica usando nf-test. Este framework de testes ajuda a garantir que seu fluxo de trabalho de chamada de variantes produza resultados consistentes e confiáveis em diferentes ambientes e à medida que você faz mudanças no código.

Você aprendeu a testar componentes críticos como:

- O processo `SAMTOOLS_INDEX` que prepara arquivos BAM para chamada de variantes
- O processo `GATK_HAPLOTYPECALLER` que identifica variantes em amostras individuais
- O processo `GATK_JOINTGENOTYPING` que combina chamadas de variantes em uma coorte

Você também implementou diferentes estratégias de teste específicas para dados genômicos:

- Verificar que arquivos VCF contêm chamadas de variantes esperadas apesar de elementos não determinísticos como timestamps
- Testar com um conjunto de dados de trio familiar para garantir identificação adequada de variantes entre amostras relacionadas
- Verificar coordenadas genômicas específicas e informações de variantes em seus arquivos de saída

Essas habilidades de teste são essenciais para desenvolver pipelines robustos de bioinformática que podem processar de forma confiável dados genômicos e produzir chamadas de variantes precisas. À medida que você continua trabalhando com Nextflow para análise genômica, essa base de testes ajudará você a manter código de alta qualidade que produz resultados científicos confiáveis.
