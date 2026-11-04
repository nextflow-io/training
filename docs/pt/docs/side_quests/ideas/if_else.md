# Parte 2: If - Else

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Tradução assistida por IA - [saiba mais e sugira melhorias](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

[TODO]

---

## 1. Faça a vaca citar cientistas famosos

Esta seção contém alguns exercícios de aprofundamento, para praticar o que você aprendeu até agora.
Fazer esses exercícios _não é obrigatório_ para entender as partes posteriores do treinamento, mas fornecem uma maneira divertida de reforçar seu aprendizado descobrindo como fazer a vaca citar cientistas famosos.

```console title="cowsay-output-Grace-Hopper.txt"
  _________________________________________________
 /                                                 \
| Humans are allergic to change. They love to       |
| say, 'We've always done it this way.' I try to fi |
| ght that. That's why I have a clock on my wall th |
| at runs counter-clockwise.                        |
| -Grace Hopper                                     |
 \                                                 /
  =================================================
                                                 \
                                                  \
                                                    ^__^
                                                    (oo)\_______
                                                    (__)\       )\/\
                                                        ||----w |
                                                        ||     ||
```

### 1.1. Modifique o script `hello-containers.nf` para usar um processo getQuote

Temos uma lista de pioneiros da computação e biologia no arquivo `containers/data/pioneers.csv`.
Em alto nível, para completar este exercício você precisará:

- Modificar o `params.input_file` padrão para apontar para o arquivo `pioneers.csv`.
- Criar um processo `getQuote` que usa o contêiner `quote` para buscar uma citação para cada entrada.
- Conectar a saída do processo `getQuote` ao processo `cowsay` para exibir a citação.

Para a imagem de contêiner `quote`, você pode usar aquela que você mesmo construiu no exercício de aprofundamento anterior ou usar a que você obteve do Seqera Containers.

!!! Hint "Dica"

    Uma boa escolha para o bloco `script` do seu processo getQuote pode ser:
        ```groovy
        script:
            def safe_author = author.tokenize(' ').join('-')
            """
            quote "$author" > quote-${safe_author}.txt
            echo "-${author}" >> quote-${safe_author}.txt
            """
        ```

Você pode encontrar uma solução para este exercício em `containers/solutions/hello-containers-4.1.nf`.

### 1.2. Modifique seu pipeline Nextflow para permitir que ele execute nos modos `quote` e `sayHello`.

Adicione alguma lógica de ramificação ao seu pipeline para permitir que ele aceite entradas destinadas tanto para `quote` quanto para `sayHello`.
Aqui está um exemplo de como usar uma instrução `if` em um fluxo de trabalho Nextflow:

```groovy title="hello-containers.nf"
workflow {
    if (params.quote) {
        ...
    }
    else {
        ...
    }
    cowSay(text_ch)
}
```

!!! Hint "Dica"

    Você pode usar `new_ch = processName.out` para atribuir um nome ao canal de saída de um processo.

Você pode encontrar uma solução para este exercício em `containers/solutions/hello-containers-4.2.nf`.

### Conclusão

Você sabe como usar contêineres no Nextflow para executar processos, e como construir alguma lógica de ramificação em seus pipelines!

### O que vem a seguir?

Comemore, faça uma pausa para se alongar e beba um pouco de água!

Quando estiver pronto, passe para a Parte 3 desta série de treinamento para aprender como aplicar o que você aprendeu até agora a um caso de uso de análise de dados mais realista.
