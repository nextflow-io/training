---
description: Material de treinamento básico do Nextflow
---

!!! warning

    Some of the translations on the training portal are out of date.
    The translated material may be incomplete or incorrect.
    We plan to update the translations later this year.
    In the meantime, please try to work through the English-language material if you can.

# Operadores

Operadores são métodos para conectar canais, transformar valores emitidos por canais ou executar regras próprias em canais.

Existem sete grupos de operadores descritos em detalhe na Documentação do Nextflow, estes são:

1. [Operadores de filtragem](https://www.nextflow.io/docs/latest/operator.html#filtering-operators)
2. [Operadores de transformação](https://www.nextflow.io/docs/latest/operator.html#transforming-operators)
3. [Operadores de divisão](https://www.nextflow.io/docs/latest/operator.html#splitting-operators)
4. [Operadores de combinação](https://www.nextflow.io/docs/latest/operator.html#combining-operators)
5. [Operadores de bifurcação](https://www.nextflow.io/docs/latest/operator.html#forking-operators)
6. [Operadores matemáticos](https://www.nextflow.io/docs/latest/operator.html#maths-operators)
7. [Outros operadores](https://www.nextflow.io/docs/latest/operator.html#other-operators)

## Exemplo básico

!!! info ""

    Clique no ícone :material-plus-circle: para ver explicações do código.

```groovy linenums="1"
nums = channel.of(1, 2, 3, 4) // (1)!
quadrados = nums.map { it -> it * it } // (2)!
quadrados.view() // (3)!
```

1. Cria um canal de fila que emite quatro valores
2. Cria um novo canal, transformando cada número ao quadrado
3. Imprime o conteúdo do canal

<figure class="excalidraw">
--8<-- "docs/basic_training/img/channel-map.excalidraw.pt.svg"
</figure>

Para implementar funcionalidades específicas operadores também podem ser encadeados. Então, o código anterior também pode ser escrito assim:

```groovy linenums="1"
channel
    .of(1, 2, 3, 4)
    .map { it -> it * it }
    .view()
```

## Operadores básicos

Agora iremos explorar alguns dos operadores mais comuns.

### `view()`

O operador `view` imprime os itens emitidos por um canal para o terminal, acrescentando um caractere de _quebra de linha_ após cada item. Por exemplo:

```groovy linenums="1"
channel
    .of('foo', 'bar', 'baz')
    .view()
```

```console title="Output"
foo
bar
baz
```

Você também pode especificar uma _clausura_ para personalizar como os itens são impressos. Por exemplo:

```groovy linenums="1"
channel
    .of('foo', 'bar', 'baz')
    .view { "- $it" }
```

```console title="Output"
- foo
- bar
- baz
```

### `map()`

O operador `map` aplica uma função de sua escolha em cada item emitido por um canal
e retorna os items obtidos como um novo canal. A função aplicada é chamada de função
de _mapeamento_ e é expressa com uma _clausura_, como demonstrado no exemplo abaixo:

```groovy linenums="1"
channel
    .of('olá', 'mundo')
    .map { it -> it.reverse() }
    .view()
```

Um `map` pode associar uma _tupla_ genérica a cada elemento e pode conter qualquer
tipo de dado.

```groovy linenums="1"
channel
    .of('olá', 'mundo')
    .map { palavra -> [palavra, palavra.size()] }
    .view { palavra, comprimento -> "$palavra contém $comprimento letras" }
```

!!! exercise

    Use `fromPath` para criar um canal emitindo os arquivos _fastq_ que correspondam à expressão `data/ggal/*.fq`, então use `map` para retornar um par contendo o nome e o caminho para o arquivo, e, por fim, use `view` para imprimir o canal resultante.

    ??? solution

        ```groovy linenums="1"
        channel
            .fromPath('data/ggal/*.fq')
            .map { arquivo -> [arquivo.name, arquivo] }
            .view { nome, arquivo -> "> $nome : $arquivo" }
        ```

### `mix()`

O operador `mix` combina os itens emitidos por dois (ou mais) canais em um único canal.

```groovy linenums="1"
meu_canal_1 = channel.of(1, 2, 3)
meu_canal_2 = channel.of('a', 'b')
meu_canal_3 = channel.of('z')

meu_canal_1
    .mix(meu_canal_2, meu_canal_3)
    .view()
```

```console title="Output"
1
2
a
3
b
z
```

!!! warning

    Os itens no canal resultante possuem a mesma ordem dos seus respectivos canais originais. No entanto, não há garantia que o elemento do segundo canal é acrescentado ao final dos elementos do primeiro canal. Como se pode observar acima, o elemento `a` foi impresso antes de `3`.

### `flatten()`

O operador `flatten` transforma um canal de maneira que cada _tupla_ é achatada, isto é, cada entrada é emitida como um único elemento pelo canal resultante.

```groovy linenums="1"
foo = [1, 2, 3]
bar = [4, 5, 6]

channel
    .of(foo, bar)
    .flatten()
    .view()
```

```console title="Output"
1
2
3
4
5
6
```

### `collect()`

O operador `collect` coleta todos os itens emitidos por um canal em uma lista e retorna o objeto como uma única emissão.

```groovy linenums="1"
channel
    .of(1, 2, 3, 4)
    .collect()
    .view()
```

Isto imprime o valor:

```console title="Output"
[1, 2, 3, 4]
```

!!! info

    O resultado do operador `collect` é um canal de **valor**.

### `groupTuple()`

O operador `groupTuple` coleta as tuplas (ou listas) de valores emitidos pelo canal de entrada, agrupando os elementos que possuem a mesma chave. Por fim, ele emite uma nova tupla para cada chave distinta.

Por exemplo:

```groovy linenums="1"
channel
    .of([1, 'A'], [1, 'B'], [2, 'C'], [3, 'B'], [1, 'C'], [2, 'A'], [3, 'D'])
    .groupTuple()
    .view()
```

```console title="Output"
[1, [A, B, C]]
[2, [C, A]]
[3, [B, D]]
```

Esse operador é útil para processar um grupo, juntando elementos que possuem uma propriedade ou uma chave em comum.

!!! exercise

    Use `fromPath` para criar um canal emitindo todos os arquivos no diretório `data/meta/`, então use `map` para associar o prefixo `baseName` a cada arquivo. Por fim, agrupe todo os arquivos que possuem o mesmo prefixo.

    ??? solution

        ```groovy linenums="1"
        channel
            .fromPath('data/meta/*')
            .map { arquivo -> tuple(arquivo.baseName, arquivo) }
            .groupTuple()
            .view { baseName, arquivo -> "> $baseName : $arquivo" }
        ```

### `join()`

O operador `join` cria um canal que combina os itens emitidos por dois canais que possuam uma chave em comum. Por padrão, a chave é definida como o primeiro elemento
em cada item emitido.

```groovy linenums="1"
esquerda = channel.of(['X', 1], ['Y', 2], ['Z', 3], ['P', 7])
direita = channel.of(['Z', 6], ['Y', 5], ['X', 4])
esquerda.join(direita).view()
```

```console title="Output"
[Z, 3, 6]
[Y, 2, 5]
[X, 1, 4]
```

!!! note

    Perceba como _P_ está ausente no resultado final.

### `branch()`

O operador `branch` permite que você envie os itens emitidos por um canal de entrada para um ou mais canais de saída.

O critério de seleção de cada canal de saída é definido especificando uma clausura que forneça uma ou mais expressões booleanas, cada uma das quais é identificada por um rótulo único. Para a primeira expressão verdadeira, o item é ligado a um canal nomeado com o rótulo. Por exemplo:

```groovy linenums="1"
channel
    .of(1, 2, 3, 40, 50)
    .branch {
        pequeno: it < 10
        grande: it > 10
    }
    .set { resultado }

resultado.pequeno.view { "$it é pequeno" }
resultado.grande.view { "$it é grande" }
```

!!! info

    O operador `branch` retorna um objeto multi-canal (isto é, uma variável que possui mais de um canal).

!!! note

    No exemplo acima, o que aconteceria com um valor igual a 10? Para lidar com isso, você pode usar `>=`.

## Outros recursos

Veja a [documentação de operadores](https://www.nextflow.io/docs/latest/operator.html) no site oficial do Nextflow.
