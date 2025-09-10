---
title: Introdução ao Groovy
description: Uma introdução rápida às estruturas básicas e expressões do Groovy
---

!!! warning

    Some of the translations on the training portal are out of date.
    The translated material may be incomplete or incorrect.
    We plan to update the translations later this year.
    In the meantime, please try to work through the English-language material if you can.

# Estruturas básicas e expressões do Groovy

Nextflow é uma linguagem específica de domínio (DSL) implementada sobre a linguagem de programação Groovy, que por sua vez é um superconjunto da linguagem de programação Java. Isso significa que o Nextflow pode executar qualquer código Groovy ou Java.

Aqui estão algumas sintaxes Groovy importantes que são comumente usadas no Nextflow.

## Imprimindo valores

Imprimir algo é tão fácil quanto usar um dos métodos `print` ou `println`.

```groovy linenums="1"
println("Olá, mundo!")
```

A única diferença entre os dois é que o método `println` anexa implicitamente um caractere de nova linha à string impressa.

!!! tip

    Parênteses para invocações de função são opcionais. Portanto, a seguinte sintaxe também é válida:

    ```groovy linenums="1"
    println "Olá, mundo!"
    ```

## Comentários

Os comentários usam a mesma sintaxe das linguagens de programação da família C:

```groovy linenums="1"
// comente uma única linha

/*
    um comentário abrangendo
    várias linhas
*/
```

## Variáveis

Para definir uma variável, basta atribuir um valor a ela:

```groovy linenums="1"
x = 1
println x

x = new java.util.Date()
println x

x = -3.1499392
println x

x = false
println x

x = "Oi"
println x
```

As variáveis locais são definidas usando a palavra-chave `def`:

```groovy linenums="1"
def x = 'foo'
```

O `def` deve ser sempre usado ao definir variáveis locais para uma função ou clausura.

## Listas

Um objeto List pode ser definido colocando os itens da lista entre colchetes:

```groovy linenums="1"
lista = [10, 20, 30, 40]
```

Você pode acessar um determinado item na lista com a notação de colchetes (índices começam em `0`) ou usando o método `get`:

```groovy linenums="1"
println lista[0]
println lista.get(0)
```

Para obter o comprimento de uma lista, você pode usar o método `size`:

```groovy linenums="1"
println lista.size()
```

Usamos a palavra-chave `assert` para testar se uma condição é verdadeira (semelhante a uma função `if`). Aqui, o Groovy não imprimirá nada se estiver correto, caso contrário, gerará uma mensagem AssertionError.

```groovy linenums="1"
assert lista[0] == 10
```

!!! note

    Esta afirmação deve estar correta, tente alterá-la para uma incorreta.

As listas também podem ser indexadas com índices negativos e intervalos invertidos.

```groovy linenums="1"
lista = [0, 1, 2]
assert lista[-1] == 2
assert lista[-1..0] == lista.reverse()
```

!!! info

    Na afirmação da última linha, estamos referenciando a lista inicial e convertendo-a com um intervalo "abreviado" (`..`), para executar do -1º elemento (2), o último, ao 0º elemento (0), o primeiro.

Objetos List implementam todos os métodos fornecidos pela interface [java.util.List](https://docs.oracle.com/javase/8/docs/api/java/util/List.html), mais os métodos de extensão fornecidos pelo [Groovy](http://docs.groovy-lang.org/latest/html/groovy-jdk/java/util/List.html).

```groovy linenums="1"
assert [1, 2, 3] << 1 == [1, 2, 3, 1]
assert [1, 2, 3] + [1] == [1, 2, 3, 1]
assert [1, 2, 3, 1] - [1] == [2, 3]
assert [1, 2, 3] * 2 == [1, 2, 3, 1, 2, 3]
assert [1, [2, 3]].flatten() == [1, 2, 3]
assert [1, 2, 3].reverse() == [3, 2, 1]
assert [1, 2, 3].collect { it + 3 } == [4, 5, 6]
assert [1, 2, 3, 1].unique().size() == 3
assert [1, 2, 3, 1].count(1) == 2
assert [1, 2, 3, 4].min() == 1
assert [1, 2, 3, 4].max() == 4
assert [1, 2, 3, 4].sum() == 10
assert [4, 2, 1, 3].sort() == [1, 2, 3, 4]
assert [4, 2, 1, 3].find { it % 2 == 0 } == 4
assert [4, 2, 1, 3].findAll { it % 2 == 0 } == [4, 2]
```

## Mapas

Os mapas são como listas que possuem uma chave arbitrária em vez de um número inteiro. Portanto, a sintaxe é bem parecida.

```groovy linenums="1"
mapa = [a: 0, b: 1, c: 2]
```

Os mapas podem ser acessados em uma sintaxe convencional de colchetes ou como se a chave fosse uma propriedade do mapa.

!!! info ""

    Clique no ícone :material-plus-circle: para ver explicações no código.

```groovy linenums="1"
assert mapa['a'] == 0 // (1)!
assert mapa.b == 1 // (2)!
assert mapa.get('c') == 2 // (3)!
```

1. Usando colchetes.
2. Usando a notação de ponto.
3. Usando o método `get`.

Para adicionar dados ou modificar um mapa, a sintaxe é semelhante à adição de valores a uma lista:

```groovy linenums="1"
mapa['a'] = 'x' // (1)!
mapa.b = 'y' // (2)!
mapa.put('c', 'z') // (3)!
assert mapa == [a: 'x', b: 'y', c: 'z']
```

1. Usando colchetes.
2. Usando a notação de ponto.
3. Usando o método put.

Objetos Map implementam todos os métodos fornecidos pela interface [java.util.Map](https://docs.oracle.com/javase/8/docs/api/java/util/Map.html), mais os métodos de extensão fornecidos pelo [Groovy](http://docs.groovy-lang.org/latest/html/groovy-jdk/java/util/Map.html).

## Interpolação de Strings

Strings podem ser definidas colocando-as entre aspas _simples_ ('') ou _duplas_ ("").

```groovy linenums="1"
tipoderaposa = 'rápida'
cordaraposa = ['m', 'a', 'r', 'r', 'o', 'm']
println "A $tipoderaposa raposa ${cordaraposa.join()}"

x = 'Olá'
println '$x + $y'
```

```console title="Output"
A rápida raposa marrom
$x + $y
```

!!! info

    Observe o uso diferente das sintaxes `$` e `${..}` para interpolar expressões de valor em uma string.
    A variável `$x` _não_ foi expandida, pois estava entre aspas simples.

Por fim, strings também podem ser definidas usando o caractere `/` como delimitador. Elas são conhecidas como strings **com barras** e são úteis para definir expressões regulares e padrões, pois não há necessidade de escapar as barras invertidas. Assim como as strings de aspas duplas, elas permitem interpolar variáveis prefixadas com um caractere `$`.

Tente o seguinte para ver a diferença:

```groovy linenums="1"
x = /tic\tac\toe/
y = 'tic\tac\toe'

println x
println y
```

```console title="Output"
tic\tac\toe
tic    ac    oe
```

## Strings de várias linhas

Um bloco de texto que abrange várias linhas pode ser definido delimitando-o com aspas simples ou duplas triplas:

```groovy linenums="1"
texto = """
    E aí, James.
    Como você está hoje?
    """
println texto
```

Por fim, strings de várias linhas também podem ser definidas com strings com barras. Por exemplo:

```groovy linenums="1"
texto = /
    Esta é uma string abrangendo
    várias linhas com barras!
    Super legal, né?!
    /
println texto
```

!!! info

    Como antes, strings de várias linhas dentro de aspas duplas e caracteres de barra suportam interpolação de variável, enquanto strings de várias linhas com aspas simples não.

## Declarações condicionais com if

A instrução `if` usa a mesma sintaxe comum em outras linguagens de programação, como Java, C, JavaScript, etc.

```groovy linenums="1"
if (< expressão booleana >) {
    // ramo verdadeiro
}
else {
    // ramo falso
}
```

O ramo `else` é opcional. Além disso, as chaves são opcionais quando a ramificação define apenas uma única instrução.

```groovy linenums="1"
x = 11
if (x > 10)
    println 'Olá'
```

```console title="Output"
Olá
```

!!! tip

    `null`, strings vazias e coleções (mapas e listas) vazias são avaliadas como `false`.

    Portanto, uma declaração como:

    ```groovy linenums="1"
    lista = [1, 2, 3]
    if (lista != null && lista.size() > 0) {
        println lista
    }
    else {
        println 'A lista está vazia'
    }
    ```

    Pode ser escrita como:

    ```groovy linenums="1"
    lista = [1, 2, 3]
    if (lista)
        println lista
    else
        println 'A lista está vazia'
    ```

    Veja o [Groovy-Truth](http://groovy-lang.org/semantics.html#Groovy-Truth) para mais detalhes.

!!! tip

    Em alguns casos, pode ser útil substituir a instrução `if` por uma expressão ternária (também conhecida como expressão condicional). Por exemplo:

    ```groovy linenums="1"
    println lista ? lista : 'A lista está vazia'
    ```

    A declaração anterior pode ser ainda mais simplificada usando o [operador Elvis](http://groovy-lang.org/operators.html#_elvis_operator), como mostrado abaixo:

    ```groovy linenums="1"
    println lista ?: 'A lista está vazia'
    ```

## Declarações de loop com for

A sintaxe clássica do loop `for` é suportada como mostrado aqui:

```groovy linenums="1"
for (int i = 0; i < 3; i++) {
    println("Olá mundo $i")
}
```

A iteração sobre objetos de lista também é possível usando a sintaxe abaixo:

```groovy linenums="1"
list = ['a', 'b', 'c']

for (String elem : lista) {
    println elem
}
```

## Funções

É possível definir uma função personalizada em um script, conforme mostrado aqui:

```groovy linenums="1"
def fib(int n) {
    return n < 2 ? 1 : fib(n - 1) + fib(n - 2)
}

assert fib(10)==89
```

Uma função pode receber vários argumentos, separando-os com uma vírgula. A palavra-chave `return` pode ser omitida e a função retorna implicitamente o valor da última expressão avaliada. Além disso, tipos explícitos podem ser omitidos, embora não sejam recomendados:

```groovy linenums="1"
def fact(n) {
    n > 1 ? n * fact(n - 1) : 1
}

assert fact(5) == 120
```

## Clausuras

Clausuras são o canivete suíço da programação com Nextflow/Groovy. Resumindo, uma clausura é um bloco de código que pode ser passado como um argumento para uma função. Clausuras também podem ser usadas para definir uma função anônima.

Mais formalmente, uma clausura permite a definição de funções como objetos de primeira classe.

```groovy linenums="1"
quadrado = { it * it }
```

As chaves ao redor da expressão `it * it` informam ao interpretador de scripts para tratar essa expressão como código. O identificador `it` é uma variável implícita que representa o valor que é passado para a função quando ela é invocada.

Depois de compilado, o objeto de função é atribuído à variável `quadrado` como qualquer outra atribuição de variável mostrada anteriormente. Para invocar a execução da clausura, use o método especial `call` ou simplesmente use os parênteses para especificar o(s) parâmetro(s) da clausura. Por exemplo:

```groovy linenums="1"
assert quadrado.call(5) == 25
assert quadrado(9) == 81
```

Da forma como foi mostrado, isso pode não parecer interessante, mas agora podemos passar a função `quadrado` como um argumento para outras funções ou métodos. Algumas funções embutidas aceitam uma função como esta como um argumento. Um exemplo é o método `collect` em listas:

```groovy linenums="1"
x = [1, 2, 3, 4].collect(quadrado)
println x
```

```console title="Output"
[1, 4, 9, 16]
```

Por padrão, as clausuras recebem um único parâmetro chamado `it`. Para dar a ele um nome diferente, use a sintaxe `->`. Por exemplo:

```groovy linenums="1"
quadrado = { num -> num * num }
```

Também é possível definir clausuras com vários parâmetros com nomes personalizados.

Por exemplo, quando o método `each()` é aplicado a um mapa, ele pode receber uma clausura com dois argumentos, para os quais passa o par _chave-valor_ para cada entrada no objeto `Map`. Por exemplo:

```groovy linenums="1"
imprimirMapa = { a, b -> println "$a com o valor $b" }
valores = ["Yue": "Wu", "Mark": "Williams", "Sudha": "Kumari"]
valores.each(imprimirMapa)
```

```console title="Output"
Yue com o valor Wu
Mark com o valor Williams
Sudha com o valor Kumari
```

Uma clausura tem duas outras características importantes.

Primeiro, ela pode acessar e _modificar_ variáveis no escopo em que está definida.

Em segundo lugar, uma clausura pode ser definida de maneira _anônima_, o que significa que não recebe um nome e é definida no local em que precisa ser usada.

Para um exemplo mostrando esses dois recursos, consulte o seguinte trecho de código:

```groovy linenums="1"
resultado = 0 // (1)!
valores = ["China": 1, "India": 2, "USA": 3] // (2)!
valores.keySet().each { resultado += valores[it] } // (3)!
println resultado
```

1. Define uma variável global.
2. Define um objeto de mapa.
3. Chama o método `each` passando o objeto de clausura que modifica a variável `resultado`.

Saiba mais sobre clausuras na [documentação do Groovy](http://groovy-lang.org/closures.html).

## Mais recursos

A documentação completa da linguagem Groovy está disponível [nesse link](http://groovy-lang.org/documentation.html#languagespecification).

Um ótimo recurso para dominar a sintaxe do Apache Groovy é o livro: [Groovy in Action](https://www.manning.com/books/groovy-in-action-second-edition).
